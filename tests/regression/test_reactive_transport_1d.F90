program test_reactive_transport_1d
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imy, iet
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density, &
    mixture_specific_gas_constant
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport
  use mixture_transport_mod, only: mixture_transport_coefficients
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_species_component, reactive_primitive_to_conserved, &
    reactive_conserved_to_primitive, reactive_transport_timestep, &
    advance_reactive_transport
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  real(dp) :: shear_errors(3), orders(2)
  integer :: grids(3), level
  logical :: ok

  grids = [32, 64, 128]
  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamics load")
  call load_h2o2_elementary_transport(transport, ok)
  call require(ok, "transport database loads")

  do level = 1, 3
    call run_shear_wave(grids(level), shear_errors(level))
  end do
  orders(1) = log(shear_errors(1) / shear_errors(2)) / log(2.0_dp)
  orders(2) = log(shear_errors(2) / shear_errors(3)) / log(2.0_dp)
  call require(all(orders > 1.75_dp), &
    "viscous shear wave is second-order accurate")
  call run_species_smoothing()
  call run_thermal_smoothing()

  write(*, '(a,3(es14.6,1x))') "Shear-wave L1 errors: ", shear_errors
  write(*, '(a,2(f10.6,1x))') "Observed orders: ", orders

contains

  subroutine run_shear_wave(nx, error)
    integer, intent(in) :: nx
    real(dp), intent(out) :: error
    real(dp), allocatable :: state(:, :), temperature(:), q(:), y(:)
    real(dp), allocatable :: diffusion(:)
    real(dp) :: xmol(7), rho, pressure, base_temperature, length, dx, x
    real(dp) :: amplitude, final_time, time, dt, max_diffusivity
    real(dp) :: viscosity, conductivity, nu, exact, local_t, c
    real(dp) :: initial_energy, final_energy, initial_momentum, final_momentum
    logical :: local_ok
    integer :: i, nvar, nprim

    nvar = reactive_nvar(size(species))
    nprim = reactive_nprim(size(species))
    allocate(state(nvar, 0:nx + 1), temperature(0:nx + 1))
    allocate(q(nprim), y(size(species)), diffusion(size(species)))
    xmol = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, 1.0e-5_dp, &
      0.0_dp, 0.55643_dp]
    call mass_fractions_from_mole_fractions(species, xmol, y, local_ok)
    call require(local_ok, "shear-wave composition converts")
    pressure = 101325.0_dp
    base_temperature = 1000.0_dp
    rho = mixture_density(species, y, pressure, base_temperature, local_ok)
    call require(local_ok, "shear-wave density evaluates")
    call mixture_transport_coefficients( &
      species, transport, y, base_temperature, pressure, viscosity, &
      conductivity, diffusion, local_ok)
    call require(local_ok, "shear-wave coefficients evaluate")
    nu = viscosity / rho
    length = 0.01_dp
    dx = length / real(nx, dp)
    amplitude = 1.0e-2_dp
    final_time = 2.0e-3_dp
    do i = 1, nx
      x = (real(i, dp) - 0.5_dp) * dx
      q(1:5) = [rho, 0.0_dp, amplitude * sin(2.0_dp * acos(-1.0_dp) * &
        x / length), 0.0_dp, pressure]
      q(6:) = y
      call reactive_primitive_to_conserved( &
        species, q, state(:, i), temperature(i), c, local_ok)
      call require(local_ok, "shear-wave state converts")
    end do
    initial_energy = dx * sum(state(iet, 1:nx))
    initial_momentum = dx * sum(state(imy, 1:nx))
    time = 0.0_dp
    do while (time < final_time - 1.0e-15_dp)
      call reactive_transport_timestep( &
        species, transport, state, temperature, nx, dx, 0.40_dp, .true., &
        .false., .false., dt, max_diffusivity, local_ok)
      call require(local_ok, "shear-wave transport timestep")
      dt = min(dt, final_time - time)
      call advance_reactive_transport( &
        species, transport, state, temperature, nx, dx, dt, "periodic", &
        .true., .false., .false., .false., local_ok)
      call require(local_ok, "shear-wave transport advance")
      time = time + dt
    end do
    error = 0.0_dp
    do i = 1, nx
      x = (real(i, dp) - 0.5_dp) * dx
      exact = amplitude * exp(-nu * (2.0_dp * acos(-1.0_dp) / length)**2 * &
        final_time) * sin(2.0_dp * acos(-1.0_dp) * x / length)
      error = error + abs(state(imy, i) / state(irho, i) - exact)
      call reactive_conserved_to_primitive( &
        species, state(:, i), temperature(i), q, local_t, c, local_ok)
      call require(local_ok, "shear-wave final state is admissible")
    end do
    error = error / real(nx, dp)
    final_energy = dx * sum(state(iet, 1:nx))
    final_momentum = dx * sum(state(imy, 1:nx))
    call require(relative_error(final_energy, initial_energy) < 3.0e-13_dp, &
      "viscous transport conserves total energy")
    call require(abs(final_momentum - initial_momentum) < 2.0e-14_dp, &
      "viscous transport conserves transverse momentum")
  end subroutine run_shear_wave

  subroutine run_species_smoothing()
    integer, parameter :: nx = 64
    real(dp), allocatable :: state(:, :), temperature(:), q(:), y(:)
    real(dp) :: xmol(7), local_xmol(7), length, dx, x, pressure, rho, c
    real(dp) :: time, final_time, dt, max_diffusivity
    real(dp) :: initial_variance, final_variance, mean_h2
    real(dp) :: initial_species(7), final_species(7)
    logical :: local_ok
    integer :: i, k, nvar

    nvar = reactive_nvar(size(species))
    allocate(state(nvar, 0:nx + 1), temperature(0:nx + 1))
    allocate(q(reactive_nprim(size(species))), y(size(species)))
    xmol = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, 1.0e-5_dp, &
      0.0_dp, 0.55643_dp]
    pressure = 101325.0_dp
    length = 0.012_dp
    dx = length / real(nx, dp)
    do i = 1, nx
      x = (real(i, dp) - 0.5_dp) * dx
      local_xmol = xmol
      local_xmol(1) = local_xmol(1) + 0.02_dp * &
        sin(2.0_dp * acos(-1.0_dp) * x / length)
      local_xmol(7) = local_xmol(7) - 0.02_dp * &
        sin(2.0_dp * acos(-1.0_dp) * x / length)
      call mass_fractions_from_mole_fractions( &
        species, local_xmol, y, local_ok)
      call require(local_ok, "species-wave composition converts")
      rho = mixture_density(species, y, pressure, 1000.0_dp, local_ok)
      call require(local_ok, "species-wave density evaluates")
      q(1:5) = [rho, 0.0_dp, 0.0_dp, 0.0_dp, pressure]
      q(6:) = y
      call reactive_primitive_to_conserved( &
        species, q, state(:, i), temperature(i), c, local_ok)
      call require(local_ok, "species-wave state converts")
    end do
    do k = 1, 7
      initial_species(k) = dx * sum(state(reactive_species_component(k), 1:nx))
    end do
    mean_h2 = sum(state(reactive_species_component(1), 1:nx) / &
      state(irho, 1:nx)) / real(nx, dp)
    initial_variance = sum((state(reactive_species_component(1), 1:nx) / &
      state(irho, 1:nx) - mean_h2)**2) / real(nx, dp)
    time = 0.0_dp
    final_time = 2.0e-5_dp
    do while (time < final_time - 1.0e-15_dp)
      call reactive_transport_timestep( &
        species, transport, state, temperature, nx, dx, 0.40_dp, .false., &
        .false., .true., dt, max_diffusivity, local_ok)
      call require(local_ok, "species-wave transport timestep")
      dt = min(dt, final_time - time)
      call advance_reactive_transport( &
        species, transport, state, temperature, nx, dx, dt, "periodic", &
        .false., .false., .true., .true., local_ok)
      call require(local_ok, "species-wave transport advance")
      time = time + dt
    end do
    mean_h2 = sum(state(reactive_species_component(1), 1:nx) / &
      state(irho, 1:nx)) / real(nx, dp)
    final_variance = sum((state(reactive_species_component(1), 1:nx) / &
      state(irho, 1:nx) - mean_h2)**2) / real(nx, dp)
    call require(final_variance < initial_variance, &
      "species diffusion reduces H2 variance")
    do k = 1, 7
      final_species(k) = dx * sum(state(reactive_species_component(k), 1:nx))
      call require(relative_error(final_species(k), initial_species(k)) < &
        5.0e-13_dp, "species diffusion conserves each periodic species mass")
      call require(minval(state(reactive_species_component(k), 1:nx)) >= &
        -1.0e-14_dp, "species diffusion preserves positivity")
    end do
  end subroutine run_species_smoothing

  subroutine run_thermal_smoothing()
    integer, parameter :: nx = 64
    real(dp), allocatable :: state(:, :), temperature(:), q(:), y(:)
    real(dp) :: xmol(7), length, dx, x, rho, r_mix, local_temperature
    real(dp) :: pressure, c, time, final_time, dt, max_diffusivity
    real(dp) :: initial_span, final_span, initial_energy, final_energy
    logical :: local_ok
    integer :: i, nvar

    nvar = reactive_nvar(size(species))
    allocate(state(nvar, 0:nx + 1), temperature(0:nx + 1))
    allocate(q(reactive_nprim(size(species))), y(size(species)))
    xmol = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, 1.0e-5_dp, &
      0.0_dp, 0.55643_dp]
    call mass_fractions_from_mole_fractions(species, xmol, y, local_ok)
    call require(local_ok, "thermal-wave composition converts")
    r_mix = mixture_specific_gas_constant(species, y, local_ok)
    call require(local_ok, "thermal-wave gas constant evaluates")
    rho = 0.30_dp
    length = 0.012_dp
    dx = length / real(nx, dp)
    do i = 1, nx
      x = (real(i, dp) - 0.5_dp) * dx
      local_temperature = 1000.0_dp + 20.0_dp * &
        sin(2.0_dp * acos(-1.0_dp) * x / length)
      pressure = rho * r_mix * local_temperature
      q(1:5) = [rho, 0.0_dp, 0.0_dp, 0.0_dp, pressure]
      q(6:) = y
      call reactive_primitive_to_conserved( &
        species, q, state(:, i), temperature(i), c, local_ok)
      call require(local_ok, "thermal-wave state converts")
    end do
    initial_span = maxval(temperature(1:nx)) - minval(temperature(1:nx))
    initial_energy = dx * sum(state(iet, 1:nx))
    time = 0.0_dp
    final_time = 2.0e-5_dp
    do while (time < final_time - 1.0e-15_dp)
      call reactive_transport_timestep( &
        species, transport, state, temperature, nx, dx, 0.40_dp, .false., &
        .true., .false., dt, max_diffusivity, local_ok)
      call require(local_ok, "thermal-wave transport timestep")
      dt = min(dt, final_time - time)
      call advance_reactive_transport( &
        species, transport, state, temperature, nx, dx, dt, "periodic", &
        .false., .true., .false., .false., local_ok)
      call require(local_ok, "thermal-wave transport advance")
      time = time + dt
    end do
    final_span = maxval(temperature(1:nx)) - minval(temperature(1:nx))
    final_energy = dx * sum(state(iet, 1:nx))
    call require(final_span < initial_span, &
      "thermal conduction reduces temperature span")
    call require(relative_error(final_energy, initial_energy) < 5.0e-13_dp, &
      "thermal conduction conserves periodic total energy")
  end subroutine run_thermal_smoothing

  pure real(dp) function relative_error(actual, expected) result(error)
    real(dp), intent(in) :: actual, expected
    error = abs(actual - expected) / max(1.0e-30_dp, abs(expected))
  end function relative_error

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') "FAILED: " // trim(message)
      error stop 1
    end if
  end subroutine require

end program test_reactive_transport_1d
