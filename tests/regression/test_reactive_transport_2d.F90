program test_reactive_transport_2d_regression
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imx, imy, iet
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
    reactive_species_component, reactive_primitive_to_conserved
  use reactive_transport_2d_mod, only: &
    reactive_transport_timestep_2d, advance_reactive_transport_2d
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  real(dp) :: errors(3), orders(2)
  integer :: grids(3), level
  logical :: ok

  grids = [16, 32, 64]
  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamics load")
  call load_h2o2_elementary_transport(transport, ok)
  call require(ok, "transport database load")
  do level = 1, 3
    call run_diagonal_shear(grids(level), errors(level))
  end do
  orders(1) = log(errors(1) / errors(2)) / log(2.0_dp)
  orders(2) = log(errors(2) / errors(3)) / log(2.0_dp)
  call require(all(orders > 1.75_dp), &
    "two-dimensional viscous shear wave converges at second order")
  call run_species_smoothing()
  call run_thermal_smoothing()
  write(*, '(a,3(es14.6,1x))') "2D shear-wave L1 errors: ", errors
  write(*, '(a,2(f10.6,1x))') "Observed orders: ", orders

contains

  subroutine base_composition(mass_fractions, ok_out)
    real(dp), intent(out) :: mass_fractions(:)
    logical, intent(out) :: ok_out
    real(dp) :: mole_fractions(7)
    mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
      1.0e-5_dp, 0.0_dp, 0.55643_dp]
    call mass_fractions_from_mole_fractions( &
      species, mole_fractions, mass_fractions, ok_out)
  end subroutine base_composition

  subroutine run_diagonal_shear(n, error)
    integer, intent(in) :: n
    real(dp), intent(out) :: error
    real(dp), allocatable :: state(:, :, :), temperature(:, :), primitive(:)
    real(dp), allocatable :: mass_fractions(:), diffusion(:)
    real(dp) :: length, dx, dy, x, y, phase, amplitude, pressure, rho
    real(dp) :: viscosity, conductivity, nu, wave_number, exact_factor
    real(dp) :: final_time, time, dt, maximum_diffusivity, theta, sound_speed
    real(dp) :: initial_energy, final_energy, initial_mx, final_mx
    real(dp) :: initial_my, final_my, numerical_u, numerical_v
    logical :: local_ok
    integer :: i, j, k, nvar

    nvar = reactive_nvar(size(species))
    allocate(state(nvar, n, n), temperature(n, n))
    allocate(primitive(reactive_nprim(size(species))))
    allocate(mass_fractions(size(species)), diffusion(size(species)))
    call base_composition(mass_fractions, local_ok)
    call require(local_ok, "shear composition")
    pressure = 101325.0_dp
    rho = mixture_density(species, mass_fractions, pressure, 1000.0_dp, local_ok)
    call require(local_ok, "shear density")
    call mixture_transport_coefficients( &
      species, transport, mass_fractions, 1000.0_dp, pressure, viscosity, &
      conductivity, diffusion, local_ok)
    call require(local_ok, "shear transport coefficients")
    nu = viscosity / rho
    length = 0.01_dp
    dx = length / real(n, dp)
    dy = dx
    amplitude = 1.0e-2_dp
    wave_number = 2.0_dp * acos(-1.0_dp) / length
    do j = 1, n
      y = (real(j, dp) - 0.5_dp) * dy
      do i = 1, n
        x = (real(i, dp) - 0.5_dp) * dx
        phase = wave_number * (x + y)
        primitive(1:5) = [rho, amplitude * sin(phase), &
          -amplitude * sin(phase), 0.0_dp, pressure]
        do k = 1, size(species)
          primitive(reactive_mass_fraction_component(k)) = mass_fractions(k)
        end do
        call reactive_primitive_to_conserved( &
          species, primitive, state(:, i, j), temperature(i, j), sound_speed, &
          local_ok)
        call require(local_ok, "shear state")
      end do
    end do
    initial_energy = dx * dy * sum(state(iet, :, :))
    initial_mx = dx * dy * sum(state(imx, :, :))
    initial_my = dx * dy * sum(state(imy, :, :))
    time = 0.0_dp
    final_time = 5.0e-4_dp
    do while (time < final_time - 1.0e-15_dp)
      call reactive_transport_timestep_2d( &
        species, transport, state, temperature, n, n, dx, dy, 0.40_dp, &
        .true., .false., .false., dt, maximum_diffusivity, local_ok)
      call require(local_ok, "shear timestep")
      dt = min(dt, final_time - time)
      call advance_reactive_transport_2d( &
        species, transport, state, temperature, n, n, dx, dy, dt, &
        .true., .false., .false., .false., theta, local_ok)
      call require(local_ok, "shear advance")
      call require(theta > 0.999999999_dp, "shear limiter inactive")
      time = time + dt
    end do
    exact_factor = exp(-2.0_dp * nu * wave_number**2 * final_time)
    error = 0.0_dp
    do j = 1, n
      y = (real(j, dp) - 0.5_dp) * dy
      do i = 1, n
        x = (real(i, dp) - 0.5_dp) * dx
        phase = wave_number * (x + y)
        numerical_u = state(imx, i, j) / state(irho, i, j)
        numerical_v = state(imy, i, j) / state(irho, i, j)
        error = error + 0.5_dp * &
          (abs(numerical_u - amplitude * exact_factor * sin(phase)) + &
           abs(numerical_v + amplitude * exact_factor * sin(phase)))
      end do
    end do
    error = error / real(n * n, dp)
    final_energy = dx * dy * sum(state(iet, :, :))
    final_mx = dx * dy * sum(state(imx, :, :))
    final_my = dx * dy * sum(state(imy, :, :))
    call require(relative_error(final_energy, initial_energy) < 2.0e-12_dp, &
      "2D viscosity conserves total energy")
    call require(abs(final_mx - initial_mx) < 2.0e-13_dp, &
      "2D viscosity conserves x momentum")
    call require(abs(final_my - initial_my) < 2.0e-13_dp, &
      "2D viscosity conserves y momentum")
  end subroutine run_diagonal_shear

  subroutine run_species_smoothing()
    integer, parameter :: n = 32
    real(dp), allocatable :: state(:, :, :), temperature(:, :), primitive(:)
    real(dp), allocatable :: mass_fractions(:)
    real(dp) :: base_x(7), local_x(7), length, dx, dy, x, y, phase
    real(dp) :: pressure, rho, sound_speed, time, final_time, dt
    real(dp) :: maximum_diffusivity, theta, mean_h2, initial_variance
    real(dp) :: final_variance, initial_species(7), final_species(7)
    logical :: local_ok
    integer :: i, j, k, nvar

    nvar = reactive_nvar(size(species))
    allocate(state(nvar, n, n), temperature(n, n))
    allocate(primitive(reactive_nprim(size(species))))
    allocate(mass_fractions(size(species)))
    base_x = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
      1.0e-5_dp, 0.0_dp, 0.55643_dp]
    length = 0.012_dp
    dx = length / real(n, dp)
    dy = dx
    pressure = 101325.0_dp
    do j = 1, n
      y = (real(j, dp) - 0.5_dp) * dy
      do i = 1, n
        x = (real(i, dp) - 0.5_dp) * dx
        phase = sin(2.0_dp * acos(-1.0_dp) * (x + y) / length)
        local_x = base_x
        local_x(1) = local_x(1) + 0.02_dp * phase
        local_x(7) = local_x(7) - 0.02_dp * phase
        call mass_fractions_from_mole_fractions( &
          species, local_x, mass_fractions, local_ok)
        call require(local_ok, "species-wave composition")
        rho = mixture_density( &
          species, mass_fractions, pressure, 1000.0_dp, local_ok)
        call require(local_ok, "species-wave density")
        primitive(1:5) = [rho, 0.0_dp, 0.0_dp, 0.0_dp, pressure]
        do k = 1, size(species)
          primitive(reactive_mass_fraction_component(k)) = mass_fractions(k)
        end do
        call reactive_primitive_to_conserved( &
          species, primitive, state(:, i, j), temperature(i, j), &
          sound_speed, local_ok)
        call require(local_ok, "species-wave state")
      end do
    end do
    do k = 1, 7
      initial_species(k) = dx * dy * &
        sum(state(reactive_species_component(k), :, :))
    end do
    mean_h2 = sum(state(reactive_species_component(1), :, :) / &
      state(irho, :, :)) / real(n * n, dp)
    initial_variance = sum((state(reactive_species_component(1), :, :) / &
      state(irho, :, :) - mean_h2)**2) / real(n * n, dp)
    time = 0.0_dp
    final_time = 5.0e-5_dp
    do while (time < final_time - 1.0e-15_dp)
      call reactive_transport_timestep_2d( &
        species, transport, state, temperature, n, n, dx, dy, 0.40_dp, &
        .false., .false., .true., dt, maximum_diffusivity, local_ok)
      call require(local_ok, "species-wave timestep")
      dt = min(dt, final_time - time)
      call advance_reactive_transport_2d( &
        species, transport, state, temperature, n, n, dx, dy, dt, &
        .false., .false., .true., .true., theta, local_ok)
      call require(local_ok, "species-wave advance")
      time = time + dt
    end do
    mean_h2 = sum(state(reactive_species_component(1), :, :) / &
      state(irho, :, :)) / real(n * n, dp)
    final_variance = sum((state(reactive_species_component(1), :, :) / &
      state(irho, :, :) - mean_h2)**2) / real(n * n, dp)
    call require(final_variance < initial_variance, &
      "2D species diffusion reduces variance")
    do k = 1, 7
      final_species(k) = dx * dy * &
        sum(state(reactive_species_component(k), :, :))
      call require(relative_error(final_species(k), initial_species(k)) < &
        2.0e-12_dp, "2D species diffusion conserves species mass")
      call require(minval(state(reactive_species_component(k), :, :)) >= &
        -1.0e-14_dp, "2D species diffusion preserves positivity")
    end do
  end subroutine run_species_smoothing

  subroutine run_thermal_smoothing()
    integer, parameter :: n = 32
    real(dp), allocatable :: state(:, :, :), temperature(:, :), primitive(:)
    real(dp), allocatable :: mass_fractions(:)
    real(dp) :: length, dx, dy, x, y, local_temperature, rho, r_mix
    real(dp) :: pressure, sound_speed, time, final_time, dt
    real(dp) :: maximum_diffusivity, theta, initial_span, final_span
    real(dp) :: initial_energy, final_energy
    logical :: local_ok
    integer :: i, j, k, nvar

    nvar = reactive_nvar(size(species))
    allocate(state(nvar, n, n), temperature(n, n))
    allocate(primitive(reactive_nprim(size(species))))
    allocate(mass_fractions(size(species)))
    call base_composition(mass_fractions, local_ok)
    call require(local_ok, "thermal composition")
    r_mix = mixture_specific_gas_constant(species, mass_fractions, local_ok)
    call require(local_ok, "thermal gas constant")
    rho = 0.30_dp
    length = 0.012_dp
    dx = length / real(n, dp)
    dy = dx
    do j = 1, n
      y = (real(j, dp) - 0.5_dp) * dy
      do i = 1, n
        x = (real(i, dp) - 0.5_dp) * dx
        local_temperature = 1000.0_dp + 20.0_dp * &
          sin(2.0_dp * acos(-1.0_dp) * (x + y) / length)
        pressure = rho * r_mix * local_temperature
        primitive(1:5) = [rho, 0.0_dp, 0.0_dp, 0.0_dp, pressure]
        do k = 1, size(species)
          primitive(reactive_mass_fraction_component(k)) = mass_fractions(k)
        end do
        call reactive_primitive_to_conserved( &
          species, primitive, state(:, i, j), temperature(i, j), &
          sound_speed, local_ok)
        call require(local_ok, "thermal state")
      end do
    end do
    initial_span = maxval(temperature) - minval(temperature)
    initial_energy = dx * dy * sum(state(iet, :, :))
    time = 0.0_dp
    final_time = 5.0e-5_dp
    do while (time < final_time - 1.0e-15_dp)
      call reactive_transport_timestep_2d( &
        species, transport, state, temperature, n, n, dx, dy, 0.40_dp, &
        .false., .true., .false., dt, maximum_diffusivity, local_ok)
      call require(local_ok, "thermal timestep")
      dt = min(dt, final_time - time)
      call advance_reactive_transport_2d( &
        species, transport, state, temperature, n, n, dx, dy, dt, &
        .false., .true., .false., .false., theta, local_ok)
      call require(local_ok, "thermal advance")
      time = time + dt
    end do
    final_span = maxval(temperature) - minval(temperature)
    final_energy = dx * dy * sum(state(iet, :, :))
    call require(final_span < initial_span, &
      "2D conduction reduces temperature span")
    call require(relative_error(final_energy, initial_energy) < 2.0e-12_dp, &
      "2D conduction conserves total energy")
  end subroutine run_thermal_smoothing

  pure real(dp) function relative_error(value, reference) result(error)
    real(dp), intent(in) :: value, reference
    error = abs(value - reference) / max(1.0_dp, abs(reference))
  end function relative_error

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') "FAIL: " // trim(message)
      error stop
    end if
  end subroutine require
end program test_reactive_transport_2d_regression
