program test_reactive_contact_advection
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved, reactive_conserved_to_primitive, &
    reactive_cfl_timestep, advance_reactive_hydro, reactive_integrals
  implicit none

  type(nasa7_species), allocatable :: species(:)
  real(dp) :: hllc_error, pelec_error, rusanov_error
  real(dp) :: hllc_conservation, pelec_conservation, rusanov_conservation
  logical :: ok

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "Failed to load contact thermodynamics"

  call run_solver("hllc", hllc_error, hllc_conservation)
  call run_solver("pelec", pelec_error, pelec_conservation)
  call run_solver("rusanov", rusanov_error, rusanov_conservation)
  write(*, '(a,1x,es16.8)') "reactive contact HLLC error:", hllc_error
  write(*, '(a,1x,es16.8)') "reactive contact PeleC error:", pelec_error
  write(*, '(a,1x,es16.8)') "reactive contact Rusanov error:", rusanov_error
  write(*, '(a,3(1x,es16.8))') "contact conservation:", &
    hllc_conservation, pelec_conservation, rusanov_conservation
  if (hllc_error >= 0.75_dp * rusanov_error) then
    error stop "Reactive HLLC did not materially improve contact transport"
  end if
  if (pelec_error >= 0.75_dp * rusanov_error) then
    error stop "Reactive PeleC flux did not materially improve contact transport"
  end if
  if (max(hllc_conservation, pelec_conservation, rusanov_conservation) > 2.0e-11_dp) then
    error stop "Reactive contact transport lost conservation"
  end if
  write(*, '(a)') "test_reactive_contact_advection: PASS"

contains

  subroutine run_solver(solver, error, conservation)
    character(len=*), intent(in) :: solver
    real(dp), intent(out) :: error, conservation
    integer, parameter :: nx = 200
    real(dp), parameter :: x_lower = 0.0_dp, x_upper = 1.0_dp
    real(dp), parameter :: pressure = 101325.0_dp
    real(dp), parameter :: temperature_value = 1000.0_dp
    real(dp), parameter :: velocity = 100.0_dp
    real(dp), parameter :: final_time = 1.0e-3_dp
    real(dp), parameter :: cfl = 0.45_dp
    real(dp), allocatable :: state(:, :), temperature(:), q(:)
    real(dp) :: xl(7), xr(7), yl(7), yr(7)
    real(dp) :: rho_l, rho_r, dx, x, shifted, dt, time
    real(dp) :: local_t, sound_speed, dummy_c
    real(dp) :: initial_integrals(5), final_integrals(5)
    real(dp) :: rho_exact, y_exact, density_error, composition_error
    logical :: local_ok
    integer :: nvar, cell, k

    xl = [0.25_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.75_dp]
    xr = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
    call mass_fractions_from_mole_fractions(species, xl, yl, local_ok)
    if (.not. local_ok) error stop "Failed to convert contact left composition"
    call mass_fractions_from_mole_fractions(species, xr, yr, local_ok)
    if (.not. local_ok) error stop "Failed to convert contact right composition"
    rho_l = mixture_density(species, yl, pressure, temperature_value, local_ok)
    if (.not. local_ok) error stop "Failed to evaluate contact left density"
    rho_r = mixture_density(species, yr, pressure, temperature_value, local_ok)
    if (.not. local_ok) error stop "Failed to evaluate contact right density"

    nvar = reactive_nvar(size(species))
    allocate(state(nvar, 0:nx + 1), temperature(0:nx + 1))
    allocate(q(reactive_nprim(size(species))))
    dx = (x_upper - x_lower) / real(nx, dp)
    do cell = 1, nx
      x = x_lower + (real(cell, dp) - 0.5_dp) * dx
      if (x >= 0.25_dp .and. x < 0.50_dp) then
        call set_state(rho_l, velocity, pressure, yl, q)
      else
        call set_state(rho_r, velocity, pressure, yr, q)
      end if
      call reactive_primitive_to_conserved( &
        species, q, state(:, cell), temperature(cell), dummy_c, local_ok)
      if (.not. local_ok) error stop "Failed to initialize contact state"
    end do
    state(:, 0) = state(:, nx)
    state(:, nx + 1) = state(:, 1)
    temperature(0) = temperature(nx)
    temperature(nx + 1) = temperature(1)
    call reactive_integrals(state, nx, dx, initial_integrals)

    time = 0.0_dp
    do while (time < final_time - 1.0e-14_dp)
      call reactive_cfl_timestep( &
        species, state, temperature, nx, dx, cfl, dt, local_ok)
      if (.not. local_ok) error stop "Contact CFL evaluation failed"
      dt = min(dt, final_time - time)
      call advance_reactive_hydro( &
        species, state, temperature, nx, dx, dt, "characteristic_plm", &
        "mc", solver, "periodic", local_ok)
      if (.not. local_ok) error stop "Reactive contact hydro step failed"
      time = time + dt
    end do
    call reactive_integrals(state, nx, dx, final_integrals)
    conservation = maxval(abs(final_integrals - initial_integrals) / &
      max(1.0_dp, abs(initial_integrals)))

    density_error = 0.0_dp
    composition_error = 0.0_dp
    do cell = 1, nx
      x = x_lower + (real(cell, dp) - 0.5_dp) * dx
      shifted = modulo(x - velocity * final_time, x_upper - x_lower)
      if (shifted >= 0.25_dp .and. shifted < 0.50_dp) then
        rho_exact = rho_l
        y_exact = yl(1)
      else
        rho_exact = rho_r
        y_exact = yr(1)
      end if
      call reactive_conserved_to_primitive( &
        species, state(:, cell), temperature(cell), q, local_t, sound_speed, &
        local_ok)
      if (.not. local_ok) error stop "Invalid final contact state"
      density_error = density_error + abs(q(1) - rho_exact) / max(rho_l, rho_r)
      composition_error = composition_error + &
        abs(q(reactive_mass_fraction_component(1)) - y_exact)
    end do
    error = 0.5_dp * (density_error + composition_error) / real(nx, dp)


  end subroutine run_solver

  subroutine set_state(rho, velocity, pressure, y, primitive)
    real(dp), intent(in) :: rho, velocity, pressure, y(:)
    real(dp), intent(out) :: primitive(:)
    integer :: species_index

    primitive = 0.0_dp
    primitive(1:5) = [rho, velocity, 0.0_dp, 0.0_dp, pressure]
    do species_index = 1, size(y)
      primitive(reactive_mass_fraction_component(species_index)) = &
        y(species_index)
    end do
  end subroutine set_state

end program test_reactive_contact_advection
