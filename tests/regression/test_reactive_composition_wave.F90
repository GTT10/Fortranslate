program test_reactive_composition_wave
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: &
    simulate_reactive_1d, reactive_conserved_to_primitive, reactive_nprim, &
    reactive_mass_fraction_component, reactive_composition_wave_exact
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  integer, parameter :: grids(3) = [40, 80, 160]
  real(dp) :: hllc_errors(3), hllc_density_errors(3), orders(2)
  real(dp) :: rusanov_error, rusanov_density_error
  logical :: ok
  integer :: i

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "Failed to load composition-wave thermodynamics"
  call load_h2o2_elementary_mechanism(reactions, ok)
  if (.not. ok) error stop "Failed to load composition-wave mechanism"

  do i = 1, size(grids)
    call run_grid(grids(i), "hllc", hllc_errors(i), hllc_density_errors(i))
  end do
  call run_grid(80, "rusanov", rusanov_error, rusanov_density_error)
  orders(1) = log(hllc_errors(1) / hllc_errors(2)) / log(2.0_dp)
  orders(2) = log(hllc_errors(2) / hllc_errors(3)) / log(2.0_dp)

  write(*, '(a,3(1x,es16.8))') "HLLC composition Y_H2 L1:", hllc_errors
  write(*, '(a,3(1x,es16.8))') "HLLC composition rho relative L1:", &
    hllc_density_errors
  write(*, '(a,2(1x,f10.6))') "observed orders:", orders
  write(*, '(a,1x,es16.8)') "Rusanov 80-cell Y_H2 L1:", rusanov_error

  if (minval(orders) < 1.70_dp) then
    error stop "Reactive HLLC composition wave lost second-order convergence"
  end if
  if (hllc_errors(3) > 2.0e-6_dp) then
    error stop "Reactive HLLC composition-wave error is too large"
  end if
  write(*, '(a)') "test_reactive_composition_wave: PASS"

contains

  subroutine run_grid(nx, solver, y_error, density_error)
    integer, intent(in) :: nx
    character(len=*), intent(in) :: solver
    real(dp), intent(out) :: y_error, density_error
    type(reactive_1d_config) :: config
    real(dp), allocatable :: state(:, :), temperature(:), primitive(:)
    real(dp) :: dx, time, x, exact_density, local_temperature, sound_speed
    real(dp) :: exact_y(7), initial_integrals(5), final_integrals(5)
    logical :: run_ok, local_ok
    integer :: steps, cell

    config%nx = nx
    config%x_lower = 0.0_dp
    config%x_upper = 0.012_dp
    config%final_time = 2.0e-5_dp
    config%cfl = 0.40_dp
    config%maximum_steps = 10000
    config%problem = "composition_wave"
    config%reconstruction = "characteristic_plm"
    config%riemann_solver = trim(solver)
    config%limiter = "mc"
    config%boundary_condition = "periodic"
    config%chemistry_enabled = .false.
    config%initial_temperature = 1000.0_dp
    config%initial_pressure = 101325.0_dp
    config%initial_velocity = 200.0_dp
    config%composition_wave_amplitude = 0.04_dp

    call simulate_reactive_1d(species, reactions, config, state, temperature, &
      dx, time, steps, initial_integrals, final_integrals, run_ok)
    if (.not. run_ok) error stop "Reactive composition-wave run failed"
    allocate(primitive(reactive_nprim(size(species))))
    y_error = 0.0_dp
    density_error = 0.0_dp
    do cell = 1, nx
      x = config%x_lower + (real(cell, dp) - 0.5_dp) * dx
      call reactive_composition_wave_exact( &
        species, x, time, config, exact_density, exact_y, local_ok)
      if (.not. local_ok) error stop "Composition-wave exact state failed"
      call reactive_conserved_to_primitive( &
        species, state(:, cell), temperature(cell), primitive, &
        local_temperature, sound_speed, local_ok)
      if (.not. local_ok) error stop "Composition-wave numerical state failed"
      y_error = y_error + abs( &
        primitive(reactive_mass_fraction_component(1)) - exact_y(1))
      density_error = density_error + &
        abs(state(1, cell) - exact_density) / exact_density
    end do
    y_error = y_error / real(nx, dp)
    density_error = density_error / real(nx, dp)
    if (maxval(abs(final_integrals - initial_integrals) / &
        max(1.0_dp, abs(initial_integrals))) > 5.0e-12_dp) then
      error stop "Reactive composition wave lost conservation"
    end if
  end subroutine run_grid

end program test_reactive_composition_wave
