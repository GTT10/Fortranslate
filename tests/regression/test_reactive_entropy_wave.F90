program test_reactive_entropy_wave
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use mixture_thermo_mod, only: mixture_density
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: &
    simulate_reactive_1d, reactive_entropy_wave_density
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  integer, parameter :: grids(3) = [40, 80, 160]
  real(dp) :: errors(3), orders(2), base_density
  logical :: ok
  integer :: i

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "Failed to load entropy-wave thermodynamics"
  call load_h2o2_elementary_mechanism(reactions, ok)
  if (.not. ok) error stop "Failed to load entropy-wave mechanism"
  block
    real(dp) :: y(7)
    y = 0.0_dp
    y(7) = 1.0_dp
    base_density = mixture_density(species, y, 101325.0_dp, 1200.0_dp, ok)
  end block
  if (.not. ok) error stop "Failed to compute entropy-wave density"

  do i = 1, 3
    call run_grid(grids(i), errors(i))
  end do
  orders(1) = log(errors(1) / errors(2)) / log(2.0_dp)
  orders(2) = log(errors(2) / errors(3)) / log(2.0_dp)
  write(*, '(a,3(1x,es16.8))') "reactive entropy-wave L1:", errors
  write(*, '(a,2(1x,f10.6))') "observed orders:", orders
  if (minval(orders) < 1.75_dp) then
    error stop "Reactive characteristic PLM lost second-order convergence"
  end if
  if (errors(3) > 1.0e-5_dp) then
    error stop "Reactive entropy-wave error is too large"
  end if
  write(*, '(a)') "test_reactive_entropy_wave: PASS"

contains

  subroutine run_grid(nx, l1_error)
    integer, intent(in) :: nx
    real(dp), intent(out) :: l1_error
    type(reactive_1d_config) :: config
    real(dp), allocatable :: state(:, :), temperature(:)
    real(dp) :: dx, time, x, exact
    real(dp) :: initial_integrals(5), final_integrals(5)
    integer :: steps, cell
    logical :: run_ok

    config%nx = nx
    config%x_lower = 0.0_dp
    config%x_upper = 1.0_dp
    config%final_time = 1.0e-3_dp
    config%cfl = 0.4_dp
    config%maximum_steps = 5000
    config%problem = "entropy_wave"
    config%reconstruction = "characteristic_plm"
    config%limiter = "mc"
    config%boundary_condition = "periodic"
    config%chemistry_enabled = .false.
    config%initial_temperature = 1200.0_dp
    config%initial_pressure = 101325.0_dp
    config%initial_velocity = 1000.0_dp
    config%density_wave_amplitude = 0.1_dp
    config%x_h2 = 0.0_dp
    config%x_h = 0.0_dp
    config%x_o = 0.0_dp
    config%x_o2 = 0.0_dp
    config%x_oh = 0.0_dp
    config%x_h2o = 0.0_dp
    config%x_n2 = 1.0_dp
    call simulate_reactive_1d(species, reactions, config, state, temperature, &
      dx, time, steps, initial_integrals, final_integrals, run_ok)
    if (.not. run_ok) error stop "Reactive entropy-wave run failed"
    l1_error = 0.0_dp
    do cell = 1, nx
      x = (real(cell, dp) - 0.5_dp) * dx
      exact = reactive_entropy_wave_density(x, time, config, base_density)
      l1_error = l1_error + abs(state(1, cell) - exact)
    end do
    l1_error = l1_error / real(nx, dp)
    if (maxval(abs(final_integrals - initial_integrals) / &
        max(1.0_dp, abs(initial_integrals))) > 5.0e-12_dp) then
      error stop "Reactive entropy-wave conservation failed"
    end if
  end subroutine run_grid

end program test_reactive_entropy_wave
