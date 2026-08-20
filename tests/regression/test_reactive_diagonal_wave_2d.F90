program test_reactive_diagonal_wave_2d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use simulation_config_reactive_2d_mod, only: reactive_2d_config
  use reactive_2d_mod, only: &
    simulate_reactive_2d, reactive_diagonal_wave_density
  use state_indices_mod, only: irho
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp) :: errors(3), error_without_ctu, order1, order2
  integer :: grids(3), m
  logical :: ok

  grids = [12, 24, 48]
  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "thermodynamic database load failed"
  call load_h2o2_elementary_mechanism(reactions, ok)
  if (.not. ok) error stop "mechanism load failed"
  do m = 1, 3
    call run_case(grids(m), .true., errors(m))
  end do
  order1 = log(errors(1) / errors(2)) / log(2.0_dp)
  order2 = log(errors(2) / errors(3)) / log(2.0_dp)
  if (order1 < 1.45_dp .or. order2 < 1.90_dp) &
    error stop "reactive 2D diagonal-wave convergence order too low"
  call run_case(24, .false., error_without_ctu)
  if (errors(2) > 1.02_dp * error_without_ctu) &
    error stop "transverse correction materially degraded diagonal advection"
  if (abs(error_without_ctu - errors(2)) < 5.0e-7_dp) &
    error stop "transverse correction has no measurable signature"

contains

  subroutine run_case(n, use_transverse, error)
    integer, intent(in) :: n
    logical, intent(in) :: use_transverse
    real(dp), intent(out) :: error
    type(reactive_2d_config) :: config
    real(dp), allocatable :: state(:, :, :), temperature(:, :)
    real(dp) :: dx, dy, time, base_density, theta
    real(dp) :: initial_integrals(5), final_integrals(5)
    real(dp) :: x, y, exact, conservation
    integer :: i, j, steps
    logical :: local_ok

    config = reactive_2d_config()
    config%nx = n
    config%ny = n
    config%problem = "diagonal_wave"
    config%reconstruction = "characteristic_plm"
    config%riemann_solver = "hllc"
    config%limiter = "mc"
    config%use_transverse_correction = use_transverse
    config%chemistry_enabled = .false.
    config%final_time = 1.0e-6_dp
    config%cfl = 0.30_dp
    config%initial_temperature = 1000.0_dp
    config%initial_pressure = 101325.0_dp
    config%initial_velocity_x = 300.0_dp
    config%initial_velocity_y = 200.0_dp
    config%density_wave_amplitude = 0.08_dp
    call simulate_reactive_2d( &
      species, reactions, config, state, temperature, dx, dy, time, steps, &
      initial_integrals, final_integrals, theta, base_density, local_ok)
    if (.not. local_ok) error stop "reactive 2D diagonal-wave run failed"
    conservation = maxval(abs(final_integrals - initial_integrals) / &
      max(1.0_dp, abs(initial_integrals)))
    if (conservation > 2.0e-11_dp) error stop "reactive 2D conservation failure"
    if (use_transverse .and. theta < 0.999999_dp) &
      error stop "unexpected diagonal-wave transverse limiter activation"
    error = 0.0_dp
    do j = 1, n
      y = config%y_lower + (real(j, dp) - 0.5_dp) * dy
      do i = 1, n
        x = config%x_lower + (real(i, dp) - 0.5_dp) * dx
        exact = reactive_diagonal_wave_density(x, y, time, config, base_density)
        error = error + abs(state(irho, i, j) - exact)
      end do
    end do
    error = error / real(n * n, dp)
  end subroutine run_case
end program test_reactive_diagonal_wave_2d
