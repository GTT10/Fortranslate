program test_reactive_characteristic_ppm_2d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use simulation_config_reactive_2d_mod, only: reactive_2d_config
  use reactive_1d_mod, only: &
    reactive_nprim, reactive_conserved_to_primitive, &
    reactive_mass_fraction_component
  use reactive_2d_mod, only: &
    simulate_reactive_2d, reactive_diagonal_wave_density, &
    reactive_diagonal_composition_wave_exact
  use state_indices_mod, only: irho
  implicit none

  integer, parameter :: grids(3) = [16, 32, 64]
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp) :: density_errors(3), composition_errors(3)
  real(dp) :: density_orders(2), composition_orders(2)
  real(dp) :: density_without_ctu
  logical :: ok
  integer :: m

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "thermodynamic database load failed"
  call load_h2o2_elementary_mechanism(reactions, ok)
  if (.not. ok) error stop "mechanism load failed"

  do m = 1, size(grids)
    call run_density_case(grids(m), .true., density_errors(m))
    call run_composition_case(grids(m), composition_errors(m))
  end do
  density_orders(1) = log(density_errors(1) / density_errors(2)) / log(2.0_dp)
  density_orders(2) = log(density_errors(2) / density_errors(3)) / log(2.0_dp)
  composition_orders(1) = &
    log(composition_errors(1) / composition_errors(2)) / log(2.0_dp)
  composition_orders(2) = &
    log(composition_errors(2) / composition_errors(3)) / log(2.0_dp)

  call run_density_case(32, .false., density_without_ctu)

  write(*, '(a,3(1x,es16.8))') &
    "2D characteristic-PPM density L1:", density_errors
  write(*, '(a,2(1x,f10.6))') "density orders:", density_orders
  write(*, '(a,3(1x,es16.8))') &
    "2D characteristic-PPM Y_H2 L1:", composition_errors
  write(*, '(a,2(1x,f10.6))') "composition orders:", composition_orders
  write(*, '(a,1x,es16.8)') &
    "32-grid density error without CTU:", density_without_ctu

  if (density_orders(1) < 1.80_dp .or. density_orders(2) < 1.90_dp) &
    error stop "2D characteristic-PPM density convergence order too low"
  if (composition_orders(1) < 1.80_dp .or. composition_orders(2) < 1.90_dp) &
    error stop "2D characteristic-PPM composition convergence order too low"
  if (density_errors(3) > 1.1e-5_dp) &
    error stop "2D characteristic-PPM density error too large"
  if (composition_errors(3) > 3.0e-6_dp) &
    error stop "2D characteristic-PPM composition error too large"
  if (density_errors(2) > 1.03_dp * density_without_ctu) &
    error stop "CTU correction materially degraded characteristic PPM"
  if (abs(density_without_ctu - density_errors(2)) < 1.0e-7_dp) &
    error stop "CTU correction lacks a characteristic-PPM signature"

contains

  subroutine configure_common(config, n)
    type(reactive_2d_config), intent(out) :: config
    integer, intent(in) :: n

    config = reactive_2d_config()
    config%nx = n
    config%ny = n
    config%reconstruction = "characteristic_ppm"
    config%riemann_solver = "hllc"
    config%limiter = "mc"
    config%chemistry_enabled = .false.
    config%final_time = 1.0e-6_dp
    config%cfl = 0.30_dp
    config%initial_temperature = 1000.0_dp
    config%initial_pressure = 101325.0_dp
    config%initial_velocity_x = 300.0_dp
    config%initial_velocity_y = 200.0_dp
    config%ppm_contact_steepening = .false.
    config%ppm_shock_flattening = .false.
  end subroutine configure_common

  subroutine run_density_case(n, use_transverse, error)
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

    call configure_common(config, n)
    config%problem = "diagonal_wave"
    config%density_wave_amplitude = 0.08_dp
    config%use_transverse_correction = use_transverse
    call simulate_reactive_2d( &
      species, reactions, config, state, temperature, dx, dy, time, steps, &
      initial_integrals, final_integrals, theta, base_density, local_ok)
    if (.not. local_ok) error stop "2D characteristic-PPM density run failed"
    conservation = maxval(abs(final_integrals - initial_integrals) / &
      max(1.0_dp, abs(initial_integrals)))
    if (conservation > 2.0e-11_dp) error stop "2D PPM density conservation failure"
    if (use_transverse .and. theta < 0.999999_dp) &
      error stop "unexpected 2D PPM transverse limiter activation"
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
  end subroutine run_density_case

  subroutine run_composition_case(n, error)
    integer, intent(in) :: n
    real(dp), intent(out) :: error
    type(reactive_2d_config) :: config
    real(dp), allocatable :: state(:, :, :), temperature(:, :)
    real(dp), allocatable :: primitive(:), exact_mass_fractions(:)
    real(dp) :: dx, dy, time, base_density, theta
    real(dp) :: initial_integrals(5), final_integrals(5)
    real(dp) :: x, y, exact_density, local_temperature, sound_speed
    real(dp) :: conservation
    integer :: i, j, steps
    logical :: local_ok

    call configure_common(config, n)
    config%problem = "diagonal_composition_wave"
    config%composition_wave_amplitude = 0.04_dp
    config%use_transverse_correction = .true.
    call simulate_reactive_2d( &
      species, reactions, config, state, temperature, dx, dy, time, steps, &
      initial_integrals, final_integrals, theta, base_density, local_ok)
    if (.not. local_ok) error stop "2D characteristic-PPM composition run failed"
    conservation = maxval(abs(final_integrals - initial_integrals) / &
      max(1.0_dp, abs(initial_integrals)))
    if (conservation > 2.0e-11_dp) &
      error stop "2D PPM composition conservation failure"
    if (theta < 0.999999_dp) &
      error stop "unexpected 2D composition transverse limiter activation"

    allocate(primitive(reactive_nprim(size(species))))
    allocate(exact_mass_fractions(size(species)))
    error = 0.0_dp
    do j = 1, n
      y = config%y_lower + (real(j, dp) - 0.5_dp) * dy
      do i = 1, n
        x = config%x_lower + (real(i, dp) - 0.5_dp) * dx
        call reactive_diagonal_composition_wave_exact( &
          species, x, y, time, config, exact_density, &
          exact_mass_fractions, local_ok)
        if (.not. local_ok) error stop "2D composition exact state failed"
        call reactive_conserved_to_primitive( &
          species, state(:, i, j), temperature(i, j), primitive, &
          local_temperature, sound_speed, local_ok)
        if (.not. local_ok) error stop "2D composition numerical state failed"
        error = error + abs( &
          primitive(reactive_mass_fraction_component(1)) - &
          exact_mass_fractions(1))
      end do
    end do
    error = error / real(n * n, dp)
  end subroutine run_composition_case
end program test_reactive_characteristic_ppm_2d
