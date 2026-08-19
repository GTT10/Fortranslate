program test_reactive_ppm_convergence
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use mixture_thermo_mod, only: mixture_density
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: &
    simulate_reactive_1d, reactive_entropy_wave_density, &
    reactive_conserved_to_primitive, reactive_nprim, &
    reactive_mass_fraction_component, reactive_composition_wave_exact
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  integer, parameter :: grids(3) = [32, 64, 128]
  real(dp) :: entropy_error(3), composition_error(3)
  real(dp) :: entropy_order(2), composition_order(2)
  logical :: ok
  integer :: i

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "Failed to load PPM thermodynamics"
  call load_h2o2_elementary_mechanism(reactions, ok)
  if (.not. ok) error stop "Failed to load PPM mechanism"

  do i = 1, size(grids)
    call run_entropy(grids(i), entropy_error(i))
    call run_composition(grids(i), composition_error(i))
  end do
  entropy_order(1) = log(entropy_error(1) / entropy_error(2)) / log(2.0_dp)
  entropy_order(2) = log(entropy_error(2) / entropy_error(3)) / log(2.0_dp)
  composition_order(1) = &
    log(composition_error(1) / composition_error(2)) / log(2.0_dp)
  composition_order(2) = &
    log(composition_error(2) / composition_error(3)) / log(2.0_dp)

  write(*, '(a,3(1x,es16.8))') "PPM entropy-wave density L1:", entropy_error
  write(*, '(a,2(1x,f10.6))') "PPM entropy-wave orders:", entropy_order
  write(*, '(a,3(1x,es16.8))') &
    "PPM composition-wave H2 L1:", composition_error
  write(*, '(a,2(1x,f10.6))') &
    "PPM composition-wave orders:", composition_order

  if (minval(entropy_order) < 1.80_dp .or. entropy_error(3) > 2.2e-5_dp) then
    error stop "Reactive PPM entropy-wave convergence failed"
  end if
  if (minval(composition_order) < 1.60_dp .or. &
      composition_error(3) > 2.6e-6_dp) then
    error stop "Reactive PPM composition-wave convergence failed"
  end if
  write(*, '(a)') "test_reactive_ppm_convergence: PASS"

contains

  subroutine run_entropy(nx, l1_error)
    integer, intent(in) :: nx
    real(dp), intent(out) :: l1_error
    type(reactive_1d_config) :: config
    real(dp), allocatable :: state(:, :), temperature(:)
    real(dp) :: dx, time, x, exact, base_density
    real(dp) :: initial_integrals(5), final_integrals(5), y(7)
    integer :: steps, cell
    logical :: run_ok

    config = reactive_1d_config()
    config%nx = nx
    config%x_lower = 0.0_dp
    config%x_upper = 1.0_dp
    config%final_time = 1.0e-3_dp
    config%cfl = 0.4_dp
    config%maximum_steps = 10000
    config%problem = "entropy_wave"
    config%reconstruction = "ppm"
    config%riemann_solver = "hllc"
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
    y = 0.0_dp
    y(7) = 1.0_dp
    base_density = mixture_density( &
      species, y, config%initial_pressure, config%initial_temperature, run_ok)
    if (.not. run_ok) error stop "Failed to build PPM entropy density"
    call simulate_reactive_1d( &
      species, reactions, config, state, temperature, dx, time, steps, &
      initial_integrals, final_integrals, run_ok)
    if (.not. run_ok) error stop "Reactive PPM entropy-wave run failed"

    l1_error = 0.0_dp
    do cell = 1, nx
      x = (real(cell, dp) - 0.5_dp) * dx
      exact = reactive_entropy_wave_density(x, time, config, base_density)
      l1_error = l1_error + abs(state(1, cell) - exact)
    end do
    l1_error = l1_error / real(nx, dp)
    call require_conservation(initial_integrals, final_integrals)
  end subroutine run_entropy

  subroutine run_composition(nx, l1_error)
    integer, intent(in) :: nx
    real(dp), intent(out) :: l1_error
    type(reactive_1d_config) :: config
    real(dp), allocatable :: state(:, :), temperature(:), primitive(:)
    real(dp) :: dx, time, x, exact_density, local_temperature, sound_speed
    real(dp) :: initial_integrals(5), final_integrals(5)
    real(dp) :: exact_y(7)
    integer :: steps, cell
    logical :: run_ok

    config = reactive_1d_config()
    config%nx = nx
    config%x_lower = 0.0_dp
    config%x_upper = 0.012_dp
    config%final_time = 2.0e-5_dp
    config%cfl = 0.4_dp
    config%maximum_steps = 10000
    config%problem = "composition_wave"
    config%reconstruction = "ppm"
    config%riemann_solver = "hllc"
    config%limiter = "mc"
    config%boundary_condition = "periodic"
    config%chemistry_enabled = .false.
    config%initial_temperature = 1000.0_dp
    config%initial_pressure = 101325.0_dp
    config%initial_velocity = 200.0_dp
    config%composition_wave_amplitude = 0.04_dp
    config%x_h2 = 0.29570_dp
    config%x_h = 1.0e-5_dp
    config%x_o = 1.0e-5_dp
    config%x_o2 = 0.14784_dp
    config%x_oh = 1.0e-5_dp
    config%x_h2o = 0.0_dp
    config%x_n2 = 0.55643_dp
    call simulate_reactive_1d( &
      species, reactions, config, state, temperature, dx, time, steps, &
      initial_integrals, final_integrals, run_ok)
    if (.not. run_ok) error stop "Reactive PPM composition-wave run failed"

    allocate(primitive(reactive_nprim(size(species))))
    l1_error = 0.0_dp
    do cell = 1, nx
      x = config%x_lower + (real(cell, dp) - 0.5_dp) * dx
      call reactive_composition_wave_exact( &
        species, x, time, config, exact_density, exact_y, run_ok)
      if (.not. run_ok) error stop "Failed to evaluate PPM composition exact state"
      call reactive_conserved_to_primitive( &
        species, state(:, cell), temperature(cell), primitive, &
        local_temperature, sound_speed, run_ok)
      if (.not. run_ok) error stop "Invalid PPM composition-wave state"
      l1_error = l1_error + abs( &
        primitive(reactive_mass_fraction_component(1)) - exact_y(1))
    end do
    l1_error = l1_error / real(nx, dp)
    call require_conservation(initial_integrals, final_integrals)
  end subroutine run_composition

  subroutine require_conservation(initial_integrals, final_integrals)
    real(dp), intent(in) :: initial_integrals(5), final_integrals(5)
    real(dp) :: error

    error = maxval(abs(final_integrals - initial_integrals) / &
      max(1.0_dp, abs(initial_integrals)))
    if (error > 5.0e-12_dp) then
      error stop "Reactive PPM conservation failed"
    end if
  end subroutine require_conservation

end program test_reactive_ppm_convergence
