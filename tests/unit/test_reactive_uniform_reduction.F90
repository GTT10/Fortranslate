program test_reactive_uniform_reduction
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: &
    initialize_reactive_1d, advance_reactive_strang, &
    advance_reactive_chemistry
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(reactive_1d_config) :: config
  real(dp), allocatable :: state(:, :), reference(:, :)
  real(dp), allocatable :: temperature(:), reference_temperature(:)
  real(dp) :: dx, dt, error
  logical :: ok

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "Failed to load uniform reactive thermodynamics"
  call load_h2o2_elementary_mechanism(reactions, ok)
  if (.not. ok) error stop "Failed to load uniform reactive mechanism"
  config%nx = 8
  config%x_lower = 0.0_dp
  config%x_upper = 0.08_dp
  config%problem = "uniform_reactor"
  config%boundary_condition = "periodic"
  config%reconstruction = "characteristic_plm"
  config%limiter = "mc"
  config%riemann_solver = "rusanov"
  config%initial_temperature = 1200.0_dp
  config%initial_pressure = 101325.0_dp
  config%initial_velocity = 0.0_dp
  call initialize_reactive_1d(species, config, state, temperature, dx, ok)
  if (.not. ok) error stop "Failed to initialize uniform reactive field"
  allocate(reference(lbound(state, 1):ubound(state, 1), &
    lbound(state, 2):ubound(state, 2)))
  allocate(reference_temperature(lbound(temperature, 1):ubound(temperature, 1)))
  reference = state
  reference_temperature = temperature
  dt = 1.0e-5_dp

  call advance_reactive_strang(species, reactions, state, temperature, &
    config%nx, dx, dt, "characteristic_plm", "mc", "rusanov", &
    "periodic", .true., &
    2.0e-7_dp, 1.0e-12_dp, ok)
  if (.not. ok) error stop "Uniform Strang step failed"
  call advance_reactive_chemistry(species, reactions, reference, &
    reference_temperature, config%nx, 0.5_dp * dt, 2.0e-7_dp, 1.0e-12_dp, &
    "periodic", ok)
  if (.not. ok) error stop "First reference chemistry half-step failed"
  call advance_reactive_chemistry(species, reactions, reference, &
    reference_temperature, config%nx, 0.5_dp * dt, 2.0e-7_dp, 1.0e-12_dp, &
    "periodic", ok)
  if (.not. ok) error stop "Second reference chemistry half-step failed"

  error = maxval(abs(state(:, 1:config%nx) - &
    reference(:, 1:config%nx)))
  if (error > 5.0e-9_dp) then
    write(*, '(a,1x,es24.16)') "uniform reduction error", error
    error stop "Uniform reactive field did not reduce to cell chemistry"
  end if
  error = maxval(abs(temperature(1:config%nx) - &
    reference_temperature(1:config%nx)))
  if (error > 5.0e-8_dp) then
    write(*, '(a,1x,es24.16)') "uniform temperature error", error
    error stop "Uniform reactive temperature reduction failed"
  end if
  write(*, '(a)') "test_reactive_uniform_reduction: PASS"
end program test_reactive_uniform_reduction
