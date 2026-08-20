program test_reactive_vortex_2d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use simulation_config_reactive_2d_mod, only: reactive_2d_config
  use reactive_2d_mod, only: simulate_reactive_2d, reactive_extrema_2d
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(reactive_2d_config) :: config
  real(dp), allocatable :: state(:, :, :), temperature(:, :)
  real(dp) :: dx, dy, time, base_density, theta
  real(dp) :: initial_integrals(5), final_integrals(5), conservation
  real(dp) :: min_rho, max_rho, min_p, max_p, min_t, max_t, max_speed, closure
  integer :: steps
  logical :: ok

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "thermodynamic database load failed"
  call load_h2o2_elementary_mechanism(reactions, ok)
  if (.not. ok) error stop "mechanism load failed"
  config = reactive_2d_config()
  config%nx = 24
  config%ny = 24
  config%problem = "reactive_vortex"
  config%reconstruction = "characteristic_plm"
  config%riemann_solver = "hllc"
  config%chemistry_enabled = .false.
  config%use_transverse_correction = .true.
  config%initial_velocity_x = 120.0_dp
  config%initial_velocity_y = 80.0_dp
  config%vortex_strength = 45.0_dp
  config%final_time = 2.0e-6_dp
  call simulate_reactive_2d( &
    species, reactions, config, state, temperature, dx, dy, time, steps, &
    initial_integrals, final_integrals, theta, base_density, ok)
  if (.not. ok) error stop "reactive 2D vortex run failed"
  call reactive_extrema_2d( &
    species, state, temperature, config%nx, config%ny, min_rho, max_rho, &
    min_p, max_p, min_t, max_t, max_speed, closure, ok)
  if (.not. ok) error stop "reactive 2D vortex diagnostics failed"
  conservation = maxval(abs(final_integrals - initial_integrals) / &
    max(1.0_dp, abs(initial_integrals)))
  if (conservation > 2.0e-11_dp) error stop "vortex conservation failure"
  if (min_rho <= 0.0_dp .or. min_p <= 0.0_dp .or. min_t <= 0.0_dp) &
    error stop "vortex positivity failure"
  if (max_speed < 145.0_dp) error stop "vortex velocity signature missing"
  if (max_p - min_p < 20.0_dp) error stop "vortex pressure response missing"
  if (closure > 2.0e-12_dp) error stop "vortex composition closure failure"
  if (theta < 0.999_dp) error stop "unexpected vortex transverse limiting"
end program test_reactive_vortex_2d
