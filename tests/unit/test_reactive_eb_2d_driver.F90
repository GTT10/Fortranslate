program test_reactive_eb_2d_driver
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nprim, reactive_conserved_to_primitive
  use reactive_2d_mod, only: initialize_reactive_2d
  use eb_geometry_2d_mod, only: &
    eb_geometry_2d, eb_covered_cell, eb_cut_cell, eb_regular_cell
  use simulation_config_reactive_eb_2d_mod, only: reactive_eb_2d_config
  use reactive_eb_2d_driver_mod, only: &
    build_configured_eb_geometry_2d, &
    compute_reactive_eb_cfl_timestep_2d, reactive_eb_integrals_2d, &
    reactive_eb_extrema_2d, simulate_reactive_eb_2d
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(reactive_eb_2d_config) :: config
  type(eb_geometry_2d) :: geometry, simulated_geometry
  real(dp), allocatable :: initial_state(:, :, :), initial_temperature(:, :)
  real(dp), allocatable :: state(:, :, :), temperature(:, :)
  real(dp), allocatable :: primitive(:), integrals(:), expected_integrals(:)
  real(dp), allocatable :: initial_integrals(:), final_integrals(:)
  real(dp) :: dx, dy, base_density, sound_speed, local_temperature
  real(dp) :: dt, expected_dt, time, minimum_dt
  real(dp) :: minimum_density, maximum_density
  real(dp) :: minimum_pressure, maximum_pressure
  real(dp) :: minimum_temperature, maximum_temperature
  real(dp) :: maximum_speed, maximum_closure_error, scale
  logical :: ok
  integer :: i, j, steps, first_i, first_j

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamic database load")
  config%flow%nx = 20
  config%flow%ny = 20
  config%flow%x_lower = 0.0_dp
  config%flow%x_upper = 1.0_dp
  config%flow%y_lower = 0.0_dp
  config%flow%y_upper = 1.0_dp
  config%flow%final_time = 1.0e-6_dp
  config%flow%cfl = 0.2_dp
  config%flow%maximum_steps = 20
  config%flow%problem = "uniform_reactor"
  config%flow%reconstruction = "pcm"
  config%flow%riemann_solver = "hllc"
  config%flow%use_transverse_correction = .false.
  config%flow%chemistry_enabled = .false.
  config%flow%transport_enabled = .false.
  config%flow%initial_velocity_x = 0.0_dp
  config%flow%initial_velocity_y = 0.0_dp
  config%flow%boundary_x_lower = "outflow"
  config%flow%boundary_x_upper = "outflow"
  config%flow%boundary_y_lower = "outflow"
  config%flow%boundary_y_upper = "outflow"
  config%geometry = "circle"
  config%circle_center_x = 0.5_dp
  config%circle_center_y = 0.5_dp
  config%circle_radius = 0.23_dp
  config%circle_fluid_inside = .false.

  call build_configured_eb_geometry_2d(config, geometry, ok)
  call require(ok .and. geometry%is_valid(), "configured circle geometry")
  call require(count(geometry%cell_type == eb_cut_cell) > 0, "cut cells")
  call require(count(geometry%cell_type == eb_covered_cell) > 0, &
    "covered cells")
  call require(count(geometry%cell_type == eb_regular_cell) > 0, &
    "regular cells")

  call initialize_reactive_2d( &
    species, config%flow, initial_state, initial_temperature, &
    dx, dy, base_density, ok)
  call require(ok, "uniform reactive initialization")
  allocate(primitive(reactive_nprim(size(species))))
  first_i = 0
  first_j = 0
  do j = 1, geometry%ny
    do i = 1, geometry%nx
      if (geometry%cell_type(i, j) /= eb_covered_cell) then
        first_i = i
        first_j = j
        exit
      end if
    end do
    if (first_i > 0) exit
  end do
  call require(first_i > 0, "active cell search")
  call reactive_conserved_to_primitive( &
    species, initial_state(:, first_i, first_j), &
    initial_temperature(first_i, first_j), primitive, local_temperature, &
    sound_speed, ok)
  call require(ok, "active primitive recovery")
  call compute_reactive_eb_cfl_timestep_2d( &
    species, initial_state, initial_temperature, geometry, &
    config%flow%cfl, dt, ok)
  call require(ok, "EB CFL")
  expected_dt = config%flow%cfl / &
    (sound_speed / geometry%dx + sound_speed / geometry%dy)
  call assert_close(dt, expected_dt, 2.0e-13_dp, "active-cell CFL value")
  call compute_reactive_eb_cfl_timestep_2d( &
    species, initial_state, initial_temperature, geometry, -1.0_dp, dt, ok)
  call require(.not. ok .and. dt == 0.0_dp, "invalid CFL transaction")

  allocate(integrals(size(initial_state, 1)))
  allocate(expected_integrals(size(initial_state, 1)))
  call reactive_eb_integrals_2d(initial_state, geometry, integrals, ok)
  call require(ok, "volume-weighted integrals")
  expected_integrals = 0.0_dp
  do j = 1, geometry%ny
    do i = 1, geometry%nx
      expected_integrals = expected_integrals + &
        geometry%volume_fraction(i, j) * initial_state(:, i, j)
    end do
  end do
  expected_integrals = expected_integrals * geometry%dx * geometry%dy
  scale = max(1.0_dp, maxval(abs(expected_integrals)))
  call require(maxval(abs(integrals - expected_integrals)) <= &
    5.0e-14_dp * scale, "volume-weighted integral value")

  call simulate_reactive_eb_2d( &
    species, config, state, temperature, simulated_geometry, time, steps, &
    initial_integrals, final_integrals, minimum_dt, base_density, ok)
  call require(ok, "runnable EB simulation")
  call require(steps == 1 .and. time == config%flow%final_time, &
    "time loop completion")
  call require(minimum_dt == config%flow%final_time, "accepted dt diagnostic")
  scale = max(1.0_dp, maxval(abs(initial_state)))
  call require(maxval(abs(state - initial_state)) <= 2.0e-11_dp * scale, &
    "stationary state preservation")
  call require(maxval(abs(temperature - initial_temperature)) <= 2.0e-8_dp, &
    "stationary temperature preservation")
  scale = max(1.0_dp, maxval(abs(initial_integrals)))
  call require(maxval(abs(final_integrals - initial_integrals)) <= &
    2.0e-12_dp * scale, "simulation conservation")

  call reactive_eb_extrema_2d( &
    species, state, temperature, simulated_geometry, minimum_density, &
    maximum_density, minimum_pressure, maximum_pressure, &
    minimum_temperature, maximum_temperature, maximum_speed, &
    maximum_closure_error, ok)
  call require(ok, "active-cell diagnostics")
  call require(minimum_density > 0.0_dp .and. minimum_pressure > 0.0_dp, &
    "positive active extrema")
  call require(maximum_speed <= 2.0e-9_dp, "stationary maximum speed")
  call require(maximum_closure_error <= 2.0e-13_dp, &
    "composition closure")

  config%flow%chemistry_enabled = .true.
  call simulate_reactive_eb_2d( &
    species, config, state, temperature, simulated_geometry, time, steps, &
    initial_integrals, final_integrals, minimum_dt, base_density, ok)
  call require(.not. ok .and. steps == 0 .and. time == 0.0_dp, &
    "direct API rejects unsupported chemistry")
  config%flow%chemistry_enabled = .false.

  config%geometry = "plane"
  config%plane_normal_x = 1.0_dp
  config%plane_normal_y = 0.0_dp
  config%plane_offset = 0.0_dp
  call build_configured_eb_geometry_2d(config, geometry, ok)
  call require(.not. ok, "configuration must create cut cells")

  write(*, '(a)') "test_reactive_eb_2d_driver: PASS"

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

  subroutine assert_close(actual, expected, relative_tolerance, message)
    real(dp), intent(in) :: actual, expected, relative_tolerance
    character(len=*), intent(in) :: message
    real(dp) :: tolerance

    tolerance = relative_tolerance * max(1.0_dp, abs(expected))
    call require(abs(actual - expected) <= tolerance, message)
  end subroutine assert_close

end program test_reactive_eb_2d_driver
