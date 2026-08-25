program test_reactive_eb_transport_2d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_2d_mod, only: initialize_reactive_2d
  use reactive_boundary_2d_mod, only: &
    reactive_boundary_set_2d, build_reactive_boundary_set_2d
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_covered_cell, eb_cut_cell
  use simulation_config_reactive_eb_2d_mod, only: reactive_eb_2d_config
  use reactive_eb_2d_driver_mod, only: &
    build_configured_eb_geometry_2d, reactive_eb_integrals_2d
  use eb_reactive_transport_2d_mod, only: &
    reactive_eb_transport_timestep_2d, advance_reactive_eb_transport_2d
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(reactive_eb_2d_config) :: config
  type(reactive_boundary_set_2d) :: boundaries
  type(eb_geometry_2d) :: geometry
  real(dp), allocatable :: state(:, :, :), temperature(:, :)
  real(dp), allocatable :: saved_state(:, :, :), saved_temperature(:, :)
  real(dp), allocatable :: initial_integrals(:), final_integrals(:)
  real(dp) :: dx, dy, base_density, dt, maximum_diffusivity, interval
  real(dp) :: initial_span, final_span, minimum_theta, scale
  logical :: ok
  integer :: i, j

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamic database load")
  call load_h2o2_elementary_transport(transport, ok)
  call require(ok, "transport database load")

  config%flow%nx = 16
  config%flow%ny = 16
  config%flow%x_lower = 0.0_dp
  config%flow%x_upper = 0.012_dp
  config%flow%y_lower = 0.0_dp
  config%flow%y_upper = 0.012_dp
  config%flow%problem = "reactive_hotspot"
  config%flow%initial_temperature = 1000.0_dp
  config%flow%initial_velocity_x = 0.0_dp
  config%flow%initial_velocity_y = 0.0_dp
  config%flow%hotspot_temperature_rise = 300.0_dp
  config%flow%hotspot_center_x = 0.006_dp
  config%flow%hotspot_center_y = 0.006_dp
  config%flow%hotspot_width = 0.001_dp
  config%flow%boundary_x_lower = "outflow"
  config%flow%boundary_x_upper = "outflow"
  config%flow%boundary_y_lower = "outflow"
  config%flow%boundary_y_upper = "outflow"
  config%geometry = "plane"
  config%plane_normal_x = 1.0_dp
  config%plane_normal_y = 0.0_dp
  config%plane_offset = 0.00437_dp
  config%state_redist_target_volume_fraction = 0.5_dp
  config%state_redist_max_order = 2

  call build_configured_eb_geometry_2d(config, geometry, ok)
  call require(ok .and. count(geometry%cell_type == eb_cut_cell) > 0, &
    "cut-cell geometry")
  call initialize_reactive_2d( &
    species, config%flow, state, temperature, dx, dy, base_density, ok)
  call require(ok .and. &
    abs(dx - geometry%dx) <= 8.0_dp * epsilon(1.0_dp) * geometry%dx .and. &
    abs(dy - geometry%dy) <= 8.0_dp * epsilon(1.0_dp) * geometry%dy, &
    "hotspot initialization")
  call build_reactive_boundary_set_2d( &
    species, config%flow, boundaries, ok)
  call require(ok, "outflow boundary construction")

  allocate(initial_integrals(size(state, 1)))
  allocate(final_integrals(size(state, 1)))
  call reactive_eb_integrals_2d(state, geometry, initial_integrals, ok)
  call require(ok, "initial EB integral")
  initial_span = maxval(temperature, &
    mask=geometry%cell_type /= eb_covered_cell) - &
    minval(temperature, mask=geometry%cell_type /= eb_covered_cell)
  allocate(saved_state, source=state)
  allocate(saved_temperature, source=temperature)

  call reactive_eb_transport_timestep_2d( &
    species, transport, state, temperature, geometry, 0.35_dp, .false., &
    .true., .false., dt, maximum_diffusivity, ok)
  call require(ok .and. dt > 0.0_dp .and. maximum_diffusivity > 0.0_dp, &
    "EB transport timestep")
  interval = min(dt, 1.0e-5_dp)
  call advance_reactive_eb_transport_2d( &
    species, transport, state, temperature, geometry, interval, .false., &
    .true., .false., .false., boundaries, &
    config%state_redist_target_volume_fraction, config%state_redist_max_order, &
    minimum_theta, ok)
  call require(ok .and. minimum_theta > 0.999999999_dp, &
    "EB conduction transaction")
  final_span = maxval(temperature, &
    mask=geometry%cell_type /= eb_covered_cell) - &
    minval(temperature, mask=geometry%cell_type /= eb_covered_cell)
  call require(final_span < initial_span, "EB conduction reduces hotspot span")
  call reactive_eb_integrals_2d(state, geometry, final_integrals, ok)
  call require(ok, "final EB integral")
  scale = max(1.0_dp, maxval(abs(initial_integrals)))
  call require(maxval(abs(final_integrals - initial_integrals)) <= &
    8.0e-12_dp * scale, "EB conduction conservation")
  do j = 1, geometry%ny
    do i = 1, geometry%nx
      if (geometry%cell_type(i, j) == eb_covered_cell) then
        call require(all(state(:, i, j) == saved_state(:, i, j)) .and. &
          temperature(i, j) == saved_temperature(i, j), &
          "covered cell remains unchanged")
      end if
    end do
  end do

  saved_state = state
  saved_temperature = temperature
  call advance_reactive_eb_transport_2d( &
    species, transport, state, temperature, geometry, -interval, .false., &
    .true., .false., .false., boundaries, &
    config%state_redist_target_volume_fraction, config%state_redist_max_order, &
    minimum_theta, ok)
  call require(.not. ok .and. all(state == saved_state) .and. &
    all(temperature == saved_temperature), "invalid interval rollback")

  write(*, '(a)') "test_reactive_eb_transport_2d: PASS"

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) error stop message
  end subroutine require

end program test_reactive_eb_transport_2d
