program test_amr_eb_multilevel_transport_2d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport
  use reactive_2d_mod, only: initialize_reactive_2d
  use reactive_boundary_2d_mod, only: &
    reactive_boundary_set_2d, build_reactive_boundary_set_2d
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_covered_cell, eb_cut_cell
  use simulation_config_reactive_eb_2d_mod, only: reactive_eb_2d_config
  use reactive_eb_2d_driver_mod, only: &
    build_configured_eb_geometry_2d, &
    build_configured_eb_geometry_region_2d
  use eb_reactive_transport_2d_mod, only: &
    reactive_eb_transport_timestep_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, build_amr_eb_patch_2d
  use amr_eb_multilevel_2d_mod, only: composite_three_level_eb_integral_2d
  use amr_eb_reactive_2d_mod, only: prolong_reactive_eb_patch_pcm_2d
  use amr_eb_multilevel_transport_2d_mod, only: &
    advance_three_level_reactive_eb_transport_2d
  implicit none

  integer, parameter :: root_nx = 8, root_ny = 8
  integer, parameter :: root_i_lower = 2, root_i_upper = 7
  integer, parameter :: root_j_lower = 2, root_j_upper = 7
  integer, parameter :: level_one_nx = 12, level_one_ny = 12
  integer, parameter :: level_one_i_lower = 3, level_one_i_upper = 10
  integer, parameter :: level_one_j_lower = 3, level_one_j_upper = 10
  integer, parameter :: level_two_nx = 16, level_two_ny = 16
  integer, parameter :: ratio = 2
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(reactive_eb_2d_config) :: config
  type(reactive_boundary_set_2d) :: boundaries
  type(eb_geometry_2d) :: root_geometry, level_one_geometry
  type(eb_geometry_2d) :: level_two_geometry
  type(amr_eb_patch_2d) :: root_patch, level_one_patch
  real(dp), allocatable :: root_state(:, :, :), root_temperature(:, :)
  real(dp), allocatable :: level_one_state(:, :, :)
  real(dp), allocatable :: level_one_temperature(:, :)
  real(dp), allocatable :: level_two_state(:, :, :)
  real(dp), allocatable :: level_two_temperature(:, :)
  real(dp), allocatable :: new_root_state(:, :, :)
  real(dp), allocatable :: new_root_temperature(:, :)
  real(dp), allocatable :: new_level_one_state(:, :, :)
  real(dp), allocatable :: new_level_one_temperature(:, :)
  real(dp), allocatable :: new_level_two_state(:, :, :)
  real(dp), allocatable :: new_level_two_temperature(:, :)
  real(dp), allocatable :: integral_before(:), integral_after(:)
  real(dp) :: root_dx, root_dy, level_one_dx, level_one_dy
  real(dp) :: root_x_lower, root_x_upper, root_y_lower, root_y_upper
  real(dp) :: one_x_lower, one_x_upper, one_y_lower, one_y_upper
  real(dp) :: two_x_lower, two_x_upper, two_y_lower, two_y_upper
  real(dp) :: base_density, root_dt, one_dt, two_dt, diffusivity
  real(dp) :: interval, initial_span, final_span, minimum_theta, scale
  logical :: ok

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamic database load")
  call load_h2o2_elementary_transport(transport, ok)
  call require(ok, "transport database load")

  config%flow%nx = root_nx
  config%flow%ny = root_ny
  config%flow%x_lower = 0.0_dp
  config%flow%x_upper = 0.012_dp
  config%flow%y_lower = 0.0_dp
  config%flow%y_upper = 0.012_dp
  config%flow%problem = "reactive_hotspot"
  config%flow%initial_temperature = 1000.0_dp
  config%flow%initial_velocity_x = 0.0_dp
  config%flow%initial_velocity_y = 0.0_dp
  config%flow%hotspot_temperature_rise = 350.0_dp
  config%flow%hotspot_center_x = 0.006_dp
  config%flow%hotspot_center_y = 0.006_dp
  config%flow%hotspot_width = 0.0015_dp
  config%flow%boundary_x_lower = "outflow"
  config%flow%boundary_x_upper = "outflow"
  config%flow%boundary_y_lower = "outflow"
  config%flow%boundary_y_upper = "outflow"
  config%geometry = "plane"
  config%plane_normal_x = 1.0_dp
  config%plane_normal_y = 0.0_dp
  config%plane_offset = 0.00337_dp
  config%state_redist_target_volume_fraction = 0.5_dp
  config%state_redist_max_order = 2

  root_x_lower = config%flow%x_lower
  root_x_upper = config%flow%x_upper
  root_y_lower = config%flow%y_lower
  root_y_upper = config%flow%y_upper
  call build_configured_eb_geometry_2d(config, root_geometry, ok)
  call require(ok .and. count(root_geometry%cell_type == eb_cut_cell) > 0, &
    "root cut-cell geometry")
  root_dx = (root_x_upper - root_x_lower) / real(root_nx, dp)
  root_dy = (root_y_upper - root_y_lower) / real(root_ny, dp)
  one_x_lower = root_x_lower + real(root_i_lower - 1, dp) * root_dx
  one_x_upper = root_x_lower + real(root_i_upper, dp) * root_dx
  one_y_lower = root_y_lower + real(root_j_lower - 1, dp) * root_dy
  one_y_upper = root_y_lower + real(root_j_upper, dp) * root_dy
  call build_configured_eb_geometry_region_2d( &
    config, level_one_nx, level_one_ny, one_x_lower, one_x_upper, &
    one_y_lower, one_y_upper, level_one_geometry, ok)
  call require(ok .and. &
    count(level_one_geometry%cell_type == eb_cut_cell) > 0, &
    "level-one cut-cell geometry")
  call build_amr_eb_patch_2d( &
    root_geometry, level_one_geometry, root_i_lower, root_i_upper, &
    root_j_lower, root_j_upper, ratio, root_patch, ok)
  call require(ok, "root patch")

  level_one_dx = (one_x_upper - one_x_lower) / real(level_one_nx, dp)
  level_one_dy = (one_y_upper - one_y_lower) / real(level_one_ny, dp)
  two_x_lower = one_x_lower + &
    real(level_one_i_lower - 1, dp) * level_one_dx
  two_x_upper = one_x_lower + real(level_one_i_upper, dp) * level_one_dx
  two_y_lower = one_y_lower + &
    real(level_one_j_lower - 1, dp) * level_one_dy
  two_y_upper = one_y_lower + real(level_one_j_upper, dp) * level_one_dy
  call build_configured_eb_geometry_region_2d( &
    config, level_two_nx, level_two_ny, two_x_lower, two_x_upper, &
    two_y_lower, two_y_upper, level_two_geometry, ok)
  call require(ok .and. &
    count(level_two_geometry%cell_type == eb_cut_cell) > 0, &
    "level-two cut-cell geometry")
  call build_amr_eb_patch_2d( &
    level_one_geometry, level_two_geometry, level_one_i_lower, &
    level_one_i_upper, level_one_j_lower, level_one_j_upper, ratio, &
    level_one_patch, ok)
  call require(ok, "level-one patch")

  call initialize_reactive_2d( &
    species, config%flow, root_state, root_temperature, root_dx, root_dy, &
    base_density, ok)
  call require(ok, "root hotspot initialization")
  allocate(level_one_state(size(root_state, 1), level_one_nx, level_one_ny))
  allocate(level_one_temperature(level_one_nx, level_one_ny))
  call prolong_reactive_eb_patch_pcm_2d( &
    species, root_state, root_temperature, root_geometry, level_one_geometry, &
    root_patch, level_one_state, level_one_temperature, ok)
  call require(ok, "level-one prolongation")
  allocate(level_two_state(size(root_state, 1), level_two_nx, level_two_ny))
  allocate(level_two_temperature(level_two_nx, level_two_ny))
  call prolong_reactive_eb_patch_pcm_2d( &
    species, level_one_state, level_one_temperature, level_one_geometry, &
    level_two_geometry, level_one_patch, level_two_state, &
    level_two_temperature, ok)
  call require(ok, "level-two prolongation")
  call build_reactive_boundary_set_2d( &
    species, config%flow, boundaries, ok)
  call require(ok, "outflow boundary construction")

  allocate(integral_before(size(root_state, 1)))
  allocate(integral_after(size(root_state, 1)))
  call composite_three_level_eb_integral_2d( &
    root_state, root_geometry, level_one_state, level_one_geometry, &
    root_patch, level_two_state, level_two_geometry, level_one_patch, &
    integral_before, ok)
  call require(ok, "initial three-level integral")
  initial_span = temperature_span( &
    root_temperature, root_geometry, level_one_temperature, &
    level_one_geometry, level_two_temperature, level_two_geometry)

  call reactive_eb_transport_timestep_2d( &
    species, transport, root_state, root_temperature, root_geometry, 0.30_dp, &
    .false., .true., .false., root_dt, diffusivity, ok)
  call require(ok, "root transport timestep")
  call reactive_eb_transport_timestep_2d( &
    species, transport, level_one_state, level_one_temperature, &
    level_one_geometry, 0.30_dp, .false., .true., .false., one_dt, &
    diffusivity, ok)
  call require(ok, "level-one transport timestep")
  call reactive_eb_transport_timestep_2d( &
    species, transport, level_two_state, level_two_temperature, &
    level_two_geometry, 0.30_dp, .false., .true., .false., two_dt, &
    diffusivity, ok)
  call require(ok, "level-two transport timestep")
  interval = min(2.0e-6_dp, min(root_dt, &
    min(real(ratio, dp) * one_dt, real(ratio * ratio, dp) * two_dt)))

  allocate(new_root_state, mold=root_state)
  allocate(new_root_temperature, mold=root_temperature)
  allocate(new_level_one_state, mold=level_one_state)
  allocate(new_level_one_temperature, mold=level_one_temperature)
  allocate(new_level_two_state, mold=level_two_state)
  allocate(new_level_two_temperature, mold=level_two_temperature)
  call advance_three_level_reactive_eb_transport_2d( &
    species, transport, root_state, root_temperature, root_geometry, &
    level_one_state, level_one_temperature, level_one_geometry, root_patch, &
    level_two_state, level_two_temperature, level_two_geometry, &
    level_one_patch, interval, .false., .true., .false., .false., &
    boundaries, config%state_redist_target_volume_fraction, &
    config%state_redist_max_order, new_root_state, new_root_temperature, &
    new_level_one_state, new_level_one_temperature, new_level_two_state, &
    new_level_two_temperature, minimum_theta, ok)
  call require(ok .and. minimum_theta > 0.999999999_dp, &
    "three-level transport transaction")
  call composite_three_level_eb_integral_2d( &
    new_root_state, root_geometry, new_level_one_state, level_one_geometry, &
    root_patch, new_level_two_state, level_two_geometry, level_one_patch, &
    integral_after, ok)
  scale = max(1.0_dp, maxval(abs(integral_before)))
  call require(ok .and. maxval(abs(integral_after - integral_before)) <= &
    3.0e-10_dp * scale, "three-level transport conservation")
  final_span = temperature_span( &
    new_root_temperature, root_geometry, new_level_one_temperature, &
    level_one_geometry, new_level_two_temperature, level_two_geometry)
  call require(final_span < initial_span, &
    "three-level conduction reduces temperature span")
  call assert_covered_unchanged( &
    root_state, new_root_state, root_geometry, "root covered state")
  call assert_covered_unchanged( &
    level_one_state, new_level_one_state, level_one_geometry, &
    "level-one covered state")
  call assert_covered_unchanged( &
    level_two_state, new_level_two_state, level_two_geometry, &
    "level-two covered state")

  call advance_three_level_reactive_eb_transport_2d( &
    species, transport, root_state, root_temperature, root_geometry, &
    level_one_state, level_one_temperature, level_one_geometry, root_patch, &
    level_two_state, level_two_temperature, level_two_geometry, &
    level_one_patch, -interval, .false., .true., .false., .false., &
    boundaries, config%state_redist_target_volume_fraction, &
    config%state_redist_max_order, new_root_state, new_root_temperature, &
    new_level_one_state, new_level_one_temperature, new_level_two_state, &
    new_level_two_temperature, minimum_theta, ok)
  call require(.not. ok .and. all(new_root_state == root_state) .and. &
    all(new_level_one_state == level_one_state) .and. &
    all(new_level_two_state == level_two_state), &
    "three-level transport rollback")

  write(*, '(a)') "test_amr_eb_multilevel_transport_2d: PASS"

contains

  real(dp) function temperature_span( &
      root_values, root_grid, one_values, one_grid, two_values, two_grid) &
      result(span)
    real(dp), intent(in) :: root_values(:, :), one_values(:, :)
    real(dp), intent(in) :: two_values(:, :)
    type(eb_geometry_2d), intent(in) :: root_grid, one_grid, two_grid

    span = max( &
      maxval(root_values, mask=root_grid%cell_type /= eb_covered_cell), &
      maxval(one_values, mask=one_grid%cell_type /= eb_covered_cell), &
      maxval(two_values, mask=two_grid%cell_type /= eb_covered_cell)) - &
      min( &
      minval(root_values, mask=root_grid%cell_type /= eb_covered_cell), &
      minval(one_values, mask=one_grid%cell_type /= eb_covered_cell), &
      minval(two_values, mask=two_grid%cell_type /= eb_covered_cell))
  end function temperature_span

  subroutine assert_covered_unchanged(original, candidate, geometry, message)
    real(dp), intent(in) :: original(:, :, :), candidate(:, :, :)
    type(eb_geometry_2d), intent(in) :: geometry
    character(len=*), intent(in) :: message
    integer :: i, j

    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) /= eb_covered_cell) cycle
        call require(all(candidate(:, i, j) == original(:, i, j)), message)
      end do
    end do
  end subroutine assert_covered_unchanged

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) error stop message
  end subroutine require

end program test_amr_eb_multilevel_transport_2d
