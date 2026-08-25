program test_amr_eb_transport_2d
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
    amr_eb_patch_2d, build_amr_eb_patch_2d, composite_eb_integral_2d
  use amr_eb_reactive_2d_mod, only: prolong_reactive_eb_patch_pcm_2d
  use amr_eb_transport_2d_mod, only: &
    advance_two_level_reactive_eb_transport_2d
  implicit none

  integer, parameter :: coarse_nx = 12, coarse_ny = 12
  integer, parameter :: coarse_i_lower = 3, coarse_i_upper = 10
  integer, parameter :: coarse_j_lower = 3, coarse_j_upper = 10
  integer, parameter :: ratio = 2
  integer, parameter :: fine_nx = &
    (coarse_i_upper - coarse_i_lower + 1) * ratio
  integer, parameter :: fine_ny = &
    (coarse_j_upper - coarse_j_lower + 1) * ratio
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(reactive_eb_2d_config) :: config
  type(reactive_boundary_set_2d) :: boundaries
  type(eb_geometry_2d) :: coarse_geometry, fine_geometry
  type(amr_eb_patch_2d) :: patch
  real(dp), allocatable :: coarse_state(:, :, :), fine_state(:, :, :)
  real(dp), allocatable :: coarse_temperature(:, :), fine_temperature(:, :)
  real(dp), allocatable :: new_coarse_state(:, :, :)
  real(dp), allocatable :: new_fine_state(:, :, :)
  real(dp), allocatable :: new_coarse_temperature(:, :)
  real(dp), allocatable :: new_fine_temperature(:, :)
  real(dp), allocatable :: saved_coarse_state(:, :, :)
  real(dp), allocatable :: saved_fine_state(:, :, :)
  real(dp), allocatable :: integral_before(:), integral_after(:)
  real(dp) :: base_density, coarse_dx, coarse_dy
  real(dp) :: fine_x_lower, fine_x_upper, fine_y_lower, fine_y_upper
  real(dp) :: coarse_dt, fine_dt, coarse_diffusivity, fine_diffusivity
  real(dp) :: interval, initial_span, final_span, minimum_theta, scale
  logical :: ok

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamic database load")
  call load_h2o2_elementary_transport(transport, ok)
  call require(ok, "transport database load")

  config%flow%nx = coarse_nx
  config%flow%ny = coarse_ny
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
  config%plane_offset = 0.00237_dp
  config%state_redist_target_volume_fraction = 0.5_dp
  config%state_redist_max_order = 2

  call build_configured_eb_geometry_2d(config, coarse_geometry, ok)
  call require(ok .and. count(coarse_geometry%cell_type == eb_cut_cell) > 0, &
    "coarse cut-cell geometry")
  coarse_dx = (config%flow%x_upper - config%flow%x_lower) / &
    real(coarse_nx, dp)
  coarse_dy = (config%flow%y_upper - config%flow%y_lower) / &
    real(coarse_ny, dp)
  fine_x_lower = config%flow%x_lower + &
    real(coarse_i_lower - 1, dp) * coarse_dx
  fine_x_upper = config%flow%x_lower + real(coarse_i_upper, dp) * coarse_dx
  fine_y_lower = config%flow%y_lower + &
    real(coarse_j_lower - 1, dp) * coarse_dy
  fine_y_upper = config%flow%y_lower + real(coarse_j_upper, dp) * coarse_dy
  call build_configured_eb_geometry_region_2d( &
    config, fine_nx, fine_ny, fine_x_lower, fine_x_upper, fine_y_lower, &
    fine_y_upper, fine_geometry, ok)
  call require(ok .and. count(fine_geometry%cell_type == eb_cut_cell) > 0, &
    "fine cut-cell geometry")
  call build_amr_eb_patch_2d( &
    coarse_geometry, fine_geometry, coarse_i_lower, coarse_i_upper, &
    coarse_j_lower, coarse_j_upper, ratio, patch, ok)
  call require(ok, "two-level EB patch")

  call initialize_reactive_2d( &
    species, config%flow, coarse_state, coarse_temperature, coarse_dx, &
    coarse_dy, base_density, ok)
  call require(ok, "coarse hotspot initialization")
  allocate(fine_state(size(coarse_state, 1), fine_nx, fine_ny))
  allocate(fine_temperature(fine_nx, fine_ny))
  call prolong_reactive_eb_patch_pcm_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, &
    fine_geometry, patch, fine_state, fine_temperature, ok)
  call require(ok, "fine hotspot prolongation")
  call build_reactive_boundary_set_2d( &
    species, config%flow, boundaries, ok)
  call require(ok, "outflow boundary construction")

  allocate(integral_before(size(coarse_state, 1)))
  allocate(integral_after(size(coarse_state, 1)))
  call composite_eb_integral_2d( &
    coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
    integral_before, ok)
  call require(ok, "initial composite integral")
  initial_span = composite_temperature_span( &
    coarse_temperature, coarse_geometry, fine_temperature, fine_geometry)
  allocate(saved_coarse_state, source=coarse_state)
  allocate(saved_fine_state, source=fine_state)

  call reactive_eb_transport_timestep_2d( &
    species, transport, coarse_state, coarse_temperature, coarse_geometry, &
    0.30_dp, .false., .true., .false., coarse_dt, coarse_diffusivity, ok)
  call require(ok .and. coarse_dt > 0.0_dp, "coarse transport timestep")
  call reactive_eb_transport_timestep_2d( &
    species, transport, fine_state, fine_temperature, fine_geometry, &
    0.30_dp, .false., .true., .false., fine_dt, fine_diffusivity, ok)
  call require(ok .and. fine_dt > 0.0_dp, "fine transport timestep")
  interval = min(2.0e-6_dp, min(coarse_dt, real(ratio, dp) * fine_dt))

  allocate(new_coarse_state, mold=coarse_state)
  allocate(new_coarse_temperature, mold=coarse_temperature)
  allocate(new_fine_state, mold=fine_state)
  allocate(new_fine_temperature, mold=fine_temperature)
  call advance_two_level_reactive_eb_transport_2d( &
    species, transport, coarse_state, coarse_temperature, coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, patch, interval, .false., &
    .true., .false., .false., boundaries, &
    config%state_redist_target_volume_fraction, config%state_redist_max_order, &
    new_coarse_state, new_coarse_temperature, new_fine_state, &
    new_fine_temperature, minimum_theta, ok)
  call require(ok .and. minimum_theta > 0.999999999_dp, &
    "two-level EB transport transaction")
  call composite_eb_integral_2d( &
    new_coarse_state, coarse_geometry, new_fine_state, fine_geometry, patch, &
    integral_after, ok)
  call require(ok, "final composite integral")
  scale = max(1.0_dp, maxval(abs(integral_before)))
  call require(maxval(abs(integral_after - integral_before)) <= &
    2.0e-10_dp * scale, "two-level EB transport conservation")
  final_span = composite_temperature_span( &
    new_coarse_temperature, coarse_geometry, new_fine_temperature, &
    fine_geometry)
  call require(final_span < initial_span, &
    "two-level EB conduction reduces hotspot span")
  call assert_covered_unchanged( &
    saved_coarse_state, new_coarse_state, coarse_geometry, &
    "coarse covered cells unchanged")
  call assert_covered_unchanged( &
    saved_fine_state, new_fine_state, fine_geometry, &
    "fine covered cells unchanged")

  call advance_two_level_reactive_eb_transport_2d( &
    species, transport, coarse_state, coarse_temperature, coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, patch, -interval, .false., &
    .true., .false., .false., boundaries, &
    config%state_redist_target_volume_fraction, config%state_redist_max_order, &
    new_coarse_state, new_coarse_temperature, new_fine_state, &
    new_fine_temperature, minimum_theta, ok)
  call require(.not. ok .and. all(new_coarse_state == coarse_state) .and. &
    all(new_fine_state == fine_state), "invalid interval rollback")

  write(*, '(a)') "test_amr_eb_transport_2d: PASS"

contains

  real(dp) function composite_temperature_span( &
      coarse_values, coarse_grid, fine_values, fine_grid) result(span)
    real(dp), intent(in) :: coarse_values(:, :), fine_values(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_grid, fine_grid

    span = max( &
      maxval(coarse_values, mask=coarse_grid%cell_type /= eb_covered_cell), &
      maxval(fine_values, mask=fine_grid%cell_type /= eb_covered_cell)) - &
      min( &
      minval(coarse_values, mask=coarse_grid%cell_type /= eb_covered_cell), &
      minval(fine_values, mask=fine_grid%cell_type /= eb_covered_cell))
  end function composite_temperature_span

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

end program test_amr_eb_transport_2d
