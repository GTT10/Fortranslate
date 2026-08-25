program test_reactive_eb_amr_2d_driver
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_cut_cell
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, composite_eb_integral_2d
  use simulation_config_reactive_eb_amr_2d_mod, only: &
    reactive_eb_amr_2d_config
  use reactive_eb_amr_2d_driver_mod, only: &
    compute_reactive_eb_amr_cfl_timestep_2d, &
    regrid_reactive_eb_amr_hierarchy_2d, &
    simulate_reactive_eb_amr_2d
  use reactive_eb_2d_driver_mod, only: reactive_eb_integrals_2d
  implicit none

  type(reactive_eb_amr_2d_config) :: config
  type(eb_geometry_2d) :: coarse_geometry, fine_geometry
  type(amr_eb_patch_2d) :: patch
  type(nasa7_species), allocatable :: species(:)
  real(dp), allocatable :: coarse_state(:, :, :), coarse_temperature(:, :)
  real(dp), allocatable :: fine_state(:, :, :), fine_temperature(:, :)
  real(dp), allocatable :: initial_integrals(:), final_integrals(:)
  real(dp), allocatable :: lifecycle_integrals(:)
  real(dp), allocatable :: reference_state(:)
  real(dp) :: time, minimum_dt, base_density, cfl_dt, conservation_error, scale
  logical :: changed, fine_active, ok
  integer :: initial_i_lower, initial_i_upper
  integer :: initial_j_lower, initial_j_upper, regrids, steps

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamic database load")
  config%eb%flow%nx = 8
  config%eb%flow%ny = 8
  config%eb%flow%x_lower = 0.0_dp
  config%eb%flow%x_upper = 1.0_dp
  config%eb%flow%y_lower = 0.0_dp
  config%eb%flow%y_upper = 1.0_dp
  config%eb%flow%final_time = 1.0e-6_dp
  config%eb%flow%cfl = 0.2_dp
  config%eb%flow%maximum_steps = 20
  config%eb%flow%problem = "uniform_reactor"
  config%eb%flow%reconstruction = "characteristic_plm"
  config%eb%flow%limiter = "mc"
  config%eb%flow%riemann_solver = "hllc"
  config%eb%flow%use_transverse_correction = .false.
  config%eb%flow%chemistry_enabled = .false.
  config%eb%flow%transport_enabled = .false.
  config%eb%flow%boundary_x_lower = "outflow"
  config%eb%flow%boundary_x_upper = "outflow"
  config%eb%flow%boundary_y_lower = "outflow"
  config%eb%flow%boundary_y_upper = "outflow"
  config%eb%flow%initial_temperature = 1000.0_dp
  config%eb%flow%initial_pressure = 101325.0_dp
  config%eb%flow%initial_velocity_x = 0.0_dp
  config%eb%flow%initial_velocity_y = 0.0_dp
  config%eb%geometry = "plane"
  config%eb%plane_normal_x = 1.0_dp
  config%eb%plane_normal_y = 1.0_dp
  config%eb%plane_offset = 0.78_dp
  config%eb%state_redist_target_volume_fraction = 0.5_dp
  config%eb%state_redist_max_order = 2
  config%coarse_i_lower = 2
  config%coarse_i_upper = 6
  config%coarse_j_lower = 2
  config%coarse_j_upper = 6
  config%refinement_ratio = 2

  call simulate_reactive_eb_amr_2d( &
    species, config, coarse_state, coarse_temperature, coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, patch, fine_active, time, &
    steps, regrids, initial_integrals, final_integrals, minimum_dt, &
    base_density, ok)
  call require(ok, "runnable static EB AMR simulation")
  call require(fine_active .and. steps == 1 .and. regrids == 0 .and. &
    time == config%eb%flow%final_time .and. &
    minimum_dt == config%eb%flow%final_time, "time-loop completion")
  call require(coarse_geometry%nx == 8 .and. coarse_geometry%ny == 8 .and. &
    fine_geometry%nx == 10 .and. fine_geometry%ny == 10, &
    "two-level dimensions")
  call require(patch%is_valid(coarse_geometry, fine_geometry), &
    "qualified static patch")
  call require(count(coarse_geometry%cell_type == eb_cut_cell) > 0 .and. &
    count(fine_geometry%cell_type == eb_cut_cell) > 0, &
    "two-level cut-cell coverage")
  scale = max(1.0_dp, maxval(abs(initial_integrals)))
  call require(maxval(abs(final_integrals - initial_integrals)) <= &
    3.0e-12_dp * scale, "static EB AMR composite conservation")
  allocate(reference_state, source=coarse_state(:, 8, 8))
  scale = max(1.0_dp, maxval(abs(reference_state)))
  call require(maxval(abs(coarse_state - &
    spread(spread(reference_state, 2, 8), 3, 8))) <= &
    3.0e-12_dp * scale, "coarse stationary state")
  call require(maxval(abs(fine_state - &
    spread(spread(reference_state, 2, 10), 3, 10))) <= &
    3.0e-12_dp * scale, "fine stationary state")
  call require(maxval(abs(coarse_temperature - 1000.0_dp)) <= 3.0e-8_dp .and. &
    maxval(abs(fine_temperature - 1000.0_dp)) <= 3.0e-8_dp, &
    "two-level stationary temperature")
  call compute_reactive_eb_amr_cfl_timestep_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, config%refinement_ratio, &
    config%eb%flow%cfl, cfl_dt, ok)
  call require(ok .and. cfl_dt > config%eb%flow%final_time, &
    "two-level CFL selection")

  config%eb%flow%nx = 12
  config%eb%flow%ny = 12
  config%eb%flow%final_time = 1.0e-8_dp
  config%eb%flow%maximum_steps = 10
  config%eb%flow%problem = "reactive_hotspot"
  config%eb%flow%hotspot_temperature_rise = 350.0_dp
  config%eb%flow%hotspot_center_x = 0.72_dp
  config%eb%flow%hotspot_center_y = 0.62_dp
  config%eb%flow%hotspot_width = 0.08_dp
  config%eb%plane_offset = 0.30_dp
  config%coarse_i_lower = 2
  config%coarse_i_upper = 5
  config%coarse_j_lower = 2
  config%coarse_j_upper = 5
  config%dynamic_regridding = .true.
  config%regrid_interval = 1
  config%regrid_relative_temperature_gradient = 0.02_dp
  config%regrid_absolute_temperature_gradient = 5.0_dp
  config%regrid_temperature_scale_floor = 1.0_dp
  config%regrid_buffer_cells = 1
  config%regrid_minimum_patch_cells_x = 4
  config%regrid_minimum_patch_cells_y = 4
  initial_i_lower = config%coarse_i_lower
  initial_i_upper = config%coarse_i_upper
  initial_j_lower = config%coarse_j_lower
  initial_j_upper = config%coarse_j_upper
  call simulate_reactive_eb_amr_2d( &
    species, config, coarse_state, coarse_temperature, coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, patch, fine_active, time, &
    steps, regrids, initial_integrals, final_integrals, minimum_dt, &
    base_density, ok)
  call require(ok .and. fine_active .and. regrids >= 1, &
    "dynamic EB AMR simulation")
  call require(patch%is_valid(coarse_geometry, fine_geometry), &
    "dynamically selected patch")
  call require(patch%coarse_i_lower /= initial_i_lower .or. &
    patch%coarse_i_upper /= initial_i_upper .or. &
    patch%coarse_j_lower /= initial_j_lower .or. &
    patch%coarse_j_upper /= initial_j_upper, "hotspot moves static patch")
  scale = max(1.0_dp, maxval(abs(initial_integrals)))
  call require(maxval(abs(final_integrals - initial_integrals)) <= &
    2.0e-10_dp * scale, "dynamic EB AMR composite conservation")

  config%eb%flow%problem = "uniform_reactor"
  config%eb%flow%initial_temperature = 1000.0_dp
  config%eb%flow%initial_velocity_x = 0.0_dp
  config%eb%flow%initial_velocity_y = 0.0_dp
  config%eb%flow%final_time = 1.0e-8_dp
  config%regrid_at_initialization = .false.
  config%remove_fine_patch_when_untagged = .true.
  call simulate_reactive_eb_amr_2d( &
    species, config, coarse_state, coarse_temperature, coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, patch, fine_active, time, &
    steps, regrids, initial_integrals, final_integrals, minimum_dt, &
    base_density, ok)
  call require(ok .and. .not. fine_active .and. steps == 1 .and. &
    regrids == 1, "time-loop fine-patch removal")
  call require(.not. allocated(fine_state) .and. &
    .not. allocated(fine_temperature) .and. &
    .not. fine_geometry%is_valid() .and. patch%refinement_ratio == 0, &
    "inactive fine storage released")
  scale = max(1.0_dp, maxval(abs(initial_integrals)))
  conservation_error = maxval(abs(final_integrals - initial_integrals)) / scale
  write(*, '(a,1x,es16.8)') &
    "Fine-patch removal conservation error:", conservation_error
  call require(conservation_error <= 3.0e-12_dp, &
    "fine-patch removal conservation")

  allocate(lifecycle_integrals(size(final_integrals)))
  coarse_temperature = 1000.0_dp
  coarse_temperature(9, 8) = 2000.0_dp
  call regrid_reactive_eb_amr_hierarchy_2d( &
    species, config, coarse_state, coarse_temperature, coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, patch, fine_active, &
    changed, ok)
  call require(ok .and. changed .and. fine_active .and. &
    allocated(fine_state) .and. allocated(fine_temperature), &
    "fine-patch re-creation from root-only state")
  call composite_eb_integral_2d( &
    coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
    lifecycle_integrals, ok)
  call require(ok .and. maxval(abs(lifecycle_integrals - &
    final_integrals)) <= 3.0e-12_dp * scale, &
    "fine-patch creation conservation")

  coarse_temperature = 1000.0_dp
  call regrid_reactive_eb_amr_hierarchy_2d( &
    species, config, coarse_state, coarse_temperature, coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, patch, fine_active, &
    changed, ok)
  call require(ok .and. changed .and. .not. fine_active, &
    "re-created fine patch collapses on empty tags")
  call reactive_eb_integrals_2d( &
    coarse_state, coarse_geometry, lifecycle_integrals, ok)
  lifecycle_integrals = lifecycle_integrals * &
    coarse_geometry%dx * coarse_geometry%dy
  call require(ok .and. maxval(abs(lifecycle_integrals - &
    final_integrals)) <= 3.0e-12_dp * scale, &
    "re-created patch collapse conservation")

  config%eb%flow%chemistry_enabled = .true.
  call simulate_reactive_eb_amr_2d( &
    species, config, coarse_state, coarse_temperature, coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, patch, fine_active, time, &
    steps, regrids, initial_integrals, final_integrals, minimum_dt, &
    base_density, ok)
  call require(.not. ok .and. steps == 0 .and. regrids == 0 .and. &
    time == 0.0_dp, &
    "unsupported AMR chemistry rejection")

  write(*, '(a)') "test_reactive_eb_amr_2d_driver: PASS"

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_reactive_eb_amr_2d_driver
