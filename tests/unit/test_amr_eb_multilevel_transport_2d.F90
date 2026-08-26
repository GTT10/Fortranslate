program test_amr_eb_multilevel_transport_2d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
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
  use amr_eb_patch_tree_2d_mod, only: &
    amr_eb_patch_tree_level_plan_2d, amr_eb_patch_tree_topology_2d, &
    initialize_amr_eb_patch_tree_topology_2d
  use amr_eb_patch_tree_reactive_2d_mod, only: &
    reactive_amr_eb_patch_tree_2d, &
    initialize_reactive_amr_eb_patch_tree_2d, &
    compute_reactive_amr_eb_patch_tree_cfl_timestep_2d, &
    compute_reactive_amr_eb_patch_tree_timestep_2d, &
    advance_reactive_amr_eb_patch_tree_transport_2d, &
    advance_reactive_amr_eb_patch_tree_full_physics_2d, &
    advance_reactive_amr_eb_patch_tree_to_time_2d, &
    composite_integral_reactive_amr_eb_patch_tree_2d
  use amr_eb_multilevel_transport_2d_mod, only: &
    advance_three_level_reactive_eb_transport_2d
  use reactive_eb_amr_2d_driver_mod, only: &
    advance_three_level_reactive_eb_strang_2d
  implicit none

  integer, parameter :: root_nx = 8, root_ny = 8
  integer, parameter :: root_i_lower = 2, root_i_upper = 7
  integer, parameter :: root_j_lower = 2, root_j_upper = 7
  integer, parameter :: level_one_nx = 12, level_one_ny = 12
  integer, parameter :: level_one_i_lower = 3, level_one_i_upper = 10
  integer, parameter :: level_one_j_lower = 3, level_one_j_upper = 10
  integer, parameter :: level_two_nx = 16, level_two_ny = 16
  integer, parameter :: branch_a_nx = 8, branch_a_ny = 8
  integer, parameter :: branch_b_nx = 4, branch_b_ny = 4
  integer, parameter :: branch_deep_nx = 8, branch_deep_ny = 8
  integer, parameter :: ratio = 2
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(reactive_eb_2d_config) :: config
  type(reactive_boundary_set_2d) :: boundaries
  type(eb_geometry_2d) :: root_geometry, level_one_geometry
  type(eb_geometry_2d) :: level_two_geometry
  type(eb_geometry_2d) :: branch_a_geometry, branch_b_geometry
  type(eb_geometry_2d) :: branch_deep_geometry
  type(amr_eb_patch_2d) :: root_patch, level_one_patch
  type(amr_eb_patch_tree_level_plan_2d), allocatable :: tree_plans(:)
  type(amr_eb_patch_tree_level_plan_2d), allocatable :: branch_plans(:)
  type(amr_eb_patch_tree_topology_2d) :: tree_topology
  type(amr_eb_patch_tree_topology_2d) :: branch_topology
  type(reactive_amr_eb_patch_tree_2d) :: tree, tree_snapshot
  type(reactive_amr_eb_patch_tree_2d) :: branch_tree, branch_snapshot
  type(reactive_amr_eb_patch_tree_2d) :: full_tree, full_snapshot
  type(reactive_amr_eb_patch_tree_2d) :: branch_full_tree
  type(reactive_amr_eb_patch_tree_2d) :: time_tree, time_reference_tree
  type(reactive_amr_eb_patch_tree_2d) :: limited_tree, limited_snapshot
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
  real(dp), allocatable :: branch_integral_before(:), branch_integral_after(:)
  real(dp) :: root_dx, root_dy, level_one_dx, level_one_dy
  real(dp) :: root_x_lower, root_x_upper, root_y_lower, root_y_upper
  real(dp) :: one_x_lower, one_x_upper, one_y_lower, one_y_upper
  real(dp) :: two_x_lower, two_x_upper, two_y_lower, two_y_upper
  real(dp) :: branch_x_lower, branch_x_upper
  real(dp) :: branch_y_lower, branch_y_upper
  real(dp) :: branch_dx, branch_dy, branch_interval
  real(dp) :: full_interval
  real(dp) :: hydro_dt, time_dt, time_value, time_final
  real(dp) :: time_minimum_dt, time_reference, reference_minimum_dt
  real(dp) :: reference_theta, step_dt, step_theta, time_tolerance
  real(dp) :: base_density, root_dt, one_dt, two_dt, diffusivity
  real(dp) :: interval, initial_span, final_span, minimum_theta, scale
  real(dp) :: tree_minimum_theta
  logical :: ok
  character(len=128) :: failure_context
  integer, allocatable :: level_advances(:)
  integer, allocatable :: chemistry_advances(:), transport_advances(:)
  integer, allocatable :: hydro_advances(:)
  integer :: time_steps, reference_steps, advanced_steps

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamic database load")
  call load_h2o2_elementary_mechanism(reactions, ok)
  call require(ok, "chemistry mechanism load")
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

  branch_x_lower = one_x_lower + 2.0_dp * level_one_dx
  branch_x_upper = one_x_lower + 6.0_dp * level_one_dx
  branch_y_lower = one_y_lower + 2.0_dp * level_one_dy
  branch_y_upper = one_y_lower + 6.0_dp * level_one_dy
  call build_configured_eb_geometry_region_2d( &
    config, branch_a_nx, branch_a_ny, branch_x_lower, branch_x_upper, &
    branch_y_lower, branch_y_upper, branch_a_geometry, ok)
  call require(ok, "branch-a transport geometry")
  call build_configured_eb_geometry_region_2d( &
    config, branch_b_nx, branch_b_ny, &
    one_x_lower + 8.0_dp * level_one_dx, &
    one_x_lower + 10.0_dp * level_one_dx, &
    one_y_lower + 8.0_dp * level_one_dy, &
    one_y_lower + 10.0_dp * level_one_dy, branch_b_geometry, ok)
  call require(ok, "branch-b transport geometry")
  branch_dx = (branch_x_upper - branch_x_lower) / real(branch_a_nx, dp)
  branch_dy = (branch_y_upper - branch_y_lower) / real(branch_a_ny, dp)
  call build_configured_eb_geometry_region_2d( &
    config, branch_deep_nx, branch_deep_ny, &
    branch_x_lower + 2.0_dp * branch_dx, &
    branch_x_lower + 6.0_dp * branch_dx, &
    branch_y_lower + 2.0_dp * branch_dy, &
    branch_y_lower + 6.0_dp * branch_dy, branch_deep_geometry, ok)
  call require(ok, "deep transport geometry")

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

  allocate(branch_plans(3))
  branch_plans%refinement_ratio = ratio
  allocate(branch_plans(1)%children(1))
  branch_plans(1)%children(1)%parent_patch = 1
  branch_plans(1)%children(1)%coarse_i_lower = root_i_lower
  branch_plans(1)%children(1)%coarse_i_upper = root_i_upper
  branch_plans(1)%children(1)%coarse_j_lower = root_j_lower
  branch_plans(1)%children(1)%coarse_j_upper = root_j_upper
  branch_plans(1)%children(1)%geometry = level_one_geometry
  allocate(branch_plans(2)%children(2))
  branch_plans(2)%children(1)%parent_patch = 1
  branch_plans(2)%children(1)%coarse_i_lower = 3
  branch_plans(2)%children(1)%coarse_i_upper = 6
  branch_plans(2)%children(1)%coarse_j_lower = 3
  branch_plans(2)%children(1)%coarse_j_upper = 6
  branch_plans(2)%children(1)%geometry = branch_a_geometry
  branch_plans(2)%children(2)%parent_patch = 1
  branch_plans(2)%children(2)%coarse_i_lower = 9
  branch_plans(2)%children(2)%coarse_i_upper = 10
  branch_plans(2)%children(2)%coarse_j_lower = 9
  branch_plans(2)%children(2)%coarse_j_upper = 10
  branch_plans(2)%children(2)%geometry = branch_b_geometry
  allocate(branch_plans(3)%children(1))
  branch_plans(3)%children(1)%parent_patch = 1
  branch_plans(3)%children(1)%coarse_i_lower = 3
  branch_plans(3)%children(1)%coarse_i_upper = 6
  branch_plans(3)%children(1)%coarse_j_lower = 3
  branch_plans(3)%children(1)%coarse_j_upper = 6
  branch_plans(3)%children(1)%geometry = branch_deep_geometry
  call initialize_amr_eb_patch_tree_topology_2d( &
    root_geometry, branch_plans, branch_topology, ok)
  call require(ok, "four-level branching transport topology")
  call initialize_reactive_amr_eb_patch_tree_2d( &
    species, root_state, root_temperature, branch_topology, branch_tree, ok)
  call require(ok, "four-level branching transport state")

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

  allocate(tree_plans(2))
  tree_plans%refinement_ratio = ratio
  allocate(tree_plans(1)%children(1))
  tree_plans(1)%children(1)%parent_patch = 1
  tree_plans(1)%children(1)%coarse_i_lower = root_i_lower
  tree_plans(1)%children(1)%coarse_i_upper = root_i_upper
  tree_plans(1)%children(1)%coarse_j_lower = root_j_lower
  tree_plans(1)%children(1)%coarse_j_upper = root_j_upper
  tree_plans(1)%children(1)%geometry = level_one_geometry
  allocate(tree_plans(2)%children(1))
  tree_plans(2)%children(1)%parent_patch = 1
  tree_plans(2)%children(1)%coarse_i_lower = level_one_i_lower
  tree_plans(2)%children(1)%coarse_i_upper = level_one_i_upper
  tree_plans(2)%children(1)%coarse_j_lower = level_one_j_lower
  tree_plans(2)%children(1)%coarse_j_upper = level_one_j_upper
  tree_plans(2)%children(1)%geometry = level_two_geometry
  call initialize_amr_eb_patch_tree_topology_2d( &
    root_geometry, tree_plans, tree_topology, ok)
  call require(ok, "three-level transport tree topology")
  call initialize_reactive_amr_eb_patch_tree_2d( &
    species, root_state, root_temperature, tree_topology, tree, ok)
  call require(ok, "three-level transport tree state")
  tree%levels(2)%patches(1)%state = level_one_state
  tree%levels(2)%patches(1)%temperature = level_one_temperature
  tree%levels(3)%patches(1)%state = level_two_state
  tree%levels(3)%patches(1)%temperature = level_two_temperature
  allocate(level_advances(tree%level_count()))
  call advance_reactive_amr_eb_patch_tree_transport_2d( &
    species, transport, tree, interval, .false., .true., .false., .false., &
    boundaries, config%state_redist_target_volume_fraction, &
    config%state_redist_max_order, tree_minimum_theta, ok, failure_context, &
    level_advances)
  call require(ok .and. all(level_advances == [2, 4, 8]) .and. &
    tree_minimum_theta > 0.999999999_dp, &
    "three-level patch-tree transport schedule: " // trim(failure_context))

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
  scale = max(1.0_dp, maxval(abs(new_root_state)), &
    maxval(abs(new_level_one_state)), maxval(abs(new_level_two_state)))
  call require(maxval(abs(tree%levels(1)%patches(1)%state - &
      new_root_state)) <= 5.0e-7_dp * scale .and. &
    maxval(abs(tree%levels(2)%patches(1)%state - &
      new_level_one_state)) <= 5.0e-7_dp * scale .and. &
    maxval(abs(tree%levels(3)%patches(1)%state - &
      new_level_two_state)) <= 5.0e-7_dp * scale, &
    "three-level patch-tree transport field parity")
  scale = max(1.0_dp, maxval(new_root_temperature), &
    maxval(new_level_one_temperature), maxval(new_level_two_temperature))
  call require(maxval(abs(tree%levels(1)%patches(1)%temperature - &
      new_root_temperature)) <= 5.0e-7_dp * scale .and. &
    maxval(abs(tree%levels(2)%patches(1)%temperature - &
      new_level_one_temperature)) <= 5.0e-7_dp * scale .and. &
    maxval(abs(tree%levels(3)%patches(1)%temperature - &
      new_level_two_temperature)) <= 5.0e-7_dp * scale, &
    "three-level patch-tree transport temperature parity")
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

  tree_snapshot = tree
  call advance_reactive_amr_eb_patch_tree_transport_2d( &
    species, transport, tree, -interval, .false., .true., .false., .false., &
    boundaries, config%state_redist_target_volume_fraction, &
    config%state_redist_max_order, tree_minimum_theta, ok, failure_context, &
    level_advances)
  call require(.not. ok .and. all(level_advances == 0) .and. &
    tree_minimum_theta == 1.0_dp .and. tree_solutions_match(tree, tree_snapshot), &
    "three-level patch-tree transport rollback")

  allocate(branch_integral_before(size(root_state, 1)))
  allocate(branch_integral_after(size(root_state, 1)))
  call composite_integral_reactive_amr_eb_patch_tree_2d( &
    branch_tree, branch_integral_before, ok)
  call require(ok, "branching transport initial composite integral")
  branch_snapshot = branch_tree
  if (allocated(level_advances)) deallocate(level_advances)
  allocate(level_advances(branch_tree%level_count()))
  branch_interval = min(interval, 1.0e-8_dp)
  call advance_reactive_amr_eb_patch_tree_transport_2d( &
    species, transport, branch_tree, branch_interval, .false., .true., &
    .false., .false., boundaries, &
    config%state_redist_target_volume_fraction, &
    config%state_redist_max_order, tree_minimum_theta, ok, failure_context, &
    level_advances)
  call require(ok .and. all(level_advances == [2, 4, 16, 16]) .and. &
    tree_minimum_theta > 0.0_dp .and. &
    .not. tree_solutions_match(branch_tree, branch_snapshot), &
    "four-level branching patch-tree transport schedule: " // &
      trim(failure_context))
  call composite_integral_reactive_amr_eb_patch_tree_2d( &
    branch_tree, branch_integral_after, ok)
  scale = max(1.0_dp, maxval(abs(branch_integral_before)))
  call require(ok .and. maxval(abs( &
      branch_integral_after - branch_integral_before)) <= &
      3.0e-10_dp * scale .and. branch_tree%is_valid(), &
    "four-level branching patch-tree transport conservation")

  branch_snapshot = branch_tree
  call advance_reactive_amr_eb_patch_tree_transport_2d( &
    species, transport, branch_tree, -branch_interval, .false., .true., &
    .false., .false., boundaries, &
    config%state_redist_target_volume_fraction, &
    config%state_redist_max_order, tree_minimum_theta, ok, failure_context, &
    level_advances)
  call require(.not. ok .and. all(level_advances == 0) .and. &
    tree_minimum_theta == 1.0_dp .and. &
    tree_solutions_match(branch_tree, branch_snapshot), &
    "four-level branching patch-tree transport rollback")

  call initialize_reactive_amr_eb_patch_tree_2d( &
    species, root_state, root_temperature, tree_topology, full_tree, ok)
  call require(ok, "three-level full-physics tree state")
  full_tree%levels(2)%patches(1)%state = level_one_state
  full_tree%levels(2)%patches(1)%temperature = level_one_temperature
  full_tree%levels(3)%patches(1)%state = level_two_state
  full_tree%levels(3)%patches(1)%temperature = level_two_temperature
  if (allocated(chemistry_advances)) deallocate(chemistry_advances)
  if (allocated(transport_advances)) deallocate(transport_advances)
  if (allocated(hydro_advances)) deallocate(hydro_advances)
  allocate(chemistry_advances(full_tree%level_count()))
  allocate(transport_advances(full_tree%level_count()))
  allocate(hydro_advances(full_tree%level_count()))
  full_interval = min(interval, 1.0e-8_dp)
  call advance_reactive_amr_eb_patch_tree_full_physics_2d( &
    species, reactions, transport, full_tree, "hllc", "pcm", "mc", &
    config%state_redist_max_order, full_interval, .true., 1.0e-7_dp, &
    1.0e-13_dp, .false., .true., .false., .false., boundaries, &
    config%state_redist_target_volume_fraction, tree_minimum_theta, ok, &
    failure_context, chemistry_advances, transport_advances, hydro_advances)
  call require(ok .and. all(chemistry_advances == [2, 2, 2]) .and. &
    all(transport_advances == [4, 8, 16]) .and. &
    all(hydro_advances == [1, 2, 4]) .and. tree_minimum_theta > 0.0_dp, &
    "three-level patch-tree full-physics schedule: " // &
      trim(failure_context))
  call advance_three_level_reactive_eb_strang_2d( &
    species, reactions, root_state, root_temperature, root_geometry, &
    level_one_state, level_one_temperature, level_one_geometry, root_patch, &
    level_two_state, level_two_temperature, level_two_geometry, &
    level_one_patch, "hllc", "pcm", "mc", &
    config%state_redist_max_order, full_interval, .true., 1.0e-7_dp, &
    1.0e-13_dp, new_root_state, new_root_temperature, new_level_one_state, &
    new_level_one_temperature, new_level_two_state, &
    new_level_two_temperature, ok, &
    config%state_redist_target_volume_fraction, transport, .true., .false., &
    .true., .false., .false., minimum_theta, boundaries)
  call require(ok .and. minimum_theta > 0.0_dp, &
    "fixed three-level full-physics transaction")
  scale = max(1.0_dp, maxval(abs(new_root_state)), &
    maxval(abs(new_level_one_state)), maxval(abs(new_level_two_state)))
  call require(maxval(abs(full_tree%levels(1)%patches(1)%state - &
      new_root_state)) <= 5.0e-7_dp * scale .and. &
    maxval(abs(full_tree%levels(2)%patches(1)%state - &
      new_level_one_state)) <= 5.0e-7_dp * scale .and. &
    maxval(abs(full_tree%levels(3)%patches(1)%state - &
      new_level_two_state)) <= 5.0e-7_dp * scale, &
    "three-level patch-tree full-physics field parity")
  scale = max(1.0_dp, maxval(new_root_temperature), &
    maxval(new_level_one_temperature), maxval(new_level_two_temperature))
  call require(maxval(abs(full_tree%levels(1)%patches(1)%temperature - &
      new_root_temperature)) <= 5.0e-7_dp * scale .and. &
    maxval(abs(full_tree%levels(2)%patches(1)%temperature - &
      new_level_one_temperature)) <= 5.0e-7_dp * scale .and. &
    maxval(abs(full_tree%levels(3)%patches(1)%temperature - &
      new_level_two_temperature)) <= 5.0e-7_dp * scale, &
    "three-level patch-tree full-physics temperature parity")

  full_snapshot = full_tree
  call advance_reactive_amr_eb_patch_tree_full_physics_2d( &
    species, reactions, transport, full_tree, "unknown", "pcm", "mc", &
    config%state_redist_max_order, full_interval, .true., 1.0e-7_dp, &
    1.0e-13_dp, .false., .true., .false., .false., boundaries, &
    config%state_redist_target_volume_fraction, tree_minimum_theta, ok, &
    failure_context, chemistry_advances, transport_advances, hydro_advances)
  call require(.not. ok .and. all(chemistry_advances == 0) .and. &
    all(transport_advances == 0) .and. all(hydro_advances == 0) .and. &
    tree_minimum_theta == 1.0_dp .and. &
    tree_solutions_match(full_tree, full_snapshot), &
    "three-level patch-tree full-physics rollback")

  call initialize_reactive_amr_eb_patch_tree_2d( &
    species, root_state, root_temperature, branch_topology, &
    branch_full_tree, ok)
  call require(ok, "four-level branching full-physics tree state")
  deallocate(chemistry_advances, transport_advances, hydro_advances)
  allocate(chemistry_advances(branch_full_tree%level_count()))
  allocate(transport_advances(branch_full_tree%level_count()))
  allocate(hydro_advances(branch_full_tree%level_count()))
  call composite_integral_reactive_amr_eb_patch_tree_2d( &
    branch_full_tree, branch_integral_before, ok)
  call require(ok, "branching full-physics initial integral")
  call advance_reactive_amr_eb_patch_tree_full_physics_2d( &
    species, reactions, transport, branch_full_tree, "hllc", "pcm", "mc", &
    config%state_redist_max_order, full_interval, .true., 1.0e-7_dp, &
    1.0e-13_dp, .false., .true., .false., .false., boundaries, &
    config%state_redist_target_volume_fraction, tree_minimum_theta, ok, &
    failure_context, chemistry_advances, transport_advances, hydro_advances)
  call require(ok .and. all(chemistry_advances == [2, 2, 4, 2]) .and. &
    all(transport_advances == [4, 8, 32, 32]) .and. &
    all(hydro_advances == [1, 2, 8, 8]) .and. &
    tree_minimum_theta > 0.0_dp .and. branch_full_tree%is_valid(), &
    "four-level branching patch-tree full-physics schedule: " // &
      trim(failure_context))
  call composite_integral_reactive_amr_eb_patch_tree_2d( &
    branch_full_tree, branch_integral_after, ok)
  scale = max(1.0_dp, maxval(abs(branch_integral_before)))
  call require(ok .and. maxval(abs( &
      branch_integral_after - branch_integral_before)) <= &
      3.0e-8_dp * scale, &
    "four-level branching patch-tree full-physics conservation")

  call initialize_reactive_amr_eb_patch_tree_2d( &
    species, root_state, root_temperature, tree_topology, &
    time_reference_tree, ok)
  call require(ok, "three-level time-loop reference state")
  time_reference_tree%levels(2)%patches(1)%state = level_one_state
  time_reference_tree%levels(2)%patches(1)%temperature = &
    level_one_temperature
  time_reference_tree%levels(3)%patches(1)%state = level_two_state
  time_reference_tree%levels(3)%patches(1)%temperature = &
    level_two_temperature
  time_tree = time_reference_tree
  limited_snapshot = time_reference_tree
  call compute_reactive_amr_eb_patch_tree_cfl_timestep_2d( &
    species, time_reference_tree, 0.004_dp, hydro_dt, ok)
  call require(ok, "three-level hydro CFL timestep")
  call compute_reactive_amr_eb_patch_tree_timestep_2d( &
    species, transport, time_reference_tree, 0.004_dp, 0.30_dp, &
    .false., .true., .false., time_dt, ok)
  scale = max(1.0_dp, abs(hydro_dt), abs(root_dt), &
    abs(real(ratio, dp) * one_dt), &
    abs(real(ratio * ratio, dp) * two_dt))
  call require(ok .and. abs(time_dt - min(hydro_dt, root_dt, &
      real(ratio, dp) * one_dt, &
      real(ratio * ratio, dp) * two_dt)) <= &
      128.0_dp * epsilon(1.0_dp) * scale, &
    "three-level combined hydro-transport timestep")

  if (allocated(chemistry_advances)) deallocate(chemistry_advances)
  if (allocated(transport_advances)) deallocate(transport_advances)
  if (allocated(hydro_advances)) deallocate(hydro_advances)
  allocate(chemistry_advances(time_reference_tree%level_count()))
  allocate(transport_advances(time_reference_tree%level_count()))
  allocate(hydro_advances(time_reference_tree%level_count()))
  time_final = 1.25_dp * time_dt
  time_reference = 0.0_dp
  reference_steps = 0
  reference_minimum_dt = 0.0_dp
  reference_theta = 1.0_dp
  time_tolerance = 16.0_dp * epsilon(1.0_dp) * &
    max(tiny(1.0_dp), abs(time_final))
  do while (time_final - time_reference > time_tolerance)
    call compute_reactive_amr_eb_patch_tree_timestep_2d( &
      species, transport, time_reference_tree, 0.004_dp, 0.30_dp, &
      .false., .true., .false., step_dt, ok)
    call require(ok, "manual time-loop timestep")
    step_dt = min(step_dt, time_final - time_reference)
    call advance_reactive_amr_eb_patch_tree_full_physics_2d( &
      species, reactions, transport, time_reference_tree, "hllc", "pcm", &
      "mc", config%state_redist_max_order, step_dt, .true., 1.0e-7_dp, &
      1.0e-13_dp, .false., .true., .false., .false., boundaries, &
      config%state_redist_target_volume_fraction, step_theta, ok, &
      failure_context, chemistry_advances, transport_advances, &
      hydro_advances)
    call require(ok, "manual full-physics time-loop step: " // &
      trim(failure_context))
    time_reference = time_reference + step_dt
    reference_steps = reference_steps + 1
    if (reference_minimum_dt == 0.0_dp) then
      reference_minimum_dt = step_dt
    else
      reference_minimum_dt = min(reference_minimum_dt, step_dt)
    end if
    reference_theta = min(reference_theta, step_theta)
  end do
  time_reference = time_final
  call require(reference_steps == 2, "manual two-step reference schedule")

  time_value = 0.0_dp
  time_steps = 0
  call advance_reactive_amr_eb_patch_tree_to_time_2d( &
    species, reactions, transport, time_tree, "hllc", "pcm", "mc", &
    config%state_redist_max_order, time_value, time_final, time_steps, 8, &
    0.004_dp, 0.30_dp, .true., 1.0e-7_dp, 1.0e-13_dp, .false., .true., &
    .false., .false., boundaries, &
    config%state_redist_target_volume_fraction, time_minimum_dt, &
    tree_minimum_theta, ok, failure_context, advanced_steps, &
    chemistry_advances, transport_advances, hydro_advances)
  call require(ok .and. time_value == time_final .and. &
    time_steps == reference_steps .and. advanced_steps == reference_steps &
    .and. all(chemistry_advances == reference_steps * [2, 2, 2]) .and. &
    all(transport_advances == reference_steps * [4, 8, 16]) .and. &
    all(hydro_advances == reference_steps * [1, 2, 4]) .and. &
    abs(time_minimum_dt - reference_minimum_dt) <= &
      128.0_dp * epsilon(1.0_dp) * &
        max(1.0_dp, abs(reference_minimum_dt)) .and. &
    tree_minimum_theta == reference_theta .and. &
    tree_solutions_match(time_tree, time_reference_tree), &
    "three-level adaptive time-loop parity: " // trim(failure_context))

  limited_tree = limited_snapshot
  time_value = 0.0_dp
  time_steps = 0
  call advance_reactive_amr_eb_patch_tree_to_time_2d( &
    species, reactions, transport, limited_tree, "hllc", "pcm", "mc", &
    config%state_redist_max_order, time_value, time_final, time_steps, 1, &
    0.004_dp, 0.30_dp, .true., 1.0e-7_dp, 1.0e-13_dp, .false., .true., &
    .false., .false., boundaries, &
    config%state_redist_target_volume_fraction, time_minimum_dt, &
    tree_minimum_theta, ok, failure_context, advanced_steps, &
    chemistry_advances, transport_advances, hydro_advances)
  call require(.not. ok .and. trim(failure_context) == "maximum steps" &
    .and. time_steps == 1 .and. advanced_steps == 1 .and. &
    time_value > 0.0_dp .and. time_value < time_final .and. &
    time_minimum_dt == time_value .and. &
    all(chemistry_advances == [2, 2, 2]) .and. &
    all(transport_advances == [4, 8, 16]) .and. &
    all(hydro_advances == [1, 2, 4]) .and. &
    .not. tree_solutions_match(limited_tree, limited_snapshot), &
    "time-loop maximum-step committed prefix")

  limited_tree = limited_snapshot
  time_value = 0.0_dp
  time_steps = 0
  call advance_reactive_amr_eb_patch_tree_to_time_2d( &
    species, reactions, transport, limited_tree, "unknown", "pcm", "mc", &
    config%state_redist_max_order, time_value, time_final, time_steps, 8, &
    0.004_dp, 0.30_dp, .true., 1.0e-7_dp, 1.0e-13_dp, .false., .true., &
    .false., .false., boundaries, &
    config%state_redist_target_volume_fraction, time_minimum_dt, &
    tree_minimum_theta, ok, failure_context, advanced_steps, &
    chemistry_advances, transport_advances, hydro_advances)
  call require(.not. ok .and. time_value == 0.0_dp .and. time_steps == 0 &
    .and. advanced_steps == 0 .and. time_minimum_dt == 0.0_dp .and. &
    tree_minimum_theta == 1.0_dp .and. all(chemistry_advances == 0) .and. &
    all(transport_advances == 0) .and. all(hydro_advances == 0) .and. &
    tree_solutions_match(limited_tree, limited_snapshot), &
    "time-loop failed first-step rollback")

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

  logical function tree_solutions_match(first, second) result(match)
    type(reactive_amr_eb_patch_tree_2d), intent(in) :: first, second
    integer :: level, patch

    match = first%level_count() == second%level_count()
    if (.not. match) return
    do level = 1, first%level_count()
      match = first%levels(level)%patch_count() == &
        second%levels(level)%patch_count()
      if (.not. match) return
      do patch = 1, first%levels(level)%patch_count()
        match = all(first%levels(level)%patches(patch)%state == &
            second%levels(level)%patches(patch)%state) .and. &
          all(first%levels(level)%patches(patch)%temperature == &
            second%levels(level)%patches(patch)%temperature)
        if (.not. match) return
      end do
    end do
  end function tree_solutions_match

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) error stop message
  end subroutine require

end program test_amr_eb_multilevel_transport_2d
