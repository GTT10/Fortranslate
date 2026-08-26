program test_amr_eb_multilevel_2d
  use, intrinsic :: ieee_arithmetic, only: &
    ieee_is_finite, ieee_value, ieee_quiet_nan
  use precision_mod, only: dp
  use state_indices_mod, only: irho, iet
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_mass_fraction_component, reactive_primitive_to_conserved
  use eb_geometry_2d_mod, only: &
    eb_geometry_2d, eb_covered_cell, build_eb_geometry_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, build_amr_eb_patch_2d
  use amr_eb_patch_tree_2d_mod, only: &
    amr_eb_patch_tree_child_plan_2d, &
    amr_eb_patch_tree_level_plan_2d, &
    amr_eb_patch_tree_topology_2d, &
    initialize_amr_eb_patch_tree_topology_2d, &
    rebuild_amr_eb_patch_tree_topology_2d
  use amr_eb_patch_tree_reactive_2d_mod, only: &
    reactive_amr_eb_patch_tree_2d, &
    initialize_reactive_amr_eb_patch_tree_2d, &
    synchronize_reactive_amr_eb_patch_tree_2d, &
    rebuild_reactive_amr_eb_patch_tree_2d, &
    compute_reactive_amr_eb_patch_tree_cfl_timestep_2d, &
    advance_reactive_amr_eb_patch_tree_hydro_2d, &
    composite_integral_reactive_amr_eb_patch_tree_2d
  use amr_eb_multilevel_2d_mod, only: &
    average_down_three_level_eb_state_2d, &
    average_down_three_level_reactive_eb_state_2d, &
    composite_three_level_eb_integral_2d
  use amr_eb_multilevel_reactive_2d_mod, only: &
    advance_three_level_reactive_eb_hydro_2d
  use reactive_eb_amr_2d_driver_mod, only: &
    advance_three_level_reactive_eb_strang_2d
  use reactive_eb_2d_driver_mod, only: &
    compute_reactive_eb_cfl_timestep_2d
  implicit none

  integer, parameter :: root_nx = 8, root_ny = 8, ratio = 2
  integer, parameter :: root_i_lower = 2, root_i_upper = 7
  integer, parameter :: root_j_lower = 2, root_j_upper = 7
  integer, parameter :: level_one_nx = &
    (root_i_upper - root_i_lower + 1) * ratio
  integer, parameter :: level_one_ny = &
    (root_j_upper - root_j_lower + 1) * ratio
  integer, parameter :: level_one_i_lower = 3, level_one_i_upper = 10
  integer, parameter :: level_one_j_lower = 3, level_one_j_upper = 10
  integer, parameter :: level_two_nx = &
    (level_one_i_upper - level_one_i_lower + 1) * ratio
  integer, parameter :: level_two_ny = &
    (level_one_j_upper - level_one_j_lower + 1) * ratio
  type(eb_geometry_2d) :: root_geometry
  type(eb_geometry_2d) :: level_one_geometry, level_two_geometry
  type(eb_geometry_2d) :: tree_level_one_geometry_a
  type(eb_geometry_2d) :: tree_level_one_geometry_b
  type(eb_geometry_2d) :: tree_level_two_geometry_a
  type(eb_geometry_2d) :: tree_level_two_geometry_b
  type(eb_geometry_2d) :: tree_level_three_geometry
  type(eb_geometry_2d) :: tree_level_three_shifted_geometry
  type(amr_eb_patch_2d) :: root_patch, level_one_patch
  type(amr_eb_patch_2d) :: tree_scratch_patch
  type(amr_eb_patch_tree_level_plan_2d), allocatable :: tree_plans(:)
  type(amr_eb_patch_tree_level_plan_2d), allocatable :: extended_tree_plans(:)
  type(amr_eb_patch_tree_level_plan_2d), allocatable :: shifted_tree_plans(:)
  type(amr_eb_patch_tree_level_plan_2d), allocatable :: invalid_tree_plans(:)
  type(amr_eb_patch_tree_level_plan_2d), allocatable :: chain_tree_plans(:)
  type(amr_eb_patch_tree_topology_2d) :: tree_topology
  type(amr_eb_patch_tree_topology_2d) :: chain_tree_topology
  type(reactive_amr_eb_patch_tree_2d) :: reactive_tree
  type(reactive_amr_eb_patch_tree_2d) :: reactive_tree_snapshot
  type(reactive_amr_eb_patch_tree_2d) :: chain_tree
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp) :: root_level_set(0:root_nx, 0:root_ny)
  real(dp), allocatable :: root_state(:, :, :)
  real(dp), allocatable :: level_one_state(:, :, :)
  real(dp), allocatable :: level_two_state(:, :, :)
  real(dp), allocatable :: synchronized_root(:, :, :)
  real(dp), allocatable :: synchronized_level_one(:, :, :)
  real(dp), allocatable :: integral_before(:), integral_after(:)
  real(dp), allocatable :: primitive(:), mass_fractions(:), state_cell(:)
  real(dp), allocatable :: tree_integral_before(:), tree_integral_after(:)
  real(dp), allocatable :: tree_level_two_saved(:, :, :)
  real(dp), allocatable :: tree_deepest_saved(:, :, :)
  real(dp), allocatable :: reactive_root(:, :, :)
  real(dp), allocatable :: reactive_level_one(:, :, :)
  real(dp), allocatable :: reactive_level_two(:, :, :)
  real(dp), allocatable :: reactive_root_sync(:, :, :)
  real(dp), allocatable :: reactive_level_one_sync(:, :, :)
  real(dp), allocatable :: reactive_level_two_sync(:, :, :)
  real(dp), allocatable :: root_temperature(:, :)
  real(dp), allocatable :: level_one_temperature(:, :)
  real(dp), allocatable :: level_two_temperature(:, :)
  real(dp), allocatable :: root_temperature_sync(:, :)
  real(dp), allocatable :: level_one_temperature_sync(:, :)
  real(dp), allocatable :: level_two_temperature_sync(:, :)
  real(dp) :: mole_fractions(7), x, y, temperature_cell, sound_speed
  real(dp) :: scale, dt, species_integral_sum, species_change
  real(dp) :: tree_dt, reference_tree_dt, initial_tree_dt, node_dt
  logical :: ok, topology_changed, reference_ok, node_ok
  character(len=128) :: tree_failure_context
  integer, allocatable :: tree_level_advances(:), chain_level_advances(:)
  integer :: i, j, k, level, patch, nvar

  do j = 0, root_ny
    y = real(j, dp) / real(root_ny, dp)
    do i = 0, root_nx
      x = real(i, dp) / real(root_nx, dp)
      root_level_set(i, j) = x + y - 0.78_dp
    end do
  end do
  call build_eb_geometry_2d( &
    root_level_set, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, root_geometry, ok)
  call require(ok, "three-level root geometry")
  call build_patch_geometry( &
    root_geometry, root_i_lower, root_i_upper, root_j_lower, root_j_upper, &
    ratio, level_one_geometry, root_patch, ok)
  call require(ok, "three-level middle geometry")
  call build_patch_geometry( &
    level_one_geometry, level_one_i_lower, level_one_i_upper, &
    level_one_j_lower, level_one_j_upper, ratio, level_two_geometry, &
    level_one_patch, ok)
  call require(ok, "three-level finest geometry")

  call build_patch_geometry( &
    root_geometry, 1, 2, 1, 2, ratio, tree_level_one_geometry_a, &
    tree_scratch_patch, ok)
  call require(ok, "patch-tree first root child geometry")
  call build_patch_geometry( &
    root_geometry, 6, 8, 6, 8, ratio, tree_level_one_geometry_b, &
    tree_scratch_patch, ok)
  call require(ok, "patch-tree second root child geometry")
  call build_patch_geometry( &
    tree_level_one_geometry_a, 1, 4, 1, 4, ratio, &
    tree_level_two_geometry_a, tree_scratch_patch, ok)
  call require(ok, "patch-tree first branch geometry")
  call build_patch_geometry( &
    tree_level_one_geometry_b, 2, 5, 2, 5, ratio, &
    tree_level_two_geometry_b, tree_scratch_patch, ok)
  call require(ok, "patch-tree second branch geometry")
  call build_patch_geometry( &
    tree_level_two_geometry_b, 2, 7, 2, 7, ratio, &
    tree_level_three_geometry, tree_scratch_patch, ok)
  call require(ok, "patch-tree fourth-level geometry")
  call build_patch_geometry( &
    tree_level_two_geometry_b, 3, 8, 2, 7, ratio, &
    tree_level_three_shifted_geometry, tree_scratch_patch, ok)
  call require(ok, "shifted patch-tree fourth-level geometry")

  allocate(tree_plans(2))
  tree_plans%refinement_ratio = ratio
  allocate(tree_plans(1)%children(2), tree_plans(2)%children(2))
  call set_tree_child_plan( &
    tree_plans(1)%children(1), 1, 1, 2, 1, 2, &
    tree_level_one_geometry_a)
  call set_tree_child_plan( &
    tree_plans(1)%children(2), 1, 6, 8, 6, 8, &
    tree_level_one_geometry_b)
  call set_tree_child_plan( &
    tree_plans(2)%children(1), 1, 1, 4, 1, 4, &
    tree_level_two_geometry_a)
  call set_tree_child_plan( &
    tree_plans(2)%children(2), 2, 2, 5, 2, 5, &
    tree_level_two_geometry_b)
  call initialize_amr_eb_patch_tree_topology_2d( &
    root_geometry, tree_plans, tree_topology, ok)
  call require(ok .and. tree_topology%is_valid() .and. &
    tree_topology%level_count() == 3 .and. &
    tree_topology%level_patch_count(0) == 1 .and. &
    tree_topology%level_patch_count(1) == 2 .and. &
    tree_topology%level_patch_count(2) == 2 .and. &
    tree_topology%relations(1)%child_index(1, 2) == 2 .and. &
    tree_topology%relations(2)%child_index(2, 1) == 2, &
    "branching EB patch-tree topology")

  allocate(extended_tree_plans(3))
  extended_tree_plans(1:2) = tree_plans
  extended_tree_plans(3)%refinement_ratio = ratio
  allocate(extended_tree_plans(3)%children(1))
  call set_tree_child_plan( &
    extended_tree_plans(3)%children(1), 2, 2, 7, 2, 7, &
    tree_level_three_geometry)
  shifted_tree_plans = extended_tree_plans
  call set_tree_child_plan( &
    shifted_tree_plans(3)%children(1), 2, 3, 8, 2, 7, &
    tree_level_three_shifted_geometry)
  call rebuild_amr_eb_patch_tree_topology_2d( &
    tree_topology, extended_tree_plans, ok, topology_changed)
  call require(ok .and. topology_changed .and. &
    tree_topology%is_valid() .and. tree_topology%level_count() == 4 .and. &
    tree_topology%level_patch_count(3) == 1 .and. &
    tree_topology%relations(3)%child_index(1, 1) == 0 .and. &
    tree_topology%relations(3)%child_index(2, 1) == 1, &
    "arbitrary-depth EB patch-tree rebuild")
  call rebuild_amr_eb_patch_tree_topology_2d( &
    tree_topology, extended_tree_plans, ok, topology_changed)
  call require(ok .and. .not. topology_changed, &
    "unchanged EB patch-tree rebuild")

  invalid_tree_plans = shifted_tree_plans
  invalid_tree_plans(3)%children(1)%parent_patch = 3
  call rebuild_amr_eb_patch_tree_topology_2d( &
    tree_topology, invalid_tree_plans, ok, topology_changed)
  call require(.not. ok .and. .not. topology_changed .and. &
    tree_topology%is_valid() .and. tree_topology%level_count() == 4 .and. &
    tree_topology%relations(3)%children(1)%parent_patch == 2, &
    "invalid EB patch-tree rebuild rollback")
  call require( &
    any(level_two_geometry%x_face_fraction(0, :) < 1.0_dp) .or. &
    any(level_two_geometry%x_face_fraction(level_two_geometry%nx, :) < &
      1.0_dp) .or. &
    any(level_two_geometry%y_face_fraction(:, 0) < 1.0_dp) .or. &
    any(level_two_geometry%y_face_fraction(:, level_two_geometry%ny) < &
      1.0_dp), "cut finest coarse/fine interface")

  allocate(root_state(1, root_nx, root_ny), source=1.0_dp)
  allocate(level_one_state(1, level_one_nx, level_one_ny), source=2.0_dp)
  allocate(level_two_state(1, level_two_nx, level_two_ny), source=3.0_dp)
  allocate(synchronized_root(1, root_nx, root_ny))
  allocate(synchronized_level_one(1, level_one_nx, level_one_ny))
  allocate(integral_before(1), integral_after(1))
  call composite_three_level_eb_integral_2d( &
    root_state, root_geometry, level_one_state, level_one_geometry, &
    root_patch, level_two_state, level_two_geometry, level_one_patch, &
    integral_before, ok)
  call require(ok, "three-level composite integral")
  call average_down_three_level_eb_state_2d( &
    root_state, root_geometry, level_one_state, level_one_geometry, &
    root_patch, level_two_state, level_two_geometry, level_one_patch, &
    synchronized_root, synchronized_level_one, ok)
  call require(ok, "deepest-to-root EB synchronization")
  integral_after(1) = sum(root_geometry%volume_fraction * &
    synchronized_root(1, :, :)) * root_geometry%dx * root_geometry%dy
  call assert_close(integral_after(1), integral_before(1), 8.0e-13_dp, &
    "three-level synchronization conservation")
  call require(maxval(abs(synchronized_root(:, 1, :) - &
    root_state(:, 1, :))) == 0.0_dp .and. &
    maxval(abs(synchronized_root(:, 8, :) - &
      root_state(:, 8, :))) == 0.0_dp, &
    "root cells outside middle patch unchanged")
  call require(maxval(abs(synchronized_level_one(:, &
    level_one_i_lower:level_one_i_upper, &
    level_one_j_lower:level_one_j_upper) - 3.0_dp)) == 0.0_dp, &
    "finest constant restricted into middle")

  level_two_state(1, level_two_nx, level_two_ny) = &
    ieee_value(0.0_dp, ieee_quiet_nan)
  call average_down_three_level_eb_state_2d( &
    root_state, root_geometry, level_one_state, level_one_geometry, &
    root_patch, level_two_state, level_two_geometry, level_one_patch, &
    synchronized_root, synchronized_level_one, ok)
  call require(.not. ok .and. all(synchronized_root == root_state) .and. &
    all(synchronized_level_one == level_one_state), &
    "three-level nonfinite rollback")
  level_two_state(1, level_two_nx, level_two_ny) = 3.0_dp

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "three-level thermodynamic database")
  call load_h2o2_elementary_mechanism(reactions, ok)
  call require(ok, "three-level elementary mechanism")
  nvar = reactive_nvar(size(species))
  allocate(primitive(reactive_nprim(size(species))))
  allocate(mass_fractions(size(species)), state_cell(nvar))
  mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.55643_dp]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  call require(ok, "three-level composition conversion")
  primitive(1:5) = [0.31_dp, 2.0_dp, -1.0_dp, 0.0_dp, 135000.0_dp]
  do i = 1, size(species)
    primitive(reactive_mass_fraction_component(i)) = mass_fractions(i)
  end do
  call reactive_primitive_to_conserved( &
    species, primitive, state_cell, temperature_cell, sound_speed, ok)
  call require(ok, "three-level reference reactive state")
  allocate(reactive_root(nvar, root_nx, root_ny))
  allocate(reactive_level_one(nvar, level_one_nx, level_one_ny))
  allocate(reactive_level_two(nvar, level_two_nx, level_two_ny))
  allocate(reactive_root_sync(nvar, root_nx, root_ny))
  allocate(reactive_level_one_sync(nvar, level_one_nx, level_one_ny))
  allocate(reactive_level_two_sync(nvar, level_two_nx, level_two_ny))
  allocate(root_temperature(root_nx, root_ny), source=temperature_cell)
  allocate(level_one_temperature( &
    level_one_nx, level_one_ny), source=temperature_cell)
  allocate(level_two_temperature( &
    level_two_nx, level_two_ny), source=temperature_cell)
  allocate(root_temperature_sync(root_nx, root_ny))
  allocate(level_one_temperature_sync(level_one_nx, level_one_ny))
  allocate(level_two_temperature_sync(level_two_nx, level_two_ny))
  reactive_root = spread(spread(state_cell, 2, root_nx), 3, root_ny)
  reactive_level_one = 1.01_dp * &
    spread(spread(state_cell, 2, level_one_nx), 3, level_one_ny)
  reactive_level_two = 0.99_dp * &
    spread(spread(state_cell, 2, level_two_nx), 3, level_two_ny)

  call initialize_amr_eb_patch_tree_topology_2d( &
    root_geometry, tree_plans, tree_topology, ok)
  call require(ok, "reactive patch-tree source topology")
  call initialize_reactive_amr_eb_patch_tree_2d( &
    species, reactive_root, root_temperature, tree_topology, &
    reactive_tree, ok)
  call require(ok .and. reactive_tree%is_valid() .and. &
    reactive_tree%level_count() == 3 .and. &
    reactive_tree%level_patch_count(1) == 2 .and. &
    reactive_tree%level_patch_count(2) == 2, &
    "reactive arbitrary-depth patch-tree initialization")
  reactive_tree%levels(2)%patches(1)%state = &
    1.01_dp * reactive_tree%levels(2)%patches(1)%state
  reactive_tree%levels(2)%patches(2)%state = &
    0.98_dp * reactive_tree%levels(2)%patches(2)%state
  reactive_tree%levels(3)%patches(1)%state = &
    1.02_dp * reactive_tree%levels(3)%patches(1)%state
  reactive_tree%levels(3)%patches(2)%state = &
    0.99_dp * reactive_tree%levels(3)%patches(2)%state
  allocate(tree_integral_before(nvar), tree_integral_after(nvar))
  call composite_integral_reactive_amr_eb_patch_tree_2d( &
    reactive_tree, tree_integral_before, ok)
  call require(ok, "reactive patch-tree composite integral")
  call synchronize_reactive_amr_eb_patch_tree_2d( &
    species, reactive_tree, ok)
  call composite_integral_reactive_amr_eb_patch_tree_2d( &
    reactive_tree, tree_integral_after, ok)
  scale = max(1.0_dp, maxval(abs(tree_integral_before)))
  call require(ok .and. maxval(abs(tree_integral_after - &
    tree_integral_before)) <= 8.0e-12_dp * scale, &
    "deepest-first reactive patch-tree synchronization")
  tree_integral_before = tree_integral_after
  allocate(tree_level_two_saved, source= &
    reactive_tree%levels(3)%patches(2)%state)
  call rebuild_reactive_amr_eb_patch_tree_2d( &
    species, reactive_tree, extended_tree_plans, ok, topology_changed, &
    tree_failure_context)
  call require(ok .and. topology_changed .and. &
    reactive_tree%is_valid() .and. reactive_tree%level_count() == 4 .and. &
    reactive_tree%level_patch_count(3) == 1, &
    "reactive arbitrary-depth patch-tree rebuild: " // &
      trim(tree_failure_context))
  call composite_integral_reactive_amr_eb_patch_tree_2d( &
    reactive_tree, tree_integral_after, ok)
  scale = max(1.0_dp, maxval(abs(tree_integral_before)))
  call require(ok .and. maxval(abs(tree_integral_after - &
    tree_integral_before)) <= 8.0e-12_dp * scale, &
    "reactive patch-tree rebuild conservation")
  call require(maxval(abs(reactive_tree%levels(3)%patches(2)%state - &
    tree_level_two_saved)) <= 8.0e-12_dp * scale, &
    "same-resolution patch-tree overlap retention")
  call require(maxval(abs( &
    reactive_tree%levels(4)%patches(1)%state(:, 1, 1) - &
    reactive_tree%levels(3)%patches(2)%state(:, 2, 2))) <= &
    8.0e-12_dp * scale, "new deepest patch PCM initialization")

  do j = 1, size(reactive_tree%levels(4)%patches(1)%state, 3)
    do i = 1, size(reactive_tree%levels(4)%patches(1)%state, 2)
      reactive_tree%levels(4)%patches(1)%state(:, i, j) = &
        (1.0_dp + 1.0e-4_dp * real(i, dp) + &
          2.0e-4_dp * real(j, dp)) * &
        reactive_tree%levels(4)%patches(1)%state(:, i, j)
    end do
  end do
  call composite_integral_reactive_amr_eb_patch_tree_2d( &
    reactive_tree, tree_integral_before, ok)
  call require(ok, "moving reactive patch-tree source integral")
  allocate(tree_deepest_saved, source= &
    reactive_tree%levels(4)%patches(1)%state)
  call rebuild_reactive_amr_eb_patch_tree_2d( &
    species, reactive_tree, shifted_tree_plans, ok, topology_changed, &
    tree_failure_context)
  call require(ok .and. topology_changed .and. reactive_tree%is_valid(), &
    "moving reactive patch-tree rebuild: " // trim(tree_failure_context))
  call composite_integral_reactive_amr_eb_patch_tree_2d( &
    reactive_tree, tree_integral_after, ok)
  scale = max(1.0_dp, maxval(abs(tree_integral_before)))
  call require(ok .and. maxval(abs(tree_integral_after - &
    tree_integral_before)) <= 8.0e-12_dp * scale, &
    "moving reactive patch-tree conservation")
  call require(maxval(abs( &
    reactive_tree%levels(4)%patches(1)%state(:, 1, 1) - &
    tree_deepest_saved(:, 3, 1))) <= 8.0e-12_dp * scale, &
    "moving reactive patch-tree physical overlap retention")

  reactive_tree_snapshot = reactive_tree
  call compute_reactive_amr_eb_patch_tree_cfl_timestep_2d( &
    species, reactive_tree, 0.4_dp, tree_dt, ok)
  call reference_reactive_tree_cfl_timestep( &
    species, reactive_tree, 0.4_dp, reference_tree_dt, reference_ok)
  call require(ok, "arbitrary-depth patch-tree CFL acceptance")
  call require(reference_ok, "reference patch-tree CFL acceptance")
  call assert_close( &
    tree_dt, reference_tree_dt, 2.0e-12_dp, &
    "arbitrary-depth reactive patch-tree CFL traversal")
  call require( &
    reactive_tree_solutions_match(reactive_tree, reactive_tree_snapshot), &
    "read-only patch-tree CFL traversal")
  initial_tree_dt = tree_dt

  primitive(1:5) = [0.31_dp, 3000.0_dp, -1.0_dp, 0.0_dp, 135000.0_dp]
  call reactive_primitive_to_conserved( &
    species, primitive, state_cell, temperature_cell, sound_speed, ok)
  call require(ok, "deepest patch-tree CFL limiting state")
  reactive_tree%levels(4)%patches(1)%state = spread(spread( &
    state_cell, 2, size(reactive_tree%levels(4)%patches(1)%state, 2)), &
    3, size(reactive_tree%levels(4)%patches(1)%state, 3))
  reactive_tree%levels(4)%patches(1)%temperature = temperature_cell
  call compute_reactive_amr_eb_patch_tree_cfl_timestep_2d( &
    species, reactive_tree, 0.4_dp, tree_dt, ok)
  call reference_reactive_tree_cfl_timestep( &
    species, reactive_tree, 0.4_dp, reference_tree_dt, reference_ok)
  call compute_reactive_eb_cfl_timestep_2d( &
    species, reactive_tree%levels(4)%patches(1)%state, &
    reactive_tree%levels(4)%patches(1)%temperature, &
    reactive_tree%topology%relations(3)%children(1)%geometry, &
    0.4_dp, node_dt, node_ok)
  call require(ok .and. reference_ok .and. node_ok .and. &
    tree_dt < initial_tree_dt .and. &
    abs(tree_dt - reference_tree_dt) <= &
      2.0e-12_dp * max(1.0_dp, abs(reference_tree_dt)) .and. &
    abs(tree_dt - real(ratio**3, dp) * node_dt) <= &
      2.0e-12_dp * max(1.0_dp, abs(tree_dt)), &
    "cumulative patch-tree subcycle scaling")

  reactive_tree_snapshot = reactive_tree
  call compute_reactive_amr_eb_patch_tree_cfl_timestep_2d( &
    species, reactive_tree, ieee_value(0.0_dp, ieee_quiet_nan), tree_dt, ok)
  call require(.not. ok .and. tree_dt == 0.0_dp .and. &
    reactive_tree_solutions_match(reactive_tree, reactive_tree_snapshot), &
    "invalid patch-tree CFL rollback")

  primitive(1:5) = [0.31_dp, 2.0_dp, -1.0_dp, 0.0_dp, 135000.0_dp]
  call reactive_primitive_to_conserved( &
    species, primitive, state_cell, temperature_cell, sound_speed, ok)
  call require(ok, "restore three-level reference reactive state")

  reactive_tree_snapshot = reactive_tree
  call rebuild_reactive_amr_eb_patch_tree_2d( &
    species, reactive_tree, shifted_tree_plans, ok, topology_changed)
  call require(ok .and. .not. topology_changed .and. &
    reactive_tree_solutions_match(reactive_tree, reactive_tree_snapshot), &
    "unchanged reactive patch-tree rebuild")
  call rebuild_reactive_amr_eb_patch_tree_2d( &
    species, reactive_tree, invalid_tree_plans, ok, topology_changed)
  call require(.not. ok .and. .not. topology_changed .and. &
    reactive_tree_solutions_match(reactive_tree, reactive_tree_snapshot) &
    .and. reactive_tree%topology%relations(3)%children(1)%parent_patch == 2, &
    "invalid reactive patch-tree rebuild rollback")

  do level = 1, reactive_tree%level_count()
    do patch = 1, reactive_tree%levels(level)%patch_count()
      reactive_tree%levels(level)%patches(patch)%state = spread(spread( &
        state_cell, 2, &
        size(reactive_tree%levels(level)%patches(patch)%state, 2)), 3, &
        size(reactive_tree%levels(level)%patches(patch)%state, 3))
      reactive_tree%levels(level)%patches(patch)%temperature = temperature_cell
    end do
  end do
  reactive_tree%levels(2)%patches(2)%state = &
    1.01_dp * reactive_tree%levels(2)%patches(2)%state
  reactive_tree%levels(3)%patches(2)%state = &
    0.99_dp * reactive_tree%levels(3)%patches(2)%state
  reactive_tree%levels(4)%patches(1)%state = &
    1.02_dp * reactive_tree%levels(4)%patches(1)%state
  reactive_tree_snapshot = reactive_tree
  call composite_integral_reactive_amr_eb_patch_tree_2d( &
    reactive_tree, tree_integral_before, ok)
  call require(ok, "patch-tree hydro initial composite integral")
  call compute_reactive_amr_eb_patch_tree_cfl_timestep_2d( &
    species, reactive_tree, 0.01_dp, tree_dt, ok)
  call require(ok, "patch-tree hydro stable timestep")
  allocate(tree_level_advances(reactive_tree%level_count()))
  call advance_reactive_amr_eb_patch_tree_hydro_2d( &
    species, reactive_tree, "hllc", "pcm", "mc", 2, tree_dt, ok, &
    0.5_dp, tree_failure_context, tree_level_advances)
  call require(ok .and. all(tree_level_advances == [1, 4, 8, 8]), &
    "arbitrary-depth branching patch-tree hydro schedule: " // &
      trim(tree_failure_context))
  call composite_integral_reactive_amr_eb_patch_tree_2d( &
    reactive_tree, tree_integral_after, ok)
  scale = max(1.0_dp, maxval(abs(tree_integral_before)))
  call require(ok .and. &
    abs(tree_integral_after(irho) - tree_integral_before(irho)) <= &
      2.0e-8_dp * scale .and. &
    abs(tree_integral_after(iet) - tree_integral_before(iet)) <= &
      2.0e-8_dp * scale, &
    "arbitrary-depth patch-tree hydro mass and energy conservation")
  do k = 1, size(species)
    call require(abs(tree_integral_after(reactive_species_component(k)) - &
      tree_integral_before(reactive_species_component(k))) <= &
      2.0e-8_dp * scale, &
      "arbitrary-depth patch-tree hydro species conservation")
  end do
  call require(.not. reactive_tree_solutions_match( &
      reactive_tree, reactive_tree_snapshot), &
    "arbitrary-depth patch-tree hydro state advance")
  call require(reactive_tree%is_valid(), &
    "arbitrary-depth patch-tree hydro thermodynamics")

  reactive_tree_snapshot = reactive_tree
  call advance_reactive_amr_eb_patch_tree_hydro_2d( &
    species, reactive_tree, "unknown", "pcm", "mc", 2, tree_dt, ok, &
    0.5_dp, tree_failure_context, tree_level_advances)
  call require(.not. ok .and. all(tree_level_advances == 0) .and. &
    reactive_tree_solutions_match(reactive_tree, reactive_tree_snapshot), &
    "arbitrary-depth patch-tree hydro rollback")

  deallocate(integral_before, integral_after)
  allocate(integral_before(nvar), integral_after(nvar))
  call composite_three_level_eb_integral_2d( &
    reactive_root, root_geometry, reactive_level_one, level_one_geometry, &
    root_patch, reactive_level_two, level_two_geometry, level_one_patch, &
    integral_before, ok)
  call require(ok, "three-level reactive composite integral")
  call average_down_three_level_reactive_eb_state_2d( &
    species, reactive_root, root_temperature, root_geometry, &
    reactive_level_one, level_one_temperature, level_one_geometry, &
    root_patch, reactive_level_two, level_two_temperature, &
    level_two_geometry, level_one_patch, reactive_root_sync, &
    root_temperature_sync, reactive_level_one_sync, &
    level_one_temperature_sync, ok)
  call require(ok, "three-level reactive synchronization")
  do i = 1, nvar
    integral_after(i) = sum(root_geometry%volume_fraction * &
      reactive_root_sync(i, :, :)) * root_geometry%dx * root_geometry%dy
  end do
  scale = max(1.0_dp, maxval(abs(integral_before)))
  call require(maxval(abs(integral_after - integral_before)) <= &
    8.0e-12_dp * scale, "three-level reactive conservation")
  call require(all(ieee_is_finite(root_temperature_sync)) .and. &
    all(ieee_is_finite(level_one_temperature_sync)) .and. &
    all(root_temperature_sync > 0.0_dp) .and. &
    all(level_one_temperature_sync > 0.0_dp), &
    "three-level reactive thermodynamics")
  call require(all(reactive_root_sync(:, 1, :) == reactive_root(:, 1, :)) &
    .and. all(reactive_level_one_sync(:, 1:2, :) == &
      reactive_level_one(:, 1:2, :)), &
    "three-level reactive unrefined ownership")

  level_two_temperature(1, 1) = -1.0_dp
  call average_down_three_level_reactive_eb_state_2d( &
    species, reactive_root, root_temperature, root_geometry, &
    reactive_level_one, level_one_temperature, level_one_geometry, &
    root_patch, reactive_level_two, level_two_temperature, &
    level_two_geometry, level_one_patch, reactive_root_sync, &
    root_temperature_sync, reactive_level_one_sync, &
    level_one_temperature_sync, ok)
  call require(.not. ok .and. all(reactive_root_sync == reactive_root) .and. &
    all(root_temperature_sync == root_temperature) .and. &
    all(reactive_level_one_sync == reactive_level_one) .and. &
    all(level_one_temperature_sync == level_one_temperature), &
    "three-level reactive rollback")

  level_two_temperature(1, 1) = temperature_cell
  primitive(2:4) = 0.0_dp
  call reactive_primitive_to_conserved( &
    species, primitive, state_cell, temperature_cell, sound_speed, ok)
  call require(ok, "stationary three-level hydro state")
  reactive_root = spread(spread(state_cell, 2, root_nx), 3, root_ny)
  reactive_level_one = 1.01_dp * &
    spread(spread(state_cell, 2, level_one_nx), 3, level_one_ny)
  reactive_level_two = 0.99_dp * &
    spread(spread(state_cell, 2, level_two_nx), 3, level_two_ny)
  root_temperature = temperature_cell
  level_one_temperature = temperature_cell
  level_two_temperature = temperature_cell
  call composite_three_level_eb_integral_2d( &
    reactive_root, root_geometry, reactive_level_one, level_one_geometry, &
    root_patch, reactive_level_two, level_two_geometry, level_one_patch, &
    integral_before, ok)
  call require(ok, "three-level hydro initial composite integral")
  dt = 0.015_dp * min(root_geometry%dx, root_geometry%dy) / sound_speed
  call advance_three_level_reactive_eb_hydro_2d( &
    species, reactive_root, root_temperature, root_geometry, &
    reactive_level_one, level_one_temperature, level_one_geometry, &
    root_patch, reactive_level_two, level_two_temperature, &
    level_two_geometry, level_one_patch, "hllc", "pcm", "mc", 2, dt, &
    reactive_root_sync, root_temperature_sync, reactive_level_one_sync, &
    level_one_temperature_sync, reactive_level_two_sync, &
    level_two_temperature_sync, ok)
  call require(ok, "recursive three-level reactive EB hydro")
  call composite_three_level_eb_integral_2d( &
    reactive_root_sync, root_geometry, reactive_level_one_sync, &
    level_one_geometry, root_patch, reactive_level_two_sync, &
    level_two_geometry, level_one_patch, integral_after, ok)
  scale = max(1.0_dp, maxval(abs(integral_before)))
  call require(ok .and. &
    abs(integral_after(irho) - integral_before(irho)) <= &
      8.0e-10_dp * scale .and. &
    abs(integral_after(iet) - integral_before(iet)) <= &
      8.0e-10_dp * scale, &
    "three-level hydro mass and energy conservation")
  do k = 1, size(species)
    call require(abs(integral_after(reactive_species_component(k)) - &
      integral_before(reactive_species_component(k))) <= &
      8.0e-10_dp * scale, "three-level hydro species conservation")
  end do
  call require(all(ieee_is_finite(root_temperature_sync)) .and. &
    all(ieee_is_finite(level_one_temperature_sync)) .and. &
    all(ieee_is_finite(level_two_temperature_sync)) .and. &
    all(root_temperature_sync > 0.0_dp) .and. &
    all(level_one_temperature_sync > 0.0_dp) .and. &
    all(level_two_temperature_sync > 0.0_dp), &
    "three-level hydro thermodynamics")

  call average_down_three_level_reactive_eb_state_2d( &
    species, reactive_root_sync, root_temperature_sync, root_geometry, &
    reactive_level_one_sync, level_one_temperature_sync, &
    level_one_geometry, root_patch, reactive_level_two_sync, &
    level_two_temperature_sync, level_two_geometry, level_one_patch, &
    reactive_root, root_temperature, reactive_level_one, &
    level_one_temperature, ok)
  call require(ok .and. maxval(abs(reactive_root - &
    reactive_root_sync)) <= 2.0e-12_dp * scale .and. &
    maxval(abs(reactive_level_one - reactive_level_one_sync)) <= &
      2.0e-12_dp * scale, "three-level hydro final synchronization")

  reactive_root = spread(spread(state_cell, 2, root_nx), 3, root_ny)
  reactive_level_one = 1.01_dp * &
    spread(spread(state_cell, 2, level_one_nx), 3, level_one_ny)
  reactive_level_two = 0.99_dp * &
    spread(spread(state_cell, 2, level_two_nx), 3, level_two_ny)
  root_temperature = temperature_cell
  level_one_temperature = temperature_cell
  level_two_temperature = temperature_cell
  call advance_three_level_reactive_eb_hydro_2d( &
    species, reactive_root, root_temperature, root_geometry, &
    reactive_level_one, level_one_temperature, level_one_geometry, &
    root_patch, reactive_level_two, level_two_temperature, &
    level_two_geometry, level_one_patch, "unknown", "pcm", "mc", 2, dt, &
    reactive_root_sync, root_temperature_sync, reactive_level_one_sync, &
    level_one_temperature_sync, reactive_level_two_sync, &
    level_two_temperature_sync, ok)
  call require(.not. ok .and. all(reactive_root_sync == reactive_root) .and. &
    all(root_temperature_sync == root_temperature) .and. &
    all(reactive_level_one_sync == reactive_level_one) .and. &
    all(level_one_temperature_sync == level_one_temperature) .and. &
    all(reactive_level_two_sync == reactive_level_two) .and. &
    all(level_two_temperature_sync == level_two_temperature), &
    "three-level hydro rollback")

  call composite_three_level_eb_integral_2d( &
    reactive_root, root_geometry, reactive_level_one, level_one_geometry, &
    root_patch, reactive_level_two, level_two_geometry, level_one_patch, &
    integral_before, ok)
  call require(ok, "three-level Strang initial composite integral")
  dt = 1.0e-8_dp
  call advance_three_level_reactive_eb_strang_2d( &
    species, reactions, reactive_root, root_temperature, root_geometry, &
    reactive_level_one, level_one_temperature, level_one_geometry, &
    root_patch, reactive_level_two, level_two_temperature, &
    level_two_geometry, level_one_patch, "hllc", "pcm", "mc", 2, dt, &
    .true., 1.0e-7_dp, 1.0e-13_dp, reactive_root_sync, &
    root_temperature_sync, reactive_level_one_sync, &
    level_one_temperature_sync, reactive_level_two_sync, &
    level_two_temperature_sync, ok, 0.5_dp)
  call require(ok, "three-level reactive EB Strang transaction")
  call composite_three_level_eb_integral_2d( &
    reactive_root_sync, root_geometry, reactive_level_one_sync, &
    level_one_geometry, root_patch, reactive_level_two_sync, &
    level_two_geometry, level_one_patch, integral_after, ok)
  scale = max(1.0_dp, maxval(abs(integral_before)))
  call require(ok .and. &
    abs(integral_after(irho) - integral_before(irho)) <= &
      2.0e-8_dp * scale .and. &
    abs(integral_after(iet) - integral_before(iet)) <= &
      2.0e-8_dp * scale, "three-level Strang mass and energy conservation")
  species_integral_sum = 0.0_dp
  species_change = 0.0_dp
  do k = 1, size(species)
    species_integral_sum = species_integral_sum + &
      integral_after(reactive_species_component(k))
    species_change = max(species_change, abs( &
      integral_after(reactive_species_component(k)) - &
      integral_before(reactive_species_component(k))))
  end do
  call require(abs(species_integral_sum - integral_after(irho)) <= &
    2.0e-10_dp * scale .and. species_change > 1.0e-15_dp * scale, &
    "three-level chemistry activity and species closure")
  call require(all(ieee_is_finite(root_temperature_sync)) .and. &
    all(ieee_is_finite(level_one_temperature_sync)) .and. &
    all(ieee_is_finite(level_two_temperature_sync)) .and. &
    all(root_temperature_sync > 0.0_dp) .and. &
    all(level_one_temperature_sync > 0.0_dp) .and. &
    all(level_two_temperature_sync > 0.0_dp), &
    "three-level Strang thermodynamics")

  call average_down_three_level_reactive_eb_state_2d( &
    species, reactive_root_sync, root_temperature_sync, root_geometry, &
    reactive_level_one_sync, level_one_temperature_sync, &
    level_one_geometry, root_patch, reactive_level_two_sync, &
    level_two_temperature_sync, level_two_geometry, level_one_patch, &
    reactive_root, root_temperature, reactive_level_one, &
    level_one_temperature, ok)
  call require(ok .and. maxval(abs(reactive_root - &
    reactive_root_sync)) <= 2.0e-12_dp * scale .and. &
    maxval(abs(reactive_level_one - reactive_level_one_sync)) <= &
      2.0e-12_dp * scale, "three-level Strang final synchronization")

  reactive_root = spread(spread(state_cell, 2, root_nx), 3, root_ny)
  reactive_level_one = 1.01_dp * &
    spread(spread(state_cell, 2, level_one_nx), 3, level_one_ny)
  reactive_level_two = 0.99_dp * &
    spread(spread(state_cell, 2, level_two_nx), 3, level_two_ny)
  root_temperature = temperature_cell
  level_one_temperature = temperature_cell
  level_two_temperature = temperature_cell
  call advance_three_level_reactive_eb_strang_2d( &
    species, reactions, reactive_root, root_temperature, root_geometry, &
    reactive_level_one, level_one_temperature, level_one_geometry, &
    root_patch, reactive_level_two, level_two_temperature, &
    level_two_geometry, level_one_patch, "unknown", "pcm", "mc", 2, dt, &
    .true., 1.0e-7_dp, 1.0e-13_dp, reactive_root_sync, &
    root_temperature_sync, reactive_level_one_sync, &
    level_one_temperature_sync, reactive_level_two_sync, &
    level_two_temperature_sync, ok, 0.5_dp)
  call require(.not. ok .and. all(reactive_root_sync == reactive_root) .and. &
    all(root_temperature_sync == root_temperature) .and. &
    all(reactive_level_one_sync == reactive_level_one) .and. &
    all(level_one_temperature_sync == level_one_temperature) .and. &
    all(reactive_level_two_sync == reactive_level_two) .and. &
    all(level_two_temperature_sync == level_two_temperature), &
    "three-level Strang rollback after chemistry")

  write(*, '(a)') "test_amr_eb_multilevel_2d: PASS"

contains

  subroutine reference_reactive_tree_cfl_timestep( &
      local_species, solution, cfl, reference_dt, valid)
    type(nasa7_species), intent(in) :: local_species(:)
    type(reactive_amr_eb_patch_tree_2d), intent(in) :: solution
    real(dp), intent(in) :: cfl
    real(dp), intent(out) :: reference_dt
    logical, intent(out) :: valid

    type(eb_geometry_2d) :: geometry
    real(dp) :: level_scale, local_dt
    logical :: local_ok
    integer :: level, patch, active_nodes

    reference_dt = huge(1.0_dp)
    valid = .false.
    level_scale = 1.0_dp
    active_nodes = 0
    do level = 1, solution%level_count()
      if (level > 1) level_scale = level_scale * real( &
        solution%topology%relations(level - 1)%refinement_ratio, dp)
      do patch = 1, solution%level_patch_count(level - 1)
        if (level == 1) then
          geometry = solution%topology%root_geometry
        else
          geometry = solution%topology%relations(level - 1)% &
            children(patch)%geometry
        end if
        if (count(geometry%cell_type /= eb_covered_cell) == 0) cycle
        active_nodes = active_nodes + 1
        call compute_reactive_eb_cfl_timestep_2d( &
          local_species, solution%levels(level)%patches(patch)%state, &
          solution%levels(level)%patches(patch)%temperature, geometry, &
          cfl, local_dt, local_ok)
        if (.not. local_ok) return
        reference_dt = min(reference_dt, level_scale * local_dt)
      end do
    end do
    valid = active_nodes > 0 .and. ieee_is_finite(reference_dt) .and. &
      reference_dt > 0.0_dp
  end subroutine reference_reactive_tree_cfl_timestep

  subroutine set_tree_child_plan( &
      plan, parent_patch, i_lower, i_upper, j_lower, j_upper, geometry)
    type(amr_eb_patch_tree_child_plan_2d), intent(out) :: plan
    integer, intent(in) :: parent_patch
    integer, intent(in) :: i_lower, i_upper, j_lower, j_upper
    type(eb_geometry_2d), intent(in) :: geometry

    plan%parent_patch = parent_patch
    plan%coarse_i_lower = i_lower
    plan%coarse_i_upper = i_upper
    plan%coarse_j_lower = j_lower
    plan%coarse_j_upper = j_upper
    plan%geometry = geometry
  end subroutine set_tree_child_plan

  subroutine build_patch_geometry( &
      parent_geometry, i_lower, i_upper, j_lower, j_upper, &
      refinement_ratio, child_geometry, patch, valid)
    type(eb_geometry_2d), intent(in) :: parent_geometry
    integer, intent(in) :: i_lower, i_upper, j_lower, j_upper
    integer, intent(in) :: refinement_ratio
    type(eb_geometry_2d), intent(out) :: child_geometry
    type(amr_eb_patch_2d), intent(out) :: patch
    logical, intent(out) :: valid

    real(dp), allocatable :: level_set(:, :)
    real(dp) :: x_lower, x_upper, y_lower, y_upper, local_x, local_y
    integer :: nx, ny, local_i, local_j

    nx = (i_upper - i_lower + 1) * refinement_ratio
    ny = (j_upper - j_lower + 1) * refinement_ratio
    x_lower = parent_geometry%x_lower + real(i_lower - 1, dp) * &
      parent_geometry%dx
    x_upper = parent_geometry%x_lower + real(i_upper, dp) * &
      parent_geometry%dx
    y_lower = parent_geometry%y_lower + real(j_lower - 1, dp) * &
      parent_geometry%dy
    y_upper = parent_geometry%y_lower + real(j_upper, dp) * &
      parent_geometry%dy
    allocate(level_set(0:nx, 0:ny))
    do local_j = 0, ny
      local_y = y_lower + real(local_j, dp) * &
        (y_upper - y_lower) / real(ny, dp)
      do local_i = 0, nx
        local_x = x_lower + real(local_i, dp) * &
          (x_upper - x_lower) / real(nx, dp)
        level_set(local_i, local_j) = local_x + local_y - 0.78_dp
      end do
    end do
    call build_eb_geometry_2d( &
      level_set, x_lower, x_upper, y_lower, y_upper, child_geometry, valid)
    if (.not. valid) return
    call build_amr_eb_patch_2d( &
      parent_geometry, child_geometry, i_lower, i_upper, j_lower, j_upper, &
      refinement_ratio, patch, valid)
  end subroutine build_patch_geometry

  logical function reactive_tree_solutions_match(first, second) &
      result(matches)
    type(reactive_amr_eb_patch_tree_2d), intent(in) :: first, second

    integer :: level, patch

    matches = first%is_valid() .and. second%is_valid() .and. &
      first%nvar == second%nvar .and. &
      first%level_count() == second%level_count()
    if (.not. matches) return
    do level = 1, first%level_count()
      matches = first%levels(level)%patch_count() == &
        second%levels(level)%patch_count()
      if (.not. matches) return
      do patch = 1, first%levels(level)%patch_count()
        matches = all(first%levels(level)%patches(patch)%state == &
            second%levels(level)%patches(patch)%state) .and. &
          all(first%levels(level)%patches(patch)%temperature == &
            second%levels(level)%patches(patch)%temperature)
        if (.not. matches) return
      end do
    end do
  end function reactive_tree_solutions_match

  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message

    call require(abs(actual - expected) <= &
      tolerance * max(1.0_dp, abs(expected)), message)
  end subroutine assert_close

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_amr_eb_multilevel_2d
