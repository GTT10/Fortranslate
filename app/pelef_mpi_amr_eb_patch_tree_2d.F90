program pelef_mpi_amr_eb_patch_tree_2d
  use, intrinsic :: iso_fortran_env, only: int64
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport
  use simulation_config_reactive_2d_mod, only: reactive_2d_config
  use reactive_1d_mod, only: reactive_nprim
  use reactive_2d_mod, only: initialize_reactive_2d
  use reactive_boundary_2d_mod, only: &
    reactive_boundary_set_2d, initialize_periodic_boundary_set_2d
  use eb_geometry_2d_mod, only: eb_geometry_2d, build_eb_geometry_2d
  use amr_eb_regrid_2d_mod, only: amr_eb_tagging_criteria_2d
  use amr_eb_patch_tree_2d_mod, only: &
    amr_eb_patch_tree_level_plan_2d, amr_eb_patch_tree_topology_2d, &
    initialize_amr_eb_patch_tree_topology_2d
  use amr_eb_patch_tree_reactive_2d_mod, only: &
    reactive_amr_eb_patch_tree_2d, &
    initialize_reactive_amr_eb_patch_tree_2d, &
    plan_tagged_reactive_amr_eb_patch_tree_2d, &
    regrid_tagged_reactive_amr_eb_patch_tree_2d, &
    compute_reactive_amr_eb_patch_tree_timestep_2d, &
    advance_reactive_amr_eb_patch_tree_chemistry_2d, &
    advance_reactive_amr_eb_patch_tree_hydro_2d, &
    advance_reactive_amr_eb_patch_tree_transport_2d, &
    advance_reactive_amr_eb_patch_tree_full_physics_2d, &
    advance_reactive_amr_eb_patch_tree_to_time_2d, &
    composite_integral_reactive_amr_eb_patch_tree_2d, &
    composite_reactive_amr_eb_patch_subtree_integral_2d
  use mpi_amr_eb_patch_tree_2d_mod, only: &
    mpi_amr_eb_patch_tree_distribution_2d, &
    mpi_sparse_reactive_amr_eb_patch_tree_2d, &
    initialize_mpi_amr_eb_patch_tree_distribution_2d, &
    mpi_amr_eb_patch_tree_distribution_matches_2d, &
    synchronize_owned_reactive_amr_eb_patch_tree_2d, &
    initialize_sparse_owned_reactive_amr_eb_patch_tree_2d, &
    materialize_sparse_owned_reactive_amr_eb_patch_tree_2d, &
    migrate_sparse_owned_reactive_amr_eb_patch_tree_2d, &
    plan_tagged_sparse_owned_reactive_amr_eb_patch_tree_2d, &
    regrid_sparse_owned_reactive_amr_eb_patch_tree_2d, &
    regrid_tagged_sparse_owned_reactive_amr_eb_patch_tree_2d, &
    compute_sparse_owned_reactive_amr_eb_patch_tree_timestep_2d, &
    advance_sparse_owned_reactive_amr_eb_patch_tree_chemistry_2d, &
    advance_sparse_owned_reactive_amr_eb_patch_tree_hydro_2d, &
    advance_sparse_owned_reactive_amr_eb_patch_tree_transport_2d, &
    advance_sparse_owned_reactive_amr_eb_patch_tree_full_physics_2d, &
    advance_sparse_owned_reactive_amr_eb_patch_tree_to_time_2d, &
    composite_sparse_amr_eb_patch_tree_integral_2d, &
    composite_sparse_amr_eb_patch_subtree_integral_2d
  implicit none

  type(MPI_Comm) :: comm
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(reactive_boundary_set_2d) :: boundaries
  type(reactive_2d_config) :: config
  type(eb_geometry_2d) :: root_geometry, level_one_geometry
  type(eb_geometry_2d) :: branch_a_geometry, branch_b_geometry
  type(eb_geometry_2d) :: deep_geometry
  type(amr_eb_patch_tree_level_plan_2d), allocatable :: plans(:)
  type(amr_eb_patch_tree_level_plan_2d), allocatable :: empty_plans(:)
  type(amr_eb_patch_tree_level_plan_2d), allocatable :: tagged_plans(:)
  type(amr_eb_patch_tree_level_plan_2d), allocatable :: serial_tagged_plans(:)
  type(amr_eb_patch_tree_topology_2d) :: topology
  type(amr_eb_patch_tree_topology_2d) :: root_only_topology
  type(reactive_amr_eb_patch_tree_2d) :: solution, accepted, failed
  type(reactive_amr_eb_patch_tree_2d) :: physical_solution
  type(reactive_amr_eb_patch_tree_2d) :: serial_chemistry
  type(reactive_amr_eb_patch_tree_2d) :: serial_hydro
  type(reactive_amr_eb_patch_tree_2d) :: serial_transport
  type(reactive_amr_eb_patch_tree_2d) :: serial_full_physics
  type(reactive_amr_eb_patch_tree_2d) :: serial_clock
  type(reactive_amr_eb_patch_tree_2d) :: tagged_serial
  type(reactive_amr_eb_patch_tree_2d) :: tagged_root
  type(reactive_amr_eb_patch_tree_2d) :: materialized
  type(mpi_amr_eb_patch_tree_distribution_2d) :: distribution
  type(mpi_amr_eb_patch_tree_distribution_2d) :: unweighted_distribution
  type(mpi_amr_eb_patch_tree_distribution_2d) :: rejected_distribution
  type(mpi_amr_eb_patch_tree_distribution_2d) :: migrated_distribution
  type(mpi_amr_eb_patch_tree_distribution_2d) :: invalid_distribution
  type(mpi_amr_eb_patch_tree_distribution_2d) :: tagged_distribution
  type(mpi_amr_eb_patch_tree_distribution_2d) :: tagged_new_distribution
  type(mpi_sparse_reactive_amr_eb_patch_tree_2d) :: sparse, sparse_snapshot
  type(mpi_sparse_reactive_amr_eb_patch_tree_2d) :: physical_sparse
  type(mpi_sparse_reactive_amr_eb_patch_tree_2d) :: tagged_sparse
  type(mpi_sparse_reactive_amr_eb_patch_tree_2d) :: tagged_sparse_snapshot
  type(amr_eb_tagging_criteria_2d) :: tagged_criteria
  type(amr_eb_tagging_criteria_2d) :: expanded_tagged_criteria
  type(amr_eb_tagging_criteria_2d) :: rejected_tagged_criteria
  real(dp), allocatable :: root_state(:, :, :), root_temperature(:, :)
  real(dp), allocatable :: serial_integral(:), sparse_integral(:)
  real(dp) :: base_density, expected_state, expected_temperature
  real(dp) :: hydro_cfl, mismatched_hydro_cfl, serial_dt, sparse_dt
  real(dp) :: transport_cfl
  real(dp) :: chemistry_interval, hydro_interval, mismatched_interval
  real(dp) :: mismatched_transport_interval, serial_transport_theta
  real(dp) :: sparse_transport_theta, transport_interval
  real(dp) :: full_physics_interval, serial_full_physics_theta
  real(dp) :: sparse_full_physics_theta
  real(dp) :: clock_final_time, serial_clock_minimum_dt
  real(dp) :: serial_clock_theta, serial_clock_time
  real(dp) :: sparse_clock_minimum_dt, sparse_clock_theta, sparse_clock_time
  real(dp) :: integral_scale
  real(dp) :: level_one_dx, level_one_dy
  integer :: ierr, rank, nranks, level, patch
  integer :: rejected_exponent
  integer :: local_publications, global_publications
  integer :: expected_transfers, global_transfers, local_allocated_cells
  integer :: local_nodes, local_transfers, new_owner, old_owner
  integer :: global_timestep_nodes, local_timestep_nodes
  integer :: child, expected_restriction_transfers, parent, relation
  integer :: global_restriction_transfers, local_restriction_transfers
  integer :: expected_integral_nodes, global_integral_nodes
  integer :: local_integral_nodes, selected_integral_patch
  integer :: expected_hydro_transfers, global_hydro_transfers
  integer :: local_hydro_transfers
  integer :: expected_transport_transfers, global_transport_transfers
  integer :: local_transport_transfers
  integer :: expected_regrid_edge_transfers
  integer :: global_candidate_transfers, local_candidate_transfers
  integer :: global_overlap_transfers, local_overlap_transfers
  integer :: global_prolongation_transfers, local_prolongation_transfers
  integer :: global_regrid_restriction_transfers
  integer :: local_regrid_restriction_transfers
  integer :: global_tagging_evaluations, local_tagging_evaluations
  integer :: serial_tagged_cells, tagged_cells, transferred_cells
  integer :: global_clock_timestep_evaluations
  integer :: local_clock_timestep_evaluations
  integer :: serial_clock_advanced_steps, serial_clock_steps
  integer :: sparse_clock_advanced_steps, sparse_clock_steps
  integer, allocatable :: global_chemistry_advances(:)
  integer, allocatable :: local_chemistry_advances(:)
  integer, allocatable :: global_hydro_advances(:)
  integer, allocatable :: local_hydro_advances(:)
  integer, allocatable :: global_transport_advances(:)
  integer, allocatable :: local_transport_advances(:)
  integer, allocatable :: serial_clock_chemistry_advances(:)
  integer, allocatable :: serial_clock_hydro_advances(:)
  integer, allocatable :: serial_clock_transport_advances(:)
  integer(int64) :: unweighted_work, weighted_work
  logical :: ok, local_ok, topology_changed

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI initialization failed"
  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank, ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI rank query failed"
  call MPI_Comm_size(comm, nranks, ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI size query failed"

  call load_h2o2_elementary_thermo(species, ok)
  call assert_all(ok, "MPI EB patch-tree thermodynamic database", comm)
  call load_h2o2_elementary_mechanism(reactions, ok)
  call assert_all(ok, "MPI EB patch-tree chemistry mechanism", comm)
  call load_h2o2_elementary_transport(transport, ok)
  call assert_all(ok, "MPI EB patch-tree transport database", comm)
  call initialize_periodic_boundary_set_2d( &
    reactive_nprim(size(species)), boundaries)
  call build_regular_geometry(8, 8, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    root_geometry, ok)
  call assert_all(ok, "MPI EB patch-tree root geometry", comm)
  call build_regular_geometry(12, 12, 0.125_dp, 0.875_dp, 0.125_dp, &
    0.875_dp, level_one_geometry, ok)
  call assert_all(ok, "MPI EB patch-tree level-one geometry", comm)
  level_one_dx = level_one_geometry%dx
  level_one_dy = level_one_geometry%dy
  call build_regular_geometry( &
    8, 8, level_one_geometry%x_lower + level_one_dx, &
    level_one_geometry%x_lower + 5.0_dp * level_one_dx, &
    level_one_geometry%y_lower + level_one_dy, &
    level_one_geometry%y_lower + 5.0_dp * level_one_dy, &
    branch_a_geometry, ok)
  call assert_all(ok, "MPI EB patch-tree branch-a geometry", comm)
  call build_regular_geometry( &
    4, 4, level_one_geometry%x_lower + 7.0_dp * level_one_dx, &
    level_one_geometry%x_lower + 9.0_dp * level_one_dx, &
    level_one_geometry%y_lower + 7.0_dp * level_one_dy, &
    level_one_geometry%y_lower + 9.0_dp * level_one_dy, &
    branch_b_geometry, ok)
  call assert_all(ok, "MPI EB patch-tree branch-b geometry", comm)
  call build_regular_geometry( &
    8, 8, branch_a_geometry%x_lower + branch_a_geometry%dx, &
    branch_a_geometry%x_lower + 5.0_dp * branch_a_geometry%dx, &
    branch_a_geometry%y_lower + branch_a_geometry%dy, &
    branch_a_geometry%y_lower + 5.0_dp * branch_a_geometry%dy, &
    deep_geometry, ok)
  call assert_all(ok, "MPI EB patch-tree deep geometry", comm)

  allocate(plans(3))
  plans%refinement_ratio = 2
  allocate(plans(1)%children(1))
  plans(1)%children(1)%parent_patch = 1
  plans(1)%children(1)%coarse_i_lower = 2
  plans(1)%children(1)%coarse_i_upper = 7
  plans(1)%children(1)%coarse_j_lower = 2
  plans(1)%children(1)%coarse_j_upper = 7
  plans(1)%children(1)%geometry = level_one_geometry
  allocate(plans(2)%children(2))
  plans(2)%children(1)%parent_patch = 1
  plans(2)%children(1)%coarse_i_lower = 2
  plans(2)%children(1)%coarse_i_upper = 5
  plans(2)%children(1)%coarse_j_lower = 2
  plans(2)%children(1)%coarse_j_upper = 5
  plans(2)%children(1)%geometry = branch_a_geometry
  plans(2)%children(2)%parent_patch = 1
  plans(2)%children(2)%coarse_i_lower = 8
  plans(2)%children(2)%coarse_i_upper = 9
  plans(2)%children(2)%coarse_j_lower = 8
  plans(2)%children(2)%coarse_j_upper = 9
  plans(2)%children(2)%geometry = branch_b_geometry
  allocate(plans(3)%children(1))
  plans(3)%children(1)%parent_patch = 1
  plans(3)%children(1)%coarse_i_lower = 2
  plans(3)%children(1)%coarse_i_upper = 5
  plans(3)%children(1)%coarse_j_lower = 2
  plans(3)%children(1)%coarse_j_upper = 5
  plans(3)%children(1)%geometry = deep_geometry
  call initialize_amr_eb_patch_tree_topology_2d( &
    root_geometry, plans, topology, ok)
  call assert_all(ok .and. topology%level_count() == 4 .and. &
    all([topology%level_patch_count(0), topology%level_patch_count(1), &
      topology%level_patch_count(2), topology%level_patch_count(3)] == &
      [1, 1, 2, 1]), "MPI EB patch-tree branching topology", comm)

  config%nx = root_geometry%nx
  config%ny = root_geometry%ny
  config%x_lower = root_geometry%x_lower
  config%x_upper = root_geometry%x_upper
  config%y_lower = root_geometry%y_lower
  config%y_upper = root_geometry%y_upper
  config%problem = "reactive_hotspot"
  config%initial_temperature = 900.0_dp
  config%initial_velocity_x = 0.0_dp
  config%initial_velocity_y = 0.0_dp
  config%hotspot_temperature_rise = 100.0_dp
  config%hotspot_center_x = 0.5_dp
  config%hotspot_center_y = 0.5_dp
  config%hotspot_width = 0.1_dp
  call initialize_reactive_2d( &
    species, config, root_state, root_temperature, root_geometry%dx, &
    root_geometry%dy, base_density, ok)
  call assert_all(ok, "MPI EB patch-tree root state", comm)
  call initialize_reactive_amr_eb_patch_tree_2d( &
    species, root_state, root_temperature, topology, solution, ok)
  call assert_all(ok, "MPI EB reactive patch-tree state", comm)
  physical_solution = solution

  call initialize_mpi_amr_eb_patch_tree_distribution_2d( &
    topology, comm, unweighted_distribution, ok, 0)
  call assert_all(ok .and. unweighted_distribution%is_valid() .and. &
    mpi_amr_eb_patch_tree_distribution_matches_2d( &
      unweighted_distribution, topology), &
    "MPI EB patch-tree unweighted ownership", comm)
  call initialize_mpi_amr_eb_patch_tree_distribution_2d( &
    topology, comm, distribution, ok, 2)
  call assert_all(ok .and. distribution%is_valid() .and. &
    mpi_amr_eb_patch_tree_distribution_matches_2d(distribution, topology), &
    "MPI EB patch-tree subcycle-weighted ownership", comm)
  unweighted_work = sum(unweighted_distribution%rank_work_counts)
  weighted_work = sum(distribution%rank_work_counts)
  call assert_all(unweighted_work == 352_int64 .and. &
    weighted_work == 6016_int64 .and. &
    sum(distribution%rank_patch_counts) == 5, &
    "MPI EB patch-tree ownership accounting", comm)

  do level = 1, solution%level_count()
    do patch = 1, solution%levels(level)%patch_count()
      expected_state = real(1000 * level + 10 * patch, dp)
      expected_temperature = 800.0_dp + real(10 * level + patch, dp)
      if (distribution%is_local(level - 1, patch)) then
        solution%levels(level)%patches(patch)%state = expected_state
        solution%levels(level)%patches(patch)%temperature = &
          expected_temperature
      else
        solution%levels(level)%patches(patch)%state = -real(rank + 1, dp)
        solution%levels(level)%patches(patch)%temperature = &
          700.0_dp + real(rank, dp)
      end if
    end do
  end do
  call synchronize_owned_reactive_amr_eb_patch_tree_2d( &
    distribution, solution, ok, local_publications)
  local_ok = ok
  do level = 1, solution%level_count()
    do patch = 1, solution%levels(level)%patch_count()
      expected_state = real(1000 * level + 10 * patch, dp)
      expected_temperature = 800.0_dp + real(10 * level + patch, dp)
      local_ok = local_ok .and. &
        all(solution%levels(level)%patches(patch)%state == expected_state) &
        .and. all(solution%levels(level)%patches(patch)%temperature == &
          expected_temperature)
    end do
  end do
  call assert_all(local_ok, "MPI EB patch-tree owner publication", comm)
  call MPI_Allreduce( &
    local_publications, global_publications, 1, MPI_INTEGER, MPI_SUM, &
    comm, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. global_publications == 5, &
    "MPI EB patch-tree publication accounting", comm)

  accepted = solution
  call initialize_sparse_owned_reactive_amr_eb_patch_tree_2d( &
    distribution, accepted, sparse, ok, local_allocated_cells)
  local_nodes = 0
  do level = 1, sparse%level_count()
    do patch = 1, sparse%levels(level)%patch_count()
      if (sparse%levels(level)%patches(patch)%has_data()) &
        local_nodes = local_nodes + 1
    end do
  end do
  call assert_all(ok .and. sparse%is_valid(distribution) .and. &
    local_allocated_cells == &
      distribution%rank_cell_counts(rank + 1) .and. &
    local_nodes == distribution%rank_patch_counts(rank + 1), &
    "MPI EB patch-tree sparse owner allocation", comm)
  call materialize_sparse_owned_reactive_amr_eb_patch_tree_2d( &
    distribution, sparse, materialized, ok, local_publications)
  call MPI_Allreduce( &
    local_publications, global_publications, 1, MPI_INTEGER, MPI_SUM, &
    comm, ierr)
  call assert_all(ok .and. ierr == MPI_SUCCESS .and. &
    global_publications == 5 .and. &
    tree_solutions_match(materialized, accepted), &
    "MPI EB patch-tree sparse materialization", comm)

  migrated_distribution = distribution
  migrated_distribution%rank_cell_counts = 0
  migrated_distribution%rank_patch_counts = 0
  migrated_distribution%rank_work_counts = 0_int64
  expected_transfers = 0
  do level = 1, migrated_distribution%level_count()
    do patch = 1, migrated_distribution%levels(level)%patch_count()
      old_owner = distribution%owner_of(level - 1, patch)
      new_owner = old_owner
      if (nranks > 1) new_owner = modulo(old_owner + 1, nranks)
      migrated_distribution%levels(level)%owners(patch) = new_owner
      migrated_distribution%rank_cell_counts(new_owner + 1) = &
        migrated_distribution%rank_cell_counts(new_owner + 1) + &
          migrated_distribution%levels(level)%cell_counts(patch)
      migrated_distribution%rank_patch_counts(new_owner + 1) = &
        migrated_distribution%rank_patch_counts(new_owner + 1) + 1
      migrated_distribution%rank_work_counts(new_owner + 1) = &
        migrated_distribution%rank_work_counts(new_owner + 1) + &
          migrated_distribution%levels(level)%work_counts(patch)
      if (new_owner /= old_owner) &
        expected_transfers = expected_transfers + 1
    end do
  end do
  call assert_all(migrated_distribution%is_valid() .and. &
    mpi_amr_eb_patch_tree_distribution_matches_2d( &
      migrated_distribution, topology), &
    "MPI EB patch-tree rotated ownership", comm)
  call migrate_sparse_owned_reactive_amr_eb_patch_tree_2d( &
    distribution, migrated_distribution, sparse, ok, local_transfers)
  call MPI_Allreduce( &
    local_transfers, global_transfers, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
  call assert_all(ok .and. ierr == MPI_SUCCESS .and. &
    global_transfers == expected_transfers .and. &
    sparse%is_valid(migrated_distribution), &
    "MPI EB patch-tree direct sparse migration", comm)
  call materialize_sparse_owned_reactive_amr_eb_patch_tree_2d( &
    migrated_distribution, sparse, materialized, ok, local_publications)
  call assert_all(ok .and. tree_solutions_match(materialized, accepted), &
    "MPI EB patch-tree migrated field parity", comm)

  sparse_snapshot = sparse
  invalid_distribution = migrated_distribution
  if (rank == 0) invalid_distribution%levels(1)%owners(1) = nranks
  call migrate_sparse_owned_reactive_amr_eb_patch_tree_2d( &
    migrated_distribution, invalid_distribution, sparse, ok, &
    local_transfers)
  call assert_all(.not. ok .and. local_transfers == 0 .and. &
    sparse_trees_match(sparse, sparse_snapshot), &
    "MPI EB patch-tree sparse migration rollback", comm)

  hydro_cfl = 0.4_dp
  transport_cfl = 0.2_dp
  call initialize_sparse_owned_reactive_amr_eb_patch_tree_2d( &
    migrated_distribution, physical_solution, physical_sparse, ok)
  call assert_all(ok, "MPI EB patch-tree physical sparse state", comm)
  call compute_reactive_amr_eb_patch_tree_timestep_2d( &
    species, transport, physical_solution, hydro_cfl, transport_cfl, &
    .true., .true., .true., serial_dt, ok)
  call assert_all(ok, "MPI EB patch-tree serial timestep reference", comm)
  call compute_sparse_owned_reactive_amr_eb_patch_tree_timestep_2d( &
    species, transport, migrated_distribution, physical_sparse, hydro_cfl, &
    transport_cfl, .true., .true., .true., sparse_dt, ok, &
    local_timestep_nodes)
  call MPI_Allreduce( &
    local_timestep_nodes, global_timestep_nodes, 1, MPI_INTEGER, MPI_SUM, &
    comm, ierr)
  call assert_all(ok .and. ierr == MPI_SUCCESS .and. &
    sparse_dt == serial_dt .and. global_timestep_nodes == 5, &
    "MPI EB patch-tree sparse timestep parity", comm)

  mismatched_hydro_cfl = hydro_cfl
  if (nranks > 1) then
    if (rank == 0) mismatched_hydro_cfl = 0.5_dp * hydro_cfl
  else
    mismatched_hydro_cfl = -hydro_cfl
  end if
  call compute_sparse_owned_reactive_amr_eb_patch_tree_timestep_2d( &
    species, transport, migrated_distribution, physical_sparse, &
    mismatched_hydro_cfl, transport_cfl, .true., .true., .true., &
    sparse_dt, ok, local_timestep_nodes)
  call assert_all(.not. ok .and. sparse_dt == 0.0_dp .and. &
    local_timestep_nodes == 0, &
    "MPI EB patch-tree timestep control consensus", comm)

  chemistry_interval = 5.0e-9_dp
  serial_chemistry = physical_solution
  call advance_reactive_amr_eb_patch_tree_chemistry_2d( &
    species, reactions, serial_chemistry, chemistry_interval, 1.0e-7_dp, &
    1.0e-13_dp, ok)
  call assert_all(ok, "MPI EB patch-tree serial chemistry reference", comm)
  expected_restriction_transfers = 0
  do relation = 1, size(topology%relations)
    do child = 1, topology%relations(relation)%child_patch_count()
      parent = topology%relations(relation)%children(child)%parent_patch
      if (migrated_distribution%owner_of(relation - 1, parent) /= &
          migrated_distribution%owner_of(relation, child)) &
        expected_restriction_transfers = expected_restriction_transfers + 1
    end do
  end do
  allocate(local_chemistry_advances(topology%level_count()))
  allocate(global_chemistry_advances(topology%level_count()))
  call advance_sparse_owned_reactive_amr_eb_patch_tree_chemistry_2d( &
    species, reactions, migrated_distribution, physical_sparse, &
    chemistry_interval, 1.0e-7_dp, 1.0e-13_dp, ok, &
    local_chemistry_advances, local_restriction_transfers)
  call MPI_Allreduce( &
    local_chemistry_advances, global_chemistry_advances, &
    size(local_chemistry_advances), MPI_INTEGER, MPI_SUM, comm, ierr)
  call MPI_Allreduce( &
    local_restriction_transfers, global_restriction_transfers, 1, &
    MPI_INTEGER, MPI_SUM, comm, ierr)
  call materialize_sparse_owned_reactive_amr_eb_patch_tree_2d( &
    migrated_distribution, physical_sparse, materialized, local_ok)
  call assert_all(ok .and. local_ok .and. ierr == MPI_SUCCESS .and. &
    all(global_chemistry_advances == [1, 1, 2, 1]) .and. &
    global_restriction_transfers == expected_restriction_transfers .and. &
    tree_solutions_match(materialized, serial_chemistry), &
    "MPI EB patch-tree owner-local chemistry parity", comm)

  allocate(serial_integral(serial_chemistry%nvar))
  allocate(sparse_integral(serial_chemistry%nvar))
  call composite_integral_reactive_amr_eb_patch_tree_2d( &
    serial_chemistry, serial_integral, ok)
  call composite_sparse_amr_eb_patch_tree_integral_2d( &
    migrated_distribution, physical_sparse, sparse_integral, local_ok, &
    local_integral_nodes)
  call MPI_Allreduce( &
    local_integral_nodes, global_integral_nodes, 1, MPI_INTEGER, MPI_SUM, &
    comm, ierr)
  integral_scale = max(1.0_dp, maxval(abs(serial_integral)))
  call assert_all(ok .and. local_ok .and. ierr == MPI_SUCCESS .and. &
    global_integral_nodes == 5 .and. &
    maxval(abs(sparse_integral - serial_integral)) <= &
      256.0_dp * epsilon(1.0_dp) * integral_scale, &
    "MPI EB patch-tree sparse composite integral", comm)
  do level = 1, serial_chemistry%level_count()
    do patch = 1, serial_chemistry%levels(level)%patch_count()
      call composite_reactive_amr_eb_patch_subtree_integral_2d( &
        serial_chemistry, level, patch, serial_integral, ok)
      call &
          composite_sparse_amr_eb_patch_subtree_integral_2d( &
        migrated_distribution, physical_sparse, level, patch, &
        sparse_integral, local_ok, local_integral_nodes)
      call MPI_Allreduce( &
        local_integral_nodes, global_integral_nodes, 1, MPI_INTEGER, MPI_SUM, &
        comm, ierr)
      expected_integral_nodes = subtree_node_count(topology, level, patch)
      integral_scale = max(1.0_dp, maxval(abs(serial_integral)))
      call assert_all(ok .and. local_ok .and. ierr == MPI_SUCCESS .and. &
        global_integral_nodes == expected_integral_nodes .and. &
        maxval(abs(sparse_integral - serial_integral)) <= &
          256.0_dp * epsilon(1.0_dp) * integral_scale, &
        "MPI EB patch-tree sparse subtree integral", comm)
    end do
  end do

  selected_integral_patch = 1
  if (nranks > 1) then
    if (rank == 0) selected_integral_patch = 2
  else
    selected_integral_patch = 0
  end if
  call composite_sparse_amr_eb_patch_subtree_integral_2d( &
    migrated_distribution, physical_sparse, 3, selected_integral_patch, &
    sparse_integral, ok, local_integral_nodes)
  call assert_all(.not. ok .and. all(sparse_integral == 0.0_dp) .and. &
    local_integral_nodes == 0, &
    "MPI EB patch-tree subtree selector consensus", comm)

  hydro_interval = 0.05_dp * serial_dt
  serial_hydro = serial_chemistry
  call advance_reactive_amr_eb_patch_tree_hydro_2d( &
    species, serial_hydro, "hllc", "pcm", "mc", 2, hydro_interval, ok)
  call assert_all(ok, "MPI EB patch-tree serial hydro reference", comm)
  allocate(local_hydro_advances(topology%level_count()))
  allocate(global_hydro_advances(topology%level_count()))
  expected_hydro_transfers = expected_sparse_hydro_transfers( &
    topology, migrated_distribution)
  call advance_sparse_owned_reactive_amr_eb_patch_tree_hydro_2d( &
    species, migrated_distribution, physical_sparse, "hllc", "pcm", &
    "mc", 2, hydro_interval, ok, 0.5_dp, &
    local_level_advances=local_hydro_advances, &
    local_entity_transfers=local_hydro_transfers)
  call MPI_Allreduce( &
    local_hydro_advances, global_hydro_advances, &
    size(local_hydro_advances), MPI_INTEGER, MPI_SUM, comm, ierr)
  call MPI_Allreduce( &
    local_hydro_transfers, global_hydro_transfers, 1, MPI_INTEGER, MPI_SUM, &
    comm, ierr)
  call materialize_sparse_owned_reactive_amr_eb_patch_tree_2d( &
    migrated_distribution, physical_sparse, materialized, local_ok)
  call assert_all(ok .and. local_ok .and. ierr == MPI_SUCCESS .and. &
    all(global_hydro_advances == [1, 2, 8, 8]) .and. &
    global_hydro_transfers == expected_hydro_transfers .and. &
    tree_solutions_close(materialized, serial_hydro, 5.0e-10_dp), &
    "MPI EB patch-tree owner-local hydro parity", comm)
  call composite_integral_reactive_amr_eb_patch_tree_2d( &
    serial_hydro, serial_integral, ok)
  call composite_sparse_amr_eb_patch_tree_integral_2d( &
    migrated_distribution, physical_sparse, sparse_integral, local_ok)
  integral_scale = max(1.0_dp, maxval(abs(serial_integral)))
  call assert_all(ok .and. local_ok .and. &
    maxval(abs(sparse_integral - serial_integral)) <= &
      5.0e-10_dp * integral_scale, &
    "MPI EB patch-tree owner-local hydro conservation", comm)

  sparse_snapshot = physical_sparse
  mismatched_interval = hydro_interval
  if (nranks > 1) then
    if (rank == 0) mismatched_interval = 0.5_dp * hydro_interval
  else
    mismatched_interval = -hydro_interval
  end if
  call advance_sparse_owned_reactive_amr_eb_patch_tree_hydro_2d( &
    species, migrated_distribution, physical_sparse, "hllc", "pcm", &
    "mc", 2, mismatched_interval, ok, 0.5_dp, &
    local_level_advances=local_hydro_advances, &
    local_entity_transfers=local_hydro_transfers)
  call assert_all(.not. ok .and. all(local_hydro_advances == 0) .and. &
    local_hydro_transfers == 0 .and. &
    sparse_trees_match(physical_sparse, sparse_snapshot), &
    "MPI EB patch-tree hydro control rollback", comm)

  transport_interval = 0.02_dp * serial_dt
  serial_transport = serial_hydro
  call advance_reactive_amr_eb_patch_tree_transport_2d( &
    species, transport, serial_transport, transport_interval, .true., &
    .true., .true., .true., boundaries, 0.5_dp, 2, &
    serial_transport_theta, ok)
  call assert_all(ok, "MPI EB patch-tree serial transport reference", comm)
  allocate(local_transport_advances(topology%level_count()))
  allocate(global_transport_advances(topology%level_count()))
  expected_transport_transfers = &
    2 * expected_hydro_transfers + expected_restriction_transfers
  call advance_sparse_owned_reactive_amr_eb_patch_tree_transport_2d( &
    species, transport, migrated_distribution, physical_sparse, &
    transport_interval, .true., .true., .true., .true., boundaries, &
    0.5_dp, 2, sparse_transport_theta, ok, &
    local_level_advances=local_transport_advances, &
    local_entity_transfers=local_transport_transfers)
  call MPI_Allreduce( &
    local_transport_advances, global_transport_advances, &
    size(local_transport_advances), MPI_INTEGER, MPI_SUM, comm, ierr)
  call MPI_Allreduce( &
    local_transport_transfers, global_transport_transfers, 1, MPI_INTEGER, &
    MPI_SUM, comm, ierr)
  call materialize_sparse_owned_reactive_amr_eb_patch_tree_2d( &
    migrated_distribution, physical_sparse, materialized, local_ok)
  call assert_all(ok .and. local_ok .and. ierr == MPI_SUCCESS .and. &
    all(global_transport_advances == [2, 4, 16, 16]) .and. &
    global_transport_transfers == expected_transport_transfers .and. &
    abs(sparse_transport_theta - serial_transport_theta) <= &
      256.0_dp * epsilon(1.0_dp) * &
        max(1.0_dp, abs(serial_transport_theta)) .and. &
    tree_solutions_close(materialized, serial_transport, 1.0e-9_dp), &
    "MPI EB patch-tree owner-local transport parity", comm)
  call composite_integral_reactive_amr_eb_patch_tree_2d( &
    serial_transport, serial_integral, ok)
  call composite_sparse_amr_eb_patch_tree_integral_2d( &
    migrated_distribution, physical_sparse, sparse_integral, local_ok)
  integral_scale = max(1.0_dp, maxval(abs(serial_integral)))
  call assert_all(ok .and. local_ok .and. &
    maxval(abs(sparse_integral - serial_integral)) <= &
      1.0e-9_dp * integral_scale, &
    "MPI EB patch-tree owner-local transport conservation", comm)

  sparse_snapshot = physical_sparse
  mismatched_transport_interval = transport_interval
  if (nranks > 1) then
    if (rank == 0) &
      mismatched_transport_interval = 0.5_dp * transport_interval
  else
    mismatched_transport_interval = -transport_interval
  end if
  call advance_sparse_owned_reactive_amr_eb_patch_tree_transport_2d( &
    species, transport, migrated_distribution, physical_sparse, &
    mismatched_transport_interval, .true., .true., .true., .true., &
    boundaries, 0.5_dp, 2, sparse_transport_theta, ok, &
    local_level_advances=local_transport_advances, &
    local_entity_transfers=local_transport_transfers)
  call assert_all(.not. ok .and. sparse_transport_theta == 1.0_dp .and. &
    all(local_transport_advances == 0) .and. &
    local_transport_transfers == 0 .and. &
    sparse_trees_match(physical_sparse, sparse_snapshot), &
    "MPI EB patch-tree transport control rollback", comm)

  full_physics_interval = chemistry_interval
  serial_full_physics = serial_transport
  call advance_reactive_amr_eb_patch_tree_full_physics_2d( &
    species, reactions, transport, serial_full_physics, "hllc", "pcm", &
    "mc", 2, full_physics_interval, .true., 1.0e-7_dp, 1.0e-13_dp, &
    .true., .true., .true., .true., boundaries, 0.5_dp, &
    serial_full_physics_theta, ok)
  call assert_all(ok, "MPI EB patch-tree serial full-physics reference", comm)
  call advance_sparse_owned_reactive_amr_eb_patch_tree_full_physics_2d( &
    species, reactions, transport, migrated_distribution, physical_sparse, &
    "hllc", "pcm", "mc", 2, full_physics_interval, .true., 1.0e-7_dp, &
    1.0e-13_dp, .true., .true., .true., .true., boundaries, 0.5_dp, &
    sparse_full_physics_theta, ok, &
    local_chemistry_level_advances=local_chemistry_advances, &
    local_transport_level_advances=local_transport_advances, &
    local_hydro_level_advances=local_hydro_advances, &
    local_chemistry_transfers=local_restriction_transfers, &
    local_transport_transfers=local_transport_transfers, &
    local_hydro_transfers=local_hydro_transfers)
  call MPI_Allreduce( &
    local_chemistry_advances, global_chemistry_advances, &
    size(local_chemistry_advances), MPI_INTEGER, MPI_SUM, comm, ierr)
  call MPI_Allreduce( &
    local_transport_advances, global_transport_advances, &
    size(local_transport_advances), MPI_INTEGER, MPI_SUM, comm, ierr)
  call MPI_Allreduce( &
    local_hydro_advances, global_hydro_advances, &
    size(local_hydro_advances), MPI_INTEGER, MPI_SUM, comm, ierr)
  call MPI_Allreduce( &
    local_restriction_transfers, global_restriction_transfers, 1, &
    MPI_INTEGER, MPI_SUM, comm, ierr)
  call MPI_Allreduce( &
    local_transport_transfers, global_transport_transfers, 1, MPI_INTEGER, &
    MPI_SUM, comm, ierr)
  call MPI_Allreduce( &
    local_hydro_transfers, global_hydro_transfers, 1, MPI_INTEGER, MPI_SUM, &
    comm, ierr)
  call materialize_sparse_owned_reactive_amr_eb_patch_tree_2d( &
    migrated_distribution, physical_sparse, materialized, local_ok)
  call assert_all(ok .and. local_ok .and. ierr == MPI_SUCCESS .and. &
    all(global_chemistry_advances == [2, 2, 4, 2]) .and. &
    all(global_transport_advances == [4, 8, 32, 32]) .and. &
    all(global_hydro_advances == [1, 2, 8, 8]) .and. &
    global_restriction_transfers == 2 * expected_restriction_transfers .and. &
    global_transport_transfers == 2 * expected_transport_transfers .and. &
    global_hydro_transfers == expected_hydro_transfers .and. &
    abs(sparse_full_physics_theta - serial_full_physics_theta) <= &
      512.0_dp * epsilon(1.0_dp) * &
        max(1.0_dp, abs(serial_full_physics_theta)) .and. &
    tree_solutions_close(materialized, serial_full_physics, 3.0e-9_dp), &
    "MPI EB patch-tree owner-local full-physics parity", comm)
  call composite_integral_reactive_amr_eb_patch_tree_2d( &
    serial_full_physics, serial_integral, ok)
  call composite_sparse_amr_eb_patch_tree_integral_2d( &
    migrated_distribution, physical_sparse, sparse_integral, local_ok)
  integral_scale = max(1.0_dp, maxval(abs(serial_integral)))
  call assert_all(ok .and. local_ok .and. &
    maxval(abs(sparse_integral - serial_integral)) <= &
      3.0e-9_dp * integral_scale, &
    "MPI EB patch-tree owner-local full-physics conservation", comm)

  sparse_snapshot = physical_sparse
  mismatched_interval = full_physics_interval
  if (nranks > 1) then
    if (rank == 0) mismatched_interval = 0.5_dp * full_physics_interval
  else
    mismatched_interval = -full_physics_interval
  end if
  call advance_sparse_owned_reactive_amr_eb_patch_tree_full_physics_2d( &
    species, reactions, transport, migrated_distribution, physical_sparse, &
    "hllc", "pcm", "mc", 2, mismatched_interval, .true., 1.0e-7_dp, &
    1.0e-13_dp, .true., .true., .true., .true., boundaries, 0.5_dp, &
    sparse_full_physics_theta, ok, &
    local_chemistry_level_advances=local_chemistry_advances, &
    local_transport_level_advances=local_transport_advances, &
    local_hydro_level_advances=local_hydro_advances, &
    local_chemistry_transfers=local_restriction_transfers, &
    local_transport_transfers=local_transport_transfers, &
    local_hydro_transfers=local_hydro_transfers)
  call assert_all(.not. ok .and. sparse_full_physics_theta == 1.0_dp .and. &
    all(local_chemistry_advances == 0) .and. &
    all(local_transport_advances == 0) .and. &
    all(local_hydro_advances == 0) .and. &
    local_restriction_transfers == 0 .and. &
    local_transport_transfers == 0 .and. local_hydro_transfers == 0 .and. &
    sparse_trees_match(physical_sparse, sparse_snapshot), &
    "MPI EB patch-tree full-physics control rollback", comm)

  clock_final_time = 0.5_dp * full_physics_interval
  serial_clock = serial_full_physics
  serial_clock_time = 0.0_dp
  serial_clock_steps = 0
  allocate(serial_clock_chemistry_advances(topology%level_count()))
  allocate(serial_clock_transport_advances(topology%level_count()))
  allocate(serial_clock_hydro_advances(topology%level_count()))
  call advance_reactive_amr_eb_patch_tree_to_time_2d( &
    species, reactions, transport, serial_clock, "hllc", "pcm", "mc", 2, &
    serial_clock_time, clock_final_time, serial_clock_steps, 2, hydro_cfl, &
    transport_cfl, .true., 1.0e-7_dp, 1.0e-13_dp, .true., .true., .true., &
    .true., boundaries, 0.5_dp, serial_clock_minimum_dt, &
    serial_clock_theta, ok, advanced_steps=serial_clock_advanced_steps, &
    chemistry_level_advances=serial_clock_chemistry_advances, &
    transport_level_advances=serial_clock_transport_advances, &
    hydro_level_advances=serial_clock_hydro_advances)
  call assert_all(ok .and. serial_clock_time == clock_final_time .and. &
    serial_clock_steps == 1 .and. serial_clock_advanced_steps == 1, &
    "MPI EB patch-tree serial public-clock reference", comm)

  sparse_clock_time = 0.0_dp
  sparse_clock_steps = 0
  call advance_sparse_owned_reactive_amr_eb_patch_tree_to_time_2d( &
    species, reactions, transport, migrated_distribution, physical_sparse, &
    "hllc", "pcm", "mc", 2, sparse_clock_time, clock_final_time, &
    sparse_clock_steps, 2, hydro_cfl, transport_cfl, .true., 1.0e-7_dp, &
    1.0e-13_dp, .true., .true., .true., .true., boundaries, 0.5_dp, &
    sparse_clock_minimum_dt, sparse_clock_theta, ok, &
    advanced_steps=sparse_clock_advanced_steps, &
    local_timestep_evaluations=local_clock_timestep_evaluations, &
    local_chemistry_level_advances=local_chemistry_advances, &
    local_transport_level_advances=local_transport_advances, &
    local_hydro_level_advances=local_hydro_advances, &
    local_chemistry_transfers=local_restriction_transfers, &
    local_transport_transfers=local_transport_transfers, &
    local_hydro_transfers=local_hydro_transfers)
  call MPI_Allreduce( &
    local_clock_timestep_evaluations, global_clock_timestep_evaluations, 1, &
    MPI_INTEGER, MPI_SUM, comm, ierr)
  call MPI_Allreduce( &
    local_chemistry_advances, global_chemistry_advances, &
    size(local_chemistry_advances), MPI_INTEGER, MPI_SUM, comm, ierr)
  call MPI_Allreduce( &
    local_transport_advances, global_transport_advances, &
    size(local_transport_advances), MPI_INTEGER, MPI_SUM, comm, ierr)
  call MPI_Allreduce( &
    local_hydro_advances, global_hydro_advances, &
    size(local_hydro_advances), MPI_INTEGER, MPI_SUM, comm, ierr)
  call MPI_Allreduce( &
    local_restriction_transfers, global_restriction_transfers, 1, &
    MPI_INTEGER, MPI_SUM, comm, ierr)
  call MPI_Allreduce( &
    local_transport_transfers, global_transport_transfers, 1, MPI_INTEGER, &
    MPI_SUM, comm, ierr)
  call MPI_Allreduce( &
    local_hydro_transfers, global_hydro_transfers, 1, MPI_INTEGER, MPI_SUM, &
    comm, ierr)
  call materialize_sparse_owned_reactive_amr_eb_patch_tree_2d( &
    migrated_distribution, physical_sparse, materialized, local_ok)
  call assert_all(ok .and. local_ok .and. ierr == MPI_SUCCESS .and. &
    sparse_clock_time == clock_final_time .and. sparse_clock_steps == 1 .and. &
    sparse_clock_advanced_steps == 1 .and. &
    sparse_clock_minimum_dt == serial_clock_minimum_dt .and. &
    global_clock_timestep_evaluations == 5 .and. &
    all(global_chemistry_advances == serial_clock_chemistry_advances) .and. &
    all(global_transport_advances == serial_clock_transport_advances) .and. &
    all(global_hydro_advances == serial_clock_hydro_advances) .and. &
    global_restriction_transfers == 2 * expected_restriction_transfers .and. &
    global_transport_transfers == 2 * expected_transport_transfers .and. &
    global_hydro_transfers == expected_hydro_transfers .and. &
    abs(sparse_clock_theta - serial_clock_theta) <= &
      1024.0_dp * epsilon(1.0_dp) * &
        max(1.0_dp, abs(serial_clock_theta)) .and. &
    tree_solutions_close(materialized, serial_clock, 6.0e-9_dp), &
    "MPI EB patch-tree owner-local public-clock parity", comm)
  call composite_integral_reactive_amr_eb_patch_tree_2d( &
    serial_clock, serial_integral, ok)
  call composite_sparse_amr_eb_patch_tree_integral_2d( &
    migrated_distribution, physical_sparse, sparse_integral, local_ok)
  integral_scale = max(1.0_dp, maxval(abs(serial_integral)))
  call assert_all(ok .and. local_ok .and. &
    maxval(abs(sparse_integral - serial_integral)) <= &
      6.0e-9_dp * integral_scale, &
    "MPI EB patch-tree owner-local public-clock conservation", comm)

  sparse_snapshot = physical_sparse
  sparse_clock_time = 0.0_dp
  sparse_clock_steps = 0
  mismatched_interval = clock_final_time
  if (nranks > 1) then
    if (rank == 0) mismatched_interval = 0.5_dp * clock_final_time
  else
    mismatched_interval = -clock_final_time
  end if
  call advance_sparse_owned_reactive_amr_eb_patch_tree_to_time_2d( &
    species, reactions, transport, migrated_distribution, physical_sparse, &
    "hllc", "pcm", "mc", 2, sparse_clock_time, mismatched_interval, &
    sparse_clock_steps, 2, hydro_cfl, transport_cfl, .true., 1.0e-7_dp, &
    1.0e-13_dp, .true., .true., .true., .true., boundaries, 0.5_dp, &
    sparse_clock_minimum_dt, sparse_clock_theta, ok, &
    advanced_steps=sparse_clock_advanced_steps, &
    local_timestep_evaluations=local_clock_timestep_evaluations, &
    local_chemistry_level_advances=local_chemistry_advances, &
    local_transport_level_advances=local_transport_advances, &
    local_hydro_level_advances=local_hydro_advances, &
    local_chemistry_transfers=local_restriction_transfers, &
    local_transport_transfers=local_transport_transfers, &
    local_hydro_transfers=local_hydro_transfers)
  call assert_all(.not. ok .and. sparse_clock_time == 0.0_dp .and. &
    sparse_clock_steps == 0 .and. sparse_clock_advanced_steps == 0 .and. &
    sparse_clock_minimum_dt == 0.0_dp .and. &
    sparse_clock_theta == 1.0_dp .and. &
    local_clock_timestep_evaluations == 0 .and. &
    all(local_chemistry_advances == 0) .and. &
    all(local_transport_advances == 0) .and. &
    all(local_hydro_advances == 0) .and. &
    local_restriction_transfers == 0 .and. &
    local_transport_transfers == 0 .and. local_hydro_transfers == 0 .and. &
    sparse_trees_match(physical_sparse, sparse_snapshot), &
    "MPI EB patch-tree public-clock control rollback", comm)

  sparse_clock_time = 0.0_dp
  sparse_clock_steps = 0
  call advance_sparse_owned_reactive_amr_eb_patch_tree_to_time_2d( &
    species, reactions, transport, migrated_distribution, physical_sparse, &
    "hllc", "pcm", "mc", 2, sparse_clock_time, clock_final_time, &
    sparse_clock_steps, 0, hydro_cfl, transport_cfl, .true., 1.0e-7_dp, &
    1.0e-13_dp, .true., .true., .true., .true., boundaries, 0.5_dp, &
    sparse_clock_minimum_dt, sparse_clock_theta, ok, &
    advanced_steps=sparse_clock_advanced_steps, &
    local_timestep_evaluations=local_clock_timestep_evaluations, &
    local_chemistry_level_advances=local_chemistry_advances, &
    local_transport_level_advances=local_transport_advances, &
    local_hydro_level_advances=local_hydro_advances, &
    local_chemistry_transfers=local_restriction_transfers, &
    local_transport_transfers=local_transport_transfers, &
    local_hydro_transfers=local_hydro_transfers)
  call assert_all(.not. ok .and. sparse_clock_time == 0.0_dp .and. &
    sparse_clock_steps == 0 .and. sparse_clock_advanced_steps == 0 .and. &
    sparse_clock_minimum_dt == 0.0_dp .and. &
    sparse_clock_theta == 1.0_dp .and. &
    local_clock_timestep_evaluations == 0 .and. &
    all(local_chemistry_advances == 0) .and. &
    all(local_transport_advances == 0) .and. &
    all(local_hydro_advances == 0) .and. &
    local_restriction_transfers == 0 .and. &
    local_transport_transfers == 0 .and. local_hydro_transfers == 0 .and. &
    sparse_trees_match(physical_sparse, sparse_snapshot), &
    "MPI EB patch-tree public-clock step-limit rollback", comm)

  sparse_snapshot = physical_sparse
  mismatched_interval = chemistry_interval
  if (nranks > 1) then
    if (rank == 0) mismatched_interval = 0.5_dp * chemistry_interval
  else
    mismatched_interval = -chemistry_interval
  end if
  call advance_sparse_owned_reactive_amr_eb_patch_tree_chemistry_2d( &
    species, reactions, migrated_distribution, physical_sparse, &
    mismatched_interval, 1.0e-7_dp, 1.0e-13_dp, ok, &
    local_chemistry_advances, local_restriction_transfers)
  call assert_all(.not. ok .and. all(local_chemistry_advances == 0) .and. &
    local_restriction_transfers == 0 .and. &
    sparse_trees_match(physical_sparse, sparse_snapshot), &
    "MPI EB patch-tree chemistry control rollback", comm)

  allocate(empty_plans(0))
  call initialize_amr_eb_patch_tree_topology_2d( &
    root_geometry, empty_plans, root_only_topology, ok)
  call assert_all(ok .and. root_only_topology%level_count() == 1, &
    "MPI tagged EB root-only topology", comm)
  call initialize_reactive_amr_eb_patch_tree_2d( &
    species, root_state, root_temperature, root_only_topology, &
    tagged_root, ok)
  call assert_all(ok .and. tagged_root%is_valid(), &
    "MPI tagged EB root-only fields", comm)
  tagged_serial = tagged_root
  call initialize_mpi_amr_eb_patch_tree_distribution_2d( &
    root_only_topology, comm, tagged_distribution, ok, 2)
  call assert_all(ok, "MPI tagged EB root-only ownership", comm)
  call initialize_sparse_owned_reactive_amr_eb_patch_tree_2d( &
    tagged_distribution, tagged_root, tagged_sparse, ok)
  call assert_all(ok .and. tagged_sparse%is_valid(tagged_distribution), &
    "MPI tagged EB root-only sparse state", comm)

  tagged_criteria%relative_gradient_threshold = 0.01_dp
  tagged_criteria%absolute_gradient_threshold = 1.0_dp
  tagged_criteria%scale_floor = 1.0_dp
  tagged_criteria%buffer_cells = 0
  tagged_criteria%minimum_patch_cells_x = 4
  tagged_criteria%minimum_patch_cells_y = 4
  tagged_criteria%maximum_patch_gap_cells = 0
  call plan_tagged_reactive_amr_eb_patch_tree_2d( &
    species, tagged_serial, tagged_criteria, 3, 2, &
    build_tagged_regular_geometry, serial_tagged_plans, &
    serial_tagged_cells, ok)
  call assert_all(ok .and. size(serial_tagged_plans) == 2 .and. &
    serial_tagged_cells > 0, "MPI tagged EB serial plan reference", comm)
  call plan_tagged_sparse_owned_reactive_amr_eb_patch_tree_2d( &
    species, tagged_distribution, tagged_sparse, tagged_criteria, 3, 2, &
    build_tagged_regular_geometry, tagged_plans, tagged_cells, ok, &
    local_tagging_evaluations, local_candidate_transfers, &
    local_regrid_restriction_transfers)
  call MPI_Allreduce( &
    local_tagging_evaluations, global_tagging_evaluations, 1, MPI_INTEGER, &
    MPI_SUM, comm, ierr)
  call MPI_Allreduce( &
    local_candidate_transfers, global_candidate_transfers, 1, MPI_INTEGER, &
    MPI_SUM, comm, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. ok .and. &
    tagged_cells == serial_tagged_cells .and. size(tagged_plans) == 2 .and. &
    all([tagged_plans(1)%patch_count(), &
      tagged_plans(2)%patch_count()] == [ &
        serial_tagged_plans(1)%patch_count(), &
        serial_tagged_plans(2)%patch_count()]) .and. &
    global_tagging_evaluations == 1 + tagged_plans(1)%patch_count() .and. &
    global_candidate_transfers >= 0 .and. &
    local_regrid_restriction_transfers == 0, &
    "MPI owner-local arbitrary-depth EB tag planning", comm)

  call regrid_tagged_reactive_amr_eb_patch_tree_2d( &
    species, tagged_serial, tagged_criteria, 3, 2, &
    build_tagged_regular_geometry, ok, topology_changed, &
    serial_tagged_cells)
  call assert_all(ok .and. topology_changed .and. &
    tagged_serial%level_count() == 3, &
    "MPI tagged EB serial first rebuild", comm)
  call regrid_sparse_owned_reactive_amr_eb_patch_tree_2d( &
    species, tagged_distribution, tagged_sparse, tagged_plans, &
    tagged_new_distribution, ok, topology_changed, transferred_cells, &
    local_regrid_restriction_transfers, local_prolongation_transfers, &
    local_overlap_transfers)
  call MPI_Allreduce( &
    local_regrid_restriction_transfers, &
    global_regrid_restriction_transfers, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
  call MPI_Allreduce( &
    local_prolongation_transfers, global_prolongation_transfers, 1, &
    MPI_INTEGER, MPI_SUM, comm, ierr)
  call MPI_Allreduce( &
    local_overlap_transfers, global_overlap_transfers, 1, MPI_INTEGER, &
    MPI_SUM, comm, ierr)
  expected_regrid_edge_transfers = expected_remote_tree_edges( &
    tagged_sparse%topology, tagged_new_distribution)
  call assert_all(ierr == MPI_SUCCESS .and. ok .and. topology_changed .and. &
    transferred_cells == 0 .and. global_overlap_transfers == 0 .and. &
    global_regrid_restriction_transfers == &
      expected_regrid_edge_transfers .and. &
    global_prolongation_transfers == expected_regrid_edge_transfers .and. &
    tagged_sparse%is_valid(tagged_new_distribution), &
    "MPI direct sparse EB tagged-tree creation traffic", comm)
  tagged_distribution = tagged_new_distribution
  call materialize_sparse_owned_reactive_amr_eb_patch_tree_2d( &
    tagged_distribution, tagged_sparse, materialized, ok)
  call assert_all(ok .and. &
    tree_solutions_close(materialized, tagged_serial, 2.0e-11_dp), &
    "MPI direct sparse EB tagged-tree creation parity", comm)

  expanded_tagged_criteria = tagged_criteria
  expanded_tagged_criteria%buffer_cells = 1
  expanded_tagged_criteria%minimum_patch_cells_x = 6
  expanded_tagged_criteria%minimum_patch_cells_y = 6
  call regrid_tagged_reactive_amr_eb_patch_tree_2d( &
    species, tagged_serial, expanded_tagged_criteria, 3, 2, &
    build_tagged_regular_geometry, ok, topology_changed, &
    serial_tagged_cells)
  call assert_all(ok .and. topology_changed, &
    "MPI tagged EB serial expanded rebuild", comm)
  call regrid_tagged_sparse_owned_reactive_amr_eb_patch_tree_2d( &
    species, tagged_distribution, tagged_sparse, expanded_tagged_criteria, &
    3, 2, build_tagged_regular_geometry, tagged_new_distribution, ok, &
    topology_changed, tagged_cells, transferred_cells, &
    local_tagging_evaluations, local_candidate_transfers, &
    local_regrid_restriction_transfers, local_prolongation_transfers, &
    local_overlap_transfers)
  call MPI_Allreduce( &
    local_overlap_transfers, global_overlap_transfers, 1, MPI_INTEGER, &
    MPI_SUM, comm, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. ok .and. topology_changed .and. &
    tagged_cells == serial_tagged_cells .and. transferred_cells > 0 .and. &
    global_overlap_transfers >= 0 .and. &
    tagged_sparse%is_valid(tagged_new_distribution), &
    "MPI owner-local EB tagged-tree overlap rebuild", comm)
  tagged_distribution = tagged_new_distribution
  call materialize_sparse_owned_reactive_amr_eb_patch_tree_2d( &
    tagged_distribution, tagged_sparse, materialized, ok)
  call assert_all(ok .and. &
    tree_solutions_close(materialized, tagged_serial, 3.0e-11_dp), &
    "MPI owner-local EB tagged-tree overlap parity", comm)

  tagged_sparse_snapshot = tagged_sparse
  call regrid_tagged_sparse_owned_reactive_amr_eb_patch_tree_2d( &
    species, tagged_distribution, tagged_sparse, expanded_tagged_criteria, &
    3, 2, build_tagged_regular_geometry, tagged_new_distribution, ok, &
    topology_changed, tagged_cells, transferred_cells)
  call assert_all(ok .and. .not. topology_changed .and. &
    tagged_cells > 0 .and. transferred_cells == 0 .and. &
    sparse_trees_match(tagged_sparse, tagged_sparse_snapshot) .and. &
    mpi_amr_eb_patch_tree_distribution_matches_2d( &
      tagged_new_distribution, tagged_sparse%topology), &
    "MPI unchanged owner-local EB tagged-tree plan", comm)

  rejected_tagged_criteria = expanded_tagged_criteria
  if (nranks > 1) then
    if (rank == 0) rejected_tagged_criteria%relative_gradient_threshold = &
      2.0_dp * rejected_tagged_criteria%relative_gradient_threshold
  else
    rejected_tagged_criteria%relative_gradient_threshold = -1.0_dp
  end if
  tagged_sparse_snapshot = tagged_sparse
  call regrid_tagged_sparse_owned_reactive_amr_eb_patch_tree_2d( &
    species, tagged_distribution, tagged_sparse, rejected_tagged_criteria, &
    3, 2, build_tagged_regular_geometry, tagged_new_distribution, ok, &
    topology_changed, tagged_cells, transferred_cells, &
    local_tagging_evaluations, local_candidate_transfers, &
    local_regrid_restriction_transfers, local_prolongation_transfers, &
    local_overlap_transfers)
  call assert_all(.not. ok .and. .not. topology_changed .and. &
    tagged_cells == 0 .and. transferred_cells == 0 .and. &
    local_tagging_evaluations == 0 .and. local_candidate_transfers == 0 .and. &
    local_regrid_restriction_transfers == 0 .and. &
    local_prolongation_transfers == 0 .and. &
    local_overlap_transfers == 0 .and. &
    sparse_trees_match(tagged_sparse, tagged_sparse_snapshot), &
    "MPI owner-local EB tagged-tree control rollback", comm)

  do level = 1, tagged_serial%level_count()
    do patch = 1, tagged_serial%levels(level)%patch_count()
      tagged_serial%levels(level)%patches(patch)%state = spread(spread( &
        root_state(:, 1, 1), 2, &
        size(tagged_serial%levels(level)%patches(patch)%state, 2)), 3, &
        size(tagged_serial%levels(level)%patches(patch)%state, 3))
      tagged_serial%levels(level)%patches(patch)%temperature = &
        root_temperature(1, 1)
    end do
  end do
  do level = 1, tagged_sparse%level_count()
    do patch = 1, tagged_sparse%levels(level)%patch_count()
      if (.not. tagged_distribution%is_local(level - 1, patch)) cycle
      tagged_sparse%levels(level)%patches(patch)%state = spread(spread( &
        root_state(:, 1, 1), 2, &
        size(tagged_sparse%levels(level)%patches(patch)%state, 2)), 3, &
        size(tagged_sparse%levels(level)%patches(patch)%state, 3))
      tagged_sparse%levels(level)%patches(patch)%temperature = &
        root_temperature(1, 1)
    end do
  end do
  call regrid_tagged_reactive_amr_eb_patch_tree_2d( &
    species, tagged_serial, expanded_tagged_criteria, 3, 2, &
    build_tagged_regular_geometry, ok, topology_changed, &
    serial_tagged_cells)
  call assert_all(ok .and. topology_changed .and. &
    serial_tagged_cells == 0 .and. tagged_serial%level_count() == 1, &
    "MPI tagged EB serial collapse", comm)
  call regrid_tagged_sparse_owned_reactive_amr_eb_patch_tree_2d( &
    species, tagged_distribution, tagged_sparse, expanded_tagged_criteria, &
    3, 2, build_tagged_regular_geometry, tagged_new_distribution, ok, &
    topology_changed, tagged_cells, transferred_cells)
  call assert_all(ok .and. topology_changed .and. tagged_cells == 0 .and. &
    tagged_sparse%level_count() == 1 .and. &
    tagged_sparse%is_valid(tagged_new_distribution), &
    "MPI owner-local EB tagged-tree collapse", comm)
  tagged_distribution = tagged_new_distribution
  call materialize_sparse_owned_reactive_amr_eb_patch_tree_2d( &
    tagged_distribution, tagged_sparse, materialized, ok)
  call assert_all(ok .and. &
    tree_solutions_close(materialized, tagged_serial, 3.0e-11_dp), &
    "MPI owner-local EB tagged-tree collapse parity", comm)

  if (rank == 0) solution%levels(1)%patches(1)%temperature(1, 1) = -1.0_dp
  failed = solution
  call synchronize_owned_reactive_amr_eb_patch_tree_2d( &
    distribution, solution, ok, local_publications)
  call assert_all(.not. ok .and. local_publications == 0 .and. &
    tree_solutions_match(solution, failed), &
    "MPI EB patch-tree collective rejection rollback", comm)
  solution = accepted

  if (nranks > 1) then
    rejected_exponent = modulo(rank, 2)
  else
    rejected_exponent = 3
  end if
  call initialize_mpi_amr_eb_patch_tree_distribution_2d( &
    topology, comm, rejected_distribution, ok, rejected_exponent)
  call assert_all(.not. ok, &
    "MPI EB patch-tree exponent consensus rejection", comm)

  if (rank == 0) &
    write(*, '(a,i0,a)') &
      "pelef_mpi_amr_eb_patch_tree_2d: PASS (", nranks, " ranks)"
  call MPI_Finalize(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI finalization failed"

contains

  subroutine build_regular_geometry( &
      nx, ny, x_lower, x_upper, y_lower, y_upper, geometry, geometry_ok)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: x_lower, x_upper, y_lower, y_upper
    type(eb_geometry_2d), intent(out) :: geometry
    logical, intent(out) :: geometry_ok

    real(dp), allocatable :: level_set(:, :)

    allocate(level_set(0:nx, 0:ny), source=1.0_dp)
    call build_eb_geometry_2d( &
      level_set, x_lower, x_upper, y_lower, y_upper, geometry, geometry_ok)
  end subroutine build_regular_geometry

  subroutine build_tagged_regular_geometry( &
      parent_geometry, i_lower, i_upper, j_lower, j_upper, &
      refinement_ratio, child_geometry, geometry_ok)
    type(eb_geometry_2d), intent(in) :: parent_geometry
    integer, intent(in) :: i_lower, i_upper, j_lower, j_upper
    integer, intent(in) :: refinement_ratio
    type(eb_geometry_2d), intent(out) :: child_geometry
    logical, intent(out) :: geometry_ok

    real(dp) :: x_lower, x_upper, y_lower, y_upper
    integer :: nx, ny

    nx = (i_upper - i_lower + 1) * refinement_ratio
    ny = (j_upper - j_lower + 1) * refinement_ratio
    x_lower = parent_geometry%x_lower + &
      real(i_lower - 1, dp) * parent_geometry%dx
    x_upper = parent_geometry%x_lower + &
      real(i_upper, dp) * parent_geometry%dx
    y_lower = parent_geometry%y_lower + &
      real(j_lower - 1, dp) * parent_geometry%dy
    y_upper = parent_geometry%y_lower + &
      real(j_upper, dp) * parent_geometry%dy
    call build_regular_geometry( &
      nx, ny, x_lower, x_upper, y_lower, y_upper, child_geometry, &
      geometry_ok)
  end subroutine build_tagged_regular_geometry

  integer function expected_remote_tree_edges( &
      tree_topology, tree_distribution) result(count)
    type(amr_eb_patch_tree_topology_2d), intent(in) :: tree_topology
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: &
      tree_distribution

    integer :: child_index, parent_index, relation_index

    count = 0
    do relation_index = 1, size(tree_topology%relations)
      do child_index = 1, tree_topology%relations(relation_index)% &
          child_patch_count()
        parent_index = tree_topology%relations(relation_index)% &
          children(child_index)%parent_patch
        if (tree_distribution%owner_of( &
              relation_index - 1, parent_index) /= &
            tree_distribution%owner_of(relation_index, child_index)) &
          count = count + 1
      end do
    end do
  end function expected_remote_tree_edges

  logical function tree_solutions_match(first, second) result(matches)
    type(reactive_amr_eb_patch_tree_2d), intent(in) :: first, second
    integer :: candidate_level, candidate_patch

    matches = first%level_count() == second%level_count()
    if (.not. matches) return
    do candidate_level = 1, first%level_count()
      matches = first%levels(candidate_level)%patch_count() == &
        second%levels(candidate_level)%patch_count()
      if (.not. matches) return
      do candidate_patch = 1, &
          first%levels(candidate_level)%patch_count()
        matches = all(first%levels(candidate_level)% &
            patches(candidate_patch)%state == second%levels(candidate_level)% &
              patches(candidate_patch)%state) .and. &
          all(first%levels(candidate_level)%patches(candidate_patch)% &
            temperature == second%levels(candidate_level)% &
              patches(candidate_patch)%temperature)
        if (.not. matches) return
      end do
    end do
  end function tree_solutions_match

  logical function tree_solutions_close( &
      first, second, relative_tolerance) result(matches)
    type(reactive_amr_eb_patch_tree_2d), intent(in) :: first, second
    real(dp), intent(in) :: relative_tolerance

    real(dp) :: scale
    integer :: candidate_level, candidate_patch

    matches = relative_tolerance >= 0.0_dp .and. &
      first%level_count() == second%level_count()
    if (.not. matches) return
    do candidate_level = 1, first%level_count()
      matches = first%levels(candidate_level)%patch_count() == &
        second%levels(candidate_level)%patch_count()
      if (.not. matches) return
      do candidate_patch = 1, &
          first%levels(candidate_level)%patch_count()
        scale = max(1.0_dp, maxval(abs(first%levels(candidate_level)% &
          patches(candidate_patch)%state)), &
          maxval(abs(second%levels(candidate_level)% &
            patches(candidate_patch)%state)))
        matches = maxval(abs(first%levels(candidate_level)% &
            patches(candidate_patch)%state - second%levels(candidate_level)% &
              patches(candidate_patch)%state)) <= &
          relative_tolerance * scale
        if (.not. matches) return
        scale = max(1.0_dp, maxval(first%levels(candidate_level)% &
          patches(candidate_patch)%temperature), &
          maxval(second%levels(candidate_level)% &
            patches(candidate_patch)%temperature))
        matches = maxval(abs(first%levels(candidate_level)% &
            patches(candidate_patch)%temperature - &
          second%levels(candidate_level)%patches(candidate_patch)% &
            temperature)) <= relative_tolerance * scale
        if (.not. matches) return
      end do
    end do
  end function tree_solutions_close

  integer function expected_sparse_hydro_transfers( &
      tree_topology, tree_distribution) result(count)
    type(amr_eb_patch_tree_topology_2d), intent(in) :: tree_topology
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: &
      tree_distribution

    integer :: child_index, parent_index, parent_invocations
    integer :: ratio, relation_index

    count = 0
    parent_invocations = 1
    do relation_index = 1, size(tree_topology%relations)
      ratio = tree_topology%relations(relation_index)%refinement_ratio
      do child_index = 1, tree_topology%relations(relation_index)% &
          child_patch_count()
        parent_index = tree_topology%relations(relation_index)% &
          children(child_index)%parent_patch
        if (tree_distribution%owner_of( &
              relation_index - 1, parent_index) /= &
            tree_distribution%owner_of(relation_index, child_index)) &
          count = count + parent_invocations * (ratio + 4)
      end do
      parent_invocations = parent_invocations * ratio
    end do
  end function expected_sparse_hydro_transfers

  logical function sparse_trees_match(first, second) result(matches)
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(in) :: &
      first, second
    integer :: candidate_level, candidate_patch

    matches = first%nvar == second%nvar .and. &
      first%level_count() == second%level_count()
    if (.not. matches) return
    do candidate_level = 1, first%level_count()
      matches = first%levels(candidate_level)%patch_count() == &
        second%levels(candidate_level)%patch_count()
      if (.not. matches) return
      do candidate_patch = 1, &
          first%levels(candidate_level)%patch_count()
        matches = allocated(first%levels(candidate_level)% &
            patches(candidate_patch)%state) .eqv. &
          allocated(second%levels(candidate_level)% &
            patches(candidate_patch)%state)
        matches = matches .and. &
          (allocated(first%levels(candidate_level)% &
            patches(candidate_patch)%temperature) .eqv. &
          allocated(second%levels(candidate_level)% &
            patches(candidate_patch)%temperature))
        if (.not. matches) return
        if (.not. first%levels(candidate_level)% &
            patches(candidate_patch)%has_data()) cycle
        matches = all(first%levels(candidate_level)% &
            patches(candidate_patch)%state == second%levels(candidate_level)% &
              patches(candidate_patch)%state) .and. &
          all(first%levels(candidate_level)%patches(candidate_patch)% &
            temperature == second%levels(candidate_level)% &
              patches(candidate_patch)%temperature)
        if (.not. matches) return
      end do
    end do
  end function sparse_trees_match

  recursive integer function subtree_node_count( &
      tree_topology, tree_level, tree_patch) result(count)
    type(amr_eb_patch_tree_topology_2d), intent(in) :: tree_topology
    integer, intent(in) :: tree_level, tree_patch

    integer :: first_child, last_child, tree_child

    count = 1
    if (tree_level >= tree_topology%level_count()) return
    first_child = tree_topology%relations(tree_level)% &
      child_offsets(tree_patch) + 1
    last_child = tree_topology%relations(tree_level)% &
      child_offsets(tree_patch + 1)
    do tree_child = first_child, last_child
      count = count + subtree_node_count( &
        tree_topology, tree_level + 1, tree_child)
    end do
  end function subtree_node_count

  subroutine assert_all(condition, message, communicator)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    type(MPI_Comm), intent(in) :: communicator

    logical :: accepted_condition
    integer :: status

    call MPI_Allreduce( &
      condition, accepted_condition, 1, MPI_LOGICAL, MPI_LAND, &
      communicator, status)
    if (status /= MPI_SUCCESS .or. .not. accepted_condition) then
      if (rank == 0) write(*, '(a)') trim(message)
      call MPI_Abort(communicator, 1, status)
      error stop trim(message)
    end if
  end subroutine assert_all

end program pelef_mpi_amr_eb_patch_tree_2d
