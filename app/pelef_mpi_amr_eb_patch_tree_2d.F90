program pelef_mpi_amr_eb_patch_tree_2d
  use, intrinsic :: iso_fortran_env, only: int64
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport
  use simulation_config_reactive_2d_mod, only: reactive_2d_config
  use reactive_2d_mod, only: initialize_reactive_2d
  use eb_geometry_2d_mod, only: eb_geometry_2d, build_eb_geometry_2d
  use amr_eb_patch_tree_2d_mod, only: &
    amr_eb_patch_tree_level_plan_2d, amr_eb_patch_tree_topology_2d, &
    initialize_amr_eb_patch_tree_topology_2d
  use amr_eb_patch_tree_reactive_2d_mod, only: &
    reactive_amr_eb_patch_tree_2d, &
    initialize_reactive_amr_eb_patch_tree_2d, &
    compute_reactive_amr_eb_patch_tree_timestep_2d
  use mpi_amr_eb_patch_tree_2d_mod, only: &
    mpi_amr_eb_patch_tree_distribution_2d, &
    mpi_sparse_reactive_amr_eb_patch_tree_2d, &
    initialize_mpi_amr_eb_patch_tree_distribution_2d, &
    mpi_amr_eb_patch_tree_distribution_matches_2d, &
    synchronize_owned_reactive_amr_eb_patch_tree_2d, &
    initialize_sparse_owned_reactive_amr_eb_patch_tree_2d, &
    materialize_sparse_owned_reactive_amr_eb_patch_tree_2d, &
    migrate_sparse_owned_reactive_amr_eb_patch_tree_2d, &
    compute_sparse_owned_reactive_amr_eb_patch_tree_timestep_2d
  implicit none

  type(MPI_Comm) :: comm
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(reactive_2d_config) :: config
  type(eb_geometry_2d) :: root_geometry, level_one_geometry
  type(eb_geometry_2d) :: branch_a_geometry, branch_b_geometry
  type(eb_geometry_2d) :: deep_geometry
  type(amr_eb_patch_tree_level_plan_2d), allocatable :: plans(:)
  type(amr_eb_patch_tree_topology_2d) :: topology
  type(reactive_amr_eb_patch_tree_2d) :: solution, accepted, failed
  type(reactive_amr_eb_patch_tree_2d) :: physical_solution
  type(reactive_amr_eb_patch_tree_2d) :: materialized
  type(mpi_amr_eb_patch_tree_distribution_2d) :: distribution
  type(mpi_amr_eb_patch_tree_distribution_2d) :: unweighted_distribution
  type(mpi_amr_eb_patch_tree_distribution_2d) :: rejected_distribution
  type(mpi_amr_eb_patch_tree_distribution_2d) :: migrated_distribution
  type(mpi_amr_eb_patch_tree_distribution_2d) :: invalid_distribution
  type(mpi_sparse_reactive_amr_eb_patch_tree_2d) :: sparse, sparse_snapshot
  type(mpi_sparse_reactive_amr_eb_patch_tree_2d) :: physical_sparse
  real(dp), allocatable :: root_state(:, :, :), root_temperature(:, :)
  real(dp) :: base_density, expected_state, expected_temperature
  real(dp) :: hydro_cfl, mismatched_hydro_cfl, serial_dt, sparse_dt
  real(dp) :: transport_cfl
  real(dp) :: level_one_dx, level_one_dy
  integer :: ierr, rank, nranks, level, patch
  integer :: rejected_exponent
  integer :: local_publications, global_publications
  integer :: expected_transfers, global_transfers, local_allocated_cells
  integer :: local_nodes, local_transfers, new_owner, old_owner
  integer :: global_timestep_nodes, local_timestep_nodes
  integer(int64) :: unweighted_work, weighted_work
  logical :: ok, local_ok

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI initialization failed"
  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank, ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI rank query failed"
  call MPI_Comm_size(comm, nranks, ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI size query failed"

  call load_h2o2_elementary_thermo(species, ok)
  call assert_all(ok, "MPI EB patch-tree thermodynamic database", comm)
  call load_h2o2_elementary_transport(transport, ok)
  call assert_all(ok, "MPI EB patch-tree transport database", comm)
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
