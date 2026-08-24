program pelef_mpi_amr_patch_1d
  use, intrinsic :: iso_fortran_env, only: int64
  use mpi_f08
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imx
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use amr_hierarchy_1d_mod, only: amr_two_level_hierarchy_1d
  use amr_patch_tree_1d_mod, only: &
    amr_patch_level_plan_1d, amr_patch_tree_hierarchy_1d, &
    amr_patch_tree_level_fields_1d, initialize_patch_tree_1d, &
    prolong_patch_tree_1d, patch_tree_child_geometry_1d
  use amr_patch_tree_reactive_1d_mod, only: &
    amr_patch_tree_reactive_solution_1d, &
    initialize_patch_tree_reactive_1d, advance_patch_tree_chemistry, &
    patch_tree_reactive_timestep_1d, advance_patch_tree_reactive_hydro_1d, &
    advance_patch_tree_transport, advance_patch_tree_reactive_1d, &
    regrid_patch_tree_reactive_1d, &
    regrid_tagged_patch_tree_reactive_1d, &
    patch_tree_reactive_integrals_1d
  use mpi_amr_patch_1d_mod, only: &
    mpi_amr_patch_distribution_1d, mpi_amr_level_halos_1d, &
    initialize_mpi_amr_patch_distribution_1d, &
    synchronize_owned_patch_tree_fields_1d, &
    exchange_owned_adjacent_patch_halos_1d, &
    synchronize_owned_patch_tree_reactive_1d, &
    advance_owned_patch_tree_chemistry_1d, &
    advance_owned_patch_tree_hydro_1d, &
    advance_owned_patch_tree_transport_1d, &
    advance_owned_patch_tree_reactive_1d
  use mpi_amr_sparse_patch_1d_mod, only: &
    mpi_amr_sparse_reactive_solution_1d, &
    mpi_amr_sparse_communication_counts_1d, &
    scatter_owned_patch_tree_reactive_1d, &
    gather_owned_patch_tree_reactive_1d, &
    migrate_owned_patch_tree_reactive_1d, &
    advance_sparse_patch_tree_chemistry_1d, &
    advance_sparse_patch_tree_hydro_1d, &
    advance_sparse_patch_tree_transport_1d, &
    advance_sparse_patch_tree_reactive_1d, &
    regrid_sparse_patch_tree_reactive_1d, &
    regrid_tagged_sparse_patch_tree_reactive_1d
  implicit none

  integer, parameter :: variable_count = 3
  integer, parameter :: halo_width = 4
  real(dp), parameter :: stale_value = -huge(1.0_dp)

  type(amr_patch_tree_hierarchy_1d) :: hierarchy, comparison_hierarchy
  type(amr_patch_level_plan_1d), allocatable :: plans(:)
  type(amr_patch_tree_level_fields_1d), allocatable :: fields(:)
  type(mpi_amr_level_halos_1d), allocatable :: halos(:)
  type(mpi_amr_patch_distribution_1d) :: distribution
  type(mpi_amr_patch_distribution_1d) :: work_distribution
  type(mpi_amr_patch_distribution_1d) :: comparison_distribution
  type(mpi_amr_patch_distribution_1d) :: reactive_distribution
  type(mpi_amr_patch_distribution_1d) :: migrated_distribution
  type(mpi_amr_patch_distribution_1d) :: adjacent_distribution
  type(mpi_amr_patch_distribution_1d) :: regridded_distribution
  type(mpi_amr_patch_distribution_1d) :: tagged_distribution
  type(mpi_amr_patch_distribution_1d) :: tagged_regridded_distribution
  type(mpi_amr_sparse_reactive_solution_1d) :: sparse_reactive
  type(mpi_amr_sparse_reactive_solution_1d) :: migrated_sparse
  type(mpi_amr_sparse_reactive_solution_1d) :: sparse_chemistry
  type(mpi_amr_sparse_reactive_solution_1d) :: sparse_hydro
  type(mpi_amr_sparse_reactive_solution_1d) :: sparse_transport
  type(mpi_amr_sparse_reactive_solution_1d) :: sparse_split
  type(mpi_amr_sparse_reactive_solution_1d) :: sparse_regrid
  type(mpi_amr_sparse_reactive_solution_1d) :: tagged_sparse
  type(mpi_amr_sparse_reactive_solution_1d) :: rejected_sparse
  type(mpi_amr_sparse_reactive_solution_1d) :: sparse_backup
  type(mpi_amr_sparse_communication_counts_1d) :: sparse_communication
  type(amr_patch_tree_reactive_solution_1d) :: initial_reactive
  type(amr_patch_tree_reactive_solution_1d) :: serial_reactive
  type(amr_patch_tree_reactive_solution_1d) :: distributed_reactive
  type(amr_patch_tree_reactive_solution_1d) :: rejected_reactive
  type(amr_patch_tree_reactive_solution_1d) :: rejected_backup
  type(amr_patch_tree_reactive_solution_1d) :: serial_hydro
  type(amr_patch_tree_reactive_solution_1d) :: distributed_hydro
  type(amr_patch_tree_reactive_solution_1d) :: serial_transport
  type(amr_patch_tree_reactive_solution_1d) :: distributed_transport
  type(amr_patch_tree_reactive_solution_1d) :: synchronized_reactive
  type(amr_patch_tree_reactive_solution_1d) :: gathered_reactive
  type(amr_patch_tree_reactive_solution_1d) :: serial_split
  type(amr_patch_tree_reactive_solution_1d) :: distributed_split
  type(amr_patch_tree_reactive_solution_1d) :: serial_regrid
  type(amr_patch_tree_reactive_solution_1d) :: tagged_initial
  type(amr_patch_tree_reactive_solution_1d) :: tagged_serial
  type(amr_patch_tree_reactive_solution_1d) :: adjacent_initial
  type(amr_patch_tree_reactive_solution_1d) :: adjacent_serial
  type(amr_patch_tree_reactive_solution_1d) :: adjacent_distributed
  type(amr_patch_level_plan_1d), allocatable :: reactive_plans(:)
  type(amr_patch_level_plan_1d), allocatable :: adjacent_reactive_plans(:)
  type(amr_patch_level_plan_1d), allocatable :: regrid_reactive_plans(:)
  type(amr_patch_level_plan_1d), allocatable :: invalid_regrid_plans(:)
  type(amr_patch_level_plan_1d), allocatable :: empty_reactive_plans(:)
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(reactive_1d_config) :: reactive_config
  type(reactive_1d_config) :: adjacent_reactive_config
  type(reactive_1d_config) :: transport_config
  type(reactive_1d_config) :: adjacent_transport_config
  type(reactive_1d_config) :: invalid_split_config
  type(reactive_1d_config) :: tagged_reactive_config
  type(reactive_1d_config) :: invalid_tagged_config
  real(dp) :: root(variable_count, 64)
  real(dp), allocatable :: initial_integral(:), final_integral(:)
  integer(int64), allocatable :: cell_distribution_work(:)
  real(dp) :: reactive_difference, conservation_error, hydro_dt, adjacent_dt
  real(dp) :: transport_dt, adjacent_transport_dt
  real(dp) :: split_dt
  logical :: ok, changed
  integer :: ierr, rank, nranks, level, patch, variable, cell
  integer :: parent, child, left_patch, right_patch, layer, cross_rank_faces
  integer :: local_chemistry_advances, global_chemistry_advances
  integer :: expected_patch_advances, corrupt_owner
  integer :: local_hydro_advances, global_hydro_advances
  integer :: expected_hydro_advances, cross_owner_hydro_faces
  integer :: local_transport_advances, global_transport_advances
  integer :: expected_transport_advances
  integer :: local_halo_transfers, global_halo_transfers
  integer :: local_parent_transfers, global_parent_transfers
  integer :: expected_parent_transfers
  integer :: local_parent_state_transfers, global_parent_state_transfers
  integer :: expected_parent_state_transfers
  integer :: root_owner, old_owner, new_owner, owner_changes
  integer :: hierarchy_difference
  integer :: local_sparse_patches, global_sparse_patches
  integer :: local_sparse_cells, global_sparse_cells
  integer :: local_sparse_values, global_sparse_values
  integer :: replicated_value_count
  integer :: transferred_cells, serial_transferred_cells
  integer :: tagged_cells, serial_tagged_cells
  integer :: local_patch_transfers, global_patch_transfers
  integer :: local_sparse_communication(3), global_sparse_communication(3)
  integer :: expected_sparse_communication(3)
  integer :: local_regrid_communication(2), global_regrid_communication(2)
  integer :: expected_regrid_communication(2)
  integer :: local_tagged_communication(4), global_tagged_communication(4)
  integer :: expected_tagged_communication(4)

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Init failed"
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Comm_rank failed"
  call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Comm_size failed"

  call build_test_hierarchy(.false., plans, hierarchy, ok)
  call assert_all(ok, "valid AMR patch tree", rank)
  call initialize_mpi_amr_patch_distribution_1d( &
    hierarchy, MPI_COMM_WORLD, distribution, ok)
  call assert_all(ok, "deterministic MPI AMR patch distribution", rank)
  call assert_all(distribution%rank == rank, "local rank metadata", rank)
  call assert_all(distribution%nranks == nranks, "rank-count metadata", rank)
  call assert_all( &
    sum(distribution%rank_patch_counts) == 9, &
    "every root/fine patch has exactly one owner", rank)
  call assert_all( &
    sum(distribution%rank_cell_counts) == 152, &
    "owned cell work is conserved", rank)
  call assert_all( &
    distribution%subcycle_exponent == 0 .and. &
      sum(distribution%rank_work_counts) == 152_int64, &
    "default distribution retains cell-count weighting", rank)
  call initialize_mpi_amr_patch_distribution_1d( &
    hierarchy, MPI_COMM_WORLD, work_distribution, ok, 2)
  call assert_all(ok .and. work_distribution%subcycle_exponent == 2, &
    "parabolic subcycle-weighted AMR distribution", rank)
  call assert_all(sum(work_distribution%rank_work_counts) == 416_int64, &
    "parabolic AMR work is conserved", rank)
  call effective_distribution_work_counts( &
    hierarchy, distribution, 2, cell_distribution_work, ok)
  call assert_all(ok .and. &
    maxval(work_distribution%rank_work_counts) <= &
      maxval(cell_distribution_work), &
    "subcycle weighting does not increase maximum AMR work", rank)
  if (nranks == 2) then
    call assert_all( &
      maxval(work_distribution%rank_work_counts) < &
        maxval(cell_distribution_work), &
      "subcycle weighting reduces two-rank AMR work imbalance", rank)
  end if
  comparison_distribution = work_distribution
  comparison_distribution%rank_work_counts(1) = &
    comparison_distribution%rank_work_counts(1) + 1_int64
  call assert_all(.not. comparison_distribution%is_valid(), &
    "inconsistent AMR rank work metadata is rejected", rank)
  call initialize_mpi_amr_patch_distribution_1d( &
    hierarchy, MPI_COMM_WORLD, comparison_distribution, ok, 3)
  call assert_all(.not. ok, "invalid AMR work exponent is rejected", rank)

  cross_rank_faces = 0
  do parent = 1, hierarchy%relations(1)%parent_patch_count()
    do child = 1, hierarchy%relations(1)%child_sets(parent)%patch_count() - 1
      if (hierarchy%relations(1)%child_sets(parent)%patches(child)%fine%upper + &
          1 /= hierarchy%relations(1)%child_sets(parent)% &
          patches(child + 1)%fine%lower) cycle
      left_patch = hierarchy%relations(1)%child_index(parent, child)
      right_patch = hierarchy%relations(1)%child_index(parent, child + 1)
      if (distribution%owner_of(1, left_patch) /= &
          distribution%owner_of(1, right_patch)) &
        cross_rank_faces = cross_rank_faces + 1
    end do
  end do
  if (nranks > 1) then
    call assert_all(cross_rank_faces >= 1, &
      "at least one adjacent face crosses ranks", rank)
  end if

  root = 0.0_dp
  call prolong_patch_tree_1d(root, hierarchy, fields, ok)
  call assert_all(ok, "patch-field allocation", rank)
  do level = 0, hierarchy%level_count() - 1
    do patch = 1, hierarchy%level_patch_count(level)
      fields(level + 1)%patches(patch)%values = stale_value
      if (.not. distribution%is_local(level, patch)) cycle
      do cell = 1, size(fields(level + 1)%patches(patch)%values, 2)
        do variable = 1, variable_count
          fields(level + 1)%patches(patch)%values(variable, cell) = &
            expected_value(level, patch, variable, cell)
        end do
      end do
    end do
  end do

  call exchange_owned_adjacent_patch_halos_1d( &
    distribution, hierarchy, fields, halo_width, halos, ok)
  call assert_all(ok, "owner-authoritative adjacent halo exchange", rank)
  do parent = 1, hierarchy%relations(1)%parent_patch_count()
    do child = 1, hierarchy%relations(1)%child_sets(parent)%patch_count() - 1
      if (hierarchy%relations(1)%child_sets(parent)%patches(child)%fine%upper + &
          1 /= hierarchy%relations(1)%child_sets(parent)% &
          patches(child + 1)%fine%lower) cycle
      left_patch = hierarchy%relations(1)%child_index(parent, child)
      right_patch = hierarchy%relations(1)%child_index(parent, child + 1)
      call assert_all(halos(2)%patches(left_patch)%has_right, &
        "left sibling receives a right halo", rank)
      call assert_all(halos(2)%patches(right_patch)%has_left, &
        "right sibling receives a left halo", rank)
      do layer = 1, halo_width
        do variable = 1, variable_count
          call assert_all( &
            halos(2)%patches(left_patch)%right(variable, layer) == &
              expected_value(1, right_patch, variable, layer), &
            "right-source halo value", rank)
          cell = size(fields(2)%patches(left_patch)%values, 2) - layer + 1
          call assert_all( &
            halos(2)%patches(right_patch)%left(variable, layer) == &
              expected_value(1, left_patch, variable, cell), &
            "left-source halo value", rank)
        end do
      end do
    end do
  end do

  call synchronize_owned_patch_tree_fields_1d( &
    distribution, hierarchy, fields, ok)
  call assert_all(ok, "owner-authoritative full-patch synchronization", rank)
  do level = 0, hierarchy%level_count() - 1
    do patch = 1, hierarchy%level_patch_count(level)
      do cell = 1, size(fields(level + 1)%patches(patch)%values, 2)
        do variable = 1, variable_count
          call assert_all( &
            fields(level + 1)%patches(patch)%values(variable, cell) == &
              expected_value(level, patch, variable, cell), &
            "replicated patch equals owner state", rank)
        end do
      end do
    end do
  end do

  if (nranks > 1) then
    call initialize_mpi_amr_patch_distribution_1d( &
      hierarchy, MPI_COMM_WORLD, comparison_distribution, ok, &
      merge(1, 2, rank == nranks - 1))
    call assert_all(.not. ok, &
      "rank-inconsistent AMR work exponent is rejected collectively", rank)
    call build_test_hierarchy( &
      rank == nranks - 1, plans, comparison_hierarchy, ok)
    call assert_all(ok, "comparison AMR patch tree", rank)
    call initialize_mpi_amr_patch_distribution_1d( &
      comparison_hierarchy, MPI_COMM_WORLD, comparison_distribution, ok)
    call assert_all(.not. ok, &
      "rank-inconsistent hierarchy is rejected collectively", rank)
  end if

  call load_h2o2_elementary_thermo(species, ok)
  call assert_all(ok, "reactive AMR thermodynamics", rank)
  call load_h2o2_elementary_mechanism(reactions, ok)
  call assert_all(ok, "reactive AMR chemistry mechanism", rank)
  call load_h2o2_elementary_transport(transport, ok)
  call assert_all(ok, "reactive AMR transport database", rank)
  call configure_reactive_case(reactive_config)
  call build_reactive_plans(reactive_plans)
  call initialize_patch_tree_reactive_1d( &
    species, reactive_config, reactive_plans, initial_reactive, ok)
  call assert_all(ok .and. initial_reactive%is_valid(), &
    "four-level reactive AMR initialization", rank)
  call initialize_mpi_amr_patch_distribution_1d( &
    initial_reactive%hierarchy, MPI_COMM_WORLD, work_distribution, ok)
  call assert_all(ok, "four-level cell-weighted AMR distribution", rank)
  call effective_distribution_work_counts( &
    initial_reactive%hierarchy, work_distribution, 2, &
    cell_distribution_work, ok)
  call assert_all(ok, "four-level cell-distribution work accounting", rank)
  call initialize_mpi_amr_patch_distribution_1d( &
    initial_reactive%hierarchy, MPI_COMM_WORLD, reactive_distribution, ok, 2)
  call assert_all(ok .and. reactive_distribution%subcycle_exponent == 2, &
    "reactive AMR owner distribution", rank)
  call assert_all(maxval(reactive_distribution%rank_work_counts) <= &
    maxval(cell_distribution_work), &
    "four-level work weighting does not increase maximum rank work", rank)
  if (nranks == 2 .or. nranks == 4) then
    call assert_all(maxval(reactive_distribution%rank_work_counts) < &
      maxval(cell_distribution_work), &
      "four-level work weighting reduces rank imbalance", rank)
  end if

  synchronized_reactive = initial_reactive
  root_owner = reactive_distribution%owner_of(0, 1)
  if (rank /= root_owner) then
    synchronized_reactive%level_advances = rank + 1
    synchronized_reactive%transport_level_advances = rank + 2
    synchronized_reactive%time = real(rank + 1, dp)
    synchronized_reactive%steps = rank + 3
    synchronized_reactive%regrid_evaluations = rank + 4
    synchronized_reactive%regrids = rank + 5
    synchronized_reactive%overlap_cells_transferred = rank + 6
  end if
  call synchronize_owned_patch_tree_reactive_1d( &
    reactive_distribution, synchronized_reactive, ok)
  call assert_all(ok .and. reactive_solution_difference( &
    synchronized_reactive, initial_reactive) == 0.0_dp, &
    "root-owner AMR bookkeeping synchronization", rank)

  gathered_reactive = initial_reactive
  call scatter_owned_patch_tree_reactive_1d( &
    reactive_distribution, gathered_reactive, sparse_reactive, ok)
  call assert_all(ok .and. sparse_reactive%is_valid(reactive_distribution), &
    "replicated-to-sparse owner scatter", rank)
  local_sparse_patches = sparse_reactive%local_patch_count()
  local_sparse_cells = sparse_reactive%local_cell_count()
  local_sparse_values = sparse_reactive%local_value_count()
  call assert_all( &
    local_sparse_patches == reactive_distribution%rank_patch_counts(rank + 1), &
    "sparse storage contains only locally owned patches", rank)
  call assert_all( &
    local_sparse_cells == reactive_distribution%rank_cell_counts(rank + 1), &
    "sparse storage contains only locally owned cells", rank)
  call MPI_Allreduce(local_sparse_patches, global_sparse_patches, 1, &
    MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS, "sparse patch-count reduction", rank)
  call MPI_Allreduce(local_sparse_cells, global_sparse_cells, 1, &
    MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS, "sparse cell-count reduction", rank)
  call MPI_Allreduce(local_sparse_values, global_sparse_values, 1, &
    MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS, "sparse value-count reduction", rank)
  replicated_value_count = reactive_solution_value_count(initial_reactive)
  call assert_all( &
    global_sparse_patches == sum(reactive_distribution%rank_patch_counts), &
    "each sparse patch is stored on exactly one rank", rank)
  call assert_all( &
    global_sparse_cells == sum(reactive_distribution%rank_cell_counts), &
    "each sparse cell is stored on exactly one rank", rank)
  call assert_all(global_sparse_values == replicated_value_count, &
    "sparse field storage has no replicated patch payload", rank)

  call poison_reactive_solution(gathered_reactive)
  call gather_owned_patch_tree_reactive_1d( &
    reactive_distribution, sparse_reactive, gathered_reactive, ok)
  call assert_all(ok .and. reactive_solution_difference( &
    gathered_reactive, initial_reactive) == 0.0_dp, &
    "sparse-to-replicated owner gather", rank)

  migrated_distribution = reactive_distribution
  migrated_distribution%rank_patch_counts = 0
  migrated_distribution%rank_cell_counts = 0
  migrated_distribution%rank_work_counts = 0_int64
  owner_changes = 0
  do level = 0, initial_reactive%hierarchy%level_count() - 1
    do patch = 1, initial_reactive%hierarchy%level_patch_count(level)
      old_owner = reactive_distribution%owner_of(level, patch)
      new_owner = mod(old_owner + 1, nranks)
      migrated_distribution%levels(level + 1)%owners(patch) = new_owner
      migrated_distribution%rank_patch_counts(new_owner + 1) = &
        migrated_distribution%rank_patch_counts(new_owner + 1) + 1
      migrated_distribution%rank_cell_counts(new_owner + 1) = &
        migrated_distribution%rank_cell_counts(new_owner + 1) + &
          migrated_distribution%levels(level + 1)%cell_counts(patch)
      migrated_distribution%rank_work_counts(new_owner + 1) = &
        migrated_distribution%rank_work_counts(new_owner + 1) + &
          migrated_distribution%levels(level + 1)%work_counts(patch)
      if (new_owner /= old_owner) owner_changes = owner_changes + 1
    end do
  end do
  call assert_all(migrated_distribution%is_valid(), &
    "rotated sparse owner distribution", rank)
  if (nranks > 1) call assert_all(owner_changes == global_sparse_patches, &
    "every sparse patch changes owner", rank)
  call migrate_owned_patch_tree_reactive_1d( &
    reactive_distribution, migrated_distribution, sparse_reactive, &
    migrated_sparse, ok, local_patch_transfers)
  call assert_all(ok .and. migrated_sparse%is_valid(migrated_distribution), &
    "same-hierarchy sparse owner migration", rank)
  call MPI_Allreduce( &
    local_patch_transfers, global_patch_transfers, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    global_patch_transfers == owner_changes, &
    "one point-to-point payload per changed patch owner", rank)
  call assert_all( &
    migrated_sparse%local_patch_count() == &
      migrated_distribution%rank_patch_counts(rank + 1), &
    "migrated sparse patch ownership", rank)
  call assert_all( &
    migrated_sparse%local_cell_count() == &
      migrated_distribution%rank_cell_counts(rank + 1), &
    "migrated sparse cell ownership", rank)
  local_sparse_values = migrated_sparse%local_value_count()
  call MPI_Allreduce(local_sparse_values, &
    global_sparse_values, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    global_sparse_values == replicated_value_count, &
    "migrated sparse field storage remains unique", rank)
  call poison_reactive_solution(gathered_reactive)
  call gather_owned_patch_tree_reactive_1d( &
    migrated_distribution, migrated_sparse, gathered_reactive, ok)
  call assert_all(ok .and. reactive_solution_difference( &
    gathered_reactive, initial_reactive) == 0.0_dp, &
    "migrated sparse owner gather", rank)

  serial_reactive = initial_reactive
  distributed_reactive = initial_reactive
  allocate(initial_integral( &
    size(initial_reactive%levels(1)%patches(1)%state, 1)))
  allocate(final_integral(size(initial_integral)))
  call patch_tree_reactive_integrals_1d( &
    initial_reactive, initial_integral, ok)
  call assert_all(ok, "initial reactive composite integral", rank)
  call advance_patch_tree_chemistry( &
    species, reactions, reactive_config, 1.0e-10_dp, serial_reactive, ok)
  call assert_all(ok .and. serial_reactive%is_valid(), &
    "serial patch-tree chemistry reference", rank)
  sparse_chemistry = sparse_reactive
  call advance_sparse_patch_tree_chemistry_1d( &
    species, reactions, reactive_config, 1.0e-10_dp, &
    reactive_distribution, sparse_chemistry, ok, local_chemistry_advances)
  call assert_all(ok .and. &
    sparse_chemistry%is_valid(reactive_distribution), &
    "direct sparse patch-tree chemistry", rank)
  call assert_all( &
    local_chemistry_advances == &
      reactive_distribution%rank_patch_counts(rank + 1), &
    "only local sparse patches execute chemistry", rank)
  call MPI_Allreduce( &
    local_chemistry_advances, global_chemistry_advances, 1, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    global_chemistry_advances == &
      sum(reactive_distribution%rank_patch_counts), &
    "every sparse patch advances exactly once globally", rank)
  gathered_reactive = initial_reactive
  call poison_reactive_solution(gathered_reactive)
  call gather_owned_patch_tree_reactive_1d( &
    reactive_distribution, sparse_chemistry, gathered_reactive, ok)
  call assert_all(ok .and. reactive_solution_difference( &
    gathered_reactive, serial_reactive) <= 5.0e-13_dp, &
    "sparse chemistry matches serial patch tree", rank)

  rejected_sparse = sparse_reactive
  corrupt_owner = reactive_distribution%owner_of(3, 1)
  if (rank == corrupt_owner) &
    rejected_sparse%levels(4)%patches(1)%state(irho, 1) = -1.0_dp
  sparse_backup = rejected_sparse
  call advance_sparse_patch_tree_chemistry_1d( &
    species, reactions, reactive_config, 1.0e-10_dp, &
    reactive_distribution, rejected_sparse, ok, local_chemistry_advances)
  call assert_all(.not. ok .and. local_chemistry_advances == 0, &
    "sparse chemistry failure is rejected globally", rank)
  rejected_reactive = initial_reactive
  call gather_owned_patch_tree_reactive_1d( &
    reactive_distribution, rejected_sparse, rejected_reactive, ok)
  call assert_all(ok, "rejected sparse chemistry gather", rank)
  rejected_backup = initial_reactive
  call gather_owned_patch_tree_reactive_1d( &
    reactive_distribution, sparse_backup, rejected_backup, ok)
  call assert_all(ok .and. reactive_solution_difference( &
    rejected_reactive, rejected_backup) == 0.0_dp, &
    "sparse chemistry rollback is exact", rank)

  call advance_owned_patch_tree_chemistry_1d( &
    species, reactions, reactive_config, 1.0e-10_dp, &
    reactive_distribution, distributed_reactive, ok, &
    local_chemistry_advances)
  call assert_all(ok .and. distributed_reactive%is_valid(), &
    "owner-only distributed patch-tree chemistry", rank)
  call assert_all( &
    local_chemistry_advances == &
      reactive_distribution%rank_patch_counts(rank + 1), &
    "only locally owned patches execute chemistry", rank)
  call MPI_Allreduce( &
    local_chemistry_advances, global_chemistry_advances, 1, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS, &
    "owner-only chemistry execution reduction", rank)
  expected_patch_advances = sum(reactive_distribution%rank_patch_counts)
  call assert_all(global_chemistry_advances == expected_patch_advances, &
    "every reactive patch advances exactly once globally", rank)
  reactive_difference = reactive_solution_difference( &
    distributed_reactive, serial_reactive)
  call assert_all(reactive_difference <= 5.0e-13_dp, &
    "distributed chemistry matches serial patch tree", rank)
  call assert_all( &
    reactive_solution_difference(distributed_reactive, initial_reactive) > &
      100.0_dp * epsilon(1.0_dp), &
    "owner-only chemistry changes the reactive state", rank)
  call patch_tree_reactive_integrals_1d( &
    distributed_reactive, final_integral, ok)
  call assert_all(ok, "distributed reactive composite integral", rank)
  conservation_error = maxval(abs( &
    final_integral(1:5) - initial_integral(1:5)) / &
    max(1.0_dp, abs(initial_integral(1:5))))
  call assert_all(conservation_error <= 3.0e-10_dp, &
    "owner-only chemistry conserves mass momentum energy", rank)

  rejected_reactive = initial_reactive
  corrupt_owner = reactive_distribution%owner_of(3, 1)
  if (rank == corrupt_owner) &
    rejected_reactive%levels(4)%patches(1)%state(irho, 1) = -1.0_dp
  call synchronize_owned_patch_tree_reactive_1d( &
    reactive_distribution, rejected_reactive, ok)
  call assert_all(ok, "corrupt owner state synchronization", rank)
  rejected_backup = rejected_reactive
  call advance_owned_patch_tree_chemistry_1d( &
    species, reactions, reactive_config, 1.0e-10_dp, &
    reactive_distribution, rejected_reactive, ok, &
    local_chemistry_advances)
  call assert_all(.not. ok .and. local_chemistry_advances == 0, &
    "owner chemistry failure is rejected globally", rank)
  call assert_all( &
    reactive_solution_difference(rejected_reactive, rejected_backup) == &
      0.0_dp, "global chemistry rollback is exact", rank)

  serial_hydro = initial_reactive
  distributed_hydro = initial_reactive
  call patch_tree_reactive_timestep_1d( &
    species, reactive_config, initial_reactive, hydro_dt, ok)
  call assert_all(ok .and. hydro_dt > 0.0_dp, &
    "four-level reactive hydro timestep", rank)
  hydro_dt = min(0.10_dp * hydro_dt, 2.0e-8_dp)
  call advance_patch_tree_reactive_hydro_1d( &
    species, reactive_config, hydro_dt, serial_hydro, ok)
  call assert_all(ok .and. serial_hydro%is_valid(), &
    "serial four-level hydro reference", rank)
  gathered_reactive = initial_reactive
  call scatter_owned_patch_tree_reactive_1d( &
    reactive_distribution, gathered_reactive, sparse_hydro, ok)
  call assert_all(ok, "four-level sparse hydro owner scatter", rank)
  call advance_sparse_patch_tree_hydro_1d( &
    species, reactive_config, hydro_dt, reactive_distribution, &
    sparse_hydro, ok, local_hydro_advances)
  expected_hydro_advances = expected_owned_hydro_advances( &
    reactive_distribution, initial_reactive%hierarchy, rank)
  call assert_all(ok .and. &
    local_hydro_advances == expected_hydro_advances .and. &
    all(sparse_hydro%level_advances == [1, 4, 12, 16]), &
    "direct sparse four-level hydro accounting", rank)
  call MPI_Allreduce( &
    local_hydro_advances, global_hydro_advances, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. global_hydro_advances == 33, &
    "direct sparse four-level hydro global count", rank)
  gathered_reactive = initial_reactive
  call poison_reactive_solution(gathered_reactive)
  call gather_owned_patch_tree_reactive_1d( &
    reactive_distribution, sparse_hydro, gathered_reactive, ok)
  call assert_all(ok .and. reactive_solution_difference( &
    gathered_reactive, serial_hydro) <= 5.0e-13_dp, &
    "sparse four-level hydro matches serial", rank)
  call advance_owned_patch_tree_hydro_1d( &
    species, reactive_config, hydro_dt, reactive_distribution, &
    distributed_hydro, ok, local_hydro_advances)
  call assert_all(ok .and. distributed_hydro%is_valid(), &
    "owner-only four-level hydro", rank)
  expected_hydro_advances = expected_owned_hydro_advances( &
    reactive_distribution, initial_reactive%hierarchy, rank)
  call assert_all(local_hydro_advances == expected_hydro_advances, &
    "four-level hydro executes on owners only", rank)
  call MPI_Allreduce( &
    local_hydro_advances, global_hydro_advances, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    global_hydro_advances == 33, &
    "four-level hydro global subcycle count", rank)
  call assert_all(all(distributed_hydro%level_advances == [1, 4, 12, 16]), &
    "four-level distributed hydro level accounting", rank)
  reactive_difference = reactive_solution_difference( &
    distributed_hydro, serial_hydro)
  call assert_all(reactive_difference <= 5.0e-13_dp, &
    "distributed four-level hydro matches serial", rank)
  call patch_tree_reactive_integrals_1d( &
    distributed_hydro, final_integral, ok)
  call assert_all(ok, "distributed hydro composite integral", rank)
  conservation_error = maxval(abs( &
    final_integral - initial_integral) / &
    max(1.0_dp, abs(initial_integral)))
  call assert_all(conservation_error <= 3.0e-10_dp, &
    "owner-only four-level hydro conservation", rank)

  rejected_sparse = sparse_reactive
  corrupt_owner = reactive_distribution%owner_of(3, 1)
  if (rank == corrupt_owner) &
    rejected_sparse%levels(4)%patches(1)%state(irho, 1) = -1.0_dp
  sparse_backup = rejected_sparse
  call advance_sparse_patch_tree_hydro_1d( &
    species, reactive_config, hydro_dt, reactive_distribution, &
    rejected_sparse, ok, local_hydro_advances)
  call assert_all(.not. ok .and. local_hydro_advances == 0, &
    "sparse hydro failure is rejected globally", rank)
  rejected_reactive = initial_reactive
  call gather_owned_patch_tree_reactive_1d( &
    reactive_distribution, rejected_sparse, rejected_reactive, ok)
  call assert_all(ok, "rejected sparse hydro gather", rank)
  rejected_backup = initial_reactive
  call gather_owned_patch_tree_reactive_1d( &
    reactive_distribution, sparse_backup, rejected_backup, ok)
  call assert_all(ok .and. reactive_solution_difference( &
    rejected_reactive, rejected_backup) == 0.0_dp, &
    "sparse hydro rollback is exact", rank)

  rejected_reactive = initial_reactive
  corrupt_owner = reactive_distribution%owner_of(3, 1)
  if (rank == corrupt_owner) &
    rejected_reactive%levels(4)%patches(1)%state(irho, 1) = -1.0_dp
  call synchronize_owned_patch_tree_reactive_1d( &
    reactive_distribution, rejected_reactive, ok)
  call assert_all(ok, "corrupt hydro owner state synchronization", rank)
  rejected_backup = rejected_reactive
  call advance_owned_patch_tree_hydro_1d( &
    species, reactive_config, hydro_dt, reactive_distribution, &
    rejected_reactive, ok, local_hydro_advances)
  call assert_all(.not. ok .and. local_hydro_advances == 0, &
    "owner hydro failure is rejected globally", rank)
  call assert_all( &
    reactive_solution_difference(rejected_reactive, rejected_backup) == &
      0.0_dp, "global hydro rollback is exact", rank)

  call configure_transport_case(reactive_config, transport_config)
  serial_transport = initial_reactive
  distributed_transport = initial_reactive
  call patch_tree_reactive_timestep_1d( &
    species, transport_config, initial_reactive, transport_dt, ok, transport)
  call assert_all(ok .and. transport_dt > 0.0_dp, &
    "four-level reactive transport timestep", rank)
  transport_dt = min(transport_dt, 1.0e-10_dp)
  call advance_patch_tree_transport( &
    species, transport, transport_config, transport_dt, serial_transport, ok)
  call assert_all(ok .and. serial_transport%is_valid(), &
    "serial four-level transport reference", rank)
  gathered_reactive = initial_reactive
  call scatter_owned_patch_tree_reactive_1d( &
    reactive_distribution, gathered_reactive, sparse_transport, ok)
  call assert_all(ok, "four-level sparse transport owner scatter", rank)
  call advance_sparse_patch_tree_transport_1d( &
    species, transport, transport_config, transport_dt, &
    reactive_distribution, sparse_transport, ok, local_transport_advances)
  expected_transport_advances = expected_owned_transport_advances( &
    reactive_distribution, initial_reactive%hierarchy, rank)
  call assert_all(ok .and. &
    local_transport_advances == expected_transport_advances .and. &
    all(sparse_transport%transport_level_advances == [1, 8, 48, 128]), &
    "direct sparse four-level transport accounting", rank)
  call MPI_Allreduce( &
    local_transport_advances, global_transport_advances, 1, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    global_transport_advances == 185, &
    "direct sparse four-level transport global count", rank)
  gathered_reactive = initial_reactive
  call poison_reactive_solution(gathered_reactive)
  call gather_owned_patch_tree_reactive_1d( &
    reactive_distribution, sparse_transport, gathered_reactive, ok)
  call assert_all(ok .and. reactive_solution_difference( &
    gathered_reactive, serial_transport) <= 5.0e-13_dp, &
    "sparse four-level transport matches serial", rank)
  call advance_owned_patch_tree_transport_1d( &
    species, transport, transport_config, transport_dt, &
    reactive_distribution, distributed_transport, ok, &
    local_transport_advances)
  call assert_all(ok .and. distributed_transport%is_valid(), &
    "owner-only four-level transport", rank)
  expected_transport_advances = expected_owned_transport_advances( &
    reactive_distribution, initial_reactive%hierarchy, rank)
  call assert_all(local_transport_advances == expected_transport_advances, &
    "four-level transport executes on owners only", rank)
  call MPI_Allreduce( &
    local_transport_advances, global_transport_advances, 1, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    global_transport_advances == 185, &
    "four-level transport global subcycle count", rank)
  call assert_all(all(distributed_transport%transport_level_advances == &
    [1, 8, 48, 128]), &
    "four-level distributed transport level accounting", rank)
  reactive_difference = reactive_solution_difference( &
    distributed_transport, serial_transport)
  call assert_all(reactive_difference <= 5.0e-13_dp, &
    "distributed four-level transport matches serial", rank)
  call assert_all( &
    reactive_solution_difference(distributed_transport, initial_reactive) > &
      100.0_dp * epsilon(1.0_dp), &
    "owner-only transport changes the reactive state", rank)
  call patch_tree_reactive_integrals_1d( &
    distributed_transport, final_integral, ok)
  call assert_all(ok, "distributed transport composite integral", rank)
  conservation_error = maxval(abs( &
    final_integral(1:5) - initial_integral(1:5)) / &
    max(1.0_dp, abs(initial_integral(1:5))))
  call assert_all(conservation_error <= 2.0e-9_dp, &
    "owner-only four-level transport conservation", rank)

  rejected_sparse = sparse_reactive
  corrupt_owner = reactive_distribution%owner_of(3, 1)
  if (rank == corrupt_owner) &
    rejected_sparse%levels(4)%patches(1)%state(irho, 1) = -1.0_dp
  sparse_backup = rejected_sparse
  call advance_sparse_patch_tree_transport_1d( &
    species, transport, transport_config, transport_dt, &
    reactive_distribution, rejected_sparse, ok, local_transport_advances)
  call assert_all(.not. ok .and. local_transport_advances == 0, &
    "sparse transport failure is rejected globally", rank)
  rejected_reactive = initial_reactive
  call gather_owned_patch_tree_reactive_1d( &
    reactive_distribution, rejected_sparse, rejected_reactive, ok)
  call assert_all(ok, "rejected sparse transport gather", rank)
  rejected_backup = initial_reactive
  call gather_owned_patch_tree_reactive_1d( &
    reactive_distribution, sparse_backup, rejected_backup, ok)
  call assert_all(ok .and. reactive_solution_difference( &
    rejected_reactive, rejected_backup) == 0.0_dp, &
    "sparse transport rollback is exact", rank)

  rejected_reactive = initial_reactive
  corrupt_owner = reactive_distribution%owner_of(3, 1)
  if (rank == corrupt_owner) &
    rejected_reactive%levels(4)%patches(1)%state(irho, 1) = -1.0_dp
  call synchronize_owned_patch_tree_reactive_1d( &
    reactive_distribution, rejected_reactive, ok)
  call assert_all(ok, "corrupt transport owner state synchronization", rank)
  rejected_backup = rejected_reactive
  call advance_owned_patch_tree_transport_1d( &
    species, transport, transport_config, transport_dt, &
    reactive_distribution, rejected_reactive, ok, local_transport_advances)
  call assert_all(.not. ok .and. local_transport_advances == 0, &
    "owner transport failure is rejected globally", rank)
  call assert_all( &
    reactive_solution_difference(rejected_reactive, rejected_backup) == &
      0.0_dp, "global transport rollback is exact", rank)

  split_dt = transport_dt
  serial_split = initial_reactive
  distributed_split = initial_reactive
  call advance_patch_tree_reactive_1d( &
    species, reactions, transport_config, split_dt, serial_split, ok, &
    transport)
  call assert_all(ok .and. serial_split%is_valid(), &
    "serial four-level full-physics reference", rank)
  call scatter_owned_patch_tree_reactive_1d( &
    reactive_distribution, initial_reactive, sparse_split, ok)
  call assert_all(ok .and. sparse_split%is_valid(reactive_distribution), &
    "four-level full-physics sparse owner scatter", rank)
  call advance_sparse_patch_tree_reactive_1d( &
    species, reactions, transport_config, split_dt, reactive_distribution, &
    sparse_split, ok, transport, local_chemistry_advances, &
    local_hydro_advances, local_transport_advances)
  call assert_all(ok .and. sparse_split%is_valid(reactive_distribution), &
    "sparse four-level full-physics transaction", rank)
  call assert_all(local_chemistry_advances == 2 * &
    reactive_distribution%rank_patch_counts(rank + 1), &
    "sparse full split chemistry executes on owners only", rank)
  expected_hydro_advances = expected_owned_hydro_advances( &
    reactive_distribution, initial_reactive%hierarchy, rank)
  call assert_all(local_hydro_advances == expected_hydro_advances, &
    "sparse full split hydro executes on owners only", rank)
  expected_transport_advances = expected_owned_transport_advances( &
    reactive_distribution, initial_reactive%hierarchy, rank)
  call assert_all(local_transport_advances == &
    2 * expected_transport_advances, &
    "sparse full split transport executes on owners only", rank)
  call MPI_Allreduce( &
    local_chemistry_advances, global_chemistry_advances, 1, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    global_chemistry_advances == 16, &
    "sparse full split global chemistry call count", rank)
  call MPI_Allreduce( &
    local_hydro_advances, global_hydro_advances, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. global_hydro_advances == 33, &
    "sparse full split global hydro call count", rank)
  call MPI_Allreduce( &
    local_transport_advances, global_transport_advances, 1, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    global_transport_advances == 370, &
    "sparse full split global transport call count", rank)
  call assert_all(all(sparse_split%level_advances == [1, 4, 12, 16]) &
    .and. all(sparse_split%transport_level_advances == &
      [2, 16, 96, 256]) .and. sparse_split%steps == 1 .and. &
    abs(sparse_split%time - split_dt) <= &
      16.0_dp * epsilon(1.0_dp) * split_dt, &
    "sparse full split time and subcycle accounting", rank)
  gathered_reactive = initial_reactive
  call poison_reactive_solution(gathered_reactive)
  call gather_owned_patch_tree_reactive_1d( &
    reactive_distribution, sparse_split, gathered_reactive, ok)
  call assert_all(ok .and. gathered_reactive%is_valid(), &
    "sparse full-physics owner gather", rank)
  reactive_difference = reactive_solution_difference( &
    gathered_reactive, serial_split)
  call assert_all(reactive_difference <= 5.0e-13_dp, &
    "sparse full-physics transaction matches serial", rank)
  call assert_all( &
    reactive_solution_difference(gathered_reactive, initial_reactive) > &
      100.0_dp * epsilon(1.0_dp), &
    "sparse full-physics transaction changes state", rank)
  call patch_tree_reactive_integrals_1d( &
    gathered_reactive, final_integral, ok)
  call assert_all(ok, "sparse full-physics composite integral", rank)
  conservation_error = maxval(abs( &
    final_integral(1:5) - initial_integral(1:5)) / &
    max(1.0_dp, abs(initial_integral(1:5))))
  call assert_all(conservation_error <= 2.0e-9_dp, &
    "sparse full-physics conservation", rank)
  call advance_owned_patch_tree_reactive_1d( &
    species, reactions, transport_config, split_dt, reactive_distribution, &
    distributed_split, ok, transport, local_chemistry_advances, &
    local_hydro_advances, local_transport_advances)
  call assert_all(ok .and. distributed_split%is_valid(), &
    "owner-only four-level full-physics transaction", rank)
  call assert_all(local_chemistry_advances == 2 * &
    reactive_distribution%rank_patch_counts(rank + 1), &
    "full split chemistry executes on owners only", rank)
  expected_hydro_advances = expected_owned_hydro_advances( &
    reactive_distribution, initial_reactive%hierarchy, rank)
  call assert_all(local_hydro_advances == expected_hydro_advances, &
    "full split hydro executes on owners only", rank)
  expected_transport_advances = expected_owned_transport_advances( &
    reactive_distribution, initial_reactive%hierarchy, rank)
  call assert_all(local_transport_advances == &
    2 * expected_transport_advances, &
    "full split transport executes on owners only", rank)
  call MPI_Allreduce( &
    local_chemistry_advances, global_chemistry_advances, 1, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    global_chemistry_advances == 16, &
    "full split global chemistry call count", rank)
  call MPI_Allreduce( &
    local_hydro_advances, global_hydro_advances, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. global_hydro_advances == 33, &
    "full split global hydro call count", rank)
  call MPI_Allreduce( &
    local_transport_advances, global_transport_advances, 1, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    global_transport_advances == 370, &
    "full split global transport call count", rank)
  call assert_all(all(distributed_split%level_advances == [1, 4, 12, 16]) &
    .and. all(distributed_split%transport_level_advances == &
      [2, 16, 96, 256]) .and. distributed_split%steps == 1 .and. &
    abs(distributed_split%time - split_dt) <= &
      16.0_dp * epsilon(1.0_dp) * split_dt, &
    "full split distributed time and subcycle accounting", rank)
  reactive_difference = reactive_solution_difference( &
    distributed_split, serial_split)
  call assert_all(reactive_difference <= 5.0e-13_dp, &
    "distributed full-physics transaction matches serial", rank)
  call assert_all( &
    reactive_solution_difference(distributed_split, initial_reactive) > &
      100.0_dp * epsilon(1.0_dp), &
    "distributed full-physics transaction changes state", rank)
  call patch_tree_reactive_integrals_1d( &
    distributed_split, final_integral, ok)
  call assert_all(ok, "distributed full-physics composite integral", rank)
  conservation_error = maxval(abs( &
    final_integral(1:5) - initial_integral(1:5)) / &
    max(1.0_dp, abs(initial_integral(1:5))))
  call assert_all(conservation_error <= 2.0e-9_dp, &
    "distributed full-physics conservation", rank)

  rejected_reactive = initial_reactive
  rejected_backup = initial_reactive
  call advance_owned_patch_tree_reactive_1d( &
    species, reactions, transport_config, split_dt, reactive_distribution, &
    rejected_reactive, ok, local_chemistry_advances = &
      local_chemistry_advances, &
    local_hydro_advances = local_hydro_advances, &
    local_transport_advances = local_transport_advances)
  call assert_all(.not. ok .and. local_chemistry_advances == 0 .and. &
    local_hydro_advances == 0 .and. local_transport_advances == 0 .and. &
    reactive_solution_difference(rejected_reactive, rejected_backup) == &
      0.0_dp, "missing split transport database is rejected", rank)

  invalid_split_config = transport_config
  invalid_split_config%amr_reconstruction = "invalid"
  rejected_reactive = initial_reactive
  rejected_backup = initial_reactive
  call advance_owned_patch_tree_reactive_1d( &
    species, reactions, invalid_split_config, split_dt, &
    reactive_distribution, rejected_reactive, ok, transport, &
    local_chemistry_advances, local_hydro_advances, &
    local_transport_advances)
  call assert_all(.not. ok .and. local_chemistry_advances == 0 .and. &
    local_hydro_advances == 0 .and. local_transport_advances == 0, &
    "post-prefix split hydro failure is rejected globally", rank)
  call assert_all( &
    reactive_solution_difference(rejected_reactive, rejected_backup) == &
      0.0_dp, "outer full-physics rollback is exact", rank)

  call scatter_owned_patch_tree_reactive_1d( &
    reactive_distribution, initial_reactive, rejected_sparse, ok)
  call assert_all(ok, "missing-transport sparse split scatter", rank)
  sparse_backup = rejected_sparse
  call advance_sparse_patch_tree_reactive_1d( &
    species, reactions, transport_config, split_dt, reactive_distribution, &
    rejected_sparse, ok, local_chemistry_advances = &
      local_chemistry_advances, &
    local_hydro_advances = local_hydro_advances, &
    local_transport_advances = local_transport_advances)
  call assert_all(.not. ok .and. local_chemistry_advances == 0 .and. &
    local_hydro_advances == 0 .and. local_transport_advances == 0, &
    "missing sparse split transport database is rejected", rank)
  rejected_reactive = initial_reactive
  call gather_owned_patch_tree_reactive_1d( &
    reactive_distribution, rejected_sparse, rejected_reactive, ok)
  call assert_all(ok, "missing-transport sparse split gather", rank)
  rejected_backup = initial_reactive
  call gather_owned_patch_tree_reactive_1d( &
    reactive_distribution, sparse_backup, rejected_backup, ok)
  call assert_all(ok .and. reactive_solution_difference( &
    rejected_reactive, rejected_backup) == 0.0_dp, &
    "missing-transport sparse split rollback is exact", rank)

  call scatter_owned_patch_tree_reactive_1d( &
    reactive_distribution, initial_reactive, rejected_sparse, ok)
  call assert_all(ok, "post-prefix sparse split scatter", rank)
  sparse_backup = rejected_sparse
  call advance_sparse_patch_tree_reactive_1d( &
    species, reactions, invalid_split_config, split_dt, &
    reactive_distribution, rejected_sparse, ok, transport, &
    local_chemistry_advances, local_hydro_advances, &
    local_transport_advances)
  call assert_all(.not. ok .and. local_chemistry_advances == 0 .and. &
    local_hydro_advances == 0 .and. local_transport_advances == 0, &
    "post-prefix sparse split hydro failure is rejected globally", rank)
  rejected_reactive = initial_reactive
  call gather_owned_patch_tree_reactive_1d( &
    reactive_distribution, rejected_sparse, rejected_reactive, ok)
  call assert_all(ok, "post-prefix sparse split gather", rank)
  rejected_backup = initial_reactive
  call gather_owned_patch_tree_reactive_1d( &
    reactive_distribution, sparse_backup, rejected_backup, ok)
  call assert_all(ok .and. reactive_solution_difference( &
    rejected_reactive, rejected_backup) == 0.0_dp, &
    "outer sparse full-physics rollback is exact", rank)

  serial_regrid = initial_reactive
  call regrid_patch_tree_reactive_1d( &
    species, transport_config, reactive_plans, serial_regrid, changed, &
    serial_transferred_cells, ok)
  call assert_all(ok .and. .not. changed .and. &
    serial_transferred_cells == 0 .and. &
    serial_regrid%regrid_evaluations == 1, &
    "serial identical regrid reference", rank)
  call scatter_owned_patch_tree_reactive_1d( &
    reactive_distribution, initial_reactive, sparse_regrid, ok)
  call assert_all(ok, "sparse identical regrid scatter", rank)
  call regrid_sparse_patch_tree_reactive_1d( &
    species, transport_config, reactive_plans, reactive_distribution, &
    sparse_regrid, regridded_distribution, changed, transferred_cells, ok, &
    local_regrid_communication(1), local_regrid_communication(2))
  call assert_all(ok .and. .not. changed .and. transferred_cells == 0 .and. &
    all(local_regrid_communication == 0) .and. &
    sparse_regrid%regrid_evaluations == 1 .and. &
    sparse_regrid%regrids == 0 .and. &
    regridded_distribution%subcycle_exponent == 2 .and. &
    sparse_regrid%is_valid(regridded_distribution), &
    "identical sparse regrid is a distributed no-op", rank)
  gathered_reactive = serial_regrid
  call poison_reactive_solution(gathered_reactive)
  call gather_owned_patch_tree_reactive_1d( &
    regridded_distribution, sparse_regrid, gathered_reactive, ok)
  call assert_all(ok .and. &
    reactive_solution_difference(gathered_reactive, serial_regrid) == &
      0.0_dp, "identical sparse regrid matches serial", rank)

  call build_regrid_reactive_plans(regrid_reactive_plans)
  call regrid_patch_tree_reactive_1d( &
    species, transport_config, regrid_reactive_plans, serial_regrid, &
    changed, serial_transferred_cells, ok)
  call assert_all(ok .and. changed .and. serial_transferred_cells > 0 .and. &
    all([ &
      size(serial_regrid%levels(1)%patches), &
      size(serial_regrid%levels(2)%patches), &
      size(serial_regrid%levels(3)%patches), &
      size(serial_regrid%levels(4)%patches)] == [1, 3, 3, 2]), &
    "serial topology-changing regrid reference", rank)
  call regrid_sparse_patch_tree_reactive_1d( &
    species, transport_config, regrid_reactive_plans, &
    regridded_distribution, sparse_regrid, migrated_distribution, changed, &
    transferred_cells, ok, local_regrid_communication(1), &
    local_regrid_communication(2))
  call assert_all(ok .and. &
    migrated_distribution%subcycle_exponent == 2, &
    "direct sparse topology regrid preserves work weighting", rank)
  call MPI_Allreduce( &
    local_regrid_communication, global_regrid_communication, 2, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call expected_direct_sparse_regrid_communication_1d( &
    initial_reactive%hierarchy, regridded_distribution, &
    serial_regrid%hierarchy, migrated_distribution, &
    expected_regrid_communication, ok)
  call assert_all(ierr == MPI_SUCCESS .and. ok .and. all( &
    global_regrid_communication == expected_regrid_communication), &
    "direct sparse regrid communication accounting", rank)
  regridded_distribution = migrated_distribution
  call assert_all(ok .and. changed .and. &
    transferred_cells == serial_transferred_cells .and. &
    sparse_regrid%is_valid(regridded_distribution), &
    "topology-changing sparse regrid accepted", rank)
  call assert_all(sparse_regrid%regrid_evaluations == 2 .and. &
    sparse_regrid%regrids == 1 .and. &
    sparse_regrid%overlap_cells_transferred == transferred_cells, &
    "sparse regrid statistics", rank)
  owner_changes = 0
  do level = 0, min(reactive_distribution%level_count(), &
      regridded_distribution%level_count()) - 1
    do patch = 1, min( &
        reactive_distribution%levels(level + 1)%patch_count(), &
        regridded_distribution%levels(level + 1)%patch_count())
      if (reactive_distribution%owner_of(level, patch) /= &
          regridded_distribution%owner_of(level, patch)) &
        owner_changes = owner_changes + 1
    end do
  end do
  call assert_all(nranks < 4 .or. owner_changes > 0, &
    "topology regrid recomputes patch ownership", rank)
  local_sparse_patches = sparse_regrid%local_patch_count()
  local_sparse_cells = sparse_regrid%local_cell_count()
  local_sparse_values = sparse_regrid%local_value_count()
  call MPI_Allreduce( &
    local_sparse_patches, global_sparse_patches, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call MPI_Allreduce( &
    local_sparse_cells, global_sparse_cells, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call MPI_Allreduce( &
    local_sparse_values, global_sparse_values, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  replicated_value_count = reactive_solution_value_count(serial_regrid)
  call assert_all(ierr == MPI_SUCCESS .and. &
    global_sparse_patches == sum(regridded_distribution%rank_patch_counts) &
    .and. global_sparse_cells == &
      sum(regridded_distribution%rank_cell_counts) .and. &
    global_sparse_values == replicated_value_count, &
    "regridded sparse payload remains globally single-copy", rank)
  gathered_reactive = serial_regrid
  call poison_reactive_solution(gathered_reactive)
  call gather_owned_patch_tree_reactive_1d( &
    regridded_distribution, sparse_regrid, gathered_reactive, ok)
  call assert_all(ok .and. &
    reactive_solution_difference(gathered_reactive, serial_regrid) == &
      0.0_dp, "topology-changing sparse regrid matches serial", rank)
  call patch_tree_reactive_integrals_1d( &
    gathered_reactive, final_integral, ok)
  call assert_all(ok, "sparse regrid composite integral", rank)
  conservation_error = maxval(abs( &
    final_integral - initial_integral) / &
    max(1.0_dp, abs(initial_integral)))
  call assert_all(conservation_error <= 2.0e-9_dp, &
    "topology-changing sparse regrid conservation", rank)

  invalid_regrid_plans = regrid_reactive_plans
  invalid_regrid_plans(2)%patches(3)%parent_patch = 4
  rejected_sparse = sparse_regrid
  sparse_backup = sparse_regrid
  call regrid_sparse_patch_tree_reactive_1d( &
    species, transport_config, invalid_regrid_plans, &
    regridded_distribution, rejected_sparse, migrated_distribution, &
    changed, transferred_cells, ok, local_regrid_communication(1), &
    local_regrid_communication(2))
  call assert_all(.not. ok .and. .not. changed .and. &
    transferred_cells == 0 .and. all(local_regrid_communication == 0) .and. &
    rejected_sparse%is_valid(migrated_distribution), &
    "invalid topology-changing sparse regrid is rejected", rank)
  rejected_reactive = serial_regrid
  call gather_owned_patch_tree_reactive_1d( &
    migrated_distribution, rejected_sparse, rejected_reactive, ok)
  call assert_all(ok, "rejected sparse regrid gather", rank)
  rejected_backup = serial_regrid
  call gather_owned_patch_tree_reactive_1d( &
    regridded_distribution, sparse_backup, rejected_backup, ok)
  call assert_all(ok .and. reactive_solution_difference( &
    rejected_reactive, rejected_backup) == 0.0_dp, &
    "topology-changing sparse regrid rollback is exact", rank)

  call configure_tagged_reactive_case( &
    reactive_config, tagged_reactive_config)
  allocate(empty_reactive_plans(0))
  call initialize_patch_tree_reactive_1d( &
    species, tagged_reactive_config, empty_reactive_plans, tagged_initial, ok)
  call assert_all(ok .and. tagged_initial%level_count() == 1, &
    "root-only tagged sparse regrid initialization", rank)
  tagged_initial%levels(1)%patches(1)%state(imx, &
    1:tagged_reactive_config%nx) = 0.0_dp
  tagged_initial%levels(1)%patches(1)%state(imx, 8) = 10.0_dp
  tagged_initial%levels(1)%patches(1)%state(imx, 24) = -10.0_dp
  call initialize_mpi_amr_patch_distribution_1d( &
    tagged_initial%hierarchy, MPI_COMM_WORLD, tagged_distribution, ok, 1)
  call assert_all(ok .and. tagged_distribution%subcycle_exponent == 1, &
    "tagged sparse owner distribution", rank)
  tagged_serial = tagged_initial
  call regrid_tagged_patch_tree_reactive_1d( &
    species, tagged_reactive_config, tagged_serial, changed, &
    serial_tagged_cells, serial_transferred_cells, ok)
  call assert_all(ok .and. changed .and. serial_tagged_cells > 0 .and. &
    serial_transferred_cells == 0 .and. &
    all([ &
      size(tagged_serial%levels(1)%patches), &
      size(tagged_serial%levels(2)%patches), &
      size(tagged_serial%levels(3)%patches), &
      size(tagged_serial%levels(4)%patches)] == [1, 2, 2, 2]), &
    "serial tag-driven regrid reference", rank)
  call scatter_owned_patch_tree_reactive_1d( &
    tagged_distribution, tagged_initial, tagged_sparse, ok)
  call assert_all(ok, "tag-driven sparse regrid scatter", rank)
  call regrid_tagged_sparse_patch_tree_reactive_1d( &
    species, tagged_reactive_config, tagged_distribution, tagged_sparse, &
    tagged_regridded_distribution, changed, tagged_cells, &
    transferred_cells, ok, local_tagged_communication(1), &
    local_tagged_communication(2), local_tagged_communication(3), &
    local_tagged_communication(4))
  call assert_all(ok .and. &
    tagged_regridded_distribution%subcycle_exponent == 1, &
    "distributed sparse tag planning preserves work weighting", rank)
  call MPI_Allreduce( &
    local_tagged_communication, global_tagged_communication, 4, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call expected_tagged_sparse_communication_1d( &
    tagged_reactive_config, tagged_initial%hierarchy, tagged_distribution, &
    tagged_serial%hierarchy, tagged_regridded_distribution, &
    expected_tagged_communication, ok)
  call assert_all(ierr == MPI_SUCCESS .and. ok .and. all( &
    global_tagged_communication == expected_tagged_communication), &
    "distributed sparse tag communication accounting", rank)
  hierarchy_difference = patch_tree_hierarchy_extent_difference( &
    tagged_sparse%hierarchy, tagged_serial%hierarchy)
  if (rank == 0 .and. hierarchy_difference /= 0) &
    write(*, '(a,1x,i0)') "Tagged hierarchy difference:", &
      hierarchy_difference
  call assert_all(hierarchy_difference == 0, &
    "distributed sparse tag hierarchy matches serial", rank)
  call assert_all(ok .and. changed .and. &
    tagged_cells == serial_tagged_cells .and. &
    transferred_cells == serial_transferred_cells .and. &
    tagged_sparse%is_valid(tagged_regridded_distribution), &
    "tag-driven sparse topology regrid accepted", rank)
  call assert_all(tagged_sparse%regrid_evaluations == 1 .and. &
    tagged_sparse%regrids == 1 .and. &
    tagged_sparse%overlap_cells_transferred == 0, &
    "tag-driven sparse regrid statistics", rank)
  local_sparse_patches = tagged_sparse%local_patch_count()
  local_sparse_cells = tagged_sparse%local_cell_count()
  local_sparse_values = tagged_sparse%local_value_count()
  call MPI_Allreduce( &
    local_sparse_patches, global_sparse_patches, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call MPI_Allreduce( &
    local_sparse_cells, global_sparse_cells, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call MPI_Allreduce( &
    local_sparse_values, global_sparse_values, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  replicated_value_count = reactive_solution_value_count(tagged_serial)
  call assert_all(ierr == MPI_SUCCESS .and. global_sparse_patches == 7 .and. &
    global_sparse_cells == &
      sum(tagged_regridded_distribution%rank_cell_counts) .and. &
    global_sparse_values == replicated_value_count, &
    "tag-driven sparse payload remains globally single-copy", rank)
  gathered_reactive = tagged_serial
  call poison_reactive_solution(gathered_reactive)
  call gather_owned_patch_tree_reactive_1d( &
    tagged_regridded_distribution, tagged_sparse, gathered_reactive, ok)
  reactive_difference = reactive_solution_difference( &
    gathered_reactive, tagged_serial)
  if (rank == 0 .and. reactive_difference /= 0.0_dp) &
    write(*, '(a,1x,es12.5)') "Tagged field difference:", &
      reactive_difference
  if (rank == 0 .and. reactive_difference /= 0.0_dp) &
    call report_reactive_field_differences(gathered_reactive, tagged_serial)
  call assert_all(ok .and. reactive_difference == 0.0_dp, &
    "tag-driven sparse regrid matches serial", rank)
  call patch_tree_reactive_integrals_1d( &
    tagged_initial, initial_integral, ok)
  call assert_all(ok, "pre-tag-driven sparse regrid integral", rank)
  call patch_tree_reactive_integrals_1d( &
    gathered_reactive, final_integral, ok)
  call assert_all(ok, "post-tag-driven sparse regrid integral", rank)
  conservation_error = maxval(abs( &
    final_integral - initial_integral) / &
    max(1.0_dp, abs(initial_integral)))
  call assert_all(conservation_error <= 2.0e-9_dp, &
    "tag-driven sparse regrid conservation", rank)

  call regrid_tagged_patch_tree_reactive_1d( &
    species, tagged_reactive_config, tagged_serial, changed, &
    serial_tagged_cells, serial_transferred_cells, ok)
  call assert_all(ok .and. .not. changed .and. &
    serial_transferred_cells == 0 .and. &
    tagged_serial%regrid_evaluations == 2, &
    "serial unchanged tag-driven regrid reference", rank)
  call regrid_tagged_sparse_patch_tree_reactive_1d( &
    species, tagged_reactive_config, tagged_regridded_distribution, &
    tagged_sparse, migrated_distribution, changed, tagged_cells, &
    transferred_cells, ok, local_tagged_communication(1), &
    local_tagged_communication(2), local_tagged_communication(3), &
    local_tagged_communication(4))
  tagged_regridded_distribution = migrated_distribution
  call assert_all(ok .and. .not. changed .and. &
    tagged_cells == serial_tagged_cells .and. transferred_cells == 0 .and. &
    all(local_tagged_communication(3:4) == 0) .and. &
    tagged_sparse%regrid_evaluations == 2 .and. &
    tagged_sparse%regrids == 1, &
    "unchanged tag-driven sparse regrid is a no-op", rank)
  gathered_reactive = tagged_serial
  call poison_reactive_solution(gathered_reactive)
  call gather_owned_patch_tree_reactive_1d( &
    tagged_regridded_distribution, tagged_sparse, gathered_reactive, ok)
  call assert_all(ok .and. reactive_solution_difference( &
    gathered_reactive, tagged_serial) == 0.0_dp, &
    "unchanged tag-driven sparse regrid matches serial", rank)

  invalid_tagged_config = tagged_reactive_config
  invalid_tagged_config%amr_tag_component = tagged_sparse%nvar + 1
  rejected_sparse = tagged_sparse
  sparse_backup = tagged_sparse
  call regrid_tagged_sparse_patch_tree_reactive_1d( &
    species, invalid_tagged_config, tagged_regridded_distribution, &
    rejected_sparse, migrated_distribution, changed, tagged_cells, &
    transferred_cells, ok, local_tagged_communication(1), &
    local_tagged_communication(2), local_tagged_communication(3), &
    local_tagged_communication(4))
  call assert_all(.not. ok .and. .not. changed .and. tagged_cells == 0 .and. &
    transferred_cells == 0 .and. all(local_tagged_communication == 0) .and. &
    rejected_sparse%is_valid(migrated_distribution), &
    "invalid tag-driven sparse regrid is rejected", rank)
  rejected_reactive = tagged_serial
  call gather_owned_patch_tree_reactive_1d( &
    migrated_distribution, rejected_sparse, rejected_reactive, ok)
  call assert_all(ok, "rejected tag-driven sparse regrid gather", rank)
  rejected_backup = tagged_serial
  call gather_owned_patch_tree_reactive_1d( &
    tagged_regridded_distribution, sparse_backup, rejected_backup, ok)
  call assert_all(ok .and. reactive_solution_difference( &
    rejected_reactive, rejected_backup) == 0.0_dp, &
    "tag-driven sparse regrid rollback is exact", rank)

  adjacent_reactive_config = reactive_config
  adjacent_reactive_config%problem = "entropy_wave"
  adjacent_reactive_config%amr_reconstruction = "ppm"
  adjacent_reactive_config%ppm_contact_steepening = .false.
  adjacent_reactive_config%ppm_shock_flattening = .false.
  adjacent_reactive_config%amr_hybrid_weno = .false.
  call build_adjacent_reactive_plans(adjacent_reactive_plans)
  call initialize_patch_tree_reactive_1d( &
    species, adjacent_reactive_config, adjacent_reactive_plans, &
    adjacent_initial, ok)
  call assert_all(ok .and. adjacent_initial%is_valid(), &
    "adjacent reactive PPM initialization", rank)
  call initialize_mpi_amr_patch_distribution_1d( &
    adjacent_initial%hierarchy, MPI_COMM_WORLD, adjacent_distribution, ok, 2)
  call assert_all(ok .and. adjacent_distribution%subcycle_exponent == 2, &
    "adjacent reactive owner distribution", rank)
  cross_owner_hydro_faces = 0
  do child = 1, adjacent_initial%hierarchy%relations(1)% &
      child_sets(1)%patch_count() - 1
    left_patch = adjacent_initial%hierarchy%relations(1)%child_index(1, child)
    right_patch = adjacent_initial%hierarchy%relations(1)% &
      child_index(1, child + 1)
    if (adjacent_distribution%owner_of(1, left_patch) /= &
        adjacent_distribution%owner_of(1, right_patch)) &
      cross_owner_hydro_faces = cross_owner_hydro_faces + 1
  end do
  if (nranks > 1) call assert_all(cross_owner_hydro_faces >= 1, &
    "adjacent reactive face crosses MPI owners", rank)
  expected_parent_transfers = 0
  root_owner = adjacent_distribution%owner_of(0, 1)
  do child = 1, adjacent_initial%hierarchy%relations(1)% &
      child_sets(1)%patch_count()
    left_patch = adjacent_initial%hierarchy%relations(1)%child_index(1, child)
    if (adjacent_distribution%owner_of(1, left_patch) /= root_owner) &
      expected_parent_transfers = expected_parent_transfers + 1
  end do
  expected_parent_state_transfers = 0
  do new_owner = 0, nranks - 1
    if (new_owner == root_owner) cycle
    do child = 1, adjacent_initial%hierarchy%relations(1)% &
        child_sets(1)%patch_count()
      left_patch = adjacent_initial%hierarchy%relations(1)% &
        child_index(1, child)
      if (adjacent_distribution%owner_of(1, left_patch) == new_owner) then
        expected_parent_state_transfers = &
          expected_parent_state_transfers + 1
        exit
      end if
    end do
  end do
  serial_reactive = adjacent_initial
  call advance_patch_tree_chemistry( &
    species, reactions, adjacent_reactive_config, 1.0e-10_dp, &
    serial_reactive, ok)
  call assert_all(ok, "serial adjacent PPM chemistry reference", rank)
  gathered_reactive = adjacent_initial
  call scatter_owned_patch_tree_reactive_1d( &
    adjacent_distribution, gathered_reactive, sparse_reactive, ok)
  call assert_all(ok, "adjacent PPM sparse owner scatter", rank)
  call advance_sparse_patch_tree_chemistry_1d( &
    species, reactions, adjacent_reactive_config, 1.0e-10_dp, &
    adjacent_distribution, sparse_reactive, ok, local_chemistry_advances, &
    local_halo_transfers, local_parent_transfers, &
    local_parent_state_transfers)
  call assert_all(ok .and. &
    local_chemistry_advances == &
      adjacent_distribution%rank_patch_counts(rank + 1), &
    "adjacent PPM sparse chemistry ownership", rank)
  call MPI_Allreduce( &
    local_halo_transfers, global_halo_transfers, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    global_halo_transfers == cross_owner_hydro_faces, &
    "adjacent sparse chemistry sends one payload per cross-owner face", rank)
  call MPI_Allreduce( &
    local_parent_transfers, global_parent_transfers, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    global_parent_transfers == expected_parent_transfers, &
    "sparse chemistry sends one payload per cross-owner child", rank)
  call MPI_Allreduce( &
    local_parent_state_transfers, global_parent_state_transfers, 1, &
    MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    global_parent_state_transfers == expected_parent_state_transfers, &
    "sparse chemistry sends one parent state per remote child owner", rank)
  gathered_reactive = adjacent_initial
  call poison_reactive_solution(gathered_reactive)
  call gather_owned_patch_tree_reactive_1d( &
    adjacent_distribution, sparse_reactive, gathered_reactive, ok)
  call assert_all(ok .and. reactive_solution_difference( &
    gathered_reactive, serial_reactive) <= 5.0e-13_dp, &
    "cross-owner adjacent sparse chemistry matches serial", rank)
  serial_hydro = adjacent_initial
  distributed_hydro = adjacent_initial
  call patch_tree_reactive_integrals_1d( &
    adjacent_initial, initial_integral, ok)
  call assert_all(ok, "adjacent initial composite integral", rank)
  call patch_tree_reactive_timestep_1d( &
    species, adjacent_reactive_config, adjacent_initial, adjacent_dt, ok)
  call assert_all(ok .and. adjacent_dt > 0.0_dp, &
    "adjacent reactive hydro timestep", rank)
  adjacent_dt = min(0.10_dp * adjacent_dt, 2.0e-8_dp)
  call advance_patch_tree_reactive_hydro_1d( &
    species, adjacent_reactive_config, adjacent_dt, serial_hydro, ok)
  call assert_all(ok, "serial adjacent PPM hydro reference", rank)
  gathered_reactive = adjacent_initial
  call scatter_owned_patch_tree_reactive_1d( &
    adjacent_distribution, gathered_reactive, sparse_hydro, ok)
  call assert_all(ok, "adjacent PPM sparse hydro owner scatter", rank)
  call advance_sparse_patch_tree_hydro_1d( &
    species, adjacent_reactive_config, adjacent_dt, adjacent_distribution, &
    sparse_hydro, ok, local_hydro_advances, sparse_communication)
  expected_hydro_advances = expected_owned_hydro_advances( &
    adjacent_distribution, adjacent_initial%hierarchy, rank)
  call assert_all(ok .and. &
    local_hydro_advances == expected_hydro_advances .and. &
    all(sparse_hydro%level_advances == [1, 12]), &
    "direct sparse adjacent PPM hydro accounting", rank)
  local_sparse_communication = [ &
    sparse_communication%interval_state_transfers, &
    sparse_communication%boundary_flux_transfers, &
    sparse_communication%shared_flux_correction_transfers]
  call MPI_Allreduce( &
    local_sparse_communication, global_sparse_communication, 3, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call expected_sparse_physics_communication_counts( &
    adjacent_distribution, adjacent_initial%hierarchy, .false., &
    expected_sparse_communication)
  call assert_all(ierr == MPI_SUCCESS .and. &
    all(global_sparse_communication == expected_sparse_communication), &
    "direct sparse hydro point-to-point communication accounting", rank)
  gathered_reactive = adjacent_initial
  call poison_reactive_solution(gathered_reactive)
  call gather_owned_patch_tree_reactive_1d( &
    adjacent_distribution, sparse_hydro, gathered_reactive, ok)
  call assert_all(ok .and. reactive_solution_difference( &
    gathered_reactive, serial_hydro) <= 5.0e-13_dp, &
    "cross-owner adjacent sparse PPM hydro matches serial", rank)
  call advance_owned_patch_tree_hydro_1d( &
    species, adjacent_reactive_config, adjacent_dt, adjacent_distribution, &
    distributed_hydro, ok, local_hydro_advances)
  call assert_all(ok .and. distributed_hydro%is_valid(), &
    "owner-only adjacent PPM hydro", rank)
  expected_hydro_advances = expected_owned_hydro_advances( &
    adjacent_distribution, adjacent_initial%hierarchy, rank)
  call assert_all(local_hydro_advances == expected_hydro_advances, &
    "adjacent PPM hydro executes on owners only", rank)
  call MPI_Allreduce( &
    local_hydro_advances, global_hydro_advances, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    global_hydro_advances == 13, &
    "adjacent PPM global subcycle count", rank)
  call assert_all(all(distributed_hydro%level_advances == [1, 12]), &
    "adjacent PPM distributed level accounting", rank)
  reactive_difference = reactive_solution_difference( &
    distributed_hydro, serial_hydro)
  call assert_all(reactive_difference <= 5.0e-13_dp, &
    "cross-owner adjacent PPM matches serial", rank)
  call patch_tree_reactive_integrals_1d( &
    distributed_hydro, final_integral, ok)
  call assert_all(ok, "adjacent distributed composite integral", rank)
  conservation_error = maxval(abs( &
    final_integral - initial_integral) / &
    max(1.0_dp, abs(initial_integral)))
  call assert_all(conservation_error <= 3.0e-10_dp, &
    "cross-owner adjacent PPM conservation", rank)

  call configure_transport_case( &
    adjacent_reactive_config, adjacent_transport_config)
  adjacent_serial = adjacent_initial
  adjacent_distributed = adjacent_initial
  call patch_tree_reactive_timestep_1d( &
    species, adjacent_transport_config, adjacent_initial, &
    adjacent_transport_dt, ok, transport)
  call assert_all(ok .and. adjacent_transport_dt > 0.0_dp, &
    "adjacent reactive transport timestep", rank)
  adjacent_transport_dt = min(adjacent_transport_dt, 1.0e-10_dp)
  call advance_patch_tree_transport( &
    species, transport, adjacent_transport_config, adjacent_transport_dt, &
    adjacent_serial, ok)
  call assert_all(ok, "serial adjacent transport reference", rank)
  gathered_reactive = adjacent_initial
  call scatter_owned_patch_tree_reactive_1d( &
    adjacent_distribution, gathered_reactive, sparse_transport, ok)
  call assert_all(ok, "adjacent sparse transport owner scatter", rank)
  call advance_sparse_patch_tree_transport_1d( &
    species, transport, adjacent_transport_config, adjacent_transport_dt, &
    adjacent_distribution, sparse_transport, ok, local_transport_advances, &
    sparse_communication)
  expected_transport_advances = expected_owned_transport_advances( &
    adjacent_distribution, adjacent_initial%hierarchy, rank)
  call assert_all(ok .and. &
    local_transport_advances == expected_transport_advances .and. &
    all(sparse_transport%transport_level_advances == [1, 24]), &
    "direct sparse adjacent transport accounting", rank)
  local_sparse_communication = [ &
    sparse_communication%interval_state_transfers, &
    sparse_communication%boundary_flux_transfers, &
    sparse_communication%shared_flux_correction_transfers]
  call MPI_Allreduce( &
    local_sparse_communication, global_sparse_communication, 3, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call expected_sparse_physics_communication_counts( &
    adjacent_distribution, adjacent_initial%hierarchy, .true., &
    expected_sparse_communication)
  call assert_all(ierr == MPI_SUCCESS .and. &
    all(global_sparse_communication == expected_sparse_communication), &
    "direct sparse transport point-to-point communication accounting", rank)
  gathered_reactive = adjacent_initial
  call poison_reactive_solution(gathered_reactive)
  call gather_owned_patch_tree_reactive_1d( &
    adjacent_distribution, sparse_transport, gathered_reactive, ok)
  call assert_all(ok .and. reactive_solution_difference( &
    gathered_reactive, adjacent_serial) <= 5.0e-13_dp, &
    "cross-owner adjacent sparse transport matches serial", rank)
  call advance_owned_patch_tree_transport_1d( &
    species, transport, adjacent_transport_config, adjacent_transport_dt, &
    adjacent_distribution, adjacent_distributed, ok, &
    local_transport_advances)
  call assert_all(ok .and. adjacent_distributed%is_valid(), &
    "owner-only adjacent transport", rank)
  expected_transport_advances = expected_owned_transport_advances( &
    adjacent_distribution, adjacent_initial%hierarchy, rank)
  call assert_all(local_transport_advances == expected_transport_advances, &
    "adjacent transport executes on owners only", rank)
  call MPI_Allreduce( &
    local_transport_advances, global_transport_advances, 1, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    global_transport_advances == 25, &
    "adjacent transport global subcycle count", rank)
  call assert_all(all(adjacent_distributed%transport_level_advances == &
    [1, 24]), "adjacent distributed transport level accounting", rank)
  reactive_difference = reactive_solution_difference( &
    adjacent_distributed, adjacent_serial)
  call assert_all(reactive_difference <= 5.0e-13_dp, &
    "cross-owner adjacent transport matches serial", rank)
  call patch_tree_reactive_integrals_1d( &
    adjacent_distributed, final_integral, ok)
  call assert_all(ok, "adjacent transport composite integral", rank)
  conservation_error = maxval(abs( &
    final_integral(1:5) - initial_integral(1:5)) / &
    max(1.0_dp, abs(initial_integral(1:5))))
  call assert_all(conservation_error <= 2.0e-9_dp, &
    "cross-owner adjacent transport conservation", rank)

  if (rank == 0) write(*, '(a,i0,a)') &
    "pelef_mpi_amr_patch_1d: PASS (", nranks, " ranks)"
  call MPI_Finalize(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Finalize failed"

contains

  subroutine effective_distribution_work_counts( &
      local_hierarchy, local_distribution, exponent, counts, local_ok)
    type(amr_patch_tree_hierarchy_1d), intent(in) :: local_hierarchy
    type(mpi_amr_patch_distribution_1d), intent(in) :: local_distribution
    integer, intent(in) :: exponent
    integer(int64), allocatable, intent(out) :: counts(:)
    logical, intent(out) :: local_ok

    integer(int64) :: level_scale, patch_work
    integer :: local_level, local_patch, owner, power, ratio

    local_ok = local_hierarchy%is_valid() .and. &
      local_distribution%is_valid() .and. exponent >= 0 .and. exponent <= 2
    if (.not. local_ok) return
    allocate(counts(local_distribution%nranks))
    counts = 0_int64
    level_scale = 1_int64
    do local_level = 0, local_hierarchy%level_count() - 1
      if (local_level > 0) then
        ratio = local_hierarchy%relations(local_level)%refinement_ratio
        do power = 1, exponent
          if (level_scale > huge(level_scale) / int(ratio, int64)) then
            local_ok = .false.
            return
          end if
          level_scale = level_scale * int(ratio, int64)
        end do
      end if
      do local_patch = 1, local_hierarchy%level_patch_count(local_level)
        if (int(local_distribution%levels(local_level + 1)% &
              cell_counts(local_patch), int64) > &
            huge(patch_work) / level_scale) then
          local_ok = .false.
          return
        end if
        patch_work = int(local_distribution%levels(local_level + 1)% &
          cell_counts(local_patch), int64) * level_scale
        owner = local_distribution%owner_of(local_level, local_patch) + 1
        if (owner < 1 .or. owner > size(counts)) then
          local_ok = .false.
          return
        end if
        if (counts(owner) > huge(patch_work) - patch_work) then
          local_ok = .false.
          return
        end if
        counts(owner) = counts(owner) + patch_work
      end do
    end do
  end subroutine effective_distribution_work_counts

  subroutine build_test_hierarchy( &
      shift_last_upper, local_plans, local_hierarchy, local_ok)
    logical, intent(in) :: shift_last_upper
    type(amr_patch_level_plan_1d), allocatable, intent(out) :: local_plans(:)
    type(amr_patch_tree_hierarchy_1d), intent(out) :: local_hierarchy
    logical, intent(out) :: local_ok

    integer, parameter :: lowers(8) = [5, 9, 17, 25, 33, 37, 49, 53]
    integer :: uppers(8)
    integer :: entry

    uppers = [8, 12, 24, 28, 36, 44, 52, 60]
    if (shift_last_upper) uppers(8) = 59
    allocate(local_plans(1))
    local_plans(1)%refinement_ratio = 2
    allocate(local_plans(1)%patches(8))
    do entry = 1, 8
      local_plans(1)%patches(entry)%parent_patch = 1
      local_plans(1)%patches(entry)%lower = lowers(entry)
      local_plans(1)%patches(entry)%upper = uppers(entry)
    end do
    call initialize_patch_tree_1d( &
      64, 0.0_dp, 1.0_dp, local_plans, local_hierarchy, local_ok)
  end subroutine build_test_hierarchy

  subroutine configure_reactive_case(local_config)
    type(reactive_1d_config), intent(out) :: local_config

    local_config = reactive_1d_config()
    local_config%nx = 32
    local_config%x_lower = 0.0_dp
    local_config%x_upper = 0.012_dp
    local_config%cfl = 0.20_dp
    local_config%problem = "reactive_hotspot"
    local_config%riemann_solver = "rusanov"
    local_config%limiter = "mc"
    local_config%boundary_condition = "periodic"
    local_config%chemistry_enabled = .true.
    local_config%transport_enabled = .false.
    local_config%chemistry_relative_tolerance = 1.0e-8_dp
    local_config%chemistry_absolute_tolerance = 1.0e-14_dp
    local_config%initial_temperature = 1200.0_dp
    local_config%initial_pressure = 101325.0_dp
    local_config%initial_velocity = 0.0_dp
    local_config%hotspot_temperature_rise = 200.0_dp
    local_config%hotspot_center = 0.006_dp
    local_config%hotspot_width = 0.0012_dp
    local_config%amr_enabled = .true.
    local_config%amr_reconstruction = "pcm"
  end subroutine configure_reactive_case

  subroutine configure_transport_case(base_config, local_config)
    type(reactive_1d_config), intent(in) :: base_config
    type(reactive_1d_config), intent(out) :: local_config

    local_config = base_config
    local_config%transport_enabled = .true.
    local_config%viscosity_enabled = .true.
    local_config%thermal_conduction_enabled = .true.
    local_config%species_diffusion_enabled = .true.
    local_config%barodiffusion_enabled = .true.
    local_config%transport_cfl = 0.30_dp
  end subroutine configure_transport_case

  subroutine configure_tagged_reactive_case(base_config, local_config)
    type(reactive_1d_config), intent(in) :: base_config
    type(reactive_1d_config), intent(out) :: local_config

    local_config = base_config
    local_config%nx = 32
    local_config%problem = "entropy_wave"
    local_config%chemistry_enabled = .false.
    local_config%transport_enabled = .false.
    local_config%initial_velocity = 25.0_dp
    local_config%density_wave_amplitude = 0.08_dp
    local_config%amr_multipatch_enabled = .true.
    local_config%amr_max_levels = 4
    local_config%amr_refinement_ratio = 2
    local_config%amr_tag_component = imx
    local_config%amr_relative_gradient_threshold = 0.20_dp
    local_config%amr_absolute_gradient_threshold = 1.0_dp
    local_config%amr_scale_floor = 1.0_dp
    local_config%amr_buffer_cells = 4
    local_config%amr_minimum_patch_cells = 8
    local_config%amr_maximum_patch_gap_cells = 4
  end subroutine configure_tagged_reactive_case

  subroutine build_reactive_plans(local_plans)
    type(amr_patch_level_plan_1d), allocatable, intent(out) :: local_plans(:)

    allocate(local_plans(3))
    local_plans(1)%refinement_ratio = 2
    allocate(local_plans(1)%patches(2))
    local_plans(1)%patches(1)%parent_patch = 1
    local_plans(1)%patches(1)%lower = 4
    local_plans(1)%patches(1)%upper = 11
    local_plans(1)%patches(2)%parent_patch = 1
    local_plans(1)%patches(2)%lower = 20
    local_plans(1)%patches(2)%upper = 27

    local_plans(2)%refinement_ratio = 2
    allocate(local_plans(2)%patches(3))
    local_plans(2)%patches(1)%parent_patch = 1
    local_plans(2)%patches(1)%lower = 3
    local_plans(2)%patches(1)%upper = 8
    local_plans(2)%patches(2)%parent_patch = 1
    local_plans(2)%patches(2)%lower = 11
    local_plans(2)%patches(2)%upper = 14
    local_plans(2)%patches(3)%parent_patch = 2
    local_plans(2)%patches(3)%lower = 5
    local_plans(2)%patches(3)%upper = 12

    local_plans(3)%refinement_ratio = 2
    allocate(local_plans(3)%patches(2))
    local_plans(3)%patches(1)%parent_patch = 1
    local_plans(3)%patches(1)%lower = 3
    local_plans(3)%patches(1)%upper = 10
    local_plans(3)%patches(2)%parent_patch = 3
    local_plans(3)%patches(2)%lower = 4
    local_plans(3)%patches(2)%upper = 13
  end subroutine build_reactive_plans

  subroutine build_regrid_reactive_plans(local_plans)
    type(amr_patch_level_plan_1d), allocatable, intent(out) :: local_plans(:)

    call build_reactive_plans(local_plans)
    deallocate(local_plans(1)%patches)
    allocate(local_plans(1)%patches(3))
    local_plans(1)%patches(1)%parent_patch = 1
    local_plans(1)%patches(1)%lower = 6
    local_plans(1)%patches(1)%upper = 13
    local_plans(1)%patches(2)%parent_patch = 1
    local_plans(1)%patches(2)%lower = 14
    local_plans(1)%patches(2)%upper = 16
    local_plans(1)%patches(3)%parent_patch = 1
    local_plans(1)%patches(3)%lower = 18
    local_plans(1)%patches(3)%upper = 25
    local_plans(2)%patches(3)%parent_patch = 3
  end subroutine build_regrid_reactive_plans

  subroutine build_adjacent_reactive_plans(local_plans)
    type(amr_patch_level_plan_1d), allocatable, intent(out) :: local_plans(:)

    integer :: entry

    allocate(local_plans(1))
    local_plans(1)%refinement_ratio = 2
    allocate(local_plans(1)%patches(6))
    do entry = 1, 6
      local_plans(1)%patches(entry)%parent_patch = 1
      local_plans(1)%patches(entry)%lower = 4 * entry
      local_plans(1)%patches(entry)%upper = 4 * entry + 3
    end do
  end subroutine build_adjacent_reactive_plans

  integer function expected_owned_hydro_advances( &
      local_distribution, local_hierarchy, local_rank) result(count)
    type(mpi_amr_patch_distribution_1d), intent(in) :: local_distribution
    type(amr_patch_tree_hierarchy_1d), intent(in) :: local_hierarchy
    integer, intent(in) :: local_rank

    integer :: local_level, local_patch, multiplier

    count = 0
    multiplier = 1
    do local_level = 0, local_hierarchy%level_count() - 1
      if (local_level > 0) multiplier = multiplier * &
        local_hierarchy%relations(local_level)%refinement_ratio
      do local_patch = 1, local_hierarchy%level_patch_count(local_level)
        if (local_distribution%owner_of(local_level, local_patch) == &
            local_rank) count = count + multiplier
      end do
    end do
  end function expected_owned_hydro_advances

  integer function expected_owned_transport_advances( &
      local_distribution, local_hierarchy, local_rank) result(count)
    type(mpi_amr_patch_distribution_1d), intent(in) :: local_distribution
    type(amr_patch_tree_hierarchy_1d), intent(in) :: local_hierarchy
    integer, intent(in) :: local_rank

    integer :: local_level, local_patch, multiplier, ratio

    count = 0
    multiplier = 1
    do local_level = 0, local_hierarchy%level_count() - 1
      if (local_level > 0) then
        ratio = local_hierarchy%relations(local_level)%refinement_ratio
        multiplier = multiplier * ratio * ratio
      end if
      do local_patch = 1, local_hierarchy%level_patch_count(local_level)
        if (local_distribution%owner_of(local_level, local_patch) == &
            local_rank) count = count + multiplier
      end do
    end do
  end function expected_owned_transport_advances

  subroutine expected_sparse_physics_communication_counts( &
      local_distribution, local_hierarchy, parabolic, counts)
    type(mpi_amr_patch_distribution_1d), intent(in) :: local_distribution
    type(amr_patch_tree_hierarchy_1d), intent(in) :: local_hierarchy
    logical, intent(in) :: parabolic
    integer, intent(out) :: counts(3)

    counts = 0
    call accumulate_expected_sparse_communication( &
      local_distribution, local_hierarchy, 1, 1, 1, parabolic, counts)
  end subroutine expected_sparse_physics_communication_counts

  recursive subroutine accumulate_expected_sparse_communication( &
      local_distribution, local_hierarchy, local_level, local_patch, &
      multiplier, parabolic, counts)
    type(mpi_amr_patch_distribution_1d), intent(in) :: local_distribution
    type(amr_patch_tree_hierarchy_1d), intent(in) :: local_hierarchy
    integer, intent(in) :: local_level, local_patch, multiplier
    logical, intent(in) :: parabolic
    integer, intent(inout) :: counts(3)

    logical, allocatable :: recipients(:)
    integer :: child_count, local_child, child_index, child_owner
    integer :: left_index, right_index, owner, ratio, subcycles
    integer :: child_multiplier, recipient

    if (local_level > size(local_hierarchy%relations)) return
    child_count = local_hierarchy%relations(local_level)% &
      child_sets(local_patch)%patch_count()
    if (child_count == 0) return
    owner = local_distribution%owner_of(local_level - 1, local_patch)
    ratio = local_hierarchy%relations(local_level)%refinement_ratio
    subcycles = ratio
    if (parabolic) subcycles = ratio * ratio

    allocate(recipients(local_distribution%nranks))
    recipients = .false.
    do local_child = 1, child_count
      child_index = local_hierarchy%relations(local_level)% &
        child_index(local_patch, local_child)
      child_owner = local_distribution%owner_of(local_level, child_index)
      recipients(child_owner + 1) = .true.
    end do
    do recipient = 0, local_distribution%nranks - 1
      if (recipient /= owner .and. recipients(recipient + 1)) &
        counts(1) = counts(1) + multiplier
    end do

    do local_child = 1, child_count - 1
      if (local_hierarchy%relations(local_level)% &
            child_sets(local_patch)%patches(local_child)% &
              fine_coarse_upper + 1 /= &
          local_hierarchy%relations(local_level)% &
            child_sets(local_patch)%patches(local_child + 1)% &
              fine_coarse_lower) cycle
      left_index = local_hierarchy%relations(local_level)% &
        child_index(local_patch, local_child)
      right_index = local_hierarchy%relations(local_level)% &
        child_index(local_patch, local_child + 1)
      if (local_distribution%owner_of(local_level, left_index) /= owner) &
        counts(3) = counts(3) + multiplier * subcycles
      if (local_distribution%owner_of(local_level, right_index) /= owner) &
        counts(3) = counts(3) + multiplier * subcycles
    end do

    child_multiplier = multiplier * subcycles
    do local_child = 1, child_count
      child_index = local_hierarchy%relations(local_level)% &
        child_index(local_patch, local_child)
      child_owner = local_distribution%owner_of(local_level, child_index)
      if (child_owner /= owner) &
        counts(2) = counts(2) + child_multiplier
      call accumulate_expected_sparse_communication( &
        local_distribution, local_hierarchy, local_level + 1, child_index, &
        child_multiplier, parabolic, counts)
    end do
  end subroutine accumulate_expected_sparse_communication

  subroutine expected_tagged_sparse_communication_1d( &
      local_config, old_hierarchy, old_distribution, new_hierarchy, &
      new_distribution, counts, local_ok)
    type(reactive_1d_config), intent(in) :: local_config
    type(amr_patch_tree_hierarchy_1d), intent(in) :: old_hierarchy
    type(mpi_amr_patch_distribution_1d), intent(in) :: old_distribution
    type(amr_patch_tree_hierarchy_1d), intent(in) :: new_hierarchy
    type(mpi_amr_patch_distribution_1d), intent(in) :: new_distribution
    integer, intent(out) :: counts(4)
    logical, intent(out) :: local_ok

    type(amr_patch_tree_hierarchy_1d) :: prefix_hierarchy
    type(mpi_amr_patch_distribution_1d) :: prefix_distribution
    integer :: direct_counts(2)
    integer :: relation, depth, final_relations, evaluated_relations

    counts = 0
    local_ok = old_hierarchy%is_valid() .and. new_hierarchy%is_valid() .and. &
      local_config%amr_max_levels >= 2
    if (.not. local_ok) return
    final_relations = size(new_hierarchy%relations)
    evaluated_relations = min( &
      local_config%amr_max_levels - 1, final_relations + 1)
    do relation = 1, evaluated_relations
      counts(1) = counts(1) + &
        new_hierarchy%level_patch_count(relation - 1)
    end do

    do depth = 1, final_relations
      prefix_hierarchy%base_cells = new_hierarchy%base_cells
      prefix_hierarchy%x_lower = new_hierarchy%x_lower
      prefix_hierarchy%x_upper = new_hierarchy%x_upper
      if (allocated(prefix_hierarchy%relations)) &
        deallocate(prefix_hierarchy%relations)
      allocate(prefix_hierarchy%relations(depth))
      prefix_hierarchy%relations = new_hierarchy%relations(1:depth)
      call initialize_mpi_amr_patch_distribution_1d( &
        prefix_hierarchy, MPI_COMM_WORLD, prefix_distribution, local_ok, &
        new_distribution%subcycle_exponent)
      if (.not. local_ok) return
      counts(2) = counts(2) + &
        expected_sparse_prolongation_transfers_1d( &
          prefix_hierarchy, prefix_distribution)
    end do

    call expected_direct_sparse_regrid_communication_1d( &
      old_hierarchy, old_distribution, new_hierarchy, new_distribution, &
      direct_counts, local_ok)
    if (.not. local_ok) return
    counts(3:4) = direct_counts
  end subroutine expected_tagged_sparse_communication_1d

  integer function expected_sparse_prolongation_transfers_1d( &
      local_hierarchy, local_distribution) result(count)
    type(amr_patch_tree_hierarchy_1d), intent(in) :: local_hierarchy
    type(mpi_amr_patch_distribution_1d), intent(in) :: local_distribution

    integer :: relation, parent, local_child, child_index

    count = 0
    do relation = 1, size(local_hierarchy%relations)
      do parent = 1, local_hierarchy%relations(relation)%parent_patch_count()
        do local_child = 1, local_hierarchy%relations(relation)% &
            child_sets(parent)%patch_count()
          child_index = local_hierarchy%relations(relation)% &
            child_index(parent, local_child)
          if (local_distribution%owner_of(relation - 1, parent) /= &
              local_distribution%owner_of(relation, child_index)) &
            count = count + 1
        end do
      end do
    end do
  end function expected_sparse_prolongation_transfers_1d

  subroutine expected_direct_sparse_regrid_communication_1d( &
      old_hierarchy, old_distribution, new_hierarchy, new_distribution, &
      counts, local_ok)
    type(amr_patch_tree_hierarchy_1d), intent(in) :: old_hierarchy
    type(mpi_amr_patch_distribution_1d), intent(in) :: old_distribution
    type(amr_patch_tree_hierarchy_1d), intent(in) :: new_hierarchy
    type(mpi_amr_patch_distribution_1d), intent(in) :: new_distribution
    integer, intent(out) :: counts(2)
    logical, intent(out) :: local_ok

    type(amr_two_level_hierarchy_1d) :: old_geometry, new_geometry
    real(dp) :: old_lower, old_upper, new_lower, new_upper
    real(dp) :: overlap_lower, overlap_upper, old_dx, new_dx, tolerance
    integer :: relation, parent, local_child, child_index
    integer :: local_level, old_patch, new_patch, common_levels
    integer :: old_owner, new_owner

    counts = 0
    local_ok = old_hierarchy%is_valid() .and. new_hierarchy%is_valid()
    if (.not. local_ok) return
    do relation = 1, size(new_hierarchy%relations)
      do parent = 1, new_hierarchy%relations(relation)%parent_patch_count()
        do local_child = 1, new_hierarchy%relations(relation)% &
            child_sets(parent)%patch_count()
          child_index = new_hierarchy%relations(relation)% &
            child_index(parent, local_child)
          if (new_distribution%owner_of(relation - 1, parent) /= &
              new_distribution%owner_of(relation, child_index)) &
            counts(1) = counts(1) + 1
        end do
      end do
    end do

    common_levels = min(old_hierarchy%level_count(), &
      new_hierarchy%level_count())
    do local_level = 2, common_levels
      old_dx = old_hierarchy%level_dx(local_level - 1)
      new_dx = new_hierarchy%level_dx(local_level - 1)
      tolerance = 128.0_dp * epsilon(1.0_dp) * &
        max(1.0_dp, abs(old_dx), abs(new_dx))
      if (abs(old_dx - new_dx) > tolerance) cycle
      do old_patch = 1, old_hierarchy%level_patch_count(local_level - 1)
        call patch_tree_child_geometry_1d( &
          old_hierarchy%relations(local_level - 1), old_patch, &
          old_geometry, local_ok)
        if (.not. local_ok) return
        call application_patch_physical_bounds_1d( &
          old_geometry, old_lower, old_upper)
        old_owner = old_distribution%owner_of(local_level - 1, old_patch)
        do new_patch = 1, new_hierarchy%level_patch_count(local_level - 1)
          call patch_tree_child_geometry_1d( &
            new_hierarchy%relations(local_level - 1), new_patch, &
            new_geometry, local_ok)
          if (.not. local_ok) return
          call application_patch_physical_bounds_1d( &
            new_geometry, new_lower, new_upper)
          overlap_lower = max(old_lower, new_lower)
          overlap_upper = min(old_upper, new_upper)
          if (overlap_upper <= overlap_lower + tolerance) cycle
          new_owner = new_distribution%owner_of(local_level - 1, new_patch)
          if (old_owner /= new_owner) counts(2) = counts(2) + 1
        end do
      end do
    end do
    local_ok = .true.
  end subroutine expected_direct_sparse_regrid_communication_1d

  pure subroutine application_patch_physical_bounds_1d( &
      geometry, lower, upper)
    type(amr_two_level_hierarchy_1d), intent(in) :: geometry
    real(dp), intent(out) :: lower, upper

    lower = geometry%x_lower + &
      real(geometry%fine_coarse_lower - 1, dp) * geometry%coarse_dx
    upper = geometry%x_lower + &
      real(geometry%fine_coarse_upper, dp) * geometry%coarse_dx
  end subroutine application_patch_physical_bounds_1d

  integer function patch_tree_hierarchy_extent_difference(first, second) &
      result(difference)
    type(amr_patch_tree_hierarchy_1d), intent(in) :: first, second

    integer :: relation, parent, local_child

    difference = huge(0)
    if (first%base_cells /= second%base_cells) return
    if (size(first%relations) /= size(second%relations)) return
    difference = 0
    do relation = 1, size(first%relations)
      if (first%relations(relation)%refinement_ratio /= &
          second%relations(relation)%refinement_ratio) then
        difference = huge(0)
        return
      end if
      if (first%relations(relation)%parent_patch_count() /= &
          second%relations(relation)%parent_patch_count()) then
        difference = huge(0)
        return
      end if
      do parent = 1, first%relations(relation)%parent_patch_count()
        if (first%relations(relation)%child_sets(parent)%patch_count() /= &
            second%relations(relation)%child_sets(parent)%patch_count()) then
          difference = huge(0)
          return
        end if
        do local_child = 1, first%relations(relation)% &
            child_sets(parent)%patch_count()
          difference = max(difference, abs( &
            first%relations(relation)%child_sets(parent)% &
              patches(local_child)%fine_coarse_lower - &
            second%relations(relation)%child_sets(parent)% &
              patches(local_child)%fine_coarse_lower))
          difference = max(difference, abs( &
            first%relations(relation)%child_sets(parent)% &
              patches(local_child)%fine_coarse_upper - &
            second%relations(relation)%child_sets(parent)% &
              patches(local_child)%fine_coarse_upper))
        end do
      end do
    end do
  end function patch_tree_hierarchy_extent_difference

  subroutine report_reactive_field_differences(first, second)
    type(amr_patch_tree_reactive_solution_1d), intent(in) :: first, second

    real(dp) :: state_error, temperature_error
    real(dp) :: ghost_state_error, ghost_temperature_error
    integer :: local_level, local_patch

    state_error = 0.0_dp
    temperature_error = 0.0_dp
    ghost_state_error = 0.0_dp
    ghost_temperature_error = 0.0_dp
    do local_level = 1, first%level_count()
      do local_patch = 1, size(first%levels(local_level)%patches)
        state_error = max(state_error, maxval(abs( &
          first%levels(local_level)%patches(local_patch)%state - &
          second%levels(local_level)%patches(local_patch)%state)))
        temperature_error = max(temperature_error, maxval(abs( &
          first%levels(local_level)%patches(local_patch)%temperature - &
          second%levels(local_level)%patches(local_patch)%temperature)))
        ghost_state_error = max(ghost_state_error, maxval(abs( &
          first%levels(local_level)%patches(local_patch)%left_ghost_state - &
          second%levels(local_level)%patches(local_patch)%left_ghost_state)))
        ghost_state_error = max(ghost_state_error, maxval(abs( &
          first%levels(local_level)%patches(local_patch)%right_ghost_state - &
          second%levels(local_level)%patches(local_patch)%right_ghost_state)))
        ghost_temperature_error = max(ghost_temperature_error, maxval(abs( &
          first%levels(local_level)%patches(local_patch)% &
            left_ghost_temperature - &
          second%levels(local_level)%patches(local_patch)% &
            left_ghost_temperature)))
        ghost_temperature_error = max(ghost_temperature_error, maxval(abs( &
          first%levels(local_level)%patches(local_patch)% &
            right_ghost_temperature - &
          second%levels(local_level)%patches(local_patch)% &
            right_ghost_temperature)))
      end do
    end do
    write(*, '(a,4(1x,es12.5))') "Tagged component differences:", &
      state_error, temperature_error, ghost_state_error, &
      ghost_temperature_error
  end subroutine report_reactive_field_differences

  real(dp) function reactive_solution_difference(first, second) result(error)
    type(amr_patch_tree_reactive_solution_1d), intent(in) :: first, second

    integer :: local_level, local_patch

    error = huge(1.0_dp)
    if (first%level_count() /= second%level_count()) return
    if (size(first%level_advances) /= size(second%level_advances)) return
    if (size(first%transport_level_advances) /= &
        size(second%transport_level_advances)) return
    error = abs(first%time - second%time) / max(1.0_dp, abs(second%time))
    error = max(error, real(abs(first%steps - second%steps), dp))
    error = max(error, real(maxval(abs( &
      first%level_advances - second%level_advances)), dp))
    error = max(error, real(maxval(abs( &
      first%transport_level_advances - &
      second%transport_level_advances)), dp))
    error = max(error, real(abs( &
      first%regrid_evaluations - second%regrid_evaluations), dp))
    error = max(error, real(abs(first%regrids - second%regrids), dp))
    error = max(error, real(abs( &
      first%overlap_cells_transferred - &
      second%overlap_cells_transferred), dp))
    do local_level = 1, first%level_count()
      if (size(first%levels(local_level)%patches) /= &
          size(second%levels(local_level)%patches)) then
        error = huge(1.0_dp)
        return
      end if
      do local_patch = 1, size(first%levels(local_level)%patches)
        error = max(error, maxval(abs( &
          first%levels(local_level)%patches(local_patch)%state - &
          second%levels(local_level)%patches(local_patch)%state) / &
          max(1.0_dp, abs(second%levels(local_level)% &
            patches(local_patch)%state))))
        error = max(error, maxval(abs( &
          first%levels(local_level)%patches(local_patch)%temperature - &
          second%levels(local_level)%patches(local_patch)%temperature) / &
          max(1.0_dp, abs(second%levels(local_level)% &
            patches(local_patch)%temperature))))
        error = max(error, maxval(abs( &
          first%levels(local_level)%patches(local_patch)%left_ghost_state - &
          second%levels(local_level)%patches(local_patch)%left_ghost_state) / &
          max(1.0_dp, abs(second%levels(local_level)% &
            patches(local_patch)%left_ghost_state))))
        error = max(error, maxval(abs( &
          first%levels(local_level)%patches(local_patch)%right_ghost_state - &
          second%levels(local_level)%patches(local_patch)% &
            right_ghost_state) / max(1.0_dp, abs(second%levels(local_level)% &
              patches(local_patch)%right_ghost_state))))
        error = max(error, maxval(abs( &
          first%levels(local_level)%patches(local_patch)% &
            left_ghost_temperature - second%levels(local_level)% &
            patches(local_patch)%left_ghost_temperature) / &
          max(1.0_dp, abs(second%levels(local_level)%patches(local_patch)% &
            left_ghost_temperature))))
        error = max(error, maxval(abs( &
          first%levels(local_level)%patches(local_patch)% &
            right_ghost_temperature - second%levels(local_level)% &
            patches(local_patch)%right_ghost_temperature) / &
          max(1.0_dp, abs(second%levels(local_level)%patches(local_patch)% &
            right_ghost_temperature))))
      end do
    end do
  end function reactive_solution_difference

  integer function reactive_solution_value_count(solution) result(count)
    type(amr_patch_tree_reactive_solution_1d), intent(in) :: solution

    integer :: local_level, local_patch

    count = 0
    do local_level = 1, solution%level_count()
      do local_patch = 1, size(solution%levels(local_level)%patches)
        count = count + &
          size(solution%levels(local_level)%patches(local_patch)%state) + &
          size(solution%levels(local_level)%patches(local_patch)%temperature) + &
          size(solution%levels(local_level)%patches(local_patch)% &
            left_ghost_state) + &
          size(solution%levels(local_level)%patches(local_patch)% &
            right_ghost_state) + &
          size(solution%levels(local_level)%patches(local_patch)% &
            left_ghost_temperature) + &
          size(solution%levels(local_level)%patches(local_patch)% &
            right_ghost_temperature)
      end do
    end do
  end function reactive_solution_value_count

  subroutine poison_reactive_solution(solution)
    type(amr_patch_tree_reactive_solution_1d), intent(inout) :: solution

    integer :: local_level, local_patch

    solution%level_advances = 101
    solution%transport_level_advances = 102
    solution%time = 103.0_dp
    solution%steps = 104
    solution%regrid_evaluations = 105
    solution%regrids = 106
    solution%overlap_cells_transferred = 107
    do local_level = 1, solution%level_count()
      do local_patch = 1, size(solution%levels(local_level)%patches)
        solution%levels(local_level)%patches(local_patch)%state = stale_value
        solution%levels(local_level)%patches(local_patch)%temperature = &
          stale_value
        solution%levels(local_level)%patches(local_patch)%left_ghost_state = &
          stale_value
        solution%levels(local_level)%patches(local_patch)%right_ghost_state = &
          stale_value
        solution%levels(local_level)%patches(local_patch)% &
          left_ghost_temperature = stale_value
        solution%levels(local_level)%patches(local_patch)% &
          right_ghost_temperature = stale_value
      end do
    end do
  end subroutine poison_reactive_solution

  pure real(dp) function expected_value( &
      local_level, local_patch, local_variable, local_cell) result(value)
    integer, intent(in) :: local_level, local_patch, local_variable, local_cell

    value = real( &
      100000 * local_level + 10000 * local_patch + &
      100 * local_variable + local_cell, dp)
  end function expected_value

  subroutine assert_all(condition, message, local_rank)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    integer, intent(in) :: local_rank

    logical :: global_condition
    integer :: local_ierr

    call MPI_Allreduce(condition, global_condition, 1, MPI_LOGICAL, MPI_LAND, &
      MPI_COMM_WORLD, local_ierr)
    if (local_ierr == MPI_SUCCESS .and. global_condition) return
    if (local_rank == 0) write(*, '(a)') "FAIL: " // trim(message)
    call MPI_Abort(MPI_COMM_WORLD, 1, local_ierr)
  end subroutine assert_all

end program pelef_mpi_amr_patch_1d
