program test_mpi_amr_eb_selected_checkpoint_presence
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use amr_eb_patch_tree_reactive_2d_mod, only: &
    reactive_amr_eb_patch_tree_checkpoint_fingerprint_2d
  use mpi_amr_eb_patch_tree_2d_mod, only: &
    mpi_amr_eb_patch_tree_distribution_2d, &
    mpi_sparse_reactive_amr_eb_patch_tree_2d
  use mpi_amr_eb_patch_tree_io_2d_mod, only: &
    read_sparse_owned_reactive_amr_eb_patch_tree_2d_checkpoint, &
    write_sparse_owned_reactive_amr_eb_patch_tree_2d_checkpoint
  implicit none

  character(len=*), parameter :: bundle_sha256 = &
    "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef"

  type(nasa7_species) :: species(1)
  type(reactive_amr_eb_patch_tree_checkpoint_fingerprint_2d) :: fingerprint
  type(mpi_amr_eb_patch_tree_distribution_2d) :: distribution
  type(mpi_sparse_reactive_amr_eb_patch_tree_2d) :: sparse
  real(dp), allocatable :: initial_integrals(:)
  real(dp) :: composition(1), minimum_dt, time
  integer :: chemistry_advances(1), hydro_advances(1), transport_advances(1)
  integer :: ierr, nranks, rank, regrids, steps, transfers
  logical :: all_rejected, any_file_exists, local_file_exists
  logical :: local_rejected, ok

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Init failed"
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Comm_rank failed"
  call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)
  if (ierr /= MPI_SUCCESS .or. nranks /= 2) &
    error stop "This test requires exactly two MPI ranks"

  species(1)%name = "X"
  composition = 1.0_dp
  fingerprint = reactive_amr_eb_patch_tree_checkpoint_fingerprint_2d()
  if (rank == 0) then
    call read_sparse_owned_reactive_amr_eb_patch_tree_2d_checkpoint( &
      "selected_checkpoint_presence_should_not_exist.dat", species, &
      MPI_COMM_WORLD, 0, 4, 1, distribution, sparse, time, steps, regrids, &
      minimum_dt, ok, fingerprint=fingerprint, &
      bundle_sha256=bundle_sha256, chemistry_integrator="explicit", &
      base_mole_fractions=composition)
  else
    call read_sparse_owned_reactive_amr_eb_patch_tree_2d_checkpoint( &
      "selected_checkpoint_presence_should_not_exist.dat", species, &
      MPI_COMM_WORLD, 0, 4, 1, distribution, sparse, time, steps, regrids, &
      minimum_dt, ok, bundle_sha256=bundle_sha256, &
      chemistry_integrator="explicit", base_mole_fractions=composition)
  end if
  local_rejected = .not. ok
  call MPI_Allreduce( &
    local_rejected, all_rejected, 1, MPI_LOGICAL, MPI_LAND, &
    MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS .or. .not. all_rejected) &
    error stop "Rank-dependent selected fingerprint presence was accepted"

  distribution = mpi_amr_eb_patch_tree_distribution_2d()
  distribution%comm = MPI_COMM_WORLD
  distribution%rank = rank
  distribution%nranks = nranks
  sparse = mpi_sparse_reactive_amr_eb_patch_tree_2d()
  sparse%nvar = reactive_nvar(size(species))
  allocate(sparse%levels(1))
  allocate(initial_integrals(sparse%nvar), source=0.0_dp)
  chemistry_advances = 0
  transport_advances = 0
  hydro_advances = 0
  if (rank == 0) then
    call write_sparse_owned_reactive_amr_eb_patch_tree_2d_checkpoint( &
      "selected_checkpoint_presence_should_not_be_written.dat", species, &
      distribution, sparse, 0, 0.0_dp, 0, 0, 1.0e-6_dp, ok, &
      local_entity_transfers=transfers, fingerprint=fingerprint, &
      minimum_transport_theta=1.0_dp, &
      initial_integrals=initial_integrals, &
      chemistry_level_advances=chemistry_advances, &
      transport_level_advances=transport_advances, &
      hydro_level_advances=hydro_advances, regrid_evaluations=0, &
      cumulative_tagged_cells=0, bundle_sha256=bundle_sha256, &
      chemistry_integrator="explicit", base_mole_fractions=composition)
  else
    call write_sparse_owned_reactive_amr_eb_patch_tree_2d_checkpoint( &
      "selected_checkpoint_presence_should_not_be_written.dat", species, &
      distribution, sparse, 0, 0.0_dp, 0, 0, 1.0e-6_dp, ok, &
      local_entity_transfers=transfers, &
      minimum_transport_theta=1.0_dp, &
      initial_integrals=initial_integrals, &
      chemistry_level_advances=chemistry_advances, &
      transport_level_advances=transport_advances, &
      hydro_level_advances=hydro_advances, regrid_evaluations=0, &
      cumulative_tagged_cells=0, bundle_sha256=bundle_sha256, &
      chemistry_integrator="explicit", base_mole_fractions=composition)
  end if
  local_rejected = .not. ok .and. transfers == 0
  call MPI_Allreduce( &
    local_rejected, all_rejected, 1, MPI_LOGICAL, MPI_LAND, &
    MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS .or. .not. all_rejected) &
    error stop "Rank-dependent selected write metadata was accepted"
  inquire( &
    file="selected_checkpoint_presence_should_not_be_written.dat", &
    exist=local_file_exists)
  call MPI_Allreduce( &
    local_file_exists, any_file_exists, 1, MPI_LOGICAL, MPI_LOR, &
    MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS .or. any_file_exists) &
    error stop "Rejected selected write published a checkpoint"

  if (rank == 0) &
    print '(a)', "test_mpi_amr_eb_selected_checkpoint_presence: PASS"
  call MPI_Finalize(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Finalize failed"

end program test_mpi_amr_eb_selected_checkpoint_presence
