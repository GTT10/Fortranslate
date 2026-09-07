program test_mpi_amr_eb_selected_context_mismatch
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_elementary_transport
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use simulation_config_reactive_eb_amr_2d_mod, only: &
    reactive_eb_amr_2d_config
  use mpi_reactive_eb_patch_tree_2d_application_mod, only: &
    run_mpi_reactive_eb_patch_tree_2d_application
  implicit none

  character(len=*), parameter :: bundle_sha256 = &
    "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef"

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(reactive_eb_amr_2d_config) :: config
  real(dp), allocatable :: composition(:)
  logical :: ok
  integer :: ierr, rank

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Init failed"
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Comm_rank failed"
  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "MPI EB context thermodynamics load failed"
  call load_h2o2_elementary_mechanism(reactions, ok)
  if (.not. ok) error stop "MPI EB context mechanism load failed"
  call load_h2o2_elementary_transport(transport, ok)
  if (.not. ok) error stop "MPI EB context transport load failed"

  config = reactive_eb_amr_2d_config()
  config%eb%flow%chemistry_model = "selected"
  config%three_level_enabled = .false.
  config%multipatch_enabled = .false.
  allocate(composition(size(species)))
  composition = [ &
    0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, 1.0e-5_dp, &
    0.0_dp, 0.55643_dp]
  if (rank == 1) then
    composition(1) = composition(1) + 1.0e-6_dp
    composition(2) = composition(2) - 1.0e-6_dp
  end if

  call run_mpi_reactive_eb_patch_tree_2d_application( &
    MPI_COMM_WORLD, "selected MPI EB context mismatch test", config, &
    species, reactions, transport, "mpi_amr_eb_context_mismatch.csv", &
    bundle_sha256, composition, "explicit")
  error stop "Rank-dependent selected MPI EB context was unexpectedly accepted"

end program test_mpi_amr_eb_selected_context_mismatch
