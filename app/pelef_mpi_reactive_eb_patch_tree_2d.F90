program pelef_mpi_reactive_eb_patch_tree_2d
  use, intrinsic :: iso_fortran_env, only: error_unit
  use mpi_f08
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: &
    load_h2o2_elementary_transport, load_h2o2_full_transport
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use simulation_config_reactive_eb_amr_2d_mod, only: &
    reactive_eb_amr_2d_config, read_reactive_eb_amr_2d_configuration
  use mpi_reactive_eb_patch_tree_2d_application_mod, only: &
    run_mpi_reactive_eb_patch_tree_2d_application
  implicit none

  type(reactive_eb_amr_2d_config) :: config
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  character(len=1024) :: input_path, message, output_path
  integer :: argument_count, ierr, rank
  logical :: ok

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Init failed"
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Comm_rank failed"

  argument_count = command_argument_count()
  if (argument_count < 1 .or. argument_count > 2) call abort_run( &
    "Usage: pelef_mpi_reactive_eb_patch_tree_2d <input.nml> [output.csv]", 2)
  call get_command_argument(1, input_path)
  call read_reactive_eb_amr_2d_configuration( &
    trim(input_path), config, ok, message)
  if (.not. ok) call abort_run(trim(message), 2)
  output_path = config%eb%flow%output_file
  if (argument_count == 2) call get_command_argument(2, output_path)
  if (len_trim(output_path) == 0) call abort_run("Output path is empty", 2)

  select case (trim(config%eb%flow%chemistry_model))
  case ("elementary")
    call load_h2o2_elementary_thermo(species, ok)
    if (.not. ok) call abort_run("Failed to load elementary thermodynamics", 3)
    call load_h2o2_elementary_mechanism(reactions, ok)
    if (.not. ok) call abort_run("Failed to load elementary mechanism", 3)
    call load_h2o2_elementary_transport(transport, ok)
    if (.not. ok) call abort_run("Failed to load elementary transport", 3)
  case ("full_h2o2")
    call load_h2o2_full_thermo(species, ok)
    if (.not. ok) call abort_run("Failed to load full H2/O2 thermodynamics", 3)
    call load_h2o2_full_mechanism(reactions, ok)
    if (.not. ok) call abort_run("Failed to load full H2/O2 mechanism", 3)
    call load_h2o2_full_transport(transport, ok)
    if (.not. ok) call abort_run("Failed to load full H2/O2 transport", 3)
  case default
    call abort_run("Unknown chemistry model", 3)
  end select

  call run_mpi_reactive_eb_patch_tree_2d_application( &
    MPI_COMM_WORLD, "sparse MPI reactive EB patch-tree 2D", config, species, &
    reactions, transport, trim(output_path))

  call MPI_Finalize(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Finalize failed"

contains

  subroutine abort_run(reason, code)
    character(len=*), intent(in) :: reason
    integer, intent(in) :: code

    integer :: abort_ierr

    if (rank == 0) write(error_unit, '(a)') trim(reason)
    call MPI_Abort(MPI_COMM_WORLD, code, abort_ierr)
    error stop code
  end subroutine abort_run

end program pelef_mpi_reactive_eb_patch_tree_2d
