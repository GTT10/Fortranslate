program pelef_mpi_amr_reactive_3d
  use, intrinsic :: iso_fortran_env, only: error_unit
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use transport_database_mod, only: &
    load_h2o2_elementary_transport, load_h2o2_full_transport
  use simulation_config_reactive_3d_mod, only: &
    reactive_3d_config, read_reactive_3d_configuration, &
    reactive_3d_mole_fractions
  use simulation_config_amr_reactive_3d_mod, only: &
    amr_reactive_3d_config, read_amr_reactive_3d_configuration
  use mpi_amr_reactive_3d_application_mod, only: &
    run_mpi_amr_reactive_3d_application
  implicit none

  type(reactive_3d_config) :: config
  type(amr_reactive_3d_config) :: amr_config
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  real(dp), allocatable :: mole_fractions(:)
  character(len=1024) :: input_path, message
  character(len=512) :: output_prefix, coarse_output, fine_output
  logical :: ok
  integer :: argument_count, ierr, rank

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Init failed"
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Comm_rank failed"

  argument_count = command_argument_count()
  if (argument_count < 1 .or. argument_count > 2) then
    call abort_run( &
      "Usage: pelef_mpi_amr_reactive_3d <input.nml> [output-prefix]", 2)
  end if
  call get_command_argument(1, input_path)
  call read_reactive_3d_configuration(trim(input_path), config, ok, message)
  if (.not. ok) call abort_run(trim(message), 2)
  call read_amr_reactive_3d_configuration( &
    trim(input_path), config, amr_config, ok, message)
  if (.not. ok) call abort_run(trim(message), 2)
  if (argument_count == 2) then
    call get_command_argument(2, output_prefix)
    if (len_trim(output_prefix) == 0) then
      call abort_run("MPI output prefix is empty", 2)
    end if
    coarse_output = trim(output_prefix) // "_coarse.csv"
    fine_output = trim(output_prefix) // "_fine.csv"
  else
    coarse_output = amr_config%coarse_output_file
    fine_output = amr_config%fine_output_file
  end if

  select case (trim(config%thermo_model))
  case ("elementary")
    call load_h2o2_elementary_thermo(species, ok)
    if (ok) call load_h2o2_elementary_mechanism(reactions, ok)
    if (ok) call load_h2o2_elementary_transport(transport, ok)
  case ("full_h2o2")
    call load_h2o2_full_thermo(species, ok)
    if (ok) call load_h2o2_full_mechanism(reactions, ok)
    if (ok) call load_h2o2_full_transport(transport, ok)
  case default
    ok = .false.
  end select
  if (.not. ok .or. .not. allocated(species) .or. &
      .not. allocated(reactions) .or. .not. allocated(transport)) then
    call abort_run("Failed to load MPI 3D AMR chemistry model", 2)
  end if

  allocate(mole_fractions(size(species)))
  call reactive_3d_mole_fractions( &
    config, size(species), mole_fractions, ok)
  if (.not. ok) call abort_run("Failed to resolve MPI 3D AMR composition", 2)
  call run_mpi_amr_reactive_3d_application( &
    MPI_COMM_WORLD, trim(input_path), &
    "rank-local sparse distributed-slab static two-level 3D reactive AMR", &
    config, amr_config, species, reactions, transport, mole_fractions, &
    trim(coarse_output), trim(fine_output))

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

end program pelef_mpi_amr_reactive_3d
