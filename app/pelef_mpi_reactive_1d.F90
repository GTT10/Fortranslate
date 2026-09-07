program pelef_mpi_reactive_1d
  use mpi_f08
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use transport_database_mod, only: load_h2o2_full_transport
  use mpi_reactive_1d_application_mod, only: &
    run_mpi_reactive_1d_application
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  character(len=256) :: output_file
  integer :: ierr
  logical :: ok

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Init failed"

  if (command_argument_count() >= 1) then
    call get_command_argument(1, output_file)
  else
    output_file = "mpi_reactive_1d.csv"
  end if

  call load_h2o2_full_thermo(species, ok)
  if (.not. ok) error stop "Failed to load full H2/O2 thermodynamics"
  call load_h2o2_full_mechanism(reactions, ok)
  if (.not. ok) error stop "Failed to load full H2/O2 mechanism"
  call load_h2o2_full_transport(transport, ok)
  if (.not. ok) error stop "Failed to load full H2/O2 transport data"

  call run_mpi_reactive_1d_application( &
    MPI_COMM_WORLD, species, reactions, transport, trim(output_file))

  call MPI_Finalize(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Finalize failed"
end program pelef_mpi_reactive_1d
