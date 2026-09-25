program pelef_mpi_amr_spray_3d
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
  use spray_mpi_amr_application_mod, only: run_mpi_amr_spray_application
  implicit none

  type(reactive_3d_config) :: config
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  real(dp), allocatable :: mole_fractions(:)
  character(len=1024) :: input_path, topology_path, schedule_path, message
  logical :: ok
  integer :: ierr

  call MPI_Init(ierr)

  if (command_argument_count()/=3) then
    write(*, '(a)') "Usage: pelef_mpi_amr_spray_3d <input.nml> <topology.dat> <spray.nml>"
    call fatal('Invalid MPI tree command or input')
  end if
  call get_command_argument(1, input_path)
  call get_command_argument(2, topology_path)
  call read_reactive_3d_configuration( &
    trim(input_path), config, ok, message,&
    allow_pressure_outlet=.true.,allow_mach_inlet=.true.,allow_characteristic_ppm=.true.,allow_difmag=.true.)
  if (.not. ok) then
    write(*, '(a)') trim(message)
    call fatal('Invalid MPI tree command or input')
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
    call fatal('Failed to load fixed 3D AMR chemistry model')
  end if

  allocate(mole_fractions(size(species)))
  call reactive_3d_mole_fractions( &
    config, size(species), mole_fractions, ok)
  if (.not. ok) call fatal('Failed to resolve fixed 3D AMR composition')
  call get_command_argument(3,schedule_path)
  call run_mpi_amr_spray_application(trim(input_path),trim(topology_path),trim(schedule_path), &
    config,species,reactions,transport,mole_fractions,MPI_COMM_WORLD)
  call MPI_Finalize(ierr)
contains
  subroutine fatal(reason)
    character(len=*), intent(in) :: reason
    integer :: status
    write(*,'(a)') reason
    call MPI_Abort(MPI_COMM_WORLD,2,status)
    error stop 2
  end subroutine
end program pelef_mpi_amr_spray_3d
