program pelef_eb_spray_3d
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_elementary_transport
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use simulation_config_reactive_eb_3d_mod, only: &
    reactive_eb_3d_config, read_reactive_eb_3d_configuration
  use spray_eb_application_mod, only: &
    run_eb_spray_application
  implicit none

  type(reactive_eb_3d_config) :: config
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  character(len=1024) :: input_path, spray_path, message
  logical :: ok

  if (command_argument_count() /= 2) then
    write(*, '(a)') "Usage: pelef_eb_spray_3d <gas.nml> <spray.nml>"
    error stop 2
  end if
  call get_command_argument(1, input_path)
  call read_reactive_eb_3d_configuration( &
    trim(input_path), config, ok, message, allow_difmag=.true.,allow_ctu=.true.)
  if (.not. ok) then
    write(*, '(a)') trim(message)
    error stop 2
  end if
  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "Failed to load elementary thermodynamics"
  call load_h2o2_elementary_mechanism(reactions, ok)
  if (.not. ok) error stop "Failed to load elementary mechanism"
  call load_h2o2_elementary_transport(transport, ok)
  if (.not. ok) error stop "Failed to load elementary transport"
  call get_command_argument(2,spray_path)
  call run_eb_spray_application(trim(input_path),trim(spray_path),config,species,reactions,transport)
end program pelef_eb_spray_3d
