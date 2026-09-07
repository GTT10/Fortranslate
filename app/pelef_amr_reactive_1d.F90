program pelef_amr_reactive_1d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: &
    load_h2o2_elementary_transport, &
    load_h2o2_full_transport
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use simulation_config_reactive_1d_mod, only: &
    reactive_1d_config, read_reactive_1d_configuration, &
    reactive_1d_mole_fractions
  use amr_reactive_1d_application_mod, only: &
    run_amr_reactive_1d_application
  implicit none

  type(reactive_1d_config) :: config
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  real(dp), allocatable :: base_mole_fractions(:)
  character(len=1024) :: input_path, message
  logical :: ok

  if (command_argument_count() /= 1) then
    write(*, '(a)') "Usage: pelef_amr_reactive_1d <input.nml>"
    error stop 2
  end if
  call get_command_argument(1, input_path)
  call read_reactive_1d_configuration(trim(input_path), config, ok, message)
  if (.not. ok) then
    write(*, '(a)') trim(message)
    error stop 2
  end if
  if (.not. config%amr_enabled) &
    error stop "AMR reactive application requires amr_enabled"

  select case (trim(config%chemistry_model))
  case ("elementary")
    call load_h2o2_elementary_thermo(species, ok)
    if (.not. ok) error stop "Failed to load elementary thermodynamics"
    call load_h2o2_elementary_mechanism(reactions, ok)
    if (.not. ok) error stop "Failed to load elementary mechanism"
    call load_h2o2_elementary_transport(transport, ok)
    if (.not. ok) error stop "Failed to load elementary transport"
  case ("full_h2o2")
    call load_h2o2_full_thermo(species, ok)
    if (.not. ok) error stop "Failed to load full H2/O2 thermodynamics"
    call load_h2o2_full_mechanism(reactions, ok)
    if (.not. ok) error stop "Failed to load full H2/O2 mechanism"
    call load_h2o2_full_transport(transport, ok)
    if (.not. ok) error stop "Failed to load full H2/O2 transport"
  case default
    error stop "Unknown chemistry model"
  end select

  allocate(base_mole_fractions(size(species)))
  call reactive_1d_mole_fractions( &
    config, size(species), base_mole_fractions, ok)
  if (.not. ok) error stop "Failed to resolve AMR reactive composition"
  call run_amr_reactive_1d_application( &
    trim(input_path), "AMR reactive 1D", "", config, species, reactions, &
    transport, base_mole_fractions)
end program pelef_amr_reactive_1d
