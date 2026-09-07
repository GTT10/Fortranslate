program pelef_amr_reactive_3d
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
  use amr_reactive_3d_application_mod, only: &
    run_amr_reactive_3d_application
  implicit none

  type(reactive_3d_config) :: config
  type(amr_reactive_3d_config) :: amr_config
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  real(dp), allocatable :: mole_fractions(:)
  character(len=1024) :: input_path, message
  logical :: ok

  if (command_argument_count() /= 1) then
    write(*, '(a)') "Usage: pelef_amr_reactive_3d <input.nml>"
    error stop 2
  end if
  call get_command_argument(1, input_path)
  call read_reactive_3d_configuration(trim(input_path), config, ok, message)
  if (.not. ok) then
    write(*, '(a)') trim(message)
    error stop 2
  end if
  call read_amr_reactive_3d_configuration( &
    trim(input_path), config, amr_config, ok, message)
  if (.not. ok) then
    write(*, '(a)') trim(message)
    error stop 2
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
    error stop "Failed to load fixed 3D AMR chemistry model"
  end if

  allocate(mole_fractions(size(species)))
  call reactive_3d_mole_fractions( &
    config, size(species), mole_fractions, ok)
  if (.not. ok) error stop "Failed to resolve fixed 3D AMR composition"
  call run_amr_reactive_3d_application( &
    trim(input_path), "static two-level 3D reactive AMR", "", &
    config, amr_config, species, reactions, transport, mole_fractions)
end program pelef_amr_reactive_3d
