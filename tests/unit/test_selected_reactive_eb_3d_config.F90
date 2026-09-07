program test_selected_reactive_eb_3d_config
  use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use simulation_config_reactive_eb_3d_mod, only: &
    reactive_eb_3d_config, read_reactive_eb_3d_configuration, &
    resolve_reactive_eb_3d_selected_composition
  implicit none

  type(reactive_eb_3d_config) :: config, valid_config
  type(nasa7_species) :: species(2)
  real(dp) :: mole_fractions(2)
  character(len=1024) :: input_path, message
  logical :: ok

  if (command_argument_count() /= 1) then
    error stop "test_selected_reactive_eb_3d_config requires an input path"
  end if
  call get_command_argument(1, input_path)
  call read_reactive_eb_3d_configuration( &
    trim(input_path), config, ok, message, allow_selected=.true.)
  if (.not. ok) error stop trim(message)
  if (trim(config%thermo_model) /= "selected" .or. &
      config%composition_count /= 2) then
    error stop "selected reactive EB 3D configuration was not read"
  end if
  valid_config = config

  call read_reactive_eb_3d_configuration( &
    trim(input_path), config, ok, message)
  if (ok .or. index(message, "explicit opt-in") == 0) then
    error stop "ordinary reactive EB 3D parser accepted selected thermo"
  end if

  species(1)%name = "H2"
  species(2)%name = "H"
  config = valid_config
  call resolve_reactive_eb_3d_selected_composition( &
    config, species, mole_fractions, ok, message)
  if (.not. ok) error stop trim(message)
  if (abs(mole_fractions(1) - 0.8_dp) > 1.0e-14_dp .or. &
      abs(mole_fractions(2) - 0.2_dp) > 1.0e-14_dp) then
    error stop "selected reactive EB 3D composition was not mapped by name"
  end if

  config = valid_config
  config%composition_species(2) = "H"
  call resolve_reactive_eb_3d_selected_composition( &
    config, species, mole_fractions, ok, message)
  if (ok .or. index(message, "duplicate species") == 0) then
    error stop "duplicate selected reactive EB 3D species was not rejected"
  end if

  config = valid_config
  config%composition_species(2) = "UNKNOWN"
  call resolve_reactive_eb_3d_selected_composition( &
    config, species, mole_fractions, ok, message)
  if (ok .or. index(message, "not found exactly once") == 0) then
    error stop "unknown selected reactive EB 3D species was not rejected"
  end if

  config = valid_config
  config%composition_mole_fractions(3) = &
    ieee_value(0.0_dp, ieee_quiet_nan)
  call resolve_reactive_eb_3d_selected_composition( &
    config, species, mole_fractions, ok, message)
  if (ok .or. index(message, "nonfinite mole fraction") == 0) then
    error stop "nonfinite selected reactive EB 3D tail was not rejected"
  end if

  config = valid_config
  config%composition_mole_fractions(1:2) = huge(1.0_dp)
  call resolve_reactive_eb_3d_selected_composition( &
    config, species, mole_fractions, ok, message)
  if (.not. ok .or. any(abs(mole_fractions - 0.5_dp) > 1.0e-14_dp)) then
    error stop "large selected reactive EB 3D composition was not scaled"
  end if

  write(*, '(a)') "test_selected_reactive_eb_3d_config: PASS"
end program test_selected_reactive_eb_3d_config
