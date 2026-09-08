program test_selected_reactive_2d_config
  use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use simulation_config_reactive_2d_mod, only: &
    reactive_2d_config, read_reactive_2d_configuration, &
    resolve_reactive_2d_selected_composition
  implicit none

  type(reactive_2d_config) :: config, valid_config
  type(nasa7_species) :: species(2)
  real(dp) :: mole_fractions(2)
  character(len=1024) :: input_path, message
  logical :: ok

  if (command_argument_count() /= 1) then
    error stop "test_selected_reactive_2d_config requires an input path"
  end if
  call get_command_argument(1, input_path)
  call read_reactive_2d_configuration( &
    trim(input_path), config, ok, message, allow_selected=.true.)
  if (.not. ok) error stop trim(message)
  if (trim(config%chemistry_model) /= "selected") then
    error stop "selected reactive 2D chemistry model was not read"
  end if
  valid_config = config

  call read_reactive_2d_configuration( &
    trim(input_path), config, ok, message)
  if (ok .or. index(message, "Unknown reactive 2D chemistry model") == 0) then
    error stop "ordinary reactive 2D parser accepted selected chemistry"
  end if

  species(1)%name = "H2"
  species(2)%name = "H"
  config = valid_config
  call resolve_reactive_2d_selected_composition( &
    config, species, mole_fractions, ok, message)
  if (.not. ok) error stop trim(message)
  if (abs(mole_fractions(1) - 0.8_dp) > 1.0e-14_dp .or. &
      abs(mole_fractions(2) - 0.2_dp) > 1.0e-14_dp) then
    error stop "selected reactive 2D composition was not mapped by name"
  end if

  config = valid_config
  config%composition_species(2) = "H"
  call resolve_reactive_2d_selected_composition( &
    config, species, mole_fractions, ok, message)
  if (ok .or. index(message, "duplicate species") == 0) then
    error stop "duplicate selected reactive 2D species was not rejected"
  end if

  config = valid_config
  config%composition_species(2) = "UNKNOWN"
  call resolve_reactive_2d_selected_composition( &
    config, species, mole_fractions, ok, message)
  if (ok .or. index(message, "not found exactly once") == 0) then
    error stop "unknown selected reactive 2D species was not rejected"
  end if

  config = valid_config
  config%composition_mole_fractions(3) = &
    ieee_value(0.0_dp, ieee_quiet_nan)
  call resolve_reactive_2d_selected_composition( &
    config, species, mole_fractions, ok, message)
  if (ok .or. index(message, "nonfinite mole fraction") == 0) then
    error stop "nonfinite selected reactive 2D tail was not rejected"
  end if

  config = valid_config
  config%problem = "diagonal_composition_wave"
  config%composition_wave_amplitude = 0.21_dp
  call resolve_reactive_2d_selected_composition( &
    config, species, mole_fractions, ok, message)
  if (ok .or. index(message, "endpoint fractions") == 0) then
    error stop "invalid selected reactive 2D composition wave was not rejected"
  end if

  config = valid_config
  config%boundary_x_lower = "slip_wall"
  config%wall_species_x_lower = "prescribed"
  config%prescribed_species_flux_x_lower(1:2) = [0.05_dp, -0.05_dp]
  call resolve_reactive_2d_selected_composition( &
    config, species, mole_fractions, ok, message)
  if (.not. ok) error stop "valid selected reactive 2D wall flux was rejected"
  config%prescribed_species_flux_x_lower(3) = 0.01_dp
  call resolve_reactive_2d_selected_composition( &
    config, species, mole_fractions, ok, message)
  if (ok .or. index(message, "wall species flux") == 0) then
    error stop "selected reactive 2D wall-flux tail was not rejected"
  end if

  config = valid_config
  config%composition_mole_fractions(1:2) = huge(1.0_dp)
  call resolve_reactive_2d_selected_composition( &
    config, species, mole_fractions, ok, message)
  if (.not. ok .or. any(abs(mole_fractions - 0.5_dp) > 1.0e-14_dp)) then
    error stop "large selected reactive 2D composition was not scaled safely"
  end if

  write(*, '(a)') "test_selected_reactive_2d_config: PASS"
end program test_selected_reactive_2d_config
