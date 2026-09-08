program test_selected_reactive_1d_config
  use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use simulation_config_reactive_1d_mod, only: &
    reactive_1d_config, read_reactive_1d_configuration, &
    resolve_reactive_1d_selected_composition
  implicit none

  type(reactive_1d_config) :: config, valid_config
  type(nasa7_species) :: species(2)
  real(dp) :: mole_fractions(2)
  character(len=1024) :: input_path, amr_input_path, message
  logical :: ok

  if (command_argument_count() /= 2) then
    error stop "test_selected_reactive_1d_config requires regular and AMR inputs"
  end if
  call get_command_argument(1, input_path)
  call get_command_argument(2, amr_input_path)
  call read_reactive_1d_configuration( &
    trim(input_path), config, ok, message, allow_selected=.true.)
  if (.not. ok) error stop trim(message)
  if (trim(config%chemistry_model) /= "selected") then
    error stop "selected reactive 1D chemistry model was not read"
  end if
  valid_config = config

  call read_reactive_1d_configuration( &
    trim(amr_input_path), config, ok, message, allow_selected=.true.)
  if (.not. ok) error stop trim(message)
  if (.not. config%amr_enabled .or. config%amr_multipatch_enabled .or. &
      trim(config%amr_reconstruction) /= "plm" .or. &
      config%amr_refinement_ratio /= 2 .or. config%amr_max_levels /= 2 .or. &
      config%amr_regrid_interval /= 1 .or. &
      config%amr_minimum_patch_cells /= 4) then
    error stop "selected reactive 1D AMR controls were not read"
  end if
  call read_reactive_1d_configuration( &
    trim(amr_input_path), config, ok, message)
  if (ok .or. index(message, "Unknown reactive 1D chemistry model") == 0) then
    error stop "ordinary reactive 1D parser accepted selected AMR chemistry"
  end if

  call read_reactive_1d_configuration( &
    trim(input_path), config, ok, message)
  if (ok .or. index(message, "Unknown reactive 1D chemistry model") == 0) then
    error stop "ordinary reactive 1D parser accepted selected chemistry"
  end if

  species(1)%name = "H2"
  species(2)%name = "H"
  config = valid_config
  call resolve_reactive_1d_selected_composition( &
    config, species, mole_fractions, ok, message)
  if (.not. ok) error stop trim(message)
  if (abs(mole_fractions(1) - 0.8_dp) > 1.0e-14_dp .or. &
      abs(mole_fractions(2) - 0.2_dp) > 1.0e-14_dp) then
    error stop "selected reactive 1D composition was not mapped by name"
  end if

  config = valid_config
  config%composition_species(2) = "H"
  call resolve_reactive_1d_selected_composition( &
    config, species, mole_fractions, ok, message)
  if (ok .or. index(message, "duplicate species") == 0) then
    error stop "duplicate selected reactive 1D species was not rejected"
  end if

  config = valid_config
  config%composition_species(2) = "UNKNOWN"
  call resolve_reactive_1d_selected_composition( &
    config, species, mole_fractions, ok, message)
  if (ok .or. index(message, "not found exactly once") == 0) then
    error stop "unknown selected reactive 1D species was not rejected"
  end if

  config = valid_config
  config%composition_mole_fractions(3) = &
    ieee_value(0.0_dp, ieee_quiet_nan)
  call resolve_reactive_1d_selected_composition( &
    config, species, mole_fractions, ok, message)
  if (ok .or. index(message, "nonfinite mole fraction") == 0) then
    error stop "nonfinite selected reactive 1D tail was not rejected"
  end if

  config = valid_config
  config%problem = "composition_wave"
  config%composition_wave_amplitude = 0.21_dp
  call resolve_reactive_1d_selected_composition( &
    config, species, mole_fractions, ok, message)
  if (ok .or. index(message, "endpoint fractions") == 0) then
    error stop "invalid selected reactive 1D composition wave was not rejected"
  end if

  config = valid_config
  config%composition_mole_fractions(1:2) = huge(1.0_dp)
  call resolve_reactive_1d_selected_composition( &
    config, species, mole_fractions, ok, message)
  if (.not. ok .or. any(abs(mole_fractions - 0.5_dp) > 1.0e-14_dp)) then
    error stop "large selected reactive 1D composition was not scaled safely"
  end if

  write(*, '(a)') "test_selected_reactive_1d_config: PASS"
end program test_selected_reactive_1d_config
