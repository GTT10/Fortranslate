program test_selected_reactor_config
  use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use simulation_config_selected_reactor_mod, only: &
    selected_reactor_config, selected_reactor_max_species, &
    read_selected_reactor_configuration, &
    validate_selected_reactor_configuration, &
    resolve_selected_reactor_composition, &
    selected_reactor_paths_alias_lexically
  implicit none

  type(selected_reactor_config) :: config, valid_config
  type(nasa7_species) :: species(2)
  real(dp) :: mole_fractions(2)
  character(len=1024) :: input_path, message
  logical :: ok
  integer :: species_index

  if (command_argument_count() /= 1) then
    error stop "test_selected_reactor_config requires an input path"
  end if
  call get_command_argument(1, input_path)
  call read_selected_reactor_configuration(trim(input_path), config, ok, message)
  if (.not. ok) error stop trim(message)
  if (config%composition_count /= 2) then
    error stop "selected reactor composition count was not read"
  end if
  if (trim(config%integrator) /= "implicit") then
    error stop "selected reactor default integrator was not implicit"
  end if
  valid_config = config

  config = valid_config
  config%integrator = "invalid"
  call validate_selected_reactor_configuration(config, ok, message)
  if (ok .or. index(message, "integrator must be") == 0) then
    error stop "unknown selected reactor integrator was not rejected"
  end if
  config = valid_config
  config%integrator = "cvode"
  call validate_selected_reactor_configuration(config, ok, message)
  if (.not. ok) then
    error stop "valid CVODE selected reactor configuration was rejected"
  end if

  species(1)%name = "H"
  species(2)%name = "H2"
  call resolve_selected_reactor_composition( &
    config, species, mole_fractions, ok, message)
  if (.not. ok) error stop trim(message)
  if (abs(mole_fractions(1) - 0.2_dp) > 1.0e-14_dp .or. &
      abs(mole_fractions(2) - 0.8_dp) > 1.0e-14_dp) then
    error stop "selected reactor composition was not mapped by species name"
  end if

  config = valid_config
  config%composition_species(2) = "H"
  call resolve_selected_reactor_composition( &
    config, species, mole_fractions, ok, message)
  if (ok .or. index(message, "duplicate species") == 0) then
    error stop "duplicate selected reactor species was not rejected"
  end if

  config = valid_config
  config%composition_species(2) = "UNKNOWN"
  call resolve_selected_reactor_composition( &
    config, species, mole_fractions, ok, message)
  if (ok .or. index(message, "not found exactly once") == 0) then
    error stop "unknown selected reactor species was not rejected"
  end if

  config = valid_config
  config%final_time = ieee_value(0.0_dp, ieee_quiet_nan)
  call validate_selected_reactor_configuration(config, ok, message)
  if (ok .or. index(message, "must be finite") == 0) then
    error stop "nonfinite selected reactor time was not rejected"
  end if

  config = valid_config
  config%initial_time_step = 0.5_dp * config%minimum_time_step
  call validate_selected_reactor_configuration(config, ok, message)
  if (ok .or. index(message, "within the time-step bounds") == 0) then
    error stop "selected reactor initial step below the minimum was not rejected"
  end if

  config = valid_config
  config%initial_time_step = 2.0_dp * config%maximum_time_step
  call validate_selected_reactor_configuration(config, ok, message)
  if (ok .or. index(message, "within the time-step bounds") == 0) then
    error stop "selected reactor initial step above the maximum was not rejected"
  end if

  config = valid_config
  config%cvode_max_internal_steps = 0
  call validate_selected_reactor_configuration(config, ok, message)
  if (ok .or. index(message, "cvode_max_internal_steps") == 0) then
    error stop "nonpositive CVODE internal-step limit was not rejected"
  end if

  config = valid_config
  config%output_interval = 0.5_dp * config%minimum_time_step
  call validate_selected_reactor_configuration(config, ok, message)
  if (ok .or. index(message, "output_interval is below") == 0) then
    error stop "selected reactor output interval below the minimum was not rejected"
  end if

  config = valid_config
  config%composition_mole_fractions(3) = &
    ieee_value(0.0_dp, ieee_quiet_nan)
  call validate_selected_reactor_configuration(config, ok, message)
  if (ok .or. index(message, "nonfinite mole fraction") == 0) then
    error stop "nonfinite selected reactor composition tail was not rejected"
  end if

  config = valid_config
  config%composition_species(3) = "EXTRA"
  call validate_selected_reactor_configuration(config, ok, message)
  if (ok .or. index(message, "names beyond") == 0) then
    error stop "selected reactor composition tail was not rejected"
  end if

  config = selected_reactor_config()
  config%composition_count = selected_reactor_max_species
  do species_index = 1, selected_reactor_max_species
    write(config%composition_species(species_index), '("S",i0)') species_index
  end do
  config%composition_mole_fractions(1) = 1.0_dp
  call validate_selected_reactor_configuration(config, ok, message)
  if (.not. ok) then
    error stop "maximum-size selected reactor configuration was rejected"
  end if

  if (.not. selected_reactor_paths_alias_lexically( &
      "reactor.nml", "./reactor.nml")) then
    error stop "selected reactor dot-component alias was not detected"
  end if
  if (.not. selected_reactor_paths_alias_lexically( &
      "reactor.nml", "scratch/../reactor.nml")) then
    error stop "selected reactor parent-component alias was not detected"
  end if
  if (selected_reactor_paths_alias_lexically( &
      "input/reactor.nml", "output/reactor.nml")) then
    error stop "distinct selected reactor paths were rejected"
  end if

  config = valid_config
  config%composition_mole_fractions(1:2) = huge(1.0_dp)
  call resolve_selected_reactor_composition( &
    config, species, mole_fractions, ok, message)
  if (.not. ok .or. any(abs(mole_fractions - 0.5_dp) > 1.0e-14_dp)) then
    error stop "large selected reactor composition was not scaled safely"
  end if

  write(*, '(a)') "test_selected_reactor_config: PASS"
end program test_selected_reactor_config
