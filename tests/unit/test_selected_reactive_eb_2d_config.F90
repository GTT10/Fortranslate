program test_selected_reactive_eb_2d_config
  use simulation_config_reactive_eb_2d_mod, only: &
    reactive_eb_2d_config, read_reactive_eb_2d_configuration
  implicit none

  type(reactive_eb_2d_config) :: config
  character(len=1024) :: path, message
  logical :: ok

  if (command_argument_count() /= 1) then
    error stop "Expected selected reactive EB 2D fixture path"
  end if
  call get_command_argument(1, path)
  call read_reactive_eb_2d_configuration( &
    trim(path), config, ok, message, allow_selected=.true.)
  if (.not. ok) then
    write(*, '(a)') trim(message)
    error stop "Selected reactive EB 2D configuration was rejected"
  end if
  if (trim(config%flow%chemistry_model) /= "selected" .or. &
      config%flow%composition_count /= 2 .or. &
      trim(config%geometry) /= "plane" .or. &
      config%state_redist_max_order /= 2 .or. &
      trim(config%flow%boundary_x_lower) /= "outflow" .or. &
      trim(config%flow%boundary_y_upper) /= "outflow") then
    error stop "Selected reactive EB 2D configuration fields are incorrect"
  end if

  call read_reactive_eb_2d_configuration( &
    trim(path), config, ok, message)
  if (ok .or. index(message, "Unknown reactive 2D chemistry model") == 0) then
    error stop "Fixed reactive EB 2D parser accepted selected chemistry"
  end if

  write(*, '(a)') "selected reactive EB 2D configuration: PASS"
end program test_selected_reactive_eb_2d_config
