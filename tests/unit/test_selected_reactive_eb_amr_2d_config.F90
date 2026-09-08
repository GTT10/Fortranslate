program test_selected_reactive_eb_amr_2d_config
  use simulation_config_reactive_eb_amr_2d_mod, only: &
    reactive_eb_amr_2d_config, read_reactive_eb_amr_2d_configuration
  implicit none

  type(reactive_eb_amr_2d_config) :: config
  character(len=1024) :: path, message
  logical :: ok

  if (command_argument_count() /= 1) then
    error stop "Expected selected reactive EB AMR 2D fixture path"
  end if
  call get_command_argument(1, path)
  call read_reactive_eb_amr_2d_configuration( &
    trim(path), config, ok, message, allow_selected=.true.)
  if (.not. ok) then
    write(*, '(a)') trim(message)
    error stop "Selected reactive EB AMR 2D configuration was rejected"
  end if
  if (trim(config%eb%flow%chemistry_model) /= "selected" .or. &
      config%eb%flow%composition_count /= 2 .or. &
      trim(config%eb%geometry) /= "plane" .or. &
      config%coarse_i_lower /= 2 .or. config%coarse_i_upper /= 7 .or. &
      config%coarse_j_lower /= 2 .or. config%coarse_j_upper /= 7 .or. &
      config%refinement_ratio /= 2 .or. &
      trim(config%prolongation_method) /= "linear" .or. &
      config%dynamic_regridding .or. config%three_level_enabled .or. &
      config%multipatch_enabled .or. &
      trim(config%fine_output_file) /= &
        "selected_fixture_reactive_eb_amr_2d_fine.csv") then
    error stop &
      "Selected reactive EB AMR 2D configuration fields are incorrect"
  end if

  call read_reactive_eb_amr_2d_configuration( &
    trim(path), config, ok, message)
  if (ok .or. index(message, "Unknown reactive 2D chemistry model") == 0) then
    error stop "Fixed reactive EB AMR 2D parser accepted selected chemistry"
  end if

  write(*, '(a)') "selected reactive EB AMR 2D configuration: PASS"
end program test_selected_reactive_eb_amr_2d_config
