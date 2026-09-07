program test_reactive_eb_3d_config
  use, intrinsic :: ieee_arithmetic, only: &
    ieee_value, ieee_quiet_nan
  use precision_mod, only: dp
  use simulation_config_reactive_eb_3d_mod, only: &
    reactive_eb_3d_config, read_reactive_eb_3d_configuration, &
    validate_reactive_eb_3d_configuration
  implicit none

  character(len=*), parameter :: path = "reactive_eb_3d_config_test.nml"
  type(reactive_eb_3d_config) :: config
  character(len=256) :: message
  logical :: ok
  integer :: unit

  config = reactive_eb_3d_config()
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(ok, "default public 3D EB configuration")
  config%chemistry_enabled = .true.
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(ok, "reacting public 3D EB configuration")

  config%plane_axis = "y"
  config%plane_position = 0.49375_dp
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(ok, "y-normal public 3D EB configuration")
  config%plane_axis = "z"
  config%plane_position = 0.4916666666666667_dp
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(ok, "z-normal public 3D EB configuration")

  config = reactive_eb_3d_config()
  config%state_redist_target_volume_fraction = 0.0_dp
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "invalid public StateRedist target")
  config = reactive_eb_3d_config()
  config%state_redist_target_volume_fraction = 0.05_dp
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "inactive public StateRedist target")
  config%redistribution = "flux_redist"
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(ok, "FluxRedist ignores StateRedist target activation")
  config = reactive_eb_3d_config()
  config%redistribution = "raw"
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "unstabilized public route rejection")
  config = reactive_eb_3d_config()
  config%plane_position = 0.4_dp
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "face-aligned public plane rejection")
  config = reactive_eb_3d_config()
  config%plane_position = 0.95_dp
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "receiver-free public plane rejection")
  config = reactive_eb_3d_config()
  config%initial_pressure = ieee_value(0.0_dp, ieee_quiet_nan)
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "nonfinite public configuration rejection")
  config = reactive_eb_3d_config()
  config%chemistry_relative_tolerance = 0.0_dp
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "nonpositive chemistry tolerance rejection")
  config = reactive_eb_3d_config()
  config%chemistry_absolute_tolerance = &
    ieee_value(0.0_dp, ieee_quiet_nan)
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "nonfinite chemistry tolerance rejection")
  config = reactive_eb_3d_config()
  config%transport_enabled = .true.
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(ok, "transport public 3D EB configuration")
  config%redistribution = "flux_redist"
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "transport requires StateRedist")
  config = reactive_eb_3d_config()
  config%transport_enabled = .true.
  config%viscosity_enabled = .false.
  config%thermal_conduction_enabled = .false.
  config%species_diffusion_enabled = .false.
  config%barodiffusion_enabled = .false.
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "transport process rejection")
  config%thermal_conduction_enabled = .true.
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(ok, "thermal-only transport configuration")
  config%barodiffusion_enabled = .true.
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "barodiffusion requires species diffusion")
  config = reactive_eb_3d_config()
  config%transport_cfl = 0.0_dp
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "nonpositive transport CFL rejection")
  config%transport_cfl = ieee_value(0.0_dp, ieee_quiet_nan)
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "nonfinite transport CFL rejection")
  config = reactive_eb_3d_config()
  config%checkpoint_interval_steps = -1
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "negative checkpoint interval rejection")
  config = reactive_eb_3d_config()
  config%stop_after_checkpoint = .true.
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "checkpoint stop without interval rejection")
  config%checkpoint_interval_steps = 1
  config%checkpoint_file = ""
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "empty checkpoint path rejection")
  config%checkpoint_file = config%output_file
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "checkpoint output collision rejection")
  config%checkpoint_file = "./reactive_eb_3d.csv"
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "lexical checkpoint output alias rejection")
  config%checkpoint_file = "scratch/../reactive_eb_3d.csv"
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "parent checkpoint output alias rejection")
  config = reactive_eb_3d_config()
  config%restart_file = config%output_file
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "restart output collision rejection")
  config%restart_file = "./reactive_eb_3d.csv"
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "lexical restart output alias rejection")
  config = reactive_eb_3d_config()
  config%checkpoint_interval_steps = 2
  config%checkpoint_file = "resume.chk"
  config%restart_file = "resume.chk"
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(ok, "same checkpoint and restart path")

  open(newunit=unit, file=path, status="replace", action="write")
  write(unit, '(a)') "&reactive_eb_3d"
  write(unit, '(a)') " nx = 12"
  write(unit, '(a)') " plane_position = 0.3291666666666667"
  write(unit, '(a)') " redistribution = 'state_redist'"
  write(unit, '(a)') " chemistry_enabled = .true."
  write(unit, '(a)') " chemistry_relative_tolerance = 1.0e-8"
  write(unit, '(a)') " chemistry_absolute_tolerance = 1.0e-13"
  write(unit, '(a)') " transport_enabled = .true."
  write(unit, '(a)') " viscosity_enabled = .false."
  write(unit, '(a)') " thermal_conduction_enabled = .true."
  write(unit, '(a)') " species_diffusion_enabled = .false."
  write(unit, '(a)') " barodiffusion_enabled = .false."
  write(unit, '(a)') " transport_cfl = 0.2"
  write(unit, '(a)') " checkpoint_interval_steps = 3"
  write(unit, '(a)') " stop_after_checkpoint = .true."
  write(unit, '(a)') " checkpoint_file = 'configured_eb_3d.chk'"
  write(unit, '(a)') " restart_file = 'restart_eb_3d.chk'"
  write(unit, '(a)') " output_file = 'configured_eb_3d.csv'"
  write(unit, '(a)') "/"
  close(unit)
  call read_reactive_eb_3d_configuration(path, config, ok, message)
  call require(ok .and. config%nx == 12 .and. &
    trim(config%redistribution) == "state_redist" .and. &
    config%chemistry_enabled .and. &
    config%chemistry_relative_tolerance == 1.0e-8_dp .and. &
    config%chemistry_absolute_tolerance == 1.0e-13_dp .and. &
    config%transport_enabled .and. .not. config%viscosity_enabled .and. &
    config%thermal_conduction_enabled .and. &
    .not. config%species_diffusion_enabled .and. &
    .not. config%barodiffusion_enabled .and. &
    config%transport_cfl == 0.2_dp .and. &
    config%checkpoint_interval_steps == 3 .and. &
    config%stop_after_checkpoint .and. &
    trim(config%checkpoint_file) == "configured_eb_3d.chk" .and. &
    trim(config%restart_file) == "restart_eb_3d.chk" .and. &
    trim(config%output_file) == "configured_eb_3d.csv", &
    "public 3D EB namelist read")
  open(newunit=unit, file=path, status="old")
  close(unit, status="delete")

  write(*, '(a)') "test_reactive_eb_3d_config: PASS"

contains

  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      write(*, '(a,1x,a)') trim(label), trim(message)
      error stop label
    end if
  end subroutine require

end program test_reactive_eb_3d_config
