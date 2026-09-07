program test_amr_reactive_3d_config
  use simulation_config_reactive_3d_mod, only: &
    reactive_3d_config, validate_reactive_3d_configuration
  use simulation_config_amr_reactive_3d_mod, only: &
    amr_reactive_3d_config, validate_amr_reactive_3d_configuration
  implicit none

  type(reactive_3d_config) :: base
  type(amr_reactive_3d_config) :: amr
  character(len=256) :: message
  logical :: ok

  base = reactive_3d_config()
  base%nx = 8
  base%ny = 8
  base%nz = 8
  amr = amr_reactive_3d_config()
  call validate_amr_reactive_3d_configuration(base, amr, ok, message)
  call require(ok, "default ratio-two interior patch")
  base%reconstruction = "characteristic_plm"
  base%limiter = "minmod"
  call validate_amr_reactive_3d_configuration(base, amr, ok, message)
  call require(ok, "characteristic PLM AMR configuration")
  base%reconstruction = "invalid"
  call validate_amr_reactive_3d_configuration(base, amr, ok, message)
  call require(.not. ok, "invalid reconstruction rejection")
  base%reconstruction = "characteristic_plm"
  base%limiter = "invalid"
  call validate_amr_reactive_3d_configuration(base, amr, ok, message)
  call require(.not. ok, "invalid limiter rejection")
  base%reconstruction = "pcm"
  base%limiter = "mc"

  base%chemistry_enabled = .true.
  call validate_amr_reactive_3d_configuration(base, amr, ok, message)
  call require(ok, "fixed-mechanism chemistry configuration")
  amr%checkpoint_interval_steps = 4
  call validate_amr_reactive_3d_configuration(base, amr, ok, message)
  call require(.not. ok, "chemistry checkpoint rejection")
  amr = amr_reactive_3d_config()
  amr%restart_file = "chemistry.chk"
  call validate_amr_reactive_3d_configuration(base, amr, ok, message)
  call require(.not. ok, "chemistry restart rejection")
  amr = amr_reactive_3d_config()
  base%chemistry_enabled = .false.
  base%transport_enabled = .true.
  call validate_amr_reactive_3d_configuration(base, amr, ok, message)
  call require(ok, "static periodic transport support")
  base%barodiffusion_enabled = .true.
  base%species_diffusion_enabled = .false.
  call validate_reactive_3d_configuration(base, ok, message)
  call require(.not. ok, "barodiffusion without species diffusion rejection")
  call require(trim(message) == &
    "Reactive 3D barodiffusion requires species diffusion", &
    "barodiffusion dependency diagnostic")
  base%species_diffusion_enabled = .true.
  amr%checkpoint_interval_steps = 1
  call validate_amr_reactive_3d_configuration(base, amr, ok, message)
  call require(ok, "transport checkpoint acceptance")
  amr = amr_reactive_3d_config()
  amr%restart_file = "transport.chk"
  call validate_amr_reactive_3d_configuration(base, amr, ok, message)
  call require(ok, "transport restart acceptance")
  amr = amr_reactive_3d_config()
  base%transport_enabled = .false.
  base%thermo_model = "selected"
  call validate_amr_reactive_3d_configuration(base, amr, ok, message)
  call require(.not. ok, "selected mechanism default rejection")
  call validate_amr_reactive_3d_configuration( &
    base, amr, ok, message, allow_selected=.true.)
  call require(ok, "selected mechanism explicit opt-in")
  base%chemistry_enabled = .true.
  amr%checkpoint_interval_steps = 1
  call validate_amr_reactive_3d_configuration( &
    base, amr, ok, message, allow_selected=.true.)
  call require(ok, "selected chemistry checkpoint support")
  amr = amr_reactive_3d_config()
  amr%checkpoint_interval_steps = 1
  amr%stop_after_checkpoint = .true.
  call validate_amr_reactive_3d_configuration( &
    base, amr, ok, message, allow_selected=.true.)
  call require(ok, "selected chemistry checkpoint stop support")
  amr = amr_reactive_3d_config()
  amr%restart_file = "selected_restart.chk"
  call validate_amr_reactive_3d_configuration( &
    base, amr, ok, message, allow_selected=.true.)
  call require(ok, "selected chemistry restart support")

  base%transport_enabled = .true.
  amr = amr_reactive_3d_config()
  amr%checkpoint_interval_steps = 1
  call validate_amr_reactive_3d_configuration( &
    base, amr, ok, message, allow_selected=.true.)
  call require(ok, "selected chemistry and transport checkpoint acceptance")
  amr = amr_reactive_3d_config()
  amr%restart_file = "selected_transport_restart.chk"
  call validate_amr_reactive_3d_configuration( &
    base, amr, ok, message, allow_selected=.true.)
  call require(ok, "selected chemistry and transport restart acceptance")

  amr = amr_reactive_3d_config()
  base%chemistry_enabled = .false.
  base%transport_enabled = .false.
  base%thermo_model = "full_h2o2"
  base%boundary_condition = "outflow"
  call validate_amr_reactive_3d_configuration(base, amr, ok, message)
  call require(.not. ok, "nonperiodic rejection")
  base%boundary_condition = "periodic"

  amr%coarse_i_lower = 1
  call validate_amr_reactive_3d_configuration(base, amr, ok, message)
  call require(.not. ok, "boundary-touching patch rejection")
  amr = amr_reactive_3d_config()
  amr%refinement_ratio = 1
  call validate_amr_reactive_3d_configuration(base, amr, ok, message)
  call require(.not. ok, "unit refinement rejection")
  amr = amr_reactive_3d_config()
  amr%fine_output_file = amr%coarse_output_file
  call validate_amr_reactive_3d_configuration(base, amr, ok, message)
  call require(.not. ok, "duplicate output rejection")
  amr = amr_reactive_3d_config()
  amr%checkpoint_interval_steps = -1
  call validate_amr_reactive_3d_configuration(base, amr, ok, message)
  call require(.not. ok, "negative checkpoint interval rejection")
  amr = amr_reactive_3d_config()
  amr%stop_after_checkpoint = .true.
  call validate_amr_reactive_3d_configuration(base, amr, ok, message)
  call require(.not. ok, "checkpoint stop without interval rejection")
  amr = amr_reactive_3d_config()
  amr%checkpoint_interval_steps = 4
  amr%checkpoint_file = amr%coarse_output_file
  call validate_amr_reactive_3d_configuration(base, amr, ok, message)
  call require(.not. ok, "checkpoint and output collision rejection")
  amr = amr_reactive_3d_config()
  amr%restart_file = amr%fine_output_file
  call validate_amr_reactive_3d_configuration(base, amr, ok, message)
  call require(.not. ok, "restart input and output collision rejection")
  write(*, '(a)') "test_amr_reactive_3d_config: PASS"

contains

  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      write(*, '(a)') "FAIL: " // trim(label)
      error stop 1
    end if
  end subroutine require

end program test_amr_reactive_3d_config
