module simulation_config_reactive_eb_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use simulation_config_reactive_2d_mod, only: &
    reactive_2d_config, read_reactive_2d_configuration
  implicit none
  private

  type, public :: reactive_eb_2d_config
    type(reactive_2d_config) :: flow
    character(len=32) :: geometry = "plane"
    real(dp) :: plane_normal_x = 1.0_dp
    real(dp) :: plane_normal_y = 0.0_dp
    real(dp) :: plane_offset = 0.0_dp
    real(dp) :: circle_center_x = 0.0_dp
    real(dp) :: circle_center_y = 0.0_dp
    real(dp) :: circle_radius = 1.0_dp
    logical :: circle_fluid_inside = .false.
    real(dp) :: state_redist_target_volume_fraction = 0.5_dp
    integer :: state_redist_max_order = 0
    character(len=24) :: embedded_wall_kind = "slip_wall"
    character(len=24) :: embedded_wall_thermal = "adiabatic"
    real(dp) :: embedded_wall_temperature = 300.0_dp
    real(dp) :: embedded_wall_velocity(3) = 0.0_dp
  end type reactive_eb_2d_config

  public :: read_reactive_eb_2d_configuration
  public :: reactive_eb_wall_configuration_is_valid
  public :: reactive_eb_wall_transport_is_active

contains

  pure logical function reactive_eb_wall_transport_is_active(config) &
      result(active)
    type(reactive_eb_2d_config), intent(in) :: config

    active = trim(config%embedded_wall_thermal) == "isothermal" .or. &
      trim(config%embedded_wall_kind) == "no_slip_wall"
  end function reactive_eb_wall_transport_is_active

  pure logical function reactive_eb_wall_configuration_is_valid(config) &
      result(valid)
    type(reactive_eb_2d_config), intent(in) :: config

    valid = &
      (trim(config%embedded_wall_kind) == "slip_wall" .or. &
       trim(config%embedded_wall_kind) == "no_slip_wall") .and. &
      (trim(config%embedded_wall_thermal) == "adiabatic" .or. &
       trim(config%embedded_wall_thermal) == "isothermal") .and. &
      ieee_is_finite(config%embedded_wall_temperature) .and. &
      config%embedded_wall_temperature > 0.0_dp .and. &
      all(ieee_is_finite(config%embedded_wall_velocity))
    if (.not. valid) return
    if (trim(config%embedded_wall_thermal) == "isothermal") then
      valid = config%flow%transport_enabled .and. &
        config%flow%thermal_conduction_enabled
      if (.not. valid) return
    end if
    if (trim(config%embedded_wall_kind) == "no_slip_wall") then
      valid = config%flow%transport_enabled .and. &
        config%flow%viscosity_enabled
      if (.not. valid) return
    end if
    if (maxval(abs(config%embedded_wall_velocity)) > 0.0_dp) then
      valid = trim(config%embedded_wall_kind) == "no_slip_wall"
    end if
  end function reactive_eb_wall_configuration_is_valid

  subroutine read_reactive_eb_2d_configuration(path, config, ok, message)
    character(len=*), intent(in) :: path
    type(reactive_eb_2d_config), intent(out) :: config
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    character(len=32) :: geometry
    real(dp) :: plane_normal_x, plane_normal_y, plane_offset
    real(dp) :: circle_center_x, circle_center_y, circle_radius
    real(dp) :: state_redist_target_volume_fraction
    real(dp) :: embedded_wall_temperature, embedded_wall_velocity(3)
    character(len=24) :: embedded_wall_kind, embedded_wall_thermal
    logical :: circle_fluid_inside
    integer :: state_redist_max_order, unit, status
    namelist /embedded_boundary/ geometry, plane_normal_x, plane_normal_y, &
      plane_offset, circle_center_x, circle_center_y, circle_radius, &
      circle_fluid_inside, state_redist_target_volume_fraction, &
      state_redist_max_order, embedded_wall_kind, embedded_wall_thermal, &
      embedded_wall_temperature, embedded_wall_velocity

    config = reactive_eb_2d_config()
    call read_reactive_2d_configuration(path, config%flow, ok, message)
    if (.not. ok) return

    geometry = config%geometry
    plane_normal_x = config%plane_normal_x
    plane_normal_y = config%plane_normal_y
    plane_offset = config%plane_offset
    circle_center_x = config%circle_center_x
    circle_center_y = config%circle_center_y
    circle_radius = config%circle_radius
    circle_fluid_inside = config%circle_fluid_inside
    state_redist_target_volume_fraction = &
      config%state_redist_target_volume_fraction
    state_redist_max_order = config%state_redist_max_order
    embedded_wall_kind = config%embedded_wall_kind
    embedded_wall_thermal = config%embedded_wall_thermal
    embedded_wall_temperature = config%embedded_wall_temperature
    embedded_wall_velocity = config%embedded_wall_velocity

    open(newunit=unit, file=trim(path), status="old", action="read", &
      iostat=status)
    if (status /= 0) then
      ok = .false.
      message = "Could not open reactive EB 2D input"
      return
    end if
    read(unit, nml=embedded_boundary, iostat=status)
    close(unit)
    if (status /= 0) then
      ok = .false.
      message = "Could not parse embedded_boundary namelist"
      return
    end if

    if (.not. all(ieee_is_finite([ &
          plane_normal_x, plane_normal_y, plane_offset, &
          circle_center_x, circle_center_y, circle_radius, &
          state_redist_target_volume_fraction, embedded_wall_temperature, &
          embedded_wall_velocity]))) then
      ok = .false.
      message = "Embedded-boundary parameters must be finite"
      return
    end if
    if (trim(geometry) /= "plane" .and. trim(geometry) /= "circle") then
      ok = .false.
      message = "Embedded-boundary geometry must be plane or circle"
      return
    end if
    if (trim(geometry) == "plane" .and. &
        hypot(plane_normal_x, plane_normal_y) <= 0.0_dp) then
      ok = .false.
      message = "Embedded-boundary plane normal must be nonzero"
      return
    end if
    if (trim(geometry) == "circle" .and. circle_radius <= 0.0_dp) then
      ok = .false.
      message = "Embedded-boundary circle radius must be positive"
      return
    end if
    if (state_redist_target_volume_fraction <= 0.0_dp .or. &
        state_redist_target_volume_fraction > 1.0_dp) then
      ok = .false.
      message = "StateRedist target volume fraction must be in (0,1]"
      return
    end if
    if (state_redist_max_order /= 0 .and. state_redist_max_order /= 2) then
      ok = .false.
      message = "StateRedist max order must be 0 or 2"
      return
    end if
    if ((trim(embedded_wall_kind) /= "slip_wall" .and. &
         trim(embedded_wall_kind) /= "no_slip_wall") .or. &
        (trim(embedded_wall_thermal) /= "adiabatic" .and. &
         trim(embedded_wall_thermal) /= "isothermal") .or. &
        embedded_wall_temperature <= 0.0_dp) then
      ok = .false.
      message = "Embedded-wall kind, thermal mode, or temperature is invalid"
      return
    end if
    if (trim(embedded_wall_thermal) == "isothermal" .and. &
        (.not. config%flow%transport_enabled .or. &
         .not. config%flow%thermal_conduction_enabled)) then
      ok = .false.
      message = "Isothermal embedded walls require thermal transport"
      return
    end if
    if (trim(embedded_wall_kind) == "no_slip_wall" .and. &
        (.not. config%flow%transport_enabled .or. &
         .not. config%flow%viscosity_enabled)) then
      ok = .false.
      message = "No-slip embedded walls require viscous transport"
      return
    end if
    if (maxval(abs(embedded_wall_velocity)) > 0.0_dp .and. &
        trim(embedded_wall_kind) /= "no_slip_wall") then
      ok = .false.
      message = "Moving embedded walls require no-slip transport"
      return
    end if
    if ((trim(config%flow%reconstruction) /= "pcm" .and. &
         trim(config%flow%reconstruction) /= "characteristic_plm") .or. &
        config%flow%use_transverse_correction) then
      ok = .false.
      message = "Reactive EB 2D supports PCM or characteristic PLM without transverse correction"
      return
    end if
    if (trim(config%flow%boundary_x_lower) /= "outflow" .or. &
        trim(config%flow%boundary_x_upper) /= "outflow" .or. &
        trim(config%flow%boundary_y_lower) /= "outflow" .or. &
        trim(config%flow%boundary_y_upper) /= "outflow") then
      ok = .false.
      message = "Reactive EB 2D currently requires outflow domain boundaries"
      return
    end if

    config%geometry = trim(geometry)
    config%plane_normal_x = plane_normal_x
    config%plane_normal_y = plane_normal_y
    config%plane_offset = plane_offset
    config%circle_center_x = circle_center_x
    config%circle_center_y = circle_center_y
    config%circle_radius = circle_radius
    config%circle_fluid_inside = circle_fluid_inside
    config%state_redist_target_volume_fraction = &
      state_redist_target_volume_fraction
    config%state_redist_max_order = state_redist_max_order
    config%embedded_wall_kind = trim(embedded_wall_kind)
    config%embedded_wall_thermal = trim(embedded_wall_thermal)
    config%embedded_wall_temperature = embedded_wall_temperature
    config%embedded_wall_velocity = embedded_wall_velocity
    message = ""
    ok = .true.
  end subroutine read_reactive_eb_2d_configuration

end module simulation_config_reactive_eb_2d_mod
