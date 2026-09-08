module simulation_config_reactive_eb_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use selected_composition_mod, only: &
    selected_composition_max_species, selected_composition_name_length, &
    validate_selected_composition_fields, resolve_selected_composition
  implicit none
  private

  integer, parameter, public :: reactive_eb_3d_max_species = &
    selected_composition_max_species

  type, public :: reactive_eb_3d_config
    integer :: nx = 10
    integer :: ny = 8
    integer :: nz = 6
    integer :: maximum_steps = 100
    real(dp) :: x_lower = 0.0_dp
    real(dp) :: x_upper = 1.0_dp
    real(dp) :: y_lower = 0.0_dp
    real(dp) :: y_upper = 1.0_dp
    real(dp) :: z_lower = 0.0_dp
    real(dp) :: z_upper = 1.0_dp
    real(dp) :: final_time = 5.0e-5_dp
    real(dp) :: cfl = 0.5_dp
    character(len=8) :: plane_axis = "x"
    real(dp) :: plane_position = 0.395_dp
    character(len=32) :: riemann_solver = "rusanov"
    character(len=32) :: redistribution = "state_redist"
    character(len=32) :: thermo_model = "elementary"
    real(dp) :: state_redist_target_volume_fraction = 0.5_dp
    logical :: chemistry_enabled = .false.
    real(dp) :: chemistry_relative_tolerance = 2.0e-7_dp
    real(dp) :: chemistry_absolute_tolerance = 1.0e-12_dp
    logical :: transport_enabled = .false.
    logical :: viscosity_enabled = .true.
    logical :: thermal_conduction_enabled = .true.
    logical :: species_diffusion_enabled = .true.
    logical :: barodiffusion_enabled = .true.
    real(dp) :: transport_cfl = 0.35_dp
    real(dp) :: regular_density = 0.31_dp
    real(dp) :: cut_density = 0.93_dp
    real(dp) :: initial_pressure = 135000.0_dp
    real(dp) :: initial_velocity_x = 0.0_dp
    real(dp) :: initial_velocity_y = 0.0_dp
    real(dp) :: initial_velocity_z = 0.0_dp
    integer :: composition_count = 0
    character(len=selected_composition_name_length) :: &
      composition_species(reactive_eb_3d_max_species) = ""
    real(dp) :: &
      composition_mole_fractions(reactive_eb_3d_max_species) = 0.0_dp
    integer :: checkpoint_interval_steps = 0
    logical :: stop_after_checkpoint = .false.
    character(len=512) :: checkpoint_file = "reactive_eb_3d.chk"
    character(len=512) :: restart_file = ""
    character(len=512) :: output_file = "reactive_eb_3d.csv"
  end type reactive_eb_3d_config

  public :: read_reactive_eb_3d_configuration
  public :: validate_reactive_eb_3d_configuration
  public :: resolve_reactive_eb_3d_selected_composition

contains

  subroutine read_reactive_eb_3d_configuration( &
      path, config, ok, message, allow_selected)
    character(len=*), intent(in) :: path
    type(reactive_eb_3d_config), intent(out) :: config
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    logical, intent(in), optional :: allow_selected

    integer :: nx, ny, nz, maximum_steps, checkpoint_interval_steps
    integer :: unit, io_status
    real(dp) :: x_lower, x_upper, y_lower, y_upper, z_lower, z_upper
    real(dp) :: final_time, cfl, plane_position
    real(dp) :: state_redist_target_volume_fraction
    real(dp) :: chemistry_relative_tolerance, chemistry_absolute_tolerance
    real(dp) :: transport_cfl
    real(dp) :: regular_density, cut_density, initial_pressure
    real(dp) :: initial_velocity_x, initial_velocity_y, initial_velocity_z
    character(len=8) :: plane_axis
    character(len=32) :: riemann_solver, redistribution, thermo_model
    character(len=512) :: checkpoint_file, restart_file, output_file
    integer :: composition_count
    character(len=selected_composition_name_length) :: &
      composition_species(reactive_eb_3d_max_species)
    real(dp) :: composition_mole_fractions(reactive_eb_3d_max_species)
    logical :: chemistry_enabled, transport_enabled, viscosity_enabled
    logical :: thermal_conduction_enabled, species_diffusion_enabled
    logical :: barodiffusion_enabled, stop_after_checkpoint
    logical :: selected_allowed

    namelist /reactive_eb_3d/ &
      nx, ny, nz, maximum_steps, &
      x_lower, x_upper, y_lower, y_upper, z_lower, z_upper, &
      final_time, cfl, plane_axis, plane_position, riemann_solver, &
      redistribution, thermo_model, state_redist_target_volume_fraction, &
      chemistry_enabled, chemistry_relative_tolerance, &
      chemistry_absolute_tolerance, &
      transport_enabled, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, transport_cfl, &
      regular_density, cut_density, initial_pressure, &
      initial_velocity_x, initial_velocity_y, initial_velocity_z, &
      composition_count, composition_species, composition_mole_fractions, &
      checkpoint_interval_steps, stop_after_checkpoint, &
      checkpoint_file, restart_file, &
      output_file

    config = reactive_eb_3d_config()
    nx = config%nx
    ny = config%ny
    nz = config%nz
    maximum_steps = config%maximum_steps
    x_lower = config%x_lower
    x_upper = config%x_upper
    y_lower = config%y_lower
    y_upper = config%y_upper
    z_lower = config%z_lower
    z_upper = config%z_upper
    final_time = config%final_time
    cfl = config%cfl
    plane_axis = config%plane_axis
    plane_position = config%plane_position
    riemann_solver = config%riemann_solver
    redistribution = config%redistribution
    thermo_model = config%thermo_model
    state_redist_target_volume_fraction = &
      config%state_redist_target_volume_fraction
    chemistry_enabled = config%chemistry_enabled
    chemistry_relative_tolerance = config%chemistry_relative_tolerance
    chemistry_absolute_tolerance = config%chemistry_absolute_tolerance
    transport_enabled = config%transport_enabled
    viscosity_enabled = config%viscosity_enabled
    thermal_conduction_enabled = config%thermal_conduction_enabled
    species_diffusion_enabled = config%species_diffusion_enabled
    barodiffusion_enabled = config%barodiffusion_enabled
    transport_cfl = config%transport_cfl
    regular_density = config%regular_density
    cut_density = config%cut_density
    initial_pressure = config%initial_pressure
    initial_velocity_x = config%initial_velocity_x
    initial_velocity_y = config%initial_velocity_y
    initial_velocity_z = config%initial_velocity_z
    composition_count = config%composition_count
    composition_species = config%composition_species
    composition_mole_fractions = config%composition_mole_fractions
    checkpoint_interval_steps = config%checkpoint_interval_steps
    stop_after_checkpoint = config%stop_after_checkpoint
    checkpoint_file = config%checkpoint_file
    restart_file = config%restart_file
    output_file = config%output_file
    selected_allowed = .false.
    if (present(allow_selected)) selected_allowed = allow_selected

    open(newunit=unit, file=trim(path), status="old", action="read", &
      iostat=io_status)
    if (io_status /= 0) then
      ok = .false.
      write(message, '(a,1x,a)') "Could not open input file:", trim(path)
      return
    end if
    read(unit, nml=reactive_eb_3d, iostat=io_status)
    close(unit)
    if (io_status /= 0) then
      ok = .false.
      write(message, '(a,1x,a)') &
        "Could not read &reactive_eb_3d from:", trim(path)
      return
    end if

    config%nx = nx
    config%ny = ny
    config%nz = nz
    config%maximum_steps = maximum_steps
    config%x_lower = x_lower
    config%x_upper = x_upper
    config%y_lower = y_lower
    config%y_upper = y_upper
    config%z_lower = z_lower
    config%z_upper = z_upper
    config%final_time = final_time
    config%cfl = cfl
    config%plane_axis = trim(plane_axis)
    config%plane_position = plane_position
    config%riemann_solver = trim(riemann_solver)
    config%redistribution = trim(redistribution)
    config%thermo_model = trim(thermo_model)
    config%state_redist_target_volume_fraction = &
      state_redist_target_volume_fraction
    config%chemistry_enabled = chemistry_enabled
    config%chemistry_relative_tolerance = chemistry_relative_tolerance
    config%chemistry_absolute_tolerance = chemistry_absolute_tolerance
    config%transport_enabled = transport_enabled
    config%viscosity_enabled = viscosity_enabled
    config%thermal_conduction_enabled = thermal_conduction_enabled
    config%species_diffusion_enabled = species_diffusion_enabled
    config%barodiffusion_enabled = barodiffusion_enabled
    config%transport_cfl = transport_cfl
    config%regular_density = regular_density
    config%cut_density = cut_density
    config%initial_pressure = initial_pressure
    config%initial_velocity_x = initial_velocity_x
    config%initial_velocity_y = initial_velocity_y
    config%initial_velocity_z = initial_velocity_z
    config%composition_count = composition_count
    config%composition_species = composition_species
    config%composition_mole_fractions = composition_mole_fractions
    config%checkpoint_interval_steps = checkpoint_interval_steps
    config%stop_after_checkpoint = stop_after_checkpoint
    config%checkpoint_file = trim(checkpoint_file)
    config%restart_file = trim(restart_file)
    config%output_file = trim(output_file)
    call validate_reactive_eb_3d_configuration( &
      config, ok, message, allow_selected=selected_allowed)
  end subroutine read_reactive_eb_3d_configuration

  pure subroutine validate_reactive_eb_3d_configuration( &
      config, ok, message, allow_selected)
    type(reactive_eb_3d_config), intent(in) :: config
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    logical, intent(in), optional :: allow_selected

    real(dp) :: axis_lower, axis_upper, spacing, coordinate, tolerance
    real(dp) :: cut_volume_fraction
    logical :: selected_allowed, selected_model

    ok = .false.
    message = ""
    selected_allowed = .false.
    if (present(allow_selected)) selected_allowed = allow_selected
    selected_model = trim(config%thermo_model) == "selected"
    if (selected_model .and. .not. selected_allowed) then
      message = "Reactive EB 3D selected thermo requires explicit opt-in"
      return
    else if (.not. selected_model .and. &
             trim(config%thermo_model) /= "elementary") then
      message = "Reactive EB 3D thermo model is unsupported"
      return
    else if (.not. selected_model .and. config%composition_count /= 0) then
      message = "Reactive EB 3D composition fields require selected thermo"
      return
    else if (selected_model) then
      call validate_selected_composition_fields( &
        "Reactive EB 3D", config%composition_count, &
        config%composition_species, config%composition_mole_fractions, &
        ok, message)
      if (.not. ok) return
      ok = .false.
    end if
    if (min(config%nx, config%ny, config%nz) < 2) then
      message = "Reactive EB 3D grid extents must all be at least 2"
      return
    else if (config%maximum_steps <= 0) then
      message = "Reactive EB 3D maximum_steps must be positive"
      return
    else if (.not. all(ieee_is_finite([ &
        config%x_lower, config%x_upper, config%y_lower, config%y_upper, &
        config%z_lower, config%z_upper, config%final_time, config%cfl, &
        config%plane_position, &
        config%state_redist_target_volume_fraction, &
        config%chemistry_relative_tolerance, &
        config%chemistry_absolute_tolerance, &
        config%transport_cfl, &
        config%regular_density, config%cut_density, &
        config%initial_pressure, config%initial_velocity_x, &
        config%initial_velocity_y, config%initial_velocity_z]))) then
      message = "Reactive EB 3D scalar inputs must be finite"
      return
    else if (config%x_upper <= config%x_lower .or. &
             config%y_upper <= config%y_lower .or. &
             config%z_upper <= config%z_lower) then
      message = "Reactive EB 3D domain maxima must exceed minima"
      return
    else if (config%final_time <= 0.0_dp) then
      message = "Reactive EB 3D final_time must be positive"
      return
    else if (config%cfl <= 0.0_dp .or. config%cfl > 1.0_dp) then
      message = "Reactive EB 3D cfl must be in (0,1]"
      return
    else if (config%state_redist_target_volume_fraction <= 0.0_dp .or. &
             config%state_redist_target_volume_fraction > 1.0_dp) then
      message = "Reactive EB 3D StateRedist target must be in (0,1]"
      return
    else if (config%chemistry_relative_tolerance <= 0.0_dp .or. &
             config%chemistry_absolute_tolerance <= 0.0_dp) then
      message = "Reactive EB 3D chemistry tolerances must be positive"
      return
    else if (config%transport_cfl <= 0.0_dp .or. &
             config%transport_cfl > 0.5_dp) then
      message = "Reactive EB 3D transport_cfl must be in (0,0.5]"
      return
    else if (config%transport_enabled .and. &
             trim(config%redistribution) /= "state_redist") then
      message = "Reactive EB 3D transport requires StateRedist"
      return
    else if (config%transport_enabled .and. &
             .not. (config%viscosity_enabled .or. &
               config%thermal_conduction_enabled .or. &
               config%species_diffusion_enabled)) then
      message = "Reactive EB 3D transport requires an enabled process"
      return
    else if (config%transport_enabled .and. &
             config%barodiffusion_enabled .and. &
             .not. config%species_diffusion_enabled) then
      message = "Reactive EB 3D barodiffusion requires species diffusion"
      return
    else if (config%regular_density <= 0.0_dp .or. &
             config%cut_density <= 0.0_dp .or. &
             config%initial_pressure <= 0.0_dp) then
      message = "Reactive EB 3D densities and pressure must be positive"
      return
    else if (config%checkpoint_interval_steps < 0) then
      message = "Reactive EB 3D checkpoint interval must be nonnegative"
      return
    else if (config%stop_after_checkpoint .and. &
             config%checkpoint_interval_steps <= 0) then
      message = "Reactive EB 3D checkpoint stop requires an interval"
      return
    else if (config%checkpoint_interval_steps > 0 .and. &
             len_trim(config%checkpoint_file) == 0) then
      message = "Reactive EB 3D checkpoint_file must not be empty"
      return
    else if (trim(config%riemann_solver) /= "rusanov" .and. &
             trim(config%riemann_solver) /= "hllc" .and. &
             trim(config%riemann_solver) /= "pelec") then
      message = "Reactive EB 3D Riemann solver is unsupported"
      return
    else if (trim(config%redistribution) /= "flux_redist" .and. &
             trim(config%redistribution) /= "state_redist") then
      message = "Reactive EB 3D redistribution is unsupported"
      return
    else if (len_trim(config%output_file) == 0) then
      message = "Reactive EB 3D output_file must not be empty"
      return
    else if (config%checkpoint_interval_steps > 0 .and. &
             paths_alias_lexically( &
               config%checkpoint_file, config%output_file)) then
      message = "Reactive EB 3D checkpoint and output paths must differ"
      return
    else if (len_trim(config%restart_file) > 0 .and. &
             paths_alias_lexically( &
               config%restart_file, config%output_file)) then
      message = "Reactive EB 3D restart and output paths must differ"
      return
    end if

    select case (trim(config%plane_axis))
    case ("x")
      axis_lower = config%x_lower
      axis_upper = config%x_upper
      spacing = (axis_upper - axis_lower) / real(config%nx, dp)
    case ("y")
      axis_lower = config%y_lower
      axis_upper = config%y_upper
      spacing = (axis_upper - axis_lower) / real(config%ny, dp)
    case ("z")
      axis_lower = config%z_lower
      axis_upper = config%z_upper
      spacing = (axis_upper - axis_lower) / real(config%nz, dp)
    case default
      message = "Reactive EB 3D plane_axis must be x, y, or z"
      return
    end select
    if (config%plane_position <= axis_lower .or. &
        config%plane_position >= axis_upper - spacing) then
      message = "Reactive EB 3D plane needs a cut cell and regular receiver"
      return
    end if
    coordinate = (config%plane_position - axis_lower) / spacing
    tolerance = 512.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, abs(coordinate))
    if (abs(coordinate - anint(coordinate)) <= tolerance) then
      message = "Reactive EB 3D plane must not align with a grid face"
      return
    end if
    cut_volume_fraction = real(ceiling(coordinate), dp) - coordinate
    if (trim(config%redistribution) == "state_redist" .and. &
        config%state_redist_target_volume_fraction <= &
          cut_volume_fraction + tolerance) then
      message = "Reactive EB 3D StateRedist target must exceed cut fraction"
      return
    end if
    ok = .true.
  end subroutine validate_reactive_eb_3d_configuration

  subroutine resolve_reactive_eb_3d_selected_composition( &
      config, species, mole_fractions, ok, message)
    type(reactive_eb_3d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(out) :: mole_fractions(:)
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    mole_fractions = 0.0_dp
    if (trim(config%thermo_model) /= "selected") then
      ok = .false.
      message = &
        "Reactive EB 3D selected composition requires thermo_model='selected'"
      return
    end if
    call resolve_selected_composition( &
      "Reactive EB 3D", config%composition_count, &
      config%composition_species, config%composition_mole_fractions, &
      species, mole_fractions, ok, message)
  end subroutine resolve_reactive_eb_3d_selected_composition

  pure logical function paths_alias_lexically(left, right) result(aliases)
    character(len=*), intent(in) :: left, right

    aliases = normalized_path(left) == normalized_path(right)
  end function paths_alias_lexically

  pure function normalized_path(path) result(normalized)
    character(len=*), intent(in) :: path
    character(len=1024) :: normalized
    character(len=1024) :: input, component
    integer :: component_starts(256), component_ends(256)
    integer :: input_length, start, finish, component_count, index
    logical :: absolute

    normalized = ""
    input = trim(path)
    input_length = len_trim(input)
    if (input_length == 0) return
    absolute = input(1:1) == "/"
    component_count = 0
    start = 1
    do while (start <= input_length)
      do while (start <= input_length)
        if (input(start:start) /= "/") exit
        start = start + 1
      end do
      if (start > input_length) exit
      finish = start
      do while (finish <= input_length)
        if (input(finish:finish) == "/") exit
        finish = finish + 1
      end do
      component = input(start:finish - 1)
      select case (trim(component))
      case ("", ".")
        continue
      case ("..")
        if (component_count > 0) then
          if (trim(input(component_starts(component_count): &
              component_ends(component_count))) /= "..") then
            component_count = component_count - 1
          else if (.not. absolute) then
            component_count = component_count + 1
            component_starts(component_count) = start
            component_ends(component_count) = finish - 1
          end if
        else if (.not. absolute) then
          component_count = component_count + 1
          component_starts(component_count) = start
          component_ends(component_count) = finish - 1
        end if
      case default
        component_count = component_count + 1
        component_starts(component_count) = start
        component_ends(component_count) = finish - 1
      end select
      start = finish + 1
    end do

    if (absolute) normalized = "/"
    do index = 1, component_count
      if (len_trim(normalized) > 0 .and. &
          trim(normalized) /= "/") normalized = trim(normalized) // "/"
      normalized = trim(normalized) // &
        trim(input(component_starts(index):component_ends(index)))
    end do
    if (len_trim(normalized) == 0) then
      if (absolute) then
        normalized = "/"
      else
        normalized = "."
      end if
    end if
  end function normalized_path

end module simulation_config_reactive_eb_3d_mod
