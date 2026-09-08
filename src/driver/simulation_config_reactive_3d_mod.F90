module simulation_config_reactive_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use selected_composition_mod, only: &
    selected_composition_max_species, selected_composition_name_length, &
    validate_selected_composition_fields, resolve_selected_composition
  implicit none
  private

  integer, parameter, public :: reactive_3d_max_species = &
    selected_composition_max_species

  type, public :: reactive_3d_config
    integer :: nx = 12
    integer :: ny = 12
    integer :: nz = 12
    integer :: maximum_steps = 100000
    real(dp) :: x_lower = 0.0_dp
    real(dp) :: x_upper = 0.01_dp
    real(dp) :: y_lower = 0.0_dp
    real(dp) :: y_upper = 0.01_dp
    real(dp) :: z_lower = 0.0_dp
    real(dp) :: z_upper = 0.01_dp
    real(dp) :: final_time = 1.0e-6_dp
    real(dp) :: cfl = 0.30_dp
    character(len=32) :: problem = "entropy_wave"
    character(len=32) :: reconstruction = "pcm"
    character(len=32) :: limiter = "mc"
    character(len=32) :: thermo_model = "full_h2o2"
    character(len=32) :: riemann_solver = "rusanov"
    character(len=32) :: boundary_condition = "periodic"
    logical :: chemistry_enabled = .false.
    real(dp) :: chemistry_relative_tolerance = 2.0e-7_dp
    real(dp) :: chemistry_absolute_tolerance = 1.0e-12_dp
    logical :: transport_enabled = .false.
    logical :: viscosity_enabled = .true.
    logical :: thermal_conduction_enabled = .true.
    logical :: species_diffusion_enabled = .true.
    logical :: barodiffusion_enabled = .true.
    real(dp) :: transport_cfl = 0.35_dp
    real(dp) :: initial_temperature = 1000.0_dp
    real(dp) :: initial_pressure = 101325.0_dp
    real(dp) :: initial_velocity_x = 300.0_dp
    real(dp) :: initial_velocity_y = 200.0_dp
    real(dp) :: initial_velocity_z = -100.0_dp
    real(dp) :: density_wave_amplitude = 0.08_dp
    integer :: wave_number_x = 1
    integer :: wave_number_y = 1
    integer :: wave_number_z = 1
    real(dp) :: hotspot_temperature_rise = 250.0_dp
    real(dp) :: hotspot_center_x = 0.005_dp
    real(dp) :: hotspot_center_y = 0.005_dp
    real(dp) :: hotspot_center_z = 0.005_dp
    real(dp) :: hotspot_width = 0.0012_dp
    real(dp) :: x_h2 = 0.29570_dp
    real(dp) :: x_h = 1.0e-5_dp
    real(dp) :: x_o = 1.0e-5_dp
    real(dp) :: x_o2 = 0.14784_dp
    real(dp) :: x_oh = 1.0e-5_dp
    real(dp) :: x_h2o = 0.0_dp
    real(dp) :: x_ho2 = 0.0_dp
    real(dp) :: x_h2o2 = 0.0_dp
    real(dp) :: x_ar = 0.0_dp
    real(dp) :: x_n2 = 0.55643_dp
    integer :: composition_count = 0
    character(len=selected_composition_name_length) :: &
      composition_species(reactive_3d_max_species) = ""
    real(dp) :: &
      composition_mole_fractions(reactive_3d_max_species) = 0.0_dp
    character(len=512) :: output_file = "reactive_entropy_wave_3d.csv"
  end type reactive_3d_config

  public :: read_reactive_3d_configuration
  public :: validate_reactive_3d_configuration
  public :: reactive_3d_mole_fractions
  public :: resolve_reactive_3d_selected_composition

contains

  subroutine read_reactive_3d_configuration( &
      path, config, ok, message, allow_selected)
    character(len=*), intent(in) :: path
    type(reactive_3d_config), intent(out) :: config
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    logical, intent(in), optional :: allow_selected

    integer :: nx, ny, nz, maximum_steps, unit, io_status
    integer :: wave_number_x, wave_number_y, wave_number_z
    real(dp) :: x_lower, x_upper, y_lower, y_upper, z_lower, z_upper
    real(dp) :: final_time, cfl, initial_temperature, initial_pressure
    real(dp) :: chemistry_relative_tolerance, chemistry_absolute_tolerance
    real(dp) :: transport_cfl
    real(dp) :: initial_velocity_x, initial_velocity_y, initial_velocity_z
    real(dp) :: density_wave_amplitude
    real(dp) :: hotspot_temperature_rise, hotspot_center_x
    real(dp) :: hotspot_center_y, hotspot_center_z, hotspot_width
    real(dp) :: x_h2, x_h, x_o, x_o2, x_oh, x_h2o
    real(dp) :: x_ho2, x_h2o2, x_ar, x_n2
    real(dp) :: composition_mole_fractions(reactive_3d_max_species)
    integer :: composition_count
    character(len=selected_composition_name_length) :: &
      composition_species(reactive_3d_max_species)
    character(len=32) :: problem, reconstruction, limiter
    character(len=32) :: thermo_model, riemann_solver
    character(len=32) :: boundary_condition
    character(len=512) :: output_file
    logical :: chemistry_enabled, transport_enabled, viscosity_enabled
    logical :: thermal_conduction_enabled, species_diffusion_enabled
    logical :: barodiffusion_enabled
    logical :: selected_allowed

    namelist /reactive_3d/ &
      nx, ny, nz, maximum_steps, &
      x_lower, x_upper, y_lower, y_upper, z_lower, z_upper, &
      final_time, cfl, problem, reconstruction, limiter, &
      thermo_model, riemann_solver, &
      boundary_condition, chemistry_enabled, chemistry_relative_tolerance, &
      chemistry_absolute_tolerance, transport_enabled, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, transport_cfl, &
      initial_temperature, initial_pressure, &
      initial_velocity_x, initial_velocity_y, initial_velocity_z, &
      density_wave_amplitude, wave_number_x, wave_number_y, wave_number_z, &
      hotspot_temperature_rise, hotspot_center_x, hotspot_center_y, &
      hotspot_center_z, hotspot_width, &
      x_h2, x_h, x_o, x_o2, x_oh, x_h2o, x_ho2, x_h2o2, x_ar, x_n2, &
      composition_count, composition_species, composition_mole_fractions, &
      output_file

    config = reactive_3d_config()
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
    problem = config%problem
    reconstruction = config%reconstruction
    limiter = config%limiter
    thermo_model = config%thermo_model
    riemann_solver = config%riemann_solver
    boundary_condition = config%boundary_condition
    chemistry_enabled = config%chemistry_enabled
    chemistry_relative_tolerance = config%chemistry_relative_tolerance
    chemistry_absolute_tolerance = config%chemistry_absolute_tolerance
    transport_enabled = config%transport_enabled
    viscosity_enabled = config%viscosity_enabled
    thermal_conduction_enabled = config%thermal_conduction_enabled
    species_diffusion_enabled = config%species_diffusion_enabled
    barodiffusion_enabled = config%barodiffusion_enabled
    transport_cfl = config%transport_cfl
    initial_temperature = config%initial_temperature
    initial_pressure = config%initial_pressure
    initial_velocity_x = config%initial_velocity_x
    initial_velocity_y = config%initial_velocity_y
    initial_velocity_z = config%initial_velocity_z
    density_wave_amplitude = config%density_wave_amplitude
    wave_number_x = config%wave_number_x
    wave_number_y = config%wave_number_y
    wave_number_z = config%wave_number_z
    hotspot_temperature_rise = config%hotspot_temperature_rise
    hotspot_center_x = config%hotspot_center_x
    hotspot_center_y = config%hotspot_center_y
    hotspot_center_z = config%hotspot_center_z
    hotspot_width = config%hotspot_width
    x_h2 = config%x_h2
    x_h = config%x_h
    x_o = config%x_o
    x_o2 = config%x_o2
    x_oh = config%x_oh
    x_h2o = config%x_h2o
    x_ho2 = config%x_ho2
    x_h2o2 = config%x_h2o2
    x_ar = config%x_ar
    x_n2 = config%x_n2
    composition_count = config%composition_count
    composition_species = config%composition_species
    composition_mole_fractions = config%composition_mole_fractions
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
    read(unit, nml=reactive_3d, iostat=io_status)
    close(unit)
    if (io_status /= 0) then
      ok = .false.
      write(message, '(a,1x,a)') &
        "Could not read &reactive_3d from:", trim(path)
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
    config%problem = trim(problem)
    config%reconstruction = trim(reconstruction)
    config%limiter = trim(limiter)
    config%thermo_model = trim(thermo_model)
    config%riemann_solver = trim(riemann_solver)
    config%boundary_condition = trim(boundary_condition)
    config%chemistry_enabled = chemistry_enabled
    config%chemistry_relative_tolerance = chemistry_relative_tolerance
    config%chemistry_absolute_tolerance = chemistry_absolute_tolerance
    config%transport_enabled = transport_enabled
    config%viscosity_enabled = viscosity_enabled
    config%thermal_conduction_enabled = thermal_conduction_enabled
    config%species_diffusion_enabled = species_diffusion_enabled
    config%barodiffusion_enabled = barodiffusion_enabled
    config%transport_cfl = transport_cfl
    config%initial_temperature = initial_temperature
    config%initial_pressure = initial_pressure
    config%initial_velocity_x = initial_velocity_x
    config%initial_velocity_y = initial_velocity_y
    config%initial_velocity_z = initial_velocity_z
    config%density_wave_amplitude = density_wave_amplitude
    config%wave_number_x = wave_number_x
    config%wave_number_y = wave_number_y
    config%wave_number_z = wave_number_z
    config%hotspot_temperature_rise = hotspot_temperature_rise
    config%hotspot_center_x = hotspot_center_x
    config%hotspot_center_y = hotspot_center_y
    config%hotspot_center_z = hotspot_center_z
    config%hotspot_width = hotspot_width
    config%x_h2 = x_h2
    config%x_h = x_h
    config%x_o = x_o
    config%x_o2 = x_o2
    config%x_oh = x_oh
    config%x_h2o = x_h2o
    config%x_ho2 = x_ho2
    config%x_h2o2 = x_h2o2
    config%x_ar = x_ar
    config%x_n2 = x_n2
    config%composition_count = composition_count
    config%composition_species = composition_species
    config%composition_mole_fractions = composition_mole_fractions
    config%output_file = trim(output_file)
    call validate_reactive_3d_configuration( &
      config, ok, message, allow_selected=selected_allowed)
  end subroutine read_reactive_3d_configuration

  pure subroutine validate_reactive_3d_configuration( &
      config, ok, message, allow_selected)
    type(reactive_3d_config), intent(in) :: config
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    logical, intent(in), optional :: allow_selected

    real(dp) :: mole_sum
    logical :: selected_allowed, selected_model

    ok = .false.
    message = ""
    selected_allowed = .false.
    if (present(allow_selected)) selected_allowed = allow_selected
    selected_model = trim(config%thermo_model) == "selected"
    if (config%nx < 4 .or. config%ny < 4 .or. config%nz < 4) then
      message = "Reactive 3D grid extents must all be at least 4"
    else if (config%maximum_steps <= 0) then
      message = "Reactive 3D maximum_steps must be positive"
    else if (config%x_upper <= config%x_lower .or. &
             config%y_upper <= config%y_lower .or. &
             config%z_upper <= config%z_lower) then
      message = "Reactive 3D domain maxima must exceed minima"
    else if (.not. all(ieee_is_finite([ &
        config%x_lower, config%x_upper, config%y_lower, config%y_upper, &
        config%z_lower, config%z_upper, config%final_time, config%cfl, &
        config%transport_cfl, &
        config%chemistry_relative_tolerance, &
        config%chemistry_absolute_tolerance, &
        config%initial_temperature, config%initial_pressure, &
        config%initial_velocity_x, config%initial_velocity_y, &
        config%initial_velocity_z, config%density_wave_amplitude, &
        config%hotspot_temperature_rise, config%hotspot_center_x, &
        config%hotspot_center_y, config%hotspot_center_z, &
        config%hotspot_width]))) then
      message = "Reactive 3D scalar inputs must be finite"
    else if (config%final_time <= 0.0_dp) then
      message = "Reactive 3D final_time must be positive"
    else if (config%cfl <= 0.0_dp .or. config%cfl > 1.0_dp) then
      message = "Reactive 3D cfl must be in (0, 1]"
    else if (config%transport_cfl <= 0.0_dp .or. &
             config%transport_cfl > 0.5_dp) then
      message = "Reactive 3D transport_cfl must be in (0, 0.5]"
    else if (config%transport_enabled .and. &
             .not. (config%viscosity_enabled .or. &
               config%thermal_conduction_enabled .or. &
               config%species_diffusion_enabled)) then
      message = "Reactive 3D transport requires an enabled process"
    else if (config%transport_enabled .and. &
             config%barodiffusion_enabled .and. &
             .not. config%species_diffusion_enabled) then
      message = "Reactive 3D barodiffusion requires species diffusion"
    else if (config%chemistry_relative_tolerance <= 0.0_dp .or. &
             config%chemistry_absolute_tolerance <= 0.0_dp) then
      message = "Reactive 3D chemistry tolerances must be positive"
    else if (config%initial_temperature <= 0.0_dp .or. &
             config%initial_pressure <= 0.0_dp) then
      message = "Reactive 3D initial temperature and pressure must be positive"
    else if (trim(config%problem) == "entropy_wave" .and. &
             abs(config%density_wave_amplitude) >= 1.0_dp) then
      message = "Reactive 3D relative density amplitude must be below one"
    else if (trim(config%problem) == "entropy_wave" .and. &
             config%wave_number_x == 0 .and. &
             config%wave_number_y == 0 .and. &
             config%wave_number_z == 0) then
      message = "Reactive 3D requires at least one nonzero wave number"
    else if (trim(config%problem) /= "entropy_wave" .and. &
             trim(config%problem) /= "uniform_reactor" .and. &
             trim(config%problem) /= "reactive_hotspot") then
      message = "Unknown reactive 3D problem"
    else if (trim(config%problem) == "reactive_hotspot" .and. &
             (config%hotspot_temperature_rise < 0.0_dp .or. &
              config%hotspot_width <= 0.0_dp .or. &
              config%hotspot_center_x < config%x_lower .or. &
              config%hotspot_center_x > config%x_upper .or. &
              config%hotspot_center_y < config%y_lower .or. &
              config%hotspot_center_y > config%y_upper .or. &
              config%hotspot_center_z < config%z_lower .or. &
              config%hotspot_center_z > config%z_upper)) then
      message = "Invalid reactive 3D hotspot geometry"
    else if (trim(config%thermo_model) /= "elementary" .and. &
             trim(config%thermo_model) /= "full_h2o2" .and. &
             .not. (selected_model .and. selected_allowed)) then
      message = &
        "Reactive 3D thermo_model must be elementary, full_h2o2, or selected"
    else if (trim(config%reconstruction) /= "pcm" .and. &
             trim(config%reconstruction) /= "characteristic_plm") then
      message = "Reactive 3D reconstruction must be pcm or characteristic_plm"
    else if (trim(config%reconstruction) == "characteristic_plm" .and. &
             min(config%nx, config%ny, config%nz) < 3) then
      message = "Reactive 3D characteristic_plm requires nx, ny, nz >= 3"
    else if (trim(config%limiter) /= "minmod" .and. &
             trim(config%limiter) /= "mc") then
      message = "Reactive 3D limiter must be minmod or mc"
    else if (trim(config%riemann_solver) /= "rusanov" .and. &
             trim(config%riemann_solver) /= "hllc" .and. &
             trim(config%riemann_solver) /= "pelec") then
      message = "Reactive 3D Riemann solver must be rusanov, hllc, or pelec"
    else if (trim(config%boundary_condition) /= "periodic") then
      message = "Reactive 3D currently requires periodic boundaries"
    else if (len_trim(config%output_file) == 0) then
      message = "Reactive 3D output_file must not be empty"
    else
      if (selected_model) then
        call validate_selected_composition_fields( &
          "Reactive 3D", config%composition_count, &
          config%composition_species, config%composition_mole_fractions, &
          ok, message)
        return
      end if
      if (.not. all(ieee_is_finite([ &
          config%x_h2, config%x_h, config%x_o, config%x_o2, config%x_oh, &
          config%x_h2o, config%x_ho2, config%x_h2o2, config%x_ar, &
          config%x_n2]))) then
        message = "Reactive 3D mole fractions must be finite"
        return
      end if
      if (min(config%x_h2, config%x_h, config%x_o, config%x_o2, &
              config%x_oh, config%x_h2o, config%x_ho2, &
              config%x_h2o2, config%x_ar, config%x_n2) < 0.0_dp) then
        message = "Reactive 3D mole fractions must be nonnegative"
        return
      end if
      mole_sum = config%x_h2 + config%x_h + config%x_o + config%x_o2 + &
        config%x_oh + config%x_h2o + config%x_ho2 + config%x_h2o2 + &
        config%x_ar + config%x_n2
      if (abs(mole_sum - 1.0_dp) > 5.0e-10_dp) then
        message = "Reactive 3D mole fractions must sum to one"
      else if (trim(config%thermo_model) == "elementary" .and. &
               max(config%x_ho2, config%x_h2o2, config%x_ar) > &
                 5.0e-14_dp) then
        message = "Elementary thermo requires zero HO2/H2O2/AR fractions"
      else
        ok = .true.
      end if
    end if
  end subroutine validate_reactive_3d_configuration

  subroutine reactive_3d_mole_fractions(config, nspecies, mole_fractions, ok)
    type(reactive_3d_config), intent(in) :: config
    integer, intent(in) :: nspecies
    real(dp), intent(out) :: mole_fractions(:)
    logical, intent(out) :: ok

    mole_fractions = 0.0_dp
    ok = size(mole_fractions) == nspecies
    if (.not. ok) return
    select case (trim(config%thermo_model))
    case ("elementary")
      if (nspecies /= 7) then
        ok = .false.
        return
      end if
      mole_fractions = [ &
        config%x_h2, config%x_h, config%x_o, config%x_o2, &
        config%x_oh, config%x_h2o, config%x_n2]
    case ("full_h2o2")
      if (nspecies /= 10) then
        ok = .false.
        return
      end if
      mole_fractions = [ &
        config%x_h2, config%x_h, config%x_o, config%x_o2, &
        config%x_oh, config%x_h2o, config%x_ho2, config%x_h2o2, &
        config%x_ar, config%x_n2]
    case default
      ok = .false.
      return
    end select
    ok = minval(mole_fractions) >= 0.0_dp .and. &
      abs(sum(mole_fractions) - 1.0_dp) <= 5.0e-10_dp
  end subroutine reactive_3d_mole_fractions


  subroutine resolve_reactive_3d_selected_composition( &
      config, species, mole_fractions, ok, message)
    type(reactive_3d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(out) :: mole_fractions(:)
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    mole_fractions = 0.0_dp
    if (trim(config%thermo_model) /= "selected") then
      ok = .false.
      message = &
        "Reactive 3D selected composition requires thermo_model='selected'"
      return
    end if
    call resolve_selected_composition( &
      "Reactive 3D", config%composition_count, config%composition_species, &
      config%composition_mole_fractions, species, mole_fractions, ok, message)
  end subroutine resolve_reactive_3d_selected_composition

end module simulation_config_reactive_3d_mod
