module simulation_config_reactive_2d_mod
  use precision_mod, only: dp
  implicit none
  private

  type, public :: reactive_2d_config
    integer :: nx = 32
    integer :: ny = 32
    integer :: maximum_steps = 200000
    real(dp) :: x_lower = 0.0_dp
    real(dp) :: x_upper = 0.01_dp
    real(dp) :: y_lower = 0.0_dp
    real(dp) :: y_upper = 0.01_dp
    real(dp) :: final_time = 5.0e-6_dp
    real(dp) :: cfl = 0.30_dp
    character(len=32) :: problem = "diagonal_wave"
    character(len=32) :: reconstruction = "characteristic_plm"
    character(len=32) :: riemann_solver = "hllc"
    character(len=32) :: limiter = "mc"
    logical :: use_transverse_correction = .true.
    logical :: chemistry_enabled = .false.
    character(len=32) :: chemistry_model = "elementary"
    logical :: transport_enabled = .false.
    logical :: viscosity_enabled = .true.
    logical :: thermal_conduction_enabled = .true.
    logical :: species_diffusion_enabled = .true.
    logical :: barodiffusion_enabled = .true.
    real(dp) :: transport_cfl = 0.35_dp
    character(len=24) :: boundary_x_lower = "periodic"
    character(len=24) :: boundary_x_upper = "periodic"
    character(len=24) :: boundary_y_lower = "periodic"
    character(len=24) :: boundary_y_upper = "periodic"
    character(len=24) :: thermal_x_lower = "adiabatic"
    character(len=24) :: thermal_x_upper = "adiabatic"
    character(len=24) :: thermal_y_lower = "adiabatic"
    character(len=24) :: thermal_y_upper = "adiabatic"
    real(dp) :: wall_temperature_x_lower = 300.0_dp
    real(dp) :: wall_temperature_x_upper = 300.0_dp
    real(dp) :: wall_temperature_y_lower = 300.0_dp
    real(dp) :: wall_temperature_y_upper = 300.0_dp
    real(dp) :: wall_velocity_x_lower(3) = 0.0_dp
    real(dp) :: wall_velocity_x_upper(3) = 0.0_dp
    real(dp) :: wall_velocity_y_lower(3) = 0.0_dp
    real(dp) :: wall_velocity_y_upper(3) = 0.0_dp
    logical :: ppm_contact_steepening = .false.
    logical :: ppm_shock_flattening = .false.
    real(dp) :: chemistry_relative_tolerance = 2.0e-7_dp
    real(dp) :: chemistry_absolute_tolerance = 1.0e-12_dp
    real(dp) :: initial_temperature = 1000.0_dp
    real(dp) :: initial_pressure = 101325.0_dp
    real(dp) :: initial_velocity_x = 300.0_dp
    real(dp) :: initial_velocity_y = 200.0_dp
    real(dp) :: density_wave_amplitude = 0.08_dp
    real(dp) :: composition_wave_amplitude = 0.04_dp
    real(dp) :: vortex_strength = 45.0_dp
    real(dp) :: vortex_center_x = 0.005_dp
    real(dp) :: vortex_center_y = 0.005_dp
    real(dp) :: vortex_radius = 0.0015_dp
    real(dp) :: hotspot_temperature_rise = 250.0_dp
    real(dp) :: hotspot_center_x = 0.005_dp
    real(dp) :: hotspot_center_y = 0.005_dp
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
    character(len=256) :: output_file = "reactive_2d.csv"
  end type reactive_2d_config

  public :: read_reactive_2d_configuration
  public :: reactive_2d_mole_fractions

contains

  pure logical function valid_boundary_kind(kind) result(valid)
    character(len=*), intent(in) :: kind
    valid = trim(kind) == "periodic" .or. trim(kind) == "slip_wall" .or. &
      trim(kind) == "no_slip_wall" .or. trim(kind) == "inflow" .or. &
      trim(kind) == "outflow"
  end function valid_boundary_kind

  pure logical function valid_boundary_pair(lower, upper) result(valid)
    character(len=*), intent(in) :: lower, upper
    valid = valid_boundary_kind(lower) .and. valid_boundary_kind(upper)
    if (.not. valid) return
    valid = (trim(lower) == "periodic") .eqv. (trim(upper) == "periodic")
  end function valid_boundary_pair

  pure logical function valid_thermal_boundary(kind) result(valid)
    character(len=*), intent(in) :: kind
    valid = trim(kind) == "adiabatic" .or. trim(kind) == "isothermal"
  end function valid_thermal_boundary

  subroutine read_reactive_2d_configuration(path, config, ok, message)
    character(len=*), intent(in) :: path
    type(reactive_2d_config), intent(out) :: config
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    integer :: nx, ny, maximum_steps, unit, status
    real(dp) :: x_lower, x_upper, y_lower, y_upper, final_time, cfl
    real(dp) :: chemistry_relative_tolerance, chemistry_absolute_tolerance
    real(dp) :: transport_cfl
    real(dp) :: wall_temperature_x_lower, wall_temperature_x_upper
    real(dp) :: wall_temperature_y_lower, wall_temperature_y_upper
    real(dp) :: wall_velocity_x_lower(3), wall_velocity_x_upper(3)
    real(dp) :: wall_velocity_y_lower(3), wall_velocity_y_upper(3)
    real(dp) :: initial_temperature, initial_pressure
    real(dp) :: initial_velocity_x, initial_velocity_y
    real(dp) :: density_wave_amplitude, composition_wave_amplitude
    real(dp) :: vortex_strength
    real(dp) :: vortex_center_x, vortex_center_y, vortex_radius
    real(dp) :: hotspot_temperature_rise, hotspot_center_x
    real(dp) :: hotspot_center_y, hotspot_width
    real(dp) :: x_h2, x_h, x_o, x_o2, x_oh, x_h2o, x_ho2, x_h2o2, x_ar, x_n2, mole_sum
    character(len=32) :: problem, reconstruction, riemann_solver, limiter
    character(len=32) :: chemistry_model
    character(len=24) :: boundary_x_lower, boundary_x_upper
    character(len=24) :: boundary_y_lower, boundary_y_upper
    character(len=24) :: thermal_x_lower, thermal_x_upper
    character(len=24) :: thermal_y_lower, thermal_y_upper
    character(len=256) :: output_file
    logical :: use_transverse_correction, chemistry_enabled
    logical :: transport_enabled, viscosity_enabled
    logical :: thermal_conduction_enabled, species_diffusion_enabled
    logical :: barodiffusion_enabled
    logical :: ppm_contact_steepening, ppm_shock_flattening
    namelist /reactive_2d/ &
      nx, ny, maximum_steps, x_lower, x_upper, y_lower, y_upper, &
      final_time, cfl, problem, reconstruction, riemann_solver, limiter, &
      use_transverse_correction, chemistry_enabled, chemistry_model, transport_enabled, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, transport_cfl, &
      boundary_x_lower, boundary_x_upper, boundary_y_lower, boundary_y_upper, &
      thermal_x_lower, thermal_x_upper, thermal_y_lower, thermal_y_upper, &
      wall_temperature_x_lower, wall_temperature_x_upper, &
      wall_temperature_y_lower, wall_temperature_y_upper, &
      wall_velocity_x_lower, wall_velocity_x_upper, &
      wall_velocity_y_lower, wall_velocity_y_upper, &
      ppm_contact_steepening, ppm_shock_flattening, &
      chemistry_relative_tolerance, chemistry_absolute_tolerance, &
      initial_temperature, initial_pressure, initial_velocity_x, &
      initial_velocity_y, density_wave_amplitude, &
      composition_wave_amplitude, vortex_strength, &
      vortex_center_x, vortex_center_y, vortex_radius, &
      hotspot_temperature_rise, hotspot_center_x, hotspot_center_y, &
      hotspot_width, x_h2, x_h, x_o, x_o2, x_oh, x_h2o, x_ho2, x_h2o2, &
      x_ar, x_n2, output_file

    config = reactive_2d_config()
    nx = config%nx
    ny = config%ny
    maximum_steps = config%maximum_steps
    x_lower = config%x_lower
    x_upper = config%x_upper
    y_lower = config%y_lower
    y_upper = config%y_upper
    final_time = config%final_time
    cfl = config%cfl
    problem = config%problem
    reconstruction = config%reconstruction
    riemann_solver = config%riemann_solver
    limiter = config%limiter
    use_transverse_correction = config%use_transverse_correction
    chemistry_enabled = config%chemistry_enabled
    chemistry_model = config%chemistry_model
    transport_enabled = config%transport_enabled
    viscosity_enabled = config%viscosity_enabled
    thermal_conduction_enabled = config%thermal_conduction_enabled
    species_diffusion_enabled = config%species_diffusion_enabled
    barodiffusion_enabled = config%barodiffusion_enabled
    transport_cfl = config%transport_cfl
    boundary_x_lower = config%boundary_x_lower
    boundary_x_upper = config%boundary_x_upper
    boundary_y_lower = config%boundary_y_lower
    boundary_y_upper = config%boundary_y_upper
    thermal_x_lower = config%thermal_x_lower
    thermal_x_upper = config%thermal_x_upper
    thermal_y_lower = config%thermal_y_lower
    thermal_y_upper = config%thermal_y_upper
    wall_temperature_x_lower = config%wall_temperature_x_lower
    wall_temperature_x_upper = config%wall_temperature_x_upper
    wall_temperature_y_lower = config%wall_temperature_y_lower
    wall_temperature_y_upper = config%wall_temperature_y_upper
    wall_velocity_x_lower = config%wall_velocity_x_lower
    wall_velocity_x_upper = config%wall_velocity_x_upper
    wall_velocity_y_lower = config%wall_velocity_y_lower
    wall_velocity_y_upper = config%wall_velocity_y_upper
    ppm_contact_steepening = config%ppm_contact_steepening
    ppm_shock_flattening = config%ppm_shock_flattening
    chemistry_relative_tolerance = config%chemistry_relative_tolerance
    chemistry_absolute_tolerance = config%chemistry_absolute_tolerance
    initial_temperature = config%initial_temperature
    initial_pressure = config%initial_pressure
    initial_velocity_x = config%initial_velocity_x
    initial_velocity_y = config%initial_velocity_y
    density_wave_amplitude = config%density_wave_amplitude
    composition_wave_amplitude = config%composition_wave_amplitude
    vortex_strength = config%vortex_strength
    vortex_center_x = config%vortex_center_x
    vortex_center_y = config%vortex_center_y
    vortex_radius = config%vortex_radius
    hotspot_temperature_rise = config%hotspot_temperature_rise
    hotspot_center_x = config%hotspot_center_x
    hotspot_center_y = config%hotspot_center_y
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
    output_file = config%output_file

    message = ""
    open(newunit=unit, file=trim(path), status="old", action="read", &
      iostat=status)
    if (status /= 0) then
      ok = .false.
      message = "Could not open reactive 2D input"
      return
    end if
    read(unit, nml=reactive_2d, iostat=status)
    close(unit)
    if (status /= 0) then
      ok = .false.
      message = "Could not parse reactive 2D namelist"
      return
    end if

    mole_sum = x_h2 + x_h + x_o + x_o2 + x_oh + x_h2o + x_ho2 + x_h2o2 + x_ar + x_n2
    ok = nx >= 4 .and. ny >= 4 .and. maximum_steps >= 1 .and. &
      x_upper > x_lower .and. y_upper > y_lower .and. final_time > 0.0_dp .and. &
      cfl > 0.0_dp .and. cfl <= 0.8_dp .and. initial_temperature > 0.0_dp .and. &
      initial_pressure > 0.0_dp .and. density_wave_amplitude >= 0.0_dp .and. &
      density_wave_amplitude < 0.9_dp .and. vortex_radius > 0.0_dp .and. &
      composition_wave_amplitude >= 0.0_dp .and. &
      hotspot_width > 0.0_dp .and. transport_cfl > 0.0_dp .and. &
      min(wall_temperature_x_lower, wall_temperature_x_upper, &
        wall_temperature_y_lower, wall_temperature_y_upper) > 0.0_dp .and. &
      transport_cfl <= 0.5_dp .and. chemistry_relative_tolerance > 0.0_dp .and. &
      chemistry_absolute_tolerance > 0.0_dp .and. &
      min(x_h2, x_h, x_o, x_o2, x_oh, x_h2o, x_ho2, x_h2o2, x_ar, x_n2) >= 0.0_dp .and. &
      abs(mole_sum - 1.0_dp) <= 5.0e-10_dp
    if (.not. ok) then
      message = "Invalid reactive 2D configuration"
      return
    end if
    if (trim(chemistry_model) /= "elementary" .and. &
        trim(chemistry_model) /= "full_h2o2") then
      ok = .false.
      message = "Unknown reactive 2D chemistry model"
      return
    end if
    if (trim(chemistry_model) == "elementary" .and. &
        max(x_ho2, x_h2o2, x_ar) > 5.0e-14_dp) then
      ok = .false.
      message = "Elementary chemistry requires zero HO2/H2O2/AR mole fractions"
      return
    end if
    if (transport_enabled .and. .not. (viscosity_enabled .or. &
        thermal_conduction_enabled .or. species_diffusion_enabled)) then
      ok = .false.
      message = "Reactive 2D transport requires an enabled process"
      return
    end if
    if (trim(problem) /= "diagonal_wave" .and. &
        trim(problem) /= "diagonal_composition_wave" .and. &
        trim(problem) /= "reactive_vortex" .and. &
        trim(problem) /= "reactive_hotspot" .and. &
        trim(problem) /= "uniform_reactor" .and. &
        trim(problem) /= "couette_channel" .and. &
        trim(problem) /= "thermal_channel" .and. &
        trim(problem) /= "inflow_outflow") then
      ok = .false.
      message = "Unknown reactive 2D problem"
      return
    end if
    if (trim(problem) == "diagonal_composition_wave" .and. &
        composition_wave_amplitude > min(x_h2, x_n2)) then
      ok = .false.
      message = "Reactive 2D composition-wave amplitude exceeds H2/N2 base fraction"
      return
    end if
    if (.not. valid_boundary_pair(boundary_x_lower, boundary_x_upper) .or. &
        .not. valid_boundary_pair(boundary_y_lower, boundary_y_upper)) then
      ok = .false.
      message = "Invalid or unmatched reactive 2D boundary pair"
      return
    end if
    if (.not. valid_thermal_boundary(thermal_x_lower) .or. &
        .not. valid_thermal_boundary(thermal_x_upper) .or. &
        .not. valid_thermal_boundary(thermal_y_lower) .or. &
        .not. valid_thermal_boundary(thermal_y_upper)) then
      ok = .false.
      message = "Unknown reactive 2D thermal boundary"
      return
    end if
    if (trim(reconstruction) /= "pcm" .and. &
        trim(reconstruction) /= "characteristic_plm" .and. &
        trim(reconstruction) /= "characteristic_ppm") then
      ok = .false.
      message = "Reactive 2D supports pcm, characteristic_plm, or characteristic_ppm"
      return
    end if
    if ((ppm_contact_steepening .or. ppm_shock_flattening) .and. &
        trim(reconstruction) /= "characteristic_ppm") then
      ok = .false.
      message = "Reactive 2D PPM controls require characteristic_ppm"
      return
    end if
    if (trim(riemann_solver) /= "rusanov" .and. &
        trim(riemann_solver) /= "hllc" .and. &
        trim(riemann_solver) /= "pelec") then
      ok = .false.
      message = "Unknown reactive 2D Riemann solver"
      return
    end if
    if (trim(limiter) /= "minmod" .and. trim(limiter) /= "mc") then
      ok = .false.
      message = "Unknown reactive 2D limiter"
      return
    end if

    config%nx = nx
    config%ny = ny
    config%maximum_steps = maximum_steps
    config%x_lower = x_lower
    config%x_upper = x_upper
    config%y_lower = y_lower
    config%y_upper = y_upper
    config%final_time = final_time
    config%cfl = cfl
    config%problem = trim(problem)
    config%reconstruction = trim(reconstruction)
    config%riemann_solver = trim(riemann_solver)
    config%limiter = trim(limiter)
    config%use_transverse_correction = use_transverse_correction
    config%chemistry_enabled = chemistry_enabled
    config%chemistry_model = trim(chemistry_model)
    config%transport_enabled = transport_enabled
    config%viscosity_enabled = viscosity_enabled
    config%thermal_conduction_enabled = thermal_conduction_enabled
    config%species_diffusion_enabled = species_diffusion_enabled
    config%barodiffusion_enabled = barodiffusion_enabled
    config%transport_cfl = transport_cfl
    config%boundary_x_lower = trim(boundary_x_lower)
    config%boundary_x_upper = trim(boundary_x_upper)
    config%boundary_y_lower = trim(boundary_y_lower)
    config%boundary_y_upper = trim(boundary_y_upper)
    config%thermal_x_lower = trim(thermal_x_lower)
    config%thermal_x_upper = trim(thermal_x_upper)
    config%thermal_y_lower = trim(thermal_y_lower)
    config%thermal_y_upper = trim(thermal_y_upper)
    config%wall_temperature_x_lower = wall_temperature_x_lower
    config%wall_temperature_x_upper = wall_temperature_x_upper
    config%wall_temperature_y_lower = wall_temperature_y_lower
    config%wall_temperature_y_upper = wall_temperature_y_upper
    config%wall_velocity_x_lower = wall_velocity_x_lower
    config%wall_velocity_x_upper = wall_velocity_x_upper
    config%wall_velocity_y_lower = wall_velocity_y_lower
    config%wall_velocity_y_upper = wall_velocity_y_upper
    config%ppm_contact_steepening = ppm_contact_steepening
    config%ppm_shock_flattening = ppm_shock_flattening
    config%chemistry_relative_tolerance = chemistry_relative_tolerance
    config%chemistry_absolute_tolerance = chemistry_absolute_tolerance
    config%initial_temperature = initial_temperature
    config%initial_pressure = initial_pressure
    config%initial_velocity_x = initial_velocity_x
    config%initial_velocity_y = initial_velocity_y
    config%density_wave_amplitude = density_wave_amplitude
    config%composition_wave_amplitude = composition_wave_amplitude
    config%vortex_strength = vortex_strength
    config%vortex_center_x = vortex_center_x
    config%vortex_center_y = vortex_center_y
    config%vortex_radius = vortex_radius
    config%hotspot_temperature_rise = hotspot_temperature_rise
    config%hotspot_center_x = hotspot_center_x
    config%hotspot_center_y = hotspot_center_y
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
    config%output_file = trim(output_file)
  end subroutine read_reactive_2d_configuration


  subroutine reactive_2d_mole_fractions(config, nspecies, mole_fractions, ok)
    type(reactive_2d_config), intent(in) :: config
    integer, intent(in) :: nspecies
    real(dp), intent(out) :: mole_fractions(:)
    logical, intent(out) :: ok

    ok = .false.
    if (nspecies < 1 .or. size(mole_fractions) /= nspecies) return
    select case (trim(config%chemistry_model))
    case ("elementary")
      if (nspecies /= 7) return
      mole_fractions = [config%x_h2, config%x_h, config%x_o, config%x_o2, &
        config%x_oh, config%x_h2o, config%x_n2]
    case ("full_h2o2")
      if (nspecies /= 10) return
      mole_fractions = [config%x_h2, config%x_h, config%x_o, config%x_o2, &
        config%x_oh, config%x_h2o, config%x_ho2, config%x_h2o2, &
        config%x_ar, config%x_n2]
    case default
      return
    end select
    ok = minval(mole_fractions) >= 0.0_dp .and. &
      abs(sum(mole_fractions) - 1.0_dp) <= 5.0e-10_dp
  end subroutine reactive_2d_mole_fractions

end module simulation_config_reactive_2d_mod
