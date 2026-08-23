module simulation_config_reactive_1d_mod
  use precision_mod, only: dp
  implicit none
  private

  type, public :: reactive_1d_config
    integer :: nx = 96
    integer :: maximum_steps = 200000
    real(dp) :: x_lower = 0.0_dp
    real(dp) :: x_upper = 0.012_dp
    real(dp) :: final_time = 4.0e-5_dp
    real(dp) :: cfl = 0.35_dp
    character(len=32) :: problem = "reactive_hotspot"
    character(len=32) :: reconstruction = "characteristic_plm"
    character(len=32) :: riemann_solver = "rusanov"
    character(len=32) :: limiter = "mc"
    character(len=32) :: boundary_condition = "periodic"
    logical :: chemistry_enabled = .true.
    character(len=32) :: chemistry_model = "elementary"
    logical :: transport_enabled = .false.
    logical :: viscosity_enabled = .true.
    logical :: thermal_conduction_enabled = .true.
    logical :: species_diffusion_enabled = .true.
    logical :: barodiffusion_enabled = .true.
    real(dp) :: transport_cfl = 0.40_dp
    logical :: ppm_contact_steepening = .false.
    logical :: ppm_shock_flattening = .false.
    logical :: amr_enabled = .false.
    character(len=32) :: amr_reconstruction = "pcm"
    integer :: amr_refinement_ratio = 2
    integer :: amr_regrid_interval = 4
    integer :: amr_tag_component = 1
    integer :: amr_buffer_cells = 2
    integer :: amr_minimum_patch_cells = 8
    real(dp) :: amr_relative_gradient_threshold = 0.02_dp
    real(dp) :: amr_absolute_gradient_threshold = 0.0_dp
    real(dp) :: amr_scale_floor = 1.0e-12_dp
    real(dp) :: chemistry_relative_tolerance = 2.0e-7_dp
    real(dp) :: chemistry_absolute_tolerance = 1.0e-12_dp
    real(dp) :: initial_temperature = 1200.0_dp
    real(dp) :: initial_pressure = 101325.0_dp
    real(dp) :: initial_velocity = 0.0_dp
    real(dp) :: density_wave_amplitude = 0.1_dp
    real(dp) :: composition_wave_amplitude = 0.04_dp
    real(dp) :: hotspot_temperature_rise = 250.0_dp
    real(dp) :: hotspot_center = 0.006_dp
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
    character(len=256) :: output_file = "reactive_1d.csv"
  end type reactive_1d_config

  public :: read_reactive_1d_configuration
  public :: reactive_1d_mole_fractions

contains

  subroutine read_reactive_1d_configuration(path, config, ok, message)
    character(len=*), intent(in) :: path
    type(reactive_1d_config), intent(out) :: config
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    integer :: nx, maximum_steps, unit, status
    integer :: amr_refinement_ratio, amr_regrid_interval
    integer :: amr_tag_component, amr_buffer_cells
    integer :: amr_minimum_patch_cells
    real(dp) :: x_lower, x_upper, final_time, cfl
    real(dp) :: chemistry_relative_tolerance, chemistry_absolute_tolerance
    real(dp) :: transport_cfl
    real(dp) :: amr_relative_gradient_threshold
    real(dp) :: amr_absolute_gradient_threshold, amr_scale_floor
    real(dp) :: initial_temperature, initial_pressure, initial_velocity
    real(dp) :: density_wave_amplitude, composition_wave_amplitude
    real(dp) :: hotspot_temperature_rise
    real(dp) :: hotspot_center, hotspot_width
    real(dp) :: x_h2, x_h, x_o, x_o2, x_oh, x_h2o, x_ho2, x_h2o2, x_ar, x_n2
    character(len=32) :: problem, reconstruction, riemann_solver, limiter
    character(len=32) :: chemistry_model, amr_reconstruction
    character(len=32) :: boundary_condition
    character(len=256) :: output_file
    logical :: chemistry_enabled
    logical :: transport_enabled, viscosity_enabled
    logical :: thermal_conduction_enabled, species_diffusion_enabled
    logical :: barodiffusion_enabled
    logical :: ppm_contact_steepening, ppm_shock_flattening
    logical :: amr_enabled
    real(dp) :: mole_sum
    namelist /reactive_1d/ &
      nx, maximum_steps, x_lower, x_upper, final_time, cfl, problem, &
      reconstruction, riemann_solver, limiter, boundary_condition, &
      chemistry_enabled, chemistry_model, transport_enabled, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, transport_cfl, ppm_contact_steepening, &
      ppm_shock_flattening, chemistry_relative_tolerance, &
      chemistry_absolute_tolerance, &
      amr_enabled, amr_reconstruction, &
      amr_refinement_ratio, amr_regrid_interval, &
      amr_tag_component, amr_buffer_cells, amr_minimum_patch_cells, &
      amr_relative_gradient_threshold, amr_absolute_gradient_threshold, &
      amr_scale_floor, &
      initial_temperature, initial_pressure, initial_velocity, &
      density_wave_amplitude, composition_wave_amplitude, &
      hotspot_temperature_rise, hotspot_center, hotspot_width, x_h2, x_h, &
      x_o, x_o2, x_oh, x_h2o, x_ho2, x_h2o2, x_ar, x_n2, output_file

    config = reactive_1d_config()
    nx = config%nx
    maximum_steps = config%maximum_steps
    x_lower = config%x_lower
    x_upper = config%x_upper
    final_time = config%final_time
    cfl = config%cfl
    problem = config%problem
    reconstruction = config%reconstruction
    riemann_solver = config%riemann_solver
    limiter = config%limiter
    boundary_condition = config%boundary_condition
    chemistry_enabled = config%chemistry_enabled
    chemistry_model = config%chemistry_model
    transport_enabled = config%transport_enabled
    viscosity_enabled = config%viscosity_enabled
    thermal_conduction_enabled = config%thermal_conduction_enabled
    species_diffusion_enabled = config%species_diffusion_enabled
    barodiffusion_enabled = config%barodiffusion_enabled
    transport_cfl = config%transport_cfl
    ppm_contact_steepening = config%ppm_contact_steepening
    ppm_shock_flattening = config%ppm_shock_flattening
    amr_enabled = config%amr_enabled
    amr_reconstruction = config%amr_reconstruction
    amr_refinement_ratio = config%amr_refinement_ratio
    amr_regrid_interval = config%amr_regrid_interval
    amr_tag_component = config%amr_tag_component
    amr_buffer_cells = config%amr_buffer_cells
    amr_minimum_patch_cells = config%amr_minimum_patch_cells
    amr_relative_gradient_threshold = &
      config%amr_relative_gradient_threshold
    amr_absolute_gradient_threshold = &
      config%amr_absolute_gradient_threshold
    amr_scale_floor = config%amr_scale_floor
    chemistry_relative_tolerance = config%chemistry_relative_tolerance
    chemistry_absolute_tolerance = config%chemistry_absolute_tolerance
    initial_temperature = config%initial_temperature
    initial_pressure = config%initial_pressure
    initial_velocity = config%initial_velocity
    density_wave_amplitude = config%density_wave_amplitude
    composition_wave_amplitude = config%composition_wave_amplitude
    hotspot_temperature_rise = config%hotspot_temperature_rise
    hotspot_center = config%hotspot_center
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
      message = "Could not open reactive 1D input"
      return
    end if
    read(unit, nml=reactive_1d, iostat=status)
    close(unit)
    if (status /= 0) then
      ok = .false.
      message = "Could not parse reactive 1D namelist"
      return
    end if

    mole_sum = x_h2 + x_h + x_o + x_o2 + x_oh + x_h2o + x_ho2 + x_h2o2 + x_ar + x_n2
    ok = nx >= 8 .and. maximum_steps >= 1 .and. x_upper > x_lower .and. &
      final_time > 0.0_dp .and. cfl > 0.0_dp .and. cfl <= 0.9_dp .and. &
      initial_temperature > 0.0_dp .and. initial_pressure > 0.0_dp .and. &
      density_wave_amplitude >= 0.0_dp .and. &
      density_wave_amplitude < 0.9_dp .and. hotspot_width > 0.0_dp .and. &
      composition_wave_amplitude >= 0.0_dp .and. &
      chemistry_relative_tolerance > 0.0_dp .and. &
      transport_cfl > 0.0_dp .and. transport_cfl <= 0.5_dp .and. &
      chemistry_absolute_tolerance > 0.0_dp .and. &
      min(x_h2, x_h, x_o, x_o2, x_oh, x_h2o, x_ho2, x_h2o2, x_ar, x_n2) >= 0.0_dp .and. &
      abs(mole_sum - 1.0_dp) <= 5.0e-10_dp
    if (ok .and. amr_enabled) then
      ok = amr_refinement_ratio >= 2 .and. amr_regrid_interval >= 1 .and. &
        amr_tag_component >= 1 .and. amr_buffer_cells >= 0 .and. &
        amr_minimum_patch_cells >= 1 .and. &
        amr_minimum_patch_cells <= nx - 2 .and. &
        amr_relative_gradient_threshold >= 0.0_dp .and. &
        amr_absolute_gradient_threshold >= 0.0_dp .and. &
        amr_scale_floor > 0.0_dp
      ok = ok .and. (trim(amr_reconstruction) == "pcm" .or. &
        trim(amr_reconstruction) == "plm")
    end if
    if (.not. ok) then
      message = "Invalid reactive 1D configuration"
      return
    end if
    if (trim(chemistry_model) /= "elementary" .and. &
        trim(chemistry_model) /= "full_h2o2") then
      ok = .false.
      message = "Unknown reactive 1D chemistry model"
      return
    end if
    if (trim(chemistry_model) == "elementary" .and. &
        max(x_ho2, x_h2o2, x_ar) > 5.0e-14_dp) then
      ok = .false.
      message = "Elementary chemistry requires zero HO2/H2O2/AR mole fractions"
      return
    end if
    if (trim(problem) /= "entropy_wave" .and. &
        trim(problem) /= "composition_wave" .and. &
        trim(problem) /= "reactive_hotspot" .and. &
        trim(problem) /= "uniform_reactor") then
      ok = .false.
      message = "Unknown reactive 1D problem"
      return
    end if
    if (trim(problem) == "composition_wave" .and. &
        composition_wave_amplitude > min(x_h2, x_n2)) then
      ok = .false.
      message = "Composition-wave amplitude exceeds the H2/N2 base fraction"
      return
    end if
    if (trim(reconstruction) /= "pcm" .and. &
        trim(reconstruction) /= "characteristic_plm" .and. &
        trim(reconstruction) /= "ppm" .and. &
        trim(reconstruction) /= "characteristic_ppm") then
      ok = .false.
      message = "Unknown reactive 1D reconstruction"
      return
    end if
    if ((ppm_contact_steepening .or. ppm_shock_flattening) .and. &
        trim(reconstruction) /= "characteristic_ppm") then
      ok = .false.
      message = "Reactive PPM steepening/flattening requires characteristic_ppm"
      return
    end if
    if (trim(riemann_solver) /= "rusanov" .and. &
        trim(riemann_solver) /= "hllc") then
      ok = .false.
      message = "Unknown reactive 1D Riemann solver"
      return
    end if
    if (trim(limiter) /= "minmod" .and. trim(limiter) /= "mc") then
      ok = .false.
      message = "Unknown reactive 1D limiter"
      return
    end if
    if (trim(boundary_condition) /= "periodic" .and. &
        trim(boundary_condition) /= "outflow") then
      ok = .false.
      message = "Unknown reactive 1D boundary condition"
      return
    end if

    config%nx = nx
    config%maximum_steps = maximum_steps
    config%x_lower = x_lower
    config%x_upper = x_upper
    config%final_time = final_time
    config%cfl = cfl
    config%problem = trim(problem)
    config%reconstruction = trim(reconstruction)
    config%riemann_solver = trim(riemann_solver)
    config%limiter = trim(limiter)
    config%boundary_condition = trim(boundary_condition)
    config%chemistry_enabled = chemistry_enabled
    config%chemistry_model = trim(chemistry_model)
    config%transport_enabled = transport_enabled
    config%viscosity_enabled = viscosity_enabled
    config%thermal_conduction_enabled = thermal_conduction_enabled
    config%species_diffusion_enabled = species_diffusion_enabled
    config%barodiffusion_enabled = barodiffusion_enabled
    config%transport_cfl = transport_cfl
    config%ppm_contact_steepening = ppm_contact_steepening
    config%ppm_shock_flattening = ppm_shock_flattening
    config%amr_enabled = amr_enabled
    config%amr_reconstruction = trim(amr_reconstruction)
    config%amr_refinement_ratio = amr_refinement_ratio
    config%amr_regrid_interval = amr_regrid_interval
    config%amr_tag_component = amr_tag_component
    config%amr_buffer_cells = amr_buffer_cells
    config%amr_minimum_patch_cells = amr_minimum_patch_cells
    config%amr_relative_gradient_threshold = &
      amr_relative_gradient_threshold
    config%amr_absolute_gradient_threshold = &
      amr_absolute_gradient_threshold
    config%amr_scale_floor = amr_scale_floor
    config%chemistry_relative_tolerance = chemistry_relative_tolerance
    config%chemistry_absolute_tolerance = chemistry_absolute_tolerance
    config%initial_temperature = initial_temperature
    config%initial_pressure = initial_pressure
    config%initial_velocity = initial_velocity
    config%density_wave_amplitude = density_wave_amplitude
    config%composition_wave_amplitude = composition_wave_amplitude
    config%hotspot_temperature_rise = hotspot_temperature_rise
    config%hotspot_center = hotspot_center
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
  end subroutine read_reactive_1d_configuration


  subroutine reactive_1d_mole_fractions(config, nspecies, mole_fractions, ok)
    type(reactive_1d_config), intent(in) :: config
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
  end subroutine reactive_1d_mole_fractions

end module simulation_config_reactive_1d_mod
