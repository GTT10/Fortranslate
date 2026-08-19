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
    real(dp) :: x_n2 = 0.55643_dp
    character(len=256) :: output_file = "reactive_1d.csv"
  end type reactive_1d_config

  public :: read_reactive_1d_configuration

contains

  subroutine read_reactive_1d_configuration(path, config, ok, message)
    character(len=*), intent(in) :: path
    type(reactive_1d_config), intent(out) :: config
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    integer :: nx, maximum_steps, unit, status
    real(dp) :: x_lower, x_upper, final_time, cfl
    real(dp) :: chemistry_relative_tolerance, chemistry_absolute_tolerance
    real(dp) :: initial_temperature, initial_pressure, initial_velocity
    real(dp) :: density_wave_amplitude, composition_wave_amplitude
    real(dp) :: hotspot_temperature_rise
    real(dp) :: hotspot_center, hotspot_width
    real(dp) :: x_h2, x_h, x_o, x_o2, x_oh, x_h2o, x_n2
    character(len=32) :: problem, reconstruction, riemann_solver, limiter
    character(len=32) :: boundary_condition
    character(len=256) :: output_file
    logical :: chemistry_enabled
    real(dp) :: mole_sum
    namelist /reactive_1d/ &
      nx, maximum_steps, x_lower, x_upper, final_time, cfl, problem, &
      reconstruction, riemann_solver, limiter, boundary_condition, &
      chemistry_enabled, &
      chemistry_relative_tolerance, chemistry_absolute_tolerance, &
      initial_temperature, initial_pressure, initial_velocity, &
      density_wave_amplitude, composition_wave_amplitude, &
      hotspot_temperature_rise, hotspot_center, hotspot_width, x_h2, x_h, &
      x_o, x_o2, x_oh, x_h2o, x_n2, output_file

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

    mole_sum = x_h2 + x_h + x_o + x_o2 + x_oh + x_h2o + x_n2
    ok = nx >= 8 .and. maximum_steps >= 1 .and. x_upper > x_lower .and. &
      final_time > 0.0_dp .and. cfl > 0.0_dp .and. cfl <= 0.9_dp .and. &
      initial_temperature > 0.0_dp .and. initial_pressure > 0.0_dp .and. &
      density_wave_amplitude >= 0.0_dp .and. &
      density_wave_amplitude < 0.9_dp .and. hotspot_width > 0.0_dp .and. &
      composition_wave_amplitude >= 0.0_dp .and. &
      chemistry_relative_tolerance > 0.0_dp .and. &
      chemistry_absolute_tolerance > 0.0_dp .and. &
      min(x_h2, x_h, x_o, x_o2, x_oh, x_h2o, x_n2) >= 0.0_dp .and. &
      abs(mole_sum - 1.0_dp) <= 5.0e-10_dp
    if (.not. ok) then
      message = "Invalid reactive 1D configuration"
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
        trim(reconstruction) /= "ppm") then
      ok = .false.
      message = "Unknown reactive 1D reconstruction"
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
    config%x_n2 = x_n2
    config%output_file = trim(output_file)
  end subroutine read_reactive_1d_configuration

end module simulation_config_reactive_1d_mod
