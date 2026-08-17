module simulation_config_2d_mod
  use precision_mod, only: dp
  use constants_mod, only: default_gamma
  implicit none
  private

  type, public :: simulation_config_2d
    integer :: nx = 64
    integer :: ny = 64
    integer :: max_steps = 100000
    real(dp) :: x_min = 0.0_dp
    real(dp) :: x_max = 10.0_dp
    real(dp) :: y_min = 0.0_dp
    real(dp) :: y_max = 10.0_dp
    real(dp) :: final_time = 1.0_dp
    real(dp) :: cfl = 0.4_dp
    real(dp) :: gamma = default_gamma
    character(len=512) :: output_file = "isentropic_vortex.csv"
    character(len=32) :: limiter = "mc"
    character(len=32) :: boundary_condition = "periodic"
    character(len=32) :: riemann_solver = "pelec"
    logical :: use_transverse_correction = .true.
  end type simulation_config_2d

  type, public :: isentropic_vortex_config
    real(dp) :: center_x = 5.0_dp
    real(dp) :: center_y = 5.0_dp
    real(dp) :: strength = 5.0_dp
    real(dp) :: base_density = 1.0_dp
    real(dp) :: base_pressure = 1.0_dp
    real(dp) :: base_velocity_x = 1.0_dp
    real(dp) :: base_velocity_y = 1.0_dp
  end type isentropic_vortex_config

  public :: read_configuration_2d
  public :: validate_configuration_2d

contains

  subroutine read_configuration_2d(path, config, vortex, ok, message)
    character(len=*), intent(in) :: path
    type(simulation_config_2d), intent(out) :: config
    type(isentropic_vortex_config), intent(out) :: vortex
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    integer :: nx, ny, max_steps, unit, io_status
    real(dp) :: x_min, x_max, y_min, y_max, final_time, cfl, gamma
    real(dp) :: center_x, center_y, strength
    real(dp) :: base_density, base_pressure
    real(dp) :: base_velocity_x, base_velocity_y
    character(len=512) :: output_file
    character(len=32) :: limiter, boundary_condition, riemann_solver
    logical :: use_transverse_correction

    namelist /simulation_2d/ nx, ny, max_steps, x_min, x_max, y_min, &
      y_max, final_time, cfl, gamma, output_file, limiter, &
      boundary_condition, riemann_solver, use_transverse_correction
    namelist /isentropic_vortex/ center_x, center_y, strength, &
      base_density, base_pressure, base_velocity_x, base_velocity_y

    config = simulation_config_2d()
    vortex = isentropic_vortex_config()

    nx = config%nx
    ny = config%ny
    max_steps = config%max_steps
    x_min = config%x_min
    x_max = config%x_max
    y_min = config%y_min
    y_max = config%y_max
    final_time = config%final_time
    cfl = config%cfl
    gamma = config%gamma
    output_file = config%output_file
    limiter = config%limiter
    boundary_condition = config%boundary_condition
    riemann_solver = config%riemann_solver
    use_transverse_correction = config%use_transverse_correction

    center_x = vortex%center_x
    center_y = vortex%center_y
    strength = vortex%strength
    base_density = vortex%base_density
    base_pressure = vortex%base_pressure
    base_velocity_x = vortex%base_velocity_x
    base_velocity_y = vortex%base_velocity_y

    open(newunit=unit, file=trim(path), status="old", action="read", &
      iostat=io_status)
    if (io_status /= 0) then
      ok = .false.
      write(message, '(a,1x,a)') "Could not open input file:", trim(path)
      return
    end if

    read(unit, nml=simulation_2d, iostat=io_status)
    if (io_status /= 0) then
      close(unit)
      ok = .false.
      write(message, '(a,1x,a)') &
        "Could not read &simulation_2d from:", trim(path)
      return
    end if

    read(unit, nml=isentropic_vortex, iostat=io_status)
    close(unit)
    if (io_status /= 0) then
      ok = .false.
      write(message, '(a,1x,a)') &
        "Could not read &isentropic_vortex from:", trim(path)
      return
    end if

    config%nx = nx
    config%ny = ny
    config%max_steps = max_steps
    config%x_min = x_min
    config%x_max = x_max
    config%y_min = y_min
    config%y_max = y_max
    config%final_time = final_time
    config%cfl = cfl
    config%gamma = gamma
    config%output_file = trim(output_file)
    config%limiter = trim(limiter)
    config%boundary_condition = trim(boundary_condition)
    config%riemann_solver = trim(riemann_solver)
    config%use_transverse_correction = use_transverse_correction

    vortex%center_x = center_x
    vortex%center_y = center_y
    vortex%strength = strength
    vortex%base_density = base_density
    vortex%base_pressure = base_pressure
    vortex%base_velocity_x = base_velocity_x
    vortex%base_velocity_y = base_velocity_y

    call validate_configuration_2d(config, vortex, ok, message)
  end subroutine read_configuration_2d

  pure subroutine validate_configuration_2d(config, vortex, ok, message)
    type(simulation_config_2d), intent(in) :: config
    type(isentropic_vortex_config), intent(in) :: vortex
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    real(dp) :: base_temperature, minimum_temperature, pi

    ok = .false.
    message = ""

    pi = acos(-1.0_dp)
    base_temperature = vortex%base_pressure / vortex%base_density
    minimum_temperature = base_temperature - &
      (config%gamma - 1.0_dp) * vortex%strength**2 * exp(1.0_dp) / &
      (8.0_dp * config%gamma * pi**2)

    if (config%nx < 8 .or. config%ny < 8) then
      message = "nx and ny must both be at least 8"
    else if (config%x_max <= config%x_min .or. &
             config%y_max <= config%y_min) then
      message = "2D domain maxima must exceed minima"
    else if (config%final_time <= 0.0_dp) then
      message = "final_time must be positive"
    else if (config%cfl <= 0.0_dp .or. config%cfl > 1.0_dp) then
      message = "cfl must be in (0, 1]"
    else if (config%gamma <= 1.0_dp) then
      message = "gamma must be greater than 1"
    else if (config%max_steps <= 0) then
      message = "max_steps must be positive"
    else if (.not. valid_limiter_2d(config%limiter)) then
      message = "limiter must be minmod or mc"
    else if (trim(config%boundary_condition) /= "periodic") then
      message = "the current 2D solver requires periodic boundaries"
    else if (.not. valid_riemann_solver_2d(config%riemann_solver)) then
      message = "riemann_solver must be rusanov or pelec"
    else if (vortex%base_density <= 0.0_dp .or. &
             vortex%base_pressure <= 0.0_dp) then
      message = "vortex base density and pressure must be positive"
    else if (vortex%center_x < config%x_min .or. &
             vortex%center_x >= config%x_max .or. &
             vortex%center_y < config%y_min .or. &
             vortex%center_y >= config%y_max) then
      message = "vortex center must lie inside the periodic domain"
    else if (minimum_temperature <= 0.0_dp) then
      message = "vortex strength produces non-positive temperature"
    else
      ok = .true.
    end if
  end subroutine validate_configuration_2d

  pure logical function valid_limiter_2d(name) result(valid)
    character(len=*), intent(in) :: name

    select case (trim(name))
    case ("minmod", "mc")
      valid = .true.
    case default
      valid = .false.
    end select
  end function valid_limiter_2d

  pure logical function valid_riemann_solver_2d(name) result(valid)
    character(len=*), intent(in) :: name

    select case (trim(name))
    case ("rusanov", "pelec")
      valid = .true.
    case default
      valid = .false.
    end select
  end function valid_riemann_solver_2d

end module simulation_config_2d_mod
