module simulation_config_mod
  use precision_mod, only: dp
  use constants_mod, only: default_gamma
  implicit none
  private

  type, public :: simulation_config
    integer :: nx = 400
    integer :: max_steps = 100000
    real(dp) :: x_min = 0.0_dp
    real(dp) :: x_max = 1.0_dp
    real(dp) :: final_time = 0.2_dp
    real(dp) :: cfl = 0.45_dp
    real(dp) :: gamma = default_gamma
    character(len=512) :: output_file = "sod.csv"
    character(len=32) :: problem = "sod"
    character(len=32) :: reconstruction = "pcm"
    character(len=32) :: limiter = "mc"
    character(len=32) :: boundary_condition = "outflow"
    character(len=32) :: riemann_solver = "rusanov"
  end type simulation_config

  type, public :: sod_config
    real(dp) :: discontinuity = 0.5_dp
    real(dp) :: rho_left = 1.0_dp
    real(dp) :: velocity_left = 0.0_dp
    real(dp) :: pressure_left = 1.0_dp
    real(dp) :: rho_right = 0.125_dp
    real(dp) :: velocity_right = 0.0_dp
    real(dp) :: pressure_right = 0.1_dp
  end type sod_config

  type, public :: shu_osher_config
    real(dp) :: shock_location = -4.0_dp
    real(dp) :: left_density = 3.857143_dp
    real(dp) :: left_velocity = 2.629369_dp
    real(dp) :: left_pressure = 10.33333_dp
    real(dp) :: density_base = 1.0_dp
    real(dp) :: density_amplitude = 0.2_dp
    real(dp) :: density_wavenumber = 5.0_dp
    real(dp) :: right_velocity = 0.0_dp
    real(dp) :: right_pressure = 1.0_dp
  end type shu_osher_config

  public :: read_configuration
  public :: read_sod_configuration
  public :: validate_configuration

  interface validate_configuration
    module procedure validate_configuration_full
    module procedure validate_configuration_sod
  end interface validate_configuration

contains

  subroutine read_configuration(path, config, sod, shu_osher, ok, message)
    character(len=*), intent(in) :: path
    type(simulation_config), intent(out) :: config
    type(sod_config), intent(out) :: sod
    type(shu_osher_config), intent(out) :: shu_osher
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    integer :: nx, max_steps, unit, io_status
    real(dp) :: x_min, x_max, final_time, cfl, gamma
    real(dp) :: discontinuity, rho_left, velocity_left, pressure_left
    real(dp) :: rho_right, velocity_right, pressure_right
    real(dp) :: shock_location, left_density, left_velocity, left_pressure
    real(dp) :: density_base, density_amplitude, density_wavenumber
    character(len=512) :: output_file
    character(len=32) :: problem, reconstruction, limiter
    character(len=32) :: boundary_condition, riemann_solver

    namelist /simulation/ nx, max_steps, x_min, x_max, final_time, cfl, &
      gamma, output_file, problem, reconstruction, limiter, &
      boundary_condition, riemann_solver
    namelist /sod_problem/ discontinuity, rho_left, velocity_left, &
      pressure_left, rho_right, velocity_right, pressure_right
    namelist /shu_osher_problem/ shock_location, left_density, &
      left_velocity, left_pressure, density_base, density_amplitude, &
      density_wavenumber, velocity_right, pressure_right

    config = simulation_config()
    sod = sod_config()
    shu_osher = shu_osher_config()

    nx = config%nx
    max_steps = config%max_steps
    x_min = config%x_min
    x_max = config%x_max
    final_time = config%final_time
    cfl = config%cfl
    gamma = config%gamma
    output_file = config%output_file
    problem = config%problem
    reconstruction = config%reconstruction
    limiter = config%limiter
    boundary_condition = config%boundary_condition
    riemann_solver = config%riemann_solver

    discontinuity = sod%discontinuity
    rho_left = sod%rho_left
    velocity_left = sod%velocity_left
    pressure_left = sod%pressure_left
    rho_right = sod%rho_right
    velocity_right = sod%velocity_right
    pressure_right = sod%pressure_right

    shock_location = shu_osher%shock_location
    left_density = shu_osher%left_density
    left_velocity = shu_osher%left_velocity
    left_pressure = shu_osher%left_pressure
    density_base = shu_osher%density_base
    density_amplitude = shu_osher%density_amplitude
    density_wavenumber = shu_osher%density_wavenumber

    open(newunit=unit, file=trim(path), status="old", action="read", &
      iostat=io_status)
    if (io_status /= 0) then
      ok = .false.
      write(message, '(a,1x,a)') "Could not open input file:", trim(path)
      return
    end if

    read(unit, nml=simulation, iostat=io_status)
    if (io_status /= 0) then
      close(unit)
      ok = .false.
      write(message, '(a,1x,a)') &
        "Could not read &simulation from:", trim(path)
      return
    end if

    select case (trim(problem))
    case ("sod")
      read(unit, nml=sod_problem, iostat=io_status)
      if (io_status /= 0) then
        close(unit)
        ok = .false.
        write(message, '(a,1x,a)') &
          "Could not read &sod_problem from:", trim(path)
        return
      end if

    case ("shu_osher")
      velocity_right = shu_osher%right_velocity
      pressure_right = shu_osher%right_pressure
      read(unit, nml=shu_osher_problem, iostat=io_status)
      if (io_status /= 0) then
        close(unit)
        ok = .false.
        write(message, '(a,1x,a)') &
          "Could not read &shu_osher_problem from:", trim(path)
        return
      end if

    case default
      close(unit)
      ok = .false.
      write(message, '(a,1x,a)') "Unsupported problem:", trim(problem)
      return
    end select
    close(unit)

    config%nx = nx
    config%max_steps = max_steps
    config%x_min = x_min
    config%x_max = x_max
    config%final_time = final_time
    config%cfl = cfl
    config%gamma = gamma
    config%output_file = trim(output_file)
    config%problem = trim(problem)
    config%reconstruction = trim(reconstruction)
    config%limiter = trim(limiter)
    config%boundary_condition = trim(boundary_condition)
    config%riemann_solver = trim(riemann_solver)

    sod%discontinuity = discontinuity
    sod%rho_left = rho_left
    sod%velocity_left = velocity_left
    sod%pressure_left = pressure_left
    sod%rho_right = rho_right
    sod%velocity_right = velocity_right
    sod%pressure_right = pressure_right

    shu_osher%shock_location = shock_location
    shu_osher%left_density = left_density
    shu_osher%left_velocity = left_velocity
    shu_osher%left_pressure = left_pressure
    shu_osher%density_base = density_base
    shu_osher%density_amplitude = density_amplitude
    shu_osher%density_wavenumber = density_wavenumber
    shu_osher%right_velocity = velocity_right
    shu_osher%right_pressure = pressure_right

    call validate_configuration_full(config, sod, shu_osher, ok, message)
  end subroutine read_configuration

  subroutine read_sod_configuration(path, config, sod, ok, message)
    character(len=*), intent(in) :: path
    type(simulation_config), intent(out) :: config
    type(sod_config), intent(out) :: sod
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    type(shu_osher_config) :: shu_osher

    call read_configuration(path, config, sod, shu_osher, ok, message)
    if (ok .and. trim(config%problem) /= "sod") then
      ok = .false.
      message = "Input does not select problem = sod"
    end if
  end subroutine read_sod_configuration

  pure subroutine validate_configuration_full( &
      config, sod, shu_osher, ok, message)
    type(simulation_config), intent(in) :: config
    type(sod_config), intent(in) :: sod
    type(shu_osher_config), intent(in) :: shu_osher
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    ok = .false.
    message = ""

    if (config%nx < 10) then
      message = "nx must be at least 10"
    else if (config%x_max <= config%x_min) then
      message = "x_max must be greater than x_min"
    else if (config%final_time <= 0.0_dp) then
      message = "final_time must be positive"
    else if (config%cfl <= 0.0_dp .or. config%cfl > 1.0_dp) then
      message = "cfl must be in (0, 1]"
    else if (config%gamma <= 1.0_dp) then
      message = "gamma must be greater than 1"
    else if (config%max_steps <= 0) then
      message = "max_steps must be positive"
    else if (.not. valid_problem(config%problem)) then
      message = "problem must be sod or shu_osher"
    else if (.not. valid_reconstruction(config%reconstruction)) then
      message = "reconstruction must be pcm, plm, or pelec_plm"
    else if (.not. valid_limiter(config%limiter)) then
      message = "limiter must be minmod or mc"
    else if (.not. valid_boundary_condition(config%boundary_condition)) then
      message = "boundary_condition must be outflow or periodic"
    else if (.not. valid_riemann_solver(config%riemann_solver)) then
      message = "riemann_solver must be rusanov or pelec"
    else
      select case (trim(config%problem))
      case ("sod")
        call validate_sod(config, sod, ok, message)
      case ("shu_osher")
        call validate_shu_osher(config, shu_osher, ok, message)
      end select
    end if
  end subroutine validate_configuration_full

  pure subroutine validate_configuration_sod(config, sod, ok, message)
    type(simulation_config), intent(in) :: config
    type(sod_config), intent(in) :: sod
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    type(shu_osher_config) :: shu_osher

    shu_osher = shu_osher_config()
    call validate_configuration_full(config, sod, shu_osher, ok, message)
  end subroutine validate_configuration_sod

  pure subroutine validate_sod(config, sod, ok, message)
    type(simulation_config), intent(in) :: config
    type(sod_config), intent(in) :: sod
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    ok = .false.
    if (sod%discontinuity <= config%x_min .or. &
        sod%discontinuity >= config%x_max) then
      message = "discontinuity must lie inside the domain"
    else if (sod%rho_left <= 0.0_dp .or. sod%rho_right <= 0.0_dp) then
      message = "left and right densities must be positive"
    else if (sod%pressure_left <= 0.0_dp .or. &
             sod%pressure_right <= 0.0_dp) then
      message = "left and right pressures must be positive"
    else
      ok = .true.
      message = ""
    end if
  end subroutine validate_sod

  pure subroutine validate_shu_osher(config, shu_osher, ok, message)
    type(simulation_config), intent(in) :: config
    type(shu_osher_config), intent(in) :: shu_osher
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    ok = .false.
    if (shu_osher%shock_location <= config%x_min .or. &
        shu_osher%shock_location >= config%x_max) then
      message = "shock_location must lie inside the domain"
    else if (shu_osher%left_density <= 0.0_dp .or. &
             shu_osher%left_pressure <= 0.0_dp) then
      message = "Shu-Osher left density and pressure must be positive"
    else if (shu_osher%density_base <= &
             abs(shu_osher%density_amplitude)) then
      message = "Shu-Osher density_base must exceed abs(density_amplitude)"
    else if (shu_osher%density_wavenumber <= 0.0_dp) then
      message = "Shu-Osher density_wavenumber must be positive"
    else if (shu_osher%right_pressure <= 0.0_dp) then
      message = "Shu-Osher right pressure must be positive"
    else
      ok = .true.
      message = ""
    end if
  end subroutine validate_shu_osher

  pure logical function valid_problem(name) result(valid)
    character(len=*), intent(in) :: name

    select case (trim(name))
    case ("sod", "shu_osher")
      valid = .true.
    case default
      valid = .false.
    end select
  end function valid_problem

  pure logical function valid_reconstruction(name) result(valid)
    character(len=*), intent(in) :: name

    select case (trim(name))
    case ("pcm", "plm", "pelec_plm")
      valid = .true.
    case default
      valid = .false.
    end select
  end function valid_reconstruction

  pure logical function valid_limiter(name) result(valid)
    character(len=*), intent(in) :: name

    select case (trim(name))
    case ("minmod", "mc")
      valid = .true.
    case default
      valid = .false.
    end select
  end function valid_limiter

  pure logical function valid_boundary_condition(name) result(valid)
    character(len=*), intent(in) :: name

    select case (trim(name))
    case ("outflow", "periodic")
      valid = .true.
    case default
      valid = .false.
    end select
  end function valid_boundary_condition

  pure logical function valid_riemann_solver(name) result(valid)
    character(len=*), intent(in) :: name

    select case (trim(name))
    case ("rusanov", "pelec")
      valid = .true.
    case default
      valid = .false.
    end select
  end function valid_riemann_solver

end module simulation_config_mod
