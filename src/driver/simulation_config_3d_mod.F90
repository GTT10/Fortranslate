module simulation_config_3d_mod
  use precision_mod, only: dp
  use constants_mod, only: default_gamma
  implicit none
  private

  type, public :: simulation_config_3d
    integer :: nx = 16
    integer :: ny = 16
    integer :: nz = 16
    integer :: max_steps = 100000
    real(dp) :: x_min = 0.0_dp
    real(dp) :: x_max = 1.0_dp
    real(dp) :: y_min = 0.0_dp
    real(dp) :: y_max = 1.0_dp
    real(dp) :: z_min = 0.0_dp
    real(dp) :: z_max = 1.0_dp
    real(dp) :: final_time = 0.05_dp
    real(dp) :: cfl = 0.4_dp
    real(dp) :: gamma = default_gamma
    character(len=512) :: output_file = "entropy_wave_3d.csv"
    character(len=32) :: boundary_condition = "periodic"
    character(len=32) :: riemann_solver = "rusanov"
  end type simulation_config_3d

  type, public :: entropy_wave_3d_config
    real(dp) :: base_density = 1.0_dp
    real(dp) :: density_amplitude = 0.1_dp
    real(dp) :: base_pressure = 1.0_dp
    real(dp) :: velocity_x = 0.7_dp
    real(dp) :: velocity_y = 0.2_dp
    real(dp) :: velocity_z = -0.1_dp
    integer :: wave_number_x = 1
    integer :: wave_number_y = 1
    integer :: wave_number_z = 1
  end type entropy_wave_3d_config

  public :: read_configuration_3d
  public :: validate_configuration_3d

contains

  subroutine read_configuration_3d(path, config, wave, ok, message)
    character(len=*), intent(in) :: path
    type(simulation_config_3d), intent(out) :: config
    type(entropy_wave_3d_config), intent(out) :: wave
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    integer :: nx, ny, nz, max_steps, unit, io_status
    integer :: wave_number_x, wave_number_y, wave_number_z
    real(dp) :: x_min, x_max, y_min, y_max, z_min, z_max
    real(dp) :: final_time, cfl, gamma
    real(dp) :: base_density, density_amplitude, base_pressure
    real(dp) :: velocity_x, velocity_y, velocity_z
    character(len=512) :: output_file
    character(len=32) :: boundary_condition, riemann_solver

    namelist /simulation_3d/ &
      nx, ny, nz, max_steps, x_min, x_max, y_min, y_max, z_min, z_max, &
      final_time, cfl, gamma, output_file, boundary_condition, riemann_solver
    namelist /entropy_wave_3d/ &
      base_density, density_amplitude, base_pressure, &
      velocity_x, velocity_y, velocity_z, &
      wave_number_x, wave_number_y, wave_number_z

    config = simulation_config_3d()
    wave = entropy_wave_3d_config()

    nx = config%nx
    ny = config%ny
    nz = config%nz
    max_steps = config%max_steps
    x_min = config%x_min
    x_max = config%x_max
    y_min = config%y_min
    y_max = config%y_max
    z_min = config%z_min
    z_max = config%z_max
    final_time = config%final_time
    cfl = config%cfl
    gamma = config%gamma
    output_file = config%output_file
    boundary_condition = config%boundary_condition
    riemann_solver = config%riemann_solver

    base_density = wave%base_density
    density_amplitude = wave%density_amplitude
    base_pressure = wave%base_pressure
    velocity_x = wave%velocity_x
    velocity_y = wave%velocity_y
    velocity_z = wave%velocity_z
    wave_number_x = wave%wave_number_x
    wave_number_y = wave%wave_number_y
    wave_number_z = wave%wave_number_z

    open(newunit=unit, file=trim(path), status="old", action="read", &
      iostat=io_status)
    if (io_status /= 0) then
      ok = .false.
      write(message, '(a,1x,a)') "Could not open input file:", trim(path)
      return
    end if

    read(unit, nml=simulation_3d, iostat=io_status)
    if (io_status /= 0) then
      close(unit)
      ok = .false.
      write(message, '(a,1x,a)') &
        "Could not read &simulation_3d from:", trim(path)
      return
    end if
    read(unit, nml=entropy_wave_3d, iostat=io_status)
    close(unit)
    if (io_status /= 0) then
      ok = .false.
      write(message, '(a,1x,a)') &
        "Could not read &entropy_wave_3d from:", trim(path)
      return
    end if

    config%nx = nx
    config%ny = ny
    config%nz = nz
    config%max_steps = max_steps
    config%x_min = x_min
    config%x_max = x_max
    config%y_min = y_min
    config%y_max = y_max
    config%z_min = z_min
    config%z_max = z_max
    config%final_time = final_time
    config%cfl = cfl
    config%gamma = gamma
    config%output_file = trim(output_file)
    config%boundary_condition = trim(boundary_condition)
    config%riemann_solver = trim(riemann_solver)

    wave%base_density = base_density
    wave%density_amplitude = density_amplitude
    wave%base_pressure = base_pressure
    wave%velocity_x = velocity_x
    wave%velocity_y = velocity_y
    wave%velocity_z = velocity_z
    wave%wave_number_x = wave_number_x
    wave%wave_number_y = wave_number_y
    wave%wave_number_z = wave_number_z

    call validate_configuration_3d(config, wave, ok, message)
  end subroutine read_configuration_3d

  pure subroutine validate_configuration_3d(config, wave, ok, message)
    type(simulation_config_3d), intent(in) :: config
    type(entropy_wave_3d_config), intent(in) :: wave
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    ok = .false.
    message = ""
    if (config%nx < 4 .or. config%ny < 4 .or. config%nz < 4) then
      message = "nx, ny, and nz must all be at least 4"
    else if (config%x_max <= config%x_min .or. &
             config%y_max <= config%y_min .or. &
             config%z_max <= config%z_min) then
      message = "3D domain maxima must exceed minima"
    else if (config%final_time <= 0.0_dp) then
      message = "final_time must be positive"
    else if (config%cfl <= 0.0_dp .or. config%cfl > 1.0_dp) then
      message = "cfl must be in (0, 1]"
    else if (config%gamma <= 1.0_dp) then
      message = "gamma must be greater than 1"
    else if (config%max_steps <= 0) then
      message = "max_steps must be positive"
    else if (len_trim(config%output_file) == 0) then
      message = "output_file must not be empty"
    else if (trim(config%boundary_condition) /= "periodic") then
      message = "the current 3D solver requires periodic boundaries"
    else if (.not. valid_riemann_solver_3d(config%riemann_solver)) then
      message = "riemann_solver must be rusanov or pelec"
    else if (wave%base_density <= 0.0_dp .or. &
             wave%base_pressure <= 0.0_dp) then
      message = "wave base density and pressure must be positive"
    else if (abs(wave%density_amplitude) >= wave%base_density) then
      message = "density amplitude must be smaller than base density"
    else if (wave%wave_number_x == 0 .and. &
             wave%wave_number_y == 0 .and. &
             wave%wave_number_z == 0) then
      message = "at least one wave number must be nonzero"
    else
      ok = .true.
    end if
  end subroutine validate_configuration_3d

  pure logical function valid_riemann_solver_3d(name) result(valid)
    character(len=*), intent(in) :: name

    select case (trim(name))
    case ("rusanov", "pelec")
      valid = .true.
    case default
      valid = .false.
    end select
  end function valid_riemann_solver_3d

end module simulation_config_3d_mod
