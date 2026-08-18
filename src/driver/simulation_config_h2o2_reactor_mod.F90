module simulation_config_h2o2_reactor_mod
  use precision_mod, only: dp
  implicit none
  private

  type, public :: h2o2_reactor_config
    real(dp) :: final_time = 2.0e-4_dp
    real(dp) :: initial_time_step = 1.0e-9_dp
    real(dp) :: maximum_time_step = 2.0e-6_dp
    real(dp) :: output_interval = 2.0e-6_dp
    real(dp) :: relative_tolerance = 1.0e-8_dp
    real(dp) :: absolute_tolerance = 1.0e-14_dp
    real(dp) :: initial_temperature = 1200.0_dp
    real(dp) :: initial_pressure = 101325.0_dp
    real(dp) :: x_h2 = 0.295_dp
    real(dp) :: x_h = 1.0e-6_dp
    real(dp) :: x_o = 1.0e-6_dp
    real(dp) :: x_o2 = 0.1475_dp
    real(dp) :: x_oh = 1.0e-6_dp
    real(dp) :: x_h2o = 0.0_dp
    real(dp) :: x_n2 = 0.557497_dp
    integer :: maximum_steps = 2000000
    character(len=512) :: output_file = "zero_d_h2o2.csv"
  end type h2o2_reactor_config

  public :: read_h2o2_reactor_configuration
  public :: validate_h2o2_reactor_configuration

contains

  subroutine read_h2o2_reactor_configuration(path, config, ok, message)
    character(len=*), intent(in) :: path
    type(h2o2_reactor_config), intent(out) :: config
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    real(dp) :: final_time, initial_time_step, maximum_time_step
    real(dp) :: output_interval, relative_tolerance, absolute_tolerance
    real(dp) :: initial_temperature, initial_pressure
    real(dp) :: x_h2, x_h, x_o, x_o2, x_oh, x_h2o, x_n2
    integer :: maximum_steps, unit, io_status
    character(len=512) :: output_file

    namelist /h2o2_reactor/ final_time, initial_time_step, maximum_time_step, &
      output_interval, relative_tolerance, absolute_tolerance, &
      initial_temperature, initial_pressure, x_h2, x_h, x_o, x_o2, x_oh, &
      x_h2o, x_n2, maximum_steps, output_file

    config = h2o2_reactor_config()
    final_time = config%final_time
    initial_time_step = config%initial_time_step
    maximum_time_step = config%maximum_time_step
    output_interval = config%output_interval
    relative_tolerance = config%relative_tolerance
    absolute_tolerance = config%absolute_tolerance
    initial_temperature = config%initial_temperature
    initial_pressure = config%initial_pressure
    x_h2 = config%x_h2
    x_h = config%x_h
    x_o = config%x_o
    x_o2 = config%x_o2
    x_oh = config%x_oh
    x_h2o = config%x_h2o
    x_n2 = config%x_n2
    maximum_steps = config%maximum_steps
    output_file = config%output_file

    open(newunit=unit, file=trim(path), status="old", action="read", &
      iostat=io_status)
    if (io_status /= 0) then
      ok = .false.
      write(message, '(a,1x,a)') "Could not open H2/O2 input:", trim(path)
      return
    end if
    read(unit, nml=h2o2_reactor, iostat=io_status)
    close(unit)
    if (io_status /= 0) then
      ok = .false.
      write(message, '(a,1x,a)') "Could not read &h2o2_reactor from:", &
        trim(path)
      return
    end if

    config%final_time = final_time
    config%initial_time_step = initial_time_step
    config%maximum_time_step = maximum_time_step
    config%output_interval = output_interval
    config%relative_tolerance = relative_tolerance
    config%absolute_tolerance = absolute_tolerance
    config%initial_temperature = initial_temperature
    config%initial_pressure = initial_pressure
    config%x_h2 = x_h2
    config%x_h = x_h
    config%x_o = x_o
    config%x_o2 = x_o2
    config%x_oh = x_oh
    config%x_h2o = x_h2o
    config%x_n2 = x_n2
    config%maximum_steps = maximum_steps
    config%output_file = trim(output_file)
    call validate_h2o2_reactor_configuration(config, ok, message)
  end subroutine read_h2o2_reactor_configuration

  subroutine validate_h2o2_reactor_configuration(config, ok, message)
    type(h2o2_reactor_config), intent(in) :: config
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    real(dp) :: mole_fraction_sum

    ok = .false.
    message = ""
    mole_fraction_sum = config%x_h2 + config%x_h + config%x_o + &
      config%x_o2 + config%x_oh + config%x_h2o + config%x_n2
    if (config%final_time <= 0.0_dp) then
      message = "final_time must be positive"
    else if (config%initial_time_step <= 0.0_dp .or. &
             config%maximum_time_step < config%initial_time_step) then
      message = "time-step bounds are invalid"
    else if (config%output_interval <= 0.0_dp .or. &
             config%output_interval > config%final_time) then
      message = "output_interval must be within the simulation interval"
    else if (config%relative_tolerance <= 0.0_dp .or. &
             config%absolute_tolerance <= 0.0_dp) then
      message = "reactor tolerances must be positive"
    else if (config%initial_temperature < 300.0_dp .or. &
             config%initial_temperature > 3500.0_dp) then
      message = "initial_temperature must be within [300, 3500] K"
    else if (config%initial_pressure <= 0.0_dp) then
      message = "initial_pressure must be positive"
    else if (min(config%x_h2, config%x_h, config%x_o, config%x_o2, &
             config%x_oh, config%x_h2o, config%x_n2) < 0.0_dp) then
      message = "mole fractions must be non-negative"
    else if (abs(mole_fraction_sum - 1.0_dp) > 1.0e-12_dp) then
      message = "mole fractions must sum to one"
    else if (config%maximum_steps < 1) then
      message = "maximum_steps must be positive"
    else
      ok = .true.
    end if
  end subroutine validate_h2o2_reactor_configuration

end module simulation_config_h2o2_reactor_mod
