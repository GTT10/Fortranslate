module simulation_config_h2o2_full_mod
  use precision_mod, only: dp
  implicit none
  private

  type, public :: h2o2_full_reactor_config
    real(dp) :: final_time = 2.0e-3_dp
    real(dp) :: output_interval = 2.0e-5_dp
    real(dp) :: initial_time_step = 1.0e-9_dp
    real(dp) :: minimum_time_step = 1.0e-14_dp
    real(dp) :: maximum_time_step = 2.0e-5_dp
    real(dp) :: relative_tolerance = 2.0e-7_dp
    real(dp) :: absolute_tolerance = 1.0e-12_dp
    integer :: maximum_steps = 500000
    integer :: maximum_newton_iterations = 16
    real(dp) :: initial_temperature = 1000.0_dp
    real(dp) :: initial_pressure = 101325.0_dp
    real(dp) :: x_h2 = 2.0_dp
    real(dp) :: x_h = 1.0e-6_dp
    real(dp) :: x_o = 1.0e-12_dp
    real(dp) :: x_o2 = 1.0_dp
    real(dp) :: x_oh = 1.0e-12_dp
    real(dp) :: x_h2o = 0.0_dp
    real(dp) :: x_ho2 = 0.0_dp
    real(dp) :: x_h2o2 = 0.0_dp
    real(dp) :: x_ar = 0.0_dp
    real(dp) :: x_n2 = 3.0_dp
    character(len=256) :: output_file = "zero_d_h2o2_full.csv"
  end type h2o2_full_reactor_config

  public :: read_h2o2_full_reactor_configuration

contains

  subroutine read_h2o2_full_reactor_configuration(path, config, ok, message)
    character(len=*), intent(in) :: path
    type(h2o2_full_reactor_config), intent(out) :: config
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    real(dp) :: final_time, output_interval
    real(dp) :: initial_time_step, minimum_time_step, maximum_time_step
    real(dp) :: relative_tolerance, absolute_tolerance
    real(dp) :: initial_temperature, initial_pressure
    real(dp) :: x_h2, x_h, x_o, x_o2, x_oh, x_h2o
    real(dp) :: x_ho2, x_h2o2, x_ar, x_n2
    integer :: maximum_steps, maximum_newton_iterations
    character(len=256) :: output_file
    integer :: unit, io_status
    namelist /h2o2_full_reactor/ &
      final_time, output_interval, initial_time_step, minimum_time_step, &
      maximum_time_step, relative_tolerance, absolute_tolerance, &
      maximum_steps, maximum_newton_iterations, initial_temperature, &
      initial_pressure, x_h2, x_h, x_o, x_o2, x_oh, x_h2o, x_ho2, &
      x_h2o2, x_ar, x_n2, output_file

    final_time = config%final_time
    output_interval = config%output_interval
    initial_time_step = config%initial_time_step
    minimum_time_step = config%minimum_time_step
    maximum_time_step = config%maximum_time_step
    relative_tolerance = config%relative_tolerance
    absolute_tolerance = config%absolute_tolerance
    maximum_steps = config%maximum_steps
    maximum_newton_iterations = config%maximum_newton_iterations
    initial_temperature = config%initial_temperature
    initial_pressure = config%initial_pressure
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

    open(newunit=unit, file=trim(path), status="old", action="read", &
      iostat=io_status)
    if (io_status /= 0) then
      ok = .false.
      message = "Could not open full H2/O2 reactor input"
      return
    end if
    read(unit, nml=h2o2_full_reactor, iostat=io_status)
    close(unit)
    if (io_status /= 0) then
      ok = .false.
      message = "Could not parse &h2o2_full_reactor namelist"
      return
    end if

    config%final_time = final_time
    config%output_interval = output_interval
    config%initial_time_step = initial_time_step
    config%minimum_time_step = minimum_time_step
    config%maximum_time_step = maximum_time_step
    config%relative_tolerance = relative_tolerance
    config%absolute_tolerance = absolute_tolerance
    config%maximum_steps = maximum_steps
    config%maximum_newton_iterations = maximum_newton_iterations
    config%initial_temperature = initial_temperature
    config%initial_pressure = initial_pressure
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
    config%output_file = output_file

    ok = final_time > 0.0_dp .and. output_interval > 0.0_dp
    ok = ok .and. initial_time_step > 0.0_dp
    ok = ok .and. minimum_time_step > 0.0_dp
    ok = ok .and. maximum_time_step >= minimum_time_step
    ok = ok .and. relative_tolerance > 0.0_dp
    ok = ok .and. absolute_tolerance > 0.0_dp
    ok = ok .and. maximum_steps > 0 .and. maximum_newton_iterations > 0
    ok = ok .and. initial_temperature >= 300.0_dp
    ok = ok .and. initial_temperature <= 3500.0_dp
    ok = ok .and. initial_pressure > 0.0_dp
    ok = ok .and. min( &
      x_h2, x_h, x_o, x_o2, x_oh, x_h2o, x_ho2, x_h2o2, x_ar, x_n2) &
      >= 0.0_dp
    ok = ok .and. x_h2 + x_h + x_o + x_o2 + x_oh + x_h2o + &
      x_ho2 + x_h2o2 + x_ar + x_n2 > 0.0_dp
    if (ok) then
      message = ""
    else
      message = "Invalid full H2/O2 reactor configuration"
    end if
  end subroutine read_h2o2_full_reactor_configuration

end module simulation_config_h2o2_full_mod
