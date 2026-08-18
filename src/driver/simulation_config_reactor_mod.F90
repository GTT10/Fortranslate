module simulation_config_reactor_mod
  use precision_mod, only: dp
  implicit none
  private

  type, public :: reactor_config
    real(dp) :: final_time = 0.2_dp
    real(dp) :: time_step = 5.0e-5_dp
    real(dp) :: density = 1.0_dp
    real(dp) :: initial_temperature = 700.0_dp
    real(dp) :: initial_reactant_fraction = 0.95_dp
    real(dp) :: pre_exponential = 2.0e4_dp
    real(dp) :: temperature_exponent = 0.0_dp
    real(dp) :: activation_temperature = 5000.0_dp
    logical :: adiabatic = .true.
    integer :: output_stride = 20
    character(len=512) :: output_file = "zero_d_isomerization.csv"
  end type reactor_config

  public :: read_reactor_configuration
  public :: validate_reactor_configuration

contains

  subroutine read_reactor_configuration(path, config, ok, message)
    character(len=*), intent(in) :: path
    type(reactor_config), intent(out) :: config
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    real(dp) :: final_time, time_step, density, initial_temperature
    real(dp) :: initial_reactant_fraction, pre_exponential
    real(dp) :: temperature_exponent, activation_temperature
    logical :: adiabatic
    integer :: output_stride, unit, io_status
    character(len=512) :: output_file

    namelist /reactor/ final_time, time_step, density, initial_temperature, &
      initial_reactant_fraction, pre_exponential, temperature_exponent, &
      activation_temperature, adiabatic, output_stride, output_file

    config = reactor_config()
    final_time = config%final_time
    time_step = config%time_step
    density = config%density
    initial_temperature = config%initial_temperature
    initial_reactant_fraction = config%initial_reactant_fraction
    pre_exponential = config%pre_exponential
    temperature_exponent = config%temperature_exponent
    activation_temperature = config%activation_temperature
    adiabatic = config%adiabatic
    output_stride = config%output_stride
    output_file = config%output_file

    open(newunit=unit, file=trim(path), status="old", action="read", &
      iostat=io_status)
    if (io_status /= 0) then
      ok = .false.
      write(message, '(a,1x,a)') "Could not open reactor input:", trim(path)
      return
    end if
    read(unit, nml=reactor, iostat=io_status)
    close(unit)
    if (io_status /= 0) then
      ok = .false.
      write(message, '(a,1x,a)') "Could not read &reactor from:", trim(path)
      return
    end if

    config%final_time = final_time
    config%time_step = time_step
    config%density = density
    config%initial_temperature = initial_temperature
    config%initial_reactant_fraction = initial_reactant_fraction
    config%pre_exponential = pre_exponential
    config%temperature_exponent = temperature_exponent
    config%activation_temperature = activation_temperature
    config%adiabatic = adiabatic
    config%output_stride = output_stride
    config%output_file = trim(output_file)

    call validate_reactor_configuration(config, ok, message)
  end subroutine read_reactor_configuration

  subroutine validate_reactor_configuration(config, ok, message)
    type(reactor_config), intent(in) :: config
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    ok = .false.
    message = ""
    if (config%final_time <= 0.0_dp) then
      message = "final_time must be positive"
    else if (config%time_step <= 0.0_dp .or. &
             config%time_step > config%final_time) then
      message = "time_step must be positive and no larger than final_time"
    else if (config%density <= 0.0_dp) then
      message = "density must be positive"
    else if (config%initial_temperature < 200.0_dp .or. &
             config%initial_temperature > 5000.0_dp) then
      message = "initial_temperature must be within [200, 5000] K"
    else if (config%initial_reactant_fraction < 0.0_dp .or. &
             config%initial_reactant_fraction > 1.0_dp) then
      message = "initial_reactant_fraction must be within [0, 1]"
    else if (config%pre_exponential < 0.0_dp) then
      message = "pre_exponential must be non-negative"
    else if (config%activation_temperature < 0.0_dp) then
      message = "activation_temperature must be non-negative"
    else if (config%output_stride < 1) then
      message = "output_stride must be positive"
    else
      ok = .true.
    end if
  end subroutine validate_reactor_configuration

end module simulation_config_reactor_mod
