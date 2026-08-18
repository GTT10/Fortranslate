program pelef0d
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_toy_isomerization_thermo
  use mixture_thermo_mod, only: mixture_pressure
  use isomerization_reactor_mod, only: &
    isomerization_reaction, advance_isomerization_rk4, &
    reactor_internal_energy, isomerization_rate_constant
  use simulation_config_reactor_mod, only: &
    reactor_config, read_reactor_configuration
  implicit none

  type(reactor_config) :: config
  type(isomerization_reaction) :: reaction
  type(nasa7_species), allocatable :: species(:)
  real(dp) :: mass_fractions(2), temperature, target_energy
  real(dp) :: current_energy, pressure, rate_constant
  real(dp) :: time, dt, relative_energy_error
  character(len=1024) :: input_path, message
  logical :: ok
  integer :: argument_count, unit, io_status, step

  argument_count = command_argument_count()
  if (argument_count /= 1) then
    write(*, '(a)') "Usage: pelef0d <input.nml>"
    error stop 2
  end if

  call get_command_argument(1, input_path)
  call read_reactor_configuration(trim(input_path), config, ok, message)
  if (.not. ok) then
    write(*, '(a)') trim(message)
    error stop 2
  end if

  call load_toy_isomerization_thermo(species, ok)
  if (.not. ok) error stop "Failed to initialize toy NASA7 species"

  reaction%reactant = 1
  reaction%product = 2
  reaction%pre_exponential = config%pre_exponential
  reaction%temperature_exponent = config%temperature_exponent
  reaction%activation_temperature = config%activation_temperature

  mass_fractions(1) = config%initial_reactant_fraction
  mass_fractions(2) = 1.0_dp - mass_fractions(1)
  temperature = config%initial_temperature
  target_energy = reactor_internal_energy( &
    species, mass_fractions, temperature, ok)
  if (.not. ok) error stop "Failed to evaluate initial reactor energy"

  open(newunit=unit, file=trim(config%output_file), status="replace", &
    action="write", iostat=io_status)
  if (io_status /= 0) error stop "Could not create reactor CSV output"
  write(unit, '(a)') &
    "time,temperature,pressure,Y_A,Y_B,specific_internal_energy," // &
    "relative_energy_error,rate_constant"

  write(*, '(a)') "PeleF " // pelef_version
  write(*, '(a,1x,a)') "Input:", trim(input_path)
  write(*, '(a,l1)') "Adiabatic: ", config%adiabatic

  time = 0.0_dp
  step = 0
  call write_reactor_row()
  do while (time < config%final_time)
    dt = min(config%time_step, config%final_time - time)
    call advance_isomerization_rk4( &
      species, reaction, dt, config%adiabatic, target_energy, &
      mass_fractions, temperature, ok)
    if (.not. ok) error stop "Zero-dimensional reactor step failed"
    time = time + dt
    step = step + 1
    if (mod(step, config%output_stride) == 0 .or. &
        time >= config%final_time) call write_reactor_row()
  end do
  close(unit)

  current_energy = reactor_internal_energy( &
    species, mass_fractions, temperature, ok)
  if (.not. ok) error stop "Failed to evaluate final reactor energy"
  relative_energy_error = abs(current_energy - target_energy) / &
    max(1.0_dp, abs(target_energy))

  write(*, '(a,i0)') "Completed steps: ", step
  write(*, '(a,es24.16)') "Final time: ", time
  write(*, '(a,es24.16)') "Final temperature: ", temperature
  write(*, '(a,es24.16)') "Final reactant fraction: ", mass_fractions(1)
  write(*, '(a,es24.16)') "Relative energy error: ", relative_energy_error
  write(*, '(a,1x,a)') "Output:", trim(config%output_file)

contains

  subroutine write_reactor_row()
    current_energy = reactor_internal_energy( &
      species, mass_fractions, temperature, ok)
    if (.not. ok) error stop "Failed to evaluate reactor energy"
    pressure = mixture_pressure( &
      species, mass_fractions, config%density, temperature, ok)
    if (.not. ok) error stop "Failed to evaluate reactor pressure"
    rate_constant = isomerization_rate_constant(reaction, temperature, ok)
    if (.not. ok) error stop "Failed to evaluate reactor rate constant"
    relative_energy_error = abs(current_energy - target_energy) / &
      max(1.0_dp, abs(target_energy))
    write(unit, '(*(es25.16e3,:,","))') &
      time, temperature, pressure, mass_fractions(1), mass_fractions(2), &
      current_energy, relative_energy_error, rate_constant
  end subroutine write_reactor_row

end program pelef0d
