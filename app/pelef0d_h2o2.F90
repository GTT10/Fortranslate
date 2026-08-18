program pelef0d_h2o2
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density, mixture_pressure
  use elementary_kinetics_mod, only: &
    elementary_reaction, elementary_production_rates
  use h2o2_elementary_mechanism_mod, only: &
    h2o2_nspecies, h2o2_h2_index, h2o2_h_index, h2o2_o_index, &
    h2o2_o2_index, h2o2_oh_index, h2o2_h2o_index, h2o2_n2_index, &
    load_h2o2_elementary_mechanism
  use constant_volume_reactor_mod, only: &
    reactor_specific_internal_energy, advance_constant_volume_adaptive
  use simulation_config_h2o2_reactor_mod, only: &
    h2o2_reactor_config, read_h2o2_reactor_configuration
  implicit none

  type(h2o2_reactor_config) :: config
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp) :: mole_fractions(h2o2_nspecies)
  real(dp) :: mass_fractions(h2o2_nspecies)
  real(dp) :: molar_production_rates(h2o2_nspecies)
  real(dp) :: temperature, density, target_internal_energy, current_energy
  real(dp) :: pressure, relative_energy_error, closure_error
  real(dp) :: time, requested_dt, accepted_dt, next_dt, next_output
  real(dp) :: time_tolerance
  character(len=1024) :: input_path, message
  logical :: ok
  integer :: argument_count, unit, io_status, step, output_count

  argument_count = command_argument_count()
  if (argument_count /= 1) then
    write(*, '(a)') "Usage: pelef0d_h2o2 <input.nml>"
    error stop 2
  end if

  call get_command_argument(1, input_path)
  call read_h2o2_reactor_configuration(trim(input_path), config, ok, message)
  if (.not. ok) then
    write(*, '(a)') trim(message)
    error stop 2
  end if
  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "Failed to load H2/O2 thermodynamics"
  call load_h2o2_elementary_mechanism(reactions, ok)
  if (.not. ok) error stop "Failed to load generated H2/O2 mechanism"

  mole_fractions = [ &
    config%x_h2, config%x_h, config%x_o, config%x_o2, config%x_oh, &
    config%x_h2o, config%x_n2 ]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  if (.not. ok) error stop "Failed to convert initial mole fractions"

  temperature = config%initial_temperature
  density = mixture_density( &
    species, mass_fractions, config%initial_pressure, temperature, ok)
  if (.not. ok) error stop "Failed to evaluate initial density"
  target_internal_energy = reactor_specific_internal_energy( &
    species, mass_fractions, temperature, ok)
  if (.not. ok) error stop "Failed to evaluate initial internal energy"

  open(newunit=unit, file=trim(config%output_file), status="replace", &
    action="write", iostat=io_status)
  if (io_status /= 0) error stop "Could not create H2/O2 CSV output"
  write(unit, '(a)') &
    "time,temperature,pressure,density,specific_internal_energy," // &
    "relative_energy_error,closure_error,Y_H2,Y_H,Y_O,Y_O2,Y_OH," // &
    "Y_H2O,Y_N2,wdot_H2,wdot_H,wdot_O,wdot_O2,wdot_OH,wdot_H2O," // &
    "wdot_N2"

  write(*, '(a)') "PeleF " // pelef_version // " elementary H2/O2"
  write(*, '(a,1x,a)') "Input:", trim(input_path)
  write(*, '(a,es24.16)') "Constant density: ", density

  time = 0.0_dp
  requested_dt = config%initial_time_step
  next_output = 0.0_dp
  step = 0
  output_count = 0
  time_tolerance = 50.0_dp * epsilon(1.0_dp) * &
    max(1.0_dp, config%final_time)

  call write_reactor_row()
  next_output = min(config%output_interval, config%final_time)
  do while (time < config%final_time - time_tolerance)
    if (step >= config%maximum_steps) then
      error stop "Maximum H2/O2 reactor step count reached"
    end if

    requested_dt = min(requested_dt, config%maximum_time_step)
    requested_dt = min(requested_dt, config%final_time - time)
    requested_dt = min(requested_dt, next_output - time)
    if (requested_dt <= 0.0_dp) then
      time = next_output
      call write_reactor_row()
      if (next_output >= config%final_time - time_tolerance) exit
      next_output = min(next_output + config%output_interval, &
        config%final_time)
      cycle
    end if

    call advance_constant_volume_adaptive( &
      species, reactions, density, target_internal_energy, requested_dt, &
      config%relative_tolerance, config%absolute_tolerance, mass_fractions, &
      temperature, accepted_dt, next_dt, ok)
    if (.not. ok) error stop "Adaptive H2/O2 reactor step failed"
    time = time + accepted_dt
    requested_dt = min(config%maximum_time_step, next_dt)
    step = step + 1

    if (time >= next_output - time_tolerance) then
      time = next_output
      call write_reactor_row()
      if (next_output >= config%final_time - time_tolerance) exit
      next_output = min(next_output + config%output_interval, &
        config%final_time)
    end if
  end do
  close(unit)

  current_energy = reactor_specific_internal_energy( &
    species, mass_fractions, temperature, ok)
  if (.not. ok) error stop "Failed to evaluate final H2/O2 energy"
  relative_energy_error = abs(current_energy - target_internal_energy) / &
    max(1.0_dp, abs(target_internal_energy))

  write(*, '(a,i0)') "Completed steps: ", step
  write(*, '(a,i0)') "Output rows: ", output_count
  write(*, '(a,es24.16)') "Final time: ", time
  write(*, '(a,es24.16)') "Final temperature: ", temperature
  write(*, '(a,es24.16)') "Final pressure: ", pressure
  write(*, '(a,es24.16)') "Relative energy error: ", relative_energy_error
  write(*, '(a,1x,a)') "Output:", trim(config%output_file)

contains

  subroutine write_reactor_row()
    current_energy = reactor_specific_internal_energy( &
      species, mass_fractions, temperature, ok)
    if (.not. ok) error stop "Failed to evaluate H2/O2 reactor energy"
    pressure = mixture_pressure( &
      species, mass_fractions, density, temperature, ok)
    if (.not. ok) error stop "Failed to evaluate H2/O2 reactor pressure"
    call elementary_production_rates( &
      species, reactions, temperature, density, mass_fractions, &
      molar_production_rates, ok)
    if (.not. ok) error stop "Failed to evaluate H2/O2 production rates"
    relative_energy_error = abs(current_energy - target_internal_energy) / &
      max(1.0_dp, abs(target_internal_energy))
    closure_error = abs(sum(mass_fractions) - 1.0_dp)
    write(unit, '(*(es25.16e3,:,","))') &
      time, temperature, pressure, density, current_energy, &
      relative_energy_error, closure_error, &
      mass_fractions(h2o2_h2_index), mass_fractions(h2o2_h_index), &
      mass_fractions(h2o2_o_index), mass_fractions(h2o2_o2_index), &
      mass_fractions(h2o2_oh_index), mass_fractions(h2o2_h2o_index), &
      mass_fractions(h2o2_n2_index), &
      molar_production_rates(h2o2_h2_index), &
      molar_production_rates(h2o2_h_index), &
      molar_production_rates(h2o2_o_index), &
      molar_production_rates(h2o2_o2_index), &
      molar_production_rates(h2o2_oh_index), &
      molar_production_rates(h2o2_h2o_index), &
      molar_production_rates(h2o2_n2_index)
    output_count = output_count + 1
  end subroutine write_reactor_row

end program pelef0d_h2o2
