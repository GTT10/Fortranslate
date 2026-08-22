program pelef0d_h2o2_full
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
  use nasa7_thermo_mod, only: nasa7_species
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density, mixture_pressure, &
    mixture_mass_properties
  use elementary_kinetics_mod, only: &
    elementary_reaction, elementary_production_rates
  use h2o2_full_thermo_mod, only: full_nspecies, load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use constant_volume_reactor_mod, only: &
    advance_constant_volume_implicit_adaptive
  use simulation_config_h2o2_full_mod, only: &
    h2o2_full_reactor_config, read_h2o2_full_reactor_configuration
  implicit none

  type(h2o2_full_reactor_config) :: config
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp) :: mole_fractions(full_nspecies), mass_fractions(full_nspecies)
  real(dp) :: production_rates(full_nspecies)
  real(dp) :: density, temperature, target_energy, current_energy, pressure
  real(dp) :: molecular_weight, gas_constant, cp, cv, gamma, enthalpy, entropy
  real(dp) :: time, requested_step, accepted_step, suggested_step, next_output
  real(dp) :: relative_energy_error, closure_error, time_tolerance
  integer :: unit, io_status, steps, total_newton, total_rejected
  integer :: newton_iterations, rejected_steps
  character(len=1024) :: input_path, message
  logical :: ok

  if (command_argument_count() /= 1) then
    write(*, '(a)') "Usage: pelef0d_h2o2_full <input.nml>"
    error stop 2
  end if
  call get_command_argument(1, input_path)
  call read_h2o2_full_reactor_configuration( &
    trim(input_path), config, ok, message)
  if (.not. ok) then
    write(*, '(a)') trim(message)
    error stop 2
  end if
  call load_h2o2_full_thermo(species, ok)
  if (.not. ok) error stop "Failed to load full H2/O2 thermodynamics"
  call load_h2o2_full_mechanism(reactions, ok)
  if (.not. ok) error stop "Failed to load full H2/O2 mechanism"

  mole_fractions = [config%x_h2, config%x_h, config%x_o, config%x_o2, &
    config%x_oh, config%x_h2o, config%x_ho2, config%x_h2o2, &
    config%x_ar, config%x_n2]
  mole_fractions = mole_fractions / sum(mole_fractions)
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  if (.not. ok) error stop "Failed to convert full H2/O2 composition"
  temperature = config%initial_temperature
  density = mixture_density( &
    species, mass_fractions, config%initial_pressure, temperature, ok)
  if (.not. ok) error stop "Failed to evaluate initial full H2/O2 density"
  call mixture_mass_properties( &
    species, mass_fractions, temperature, molecular_weight, gas_constant, &
    cp, cv, gamma, enthalpy, target_energy, entropy, ok)
  if (.not. ok) error stop "Failed to evaluate initial full H2/O2 energy"

  open(newunit=unit, file=trim(config%output_file), status="replace", &
    action="write", iostat=io_status)
  if (io_status /= 0) error stop "Could not create full H2/O2 CSV"
  write(unit, '(a)') &
    "time,temperature,pressure,density,specific_internal_energy," // &
    "relative_energy_error,closure_error,Y_H2,Y_H,Y_O,Y_O2,Y_OH,Y_H2O," // &
    "Y_HO2,Y_H2O2,Y_AR,Y_N2,wdot_H2,wdot_H,wdot_O,wdot_O2,wdot_OH," // &
    "wdot_H2O,wdot_HO2,wdot_H2O2,wdot_AR,wdot_N2"

  write(*, '(a)') "PeleF " // pelef_version // " full H2/O2"
  write(*, '(a,1x,a)') "Input:", trim(input_path)
  write(*, '(a,es24.16)') "Constant density: ", density

  time = 0.0_dp
  requested_step = config%initial_time_step
  next_output = 0.0_dp
  steps = 0
  total_newton = 0
  total_rejected = 0
  time_tolerance = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, config%final_time)
  call write_row()
  next_output = min(config%output_interval, config%final_time)

  do while (time < config%final_time - time_tolerance)
    if (steps >= config%maximum_steps) error stop "Full H2/O2 step limit reached"
    requested_step = min(requested_step, config%maximum_time_step)
    requested_step = min(requested_step, config%final_time - time)
    requested_step = min(requested_step, next_output - time)
    if (requested_step <= time_tolerance) then
      time = next_output
      call write_row()
      if (next_output >= config%final_time - time_tolerance) exit
      next_output = min(next_output + config%output_interval, config%final_time)
      cycle
    end if
    call advance_constant_volume_implicit_adaptive( &
      species, reactions, density, target_energy, requested_step, &
      config%relative_tolerance, config%absolute_tolerance, mass_fractions, &
      temperature, accepted_step, suggested_step, newton_iterations, &
      rejected_steps, ok)
    if (.not. ok .or. accepted_step < config%minimum_time_step) then
      error stop "Implicit full H2/O2 step failed"
    end if
    time = time + accepted_step
    requested_step = min(config%maximum_time_step, suggested_step)
    steps = steps + 1
    total_newton = total_newton + newton_iterations
    total_rejected = total_rejected + rejected_steps
    if (time >= next_output - time_tolerance) then
      time = next_output
      call write_row()
      if (next_output >= config%final_time - time_tolerance) exit
      next_output = min(next_output + config%output_interval, config%final_time)
    end if
  end do
  close(unit)

  call mixture_mass_properties( &
    species, mass_fractions, temperature, molecular_weight, gas_constant, &
    cp, cv, gamma, enthalpy, current_energy, entropy, ok)
  if (.not. ok) error stop "Failed to evaluate final full H2/O2 energy"
  relative_energy_error = abs(current_energy - target_energy) / &
    max(1.0_dp, abs(target_energy))
  pressure = mixture_pressure(species, mass_fractions, density, temperature, ok)
  if (.not. ok) error stop "Failed to evaluate final full H2/O2 pressure"
  write(*, '(a,i0)') "Accepted implicit steps: ", steps
  write(*, '(a,i0)') "Newton iterations: ", total_newton
  write(*, '(a,i0)') "Rejected trial steps: ", total_rejected
  write(*, '(a,es24.16)') "Final temperature: ", temperature
  write(*, '(a,es24.16)') "Final pressure: ", pressure
  write(*, '(a,es24.16)') "Relative energy error: ", relative_energy_error
  write(*, '(a,1x,a)') "Output:", trim(config%output_file)

contains

  subroutine write_row()
    call mixture_mass_properties( &
      species, mass_fractions, temperature, molecular_weight, gas_constant, &
      cp, cv, gamma, enthalpy, current_energy, entropy, ok)
    if (.not. ok) error stop "Failed to evaluate full H2/O2 row energy"
    pressure = mixture_pressure(species, mass_fractions, density, temperature, ok)
    if (.not. ok) error stop "Failed to evaluate full H2/O2 row pressure"
    call elementary_production_rates( &
      species, reactions, temperature, density, mass_fractions, &
      production_rates, ok)
    if (.not. ok) error stop "Failed to evaluate full H2/O2 production rates"
    relative_energy_error = abs(current_energy - target_energy) / &
      max(1.0_dp, abs(target_energy))
    closure_error = abs(sum(mass_fractions) - 1.0_dp)
    write(unit, '(*(es25.16e3,:,","))') time, temperature, pressure, density, &
      current_energy, relative_energy_error, closure_error, mass_fractions, &
      production_rates
  end subroutine write_row

end program pelef0d_h2o2_full
