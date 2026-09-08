program test_sundials_constant_volume_reactor
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density, &
    mixture_mass_properties
  use elementary_kinetics_mod, only: elementary_reaction
  use sundials_constant_volume_reactor_mod, only: &
    cvode_reactor_context, cvode_reactor_statistics, &
    initialize_constant_volume_cvode, &
    advance_constant_volume_cvode, &
    get_constant_volume_cvode_statistics, finalize_constant_volume_cvode
  use fixture_mechanism_mod, only: &
    load_fixture_thermo_data, load_fixture_mechanism
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(cvode_reactor_context) :: context
  type(cvode_reactor_statistics) :: statistics, calibration_statistics
  real(dp) :: mole_fractions(2), mass_fractions(2), saved_state(2)
  real(dp) :: initial_state(2), density, temperature, saved_temperature
  real(dp) :: initial_temperature
  real(dp) :: molecular_weight, gas_constant, cp, cv, gamma
  real(dp) :: enthalpy, target_energy, final_energy, entropy
  character(len=1024) :: message
  logical :: ok
  integer :: cumulative_limit

  call load_fixture_thermo_data(species, ok)
  if (.not. ok) error stop "Could not load fixture thermodynamics"
  call load_fixture_mechanism(reactions, ok)
  if (.not. ok) error stop "Could not load fixture reactions"
  mole_fractions = [0.8_dp, 0.2_dp]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  if (.not. ok) error stop "Could not create fixture mass fractions"
  temperature = 1000.0_dp
  density = mixture_density( &
    species, mass_fractions, 101325.0_dp, temperature, ok)
  if (.not. ok) error stop "Could not create fixture density"
  call mixture_mass_properties( &
    species, mass_fractions, temperature, molecular_weight, gas_constant, &
    cp, cv, gamma, enthalpy, target_energy, entropy, ok)
  if (.not. ok) error stop "Could not create fixture energy"
  initial_state = mass_fractions
  initial_temperature = temperature

  call initialize_constant_volume_cvode( &
    context, species, reactions, density, target_energy, 0.0_dp, 1.0e-12_dp, &
    1.0e-16_dp, 2.5e-8_dp, 1.0e-7_dp, 1.0e-13_dp, 1, &
    mass_fractions, temperature, ok, message)
  if (.not. ok) error stop trim(message)
  saved_state = mass_fractions
  saved_temperature = temperature
  call advance_constant_volume_cvode( &
    context, 1.0e-7_dp, mass_fractions, temperature, ok, message)
  if (ok .or. index(message, "FCVode failed") == 0) then
    error stop "CVODE maximum-step failure was not reported"
  end if
  if (any(mass_fractions /= saved_state) .or. &
      temperature /= saved_temperature) then
    error stop "Failed CVODE advance modified the caller state"
  end if
  call finalize_constant_volume_cvode(context, ok, message)
  if (.not. ok) error stop trim(message)

  mass_fractions = initial_state
  temperature = initial_temperature
  call initialize_constant_volume_cvode( &
    context, species, reactions, density, target_energy, 0.0_dp, 1.0e-12_dp, &
    1.0e-16_dp, 2.5e-8_dp, 1.0e-7_dp, 1.0e-13_dp, 100000, &
    mass_fractions, temperature, ok, message)
  if (.not. ok) error stop trim(message)
  call advance_constant_volume_cvode( &
    context, 2.5e-8_dp, mass_fractions, temperature, ok, message)
  if (.not. ok) error stop trim(message)
  call get_constant_volume_cvode_statistics( &
    context, calibration_statistics, ok, message)
  if (.not. ok) error stop trim(message)
  cumulative_limit = int(calibration_statistics%internal_steps)
  if (cumulative_limit <= 0) then
    error stop "CVODE cumulative-step calibration did not advance"
  end if
  call finalize_constant_volume_cvode(context, ok, message)
  if (.not. ok) error stop trim(message)

  mass_fractions = initial_state
  temperature = initial_temperature
  call initialize_constant_volume_cvode( &
    context, species, reactions, density, target_energy, 0.0_dp, 1.0e-12_dp, &
    1.0e-16_dp, 2.5e-8_dp, 1.0e-7_dp, 1.0e-13_dp, cumulative_limit, &
    mass_fractions, temperature, ok, message)
  if (.not. ok) error stop trim(message)
  call advance_constant_volume_cvode( &
    context, 2.5e-8_dp, mass_fractions, temperature, ok, message)
  if (.not. ok) error stop trim(message)
  call get_constant_volume_cvode_statistics( &
    context, statistics, ok, message)
  if (.not. ok) error stop trim(message)
  if (statistics%internal_steps /= calibration_statistics%internal_steps) then
    error stop "CVODE cumulative-step calibration was not deterministic"
  end if
  saved_state = mass_fractions
  saved_temperature = temperature
  call advance_constant_volume_cvode( &
    context, 5.0e-8_dp, mass_fractions, temperature, ok, message)
  if (ok .or. index(message, "cumulative step limit reached") == 0) then
    error stop "CVODE cumulative step limit was not enforced before advance"
  end if
  if (any(mass_fractions /= saved_state) .or. &
      temperature /= saved_temperature) then
    error stop "Cumulative-limit failure modified the caller state"
  end if
  call get_constant_volume_cvode_statistics( &
    context, statistics, ok, message)
  if (.not. ok) error stop trim(message)
  if (statistics%internal_steps /= calibration_statistics%internal_steps) then
    error stop "Cumulative-limit rejection performed an internal step"
  end if
  call finalize_constant_volume_cvode(context, ok, message)
  if (.not. ok) error stop trim(message)

  mass_fractions = initial_state
  temperature = initial_temperature
  call initialize_constant_volume_cvode( &
    context, species, reactions, density, target_energy, 0.0_dp, 1.0e-12_dp, &
    1.0e-16_dp, 2.5e-8_dp, 1.0e-7_dp, 1.0e-13_dp, 100000, &
    mass_fractions, temperature, ok, message)
  if (.not. ok) error stop trim(message)
  call initialize_constant_volume_cvode( &
    context, species, reactions, density, target_energy, 0.0_dp, 1.0e-12_dp, &
    1.0e-16_dp, 2.5e-8_dp, 1.0e-7_dp, 1.0e-13_dp, 100000, &
    mass_fractions, temperature, ok, message)
  if (ok .or. index(message, "already assigned") == 0) then
    error stop "Assigned CVODE handle reinitialization was not rejected"
  end if

  call advance_constant_volume_cvode( &
    context, 2.5e-8_dp, mass_fractions, temperature, ok, message)
  if (.not. ok) error stop trim(message)
  saved_state = mass_fractions
  saved_temperature = temperature
  call advance_constant_volume_cvode( &
    context, 2.5e-8_dp, mass_fractions, temperature, ok, message)
  if (ok .or. index(message, "monotonically") == 0) then
    error stop "Nonmonotonic CVODE output time was not rejected"
  end if
  if (any(mass_fractions /= saved_state) .or. &
      temperature /= saved_temperature) then
    error stop "Rejected CVODE target time modified the caller state"
  end if

  call get_constant_volume_cvode_statistics( &
    context, statistics, ok, message)
  if (.not. ok) error stop trim(message)
  if (statistics%internal_steps <= 0 .or. &
      statistics%rhs_evaluations <= 0 .or. &
      statistics%jacobian_evaluations <= 0 .or. &
      statistics%nonlinear_iterations <= 0) then
    error stop "CVODE statistics did not record the integrated trajectory"
  end if
  if (minval(mass_fractions) < 0.0_dp .or. &
      abs(sum(mass_fractions) - 1.0_dp) > 1.0e-14_dp) then
    error stop "CVODE fixture state violated composition closure"
  end if
  call mixture_mass_properties( &
    species, mass_fractions, temperature, molecular_weight, gas_constant, &
    cp, cv, gamma, enthalpy, final_energy, entropy, ok)
  if (.not. ok) error stop "Could not evaluate the CVODE final energy"
  if (abs(final_energy - target_energy) / max(1.0_dp, abs(target_energy)) &
      > 5.0e-10_dp) then
    error stop "CVODE fixture state violated constant-volume energy"
  end if
  call finalize_constant_volume_cvode(context, ok, message)
  if (.not. ok) error stop trim(message)
  call finalize_constant_volume_cvode(context, ok, message)
  if (.not. ok) error stop "Idempotent CVODE finalization failed"

  write(*, '(a)') "test_sundials_constant_volume_reactor: PASS"
end program test_sundials_constant_volume_reactor
