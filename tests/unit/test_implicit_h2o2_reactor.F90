program test_implicit_h2o2_reactor
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density, mixture_mass_properties
  use pressure_dependent_kinetics_mod, only: pressure_dependent_reaction
  use h2o2_full_thermo_mod, only: full_nspecies, load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use implicit_reactor_mod, only: advance_implicit_adaptive
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(pressure_dependent_reaction), allocatable :: reactions(:)
  real(dp) :: mole_fractions(full_nspecies), mass_fractions(full_nspecies)
  real(dp) :: density, temperature, target_energy, current_energy
  real(dp) :: molecular_weight, gas_constant, cp, cv, gamma, enthalpy, entropy
  real(dp) :: accepted_step, suggested_step
  integer :: newton_iterations, rejected_steps
  logical :: ok

  call load_h2o2_full_thermo(species, ok)
  if (.not. ok) error stop "Failed to load full thermo"
  call load_h2o2_full_mechanism(reactions, ok)
  if (.not. ok) error stop "Failed to load full mechanism"
  mole_fractions = [2.0_dp, 1.0e-6_dp, 1.0e-12_dp, 1.0_dp, &
    1.0e-12_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 3.0_dp]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  if (.not. ok) error stop "Composition conversion failed"
  temperature = 1000.0_dp
  density = mixture_density( &
    species, mass_fractions, 101325.0_dp, temperature, ok)
  if (.not. ok) error stop "Density evaluation failed"
  call mixture_mass_properties( &
    species, mass_fractions, temperature, molecular_weight, gas_constant, &
    cp, cv, gamma, enthalpy, target_energy, entropy, ok)
  if (.not. ok) error stop "Initial energy evaluation failed"

  call advance_implicit_adaptive( &
    species, reactions, density, target_energy, 1.0e-7_dp, 1.0e-14_dp, &
    1.0e-5_dp, 1.0e-6_dp, 1.0e-12_dp, 18, mass_fractions, temperature, &
    accepted_step, suggested_step, newton_iterations, rejected_steps, ok)
  if (.not. ok) error stop "Implicit adaptive step failed"
  if (accepted_step <= 0.0_dp .or. suggested_step <= 0.0_dp) then
    error stop "Implicit adaptive step returned invalid step sizes"
  end if
  if (minval(mass_fractions) < -1.0e-13_dp) then
    error stop "Implicit step generated negative composition"
  end if
  if (abs(sum(mass_fractions) - 1.0_dp) > 1.0e-12_dp) then
    error stop "Implicit step violated composition closure"
  end if
  call mixture_mass_properties( &
    species, mass_fractions, temperature, molecular_weight, gas_constant, &
    cp, cv, gamma, enthalpy, current_energy, entropy, ok)
  if (.not. ok) error stop "Final energy evaluation failed"
  if (abs(current_energy - target_energy) / max(1.0_dp, abs(target_energy)) &
      > 5.0e-10_dp) then
    error stop "Implicit step violated constant-volume energy"
  end if

  write(*, '(a)') "test_implicit_h2o2_reactor: PASS"
end program test_implicit_h2o2_reactor
