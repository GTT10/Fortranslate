program test_implicit_h2o2_reactor
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density, mixture_mass_properties
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_full_thermo_mod, only: full_nspecies, load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use constant_volume_reactor_mod, only: &
    advance_constant_volume_implicit_adaptive
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp) :: x(full_nspecies), y(full_nspecies), density, temperature
  real(dp) :: target_energy, final_energy, molecular_weight, gas_constant
  real(dp) :: cp, cv, gamma, enthalpy, entropy, accepted, next_step
  integer :: newton_iterations, rejected
  logical :: ok

  call load_h2o2_full_thermo(species, ok)
  if (.not. ok) error stop "Failed to load full thermo"
  call load_h2o2_full_mechanism(reactions, ok)
  if (.not. ok) error stop "Failed to load full mechanism"
  x = [2.0_dp, 1.0e-6_dp, 1.0e-12_dp, 1.0_dp, 1.0e-12_dp, &
    0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 3.0_dp]
  x = x / sum(x)
  call mass_fractions_from_mole_fractions(species, x, y, ok)
  if (.not. ok) error stop "Composition conversion failed"
  temperature = 1000.0_dp
  density = mixture_density(species, y, 101325.0_dp, temperature, ok)
  if (.not. ok) error stop "Density evaluation failed"
  call mixture_mass_properties(species, y, temperature, molecular_weight, &
    gas_constant, cp, cv, gamma, enthalpy, target_energy, entropy, ok)
  if (.not. ok) error stop "Initial energy evaluation failed"
  call advance_constant_volume_implicit_adaptive( &
    species, reactions, density, target_energy, 1.0e-7_dp, 1.0e-6_dp, &
    1.0e-12_dp, y, temperature, accepted, next_step, newton_iterations, &
    rejected, ok)
  if (.not. ok) error stop "Implicit adaptive step failed"
  if (accepted <= 0.0_dp .or. next_step <= 0.0_dp) then
    error stop "Implicit adaptive step returned invalid time steps"
  end if
  if (minval(y) < -1.0e-13_dp .or. abs(sum(y) - 1.0_dp) > 1.0e-12_dp) then
    error stop "Implicit step violated composition bounds"
  end if
  call mixture_mass_properties(species, y, temperature, molecular_weight, &
    gas_constant, cp, cv, gamma, enthalpy, final_energy, entropy, ok)
  if (.not. ok) error stop "Final energy evaluation failed"
  if (abs(final_energy - target_energy) / max(1.0_dp, abs(target_energy)) &
      > 5.0e-10_dp) then
    error stop "Implicit step violated constant-volume energy"
  end if
  write(*, '(a)') "test_implicit_h2o2_reactor: PASS"
end program test_implicit_h2o2_reactor
