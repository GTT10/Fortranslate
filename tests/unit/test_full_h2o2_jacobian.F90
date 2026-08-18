program test_full_h2o2_jacobian
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density, mixture_mass_properties
  use pressure_dependent_kinetics_mod, only: pressure_dependent_reaction
  use h2o2_full_thermo_mod, only: full_nspecies, load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use h2o2_full_jacobian_mod, only: &
    h2o2_full_active_species, evaluate_h2o2_full_jacobian
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(pressure_dependent_reaction), allocatable :: reactions(:)
  real(dp) :: mole_fractions(full_nspecies), mass_fractions(full_nspecies)
  real(dp) :: derivative(full_nspecies)
  real(dp) :: jacobian(h2o2_full_active_species, h2o2_full_active_species)
  real(dp) :: density, temperature, pressure, target_energy
  real(dp) :: molecular_weight, gas_constant, cp, cv, gamma, enthalpy, entropy
  real(dp) :: maximum_magnitude
  logical :: ok

  call load_h2o2_full_thermo(species, ok)
  if (.not. ok) error stop "Failed to load full thermo"
  call load_h2o2_full_mechanism(reactions, ok)
  if (.not. ok) error stop "Failed to load full mechanism"

  mole_fractions = [2.0_dp, 1.0e-5_dp, 2.0e-6_dp, 1.0_dp, &
    3.0e-6_dp, 1.0e-4_dp, 1.0e-5_dp, 1.0e-6_dp, 0.1_dp, 3.0_dp]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  if (.not. ok) error stop "Mole-to-mass conversion failed"
  temperature = 1100.0_dp
  pressure = 202650.0_dp
  density = mixture_density(species, mass_fractions, pressure, temperature, ok)
  if (.not. ok) error stop "Density evaluation failed"
  call mixture_mass_properties( &
    species, mass_fractions, temperature, molecular_weight, gas_constant, &
    cp, cv, gamma, enthalpy, target_energy, entropy, ok)
  if (.not. ok) error stop "Energy evaluation failed"

  call evaluate_h2o2_full_jacobian( &
    species, reactions, density, target_energy, mass_fractions, temperature, &
    derivative, jacobian, temperature, ok)
  if (.not. ok) error stop "Generated Jacobian evaluation failed"
  maximum_magnitude = maxval(abs(jacobian))
  if (maximum_magnitude <= 1.0_dp) then
    error stop "Generated Jacobian did not capture stiff chemistry"
  end if
  if (any(jacobian /= jacobian)) error stop "Generated Jacobian contains NaN"
  if (maxval(abs(derivative)) <= 0.0_dp) error stop "Full RHS is identically zero"

  write(*, '(a,es24.16)') &
    "test_full_h2o2_jacobian: PASS, max|J|=", maximum_magnitude
end program test_full_h2o2_jacobian
