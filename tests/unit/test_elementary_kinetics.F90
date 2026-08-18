program test_elementary_kinetics
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use elementary_kinetics_mod, only: &
    elementary_reaction, valid_elementary_reaction, arrhenius_rate_constant, &
    reaction_equilibrium_constant, elementary_production_rates
  use h2o2_elementary_mechanism_mod, only: &
    h2o2_nspecies, h2o2_nreactions, h2o2_h2_index, h2o2_h_index, &
    h2o2_o_index, h2o2_o2_index, h2o2_oh_index, h2o2_h2o_index, &
    h2o2_n2_index, load_h2o2_elementary_mechanism, &
    h2o2_elementary_production_rates
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp) :: mole_fractions(h2o2_nspecies), mass_fractions(h2o2_nspecies)
  real(dp) :: generic_rates(h2o2_nspecies), generated_rates(h2o2_nspecies)
  real(dp) :: temperature, density, rate_constant, equilibrium_constant
  real(dp) :: mass_source, hydrogen_source, oxygen_source
  logical :: ok
  integer :: i

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok .or. size(species) /= h2o2_nspecies) &
    error stop "H2/O2 thermo loading failed"
  call load_h2o2_elementary_mechanism(reactions, ok)
  if (.not. ok .or. size(reactions) /= h2o2_nreactions) &
    error stop "H2/O2 mechanism loading failed"
  do i = 1, size(reactions)
    if (.not. valid_elementary_reaction(reactions(i), size(species))) &
      error stop "invalid generated elementary reaction"
  end do

  temperature = 1200.0_dp
  rate_constant = arrhenius_rate_constant( &
    reactions(1)%forward_rate, temperature, ok)
  if (.not. ok) error stop "Arrhenius evaluation failed"
  if (abs(rate_constant - 38.7_dp * temperature**2.7_dp * &
      exp(-26191840.0_dp / (8314.46261815324_dp * temperature))) > &
      1.0e-12_dp * rate_constant) error stop "Arrhenius mismatch"

  equilibrium_constant = reaction_equilibrium_constant( &
    species, reactions(1), temperature, ok)
  if (.not. ok .or. equilibrium_constant <= 0.0_dp) &
    error stop "equilibrium-constant evaluation failed"

  mole_fractions = [ &
    0.295_dp, 1.0e-6_dp, 1.0e-6_dp, 0.1475_dp, 1.0e-6_dp, &
    0.0_dp, 0.557497_dp ]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  if (.not. ok) error stop "mole-to-mass conversion failed"
  density = 0.25_dp

  call elementary_production_rates( &
    species, reactions, temperature, density, mass_fractions, &
    generic_rates, ok)
  if (.not. ok) error stop "generic production-rate evaluation failed"
  call h2o2_elementary_production_rates( &
    species, reactions, temperature, density, mass_fractions, &
    generated_rates, ok)
  if (.not. ok) error stop "generated production-rate evaluation failed"
  if (maxval(abs(generic_rates - generated_rates)) > &
      1.0e-13_dp * max(1.0_dp, maxval(abs(generic_rates)))) &
    error stop "generated and generic kernels disagree"

  mass_source = 0.0_dp
  do i = 1, size(species)
    mass_source = mass_source + species(i)%molecular_weight * generic_rates(i)
  end do
  hydrogen_source = 2.0_dp * generic_rates(h2o2_h2_index) + &
    generic_rates(h2o2_h_index) + generic_rates(h2o2_oh_index) + &
    2.0_dp * generic_rates(h2o2_h2o_index)
  oxygen_source = generic_rates(h2o2_o_index) + &
    2.0_dp * generic_rates(h2o2_o2_index) + &
    generic_rates(h2o2_oh_index) + generic_rates(h2o2_h2o_index)

  if (abs(mass_source) > 1.0e-9_dp * &
      max(1.0_dp, maxval(abs(generic_rates)))) &
    error stop "elementary mechanism violates mass conservation"
  if (abs(hydrogen_source) > 1.0e-10_dp * &
      max(1.0_dp, maxval(abs(generic_rates)))) &
    error stop "elementary mechanism violates H conservation"
  if (abs(oxygen_source) > 1.0e-10_dp * &
      max(1.0_dp, maxval(abs(generic_rates)))) &
    error stop "elementary mechanism violates O conservation"
  if (abs(generic_rates(h2o2_n2_index)) > tiny(1.0_dp)) &
    error stop "inert N2 has a production rate"

  write(*, '(a,es24.16)') "reaction_1_rate_constant=", rate_constant
  write(*, '(a,es24.16)') "reaction_1_equilibrium_constant=", &
    equilibrium_constant
  write(*, '(a,es24.16)') "maximum_production_rate=", &
    maxval(abs(generic_rates))
  write(*, '(a)') "test_elementary_kinetics: PASS"
end program test_elementary_kinetics
