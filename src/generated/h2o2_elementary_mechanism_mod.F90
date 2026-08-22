module h2o2_elementary_mechanism_mod
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: &
    elementary_reaction, elementary_production_rates, &
    elementary_mass_fraction_jacobian, reaction_kind_elementary, &
    reaction_kind_three_body, reaction_kind_falloff
  implicit none
  private

  integer, parameter, public :: h2o2_nspecies = 7
  integer, parameter, public :: h2o2_nreactions = 4
  integer, parameter, public :: h2o2_h2_index = 1
  integer, parameter, public :: h2o2_h_index = 2
  integer, parameter, public :: h2o2_o_index = 3
  integer, parameter, public :: h2o2_o2_index = 4
  integer, parameter, public :: h2o2_oh_index = 5
  integer, parameter, public :: h2o2_h2o_index = 6
  integer, parameter, public :: h2o2_n2_index = 7

  public :: load_h2o2_elementary_mechanism
  public :: h2o2_elementary_production_rates
  public :: h2o2_elementary_mass_fraction_jacobian

contains

  subroutine load_h2o2_elementary_mechanism(reactions, ok)
    type(elementary_reaction), allocatable, intent(out) :: reactions(:)
    logical, intent(out) :: ok

    allocate(reactions(4))

    reactions(1)%equation = "O + H2 <=> H + OH"
    reactions(1)%kind = reaction_kind_elementary
    allocate(reactions(1)%reactant_stoich(h2o2_nspecies))
    allocate(reactions(1)%product_stoich(h2o2_nspecies))
    reactions(1)%reactant_stoich = 0.0_dp
    reactions(1)%product_stoich = 0.0_dp
    reactions(1)%reactant_stoich(1) = 1.000000000000e+00_dp
    reactions(1)%reactant_stoich(3) = 1.000000000000e+00_dp
    reactions(1)%product_stoich(2) = 1.000000000000e+00_dp
    reactions(1)%product_stoich(5) = 1.000000000000e+00_dp
    reactions(1)%forward_rate%pre_exponential = 3.870000000000e+01_dp
    reactions(1)%forward_rate%temperature_exponent = 2.700000000000e+00_dp
    reactions(1)%forward_rate%activation_energy = 2.619184000000e+07_dp
    reactions(1)%reversible = .true.

    reactions(2)%equation = "H + O2 <=> O + OH"
    reactions(2)%kind = reaction_kind_elementary
    allocate(reactions(2)%reactant_stoich(h2o2_nspecies))
    allocate(reactions(2)%product_stoich(h2o2_nspecies))
    reactions(2)%reactant_stoich = 0.0_dp
    reactions(2)%product_stoich = 0.0_dp
    reactions(2)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(2)%reactant_stoich(4) = 1.000000000000e+00_dp
    reactions(2)%product_stoich(3) = 1.000000000000e+00_dp
    reactions(2)%product_stoich(5) = 1.000000000000e+00_dp
    reactions(2)%forward_rate%pre_exponential = 2.650000000000e+13_dp
    reactions(2)%forward_rate%temperature_exponent = -6.707000000000e-01_dp
    reactions(2)%forward_rate%activation_energy = 7.129954400000e+07_dp
    reactions(2)%reversible = .true.

    reactions(3)%equation = "OH + H2 <=> H + H2O"
    reactions(3)%kind = reaction_kind_elementary
    allocate(reactions(3)%reactant_stoich(h2o2_nspecies))
    allocate(reactions(3)%product_stoich(h2o2_nspecies))
    reactions(3)%reactant_stoich = 0.0_dp
    reactions(3)%product_stoich = 0.0_dp
    reactions(3)%reactant_stoich(1) = 1.000000000000e+00_dp
    reactions(3)%reactant_stoich(5) = 1.000000000000e+00_dp
    reactions(3)%product_stoich(2) = 1.000000000000e+00_dp
    reactions(3)%product_stoich(6) = 1.000000000000e+00_dp
    reactions(3)%forward_rate%pre_exponential = 2.160000000000e+05_dp
    reactions(3)%forward_rate%temperature_exponent = 1.510000000000e+00_dp
    reactions(3)%forward_rate%activation_energy = 1.435112000000e+07_dp
    reactions(3)%reversible = .true.

    reactions(4)%equation = "2 OH <=> O + H2O"
    reactions(4)%kind = reaction_kind_elementary
    allocate(reactions(4)%reactant_stoich(h2o2_nspecies))
    allocate(reactions(4)%product_stoich(h2o2_nspecies))
    reactions(4)%reactant_stoich = 0.0_dp
    reactions(4)%product_stoich = 0.0_dp
    reactions(4)%reactant_stoich(5) = 2.000000000000e+00_dp
    reactions(4)%product_stoich(3) = 1.000000000000e+00_dp
    reactions(4)%product_stoich(6) = 1.000000000000e+00_dp
    reactions(4)%forward_rate%pre_exponential = 3.570000000000e+01_dp
    reactions(4)%forward_rate%temperature_exponent = 2.400000000000e+00_dp
    reactions(4)%forward_rate%activation_energy = -8.828240000000e+06_dp
    reactions(4)%reversible = .true.

    ok = .true.
  end subroutine load_h2o2_elementary_mechanism

  subroutine h2o2_elementary_production_rates( &
      species, reactions, temperature, density, mass_fractions, &
      molar_production_rates, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: temperature, density, mass_fractions(:)
    real(dp), intent(out) :: molar_production_rates(:)
    logical, intent(out) :: ok

    ok = size(species) == h2o2_nspecies .and. &
      size(reactions) == h2o2_nreactions
    if (.not. ok) then
      molar_production_rates = 0.0_dp
      return
    end if
    call elementary_production_rates( &
      species, reactions, temperature, density, mass_fractions, &
      molar_production_rates, ok)
  end subroutine h2o2_elementary_production_rates

  subroutine h2o2_elementary_mass_fraction_jacobian( &
      species, reactions, temperature, density, mass_fractions, &
      jacobian, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: temperature, density, mass_fractions(:)
    real(dp), intent(out) :: jacobian(:, :)
    logical, intent(out) :: ok

    ok = size(species) == h2o2_nspecies .and. &
      size(reactions) == h2o2_nreactions
    if (.not. ok) then
      jacobian = 0.0_dp
      return
    end if
    call elementary_mass_fraction_jacobian( &
      species, reactions, temperature, density, mass_fractions, &
      jacobian, ok)
  end subroutine h2o2_elementary_mass_fraction_jacobian

end module h2o2_elementary_mechanism_mod
