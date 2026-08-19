module h2o2_full_mechanism_mod
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: &
    elementary_reaction, elementary_production_rates, &
    elementary_mass_fraction_jacobian, reaction_kind_elementary, &
    reaction_kind_three_body, reaction_kind_falloff
  implicit none
  private

  integer, parameter, public :: h2o2_full_nspecies = 10
  integer, parameter, public :: h2o2_full_nreactions = 29
  integer, parameter, public :: h2o2_full_h2_index = 1
  integer, parameter, public :: h2o2_full_h_index = 2
  integer, parameter, public :: h2o2_full_o_index = 3
  integer, parameter, public :: h2o2_full_o2_index = 4
  integer, parameter, public :: h2o2_full_oh_index = 5
  integer, parameter, public :: h2o2_full_h2o_index = 6
  integer, parameter, public :: h2o2_full_ho2_index = 7
  integer, parameter, public :: h2o2_full_h2o2_index = 8
  integer, parameter, public :: h2o2_full_ar_index = 9
  integer, parameter, public :: h2o2_full_n2_index = 10

  public :: load_h2o2_full_mechanism
  public :: h2o2_full_production_rates
  public :: h2o2_full_mass_fraction_jacobian

contains

  subroutine load_h2o2_full_mechanism(reactions, ok)
    type(elementary_reaction), allocatable, intent(out) :: reactions(:)
    logical, intent(out) :: ok

    allocate(reactions(29))

    reactions(1)%equation = "2 O + M <=> O2 + M"
    reactions(1)%kind = reaction_kind_three_body
    allocate(reactions(1)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(1)%product_stoich(h2o2_full_nspecies))
    reactions(1)%reactant_stoich = 0.0_dp
    reactions(1)%product_stoich = 0.0_dp
    reactions(1)%reactant_stoich(3) = 2.000000000000e+00_dp
    reactions(1)%product_stoich(4) = 1.000000000000e+00_dp
    reactions(1)%forward_rate%pre_exponential = 1.200000000000e+11_dp
    reactions(1)%forward_rate%temperature_exponent = -1.000000000000e+00_dp
    reactions(1)%forward_rate%activation_energy = 0.000000000000e+00_dp
    allocate(reactions(1)%third_body_efficiencies(h2o2_full_nspecies))
    reactions(1)%third_body_efficiencies = 1.000000000000e+00_dp
    reactions(1)%third_body_efficiencies(1) = 2.400000000000e+00_dp
    reactions(1)%third_body_efficiencies(6) = 1.540000000000e+01_dp
    reactions(1)%third_body_efficiencies(9) = 8.300000000000e-01_dp
    reactions(1)%reversible = .true.

    reactions(2)%equation = "O + H + M <=> OH + M"
    reactions(2)%kind = reaction_kind_three_body
    allocate(reactions(2)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(2)%product_stoich(h2o2_full_nspecies))
    reactions(2)%reactant_stoich = 0.0_dp
    reactions(2)%product_stoich = 0.0_dp
    reactions(2)%reactant_stoich(3) = 1.000000000000e+00_dp
    reactions(2)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(2)%product_stoich(5) = 1.000000000000e+00_dp
    reactions(2)%forward_rate%pre_exponential = 5.000000000000e+11_dp
    reactions(2)%forward_rate%temperature_exponent = -1.000000000000e+00_dp
    reactions(2)%forward_rate%activation_energy = 0.000000000000e+00_dp
    allocate(reactions(2)%third_body_efficiencies(h2o2_full_nspecies))
    reactions(2)%third_body_efficiencies = 1.000000000000e+00_dp
    reactions(2)%third_body_efficiencies(1) = 2.000000000000e+00_dp
    reactions(2)%third_body_efficiencies(6) = 6.000000000000e+00_dp
    reactions(2)%third_body_efficiencies(9) = 7.000000000000e-01_dp
    reactions(2)%reversible = .true.

    reactions(3)%equation = "O + H2 <=> H + OH"
    reactions(3)%kind = reaction_kind_elementary
    allocate(reactions(3)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(3)%product_stoich(h2o2_full_nspecies))
    reactions(3)%reactant_stoich = 0.0_dp
    reactions(3)%product_stoich = 0.0_dp
    reactions(3)%reactant_stoich(3) = 1.000000000000e+00_dp
    reactions(3)%reactant_stoich(1) = 1.000000000000e+00_dp
    reactions(3)%product_stoich(2) = 1.000000000000e+00_dp
    reactions(3)%product_stoich(5) = 1.000000000000e+00_dp
    reactions(3)%forward_rate%pre_exponential = 3.870000000000e+01_dp
    reactions(3)%forward_rate%temperature_exponent = 2.700000000000e+00_dp
    reactions(3)%forward_rate%activation_energy = 2.619184000000e+07_dp
    reactions(3)%reversible = .true.

    reactions(4)%equation = "O + HO2 <=> OH + O2"
    reactions(4)%kind = reaction_kind_elementary
    allocate(reactions(4)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(4)%product_stoich(h2o2_full_nspecies))
    reactions(4)%reactant_stoich = 0.0_dp
    reactions(4)%product_stoich = 0.0_dp
    reactions(4)%reactant_stoich(3) = 1.000000000000e+00_dp
    reactions(4)%reactant_stoich(7) = 1.000000000000e+00_dp
    reactions(4)%product_stoich(5) = 1.000000000000e+00_dp
    reactions(4)%product_stoich(4) = 1.000000000000e+00_dp
    reactions(4)%forward_rate%pre_exponential = 2.000000000000e+10_dp
    reactions(4)%forward_rate%temperature_exponent = 0.000000000000e+00_dp
    reactions(4)%forward_rate%activation_energy = 0.000000000000e+00_dp
    reactions(4)%reversible = .true.

    reactions(5)%equation = "O + H2O2 <=> OH + HO2"
    reactions(5)%kind = reaction_kind_elementary
    allocate(reactions(5)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(5)%product_stoich(h2o2_full_nspecies))
    reactions(5)%reactant_stoich = 0.0_dp
    reactions(5)%product_stoich = 0.0_dp
    reactions(5)%reactant_stoich(3) = 1.000000000000e+00_dp
    reactions(5)%reactant_stoich(8) = 1.000000000000e+00_dp
    reactions(5)%product_stoich(5) = 1.000000000000e+00_dp
    reactions(5)%product_stoich(7) = 1.000000000000e+00_dp
    reactions(5)%forward_rate%pre_exponential = 9.630000000000e+03_dp
    reactions(5)%forward_rate%temperature_exponent = 2.000000000000e+00_dp
    reactions(5)%forward_rate%activation_energy = 1.673600000000e+07_dp
    reactions(5)%reversible = .true.

    reactions(6)%equation = "H + O2 + M <=> HO2 + M"
    reactions(6)%kind = reaction_kind_three_body
    allocate(reactions(6)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(6)%product_stoich(h2o2_full_nspecies))
    reactions(6)%reactant_stoich = 0.0_dp
    reactions(6)%product_stoich = 0.0_dp
    reactions(6)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(6)%reactant_stoich(4) = 1.000000000000e+00_dp
    reactions(6)%product_stoich(7) = 1.000000000000e+00_dp
    reactions(6)%forward_rate%pre_exponential = 2.800000000000e+12_dp
    reactions(6)%forward_rate%temperature_exponent = -8.600000000000e-01_dp
    reactions(6)%forward_rate%activation_energy = 0.000000000000e+00_dp
    allocate(reactions(6)%third_body_efficiencies(h2o2_full_nspecies))
    reactions(6)%third_body_efficiencies = 1.000000000000e+00_dp
    reactions(6)%third_body_efficiencies(4) = 0.000000000000e+00_dp
    reactions(6)%third_body_efficiencies(6) = 0.000000000000e+00_dp
    reactions(6)%third_body_efficiencies(10) = 0.000000000000e+00_dp
    reactions(6)%third_body_efficiencies(9) = 0.000000000000e+00_dp
    reactions(6)%reversible = .true.

    reactions(7)%equation = "H + 2 O2 <=> HO2 + O2"
    reactions(7)%kind = reaction_kind_elementary
    allocate(reactions(7)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(7)%product_stoich(h2o2_full_nspecies))
    reactions(7)%reactant_stoich = 0.0_dp
    reactions(7)%product_stoich = 0.0_dp
    reactions(7)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(7)%reactant_stoich(4) = 2.000000000000e+00_dp
    reactions(7)%product_stoich(7) = 1.000000000000e+00_dp
    reactions(7)%product_stoich(4) = 1.000000000000e+00_dp
    reactions(7)%forward_rate%pre_exponential = 2.080000000000e+13_dp
    reactions(7)%forward_rate%temperature_exponent = -1.240000000000e+00_dp
    reactions(7)%forward_rate%activation_energy = 0.000000000000e+00_dp
    reactions(7)%reversible = .true.

    reactions(8)%equation = "H + O2 + H2O <=> HO2 + H2O"
    reactions(8)%kind = reaction_kind_elementary
    allocate(reactions(8)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(8)%product_stoich(h2o2_full_nspecies))
    reactions(8)%reactant_stoich = 0.0_dp
    reactions(8)%product_stoich = 0.0_dp
    reactions(8)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(8)%reactant_stoich(4) = 1.000000000000e+00_dp
    reactions(8)%reactant_stoich(6) = 1.000000000000e+00_dp
    reactions(8)%product_stoich(7) = 1.000000000000e+00_dp
    reactions(8)%product_stoich(6) = 1.000000000000e+00_dp
    reactions(8)%forward_rate%pre_exponential = 1.126000000000e+13_dp
    reactions(8)%forward_rate%temperature_exponent = -7.600000000000e-01_dp
    reactions(8)%forward_rate%activation_energy = 0.000000000000e+00_dp
    reactions(8)%reversible = .true.

    reactions(9)%equation = "H + O2 + N2 <=> HO2 + N2"
    reactions(9)%kind = reaction_kind_elementary
    allocate(reactions(9)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(9)%product_stoich(h2o2_full_nspecies))
    reactions(9)%reactant_stoich = 0.0_dp
    reactions(9)%product_stoich = 0.0_dp
    reactions(9)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(9)%reactant_stoich(4) = 1.000000000000e+00_dp
    reactions(9)%reactant_stoich(10) = 1.000000000000e+00_dp
    reactions(9)%product_stoich(7) = 1.000000000000e+00_dp
    reactions(9)%product_stoich(10) = 1.000000000000e+00_dp
    reactions(9)%forward_rate%pre_exponential = 2.600000000000e+13_dp
    reactions(9)%forward_rate%temperature_exponent = -1.240000000000e+00_dp
    reactions(9)%forward_rate%activation_energy = 0.000000000000e+00_dp
    reactions(9)%reversible = .true.

    reactions(10)%equation = "H + O2 + AR <=> HO2 + AR"
    reactions(10)%kind = reaction_kind_elementary
    allocate(reactions(10)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(10)%product_stoich(h2o2_full_nspecies))
    reactions(10)%reactant_stoich = 0.0_dp
    reactions(10)%product_stoich = 0.0_dp
    reactions(10)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(10)%reactant_stoich(4) = 1.000000000000e+00_dp
    reactions(10)%reactant_stoich(9) = 1.000000000000e+00_dp
    reactions(10)%product_stoich(7) = 1.000000000000e+00_dp
    reactions(10)%product_stoich(9) = 1.000000000000e+00_dp
    reactions(10)%forward_rate%pre_exponential = 7.000000000000e+11_dp
    reactions(10)%forward_rate%temperature_exponent = -8.000000000000e-01_dp
    reactions(10)%forward_rate%activation_energy = 0.000000000000e+00_dp
    reactions(10)%reversible = .true.

    reactions(11)%equation = "H + O2 <=> O + OH"
    reactions(11)%kind = reaction_kind_elementary
    allocate(reactions(11)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(11)%product_stoich(h2o2_full_nspecies))
    reactions(11)%reactant_stoich = 0.0_dp
    reactions(11)%product_stoich = 0.0_dp
    reactions(11)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(11)%reactant_stoich(4) = 1.000000000000e+00_dp
    reactions(11)%product_stoich(3) = 1.000000000000e+00_dp
    reactions(11)%product_stoich(5) = 1.000000000000e+00_dp
    reactions(11)%forward_rate%pre_exponential = 2.650000000000e+13_dp
    reactions(11)%forward_rate%temperature_exponent = -6.707000000000e-01_dp
    reactions(11)%forward_rate%activation_energy = 7.129954400000e+07_dp
    reactions(11)%reversible = .true.

    reactions(12)%equation = "2 H + M <=> H2 + M"
    reactions(12)%kind = reaction_kind_three_body
    allocate(reactions(12)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(12)%product_stoich(h2o2_full_nspecies))
    reactions(12)%reactant_stoich = 0.0_dp
    reactions(12)%product_stoich = 0.0_dp
    reactions(12)%reactant_stoich(2) = 2.000000000000e+00_dp
    reactions(12)%product_stoich(1) = 1.000000000000e+00_dp
    reactions(12)%forward_rate%pre_exponential = 1.000000000000e+12_dp
    reactions(12)%forward_rate%temperature_exponent = -1.000000000000e+00_dp
    reactions(12)%forward_rate%activation_energy = 0.000000000000e+00_dp
    allocate(reactions(12)%third_body_efficiencies(h2o2_full_nspecies))
    reactions(12)%third_body_efficiencies = 1.000000000000e+00_dp
    reactions(12)%third_body_efficiencies(1) = 0.000000000000e+00_dp
    reactions(12)%third_body_efficiencies(6) = 0.000000000000e+00_dp
    reactions(12)%third_body_efficiencies(9) = 6.300000000000e-01_dp
    reactions(12)%reversible = .true.

    reactions(13)%equation = "2 H + H2 <=> 2 H2"
    reactions(13)%kind = reaction_kind_elementary
    allocate(reactions(13)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(13)%product_stoich(h2o2_full_nspecies))
    reactions(13)%reactant_stoich = 0.0_dp
    reactions(13)%product_stoich = 0.0_dp
    reactions(13)%reactant_stoich(2) = 2.000000000000e+00_dp
    reactions(13)%reactant_stoich(1) = 1.000000000000e+00_dp
    reactions(13)%product_stoich(1) = 2.000000000000e+00_dp
    reactions(13)%forward_rate%pre_exponential = 9.000000000000e+10_dp
    reactions(13)%forward_rate%temperature_exponent = -6.000000000000e-01_dp
    reactions(13)%forward_rate%activation_energy = 0.000000000000e+00_dp
    reactions(13)%reversible = .true.

    reactions(14)%equation = "2 H + H2O <=> H2 + H2O"
    reactions(14)%kind = reaction_kind_elementary
    allocate(reactions(14)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(14)%product_stoich(h2o2_full_nspecies))
    reactions(14)%reactant_stoich = 0.0_dp
    reactions(14)%product_stoich = 0.0_dp
    reactions(14)%reactant_stoich(2) = 2.000000000000e+00_dp
    reactions(14)%reactant_stoich(6) = 1.000000000000e+00_dp
    reactions(14)%product_stoich(1) = 1.000000000000e+00_dp
    reactions(14)%product_stoich(6) = 1.000000000000e+00_dp
    reactions(14)%forward_rate%pre_exponential = 6.000000000000e+13_dp
    reactions(14)%forward_rate%temperature_exponent = -1.250000000000e+00_dp
    reactions(14)%forward_rate%activation_energy = 0.000000000000e+00_dp
    reactions(14)%reversible = .true.

    reactions(15)%equation = "H + OH + M <=> H2O + M"
    reactions(15)%kind = reaction_kind_three_body
    allocate(reactions(15)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(15)%product_stoich(h2o2_full_nspecies))
    reactions(15)%reactant_stoich = 0.0_dp
    reactions(15)%product_stoich = 0.0_dp
    reactions(15)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(15)%reactant_stoich(5) = 1.000000000000e+00_dp
    reactions(15)%product_stoich(6) = 1.000000000000e+00_dp
    reactions(15)%forward_rate%pre_exponential = 2.200000000000e+16_dp
    reactions(15)%forward_rate%temperature_exponent = -2.000000000000e+00_dp
    reactions(15)%forward_rate%activation_energy = 0.000000000000e+00_dp
    allocate(reactions(15)%third_body_efficiencies(h2o2_full_nspecies))
    reactions(15)%third_body_efficiencies = 1.000000000000e+00_dp
    reactions(15)%third_body_efficiencies(1) = 7.300000000000e-01_dp
    reactions(15)%third_body_efficiencies(6) = 3.650000000000e+00_dp
    reactions(15)%third_body_efficiencies(9) = 3.800000000000e-01_dp
    reactions(15)%reversible = .true.

    reactions(16)%equation = "H + HO2 <=> O + H2O"
    reactions(16)%kind = reaction_kind_elementary
    allocate(reactions(16)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(16)%product_stoich(h2o2_full_nspecies))
    reactions(16)%reactant_stoich = 0.0_dp
    reactions(16)%product_stoich = 0.0_dp
    reactions(16)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(16)%reactant_stoich(7) = 1.000000000000e+00_dp
    reactions(16)%product_stoich(3) = 1.000000000000e+00_dp
    reactions(16)%product_stoich(6) = 1.000000000000e+00_dp
    reactions(16)%forward_rate%pre_exponential = 3.970000000000e+09_dp
    reactions(16)%forward_rate%temperature_exponent = 0.000000000000e+00_dp
    reactions(16)%forward_rate%activation_energy = 2.807464000000e+06_dp
    reactions(16)%reversible = .true.

    reactions(17)%equation = "H + HO2 <=> O2 + H2"
    reactions(17)%kind = reaction_kind_elementary
    allocate(reactions(17)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(17)%product_stoich(h2o2_full_nspecies))
    reactions(17)%reactant_stoich = 0.0_dp
    reactions(17)%product_stoich = 0.0_dp
    reactions(17)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(17)%reactant_stoich(7) = 1.000000000000e+00_dp
    reactions(17)%product_stoich(4) = 1.000000000000e+00_dp
    reactions(17)%product_stoich(1) = 1.000000000000e+00_dp
    reactions(17)%forward_rate%pre_exponential = 4.480000000000e+10_dp
    reactions(17)%forward_rate%temperature_exponent = 0.000000000000e+00_dp
    reactions(17)%forward_rate%activation_energy = 4.468512000000e+06_dp
    reactions(17)%reversible = .true.

    reactions(18)%equation = "H + HO2 <=> 2 OH"
    reactions(18)%kind = reaction_kind_elementary
    allocate(reactions(18)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(18)%product_stoich(h2o2_full_nspecies))
    reactions(18)%reactant_stoich = 0.0_dp
    reactions(18)%product_stoich = 0.0_dp
    reactions(18)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(18)%reactant_stoich(7) = 1.000000000000e+00_dp
    reactions(18)%product_stoich(5) = 2.000000000000e+00_dp
    reactions(18)%forward_rate%pre_exponential = 8.400000000000e+10_dp
    reactions(18)%forward_rate%temperature_exponent = 0.000000000000e+00_dp
    reactions(18)%forward_rate%activation_energy = 2.656840000000e+06_dp
    reactions(18)%reversible = .true.

    reactions(19)%equation = "H + H2O2 <=> HO2 + H2"
    reactions(19)%kind = reaction_kind_elementary
    allocate(reactions(19)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(19)%product_stoich(h2o2_full_nspecies))
    reactions(19)%reactant_stoich = 0.0_dp
    reactions(19)%product_stoich = 0.0_dp
    reactions(19)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(19)%reactant_stoich(8) = 1.000000000000e+00_dp
    reactions(19)%product_stoich(7) = 1.000000000000e+00_dp
    reactions(19)%product_stoich(1) = 1.000000000000e+00_dp
    reactions(19)%forward_rate%pre_exponential = 1.210000000000e+04_dp
    reactions(19)%forward_rate%temperature_exponent = 2.000000000000e+00_dp
    reactions(19)%forward_rate%activation_energy = 2.175680000000e+07_dp
    reactions(19)%reversible = .true.

    reactions(20)%equation = "H + H2O2 <=> OH + H2O"
    reactions(20)%kind = reaction_kind_elementary
    allocate(reactions(20)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(20)%product_stoich(h2o2_full_nspecies))
    reactions(20)%reactant_stoich = 0.0_dp
    reactions(20)%product_stoich = 0.0_dp
    reactions(20)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(20)%reactant_stoich(8) = 1.000000000000e+00_dp
    reactions(20)%product_stoich(5) = 1.000000000000e+00_dp
    reactions(20)%product_stoich(6) = 1.000000000000e+00_dp
    reactions(20)%forward_rate%pre_exponential = 1.000000000000e+10_dp
    reactions(20)%forward_rate%temperature_exponent = 0.000000000000e+00_dp
    reactions(20)%forward_rate%activation_energy = 1.506240000000e+07_dp
    reactions(20)%reversible = .true.

    reactions(21)%equation = "OH + H2 <=> H + H2O"
    reactions(21)%kind = reaction_kind_elementary
    allocate(reactions(21)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(21)%product_stoich(h2o2_full_nspecies))
    reactions(21)%reactant_stoich = 0.0_dp
    reactions(21)%product_stoich = 0.0_dp
    reactions(21)%reactant_stoich(5) = 1.000000000000e+00_dp
    reactions(21)%reactant_stoich(1) = 1.000000000000e+00_dp
    reactions(21)%product_stoich(2) = 1.000000000000e+00_dp
    reactions(21)%product_stoich(6) = 1.000000000000e+00_dp
    reactions(21)%forward_rate%pre_exponential = 2.160000000000e+05_dp
    reactions(21)%forward_rate%temperature_exponent = 1.510000000000e+00_dp
    reactions(21)%forward_rate%activation_energy = 1.435112000000e+07_dp
    reactions(21)%reversible = .true.

    reactions(22)%equation = "2 OH (+M) <=> H2O2 (+M)"
    reactions(22)%kind = reaction_kind_falloff
    allocate(reactions(22)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(22)%product_stoich(h2o2_full_nspecies))
    reactions(22)%reactant_stoich = 0.0_dp
    reactions(22)%product_stoich = 0.0_dp
    reactions(22)%reactant_stoich(5) = 2.000000000000e+00_dp
    reactions(22)%product_stoich(8) = 1.000000000000e+00_dp
    reactions(22)%low_pressure_rate%pre_exponential = 2.300000000000e+12_dp
    reactions(22)%low_pressure_rate%temperature_exponent = -9.000000000000e-01_dp
    reactions(22)%low_pressure_rate%activation_energy = -7.112800000000e+06_dp
    reactions(22)%high_pressure_rate%pre_exponential = 7.400000000000e+10_dp
    reactions(22)%high_pressure_rate%temperature_exponent = -3.700000000000e-01_dp
    reactions(22)%high_pressure_rate%activation_energy = 0.000000000000e+00_dp
    allocate(reactions(22)%third_body_efficiencies(h2o2_full_nspecies))
    reactions(22)%third_body_efficiencies = 1.000000000000e+00_dp
    reactions(22)%third_body_efficiencies(1) = 2.000000000000e+00_dp
    reactions(22)%third_body_efficiencies(6) = 6.000000000000e+00_dp
    reactions(22)%third_body_efficiencies(9) = 7.000000000000e-01_dp
    reactions(22)%troe%enabled = .true.
    reactions(22)%troe%alpha = 7.346000000000e-01_dp
    reactions(22)%troe%temperature_3 = 9.400000000000e+01_dp
    reactions(22)%troe%temperature_1 = 1.756000000000e+03_dp
    reactions(22)%troe%temperature_2 = 5.182000000000e+03_dp
    reactions(22)%reversible = .true.

    reactions(23)%equation = "2 OH <=> O + H2O"
    reactions(23)%kind = reaction_kind_elementary
    allocate(reactions(23)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(23)%product_stoich(h2o2_full_nspecies))
    reactions(23)%reactant_stoich = 0.0_dp
    reactions(23)%product_stoich = 0.0_dp
    reactions(23)%reactant_stoich(5) = 2.000000000000e+00_dp
    reactions(23)%product_stoich(3) = 1.000000000000e+00_dp
    reactions(23)%product_stoich(6) = 1.000000000000e+00_dp
    reactions(23)%forward_rate%pre_exponential = 3.570000000000e+01_dp
    reactions(23)%forward_rate%temperature_exponent = 2.400000000000e+00_dp
    reactions(23)%forward_rate%activation_energy = -8.828240000000e+06_dp
    reactions(23)%reversible = .true.

    reactions(24)%equation = "OH + HO2 <=> O2 + H2O [1]"
    reactions(24)%kind = reaction_kind_elementary
    allocate(reactions(24)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(24)%product_stoich(h2o2_full_nspecies))
    reactions(24)%reactant_stoich = 0.0_dp
    reactions(24)%product_stoich = 0.0_dp
    reactions(24)%reactant_stoich(5) = 1.000000000000e+00_dp
    reactions(24)%reactant_stoich(7) = 1.000000000000e+00_dp
    reactions(24)%product_stoich(4) = 1.000000000000e+00_dp
    reactions(24)%product_stoich(6) = 1.000000000000e+00_dp
    reactions(24)%forward_rate%pre_exponential = 1.450000000000e+10_dp
    reactions(24)%forward_rate%temperature_exponent = 0.000000000000e+00_dp
    reactions(24)%forward_rate%activation_energy = -2.092000000000e+06_dp
    reactions(24)%reversible = .true.

    reactions(25)%equation = "OH + H2O2 <=> HO2 + H2O [1]"
    reactions(25)%kind = reaction_kind_elementary
    allocate(reactions(25)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(25)%product_stoich(h2o2_full_nspecies))
    reactions(25)%reactant_stoich = 0.0_dp
    reactions(25)%product_stoich = 0.0_dp
    reactions(25)%reactant_stoich(5) = 1.000000000000e+00_dp
    reactions(25)%reactant_stoich(8) = 1.000000000000e+00_dp
    reactions(25)%product_stoich(7) = 1.000000000000e+00_dp
    reactions(25)%product_stoich(6) = 1.000000000000e+00_dp
    reactions(25)%forward_rate%pre_exponential = 2.000000000000e+09_dp
    reactions(25)%forward_rate%temperature_exponent = 0.000000000000e+00_dp
    reactions(25)%forward_rate%activation_energy = 1.786568000000e+06_dp
    reactions(25)%reversible = .true.

    reactions(26)%equation = "OH + H2O2 <=> HO2 + H2O [2]"
    reactions(26)%kind = reaction_kind_elementary
    allocate(reactions(26)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(26)%product_stoich(h2o2_full_nspecies))
    reactions(26)%reactant_stoich = 0.0_dp
    reactions(26)%product_stoich = 0.0_dp
    reactions(26)%reactant_stoich(5) = 1.000000000000e+00_dp
    reactions(26)%reactant_stoich(8) = 1.000000000000e+00_dp
    reactions(26)%product_stoich(7) = 1.000000000000e+00_dp
    reactions(26)%product_stoich(6) = 1.000000000000e+00_dp
    reactions(26)%forward_rate%pre_exponential = 1.700000000000e+15_dp
    reactions(26)%forward_rate%temperature_exponent = 0.000000000000e+00_dp
    reactions(26)%forward_rate%activation_energy = 1.230514400000e+08_dp
    reactions(26)%reversible = .true.

    reactions(27)%equation = "2 HO2 <=> O2 + H2O2 [1]"
    reactions(27)%kind = reaction_kind_elementary
    allocate(reactions(27)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(27)%product_stoich(h2o2_full_nspecies))
    reactions(27)%reactant_stoich = 0.0_dp
    reactions(27)%product_stoich = 0.0_dp
    reactions(27)%reactant_stoich(7) = 2.000000000000e+00_dp
    reactions(27)%product_stoich(4) = 1.000000000000e+00_dp
    reactions(27)%product_stoich(8) = 1.000000000000e+00_dp
    reactions(27)%forward_rate%pre_exponential = 1.300000000000e+08_dp
    reactions(27)%forward_rate%temperature_exponent = 0.000000000000e+00_dp
    reactions(27)%forward_rate%activation_energy = -6.819920000000e+06_dp
    reactions(27)%reversible = .true.

    reactions(28)%equation = "2 HO2 <=> O2 + H2O2 [2]"
    reactions(28)%kind = reaction_kind_elementary
    allocate(reactions(28)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(28)%product_stoich(h2o2_full_nspecies))
    reactions(28)%reactant_stoich = 0.0_dp
    reactions(28)%product_stoich = 0.0_dp
    reactions(28)%reactant_stoich(7) = 2.000000000000e+00_dp
    reactions(28)%product_stoich(4) = 1.000000000000e+00_dp
    reactions(28)%product_stoich(8) = 1.000000000000e+00_dp
    reactions(28)%forward_rate%pre_exponential = 4.200000000000e+11_dp
    reactions(28)%forward_rate%temperature_exponent = 0.000000000000e+00_dp
    reactions(28)%forward_rate%activation_energy = 5.020800000000e+07_dp
    reactions(28)%reversible = .true.

    reactions(29)%equation = "OH + HO2 <=> O2 + H2O [2]"
    reactions(29)%kind = reaction_kind_elementary
    allocate(reactions(29)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(29)%product_stoich(h2o2_full_nspecies))
    reactions(29)%reactant_stoich = 0.0_dp
    reactions(29)%product_stoich = 0.0_dp
    reactions(29)%reactant_stoich(5) = 1.000000000000e+00_dp
    reactions(29)%reactant_stoich(7) = 1.000000000000e+00_dp
    reactions(29)%product_stoich(4) = 1.000000000000e+00_dp
    reactions(29)%product_stoich(6) = 1.000000000000e+00_dp
    reactions(29)%forward_rate%pre_exponential = 5.000000000000e+12_dp
    reactions(29)%forward_rate%temperature_exponent = 0.000000000000e+00_dp
    reactions(29)%forward_rate%activation_energy = 7.250872000000e+07_dp
    reactions(29)%reversible = .true.

    ok = .true.
  end subroutine load_h2o2_full_mechanism

  subroutine h2o2_full_production_rates( &
      species, reactions, temperature, density, mass_fractions, &
      molar_production_rates, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: temperature, density, mass_fractions(:)
    real(dp), intent(out) :: molar_production_rates(:)
    logical, intent(out) :: ok

    ok = size(species) == h2o2_full_nspecies .and. &
      size(reactions) == h2o2_full_nreactions
    if (.not. ok) then
      molar_production_rates = 0.0_dp
      return
    end if
    call elementary_production_rates( &
      species, reactions, temperature, density, mass_fractions, &
      molar_production_rates, ok)
  end subroutine h2o2_full_production_rates

  subroutine h2o2_full_mass_fraction_jacobian( &
      species, reactions, temperature, density, mass_fractions, &
      jacobian, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: temperature, density, mass_fractions(:)
    real(dp), intent(out) :: jacobian(:, :)
    logical, intent(out) :: ok

    ok = size(species) == h2o2_full_nspecies .and. &
      size(reactions) == h2o2_full_nreactions
    if (.not. ok) then
      jacobian = 0.0_dp
      return
    end if
    call elementary_mass_fraction_jacobian( &
      species, reactions, temperature, density, mass_fractions, &
      jacobian, ok)
  end subroutine h2o2_full_mass_fraction_jacobian

end module h2o2_full_mechanism_mod
