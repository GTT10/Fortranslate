module h2o2_full_mechanism_mod
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species, valid_nasa7_species
  use elementary_kinetics_mod, only: &
    elementary_reaction, elementary_production_rates, &
    elementary_mass_fraction_jacobian, reaction_kind_elementary, &
    reaction_kind_three_body, reaction_kind_falloff
  implicit none
  private

  integer, parameter, public :: h2o2_full_nspecies = 10
  integer, parameter, public :: h2o2_full_nreactions = 29
  character(len=*), parameter, public :: h2o2_full_chemistry_integrator = &
    "implicit"
  character(len=*), parameter, public :: h2o2_full_source_sha256 = &
    "0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4"
  character(len=*), parameter, public :: h2o2_full_source_phase = &
    "ohmech"
  character(len=*), parameter, public :: h2o2_full_source_cantera_version = &
    "2.5.0"
  character(len=*), parameter, public :: h2o2_full_runtime_cantera_version = &
    "3.2.0"
  character(len=*), parameter, public :: h2o2_full_runtime_cantera_git_commit = &
    "4a8358e"
  integer, parameter, public :: h2o2_full_source_indices(29) = [ &
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, &
    11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
    21, 22, 23, 24, 25, 26, 27, 28, 29 ]
  logical, parameter, public :: h2o2_full_duplicate_reactions(29) = [ &
    .false., .false., .false., .false., .false., .false., .false., &
    .false., .false., .false., .false., .false., .false., .false., &
    .false., .false., .false., .false., .false., .false., .false., &
    .false., .false., .true., .true., .true., .true., .true., &
    .true. ]
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
  public :: load_h2o2_full_thermo_data
  public :: load_h2o2_full_transport_data

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
    reactions(1)%third_body_efficiencies(9) = 8.300000000000e-01_dp
    reactions(1)%third_body_efficiencies(1) = 2.400000000000e+00_dp
    reactions(1)%third_body_efficiencies(6) = 1.540000000000e+01_dp
    reactions(1)%reversible = .true.

    reactions(2)%equation = "H + O + M <=> OH + M"
    reactions(2)%kind = reaction_kind_three_body
    allocate(reactions(2)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(2)%product_stoich(h2o2_full_nspecies))
    reactions(2)%reactant_stoich = 0.0_dp
    reactions(2)%product_stoich = 0.0_dp
    reactions(2)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(2)%reactant_stoich(3) = 1.000000000000e+00_dp
    reactions(2)%product_stoich(5) = 1.000000000000e+00_dp
    reactions(2)%forward_rate%pre_exponential = 5.000000000000e+11_dp
    reactions(2)%forward_rate%temperature_exponent = -1.000000000000e+00_dp
    reactions(2)%forward_rate%activation_energy = 0.000000000000e+00_dp
    allocate(reactions(2)%third_body_efficiencies(h2o2_full_nspecies))
    reactions(2)%third_body_efficiencies = 1.000000000000e+00_dp
    reactions(2)%third_body_efficiencies(9) = 7.000000000000e-01_dp
    reactions(2)%third_body_efficiencies(1) = 2.000000000000e+00_dp
    reactions(2)%third_body_efficiencies(6) = 6.000000000000e+00_dp
    reactions(2)%reversible = .true.

    reactions(3)%equation = "H2 + O <=> H + OH"
    reactions(3)%kind = reaction_kind_elementary
    allocate(reactions(3)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(3)%product_stoich(h2o2_full_nspecies))
    reactions(3)%reactant_stoich = 0.0_dp
    reactions(3)%product_stoich = 0.0_dp
    reactions(3)%reactant_stoich(1) = 1.000000000000e+00_dp
    reactions(3)%reactant_stoich(3) = 1.000000000000e+00_dp
    reactions(3)%product_stoich(2) = 1.000000000000e+00_dp
    reactions(3)%product_stoich(5) = 1.000000000000e+00_dp
    reactions(3)%forward_rate%pre_exponential = 3.870000000000e+01_dp
    reactions(3)%forward_rate%temperature_exponent = 2.700000000000e+00_dp
    reactions(3)%forward_rate%activation_energy = 2.619184000000e+07_dp
    reactions(3)%reversible = .true.

    reactions(4)%equation = "HO2 + O <=> O2 + OH"
    reactions(4)%kind = reaction_kind_elementary
    allocate(reactions(4)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(4)%product_stoich(h2o2_full_nspecies))
    reactions(4)%reactant_stoich = 0.0_dp
    reactions(4)%product_stoich = 0.0_dp
    reactions(4)%reactant_stoich(7) = 1.000000000000e+00_dp
    reactions(4)%reactant_stoich(3) = 1.000000000000e+00_dp
    reactions(4)%product_stoich(4) = 1.000000000000e+00_dp
    reactions(4)%product_stoich(5) = 1.000000000000e+00_dp
    reactions(4)%forward_rate%pre_exponential = 2.000000000000e+10_dp
    reactions(4)%forward_rate%temperature_exponent = 0.000000000000e+00_dp
    reactions(4)%forward_rate%activation_energy = 0.000000000000e+00_dp
    reactions(4)%reversible = .true.

    reactions(5)%equation = "H2O2 + O <=> HO2 + OH"
    reactions(5)%kind = reaction_kind_elementary
    allocate(reactions(5)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(5)%product_stoich(h2o2_full_nspecies))
    reactions(5)%reactant_stoich = 0.0_dp
    reactions(5)%product_stoich = 0.0_dp
    reactions(5)%reactant_stoich(8) = 1.000000000000e+00_dp
    reactions(5)%reactant_stoich(3) = 1.000000000000e+00_dp
    reactions(5)%product_stoich(7) = 1.000000000000e+00_dp
    reactions(5)%product_stoich(5) = 1.000000000000e+00_dp
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
    reactions(6)%third_body_efficiencies(9) = 0.000000000000e+00_dp
    reactions(6)%third_body_efficiencies(6) = 0.000000000000e+00_dp
    reactions(6)%third_body_efficiencies(10) = 0.000000000000e+00_dp
    reactions(6)%third_body_efficiencies(4) = 0.000000000000e+00_dp
    reactions(6)%reversible = .true.

    reactions(7)%equation = "H + O2 + O2 <=> HO2 + O2"
    reactions(7)%kind = reaction_kind_three_body
    allocate(reactions(7)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(7)%product_stoich(h2o2_full_nspecies))
    reactions(7)%reactant_stoich = 0.0_dp
    reactions(7)%product_stoich = 0.0_dp
    reactions(7)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(7)%reactant_stoich(4) = 1.000000000000e+00_dp
    reactions(7)%product_stoich(7) = 1.000000000000e+00_dp
    reactions(7)%forward_rate%pre_exponential = 2.080000000000e+13_dp
    reactions(7)%forward_rate%temperature_exponent = -1.240000000000e+00_dp
    reactions(7)%forward_rate%activation_energy = 0.000000000000e+00_dp
    allocate(reactions(7)%third_body_efficiencies(h2o2_full_nspecies))
    reactions(7)%third_body_efficiencies = 0.000000000000e+00_dp
    reactions(7)%third_body_efficiencies(4) = 1.000000000000e+00_dp
    reactions(7)%reversible = .true.

    reactions(8)%equation = "H + O2 + H2O <=> HO2 + H2O"
    reactions(8)%kind = reaction_kind_three_body
    allocate(reactions(8)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(8)%product_stoich(h2o2_full_nspecies))
    reactions(8)%reactant_stoich = 0.0_dp
    reactions(8)%product_stoich = 0.0_dp
    reactions(8)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(8)%reactant_stoich(4) = 1.000000000000e+00_dp
    reactions(8)%product_stoich(7) = 1.000000000000e+00_dp
    reactions(8)%forward_rate%pre_exponential = 1.126000000000e+13_dp
    reactions(8)%forward_rate%temperature_exponent = -7.600000000000e-01_dp
    reactions(8)%forward_rate%activation_energy = 0.000000000000e+00_dp
    allocate(reactions(8)%third_body_efficiencies(h2o2_full_nspecies))
    reactions(8)%third_body_efficiencies = 0.000000000000e+00_dp
    reactions(8)%third_body_efficiencies(6) = 1.000000000000e+00_dp
    reactions(8)%reversible = .true.

    reactions(9)%equation = "H + O2 + N2 <=> HO2 + N2"
    reactions(9)%kind = reaction_kind_three_body
    allocate(reactions(9)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(9)%product_stoich(h2o2_full_nspecies))
    reactions(9)%reactant_stoich = 0.0_dp
    reactions(9)%product_stoich = 0.0_dp
    reactions(9)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(9)%reactant_stoich(4) = 1.000000000000e+00_dp
    reactions(9)%product_stoich(7) = 1.000000000000e+00_dp
    reactions(9)%forward_rate%pre_exponential = 2.600000000000e+13_dp
    reactions(9)%forward_rate%temperature_exponent = -1.240000000000e+00_dp
    reactions(9)%forward_rate%activation_energy = 0.000000000000e+00_dp
    allocate(reactions(9)%third_body_efficiencies(h2o2_full_nspecies))
    reactions(9)%third_body_efficiencies = 0.000000000000e+00_dp
    reactions(9)%third_body_efficiencies(10) = 1.000000000000e+00_dp
    reactions(9)%reversible = .true.

    reactions(10)%equation = "H + O2 + AR <=> HO2 + AR"
    reactions(10)%kind = reaction_kind_three_body
    allocate(reactions(10)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(10)%product_stoich(h2o2_full_nspecies))
    reactions(10)%reactant_stoich = 0.0_dp
    reactions(10)%product_stoich = 0.0_dp
    reactions(10)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(10)%reactant_stoich(4) = 1.000000000000e+00_dp
    reactions(10)%product_stoich(7) = 1.000000000000e+00_dp
    reactions(10)%forward_rate%pre_exponential = 7.000000000000e+11_dp
    reactions(10)%forward_rate%temperature_exponent = -8.000000000000e-01_dp
    reactions(10)%forward_rate%activation_energy = 0.000000000000e+00_dp
    allocate(reactions(10)%third_body_efficiencies(h2o2_full_nspecies))
    reactions(10)%third_body_efficiencies = 0.000000000000e+00_dp
    reactions(10)%third_body_efficiencies(9) = 1.000000000000e+00_dp
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
    reactions(12)%third_body_efficiencies(9) = 6.300000000000e-01_dp
    reactions(12)%third_body_efficiencies(1) = 0.000000000000e+00_dp
    reactions(12)%third_body_efficiencies(6) = 0.000000000000e+00_dp
    reactions(12)%reversible = .true.

    reactions(13)%equation = "2 H + H2 <=> H2 + H2"
    reactions(13)%kind = reaction_kind_three_body
    allocate(reactions(13)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(13)%product_stoich(h2o2_full_nspecies))
    reactions(13)%reactant_stoich = 0.0_dp
    reactions(13)%product_stoich = 0.0_dp
    reactions(13)%reactant_stoich(2) = 2.000000000000e+00_dp
    reactions(13)%product_stoich(1) = 1.000000000000e+00_dp
    reactions(13)%forward_rate%pre_exponential = 9.000000000000e+10_dp
    reactions(13)%forward_rate%temperature_exponent = -6.000000000000e-01_dp
    reactions(13)%forward_rate%activation_energy = 0.000000000000e+00_dp
    allocate(reactions(13)%third_body_efficiencies(h2o2_full_nspecies))
    reactions(13)%third_body_efficiencies = 0.000000000000e+00_dp
    reactions(13)%third_body_efficiencies(1) = 1.000000000000e+00_dp
    reactions(13)%reversible = .true.

    reactions(14)%equation = "2 H + H2O <=> H2 + H2O"
    reactions(14)%kind = reaction_kind_three_body
    allocate(reactions(14)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(14)%product_stoich(h2o2_full_nspecies))
    reactions(14)%reactant_stoich = 0.0_dp
    reactions(14)%product_stoich = 0.0_dp
    reactions(14)%reactant_stoich(2) = 2.000000000000e+00_dp
    reactions(14)%product_stoich(1) = 1.000000000000e+00_dp
    reactions(14)%forward_rate%pre_exponential = 6.000000000000e+13_dp
    reactions(14)%forward_rate%temperature_exponent = -1.250000000000e+00_dp
    reactions(14)%forward_rate%activation_energy = 0.000000000000e+00_dp
    allocate(reactions(14)%third_body_efficiencies(h2o2_full_nspecies))
    reactions(14)%third_body_efficiencies = 0.000000000000e+00_dp
    reactions(14)%third_body_efficiencies(6) = 1.000000000000e+00_dp
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
    reactions(15)%third_body_efficiencies(9) = 3.800000000000e-01_dp
    reactions(15)%third_body_efficiencies(1) = 7.300000000000e-01_dp
    reactions(15)%third_body_efficiencies(6) = 3.650000000000e+00_dp
    reactions(15)%reversible = .true.

    reactions(16)%equation = "H + HO2 <=> H2O + O"
    reactions(16)%kind = reaction_kind_elementary
    allocate(reactions(16)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(16)%product_stoich(h2o2_full_nspecies))
    reactions(16)%reactant_stoich = 0.0_dp
    reactions(16)%product_stoich = 0.0_dp
    reactions(16)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(16)%reactant_stoich(7) = 1.000000000000e+00_dp
    reactions(16)%product_stoich(6) = 1.000000000000e+00_dp
    reactions(16)%product_stoich(3) = 1.000000000000e+00_dp
    reactions(16)%forward_rate%pre_exponential = 3.970000000000e+09_dp
    reactions(16)%forward_rate%temperature_exponent = 0.000000000000e+00_dp
    reactions(16)%forward_rate%activation_energy = 2.807464000000e+06_dp
    reactions(16)%reversible = .true.

    reactions(17)%equation = "H + HO2 <=> H2 + O2"
    reactions(17)%kind = reaction_kind_elementary
    allocate(reactions(17)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(17)%product_stoich(h2o2_full_nspecies))
    reactions(17)%reactant_stoich = 0.0_dp
    reactions(17)%product_stoich = 0.0_dp
    reactions(17)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(17)%reactant_stoich(7) = 1.000000000000e+00_dp
    reactions(17)%product_stoich(1) = 1.000000000000e+00_dp
    reactions(17)%product_stoich(4) = 1.000000000000e+00_dp
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

    reactions(19)%equation = "H + H2O2 <=> H2 + HO2"
    reactions(19)%kind = reaction_kind_elementary
    allocate(reactions(19)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(19)%product_stoich(h2o2_full_nspecies))
    reactions(19)%reactant_stoich = 0.0_dp
    reactions(19)%product_stoich = 0.0_dp
    reactions(19)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(19)%reactant_stoich(8) = 1.000000000000e+00_dp
    reactions(19)%product_stoich(1) = 1.000000000000e+00_dp
    reactions(19)%product_stoich(7) = 1.000000000000e+00_dp
    reactions(19)%forward_rate%pre_exponential = 1.210000000000e+04_dp
    reactions(19)%forward_rate%temperature_exponent = 2.000000000000e+00_dp
    reactions(19)%forward_rate%activation_energy = 2.175680000000e+07_dp
    reactions(19)%reversible = .true.

    reactions(20)%equation = "H + H2O2 <=> H2O + OH"
    reactions(20)%kind = reaction_kind_elementary
    allocate(reactions(20)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(20)%product_stoich(h2o2_full_nspecies))
    reactions(20)%reactant_stoich = 0.0_dp
    reactions(20)%product_stoich = 0.0_dp
    reactions(20)%reactant_stoich(2) = 1.000000000000e+00_dp
    reactions(20)%reactant_stoich(8) = 1.000000000000e+00_dp
    reactions(20)%product_stoich(6) = 1.000000000000e+00_dp
    reactions(20)%product_stoich(5) = 1.000000000000e+00_dp
    reactions(20)%forward_rate%pre_exponential = 1.000000000000e+10_dp
    reactions(20)%forward_rate%temperature_exponent = 0.000000000000e+00_dp
    reactions(20)%forward_rate%activation_energy = 1.506240000000e+07_dp
    reactions(20)%reversible = .true.

    reactions(21)%equation = "H2 + OH <=> H + H2O"
    reactions(21)%kind = reaction_kind_elementary
    allocate(reactions(21)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(21)%product_stoich(h2o2_full_nspecies))
    reactions(21)%reactant_stoich = 0.0_dp
    reactions(21)%product_stoich = 0.0_dp
    reactions(21)%reactant_stoich(1) = 1.000000000000e+00_dp
    reactions(21)%reactant_stoich(5) = 1.000000000000e+00_dp
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
    reactions(22)%third_body_efficiencies(9) = 7.000000000000e-01_dp
    reactions(22)%third_body_efficiencies(1) = 2.000000000000e+00_dp
    reactions(22)%third_body_efficiencies(6) = 6.000000000000e+00_dp
    reactions(22)%troe%enabled = .true.
    reactions(22)%troe%alpha = 7.346000000000e-01_dp
    reactions(22)%troe%temperature_3 = 9.400000000000e+01_dp
    reactions(22)%troe%temperature_1 = 1.756000000000e+03_dp
    reactions(22)%troe%temperature_2 = 5.182000000000e+03_dp
    reactions(22)%reversible = .true.

    reactions(23)%equation = "2 OH <=> H2O + O"
    reactions(23)%kind = reaction_kind_elementary
    allocate(reactions(23)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(23)%product_stoich(h2o2_full_nspecies))
    reactions(23)%reactant_stoich = 0.0_dp
    reactions(23)%product_stoich = 0.0_dp
    reactions(23)%reactant_stoich(5) = 2.000000000000e+00_dp
    reactions(23)%product_stoich(6) = 1.000000000000e+00_dp
    reactions(23)%product_stoich(3) = 1.000000000000e+00_dp
    reactions(23)%forward_rate%pre_exponential = 3.570000000000e+01_dp
    reactions(23)%forward_rate%temperature_exponent = 2.400000000000e+00_dp
    reactions(23)%forward_rate%activation_energy = -8.828240000000e+06_dp
    reactions(23)%reversible = .true.

    reactions(24)%equation = "HO2 + OH <=> H2O + O2"
    reactions(24)%kind = reaction_kind_elementary
    allocate(reactions(24)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(24)%product_stoich(h2o2_full_nspecies))
    reactions(24)%reactant_stoich = 0.0_dp
    reactions(24)%product_stoich = 0.0_dp
    reactions(24)%reactant_stoich(7) = 1.000000000000e+00_dp
    reactions(24)%reactant_stoich(5) = 1.000000000000e+00_dp
    reactions(24)%product_stoich(6) = 1.000000000000e+00_dp
    reactions(24)%product_stoich(4) = 1.000000000000e+00_dp
    reactions(24)%forward_rate%pre_exponential = 1.450000000000e+10_dp
    reactions(24)%forward_rate%temperature_exponent = 0.000000000000e+00_dp
    reactions(24)%forward_rate%activation_energy = -2.092000000000e+06_dp
    reactions(24)%reversible = .true.

    reactions(25)%equation = "H2O2 + OH <=> H2O + HO2"
    reactions(25)%kind = reaction_kind_elementary
    allocate(reactions(25)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(25)%product_stoich(h2o2_full_nspecies))
    reactions(25)%reactant_stoich = 0.0_dp
    reactions(25)%product_stoich = 0.0_dp
    reactions(25)%reactant_stoich(8) = 1.000000000000e+00_dp
    reactions(25)%reactant_stoich(5) = 1.000000000000e+00_dp
    reactions(25)%product_stoich(6) = 1.000000000000e+00_dp
    reactions(25)%product_stoich(7) = 1.000000000000e+00_dp
    reactions(25)%forward_rate%pre_exponential = 2.000000000000e+09_dp
    reactions(25)%forward_rate%temperature_exponent = 0.000000000000e+00_dp
    reactions(25)%forward_rate%activation_energy = 1.786568000000e+06_dp
    reactions(25)%reversible = .true.

    reactions(26)%equation = "H2O2 + OH <=> H2O + HO2"
    reactions(26)%kind = reaction_kind_elementary
    allocate(reactions(26)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(26)%product_stoich(h2o2_full_nspecies))
    reactions(26)%reactant_stoich = 0.0_dp
    reactions(26)%product_stoich = 0.0_dp
    reactions(26)%reactant_stoich(8) = 1.000000000000e+00_dp
    reactions(26)%reactant_stoich(5) = 1.000000000000e+00_dp
    reactions(26)%product_stoich(6) = 1.000000000000e+00_dp
    reactions(26)%product_stoich(7) = 1.000000000000e+00_dp
    reactions(26)%forward_rate%pre_exponential = 1.700000000000e+15_dp
    reactions(26)%forward_rate%temperature_exponent = 0.000000000000e+00_dp
    reactions(26)%forward_rate%activation_energy = 1.230514400000e+08_dp
    reactions(26)%reversible = .true.

    reactions(27)%equation = "2 HO2 <=> H2O2 + O2"
    reactions(27)%kind = reaction_kind_elementary
    allocate(reactions(27)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(27)%product_stoich(h2o2_full_nspecies))
    reactions(27)%reactant_stoich = 0.0_dp
    reactions(27)%product_stoich = 0.0_dp
    reactions(27)%reactant_stoich(7) = 2.000000000000e+00_dp
    reactions(27)%product_stoich(8) = 1.000000000000e+00_dp
    reactions(27)%product_stoich(4) = 1.000000000000e+00_dp
    reactions(27)%forward_rate%pre_exponential = 1.300000000000e+08_dp
    reactions(27)%forward_rate%temperature_exponent = 0.000000000000e+00_dp
    reactions(27)%forward_rate%activation_energy = -6.819920000000e+06_dp
    reactions(27)%reversible = .true.

    reactions(28)%equation = "2 HO2 <=> H2O2 + O2"
    reactions(28)%kind = reaction_kind_elementary
    allocate(reactions(28)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(28)%product_stoich(h2o2_full_nspecies))
    reactions(28)%reactant_stoich = 0.0_dp
    reactions(28)%product_stoich = 0.0_dp
    reactions(28)%reactant_stoich(7) = 2.000000000000e+00_dp
    reactions(28)%product_stoich(8) = 1.000000000000e+00_dp
    reactions(28)%product_stoich(4) = 1.000000000000e+00_dp
    reactions(28)%forward_rate%pre_exponential = 4.200000000000e+11_dp
    reactions(28)%forward_rate%temperature_exponent = 0.000000000000e+00_dp
    reactions(28)%forward_rate%activation_energy = 5.020800000000e+07_dp
    reactions(28)%reversible = .true.

    reactions(29)%equation = "HO2 + OH <=> H2O + O2"
    reactions(29)%kind = reaction_kind_elementary
    allocate(reactions(29)%reactant_stoich(h2o2_full_nspecies))
    allocate(reactions(29)%product_stoich(h2o2_full_nspecies))
    reactions(29)%reactant_stoich = 0.0_dp
    reactions(29)%product_stoich = 0.0_dp
    reactions(29)%reactant_stoich(7) = 1.000000000000e+00_dp
    reactions(29)%reactant_stoich(5) = 1.000000000000e+00_dp
    reactions(29)%product_stoich(6) = 1.000000000000e+00_dp
    reactions(29)%product_stoich(4) = 1.000000000000e+00_dp
    reactions(29)%forward_rate%pre_exponential = 5.000000000000e+12_dp
    reactions(29)%forward_rate%temperature_exponent = 0.000000000000e+00_dp
    reactions(29)%forward_rate%activation_energy = 7.250872000000e+07_dp
    reactions(29)%reversible = .true.

    ok = .true.
  end subroutine load_h2o2_full_mechanism

  subroutine load_h2o2_full_thermo_data(species, ok)
    type(nasa7_species), allocatable, intent(out) :: species(:)
    logical, intent(out) :: ok
    integer :: species_index

    allocate(species(10))

    species(1)%name = "H2"
    species(1)%molecular_weight = 2.016000000000e+00_dp
    species(1)%temperature_min = 2.000000000000e+02_dp
    species(1)%temperature_mid = 1.000000000000e+03_dp
    species(1)%temperature_max = 3.500000000000e+03_dp
    species(1)%low_coefficients = [ &
      2.344331120000e+00_dp, 7.980520750000e-03_dp, -1.947815100000e-05_dp, &
      2.015720940000e-08_dp, -7.376117610000e-12_dp, -9.179351730000e+02_dp, &
      6.830102380000e-01_dp ]
    species(1)%high_coefficients = [ &
      3.337279200000e+00_dp, -4.940247310000e-05_dp, 4.994567780000e-07_dp, &
      -1.795663940000e-10_dp, 2.002553760000e-14_dp, -9.501589220000e+02_dp, &
      -3.205023310000e+00_dp ]

    species(2)%name = "H"
    species(2)%molecular_weight = 1.008000000000e+00_dp
    species(2)%temperature_min = 2.000000000000e+02_dp
    species(2)%temperature_mid = 1.000000000000e+03_dp
    species(2)%temperature_max = 3.500000000000e+03_dp
    species(2)%low_coefficients = [ &
      2.500000000000e+00_dp, 7.053328190000e-13_dp, -1.995919640000e-15_dp, &
      2.300816320000e-18_dp, -9.277323320000e-22_dp, 2.547365990000e+04_dp, &
      -4.466828530000e-01_dp ]
    species(2)%high_coefficients = [ &
      2.500000010000e+00_dp, -2.308429730000e-11_dp, 1.615619480000e-14_dp, &
      -4.735152350000e-18_dp, 4.981973570000e-22_dp, 2.547365990000e+04_dp, &
      -4.466829140000e-01_dp ]

    species(3)%name = "O"
    species(3)%molecular_weight = 1.599900000000e+01_dp
    species(3)%temperature_min = 2.000000000000e+02_dp
    species(3)%temperature_mid = 1.000000000000e+03_dp
    species(3)%temperature_max = 3.500000000000e+03_dp
    species(3)%low_coefficients = [ &
      3.168267100000e+00_dp, -3.279318840000e-03_dp, 6.643063960000e-06_dp, &
      -6.128066240000e-09_dp, 2.112659710000e-12_dp, 2.912225920000e+04_dp, &
      2.051933460000e+00_dp ]
    species(3)%high_coefficients = [ &
      2.569420780000e+00_dp, -8.597411370000e-05_dp, 4.194845890000e-08_dp, &
      -1.001777990000e-11_dp, 1.228336910000e-15_dp, 2.921757910000e+04_dp, &
      4.784338640000e+00_dp ]

    species(4)%name = "O2"
    species(4)%molecular_weight = 3.199800000000e+01_dp
    species(4)%temperature_min = 2.000000000000e+02_dp
    species(4)%temperature_mid = 1.000000000000e+03_dp
    species(4)%temperature_max = 3.500000000000e+03_dp
    species(4)%low_coefficients = [ &
      3.782456360000e+00_dp, -2.996734160000e-03_dp, 9.847302010000e-06_dp, &
      -9.681295090000e-09_dp, 3.243728370000e-12_dp, -1.063943560000e+03_dp, &
      3.657675730000e+00_dp ]
    species(4)%high_coefficients = [ &
      3.282537840000e+00_dp, 1.483087540000e-03_dp, -7.579666690000e-07_dp, &
      2.094705550000e-10_dp, -2.167177940000e-14_dp, -1.088457720000e+03_dp, &
      5.453231290000e+00_dp ]

    species(5)%name = "OH"
    species(5)%molecular_weight = 1.700700000000e+01_dp
    species(5)%temperature_min = 2.000000000000e+02_dp
    species(5)%temperature_mid = 1.000000000000e+03_dp
    species(5)%temperature_max = 3.500000000000e+03_dp
    species(5)%low_coefficients = [ &
      3.992015430000e+00_dp, -2.401317520000e-03_dp, 4.617938410000e-06_dp, &
      -3.881133330000e-09_dp, 1.364114700000e-12_dp, 3.615080560000e+03_dp, &
      -1.039254580000e-01_dp ]
    species(5)%high_coefficients = [ &
      3.092887670000e+00_dp, 5.484297160000e-04_dp, 1.265052280000e-07_dp, &
      -8.794615560000e-11_dp, 1.174123760000e-14_dp, 3.858657000000e+03_dp, &
      4.476696100000e+00_dp ]

    species(6)%name = "H2O"
    species(6)%molecular_weight = 1.801500000000e+01_dp
    species(6)%temperature_min = 2.000000000000e+02_dp
    species(6)%temperature_mid = 1.000000000000e+03_dp
    species(6)%temperature_max = 3.500000000000e+03_dp
    species(6)%low_coefficients = [ &
      4.198640560000e+00_dp, -2.036434100000e-03_dp, 6.520402110000e-06_dp, &
      -5.487970620000e-09_dp, 1.771978170000e-12_dp, -3.029372670000e+04_dp, &
      -8.490322080000e-01_dp ]
    species(6)%high_coefficients = [ &
      3.033992490000e+00_dp, 2.176918040000e-03_dp, -1.640725180000e-07_dp, &
      -9.704198700000e-11_dp, 1.682009920000e-14_dp, -3.000429710000e+04_dp, &
      4.966770100000e+00_dp ]

    species(7)%name = "HO2"
    species(7)%molecular_weight = 3.300600000000e+01_dp
    species(7)%temperature_min = 2.000000000000e+02_dp
    species(7)%temperature_mid = 1.000000000000e+03_dp
    species(7)%temperature_max = 3.500000000000e+03_dp
    species(7)%low_coefficients = [ &
      4.301798010000e+00_dp, -4.749120510000e-03_dp, 2.115828910000e-05_dp, &
      -2.427638940000e-08_dp, 9.292251240000e-12_dp, 2.948080400000e+02_dp, &
      3.716662450000e+00_dp ]
    species(7)%high_coefficients = [ &
      4.017210900000e+00_dp, 2.239820130000e-03_dp, -6.336581500000e-07_dp, &
      1.142463700000e-10_dp, -1.079085350000e-14_dp, 1.118567130000e+02_dp, &
      3.785102150000e+00_dp ]

    species(8)%name = "H2O2"
    species(8)%molecular_weight = 3.401400000000e+01_dp
    species(8)%temperature_min = 2.000000000000e+02_dp
    species(8)%temperature_mid = 1.000000000000e+03_dp
    species(8)%temperature_max = 3.500000000000e+03_dp
    species(8)%low_coefficients = [ &
      4.276112690000e+00_dp, -5.428224170000e-04_dp, 1.673357010000e-05_dp, &
      -2.157708130000e-08_dp, 8.624543630000e-12_dp, -1.770258210000e+04_dp, &
      3.435050740000e+00_dp ]
    species(8)%high_coefficients = [ &
      4.165002850000e+00_dp, 4.908316940000e-03_dp, -1.901392250000e-06_dp, &
      3.711859860000e-10_dp, -2.879083050000e-14_dp, -1.786178770000e+04_dp, &
      2.916156620000e+00_dp ]

    species(9)%name = "AR"
    species(9)%molecular_weight = 3.995000000000e+01_dp
    species(9)%temperature_min = 3.000000000000e+02_dp
    species(9)%temperature_mid = 1.000000000000e+03_dp
    species(9)%temperature_max = 5.000000000000e+03_dp
    species(9)%low_coefficients = [ &
      2.500000000000e+00_dp, 0.000000000000e+00_dp, 0.000000000000e+00_dp, &
      0.000000000000e+00_dp, 0.000000000000e+00_dp, -7.453750000000e+02_dp, &
      4.366000000000e+00_dp ]
    species(9)%high_coefficients = [ &
      2.500000000000e+00_dp, 0.000000000000e+00_dp, 0.000000000000e+00_dp, &
      0.000000000000e+00_dp, 0.000000000000e+00_dp, -7.453750000000e+02_dp, &
      4.366000000000e+00_dp ]

    species(10)%name = "N2"
    species(10)%molecular_weight = 2.801400000000e+01_dp
    species(10)%temperature_min = 3.000000000000e+02_dp
    species(10)%temperature_mid = 1.000000000000e+03_dp
    species(10)%temperature_max = 5.000000000000e+03_dp
    species(10)%low_coefficients = [ &
      3.298677000000e+00_dp, 1.408240400000e-03_dp, -3.963222000000e-06_dp, &
      5.641515000000e-09_dp, -2.444854000000e-12_dp, -1.020899900000e+03_dp, &
      3.950372000000e+00_dp ]
    species(10)%high_coefficients = [ &
      2.926640000000e+00_dp, 1.487976800000e-03_dp, -5.684760000000e-07_dp, &
      1.009703800000e-10_dp, -6.753351000000e-15_dp, -9.227977000000e+02_dp, &
      5.980528000000e+00_dp ]

    ok = .true.
    do species_index = 1, size(species)
      if (.not. valid_nasa7_species(species(species_index))) then
        ok = .false.
        return
      end if
    end do
  end subroutine load_h2o2_full_thermo_data

  subroutine load_h2o2_full_transport_data( &
      names, geometries, values, ok)
    character(len=*), intent(out) :: names(:)
    integer, intent(out) :: geometries(:)
    real(dp), intent(out) :: values(:, :)
    logical, intent(out) :: ok

    ok = size(names) == 10 .and. &
      size(geometries) == 10 .and. &
      all(shape(values) == [5, 10])
    if (.not. ok) return
    names = ""
    geometries = 0
    values = 0.0_dp
    names(1) = "H2"
    geometries(1) = 1
    values(:, 1) = [ &
      3.800000000000e+01_dp, 2.920000000000e+00_dp, 0.000000000000e+00_dp, &
      7.900000000000e-01_dp, 2.800000000000e+02_dp ]
    names(2) = "H"
    geometries(2) = 0
    values(:, 2) = [ &
      1.450000000000e+02_dp, 2.050000000000e+00_dp, 0.000000000000e+00_dp, &
      0.000000000000e+00_dp, 0.000000000000e+00_dp ]
    names(3) = "O"
    geometries(3) = 0
    values(:, 3) = [ &
      8.000000000000e+01_dp, 2.750000000000e+00_dp, 0.000000000000e+00_dp, &
      0.000000000000e+00_dp, 0.000000000000e+00_dp ]
    names(4) = "O2"
    geometries(4) = 1
    values(:, 4) = [ &
      1.074000000000e+02_dp, 3.458000000000e+00_dp, 0.000000000000e+00_dp, &
      1.600000000000e+00_dp, 3.800000000000e+00_dp ]
    names(5) = "OH"
    geometries(5) = 1
    values(:, 5) = [ &
      8.000000000000e+01_dp, 2.750000000000e+00_dp, 0.000000000000e+00_dp, &
      0.000000000000e+00_dp, 0.000000000000e+00_dp ]
    names(6) = "H2O"
    geometries(6) = 2
    values(:, 6) = [ &
      5.724000000000e+02_dp, 2.605000000000e+00_dp, 1.844000000000e+00_dp, &
      0.000000000000e+00_dp, 4.000000000000e+00_dp ]
    names(7) = "HO2"
    geometries(7) = 2
    values(:, 7) = [ &
      1.074000000000e+02_dp, 3.458000000000e+00_dp, 0.000000000000e+00_dp, &
      0.000000000000e+00_dp, 1.000000000000e+00_dp ]
    names(8) = "H2O2"
    geometries(8) = 2
    values(:, 8) = [ &
      1.074000000000e+02_dp, 3.458000000000e+00_dp, 0.000000000000e+00_dp, &
      0.000000000000e+00_dp, 3.800000000000e+00_dp ]
    names(9) = "AR"
    geometries(9) = 0
    values(:, 9) = [ &
      1.365000000000e+02_dp, 3.330000000000e+00_dp, 0.000000000000e+00_dp, &
      0.000000000000e+00_dp, 0.000000000000e+00_dp ]
    names(10) = "N2"
    geometries(10) = 1
    values(:, 10) = [ &
      9.753000000000e+01_dp, 3.621000000000e+00_dp, 0.000000000000e+00_dp, &
      1.760000000000e+00_dp, 4.000000000000e+00_dp ]
  end subroutine load_h2o2_full_transport_data

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
