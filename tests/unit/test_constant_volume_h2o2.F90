program test_constant_volume_h2o2
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: &
    h2o2_nspecies, h2o2_h2_index, h2o2_h_index, h2o2_o_index, &
    h2o2_o2_index, h2o2_oh_index, h2o2_h2o_index, &
    load_h2o2_elementary_mechanism
  use constant_volume_reactor_mod, only: &
    reactor_specific_internal_energy, advance_constant_volume_adaptive
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp) :: mole_fractions(h2o2_nspecies), mass_fractions(h2o2_nspecies)
  real(dp) :: initial_mass_fractions(h2o2_nspecies)
  real(dp) :: density, temperature, target_energy, current_energy
  real(dp) :: accepted_dt, next_dt, requested_dt, time
  real(dp) :: initial_h_atoms, initial_o_atoms, h_atoms, o_atoms
  logical :: ok
  integer :: steps

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "thermo load failed"
  call load_h2o2_elementary_mechanism(reactions, ok)
  if (.not. ok) error stop "mechanism load failed"
  mole_fractions = [ &
    0.295_dp, 1.0e-6_dp, 1.0e-6_dp, 0.1475_dp, 1.0e-6_dp, &
    0.0_dp, 0.557497_dp ]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  if (.not. ok) error stop "composition conversion failed"
  initial_mass_fractions = mass_fractions
  temperature = 1200.0_dp
  density = mixture_density( &
    species, mass_fractions, 101325.0_dp, temperature, ok)
  if (.not. ok) error stop "density evaluation failed"
  target_energy = reactor_specific_internal_energy( &
    species, mass_fractions, temperature, ok)
  if (.not. ok) error stop "initial energy evaluation failed"

  initial_h_atoms = atom_inventory_h(mass_fractions)
  initial_o_atoms = atom_inventory_o(mass_fractions)
  requested_dt = 1.0e-9_dp
  time = 0.0_dp
  steps = 0
  do while (time < 2.0e-6_dp)
    requested_dt = min(requested_dt, 2.0e-6_dp - time)
    call advance_constant_volume_adaptive( &
      species, reactions, density, target_energy, requested_dt, 1.0e-8_dp, &
      1.0e-14_dp, mass_fractions, temperature, accepted_dt, next_dt, ok)
    if (.not. ok) error stop "adaptive reactor step failed"
    time = time + accepted_dt
    requested_dt = min(2.0e-7_dp, next_dt)
    steps = steps + 1
    if (steps > 100000) error stop "adaptive reactor did not finish"
  end do

  current_energy = reactor_specific_internal_energy( &
    species, mass_fractions, temperature, ok)
  if (.not. ok) error stop "final energy evaluation failed"
  h_atoms = atom_inventory_h(mass_fractions)
  o_atoms = atom_inventory_o(mass_fractions)

  if (abs(sum(mass_fractions) - 1.0_dp) > 2.0e-12_dp) &
    error stop "reactor composition closure failed"
  if (any(mass_fractions < 0.0_dp)) &
    error stop "reactor generated a negative mass fraction"
  if (abs(current_energy - target_energy) > &
      2.0e-9_dp * max(1.0_dp, abs(target_energy))) &
    error stop "reactor energy conservation failed"
  if (abs(h_atoms - initial_h_atoms) > 2.0e-10_dp) &
    error stop "reactor H inventory changed"
  if (abs(o_atoms - initial_o_atoms) > 2.0e-10_dp) &
    error stop "reactor O inventory changed"
  if (maxval(abs(mass_fractions - initial_mass_fractions)) < 1.0e-10_dp) &
    error stop "reactor did not evolve"

  write(*, '(a,i0)') "steps=", steps
  write(*, '(a,es24.16)') "temperature=", temperature
  write(*, '(a,es24.16)') "maximum_composition_change=", &
    maxval(abs(mass_fractions - initial_mass_fractions))
  write(*, '(a)') "test_constant_volume_h2o2: PASS"

contains

  real(dp) function atom_inventory_h(y) result(value)
    real(dp), intent(in) :: y(:)
    value = 2.0_dp * y(h2o2_h2_index) / species(h2o2_h2_index)%molecular_weight + &
      y(h2o2_h_index) / species(h2o2_h_index)%molecular_weight + &
      y(h2o2_oh_index) / species(h2o2_oh_index)%molecular_weight + &
      2.0_dp * y(h2o2_h2o_index) / species(h2o2_h2o_index)%molecular_weight
  end function atom_inventory_h

  real(dp) function atom_inventory_o(y) result(value)
    real(dp), intent(in) :: y(:)
    value = y(h2o2_o_index) / species(h2o2_o_index)%molecular_weight + &
      2.0_dp * y(h2o2_o2_index) / species(h2o2_o2_index)%molecular_weight + &
      y(h2o2_oh_index) / species(h2o2_oh_index)%molecular_weight + &
      y(h2o2_h2o_index) / species(h2o2_h2o_index)%molecular_weight
  end function atom_inventory_o

end program test_constant_volume_h2o2
