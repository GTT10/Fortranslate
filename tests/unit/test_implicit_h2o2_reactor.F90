program test_implicit_h2o2_reactor
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density, mixture_mass_properties
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_full_thermo_mod, only: full_nspecies, load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use constant_volume_reactor_mod, only: &
    advance_constant_volume_implicit_adaptive,forced_reactor_rhs,reactor_rhs
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp) :: x(full_nspecies), y(full_nspecies), density, temperature
  real(dp) :: target_energy, final_energy, molecular_weight, gas_constant
  real(dp) :: cp, cv, gamma, enthalpy, entropy, accepted, next_step
  integer :: newton_iterations, rejected
  logical :: ok
  type(elementary_reaction), allocatable :: inert(:)
  real(dp) :: rhoy(full_nspecies),source(full_nspecies),rates(full_nspecies),baseline(full_nspecies),source_t

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
  allocate(inert(0));rhoy=density*y;source=0
  call reactor_rhs(species,reactions,density,target_energy,y,temperature,baseline,source_t,ok)
  if(.not.ok) error stop 'unforced RHS baseline'
  call forced_reactor_rhs(species,reactions,density*target_energy,0._dp,rhoy,source,0._dp, &
    temperature,rates,source_t,ok)
  if(.not.ok.or.maxval(abs(rates-density*baseline)/max(1._dp,abs(density*baseline)))>1e-12_dp) &
    error stop 'forced RHS zero-source chemistry limit'
  source=.2_dp*y;rhoy=(density+.05_dp)*y
  call forced_reactor_rhs(species,inert,density*target_energy,.25_dp,rhoy,source,.2_dp*target_energy, &
    temperature,rates,source_t,ok)
  if(.not.ok.or.any(rates/=source).or.abs(source_t-temperature)>1e-8_dp) &
    error stop 'forced RHS variable density and energy'
  call reactor_rhs(species,reactions,density+.05_dp,target_energy,y,temperature,baseline,source_t,ok)
  if(.not.ok) error stop 'forced chemistry comparison baseline'
  call forced_reactor_rhs(species,reactions,density*target_energy,.25_dp,rhoy,source,.2_dp*target_energy, &
    temperature,rates,source_t,ok)
  if(.not.ok.or.maxval(abs(rates-((density+.05_dp)*baseline+source))/ &
    max(1._dp,abs((density+.05_dp)*baseline+source)))>1e-12_dp) error stop 'chemistry and external species rates'
  rhoy(1)=-1
  call forced_reactor_rhs(species,reactions,density*target_energy,.25_dp,rhoy,source,.2_dp*target_energy, &
    temperature,rates,source_t,ok)
  if(ok.or.any(rates/=0).or.source_t/=0) error stop 'forced RHS invalid species zero outputs'
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
