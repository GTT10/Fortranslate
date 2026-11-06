program test_sundials_forced_reactor
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use mixture_thermo_mod, only: mixture_mass_properties
  use constant_volume_reactor_mod, only: advance_forced_reactor_interval
  use sundials_forced_reactor_mod, only: advance_forced_reactor_cvode
  use fixture_mechanism_mod, only: load_fixture_thermo_data,load_fixture_mechanism
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(elementary_reaction) :: inert(0)
  type(nasa7_species) :: analytic_species(3)
  type(elementary_reaction) :: analytic_reactions(1)
  real(dp) :: analytic_initial(3),analytic_state(3),analytic_expected(3),analytic_source(3)
  real(dp) :: initial(2),rho_y(2),native(2),source(2),temperature,nt,energy,mw,r,cp,cv,gamma,h,s,e
  real(dp) :: energy_source,duration,expected(2),relative
  integer :: steps,native_steps,i
  logical :: ok
  call load_fixture_thermo_data(species,ok)
  call require(ok,'thermodynamics')
  call load_fixture_mechanism(reactions,ok)
  call require(ok,'mechanism')
  initial=[.8_dp,.2_dp];temperature=1000
  call mixture_mass_properties(species,initial,temperature,mw,r,cp,cv,gamma,h,energy,s,ok)
  call require(ok,'initial internal energy')
  source=[.1_dp,-.02_dp];energy_source=1e6_dp;duration=1e-4_dp
  do i=1,70
    rho_y=initial;temperature=1000
    call advance_forced_reactor_cvode(species,inert,energy,source,energy_source,duration, &
      1e-10_dp,1e-13_dp,rho_y,temperature,steps,ok)
    call require(ok.and.steps>0,'inert forced integration and context reuse')
    expected=initial+duration*source
    call require(maxval(abs(rho_y-expected))<1e-12_dp,'linear density-source oracle')
    call mixture_mass_properties(species,rho_y/sum(rho_y),temperature,mw,r,cp,cv,gamma,h,e,s,ok)
    call require(ok,'final forced thermodynamics')
    call require(abs(sum(rho_y)*e-energy-duration*energy_source)<1e-8_dp*max(1._dp,abs(energy)), &
      'linear energy-source oracle')
  enddo
  duration=1e-7_dp;rho_y=initial;native=initial;temperature=1000;nt=1000
  call advance_forced_reactor_cvode(species,reactions,energy,source,energy_source,duration, &
    1e-10_dp,1e-13_dp,rho_y,temperature,steps,ok)
  call require(ok,'active forced CVODE reaction')
  call require(abs(sum(rho_y)-sum(initial)-duration*sum(source))<1e-11_dp,'reacting total mass source')
  call mixture_mass_properties(species,rho_y/sum(rho_y),temperature,mw,r,cp,cv,gamma,h,e,s,ok)
  call require(ok.and.abs(sum(rho_y)*e-energy-duration*energy_source)<1e-8_dp*max(1._dp,abs(energy)), &
    'reacting total internal energy source')
  ! Loose integration accuracy must not relax the exact total-density identity.
  rho_y=initial;temperature=1000
  call advance_forced_reactor_cvode(species,reactions,energy,source,energy_source,duration, &
    2e-7_dp,1e-12_dp,rho_y,temperature,steps,ok)
  call require(ok,'active forced CVODE at application tolerances')
  call require(abs(sum(rho_y)-(sum(initial)+duration*sum(source)))<4*epsilon(1._dp), &
    'analytic total density independent of integration tolerance')
  ! The native low-order solver is used only over an early, resolved interval.
  duration=1e-10_dp;rho_y=initial;temperature=1000
  call advance_forced_reactor_cvode(species,reactions,energy,source,energy_source,duration, &
    1e-10_dp,1e-13_dp,rho_y,temperature,steps,ok)
  call require(ok,'early forced CVODE reaction')
  call advance_forced_reactor_interval(species,reactions,energy,source,energy_source,duration, &
    1e-10_dp,1e-13_dp,native,nt,native_steps,ok)
  call require(ok,'independent native forced reaction')
  relative=maxval(abs(rho_y-native))/maxval(abs(native))
  call require(relative<2e-7_dp.and.abs(temperature-nt)<2e-4_dp,'active native/CVODE agreement')
  ! A -> B with a carrier and frozen sources 1e-10*(1,-0.1,0):
  ! A=1e-10*(1-exp(-10t))/10, B=0.9e-10*t-A, C=1. B is initially negative but
  ! positive at t=0.1. Per-substep nonnegativity wrongly stalls this valid
  ! forced endpoint. The transient deficit stays inside the existing signed
  ! trial bound. Identical thermodynamics preserve exactly constant T.
  analytic_species=species(1)
  analytic_species(1)%name='A';analytic_species(2)%name='B';analytic_species(3)%name='C'
  analytic_reactions(1)%reactant_stoich=[1._dp,0._dp,0._dp]
  analytic_reactions(1)%product_stoich=[0._dp,1._dp,0._dp]
  analytic_reactions(1)%forward_rate%pre_exponential=10._dp
  analytic_reactions(1)%reversible=.false.
  analytic_initial=[0._dp,0._dp,1._dp];analytic_source=[1e-10_dp,-1e-11_dp,0._dp]
  call mixture_mass_properties(analytic_species,analytic_initial,1000._dp,mw,r,cp,cv,gamma,h,e,s,ok)
  call require(ok,'analytic forced thermodynamics')
  analytic_state=analytic_initial;temperature=1000
  call advance_forced_reactor_cvode(analytic_species,analytic_reactions,e,analytic_source,.9e-10_dp*e,.1_dp, &
    1e-10_dp,1e-19_dp,analytic_state,temperature,steps,ok,maximum_steps=5000)
  analytic_expected(1)=1e-10_dp*(1-exp(-1._dp))/10
  analytic_expected(2)=.09e-10_dp-analytic_expected(1);analytic_expected(3)=1
  call require(ok.and.steps>0.and.maxval(abs(analytic_state(:2)-analytic_expected(:2)))<1e-17_dp.and. &
    abs(analytic_state(3)-1._dp)<8*epsilon(1._dp), &
    'signed forced intermediate with analytic nonnegative endpoint')
  call require(abs(sum(analytic_state)-(1+.09e-10_dp))<8*epsilon(1._dp).and.abs(temperature-1000)<1e-7_dp, &
    'forced analytic mass and energy identities')
  analytic_state=analytic_initial;temperature=1000
  call advance_forced_reactor_cvode(analytic_species,analytic_reactions,e,analytic_source,.9e-10_dp*e,1e-3_dp, &
    1e-10_dp,1e-19_dp,analytic_state,temperature,steps,ok,maximum_steps=5000)
  call require(.not.ok.and.steps==0.and.all(analytic_state==analytic_initial).and.temperature==1000, &
    'negative forced endpoint still refuses with exact rollback')
  call require(maxval(abs(rho_y-initial-duration*source))>1e-9_dp,'reaction must be active')
  print *, 'Forced native/CVODE relative difference:',relative
  rho_y=initial;temperature=1000
  call advance_forced_reactor_cvode(species,reactions,energy,source,energy_source,1e-3_dp, &
    1e-12_dp,1e-15_dp,rho_y,temperature,steps,ok,maximum_steps=1)
  call require(.not.ok.and.all(rho_y==initial).and.temperature==1000.and.steps==0,'step-limit rollback')
  call advance_forced_reactor_cvode(species,inert,energy,source,energy_source,1e-4_dp, &
    1e-10_dp,1e-13_dp,rho_y,temperature,steps,ok)
  call require(ok,'context released after failed solve')
  rho_y=[1._dp,0._dp];temperature=1000
  call mixture_mass_properties(species,rho_y,temperature,mw,r,cp,cv,gamma,h,e,s,ok)
  call require(ok,'zero-species fixture energy')
  call advance_forced_reactor_cvode(species,reactions,e,[0._dp,0._dp],0._dp,1e-7_dp, &
    2e-7_dp,1e-12_dp,rho_y,temperature,steps,ok)
  call require(ok.and.all(rho_y>=0).and.abs(sum(rho_y)-1)<4*epsilon(1._dp), &
    'absent species constrained endpoint and exact mass')
  call mixture_mass_properties(species(:1),[1._dp],1000._dp,mw,r,cp,cv,gamma,h,e,s,ok)
  call require(ok,'single-species thermodynamics')
  rho_y(1)=.8_dp;temperature=1000
  call advance_forced_reactor_cvode(species(:1),inert,.8_dp*e,[.1_dp],1e6_dp,1e-4_dp, &
    1e-10_dp,1e-13_dp,rho_y(:1),temperature,steps,ok)
  call require(ok.and.abs(rho_y(1)-.80001_dp)<epsilon(1._dp),'single-species analytic density')
contains
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if(.not.condition) then
      print *, 'FAIL: ',message
      error stop 1
    endif
  end subroutine
end program
