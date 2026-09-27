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
  call require(maxval(abs(rho_y-initial-duration*source))>1e-9_dp,'reaction must be active')
  print *, 'Forced native/CVODE relative difference:',relative
  rho_y=initial;temperature=1000
  call advance_forced_reactor_cvode(species,reactions,energy,source,energy_source,1e-3_dp, &
    1e-12_dp,1e-15_dp,rho_y,temperature,steps,ok,maximum_steps=1)
  call require(.not.ok.and.all(rho_y==initial).and.temperature==1000.and.steps==0,'step-limit rollback')
  call advance_forced_reactor_cvode(species,inert,energy,source,energy_source,1e-4_dp, &
    1e-10_dp,1e-13_dp,rho_y,temperature,steps,ok)
  call require(ok,'context released after failed solve')
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
