program test_rotating_reaction_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species,nasa7_mass_properties
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_primitive_to_conserved
  use reactive_3d_mod, only: react_sdc_state_3d
  use reactive_eb_3d_driver_mod, only: react_sdc_eb_state_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp), allocatable :: q(:),old(:,:,:,:),rot(:,:,:,:),rates(:,:,:,:),rr(:,:,:,:),u(:,:,:,:),v(:,:,:,:), &
    ir(:,:,:,:),jr(:,:,:,:)
  real(dp) :: heat(1,1,1),rot_heat(1,1,1),expected,cp,cv,h,e,entropy
  real(dp) :: t(1,1,1),tr(1,1,1),ts(1,1,1),sound,phi(1,1,1),error,dt
  integer :: nv,k,h2,n2
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok);call require(ok,'thermo')
  call load_h2o2_elementary_mechanism(reactions,ok);call require(ok,'reactions')
  nv=reactive_nvar(size(species));dt=1e-7_dp
  allocate(q(reactive_nprim(size(species))),old(nv,1,1,1),rot(nv,1,1,1),rates(nv,1,1,1),rr(nv,1,1,1), &
    u(nv,1,1,1),v(nv,1,1,1),ir(nv-4,1,1,1),jr(nv-4,1,1,1))
  q=0;q(1)=1;q(2:4)=[10._dp,-20._dp,30._dp];q(5)=6e5_dp
  do k=1,size(species)
    select case(species(k)%name)
    case('H2');q(5+k)=.028_dp;h2=5+k
    case('O2');q(5+k)=.222_dp
    case('OH');q(5+k)=1e-4_dp
    case('N2');q(5+k)=.7499_dp;n2=5+k
    end select
  end do
  call reactive_primitive_to_conserved(species,q,old(:,1,1,1),t(1,1,1),sound,ok)
  call require(ok,'initial state')
  phi=2e6_dp;rot=old;rot(5,1,1,1)=old(5,1,1,1)-old(1,1,1,1)*phi(1,1,1)
  call require(rot(5,1,1,1)<0,'negative rotating energy')
  rr=0;rr(1,1,1,1)=3;rr(n2,1,1,1)=3;rr(2:4,1,1,1)=[2000._dp,-1000._dp,300._dp]
  rr(5,1,1,1)=3e5_dp
  rates=rr;rates(5,1,1,1)=rr(5,1,1,1)+phi(1,1,1)*rr(1,1,1,1)
  call react_sdc_state_3d(species,reactions,old,t,rates,dt,1e-8_dp,1e-12_dp,u,ts,ir,ok,heat_release=heat)
  call require(ok,'ordinary forced reaction')
  call react_sdc_state_3d(species,reactions,rot,t,rr,dt,1e-8_dp,1e-12_dp,v,tr,jr,ok,rotating_potential=phi,heat_release=rot_heat)
  call require(ok,'rotating forced reaction')
  expected=0
  do k=1,size(species)
    call nasa7_mass_properties(species(k),ts(1,1,1),cp,cv,h,e,entropy,ok)
    call require(ok,'accepted-time heat-release enthalpy')
    expected=expected-h*ir(k,1,1,1)
  enddo
  call require(abs(heat(1,1,1)-expected)<=16*epsilon(1._dp)*max(1._dp,abs(expected)), &
    'stage heat release excludes external species forcing')
  call require(abs(rot_heat(1,1,1)-heat(1,1,1))<=3e-10_dp*max(1._dp,abs(heat(1,1,1))), &
    'heat release uses thermal enthalpy independent of rotating energy')

  v(5,1,1,1)=v(5,1,1,1)+v(1,1,1,1)*phi(1,1,1)
  error=maxval(abs(u-v)/max(1._dp,abs(u)))
  call require(error<3e-11_dp.and.abs(tr(1,1,1)/ts(1,1,1)-1)<3e-11_dp,'energy convention invariance')
  call require(abs(u(h2,1,1,1)-old(h2,1,1,1))>1e-12_dp,'actual chemical conversion')
  phi=0
  call react_sdc_state_3d(species,reactions,old,t,rates,dt,1e-8_dp,1e-12_dp,v,tr,jr,ok,rotating_potential=phi)
  call require(ok.and.all(v==u).and.all(tr==ts).and.all(jr==ir),'exact zero potential identity')
  call check_old_internal_energy()
  phi=ieee_value(0._dp,ieee_quiet_nan);v=17;tr=19;jr=23;rot_heat=29
  call react_sdc_state_3d(species,reactions,old,t,rates,dt,1e-8_dp,1e-12_dp,v,tr,jr,ok,rotating_potential=phi,heat_release=rot_heat)
  call require(.not.ok.and.all(v==17).and.all(tr==19).and.all(jr==23).and.all(rot_heat==29),'invalid potential rollback')
  block
    type(eb_geometry_3d) :: geometry
    real(dp) :: eb_old(nv,2,1,1),eb_rates(nv,2,1,1),eb_state(nv,2,1,1),eb_ir(nv-4,2,1,1)
    real(dp) :: eb_t(2,1,1),eb_heat(2,1,1),eb_initial_t(2,1,1),eb_held(2,1,1),energy_correction
    call build_axis_plane_eb_geometry_3d(2,1,1,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
      'x',.6_dp,geometry,ok)
    call require(ok,'heat-release EB geometry')
    call require(geometry%volume_fraction(1,1,1)==0.and.geometry%volume_fraction(2,1,1)>0, &
      'heat-release fixture contains solid and physical cells')
    eb_old(:,1,1,1)=ieee_value(0._dp,ieee_quiet_nan);eb_old(:,2,1,1)=old(:,1,1,1)
    eb_rates(:,1,1,1)=ieee_value(0._dp,ieee_quiet_nan);eb_rates(:,2,1,1)=rates(:,1,1,1)
    eb_initial_t(1,1,1)=ieee_value(0._dp,ieee_quiet_nan);eb_initial_t(2,1,1)=t(1,1,1)
    eb_state=17;eb_t=19;eb_ir=23;eb_heat=29
    call react_sdc_eb_state_3d(species,reactions,eb_old,eb_initial_t,geometry,eb_rates,dt,1e-8_dp,1e-12_dp, &
      eb_state,eb_t,eb_ir,ok,heat_release=eb_heat)
    call require(ok,'EB heat release ignores poisoned covered-cell thermo')
    call require(eb_heat(1,1,1)==0.and.eb_heat(2,1,1)==heat(1,1,1),'EB physical heat matches cell stage')
    call require(all(eb_state(:,2,1,1)==u(:,1,1,1)).and.all(eb_ir(:,2,1,1)==ir(:,1,1,1)), &
      'optional EB heat leaves physical gas and chemical source unchanged')
    eb_state=17;eb_t=19;eb_ir=23;eb_held=29
    call react_sdc_eb_state_3d(species,reactions,eb_old,eb_initial_t,geometry,eb_rates,dt,1e-8_dp,1e-12_dp, &
      eb_state,eb_t,eb_ir,ok,internal_energy_density=eb_held)
    call require(ok,'EB forced-reactor retained endpoint')
    call require(eb_held(1,1,1)==29,'covered retained scalar preserved')
    call require(abs(eb_held(2,1,1)-(eb_state(5,2,1,1)- &
      .5_dp*sum(eb_state(2:4,2,1,1)**2)/eb_state(1,2,1,1)))< &
      2e-12_dp*max(1._dp,abs(eb_held(2,1,1))),'retained forced endpoint matches accepted gas')
    energy_correction=7
    call react_sdc_eb_state_3d(species,reactions,eb_old,eb_initial_t,geometry,eb_rates,dt,1e-8_dp,1e-12_dp, &
      eb_state,eb_t,eb_ir,ok,internal_energy_density=eb_held,energy_correction=energy_correction)
    call require(ok.and.energy_correction==7.and.eb_held(1,1,1)==29, &
      'EB retained EOS policy preserves covered scalar and correction ledger')
    eb_state=17;eb_t=19;eb_ir=23;eb_held=29
    call react_sdc_eb_state_3d(species,reactions,eb_old,eb_initial_t,geometry,eb_rates,-dt,1e-8_dp,1e-12_dp, &
      eb_state,eb_t,eb_ir,ok,internal_energy_density=eb_held)
    call require(.not.ok.and.all(eb_held==29),'EB retained endpoint rollback')
    eb_state=17;eb_t=19;eb_ir=23;eb_heat=29
    call react_sdc_eb_state_3d(species,reactions,eb_old,eb_initial_t,geometry,eb_rates,-dt,1e-8_dp,1e-12_dp, &
      eb_state,eb_t,eb_ir,ok,heat_release=eb_heat)
    call require(.not.ok.and.all(eb_state==17).and.all(eb_t==19).and.all(eb_ir==23).and.all(eb_heat==29), &
      'EB heat-release rejection preserves all output arrays')
  end block
  print *, 'PASS: rotating forced chemistry; scaled state difference, H2 change:',error,u(h2,1,1,1)-old(h2,1,1,1)
contains
  subroutine check_old_internal_energy()
    type(eb_geometry_3d) :: geometry
    real(dp) :: shadow(nv,1,1,1),ref(nv,1,1,1),actual(nv,1,1,1),ri(nv-4,1,1,1),ai(nv-4,1,1,1)
    real(dp) :: rt(1,1,1),at(1,1,1),initial(1,1,1),endpoint(1,1,1),trial(nv,1,1,1),projected,offset
    real(dp) :: eb_old(nv,2,1,1),eb_rates(nv,2,1,1),eb_state(nv,2,1,1),eb_ir(nv-4,2,1,1)
    real(dp) :: eb_t(2,1,1),eb_initial(2,1,1),eb_endpoint(2,1,1),eb_old_t(2,1,1)
    ! A chemically equivalent shifted-total state provides an independent
    ! native trajectory for the held initial energy. Publication still uses
    ! the original total/momentum-derived endpoint, as pinned React.cpp does.
    offset=16384._dp;shadow=old;shadow(5,1,1,1)=shadow(5,1,1,1)+offset
    initial=old(5,1,1,1)-.5_dp*sum(old(2:4,1,1,1)**2)/old(1,1,1,1)+offset
    call react_sdc_state_3d(species,reactions,shadow,t,rates,dt,1e-8_dp,1e-12_dp,ref,rt,ri,ok)
    call require(ok,'shifted-energy native reference')
    endpoint=29
    call react_sdc_state_3d(species,reactions,old,t,rates,dt,1e-8_dp,1e-12_dp,actual,at,ai,ok, &
      internal_energy_density=endpoint,old_internal_energy_density=initial)
    call require(ok,'held old energy feeds native reaction')
    call require(all(abs(actual(6:,:,:,:)-ref(6:,:,:,:))<=64*epsilon(1._dp)*max(1._dp,abs(ref(6:,:,:,:)))), &
      'held initial energy matches chemically equivalent trajectory')
    call require(any(actual(6:,:,:,:)/=u(6:,:,:,:)),'held old energy changes actual chemical conversion')
    trial=old+dt*rates
    projected=trial(5,1,1,1)-.5_dp*sum(trial(2:4,1,1,1)**2)/trial(1,1,1,1)
    call require(endpoint(1,1,1)==projected,'held output retains projected forced endpoint')
    call require(initial(1,1,1)==old(5,1,1,1)-.5_dp*sum(old(2:4,1,1,1)**2)/old(1,1,1,1)+offset, &
      'held old input remains unchanged')
    call require(all(abs(ai(:nv-5,:,:,:)-ri(:nv-5,:,:,:))<= &
      64*epsilon(1._dp)*max(1._dp,abs(ri(:nv-5,:,:,:)))),'held reaction-source trajectory')
    actual=17;at=19;ai=23;endpoint=29
    call react_sdc_state_3d(species,reactions,old,t,rates,dt,1e-8_dp,1e-12_dp,actual,at,ai,ok, &
      evaluate_reaction=refuse_reactor,internal_energy_density=endpoint,old_internal_energy_density=initial)
    call require(.not.ok.and.all(actual==17).and.all(at==19).and.all(ai==23).and.all(endpoint==29), &
      'held supplied-reactor rejection is atomic')
    call build_axis_plane_eb_geometry_3d(2,1,1,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',.6_dp,geometry,ok)
    call require(ok,'held initial EB geometry')
    eb_old(:,1,1,1)=ieee_value(0._dp,ieee_quiet_nan);eb_old(:,2,1,1)=old(:,1,1,1)
    eb_rates(:,1,1,1)=ieee_value(0._dp,ieee_quiet_nan);eb_rates(:,2,1,1)=rates(:,1,1,1)
    eb_old_t(1,1,1)=ieee_value(0._dp,ieee_quiet_nan);eb_old_t(2,1,1)=t(1,1,1)
    eb_initial(1,1,1)=ieee_value(0._dp,ieee_quiet_nan);eb_initial(2,1,1)=initial(1,1,1)
    eb_state=17;eb_t=19;eb_ir=23;eb_endpoint=29
    call react_sdc_eb_state_3d(species,reactions,eb_old,eb_old_t,geometry,eb_rates,dt,1e-8_dp,1e-12_dp, &
      eb_state,eb_t,eb_ir,ok,internal_energy_density=eb_endpoint,old_internal_energy_density=eb_initial)
    call require(ok.and.eb_endpoint(1,1,1)==29,'held EB input ignores poisoned covered scalar')
    call require(all(eb_state(6:,2,1,1)==ref(6:,1,1,1)).and.eb_endpoint(2,1,1)==projected, &
      'held EB input preserves native trajectory and projected endpoint')
  end subroutine check_old_internal_energy

  subroutine refuse_reactor(species,reactions,initial_energy_density,external_species,external_energy, &
      duration,rtol,atol,species_density,temperature,completed_steps,ok,maximum_steps)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: initial_energy_density,external_species(:),external_energy,duration,rtol,atol
    real(dp), intent(inout) :: species_density(:),temperature
    integer, intent(out) :: completed_steps
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_steps
    species_density=-1;temperature=-1;completed_steps=0;ok=.false.
  end subroutine refuse_reactor

  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
