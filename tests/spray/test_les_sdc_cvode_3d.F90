program test_les_sdc_cvode_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use mixture_thermo_mod, only: mixture_density
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use les_sdc_3d_mod, only: advance_les_sdc_3d
  use sundials_forced_reactor_mod, only: advance_forced_reactor_cvode
  implicit none
  integer, parameter :: n=4,nv=12
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  real(dp) :: initial(nv,n,n,n),it(n,n,n),u(nv,n,n,n),t(n,n,n),native(nv,n,n,n)
  real(dp) :: ir(nv-4,n,n,n),q(nv),sound,h(3),boundary(nv),source(nv),diffusion
  real(dp) :: fx(nv,n+1,n,n),fy(nv,n,n+1,n),fz(nv,n,n,n+1),estimate(nv,n,n,n,3)
  real(dp) :: before(nv),after(nv),residual,delta,maximum
  real(dp), parameter :: dt=1e-10_dp
  character(len=16) :: models(4)
  logical :: ok,fail
  integer :: i,j,k,c,mode,calls
  models=[character(len=16)::'smagorinsky','wale','vreman','dynamic']
  call load_h2o2_elementary_thermo(species,ok);call require(ok,'thermo')
  call load_h2o2_elementary_mechanism(reactions,ok);call require(ok,'reactions')
  call load_h2o2_elementary_transport(transport,ok);call require(ok,'transport')
  q=0
  do c=1,size(species)
    select case(trim(species(c)%name))
    case('H2');q(5+c)=.1_dp
    case('O2');q(5+c)=.2_dp
    case('N2');q(5+c)=.7_dp
    end select
  end do
  call require(abs(sum(q(6:))-1)<1e-15_dp,'composition')
  do k=1,n
  do j=1,n
  do i=1,n
    q(2:4)=[real(2*j,dp),real(3*i,dp),real(k,dp)];q(5)=135000+10*j
    q(1)=mixture_density(species,q(6:),q(5),1100._dp,ok);call require(ok,'density')
    call reactive_primitive_to_conserved(species,q,initial(:,i,j,k),it(i,j,k),sound,ok)
    call require(ok,'state')
  end do
  end do
  end do
  h=.001_dp;maximum=0;delta=0
  do c=1,nv
    before(c)=sum(initial(c,:,:,:))*product(h)
  end do
  do mode=1,4
    u=initial;t=it;ir=0
    call advance_les_sdc_3d(species,u,t,h,dt,0._dp,2,.17_dp,.7_dp,.7_dp,boundary,source,diffusion,ok, &
      model=trim(models(mode)),transport=transport,reactions=reactions,reaction_source=ir, &
      external_force=[10._dp,-20._dp,30._dp], &
      chemistry_relative_tolerance=1e-10_dp,chemistry_absolute_tolerance=1e-13_dp)
    call require(ok,'native reacting LES reference');native=u
    u=initial;t=it;ir=0;estimate=0;calls=0;fail=.false.
    call run_step()
    call require(ok.and.calls==2*n**3,'CVODE visits every cell twice')
    call require(diffusion>0.and.maxval(abs(ir))>0,'active LES and reaction history')
    delta=max(delta,maxval(abs(u-native)/max(1._dp,abs(native))))
    call require(maxval(abs(u-native)/max(1._dp,abs(native)))<2e-7_dp,'early-time native comparison')
    do c=1,nv
      after(c)=sum(u(c,:,:,:))*product(h)
    end do
    residual=maxval(abs(after(:5)-before(:5)-boundary(:5)-source(:5))/max(abs(before(:5)),before(1)))
    maximum=max(maximum,residual)
    call require(residual<3e-12_dp,'Euler inventory including forcing')
    u=initial;t=it;ir=0;estimate=0;calls=0;fail=.true.
    boundary=42;source=43;diffusion=44;fx=45;fy=46;fz=47
    call run_step()
    call require(.not.ok.and.calls==n**3+1,'second-iteration reactor rejection')
    call require(all(u==initial).and.all(t==it).and.all(ir==0).and.all(estimate==0),'state and history rollback')
    call require(all(boundary==42).and.all(source==43).and.diffusion==44,'inventory rollback')
    call require(all(fx==45).and.all(fy==46).and.all(fz==47),'face flux rollback')
  end do
  print *, 'PASS: four LES closures with CVODE and rollback; native difference, conservation:',delta,maximum
contains
  subroutine run_step()
    call advance_les_sdc_3d(species,u,t,h,dt,0._dp,2,.17_dp,.7_dp,.7_dp,boundary,source,diffusion,ok, &
      model=trim(models(mode)),transport=transport,reactions=reactions,reaction_source=ir, &
      external_force=[10._dp,-20._dp,30._dp], &
      chemistry_relative_tolerance=1e-10_dp,chemistry_absolute_tolerance=1e-13_dp, &
      flux_x=fx,flux_y=fy,flux_z=fz,new_source_estimate=estimate,evaluate_reaction=reactor)
  end subroutine
  subroutine reactor(sp,rx,e0,external_y,external_e,duration,rtol,atol,densities,temp,count,good,maximum_steps)
    type(nasa7_species), intent(in) :: sp(:)
    type(elementary_reaction), intent(in) :: rx(:)
    real(dp), intent(in) :: e0,external_y(:),external_e,duration,rtol,atol
    real(dp), intent(inout) :: densities(:),temp
    integer, intent(out) :: count
    logical, intent(out) :: good
    integer, optional, intent(in) :: maximum_steps
    calls=calls+1
    if(fail.and.calls==n**3+1) then
      densities=-42;temp=-7;count=0;good=.false.;return
    end if
    call advance_forced_reactor_cvode(sp,rx,e0,external_y,external_e,duration,rtol,atol, &
      densities,temp,count,good,maximum_steps)
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
