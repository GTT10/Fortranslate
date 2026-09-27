program test_rotating_reaction_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_primitive_to_conserved
  use reactive_3d_mod, only: react_sdc_state_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp), allocatable :: q(:),old(:,:,:,:),rot(:,:,:,:),rates(:,:,:,:),rr(:,:,:,:),u(:,:,:,:),v(:,:,:,:), &
    ir(:,:,:,:),jr(:,:,:,:)
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
  call react_sdc_state_3d(species,reactions,old,t,rates,dt,1e-8_dp,1e-12_dp,u,ts,ir,ok)
  call require(ok,'ordinary forced reaction')
  call react_sdc_state_3d(species,reactions,rot,t,rr,dt,1e-8_dp,1e-12_dp,v,tr,jr,ok,rotating_potential=phi)
  call require(ok,'rotating forced reaction')
  v(5,1,1,1)=v(5,1,1,1)+v(1,1,1,1)*phi(1,1,1)
  error=maxval(abs(u-v)/max(1._dp,abs(u)))
  call require(error<3e-11_dp.and.abs(tr(1,1,1)/ts(1,1,1)-1)<3e-11_dp,'energy convention invariance')
  call require(abs(u(h2,1,1,1)-old(h2,1,1,1))>1e-12_dp,'actual chemical conversion')
  phi=0
  call react_sdc_state_3d(species,reactions,old,t,rates,dt,1e-8_dp,1e-12_dp,v,tr,jr,ok,rotating_potential=phi)
  call require(ok.and.all(v==u).and.all(tr==ts).and.all(jr==ir),'exact zero potential identity')
  phi=ieee_value(0._dp,ieee_quiet_nan);v=17;tr=19;jr=23
  call react_sdc_state_3d(species,reactions,old,t,rates,dt,1e-8_dp,1e-12_dp,v,tr,jr,ok,rotating_potential=phi)
  call require(.not.ok.and.all(v==17).and.all(tr==19).and.all(jr==23),'invalid potential rollback')
  print *, 'PASS: rotating forced chemistry; scaled state difference, H2 change:',error,u(h2,1,1,1)-old(h2,1,1,1)
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
