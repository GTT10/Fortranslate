program test_spray_eb_field
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d
  use spray_eb_field_mod, only: eb_spray_stencil,sample_eb_spray_field,deposit_eb_spray_gas
  implicit none
  type(eb_geometry_3d) :: g
  type(nasa7_species), allocatable :: species(:)
  real(dp) :: phi(0:4,0:4,0:4),positions(3,2),w(8),v(8),t,sound,r,temperature(4,4,4),saved_t(4,4,4)
  real(dp) :: field(1,4,4,4)
  real(dp), allocatable :: primitive(:),cell(:),y(:),state(:,:,:,:),saved(:,:,:,:),feedback(:,:),samples(:,:),before(:),after(:)
  integer :: i,j,k,nv,nitrogen,c(3,8),q
  logical :: ok
  do i=0,4
    phi(i,:,:)=real(i,dp)-1.25_dp
  end do
  call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[4._dp,4._dp,4._dp],g,ok)
  call require(ok,'geometry')
  positions(:,1)=[1.4_dp,1.5_dp,1.5_dp];positions(:,2)=[2._dp,2._dp,2._dp]
  c=0;w=0;v=0
  call eb_spray_stencil(g,phi,positions(:,1),c,w,v,ok)
  call require(ok.and.abs(sum(w)-1)<1e-14_dp.and.all(w>=0),'normalized positive stencil')
  do q=1,8
    if(w(q)<=0) cycle
    call require(c(1,q)>=2.and.v(q)>0,'covered cells excluded')
  end do
  field=3;field(:,1,:,:)=999
  call sample_eb_spray_field(g,phi,field,positions,samples,ok)
  call require(ok.and.maxval(abs(samples-3))<1e-14_dp,'constant fluid field sampled')
  call load_h2o2_full_thermo(species,ok)
  call require(ok,'thermo')
  nv=reactive_nvar(size(species))
  allocate(primitive(nv),cell(nv),y(size(species)),state(nv,4,4,4),feedback(nv,2),before(nv),after(nv))
  nitrogen=0
  do i=1,size(species)
    if(species(i)%name=='N2') nitrogen=i
  end do
  call require(nitrogen>0,'nitrogen identity')
  y=0;y(nitrogen)=1;r=mixture_specific_gas_constant(species,y,ok)
  call require(ok,'gas constant')
  primitive=0;primitive(1)=1;primitive(2:4)=[1._dp,2._dp,3._dp];primitive(5)=r*1000;primitive(6:)=y
  call reactive_primitive_to_conserved(species,primitive,cell,t,sound,ok)
  call require(ok,'initial gas')
  do k=1,4
  do j=1,4
  do i=1,4
    state(:,i,j,k)=cell
  end do
  end do
  end do
  state(:,1,:,:)=0;temperature=t;temperature(1,:,:)=42
  saved=state;saved_t=temperature
  feedback(:,1)=.01_dp*cell;feedback(:,2)=.02_dp*cell
  call inventory(before)
  call deposit_eb_spray_gas(species,g,phi,positions,feedback,state,temperature,ok)
  call require(ok,'gas deposition and EOS')
  call inventory(after)
  call require(maxval(abs(after-before-sum(feedback,dim=2))/max(1._dp,abs(before)))<1e-14_dp,'fluid-volume inventory')
  call require(all(state(:,1,:,:)==0).and.all(temperature(1,:,:)==42),'covered state untouched')
  call require(maxval(abs(temperature(2:,:,:)-t))<1e-8_dp,'recovered temperature')
  saved=state;saved_t=temperature
  feedback=-1e6_dp*feedback
  call deposit_eb_spray_gas(species,g,phi,positions,feedback,state,temperature,ok)
  call require(.not.ok.and.all(state==saved).and.all(temperature==saved_t),'nonphysical feedback rollback')
  positions(:,2)=.1_dp
  call sample_eb_spray_field(g,phi,field,positions,samples,ok)
  call require(.not.ok.and.maxval(abs(samples-3))<1e-14_dp,'solid query sampling rollback')
  print *, 'PASS: EB particle field exchange and EOS transaction'
contains
  subroutine inventory(total)
    real(dp), intent(out) :: total(:)
    integer :: ii,jj,kk
    total=0
    do kk=1,4
    do jj=1,4
    do ii=1,4
      total=total+state(:,ii,jj,kk)*g%volume_fraction(ii,jj,kk)*g%dx*g%dy*g%dz
    end do
    end do
    end do
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
