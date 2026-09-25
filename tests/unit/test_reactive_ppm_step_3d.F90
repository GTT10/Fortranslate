program test_reactive_ppm_step_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use elementary_kinetics_mod, only: elementary_reaction
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved
  use reactive_3d_mod, only: advance_reactive_euler_ppm_ctu_3d,advance_reactive_full_3d,advance_reactive_boundary_hydro_3d
  use reactive_transport_3d_mod, only: advance_reactive_transport_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:),saved(:,:,:,:),st(:,:,:),reference(:,:,:,:),rt(:,:,:),q(:),before(:),change(:)
  real(dp) :: theta,stage_theta,speed,phase,dt,residual,after
  integer :: nv,i,j,k,m,step,policy
  character(len=5) :: solver
  logical :: ok
  character(len=8) :: faces(6)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  allocate(reactions(0));nv=reactive_nvar(size(species))
  allocate(u(nv,4,4,4),t(4,4,4),q(nv),before(nv),change(nv));dt=1e-5_dp
  do k=1,4
  do j=1,4
  do i=1,4
    phase=2*acos(-1._dp)*real(i+j+k,dp)/4
    q=0;q(1)=1+.05_dp*sin(phase);q(2:4)=[20._dp,10._dp,-5._dp];q(5)=500000
    q(6)=.05_dp+.01_dp*cos(phase);q(9)=.1_dp+.007_dp*sin(phase);q(nv)=1-q(6)-q(9)
    call reactive_primitive_to_conserved(species,q,u(:,i,j,k),t(i,j,k),speed,ok)
    call require(ok,'initial state')
  end do
  end do
  end do
  saved=u;st=t
  do m=1,nv
    before(m)=sum(u(m,:,:,:))
  end do
  do policy=1,2
  solver="hllc"
  if (policy==2) solver="pelec"
  u=saved;t=st
  do step=1,4
    call advance_reactive_euler_ppm_ctu_3d(species,u,t,[.25_dp,.25_dp,.25_dp],dt, &
      .false.,.true.,solver,theta,ok,difmag=merge(0._dp,.1_dp,mod(step,2)==0))
    call require(ok.and.theta>=0.and.theta<=1.and.all(t>0),'accepted CTU step')
    do m=1,nv
      after=sum(u(m,:,:,:));residual=abs(after-before(m))/max(1._dp,abs(before(m)))
      call require(residual<3e-14_dp,'multistep conserved inventory')
    end do
    call require(maxval(abs(sum(u(6:,:,:,:),dim=1)-u(1,:,:,:)))<3e-14_dp,'species mass closure')
  end do
  call require(any(u/=saved),'nontrivial time evolution')
  u=saved;t=st;theta=17
  call advance_reactive_euler_ppm_ctu_3d(species,u,t,[.25_dp,.25_dp,.25_dp],1._dp, &
    .false.,.true.,solver,theta,ok)
  call require(.not.ok.and.all(u==saved).and.all(t==st).and.theta==17,'excessive dt complete rollback')
  call advance_reactive_euler_ppm_ctu_3d(species,u,t,[.25_dp,.25_dp,.25_dp],0._dp, &
    .false.,.true.,solver,theta,ok)
  call require(ok.and.all(u==saved).and.all(t==st).and.theta==1,'zero step exact identity')
  reference=saved;rt=st
  call advance_reactive_transport_3d(species,transport,reference,rt,4,4,4,.25_dp,.25_dp,.25_dp,.5_dp*dt, &
    .true.,.true.,.true.,.true.,stage_theta,ok)
  call require(ok,'first transport half')
  call advance_reactive_euler_ppm_ctu_3d(species,reference,rt,[.25_dp,.25_dp,.25_dp],dt, &
    .false.,.true.,solver,theta,ok)
  call require(ok,'reference PPM hydro')
  call advance_reactive_transport_3d(species,transport,reference,rt,4,4,4,.25_dp,.25_dp,.25_dp,.5_dp*dt, &
    .true.,.true.,.true.,.true.,stage_theta,ok)
  call require(ok,'second transport half')
  call advance_reactive_full_3d(species,reactions,transport,u,t,4,4,4,.25_dp,.25_dp,.25_dp,dt,solver, &
    .false.,1e-7_dp,1e-12_dp,.true.,.true.,.true.,.true.,.true.,theta,ok,reconstruction='characteristic_ppm')
  call require(ok.and.all(u==reference).and.all(t==rt),'exact transport PPM split dispatch')
  u=saved;t=st;faces='outflow'
  reference=saved;rt=st
  call advance_reactive_transport_3d(species,transport,reference,rt,4,4,4,.25_dp,.25_dp,.25_dp,.5_dp*dt, &
    .true.,.true.,.true.,.true.,stage_theta,ok,faces)
  call require(ok,'physical first transport half')
  call advance_reactive_boundary_hydro_3d(species,reference,rt,[.25_dp,.25_dp,.25_dp],dt, &
    faces,'characteristic_ppm','mc',solver,change,ok)
  call require(ok,'physical PPM hydro')
  call advance_reactive_transport_3d(species,transport,reference,rt,4,4,4,.25_dp,.25_dp,.25_dp,.5_dp*dt, &
    .true.,.true.,.true.,.true.,stage_theta,ok,faces)
  call require(ok,'physical second transport half')
  call advance_reactive_full_3d(species,reactions,transport,u,t,4,4,4,.25_dp,.25_dp,.25_dp,dt,solver, &
    .false.,1e-7_dp,1e-12_dp,.true.,.true.,.true.,.true.,.true.,theta,ok, &
    reconstruction='characteristic_ppm',boundary_faces=faces)
  call require(ok.and.all(u==reference).and.all(t==rt),'exact physical transport PPM split dispatch')
  end do
  print *, 'test_reactive_ppm_step_3d: PASS'
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
