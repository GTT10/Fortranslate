program test_reactive_ppm_step_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use elementary_kinetics_mod, only: elementary_reaction
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved
  use reactive_3d_mod, only: advance_reactive_euler_ppm_ctu_3d,advance_reactive_full_3d,advance_reactive_boundary_hydro_3d
  use reactive_3d_mod, only: compute_reactive_external_source_3d
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
  call check_source_step()
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
  subroutine check_source_step()
    real(dp) :: state(nv,4,4,4),temp(4,4,4),source(nv,4,4,4),old(nv,4,4,4),ot(4,4,4),lim
    real(dp) :: halo_source(nv,-3:8,-3:8,-3:8),boundary(nv)
    character(len=8) :: walls(6)
    real(dp) :: composed(nv,4,4,4),composed_t(4,4,4),boundary_ref(nv),transport_limit
    integer :: component,bc
    state=u;temp=t
    call compute_reactive_external_source_3d(state,[2._dp,-3._dp,4._dp],source,ok)
    call require(ok,'external force source generation')
    call require(all(source(1,:,:,:)==0).and.all(source(6:,:,:,:)==0),'force leaves mass and species unchanged')
    call require(all(source(2,:,:,:)==2).and.all(source(3,:,:,:)==-3).and.all(source(4,:,:,:)==4), &
      'force is per volume independent of density')
    call require(maxval(abs(source(5,:,:,:)+10._dp))<1e-13_dp,'external source mechanical work')
    state(1,4,4,4)=-1
    call compute_reactive_external_source_3d(state,[2._dp,-3._dp,4._dp],source,ok)
    call require(.not.ok.and.all(source==0),'external source rejects invalid density transactionally')
    ! Nonuniform fixture: periodic flux divergence must leave only source inventory.
    state=u;temp=t;source=0
    source(2,:,:,:)=.2_dp;source(5,:,:,:)=100._dp
    source(6,:,:,:)=.01_dp;source(nv,:,:,:)=-.01_dp
    call advance_reactive_euler_ppm_ctu_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt, &
      .false.,.true.,'pelec',lim,ok,conserved_source=source)
    call require(ok,'source CTU step')
    do component=1,nv
      call require(abs(sum(state(component,:,:,:))-sum(u(component,:,:,:))- &
        dt*sum(source(component,:,:,:)))<3e-14_dp*max(1._dp,abs(sum(u(component,:,:,:)))), &
        'periodic source inventory')
    end do
    call require(maxval(abs(sum(state(6:,:,:,:),dim=1)-state(1,:,:,:)))<3e-14_dp,'source species closure')
    state=u;temp=t;walls='outflow'
    do k=-3,8
    do j=-3,8
    do i=-3,8
      halo_source(:,i,j,k)=source(:,1,1,1)
    end do
    end do
    end do
    call advance_reactive_boundary_hydro_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt,walls, &
      'characteristic_ppm','mc','pelec',boundary,ok,conserved_source=halo_source)
    call require(ok,'boundary source CTU step')
    do component=1,nv
      call require(abs((sum(state(component,:,:,:))-sum(u(component,:,:,:)))*.25_dp**3- &
        boundary(component)-dt*sum(source(component,:,:,:))*.25_dp**3)< &
        3e-14_dp*max(1._dp,abs(sum(u(component,:,:,:)))*.25_dp**3),'boundary and source inventory separation')
    end do
    do bc=1,2
      walls='periodic'
      if(bc==2) walls='outflow'
      state=u;temp=t
      call advance_reactive_transport_3d(species,transport,state,temp,4,4,4,.25_dp,.25_dp,.25_dp,.5_dp*dt, &
        .true.,.true.,.true.,.true.,transport_limit,ok,walls)
      call require(ok,'source split first transport')
      if(bc==1) then
        call advance_reactive_euler_ppm_ctu_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt, &
          .false.,.true.,'pelec',lim,ok,conserved_source=source)
      else
        call advance_reactive_boundary_hydro_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt,walls, &
          'characteristic_ppm','mc','pelec',boundary_ref,ok,conserved_source=halo_source)
      end if
      call require(ok,'source split reference hydro')
      call advance_reactive_transport_3d(species,transport,state,temp,4,4,4,.25_dp,.25_dp,.25_dp,.5_dp*dt, &
        .true.,.true.,.true.,.true.,transport_limit,ok,walls)
      call require(ok,'source split second transport')
      composed=u;composed_t=t
      call advance_reactive_full_3d(species,reactions,transport,composed,composed_t,4,4,4, &
        .25_dp,.25_dp,.25_dp,dt,'pelec',.false.,1e-7_dp,1e-12_dp,.true.,.true.,.true.,.true.,.true.,lim,ok, &
        reconstruction='characteristic_ppm',boundary_faces=walls,conserved_source=halo_source)
      call require(ok.and.all(composed==state).and.all(composed_t==temp),'exact source transport split')
    end do
    state=u;temp=t;halo_source=0;halo_source(5,:,:,:)=-1e14_dp
    call advance_reactive_boundary_hydro_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt,walls, &
      'characteristic_ppm','mc','pelec',boundary,ok,conserved_source=halo_source)
    call require(.not.ok.and.all(state==u).and.all(temp==t).and.all(boundary==0),'boundary source rollback')
    call advance_reactive_full_3d(species,reactions,transport,state,temp,4,4,4, &
      .25_dp,.25_dp,.25_dp,dt,'pelec',.false.,1e-7_dp,1e-12_dp,.true.,.true.,.true.,.true.,.true.,lim,ok, &
      reconstruction='characteristic_ppm',boundary_faces=walls,conserved_source=halo_source,boundary_change=boundary)
    call require(.not.ok.and.all(state==u).and.all(temp==t).and.all(boundary==0),'source split complete rollback')
    ! Uniform fields isolate the source increment, including density injection.
    do k=1,4
    do j=1,4
    do i=1,4
      old(:,i,j,k)=u(:,1,1,1);ot(i,j,k)=t(1,1,1)
      source(:,i,j,k)=.1_dp*old(:,i,j,k)
    end do
    end do
    end do
    state=old;temp=ot
    call advance_reactive_euler_ppm_ctu_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt, &
      .false.,.true.,'pelec',lim,ok,conserved_source=source)
    call require(ok.and.maxval(abs(state-old-dt*source)/max(1._dp,abs(old)))<3e-14_dp, &
      'uniform source applied exactly once')
    state=old;temp=ot;lim=17
    source=0;source(5,:,:,:)=-1e14_dp
    call advance_reactive_euler_ppm_ctu_3d(species,state,temp,[.25_dp,.25_dp,.25_dp],dt, &
      .false.,.true.,'pelec',lim,ok,conserved_source=source)
    call require(.not.ok.and.all(state==old).and.all(temp==ot).and.lim==17,'source failure complete rollback')
  end subroutine check_source_step
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
