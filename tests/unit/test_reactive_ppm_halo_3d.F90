program test_reactive_ppm_halo_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved,reactive_ppm_flattening_coefficient
  use reactive_2d_mod, only: reconstruct_reactive_characteristic_ppm_cell_2d
  use reactive_3d_mod, only: compute_reactive_ppm_ctu_fluxes_3d,compute_reactive_ppm_halo_ctu_fluxes_3d
  use reactive_3d_mod, only: advance_reactive_boundary_hydro_3d,reactive_ppm_transverse_pressure_3d
  use simulation_config_reactive_3d_mod, only: fill_reactive_boundary_ghosts_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:),g(:,:,:,:),gt(:,:,:),q(:)
  real(dp), allocatable :: fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:),hx(:,:,:,:),hy(:,:,:,:),hz(:,:,:,:)
  real(dp), allocatable :: saved(:,:,:,:),st(:,:,:),change(:),before(:)
  real(dp) :: speed,theta,halo_theta,phase,dt,mach(6),tin(6),pin(6),residual
  integer :: nv,i,j,k,mode,n,face,m,policy
  character(len=5) :: solver
  logical :: ok
  character(len=32) :: faces(6)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  nv=reactive_nvar(size(species));dt=1e-5_dp;faces='periodic'
  allocate(g(nv,12,12,12),gt(12,12,12),q(nv))
  allocate(hx(nv,0:4,4,4),hy(nv,4,0:4,4),hz(nv,4,4,0:4))
  call check_pressure_work()
  call check_transverse_flattening()
  do policy=1,2
  solver="hllc"
  if (policy==2) solver="pelec"
  do mode=1,2
    n=4
    if (mode==2) n=12
    allocate(u(nv,n,n,n),t(n,n,n),fx(nv,0:n,n,n),fy(nv,n,0:n,n),fz(nv,n,n,0:n))
    do k=1,n
    do j=1,n
    do i=1,n
      phase=2*acos(-1._dp)*real(i+2*j+3*k,dp)/n
      q=0;q(1)=1+.05_dp*sin(phase);q(2:4)=[20._dp,10._dp,-5._dp]*(1+.02_dp*sin(phase))
      q(5)=500000*(1+.01_dp*cos(phase));q(6)=.05_dp+.01_dp*cos(phase)
      if (policy==2) then
        ! Transverse compression must flatten normal pencils too.
        if (j>n/2) q(5)=2*q(5)
        q(3)=merge(40._dp,-40._dp,j<=n/2)
      end if
      q(9)=.1_dp;q(nv)=1-q(6)-q(9)
      call reactive_primitive_to_conserved(species,q,u(:,i,j,k),t(i,j,k),speed,ok)
      call require(ok,'initial state')
    end do
    end do
    end do
    call compute_reactive_ppm_ctu_fluxes_3d(species,u,t,[.25_dp,.25_dp,.25_dp],dt, &
      .true.,.true.,solver,fx,fy,fz,theta,ok,pressure_tracing=policy==2)
    call require(ok,'periodic reference')
    if (mode==1) then
      call fill_reactive_boundary_ghosts_3d(u,t,faces,g,gt,ok,ghost_width=4)
      call require(ok,'periodic four-layer ghosts')
    else
      ! A nonperiodic subdomain of a larger independent periodic calculation.
      ! The supplied donors differ from wrapping the four physical cells.
      g=u;gt=t
    end if
    call compute_reactive_ppm_halo_ctu_fluxes_3d(species,g,gt,[4,4,4],[.25_dp,.25_dp,.25_dp],dt, &
      .true.,.true.,solver,hx,hy,hz,halo_theta,ok,pressure_tracing=policy==2)
    call require(ok.and.halo_theta==1,'halo CTU')
    if (mode==1) then
      call require(all(hx==fx).and.all(hy==fy).and.all(hz==fz),'exact periodic flux parity')
    else
      call require(all(hx==fx(:,4:8,5:8,5:8)).and.all(hy==fy(:,5:8,4:8,5:8)).and. &
        all(hz==fz(:,5:8,5:8,4:8)),'exact interior subdomain flux parity')
    end if
    if (policy==2.and.mode==1) then
      call compute_reactive_ppm_ctu_fluxes_3d(species,u,t,[.25_dp,.25_dp,.25_dp],dt, &
        .true.,.true.,solver,fx,fy,fz,theta,ok)
      call require(ok.and.all(hx==fx).and.all(hy==fy).and.all(hz==fz),'default pelec uses pressure tracing')
      call compute_reactive_ppm_halo_ctu_fluxes_3d(species,g,gt,[4,4,4],[.25_dp,.25_dp,.25_dp],dt, &
        .true.,.true.,solver,hx,hy,hz,halo_theta,ok)
      call require(ok.and.all(hx==fx).and.all(hy==fy).and.all(hz==fz),'halo default uses pressure tracing')
    end if
    deallocate(u,t,fx,fy,fz)
  end do
  end do
  call compute_reactive_ppm_halo_ctu_fluxes_3d(species,g,gt,[4,4,4],[.25_dp,.25_dp,.25_dp],dt, &
    .true.,.true.,"hllc",hx,hy,hz,halo_theta,ok,pressure_tracing=.true.)
  call require(.not.ok.and.all(hx==0).and.all(hy==0).and.all(hz==0),"pressure solver mismatch rollback")
  call compute_reactive_ppm_halo_ctu_fluxes_3d(species,g,gt,[4,4,4],[.25_dp,.25_dp,.25_dp],1._dp, &
    .true.,.true.,"pelec",hx,hy,hz,halo_theta,ok,pressure_tracing=.true.)
  call require(.not.ok.and.all(hx==0).and.all(hy==0).and.all(hz==0),"pressure CFL rollback")
  ! Unused far-corner padding must not participate in EOS/CFL checks.
  g(:,1,1,1)=0;gt(1,1,1)=-1
  call compute_reactive_ppm_halo_ctu_fluxes_3d(species,g,gt,[4,4,4],[.25_dp,.25_dp,.25_dp],dt, &
    .true.,.true.,'hllc',hx,hy,hz,halo_theta,ok)
  call require(ok,'unused corner ignored')
  call compute_reactive_ppm_halo_ctu_fluxes_3d(species,g,gt,[4,4,4],[.25_dp,.25_dp,.25_dp],1._dp, &
    .true.,.true.,'hllc',hx,hy,hz,halo_theta,ok)
  call require(.not.ok.and.all(hx==0).and.all(hy==0).and.all(hz==0).and.halo_theta==0,'CFL rollback')
  allocate(u(nv,4,4,4),t(4,4,4),change(nv),before(nv))
  saved=g(:,5:8,5:8,5:8);st=gt(5:8,5:8,5:8)
  do m=1,nv
    before(m)=sum(saved(m,:,:,:))*.25_dp**3
  end do
  do face=1,6
    u=saved;t=st;faces='outflow';faces(face)='inflow'
    mach=-1;mach(face)=.1_dp;tin=0;tin(face)=400;pin=0;pin(face)=500000
    call advance_reactive_boundary_hydro_3d(species,u,t,[.25_dp,.25_dp,.25_dp],dt, &
      faces,'characteristic_ppm','mc','hllc',change,ok,inflow_temperature=tin, &
      flux_x=hx,flux_y=hy,flux_z=hz,inflow_mach=mach,inflow_pressure=pin)
    call require(ok,'six-face Mach inlet PPM step')
    do m=1,nv
      residual=abs(sum(u(m,:,:,:))*.25_dp**3-before(m)-change(m))/max(1._dp,abs(before(m)))
      call require(residual<3e-14_dp,'open boundary conserved inventory ledger')
    end do
    call require(any(u/=saved),'nontrivial inlet step')
    u=saved;t=st
    call advance_reactive_boundary_hydro_3d(species,u,t,[.25_dp,.25_dp,.25_dp],1._dp, &
      faces,'characteristic_ppm','mc','hllc',change,ok,inflow_temperature=tin, &
      flux_x=hx,flux_y=hy,flux_z=hz,inflow_mach=mach,inflow_pressure=pin)
    call require(.not.ok.and.all(u==saved).and.all(t==st).and.all(change==0).and. &
      all(hx==0).and.all(hy==0).and.all(hz==0),'boundary step rollback and unpublished fluxes')
  end do
  q=0;q(1)=1;q(5)=500000;q(nv)=1
  call reactive_primitive_to_conserved(species,q,u(:,1,1,1),t(1,1,1),speed,ok)
  call require(ok,'wall initial state')
  do k=1,4
  do j=1,4
  do i=1,4
    u(:,i,j,k)=u(:,1,1,1);t(i,j,k)=t(1,1,1)
  end do
  end do
  end do
  saved=u;st=t;faces='slip_wall'
  call advance_reactive_boundary_hydro_3d(species,u,t,[.25_dp,.25_dp,.25_dp],dt, &
    faces,'characteristic_ppm','mc','hllc',change,ok,flux_x=hx,flux_y=hy,flux_z=hz)
  call require(ok.and.all(u==saved).and.all(change==0),'stationary wall equilibrium')
  call require(all(hx(1,:,:,:)==0).and.all(hy(1,:,:,:)==0).and.all(hz(1,:,:,:)==0),'impermeable wall flux')
  print *, 'test_reactive_ppm_halo_3d: PASS' 
contains
  subroutine check_transverse_flattening()
    real(dp) :: stencil(nv,-3:3),qm(nv),qp(nv),dummy(nv),temp,c,flat
    real(dp) :: pressure(-3:3),velocity(-3:3)
    integer :: offset
    logical :: valid
    do offset=-3,3
      stencil(:,offset)=0
      stencil(1,offset)=1+.02_dp*offset
      stencil(2,offset)=20
      stencil(5,offset)=500000
      stencil(nv,offset)=1
      pressure(offset)=merge(500000._dp,1000000._dp,offset<=0)
      velocity(offset)=merge(40._dp,-40._dp,offset<=0)
    end do
    call reactive_primitive_to_conserved(species,stencil(:,0),dummy,temp,c,valid)
    call require(valid,'flattening center EOS')
    flat=reactive_ppm_flattening_coefficient(pressure,velocity)
    call require(flat==0,'transverse compressive shock detected')
    call reconstruct_reactive_characteristic_ppm_cell_2d(species,stencil,c,1e-5_dp, &
      .false.,.true.,qm,qp,valid)
    call require(valid.and.any(qm/=stencil(:,0)).and.any(qp/=stencil(:,0)),'unflattened normal gradient')
    call reconstruct_reactive_characteristic_ppm_cell_2d(species,stencil,c,1e-5_dp, &
      .false.,.true.,qm,qp,valid,flattening_override=flat)
    call require(valid.and.all(qm==stencil(:,0)).and.all(qp==stencil(:,0)),'transverse shock flattens normal states')
    call reconstruct_reactive_characteristic_ppm_cell_2d(species,stencil,c,1e-5_dp, &
      .false.,.true.,qm,qp,valid,flattening_override=-1._dp)
    call require(.not.valid.and.all(qm==0).and.all(qp==0),'invalid flattening rollback')
  end subroutine
  subroutine check_pressure_work()
    real(dp) :: normal(7),delta(7),primitive(7),e
    ! Corrected rho=2, momentum=(2,4,0), rhoE=30; therefore rhoe=25.
    normal=[3._dp,3._dp,6._dp,0._dp,40._dp,1._dp,2._dp]
    delta=[1._dp,1._dp,2._dp,0._dp,10._dp,.25_dp,.75_dp]
    call reactive_ppm_transverse_pressure_3d(normal,100._dp,delta,1.4_dp, &
      [10._dp],[20._dp],[2._dp],[3._dp],[.1_dp],primitive,e,ok)
    call require(ok.and.abs(primitive(5)-95.4_dp)<1e-13_dp,'single transverse pressure work')
    call require(e==25.and.all(primitive(1:4)==[2._dp,1._dp,2._dp,0._dp]),'independent energy and momentum')
    call require(all(primitive(6:)==[.375_dp,.625_dp]),'conserved species correction')
    call reactive_ppm_transverse_pressure_3d(normal,100._dp,delta,1.4_dp, &
      [10._dp,20._dp],[20._dp,10._dp],[2._dp,3._dp],[3._dp,2._dp], &
      [.1_dp,.2_dp],primitive,e,ok)
    call require(ok.and.abs(primitive(5)-104.6_dp)<1e-13_dp.and.e==25,'double transverse pressure work')
    delta(1)=4
    call reactive_ppm_transverse_pressure_3d(normal,100._dp,delta,1.4_dp, &
      [10._dp],[20._dp],[2._dp],[3._dp],[.1_dp],primitive,e,ok)
    call require(.not.ok.and.all(primitive==0).and.e==0,'pressure correction rollback')
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
