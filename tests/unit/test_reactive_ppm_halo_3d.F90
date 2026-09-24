program test_reactive_ppm_halo_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved
  use reactive_3d_mod, only: compute_reactive_ppm_ctu_fluxes_3d,compute_reactive_ppm_halo_ctu_fluxes_3d
  use reactive_3d_mod, only: advance_reactive_boundary_hydro_3d
  use simulation_config_reactive_3d_mod, only: fill_reactive_boundary_ghosts_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:),g(:,:,:,:),gt(:,:,:),q(:)
  real(dp), allocatable :: fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:),hx(:,:,:,:),hy(:,:,:,:),hz(:,:,:,:)
  real(dp), allocatable :: saved(:,:,:,:),st(:,:,:),change(:),before(:)
  real(dp) :: speed,theta,halo_theta,phase,dt,mach(6),tin(6),pin(6),residual
  integer :: nv,i,j,k,mode,n,face,m
  logical :: ok
  character(len=32) :: faces(6)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  nv=reactive_nvar(size(species));dt=1e-5_dp;faces='periodic'
  allocate(g(nv,12,12,12),gt(12,12,12),q(nv))
  allocate(hx(nv,0:4,4,4),hy(nv,4,0:4,4),hz(nv,4,4,0:4))
  do mode=1,2
    n=4
    if (mode==2) n=12
    allocate(u(nv,n,n,n),t(n,n,n),fx(nv,0:n,n,n),fy(nv,n,0:n,n),fz(nv,n,n,0:n))
    do k=1,n
    do j=1,n
    do i=1,n
      phase=2*acos(-1._dp)*real(i+2*j+3*k,dp)/n
      q=0;q(1)=1+.05_dp*sin(phase);q(2:4)=[20._dp,10._dp,-5._dp]
      q(5)=500000*(1+.01_dp*cos(phase));q(6)=.05_dp+.01_dp*cos(phase)
      q(9)=.1_dp;q(nv)=1-q(6)-q(9)
      call reactive_primitive_to_conserved(species,q,u(:,i,j,k),t(i,j,k),speed,ok)
      call require(ok,'initial state')
    end do
    end do
    end do
    call compute_reactive_ppm_ctu_fluxes_3d(species,u,t,[.25_dp,.25_dp,.25_dp],dt, &
      .true.,.true.,'hllc',fx,fy,fz,theta,ok)
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
      .true.,.true.,'hllc',hx,hy,hz,halo_theta,ok)
    call require(ok.and.halo_theta==1,'halo CTU')
    if (mode==1) then
      call require(all(hx==fx).and.all(hy==fy).and.all(hz==fz),'exact periodic flux parity')
    else
      call require(all(hx==fx(:,4:8,5:8,5:8)).and.all(hy==fy(:,5:8,4:8,5:8)).and. &
        all(hz==fz(:,5:8,5:8,4:8)),'exact interior subdomain flux parity')
    end if
    deallocate(u,t,fx,fy,fz)
  end do
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
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
