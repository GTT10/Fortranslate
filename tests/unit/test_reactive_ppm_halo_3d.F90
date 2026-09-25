program test_reactive_ppm_halo_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved,reactive_ppm_flattening_coefficient, &
    reconstruct_pelec_plm_cell,compute_pelec_plm_slopes,reactive_pelec_primitive_flux_x
  use reactive_2d_mod, only: reconstruct_reactive_characteristic_ppm_cell_2d
  use reactive_3d_mod, only: compute_reactive_ppm_ctu_fluxes_3d,compute_reactive_ppm_halo_ctu_fluxes_3d
  use reactive_3d_mod, only: advance_reactive_boundary_hydro_3d,reactive_ppm_transverse_pressure_3d, &
    adjust_reactive_artificial_viscosity_3d
  use simulation_config_reactive_3d_mod, only: fill_reactive_boundary_ghosts_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use, intrinsic :: ieee_exceptions, only: ieee_set_flag,ieee_get_flag,ieee_overflow
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
  call check_artificial_viscosity()
  call check_pressure_work()
  call check_flattening_denominator()
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
      call compute_reactive_ppm_ctu_fluxes_3d(species,u,t,[.25_dp,.25_dp,.25_dp],dt, &
        .true.,.true.,solver,fx,fy,fz,theta,ok,difmag=0._dp)
      call require(ok.and.(any(hx/=fx).or.any(hy/=fy).or.any(hz/=fz)),'default viscosity changes compressive flux')
      call adjust_reactive_artificial_viscosity_3d(g(:,4:9,4:9,4:9),[4,4,4],[.25_dp,.25_dp,.25_dp], &
        .1_dp,[(.false.,m=1,6)],fx,fy,fz,ok)
      call require(ok.and.all(hx==fx).and.all(hy==fy).and.all(hz==fz),'CTU then viscosity exact composition')
      call compute_reactive_ppm_ctu_fluxes_3d(species,u,t,[.25_dp,.25_dp,.25_dp],dt, &
        .true.,.true.,solver,fx,fy,fz,theta,ok,difmag=-1._dp)
      call require(.not.ok.and.all(fx==0).and.all(fy==0).and.all(fz==0),'invalid viscosity unpublished flux')
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
  call compute_reactive_ppm_halo_ctu_fluxes_3d(species,g,gt,[4,4,4],[.25_dp,.25_dp,.25_dp],dt, &
    .true.,.true.,'pelec',hx,hy,hz,halo_theta,ok)
  call require(ok,'unmarked CTU boundary reference')
  allocate(fx(nv,0:4,4,4),fy(nv,4,0:4,4),fz(nv,4,4,0:4))
  fx=hx;fy=hy;fz=hz
  call compute_reactive_ppm_halo_ctu_fluxes_3d(species,g,gt,[4,4,4],[.25_dp,.25_dp,.25_dp],dt, &
    .true.,.true.,'pelec',hx,hy,hz,halo_theta,ok,inflow_faces=[(.true.,m=1,6)])
  call require(ok,'six PLM boundary replacements')
  call require(all(hx(:,1:3,:,:)==fx(:,1:3,:,:)).and.all(hy(:,:,1:3,:)==fy(:,:,1:3,:)).and. &
    all(hz(:,:,:,1:3)==fz(:,:,:,1:3)),'interior CTU fluxes preserved')
  call require(any(hx(:,0,:,:)/=fx(:,0,:,:)).and.any(hx(:,4,:,:)/=fx(:,4,:,:)).and. &
    any(hy(:,:,0,:)/=fy(:,:,0,:)).and.any(hy(:,:,4,:)/=fy(:,:,4,:)).and. &
    any(hz(:,:,:,0)/=fz(:,:,:,0)).and.any(hz(:,:,:,4)/=fz(:,:,:,4)),'each boundary plane replaced')
  deallocate(fx,fy,fz)
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
  do policy=1,2
  solver="hllc"
  if (policy==2) solver="pelec"
  do face=1,6
    u=saved;t=st;faces='outflow';faces(face)='inflow'
    mach=-1;mach(face)=.1_dp;tin=0;tin(face)=400;pin=0;pin(face)=500000
    call advance_reactive_boundary_hydro_3d(species,u,t,[.25_dp,.25_dp,.25_dp],dt, &
      faces,'characteristic_ppm','mc',solver,change,ok,inflow_temperature=tin, &
      flux_x=hx,flux_y=hy,flux_z=hz,inflow_mach=mach,inflow_pressure=pin)
    call require(ok,'six-face Mach inlet PPM step')
    do m=1,nv
      residual=abs(sum(u(m,:,:,:))*.25_dp**3-before(m)-change(m))/max(1._dp,abs(before(m)))
      call require(residual<3e-14_dp,'open boundary conserved inventory ledger')
    end do
    call require(any(u/=saved),'nontrivial inlet step')
    u=saved;t=st
    call advance_reactive_boundary_hydro_3d(species,u,t,[.25_dp,.25_dp,.25_dp],1._dp, &
      faces,'characteristic_ppm','mc',solver,change,ok,inflow_temperature=tin, &
      flux_x=hx,flux_y=hy,flux_z=hz,inflow_mach=mach,inflow_pressure=pin)
    call require(.not.ok.and.all(u==saved).and.all(t==st).and.all(change==0).and. &
      all(hx==0).and.all(hy==0).and.all(hz==0),'boundary step rollback and unpublished fluxes')
  end do
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
  subroutine check_artificial_viscosity()
    real(dp) :: cells(7,0:3,0:3,0:3),ax(7,0:2,2,2),ay(7,2,0:2,2),az(7,2,2,0:2)
    real(dp) :: keepx(7,0:2,2,2),keepy(7,2,0:2,2),keepz(7,2,2,0:2),rho,expected_flux(7)
    integer :: x,y,z
    logical :: valid,walls(6)
    do z=0,3
    do y=0,3
    do x=0,3
      rho=1+.1_dp*x+.2_dp*y+.3_dp*z
      cells(:,x,y,z)=[rho,-rho*x,-rho*2*y,-rho*3*z,10*rho,.25_dp*rho,.75_dp*rho]
    end do
    end do
    end do
    ax=0;ay=0;az=0;walls=.false.
    call adjust_reactive_artificial_viscosity_3d(cells,[2,2,2],[1._dp,1._dp,1._dp],.1_dp,walls,ax,ay,az,valid)
    call require(valid,'compressive artificial viscosity')
    call require(maxval(abs(ax(1,:,:,:)+.06_dp))<1e-15_dp.and. &
      maxval(abs(ay(1,:,:,:)+.12_dp))<1e-15_dp.and.maxval(abs(az(1,:,:,:)+.18_dp))<1e-15_dp, &
      'analytic divergence minus six and directional density diffusion')
    expected_flux=-.6_dp*(cells(:,1,1,1)-cells(:,0,1,1))
    call require(maxval(abs(ax(:,0,1,1)-expected_flux))<2e-15_dp,'conserved jump diffusion')
    call require(maxval(abs(sum(ax(6:,:,:,:),dim=1)-ax(1,:,:,:)))<1e-16_dp,'diffusive species mass closure')
    keepx=ax;keepy=ay;keepz=az
    call check_eb_viscosity(cells,keepx,keepy,keepz)
    ax=0;ay=0;az=0;walls(1)=.true.
    call adjust_reactive_artificial_viscosity_3d(cells,[2,2,2],[1._dp,1._dp,1._dp],.1_dp,walls,ax,ay,az,valid)
    call require(valid.and.all(ax(3:4,0,:,:)==0).and.all(ax(:,1:2,:,:)==keepx(:,1:2,:,:)).and. &
      all(ay==keepy).and.all(az==keepz),'no-slip tangential momentum exclusion only at wall')
    cells(2:4,:,:,:)=-cells(2:4,:,:,:);ax=0;ay=0;az=0
    call adjust_reactive_artificial_viscosity_3d(cells,[2,2,2],[1._dp,1._dp,1._dp],.1_dp,walls,ax,ay,az,valid)
    call require(valid.and.all(ax==0).and.all(ay==0).and.all(az==0),'expansion has no artificial viscosity')
    ax=keepx;ay=keepy;az=keepz;cells(1,0,0,0)=0
    call adjust_reactive_artificial_viscosity_3d(cells,[2,2,2],[1._dp,1._dp,1._dp],.1_dp,walls,ax,ay,az,valid)
    call require(.not.valid.and.all(ax==keepx).and.all(ay==keepy).and.all(az==keepz),'viscosity rollback')
  end subroutine
  subroutine check_eb_viscosity(input,rx,ry,rz)
    real(dp), intent(in) :: input(7,0:3,0:3,0:3),rx(7,0:2,2,2),ry(7,2,0:2,2),rz(7,2,2,0:2)
    real(dp) :: cells(7,0:3,0:3,0:3),vf(0:3,0:3,0:3)
    real(dp) :: ax(7,0:2,2,2),ay(7,2,0:2,2),az(7,2,2,0:2)
    logical :: ox(0:2,2,2),oy(2,0:2,2),oz(2,2,0:2),walls(6),valid
    cells=input;vf=1;ox=.true.;oy=.true.;oz=.true.;walls=.false.
    ax=0;ay=0;az=0
    call adjust_reactive_artificial_viscosity_3d(cells,[2,2,2],[1._dp,1._dp,1._dp],.1_dp, &
      walls,ax,ay,az,valid,vf,ox,oy,oz)
    call require(valid.and.all(ax==rx).and.all(ay==ry).and.all(az==rz),'EB regular limit exact')
    vf(0,0,0)=1-4e-12_dp;ax=0;ay=0;az=0
    call adjust_reactive_artificial_viscosity_3d(cells,[2,2,2],[1._dp,1._dp,1._dp],.1_dp, &
      walls,ax,ay,az,valid,vf,ox,oy,oz)
    call require(valid.and.all(ax==rx).and.all(ay==ry).and.all(az==rz),'EB reference near-unit cutoff')
    vf(0,0,0)=.5_dp;ax=0;ay=0;az=0
    call adjust_reactive_artificial_viscosity_3d(cells,[2,2,2],[1._dp,1._dp,1._dp],.1_dp, &
      walls,ax,ay,az,valid,vf,ox,oy,oz)
    call require(valid.and.abs(ax(1,0,1,1)-.75_dp*rx(1,0,1,1))<1e-15_dp, &
      'cut cell suppresses one of four face nodes')
    call require(all(ax(:,1:,:,:)==rx(:,1:,:,:)),'distant regular faces unchanged')
    vf(0,0,0)=0;cells(:,0,0,0)=0;ax=0;ay=0;az=0
    call adjust_reactive_artificial_viscosity_3d(cells,[2,2,2],[1._dp,1._dp,1._dp],.1_dp, &
      walls,ax,ay,az,valid,vf,ox,oy,oz)
    call require(valid.and.abs(ax(1,0,1,1)-.75_dp*rx(1,0,1,1))<1e-15_dp,'covered corner state not read')
    vf=0;cells=0;ox=.false.;oy=.false.;oz=.false.;ax=rx;ay=ry;az=rz
    call adjust_reactive_artificial_viscosity_3d(cells,[2,2,2],[1._dp,1._dp,1._dp],.1_dp, &
      walls,ax,ay,az,valid,vf,ox,oy,oz)
    call require(valid.and.all(ax==0).and.all(ay==0).and.all(az==0),'closed faces zero without covered EOS')
    ox(0,1,1)=.true.;ax=rx;ay=ry;az=rz
    call adjust_reactive_artificial_viscosity_3d(cells,[2,2,2],[1._dp,1._dp,1._dp],.1_dp, &
      walls,ax,ay,az,valid,vf,ox,oy,oz)
    call require(.not.valid.and.all(ax==rx).and.all(ay==ry).and.all(az==rz),'inconsistent EB masks rollback')
  end subroutine check_eb_viscosity

  subroutine check_flattening_denominator()
    real(dp) :: pressure(-3:3),velocity(-3:3),flat
    integer :: direction,shifted,point
    logical :: overflowed
    do direction=-1,1,2
    do shifted=0,1
      pressure=1e5_dp
      if(shifted==0) then
        pressure(direction)=2e5_dp
      else
        pressure=1e5_dp*[2._dp,1._dp,1._dp,2._dp,2._dp,2._dp,1._dp]
        if(direction<0) pressure=pressure(3:-3:-1)
      end if
      do point=-3,3
        velocity(point)=-real(point,dp)
      end do
      ! The direct or shifted pressure denominator vanishes while its numerator
      ! is nonzero. Saturation must not depend on an overflowing division.
      call ieee_set_flag(ieee_overflow,.false.)
      flat=reactive_ppm_flattening_coefficient(pressure,velocity)
      call ieee_get_flag(ieee_overflow,overflowed)
      call require(flat==0.and..not.overflowed,'flattening handles zero denominator without overflow')
      velocity=-velocity
      call ieee_set_flag(ieee_overflow,.false.)
      flat=reactive_ppm_flattening_coefficient(pressure,velocity)
      call ieee_get_flag(ieee_overflow,overflowed)
      call require(flat==1.and..not.overflowed,'expansion sensor avoids zero-denominator overflow')
    end do
    end do
  end subroutine check_flattening_denominator

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
    call reconstruct_pelec_plm_cell(stencil(:,-2:2),c,0._dp,1._dp,qm,qp,valid)
    call require(valid.and.abs(qm(1)-.99_dp)<1e-14_dp.and.abs(qp(1)-1.01_dp)<1e-14_dp, &
      'fourth-order PLM linear spatial edges')
    call reconstruct_pelec_plm_cell(stencil(:,-2:2),c,1e-5_dp,0._dp,qm,qp,valid)
    call require(valid.and.all(qm==stencil(:,0)).and.all(qp==stencil(:,0)),'PLM fully flattened predictor')
    call reconstruct_pelec_plm_cell(stencil(:,-2:2),c,1._dp,1._dp,qm,qp,valid)
    call require(.not.valid.and.all(qm==0).and.all(qp==0),'PLM CFL rollback')
    block
      real(dp) :: mixture(nv,-2:2)
      mixture=stencil(:,-2:2)
      mixture(6,:)=[.1_dp,.1_dp,.2_dp,.3_dp,.3_dp]
      mixture(7,:)=[.2_dp,.25_dp,.2_dp,.1_dp,.05_dp]
      mixture(nv,:)=1-mixture(6,:)-mixture(7,:)
      call reconstruct_pelec_plm_cell(mixture,c,1e-5_dp,1._dp,qm,qp,valid)
      call require(valid.and.abs(sum(qm(6:))-1)<3e-16_dp.and.abs(sum(qp(6:))-1)<3e-16_dp, &
        'independently limited PLM species remain normalized')
    end block
    block
      real(dp) :: line(1,-2:2),gradient(1),expected,eb_stencil(nv,-2:2)
      logical :: covered(-2:2)
      integer :: mask,order,orders(3),o,m
      orders=[1,2,4]
      do mask=0,31
        do m=-2,2
          covered(m)=btest(mask,m+2)
          line(1,m)=1+.02_dp*m
          if(covered(m)) line(1,m)=ieee_value(0._dp,ieee_quiet_nan)
        end do
        do o=1,3
          order=orders(o);expected=0
          if(order>1.and..not.any(covered(-1:1))) then
            expected=.02_dp*4/3
            if(order==4.and..not.any(covered)) expected=.02_dp
          end if
          call compute_pelec_plm_slopes(line,1._dp,order,gradient,valid,covered)
          call require(valid.and.abs(gradient(1)-expected)<1e-14_dp,'EB PLM covered-stencil order fallback')
        end do
      end do
      covered=.false.;covered(-2)=.true.
      eb_stencil=stencil(:,-2:2)
      eb_stencil(:,-2)=ieee_value(0._dp,ieee_quiet_nan)
      call reconstruct_pelec_plm_cell(eb_stencil,c,0._dp,1._dp,qm,qp,valid,covered=covered)
      call require(valid.and.abs(qm(1)-(1._dp-.01_dp*4/3))<1e-14_dp.and. &
        abs(qp(1)-(1._dp+.01_dp*4/3))<1e-14_dp,'EB second-order fallback normal predictor')
      covered(-1)=.true.;eb_stencil(:,-1)=ieee_value(0._dp,ieee_quiet_nan)
      call reconstruct_pelec_plm_cell(eb_stencil,c,0._dp,1._dp,qm,qp,valid,covered=covered)
      call require(valid.and.all(qm==stencil(:,0)).and.all(qp==stencil(:,0)),'EB blocked stencil zero slope')
      call compute_pelec_plm_slopes(line,1._dp,3,gradient,valid,covered)
      call require(.not.valid.and.all(gradient==0),'invalid PLM order rejected')
    end block
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
    real(dp) :: cell(nv),faceq(nv),flux(nv),difference(nv),base(nv),restored(nv),temp,c
    integer :: directions
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
    ! The pinned passive transverse kernels preserve signed face weights.
    ! Accepted cell compositions still use the ordinary nonnegative EOS gate.
    base=0;base(1)=1;base(2)=100;base(5)=300000;base(nv)=1
    call reactive_primitive_to_conserved(species,base,cell,temp,c,ok)
    call require(ok,'signed transverse fixture')
    difference=0;difference(6)=1e-3_dp;difference(nv)=-1e-3_dp
    do directions=1,2
      call reactive_ppm_transverse_pressure_3d(cell,base(5),difference,1.4_dp, &
        spread(100._dp,1,directions),spread(100._dp,1,directions), &
        spread(0._dp,1,directions),spread(0._dp,1,directions), &
        spread(.1_dp,1,directions),faceq,e,ok)
      call require(ok,'signed transient species after transverse correction')
      call require(faceq(6)==-1e-3_dp.and.faceq(nv)==1.001_dp,'unclipped transverse species')
      call require(abs(sum(faceq(6:))-1)<1e-14_dp,'transverse species closure')
      call reactive_pelec_primitive_flux_x(species,faceq,faceq,c,flux,ok)
      call require(ok,'signed face reaches Riemann solver')
      call require(abs(sum(flux(6:))-flux(1))<1e-12_dp,'signed face mass flux closure')
      call reactive_primitive_to_conserved(species,faceq,restored,temp,c,ok)
      call require(.not.ok,'signed face is not an admissible cell composition')
      ! Refresh sound speed after the intentionally failed ordinary EOS call.
      call reactive_primitive_to_conserved(species,base,restored,temp,c,ok)
      call require(ok,'restore fixture sound speed')
    end do
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
