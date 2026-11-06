program test_constant_transport_3d
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use nasa7_thermo_mod, only: nasa7_species,nasa7_mass_properties
  use gas_transport_mod, only: gas_transport_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: load_h2o2_elementary_transport
  use mixture_thermo_mod, only: mixture_density
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_transport_3d_mod, only: reactive_transport_face_flux_3d,reactive_transport_fluxes_3d, &
    reactive_transport_timestep_3d,physical_transport_fluxes_3d,reactive_isothermal_timestep_3d
  use reactive_3d_mod, only: advance_reactive_sdc_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  real(dp) :: left(12),right(12),g(3,3),tau(3,3),flux(12),expected(12),es,coef(2),h,sound
  real(dp) :: u(12,4,4,4),t(4,4,4),x(12,4,4,4),y(12,4,4,4),z(12,4,4,4)
  real(dp) :: refx(12,4,4,4),refy(12,4,4,4),refz(12,4,4,4),theta,dt,diff,dt2,diff2,error
  logical :: ok
  integer :: i,j,k,axis
  call load_h2o2_elementary_thermo(species,ok);call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok);call require(ok,'transport')
  left=0;right=0;left(12)=1;right(12)=1;left(5)=101325;right(5)=101325
  left(1)=mixture_density(species,left(6:),left(5),500._dp,ok);call require(ok,'left density')
  right(1)=mixture_density(species,right(6:),right(5),600._dp,ok);call require(ok,'right density')
  left(2:4)=[2._dp,3._dp,4._dp];right(2:4)=[5._dp,6._dp,7._dp]
  g=reshape([1._dp,2._dp,3._dp,4._dp,5._dp,6._dp,7._dp,8._dp,9._dp],[3,3])
  coef=[1e-5_dp,.02_dp];h=.01_dp;error=0
  tau=coef(1)*(g+transpose(g))
  do i=1,3
    tau(i,i)=tau(i,i)-(2._dp/3)*coef(1)*(g(1,1)+g(2,2)+g(3,3))
  end do
  do axis=1,3
    call reactive_transport_face_flux_3d(species,transport,left,right,500._dp,600._dp,h,axis,g, &
      .true.,.true.,.false.,.false.,flux,es,ok,constant_coefficients=coef)
    call require(ok,'constant face')
    expected=0;expected(2:4)=-tau(:,axis)
    expected(5)=-dot_product(tau(:,axis),.5_dp*(left(2:4)+right(2:4)))-coef(2)*100/h
    error=max(error,maxval(abs(flux-expected)))
    call require(maxval(abs(flux-expected))<1e-12_dp.and.es==0,'analytic stress, work and heat flux')
  end do
  do k=1,4
  do j=1,4
  do i=1,4
    left(2:4)=[real(2*j,dp),real(3*k,dp),real(i,dp)]
    left(1)=mixture_density(species,left(6:),left(5),500._dp+2*i,ok)
    call require(ok,'grid density')
    call reactive_primitive_to_conserved(species,left,u(:,i,j,k),t(i,j,k),sound,ok)
    call require(ok,'grid state')
  end do
  end do
  end do
  call fluxes(coef)
  call require(ok.and.theta==1,'periodic fluxes')
  refx=x;refy=y;refz=z
  call fluxes(2*coef)
  call require(ok,'scaled fluxes')
  call require(all(x==2*refx).and.all(y==2*refy).and.all(z==2*refz),'fluxes scale with prescribed coefficients')
  call timestep(coef,dt,diff)
  call require(ok.and.diff>0,'constant timestep')
  call timestep(2*coef,dt2,diff2)
  call require(ok.and.diff2==2*diff.and.dt2==.5_dp*dt,'timestep uses the same coefficients')
  call fluxes([0._dp,0._dp])
  call require(ok.and.all(x==0).and.all(y==0).and.all(z==0),'zero coefficients disable stress and heat')
  call fluxes([-1._dp,0._dp])
  call require(.not.ok.and.all(x==0).and.all(y==0).and.all(z==0),'negative coefficient rejection')
  call timestep([ieee_value(0._dp,ieee_quiet_nan),0._dp],dt,diff)
  call require(.not.ok.and.dt==0.and.diff==0,'nonfinite coefficient rejection')
  ! An intentionally unusable molecular database must be irrelevant when the
  ! only enabled constitutive inputs are prescribed viscosity/conductivity.
  transport%diameter=0;transport%well_depth=0
  call fluxes(coef)
  call require(ok.and.all(x==refx).and.all(y==refy).and.all(z==refz),'prescribed flux needs no molecular data')
  call timestep(coef,dt2,diff2)
  call require(ok.and.diff2>0,'prescribed timestep needs no molecular data')
  call reactive_transport_fluxes_3d(species,transport,u,t,4,4,4,h,h,h,0._dp, &
    .true.,.true.,.false.,.false.,x,y,z,theta,ok)
  call require(.not.ok.and.all(x==0).and.all(y==0).and.all(z==0),'molecular path still rejects invalid database')
  call reactive_transport_fluxes_3d(species,transport,u,t,4,4,4,h,h,h,0._dp, &
    .true.,.true.,.true.,.false.,x,y,z,theta,ok,constant_coefficients=coef)
  call require(.not.ok,'species diffusion still requires molecular data')
  call check_prescribed_bulk_diffusion()
  call check_physical_constant_transport()
  u(1,1,1,1)=-1
  call fluxes(coef)
  call require(.not.ok,'prescribed transport still validates thermodynamic state')
  print *, 'PASS: constant stress/heat, periodic flux scaling and timestep; face error',error
contains
  subroutine check_prescribed_bulk_diffusion()
    real(dp) :: a(12),b(12),f(12),sf,full(4),grad(3,3),divergence
    real(dp) :: cp,cv,ha,hb,e,entropy,step,maxdiff
    integer :: axis
    a=0;b=0;a(6)=1;b(12)=1;a(5)=101325;b(5)=101325
    a(1)=mixture_density(species,a(6:),a(5),500._dp,ok);call require(ok,'constant diffusion left EOS')
    b(1)=mixture_density(species,b(6:),b(5),500._dp,ok);call require(ok,'constant diffusion right EOS')
    grad=0;full=[0._dp,0._dp,0._dp,.02_dp]
    do axis=1,3
      call reactive_transport_face_flux_3d(species,transport,a,b,500._dp,500._dp,h,axis,grad, &
        .false.,.false.,.true.,.true.,f,sf,ok,constant_coefficients=full)
      call require(ok,'prescribed rhoD bypasses invalid molecular database')
      call require(abs(f(6)-full(4)/h)<1e-12_dp.and.abs(f(12)+full(4)/h)<1e-12_dp, &
        'prescribed rhoD has no additional density factor')
      call require(sum(f(6:))==0,'prescribed diffusion mass closure')
      call nasa7_mass_properties(species(1),500._dp,cp,cv,ha,e,entropy,ok)
      call require(ok,'diffusion H2 enthalpy')
      call nasa7_mass_properties(species(7),500._dp,cp,cv,hb,e,entropy,ok)
      call require(ok,'diffusion N2 enthalpy')
      call require(abs(sf-(full(4)/h)*(ha-hb))/max(1._dp,abs(sf))<3e-14_dp, &
        'prescribed diffusion transports species enthalpy')
    enddo
    a(2:4)=[2._dp,3._dp,4._dp];b(2:4)=a(2:4)
    grad=0;grad(1,1)=1;grad(2,2)=2;grad(3,3)=3;divergence=6
    full=[0._dp,0._dp,3e-5_dp,0._dp]
    do axis=1,3
      call reactive_transport_face_flux_3d(species,transport,a,b,500._dp,500._dp,h,axis,grad, &
        .true.,.false.,.false.,.false.,f,sf,ok,constant_coefficients=full)
      call require(ok,'prescribed bulk viscosity')
      call require(abs(f(axis+1)+full(3)*divergence)<1e-18_dp,'bulk normal stress')
      call require(abs(f(5)+full(3)*divergence*a(axis+1))<1e-18_dp,'bulk viscous work')
    enddo
    full=[0._dp,0._dp,0._dp,.02_dp]
    call reactive_transport_timestep_3d(species,transport,u,t,4,4,4,h,h,h,.35_dp, &
      .false.,.false.,.true.,step,maxdiff,ok,constant_coefficients=full)
    call require(ok,'prescribed diffusion timestep with invalid molecular database')
    call require(abs(maxdiff-full(4)/minval(u(1,:,:,:)))<1e-14_dp,'rhoD timestep divides by cell density')
    full=[0._dp,0._dp,3e-5_dp,-1._dp]
    call reactive_transport_timestep_3d(species,transport,u,t,4,4,4,h,h,h,.35_dp, &
      .true.,.false.,.false.,step,maxdiff,ok,constant_coefficients=full)
    call require(ok.and.abs(maxdiff-full(3)/minval(u(1,:,:,:)))<1e-14_dp,'bulk viscosity timestep')
    call reactive_transport_face_flux_3d(species,transport,a,b,500._dp,500._dp,h,1,grad, &
      .true.,.true.,.true.,.true.,f,sf,ok,constant_coefficients=[0._dp,0._dp,0._dp])
    call require(.not.ok.and.all(f==0),'malformed coefficient vector rejected')
  end subroutine check_prescribed_bulk_diffusion
  subroutine check_physical_constant_transport()
    real(dp) :: fx(12,0:4,4,4),fy(12,4,0:4,4),fz(12,4,4,0:4)
    real(dp) :: bx(12,0:4,4,4),by(12,4,0:4,4),bz(12,4,4,0:4)
    real(dp) :: wall_t(6),wall_speed(3,6),wall_dt,wall_dt2,old_totals(12),boundary(12),source(12)
    real(dp) :: candidate(12,4,4,4),ct(4,4,4),residual
    character(len=32) :: faces(6)
    integer :: c
    faces='outflow';faces(1)='isothermal_wall'
    wall_t=600;wall_speed=0
    call physical_transport_fluxes_3d(species,transport,u,t,[4,4,4],[h,h,h],0._dp,faces, &
      .true.,.true.,.false.,.false.,fx,fy,fz,theta,ok,wall_temperature=wall_t,wall_velocity=wall_speed, &
      constant_coefficients=coef)
    call require(ok,'physical prescribed transport ignores invalid molecular database')
    call require(maxval(abs(fx(5,0,:,:)-coef(2)*(wall_t(1)-t(1,:,:))/(.5_dp*h)))<1e-12_dp, &
      'isothermal wall uses prescribed conductivity')
    bx=fx;by=fy;bz=fz
    call physical_transport_fluxes_3d(species,transport,u,t,[4,4,4],[h,h,h],0._dp,faces, &
      .true.,.true.,.false.,.false.,fx,fy,fz,theta,ok,wall_temperature=wall_t,wall_velocity=wall_speed, &
      constant_coefficients=2*coef)
    call require(ok.and.all(fx==2*bx).and.all(fy==2*by).and.all(fz==2*bz),'physical coefficient scaling')
    call reactive_isothermal_timestep_3d(species,transport,u,t,[h,h,h],faces,wall_t,.35_dp,wall_dt,ok,coef)
    call require(ok.and.wall_dt>0,'prescribed isothermal timestep')
    call reactive_isothermal_timestep_3d(species,transport,u,t,[h,h,h],faces,wall_t,.35_dp,wall_dt2,ok,2*coef)
    call require(ok.and.wall_dt2==.5_dp*wall_dt,'isothermal timestep uses prescribed conductivity')
    candidate=u;ct=t
    do c=1,12
      old_totals(c)=h**3*sum(u(c,:,:,:))
    enddo
    call advance_reactive_sdc_3d(species,candidate,ct,[h,h,h],1.e-9_dp,0._dp,2,2,ok=ok,transport=transport, &
      external_force=[0._dp,0._dp,0._dp], &
      transport_flags=[.true.,.true.,.false.,.false.],boundary_faces=faces,boundary_change=boundary, &
      source_change=source,wall_temperature=wall_t,wall_velocity=wall_speed,constant_coefficients=coef)
    call require(ok,'nonperiodic SDC prescribed transport')
    residual=0
    do c=1,12
      residual=max(residual,abs(h**3*sum(candidate(c,:,:,:))-old_totals(c)-boundary(c)-source(c))/ &
        max(1._dp,abs(old_totals(c))))
    enddo
    call require(residual<3e-14_dp,'nonperiodic SDC boundary and source inventory')
  end subroutine check_physical_constant_transport
  subroutine fluxes(coefficients)
    real(dp), intent(in) :: coefficients(2)
    call reactive_transport_fluxes_3d(species,transport,u,t,4,4,4,h,h,h,0._dp, &
      .true.,.true.,.false.,.false.,x,y,z,theta,ok,constant_coefficients=coefficients)
  end subroutine
  subroutine timestep(coefficients,step,diffusivity)
    real(dp), intent(in) :: coefficients(2)
    real(dp), intent(out) :: step,diffusivity
    call reactive_transport_timestep_3d(species,transport,u,t,4,4,4,h,h,h,.35_dp, &
      .true.,.true.,.false.,step,diffusivity,ok,constant_coefficients=coefficients)
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
