program test_constant_transport_3d
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: load_h2o2_elementary_transport
  use mixture_thermo_mod, only: mixture_density
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_transport_3d_mod, only: reactive_transport_face_flux_3d,reactive_transport_fluxes_3d, &
    reactive_transport_timestep_3d
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
  print *, 'PASS: constant stress/heat, periodic flux scaling and timestep; face error',error
contains
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
