program test_transport_cell_reuse_3d
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use gas_transport_mod, only: gas_transport_species,set_simple_transport_fits,set_constant_transport_coefficients
  use transport_database_mod, only: load_h2o2_elementary_transport
  use mixture_thermo_mod, only: mixture_density
  use reactive_1d_mod, only: reactive_primitive_to_conserved,reactive_conserved_to_primitive
  use reactive_transport_3d_mod, only: reactive_transport_fluxes_3d,reactive_transport_ghosted_fluxes_3d, &
    reactive_transport_face_flux_3d,physical_transport_fluxes_3d
  implicit none
  integer, parameter :: nx=4,ny=3,nz=2,nv=12,ns=7
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: ordinary(:),transport(:)
  real(dp), parameter :: h(3)=[.001_dp,.0013_dp,.0009_dp]
  real(dp) :: u(nv,nx,ny,nz),t(nx,ny,nz),q(nv,nx,ny,nz),ct(nx,ny,nz)
  real(dp) :: ug(nv,0:nx+1,0:ny+1,0:nz+1),tg(0:nx+1,0:ny+1,0:nz+1)
  real(dp) :: x(nv,nx,ny,nz),y(nv,nx,ny,nz),z(nv,nx,ny,nz)
  real(dp) :: gx(nv,0:nx,ny,nz),gy(nv,nx,0:ny,nz),gz(nv,nx,ny,0:nz)
  real(dp) :: fits(4,2*ns+ns*(ns+1)/2),coef(4),theta,sound
  real(dp) :: p(nv),gradient(3,3),f(nv),energy
  integer :: model,variant,mode,i,j,k,axis,l(3),r(3),unit,checks
  logical :: ok,flags(4),dump
  character(len=1024) :: output

  checks=0;call get_command_argument(1,output);dump=len_trim(output)>0
  if(dump) open(newunit=unit,file=trim(output),access='stream',form='unformatted',status='replace')
  call load_h2o2_elementary_thermo(species,ok);call require(ok,'thermo load')
  call load_h2o2_elementary_transport(ordinary,ok);call require(ok,'transport load')
  fits=0;fits(1,:ns)=log(3.e-4_dp);fits(1,ns+1:2*ns)=log(2.e4_dp)
  fits(1,2*ns+1:)=log(.8_dp)
  do model=1,6
    transport=ordinary
    if(model==2.or.model==6) then
      call set_simple_transport_fits(transport,fits,ok);call require(ok,'synthetic Simple fits')
    endif
    if(model==5.or.model==6) then
      call set_constant_transport_coefficients(transport,[1.e-5_dp,.02_dp,.5e-5_dp,.8e-5_dp],ok)
      call require(ok,'record constant transport')
    endif
    do variant=1,3
      ! New calls change the state and then the transport data. No persistent
      ! cache may carry the preceding values or the preceding transport model.
      call build_state(variant)
      if(variant==3) then
        if(model==2.or.model==6) then
          do i=1,ns
            transport(i)%simple_viscosity(1)=transport(i)%simple_viscosity(1)+.03_dp
            transport(i)%simple_conductivity(1)=transport(i)%simple_conductivity(1)-.02_dp
          enddo
        else
          transport(1)%diameter=transport(1)%diameter*1.01_dp
        endif
      endif
      do mode=0,15
        if(btest(mode,3).and..not.btest(mode,2)) cycle
        flags=[btest(mode,0),btest(mode,1),btest(mode,2),btest(mode,3)]
        coef=[1.e-5_dp,.02_dp,0._dp,-1._dp]
        if(model==4) coef=[1.e-5_dp,.02_dp,.5e-5_dp,.8e-5_dp]
        if(model==3) then
          call check_case(coef(:2))
        else if(model==4) then
          call check_case(coef)
        else
          call check_case()
        endif
      enddo
      if(model==2) call check_physical()
    enddo
  enddo
  ! Invalid data after successful calls must still reject and reset outputs.
  transport=ordinary;call set_simple_transport_fits(transport,fits,ok)
  flags=.true.;transport(ns)%simple_binary(1,ns)=ieee_value(0._dp,ieee_quiet_nan)
  call reject_case('invalid fit after valid calls')
  transport=ordinary;u(1,nx,ny,nz)=-1
  call reject_case('invalid conserved state after valid calls')
  call build_state(1);flags=[.true.,.false.,.false.,.true.]
  call reject_case('barodiffusion without species diffusion')
  flags=.true.;tg(0,0,0)=ieee_value(0._dp,ieee_quiet_nan)
  gx=17;gy=19;gz=23
  call reactive_transport_ghosted_fluxes_3d(species,transport,ug,tg,nx,ny,nz,h(1),h(2),h(3),0._dp, &
    flags(1),flags(2),flags(3),flags(4),gx,gy,gz,theta,ok)
  call require(.not.ok.and.all(gx==0).and.all(gy==0).and.all(gz==0).and.theta==1,'invalid ghost reset')
  if(dump) close(unit)
  print *, 'PASS: cell reuse against uncached public faces; exact checks=',checks
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL: ',trim(label),' model/variant/mode=',model,variant,mode
      error stop 1
    endif
    checks=checks+1
  end subroutine

  integer function wrap(value,n) result(index)
    integer, intent(in) :: value,n
    index=modulo(value-1,n)+1
  end function

  subroutine build_state(which)
    integer, intent(in) :: which
    real(dp) :: temperature,phase
    integer :: ii,jj,kk
    do kk=1,nz
      do jj=1,ny
        do ii=1,nx
          phase=real(ii*jj+jj*kk+ii*kk,dp)
          temperature=550+3*phase+which
          p=0;p(5)=101325*(1+.002_dp*phase)
          p(2:4)=[sin(phase),cos(.2_dp*phase),.1_dp*phase]
          p(6)=.035_dp+.001_dp*ii;p(9)=.21_dp+.001_dp*jj;p(12)=1-p(6)-p(9)
          p(1)=mixture_density(species,p(6:),p(5),temperature,ok);call require(ok,'state density')
          call reactive_primitive_to_conserved(species,p,u(:,ii,jj,kk),temperature,sound,ok)
          call require(ok,'state construction');t(ii,jj,kk)=temperature
          call reactive_conserved_to_primitive(species,u(:,ii,jj,kk),t(ii,jj,kk),q(:,ii,jj,kk), &
            ct(ii,jj,kk),sound,ok)
          call require(ok,'uncached primitive recovery')
        enddo
      enddo
    enddo
    do kk=0,nz+1
      do jj=0,ny+1
        do ii=0,nx+1
          ug(:,ii,jj,kk)=u(:,wrap(ii,nx),wrap(jj,ny),wrap(kk,nz))
          tg(ii,jj,kk)=t(wrap(ii,nx),wrap(jj,ny),wrap(kk,nz))
        enddo
      enddo
    enddo
  end subroutine

  subroutine reference_face(left,right,normal,coefficients)
    integer, intent(in) :: left(3),right(3),normal
    real(dp), optional, intent(in) :: coefficients(:)
    integer :: v,dir,lp(3),lm(3),rp(3),rm(3)
    real(dp) :: dl,dr
    do v=1,3
      gradient(v,normal)=(q(v+1,right(1),right(2),right(3))- &
        q(v+1,left(1),left(2),left(3)))/h(normal)
      do dir=1,3
        if(dir==normal) cycle
        lp=left;lm=left;rp=right;rm=right
        lp(dir)=lp(dir)+1;lm(dir)=lm(dir)-1;rp(dir)=rp(dir)+1;rm(dir)=rm(dir)-1
        lp=[wrap(lp(1),nx),wrap(lp(2),ny),wrap(lp(3),nz)]
        lm=[wrap(lm(1),nx),wrap(lm(2),ny),wrap(lm(3),nz)]
        rp=[wrap(rp(1),nx),wrap(rp(2),ny),wrap(rp(3),nz)]
        rm=[wrap(rm(1),nx),wrap(rm(2),ny),wrap(rm(3),nz)]
        dl=(q(v+1,lp(1),lp(2),lp(3))-q(v+1,lm(1),lm(2),lm(3)))/(2*h(dir))
        dr=(q(v+1,rp(1),rp(2),rp(3))-q(v+1,rm(1),rm(2),rm(3)))/(2*h(dir))
        gradient(v,dir)=.5_dp*(dl+dr)
      enddo
    enddo
    f=0;energy=0
    if(.not.any(flags(:3))) return
    call reactive_transport_face_flux_3d(species,transport,q(:,left(1),left(2),left(3)), &
      q(:,right(1),right(2),right(3)),ct(left(1),left(2),left(3)),ct(right(1),right(2),right(3)), &
      h(normal),normal,gradient,flags(1),flags(2),flags(3),flags(4),f,energy,ok,coefficients)
    call require(ok,'uncached public face')
  end subroutine

  subroutine check_case(coefficients)
    real(dp), optional, intent(in) :: coefficients(:)
    call reactive_transport_fluxes_3d(species,transport,u,t,nx,ny,nz,h(1),h(2),h(3),0._dp, &
      flags(1),flags(2),flags(3),flags(4),x,y,z,theta,ok,coefficients)
    call require(ok.and.theta==1,'periodic operator')
    call reactive_transport_ghosted_fluxes_3d(species,transport,ug,tg,nx,ny,nz,h(1),h(2),h(3),0._dp, &
      flags(1),flags(2),flags(3),flags(4),gx,gy,gz,theta,ok,constant_coefficients=coefficients)
    call require(ok.and.theta==1,'ghosted operator')
    do k=1,nz
      do j=1,ny
        do i=1,nx
          l=[i,j,k]
          do axis=1,3
            r=l;r(axis)=r(axis)+1
            r=[wrap(r(1),nx),wrap(r(2),ny),wrap(r(3),nz)]
            call reference_face(l,r,axis,coefficients)
            select case(axis)
            case(1)
              call require(all(x(:,i,j,k)==f),'periodic x exact public face')
              call require(all(gx(:,i,j,k)==f),'ghost x exact public face')
              if(i==nx) call require(all(gx(:,0,j,k)==f),'ghost x lower wrap')
            case(2)
              call require(all(y(:,i,j,k)==f),'periodic y exact public face')
              call require(all(gy(:,i,j,k)==f),'ghost y exact public face')
              if(j==ny) call require(all(gy(:,i,0,k)==f),'ghost y lower wrap')
            case(3)
              call require(all(z(:,i,j,k)==f),'periodic z exact public face')
              call require(all(gz(:,i,j,k)==f),'ghost z exact public face')
              if(k==nz) call require(all(gz(:,i,j,0)==f),'ghost z lower wrap')
            end select
          enddo
        enddo
      enddo
    enddo
    if(dump) write(unit) x,y,z,gx,gy,gz
  end subroutine

  subroutine check_physical()
    character(len=32) :: faces(6)
    character(len=32), parameter :: names(8)=[character(len=32) :: 'outflow','slip walls', &
      'moving no-slip walls','isothermal walls','pressure outlet','reservoir inflow','Mach inflow','held energy']
    real(dp) :: wall_t(6),wall_v(3,6),outlet(6),inlet(nv,6),inlet_t(6),mach(6),target_p(6)
    real(dp) :: held(nx,ny,nz),saved_fit
    integer :: boundary
    flags=.true.;wall_t=650;wall_v=0
    inlet=spread(u(:,1,1,1),2,6);inlet_t=t(1,1,1)
    do boundary=1,8
      faces='outflow';outlet=0;mach=-1;target_p=0;wall_v=0
      select case(boundary)
      case(2)
        faces='slip_wall'
      case(3)
        faces='no_slip_wall';wall_v(2,1)=.3_dp
      case(4)
        faces='isothermal_wall'
      case(5)
        outlet(1)=102000
      case(6,7)
        faces(1)='inflow'
      end select
      if(boundary==7) then
        mach(1)=.2_dp;target_p(1)=103000
        call physical_transport_fluxes_3d(species,transport,u,t,[nx,ny,nz],h,0._dp,faces, &
          .true.,.true.,.true.,.true.,gx,gy,gz,theta,ok, &
          inflow_temperature=wall_t,inflow_mach=mach,inflow_pressure=target_p)
      else if(boundary==8) then
        held=u(5,:,:,:)-.5_dp*sum(u(2:4,:,:,:)**2,dim=1)/u(1,:,:,:)
        held=held+100._dp
        call physical_transport_fluxes_3d(species,transport,u,t,[nx,ny,nz],h,0._dp,faces, &
          .true.,.true.,.true.,.true.,gx,gy,gz,theta,ok,internal_energy_density=held)
      else
        call physical_transport_fluxes_3d(species,transport,u,t,[nx,ny,nz],h,0._dp,faces, &
          .true.,.true.,.true.,.true.,gx,gy,gz,theta,ok, &
          inflow_state=inlet,inflow_temperature=inlet_t,wall_temperature=wall_t, &
          wall_velocity=wall_v,outlet_pressure=outlet)
      endif
      if(.not.ok.or.theta/=1) print *, 'physical boundary case=',trim(names(boundary)),' ok=',ok,' theta=',theta
      call require(ok.and.theta==1,'Simple physical '//trim(names(boundary)))
      call require(any(gx/=0).or.any(gy/=0).or.any(gz/=0),'nonuniform physical transport active')
      if(dump) write(unit) gx,gy,gz
    enddo
    saved_fit=transport(ns)%simple_binary(1,ns)
    transport(ns)%simple_binary(1,ns)=ieee_value(0._dp,ieee_quiet_nan)
    faces='outflow'
    call physical_transport_fluxes_3d(species,transport,u,t,[nx,ny,nz],h,0._dp,faces, &
      .true.,.true.,.true.,.true.,gx,gy,gz,theta,ok)
    call require(.not.ok,'Simple physical invalid fit rejection')
    transport(ns)%simple_binary(1,ns)=saved_fit
    call physical_transport_fluxes_3d(species,transport,u,t,[nx,ny,nz],h,0._dp,faces, &
      .true.,.true.,.true.,.true.,gx,gy,gz,theta,ok)
    call require(ok,'Simple physical valid data after failed call')
    if(dump) write(unit) gx,gy,gz
  end subroutine

  subroutine reject_case(label)
    character(len=*), intent(in) :: label
    integer :: ii,jj,kk
    x=17;y=19;z=23
    call reactive_transport_fluxes_3d(species,transport,u,t,nx,ny,nz,h(1),h(2),h(3),0._dp, &
      flags(1),flags(2),flags(3),flags(4),x,y,z,theta,ok)
    call require(.not.ok.and.all(x==0).and.all(y==0).and.all(z==0).and.theta==1,label)
    do kk=0,nz+1
      do jj=0,ny+1
        do ii=0,nx+1
          ug(:,ii,jj,kk)=u(:,wrap(ii,nx),wrap(jj,ny),wrap(kk,nz))
          tg(ii,jj,kk)=t(wrap(ii,nx),wrap(jj,ny),wrap(kk,nz))
        enddo
      enddo
    enddo
    gx=17;gy=19;gz=23
    call reactive_transport_ghosted_fluxes_3d(species,transport,ug,tg,nx,ny,nz,h(1),h(2),h(3),0._dp, &
      flags(1),flags(2),flags(3),flags(4),gx,gy,gz,theta,ok)
    call require(.not.ok.and.all(gx==0).and.all(gy==0).and.all(gz==0).and.theta==1,label//' ghosted')
  end subroutine
end program
