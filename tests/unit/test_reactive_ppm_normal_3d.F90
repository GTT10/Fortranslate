program test_reactive_ppm_normal_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved
  use reactive_2d_mod, only: reconstruct_reactive_characteristic_ppm_cell_2d
  use reactive_3d_mod, only: compute_reactive_ppm_normal_fluxes_3d,compute_reactive_ppm_ctu_fluxes_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:),fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:),reference(:,:),q(:),v(:),f(:)
  real(dp) :: speed,swap,dt,pi,phase,residual,theta,expected,error_normal,error_ctu,value
  integer :: nv,i,j,k,axis,mode,index(3),m
  real(dp) :: stencil(12,-3:3),qm(12),qp(12),scratch(12),flat,ratio,temperature
  integer :: status
  character(len=32) :: argument
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  if(command_argument_count()>0) then
    call get_command_argument(1,argument)
    select case(trim(argument))
    case('--thermo-stream')
      do k=1,size(species)
        write(*,'(16(es26.17e3,1x))') species(k)%molecular_weight,species(k)%temperature_mid, &
          species(k)%low_coefficients,species(k)%high_coefficients
      end do
    case('--ctu-stream')
      block
        real(dp) :: cu(12,4,4,4),ct(4,4,4),cx(12,0:4,4,4),cy(12,4,0:4,4),cz(12,4,4,0:4),cq(12)
        integer :: ii,jj,kk,flatten_flag
        real(dp) :: viscosity
        do
          read(*,*,iostat=status) dt,flatten_flag,viscosity
          if(status<0) exit
          if(status/=0) error stop 'invalid CTU interval'
          do kk=1,4
          do jj=1,4
          do ii=1,4
            read(*,*,iostat=status) cq
            if(status/=0) error stop 'invalid CTU cell'
            call reactive_primitive_to_conserved(species,cq,cu(:,ii,jj,kk),ct(ii,jj,kk),speed,ok)
            call require(ok,'CTU input EOS')
          end do
          end do
          end do
          call compute_reactive_ppm_ctu_fluxes_3d(species,cu,ct,[.25_dp,.25_dp,.25_dp],dt, &
            .false.,flatten_flag/=0,'pelec',cx,cy,cz,theta,ok,difmag=viscosity)
          call require(ok,'CTU reference fluxes')
          do kk=1,4
          do jj=1,4
          do ii=0,4
            write(*,'(12(es26.17e3,1x))') cx(:,ii,jj,kk)
          end do
          end do
          end do
          do kk=1,4
          do jj=0,4
          do ii=1,4
            write(*,'(12(es26.17e3,1x))') cy(:,ii,jj,kk)
          end do
          end do
          end do
          do kk=0,4
          do jj=1,4
          do ii=1,4
            write(*,'(12(es26.17e3,1x))') cz(:,ii,jj,kk)
          end do
          end do
          end do
        end do
      end block
    case('--reference-stream')
      do
        read(*,*,iostat=status) stencil,flat,ratio
        if(status<0) exit
        if(status/=0) error stop 'invalid reference input'
        call reactive_primitive_to_conserved(species,stencil(:,0),scratch,temperature,speed,ok)
        call require(ok,'reference center EOS')
        call reconstruct_reactive_characteristic_ppm_cell_2d(species,stencil,speed,ratio, &
          .false.,.false.,qm,qp,ok,flat,pelec_projection=.true.)
        call require(ok,'reference normal prediction')
        write(*,'(24(es26.17e3,1x))') qm,qp
      end do
    case default
      error stop 'unknown mode'
    end select
    stop
  end if
  nv=reactive_nvar(size(species));pi=acos(-1._dp);dt=1e-5_dp
  allocate(u(nv,8,8,8),t(8,8,8),fx(nv,0:8,8,8),fy(nv,8,0:8,8),fz(nv,8,8,0:8))
  allocate(reference(nv,0:8),q(nv),v(nv),f(nv))
  do mode=0,1
  do axis=1,3
    do k=1,8
    do j=1,8
    do i=1,8
      index=[i,j,k];phase=2*pi*(real(index(axis),dp)-.5_dp)/8
      q=0;q(1)=1+.1_dp*mode*sin(phase);q(2:4)=[20._dp,5._dp,-7._dp]
      q(5)=101325;q(nv)=1
      swap=q(2);q(2)=q(axis+1);q(axis+1)=swap
      call reactive_primitive_to_conserved(species,q,u(:,i,j,k),t(i,j,k),speed,ok)
      call require(ok,'initial EOS')
    end do
    end do
    end do
    call compute_reactive_ppm_normal_fluxes_3d(species,u,t,[.125_dp,.125_dp,.125_dp],dt, &
      .true.,.true.,'hllc',fx,fy,fz,ok)
    call require(ok,'normal PPM flux')
    call require(all(fx(:,0,:,:)==fx(:,8,:,:)).and.all(fy(:,:,0,:)==fy(:,:,8,:)).and. &
      all(fz(:,:,:,0)==fz(:,:,:,8)),'single periodic face identity')
    do m=1,nv
      residual=sum(fx(m,1:8,:,:)-fx(m,0:7,:,:))+sum(fy(m,:,1:8,:)-fy(m,:,0:7,:))+ &
        sum(fz(m,:,:,1:8)-fz(m,:,:,0:7))
      call require(abs(residual)<2e-13_dp*max(1._dp,maxval(abs(fx(m,:,:,:))), &
        maxval(abs(fy(m,:,:,:))),maxval(abs(fz(m,:,:,:)))),'periodic flux conservation')
    end do
    do i=0,8
      select case(axis)
      case(1);f=fx(:,i,3,3)
      case(2);f=fy(:,3,i,3)
      case(3);f=fz(:,3,3,i)
      end select
      swap=f(2);f(2)=f(axis+1);f(axis+1)=swap
      if(axis==1) reference(:,i)=f
      call require(all(abs(f-reference(:,i))<=2e-13_dp*max(1._dp,abs(reference(:,i)))), &
        'three-axis dimensional reduction with both transverse velocities')
      if (mode==0) then
        v=u(:,3,3,3);swap=v(2);v(2)=v(axis+1);v(axis+1)=swap
        v=20*v;v(2)=v(2)+101325;v(5)=v(5)+20*101325
        call require(all(abs(f-v)<=2e-13_dp*max(1._dp,abs(v))),'analytic uniform Euler flux')
      end if
    end do
    call compute_reactive_ppm_ctu_fluxes_3d(species,u,t,[.125_dp,.125_dp,.125_dp],dt, &
      .true.,.true.,'hllc',fx,fy,fz,theta,ok)
    call require(ok.and.theta==1,'CTU uniform/one-dimensional predictor')
    do i=0,8
      select case(axis)
      case(1);f=fx(:,i,3,3)
      case(2);f=fy(:,3,i,3)
      case(3);f=fz(:,3,3,i)
      end select
      swap=f(2);f(2)=f(axis+1);f(axis+1)=swap
      call require(all(abs(f-reference(:,i))<=2e-13_dp*max(1._dp,abs(reference(:,i)))), &
        'CTU reduces to normal prediction in one dimension')
    end do
  end do
  end do
  call compute_reactive_ppm_normal_fluxes_3d(species,u,t,[.125_dp,.125_dp,.125_dp],0._dp, &
    .true.,.true.,'hllc',fx,fy,fz,ok)
  f=fz(:,3,3,0);swap=f(2);f(2)=f(4);f(4)=swap
  call require(ok.and.any(f/=reference(:,0)),'profile tracing depends on time interval')
  call compute_reactive_ppm_normal_fluxes_3d(species,u,t,[.125_dp,.125_dp,.125_dp],1._dp, &
    .true.,.true.,'hllc',fx,fy,fz,ok)
  call require(.not.ok.and.all(fx==0).and.all(fy==0).and.all(fz==0),'CFL rejection clears all flux outputs')
  dt=1e-4_dp
  do k=1,8
  do j=1,8
  do i=1,8
    phase=2*pi*(real(i+j+k,dp)-1.5_dp)/8
    q=0;q(1)=1+.1_dp*sin(phase);q(2:4)=20;q(5)=101325;q(nv)=1
    call reactive_primitive_to_conserved(species,q,u(:,i,j,k),t(i,j,k),speed,ok)
    call require(ok,'diagonal contact initial state')
  end do
  end do
  end do
  do mode=0,1
    if (mode==0) then
      call compute_reactive_ppm_normal_fluxes_3d(species,u,t,[.125_dp,.125_dp,.125_dp],dt, &
        .false.,.false.,'hllc',fx,fy,fz,ok)
    else
      call compute_reactive_ppm_ctu_fluxes_3d(species,u,t,[.125_dp,.125_dp,.125_dp],dt, &
        .false.,.false.,'hllc',fx,fy,fz,theta,ok)
      call require(theta==1,'unlimited smooth CTU')
    end if
    call require(ok,'diagonal flux')
    residual=0
    do k=1,8
    do j=1,8
    do i=1,8
      expected=1+.1_dp*sin(2*pi*((real(i+j+k,dp)-1.5_dp)/8-60*dt))
      value=u(1,i,j,k)-dt/.125_dp*(fx(1,i,j,k)-fx(1,i-1,j,k)+ &
        fy(1,i,j,k)-fy(1,i,j-1,k)+fz(1,i,j,k)-fz(1,i,j,k-1))
      residual=residual+abs(value-expected)/512
    end do
    end do
    end do
    if (mode==0) error_normal=residual
    if (mode==1) error_ctu=residual
  end do
  print *, 'diagonal contact normal/CTU density L1: ',error_normal,error_ctu
  call require(error_ctu<error_normal,'corner transport improves diagonal contact prediction')
  print *, 'test_reactive_ppm_normal_3d: PASS'
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
