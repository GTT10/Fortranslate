program test_reactive_boundary_transport_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species, load_h2o2_elementary_transport
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: reactive_integrals_3d, advance_reactive_boundary_hydro_3d
  use reactive_transport_3d_mod, only: advance_reactive_transport_3d, reactive_transport_euler_update_3d, &
    reactive_isothermal_timestep_3d
  use mixture_thermo_mod, only: mixture_mass_properties
  use mixture_transport_mod, only: mixture_transport_coefficients
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  real(dp) :: u(12,4,4,4),t(4,4,4),before(12,4,4,4),bt(4,4,4)
  real(dp) :: q(12),sound,change(12),total0(12),total1(12),theta
  real(dp), parameter :: dx=0.001_dp,dt=1e-6_dp
  character(len=32) :: faces(6)
  logical :: ok,visc,heat,diff
  integer :: i,j,k,mode,b
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  call check_wall_shear()
  call check_wall_heat()
  do mode=1,4
    visc=mode==1.or.mode==4; heat=mode==2.or.mode==4; diff=mode>=3
    do b=1,3
      faces='slip_wall'
      if (b==3) faces='no_slip_wall'
      if (b==2) then
        faces(1:4)='outflow'; faces(5:6)='periodic'
      end if
      do k=1,4
        do j=1,4
          do i=1,4
            q=0.0_dp; q(1)=0.3_dp*(1+0.02_dp*(i+j+k)); q(5)=101325.0_dp*(1+0.02_dp*k)
            q(2:4)=[10.0_dp+i,20.0_dp+j,-4.0_dp+0.5_dp*k]
            q(6)=0.025_dp+0.003_dp*i; q(9)=0.22_dp; q(12)=1-q(6)-q(9)
            call reactive_primitive_to_conserved(species,q,u(:,i,j,k),t(i,j,k),sound,ok)
            call require(ok,'initial mixture')
          end do
        end do
      end do
      before=u; bt=t
      call reactive_integrals_3d(u,4,4,4,dx,dx,dx,total0,ok)
      call require(ok,'initial integrals')
      call advance_reactive_transport_3d(species,transport,u,t,4,4,4,dx,dx,dx,dt, &
        visc,heat,diff,diff,theta,ok,faces,change)
      call require(ok,'physical transport')
      call require(any(u/=before),'active transport')
      call require(all(u(6:,:,:,:)>=0),'positive species')
      call require(maxval(abs(sum(u(6:,:,:,:),dim=1)-u(1,:,:,:)))<1e-12_dp,'species closure')
      call reactive_integrals_3d(u,4,4,4,dx,dx,dx,total1,ok)
      call require(ok,'final integrals')
      call require(maxval(abs(total1-total0-change)/max(abs(total0),1e-12_dp))<2e-12_dp,'surface budget')
      if (b/=2) then
        call require(maxval(abs(change([1,5,6,7,8,9,10,11,12])))<1e-18_dp,'adiabatic impermeable wall')
      else if (visc) then
        write(*,'(a,i0,a,es14.6)') 'mode=',mode,' open viscous energy=',change(5)
        call require(abs(change(5))>1e-12_dp,'open viscous work accounted')
      end if
      if (mode==3.and.b==1) then
        call advance_reactive_transport_3d(species,transport,u,t,4,4,4,dx,dx,dx,0.01_dp, &
          .false.,.false.,.true.,.true.,theta,ok,faces,change)
        call require(ok.and.theta<1.0_dp,'active physical species limiter')
        call require(all(u(6:,:,:,:)>=0),'limited species remain nonnegative')
        call require(maxval(abs(sum(u(6:,:,:,:),dim=1)-u(1,:,:,:)))<1e-12_dp,'limited species closure')
      end if
      before=u; bt=t
      call advance_reactive_transport_3d(species,transport,u,t,4,4,4,dx,dx,dx,-dt, &
        visc,heat,diff,diff,theta,ok,faces,change)
      call require(.not.ok.and.all(u==before).and.all(t==bt).and.all(change==0),'rollback')
    end do
  end do
  print *, 'PASS: physical viscosity/conduction/diffusion, surface ledger, walls and rollback'
contains
  subroutine check_wall_heat()
    real(dp) :: mu,kappa,diffusion(7),wt(6),expected,wall_dt,rate,mw,gas_r,cp,cv,gamma,h,e,entropy
    integer :: face,thermal_case
    q=0.0_dp; q(1)=0.3_dp; q(5)=101325.0_dp
    q(6)=0.03_dp; q(9)=0.22_dp; q(12)=0.75_dp
    do k=1,4
      do j=1,4
        do i=1,4
          call reactive_primitive_to_conserved(species,q,before(:,i,j,k),bt(i,j,k),sound,ok)
          call require(ok,'uniform heat-transfer state')
        end do
      end do
    end do
    do thermal_case=1,2
      wt=bt(1,1,1)*0.5_dp
      if (thermal_case==2) wt=bt(1,1,1)*1.5_dp
      call mixture_transport_coefficients(species,transport,q(6:),wt(1),q(5),mu,kappa,diffusion,ok)
      call require(ok,'wall conductivity')
      call mixture_mass_properties(species,q(6:),bt(1,1,1),mw,gas_r,cp,cv,gamma,h,e,entropy,ok)
      call require(ok,'cell thermal capacity')
      expected=2*kappa*(wt(1)-bt(1,1,1))/dx*(4*dx)**2*dt
      rate=2*kappa/(q(1)*cv*dx**2)
      do face=1,6
        faces='no_slip_wall'; faces(face)='isothermal_wall'
        call reactive_isothermal_timestep_3d(species,transport,before,bt,[dx,dx,dx],faces,wt,.35_dp,wall_dt,ok)
        call require(ok.and.abs(wall_dt/(.35_dp/rate)-1)<1e-12_dp,'wall thermal timestep')
        call reactive_transport_euler_update_3d(species,transport,before,bt,4,4,4,dx,dx,dx,dt, &
          .false.,.true.,.false.,.false.,u,t,theta,ok,faces,change,wall_temperature=wt)
        call require(ok,'wall heat Euler update')
        call require(abs(change(5)/expected-1)<1e-12_dp,'analytical Fourier wall flux')
        call require(maxval(abs(change([1,2,3,4,6,7,8,9,10,11,12])))<1e-18_dp,'heat-only boundary exchange')
        call require(abs((sum(u(5,:,:,:))-sum(before(5,:,:,:)))*dx**3/change(5)-1)<1e-9_dp, &
          'heat enters gas energy')
      end do
    end do
    u=before; t=bt; wt=-1
    call advance_reactive_transport_3d(species,transport,u,t,4,4,4,dx,dx,dx,dt, &
      .true.,.true.,.true.,.true.,theta,ok,faces,change,wall_temperature=wt)
    call require(.not.ok.and.all(u==before).and.all(t==bt).and.all(change==0),'invalid wall rollback')
    print *, 'PASS: six-face hot/cold Fourier heat flux, thermal timestep and rejection'
  end subroutine check_wall_heat
  subroutine check_wall_shear()
    real(dp) :: mu,kappa,diffusion(7),expected,kinetic0,kinetic1
    ! Uniform tangential flow between stationary y walls. The face derivative
    ! is U/(dy/2). Both walls exert -mu*2U/dy over area (4dx)^2.
    q=0.0_dp; q(1)=0.3_dp; q(2)=20.0_dp; q(5)=101325.0_dp
    q(6)=0.03_dp; q(9)=0.22_dp; q(12)=0.75_dp
    do k=1,4
      do j=1,4
        do i=1,4
          call reactive_primitive_to_conserved(species,q,u(:,i,j,k),t(i,j,k),sound,ok)
          call require(ok,'uniform wall shear state')
        end do
      end do
    end do
    call mixture_transport_coefficients(species,transport,q(6:),t(1,1,1),q(5),mu,kappa,diffusion,ok)
    call require(ok,'uniform viscosity')
    faces='periodic'; faces(3:4)='no_slip_wall'
    expected=-4.0_dp*mu*q(2)/dx*(4*dx)**2*dt
    before=u; bt=t
    call advance_reactive_boundary_hydro_3d(species,u,t,[dx,dx,dx],dt,faces,'pcm','mc','rusanov',change,ok)
    call require(ok,'impermeable hydro step')
    call require(maxval(abs(u-before)/max(1.0_dp,abs(before)))<1e-13_dp,'no artificial hydro wall drag')
    call reactive_transport_euler_update_3d(species,transport,before,bt,4,4,4,dx,dx,dx,dt, &
      .true.,.false.,.false.,.false.,u,t,theta,ok,faces,change)
    call require(ok,'wall shear Euler step')
    call require(abs(change(2)/expected-1.0_dp)<1e-12_dp,'analytical wall shear')
    call require(abs(change(5))<1e-18_dp,'stationary wall does no work')
    kinetic0=sum(0.5_dp*sum(before(2:4,:,:,:)**2,dim=1)/before(1,:,:,:))
    kinetic1=sum(0.5_dp*sum(u(2:4,:,:,:)**2,dim=1)/u(1,:,:,:))
    call require(kinetic1<kinetic0.and.maxval(t)>maxval(bt),'viscous heating')
    call require(abs(sum(u(5,:,:,:))-sum(before(5,:,:,:)))<1e-9_dp,'closed total energy')
    print *, 'PASS: analytical no-slip shear and kinetic-to-internal energy conversion'
  end subroutine check_wall_shear
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *,label
      error stop 1
    end if
  end subroutine
end program
