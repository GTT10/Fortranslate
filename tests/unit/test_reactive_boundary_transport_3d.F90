program test_reactive_boundary_transport_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species, load_h2o2_elementary_transport
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: reactive_integrals_3d
  use reactive_transport_3d_mod, only: advance_reactive_transport_3d
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
  do mode=1,4
    visc=mode==1.or.mode==4; heat=mode==2.or.mode==4; diff=mode>=3
    do b=1,2
      faces='slip_wall'
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
      if (b==1) then
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
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *,label
      error stop 1
    end if
  end subroutine
end program
