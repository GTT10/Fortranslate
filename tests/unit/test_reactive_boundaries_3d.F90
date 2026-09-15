program test_reactive_boundaries_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: advance_reactive_boundary_hydro_3d, reactive_integrals_3d
  use simulation_config_reactive_3d_mod, only: reactive_3d_config, &
    validate_reactive_3d_configuration, resolved_reactive_boundaries_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(reactive_3d_config) :: config
  real(dp) :: u(12,4,4,4),t(4,4,4),before(12,4,4,4),bt(4,4,4)
  real(dp) :: q(12),sound,change(12),total0(12),total1(12)
  character(len=32) :: faces(6),reconstruction,solver
  character(len=32), parameter :: solvers(3)=[character(len=32)::'rusanov','hllc','pelec']
  character(len=256) :: message
  logical :: ok
  integer :: i,j,k,r,b,s
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call check_inflow()
  config%boundary_condition='isothermal_wall'
  call validate_reactive_3d_configuration(config,ok,message)
  call require(.not.ok,'isothermal wall requires conduction')
  config%transport_enabled=.true.
  call validate_reactive_3d_configuration(config,ok,message)
  call require(ok,'valid isothermal config')
  config%wall_temperature(1)=-1
  call validate_reactive_3d_configuration(config,ok,message)
  call require(.not.ok,'invalid wall temperature')
  config=reactive_3d_config()
  config%boundary_condition='outflow'
  config%boundary_faces(3:6)='slip_wall'
  call validate_reactive_3d_configuration(config,ok,message)
  call require(ok,'mixed physical config')
  faces=resolved_reactive_boundaries_3d(config)
  call require(all(faces(1:2)=='outflow'),'face defaults')
  config%boundary_faces(1)='periodic'
  call validate_reactive_3d_configuration(config,ok,message)
  call require(.not.ok,'unpaired periodic config')
  config%boundary_faces=''; config%transport_enabled=.true.
  call validate_reactive_3d_configuration(config,ok,message)
  call require(ok,'physical transport config')
  do r=1,2
    reconstruction='pcm'
    if (r==2) reconstruction='characteristic_plm'
    do s=1,3
      solver=solvers(s)
      do b=1,4
        faces='outflow'
        if (b==2) faces='slip_wall'
        if (b==3) faces(3:6)='periodic'
        if (b==4) faces='no_slip_wall'
        do k=1,4
          do j=1,4
            do i=1,4
              q=0.0_dp; q(1)=0.3_dp*(1.0_dp+0.02_dp*(i+j+k)); q(5)=101325.0_dp
              q(2:4)=[12.0_dp+i,-8.0_dp+j,4.0_dp-k]
              q(6)=0.03_dp; q(9)=0.22_dp; q(12)=0.75_dp
              call reactive_primitive_to_conserved(species,q,u(:,i,j,k),t(i,j,k),sound,ok)
              call require(ok,'initial mixture')
            end do
          end do
        end do
        before=u; bt=t
        call reactive_integrals_3d(u,4,4,4,1.0_dp,1.0_dp,1.0_dp,total0,ok)
        call require(ok,'initial integrals')
        call advance_reactive_boundary_hydro_3d(species,u,t,[1.0_dp,1.0_dp,1.0_dp],1e-5_dp, &
          faces,reconstruction,'mc',solver,change,ok)
        call require(ok,'physical hydro step')
        call reactive_integrals_3d(u,4,4,4,1.0_dp,1.0_dp,1.0_dp,total1,ok)
        call require(ok,'final integrals')
        call require(maxval(abs(total1-total0-change)/max(1.0_dp,abs(total0)))<2e-12_dp,'surface budget')
        call require(any(u/=before),'nontrivial evolution')
        if (b==2.or.b==4) then
          call require(maxval(abs(change([1,5,6,7,8,9,10,11,12])))<1e-10_dp,'impermeable adiabatic walls')
          call require(maxval(abs(change(2:4)))>1e-5_dp,'wall pressure changes momentum')
        end if
        before=u; bt=t; faces(1)='invalid'
        call advance_reactive_boundary_hydro_3d(species,u,t,[1.0_dp,1.0_dp,1.0_dp],1e-5_dp, &
          faces,reconstruction,'mc',solver,change,ok)
        call require(.not.ok.and.all(u==before).and.all(t==bt).and.all(change==0.0_dp),'failure rollback')
      end do
    end do
  end do
  print *, 'PASS: six-face PCM/PLM, three solvers, surface balance, wall flux and rollback'
contains
  subroutine check_inflow()
    real(dp) :: inlet(12,6),it(6),base(12,4,4,4),base_t(4,4,4)
    integer :: face,axis
    q=0.0_dp; q(1)=0.3_dp; q(5)=101325.0_dp
    q(6)=0.03_dp; q(9)=0.22_dp; q(12)=0.75_dp
    do k=1,4
      do j=1,4
        do i=1,4
          call reactive_primitive_to_conserved(species,q,base(:,i,j,k),base_t(i,j,k),sound,ok)
          call require(ok,'base inlet state')
        end do
      end do
    end do
    do face=1,6
      axis=(face+1)/2
      faces='outflow'; faces(face)='inflow'
      u=base; t=base_t
      inlet=0.0_dp; it=0.0_dp; q(2:4)=0.0_dp
      q(axis+1)=50.0_dp
      if (modulo(face,2)==0) q(axis+1)=-50.0_dp
      call reactive_primitive_to_conserved(species,q,inlet(:,face),it(face),sound,ok)
      call require(ok,'prescribed inlet state')
      call reactive_integrals_3d(u,4,4,4,1.0_dp,1.0_dp,1.0_dp,total0,ok)
      call require(ok,'inlet initial integrals')
      call advance_reactive_boundary_hydro_3d(species,u,t,[1.0_dp,1.0_dp,1.0_dp],1e-5_dp, &
        faces,'characteristic_plm','mc','hllc',change,ok,inlet,it)
      call require(ok,'six-face inlet step')
      call require(change(1)>0,'inward inlet mass')
      call require(change(axis+1)*q(axis+1)>0,'inward inlet momentum')
      call reactive_integrals_3d(u,4,4,4,1.0_dp,1.0_dp,1.0_dp,total1,ok)
      call require(ok,'inlet final integrals')
      call require(maxval(abs(total1-total0-change)/max(abs(total0),1.0_dp))<2e-12_dp,'inlet budget')
      before=u; bt=t
      call advance_reactive_boundary_hydro_3d(species,u,t,[1.0_dp,1.0_dp,1.0_dp],1e-5_dp, &
        faces,'characteristic_plm','mc','hllc',change,ok)
      call require(.not.ok.and.all(u==before).and.all(t==bt).and.all(change==0),'missing inlet rollback')
    end do
    config=reactive_3d_config(); config%boundary_faces(1)='inflow'
    config%boundary_faces(2)='outflow'; config%inflow_temperature(1)=-1
    call validate_reactive_3d_configuration(config,ok,message)
    call require(.not.ok,'invalid inlet temperature')
    config=reactive_3d_config()
  end subroutine check_inflow
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *,label
      error stop 1
    end if
  end subroutine
end program
