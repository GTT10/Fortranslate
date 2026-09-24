program test_reactive_boundaries_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved,reactive_conserved_to_primitive
  use reactive_3d_mod, only: advance_reactive_boundary_hydro_3d, reactive_integrals_3d
  use simulation_config_reactive_3d_mod, only: reactive_3d_config, &
    validate_reactive_3d_configuration, resolved_reactive_boundaries_3d,characteristic_pressure_ghost_3d, &
    read_reactive_3d_configuration,mach_inlet_ghost_3d
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
  call check_pressure_input()
  call check_pressure_outlet()
  call check_pressure_ghost()
  call check_mach_inlet()
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
  config%outlet_pressure(2)=102325
  call validate_reactive_3d_configuration(config,ok,message)
  call require(ok,'pressure outlet config')
  config%outlet_pressure(3)=102325
  call validate_reactive_3d_configuration(config,ok,message)
  call require(.not.ok,'pressure on wall rejected')
  config%outlet_pressure=0;config%outlet_pressure(2)=-1
  call validate_reactive_3d_configuration(config,ok,message)
  call require(.not.ok,'negative outlet pressure rejected')
  config%outlet_pressure=0
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
  subroutine check_pressure_input()
    integer :: unit
    type(reactive_3d_config) :: input_config
    open(newunit=unit,file='pressure_outlet_input_test.nml',status='replace',action='write')
    write(unit,'(a)') "&reactive_3d boundary_condition='outflow', outlet_pressure(2)=102325 /"
    close(unit)
    call read_reactive_3d_configuration('pressure_outlet_input_test.nml',input_config,ok,message)
    call require(.not.ok.and.index(message,'not connected')>0,'unconnected reader rejects pressure outlet')
    call read_reactive_3d_configuration('pressure_outlet_input_test.nml',input_config,ok,message, &
      allow_pressure_outlet=.true.)
    call require(ok.and.input_config%outlet_pressure(2)==102325,'connected reader accepts pressure outlet')
    open(newunit=unit,file='pressure_outlet_input_test.nml',status='old')
    close(unit,status='delete')
  end subroutine

  subroutine check_pressure_outlet()
    real(dp) :: outlet(6)
    integer :: face,order
    q=0;q(1)=.3_dp;q(5)=101325;q(6)=.03_dp;q(9)=.22_dp;q(12)=.75_dp
    do k=1,4
    do j=1,4
    do i=1,4
      call reactive_primitive_to_conserved(species,q,before(:,i,j,k),bt(i,j,k),sound,ok)
      call require(ok,'pressure outlet base state')
    end do
    end do
    end do
    faces='outflow'
    do order=1,2
      reconstruction='pcm'
      if(order==2) reconstruction='characteristic_plm'
      do face=1,6
        outlet=0;outlet(face)=102325;u=before;t=bt
        call reactive_integrals_3d(u,4,4,4,1._dp,1._dp,1._dp,total0,ok)
        call advance_reactive_boundary_hydro_3d(species,u,t,[1._dp,1._dp,1._dp],1e-5_dp, &
          faces,reconstruction,'mc','hllc',change,ok,outlet_pressure=outlet)
        call require(ok.and.any(change/=0),'pressure outlet PCM PLM six faces')
        call reactive_integrals_3d(u,4,4,4,1._dp,1._dp,1._dp,total1,ok)
        call require(all(abs(total1-total0-change)<=3e-12_dp*max(abs(total0),abs(total1),1._dp)), &
          'pressure outlet flux inventory')
      end do
    end do
    outlet=0;outlet(2)=1;u=before;t=bt
    call advance_reactive_boundary_hydro_3d(species,u,t,[1._dp,1._dp,1._dp],1e-5_dp, &
      faces,'pcm','mc','hllc',change,ok,outlet_pressure=outlet)
    call require(.not.ok.and.all(u==before).and.all(t==bt),'pressure ghost failure preserves fluid')
  end subroutine
  subroutine check_mach_inlet()
    real(dp) :: qi(12),ui(12),ug(12),qg(12),ti,tg,ci,cg,recovered,target_p,expected_v(3),normal
    integer :: face,axis,mode
    qi=0;qi(1)=.3_dp;qi(2:4)=[20._dp,30._dp,40._dp];qi(5)=101325
    qi(6)=.03_dp;qi(9)=.22_dp;qi(12)=.75_dp
    call reactive_primitive_to_conserved(species,qi,ui,ti,ci,ok)
    call require(ok,'Mach inlet interior')
    do face=1,6
      axis=(face+1)/2;normal=merge(1._dp,-1._dp,modulo(face,2)==1)
      do mode=1,3
        target_p=qi(5)
        select case(mode)
        case(1)
          call mach_inlet_ghost_3d(species,ui,ti,face,700._dp,.2_dp,ug,tg,ok)
        case(2)
          target_p=120000
          call mach_inlet_ghost_3d(species,ui,ti,face,700._dp,.2_dp,ug,tg,ok, &
            inlet_pressure=target_p,velocity_fluctuation=[1._dp,2._dp,3._dp])
        case(3)
          call mach_inlet_ghost_3d(species,ui,ti,face,700._dp,.2_dp,ug,tg,ok,total_velocity=[4._dp,5._dp,6._dp])
        end select
        call require(ok,'Mach inlet modes and six faces')
        call reactive_conserved_to_primitive(species,ug,tg,qg,recovered,cg,ok)
        call require(ok,'Mach inlet recovery')
        expected_v=0;expected_v(axis)=normal*.2_dp*cg
        if(mode==2) expected_v=expected_v+[1._dp,2._dp,3._dp]
        if(mode==3) expected_v=[4._dp,5._dp,6._dp]
        call require(maxval(abs(qg(2:4)-expected_v))<=3e-12_dp*maxval(abs(expected_v)),'Mach inlet velocity')
        call require(abs(qg(5)-target_p)<=3e-12_dp*target_p.and.abs(tg-700)<=3e-12_dp*700,'Mach inlet pressure and T')
        call require(all(abs(qg(6:)-qi(6:))<=3e-12_dp*max(qi(6:),1e-30_dp)),'Mach inlet frozen composition')
        call require(abs(qg(1)/qi(1)-(target_p/qi(5))*ti/700)<=3e-12_dp,'Mach inlet ideal-gas density')
      end do
    end do
    ug=42;tg=43
    call mach_inlet_ghost_3d(species,ui,ti,1,700._dp,.2_dp,ug,tg,ok,inlet_pressure=-1._dp)
    call require(.not.ok.and.all(ug==42).and.tg==43,'Mach inlet bad pressure rollback')
    call mach_inlet_ghost_3d(species,ui,ti,1,700._dp,.2_dp,ug,tg,ok, &
      velocity_fluctuation=[1._dp,2._dp,3._dp],total_velocity=[4._dp,5._dp,6._dp])
    call require(.not.ok.and.all(ug==42).and.tg==43,'Mach inlet conflicting velocity rollback')
    call mach_inlet_ghost_3d(species,ui,ti,1,-1._dp,.2_dp,ug,tg,ok)
    call require(.not.ok.and.all(ug==42).and.tg==43,'Mach inlet bad temperature rollback')
  end subroutine
  subroutine check_pressure_ghost()
    real(dp) :: qi(12),ui(12),ug(12),qg(12),ti,tg,ci,cg,expected(12),normal
    integer :: face,axis
    qi=0;qi(1)=.3_dp;qi(2:4)=[20._dp,30._dp,40._dp];qi(5)=101325
    qi(6)=.03_dp;qi(9)=.22_dp;qi(12)=.75_dp
    call reactive_primitive_to_conserved(species,qi,ui,ti,ci,ok)
    call require(ok,'pressure ghost input')
    do face=1,6
      axis=(face+1)/2;normal=merge(-1._dp,1._dp,modulo(face,2)==1)
      expected=qi;expected(1)=qi(1)+2000/(ci*ci)
      expected(axis+1)=qi(axis+1)-normal*2000/(qi(1)*ci);expected(5)=qi(5)+2000
      ug=42;tg=43
      call characteristic_pressure_ghost_3d(species,ui,ti,face,102325._dp,2._dp,ug,tg,ok)
      call require(ok,'pressure ghost six faces')
      call reactive_conserved_to_primitive(species,ug,tg,qg,ti,cg,ok)
      call require(ok.and.all(abs(qg-expected)<=3e-12_dp*max(abs(expected),1e-30_dp)), &
        'pressure density velocity and frozen composition')
      ! Restore the interior temperature guess after the output conversion.
      call reactive_primitive_to_conserved(species,qi,ui,ti,ci,ok)
    end do
    ug=42;tg=43
    call characteristic_pressure_ghost_3d(species,ui,ti,2,1._dp,2._dp,ug,tg,ok)
    call require(.not.ok.and.all(ug==42).and.tg==43,'nonpositive extrapolated pressure rollback')
  end subroutine
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
