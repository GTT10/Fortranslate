program test_eb_redistributed_exchange_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use reactive_1d_mod, only: reactive_nvar
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d,reactive_eb_integrals_3d, &
    advance_reactive_eb_full_3d
  use eb_reactive_hydro_3d_mod, only: advance_reactive_eb_redistributed_euler_3d, &
    advance_reactive_eb_state_redistributed_euler_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(eb_geometry_3d) :: geometry
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:),next(:,:,:,:),nt(:,:,:),reference(:,:,:,:),rt(:,:,:)
  real(dp), allocatable :: before(:),after(:),boundary(:),wall(:),chemical(:)
  real(dp), parameter :: dt=1.e-8_dp,velocity=10._dp
  real(dp) :: theta
  character(len=12), parameter :: methods(2)=[character(len=12)::'flux_redist','state_redist']
  integer :: nvar,axis,mode
  logical :: ok
  character(len=1), parameter :: axes(3)=['x','y','z']
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_mechanism(reactions,ok)
  call require(ok,'reactions')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  nvar=reactive_nvar(size(species))
  allocate(u(nvar,4,4,4),t(4,4,4),next(nvar,4,4,4),nt(4,4,4),reference(nvar,4,4,4),rt(4,4,4))
  allocate(before(nvar),after(nvar),boundary(nvar),wall(nvar),chemical(nvar))
  config%nx=4;config%ny=4;config%nz=4
  do axis=1,3
    call build_axis_plane_eb_geometry_3d(4,4,4,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
      axes(axis),.39_dp,geometry,ok)
    call require(ok,'geometry')
    config%initial_velocity_x=0._dp;config%initial_velocity_y=0._dp;config%initial_velocity_z=0._dp
    select case(axis)
    case(1);config%initial_velocity_x=velocity
    case(2);config%initial_velocity_y=velocity
    case(3);config%initial_velocity_z=velocity
    end select
    call initialize_reactive_eb_density_sheet_3d(species,config,geometry,u,t,ok)
    call require(ok,'initialization')
    call reactive_eb_integrals_3d(u,geometry,before,ok)
    call require(ok,'initial inventory')
    do mode=1,2
      boundary=17._dp;wall=19._dp
      if (mode==1) then
        call advance_reactive_eb_redistributed_euler_3d(species,u,t,geometry,'rusanov',dt,next,nt,ok,boundary,wall)
      else
        call advance_reactive_eb_state_redistributed_euler_3d(species,u,t,geometry,'rusanov',dt,next,nt,ok, &
          boundary_change=boundary,wall_change=wall)
      end if
      call require(ok,'redistributed exchange step')
      call require(abs(boundary(1)+dt*config%regular_density*velocity)<1.e-20_dp,'analytic outflow mass')
      call require(abs(wall(axis+1)-dt*config%initial_pressure)<1.e-15_dp,'analytic wall pressure impulse')
      call require(wall(1)==0._dp.and.all(wall(5:)==0._dp),'impermeable adiabatic wall')
      call reactive_eb_integrals_3d(next,geometry,after,ok)
      call require(ok,'final inventory')
      call require(all(abs(after-before-boundary-wall)<=128._dp*epsilon(1._dp)*max(1._dp,abs(before))), &
        'every component balances actual face and wall exchange')
      if (mode==1) then
        call advance_reactive_eb_redistributed_euler_3d(species,u,t,geometry,'rusanov',dt,reference,rt,ok)
      else
        call advance_reactive_eb_state_redistributed_euler_3d(species,u,t,geometry,'rusanov',dt,reference,rt,ok)
      end if
      call require(ok.and.all(reference==next).and.all(rt==nt),'diagnostics preserve exact update')
      boundary=17._dp;wall=19._dp
      if (mode==1) then
        call advance_reactive_eb_redistributed_euler_3d(species,u,t,geometry,'rusanov',1.e8_dp,next,nt,ok,boundary,wall)
      else
        call advance_reactive_eb_state_redistributed_euler_3d(species,u,t,geometry,'rusanov',1.e8_dp,next,nt,ok, &
          boundary_change=boundary,wall_change=wall)
      end if
      call require(.not.ok.and.all(boundary==17._dp).and.all(wall==19._dp),'failed update retains ledgers')
      call advance_reactive_eb_full_3d(species,reactions,transport,u,t,geometry,'rusanov',methods(mode),dt, &
        .true.,2.e-7_dp,1.e-12_dp,mode==2,.true.,.true.,.true.,.true.,next,nt,theta,ok, &
        boundary_change=boundary,wall_change=wall,chemistry_change=chemical)
      call require(ok,'reacting redistributed exchange')
      call reactive_eb_integrals_3d(next,geometry,after,ok)
      call require(ok,'reacting inventory')
      call require(all(abs(after-before-boundary-wall-chemical)<=128._dp*epsilon(1._dp)*max(1._dp,abs(before))), &
        'reacting transport inventory includes every exchange')
      call require(maxval(abs(chemical(6:)))>1.e-12_dp,'chemical species contribution is nonzero')
      call advance_reactive_eb_full_3d(species,reactions,transport,u,t,geometry,'rusanov',methods(mode),dt, &
        .true.,2.e-7_dp,1.e-12_dp,mode==2,.true.,.true.,.true.,.true.,reference,rt,theta,ok)
      call require(ok.and.all(reference==next).and.all(rt==nt),'full diagnostic update remains exact')
      boundary=17._dp;wall=19._dp;chemical=23._dp
      call advance_reactive_eb_full_3d(species,reactions,transport,u,t,geometry,'invalid',methods(mode),dt, &
        .true.,2.e-7_dp,1.e-12_dp,mode==2,.true.,.true.,.true.,.true.,next,nt,theta,ok, &
        boundary_change=boundary,wall_change=wall,chemistry_change=chemical)
      call require(.not.ok.and.all(boundary==17._dp).and.all(wall==19._dp).and.all(chemical==23._dp), &
        'failure after reaction and transport preserves all ledgers')

    end do
  end do
  write(*,'(a)') 'test_eb_redistributed_exchange_3d: PASS'
contains
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not.condition) error stop message
  end subroutine require
end program test_eb_redistributed_exchange_3d
