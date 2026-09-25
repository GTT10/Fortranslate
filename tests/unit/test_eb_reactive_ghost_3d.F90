program test_eb_reactive_ghost_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar,reactive_riemann_flux_x
  use reactive_2d_mod, only: reactive_riemann_flux_y
  use reactive_directional_flux_3d_mod, only: reactive_riemann_flux_z
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d
  use eb_reactive_hydro_3d_mod, only: reactive_eb_riemann_fluxes_3d,advance_reactive_eb_euler_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_reactive_wall_flux_3d_mod, only: reactive_eb_slip_wall_source_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(eb_geometry_3d) :: geometry
  type(reactive_eb_3d_config) :: config
  type(eb_integrated_flux_3d) :: base,flux,integral,saved_flux
  real(dp), allocatable :: state(:,:,:,:),temp(:,:,:),ghost(:,:,:,:),gt(:,:,:),next(:,:,:,:),nt(:,:,:)
  real(dp), allocatable :: expected(:),left(:),right(:),wall(:,:,:,:)
  real(dp) :: tl,tr,budget,change,volume,position
  real(dp), parameter :: dt=1.e-9_dp
  integer :: nvar,i,j,k,axis,side,idx(3),inside(3),outside(3),mode,solver_id,c
  character(len=7), parameter :: solvers(3)=[character(len=7)::'rusanov','hllc','pelec']
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo load')
  nvar=reactive_nvar(size(species));config%nx=4;config%ny=4;config%nz=4
  allocate(state(nvar,4,4,4),temp(4,4,4),next(nvar,4,4,4),nt(4,4,4),wall(nvar,4,4,4))
  allocate(ghost(nvar,0:5,0:5,0:5),gt(0:5,0:5,0:5),expected(nvar),left(nvar),right(nvar))
  allocate(base%x(nvar,0:4,4,4),base%y(nvar,4,0:4,4),base%z(nvar,4,4,0:4))
  base%x=0._dp;base%y=0._dp;base%z=0._dp
  flux=base
  do mode=1,2
    position=-.1_dp
    if(mode==2) position=.37_dp
    call build_axis_plane_eb_geometry_3d(4,4,4,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',position,geometry,ok)
    call require(ok,'geometry')
    call initialize_reactive_eb_density_sheet_3d(species,config,geometry,state,temp,ok)
    call require(ok,'initial state')
    do solver_id=1,3
      ghost=ieee_value(0._dp,ieee_quiet_nan);gt=ieee_value(0._dp,ieee_quiet_nan)
      do k=0,5
        do j=0,5
          do i=0,5
            idx=[i,j,k]
            if(count(idx==0.or.idx==5)/=1) cycle
            inside=max(1,min(4,idx))
            ghost(:,i,j,k)=state(:,inside(1),inside(2),inside(3))
            gt(i,j,k)=temp(inside(1),inside(2),inside(3))
          end do
        end do
      end do
      if (mode==2) then
        ghost(:,0,:,:)=ieee_value(0._dp,ieee_quiet_nan)
        gt(0,:,:)=ieee_value(0._dp,ieee_quiet_nan)
      end if
      call reactive_eb_riemann_fluxes_3d(species,state,temp,geometry,solvers(solver_id),&
        base%x,base%y,base%z,ok,difmag=0._dp)
      call require(ok,'legacy outflow flux')
      call reactive_eb_riemann_fluxes_3d(species,state,temp,geometry,solvers(solver_id),flux%x,flux%y,flux%z, &
        ok,ghost,gt,difmag=0._dp)
      call require(ok.and.all(flux%x==base%x).and.all(flux%y==base%y).and.all(flux%z==base%z), &
        'outflow halo exact legacy parity including unused NaNs')
      ! The default positive viscosity requires volume fractions and complete halos.
      call reactive_eb_riemann_fluxes_3d(species,state,temp,geometry,solvers(solver_id),flux%x,flux%y,flux%z, &
        ok,ghost,gt)
      call require(.not.ok.and.all(flux%x==0._dp).and.all(flux%y==0._dp).and.all(flux%z==0._dp), &
        'default viscosity rejects missing volume fraction halo')
      ghost=1.1_dp*ghost
      call reactive_eb_riemann_fluxes_3d(species,state,temp,geometry,solvers(solver_id),flux%x,flux%y,flux%z, &
        ok,ghost,gt,difmag=0._dp)
      call require(ok,'nonuniform boundary donor flux')
      do axis=1,3
        do side=0,1
          idx=2;idx(axis)=4*side
          inside=idx;inside(axis)=1+3*side
          outside=idx;outside(axis)=5*side
          if(axis==1.and.side==0.and.mode==2) cycle
          left=state(:,inside(1),inside(2),inside(3));right=left
          tl=temp(inside(1),inside(2),inside(3));tr=tl
          if(side==0) then
            left=ghost(:,outside(1),outside(2),outside(3));tl=gt(outside(1),outside(2),outside(3))
          else
            right=ghost(:,outside(1),outside(2),outside(3));tr=gt(outside(1),outside(2),outside(3))
          end if
          select case(axis)
          case(1)
            call reactive_riemann_flux_x(species,left,right,tl,tr,solvers(solver_id),expected,ok)
            call require(ok.and.all(expected==flux%x(:,idx(1),idx(2),idx(3))),'x donor orientation')
          case(2)
            call reactive_riemann_flux_y(species,left,right,tl,tr,solvers(solver_id),expected,ok)
            call require(ok.and.all(expected==flux%y(:,idx(1),idx(2),idx(3))),'y donor orientation')
          case(3)
            call reactive_riemann_flux_z(species,left,right,tl,tr,solvers(solver_id),expected,ok)
            call require(ok.and.all(expected==flux%z(:,idx(1),idx(2),idx(3))),'z donor orientation')
          end select
        end do
      end do
      call advance_reactive_eb_euler_3d(species,state,temp,geometry,solvers(solver_id),&
        dt,next,nt,ok,ghost,gt,integral,difmag=0._dp)
      call require(ok.and.all(nt>0._dp),'halo Euler update')
      call require(all(integral%x==dt*flux%x).and.all(integral%y==dt*flux%y).and. &
        all(integral%z==dt*flux%z),'accepted integrated flux matches update')
      call require(all(lbound(integral%x)==[1,0,1,1]),'integrated flux normal lower bound')
      call reactive_eb_slip_wall_source_3d(species,state,temp,geometry,wall,ok)
      call require(ok,'wall source')
      volume=geometry%dx*geometry%dy*geometry%dz
      do c=1,nvar
        budget=dt*volume*sum(wall(c,:,:,:)*geometry%volume_fraction) &
          +geometry%dy*geometry%dz*sum(geometry%x_face_fraction(0,:,:)*integral%x(c,0,:,:)- &
            geometry%x_face_fraction(4,:,:)*integral%x(c,4,:,:)) &
          +geometry%dx*geometry%dz*sum(geometry%y_face_fraction(:,0,:)*integral%y(c,:,0,:)- &
            geometry%y_face_fraction(:,4,:)*integral%y(c,:,4,:)) &
          +geometry%dx*geometry%dy*sum(geometry%z_face_fraction(:,:,0)*integral%z(c,:,:,0)- &
            geometry%z_face_fraction(:,:,4)*integral%z(c,:,:,4))
        change=volume*sum((next(c,:,:,:)-state(c,:,:,:))*geometry%volume_fraction)
        call require(abs(change-budget)<=2.e-12_dp*max(1._dp,volume*sum(abs(state(c,:,:,:))* &
          geometry%volume_fraction)),'Euler boundary and wall inventory')
      end do
      saved_flux=integral
      ghost(1,5,4,4)=-1._dp
      call advance_reactive_eb_euler_3d(species,state,temp,geometry,solvers(solver_id),&
        dt,next,nt,ok,ghost,gt,integral,difmag=0._dp)
      call require(.not.ok.and.all(next==state).and.all(nt==temp).and.all(integral%x==saved_flux%x).and. &
        all(integral%y==saved_flux%y).and.all(integral%z==saved_flux%z),'bad donor state and flux rollback')
      call reactive_eb_riemann_fluxes_3d(species,state,temp,geometry,solvers(solver_id),flux%x,flux%y,flux%z, &
        ok,ghost_state=ghost,difmag=0._dp)
      call require(.not.ok.and.all(flux%x==0._dp).and.all(flux%y==0._dp).and.all(flux%z==0._dp), &
        'incomplete halo rejected with zero flux')
      call reactive_eb_riemann_fluxes_3d(species,state,temp,geometry,solvers(solver_id),flux%x,flux%y,flux%z, &
        ok,ghost(:,:4,:,:),gt,difmag=0._dp)
      call require(.not.ok.and.all(flux%x==0._dp).and.all(flux%y==0._dp).and.all(flux%z==0._dp), &
        'wrong halo dimensions rejected')
    end do
  end do
  write(*,'(a)') 'test_eb_reactive_ghost_3d: PASS'
contains
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not.condition) error stop message
  end subroutine require
end program test_eb_reactive_ghost_3d
