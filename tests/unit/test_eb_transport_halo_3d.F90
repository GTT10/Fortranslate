program test_eb_transport_halo_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_species_component, &
    reactive_mass_fraction_component,reactive_primitive_to_conserved,reactive_conserved_to_primitive
  use reactive_transport_3d_mod, only: reactive_transport_face_flux_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d
  use eb_reactive_transport_3d_mod, only: reactive_eb_transport_fluxes_rhs_3d, &
    reactive_eb_transport_euler_update_3d,advance_reactive_eb_transport_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use simulation_config_reactive_3d_mod, only: fill_reactive_boundary_ghosts_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(eb_geometry_3d) :: geometry
  type(reactive_eb_3d_config) :: config
  type(eb_integrated_flux_3d) :: integrated,saved_flux,first_flux,second_flux
  real(dp), allocatable :: stage(:,:,:,:),stage_t(:,:,:),twice(:,:,:,:),twice_t(:,:,:)
  integer :: provider_calls,fail_call
  logical :: observed_stage_change
  real(dp), allocatable :: advanced(:,:,:,:),advanced_t(:,:,:)
  real(dp), parameter :: dt=1.e-3_dp
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:),g(:,:,:,:),gt(:,:,:),rhs(:,:,:,:),x(:,:,:,:),y(:,:,:,:),z(:,:,:,:)
  real(dp), allocatable :: p(:),donor(:),gp(:),expected(:),left(:),right(:),uniform(:)
  real(dp) :: temp,speed,donor_t,checked_t,gradient(3,3),species_energy,theta,volume,budget,change,position,guess_t
  integer :: nvar,i,j,k,mode,axis,side,idx(3),c,s
  logical :: ok,faces(6)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  nvar=reactive_nvar(size(species));config%nx=4;config%ny=4;config%nz=4
  allocate(u(nvar,4,4,4),t(4,4,4),rhs(nvar,4,4,4),g(nvar,0:5,0:5,0:5),gt(0:5,0:5,0:5))
  allocate(advanced(nvar,4,4,4),advanced_t(4,4,4))
  allocate(stage(nvar,4,4,4),stage_t(4,4,4),twice(nvar,4,4,4),twice_t(4,4,4))
  allocate(x(nvar,0:4,4,4),y(nvar,4,0:4,4),z(nvar,4,4,0:4))
  allocate(p(reactive_nprim(size(species))))
  allocate(gp(size(p)),left(size(p)),right(size(p)), &
    donor(nvar),expected(nvar),uniform(nvar))
  do mode=1,2
    position=-.1_dp
    if(mode==2) position=.31_dp
    call build_axis_plane_eb_geometry_3d(4,4,4,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'z',position,geometry,ok)
    call require(ok,'geometry')
    call initialize_reactive_eb_density_sheet_3d(species,config,geometry,u,t,ok)
    call require(ok,'initial')
    uniform=u(:,4,4,4);temp=t(4,4,4)
    do k=1,4
      do j=1,4
        do i=1,4
          u(:,i,j,k)=uniform;t(i,j,k)=temp
        end do
      end do
    end do
    call reactive_conserved_to_primitive(species,uniform,temp,p,checked_t,speed,ok)
    call require(ok,'uniform primitive')
    gp=p;gp(1)=gp(1)*1.1_dp;gp(5)=gp(5)*1.2_dp;gp(2:4)=[3._dp,4._dp,5._dp]
    gp(reactive_mass_fraction_component(1))=gp(reactive_mass_fraction_component(1))+.01_dp
    gp(reactive_mass_fraction_component(size(species)))=gp(reactive_mass_fraction_component(size(species)))-.01_dp
    call reactive_primitive_to_conserved(species,gp,donor,donor_t,speed,ok)
    call require(ok,'nonuniform donor')
    guess_t=donor_t
    call reactive_conserved_to_primitive(species,donor,guess_t,gp,donor_t,speed,ok)
    call require(ok,'checked donor primitive')
    do k=0,5
      do j=0,5
        do i=0,5
          g(:,i,j,k)=donor;gt(i,j,k)=donor_t
        end do
      end do
    end do
    faces=.false.
    call reactive_eb_transport_fluxes_rhs_3d(species,transport,u,t,geometry,0._dp, &
      .true.,.true.,.true.,.true.,rhs,x,y,z,theta,ok,g,gt,faces)
    call require(ok.and.maxval(abs(rhs))==0._dp,'disabled boundary faces preserve uniform legacy result')
    do axis=1,3
      do side=0,1
        if(mode==2.and.axis==3.and.side==0) cycle
        faces=.false.;faces(2*axis-1+side)=.true.
        call reactive_eb_transport_fluxes_rhs_3d(species,transport,u,t,geometry,0._dp, &
          .true.,.true.,.true.,.true.,rhs,x,y,z,theta,ok,g,gt,faces)
        call require(ok.and.theta==1._dp,'halo transport flux')
        gradient=0._dp;left=p;right=gp
        if(side==0) then
          left=gp;right=p
        end if
        gradient(:,axis)=(right(2:4)-left(2:4))/.25_dp
        if(side==0) then
          call reactive_transport_face_flux_3d(species,transport,left,right,donor_t,checked_t,.25_dp,axis, &
            gradient,.true.,.true.,.true.,.true.,expected,species_energy,ok)
        else
          call reactive_transport_face_flux_3d(species,transport,left,right,checked_t,donor_t,.25_dp,axis, &
            gradient,.true.,.true.,.true.,.true.,expected,species_energy,ok)
        end if
        call require(ok,'independent face transport')
        idx=2;idx(axis)=4*side
        select case(axis)
        case(1);uniform=x(:,idx(1),idx(2),idx(3))
        case(2);uniform=y(:,idx(1),idx(2),idx(3))
        case(3);uniform=z(:,idx(1),idx(2),idx(3))
        end select
        call require(maxval(abs(uniform-expected))<=2.e-12_dp*max(1._dp,maxval(abs(expected))), &
          'boundary transport matches directional face kernel')
        call require(maxval(abs(uniform))>0._dp,'nonzero boundary transport')
        volume=geometry%dx*geometry%dy*geometry%dz
        do c=1,nvar
          budget=.25_dp**2*(sum(geometry%x_face_fraction(0,:,:)*x(c,0,:,:)- &
            geometry%x_face_fraction(4,:,:)*x(c,4,:,:))+sum(geometry%y_face_fraction(:,0,:)*y(c,:,0,:)- &
            geometry%y_face_fraction(:,4,:)*y(c,:,4,:))+sum(geometry%z_face_fraction(:,:,0)*z(c,:,:,0)- &
            geometry%z_face_fraction(:,:,4)*z(c,:,:,4)))
          change=volume*sum(rhs(c,:,:,:)*geometry%volume_fraction)
          call require(abs(change-budget)<=2.e-12_dp*max(1._dp,abs(budget)),'transport aperture inventory')
        end do
        call reactive_eb_transport_euler_update_3d(species,transport,u,t,geometry,dt, &
          .true.,.true.,.true.,.true.,.5_dp,advanced,advanced_t,theta,ok,g,gt,faces,integrated)
        call require(ok,'redistributed halo transport update')
        call require(maxval(abs(advanced-u))>0._dp,'halo exchange changes state')
        call require(maxval(abs(integrated%x-dt*x))<=2.e-12_dp*max(1._dp,maxval(abs(dt*x))).and. &
          maxval(abs(integrated%y-dt*y))<=2.e-12_dp*max(1._dp,maxval(abs(dt*y))).and. &
          maxval(abs(integrated%z-dt*z))<=2.e-12_dp*max(1._dp,maxval(abs(dt*z))), &
          'redistributed update returns actual time-integrated face fluxes')
        do c=1,nvar
          budget=.25_dp**2*(sum(geometry%x_face_fraction(0,:,:)*integrated%x(c,0,:,:)- &
            geometry%x_face_fraction(4,:,:)*integrated%x(c,4,:,:))+ &
            sum(geometry%y_face_fraction(:,0,:)*integrated%y(c,:,0,:)- &
            geometry%y_face_fraction(:,4,:)*integrated%y(c,:,4,:))+ &
            sum(geometry%z_face_fraction(:,:,0)*integrated%z(c,:,:,0)- &
            geometry%z_face_fraction(:,:,4)*integrated%z(c,:,:,4)))
          change=volume*sum((advanced(c,:,:,:)-u(c,:,:,:))*geometry%volume_fraction)
          call require(abs(change-budget)<=2.e-12_dp*max(1._dp, &
            volume*sum(abs(u(c,:,:,:))*geometry%volume_fraction)), 'redistributed halo inventory')
          call require(all(advanced(c,:,:,:)==u(c,:,:,:).or.geometry%volume_fraction>0), &
            'covered storage preserved')
        end do
        saved_flux=integrated
        call reactive_eb_transport_euler_update_3d(species,transport,u,t,geometry,dt, &
          .true.,.true.,.true.,.true.,-1._dp,advanced,advanced_t,theta,ok,g,gt,faces,integrated)
        call require(.not.ok.and.all(advanced==u).and.all(advanced_t==t), &
          'redistribution failure restores input state and temperature')
        call require(all(integrated%x==saved_flux%x).and.all(integrated%y==saved_flux%y).and. &
          all(integrated%z==saved_flux%z),'failed redistribution preserves published flux')
        call reactive_eb_transport_fluxes_rhs_3d(species,transport,u,t,geometry,1.e8_dp, &
          .true.,.true.,.true.,.true.,rhs,x,y,z,theta,ok,g,gt,faces)
        call require(ok.and.theta<1._dp,'boundary species limiter active')
        do s=1,size(species)
          c=reactive_species_component(s)
          call require(minval(u(c,:,:,:)+1.e8_dp*rhs(c,:,:,:),mask=geometry%volume_fraction>0)>=-1.e-12_dp, &
            'limited boundary species inventory nonnegative')
        end do
      end do
    end do
    call reactive_eb_transport_fluxes_rhs_3d(species,transport,u,t,geometry,0._dp, &
      .true.,.true.,.true.,.true.,rhs,x,y,z,theta,ok,ghost_state=g)
    call require(.not.ok,'incomplete transport halo rejected')
    call reactive_eb_transport_euler_update_3d(species,transport,u,t,geometry,dt, &
      .true.,.true.,.true.,.true.,.5_dp,advanced,advanced_t,theta,ok, &
      ghost_state=g,integrated_flux=integrated)
    call require(.not.ok.and.all(advanced==u).and.all(advanced_t==t),'incomplete halo update rejected')
    call require(all(integrated%x==saved_flux%x).and.all(integrated%y==saved_flux%y).and. &
      all(integrated%z==saved_flux%z),'failed halo preserves published flux')
    call reactive_eb_transport_euler_update_3d(species,transport,u,t,geometry,0._dp, &
      .true.,.true.,.true.,.true.,.5_dp,advanced,advanced_t,theta,ok,integrated_flux=integrated)
    call require(ok.and.all(advanced==u).and.all(advanced_t==t),'zero interval identity')
    call require(all(integrated%x==0).and.all(integrated%y==0).and.all(integrated%z==0), &
      'zero interval clears exchange')
    provider_calls=0;fail_call=0;observed_stage_change=.false.
    advanced=u;advanced_t=t
    call advance_reactive_eb_transport_3d(species,transport,advanced,advanced_t,geometry,dt, &
      .true.,.true.,.true.,.true.,.5_dp,theta,ok,prepare_stage_halo,integrated)
    call require(ok.and.provider_calls==2.and.observed_stage_change,'both stages rebuild current donors')
    saved_flux=integrated
    call prepare_stage_halo(u,t,geometry,g,gt,faces,ok)
    call require(ok,'manual initial halo')
    call reactive_eb_transport_euler_update_3d(species,transport,u,t,geometry,dt, &
      .true.,.true.,.true.,.true.,.5_dp,stage,stage_t,theta,ok,g,gt,faces,first_flux)
    call require(ok,'manual first stage')
    call prepare_stage_halo(stage,stage_t,geometry,g,gt,faces,ok)
    call require(ok,'manual second halo')
    call reactive_eb_transport_euler_update_3d(species,transport,stage,stage_t,geometry,dt, &
      .true.,.true.,.true.,.true.,.5_dp,twice,twice_t,theta,ok,g,gt,faces,second_flux)
    call require(ok,'manual second stage')
    call require(maxval(abs(advanced-.5_dp*(u+twice)))<2.e-12_dp*max(1._dp,maxval(abs(u))), &
      'SSPRK2 agrees with independently assembled stages')
    call require(all(integrated%x==.5_dp*first_flux%x+.5_dp*second_flux%x).and. &
      all(integrated%y==.5_dp*first_flux%y+.5_dp*second_flux%y).and. &
      all(integrated%z==.5_dp*first_flux%z+.5_dp*second_flux%z),'SSPRK2 integrated flux weights')
    do c=1,nvar
      budget=.25_dp**2*(sum(geometry%x_face_fraction(0,:,:)*integrated%x(c,0,:,:)- &
        geometry%x_face_fraction(4,:,:)*integrated%x(c,4,:,:))+ &
        sum(geometry%y_face_fraction(:,0,:)*integrated%y(c,:,0,:)- &
        geometry%y_face_fraction(:,4,:)*integrated%y(c,:,4,:))+ &
        sum(geometry%z_face_fraction(:,:,0)*integrated%z(c,:,:,0)- &
        geometry%z_face_fraction(:,:,4)*integrated%z(c,:,:,4)))
      change=volume*sum((advanced(c,:,:,:)-u(c,:,:,:))*geometry%volume_fraction)
      call require(abs(change-budget)<=2.e-12_dp*max(1._dp, &
        volume*sum(abs(u(c,:,:,:))*geometry%volume_fraction)), 'SSPRK2 halo inventory')
    end do
    provider_calls=0;fail_call=2;advanced=u;advanced_t=t
    call advance_reactive_eb_transport_3d(species,transport,advanced,advanced_t,geometry,dt, &
      .true.,.true.,.true.,.true.,.5_dp,theta,ok,prepare_stage_halo,integrated)
    call require(.not.ok.and.provider_calls==2.and.all(advanced==u).and.all(advanced_t==t), &
      'second-stage halo failure rolls back whole interval')
    call require(all(integrated%x==saved_flux%x).and.all(integrated%y==saved_flux%y).and. &
      all(integrated%z==saved_flux%z),'second-stage halo failure preserves flux')
  end do
  write(*,'(a)') 'test_eb_transport_halo_3d: PASS'
contains
  subroutine prepare_stage_halo(state,temperature,metrics,ghost,ghost_temperature,enabled,valid)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:)
    type(eb_geometry_3d), intent(in) :: metrics
    real(dp), intent(out) :: ghost(:,0:,0:,0:),ghost_temperature(0:,0:,0:)
    logical, intent(out) :: enabled(6),valid
    character(len=16) :: labels(6)
    real(dp) :: inlet(size(state,1),6),inlet_t(6)
    integer :: face
    provider_calls=provider_calls+1;valid=.false.;enabled=.true.
    if(provider_calls==fail_call) return
    if(any(shape(temperature)/=[metrics%nx,metrics%ny,metrics%nz])) return
    if(provider_calls==2) observed_stage_change=any(state/=u)
    labels='outflow';labels(1)='inflow'
    do face=1,6
      inlet(:,face)=donor;inlet_t(face)=donor_t
    end do
    call fill_reactive_boundary_ghosts_3d(state,temperature,labels,ghost,ghost_temperature,valid, &
      viscous_wall=.true.,inflow_state=inlet,inflow_temperature=inlet_t,ghost_width=1)
  end subroutine prepare_stage_halo
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if(.not.condition) error stop message
  end subroutine require
end program test_eb_transport_halo_3d
