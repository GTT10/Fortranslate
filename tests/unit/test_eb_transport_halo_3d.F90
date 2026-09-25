program test_eb_transport_halo_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_species_component, &
    reactive_mass_fraction_component,reactive_primitive_to_conserved,reactive_conserved_to_primitive
  use reactive_transport_3d_mod, only: reactive_transport_face_flux_3d,physical_transport_fluxes_3d,reactive_isothermal_timestep_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d
  use eb_reactive_transport_3d_mod, only: reactive_eb_transport_fluxes_rhs_3d, &
    reactive_eb_transport_euler_update_3d,advance_reactive_eb_transport_3d,reactive_eb_isothermal_timestep_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config,build_configured_eb_geometry_3d
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d, &
    advance_reactive_eb_full_3d,reactive_eb_integrals_3d
  use simulation_config_reactive_3d_mod, only: fill_reactive_boundary_ghosts_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(eb_geometry_3d) :: geometry,extended
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp), allocatable :: inflow(:,:),inflow_t(:),boundary(:),wall(:),chemical(:),before(:),after(:)
  real(dp), allocatable :: reference_boundary(:),reference_wall(:),reference_chemical(:)
  character(len=16) :: labels(6)
  integer :: reacting
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
  call load_h2o2_elementary_mechanism(reactions,ok)
  call require(ok,'reactions')
  nvar=reactive_nvar(size(species));config%nx=4;config%ny=4;config%nz=4
  allocate(u(nvar,4,4,4),t(4,4,4),rhs(nvar,4,4,4),g(nvar,0:5,0:5,0:5),gt(0:5,0:5,0:5))
  allocate(advanced(nvar,4,4,4),advanced_t(4,4,4))
  allocate(stage(nvar,4,4,4),stage_t(4,4,4),twice(nvar,4,4,4),twice_t(4,4,4))
  allocate(reference_boundary(nvar),reference_wall(nvar),reference_chemical(nvar))
  allocate(inflow(nvar,6),inflow_t(6),boundary(nvar),wall(nvar),chemical(nvar),before(nvar),after(nvar))
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
    position=-10._dp
    if(mode==2) position=.31_dp
    call build_axis_plane_eb_geometry_3d(16,16,16,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
      -1.5_dp,2.5_dp,'z',position,extended,ok)
    call require(ok,'coupled CTU extended geometry')
    labels='outflow';labels(1)='inflow'
    do c=1,6
      inflow(:,c)=donor;inflow_t(c)=donor_t
    end do
    call reactive_eb_integrals_3d(u,geometry,before,ok)
    call require(ok,'coupled initial inventory')
    do reacting=0,1
      provider_calls=0;fail_call=0
      call advance_reactive_eb_full_3d(species,reactions,transport,u,t,geometry,'pelec','state_redist', &
        1.e-7_dp,reacting==1,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true., &
        advanced,advanced_t,theta,ok,boundary_change=boundary,wall_change=wall,chemistry_change=chemical, &
        extended_geometry=extended,boundary_faces=labels,inflow_state=inflow,inflow_temperature=inflow_t, &
        prepare_transport_halo=prepare_stage_halo)
      call require(ok.and.provider_calls==4,'coupled CTU rebuilds four transport stage halos')
      call reactive_eb_integrals_3d(advanced,geometry,after,ok)
      call require(ok,'coupled final inventory')
      call require(maxval(abs(after-before-boundary-wall-chemical)/max(1._dp,abs(before)))<=2.e-12_dp, &
        'coupled CTU hydro transport and chemistry boundary inventory')
      if(reacting==1) call require(maxval(abs(chemical))>1.e-15_dp,'coupled chemistry has resolved exchange')
      twice=advanced;twice_t=advanced_t
      reference_boundary=boundary;reference_wall=wall;reference_chemical=chemical
      provider_calls=0
      call advance_reactive_eb_full_3d(species,reactions,transport,u,t,geometry,'pelec','state_redist', &
        1.e-7_dp,reacting==1,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true., &
        advanced,advanced_t,theta,ok,boundary_change=boundary,wall_change=wall,chemistry_change=chemical, &
        extended_geometry=extended,boundary_faces=labels,inflow_state=inflow,inflow_temperature=inflow_t)
      call require(ok.and.provider_calls==0,'built-in coupled inflow donor preparation')
      call require(all(advanced==twice).and.all(advanced_t==twice_t), &
        'built-in donors exactly match independently supplied stage donors')
      call require(all(boundary==reference_boundary).and.all(wall==reference_wall).and. &
        all(chemical==reference_chemical),'built-in donor ledgers match custom provider')
    end do
    provider_calls=0
    call advance_reactive_eb_full_3d(species,reactions,transport,u,t,geometry,'pelec','state_redist', &
      1.e-7_dp,.true.,2.e-7_dp,1.e-12_dp,.false.,.false.,.false.,.false.,.false., &
      advanced,advanced_t,theta,ok,boundary_change=boundary,wall_change=wall,chemistry_change=chemical, &
      extended_geometry=extended,boundary_faces=labels,inflow_state=inflow,inflow_temperature=inflow_t)
    call require(ok.and.provider_calls==0,'transport-disabled route forwards CTU boundaries')
    call reactive_eb_integrals_3d(advanced,geometry,after,ok)
    call require(ok,'Strang final inventory')
    call require(maxval(abs(after-before-boundary-wall-chemical)/max(1._dp,abs(before)))<=2.e-12_dp, &
      'Strang CTU boundary and chemistry inventory')
    provider_calls=0;fail_call=4;boundary=-81;wall=-82;chemical=-83
    call advance_reactive_eb_full_3d(species,reactions,transport,u,t,geometry,'pelec','state_redist', &
      1.e-7_dp,.true.,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true., &
      advanced,advanced_t,theta,ok,boundary_change=boundary,wall_change=wall,chemistry_change=chemical, &
      extended_geometry=extended,boundary_faces=labels,inflow_state=inflow,inflow_temperature=inflow_t, &
      prepare_transport_halo=prepare_stage_halo)
    call require(.not.ok.and.provider_calls==4.and.all(advanced==u).and.all(advanced_t==t), &
      'last coupled transport halo failure restores complete input')
    call require(all(boundary==-81).and.all(wall==-82).and.all(chemical==-83), &
      'last coupled transport halo failure preserves ledgers')
    call advance_reactive_eb_full_3d(species,reactions,transport,u,t,geometry,'pelec','state_redist', &
      1.e-7_dp,.false.,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true., &
      advanced,advanced_t,theta,ok,extended_geometry=extended,boundary_faces=labels)
    call require(.not.ok.and.all(advanced==u).and.all(advanced_t==t),'built-in inflow requires prescribed state and temperature')
    labels='slip_wall';provider_calls=0
    call advance_reactive_eb_full_3d(species,reactions,transport,u,t,geometry,'pelec','state_redist', &
      1.e-7_dp,.false.,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true., &
      advanced,advanced_t,theta,ok,extended_geometry=extended,boundary_faces=labels, &
      prepare_transport_halo=prepare_stage_halo)
    call require(.not.ok.and.provider_calls==0,'unsupported physical wall transport rejects before stepping')
    call verify_wall_transport()
  end do
  write(*,'(a)') 'test_eb_transport_halo_3d: PASS'
contains
  subroutine verify_wall_transport()
    character(len=16), parameter :: kinds(3)=[character(len=16)::'slip_wall','no_slip_wall','isothermal_wall']
    character(len=16) :: wall_faces(6)
    real(dp) :: wall_t(6),weights(3),rho,theta_ref
    type(reactive_eb_3d_config) :: wall_config
    type(eb_geometry_3d) :: wall_extended,thin
    real(dp) :: wall_dt,regular_dt,thin_dt
    real(dp), allocatable :: moving(:,:,:,:),rx(:,:,:,:),ry(:,:,:,:),rz(:,:,:,:)
    integer :: kind,normal,which,point(3),component
    wall_config=config;wall_config%plane_axis='z';wall_config%plane_position=-10._dp
    if(mode==2) wall_config%plane_position=.31_dp
    wall_config%boundary_faces='slip_wall'
    call build_configured_eb_geometry_3d(wall_config,wall_extended,ok,halo=6)
    call require(ok,'configured reflected wall geometry')
    allocate(moving,source=u)
    allocate(rx(nvar,0:4,4,4),ry(nvar,4,0:4,4),rz(nvar,4,4,0:4))
    do k=1,4
    do j=1,4
    do i=1,4
      rho=moving(1,i,j,k)
      moving(5,i,j,k)=moving(5,i,j,k)-sum(moving(2:4,i,j,k)**2)/(2*rho)+25*rho
      moving(2:4,i,j,k)=rho*[3._dp,4._dp,5._dp]
    end do
    end do
    end do
    wall_t=temp+100
    weights=[geometry%dx,geometry%dy,geometry%dz]
    do kind=1,3
    do normal=1,3
    do which=0,1
      if(mode==2.and.normal==3.and.which==0) cycle
      wall_faces='outflow';wall_faces(2*normal-1+which)=kinds(kind)
      call fill_reactive_boundary_ghosts_3d(moving,t,wall_faces,g,gt,ok,viscous_wall=.true.,ghost_width=1)
      call require(ok,'wall donor fill')
      faces=.true.
      call reactive_eb_transport_fluxes_rhs_3d(species,transport,moving,t,geometry,0._dp, &
        .true.,.true.,.true.,.true.,rhs,x,y,z,theta,ok,g,gt,faces,wall_faces,wall_t)
      call require(ok,'EB wall transport flux')
      call physical_transport_fluxes_3d(species,transport,moving,t,[4,4,4],weights,0._dp,wall_faces, &
        .true.,.true.,.true.,.true.,rx,ry,rz,theta_ref,ok,wall_temperature=wall_t)
      call require(ok,'regular wall transport reference')
      point=3;point(normal)=4*which
      select case(normal)
      case(1);expected=x(:,point(1),point(2),point(3));uniform=rx(:,point(1),point(2),point(3))
      case(2);expected=y(:,point(1),point(2),point(3));uniform=ry(:,point(1),point(2),point(3))
      case(3);expected=z(:,point(1),point(2),point(3));uniform=rz(:,point(1),point(2),point(3))
      end select
      call require(maxval(abs(expected-uniform))<=3.e-12_dp*max(1._dp,maxval(abs(uniform))), &
        'wall flux agrees with regular transport reference for uniform moving state')
      call require(expected(1)==0.and.all(expected(6:)==0),'impermeable wall species and mass flux')
      if(kind<3) call require(expected(5)==0,'adiabatic stationary wall does no work')
      if(kind==3) call require(real(1-2*which,dp)*expected(5)>0,'hot wall heat enters fluid')
    end do
    end do
    wall_faces=kinds(kind)
    call advance_reactive_eb_full_3d(species,reactions,transport,moving,t,geometry,'pelec','state_redist', &
      1.e-7_dp,.true.,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true., &
      advanced,advanced_t,theta,ok,boundary_change=boundary,wall_change=wall,chemistry_change=chemical, &
      extended_geometry=wall_extended,boundary_faces=wall_faces,wall_temperature=wall_t)
    call require(ok,'coupled wall reaction transport CTU step')
    call reactive_eb_integrals_3d(moving,geometry,before,ok)
    call require(ok,'wall initial inventory')
    call reactive_eb_integrals_3d(advanced,geometry,after,ok)
    call require(ok,'wall final inventory')
    call require(maxval(abs(after-before-boundary-wall-chemical)/max(1._dp,abs(before)))<=2.e-12_dp, &
      'coupled wall inventory closes')
    do component=6,nvar
      call require(abs(boundary(component))<=1.e-25_dp,'no wall species exchange')
    end do
    end do
    wall_faces='isothermal_wall'
    call reactive_eb_isothermal_timestep_3d(species,transport,moving,t,geometry,wall_faces,wall_t,.35_dp,wall_dt,ok)
    call require(ok.and.wall_dt>0,'EB wall thermal bound')
    call reactive_isothermal_timestep_3d(species,transport,moving,t,weights,wall_faces,wall_t,.35_dp,regular_dt,ok)
    call require(ok.and.abs(wall_dt-regular_dt)<=3.e-12_dp*regular_dt, &
      'full-fluid and tangential-cut wall bounds match regular corner rate')
    call build_axis_plane_eb_geometry_3d(4,4,4,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',.95_dp,thin,ok)
    call require(ok,'thin fluid wall geometry')
    wall_faces='outflow';wall_faces(2)='isothermal_wall'
    call reactive_isothermal_timestep_3d(species,transport,moving,t,weights,wall_faces,wall_t,.35_dp,regular_dt,ok)
    call require(ok,'regular single-wall rate')
    wall_faces(1)='isothermal_wall';wall_t(1)=-1
    call reactive_eb_isothermal_timestep_3d(species,transport,moving,t,thin,wall_faces,wall_t,.35_dp,thin_dt,ok)
    call require(ok.and.abs(thin_dt-.2_dp*regular_dt)<=3.e-12_dp*regular_dt, &
      'thin-cell aperture over fluid volume strengthens wall bound; covered wall ignored')
    wall_t(2)=-1
    call reactive_eb_isothermal_timestep_3d(species,transport,moving,t,thin,wall_faces,wall_t,.35_dp,thin_dt,ok)
    call require(.not.ok.and.thin_dt==0,'invalid active wall resets thermal timestep')
    wall_t=-1;wall_faces='isothermal_wall';boundary=-81
    call advance_reactive_eb_full_3d(species,reactions,transport,moving,t,geometry,'pelec','state_redist', &
      1.e-7_dp,.false.,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true., &
      advanced,advanced_t,theta,ok,boundary_change=boundary,extended_geometry=wall_extended, &
      boundary_faces=wall_faces,wall_temperature=wall_t)
    call require(.not.ok.and.all(advanced==moving).and.all(advanced_t==t).and.all(boundary==-81), &
      'invalid wall temperature rolls back state and boundary ledger')
  end subroutine verify_wall_transport
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
