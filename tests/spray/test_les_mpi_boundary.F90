program test_les_mpi_boundary
  use mpi_f08
  use spray_mpi_mod, only: fill_mpi_spray_ghosts,deposit_mpi_spray_field,migrate_spray_parcels,migrate_mpi_wall_films
  use rotating_frame_3d_mod, only: rotating_frame_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved,reactive_conserved_to_primitive
  use les_sdc_3d_mod, only: advance_les_sdc_3d
  use reactive_mol_3d_mod, only: advance_reactive_mol_3d,convert_regular_energy_3d,reactive_mol_blend_agreement
  use reactive_3d_mod, only: advance_reactive_sdc_3d
  use constant_volume_reactor_mod, only: forced_reactor_integrator,advance_forced_reactor_interval
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use spray_mpi_sdc_mod, only: advance_mpi_cartesian_gas_sdc,sample_mpi_sdc_primitive,kick_mpi_sdc_transfer, &
    deposit_mpi_sdc_cell_source,kick_mpi_sdc_films,mpi_sdc_film_stage, &
    advance_mpi_sdc_film_stage,append_mpi_sdc_film_births,kick_mpi_sdc_breakup,drift_mpi_sdc_particles, &
    mpi_sdc_particle_stage,advance_mpi_sdc_particle_stage,advance_mpi_cartesian_spray_sdc
  use spray_reactive_3d_mod, only: advance_spray_reactive_sdc_3d
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options,wall_film_totals
  use spray_sdc_3d_mod, only: kick_serial_sdc_films,spray_sdc_stage_3d,advance_spray_sdc_stage_3d
  use spray_parcel_mod, only: spray_parcel,spray_liquid,spray_options,initialize_spray_liquid,parcel_totals
  use spray_mixture_mod, only: advance_mixture_frozen_gas
  use spray_coupling_3d_mod, only: particle_boundary_stencil
  use mixture_thermo_mod, only: mixture_pressure
  use, intrinsic :: iso_fortran_env, only: int64
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_elementary_transport
  use spray_partition_mod, only: spray_cartesian_block,spray_cartesian_owner
  use les_mpi_dynamic_mod, only: compute_mpi_dynamic_les_boundary_rhs_3d,advance_mpi_dynamic_les_ssprk2_3d
  use les_dynamic_smagorinsky_mod, only: compute_dynamic_les_boundary_rhs_3d,advance_dynamic_les_ssprk2_3d
  use les_mpi_boundary_mod, only: fill_mpi_les_boundary_ghosts_3d,fill_mpi_sdc_source_ghosts_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use simulation_config_reactive_3d_mod, only: fill_reactive_boundary_ghosts_3d,fill_pressure_outlet_ghosts_3d, &
    fill_mach_inlet_boundary_ghosts_3d
  implicit none
  integer, parameter :: cells(3)=[5,3,2],width=7,nv=12
  type(nasa7_species), allocatable :: species(:)
  real(dp) :: full(nv,5,3,2),temperature(5,3,2),q(nv),sound,walls(3,6),pressure(6),inlet(nv,6),inlet_t(6)
  real(dp) :: mach(6),pin(6)
  real(dp) :: reference(nv,19,17,16),reference_t(19,17,16)
  real(dp), allocatable :: local(:,:,:,:),local_t(:,:,:),ghost(:,:,:,:),gt(:,:,:),saved(:,:,:,:),saved_t(:,:,:)
  real(dp) :: rhs_full(nv,5,3,2),serial_rate(nv),rate(nv),diffusion,serial_diffusion,h(3),wall_t(6)
  real(dp) :: serial_u(nv,5,3,2),serial_t(5,3,2),serial_budget(nv),budget(nv),saved_budget(nv),scale
  real(dp), allocatable :: rhs(:,:,:,:),u(:,:,:,:),t(:,:,:)
  character(len=32) :: faces(6)
  integer :: rank,nranks,ierr,first(3),extent(3),i,j,k,mode,lo,hi,requested_width,iv
  logical :: ok
  character(len=32) :: test_mode
  integer :: agreement_calls=0,source_calls=0,reactor_calls=0
  logical :: fail_les_halo=.false.,fail_les_auxiliary=.false.
  integer :: profile_unit,profile_calls=0,profile_phase=0,profile_io
  character(len=1024) :: profile_prefix,profile_path
  character(len=32) :: profile_mode,profile_rank,profile_index
  call get_command_argument(1,test_mode)
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call spray_cartesian_block(cells,[nranks,1,1],rank,first,extent,ok)
  call require(ok,'partition')
  if(test_mode=='source-halo') then
    call check_source_halo()
    call MPI_Finalize(ierr)
    stop
  endif
  do k=1,2
    do j=1,3
      do i=1,5
        q=0;q(1)=.3_dp;q(2:4)=[real(i,dp),real(j,dp),real(k,dp)];q(5)=101325+100*i;q(12)=1
        call reactive_primitive_to_conserved(species,q,full(:,i,j,k),temperature(i,j,k),sound,ok)
        call require(ok,'initial state')
      end do
    end do
  end do
  lo=first(1)+1;hi=first(1)+extent(1)
  if(test_mode=='sdc-agreement') then
    call check_sdc_agreement()
    call MPI_Finalize(ierr)
    stop
  endif
  local=full(:,lo:hi,:,:);local_t=temperature(lo:hi,:,:)
  if(test_mode=='mol-blend-agreement') then
    call check_mol_blend_agreement()
    call MPI_Finalize(ierr)
    stop
  endif
  if(test_mode=='sdc-coupled'.or.test_mode=='mol-coupled'.or.test_mode=='mol-coupled-rotation') then
    call get_command_argument(3,profile_mode)
    if(profile_mode=='profile') then
      call get_command_argument(4,profile_prefix)
      call require(len_trim(profile_prefix)>0,'reactor profile prefix')
      write(profile_rank,'(i0)') rank
      profile_path=trim(profile_prefix)//'.rank'//trim(profile_rank)//'.csv'
      open(newunit=profile_unit,file=trim(profile_path),status='new',action='write',iostat=profile_io)
      call require(profile_io==0,'reactor profile log')
      write(profile_unit,'(a)') 'call,phase,event,seconds,steps,ok,minimum_rhoY,temperature'
      flush(profile_unit)
      call check_coupled_sdc(observe_native_reactor)
      close(profile_unit)
    else
      call check_coupled_sdc()
    endif
    call MPI_Finalize(ierr)
    stop
  endif
  if(test_mode=='sdc-particle-stage') then
    call check_sdc_particle_stage()
    call MPI_Finalize(ierr)
    stop
  endif
  if(test_mode=='sdc-drift') then
    call check_sdc_drift()
    call MPI_Finalize(ierr)
    stop
  endif
  if(test_mode=='sdc-breakup') then
    call check_sdc_breakup()
    call MPI_Finalize(ierr)
    stop
  endif
  if(test_mode=='sdc-film') then
    call check_sdc_film()
    call MPI_Finalize(ierr)
    stop
  endif
  if(test_mode=='sdc-transfer') then
    call check_sdc_transfer()
    call MPI_Finalize(ierr)
    stop
  endif
  if(test_mode=='gas-sdc-les') then
    call check_cartesian_les_sdc()
    call MPI_Finalize(ierr)
    stop
  endif
  if(test_mode=='gas-sdc'.or.test_mode=='gas-mol'.or.test_mode=='gas-mol-rotation') then
    call check_cartesian_gas_sdc()
    call MPI_Finalize(ierr)
    stop
  endif
  inlet=spread(full(:,1,1,1),2,6);inlet_t=temperature(1,1,1)
  allocate(rhs,mold=local)
  h=[.001_dp,.0015_dp,.002_dp];wall_t=1100
  do mode=1,8
    walls=0;pressure=0;faces='periodic';mach=-1;pin=0
    select case(mode)
    case(2)
      faces='outflow';faces(1)='inflow'
    case(3)
      faces='slip_wall'
    case(4)
      faces='no_slip_wall';walls(2,1)=1;walls(2,2)=3;walls(1,3)=2
    case(5)
      faces='isothermal_wall';walls(2,1)=1;walls(2,2)=3;walls(1,3)=2
    case(6)
      faces='outflow';faces(3:4)='no_slip_wall';walls(1,3)=2;pressure(2)=101000
    case(7)
      faces='outflow';faces(1)='inflow';faces(3:4)='no_slip_wall';walls(1,3)=2;mach(1)=.2_dp
    case(8)
      faces='outflow';faces(1)='inflow';faces(6)='inflow';mach(6)=.3_dp;pin(6)=120000
    end select
    call fill_mach_inlet_boundary_ghosts_3d(species,full,temperature,faces,mach,reference,reference_t,ok,viscous_wall=.true., &
      inflow_state=inlet,inflow_temperature=inlet_t,wall_velocity=walls,ghost_width=width,inflow_pressure=pin)
    call require(ok,'serial ordinary donors')
    call fill_pressure_outlet_ghosts_3d(species,cells,faces,pressure,reference,reference_t,ok,ghost_width=width)
    call require(ok,'serial pressure donors')
    call fill_mpi_les_boundary_ghosts_3d(MPI_COMM_WORLD,species,local,local_t,cells,faces,width,ghost,gt,ok, &
      inlet,inlet_t,walls,pressure,mach,pin)
    call require(ok,'distributed boundary donors')
    call require(all(ghost==reference(:,lo:hi+2*width,:,:)),'exact full-width serial states')
    call require(all(gt==reference_t(lo:hi+2*width,:,:)),'exact full-width serial temperatures')
    if(mode>=7) then
      do iv=0,1
        call fill_mach_inlet_boundary_ghosts_3d(species,full,temperature,faces,mach, &
          reference(:,1:9,1:7,1:6),reference_t(1:9,1:7,1:6),ok,viscous_wall=iv==1, &
          inflow_state=inlet,inflow_temperature=inlet_t,inflow_pressure=pin,wall_velocity=walls)
        call require(ok,'serial two-layer Mach donors')
        call fill_mpi_spray_ghosts(MPI_COMM_WORLD,local,local_t,cells,faces,ghost,gt,ok,viscous_wall=iv==1, &
          inflow_state=inlet,inflow_temperature=inlet_t,wall_velocity=walls,species=species, &
          inflow_mach=mach,inflow_pressure=pin)
        call require(ok,'MPI spray Mach donors')
        call require(all(ghost==reference(:,lo:hi+4,1:7,1:6)),'exact inviscid/viscous spray Mach states')
        call require(all(gt==reference_t(lo:hi+4,1:7,1:6)),'exact inviscid/viscous spray Mach temperatures')
      end do
    end if
    call compute_dynamic_les_boundary_rhs_3d(species,full,temperature,h,'box',2,rhs_full,serial_diffusion,ok, &
      boundary_faces=faces,boundary_rate=serial_rate,inflow_state=inlet,inflow_temperature=inlet_t, &
      wall_temperature=wall_t,wall_velocity=walls,outlet_pressure=pressure,inflow_mach=mach,inflow_pressure=pin)
    call require(ok,'serial physical RHS')
    call compute_mpi_dynamic_les_boundary_rhs_3d(MPI_COMM_WORLD,species,local,local_t,h,cells,'box',2, &
      rhs,diffusion,ok,boundary_faces=faces,boundary_rate=rate,inflow_state=inlet,inflow_temperature=inlet_t, &
      wall_temperature=wall_t,wall_velocity=walls,outlet_pressure=pressure,inflow_mach=mach,inflow_pressure=pin)
    call require(ok,'distributed physical RHS')
    if(.not.all(rhs==rhs_full(:,lo:hi,:,:))) print *, 'RHS mismatch mode/rank/gap:',mode,rank, &
      maxval(abs(rhs-rhs_full(:,lo:hi,:,:)))
    call require(all(rhs==rhs_full(:,lo:hi,:,:)),'exact serial RHS')
    scale=max(maxval(abs(serial_rate)),1e-30_dp)
    call require(maxval(abs(rate-serial_rate))<=1e-12_dp*scale,'global boundary rate')
    if(mode==2.or.mode>=5) then
      serial_u=full;serial_t=temperature;u=local;t=local_t
      call advance_dynamic_les_ssprk2_3d(species,serial_u,serial_t,h,'box',2,1e-9_dp,ok, &
        boundary_faces=faces,boundary_change=serial_budget,inflow_state=inlet,inflow_temperature=inlet_t, &
        wall_temperature=wall_t,wall_velocity=walls,outlet_pressure=pressure,inflow_mach=mach,inflow_pressure=pin)
      call require(ok,'serial physical step')
      call advance_mpi_dynamic_les_ssprk2_3d(MPI_COMM_WORLD,species,u,t,h,cells,'box',2,1e-9_dp,ok, &
        boundary_faces=faces,boundary_change=budget,inflow_state=inlet,inflow_temperature=inlet_t, &
        wall_temperature=wall_t,wall_velocity=walls,outlet_pressure=pressure,inflow_mach=mach,inflow_pressure=pin)
      call require(ok,'distributed physical step')
      call require(all(u==serial_u(:,lo:hi,:,:)).and.all(t==serial_t(lo:hi,:,:)),'exact serial SSPRK2')
      scale=max(maxval(abs(serial_budget)),1e-30_dp)
      call require(maxval(abs(budget-serial_budget))<=1e-12_dp*scale,'integrated boundary budget')
      saved_budget=budget
      call advance_mpi_dynamic_les_ssprk2_3d(MPI_COMM_WORLD,species,u,t,h,cells,'box',2,1._dp,ok, &
        boundary_faces=faces,boundary_change=budget,inflow_state=inlet,inflow_temperature=inlet_t, &
        wall_temperature=wall_t,wall_velocity=walls,outlet_pressure=pressure,inflow_mach=mach,inflow_pressure=pin)
      call require(.not.ok,'excessive dt rejects')
      call require(all(u==serial_u(:,lo:hi,:,:)).and.all(t==serial_t(lo:hi,:,:)).and. &
        all(budget==saved_budget),'step rollback')
    end if
  end do
  saved=ghost;saved_t=gt
  if(rank==nranks-1) mach(6)=.4_dp
  call fill_mpi_les_boundary_ghosts_3d(MPI_COMM_WORLD,species,local,local_t,cells,faces,width,ghost,gt,ok, &
    inlet,inlet_t,walls,pressure,mach,pin)
  call require(.not.ok.and.all(ghost==saved).and.all(gt==saved_t),'Mach disagreement rollback')
  mach(6)=.3_dp
  if(rank==nranks-1) pin(6)=121000
  rhs=123;diffusion=456;rate=789
  call compute_mpi_dynamic_les_boundary_rhs_3d(MPI_COMM_WORLD,species,local,local_t,h,cells,'box',2, &
    rhs,diffusion,ok,boundary_faces=faces,boundary_rate=rate,inflow_state=inlet,inflow_temperature=inlet_t, &
    inflow_mach=mach,inflow_pressure=pin)
  call require(.not.ok.and.all(rhs==123).and.diffusion==456.and.all(rate==789),'Mach pressure disagreement rollback')
  pin=0;mach=-1;mach(1)=.2_dp;faces='outflow';faces(1)='inflow'
  call fill_mpi_les_boundary_ghosts_3d(MPI_COMM_WORLD,species,local,local_t,cells,faces,width,ghost,gt,ok, &
    inflow_temperature=inlet_t,inflow_mach=mach)
  call require(ok,'Mach-only inlet needs no reservoir table')
  saved=ghost;saved_t=gt
  if(rank==nranks-1) inlet_t(1)=-1
  call fill_mpi_les_boundary_ghosts_3d(MPI_COMM_WORLD,species,local,local_t,cells,faces,width,ghost,gt,ok, &
    inflow_temperature=inlet_t,inflow_mach=mach)
  call require(.not.ok.and.all(ghost==saved).and.all(gt==saved_t),'Mach temperature disagreement rollback')
  inlet_t=temperature(1,1,1)
  faces='outflow';faces(3:4)='no_slip_wall';walls=0;walls(1,3)=2;pressure=0;pressure(2)=101000
  requested_width=width
  if(rank==nranks-1) requested_width=width-1
  call fill_mpi_les_boundary_ghosts_3d(MPI_COMM_WORLD,species,local,local_t,cells,faces,requested_width,ghost,gt,ok, &
    inlet,inlet_t,walls,pressure)
  call require(.not.ok.and.all(ghost==saved).and.all(gt==saved_t),'width disagreement rollback')
  if(rank==nranks-1) walls(1,3)=3
  call fill_mpi_les_boundary_ghosts_3d(MPI_COMM_WORLD,species,local,local_t,cells,faces,width,ghost,gt,ok, &
    inlet,inlet_t,walls,pressure)
  call require(.not.ok.and.all(ghost==saved).and.all(gt==saved_t),'wall disagreement rollback')
  walls(1,3)=2
  if(rank==nranks-1) pressure(2)=102000
  call fill_mpi_les_boundary_ghosts_3d(MPI_COMM_WORLD,species,local,local_t,cells,faces,width,ghost,gt,ok, &
    inlet,inlet_t,walls,pressure)
  call require(.not.ok.and.all(ghost==saved).and.all(gt==saved_t),'pressure disagreement rollback')
  rhs=123;diffusion=456;rate=789
  pressure=0
  if(rank==nranks-1) wall_t(3)=1200
  call compute_mpi_dynamic_les_boundary_rhs_3d(MPI_COMM_WORLD,species,local,local_t,h,cells,'box',2, &
    rhs,diffusion,ok,boundary_faces=faces,boundary_rate=rate,wall_temperature=wall_t,wall_velocity=walls)
  call require(.not.ok.and.all(rhs==123).and.diffusion==456.and.all(rate==789),'RHS disagreement rollback')
  if(rank==0) print *, 'PASS: distributed physical donors, Dynamic fluxes and SSPRK2 match serial with rollback'
  call MPI_Finalize(ierr)
contains
  subroutine check_coupled_sdc(observer)
    procedure(forced_reactor_integrator), optional :: observer
    type(gas_transport_species), allocatable :: properties(:)
    type(elementary_reaction), allocatable :: reactions(:)
    type(spray_liquid) :: liquid(1)
    type(spray_options) :: options
    type(spray_wall_options) :: settings
    type(spray_parcel), allocatable :: all_initial(:),drops(:),reference(:),expected(:)
    type(spray_wall_film), allocatable :: films(:),ref_films(:),initial_films(:),expected_films(:)
    real(dp), allocatable :: mass(:,:),ref_mass(:,:),film_mass(:,:),ref_film_mass(:,:),initial_film_mass(:,:)
    real(dp), allocatable :: state(:,:,:,:),temp(:,:,:),ir(:,:,:,:),heat(:,:,:)
    real(dp) :: ref_state(nv,5,3,2),ref_temp(5,3,2),ref_heat(5,3,2),ref_ir(nv-4,5,3,2)
    real(dp) :: spacing(3)=[.001_dp,.0015_dp,.002_dp],zero(3)=0,force(3)=[10._dp,-20._dp,30._dp]
    real(dp) :: boundary(nv),source(nv),exchange(nv),wall(5),outside(nv),forced(5)
    real(dp) :: rb(nv),rs(nv),rx(nv),rw(5),ro(nv),rf(5),wall_t(6)=400
    real(dp), parameter :: interval=1e-8_dp
    integer(int64) :: watermark,ref_id
    logical, allocatable :: owned(:)
    type(rotating_frame_3d), allocatable :: frame
    real(dp) :: base_state(nv,5,3,2),base_t(5,3,2),delta(5),global_delta(5),energy_error
    logical :: flags(4),mol,rotating
    integer :: scenario,p,block_first(3),block_extent(3)
    integer :: selected_scenario,parse_status
    character(len=16) :: scenario_arg
    real(dp) :: phase_clock
    selected_scenario=0
    call get_command_argument(2,scenario_arg)
    if(len_trim(scenario_arg)>0) then
      read(scenario_arg,*,iostat=parse_status) selected_scenario
      call require(parse_status==0.and.selected_scenario>=1.and.selected_scenario<=3,'selected coupled scenario')
    endif
    rotating=test_mode=='mol-coupled-rotation';mol=test_mode=='mol-coupled'.or.rotating
    base_state=full;base_t=temperature
    if(rotating) then
      allocate(frame);frame%axis=3;frame%angular_speed=3000;frame%axis_point=[.03_dp,-.04_dp,.02_dp]
    endif
    call load_h2o2_elementary_transport(properties,ok);call require(ok,'coupled transport')
    call initialize_spray_liquid(species,'H2O',1000._dp,4184._dp,373.15_dp,101325._dp,2.257e6_dp,liquid(1),ok)
    call require(ok,'coupled liquid')
    allocate(reactions(0),all_initial(nranks),ir(nv-4,extent(1),extent(2),extent(3)), &
      heat(extent(1),extent(2),extent(3)))
    if(rotating) then
      deallocate(reactions);call load_h2o2_elementary_mechanism(reactions,ok)
      call require(ok,'rotating coupled real chemistry')
    endif
    do scenario=1,3
      if(selected_scenario>0.and.scenario/=selected_scenario) cycle
      phase_clock=MPI_Wtime()
      if(rank==0) print *, 'Coupled scenario start:',scenario
      full=base_state;temperature=base_t
      options=spray_options();faces='periodic';flags=.false.;settings=spray_wall_options()
      options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
      if(scenario==2) options%mixture_time_scheme='euler'
      if(scenario==3) then
        faces(3:4)='isothermal_wall';settings%enabled(3)=.true.;settings%evolve=.true.
        flags(2)=.true.
        settings%temperature=400;settings%liquid_conductivity=.6_dp;settings%minimum_height=1e-15_dp
        options%drag=.false.;options%heat=.false.;options%evaporation=.false.
        ! A common 400 K wall heats the gas and the liquid films.
        q=0;q(1)=1;q(2:4)=[1._dp,1._dp,1._dp];q(5)=101325;q(12)=1
        call reactive_primitive_to_conserved(species,q,full(:,1,1,1),temperature(1,1,1),sound,ok)
        call require(ok,'wall/film initial gas thermodynamics')
        do k=1,2
          do j=1,3
            do i=1,5
              full(:,i,j,k)=full(:,1,1,1);temperature(i,j,k)=temperature(1,1,1)
            enddo
          enddo
        enddo
        local=full(:,lo:hi,:,:);local_t=temperature(lo:hi,:,:)
      endif
      if(rotating) then
        call convert_regular_energy_3d(frame,zero,spacing,full,0,.true.,ok)
        call require(ok,'coupled rotating initial energy')
      endif
      local=full(:,lo:hi,:,:);local_t=temperature(lo:hi,:,:)
      do p=1,nranks
        call spray_cartesian_block(cells,[nranks,1,1],p-1,block_first,block_extent,ok)
        call require(ok,'global serial oracle particle origin')
        all_initial(p)=spray_parcel(int(p,int64), &
          [(real(block_first(1)+block_extent(1),dp)-.00005_dp)*spacing(1),.5_dp*spacing(2),.5_dp*spacing(3)], &
          [20._dp,1._dp,0._dp],1e-12_dp,330._dp,10._dp)
        if(scenario==3) then
          all_initial(p)%position(2)=.1_dp*interval;all_initial(p)%velocity=[0._dp,-1._dp,0._dp]
        endif
      enddo
      reference=all_initial;ref_mass=reshape(reference%mass,[1,nranks])
      drops=[all_initial(rank+1)];mass=reshape(drops%mass,[1,1]);state=local;temp=local_t
      ref_state=full;ref_temp=temperature;ref_heat=53;ref_ir=0;ir=0;heat=53;watermark=int(nranks,int64);ref_id=watermark
      profile_phase=1
      if(scenario<3) then
        if(mol) then
          call advance_mpi_cartesian_spray_sdc(MPI_COMM_SELF,species,reactions,properties,liquid,options, &
            ref_state,ref_temp,reference,ref_mass,zero,spacing,cells,interval,0._dp,2,force,flags,1e-6_dp,1e-12_dp,.1_dp, &
            faces,ref_ir,ref_heat,ref_id,rb,rs,rx,rw,ro,rf,ok,acceleration=force,mol_mode=.true., &
            mol_lax_friedrichs=scenario==2,frame=frame,evaluate_reaction=observer)
        else
        call advance_spray_reactive_sdc_3d(species,reactions,properties,liquid,options,ref_state,ref_temp, &
          reference,ref_mass,zero,spacing,interval,0._dp,2,.false.,flags,force,0._dp,.7_dp,.7_dp,1e-6_dp,1e-12_dp, &
          faces,rb,rs,rx,rw,ro,rf,ok,acceleration=force,last_id=ref_id,difmag=.1_dp,heat_release=ref_heat)
        endif
        if(.not.ok) print *, 'Serial coupled oracle mode/scenario: ',test_mode,scenario
        call require(ok,'serial fully coupled SDC oracle')
        if(rank==0) print *, 'Coupled serial oracle seconds:',scenario,MPI_Wtime()-phase_clock
        phase_clock=MPI_Wtime();profile_phase=2
        call advance_mpi_cartesian_spray_sdc(MPI_COMM_WORLD,species,reactions,properties,liquid,options, &
          state,temp,drops,mass,zero,spacing,cells,interval,0._dp,2,force,flags,1e-6_dp,1e-12_dp,.1_dp, &
          faces,ir,heat,watermark,boundary,source,exchange,wall,outside,forced,ok,acceleration=force, &
          mol_mode=mol,mol_lax_friedrichs=mol.and.scenario==2,frame=frame,evaluate_reaction=observer)
      else
        allocate(ref_films(nranks),ref_film_mass(1,nranks))
        do p=1,nranks
          call spray_cartesian_block(cells,[nranks,1,1],p-1,block_first,block_extent,ok)
          ref_films(p)=spray_wall_film([(real(block_first(1),dp)+.5_dp)*spacing(1),0._dp,.5_dp*spacing(3)], &
            [0._dp,1._dp,0._dp],1e-12_dp,330._dp,1e-6_dp,1e-9_dp,1)
          ref_film_mass(1,p)=ref_films(p)%mass
        enddo
        films=[ref_films(rank+1)];film_mass=reshape(films%mass,[1,1])
        initial_films=films;initial_film_mass=film_mass
        if(mol) then
          call advance_mpi_cartesian_spray_sdc(MPI_COMM_SELF,species,reactions,properties,liquid,options, &
            ref_state,ref_temp,reference,ref_mass,zero,spacing,cells,interval,0._dp,2,force,flags,1e-6_dp,1e-12_dp,.1_dp, &
            faces,ref_ir,ref_heat,ref_id,rb,rs,rx,rw,ro,rf,ok,acceleration=force,mol_mode=.true., &
            wall_options=settings,wall_films=ref_films,film_component_mass=ref_film_mass,wall_temperature=wall_t,frame=frame)
        else
        call advance_spray_reactive_sdc_3d(species,reactions,properties,liquid,options,ref_state,ref_temp, &
          reference,ref_mass,zero,spacing,interval,0._dp,2,.false.,flags,force,0._dp,.7_dp,.7_dp,1e-6_dp,1e-12_dp, &
          faces,rb,rs,rx,rw,ro,rf,ok,acceleration=force,last_id=ref_id,difmag=.1_dp,heat_release=ref_heat, &
          wall_options=settings,wall_films=ref_films,film_component_mass=ref_film_mass,wall_temperature=wall_t)
        endif
        call require(ok,'serial fully coupled wall/film SDC oracle')
        call advance_mpi_cartesian_spray_sdc(MPI_COMM_WORLD,species,reactions,properties,liquid,options, &
          state,temp,drops,mass,zero,spacing,cells,interval,0._dp,2,force,flags,1e-6_dp,1e-12_dp,.1_dp, &
          faces,ir,heat,watermark,boundary,source,exchange,wall,outside,forced,ok,acceleration=force, &
          wall_options=settings,wall_films=films,film_component_mass=film_mass,wall_temperature=wall_t,mol_mode=mol,frame=frame)
      endif
      call require(ok,'distributed fully coupled SDC accepted')
      call require(maxval(abs(state-ref_state(:,lo:hi,:,:)))<2e-12_dp*maxval(abs(ref_state)).and. &
        maxval(abs(temp-ref_temp(lo:hi,:,:)))<2e-8_dp,'coupled gas and temperature match global serial oracle')
      if(.not.rotating) call require(all(ir==0).and.all(heat==0),'coupled nonreacting IR and heat reset')
      if(rotating) then
        call require(maxval(abs(heat-ref_heat(lo:hi,:,:)))<1e-12_dp,'coupled rotating reaction heat agreement')
        delta=0
        do p=1,5
          delta(p)=product(spacing)*sum(state(p,:,:,:)-local(p,:,:,:))
        enddo
        delta=delta-parcel_totals(all_initial(rank+1),liquid(1))
        do p=1,size(drops)
          delta=delta+parcel_totals(drops(p),liquid(1))
        enddo
        if(scenario<3) then
          call MPI_Allreduce(delta,global_delta,5,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
          global_delta=global_delta-boundary(:5)-source(:5)-forced+outside(:5)
          energy_error=abs(global_delta(5))
          call require(energy_error<2e-16_dp,'rotating gas/particle mixed-energy inventory including Phi mass source')
          if(rank==0) print *, 'Rotating coupled MOL scenario, energy inventory residual:',scenario,energy_error
        endif
      endif
      call require(maxval(abs(boundary-rb))<max(1e-20_dp,2e-10_dp*maxval(abs(rb))).and. &
        maxval(abs(source-rs))<max(1e-20_dp,2e-10_dp*maxval(abs(rs))).and. &
        maxval(abs(exchange-rx))<max(1e-20_dp,2e-10_dp*maxval(abs(rx))), &
        'coupled hydro, external and particle exchange ledgers match serial')
      call require(maxval(abs(wall-rw))<max(1e-20_dp,2e-10_dp*maxval(abs(rw))).and. &
        maxval(abs(outside-ro))<max(1e-20_dp,2e-10_dp*maxval(abs(ro))).and. &
        maxval(abs(forced-rf))<max(1e-20_dp,2e-10_dp*maxval(abs(rf))), &
        'coupled particle wall, outflow and force ledgers match serial')
      owned=spread(.false.,1,size(reference))
      do p=1,size(reference)
        owned(p)=spray_cartesian_owner(reference(p)%position,zero,spacing,cells,[nranks,1,1], &
          [.true.,.true.,.true.])==rank
      enddo
      expected=pack(reference,owned)
      call require(size(drops)==size(expected).and.watermark==ref_id,'coupled migrated population and ID oracle')
      do p=1,size(expected)
        call require(drops(p)%id==expected(p)%id.and.abs(drops(p)%mass-expected(p)%mass)<1e-24_dp.and. &
          abs(drops(p)%temperature-expected(p)%temperature)<1e-10_dp.and. &
          maxval(abs(drops(p)%velocity-expected(p)%velocity))<1e-10_dp.and. &
          maxval(abs(drops(p)%position-expected(p)%position))<1e-14_dp,'coupled particle physical oracle')
      enddo
      if(scenario==3) then
        owned=spread(.false.,1,size(ref_films))
        do p=1,size(ref_films)
          owned(p)=spray_cartesian_owner(ref_films(p)%position,zero,spacing,cells,[nranks,1,1], &
            [.true.,.true.,.true.])==rank
        enddo
        expected_films=pack(ref_films,owned)
        call require(size(films)==size(expected_films),'coupled film count oracle')
        do p=1,size(films)
          call require(abs(films(p)%mass-expected_films(p)%mass)<1e-24_dp.and. &
            abs(films(p)%temperature-expected_films(p)%temperature)<1e-10_dp,'film old/final kicks match serial')
        enddo
        ! One iteration reaches final particle/film candidates before the
        ! rank-zero chemistry callback rejects gas publication.
        deallocate(reactions);allocate(reactions(1));drops=[all_initial(rank+1)]
        mass=reshape(drops%mass,[1,1]);films=initial_films;film_mass=initial_film_mass
        state=local;temp=local_t;ir=0;heat=53;watermark=int(nranks,int64)
        boundary=41;source=42;exchange=43;wall=44;outside=45;forced=46;reactor_calls=0
        call advance_mpi_cartesian_spray_sdc(MPI_COMM_WORLD,species,reactions,properties,liquid,options, &
          state,temp,drops,mass,zero,spacing,cells,interval,0._dp,1,force,flags,1e-6_dp,1e-12_dp,.1_dp, &
          faces,ir,heat,watermark,boundary,source,exchange,wall,outside,forced,ok,acceleration=force, &
          wall_options=settings,wall_films=films,film_component_mass=film_mass,wall_temperature=wall_t, &
          evaluate_reaction=rank_reactor,mol_mode=mol,frame=frame)
        call require(.not.ok.and.reactor_calls>0,'late reacting gas rejection after final particle/film stages')
        call require(all(state==local).and.all(temp==local_t).and.all(ir==0).and.all(heat==53).and. &
          watermark==nranks.and.size(drops)==1.and.all(mass==all_initial(rank+1)%mass).and. &
          drops(1)%mass==all_initial(rank+1)%mass.and. &
          all(drops(1)%position==all_initial(rank+1)%position),'whole coupled gas/particle/ID rollback')
        call require(all(film_mass==initial_film_mass).and.films(1)%temperature==initial_films(1)%temperature.and. &
          all(boundary==41).and.all(source==42).and.all(exchange==43).and.all(wall==44).and. &
          all(outside==45).and.all(forced==46),'whole coupled film and ledger rollback')
      endif
    enddo
    if(rank==0) print *, 'PASS: coupled MPI gas/particle/film and late chemistry rollback, MOL=',mol
  end subroutine check_coupled_sdc

  subroutine check_sdc_particle_stage()
    type(gas_transport_species), allocatable :: properties(:)
    type(spray_liquid) :: liquid(1)
    type(spray_options) :: options
    type(spray_wall_options) :: settings
    type(spray_wall_film), allocatable :: films(:),reference_films(:)
    real(dp), allocatable :: film_mass(:,:),reference_film_mass(:,:)
    type(spray_parcel), allocatable :: drops(:),reference(:),oracle(:),initial(:)
    type(spray_sdc_stage_3d) :: serial_work
    type(mpi_sdc_particle_stage) :: work,kept_work
    real(dp), allocatable :: mass(:,:),ref_mass(:,:),oracle_mass(:,:),initial_mass(:,:),rate(:,:,:,:),kept_rate(:,:,:,:)
    real(dp) :: spacing(3)=[.001_dp,.0015_dp,.002_dp],zero(3)=0,accel(3)=[100._dp,-200._dp,300._dp]
    real(dp) :: wall(5),outside(nv),forced(5),ref_wall(5),ref_out(nv),ref_force(5)
    real(dp) :: global_wall(5),global_out(nv),global_force(5),serial_rhs(nv,5,3,2),global_rhs(nv,5,3,2)
    real(dp), parameter :: interval=1e-7_dp
    integer(int64) :: watermark,serial_watermark
    integer :: scenario,stage,p
    call load_h2o2_elementary_transport(properties,ok);call require(ok,'stage transport')
    call initialize_spray_liquid(species,'H2O',1000._dp,4184._dp,373.15_dp,101325._dp,2.257e6_dp,liquid(1),ok)
    call require(ok,'stage liquid')
    allocate(initial(1),initial_mass(1,1),rate(nv,extent(1),extent(2),extent(3)))
    do scenario=1,4
      faces='periodic';options=spray_options()
      options%drag=.false.;options%heat=.false.;options%evaporation=.false.
      if(scenario==2) options%mixture_time_scheme='euler'
      if(scenario==3) options%fixed_parts=.true.
      if(scenario==4) then
        options%maximum_multiplicity=60;options%breakup_model='etab'
        options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
      endif
      if(scenario<4) then
        options%drag=.true.;options%heat=.true.;options%evaporation=.true.
      endif
      initial(1)=spray_parcel(int(rank+1,int64),[(real(first(1)+extent(1),dp)-.1_dp)*spacing(1), &
        .5_dp*spacing(2),.5_dp*spacing(3)],[4000._dp,1._dp,0._dp],1e-12_dp,330._dp,100._dp)
      if(scenario==4) initial(1)%velocity=[20._dp,1._dp,0._dp]
      initial_mass(1,1)=initial(1)%mass
      drops=initial;mass=initial_mass;reference=initial;ref_mass=initial_mass
      serial_work=spray_sdc_stage_3d();work=mpi_sdc_particle_stage()
      watermark=int(nranks,int64);serial_watermark=watermark
      do stage=0,2
        call advance_spray_sdc_stage_3d(species,properties,liquid,options,full,temperature,zero,spacing, &
          0._dp,interval,stage,2,faces,reference,ref_mass,serial_work,serial_rhs,ref_wall,ref_out,ref_force,ok, &
          accel,serial_watermark,native_flow_cfl=1._dp)
        if(.not.ok) print *, 'serial stage rejected scenario/stage/rank:',scenario,stage,rank
        call require(ok,'serial complete particle stage oracle')
        call MPI_Allreduce(serial_rhs,global_rhs,size(serial_rhs),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(ref_wall,global_wall,5,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(ref_out,global_out,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(ref_force,global_force,5,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call advance_mpi_sdc_particle_stage(MPI_COMM_WORLD,species,properties,liquid,options,local,local_t, &
          zero,spacing,cells,[nranks,1,1],0._dp,interval,stage,2,faces,drops,mass,work,rate,wall,outside, &
          forced,watermark,ok,accel,native_flow_cfl=1._dp)
        call require(ok,'complete collective particle stage accepted')
        call require(maxval(abs(rate-global_rhs(:,lo:hi,:,:)))<1e-12_dp*max(1._dp,maxval(abs(global_rhs))), &
          'stage source matches serial predictor/final kick')
        call require(maxval(abs(wall-global_wall))<1e-18_dp.and.maxval(abs(outside-global_out))<1e-18_dp.and. &
          maxval(abs(forced-global_force))<1e-17_dp,'stage global ledgers match serial')
        if(stage<2) then
          call require(all(mass==initial_mass).and.watermark==nranks.and.all(drops(1)%position==initial(1)%position), &
            'predictor and repeated correction keep public particles private')
          oracle=serial_work%midpoint;oracle_mass=serial_work%midpoint_mass
        else
          oracle=reference;oracle_mass=ref_mass
        endif
        ! Computational split IDs use global reservation in MPI and local
        ! reservation in the serial oracle. Compare physical parcels below.
        if(scenario<4.or.stage<2) then
          call migrate_spray_parcels(MPI_COMM_WORLD,oracle,zero,spacing,cells,[nranks,1,1], &
            [.true.,.true.,.true.],100,ok,oracle_mass)
          call require(ok,'stage serial owner redistribution')
          if(stage<2) then
            call require(size(work%midpoint)==size(oracle),'midpoint owner count')
            do p=1,size(oracle)
              call require(work%midpoint(p)%id==oracle(p)%id.and. &
                maxval(abs(work%midpoint(p)%position-oracle(p)%position))<1e-14_dp.and. &
                maxval(abs(work%midpoint(p)%velocity-oracle(p)%velocity))<1e-10_dp,'midpoint physical oracle')
            enddo
          else
            call require(size(drops)==size(oracle),'final owner count')
            do p=1,size(oracle)
              call require(drops(p)%id==oracle(p)%id.and.maxval(abs(drops(p)%position-oracle(p)%position))<1e-14_dp.and. &
                maxval(abs(drops(p)%velocity-oracle(p)%velocity))<1e-10_dp,'final physical oracle')
            enddo
          endif
        else
          call require(size(drops)==2.and.watermark==2*nranks,'global final computational split count and IDs')
          call require(all(drops%multiplicity==50._dp).and.all(drops%mass==1e-12_dp), &
            'computational split preserves physical parcel mass')
        endif
        if(stage==1) then
          kept_work=work;kept_rate=rate
          options%heat=.not.options%heat
          call advance_mpi_sdc_particle_stage(MPI_COMM_WORLD,species,properties,liquid,options,local,local_t, &
            zero,spacing,cells,[nranks,1,1],0._dp,interval,2,2,faces,drops,mass,work,rate,wall,outside, &
            forced,watermark,ok,accel,native_flow_cfl=1._dp)
          call require(.not.ok.and.work%next_stage==kept_work%next_stage.and.all(rate==kept_rate).and. &
            all(mass==initial_mass).and.watermark==nranks,'mid-interval material change rollback')
          options%heat=.not.options%heat
          local(1,1,1,1)=-1
          call advance_mpi_sdc_particle_stage(MPI_COMM_WORLD,species,properties,liquid,options,local,local_t, &
            zero,spacing,cells,[nranks,1,1],0._dp,interval,2,2,faces,drops,mass,work,rate,wall,outside, &
            forced,watermark,ok,accel,native_flow_cfl=1._dp)
          call require(.not.ok.and.work%next_stage==2.and.all(rate==kept_rate).and.all(mass==initial_mass), &
            'final sampled gas failure preserves stage and public state')
          local=full(:,lo:hi,:,:)
        endif
      enddo
      call require(work%next_stage==0,'accepted final stage resets lifecycle')
    enddo
    ! Dry wall deposits remain private together with the predictor particles.
    options=spray_options();options%drag=.false.;options%heat=.false.;options%evaporation=.false.
    options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
    faces='periodic';faces(3:4)='slip_wall';settings%enabled(3)=.true.
    initial(1)=spray_parcel(int(rank+1,int64),[(real(first(1),dp)+.5_dp)*spacing(1), &
      .1_dp*interval,.5_dp*spacing(3)],[0._dp,-1._dp,0._dp],1e-12_dp,330._dp,100._dp)
    initial_mass(1,1)=initial(1)%mass;drops=initial;mass=initial_mass;reference=initial;ref_mass=initial_mass
    allocate(films(0),film_mass(1,0),reference_films(0),reference_film_mass(1,0))
    work=mpi_sdc_particle_stage();serial_work=spray_sdc_stage_3d()
    watermark=int(nranks,int64);serial_watermark=watermark
    do stage=0,2
      call advance_spray_sdc_stage_3d(species,properties,liquid,options,full,temperature,zero,spacing, &
        0._dp,interval,stage,2,faces,reference,ref_mass,serial_work,serial_rhs,ref_wall,ref_out,ref_force,ok, &
        last_id=serial_watermark,wall_options=settings,wall_films=reference_films,film_component_mass=reference_film_mass)
      if(.not.ok) print *, 'serial deposited film rejected stage/rank:',stage,rank
      call require(ok,'serial deposited film stage oracle')
      call advance_mpi_sdc_particle_stage(MPI_COMM_WORLD,species,properties,liquid,options,local,local_t, &
        zero,spacing,cells,[nranks,1,1],0._dp,interval,stage,2,faces,drops,mass,work,rate,wall,outside,forced, &
        watermark,ok,wall_options=settings,wall_films=films,film_component_mass=film_mass)
      call require(ok,'collective deposited film stage accepted')
      call MPI_Allreduce(ref_wall,global_wall,5,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call require(maxval(abs(wall-global_wall))<1e-18_dp,'whole-stage deposit wall ledger')
      if(stage<2) then
        call require(size(films)==0.and.size(drops)==1.and.size(work%midpoint_films)==1.and. &
          size(work%midpoint)==0,'deposits and removed parent stay private through correction')
      else
        call require(size(films)==size(reference_films).and.size(drops)==0.and.all(film_mass==reference_film_mass), &
          'final deposit and removed parent commit together')
        call require(all(films(1)%position==reference_films(1)%position).and. &
          films(1)%mass==reference_films(1)%mass,'committed film physical oracle')
      endif
    enddo
    if(rank==0) print *, 'PASS: MPI particle MKD/MK lifecycle, force, split and whole-stage rollback'
  end subroutine check_sdc_particle_stage

  subroutine check_sdc_drift()
    type(gas_transport_species), allocatable :: properties(:)
    type(spray_liquid) :: liquid(1)
    type(spray_options) :: options,materials(2)
    type(spray_wall_options) :: settings
    type(spray_parcel), allocatable :: drops(:),reference(:),kept(:)
    type(spray_wall_film), allocatable :: films(:),ref_films(:)
    type(spray_sdc_stage_3d) :: serial_work
    real(dp), allocatable :: mass(:,:),ref_mass(:,:),history(:),kept_mass(:,:),kept_history(:),fm(:,:),ref_fm(:,:)
    real(dp) :: spacing(3)=[.001_dp,.0015_dp,.002_dp],zero(3)=0,wall(5),outside(nv),ref_wall(5),ref_out(nv),forced(5)
    real(dp) :: global_wall(5),global_out(nv),serial_rhs(nv,5,3,2)
    real(dp), parameter :: interval=1e-7_dp
    integer(int64) :: watermark,serial_watermark
    integer :: scenario,p
    call load_h2o2_elementary_transport(properties,ok);call require(ok,'drift transport')
    call initialize_spray_liquid(species,'H2O',1000._dp,4184._dp,373.15_dp,101325._dp,2.257e6_dp,liquid(1),ok)
    call require(ok,'drift liquid')
    options%drag=.false.;options%heat=.false.;options%evaporation=.false.
    options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
    materials=options;materials(2)%surface_tension=.05_dp;materials(2)%liquid_viscosity=.002_dp
    allocate(drops(1),mass(1,1),history(1),films(0),fm(1,0),ref_films(0),ref_fm(1,0))
    do scenario=1,4
      faces='periodic';settings%enabled=.false.
      if(scenario==2) faces(1:2)='slip_wall'
      if(scenario==3) faces(1:2)='outflow'
      if(scenario==4) then
        faces(3:4)='slip_wall';settings%enabled(3)=.true.
      endif
      drops=[spray_parcel(int(rank+1,int64),[(real(first(1)+extent(1),dp)-.1_dp)*spacing(1), &
        .5_dp*spacing(2),.5_dp*spacing(3)],[4000._dp,0._dp,0._dp],1e-12_dp,330._dp,100._dp)]
      drops(1)%material_id=1+modulo(rank,2)
      if(scenario==4) then
        drops(1)%position(2)=.1_dp*interval;drops(1)%velocity=[0._dp,-1._dp,0._dp]
      endif
      mass=reshape([drops(1)%mass],[1,1]);history=[.25_dp*real(drops(1)%id,dp)]
      reference=drops;ref_mass=mass;serial_work=spray_sdc_stage_3d()
      films=films(:0);fm=fm(:,:0);ref_films=ref_films(:0);ref_fm=ref_fm(:,:0)
      watermark=int(nranks,int64);serial_watermark=watermark
      if(scenario==4) then
        call advance_spray_sdc_stage_3d(species,properties,liquid,options,full,temperature,zero,spacing,0._dp,interval, &
          0,2,faces,reference,ref_mass,serial_work,serial_rhs,ref_wall,ref_out,forced,ok,last_id=serial_watermark, &
          breakup_properties=materials,wall_options=settings,wall_films=ref_films,film_component_mass=ref_fm)
      else
        call advance_spray_sdc_stage_3d(species,properties,liquid,options,full,temperature,zero,spacing,0._dp,interval, &
          0,2,faces,reference,ref_mass,serial_work,serial_rhs,ref_wall,ref_out,forced,ok,last_id=serial_watermark)
      endif
      call require(ok,'serial ballistic/impact predictor oracle')
      reference=pack(serial_work%midpoint,serial_work%midpoint%mass>0)
      ref_mass=reshape(pack(serial_work%midpoint_mass(1,:),serial_work%midpoint%mass>0),[1,size(reference)])
      call migrate_spray_parcels(MPI_COMM_WORLD,reference,zero,spacing,cells,[nranks,1,1], &
        [faces(1)=='periodic',faces(3)=='periodic',faces(5)=='periodic'],100,ok,ref_mass)
      call require(ok,'serial oracle owner redistribution')
      call MPI_Allreduce(ref_wall,global_wall,5,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(ref_out,global_out,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      wall=71;outside=72
      if(scenario==4) then
        call drift_mpi_sdc_particles(MPI_COMM_WORLD,species,liquid,options,local,local_t,drops,mass,zero,spacing, &
          cells,[nranks,1,1],faces,interval,watermark,wall,outside,ok,history,settings,films,fm,materials)
        call require(ok,'transactional wall impact drift accepted')
        call migrate_mpi_wall_films(MPI_COMM_WORLD,ref_films,ref_fm,zero,spacing,cells,[nranks,1,1],100,ok)
        call require(ok,'serial impact film owner redistribution')
        call require(size(films)==size(ref_films).and.all(shape(fm)==shape(ref_fm)),'film migration oracle shape')
        do p=1,size(films)
          call require(abs(films(p)%mass-ref_films(p)%mass)<1e-24_dp.and. &
            maxval(abs(films(p)%position-ref_films(p)%position))<1e-14_dp,'owner-local impact film oracle')
        enddo
      else
        call drift_mpi_sdc_particles(MPI_COMM_WORLD,species,liquid,options,local,local_t,drops,mass,zero,spacing, &
          cells,[nranks,1,1],faces,interval,watermark,wall,outside,ok,history)
      endif
      call require(ok.and.watermark==nranks.and.size(drops)==size(reference),'transactional drift owner population')
      do p=1,size(drops)
        call require(drops(p)%id==reference(p)%id.and.maxval(abs(drops(p)%position-reference(p)%position))<1e-14_dp.and. &
          all(drops(p)%velocity==reference(p)%velocity),'drift/reflection/outflow matches serial predictor')
        call require(history(p)==.25_dp*real(drops(p)%id,dp),'TAB tangent follows ID through partition and sorted exchange')
      enddo
      call require(all(mass==ref_mass),'component mass follows moved and surviving particles')
      call require(maxval(abs(wall-global_wall))<1e-17_dp.and.maxval(abs(outside-global_out))<1e-17_dp, &
        'global drift wall and component outflow ledgers match serial')
    enddo
    ! Every rank starts within its local capacity; all particles arriving on
    ! rank zero exceed that rank's capacity after the completed ballistic step.
    faces='periodic';options%maximum_parcels=1
    drops=[spray_parcel(int(rank+1,int64),[.0005_dp,.0005_dp,.0005_dp],zero,1e-12_dp,330._dp,100._dp)]
    mass=reshape([drops(1)%mass],[1,1]);history=[real(rank+1,dp)]
    kept=drops;kept_mass=mass;kept_history=history;wall=81;outside=82;watermark=int(nranks,int64)
    call drift_mpi_sdc_particles(MPI_COMM_WORLD,species,liquid,options,local,local_t,drops,mass,zero,spacing, &
      cells,[nranks,1,1],faces,interval,watermark,wall,outside,ok,history)
    call require(.not.ok.and.watermark==nranks.and.all(wall==81).and.all(outside==82).and.all(mass==kept_mass).and. &
      all(history==kept_history).and.all(drops(1)%position==kept(1)%position),'late migration capacity rejects whole drift')
    call require(all(local==full(:,lo:hi,:,:)).and.all(local_t==temperature(lo:hi,:,:)),'drift never modifies gas stage')
    if(rank==0) print *, 'PASS: MPI drift, wall/outflow budgets, film migration and TAB payload transaction'
  end subroutine check_sdc_drift

  subroutine check_sdc_breakup()
    type(gas_transport_species), allocatable :: properties(:)
    type(spray_liquid) :: liquid(1)
    type(spray_options) :: options,materials(2)
    type(spray_parcel) :: initial
    type(spray_parcel), allocatable :: drops(:),reference(:),kept(:)
    type(spray_sdc_stage_3d) :: serial_work
    real(dp), allocatable :: mass(:,:),reference_mass(:,:),sample(:,:),sample_t(:),reynolds(:),rhs(:,:,:,:)
    real(dp), allocatable :: kept_mass(:,:),kept_rhs(:,:,:,:),tangent(:),cycle_rhs(:,:,:,:)
    real(dp) :: spacing(3)=[.001_dp,.0015_dp,.002_dp],zero(3)=0,force(5),radius
    real(dp) :: serial_rhs(nv,5,3,2),global_rhs(nv,5,3,2),wall(5),outside(nv),serial_force(5),position(3)
    real(dp) :: before(5),after(5),delta(nv),global_delta(nv)
    real(dp), parameter :: interval=5e-8_dp
    integer(int64) :: watermark,serial_watermark,births,total_births,kept_watermark
    integer :: variant,p,c,flow_cycle
    call load_h2o2_elementary_transport(properties,ok);call require(ok,'breakup transport')
    call initialize_spray_liquid(species,'H2O',1000._dp,4184._dp,373.15_dp,101325._dp,2.257e6_dp,liquid(1),ok)
    call require(ok,'breakup liquid')
    faces='periodic'
    allocate(drops(1),mass(1,1),reference(1),reference_mass(1,1),rhs(nv,extent(1),extent(2),extent(3)))
    do variant=1,4
      options=spray_options();options%heat=.false.;options%evaporation=.false.
      options%breakup_model='etab';options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
      if(variant==2) then
        options%etab_parcel_factor=1;options%maximum_multiplicity=3
      endif
      if(variant>=3) then
        options%breakup_model='khrt';options%khrt_children=.true.
        options%khrt_reference_diameter_cubed=1e-25_dp;options%drag=.false.
        if(variant==4) options%khrt_parcel_factor=1
      endif
      materials=options;materials(2)%surface_tension=.05_dp;materials(2)%liquid_viscosity=.002_dp
      initial=spray_parcel(int(rank+1,int64),[(real(first(1)+extent(1),dp)-.25_dp)*spacing(1), &
        1.5_dp*spacing(2),spacing(3)],[1000._dp,0._dp,0._dp],1e-12_dp,350._dp,10._dp)
      initial%deformation=1.1_dp;initial%material_id=1+modulo(rank,2)
      if(variant>=3) then
        initial%mass=1e-9_dp;initial%multiplicity=.01_dp;initial%velocity=[200._dp,0._dp,0._dp]
        radius=(3*initial%mass/(4*acos(-1._dp)*liquid(1)%density))**(1._dp/3._dp)
        initial%khrt_initial_count=initial%multiplicity
        initial%khrt_shed_volume=.95_dp*initial%multiplicity*radius**3
      endif
      drops=[initial];mass=reshape([initial%mass],[1,1]);reference=drops;reference_mass=mass
      serial_work=spray_sdc_stage_3d();watermark=int(nranks,int64);serial_watermark=watermark
      call advance_spray_sdc_stage_3d(species,properties,liquid,options,full,temperature,zero,spacing,0._dp,interval, &
        0,2,faces,reference,reference_mass,serial_work,serial_rhs,wall,outside,serial_force,ok, &
        last_id=serial_watermark,breakup_properties=materials)
      call require(ok,'serial transfer breakup predictor oracle')
      births=serial_work%midpoint_watermark-watermark
      call MPI_Allreduce(births,total_births,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(serial_rhs,global_rhs,size(serial_rhs),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call kick_mpi_sdc_transfer(MPI_COMM_WORLD,species,properties,liquid,options,local,local_t,drops,mass,zero, &
        spacing,cells,[nranks,1,1],faces,interval,zero,rhs,force,ok,sample,sample_t,reynolds,materials)
      call require(ok,'half transfer exports its sample and final internal Reynolds')
      kept=drops;kept_mass=mass;kept_rhs=rhs;kept_watermark=watermark
      if(rank==nranks-1) sample(1,1)=-1
      call kick_mpi_sdc_breakup(MPI_COMM_WORLD,species,properties,liquid,options,drops,mass,sample,sample_t,reynolds, &
        zero,spacing,cells,[nranks,1,1],faces,interval,.5_dp*interval,.5_dp*interval,.true.,.true., &
        watermark,rhs,ok,materials)
      call require(.not.ok.and.watermark==kept_watermark.and.all(rhs==kept_rhs).and.all(mass==kept_mass).and. &
        all(drops(1)%velocity==kept(1)%velocity).and.drops(1)%mass==kept(1)%mass, &
        'one-rank breakup rejection retains transfer candidate and all birth IDs')
      sample(1,1)=.3_dp
      call kick_mpi_sdc_breakup(MPI_COMM_WORLD,species,properties,liquid,options,drops,mass,sample,sample_t,reynolds, &
        zero,spacing,cells,[nranks,1,1],faces,interval,.5_dp*interval,.5_dp*interval,.true.,.true., &
        watermark,rhs,ok,materials)
      call require(ok.and.watermark==int(nranks,int64)+total_births,'collective breakup ID watermark')
      call require(size(drops)==size(serial_work%midpoint),'breakup retained and born population matches serial')
      if(variant==2) call require(total_births>0,'ETAB replacement creates children')
      if(variant>=3) call require(size(drops)>1.and.drops(1)%mass==0,'KHRT release creates children')
      do p=1,size(drops)
        call require(abs(drops(p)%mass-serial_work%midpoint(p)%mass)<1e-24_dp.and. &
          abs(drops(p)%multiplicity-serial_work%midpoint(p)%multiplicity)<1e-12_dp.and. &
          maxval(abs(drops(p)%velocity-serial_work%midpoint(p)%velocity))<1e-10_dp, &
          'material-specific breakup and sampled Reynolds match serial predictor')
        position=serial_work%midpoint(p)%position
        if(serial_work%midpoint(p)%mass>0) position=position-interval*serial_work%midpoint(p)%velocity
        call require(maxval(abs(drops(p)%position-position))<1e-14_dp,'breakup children remain at pre-drift birth location')
        if(drops(p)%parent_id/=0) call require(drops(p)%parent_id==initial%id.and.drops(p)%id>nranks, &
          'replacement children retain parent identity with global IDs')
      enddo
      call require(maxval(abs(rhs-global_rhs(:,lo:hi,:,:)))<1e-12_dp*max(1._dp,maxval(abs(global_rhs))), &
        'transfer plus breakup containing-cell source matches serial predictor')
      before=parcel_totals(initial,liquid(1));after=0
      do p=1,size(drops)
        after=after+parcel_totals(drops(p),liquid(1))
      enddo
      do c=1,nv
        delta(c)=sum(rhs(c,:,:,:))*product(spacing)*.5_dp*interval
      enddo
      delta(:5)=delta(:5)+after-before
      call MPI_Allreduce(delta,global_delta,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call require(maxval(abs(global_delta))<1e-17_dp,'breakup gas compensation closes physical inventory')
      ! A valid private midpoint with the same release trigger exercises MK:
      ! size evolves, but neither ETAB replacement nor KHRT births may occur.
      drops=[initial];mass=reshape([initial%mass],[1,1]);reference=drops;reference_mass=mass
      watermark=int(nranks,int64);serial_watermark=watermark;serial_work=spray_sdc_stage_3d()
      serial_work%next_stage=2;serial_work%iterations=2;serial_work%interval=interval
      serial_work%midpoint=reference;serial_work%midpoint_mass=reference_mass
      serial_work%old_rate=serial_rhs;serial_work%old_rate=0
      serial_work%boundary_change=spread(0._dp,1,nv)
      serial_work%has_watermark=.true.;serial_work%initial_watermark=watermark;serial_work%midpoint_watermark=watermark
      serial_work%breakup_properties=materials
      call advance_spray_sdc_stage_3d(species,properties,liquid,options,full,temperature,zero,spacing,0._dp,interval, &
        2,2,faces,reference,reference_mass,serial_work,serial_rhs,wall,outside,serial_force,ok, &
        last_id=serial_watermark,breakup_properties=materials)
      call require(ok.and.serial_watermark==watermark.and.size(reference)==1,'serial no-birth MK oracle')
      call MPI_Allreduce(serial_rhs,global_rhs,size(serial_rhs),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call kick_mpi_sdc_transfer(MPI_COMM_WORLD,species,properties,liquid,options,local,local_t,drops,mass,zero, &
        spacing,cells,[nranks,1,1],faces,interval,zero,rhs,force,ok,sample,sample_t,reynolds,materials)
      call require(ok,'MK transfer sample export')
      call kick_mpi_sdc_breakup(MPI_COMM_WORLD,species,properties,liquid,options,drops,mass,sample,sample_t,reynolds, &
        zero,spacing,cells,[nranks,1,1],faces,interval,.5_dp*interval,.5_dp*interval,.false.,.true., &
        watermark,rhs,ok,materials)
      call require(ok.and.size(drops)==1.and.watermark==nranks,'MK suppresses replacement IDs and births')
      call require(abs(drops(1)%mass-reference(1)%mass)<1e-24_dp.and. &
        abs(drops(1)%multiplicity-reference(1)%multiplicity)<1e-12_dp.and. &
        maxval(abs(drops(1)%velocity-reference(1)%velocity))<1e-10_dp.and. &
        all(drops(1)%position==initial%position),'MK breakup updates size without drift')
      call require(maxval(abs(rhs-global_rhs(:,lo:hi,:,:)))<1e-12_dp*max(1._dp,maxval(abs(global_rhs))), &
        'no-birth MK compensation matches serial final half')
    enddo
    ! Two flow subcycles carry TAB's cumulative tangent along the trajectory.
    ! The reference stage sums both TAB contributions
    ! before replacement; the MPI payload must survive the first subcycle.
    options=spray_options();options%drag=.false.;options%heat=.false.;options%evaporation=.false.
    options%breakup_model='etab';options%etab_parcel_factor=1
    options%maximum_multiplicity=3;options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
    materials=options;materials(2)%surface_tension=.05_dp;materials(2)%liquid_viscosity=.002_dp
    initial%mass=1e-12_dp;initial%multiplicity=10;initial%deformation=1.1_dp
    initial%deformation_rate=0;initial%khrt_initial_count=0;initial%khrt_shed_volume=0
    initial%rt_time=-1;initial%velocity=[1000._dp,0._dp,0._dp]
    drops=[initial];mass=reshape([initial%mass],[1,1]);reference=drops;reference_mass=mass
    watermark=int(nranks,int64);serial_watermark=watermark;serial_work=spray_sdc_stage_3d()
    call advance_spray_sdc_stage_3d(species,properties,liquid,options,full,temperature,zero,spacing,0._dp,interval, &
      0,2,faces,reference,reference_mass,serial_work,serial_rhs,wall,outside,serial_force,ok, &
      last_id=serial_watermark,breakup_properties=materials,native_flow_cfl=1._dp)
    call require(ok.and.serial_work%flow_subcycles==2,'serial two-flow TAB predictor')
    call MPI_Allreduce(serial_rhs,global_rhs,size(serial_rhs),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    allocate(tangent(1));tangent=0;rhs=0;cycle_rhs=rhs
    do flow_cycle=1,2
      call kick_mpi_sdc_transfer(MPI_COMM_WORLD,species,properties,liquid,options,local,local_t,drops,mass,zero, &
        spacing,cells,[nranks,1,1],faces,.5_dp*interval,zero,cycle_rhs,force,ok,sample,sample_t,reynolds,materials)
      call require(ok,'TAB flow-subcycle transfer')
      rhs=rhs+.5_dp*cycle_rhs
      call kick_mpi_sdc_breakup(MPI_COMM_WORLD,species,properties,liquid,options,drops,mass,sample,sample_t,reynolds, &
        zero,spacing,cells,[nranks,1,1],faces,.5_dp*interval,.5_dp*interval,.5_dp*interval,.true.,flow_cycle==2, &
        watermark,rhs,ok,materials,tangent)
      call require(ok.and.size(tangent)==1.and.tangent(1)>0,'TAB tangent survives flow subcycles')
      if(flow_cycle==1) call require(size(drops)==1.and.watermark==nranks,'only last flow subcycle creates TAB children')
      do p=1,size(drops)
        drops(p)%position=drops(p)%position+.5_dp*interval*drops(p)%velocity
      enddo
    enddo
    call require(size(drops)==size(serial_work%midpoint),'cumulative tangent replacement population')
    do p=1,size(drops)
      call require(maxval(abs(drops(p)%velocity-serial_work%midpoint(p)%velocity))<1e-10_dp.and. &
        abs(drops(p)%mass-serial_work%midpoint(p)%mass)<1e-24_dp,'cumulative tangent matches serial child velocities')
    enddo
    call require(maxval(abs(rhs-global_rhs(:,lo:hi,:,:)))<1e-12_dp*max(1._dp,maxval(abs(global_rhs))), &
      'cumulative TAB tangent compensation matches serial flow subcycles')
    if(rank==0) print *, 'PASS: MPI transfer Reynolds, material-specific ETAB/KHRT births and gas compensation'
  end subroutine check_sdc_breakup

  subroutine check_sdc_film()
    type(gas_transport_species), allocatable :: properties(:)
    type(spray_liquid) :: liquid(1)
    type(spray_options) :: options
    type(spray_wall_options) :: settings
    type(spray_wall_film), allocatable :: films(:),expected(:),initial(:),kept(:),births(:)
    type(mpi_sdc_film_stage) :: work,kept_work
    real(dp), allocatable :: mass(:,:),expected_mass(:,:),initial_mass(:,:),kept_mass(:,:),rhs(:,:,:,:),kept_rhs(:,:,:,:)
    real(dp) :: spacing(3)=[.001_dp,.0015_dp,.002_dp],zero(3)=0,wall(5),expected_wall(5),kept_wall(5)
    real(dp) :: serial_rhs(nv,5,3,2),global_rhs(nv,5,3,2),global_wall(5),delta(nv),total_delta(nv),gas_change(nv)
    real(dp), allocatable :: birth_mass(:,:)
    real(dp) :: old_serial_rhs(nv,5,3,2),old_wall(5),old_serial_wall(5)
    real(dp), parameter :: duration=5e-8_dp
    integer :: np,scheme,c
    call load_h2o2_elementary_transport(properties,ok);call require(ok,'film transport')
    call initialize_spray_liquid(species,'H2O',1000._dp,4184._dp,373.15_dp,101325._dp,2.257e6_dp,liquid(1),ok)
    call require(ok,'film liquid')
    faces='periodic';faces(3:4)='isothermal_wall'
    settings%enabled(3)=.true.;settings%evolve=.true.;settings%temperature=400
    settings%liquid_conductivity=.6_dp;settings%minimum_height=1e-15_dp
    ! Empty rank zero and samples crossing each uneven x-slab seam.
    np=merge(0,1,rank==0)
    allocate(initial(np),initial_mass(1,np),rhs(nv,extent(1),extent(2),extent(3)))
    if(np==1) then
      initial(1)=spray_wall_film([(real(first(1)+extent(1),dp)-.25_dp)*spacing(1),0._dp,spacing(3)], &
        [0._dp,1._dp,0._dp],1e-10_dp,330._dp,1e-4_dp,1e-5_dp,1)
      initial_mass(1,1)=initial(1)%mass
    endif
    do scheme=1,2
      options%mixture_time_scheme='conservative'
      if(scheme==2) options%mixture_time_scheme='euler'
      films=initial;mass=initial_mass;expected=initial;expected_mass=initial_mass
      serial_rhs=71;expected_wall=72
      call kick_serial_sdc_films(species,properties,liquid,options,settings,full,temperature,zero,spacing, &
        duration,faces,expected,expected_mass,serial_rhs,expected_wall,ok)
      call require(ok,'serial immutable film oracle')
      call MPI_Allreduce(serial_rhs,global_rhs,size(serial_rhs),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(expected_wall,global_wall,5,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      rhs=73;wall=74
      call kick_mpi_sdc_films(MPI_COMM_WORLD,species,properties,liquid,options,settings,local,local_t,zero,spacing, &
        cells,[nranks,1,1],duration,faces,films,mass,rhs,wall,ok)
      call require(ok,'distributed immutable-gas film kick')
      call require(size(films)==size(expected).and.all(shape(mass)==shape(expected_mass)),'film population oracle')
      if(np==1) then
        call require(films(1)%temperature>initial(1)%temperature.and.films(1)%mass<initial(1)%mass, &
          'film heats and evaporates')
        call require(abs(films(1)%temperature-expected(1)%temperature)<1e-10_dp.and. &
          abs(films(1)%mass-expected(1)%mass)<1e-24_dp.and. &
          abs(films(1)%height-expected(1)%height)<1e-17_dp.and. &
          maxval(abs(mass-expected_mass))<1e-24_dp,'film thermodynamics matches serial primitive sampling')
      endif
      call require(maxval(abs(rhs-global_rhs(:,lo:hi,:,:)))<1e-12_dp*max(1._dp,maxval(abs(global_rhs))), &
        'containing-cell film source matches serial oracle')
      call require(maxval(abs(wall-global_wall))<1e-12_dp*maxval(abs(global_wall)),'global film wall ledger')
      delta=0
      if(np==1) then
        delta(:5)=wall_film_totals(films(1),liquid(1))-wall_film_totals(initial(1),liquid(1))
        delta(5+liquid(1)%vapor_index)=films(1)%mass-initial(1)%mass
      endif
      do c=1,nv
        gas_change(c)=sum(rhs(c,:,:,:))*product(spacing)*duration
      enddo
      delta=delta+gas_change
      call MPI_Allreduce(delta,total_delta,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      total_delta(:5)=total_delta(:5)-wall
      call require(maxval(abs(total_delta))<1e-18_dp,'film gas species and wall energy conservation')
      kept=films;kept_mass=mass;kept_rhs=rhs;kept_wall=wall
      if(rank==0) settings%temperature(3)=401
      call kick_mpi_sdc_films(MPI_COMM_WORLD,species,properties,liquid,options,settings,local,local_t,zero,spacing, &
        cells,[nranks,1,1],duration,faces,films,mass,rhs,wall,ok)
      call require(.not.ok.and.all(rhs==kept_rhs).and.all(wall==kept_wall).and.all(mass==kept_mass), &
        'rank-disagreeing film settings rollback')
      settings%temperature(3)=400
      if(rank==0) liquid(1)%density=1001
      call kick_mpi_sdc_films(MPI_COMM_WORLD,species,properties,liquid,options,settings,local,local_t,zero,spacing, &
        cells,[nranks,1,1],duration,faces,films,mass,rhs,wall,ok)
      call require(.not.ok.and.all(rhs==kept_rhs).and.all(wall==kept_wall),'rank-disagreeing liquid table rollback')
      liquid(1)%density=1000
      ! All controls agree, but only the last rank owns a material whose
      ! valid conductivity table cannot evaluate this film's temperature.
      settings%port_conductivity_models(2)='polynomial'
      settings%port_conductivity_coefficients(:,2)=[.6_dp,0._dp,0._dp,0._dp]
      settings%port_conductivity_temperature_ranges(:,2)=[340._dp,350._dp]
      if(rank==nranks-1) films(1)%material_id=2
      kept=films
      call kick_mpi_sdc_films(MPI_COMM_WORLD,species,properties,liquid,options,settings,local,local_t,zero,spacing, &
        cells,[nranks,1,1],duration,faces,films,mass,rhs,wall,ok)
      call require(.not.ok.and.all(rhs==kept_rhs).and.all(wall==kept_wall).and.all(mass==kept_mass), &
        'one-rank adaptive film failure rejects all private candidates')
      if(np==1) call require(films(1)%mass==kept(1)%mass.and.films(1)%temperature==kept(1)%temperature.and. &
        films(1)%height==kept(1)%height.and.all(films(1)%position==kept(1)%position),'film rollback retains payload')
      settings%port_conductivity_models(2)='';settings%port_conductivity_coefficients(:,2)=0
      settings%port_conductivity_temperature_ranges(:,2)=0
      call require(all(local==full(:,lo:hi,:,:)).and.all(local_t==temperature(lo:hi,:,:)), &
        'film kick never changes gas stage')
    enddo
    ! Full old/cached/final lifecycle, including a film born after the old
    ! kick on the previously empty rank. Births receive only the final kick.
    options%mixture_time_scheme='conservative'
    films=initial;mass=initial_mass;expected=initial;expected_mass=initial_mass
    call kick_serial_sdc_films(species,properties,liquid,options,settings,full,temperature,zero,spacing, &
      duration,faces,expected,expected_mass,old_serial_rhs,old_serial_wall,ok)
    call require(ok,'serial film predictor oracle')
    call MPI_Allreduce(old_serial_rhs,global_rhs,size(old_serial_rhs),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(old_serial_wall,old_wall,5,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call advance_mpi_sdc_film_stage(MPI_COMM_WORLD,species,properties,liquid,options,settings,local,local_t, &
      zero,spacing,cells,[nranks,1,1],0._dp,2*duration,0,2,faces,films,mass,work,rhs,wall,ok)
    call require(ok.and.work%next_stage==1,'film predictor lifecycle')
    call require(all(mass==initial_mass),'predictor retains physical component inventory')
    if(np==1) call require(films(1)%temperature==initial(1)%temperature.and. &
      films(1)%mass==initial(1)%mass,'predictor film payload stays private')
    kept_rhs=rhs;kept_wall=wall
    allocate(births(merge(1,0,rank==0)),birth_mass(1,merge(1,0,rank==0)))
    if(rank==0) then
      births(1)=spray_wall_film([.0009_dp,0._dp,.002_dp],[0._dp,1._dp,0._dp], &
        2e-10_dp,335._dp,1e-4_dp,1e-5_dp,1)
      birth_mass(1,1)=births(1)%mass
    endif
    call append_mpi_sdc_film_births(MPI_COMM_WORLD,births,birth_mass,settings%maximum_films,work,ok)
    call require(ok.and.work%births_attached,'predictor film births attach once')
    expected=[expected,births]
    expected_mass=reshape([expected_mass,birth_mass],[1,size(expected)])
    kept_work=work
    call append_mpi_sdc_film_births(MPI_COMM_WORLD,births,birth_mass,settings%maximum_films,work,ok)
    call require(.not.ok.and.size(work%midpoint)==size(kept_work%midpoint).and. &
      all(work%midpoint_mass==kept_work%midpoint_mass),'duplicate film births reject privately')
    settings%temperature(3)=401
    call advance_mpi_sdc_film_stage(MPI_COMM_WORLD,species,properties,liquid,options,settings,local,local_t, &
      zero,spacing,cells,[nranks,1,1],0._dp,2*duration,1,2,faces,films,mass,work,rhs,wall,ok)
    call require(.not.ok.and.work%next_stage==1.and.all(rhs==kept_rhs).and.all(wall==kept_wall), &
      'collectively changed mid-interval wall settings reject cached stage')
    settings%temperature(3)=400
    call advance_mpi_sdc_film_stage(MPI_COMM_WORLD,species,properties,liquid,options,settings,local,local_t, &
      zero,spacing,cells,[nranks,1,1],0._dp,2*duration,1,2,faces,films,mass,work,rhs,wall,ok)
    call require(ok.and.work%next_stage==2.and.all(rhs==kept_rhs).and.all(wall==kept_wall), &
      'intermediate film stage reuses old source exactly')
    kept_work=work
    if(rank==nranks-1) work%midpoint(1)%temperature=0
    call advance_mpi_sdc_film_stage(MPI_COMM_WORLD,species,properties,liquid,options,settings,local,local_t, &
      zero,spacing,cells,[nranks,1,1],0._dp,2*duration,2,2,faces,films,mass,work,rhs,wall,ok)
    call require(.not.ok.and.work%next_stage==2.and.all(rhs==kept_rhs).and.all(wall==kept_wall).and. &
      all(mass==initial_mass),'one-rank final film rejection retains stage and physical inventory')
    work=kept_work
    call kick_serial_sdc_films(species,properties,liquid,options,settings,full,temperature,zero,spacing, &
      duration,faces,expected,expected_mass,serial_rhs,expected_wall,ok)
    call require(ok,'serial final film and birth oracle')
    call MPI_Allreduce(serial_rhs,global_rhs,size(serial_rhs),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(expected_wall,global_wall,5,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call advance_mpi_sdc_film_stage(MPI_COMM_WORLD,species,properties,liquid,options,settings,local,local_t, &
      zero,spacing,cells,[nranks,1,1],0._dp,2*duration,2,2,faces,films,mass,work,rhs,wall,ok)
    call require(ok.and.work%next_stage==0.and.size(films)==size(expected),'final film stage commits and closes')
    do c=1,size(films)
      call require(abs(films(c)%temperature-expected(c)%temperature)<1e-10_dp.and. &
        abs(films(c)%mass-expected(c)%mass)<1e-24_dp,'committed midpoint and born films match serial oracle')
    enddo
    call require(maxval(abs(rhs-global_rhs(:,lo:hi,:,:)))<1e-12_dp*max(1._dp,maxval(abs(global_rhs))), &
      'final film source excludes predictor births from old half')
    call require(maxval(abs(wall-global_wall-old_wall))<1e-18_dp,'final wall ledger includes both film halves')
    if(rank==0) print *, 'PASS: MPI immutable film transfer, empty rank, wall accounting and collective rollback'
  end subroutine check_sdc_film

  subroutine check_sdc_transfer()
    type(gas_transport_species), allocatable :: properties(:)
    type(spray_liquid) :: liquid(1)
    type(spray_options) :: options
    type(spray_parcel), allocatable :: drops(:),kept(:)
    type(spray_parcel) :: expected
    real(dp), allocatable :: mass(:,:),kept_mass(:,:),sample(:,:),sample_t(:),rhs(:,:,:,:),kept_rhs(:,:,:,:)
    real(dp) :: position(3,1),spacing(3)=[.001_dp,.0015_dp,.002_dp],zero(3)=0
    real(dp) :: oracle(nv),cell(nv),tg,tc,speed,weights(8),expected_mass(1),change(nv),force(5),expected_force(5)
    real(dp) :: local_exchange(nv),global_exchange(nv),local_expected(nv),global_expected(nv)
    real(dp) :: expected_rhs(nv,5,3,2),global_rhs(nv,5,3,2),kept_force(5),global_force(5)
    real(dp) :: deposit_positions(3,3),deposit_changes(nv,3)
    integer :: stencil(3,8),c,scheme,index(3),p,a
    call load_h2o2_elementary_transport(properties,ok);call require(ok,'transfer transport')
    call initialize_spray_liquid(species,'H2O',1000._dp,4184._dp,373.15_dp,101325._dp,2.257e6_dp,liquid(1),ok)
    call require(ok,'transfer liquid')
    faces='periodic'
    position(:,1)=[(real(first(1)+extent(1),dp)-.25_dp)*spacing(1),1.5_dp*spacing(2),spacing(3)]
    call sample_mpi_sdc_primitive(MPI_COMM_WORLD,species,local,local_t,position,zero,spacing,cells, &
      [nranks,1,1],faces,sample,sample_t,ok)
    call require(ok,'distributed primitive sampling')
    call particle_boundary_stencil(position(:,1),zero,spacing,cells,faces,stencil,weights,ok)
    call require(ok,'serial sample stencil')
    oracle=0;tg=0
    do c=1,8
      call reactive_conserved_to_primitive(species,full(:,stencil(1,c),stencil(2,c),stencil(3,c)), &
        temperature(stencil(1,c),stencil(2,c),stencil(3,c)),cell,tc,speed,ok)
      call require(ok,'serial cell EOS')
      oracle=oracle+weights(c)*cell;tg=tg+weights(c)*tc
    enddo
    oracle(5)=mixture_pressure(species,oracle(6:),oracle(1),tg,ok)
    call require(ok.and.maxval(abs(sample(:,1)-oracle))<2e-12_dp*maxval(abs(oracle)), &
      'interpolate primitive cells before EOS pressure recovery')
    call require(abs(sample_t(1)-tg)<1e-10_dp,'weighted cell temperature retains no velocity variance')
    allocate(drops(1),mass(1,1),rhs(nv,extent(1),extent(2),extent(3)))
    ! Signed containing-cell packets cross rank seams and periodic endpoints;
    ! every other cell must stay exactly zero (CIC roundoff is unacceptable).
    deposit_positions(:,1)=position(:,1)+[spacing(1),0._dp,0._dp]
    deposit_positions(:,2)=position(:,1)-[real(cells(1),dp)*spacing(1),0._dp,0._dp]
    deposit_positions(:,3)=real(cells,dp)*spacing
    do p=1,3
      do c=1,nv
        deposit_changes(c,p)=real((rank+1)*p*c,dp)*(-1._dp)**c*product(spacing)
      enddo
    enddo
    rhs=0;expected_rhs=0
    do p=1,3
      do a=1,3
        index(a)=floor(modulo(deposit_positions(a,p),real(cells(a),dp)*spacing(a))/spacing(a))+1
      enddo
      expected_rhs(:,index(1),index(2),index(3))=expected_rhs(:,index(1),index(2),index(3))+ &
        deposit_changes(:,p)/product(spacing)
    enddo
    call MPI_Allreduce(expected_rhs,global_rhs,size(expected_rhs),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call deposit_mpi_sdc_cell_source(MPI_COMM_WORLD,rhs,deposit_positions,deposit_changes,zero,spacing,cells, &
      [nranks,1,1],faces,ok)
    call require(ok.and.all(rhs==global_rhs(:,lo:hi,:,:)),'integer containing-cell routing across ranks and wrap')
    kept_rhs=rhs
    call deposit_mpi_spray_field(MPI_COMM_WORLD,rhs,deposit_positions,deposit_changes,zero,spacing,cells, &
      [nranks,1,1],faces,ok,containing_cell=rank/=0)
    call require(.not.ok.and.all(rhs==kept_rhs),'rank-disagreeing deposition mode rejects without communication mismatch')
    do scheme=1,2
      options%mixture_time_scheme='conservative'
      if(scheme==2) options%mixture_time_scheme='euler'
      drops(1)=spray_parcel(int(rank+1,int64),position(:,1),[20._dp,-2._dp,1._dp],1e-12_dp,350._dp,100._dp)
      mass=drops(1)%mass;kept=drops;kept_mass=mass;expected=drops(1);expected_mass=mass(:,1)
      if(scheme==2) then
        call advance_mixture_frozen_gas(species,properties,oracle,tg,expected,liquid,expected_mass,options, &
          1e-8_dp,change,ok,time_modifier=.5_dp,force_change=expected_force)
      else
        call advance_mixture_frozen_gas(species,properties,oracle,tg,expected,liquid,expected_mass,options, &
          5e-9_dp,change,ok,force_change=expected_force)
      endif
      call require(ok,'serial transfer oracle')
      rhs=51;force=52
      call kick_mpi_sdc_transfer(MPI_COMM_WORLD,species,properties,liquid,options,local,local_t,drops,mass, &
        zero,spacing,cells,[nranks,1,1],faces,1e-8_dp,zero,rhs,force,ok)
      call require(ok,'distributed immutable-gas half transfer')
      call require(all(drops(1)%position==kept(1)%position),'half transfer never drifts')
      call require(maxval(abs(drops(1)%velocity-expected%velocity))<1e-12_dp.and. &
        abs(drops(1)%temperature-expected%temperature)<1e-10_dp.and. &
        abs(mass(1,1)-expected_mass(1))<1e-24_dp,'distributed transfer matches serial sample oracle')
      do c=1,nv
        local_exchange(c)=sum(rhs(c,:,:,:))*product(spacing)*5e-9_dp
      enddo
      local_expected=change
      call MPI_Allreduce(local_exchange,global_exchange,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(local_expected,global_expected,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call require(maxval(abs(global_exchange-global_expected))<1e-12_dp*maxval(abs(global_expected)), &
        'distributed half-source integrates to physical particle exchange')
      call MPI_Allreduce(expected_force,global_force,5,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call require(all(force==global_force),'physical force ledger matches serial transfer')
      kept=drops;kept_rhs=rhs;kept_force=force
      if(rank==0) mass(1,1)=-1
      kept_mass=mass
      call kick_mpi_sdc_transfer(MPI_COMM_WORLD,species,properties,liquid,options,local,local_t,drops,mass, &
        zero,spacing,cells,[nranks,1,1],faces,1e-8_dp,zero,rhs,force,ok)
      call require(.not.ok.and.all(rhs==kept_rhs).and.all(force==kept_force).and.all(mass==kept_mass).and. &
        all(drops(1)%velocity==kept(1)%velocity).and.all(drops(1)%position==kept(1)%position).and. &
        drops(1)%mass==kept(1)%mass.and.drops(1)%temperature==kept(1)%temperature, &
        'one-rank invalid inventory rejects transfer collectively without publishing')
      call require(all(local==full(:,lo:hi,:,:)).and.all(local_t==temperature(lo:hi,:,:)), &
        'half transfer leaves gas immutable')
    enddo
    if(rank==0) print *, 'PASS: MPI primitive interpolation and immutable-gas half transfer'
  end subroutine

  subroutine check_cartesian_les_sdc()
    type(gas_transport_species), allocatable :: properties(:)
    type(elementary_reaction), allocatable :: reactions(:)
    real(dp) :: reference_state(nv,5,3,2),reference_temp(5,3,2),ref_heat(5,3,2)
    real(dp) :: serial_boundary(nv),serial_source(nv),boundary(nv),source(nv)
    real(dp) :: spacing(3)=[.001_dp,.0015_dp,.002_dp],force(3)=[10._dp,-20._dp,30._dp],wall_t(6)=1100
    real(dp) :: diffusion,isotropic,cs,old_boundary(nv),old_source(nv)
    real(dp), allocatable :: state(:,:,:,:),temp(:,:,:),ir(:,:,:,:),heat(:,:,:)
    logical :: flags(4)
    integer :: mode,closure_index
    character(len=24), parameter :: closures(4)=[character(len=24)::'smagorinsky','wale','vreman','dynamic']
    call load_h2o2_elementary_transport(properties,ok);call require(ok,'Cartesian LES SDC transport')
    allocate(reactions(0),ir(nv-4,extent(1),extent(2),extent(3)),heat(extent(1),extent(2),extent(3)))
    do closure_index=1,4
    do mode=1,2
      faces='periodic';flags=.false.;isotropic=.1_dp;cs=.2_dp
      if(closure_index==4) isotropic=0
      if(mode==2) then
        faces(1:2)='isothermal_wall';flags(:2)=.true.
      endif
      reference_state=full;reference_temp=temperature;ref_heat=53
      call advance_les_sdc_3d(species,reference_state,reference_temp,spacing,1e-9_dp,0._dp,2,cs,.7_dp,.7_dp, &
        serial_boundary,serial_source,diffusion,ok,model=trim(closures(closure_index)),ci=isotropic, &
        external_force=force,transport=properties,transport_flags=flags,boundary_faces=faces, &
        wall_temperature=wall_t,heat_release=ref_heat,test_filter='gaussian',test_ratio=4,grid_ratio=2, &
        harmonic_mean=closure_index==4,difmag=.1_dp)
      call require(ok,'serial LES gas SDC oracle')
      state=local;temp=local_t;ir=0;heat=53;boundary=41;source=42
      call advance_mpi_cartesian_gas_sdc(MPI_COMM_WORLD,species,reactions,properties,state,temp,cells,spacing, &
        1e-9_dp,0._dp,2,force,flags,1e-6_dp,1e-12_dp,.1_dp,faces,ir,heat,boundary,source,ok, &
        wall_temperature=wall_t,les_parameters=[cs,.7_dp,.7_dp],les_model=trim(closures(closure_index)), &
        les_ci=isotropic,les_test_filter='gaussian',les_test_ratio=4,les_grid_ratio=2,les_harmonic_mean=closure_index==4)
      call require(ok,'distributed LES gas SDC accepted')
      call require(maxval(abs(state-reference_state(:,lo:hi,:,:)))<2e-12_dp*maxval(abs(reference_state)).and. &
        maxval(abs(temp-reference_temp(lo:hi,:,:)))<2e-8_dp,'owned LES gas and temperature serial oracle')
      call require(all(ir==0).and.all(heat==0),'nonreacting LES rates and heat reset')
      call require(maxval(abs(boundary-serial_boundary))<max(1e-20_dp,2e-10_dp*maxval(abs(serial_boundary))).and. &
        maxval(abs(source-serial_source))<max(1e-20_dp,2e-10_dp*maxval(abs(serial_source))), &
        'LES wall flux is boundary inventory and not external forcing')
      old_boundary=boundary;old_source=source
      cs=.2_dp
      if(rank==0) cs=.3_dp
      call advance_mpi_cartesian_gas_sdc(MPI_COMM_WORLD,species,reactions,properties,state,temp,cells,spacing, &
        1e-9_dp,0._dp,2,force,flags,1e-6_dp,1e-12_dp,.1_dp,faces,ir,heat,boundary,source,ok, &
        wall_temperature=wall_t,les_parameters=[cs,.7_dp,.7_dp],les_model=trim(closures(closure_index)), &
        les_ci=isotropic,les_test_filter='gaussian',les_test_ratio=4,les_grid_ratio=2,les_harmonic_mean=closure_index==4)
      call require(.not.ok.and.all(boundary==old_boundary).and.all(source==old_source), &
        'rank-disagreeing LES controls preserve accepted gas/ledgers')
    enddo
    enddo
    if(rank==0) print *, 'PASS: MPI SDC four LES closures, wide halos, physical ledgers and control rollback'
  end subroutine check_cartesian_les_sdc

  subroutine check_cartesian_gas_sdc()
    type(gas_transport_species), allocatable :: properties(:)
    type(elementary_reaction), allocatable :: reactions(:)
    real(dp) :: reference_state(nv,5,3,2),reference_temp(5,3,2),ref_heat(5,3,2)
    real(dp) :: serial_boundary(nv),serial_source(nv),boundary(nv),source(nv)
    real(dp) :: spacing(3)=[.001_dp,.0015_dp,.002_dp],force(3)=[10._dp,-20._dp,30._dp],wall_t(6)=1100
    real(dp), allocatable :: state(:,:,:,:),temp(:,:,:),ir(:,:,:,:),heat(:,:,:)
    type(rotating_frame_3d), allocatable :: frame
    type(rotating_frame_3d) :: reference_frame
    real(dp) :: lower(3),kept_lower(3),error
    logical :: flags(4),mol,rotating
    integer :: mode,mismatch
    rotating=test_mode=='gas-mol-rotation';mol=test_mode=='gas-mol'.or.rotating
    lower=0
    if(rotating) then
      allocate(frame);frame%axis=3;frame%angular_speed=300;frame%axis_point=[.003_dp,.002_dp,.001_dp]
      lower=[.02_dp,-.03_dp,.01_dp]
      call convert_regular_energy_3d(frame,lower,spacing,full,0,.true.,ok)
      call require(ok,'rotating Cartesian initial energy');local=full(:,lo:hi,:,:)
    endif
    call load_h2o2_elementary_transport(properties,ok);call require(ok,'Cartesian SDC transport')
    allocate(reactions(0),ir(nv-4,extent(1),extent(2),extent(3)),heat(extent(1),extent(2),extent(3)))
    do mode=1,merge(3,2,rotating)
      faces='periodic';flags=.false.
      if(mode==2) then
        faces(1:2)='isothermal_wall';flags(:2)=.true.
      endif
      reference_state=full;reference_temp=temperature;ref_heat=53
      if(mode==3) then
        ir=0
        block
          real(dp) :: reference_ir(nv-4,5,3,2)
          reference_ir=0
          call advance_mpi_cartesian_gas_sdc(MPI_COMM_SELF,species,reactions,properties,reference_state, &
            reference_temp,cells,spacing,1e-9_dp,0._dp,2,force,flags,1e-6_dp,1e-12_dp,.1_dp,faces, &
            reference_ir,ref_heat,serial_boundary,serial_source,ok,mol_mode=.true.,domain_lower=lower,frame=frame, &
            les_parameters=[.17_dp,.7_dp,.7_dp],les_model='dynamic',les_test_filter='gaussian', &
            les_test_ratio=4,les_grid_ratio=2)
        end block
      else if(mol) then
        call advance_reactive_mol_3d(species,reference_state,reference_temp,spacing,1e-9_dp,0._dp,2,.false., &
          1e-6_dp,1e-12_dp,ok,external_force=force,transport=properties,transport_flags=flags,boundary_faces=faces, &
          wall_temperature=wall_t,boundary_change=serial_boundary,source_change=serial_source,heat_release=ref_heat, &
          lax_friedrichs=mode==2,domain_lower=lower,frame=frame)
      else
      call advance_reactive_sdc_3d(species,reference_state,reference_temp,spacing,1e-9_dp,0._dp,2,2, &
        ok=ok,external_force=force,transport=properties,transport_flags=flags,boundary_faces=faces, &
        wall_temperature=wall_t,boundary_change=serial_boundary,source_change=serial_source,heat_release=ref_heat)
      endif
      call require(ok,'serial gas SDC reference')
      state=local;temp=local_t;ir=0;heat=53;boundary=41;source=42
      if(mode==3) then
        call advance_mpi_cartesian_gas_sdc(MPI_COMM_WORLD,species,reactions,properties,state,temp,cells,spacing, &
          1e-9_dp,0._dp,2,force,flags,1e-6_dp,1e-12_dp,.1_dp,faces,ir,heat,boundary,source,ok, &
          mol_mode=.true.,domain_lower=lower,frame=frame,les_parameters=[.17_dp,.7_dp,.7_dp],les_model='dynamic', &
          les_test_filter='gaussian',les_test_ratio=4,les_grid_ratio=2)
      else
      call advance_mpi_cartesian_gas_sdc(MPI_COMM_WORLD,species,reactions,properties,state,temp,cells,spacing, &
        1e-9_dp,0._dp,2,force,flags,1e-6_dp,1e-12_dp,.1_dp,faces,ir,heat,boundary,source,ok, &
        wall_temperature=wall_t,mol_mode=mol,mol_lax_friedrichs=mol.and.mode==2,domain_lower=lower,frame=frame)
      endif
      call require(ok,'distributed Cartesian gas SDC')
      if(rotating) then
        error=maxval(abs(state-reference_state(:,lo:hi,:,:))/max(1e-14_dp,abs(reference_state(:,lo:hi,:,:))))
        call require(error<2e-11_dp,'rotating Cartesian owned component agreement')
        if(rank==0) print *, 'Rotating MPI MOL mode, component error:',mode,error
      endif
      call require(maxval(abs(state-reference_state(:,lo:hi,:,:)))/maxval(abs(reference_state))<2e-12_dp, &
        'owned gas SDC matches serial global cells')
      call require(maxval(abs(temp-reference_temp(lo:hi,:,:)))<2e-8_dp,'owned gas SDC temperature matches serial')
      call require(all(heat==0).and.all(ir==0),'nonreacting MPI SDC zero rates and heat')
      call require(maxval(abs(boundary-serial_boundary))<2e-12_dp*max(1._dp,maxval(abs(serial_boundary))), &
        'MPI global physical boundary ledger matches serial')
      call require(maxval(abs(source-serial_source))<2e-12_dp*max(1._dp,maxval(abs(serial_source))), &
        'MPI global volumetric source ledger matches serial')
    enddo
    if(rotating) then
      reference_frame=frame;kept_lower=lower
      do mismatch=1,5
        state=local;temp=local_t;ir=0;heat=53;boundary=41;source=42
        if(rank==0) then
          select case(mismatch)
          case(1);frame%axis=2
          case(2);frame%angular_speed=301
          case(3);frame%axis_point(1)=frame%axis_point(1)+.001_dp
          case(4);lower(1)=lower(1)+.001_dp
          case(5);deallocate(frame)
          end select
        endif
        call advance_mpi_cartesian_gas_sdc(MPI_COMM_WORLD,species,reactions,properties,state,temp,cells,spacing, &
          1e-9_dp,0._dp,2,force,flags,1e-6_dp,1e-12_dp,.1_dp,faces,ir,heat,boundary,source,ok, &
          mol_mode=.true.,domain_lower=lower,frame=frame)
        call require(.not.ok.and.all(state==local).and.all(temp==local_t).and.all(ir==0).and.all(heat==53), &
          'rank-disagreeing frame rejects without publication')
        call require(all(boundary==41).and.all(source==42),'rank-disagreeing frame retains ledgers')
        frame=reference_frame;lower=kept_lower
      enddo
    endif
    deallocate(reactions);allocate(reactions(1))
    state=local;temp=local_t;ir=0;heat=53;boundary=41;source=42;reactor_calls=0
    call advance_mpi_cartesian_gas_sdc(MPI_COMM_WORLD,species,reactions,properties,state,temp,cells,spacing, &
      1e-9_dp,0._dp,2,force,flags,1e-6_dp,1e-12_dp,.1_dp,faces,ir,heat,boundary,source,ok, &
      wall_temperature=wall_t,evaluate_reaction=rank_reactor,mol_mode=mol,domain_lower=lower,frame=frame)
    call require(.not.ok.and.reactor_calls>0,'distributed gas reactor failure is collective')
    call require(all(state==local).and.all(temp==local_t).and.all(ir==0).and.all(heat==53), &
      'distributed gas failure preserves state, temperature, IR and heat')
    call require(all(boundary==41).and.all(source==42),'distributed gas failure preserves global budgets')
    if(rank==0) print *, 'PASS: Cartesian gas SDC, physical masks and collective rollback'
  end subroutine

  subroutine check_mol_blend_agreement()
    procedure(reactive_mol_blend_agreement), pointer :: blend_proc
    real(dp), allocatable :: gas(:,:,:,:),temp(:,:,:),held(:,:,:),initial_held(:,:,:)
    real(dp) :: source(nv),boundary(nv),prefix,blend
    nullify(blend_proc)
    if(rank==0) blend_proc=>unused_blend
    gas=local;temp=local_t
    held=gas(5,:,:,:)-.5_dp*sum(gas(2:4,:,:,:)**2,dim=1)/gas(1,:,:,:);initial_held=held
    source=41;boundary=43;prefix=47;blend=53
    call advance_reactive_mol_3d(species,gas,temp,[.001_dp,.0015_dp,.002_dp],1e-9_dp,0._dp,1,.false., &
      1e-6_dp,1e-12_dp,ok,source_change=source,boundary_change=boundary,internal_energy_density=held, &
      energy_correction=prefix,minimum_flux_blend=blend,agree_stage=blend_stage_agreement,agree_flux_blend=blend_proc)
    call require(.not.ok,'rank-disagreeing window callback presence refuses before recovery collectives')
    call require(all(gas==local).and.all(temp==local_t).and.all(held==initial_held).and.prefix==47.and.blend==53, &
      'window callback-presence gas/UEINT/prefix rollback')
    call require(all(source==41).and.all(boundary==43),'window callback-presence inventory rollback')
    if(rank==0) print *, 'PASS: collective MOL blend callback presence disagreement and full rollback'
  end subroutine check_mol_blend_agreement
  subroutine blend_stage_agreement(local_ok,global_ok)
    logical, intent(in) :: local_ok
    logical, intent(out) :: global_ok
    integer :: status
    call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,status)
    if(status/=MPI_SUCCESS) global_ok=.false.
  end subroutine blend_stage_agreement
  subroutine unused_blend(local_blend,global_blend,good)
    real(dp), intent(in) :: local_blend
    real(dp), intent(out) :: global_blend
    logical, intent(out) :: good
    global_blend=local_blend;good=.false.
    call require(.false.,'blend collective entered despite callback presence mismatch')
  end subroutine unused_blend
  subroutine check_sdc_agreement()
    type(elementary_reaction) :: reaction(1)
    real(dp) :: reference_state(nv,5,3,2),reference_temp(5,3,2)
    real(dp) :: state(nv,5,3,2),temp(5,3,2),ir(nv-4,5,3,2),heat(5,3,2)
    real(dp) :: zero_force(3)=0,spacing(3)=[.001_dp,.0015_dp,.002_dp]
    integer :: minimum_calls,maximum_calls
    real(dp) :: les_boundary(nv),les_source(nv),les_diffusion,coefficient,auxiliary(nv)
    reference_state=full;reference_temp=temperature
    call advance_reactive_sdc_3d(species,reference_state,reference_temp,spacing,1e-8_dp,0._dp,2,1, &
      ok=ok,external_force=zero_force)
    call require(ok,'serial SDC agreement reference')
    state=full;temp=temperature;agreement_calls=0
    call advance_reactive_sdc_3d(species,state,temp,spacing,1e-8_dp,0._dp,2,1, &
      ok=ok,external_force=zero_force,agree_stage=agree_sdc_stage)
    call require(ok.and.agreement_calls>0,'collective SDC accepted stage path')
    call require(all(state==reference_state).and.all(temp==reference_temp),'agreement retains exact serial SDC result')
    state=full;temp=temperature;ir=0;heat=53;reactor_calls=0;agreement_calls=0
    call advance_reactive_sdc_3d(species,state,temp,spacing,1e-8_dp,0._dp,2,1, &
      ok=ok,external_force=zero_force,reactions=reaction,reaction_source=ir,heat_release=heat, &
      evaluate_reaction=rank_reactor,agree_stage=agree_sdc_stage)
    call require(.not.ok.and.reactor_calls>0,'single-rank reactor rejection is collective')
    call require(all(state==full).and.all(temp==temperature).and.all(ir==0).and.all(heat==53), &
      'collective reactor rejection retains gas, rates and heat')
    call MPI_Allreduce(agreement_calls,minimum_calls,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(agreement_calls,maximum_calls,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    call require(minimum_calls==maximum_calls,'failed chemistry stops peers at same stage')
    state=full;temp=temperature;source_calls=0;agreement_calls=0
    call advance_reactive_sdc_3d(species,state,temp,spacing,1e-8_dp,0._dp,2,1,rank_source,ok, &
      agree_stage=agree_sdc_stage)
    call require(.not.ok.and.source_calls==2,'single-rank new source failure stops before next stage')
    call require(all(state==full).and.all(temp==temperature),'source failure preserves interval input')
    call MPI_Allreduce(agreement_calls,minimum_calls,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(agreement_calls,maximum_calls,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    call require(minimum_calls==maximum_calls,'failed source stops peers at same stage')
    reference_state=full;reference_temp=temperature
    call advance_les_sdc_3d(species,reference_state,reference_temp,spacing,1e-8_dp,0._dp,2,.1_dp,.7_dp,.7_dp, &
      les_boundary,les_source,les_diffusion,ok)
    call require(ok,'serial LES SDC agreement reference')
    state=full;temp=temperature;agreement_calls=0
    call advance_les_sdc_3d(species,state,temp,spacing,1e-8_dp,0._dp,2,.1_dp,.7_dp,.7_dp, &
      les_boundary,les_source,les_diffusion,ok,agree_stage=agree_sdc_stage)
    call require(ok.and.agreement_calls>0,'LES wrapper forwards collective stage agreement')
    call require(all(state==reference_state).and.all(temp==reference_temp),'agreement retains exact LES SDC result')
    state=full;temp=temperature;ir=0;heat=53;reactor_calls=0;agreement_calls=0
    les_boundary=41;les_source=42;les_diffusion=43
    call advance_les_sdc_3d(species,state,temp,spacing,1e-8_dp,0._dp,2,.1_dp,.7_dp,.7_dp, &
      les_boundary,les_source,les_diffusion,ok,reactions=reaction,reaction_source=ir,heat_release=heat, &
      evaluate_reaction=rank_reactor,agree_stage=agree_sdc_stage)
    call require(.not.ok.and.reactor_calls>0,'LES single-rank reactor failure is collective')
    call require(all(state==full).and.all(temp==temperature).and.all(ir==0).and.all(heat==53), &
      'LES reactor rejection retains gas, IR and heat')
    call require(all(les_boundary==41).and.all(les_source==42).and.les_diffusion==43,'LES budget rollback')
    coefficient=.1_dp
    if(rank==0) coefficient=-1
    agreement_calls=0
    call advance_les_sdc_3d(species,state,temp,spacing,1e-8_dp,0._dp,2,coefficient,.7_dp,.7_dp, &
      les_boundary,les_source,les_diffusion,ok,agree_stage=agree_sdc_stage)
    call require(.not.ok.and.agreement_calls==1,'LES local preflight rejection is collective')
    fail_les_halo=.true.;source_calls=0;auxiliary=44;agreement_calls=0
    call advance_les_sdc_3d(species,state,temp,spacing,1e-8_dp,0._dp,2,.1_dp,.7_dp,.7_dp, &
      les_boundary,les_source,les_diffusion,ok,agree_stage=agree_sdc_stage,fill_halos=les_failure_halo, &
      fill_source_halos=les_local_source_halo,auxiliary_source=collective_les_auxiliary,auxiliary_change=auxiliary)
    call require(.not.ok.and.source_calls==0,'LES local kernel failure stops before collective auxiliary callback')
    call require(all(state==full).and.all(temp==temperature).and.all(auxiliary==44),'LES kernel failure rollback')
    call require(all(les_boundary==41).and.all(les_source==42).and.les_diffusion==43,'LES kernel budget rollback')
    call MPI_Allreduce(agreement_calls,minimum_calls,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(agreement_calls,maximum_calls,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    call require(minimum_calls==maximum_calls,'LES failure stops peers at same stage')
    fail_les_halo=.false.;fail_les_auxiliary=.true.;source_calls=0;agreement_calls=0
    call advance_les_sdc_3d(species,state,temp,spacing,1e-8_dp,0._dp,2,.1_dp,.7_dp,.7_dp, &
      les_boundary,les_source,les_diffusion,ok,agree_stage=agree_sdc_stage,fill_halos=les_failure_halo, &
      fill_source_halos=les_local_source_halo,auxiliary_source=collective_les_auxiliary,auxiliary_change=auxiliary)
    call require(.not.ok.and.source_calls==2,'LES rank-local auxiliary rejection stops peers at new source')
    call require(all(state==full).and.all(temp==temperature).and.all(auxiliary==44),'LES auxiliary failure rollback')
    call require(all(les_boundary==41).and.all(les_source==42).and.les_diffusion==43,'LES auxiliary budget rollback')
    fail_les_auxiliary=.false.
    if(rank==0) print *, 'PASS: SDC collective stage agreement and local failure rollback'
  end subroutine

  subroutine observe_native_reactor(records,reactions,energy,forcing,energy_forcing,duration,rtol,atol, &
      densities,temp,steps,good,maximum_steps)
    type(nasa7_species), intent(in) :: records(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: energy,forcing(:),energy_forcing,duration,rtol,atol
    real(dp), intent(inout) :: densities(:),temp
    integer, intent(out) :: steps
    logical, intent(out) :: good
    integer, optional, intent(in) :: maximum_steps
    real(dp) :: started
    integer :: input_unit,limit
    profile_calls=profile_calls+1
    write(profile_index,'(i0)') profile_calls
    profile_path=trim(profile_prefix)//'.rank'//trim(profile_rank)//'.input'//trim(profile_index)//'.dat'
    open(newunit=input_unit,file=trim(profile_path),status='new',action='write',iostat=profile_io)
    call require(profile_io==0,'native reactor exact input capture')
    limit=-1
    if(present(maximum_steps)) limit=maximum_steps
    write(input_unit,'(3(i0,1x))') size(records),size(reactions),limit
    write(input_unit,'(*(es27.18e3,1x))') energy,energy_forcing,duration,rtol,atol,temp
    write(input_unit,'(*(es27.18e3,1x))') densities
    write(input_unit,'(*(es27.18e3,1x))') forcing
    close(input_unit)
    write(profile_unit,'(i0,",",i0,a)') profile_calls,profile_phase,',start,0,0,F,0,0'
    flush(profile_unit)
    started=MPI_Wtime()
    call advance_forced_reactor_interval(records,reactions,energy,forcing,energy_forcing,duration,rtol,atol, &
      densities,temp,steps,good,maximum_steps=maximum_steps)
    write(profile_unit,'(i0,",",i0,a,es24.16,",",i0,",",l1,",",es24.16,",",es24.16)') &
      profile_calls,profile_phase,',end,',MPI_Wtime()-started,steps,good,minval(densities),temp
    flush(profile_unit)
  end subroutine observe_native_reactor

  subroutine agree_sdc_stage(local_ok,all_ok)
    logical, intent(in) :: local_ok
    logical, intent(out) :: all_ok
    integer :: status
    agreement_calls=agreement_calls+1
    call MPI_Allreduce(local_ok,all_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,status)
    if(status/=MPI_SUCCESS) all_ok=.false.
  end subroutine

  subroutine les_failure_halo(state,temp,time,width,viscous,ghost,gt,good, &
      internal_energy_density,internal_energy_ghost)
    real(dp), optional, intent(in) :: internal_energy_density(:,:,:)
    real(dp), optional, intent(out) :: internal_energy_ghost(:,:,:)
    real(dp), intent(in) :: state(:,:,:,:),temp(:,:,:),time
    integer, intent(in) :: width
    logical, intent(in) :: viscous
    real(dp), intent(out) :: ghost(:,:,:,:),gt(:,:,:)
    logical, intent(out) :: good
    character(len=8) :: periodic_faces(6)
    periodic_faces='periodic'
    call fill_reactive_boundary_ghosts_3d(state,temp,periodic_faces,ghost,gt,good, &
      ghost_width=width,viscous_wall=viscous, &
      internal_energy_density=internal_energy_density,internal_energy_ghost=internal_energy_ghost)
    ! Return valid communication data but force only rank zero's local LES
    ! physics kernel to fail, before the auxiliary source enters a collective.
    if(fail_les_halo.and.rank==0.and.width==2) ghost(1,:,:,:)=-1
  end subroutine

  subroutine les_local_source_halo(source,time,ghost,good)
    real(dp), intent(in) :: source(:,:,:,:),time
    real(dp), intent(out) :: ghost(:,:,:,:)
    logical, intent(out) :: good
    integer :: n(3),a,b,c,index(3)
    n=shape(source(1,:,:,:))
    do c=1,n(3)+8
    do b=1,n(2)+8
    do a=1,n(1)+8
      index=modulo([a,b,c]-5,n)+1
      ghost(:,a,b,c)=source(:,index(1),index(2),index(3))
    enddo
    enddo
    enddo
    good=time>=0
  end subroutine

  subroutine collective_les_auxiliary(state,temp,time,interval,stage,rate,good)
    real(dp), intent(in) :: state(:,:,:,:),temp(:,:,:),time,interval
    integer, intent(in) :: stage
    real(dp), intent(out) :: rate(:,:,:,:)
    logical, intent(out) :: good
    integer :: total,status
    source_calls=source_calls+1;rate=0
    call MPI_Allreduce(stage,total,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,status)
    good=status==MPI_SUCCESS.and.total==nranks*stage
    if(fail_les_auxiliary.and.rank==0.and.stage==1) good=.false.
  end subroutine

  subroutine rank_source(state,temp,time,rates,good)
    real(dp), intent(in) :: state(:,:,:,:),temp(:,:,:),time
    real(dp), intent(out) :: rates(:,:,:,:,:)
    logical, intent(out) :: good
    source_calls=source_calls+1;rates=0
    good=rank/=0.or.source_calls<2
  end subroutine

  subroutine rank_reactor(sp,rx,e0,external_y,external_e,duration,rtol,atol,densities,temp,count,good,maximum_steps)
    type(nasa7_species), intent(in) :: sp(:)
    type(elementary_reaction), intent(in) :: rx(:)
    real(dp), intent(in) :: e0,external_y(:),external_e,duration,rtol,atol
    real(dp), intent(inout) :: densities(:),temp
    integer, intent(out) :: count
    logical, intent(out) :: good
    integer, optional, intent(in) :: maximum_steps
    reactor_calls=reactor_calls+1;count=1;good=rank/=0
  end subroutine

  subroutine check_source_halo()
    real(dp) :: global_source(3,5,3,2)
    real(dp), allocatable :: source(:,:,:,:),halo(:,:,:,:),kept(:,:,:,:)
    integer :: x,y,z,m,w,a,donor(3),index(3)
    do z=1,2
    do y=1,3
    do x=1,5
      global_source(:,x,y,z)=[real(x+10*y+100*z,dp),-real(2*x+y+z,dp),0._dp]
    enddo
    enddo
    enddo
    source=global_source(:,first(1)+1:first(1)+extent(1),:,:)
    do m=1,3
      faces='periodic'
      if(m==2) faces='outflow'
      if(m==3) then
        faces(1:2)='inflow';faces(3:4)='slip_wall';faces(5:6)='isothermal_wall'
      endif
      do w=1,7,3
        call fill_mpi_sdc_source_ghosts_3d(MPI_COMM_WORLD,source,cells,faces,w,halo,ok)
        call require(ok,'signed distributed source halo')
        do z=1,extent(3)+2*w
        do y=1,extent(2)+2*w
        do x=1,extent(1)+2*w
          index=[first(1)+x-w,y-w,z-w]
          do a=1,3
            donor(a)=max(1,min(cells(a),index(a)))
            if(faces(2*a-1)=='periodic') donor(a)=modulo(index(a)-1,cells(a))+1
          enddo
          call require(all(halo(:,x,y,z)==global_source(:,donor(1),donor(2),donor(3))), &
            'source halo matches global signed source exactly')
        enddo
        enddo
        enddo
      enddo
    enddo
    kept=halo
    if(rank==0) source(1,1,1,1)=ieee_value(0._dp,ieee_quiet_nan)
    call fill_mpi_sdc_source_ghosts_3d(MPI_COMM_WORLD,source,cells,faces,7,halo,ok)
    call require(.not.ok.and.all(halo==kept),'nonfinite source rejects collectively and preserves halo')
    source=global_source(:,first(1)+1:first(1)+extent(1),:,:)
    w=7;if(rank==0) w=4
    call fill_mpi_sdc_source_ghosts_3d(MPI_COMM_WORLD,source,cells,faces,w,halo,ok)
    call require(.not.ok.and.all(halo==kept),'source width disagreement preserves halo')
    faces='periodic';if(rank==0) faces='outflow'
    call fill_mpi_sdc_source_ghosts_3d(MPI_COMM_WORLD,source,cells,faces,7,halo,ok)
    call require(.not.ok.and.all(halo==kept),'source boundary disagreement preserves halo')
    if(rank==0) print *, 'PASS: signed SDC source halos and collective rollback'
  end subroutine
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    integer :: status
    if(condition) return
    print *,rank,'FAIL: '//message
    call MPI_Abort(MPI_COMM_WORLD,1,status)
  end subroutine
end program
