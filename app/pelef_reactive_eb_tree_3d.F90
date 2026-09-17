program pelef_reactive_eb_tree_3d
#ifdef PELEF_MPI_EB_TREE
  use mpi_f08
  use mpi_eb_tree_initialize_3d_mod, only: initialize_mpi_eb_tree_3d
  use mpi_eb_tree_checkpoint_read_3d_mod, only: read_mpi_eb_tree_checkpoint_3d
  use mpi_eb_tree_checkpoint_3d_mod, only: write_mpi_eb_tree_checkpoint_3d
  use mpi_eb_tree_csv_3d_mod, only: write_mpi_eb_tree_csv_3d
  use amr_patch_tree_3d_mod, only: amr_tree_state_3d
  use mpi_amr_patch_tree_distribution_3d_mod, only: partition_amr_tree_3d,rebalance_mpi_amr_tree_3d
  use mpi_amr_tree_work_mod, only: begin_tree_work,end_tree_work
  use mpi_eb_tree_tagging_3d_mod, only: tag_mpi_eb_tree_temperature_3d
  use mpi_eb_tree_remap_3d_mod, only: remap_mpi_eb_tree_3d
  use mpi_eb_tree_timestep_3d_mod, only: mpi_eb_tree_timestep_3d
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_full_3d
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d
#endif
#ifdef PELEF_MPI_EB_SPRAY
  use spray_parcel_mod
  use spray_mixture_mod, only: liquid_mixture_properties
  use spray_mpi_amr_tree_mod, only: migrate_eb_tree_spray_parcels
  use spray_mpi_eb_reactive_mod, only: advance_mpi_eb_spray_reactive
#endif
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64,iostat_end,iostat_eor
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
#ifndef PELEF_SELECTED_EB_TREE
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use transport_database_mod, only: load_h2o2_elementary_transport
#else
  use @PELEF_BUNDLE_MODULE_NAME@, only: &
    load_selected_thermo => @PELEF_BUNDLE_THERMO_LOADER_NAME@, &
    load_selected_reactions => @PELEF_BUNDLE_LOADER_NAME@, &
    load_selected_transport => @PELEF_BUNDLE_TRANSPORT_LOADER_NAME@, &
    selected_nspecies => @PELEF_BUNDLE_NSPECIES_SYMBOL@, &
    selected_nreactions => @PELEF_BUNDLE_NREACTIONS_SYMBOL@, &
    selected_integrator => @PELEF_BUNDLE_CHEMISTRY_INTEGRATOR_SYMBOL@, &
    element_count => @PELEF_BUNDLE_NELEMENTS_SYMBOL@, &
    selected_element_names => @PELEF_BUNDLE_ELEMENT_NAMES_SYMBOL@, &
    selected_element_atoms => @PELEF_BUNDLE_ELEMENT_ATOMS_SYMBOL@
  use selected_mechanism_runtime_mod, only: prepare_selected_mechanism,validate_selected_temperature_range
  use simulation_config_reactive_eb_3d_mod, only: resolve_reactive_eb_3d_selected_composition
  use reactive_eb_3d_driver_mod, only: reactive_eb_composition_integrals_3d
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions,mixture_specific_gas_constant
#endif
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,initialize_amr_patch_tree_3d,amr_tree_integer_layout_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config,build_configured_eb_geometry_3d
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    advance_reactive_eb_tree_full_3d,integrate_reactive_eb_tree_3d,reactive_eb_tree_timestep_3d, &
    remap_reactive_eb_tree_3d,tag_reactive_eb_tree_temperature_3d
  use eb_reactive_tree_checkpoint_3d_mod, only: eb_tree_checkpoint_3d,reactive_eb_tree_checkpoint_record_3d
  use selected_composition_mod, only: selected_composition_max_species,selected_composition_name_length
  implicit none
  character(len=:), allocatable :: element_names(:),bundle_identity,checkpoint_integrator
  real(dp), allocatable :: element_atoms(:,:),base_mole_fractions(:)
#ifdef PELEF_MPI_EB_SPRAY
  type(spray_parcel), allocatable :: drops(:)
  type(spray_liquid), allocatable :: liquids(:)
  type(spray_options) :: particle_options
  real(dp), allocatable :: nodal_phi(:,:,:),component_mass(:,:)
  character(len=1024) :: spray_input
  integer :: particle_steps
#endif
  integer :: composition_count=0
  character(len=selected_composition_name_length) :: composition_species(selected_composition_max_species)=''
  real(dp) :: composition_mole_fractions(selected_composition_max_species)=0._dp
#ifdef PELEF_SELECTED_EB_TREE
  character(len=24) :: transport_names(selected_nspecies)
  integer :: transport_geometries(selected_nspecies)
  real(dp) :: transport_values(5,selected_nspecies),mass_fractions(selected_nspecies)
  real(dp) :: common_tmin,common_tmax,gas_constant,regular_t,cut_t
#endif
  integer :: root_shape(3)=[4,4,4],child_count=0,parent(32)=1,ratio(32)=2,lower(3,32)=1,upper(3,32)=2
  logical :: temperature_adaptive=.false.
  integer :: tag_depth=2,tag_ratios(8)=2,tag_tile_shape(3)=[2,2,2],tag_buffer=0,tag_maximum_nodes=33
  real(dp) :: tag_thresholds(8)=1000._dp
  integer :: regrid_count=0,regrid_stage=0,regrid_child_count(16)=0
  integer :: regrid_parent(32,16)=1,regrid_ratio(32,16)=2,regrid_lower(3,32,16)=1,regrid_upper(3,32,16)=2
  real(dp) :: regrid_times(16)=0._dp,next_time
  integer :: maximum_steps=100,maximum_step_retries=12,checkpoint_interval_steps=0
  real(dp) :: domain_lower(3)=0._dp,domain_upper(3)=1._dp,plane_position=.31_dp
  logical :: automatic_timestep=.false.,subcycling=.false.
  real(dp) :: hydro_cfl=.5_dp,transport_cfl=.35_dp
  real(dp) :: time_step=1.e-9_dp,final_time=2.e-9_dp,rtol=2.e-7_dp,atol=1.e-12_dp
  logical :: chemistry_enabled=.true.,transport_enabled=.true.,stop_after_checkpoint=.false.
  character(len=16) :: geometry_type='plane',plane_axis='z',extrusion_axis='z'
  character(len=16) :: riemann_solver='rusanov'
#ifdef PELEF_SELECTED_EB_TREE
  character(len=16) :: chemistry_integrator=selected_integrator
#else
  character(len=16) :: chemistry_integrator='explicit'
#endif
  character(len=512) :: level_set_file='',section_file='',output_file='eb_tree.csv',checkpoint_file='',restart_file=''
  integer :: rebalance_interval_steps=0
  real(dp) :: rebalance_minimum_gain=0.1_dp
  namelist /eb_tree/ root_shape,child_count,parent,ratio,lower,upper,domain_lower,domain_upper, &
    geometry_type,plane_axis,plane_position,extrusion_axis,level_set_file,section_file,time_step,final_time, &
    maximum_steps,maximum_step_retries,rtol,atol,chemistry_enabled,transport_enabled,riemann_solver, &
    automatic_timestep,subcycling,hydro_cfl,transport_cfl,chemistry_integrator,output_file,checkpoint_file,restart_file, &
    checkpoint_interval_steps,stop_after_checkpoint,regrid_count,regrid_times,regrid_child_count, &
    regrid_parent,regrid_ratio,regrid_lower,regrid_upper,temperature_adaptive,tag_depth,tag_ratios, &
    tag_tile_shape,tag_buffer,tag_maximum_nodes,tag_thresholds, &
    composition_count,composition_species,composition_mole_fractions,rebalance_interval_steps,rebalance_minimum_gain
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_tree_3d) :: topology
  type(amr_patch_tree_3d), allocatable :: layouts(:)
  type(eb_geometry_3d), allocatable :: geometry_levels(:)
  type(eb_reactive_tree_3d) :: tree
  type(eb_tree_checkpoint_3d) :: record
  type(eb_geometry_3d) :: finest,coarse
  type(reactive_eb_3d_config) :: config,root_config
  integer(int64), allocatable :: origin(:,:),domain(:,:)
  integer, allocatable :: ratios(:),shared_ratios(:)
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:),current(:),boundary(:),wall(:),chemical(:)
  real(dp) :: dt,theta,error,estimated_dt
  integer :: unit,ios,status,node,n(3),nv,retry,depth,event,common
  logical :: ok,stopped=.false.
  character(len=1024) :: input_path,message,budget_path
#ifdef PELEF_MPI_EB_TREE
  type(amr_tree_state_3d) :: mpi_state,mpi_temperature
  type(eb_geometry_3d), allocatable :: mpi_geometry(:)
  integer, allocatable :: mpi_owners(:)
  integer :: mpi_rank,mpi_size,mpi_ierr
  logical :: distributed=.false.
  logical :: sample_work,redistributed,work_ready
  real(dp), allocatable :: measured_cost(:)
  call MPI_Init(mpi_ierr)
  if (mpi_ierr/=MPI_SUCCESS) error stop 'MPI initialization failed'
  call MPI_Comm_rank(MPI_COMM_WORLD,mpi_rank,mpi_ierr)
  if (mpi_ierr/=MPI_SUCCESS) error stop 'MPI rank failed'
  call MPI_Comm_size(MPI_COMM_WORLD,mpi_size,mpi_ierr)
  if (mpi_ierr/=MPI_SUCCESS) error stop 'MPI size failed'
#endif
#ifdef PELEF_MPI_EB_SPRAY
  if(command_argument_count()/=2) error stop 'Usage: MPI EB spray <tree.nml> <spray.nml>'
  call get_command_argument(2,spray_input)
#else
  if (command_argument_count()/=1) error stop 'Usage: pelef_reactive_eb_tree_3d input.nml'
#endif
  call get_command_argument(1,input_path)
  open(newunit=unit,file=trim(input_path),status='old',action='read',iostat=ios)
  if (ios/=0) error stop 'Cannot open input'
  read(unit,nml=eb_tree,iostat=ios)
  close(unit,iostat=status)
  if (ios/=0.or.status/=0) error stop 'Invalid eb_tree input'
#ifdef PELEF_MPI_EB_SPRAY
  if(geometry_type/='level_set') error stop 'MPI EB spray requires original nodal geometry'
  if(checkpoint_interval_steps/=0.or.len_trim(restart_file)>0.or.regrid_count/=0.or.temperature_adaptive) &
    error stop 'MPI EB spray checkpoint and regridding are not connected'
  if(rebalance_interval_steps/=0) error stop 'MPI EB spray rebalancing is not connected'
  if(trim(output_file)==trim(spray_input).or.trim(output_file)//'.budget.csv'==trim(spray_input)) &
    error stop 'Spray input aliases output'
#endif
  if (rebalance_interval_steps<0) error stop 'Invalid rebalance interval'
  if (.not.ieee_is_finite(rebalance_minimum_gain)) error stop 'Invalid rebalance gain'
  if (rebalance_minimum_gain<0._dp.or.rebalance_minimum_gain>=1._dp) error stop 'Invalid rebalance gain'
#ifndef PELEF_MPI_EB_TREE
  if (rebalance_interval_steps/=0) error stop 'Rebalancing requires MPI executable'
#endif
  if (child_count<0.or.child_count>32.or.any(root_shape<2)) error stop 'Invalid tree dimensions'
  if (maximum_steps<1.or.maximum_step_retries<0.or.maximum_step_retries>60) error stop 'Invalid step limits'
  if (.not.all(ieee_is_finite([domain_lower,domain_upper,plane_position,time_step,final_time,rtol,atol,hydro_cfl,transport_cfl]))) &
    error stop 'Nonfinite controls'
  if (any(domain_upper<=domain_lower).or.min(time_step,final_time,rtol,atol)<=0._dp) error stop 'Invalid controls'
  if (hydro_cfl<=0._dp.or.hydro_cfl>1._dp.or.transport_cfl<=0._dp.or.transport_cfl>.5_dp) &
    error stop 'Invalid CFL controls'
  if (chemistry_integrator/='explicit'.and.chemistry_integrator/='implicit') error stop 'Invalid chemistry integrator'
  if (riemann_solver/='rusanov'.and.riemann_solver/='hllc'.and.riemann_solver/='pelec') error stop 'Invalid solver'
  if (len_trim(output_file)==0.or.checkpoint_interval_steps<0) error stop 'Invalid output policy'
  if (checkpoint_interval_steps>0.and.len_trim(checkpoint_file)==0) error stop 'Missing checkpoint path'
  if (stop_after_checkpoint.and.checkpoint_interval_steps==0) error stop 'Missing checkpoint interval'
  budget_path=trim(output_file)//'.budget.csv'
  if (output_file==input_path.or.budget_path==input_path.or.output_file==restart_file.or.budget_path==restart_file.or. &
    output_file==level_set_file.or.budget_path==level_set_file.or.output_file==section_file.or.budget_path==section_file) &
    error stop 'Output aliases input'
  if (len_trim(checkpoint_file)>0) then
    if (checkpoint_file==input_path.or.checkpoint_file==output_file.or.checkpoint_file==budget_path.or. &
      checkpoint_file==restart_file.or.checkpoint_file==level_set_file.or.checkpoint_file==section_file) &
      error stop 'Checkpoint aliases another file'
  end if
  if (regrid_count<0.or.regrid_count>16) error stop 'Invalid regrid count'
  if (temperature_adaptive) then
    if (regrid_count==0.or.tag_depth<1.or.tag_depth>8) error stop 'Invalid adaptation depth or schedule'
    if (any(regrid_child_count/=0).or.any(regrid_parent/=1).or.any(regrid_ratio/=2).or. &
      any(regrid_lower/=1).or.any(regrid_upper/=2)) error stop 'Adaptive mode cannot prescribe event boxes'
    if (any(tag_ratios(:tag_depth)<2).or.any(tag_tile_shape<1).or.tag_buffer<0.or. &
      tag_maximum_nodes<1.or.tag_maximum_nodes>33) error stop 'Invalid tagging controls'
    if (.not.all(ieee_is_finite(tag_thresholds(:tag_depth))).or.any(tag_thresholds(:tag_depth)<=0._dp)) &
      error stop 'Invalid tagging thresholds'
    do depth=2,tag_depth
      if (tag_thresholds(depth)<tag_thresholds(depth-1)) error stop 'Tag thresholds must be nondecreasing'
    end do
  end if
  allocate(layouts(0:regrid_count))
  call make_layout(child_count,parent,ratio,lower,upper,layouts(0))
  do event=1,regrid_count
    if (.not.ieee_is_finite(regrid_times(event))) error stop 'Nonfinite regrid time'
    if (regrid_times(event)<=0._dp) error stop 'Invalid regrid time'
    if (event>1) then
      if (regrid_times(event)<=regrid_times(event-1)) error stop 'Regrid times must increase'
    end if
    call make_layout(regrid_child_count(event),regrid_parent(:,event),regrid_ratio(:,event), &
      regrid_lower(:,:,event),regrid_upper(:,:,event),layouts(event))
  end do
  do event=0,regrid_count
    call amr_tree_integer_layout_3d(layouts(event),origin,domain,ratios,ok)
    if (.not.ok) error stop 'Invalid integer layout'
    if (.not.allocated(shared_ratios)) shared_ratios=ratios
    common=min(size(shared_ratios),size(ratios))
    if (any(shared_ratios(:common)/=ratios(:common))) error stop 'Schedule changes refinement lattice'
    if (size(ratios)>size(shared_ratios)) shared_ratios=ratios
  end do
  if (temperature_adaptive) then
    common=min(size(shared_ratios)-1,tag_depth)
    if (any(shared_ratios(2:common+1)/=tag_ratios(:common))) error stop 'Tagging changes refinement lattice'
    if (tag_depth+1>size(shared_ratios)) shared_ratios=[1,tag_ratios(:tag_depth)]
  end if
  n=root_shape
  do depth=2,size(shared_ratios)
    if (any(n>huge(1)/shared_ratios(depth))) error stop 'Finest lattice too large'
    n=n*shared_ratios(depth)
  end do
  topology=layouts(0)
  config%nx=n(1);config%ny=n(2);config%nz=n(3)
  config%x_lower=domain_lower(1);config%y_lower=domain_lower(2);config%z_lower=domain_lower(3)
  config%x_upper=domain_upper(1);config%y_upper=domain_upper(2);config%z_upper=domain_upper(3)
  config%geometry_type=geometry_type;config%plane_axis=plane_axis;config%plane_position=plane_position
  config%extrusion_axis=extrusion_axis;config%section_file=section_file;config%level_set_file=level_set_file
  config%cfl=hydro_cfl;config%transport_cfl=transport_cfl
  config%redistribution='none';config%riemann_solver=riemann_solver
  config%chemistry_enabled=chemistry_enabled;config%transport_enabled=transport_enabled
  config%chemistry_relative_tolerance=rtol;config%chemistry_absolute_tolerance=atol
  config%viscosity_enabled=.true.;config%thermal_conduction_enabled=.true.
  config%species_diffusion_enabled=.true.;config%barodiffusion_enabled=.true.
  config%final_time=final_time;config%maximum_steps=maximum_steps
#ifdef PELEF_MPI_EB_SPRAY
  call build_configured_eb_geometry_3d(config,finest,ok,nodal_phi)
#else
  call build_configured_eb_geometry_3d(config,finest,ok)
#endif
  if (.not.ok) error stop 'Cannot build finest geometry'
  allocate(geometry_levels(size(shared_ratios)));geometry_levels(size(shared_ratios))=finest
  do depth=size(shared_ratios),2,-1
    call coarsen_eb_geometry_3d(geometry_levels(depth),shared_ratios(depth),geometry_levels(depth-1),ok)
    if (.not.ok) error stop 'Cannot coarsen geometry'
  end do
  coarse=geometry_levels(1)
#ifdef PELEF_SELECTED_EB_TREE
  call load_selected_thermo(species,ok)
  if (.not.ok) error stop 'Selected thermodynamic load failed'
  call load_selected_reactions(reactions,ok)
  if (.not.ok) error stop 'Selected reaction load failed'
  call load_selected_transport(transport_names,transport_geometries,transport_values,ok)
  if (.not.ok) error stop 'Selected transport load failed'
  call prepare_selected_mechanism('EB tree',selected_nspecies,selected_nreactions,species,reactions, &
    transport_names,transport_geometries,transport_values,transport,common_tmin,common_tmax,ok,message)
  if (.not.ok) then
    write(*,'(a)') trim(message)
    error stop 'Invalid selected mechanism'
  end if
  element_names=selected_element_names; element_atoms=selected_element_atoms
  config%thermo_model='selected'
  config%composition_count=composition_count;config%composition_species=composition_species
  config%composition_mole_fractions=composition_mole_fractions
  allocate(base_mole_fractions(size(species)))
  call resolve_reactive_eb_3d_selected_composition(config,species,base_mole_fractions,ok,message)
  if (.not.ok) then
    write(*,'(a)') trim(message)
    error stop 'Invalid selected composition'
  end if
  call mass_fractions_from_mole_fractions(species,base_mole_fractions,mass_fractions,ok)
  if (.not.ok) error stop 'Invalid selected mass fractions'
  gas_constant=mixture_specific_gas_constant(species,mass_fractions,ok)
  if (.not.ok) error stop 'Invalid selected gas constant'
  regular_t=config%initial_pressure/(config%regular_density*gas_constant)
  cut_t=config%initial_pressure/(config%cut_density*gas_constant)
  call validate_selected_temperature_range('EB tree',min(regular_t,cut_t),max(regular_t,cut_t), &
    common_tmin,common_tmax,ok,message)
  if (.not.ok) then
    write(*,'(a)') trim(message)
    error stop 'Selected initial temperature outside mechanism range'
  end if
  bundle_identity='@PELEF_BUNDLE_BUNDLE_SHA256@'
  checkpoint_integrator=chemistry_integrator
#else
  if (composition_count/=0.or.any(composition_species/='').or.any(composition_mole_fractions/=0._dp)) &
    error stop 'Composition inputs require selected executable'
  call load_h2o2_elementary_thermo(species,ok)
  if (.not.ok) error stop 'Thermodynamic load failed'
  call load_h2o2_elementary_mechanism(reactions,ok)
  if (.not.ok) error stop 'Mechanism load failed'
  call load_h2o2_elementary_transport(transport,ok)
  if (.not.ok) error stop 'Transport load failed'
#endif
  nv=reactive_nvar(size(species));n=root_shape
  allocate(current(nv),boundary(nv),wall(nv),chemical(nv))
  root_config=config;root_config%nx=n(1);root_config%ny=n(2);root_config%nz=n(3)
#ifdef PELEF_MPI_EB_TREE
  call initialize_owned()
  call integrate_mpi_eb_tree_3d(topology,mpi_owners,mpi_geometry,mpi_state,current,MPI_COMM_WORLD,ok)
#else
  allocate(u(nv,n(1),n(2),n(3)),t(n(1),n(2),n(3)))
  call initialize_reactive_eb_density_sheet_3d(species,root_config,coarse,u,t,ok,base_mole_fractions)
  if (.not.ok) error stop 'Initial state failed'
  call initialize_reactive_eb_tree_3d(species,topology,geometry_levels(topology%level_count()),u,t,tree,ok)
  if (.not.ok) error stop 'Tree initialization failed'
  call integrate_reactive_eb_tree_3d(tree,current,ok)
#endif
  if (.not.ok) error stop 'Initial inventory failed'
#ifdef PELEF_MPI_EB_SPRAY
  call initialize_particles()
  call add_particle_inventory(current)
#endif
  record%initial=current;record%scale=abs(current);record%boundary=0._dp*current
  record%wall=0._dp*current;record%chemical=0._dp*current
  if (len_trim(restart_file)>0) call restore()
  do while(record%time<final_time.and.record%steps<maximum_steps)
    dt=min(time_step,final_time-record%time)
    if (automatic_timestep) then
#ifdef PELEF_MPI_EB_TREE
      call mpi_eb_tree_timestep_3d(species,transport,topology,mpi_owners,mpi_geometry,mpi_state,mpi_temperature, &
        hydro_cfl,transport_cfl,transport_enabled,.true.,.true.,.true.,estimated_dt,MPI_COMM_WORLD,ok,subcycling)
#else
      call reactive_eb_tree_timestep_3d(species,transport,tree,hydro_cfl,transport_cfl,transport_enabled, &
        .true.,.true.,.true.,estimated_dt,ok,subcycling)
#endif
      if (.not.ok) error stop 'Tree timestep estimate failed'
      dt=min(dt,estimated_dt)
    end if
    if (regrid_stage<regrid_count) dt=min(dt,regrid_times(regrid_stage+1)-record%time)
#ifdef PELEF_MPI_EB_TREE
    sample_work=.false.
    if (rebalance_interval_steps>0) sample_work=mod(record%steps+1,rebalance_interval_steps)==0
#endif
    do retry=0,maximum_step_retries
      if (record%time+dt<=record%time) error stop 'Time increment underflow'
#ifdef PELEF_MPI_EB_TREE
      if (sample_work) then
        call begin_tree_work(size(topology%nodes),ok)
        call MPI_Allreduce(ok,work_ready,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,mpi_ierr)
        if (mpi_ierr/=MPI_SUCCESS.or..not.work_ready) error stop 'Cannot allocate tree work counters'
      end if
#ifdef PELEF_MPI_EB_SPRAY
      call advance_mpi_eb_spray_reactive(MPI_COMM_WORLD,topology,mpi_owners,mpi_geometry,nodal_phi, &
        species,reactions,transport,liquids,particle_options,mpi_state,mpi_temperature,drops,component_mass, &
        domain_lower,(domain_upper-domain_lower)/real(root_shape,dp),dt,riemann_solver, &
        chemistry_enabled,transport_enabled,.true.,.true.,.true.,.true.,rtol,atol,particle_steps, &
        theta,boundary,wall,chemical,ok,chemistry_integrator,subcycling)
#else
      call advance_mpi_eb_tree_full_3d(species,reactions,transport,topology,mpi_owners,mpi_geometry, &
        riemann_solver,dt,chemistry_enabled,rtol,atol,transport_enabled,.true.,.true.,.true.,.true., &
        mpi_state,mpi_temperature,boundary,wall,chemical,theta,MPI_COMM_WORLD,ok,chemistry_integrator,subcycling)
#endif
      if (sample_work) call end_tree_work(measured_cost)
#else
      call advance_reactive_eb_tree_full_3d(species,reactions,transport,riemann_solver,dt,chemistry_enabled, &
        rtol,atol,transport_enabled,.true.,.true.,.true.,.true.,tree,boundary,wall,chemical,theta,ok, &
        chemistry_integrator,subcycling)
#endif
      if (ok) exit
      if (retry==maximum_step_retries) error stop 'Tree step retry limit exhausted'
      dt=.5_dp*dt
    end do
    record%boundary=record%boundary+boundary;record%wall=record%wall+wall;record%chemical=record%chemical+chemical
    next_time=record%time+dt
    if (regrid_stage<regrid_count) then
      if (dt==regrid_times(regrid_stage+1)-record%time) next_time=regrid_times(regrid_stage+1)
    end if
    record%time=next_time;record%steps=record%steps+1
    if (record%steps==1) record%minimum_dt=dt
    record%minimum_dt=min(record%minimum_dt,dt);record%minimum_theta=min(record%minimum_theta,theta)
#ifdef PELEF_MPI_EB_TREE
    if (sample_work) then
      call rebalance_mpi_amr_tree_3d(topology,mpi_owners,nv,mpi_state,mpi_temperature, &
        measured_cost,rebalance_minimum_gain,MPI_COMM_WORLD,redistributed,ok)
      if (.not.ok) error stop 'Cannot rebalance MPI EB tree'
      if (redistributed.and.mpi_rank==0) write(*,'(a,i0)') 'Rebalanced at step: ',record%steps
    end if
#endif
    if (regrid_stage<regrid_count) then
      if (record%time>=regrid_times(regrid_stage+1)) then
        if (temperature_adaptive) then
          depth=max(topology%level_count(),tag_depth+1)
#ifdef PELEF_MPI_EB_TREE
          call tag_mpi_eb_tree_temperature_3d(topology,mpi_owners,mpi_temperature,mpi_geometry,geometry_levels(depth), &
            tag_thresholds(:tag_depth),tag_ratios(:tag_depth),tag_tile_shape,tag_buffer,tag_maximum_nodes, &
            layouts(regrid_stage+1),MPI_COMM_WORLD,ok)
#else
          call tag_reactive_eb_tree_temperature_3d(tree,geometry_levels(depth),tag_thresholds(:tag_depth), &
            tag_ratios(:tag_depth),tag_tile_shape,tag_buffer,tag_maximum_nodes,layouts(regrid_stage+1),ok)
#endif
          if (.not.ok) error stop 'Temperature adaptation failed or exceeded capacity'
        end if
        call activate(regrid_stage+1,.true.)
#ifdef PELEF_MPI_EB_TREE
        if (mpi_rank==0) then
#endif
        write(*,'(a,i0,a,i0)') 'Regrid stage: ',regrid_stage,' nodes: ',size(tree%nodes)
#ifdef PELEF_MPI_EB_TREE
        end if
#endif
      end if
    end if
    call balance()
    if (checkpoint_interval_steps>0) then
      if (mod(record%steps,checkpoint_interval_steps)==0) then
        call save()
        if (stop_after_checkpoint) then
          stopped=.true.
          exit
        end if
      end if
    end if
  end do
  if (record%time<final_time.and..not.stopped) error stop 'Maximum steps reached before final time'
  call balance()
#ifdef PELEF_MPI_EB_SPRAY
  call output_particles()
#endif
  call output()
#ifdef PELEF_MPI_EB_TREE
  if (mpi_rank==0) then
#endif
  write(*,'(a,i0,a,es24.16)') 'Completed steps: ',record%steps,' time: ',record%time
  write(*,'(a,es24.16)') 'Composite inventory error: ',error
  write(*,'(a,es24.16)') 'Minimum accepted interval: ',record%minimum_dt
#ifdef PELEF_MPI_EB_TREE
  end if
  call MPI_Finalize(mpi_ierr)
  if (mpi_ierr/=MPI_SUCCESS) error stop 'MPI finalize failed'
#endif
contains
#ifdef PELEF_MPI_EB_SPRAY
  subroutine initialize_particles()
    integer :: nc,number,q,c,stream,code,closed
    character(len=24) :: vapor_names(16)
    real(dp) :: densities(16),heat_capacities(16),reference_temperatures(16)
    real(dp) :: reference_pressures(16),latent_heats(16),fractions(16)
    real(dp) :: total_mass,center(3),axis(3),speed,half_angle,diameter,liquid_temperature,rho,cp,e
    type(spray_liquid) :: effective
    namelist /eb_spray/ nc,number,vapor_names,densities,heat_capacities,reference_temperatures, &
      reference_pressures,latent_heats,fractions,total_mass,center,axis,speed,half_angle,diameter, &
      liquid_temperature,particle_options
    nc=1;number=1;vapor_names='H2O';densities=1000;heat_capacities=4184;reference_temperatures=373.15_dp
    reference_pressures=101325;latent_heats=2.257e6_dp;fractions=0;fractions(1)=1
    total_mass=1e-10_dp;center=0;axis=[1._dp,0._dp,0._dp];speed=1;half_angle=0
    diameter=1e-5_dp;liquid_temperature=350;particle_options=spray_options()
    open(newunit=stream,file=trim(spray_input),status='old',action='read',iostat=code)
    if(code/=0) error stop 'Cannot read MPI EB spray input'
    read(stream,nml=eb_spray,iostat=code)
    close(stream,iostat=closed)
    if(code/=0.or.closed/=0) error stop 'Invalid MPI EB spray input'
    if(nc<1.or.nc>16.or.number<1.or.number>125000) error stop 'Invalid MPI EB spray dimensions'
    if(.not.all(ieee_is_finite(fractions)).or.any(fractions<0)) error stop 'Invalid liquid fractions'
    if(abs(sum(fractions(:nc))-1)>1e-12_dp.or.any(fractions(nc+1:)/=0)) error stop 'Invalid fraction sum'
    if(.not.valid_spray_options(particle_options).or.particle_options%khrt_children) &
      error stop 'Unsupported particle options'
    if(number>particle_options%maximum_parcels) error stop 'Particle capacity exceeded'
    allocate(liquids(nc))
    do c=1,nc
      call initialize_spray_liquid(species,vapor_names(c),densities(c),heat_capacities(c),reference_temperatures(c), &
        reference_pressures(c),latent_heats(c),liquids(c),ok)
      if(.not.ok) error stop 'Invalid liquid properties'
    end do
    call liquid_mixture_properties(liquids,fractions(:nc),liquid_temperature,rho,cp,e,ok)
    if(.not.ok) error stop 'Invalid liquid mixture'
    effective=liquids(1);effective%density=rho;effective%heat_capacity=cp
    effective%reference_temperature=liquid_temperature;effective%reference_energy=e
    ! Generate deterministically on every rank to validate identical input semantics;
    ! rank zero retains the initial population for the collective owner migration.
    call inject_cone_parcels(effective,total_mass,number,1_int64,center,axis,speed,half_angle, &
      diameter,liquid_temperature,drops,ok)
    if(.not.ok) error stop 'Invalid initial injection'
    if(mpi_rank/=0) then
      deallocate(drops)
      allocate(drops(0))
    end if
    allocate(component_mass(nc,size(drops)))
    do q=1,size(drops)
      component_mass(:,q)=fractions(:nc)*drops(q)%mass
    end do
    call migrate_eb_tree_spray_parcels(MPI_COMM_WORLD,topology,mpi_owners,nodal_phi,drops, &
      domain_lower,(domain_upper-domain_lower)/real(root_shape,dp),particle_options%maximum_parcels,ok,component_mass)
    if(.not.ok) error stop 'Initial EB particle ownership failed'
  end subroutine

  subroutine add_particle_inventory(inventory)
    real(dp), intent(inout) :: inventory(:)
    real(dp) :: local(size(inventory)),global(size(inventory)),rho,cp,e,amount
    integer :: q,c
    local=0
    do q=1,size(drops)
      if(drops(q)%mass<=0) cycle
      call liquid_mixture_properties(liquids,component_mass(:,q)/drops(q)%mass,drops(q)%temperature,rho,cp,e,ok)
      if(.not.ok) error stop 'Invalid particle inventory'
      amount=drops(q)%mass*drops(q)%multiplicity
      local(1)=local(1)+amount;local(2:4)=local(2:4)+amount*drops(q)%velocity
      local(5)=local(5)+amount*(e+.5_dp*sum(drops(q)%velocity**2))
      do c=1,size(liquids)
        local(5+liquids(c)%vapor_index)=local(5+liquids(c)%vapor_index)+component_mass(c,q)*drops(q)%multiplicity
      end do
    end do
    call MPI_Allreduce(local,global,size(local),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpi_ierr)
    if(mpi_ierr/=MPI_SUCCESS) error stop 'Particle inventory reduction failed'
    inventory=inventory+global
  end subroutine

  subroutine output_particles()
    integer :: stream,code,closed,q,c
    character(len=32) :: rank_label,label
    character(len=:), allocatable :: path
    write(rank_label,'(i0)') mpi_rank
    path=trim(output_file)//'.particles.rank'//trim(rank_label)//'.csv'
    open(newunit=stream,file=path,status='new',action='write',iostat=code)
    if(code/=0) error stop 'Cannot create particle CSV (existing files are protected)'
    write(stream,'(a)',advance='no') 'id,x,y,z,u,v,w,mass,temperature,multiplicity'
    do c=1,size(liquids)
      write(label,'(i0)') c
      write(stream,'(a)',advance='no') ',component_mass_'//trim(label)
    end do
    write(stream,*)
    do q=1,size(drops)
      write(stream,'(i0,",",*(es25.16,:,","))',iostat=code) drops(q)%id,drops(q)%position,drops(q)%velocity, &
        drops(q)%mass,drops(q)%temperature,drops(q)%multiplicity,component_mass(:,q)
      if(code/=0) error stop 'Cannot write particle CSV'
    end do
    close(stream,iostat=closed)
    if(closed/=0) error stop 'Cannot close particle CSV'
  end subroutine
#endif

#ifdef PELEF_MPI_EB_TREE
  subroutine initialize_owned()
    integer :: p,ios,ierr
    logical :: valid,all_valid
    real(dp), allocatable :: seed_u(:,:,:,:),seed_t(:,:,:)
    call partition_amr_tree_3d(topology,mpi_size,0,mpi_owners,valid)
    if (.not.valid) error stop 'MPI tree partition failed'
    if (mpi_rank==mpi_owners(1)) then
      allocate(seed_u(nv,n(1),n(2),n(3)),seed_t(n(1),n(2),n(3)),stat=ios)
      valid=ios==0
      if (valid) call initialize_reactive_eb_density_sheet_3d( &
        species,root_config,coarse,seed_u,seed_t,valid,base_mole_fractions)
    else
      ! The collective initializer reads the seed only on the root patch owner.
      allocate(seed_u(0,0,0,0),seed_t(0,0,0),stat=ios)
      valid=ios==0
    end if
    call MPI_Allreduce(valid,all_valid,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
    if (ierr/=MPI_SUCCESS) error stop 'MPI seed agreement failed'
    if (.not.all_valid) error stop 'Initial state failed'
    call initialize_mpi_eb_tree_3d(species,topology,mpi_owners,geometry_levels(topology%level_count()),seed_u,seed_t, &
      mpi_geometry,mpi_state,mpi_temperature,MPI_COMM_WORLD,valid)
    if (.not.valid) error stop 'MPI tree initialization failed'
    tree%topology=topology
    if (allocated(tree%nodes)) deallocate(tree%nodes)
    allocate(tree%nodes(size(topology%nodes)))
    do p=1,size(tree%nodes)
      tree%nodes(p)%geometry=mpi_geometry(p)
    end do
    distributed=.true.
  end subroutine
#endif
  subroutine make_layout(number,parents,refinement,lo,hi,result)
    integer, intent(in) :: number,parents(:),refinement(:),lo(:,:),hi(:,:)
    type(amr_patch_tree_3d), intent(out) :: result
    type(amr_patch_3d), allocatable :: patches(:)
    integer :: node,parent_node,dims(3)
    logical :: valid
    if (number<0.or.number>32) error stop 'Invalid scheduled child count'
    allocate(patches(number))
    do node=1,number
      parent_node=parents(node)
      if (parent_node<1.or.parent_node>node) error stop 'Parents must precede children'
      dims=root_shape
      if (parent_node>1) dims=(hi(:,parent_node-1)-lo(:,parent_node-1)+1)*refinement(parent_node-1)
      call initialize_amr_patch_3d(dims(1),dims(2),dims(3),lo(1,node),hi(1,node),lo(2,node),hi(2,node), &
        lo(3,node),hi(3,node),refinement(node),patches(node),valid)
      if (.not.valid) error stop 'Invalid scheduled box'
    end do
    call initialize_amr_patch_tree_3d(root_shape,parents(:number),patches,result,valid)
    if (.not.valid) error stop 'Invalid scheduled topology'
  end subroutine
  subroutine activate(stage,remapping)
    integer, intent(in) :: stage
    logical, intent(in) :: remapping
    integer :: level
    logical :: valid
#ifdef PELEF_MPI_EB_TREE
    integer, allocatable :: next_owners(:)
    integer :: p
    if (.not.remapping) then
      topology=layouts(stage);regrid_stage=stage
      call initialize_owned()
      return
    end if
    if (distributed.and.remapping) then
      level=max(topology%level_count(),layouts(stage)%level_count())
      call partition_amr_tree_3d(layouts(stage),mpi_size,0,next_owners,valid)
      if (.not.valid) error stop 'MPI regrid partition failed'
      call remap_mpi_eb_tree_3d(species,topology,layouts(stage),mpi_owners,next_owners,geometry_levels(level), &
        mpi_geometry,mpi_state,mpi_temperature,MPI_COMM_WORLD,valid)
      if (.not.valid) error stop 'MPI conservative regrid failed'
      mpi_owners=next_owners;topology=layouts(stage);regrid_stage=stage
      tree%topology=topology
      deallocate(tree%nodes)
      allocate(tree%nodes(size(topology%nodes)))
      do p=1,size(tree%nodes)
        tree%nodes(p)%geometry=mpi_geometry(p)
      end do
      return
    end if
#endif
    if (remapping) then
      level=max(topology%level_count(),layouts(stage)%level_count())
      call remap_reactive_eb_tree_3d(species,layouts(stage),geometry_levels(level),tree,valid)
    else
      level=layouts(stage)%level_count()
      call initialize_reactive_eb_tree_3d(species,layouts(stage),geometry_levels(level),u,t,tree,valid)
    end if
    if (.not.valid) error stop 'Cannot activate scheduled tree'
    topology=layouts(stage);regrid_stage=stage
  end subroutine
  function descriptor_for(layout,index) result(descriptor)
    type(amr_patch_tree_3d), intent(in) :: layout
    integer, intent(in) :: index
    integer :: descriptor(12)
    associate(n=>layout%nodes(index),p=>layout%nodes(index)%patch)
      descriptor=[n%parent,n%level,n%shape,p%refinement_ratio,p%coarse_i_lower,p%coarse_i_upper, &
        p%coarse_j_lower,p%coarse_j_upper,p%coarse_k_lower,p%coarse_k_upper]
    end associate
  end function
  subroutine adaptive_record(stream,reading)
    integer, intent(in) :: stream
    logical, intent(in) :: reading
    integer :: controls(6),got(6),saved_ratios(8),stored_stage,stored_count,event,node,node_count,axis,p,level
    integer :: descriptors(12,33),parents(32),refinement(32),lo(3,32),hi(3,32),previous(3),dims(3)
    real(dp) :: thresholds(8),event_time
    controls=[tag_depth,tag_tile_shape,tag_buffer,tag_maximum_nodes]
    if (reading) then
      read(stream,*,iostat=ios) got
      if (ios/=0) error stop 'Missing tagging policy'
      if (any(got/=controls)) error stop 'Tagging policy mismatch'
      read(stream,*,iostat=ios) saved_ratios(:tag_depth)
      if (ios/=0) error stop 'Missing tagging ratios'
      read(stream,*,iostat=ios) thresholds(:tag_depth)
      if (ios/=0) error stop 'Missing tagging thresholds'
      if (any(saved_ratios(:tag_depth)/=tag_ratios(:tag_depth)).or. &
        any(thresholds(:tag_depth)/=tag_thresholds(:tag_depth))) error stop 'Tagging thresholds or ratios mismatch'
      read(stream,*,iostat=ios) stored_count,stored_stage
      if (ios/=0) error stop 'Missing adaptive schedule'
      if (stored_count/=regrid_count.or.stored_stage<0.or.stored_stage>regrid_count) error stop 'Invalid adaptive stage'
    else
      stored_stage=regrid_stage
      write(stream,*,iostat=ios) controls
      if (ios/=0) error stop 'Cannot write tagging policy'
      write(stream,*,iostat=ios) tag_ratios(:tag_depth)
      if (ios/=0) error stop 'Cannot write tagging ratios'
      write(stream,'(*(es27.18e3,1x))',iostat=ios) tag_thresholds(:tag_depth)
      if (ios/=0) error stop 'Cannot write tagging thresholds'
      write(stream,*,iostat=ios) regrid_count,regrid_stage
      if (ios/=0) error stop 'Cannot write adaptive schedule'
    end if
    do event=0,regrid_count
      if (event>0) then
        if (reading) then
          read(stream,*,iostat=ios) event_time
          if (ios/=0) error stop 'Missing adaptive time'
          if (event_time/=regrid_times(event)) error stop 'Adaptive event time mismatch'
        else
          write(stream,'(es27.18e3)',iostat=ios) regrid_times(event)
          if (ios/=0) error stop 'Cannot write adaptive time'
        end if
      end if
      node_count=0
      if (event<=stored_stage) node_count=size(layouts(event)%nodes)
      if (reading) then
        read(stream,*,iostat=ios) node_count
        if (ios/=0) error stop 'Missing adaptive topology'
        if (event>stored_stage) then
          if (node_count/=0) error stop 'Unexpected future adaptive topology'
          cycle
        end if
        if (node_count<1.or.node_count>33) error stop 'Invalid adaptive node count'
        if (event>0.and.node_count>tag_maximum_nodes) error stop 'Adaptive node capacity exceeded'
        do node=1,node_count
          read(stream,*,iostat=ios) descriptors(:,node)
          if (ios/=0) error stop 'Missing adaptive node'
        end do
        if (event==0) then
          if (node_count/=size(layouts(0)%nodes)) error stop 'Initial adaptive topology mismatch'
          do node=1,node_count
            if (any(descriptors(:,node)/=descriptor_for(layouts(0),node))) error stop 'Initial adaptive node mismatch'
          end do
          cycle
        end if
        if (any(descriptors(:,1)/=descriptor_for(layouts(0),1))) error stop 'Adaptive root mismatch'
        do node=2,node_count
          p=node-1;parents(p)=descriptors(1,node);refinement(p)=descriptors(6,node)
          lo(:,p)=descriptors([7,9,11],node);hi(:,p)=descriptors([8,10,12],node)
        end do
        call make_layout(node_count-1,parents,refinement,lo,hi,layouts(event))
        do node=2,node_count
          if (any(descriptors(:,node)/=descriptor_for(layouts(event),node))) error stop 'Adaptive descriptor mismatch'
          level=layouts(event)%nodes(node)%level;p=layouts(event)%nodes(node)%parent
          if (level>tag_depth) error stop 'Adaptive depth exceeds policy'
          if (refinement(node-1)/=tag_ratios(level)) error stop 'Adaptive ratio differs from policy'
          dims=layouts(event)%nodes(p)%shape
          if (any(mod(lo(:,node-1)-1,tag_tile_shape)/=0)) error stop 'Adaptive tile alignment mismatch'
          if (any(hi(:,node-1)/=lo(:,node-1)+min(tag_tile_shape-1,dims-lo(:,node-1)))) &
            error stop 'Adaptive tile extent mismatch'
          if (node>2) then
            if (p<layouts(event)%nodes(node-1)%parent) error stop 'Adaptive parent order mismatch'
            if (p==layouts(event)%nodes(node-1)%parent) then
              previous=lo(:,node-2)
              do axis=3,1,-1
                if (previous(axis)==lo(axis,node-1)) cycle
                if (previous(axis)>lo(axis,node-1)) error stop 'Adaptive tile order mismatch'
                exit
              end do
            end if
          end if
        end do
      else
        write(stream,*,iostat=ios) node_count
        if (ios/=0) error stop 'Cannot write adaptive node count'
        do node=1,node_count
          write(stream,*,iostat=ios) descriptor_for(layouts(event),node)
          if (ios/=0) error stop 'Cannot write adaptive node'
        end do
      end if
    end do
    if (reading) call activate(stored_stage,.false.)
  end subroutine
  subroutine schedule_record(stream,reading)
    integer, intent(in) :: stream
    logical, intent(in) :: reading
    integer :: event,node,stored_count,stored_stage,node_count,descriptor(12),stored(12)
    real(dp) :: event_time
    if (reading) then
      read(stream,*,iostat=ios) stored_count,stored_stage
      if (ios/=0) error stop 'Missing restart schedule'
      if (stored_count/=regrid_count.or.stored_stage<0.or.stored_stage>regrid_count) &
        error stop 'Invalid restart schedule count or stage'
    else
      write(stream,*,iostat=ios) regrid_count,regrid_stage
      if (ios/=0) error stop 'Cannot write schedule count'
    end if
    do event=0,regrid_count
      if (event>0) then
        if (reading) then
          read(stream,*,iostat=ios) event_time
          if (ios/=0) error stop 'Missing event time'
          if (event_time/=regrid_times(event)) error stop 'Restart event time mismatch'
        else
          write(stream,'(es27.18e3)',iostat=ios) regrid_times(event)
          if (ios/=0) error stop 'Cannot write event time'
        end if
      end if
      node_count=size(layouts(event)%nodes)
      if (reading) then
        read(stream,*,iostat=ios) stored_count
        if (ios/=0) error stop 'Missing event topology'
        if (stored_count/=node_count) error stop 'Restart event node count mismatch'
      else
        write(stream,*,iostat=ios) node_count
        if (ios/=0) error stop 'Cannot write node count'
      end if
      do node=1,node_count
        associate(n=>layouts(event)%nodes(node),p=>layouts(event)%nodes(node)%patch)
          descriptor=[n%parent,n%level,n%shape,p%refinement_ratio,p%coarse_i_lower,p%coarse_i_upper, &
            p%coarse_j_lower,p%coarse_j_upper,p%coarse_k_lower,p%coarse_k_upper]
        end associate
        if (reading) then
          read(stream,*,iostat=ios) stored
          if (ios/=0) error stop 'Missing event node'
          if (any(stored/=descriptor)) error stop 'Restart event topology mismatch'
        else
          write(stream,*,iostat=ios) descriptor
          if (ios/=0) error stop 'Cannot write event node'
        end if
      end do
    end do
    if (reading) call activate(stored_stage,.false.)
  end subroutine
#ifdef PELEF_SELECTED_EB_TREE
  logical function chemical_elements_valid(initial_values,changes,accepted_steps) result(valid)
    real(dp), intent(in) :: initial_values(:),changes(:)
    integer, intent(in) :: accepted_steps
    real(dp) :: element_scale(element_count),element_change(element_count),roundoff(element_count)
    integer :: e
    element_scale=matmul(element_atoms,abs(initial_values(6:))/species%molecular_weight)
    element_change=matmul(element_atoms,changes(6:)/species%molecular_weight)
    ! A zero initial element still sees floating-point species mass closure.
    ! Bound that absolute error by total species mass and atoms per unit mass;
    ! retain the relative tolerance for each populated elemental inventory.
    do e=1,element_count
      roundoff(e)=64._dp*epsilon(1._dp)*real(max(1,accepted_steps),dp)*sum(abs(initial_values(6:)))* &
        maxval(element_atoms(e,:)/species%molecular_weight)
    end do
    valid=all(ieee_is_finite(element_change)).and.all(ieee_is_finite(element_scale)).and.all(ieee_is_finite(roundoff))
    if (valid) valid=all(abs(element_change)<=2.e-10_dp*element_scale+roundoff)
    if (.not.valid) then
      do e=1,element_count
        write(*,'(a,1x,2(es24.16,1x))') trim(element_names(e)),element_change(e),element_scale(e)
      end do
    end if
  end function chemical_elements_valid
#endif

  subroutine balance()
#ifdef PELEF_SELECTED_EB_TREE
    if (.not.chemical_elements_valid(record%initial,record%chemical,record%steps)) &
      error stop 'Chemical element conservation failed'
#endif
#ifdef PELEF_MPI_EB_TREE
    if (distributed) then
      call integrate_mpi_eb_tree_3d(topology,mpi_owners,mpi_geometry,mpi_state,current,MPI_COMM_WORLD,ok)
    else
      call integrate_reactive_eb_tree_3d(tree,current,ok)
    end if
#else
    call integrate_reactive_eb_tree_3d(tree,current,ok)
#endif
    if (.not.ok) error stop 'Inventory failed'
#ifdef PELEF_MPI_EB_SPRAY
    call add_particle_inventory(current)
#endif
    error=maxval(abs(current-record%initial-record%boundary-record%wall-record%chemical)/ &
      max(1._dp,record%scale,abs(current),abs(record%boundary)+abs(record%wall)+abs(record%chemical)))
    if (.not.ieee_is_finite(error).or.error>5.e-10_dp) error stop 'Composite budget failed'
  end subroutine
  subroutine envelope(stream,reading,header_only)
    integer, intent(in) :: stream
    logical, intent(in) :: reading
    logical, optional, intent(in) :: header_only
    character(len=64) :: magic,integrator
    logical :: automatic
    real(dp) :: interval,stored_hydro_cfl,stored_transport_cfl
    if (reading) then
      read(stream,'(a)',iostat=ios) magic
      if (ios/=0) error stop 'Missing tree application header'
      if (magic=='PELEF_EB_TREE_SUBCYCLE_V1') then
        if (.not.subcycling) error stop 'Restart subcycling mode mismatch'
        read(stream,'(a)',iostat=ios) magic
        if (ios/=0) error stop 'Missing subcycled application header'
      else if (subcycling) then
        error stop 'Restart subcycling mode mismatch'
      end if
      if (magic/='PELEF_EB_TREE_APPLICATION_V1'.and.magic/='PELEF_EB_TREE_APPLICATION_V2'.and. &
        magic/='PELEF_EB_TREE_APPLICATION_V3'.and.magic/='PELEF_EB_TREE_APPLICATION_V4') &
        error stop 'Wrong tree application header'
      if ((magic=='PELEF_EB_TREE_APPLICATION_V4').neqv.temperature_adaptive) error stop 'Adaptation mode mismatch'
      if (magic=='PELEF_EB_TREE_APPLICATION_V1'.and.automatic_timestep) error stop 'Restart policy mismatch'
      read(stream,*,iostat=ios) interval,integrator
      if (ios/=0) error stop 'Missing restart policy'
      if (interval/=time_step.or.integrator/=chemistry_integrator) error stop 'Restart policy mismatch'
      if (magic/='PELEF_EB_TREE_APPLICATION_V1') then
        read(stream,*,iostat=ios) automatic,stored_hydro_cfl,stored_transport_cfl
        if (ios/=0) error stop 'Missing automatic interval policy'
        if (automatic.neqv.automatic_timestep) error stop 'Restart interval mode mismatch'
        if (stored_hydro_cfl/=hydro_cfl.or.stored_transport_cfl/=transport_cfl) error stop 'Restart CFL mismatch'
      end if
    else
      if (subcycling) then
        write(stream,'(a)',iostat=ios) 'PELEF_EB_TREE_SUBCYCLE_V1'
        if (ios/=0) error stop 'Cannot write subcycling mode'
      end if
      magic='PELEF_EB_TREE_APPLICATION_V2'
      if (regrid_count>0) magic='PELEF_EB_TREE_APPLICATION_V3'
      if (temperature_adaptive) magic='PELEF_EB_TREE_APPLICATION_V4'
      write(stream,'(a)',iostat=ios) trim(magic)
      if (ios/=0) error stop 'Cannot write application header'
      write(stream,'(es27.18e3,1x,a)',iostat=ios) time_step,trim(chemistry_integrator)
      if (ios/=0) error stop 'Cannot write policy'
      write(stream,'(l1,2(1x,es27.18e3))',iostat=ios) automatic_timestep,hydro_cfl,transport_cfl
      if (ios/=0) error stop 'Cannot write automatic interval policy'
    end if
    if (temperature_adaptive) then
      call adaptive_record(stream,reading)
    else if (magic=='PELEF_EB_TREE_APPLICATION_V3') then
      call schedule_record(stream,reading)
    else if (regrid_count/=0) then
      error stop 'Restart schedule missing'
    end if
    if (present(header_only)) then
      if (header_only) return
    end if
    call reactive_eb_tree_checkpoint_record_3d(stream,reading,species,reactions,transport,config,tree,record,ok,message, &
      bundle_identity,checkpoint_integrator,base_mole_fractions,element_names,element_atoms)
    if (.not.ok) then
      print *,trim(message)
      error stop 'Invalid tree checkpoint record'
    end if
    if (reading) then
      if (regrid_stage/=count(regrid_times(:regrid_count)<=record%time)) error stop 'Restart stage/time mismatch'
    end if
  end subroutine
  subroutine restore()
    character(len=1) :: extra
    open(newunit=unit,file=trim(restart_file),status='old',action='read',access='stream',form='formatted',iostat=ios)
    if (ios/=0) error stop 'Cannot open restart'
#ifdef PELEF_MPI_EB_TREE
    call envelope(unit,.true.,header_only=.true.)
    call read_mpi_eb_tree_checkpoint_3d(trim(restart_file),unit,species,reactions,transport,config,topology,mpi_owners, &
      mpi_geometry,mpi_state,mpi_temperature,record,MPI_COMM_WORLD,ok, &
      bundle_identity,checkpoint_integrator,base_mole_fractions,element_names,element_atoms)
    if (.not.ok) error stop 'Invalid distributed tree checkpoint record'
    if (regrid_stage/=count(regrid_times(:regrid_count)<=record%time)) error stop 'Restart stage/time mismatch'
#else
    call envelope(unit,.true.)
    read(unit,'(a)',iostat=ios) extra
    if (ios/=iostat_end) error stop 'Trailing restart data'
#endif
    close(unit,iostat=ios)
    if (ios/=0) error stop 'Restart close failed'
    if (record%time>final_time.or.record%steps>maximum_steps) error stop 'Restart beyond requested limits'
    call balance()
  end subroutine
  subroutine save()
    integer :: scratch,destination,amount
    character(len=4096) :: chunk
    open(newunit=scratch,status='scratch',form='formatted',action='readwrite',iostat=ios)
    if (ios/=0) error stop 'Cannot stage checkpoint'
#ifdef PELEF_MPI_EB_TREE
    call envelope(scratch,.false.,header_only=.true.)
    call write_mpi_eb_tree_checkpoint_3d(trim(checkpoint_file),scratch,species,reactions,transport,config,topology, &
      mpi_owners,mpi_geometry,mpi_state,mpi_temperature,record,MPI_COMM_WORLD,ok, &
      bundle_identity,checkpoint_integrator,base_mole_fractions,element_names,element_atoms)
    close(scratch,iostat=ios)
    if (.not.ok.or.ios/=0) error stop 'MPI checkpoint save failed'
#else
    call envelope(scratch,.false.)
    rewind(scratch,iostat=ios)
    if (ios/=0) error stop 'Checkpoint rewind failed'
    open(newunit=destination,file=trim(checkpoint_file),status='replace',action='write',iostat=ios)
    if (ios/=0) error stop 'Cannot open checkpoint destination'
    do
      read(scratch,'(a)',advance='no',size=amount,iostat=ios) chunk
      if (ios==iostat_end) exit
      if (ios/=0.and.ios/=iostat_eor) error stop 'Checkpoint staging read failed'
      write(destination,'(a)',advance='no',iostat=status) chunk(:amount)
      if (status/=0) error stop 'Checkpoint publication failed'
      if (ios==iostat_eor) then
        write(destination,'(a)',iostat=status) ''
        if (status/=0) error stop 'Checkpoint newline failed'
      end if
    end do
    close(destination,iostat=ios)
    if (ios/=0) error stop 'Checkpoint close failed'
    close(scratch,iostat=ios)
    if (ios/=0) error stop 'Checkpoint staging close failed'
#endif
  end subroutine
  subroutine output()
    integer :: i,j,k,c,child,lo(3),hi(3)
    logical, allocatable :: exposed(:,:,:)
#ifdef PELEF_MPI_EB_TREE
    call write_mpi_eb_tree_csv_3d(trim(output_file),topology,mpi_owners,species,mpi_state,mpi_temperature, &
      mpi_geometry,MPI_COMM_WORLD,ok)
    if (.not.ok) error stop 'MPI CSV write failed'
    if (mpi_rank/=0) return
#else
    open(newunit=unit,file=trim(output_file),status='replace',action='write',iostat=ios)
    if (ios/=0) error stop 'Cannot open CSV'
    write(unit,'(a)',advance='no') 'node,level,i,j,k,fluid_volume,temperature,rho,rhou,rhov,rhow,rhoE'
    do c=1,size(species)
      write(unit,'(a)',advance='no') ',rhoY_'//trim(species(c)%name)
    end do
    write(unit,*)
    do node=1,size(tree%nodes)
      n=topology%nodes(node)%shape
      allocate(exposed(n(1),n(2),n(3)));exposed=tree%nodes(node)%geometry%volume_fraction>0._dp
      do child=node+1,size(tree%nodes)
        if (topology%nodes(child)%parent/=node) cycle
        associate(box=>topology%nodes(child)%patch)
          lo=[box%coarse_i_lower,box%coarse_j_lower,box%coarse_k_lower]
          hi=[box%coarse_i_upper,box%coarse_j_upper,box%coarse_k_upper]
        end associate
        exposed(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.false.
      end do
      associate(g=>tree%nodes(node)%geometry)
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if (.not.exposed(i,j,k)) cycle
        write(unit,'(i0,4(",",i0),*(:,",",es26.17e3))') node,topology%nodes(node)%level,i,j,k, &
          g%volume_fraction(i,j,k)*g%dx*g%dy*g%dz,tree%nodes(node)%temperature(i,j,k),tree%nodes(node)%state(:,i,j,k)
      end do
      end do
      end do
      end associate
      deallocate(exposed)
    end do
    close(unit,iostat=ios)
    if (ios/=0) error stop 'CSV close failed'
#endif
    open(newunit=unit,file=trim(budget_path),status='replace',action='write',iostat=ios)
    if (ios/=0) error stop 'Cannot open budget CSV'
    write(unit,'(a)') 'component,initial,final,boundary,wall,chemistry,residual'
    do c=1,nv
      write(unit,'(i0,6(",",es26.17e3))') c,record%initial(c),current(c),record%boundary(c), &
        record%wall(c),record%chemical(c),current(c)-record%initial(c)-record%boundary(c)-record%wall(c)-record%chemical(c)
    end do
    close(unit,iostat=ios)
    if (ios/=0) error stop 'Budget close failed'
  end subroutine
end program pelef_reactive_eb_tree_3d
