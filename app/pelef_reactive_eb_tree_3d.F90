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
  use, intrinsic :: iso_c_binding, only: c_int,c_char,c_null_char
  use spray_mpi_eb_particle_checkpoint_mod, only: mpi_eb_particle_checkpoint
  use spray_injection_profile_mod, only: integrate_spray_injection_profile,read_spray_injection_profile
  use spray_mpi_eb_remap_mod, only: remap_mpi_eb_spray
  use spray_mpi_amr_balance_mod, only: rebalance_mpi_amr_spray
  use spray_parcel_mod
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options,valid_wall_options
  use spray_mixture_mod, only: liquid_mixture_properties,mixture_evolution_id
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
  use mixture_thermo_mod, only: mixture_specific_gas_constant
#endif
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar,reactive_conserved_to_primitive,reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,initialize_amr_patch_tree_3d,amr_tree_integer_layout_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config,build_configured_eb_geometry_3d
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    advance_reactive_eb_tree_full_3d,integrate_reactive_eb_tree_3d,reactive_eb_tree_timestep_3d, &
    build_eb_tree_extended_geometry_3d, &
    remap_reactive_eb_tree_3d,tag_reactive_eb_tree_temperature_3d
  use eb_reactive_tree_checkpoint_3d_mod, only: eb_tree_checkpoint_3d,reactive_eb_tree_checkpoint_record_3d
  use selected_composition_mod, only: selected_composition_max_species,selected_composition_name_length, &
    resolve_selected_composition
  implicit none
#ifdef PELEF_MPI_EB_SPRAY
  interface
    function publish_link(old,new) bind(C,name='link') result(code)
      import c_int,c_char
      character(c_char), intent(in) :: old(*),new(*)
      integer(c_int) :: code
    end function
  end interface
#endif
  character(len=:), allocatable :: element_names(:),bundle_identity,checkpoint_integrator
  real(dp), allocatable :: element_atoms(:,:),base_mole_fractions(:)
#ifdef PELEF_MPI_EB_SPRAY
  type eb_injection_port
    type(spray_parcel), allocatable :: parcels(:)
    real(dp), allocatable :: fractions(:)
  end type
  type(eb_injection_port), allocatable :: injection_ports(:)
  type(spray_parcel), allocatable :: drops(:)
  integer :: port_count=1,event_ports(1024)=1
  integer :: injection_events=1,next_injection=1
  integer(int64) :: last_particle_id=0,split_created_total=0,step_split_created=0
  real(dp) :: injection_times(1024)=0,injection_scales(1024)=1,event_factors(1024)=1
  real(dp) :: event_velocity_scales(1024)=1
  type(spray_liquid), allocatable :: liquids(:)
  type(spray_options) :: particle_options
  type(spray_wall_options) :: walls
  type(spray_wall_film), allocatable :: films(:)
  real(dp), allocatable :: film_mass(:,:)
  real(dp) :: acceleration(3)=0
  logical :: wall_impacts=.false.,particle_outflow=.false.,model_rejected
  real(dp), allocatable :: nodal_phi(:,:,:),component_mass(:,:)
  character(len=:), allocatable :: spray_context
  character(len=1024) :: spray_input
  integer :: particle_steps,restart_layout_stage=0
#endif
  integer :: inflow_composition_count(6)=0
  character(len=selected_composition_name_length) :: &
    inflow_composition_species(selected_composition_max_species,6)=''
  real(dp) :: inflow_composition_mole_fractions(selected_composition_max_species,6)=0._dp
  real(dp), allocatable :: inlet_mass_fractions(:,:)
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
  real(dp) :: inflow_mach(6)=-1._dp,inflow_temperature(6)=0._dp,inflow_pressure(6)=0._dp
  real(dp), allocatable :: active_mach(:),active_temperature(:),active_pressure(:),inlet_state(:,:),inlet_primitive(:)
  real(dp) :: wall_temperature(6)=0
  real(dp), allocatable :: active_wall_temperature(:)
  character(len=16) :: boundary_faces(6)='legacy'
  character(len=16), allocatable :: active_faces(:)
  real(dp) :: inflow_density(6)=.31_dp,inflow_velocity(3,6)=0,inlet_t,inlet_sound
  integer :: inlet_face
  real(dp) :: regular_density=.31_dp,cut_density=.93_dp,initial_pressure=135000._dp
  real(dp) :: hydro_cfl=.5_dp,transport_cfl=.35_dp,difmag=.1_dp
  real(dp) :: time_step=1.e-9_dp,final_time=2.e-9_dp,rtol=2.e-7_dp,atol=1.e-12_dp
  logical :: chemistry_enabled=.true.,transport_enabled=.true.,stop_after_checkpoint=.false.
  character(len=16) :: geometry_type='plane',plane_axis='z',extrusion_axis='z'
  character(len=16) :: riemann_solver='rusanov',hydro_method='pcm'
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
    inflow_mach,inflow_temperature,inflow_pressure,difmag,hydro_method,boundary_faces,inflow_density,inflow_velocity, &
    inflow_composition_count,inflow_composition_species,inflow_composition_mole_fractions, &
    regular_density,cut_density,initial_pressure,wall_temperature, &
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
  type(eb_geometry_3d), allocatable :: geometry_levels(:),extended_levels(:),extended_geometry(:)
  integer, allocatable :: cached_descriptors(:,:)
  type(eb_reactive_tree_3d) :: tree
  type(eb_tree_checkpoint_3d) :: record
  type(eb_geometry_3d) :: finest,coarse
  type(reactive_eb_3d_config) :: config,root_config
  integer(int64), allocatable :: origin(:,:),domain(:,:)
  integer, allocatable :: ratios(:),shared_ratios(:)
  real(dp), allocatable :: force_total(:),step_force(:)
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:),current(:),boundary(:),wall(:),chemical(:),injected(:)
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
  if(.not.all(ieee_is_finite([regular_density,cut_density,initial_pressure]))) error stop 'Invalid initial state'
  if(min(regular_density,cut_density,initial_pressure)<=0) error stop 'Nonpositive initial state'
  if(hydro_method/='pcm'.and.hydro_method/='plm_ctu') error stop 'Invalid hydro method'
  if(hydro_method=='plm_ctu'.and.riemann_solver/='pelec') error stop 'CTU requires pelec solver'
  if (.not.ieee_is_finite(difmag).or.difmag<0) error stop 'Invalid artificial-viscosity coefficient'
  if (.not.all(ieee_is_finite(inflow_mach)).or.any(inflow_mach>100)) error stop 'Invalid inlet Mach'
  if (.not.all(ieee_is_finite(inflow_temperature))) error stop 'Invalid inlet temperature'
  if (.not.all(ieee_is_finite(inflow_pressure))) error stop 'Invalid inlet pressure'
  if (any(inflow_pressure<0).or.any(inflow_pressure>1e15_dp)) error stop 'Invalid inlet pressure'
  if(any(boundary_faces/='legacy')) then
    if(any(boundary_faces/='inflow'.and.boundary_faces/='outflow'.and.boundary_faces/='periodic'.and. &
      boundary_faces/='slip_wall'.and.boundary_faces/='no_slip_wall'.and.boundary_faces/='isothermal_wall')) &
      error stop 'Invalid tree boundary faces'
    do inlet_face=1,5,2
      if((boundary_faces(inlet_face)=='periodic').neqv.(boundary_faces(inlet_face+1)=='periodic')) &
        error stop 'Unpaired periodic tree faces'
    end do
    if(any(inflow_mach>=0).or.any(inflow_temperature/=0)) error stop 'Conflicting tree inlet modes'
    if(.not.all(ieee_is_finite(inflow_density)).or..not.all(ieee_is_finite(inflow_velocity))) &
      error stop 'Invalid prescribed tree inlet'
    if(any(inflow_density<=0).or.any(inflow_pressure<=0.and.boundary_faces=='inflow')) &
      error stop 'Nonpositive prescribed tree inlet'
    active_faces=boundary_faces
  end if
  if(.not.all(ieee_is_finite(wall_temperature))) error stop 'Invalid tree wall temperature'
  if(any(boundary_faces=='isothermal_wall')) then
    if(any(wall_temperature<=0)) error stop 'Isothermal tree walls require six positive temperatures'
    if(.not.transport_enabled) error stop 'Isothermal tree walls require transport'
    active_wall_temperature=wall_temperature
  else if(any(wall_temperature/=0)) then
    error stop 'Unused tree wall temperatures'
  end if
  if (any(inflow_temperature<=0.and.inflow_mach>=0)) error stop 'Active inlet requires positive temperature'
  if (.not.allocated(active_faces).and.any((inflow_temperature/=0.or.inflow_pressure/=0).and.inflow_mach<0)) &
    error stop 'Inactive inlet has unused temperature or pressure'
  if (any(inflow_mach>=0)) then
    active_mach=inflow_mach;active_temperature=inflow_temperature;active_pressure=inflow_pressure
  end if
#ifdef PELEF_MPI_EB_SPRAY
  if(geometry_type/='level_set') error stop 'MPI EB spray requires original nodal geometry'
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
  config%redistribution='none';config%riemann_solver=riemann_solver;config%difmag=difmag
  config%hydro_method=hydro_method
  if(allocated(active_faces)) config%boundary_faces=active_faces
  config%regular_density=regular_density;config%cut_density=cut_density;config%initial_pressure=initial_pressure
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
  if(hydro_method=='plm_ctu') then
    common=6
    do depth=2,size(shared_ratios)
      if(common>huge(1)/shared_ratios(depth)) error stop 'CTU padding too large'
      common=common*shared_ratios(depth)
    end do
    allocate(extended_levels(size(shared_ratios)))
    call build_configured_eb_geometry_3d(config,extended_levels(size(shared_ratios)),ok,halo=common)
    if(.not.ok) error stop 'Cannot build extended finest geometry'
    do depth=size(shared_ratios),2,-1
      call coarsen_eb_geometry_3d(extended_levels(depth),shared_ratios(depth),extended_levels(depth-1),ok)
      if(.not.ok) error stop 'Cannot coarsen extended geometry'
    end do
  end if
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
  call prepare_inlet_compositions()
  nv=reactive_nvar(size(species));n=root_shape
  allocate(current(nv),boundary(nv),wall(nv),chemical(nv),injected(nv))
  allocate(force_total(nv),step_force(nv));force_total=0;step_force=0
  injected=0
  root_config=config;root_config%nx=n(1);root_config%ny=n(2);root_config%nz=n(3)
#ifdef PELEF_MPI_EB_TREE
  call initialize_owned()
  call integrate_mpi_eb_tree_3d(topology,mpi_owners,mpi_geometry,mpi_state,current,MPI_COMM_WORLD,ok)
#else
  allocate(u(nv,n(1),n(2),n(3)),t(n(1),n(2),n(3)))
  call initialize_reactive_eb_density_sheet_3d(species,root_config,coarse,u,t,ok,base_mole_fractions)
  if (.not.ok) error stop 'Initial state failed'
  if(allocated(active_faces)) then
    call prepare_prescribed_inlet(u(:,1,1,1),t(1,1,1),ok)
    if(.not.ok) error stop 'Invalid prescribed inlet EOS'
  end if
  call initialize_reactive_eb_tree_3d(species,topology,geometry_levels(topology%level_count()),u,t,tree,ok)
  if (.not.ok) error stop 'Tree initialization failed'
  call integrate_reactive_eb_tree_3d(tree,current,ok)
#endif
  if (.not.ok) error stop 'Initial inventory failed'
#ifdef PELEF_MPI_EB_SPRAY
  call initialize_particles()
  call add_particle_inventory(current,drops,component_mass)
#endif
  record%initial=current;record%scale=abs(current);record%boundary=0._dp*current
  record%wall=0._dp*current;record%chemical=0._dp*current
#ifdef PELEF_MPI_EB_SPRAY
  if(len_trim(restart_file)>0) call spray_checkpoint(.false.)
#else
  if (len_trim(restart_file)>0) call restore()
#endif
  do
#ifdef PELEF_MPI_EB_SPRAY
    call inject_due_events()
#endif
    if(record%time>=final_time.or.record%steps>=maximum_steps) exit
    call refresh_ctu_geometry()
    dt=min(time_step,final_time-record%time)
#ifdef PELEF_MPI_EB_SPRAY
    if(next_injection<=injection_events) dt=min(dt,injection_times(next_injection)-record%time)
#endif
    if (automatic_timestep) then
#ifdef PELEF_MPI_EB_TREE
      call mpi_eb_tree_timestep_3d(species,transport,topology,mpi_owners,mpi_geometry,mpi_state,mpi_temperature, &
        hydro_cfl,transport_cfl,transport_enabled,.true.,.true.,.true.,estimated_dt,MPI_COMM_WORLD,ok,subcycling, &
        inflow_mach=active_mach,inflow_temperature=active_temperature,inflow_pressure=active_pressure, &
        boundary_faces=active_faces,inflow_state=inlet_state,wall_temperature=active_wall_temperature)
#else
      call reactive_eb_tree_timestep_3d(species,transport,tree,hydro_cfl,transport_cfl,transport_enabled, &
        .true.,.true.,.true.,estimated_dt,ok,subcycling, &
        inflow_mach=active_mach,inflow_temperature=active_temperature,inflow_pressure=active_pressure, &
        boundary_faces=active_faces,inflow_state=inlet_state,wall_temperature=active_wall_temperature)
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
      if(wall_impacts) then
      call advance_mpi_eb_spray_reactive(MPI_COMM_WORLD,topology,mpi_owners,mpi_geometry,nodal_phi, &
        species,reactions,transport,liquids,particle_options,mpi_state,mpi_temperature,drops,component_mass, &
        domain_lower,(domain_upper-domain_lower)/real(root_shape,dp),dt,riemann_solver, &
        chemistry_enabled,transport_enabled,.true.,.true.,.true.,.true.,rtol,atol,particle_steps, &
        theta,boundary,wall,chemical,ok,chemistry_integrator,subcycling, &
        films=films,film_mass=film_mass,film_settings=walls,last_id=last_particle_id,model_rejected=model_rejected, &
        split_created=step_split_created,particle_outflow=particle_outflow, &
        acceleration=acceleration,source_change=step_force, &
        inflow_mach=active_mach,inflow_temperature=active_temperature,inflow_pressure=active_pressure, &
        difmag=difmag,extended_geometry=extended_geometry, &
        boundary_faces=active_faces,inflow_state=inlet_state,wall_temperature=active_wall_temperature)
      if(model_rejected) error stop 'EB wall impact outside admissible empirical model range'
      else
      call advance_mpi_eb_spray_reactive(MPI_COMM_WORLD,topology,mpi_owners,mpi_geometry,nodal_phi, &
        species,reactions,transport,liquids,particle_options,mpi_state,mpi_temperature,drops,component_mass, &
        domain_lower,(domain_upper-domain_lower)/real(root_shape,dp),dt,riemann_solver, &
        chemistry_enabled,transport_enabled,.true.,.true.,.true.,.true.,rtol,atol,particle_steps, &
        theta,boundary,wall,chemical,ok,chemistry_integrator,subcycling,last_id=last_particle_id, &
        split_created=step_split_created,particle_outflow=particle_outflow, &
        acceleration=acceleration,source_change=step_force, &
        inflow_mach=active_mach,inflow_temperature=active_temperature,inflow_pressure=active_pressure, &
        difmag=difmag,extended_geometry=extended_geometry, &
        boundary_faces=active_faces,inflow_state=inlet_state,wall_temperature=active_wall_temperature)
      end if
#else
      call advance_mpi_eb_tree_full_3d(species,reactions,transport,topology,mpi_owners,mpi_geometry, &
        riemann_solver,dt,chemistry_enabled,rtol,atol,transport_enabled,.true.,.true.,.true.,.true., &
        mpi_state,mpi_temperature,boundary,wall,chemical,theta,MPI_COMM_WORLD,ok,chemistry_integrator,subcycling, &
        inflow_mach=active_mach,inflow_temperature=active_temperature,inflow_pressure=active_pressure, &
        difmag=difmag,extended_geometry=extended_geometry, &
        boundary_faces=active_faces,inflow_state=inlet_state,wall_temperature=active_wall_temperature)
#endif
      if (sample_work) call end_tree_work(measured_cost)
#else
      call advance_reactive_eb_tree_full_3d(species,reactions,transport,riemann_solver,dt,chemistry_enabled, &
        rtol,atol,transport_enabled,.true.,.true.,.true.,.true.,tree,boundary,wall,chemical,theta,ok, &
        chemistry_integrator,subcycling, &
        inflow_mach=active_mach,inflow_temperature=active_temperature,inflow_pressure=active_pressure, &
        difmag=difmag,extended_geometry=extended_geometry, &
        boundary_faces=active_faces,inflow_state=inlet_state,wall_temperature=active_wall_temperature)
#endif
      if (ok) exit
      if (retry==maximum_step_retries) error stop 'Tree step retry limit exhausted'
      dt=.5_dp*dt
    end do
#ifdef PELEF_MPI_EB_SPRAY
    split_created_total=split_created_total+step_split_created
#endif
    force_total=force_total+step_force
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
#ifdef PELEF_MPI_EB_SPRAY
      call rebalance_mpi_amr_spray(MPI_COMM_WORLD,topology,mpi_owners,nv,mpi_state,mpi_temperature, &
        drops,component_mass,domain_lower,(domain_upper-domain_lower)/real(root_shape,dp), &
        [.false.,.false.,.false.],measured_cost,rebalance_minimum_gain,particle_options%maximum_parcels, &
        redistributed,ok,wall_films=films,film_mass=film_mass,maximum_films=walls%maximum_films,nodal_phi=nodal_phi)
#else
      call rebalance_mpi_amr_tree_3d(topology,mpi_owners,nv,mpi_state,mpi_temperature, &
        measured_cost,rebalance_minimum_gain,MPI_COMM_WORLD,redistributed,ok)
#endif
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
#ifdef PELEF_MPI_EB_SPRAY
        call spray_checkpoint(.true.)
#else
        call save()
#endif
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
    integer :: nc,number,q,c,stream,code,closed,port,port_numbers(16),port_profile_points(16)
    real(dp) :: port_profile_times(64,16),port_profile_rates(64,16),injection_end_times(1024),profile_mass
    real(dp) :: port_centers(3,16),port_axes(3,16),port_masses(16),port_speeds(16)
    real(dp) :: port_angles(16),port_diameters(16),port_temperatures(16),port_fractions(16,16)
    real(dp) :: port_radii(16)=0
    type(spray_nozzle_model) :: port_nozzles(16)
    type(spray_size_options) :: port_sizes(16)
    logical :: port_velocity_from_rate(16)=.false.
    real(dp) :: port_discharge_coefficients(16)=1,template_speed,denominator
    character(len=16) :: mixture_density_models(16)='constant',mixture_cp_models(16)='constant'
    character(len=16) :: mixture_vapor_pressure_models(16)='clausius',mixture_boiling_models(16)='reject'
    real(dp) :: mixture_density_coefficients(4,16)=0,mixture_density_temperature_ranges(2,16)=0
    real(dp) :: mixture_cp_coefficients(4,16)=0,mixture_cp_temperature_ranges(2,16)=0
    real(dp) :: mixture_antoine_coefficients(4,16)=0,mixture_antoine_temperature_ranges(2,16)=0
    real(dp) :: mixture_molecular_weights(16)=0
    real(dp) :: mixture_critical_temperatures(16)=0,mixture_normal_boiling_temperatures(16)=0
    character(len=4096) :: context_line
    character(len=1024) :: profile_files(16)
    namelist /eb_spray_rate_files/ profile_files
    character(len=24) :: vapor_names(16)
    real(dp) :: densities(16),heat_capacities(16),reference_temperatures(16)
    real(dp) :: reference_pressures(16),latent_heats(16),fractions(16)
    real(dp) :: total_mass,center(3),axis(3),speed,half_angle,diameter,liquid_temperature,rho,cp,e
    type(spray_liquid) :: effective
    namelist /eb_spray/ nc,number,vapor_names,densities,heat_capacities,reference_temperatures, &
      reference_pressures,latent_heats,fractions,total_mass,center,axis,speed,half_angle,diameter, &
      liquid_temperature,particle_options,injection_events,injection_times,injection_scales,port_count,event_ports, &
      port_numbers,port_centers,port_axes,port_masses,port_speeds,port_angles,port_diameters,port_temperatures,port_fractions, &
      port_profile_points,port_profile_times,port_profile_rates,injection_end_times,wall_impacts,walls, &
      particle_outflow,acceleration, &
      port_radii,port_nozzles,port_sizes,port_velocity_from_rate,port_discharge_coefficients, &
      mixture_density_models,mixture_density_coefficients,mixture_density_temperature_ranges, &
      mixture_cp_models,mixture_cp_coefficients,mixture_cp_temperature_ranges, &
      mixture_vapor_pressure_models,mixture_antoine_coefficients,mixture_antoine_temperature_ranges, &
      mixture_boiling_models,mixture_critical_temperatures,mixture_molecular_weights, &
      mixture_normal_boiling_temperatures
    port_profile_points=0;port_profile_times=0;port_profile_rates=0;injection_end_times=-1
    port_count=1;event_ports=1;port_numbers=-1
    port_centers=huge(1._dp);port_axes=huge(1._dp)
    port_masses=-1;port_speeds=-1;port_angles=-1;port_diameters=-1;port_temperatures=-1;port_fractions=-1
    nc=1;number=1;vapor_names='H2O';densities=1000;heat_capacities=4184;reference_temperatures=373.15_dp
    reference_pressures=101325;latent_heats=2.257e6_dp;fractions=0;fractions(1)=1
    total_mass=1e-10_dp;center=0;axis=[1._dp,0._dp,0._dp];speed=1;half_angle=0
    diameter=1e-5_dp;liquid_temperature=350;particle_options=spray_options()
    open(newunit=stream,file=trim(spray_input),status='old',action='read',iostat=code)
    if(code/=0) error stop 'Cannot read MPI EB spray input'
    read(stream,nml=eb_spray,iostat=code)
    if(code/=0) error stop 'Invalid MPI EB spray input'
    if(allocated(active_faces).and..not.particle_outflow) &
      error stop 'Explicit MPI EB spray boundaries require particle_outflow'
    profile_files=''
    rewind(stream)
    read(stream,nml=eb_spray_rate_files,iostat=code)
    if(code>0) error stop 'Invalid EB spray rate files input'
    close(stream,iostat=closed)
    if(closed/=0) error stop 'Invalid MPI EB spray input'
    if(port_count<1.or.port_count>16) error stop 'Invalid port count'
    if(any(profile_files(port_count+1:)/='')) error stop 'Unused EB rate files'
    do port=1,port_count
      if(len_trim(profile_files(port))==0) cycle
      if(port_profile_points(port)/=0.or.any(port_profile_times(:,port)/=0).or. &
        any(port_profile_rates(:,port)/=0)) error stop 'Rate file conflicts with inline profile'
      call read_spray_injection_profile(trim(profile_files(port)),port_profile_times(:,port), &
        port_profile_rates(:,port),port_profile_points(port),ok)
      if(.not.ok) error stop 'Invalid EB spray rate file'
    end do
    if(nc<1.or.nc>16.or.number<1.or.number>125000) error stop 'Invalid MPI EB spray dimensions'
    if(.not.all(ieee_is_finite(fractions)).or.any(fractions<0)) error stop 'Invalid liquid fractions'
    if(abs(sum(fractions(:nc))-1)>1e-12_dp.or.any(fractions(nc+1:)/=0)) error stop 'Invalid fraction sum'
    if(.not.valid_spray_options(particle_options)) &
      error stop 'Unsupported particle options'
    if(number>particle_options%maximum_parcels) error stop 'Particle capacity exceeded'
    if(.not.all(ieee_is_finite(acceleration))) error stop 'Invalid MPI EB acceleration'
    if(any(abs(acceleration)>1e6_dp)) error stop 'Invalid MPI EB acceleration'
    if(wall_impacts) then
      if(.not.valid_wall_options(walls,[character(len=16)::'slip_wall','slip_wall','slip_wall', &
        'slip_wall','slip_wall','slip_wall'])) error stop 'Invalid EB wall options'
      if(any(walls%temperature/=walls%temperature(1))) error stop 'EB wall temperature must be uniform'
      if(particle_options%liquid_viscosity<=0.or.particle_options%surface_tension<=0) &
        error stop 'EB impacts require viscosity and surface tension'
    else if(walls%evolve) then
      error stop 'EB film evolution requires wall impacts'
    end if
    if(injection_events<1.or.injection_events>1024) error stop 'Invalid injection event count'
    if(.not.all(ieee_is_finite(injection_times)).or..not.all(ieee_is_finite(injection_scales))) &
      error stop 'Nonfinite injection schedule'
    if(any(injection_times<0).or.any(injection_times>1e6_dp).or.any(injection_scales<0).or. &
      any(injection_scales>1e6_dp)) error stop 'Invalid injection schedule'
    do q=2,injection_events
      if(injection_times(q)<injection_times(q-1)) error stop 'Injection times must be nondecreasing'
    end do
    if(any(mixture_density_models(nc+1:)/='constant')) error stop 'Unused mixture property entries'
    if(any(mixture_cp_models(nc+1:)/='constant')) error stop 'Unused mixture property entries'
    if(any(mixture_vapor_pressure_models(nc+1:)/='clausius')) error stop 'Unused mixture property entries'
    if(any(mixture_boiling_models(nc+1:)/='reject')) error stop 'Unused mixture property entries'
    if(any(mixture_critical_temperatures(nc+1:)/=0)) error stop 'Unused mixture property entries'
    if(any(mixture_normal_boiling_temperatures(nc+1:)/=0)) error stop 'Unused mixture property entries'
    if(any(mixture_molecular_weights(nc+1:)/=0)) error stop 'Unused mixture property entries'
    if(any(mixture_density_coefficients(:,nc+1:)/=0)) error stop 'Unused mixture property entries'
    if(any(mixture_density_temperature_ranges(:,nc+1:)/=0)) error stop 'Unused mixture property entries'
    if(any(mixture_cp_coefficients(:,nc+1:)/=0)) error stop 'Unused mixture property entries'
    if(any(mixture_cp_temperature_ranges(:,nc+1:)/=0)) error stop 'Unused mixture property entries'
    if(any(mixture_antoine_coefficients(:,nc+1:)/=0)) error stop 'Unused mixture property entries'
    if(any(mixture_antoine_temperature_ranges(:,nc+1:)/=0)) error stop 'Unused mixture property entries'
    allocate(liquids(nc))
    do c=1,nc
      call initialize_spray_liquid(species,vapor_names(c),densities(c),heat_capacities(c),reference_temperatures(c), &
        reference_pressures(c),latent_heats(c),liquids(c),ok,mixture_vapor_pressure_models(c), &
        mixture_antoine_coefficients(:,c),mixture_antoine_temperature_ranges(:,c))
      if(.not.ok) error stop 'Invalid liquid properties'
      liquids(c)%density_model=mixture_density_models(c)
      liquids(c)%density_coefficients=mixture_density_coefficients(:,c)
      liquids(c)%density_temperature_range=mixture_density_temperature_ranges(:,c)
      liquids(c)%heat_capacity_model=mixture_cp_models(c)
      liquids(c)%heat_capacity_coefficients=mixture_cp_coefficients(:,c)
      liquids(c)%heat_capacity_temperature_range=mixture_cp_temperature_ranges(:,c)
      liquids(c)%boiling_model=mixture_boiling_models(c)
      liquids(c)%critical_temperature=mixture_critical_temperatures(c)
      liquids(c)%normal_boiling_temperature=mixture_normal_boiling_temperatures(c)
      liquids(c)%molecular_weight=mixture_molecular_weights(c)
      if(.not.valid_spray_liquid(liquids(c))) error stop 'Invalid mixture component property model'
    end do
    if(port_count<1.or.port_count>16) error stop 'Invalid injection port count'
    if(any(event_ports(:injection_events)<1).or.any(event_ports(:injection_events)>port_count)) &
      error stop 'Invalid injection event port'
    allocate(injection_ports(port_count))
    do port=1,port_count
      if(port_numbers(port)==-1) port_numbers(port)=number
      if(all(port_centers(:,port)==huge(1._dp))) port_centers(:,port)=center
      if(all(port_axes(:,port)==huge(1._dp))) port_axes(:,port)=axis
      if(port_masses(port)==-1) port_masses(port)=total_mass
      if(port_speeds(port)==-1) port_speeds(port)=speed
      if(port_angles(port)==-1) port_angles(port)=half_angle
      if(port_diameters(port)==-1) port_diameters(port)=diameter
      if(port_temperatures(port)==-1) port_temperatures(port)=liquid_temperature
      if(all(port_fractions(:,port)==-1)) port_fractions(:,port)=fractions
      if(port_numbers(port)<1.or.port_numbers(port)>min(125000,particle_options%maximum_parcels)) &
        error stop 'Invalid port population'
      if(.not.all(ieee_is_finite(port_fractions(:,port))).or.any(port_fractions(:,port)<0)) &
        error stop 'Invalid port fractions'
      if(abs(sum(port_fractions(:nc,port))-1)>1e-12_dp.or.any(port_fractions(nc+1:,port)/=0)) &
        error stop 'Invalid port fraction sum'
      call liquid_mixture_properties(liquids,port_fractions(:nc,port),port_temperatures(port),rho,cp,e,ok)
      if(.not.ok) error stop 'Invalid injected port mixture'
      effective=spray_liquid(vapor_name=liquids(1)%vapor_name,vapor_index=liquids(1)%vapor_index, &
        density=rho,heat_capacity=cp,reference_temperature=port_temperatures(port), &
        reference_pressure=liquids(1)%reference_pressure,latent_heat=liquids(1)%latent_heat,reference_energy=e)
      template_speed=port_speeds(port)
      if(port_velocity_from_rate(port)) then
        if(.not.ieee_is_finite(port_discharge_coefficients(port))) error stop 'Nonfinite discharge coefficient'
        if(port_discharge_coefficients(port)<1e-6_dp.or.port_discharge_coefficients(port)>1) &
          error stop 'Invalid discharge coefficient'
        if(.not.ieee_is_finite(port_radii(port))) error stop 'Nonfinite rate-driven outlet radius'
        if(port_radii(port)<1e-12_dp) error stop 'Rate-driven outlet requires positive radius'
        if(port_profile_points(port)<2) error stop 'Rate-driven outlet requires a mass-flow profile'
        template_speed=1
      else if(port_discharge_coefficients(port)/=1) then
        error stop 'Discharge coefficient requires rate-driven velocity'
      end if
      call inject_cone_parcels(effective,port_masses(port),port_numbers(port),1_int64, &
        port_centers(:,port),port_axes(:,port),template_speed,port_angles(port),port_diameters(port), &
        port_temperatures(port),injection_ports(port)%parcels,ok, &
        diameter_shape=port_sizes(port)%shape,diameter_min=port_sizes(port)%minimum, &
        diameter_max=port_sizes(port)%maximum,diameter_weighting=port_sizes(port)%weighting, &
        diameter_distribution=port_sizes(port)%distribution,diameter_std=port_sizes(port)%std, &
        nozzle_radius=port_radii(port),nozzle_model=port_nozzles(port))
      if(.not.ok) error stop 'Invalid injection port template'
      injection_ports(port)%fractions=port_fractions(:nc,port)
      injection_ports(port)%parcels%material_id=port
    end do
    event_factors=injection_scales;event_velocity_scales=1
    do port=1,port_count
      q=port_profile_points(port)
      if(q/=0.and.(q<2.or.q>64)) error stop 'Invalid port profile size'
      if(q==0) cycle
      call integrate_spray_injection_profile(port_profile_times(:q,port),port_profile_rates(:q,port), &
        port_profile_times(1,port),port_profile_times(q,port),profile_mass,ok)
      if(.not.ok) error stop 'Invalid port rate profile'
    end do
    do c=1,injection_events
      port=event_ports(c);q=port_profile_points(port)
      if(q==0) cycle
      call integrate_spray_injection_profile(port_profile_times(:q,port),port_profile_rates(:q,port), &
        injection_times(c),injection_end_times(c),profile_mass,ok)
      if(.not.ok) error stop 'Invalid injection integration interval'
      profile_mass=profile_mass*injection_scales(c)
      if(profile_mass>1e6_dp.or.profile_mass>1e6_dp*port_masses(port)) error stop 'Integrated injection exceeds bounds'
      event_factors(c)=profile_mass/port_masses(port)
      if(port_velocity_from_rate(port)) then
        event_velocity_scales(c)=0
        if(profile_mass>0) then
          call liquid_mixture_properties(liquids,port_fractions(:nc,port),port_temperatures(port),rho,cp,e,ok)
          if(.not.ok) error stop 'Invalid rate-driven mixture density'
          denominator=(injection_end_times(c)-injection_times(c))*rho*acos(-1._dp)* &
            port_radii(port)**2*port_discharge_coefficients(port)
          if(denominator<=0) error stop 'Invalid rate-driven injection interval'
          if(profile_mass>1e7_dp*denominator) error stop 'Rate-driven injection speed exceeds bounds'
          event_velocity_scales(c)=profile_mass/denominator
        end if
      end if
    end do
    allocate(drops(0),component_mass(nc,0),films(0),film_mass(nc,0))
    if(checkpoint_interval_steps>0.or.len_trim(restart_file)>0) then
      open(newunit=stream,status='scratch',form='formatted',iostat=code)
      if(code/=0) error stop 'Cannot stage spray context'
      write(stream,nml=eb_spray,iostat=code)
      if(code==0) write(stream,'(es27.18e3)',iostat=code) nodal_phi
      if(code/=0) error stop 'Cannot serialize spray context'
      rewind(stream)
      spray_context=mixture_evolution_id//'|'
      do
        read(stream,'(a)',iostat=code) context_line
        if(code==iostat_end) exit
        if(code/=0) error stop 'Cannot read spray context'
        if(len(spray_context)+len_trim(context_line)+1>130000) error stop 'Spray context capacity exceeded'
        spray_context=spray_context//trim(context_line)//'|'
      end do
      close(stream)
    end if
  end subroutine

  subroutine inject_due_events()
    type(spray_parcel) :: injected_parcel
    type(spray_parcel), allocatable :: candidate(:)
    real(dp), allocatable :: masses(:,:)
    real(dp) :: change(nv),totals(5),rho,cp,e,scale
    type(spray_liquid) :: effective
    integer :: q,c,number,old_count,local_count,global_count
    do while(next_injection<=injection_events)
      if(injection_times(next_injection)>record%time) exit
      scale=event_factors(next_injection)
      if(scale>0) then
        associate(injection_template=>injection_ports(event_ports(next_injection))%parcels, &
          injection_fractions=>injection_ports(event_ports(next_injection))%fractions)
        number=size(injection_template);local_count=size(drops)
        call MPI_Allreduce(local_count,global_count,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpi_ierr)
        if(mpi_ierr/=MPI_SUCCESS) error stop 'Injection population reduction failed'
        if(global_count>min(125000,particle_options%maximum_parcels)-number) error stop 'Injection capacity exceeded'
        if(last_particle_id>huge(last_particle_id)-number) error stop 'Particle ID overflow'
        candidate=drops;masses=component_mass;change=0
        ! The replicated template defines the same injected inventory on every rank.
        do q=1,number
          call liquid_mixture_properties(liquids,injection_fractions,injection_template(q)%temperature,rho,cp,e,ok)
          if(.not.ok) error stop 'Invalid injection properties'
          effective=spray_liquid(vapor_name=liquids(1)%vapor_name,vapor_index=liquids(1)%vapor_index, &
            density=rho,heat_capacity=cp,reference_temperature=injection_template(q)%temperature, &
            reference_pressure=liquids(1)%reference_pressure,latent_heat=liquids(1)%latent_heat,reference_energy=e)
          injected_parcel=injection_template(q)
          injected_parcel%velocity=event_velocity_scales(next_injection)*injected_parcel%velocity
          totals=scale*parcel_totals(injected_parcel,effective)
          change(:5)=change(:5)+totals
          do c=1,size(liquids)
            change(5+liquids(c)%vapor_index)=change(5+liquids(c)%vapor_index)+totals(1)*injection_fractions(c)
          end do
        end do
        if(mpi_rank==0) then
          old_count=size(drops)
          deallocate(candidate,masses)
          allocate(candidate(old_count+number),masses(size(liquids),old_count+number))
          candidate(:old_count)=drops;masses(:,:old_count)=component_mass
          do q=1,number
            candidate(old_count+q)=injection_template(q)
            candidate(old_count+q)%velocity=event_velocity_scales(next_injection)*injection_template(q)%velocity
            candidate(old_count+q)%id=last_particle_id+q
            candidate(old_count+q)%multiplicity=scale*injection_template(q)%multiplicity
            masses(:,old_count+q)=injection_fractions*candidate(old_count+q)%mass
          end do
        end if
        call migrate_eb_tree_spray_parcels(MPI_COMM_WORLD,topology,mpi_owners,nodal_phi,candidate, &
          domain_lower,(domain_upper-domain_lower)/real(root_shape,dp),particle_options%maximum_parcels,ok,masses)
        if(.not.ok) error stop 'Scheduled EB injection failed'
        call move_alloc(candidate,drops)
        call move_alloc(masses,component_mass)
        last_particle_id=last_particle_id+number;injected=injected+change
        end associate
      end if
      next_injection=next_injection+1
    end do
  end subroutine

  subroutine add_particle_inventory(inventory,particles,components)
    real(dp), intent(inout) :: inventory(:)
    type(spray_parcel), intent(in) :: particles(:)
    real(dp), intent(in) :: components(:,:)
    real(dp) :: local(size(inventory)),global(size(inventory)),rho,cp,e,amount
    integer :: q,c
    local=0
    do q=1,size(particles)
      if(particles(q)%mass<=0) cycle
      call liquid_mixture_properties(liquids,components(:,q)/particles(q)%mass,particles(q)%temperature,rho,cp,e,ok)
      if(.not.ok) error stop 'Invalid particle inventory'
      amount=particles(q)%mass*particles(q)%multiplicity
      local(1)=local(1)+amount;local(2:4)=local(2:4)+amount*particles(q)%velocity
      local(5)=local(5)+amount*(e+.5_dp*sum(particles(q)%velocity**2))
      do c=1,size(liquids)
        local(5+liquids(c)%vapor_index)=local(5+liquids(c)%vapor_index)+components(c,q)*particles(q)%multiplicity
      end do
    end do
    call MPI_Allreduce(local,global,size(local),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpi_ierr)
    if(mpi_ierr/=MPI_SUCCESS) error stop 'Particle inventory reduction failed'
    inventory=inventory+global
  end subroutine

  subroutine add_film_inventory(inventory,patches,components)
    real(dp), intent(inout) :: inventory(:)
    type(spray_wall_film), intent(in) :: patches(:)
    real(dp), intent(in) :: components(:,:)
    real(dp) :: local(size(inventory)),rho,cp,e
    integer :: q,c
    local=0
    do q=1,size(patches)
      call liquid_mixture_properties(liquids,components(:,q)/patches(q)%mass,patches(q)%temperature,rho,cp,e,ok)
      if(.not.ok) error stop 'Invalid film inventory'
      local(1)=local(1)+patches(q)%mass;local(5)=local(5)+patches(q)%mass*e
      do c=1,size(liquids)
        local(5+liquids(c)%vapor_index)=local(5+liquids(c)%vapor_index)+components(c,q)
      end do
    end do
    call MPI_Allreduce(MPI_IN_PLACE,local,size(local),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpi_ierr)
    if(mpi_ierr/=MPI_SUCCESS) error stop 'Film inventory reduction failed'
    inventory=inventory+local
  end subroutine

  subroutine spray_checkpoint(writing)
    logical, intent(in) :: writing
    type(amr_tree_state_3d) :: candidate,ct
    type(amr_patch_tree_3d) :: checkpoint_tree
    type(eb_geometry_3d), allocatable :: checkpoint_geometry(:)
    integer, allocatable :: checkpoint_owners(:)
    integer :: p,level
    type(eb_tree_checkpoint_3d) :: gas_record,coupled
    type(spray_parcel), allocatable :: particles(:)
    real(dp), allocatable :: components(:,:),patch_mass(:,:)
    type(spray_wall_film), allocatable :: patches(:)
    real(dp) :: source(nv),expected_source(nv),inventory(nv),saved_time,saved_force(nv)
    integer :: stream,code,closed,saved_steps,saved_event,expected_event,date(8)
    integer(int64) :: watermark,expected_mark,tick,expected_ids,expected_splits
    character(len=256) :: token,step_label
    character(len=:), allocatable :: path,bound_context,line
    logical :: valid,agreed,exists
    allocate(character(len=131072)::line)
    candidate=mpi_state;ct=mpi_temperature;coupled=record;gas_record=record
    checkpoint_tree=topology;checkpoint_geometry=mpi_geometry;checkpoint_owners=mpi_owners
    particles=drops;components=component_mass;patches=films;patch_mass=film_mass
    saved_time=record%time;saved_steps=record%steps;saved_event=next_injection
    watermark=last_particle_id;source=injected
    if(writing) then
      write(step_label,'(i0)') record%steps
      path=trim(checkpoint_file)//'.step'//trim(step_label)
      valid=.true.;token=''
      if(mpi_rank==0) then
        inquire(file=path,exist=exists);valid=.not.exists
        inquire(file=path//'.gas',exist=exists);valid=valid.and..not.exists
        inquire(file=path//'.particles',exist=exists);valid=valid.and..not.exists
        inquire(file=path//'.particles.pending',exist=exists);valid=valid.and..not.exists
        code=1
        if(valid) open(newunit=stream,file=path//'.pending',status='new',action='write',iostat=code)
        valid=code==0
        if(valid) then
          call date_and_time(values=date);call system_clock(tick)
          write(token,'(8(i0,":"),i0)') date,tick
          write(stream,'(a)',iostat=code) 'PELEF_MPI_EB_SPRAY_V1',trim(token),'END_EB_SPRAY_MANIFEST'
          close(stream,iostat=closed);valid=code==0.and.closed==0
        end if
      end if
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,mpi_ierr)
      if(mpi_ierr/=MPI_SUCCESS.or..not.agreed) error stop 'Cannot reserve EB spray checkpoint'
      call MPI_Bcast(token,len(token),MPI_CHARACTER,0,MPI_COMM_WORLD,mpi_ierr)
      if(mpi_ierr/=MPI_SUCCESS) error stop 'Checkpoint generation broadcast failed'
      bound_context=spray_context//'GENERATION:'//trim(token)
      open(newunit=stream,status='scratch',form='formatted',iostat=code)
      if(code/=0) error stop 'Cannot stage coupled checkpoint'
      call envelope(stream,.false.,header_only=.true.)
      write(stream,'(a)',iostat=code) bound_context
      if(code==0) then
        if(particle_options%maximum_multiplicity>0) then
          write(stream,*,iostat=code) next_injection,last_particle_id,split_created_total
        else
          write(stream,*,iostat=code) next_injection,last_particle_id
        end if
      end if
      if(code==0) write(stream,'(*(es27.18e3,1x))',iostat=code) &
        record%initial,record%scale,record%boundary,record%wall,record%chemical,injected
      if(any(acceleration/=0).and.code==0) write(stream,'(*(es27.18e3,1x))',iostat=code) force_total
      if(code/=0) error stop 'Cannot write coupled checkpoint prefix'
      ! The gas archive is a snapshot with its own zero-exchange reference.
      ! The complete coupled history is retained separately in the prefix above.
      call integrate_mpi_eb_tree_3d(checkpoint_tree,checkpoint_owners,checkpoint_geometry,candidate,inventory,MPI_COMM_WORLD,valid)
      if(.not.valid) error stop 'Cannot integrate gas snapshot'
      gas_record%initial=inventory;gas_record%scale=abs(inventory)
      gas_record%boundary=0;gas_record%wall=0;gas_record%chemical=0
      call write_mpi_eb_tree_checkpoint_3d(path//'.gas',stream,species,reactions,transport,config,topology, &
        mpi_owners,mpi_geometry,candidate,ct,gas_record,MPI_COMM_WORLD,valid, &
        bundle_identity,checkpoint_integrator,base_mole_fractions,element_names,element_atoms)
      close(stream,iostat=closed)
      if(.not.valid.or.closed/=0) error stop 'Gas snapshot write failed'
    else
      path=trim(restart_file)
      open(newunit=stream,file=path,status='old',action='read',iostat=code)
      if(code/=0) error stop 'Cannot open EB spray manifest'
      read(stream,'(a)',iostat=code) line
      if(code/=0.or.line/='PELEF_MPI_EB_SPRAY_V1') error stop 'Invalid EB spray manifest'
      read(stream,'(a)',iostat=code) token
      if(code/=0.or.len_trim(token)==0) error stop 'Missing checkpoint generation'
      read(stream,'(a)',iostat=code) line
      if(code/=0.or.line/='END_EB_SPRAY_MANIFEST') error stop 'Invalid manifest footer'
      read(stream,'(a)',iostat=code) line
      if(code/=iostat_end) error stop 'Trailing manifest data'
      close(stream)
      bound_context=spray_context//'GENERATION:'//trim(token)
      open(newunit=stream,file=path//'.gas',status='old',action='read',access='stream',form='formatted',iostat=code)
      if(code/=0) error stop 'Cannot open EB spray gas snapshot'
      restart_layout_stage=0
      call envelope(stream,.true.,header_only=.true.)
      checkpoint_tree=layouts(restart_layout_stage)
      call partition_amr_tree_3d(checkpoint_tree,mpi_size,0,checkpoint_owners,valid)
      if(.not.valid) error stop 'Restart partition failed'
      level=max(topology%level_count(),checkpoint_tree%level_count())
      call remap_mpi_eb_tree_3d(species,topology,checkpoint_tree,mpi_owners,checkpoint_owners,geometry_levels(level), &
        checkpoint_geometry,candidate,ct,MPI_COMM_WORLD,valid)
      if(.not.valid) error stop 'Cannot prepare restart hierarchy'

      read(stream,'(a)',iostat=code) line
      if(code/=0.or.line/=bound_context) error stop 'EB spray context mismatch'
      expected_splits=0
      if(particle_options%maximum_multiplicity>0) then
        read(stream,*,iostat=code) expected_event,expected_mark,expected_splits
      else
        read(stream,*,iostat=code) expected_event,expected_mark
      end if
      if(code/=0) error stop 'Missing injection metadata'
      read(stream,*,iostat=code) coupled%initial,coupled%scale,coupled%boundary,coupled%wall,coupled%chemical,expected_source
      if(code/=0) error stop 'Missing coupled history'
      if(.not.all(ieee_is_finite([coupled%initial,coupled%scale,coupled%boundary,coupled%wall, &
        coupled%chemical,expected_source])).or.any(coupled%scale<0)) error stop 'Invalid coupled history'
      saved_force=0
      if(any(acceleration/=0)) then
        read(stream,*,iostat=code) saved_force
        if(code/=0) error stop 'Missing external force history'
        if(.not.all(ieee_is_finite(saved_force))) error stop 'Invalid external force history'
      end if
      call read_mpi_eb_tree_checkpoint_3d(path//'.gas',stream,species,reactions,transport,config, &
        checkpoint_tree,checkpoint_owners, &
        checkpoint_geometry,candidate,ct,gas_record,MPI_COMM_WORLD,valid, &
        bundle_identity,checkpoint_integrator,base_mole_fractions,element_names,element_atoms)
      close(stream,iostat=closed)
      if(.not.valid.or.closed/=0) error stop 'Gas snapshot restore failed'
    end if
    if(wall_impacts) then
    call mpi_eb_particle_checkpoint(MPI_COMM_WORLD,writing,path//'.particles',bound_context,checkpoint_tree,checkpoint_owners, &
      nodal_phi,domain_lower,(domain_upper-domain_lower)/real(root_shape,dp), &
      min(125000,particle_options%maximum_parcels),particles,components,saved_time,saved_steps,saved_event, &
      watermark,source,valid,patches,patch_mass,walls%maximum_films)
    else
    call mpi_eb_particle_checkpoint(MPI_COMM_WORLD,writing,path//'.particles',bound_context,checkpoint_tree,checkpoint_owners, &
      nodal_phi,domain_lower,(domain_upper-domain_lower)/real(root_shape,dp), &
      min(125000,particle_options%maximum_parcels),particles,components,saved_time,saved_steps,saved_event, &
      watermark,source,valid)
    end if
    if(.not.valid) error stop 'Particle snapshot failed'
    if(writing) then
      if(mpi_rank==0) valid=publish_link(path//'.pending'//c_null_char,path//c_null_char)==0
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,mpi_ierr)
      if(mpi_ierr/=MPI_SUCCESS.or..not.agreed) error stop 'Checkpoint manifest publication failed'
    else
      if(saved_time/=gas_record%time.or.saved_steps/=gas_record%steps.or.saved_event/=expected_event.or. &
        watermark/=expected_mark.or.any(source/=expected_source)) error stop 'Checkpoint payload mismatch'
      if(restart_layout_stage/=count(regrid_times(:regrid_count)<=saved_time)) error stop 'Restart regrid stage mismatch'
      if(saved_time>final_time.or.saved_steps>maximum_steps) error stop 'Restart beyond requested limits'
      if(saved_event<1.or.saved_event>injection_events+1) error stop 'Invalid restart event'
      if(any(injection_times(:saved_event-1)>saved_time)) error stop 'Processed future injection'
      if(saved_event<=injection_events) then
        if(injection_times(saved_event)<saved_time) error stop 'Missed restart injection'
      end if
      expected_ids=0
      do p=1,saved_event-1
        if(event_factors(p)>0) expected_ids=expected_ids+size(injection_ports(event_ports(p))%parcels)
      end do
      if(watermark<expected_ids) error stop 'Invalid restart ID watermark'
      if(expected_splits<0.or.expected_splits>watermark-expected_ids) error stop 'Invalid restart split ID count'
      if(wall_impacts.or.particle_options%khrt_children) then
        if(modulo(watermark-expected_ids-expected_splits,4_int64)/=0) error stop 'Invalid restart child ID watermark'
      else
        if(watermark-expected_ids/=expected_splits) error stop 'Invalid restart ID watermark'
      end if
      do p=1,size(patches)
        if(patches(p)%material_id>port_count) error stop 'Invalid restart film port'
      end do
      do p=1,size(particles)
        if(particles(p)%material_id>port_count) error stop 'Invalid restart particle port'
      end do
      coupled%time=saved_time;coupled%steps=saved_steps
      coupled%minimum_dt=gas_record%minimum_dt;coupled%minimum_theta=gas_record%minimum_theta
      call integrate_mpi_eb_tree_3d(checkpoint_tree,checkpoint_owners,checkpoint_geometry,candidate,inventory,MPI_COMM_WORLD,valid)
      if(.not.valid) error stop 'Cannot integrate restored gas'
      call add_particle_inventory(inventory,particles,components)
      call add_film_inventory(inventory,patches,patch_mass)
      if(any(abs(inventory-coupled%initial-coupled%boundary-coupled%wall-coupled%chemical-source-saved_force)> &
        5e-10_dp*max(1._dp,coupled%scale,abs(inventory), &
        abs(coupled%boundary)+abs(coupled%wall)+abs(coupled%chemical)+abs(source)+abs(saved_force)))) &
        error stop 'Restored budget failed'
      mpi_state=candidate;mpi_temperature=ct;record=coupled
      topology=checkpoint_tree;mpi_owners=checkpoint_owners;mpi_geometry=checkpoint_geometry
      regrid_stage=restart_layout_stage;tree%topology=topology
      if(allocated(tree%nodes)) deallocate(tree%nodes)
      allocate(tree%nodes(size(topology%nodes)))
      do p=1,size(tree%nodes)
        tree%nodes(p)%geometry=mpi_geometry(p)
      end do
      call move_alloc(particles,drops);call move_alloc(components,component_mass)
      call move_alloc(patches,films);call move_alloc(patch_mass,film_mass)
      force_total=saved_force
      injected=source;next_injection=saved_event;last_particle_id=watermark;split_created_total=expected_splits
    end if
  end subroutine

  subroutine output_particles()
    integer :: stream,code,closed,q,c
    character(len=32) :: rank_label,label
    character(len=:), allocatable :: path
    write(rank_label,'(i0)') mpi_rank
    path=trim(output_file)//'.particles.rank'//trim(rank_label)//'.csv'
    open(newunit=stream,file=path,status='new',action='write',iostat=code)
    if(code/=0) error stop 'Cannot create particle CSV (existing files are protected)'
    write(stream,'(a)',advance='no') 'id,material_id,x,y,z,u,v,w,mass,temperature,multiplicity'
    do c=1,size(liquids)
      write(label,'(i0)') c
      write(stream,'(a)',advance='no') ',component_mass_'//trim(label)
    end do
    write(stream,*)
    do q=1,size(drops)
      write(stream,'(i0,",",i0,",",*(es25.16,:,","))',iostat=code) &
        drops(q)%id,drops(q)%material_id,drops(q)%position,drops(q)%velocity, &
        drops(q)%mass,drops(q)%temperature,drops(q)%multiplicity,component_mass(:,q)
      if(code/=0) error stop 'Cannot write particle CSV'
    end do
    close(stream,iostat=closed)
    if(closed/=0) error stop 'Cannot close particle CSV'
    path=trim(output_file)//'.films.rank'//trim(rank_label)//'.csv'
    open(newunit=stream,file=path,status='new',action='write',iostat=code)
    if(code/=0) error stop 'Cannot create film CSV (existing files are protected)'
    write(stream,'(a)',advance='no') 'material_id,x,y,z,nx,ny,nz,mass,temperature,diameter,height'
    do c=1,size(liquids)
      write(label,'(i0)') c
      write(stream,'(a)',advance='no') ',component_mass_'//trim(label)
    end do
    write(stream,*)
    do q=1,size(films)
      write(stream,'(i0,",",*(es25.16,:,","))',iostat=code) films(q)%material_id, &
        films(q)%position,films(q)%normal,films(q)%mass,films(q)%temperature, &
        films(q)%diameter,films(q)%height,film_mass(:,q)
      if(code/=0) error stop 'Cannot write film CSV'
    end do
    close(stream,iostat=closed)
    if(closed/=0) error stop 'Cannot close film CSV'
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
    if(allocated(active_faces).and..not.allocated(inlet_state)) then
      if(mpi_rank==mpi_owners(1)) then
        call prepare_prescribed_inlet(seed_u(:,1,1,1),seed_t(1,1,1),valid)
      else
        allocate(inlet_state(nv,6),active_temperature(6))
      end if
      call MPI_Allreduce(valid,all_valid,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
      if(ierr/=MPI_SUCCESS.or..not.all_valid) error stop 'Invalid prescribed inlet EOS'
      call MPI_Bcast(inlet_state,nv*6,MPI_DOUBLE_PRECISION,mpi_owners(1),MPI_COMM_WORLD,ierr)
      if(ierr/=MPI_SUCCESS) error stop 'Cannot broadcast inlet states'
      call MPI_Bcast(active_temperature,6,MPI_DOUBLE_PRECISION,mpi_owners(1),MPI_COMM_WORLD,ierr)
      if(ierr/=MPI_SUCCESS) error stop 'Cannot broadcast inlet temperatures'
    end if
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
#ifdef PELEF_MPI_EB_SPRAY
      call remap_mpi_eb_spray(MPI_COMM_WORLD,species,topology,mpi_owners,mpi_geometry,mpi_state,mpi_temperature, &
        drops,component_mass,layouts(stage),next_owners,geometry_levels(level),nodal_phi,domain_lower, &
        (domain_upper-domain_lower)/real(root_shape,dp),particle_options%maximum_parcels,valid,films,film_mass,walls%maximum_films)
#else
      call remap_mpi_eb_tree_3d(species,topology,layouts(stage),mpi_owners,next_owners,geometry_levels(level), &
        mpi_geometry,mpi_state,mpi_temperature,MPI_COMM_WORLD,valid)
#endif
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
  subroutine refresh_ctu_geometry()
    type(eb_geometry_3d), allocatable :: physical(:)
    integer, allocatable :: descriptors(:,:)
    integer :: p
    logical :: valid
    if(hydro_method/='plm_ctu') return
    allocate(descriptors(12,size(topology%nodes)))
    do p=1,size(topology%nodes)
      descriptors(:,p)=descriptor_for(topology,p)
    end do
    if(allocated(cached_descriptors)) then
      if(all(shape(cached_descriptors)==shape(descriptors))) then
        if(all(cached_descriptors==descriptors)) return
      end if
    end if
    allocate(physical(size(topology%nodes)))
    do p=1,size(topology%nodes)
#ifdef PELEF_MPI_EB_TREE
      physical(p)=mpi_geometry(p)
#else
      physical(p)=tree%nodes(p)%geometry
#endif
    end do
    call build_eb_tree_extended_geometry_3d(topology,physical,extended_levels(topology%level_count()), &
      extended_geometry,valid)
    if(.not.valid) error stop 'Cannot build CTU patch geometry cache'
    call move_alloc(descriptors,cached_descriptors)
  end subroutine refresh_ctu_geometry
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
#ifdef PELEF_MPI_EB_SPRAY
    if(reading) restart_layout_stage=stored_stage
#else
    if (reading) call activate(stored_stage,.false.)
#endif
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
#ifdef PELEF_MPI_EB_SPRAY
    if(reading) restart_layout_stage=stored_stage
#else
    if (reading) call activate(stored_stage,.false.)
#endif
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
    call add_particle_inventory(current,drops,component_mass)
    call add_film_inventory(current,films,film_mass)
#endif
    error=maxval(abs(current-record%initial-record%boundary-record%wall-record%chemical-injected-force_total)/ &
      max(1._dp,record%scale,abs(current),abs(record%boundary)+abs(record%wall)+ &
      abs(record%chemical)+abs(injected)+abs(force_total)))
    if (.not.ieee_is_finite(error).or.error>5.e-10_dp) error stop 'Composite budget failed'
  end subroutine
  subroutine prepare_inlet_compositions()
    real(dp) :: mole(size(species))
    integer :: face
    allocate(inlet_mass_fractions(size(species),6));inlet_mass_fractions=0._dp
    do face=1,6
      if(inflow_composition_count(face)==0) then
        if(any(inflow_composition_species(:,face)/='').or. &
          any(inflow_composition_mole_fractions(:,face)/=0._dp)) &
          error stop 'Inlet composition fields require a positive count'
        cycle
      end if
      if(boundary_faces(face)/='inflow') error stop 'Inlet composition requires an inflow face'
      call resolve_selected_composition('Tree inlet',inflow_composition_count(face), &
        inflow_composition_species(:,face),inflow_composition_mole_fractions(:,face),species,mole,ok,message)
      if(.not.ok) then
        write(*,'(a)') trim(message)
        error stop 'Invalid inlet composition'
      end if
      call mass_fractions_from_mole_fractions(species,mole,inlet_mass_fractions(:,face),ok)
      if(.not.ok) error stop 'Invalid inlet mass fractions'
    end do
  end subroutine
  subroutine prepare_prescribed_inlet(seed,seed_temperature,valid)
    real(dp), intent(in) :: seed(:),seed_temperature
    logical, intent(out) :: valid
    real(dp) :: seed_primitive(size(seed))
    allocate(inlet_state(size(seed),6),inlet_primitive(size(seed)),active_temperature(6))
    call reactive_conserved_to_primitive(species,seed,seed_temperature,inlet_primitive,inlet_t,inlet_sound,valid)
    if(.not.valid) return
    seed_primitive=inlet_primitive
    do inlet_face=1,6
      inlet_state(:,inlet_face)=seed;active_temperature(inlet_face)=seed_temperature
      if(active_faces(inlet_face)/='inflow') cycle
      inlet_primitive=seed_primitive
      if(inflow_composition_count(inlet_face)>0) inlet_primitive(6:)=inlet_mass_fractions(:,inlet_face)
      inlet_primitive(1)=inflow_density(inlet_face);inlet_primitive(2:4)=inflow_velocity(:,inlet_face)
      inlet_primitive(5)=inflow_pressure(inlet_face)
      call reactive_primitive_to_conserved(species,inlet_primitive,inlet_state(:,inlet_face), &
        active_temperature(inlet_face),inlet_sound,valid)
      if(.not.valid) return
    end do
  end subroutine
  subroutine envelope(stream,reading,header_only)
    integer, intent(in) :: stream
    logical, intent(in) :: reading
    logical, optional, intent(in) :: header_only
    character(len=64) :: magic,integrator
    logical :: automatic
    character(len=16) :: stored_faces(6)
    integer :: stored_counts(6),stored_nspecies
    real(dp) :: stored_composition(size(species),6)
    real(dp) :: stored_prescribed(30),stored_wall_temperature(6)
    real(dp) :: interval,stored_hydro_cfl,stored_transport_cfl,stored_inlet(18)
    if (reading) then
      read(stream,'(a)',iostat=ios) magic
      if (ios/=0) error stop 'Missing tree application header'
      if(magic=='PELEF_EB_TREE_ISOTHERMAL_V1') then
        if(.not.allocated(active_wall_temperature)) error stop 'Restart wall temperature mismatch'
        read(stream,*,iostat=ios) stored_wall_temperature
        if(ios/=0) error stop 'Missing restart wall temperatures'
        if(.not.all(ieee_is_finite(stored_wall_temperature))) error stop 'Invalid restart wall temperatures'
        if(any(stored_wall_temperature/=wall_temperature)) error stop 'Restart wall temperature mismatch'
        read(stream,'(a)',iostat=ios) magic
        if(ios/=0) error stop 'Missing isothermal application header'
      else if(allocated(active_wall_temperature)) then
        error stop 'Restart wall temperature mismatch'
      end if
      if(magic=='PELEF_EB_TREE_INLET_COMPOSITION_V1') then
        read(stream,*,iostat=ios) stored_nspecies,stored_counts
        if(ios/=0) error stop 'Missing inlet composition dimensions'
        if(stored_nspecies/=size(species).or.any(stored_counts/=inflow_composition_count)) &
          error stop 'Restart inlet composition mismatch'
        read(stream,*,iostat=ios) stored_composition
        if(ios/=0) error stop 'Missing inlet composition values'
        if(.not.all(ieee_is_finite(stored_composition))) error stop 'Invalid restart inlet composition'
        if(any(stored_composition/=inlet_mass_fractions)) error stop 'Restart inlet composition mismatch'
        read(stream,'(a)',iostat=ios) magic
        if(ios/=0) error stop 'Missing inlet composition application header'
      else if(any(inflow_composition_count/=0)) then
        error stop 'Restart inlet composition mismatch'
      end if
      if(magic=='PELEF_EB_TREE_PRESCRIBED_INLET_V1') then
        if(.not.allocated(active_faces)) error stop 'Restart prescribed inlet mismatch'
        read(stream,*,iostat=ios) stored_faces
        if(ios/=0) error stop 'Missing restart boundary labels'
        if(any(stored_faces/=boundary_faces)) error stop 'Restart boundary labels mismatch'
        read(stream,*,iostat=ios) stored_prescribed
        if(ios/=0) error stop 'Missing prescribed inlet controls'
        if(.not.all(ieee_is_finite(stored_prescribed))) error stop 'Invalid prescribed inlet restart controls'
        if(any(stored_prescribed/=[inflow_density,inflow_pressure,reshape(inflow_velocity,[18])])) &
          error stop 'Restart prescribed inlet mismatch'
        read(stream,'(a)',iostat=ios) magic
        if(ios/=0) error stop 'Missing prescribed inlet application header'
      else if(allocated(active_faces)) then
        error stop 'Restart prescribed inlet mismatch'
      end if
      if (magic=='PELEF_EB_TREE_MACH_INLET_V1') then
        read(stream,*,iostat=ios) stored_inlet
        if (ios/=0) error stop 'Missing inlet restart controls'
        if (.not.all(ieee_is_finite(stored_inlet))) error stop 'Invalid inlet restart controls'
        if (any(stored_inlet/=[inflow_mach,inflow_temperature,inflow_pressure])) &
          error stop 'Restart inlet mismatch'
        read(stream,'(a)',iostat=ios) magic
        if (ios/=0) error stop 'Missing inlet application header'
      else if (any(inflow_mach>=0)) then
        error stop 'Restart inlet mismatch'
      end if
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
      if(allocated(active_wall_temperature)) then
        write(stream,'(a)',iostat=ios) 'PELEF_EB_TREE_ISOTHERMAL_V1'
        if(ios==0) write(stream,'(6(es27.18e3,1x))',iostat=ios) wall_temperature
        if(ios/=0) error stop 'Cannot write wall temperature context'
      end if
      if(any(inflow_composition_count/=0)) then
        write(stream,'(a)',iostat=ios) 'PELEF_EB_TREE_INLET_COMPOSITION_V1'
        if(ios==0) write(stream,*,iostat=ios) size(species),inflow_composition_count
        if(ios==0) write(stream,'(*(es27.18e3,1x))',iostat=ios) inlet_mass_fractions
        if(ios/=0) error stop 'Cannot write inlet composition context'
      end if
      if(allocated(active_faces)) then
        write(stream,'(a)',iostat=ios) 'PELEF_EB_TREE_PRESCRIBED_INLET_V1'
        if(ios==0) write(stream,'(6(a,1x))',iostat=ios) boundary_faces
        if(ios==0) write(stream,'(30(es27.18e3,1x))',iostat=ios) &
          inflow_density,inflow_pressure,inflow_velocity
        if(ios/=0) error stop 'Cannot write prescribed inlet context'
      end if
      if (any(inflow_mach>=0)) then
        write(stream,'(a)',iostat=ios) 'PELEF_EB_TREE_MACH_INLET_V1'
        if (ios/=0) error stop 'Cannot write inlet header'
        write(stream,'(18(es27.18e3,1x))',iostat=ios) inflow_mach,inflow_temperature,inflow_pressure
        if (ios/=0) error stop 'Cannot write inlet restart controls'
      end if
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
#ifdef PELEF_MPI_EB_SPRAY
    write(unit,'(a)') 'component,initial,final,boundary,wall,chemistry,injected,force,residual'
    do c=1,nv
      write(unit,'(i0,8(",",es26.17e3))') c,record%initial(c),current(c),record%boundary(c), &
        record%wall(c),record%chemical(c),injected(c),force_total(c), &
        current(c)-record%initial(c)-record%boundary(c)-record%wall(c)-record%chemical(c)-injected(c)-force_total(c)
    end do
    write(*,'(a,i0)') 'Processed injection events: ',next_injection-1
#else
    write(unit,'(a)') 'component,initial,final,boundary,wall,chemistry,residual'
    do c=1,nv
      write(unit,'(i0,6(",",es26.17e3))') c,record%initial(c),current(c),record%boundary(c), &
        record%wall(c),record%chemical(c),current(c)-record%initial(c)-record%boundary(c)-record%wall(c)-record%chemical(c)
    end do
#endif
    close(unit,iostat=ios)
    if (ios/=0) error stop 'Budget close failed'
  end subroutine
end program pelef_reactive_eb_tree_3d
