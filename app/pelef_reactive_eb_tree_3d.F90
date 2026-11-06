program pelef_reactive_eb_tree_3d
#ifdef PELEF_MPI_EB_TREE
  use mpi_f08
  use mpi_reactive_eos_policy_mod, only: agree_mpi_reactive_eos_1d
  use mpi_eb_wall_heat_3d_mod, only: compute_mpi_eb_wall_transport_3d
  use mpi_eb_tree_initialize_3d_mod, only: initialize_mpi_eb_tree_3d
  use mpi_eb_tree_checkpoint_read_3d_mod, only: read_mpi_eb_tree_checkpoint_3d
  use mpi_eb_tree_checkpoint_3d_mod, only: write_mpi_eb_tree_checkpoint_3d
  use mpi_eb_tree_csv_3d_mod, only: write_mpi_eb_tree_csv_3d
  use amr_patch_tree_3d_mod, only: amr_tree_state_3d
  use mpi_amr_patch_tree_distribution_3d_mod, only: partition_amr_tree_3d,rebalance_mpi_amr_tree_3d, &
    rebalance_mpi_amr_tree_with_field_3d
  use mpi_amr_tree_work_mod, only: begin_tree_work,end_tree_work
  use mpi_eb_tree_tagging_3d_mod, only: tag_mpi_eb_tree_temperature_3d
  use mpi_eb_tree_remap_3d_mod, only: remap_mpi_eb_tree_3d
  use mpi_eb_tree_timestep_3d_mod, only: mpi_eb_tree_timestep_3d
  use mpi_eb_tree_full_3d_mod, only: initialize_mpi_eb_tree_reaction_source_3d, &
    advance_mpi_eb_tree_full_3d,advance_mpi_eb_tree_sdc_3d, &
    advance_mpi_eb_tree_sdc_subcycle_3d,mpi_eb_sdc_diffusion_source,mpi_eb_sdc_spatial_diffusion_source
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d
#endif
  use eb_mol_hydro_3d_mod, only: eb_mol_diagnostics_3d
  use gas_transport_mod, only: set_constant_transport_coefficients
  use reactive_3d_application_mod, only: read_constant_transport
  use reactive_transport_3d_mod, only: valid_constant_transport_coefficients
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
#ifdef PELEF_MPI_EB_LES
  use les_mpi_eb_sdc_mod, only: compute_mpi_regular_eb_les_sdc_3d,compute_mpi_regular_eb_les_spatial_sdc_3d
#endif
#if defined(PELEF_MPI_EB_LES) || defined(PELEF_SERIAL_EB_LES)
  use les_filter_3d_mod, only: les_filter_weights
#endif
#ifdef PELEF_SERIAL_EB_LES
  use les_eb_sdc_mod, only: compute_regular_eb_les_spatial_sdc_3d
#endif
#ifdef PELEF_MPI_EB_SPRAY
  use, intrinsic :: iso_c_binding, only: c_int,c_char,c_null_char
  use spray_mpi_eb_particle_checkpoint_mod, only: mpi_eb_particle_checkpoint
  use spray_injection_profile_mod, only: integrate_spray_injection_profile,read_spray_rate_profile, &
    spray_rate_profile,spray_event_schedule,read_spray_event_schedule
  use spray_mpi_amr_step_mod, only: compute_mpi_amr_spray_timestep
  use spray_mpi_eb_remap_mod, only: remap_mpi_eb_spray
  use spray_mpi_amr_balance_mod, only: rebalance_mpi_amr_spray
  use spray_parcel_mod
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options,valid_wall_options
  use spray_mixture_mod, only: apply_mixture_breakup_viscosity, liquid_mixture_properties,mixture_evolution_id
  use spray_mpi_amr_tree_mod, only: migrate_eb_tree_spray_parcels
  use spray_mpi_eb_reactive_mod, only: advance_mpi_eb_spray_reactive,advance_mpi_eb_spray_sdc, &
    advance_mpi_eb_spray_sdc_subcycle
#endif
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64,iostat_end,iostat_eor
  use precision_mod, only: dp
  use simulation_config_selected_reactor_mod, only: paths_alias => selected_reactor_paths_alias_lexically
  use rotating_frame_3d_mod, only: rotating_frame_3d,convert_rotating_energy_3d,rotating_potential_3d
  use dual_energy_mod, only: dual_energy_policy
  use constant_volume_reactor_mod, only: forced_reactor_integrator
#ifdef PELEF_SDC_CVODE
  use sundials_forced_reactor_mod, only: advance_forced_reactor_cvode,advance_forced_reactor_cvode_full
#endif
  use nasa7_thermo_mod, only: nasa7_species,configure_nasa7_eos_policy
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
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,initialize_amr_patch_tree_3d,amr_tree_integer_layout_3d,amr_tree_state_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config,build_configured_eb_geometry_3d
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use eb_reactive_patch_tree_3d_mod, only: initialize_eb_tree_reaction_source_3d, &
    eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    advance_reactive_eb_tree_full_3d,integrate_reactive_eb_tree_3d,reactive_eb_tree_timestep_3d, &
    build_eb_tree_extended_geometry_3d,advance_reactive_eb_tree_sdc_3d,advance_reactive_eb_tree_sdc_subcycle_3d, &
    eb_sdc_spatial_diffusion_source,fill_reactive_eb_tree_halo_3d, &
    remap_reactive_eb_tree_3d,tag_reactive_eb_tree_temperature_3d,compute_eb_tree_wall_transport_3d
  use eb_reactive_tree_checkpoint_3d_mod, only: eb_tree_checkpoint_3d,reactive_eb_tree_checkpoint_record_3d
  use selected_composition_mod, only: selected_composition_max_species,selected_composition_name_length, &
    resolve_selected_composition
#ifdef PELEF_SELECTED_SIMPLE_TRANSPORT
  use @PELEF_BUNDLE_MODULE_NAME@, only: load_selected_simple => @PELEF_BUNDLE_TRANSPORT_LOADER_NAME@_simple
  use gas_transport_mod, only: set_simple_transport_fits
#endif
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
  logical :: initial_reactions_pending=.false.
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
  real(dp) :: injection_times(1024)=0,injection_scales(1024)=1
  type(spray_event_schedule) :: event_schedule
  real(dp), allocatable :: event_factors(:),event_velocity_scales(:)
  type(spray_liquid), allocatable :: liquids(:)
  type(spray_options) :: particle_options
  type(spray_wall_options) :: walls
  type(spray_wall_film), allocatable :: films(:)
  real(dp), allocatable :: film_mass(:,:)
  real(dp) :: maximum_injection_speed,particle_dt
  real(dp) :: acceleration(3)=0
  logical :: wall_impacts=.false.,particle_outflow=.false.,model_rejected
  real(dp), allocatable :: nodal_phi(:,:,:),component_mass(:,:)
  character(len=:), allocatable :: spray_context
  character(len=1024) :: spray_input
  integer :: particle_steps,restart_layout_stage=0
#endif
  logical :: les_enabled=.false.,les_harmonic_mean=.false.
  real(dp) :: les_cs=.17_dp,les_prt=.7_dp,les_sct=.7_dp,les_ci=0,les_cfl=.1_dp
  real(dp) :: les_diffusivity,les_theta,les_spacing(3)
  character(len=32) :: les_model='smagorinsky',les_test_model='box'
  integer :: les_test_ratio=2,les_grid_ratio=1
  character(len=:), allocatable :: les_filter
  integer, allocatable :: les_ratio,les_grid
  logical, allocatable :: les_harmonic
  type(amr_tree_state_3d) :: les_rates
  type(eb_integrated_flux_3d), allocatable :: les_flux(:)
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
  logical :: automatic_timestep=.false.,subcycling=.false.,timestep_after_regrid=.false.
  real(dp) :: initial_timestep_factor=1._dp,maximum_timestep_growth=1.1_dp,previous_timestep=0._dp
  real(dp), allocatable :: constant_coefficients(:),hit_transport_parameters(:)
  character(len=16) :: time_integrator="split",sdc_reactor="native"
  integer :: sdc_iterations=1 ! Pinned SDC/MOL default; explicit iterations remain configurable.
  real(dp) :: external_forcing(3)=0
  type(amr_tree_state_3d), allocatable :: sdc_ir,sdc_heat,held_energy
  type(eb_integrated_flux_3d), allocatable :: held_flux(:)
  type(dual_energy_policy), allocatable :: reset_policy
  logical :: dual_energy_enabled=.false.,allow_small_energy=.true.,allow_negative_energy=.true.
  logical :: update_total_from_internal=.true.
  real(dp) :: dual_energy_eta2=1.e-4_dp,numerical_energy=0._dp
  real(dp) :: internal_energy_pressure_work=0._dp,step_internal_energy_pressure_work=0._dp
  real(dp) :: inflow_mach(6)=-1._dp,inflow_temperature(6)=0._dp,inflow_pressure(6)=0._dp
  real(dp), allocatable :: active_mach(:),active_temperature(:),active_pressure(:),inlet_state(:,:),inlet_primitive(:)
  real(dp) :: wall_temperature(6)=0,embedded_wall_temperature=0
  logical :: embedded_noslip=.false.,mol_mode=.false.,mol_lax_friedrichs=.false.
  logical :: rotating_frame_enabled=.false.
  integer :: rotation_axis=2
  real(dp) :: rotation_angular_speed=0._dp,rotation_axis_point(3)=0._dp
  type(rotating_frame_3d), allocatable :: frame
  real(dp), allocatable :: active_embedded_temperature
  type(eb_geometry_3d), allocatable :: heat_geometry(:)
  type(amr_tree_state_3d) :: heat_rates
  real(dp) :: heat_frequency
  real(dp), allocatable :: wall_exchange(:)
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
  integer :: state_redist_max_order=0
  real(dp) :: state_redist_target_volume_fraction=.5_dp
  character(len=16) :: riemann_solver='rusanov',hydro_method='pcm',redistribution='none'
#ifdef PELEF_SELECTED_EB_TREE
  character(len=16) :: chemistry_integrator=selected_integrator
#else
  character(len=16) :: chemistry_integrator='explicit'
#endif
  integer :: mol_trace_count=0,mol_trace_nodes(16)=1,mol_trace_cells(3,16)=1
#ifndef PELEF_MPI_EB_TREE
  integer :: cell_stage_unit,cell_sync_unit,cell_reaction_unit,cell_stencil_unit,cell_trace_event=0
  procedure(write_mol_cell_stage), pointer :: cell_stage_observer=>null()
  procedure(write_mol_cell_sync), pointer :: cell_sync_observer=>null()
  procedure(write_mol_cell_reaction), pointer :: cell_reaction_observer=>null()
#endif
  character(len=512) :: initial_state_file=''
  character(len=512) :: level_set_file='',section_file='',output_file='eb_tree.csv',checkpoint_file='',restart_file=''
  integer :: rebalance_interval_steps=0
  integer :: output_interval_steps=0,last_output_step=-1
  real(dp) :: rebalance_minimum_gain=0.1_dp
  logical :: simple_transport_enabled=.false.
  logical :: viscosity_enabled=.true.,thermal_conduction_enabled=.true.
  logical :: species_diffusion_enabled=.true.,barodiffusion_enabled=.true.
  logical :: transport_flags(4)
#ifdef PELEF_SELECTED_SIMPLE_TRANSPORT
  real(dp) :: simple_coefficients(4,2*selected_nspecies+selected_nspecies*(selected_nspecies+1)/2)
#endif
  logical :: pinned_eos_extrapolation=.false.
  logical :: pinned_eos_convergence=.false.
  logical :: separate_momentum=.false.,consistent_velocity=.false.,paired_pressure=.false.,connected_directions=.false.
  real(dp) :: pinned_eos_temperature_min=0._dp,pinned_eos_temperature_max=0._dp
  namelist /eb_tree/ root_shape,child_count,parent,ratio,lower,upper,domain_lower,domain_upper, &
    geometry_type,plane_axis,plane_position,extrusion_axis,level_set_file,section_file,initial_state_file, &
    mol_trace_count,mol_trace_nodes,mol_trace_cells, &
    time_step,final_time, &
    maximum_steps,maximum_step_retries,rtol,atol,chemistry_enabled,transport_enabled,simple_transport_enabled,riemann_solver, &
    viscosity_enabled,thermal_conduction_enabled,species_diffusion_enabled,barodiffusion_enabled, &
    inflow_mach,inflow_temperature,inflow_pressure,difmag,hydro_method,boundary_faces,inflow_density,inflow_velocity, &
    inflow_composition_count,inflow_composition_species,inflow_composition_mole_fractions, &
    regular_density,cut_density,initial_pressure,wall_temperature,embedded_wall_temperature,embedded_noslip, &
    time_integrator,sdc_iterations,sdc_reactor,mol_mode,mol_lax_friedrichs,separate_momentum,consistent_velocity, &
      paired_pressure,connected_directions, &
    rotating_frame_enabled, &
    dual_energy_enabled,allow_small_energy,allow_negative_energy,update_total_from_internal,dual_energy_eta2, &
    pinned_eos_extrapolation,pinned_eos_temperature_min,pinned_eos_temperature_max,pinned_eos_convergence, &
    rotation_axis,rotation_angular_speed,rotation_axis_point,external_forcing,redistribution, &
        state_redist_target_volume_fraction,state_redist_max_order, &
    automatic_timestep,initial_timestep_factor,maximum_timestep_growth,subcycling, &
    hydro_cfl,transport_cfl,chemistry_integrator,output_file,output_interval_steps,checkpoint_file,restart_file, &
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
  real(dp), allocatable :: force_total(:),step_force(:),coordinate_total(:),numerical_total(:)
  real(dp) :: frame_energy_step=0._dp
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:),current(:),boundary(:),wall(:),chemical(:),injected(:)
  real(dp) :: dt,theta,error,estimated_dt
  integer :: unit,ios,status,node,n(3),nv,retry,depth,event,common
  integer :: ctu_halo=6
  logical :: ok,stopped=.false.
  character(len=1024) :: input_path,message,budget_path
#ifdef PELEF_MPI_EB_TREE
  type(amr_tree_state_3d) :: mpi_state,mpi_temperature
  type(eb_geometry_3d), allocatable :: mpi_geometry(:)
  integer, allocatable :: mpi_owners(:)
  integer :: mpi_rank,mpi_size,mpi_ierr,momentum_min,momentum_max
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
  if(mol_trace_count<0.or.mol_trace_count>16) error stop 'Invalid MOL trace cell count'
  if(mol_trace_count>0) then
#ifdef PELEF_MPI_EB_TREE
    error stop 'MOL cell trace currently requires serial gas'
#else
    if(.not.mol_mode.or..not.dual_energy_enabled.or.time_integrator/='sdc') &
      error stop 'MOL cell trace requires held SDC MOL'
    if(.not.subcycling) error stop 'MOL cell trace currently requires subcycling'
    if(redistribution/='none') error stop 'MOL cell trace currently requires redistribution none'
    cell_stage_observer=>write_mol_cell_stage;cell_sync_observer=>write_mol_cell_sync
    cell_reaction_observer=>write_mol_cell_reaction
#endif
  endif
  if(len_trim(initial_state_file)>0) then
    if(len_trim(restart_file)>0) error stop 'Initial state file conflicts with restart'
    if(paths_alias(initial_state_file,output_file).or.paths_alias(initial_state_file,checkpoint_file)) &
      error stop 'Initial state file aliases output'
    if(rotating_frame_enabled.or..not.dual_energy_enabled) &
      error stop 'Initial state file requires nonrotating held energy'
#ifdef PELEF_MPI_EB_SPRAY
    error stop 'Initial state file is currently gas only'
#endif
  endif
#ifdef PELEF_MPI_EB_TREE
  call MPI_Allreduce(merge(1,0,separate_momentum)+merge(2,0,consistent_velocity)+merge(4,0, &
      paired_pressure)+merge(8,0,connected_directions), &
    momentum_min,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,mpi_ierr)
  if(mpi_ierr/=MPI_SUCCESS) error stop 'MPI momentum policy agreement failed'
  call MPI_Allreduce(merge(1,0,separate_momentum)+merge(2,0,consistent_velocity)+merge(4,0, &
      paired_pressure)+merge(8,0,connected_directions), &
    momentum_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpi_ierr)
  if(mpi_ierr/=MPI_SUCCESS.or.momentum_min/=momentum_max) error stop 'MPI EB-tree momentum policy mismatch'
#endif
  if(connected_directions.and..not.separate_momentum) error stop 'Connected directions requires separate momentum'
  if(consistent_velocity.and..not.separate_momentum) error stop 'Consistent velocity requires separate momentum'
  if(paired_pressure.and..not.consistent_velocity) error stop 'Paired pressure requires consistent velocity'
  if(paired_pressure.and.rotating_frame_enabled) error stop 'Paired pressure requires a nonrotating frame'
  if(separate_momentum) then
    if(time_integrator/='sdc'.or..not.mol_mode.or..not.dual_energy_enabled.or.redistribution/='none') &
      error stop 'Separate momentum requires held MOL without redistribution'
#ifdef PELEF_MPI_EB_SPRAY
    error stop 'Separate momentum is unavailable for EB spray'
#endif
  endif
#ifndef PELEF_MPI_EB_SPRAY
  call read_gas_les_input()
#if defined(PELEF_MPI_EB_LES) || defined(PELEF_SERIAL_EB_LES)
  if(separate_momentum.and.les_enabled) error stop 'Separate momentum is unavailable for EB LES'
#endif
#endif
  if(time_integrator/='split'.and.time_integrator/='sdc') error stop 'Invalid time integrator'
  call read_constant_transport(trim(input_path),constant_coefficients,hit_transport_parameters)
  if(allocated(hit_transport_parameters)) error stop 'HIT transport scales require initialized regular HIT input'
  if(allocated(constant_coefficients)) then
    if(time_integrator/='sdc'.or..not.transport_enabled) error stop 'EB constant transport requires enabled SDC transport'
  endif
  if(mol_lax_friedrichs.and..not.mol_mode) error stop 'MOL LF flux requires MOL mode'
  if(rotating_frame_enabled) then
    if(.not.mol_mode) error stop 'Rotating frame requires MOL mode'
    allocate(frame)
    frame%axis=rotation_axis+1
    frame%angular_speed=rotation_angular_speed
    frame%axis_point=rotation_axis_point
    if(.not.frame%is_valid()) error stop 'Invalid rotating frame'
  else
    if(rotation_axis/=2.or.rotation_angular_speed/=0._dp.or.any(rotation_axis_point/=0._dp)) &
      error stop 'Rotating frame controls require enable flag'
  end if
  if(mol_mode) then
    if(time_integrator/='sdc') error stop 'MOL requires the coupled time integrator'
  end if
  if(dual_energy_enabled) then
    if(.not.ieee_is_finite(dual_energy_eta2)) error stop 'Invalid dual-energy eta2'
    if(time_integrator/='sdc') error stop 'Held UEINT requires SDC'
    if(redistribution/='none'.and.redistribution/='state_redist') &
      error stop 'Held UEINT currently requires no redistribution or StateRedist'
    allocate(held_energy,reset_policy)
    reset_policy=dual_energy_policy(allow_small_energy,allow_negative_energy, &
      update_total_from_internal,dual_energy_eta2)
  else
    if(.not.allow_small_energy.or..not.allow_negative_energy.or..not.update_total_from_internal.or. &
      dual_energy_eta2/=1.e-4_dp) error stop 'Dual-energy reset controls require enable flag'
  endif
  if(sdc_iterations<1.or..not.all(ieee_is_finite(external_forcing))) error stop 'Invalid SDC controls'
  if(sdc_reactor/='native'.and.sdc_reactor/='cvode'.and.sdc_reactor/='cvode_full') &
    error stop 'Unknown EB tree SDC reactor'
  if(sdc_reactor=='cvode'.or.sdc_reactor=='cvode_full') then
#ifndef PELEF_SDC_CVODE
    error stop 'SDC CVODE backend is unavailable in this application'
#endif
    if(time_integrator/='sdc'.or..not.chemistry_enabled) error stop 'EB tree CVODE requires reacting SDC'
  end if
  if(redistribution/='none'.and.redistribution/='flux_redist'.and.redistribution/='state_redist') &
    error stop 'Invalid tree redistribution'
  if(state_redist_max_order/=0.and.state_redist_max_order/=2.and.state_redist_max_order/=4) &
    error stop 'StateRedist order must be 0, 2 or 4'
  if(.not.ieee_is_finite(state_redist_target_volume_fraction)) error stop 'Invalid redistribution target'
  if(state_redist_target_volume_fraction<=0.or.state_redist_target_volume_fraction>1) &
    error stop 'Invalid redistribution target'
  if(redistribution/='none') then
    if(time_integrator/='sdc') error stop 'Tree redistribution requires SDC'
  end if
  if(time_integrator=='sdc') then
    if(hydro_method/='plm_ctu'.or.riemann_solver/='pelec') error stop 'EB tree SDC requires PeleC CTU'
    if(any(boundary_faces=='legacy')) error stop 'EB tree SDC requires explicit boundary faces'
  else
    if(any(external_forcing/=0)) error stop 'External forcing requires SDC'
  end if
  if(.not.ieee_is_finite(embedded_wall_temperature).or.embedded_wall_temperature<0) &
    error stop 'Invalid embedded wall temperature'
  if(embedded_wall_temperature>0) then
    if(time_integrator/='sdc'.or..not.transport_enabled.or..not.thermal_conduction_enabled) &
      error stop 'Embedded wall heating requires SDC and thermal conduction'
    allocate(active_embedded_temperature,source=embedded_wall_temperature)
  end if
  if(embedded_noslip) then
    if(time_integrator/='sdc'.or..not.transport_enabled.or..not.viscosity_enabled) &
      error stop 'Embedded no-slip wall requires SDC and viscosity'
  end if
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
    if(.not.transport_enabled.or..not.thermal_conduction_enabled) &
      error stop 'Isothermal tree walls require thermal conduction'
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
  if (.not.all(ieee_is_finite([domain_lower,domain_upper,plane_position,time_step,final_time,rtol,atol,hydro_cfl, &
      transport_cfl]))) &
    error stop 'Nonfinite controls'
  if (any(domain_upper<=domain_lower).or.min(time_step,final_time,rtol,atol)<=0._dp) error stop 'Invalid controls'
  if (hydro_cfl<=0._dp.or.hydro_cfl>1._dp.or.transport_cfl<=0._dp.or.transport_cfl>.5_dp) &
    error stop 'Invalid CFL controls'
  if(.not.all(ieee_is_finite([initial_timestep_factor,maximum_timestep_growth]))) &
    error stop 'Nonfinite automatic timestep controls'
  if(initial_timestep_factor<=0.or.initial_timestep_factor>1.or.maximum_timestep_growth<1) &
    error stop 'Invalid automatic timestep controls'
  if (chemistry_integrator/='explicit'.and.chemistry_integrator/='implicit') error stop 'Invalid chemistry integrator'
  if (riemann_solver/='rusanov'.and.riemann_solver/='hllc'.and.riemann_solver/='pelec') error stop 'Invalid solver'
  if (len_trim(output_file)==0.or.checkpoint_interval_steps<0.or.output_interval_steps<0) &
    error stop 'Invalid output policy'
#ifdef PELEF_MPI_EB_SPRAY
  if(output_interval_steps>0) error stop 'Intermediate snapshots currently require the gas executable'
#endif
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
  config%state_redist_max_order=state_redist_max_order
  config%state_redist_target_volume_fraction=state_redist_target_volume_fraction
  config%redistribution=redistribution;config%riemann_solver=riemann_solver;config%difmag=difmag
  config%hydro_method=hydro_method
  config%sdc_reactor=sdc_reactor
  config%pinned_eos_extrapolation=pinned_eos_extrapolation
  config%pinned_eos_temperature_min=pinned_eos_temperature_min
  config%pinned_eos_temperature_max=pinned_eos_temperature_max
  config%time_integrator=time_integrator;config%sdc_iterations=sdc_iterations;config%external_forcing=external_forcing
  if(allocated(active_faces)) config%boundary_faces=active_faces
  config%regular_density=regular_density;config%cut_density=cut_density;config%initial_pressure=initial_pressure
  config%chemistry_enabled=chemistry_enabled;config%transport_enabled=transport_enabled
  config%chemistry_relative_tolerance=rtol;config%chemistry_absolute_tolerance=atol
  barodiffusion_enabled=barodiffusion_enabled.and.species_diffusion_enabled
  transport_flags=transport_enabled.and.[viscosity_enabled,thermal_conduction_enabled, &
    species_diffusion_enabled,barodiffusion_enabled]
  config%viscosity_enabled=viscosity_enabled;config%thermal_conduction_enabled=thermal_conduction_enabled
  config%species_diffusion_enabled=species_diffusion_enabled;config%barodiffusion_enabled=barodiffusion_enabled
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
#ifdef PELEF_MPI_EB_SPRAY
    if(time_integrator=='sdc'.and..not.mol_mode.and.redistribution=='state_redist') ctu_halo=10
#else
    if(time_integrator=='sdc'.and..not.mol_mode.and.redistribution=='state_redist') ctu_halo=10
#endif
    common=ctu_halo
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
#ifndef PELEF_SELECTED_SIMPLE_TRANSPORT
  if(simple_transport_enabled) error stop 'Simple transport requires a bundle with pinned fits'
#endif
#ifdef PELEF_SELECTED_EB_TREE
  call load_selected_thermo(species,ok)
  if (.not.ok) error stop 'Selected thermodynamic load failed'
  call load_selected_reactions(reactions,ok)
  if (.not.ok) error stop 'Selected reaction load failed'
  call load_selected_transport(transport_names,transport_geometries,transport_values,ok)
  if (.not.ok) error stop 'Selected transport load failed'
  if(selected_nreactions==0.and.chemistry_enabled) error stop 'Nonreacting EB tree requires chemistry disabled'
  call prepare_selected_mechanism('EB tree',selected_nspecies,selected_nreactions,species,reactions, &
    transport_names,transport_geometries,transport_values,transport,common_tmin,common_tmax,ok,message, &
    allow_nonreacting=selected_nreactions==0)
  if (.not.ok) then
    write(*,'(a)') trim(message)
    error stop 'Invalid selected mechanism'
  end if
#ifdef PELEF_SELECTED_SIMPLE_TRANSPORT
  if(simple_transport_enabled) then
    call load_selected_simple(simple_coefficients,ok)
    if(.not.ok) error stop 'Simple transport coefficient load failed'
    call set_simple_transport_fits(transport,simple_coefficients,ok)
    if(.not.ok) error stop 'Invalid Simple transport coefficients'
  end if
#endif
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
  if(.not.pinned_eos_extrapolation) then
    call validate_selected_temperature_range('EB tree',min(regular_t,cut_t),max(regular_t,cut_t), &
      common_tmin,common_tmax,ok,message)
    if (.not.ok) then
      write(*,'(a)') trim(message)
      error stop 'Selected initial temperature outside mechanism range'
    end if
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
  call configure_nasa7_eos_policy(species,pinned_eos_extrapolation,ok, &
    pinned_eos_temperature_min,pinned_eos_temperature_max,converged_inverse=pinned_eos_convergence)
  if(.not.ok) error stop 'Invalid EB-tree EOS policy'
#ifdef PELEF_MPI_EB_TREE
  call agree_mpi_reactive_eos_1d(species,MPI_COMM_WORLD,ok)
  if(.not.ok) error stop 'MPI EB-tree EOS policy mismatch'
#endif
  call prepare_inlet_compositions()
  nv=reactive_nvar(size(species));n=root_shape
  allocate(current(nv),boundary(nv),wall(nv),chemical(nv),injected(nv))
  allocate(force_total(nv),step_force(nv),coordinate_total(nv));force_total=0;step_force=0;coordinate_total=0
  allocate(numerical_total(nv));numerical_total=0
  injected=0
  root_config=config;root_config%nx=n(1);root_config%ny=n(2);root_config%nz=n(3)
#ifdef PELEF_MPI_EB_TREE
  call initialize_owned()
  if(len_trim(initial_state_file)>0) call import_initial_state()
  call integrate_mpi_eb_tree_3d(topology,mpi_owners,mpi_geometry,mpi_state,current,MPI_COMM_WORLD,ok)
#else
  allocate(u(nv,n(1),n(2),n(3)),t(n(1),n(2),n(3)))
  call initialize_reactive_eb_density_sheet_3d(species,root_config,coarse,u,t,ok,base_mole_fractions)
  if (.not.ok) error stop 'Initial state failed'
  if(allocated(active_faces)) then
    call prepare_prescribed_inlet(u(:,1,1,1),t(1,1,1),ok)
    if(.not.ok) error stop 'Invalid prescribed inlet EOS'
  end if
  call rotate_seed(u,coarse,ok)
  if(.not.ok) error stop 'Initial rotating energy failed'
  call initialize_reactive_eb_tree_3d(species,topology,geometry_levels(topology%level_count()),u,t,tree,ok,frame,held_energy)
  if (.not.ok) error stop 'Tree initialization failed'
  if(len_trim(initial_state_file)>0) call import_initial_state()
  call integrate_reactive_eb_tree_3d(tree,current,ok)
#endif
  if (.not.ok) error stop 'Initial inventory failed'
#ifdef PELEF_MPI_EB_SPRAY
  if(allocated(constant_coefficients)) then
    call set_constant_transport_coefficients(transport,constant_coefficients,ok)
    if(.not.ok) error stop 'Invalid particle constant transport controls'
  endif
  call initialize_particles()
  call add_particle_inventory(current,drops,component_mass)
#endif
  if(allocated(held_energy)) then
    record%internal_energy=held_energy;record%energy_policy=reset_policy
    record%energy_correction=numerical_energy
    record%internal_energy_pressure_work=internal_energy_pressure_work
  endif
  record%initial=current;record%scale=abs(current);record%boundary=0._dp*current
  record%wall=0._dp*current;record%chemical=0._dp*current
  if(time_integrator=='sdc') then
    call initialize_heat_snapshot()
    record%heat_release=sdc_heat
    record%forcing=force_total
    if(chemistry_enabled) then
      allocate(sdc_ir)
      allocate(sdc_ir%fields(size(tree%nodes)))
      do node=1,size(tree%nodes)
#ifdef PELEF_MPI_EB_TREE
        if(mpi_owners(node)/=mpi_rank) cycle
#endif
        n=tree%topology%nodes(node)%shape
        allocate(sdc_ir%fields(node)%values(nv-4,n(1),n(2),n(3)),source=0._dp)
      end do
      record%reaction_source=sdc_ir
    end if
  end if
#ifdef PELEF_MPI_EB_SPRAY
  if(len_trim(restart_file)>0) call spray_checkpoint(.false.)
#else
  if (len_trim(restart_file)>0) call restore()
#endif
#ifndef PELEF_MPI_EB_SPRAY
  initial_reactions_pending=mol_mode.and.chemistry_enabled.and.len_trim(restart_file)==0
#endif
#ifndef PELEF_MPI_EB_TREE
  if(mol_trace_count>0) call open_mol_cell_trace()
#endif
  if(output_interval_steps>0) then
    call balance()
    call output_snapshot()
  endif
  do
#ifdef PELEF_MPI_EB_SPRAY
    call inject_due_events()
#endif
    if(final_time_reached().or.record%steps>=maximum_steps) exit
    call refresh_ctu_geometry()
    dt=time_step
#ifdef PELEF_MPI_EB_SPRAY
    call compute_mpi_amr_spray_timestep(MPI_COMM_WORLD,topology,drops,domain_lower, &
      (domain_upper-domain_lower)/real(root_shape,dp),boundary_faces([1,3,5])=="periodic", &
      particle_options,maximum_injection_speed,subcycling,particle_dt,ok)
    if(.not.ok) error stop "EB particle timestep calculation failed"
    dt=min(dt,particle_dt)
#endif
    if (automatic_timestep) then
#ifdef PELEF_MPI_EB_TREE
      call mpi_eb_tree_timestep_3d(species,transport,topology,mpi_owners,mpi_geometry,mpi_state,mpi_temperature, &
        hydro_cfl,transport_cfl,transport_enabled, &
        viscosity_enabled,thermal_conduction_enabled,species_diffusion_enabled,estimated_dt,MPI_COMM_WORLD,ok,subcycling, &
        inflow_mach=active_mach,inflow_temperature=active_temperature,inflow_pressure=active_pressure, &
        boundary_faces=active_faces,inflow_state=inlet_state,wall_temperature=active_wall_temperature, &
        coupled_subcycling=time_integrator=="sdc",frame=frame, &
        constant_coefficients=constant_coefficients,internal_energy=held_energy)
#else
      call reactive_eb_tree_timestep_3d(species,transport,tree,hydro_cfl,transport_cfl,transport_enabled, &
        viscosity_enabled,thermal_conduction_enabled,species_diffusion_enabled,estimated_dt,ok,subcycling, &
        inflow_mach=active_mach,inflow_temperature=active_temperature,inflow_pressure=active_pressure, &
        boundary_faces=active_faces,inflow_state=inlet_state,wall_temperature=active_wall_temperature, &
        coupled_subcycling=time_integrator=="sdc",frame=frame, &
        constant_coefficients=constant_coefficients,internal_energy=held_energy)
#endif
      if (.not.ok) error stop 'Tree timestep estimate failed'
      dt=min(dt,estimated_dt)
    end if
    if(allocated(active_embedded_temperature).or.embedded_noslip) then
      if(.not.allocated(wall_exchange)) allocate(wall_exchange(reactive_nvar(size(species))))
#ifdef PELEF_MPI_EB_TREE
      call compute_mpi_eb_wall_transport_3d(species,transport,topology,mpi_owners,mpi_geometry,heat_geometry, &
        mpi_state,mpi_temperature,allocated(active_embedded_temperature),embedded_noslip,embedded_wall_temperature, &
        heat_rates,wall_exchange,heat_frequency,MPI_COMM_WORLD,ok, &
        boundary_faces=active_faces,inflow_state=inlet_state,inflow_mach=active_mach, &
        inflow_temperature=active_temperature,inflow_pressure=active_pressure,frame=frame, &
        constant_coefficients=constant_coefficients,internal_energy=held_energy)
#else
      call compute_eb_tree_wall_transport_3d(species,transport,tree,heat_geometry, &
        allocated(active_embedded_temperature),embedded_noslip,embedded_wall_temperature, &
        heat_rates,wall_exchange,heat_frequency,ok,boundary_faces=active_faces, &
        inflow_state=inlet_state,inflow_mach=active_mach,inflow_temperature=active_temperature, &
        inflow_pressure=active_pressure,frame=frame, &
        constant_coefficients=constant_coefficients,internal_energy=held_energy)
#endif
      if(.not.ok) error stop 'Embedded wall timestep estimate failed'
      if(heat_frequency>0) dt=min(dt,transport_cfl/heat_frequency)
    end if
#ifdef PELEF_MPI_EB_LES
    if(les_enabled) then
      call les_source(mpi_state,mpi_temperature,dt,-1,les_rates,les_flux,les_theta,MPI_COMM_WORLD,ok)
#ifdef PELEF_MPI_EB_SPRAY
      if(.not.ok) error stop 'Invalid particle LES state or non-regular geometry'
#else
      if(.not.ok) error stop 'Invalid LES state or non-regular geometry'
#endif
      les_spacing=[minval(mpi_geometry%dx),minval(mpi_geometry%dy),minval(mpi_geometry%dz)]
      if(les_diffusivity>0) dt=min(dt,les_cfl/(les_diffusivity*sum(1._dp/les_spacing**2)))
    end if
#endif
#ifdef PELEF_SERIAL_EB_LES
    if(les_enabled) then
      call estimate_serial_les_diffusivity()
      les_spacing=[minval(tree%nodes%geometry%dx),minval(tree%nodes%geometry%dy),minval(tree%nodes%geometry%dz)]
      if(les_diffusivity>0) dt=min(dt,les_cfl/(les_diffusivity*sum(1._dp/les_spacing**2)))
    endif
#endif
    if(automatic_timestep) then
      if(record%steps==0) then
        dt=initial_timestep_factor*dt
      else if(timestep_after_regrid) then
        dt=min(dt,previous_timestep)
      else
        ! Compare before multiplication to avoid overflowing an inactive cap.
        if(previous_timestep<dt/maximum_timestep_growth) dt=maximum_timestep_growth*previous_timestep
      endif
    endif
    ! Apply event/stop caps after initial shrink and accepted-step growth.
    if(automatic_timestep) then
      if(final_time-record%time-dt<.001_dp*dt) dt=final_time-record%time
    endif
    dt=min(dt,final_time-record%time)
#ifdef PELEF_MPI_EB_SPRAY
    if(next_injection<=injection_events) dt=min(dt,event_schedule%times(next_injection)-record%time)
#endif
    if (regrid_stage<regrid_count) dt=min(dt,regrid_times(regrid_stage+1)-record%time)
#ifdef PELEF_MPI_EB_TREE
    sample_work=.false.
    if (rebalance_interval_steps>0) sample_work=mod(record%steps+1,rebalance_interval_steps)==0
#endif
    do retry=0,maximum_step_retries
      if (record%time+dt<=record%time) error stop 'Time increment underflow'
#ifdef PELEF_SDC_CVODE
      if(sdc_reactor=='cvode') then
        call advance_step(advance_forced_reactor_cvode)
      else if(sdc_reactor=='cvode_full') then
        call advance_step(advance_forced_reactor_cvode_full)
      else
        call advance_step()
      end if
#else
      call advance_step()
#endif
      if (ok) exit
#ifdef PELEF_MPI_EB_TREE
      if(mpi_rank==0) then
#endif
        write(*,'(a,i0,a,es24.16,a,es24.16,a,i0)') 'Rejected tree step: ',record%steps+1, &
          ' time: ',record%time,' interval: ',dt,' retry: ',retry
#ifdef PELEF_MPI_EB_SPRAY
        write(*,'(a,i0)') 'Processed injection events at rejection: ',next_injection-1
#endif
#ifdef PELEF_MPI_EB_TREE
      end if
#endif
      if (retry==maximum_step_retries) error stop 'Tree step retry limit exhausted'
      dt=.5_dp*dt
    end do
    if(initial_reactions_pending) then
#ifdef PELEF_MPI_EB_TREE
      if(mpi_rank==0) &
#endif
        print *, 'Initialized MOL reaction history for first accepted step'
      initial_reactions_pending=.false.
    endif
#ifdef PELEF_MPI_EB_SPRAY
    split_created_total=split_created_total+step_split_created
#endif
    if(allocated(held_energy)) then
      ! This is UEINT volume work, distinct from the UEDEN conservation
      ! residual. A rejected native step never charges the history prefix.
      if(step_internal_energy_pressure_work>0) then
        if(internal_energy_pressure_work>huge(1._dp)-step_internal_energy_pressure_work) &
          error stop 'Internal-energy pressure-work history overflow'
      else if(step_internal_energy_pressure_work<0) then
        if(internal_energy_pressure_work< -huge(1._dp)-step_internal_energy_pressure_work) &
          error stop 'Internal-energy pressure-work history overflow'
      endif
      internal_energy_pressure_work=internal_energy_pressure_work+step_internal_energy_pressure_work
      if(.not.ieee_is_finite(internal_energy_pressure_work)) error stop 'Nonfinite internal-energy pressure work'
      record%internal_energy=held_energy;record%energy_correction=numerical_energy
      record%internal_energy_pressure_work=internal_energy_pressure_work
      numerical_total(5)=numerical_energy
    endif
    force_total=force_total+step_force
    coordinate_total(5)=coordinate_total(5)+frame_energy_step
    if(.not.ieee_is_finite(coordinate_total(5))) error stop "Nonfinite coordinate energy history"
    if(time_integrator=='sdc') then
      record%forcing=force_total
      if(allocated(sdc_ir)) record%reaction_source=sdc_ir
      if(allocated(sdc_heat)) record%heat_release=sdc_heat
    end if
    record%boundary=record%boundary+boundary;record%wall=record%wall+wall;record%chemical=record%chemical+chemical
    next_time=record%time+dt
    if (regrid_stage<regrid_count) then
      if (dt==regrid_times(regrid_stage+1)-record%time) next_time=regrid_times(regrid_stage+1)
    end if
    previous_timestep=dt;timestep_after_regrid=.false.
    record%time=next_time;record%steps=record%steps+1
    if (record%steps==1) record%minimum_dt=dt
    record%minimum_dt=min(record%minimum_dt,dt);record%minimum_theta=min(record%minimum_theta,theta)
#ifdef PELEF_MPI_EB_TREE
    if (sample_work) then
#ifdef PELEF_MPI_EB_SPRAY
      call rebalance_mpi_amr_spray(MPI_COMM_WORLD,topology,mpi_owners,nv,mpi_state,mpi_temperature, &
        drops,component_mass,domain_lower,(domain_upper-domain_lower)/real(root_shape,dp), &
        boundary_faces([1,3,5])=='periodic',measured_cost,rebalance_minimum_gain,particle_options%maximum_parcels, &
        redistributed,ok,wall_films=films,film_mass=film_mass,maximum_films=walls%maximum_films,nodal_phi=nodal_phi, &
        reaction_source=sdc_ir,heat_release=sdc_heat)
      if(ok.and.allocated(sdc_ir)) record%reaction_source=sdc_ir
      if(ok.and.allocated(sdc_heat)) record%heat_release=sdc_heat
#else
      if(allocated(sdc_ir)) then
        call rebalance_mpi_amr_tree_with_field_3d(topology,mpi_owners,nv,mpi_state,mpi_temperature,nv-4,sdc_ir, &
          measured_cost,rebalance_minimum_gain,MPI_COMM_WORLD,redistributed,ok,sdc_heat,held_energy)
        if(ok) then
          record%reaction_source=sdc_ir
          if(allocated(sdc_heat)) record%heat_release=sdc_heat
        end if
      else if(allocated(sdc_heat)) then
        call rebalance_mpi_amr_tree_with_field_3d(topology,mpi_owners,nv,mpi_state,mpi_temperature,1,sdc_heat, &
          measured_cost,rebalance_minimum_gain,MPI_COMM_WORLD,redistributed,ok,internal_energy=held_energy)
        if(ok) record%heat_release=sdc_heat
      else
      call rebalance_mpi_amr_tree_3d(topology,mpi_owners,nv,mpi_state,mpi_temperature, &
        measured_cost,rebalance_minimum_gain,MPI_COMM_WORLD,redistributed,ok)
      end if
#endif
      if(ok.and.allocated(held_energy)) record%internal_energy=held_energy
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
        timestep_after_regrid=.true.
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
    if(output_interval_steps>0) then
      if(mod(record%steps,output_interval_steps)==0) call output_snapshot()
    endif
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
  if (.not.final_time_reached().and..not.stopped) error stop 'Maximum steps reached before final time'
  call balance()
#ifdef PELEF_MPI_EB_SPRAY
  call output_particles()
#endif
  if(output_interval_steps>0.and.last_output_step/=record%steps) call output_snapshot()
  call output(trim(output_file))
#ifdef PELEF_MPI_EB_TREE
  if (mpi_rank==0) then
#endif
  write(*,'(a,i0,a,es24.16)') 'Completed steps: ',record%steps,' time: ',record%time
  if(time_integrator=='sdc') write(*,'(a)') 'SDC reactor: '//trim(sdc_reactor)
  if(mol_mode) write(*,'(a)') 'Hydrodynamic update: MOL'
  if(mol_lax_friedrichs) write(*,'(a)') 'MOL interface flux: Lax-Friedrichs'
  write(*,'(a,es24.16)') 'Composite inventory error: ',error
  write(*,'(a,es24.16)') 'Minimum accepted interval: ',record%minimum_dt
#ifdef PELEF_MPI_EB_SPRAY
  write(*,'(a,i0)') 'Particle ID watermark: ',last_particle_id
  write(*,'(a,i0)') 'Computational split parcels: ',split_created_total
#endif
#ifdef PELEF_MPI_EB_TREE
  end if
  call MPI_Finalize(mpi_ierr)
  if (mpi_ierr/=MPI_SUCCESS) error stop 'MPI finalize failed'
#endif
contains
  logical function final_time_reached()
    ! Fixed-step accumulation can finish a few ulps below the requested clock.
    ! Scale to physical time (not unity), and preserve the actual saved clock.
    final_time_reached=record%time>=final_time.or. &
      final_time-record%time<=8.0_dp*epsilon(1.0_dp)*max(abs(final_time),abs(record%time))
  end function final_time_reached
#ifndef PELEF_MPI_EB_TREE
  subroutine open_mol_cell_trace()
    integer :: q,p,c(3)
    character(len=1024) :: stage_path,sync_path,reaction_path,stencil_path
    stage_path=trim(output_file)//'.mol_stage.csv';sync_path=trim(output_file)//'.mol_sync.csv'
    reaction_path=trim(output_file)//'.mol_reaction.csv'
    stencil_path=trim(output_file)//'.mol_stencil.csv'
    if(paths_alias(stage_path,input_path).or.paths_alias(stage_path,initial_state_file).or. &
      paths_alias(stage_path,restart_file).or.paths_alias(stage_path,checkpoint_file).or. &
      paths_alias(sync_path,input_path).or.paths_alias(sync_path,initial_state_file).or. &
      paths_alias(sync_path,restart_file).or.paths_alias(sync_path,checkpoint_file).or. &
      paths_alias(reaction_path,input_path).or.paths_alias(reaction_path,initial_state_file).or. &
      paths_alias(reaction_path,restart_file).or.paths_alias(reaction_path,checkpoint_file).or. &
      paths_alias(stencil_path,input_path).or.paths_alias(stencil_path,initial_state_file).or. &
      paths_alias(stencil_path,restart_file).or.paths_alias(stencil_path,checkpoint_file)) &
      error stop 'MOL trace path aliases input or checkpoint'
    do q=1,mol_trace_count
      p=mol_trace_nodes(q);c=mol_trace_cells(:,q)
      if(p<1.or.p>size(tree%nodes)) error stop 'Invalid MOL trace node'
      if(any(c<1).or.any(c>tree%topology%nodes(p)%shape)) error stop 'Invalid MOL trace cell'
    enddo
    open(newunit=cell_stage_unit,file=trim(stage_path),status='replace',action='write',iostat=ios)
    if(ios/=0) error stop 'Cannot open MOL stage trace'
    write(cell_stage_unit,'(a)') 'step,event,stage,depth,time_start_s,interval_s,node,probe,i,j,k,component,axis,side,'// &
      'state,pressure_cell_Pa,pressure_left_Pa,pressure_right_Pa,raw_hydro_flux,centroid_hydro_flux,'// &
      'tree_flux_rate,hydro_rate,wall_rate,vfrac,aperture,'// &
      'rho_left,rho_right,normal_velocity_left,normal_velocity_right,sound_left,sound_right'
    open(newunit=cell_stencil_unit,file=trim(stencil_path),status='replace',action='write',iostat=ios)
    if(ios/=0) error stop 'Cannot open MOL stencil trace'
    write(cell_stencil_unit,'(a)',advance='no') &
      'step,event,stage,depth,time_start_s,interval_s,node,probe,i,j,k,axis,side,offset,sound,rho,u,v,w,p'
    do q=1,size(species)
      write(cell_stencil_unit,'(a,i0)',advance='no') ',Y_',q
    enddo
    write(cell_stencil_unit,*)
    open(newunit=cell_sync_unit,file=trim(sync_path),status='replace',action='write',iostat=ios)
    if(ios/=0) error stop 'Cannot open MOL sync trace'
    write(cell_sync_unit,'(a)') 'step,event,phase,depth,time_start_s,interval_s,node,probe,i,j,k,component,axis,side,'// &
      'state,integrated_flux,raw_increment,applied_increment,vfrac,aperture'
    open(newunit=cell_reaction_unit,file=trim(reaction_path),status='replace',action='write',iostat=ios)
    if(ios/=0) error stop 'Cannot open MOL reaction trace'
    write(cell_reaction_unit,'(a)') 'step,event,stage,depth,time_start_s,interval_s,node,probe,i,j,k,component,'// &
      'old_state,updated_state,forcing_rate,spatial_old,spatial_new,external_old,external_new,reaction_rate'
  end subroutine

  ! Attempt diagnostics, not accepted-step receipts. Component order is
  ! rho,momentum(3),rhoE,rhoY(:). Hydro rates exclude reaction/transport;
  ! tree_flux_rate includes the flux contributions supplied by the integrator.
  subroutine write_mol_cell_stage(stage,depth,time,interval,snapshot,diagnostics,flux)
    integer, intent(in) :: stage,depth
    real(dp), intent(in) :: time,interval
    type(eb_reactive_tree_3d), intent(in) :: snapshot
    type(eb_mol_diagnostics_3d), intent(in) :: diagnostics(:)
    type(eb_integrated_flux_3d), intent(in) :: flux(:)
    integer :: q,p,c(3),v(3),f(3),r(3),h,g,a,side,m,offset,donor(3)
    real(dp) :: pair(2),kinematics(6),raw,centroid,face,aperture
    cell_trace_event=cell_trace_event+1
    do q=1,mol_trace_count
      p=mol_trace_nodes(q);c=mol_trace_cells(:,q)
      if(p>size(snapshot%nodes)) error stop 'MOL trace node removed'
      if(any(c>snapshot%topology%nodes(p)%shape)) error stop 'MOL trace cell removed'
      if(.not.allocated(diagnostics(p)%rate)) cycle
      if(snapshot%nodes(p)%geometry%volume_fraction(c(1),c(2),c(3))<=0) cycle
      h=(size(diagnostics(p)%cell_pressure,1)-snapshot%topology%nodes(p)%shape(1))/2;g=h-3
      do a=1,3
      do side=0,1
        v=c;v(a)=v(a)-1+side;f=v+g;r=f+1
        ! Offset is relative to the right donor of this face (-2,-1,0,1).
        ! Input halo indexing is physical cell index plus h in each direction.
        do offset=-2,1
          donor=c+h;donor(a)=donor(a)+side+offset
          write(cell_stencil_unit,'(*(g0,:,","))') record%steps+1,cell_trace_event,stage,depth, &
            record%time+time,interval,p,q,c,a,side,offset, &
            diagnostics(p)%cell_sound(donor(1),donor(2),donor(3)), &
            diagnostics(p)%cell_primitive(:,donor(1),donor(2),donor(3))
        enddo
        do m=1,size(snapshot%nodes(p)%state,1)
          select case(a)
          case(1)
            aperture=snapshot%nodes(p)%geometry%x_face_fraction(v(1),v(2),v(3))
            pair=diagnostics(p)%pressure%x(:,r(1),r(2),r(3))
            kinematics=diagnostics(p)%kinematics%x(:,r(1),r(2),r(3))
            raw=diagnostics(p)%raw_flux%x(m,r(1),r(2),r(3))
            centroid=diagnostics(p)%centroid_flux%x(m,f(1),f(2),f(3));face=flux(p)%x(m,v(1),v(2),v(3))
          case(2)
            aperture=snapshot%nodes(p)%geometry%y_face_fraction(v(1),v(2),v(3))
            pair=diagnostics(p)%pressure%y(:,r(1),r(2),r(3))
            kinematics=diagnostics(p)%kinematics%y(:,r(1),r(2),r(3))
            raw=diagnostics(p)%raw_flux%y(m,r(1),r(2),r(3))
            centroid=diagnostics(p)%centroid_flux%y(m,f(1),f(2),f(3));face=flux(p)%y(m,v(1),v(2),v(3))
          case(3)
            aperture=snapshot%nodes(p)%geometry%z_face_fraction(v(1),v(2),v(3))
            pair=diagnostics(p)%pressure%z(:,r(1),r(2),r(3))
            kinematics=diagnostics(p)%kinematics%z(:,r(1),r(2),r(3))
            raw=diagnostics(p)%raw_flux%z(m,r(1),r(2),r(3))
            centroid=diagnostics(p)%centroid_flux%z(m,f(1),f(2),f(3));face=flux(p)%z(m,v(1),v(2),v(3))
          end select
          write(cell_stage_unit,'(*(g0,:,","))') record%steps+1,cell_trace_event,stage,depth, &
            record%time+time,interval,p,q,c,m,a,side,snapshot%nodes(p)%state(m,c(1),c(2),c(3)), &
            diagnostics(p)%cell_pressure(c(1)+h,c(2)+h,c(3)+h),pair,raw,centroid,face/interval, &
            diagnostics(p)%rate(m,c(1)+g,c(2)+g,c(3)+g),diagnostics(p)%wall_rate(m,c(1)+g,c(2)+g,c(3)+g), &
            snapshot%nodes(p)%geometry%volume_fraction(c(1),c(2),c(3)),aperture,kinematics
        enddo
      enddo
      enddo
    enddo
    flush(cell_stage_unit)
    flush(cell_stencil_unit)
  end subroutine

  subroutine write_mol_cell_reaction(stage,depth,time,interval,old,updated,forcing, &
      spatial_old,spatial_new,external_old,external_new,reaction)
    integer, intent(in) :: stage,depth
    real(dp), intent(in) :: time,interval
    type(eb_reactive_tree_3d), intent(in) :: old
    type(amr_tree_state_3d), intent(in) :: updated,forcing,spatial_old,spatial_new,external_old,external_new,reaction
    integer :: q,p,c(3),m,nvar
    real(dp) :: chemical_rate
    cell_trace_event=cell_trace_event+1
    do q=1,mol_trace_count
      p=mol_trace_nodes(q);c=mol_trace_cells(:,q)
      if(p>size(old%nodes)) error stop 'MOL trace node removed'
      if(depth>=0.and.old%topology%nodes(p)%level/=depth) cycle
      if(any(c>old%topology%nodes(p)%shape)) error stop 'MOL trace cell removed'
      if(old%nodes(p)%geometry%volume_fraction(c(1),c(2),c(3))<=0) cycle
      nvar=size(old%nodes(p)%state,1)
      do m=1,nvar
        chemical_rate=0
        if(m>=6) chemical_rate=reaction%fields(p)%values(m-5,c(1),c(2),c(3))
        if(m==5) chemical_rate=reaction%fields(p)%values(nvar-4,c(1),c(2),c(3))
        write(cell_reaction_unit,'(*(g0,:,","))') record%steps+1,cell_trace_event,stage,depth, &
          record%time+time,interval,p,q,c,m,old%nodes(p)%state(m,c(1),c(2),c(3)), &
          updated%fields(p)%values(m,c(1),c(2),c(3)),forcing%fields(p)%values(m,c(1),c(2),c(3)), &
          spatial_old%fields(p)%values(m,c(1),c(2),c(3)),spatial_new%fields(p)%values(m,c(1),c(2),c(3)), &
          external_old%fields(p)%values(m,c(1),c(2),c(3)),external_new%fields(p)%values(m,c(1),c(2),c(3)), &
          chemical_rate
      enddo
    enddo
    flush(cell_reaction_unit)
  end subroutine

  subroutine write_mol_cell_sync(phase,depth,time,interval,snapshot,stored,flux,scalar_flux, &
      raw_increment,applied_increment)
    integer, intent(in) :: phase,depth
    real(dp), intent(in) :: time,interval
    type(eb_reactive_tree_3d), intent(in) :: snapshot
    type(amr_tree_state_3d), intent(in) :: stored
    type(eb_integrated_flux_3d), intent(in) :: flux(:),scalar_flux(:)
    type(amr_tree_state_3d), optional, intent(in) :: raw_increment,applied_increment
    integer :: q,p,c(3),f(3),a,side,m
    real(dp) :: face,raw,applied,aperture
    cell_trace_event=cell_trace_event+1
    do q=1,mol_trace_count
      p=mol_trace_nodes(q);c=mol_trace_cells(:,q)
      if(p>size(snapshot%nodes)) error stop 'MOL trace node removed'
      if(any(c>snapshot%topology%nodes(p)%shape)) error stop 'MOL trace cell removed'
      if(snapshot%nodes(p)%geometry%volume_fraction(c(1),c(2),c(3))<=0) cycle
      do a=1,3
      do side=0,1
        f=c;f(a)=f(a)-1+side
        do m=1,size(snapshot%nodes(p)%state,1)
          raw=0;applied=0
          if(present(raw_increment)) raw=raw_increment%fields(p)%values(m,c(1),c(2),c(3))
          if(present(applied_increment)) applied=applied_increment%fields(p)%values(m,c(1),c(2),c(3))
          select case(a)
          case(1)
            face=flux(p)%x(m,f(1),f(2),f(3));aperture=snapshot%nodes(p)%geometry%x_face_fraction(f(1),f(2),f(3))
          case(2)
            face=flux(p)%y(m,f(1),f(2),f(3));aperture=snapshot%nodes(p)%geometry%y_face_fraction(f(1),f(2),f(3))
          case(3)
            face=flux(p)%z(m,f(1),f(2),f(3));aperture=snapshot%nodes(p)%geometry%z_face_fraction(f(1),f(2),f(3))
          end select
          write(cell_sync_unit,'(*(g0,:,","))') record%steps+1,cell_trace_event,phase,depth, &
            record%time+time,interval,p,q,c,m,a,side,snapshot%nodes(p)%state(m,c(1),c(2),c(3)), &
            face,raw,applied,snapshot%nodes(p)%geometry%volume_fraction(c(1),c(2),c(3)),aperture
        enddo
      enddo
      enddo
    enddo
    flush(cell_sync_unit)
  end subroutine
#endif
#if defined(PELEF_MPI_EB_TREE) && !defined(PELEF_MPI_EB_SPRAY)
  subroutine advance_gas_sdc(evaluator,evaluate_reaction,spatial_evaluator)
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    procedure(mpi_eb_sdc_diffusion_source), optional :: evaluator
    procedure(mpi_eb_sdc_spatial_diffusion_source), optional :: spatial_evaluator
    type(amr_tree_state_3d) :: saved_source
    if(initial_reactions_pending) then
      saved_source=sdc_ir
      call initialize_mpi_eb_tree_reaction_source_3d(species,reactions,topology,mpi_owners,mpi_geometry,mpi_state, &
        mpi_temperature,dt,subcycling,rtol,atol,sdc_ir,MPI_COMM_WORLD,ok, &
        evaluate_reaction=evaluate_reaction,internal_energy=held_energy,frame=frame)
      if(.not.ok) return
    endif
    if(allocated(held_energy)) then
      call advance_gas_sdc_fields(evaluator,evaluate_reaction,spatial_evaluator,held_energy,held_flux, &
        reset_policy,numerical_energy, &
        step_internal_energy_pressure_work)
    else
      call advance_gas_sdc_fields(evaluator,evaluate_reaction,spatial_evaluator)
    endif
    if(initial_reactions_pending.and..not.ok) call move_alloc(saved_source%fields,sdc_ir%fields)
  end subroutine advance_gas_sdc
  subroutine advance_gas_sdc_fields(evaluator,evaluate_reaction,spatial_evaluator,held,scalar_flux,policy,prefix, &
      pressure_work)
    type(amr_tree_state_3d), optional, intent(inout) :: held
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: scalar_flux(:)
    type(dual_energy_policy), optional, intent(in) :: policy
    real(dp), optional, intent(inout) :: prefix,pressure_work
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    procedure(mpi_eb_sdc_diffusion_source), optional :: evaluator
    procedure(mpi_eb_sdc_spatial_diffusion_source), optional :: spatial_evaluator
    integer :: callback_width
#ifdef PELEF_MPI_EB_LES
    real(dp), allocatable :: filter_weights(:)
    logical :: filter_valid
#endif
    callback_width=6
#ifdef PELEF_MPI_EB_LES
    if(present(spatial_evaluator).and.les_model=='dynamic') then
      call les_filter_weights(les_filter,les_ratio,filter_weights,filter_valid)
      if(.not.filter_valid) error stop 'Invalid Dynamic LES filter stencil'
      callback_width=(size(filter_weights)-1)/2+8
    end if
#endif
        if(subcycling) then
        call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,mpi_owners,mpi_geometry, &
          mpi_state,mpi_temperature,dt,sdc_iterations,chemistry_enabled, &
          transport_flags,external_forcing,rtol,atol, &
          extended_geometry,active_faces,boundary,wall,chemical,step_force,theta,MPI_COMM_WORLD,ok,sdc_ir, &
          difmag,inlet_state,active_mach,active_temperature,active_pressure,active_wall_temperature, &
          redistribution=redistribution,target_volume_fraction=state_redist_target_volume_fraction, &
          evaluate_reaction=evaluate_reaction,embedded_wall_temperature=active_embedded_temperature, &
          embedded_noslip=embedded_noslip,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
      state_redist_max_order=state_redist_max_order,evaluate_diffusion=evaluator, &
          evaluate_spatial_diffusion=spatial_evaluator,spatial_diffusion_halo_width=callback_width, &
          constant_coefficients=constant_coefficients,heat_release=sdc_heat, &
          internal_energy=held,internal_energy_flux=scalar_flux,energy_policy=policy,energy_correction=prefix, &
          internal_energy_pressure_work_change=pressure_work, &
        spatial_diffusion_is_source=les_enabled.and.mol_mode, &
        separate_momentum=separate_momentum,consistent_velocity=consistent_velocity, &
      paired_pressure=paired_pressure,connected_directions=connected_directions)
        else
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,mpi_owners,mpi_geometry, &
          mpi_state,mpi_temperature,extended_geometry,dt,sdc_iterations,chemistry_enabled, &
          transport_flags,external_forcing,rtol,atol, &
          active_faces,boundary,wall,chemical,step_force,theta,MPI_COMM_WORLD,ok,sdc_ir, &
          difmag,inlet_state,active_mach,active_temperature,active_pressure,active_wall_temperature, &
          redistribution=redistribution,target_volume_fraction=state_redist_target_volume_fraction, &
          evaluate_reaction=evaluate_reaction,embedded_wall_temperature=active_embedded_temperature, &
          embedded_noslip=embedded_noslip,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
      state_redist_max_order=state_redist_max_order,evaluate_diffusion=evaluator, &
          evaluate_spatial_diffusion=spatial_evaluator,spatial_diffusion_halo_width=callback_width, &
          constant_coefficients=constant_coefficients,heat_release=sdc_heat, &
          internal_energy=held,internal_energy_flux=scalar_flux,energy_policy=policy,energy_correction=prefix, &
          internal_energy_pressure_work_change=pressure_work, &
        spatial_diffusion_is_source=les_enabled.and.mol_mode, &
        separate_momentum=separate_momentum,consistent_velocity=consistent_velocity, &
      paired_pressure=paired_pressure,connected_directions=connected_directions)
        end if
  end subroutine advance_gas_sdc_fields
#endif
#ifndef PELEF_MPI_EB_SPRAY
  subroutine read_gas_les_input()
    integer :: stream,code,closed
    namelist /eb_tree_les/ les_enabled,les_model,les_cs,les_prt,les_sct,les_ci,les_cfl, &
      les_test_model,les_test_ratio,les_grid_ratio,les_harmonic_mean
    open(newunit=stream,file=trim(input_path),status='old',action='read',iostat=code)
    if(code/=0) error stop 'Cannot open gas LES input'
    read(stream,nml=eb_tree_les,iostat=code)
    close(stream,iostat=closed)
    if(code>0.or.closed/=0) error stop 'Invalid gas LES input'
    if(.not.les_enabled) return
#if !defined(PELEF_MPI_EB_LES) && !defined(PELEF_SERIAL_EB_LES)
    error stop 'Gas LES requires PELEF_ENABLE_SPRAY'
#endif
    if(time_integrator/='sdc') error stop 'Gas LES requires SDC'
    if(.not.all(ieee_is_finite([les_cs,les_prt,les_sct,les_ci,les_cfl]))) error stop 'Invalid gas LES controls'
    if(les_cs<0.or.les_prt<=0.or.les_sct<=0.or.les_ci<0.or.les_cfl<=0.or.les_cfl>.2_dp) &
      error stop 'Invalid gas LES controls'
    select case(les_model)
    case('smagorinsky','wale','vreman','dynamic')
    case default
      error stop 'Invalid gas LES model'
    end select
    if(les_model=='dynamic') then
      les_filter=trim(les_test_model);les_ratio=les_test_ratio;les_grid=les_grid_ratio;les_harmonic=les_harmonic_mean
    else
      if(les_test_model/='box'.or.les_test_ratio/=2.or.les_grid_ratio/=1.or.les_harmonic_mean) &
        error stop 'Gas LES filter controls require dynamic model'
    end if
  end subroutine read_gas_les_input
#endif
#ifndef PELEF_MPI_EB_TREE
  subroutine advance_serial_gas_sdc(evaluate_reaction,spatial_evaluator)
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    procedure(eb_sdc_spatial_diffusion_source), optional :: spatial_evaluator
    type(amr_tree_state_3d) :: saved_source
    if(initial_reactions_pending) then
      saved_source=sdc_ir
      call initialize_eb_tree_reaction_source_3d(species,reactions,tree,dt,subcycling,rtol,atol,sdc_ir,ok, &
        evaluate_reaction=evaluate_reaction,internal_energy=held_energy,frame=frame)
      if(.not.ok) return
    endif
    if(allocated(held_energy)) then
      call advance_serial_gas_sdc_fields(evaluate_reaction,spatial_evaluator,held_energy,held_flux,reset_policy, &
      numerical_energy, &
        step_internal_energy_pressure_work)
    else
      call advance_serial_gas_sdc_fields(evaluate_reaction,spatial_evaluator)
    endif
    if(initial_reactions_pending.and..not.ok) call move_alloc(saved_source%fields,sdc_ir%fields)
  end subroutine advance_serial_gas_sdc
  subroutine advance_serial_gas_sdc_fields(evaluate_reaction,spatial_evaluator,held,scalar_flux,policy,prefix,pressure_work)
    type(amr_tree_state_3d), optional, intent(inout) :: held
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: scalar_flux(:)
    type(dual_energy_policy), optional, intent(in) :: policy
    real(dp), optional, intent(inout) :: prefix,pressure_work
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    procedure(eb_sdc_spatial_diffusion_source), optional :: spatial_evaluator
    integer :: callback_width
#ifdef PELEF_SERIAL_EB_LES
    real(dp), allocatable :: weights(:)
    logical :: valid
#endif
    callback_width=6
#ifdef PELEF_SERIAL_EB_LES
    if(present(spatial_evaluator).and.les_model=='dynamic') then
      call les_filter_weights(les_filter,les_ratio,weights,valid)
      if(.not.valid) error stop 'Invalid Dynamic LES filter stencil'
      callback_width=(size(weights)-1)/2+8
    endif
#endif
        if(subcycling) then
        call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,dt,sdc_iterations,chemistry_enabled, &
          transport_flags,external_forcing,rtol,atol, &
          tree,extended_geometry,active_faces,boundary,wall,chemical,step_force,theta,ok,sdc_ir, &
          difmag,inlet_state,active_mach,active_temperature,active_pressure,active_wall_temperature, &
          redistribution=redistribution,target_volume_fraction=state_redist_target_volume_fraction, &
          evaluate_reaction=evaluate_reaction,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs, &
          embedded_wall_temperature=active_embedded_temperature,embedded_noslip=embedded_noslip,frame=frame, &
      state_redist_max_order=state_redist_max_order,evaluate_spatial_diffusion=spatial_evaluator, &
          spatial_diffusion_halo_width=callback_width, &
          constant_coefficients=constant_coefficients,heat_release=sdc_heat, &
          internal_energy=held,internal_energy_flux=scalar_flux,energy_policy=policy,energy_correction=prefix, &
          internal_energy_pressure_work_change=pressure_work, &
        spatial_diffusion_is_source=les_enabled.and.mol_mode, &
        separate_momentum=separate_momentum,consistent_velocity=consistent_velocity, &
      paired_pressure=paired_pressure,connected_directions=connected_directions, &
          observe_mol=cell_stage_observer,observe_mol_sync=cell_sync_observer, &
          observe_mol_reaction=cell_reaction_observer)
        else
        call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,sdc_iterations,chemistry_enabled, &
          transport_flags,external_forcing,rtol,atol, &
          tree,extended_geometry,active_faces,boundary,wall,chemical,step_force,theta,ok,sdc_ir, &
          difmag,inlet_state,active_mach,active_temperature,active_pressure,active_wall_temperature, &
          redistribution=redistribution,target_volume_fraction=state_redist_target_volume_fraction, &
          evaluate_reaction=evaluate_reaction,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs, &
          embedded_wall_temperature=active_embedded_temperature,embedded_noslip=embedded_noslip,frame=frame, &
      state_redist_max_order=state_redist_max_order,evaluate_spatial_diffusion=spatial_evaluator, &
          spatial_diffusion_halo_width=callback_width, &
          constant_coefficients=constant_coefficients,heat_release=sdc_heat, &
          internal_energy=held,internal_energy_flux=scalar_flux,energy_policy=policy,energy_correction=prefix, &
          internal_energy_pressure_work_change=pressure_work, &
        spatial_diffusion_is_source=les_enabled.and.mol_mode, &
        separate_momentum=separate_momentum,consistent_velocity=consistent_velocity, &
      paired_pressure=paired_pressure,connected_directions=connected_directions)
        end if
  end subroutine advance_serial_gas_sdc_fields
#endif
#ifdef PELEF_SERIAL_EB_LES
  subroutine serial_les_source(state,geometry,old_halos,interval,depth,rate,spatial_rate,flux,limiter,good,frame)
    type(eb_reactive_tree_3d), intent(in) :: state
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: old_halos
    real(dp), intent(in) :: interval
    integer, intent(in) :: depth
    type(amr_tree_state_3d), intent(inout) :: rate,spatial_rate
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    real(dp), intent(inout) :: limiter
    logical, intent(out) :: good
    type(rotating_frame_3d), optional, intent(in) :: frame
    call compute_regular_eb_les_spatial_sdc_3d(state,geometry,species,old_halos,interval,depth, &
      les_cs,les_prt,les_sct,rate,spatial_rate,flux,les_diffusivity,good,model=les_model,ci=les_ci, &
      test_model=les_filter,test_ratio=les_ratio,grid_ratio=les_grid,harmonic_mean=les_harmonic,frame=frame)
    if(good) limiter=1
  end subroutine serial_les_source
  subroutine estimate_serial_les_diffusivity()
    type(amr_tree_state_3d) :: old_halos,spatial_rates
    type(eb_geometry_3d), allocatable :: padded(:)
    real(dp), allocatable :: weights(:),gt(:,:,:),vf(:,:,:),metrics(:,:,:)
    integer :: p,m(3),width,padding
    logical :: valid
    width=6
    if(les_model=='dynamic') then
      call les_filter_weights(les_filter,les_ratio,weights,valid)
      if(.not.valid) error stop 'Invalid Dynamic LES filter stencil'
      width=(size(weights)-1)/2+8
    endif
    allocate(old_halos%fields(size(tree%nodes)),padded(size(tree%nodes)))
    do p=1,size(tree%nodes)
      m=tree%topology%nodes(p)%shape
      padding=(extended_geometry(p)%nx-m(1))/2
      call extract_eb_geometry_3d(extended_geometry(p),[padding-5,padding-5,padding-5], &
        m+padding+6,padded(p),valid)
      if(.not.valid) error stop 'Invalid LES prepared geometry'
      if(any(padded(p)%volume_fraction/=1)) error stop 'Invalid LES state or non-regular geometry'
      allocate(old_halos%fields(p)%values(nv,1-width:m(1)+width,1-width:m(2)+width,1-width:m(3)+width), &
        gt(1-width:m(1)+width,1-width:m(2)+width,1-width:m(3)+width), &
        vf(1-width:m(1)+width,1-width:m(2)+width,1-width:m(3)+width), &
        metrics(m(1)+2*width,m(2)+2*width,m(3)+2*width),source=1._dp)
      call fill_reactive_eb_tree_halo_3d(species,tree,p,old_halos%fields(p)%values,gt,valid, &
        boundary_faces=active_faces,inflow_state=inlet_state,inflow_mach=active_mach, &
        inflow_temperature=active_temperature,inflow_pressure=active_pressure,ghost_volume_fraction=vf, &
        halo_width=width,target_volume_fraction=metrics,restrict_snapshot=.false.,frame=frame)
      if(.not.valid) error stop 'Invalid LES timestep donors'
      deallocate(gt,vf,metrics)
    enddo
    call compute_regular_eb_les_spatial_sdc_3d(tree,padded,species,old_halos,dt,-1,les_cs,les_prt,les_sct, &
      les_rates,spatial_rates,les_flux,les_diffusivity,valid,model=les_model,ci=les_ci, &
      test_model=les_filter,test_ratio=les_ratio,grid_ratio=les_grid,harmonic_mean=les_harmonic,frame=frame)
    if(.not.valid) error stop 'Invalid LES state or non-regular geometry'
  end subroutine estimate_serial_les_diffusivity
#endif
  subroutine advance_step(evaluate_reaction)
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    frame_energy_step=0;step_internal_energy_pressure_work=0
#ifdef PELEF_MPI_EB_TREE
      if (sample_work) then
        call begin_tree_work(size(topology%nodes),ok)
        call MPI_Allreduce(ok,work_ready,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,mpi_ierr)
        if (mpi_ierr/=MPI_SUCCESS.or..not.work_ready) error stop 'Cannot allocate tree work counters'
      end if
#ifdef PELEF_MPI_EB_SPRAY
      if(time_integrator=='sdc') then
        if(les_enabled) then
          if(redistribution=='state_redist') then
            call advance_particle_sdc(evaluate_reaction=evaluate_reaction,spatial_evaluator=les_spatial_source)
          else
            call advance_particle_sdc(les_source,evaluate_reaction)
          endif
        else
          call advance_particle_sdc(evaluate_reaction=evaluate_reaction)
        end if
      else if(wall_impacts) then
      call advance_mpi_eb_spray_reactive(MPI_COMM_WORLD,topology,mpi_owners,mpi_geometry,nodal_phi, &
        species,reactions,transport,liquids,particle_options,mpi_state,mpi_temperature,drops,component_mass, &
        domain_lower,(domain_upper-domain_lower)/real(root_shape,dp),dt,riemann_solver, &
        chemistry_enabled,transport_enabled,viscosity_enabled,thermal_conduction_enabled, &
        species_diffusion_enabled,barodiffusion_enabled,rtol,atol,particle_steps, &
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
        chemistry_enabled,transport_enabled,viscosity_enabled,thermal_conduction_enabled, &
        species_diffusion_enabled,barodiffusion_enabled,rtol,atol,particle_steps, &
        theta,boundary,wall,chemical,ok,chemistry_integrator,subcycling,last_id=last_particle_id, &
        split_created=step_split_created,particle_outflow=particle_outflow, &
        acceleration=acceleration,source_change=step_force, &
        inflow_mach=active_mach,inflow_temperature=active_temperature,inflow_pressure=active_pressure, &
        difmag=difmag,extended_geometry=extended_geometry, &
        boundary_faces=active_faces,inflow_state=inlet_state,wall_temperature=active_wall_temperature)
      end if
#else
      if(time_integrator=='sdc') then
#ifdef PELEF_MPI_EB_LES
        if(les_enabled) then
          if(redistribution=='state_redist'.or.mol_mode) then
            call advance_gas_sdc(evaluate_reaction=evaluate_reaction,spatial_evaluator=les_spatial_source)
          else
            call advance_gas_sdc(evaluator=les_source,evaluate_reaction=evaluate_reaction)
          endif
        else
#endif
          call advance_gas_sdc(evaluate_reaction=evaluate_reaction)
#ifdef PELEF_MPI_EB_LES
        endif
#endif
      else
        call advance_mpi_eb_tree_full_3d(species,reactions,transport,topology,mpi_owners,mpi_geometry, &
          riemann_solver,dt,chemistry_enabled,rtol,atol,transport_enabled,viscosity_enabled,thermal_conduction_enabled, &
          species_diffusion_enabled,barodiffusion_enabled,mpi_state,mpi_temperature,boundary,wall,chemical, &
          theta,MPI_COMM_WORLD,ok,chemistry_integrator,subcycling,inflow_mach=active_mach, &
          inflow_temperature=active_temperature,inflow_pressure=active_pressure,difmag=difmag, &
          extended_geometry=extended_geometry,boundary_faces=active_faces,inflow_state=inlet_state, &
          wall_temperature=active_wall_temperature)
      endif
#endif
      if(sample_work) call end_tree_work(measured_cost)
#else
      if(time_integrator=='sdc') then
#ifdef PELEF_SERIAL_EB_LES
        if(les_enabled) then
          call advance_serial_gas_sdc(evaluate_reaction,serial_les_source)
        else
#endif
          call advance_serial_gas_sdc(evaluate_reaction)
#ifdef PELEF_SERIAL_EB_LES
        endif
#endif
      else
      call advance_reactive_eb_tree_full_3d(species,reactions,transport,riemann_solver,dt,chemistry_enabled, &
        rtol,atol,transport_enabled,viscosity_enabled,thermal_conduction_enabled, &
        species_diffusion_enabled,barodiffusion_enabled,tree,boundary,wall,chemical,theta,ok, &
        chemistry_integrator,subcycling, &
        inflow_mach=active_mach,inflow_temperature=active_temperature,inflow_pressure=active_pressure, &
        difmag=difmag,extended_geometry=extended_geometry, &
        boundary_faces=active_faces,inflow_state=inlet_state,wall_temperature=active_wall_temperature)
      end if
#endif
  end subroutine advance_step
#ifdef PELEF_MPI_EB_SPRAY
  subroutine advance_particle_sdc(evaluator,evaluate_reaction,spatial_evaluator)
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    procedure(mpi_eb_sdc_diffusion_source), optional :: evaluator
    procedure(mpi_eb_sdc_spatial_diffusion_source), optional :: spatial_evaluator
    real(dp), allocatable :: filter_weights(:),held_prefix,held_pressure_work
    type(dual_energy_policy), allocatable :: held_policy
    integer :: callback_width
    logical :: filter_valid
    if(allocated(held_energy)) then
      allocate(held_prefix,held_pressure_work,held_policy)
      held_policy=reset_policy
      held_prefix=numerical_energy;held_pressure_work=0
    endif
    callback_width=6
    if(present(spatial_evaluator).and.les_model=='dynamic') then
      call les_filter_weights(les_filter,les_ratio,filter_weights,filter_valid)
      if(.not.filter_valid) error stop 'Invalid Dynamic LES filter stencil'
      callback_width=(size(filter_weights)-1)/2+8
    endif
      if(subcycling) then
      if(wall_impacts) then
      call advance_mpi_eb_spray_sdc_subcycle(MPI_COMM_WORLD,topology,mpi_owners,mpi_geometry,extended_geometry,nodal_phi, &
        species,reactions,transport,liquids,particle_options,mpi_state,mpi_temperature,drops,component_mass, &
        domain_lower,(domain_upper-domain_lower)/real(root_shape,dp),dt,sdc_iterations,chemistry_enabled, &
        transport_flags,external_forcing, &
        rtol,atol,active_faces,boundary,wall,chemical,step_force,theta,particle_steps,last_particle_id,ok, &
        reaction_source=sdc_ir,difmag=difmag,redistribution=redistribution, &
        target_volume_fraction=state_redist_target_volume_fraction, &
        inflow_state=inlet_state,inflow_mach=active_mach,inflow_temperature=active_temperature, &
        inflow_pressure=active_pressure,wall_temperature=active_wall_temperature, &
        particle_acceleration=acceleration,gas_acceleration=acceleration,split_created=step_split_created, &
        films=films,film_mass=film_mass,film_settings=walls,film_impacts=.true.,model_rejected=model_rejected, &
        evaluate_diffusion=evaluator,evaluate_reaction=evaluate_reaction, &
        embedded_wall_temperature=active_embedded_temperature,embedded_noslip=embedded_noslip, &
        mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
      state_redist_max_order=state_redist_max_order,evaluate_spatial_diffusion=spatial_evaluator, &
      spatial_diffusion_halo_width=callback_width, &
        constant_coefficients=constant_coefficients,heat_release=sdc_heat,frame_energy_change=frame_energy_step, &
        internal_energy=held_energy,energy_policy=held_policy,energy_correction=held_prefix, &
        internal_energy_pressure_work_change=held_pressure_work)
      if(model_rejected) error stop 'EB wall impact outside admissible empirical model range'
      else
      call advance_mpi_eb_spray_sdc_subcycle(MPI_COMM_WORLD,topology,mpi_owners,mpi_geometry,extended_geometry,nodal_phi, &
        species,reactions,transport,liquids,particle_options,mpi_state,mpi_temperature,drops,component_mass, &
        domain_lower,(domain_upper-domain_lower)/real(root_shape,dp),dt,sdc_iterations,chemistry_enabled, &
        transport_flags,external_forcing, &
        rtol,atol,active_faces,boundary,wall,chemical,step_force,theta,particle_steps,last_particle_id,ok, &
        reaction_source=sdc_ir,difmag=difmag,redistribution=redistribution, &
        target_volume_fraction=state_redist_target_volume_fraction, &
        inflow_state=inlet_state,inflow_mach=active_mach,inflow_temperature=active_temperature, &
        inflow_pressure=active_pressure,wall_temperature=active_wall_temperature, &
        particle_acceleration=acceleration,gas_acceleration=acceleration,split_created=step_split_created, &
        evaluate_diffusion=evaluator,evaluate_reaction=evaluate_reaction, &
        embedded_wall_temperature=active_embedded_temperature,embedded_noslip=embedded_noslip, &
        mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
      state_redist_max_order=state_redist_max_order,evaluate_spatial_diffusion=spatial_evaluator, &
      spatial_diffusion_halo_width=callback_width, &
        constant_coefficients=constant_coefficients,heat_release=sdc_heat,frame_energy_change=frame_energy_step, &
        internal_energy=held_energy,energy_policy=held_policy,energy_correction=held_prefix, &
        internal_energy_pressure_work_change=held_pressure_work)
      end if
      else if(wall_impacts) then
      call advance_mpi_eb_spray_sdc(MPI_COMM_WORLD,topology,mpi_owners,mpi_geometry,extended_geometry,nodal_phi, &
        species,reactions,transport,liquids,particle_options,mpi_state,mpi_temperature,drops,component_mass, &
        domain_lower,(domain_upper-domain_lower)/real(root_shape,dp),dt,sdc_iterations,chemistry_enabled, &
        transport_flags,external_forcing, &
        rtol,atol,active_faces,boundary,wall,chemical,step_force,theta,particle_steps,last_particle_id,ok, &
        reaction_source=sdc_ir,difmag=difmag,redistribution=redistribution, &
        target_volume_fraction=state_redist_target_volume_fraction, &
        inflow_state=inlet_state,inflow_mach=active_mach,inflow_temperature=active_temperature, &
        inflow_pressure=active_pressure,wall_temperature=active_wall_temperature, &
        split_created=step_split_created,particle_acceleration=acceleration,gas_acceleration=acceleration, &
        films=films,film_mass=film_mass,film_settings=walls,film_impacts=.true.,model_rejected=model_rejected, &
        evaluate_diffusion=evaluator,evaluate_reaction=evaluate_reaction, &
        embedded_wall_temperature=active_embedded_temperature,embedded_noslip=embedded_noslip, &
        mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
      state_redist_max_order=state_redist_max_order,evaluate_spatial_diffusion=spatial_evaluator, &
      spatial_diffusion_halo_width=callback_width, &
        constant_coefficients=constant_coefficients,heat_release=sdc_heat,frame_energy_change=frame_energy_step, &
        internal_energy=held_energy,energy_policy=held_policy,energy_correction=held_prefix, &
        internal_energy_pressure_work_change=held_pressure_work)
      if(model_rejected) error stop 'EB wall impact outside admissible empirical model range'
      else
      call advance_mpi_eb_spray_sdc(MPI_COMM_WORLD,topology,mpi_owners,mpi_geometry,extended_geometry,nodal_phi, &
        species,reactions,transport,liquids,particle_options,mpi_state,mpi_temperature,drops,component_mass, &
        domain_lower,(domain_upper-domain_lower)/real(root_shape,dp),dt,sdc_iterations,chemistry_enabled, &
        transport_flags,external_forcing, &
        rtol,atol,active_faces,boundary,wall,chemical,step_force,theta,particle_steps,last_particle_id,ok, &
        reaction_source=sdc_ir,difmag=difmag,redistribution=redistribution, &
        target_volume_fraction=state_redist_target_volume_fraction, &
        inflow_state=inlet_state,inflow_mach=active_mach,inflow_temperature=active_temperature, &
        inflow_pressure=active_pressure,wall_temperature=active_wall_temperature, &
        split_created=step_split_created,particle_acceleration=acceleration,gas_acceleration=acceleration, &
        evaluate_diffusion=evaluator,evaluate_reaction=evaluate_reaction, &
        embedded_wall_temperature=active_embedded_temperature,embedded_noslip=embedded_noslip, &
        mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
      state_redist_max_order=state_redist_max_order,evaluate_spatial_diffusion=spatial_evaluator, &
      spatial_diffusion_halo_width=callback_width, &
        constant_coefficients=constant_coefficients,heat_release=sdc_heat,frame_energy_change=frame_energy_step, &
        internal_energy=held_energy,energy_policy=held_policy,energy_correction=held_prefix, &
        internal_energy_pressure_work_change=held_pressure_work)
      end if
    if(ok.and.allocated(held_energy)) then
      numerical_energy=held_prefix
      step_internal_energy_pressure_work=held_pressure_work
    endif
  end subroutine advance_particle_sdc
#endif
#ifdef PELEF_MPI_EB_LES
  subroutine les_spatial_source(state,temperature,geometry,old_halos,interval,depth, &
      rate,spatial_rate,flux,limiter,comm,good,frame)
    type(amr_tree_state_3d), intent(in) :: state,temperature,old_halos
    type(eb_geometry_3d), intent(in) :: geometry(:)
    real(dp), intent(in) :: interval
    integer, intent(in) :: depth
    type(amr_tree_state_3d), intent(inout) :: rate,spatial_rate
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    real(dp), intent(inout) :: limiter
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    type(rotating_frame_3d), optional, intent(in) :: frame
    call compute_mpi_regular_eb_les_spatial_sdc_3d(topology,mpi_owners,geometry,species,old_halos, &
      interval,depth,les_cs,les_prt,les_sct,rate,spatial_rate,flux,les_diffusivity,comm,good, &
      model=les_model,ci=les_ci,test_model=les_filter,test_ratio=les_ratio,grid_ratio=les_grid, &
      harmonic_mean=les_harmonic,frame=frame)
    if(good) limiter=1
  end subroutine les_spatial_source
  subroutine les_source(state,temperature,interval,depth,rate,flux,limiter,comm,good, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,fractions)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: interval
    integer, intent(in) :: depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    real(dp), intent(inout) :: limiter
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    real(dp), optional, intent(in) :: fractions(:)
    if(allocated(frame)) then
    call compute_mpi_regular_eb_les_sdc_3d(topology,mpi_owners,mpi_geometry,species,state,temperature, &
      (domain_upper-domain_lower)/real(root_shape,dp),active_faces,interval,depth,les_cs,les_prt,les_sct, &
      rate,flux,les_diffusivity,comm,good,ancestor_start,ancestor_end,temperature_start,temperature_end,fractions, &
      model=les_model,ci=les_ci,test_model=les_filter,test_ratio=les_ratio,grid_ratio=les_grid,harmonic_mean=les_harmonic, &
      inflow_state=inlet_state,inflow_mach=active_mach,inflow_temperature=active_temperature, &
      inflow_pressure=active_pressure,wall_temperature=active_wall_temperature, &
      frame=frame,domain_lower=domain_lower)
    else
    call compute_mpi_regular_eb_les_sdc_3d(topology,mpi_owners,mpi_geometry,species,state,temperature, &
      (domain_upper-domain_lower)/real(root_shape,dp),active_faces,interval,depth,les_cs,les_prt,les_sct, &
      rate,flux,les_diffusivity,comm,good,ancestor_start,ancestor_end,temperature_start,temperature_end,fractions, &
      model=les_model,ci=les_ci,test_model=les_filter,test_ratio=les_ratio,grid_ratio=les_grid,harmonic_mean=les_harmonic, &
      inflow_state=inlet_state,inflow_mach=active_mach,inflow_temperature=active_temperature, &
      inflow_pressure=active_pressure,wall_temperature=active_wall_temperature)
    end if
    if(good) limiter=1
  end subroutine les_source
#endif
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
    character(len=16) :: mixture_conductivity_models(16)='supplied'
    character(len=24) :: mixture_viscosity_models(16)='supplied'
    real(dp) :: mixture_conductivities(16)=0,mixture_conductivity_coefficients(4,16)=0
    real(dp) :: mixture_viscosities(16)=0,mixture_viscosity_coefficients(4,16)=0
    real(dp) :: mixture_conductivity_temperature_ranges(2,16)=0
    real(dp) :: mixture_viscosity_temperature_ranges(2,16)=0
    character(len=16) :: mixture_vapor_pressure_models(16)='clausius',mixture_boiling_models(16)='reject'
    real(dp) :: mixture_density_coefficients(4,16)=0,mixture_density_temperature_ranges(2,16)=0
    real(dp) :: mixture_cp_coefficients(4,16)=0,mixture_cp_temperature_ranges(2,16)=0
    real(dp) :: mixture_antoine_coefficients(4,16)=0,mixture_antoine_temperature_ranges(2,16)=0
    real(dp) :: mixture_molecular_weights(16)=0
    real(dp) :: mixture_critical_temperatures(16)=0,mixture_normal_boiling_temperatures(16)=0
    type(spray_rate_profile) :: rate_profiles(16)
    character(len=:), allocatable :: profile_context
    character(len=4096) :: context_line
    character(len=1024) :: profile_files(16)
    character(len=1024) :: event_file
    namelist /eb_spray_events_file/ event_file
    namelist /eb_spray_rate_files/ profile_files
    namelist /eb_spray_les/ les_enabled,les_model,les_cs,les_prt,les_sct,les_ci,les_cfl, &
      les_test_model,les_test_ratio,les_grid_ratio,les_harmonic_mean
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
      mixture_cp_models, &
      mixture_viscosity_models,mixture_viscosities,mixture_viscosity_coefficients, &
      mixture_viscosity_temperature_ranges, &
      mixture_conductivity_models,mixture_conductivities,mixture_conductivity_coefficients, &
      mixture_conductivity_temperature_ranges,mixture_cp_coefficients,mixture_cp_temperature_ranges, &
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
    particle_options%particle_cfl=.5_dp
    open(newunit=stream,file=trim(spray_input),status='old',action='read',iostat=code)
    if(code/=0) error stop 'Cannot read MPI EB spray input'
    read(stream,nml=eb_spray,iostat=code)
    if(code/=0) error stop 'Invalid MPI EB spray input'
    if(allocated(active_faces).and..not.particle_outflow) &
      error stop 'Explicit MPI EB spray boundaries require particle_outflow'
    event_file=''
    rewind(stream)
    read(stream,nml=eb_spray_events_file,iostat=code)
    if(code>0) error stop 'Invalid EB spray event file input'
    profile_files=''
    rewind(stream)
    read(stream,nml=eb_spray_rate_files,iostat=code)
    if(code>0) error stop 'Invalid EB spray rate files input'
    rewind(stream)
    read(stream,nml=eb_spray_les,iostat=code)
    if(code>0) error stop 'Invalid EB spray LES input'
    if(les_enabled) then
      if(time_integrator/='sdc') error stop 'Particle LES requires SDC'
      if(.not.ieee_is_finite(les_cfl)) error stop 'Invalid particle LES CFL'
      if(les_cfl<=0.or.les_cfl>.2_dp) error stop 'Invalid particle LES CFL'
      if(les_model=='dynamic') then
        les_filter=trim(les_test_model);les_ratio=les_test_ratio;les_grid=les_grid_ratio;les_harmonic=les_harmonic_mean
      else
        if(les_test_model/='box'.or.les_test_ratio/=2.or.les_grid_ratio/=1.or.les_harmonic_mean) &
          error stop 'Particle LES filter controls require dynamic model'
      end if
    end if
    close(stream,iostat=closed)
    if(closed/=0) error stop 'Invalid MPI EB spray input'
    if(port_count<1.or.port_count>16) error stop 'Invalid port count'
    if(any(profile_files(port_count+1:)/='')) error stop 'Unused EB rate files'
    do port=1,port_count
      if(len_trim(profile_files(port))==0) cycle
      if(port_profile_points(port)/=0.or.any(port_profile_times(:,port)/=0).or. &
        any(port_profile_rates(:,port)/=0)) error stop 'Rate file conflicts with inline profile'
      call read_spray_rate_profile(trim(profile_files(port)),rate_profiles(port),ok)
      if(.not.ok) error stop 'Invalid EB spray rate file'
      port_profile_points(port)=size(rate_profiles(port)%times)
      if(port_profile_points(port)<=64) then
        port_profile_times(:port_profile_points(port),port)=rate_profiles(port)%times
        port_profile_rates(:port_profile_points(port),port)=rate_profiles(port)%rates
      end if
    end do
    do port=1,port_count
      if(allocated(rate_profiles(port)%times)) cycle
      q=port_profile_points(port)
      if(q/=0.and.(q<2.or.q>64)) error stop 'Invalid inline port profile size'
      rate_profiles(port)%times=port_profile_times(:q,port)
      rate_profiles(port)%rates=port_profile_rates(:q,port)
    end do
    if(nc<1.or.nc>16.or.number<1.or.number>125000) error stop 'Invalid MPI EB spray dimensions'
    if(.not.all(ieee_is_finite(fractions)).or.any(fractions<0)) error stop 'Invalid liquid fractions'
    if(abs(sum(fractions(:nc))-1)>1e-12_dp.or.any(fractions(nc+1:)/=0)) error stop 'Invalid fraction sum'
    if(number>particle_options%maximum_parcels) error stop 'Particle capacity exceeded'
    if(.not.all(ieee_is_finite(acceleration))) error stop 'Invalid MPI EB acceleration'
    if(any(abs(acceleration)>1e6_dp)) error stop 'Invalid MPI EB acceleration'
    if(wall_impacts) then
      if(.not.valid_wall_options(walls,[character(len=16)::'slip_wall','slip_wall','slip_wall', &
        'slip_wall','slip_wall','slip_wall'])) error stop 'Invalid EB wall options'
      if(any(walls%temperature/=walls%temperature(1))) error stop 'EB wall temperature must be uniform'
    else if(walls%evolve) then
      error stop 'EB film evolution requires wall impacts'
    end if
    if(len_trim(event_file)>0) then
      if(injection_events/=1.or.any(injection_times/=0).or.any(injection_scales/=1).or. &
          any(event_ports/=1).or.any(injection_end_times/=-1)) error stop 'Event file conflicts with inline schedule'
      call read_spray_event_schedule(trim(event_file),event_schedule,ok)
      if(.not.ok) error stop 'Invalid EB spray event file'
      injection_events=size(event_schedule%times)
    else
      if(injection_events<1.or.injection_events>1024) error stop 'Invalid injection event count'
      if(.not.all(ieee_is_finite(injection_times)).or..not.all(ieee_is_finite(injection_scales))) &
        error stop 'Nonfinite injection schedule'
      if(any(injection_times<0).or.any(injection_times>1e6_dp).or.any(injection_scales<0).or. &
        any(injection_scales>1e6_dp)) error stop 'Invalid injection schedule'
      do q=2,injection_events
        if(injection_times(q)<injection_times(q-1)) error stop 'Injection times must be nondecreasing'
      end do
      event_schedule%times=injection_times(:injection_events)
      event_schedule%ends=injection_end_times(:injection_events)
      event_schedule%ports=event_ports(:injection_events)
      event_schedule%scales=injection_scales(:injection_events)
    end if
    if(any(mixture_density_models(nc+1:)/='constant')) error stop 'Unused mixture property entries'
    if(any(mixture_cp_models(nc+1:)/='constant')) error stop 'Unused mixture property entries'
    if(any(mixture_conductivity_models(nc+1:)/='supplied')) error stop 'Unused mixture conductivity entries'
    if(any(mixture_viscosity_models(nc+1:)/='supplied')) error stop 'Unused mixture viscosity entries'
    if(any(mixture_conductivities(nc+1:)/=0)) error stop 'Unused mixture conductivity entries'
    if(any(mixture_viscosities(nc+1:)/=0)) error stop 'Unused mixture viscosity entries'
    if(any(mixture_conductivity_coefficients(:,nc+1:)/=0)) error stop 'Unused mixture conductivity entries'
    if(any(mixture_viscosity_coefficients(:,nc+1:)/=0)) error stop 'Unused mixture viscosity entries'
    if(any(mixture_conductivity_temperature_ranges(:,nc+1:)/=0)) error stop 'Unused mixture conductivity entries'
    if(any(mixture_viscosity_temperature_ranges(:,nc+1:)/=0)) error stop 'Unused mixture viscosity entries'
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
      liquids(c)%conductivity_model=mixture_conductivity_models(c)
      liquids(c)%viscosity_model=mixture_viscosity_models(c)
      liquids(c)%conductivity=mixture_conductivities(c)
      liquids(c)%viscosity=mixture_viscosities(c)
      liquids(c)%conductivity_coefficients=mixture_conductivity_coefficients(:,c)
      liquids(c)%viscosity_coefficients=mixture_viscosity_coefficients(:,c)
      liquids(c)%conductivity_temperature_range=mixture_conductivity_temperature_ranges(:,c)
      liquids(c)%viscosity_temperature_range=mixture_viscosity_temperature_ranges(:,c)
      liquids(c)%boiling_model=mixture_boiling_models(c)
      liquids(c)%critical_temperature=mixture_critical_temperatures(c)
      liquids(c)%normal_boiling_temperature=mixture_normal_boiling_temperatures(c)
      liquids(c)%molecular_weight=mixture_molecular_weights(c)
      if(.not.valid_spray_liquid(liquids(c))) error stop 'Invalid mixture component property model'
    end do
    if(all(liquids%viscosity_model/='supplied')) then
      call apply_mixture_breakup_viscosity(liquids,fractions(:nc),liquid_temperature,particle_options,ok, &
        initial_pressure,species(liquids%vapor_index)%molecular_weight)
      if(.not.ok) error stop 'Invalid component viscosity at injection'
    end if
    if(.not.valid_spray_options(particle_options)) error stop 'Unsupported particle options'
    if(wall_impacts) then
      if(particle_options%liquid_viscosity<=0.or.particle_options%surface_tension<=0) &
        error stop 'EB impacts require viscosity and surface tension'
    end if
    if(port_count<1.or.port_count>16) error stop 'Invalid injection port count'
    if(any(event_schedule%ports<1).or.any(event_schedule%ports>port_count)) &
      error stop 'Invalid injection event port'
    allocate(injection_ports(port_count))
    maximum_injection_speed=0
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
        denominator=rho*acos(-1._dp)*port_radii(port)**2*port_discharge_coefficients(port)
        maximum_injection_speed=max(maximum_injection_speed,maxval(rate_profiles(port)%rates)/denominator)
        template_speed=1
      else if(port_discharge_coefficients(port)/=1) then
        error stop 'Discharge coefficient requires rate-driven velocity'
      end if
      if(.not.port_velocity_from_rate(port)) &
        maximum_injection_speed=max(maximum_injection_speed,port_speeds(port))
      call inject_cone_parcels(effective,port_masses(port),port_numbers(port),1_int64, &
        port_centers(:,port),port_axes(:,port),template_speed,port_angles(port),port_diameters(port), &
        port_temperatures(port),injection_ports(port)%parcels,ok, &
        diameter_shape=port_sizes(port)%shape,diameter_min=port_sizes(port)%minimum, &
        diameter_max=port_sizes(port)%maximum,diameter_weighting=port_sizes(port)%weighting, &
        diameter_distribution=port_sizes(port)%distribution,diameter_std=port_sizes(port)%std, &
        nozzle_radius=port_radii(port),nozzle_model=port_nozzles(port), &
        reference_breakup=particle_options%etab_parcel_factor>=0.or.particle_options%khrt_parcel_factor>=0)
      if(.not.ok) error stop 'Invalid injection port template'
      injection_ports(port)%fractions=port_fractions(:nc,port)
      injection_ports(port)%parcels%material_id=port
    end do
    event_factors=event_schedule%scales
    allocate(event_velocity_scales(injection_events));event_velocity_scales=1
    do port=1,port_count
      q=port_profile_points(port)
      if(q/=0.and.q<2) error stop 'Invalid port profile size'
      if(q==0) cycle
      call integrate_spray_injection_profile(rate_profiles(port)%times,rate_profiles(port)%rates, &
        rate_profiles(port)%times(1),rate_profiles(port)%times(q),profile_mass,ok)
      if(.not.ok) error stop 'Invalid port rate profile'
    end do
    do c=1,injection_events
      port=event_schedule%ports(c);q=port_profile_points(port)
      if(q==0) cycle
      call integrate_spray_injection_profile(rate_profiles(port)%times,rate_profiles(port)%rates, &
        event_schedule%times(c),event_schedule%ends(c),profile_mass,ok)
      if(.not.ok) error stop 'Invalid injection integration interval'
      profile_mass=profile_mass*event_schedule%scales(c)
      if(profile_mass>1e6_dp.or.profile_mass>1e6_dp*port_masses(port)) error stop 'Integrated injection exceeds bounds'
      event_factors(c)=profile_mass/port_masses(port)
      if(port_velocity_from_rate(port)) then
        event_velocity_scales(c)=0
        if(profile_mass>0) then
          call liquid_mixture_properties(liquids,port_fractions(:nc,port),port_temperatures(port),rho,cp,e,ok)
          if(.not.ok) error stop 'Invalid rate-driven mixture density'
          denominator=(event_schedule%ends(c)-event_schedule%times(c))*rho*acos(-1._dp)* &
            port_radii(port)**2*port_discharge_coefficients(port)
          if(denominator<=0) error stop 'Invalid rate-driven injection interval'
          if(profile_mass>1e7_dp*denominator) error stop 'Rate-driven injection speed exceeds bounds'
          event_velocity_scales(c)=profile_mass/denominator
          maximum_injection_speed=max(maximum_injection_speed,event_velocity_scales(c), &
            maxval(rate_profiles(port)%rates)*event_schedule%scales(c)/ &
            (rho*acos(-1._dp)*port_radii(port)**2*port_discharge_coefficients(port)))
        end if
      end if
    end do
    allocate(drops(0),component_mass(nc,0),films(0),film_mass(nc,0))
    if(checkpoint_interval_steps>0.or.len_trim(restart_file)>0) then
      open(newunit=stream,status='scratch',form='formatted',iostat=code)
      if(code/=0) error stop 'Cannot stage spray context'
      write(stream,nml=eb_spray,iostat=code)
      if(code==0.and.les_enabled) write(stream,nml=eb_spray_les,iostat=code)
      if(code==0) write(stream,'(es27.18e3)',iostat=code) nodal_phi
      if(code/=0) error stop 'Cannot serialize spray context'
      rewind(stream)
      spray_context=mixture_evolution_id//'|'
      spray_context=spray_context//'PARTICLE_TIMESTEP_V1|'
      if(allocated(frame)) spray_context=spray_context//'ROTATING_PARTICLE_ENERGY_V1|'
      if(time_integrator=='sdc') spray_context=spray_context//'POPULATION_PRIMITIVE_V1|POPULATION_NATIVE_CFL_V1|'
      if(particle_options%etab_parcel_factor>=0.or.particle_options%khrt_parcel_factor>=0) &
        spray_context=spray_context//'REF_BREAKUP_N0_V1|'
      if(wall_impacts) spray_context=spray_context//'SPLASH_N0_V2|'
      if(allocated(constant_coefficients)) spray_context=spray_context//'CONSTANT_PARTICLE_TRANSPORT_V1|'
      if(particle_options%breakup_model=='etab') spray_context=spray_context//'ETAB_MINMASS_V2|'
      if(particle_options%khrt_parcel_factor>=0) spray_context=spray_context//'KHRT_MINMASS_V2|'
      if(time_integrator=='sdc'.and..not.subcycling) spray_context=spray_context//'COMMON_SDC_CELL_SOURCE_V1|'
      if(time_integrator=='sdc'.and.subcycling) &
        spray_context=spray_context//'RECURSIVE_SDC_PHYSICAL_SYNC_V1|'
      do
        read(stream,'(a)',iostat=code) context_line
        if(code==iostat_end) exit
        if(code/=0) error stop 'Cannot read spray context'
        if(particle_options%khrt_parcel_factor==-1.and.index(context_line,'%KHRT_PARCEL_FACTOR=')>0) cycle
        if(particle_options%mixture_time_scheme=='conservative'.and. &
            index(context_line,'%MIXTURE_TIME_SCHEME=')>0) cycle
        ! Use the same bounded context capacity as profiles and the reader.
        if(len(spray_context)+len_trim(context_line)+1>16000000) error stop 'Spray context capacity exceeded'
        spray_context=spray_context//trim(context_line)//'|'
      end do
      close(stream)
      do port=1,port_count
        q=port_profile_points(port)
        if(q<=64) cycle
        allocate(character(len=56*q)::profile_context)
        do c=1,q
          write(profile_context(56*(c-1)+1:56*c),'(2(es27.18e3,1x))') &
            rate_profiles(port)%times(c),rate_profiles(port)%rates(c)
        end do
        write(context_line,'(a,i0,a,i0,a)') 'LONG_RATE_V1:',port,':',q,':'
        if(len(spray_context)+len_trim(context_line)+len(profile_context)+1>16000000) &
          error stop 'Spray context capacity exceeded'
        spray_context=spray_context//trim(context_line)//profile_context//'|'
        deallocate(profile_context)
      end do
      if(len_trim(event_file)>0) then
        allocate(character(len=96*injection_events)::profile_context)
        do c=1,injection_events
          write(profile_context(96*(c-1)+1:96*c),'(2(es27.18e3,1x),i11,1x,es27.18e3,1x)') &
            event_schedule%times(c),event_schedule%ends(c),event_schedule%ports(c),event_schedule%scales(c)
        end do
        if(len(spray_context)+len(profile_context)+16>16000000) error stop 'Spray context capacity exceeded'
        spray_context=spray_context//'EVENT_FILE_V1:'//profile_context//'|'
        deallocate(profile_context)
      end if
    end if
  end subroutine

  subroutine inject_due_events()
    type(spray_parcel) :: injected_parcel
    type(spray_parcel), allocatable :: candidate(:)
    real(dp), allocatable :: masses(:,:)
    real(dp) :: change(nv),totals(5),rho,cp,e,scale,potential
    type(spray_liquid) :: effective
    integer :: q,c,number,old_count,local_count,global_count
    do while(next_injection<=injection_events)
      if(event_schedule%times(next_injection)>record%time) exit
      scale=event_factors(next_injection)
      if(scale>0) then
        associate(injection_template=>injection_ports(event_schedule%ports(next_injection))%parcels, &
          injection_fractions=>injection_ports(event_schedule%ports(next_injection))%fractions)
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
          if(allocated(frame)) then
            call rotating_potential_3d(frame,injected_parcel%position,potential,ok)
            if(.not.ok) error stop 'Invalid rotating injection position'
            totals(5)=totals(5)-totals(1)*potential
          end if
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
            candidate(old_count+q)%khrt_initial_count=scale*injection_template(q)%khrt_initial_count
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
    real(dp) :: local(size(inventory)),global(size(inventory)),rho,cp,e,amount,potential
    integer :: q,c
    local=0
    do q=1,size(particles)
      if(particles(q)%mass<=0) cycle
      call liquid_mixture_properties(liquids,components(:,q)/particles(q)%mass,particles(q)%temperature,rho,cp,e,ok)
      if(.not.ok) error stop 'Invalid particle inventory'
      amount=particles(q)%mass*particles(q)%multiplicity
      local(1)=local(1)+amount;local(2:4)=local(2:4)+amount*particles(q)%velocity
      potential=0
      if(allocated(frame)) then
        call rotating_potential_3d(frame,particles(q)%position,potential,ok)
        if(.not.ok) error stop 'Invalid rotating particle inventory'
      end if
      local(5)=local(5)+amount*(e+.5_dp*sum(particles(q)%velocity**2)-potential)
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
    real(dp) :: local(size(inventory)),rho,cp,e,potential
    integer :: q,c
    local=0
    do q=1,size(patches)
      call liquid_mixture_properties(liquids,components(:,q)/patches(q)%mass,patches(q)%temperature,rho,cp,e,ok)
      if(.not.ok) error stop 'Invalid film inventory'
      potential=0
      if(allocated(frame)) then
        call rotating_potential_3d(frame,patches(q)%position,potential,ok)
        if(.not.ok) error stop 'Invalid rotating film inventory'
      end if
      local(1)=local(1)+patches(q)%mass;local(5)=local(5)+patches(q)%mass*(e-potential)
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
    real(dp) :: source(nv),expected_source(nv),inventory(nv),saved_time,saved_force(nv),saved_coordinate
    real(dp) :: saved_numerical,saved_pressure_work
    integer :: stream,code,closed,saved_steps,saved_event,expected_event,date(8)
    integer(int64) :: watermark,expected_mark,tick,expected_ids,expected_splits
    character(len=256) :: token,step_label
    character(len=:), allocatable :: path,bound_context,line
    logical :: valid,agreed,exists
    allocate(character(len=max(131072,len(spray_context)+8192))::line)
    candidate=mpi_state;ct=mpi_temperature;coupled=record;gas_record=record
    checkpoint_tree=topology;checkpoint_geometry=mpi_geometry;checkpoint_owners=mpi_owners
    particles=drops;components=component_mass;patches=films;patch_mass=film_mass
    saved_numerical=0;saved_pressure_work=0
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
      if((any(acceleration/=0).or.time_integrator=='sdc').and.code==0) write(stream,'(*(es27.18e3,1x))',iostat=code) force_total
      if(allocated(frame).and.code==0) write(stream,'(es27.18e3)',iostat=code) coordinate_total(5)
      if(allocated(held_energy).and.code==0) then
        write(stream,'(a)',iostat=code) 'PELEF_MPI_EB_SPRAY_HELD_LEDGER_V1'
        if(code==0) write(stream,'(2(es27.18e3,1x))',iostat=code) numerical_energy,internal_energy_pressure_work
      endif
      if(code/=0) error stop 'Cannot write coupled checkpoint prefix'
      ! The gas archive is a snapshot with its own zero-exchange reference.
      ! The complete coupled history is retained separately in the prefix above.
      call integrate_mpi_eb_tree_3d(checkpoint_tree,checkpoint_owners,checkpoint_geometry,candidate,inventory, &
      MPI_COMM_WORLD,valid)
      if(.not.valid) error stop 'Cannot integrate gas snapshot'
      gas_record%initial=inventory;gas_record%scale=abs(inventory)
      gas_record%boundary=0;gas_record%wall=0;gas_record%chemical=0
      ! Snapshot initial inventory is already corrected. Its gas-only prefix
      ! starts at zero; the coupled history above retains the cumulative reset.
      if(allocated(held_energy)) gas_record%energy_correction=0
      if(time_integrator=='sdc') gas_record%forcing=0._dp*inventory
      call write_mpi_eb_tree_checkpoint_3d(path//'.gas',stream,species,reactions,transport,config,topology, &
        mpi_owners,mpi_geometry,candidate,ct,gas_record,MPI_COMM_WORLD,valid, &
        bundle_identity,checkpoint_integrator,base_mole_fractions,element_names,element_atoms,frame)
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
      if(any(acceleration/=0).or.time_integrator=='sdc') then
        read(stream,*,iostat=code) saved_force
        if(code/=0) error stop 'Missing external force history'
        if(.not.all(ieee_is_finite(saved_force))) error stop 'Invalid external force history'
      end if
      saved_coordinate=0
      if(allocated(frame)) then
        read(stream,*,iostat=code) saved_coordinate
        if(code/=0) error stop 'Missing coordinate energy history'
        if(.not.ieee_is_finite(saved_coordinate)) error stop 'Invalid coordinate energy history'
      endif
      if(allocated(held_energy)) then
        read(stream,'(a)',iostat=code) line
        if(code/=0.or.line/='PELEF_MPI_EB_SPRAY_HELD_LEDGER_V1') error stop 'Missing coupled held-energy ledger'
        read(stream,*,iostat=code) saved_numerical,saved_pressure_work
        if(code/=0) error stop 'Invalid coupled held-energy ledger'
        if(.not.all(ieee_is_finite([saved_numerical,saved_pressure_work]))) &
          error stop 'Nonfinite coupled held-energy ledger'
      endif
      call read_mpi_eb_tree_checkpoint_3d(path//'.gas',stream,species,reactions,transport,config, &
        checkpoint_tree,checkpoint_owners, &
        checkpoint_geometry,candidate,ct,gas_record,MPI_COMM_WORLD,valid, &
        bundle_identity,checkpoint_integrator,base_mole_fractions,element_names,element_atoms,frame)
      close(stream,iostat=closed)
      if(.not.valid.or.closed/=0) error stop 'Gas snapshot restore failed'
    end if
    if(wall_impacts) then
    call mpi_eb_particle_checkpoint(MPI_COMM_WORLD,writing,path//'.particles',bound_context,checkpoint_tree,checkpoint_owners, &
      nodal_phi,domain_lower,(domain_upper-domain_lower)/real(root_shape,dp), &
      min(125000,particle_options%maximum_parcels),particles,components,saved_time,saved_steps,saved_event, &
      watermark,source,valid,patches,patch_mass,walls%maximum_films,periodic=boundary_faces([1,3,5])=='periodic')
    else
    call mpi_eb_particle_checkpoint(MPI_COMM_WORLD,writing,path//'.particles',bound_context,checkpoint_tree,checkpoint_owners, &
      nodal_phi,domain_lower,(domain_upper-domain_lower)/real(root_shape,dp), &
      min(125000,particle_options%maximum_parcels),particles,components,saved_time,saved_steps,saved_event, &
      watermark,source,valid,periodic=boundary_faces([1,3,5])=='periodic')
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
      if(automatic_timestep.and.saved_steps>0.and.previous_timestep<=0) &
        error stop 'Restart accepted timestep missing'
      if(saved_time>final_time.or.saved_steps>maximum_steps) error stop 'Restart beyond requested limits'
      if(saved_event<1.or.saved_event>injection_events+1) error stop 'Invalid restart event'
      if(any(event_schedule%times(:saved_event-1)>saved_time)) error stop 'Processed future injection'
      if(saved_event<=injection_events) then
        if(event_schedule%times(saved_event)<saved_time) error stop 'Missed restart injection'
      end if
      expected_ids=0
      do p=1,saved_event-1
        if(event_factors(p)>0) expected_ids=expected_ids+size(injection_ports(event_schedule%ports(p))%parcels)
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
      if(allocated(held_energy)) then
        if(.not.allocated(gas_record%internal_energy%fields)) error stop 'Held spray restart lacks UEINT fields'
        if(.not.mol_mode.and..not.gas_record%has_internal_energy_pressure_work) &
          error stop 'Held spray restart lacks pressure-work history'
        if(gas_record%internal_energy_pressure_work/=saved_pressure_work) &
          error stop 'Held spray pressure-work history mismatch'
        if(gas_record%energy_correction/=0) error stop 'Held spray gas snapshot has a nonzero reset prefix'
        coupled%internal_energy=gas_record%internal_energy;coupled%energy_policy=gas_record%energy_policy
        coupled%energy_correction=saved_numerical;coupled%internal_energy_pressure_work=saved_pressure_work
        coupled%has_internal_energy_pressure_work=gas_record%has_internal_energy_pressure_work
      endif
      call integrate_mpi_eb_tree_3d(checkpoint_tree,checkpoint_owners,checkpoint_geometry,candidate,inventory, &
      MPI_COMM_WORLD,valid)
      if(.not.valid) error stop 'Cannot integrate restored gas'
      call add_particle_inventory(inventory,particles,components)
      call add_film_inventory(inventory,patches,patch_mass)
      inventory(5)=inventory(5)-saved_numerical
      if(any(abs(inventory-coupled%initial-coupled%boundary-coupled%wall-coupled%chemical-source-saved_force)> &
        5e-10_dp*max(1._dp,coupled%scale,abs(inventory), &
        abs(coupled%boundary)+abs(coupled%wall)+abs(coupled%chemical)+abs(source)+abs(saved_force)))) &
        error stop 'Restored budget failed'
      if(time_integrator=='sdc') then
        coupled%forcing=saved_force
        if(allocated(sdc_ir)) then
          sdc_ir=gas_record%reaction_source;coupled%reaction_source=sdc_ir
        end if
      end if
      mpi_state=candidate;mpi_temperature=ct;record=coupled
      if(allocated(held_energy)) then
        held_energy=coupled%internal_energy;numerical_energy=saved_numerical
        internal_energy_pressure_work=saved_pressure_work;numerical_total(5)=saved_numerical
      endif
      topology=checkpoint_tree;mpi_owners=checkpoint_owners;mpi_geometry=checkpoint_geometry
      regrid_stage=restart_layout_stage;tree%topology=topology
      if(time_integrator=='sdc') then
        if(allocated(gas_record%heat_release%fields)) then
          sdc_heat=gas_record%heat_release;record%heat_release=sdc_heat
        else
          if(allocated(record%heat_release%fields)) deallocate(record%heat_release%fields)
          call initialize_heat_snapshot()
        end if
      end if
      if(allocated(tree%nodes)) deallocate(tree%nodes)
      allocate(tree%nodes(size(topology%nodes)))
      do p=1,size(tree%nodes)
        tree%nodes(p)%geometry=mpi_geometry(p)
      end do
      call move_alloc(particles,drops);call move_alloc(components,component_mass)
      call move_alloc(patches,films);call move_alloc(patch_mass,film_mass)
      force_total=saved_force;coordinate_total=0;coordinate_total(5)=saved_coordinate
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
    if(mpi_rank==mpi_owners(1)) call rotate_seed(seed_u,coarse,valid)
    call MPI_Allreduce(valid,all_valid,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS.or..not.all_valid) error stop 'Initial rotating energy failed'
    call initialize_mpi_eb_tree_3d(species,topology,mpi_owners,geometry_levels(topology%level_count()),seed_u,seed_t, &
      mpi_geometry,mpi_state,mpi_temperature,MPI_COMM_WORLD,valid,frame,held_energy)
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
  ! Ordered full-node records, including covered coarse donors. Commit only
  ! after all records pass admission; MPI reads one authoritative byte stream.
  subroutine import_initial_state()
    type(amr_tree_state_3d) :: candidate,ct,ch
    integer :: stream,code,p,c,i,j,k,nn,ns,idx(4),dims(3),child,lo(3),hi(3)
    real(dp) :: values(nv+2),q(nv),temperature_value,sound,kinetic,scale
    character(len=8192) :: line
    character(len=64) :: name
    logical :: reader,owned,valid
    logical, allocatable :: covered(:,:,:)
    reader=.true.
#ifdef PELEF_MPI_EB_TREE
    reader=mpi_rank==0
#endif
    stream=0;code=0
    if(reader) open(newunit=stream,file=trim(initial_state_file),status='old',action='read',iostat=code)
    call initial_read_status(code)
    if(code/=0) error stop 'Cannot open initial state file'
    call initial_read_line(stream,line,code)
    if(code/=0.or.trim(line)/='PELEF_EB_TREE_STATE_V1') error stop 'Invalid initial state header'
    call initial_read_line(stream,line,code)
    if(code/=0) error stop 'Missing initial state dimensions'
    read(line,*,iostat=code) nn,ns
    if(code/=0) error stop 'Invalid initial state dimensions'
    if(nn/=size(topology%nodes).or.ns/=size(species)) error stop 'Initial state dimensions mismatch'
    do c=1,ns
      call initial_read_line(stream,line,code)
      if(code/=0) error stop 'Missing initial state species'
      read(line,*,iostat=code) name
      if(code/=0) error stop 'Invalid initial state species'
      if(trim(name)/=trim(species(c)%name)) error stop 'Initial state species order mismatch'
    enddo
    allocate(candidate%fields(nn),ct%fields(nn),ch%fields(nn))
    do p=1,nn
      call initial_read_line(stream,line,code)
      if(code/=0) error stop 'Missing initial node shape'
      read(line,*,iostat=code) idx(1),dims
      if(code/=0) error stop 'Invalid initial node shape'
      if(idx(1)/=p.or.any(dims/=topology%nodes(p)%shape)) error stop 'Initial node shape mismatch'
      allocate(covered(dims(1),dims(2),dims(3)))
      covered=.false.
      do child=p+1,nn
        if(topology%nodes(child)%parent/=p) cycle
        lo=[topology%nodes(child)%patch%coarse_i_lower,topology%nodes(child)%patch%coarse_j_lower, &
          topology%nodes(child)%patch%coarse_k_lower]
        hi=[topology%nodes(child)%patch%coarse_i_upper,topology%nodes(child)%patch%coarse_j_upper, &
          topology%nodes(child)%patch%coarse_k_upper]
        covered(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.true.
      enddo
      owned=.true.
#ifdef PELEF_MPI_EB_TREE
      owned=mpi_owners(p)==mpi_rank
#endif
      if(owned) then
        allocate(candidate%fields(p)%values(nv,dims(1),dims(2),dims(3)))
        allocate(ct%fields(p)%values(1,dims(1),dims(2),dims(3)))
        allocate(ch%fields(p)%values(1,dims(1),dims(2),dims(3)))
      endif
      do k=1,dims(3)
      do j=1,dims(2)
      do i=1,dims(1)
        call initial_read_line(stream,line,code)
        if(code/=0) error stop 'Missing initial cell'
        read(line,*,iostat=code) idx,values
        if(code/=0) error stop 'Invalid initial cell'
        if(any(idx/=[p,i,j,k])) error stop 'Initial cell ordering mismatch'
        if(.not.all(ieee_is_finite(values))) error stop 'Nonfinite initial cell'
        if(values(1)<0.or.any(values(6:nv)<0)) error stop 'Negative initial density'
        if(values(1)>0) then
          if(values(nv+2)<=0) error stop 'Invalid initial temperature'
          if(abs(sum(values(6:nv))-values(1))>5.e-12_dp*values(1)) &
            error stop 'Initial species mass mismatch'
          kinetic=.5_dp*sum(values(2:4)**2)/values(1)
          scale=max(1._dp,abs(values(5)),abs(values(nv+1)),abs(kinetic))
          if(covered(i,j,k)) then
            ! Independent restriction of E and H retains nonnegative unresolved
            ! kinetic energy: average(KE) - KE(average(rho,momentum)).
            if(values(5)-kinetic-values(nv+1)<-1.e-10_dp*scale) &
              error stop 'Initial covered energy deficit'
          else
            if(abs(values(5)-kinetic-values(nv+1))>1.e-10_dp*scale) &
              error stop 'Initial held energy mismatch'
          endif
          call reactive_conserved_to_primitive(species,values(:nv),values(nv+2),q, &
            temperature_value,sound,valid,internal_energy_density=values(nv+1))
          if(.not.valid) error stop 'Invalid initial cell EOS'
          ! A covered donor can carry independently restricted temperature.
          ! Preserve that inverse-EOS seed; only leaves require EOS identity.
          if(.not.covered(i,j,k)) then
            if(abs(temperature_value-values(nv+2))>1.e-8_dp*max(1._dp,abs(temperature_value))) &
              error stop 'Initial temperature mismatch'
          endif
        else
          if(any(values(:nv+1)/=0).or.values(nv+2)<0) error stop 'Invalid empty initial cell'
        endif
        if(owned) then
          if(tree%nodes(p)%geometry%volume_fraction(i,j,k)>0.and.values(1)<=0) &
            error stop 'Empty initial fluid cell'
          candidate%fields(p)%values(:,i,j,k)=values(:nv)
          ch%fields(p)%values(1,i,j,k)=values(nv+1)
          ct%fields(p)%values(1,i,j,k)=values(nv+2)
        endif
      enddo
      enddo
      enddo
      deallocate(covered)
    enddo
    call initial_read_line(stream,line,code)
    if(code>=0) error stop 'Unexpected initial state trailing record'
    code=0
    if(reader) close(stream,iostat=code)
    call initial_read_status(code)
    if(code/=0) error stop 'Cannot close initial state file'
#ifdef PELEF_MPI_EB_TREE
    mpi_state=candidate;mpi_temperature=ct
#else
    do p=1,nn
      tree%nodes(p)%state=candidate%fields(p)%values
      tree%nodes(p)%temperature=ct%fields(p)%values(1,:,:,:)
    enddo
#endif
    held_energy=ch
  end subroutine import_initial_state

  subroutine initial_read_status(code)
    integer, intent(inout) :: code
#ifdef PELEF_MPI_EB_TREE
    integer :: ierr
    call MPI_Bcast(code,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'Initial state status broadcast failed'
#endif
  end subroutine initial_read_status

  subroutine initial_read_line(stream,line,code)
    integer, intent(in) :: stream
    character(len=*), intent(out) :: line
    integer, intent(out) :: code
#ifdef PELEF_MPI_EB_TREE
    integer :: ierr
#endif
    line='';code=0
#ifdef PELEF_MPI_EB_TREE
    if(mpi_rank==0) then
#endif
      read(stream,'(a)',iostat=code) line
#ifdef PELEF_MPI_EB_TREE
    endif
#endif
    call initial_read_status(code)
#ifdef PELEF_MPI_EB_TREE
    if(code==0) then
      call MPI_Bcast(line,len(line),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      if(ierr/=MPI_SUCCESS) error stop 'Initial state row broadcast failed'
    endif
#endif
  end subroutine initial_read_line

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
  subroutine initialize_heat_snapshot()
    integer :: p,n(3)
    if(allocated(sdc_heat)) deallocate(sdc_heat)
    allocate(sdc_heat)
    allocate(sdc_heat%fields(size(topology%nodes)))
    do p=1,size(topology%nodes)
#ifdef PELEF_MPI_EB_TREE
      if(mpi_owners(p)/=mpi_rank) cycle
#endif
      n=topology%nodes(p)%shape
      allocate(sdc_heat%fields(p)%values(1,n(1),n(2),n(3)),source=0._dp)
    end do
  end subroutine initialize_heat_snapshot
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
      if(allocated(held_energy)) record%internal_energy=held_energy
      return
    end if
    if (distributed.and.remapping) then
      level=max(topology%level_count(),layouts(stage)%level_count())
      call partition_amr_tree_3d(layouts(stage),mpi_size,0,next_owners,valid)
      if (.not.valid) error stop 'MPI regrid partition failed'
#ifdef PELEF_MPI_EB_SPRAY
      call remap_mpi_eb_spray(MPI_COMM_WORLD,species,topology,mpi_owners,mpi_geometry,mpi_state,mpi_temperature, &
        drops,component_mass,layouts(stage),next_owners,geometry_levels(level),nodal_phi,domain_lower, &
        (domain_upper-domain_lower)/real(root_shape,dp),particle_options%maximum_parcels,valid, &
        films,film_mass,walls%maximum_films,sdc_ir,periodic=boundary_faces([1,3,5])=='periodic', &
        frame=frame,heat_release=sdc_heat)
      if(valid.and.allocated(sdc_ir)) record%reaction_source=sdc_ir
      if(valid.and.allocated(sdc_heat)) record%heat_release=sdc_heat
#else
      call remap_mpi_eb_tree_3d(species,topology,layouts(stage),mpi_owners,next_owners,geometry_levels(level), &
        mpi_geometry,mpi_state,mpi_temperature,MPI_COMM_WORLD,valid,sdc_ir,frame,sdc_heat,held_energy)
      if(valid.and.allocated(sdc_ir)) record%reaction_source=sdc_ir
      if(valid.and.allocated(sdc_heat)) record%heat_release=sdc_heat
#endif
      if (.not.valid) error stop 'MPI conservative regrid failed'
      if(allocated(held_energy)) record%internal_energy=held_energy
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
      call remap_reactive_eb_tree_3d(species,layouts(stage),geometry_levels(level),tree,valid,sdc_ir,frame,sdc_heat,held_energy)
      if(valid.and.allocated(sdc_ir)) record%reaction_source=sdc_ir
      if(valid.and.allocated(sdc_heat)) record%heat_release=sdc_heat
    else
      level=layouts(stage)%level_count()
      call initialize_reactive_eb_tree_3d(species,layouts(stage),geometry_levels(level),u,t,tree,valid,frame,held_energy)
    end if
    if (.not.valid) error stop 'Cannot activate scheduled tree'
    if(allocated(held_energy)) record%internal_energy=held_energy
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
      extended_geometry,valid,halo_width=ctu_halo)
    if(.not.valid) error stop 'Cannot build CTU patch geometry cache'
    if(allocated(active_embedded_temperature).or.embedded_noslip) then
      if(allocated(heat_geometry)) deallocate(heat_geometry)
      allocate(heat_geometry(size(topology%nodes)))
      do p=1,size(topology%nodes)
        call extract_eb_geometry_3d(extended_geometry(p),[ctu_halo-1,ctu_halo-1,ctu_halo-1], &
          topology%nodes(p)%shape+ctu_halo+2,heat_geometry(p),valid)
        if(.not.valid) error stop 'Cannot build embedded heat geometry cache'
      end do
    end if
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
      write(*,'(a,*(es24.16,1x))') 'Rejected chemical initial species masses: ',initial_values(6:)
      write(*,'(a,*(es24.16,1x))') 'Rejected chemical species changes: ',changes(6:)
      do e=1,element_count
        write(*,'(a,1x,2(es24.16,1x))') trim(element_names(e)),element_change(e),element_scale(e)
      end do
    end if
  end function chemical_elements_valid
#endif

  subroutine balance()
    integer :: component
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
    error=maxval(abs(current-record%initial-record%boundary-record%wall-record%chemical- &
      injected-force_total-coordinate_total-numerical_total)/ &
      max(1._dp,record%scale,abs(current),abs(record%boundary)+abs(record%wall)+ &
      abs(record%chemical)+abs(injected)+abs(force_total)+abs(coordinate_total)+abs(numerical_total)))
    if (.not.ieee_is_finite(error).or.error>5.e-10_dp) then
      do component=1,size(current)
        write(*,'(a,i0,10(1x,es24.16))') 'Rejected composite component: ',component,current(component), &
          record%initial(component),record%boundary(component),record%wall(component),record%chemical(component), &
          injected(component),force_total(component),coordinate_total(component),numerical_total(component), &
          current(component)-record%initial(component)- &
          record%boundary(component)-record%wall(component)-record%chemical(component)-injected(component)- &
          force_total(component)-coordinate_total(component)-numerical_total(component)
      enddo
      error stop 'Composite budget failed'
    endif
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
  subroutine rotate_seed(state,grid,valid)
    real(dp), intent(inout) :: state(:,:,:,:)
    type(eb_geometry_3d), intent(in) :: grid
    logical, intent(out) :: valid
    real(dp) :: position(3),energy
    integer :: x,y,z
    valid=.true.
    if(.not.allocated(frame)) return
    valid=.false.
    if(any(shape(state)/=[nv,grid%nx,grid%ny,grid%nz])) return
    do z=1,grid%nz
    do y=1,grid%ny
    do x=1,grid%nx
      if(grid%volume_fraction(x,y,z)<=0) cycle
      position=[grid%x_lower,grid%y_lower,grid%z_lower]+ &
        (real([x,y,z],dp)-.5_dp)*[grid%dx,grid%dy,grid%dz]
      energy=state(5,x,y,z)
      call convert_rotating_energy_3d(frame,position,state(1,x,y,z),state(5,x,y,z),.true.,energy,valid)
      if(.not.valid) return
      state(5,x,y,z)=energy
    end do
    end do
    end do
    valid=.true.
  end subroutine rotate_seed
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
    real(dp), allocatable :: stored_constants(:)
    integer :: stored_constant_count
    character(len=64) :: magic,integrator
    logical :: automatic
#if (defined(PELEF_MPI_EB_LES) || defined(PELEF_SERIAL_EB_LES)) && !defined(PELEF_MPI_EB_SPRAY)
    character(len=32) :: stored_les_model,stored_les_filter
    integer :: stored_les_ratios(2)
    real(dp) :: stored_les_controls(5)
    logical :: stored_les_harmonic
#endif
    character(len=16) :: stored_faces(6)
    integer :: stored_counts(6),stored_nspecies
    integer :: stored_rotation_axis
    real(dp) :: stored_rotation(4)
    real(dp) :: stored_composition(size(species),6)
    real(dp) :: stored_prescribed(30),stored_wall_temperature(6),stored_embedded_temperature
    real(dp) :: interval,stored_hydro_cfl,stored_transport_cfl,stored_inlet(18),stored_timestep_controls(3)
    if (reading) then
      read(stream,'(a)',iostat=ios) magic
      if (ios/=0) error stop 'Missing tree application header'
      if(magic=='PELEF_EB_TREE_CONNECTED_DIRECTIONS_V1') then
        if(.not.connected_directions) error stop 'Restart connected-direction policy mismatch'
        read(stream,'(a)',iostat=ios) magic
        if(ios/=0) error stop 'Missing connected-direction-policy application header'
      else if(connected_directions) then
        error stop 'Restart connected-direction policy mismatch'
      endif
      if(magic=='PELEF_EB_TREE_PAIRED_PRESSURE_V1') then
        if(.not.paired_pressure) error stop 'Restart pressure policy mismatch'
        read(stream,'(a)',iostat=ios) magic
        if(ios/=0) error stop 'Missing pressure-policy application header'
      else if(paired_pressure) then
        error stop 'Restart pressure policy mismatch'
      endif
      if(magic=='PELEF_EB_TREE_CONSISTENT_VELOCITY_V1') then
        if(.not.consistent_velocity) error stop 'Restart velocity policy mismatch'
        read(stream,'(a)',iostat=ios) magic
        if(ios/=0) error stop 'Missing velocity-policy application header'
      else if(consistent_velocity) then
        error stop 'Restart velocity policy mismatch'
      endif
      if(magic=='PELEF_EB_TREE_SEPARATE_MOMENTUM_V1') then
        if(.not.separate_momentum) error stop 'Restart momentum policy mismatch'
        read(stream,'(a)',iostat=ios) magic
        if(ios/=0) error stop 'Missing momentum-policy application header'
      else if(separate_momentum) then
        error stop 'Restart momentum policy mismatch'
      endif
      if(magic=='PELEF_EB_TREE_CONSTANT_TRANSPORT_V1') then
        if(.not.allocated(constant_coefficients)) error stop 'Restart constant transport enable mismatch'
        read(stream,*,iostat=ios) stored_constant_count
        if(ios/=0) error stop 'Missing restart constant transport size'
        if(stored_constant_count/=size(constant_coefficients)) error stop 'Restart constant transport size mismatch'
        allocate(stored_constants(stored_constant_count))
        read(stream,*,iostat=ios) stored_constants
        if(ios/=0) error stop 'Missing restart constant transport coefficients'
        if(.not.valid_constant_transport_coefficients(stored_constants)) &
          error stop 'Invalid restart constant transport coefficients'
        if(any(stored_constants/=constant_coefficients)) error stop 'Restart constant transport coefficients mismatch'
        read(stream,'(a)',iostat=ios) magic
        if(ios/=0) error stop 'Missing constant transport application header'
      else if(allocated(constant_coefficients)) then
        error stop 'Restart constant transport enable mismatch'
      endif
#if (defined(PELEF_MPI_EB_LES) || defined(PELEF_SERIAL_EB_LES)) && !defined(PELEF_MPI_EB_SPRAY)
      if(magic=='PELEF_EB_TREE_LES_V1') then
        if(.not.les_enabled) error stop 'Restart LES enable mismatch'
        read(stream,*,iostat=ios) stored_les_model,stored_les_filter,stored_les_ratios,stored_les_harmonic
        if(ios/=0) error stop 'Missing restart LES model controls'
        if(stored_les_model/=les_model.or.stored_les_filter/=les_test_model.or. &
          any(stored_les_ratios/=[les_test_ratio,les_grid_ratio]).or. &
          (stored_les_harmonic.neqv.les_harmonic_mean)) error stop 'Restart LES model mismatch'
        read(stream,*,iostat=ios) stored_les_controls
        if(ios/=0) error stop 'Missing restart LES coefficients'
        if(.not.all(ieee_is_finite(stored_les_controls))) error stop 'Invalid restart LES coefficients'
        if(any(stored_les_controls/=[les_cs,les_prt,les_sct,les_ci,les_cfl])) &
          error stop 'Restart LES coefficients mismatch'
        read(stream,'(a)',iostat=ios) magic
        if(ios/=0) error stop 'Missing LES application header'
      else if(les_enabled) then
        error stop 'Restart LES enable mismatch'
      end if
#endif
      if(magic=='PELEF_EB_TREE_ROTATING_FRAME_V1') then
        if(.not.allocated(frame)) error stop 'Restart rotating frame mismatch'
        read(stream,*,iostat=ios) stored_rotation_axis,stored_rotation
        if(ios/=0) error stop 'Missing restart rotating frame'
        if(.not.all(ieee_is_finite(stored_rotation))) error stop 'Invalid restart rotating frame'
        if(stored_rotation_axis/=rotation_axis.or. &
          any(stored_rotation/=[rotation_angular_speed,rotation_axis_point])) &
          error stop 'Restart rotating frame mismatch'
        read(stream,'(a)',iostat=ios) magic
        if(ios/=0) error stop 'Missing rotating frame application header'
      else if(allocated(frame)) then
        error stop 'Restart rotating frame mismatch'
      end if
      if(magic=='PELEF_EB_TREE_MOL_LF_V1') then
        if(.not.mol_lax_friedrichs) error stop 'Restart MOL flux mismatch'
        read(stream,'(a)',iostat=ios) magic
        if(ios/=0) error stop 'Missing MOL LF application header'
      else if(mol_lax_friedrichs) then
        error stop 'Restart MOL flux mismatch'
      end if
      if(magic=='PELEF_EB_TREE_SPATIAL_DENSITY_RATE_V2') then
        if(time_integrator/='sdc') &
          error stop 'Restart EB density-source convention mismatch'
        read(stream,'(a)',iostat=ios) magic
        if(ios/=0) error stop 'Missing EB density-source application header'
      else if(time_integrator=='sdc') then
        error stop 'Restart EB density-source convention mismatch'
      endif
      if(magic=='PELEF_EB_TREE_MOL_REGISTER_REFLUX_V2') then
        if(.not.mol_mode) error stop 'Restart MOL mode mismatch'
        read(stream,'(a)',iostat=ios) magic
        if(ios/=0) error stop 'Missing MOL application header'
      else if(mol_mode) then
        error stop 'Restart MOL mode mismatch'
      end if
      if(magic=='PELEF_EB_TREE_EMBEDDED_NOSLIP_V1') then
        if(.not.embedded_noslip) error stop 'Restart embedded no-slip mismatch'
        read(stream,'(a)',iostat=ios) magic
        if(ios/=0) error stop 'Missing embedded no-slip application header'
      else if(embedded_noslip) then
        error stop 'Restart embedded no-slip mismatch'
      end if
      if(magic=='PELEF_EB_TREE_EMBEDDED_HEAT_V1') then
        if(.not.allocated(active_embedded_temperature)) error stop 'Restart embedded wall temperature mismatch'
        read(stream,*,iostat=ios) stored_embedded_temperature
        if(ios/=0) error stop 'Missing restart embedded wall temperature'
        if(.not.ieee_is_finite(stored_embedded_temperature)) error stop 'Invalid restart embedded wall temperature'
        if(stored_embedded_temperature/=embedded_wall_temperature) error stop 'Restart embedded wall temperature mismatch'
        read(stream,'(a)',iostat=ios) magic
        if(ios/=0) error stop 'Missing embedded heat application header'
      else if(allocated(active_embedded_temperature)) then
        error stop 'Restart embedded wall temperature mismatch'
      end if
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
      if(magic=='PELEF_EB_TREE_AUTOMATIC_TIMESTEP_V1') then
        if(.not.automatic_timestep) error stop 'Restart automatic timestep enable mismatch'
        read(stream,*,iostat=ios) stored_timestep_controls,timestep_after_regrid
        if(ios/=0) error stop 'Missing automatic timestep history'
        if(.not.all(ieee_is_finite(stored_timestep_controls))) error stop 'Invalid automatic timestep history'
        if(any(stored_timestep_controls(:2)/=[initial_timestep_factor,maximum_timestep_growth])) &
          error stop 'Restart automatic timestep controls mismatch'
        previous_timestep=stored_timestep_controls(3)
        if(previous_timestep<0.or.previous_timestep/1.001_dp>time_step) error stop 'Invalid previous timestep'
        read(stream,'(a)',iostat=ios) magic
        if(ios/=0) error stop 'Missing automatic timestep application header'
      else if(automatic_timestep) then
        error stop 'Restart automatic timestep history missing'
      endif
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
      if(connected_directions) then
        write(stream,'(a)',iostat=ios) 'PELEF_EB_TREE_CONNECTED_DIRECTIONS_V1'
        if(ios/=0) error stop 'Cannot write connected-direction policy restart context'
      endif
      if(paired_pressure) then
        write(stream,'(a)',iostat=ios) 'PELEF_EB_TREE_PAIRED_PRESSURE_V1'
        if(ios/=0) error stop 'Cannot write pressure policy restart context'
      endif
      if(consistent_velocity) then
        write(stream,'(a)',iostat=ios) 'PELEF_EB_TREE_CONSISTENT_VELOCITY_V1'
        if(ios/=0) error stop 'Cannot write velocity policy restart context'
      endif
      if(separate_momentum) then
        write(stream,'(a)',iostat=ios) 'PELEF_EB_TREE_SEPARATE_MOMENTUM_V1'
        if(ios/=0) error stop 'Cannot write momentum policy restart context'
      endif
      if(allocated(constant_coefficients)) then
        write(stream,'(a)',iostat=ios) 'PELEF_EB_TREE_CONSTANT_TRANSPORT_V1'
        if(ios==0) write(stream,'(i0)',iostat=ios) size(constant_coefficients)
        if(ios==0) write(stream,'(*(es27.18e3,1x))',iostat=ios) constant_coefficients
        if(ios/=0) error stop 'Cannot write constant transport restart context'
      endif
#if (defined(PELEF_MPI_EB_LES) || defined(PELEF_SERIAL_EB_LES)) && !defined(PELEF_MPI_EB_SPRAY)
      if(les_enabled) then
        write(stream,'(a)',iostat=ios) 'PELEF_EB_TREE_LES_V1'
        if(ios==0) write(stream,*,iostat=ios) trim(les_model),' ',trim(les_test_model), &
          les_test_ratio,les_grid_ratio,les_harmonic_mean
        if(ios==0) write(stream,'(5(es27.18e3,1x))',iostat=ios) les_cs,les_prt,les_sct,les_ci,les_cfl
        if(ios/=0) error stop 'Cannot write LES restart context'
      end if
#endif
      if(allocated(frame)) then
        write(stream,'(a)',iostat=ios) 'PELEF_EB_TREE_ROTATING_FRAME_V1'
        if(ios==0) write(stream,'(i0,1x,4(es27.18e3,1x))',iostat=ios) &
          rotation_axis,rotation_angular_speed,rotation_axis_point
        if(ios/=0) error stop 'Cannot write rotating frame context'
      end if
      if(mol_lax_friedrichs) then
        write(stream,'(a)',iostat=ios) 'PELEF_EB_TREE_MOL_LF_V1'
        if(ios/=0) error stop 'Cannot write MOL LF context'
      end if
      if(time_integrator=='sdc') then
        write(stream,'(a)',iostat=ios) 'PELEF_EB_TREE_SPATIAL_DENSITY_RATE_V2'
        if(ios/=0) error stop 'Cannot write EB density-source context'
      endif
      if(mol_mode) then
        write(stream,'(a)',iostat=ios) 'PELEF_EB_TREE_MOL_REGISTER_REFLUX_V2'
        if(ios/=0) error stop 'Cannot write MOL context'
      end if
      if(embedded_noslip) then
        write(stream,'(a)',iostat=ios) 'PELEF_EB_TREE_EMBEDDED_NOSLIP_V1'
        if(ios/=0) error stop 'Cannot write embedded no-slip context'
      end if
      if(allocated(active_embedded_temperature)) then
        write(stream,'(a)',iostat=ios) 'PELEF_EB_TREE_EMBEDDED_HEAT_V1'
        if(ios==0) write(stream,'(es27.18e3)',iostat=ios) embedded_wall_temperature
        if(ios/=0) error stop 'Cannot write embedded wall temperature context'
      end if
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
      if(automatic_timestep) then
        write(stream,'(a)',iostat=ios) 'PELEF_EB_TREE_AUTOMATIC_TIMESTEP_V1'
        if(ios==0) write(stream,'(3(es27.18e3,1x),l1)',iostat=ios) &
          initial_timestep_factor,maximum_timestep_growth,previous_timestep,timestep_after_regrid
        if(ios/=0) error stop 'Cannot write automatic timestep history'
      endif
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
      bundle_identity,checkpoint_integrator,base_mole_fractions,element_names,element_atoms,frame)
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
      bundle_identity,checkpoint_integrator,base_mole_fractions,element_names,element_atoms,frame)
    if (.not.ok) error stop 'Invalid distributed tree checkpoint record'
    if (regrid_stage/=count(regrid_times(:regrid_count)<=record%time)) error stop 'Restart stage/time mismatch'
#else
    call envelope(unit,.true.)
    read(unit,'(a)',iostat=ios) extra
    if (ios/=iostat_end) error stop 'Trailing restart data'
#endif
    close(unit,iostat=ios)
    if (ios/=0) error stop 'Restart close failed'
    if(time_integrator=='sdc') then
      force_total=record%forcing
      if(allocated(sdc_ir)) sdc_ir=record%reaction_source
      if(allocated(record%heat_release%fields)) then
        sdc_heat=record%heat_release
      else
        call initialize_heat_snapshot()
      end if
    end if
    if(allocated(held_energy)) then
      held_energy=record%internal_energy;numerical_energy=record%energy_correction
      if(.not.mol_mode.and..not.record%has_internal_energy_pressure_work) &
        error stop 'Godunov held restart lacks pressure-work history'
      internal_energy_pressure_work=record%internal_energy_pressure_work
      numerical_total(5)=numerical_energy
    endif
    if(automatic_timestep.and.record%steps>0.and.previous_timestep<=0) &
      error stop 'Restart accepted timestep missing'
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
      bundle_identity,checkpoint_integrator,base_mole_fractions,element_names,element_atoms,frame)
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
  subroutine output_snapshot()
    character(len=32) :: step_label
    character(len=:), allocatable :: path,time_path
    integer :: snapshot_unit,status
    ! Publish only accepted, synchronized states, after regridding/rebalancing.
    ! The restored checkpoint is the initial state of a continuation segment.
    write(step_label,'(i0)') record%steps
    path=trim(output_file)//'.step'//trim(step_label)//'.csv'
    time_path=path//'.time.csv'
    call check_snapshot_path(path)
    call check_snapshot_path(path//'.budget.csv')
    call check_snapshot_path(time_path)
    call output(path)
#ifdef PELEF_MPI_EB_TREE
    if(mpi_rank==0) then
#endif
    open(newunit=snapshot_unit,file=time_path,status='replace',action='write',iostat=status)
    if(status/=0) error stop 'Cannot open snapshot time CSV'
    write(snapshot_unit,'(a)',iostat=status) 'step,time_s,nodes,regrid_stage'
    if(status/=0) error stop 'Cannot write snapshot time header'
    write(snapshot_unit,'(i0,",",es26.17e3,2(",",i0))',iostat=status) &
      record%steps,record%time,size(topology%nodes),regrid_stage
    if(status/=0) error stop 'Cannot write snapshot time'
    close(snapshot_unit,iostat=status)
    if(status/=0) error stop 'Cannot close snapshot time CSV'
#ifdef PELEF_MPI_EB_TREE
    endif
#endif
    last_output_step=record%steps
  end subroutine

  subroutine check_snapshot_path(path)
    character(len=*), intent(in) :: path
    if(paths_alias(path,input_path).or.paths_alias(path,restart_file).or.paths_alias(path,checkpoint_file).or. &
      paths_alias(path,level_set_file).or.paths_alias(path,section_file).or.paths_alias(path,initial_state_file)) &
      error stop 'Snapshot aliases input or checkpoint'
  end subroutine

  subroutine output(path)
    character(len=*), intent(in) :: path
    integer :: i,j,k,c,child,lo(3),hi(3),node,n(3)
    logical, allocatable :: exposed(:,:,:)
#ifdef PELEF_MPI_EB_TREE
    if(allocated(record%heat_release%fields)) then
      call write_mpi_eb_tree_csv_3d(path,topology,mpi_owners,species,mpi_state,mpi_temperature, &
        mpi_geometry,MPI_COMM_WORLD,ok,record%heat_release,held_energy)
    else
    call write_mpi_eb_tree_csv_3d(path,topology,mpi_owners,species,mpi_state,mpi_temperature, &
      mpi_geometry,MPI_COMM_WORLD,ok,internal_energy=held_energy)
    end if
    if (.not.ok) error stop 'MPI CSV write failed'
    if (mpi_rank/=0) return
#else
    open(newunit=unit,file=path,status='replace',action='write',iostat=ios)
    if (ios/=0) error stop 'Cannot open CSV'
    write(unit,'(a)',advance='no') 'node,level,i,j,k,fluid_volume,temperature,rho,rhou,rhov,rhow,rhoE'
    do c=1,size(species)
      write(unit,'(a)',advance='no') ',rhoY_'//trim(species(c)%name)
    end do
    if(allocated(record%heat_release%fields)) write(unit,'(a)',advance='no') ',heat_release_W_m3'
    if(allocated(held_energy)) write(unit,'(a)',advance='no') ',rho_internal_energy_J_m3'
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
        write(unit,'(i0,4(",",i0),*(:,",",es26.17e3))',advance='no') node,topology%nodes(node)%level,i,j,k, &
          g%volume_fraction(i,j,k)*g%dx*g%dy*g%dz,tree%nodes(node)%temperature(i,j,k),tree%nodes(node)%state(:,i,j,k)
        if(allocated(record%heat_release%fields)) write(unit,'(",",es26.17e3)',advance='no') &
          record%heat_release%fields(node)%values(1,i,j,k)
        if(allocated(held_energy)) write(unit,'(",",es26.17e3)',advance='no') held_energy%fields(node)%values(1,i,j,k)
        write(unit,*)
      end do
      end do
      end do
      end associate
      deallocate(exposed)
    end do
    close(unit,iostat=ios)
    if (ios/=0) error stop 'CSV close failed'
#endif
    open(newunit=unit,file=path//'.budget.csv',status='replace',action='write',iostat=ios)
    if (ios/=0) error stop 'Cannot open budget CSV'
#ifdef PELEF_MPI_EB_SPRAY
    if(allocated(held_energy)) then
      write(unit,'(a)') 'component,initial,final,boundary,wall,chemistry,injected,force,coordinate_energy,'// &
        'numerical_energy,internal_energy_pressure_work,residual'
      do c=1,nv
        write(unit,'(i0,11(",",es26.17e3))') c,record%initial(c),current(c),record%boundary(c), &
          record%wall(c),record%chemical(c),injected(c),force_total(c),coordinate_total(c),numerical_total(c), &
          merge(internal_energy_pressure_work,0._dp,c==5), &
          current(c)-record%initial(c)-record%boundary(c)-record%wall(c)-record%chemical(c)- &
          injected(c)-force_total(c)-coordinate_total(c)-numerical_total(c)
      enddo
    else if(allocated(frame)) then
      write(unit,'(a)') 'component,initial,final,boundary,wall,chemistry,injected,force,coordinate_energy,residual'
      do c=1,nv
        write(unit,'(i0,9(",",es26.17e3))') c,record%initial(c),current(c),record%boundary(c), &
          record%wall(c),record%chemical(c),injected(c),force_total(c),coordinate_total(c), &
          current(c)-record%initial(c)-record%boundary(c)-record%wall(c)-record%chemical(c)- &
          injected(c)-force_total(c)-coordinate_total(c)
      enddo
    else
      write(unit,'(a)') 'component,initial,final,boundary,wall,chemistry,injected,force,residual'
      do c=1,nv
        write(unit,'(i0,8(",",es26.17e3))') c,record%initial(c),current(c),record%boundary(c), &
          record%wall(c),record%chemical(c),injected(c),force_total(c), &
          current(c)-record%initial(c)-record%boundary(c)-record%wall(c)-record%chemical(c)-injected(c)-force_total(c)
      enddo
    endif
    write(*,'(a,i0)') 'Processed injection events: ',next_injection-1
#else
    if(allocated(held_energy)) then
      write(unit,'(a)') 'component,initial,final,boundary,wall,chemistry,forcing,numerical_energy,'// &
        'internal_energy_pressure_work,residual'
      do c=1,nv
        write(unit,'(i0,9(",",es26.17e3))') c,record%initial(c),current(c),record%boundary(c), &
          record%wall(c),record%chemical(c),force_total(c),numerical_total(c), &
          merge(internal_energy_pressure_work,0._dp,c==5), &
          current(c)-record%initial(c)-record%boundary(c)-record%wall(c)-record%chemical(c)- &
          force_total(c)-numerical_total(c)
      enddo
    else if(time_integrator=='sdc') then
      write(unit,'(a)') 'component,initial,final,boundary,wall,chemistry,forcing,residual'
      do c=1,nv
        write(unit,'(i0,7(",",es26.17e3))') c,record%initial(c),current(c),record%boundary(c), &
          record%wall(c),record%chemical(c),force_total(c), &
          current(c)-record%initial(c)-record%boundary(c)-record%wall(c)-record%chemical(c)-force_total(c)
      end do
    else
    write(unit,'(a)') 'component,initial,final,boundary,wall,chemistry,residual'
    do c=1,nv
      write(unit,'(i0,6(",",es26.17e3))') c,record%initial(c),current(c),record%boundary(c), &
        record%wall(c),record%chemical(c),current(c)-record%initial(c)-record%boundary(c)-record%wall(c)-record%chemical(c)
    end do
    end if
#endif
    close(unit,iostat=ios)
    if (ios/=0) error stop 'Budget close failed'
  end subroutine
end program pelef_reactive_eb_tree_3d
