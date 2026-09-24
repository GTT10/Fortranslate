module spray_mpi_amr_application_mod
  use les_filter_3d_mod, only: les_filter_weights
  use les_mpi_amr_tree_mod, only: compute_mpi_amr_tree_les_timestep_3d
  use spray_injection_profile_mod, only: integrate_spray_injection_profile,read_spray_injection_profile
  use spray_mixture_mod, only: liquid_mixture_properties
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: iso_c_binding, only: c_char,c_size_t,c_ptr,c_associated,c_null_char
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved
  use mixture_thermo_mod, only: mixture_density,mass_fractions_from_mole_fractions
  use selected_composition_mod, only: resolve_selected_composition
  use reactive_transport_3d_mod, only: reactive_transport_timestep_3d
  use simulation_config_selected_reactor_mod, only: selected_reactor_paths_alias_lexically
  use simulation_config_reactive_3d_mod, only: reactive_3d_config,resolved_reactive_boundaries_3d
  use reactive_entropy_wave_3d_problem_mod, only: initialize_reactive_problem_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    amr_tree_integer_layout_3d
  use mpi_amr_patch_tree_thermo_3d_mod, only: compute_mpi_amr_tree_cfl_3d,recover_mpi_amr_tree_temperatures_3d
  use mpi_amr_patch_tree_transport_3d_mod, only: compute_mpi_amr_tree_transport_timestep_3d
  use mpi_amr_patch_tree_restriction_3d_mod, only: average_down_mpi_amr_tree_3d
  use mpi_amr_patch_tree_tagging_3d_mod, only: tag_mpi_amr_tree_temperature_3d
  use mpi_amr_patch_tree_csv_3d_mod, only: write_mpi_amr_tree_csv_3d
  use amr_patch_tree_checkpoint_3d_mod, only: build_amr_tree_checkpoint_context_3d
  use spray_mpi_amr_checkpoint_mod, only: mpi_amr_spray_checkpoint
  use spray_mpi_amr_balance_mod, only: rebalance_mpi_amr_spray
  use mpi_amr_tree_work_mod, only: begin_tree_work,end_tree_work
  use spray_parcel_mod
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options,valid_wall_options
  use spray_mpi_mod, only: inject_mpi_cone_parcels
  use spray_mpi_amr_tree_mod, only: migrate_amr_tree_spray_parcels
  use spray_mpi_amr_reactive_mod, only: advance_mpi_amr_spray_reactive
  use spray_mpi_amr_remap_mod, only: remap_mpi_amr_spray
  implicit none
  private
  public :: run_mpi_amr_spray_application
  interface
    function system_getcwd(buffer,length) bind(C,name='getcwd') result(pointer)
      import c_char,c_size_t,c_ptr
      character(kind=c_char), intent(out) :: buffer(*)
      integer(c_size_t), value :: length
      type(c_ptr) :: pointer
    end function
  end interface
contains
  subroutine run_mpi_amr_spray_application(input_path,topology_path,spray_path,config,species,reactions, &
      transport,mole_fractions,comm,chemistry_integrator)
    character(len=*), intent(in) :: input_path,topology_path,spray_path
    character(len=*), optional, intent(in) :: chemistry_integrator
    type(reactive_3d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: mole_fractions(:)
    type(MPI_Comm), intent(in) :: comm
    type(reactive_3d_config) :: local
    type(amr_patch_tree_3d) :: tree,target,tagged
    type(amr_tree_state_3d) :: state,temperature
    type(spray_parcel), allocatable :: parcels(:)
    type(spray_liquid), allocatable :: liquids(:)
    type(spray_options) :: options
    type(spray_wall_options) :: walls
    type(spray_wall_film), allocatable :: films(:)
    real(dp), allocatable :: film_mass(:,:)
    real(dp) :: les_cs,les_prt,les_sct,les_ci,les_cfl,les_dt
    character(len=16) :: les_model
    character(len=64) :: les_test_model
    integer :: les_test_ratio,les_grid_ratio
    logical :: les_harmonic_mean
    character(len=:), allocatable :: dynamic_filter
    integer, allocatable :: dynamic_ratio,dynamic_grid
    logical, allocatable :: dynamic_harmonic
    real(dp), allocatable :: filter_weights(:)
    namelist /amr_spray_dynamic/ les_test_model,les_test_ratio,les_grid_ratio,les_harmonic_mean
    namelist /amr_spray_les/ les_cs,les_prt,les_sct,les_ci,les_cfl,les_model
    real(dp) :: acceleration(3)
    real(dp), allocatable :: force_total(:),force_change(:)
    namelist /amr_spray_force/ acceleration
    character(len=16) :: mixture_density_models(16),mixture_cp_models(16)
    character(len=16) :: mixture_vapor_pressure_models(16),mixture_boiling_models(16)
    real(dp) :: mixture_density_coefficients(4,16),mixture_density_temperature_ranges(2,16)
    real(dp) :: mixture_cp_coefficients(4,16),mixture_cp_temperature_ranges(2,16)
    real(dp) :: mixture_antoine_coefficients(4,16),mixture_antoine_temperature_ranges(2,16)
    real(dp) :: mixture_critical_temperatures(16)
    logical :: properties_present
    namelist /amr_spray_properties/ mixture_density_models,mixture_cp_models,mixture_vapor_pressure_models, &
      mixture_boiling_models,mixture_density_coefficients,mixture_density_temperature_ranges, &
      mixture_cp_coefficients,mixture_cp_temperature_ranges,mixture_antoine_coefficients, &
      mixture_antoine_temperature_ranges,mixture_critical_temperatures
    logical :: wall_active
    real(dp), allocatable :: masses(:,:),x(:),y(:),z(:),mf(:),ledger(:),global_ledger(:)
    integer, allocatable :: owners(:),assigned(:),ratios(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer(int64) :: last_id
    real(dp) :: lower(3),spacing(3),h(3),density,time,dt,transport_dt,diffusivity,correction,theta,stop_time
    integer :: rank,nranks,ierr,unit,ios,node,n(3),i,p,nv,step,attempt,particle_steps,number
    integer :: component_count,parcel_count,maximum_parcels
    integer :: adaptation_interval,refinement_levels,refinement_ratios(16),tile_shape(3),tag_buffer,maximum_nodes
    integer :: last_adaptation_step,rebalance_interval
    real(dp) :: rebalance_minimum_gain
    real(dp), allocatable :: measured_cost(:)
    logical :: sample_work,redistributed,periodic(3)
    real(dp), allocatable :: boundary(:),change(:),inlet(:,:),inlet_x(:),q(:),inlet_cube(:,:,:,:)
    real(dp) :: inlet_t(6),sound,min_h(3),inlet_dt,face_dt,inlet_temperature(2,2,2)
    integer :: face,j,k
    character(len=512) :: message
    real(dp) :: temperature_thresholds(16)
    integer, parameter :: maximum_events=1024
    integer :: profile_points(16),profile_count
    real(dp) :: profile_times(64,16),profile_rates(64,16),event_end_times(maximum_events),integrated_mass
    logical :: rate_velocity(16),velocity_present
    real(dp) :: discharge_coefficients(16),event_speeds(maximum_events),rate_denominator,liquid_rho,liquid_heat,liquid_energy
    namelist /amr_spray_velocity/ rate_velocity,discharge_coefficients
    character(len=1024) :: profile_files(16)
    namelist /amr_spray_rate_files/ profile_files
    logical :: rates_present
    namelist /amr_spray_rates/ profile_points,profile_times,profile_rates,event_end_times
    integer :: injection_events,event,event_parcels(maximum_events),event_ports(maximum_events),port_count
    real(dp) :: port_centers(3,16),port_axes(3,16),port_speeds(16),port_half_angles(16), &
      port_diameters(16),port_temperatures(16),port_fractions(16,16),port_radii(16)
    type(spray_nozzle_model) :: port_nozzles(16)
    type(spray_size_options) :: port_sizes(16)
    logical :: sizes_present
    namelist /amr_spray_sizes/ port_sizes
    real(dp) :: injection_times(maximum_events),event_mass(maximum_events)
    character(len=24) :: vapor_names(16)
    real(dp) :: fractions(16),liquid_density(16),liquid_cp(16),reference_temperature(16), &
      reference_pressure(16),latent_heat(16)
    real(dp) :: total_mass,center(3),axis(3),speed,half_angle,diameter,drop_temperature,regrid_time,local_mass
    logical :: ok,subcycling,regridded
    character(len=32) :: faces(6)
    character(len=1024) :: regrid_topology,particle_path
    character(len=24) :: suffix
    character(len=:), allocatable :: header,context,spray_context,expected_tree,method
    character(len=1024) :: checkpoint_file,restart_file,context_line
    real(dp) :: checkpoint_time
    logical :: stop_after_checkpoint,checkpoint_written
    namelist /amr_spray_checkpoint/ checkpoint_file,restart_file,checkpoint_time,stop_after_checkpoint
    namelist /amr_spray/ component_count,vapor_names,fractions,liquid_density,liquid_cp,reference_temperature, &
      reference_pressure,latent_heat,parcel_count,total_mass,center,axis,speed,half_angle,diameter,drop_temperature, &
      maximum_parcels,options,subcycling,regrid_time,regrid_topology, &
      injection_events,injection_times,event_parcels,event_mass,event_ports,port_count,port_centers,port_axes, &
      port_speeds,port_half_angles,port_diameters,port_temperatures,port_fractions,port_radii,port_nozzles, &
      adaptation_interval,refinement_levels,refinement_ratios,tile_shape,tag_buffer,maximum_nodes,temperature_thresholds, &
      rebalance_interval,rebalance_minimum_gain,walls
    call MPI_Comm_rank(comm,rank,ierr)
    call MPI_Comm_size(comm,nranks,ierr)
    faces=resolved_reactive_boundaries_3d(config)
    periodic=faces([1,3,5])=='periodic'
    ! Inputs are shared files, as for the existing MPI tree applications.
    component_count=1;vapor_names='';vapor_names(1)='H2O';fractions=0;fractions(1)=1
    liquid_density=1000;liquid_cp=4184;reference_temperature=373.15_dp
    reference_pressure=101325;latent_heat=2.257e6_dp
    parcel_count=1;total_mass=1e-10_dp;center=0;axis=[1._dp,0._dp,0._dp]
    speed=10;half_angle=0;diameter=1e-5_dp;drop_temperature=350;maximum_parcels=100000
    subcycling=.false.;regrid_time=-1;regrid_topology='';regridded=.false.
    rebalance_interval=0;rebalance_minimum_gain=.1_dp
    adaptation_interval=0;refinement_levels=1;refinement_ratios=2;tile_shape=4;tag_buffer=1;maximum_nodes=4096
    temperature_thresholds=1200;last_adaptation_step=0
    injection_events=0;injection_times=0;event_parcels=0;event_mass=0;event_ports=1
    port_count=0;port_centers=0;port_axes=0;port_axes(1,:)=1;port_speeds=10;port_half_angles=0
    port_diameters=1e-5_dp;port_temperatures=350;port_fractions=-1;port_radii=0
    mixture_density_models='constant';mixture_cp_models='constant'
    mixture_vapor_pressure_models='clausius';mixture_boiling_models='reject'
    mixture_density_coefficients=0;mixture_density_temperature_ranges=0
    mixture_cp_coefficients=0;mixture_cp_temperature_ranges=0
    mixture_antoine_coefficients=0;mixture_antoine_temperature_ranges=0;mixture_critical_temperatures=0
    rate_velocity=.false.;discharge_coefficients=1
    profile_files=''
    profile_points=0;profile_times=0;profile_rates=0;event_end_times=-1
    port_sizes=spray_size_options()
    acceleration=0
    les_test_model='box';les_test_ratio=2;les_grid_ratio=1;les_harmonic_mean=.false.
    les_cs=0;les_prt=.7_dp;les_sct=.7_dp;les_ci=0;les_cfl=.1_dp;les_model='smagorinsky'
    open(newunit=unit,file=spray_path,status='old',action='read',iostat=ios)
    if (ios/=0) call fatal('Cannot open AMR spray input')
    read(unit,nml=amr_spray,iostat=ios)
    if (ios/=0) call fatal('Invalid amr_spray namelist')
    checkpoint_file='';restart_file='';checkpoint_time=-1;stop_after_checkpoint=.false.;checkpoint_written=.false.
    rewind(unit)
    read(unit,nml=amr_spray_checkpoint,iostat=ios)
    if (ios>0) call fatal('Invalid amr_spray_checkpoint namelist')
    rewind(unit)
    read(unit,nml=amr_spray_les,iostat=ios)
    if (ios>0) call fatal('Invalid amr_spray_les namelist')
    rewind(unit)
    read(unit,nml=amr_spray_dynamic,iostat=ios)
    if(ios>0) call fatal('Invalid amr_spray_dynamic namelist')
    rewind(unit)
    read(unit,nml=amr_spray_force,iostat=ios)
    if (ios>0) call fatal('Invalid amr_spray_force namelist')
    rewind(unit)
    read(unit,nml=amr_spray_properties,iostat=ios)
    if (ios>0) call fatal('Invalid amr_spray_properties namelist')
    properties_present=ios==0
    rewind(unit)
    read(unit,nml=amr_spray_sizes,iostat=ios)
    if (ios>0) call fatal('Invalid amr_spray_sizes namelist')
    sizes_present=ios==0
    rewind(unit)
    read(unit,nml=amr_spray_rates,iostat=ios)
    if (ios>0) call fatal('Invalid amr_spray_rates namelist')
    rates_present=ios==0
    rewind(unit)
    read(unit,nml=amr_spray_velocity,iostat=ios)
    if (ios>0) call fatal('Invalid amr_spray_velocity namelist')
    velocity_present=ios==0
    rewind(unit)
    read(unit,nml=amr_spray_rate_files,iostat=ios)
    close(unit)
    if (ios>0) call fatal('Invalid amr_spray_rate_files namelist')
    if (.not.all(ieee_is_finite(acceleration))) call fatal('Nonfinite acceleration')
    if (any(abs(acceleration)>1e6_dp)) call fatal('Invalid acceleration')
    if (.not.all(ieee_is_finite([les_cs,les_prt,les_sct,les_ci,les_cfl]))) call fatal('Nonfinite LES settings')
    if (les_cs<0.or.les_cs>2.or.les_ci<0.or.les_ci>1) call fatal('Invalid LES coefficients')
    if (min(les_prt,les_sct)<1e-6_dp.or.max(les_prt,les_sct)>1e6_dp) call fatal('Invalid LES PrT/ScT')
    if (les_cfl<=0.or.les_cfl>.2_dp) call fatal('Invalid LES diffusion CFL')
    if (les_model/='smagorinsky'.and.les_model/='wale'.and.les_model/='vreman'.and.les_model/='dynamic') &
      call fatal('Invalid LES model')
    if(les_model=='dynamic') then
      call les_filter_weights(les_test_model,les_test_ratio,filter_weights,ok)
      if(.not.ok.or.les_grid_ratio<1.or.les_grid_ratio>64) call fatal('Invalid dynamic LES filter')
      if(any(faces/='periodic'.and.faces/='inflow'.and.faces/='outflow')) &
        call fatal('Unsupported dynamic LES boundary')
      if(any(config%outlet_pressure/=0)) call fatal('Dynamic LES pressure outlet is not supported')
      dynamic_filter=trim(les_test_model)
      allocate(dynamic_ratio,dynamic_grid,dynamic_harmonic)
      dynamic_ratio=les_test_ratio;dynamic_grid=les_grid_ratio;dynamic_harmonic=les_harmonic_mean
    else
      if(les_test_model/='box'.or.les_test_ratio/=2.or.les_grid_ratio/=1.or.les_harmonic_mean) &
        call fatal('Dynamic filter controls require dynamic LES')
    end if
    if (.not.ieee_is_finite(checkpoint_time)) call fatal('Invalid checkpoint time')
    if ((checkpoint_time>=0).neqv.(len_trim(checkpoint_file)>0)) call fatal('Checkpoint requires path and time')
    if (stop_after_checkpoint.and.checkpoint_time<0) call fatal('Stop requires a checkpoint')
    if (len_trim(checkpoint_file)>0) then
      call check_checkpoint_path(trim(checkpoint_file))
      call check_checkpoint_path(trim(checkpoint_file)//'.gas')
      call check_checkpoint_path(trim(checkpoint_file)//'.particles')
      call check_checkpoint_path(trim(checkpoint_file)//'.pending')
    end if
    if (len_trim(restart_file)>0) then
      call check_checkpoint_path(trim(restart_file))
      call check_checkpoint_path(trim(restart_file)//'.gas')
      call check_checkpoint_path(trim(restart_file)//'.particles')
      call check_checkpoint_path(trim(restart_file)//'.pending')
    end if
    if (component_count<1.or.component_count>16) call fatal('Invalid liquid component count')
    if (.not.ieee_is_finite(regrid_time)) call fatal('Invalid regrid time')
    if ((regrid_time>=0).neqv.(len_trim(regrid_topology)>0)) call fatal('Regrid needs time and topology')
    if (rebalance_interval<0.or..not.ieee_is_finite(rebalance_minimum_gain)) call fatal('Invalid rebalance controls')
    if (rebalance_minimum_gain<0.or.rebalance_minimum_gain>=1) call fatal('Invalid rebalance minimum gain')
    if (adaptation_interval<0) call fatal('Invalid adaptation interval')
    if (adaptation_interval>0) then
      if (regrid_time>=0) call fatal('Automatic adaptation cannot use a prescribed regrid')
      if (refinement_levels<0.or.refinement_levels>16) call fatal('Invalid refinement depth')
      if (any(refinement_ratios(:refinement_levels)<2)) call fatal('Invalid refinement ratios')
      if (.not.all(ieee_is_finite(temperature_thresholds(:refinement_levels)))) call fatal('Nonfinite thresholds')
      if (any(temperature_thresholds(:refinement_levels)<=0)) call fatal('Invalid thresholds')
      do i=2,refinement_levels
        if (temperature_thresholds(i)<temperature_thresholds(i-1)) call fatal('Thresholds must not decrease')
      end do
      if (any(tile_shape<1).or.tag_buffer<0.or.maximum_nodes<1) call fatal('Invalid tagging bounds')
    end if
    if (port_count<0.or.port_count>16) call fatal('Invalid port count')
    if (port_count==0) then
      port_count=1;port_centers(:,1)=center;port_axes(:,1)=axis;port_speeds(1)=speed
      port_half_angles(1)=half_angle;port_diameters(1)=diameter;port_temperatures(1)=drop_temperature
    end if
    do i=1,port_count
      if (all(port_fractions(:,i)==-1)) port_fractions(:,i)=fractions
      if (.not.all(ieee_is_finite(port_fractions(:component_count,i)))) call fatal('Nonfinite port composition')
      if (any(port_fractions(:component_count,i)<0).or. &
          abs(sum(port_fractions(:component_count,i))-1._dp)>1e-12_dp) call fatal('Invalid port composition')
    end do
    if (injection_events<0.or.injection_events>maximum_events) call fatal('Invalid injection event count')
    if (injection_events==0) then
      injection_events=1;injection_times(1)=0;event_parcels(1)=parcel_count;event_mass(1)=total_mass
    end if
    if (any(len_trim(profile_files(port_count+1:))>0)) call fatal('Unused profile files')
    do i=1,port_count
      if (len_trim(profile_files(i))==0) cycle
      if (profile_points(i)/=0.or.any(profile_times(:,i)/=0).or.any(profile_rates(:,i)/=0)) &
        call fatal('Profile file conflicts with inline table')
      call read_spray_injection_profile(trim(profile_files(i)),profile_times(:,i),profile_rates(:,i),profile_points(i),ok)
      if (.not.ok) call fatal('Invalid injection rate file')
      rates_present=.true.
    end do
    if (any(profile_points(port_count+1:)/=0).or.any(profile_times(:,port_count+1:)/=0).or. &
        any(profile_rates(:,port_count+1:)/=0)) call fatal('Unused port rate profile entries')
    if (any(event_end_times(injection_events+1:)/=-1)) call fatal('Unused event integration intervals')
    do i=1,port_count
      profile_count=profile_points(i)
      if (profile_count/=0.and.(profile_count<2.or.profile_count>64)) call fatal('Invalid rate profile count')
      if (any(profile_times(profile_count+1:,i)/=0).or.any(profile_rates(profile_count+1:,i)/=0)) &
        call fatal('Unused rate profile entries')
      if (profile_count==0) cycle
      call integrate_spray_injection_profile(profile_times(:profile_count,i),profile_rates(:profile_count,i), &
        profile_times(1,i),profile_times(profile_count,i),integrated_mass,ok)
      if (.not.ok) call fatal('Invalid port rate profile')
    end do
    do i=1,injection_events
      if (event_ports(i)<1.or.event_ports(i)>port_count) call fatal('Invalid event port')
      profile_count=profile_points(event_ports(i))
      if (profile_count==0) then
        if (event_end_times(i)/=-1) call fatal('Integration interval requires rate profile')
        cycle
      end if
      if (event_mass(i)/=0) call fatal('Rate profile requires zero explicit event mass')
      call integrate_spray_injection_profile(profile_times(:profile_count,event_ports(i)), &
        profile_rates(:profile_count,event_ports(i)),injection_times(i),event_end_times(i),integrated_mass,ok)
      if (.not.ok) call fatal('Invalid injection integration interval')
      event_mass(i)=integrated_mass
    end do
    if (.not.all(ieee_is_finite(injection_times(:injection_events))).or. &
        .not.all(ieee_is_finite(event_mass(:injection_events)))) call fatal('Nonfinite injection schedule')
    if (any(injection_times(:injection_events)<0).or.any(event_mass(:injection_events)<0).or. &
        any(event_parcels(:injection_events)<0).or.any(event_parcels(:injection_events)>maximum_parcels)) &
      call fatal('Invalid injection schedule values')
    if (any(event_ports(:injection_events)<1).or.any(event_ports(:injection_events)>port_count)) &
      call fatal('Invalid event port')
    do i=1,injection_events
      if ((event_parcels(i)==0).neqv.(event_mass(i)==0)) call fatal('Injection count/mass mismatch')
      if (i>1) then
        if (injection_times(i)<injection_times(i-1)) call fatal('Injection times must not decrease')
      end if
    end do
    call check_output(trim(config%output_file))
    write(suffix,'(i0)') rank
    call check_output(trim(config%output_file)//'.particles.'//trim(suffix)//'.csv')
    options%maximum_parcels=maximum_parcels
    if (.not.valid_spray_options(options)) call fatal('Invalid spray options')
    if (.not.valid_wall_options(walls,faces)) call fatal('Invalid wall impact/film options')
    wall_active=any(walls%enabled)
    if (walls%evolve.and..not.wall_active) call fatal('Film evolution requires enabled walls')
    if (wall_active) then
      if (options%surface_tension<1e-8_dp.or.options%surface_tension>1e3_dp) call fatal('Impact needs surface tension')
      do i=1,port_count
        if (liquid_viscosity_at_temperature(options,port_temperatures(i))<=0) &
          call fatal('Impact needs valid liquid viscosity at each port temperature')
      end do
      allocate(films(0),film_mass(component_count,0))
      write(suffix,'(i0)') rank
      call check_output(trim(config%output_file)//'.films.'//trim(suffix)//'.csv')
    end if
    if(any(mixture_density_models(component_count+1:)/='constant')) call fatal('Unused mixture property entries')
    if(any(mixture_cp_models(component_count+1:)/='constant')) call fatal('Unused mixture property entries')
    if(any(mixture_vapor_pressure_models(component_count+1:)/='clausius')) call fatal('Unused mixture property entries')
    if(any(mixture_boiling_models(component_count+1:)/='reject')) call fatal('Unused mixture property entries')
    if(any(mixture_critical_temperatures(component_count+1:)/=0)) call fatal('Unused mixture property entries')
    if(any(mixture_density_coefficients(:,component_count+1:)/=0)) call fatal('Unused mixture property entries')
    if(any(mixture_density_temperature_ranges(:,component_count+1:)/=0)) call fatal('Unused mixture property entries')
    if(any(mixture_cp_coefficients(:,component_count+1:)/=0)) call fatal('Unused mixture property entries')
    if(any(mixture_cp_temperature_ranges(:,component_count+1:)/=0)) call fatal('Unused mixture property entries')
    if(any(mixture_antoine_coefficients(:,component_count+1:)/=0)) call fatal('Unused mixture property entries')
    if(any(mixture_antoine_temperature_ranges(:,component_count+1:)/=0)) call fatal('Unused mixture property entries')
    allocate(liquids(component_count),parcels(0),masses(component_count,0))
    do i=1,component_count
      call initialize_spray_liquid(species,trim(vapor_names(i)),liquid_density(i),liquid_cp(i), &
        reference_temperature(i),reference_pressure(i),latent_heat(i),liquids(i),ok,mixture_vapor_pressure_models(i), &
        mixture_antoine_coefficients(:,i),mixture_antoine_temperature_ranges(:,i))
      if (.not.ok) call fatal('Invalid liquid properties/species')
      liquids(i)%density_model=mixture_density_models(i)
      liquids(i)%density_coefficients=mixture_density_coefficients(:,i)
      liquids(i)%density_temperature_range=mixture_density_temperature_ranges(:,i)
      liquids(i)%heat_capacity_model=mixture_cp_models(i)
      liquids(i)%heat_capacity_coefficients=mixture_cp_coefficients(:,i)
      liquids(i)%heat_capacity_temperature_range=mixture_cp_temperature_ranges(:,i)
      liquids(i)%boiling_model=mixture_boiling_models(i)
      liquids(i)%critical_temperature=mixture_critical_temperatures(i)
      if(.not.valid_spray_liquid(liquids(i))) call fatal('Invalid mixture component property model')
    end do
    if (any(rate_velocity(port_count+1:)).or.any(discharge_coefficients(port_count+1:)/=1)) &
      call fatal('Unused velocity controls')
    do i=1,port_count
      if (.not.ieee_is_finite(discharge_coefficients(i))) call fatal('Nonfinite discharge coefficient')
      if (.not.rate_velocity(i)) then
        if (discharge_coefficients(i)/=1) call fatal('Discharge coefficient requires rate velocity')
        cycle
      end if
      if (profile_points(i)<2) call fatal('Rate velocity requires a rate profile')
      if (discharge_coefficients(i)<1e-6_dp.or.discharge_coefficients(i)>1) call fatal('Invalid discharge coefficient')
      if (.not.ieee_is_finite(port_radii(i))) call fatal('Nonfinite rate outlet radius')
      if (port_radii(i)<1e-12_dp.or.port_radii(i)>1) call fatal('Invalid rate outlet radius')
    end do
    event_speeds=0
    do i=1,injection_events
      p=event_ports(i);event_speeds(i)=port_speeds(p)
      if (.not.rate_velocity(p)) cycle
      event_speeds(i)=0
      if (event_mass(i)==0) cycle
      call liquid_mixture_properties(liquids,port_fractions(:component_count,p),port_temperatures(p), &
        liquid_rho,liquid_heat,liquid_energy,ok)
      if (.not.ok) call fatal('Invalid rate-driven liquid state')
      rate_denominator=(event_end_times(i)-injection_times(i))*liquid_rho*acos(-1._dp)* &
        port_radii(p)**2*discharge_coefficients(p)
      if (rate_denominator<=0) call fatal('Invalid rate-driven interval')
      if (event_mass(i)>1e7_dp*rate_denominator) call fatal('Rate-driven speed exceeds bounds')
      event_speeds(i)=event_mass(i)/rate_denominator
    end do
    call read_topology(topology_path,tree)
    if (regrid_time>=0) call read_topology(trim(regrid_topology),target)
    open(newunit=unit,status='scratch',form='formatted',iostat=ios)
    if (ios/=0) call fatal('Cannot allocate checkpoint context')
    write(unit,nml=amr_spray,iostat=ios)
    if (ios/=0) call fatal('Cannot serialize spray context')
    if (les_cs/=0.or.les_ci/=0.or.les_prt/=.7_dp.or.les_sct/=.7_dp.or. &
        les_cfl/=.1_dp.or.les_model/='smagorinsky') then
      write(unit,nml=amr_spray_les,iostat=ios)
      if (ios/=0) call fatal('Cannot serialize LES context')
    end if
    if(les_model=='dynamic') then
      write(unit,nml=amr_spray_dynamic,iostat=ios)
      if(ios/=0) call fatal('Cannot serialize dynamic LES context')
    end if
    if (any(acceleration/=0)) then
      write(unit,nml=amr_spray_force,iostat=ios)
      if (ios/=0) call fatal('Cannot serialize force context')
    end if
    if (properties_present) then
      write(unit,nml=amr_spray_properties,iostat=ios)
      if (ios/=0) call fatal('Cannot serialize liquid properties context')
    end if
    if (sizes_present) then
      write(unit,nml=amr_spray_sizes,iostat=ios)
      if (ios/=0) call fatal('Cannot serialize size distribution context')
    end if
    if (rates_present) then
      write(unit,nml=amr_spray_rates,iostat=ios)
      if (ios/=0) call fatal('Cannot serialize injection rate context')
    end if
    if (velocity_present) then
      write(unit,nml=amr_spray_velocity,iostat=ios)
      if (ios/=0) call fatal('Cannot serialize rate velocity context')
    end if
    rewind(unit);spray_context=''
    do
      read(unit,'(a)',iostat=ios) context_line
      if (ios<0) exit
      if (ios/=0) call fatal('Cannot read spray context')
      spray_context=spray_context//trim(context_line)//achar(10)
    end do
    close(unit)
    spray_context=spray_context//tree_key(tree)
    if (regrid_time>=0) spray_context=spray_context//tree_key(target)
    method='explicit'
    if (present(chemistry_integrator)) method=chemistry_integrator
    call build_amr_tree_checkpoint_context_3d(config,species,reactions,transport,mole_fractions, &
      'AMR_SPRAY_V1:'//trim(config%thermo_model),method,spray_context,context,ok)
    if (.not.ok) call fatal('Cannot bind AMR spray checkpoint context')
    owners=[(modulo(node-1,nranks),node=1,size(tree%nodes))]
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,ok)
    if (.not.ok) call fatal('Invalid tree coordinate layout')
    lower=[config%x_lower,config%y_lower,config%z_lower]
    spacing=([config%x_upper,config%y_upper,config%z_upper]-lower)/real(tree%nodes(1)%shape,dp)
    nv=reactive_nvar(size(species))
    allocate(force_change(nv));force_change=0
    if (any(acceleration/=0)) then
      allocate(force_total(nv));force_total=0
    end if
    allocate(state%fields(size(tree%nodes)),temperature%fields(size(tree%nodes)),mf(size(species)))
    allocate(ledger(nv),global_ledger(nv),boundary(nv),change(nv),inlet(nv,6),inlet_x(size(species)), &
      q(nv),inlet_cube(nv,2,2,2))
    inlet=0;inlet_t=0
    call mass_fractions_from_mole_fractions(species,mole_fractions,mf,ok)
    if (.not.ok) call fatal('Invalid default inlet composition')
    if (any(faces=='inflow').and.config%inflow_composition_count/=0) then
      call resolve_selected_composition('Inflow',config%inflow_composition_count,config%inflow_composition_species, &
        config%inflow_composition_mole_fractions,species,inlet_x,ok,message)
      if (.not.ok) call fatal('Invalid inflow composition')
      call mass_fractions_from_mole_fractions(species,inlet_x,mf,ok)
      if (.not.ok) call fatal('Invalid inflow mass fractions')
    end if
    do face=1,6
      if (faces(face)/='inflow') cycle
      q=0;q(1)=mixture_density(species,mf,config%inflow_pressure(face),config%inflow_temperature(face),ok)
      if (.not.ok) call fatal('Invalid inflow density')
      q(2:4)=config%inflow_velocity(:,face);q(5)=config%inflow_pressure(face);q(6:)=mf
      call reactive_primitive_to_conserved(species,q,inlet(:,face),inlet_t(face),sound,ok)
      if (.not.ok) call fatal('Invalid inflow state')
    end do
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      h=spacing*real(domain(:,1),dp)/real(domain(:,tree%nodes(node)%level+1),dp)
      x=[(lower(1)+(real(origin(1,node),dp)+i-.5_dp)*h(1),i=1,n(1))]
      y=[(lower(2)+(real(origin(2,node),dp)+i-.5_dp)*h(2),i=1,n(2))]
      z=[(lower(3)+(real(origin(3,node),dp)+i-.5_dp)*h(3),i=1,n(3))]
      allocate(state%fields(node)%values(nv,n(1),n(2),n(3)),temperature%fields(node)%values(1,n(1),n(2),n(3)))
      local=config;local%nx=n(1);local%ny=n(2);local%nz=n(3)
      call initialize_reactive_problem_3d(species,local,x,y,z,state%fields(node)%values, &
        temperature%fields(node)%values(1,:,:,:),density,mf,ok,mole_fractions)
      if (.not.ok) call fatal('Gas initialization failed')
    end do
    call average_down_mpi_amr_tree_3d(tree,owners,state,nv,comm,ok)
    if (.not.ok) call fatal('Initial parent synchronization failed')
    call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,state,temperature,comm,ok)
    if (.not.ok) call fatal('Initial gas recovery failed')
    time=0;step=0;event=1;last_id=0;global_ledger=0;boundary=0
    expected_tree=tree_key(tree)
    if (len_trim(restart_file)>0) then
      call checkpoint(.false.,trim(restart_file))
      if (time>config%final_time.or.step>config%maximum_steps) call fatal('Restart exceeds run limits')
      if (event/=count(injection_times(:injection_events)<=time)+1) call fatal('Restart injection event mismatch')
      if (regridded.neqv.(regrid_time>=0.and.time>=regrid_time)) call fatal('Restart regrid event mismatch')
      if (regridded) expected_tree=tree_key(target)
      if (adaptation_interval==0) then
        if (tree_key(tree)/=expected_tree) call fatal('Restart topology mismatch')
      else if (step<adaptation_interval) then
        if (tree_key(tree)/=expected_tree) call fatal('Restart initial topology mismatch')
      else
        call amr_tree_integer_layout_3d(tree,origin,domain,ratios,ok)
        if (.not.ok) call fatal('Restart adaptive layout invalid')
        if (size(tree%nodes)>maximum_nodes.or.tree%level_count()>refinement_levels+1) &
          call fatal('Restart adaptation bounds mismatch')
        if (any(ratios(2:)/=refinement_ratios(:size(ratios)-1))) call fatal('Restart refinement ratios mismatch')
      end if
      last_adaptation_step=step
      if (checkpoint_time>=0.and.checkpoint_time<time) call fatal('Checkpoint time precedes restart')
    end if
    do
      if (.not.regridded.and.regrid_time>=0) then
        if (time==regrid_time) then
          assigned=[(modulo(node-1,nranks),node=1,size(target%nodes))]
          if (wall_active) then
            call remap_mpi_amr_spray(comm,species,tree,owners,state,temperature,parcels,masses,target,assigned, &
              lower,spacing,periodic,maximum_parcels,ok, &
              wall_films=films,film_mass=film_mass,maximum_films=walls%maximum_films)
          else
            call remap_mpi_amr_spray(comm,species,tree,owners,state,temperature,parcels,masses,target,assigned, &
              lower,spacing,periodic,maximum_parcels,ok)
          end if
          if (.not.ok) call fatal('Coupled regrid failed')
          regridded=.true.
        end if
      end if
      if (adaptation_interval>0) then
        if (step>last_adaptation_step) then
          if (mod(step,adaptation_interval)==0) then
            call tag_mpi_amr_tree_temperature_3d(tree,owners,temperature,temperature_thresholds(:refinement_levels), &
              refinement_ratios(:refinement_levels),tile_shape,tag_buffer,maximum_nodes,periodic,tagged,comm,ok)
            if (.not.ok) call fatal('Automatic temperature tagging failed')
            if (tree_key(tagged)/=tree_key(tree)) then
              assigned=[(modulo(node-1,nranks),node=1,size(tagged%nodes))]
              if (wall_active) then
                call remap_mpi_amr_spray(comm,species,tree,owners,state,temperature,parcels,masses,tagged,assigned, &
                  lower,spacing,periodic,maximum_parcels,ok, &
                  wall_films=films,film_mass=film_mass,maximum_films=walls%maximum_films)
              else
                call remap_mpi_amr_spray(comm,species,tree,owners,state,temperature,parcels,masses,tagged,assigned, &
                  lower,spacing,periodic,maximum_parcels,ok)
              end if
              if (.not.ok) call fatal('Adaptive gas/particle remap failed')
            end if
            last_adaptation_step=step
            if (rank==0) write(*,'(a,i0,a,i0)') 'Adaptation check step: ',step,' Nodes: ',size(tree%nodes)
          end if
        end if
      end if
      do while (event<=injection_events)
        if (time/=injection_times(event)) exit
        call inject_event()
        event=event+1
      end do
      if (.not.checkpoint_written.and.checkpoint_time>=0) then
        if (time==checkpoint_time) then
          call checkpoint(.true.,trim(checkpoint_file))
          checkpoint_written=.true.
          if (stop_after_checkpoint) then
            if (rank==0) write(*,'(a,es24.16)') 'Stopped after checkpoint at time: ',time
            return
          end if
        end if
      end if
      if (time>=config%final_time) exit
      if (step>=config%maximum_steps) call fatal('Maximum steps reached')
      stop_time=config%final_time
      if (.not.checkpoint_written.and.checkpoint_time>time) stop_time=min(stop_time,checkpoint_time)
      if (event<=injection_events) stop_time=min(stop_time,injection_times(event))
      if (.not.regridded.and.regrid_time>time) stop_time=min(stop_time,regrid_time)
      ! Finest-level bounds also safely bound recursive gas stepping.
      call compute_mpi_amr_tree_cfl_3d(tree,owners,species,state,temperature,spacing,config%cfl,dt,comm,ok)
      if (.not.ok) call fatal('Gas CFL calculation failed')
      call inlet_timestep()
      dt=min(dt,inlet_dt,stop_time-time)
      if (config%transport_enabled) then
        call compute_mpi_amr_tree_transport_timestep_3d(tree,owners,species,transport,state,temperature,spacing, &
          config%transport_cfl,config%viscosity_enabled,config%thermal_conduction_enabled,config%species_diffusion_enabled, &
          transport_dt,diffusivity,comm,ok,faces,config%wall_temperature)
        if (.not.ok) call fatal('Transport timestep failed')
        dt=min(dt,transport_dt)
      end if
      if (les_cs>0.or.les_ci>0.or.les_model=='dynamic') then
        call compute_mpi_amr_tree_les_timestep_3d(tree,owners,species,state,temperature,spacing,les_cfl,faces, &
          les_cs,les_prt,les_sct,les_dt,diffusivity,comm,ok,inflow_state=inlet,inflow_temperature=inlet_t, &
          wall_temperature=config%wall_temperature,wall_velocity=config%wall_velocity,model=les_model,ci=les_ci, &
          outlet_pressure=config%outlet_pressure,test_model=dynamic_filter,test_ratio=dynamic_ratio, &
          grid_ratio=dynamic_grid,harmonic_mean=dynamic_harmonic)
        if (.not.ok) call fatal('LES timestep failed')
        dt=min(dt,les_dt)
      end if
      sample_work=.false.
      if (rebalance_interval>0) sample_work=mod(step+1,rebalance_interval)==0
      do attempt=1,20
        if (time+dt<=time) call fatal('Timestep cannot advance time')
        if (sample_work) then
          call begin_tree_work(size(tree%nodes),ok)
          if (.not.ok) call fatal('Cannot start gas work sampling')
        end if
        correction=0;theta=1
        if (wall_active) then
          call advance_mpi_amr_spray_reactive(comm,tree,owners,species,reactions,transport,liquids,options, &
            state,temperature,parcels,masses,lower,spacing,dt,config%reconstruction,config%limiter,config%riemann_solver, &
            config%chemistry_enabled,config%transport_enabled.and.config%viscosity_enabled, &
            config%transport_enabled.and.config%thermal_conduction_enabled, &
            config%transport_enabled.and.config%species_diffusion_enabled, &
            config%transport_enabled.and.config%barodiffusion_enabled, &
            config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance,particle_steps,correction,theta, &
            ok,last_id,chemistry_integrator,subcycling,faces,change,inlet,inlet_t, &
            config%wall_temperature,config%wall_velocity, &
            wall_films=films,film_mass=film_mass,film_settings=walls, &
            les_settings=[les_cs,les_prt,les_sct,les_ci],les_model=les_model,outlet_pressure=config%outlet_pressure, &
            acceleration=acceleration,source_change=force_change, &
            les_test_model=dynamic_filter,les_test_ratio=dynamic_ratio, &
            les_grid_ratio=dynamic_grid,les_harmonic_mean=dynamic_harmonic)
        else
          call advance_mpi_amr_spray_reactive(comm,tree,owners,species,reactions,transport,liquids,options, &
            state,temperature,parcels,masses,lower,spacing,dt,config%reconstruction,config%limiter,config%riemann_solver, &
            config%chemistry_enabled,config%transport_enabled.and.config%viscosity_enabled, &
            config%transport_enabled.and.config%thermal_conduction_enabled, &
            config%transport_enabled.and.config%species_diffusion_enabled, &
            config%transport_enabled.and.config%barodiffusion_enabled, &
            config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance,particle_steps,correction,theta, &
            ok,last_id,chemistry_integrator,subcycling,faces,change,inlet,inlet_t, &
            config%wall_temperature,config%wall_velocity, &
            les_settings=[les_cs,les_prt,les_sct,les_ci],les_model=les_model,outlet_pressure=config%outlet_pressure, &
            acceleration=acceleration,source_change=force_change, &
            les_test_model=dynamic_filter,les_test_ratio=dynamic_ratio, &
            les_grid_ratio=dynamic_grid,les_harmonic_mean=dynamic_harmonic)
        end if
        if (sample_work) call end_tree_work(measured_cost)
        if (ok) exit
        dt=.5_dp*dt
      end do
      if (.not.ok) call fatal('Coupled step failed after retries')
      if (any(abs(boundary)>huge(1._dp)-abs(change))) call fatal('Boundary ledger overflow')
      if (allocated(force_total)) then
        if (any(abs(force_total)>huge(1._dp)-abs(force_change))) call fatal('Source ledger overflow')
        force_total=force_total+force_change
      end if
      boundary=boundary+change
      if (dt==stop_time-time) then
        time=stop_time
      else
        time=time+dt
      end if
      step=step+1
      if (sample_work) then
        if (wall_active) then
          call rebalance_mpi_amr_spray(comm,tree,owners,nv,state,temperature,parcels,masses,lower,spacing, &
            periodic,measured_cost,rebalance_minimum_gain,maximum_parcels,redistributed,ok, &
            wall_films=films,film_mass=film_mass,maximum_films=walls%maximum_films)
        else
          call rebalance_mpi_amr_spray(comm,tree,owners,nv,state,temperature,parcels,masses,lower,spacing, &
            periodic,measured_cost,rebalance_minimum_gain,maximum_parcels,redistributed,ok)
        end if
        if (.not.ok) call fatal('Coupled gas/particle rebalance failed')
        if (rank==0) write(*,'(a,i0,a,l1)') 'Rebalance step: ',step,' Migrated: ',redistributed
      end if
    end do
    call write_mpi_amr_tree_csv_3d(trim(config%output_file),tree,owners,species,state,temperature,lower,spacing,time,comm,ok)
    if (.not.ok) call fatal('Gas CSV output failed')
    write(suffix,'(i0)') rank
    particle_path=trim(config%output_file)//'.particles.'//trim(suffix)//'.csv'
    open(newunit=unit,file=trim(particle_path),status='replace',action='write',iostat=ios)
    if (ios/=0) call fatal('Cannot open particle CSV')
    header='id,parent_id,material_id,x,y,z,u,v,w,mass,temperature,multiplicity'
    do i=1,component_count
      write(suffix,'(i0)') i
      header=header//',component_mass_'//trim(suffix)
    end do
    write(unit,'(a)',iostat=ios) header
    if (ios/=0) call fatal('Cannot write particle CSV header')
    do p=1,size(parcels)
      write(unit,'(i0,2(",",i0),",",*(es25.16e3,:,","))',iostat=ios) parcels(p)%id,parcels(p)%parent_id, &
        parcels(p)%material_id,parcels(p)%position,parcels(p)%velocity,parcels(p)%mass,parcels(p)%temperature, &
        parcels(p)%multiplicity,masses(:,p)
      if (ios/=0) call fatal('Cannot write particle CSV')
    end do
    close(unit,iostat=ios)
    if (ios/=0) call fatal('Cannot close particle CSV')
    if (wall_active) call write_films()
    if (rank==0) then
      write(*,'(a,i0,a,es24.16,a,i0)') 'Steps: ',step,' Final time: ',time,' Last particle ID: ',last_id
      write(*,'(a,l1)') 'Regridded: ',regridded
      write(*,'(a,i0)') 'Injection events applied: ',event-1
      write(*,'(a,*(es24.16,1x))') 'Injected inventory: ',global_ledger
      write(*,'(a,*(es24.16,1x))') 'Boundary inventory: ',boundary
      if (allocated(force_total)) write(*,'(a,*(es24.16,1x))') 'Source inventory: ',force_total
    end if
  contains
    subroutine write_films()
      integer :: f,c,out,status
      character(len=:), allocatable :: names,path
      write(suffix,'(i0)') rank
      path=trim(config%output_file)//'.films.'//trim(suffix)//'.csv'
      open(newunit=out,file=path,status='replace',action='write',iostat=status)
      if (status/=0) call fatal('Cannot open film CSV')
      names='material_id,x,y,z,nx,ny,nz,mass,temperature,diameter,height'
      do c=1,component_count
        write(suffix,'(i0)') c
        names=names//',component_mass_'//trim(suffix)
      end do
      write(out,'(a)',iostat=status) names
      if (status/=0) call fatal('Cannot write film CSV header')
      do f=1,size(films)
        write(out,'(i0,",",*(es25.16e3,:,","))',iostat=status) films(f)%material_id,films(f)%position, &
          films(f)%normal,films(f)%mass,films(f)%temperature,films(f)%diameter,films(f)%height,film_mass(:,f)
        if (status/=0) call fatal('Cannot write film CSV')
      end do
      close(out,iostat=status)
      if (status/=0) call fatal('Cannot close film CSV')
    end subroutine
    subroutine inlet_timestep()
      ! Recompute on the current tree: adaptation may change the finest spacing.
      inlet_dt=huge(1._dp)
      if (.not.any(faces=='inflow')) return
      call amr_tree_integer_layout_3d(tree,origin,domain,ratios,ok)
      if (.not.ok) call fatal('Invalid inlet timestep layout')
      min_h=spacing*real(domain(:,1),dp)/real(domain(:,size(ratios)),dp)
      do face=1,6
        if (faces(face)/='inflow') cycle
        ! Recover sound speed from the already resolved reservoir composition.
        q=0;q(1)=inlet(1,face);q(2:4)=config%inflow_velocity(:,face)
        q(5)=config%inflow_pressure(face);q(6:)=inlet(6:,face)/inlet(1,face)
        call reactive_primitive_to_conserved(species,q,change,inlet_t(face),sound,ok)
        if (.not.ok) call fatal('Invalid inlet CFL state')
        inlet_dt=min(inlet_dt,config%cfl/sum((abs(q(2:4))+sound)/min_h))
        if (.not.config%transport_enabled) cycle
        do k=1,2
          do j=1,2
            do i=1,2
              inlet_cube(:,i,j,k)=inlet(:,face)
            end do
          end do
        end do
        inlet_temperature=inlet_t(face)
        call reactive_transport_timestep_3d(species,transport,inlet_cube,inlet_temperature,2,2,2, &
          min_h(1),min_h(2),min_h(3),config%transport_cfl,config%viscosity_enabled, &
          config%thermal_conduction_enabled,config%species_diffusion_enabled,face_dt,diffusivity,ok)
        if (.not.ok) call fatal('Invalid inlet transport bound')
        inlet_dt=min(inlet_dt,face_dt)
      end do
    end subroutine
    subroutine checkpoint(writing,path)
      logical, intent(in) :: writing
      character(len=*), intent(in) :: path
      if (wall_active) then
        call mpi_amr_spray_checkpoint(comm,writing,path,context,species,tree,owners,state,temperature,parcels,masses, &
          time,step,event,regridded,last_id,global_ledger,lower,spacing,maximum_parcels,ok,periodic,boundary, &
          wall_films=films,film_mass=film_mass,maximum_films=walls%maximum_films,source=force_total)
      else
        call mpi_amr_spray_checkpoint(comm,writing,path,context,species,tree,owners,state,temperature,parcels,masses, &
          time,step,event,regridded,last_id,global_ledger,lower,spacing,maximum_parcels,ok,periodic,boundary,source=force_total)
      end if
      if (.not.ok) call fatal('AMR spray checkpoint operation failed')
    end subroutine
    function tree_key(geometry) result(key)
      type(amr_patch_tree_3d), intent(in) :: geometry
      character(len=:), allocatable :: key
      character(len=512) :: record
      integer :: peer
      write(record,'(*(i0,1x))') size(geometry%nodes),geometry%nodes(1)%shape
      key=trim(record)//achar(10)
      do peer=2,size(geometry%nodes)
        associate(item=>geometry%nodes(peer),box=>geometry%nodes(peer)%patch)
          write(record,'(*(i0,1x))') item%parent,box%coarse_i_lower,box%coarse_i_upper, &
            box%coarse_j_lower,box%coarse_j_upper,box%coarse_k_lower,box%coarse_k_upper,box%refinement_ratio
        end associate
        key=key//trim(record)//achar(10)
      end do
    end function
    subroutine inject_event()
      type(spray_parcel), allocatable :: candidate(:)
      real(dp), allocatable :: candidate_mass(:,:)
      real(dp) :: added(nv),next_ledger(nv)
      integer(int64) :: next_id
      integer :: outlet,old_count,born
      outlet=event_ports(event);old_count=size(parcels)
      candidate=parcels;candidate_mass=masses;next_id=last_id
      number=0;local_mass=0
      if (rank==0) then
        number=event_parcels(event);local_mass=event_mass(event)
      end if
      call inject_mpi_cone_parcels(comm,liquids(1),local_mass,number,port_centers(:,outlet),port_axes(:,outlet), &
        event_speeds(event),port_half_angles(outlet),port_diameters(outlet),port_temperatures(outlet), &
        maximum_parcels,next_id,candidate,ledger,ok,liquids,port_fractions(:component_count,outlet),candidate_mass, &
        nozzle_radius=port_radii(outlet),nozzle_model=port_nozzles(outlet), &
        diameter_shape=port_sizes(outlet)%shape,diameter_min=port_sizes(outlet)%minimum, &
        diameter_max=port_sizes(outlet)%maximum,diameter_weighting=port_sizes(outlet)%weighting, &
        diameter_distribution=port_sizes(outlet)%distribution,diameter_std=port_sizes(outlet)%std)
      if (.not.ok) call fatal('Scheduled injection failed')
      do born=old_count+1,size(candidate)
        candidate(born)%material_id=outlet
      end do
      call migrate_amr_tree_spray_parcels(comm,tree,owners,candidate,lower,spacing,periodic, &
        maximum_parcels,ok,candidate_mass)
      if (.not.ok) call fatal('Scheduled particle migration failed')
      call MPI_Allreduce(ledger,added,nv,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      if (ierr/=MPI_SUCCESS) call fatal('Injection ledger reduction failed')
      if (any(abs(global_ledger)>huge(1._dp)-abs(added))) call fatal('Injection ledger overflow')
      next_ledger=global_ledger+added
      if (.not.all(ieee_is_finite(next_ledger))) call fatal('Nonfinite injection ledger')
      call move_alloc(candidate,parcels);call move_alloc(candidate_mass,masses)
      last_id=next_id;global_ledger=next_ledger
    end subroutine
    logical function paths_alias(left,right) result(same)
      character(len=*), intent(in) :: left,right
      same=selected_reactor_paths_alias_lexically(absolute_path(left),absolute_path(right))
    end function
    function absolute_path(path) result(value)
      character(len=*), intent(in) :: path
      character(len=:), allocatable :: value
      character(kind=c_char) :: buffer(4096)
      type(c_ptr) :: pointer
      integer :: index
      if (len_trim(path)==0) call fatal('Empty MPI tree path')
      if (path(1:1)=='/') then
        value=trim(path)
      else
        pointer=system_getcwd(buffer,int(size(buffer),c_size_t))
        if (.not.c_associated(pointer)) call fatal('Cannot resolve MPI tree working directory')
        value=''
        do index=1,size(buffer)
          if (buffer(index)==c_null_char) exit
          value=value//buffer(index)
        end do
        value=value//'/'//trim(path)
      end if
      if (len(value)>1024) call fatal('MPI tree path exceeds lexical comparison capacity')
    end function
    subroutine check_checkpoint_path(path)
      character(len=*), intent(in) :: path
      character(len=24) :: rank_text
      character(len=:), allocatable :: output
      call check_output(path)
      if (paths_alias(path,trim(config%output_file))) &
        call fatal('Checkpoint aliases gas output')
      write(rank_text,'(i0)') rank
      output=trim(config%output_file)//'.particles.'//trim(rank_text)//'.csv'
      if (paths_alias(path,output)) call fatal('Checkpoint aliases particle output')
      output=trim(config%output_file)//'.films.'//trim(rank_text)//'.csv'
      if (paths_alias(path,output)) call fatal('Checkpoint aliases film output')
    end subroutine
    subroutine check_output(path)
      character(len=*), intent(in) :: path
      if (paths_alias(path,input_path).or. &
          paths_alias(path,topology_path).or. &
          paths_alias(path,spray_path)) call fatal('Output aliases input')
      if (len_trim(regrid_topology)>0) then
        if (paths_alias(path,trim(regrid_topology))) call fatal('Output aliases regrid input')
      end if
    end subroutine
    subroutine fatal(message)
      character(len=*), intent(in) :: message
      integer :: status
      write(*,'(a,i0,2a)') 'AMR spray rank ',rank,': ',message
      call MPI_Abort(comm,2,status)
      error stop 2
    end subroutine
    subroutine read_topology(path,result)
      character(len=*), intent(in) :: path
      type(amr_patch_tree_3d), intent(out) :: result
      type(amr_patch_3d), allocatable :: patches(:)
      integer, allocatable :: parents(:),shapes(:,:)
      integer :: input,ios,nchild,peer,parent,dims(3),lower(3),upper(3),ratio
      logical :: valid
      open(newunit=input,file=path,status='old',action='read',iostat=ios)
      if (ios/=0) call fatal('Cannot open tree topology')
      read(input,*,iostat=ios) nchild
      if (ios/=0) call fatal('Invalid tree child count')
      if (nchild<0.or.nchild==huge(nchild)) call fatal('Invalid tree child count')
      allocate(patches(nchild),parents(nchild),shapes(3,nchild+1),stat=ios)
      if (ios/=0) call fatal('Cannot allocate tree topology')
      shapes(:,1)=[config%nx,config%ny,config%nz]
      do peer=1,nchild
        read(input,*,iostat=ios) parent,lower(1),upper(1),lower(2),upper(2),lower(3),upper(3),ratio
        if (ios/=0) call fatal('Invalid tree patch record')
        if (parent<1.or.parent>peer) call fatal('Tree parents must precede children')
        parents(peer)=parent;dims=shapes(:,parent)
        call initialize_amr_patch_3d(dims(1),dims(2),dims(3),lower(1),upper(1),lower(2),upper(2), &
          lower(3),upper(3),ratio,patches(peer),valid)
        if (.not.valid) call fatal('Invalid tree patch geometry')
        if (any(int(upper-lower+1,int64)*int(ratio,int64)>int(huge(1),int64))) call fatal('Tree patch extent overflow')
        shapes(:,peer+1)=[patches(peer)%fine_nx(),patches(peer)%fine_ny(),patches(peer)%fine_nz()]
      end do
      close(input)
      call initialize_amr_patch_tree_3d(shapes(:,1),parents,patches,result,valid)
      if (.not.valid) call fatal('Invalid tree topology')
    end subroutine
  end subroutine
end module
