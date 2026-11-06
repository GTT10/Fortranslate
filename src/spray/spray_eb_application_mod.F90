module spray_eb_application_mod
  use selected_composition_mod, only: resolve_selected_composition
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use precision_mod, only: dp
  use, intrinsic :: iso_c_binding, only: c_char,c_size_t,c_ptr,c_associated,c_null_char
  use simulation_config_selected_reactor_mod, only: selected_reactor_paths_alias_lexically
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species,configure_nasa7_eos_policy,nasa7_eos_policy_identity
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim, &
    reactive_conserved_to_primitive,reactive_primitive_to_conserved
  use simulation_config_reactive_3d_mod, only: fill_reactive_boundary_ghosts_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config,build_configured_eb_geometry_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d,write_reactive_eb_3d_csv, &
    reactive_eb_integrals_3d
  use reactive_eb_cfl_3d_mod, only: compute_reactive_eb_cfl_timestep_3d
  use eb_reactive_transport_3d_mod, only: reactive_eb_transport_timestep_3d, &
    reactive_eb_halo_bounds_3d,reactive_eb_isothermal_timestep_3d
  use spray_parcel_mod
  use spray_injection_profile_mod, only: integrate_spray_injection_profile,read_spray_injection_profile
  use spray_mixture_mod, only: apply_mixture_breakup_viscosity, liquid_mixture_properties,mixture_evolution_id
  use spray_eb_reactive_mod, only: advance_eb_spray_reactive_allocated
  use spray_eb_impact_mod, only: spray_eb_wall_options
  use spray_wall_impact_mod, only: spray_wall_film
  use spray_eb_film_step_mod, only: advance_eb_film_coupling
  use spray_eb_step_mod, only: advance_eb_mixture_coupling_allocated
  use spray_eb_field_mod, only: eb_spray_stencil
  use spray_checkpoint_mod, only: write_spray_checkpoint,read_spray_checkpoint,spray_checkpoint_context_capacity
  use amr_reactive_3d_checkpoint_mod, only: write_selected_checkpoint_context
  use gas_transport_mod, only: write_transport_context
  implicit none
  private
  public :: run_eb_spray_application
  interface
    function system_getcwd(buffer,length) bind(C,name='getcwd') result(pointer)
      import c_char,c_size_t,c_ptr
      character(kind=c_char), intent(out) :: buffer(*)
      integer(c_size_t), value :: length
      type(c_ptr) :: pointer
    end function
  end interface
contains
  subroutine run_eb_spray_application(input_path,spray_path,config,species,reactions,transport,mole_fractions,integrator)
    character(len=*), intent(in) :: input_path,spray_path
    type(reactive_eb_3d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    type(nasa7_species), allocatable :: eos_species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), optional, intent(in) :: mole_fractions(:)
    character(len=*), optional, intent(in) :: integrator
    integer :: port_count=1,event_ports(1024)=1,port_numbers(16),port
    real(dp) :: port_centers(3,16),port_axes(3,16),port_masses(16),port_speeds(16)
    real(dp) :: port_angles(16),port_diameters(16),port_temperatures(16),port_fractions(16,16)
    type(spray_liquid) :: port_liquids(16)
    integer :: nc=1,number=1,injection_events=1
    real(dp) :: injection_times(1024)=0,injection_scales(1024)=1
    character(len=24) :: vapor_names(16)='H2O'
    real(dp) :: densities(16)=1000,heat_capacities(16)=4184,reference_temperatures(16)=373.15_dp
    real(dp) :: reference_pressures(16)=101325,latent_heats(16)=2.257e6_dp,fractions(16)=0
    real(dp) :: total_mass=1e-10_dp,center(3)=0,axis(3)=[1._dp,0._dp,0._dp],speed=1,half_angle=0
    real(dp) :: diameter=1e-5_dp,liquid_temperature=350,maximum_dt=1e-8_dp
    real(dp) :: nozzle_radius=0
    type(spray_nozzle_model) :: nozzle
    type(spray_size_options) :: sizes,port_sizes(16)
    type(spray_nozzle_model) :: port_nozzles(16)
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
    real(dp) :: mixture_normal_boiling_temperatures(16)
    real(dp) :: mixture_critical_temperatures(16)=0
    integer :: port_profile_points(16)
    real(dp) :: port_profile_times(64,16),port_profile_rates(64,16),port_radii(16),port_discharge_coefficients(16)
    integer :: injection_profile_points=0
    character(len=1024) :: profile_file,profile_files(16)
    namelist /eb_spray_rate_files/ profile_file,profile_files
    real(dp) :: injection_profile_times(64)=0,injection_profile_rates(64)=0,injection_end_times(1024)=-1
    logical :: velocity_from_rate=.false.
    real(dp) :: maximum_injection_speed,particle_dt
    real(dp) :: discharge_coefficient=1,event_masses(1024),event_speeds(1024),profile_mass,denominator
    type(spray_options) :: options
    type(spray_eb_wall_options) :: walls
    real(dp) :: acceleration(3)=0
    real(dp), allocatable :: force(:),sf(:)
    logical :: wall_impacts=.false.,particle_outflow=.false.
    type(spray_wall_film), allocatable :: films(:)
    real(dp), allocatable :: film_mass(:,:)
    namelist /eb_spray/ nc,number,vapor_names,densities,heat_capacities,reference_temperatures, &
      reference_pressures,latent_heats,fractions,total_mass,center,axis,speed,half_angle,diameter, &
      liquid_temperature,maximum_dt,options,injection_events,injection_times,injection_scales,wall_impacts,walls, &
      nozzle_radius,nozzle,sizes,mixture_density_models,mixture_density_coefficients,mixture_density_temperature_ranges, &
      mixture_cp_models, &
      mixture_viscosity_models,mixture_viscosities,mixture_viscosity_coefficients, &
      mixture_viscosity_temperature_ranges, &
      mixture_conductivity_models,mixture_conductivities,mixture_conductivity_coefficients, &
      mixture_conductivity_temperature_ranges,mixture_cp_coefficients,mixture_cp_temperature_ranges, &
      mixture_vapor_pressure_models,mixture_antoine_coefficients,mixture_antoine_temperature_ranges, &
      mixture_boiling_models,mixture_molecular_weights, &
      mixture_normal_boiling_temperatures, &
      mixture_critical_temperatures,injection_profile_points, &
      injection_profile_times,injection_profile_rates,injection_end_times,velocity_from_rate,discharge_coefficient, &
      port_count,event_ports,port_numbers,port_centers,port_axes,port_masses,port_speeds,port_angles, &
      port_diameters,port_temperatures,port_fractions,port_profile_points,port_profile_times,port_profile_rates, &
      port_radii,port_discharge_coefficients,port_nozzles,port_sizes,particle_outflow,acceleration
    type(eb_geometry_3d) :: g
    type(eb_geometry_3d), allocatable :: extended_geometry
    character(len=16), allocatable :: boundary_faces(:)
    real(dp) :: seed_mass(size(species)),inlet_mole(size(species))
    real(dp), allocatable :: inlet(:,:),inlet_t(:),primitive(:),ghost(:,:,:,:),ghost_t(:,:,:),wall_t(:)
    real(dp) :: inlet_sound,inlet_temp,donor_speeds(4),donor_diffusivity,donor_rate,wall_dt
    integer :: face
    type(reactive_eb_3d_config) :: canonical
    type(spray_liquid), allocatable :: liquids(:)
    type(spray_liquid) :: effective
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: phi(:,:,:),state(:,:,:,:),temperature(:,:,:),masses(:,:)
    real(dp), allocatable :: boundary(:),wall(:),chemical(:),sb(:),sw(:),sc(:),injected(:),initial(:)
    real(dp) :: time,dt,hydro_dt,transport_dt,diffusivity,theta,rho,cp,e
    integer :: unit,status,close_status,nv,n,p,step,attempt,particle_steps,next_event
    integer(int64) :: last_id,step_split_created
    integer(int64), allocatable :: split_created_total
    logical :: ok,model_rejected
    integer :: last_checkpoint_step
    character(len=:), allocatable :: context,checkpoint_path
    character(len=1024) :: message
    character(len=24) :: component_label
    character(len=:), allocatable :: header
    ! Configuration and namelist objects are local to one invocation.
    allocate(eos_species,source=species)
    call configure_nasa7_eos_policy(eos_species,config%pinned_eos_extrapolation,ok, &
      config%pinned_eos_temperature_min,config%pinned_eos_temperature_max)
    if(.not.ok) error stop 'Invalid EB spray EOS policy'
    port_profile_points=-1;port_profile_times=0;port_profile_rates=0;port_radii=-1;port_discharge_coefficients=-1
    port_nozzles=spray_nozzle_model(profile='inherit',swirl=huge(1._dp),spread=-1._dp)
    port_sizes=spray_size_options(distribution='inherit',weighting='inherit',std=-1._dp, &
      shape=-1._dp,minimum=-1._dp,maximum=-1._dp)
    port_count=1;event_ports=1;port_numbers=-1
    port_centers=huge(1._dp);port_axes=huge(1._dp);port_masses=-1;port_speeds=-1
    port_angles=-1;port_diameters=-1;port_temperatures=-1;port_fractions=-1
    nc=1;number=1;fractions=0;fractions(1)=1
    vapor_names='H2O';densities=1000;heat_capacities=4184;reference_temperatures=373.15_dp
    reference_pressures=101325;latent_heats=2.257e6_dp;total_mass=1e-10_dp;center=0
    axis=[1._dp,0._dp,0._dp];speed=1;half_angle=0;diameter=1e-5_dp;liquid_temperature=350
    maximum_dt=1e-8_dp;options=spray_options();options%particle_cfl=.5_dp;walls=spray_eb_wall_options();wall_impacts=.false.
    mixture_density_models='constant';mixture_cp_models='constant'
    mixture_conductivity_models='supplied';mixture_conductivities=0
    mixture_viscosity_models='supplied';mixture_viscosities=0
    mixture_conductivity_coefficients=0;mixture_conductivity_temperature_ranges=0
    mixture_viscosity_coefficients=0;mixture_viscosity_temperature_ranges=0
    mixture_vapor_pressure_models='clausius';mixture_boiling_models='reject'
    mixture_density_coefficients=0;mixture_density_temperature_ranges=0
    mixture_cp_coefficients=0;mixture_cp_temperature_ranges=0
    mixture_antoine_coefficients=0;mixture_antoine_temperature_ranges=0;mixture_critical_temperatures=0
    mixture_normal_boiling_temperatures=0
    mixture_molecular_weights=0
    nozzle_radius=0;nozzle=spray_nozzle_model();sizes=spray_size_options()
    injection_profile_points=0;injection_profile_times=0;injection_profile_rates=0;injection_end_times=-1
    velocity_from_rate=.false.;discharge_coefficient=1;particle_outflow=.false.;acceleration=0
    injection_events=1;injection_times=0;injection_scales=1;next_event=1;last_id=0
    call check_output(config%output_file)
    call check_output(trim(config%output_file)//'.particles.csv')
    call check_output(trim(config%output_file)//'.films.csv')
    if (config%geometry_type/='level_set') error stop 'EB spray requires nodal level-set geometry'
    if(len_trim(config%restart_file)>0) call check_checkpoint_path(config%restart_file)
    if(config%checkpoint_interval_steps>0) call check_checkpoint_path(config%checkpoint_file)
    open(newunit=unit,file=spray_path,status='old',action='read',iostat=status)
    if(status/=0) error stop 'Cannot open EB spray input'
    read(unit,nml=eb_spray,iostat=status)
    if(status/=0) error stop 'Invalid EB spray namelist'
    profile_file='';profile_files=''
    rewind(unit)
    read(unit,nml=eb_spray_rate_files,iostat=status)
    if(status>0) error stop 'Invalid EB spray rate files input'
    close(unit,iostat=close_status)
    if(close_status/=0) error stop 'Invalid EB spray namelist'
    if(port_count<1.or.port_count>16) error stop 'Invalid port count'
    if(any(profile_files(port_count+1:)/='')) error stop 'Unused EB rate files'
    if(len_trim(profile_file)>0) then
      if(injection_profile_points/=0.or.any(injection_profile_times/=0).or. &
        any(injection_profile_rates/=0)) error stop 'Rate file conflicts with inline profile'
      call read_spray_injection_profile(trim(profile_file),injection_profile_times, &
        injection_profile_rates,injection_profile_points,ok)
      if(.not.ok) error stop 'Invalid EB spray rate file'
    end if
    do port=1,port_count
      if(len_trim(profile_files(port))==0) cycle
      if((port_profile_points(port)/=-1.and.port_profile_points(port)/=0).or. &
        any(port_profile_times(:,port)/=0).or.any(port_profile_rates(:,port)/=0)) &
        error stop 'Rate file conflicts with inline profile'
      call read_spray_injection_profile(trim(profile_files(port)),port_profile_times(:,port), &
        port_profile_rates(:,port),port_profile_points(port),ok)
      if(.not.ok) error stop 'Invalid EB spray rate file'
    end do
    if(nc<1.or.nc>16.or.number<1.or.number>125000) error stop 'Invalid EB spray population or components'
    if(.not.ieee_is_finite(maximum_dt)) error stop 'Invalid EB spray maximum_dt'
    if(maximum_dt<=0.or.maximum_dt>1e6_dp) error stop 'Invalid EB spray maximum_dt'
    if(.not.all(ieee_is_finite(fractions)).or.any(fractions<0)) error stop 'Invalid liquid fractions'
    if(abs(sum(fractions(:nc))-1)>1e-12_dp.or.any(fractions(nc+1:)/=0)) error stop 'Liquid fractions must sum to one'
    if(options%maximum_multiplicity>0) then
      allocate(split_created_total);split_created_total=0
    end if
    if(injection_events<1.or.injection_events>1024) error stop 'Invalid injection_events'
    if(.not.all(ieee_is_finite(injection_times)).or..not.all(ieee_is_finite(injection_scales))) &
      error stop 'Nonfinite injection schedule'
    if(any(injection_times<0).or.any(injection_times>1e6_dp).or.any(injection_scales<0).or. &
       any(injection_scales>1e6_dp)) error stop 'Invalid injection times or scales'
    do n=2,injection_events
      if(injection_times(n)<injection_times(n-1)) error stop 'Injection times must be nondecreasing'
    end do
    if(any(config%boundary_faces/='legacy')) then
      if(.not.particle_outflow) error stop 'Explicit EB spray boundaries require particle_outflow'
    end if
    call build_configured_eb_geometry_3d(config,g,ok,phi)
    if(.not.ok) error stop 'Cannot build EB spray geometry'
    if(config%hydro_method=='plm_ctu') then
      allocate(extended_geometry)
      call build_configured_eb_geometry_3d(config,extended_geometry,ok,halo=6)
      if(.not.ok) error stop 'Cannot build extended EB spray geometry'
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
    do n=1,nc
      call initialize_spray_liquid(eos_species,vapor_names(n),densities(n),heat_capacities(n),reference_temperatures(n), &
        reference_pressures(n),latent_heats(n),liquids(n),ok,mixture_vapor_pressure_models(n), &
        mixture_antoine_coefficients(:,n),mixture_antoine_temperature_ranges(:,n))
      if(.not.ok) error stop 'Invalid liquid properties'
      liquids(n)%density_model=mixture_density_models(n)
      liquids(n)%density_coefficients=mixture_density_coefficients(:,n)
      liquids(n)%density_temperature_range=mixture_density_temperature_ranges(:,n)
      liquids(n)%heat_capacity_model=mixture_cp_models(n)
      liquids(n)%heat_capacity_coefficients=mixture_cp_coefficients(:,n)
      liquids(n)%heat_capacity_temperature_range=mixture_cp_temperature_ranges(:,n)
      liquids(n)%conductivity_model=mixture_conductivity_models(n)
      liquids(n)%viscosity_model=mixture_viscosity_models(n)
      liquids(n)%conductivity=mixture_conductivities(n)
      liquids(n)%viscosity=mixture_viscosities(n)
      liquids(n)%conductivity_coefficients=mixture_conductivity_coefficients(:,n)
      liquids(n)%viscosity_coefficients=mixture_viscosity_coefficients(:,n)
      liquids(n)%conductivity_temperature_range=mixture_conductivity_temperature_ranges(:,n)
      liquids(n)%viscosity_temperature_range=mixture_viscosity_temperature_ranges(:,n)
      liquids(n)%boiling_model=mixture_boiling_models(n)
      liquids(n)%critical_temperature=mixture_critical_temperatures(n)
      liquids(n)%normal_boiling_temperature=mixture_normal_boiling_temperatures(n)
      liquids(n)%molecular_weight=mixture_molecular_weights(n)
      if(.not.valid_spray_liquid(liquids(n))) error stop 'Invalid mixture component property model'
    end do
    if(all(liquids%viscosity_model/='supplied')) then
      call apply_mixture_breakup_viscosity(liquids,fractions(:nc),liquid_temperature,options,ok, &
        config%initial_pressure,eos_species(liquids%vapor_index)%molecular_weight)
      if(.not.ok) error stop 'Invalid component viscosity at injection'
    end if
    if(.not.valid_spray_options(options)) error stop 'Unsupported EB spray options'
    if(port_count<1.or.port_count>16) error stop 'Invalid injection port count'
    if(any(event_ports(:injection_events)<1).or.any(event_ports(:injection_events)>port_count)) &
      error stop 'Invalid injection event port'
    do port=1,port_count
      if(port_nozzles(port)%profile=='inherit') port_nozzles(port)%profile=nozzle%profile
      if(port_nozzles(port)%swirl==huge(1._dp)) port_nozzles(port)%swirl=nozzle%swirl
      if(port_nozzles(port)%spread==-1) port_nozzles(port)%spread=nozzle%spread
      if(port_sizes(port)%distribution=='inherit') port_sizes(port)%distribution=sizes%distribution
      if(port_sizes(port)%weighting=='inherit') port_sizes(port)%weighting=sizes%weighting
      if(port_sizes(port)%std==-1) port_sizes(port)%std=sizes%std
      if(port_sizes(port)%shape==-1) port_sizes(port)%shape=sizes%shape
      if(port_sizes(port)%minimum==-1) port_sizes(port)%minimum=sizes%minimum
      if(port_sizes(port)%maximum==-1) port_sizes(port)%maximum=sizes%maximum
      if(port_numbers(port)==-1) port_numbers(port)=number
      if(all(port_centers(:,port)==huge(1._dp))) port_centers(:,port)=center
      if(all(port_axes(:,port)==huge(1._dp))) port_axes(:,port)=axis
      if(port_masses(port)==-1) port_masses(port)=total_mass
      if(port_speeds(port)==-1) port_speeds(port)=speed
      if(port_angles(port)==-1) port_angles(port)=half_angle
      if(port_diameters(port)==-1) port_diameters(port)=diameter
      if(port_temperatures(port)==-1) port_temperatures(port)=liquid_temperature
      if(all(port_fractions(:,port)==-1)) port_fractions(:,port)=fractions
      if(port_numbers(port)<1.or.port_numbers(port)>min(125000,options%maximum_parcels)) &
        error stop 'Invalid port population'
      if(.not.all(ieee_is_finite(port_fractions(:,port))).or.any(port_fractions(:,port)<0)) &
        error stop 'Invalid port fractions'
      if(abs(sum(port_fractions(:nc,port))-1)>1e-12_dp.or.any(port_fractions(nc+1:,port)/=0)) &
        error stop 'Invalid port fraction sum'
      call liquid_mixture_properties(liquids,port_fractions(:nc,port),port_temperatures(port),rho,cp,e,ok)
      if(.not.ok) error stop 'Invalid injected liquid mixture'
      port_liquids(port)=spray_liquid(vapor_name=liquids(1)%vapor_name,vapor_index=liquids(1)%vapor_index, &
        density=rho,heat_capacity=cp,reference_temperature=port_temperatures(port), &
        reference_pressure=liquids(1)%reference_pressure,latent_heat=liquids(1)%latent_heat,reference_energy=e)
      if(.not.ieee_is_finite(port_masses(port))) error stop 'Invalid injection mass'
      if(port_masses(port)<=0.or.port_masses(port)>1e6_dp) error stop 'Invalid injection mass'
    end do
    if(injection_profile_points/=0.and.(injection_profile_points<2.or.injection_profile_points>64)) &
      error stop 'Invalid injection profile size'
    if(any(injection_profile_times(injection_profile_points+1:)/=0).or. &
       any(injection_profile_rates(injection_profile_points+1:)/=0)) error stop 'Unused injection profile entries'
    do port=1,port_count
      if(port_radii(port)==-1) port_radii(port)=nozzle_radius
      if(port_discharge_coefficients(port)==-1) port_discharge_coefficients(port)=discharge_coefficient
      if(port_profile_points(port)==-1) then
        if(any(port_profile_times(:,port)/=0).or.any(port_profile_rates(:,port)/=0)) &
          error stop 'Port profile data requires an explicit point count'
        port_profile_points(port)=injection_profile_points
        port_profile_times(:,port)=injection_profile_times;port_profile_rates(:,port)=injection_profile_rates
      end if
      n=port_profile_points(port)
      if(n/=0.and.(n<2.or.n>64)) error stop 'Invalid port profile size'
      if(any(port_profile_times(n+1:,port)/=0).or.any(port_profile_rates(n+1:,port)/=0)) &
        error stop 'Unused port profile entries'
      if(n>0) then
        call integrate_spray_injection_profile(port_profile_times(:n,port),port_profile_rates(:n,port), &
          port_profile_times(1,port),port_profile_times(n,port),profile_mass,ok)
        if(.not.ok) error stop 'Invalid port rate profile'
      end if
      if(.not.ieee_is_finite(port_discharge_coefficients(port))) error stop 'Invalid discharge coefficient'
      if(velocity_from_rate) then
        if(port_discharge_coefficients(port)<1e-6_dp.or.port_discharge_coefficients(port)>1) &
          error stop 'Invalid discharge coefficient'
        if(.not.ieee_is_finite(port_radii(port))) error stop 'Invalid rate-driven outlet radius'
        if(port_radii(port)<1e-12_dp.or.port_radii(port)>1) error stop 'Invalid rate-driven outlet radius'
        if(n<2) error stop 'Rate-driven velocity requires an injection profile'
      else if(port_discharge_coefficients(port)/=1) then
        error stop 'Discharge coefficient requires rate-driven velocity'
      end if
    end do
    event_masses=0;event_speeds=0
    do p=1,injection_events
      port=event_ports(p);n=port_profile_points(port);event_speeds(p)=port_speeds(port)
      if(n==0) then
        if(injection_end_times(p)/=-1) error stop 'Integration intervals require an injection profile'
        event_masses(p)=port_masses(port)*injection_scales(p)
      else
        call integrate_spray_injection_profile(port_profile_times(:n,port),port_profile_rates(:n,port), &
          injection_times(p),injection_end_times(p),profile_mass,ok)
        if(.not.ok) error stop 'Invalid injection integration interval'
        event_masses(p)=profile_mass*injection_scales(p)
      end if
      if(event_masses(p)>1e6_dp) error stop 'Integrated injection exceeds bounds'
      if(velocity_from_rate) then
        event_speeds(p)=0
        if(event_masses(p)>0) then
          denominator=(injection_end_times(p)-injection_times(p))*port_liquids(port)%density*acos(-1._dp)* &
            port_radii(port)**2*port_discharge_coefficients(port)
          if(denominator<=0) error stop 'Invalid rate-driven interval'
          if(event_masses(p)>1e7_dp*denominator) error stop 'Rate-driven speed exceeds bounds'
          event_speeds(p)=event_masses(p)/denominator
        end if
      end if
    end do
    maximum_injection_speed=maxval(event_speeds(:injection_events))
    do port=1,port_count
      if(.not.velocity_from_rate) then
        maximum_injection_speed=max(maximum_injection_speed,port_speeds(port))
        cycle
      endif
      n=port_profile_points(port)
      denominator=port_liquids(port)%density*acos(-1._dp)*port_radii(port)**2*port_discharge_coefficients(port)
      maximum_injection_speed=max(maximum_injection_speed,maxval(port_profile_rates(:n,port))/denominator)
      do p=1,injection_events
        if(event_ports(p)/=port) cycle
        maximum_injection_speed=max(maximum_injection_speed, &
          maxval(port_profile_rates(:n,port))*injection_scales(p)/denominator)
      enddo
    enddo
    allocate(drops(0),masses(nc,0),films(0),film_mass(nc,0))
    nv=reactive_nvar(size(eos_species))
    allocate(state(nv,g%nx,g%ny,g%nz),temperature(g%nx,g%ny,g%nz))
    if(.not.all(ieee_is_finite(acceleration))) error stop 'Invalid EB acceleration'
    if(any(abs(acceleration)>1e6_dp)) error stop 'Invalid EB acceleration'
    allocate(sf(nv));sf=0
    if(any(acceleration/=0).or.options%fixed_parts) then
      allocate(force(nv));force=0
    end if
    allocate(boundary(nv),wall(nv),chemical(nv),sb(nv),sw(nv),sc(nv),injected(nv),initial(nv))
    call initialize_reactive_eb_density_sheet_3d(eos_species,config,g,state,temperature,ok,mole_fractions)
    if(.not.ok) error stop 'Invalid EB initial gas'
    if(any(config%boundary_faces/='legacy')) then
      allocate(boundary_faces,source=config%boundary_faces)
      if(any(boundary_faces=='isothermal_wall')) allocate(wall_t,source=config%wall_temperature)
      allocate(inlet(nv,6),inlet_t(6),primitive(reactive_nprim(size(eos_species))))
      allocate(ghost(nv,0:config%nx+1,0:config%ny+1,0:config%nz+1), &
        ghost_t(0:config%nx+1,0:config%ny+1,0:config%nz+1))
      call reactive_conserved_to_primitive(eos_species,state(:,1,1,1),temperature(1,1,1), &
        primitive,inlet_temp,inlet_sound,ok)
      if(.not.ok) error stop 'Failed to obtain EB inlet composition'
      seed_mass=primitive(6:)
      do face=1,6
        inlet(:,face)=state(:,1,1,1);inlet_t(face)=temperature(1,1,1)
        if(boundary_faces(face)/='inflow') cycle
        primitive(6:)=seed_mass
        if(config%inflow_composition_count(face)/=0) then
          call resolve_selected_composition('EB spray inlet',config%inflow_composition_count(face), &
            config%inflow_composition_species(:,face),config%inflow_composition_mole_fractions(:,face), &
            eos_species,inlet_mole,ok,message)
          if(.not.ok) then
            write(*,'(a)') trim(message)
            error stop 'Invalid EB spray inlet composition'
          end if
          call mass_fractions_from_mole_fractions(eos_species,inlet_mole,primitive(6:),ok)
          if(.not.ok) error stop 'Invalid EB spray inlet mass fractions'
        end if
        primitive(1)=config%inflow_density(face);primitive(2:4)=config%inflow_velocity(:,face)
        primitive(5)=config%inflow_pressure(face)
        call reactive_primitive_to_conserved(eos_species,primitive,inlet(:,face),inlet_t(face),inlet_sound,ok)
        if(.not.ok) error stop 'Invalid prescribed EB inlet thermodynamic state'
      end do
    end if
    if (walls%evolve.and..not.wall_impacts) error stop 'Film evolution requires wall_impacts'
    if (wall_impacts) then
      if (options%surface_tension<1e-8_dp.or.liquid_viscosity_at_temperature(options,liquid_temperature)<=0) &
        error stop 'EB wall impacts require positive surface tension and viscosity'
      call advance_eb_mixture_coupling_allocated(eos_species,transport,liquids,options,g,phi,state,temperature, &
        drops,masses,0._dp,particle_steps,sw(:5),ok,walls,films,film_mass,last_id)
      if (.not.ok) error stop 'Invalid EB impact configuration'
    end if
    call reactive_eb_integrals_3d(state,g,initial,ok)
    if(.not.ok) error stop 'Invalid initial EB gas inventory'
    time=0;step=0;boundary=0;wall=0;chemical=0;injected=0;last_checkpoint_step=0
    if(config%checkpoint_interval_steps>0.or.len_trim(config%restart_file)>0) call build_context()
    if(len_trim(config%restart_file)>0) then
      call read_spray_checkpoint(config%restart_file,context,eos_species,state,temperature,drops,time,step,next_event, &
        initial,injected,ok,message,boundary_totals=boundary,source_totals=wall,chemistry_totals=chemical, &
        breakup_state=.true.,split_state=.true.,khrt_state=.true.,material_state=.true.,shedding_state=.true., &
        component_mass=masses,last_id=last_id,film_state=.true.,wall_films=films,film_component_mass=film_mass, &
        retired_parent_state=wall_impacts,split_created_total=split_created_total,external_totals=force)
      if(.not.ok) then
        write(*,'(a)') trim(message)
        error stop 'EB spray restart rejected'
      end if
      call validate_restart()
      last_checkpoint_step=step
    end if
    do
      do while(next_event<=injection_events)
        if(injection_times(next_event)>time) exit
        call inject_event()
        next_event=next_event+1
      end do
      if(config%checkpoint_interval_steps>0.and.step>last_checkpoint_step) then
        if(modulo(step,config%checkpoint_interval_steps)==0) then
          write(component_label,'(i0)') step
          checkpoint_path=trim(config%checkpoint_file)//'.step'//trim(component_label)
          call check_checkpoint_path(checkpoint_path)
          call write_spray_checkpoint(checkpoint_path,context,eos_species,state,temperature,drops,time,step,next_event, &
            initial,injected,ok,message,boundary_totals=boundary,source_totals=wall,chemistry_totals=chemical, &
            breakup_state=.true.,split_state=.true.,khrt_state=.true.,material_state=.true.,shedding_state=.true., &
            component_mass=masses,last_id=last_id,film_state=.true.,wall_films=films,film_component_mass=film_mass, &
        retired_parent_state=wall_impacts,split_created_total=split_created_total,external_totals=force)
          if(.not.ok) then
            write(*,'(a)') trim(message)
            error stop 'EB spray checkpoint rejected'
          end if
          last_checkpoint_step=step
          if(config%stop_after_checkpoint) exit
        end if
      end if
      if(time>=config%final_time) exit
      if(step>=config%maximum_steps) error stop 'EB spray maximum_steps exhausted'
      call compute_reactive_eb_cfl_timestep_3d(eos_species,state,temperature,g,config%cfl,hydro_dt,ok)
      if(.not.ok) error stop 'Invalid EB gas CFL bound'
      call compute_spray_timestep(drops,[g%dx,g%dy,g%dz],options%particle_cfl,options%fixed_parts, &
        maximum_injection_speed,particle_dt,ok)
      if(.not.ok) error stop "Invalid EB particle timestep"
      dt=min(hydro_dt,maximum_dt,particle_dt,config%final_time-time)
      if(next_event<=injection_events) dt=min(dt,injection_times(next_event)-time)
      if(config%transport_enabled) then
        call reactive_eb_transport_timestep_3d(eos_species,transport,state,temperature,g,config%transport_cfl, &
          config%viscosity_enabled,config%thermal_conduction_enabled,config%species_diffusion_enabled, &
          transport_dt,diffusivity,ok)
        if(.not.ok) error stop 'Invalid EB transport bound'
        dt=min(dt,transport_dt)
      end if
      if(allocated(boundary_faces)) then
        call fill_reactive_boundary_ghosts_3d(state,temperature,boundary_faces,ghost,ghost_t,ok, &
          viscous_wall=.true.,inflow_state=inlet,inflow_temperature=inlet_t,ghost_width=1)
        if(.not.ok) error stop 'Failed to prepare EB timestep boundary donors'
        call reactive_eb_halo_bounds_3d(eos_species,transport,g,ghost,ghost_t,config%transport_enabled, &
          config%viscosity_enabled,config%thermal_conduction_enabled,config%species_diffusion_enabled, &
          donor_speeds,donor_diffusivity,ok)
        if(.not.ok) error stop 'Failed to compute EB boundary timestep bounds'
        donor_rate=sum(donor_speeds(:3)/[g%dx,g%dy,g%dz])
        if(.not.ieee_is_finite(donor_rate)) error stop 'Invalid EB boundary signal rate'
        if(donor_rate>0) dt=min(dt,config%cfl/donor_rate)
        if(config%transport_enabled.and.donor_diffusivity>0) then
          donor_rate=donor_diffusivity*sum(1._dp/[g%dx,g%dy,g%dz]**2)
          if(.not.ieee_is_finite(donor_rate)) error stop 'Invalid EB boundary diffusion rate'
          dt=min(dt,config%transport_cfl/donor_rate)
        end if
      end if
      if(allocated(wall_t)) then
        call reactive_eb_isothermal_timestep_3d(eos_species,transport,state,temperature,g, &
          boundary_faces,wall_t,config%transport_cfl,wall_dt,ok)
        if(.not.ok) error stop 'Invalid EB spray wall thermal timestep'
        dt=min(dt,wall_dt)
      end if
      do attempt=1,30
        ok=.false.
        if(dt<=0.or.time+dt<=time) exit
        theta=1;sb=0;sw=0;sc=0
        if (wall_impacts) then
        call advance_eb_spray_reactive_allocated(eos_species,reactions,transport,liquids,options,g,phi,state,temperature, &
          drops,masses,dt,config%riemann_solver,config%redistribution,config%chemistry_enabled,config%transport_enabled, &
          config%viscosity_enabled,config%thermal_conduction_enabled,config%species_diffusion_enabled, &
          config%barodiffusion_enabled,config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance, &
          particle_steps,theta,sb,sw,sc,ok,config%state_redist_target_volume_fraction,integrator, &
          walls,films,film_mass,last_id,model_rejected,step_split_created,particle_outflow,acceleration,sf, &
          difmag=config%difmag,extended_geometry=extended_geometry,boundary_faces=boundary_faces, &
          inflow_state=inlet,inflow_temperature=inlet_t,wall_temperature=wall_t, &
          state_redist_max_order=config%state_redist_max_order)
          if(model_rejected) error stop 'EB impact correlation outside admissible mass/velocity range'
        else
        call advance_eb_spray_reactive_allocated(eos_species,reactions,transport,liquids,options,g,phi,state,temperature, &
          drops,masses,dt,config%riemann_solver,config%redistribution,config%chemistry_enabled,config%transport_enabled, &
          config%viscosity_enabled,config%thermal_conduction_enabled,config%species_diffusion_enabled, &
          config%barodiffusion_enabled,config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance, &
          particle_steps,theta,sb,sw,sc,ok,config%state_redist_target_volume_fraction,integrator, &
          last_id=last_id,split_created=step_split_created,particle_outflow=particle_outflow, &
          acceleration=acceleration,source_change=sf, &
          difmag=config%difmag,extended_geometry=extended_geometry,boundary_faces=boundary_faces, &
          inflow_state=inlet,inflow_temperature=inlet_t,wall_temperature=wall_t, &
          state_redist_max_order=config%state_redist_max_order)
        end if
        if(ok) exit
        dt=.5_dp*dt
      end do
      if(.not.ok) error stop 'EB spray coupled step rejected'
      if(allocated(split_created_total)) split_created_total=split_created_total+step_split_created
      if(allocated(force)) force=force+sf
      time=time+dt;step=step+1;boundary=boundary+sb;wall=wall+sw;chemical=chemical+sc
    end do
    call write_reactive_eb_3d_csv(config%output_file,eos_species,g,state,temperature,time,ok,message)
    if(.not.ok) error stop 'Cannot write EB spray gas CSV'
    open(newunit=unit,file=trim(config%output_file)//'.particles.csv',status='replace',action='write',iostat=status)
    if(status/=0) error stop 'Cannot write EB spray particle CSV'
    header='id,material_id,x,y,z,u,v,w,mass,temperature,multiplicity'
    do n=1,nc
      write(component_label,'(i0)') n
      header=header//',component_mass_'//trim(component_label)
    end do
    write(unit,'(a)',iostat=status) header
    do p=1,size(drops)
      if(status/=0) exit
      write(unit,'(i0,",",i0,*(",",es25.16,:))',iostat=status) &
        drops(p)%id,drops(p)%material_id,drops(p)%position,drops(p)%velocity, &
        drops(p)%mass,drops(p)%temperature,drops(p)%multiplicity,masses(:,p)
    end do
    close(unit,iostat=close_status)
    if(status/=0.or.close_status/=0) error stop 'EB particle CSV write failed'
    open(newunit=unit,file=trim(config%output_file)//'.films.csv',status='replace',action='write',iostat=status)
    if(status/=0) error stop 'Cannot write EB film CSV'
    header='material_id,x,y,z,nx,ny,nz,mass,temperature,diameter,height'
    do n=1,nc
      write(component_label,'(i0)') n
      header=header//',component_mass_'//trim(component_label)
    end do
    write(unit,'(a)',iostat=status) header
    do p=1,size(films)
      if(status/=0) exit
      write(unit,'(i0,*(",",es25.16,:))',iostat=status) films(p)%material_id,films(p)%position,films(p)%normal, &
        films(p)%mass,films(p)%temperature,films(p)%diameter,films(p)%height,film_mass(:,p)
    end do
    close(unit,iostat=close_status)
    if(status/=0.or.close_status/=0) error stop 'EB film CSV write failed'
    write(*,'(a,i0,a,es25.16)') 'EB spray steps: ',step,' time: ',time
    write(*,'(a,i0)') 'Injection events processed: ',next_event-1
    write(*,'(a,*(es25.16,1x))') 'Injected inventory: ',injected
    write(*,'(a,*(es25.16,1x))') 'Boundary inventory: ',boundary
    write(*,'(a,*(es25.16,1x))') 'Wall inventory: ',wall
    write(*,'(a,*(es25.16,1x))') 'Chemistry inventory: ',chemical
    if(allocated(force)) write(*,'(a,*(es25.16,1x))') 'Force inventory: ',force
  contains
    subroutine build_context()
      integer :: scratch,ios,j,k,length
      character(len=spray_checkpoint_context_capacity) :: line
      character(len=32) :: method
      real(dp) :: composition(size(eos_species))
      canonical=config
      canonical%final_time=0;canonical%maximum_steps=0;canonical%output_file=''
      canonical%checkpoint_interval_steps=0;canonical%checkpoint_file='';canonical%restart_file=''
      canonical%stop_after_checkpoint=.false.;canonical%level_set_file=''
      method='default';if(present(integrator)) method=integrator
      composition=0
      if(present(mole_fractions)) composition=mole_fractions
      open(newunit=scratch,status='scratch',form='formatted',delim='quote',iostat=ios)
      if(ios/=0) error stop 'Cannot prepare EB checkpoint context'
      write(scratch,'(a)',iostat=ios) 'EB_SPRAY_CONTEXT_V2',mixture_evolution_id
      ! Keep legacy zero-coefficient context independent of added config fields.
      if(ios==0) write(scratch,*,iostat=ios) &
        canonical%nx,canonical%ny,canonical%nz, &
        canonical%maximum_steps,canonical%x_lower,canonical%x_upper, &
        canonical%y_lower,canonical%y_upper,canonical%z_lower, &
        canonical%z_upper,canonical%final_time,canonical%cfl, &
        canonical%geometry_type,canonical%extrusion_axis,canonical%section_file, &
        canonical%level_set_file,canonical%plane_axis,canonical%plane_position, &
        canonical%riemann_solver,canonical%redistribution,canonical%thermo_model, &
        canonical%state_redist_target_volume_fraction,canonical%chemistry_enabled,canonical%chemistry_relative_tolerance, &
        canonical%chemistry_absolute_tolerance,canonical%transport_enabled,canonical%viscosity_enabled, &
        canonical%thermal_conduction_enabled,canonical%species_diffusion_enabled,canonical%barodiffusion_enabled, &
        canonical%transport_cfl,canonical%regular_density,canonical%cut_density, &
        canonical%initial_pressure,canonical%initial_velocity_x,canonical%initial_velocity_y, &
        canonical%initial_velocity_z,canonical%composition_count,canonical%composition_species, &
        canonical%composition_mole_fractions,canonical%checkpoint_interval_steps,canonical%stop_after_checkpoint, &
        canonical%checkpoint_file,canonical%restart_file,canonical%output_file
      if(ios==0.and.config%state_redist_max_order/=0) &
        write(scratch,'(a,1x,i0)',iostat=ios) 'state-redist-order-v1',config%state_redist_max_order
      if(ios==0.and.config%difmag>0) &
        write(scratch,'(a,1x,es27.18e3)',iostat=ios) 'difmag-v1',config%difmag
      if(ios==0.and.config%hydro_method/='pcm') &
        write(scratch,'(a,1x,a)',iostat=ios) 'hydro-v1',trim(config%hydro_method)
      if(ios==0.and.allocated(boundary_faces)) then
        write(scratch,'(a)',iostat=ios) 'boundary-v1'
        do j=1,6
          if(ios==0) write(scratch,'(a,5(1x,es27.18e3))',iostat=ios) trim(boundary_faces(j)), &
            config%inflow_density(j),config%inflow_pressure(j),config%inflow_velocity(:,j)
        end do
      end if
      if(ios==0.and.any(config%inflow_composition_count/=0)) then
        write(scratch,'(a,6(1x,i0))',iostat=ios) 'inlet-composition-v1',config%inflow_composition_count
        do j=1,6
          do k=1,config%inflow_composition_count(j)
            if(ios==0) write(scratch,'(a,1x,es27.18e3)',iostat=ios) &
              trim(config%inflow_composition_species(k,j)),config%inflow_composition_mole_fractions(k,j)
          end do
        end do
      end if
      if(ios==0.and.allocated(wall_t)) &
        write(scratch,'(a,6(1x,es27.18e3))',iostat=ios) 'wall-v1',wall_t
      if(ios==0) write(scratch,nml=eb_spray,iostat=ios)
      do j=1,size(eos_species)
        if(ios==0) write(scratch,*,iostat=ios) eos_species(j)%name,eos_species(j)%molecular_weight,eos_species(j)%temperature_min, &
          eos_species(j)%temperature_mid,eos_species(j)%temperature_max, &
          eos_species(j)%low_coefficients,eos_species(j)%high_coefficients
        if(ios==0.and.eos_species(j)%pinned_eos_extrapolation) &
          write(scratch,'(a)',iostat=ios) trim(nasa7_eos_policy_identity(eos_species))
        if(ios==0) call write_transport_context(scratch,transport(j),ios)
      end do
      if(ios==0) call write_selected_checkpoint_context(scratch,reactions,'EB_SPRAY_V1',method,composition, &
        config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance,ios)
      do k=0,g%nz
      do j=0,g%ny
        if(ios==0) write(scratch,'(*(es26.17e3,1x))',iostat=ios) phi(:,j,k)
      end do
      end do
      if(ios/=0) error stop 'Cannot serialize EB checkpoint context'
      rewind(scratch,iostat=ios)
      if(ios/=0) error stop 'Cannot rewind EB checkpoint context'
      context='PARTICLE_TIMESTEP_V1|'
      if(options%etab_parcel_factor>=0.or.options%khrt_parcel_factor>=0) context=context//'REF_BREAKUP_N0_V1|'
      if(wall_impacts) context=context//'SPLASH_N0_V2|'
      do
        read(scratch,'(a)',iostat=ios) line
        if(ios<0) exit
        if(ios/=0.or.len_trim(line)==len(line)) error stop 'EB checkpoint context line exceeds capacity'
        length=len(context)+len_trim(line)+1
        if(length>spray_checkpoint_context_capacity) error stop 'EB checkpoint context exceeds capacity'
        if(options%khrt_parcel_factor==-1.and.index(line,'%KHRT_PARCEL_FACTOR=')>0) cycle
        context=context//trim(line)//'|'
      end do
      close(scratch)
      if(options%breakup_model=='etab') then
        if(len(context)+len('ETAB_MINMASS_V2|')>spray_checkpoint_context_capacity) &
          error stop 'EB checkpoint context exceeds capacity'
        context=context//'ETAB_MINMASS_V2|'
      endif
      if(options%khrt_parcel_factor>=0) then
        if(len(context)+len('KHRT_MINMASS_V2|')>spray_checkpoint_context_capacity) &
          error stop 'EB checkpoint context exceeds capacity'
        context=context//'KHRT_MINMASS_V2|'
      endif
    end subroutine
    subroutine validate_restart()
      integer :: q,expected_count,cells(3,8)
      real(dp) :: weights(8),volumes(8)
      logical :: valid
      if(time>config%final_time.or.step>config%maximum_steps) error stop 'Restart exceeds continuation limits'
      if(next_event<1.or.next_event>injection_events+1) error stop 'Invalid restart event index'
      if(next_event<=injection_events) then
        if(injection_times(next_event)<=time) error stop 'Unprocessed restart injection event'
      end if
      if(next_event>1) then
        if(injection_times(next_event-1)>time) error stop 'Restart injection ahead of time'
      end if
      expected_count=0
      do q=1,next_event-1
        if(event_masses(q)>0) expected_count=expected_count+port_numbers(event_ports(q))
      end do
      step_split_created=0
      if(allocated(split_created_total)) step_split_created=split_created_total
      if(last_id<int(expected_count,int64)) error stop 'Invalid restart injection watermark'
      if(step_split_created<0.or.step_split_created>last_id-int(expected_count,int64)) &
        error stop 'Invalid restart split count'
      if (wall_impacts.or.options%khrt_children) then
        if(modulo(last_id-int(expected_count,int64)-step_split_created,4_int64)/=0) error stop 'Invalid restart child IDs'
      else
        if(last_id-int(expected_count,int64)/=step_split_created) error stop 'Invalid restart injection IDs'
      end if
      if (.not.wall_impacts) then
        if(size(films)/=0) error stop 'Unexpected restart films'
        if(int(size(drops),int64)/=last_id) error stop 'Invalid restart particle count'
      end if
      if(size(drops)>min(125000,options%maximum_parcels)) error stop 'Restart exceeds population capacity'
      if(any(shape(masses)/=[nc,size(drops)])) error stop 'Restart component shape mismatch'
      do q=1,size(films)
        if(films(q)%material_id>port_count) error stop 'Invalid restart film port'
      end do
      cells=0;weights=0;volumes=0
      do q=1,size(drops)
        if(drops(q)%material_id>port_count) error stop 'Invalid restart particle port'
        if (wall_impacts) then
          if(drops(q)%id>last_id) error stop 'Restart ID exceeds watermark'
          if(q>1) then
            if(drops(q)%id<=drops(q-1)%id) error stop 'Restart IDs not ordered'
          end if
        else
          if(drops(q)%id/=int(q,int64)) error stop 'Restart parcel ID sequence mismatch'
        end if
        call eb_spray_stencil(g,phi,drops(q)%position,cells,weights,volumes,valid, &
          config%boundary_faces(1:6:2)=='periodic')
        if(.not.valid) error stop 'Restart parcel outside fluid'
      end do
      if (wall_impacts) then
        call advance_eb_film_coupling(eos_species,transport,liquids,options,walls,g,phi,state,temperature, &
          films,film_mass,0._dp,particle_steps,sw(:5),valid,config%boundary_faces(1:6:2)=='periodic')
        if(.not.valid) error stop 'Invalid restart EB films'
      end if
    end subroutine
    subroutine check_checkpoint_path(path)
      character(len=*), intent(in) :: path
      call check_output(path)
      if(selected_reactor_paths_alias_lexically(absolute_path(path),absolute_path(config%output_file)).or. &
         selected_reactor_paths_alias_lexically(absolute_path(path), &
           absolute_path(trim(config%output_file)//'.particles.csv')).or. &
         selected_reactor_paths_alias_lexically(absolute_path(path), &
           absolute_path(trim(config%output_file)//'.films.csv'))) error stop 'Checkpoint aliases CSV output'
    end subroutine
    subroutine inject_event()
      type(spray_parcel), allocatable :: added(:),candidate(:)
      real(dp), allocatable :: candidate_mass(:,:)
      real(dp) :: change(nv),totals(5),weights(8),volumes(8)
      integer :: cells(3,8),q,j,old_count,ios,port,nadd
      logical :: valid
      if(event_masses(next_event)==0) return
      port=event_ports(next_event);nadd=port_numbers(port);effective=port_liquids(port)
      old_count=size(drops)
      if(old_count>min(125000,options%maximum_parcels)-nadd) error stop 'EB injection capacity exceeded'
      call inject_cone_parcels(effective,event_masses(next_event),nadd,last_id+1_int64, &
        port_centers(:,port),port_axes(:,port),event_speeds(next_event),port_angles(port), &
        port_diameters(port),port_temperatures(port),added,valid, &
        diameter_shape=port_sizes(port)%shape,diameter_min=port_sizes(port)%minimum,diameter_max=port_sizes(port)%maximum, &
        diameter_weighting=port_sizes(port)%weighting,diameter_distribution=port_sizes(port)%distribution, &
        diameter_std=port_sizes(port)%std, &
        nozzle_radius=port_radii(port),nozzle_model=port_nozzles(port), &
        reference_breakup=options%etab_parcel_factor>=0.or.options%khrt_parcel_factor>=0)
      if(.not.valid) error stop 'Invalid scheduled EB injection'
      change=0;cells=0;weights=0;volumes=0
      allocate(candidate(old_count+nadd),candidate_mass(nc,old_count+nadd),stat=ios)
      if(ios/=0) error stop 'Cannot allocate scheduled injection'
      candidate(:old_count)=drops;candidate_mass(:,:old_count)=masses
      do q=1,nadd
        call eb_spray_stencil(g,phi,added(q)%position,cells,weights,volumes,valid, &
          config%boundary_faces(1:6:2)=='periodic')
        if(.not.valid) error stop 'EB injection point is not in fluid'
        added(q)%material_id=port
        candidate(old_count+q)=added(q)
        candidate_mass(:,old_count+q)=added(q)%mass*port_fractions(:nc,port)
        totals=parcel_totals(added(q),effective);change(1:5)=change(1:5)+totals
        do j=1,nc
          change(5+liquids(j)%vapor_index)=change(5+liquids(j)%vapor_index)+ &
            candidate_mass(j,old_count+q)*added(q)%multiplicity
        end do
      end do
      if(.not.all(ieee_is_finite(change))) error stop 'Invalid injection inventory'
      call move_alloc(candidate,drops);call move_alloc(candidate_mass,masses)
      last_id=last_id+int(nadd,int64);injected=injected+change
    end subroutine
    function absolute_path(path) result(value)
      character(len=*), intent(in) :: path
      character(len=:), allocatable :: value
      character(kind=c_char) :: buffer(4096)
      type(c_ptr) :: pointer
      integer :: index
      if (len_trim(path)==0) error stop 'Empty EB path'
      if (path(1:1)=='/') then
        value=trim(path)
      else
        pointer=system_getcwd(buffer,int(size(buffer),c_size_t))
        if (.not.c_associated(pointer)) error stop 'Cannot resolve EB working directory'
        value=''
        do index=1,size(buffer)
          if (buffer(index)==c_null_char) exit
          value=value//buffer(index)
        end do
        value=value//'/'//trim(path)
      end if
      if (len(value)>1024) error stop 'EB path too long'
    end function
    subroutine check_output(path)
      character(len=*), intent(in) :: path
      if(selected_reactor_paths_alias_lexically(absolute_path(path),absolute_path(input_path)).or. &
         selected_reactor_paths_alias_lexically(absolute_path(path),absolute_path(spray_path)).or. &
         selected_reactor_paths_alias_lexically(absolute_path(path),absolute_path(config%level_set_file))) &
        error stop 'EB output aliases input'
    end subroutine
  end subroutine
end module
