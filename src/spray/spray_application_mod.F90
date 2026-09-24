module spray_application_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use constants_mod, only: pelef_version,mechanism_species_capacity
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species, compatible_transport_database
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions, mixture_specific_gas_constant
  use reactive_1d_mod, only: reactive_nvar, reactive_nprim, reactive_primitive_to_conserved
  use reactive_3d_mod, only: compute_reactive_cfl_timestep_3d
  use selected_composition_mod, only: resolve_selected_composition
  use reactive_transport_3d_mod, only: reactive_transport_timestep_3d,reactive_isothermal_timestep_3d
  use simulation_config_reactive_3d_mod, only: valid_wall_velocity_3d,fill_mach_inlet_boundary_ghosts_3d
  use spray_wall_impact_mod
  use spray_injection_profile_mod, only: integrate_spray_injection_profile,read_spray_injection_profile
  use spray_mixture_mod
  use spray_parcel_mod
  use spray_coupling_3d_mod, only: spray_system_totals, mixture_system_totals, valid_spray_boundary_faces
  use les_filter_3d_mod, only: les_filter_weights
  use spray_reactive_3d_mod, only: advance_spray_reactive_3d, spray_cvode_available
  use spray_checkpoint_mod, only: write_spray_checkpoint, read_spray_checkpoint
  implicit none
  private
  public :: run_spray_application
contains
  subroutine run_spray_application(input_path,species,reactions,transport,mechanism_id)
    character(len=*), intent(in) :: input_path,mechanism_id
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    integer :: mixture_count
    character(len=24) :: mixture_vapor_names(16)
    real(dp) :: mixture_fractions(16),mixture_densities(16),mixture_cps(16),mixture_reference_temperatures(16)
    real(dp) :: mixture_reference_pressures(16),mixture_latent_heats(16),port_mixture_fractions(16,16)
    character(len=16) :: mixture_density_models(16),mixture_cp_models(16)
    character(len=16) :: mixture_vapor_pressure_models(16),mixture_boiling_models(16)
    real(dp) :: mixture_density_coefficients(4,16),mixture_density_temperature_ranges(2,16)
    real(dp) :: mixture_cp_coefficients(4,16),mixture_cp_temperature_ranges(2,16)
    real(dp) :: mixture_antoine_coefficients(4,16),mixture_antoine_temperature_ranges(2,16)
    real(dp) :: mixture_molecular_weights(16)=0
    real(dp) :: mixture_normal_boiling_temperatures(16)
    real(dp) :: mixture_critical_temperatures(16)
    type(spray_liquid), allocatable :: mixture_components(:)
    real(dp), allocatable :: component_mass(:,:),film_component_mass(:,:)
    character(len=16) :: liquid_cp_model,port_cp_models(16)
    real(dp) :: liquid_cp_coefficients(4),liquid_cp_temperature_range(2)
    real(dp) :: port_cp_coefficients(4,16),port_cp_temperature_ranges(2,16)
    character(len=16) :: wall_conductivity_model,port_wall_conductivity_models(16)
    real(dp) :: wall_conductivity_coefficients(4),wall_conductivity_temperature_range(2)
    real(dp) :: port_wall_conductivities(16),port_wall_conductivity_coefficients(4,16)
    real(dp) :: port_wall_conductivity_temperature_ranges(2,16)
    logical :: port_wall_properties
    logical :: wall_film_evolution
    real(dp) :: wall_liquid_conductivity,minimum_film_height
    logical :: wall_impact_faces(6),impact_enabled,lineage_enabled
    real(dp) :: wall_contact_angle
    integer :: maximum_wall_films,film_unit
    type(spray_wall_options) :: wall_options
    type(spray_wall_film), allocatable :: wall_films(:)
    type(spray_parcel) :: film_identity
    logical :: port_cp_overrides
    logical :: port_breakup_properties
    real(dp) :: port_surface_tensions(16),port_liquid_viscosities(16)
    character(len=24) :: port_viscosity_models(16)
    real(dp) :: port_viscosity_coefficients(4,16),port_viscosity_temperature_ranges(2,16)
    type(spray_options), allocatable :: breakup_options(:)
    logical :: port_property_models
    character(len=16) :: port_density_models(16),port_vapor_pressure_models(16),port_boiling_models(16)
    real(dp) :: port_density_coefficients(4,16),port_density_temperature_ranges(2,16)
    real(dp) :: port_antoine_coefficients(4,16),port_antoine_temperature_ranges(2,16),port_critical_temperatures(16)
    logical :: multi_material
    character(len=24) :: port_vapor_names(16)
    real(dp) :: port_liquid_densities(16),port_liquid_cps(16),port_reference_temperatures(16)
    real(dp) :: port_reference_pressures(16),port_latent_heats(16)
    type(spray_liquid), allocatable :: materials(:)
    type(spray_liquid) :: material
    logical :: port_state_overrides
    real(dp) :: injection_port_speeds(16),injection_port_diameters(16),injection_port_temperatures(16)
    real(dp) :: injection_port_cone_angles(16)
    logical :: independent_port_schedules
    integer :: port_events(16),selected_port,completed,total_events
    real(dp) :: port_mass_flow_rates(16),port_injection_starts(16),port_injection_ends(16),port_injection_intervals(16)
    integer :: injection_port_count,parcels_per_event,port
    real(dp) :: injection_port_centers(3,16),injection_port_axes(3,16),injection_port_fractions(16)
    integer :: port_profile_points(16),profile_batches
    real(dp) :: port_profile_times(64,16),port_profile_rates(64,16),profile_dose
    integer :: injection_profile_points,expected_batches,profile_event
    character(len=1024) :: profile_file,profile_files(16)
    namelist /spray_rate_files/ profile_file,profile_files
    real(dp) :: injection_profile_times(64),injection_profile_rates(64),profile_start
    integer :: cells(3),parcel_count,maximum_steps,stop_after_steps,maximum_substeps
    integer :: unit,ios,i,j,k,s,steps,event,retry,history_unit,gas_unit,parcel_unit,nvar
    real(dp) :: domain_length(3),gas_temperature,gas_pressure,gas_velocity(3),shear_velocity
    real(dp) :: mole_amounts(mechanism_species_capacity),liquid_density,liquid_cp
    character(len=16) :: boiling_model
    real(dp) :: normal_boiling_temperature,port_normal_boiling_temperatures(16)
    real(dp) :: liquid_molecular_weight,port_molecular_weights(16)
    real(dp) :: critical_temperature
    character(len=16) :: liquid_density_model
    real(dp) :: liquid_density_coefficients(4),liquid_density_temperature_range(2)
    character(len=16) :: vapor_pressure_model
    real(dp) :: antoine_coefficients(4),antoine_temperature_range(2)
    real(dp) :: reference_temperature,reference_pressure,latent_heat
    character(len=16) :: nozzle_profile,port_nozzle_profiles(16)
    real(dp) :: swirl_angle,hollow_spread,port_swirl_angles(16),port_hollow_spreads(16)
    type(spray_nozzle_model) :: nozzle_model
    type(spray_nozzle_model), allocatable :: port_nozzle_models(:)
    logical :: port_nozzle_overrides
    logical :: khrt_children
    real(dp) :: khrt_reference_diameter_cubed
    logical :: velocity_from_rate
    real(dp) :: discharge_coefficient,port_discharge_coefficients(16),event_speed
    real(dp) :: nozzle_radius,port_nozzle_radii(16)
    logical :: port_radius_overrides
    character(len=16) :: diameter_weighting,diameter_distribution
    character(len=16) :: port_diameter_distributions(16),port_diameter_weightings(16)
    real(dp) :: port_diameter_stds(16),port_diameter_shapes(16),port_diameter_minima(16),port_diameter_maxima(16)
    type(spray_size_options), allocatable :: port_sizes(:)
    logical :: port_size_overrides
    real(dp) :: diameter_std
    real(dp) :: diameter_shape,diameter_min,diameter_max
    real(dp) :: drop_temperature,drop_diameter,initial_liquid_mass,injection_center(3),injection_axis(3)
    real(dp) :: injection_speed,cone_angle,mass_flow_rate,injection_start,injection_end,injection_interval
    character(len=16) :: les_model,les_test_filter
    integer :: les_test_ratio,les_grid_ratio
    real(dp) :: les_ci
    logical :: les_harmonic_mean
    real(dp), allocatable :: les_weights(:)
    real(dp) :: final_time,maximum_dt,cfl,rtol,atol,smagorinsky_constant,turbulent_prandtl,turbulent_schmidt
    real(dp) :: max_fractional_change,particle_cfl,minimum_diameter,time,dt,cfl_dt,next_event,tolerance,amount
    real(dp) :: spacing(3),sound,local_temperature,gas_constant,totals5(5),volume
    real(dp), allocatable :: state(:,:,:,:),temperature(:,:,:),x(:),y(:),primitive(:)
    real(dp), allocatable :: reference_totals(:),injected_totals(:),totals(:),expected_initial(:)
    type(spray_parcel), allocatable :: parcels(:),new_parcels(:)
    type(spray_liquid) :: liquid
    type(spray_options) :: options
    character(len=24) :: drag_model
    logical :: chemistry_enabled,transport_enabled,drag_enabled,heat_enabled,evaporation_enabled,abramzon_sirignano
    real(dp), allocatable :: boundary_totals(:),boundary_step(:)
    character(len=8) :: breakup_model
    character(len=24) :: liquid_viscosity_model
    real(dp) :: liquid_viscosity_coefficients(4),liquid_viscosity_temperature_range(2)
    real(dp) :: surface_tension,liquid_viscosity
    real(dp) :: maximum_multiplicity
    integer :: maximum_parcels
    logical :: split_enabled
    real(dp) :: khrt_b0,khrt_b1,khrt_c3
    logical :: khrt_enabled
    logical :: breakup_enabled
    real(dp) :: outlet_pressure(6),inflow_mach(6),mach_pressure(6)
    real(dp), allocatable :: inlet_ghost(:,:,:,:),inlet_ghost_t(:,:,:)
    real(dp) :: wall_temperature(6),wall_velocity(3,6),wall_dt
    real(dp) :: acceleration(3)
    real(dp), allocatable :: source_totals(:),source_step(:)
    logical :: physical_boundaries
    real(dp) :: inflow_temperature(6),inflow_pressure(6),inflow_velocity(3,6),inflow_t(6),inlet_dt,face_dt,face_diffusivity
    real(dp), allocatable :: inflow_state(:,:),inlet_x(:),inlet_y(:),inlet_cube(:,:,:,:)
    real(dp) :: inlet_temperature(2,2,2),inflow_composition_mole_fractions(mechanism_species_capacity)
    character(len=24) :: inflow_composition_species(mechanism_species_capacity)
    integer :: inflow_composition_count,face
    character(len=32) :: boundary_faces(6)
    logical :: ok,exists
    character(len=24) :: mole_names(mechanism_species_capacity),vapor_name,chemistry_method,hydro_reconstruction,riemann_solver
    character(len=1024) :: output_prefix,checkpoint_file,restart_file,message
    character(len=16384) :: context
    character(len=1050), allocatable :: paths(:)
    integer(int64) :: first_id
    namelist /spray/ port_mixture_fractions,mixture_count,mixture_vapor_names,mixture_fractions,mixture_densities,mixture_cps, &
      mixture_reference_temperatures,mixture_reference_pressures,mixture_latent_heats, &
      mixture_density_models,mixture_density_coefficients,mixture_density_temperature_ranges, &
      mixture_cp_models,mixture_cp_coefficients,mixture_cp_temperature_ranges, &
      mixture_vapor_pressure_models,mixture_antoine_coefficients,mixture_antoine_temperature_ranges, &
      mixture_boiling_models,mixture_molecular_weights, &
      mixture_normal_boiling_temperatures, &
      mixture_critical_temperatures, &
      cells,domain_length,boundary_faces,wall_temperature,wall_velocity,acceleration,outlet_pressure,inflow_mach, &
      breakup_model,surface_tension,liquid_viscosity,khrt_children,khrt_reference_diameter_cubed, &
      liquid_viscosity_model,liquid_viscosity_coefficients,liquid_viscosity_temperature_range, &
      maximum_multiplicity,maximum_parcels,khrt_b0,khrt_b1,khrt_c3, &
      inflow_temperature,inflow_pressure,inflow_velocity,inflow_composition_count, &
      inflow_composition_species,inflow_composition_mole_fractions, &
      hydro_reconstruction,riemann_solver,gas_temperature,gas_pressure, &
      gas_velocity,shear_velocity,mole_names,mole_amounts, &
      boiling_model,critical_temperature,normal_boiling_temperature,liquid_molecular_weight, &
      liquid_density_model,liquid_density_coefficients,liquid_density_temperature_range, &
      vapor_pressure_model,antoine_coefficients,antoine_temperature_range, &
      vapor_name,liquid_density,liquid_cp,reference_temperature,reference_pressure,latent_heat,drop_temperature, &
      drop_diameter,diameter_shape,diameter_min,diameter_max,diameter_weighting,initial_liquid_mass,parcel_count, &
      diameter_distribution,diameter_std,nozzle_radius,port_nozzle_radii, &
      nozzle_profile,swirl_angle,hollow_spread,port_nozzle_profiles,port_swirl_angles,port_hollow_spreads, &
      velocity_from_rate,discharge_coefficient,port_discharge_coefficients, &
      port_diameter_distributions,port_diameter_weightings,port_diameter_stds, &
      port_diameter_shapes,port_diameter_minima,port_diameter_maxima, &
      wall_film_evolution,wall_liquid_conductivity,minimum_film_height, &
      wall_conductivity_model,wall_conductivity_coefficients,wall_conductivity_temperature_range, &
      port_wall_conductivities,port_wall_conductivity_models, &
      port_wall_conductivity_coefficients,port_wall_conductivity_temperature_ranges, &
      wall_impact_faces,wall_contact_angle,maximum_wall_films, &
      liquid_cp_model,liquid_cp_coefficients,liquid_cp_temperature_range, &
      port_cp_models,port_cp_coefficients,port_cp_temperature_ranges, &
      port_surface_tensions,port_liquid_viscosities,port_viscosity_models, &
      port_viscosity_coefficients,port_viscosity_temperature_ranges, &
      port_density_models,port_vapor_pressure_models,port_boiling_models,port_critical_temperatures, &
      port_normal_boiling_temperatures,port_molecular_weights, &
      port_density_coefficients,port_density_temperature_ranges,port_antoine_coefficients,port_antoine_temperature_ranges, &
      port_vapor_names,port_liquid_densities,port_liquid_cps,port_reference_temperatures, &
      port_reference_pressures,port_latent_heats, &
      injection_port_count,injection_port_centers,injection_port_axes,injection_port_fractions, &
      injection_port_speeds,injection_port_diameters,injection_port_temperatures,injection_port_cone_angles, &
      independent_port_schedules,port_mass_flow_rates,port_injection_starts,port_injection_ends,port_injection_intervals, &
      injection_center,injection_axis,injection_speed,cone_angle, &
      injection_profile_points,injection_profile_times,injection_profile_rates, &
      port_profile_points,port_profile_times,port_profile_rates, &
      mass_flow_rate,injection_start,injection_end,injection_interval,final_time,maximum_dt,cfl,maximum_steps, &
      chemistry_enabled,transport_enabled,chemistry_method,rtol,atol,drag_model,drag_enabled,heat_enabled,evaporation_enabled, &
      abramzon_sirignano,max_fractional_change,particle_cfl,minimum_diameter,maximum_substeps,smagorinsky_constant,les_model, &
      les_test_filter,les_test_ratio,les_grid_ratio,les_harmonic_mean,les_ci, &
      turbulent_prandtl,turbulent_schmidt,output_prefix,checkpoint_file,restart_file,stop_after_steps

    wall_conductivity_model='constant';wall_conductivity_coefficients=0;wall_conductivity_temperature_range=0
    port_wall_conductivities=-1;port_wall_conductivity_models=''
    port_wall_conductivity_coefficients=0;port_wall_conductivity_temperature_ranges=0
    wall_film_evolution=.false.;wall_liquid_conductivity=0;minimum_film_height=1e-7_dp
    wall_impact_faces=.false.;wall_contact_angle=acos(-1.0_dp)/3;maximum_wall_films=1000000
    allocate(wall_films(0))
    liquid_cp_model='constant';liquid_cp_coefficients=0;liquid_cp_temperature_range=0
    port_cp_models='constant';port_cp_coefficients=0;port_cp_temperature_ranges=0
    port_surface_tensions=-1;port_liquid_viscosities=-1;port_viscosity_models=''
    port_viscosity_coefficients=0;port_viscosity_temperature_ranges=0
    port_density_models='constant';port_vapor_pressure_models='clausius';port_boiling_models='reject'
    port_density_coefficients=0;port_density_temperature_ranges=0
    port_antoine_coefficients=0;port_antoine_temperature_ranges=0;port_critical_temperatures=0
    port_vapor_names='';port_liquid_densities=-1;port_liquid_cps=-1
    port_reference_temperatures=-1;port_reference_pressures=-1;port_latent_heats=-1
    inflow_temperature=1200;inflow_pressure=101325;inflow_velocity=0
    inflow_composition_count=0;inflow_composition_species="";inflow_composition_mole_fractions=0
    khrt_b0=.61_dp;khrt_b1=7.0_dp;khrt_c3=1.0_dp
    maximum_multiplicity=0;maximum_parcels=1000000
    liquid_viscosity_model='constant';liquid_viscosity_coefficients=0;liquid_viscosity_temperature_range=0
    breakup_model="none"; surface_tension=0; liquid_viscosity=0
    acceleration=0
    inflow_mach=-1
    outlet_pressure=0
    boundary_faces="periodic";wall_temperature=300;wall_velocity=0
    cells=[4,4,4]; domain_length=0.004_dp
    hydro_reconstruction='characteristic_plm'; riemann_solver='hllc'
    gas_temperature=1200.0_dp; gas_pressure=101325.0_dp; gas_velocity=0.0_dp; shear_velocity=0.0_dp
    mole_names=''; mole_names(1:3)=[character(len=24)::'H2','O2','N2']
    mole_amounts=0.0_dp; mole_amounts(1:3)=[2.0_dp,1.0_dp,3.76_dp]
    boiling_model='reject';critical_temperature=0;normal_boiling_temperature=0
    port_normal_boiling_temperatures=0
    liquid_molecular_weight=0;port_molecular_weights=0
    liquid_density_model='constant';liquid_density_coefficients=0;liquid_density_temperature_range=0
    vapor_pressure_model='clausius';antoine_coefficients=0;antoine_temperature_range=0
    vapor_name='H2O'; liquid_density=1000.0_dp; liquid_cp=4184.0_dp
    reference_temperature=373.15_dp; reference_pressure=101325.0_dp; latent_heat=2.257e6_dp
    drop_temperature=330.0_dp; drop_diameter=20.0e-6_dp; initial_liquid_mass=1.0e-10_dp; parcel_count=4
    port_profile_points=0;port_profile_times=0;port_profile_rates=0
    injection_profile_points=0;injection_profile_times=0;injection_profile_rates=0
    mixture_density_models='constant';mixture_cp_models='constant'
    mixture_vapor_pressure_models='clausius';mixture_boiling_models='reject'
    mixture_density_coefficients=0;mixture_density_temperature_ranges=0
    mixture_cp_coefficients=0;mixture_cp_temperature_ranges=0
    mixture_antoine_coefficients=0;mixture_antoine_temperature_ranges=0;mixture_critical_temperatures=0
    mixture_normal_boiling_temperatures=0
    mixture_molecular_weights=0
    port_mixture_fractions=0
    mixture_count=0;mixture_vapor_names='';mixture_fractions=0;mixture_densities=0;mixture_cps=0
    mixture_reference_temperatures=0;mixture_reference_pressures=101325;mixture_latent_heats=0
    khrt_children=.false.;khrt_reference_diameter_cubed=0
    velocity_from_rate=.false.;discharge_coefficient=1;port_discharge_coefficients=-1
    nozzle_profile='legacy';swirl_angle=0;hollow_spread=0
    port_nozzle_profiles='';port_swirl_angles=-999;port_hollow_spreads=-1
    nozzle_radius=0;port_nozzle_radii=-1
    port_diameter_distributions='';port_diameter_weightings=''
    port_diameter_stds=-1;port_diameter_shapes=-1;port_diameter_minima=-1;port_diameter_maxima=-1
    diameter_distribution='weibull';diameter_std=0
    diameter_weighting='mass'
    diameter_shape=0.0_dp; diameter_min=1.0e-12_dp; diameter_max=1.0_dp
    injection_port_speeds=-1;injection_port_diameters=-1;injection_port_temperatures=-1;injection_port_cone_angles=-1
    independent_port_schedules=.false.;port_events=0
    port_mass_flow_rates=0;port_injection_starts=0;port_injection_ends=0;port_injection_intervals=0
    injection_port_count=0;injection_port_centers=0;injection_port_axes=0;injection_port_fractions=0
    injection_center=0.0015_dp; injection_axis=[0.0_dp,0.0_dp,1.0_dp]; injection_speed=10.0_dp; cone_angle=0.1_dp
    mass_flow_rate=0.0_dp; injection_start=0.0_dp; injection_end=1.0e-5_dp; injection_interval=1.0e-6_dp
    final_time=1.0e-5_dp; maximum_dt=1.0e-7_dp; cfl=0.3_dp; maximum_steps=100000; stop_after_steps=0
    chemistry_enabled=.true.; transport_enabled=.true.; chemistry_method='implicit'; rtol=1.0e-7_dp; atol=1.0e-12_dp
    drag_model="schiller-naumann"
    drag_enabled=.true.; heat_enabled=.true.; evaporation_enabled=.true.; abramzon_sirignano=.true.
    max_fractional_change=0.05_dp; particle_cfl=0.25_dp; minimum_diameter=1.0e-9_dp; maximum_substeps=100000
    les_ci=0
    les_model='smagorinsky';les_test_filter='box';les_test_ratio=2;les_grid_ratio=1;les_harmonic_mean=.false.
    smagorinsky_constant=0.0_dp; turbulent_prandtl=0.7_dp; turbulent_schmidt=0.7_dp
    output_prefix='spray'; checkpoint_file=''; restart_file=''
    open(newunit=unit,file=trim(input_path),status='old',action='read',iostat=ios)
    if (ios/=0) call fail('Cannot open spray input')
    read(unit,nml=spray,iostat=ios,iomsg=message)
    if (ios/=0) call fail('Invalid spray namelist: '//trim(message))
    profile_file='';profile_files=''
    rewind(unit)
    read(unit,nml=spray_rate_files,iostat=ios)
    if(ios>0) call fail('Invalid spray_rate_files namelist')
    close(unit)
    if(injection_port_count<0.or.injection_port_count>16) call fail('Invalid injection port count')
    if(any(profile_files(injection_port_count+1:)/='')) call fail('Unused rate files')
    if(len_trim(profile_file)>0) then
      if(injection_profile_points/=0.or.any(injection_profile_times/=0).or. &
        any(injection_profile_rates/=0)) call fail('Rate file conflicts with inline profile')
      call read_spray_injection_profile(trim(profile_file),injection_profile_times, &
        injection_profile_rates,injection_profile_points,ok)
      if(.not.ok) call fail('Invalid spray rate file')
    end if
    do port=1,injection_port_count
      if(len_trim(profile_files(port))==0) cycle
      if(port_profile_points(port)/=0.or.any(port_profile_times(:,port)/=0).or. &
        any(port_profile_rates(:,port)/=0)) call fail('Rate file conflicts with inline profile')
      call read_spray_injection_profile(trim(profile_files(port)),port_profile_times(:,port), &
        port_profile_rates(:,port),port_profile_points(port),ok)
      if(.not.ok) call fail('Invalid spray rate file')
    end do
    if (diameter_distribution/='weibull'.and.diameter_distribution/='normal'.and. &
        diameter_distribution/='lognormal'.and.diameter_distribution/='chi-squared') &
      call fail('Invalid diameter distribution')
    if (.not.ieee_is_finite(diameter_std)) call fail('Nonfinite diameter standard deviation')
    if (diameter_std<0.or.diameter_std>1.or.(diameter_std>0.and.diameter_std<1e-15_dp)) &
      call fail('Invalid diameter standard deviation')
    if ((diameter_distribution=='weibull'.or.diameter_distribution=='chi-squared').and.diameter_std/=0) &
      call fail('Unused diameter standard deviation')
    if (diameter_distribution/='weibull'.and.diameter_shape/=0) call fail('Unused Weibull shape')
    if (diameter_weighting/='mass'.and.diameter_weighting/='number') call fail('Invalid diameter weighting')
    if (.not.valid_spray_boundary_faces(boundary_faces)) call fail('Invalid spray boundary faces')
    if (.not.all(ieee_is_finite(acceleration))) call fail('Invalid acceleration')
    if (maxval(abs(acceleration))>1e6_dp) call fail('Excessive acceleration')
    if (.not.all(ieee_is_finite(inflow_mach))) call fail('Nonfinite inlet Mach number')
    if (any(inflow_mach>100).or.any(inflow_mach>=0.and.boundary_faces/='inflow')) &
      call fail('Active inlet Mach number requires an inflow face and value at most 100')
    if (.not.all(ieee_is_finite(outlet_pressure))) call fail('Invalid outlet pressure')
    if (any(outlet_pressure<0).or.any(outlet_pressure>1e15_dp)) call fail('Invalid outlet pressure')
    if (any(outlet_pressure>0.and.boundary_faces/='outflow')) call fail('Outlet pressure requires outflow face')
    physical_boundaries=any(boundary_faces/='periodic')
    impact_enabled=any(wall_impact_faces)
    wall_options%evolve=wall_film_evolution;wall_options%liquid_conductivity=wall_liquid_conductivity
    wall_options%conductivity_model=wall_conductivity_model
    wall_options%conductivity_coefficients=wall_conductivity_coefficients
    wall_options%conductivity_temperature_range=wall_conductivity_temperature_range
    wall_options%port_conductivities=port_wall_conductivities
    wall_options%port_conductivity_models=port_wall_conductivity_models
    wall_options%port_conductivity_coefficients=port_wall_conductivity_coefficients
    wall_options%port_conductivity_temperature_ranges=port_wall_conductivity_temperature_ranges
    wall_options%minimum_height=minimum_film_height
    wall_options%enabled=wall_impact_faces;wall_options%temperature=wall_temperature
    wall_options%contact_angle=wall_contact_angle;wall_options%maximum_films=maximum_wall_films
    if (.not.valid_wall_options(wall_options,boundary_faces)) call fail('Invalid wall impact settings')
    if (impact_enabled.and.any(wall_velocity/=0)) call fail('Wall impact requires stationary walls')
    if (.not.valid_wall_velocity_3d(boundary_faces,wall_velocity)) call fail('Invalid tangential wall velocity')
    if (maxval(abs(wall_velocity))>1e7_dp) call fail('Excessive wall velocity')
    if (any(boundary_faces=='isothermal_wall').and..not.transport_enabled) &
      call fail('Isothermal walls require molecular transport')
    if (any(abs(wall_velocity)>0).and..not.transport_enabled) call fail('Moving walls require molecular transport')
    do face=1,6
      if (boundary_faces(face)/='isothermal_wall') cycle
      if (.not.ieee_is_finite(wall_temperature(face))) call fail('Invalid wall temperature')
      if (wall_temperature(face)<=0) call fail('Invalid wall temperature')
    end do
    nvar=reactive_nvar(size(species))
    if (nvar<7) call fail('Spray requires 2--128 gas species')
    if (.not. compatible_transport_database(species,transport)) call fail('Incompatible gas transport database')
    if (trim(hydro_reconstruction)/='pcm' .and. trim(hydro_reconstruction)/='characteristic_plm') &
      call fail('Invalid hydro reconstruction')
    if (trim(hydro_reconstruction)=='characteristic_plm' .and. any(cells<3)) &
      call fail('Characteristic PLM requires at least three cells per direction')
    if (trim(riemann_solver)/='hllc' .and. trim(riemann_solver)/='rusanov' .and. trim(riemann_solver)/='pelec') &
      call fail('Invalid Riemann solver')
    if (any(cells<2) .or. any(cells>64)) call fail('Spray cell counts must be in [2,64]')
    if (.not. all(ieee_is_finite([domain_length,gas_temperature,gas_pressure,gas_velocity,shear_velocity, &
        mole_amounts,liquid_density,liquid_cp,reference_temperature,reference_pressure,latent_heat,drop_temperature, &
        drop_diameter,diameter_shape,diameter_min,diameter_max,initial_liquid_mass, &
        injection_center,injection_axis,injection_speed,cone_angle,mass_flow_rate, &
        injection_start,injection_end,injection_interval,final_time,maximum_dt,cfl,rtol,atol,smagorinsky_constant, &
        turbulent_prandtl,turbulent_schmidt,max_fractional_change,particle_cfl,minimum_diameter]))) &
      call fail('Nonfinite spray input')
    if (any(domain_length<1.0e-6_dp) .or. any(domain_length>1.0e3_dp)) call fail('Invalid domain length')
    if (gas_temperature<1.0_dp .or. gas_temperature>1.0e5_dp .or. gas_pressure<1.0_dp .or. &
        gas_pressure>1.0e10_dp .or. maxval(abs(gas_velocity))>1.0e5_dp .or. &
        abs(shear_velocity)>1.0e5_dp) call fail('Invalid gas state')
    if (maximum_dt<1.0e-15_dp .or. maximum_dt>1.0_dp .or. final_time<=0.0_dp .or. &
        final_time>1000.0_dp .or. cfl<=0.0_dp .or. cfl>0.5_dp .or. maximum_steps<1 .or. &
        maximum_steps>10000000 .or. stop_after_steps<0) call fail('Invalid time controls')
    if (rtol<=0.0_dp .or. rtol>0.01_dp .or. atol<=0.0_dp .or. atol>0.001_dp) call fail('Invalid chemistry tolerances')
    if (les_model/='smagorinsky'.and.les_model/='wale'.and.les_model/='vreman'.and.les_model/='dynamic') &
      call fail('Invalid LES model')
    if (.not.ieee_is_finite(les_ci)) call fail('Invalid LES CI')
    if (les_ci<0.or.les_ci>1) call fail('Invalid LES CI')
    if (les_ci/=0.and.les_model=='dynamic') call fail('Fixed LES CI requires a static model')
    if (les_harmonic_mean.and.les_model/='dynamic') call fail('Harmonic LES interpolation requires dynamic model')
    if (les_model=='dynamic') then
      call les_filter_weights(les_test_filter,les_test_ratio,les_weights,ok)
      if (.not.ok.or.les_grid_ratio<1.or.les_grid_ratio>64) call fail('Invalid dynamic LES filter')
    else
      if (les_test_filter/='box'.or.les_test_ratio/=2.or.les_grid_ratio/=1) &
        call fail('Dynamic LES filter controls require dynamic model')
    end if
    if (smagorinsky_constant<0.0_dp .or. smagorinsky_constant>0.5_dp .or. &
        min(turbulent_prandtl,turbulent_schmidt)<0.1_dp .or. &
        max(turbulent_prandtl,turbulent_schmidt)>10.0_dp) call fail('Invalid LES controls')
    if (trim(chemistry_method)/='implicit' .and. trim(chemistry_method)/='explicit' .and. &
        trim(chemistry_method)/='cvode') call fail('Unknown chemistry backend')
    if (trim(chemistry_method)=='cvode' .and. .not. spray_cvode_available()) &
      call fail('CVODE spray backend not built; enable PELEF_ENABLE_SUNDIALS')
    if (initial_liquid_mass<0.0_dp .or. initial_liquid_mass>1.0_dp .or. &
        parcel_count<1 .or. parcel_count>100000) call fail('Invalid initial injection')
    if (injection_profile_points<0.or.injection_profile_points>64.or.injection_profile_points==1) &
      call fail('Injection profile requires zero or 2..64 points')
    if (.not.all(ieee_is_finite([injection_profile_times,injection_profile_rates]))) &
      call fail('Nonfinite injection profile')
    if (any(injection_profile_times(injection_profile_points+1:)/=0).or. &
        any(injection_profile_rates(injection_profile_points+1:)/=0)) call fail('Unused injection profile entries')
    if (injection_profile_points>0) then
      if (mass_flow_rate/=0.or.independent_port_schedules) call fail('Profile requires zero global rate and shared schedule')
      if (any(injection_profile_times<0).or.any(injection_profile_times>1000).or. &
          any(injection_profile_rates<0).or.any(injection_profile_rates>1)) call fail('Invalid injection profile values')
      do i=2,injection_profile_points
        if (injection_profile_times(i)<=injection_profile_times(i-1)) call fail('Injection profile times must increase')
      end do
      if (injection_profile_times(1)>injection_start.or. &
          injection_profile_times(injection_profile_points)<injection_end) call fail('Profile does not cover injection window')
      mass_flow_rate=maxval(injection_profile_rates(:injection_profile_points))
    end if
    if (mass_flow_rate<0.0_dp .or. mass_flow_rate>1.0_dp .or. injection_start<0.0_dp .or. &
        injection_end<injection_start .or. injection_end>1000.0_dp .or. &
        injection_interval<1.0e-15_dp .or. injection_interval>1000.0_dp) call fail('Invalid injection schedule')
    if (mass_flow_rate>0.0_dp) then
      if (injection_end<=injection_start) call fail('Nonzero injection requires a nonempty time window')
      if ((injection_end-injection_start)/injection_interval>1000000.0_dp) call fail('Too many injection events')
    end if
    if (.not.all(ieee_is_finite([nozzle_radius,port_nozzle_radii]))) call fail('Nonfinite nozzle radius')
    if (nozzle_radius<0.or.nozzle_radius>1e9_dp) call fail('Invalid nozzle radius')
    if (injection_port_count==0) then
      if (.not.nozzle_inside_domain(injection_center,injection_axis,nozzle_radius)) &
        call fail('Nozzle disk outside domain or invalid axis')
      if (any(injection_center<0.0_dp) .or. any(injection_center>=domain_length)) &
        call fail('Injection center must lie in the domain')
    end if
    if (injection_port_count<0.or.injection_port_count>16) call fail('Invalid injection port count')
    if (any(port_nozzle_radii(injection_port_count+1:)/=-1)) call fail('Unused port nozzle radius')
    nozzle_model=spray_nozzle_model(nozzle_profile,swirl_angle,hollow_spread)
    if (.not.all(ieee_is_finite([swirl_angle,hollow_spread,port_swirl_angles,port_hollow_spreads]))) &
      call fail('Nonfinite nozzle model')
    if (any(port_nozzle_profiles(injection_port_count+1:)/='').or. &
        any(port_swirl_angles(injection_port_count+1:)/=-999).or. &
        any(port_hollow_spreads(injection_port_count+1:)/=-1)) call fail('Unused port nozzle model')
    port_nozzle_overrides=any(port_nozzle_profiles/='').or.any(port_swirl_angles/=-999).or.any(port_hollow_spreads/=-1)
    allocate(port_nozzle_models(injection_port_count))
    do port=1,injection_port_count
      if (port_nozzle_profiles(port)=='') port_nozzle_profiles(port)=nozzle_profile
      if (port_swirl_angles(port)==-999) port_swirl_angles(port)=swirl_angle
      if (port_hollow_spreads(port)==-1) port_hollow_spreads(port)=hollow_spread
      port_nozzle_models(port)=spray_nozzle_model(port_nozzle_profiles(port),port_swirl_angles(port),port_hollow_spreads(port))
    end do
    if (.not.all(ieee_is_finite([discharge_coefficient,port_discharge_coefficients]))) &
      call fail('Nonfinite discharge coefficient')
    if (any(port_discharge_coefficients(injection_port_count+1:)/=-1)) call fail('Unused port discharge coefficient')
    if (.not.velocity_from_rate.and.(discharge_coefficient/=1.or.any(port_discharge_coefficients/=-1))) &
      call fail('Discharge coefficients require rate-derived velocity')
    if (velocity_from_rate) then
      if (initial_liquid_mass/=0) call fail('Rate-derived velocity requires scheduled-only injection')
      if (discharge_coefficient<1e-6_dp.or.discharge_coefficient>1) call fail('Invalid discharge coefficient')
      if (injection_port_count==0.and.nozzle_radius<1e-12_dp) call fail('Rate-derived velocity needs finite nozzle radius')
    end if
    do port=1,injection_port_count
      if (port_discharge_coefficients(port)==-1) port_discharge_coefficients(port)=discharge_coefficient
      if (port_discharge_coefficients(port)<1e-6_dp.or.port_discharge_coefficients(port)>1) &
        call fail('Invalid port discharge coefficient')
    end do
    port_radius_overrides=any(port_nozzle_radii/=-1)
    do port=1,injection_port_count
      if (port_nozzle_radii(port)==-1) port_nozzle_radii(port)=nozzle_radius
      if (velocity_from_rate.and.port_nozzle_radii(port)<1e-12_dp) call fail('Rate-derived port needs finite radius')
      if (.not.nozzle_inside_domain(injection_port_centers(:,port),injection_port_axes(:,port),port_nozzle_radii(port))) &
        call fail('Port nozzle disk outside domain or invalid radius/axis')
    end do
    if (.not.all(ieee_is_finite(injection_port_centers)).or..not.all(ieee_is_finite(injection_port_axes)).or. &
        .not.all(ieee_is_finite(injection_port_fractions))) call fail('Nonfinite injection ports')
    if (any(injection_port_centers(:,injection_port_count+1:)/=0).or. &
        any(injection_port_axes(:,injection_port_count+1:)/=0).or. &
        any(injection_port_fractions(injection_port_count+1:)/=0)) call fail('Unused injection port values')
    if (.not.all(ieee_is_finite([injection_port_speeds,injection_port_diameters, &
        injection_port_temperatures,injection_port_cone_angles]))) call fail('Nonfinite injection port state')
    if (any(injection_port_speeds(injection_port_count+1:)/=-1).or. &
        any(injection_port_diameters(injection_port_count+1:)/=-1).or. &
        any(injection_port_temperatures(injection_port_count+1:)/=-1).or. &
        any(injection_port_cone_angles(injection_port_count+1:)/=-1)) call fail('Unused injection port state')
    if (.not.all(ieee_is_finite([port_diameter_stds,port_diameter_shapes, &
        port_diameter_minima,port_diameter_maxima]))) call fail('Nonfinite port diameter distribution')
    if (any(port_diameter_distributions(injection_port_count+1:)/='').or. &
        any(port_diameter_weightings(injection_port_count+1:)/='').or. &
        any(port_diameter_stds(injection_port_count+1:)/=-1).or. &
        any(port_diameter_shapes(injection_port_count+1:)/=-1).or. &
        any(port_diameter_minima(injection_port_count+1:)/=-1).or. &
        any(port_diameter_maxima(injection_port_count+1:)/=-1)) call fail('Unused port diameter distribution')
    port_size_overrides=any(port_diameter_distributions/='').or.any(port_diameter_weightings/='').or. &
      any(port_diameter_stds/=-1).or.any(port_diameter_shapes/=-1).or. &
      any(port_diameter_minima/=-1).or.any(port_diameter_maxima/=-1)
    allocate(port_sizes(injection_port_count))
    do port=1,injection_port_count
      if (port_diameter_distributions(port)=='') port_diameter_distributions(port)=diameter_distribution
      if (port_diameter_weightings(port)=='') port_diameter_weightings(port)=diameter_weighting
      if (port_diameter_stds(port)==-1) port_diameter_stds(port)=diameter_std
      if (port_diameter_shapes(port)==-1) port_diameter_shapes(port)=diameter_shape
      if (port_diameter_minima(port)==-1) port_diameter_minima(port)=diameter_min
      if (port_diameter_maxima(port)==-1) port_diameter_maxima(port)=diameter_max
      port_sizes(port)=spray_size_options(port_diameter_distributions(port),port_diameter_weightings(port), &
        port_diameter_stds(port),port_diameter_shapes(port),port_diameter_minima(port),port_diameter_maxima(port))
    end do
    port_state_overrides=any(injection_port_speeds/=-1).or.any(injection_port_diameters/=-1).or. &
      any(injection_port_temperatures/=-1).or.any(injection_port_cone_angles/=-1)
    do port=1,injection_port_count
      if (injection_port_speeds(port)==-1) injection_port_speeds(port)=injection_speed
      if (injection_port_diameters(port)==-1) injection_port_diameters(port)=drop_diameter
      if (injection_port_temperatures(port)==-1) injection_port_temperatures(port)=drop_temperature
      if (injection_port_cone_angles(port)==-1) injection_port_cone_angles(port)=cone_angle
      if (any(injection_port_centers(:,port)<0).or.any(injection_port_centers(:,port)>=domain_length)) &
        call fail('Injection port outside domain')
    end do
    if (.not.all(ieee_is_finite([port_mass_flow_rates,port_injection_starts,port_injection_ends, &
        port_injection_intervals]))) call fail('Nonfinite port schedule')
    if (any(port_profile_points<0).or.any(port_profile_points>64).or.any(port_profile_points==1)) &
      call fail('Port profile requires zero or 2..64 points')
    if (.not.all(ieee_is_finite(port_profile_times)).or..not.all(ieee_is_finite(port_profile_rates))) &
      call fail('Nonfinite port profile')
    if (any(port_profile_points(injection_port_count+1:)/=0)) call fail('Unused port profile')
    if (any(port_profile_points>0).and..not.independent_port_schedules) &
      call fail('Port profiles require independent schedules')
    do port=1,16
      if (any(port_profile_times(port_profile_points(port)+1:,port)/=0).or. &
          any(port_profile_rates(port_profile_points(port)+1:,port)/=0)) call fail('Unused port profile entries')
      if (port_profile_points(port)==0) cycle
      if (port_mass_flow_rates(port)/=0) call fail('Port profile conflicts with constant rate')
      if (any(port_profile_times(:,port)<0).or.any(port_profile_times(:,port)>1000).or. &
          any(port_profile_rates(:,port)<0).or.any(port_profile_rates(:,port)>1)) call fail('Invalid port profile values')
      do i=2,port_profile_points(port)
        if (port_profile_times(i,port)<=port_profile_times(i-1,port)) call fail('Port profile times must increase')
      end do
      if (port_profile_times(1,port)>port_injection_starts(port).or. &
          port_profile_times(port_profile_points(port),port)<port_injection_ends(port)) &
        call fail('Port profile does not cover injection window')
      port_mass_flow_rates(port)=maxval(port_profile_rates(:port_profile_points(port),port))
    end do
    if (independent_port_schedules) then
      if (injection_port_count==0.or.mass_flow_rate/=0) call fail('Independent schedules require ports and zero global flow')
      if (any(port_mass_flow_rates(injection_port_count+1:)/=0).or. &
          any(port_injection_starts(injection_port_count+1:)/=0).or. &
          any(port_injection_ends(injection_port_count+1:)/=0).or. &
          any(port_injection_intervals(injection_port_count+1:)/=0)) call fail('Unused port schedule values')
      total_events=0
      do port=1,injection_port_count
        if (port_mass_flow_rates(port)<0.or.port_mass_flow_rates(port)>1.or.port_injection_starts(port)<0.or. &
            port_injection_ends(port)<port_injection_starts(port).or.port_injection_ends(port)>1000.or. &
            port_injection_intervals(port)<0.or.port_injection_intervals(port)>1000) call fail('Invalid port schedule')
        if (port_mass_flow_rates(port)==0) cycle
        if (port_injection_intervals(port)<1e-15_dp.or.port_injection_ends(port)<=port_injection_starts(port)) &
          call fail('Active port needs a positive interval and time window')
        amount=(port_injection_ends(port)-port_injection_starts(port))/port_injection_intervals(port)
        if (amount>1000000) call fail('Too many independent injection events')
        total_events=total_events+ceiling(amount)
      end do
      if (total_events>1000000) call fail('Too many independent injection events')
    else
      if (any(port_mass_flow_rates/=0).or.any(port_injection_starts/=0).or.any(port_injection_ends/=0).or. &
          any(port_injection_intervals/=0)) call fail('Port schedules require independent mode')
    end if
    parcels_per_event=parcel_count*max(1,injection_port_count)
    if (mixture_count<0.or.mixture_count>min(16,size(species))) call fail('Invalid mixture component count')
    if (any(mixture_vapor_names(mixture_count+1:)/='').or.any(mixture_fractions(mixture_count+1:)/=0).or. &
        any(mixture_densities(mixture_count+1:)/=0).or.any(mixture_cps(mixture_count+1:)/=0).or. &
        any(mixture_reference_temperatures(mixture_count+1:)/=0).or. &
        any(mixture_reference_pressures(mixture_count+1:)/=101325).or.any(mixture_latent_heats(mixture_count+1:)/=0)) &
      call fail('Unused mixture component entries')
    if (any(mixture_density_models(mixture_count+1:)/='constant')) call fail('Unused mixture property entries')
    if (any(mixture_cp_models(mixture_count+1:)/='constant')) call fail('Unused mixture property entries')
    if (any(mixture_vapor_pressure_models(mixture_count+1:)/='clausius')) call fail('Unused mixture property entries')
    if (any(mixture_boiling_models(mixture_count+1:)/='reject')) call fail('Unused mixture property entries')
    if (any(mixture_critical_temperatures(mixture_count+1:)/=0)) call fail('Unused mixture property entries')
    if (any(mixture_normal_boiling_temperatures(mixture_count+1:)/=0)) call fail('Unused mixture property entries')
    if (any(mixture_molecular_weights(mixture_count+1:)/=0)) call fail('Unused mixture property entries')
    if (any(mixture_density_coefficients(:,mixture_count+1:)/=0)) call fail('Unused mixture property entries')
    if (any(mixture_density_temperature_ranges(:,mixture_count+1:)/=0)) call fail('Unused mixture property entries')
    if (any(mixture_cp_coefficients(:,mixture_count+1:)/=0)) call fail('Unused mixture property entries')
    if (any(mixture_cp_temperature_ranges(:,mixture_count+1:)/=0)) call fail('Unused mixture property entries')
    if (any(mixture_antoine_coefficients(:,mixture_count+1:)/=0)) call fail('Unused mixture property entries')
    if (any(mixture_antoine_temperature_ranges(:,mixture_count+1:)/=0)) call fail('Unused mixture property entries')
    if (.not.all(ieee_is_finite(port_mixture_fractions))) call fail('Nonfinite port mixture fractions')
    if (any(port_mixture_fractions<0).or.any(port_mixture_fractions>1)) call fail('Invalid port mixture fractions')
    if (any(port_mixture_fractions(mixture_count+1:,:)/=0).or. &
        any(port_mixture_fractions(:,injection_port_count+1:)/=0)) call fail('Unused port mixture fractions')
    do port=1,injection_port_count
      if (all(port_mixture_fractions(:,port)==0)) cycle
      if (abs(sum(port_mixture_fractions(:,port))-1.0_dp)>1e-12_dp) call fail('Port mixture fractions must sum to one')
    end do
    if (mixture_count>0) then
      if (liquid_cp_model/='constant'.or.liquid_density_model/='constant'.or.boiling_model/='reject'.or. &
          vapor_pressure_model/='clausius') call fail('Use component property models in mixture mode')
      allocate(mixture_components(mixture_count))
      do i=1,mixture_count
        call initialize_spray_liquid(species,trim(mixture_vapor_names(i)),mixture_densities(i),mixture_cps(i), &
          mixture_reference_temperatures(i),mixture_reference_pressures(i),mixture_latent_heats(i),mixture_components(i),ok, &
          mixture_vapor_pressure_models(i),mixture_antoine_coefficients(:,i),mixture_antoine_temperature_ranges(:,i))
        if (.not.ok) call fail('Invalid mixture component properties')
        mixture_components(i)%density_model=mixture_density_models(i)
        mixture_components(i)%density_coefficients=mixture_density_coefficients(:,i)
        mixture_components(i)%density_temperature_range=mixture_density_temperature_ranges(:,i)
        mixture_components(i)%heat_capacity_model=mixture_cp_models(i)
        mixture_components(i)%heat_capacity_coefficients=mixture_cp_coefficients(:,i)
        mixture_components(i)%heat_capacity_temperature_range=mixture_cp_temperature_ranges(:,i)
        mixture_components(i)%boiling_model=mixture_boiling_models(i)
        mixture_components(i)%critical_temperature=mixture_critical_temperatures(i)
        mixture_components(i)%normal_boiling_temperature=mixture_normal_boiling_temperatures(i)
        mixture_components(i)%molecular_weight=mixture_molecular_weights(i)
        if (.not.valid_spray_liquid(mixture_components(i))) call fail('Invalid mixture component property model')
      end do
      liquid=mixture_injection_liquid(drop_temperature)
    else
    call initialize_spray_liquid(species,trim(vapor_name),liquid_density,liquid_cp,reference_temperature, &
      reference_pressure,latent_heat,liquid,ok,vapor_pressure_model,antoine_coefficients,antoine_temperature_range)
    if (.not. ok) call fail('Invalid liquid/vapor thermodynamics')
    end if
    liquid%heat_capacity_model=liquid_cp_model;liquid%heat_capacity_coefficients=liquid_cp_coefficients
    liquid%heat_capacity_temperature_range=liquid_cp_temperature_range
    liquid%boiling_model=boiling_model;liquid%critical_temperature=critical_temperature
    liquid%normal_boiling_temperature=normal_boiling_temperature
    liquid%molecular_weight=liquid_molecular_weight
    liquid%density_model=liquid_density_model;liquid%density_coefficients=liquid_density_coefficients
    liquid%density_temperature_range=liquid_density_temperature_range
    if (.not.valid_spray_liquid(liquid)) call fail('Invalid liquid property model')
    if (.not.all(ieee_is_finite([port_liquid_densities,port_liquid_cps,port_reference_temperatures, &
        port_reference_pressures,port_latent_heats]))) call fail('Nonfinite port material')
    if (any(port_vapor_names(injection_port_count+1:)/='').or. &
        any(port_liquid_densities(injection_port_count+1:)/=-1).or. &
        any(port_liquid_cps(injection_port_count+1:)/=-1).or. &
        any(port_reference_temperatures(injection_port_count+1:)/=-1).or. &
        any(port_reference_pressures(injection_port_count+1:)/=-1).or. &
        any(port_latent_heats(injection_port_count+1:)/=-1)) call fail('Unused port material')
    if (.not.all(ieee_is_finite(port_density_coefficients)).or. &
        .not.all(ieee_is_finite(port_density_temperature_ranges)).or. &
        .not.all(ieee_is_finite(port_antoine_coefficients)).or. &
        .not.all(ieee_is_finite(port_antoine_temperature_ranges)).or. &
        .not.all(ieee_is_finite(port_critical_temperatures))) call fail('Nonfinite port property model')
    if (.not.all(ieee_is_finite(port_normal_boiling_temperatures))) call fail('Nonfinite port boiling point')
    if (any(port_normal_boiling_temperatures(injection_port_count+1:)/=0)) call fail('Unused port boiling point')
    if (.not.all(ieee_is_finite(port_molecular_weights))) call fail('Nonfinite port molecular weight')
    if (any(port_molecular_weights(injection_port_count+1:)/=0)) call fail('Unused port molecular weight')
    if (any(port_density_models(injection_port_count+1:)/='constant').or. &
        any(port_vapor_pressure_models(injection_port_count+1:)/='clausius').or. &
        any(port_boiling_models(injection_port_count+1:)/='reject').or. &
        any(port_density_coefficients(:,injection_port_count+1:)/=0).or. &
        any(port_density_temperature_ranges(:,injection_port_count+1:)/=0).or. &
        any(port_antoine_coefficients(:,injection_port_count+1:)/=0).or. &
        any(port_antoine_temperature_ranges(:,injection_port_count+1:)/=0).or. &
        any(port_critical_temperatures(injection_port_count+1:)/=0)) call fail('Unused port property model')
    port_property_models=any(port_density_models/='constant').or.any(port_vapor_pressure_models/='clausius').or. &
      any(port_boiling_models/='reject').or.any(port_density_coefficients/=0).or. &
      any(port_density_temperature_ranges/=0).or.any(port_antoine_coefficients/=0).or. &
      any(port_antoine_temperature_ranges/=0).or.any(port_critical_temperatures/=0).or. &
      any(port_normal_boiling_temperatures/=0).or.any(port_molecular_weights/=0)
    if (.not.all(ieee_is_finite(port_surface_tensions)).or. &
        .not.all(ieee_is_finite(port_liquid_viscosities)).or. &
        .not.all(ieee_is_finite(port_viscosity_coefficients)).or. &
        .not.all(ieee_is_finite(port_viscosity_temperature_ranges))) call fail('Nonfinite port breakup property')
    if (any(port_surface_tensions(injection_port_count+1:)/=-1).or. &
        any(port_liquid_viscosities(injection_port_count+1:)/=-1).or. &
        any(port_viscosity_models(injection_port_count+1:)/='').or. &
        any(port_viscosity_coefficients(:,injection_port_count+1:)/=0).or. &
        any(port_viscosity_temperature_ranges(:,injection_port_count+1:)/=0)) call fail('Unused port breakup property')
    port_breakup_properties=any(port_surface_tensions/=-1).or.any(port_liquid_viscosities/=-1).or. &
      any(port_viscosity_models/='').or.any(port_viscosity_coefficients/=0).or.any(port_viscosity_temperature_ranges/=0)
    if (port_breakup_properties.and.breakup_model=='none'.and..not.impact_enabled) &
      call fail('Port breakup properties require breakup or wall impact')
    if (.not.all(ieee_is_finite(port_cp_coefficients)).or. &
        .not.all(ieee_is_finite(port_cp_temperature_ranges))) call fail('Nonfinite port heat capacity')
    if (any(port_cp_models(injection_port_count+1:)/='constant').or. &
        any(port_cp_coefficients(:,injection_port_count+1:)/=0).or. &
        any(port_cp_temperature_ranges(:,injection_port_count+1:)/=0)) call fail('Unused port heat capacity')
    port_cp_overrides=any(port_cp_models/='constant').or.any(port_cp_coefficients/=0).or.any(port_cp_temperature_ranges/=0)
    if (any(port_wall_conductivities(injection_port_count+1:)/=-1).or. &
        any(port_wall_conductivity_models(injection_port_count+1:)/='').or. &
        any(port_wall_conductivity_coefficients(:,injection_port_count+1:)/=0).or. &
        any(port_wall_conductivity_temperature_ranges(:,injection_port_count+1:)/=0)) &
      call fail('Unused port wall conductivity')
    port_wall_properties=any(port_wall_conductivities/=-1).or.any(port_wall_conductivity_models/='')
    if ((port_wall_properties.or.wall_conductivity_model/='constant').and..not.wall_film_evolution) &
      call fail('Wall conductivity models require film evolution')
    multi_material=port_wall_properties.or.port_cp_overrides.or.port_breakup_properties.or.port_property_models.or. &
      any(port_vapor_names/='').or. &
      any(port_liquid_densities/=-1).or.any(port_liquid_cps/=-1).or. &
      any(port_reference_temperatures/=-1).or.any(port_reference_pressures/=-1).or.any(port_latent_heats/=-1)
    if (multi_material) then
      allocate(materials(injection_port_count))
      do port=1,injection_port_count
        if (port_vapor_names(port)=='') port_vapor_names(port)=vapor_name
        if (port_liquid_densities(port)==-1) port_liquid_densities(port)=liquid_density
        if (port_liquid_cps(port)==-1) port_liquid_cps(port)=liquid_cp
        if (port_reference_temperatures(port)==-1) port_reference_temperatures(port)=reference_temperature
        if (port_reference_pressures(port)==-1) port_reference_pressures(port)=reference_pressure
        if (port_latent_heats(port)==-1) port_latent_heats(port)=latent_heat
        call initialize_spray_liquid(species,trim(port_vapor_names(port)),port_liquid_densities(port), &
          port_liquid_cps(port),port_reference_temperatures(port),port_reference_pressures(port), &
          port_latent_heats(port),materials(port),ok,port_vapor_pressure_models(port), &
          port_antoine_coefficients(:,port),port_antoine_temperature_ranges(:,port))
        if (.not.ok) call fail('Invalid port liquid/vapor thermodynamics')
        materials(port)%heat_capacity_model=port_cp_models(port)
        materials(port)%heat_capacity_coefficients=port_cp_coefficients(:,port)
        materials(port)%heat_capacity_temperature_range=port_cp_temperature_ranges(:,port)
        materials(port)%density_model=port_density_models(port)
        materials(port)%density_coefficients=port_density_coefficients(:,port)
        materials(port)%density_temperature_range=port_density_temperature_ranges(:,port)
        materials(port)%boiling_model=port_boiling_models(port)
        materials(port)%critical_temperature=port_critical_temperatures(port)
        materials(port)%normal_boiling_temperature=port_normal_boiling_temperatures(port)
        materials(port)%molecular_weight=port_molecular_weights(port)
        if (.not.valid_spray_liquid(materials(port))) call fail('Invalid port property model')
      end do
    end if
    options%maximum_multiplicity=maximum_multiplicity;options%maximum_parcels=maximum_parcels
    options%khrt_children=khrt_children;options%khrt_reference_diameter_cubed=khrt_reference_diameter_cubed
    split_enabled=maximum_multiplicity>0.or.khrt_children
    lineage_enabled=split_enabled.or.impact_enabled
    options%khrt_b0=khrt_b0;options%khrt_b1=khrt_b1;options%khrt_c3=khrt_c3
    khrt_enabled=trim(breakup_model)=="khrt"
    options%breakup_model=breakup_model; options%surface_tension=surface_tension
    options%liquid_viscosity_model=liquid_viscosity_model
    options%liquid_viscosity_coefficients=liquid_viscosity_coefficients
    options%liquid_viscosity_temperature_range=liquid_viscosity_temperature_range
    options%liquid_viscosity=liquid_viscosity; breakup_enabled=trim(breakup_model)/="none"
    options%drag_model=drag_model
    options%drag=drag_enabled; options%heat=heat_enabled; options%evaporation=evaporation_enabled
    options%abramzon_sirignano=abramzon_sirignano; options%max_fractional_change=max_fractional_change
    options%particle_cfl=particle_cfl; options%minimum_diameter=minimum_diameter; options%maximum_substeps=maximum_substeps
    if (.not. valid_spray_options(options)) call fail('Invalid parcel integration controls')
    if (port_breakup_properties) then
      allocate(breakup_options(injection_port_count))
      do port=1,injection_port_count
        breakup_options(port)=options
        if (port_surface_tensions(port)/=-1) breakup_options(port)%surface_tension=port_surface_tensions(port)
        if (port_liquid_viscosities(port)/=-1) breakup_options(port)%liquid_viscosity=port_liquid_viscosities(port)
        if (port_viscosity_models(port)/='') then
          breakup_options(port)%liquid_viscosity_model=port_viscosity_models(port)
          breakup_options(port)%liquid_viscosity_coefficients=port_viscosity_coefficients(:,port)
          breakup_options(port)%liquid_viscosity_temperature_range=port_viscosity_temperature_ranges(:,port)
        else if (any(port_viscosity_coefficients(:,port)/=0).or.any(port_viscosity_temperature_ranges(:,port)/=0)) then
          call fail('Port viscosity coefficients require explicit model')
        end if
        if (.not.valid_spray_options(breakup_options(port))) call fail('Invalid port breakup property')
        if (liquid_viscosity_at_temperature(breakup_options(port),injection_port_temperatures(port))<=0) &
          call fail('Port injection outside viscosity range')
      end do
    end if
    if (impact_enabled) then
      block
        type(spray_options) :: checked_properties
        do port=1,max(1,injection_port_count)
          checked_properties=options
          if (allocated(breakup_options)) checked_properties=breakup_options(port)
          if (checked_properties%breakup_model=='none') checked_properties%breakup_model='etab'
          if (.not.valid_spray_options(checked_properties)) call fail('Invalid wall liquid surface tension or viscosity')
        end do
      end block
    end if
    if (initial_liquid_mass>0.and.parcels_per_event>maximum_parcels) call fail('Initial parcel count exceeds storage limit')
    if (mixture_count>0) then
      if (multi_material) call fail('Mixture mode cannot use pure-port properties')
      if (injection_port_count>0) then
        allocate(materials(injection_port_count));multi_material=.true.
        do port=1,injection_port_count
          materials(port)=mixture_injection_liquid(injection_port_temperatures(port),port)
        end do
      end if
    end if
    ! Validate injection geometry even for a zero-mass/gas-only run.
    call inject_ports(1.0e-12_dp,1,1_int64,new_parcels,ok)
    if (.not. ok) call fail('Invalid parcel injection geometry or material state')
    deallocate(new_parcels)
    allocate(x(size(species)),y(size(species)),primitive(reactive_nprim(size(species))))
    x=0.0_dp
    if (any(mole_amounts<0.0_dp) .or. maxval(mole_amounts)>1.0e10_dp) call fail('Invalid gas mole amounts')
    do i=1,size(mole_names)
      if (len_trim(mole_names(i))==0) then
        if (mole_amounts(i)/=0.0_dp) call fail('Mole amount lacks a species name')
        cycle
      end if
      s=0
      do j=1,size(species)
        if (trim(species(j)%name)==trim(mole_names(i))) s=j
      end do
      if (s==0) call fail('Unknown gas species: '//trim(mole_names(i)))
      if (count(mole_names==mole_names(i))/=1) call fail('Duplicate input species name')
      x(s)=mole_amounts(i)
    end do
    if (sum(x)<=0.0_dp) call fail('Gas mole amounts sum to zero')
    x=x/sum(x)
    call mass_fractions_from_mole_fractions(species,x,y,ok)
    if (.not. ok) call fail('Invalid gas composition')
    gas_constant=mixture_specific_gas_constant(species,y,ok)
    if (.not. ok) call fail('Invalid gas molecular weight')
    spacing=domain_length/real(cells,dp); volume=product(spacing)
    allocate(state(nvar,cells(1),cells(2),cells(3)),temperature(cells(1),cells(2),cells(3)))
    allocate(reference_totals(nvar),injected_totals(nvar),totals(nvar))
    if (any(acceleration/=0)) then
      allocate(source_totals(nvar),source_step(nvar));source_totals=0;source_step=0
    end if
    allocate(boundary_totals(nvar),boundary_step(nvar)); boundary_totals=0; boundary_step=0
    primitive(1)=gas_pressure/(gas_constant*gas_temperature); primitive(2:4)=gas_velocity
    primitive(5)=gas_pressure; primitive(6:)=y
    do k=1,cells(3)
      do j=1,cells(2)
        primitive(2)=gas_velocity(1)+shear_velocity*sin(2.0_dp*acos(-1.0_dp)*(real(j,dp)-0.5_dp)/real(cells(2),dp))
        do i=1,cells(1)
          call reactive_primitive_to_conserved(species,primitive,state(:,i,j,k),local_temperature,sound,ok)
          if (.not. ok) call fail('Initial gas state is outside the mechanism thermodynamic range')
          temperature(i,j,k)=local_temperature
        end do
      end do
    end do
    allocate(inflow_state(nvar,6),inlet_x(size(species)),inlet_y(size(species)),inlet_cube(nvar,2,2,2))
    inflow_state=0;inflow_t=inflow_temperature;inlet_dt=huge(1.0_dp);inlet_x=x
    mach_pressure=merge(inflow_pressure,0._dp,inflow_mach>=0)
    if (any(inflow_mach>=0)) then
      allocate(inlet_ghost(nvar,cells(1)+4,cells(2)+4,cells(3)+4))
      allocate(inlet_ghost_t(cells(1)+4,cells(2)+4,cells(3)+4))
    end if
    if (any(boundary_faces=='inflow')) then
      if (inflow_composition_count/=0) then
        call resolve_selected_composition('Spray inflow',inflow_composition_count,inflow_composition_species, &
          inflow_composition_mole_fractions,species,inlet_x,ok,message)
        if (.not.ok) call fail(trim(message))
      end if
      call mass_fractions_from_mole_fractions(species,inlet_x,inlet_y,ok)
      if (.not.ok) call fail('Invalid inlet composition')
      gas_constant=mixture_specific_gas_constant(species,inlet_y,ok)
      if (.not.ok) call fail('Invalid inlet gas constant')
      do face=1,6
        if (boundary_faces(face)/='inflow') cycle
        if (.not.all(ieee_is_finite([inflow_temperature(face),inflow_pressure(face),inflow_velocity(:,face)]))) &
          call fail('Nonfinite inlet state')
        if (inflow_mach(face)>=0) then
          if (inflow_temperature(face)<=0.or.inflow_pressure(face)<0.or.inflow_pressure(face)>1e15_dp) &
            call fail('Invalid Mach inlet temperature or pressure')
          if (any(inflow_velocity(:,face)/=0)) call fail('Mach inlet requires zero inflow_velocity')
          cycle
        end if
        if (min(inflow_temperature(face),inflow_pressure(face))<=0.or.maxval(abs(inflow_velocity(:,face)))>1e7_dp) &
          call fail('Invalid inlet state')
        primitive=0
        primitive(1)=inflow_pressure(face)/(gas_constant*inflow_temperature(face))
        primitive(2:4)=inflow_velocity(:,face);primitive(5)=inflow_pressure(face);primitive(6:)=inlet_y
        call reactive_primitive_to_conserved(species,primitive,inflow_state(:,face),inflow_t(face),sound,ok)
        if (.not.ok) call fail('Invalid inlet thermodynamics')
        inlet_dt=min(inlet_dt,cfl/sum((abs(primitive(2:4))+sound)/spacing))
        if (transport_enabled) then
          do k=1,2
            do j=1,2
              do i=1,2
                inlet_cube(:,i,j,k)=inflow_state(:,face)
              end do
            end do
          end do
          inlet_temperature=inflow_t(face)
          call reactive_transport_timestep_3d(species,transport,inlet_cube,inlet_temperature,2,2,2, &
            spacing(1),spacing(2),spacing(3),.2_dp,.true.,.true.,.true.,face_dt,face_diffusivity,ok)
          if (.not.ok) call fail('Invalid inlet transport timestep')
          inlet_dt=min(inlet_dt,face_dt)
        end if
      end do
    end if
    allocate(parcels(0)); time=0.0_dp; steps=0; event=0; injected_totals=0.0_dp
    if (initial_liquid_mass>0.0_dp) then
      call inject_ports(initial_liquid_mass,parcel_count,1_int64,parcels,ok)
      if (.not. ok) call fail('Initial parcel injection failed')
    end if
    if (mixture_count>0) then
      allocate(component_mass(mixture_count,size(parcels)),film_component_mass(mixture_count,0))
      do i=1,size(parcels)
        component_mass(:,i)=parcels(i)%mass*injection_composition(parcels(i)%material_id)
      end do
    end if
    call current_totals(reference_totals)
    expected_initial=reference_totals
    ! The checkpoint contract includes every physical, injection and step-control
    ! input. Paths, requested stop time/count and output scheduling are not physics.
    write(context,'(a,1x,a,1x,*(i0,1x))',iostat=ios) &
      'spray-v2-ambient-skin/'//pelef_version,trim(mechanism_id),cells,parcel_count,maximum_substeps
    if (ios/=0) call fail('Cannot encode checkpoint context')
    do i=1,size(species)
      context=trim(context)//' '//trim(species(i)%name)
    end do
    call append_reals([domain_length,gas_temperature,gas_pressure,gas_velocity,shear_velocity,x,liquid_density,liquid_cp, &
      reference_temperature,reference_pressure,latent_heat,drop_temperature,drop_diameter, &
      diameter_shape,diameter_min,diameter_max,initial_liquid_mass, &
      injection_center,injection_axis,injection_speed,cone_angle,mass_flow_rate,injection_start,injection_end, &
      injection_interval,maximum_dt,cfl,rtol,atol,smagorinsky_constant,turbulent_prandtl,turbulent_schmidt, &
      max_fractional_change,particle_cfl,minimum_diameter])
    context=trim(context)//' '//trim(vapor_name)//' '//trim(chemistry_method)// &
      ' '//trim(hydro_reconstruction)//' '//trim(riemann_solver)
    block
      character(len=32) :: flags
      write(flags,'(7l2)') chemistry_enabled,transport_enabled,drag_enabled,heat_enabled,evaporation_enabled, &
        abramzon_sirignano,spray_cvode_available()
      context=trim(context)//' '//trim(flags)
    end block
    if (les_model/='smagorinsky') then
      if (len_trim(context)+14+len_trim(les_model)>len(context)) call fail('LES model context too long')
      context=trim(context)//' les-model-v1 '//trim(les_model)
    end if
    if (les_model=='dynamic') then
      if (len_trim(context)+23+len_trim(les_test_filter)>len(context)) call fail('Dynamic LES context too long')
      context=trim(context)//' dynamic-les-filter-v1 '//trim(les_test_filter)
      call append_reals([real(les_test_ratio,dp),real(les_grid_ratio,dp)])
    end if
    if (any(inflow_mach>=0)) then
      context=trim(context)//' mach-inlet-v1'
      call append_reals(inflow_mach)
    end if
    if (any(outlet_pressure>0)) then
      context=trim(context)//' pressure-outlet-v1'
      call append_reals(outlet_pressure)
    end if
    if (les_ci/=0) then
      context=trim(context)//' les-ci-v1'
      call append_reals([les_ci])
    end if
    if (les_harmonic_mean) then
      if (len_trim(context)+len(' les-harmonic-v1')>len(context)) call fail('Harmonic LES context too long')
      context=trim(context)//' les-harmonic-v1'
    end if
    if (injection_port_count>0) then
      context=trim(context)//' injection-ports-v1'
      call append_reals([real(injection_port_count,dp), &
        reshape(injection_port_centers(:,:injection_port_count),[3*injection_port_count]), &
        reshape(injection_port_axes(:,:injection_port_count),[3*injection_port_count]), &
        injection_port_fractions(:injection_port_count)])
    end if
    if (multi_material) then
      context=trim(context)//' port-materials-v1'
      do port=1,injection_port_count
        context=trim(context)//' '//trim(materials(port)%vapor_name)
        call append_reals([materials(port)%density,materials(port)%heat_capacity, &
          materials(port)%reference_temperature,materials(port)%reference_pressure,materials(port)%latent_heat])
      end do
    end if
    if (port_breakup_properties) then
      context=trim(context)//' port-breakup-properties-v1'
      do port=1,injection_port_count
        context=trim(context)//' '//trim(breakup_options(port)%liquid_viscosity_model)
        call append_reals([breakup_options(port)%surface_tension,breakup_options(port)%liquid_viscosity, &
          breakup_options(port)%liquid_viscosity_coefficients,breakup_options(port)%liquid_viscosity_temperature_range])
      end do
    end if
    if (liquid_cp_model=='polynomial') then
      context=trim(context)//' liquid-cp-polynomial-v1'
      call append_reals([liquid_cp_coefficients,liquid_cp_temperature_range])
    end if
    if (port_cp_overrides) then
      context=trim(context)//' port-cp-models-v1'
      do port=1,injection_port_count
        context=trim(context)//' '//trim(materials(port)%heat_capacity_model)
        call append_reals([materials(port)%heat_capacity_coefficients,materials(port)%heat_capacity_temperature_range])
      end do
    end if
    if (boiling_model=='pelemp-cap'.or.any(port_boiling_models=='pelemp-cap').or. &
        any(mixture_boiling_models=='pelemp-cap')) context=trim(context)//' pelemp-impact-boiling-v1'
    if (liquid_molecular_weight/=0.or.any(port_molecular_weights/=0).or.any(mixture_molecular_weights/=0)) then
      context=trim(context)//' liquid-molecular-weights-v1'
      call append_reals([liquid_molecular_weight,port_molecular_weights(:injection_port_count), &
        mixture_molecular_weights(:mixture_count)])
    end if
    if (port_property_models) then
      context=trim(context)//' port-property-models-v1'
      do port=1,injection_port_count
        context=trim(context)//' '//trim(materials(port)%density_model)//' '// &
          trim(materials(port)%vapor_pressure_model)//' '//trim(materials(port)%boiling_model)
        call append_reals([materials(port)%density_coefficients,materials(port)%density_temperature_range, &
          materials(port)%antoine_coefficients,materials(port)%antoine_temperature_range, &
          materials(port)%critical_temperature])
        if (materials(port)%boiling_model=='pelemp-cap') &
          call append_reals([materials(port)%normal_boiling_temperature])
      end do
    end if
    if (port_state_overrides) then
      context=trim(context)//' injection-port-states-v1'
      call append_reals([injection_port_speeds(:injection_port_count),injection_port_diameters(:injection_port_count), &
        injection_port_temperatures(:injection_port_count),injection_port_cone_angles(:injection_port_count)])
    end if
    if (independent_port_schedules) then
      context=trim(context)//' independent-port-schedules-v1'
      call append_reals([port_mass_flow_rates(:injection_port_count),port_injection_starts(:injection_port_count), &
        port_injection_ends(:injection_port_count),port_injection_intervals(:injection_port_count)])
    end if
    if (trim(boiling_model)=='pelemp-cap') then
      context=trim(context)//' boiling-pelemp-cap-v1'
      call append_reals([critical_temperature,normal_boiling_temperature])
    end if
    if (trim(boiling_model)=='surface-cap') then
      context=trim(context)//' boiling-surface-cap-v1'
      call append_reals([critical_temperature])
    end if
    if (trim(liquid_density_model)=='polynomial') then
      context=trim(context)//' liquid-density-polynomial-v1'
      call append_reals([liquid_density_coefficients,liquid_density_temperature_range])
    end if
    if (trim(vapor_pressure_model)=='antoine') then
      context=trim(context)//' antoine-v1'
      call append_reals([antoine_coefficients,antoine_temperature_range])
    end if
    if (allocated(source_totals)) then
      context=trim(context)//' uniform-acceleration-v1'
      call append_reals(acceleration)
    end if
    if (breakup_enabled) then
      if (khrt_enabled) then
        context=trim(context)//' khrt-representative-v1'
        call append_reals([khrt_b0,khrt_b1,khrt_c3])
      else
      context=trim(context)//' etab-size-deformation-v1/'//trim(breakup_model)
      end if
      call append_reals([surface_tension,liquid_viscosity])
      if (trim(liquid_viscosity_model)=='inverse-temperature') then
        context=trim(context)//' liquid-viscosity-inverse-temperature-v1'
        call append_reals([liquid_viscosity_coefficients,liquid_viscosity_temperature_range])
      end if
    end if
    if (split_enabled) then
      context=trim(context)//' computational-split-lineage-v1'
      call append_reals([maximum_multiplicity,real(maximum_parcels,dp)])
    end if
    if (khrt_children) then
      context=trim(context)//' khrt-shed-children-v1'
      call append_reals([khrt_reference_diameter_cubed])
    end if
    if (mixture_count>0) then
      if (any(port_mixture_fractions/=0)) then
        context=trim(context)//' port-mixture-composition-v1'
        do port=1,injection_port_count
          call append_reals(injection_composition(port))
        end do
      end if
      context=trim(context)//' '//mixture_evolution_id
      call append_reals([real(mixture_count,dp)])
      do i=1,mixture_count
        context=trim(context)//' '//trim(mixture_vapor_names(i))
        call append_reals([mixture_fractions(i),mixture_densities(i),mixture_cps(i), &
          mixture_reference_temperatures(i),mixture_reference_pressures(i),mixture_latent_heats(i)])
        if (mixture_density_models(i)/='constant'.or.mixture_cp_models(i)/='constant'.or. &
            mixture_vapor_pressure_models(i)/='clausius'.or.mixture_boiling_models(i)/='reject') then
          context=trim(context)//' component-properties-v1 '//trim(mixture_density_models(i))//' '// &
            trim(mixture_cp_models(i))//' '//trim(mixture_vapor_pressure_models(i))//' '//trim(mixture_boiling_models(i))
          call append_reals([mixture_density_coefficients(:,i),mixture_density_temperature_ranges(:,i), &
            mixture_cp_coefficients(:,i),mixture_cp_temperature_ranges(:,i),mixture_antoine_coefficients(:,i), &
            mixture_antoine_temperature_ranges(:,i),mixture_critical_temperatures(i)])
          if (mixture_boiling_models(i)=='pelemp-cap') call append_reals([mixture_normal_boiling_temperatures(i)])
        end if
      end do
    end if
    if (velocity_from_rate) then
      context=trim(context)//' rate-derived-velocity-v1'
      call append_reals([discharge_coefficient,port_discharge_coefficients(:injection_port_count)])
    end if
    if (any(port_profile_points>0)) then
      context=trim(context)//' port-injection-profiles-v1'
      do port=1,injection_port_count
        call append_reals([real(port_profile_points(port),dp), &
          port_profile_times(:port_profile_points(port),port),port_profile_rates(:port_profile_points(port),port)])
      end do
    end if
    if (injection_profile_points>0) then
      context=trim(context)//' piecewise-linear-injection-v1'
      call append_reals([real(injection_profile_points,dp),injection_profile_times(:injection_profile_points), &
        injection_profile_rates(:injection_profile_points)])
    end if
    if (nozzle_profile/='legacy'.or.swirl_angle/=0.or.hollow_spread/=0.or.port_nozzle_overrides) then
      context=trim(context)//' nozzle-profile-v1 '//trim(nozzle_profile)
      call append_reals([swirl_angle,hollow_spread])
      do port=1,injection_port_count
        context=trim(context)//' '//trim(port_nozzle_models(port)%profile)
        call append_reals([port_nozzle_models(port)%swirl,port_nozzle_models(port)%spread])
      end do
    end if
    if (nozzle_radius>0.or.port_radius_overrides) then
      context=trim(context)//' finite-nozzle-radius-v1'
      call append_reals([nozzle_radius,port_nozzle_radii(:injection_port_count)])
    end if
    if (port_size_overrides) then
      context=trim(context)//' port-diameter-distributions-v1'
      do port=1,injection_port_count
        context=trim(context)//' '//trim(port_sizes(port)%distribution)//' '//trim(port_sizes(port)%weighting)
        call append_reals([port_sizes(port)%std,port_sizes(port)%shape,port_sizes(port)%minimum,port_sizes(port)%maximum])
      end do
    end if
    if (diameter_distribution/='weibull') then
      context=trim(context)//' diameter-distribution-v1 '//trim(diameter_distribution)
      call append_reals([diameter_std])
    end if
    if (diameter_weighting/='mass') context=trim(context)//' diameter-weighting-v1 '//trim(diameter_weighting)
    if (drag_model/='schiller-naumann') context=trim(context)//' drag-model-v1 '//trim(drag_model)
    if (wall_film_evolution) then
      context=trim(context)//' wall-film-evolution-v2-half-skin'
      call append_reals([wall_liquid_conductivity,minimum_film_height])
      if (wall_conductivity_model/='constant'.or.port_wall_properties) then
        context=trim(context)//' wall-conductivity-v1 '//trim(wall_conductivity_model)
        call append_reals([wall_conductivity_coefficients,wall_conductivity_temperature_range])
        do port=1,injection_port_count
          context=trim(context)//' ['//trim(port_wall_conductivity_models(port))//']'
          call append_reals([port_wall_conductivities(port),port_wall_conductivity_coefficients(:,port), &
            port_wall_conductivity_temperature_ranges(:,port)])
        end do
      end if
    end if
    if (impact_enabled) then
      context=trim(context)//' wall-impact-film-v1'
      if (mixture_count>0) context=trim(context)//' '//mixture_film_evolution_id
      call append_reals([merge(1.0_dp,0.0_dp,wall_impact_faces),wall_temperature,wall_contact_angle, &
        real(maximum_wall_films,dp),real(maximum_parcels,dp),surface_tension,liquid_viscosity])
      context=trim(context)//' '//trim(liquid_viscosity_model)
      call append_reals([liquid_viscosity_coefficients,liquid_viscosity_temperature_range])
    end if
    if (physical_boundaries) then
      context=trim(context)//' specular-outflow-v1'
      do i=1,6
        context=trim(context)//' '//trim(boundary_faces(i))
      end do
    end if
    if (any(boundary_faces=='isothermal_wall').or.any(abs(wall_velocity)>0)) then
      context=trim(context)//' spray-wall-thermal-motion-v1'
      call append_reals(reshape(wall_velocity,[18]))
      do face=1,6
        if (boundary_faces(face)=='isothermal_wall') call append_reals([wall_temperature(face)])
      end do
    end if
    if (any(boundary_faces=='inflow')) then
      context=trim(context)//' gas-reservoir-v1'
      call append_reals(inlet_x)
      do face=1,6
        if (boundary_faces(face)=='inflow') &
          call append_reals([inflow_temperature(face),inflow_pressure(face),inflow_velocity(:,face)])
      end do
    end if
    if (len_trim(restart_file)>0) then
      if (mixture_count>0) then
        if (impact_enabled) then
          call load_checkpoint(component_mass,film_component_mass)
        else
          call load_checkpoint(component_mass)
        end if
      else
        call load_checkpoint()
      end if
      if (.not. ok) call fail(trim(message))
      if (mixture_count>0) then
        if (any(shape(component_mass)/=[mixture_count,size(parcels)])) call fail('Restart mixture shape mismatch')
      end if
      do i=1,size(parcels)
        material=select_spray_liquid(parcels(i),liquid,materials)
        if (.not.valid_spray_liquid(material)) call fail('Restart material ID outside configured pool')
      end do
      if (size(wall_films)>maximum_wall_films) call fail('Restart exceeds film storage limit')
      do i=1,size(wall_films)
        film_identity=spray_parcel();film_identity%material_id=wall_films(i)%material_id
        material=select_spray_liquid(film_identity,liquid,materials)
        if (.not.valid_spray_liquid(material)) call fail('Invalid restart film material')
        if (mixture_count>0) material=current_film_liquid(i)
        if (liquid_heat_capacity_at_temperature(material,wall_films(i)%temperature)<=0) &
          call fail('Invalid restart film temperature')
        if (wall_film_evolution) then
          if (wall_liquid_conductivity_at_temperature(wall_options,wall_films(i)%material_id, &
              wall_films(i)%temperature)<=0) call fail('Invalid restart film conductivity temperature')
        end if
        amount=liquid_density_at_temperature(material,wall_films(i)%temperature)* &
          acos(-1.0_dp)*.25_dp*wall_films(i)%diameter**2*wall_films(i)%height
        if (abs(amount-wall_films(i)%mass)>1e-10_dp*wall_films(i)%mass) call fail('Restart film volume mismatch')
        if (any(wall_films(i)%position<0).or.any(wall_films(i)%position>domain_length)) &
          call fail('Restart film outside domain')
        exists=.false.
        do face=1,6
          if (.not.wall_impact_faces(face)) cycle
          j=(face+1)/2
          amount=merge(0.0_dp,domain_length(j),modulo(face,2)==1)
          if (abs(wall_films(i)%position(j)-amount)>1e-12_dp*domain_length(j)) cycle
          if (abs(wall_films(i)%normal(j)-merge(1.0_dp,-1.0_dp,modulo(face,2)==1))>1e-12_dp) cycle
          exists=.true.
        end do
        if (.not.exists) call fail('Restart film not attached to an enabled wall')
      end do
      if (size(parcels)>maximum_parcels) call fail('Restart exceeds parcel storage limit')
      if (time>final_time .or. steps>maximum_steps .or. event>1000000) call fail('Restart exceeds requested limits')
      if (any(abs(reference_totals-expected_initial)>1.0e-11_dp* &
          max(abs(expected_initial),1.0e-30_dp))) call fail('Restart initial inventory mismatch')
      if (independent_port_schedules) then
        amount=0;port_events=0
        tolerance=32*epsilon(1.0_dp)*max(time,maximum_dt)
        profile_batches=0
        do completed=1,event
          call next_port_event(next_event,selected_port)
          if (selected_port==0) call fail('Restart has an impossible port event')
          if (next_event>time+tolerance) call fail('Restart port event is in the future')
          profile_dose=port_mass_flow_rates(selected_port)* &
            min(port_injection_intervals(selected_port),port_injection_ends(selected_port)-next_event)
          if (port_profile_points(selected_port)>0) profile_dose=profile_mass(next_event, &
            min(port_injection_ends(selected_port),next_event+port_injection_intervals(selected_port)),selected_port)
          amount=amount+profile_dose
          if (profile_dose>0) profile_batches=profile_batches+1
          port_events(selected_port)=port_events(selected_port)+1
        end do
        if (abs(injected_totals(1)-amount)>1e-11_dp*max(amount,1e-30_dp)) call fail('Restart port mass ledger mismatch')
        call next_port_event(next_event,selected_port)
        if (next_event<time-tolerance) call fail('Restart missed a port event')
      else
      if (mass_flow_rate==0.0_dp .and. event/=0) call fail('Restart has an impossible injection event')
      amount=mass_flow_rate*min(real(event,dp)*injection_interval,injection_end-injection_start)
      if (injection_profile_points>0) amount=profile_mass(injection_start, &
        min(injection_end,injection_start+real(event,dp)*injection_interval))
      if (abs(injected_totals(1)-amount)>1.0e-11_dp*max(amount,1.0e-30_dp)) &
        call fail('Restart injection mass ledger mismatch')
      tolerance=32.0_dp*epsilon(1.0_dp)*max(time,maximum_dt,injection_interval)
      if (mass_flow_rate>0.0_dp) then
        next_event=injection_start+real(event,dp)*injection_interval
        if (event>0) then
          if (next_event-injection_interval>=injection_end .or. next_event-injection_interval>time+tolerance) &
            call fail('Restart injection event is in the future')
        end if
        if (next_event<injection_end .and. next_event<time-tolerance) call fail('Restart missed an injection event')
      end if
      end if
      expected_batches=event
      if (independent_port_schedules.and.any(port_profile_points>0)) expected_batches=profile_batches
      if (injection_profile_points>0) then
        expected_batches=0
        do profile_event=0,event-1
          profile_start=injection_start+real(profile_event,dp)*injection_interval
          if (profile_mass(profile_start,min(injection_end,profile_start+injection_interval))>0) &
            expected_batches=expected_batches+1
        end do
      end if
      if (int(count(parcels%parent_id==0_int64),int64)/=int(expected_batches,int64)*&
          int(merge(parcel_count,parcels_per_event,independent_port_schedules),int64)+ &
          int(merge(parcels_per_event,0,initial_liquid_mass>0.0_dp),int64)) call fail('Restart parcel/event count mismatch')
      call current_totals(totals)
      if (allocated(source_totals)) totals=totals-source_totals
      if (abs(totals(1)-reference_totals(1)-injected_totals(1)-boundary_totals(1))> &
          1.0e-9_dp*max(reference_totals(1)+injected_totals(1),1.0e-30_dp)) call fail('Restart mass ledger mismatch')
      if (abs(totals(5)-reference_totals(5)-injected_totals(5)-boundary_totals(5))> &
          1.0e-9_dp*max(abs(reference_totals(5))+abs(injected_totals(5)),1.0e-20_dp)) &
        call fail('Restart energy ledger mismatch')
      if (any(parcels%id<1_int64)) call fail('Invalid restart parcel IDs')
      do i=1,size(parcels)
        if (parcels(i)%id/=int(i,int64)) call fail('Restart parcel IDs are not canonical')
        if (any(parcels(i)%position<0.0_dp) .or. any(parcels(i)%position>=domain_length)) &
          call fail('Restart parcel outside domain')
      end do
    end if
    if (len_trim(output_prefix)==0) call fail('Empty output prefix')
    paths=[character(len=1050)::trim(output_prefix)//'.history.csv', &
      trim(output_prefix)//'.gas.csv',trim(output_prefix)//'.parcels.csv']
    if (impact_enabled) paths=[character(len=1050)::paths,trim(output_prefix)//'.films.csv']
    do i=1,size(paths)
      inquire(file=trim(paths(i)),exist=exists)
      if (exists) call fail('Output exists; refusing to overwrite: '//trim(paths(i)))
    end do
    if (len_trim(checkpoint_file)>0) then
      inquire(file=trim(checkpoint_file),exist=exists)
      if (exists) call fail('Checkpoint destination already exists')
      if (any(paths==checkpoint_file)) call fail('Checkpoint path aliases an output')
    end if
    open(newunit=history_unit,file=trim(paths(1)),status='new',action='write',iostat=ios)
    if (ios/=0) call fail('Cannot create history output')
    write(history_unit,'(a)',advance='no',iostat=ios) 'time,steps,gas_mass,liquid_mass,system_mass,px,py,pz,total_energy,'// &
      'injected_mass,mass_error,px_error,py_error,pz_error,energy_error,tmin,tmax,live_parcels,smd'
    if (ios==0.and.impact_enabled) write(history_unit,'(a)',advance='no',iostat=ios) ',film_mass'
    if (ios==0.and.allocated(source_totals)) write(history_unit,'(a)',advance='no',iostat=ios) &
      ',source_px,source_py,source_pz,source_energy'
    if (ios==0) then
      if (physical_boundaries) then
        write(history_unit,'(a)',iostat=ios) ',boundary_mass,boundary_px,boundary_py,boundary_pz,boundary_energy'
      else
        write(history_unit,'(a)',iostat=ios) ''
      end if
    end if
    if (ios/=0) call fail('Cannot write history header')
    call write_history()
    do while(time<final_time)
      if (steps>=maximum_steps) call fail('Maximum step count reached before final time')
      if (stop_after_steps>0) then
        if (steps>=stop_after_steps) exit
      end if
      tolerance=16.0_dp*epsilon(1.0_dp)*max(time,maximum_dt,injection_interval)
      next_event=huge(1.0_dp)
      if (independent_port_schedules) then
        tolerance=16*epsilon(1.0_dp)*max(time,maximum_dt)
        call next_port_event(next_event,selected_port)
        do while(selected_port>0.and.next_event<=time+tolerance)
          amount=port_mass_flow_rates(selected_port)* &
            min(port_injection_intervals(selected_port),port_injection_ends(selected_port)-next_event)
          if (port_profile_points(selected_port)>0) amount=profile_mass(next_event, &
            min(port_injection_ends(selected_port),next_event+port_injection_intervals(selected_port)),selected_port)
          if (amount>0) then
          if (size(parcels)+parcel_count>maximum_parcels) call fail('Parcel storage limit reached')
          first_id=1_int64
          if (size(parcels)>0) first_id=parcels(size(parcels))%id+1_int64
          material=liquid
          if (multi_material) material=materials(selected_port)
          event_speed=injection_port_speeds(selected_port)
          if (velocity_from_rate) event_speed=rate_velocity(amount/ &
            min(port_injection_intervals(selected_port),port_injection_ends(selected_port)-next_event),selected_port)
          call inject_cone_parcels(material,amount,parcel_count,first_id,injection_port_centers(:,selected_port), &
            injection_port_axes(:,selected_port),event_speed,injection_port_cone_angles(selected_port), &
            injection_port_diameters(selected_port),injection_port_temperatures(selected_port), &
            new_parcels,ok,port_sizes(selected_port)%shape,port_sizes(selected_port)%minimum, &
            port_sizes(selected_port)%maximum,port_sizes(selected_port)%weighting, &
            port_sizes(selected_port)%distribution,port_sizes(selected_port)%std,port_nozzle_radii(selected_port), &
            port_nozzle_models(selected_port))
          if (.not.ok) call fail('Independent port injection failed')
          if (multi_material) new_parcels%material_id=selected_port
          call append_injected()
          end if
          event=event+1;port_events(selected_port)=port_events(selected_port)+1
          call next_port_event(next_event,selected_port)
        end do
      else if (mass_flow_rate>0.0_dp) then
        next_event=injection_start+real(event,dp)*injection_interval
        if (next_event>=injection_end) next_event=huge(1.0_dp)
        if (next_event<=time+tolerance) then
          amount=mass_flow_rate*min(injection_interval,injection_end-next_event)
          if (injection_profile_points>0) amount=profile_mass(next_event,min(injection_end,next_event+injection_interval))
          if (amount>0) then
          if (size(parcels)+parcels_per_event>maximum_parcels) call fail('Parcel storage limit reached')
          first_id=1_int64
          if (size(parcels)>0) first_id=parcels(size(parcels))%id+1_int64
          call inject_ports(amount,parcel_count,first_id,new_parcels,ok, &
            amount/min(injection_interval,injection_end-next_event))
          if (.not. ok) call fail('Scheduled parcel injection failed')
          call append_injected()
          end if
          event=event+1
          next_event=injection_start+real(event,dp)*injection_interval
          if (next_event>=injection_end) next_event=huge(1.0_dp)
        end if
      end if
      call compute_reactive_cfl_timestep_3d(species,state,temperature,cells(1),cells(2),cells(3), &
        spacing(1),spacing(2),spacing(3),cfl,cfl_dt,ok)
      if (.not. ok) call fail('Cannot compute gas CFL step')
      if (any(inflow_mach>=0)) then
        call fill_mach_inlet_boundary_ghosts_3d(species,state,temperature,boundary_faces,inflow_mach, &
          inlet_ghost,inlet_ghost_t,ok,inflow_state=inflow_state,inflow_temperature=inflow_t, &
          inflow_pressure=mach_pressure,viscous_wall=.true.,wall_velocity=wall_velocity)
        if (.not.ok) call fail('Cannot construct Mach inlet timestep states')
        call compute_reactive_cfl_timestep_3d(species,inlet_ghost,inlet_ghost_t,cells(1)+4,cells(2)+4,cells(3)+4, &
          spacing(1),spacing(2),spacing(3),cfl,inlet_dt,ok)
        if (.not.ok) call fail('Cannot compute Mach inlet CFL step')
        if (transport_enabled) then
          call reactive_transport_timestep_3d(species,transport,inlet_ghost,inlet_ghost_t, &
            cells(1)+4,cells(2)+4,cells(3)+4,spacing(1),spacing(2),spacing(3), &
            .2_dp,.true.,.true.,.true.,face_dt,face_diffusivity,ok)
          if (.not.ok) call fail('Cannot compute Mach inlet transport step')
          inlet_dt=min(inlet_dt,face_dt)
        end if
      end if
      dt=min(maximum_dt,cfl_dt,inlet_dt,final_time-time)
      if (any(acceleration/=0)) dt=min(dt,0.5_dp*cfl_dt,sqrt(cfl*minval(spacing)/maxval(abs(acceleration))))
      if (any(boundary_faces=='isothermal_wall')) then
        call reactive_isothermal_timestep_3d(species,transport,state,temperature,spacing, &
          boundary_faces,wall_temperature,0.1_dp,wall_dt,ok)
        if (.not.ok) call fail('Invalid isothermal wall timestep')
        dt=min(dt,wall_dt)
      end if
      if (next_event>time) dt=min(dt,next_event-time)
      ! Avoid a spurious vanishing final step from accumulated roundoff.
      if (abs(final_time-(time+dt)) <= 8.0_dp*epsilon(1.0_dp)*max(abs(final_time),abs(time))) dt=final_time-time
      do retry=1,30
        if (mixture_count>0) then
        call advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,parcels, &
          [0.0_dp,0.0_dp,0.0_dp],spacing,dt,chemistry_enabled,transport_enabled,trim(chemistry_method),rtol,atol, &
          smagorinsky_constant,turbulent_prandtl,turbulent_schmidt,ok,message,hydro_reconstruction,riemann_solver, &
          boundary_faces,boundary_step,inflow_state,inflow_t,wall_temperature,wall_velocity,acceleration,source_step, &
          mixture_components=mixture_components,component_mass=component_mass, &
          wall_options=wall_options,wall_films=wall_films,film_component_mass=film_component_mass,les_model=les_model, &
          les_test_filter=les_test_filter,les_test_ratio=les_test_ratio,les_grid_ratio=les_grid_ratio, &
          les_harmonic_mean=les_harmonic_mean,les_ci=les_ci,outlet_pressure=outlet_pressure, &
          inflow_mach=inflow_mach,inflow_pressure=mach_pressure)
        else
        call advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,parcels, &
          [0.0_dp,0.0_dp,0.0_dp],spacing,dt,chemistry_enabled,transport_enabled,trim(chemistry_method),rtol,atol, &
          smagorinsky_constant,turbulent_prandtl,turbulent_schmidt,ok,message,hydro_reconstruction,riemann_solver, &
          boundary_faces,boundary_step,inflow_state,inflow_t,wall_temperature,wall_velocity,acceleration,source_step,materials, &
          breakup_options,wall_options,wall_films,les_model=les_model, &
          les_test_filter=les_test_filter,les_test_ratio=les_test_ratio,les_grid_ratio=les_grid_ratio, &
          les_harmonic_mean=les_harmonic_mean,les_ci=les_ci,outlet_pressure=outlet_pressure, &
          inflow_mach=inflow_mach,inflow_pressure=mach_pressure)
        end if
        if (ok) exit
        dt=0.5_dp*dt
        if (time+dt<=time .or. dt<1.0e-18_dp) exit
      end do
      if (.not. ok) call fail('Coupled step failed with rollback: '//trim(message))
      if (allocated(source_totals)) source_totals=source_totals+source_step
      boundary_totals=boundary_totals+boundary_step
      time=time+dt; steps=steps+1
      call write_history()
    end do
    close(history_unit,iostat=ios)
    if (ios/=0) call fail('History close failed')
    if (len_trim(checkpoint_file)>0) then
      if (mixture_count>0) then
        if (impact_enabled) then
          call save_checkpoint(component_mass,film_component_mass)
        else
          call save_checkpoint(component_mass)
        end if
      else
        call save_checkpoint()
      end if
      if (.not. ok) call fail(trim(message))
    end if
    call write_snapshots()
    write(*,'(a,es14.6,a,i0,a,i0)') 'Spray run completed: t=',time,', steps=',steps,', parcels=',size(parcels)
  contains
    subroutine load_checkpoint(checkpoint_mass,checkpoint_film_mass)
      real(dp), allocatable, intent(inout), optional :: checkpoint_mass(:,:),checkpoint_film_mass(:,:)
      if (physical_boundaries) then
      call read_spray_checkpoint(trim(restart_file),trim(context),species,state,temperature,parcels,time,steps,event, &
        reference_totals,injected_totals,ok,message,boundary_totals,breakup_state=breakup_enabled,split_state=split_enabled, &
        khrt_state=khrt_enabled,shedding_state=khrt_children,source_totals=source_totals,material_state=multi_material, &
        film_state=impact_enabled,wall_films=wall_films,component_mass=checkpoint_mass,film_component_mass=checkpoint_film_mass)
      else
      call read_spray_checkpoint(trim(restart_file),trim(context),species,state,temperature,parcels,time,steps,event, &
        reference_totals,injected_totals,ok,message,breakup_state=breakup_enabled,split_state=split_enabled, &
        khrt_state=khrt_enabled,shedding_state=khrt_children,source_totals=source_totals,material_state=multi_material, &
        film_state=impact_enabled,wall_films=wall_films,component_mass=checkpoint_mass,film_component_mass=checkpoint_film_mass)
      end if
    end subroutine
    subroutine save_checkpoint(checkpoint_mass,checkpoint_film_mass)
      real(dp), allocatable, intent(in), optional :: checkpoint_mass(:,:),checkpoint_film_mass(:,:)
      if (physical_boundaries) then
      call write_spray_checkpoint(trim(checkpoint_file),trim(context),species,state,temperature,parcels,time,steps,event, &
        reference_totals,injected_totals,ok,message,boundary_totals,breakup_state=breakup_enabled,split_state=split_enabled, &
        khrt_state=khrt_enabled,shedding_state=khrt_children,source_totals=source_totals,material_state=multi_material, &
        film_state=impact_enabled,wall_films=wall_films,component_mass=checkpoint_mass,film_component_mass=checkpoint_film_mass)
      else
      call write_spray_checkpoint(trim(checkpoint_file),trim(context),species,state,temperature,parcels,time,steps,event, &
        reference_totals,injected_totals,ok,message,breakup_state=breakup_enabled,split_state=split_enabled, &
        khrt_state=khrt_enabled,shedding_state=khrt_children,source_totals=source_totals,material_state=multi_material, &
        film_state=impact_enabled,wall_films=wall_films,component_mass=checkpoint_mass,film_component_mass=checkpoint_film_mass)
      end if
    end subroutine
    function injection_composition(port_id) result(fractions)
      integer, intent(in) :: port_id
      real(dp) :: fractions(mixture_count)
      fractions=mixture_fractions(:mixture_count)
      if (port_id<1.or.port_id>injection_port_count) return
      if (any(port_mixture_fractions(:,port_id)/=0)) fractions=port_mixture_fractions(:mixture_count,port_id)
    end function
    function mixture_injection_liquid(t,port_id) result(fluid)
      real(dp), intent(in) :: t
      integer, optional, intent(in) :: port_id
      type(spray_liquid) :: fluid
      real(dp) :: rho,cp,e,fractions(mixture_count)
      logical :: valid
      fractions=mixture_fractions(:mixture_count)
      if (present(port_id)) fractions=injection_composition(port_id)
      call liquid_mixture_properties(mixture_components,fractions,t,rho,cp,e,valid)
      if (.not.valid) call fail('Invalid mixture composition or injection temperature')
      fluid=mixture_components(1);fluid%density=rho;fluid%heat_capacity=cp
      fluid%reference_temperature=t;fluid%reference_energy=e
      fluid%density_model='constant';fluid%density_coefficients=0;fluid%density_temperature_range=0
      fluid%heat_capacity_model='constant';fluid%heat_capacity_coefficients=0;fluid%heat_capacity_temperature_range=0
      fluid%boiling_model='reject';fluid%critical_temperature=0
      fluid%normal_boiling_temperature=0
      if (.not.valid_spray_liquid(fluid)) call fail('Invalid effective injection properties')
    end function
    function current_film_liquid(index) result(fluid)
      integer, intent(in) :: index
      type(spray_liquid) :: fluid
      real(dp) :: rho,cp,e
      logical :: valid
      if (.not.valid_wall_component_inventory(wall_films(index),film_component_mass(:,index))) &
        call fail('Invalid mixed film component inventory')
      call liquid_mixture_properties(mixture_components,film_component_mass(:,index)/wall_films(index)%mass, &
        wall_films(index)%temperature,rho,cp,e,valid)
      if (.not.valid) call fail('Invalid mixed film thermodynamics')
      fluid=mixture_components(1);fluid%density=rho;fluid%heat_capacity=cp
      fluid%reference_temperature=wall_films(index)%temperature;fluid%reference_energy=e
      fluid%density_model='constant';fluid%density_coefficients=0;fluid%density_temperature_range=0
      fluid%heat_capacity_model='constant';fluid%heat_capacity_coefficients=0;fluid%heat_capacity_temperature_range=0
      fluid%boiling_model='reject';fluid%critical_temperature=0
      fluid%normal_boiling_temperature=0
    end function
    subroutine current_totals(result)
      real(dp), intent(out) :: result(:)
      logical :: valid
      if (mixture_count>0) then
        call mixture_system_totals(state,spacing,parcels,mixture_components,component_mass,result,valid, &
          wall_films,film_component_mass)
        if (.not.valid) call fail('Invalid mixed-liquid inventory')
      else
        call spray_system_totals(state,spacing,parcels,liquid,result,materials,wall_films)
      end if
    end subroutine
    real(dp) function current_diameter(index) result(diameter)
      integer, intent(in) :: index
      real(dp) :: rho,cp,e
      logical :: valid
      if (mixture_count==0) then
        diameter=parcel_diameter(parcels(index),select_spray_liquid(parcels(index),liquid,materials));return
      end if
      diameter=0
      if (parcels(index)%mass<=0) return
      call liquid_mixture_properties(mixture_components,component_mass(:,index)/sum(component_mass(:,index)), &
        parcels(index)%temperature,rho,cp,e,valid)
      if (.not.valid) call fail('Invalid mixture snapshot properties')
      diameter=(6*parcels(index)%mass/(acos(-1.0_dp)*rho))**(1.0_dp/3)
    end function
    subroutine append_injected()
      real(dp), allocatable :: updated(:,:)
      real(dp) :: fractions(mixture_count)
      integer :: n,c,old_count
      old_count=size(parcels)
      if (mixture_count>0) then
        allocate(updated(mixture_count,old_count+size(new_parcels)))
        updated(:,:old_count)=component_mass
      end if
      do n=1,size(new_parcels)
        material=select_spray_liquid(new_parcels(n),liquid,materials)
        totals5=parcel_totals(new_parcels(n),material)
        injected_totals(1:5)=injected_totals(1:5)+totals5
        if (mixture_count>0) then
          fractions=injection_composition(new_parcels(n)%material_id)
          updated(:,old_count+n)=new_parcels(n)%mass*fractions
          do c=1,mixture_count
            i=5+mixture_components(c)%vapor_index
            injected_totals(i)=injected_totals(i)+totals5(1)*fractions(c)
          end do
        else
          i=5+material%vapor_index
          injected_totals(i)=injected_totals(i)+totals5(1)
        end if
      end do
      if (mixture_count>0) call move_alloc(updated,component_mass)
      parcels=[parcels,new_parcels];deallocate(new_parcels)
    end subroutine
    subroutine next_port_event(event_time,event_port)
      real(dp), intent(out) :: event_time
      integer, intent(out) :: event_port
      real(dp) :: candidate_time
      integer :: index
      event_time=huge(1.0_dp);event_port=0
      do index=1,injection_port_count
        if (port_mass_flow_rates(index)<=0) cycle
        candidate_time=port_injection_starts(index)+real(port_events(index),dp)*port_injection_intervals(index)
        if (candidate_time>=port_injection_ends(index)) cycle
        ! Strict comparison breaks simultaneous events by increasing port index.
        if (candidate_time<event_time) then
          event_time=candidate_time;event_port=index
        end if
      end do
    end subroutine next_port_event

    real(dp) function profile_mass(begin_time,end_time,profile_port) result(mass)
      real(dp), intent(in) :: begin_time,end_time
      integer, intent(in), optional :: profile_port
      real(dp) :: profile_times(64),profile_rates(64)
      integer :: points
      logical :: profile_ok
      points=injection_profile_points;profile_times=injection_profile_times;profile_rates=injection_profile_rates
      if (present(profile_port)) then
        points=port_profile_points(profile_port)
        profile_times=port_profile_times(:,profile_port);profile_rates=port_profile_rates(:,profile_port)
      end if
      call integrate_spray_injection_profile(profile_times(:points),profile_rates(:points), &
        begin_time,end_time,mass,profile_ok)
      if (.not.profile_ok) call fail('Invalid injection profile integration interval')
    end function

    logical function nozzle_inside_domain(center,axis,radius) result(valid)
      real(dp), intent(in) :: center(3),axis(3),radius
      real(dp) :: normal(3),extent(3)
      valid=.false.
      if (.not.all(ieee_is_finite([center,axis,radius]))) return
      if (radius<0.or.radius>1e9_dp.or.maxval(abs(axis))>1e9_dp) return
      if (norm2(axis)<=1e-12_dp) return
      normal=axis/norm2(axis)
      extent=radius*sqrt(max(0.0_dp,1-normal**2))
      valid=all(center-extent>=0).and.all(center+extent<domain_length)
    end function

    real(dp) function rate_velocity(rate,index) result(speed)
      real(dp), intent(in) :: rate
      integer, intent(in) :: index
      type(spray_liquid) :: fluid
      real(dp) :: density,radius,coefficient,t
      fluid=liquid;t=drop_temperature;radius=nozzle_radius;coefficient=discharge_coefficient
      if (index>0) then
        if (allocated(materials)) fluid=materials(index)
        t=injection_port_temperatures(index);radius=port_nozzle_radii(index)
        coefficient=port_discharge_coefficients(index)
      end if
      speed=-1;density=liquid_density_at_temperature(fluid,t)
      if (density<=0) return
      speed=rate/(density*acos(-1.0_dp)*radius**2*coefficient)
    end function

    subroutine inject_ports(mass,number,id,result,success,batch_rate)
      real(dp), intent(in), optional :: batch_rate
      real(dp) :: speed,port_speeds(16),rate
      integer :: index
      real(dp), intent(in) :: mass
      integer, intent(in) :: number
      integer(int64), intent(in) :: id
      type(spray_parcel), allocatable, intent(out) :: result(:)
      logical, intent(out) :: success
      speed=injection_speed;port_speeds=injection_port_speeds
      if (velocity_from_rate) then
        rate=mass_flow_rate
        if (present(batch_rate)) rate=batch_rate
        if (injection_port_count==0) then
          speed=rate_velocity(rate,0)
        else
          do index=1,injection_port_count
            if (independent_port_schedules) then
              port_speeds(index)=rate_velocity(port_mass_flow_rates(index),index)
            else
              port_speeds(index)=rate_velocity(rate*injection_port_fractions(index)/ &
                sum(injection_port_fractions(:injection_port_count)),index)
            end if
          end do
        end if
      end if
      if (injection_port_count==0) then
        call inject_cone_parcels(liquid,mass,number,id,injection_center,injection_axis,speed, &
          cone_angle,drop_diameter,drop_temperature,result,success,diameter_shape,diameter_min,diameter_max, &
          diameter_weighting,diameter_distribution,diameter_std,nozzle_radius,nozzle_model)
      else
        call inject_cone_ports(liquid,mass,number,id,injection_port_centers(:,:injection_port_count), &
          injection_port_axes(:,:injection_port_count),injection_port_fractions(:injection_port_count), &
          speed,cone_angle,drop_diameter,drop_temperature,result,success,diameter_shape,diameter_min,diameter_max, &
          port_speeds(:injection_port_count),injection_port_diameters(:injection_port_count), &
          injection_port_temperatures(:injection_port_count),injection_port_cone_angles(:injection_port_count), &
          materials,diameter_weighting,diameter_distribution,diameter_std,port_sizes,nozzle_radius, &
          port_nozzle_radii(:injection_port_count),nozzle_model,port_nozzle_models)
      end if
    end subroutine inject_ports
    subroutine fail(text)
      character(len=*), intent(in) :: text
      write(*,'(a)') trim(text)
      error stop 2
    end subroutine
    subroutine append_reals(values)
      real(dp), intent(in) :: values(:)
      character(len=8192) :: line
      write(line,'(*(es26.17e3,1x))',iostat=ios) values
      if (ios/=0 .or. len_trim(context)+len_trim(line)+1>len(context)) call fail('Checkpoint context too long')
      context=trim(context)//' '//trim(line)
    end subroutine
    subroutine write_history()
      real(dp) :: gas_mass,liquid_mass,error5(5),d2,d3,d,multiplicity,smd
      integer :: p,live
      call current_totals(totals)
      gas_mass=sum(state(1,:,:,:))*volume; liquid_mass=0.0_dp; d2=0.0_dp; d3=0.0_dp; live=0
      do p=1,size(parcels)
        liquid_mass=liquid_mass+parcels(p)%mass*parcels(p)%multiplicity
        if (parcels(p)%mass<=0.0_dp) cycle
        live=live+1; d=current_diameter(p)
        multiplicity=parcels(p)%multiplicity
        d2=d2+multiplicity*d*d; d3=d3+multiplicity*d*d*d
      end do
      smd=0.0_dp
      if (d2>0.0_dp) smd=d3/d2
      error5=totals(1:5)-reference_totals(1:5)-injected_totals(1:5)-boundary_totals(1:5)
      if (allocated(source_totals)) error5=error5-source_totals(1:5)
      write(history_unit,'(es26.17e3,a,i0,15(a,es26.17e3),a,i0,a,es26.17e3)',advance='no',iostat=ios) time,',',steps, &
        ',',gas_mass,',',liquid_mass,',',totals(1),',',totals(2),',',totals(3),',',totals(4),',',totals(5), &
        ',',injected_totals(1),',',error5(1),',',error5(2),',',error5(3),',',error5(4),',',error5(5), &
        ',',minval(temperature),',',maxval(temperature),',',live,',',smd
      if (ios==0.and.impact_enabled) write(history_unit,'(a,es26.17e3)',advance='no',iostat=ios) &
        ',',sum(wall_films%mass)
      if (ios==0.and.allocated(source_totals)) write(history_unit,'(4(a,es26.17e3))',advance='no',iostat=ios) &
        ',',source_totals(2),',',source_totals(3),',',source_totals(4),',',source_totals(5)
      if (ios==0) then
        if (physical_boundaries) then
          write(history_unit,'(5(a,es26.17e3))',iostat=ios) ',',boundary_totals(1),',',boundary_totals(2), &
            ',',boundary_totals(3),',',boundary_totals(4),',',boundary_totals(5)
        else
          write(history_unit,'(a)',iostat=ios) ''
        end if
      end if
      if (ios/=0) call fail('History write failed')
    end subroutine
    subroutine write_snapshots()
      integer :: p
      character(len=2048) :: header
      open(newunit=gas_unit,file=trim(paths(2)),status='new',action='write',iostat=ios)
      if (ios/=0) call fail('Cannot create gas snapshot')
      header='x,y,z,temperature,rho,rhou,rhov,rhow,rhoE'
      do p=1,size(species)
        header=trim(header)//',rhoY_'//trim(species(p)%name)
      end do
      write(gas_unit,'(a)',iostat=ios) trim(header)
      if (ios/=0) call fail('Gas header write failed')
      do k=1,cells(3)
        do j=1,cells(2)
          do i=1,cells(1)
            write(gas_unit,'(*(es26.17e3,:,a))',iostat=ios) &
              (real(i,dp)-0.5_dp)*spacing(1),',',(real(j,dp)-0.5_dp)*spacing(2),',', &
              (real(k,dp)-0.5_dp)*spacing(3),',',temperature(i,j,k),',',(state(p,i,j,k),',',p=1,nvar-1),state(nvar,i,j,k)
            if (ios/=0) call fail('Gas snapshot write failed')
          end do
        end do
      end do
      close(gas_unit,iostat=ios)
      if (ios/=0) call fail('Gas snapshot close failed')
      open(newunit=parcel_unit,file=trim(paths(3)),status='new',action='write',iostat=ios)
      if (ios/=0) call fail('Cannot create parcel snapshot')
      write(parcel_unit,'(a)',advance='no',iostat=ios) 'id,x,y,z,u,v,w,mass,temperature,multiplicity,diameter'
      do i=1,mixture_count
        if (ios==0) write(parcel_unit,'(a)',advance='no',iostat=ios) ',mass_'//trim(mixture_vapor_names(i))
      end do
      if (ios==0.and.multi_material) write(parcel_unit,'(a)',advance='no',iostat=ios) ',material_id'
      if (ios==0.and.breakup_enabled) write(parcel_unit,'(a)',advance='no',iostat=ios) ',deformation,deformation_rate'
      if (ios==0.and.khrt_enabled) write(parcel_unit,'(a)',advance='no',iostat=ios) ',rt_time'
      if (ios==0.and.khrt_children) write(parcel_unit,'(a)',advance='no',iostat=ios) &
        ',khrt_initial_count,khrt_shed_volume'
      if (ios==0.and.lineage_enabled) write(parcel_unit,'(a)',advance='no',iostat=ios) ',parent_id'
      if (ios==0) write(parcel_unit,'(a)',iostat=ios) ''
      if (ios/=0) call fail('Parcel header write failed')
      do p=1,size(parcels)
        write(parcel_unit,'(i0,10(a,es26.17e3))',advance='no',iostat=ios) parcels(p)%id,',',parcels(p)%position(1), &
          ',',parcels(p)%position(2),',',parcels(p)%position(3),',',parcels(p)%velocity(1), &
          ',',parcels(p)%velocity(2),',',parcels(p)%velocity(3),',',parcels(p)%mass,',',parcels(p)%temperature, &
          ',',parcels(p)%multiplicity,',',current_diameter(p)
        do i=1,mixture_count
          if (ios==0) write(parcel_unit,'(a,es26.17e3)',advance='no',iostat=ios) ',',component_mass(i,p)
        end do
        if (ios==0.and.multi_material) write(parcel_unit,'(a,i0)',advance='no',iostat=ios) ',',parcels(p)%material_id
        if (ios==0.and.breakup_enabled) write(parcel_unit,'(2(a,es26.17e3))',advance='no',iostat=ios) &
          ',',parcels(p)%deformation,',',parcels(p)%deformation_rate
        if (ios==0.and.khrt_enabled) write(parcel_unit,'(a,es26.17e3)',advance='no',iostat=ios) ',',parcels(p)%rt_time
        if (ios==0.and.khrt_children) write(parcel_unit,'(2(a,es26.17e3))',advance='no',iostat=ios) &
          ',',parcels(p)%khrt_initial_count,',',parcels(p)%khrt_shed_volume
        if (ios==0.and.lineage_enabled) write(parcel_unit,'(a,i0)',advance='no',iostat=ios) ',',parcels(p)%parent_id
        if (ios==0) write(parcel_unit,'(a)',iostat=ios) ''
        if (ios/=0) call fail('Parcel snapshot write failed')
      end do
      close(parcel_unit,iostat=ios)
      if (ios/=0) call fail('Parcel snapshot close failed')
      if (impact_enabled) then
        open(newunit=film_unit,file=trim(paths(4)),status='new',action='write',iostat=ios)
        if (ios/=0) call fail('Cannot create film snapshot')
        write(film_unit,'(a)',advance='no',iostat=ios) 'material_id,x,y,z,nx,ny,nz,mass,temperature,diameter,height'
        do i=1,mixture_count
          if (ios==0) write(film_unit,'(a)',advance='no',iostat=ios) ',mass_'//trim(mixture_vapor_names(i))
        end do
        if (ios==0) write(film_unit,*,iostat=ios)
        do p=1,size(wall_films)
          if (ios/=0) exit
          write(film_unit,'(i0,10(a,es26.17e3))',advance='no',iostat=ios) wall_films(p)%material_id, &
            ',',wall_films(p)%position(1),',',wall_films(p)%position(2),',',wall_films(p)%position(3), &
            ',',wall_films(p)%normal(1),',',wall_films(p)%normal(2),',',wall_films(p)%normal(3), &
            ',',wall_films(p)%mass,',',wall_films(p)%temperature,',',wall_films(p)%diameter,',',wall_films(p)%height
          do i=1,mixture_count
            if (ios==0) write(film_unit,'(a,es26.17e3)',advance='no',iostat=ios) ',',film_component_mass(i,p)
          end do
          if (ios==0) write(film_unit,*,iostat=ios)
        end do
        if (ios/=0) call fail('Film snapshot write failed')
        close(film_unit,iostat=ios)
        if (ios/=0) call fail('Film snapshot close failed')
      end if
    end subroutine
  end subroutine run_spray_application
end module spray_application_mod
