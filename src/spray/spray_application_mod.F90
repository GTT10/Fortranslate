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
  use simulation_config_reactive_3d_mod, only: valid_wall_velocity_3d
  use spray_parcel_mod
  use spray_coupling_3d_mod, only: spray_system_totals, valid_spray_boundary_faces
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
    character(len=16) :: liquid_cp_model,port_cp_models(16)
    real(dp) :: liquid_cp_coefficients(4),liquid_cp_temperature_range(2)
    real(dp) :: port_cp_coefficients(4,16),port_cp_temperature_ranges(2,16)
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
    integer :: cells(3),parcel_count,maximum_steps,stop_after_steps,maximum_substeps
    integer :: unit,ios,i,j,k,s,steps,event,retry,history_unit,gas_unit,parcel_unit,nvar
    real(dp) :: domain_length(3),gas_temperature,gas_pressure,gas_velocity(3),shear_velocity
    real(dp) :: mole_amounts(mechanism_species_capacity),liquid_density,liquid_cp
    character(len=16) :: boiling_model
    real(dp) :: critical_temperature
    character(len=16) :: liquid_density_model
    real(dp) :: liquid_density_coefficients(4),liquid_density_temperature_range(2)
    character(len=16) :: vapor_pressure_model
    real(dp) :: antoine_coefficients(4),antoine_temperature_range(2)
    real(dp) :: reference_temperature,reference_pressure,latent_heat
    real(dp) :: diameter_shape,diameter_min,diameter_max
    real(dp) :: drop_temperature,drop_diameter,initial_liquid_mass,injection_center(3),injection_axis(3)
    real(dp) :: injection_speed,cone_angle,mass_flow_rate,injection_start,injection_end,injection_interval
    real(dp) :: final_time,maximum_dt,cfl,rtol,atol,smagorinsky_constant,turbulent_prandtl,turbulent_schmidt
    real(dp) :: max_fractional_change,particle_cfl,minimum_diameter,time,dt,cfl_dt,next_event,tolerance,amount
    real(dp) :: spacing(3),sound,local_temperature,gas_constant,totals5(5),volume
    real(dp), allocatable :: state(:,:,:,:),temperature(:,:,:),x(:),y(:),primitive(:)
    real(dp), allocatable :: reference_totals(:),injected_totals(:),totals(:),expected_initial(:)
    type(spray_parcel), allocatable :: parcels(:),new_parcels(:)
    type(spray_liquid) :: liquid
    type(spray_options) :: options
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
    character(len=1050) :: paths(3)
    integer(int64) :: first_id
    namelist /spray/ cells,domain_length,boundary_faces,wall_temperature,wall_velocity,acceleration, &
      breakup_model,surface_tension,liquid_viscosity, &
      liquid_viscosity_model,liquid_viscosity_coefficients,liquid_viscosity_temperature_range, &
      maximum_multiplicity,maximum_parcels,khrt_b0,khrt_b1,khrt_c3, &
      inflow_temperature,inflow_pressure,inflow_velocity,inflow_composition_count, &
      inflow_composition_species,inflow_composition_mole_fractions, &
      hydro_reconstruction,riemann_solver,gas_temperature,gas_pressure, &
      gas_velocity,shear_velocity,mole_names,mole_amounts, &
      boiling_model,critical_temperature, &
      liquid_density_model,liquid_density_coefficients,liquid_density_temperature_range, &
      vapor_pressure_model,antoine_coefficients,antoine_temperature_range, &
      vapor_name,liquid_density,liquid_cp,reference_temperature,reference_pressure,latent_heat,drop_temperature, &
      drop_diameter,diameter_shape,diameter_min,diameter_max,initial_liquid_mass,parcel_count, &
      liquid_cp_model,liquid_cp_coefficients,liquid_cp_temperature_range, &
      port_cp_models,port_cp_coefficients,port_cp_temperature_ranges, &
      port_surface_tensions,port_liquid_viscosities,port_viscosity_models, &
      port_viscosity_coefficients,port_viscosity_temperature_ranges, &
      port_density_models,port_vapor_pressure_models,port_boiling_models,port_critical_temperatures, &
      port_density_coefficients,port_density_temperature_ranges,port_antoine_coefficients,port_antoine_temperature_ranges, &
      port_vapor_names,port_liquid_densities,port_liquid_cps,port_reference_temperatures, &
      port_reference_pressures,port_latent_heats, &
      injection_port_count,injection_port_centers,injection_port_axes,injection_port_fractions, &
      injection_port_speeds,injection_port_diameters,injection_port_temperatures,injection_port_cone_angles, &
      independent_port_schedules,port_mass_flow_rates,port_injection_starts,port_injection_ends,port_injection_intervals, &
      injection_center,injection_axis,injection_speed,cone_angle, &
      mass_flow_rate,injection_start,injection_end,injection_interval,final_time,maximum_dt,cfl,maximum_steps, &
      chemistry_enabled,transport_enabled,chemistry_method,rtol,atol,drag_enabled,heat_enabled,evaporation_enabled, &
      abramzon_sirignano,max_fractional_change,particle_cfl,minimum_diameter,maximum_substeps,smagorinsky_constant, &
      turbulent_prandtl,turbulent_schmidt,output_prefix,checkpoint_file,restart_file,stop_after_steps

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
    boundary_faces="periodic";wall_temperature=300;wall_velocity=0
    cells=[4,4,4]; domain_length=0.004_dp
    hydro_reconstruction='characteristic_plm'; riemann_solver='hllc'
    gas_temperature=1200.0_dp; gas_pressure=101325.0_dp; gas_velocity=0.0_dp; shear_velocity=0.0_dp
    mole_names=''; mole_names(1:3)=[character(len=24)::'H2','O2','N2']
    mole_amounts=0.0_dp; mole_amounts(1:3)=[2.0_dp,1.0_dp,3.76_dp]
    boiling_model='reject';critical_temperature=0
    liquid_density_model='constant';liquid_density_coefficients=0;liquid_density_temperature_range=0
    vapor_pressure_model='clausius';antoine_coefficients=0;antoine_temperature_range=0
    vapor_name='H2O'; liquid_density=1000.0_dp; liquid_cp=4184.0_dp
    reference_temperature=373.15_dp; reference_pressure=101325.0_dp; latent_heat=2.257e6_dp
    drop_temperature=330.0_dp; drop_diameter=20.0e-6_dp; initial_liquid_mass=1.0e-10_dp; parcel_count=4
    diameter_shape=0.0_dp; diameter_min=1.0e-12_dp; diameter_max=1.0_dp
    injection_port_speeds=-1;injection_port_diameters=-1;injection_port_temperatures=-1;injection_port_cone_angles=-1
    independent_port_schedules=.false.;port_events=0
    port_mass_flow_rates=0;port_injection_starts=0;port_injection_ends=0;port_injection_intervals=0
    injection_port_count=0;injection_port_centers=0;injection_port_axes=0;injection_port_fractions=0
    injection_center=0.0015_dp; injection_axis=[0.0_dp,0.0_dp,1.0_dp]; injection_speed=10.0_dp; cone_angle=0.1_dp
    mass_flow_rate=0.0_dp; injection_start=0.0_dp; injection_end=1.0e-5_dp; injection_interval=1.0e-6_dp
    final_time=1.0e-5_dp; maximum_dt=1.0e-7_dp; cfl=0.3_dp; maximum_steps=100000; stop_after_steps=0
    chemistry_enabled=.true.; transport_enabled=.true.; chemistry_method='implicit'; rtol=1.0e-7_dp; atol=1.0e-12_dp
    drag_enabled=.true.; heat_enabled=.true.; evaporation_enabled=.true.; abramzon_sirignano=.true.
    max_fractional_change=0.05_dp; particle_cfl=0.25_dp; minimum_diameter=1.0e-9_dp; maximum_substeps=100000
    smagorinsky_constant=0.0_dp; turbulent_prandtl=0.7_dp; turbulent_schmidt=0.7_dp
    output_prefix='spray'; checkpoint_file=''; restart_file=''
    open(newunit=unit,file=trim(input_path),status='old',action='read',iostat=ios)
    if (ios/=0) call fail('Cannot open spray input')
    read(unit,nml=spray,iostat=ios,iomsg=message)
    close(unit)
    if (ios/=0) call fail('Invalid spray namelist: '//trim(message))
    if (.not.valid_spray_boundary_faces(boundary_faces)) call fail('Invalid spray boundary faces')
    if (.not.all(ieee_is_finite(acceleration))) call fail('Invalid acceleration')
    if (maxval(abs(acceleration))>1e6_dp) call fail('Excessive acceleration')
    physical_boundaries=any(boundary_faces/='periodic')
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
    if (smagorinsky_constant<0.0_dp .or. smagorinsky_constant>0.5_dp .or. &
        min(turbulent_prandtl,turbulent_schmidt)<0.1_dp .or. &
        max(turbulent_prandtl,turbulent_schmidt)>10.0_dp) call fail('Invalid LES controls')
    if (trim(chemistry_method)/='implicit' .and. trim(chemistry_method)/='explicit' .and. &
        trim(chemistry_method)/='cvode') call fail('Unknown chemistry backend')
    if (trim(chemistry_method)=='cvode' .and. .not. spray_cvode_available()) &
      call fail('CVODE spray backend not built; enable PELEF_ENABLE_SUNDIALS')
    if (initial_liquid_mass<0.0_dp .or. initial_liquid_mass>1.0_dp .or. &
        parcel_count<1 .or. parcel_count>100000) call fail('Invalid initial injection')
    if (mass_flow_rate<0.0_dp .or. mass_flow_rate>1.0_dp .or. injection_start<0.0_dp .or. &
        injection_end<injection_start .or. injection_end>1000.0_dp .or. &
        injection_interval<1.0e-15_dp .or. injection_interval>1000.0_dp) call fail('Invalid injection schedule')
    if (mass_flow_rate>0.0_dp) then
      if (injection_end<=injection_start) call fail('Nonzero injection requires a nonempty time window')
      if ((injection_end-injection_start)/injection_interval>1000000.0_dp) call fail('Too many injection events')
    end if
    if (injection_port_count==0) then
      if (any(injection_center<0.0_dp) .or. any(injection_center>=domain_length)) &
        call fail('Injection center must lie in the domain')
    end if
    if (injection_port_count<0.or.injection_port_count>16) call fail('Invalid injection port count')
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
    call initialize_spray_liquid(species,trim(vapor_name),liquid_density,liquid_cp,reference_temperature, &
      reference_pressure,latent_heat,liquid,ok,vapor_pressure_model,antoine_coefficients,antoine_temperature_range)
    if (.not. ok) call fail('Invalid liquid/vapor thermodynamics')
    liquid%heat_capacity_model=liquid_cp_model;liquid%heat_capacity_coefficients=liquid_cp_coefficients
    liquid%heat_capacity_temperature_range=liquid_cp_temperature_range
    liquid%boiling_model=boiling_model;liquid%critical_temperature=critical_temperature
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
      any(port_antoine_temperature_ranges/=0).or.any(port_critical_temperatures/=0)
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
    if (port_breakup_properties.and.breakup_model=='none') call fail('Port breakup properties require breakup')
    if (.not.all(ieee_is_finite(port_cp_coefficients)).or. &
        .not.all(ieee_is_finite(port_cp_temperature_ranges))) call fail('Nonfinite port heat capacity')
    if (any(port_cp_models(injection_port_count+1:)/='constant').or. &
        any(port_cp_coefficients(:,injection_port_count+1:)/=0).or. &
        any(port_cp_temperature_ranges(:,injection_port_count+1:)/=0)) call fail('Unused port heat capacity')
    port_cp_overrides=any(port_cp_models/='constant').or.any(port_cp_coefficients/=0).or.any(port_cp_temperature_ranges/=0)
    multi_material=port_cp_overrides.or.port_breakup_properties.or.port_property_models.or. &
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
        if (.not.valid_spray_liquid(materials(port))) call fail('Invalid port property model')
      end do
    end if
    options%maximum_multiplicity=maximum_multiplicity;options%maximum_parcels=maximum_parcels
    split_enabled=maximum_multiplicity>0
    options%khrt_b0=khrt_b0;options%khrt_b1=khrt_b1;options%khrt_c3=khrt_c3
    khrt_enabled=trim(breakup_model)=="khrt"
    options%breakup_model=breakup_model; options%surface_tension=surface_tension
    options%liquid_viscosity_model=liquid_viscosity_model
    options%liquid_viscosity_coefficients=liquid_viscosity_coefficients
    options%liquid_viscosity_temperature_range=liquid_viscosity_temperature_range
    options%liquid_viscosity=liquid_viscosity; breakup_enabled=trim(breakup_model)/="none"
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
    if (initial_liquid_mass>0.and.parcels_per_event>maximum_parcels) call fail('Initial parcel count exceeds storage limit')
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
    inflow_state=0;inflow_t=0;inlet_dt=huge(1.0_dp);inlet_x=x
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
    call spray_system_totals(state,spacing,parcels,liquid,reference_totals,materials)
    expected_initial=reference_totals
    ! The checkpoint contract includes every physical, injection and step-control
    ! input. Paths, requested stop time/count and output scheduling are not physics.
    write(context,'(a,1x,a,1x,*(i0,1x))',iostat=ios) &
      'spray-v1/'//pelef_version,trim(mechanism_id),cells,parcel_count,maximum_substeps
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
    if (port_property_models) then
      context=trim(context)//' port-property-models-v1'
      do port=1,injection_port_count
        context=trim(context)//' '//trim(materials(port)%density_model)//' '// &
          trim(materials(port)%vapor_pressure_model)//' '//trim(materials(port)%boiling_model)
        call append_reals([materials(port)%density_coefficients,materials(port)%density_temperature_range, &
          materials(port)%antoine_coefficients,materials(port)%antoine_temperature_range, &
          materials(port)%critical_temperature])
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
      if (physical_boundaries) then
      call read_spray_checkpoint(trim(restart_file),trim(context),species,state,temperature,parcels,time,steps,event, &
        reference_totals,injected_totals,ok,message,boundary_totals,breakup_state=breakup_enabled,split_state=split_enabled, &
        khrt_state=khrt_enabled,source_totals=source_totals,material_state=multi_material)
      else
      call read_spray_checkpoint(trim(restart_file),trim(context),species,state,temperature,parcels,time,steps,event, &
        reference_totals,injected_totals,ok,message,breakup_state=breakup_enabled,split_state=split_enabled, &
        khrt_state=khrt_enabled,source_totals=source_totals,material_state=multi_material)
      end if
      if (.not. ok) call fail(trim(message))
      do i=1,size(parcels)
        material=select_spray_liquid(parcels(i),liquid,materials)
        if (.not.valid_spray_liquid(material)) call fail('Restart material ID outside configured pool')
      end do
      if (size(parcels)>maximum_parcels) call fail('Restart exceeds parcel storage limit')
      if (time>final_time .or. steps>maximum_steps .or. event>1000000) call fail('Restart exceeds requested limits')
      if (any(abs(reference_totals-expected_initial)>1.0e-11_dp* &
          max(abs(expected_initial),1.0e-30_dp))) call fail('Restart initial inventory mismatch')
      if (independent_port_schedules) then
        amount=0;port_events=0
        tolerance=32*epsilon(1.0_dp)*max(time,maximum_dt)
        do completed=1,event
          call next_port_event(next_event,selected_port)
          if (selected_port==0) call fail('Restart has an impossible port event')
          if (next_event>time+tolerance) call fail('Restart port event is in the future')
          amount=amount+port_mass_flow_rates(selected_port)* &
            min(port_injection_intervals(selected_port),port_injection_ends(selected_port)-next_event)
          port_events(selected_port)=port_events(selected_port)+1
        end do
        if (abs(injected_totals(1)-amount)>1e-11_dp*max(amount,1e-30_dp)) call fail('Restart port mass ledger mismatch')
        call next_port_event(next_event,selected_port)
        if (next_event<time-tolerance) call fail('Restart missed a port event')
      else
      if (mass_flow_rate==0.0_dp .and. event/=0) call fail('Restart has an impossible injection event')
      amount=mass_flow_rate*min(real(event,dp)*injection_interval,injection_end-injection_start)
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
      if (int(count(parcels%parent_id==0_int64),int64)/=int(event,int64)*&
          int(merge(parcel_count,parcels_per_event,independent_port_schedules),int64)+ &
          int(merge(parcels_per_event,0,initial_liquid_mass>0.0_dp),int64)) call fail('Restart parcel/event count mismatch')
      call spray_system_totals(state,spacing,parcels,liquid,totals,materials)
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
    do i=1,3
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
          if (size(parcels)+parcel_count>maximum_parcels) call fail('Parcel storage limit reached')
          first_id=1_int64
          if (size(parcels)>0) first_id=parcels(size(parcels))%id+1_int64
          material=liquid
          if (multi_material) material=materials(selected_port)
          call inject_cone_parcels(material,amount,parcel_count,first_id,injection_port_centers(:,selected_port), &
            injection_port_axes(:,selected_port),injection_port_speeds(selected_port),injection_port_cone_angles(selected_port), &
            injection_port_diameters(selected_port),injection_port_temperatures(selected_port), &
            new_parcels,ok,diameter_shape,diameter_min,diameter_max)
          if (.not.ok) call fail('Independent port injection failed')
          if (multi_material) new_parcels%material_id=selected_port
          do i=1,size(new_parcels)
            material=select_spray_liquid(new_parcels(i),liquid,materials)
            totals5=parcel_totals(new_parcels(i),material)
            injected_totals(1:5)=injected_totals(1:5)+totals5
            injected_totals(5+material%vapor_index)=injected_totals(5+material%vapor_index)+totals5(1)
          end do
          parcels=[parcels,new_parcels];deallocate(new_parcels)
          event=event+1;port_events(selected_port)=port_events(selected_port)+1
          call next_port_event(next_event,selected_port)
        end do
      else if (mass_flow_rate>0.0_dp) then
        next_event=injection_start+real(event,dp)*injection_interval
        if (next_event>=injection_end) next_event=huge(1.0_dp)
        if (next_event<=time+tolerance) then
          amount=mass_flow_rate*min(injection_interval,injection_end-next_event)
          if (size(parcels)+parcels_per_event>maximum_parcels) call fail('Parcel storage limit reached')
          first_id=1_int64
          if (size(parcels)>0) first_id=parcels(size(parcels))%id+1_int64
          call inject_ports(amount,parcel_count,first_id,new_parcels,ok)
          if (.not. ok) call fail('Scheduled parcel injection failed')
          do i=1,size(new_parcels)
            material=select_spray_liquid(new_parcels(i),liquid,materials)
            totals5=parcel_totals(new_parcels(i),material)
            injected_totals(1:5)=injected_totals(1:5)+totals5
            injected_totals(5+material%vapor_index)=injected_totals(5+material%vapor_index)+totals5(1)
          end do
          parcels=[parcels,new_parcels]; deallocate(new_parcels); event=event+1
          next_event=injection_start+real(event,dp)*injection_interval
          if (next_event>=injection_end) next_event=huge(1.0_dp)
        end if
      end if
      call compute_reactive_cfl_timestep_3d(species,state,temperature,cells(1),cells(2),cells(3), &
        spacing(1),spacing(2),spacing(3),cfl,cfl_dt,ok)
      if (.not. ok) call fail('Cannot compute gas CFL step')
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
        call advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,parcels, &
          [0.0_dp,0.0_dp,0.0_dp],spacing,dt,chemistry_enabled,transport_enabled,trim(chemistry_method),rtol,atol, &
          smagorinsky_constant,turbulent_prandtl,turbulent_schmidt,ok,message,hydro_reconstruction,riemann_solver, &
          boundary_faces,boundary_step,inflow_state,inflow_t,wall_temperature,wall_velocity,acceleration,source_step,materials, &
          breakup_options)
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
      if (physical_boundaries) then
      call write_spray_checkpoint(trim(checkpoint_file),trim(context),species,state,temperature,parcels,time,steps,event, &
        reference_totals,injected_totals,ok,message,boundary_totals,breakup_state=breakup_enabled,split_state=split_enabled, &
        khrt_state=khrt_enabled,source_totals=source_totals,material_state=multi_material)
      else
      call write_spray_checkpoint(trim(checkpoint_file),trim(context),species,state,temperature,parcels,time,steps,event, &
        reference_totals,injected_totals,ok,message,breakup_state=breakup_enabled,split_state=split_enabled, &
        khrt_state=khrt_enabled,source_totals=source_totals,material_state=multi_material)
      end if
      if (.not. ok) call fail(trim(message))
    end if
    call write_snapshots()
    write(*,'(a,es14.6,a,i0,a,i0)') 'Spray run completed: t=',time,', steps=',steps,', parcels=',size(parcels)
  contains
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

    subroutine inject_ports(mass,number,id,result,success)
      real(dp), intent(in) :: mass
      integer, intent(in) :: number
      integer(int64), intent(in) :: id
      type(spray_parcel), allocatable, intent(out) :: result(:)
      logical, intent(out) :: success
      if (injection_port_count==0) then
        call inject_cone_parcels(liquid,mass,number,id,injection_center,injection_axis,injection_speed, &
          cone_angle,drop_diameter,drop_temperature,result,success,diameter_shape,diameter_min,diameter_max)
      else
        call inject_cone_ports(liquid,mass,number,id,injection_port_centers(:,:injection_port_count), &
          injection_port_axes(:,:injection_port_count),injection_port_fractions(:injection_port_count), &
          injection_speed,cone_angle,drop_diameter,drop_temperature,result,success,diameter_shape,diameter_min,diameter_max, &
          injection_port_speeds(:injection_port_count),injection_port_diameters(:injection_port_count), &
          injection_port_temperatures(:injection_port_count),injection_port_cone_angles(:injection_port_count),materials)
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
      call spray_system_totals(state,spacing,parcels,liquid,totals,materials)
      gas_mass=sum(state(1,:,:,:))*volume; liquid_mass=0.0_dp; d2=0.0_dp; d3=0.0_dp; live=0
      do p=1,size(parcels)
        liquid_mass=liquid_mass+parcels(p)%mass*parcels(p)%multiplicity
        if (parcels(p)%mass<=0.0_dp) cycle
        live=live+1; d=parcel_diameter(parcels(p),select_spray_liquid(parcels(p),liquid,materials))
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
      if (ios==0.and.multi_material) write(parcel_unit,'(a)',advance='no',iostat=ios) ',material_id'
      if (ios==0) then
        if (breakup_enabled) then
          write(parcel_unit,'(a)',advance='no',iostat=ios) ',deformation,deformation_rate'
          if (ios==0.and.khrt_enabled) write(parcel_unit,'(a)',advance='no',iostat=ios) ',rt_time'
          if (ios==0) then
            if (split_enabled) then
              write(parcel_unit,'(a)',iostat=ios) ',parent_id'
            else
              write(parcel_unit,'(a)',iostat=ios) ''
            end if
          end if
        else
          write(parcel_unit,'(a)',iostat=ios) ''
        end if
      end if
      if (ios/=0) call fail('Parcel header write failed')
      do p=1,size(parcels)
        write(parcel_unit,'(i0,10(a,es26.17e3))',advance='no',iostat=ios) parcels(p)%id,',',parcels(p)%position(1), &
          ',',parcels(p)%position(2),',',parcels(p)%position(3),',',parcels(p)%velocity(1), &
          ',',parcels(p)%velocity(2),',',parcels(p)%velocity(3),',',parcels(p)%mass,',',parcels(p)%temperature, &
          ',',parcels(p)%multiplicity,',',parcel_diameter(parcels(p),select_spray_liquid(parcels(p),liquid,materials))
        if (ios==0.and.multi_material) write(parcel_unit,'(a,i0)',advance='no',iostat=ios) ',',parcels(p)%material_id
        if (ios==0) then
          if (breakup_enabled) then
            write(parcel_unit,'(2(a,es26.17e3))',advance='no',iostat=ios) ',',parcels(p)%deformation,',',parcels(p)%deformation_rate
            if (ios==0.and.khrt_enabled) &
              write(parcel_unit,'(a,es26.17e3)',advance='no',iostat=ios) ',',parcels(p)%rt_time
            if (ios==0) then
              if (split_enabled) then
                write(parcel_unit,'(a,i0)',iostat=ios) ',',parcels(p)%parent_id
              else
                write(parcel_unit,'(a)',iostat=ios) ''
              end if
            end if
          else
            write(parcel_unit,'(a)',iostat=ios) ''
          end if
        end if
        if (ios/=0) call fail('Parcel snapshot write failed')
      end do
      close(parcel_unit,iostat=ios)
      if (ios/=0) call fail('Parcel snapshot close failed')
    end subroutine
  end subroutine run_spray_application
end module spray_application_mod
