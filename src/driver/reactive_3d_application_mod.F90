module reactive_3d_application_mod
  use dual_energy_mod, only: dual_energy_policy
  use reactive_3d_checkpoint_mod, only: write_regular_gas_checkpoint,read_regular_gas_checkpoint
  use precision_mod, only: dp
  use reactive_hit_3d_mod, only: initialize_reactive_hit_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use constant_volume_reactor_mod, only: forced_reactor_integrator
  use constants_mod, only: pelef_version
  use state_indices_mod, only: irho, imx, imy, imz, iet, ncons
  use nasa7_thermo_mod, only: nasa7_species,configure_nasa7_eos_policy
  use elementary_kinetics_mod, only: elementary_reaction
  use reactive_checkpoint_context_mod, only: write_selected_checkpoint_context
  use gas_transport_mod, only: gas_transport_species,write_simple_transport_fits
  use mesh_3d_mod, only: uniform_cell_centers_3d
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive, reactive_primitive_to_conserved
  use mixture_thermo_mod, only: mixture_density, mass_fractions_from_mole_fractions
  use selected_composition_mod, only: resolve_selected_composition
  use simulation_config_reactive_3d_mod, only: reactive_3d_config, resolved_reactive_boundaries_3d, &
    fill_mach_inlet_boundary_ghosts_3d,resolved_reactive_difmag_3d
  use reactive_entropy_wave_3d_problem_mod, only: &
    initialize_reactive_problem_3d, reactive_entropy_wave_density_3d
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use reactive_mol_3d_mod, only: advance_reactive_mol_3d,convert_regular_energy_3d
  use reactive_3d_mod, only: &
    compute_reactive_cfl_timestep_3d, advance_reactive_full_3d,advance_reactive_sdc_3d, &
    reactive_integrals_3d, reactive_extrema_3d,compute_reactive_external_source_3d,reactive_sdc_source_evaluator
  use reactive_transport_3d_mod, only: reactive_transport_timestep_3d, reactive_isothermal_timestep_3d
  use reactive_csv_io_3d_mod, only: write_reactive_3d_csv
  implicit none
  private

  public :: read_constant_transport
  public :: run_reactive_3d_application,read_regular_les_controls,read_reactive_rotation_controls

  abstract interface
    subroutine application_diffusion(u,t,h,faces,inlet,inlet_t,rates,boundary_rate,diffusivity,ok)
      import dp
      real(dp), intent(in) :: u(:,:,:,:),t(:,:,:),h(3),inlet(:,:),inlet_t(:)
      character(len=*), intent(in) :: faces(6)
      real(dp), intent(out) :: rates(:,:,:,:),boundary_rate(:),diffusivity
      logical, intent(out) :: ok
    end subroutine
  end interface
contains
  subroutine read_reactive_rotation_controls(input_path,time_integrator,frame)
    character(len=*), intent(in) :: input_path,time_integrator
    type(rotating_frame_3d), allocatable, intent(out) :: frame
    logical :: enabled
    integer :: axis,unit,ios
    real(dp) :: angular_speed,axis_point(3)
    namelist /reactive_rotation/ enabled,axis,angular_speed,axis_point
    enabled=.false.;axis=2;angular_speed=0;axis_point=0
    open(newunit=unit,file=input_path,status='old',action='read',iostat=ios)
    if(ios/=0) error stop 'Cannot open regular rotation input'
    read(unit,nml=reactive_rotation,iostat=ios)
    close(unit)
    if(ios>0) error stop 'Invalid reactive_rotation namelist'
    if(axis<0.or.axis>2.or..not.all(ieee_is_finite([angular_speed,axis_point]))) &
      error stop 'Invalid regular rotating frame'
    if(.not.enabled) then
      if(axis/=2.or.angular_speed/=0.or.any(axis_point/=0)) error stop 'Rotation controls require enable flag'
      return
    endif
    if(time_integrator/='mol') error stop 'Regular rotation requires MOL'
    allocate(frame);frame%axis=axis+1;frame%angular_speed=angular_speed;frame%axis_point=axis_point
    if(.not.frame%is_valid()) error stop 'Invalid regular rotating frame'
  end subroutine read_reactive_rotation_controls

  subroutine read_constant_transport(input_path,coefficients,hit_parameters)
    character(len=*), intent(in) :: input_path
    real(dp), allocatable, intent(out) :: coefficients(:),hit_parameters(:)
    real(dp) :: viscosity,conductivity,bulk_viscosity,diffusivity,thermal_diffusion_ratio
    real(dp) :: lambda0,reynolds_lambda0,prandtl
    character(len=16) :: mode
    logical :: enabled
    integer :: unit,ios
    namelist /reactive_constant_transport/ enabled,viscosity,conductivity,mode,lambda0,reynolds_lambda0,prandtl, &
      bulk_viscosity,diffusivity,thermal_diffusion_ratio
    enabled=.false.;viscosity=0;conductivity=0;bulk_viscosity=0;diffusivity=-1;thermal_diffusion_ratio=0
    mode='specified';lambda0=.5_dp;reynolds_lambda0=100;prandtl=.71_dp
    open(newunit=unit,file=input_path,status='old',action='read',iostat=ios)
    if(ios/=0) error stop 'Cannot open constant transport input'
    read(unit,nml=reactive_constant_transport,iostat=ios)
    close(unit)
    if(ios>0) error stop 'Invalid reactive_constant_transport namelist'
    if(.not.enabled) return
    if(.not.all(ieee_is_finite([viscosity,conductivity]))) error stop 'Nonfinite constant transport coefficient'
    if(min(viscosity,conductivity)<0) error stop 'Negative constant transport coefficient'
    if(.not.all(ieee_is_finite([bulk_viscosity,diffusivity,thermal_diffusion_ratio]))) &
      error stop 'Nonfinite additional constant transport coefficient'
    if(bulk_viscosity<0.or.(diffusivity<0.and.diffusivity/=-1)) error stop 'Invalid constant bulk viscosity or rhoD'
    if(abs(thermal_diffusion_ratio)>1e-24_dp) error stop 'Pinned PeleC does not support Soret effects'
    select case(mode)
    case('specified')
      if(bulk_viscosity/=0.or.diffusivity>=0) then
        coefficients=[viscosity,conductivity,bulk_viscosity,diffusivity]
      else
        coefficients=[viscosity,conductivity]
      endif
    case('hit')
      if(viscosity/=0.or.conductivity/=0.or.bulk_viscosity/=0.or.diffusivity/=-1) &
        error stop 'HIT mode cannot also specify transport coefficients'
      if(.not.all(ieee_is_finite([lambda0,reynolds_lambda0,prandtl]))) error stop 'Nonfinite HIT transport parameter'
      if(min(lambda0,reynolds_lambda0,prandtl)<=0) error stop 'HIT transport parameters must be positive'
      hit_parameters=[lambda0,reynolds_lambda0,prandtl];coefficients=[0._dp,0._dp]
    case default
      error stop 'Unknown constant transport mode'
    end select
  end subroutine

  subroutine resolve_hit_transport(properties,parameters,coefficients)
    real(dp), intent(in) :: properties(3),parameters(3)
    real(dp), intent(out) :: coefficients(2)
    real(dp) :: value
    if(.not.all(ieee_is_finite(properties)).or.any(properties<=0)) &
      error stop 'HIT transport requires positive initialized HIT scales'
    call multiply(properties(1),properties(2),value)
    call multiply(value,parameters(1),coefficients(1))
    call divide(coefficients(1),parameters(2),value)
    coefficients(1)=value
    call multiply(value,properties(3),coefficients(2))
    call divide(coefficients(2),parameters(3),value)
    coefficients(2)=value
  contains
    subroutine multiply(a,b,result)
      real(dp), intent(in) :: a,b
      real(dp), intent(out) :: result
      if(b>1) then
        if(a>huge(1._dp)/b) error stop 'HIT transport coefficient overflow'
      end if
      result=a*b
    end subroutine
    subroutine divide(a,b,result)
      real(dp), intent(in) :: a,b
      real(dp), intent(out) :: result
      if(b<1) then
        if(a>huge(1._dp)*b) error stop 'HIT transport coefficient overflow'
      end if
      result=a/b
    end subroutine
  end subroutine

  subroutine read_linear_forcing(input_path,coefficient,momentum_offset)
    character(len=*), intent(in) :: input_path
    real(dp), intent(out) :: coefficient,momentum_offset(3)
    integer :: unit,ios
    namelist /reactive_linear_forcing/ coefficient,momentum_offset
    coefficient=0;momentum_offset=0
    open(newunit=unit,file=input_path,status='old',action='read',iostat=ios)
    if(ios/=0) error stop 'Cannot open linear forcing input'
    read(unit,nml=reactive_linear_forcing,iostat=ios)
    close(unit)
    if(ios>0) error stop 'Invalid reactive_linear_forcing namelist'
    if(.not.all(ieee_is_finite([coefficient,momentum_offset]))) error stop 'Invalid linear forcing parameters'
    if(maxval(abs([coefficient,momentum_offset]))>1e30_dp) error stop 'Oversized linear forcing parameters'
  end subroutine

  subroutine read_regular_timestep_limit(input_path,maximum_dt,initial_timestep_factor,maximum_timestep_growth, &
      initial_dt,fixed_dt)
    character(len=*), intent(in) :: input_path
    real(dp), intent(out) :: maximum_dt,initial_timestep_factor,maximum_timestep_growth,initial_dt,fixed_dt
    integer :: unit,ios
    namelist /reactive_timestep/ maximum_dt,initial_timestep_factor,maximum_timestep_growth,initial_dt,fixed_dt
    maximum_dt=huge(1._dp);initial_timestep_factor=1;maximum_timestep_growth=1.1_dp;initial_dt=0;fixed_dt=0
    open(newunit=unit,file=input_path,status='old',action='read',iostat=ios)
    if(ios/=0) error stop 'Cannot open timestep input'
    read(unit,nml=reactive_timestep,iostat=ios)
    close(unit)
    if(ios>0) error stop 'Invalid reactive_timestep namelist'
    if(.not.ieee_is_finite(maximum_dt)) error stop 'Invalid maximum timestep'
    if(.not.all(ieee_is_finite([initial_timestep_factor,maximum_timestep_growth]))) &
      error stop 'Invalid automatic timestep controls'
    if(initial_timestep_factor<=0.or.initial_timestep_factor>1.or.maximum_timestep_growth<1) &
      error stop 'Invalid automatic timestep controls'
    if(maximum_dt<=0) error stop 'Invalid maximum timestep'
    if(.not.all(ieee_is_finite([initial_dt,fixed_dt]))) error stop 'Nonfinite prescribed timestep'
    if(fixed_dt>maximum_dt) error stop 'Fixed timestep exceeds maximum timestep'
  end subroutine

  subroutine read_regular_les_controls(input_path,enabled,cs,prt,sct,ci,model,test_model,test_ratio,grid_ratio, &
      harmonic_mean,reference_heat_flux)
    character(len=*), intent(in) :: input_path
    logical, intent(out) :: enabled,harmonic_mean
    logical, optional, intent(out) :: reference_heat_flux
    logical :: pinned_heat_flux
    real(dp), intent(out) :: cs,prt,sct,ci
    character(len=*), intent(out) :: model,test_model
    integer, intent(out) :: test_ratio,grid_ratio
    integer :: unit,ios
    namelist /reactive_les/ enabled,cs,prt,sct,ci,model,test_model,test_ratio,grid_ratio,harmonic_mean,pinned_heat_flux
    enabled=.false.;cs=.17_dp;prt=.7_dp;sct=.7_dp;ci=0;pinned_heat_flux=.false.
    model='smagorinsky';test_model='box';test_ratio=2;grid_ratio=1;harmonic_mean=.false.
    open(newunit=unit,file=input_path,status='old',action='read',iostat=ios)
    if(ios/=0) error stop 'Cannot open LES input'
    read(unit,nml=reactive_les,iostat=ios)
    close(unit)
    if(ios>0) error stop 'Invalid reactive_les namelist'
    if(pinned_heat_flux.and.model/='dynamic') error stop 'Pinned heat flux requires dynamic LES'
    if(present(reference_heat_flux)) reference_heat_flux=pinned_heat_flux
  end subroutine

  subroutine run_reactive_3d_application( &
      input_path, application_label, bundle_sha256, config, species, &
      reactions, transport, base_mole_fractions, chemistry_integrator,evaluate_reaction,evaluate_diffusion)
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    procedure(application_diffusion), optional :: evaluate_diffusion
    real(dp), allocatable :: constant_coefficients(:),hit_parameters(:)
    real(dp) :: hit_properties(3)
    real(dp), allocatable :: sdc_force(:),diffusion_rates(:,:,:,:),old_boundary_rate(:),new_boundary_rate(:)
    real(dp) :: extra_diffusivity,maximum_dt,linear_coefficient,momentum_offset(3)
    character(len=1024) :: restart_file,checkpoint_file
    real(dp) :: checkpoint_metrics(2),initial_timestep_factor,maximum_timestep_growth,previous_timestep,initial_dt,fixed_dt
    integer :: stop_after_steps
    integer :: source_calls
    logical :: les_enabled,les_harmonic,les_reference_heat
    real(dp), allocatable :: checkpoint_initial_state(:,:,:,:),checkpoint_initial_temperature(:,:,:)
    type(dual_energy_policy) :: energy_policy
    real(dp), allocatable :: auxiliary_internal_energy(:,:,:),energy_reset_correction
    real(dp), allocatable :: corrected_totals(:)
    real(dp) :: les_cs,les_prt,les_sct,les_ci
    integer :: les_test_ratio,les_grid_ratio
    character(len=32) :: les_model,les_filter
    character(len=*), intent(in) :: input_path, application_label
    character(len=*), intent(in) :: bundle_sha256
    type(reactive_3d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    type(nasa7_species), allocatable :: eos_species(:)
    logical :: eos_policy_ok
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: base_mole_fractions(:)
    character(len=*), intent(in), optional :: chemistry_integrator

    type(gas_transport_species), allocatable :: sdc_transport(:)
    type(elementary_reaction), allocatable :: sdc_reactions(:)
    real(dp), allocatable :: sdc_reaction_source(:,:,:,:),sdc_heat_release(:,:,:)
    real(dp), allocatable :: x(:), y(:), z(:)
    real(dp), allocatable, target :: state(:, :, :, :), temperature(:, :, :)
    real(dp), allocatable, target :: ordinary_state(:,:,:,:)
    real(dp), pointer :: thermodynamic_state(:,:,:,:)
    real(dp), allocatable :: rotating_potential(:,:,:)
    type(rotating_frame_3d), allocatable :: frame
    real(dp) :: domain_lower(3),cell_position(3)
    real(dp), allocatable :: initial_totals(:), final_totals(:), boundary_totals(:), step_boundary_change(:)
    character(len=32) :: faces(6)
    real(dp), allocatable :: inflow_state(:,:),inflow_y(:),inflow_x(:),inlet_cube(:,:,:,:)
    real(dp), allocatable :: inlet_ghost(:,:,:,:),inlet_ghost_t(:,:,:)
    real(dp), allocatable :: external_force(:),source_totals(:),step_source_change(:)
    real(dp) :: mach_pressure(6),difmag
    real(dp) :: inflow_t(6),inlet_temperature(2,2,2)
    real(dp) :: inlet_dt,inlet_transport_dt,inlet_diffusivity,face_dt,face_diffusivity,wall_dt
    real(dp), allocatable :: mass_fractions(:), primitive(:)
    real(dp) :: dx, dy, dz, time, dt, hydro_dt, transport_dt
    real(dp) :: base_density, density_l1_error, exact_density
    real(dp) :: minimum_density, maximum_density
    real(dp) :: minimum_pressure, maximum_pressure
    real(dp) :: minimum_temperature, maximum_temperature
    real(dp) :: maximum_speed, maximum_closure_error
    real(dp) :: maximum_euler_conservation_error
    real(dp) :: maximum_all_component_conservation_error
    real(dp) :: maximum_species_change, local_temperature, sound_speed
    real(dp) :: maximum_diffusivity, step_maximum_diffusivity
    real(dp) :: minimum_transport_theta, step_transport_theta
    character(len=1024) :: message
    logical :: ok, cell_ok
    integer :: step, i, j, k, nvar, face

    eos_species=species
    call configure_nasa7_eos_policy(eos_species,config%pinned_eos_extrapolation,eos_policy_ok, &
      config%pinned_eos_temperature_min,config%pinned_eos_temperature_max)
    if(.not.eos_policy_ok) error stop 'Invalid application EOS policy'
    call read_rotation_controls()
    domain_lower=[config%x_lower,config%y_lower,config%z_lower]
    call read_restart_controls()
    call read_regular_timestep_limit(input_path,maximum_dt,initial_timestep_factor,maximum_timestep_growth, &
      initial_dt,fixed_dt)
    call read_linear_forcing(input_path,linear_coefficient,momentum_offset)
    if(linear_coefficient/=0.and.(config%time_integrator/='sdc'.and.config%time_integrator/='mol')) &
      error stop 'Linear forcing requires SDC or MOL'
    call read_regular_les_controls(input_path,les_enabled,les_cs,les_prt,les_sct,les_ci, &
        les_model,les_filter,les_test_ratio,les_grid_ratio,les_harmonic,les_reference_heat)
    if(les_enabled.and..not.present(evaluate_diffusion)) error stop 'Regular LES backend is unavailable'
    if(present(evaluate_diffusion).and.(config%time_integrator/='sdc'.and.config%time_integrator/='mol')) &
      error stop 'Additional diffusion requires SDC or MOL'
    if((config%sdc_reactor=='cvode').neqv.present(evaluate_reaction)) &
      error stop 'SDC reactor selection does not match supplied backend'
    if (size(eos_species) < 1 .or. (config%chemistry_enabled.and.size(reactions)<1) .or. &
        size(transport) /= size(eos_species) .or. &
        size(base_mole_fractions) /= size(eos_species)) then
      error stop "Reactive 3D application data have incompatible dimensions"
    end if
    nvar = reactive_nvar(size(eos_species))
    allocate(x(config%nx), y(config%ny), z(config%nz))
    allocate(state(nvar, config%nx, config%ny, config%nz))
    allocate(temperature(config%nx, config%ny, config%nz))
    allocate(initial_totals(nvar), final_totals(nvar), boundary_totals(nvar), step_boundary_change(nvar))
    if((config%time_integrator=='sdc'.or.config%time_integrator=='mol')) &
      allocate(sdc_heat_release(config%nx,config%ny,config%nz),source=0._dp)
    if((config%time_integrator=='sdc'.or.config%time_integrator=='mol').and.config%transport_enabled) sdc_transport=transport
    if((config%time_integrator=='sdc'.or.config%time_integrator=='mol').and.config%chemistry_enabled) then
      sdc_reactions=reactions
      allocate(sdc_reaction_source(nvar-4,config%nx,config%ny,config%nz));sdc_reaction_source=0
    end if
    if(present(evaluate_diffusion)) then
      allocate(diffusion_rates, mold=state)
      allocate(old_boundary_rate(nvar),new_boundary_rate(nvar))
    else if(linear_coefficient==0) then
      sdc_force=config%external_forcing
    end if
    boundary_totals=0.0_dp
    allocate(source_totals(nvar),step_source_change(nvar));source_totals=0
    if(any(config%external_forcing/=0)) external_force=config%external_forcing
    faces=resolved_reactive_boundaries_3d(config)
    call read_constant_transport(input_path,constant_coefficients,hit_parameters)
    if(config%thermo_model=='gamma_law'.and.config%transport_enabled.and..not.allocated(constant_coefficients)) &
      error stop 'Gamma-law transport requires explicit constant transport'
    if(allocated(constant_coefficients)) then
      if((config%time_integrator/='sdc'.and.config%time_integrator/='mol').or..not.config%transport_enabled) &
        error stop 'Constant transport requires SDC with transport enabled'
    end if
    allocate(mass_fractions(size(eos_species)))
    allocate(primitive(reactive_nprim(size(eos_species))))
    call uniform_cell_centers_3d( &
      config%nx, config%ny, config%nz, &
      config%x_lower, config%x_upper, &
      config%y_lower, config%y_upper, &
      config%z_lower, config%z_upper, x, y, z, dx, dy, dz)
    call initialize_reactive_problem_3d( &
      eos_species, config, x, y, z, state, temperature, base_density, &
      mass_fractions, ok, base_mole_fractions)
    if (.not. ok) error stop "Failed to initialize reactive 3D problem"
    call initialize_reactive_hit_3d(input_path,config,eos_species,x,y,z,state,temperature,ok,hit_properties)
    if(.not.ok) error stop 'Invalid HIT velocity input'
    if(allocated(frame)) then
      allocate(rotating_potential(config%nx,config%ny,config%nz))
      do k=1,config%nz
      do j=1,config%ny
      do i=1,config%nx
        cell_position=[x(i),y(j),z(k)]
        call rotating_potential_3d(frame,cell_position,rotating_potential(i,j,k),ok)
        if(.not.ok) error stop 'Invalid regular rotating potential'
      enddo
      enddo
      enddo
      call convert_regular_energy_3d(frame,domain_lower,[dx,dy,dz],state,0,.true.,ok)
      if(.not.ok) error stop 'Invalid initial regular rotating energy'
    endif
    energy_policy%allow_small_energy=config%allow_small_energy
    energy_policy%allow_negative_energy=config%allow_negative_energy
    energy_policy%update_total_from_internal=config%dual_energy_update_E_from_e
    energy_policy%eta2=config%dual_energy_eta2
    if(config%time_integrator=='mol') then
      allocate(auxiliary_internal_energy(config%nx,config%ny,config%nz),energy_reset_correction)
      auxiliary_internal_energy=state(5,:,:,:)-.5_dp*sum(state(2:4,:,:,:)**2,dim=1)/state(1,:,:,:)
      if(allocated(rotating_potential)) auxiliary_internal_energy= &
        auxiliary_internal_energy+state(1,:,:,:)*rotating_potential
      energy_reset_correction=0
    endif
    call refresh_ordinary_state()
    if(len_trim(restart_file)>0.or.len_trim(checkpoint_file)>0) then
      checkpoint_initial_state=state;checkpoint_initial_temperature=temperature
    endif
    if(allocated(hit_parameters)) then
      call resolve_hit_transport(hit_properties,hit_parameters,constant_coefficients)
      write(*,'(a,3es24.16)') 'HIT density/urms/cp: ',hit_properties
    end if
    if(allocated(constant_coefficients)) &
      write(*,'(a,2es24.16)') 'Constant viscosity/conductivity: ',constant_coefficients(1:2)
    call reactive_integrals_3d( &
      state, config%nx, config%ny, config%nz, dx, dy, dz, &
      initial_totals, ok)
    if (.not. ok) error stop "Failed to integrate initial reactive 3D state"

    allocate(inflow_state(nvar,6),inflow_y(size(eos_species)),inflow_x(size(eos_species)))
    allocate(inlet_cube(nvar,2,2,2))
    inflow_state=0.0_dp; inflow_t=config%inflow_temperature
    mach_pressure=merge(config%inflow_pressure,0.0_dp,config%inflow_mach>=0)
    if (any(config%inflow_mach>=0)) then
      allocate(inlet_ghost(nvar,config%nx+4,config%ny+4,config%nz+4))
      allocate(inlet_ghost_t(config%nx+4,config%ny+4,config%nz+4))
    end if
    inlet_dt=huge(1.0_dp); inlet_transport_dt=huge(1.0_dp); inlet_diffusivity=0.0_dp
    inflow_y=mass_fractions
    if (any(faces=='inflow').and.config%inflow_composition_count/=0) then
      call resolve_selected_composition('Inflow',config%inflow_composition_count, &
        config%inflow_composition_species,config%inflow_composition_mole_fractions,eos_species,inflow_x,ok,message)
      if (.not.ok) then
        write(*,'(a)') trim(message)
        error stop 'Invalid inflow composition'
      end if
      call mass_fractions_from_mole_fractions(eos_species,inflow_x,inflow_y,ok)
      if (.not.ok) error stop 'Invalid inflow mass fractions'
    end if
    do face=1,6
      if (faces(face)/='inflow'.or.config%inflow_mach(face)>=0) cycle
      primitive=0.0_dp
      primitive(1)=mixture_density(eos_species,inflow_y,config%inflow_pressure(face),config%inflow_temperature(face),ok)
      if (.not.ok) error stop 'Invalid inflow density'
      primitive(2:4)=config%inflow_velocity(:,face); primitive(5)=config%inflow_pressure(face)
      primitive(6:)=inflow_y
      call reactive_primitive_to_conserved(eos_species,primitive,inflow_state(:,face),inflow_t(face),sound_speed,ok)
      if (.not.ok) error stop 'Invalid inflow thermodynamics'
      inlet_dt=min(inlet_dt,config%cfl/sum((abs(primitive(2:4))+sound_speed)/[dx,dy,dz]))
      if (config%transport_enabled) then
        do k=1,2
          do j=1,2
            do i=1,2
              inlet_cube(:,i,j,k)=inflow_state(:,face)
            end do
          end do
        end do
        inlet_temperature=inflow_t(face)
        call reactive_transport_timestep_3d(eos_species,transport,inlet_cube,inlet_temperature,2,2,2,dx,dy,dz, &
          config%transport_cfl,config%viscosity_enabled,config%thermal_conduction_enabled, &
          config%species_diffusion_enabled,face_dt,face_diffusivity,ok,constant_coefficients=constant_coefficients)
        if (.not.ok) error stop 'Invalid inflow transport timestep'
        inlet_transport_dt=min(inlet_transport_dt,face_dt)
        inlet_diffusivity=max(inlet_diffusivity,face_diffusivity)
      end if
    end do

    write(*, '(a)') &
      "PeleF " // pelef_version // " " // trim(application_label)
    if (len_trim(bundle_sha256) > 0) then
      write(*, '(a,1x,a)') "Bundle SHA-256:", trim(bundle_sha256)
      if((config%time_integrator=='sdc'.or.config%time_integrator=='mol').and.config%chemistry_enabled) then
        write(*,'(a,a)') 'SDC reactor: ',trim(config%sdc_reactor)
      else if (present(chemistry_integrator)) then
        write(*, '(a,1x,a)') "Chemistry integrator:", &
          trim(chemistry_integrator)
      end if
    end if
    write(*, '(a,1x,a)') "Input:", trim(input_path)
    write(*, '(a,i0,a,i0,a,i0)') &
      "Grid: nx=", config%nx, ", ny=", config%ny, ", nz=", config%nz
    write(*, '(a,1x,a)') "Problem:", trim(config%problem)
    write(*, '(a,1x,a)') "Thermodynamics:", trim(config%thermo_model)
    write(*,'(a,a)') 'Time integrator: ',trim(config%time_integrator)
    if((config%time_integrator=='sdc'.or.config%time_integrator=='mol').and.len_trim(bundle_sha256)==0) &
      write(*,'(a,a)') 'SDC reactor: ',trim(config%sdc_reactor)
    if((config%time_integrator=='sdc'.or.config%time_integrator=='mol')) write(*,'(a,i0)') 'SDC iterations: ',config%sdc_iterations
    difmag=resolved_reactive_difmag_3d(config)
    write(*,'(a,es24.16)') 'Artificial-viscosity coefficient: ',difmag
    write(*, '(a,1x,a)') "Reconstruction:", trim(config%reconstruction)
    write(*, '(a,1x,a)') "Limiter:", trim(config%limiter)
    write(*, '(a,1x,a)') "Riemann solver:", trim(config%riemann_solver)
    write(*, '(a,1x,a)') &
      "Boundary condition:", trim(config%boundary_condition)
    write(*, '(a,6(1x,a))') "Boundary faces (x-/x+/y-/y+/z-/z+):", faces
    write(*, '(a,l2)') "Chemistry: ", config%chemistry_enabled
    if (config%chemistry_enabled) then
      write(*, '(a,es12.5)') &
        "Chemistry relative tolerance: ", &
        config%chemistry_relative_tolerance
      write(*, '(a,es12.5)') &
        "Chemistry absolute tolerance: ", &
        config%chemistry_absolute_tolerance
    end if
    write(*, '(a,l2)') "Molecular transport: ", config%transport_enabled
    if (config%transport_enabled) then
      write(*, '(a,l2)') "Viscosity: ", config%viscosity_enabled
      write(*, '(a,l2)') &
        "Thermal conduction: ", config%thermal_conduction_enabled
      write(*, '(a,l2)') &
        "Species diffusion: ", config%species_diffusion_enabled
      write(*, '(a,l2)') &
        "Barodiffusion: ", config%barodiffusion_enabled
      write(*, '(a,es12.5)') "Transport CFL: ", config%transport_cfl
    end if

    time = 0.0_dp
    previous_timestep=0
    step = 0
    maximum_diffusivity = 0.0_dp
    minimum_transport_theta = 1.0_dp
    if(len_trim(restart_file)>0) then
      checkpoint_metrics=[maximum_diffusivity,minimum_transport_theta]
      call read_regular_gas_checkpoint(trim(restart_file),checkpoint_identity,eos_species,state,temperature,time,step, &
        initial_totals,boundary_totals,source_totals,checkpoint_metrics,ok,message,sdc_reaction_source,sdc_heat_release, &
        spacing=[dx,dy,dz],final_time=config%final_time,rotating_potential=rotating_potential, &
        previous_timestep=previous_timestep,internal_energy_density=auxiliary_internal_energy, &
        energy_correction=energy_reset_correction)
      if(.not.ok) then
        write(*,'(a)') trim(message);error stop 'Regular gas restart rejected'
      endif
      if(time>config%final_time) error stop 'Regular gas restart exceeds final time'
      maximum_diffusivity=checkpoint_metrics(1);minimum_transport_theta=checkpoint_metrics(2)
    endif
    do while (time < config%final_time)
      if(stop_after_steps>0.and.step>=stop_after_steps) exit
      if (step >= config%maximum_steps) then
        error stop "Maximum step count reached before reactive 3D final_time"
      end if
      call refresh_ordinary_state()
      if(step==0.and.initial_dt>0) then
        dt=initial_dt
      else if(fixed_dt>0) then
        dt=fixed_dt
      else
        call compute_reactive_cfl_timestep_3d( &
          eos_species, thermodynamic_state, temperature, config%nx, config%ny, config%nz, &
          dx, dy, dz, config%cfl, hydro_dt, ok)
        if (.not. ok) error stop "Failed to compute reactive 3D CFL timestep"
        if (any(config%inflow_mach>=0)) then
          call fill_mach_inlet_boundary_ghosts_3d(eos_species,thermodynamic_state,temperature,faces,config%inflow_mach, &
            inlet_ghost,inlet_ghost_t,ok,inflow_state=inflow_state,inflow_temperature=inflow_t, &
            inflow_pressure=mach_pressure,viscous_wall=.true.,wall_velocity=config%wall_velocity)
          if (.not.ok) error stop 'Failed to construct Mach inlet timestep states'
          call compute_reactive_cfl_timestep_3d(eos_species,inlet_ghost,inlet_ghost_t, &
            config%nx+4,config%ny+4,config%nz+4,dx,dy,dz,config%cfl,inlet_dt,ok)
          if (.not.ok) error stop 'Failed to compute Mach inlet CFL timestep'
          if (config%transport_enabled) then
            call reactive_transport_timestep_3d(eos_species,transport,inlet_ghost,inlet_ghost_t, &
              config%nx+4,config%ny+4,config%nz+4,dx,dy,dz,config%transport_cfl, &
              config%viscosity_enabled,config%thermal_conduction_enabled,config%species_diffusion_enabled, &
              inlet_transport_dt,inlet_diffusivity,ok,constant_coefficients=constant_coefficients)
            if (.not.ok) error stop 'Failed to compute Mach inlet transport timestep'
          end if
        end if
        dt = min(hydro_dt,inlet_dt)
        if (config%transport_enabled) then
          call reactive_transport_timestep_3d( &
            eos_species, transport, thermodynamic_state, temperature, &
            config%nx, config%ny, config%nz, dx, dy, dz, &
            config%transport_cfl, config%viscosity_enabled, &
            config%thermal_conduction_enabled, &
            config%species_diffusion_enabled, transport_dt, &
            step_maximum_diffusivity, ok,constant_coefficients=constant_coefficients)
          if (.not. ok) then
            error stop "Failed to compute reactive 3D transport step"
          end if
          if (any(faces=='isothermal_wall')) then
            call reactive_isothermal_timestep_3d(eos_species,transport,thermodynamic_state,temperature,[dx,dy,dz], &
              faces,config%wall_temperature,config%transport_cfl,wall_dt,ok,constant_coefficients)
            if (.not.ok) error stop 'Invalid isothermal wall timestep'
            ! Reserve half the interior diffusion allowance for boundary exchange.
            dt=min(dt,0.5_dp*transport_dt,0.5_dp*inlet_transport_dt,wall_dt)
          else
            dt = min(dt, transport_dt, inlet_transport_dt)
          end if
          maximum_diffusivity = max( &
            maximum_diffusivity, step_maximum_diffusivity, inlet_diffusivity)
        end if
        if(present(evaluate_diffusion)) then
          call evaluate_diffusion(thermodynamic_state,temperature,[dx,dy,dz],faces,inflow_state,inflow_t, &
            diffusion_rates,new_boundary_rate,extra_diffusivity,ok)
          if(.not.ok) error stop 'Failed to compute additional diffusion timestep'
          if(.not.ieee_is_finite(extra_diffusivity).or.extra_diffusivity<0) error stop 'Invalid additional diffusivity'
          if(extra_diffusivity>0) dt=min(dt,.1_dp/(extra_diffusivity*sum(1._dp/[dx,dy,dz]**2)))
          maximum_diffusivity=max(maximum_diffusivity,extra_diffusivity)
        end if
        dt=min(dt,maximum_dt)
      endif
      if(step==0) then
        if(initial_dt<=0) dt=initial_timestep_factor*dt
      else if(fixed_dt<=0) then
        if(previous_timestep<dt/maximum_timestep_growth) dt=maximum_timestep_growth*previous_timestep
      endif
      if(config%final_time-time-dt<.001_dp*dt) dt=config%final_time-time
      dt=min(dt,config%final_time-time)
      if(time+dt<=time) error stop 'Maximum timestep cannot advance time'
      if((config%time_integrator=='sdc'.or.config%time_integrator=='mol')) then
        step_boundary_change=0;step_transport_theta=1;source_calls=0
        if(present(evaluate_diffusion).or.linear_coefficient/=0) then
          call perform_sdc(evaluate_sources)
        else
          call perform_sdc()
        end if
      else
      call advance_reactive_full_3d( &
        eos_species, reactions, transport, state, temperature, &
        config%nx, config%ny, config%nz, dx, dy, dz, dt, &
        config%riemann_solver, config%chemistry_enabled, &
        config%chemistry_relative_tolerance, &
        config%chemistry_absolute_tolerance, config%transport_enabled, &
        config%viscosity_enabled, config%thermal_conduction_enabled, &
        config%species_diffusion_enabled, config%barodiffusion_enabled, &
        step_transport_theta, ok, config%reconstruction, config%limiter, &
        chemistry_integrator=chemistry_integrator, boundary_faces=faces, boundary_change=step_boundary_change, &
        inflow_state=inflow_state,inflow_temperature=inflow_t,wall_temperature=config%wall_temperature, &
        wall_velocity=config%wall_velocity,outlet_pressure=config%outlet_pressure, &
        inflow_mach=config%inflow_mach,inflow_pressure=mach_pressure,difmag=difmag, &
        external_force=external_force,source_change=step_source_change)
      end if
      if (.not. ok) then
        error stop "Reactive 3D update rejected its candidate"
      end if
      if(present(evaluate_diffusion)) then
        if(config%time_integrator=='mol'.and.config%chemistry_enabled.and.config%sdc_iterations>1) new_boundary_rate=0
        step_boundary_change=step_boundary_change+.5_dp*dt*(old_boundary_rate+new_boundary_rate)
        step_source_change=step_source_change-.5_dp*dt*(old_boundary_rate+new_boundary_rate)
      end if
      source_totals=source_totals+step_source_change
      boundary_totals = boundary_totals + step_boundary_change
      minimum_transport_theta = min( &
        minimum_transport_theta, step_transport_theta)
      previous_timestep=dt
      time = time + dt
      step = step + 1
    end do

    if(len_trim(checkpoint_file)>0) then
      checkpoint_metrics=[maximum_diffusivity,minimum_transport_theta]
      call write_regular_gas_checkpoint(trim(checkpoint_file),checkpoint_identity,eos_species,state,temperature,time,step, &
        initial_totals,boundary_totals,source_totals,checkpoint_metrics,ok,message,sdc_reaction_source,sdc_heat_release, &
        rotating_potential=rotating_potential,previous_timestep=previous_timestep, &
        internal_energy_density=auxiliary_internal_energy,energy_correction=energy_reset_correction)
      if(.not.ok) then
        write(*,'(a)') trim(message);error stop 'Regular gas checkpoint write rejected'
      endif
    endif
    call reactive_integrals_3d( &
      state, config%nx, config%ny, config%nz, dx, dy, dz, final_totals, ok)
    if (.not. ok) error stop "Failed to integrate final reactive 3D state"
    corrected_totals=final_totals
    if(allocated(energy_reset_correction)) corrected_totals(5)=corrected_totals(5)-energy_reset_correction
    maximum_euler_conservation_error = maxval( &
      abs(corrected_totals(1:ncons) - initial_totals(1:ncons) - boundary_totals(1:ncons) - source_totals(1:ncons)) / &
        max(1.0_dp, abs(initial_totals(1:ncons))))
    maximum_all_component_conservation_error = maxval( &
      abs(corrected_totals - initial_totals - boundary_totals - source_totals) / max(1.0_dp, abs(initial_totals)))
    maximum_species_change = maxval(abs( &
      final_totals(6:nvar) - initial_totals(6:nvar)))
    call refresh_ordinary_state()
    call reactive_extrema_3d( &
      eos_species, thermodynamic_state, temperature, config%nx, config%ny, config%nz, &
      minimum_density, maximum_density, minimum_pressure, maximum_pressure, &
      minimum_temperature, maximum_temperature, maximum_speed, &
      maximum_closure_error, ok,internal_energy_density=auxiliary_internal_energy)
    if (.not. ok) error stop "Final reactive 3D state is not physical"

    if (trim(config%problem) == "entropy_wave" .and. all(faces=="periodic").and. &
        .not.allocated(external_force).and..not.allocated(frame)) then
      density_l1_error = 0.0_dp
      do k = 1, config%nz
        do j = 1, config%ny
          do i = 1, config%nx
            call reactive_conserved_to_primitive( &
              eos_species, state(:, i, j, k), temperature(i, j, k), primitive, &
              local_temperature, sound_speed, cell_ok)
            if (.not. cell_ok) error stop "Invalid final reactive 3D cell"
            exact_density = reactive_entropy_wave_density_3d( &
              x(i), y(j), z(k), time, config, base_density)
            density_l1_error = density_l1_error + &
              abs(primitive(1) - exact_density)
          end do
        end do
      end do
      density_l1_error = density_l1_error / &
        real(config%nx * config%ny * config%nz, dp)
    end if

    if (maximum_euler_conservation_error > 5.0e-11_dp) then
      error stop "Reactive 3D boundary-accounted conservation gate failed"
    end if
    if (.not. config%chemistry_enabled .and. &
        maximum_all_component_conservation_error > 5.0e-11_dp) then
      error stop "Inert reactive 3D species conservation gate failed"
    end if
    if (maximum_closure_error > 5.0e-11_dp) then
      error stop "Reactive 3D species-closure gate failed"
    end if
    write(*, '(a,*(es24.16,1x))') "Integrated boundary change: ", boundary_totals
    write(*, '(a,*(es24.16,1x))') "Integrated source change: ", source_totals
    if(allocated(energy_reset_correction)) write(*,'(a,es24.16)') &
      "Integrated energy reset correction: ",energy_reset_correction
    call write_reactive_3d_csv( &
      trim(config%output_file), eos_species, x, y, z, state, temperature, &
      config%nx, config%ny, config%nz, time, ok, message,heat_release=sdc_heat_release, &
      rotating_potential=rotating_potential,internal_energy_density=auxiliary_internal_energy)
    if (.not. ok) then
      write(*, '(a)') trim(message)
      error stop 3
    end if

    write(*, '(a,i0)') "Completed steps: ", step
    write(*, '(a,es24.16)') "Final time: ", time
    write(*, '(a,es24.16)') "Base density: ", base_density
    if (trim(config%problem) == "entropy_wave" .and. all(faces=="periodic").and..not.allocated(external_force)) then
      write(*, '(a,es24.16)') "Density L1 error: ", density_l1_error
    else
      write(*, '(a)') "Density L1 error: not applicable"
    end if
    write(*, '(a,es24.16)') "Minimum density: ", minimum_density
    write(*, '(a,es24.16)') "Maximum density: ", maximum_density
    write(*, '(a,es24.16)') "Minimum pressure: ", minimum_pressure
    write(*, '(a,es24.16)') "Maximum pressure: ", maximum_pressure
    write(*, '(a,es24.16)') "Minimum temperature: ", minimum_temperature
    write(*, '(a,es24.16)') "Maximum temperature: ", maximum_temperature
    write(*, '(a,es24.16)') "Maximum speed: ", maximum_speed
    write(*, '(a,es24.16)') &
      "Maximum species closure error: ", maximum_closure_error
    write(*, '(a,es24.16)') &
      "Maximum normalized Euler conserved change: ", &
      maximum_euler_conservation_error
    write(*, '(a,es24.16)') &
      "Maximum normalized all-component change: ", &
      maximum_all_component_conservation_error
    write(*, '(a,es24.16)') &
      "Maximum absolute species integral change: ", maximum_species_change
    write(*, '(a,es24.16)') &
      "Maximum transport diffusivity: ", maximum_diffusivity
    write(*, '(a,es24.16)') &
      "Minimum transport theta: ", minimum_transport_theta
    write(*, '(a,es24.16)') &
      "Mass change: ", final_totals(irho) - initial_totals(irho)
    write(*, '(a,es24.16)') &
      "X-momentum change: ", final_totals(imx) - initial_totals(imx)
    write(*, '(a,es24.16)') &
      "Y-momentum change: ", final_totals(imy) - initial_totals(imy)
    write(*, '(a,es24.16)') &
      "Z-momentum change: ", final_totals(imz) - initial_totals(imz)
    write(*, '(a,es24.16)') &
      "Energy change: ", final_totals(iet) - initial_totals(iet)
    write(*, '(a,1x,a)') "Output:", trim(config%output_file)
  contains
    subroutine read_restart_controls()
      integer :: unit,status
      namelist /reactive_restart/ restart_file,checkpoint_file,stop_after_steps
      restart_file='';checkpoint_file='';stop_after_steps=0
      open(newunit=unit,file=input_path,status='old',action='read',iostat=status)
      if(status/=0) error stop 'Cannot open regular restart input'
      read(unit,nml=reactive_restart,iostat=status)
      close(unit)
      if(status>0.or.stop_after_steps<0) error stop 'Invalid regular restart controls'
    end subroutine
    subroutine read_rotation_controls()
      call read_reactive_rotation_controls(input_path,config%time_integrator,frame)
    end subroutine read_rotation_controls
    subroutine refresh_ordinary_state()
      thermodynamic_state=>state
      if(.not.allocated(frame)) return
      ordinary_state=state
      call convert_regular_energy_3d(frame,domain_lower,[dx,dy,dz],ordinary_state,0,.false.,ok)
      if(.not.ok) error stop 'Invalid regular rotating thermodynamic state'
      thermodynamic_state=>ordinary_state
    end subroutine refresh_ordinary_state
    subroutine checkpoint_identity(unit,status)
      integer, intent(in) :: unit
      integer, intent(out) :: status
      type(reactive_3d_config) :: physics
      real(dp) :: constant_values(4),hit_values(3)
      integer :: constant_count,hit_count,record,scratch
      character(len=32) :: integrator
      namelist /regular_checkpoint_physics/ physics,maximum_dt,initial_timestep_factor,maximum_timestep_growth, &
        initial_dt,fixed_dt, &
        linear_coefficient,momentum_offset, &
        constant_count,constant_values,hit_count,hit_values,les_enabled,les_cs,les_prt,les_sct,les_ci, &
        les_model,les_filter,les_test_ratio,les_grid_ratio,les_harmonic,les_reference_heat,integrator,bundle_sha256
      namelist /regular_checkpoint_materials/ species
      namelist /regular_checkpoint_initial/ checkpoint_initial_state,checkpoint_initial_temperature
      physics=config;physics%final_time=0;physics%maximum_steps=0;physics%output_file=''
      constant_count=0;constant_values=0;hit_count=0;hit_values=0
      if(allocated(constant_coefficients)) then
        constant_count=size(constant_coefficients);constant_values(:constant_count)=constant_coefficients
      endif
      if(allocated(hit_parameters)) then
        hit_count=size(hit_parameters);hit_values(:hit_count)=hit_parameters
      endif
      integrator='implicit';if(present(chemistry_integrator)) integrator=chemistry_integrator
      write(unit,'(a)',iostat=status) 'REGULAR_CHECKPOINT_PHYSICS_V3'
      if(config%time_integrator=='mol'.and.status==0) &
        write(unit,'(a)',iostat=status) 'REGULAR_MOL_AUXILIARY_ENERGY_V2'
      if(status==0) then
        open(newunit=scratch,status='scratch',action='readwrite',iostat=status)
        if(status/=0) return
        write(scratch,nml=regular_checkpoint_physics,iostat=status)
        call copy_legacy_namelist(scratch,unit,status)
        close(scratch)
      end if
      if(allocated(frame).and.status==0) then
        write(unit,'(a)',iostat=status) 'REGULAR_ROTATING_FRAME_V1'
        if(status==0) write(unit,'(i0,1x,4(es27.18e3,1x))',iostat=status) &
          frame%axis,frame%angular_speed,frame%axis_point
      endif
      if(status==0) then
        open(newunit=scratch,status='scratch',action='readwrite',iostat=status)
        if(status/=0) return
        write(scratch,nml=regular_checkpoint_materials,iostat=status)
        call copy_legacy_namelist(scratch,unit,status)
        close(scratch)
      end if
      if(status==0) call write_selected_checkpoint_context(unit,reactions,bundle_sha256,integrator, &
        base_mole_fractions,config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance,status, &
        fixed_model=config%thermo_model/='selected')
      do record=1,size(transport)
        if(status/=0) return
        write(unit,'(a,1x,i0)',iostat=status) trim(transport(record)%name),transport(record)%geometry
        if(status==0) write(unit,'(*(es27.18e3,1x))',iostat=status) transport(record)%well_depth, &
          transport(record)%diameter,transport(record)%dipole,transport(record)%polarizability, &
          transport(record)%rotational_relaxation,transport(record)%simple_viscosity,transport(record)%simple_conductivity
        if(status==0) write(unit,'(l1)',iostat=status) allocated(transport(record)%constant_coefficients)
        if(status==0.and.allocated(transport(record)%constant_coefficients)) &
          write(unit,'(*(es27.18e3,1x))',iostat=status) transport(record)%constant_coefficients
        if(status==0) write(unit,'(l1)',iostat=status) allocated(transport(record)%simple_binary)
        if(status==0) call write_simple_transport_fits(unit,transport(record),status)
      enddo
      if(status==0) write(unit,nml=regular_checkpoint_initial,iostat=status)
    end subroutine
    subroutine copy_legacy_namelist(scratch,unit,status)
      integer, intent(in) :: scratch,unit
      integer, intent(inout) :: status
      integer :: ios
      character(len=65536) :: line
      if(status/=0) return
      rewind(scratch,iostat=status)
      if(status/=0) return
      do
        read(scratch,'(a)',iostat=ios) line
        if(ios<0) exit
        if(ios/=0) then
          status=ios;return
        endif
        ! The EOS policy is carried by the versioned outer identity record.
        ! Preserve legacy component names and omit only the new policy field.
        if(index(line,'%PINNED_EOS_')>0) cycle
        write(unit,'(a)',iostat=status) trim(line)
        if(status/=0) return
      enddo
    end subroutine copy_legacy_namelist
    subroutine perform_sdc(source_callback)
      procedure(reactive_sdc_source_evaluator), optional :: source_callback
      if(config%time_integrator=='mol') then
        call advance_reactive_mol_3d(eos_species,state,temperature,[dx,dy,dz],dt,time,config%sdc_iterations, &
          config%chemistry_enabled,config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance,ok, &
          reactions=sdc_reactions,reaction_source=sdc_reaction_source,heat_release=sdc_heat_release, &
          source_count=merge(1+merge(1,0,present(evaluate_diffusion)),0,present(source_callback)), &
          evaluate_sources=source_callback,external_force=sdc_force,transport=sdc_transport, &
          transport_flags=[config%viscosity_enabled,config%thermal_conduction_enabled, &
            config%species_diffusion_enabled,config%barodiffusion_enabled],minimum_transport_theta=step_transport_theta, &
          boundary_change=step_boundary_change,source_change=step_source_change,boundary_faces=faces, &
          inflow_state=inflow_state,inflow_temperature=inflow_t,inflow_mach=config%inflow_mach,inflow_pressure=mach_pressure, &
          outlet_pressure=config%outlet_pressure,wall_temperature=config%wall_temperature,wall_velocity=config%wall_velocity, &
          constant_coefficients=constant_coefficients,evaluate_reaction=evaluate_reaction, &
          lax_friedrichs=config%mol_lax_friedrichs,domain_lower=domain_lower,frame=frame, &
          internal_energy_density=auxiliary_internal_energy,energy_correction=energy_reset_correction,energy_policy=energy_policy)
      else
        call advance_reactive_sdc_3d(eos_species,state,temperature,[dx,dy,dz],dt,time, &
          config%sdc_iterations,merge(2,1,config%transport_enabled)+merge(1,0,present(evaluate_diffusion)), &
          evaluate_sources=source_callback,ok=ok,difmag=difmag,external_force=sdc_force, &
          source_change=step_source_change,transport=sdc_transport, &
          transport_flags=[config%viscosity_enabled,config%thermal_conduction_enabled, &
            config%species_diffusion_enabled,config%barodiffusion_enabled],minimum_transport_theta=step_transport_theta, &
          reactions=sdc_reactions,reaction_source=sdc_reaction_source, &
          chemistry_relative_tolerance=config%chemistry_relative_tolerance, &
          chemistry_absolute_tolerance=config%chemistry_absolute_tolerance, &
          boundary_faces=faces,boundary_change=step_boundary_change,outlet_pressure=config%outlet_pressure, &
          inflow_state=inflow_state,inflow_temperature=inflow_t,inflow_mach=config%inflow_mach,inflow_pressure=mach_pressure, &
          wall_temperature=config%wall_temperature,wall_velocity=config%wall_velocity, &
          evaluate_reaction=evaluate_reaction,constant_coefficients=constant_coefficients,heat_release=sdc_heat_release)
      endif
    end subroutine
    subroutine evaluate_sources(u,t,stage_time,rates,good)
      real(dp), intent(in) :: u(:,:,:,:),t(:,:,:),stage_time
      real(dp), intent(out) :: rates(:,:,:,:,:)
      logical, intent(out) :: good
      integer :: direction
      real(dp), allocatable :: stage_ordinary(:,:,:,:)
      good=.false.;rates=0
      if(.not.ieee_is_finite(stage_time).or.size(rates,5)/=1+merge(1,0,present(evaluate_diffusion))) return
      call compute_reactive_external_source_3d(u,config%external_forcing,rates(:,:,:,:,1),good)
      if(.not.good) return
      if(linear_coefficient/=0) then
        good=.false.
        if(.not.all(ieee_is_finite(u(1:4,:,:,:)))) return
        if(any(u(1,:,:,:)<=0).or.maxval(abs(u(1:4,:,:,:)))>1e30_dp) return
        ! Pinned Forcing.cpp uses rho*(momentum-offset), not rho*(velocity-offset).
        ! It supplies momentum only; do not invent an accompanying energy source.
        do direction=1,3
          rates(direction+1,:,:,:,1)=rates(direction+1,:,:,:,1)+ &
            linear_coefficient*u(1,:,:,:)*(u(direction+1,:,:,:)-momentum_offset(direction))
        end do
        good=all(ieee_is_finite(rates))
        if(.not.good) return
      end if
      if(.not.present(evaluate_diffusion)) return
      if(allocated(frame)) then
        stage_ordinary=u
        call convert_regular_energy_3d(frame,domain_lower,[dx,dy,dz],stage_ordinary,0,.false.,good)
        if(.not.good) return
        call evaluate_diffusion(stage_ordinary,t,[dx,dy,dz],faces,inflow_state,inflow_t,rates(:,:,:,:,2), &
          new_boundary_rate,extra_diffusivity,good)
      else
        call evaluate_diffusion(u,t,[dx,dy,dz],faces,inflow_state,inflow_t,rates(:,:,:,:,2), &
          new_boundary_rate,extra_diffusivity,good)
      endif
      if(.not.good) return
      if(.not.all(ieee_is_finite(new_boundary_rate)).or..not.ieee_is_finite(extra_diffusivity)) then
        good=.false.;return
      end if
      source_calls=source_calls+1
      if(source_calls==1) old_boundary_rate=new_boundary_rate
      maximum_diffusivity=max(maximum_diffusivity,extra_diffusivity)
    end subroutine
  end subroutine run_reactive_3d_application

end module reactive_3d_application_mod
