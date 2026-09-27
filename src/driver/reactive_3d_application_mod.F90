module reactive_3d_application_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use constant_volume_reactor_mod, only: forced_reactor_integrator
  use constants_mod, only: pelef_version
  use state_indices_mod, only: irho, imx, imy, imz, iet, ncons
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use mesh_3d_mod, only: uniform_cell_centers_3d
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive, reactive_primitive_to_conserved
  use mixture_thermo_mod, only: mixture_density, mass_fractions_from_mole_fractions
  use selected_composition_mod, only: resolve_selected_composition
  use simulation_config_reactive_3d_mod, only: reactive_3d_config, resolved_reactive_boundaries_3d, &
    fill_mach_inlet_boundary_ghosts_3d,resolved_reactive_difmag_3d
  use reactive_entropy_wave_3d_problem_mod, only: &
    initialize_reactive_problem_3d, reactive_entropy_wave_density_3d
  use reactive_3d_mod, only: &
    compute_reactive_cfl_timestep_3d, advance_reactive_full_3d,advance_reactive_sdc_3d, &
    reactive_integrals_3d, reactive_extrema_3d,compute_reactive_external_source_3d,reactive_sdc_source_evaluator
  use reactive_transport_3d_mod, only: reactive_transport_timestep_3d, reactive_isothermal_timestep_3d
  use reactive_csv_io_3d_mod, only: write_reactive_3d_csv
  implicit none
  private

  public :: run_reactive_3d_application,read_regular_les_controls

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
  subroutine read_regular_les_controls(input_path,enabled,cs,prt,sct,ci,model,test_model,test_ratio,grid_ratio,harmonic_mean)
    character(len=*), intent(in) :: input_path
    logical, intent(out) :: enabled,harmonic_mean
    real(dp), intent(out) :: cs,prt,sct,ci
    character(len=*), intent(out) :: model,test_model
    integer, intent(out) :: test_ratio,grid_ratio
    integer :: unit,ios
    namelist /reactive_les/ enabled,cs,prt,sct,ci,model,test_model,test_ratio,grid_ratio,harmonic_mean
    enabled=.false.;cs=.17_dp;prt=.7_dp;sct=.7_dp;ci=0
    model='smagorinsky';test_model='box';test_ratio=2;grid_ratio=1;harmonic_mean=.false.
    open(newunit=unit,file=input_path,status='old',action='read',iostat=ios)
    if(ios/=0) error stop 'Cannot open LES input'
    read(unit,nml=reactive_les,iostat=ios)
    close(unit)
    if(ios>0) error stop 'Invalid reactive_les namelist'
  end subroutine

  subroutine run_reactive_3d_application( &
      input_path, application_label, bundle_sha256, config, species, &
      reactions, transport, base_mole_fractions, chemistry_integrator,evaluate_reaction,evaluate_diffusion)
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    procedure(application_diffusion), optional :: evaluate_diffusion
    real(dp), allocatable :: sdc_force(:),diffusion_rates(:,:,:,:),old_boundary_rate(:),new_boundary_rate(:)
    real(dp) :: extra_diffusivity
    integer :: source_calls
    logical :: les_enabled,les_harmonic
    real(dp) :: les_cs,les_prt,les_sct,les_ci
    integer :: les_test_ratio,les_grid_ratio
    character(len=32) :: les_model,les_filter
    character(len=*), intent(in) :: input_path, application_label
    character(len=*), intent(in) :: bundle_sha256
    type(reactive_3d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: base_mole_fractions(:)
    character(len=*), intent(in), optional :: chemistry_integrator

    type(gas_transport_species), allocatable :: sdc_transport(:)
    type(elementary_reaction), allocatable :: sdc_reactions(:)
    real(dp), allocatable :: sdc_reaction_source(:,:,:,:)
    real(dp), allocatable :: x(:), y(:), z(:)
    real(dp), allocatable :: state(:, :, :, :), temperature(:, :, :)
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

    if(.not.present(evaluate_diffusion)) then
      call read_regular_les_controls(input_path,les_enabled,les_cs,les_prt,les_sct,les_ci, &
        les_model,les_filter,les_test_ratio,les_grid_ratio,les_harmonic)
      if(les_enabled) error stop 'Regular LES backend is unavailable'
    end if
    if(present(evaluate_diffusion).and.config%time_integrator/='sdc') error stop 'Additional diffusion requires SDC'
    if((config%sdc_reactor=='cvode').neqv.present(evaluate_reaction)) &
      error stop 'SDC reactor selection does not match supplied backend'
    if (size(species) < 1 .or. size(reactions) < 1 .or. &
        size(transport) /= size(species) .or. &
        size(base_mole_fractions) /= size(species)) then
      error stop "Reactive 3D application data have incompatible dimensions"
    end if
    nvar = reactive_nvar(size(species))
    allocate(x(config%nx), y(config%ny), z(config%nz))
    allocate(state(nvar, config%nx, config%ny, config%nz))
    allocate(temperature(config%nx, config%ny, config%nz))
    allocate(initial_totals(nvar), final_totals(nvar), boundary_totals(nvar), step_boundary_change(nvar))
    if(config%time_integrator=='sdc'.and.config%transport_enabled) sdc_transport=transport
    if(config%time_integrator=='sdc'.and.config%chemistry_enabled) then
      sdc_reactions=reactions
      allocate(sdc_reaction_source(nvar-4,config%nx,config%ny,config%nz));sdc_reaction_source=0
    end if
    if(present(evaluate_diffusion)) then
      allocate(diffusion_rates, mold=state)
      allocate(old_boundary_rate(nvar),new_boundary_rate(nvar))
    else
      sdc_force=config%external_forcing
    end if
    boundary_totals=0.0_dp
    allocate(source_totals(nvar),step_source_change(nvar));source_totals=0
    if(any(config%external_forcing/=0)) external_force=config%external_forcing
    faces=resolved_reactive_boundaries_3d(config)
    allocate(mass_fractions(size(species)))
    allocate(primitive(reactive_nprim(size(species))))
    call uniform_cell_centers_3d( &
      config%nx, config%ny, config%nz, &
      config%x_lower, config%x_upper, &
      config%y_lower, config%y_upper, &
      config%z_lower, config%z_upper, x, y, z, dx, dy, dz)
    call initialize_reactive_problem_3d( &
      species, config, x, y, z, state, temperature, base_density, &
      mass_fractions, ok, base_mole_fractions)
    if (.not. ok) error stop "Failed to initialize reactive 3D problem"
    call reactive_integrals_3d( &
      state, config%nx, config%ny, config%nz, dx, dy, dz, &
      initial_totals, ok)
    if (.not. ok) error stop "Failed to integrate initial reactive 3D state"

    allocate(inflow_state(nvar,6),inflow_y(size(species)),inflow_x(size(species)))
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
        config%inflow_composition_species,config%inflow_composition_mole_fractions,species,inflow_x,ok,message)
      if (.not.ok) then
        write(*,'(a)') trim(message)
        error stop 'Invalid inflow composition'
      end if
      call mass_fractions_from_mole_fractions(species,inflow_x,inflow_y,ok)
      if (.not.ok) error stop 'Invalid inflow mass fractions'
    end if
    do face=1,6
      if (faces(face)/='inflow'.or.config%inflow_mach(face)>=0) cycle
      primitive=0.0_dp
      primitive(1)=mixture_density(species,inflow_y,config%inflow_pressure(face),config%inflow_temperature(face),ok)
      if (.not.ok) error stop 'Invalid inflow density'
      primitive(2:4)=config%inflow_velocity(:,face); primitive(5)=config%inflow_pressure(face)
      primitive(6:)=inflow_y
      call reactive_primitive_to_conserved(species,primitive,inflow_state(:,face),inflow_t(face),sound_speed,ok)
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
        call reactive_transport_timestep_3d(species,transport,inlet_cube,inlet_temperature,2,2,2,dx,dy,dz, &
          config%transport_cfl,config%viscosity_enabled,config%thermal_conduction_enabled, &
          config%species_diffusion_enabled,face_dt,face_diffusivity,ok)
        if (.not.ok) error stop 'Invalid inflow transport timestep'
        inlet_transport_dt=min(inlet_transport_dt,face_dt)
        inlet_diffusivity=max(inlet_diffusivity,face_diffusivity)
      end if
    end do

    write(*, '(a)') &
      "PeleF " // pelef_version // " " // trim(application_label)
    if (len_trim(bundle_sha256) > 0) then
      write(*, '(a,1x,a)') "Bundle SHA-256:", trim(bundle_sha256)
      if(config%time_integrator=='sdc'.and.config%chemistry_enabled) then
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
    if(config%time_integrator=='sdc'.and.len_trim(bundle_sha256)==0) &
      write(*,'(a,a)') 'SDC reactor: ',trim(config%sdc_reactor)
    if(config%time_integrator=='sdc') write(*,'(a,i0)') 'SDC iterations: ',config%sdc_iterations
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
    step = 0
    maximum_diffusivity = 0.0_dp
    minimum_transport_theta = 1.0_dp
    do while (time < config%final_time)
      if (step >= config%maximum_steps) then
        error stop "Maximum step count reached before reactive 3D final_time"
      end if
      call compute_reactive_cfl_timestep_3d( &
        species, state, temperature, config%nx, config%ny, config%nz, &
        dx, dy, dz, config%cfl, hydro_dt, ok)
      if (.not. ok) error stop "Failed to compute reactive 3D CFL timestep"
      if (any(config%inflow_mach>=0)) then
        call fill_mach_inlet_boundary_ghosts_3d(species,state,temperature,faces,config%inflow_mach, &
          inlet_ghost,inlet_ghost_t,ok,inflow_state=inflow_state,inflow_temperature=inflow_t, &
          inflow_pressure=mach_pressure,viscous_wall=.true.,wall_velocity=config%wall_velocity)
        if (.not.ok) error stop 'Failed to construct Mach inlet timestep states'
        call compute_reactive_cfl_timestep_3d(species,inlet_ghost,inlet_ghost_t, &
          config%nx+4,config%ny+4,config%nz+4,dx,dy,dz,config%cfl,inlet_dt,ok)
        if (.not.ok) error stop 'Failed to compute Mach inlet CFL timestep'
        if (config%transport_enabled) then
          call reactive_transport_timestep_3d(species,transport,inlet_ghost,inlet_ghost_t, &
            config%nx+4,config%ny+4,config%nz+4,dx,dy,dz,config%transport_cfl, &
            config%viscosity_enabled,config%thermal_conduction_enabled,config%species_diffusion_enabled, &
            inlet_transport_dt,inlet_diffusivity,ok)
          if (.not.ok) error stop 'Failed to compute Mach inlet transport timestep'
        end if
      end if
      dt = min(hydro_dt,inlet_dt)
      if (config%transport_enabled) then
        call reactive_transport_timestep_3d( &
          species, transport, state, temperature, &
          config%nx, config%ny, config%nz, dx, dy, dz, &
          config%transport_cfl, config%viscosity_enabled, &
          config%thermal_conduction_enabled, &
          config%species_diffusion_enabled, transport_dt, &
          step_maximum_diffusivity, ok)
        if (.not. ok) then
          error stop "Failed to compute reactive 3D transport step"
        end if
        if (any(faces=='isothermal_wall')) then
          call reactive_isothermal_timestep_3d(species,transport,state,temperature,[dx,dy,dz], &
            faces,config%wall_temperature,config%transport_cfl,wall_dt,ok)
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
        call evaluate_diffusion(state,temperature,[dx,dy,dz],faces,inflow_state,inflow_t, &
          diffusion_rates,new_boundary_rate,extra_diffusivity,ok)
        if(.not.ok) error stop 'Failed to compute additional diffusion timestep'
        if(.not.ieee_is_finite(extra_diffusivity).or.extra_diffusivity<0) error stop 'Invalid additional diffusivity'
        if(extra_diffusivity>0) dt=min(dt,.1_dp/(extra_diffusivity*sum(1._dp/[dx,dy,dz]**2)))
        maximum_diffusivity=max(maximum_diffusivity,extra_diffusivity)
      end if
      dt = min(dt, config%final_time - time)
      if(config%time_integrator=='sdc') then
        step_boundary_change=0;step_transport_theta=1;source_calls=0
        if(present(evaluate_diffusion)) then
          call perform_sdc(evaluate_sources)
        else
          call perform_sdc()
        end if
      else
      call advance_reactive_full_3d( &
        species, reactions, transport, state, temperature, &
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
        step_boundary_change=step_boundary_change+.5_dp*dt*(old_boundary_rate+new_boundary_rate)
        step_source_change=step_source_change-.5_dp*dt*(old_boundary_rate+new_boundary_rate)
      end if
      source_totals=source_totals+step_source_change
      boundary_totals = boundary_totals + step_boundary_change
      minimum_transport_theta = min( &
        minimum_transport_theta, step_transport_theta)
      time = time + dt
      step = step + 1
    end do

    call reactive_integrals_3d( &
      state, config%nx, config%ny, config%nz, dx, dy, dz, final_totals, ok)
    if (.not. ok) error stop "Failed to integrate final reactive 3D state"
    maximum_euler_conservation_error = maxval( &
      abs(final_totals(1:ncons) - initial_totals(1:ncons) - boundary_totals(1:ncons) - source_totals(1:ncons)) / &
        max(1.0_dp, abs(initial_totals(1:ncons))))
    maximum_all_component_conservation_error = maxval( &
      abs(final_totals - initial_totals - boundary_totals - source_totals) / max(1.0_dp, abs(initial_totals)))
    maximum_species_change = maxval(abs( &
      final_totals(6:nvar) - initial_totals(6:nvar)))
    call reactive_extrema_3d( &
      species, state, temperature, config%nx, config%ny, config%nz, &
      minimum_density, maximum_density, minimum_pressure, maximum_pressure, &
      minimum_temperature, maximum_temperature, maximum_speed, &
      maximum_closure_error, ok)
    if (.not. ok) error stop "Final reactive 3D state is not physical"

    if (trim(config%problem) == "entropy_wave" .and. all(faces=="periodic").and..not.allocated(external_force)) then
      density_l1_error = 0.0_dp
      do k = 1, config%nz
        do j = 1, config%ny
          do i = 1, config%nx
            call reactive_conserved_to_primitive( &
              species, state(:, i, j, k), temperature(i, j, k), primitive, &
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
    call write_reactive_3d_csv( &
      trim(config%output_file), species, x, y, z, state, temperature, &
      config%nx, config%ny, config%nz, time, ok, message)
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
    subroutine perform_sdc(source_callback)
      procedure(reactive_sdc_source_evaluator), optional :: source_callback
        call advance_reactive_sdc_3d(species,state,temperature,[dx,dy,dz],dt,time, &
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
          evaluate_reaction=evaluate_reaction)
    end subroutine
    subroutine evaluate_sources(u,t,stage_time,rates,good)
      real(dp), intent(in) :: u(:,:,:,:),t(:,:,:),stage_time
      real(dp), intent(out) :: rates(:,:,:,:,:)
      logical, intent(out) :: good
      good=.false.;rates=0
      if(.not.ieee_is_finite(stage_time).or.size(rates,5)/=2) return
      call compute_reactive_external_source_3d(u,config%external_forcing,rates(:,:,:,:,1),good)
      if(.not.good) return
      call evaluate_diffusion(u,t,[dx,dy,dz],faces,inflow_state,inflow_t,rates(:,:,:,:,2), &
        new_boundary_rate,extra_diffusivity,good)
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
