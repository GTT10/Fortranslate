module reactive_eb_3d_application_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
  use state_indices_mod, only: ncons, imx, imy, imz
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d, eb_cut_cell_3d, &
    eb_regular_cell_3d
  use reactive_eb_cfl_3d_mod, only: &
    compute_reactive_eb_cfl_timestep_3d
  use eb_reactive_transport_3d_mod, only: &
    reactive_eb_transport_timestep_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config,build_configured_eb_geometry_3d
  use reactive_eb_3d_checkpoint_mod, only: &
    write_reactive_eb_3d_checkpoint, read_reactive_eb_3d_checkpoint
  use reactive_eb_3d_driver_mod, only: &
    initialize_reactive_eb_density_sheet_3d, &
    advance_reactive_eb_full_3d, &
    reactive_eb_integrals_3d, reactive_eb_extrema_3d, &
    reactive_eb_element_integrals_3d, reactive_eb_composition_integrals_3d, &
    reactive_eb_element_species_supported_3d, reactive_eb_element_amounts_3d, &
    write_reactive_eb_3d_csv
  implicit none
  private

  public :: run_reactive_eb_3d_application

contains

  subroutine run_reactive_eb_3d_application( &
      input_path, application_label, bundle_sha256, config, species, &
      reactions, transport, base_mole_fractions, chemistry_integrator, element_names, element_atoms)
    character(len=*), intent(in) :: input_path, application_label
    character(len=*), intent(in) :: bundle_sha256
    type(reactive_eb_3d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in), optional :: base_mole_fractions(:)
    character(len=*), intent(in), optional :: chemistry_integrator
    character(len=*), intent(in), optional :: element_names(:)
    real(dp), intent(in), optional :: element_atoms(:,:)

    type(eb_geometry_3d) :: geometry
    real(dp), allocatable :: state(:, :, :, :), new_state(:, :, :, :)
    real(dp), allocatable :: temperature(:, :, :), new_temperature(:, :, :)
    real(dp), allocatable :: initial_integrals(:), final_integrals(:)
    real(dp), allocatable :: initial_l1_integrals(:), final_l1_integrals(:)
    real(dp), allocatable :: initial_element_integrals(:),final_element_integrals(:),element_roundoff(:)
    real(dp), allocatable :: boundary(:),wall(:),chemical(:),step_boundary(:),step_wall(:),step_chemical(:)
    real(dp), allocatable :: external_elements(:)
    real(dp) :: element_residual,component_roundoff
    real(dp) :: time, dt, hydro_dt, transport_dt, minimum_dt
    real(dp) :: conservation_error
    real(dp) :: component_error, normal_momentum_change
    real(dp) :: element_conservation_error
    real(dp) :: maximum_transport_diffusivity, step_transport_diffusivity
    real(dp) :: minimum_transport_theta, step_transport_theta
    real(dp) :: minimum_density, maximum_density
    real(dp) :: minimum_pressure, maximum_pressure
    real(dp) :: minimum_temperature, maximum_temperature
    real(dp) :: maximum_speed, maximum_closure_error
    character(len=1024) :: message
    logical :: ok, restarted, stopped_after_checkpoint, selected_context
    integer :: step, nvar, component, normal_momentum,element_count,e

    if (size(species) < 2 .or. size(reactions) < 1 .or. &
        size(transport) /= size(species)) then
      error stop "Reactive EB 3D application data have incompatible dimensions"
    end if
    selected_context = len_trim(bundle_sha256) > 0
    if (selected_context) then
      if (.not. (present(base_mole_fractions) .and. &
          present(chemistry_integrator))) then
        error stop "Reactive EB 3D selected context is incomplete"
      end if
    else if (present(base_mole_fractions) .or. &
             present(chemistry_integrator)) then
      error stop "Reactive EB 3D fixed runtime received selected context"
    end if
    if (present(base_mole_fractions)) then
      if (size(base_mole_fractions) /= size(species)) then
        error stop "Reactive EB 3D composition has incompatible dimensions"
      end if
    end if
    if (present(element_names).neqv.present(element_atoms)) error stop 'Incomplete EB element metadata'
    element_count=3
    if (present(element_names)) then
      if (.not.selected_context) error stop 'Element metadata requires selected EB context'
      element_count=size(element_names)
      if (element_count<1) error stop 'Empty EB element metadata'
      if (any(shape(element_atoms)/=[element_count,size(species)])) error stop 'Invalid EB element metadata shape'
    else if (config%chemistry_enabled) then
      if (.not.reactive_eb_element_species_supported_3d(species)) &
        error stop "Reactive EB 3D chemistry requires supported H/O/N species"
    end if
    allocate(initial_element_integrals(element_count),final_element_integrals(element_count),element_roundoff(element_count))

    call build_configured_eb_geometry_3d(config,geometry,ok)
    if (.not. ok) error stop "Failed to build public 3D EB geometry"

    nvar = reactive_nvar(size(species))
    allocate(state(nvar, config%nx, config%ny, config%nz))
    allocate(new_state(nvar, config%nx, config%ny, config%nz))
    allocate(temperature(config%nx, config%ny, config%nz))
    allocate(new_temperature(config%nx, config%ny, config%nz))
    allocate(initial_integrals(nvar), final_integrals(nvar))
    allocate(boundary(nvar),wall(nvar),chemical(nvar),step_boundary(nvar),step_wall(nvar),step_chemical(nvar))
    allocate(external_elements(element_count))
    boundary=0._dp;wall=0._dp;chemical=0._dp
    allocate(initial_l1_integrals(nvar), final_l1_integrals(nvar))
    call initialize_reactive_eb_density_sheet_3d( &
      species, config, geometry, state, temperature, ok, base_mole_fractions)
    if (.not. ok) error stop "Failed to initialize public 3D EB state"
    call reactive_eb_integrals_3d( &
      state, geometry, initial_integrals, ok, initial_l1_integrals)
    if (.not. ok) error stop "Failed to integrate initial 3D EB state"
    initial_element_integrals = 0.0_dp
    if (config%chemistry_enabled.or.present(element_atoms)) then
      if (present(element_atoms)) then
        call reactive_eb_composition_integrals_3d(species,state,geometry,element_atoms,initial_element_integrals,ok)
      else
        call reactive_eb_element_integrals_3d(species,state,geometry,initial_element_integrals,ok)
      end if
      if (.not. ok) error stop "Failed to integrate initial EB elements"
    end if

    time = 0.0_dp
    step = 0
    minimum_dt = huge(1.0_dp)
    maximum_transport_diffusivity = 0.0_dp
    minimum_transport_theta = 1.0_dp
    restarted = .false.
    stopped_after_checkpoint = .false.
    if (len_trim(config%restart_file) > 0) then
      if (selected_context) then
        call read_reactive_eb_3d_checkpoint( &
          trim(config%restart_file), species, reactions, transport, config, &
          geometry, state, temperature, time, step, initial_integrals, &
          initial_l1_integrals, initial_element_integrals, minimum_dt, &
          maximum_transport_diffusivity, minimum_transport_theta, ok, &
          message, bundle_sha256=bundle_sha256, &
          chemistry_integrator=chemistry_integrator, &
          base_mole_fractions=base_mole_fractions,element_names=element_names,element_atoms=element_atoms, &
          boundary_change=boundary,wall_change=wall,chemistry_change=chemical)
      else
        call read_reactive_eb_3d_checkpoint( &
          trim(config%restart_file), species, reactions, transport, config, &
          geometry, state, temperature, time, step, initial_integrals, &
          initial_l1_integrals, initial_element_integrals, minimum_dt, &
          maximum_transport_diffusivity, minimum_transport_theta, ok, &
          message,boundary_change=boundary,wall_change=wall,chemistry_change=chemical)
      end if
      if (.not. ok) then
        write(*, '(a)') trim(message)
        error stop "Failed to restart public 3D EB state"
      end if
      call check_inventory()
      restarted = .true.
      write(*, '(a,1x,a)') &
        "Restarted from checkpoint:", trim(config%restart_file)
      write(*, '(a,i0,a,es24.16)') &
        "Restored step ", step, ", time ", time
    end if
    do while (time < config%final_time)
      if (step >= config%maximum_steps) then
        error stop "Maximum step count reached before 3D EB final_time"
      end if
      call compute_reactive_eb_cfl_timestep_3d( &
        species, state, temperature, geometry, config%cfl, hydro_dt, ok)
      if (.not. ok) error stop "Failed to compute public 3D EB timestep"
      dt = hydro_dt
      if (config%transport_enabled) then
        call reactive_eb_transport_timestep_3d( &
          species, transport, state, temperature, geometry, &
          config%transport_cfl, config%viscosity_enabled, &
          config%thermal_conduction_enabled, &
          config%species_diffusion_enabled, transport_dt, &
          step_transport_diffusivity, ok)
        if (.not. ok) then
          error stop "Failed to compute public 3D EB transport step"
        end if
        dt = min(dt, transport_dt)
        maximum_transport_diffusivity = max( &
          maximum_transport_diffusivity, step_transport_diffusivity)
      end if
      dt = min(dt, config%final_time - time)
      call advance_reactive_eb_full_3d( &
        species, reactions, transport, state, temperature, geometry, &
        config%riemann_solver, config%redistribution, dt, &
        config%chemistry_enabled, config%chemistry_relative_tolerance, &
        config%chemistry_absolute_tolerance, config%transport_enabled, &
        config%viscosity_enabled, config%thermal_conduction_enabled, &
        config%species_diffusion_enabled, config%barodiffusion_enabled, &
        new_state, new_temperature, step_transport_theta, ok, &
        config%state_redist_target_volume_fraction, chemistry_integrator, &
        step_boundary,step_wall,step_chemical,difmag=config%difmag)
      if (.not. ok) error stop "Public 3D EB update rejected its candidate"
      boundary=boundary+step_boundary
      wall=wall+step_wall
      chemical=chemical+step_chemical
      minimum_transport_theta = min( &
        minimum_transport_theta, step_transport_theta)
      state = new_state
      temperature = new_temperature
      time = time + dt
      step = step + 1
      minimum_dt = min(minimum_dt, dt)
      if (config%checkpoint_interval_steps > 0) then
        if (mod(step, config%checkpoint_interval_steps) == 0) then
          call check_inventory()
          if (selected_context) then
            call write_reactive_eb_3d_checkpoint( &
              trim(config%checkpoint_file), species, reactions, transport, &
              config, geometry, state, temperature, time, step, &
              initial_integrals, initial_l1_integrals, &
              initial_element_integrals, minimum_dt, &
              maximum_transport_diffusivity, minimum_transport_theta, &
              ok, message, bundle_sha256=bundle_sha256, &
              chemistry_integrator=chemistry_integrator, &
              base_mole_fractions=base_mole_fractions,element_names=element_names,element_atoms=element_atoms, &
          boundary_change=boundary,wall_change=wall,chemistry_change=chemical)
          else
            call write_reactive_eb_3d_checkpoint( &
              trim(config%checkpoint_file), species, reactions, transport, &
              config, geometry, state, temperature, time, step, &
              initial_integrals, initial_l1_integrals, &
              initial_element_integrals, minimum_dt, &
              maximum_transport_diffusivity, minimum_transport_theta, &
              ok, message,boundary_change=boundary,wall_change=wall,chemistry_change=chemical)
          end if
          if (.not. ok) then
            write(*, '(a)') trim(message)
            error stop "Failed to write public 3D EB checkpoint"
          end if
          write(*, '(a,1x,a,a,i0,a,es24.16)') &
            "Wrote checkpoint:", trim(config%checkpoint_file), &
            ", step ", step, ", time ", time
          if (config%stop_after_checkpoint) then
            stopped_after_checkpoint = .true.
            exit
          end if
        end if
      end if
    end do

    call check_inventory()
    do component=1,nvar
      write(*,'(a,i0,6(1x,es24.16))') "EB inventory initial/final/boundary/wall/chemical/residual: ",component, &
        initial_integrals(component),final_integrals(component),boundary(component),wall(component),chemical(component), &
        final_integrals(component)-initial_integrals(component)-boundary(component)-wall(component)-chemical(component)
    end do
    call write_reactive_eb_3d_csv( &
      trim(config%output_file), species, geometry, state, temperature, &
      time, ok, message)
    if (.not. ok) then
      write(*, '(a)') trim(message)
      error stop 3
    end if

    write(*, '(a)') &
      "PeleF " // pelef_version // " " // trim(application_label)
    if (len_trim(bundle_sha256) > 0) then
      write(*, '(a,1x,a)') "Bundle SHA-256:", trim(bundle_sha256)
      if (present(chemistry_integrator)) then
        write(*, '(a,1x,a)') &
          "Chemistry integrator:", trim(chemistry_integrator)
      end if
    end if
    write(*, '(a,1x,a)') "Input:", trim(input_path)
    write(*, '(a,i0,a,i0,a,i0)') &
      "Grid: nx=", config%nx, ", ny=", config%ny, ", nz=", config%nz
    write(*,'(a,a)') 'Geometry: ',trim(config%geometry_type)
    if (config%geometry_type=='extruded') then
      write(*,'(a,a)') 'Extrusion axis: ',trim(config%extrusion_axis)
    end if
    if (config%geometry_type/='plane') then
      write(*,'(a,3es24.16)') 'Momentum changes: ',final_integrals(imx:imz)-initial_integrals(imx:imz)
    end if
    if (config%geometry_type=='plane') then
    write(*, '(a,1x,a)') "Plane axis:", trim(config%plane_axis)
    write(*, '(a,es24.16)') "Plane position: ", config%plane_position
    end if
    write(*, '(a,1x,a)') "Riemann solver:", trim(config%riemann_solver)
    write(*, '(a,1x,a)') "Redistribution:", trim(config%redistribution)
    write(*, '(a,l2)') "Chemistry: ", config%chemistry_enabled
    write(*, '(a,es24.16)') "Chemistry relative tolerance: ", &
      config%chemistry_relative_tolerance
    write(*, '(a,es24.16)') "Chemistry absolute tolerance: ", &
      config%chemistry_absolute_tolerance
    write(*, '(a,l2)') "Molecular transport: ", config%transport_enabled
    write(*, '(a,l2)') "Viscosity: ", config%viscosity_enabled
    write(*, '(a,l2)') "Thermal conduction: ", &
      config%thermal_conduction_enabled
    write(*, '(a,l2)') "Species diffusion: ", &
      config%species_diffusion_enabled
    write(*, '(a,l2)') "Barodiffusion: ", config%barodiffusion_enabled
    if (config%chemistry_enabled .and. config%transport_enabled) then
      write(*, '(a)') "Operator sequence: R-T-H-T-R"
    else if (config%chemistry_enabled) then
      write(*, '(a)') "Operator sequence: R-H-R"
    else if (config%transport_enabled) then
      write(*, '(a)') "Operator sequence: T-H-T"
    else
      write(*, '(a)') "Operator sequence: H"
    end if
    write(*, '(a,es24.16)') "Transport CFL: ", config%transport_cfl
    write(*, '(a,es24.16)') "StateRedist target: ", &
      config%state_redist_target_volume_fraction
    write(*, '(a,i0)') "Regular cells: ", &
      count(geometry%cell_type == eb_regular_cell_3d)
    write(*, '(a,i0)') "Cut cells: ", &
      count(geometry%cell_type == eb_cut_cell_3d)
    write(*, '(a,i0)') "Covered cells: ", &
      count(geometry%cell_type == eb_covered_cell_3d)
    write(*, '(a,i0)') "Completed steps: ", step
    write(*, '(a,es24.16)') "Final time: ", time
    if (restarted) write(*, '(a)') "Restart continuation: complete"
    if (stopped_after_checkpoint) write(*, '(a)') &
      "Stopped after checkpoint"
    write(*, '(a,es24.16)') "Minimum accepted dt: ", minimum_dt
    write(*, '(a,es24.16)') "Maximum transport diffusivity: ", &
      maximum_transport_diffusivity
    write(*, '(a,es24.16)') "Minimum transport flux theta: ", &
      minimum_transport_theta
    write(*, '(a,es24.16)') "Maximum invariant conservation error: ", &
      conservation_error
    write(*, '(a,es24.16)') "Maximum elemental conservation error: ", &
      element_conservation_error
    if (config%geometry_type=='plane') write(*, '(a,es24.16)') "Normal momentum change: ", &
      normal_momentum_change
    write(*, '(a,es24.16)') "Minimum density: ", minimum_density
    write(*, '(a,es24.16)') "Maximum density: ", maximum_density
    write(*, '(a,es24.16)') "Minimum pressure: ", minimum_pressure
    write(*, '(a,es24.16)') "Maximum pressure: ", maximum_pressure
    write(*, '(a,es24.16)') "Minimum temperature: ", minimum_temperature
    write(*, '(a,es24.16)') "Maximum temperature: ", maximum_temperature
    write(*, '(a,es24.16)') "Maximum speed: ", maximum_speed
    write(*, '(a,es24.16)') "Maximum species closure error: ", &
      maximum_closure_error
    write(*, '(a,1x,a)') "Output:", trim(config%output_file)
  contains

    subroutine check_inventory()
      if (.not.all(ieee_is_finite(boundary)).or..not.all(ieee_is_finite(wall)).or. &
          .not.all(ieee_is_finite(chemical))) error stop "Nonfinite cumulative EB exchanges"
      call reactive_eb_integrals_3d( &
        state, geometry, final_integrals, ok, final_l1_integrals)
      if (.not. ok) error stop "Failed to integrate final 3D EB state"
      element_conservation_error = 0.0_dp
      if (config%chemistry_enabled) then
        call reactive_eb_element_amounts_3d(species,boundary+wall,external_elements,ok,element_atoms)
        if (.not.ok) error stop "Failed to project external EB exchanges"
        if (present(element_atoms)) then
          call reactive_eb_composition_integrals_3d(species,state,geometry,element_atoms,final_element_integrals,ok)
          if (.not.ok) error stop 'Failed to integrate final EB elements'
          ! Preserve the relative gate and separately account for roundoff of
          ! initially absent elements from species mass closure and summation.
          do e=1,element_count
            element_roundoff(e)=64._dp*epsilon(1._dp)*real(max(1,step),dp)*sum(abs(initial_integrals(6:)))* &
              maxval(element_atoms(e,:)/species%molecular_weight)
            element_residual=max(0._dp,abs(final_element_integrals(e)-initial_element_integrals(e)-external_elements(e)) &
              -element_roundoff(e))
            element_conservation_error=max(element_conservation_error, &
              element_residual/max(1.e-30_dp,abs(initial_element_integrals(e))))
          end do
          if (.not.ieee_is_finite(element_conservation_error)) error stop 'Nonfinite elemental conservation error'
        else
          call reactive_eb_element_integrals_3d(species,state,geometry,final_element_integrals,ok)
          if (.not.ok) error stop 'Failed to integrate final EB elements'
          element_conservation_error=maxval(abs(final_element_integrals-initial_element_integrals-external_elements)/ &
            max(1.e-30_dp,abs(initial_element_integrals)))
        end if
      end if
      normal_momentum=imx
      if (config%geometry_type=='plane') then
      select case (trim(config%plane_axis))
      case ("x")
        normal_momentum = imx
      case ("y")
        normal_momentum = imy
      case ("z")
        normal_momentum = imz
      case default
        error stop "Unexpected public 3D EB plane axis"
      end select
      end if
      conservation_error = 0.0_dp
      do component = 1, nvar
        component_roundoff=0._dp
        if (component>ncons) component_roundoff=64._dp*epsilon(1._dp)*real(max(1,step),dp)* &
          sum(abs(initial_integrals(6:)))
        component_error = max(0._dp,abs( &
          final_integrals(component)-initial_integrals(component)-boundary(component)-wall(component)-chemical(component)) &
          -component_roundoff) / &
          max(tiny(1.0_dp), abs(initial_integrals(component)), &
            initial_l1_integrals(component), final_l1_integrals(component), &
            abs(boundary(component)),abs(wall(component)),abs(chemical(component)))
        if (component_error>5.e-11_dp) write(*,'(a,i0,3(1x,es24.16))') &
          'Invariant component initial/final/error: ',component, &
          initial_integrals(component),final_integrals(component),component_error
        conservation_error = max(conservation_error, component_error)
      end do
      if (.not.ieee_is_finite(conservation_error)) error stop "Nonfinite EB conservation error"
      normal_momentum_change = &
        final_integrals(normal_momentum) - initial_integrals(normal_momentum)
      call reactive_eb_extrema_3d( &
        species, state, temperature, geometry, minimum_density, &
        maximum_density, minimum_pressure, maximum_pressure, &
        minimum_temperature, maximum_temperature, maximum_speed, &
        maximum_closure_error, ok)
      if (.not. ok) error stop "Final public 3D EB state is not physical"
      if (conservation_error > 5.0e-11_dp) then
        write(*, '(a,es24.16)') &
          "Rejected invariant conservation error: ", conservation_error
        error stop "Public 3D EB invariant-conservation gate failed"
      end if
      if (config%chemistry_enabled .and. &
          element_conservation_error > 5.0e-10_dp) then
        write(*, '(a,es24.16)') &
          "Rejected elemental conservation error: ", element_conservation_error
        error stop "Public 3D EB elemental-conservation gate failed"
      end if
      if (maximum_closure_error > 5.0e-11_dp) then
        error stop "Public 3D EB species-closure gate failed"
      end if
    end subroutine check_inventory

  end subroutine run_reactive_eb_3d_application

end module reactive_eb_3d_application_mod
