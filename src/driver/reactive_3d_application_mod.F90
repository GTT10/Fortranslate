module reactive_3d_application_mod
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
  use state_indices_mod, only: irho, imx, imy, imz, iet, ncons
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use mesh_3d_mod, only: uniform_cell_centers_3d
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use simulation_config_reactive_3d_mod, only: reactive_3d_config
  use reactive_entropy_wave_3d_problem_mod, only: &
    initialize_reactive_problem_3d, reactive_entropy_wave_density_3d
  use reactive_3d_mod, only: &
    compute_reactive_cfl_timestep_3d, advance_reactive_full_3d, &
    reactive_integrals_3d, reactive_extrema_3d
  use reactive_transport_3d_mod, only: reactive_transport_timestep_3d
  use reactive_csv_io_3d_mod, only: write_reactive_3d_csv
  implicit none
  private

  public :: run_reactive_3d_application

contains

  subroutine run_reactive_3d_application( &
      input_path, application_label, bundle_sha256, config, species, &
      reactions, transport, base_mole_fractions, chemistry_integrator)
    character(len=*), intent(in) :: input_path, application_label
    character(len=*), intent(in) :: bundle_sha256
    type(reactive_3d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: base_mole_fractions(:)
    character(len=*), intent(in), optional :: chemistry_integrator

    real(dp), allocatable :: x(:), y(:), z(:)
    real(dp), allocatable :: state(:, :, :, :), temperature(:, :, :)
    real(dp), allocatable :: initial_totals(:), final_totals(:)
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
    integer :: step, i, j, k, nvar

    if (size(species) < 1 .or. size(reactions) < 1 .or. &
        size(transport) /= size(species) .or. &
        size(base_mole_fractions) /= size(species)) then
      error stop "Reactive 3D application data have incompatible dimensions"
    end if
    nvar = reactive_nvar(size(species))
    allocate(x(config%nx), y(config%ny), z(config%nz))
    allocate(state(nvar, config%nx, config%ny, config%nz))
    allocate(temperature(config%nx, config%ny, config%nz))
    allocate(initial_totals(nvar), final_totals(nvar))
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

    write(*, '(a)') &
      "PeleF " // pelef_version // " " // trim(application_label)
    if (len_trim(bundle_sha256) > 0) then
      write(*, '(a,1x,a)') "Bundle SHA-256:", trim(bundle_sha256)
      if (present(chemistry_integrator)) then
        write(*, '(a,1x,a)') "Chemistry integrator:", &
          trim(chemistry_integrator)
      end if
    end if
    write(*, '(a,1x,a)') "Input:", trim(input_path)
    write(*, '(a,i0,a,i0,a,i0)') &
      "Grid: nx=", config%nx, ", ny=", config%ny, ", nz=", config%nz
    write(*, '(a,1x,a)') "Problem:", trim(config%problem)
    write(*, '(a,1x,a)') "Thermodynamics:", trim(config%thermo_model)
    write(*, '(a,1x,a)') "Reconstruction:", trim(config%reconstruction)
    write(*, '(a,1x,a)') "Limiter:", trim(config%limiter)
    write(*, '(a,1x,a)') "Riemann solver:", trim(config%riemann_solver)
    write(*, '(a,1x,a)') &
      "Boundary condition:", trim(config%boundary_condition)
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
      dt = hydro_dt
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
        dt = min(dt, transport_dt)
        maximum_diffusivity = max( &
          maximum_diffusivity, step_maximum_diffusivity)
      end if
      dt = min(dt, config%final_time - time)
      call advance_reactive_full_3d( &
        species, reactions, transport, state, temperature, &
        config%nx, config%ny, config%nz, dx, dy, dz, dt, &
        config%riemann_solver, config%chemistry_enabled, &
        config%chemistry_relative_tolerance, &
        config%chemistry_absolute_tolerance, config%transport_enabled, &
        config%viscosity_enabled, config%thermal_conduction_enabled, &
        config%species_diffusion_enabled, config%barodiffusion_enabled, &
        step_transport_theta, ok, config%reconstruction, config%limiter, &
        chemistry_integrator=chemistry_integrator)
      if (.not. ok) then
        error stop "Reactive 3D split update rejected its candidate"
      end if
      minimum_transport_theta = min( &
        minimum_transport_theta, step_transport_theta)
      time = time + dt
      step = step + 1
    end do

    call reactive_integrals_3d( &
      state, config%nx, config%ny, config%nz, dx, dy, dz, final_totals, ok)
    if (.not. ok) error stop "Failed to integrate final reactive 3D state"
    maximum_euler_conservation_error = maxval( &
      abs(final_totals(1:ncons) - initial_totals(1:ncons)) / &
        max(1.0_dp, abs(initial_totals(1:ncons))))
    maximum_all_component_conservation_error = maxval( &
      abs(final_totals - initial_totals) / max(1.0_dp, abs(initial_totals)))
    maximum_species_change = maxval(abs( &
      final_totals(6:nvar) - initial_totals(6:nvar)))
    call reactive_extrema_3d( &
      species, state, temperature, config%nx, config%ny, config%nz, &
      minimum_density, maximum_density, minimum_pressure, maximum_pressure, &
      minimum_temperature, maximum_temperature, maximum_speed, &
      maximum_closure_error, ok)
    if (.not. ok) error stop "Final reactive 3D state is not physical"

    if (trim(config%problem) == "entropy_wave") then
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
      error stop "Reactive 3D periodic conservation gate failed"
    end if
    if (.not. config%chemistry_enabled .and. &
        maximum_all_component_conservation_error > 5.0e-11_dp) then
      error stop "Inert reactive 3D species conservation gate failed"
    end if
    if (maximum_closure_error > 5.0e-11_dp) then
      error stop "Reactive 3D species-closure gate failed"
    end if
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
    if (trim(config%problem) == "entropy_wave") then
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
  end subroutine run_reactive_3d_application

end module reactive_3d_application_mod
