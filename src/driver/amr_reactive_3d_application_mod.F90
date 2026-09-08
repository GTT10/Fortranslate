module amr_reactive_3d_application_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
  use state_indices_mod, only: ncons
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: &
    gas_transport_species, compatible_transport_database
  use mesh_3d_mod, only: uniform_cell_centers_3d
  use reactive_1d_mod, only: reactive_nvar
  use simulation_config_reactive_3d_mod, only: reactive_3d_config
  use simulation_config_amr_reactive_3d_mod, only: amr_reactive_3d_config
  use reactive_entropy_wave_3d_problem_mod, only: &
    initialize_reactive_problem_3d
  use reactive_3d_mod, only: &
    recover_reactive_temperatures_3d, reactive_extrema_3d
  use amr_hierarchy_3d_mod, only: &
    amr_patch_3d, initialize_amr_patch_3d, restrict_average_3d, &
    average_down_3d, composite_integrals_amr_3d
  use amr_reactive_3d_mod, only: &
    compute_amr_reactive_cfl_timestep_3d, &
    composite_element_integrals_amr_3d, &
    element_integrals_from_reactive_integrals_3d
  use amr_reactive_transport_3d_mod, only: &
    compute_amr_reactive_transport_timestep_3d, &
    advance_amr_reactive_full_3d
  use amr_reactive_3d_checkpoint_mod, only: &
    write_amr_reactive_3d_checkpoint, read_amr_reactive_3d_checkpoint
  use reactive_csv_io_3d_mod, only: write_reactive_3d_csv
  implicit none
  private

  public :: run_amr_reactive_3d_application

contains

  subroutine run_amr_reactive_3d_application( &
      input_path, application_label, bundle_sha256, config, amr_config, &
      species, reactions, transport, base_mole_fractions, &
      chemistry_integrator)
    character(len=*), intent(in) :: input_path, application_label
    character(len=*), intent(in) :: bundle_sha256
    type(reactive_3d_config), intent(in) :: config
    type(amr_reactive_3d_config), intent(in) :: amr_config
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: base_mole_fractions(:)
    character(len=*), intent(in), optional :: chemistry_integrator

    type(reactive_3d_config) :: fine_config
    type(amr_patch_3d) :: patch
    real(dp), allocatable :: coarse_state(:, :, :, :)
    real(dp), allocatable :: fine_state(:, :, :, :)
    real(dp), allocatable :: coarse_temperature(:, :, :)
    real(dp), allocatable :: fine_temperature(:, :, :)
    real(dp), allocatable :: recovered_temperature(:, :, :)
    real(dp), allocatable :: restricted(:, :, :, :)
    real(dp), allocatable :: initial_integrals(:), final_integrals(:)
    real(dp), allocatable :: mass_fractions(:)
    real(dp), allocatable :: x(:), y(:), z(:), xf(:), yf(:), zf(:)
    real(dp) :: initial_elements(3), final_elements(3)
    real(dp) :: dx, dy, dz, time, dt, hydro_dt, transport_dt
    real(dp) :: base_density, fine_base_density
    real(dp) :: step_reflux, maximum_reflux, conservation_error
    real(dp) :: step_transport_theta, minimum_transport_theta
    real(dp) :: step_maximum_diffusivity, maximum_diffusivity
    real(dp) :: elemental_conservation_error, species_integral_change
    real(dp) :: synchronization_error, coarse_closure, fine_closure
    character(len=1024) :: message
    logical :: ok, restarted, stopped_after_checkpoint
    logical :: selected_context, element_diagnostics
    integer :: step, nvar, covered_nx, covered_ny, covered_nz

    selected_context = len_trim(bundle_sha256) > 0
    if (size(species) < 1 .or. size(reactions) < 1 .or. &
        size(base_mole_fractions) /= size(species)) then
      error stop "Static 3D AMR application data have incompatible dimensions"
    end if
    if (.not. compatible_transport_database(species, transport)) then
      error stop "Static 3D AMR application data have incompatible dimensions"
    end if
    if (.not. all(ieee_is_finite(base_mole_fractions))) then
      error stop "Static 3D AMR composition is invalid"
    end if
    if (minval(base_mole_fractions) < 0.0_dp .or. &
        abs(sum(base_mole_fractions) - 1.0_dp) > 5.0e-10_dp) then
      error stop "Static 3D AMR composition is invalid"
    end if
    if (selected_context) then
      if (len_trim(bundle_sha256) /= 64 .or. &
          .not. present(chemistry_integrator)) then
        error stop "Selected static 3D AMR context is incomplete"
      end if
    end if

    call initialize_amr_patch_3d( &
      config%nx, config%ny, config%nz, &
      amr_config%coarse_i_lower, amr_config%coarse_i_upper, &
      amr_config%coarse_j_lower, amr_config%coarse_j_upper, &
      amr_config%coarse_k_lower, amr_config%coarse_k_upper, &
      amr_config%refinement_ratio, patch, ok)
    if (.not. ok .or. .not. patch%is_strictly_interior()) then
      error stop "Failed to construct strictly interior 3D AMR patch"
    end if

    nvar = reactive_nvar(size(species))
    allocate(x(config%nx), y(config%ny), z(config%nz))
    allocate(xf(patch%fine_nx()), yf(patch%fine_ny()), zf(patch%fine_nz()))
    allocate(coarse_state(nvar, config%nx, config%ny, config%nz))
    allocate(coarse_temperature(config%nx, config%ny, config%nz))
    allocate(fine_state(nvar, patch%fine_nx(), patch%fine_ny(), &
      patch%fine_nz()))
    allocate(fine_temperature(patch%fine_nx(), patch%fine_ny(), &
      patch%fine_nz()))
    allocate(mass_fractions(size(species)))
    allocate(initial_integrals(nvar), final_integrals(nvar))
    call uniform_cell_centers_3d( &
      config%nx, config%ny, config%nz, config%x_lower, config%x_upper, &
      config%y_lower, config%y_upper, config%z_lower, config%z_upper, &
      x, y, z, dx, dy, dz)
    call fine_patch_centers(patch, config, dx, dy, dz, xf, yf, zf)
    allocate(recovered_temperature(config%nx, config%ny, config%nz))
    time = 0.0_dp
    step = 0
    maximum_reflux = 0.0_dp
    maximum_diffusivity = 0.0_dp
    minimum_transport_theta = 1.0_dp
    restarted = len_trim(amr_config%restart_file) > 0
    if (restarted) then
      if (selected_context) then
        if (config%transport_enabled) then
          call read_amr_reactive_3d_checkpoint( &
            trim(amr_config%restart_file), species, config, amr_config, &
            patch, coarse_state, coarse_temperature, fine_state, &
            fine_temperature, time, step, initial_integrals, maximum_reflux, &
            ok, message, bundle_sha256=bundle_sha256, &
            chemistry_integrator=chemistry_integrator, &
            base_mole_fractions=base_mole_fractions, reactions=reactions, &
            transport=transport, &
            maximum_transport_diffusivity=maximum_diffusivity, &
            minimum_transport_theta=minimum_transport_theta)
        else
          call read_amr_reactive_3d_checkpoint( &
            trim(amr_config%restart_file), species, config, amr_config, patch, &
            coarse_state, coarse_temperature, fine_state, fine_temperature, &
            time, step, initial_integrals, maximum_reflux, ok, message, &
            bundle_sha256=bundle_sha256, &
            chemistry_integrator=chemistry_integrator, &
            base_mole_fractions=base_mole_fractions, reactions=reactions)
        end if
      else
        if (config%transport_enabled) then
          call read_amr_reactive_3d_checkpoint( &
            trim(amr_config%restart_file), species, config, amr_config, patch, &
            coarse_state, coarse_temperature, fine_state, fine_temperature, &
            time, step, initial_integrals, maximum_reflux, ok, message, &
            transport=transport, &
            maximum_transport_diffusivity=maximum_diffusivity, &
            minimum_transport_theta=minimum_transport_theta)
        else
          call read_amr_reactive_3d_checkpoint( &
            trim(amr_config%restart_file), species, config, amr_config, patch, &
            coarse_state, coarse_temperature, fine_state, fine_temperature, &
            time, step, initial_integrals, maximum_reflux, ok, message)
        end if
      end if
      if (.not. ok) then
        write(*, '(a)') trim(message)
        error stop 3
      end if
    else
      call initialize_reactive_problem_3d( &
        species, config, x, y, z, coarse_state, coarse_temperature, &
        base_density, mass_fractions, ok, base_mole_fractions)
      if (.not. ok) error stop "Failed to initialize coarse 3D state"
      fine_config = config
      fine_config%nx = patch%fine_nx()
      fine_config%ny = patch%fine_ny()
      fine_config%nz = patch%fine_nz()
      call initialize_reactive_problem_3d( &
        species, fine_config, xf, yf, zf, fine_state, fine_temperature, &
        fine_base_density, mass_fractions, ok, base_mole_fractions)
      if (.not. ok) error stop "Failed to initialize fine 3D state"
      call average_down_3d(coarse_state, fine_state, patch, ok)
      if (.not. ok) error stop "Failed to synchronize initial 3D AMR state"
      call recover_reactive_temperatures_3d( &
        species, coarse_state, coarse_temperature, &
        config%nx, config%ny, config%nz, recovered_temperature, ok)
      if (.not. ok) error stop "Failed to recover initial coarse temperature"
      coarse_temperature = recovered_temperature
      call composite_integrals_amr_3d( &
        coarse_state, fine_state, patch, dx, dy, dz, initial_integrals, ok)
      if (.not. ok) error stop "Failed to integrate initial 3D AMR state"
    end if
    element_diagnostics = supports_hon_element_diagnostics(species)
    initial_elements = 0.0_dp
    if (element_diagnostics) then
      if (restarted) then
        call element_integrals_from_reactive_integrals_3d( &
          species, initial_integrals, initial_elements, ok)
      else
        call composite_element_integrals_amr_3d( &
          species, patch, coarse_state, fine_state, dx, dy, dz, &
          initial_elements, ok)
      end if
      if (.not. ok) error stop "Failed to integrate initial 3D AMR elements"
    end if

    write(*, '(a)') "PeleF " // pelef_version // " " // &
      trim(application_label)
    if (selected_context) then
      write(*, '(a,1x,a)') "Bundle SHA-256:", trim(bundle_sha256)
      write(*, '(a,i0)') "Species: ", size(species)
      write(*, '(a,i0)') "Reactions: ", size(reactions)
      write(*, '(a,1x,a)') "Chemistry integrator:", &
        trim(chemistry_integrator)
    end if
    write(*, '(a,1x,a)') "Input:", trim(input_path)
    write(*, '(a,i0,a,i0,a,i0)') &
      "Coarse grid: nx=", config%nx, ", ny=", config%ny, ", nz=", config%nz
    write(*, '(a,6(i0,1x))') "Coarse patch bounds: ", &
      patch%coarse_i_lower, patch%coarse_i_upper, &
      patch%coarse_j_lower, patch%coarse_j_upper, &
      patch%coarse_k_lower, patch%coarse_k_upper
    write(*, '(a,i0)') "Refinement ratio: ", patch%refinement_ratio
    write(*, '(a,1x,a)') "Thermodynamics:", trim(config%thermo_model)
    write(*, '(a,1x,a)') "Reconstruction:", trim(config%reconstruction)
    write(*, '(a,1x,a)') "Limiter:", trim(config%limiter)
    write(*, '(a,1x,a)') "Riemann solver:", trim(config%riemann_solver)
    write(*, '(a,l2)') "Chemistry: ", config%chemistry_enabled
    write(*, '(a,l2)') "Molecular transport: ", config%transport_enabled
    if (config%transport_enabled) then
      write(*, '(a,l2)') "Viscosity: ", config%viscosity_enabled
      write(*, '(a,l2)') &
        "Thermal conduction: ", config%thermal_conduction_enabled
      write(*, '(a,l2)') &
        "Species diffusion: ", config%species_diffusion_enabled
      write(*, '(a,l2)') "Barodiffusion: ", config%barodiffusion_enabled
      write(*, '(a,es12.5)') "Transport CFL: ", config%transport_cfl
      write(*, '(a,i0)') "Transport fine subcycles per stage: ", &
        patch%refinement_ratio**2
    end if
    if (restarted) then
      write(*, '(a,1x,a)') "Restarted from checkpoint:", &
        trim(amr_config%restart_file)
      write(*, '(a,i0,a,es24.16)') &
        "Restored coarse steps: ", step, ", time: ", time
    end if

    stopped_after_checkpoint = .false.
    do while (time < config%final_time)
      if (step >= config%maximum_steps) then
        error stop "Maximum step count reached before 3D AMR final_time"
      end if
      call compute_amr_reactive_cfl_timestep_3d( &
        species, patch, coarse_state, coarse_temperature, &
        fine_state, fine_temperature, dx, dy, dz, config%cfl, hydro_dt, ok)
      if (.not. ok) error stop "Failed to compute 3D AMR CFL timestep"
      dt = hydro_dt
      if (config%transport_enabled) then
        call compute_amr_reactive_transport_timestep_3d( &
          species, transport, patch, coarse_state, coarse_temperature, &
          fine_state, fine_temperature, dx, dy, dz, config%transport_cfl, &
          config%viscosity_enabled, config%thermal_conduction_enabled, &
          config%species_diffusion_enabled, transport_dt, &
          step_maximum_diffusivity, ok)
        if (.not. ok) error stop "Failed to compute 3D AMR transport step"
        dt = min(dt, transport_dt)
        maximum_diffusivity = max( &
          maximum_diffusivity, step_maximum_diffusivity)
      end if
      dt = min(dt, config%final_time - time)
      call advance_amr_reactive_full_3d( &
        species, reactions, transport, patch, &
        coarse_state, coarse_temperature, &
        fine_state, fine_temperature, dx, dy, dz, dt, &
        config%riemann_solver, config%chemistry_enabled, &
        config%chemistry_relative_tolerance, &
        config%chemistry_absolute_tolerance, config%transport_enabled, &
        config%viscosity_enabled, config%thermal_conduction_enabled, &
        config%species_diffusion_enabled, config%barodiffusion_enabled, &
        step_transport_theta, step_reflux, ok, &
        config%reconstruction, config%limiter, &
        chemistry_integrator=chemistry_integrator)
      if (.not. ok) error stop "3D AMR reactive split rejected its candidate"
      maximum_reflux = max(maximum_reflux, step_reflux)
      minimum_transport_theta = min( &
        minimum_transport_theta, step_transport_theta)
      time = time + dt
      step = step + 1
      if (amr_config%checkpoint_interval_steps > 0) then
        if (mod(step, amr_config%checkpoint_interval_steps) == 0) then
          if (selected_context) then
            if (config%transport_enabled) then
              call write_amr_reactive_3d_checkpoint( &
                trim(amr_config%checkpoint_file), species, config, &
                amr_config, patch, coarse_state, coarse_temperature, &
                fine_state, fine_temperature, time, step, initial_integrals, &
                maximum_reflux, ok, message, bundle_sha256=bundle_sha256, &
                chemistry_integrator=chemistry_integrator, &
                base_mole_fractions=base_mole_fractions, reactions=reactions, &
                transport=transport, &
                maximum_transport_diffusivity=maximum_diffusivity, &
                minimum_transport_theta=minimum_transport_theta)
            else
              call write_amr_reactive_3d_checkpoint( &
                trim(amr_config%checkpoint_file), species, config, &
                amr_config, patch, coarse_state, coarse_temperature, &
                fine_state, fine_temperature, time, step, initial_integrals, &
                maximum_reflux, ok, message, bundle_sha256=bundle_sha256, &
                chemistry_integrator=chemistry_integrator, &
                base_mole_fractions=base_mole_fractions, reactions=reactions)
            end if
          else
            if (config%transport_enabled) then
              call write_amr_reactive_3d_checkpoint( &
                trim(amr_config%checkpoint_file), species, config, &
                amr_config, patch, coarse_state, coarse_temperature, &
                fine_state, fine_temperature, time, step, initial_integrals, &
                maximum_reflux, ok, message, transport=transport, &
                maximum_transport_diffusivity=maximum_diffusivity, &
                minimum_transport_theta=minimum_transport_theta)
            else
              call write_amr_reactive_3d_checkpoint( &
                trim(amr_config%checkpoint_file), species, config, &
                amr_config, patch, coarse_state, coarse_temperature, &
                fine_state, fine_temperature, time, step, initial_integrals, &
                maximum_reflux, ok, message)
            end if
          end if
          if (.not. ok) then
            write(*, '(a)') trim(message)
            error stop 3
          end if
          write(*, '(a,1x,a,a,i0,a,es24.16)') &
            "Wrote checkpoint:", trim(amr_config%checkpoint_file), &
            ", coarse step ", step, ", time ", time
          if (amr_config%stop_after_checkpoint) then
            stopped_after_checkpoint = .true.
            exit
          end if
        end if
      end if
    end do

    call composite_integrals_amr_3d( &
      coarse_state, fine_state, patch, dx, dy, dz, final_integrals, ok)
    if (.not. ok) error stop "Failed to integrate final 3D AMR state"
    if (config%chemistry_enabled) then
      conservation_error = maxval(abs( &
        final_integrals(1:ncons) - initial_integrals(1:ncons)) / &
        max(1.0_dp, abs(initial_integrals(1:ncons))))
    else
      conservation_error = maxval(abs(final_integrals - initial_integrals) / &
        max(1.0_dp, abs(initial_integrals)))
    end if
    species_integral_change = maxval(abs( &
      final_integrals(ncons + 1:nvar) - initial_integrals(ncons + 1:nvar)) / &
      max(1.0_dp, abs(initial_integrals(ncons + 1:nvar))))
    elemental_conservation_error = 0.0_dp
    final_elements = 0.0_dp
    if (element_diagnostics) then
      call composite_element_integrals_amr_3d( &
        species, patch, coarse_state, fine_state, dx, dy, dz, &
        final_elements, ok)
      if (.not. ok) error stop "Failed to integrate final 3D AMR elements"
      elemental_conservation_error = maxval(abs( &
        final_elements - initial_elements) / &
        max(1.0e-30_dp, abs(initial_elements)))
    end if
    covered_nx = patch%coarse_i_upper - patch%coarse_i_lower + 1
    covered_ny = patch%coarse_j_upper - patch%coarse_j_lower + 1
    covered_nz = patch%coarse_k_upper - patch%coarse_k_lower + 1
    allocate(restricted(nvar, covered_nx, covered_ny, covered_nz))
    call restrict_average_3d(fine_state, patch, restricted, ok)
    if (.not. ok) error stop "Failed to restrict final 3D AMR state"
    synchronization_error = maxval(abs(restricted - coarse_state(:, &
      patch%coarse_i_lower:patch%coarse_i_upper, &
      patch%coarse_j_lower:patch%coarse_j_upper, &
      patch%coarse_k_lower:patch%coarse_k_upper))) / &
      max(1.0_dp, maxval(abs(restricted)))
    call level_extrema( &
      species, coarse_state, coarse_temperature, coarse_closure, ok)
    if (.not. ok) error stop "Final coarse 3D AMR state is not physical"
    call level_extrema(species, fine_state, fine_temperature, fine_closure, ok)
    if (.not. ok) error stop "Final fine 3D AMR state is not physical"
    if (conservation_error > 5.0e-11_dp) then
      error stop "3D AMR composite conservation gate failed"
    end if
    if (config%chemistry_enabled .and. element_diagnostics .and. &
        elemental_conservation_error > 5.0e-10_dp) then
      error stop "3D AMR elemental conservation gate failed"
    end if
    if (synchronization_error > 5.0e-13_dp) then
      error stop "3D AMR average-down gate failed"
    end if
    if (max(coarse_closure, fine_closure) > 5.0e-11_dp) then
      error stop "3D AMR species-closure gate failed"
    end if
    call write_reactive_3d_csv( &
      trim(amr_config%coarse_output_file), species, x, y, z, &
      coarse_state, coarse_temperature, config%nx, config%ny, config%nz, &
      time, ok, message)
    if (.not. ok) then
      write(*, '(a)') trim(message)
      error stop 3
    end if
    call write_reactive_3d_csv( &
      trim(amr_config%fine_output_file), species, xf, yf, zf, &
      fine_state, fine_temperature, patch%fine_nx(), patch%fine_ny(), &
      patch%fine_nz(), time, ok, message)
    if (.not. ok) then
      write(*, '(a)') trim(message)
      error stop 3
    end if
    write(*, '(a,i0)') "Completed coarse steps: ", step
    write(*, '(a,i0)') "Completed hydro fine substeps: ", &
      step * patch%refinement_ratio
    write(*, '(a,es24.16)') "Final time: ", time
    write(*, '(a,es24.16)') "Maximum reflux correction: ", maximum_reflux
    if (config%transport_enabled) then
      write(*, '(a,i0)') "Completed transport fine substeps: ", &
        4 * step * patch%refinement_ratio**2
      write(*, '(a,es24.16)') "Maximum transport diffusivity: ", &
        maximum_diffusivity
      write(*, '(a,es24.16)') "Minimum transport theta: ", &
        minimum_transport_theta
    end if
    write(*, '(a,es24.16)') "Composite conservation error: ", &
      conservation_error
    if (element_diagnostics) then
      write(*, '(a,es24.16)') "Maximum elemental conservation error: ", &
        elemental_conservation_error
    else
      write(*, '(a)') "Maximum elemental conservation error: not available"
    end if
    write(*, '(a,es24.16)') "Maximum species integral change: ", &
      species_integral_change
    write(*, '(a,es24.16)') "Average-down synchronization error: ", &
      synchronization_error
    write(*, '(a,es24.16)') "Maximum species closure error: ", &
      max(coarse_closure, fine_closure)
    write(*, '(a,1x,a)') "Wrote coarse CSV:", &
      trim(amr_config%coarse_output_file)
    write(*, '(a,1x,a)') "Wrote fine CSV:", &
      trim(amr_config%fine_output_file)
    if (stopped_after_checkpoint) then
      write(*, '(a,1x,a)') "Stopped after checkpoint:", &
        trim(amr_config%checkpoint_file)
    end if
  end subroutine run_amr_reactive_3d_application

  pure logical function supports_hon_element_diagnostics(species) &
      result(supported)
    type(nasa7_species), intent(in) :: species(:)

    integer :: index

    supported = size(species) > 0
    do index = 1, size(species)
      select case (trim(species(index)%name))
      case ("H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", &
          "N2", "AR")
      case default
        supported = .false.
        return
      end select
    end do
  end function supports_hon_element_diagnostics

  subroutine fine_patch_centers(patch, config, dx, dy, dz, x, y, z)
    type(amr_patch_3d), intent(in) :: patch
    type(reactive_3d_config), intent(in) :: config
    real(dp), intent(in) :: dx, dy, dz
    real(dp), intent(out) :: x(:), y(:), z(:)
    integer :: index

    do index = 1, size(x)
      x(index) = config%x_lower + &
        real(patch%coarse_i_lower - 1, dp) * dx + &
        (real(index, dp) - 0.5_dp) * dx / &
        real(patch%refinement_ratio, dp)
    end do
    do index = 1, size(y)
      y(index) = config%y_lower + &
        real(patch%coarse_j_lower - 1, dp) * dy + &
        (real(index, dp) - 0.5_dp) * dy / &
        real(patch%refinement_ratio, dp)
    end do
    do index = 1, size(z)
      z(index) = config%z_lower + &
        real(patch%coarse_k_lower - 1, dp) * dz + &
        (real(index, dp) - 0.5_dp) * dz / &
        real(patch%refinement_ratio, dp)
    end do
  end subroutine fine_patch_centers

  subroutine level_extrema(species, state, temperature, closure, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(out) :: closure
    logical, intent(out) :: ok

    real(dp) :: minimum_density, maximum_density
    real(dp) :: minimum_pressure, maximum_pressure
    real(dp) :: minimum_temperature, maximum_temperature, maximum_speed

    call reactive_extrema_3d( &
      species, state, temperature, size(state, 2), size(state, 3), &
      size(state, 4), minimum_density, maximum_density, minimum_pressure, &
      maximum_pressure, minimum_temperature, maximum_temperature, &
      maximum_speed, closure, ok)
    if (.not. ok) return
    ok = all(ieee_is_finite([ &
      minimum_density, maximum_density, minimum_pressure, maximum_pressure, &
      minimum_temperature, maximum_temperature, maximum_speed, closure]))
    if (.not. ok) return
    ok = minimum_density > 0.0_dp .and. &
      minimum_pressure > 0.0_dp .and. minimum_temperature > 0.0_dp .and. &
      maximum_density > 0.0_dp .and. maximum_pressure > 0.0_dp .and. &
      maximum_temperature > 0.0_dp .and. maximum_speed >= 0.0_dp
  end subroutine level_extrema

end module amr_reactive_3d_application_mod
