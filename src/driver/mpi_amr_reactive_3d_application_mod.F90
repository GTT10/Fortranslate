module mpi_amr_reactive_3d_application_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: error_unit, int64
  use mpi_f08
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
  use nasa7_thermo_mod, only: nasa7_species, valid_nasa7_species
  use elementary_kinetics_mod, only: &
    elementary_reaction, valid_elementary_reaction
  use gas_transport_mod, only: &
    gas_transport_species, compatible_transport_database
  use mesh_3d_mod, only: uniform_cell_centers_3d
  use state_indices_mod, only: ncons
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
  use mpi_amr_sparse_reactive_3d_mod, only: &
    mpi_amr_sparse_distribution_3d, mpi_amr_sparse_hierarchy_3d, &
    initialize_mpi_amr_sparse_distribution_3d, &
    scatter_mpi_amr_sparse_hierarchy_3d, &
    gather_mpi_amr_sparse_hierarchy_3d, &
    compute_mpi_amr_sparse_cfl_timestep_3d, &
    compute_mpi_amr_sparse_transport_timestep_3d, &
    advance_mpi_amr_sparse_full_3d
  use amr_reactive_3d_mod, only: &
    composite_element_integrals_amr_3d, &
    element_integrals_from_reactive_integrals_3d
  use amr_reactive_3d_checkpoint_mod, only: &
    write_amr_reactive_3d_checkpoint, read_amr_reactive_3d_checkpoint, &
    amr_reactive_3d_transport_operator
  use reactive_csv_io_3d_mod, only: write_reactive_3d_csv
  implicit none
  private

  public :: run_mpi_amr_reactive_3d_application
  public :: validate_mpi_amr_selected_context_3d
  public :: validate_mpi_amr_transport_context_3d

contains

  subroutine run_mpi_amr_reactive_3d_application( &
      comm, input_path, application_label, config, amr_config, species, &
      reactions, transport, base_mole_fractions, coarse_output, fine_output, &
      bundle_sha256, chemistry_integrator)
    type(MPI_Comm), intent(in) :: comm
    character(len=*), intent(in) :: input_path, application_label
    type(reactive_3d_config), intent(in) :: config
    type(amr_reactive_3d_config), intent(in) :: amr_config
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: base_mole_fractions(:)
    character(len=*), intent(in) :: coarse_output, fine_output
    character(len=*), intent(in), optional :: bundle_sha256
    character(len=*), intent(in), optional :: chemistry_integrator

    type(reactive_3d_config) :: fine_config
    type(amr_patch_3d) :: patch
    type(mpi_amr_sparse_distribution_3d) :: distribution
    type(mpi_amr_sparse_hierarchy_3d) :: hierarchy
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
    logical :: ok, restarted, stopped_after_checkpoint, root_ok
    logical :: selected_context, element_diagnostics
    logical :: transport_enabled, transport_policy_ok
    logical :: any_selected_context, all_selected_context
    integer :: step, nvar, covered_nx, covered_ny, covered_nz
    integer :: ierr, rank, nranks

    call MPI_Comm_rank(comm, rank, ierr)
    if (ierr /= MPI_SUCCESS) error stop "MPI_Comm_rank failed"
    call MPI_Comm_size(comm, nranks, ierr)
    if (ierr /= MPI_SUCCESS) call abort_run("MPI_Comm_size failed", 2)

    selected_context = present(bundle_sha256) .or. &
      present(chemistry_integrator)
    call MPI_Allreduce( &
      selected_context, any_selected_context, 1, MPI_LOGICAL, MPI_LOR, &
      comm, ierr)
    if (ierr /= MPI_SUCCESS) call abort_run( &
      "Failed to reduce selected sparse MPI 3D AMR context presence", 2)
    call MPI_Allreduce( &
      selected_context, all_selected_context, 1, MPI_LOGICAL, MPI_LAND, &
      comm, ierr)
    if (ierr /= MPI_SUCCESS) call abort_run( &
      "Failed to reduce selected sparse MPI 3D AMR context presence", 2)
    if (any_selected_context .neqv. all_selected_context) then
      call collective_abort_run( &
        "Selected sparse MPI 3D AMR context presence differs by rank", 2)
    end if
    selected_context = all_selected_context
    if (selected_context) then
      ok = present(bundle_sha256) .and. present(chemistry_integrator)
      call collective_require( &
        ok, "Selected sparse MPI 3D AMR context is incomplete", 2)
    end if
    ok = size(species) > 0 .and. size(reactions) > 0
    if (ok) ok = size(base_mole_fractions) == size(species)
    if (ok) ok = size(transport) == size(species)
    if (ok) ok = all(ieee_is_finite(base_mole_fractions))
    if (ok) ok = minval(base_mole_fractions) >= 0.0_dp
    if (ok) ok = abs(sum(base_mole_fractions) - 1.0_dp) <= 5.0e-10_dp
    if (ok) ok = compatible_transport_database(species, transport)
    call collective_require( &
      ok, "Sparse MPI 3D AMR application data are incompatible", 2)
    transport_enabled = config%transport_enabled
    call consensus_logical_vector( &
      comm, [transport_enabled], transport_policy_ok)
    if (.not. transport_policy_ok) then
      call collective_abort_run( &
        "Sparse MPI 3D AMR transport-enabled policy differs by rank", 3)
    end if
    if (transport_enabled) then
      call validate_mpi_amr_transport_context_3d( &
        comm, config, species, transport, ok)
      call collective_require( &
        ok, "Sparse MPI 3D AMR transport context mismatch", 3)
    end if
    ok = len_trim(coarse_output) > 0 .and. len_trim(fine_output) > 0
    if (ok) ok = trim(coarse_output) /= trim(fine_output)
    if (ok .and. len_trim(amr_config%checkpoint_file) > 0) then
      ok = trim(coarse_output) /= trim(amr_config%checkpoint_file) .and. &
        trim(fine_output) /= trim(amr_config%checkpoint_file)
    end if
    if (ok .and. len_trim(amr_config%restart_file) > 0) then
      ok = trim(coarse_output) /= trim(amr_config%restart_file) .and. &
        trim(fine_output) /= trim(amr_config%restart_file)
    end if
    call collective_require(ok, "MPI AMR I/O paths collide", 2)
    if (selected_context) then
      call validate_mpi_amr_selected_context_3d( &
        comm, bundle_sha256, chemistry_integrator, &
        config%chemistry_relative_tolerance, &
        config%chemistry_absolute_tolerance, base_mole_fractions, &
        species, reactions, transport, ok)
      call collective_require( &
        ok, "Selected sparse MPI 3D AMR rank context mismatch", 3)
    end if

    call initialize_amr_patch_3d( &
      config%nx, config%ny, config%nz, &
      amr_config%coarse_i_lower, amr_config%coarse_i_upper, &
      amr_config%coarse_j_lower, amr_config%coarse_j_upper, &
      amr_config%coarse_k_lower, amr_config%coarse_k_upper, &
      amr_config%refinement_ratio, patch, ok)
    if (ok) ok = patch%is_strictly_interior()
    call collective_require(ok, "Failed to construct MPI 3D AMR patch", 2)
    call initialize_mpi_amr_sparse_distribution_3d( &
      patch, comm, distribution, ok)
    call collective_require( &
      ok, "Failed to distribute sparse MPI 3D AMR slabs", 2)

    nvar = reactive_nvar(size(species))
    allocate(x(config%nx), y(config%ny), z(config%nz))
    allocate(xf(patch%fine_nx()), yf(patch%fine_ny()), zf(patch%fine_nz()))
    if (rank == 0) then
      allocate(coarse_state(nvar, config%nx, config%ny, config%nz))
      allocate(coarse_temperature(config%nx, config%ny, config%nz))
      allocate(fine_state(nvar, patch%fine_nx(), patch%fine_ny(), &
        patch%fine_nz()))
      allocate(fine_temperature( &
        patch%fine_nx(), patch%fine_ny(), patch%fine_nz()))
      allocate(recovered_temperature(config%nx, config%ny, config%nz))
    else
      allocate(coarse_state(0, 0, 0, 0), coarse_temperature(0, 0, 0))
      allocate(fine_state(0, 0, 0, 0), fine_temperature(0, 0, 0))
      allocate(recovered_temperature(0, 0, 0))
    end if
    allocate(mass_fractions(size(species)))
    allocate(initial_integrals(nvar), final_integrals(nvar))
    if (rank == 0) then
      coarse_state = 0.0_dp
      coarse_temperature = 0.0_dp
      fine_state = 0.0_dp
      fine_temperature = 0.0_dp
    end if
    initial_integrals = 0.0_dp
    call uniform_cell_centers_3d( &
      config%nx, config%ny, config%nz, config%x_lower, config%x_upper, &
      config%y_lower, config%y_upper, config%z_lower, config%z_upper, &
      x, y, z, dx, dy, dz)
    call fine_patch_centers(patch, config, dx, dy, dz, xf, yf, zf)
    time = 0.0_dp
    step = 0
    maximum_reflux = 0.0_dp
    maximum_diffusivity = 0.0_dp
    minimum_transport_theta = 1.0_dp
    restarted = len_trim(amr_config%restart_file) > 0
    element_diagnostics = supports_hon_element_diagnostics(species)
    initial_elements = 0.0_dp
    root_ok = .true.
    if (rank == 0) then
      if (restarted) then
        if (selected_context) then
          if (transport_enabled) then
            call read_amr_reactive_3d_checkpoint( &
              trim(amr_config%restart_file), species, config, amr_config, &
              patch, coarse_state, coarse_temperature, fine_state, &
              fine_temperature, time, step, initial_integrals, &
              maximum_reflux, root_ok, message, bundle_sha256=bundle_sha256, &
              chemistry_integrator=chemistry_integrator, &
              base_mole_fractions=base_mole_fractions, reactions=reactions, &
              transport=transport, &
              maximum_transport_diffusivity=maximum_diffusivity, &
              minimum_transport_theta=minimum_transport_theta)
          else
            call read_amr_reactive_3d_checkpoint( &
              trim(amr_config%restart_file), species, config, amr_config, &
              patch, coarse_state, coarse_temperature, fine_state, &
              fine_temperature, time, step, initial_integrals, &
              maximum_reflux, root_ok, message, bundle_sha256=bundle_sha256, &
              chemistry_integrator=chemistry_integrator, &
              base_mole_fractions=base_mole_fractions, reactions=reactions)
          end if
        else
          if (transport_enabled) then
            call read_amr_reactive_3d_checkpoint( &
              trim(amr_config%restart_file), species, config, amr_config, &
              patch, coarse_state, coarse_temperature, fine_state, &
              fine_temperature, time, step, initial_integrals, &
              maximum_reflux, root_ok, message, transport=transport, &
              maximum_transport_diffusivity=maximum_diffusivity, &
              minimum_transport_theta=minimum_transport_theta)
          else
            call read_amr_reactive_3d_checkpoint( &
              trim(amr_config%restart_file), species, config, amr_config, &
              patch, coarse_state, coarse_temperature, fine_state, &
              fine_temperature, time, step, initial_integrals, &
              maximum_reflux, root_ok, message)
          end if
        end if
      else
        call initialize_reactive_problem_3d( &
          species, config, x, y, z, coarse_state, coarse_temperature, &
          base_density, mass_fractions, root_ok, base_mole_fractions)
        if (root_ok) then
          fine_config = config
          fine_config%nx = patch%fine_nx()
          fine_config%ny = patch%fine_ny()
          fine_config%nz = patch%fine_nz()
          call initialize_reactive_problem_3d( &
            species, fine_config, xf, yf, zf, fine_state, fine_temperature, &
            fine_base_density, mass_fractions, root_ok, base_mole_fractions)
        end if
        if (root_ok) call average_down_3d( &
          coarse_state, fine_state, patch, root_ok)
        if (root_ok) then
          call recover_reactive_temperatures_3d( &
            species, coarse_state, coarse_temperature, config%nx, &
            config%ny, config%nz, recovered_temperature, root_ok)
        end if
        if (root_ok) coarse_temperature = recovered_temperature
        if (root_ok) then
          call composite_integrals_amr_3d( &
            coarse_state, fine_state, patch, dx, dy, dz, &
            initial_integrals, root_ok)
        end if
        if (.not. root_ok) message = "Failed to initialize MPI 3D AMR state"
      end if
      if (root_ok .and. element_diagnostics) then
        if (restarted) then
          call element_integrals_from_reactive_integrals_3d( &
            species, initial_integrals, initial_elements, root_ok)
        else
          call composite_element_integrals_amr_3d( &
            species, patch, coarse_state, fine_state, dx, dy, dz, &
            initial_elements, root_ok)
        end if
      end if
    end if
    call root_require(root_ok, message, 3)
    call MPI_Bcast(time, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call collective_require( &
      ierr == MPI_SUCCESS, "Failed to broadcast sparse MPI 3D AMR time", 3)
    call MPI_Bcast(step, 1, MPI_INTEGER, 0, comm, ierr)
    call collective_require( &
      ierr == MPI_SUCCESS, "Failed to broadcast sparse MPI 3D AMR step", 3)
    call MPI_Bcast( &
      maximum_reflux, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call collective_require( &
      ierr == MPI_SUCCESS, "Failed to broadcast sparse MPI 3D AMR reflux", 3)
    call MPI_Bcast( &
      maximum_diffusivity, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call collective_require( &
      ierr == MPI_SUCCESS, &
      "Failed to broadcast sparse MPI 3D AMR transport diffusivity", 3)
    call MPI_Bcast( &
      minimum_transport_theta, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call collective_require( &
      ierr == MPI_SUCCESS, &
      "Failed to broadcast sparse MPI 3D AMR transport theta", 3)
    call MPI_Bcast( &
      initial_integrals, size(initial_integrals), MPI_DOUBLE_PRECISION, &
      0, comm, ierr)
    call collective_require( &
      ierr == MPI_SUCCESS, &
      "Failed to broadcast sparse MPI 3D AMR integrals", 3)
    call scatter_mpi_amr_sparse_hierarchy_3d( &
      distribution, species, patch, 0, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, hierarchy, ok)
    call collective_require( &
      ok, "Failed to scatter sparse MPI 3D AMR hierarchy", 3)
    deallocate(coarse_state, coarse_temperature, fine_state, fine_temperature)

    if (rank == 0) then
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
      write(*, '(a,i0)') "MPI ranks: ", nranks
      write(*, '(a,i0,a,i0)') "Coarse x planes per rank: min=", &
        minval(distribution%coarse_counts), ", max=", &
        maxval(distribution%coarse_counts)
      write(*, '(a,i0,a,i0)') "Fine x planes per rank: min=", &
        minval(distribution%fine_counts), ", max=", &
        maxval(distribution%fine_counts)
      write(*, '(a,1x,a)') "Reconstruction:", trim(config%reconstruction)
      write(*, '(a,1x,a)') "Limiter:", trim(config%limiter)
      write(*, '(a,l2)') "Chemistry: ", config%chemistry_enabled
      write(*, '(a,l2)') "Molecular transport: ", transport_enabled
      if (transport_enabled) then
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
    end if

    stopped_after_checkpoint = .false.
    do while (time < config%final_time)
      ok = step < config%maximum_steps
      call collective_require( &
        ok, "Maximum step count reached before MPI 3D AMR final_time", 4)
      call compute_mpi_amr_sparse_cfl_timestep_3d( &
        distribution, species, patch, hierarchy, &
        dx, dy, dz, config%cfl, hydro_dt, ok)
      call collective_require(ok, "Failed to compute MPI 3D AMR timestep", 4)
      dt = hydro_dt
      if (transport_enabled) then
        call compute_mpi_amr_sparse_transport_timestep_3d( &
          distribution, species, transport, patch, hierarchy, &
          dx, dy, dz, config%transport_cfl, config%viscosity_enabled, &
          config%thermal_conduction_enabled, &
          config%species_diffusion_enabled, transport_dt, &
          step_maximum_diffusivity, ok)
        call collective_require( &
          ok, "Failed to compute MPI 3D AMR transport timestep", 4)
        dt = min(dt, transport_dt)
        maximum_diffusivity = max( &
          maximum_diffusivity, step_maximum_diffusivity)
      end if
      dt = min(dt, config%final_time - time)
      if (present(chemistry_integrator)) then
        call advance_mpi_amr_sparse_full_3d( &
          distribution, species, reactions, transport, patch, hierarchy, &
          dx, dy, dz, dt, config%riemann_solver, config%reconstruction, &
          config%limiter, config%chemistry_enabled, &
          config%chemistry_relative_tolerance, &
          config%chemistry_absolute_tolerance, transport_enabled, &
          config%viscosity_enabled, config%thermal_conduction_enabled, &
          config%species_diffusion_enabled, config%barodiffusion_enabled, &
          step_transport_theta, step_reflux, ok, &
          chemistry_integrator=chemistry_integrator)
      else
        call advance_mpi_amr_sparse_full_3d( &
          distribution, species, reactions, transport, patch, hierarchy, &
          dx, dy, dz, dt, config%riemann_solver, config%reconstruction, &
          config%limiter, config%chemistry_enabled, &
          config%chemistry_relative_tolerance, &
          config%chemistry_absolute_tolerance, transport_enabled, &
          config%viscosity_enabled, config%thermal_conduction_enabled, &
          config%species_diffusion_enabled, config%barodiffusion_enabled, &
          step_transport_theta, step_reflux, ok)
      end if
      call collective_require( &
        ok, "MPI 3D AMR reactive split rejected its candidate", 4)
      maximum_reflux = max(maximum_reflux, step_reflux)
      minimum_transport_theta = min( &
        minimum_transport_theta, step_transport_theta)
      time = time + dt
      step = step + 1
      if (amr_config%checkpoint_interval_steps > 0) then
        if (mod(step, amr_config%checkpoint_interval_steps) == 0) then
          call gather_mpi_amr_sparse_hierarchy_3d( &
            distribution, patch, 0, hierarchy, coarse_state, &
            coarse_temperature, fine_state, fine_temperature, ok)
          call collective_require( &
            ok, "Failed to gather sparse MPI 3D AMR checkpoint state", 5)
          root_ok = .true.
          if (rank == 0) then
            if (selected_context) then
              if (transport_enabled) then
                call write_amr_reactive_3d_checkpoint( &
                  trim(amr_config%checkpoint_file), species, config, &
                  amr_config, patch, coarse_state, coarse_temperature, &
                  fine_state, fine_temperature, time, step, &
                  initial_integrals, maximum_reflux, root_ok, message, &
                  bundle_sha256=bundle_sha256, &
                  chemistry_integrator=chemistry_integrator, &
                  base_mole_fractions=base_mole_fractions, &
                  reactions=reactions, transport=transport, &
                  maximum_transport_diffusivity=maximum_diffusivity, &
                  minimum_transport_theta=minimum_transport_theta)
              else
                call write_amr_reactive_3d_checkpoint( &
                  trim(amr_config%checkpoint_file), species, config, &
                  amr_config, patch, coarse_state, coarse_temperature, &
                  fine_state, fine_temperature, time, step, &
                  initial_integrals, maximum_reflux, root_ok, message, &
                  bundle_sha256=bundle_sha256, &
                  chemistry_integrator=chemistry_integrator, &
                  base_mole_fractions=base_mole_fractions, &
                  reactions=reactions)
              end if
            else
              if (transport_enabled) then
                call write_amr_reactive_3d_checkpoint( &
                  trim(amr_config%checkpoint_file), species, config, &
                  amr_config, patch, coarse_state, coarse_temperature, &
                  fine_state, fine_temperature, time, step, &
                  initial_integrals, maximum_reflux, root_ok, message, &
                  transport=transport, &
                  maximum_transport_diffusivity=maximum_diffusivity, &
                  minimum_transport_theta=minimum_transport_theta)
              else
                call write_amr_reactive_3d_checkpoint( &
                  trim(amr_config%checkpoint_file), species, config, &
                  amr_config, patch, coarse_state, coarse_temperature, &
                  fine_state, fine_temperature, time, step, &
                  initial_integrals, maximum_reflux, root_ok, message)
              end if
            end if
          end if
          call root_require(root_ok, message, 5)
          deallocate(coarse_state, coarse_temperature, fine_state, &
            fine_temperature)
          if (rank == 0) then
            write(*, '(a,1x,a,a,i0,a,es24.16)') &
              "Wrote checkpoint:", trim(amr_config%checkpoint_file), &
              ", coarse step ", step, ", time ", time
          end if
          if (amr_config%stop_after_checkpoint) then
            stopped_after_checkpoint = .true.
            exit
          end if
        end if
      end if
    end do

    call gather_mpi_amr_sparse_hierarchy_3d( &
      distribution, patch, 0, hierarchy, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, ok)
    call collective_require( &
      ok, "Failed to gather final sparse MPI 3D AMR state", 6)
    covered_nx = patch%coarse_i_upper - patch%coarse_i_lower + 1
    covered_ny = patch%coarse_j_upper - patch%coarse_j_lower + 1
    covered_nz = patch%coarse_k_upper - patch%coarse_k_lower + 1
    root_ok = .true.
    if (rank == 0) then
      call composite_integrals_amr_3d( &
        coarse_state, fine_state, patch, dx, dy, dz, final_integrals, root_ok)
      if (root_ok) then
        if (config%chemistry_enabled) then
          conservation_error = maxval(abs( &
            final_integrals(1:ncons) - initial_integrals(1:ncons)) / &
            max(1.0_dp, abs(initial_integrals(1:ncons))))
        else
          conservation_error = maxval( &
            abs(final_integrals - initial_integrals) / &
            max(1.0_dp, abs(initial_integrals)))
        end if
        species_integral_change = maxval(abs( &
          final_integrals(ncons + 1:nvar) - &
          initial_integrals(ncons + 1:nvar)) / &
          max(1.0_dp, abs(initial_integrals(ncons + 1:nvar))))
      end if
      elemental_conservation_error = 0.0_dp
      final_elements = 0.0_dp
      if (root_ok .and. element_diagnostics) then
        call composite_element_integrals_amr_3d( &
          species, patch, coarse_state, fine_state, dx, dy, dz, &
          final_elements, root_ok)
      end if
      if (root_ok .and. element_diagnostics) then
        elemental_conservation_error = maxval(abs( &
          final_elements - initial_elements) / &
          max(1.0e-30_dp, abs(initial_elements)))
      end if
      if (root_ok) then
        allocate(restricted(nvar, covered_nx, covered_ny, covered_nz))
        call restrict_average_3d(fine_state, patch, restricted, root_ok)
      end if
      if (root_ok) then
        synchronization_error = maxval(abs(restricted - coarse_state(:, &
          patch%coarse_i_lower:patch%coarse_i_upper, &
          patch%coarse_j_lower:patch%coarse_j_upper, &
          patch%coarse_k_lower:patch%coarse_k_upper))) / &
          max(1.0_dp, maxval(abs(restricted)))
        call level_extrema( &
          species, coarse_state, coarse_temperature, coarse_closure, root_ok)
      end if
      if (root_ok) call level_extrema( &
        species, fine_state, fine_temperature, fine_closure, root_ok)
      if (root_ok) then
        root_ok = conservation_error <= 5.0e-11_dp .and. &
          synchronization_error <= 5.0e-13_dp .and. &
          max(coarse_closure, fine_closure) <= 5.0e-11_dp
      end if
      if (root_ok .and. config%chemistry_enabled .and. &
          element_diagnostics) then
        root_ok = elemental_conservation_error <= 5.0e-10_dp
      end if
      if (.not. root_ok) &
        message = "Final sparse MPI 3D AMR invariant gate failed"
    end if
    call root_require(root_ok, message, 6)

    root_ok = .true.
    if (rank == 0) then
      call write_reactive_3d_csv( &
        trim(coarse_output), species, x, y, z, coarse_state, &
        coarse_temperature, config%nx, config%ny, config%nz, time, &
        root_ok, message)
      if (root_ok) then
        call write_reactive_3d_csv( &
          trim(fine_output), species, xf, yf, zf, fine_state, &
          fine_temperature, patch%fine_nx(), patch%fine_ny(), &
          patch%fine_nz(), time, root_ok, message)
      end if
    end if
    call root_require(root_ok, message, 7)
    if (rank == 0) then
      write(*, '(a,i0)') "Completed coarse steps: ", step
      write(*, '(a,i0)') "Completed hydro fine substeps: ", &
        step * patch%refinement_ratio
      write(*, '(a,es24.16)') "Final time: ", time
      write(*, '(a,es24.16)') "Maximum reflux correction: ", maximum_reflux
      if (transport_enabled) then
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
        write(*, '(a)') &
          "Maximum elemental conservation error: not available"
      end if
      write(*, '(a,es24.16)') "Maximum species integral change: ", &
        species_integral_change
      write(*, '(a,es24.16)') "Average-down synchronization error: ", &
        synchronization_error
      write(*, '(a,es24.16)') "Maximum species closure error: ", &
        max(coarse_closure, fine_closure)
      write(*, '(a,1x,a)') "Wrote coarse CSV:", trim(coarse_output)
      write(*, '(a,1x,a)') "Wrote fine CSV:", trim(fine_output)
      if (stopped_after_checkpoint) then
        write(*, '(a,1x,a)') "Stopped after checkpoint:", &
          trim(amr_config%checkpoint_file)
      end if
    end if

  contains

    subroutine collective_require(local_ok, failure_message, code)
      logical, intent(in) :: local_ok
      character(len=*), intent(in) :: failure_message
      integer, intent(in) :: code

      logical :: global_ok
      integer :: local_ierr

      call MPI_Allreduce( &
        local_ok, global_ok, 1, MPI_LOGICAL, MPI_LAND, comm, local_ierr)
      if (local_ierr /= MPI_SUCCESS) then
        call abort_run(failure_message, code)
      end if
      if (.not. global_ok) call collective_abort_run(failure_message, code)
    end subroutine collective_require

    subroutine root_require(local_root_ok, failure_message, code)
      logical, intent(in) :: local_root_ok
      character(len=*), intent(in) :: failure_message
      integer, intent(in) :: code

      logical :: broadcast_ok
      integer :: local_ierr

      broadcast_ok = .true.
      if (rank == 0) broadcast_ok = local_root_ok
      call MPI_Bcast( &
        broadcast_ok, 1, MPI_LOGICAL, 0, comm, local_ierr)
      if (local_ierr /= MPI_SUCCESS) then
        call abort_run(failure_message, code)
      end if
      if (.not. broadcast_ok) &
        call collective_abort_run(failure_message, code)
    end subroutine root_require

    subroutine collective_abort_run(reason, code)
      character(len=*), intent(in) :: reason
      integer, intent(in) :: code

      integer :: abort_ierr, sync_ierr

      if (rank == 0) then
        write(error_unit, '(a)') trim(reason)
        flush(error_unit)
      end if
      call MPI_Barrier(comm, sync_ierr)
      call MPI_Abort(comm, code, abort_ierr)
      error stop code
    end subroutine collective_abort_run

    subroutine abort_run(reason, code)
      character(len=*), intent(in) :: reason
      integer, intent(in) :: code

      integer :: abort_ierr

      if (rank == 0) then
        write(error_unit, '(a)') trim(reason)
        flush(error_unit)
      end if
      call MPI_Abort(comm, code, abort_ierr)
      error stop code
    end subroutine abort_run

  end subroutine run_mpi_amr_reactive_3d_application

  subroutine validate_mpi_amr_transport_context_3d( &
      comm, config, species, transport, ok)
    type(MPI_Comm), intent(in) :: comm
    type(reactive_3d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    logical, intent(out) :: ok

    character(len=32) :: policy_names(3)
    character(len=24), allocatable :: species_names(:), transport_names(:)
    integer, allocatable :: transport_geometries(:)
    integer :: sizes(2), species_index, value_offset, ierr
    logical :: controls(4), local_ok, global_ok, field_ok
    real(dp) :: transport_cfl(1)
    real(dp), allocatable :: transport_values(:)

    sizes = [size(species), size(transport)]
    call consensus_integer_vector(comm, sizes, field_ok)
    if (.not. field_ok) then
      ok = .false.
      return
    end if
    local_ok = config%transport_enabled .and. size(species) > 0
    if (local_ok) local_ok = compatible_transport_database(species, transport)
    if (local_ok) local_ok = ieee_is_finite(config%transport_cfl)
    if (local_ok) local_ok = config%transport_cfl > 0.0_dp
    if (local_ok) local_ok = .not. config%barodiffusion_enabled .or. &
      config%species_diffusion_enabled
    if (local_ok) then
      local_ok = trim(config%reconstruction) == "pcm" .or. &
        trim(config%reconstruction) == "characteristic_plm"
    end if
    if (local_ok) then
      local_ok = trim(config%limiter) == "minmod" .or. &
        trim(config%limiter) == "mc"
    end if
    call MPI_Allreduce( &
      local_ok, global_ok, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    if (.not. global_ok) then
      ok = .false.
      return
    end if

    allocate(species_names(size(species)), transport_names(size(transport)))
    allocate(transport_geometries(size(transport)))
    allocate(transport_values(5 * size(transport)))
    do species_index = 1, size(species)
      species_names(species_index) = species(species_index)%name
      transport_names(species_index) = transport(species_index)%name
      transport_geometries(species_index) = transport(species_index)%geometry
      value_offset = 5 * (species_index - 1)
      transport_values(value_offset + 1:value_offset + 5) = [ &
        transport(species_index)%well_depth, &
        transport(species_index)%diameter, &
        transport(species_index)%dipole, &
        transport(species_index)%polarizability, &
        transport(species_index)%rotational_relaxation]
    end do
    controls = [ &
      config%viscosity_enabled, config%thermal_conduction_enabled, &
      config%species_diffusion_enabled, config%barodiffusion_enabled]
    transport_cfl(1) = config%transport_cfl
    policy_names = [character(len=32) :: &
      trim(config%reconstruction), trim(config%limiter), &
      amr_reactive_3d_transport_operator]

    ok = .true.
    call consensus_character_vector(comm, species_names, field_ok)
    ok = ok .and. field_ok
    call consensus_character_vector(comm, transport_names, field_ok)
    ok = ok .and. field_ok
    call consensus_integer_vector(comm, transport_geometries, field_ok)
    ok = ok .and. field_ok
    call consensus_real_vector(comm, transport_values, field_ok)
    ok = ok .and. field_ok
    call consensus_logical_vector(comm, controls, field_ok)
    ok = ok .and. field_ok
    call consensus_real_vector(comm, transport_cfl, field_ok)
    ok = ok .and. field_ok
    call consensus_character_vector(comm, policy_names, field_ok)
    ok = ok .and. field_ok
  end subroutine validate_mpi_amr_transport_context_3d

  subroutine validate_mpi_amr_selected_context_3d( &
      comm, bundle_sha256, chemistry_integrator, relative_tolerance, &
      absolute_tolerance, composition, species, reactions, transport, ok)
    type(MPI_Comm), intent(in) :: comm
    character(len=*), intent(in) :: bundle_sha256, chemistry_integrator
    real(dp), intent(in) :: relative_tolerance, absolute_tolerance
    real(dp), intent(in) :: composition(:)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    logical, intent(out) :: ok

    character(len=64) :: bundle_values(1)
    character(len=32) :: integrator_values(1)
    character(len=24), allocatable :: species_names(:), transport_names(:)
    character(len=128), allocatable :: reaction_equations(:)
    integer, allocatable :: transport_geometries(:), reaction_kinds(:)
    integer :: sizes(4), nspecies, nreactions
    real(dp), allocatable :: species_values(:), transport_values(:)
    real(dp), allocatable :: reaction_values(:)
    real(dp) :: chemistry_controls(2)
    logical, allocatable :: reversible(:), efficiencies_allocated(:)
    logical, allocatable :: troe_enabled(:)
    logical :: local_ok, global_ok, field_ok
    integer :: ierr, species_index, comparison_index, reaction_index
    integer :: value_offset
    integer :: reaction_stride

    nspecies = size(species)
    nreactions = size(reactions)
    sizes = [nspecies, nreactions, size(transport), size(composition)]
    local_ok = nspecies > 0 .and. nreactions > 0
    if (local_ok) local_ok = size(transport) == nspecies
    if (local_ok) local_ok = size(composition) == nspecies
    if (local_ok) local_ok = len_trim(bundle_sha256) == 64
    if (local_ok) local_ok = is_sha256(bundle_sha256)
    if (local_ok) then
      local_ok = trim(chemistry_integrator) == "explicit" .or. &
        trim(chemistry_integrator) == "implicit"
    end if
    if (local_ok) local_ok = all(ieee_is_finite( &
      [relative_tolerance, absolute_tolerance]))
    if (local_ok) local_ok = relative_tolerance > 0.0_dp .and. &
      absolute_tolerance > 0.0_dp
    if (local_ok) local_ok = all(ieee_is_finite(composition))
    if (local_ok) local_ok = minval(composition) >= 0.0_dp
    if (local_ok) &
      local_ok = abs(sum(composition) - 1.0_dp) <= 5.0e-10_dp
    if (local_ok) local_ok = compatible_transport_database(species, transport)
    if (local_ok) then
      do species_index = 1, nspecies
        if (.not. valid_nasa7_species(species(species_index))) then
          local_ok = .false.
          exit
        end if
        if (len_trim(species(species_index)%name) == 0) then
          local_ok = .false.
          exit
        end if
      end do
    end if
    if (local_ok) then
      do species_index = 1, nspecies
        do comparison_index = species_index + 1, nspecies
          if (trim(species(species_index)%name) == &
              trim(species(comparison_index)%name)) then
            local_ok = .false.
            exit
          end if
        end do
        if (.not. local_ok) exit
      end do
    end if
    if (local_ok) then
      do reaction_index = 1, nreactions
        if (.not. valid_elementary_reaction( &
            reactions(reaction_index), nspecies)) then
          local_ok = .false.
          exit
        end if
        if (allocated( &
            reactions(reaction_index)%third_body_efficiencies)) then
          if (size( &
              reactions(reaction_index)%third_body_efficiencies) /= &
              nspecies) then
            local_ok = .false.
            exit
          end if
        end if
      end do
    end if
    call MPI_Allreduce( &
      local_ok, global_ok, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    if (.not. global_ok) then
      ok = .false.
      return
    end if
    call consensus_integer_vector(comm, sizes, field_ok)
    if (.not. field_ok) then
      ok = .false.
      return
    end if

    allocate(species_names(nspecies), transport_names(nspecies))
    allocate(reaction_equations(nreactions))
    allocate(transport_geometries(nspecies), reaction_kinds(nreactions))
    allocate(reversible(nreactions), efficiencies_allocated(nreactions))
    allocate(troe_enabled(nreactions))
    allocate(species_values(18 * nspecies))
    allocate(transport_values(5 * nspecies))
    reaction_stride = 3 * nspecies + 13
    allocate(reaction_values(reaction_stride * nreactions))
    species_values = 0.0_dp
    transport_values = 0.0_dp
    reaction_values = 0.0_dp

    bundle_values(1) = bundle_sha256
    integrator_values(1) = chemistry_integrator
    chemistry_controls = [relative_tolerance, absolute_tolerance]
    do species_index = 1, nspecies
      species_names(species_index) = species(species_index)%name
      transport_names(species_index) = transport(species_index)%name
      transport_geometries(species_index) = transport(species_index)%geometry
      value_offset = 18 * (species_index - 1)
      species_values(value_offset + 1:value_offset + 4) = [ &
        species(species_index)%molecular_weight, &
        species(species_index)%temperature_min, &
        species(species_index)%temperature_mid, &
        species(species_index)%temperature_max]
      species_values(value_offset + 5:value_offset + 11) = &
        species(species_index)%low_coefficients
      species_values(value_offset + 12:value_offset + 18) = &
        species(species_index)%high_coefficients
      value_offset = 5 * (species_index - 1)
      transport_values(value_offset + 1:value_offset + 5) = [ &
        transport(species_index)%well_depth, &
        transport(species_index)%diameter, &
        transport(species_index)%dipole, &
        transport(species_index)%polarizability, &
        transport(species_index)%rotational_relaxation]
    end do
    do reaction_index = 1, nreactions
      reaction_equations(reaction_index) = reactions(reaction_index)%equation
      reaction_kinds(reaction_index) = reactions(reaction_index)%kind
      reversible(reaction_index) = reactions(reaction_index)%reversible
      efficiencies_allocated(reaction_index) = &
        allocated(reactions(reaction_index)%third_body_efficiencies)
      troe_enabled(reaction_index) = reactions(reaction_index)%troe%enabled
      value_offset = reaction_stride * (reaction_index - 1)
      reaction_values( &
        value_offset + 1:value_offset + nspecies) = &
        reactions(reaction_index)%reactant_stoich
      reaction_values( &
        value_offset + nspecies + 1:value_offset + 2 * nspecies) = &
        reactions(reaction_index)%product_stoich
      if (efficiencies_allocated(reaction_index)) then
        reaction_values( &
          value_offset + 2 * nspecies + 1:value_offset + 3 * nspecies) = &
          reactions(reaction_index)%third_body_efficiencies
      end if
      value_offset = value_offset + 3 * nspecies
      reaction_values(value_offset + 1:value_offset + 3) = [ &
        reactions(reaction_index)%forward_rate%pre_exponential, &
        reactions(reaction_index)%forward_rate%temperature_exponent, &
        reactions(reaction_index)%forward_rate%activation_energy]
      reaction_values(value_offset + 4:value_offset + 6) = [ &
        reactions(reaction_index)%low_pressure_rate%pre_exponential, &
        reactions(reaction_index)%low_pressure_rate%temperature_exponent, &
        reactions(reaction_index)%low_pressure_rate%activation_energy]
      reaction_values(value_offset + 7:value_offset + 9) = [ &
        reactions(reaction_index)%high_pressure_rate%pre_exponential, &
        reactions(reaction_index)%high_pressure_rate%temperature_exponent, &
        reactions(reaction_index)%high_pressure_rate%activation_energy]
      reaction_values(value_offset + 10:value_offset + 13) = [ &
        reactions(reaction_index)%troe%alpha, &
        reactions(reaction_index)%troe%temperature_3, &
        reactions(reaction_index)%troe%temperature_1, &
        reactions(reaction_index)%troe%temperature_2]
    end do
    local_ok = all(ieee_is_finite(species_values))
    if (local_ok) local_ok = all(ieee_is_finite(transport_values))
    if (local_ok) local_ok = all(ieee_is_finite(reaction_values))
    call MPI_Allreduce( &
      local_ok, global_ok, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    if (.not. global_ok) then
      ok = .false.
      return
    end if

    ok = .true.
    call consensus_character_vector(comm, bundle_values, field_ok)
    ok = ok .and. field_ok
    call consensus_character_vector(comm, integrator_values, field_ok)
    ok = ok .and. field_ok
    call consensus_real_vector(comm, chemistry_controls, field_ok)
    ok = ok .and. field_ok
    call consensus_real_vector(comm, composition, field_ok)
    ok = ok .and. field_ok
    call consensus_character_vector(comm, species_names, field_ok)
    ok = ok .and. field_ok
    call consensus_real_vector(comm, species_values, field_ok)
    ok = ok .and. field_ok
    call consensus_character_vector(comm, transport_names, field_ok)
    ok = ok .and. field_ok
    call consensus_integer_vector(comm, transport_geometries, field_ok)
    ok = ok .and. field_ok
    call consensus_real_vector(comm, transport_values, field_ok)
    ok = ok .and. field_ok
    call consensus_character_vector(comm, reaction_equations, field_ok)
    ok = ok .and. field_ok
    call consensus_integer_vector(comm, reaction_kinds, field_ok)
    ok = ok .and. field_ok
    call consensus_logical_vector(comm, reversible, field_ok)
    ok = ok .and. field_ok
    call consensus_logical_vector(comm, efficiencies_allocated, field_ok)
    ok = ok .and. field_ok
    call consensus_logical_vector(comm, troe_enabled, field_ok)
    ok = ok .and. field_ok
    call consensus_real_vector(comm, reaction_values, field_ok)
    ok = ok .and. field_ok
  end subroutine validate_mpi_amr_selected_context_3d

  subroutine consensus_character_vector(comm, values, ok)
    type(MPI_Comm), intent(in) :: comm
    character(len=*), intent(in) :: values(:)
    logical, intent(out) :: ok

    character(len=len(values(1))), allocatable :: root_values(:)
    logical :: local_ok, global_ok
    integer :: ierr, rank

    allocate(root_values(size(values)))
    root_values = ""
    call MPI_Comm_rank(comm, rank, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    if (rank == 0) root_values = values
    call MPI_Bcast( &
      root_values, len(values(1)) * size(values), MPI_CHARACTER, 0, comm, ierr)
    local_ok = ierr == MPI_SUCCESS
    if (local_ok) local_ok = all(root_values == values)
    call MPI_Allreduce( &
      local_ok, global_ok, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
    else
      ok = global_ok
    end if
  end subroutine consensus_character_vector

  subroutine consensus_integer_vector(comm, values, ok)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: values(:)
    logical, intent(out) :: ok

    integer, allocatable :: root_values(:)
    logical :: local_ok, global_ok
    integer :: ierr, rank

    allocate(root_values(size(values)))
    root_values = 0
    call MPI_Comm_rank(comm, rank, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    if (rank == 0) root_values = values
    call MPI_Bcast( &
      root_values, size(root_values), MPI_INTEGER, 0, comm, ierr)
    local_ok = ierr == MPI_SUCCESS
    if (local_ok) local_ok = all(root_values == values)
    call MPI_Allreduce( &
      local_ok, global_ok, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
    else
      ok = global_ok
    end if
  end subroutine consensus_integer_vector

  subroutine consensus_logical_vector(comm, values, ok)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(in) :: values(:)
    logical, intent(out) :: ok

    logical, allocatable :: root_values(:)
    logical :: local_ok, global_ok
    integer :: ierr, rank

    allocate(root_values(size(values)))
    root_values = .false.
    call MPI_Comm_rank(comm, rank, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    if (rank == 0) root_values = values
    call MPI_Bcast( &
      root_values, size(root_values), MPI_LOGICAL, 0, comm, ierr)
    local_ok = ierr == MPI_SUCCESS
    if (local_ok) local_ok = all(root_values .eqv. values)
    call MPI_Allreduce( &
      local_ok, global_ok, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
    else
      ok = global_ok
    end if
  end subroutine consensus_logical_vector

  subroutine consensus_real_vector(comm, values, ok)
    type(MPI_Comm), intent(in) :: comm
    real(dp), intent(in) :: values(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: root_values(:)
    integer(int64), allocatable :: root_bits(:), value_bits(:)
    logical :: local_ok, global_ok
    integer :: ierr, rank

    allocate(root_values(size(values)))
    allocate(root_bits(size(values)), value_bits(size(values)))
    root_values = 0.0_dp
    call MPI_Comm_rank(comm, rank, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    if (rank == 0) root_values = values
    call MPI_Bcast( &
      root_values, size(root_values), MPI_DOUBLE_PRECISION, 0, comm, ierr)
    local_ok = ierr == MPI_SUCCESS
    if (local_ok) then
      root_bits = transfer(root_values, root_bits)
      value_bits = transfer(values, value_bits)
      local_ok = all(root_bits == value_bits)
    end if
    call MPI_Allreduce( &
      local_ok, global_ok, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
    else
      ok = global_ok
    end if
  end subroutine consensus_real_vector

  pure logical function is_sha256(value) result(valid)
    character(len=*), intent(in) :: value

    integer :: index

    valid = len_trim(value) == 64
    if (.not. valid) return
    do index = 1, 64
      select case (value(index:index))
      case ('0':'9', 'a':'f', 'A':'F')
      case default
        valid = .false.
        return
      end select
    end do
  end function is_sha256

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
    if (ok) ok = minimum_density > 0.0_dp
    if (ok) ok = minimum_pressure > 0.0_dp
    if (ok) ok = minimum_temperature > 0.0_dp
    if (ok) ok = maximum_density > 0.0_dp
    if (ok) ok = maximum_pressure > 0.0_dp
    if (ok) ok = maximum_temperature > 0.0_dp
    if (ok) ok = maximum_speed >= 0.0_dp
  end subroutine level_extrema

end module mpi_amr_reactive_3d_application_mod
