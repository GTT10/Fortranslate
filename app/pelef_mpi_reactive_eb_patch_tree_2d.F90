program pelef_mpi_reactive_eb_patch_tree_2d
  use, intrinsic :: iso_fortran_env, only: error_unit
  use mpi_f08
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport, &
    load_h2o2_full_transport
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use reactive_1d_mod, only: reactive_nvar
  use simulation_config_reactive_eb_amr_2d_mod, only: &
    reactive_eb_amr_2d_config, read_reactive_eb_amr_2d_configuration
  use reactive_2d_mod, only: initialize_reactive_2d
  use reactive_boundary_2d_mod, only: reactive_boundary_set_2d
  use reactive_eb_2d_driver_mod, only: &
    build_configured_eb_geometry_2d, build_configured_eb_geometry_region_2d, &
    build_configured_reactive_boundary_set_2d
  use eb_geometry_2d_mod, only: eb_geometry_2d
  use amr_eb_regrid_2d_mod, only: amr_eb_tagging_criteria_2d
  use amr_eb_patch_tree_2d_mod, only: &
    amr_eb_patch_tree_level_plan_2d, amr_eb_patch_tree_topology_2d, &
    initialize_amr_eb_patch_tree_topology_2d
  use amr_eb_patch_tree_reactive_2d_mod, only: &
    reactive_amr_eb_patch_tree_checkpoint_fingerprint_2d
  use reactive_eb_amr_2d_driver_mod, only: &
    build_reactive_amr_eb_patch_tree_checkpoint_fingerprint_2d
  use mpi_amr_eb_patch_tree_2d_mod, only: &
    mpi_amr_eb_patch_tree_distribution_2d, &
    mpi_sparse_reactive_amr_eb_patch_tree_2d, &
    initialize_mpi_amr_eb_patch_tree_distribution_2d, &
    initialize_sparse_owned_reactive_amr_eb_patch_tree_root_2d, &
    regrid_tagged_sparse_owned_reactive_amr_eb_patch_tree_2d, &
    compute_sparse_owned_reactive_amr_eb_patch_tree_timestep_2d, &
    advance_sparse_owned_reactive_amr_eb_patch_tree_full_physics_2d, &
    composite_sparse_amr_eb_patch_tree_integral_2d
  use mpi_amr_eb_patch_tree_io_2d_mod, only: &
    write_sparse_owned_reactive_amr_eb_patch_tree_2d_checkpoint, &
    read_sparse_owned_reactive_amr_eb_patch_tree_2d_checkpoint, &
    write_sparse_owned_reactive_amr_eb_patch_tree_2d_csv
  implicit none

  integer, parameter :: io_root = 0
  type(reactive_eb_amr_2d_config) :: config
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(reactive_boundary_set_2d) :: boundaries
  type(eb_geometry_2d) :: root_geometry
  type(amr_eb_patch_tree_level_plan_2d), allocatable :: empty_plans(:)
  type(amr_eb_patch_tree_topology_2d) :: topology
  type(reactive_amr_eb_patch_tree_checkpoint_fingerprint_2d) :: fingerprint
  type(mpi_amr_eb_patch_tree_distribution_2d) :: distribution
  type(mpi_amr_eb_patch_tree_distribution_2d) :: new_distribution
  type(mpi_sparse_reactive_amr_eb_patch_tree_2d) :: sparse
  type(amr_eb_tagging_criteria_2d) :: criteria
  real(dp), allocatable :: root_state(:, :, :), root_temperature(:, :)
  real(dp), allocatable :: initial_integrals(:), final_integrals(:)
  integer, allocatable :: chemistry_level_advances(:)
  integer, allocatable :: transport_level_advances(:)
  integer, allocatable :: hydro_level_advances(:)
  integer, allocatable :: local_step_chemistry_advances(:)
  integer, allocatable :: local_step_transport_advances(:)
  integer, allocatable :: local_step_hydro_advances(:)
  integer, allocatable :: global_step_chemistry_advances(:)
  integer, allocatable :: global_step_transport_advances(:)
  integer, allocatable :: global_step_hydro_advances(:)
  real(dp) :: base_density, conservation_error, dt, dx, dy
  real(dp) :: minimum_dt, minimum_transport_theta, remaining
  real(dp) :: step_theta, time, time_tolerance
  character(len=1024) :: input_path, message, output_path
  character(len=160) :: physics_context
  integer :: argument_count, ierr, last_checkpoint_step, nranks, rank
  integer :: local_root_initializers, root_initializer_ranks
  integer :: cumulative_tagged_cells, regrid_evaluations
  integer :: regrids, steps, tagged_cells, transferred_cells
  logical :: barodiffusion_active, changed, history_ok, ok, restart_run
  logical :: species_diffusion_active, stopped_after_checkpoint
  logical :: thermal_conduction_active, viscosity_active

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Init failed"
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Comm_rank failed"
  call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)
  if (ierr /= MPI_SUCCESS) call abort_run("MPI_Comm_size failed", 2)

  argument_count = command_argument_count()
  if (argument_count < 1 .or. argument_count > 2) call abort_run( &
    "Usage: pelef_mpi_reactive_eb_patch_tree_2d <input.nml> [output.csv]", 2)
  call get_command_argument(1, input_path)
  call read_reactive_eb_amr_2d_configuration( &
    trim(input_path), config, ok, message)
  if (.not. ok) call abort_run(trim(message), 2)
  viscosity_active = config%eb%flow%transport_enabled .and. &
    config%eb%flow%viscosity_enabled
  thermal_conduction_active = config%eb%flow%transport_enabled .and. &
    config%eb%flow%thermal_conduction_enabled
  species_diffusion_active = config%eb%flow%transport_enabled .and. &
    config%eb%flow%species_diffusion_enabled
  barodiffusion_active = species_diffusion_active .and. &
    config%eb%flow%barodiffusion_enabled
  if (config%three_level_enabled .or. config%multipatch_enabled) &
    call abort_run("Sparse MPI patch tree excludes fixed-depth modes", 2)
  call build_reactive_amr_eb_patch_tree_checkpoint_fingerprint_2d( &
    config, fingerprint, ok)
  if (.not. ok) call abort_run("Checkpoint fingerprint failed", 2)
  output_path = config%eb%flow%output_file
  if (argument_count == 2) call get_command_argument(2, output_path)
  if (len_trim(output_path) == 0) call abort_run("Output path is empty", 2)

  select case (trim(config%eb%flow%chemistry_model))
  case ("elementary")
    call load_h2o2_elementary_thermo(species, ok)
    if (.not. ok) call abort_run("Failed to load elementary thermodynamics", 3)
    call load_h2o2_elementary_mechanism(reactions, ok)
    if (.not. ok) call abort_run("Failed to load elementary mechanism", 3)
    call load_h2o2_elementary_transport(transport, ok)
    if (.not. ok) call abort_run("Failed to load elementary transport", 3)
  case ("full_h2o2")
    call load_h2o2_full_thermo(species, ok)
    if (.not. ok) call abort_run("Failed to load full H2/O2 thermodynamics", 3)
    call load_h2o2_full_mechanism(reactions, ok)
    if (.not. ok) call abort_run("Failed to load full H2/O2 mechanism", 3)
    call load_h2o2_full_transport(transport, ok)
    if (.not. ok) call abort_run("Failed to load full H2/O2 transport", 3)
  case default
    call abort_run("Unknown chemistry model", 3)
  end select
  call build_configured_reactive_boundary_set_2d( &
    species, config%eb, boundaries, ok)
  if (.not. ok) call abort_run("Boundary initialization failed", 3)

  criteria%relative_gradient_threshold = &
    config%regrid_relative_temperature_gradient
  criteria%absolute_gradient_threshold = &
    config%regrid_absolute_temperature_gradient
  criteria%scale_floor = config%regrid_temperature_scale_floor
  criteria%buffer_cells = config%regrid_buffer_cells
  criteria%minimum_patch_cells_x = config%regrid_minimum_patch_cells_x
  criteria%minimum_patch_cells_y = config%regrid_minimum_patch_cells_y
  criteria%maximum_patch_gap_cells = config%regrid_maximum_patch_gap_cells

  time = 0.0_dp
  steps = 0
  regrids = 0
  regrid_evaluations = 0
  cumulative_tagged_cells = 0
  minimum_dt = 0.0_dp
  minimum_transport_theta = 1.0_dp
  restart_run = len_trim(config%restart_file) > 0
  if (restart_run) then
    call read_sparse_owned_reactive_amr_eb_patch_tree_2d_checkpoint( &
      config%restart_file, species, MPI_COMM_WORLD, io_root, &
      config%patch_tree_maximum_levels, config%patch_tree_mpi_work_exponent, &
      distribution, sparse, time, steps, regrids, minimum_dt, ok, &
      fingerprint=fingerprint, &
      minimum_transport_theta=minimum_transport_theta, &
      initial_integrals=initial_integrals, &
      chemistry_level_advances=chemistry_level_advances, &
      transport_level_advances=transport_level_advances, &
      hydro_level_advances=hydro_level_advances, &
      regrid_evaluations=regrid_evaluations, &
      cumulative_tagged_cells=cumulative_tagged_cells)
    if (.not. ok) call abort_run("Sparse patch-tree restart failed", 4)
    if (size(chemistry_level_advances) /= &
          config%patch_tree_maximum_levels .or. &
        size(transport_level_advances) /= &
          config%patch_tree_maximum_levels .or. &
        size(hydro_level_advances) /= &
          config%patch_tree_maximum_levels) &
      call abort_run("Sparse patch-tree counter capacity mismatch", 4)
  else
    call build_configured_eb_geometry_2d(config%eb, root_geometry, ok)
    if (.not. ok) call abort_run("Root EB geometry failed", 4)
    allocate(empty_plans(0))
    call initialize_amr_eb_patch_tree_topology_2d( &
      root_geometry, empty_plans, topology, ok)
    if (.not. ok) call abort_run("Root patch-tree topology failed", 4)
    call initialize_mpi_amr_eb_patch_tree_distribution_2d( &
      topology, MPI_COMM_WORLD, distribution, ok, &
      config%patch_tree_mpi_work_exponent)
    if (.not. ok) call abort_run("Sparse patch-tree ownership failed", 4)
    local_root_initializers = 0
    if (distribution%is_local(0, 1)) then
      call initialize_reactive_2d( &
        species, config%eb%flow, root_state, root_temperature, dx, dy, &
        base_density, ok)
      if (.not. ok) call abort_run("Root reactive state failed", 4)
      local_root_initializers = 1
    end if
    call MPI_Allreduce( &
      local_root_initializers, root_initializer_ranks, 1, MPI_INTEGER, &
      MPI_SUM, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) &
      call abort_run("Root reactive state ownership reduction failed", 4)
    if (root_initializer_ranks /= 1) &
      call abort_run("Root reactive state ownership failed", 4)
    call initialize_sparse_owned_reactive_amr_eb_patch_tree_root_2d( &
      distribution, topology, reactive_nvar(size(species)), root_state, &
      root_temperature, sparse, ok)
    if (.not. ok) call abort_run("Sparse patch-tree initialization failed", 4)
    if (allocated(root_state) .or. allocated(root_temperature)) &
      call abort_run("Root reactive state transfer failed", 4)
    if (config%dynamic_regridding .and. &
        config%regrid_at_initialization) then
      call regrid_tagged_sparse_owned_reactive_amr_eb_patch_tree_2d( &
        species, distribution, sparse, criteria, &
        config%patch_tree_maximum_levels, config%refinement_ratio, &
        build_patch_tree_geometry, new_distribution, ok, changed, &
        tagged_cells, transferred_cells, &
        prolongation_method=config%prolongation_method)
      if (.not. ok) call abort_run("Initial sparse regrid failed", 4)
      call accumulate_regrid_history(tagged_cells, history_ok)
      if (.not. history_ok) call abort_run("Regrid history overflow", 4)
      distribution = new_distribution
      if (changed) regrids = regrids + 1
    end if
    minimum_dt = huge(1.0_dp)
    allocate(initial_integrals(sparse%nvar))
    call composite_sparse_amr_eb_patch_tree_integral_2d( &
      distribution, sparse, initial_integrals, ok)
    if (.not. ok) call abort_run("Initial sparse integral failed", 4)
    allocate(chemistry_level_advances( &
      config%patch_tree_maximum_levels), source=0)
    allocate(transport_level_advances( &
      config%patch_tree_maximum_levels), source=0)
    allocate(hydro_level_advances( &
      config%patch_tree_maximum_levels), source=0)
  end if

  time_tolerance = 16.0_dp * epsilon(1.0_dp) * &
    max(tiny(1.0_dp), abs(config%eb%flow%final_time))
  if (time > config%eb%flow%final_time + time_tolerance) &
    call abort_run("Restart time exceeds configured final time", 4)
  allocate(final_integrals(sparse%nvar))

  stopped_after_checkpoint = .false.
  last_checkpoint_step = -1
  do
    remaining = config%eb%flow%final_time - time
    if (remaining <= time_tolerance) exit
    if (steps >= config%eb%flow%maximum_steps) &
      call abort_run("Sparse patch-tree step limit reached", 5)
    call compute_sparse_owned_reactive_amr_eb_patch_tree_timestep_2d( &
      species, transport, distribution, sparse, config%eb%flow%cfl, &
      config%eb%flow%transport_cfl, viscosity_active, &
      thermal_conduction_active, species_diffusion_active, dt, ok)
    if (.not. ok) call abort_run("Sparse patch-tree timestep failed", 5)
    dt = min(dt, remaining)
    if (allocated(local_step_chemistry_advances)) &
      deallocate(local_step_chemistry_advances)
    if (allocated(local_step_transport_advances)) &
      deallocate(local_step_transport_advances)
    if (allocated(local_step_hydro_advances)) &
      deallocate(local_step_hydro_advances)
    if (allocated(global_step_chemistry_advances)) &
      deallocate(global_step_chemistry_advances)
    if (allocated(global_step_transport_advances)) &
      deallocate(global_step_transport_advances)
    if (allocated(global_step_hydro_advances)) &
      deallocate(global_step_hydro_advances)
    allocate(local_step_chemistry_advances(sparse%level_count()), source=0)
    allocate(local_step_transport_advances(sparse%level_count()), source=0)
    allocate(local_step_hydro_advances(sparse%level_count()), source=0)
    allocate(global_step_chemistry_advances(sparse%level_count()), source=0)
    allocate(global_step_transport_advances(sparse%level_count()), source=0)
    allocate(global_step_hydro_advances(sparse%level_count()), source=0)
    call advance_sparse_owned_reactive_amr_eb_patch_tree_full_physics_2d( &
      species, reactions, transport, distribution, sparse, &
      config%eb%flow%riemann_solver, config%eb%flow%reconstruction, &
      config%eb%flow%limiter, config%eb%state_redist_max_order, dt, &
      config%eb%flow%chemistry_enabled, &
      config%eb%flow%chemistry_relative_tolerance, &
      config%eb%flow%chemistry_absolute_tolerance, &
      viscosity_active, thermal_conduction_active, &
      species_diffusion_active, barodiffusion_active, boundaries, &
      config%eb%state_redist_target_volume_fraction, step_theta, ok, &
      physics_context, local_step_chemistry_advances, &
      local_step_transport_advances, local_step_hydro_advances)
    if (.not. ok) call abort_run( &
      "Sparse full physics failed: " // trim(physics_context), 5)
    call MPI_Allreduce( &
      local_step_chemistry_advances, global_step_chemistry_advances, &
      size(local_step_chemistry_advances), MPI_INTEGER, MPI_SUM, &
      MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) &
      call abort_run("Sparse chemistry counter reduction failed", 5)
    call MPI_Allreduce( &
      local_step_transport_advances, global_step_transport_advances, &
      size(local_step_transport_advances), MPI_INTEGER, MPI_SUM, &
      MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) &
      call abort_run("Sparse transport counter reduction failed", 5)
    call MPI_Allreduce( &
      local_step_hydro_advances, global_step_hydro_advances, &
      size(local_step_hydro_advances), MPI_INTEGER, MPI_SUM, &
      MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) &
      call abort_run("Sparse hydro counter reduction failed", 5)
    if (any(global_step_chemistry_advances > &
          huge(1) - chemistry_level_advances( &
            1:size(global_step_chemistry_advances))) .or. &
        any(global_step_transport_advances > &
          huge(1) - transport_level_advances( &
            1:size(global_step_transport_advances))) .or. &
        any(global_step_hydro_advances > &
          huge(1) - hydro_level_advances( &
            1:size(global_step_hydro_advances)))) &
      call abort_run("Sparse operator counter overflow", 5)
    chemistry_level_advances(1:size(global_step_chemistry_advances)) = &
      chemistry_level_advances(1:size(global_step_chemistry_advances)) + &
      global_step_chemistry_advances
    transport_level_advances(1:size(global_step_transport_advances)) = &
      transport_level_advances(1:size(global_step_transport_advances)) + &
      global_step_transport_advances
    hydro_level_advances(1:size(global_step_hydro_advances)) = &
      hydro_level_advances(1:size(global_step_hydro_advances)) + &
      global_step_hydro_advances
    time = time + dt
    minimum_dt = min(minimum_dt, dt)
    minimum_transport_theta = min(minimum_transport_theta, step_theta)
    steps = steps + 1

    if (config%dynamic_regridding .and. &
        modulo(steps, config%regrid_interval) == 0) then
      call regrid_tagged_sparse_owned_reactive_amr_eb_patch_tree_2d( &
        species, distribution, sparse, criteria, &
        config%patch_tree_maximum_levels, config%refinement_ratio, &
        build_patch_tree_geometry, new_distribution, ok, changed, &
        tagged_cells, transferred_cells, &
        prolongation_method=config%prolongation_method)
      if (.not. ok) call abort_run("Periodic sparse regrid failed", 5)
      call accumulate_regrid_history(tagged_cells, history_ok)
      if (.not. history_ok) call abort_run("Regrid history overflow", 5)
      distribution = new_distribution
      if (changed) regrids = regrids + 1
    end if
    if (config%checkpoint_interval > 0) then
      if (modulo(steps, config%checkpoint_interval) == 0) then
        call write_sparse_checkpoint(ok)
        if (.not. ok) call abort_run("Sparse checkpoint write failed", 6)
        last_checkpoint_step = steps
        if (config%checkpoint_stop_after_write) then
          stopped_after_checkpoint = .true.
          exit
        end if
      end if
    end if
  end do
  if (.not. stopped_after_checkpoint) time = config%eb%flow%final_time

  if (len_trim(config%checkpoint_file) > 0 .and. &
      last_checkpoint_step /= steps) then
    call write_sparse_checkpoint(ok)
    if (.not. ok) call abort_run("Final sparse checkpoint failed", 6)
  end if
  call composite_sparse_amr_eb_patch_tree_integral_2d( &
    distribution, sparse, final_integrals, ok)
  if (.not. ok) call abort_run("Final sparse integral failed", 6)
  conservation_error = maxval(abs(final_integrals - initial_integrals) / &
    max(1.0_dp, abs(initial_integrals)))
  call write_sparse_owned_reactive_amr_eb_patch_tree_2d_csv( &
    trim(output_path), species, distribution, sparse, io_root, time, ok)
  if (.not. ok) call abort_run("Sparse composite output failed", 6)

  if (rank == io_root) then
    write(*, '(a)') "PeleF " // pelef_version // &
      " sparse MPI reactive EB patch-tree 2D"
    write(*, '(a,i0)') "MPI ranks: ", nranks
    write(*, '(a,i0)') "Levels: ", sparse%level_count()
    write(*, '(a,i0)') "Patches: ", sum(distribution%rank_patch_counts)
    write(*, '(a,i0)') "Completed root steps: ", steps
    write(*, '(a,i0)') "Completed regrids: ", regrids
    write(*, '(a,i0)') "Regrid evaluations: ", regrid_evaluations
    write(*, '(a,i0)') "Cumulative tagged cells: ", cumulative_tagged_cells
    write(*, '(a,i0)') "MPI work exponent: ", &
      config%patch_tree_mpi_work_exponent
    if (.not. restart_run) write(*, '(a,i0)') &
      "Fresh root initializer ranks: ", root_initializer_ranks
    write(*, '(a,l2)') "Restarted: ", restart_run
    write(*, '(a,l2)') "Stopped after checkpoint: ", &
      stopped_after_checkpoint
    write(*, '(a,es24.16)') "Final time: ", time
    write(*, '(a,es24.16)') "Minimum accepted root dt: ", minimum_dt
    write(*, '(a,es24.16)') "Minimum transport limiter theta: ", &
      minimum_transport_theta
    write(*, '(a,es24.16)') "Maximum composite conservation error: ", &
      conservation_error
    write(*, '(a,*(1x,i0))') "Chemistry level advances:", &
      chemistry_level_advances
    write(*, '(a,*(1x,i0))') "Transport level advances:", &
      transport_level_advances
    write(*, '(a,*(1x,i0))') "Hydro level advances:", &
      hydro_level_advances
    write(*, '(a,1x,a)') "Composite output:", trim(output_path)
  end if

  call MPI_Finalize(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Finalize failed"

contains

  subroutine accumulate_regrid_history(tagged, history_ok)
    integer, intent(in) :: tagged
    logical, intent(out) :: history_ok

    history_ok = tagged >= 0 .and. &
      regrid_evaluations < huge(regrid_evaluations) .and. &
      tagged <= huge(cumulative_tagged_cells) - cumulative_tagged_cells
    if (.not. history_ok) return
    regrid_evaluations = regrid_evaluations + 1
    cumulative_tagged_cells = cumulative_tagged_cells + tagged
  end subroutine accumulate_regrid_history

  subroutine build_patch_tree_geometry( &
      parent_geometry, coarse_i_lower, coarse_i_upper, coarse_j_lower, &
      coarse_j_upper, refinement_ratio, child_geometry, geometry_ok)
    type(eb_geometry_2d), intent(in) :: parent_geometry
    integer, intent(in) :: coarse_i_lower, coarse_i_upper
    integer, intent(in) :: coarse_j_lower, coarse_j_upper, refinement_ratio
    type(eb_geometry_2d), intent(out) :: child_geometry
    logical, intent(out) :: geometry_ok

    real(dp) :: x_lower, x_upper, y_lower, y_upper
    integer :: nx, ny

    geometry_ok = coarse_i_lower >= 1 .and. &
      coarse_i_upper <= parent_geometry%nx .and. &
      coarse_j_lower >= 1 .and. coarse_j_upper <= parent_geometry%ny .and. &
      coarse_i_upper >= coarse_i_lower .and. &
      coarse_j_upper >= coarse_j_lower .and. refinement_ratio >= 2
    if (.not. geometry_ok) return
    nx = (coarse_i_upper - coarse_i_lower + 1) * refinement_ratio
    ny = (coarse_j_upper - coarse_j_lower + 1) * refinement_ratio
    x_lower = parent_geometry%x_lower + &
      real(coarse_i_lower - 1, dp) * parent_geometry%dx
    x_upper = parent_geometry%x_lower + &
      real(coarse_i_upper, dp) * parent_geometry%dx
    y_lower = parent_geometry%y_lower + &
      real(coarse_j_lower - 1, dp) * parent_geometry%dy
    y_upper = parent_geometry%y_lower + &
      real(coarse_j_upper, dp) * parent_geometry%dy
    call build_configured_eb_geometry_region_2d( &
      config%eb, nx, ny, x_lower, x_upper, y_lower, y_upper, &
      child_geometry, geometry_ok)
  end subroutine build_patch_tree_geometry

  subroutine write_sparse_checkpoint(checkpoint_ok)
    logical, intent(out) :: checkpoint_ok

    call write_sparse_owned_reactive_amr_eb_patch_tree_2d_checkpoint( &
      config%checkpoint_file, species, distribution, sparse, io_root, time, &
      steps, regrids, minimum_dt, checkpoint_ok, fingerprint=fingerprint, &
      minimum_transport_theta=minimum_transport_theta, &
      initial_integrals=initial_integrals, &
      chemistry_level_advances=chemistry_level_advances, &
      transport_level_advances=transport_level_advances, &
      hydro_level_advances=hydro_level_advances, &
      regrid_evaluations=regrid_evaluations, &
      cumulative_tagged_cells=cumulative_tagged_cells)
  end subroutine write_sparse_checkpoint

  subroutine abort_run(reason, code)
    character(len=*), intent(in) :: reason
    integer, intent(in) :: code

    integer :: abort_ierr

    if (rank == io_root) write(error_unit, '(a)') trim(reason)
    call MPI_Abort(MPI_COMM_WORLD, code, abort_ierr)
    error stop code
  end subroutine abort_run

end program pelef_mpi_reactive_eb_patch_tree_2d
