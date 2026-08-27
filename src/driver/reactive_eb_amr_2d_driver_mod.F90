module reactive_eb_amr_2d_driver_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use state_indices_mod, only: irho
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: gas_transport_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use reactive_2d_mod, only: &
    initialize_reactive_2d, advance_reactive_chemistry_2d
  use reactive_boundary_2d_mod, only: &
    reactive_boundary_set_2d, build_reactive_boundary_set_2d
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_covered_cell
  use simulation_config_reactive_eb_amr_2d_mod, only: &
    reactive_eb_amr_2d_config
  use reactive_eb_2d_driver_mod, only: &
    build_configured_eb_geometry_2d, &
    build_configured_eb_geometry_region_2d, &
    compute_reactive_eb_cfl_timestep_2d, reactive_eb_integrals_2d, &
    advance_reactive_eb_strang_2d
  use eb_reactive_transport_2d_mod, only: reactive_eb_transport_timestep_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, build_amr_eb_patch_2d, composite_eb_integral_2d, &
    average_down_reactive_eb_state_patch_2d
  use amr_eb_reactive_2d_mod, only: &
    prolong_reactive_eb_patch_pcm_2d, &
    advance_two_level_reactive_eb_hydro_2d
  use amr_eb_transport_2d_mod, only: &
    advance_two_level_reactive_eb_transport_2d
  use amr_eb_multilevel_2d_mod, only: &
    average_down_three_level_reactive_eb_state_2d, &
    composite_three_level_eb_integral_2d
  use amr_eb_multilevel_reactive_2d_mod, only: &
    advance_three_level_reactive_eb_hydro_2d
  use amr_eb_multilevel_transport_2d_mod, only: &
    advance_three_level_reactive_eb_transport_2d
  use amr_eb_patch_tree_2d_mod, only: &
    amr_eb_patch_tree_level_plan_2d, amr_eb_patch_tree_topology_2d, &
    initialize_amr_eb_patch_tree_topology_2d
  use amr_eb_patch_tree_reactive_2d_mod, only: &
    reactive_amr_eb_patch_tree_2d, &
    initialize_reactive_amr_eb_patch_tree_2d, &
    regrid_tagged_reactive_amr_eb_patch_tree_2d, &
    write_reactive_amr_eb_patch_tree_2d_checkpoint, &
    read_reactive_amr_eb_patch_tree_2d_checkpoint, &
    compute_reactive_amr_eb_patch_tree_timestep_2d, &
    advance_reactive_amr_eb_patch_tree_full_physics_2d, &
    composite_integral_reactive_amr_eb_patch_tree_2d
  use amr_eb_multipatch_transport_2d_mod, only: &
    advance_reactive_eb_patch_set_transport_2d
  use amr_eb_regrid_2d_mod, only: &
    amr_eb_tagging_criteria_2d, amr_eb_regrid_plan_2d, &
    amr_eb_regrid_plan_collection_2d, &
    reactive_eb_patch_set_2d, &
    tag_reactive_eb_temperature_gradient_2d, &
    build_amr_eb_regrid_plan_2d, &
    plan_reactive_eb_temperature_regrid_2d, &
    plan_reactive_eb_temperature_regrid_collection_2d, &
    initialize_reactive_eb_patch_set_2d, &
    average_down_reactive_eb_patch_set_2d, &
    composite_reactive_eb_patch_set_integral_2d, &
    regrid_reactive_eb_patch_set_2d, &
    advance_reactive_eb_patch_set_hydro_2d, &
    collapse_two_level_reactive_eb_patch_2d, &
    regrid_two_level_reactive_eb_patch_2d
  implicit none
  private

  character(len=*), parameter :: reactive_eb_amr_checkpoint_magic = &
    "PELEF_REACTIVE_EB_AMR_2D_CHECKPOINT"
  integer, parameter :: reactive_eb_amr_checkpoint_schema = 1
  character(len=*), parameter :: reactive_eb_patch_set_checkpoint_magic = &
    "PELEF_REACTIVE_EB_AMR_PATCH_SET_2D_CHECKPOINT"
  integer, parameter :: reactive_eb_patch_set_checkpoint_schema = 1
  character(len=*), parameter :: reactive_eb_three_level_checkpoint_magic = &
    "PELEF_REACTIVE_EB_AMR_THREE_LEVEL_2D_CHECKPOINT"
  integer, parameter :: reactive_eb_three_level_checkpoint_schema = 1
  character(len=*), parameter :: &
    reactive_eb_dynamic_three_level_checkpoint_magic = &
      "PELEF_REACTIVE_EB_AMR_DYNAMIC_THREE_LEVEL_2D_CHECKPOINT"
  integer, parameter :: reactive_eb_dynamic_three_level_checkpoint_schema = 1

  public :: compute_reactive_eb_amr_cfl_timestep_2d
  public :: compute_three_level_reactive_eb_cfl_timestep_2d
  public :: compute_reactive_eb_patch_set_cfl_timestep_2d
  public :: advance_two_level_reactive_eb_strang_2d
  public :: advance_three_level_reactive_eb_strang_2d
  public :: advance_reactive_eb_patch_set_strang_2d
  public :: regrid_reactive_eb_amr_hierarchy_2d
  public :: regrid_three_level_reactive_eb_amr_finest_2d
  public :: regrid_reactive_eb_amr_patch_set_2d
  public :: write_reactive_eb_amr_2d_checkpoint
  public :: read_reactive_eb_amr_2d_checkpoint
  public :: write_reactive_eb_amr_patch_set_2d_checkpoint
  public :: read_reactive_eb_amr_patch_set_2d_checkpoint
  public :: write_reactive_eb_amr_three_level_2d_checkpoint
  public :: read_reactive_eb_amr_three_level_2d_checkpoint
  public :: simulate_reactive_eb_amr_2d
  public :: simulate_reactive_eb_amr_patch_set_2d
  public :: simulate_three_level_reactive_eb_amr_2d
  public :: simulate_reactive_amr_eb_patch_tree_2d

contains

  pure logical function supported_reactive_eb_amr_config(config) &
      result(supported)
    type(reactive_eb_amr_2d_config), intent(in) :: config

    supported = config%eb%flow%nx >= 4 .and. &
      config%eb%flow%ny >= 4 .and. &
      config%eb%flow%maximum_steps >= 1 .and. &
      config%eb%flow%x_upper > config%eb%flow%x_lower .and. &
      config%eb%flow%y_upper > config%eb%flow%y_lower .and. &
      ieee_is_finite(config%eb%flow%final_time) .and. &
      config%eb%flow%final_time > 0.0_dp .and. &
      ieee_is_finite(config%eb%flow%cfl) .and. &
      config%eb%flow%cfl > 0.0_dp .and. config%eb%flow%cfl <= 0.8_dp .and. &
      (.not. config%eb%flow%transport_enabled .or. &
       (len_trim(config%checkpoint_file) == 0 .and. &
        len_trim(config%restart_file) == 0)) .and. &
      ieee_is_finite(config%eb%flow%chemistry_relative_tolerance) .and. &
      config%eb%flow%chemistry_relative_tolerance > 0.0_dp .and. &
      ieee_is_finite(config%eb%flow%chemistry_absolute_tolerance) .and. &
      config%eb%flow%chemistry_absolute_tolerance > 0.0_dp .and. &
      .not. config%eb%flow%use_transverse_correction .and. &
      (trim(config%eb%flow%reconstruction) == "pcm" .or. &
       trim(config%eb%flow%reconstruction) == "characteristic_plm") .and. &
      trim(config%eb%flow%boundary_x_lower) == "outflow" .and. &
      trim(config%eb%flow%boundary_x_upper) == "outflow" .and. &
      trim(config%eb%flow%boundary_y_lower) == "outflow" .and. &
      trim(config%eb%flow%boundary_y_upper) == "outflow" .and. &
      ieee_is_finite(config%eb%state_redist_target_volume_fraction) .and. &
      config%eb%state_redist_target_volume_fraction > 0.0_dp .and. &
      config%eb%state_redist_target_volume_fraction <= 1.0_dp .and. &
      (config%eb%state_redist_max_order == 0 .or. &
       config%eb%state_redist_max_order == 2) .and. &
      config%coarse_i_lower >= 1 .and. &
      config%coarse_i_upper <= config%eb%flow%nx .and. &
      config%coarse_j_lower >= 1 .and. &
      config%coarse_j_upper <= config%eb%flow%ny .and. &
      config%coarse_i_upper >= config%coarse_i_lower .and. &
      config%coarse_j_upper >= config%coarse_j_lower .and. &
      config%refinement_ratio >= 2 .and. config%regrid_interval >= 1 .and. &
      (.not. config%remove_fine_patch_when_untagged .or. &
       config%dynamic_regridding) .and. &
      ieee_is_finite(config%regrid_relative_temperature_gradient) .and. &
      config%regrid_relative_temperature_gradient >= 0.0_dp .and. &
      ieee_is_finite(config%regrid_absolute_temperature_gradient) .and. &
      config%regrid_absolute_temperature_gradient >= 0.0_dp .and. &
      ieee_is_finite(config%regrid_temperature_scale_floor) .and. &
      config%regrid_temperature_scale_floor > 0.0_dp .and. &
      config%regrid_buffer_cells >= 0 .and. &
      config%regrid_minimum_patch_cells_x >= 1 .and. &
      config%regrid_minimum_patch_cells_x <= config%eb%flow%nx .and. &
      config%regrid_minimum_patch_cells_y >= 1 .and. &
      config%regrid_minimum_patch_cells_y <= config%eb%flow%ny .and. &
      config%regrid_maximum_patch_gap_cells >= 0 .and. &
      (.not. config%multipatch_enabled .or. config%dynamic_regridding) .and. &
      config%checkpoint_interval >= 0 .and. &
      (config%checkpoint_interval == 0 .or. &
       len_trim(config%checkpoint_file) > 0) .and. &
      (.not. config%checkpoint_stop_after_write .or. &
       (config%checkpoint_interval > 0 .and. &
        len_trim(config%checkpoint_file) > 0))
  end function supported_reactive_eb_amr_config

  pure logical function supported_three_level_reactive_eb_amr_config( &
      config) result(supported)
    type(reactive_eb_amr_2d_config), intent(in) :: config

    integer :: level_one_nx, level_one_ny

    level_one_nx = (config%coarse_i_upper - config%coarse_i_lower + 1) * &
      config%refinement_ratio
    level_one_ny = (config%coarse_j_upper - config%coarse_j_lower + 1) * &
      config%refinement_ratio
    supported = supported_reactive_eb_amr_config(config) .and. &
      config%three_level_enabled .and. .not. config%multipatch_enabled .and. &
      config%level_two_i_lower >= 3 .and. &
      config%level_two_i_upper <= level_one_nx - 2 .and. &
      config%level_two_j_lower >= 3 .and. &
      config%level_two_j_upper <= level_one_ny - 2 .and. &
      config%level_two_i_upper >= config%level_two_i_lower .and. &
      config%level_two_j_upper >= config%level_two_j_lower .and. &
      len_trim(config%fine_output_file) > 0 .and. &
      trim(config%fine_output_file) /= trim(config%eb%flow%output_file) .and. &
      len_trim(config%level_two_output_file) > 0 .and. &
      trim(config%level_two_output_file) /= &
        trim(config%eb%flow%output_file) .and. &
      trim(config%level_two_output_file) /= trim(config%fine_output_file) .and. &
      (.not. config%dynamic_regridding .or. &
       (.not. config%remove_fine_patch_when_untagged .and. &
        level_one_nx >= 8 .and. level_one_ny >= 8 .and. &
        config%regrid_minimum_patch_cells_x <= level_one_nx - 4 .and. &
        config%regrid_minimum_patch_cells_y <= level_one_ny - 4))
  end function supported_three_level_reactive_eb_amr_config

  subroutine compute_reactive_eb_amr_cfl_timestep_2d( &
      species, coarse_state, coarse_temperature, coarse_geometry, &
      fine_state, fine_temperature, fine_geometry, refinement_ratio, &
      cfl, dt, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: coarse_state(:, :, :), coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    real(dp), intent(in) :: fine_state(:, :, :), fine_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: fine_geometry
    integer, intent(in) :: refinement_ratio
    real(dp), intent(in) :: cfl
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok

    real(dp) :: coarse_dt, fine_dt
    logical :: local_ok

    dt = 0.0_dp
    ok = .false.
    if (refinement_ratio < 2) return
    call compute_reactive_eb_cfl_timestep_2d( &
      species, coarse_state, coarse_temperature, coarse_geometry, &
      cfl, coarse_dt, local_ok)
    if (.not. local_ok) return
    call compute_reactive_eb_cfl_timestep_2d( &
      species, fine_state, fine_temperature, fine_geometry, &
      cfl, fine_dt, local_ok)
    if (.not. local_ok) return
    dt = min(coarse_dt, real(refinement_ratio, dp) * fine_dt)
    ok = ieee_is_finite(dt) .and. dt > 0.0_dp
  end subroutine compute_reactive_eb_amr_cfl_timestep_2d

  subroutine compute_three_level_reactive_eb_cfl_timestep_2d( &
      species, root_state, root_temperature, root_geometry, &
      level_one_state, level_one_temperature, level_one_geometry, &
      root_patch, level_two_state, level_two_temperature, &
      level_two_geometry, level_one_patch, cfl, dt, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: root_state(:, :, :), root_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: root_geometry
    real(dp), intent(in) :: level_one_state(:, :, :)
    real(dp), intent(in) :: level_one_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: level_one_geometry
    type(amr_eb_patch_2d), intent(in) :: root_patch
    real(dp), intent(in) :: level_two_state(:, :, :)
    real(dp), intent(in) :: level_two_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: level_two_geometry
    type(amr_eb_patch_2d), intent(in) :: level_one_patch
    real(dp), intent(in) :: cfl
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok

    real(dp) :: root_dt, level_one_dt, level_two_dt
    logical :: local_ok

    dt = 0.0_dp
    ok = .false.
    if (.not. root_patch%is_valid(root_geometry, level_one_geometry) .or. &
        .not. level_one_patch%is_valid( &
          level_one_geometry, level_two_geometry)) return
    call compute_reactive_eb_cfl_timestep_2d( &
      species, root_state, root_temperature, root_geometry, cfl, &
      root_dt, local_ok)
    if (.not. local_ok) return
    call compute_reactive_eb_cfl_timestep_2d( &
      species, level_one_state, level_one_temperature, level_one_geometry, &
      cfl, level_one_dt, local_ok)
    if (.not. local_ok) return
    call compute_reactive_eb_cfl_timestep_2d( &
      species, level_two_state, level_two_temperature, level_two_geometry, &
      cfl, level_two_dt, local_ok)
    if (.not. local_ok) return
    dt = min(root_dt, &
      real(root_patch%refinement_ratio, dp) * level_one_dt, &
      real(root_patch%refinement_ratio * &
        level_one_patch%refinement_ratio, dp) * level_two_dt)
    ok = ieee_is_finite(dt) .and. dt > 0.0_dp
  end subroutine compute_three_level_reactive_eb_cfl_timestep_2d

  subroutine compute_reactive_eb_patch_set_cfl_timestep_2d( &
      species, coarse_state, coarse_temperature, coarse_geometry, &
      patch_set, cfl, dt, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: coarse_state(:, :, :), coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set
    real(dp), intent(in) :: cfl
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok

    real(dp) :: fine_dt
    logical :: local_ok
    integer :: child, nvar

    dt = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. patch_set%is_valid(coarse_geometry, nvar)) return
    call compute_reactive_eb_cfl_timestep_2d( &
      species, coarse_state, coarse_temperature, coarse_geometry, &
      cfl, dt, local_ok)
    if (.not. local_ok) return
    do child = 1, patch_set%patch_count()
      call compute_reactive_eb_cfl_timestep_2d( &
        species, patch_set%children(child)%state, &
        patch_set%children(child)%temperature, &
        patch_set%children(child)%geometry, cfl, fine_dt, local_ok)
      if (.not. local_ok) return
      dt = min(dt, real( &
        patch_set%children(child)%patch%refinement_ratio, dp) * fine_dt)
    end do
    ok = ieee_is_finite(dt) .and. dt > 0.0_dp
  end subroutine compute_reactive_eb_patch_set_cfl_timestep_2d

  subroutine advance_reactive_eb_chemistry_level_2d( &
      species, reactions, geometry, interval, rtol, atol, state, &
      temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), intent(in) :: interval, rtol, atol
    real(dp), intent(inout) :: state(:, :, :), temperature(:, :)
    logical, intent(out) :: ok

    logical, allocatable :: active_mask(:, :)

    ok = .false.
    if (.not. geometry%is_valid()) return
    allocate(active_mask(geometry%nx, geometry%ny))
    active_mask = geometry%cell_type /= eb_covered_cell
    call advance_reactive_chemistry_2d( &
      species, reactions, state, temperature, geometry%nx, geometry%ny, &
      interval, rtol, atol, ok, active_mask)
  end subroutine advance_reactive_eb_chemistry_level_2d

  subroutine advance_two_level_reactive_eb_strang_2d( &
      species, reactions, coarse_state, coarse_temperature, coarse_geometry, &
      fine_state, fine_temperature, fine_geometry, patch, solver, &
      reconstruction, limiter, state_redist_max_order, dt, &
      chemistry_enabled, rtol, atol, new_coarse_state, &
      new_coarse_temperature, new_fine_state, new_fine_temperature, ok, &
      target_volume_fraction, transport, transport_enabled, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, &
      minimum_transport_theta, boundaries)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: coarse_state(:, :, :), coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    real(dp), intent(in) :: fine_state(:, :, :), fine_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    character(len=*), intent(in) :: solver, reconstruction, limiter
    integer, intent(in) :: state_redist_max_order
    real(dp), intent(in) :: dt, rtol, atol
    logical, intent(in) :: chemistry_enabled
    real(dp), intent(out) :: new_coarse_state(:, :, :)
    real(dp), intent(out) :: new_coarse_temperature(:, :)
    real(dp), intent(out) :: new_fine_state(:, :, :)
    real(dp), intent(out) :: new_fine_temperature(:, :)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: target_volume_fraction
    type(gas_transport_species), intent(in), optional :: transport(:)
    logical, intent(in), optional :: transport_enabled, viscosity_enabled
    logical, intent(in), optional :: thermal_conduction_enabled
    logical, intent(in), optional :: species_diffusion_enabled
    logical, intent(in), optional :: barodiffusion_enabled
    real(dp), intent(out), optional :: minimum_transport_theta
    type(reactive_boundary_set_2d), intent(in), optional :: boundaries

    real(dp), allocatable :: candidate_coarse_state(:, :, :)
    real(dp), allocatable :: candidate_coarse_temperature(:, :)
    real(dp), allocatable :: candidate_fine_state(:, :, :)
    real(dp), allocatable :: candidate_fine_temperature(:, :)
    real(dp), allocatable :: hydro_coarse_state(:, :, :)
    real(dp), allocatable :: hydro_coarse_temperature(:, :)
    real(dp), allocatable :: hydro_fine_state(:, :, :)
    real(dp), allocatable :: hydro_fine_temperature(:, :)
    real(dp), allocatable :: transport_coarse_state(:, :, :)
    real(dp), allocatable :: transport_coarse_temperature(:, :)
    real(dp), allocatable :: transport_fine_state(:, :, :)
    real(dp), allocatable :: transport_fine_temperature(:, :)
    real(dp), allocatable :: synchronized_coarse_state(:, :, :)
    real(dp), allocatable :: synchronized_coarse_temperature(:, :)
    logical :: local_ok, use_transport, use_viscosity, use_conduction
    logical :: use_diffusion, use_barodiffusion
    real(dp) :: selected_target, stage_transport_theta, transport_theta

    new_coarse_state = coarse_state
    new_coarse_temperature = coarse_temperature
    new_fine_state = fine_state
    new_fine_temperature = fine_temperature
    ok = .false.
    if (present(minimum_transport_theta)) minimum_transport_theta = 1.0_dp
    if (chemistry_enabled .and. size(reactions) < 1) return
    selected_target = 0.5_dp
    if (present(target_volume_fraction)) &
      selected_target = target_volume_fraction
    use_transport = .false.
    use_viscosity = .true.
    use_conduction = .true.
    use_diffusion = .true.
    use_barodiffusion = .true.
    if (present(transport_enabled)) use_transport = transport_enabled
    if (present(viscosity_enabled)) use_viscosity = viscosity_enabled
    if (present(thermal_conduction_enabled)) &
      use_conduction = thermal_conduction_enabled
    if (present(species_diffusion_enabled)) &
      use_diffusion = species_diffusion_enabled
    if (present(barodiffusion_enabled)) &
      use_barodiffusion = barodiffusion_enabled
    if (use_transport .and. .not. present(transport)) return
    if (use_transport .and. .not. present(boundaries)) return
    transport_theta = 1.0_dp
    allocate(candidate_coarse_state, source=coarse_state)
    allocate(candidate_coarse_temperature, source=coarse_temperature)
    allocate(candidate_fine_state, source=fine_state)
    allocate(candidate_fine_temperature, source=fine_temperature)
    if (chemistry_enabled) then
      call advance_reactive_eb_chemistry_level_2d( &
        species, reactions, coarse_geometry, 0.5_dp * dt, rtol, atol, &
        candidate_coarse_state, candidate_coarse_temperature, local_ok)
      if (.not. local_ok) return
      call advance_reactive_eb_chemistry_level_2d( &
        species, reactions, fine_geometry, 0.5_dp * dt, rtol, atol, &
        candidate_fine_state, candidate_fine_temperature, local_ok)
      if (.not. local_ok) return
    end if

    if (use_transport) then
      allocate(transport_coarse_state, mold=coarse_state)
      allocate(transport_coarse_temperature, mold=coarse_temperature)
      allocate(transport_fine_state, mold=fine_state)
      allocate(transport_fine_temperature, mold=fine_temperature)
      call advance_two_level_reactive_eb_transport_2d( &
        species, transport, candidate_coarse_state, &
        candidate_coarse_temperature, coarse_geometry, candidate_fine_state, &
        candidate_fine_temperature, fine_geometry, patch, 0.5_dp * dt, &
        use_viscosity, use_conduction, use_diffusion, use_barodiffusion, &
        boundaries, selected_target, state_redist_max_order, &
        transport_coarse_state, transport_coarse_temperature, &
        transport_fine_state, transport_fine_temperature, &
        stage_transport_theta, local_ok)
      if (.not. local_ok) return
      candidate_coarse_state = transport_coarse_state
      candidate_coarse_temperature = transport_coarse_temperature
      candidate_fine_state = transport_fine_state
      candidate_fine_temperature = transport_fine_temperature
      transport_theta = min(transport_theta, stage_transport_theta)
    end if

    allocate(hydro_coarse_state, mold=coarse_state)
    allocate(hydro_coarse_temperature, mold=coarse_temperature)
    allocate(hydro_fine_state, mold=fine_state)
    allocate(hydro_fine_temperature, mold=fine_temperature)
    call advance_two_level_reactive_eb_hydro_2d( &
      species, candidate_coarse_state, candidate_coarse_temperature, &
      coarse_geometry, candidate_fine_state, candidate_fine_temperature, &
      fine_geometry, patch, solver, reconstruction, limiter, &
      state_redist_max_order, dt, hydro_coarse_state, &
      hydro_coarse_temperature, hydro_fine_state, hydro_fine_temperature, &
      local_ok, target_volume_fraction)
    if (.not. local_ok) return
    candidate_coarse_state = hydro_coarse_state
    candidate_coarse_temperature = hydro_coarse_temperature
    candidate_fine_state = hydro_fine_state
    candidate_fine_temperature = hydro_fine_temperature

    if (use_transport) then
      call advance_two_level_reactive_eb_transport_2d( &
        species, transport, candidate_coarse_state, &
        candidate_coarse_temperature, coarse_geometry, candidate_fine_state, &
        candidate_fine_temperature, fine_geometry, patch, 0.5_dp * dt, &
        use_viscosity, use_conduction, use_diffusion, use_barodiffusion, &
        boundaries, selected_target, state_redist_max_order, &
        transport_coarse_state, transport_coarse_temperature, &
        transport_fine_state, transport_fine_temperature, &
        stage_transport_theta, local_ok)
      if (.not. local_ok) return
      candidate_coarse_state = transport_coarse_state
      candidate_coarse_temperature = transport_coarse_temperature
      candidate_fine_state = transport_fine_state
      candidate_fine_temperature = transport_fine_temperature
      transport_theta = min(transport_theta, stage_transport_theta)
    end if

    if (chemistry_enabled) then
      call advance_reactive_eb_chemistry_level_2d( &
        species, reactions, coarse_geometry, 0.5_dp * dt, rtol, atol, &
        candidate_coarse_state, candidate_coarse_temperature, local_ok)
      if (.not. local_ok) return
      call advance_reactive_eb_chemistry_level_2d( &
        species, reactions, fine_geometry, 0.5_dp * dt, rtol, atol, &
        candidate_fine_state, candidate_fine_temperature, local_ok)
      if (.not. local_ok) return
      allocate(synchronized_coarse_state, mold=coarse_state)
      allocate(synchronized_coarse_temperature, mold=coarse_temperature)
      call average_down_reactive_eb_state_patch_2d( &
        species, candidate_coarse_state, candidate_coarse_temperature, &
        coarse_geometry, candidate_fine_state, fine_geometry, patch, &
        synchronized_coarse_state, synchronized_coarse_temperature, local_ok)
      if (.not. local_ok) return
      candidate_coarse_state = synchronized_coarse_state
      candidate_coarse_temperature = synchronized_coarse_temperature
    end if
    new_coarse_state = candidate_coarse_state
    new_coarse_temperature = candidate_coarse_temperature
    new_fine_state = candidate_fine_state
    new_fine_temperature = candidate_fine_temperature
    if (present(minimum_transport_theta)) &
      minimum_transport_theta = transport_theta
    ok = .true.
  end subroutine advance_two_level_reactive_eb_strang_2d

  subroutine advance_three_level_reactive_eb_strang_2d( &
      species, reactions, root_state, root_temperature, root_geometry, &
      level_one_state, level_one_temperature, level_one_geometry, &
      root_patch, level_two_state, level_two_temperature, &
      level_two_geometry, level_one_patch, solver, reconstruction, limiter, &
      state_redist_max_order, dt, chemistry_enabled, rtol, atol, &
      new_root_state, new_root_temperature, new_level_one_state, &
      new_level_one_temperature, new_level_two_state, &
      new_level_two_temperature, ok, target_volume_fraction, transport, &
      transport_enabled, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, &
      minimum_transport_theta, boundaries)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: root_state(:, :, :), root_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: root_geometry
    real(dp), intent(in) :: level_one_state(:, :, :)
    real(dp), intent(in) :: level_one_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: level_one_geometry
    type(amr_eb_patch_2d), intent(in) :: root_patch
    real(dp), intent(in) :: level_two_state(:, :, :)
    real(dp), intent(in) :: level_two_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: level_two_geometry
    type(amr_eb_patch_2d), intent(in) :: level_one_patch
    character(len=*), intent(in) :: solver, reconstruction, limiter
    integer, intent(in) :: state_redist_max_order
    real(dp), intent(in) :: dt, rtol, atol
    logical, intent(in) :: chemistry_enabled
    real(dp), intent(out) :: new_root_state(:, :, :)
    real(dp), intent(out) :: new_root_temperature(:, :)
    real(dp), intent(out) :: new_level_one_state(:, :, :)
    real(dp), intent(out) :: new_level_one_temperature(:, :)
    real(dp), intent(out) :: new_level_two_state(:, :, :)
    real(dp), intent(out) :: new_level_two_temperature(:, :)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: target_volume_fraction
    type(gas_transport_species), intent(in), optional :: transport(:)
    logical, intent(in), optional :: transport_enabled, viscosity_enabled
    logical, intent(in), optional :: thermal_conduction_enabled
    logical, intent(in), optional :: species_diffusion_enabled
    logical, intent(in), optional :: barodiffusion_enabled
    real(dp), intent(out), optional :: minimum_transport_theta
    type(reactive_boundary_set_2d), intent(in), optional :: boundaries

    real(dp), allocatable :: root_candidate(:, :, :)
    real(dp), allocatable :: root_candidate_temperature(:, :)
    real(dp), allocatable :: level_one_candidate(:, :, :)
    real(dp), allocatable :: level_one_candidate_temperature(:, :)
    real(dp), allocatable :: level_two_candidate(:, :, :)
    real(dp), allocatable :: level_two_candidate_temperature(:, :)
    real(dp), allocatable :: hydro_root(:, :, :)
    real(dp), allocatable :: hydro_root_temperature(:, :)
    real(dp), allocatable :: hydro_level_one(:, :, :)
    real(dp), allocatable :: hydro_level_one_temperature(:, :)
    real(dp), allocatable :: hydro_level_two(:, :, :)
    real(dp), allocatable :: hydro_level_two_temperature(:, :)
    real(dp), allocatable :: transport_root(:, :, :)
    real(dp), allocatable :: transport_root_temperature(:, :)
    real(dp), allocatable :: transport_level_one(:, :, :)
    real(dp), allocatable :: transport_level_one_temperature(:, :)
    real(dp), allocatable :: transport_level_two(:, :, :)
    real(dp), allocatable :: transport_level_two_temperature(:, :)
    real(dp), allocatable :: synchronized_root(:, :, :)
    real(dp), allocatable :: synchronized_root_temperature(:, :)
    real(dp), allocatable :: synchronized_level_one(:, :, :)
    real(dp), allocatable :: synchronized_level_one_temperature(:, :)
    logical :: local_ok, use_transport, use_viscosity, use_conduction
    logical :: use_diffusion, use_barodiffusion
    real(dp) :: selected_target, stage_transport_theta, transport_theta

    new_root_state = root_state
    new_root_temperature = root_temperature
    new_level_one_state = level_one_state
    new_level_one_temperature = level_one_temperature
    new_level_two_state = level_two_state
    new_level_two_temperature = level_two_temperature
    ok = .false.
    if (present(minimum_transport_theta)) minimum_transport_theta = 1.0_dp
    if (chemistry_enabled .and. size(reactions) < 1) return
    selected_target = 0.5_dp
    if (present(target_volume_fraction)) &
      selected_target = target_volume_fraction
    use_transport = .false.
    use_viscosity = .true.
    use_conduction = .true.
    use_diffusion = .true.
    use_barodiffusion = .true.
    if (present(transport_enabled)) use_transport = transport_enabled
    if (present(viscosity_enabled)) use_viscosity = viscosity_enabled
    if (present(thermal_conduction_enabled)) &
      use_conduction = thermal_conduction_enabled
    if (present(species_diffusion_enabled)) &
      use_diffusion = species_diffusion_enabled
    if (present(barodiffusion_enabled)) &
      use_barodiffusion = barodiffusion_enabled
    if (use_transport .and. .not. present(transport)) return
    if (use_transport .and. .not. present(boundaries)) return
    transport_theta = 1.0_dp
    allocate(root_candidate, source=root_state)
    allocate(root_candidate_temperature, source=root_temperature)
    allocate(level_one_candidate, source=level_one_state)
    allocate(level_one_candidate_temperature, source=level_one_temperature)
    allocate(level_two_candidate, source=level_two_state)
    allocate(level_two_candidate_temperature, source=level_two_temperature)
    if (chemistry_enabled) then
      call advance_reactive_eb_chemistry_level_2d( &
        species, reactions, root_geometry, 0.5_dp * dt, rtol, atol, &
        root_candidate, root_candidate_temperature, local_ok)
      if (.not. local_ok) return
      call advance_reactive_eb_chemistry_level_2d( &
        species, reactions, level_one_geometry, 0.5_dp * dt, rtol, atol, &
        level_one_candidate, level_one_candidate_temperature, local_ok)
      if (.not. local_ok) return
      call advance_reactive_eb_chemistry_level_2d( &
        species, reactions, level_two_geometry, 0.5_dp * dt, rtol, atol, &
        level_two_candidate, level_two_candidate_temperature, local_ok)
      if (.not. local_ok) return
    end if

    if (use_transport) then
      allocate(transport_root, mold=root_state)
      allocate(transport_root_temperature, mold=root_temperature)
      allocate(transport_level_one, mold=level_one_state)
      allocate(transport_level_one_temperature, mold=level_one_temperature)
      allocate(transport_level_two, mold=level_two_state)
      allocate(transport_level_two_temperature, mold=level_two_temperature)
      call advance_three_level_reactive_eb_transport_2d( &
        species, transport, root_candidate, root_candidate_temperature, &
        root_geometry, level_one_candidate, level_one_candidate_temperature, &
        level_one_geometry, root_patch, level_two_candidate, &
        level_two_candidate_temperature, level_two_geometry, level_one_patch, &
        0.5_dp * dt, use_viscosity, use_conduction, use_diffusion, &
        use_barodiffusion, boundaries, selected_target, &
        state_redist_max_order, transport_root, transport_root_temperature, &
        transport_level_one, transport_level_one_temperature, &
        transport_level_two, transport_level_two_temperature, &
        stage_transport_theta, local_ok)
      if (.not. local_ok) return
      root_candidate = transport_root
      root_candidate_temperature = transport_root_temperature
      level_one_candidate = transport_level_one
      level_one_candidate_temperature = transport_level_one_temperature
      level_two_candidate = transport_level_two
      level_two_candidate_temperature = transport_level_two_temperature
      transport_theta = min(transport_theta, stage_transport_theta)
    end if

    allocate(hydro_root, mold=root_state)
    allocate(hydro_root_temperature, mold=root_temperature)
    allocate(hydro_level_one, mold=level_one_state)
    allocate(hydro_level_one_temperature, mold=level_one_temperature)
    allocate(hydro_level_two, mold=level_two_state)
    allocate(hydro_level_two_temperature, mold=level_two_temperature)
    call advance_three_level_reactive_eb_hydro_2d( &
      species, root_candidate, root_candidate_temperature, root_geometry, &
      level_one_candidate, level_one_candidate_temperature, &
      level_one_geometry, root_patch, level_two_candidate, &
      level_two_candidate_temperature, level_two_geometry, level_one_patch, &
      solver, reconstruction, limiter, state_redist_max_order, dt, &
      hydro_root, hydro_root_temperature, hydro_level_one, &
      hydro_level_one_temperature, hydro_level_two, &
      hydro_level_two_temperature, local_ok, target_volume_fraction)
    if (.not. local_ok) return
    root_candidate = hydro_root
    root_candidate_temperature = hydro_root_temperature
    level_one_candidate = hydro_level_one
    level_one_candidate_temperature = hydro_level_one_temperature
    level_two_candidate = hydro_level_two
    level_two_candidate_temperature = hydro_level_two_temperature

    if (use_transport) then
      call advance_three_level_reactive_eb_transport_2d( &
        species, transport, root_candidate, root_candidate_temperature, &
        root_geometry, level_one_candidate, level_one_candidate_temperature, &
        level_one_geometry, root_patch, level_two_candidate, &
        level_two_candidate_temperature, level_two_geometry, level_one_patch, &
        0.5_dp * dt, use_viscosity, use_conduction, use_diffusion, &
        use_barodiffusion, boundaries, selected_target, &
        state_redist_max_order, transport_root, transport_root_temperature, &
        transport_level_one, transport_level_one_temperature, &
        transport_level_two, transport_level_two_temperature, &
        stage_transport_theta, local_ok)
      if (.not. local_ok) return
      root_candidate = transport_root
      root_candidate_temperature = transport_root_temperature
      level_one_candidate = transport_level_one
      level_one_candidate_temperature = transport_level_one_temperature
      level_two_candidate = transport_level_two
      level_two_candidate_temperature = transport_level_two_temperature
      transport_theta = min(transport_theta, stage_transport_theta)
    end if

    if (chemistry_enabled) then
      call advance_reactive_eb_chemistry_level_2d( &
        species, reactions, root_geometry, 0.5_dp * dt, rtol, atol, &
        root_candidate, root_candidate_temperature, local_ok)
      if (.not. local_ok) return
      call advance_reactive_eb_chemistry_level_2d( &
        species, reactions, level_one_geometry, 0.5_dp * dt, rtol, atol, &
        level_one_candidate, level_one_candidate_temperature, local_ok)
      if (.not. local_ok) return
      call advance_reactive_eb_chemistry_level_2d( &
        species, reactions, level_two_geometry, 0.5_dp * dt, rtol, atol, &
        level_two_candidate, level_two_candidate_temperature, local_ok)
      if (.not. local_ok) return
      allocate(synchronized_root, mold=root_state)
      allocate(synchronized_root_temperature, mold=root_temperature)
      allocate(synchronized_level_one, mold=level_one_state)
      allocate(synchronized_level_one_temperature, mold=level_one_temperature)
      call average_down_three_level_reactive_eb_state_2d( &
        species, root_candidate, root_candidate_temperature, root_geometry, &
        level_one_candidate, level_one_candidate_temperature, &
        level_one_geometry, root_patch, level_two_candidate, &
        level_two_candidate_temperature, level_two_geometry, level_one_patch, &
        synchronized_root, synchronized_root_temperature, &
        synchronized_level_one, synchronized_level_one_temperature, local_ok)
      if (.not. local_ok) return
      root_candidate = synchronized_root
      root_candidate_temperature = synchronized_root_temperature
      level_one_candidate = synchronized_level_one
      level_one_candidate_temperature = synchronized_level_one_temperature
    end if

    new_root_state = root_candidate
    new_root_temperature = root_candidate_temperature
    new_level_one_state = level_one_candidate
    new_level_one_temperature = level_one_candidate_temperature
    new_level_two_state = level_two_candidate
    new_level_two_temperature = level_two_candidate_temperature
    if (present(minimum_transport_theta)) &
      minimum_transport_theta = transport_theta
    ok = .true.
  end subroutine advance_three_level_reactive_eb_strang_2d

  subroutine advance_reactive_eb_patch_set_strang_2d( &
      species, reactions, coarse_state, coarse_temperature, coarse_geometry, &
      patch_set, solver, reconstruction, limiter, state_redist_max_order, &
      dt, chemistry_enabled, rtol, atol, new_coarse_state, &
      new_coarse_temperature, new_patch_set, ok, target_volume_fraction, &
      failure_context, transport, transport_enabled, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, minimum_transport_theta, boundaries)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: coarse_state(:, :, :), coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set
    character(len=*), intent(in) :: solver, reconstruction, limiter
    integer, intent(in) :: state_redist_max_order
    real(dp), intent(in) :: dt, rtol, atol
    logical, intent(in) :: chemistry_enabled
    real(dp), intent(out) :: new_coarse_state(:, :, :)
    real(dp), intent(out) :: new_coarse_temperature(:, :)
    type(reactive_eb_patch_set_2d), intent(out) :: new_patch_set
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: target_volume_fraction
    character(len=*), intent(out), optional :: failure_context
    type(gas_transport_species), intent(in), optional :: transport(:)
    logical, intent(in), optional :: transport_enabled, viscosity_enabled
    logical, intent(in), optional :: thermal_conduction_enabled
    logical, intent(in), optional :: species_diffusion_enabled
    logical, intent(in), optional :: barodiffusion_enabled
    real(dp), intent(out), optional :: minimum_transport_theta
    type(reactive_boundary_set_2d), intent(in), optional :: boundaries

    type(reactive_eb_patch_set_2d) :: candidate_set, hydro_set, transport_set
    real(dp), allocatable :: candidate_coarse(:, :, :)
    real(dp), allocatable :: candidate_coarse_temperature(:, :)
    real(dp), allocatable :: hydro_coarse(:, :, :)
    real(dp), allocatable :: hydro_coarse_temperature(:, :)
    real(dp), allocatable :: transport_coarse(:, :, :)
    real(dp), allocatable :: transport_coarse_temperature(:, :)
    real(dp), allocatable :: synchronized_coarse(:, :, :)
    real(dp), allocatable :: synchronized_coarse_temperature(:, :)
    logical :: local_ok, use_transport, use_viscosity, use_conduction
    logical :: use_diffusion, use_barodiffusion
    real(dp) :: selected_target, stage_transport_theta, transport_theta
    integer :: child, nvar

    new_coarse_state = 0.0_dp
    new_coarse_temperature = 0.0_dp
    new_patch_set = reactive_eb_patch_set_2d()
    ok = .false.
    if (present(minimum_transport_theta)) minimum_transport_theta = 1.0_dp
    if (present(failure_context)) failure_context = "input validation"
    nvar = reactive_nvar(size(species))
    if (nvar < 1 .or. &
        any(shape(new_coarse_state) /= shape(coarse_state)) .or. &
        any(shape(new_coarse_temperature) /= shape(coarse_temperature))) &
      return
    new_coarse_state = coarse_state
    new_coarse_temperature = coarse_temperature
    new_patch_set = patch_set
    if (.not. patch_set%is_valid(coarse_geometry, nvar) .or. &
        (chemistry_enabled .and. size(reactions) < 1)) return
    selected_target = 0.5_dp
    if (present(target_volume_fraction)) &
      selected_target = target_volume_fraction
    use_transport = .false.
    use_viscosity = .true.
    use_conduction = .true.
    use_diffusion = .true.
    use_barodiffusion = .true.
    if (present(transport_enabled)) use_transport = transport_enabled
    if (present(viscosity_enabled)) use_viscosity = viscosity_enabled
    if (present(thermal_conduction_enabled)) &
      use_conduction = thermal_conduction_enabled
    if (present(species_diffusion_enabled)) &
      use_diffusion = species_diffusion_enabled
    if (present(barodiffusion_enabled)) &
      use_barodiffusion = barodiffusion_enabled
    if (use_transport .and. .not. present(transport)) return
    if (use_transport .and. .not. present(boundaries)) return
    transport_theta = 1.0_dp

    allocate(candidate_coarse, source=coarse_state)
    allocate(candidate_coarse_temperature, source=coarse_temperature)
    candidate_set = patch_set
    if (chemistry_enabled) then
      if (present(failure_context)) &
        failure_context = "first coarse chemistry half-step"
      call advance_reactive_eb_chemistry_level_2d( &
        species, reactions, coarse_geometry, 0.5_dp * dt, rtol, atol, &
        candidate_coarse, candidate_coarse_temperature, local_ok)
      if (.not. local_ok) return
      do child = 1, candidate_set%patch_count()
        if (present(failure_context)) &
          failure_context = "first fine chemistry half-step"
        call advance_reactive_eb_chemistry_level_2d( &
          species, reactions, candidate_set%children(child)%geometry, &
          0.5_dp * dt, rtol, atol, candidate_set%children(child)%state, &
          candidate_set%children(child)%temperature, local_ok)
        if (.not. local_ok) return
      end do
    end if

    if (use_transport) then
      allocate(transport_coarse, mold=coarse_state)
      allocate(transport_coarse_temperature, mold=coarse_temperature)
      if (present(failure_context)) &
        failure_context = "first multipatch transport half-step"
      call advance_reactive_eb_patch_set_transport_2d( &
        species, transport, candidate_coarse, candidate_coarse_temperature, &
        coarse_geometry, candidate_set, 0.5_dp * dt, use_viscosity, &
        use_conduction, use_diffusion, use_barodiffusion, boundaries, &
        selected_target, state_redist_max_order, transport_coarse, &
        transport_coarse_temperature, transport_set, stage_transport_theta, &
        local_ok, failure_context)
      if (.not. local_ok) return
      candidate_coarse = transport_coarse
      candidate_coarse_temperature = transport_coarse_temperature
      candidate_set = transport_set
      transport_theta = min(transport_theta, stage_transport_theta)
    end if

    allocate(hydro_coarse, mold=coarse_state)
    allocate(hydro_coarse_temperature, mold=coarse_temperature)
    if (present(failure_context)) failure_context = "multipatch hydro"
    call advance_reactive_eb_patch_set_hydro_2d( &
      species, candidate_coarse, candidate_coarse_temperature, &
      coarse_geometry, candidate_set, solver, reconstruction, limiter, &
      state_redist_max_order, dt, hydro_coarse, hydro_coarse_temperature, &
      hydro_set, local_ok, target_volume_fraction, failure_context)
    if (.not. local_ok) return
    candidate_coarse = hydro_coarse
    candidate_coarse_temperature = hydro_coarse_temperature
    candidate_set = hydro_set

    if (use_transport) then
      if (present(failure_context)) &
        failure_context = "second multipatch transport half-step"
      call advance_reactive_eb_patch_set_transport_2d( &
        species, transport, candidate_coarse, candidate_coarse_temperature, &
        coarse_geometry, candidate_set, 0.5_dp * dt, use_viscosity, &
        use_conduction, use_diffusion, use_barodiffusion, boundaries, &
        selected_target, state_redist_max_order, transport_coarse, &
        transport_coarse_temperature, transport_set, stage_transport_theta, &
        local_ok, failure_context)
      if (.not. local_ok) return
      candidate_coarse = transport_coarse
      candidate_coarse_temperature = transport_coarse_temperature
      candidate_set = transport_set
      transport_theta = min(transport_theta, stage_transport_theta)
    end if

    if (chemistry_enabled) then
      if (present(failure_context)) &
        failure_context = "second coarse chemistry half-step"
      call advance_reactive_eb_chemistry_level_2d( &
        species, reactions, coarse_geometry, 0.5_dp * dt, rtol, atol, &
        candidate_coarse, candidate_coarse_temperature, local_ok)
      if (.not. local_ok) return
      do child = 1, candidate_set%patch_count()
        if (present(failure_context)) &
          failure_context = "second fine chemistry half-step"
        call advance_reactive_eb_chemistry_level_2d( &
          species, reactions, candidate_set%children(child)%geometry, &
          0.5_dp * dt, rtol, atol, candidate_set%children(child)%state, &
          candidate_set%children(child)%temperature, local_ok)
        if (.not. local_ok) return
      end do
      allocate(synchronized_coarse, mold=coarse_state)
      allocate(synchronized_coarse_temperature, mold=coarse_temperature)
      if (present(failure_context)) &
        failure_context = "post-chemistry patch-set synchronization"
      call average_down_reactive_eb_patch_set_2d( &
        species, candidate_coarse, candidate_coarse_temperature, &
        coarse_geometry, candidate_set, synchronized_coarse, &
        synchronized_coarse_temperature, local_ok)
      if (.not. local_ok) return
      candidate_coarse = synchronized_coarse
      candidate_coarse_temperature = synchronized_coarse_temperature
    end if

    new_coarse_state = candidate_coarse
    new_coarse_temperature = candidate_coarse_temperature
    new_patch_set = candidate_set
    if (present(minimum_transport_theta)) &
      minimum_transport_theta = transport_theta
    ok = .true.
    if (present(failure_context)) failure_context = "none"
  end subroutine advance_reactive_eb_patch_set_strang_2d

  subroutine compute_reactive_eb_amr_integrals_2d( &
      coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
      fine_active, integrals, ok)
    real(dp), intent(in) :: coarse_state(:, :, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    real(dp), allocatable, intent(in) :: fine_state(:, :, :)
    type(eb_geometry_2d), intent(in) :: fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    logical, intent(in) :: fine_active
    real(dp), intent(out) :: integrals(:)
    logical, intent(out) :: ok

    if (fine_active) then
      ok = allocated(fine_state)
      if (.not. ok) return
      call composite_eb_integral_2d( &
        coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
        integrals, ok)
      return
    end if
    if (allocated(fine_state)) then
      integrals = 0.0_dp
      ok = .false.
      return
    end if
    call reactive_eb_integrals_2d( &
      coarse_state, coarse_geometry, integrals, ok)
  end subroutine compute_reactive_eb_amr_integrals_2d

  subroutine build_reactive_eb_amr_patch_geometry_2d( &
      config, coarse_geometry, coarse_i_lower, coarse_i_upper, &
      coarse_j_lower, coarse_j_upper, fine_geometry, patch, ok)
    type(reactive_eb_amr_2d_config), intent(in) :: config
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    integer, intent(in) :: coarse_i_lower, coarse_i_upper
    integer, intent(in) :: coarse_j_lower, coarse_j_upper
    type(eb_geometry_2d), intent(out) :: fine_geometry
    type(amr_eb_patch_2d), intent(out) :: patch
    logical, intent(out) :: ok

    real(dp) :: fine_x_lower, fine_x_upper, fine_y_lower, fine_y_upper
    integer :: fine_nx, fine_ny

    ok = .false.
    if (coarse_i_lower < 1 .or. &
        coarse_i_upper > coarse_geometry%nx .or. &
        coarse_j_lower < 1 .or. &
        coarse_j_upper > coarse_geometry%ny .or. &
        coarse_i_upper < coarse_i_lower .or. &
        coarse_j_upper < coarse_j_lower) return
    fine_x_lower = coarse_geometry%x_lower + &
      real(coarse_i_lower - 1, dp) * coarse_geometry%dx
    fine_x_upper = coarse_geometry%x_lower + &
      real(coarse_i_upper, dp) * coarse_geometry%dx
    fine_y_lower = coarse_geometry%y_lower + &
      real(coarse_j_lower - 1, dp) * coarse_geometry%dy
    fine_y_upper = coarse_geometry%y_lower + &
      real(coarse_j_upper, dp) * coarse_geometry%dy
    fine_nx = (coarse_i_upper - coarse_i_lower + 1) * &
      config%refinement_ratio
    fine_ny = (coarse_j_upper - coarse_j_lower + 1) * &
      config%refinement_ratio
    call build_configured_eb_geometry_region_2d( &
      config%eb, fine_nx, fine_ny, fine_x_lower, fine_x_upper, &
      fine_y_lower, fine_y_upper, fine_geometry, ok)
    if (.not. ok) return
    call build_amr_eb_patch_2d( &
      coarse_geometry, fine_geometry, coarse_i_lower, coarse_i_upper, &
      coarse_j_lower, coarse_j_upper, config%refinement_ratio, patch, ok)
  end subroutine build_reactive_eb_amr_patch_geometry_2d

  pure subroutine build_initial_reactive_eb_patch_collection_2d( &
      config, collection)
    type(reactive_eb_amr_2d_config), intent(in) :: config
    type(amr_eb_regrid_plan_collection_2d), intent(out) :: collection

    integer :: tagged_count

    tagged_count = (config%coarse_i_upper - config%coarse_i_lower + 1) * &
      (config%coarse_j_upper - config%coarse_j_lower + 1)
    collection = amr_eb_regrid_plan_collection_2d( &
      coarse_nx=config%eb%flow%nx, coarse_ny=config%eb%flow%ny, &
      tagged_cell_count=tagged_count)
    allocate(collection%plans(1))
    collection%plans(1) = amr_eb_regrid_plan_2d( &
      active=.true., coarse_nx=config%eb%flow%nx, &
      coarse_ny=config%eb%flow%ny, tagged_cell_count=tagged_count, &
      tag_i_lower=config%coarse_i_lower, &
      tag_i_upper=config%coarse_i_upper, &
      tag_j_lower=config%coarse_j_lower, &
      tag_j_upper=config%coarse_j_upper, &
      coarse_i_lower=config%coarse_i_lower, &
      coarse_i_upper=config%coarse_i_upper, &
      coarse_j_lower=config%coarse_j_lower, &
      coarse_j_upper=config%coarse_j_upper)
  end subroutine build_initial_reactive_eb_patch_collection_2d

  subroutine build_reactive_eb_patch_set_geometries_2d( &
      config, coarse_geometry, collection, fine_geometries, ok)
    type(reactive_eb_amr_2d_config), intent(in) :: config
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(amr_eb_regrid_plan_collection_2d), intent(in) :: collection
    type(eb_geometry_2d), allocatable, intent(out) :: fine_geometries(:)
    logical, intent(out) :: ok

    type(amr_eb_patch_2d) :: patch
    logical :: local_ok
    integer :: child

    ok = .false.
    if (.not. collection%is_valid() .or. &
        collection%coarse_nx /= coarse_geometry%nx .or. &
        collection%coarse_ny /= coarse_geometry%ny) return
    allocate(fine_geometries(collection%patch_count()))
    do child = 1, collection%patch_count()
      call build_reactive_eb_amr_patch_geometry_2d( &
        config, coarse_geometry, &
        collection%plans(child)%coarse_i_lower, &
        collection%plans(child)%coarse_i_upper, &
        collection%plans(child)%coarse_j_lower, &
        collection%plans(child)%coarse_j_upper, &
        fine_geometries(child), patch, local_ok)
      if (.not. local_ok) return
    end do
    ok = .true.
  end subroutine build_reactive_eb_patch_set_geometries_2d

  pure logical function reactive_eb_patch_set_matches_collection_2d( &
      patch_set, collection, coarse_geometry, nvar) result(matches)
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set
    type(amr_eb_regrid_plan_collection_2d), intent(in) :: collection
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    integer, intent(in) :: nvar

    integer :: child

    matches = patch_set%is_valid(coarse_geometry, nvar) .and. &
      collection%is_valid() .and. &
      patch_set%patch_count() == collection%patch_count()
    if (.not. matches) return
    do child = 1, patch_set%patch_count()
      matches = &
        patch_set%children(child)%patch%coarse_i_lower == &
          collection%plans(child)%coarse_i_lower .and. &
        patch_set%children(child)%patch%coarse_i_upper == &
          collection%plans(child)%coarse_i_upper .and. &
        patch_set%children(child)%patch%coarse_j_lower == &
          collection%plans(child)%coarse_j_lower .and. &
        patch_set%children(child)%patch%coarse_j_upper == &
          collection%plans(child)%coarse_j_upper
      if (.not. matches) return
    end do
  end function reactive_eb_patch_set_matches_collection_2d

  subroutine regrid_reactive_eb_amr_patch_set_2d( &
      species, config, coarse_state, coarse_temperature, coarse_geometry, &
      patch_set, changed, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_eb_amr_2d_config), intent(in) :: config
    real(dp), allocatable, intent(inout) :: coarse_state(:, :, :)
    real(dp), allocatable, intent(inout) :: coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(inout) :: patch_set
    logical, intent(out) :: changed, ok

    type(amr_eb_tagging_criteria_2d) :: criteria
    type(amr_eb_regrid_plan_collection_2d) :: collection
    type(reactive_eb_patch_set_2d) :: candidate_set
    type(eb_geometry_2d), allocatable :: fine_geometries(:)
    real(dp), allocatable :: candidate_coarse(:, :, :)
    real(dp), allocatable :: candidate_coarse_temperature(:, :)
    logical, allocatable :: tags(:, :)
    logical :: local_ok
    integer :: nvar

    changed = .false.
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. config%multipatch_enabled .or. &
        .not. patch_set%is_valid(coarse_geometry, nvar)) return
    criteria%relative_gradient_threshold = &
      config%regrid_relative_temperature_gradient
    criteria%absolute_gradient_threshold = &
      config%regrid_absolute_temperature_gradient
    criteria%scale_floor = config%regrid_temperature_scale_floor
    criteria%buffer_cells = config%regrid_buffer_cells
    criteria%minimum_patch_cells_x = config%regrid_minimum_patch_cells_x
    criteria%minimum_patch_cells_y = config%regrid_minimum_patch_cells_y
    criteria%maximum_patch_gap_cells = &
      config%regrid_maximum_patch_gap_cells
    allocate(tags(coarse_geometry%nx, coarse_geometry%ny))
    call plan_reactive_eb_temperature_regrid_collection_2d( &
      coarse_temperature, coarse_geometry, criteria, tags, collection, &
      local_ok)
    if (.not. local_ok) return
    if (collection%patch_count() == 0 .and. &
        .not. config%remove_fine_patch_when_untagged) then
      ok = .true.
      return
    end if
    if (reactive_eb_patch_set_matches_collection_2d( &
        patch_set, collection, coarse_geometry, nvar)) then
      ok = .true.
      return
    end if
    call build_reactive_eb_patch_set_geometries_2d( &
      config, coarse_geometry, collection, fine_geometries, local_ok)
    if (.not. local_ok) return
    allocate(candidate_coarse, mold=coarse_state)
    allocate(candidate_coarse_temperature, mold=coarse_temperature)
    call regrid_reactive_eb_patch_set_2d( &
      species, coarse_state, coarse_temperature, coarse_geometry, patch_set, &
      fine_geometries, collection, config%refinement_ratio, &
      candidate_coarse, candidate_coarse_temperature, candidate_set, &
      local_ok)
    if (.not. local_ok) return
    call move_alloc(candidate_coarse, coarse_state)
    call move_alloc(candidate_coarse_temperature, coarse_temperature)
    patch_set = candidate_set
    changed = .true.
    ok = .true.
  end subroutine regrid_reactive_eb_amr_patch_set_2d

  subroutine regrid_reactive_eb_amr_hierarchy_2d( &
      species, config, coarse_state, coarse_temperature, coarse_geometry, &
      fine_state, fine_temperature, fine_geometry, patch, fine_active, &
      changed, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_eb_amr_2d_config), intent(in) :: config
    real(dp), allocatable, intent(inout) :: coarse_state(:, :, :)
    real(dp), allocatable, intent(inout) :: coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    real(dp), allocatable, intent(inout) :: fine_state(:, :, :)
    real(dp), allocatable, intent(inout) :: fine_temperature(:, :)
    type(eb_geometry_2d), intent(inout) :: fine_geometry
    type(amr_eb_patch_2d), intent(inout) :: patch
    logical, intent(inout) :: fine_active
    logical, intent(out) :: changed, ok

    type(amr_eb_tagging_criteria_2d) :: criteria
    type(amr_eb_regrid_plan_2d) :: plan
    type(eb_geometry_2d) :: candidate_fine_geometry
    type(amr_eb_patch_2d) :: candidate_patch
    real(dp), allocatable :: candidate_coarse_state(:, :, :)
    real(dp), allocatable :: candidate_coarse_temperature(:, :)
    real(dp), allocatable :: candidate_fine_state(:, :, :)
    real(dp), allocatable :: candidate_fine_temperature(:, :)
    logical, allocatable :: tags(:, :)
    logical :: local_ok
    integer :: nvar

    changed = .false.
    ok = .false.
    if (fine_active) then
      if (.not. allocated(fine_state) .or. &
          .not. allocated(fine_temperature) .or. &
          .not. patch%is_valid(coarse_geometry, fine_geometry)) return
    else
      if (allocated(fine_state) .or. allocated(fine_temperature) .or. &
          fine_geometry%is_valid() .or. patch%refinement_ratio /= 0) return
    end if
    criteria%relative_gradient_threshold = &
      config%regrid_relative_temperature_gradient
    criteria%absolute_gradient_threshold = &
      config%regrid_absolute_temperature_gradient
    criteria%scale_floor = config%regrid_temperature_scale_floor
    criteria%buffer_cells = config%regrid_buffer_cells
    criteria%minimum_patch_cells_x = config%regrid_minimum_patch_cells_x
    criteria%minimum_patch_cells_y = config%regrid_minimum_patch_cells_y
    allocate(tags(coarse_geometry%nx, coarse_geometry%ny))
    call plan_reactive_eb_temperature_regrid_2d( &
      coarse_temperature, coarse_geometry, criteria, tags, plan, local_ok)
    if (.not. local_ok) return
    if (.not. plan%active) then
      if (fine_active .and. config%remove_fine_patch_when_untagged) then
        allocate(candidate_coarse_state, mold=coarse_state)
        allocate(candidate_coarse_temperature, mold=coarse_temperature)
        call collapse_two_level_reactive_eb_patch_2d( &
          species, coarse_state, coarse_temperature, coarse_geometry, &
          fine_state, fine_geometry, patch, candidate_coarse_state, &
          candidate_coarse_temperature, local_ok)
        if (.not. local_ok) return
        call move_alloc(candidate_coarse_state, coarse_state)
        call move_alloc(candidate_coarse_temperature, coarse_temperature)
        deallocate(fine_state, fine_temperature)
        fine_geometry = eb_geometry_2d()
        patch = amr_eb_patch_2d()
        fine_active = .false.
        changed = .true.
      end if
      ok = .true.
      return
    end if
    if (fine_active .and. &
        plan%coarse_i_lower == patch%coarse_i_lower .and. &
        plan%coarse_i_upper == patch%coarse_i_upper .and. &
        plan%coarse_j_lower == patch%coarse_j_lower .and. &
        plan%coarse_j_upper == patch%coarse_j_upper) then
      ok = .true.
      return
    end if

    call build_reactive_eb_amr_patch_geometry_2d( &
      config, coarse_geometry, plan%coarse_i_lower, &
      plan%coarse_i_upper, plan%coarse_j_lower, plan%coarse_j_upper, &
      candidate_fine_geometry, candidate_patch, local_ok)
    if (.not. local_ok) return
    nvar = size(coarse_state, 1)
    allocate(candidate_coarse_state, mold=coarse_state)
    allocate(candidate_coarse_temperature, mold=coarse_temperature)
    allocate(candidate_fine_state( &
      nvar, candidate_fine_geometry%nx, candidate_fine_geometry%ny))
    allocate(candidate_fine_temperature( &
      candidate_fine_geometry%nx, candidate_fine_geometry%ny))
    if (fine_active) then
      call regrid_two_level_reactive_eb_patch_2d( &
        species, coarse_state, coarse_temperature, coarse_geometry, &
        fine_state, fine_temperature, fine_geometry, patch, &
        candidate_fine_geometry, candidate_patch, candidate_coarse_state, &
        candidate_coarse_temperature, candidate_fine_state, &
        candidate_fine_temperature, local_ok)
    else
      candidate_coarse_state = coarse_state
      candidate_coarse_temperature = coarse_temperature
      call prolong_reactive_eb_patch_pcm_2d( &
        species, coarse_state, coarse_temperature, coarse_geometry, &
        candidate_fine_geometry, candidate_patch, candidate_fine_state, &
        candidate_fine_temperature, local_ok)
    end if
    if (.not. local_ok) return
    call move_alloc(candidate_coarse_state, coarse_state)
    call move_alloc(candidate_coarse_temperature, coarse_temperature)
    call move_alloc(candidate_fine_state, fine_state)
    call move_alloc(candidate_fine_temperature, fine_temperature)
    fine_geometry = candidate_fine_geometry
    patch = candidate_patch
    fine_active = .true.
    changed = .true.
    ok = .true.
  end subroutine regrid_reactive_eb_amr_hierarchy_2d

  subroutine regrid_three_level_reactive_eb_amr_finest_2d( &
      species, config, level_one_state, level_one_temperature, &
      level_one_geometry, level_two_state, level_two_temperature, &
      level_two_geometry, level_one_patch, changed, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_eb_amr_2d_config), intent(in) :: config
    real(dp), allocatable, intent(inout) :: level_one_state(:, :, :)
    real(dp), allocatable, intent(inout) :: level_one_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: level_one_geometry
    real(dp), allocatable, intent(inout) :: level_two_state(:, :, :)
    real(dp), allocatable, intent(inout) :: level_two_temperature(:, :)
    type(eb_geometry_2d), intent(inout) :: level_two_geometry
    type(amr_eb_patch_2d), intent(inout) :: level_one_patch
    logical, intent(out) :: changed, ok

    type(amr_eb_tagging_criteria_2d) :: criteria
    type(amr_eb_regrid_plan_2d) :: plan
    type(eb_geometry_2d) :: candidate_level_two_geometry
    type(amr_eb_patch_2d) :: candidate_level_one_patch
    real(dp), allocatable :: candidate_level_one_state(:, :, :)
    real(dp), allocatable :: candidate_level_one_temperature(:, :)
    real(dp), allocatable :: candidate_level_two_state(:, :, :)
    real(dp), allocatable :: candidate_level_two_temperature(:, :)
    logical, allocatable :: tags(:, :), interior_tags(:, :)
    logical :: local_ok
    integer :: nvar

    changed = .false.
    ok = .false.
    if (.not. supported_three_level_reactive_eb_amr_config(config) .or. &
        .not. config%dynamic_regridding .or. &
        .not. allocated(level_one_state) .or. &
        .not. allocated(level_one_temperature) .or. &
        .not. allocated(level_two_state) .or. &
        .not. allocated(level_two_temperature) .or. &
        .not. level_one_patch%is_valid( &
          level_one_geometry, level_two_geometry)) return

    criteria%relative_gradient_threshold = &
      config%regrid_relative_temperature_gradient
    criteria%absolute_gradient_threshold = &
      config%regrid_absolute_temperature_gradient
    criteria%scale_floor = config%regrid_temperature_scale_floor
    criteria%buffer_cells = config%regrid_buffer_cells
    criteria%minimum_patch_cells_x = config%regrid_minimum_patch_cells_x
    criteria%minimum_patch_cells_y = config%regrid_minimum_patch_cells_y
    allocate(tags(level_one_geometry%nx, level_one_geometry%ny))
    call tag_reactive_eb_temperature_gradient_2d( &
      level_one_temperature, level_one_geometry, criteria, tags, local_ok)
    if (.not. local_ok) return

    allocate(interior_tags( &
      level_one_geometry%nx - 4, level_one_geometry%ny - 4))
    interior_tags = tags(3:level_one_geometry%nx - 2, &
      3:level_one_geometry%ny - 2)
    call build_amr_eb_regrid_plan_2d( &
      interior_tags, criteria, plan, local_ok)
    if (.not. local_ok) return
    if (.not. plan%active) then
      ok = .true.
      return
    end if
    plan%coarse_nx = level_one_geometry%nx
    plan%coarse_ny = level_one_geometry%ny
    plan%tag_i_lower = plan%tag_i_lower + 2
    plan%tag_i_upper = plan%tag_i_upper + 2
    plan%tag_j_lower = plan%tag_j_lower + 2
    plan%tag_j_upper = plan%tag_j_upper + 2
    plan%coarse_i_lower = plan%coarse_i_lower + 2
    plan%coarse_i_upper = plan%coarse_i_upper + 2
    plan%coarse_j_lower = plan%coarse_j_lower + 2
    plan%coarse_j_upper = plan%coarse_j_upper + 2
    if (.not. plan%is_valid() .or. &
        plan%coarse_i_lower < 3 .or. &
        plan%coarse_i_upper > level_one_geometry%nx - 2 .or. &
        plan%coarse_j_lower < 3 .or. &
        plan%coarse_j_upper > level_one_geometry%ny - 2) return
    if (plan%coarse_i_lower == level_one_patch%coarse_i_lower .and. &
        plan%coarse_i_upper == level_one_patch%coarse_i_upper .and. &
        plan%coarse_j_lower == level_one_patch%coarse_j_lower .and. &
        plan%coarse_j_upper == level_one_patch%coarse_j_upper) then
      ok = .true.
      return
    end if

    call build_reactive_eb_amr_patch_geometry_2d( &
      config, level_one_geometry, plan%coarse_i_lower, &
      plan%coarse_i_upper, plan%coarse_j_lower, plan%coarse_j_upper, &
      candidate_level_two_geometry, candidate_level_one_patch, local_ok)
    if (.not. local_ok) return
    nvar = size(level_one_state, 1)
    allocate(candidate_level_one_state, mold=level_one_state)
    allocate(candidate_level_one_temperature, mold=level_one_temperature)
    allocate(candidate_level_two_state( &
      nvar, candidate_level_two_geometry%nx, &
      candidate_level_two_geometry%ny))
    allocate(candidate_level_two_temperature( &
      candidate_level_two_geometry%nx, candidate_level_two_geometry%ny))
    call regrid_two_level_reactive_eb_patch_2d( &
      species, level_one_state, level_one_temperature, level_one_geometry, &
      level_two_state, level_two_temperature, level_two_geometry, &
      level_one_patch, candidate_level_two_geometry, &
      candidate_level_one_patch, candidate_level_one_state, &
      candidate_level_one_temperature, candidate_level_two_state, &
      candidate_level_two_temperature, local_ok)
    if (.not. local_ok) return

    call move_alloc(candidate_level_one_state, level_one_state)
    call move_alloc( &
      candidate_level_one_temperature, level_one_temperature)
    call move_alloc(candidate_level_two_state, level_two_state)
    call move_alloc( &
      candidate_level_two_temperature, level_two_temperature)
    level_two_geometry = candidate_level_two_geometry
    level_one_patch = candidate_level_one_patch
    changed = .true.
    ok = .true.
  end subroutine regrid_three_level_reactive_eb_amr_finest_2d

  pure elemental logical function checkpoint_real_matches(actual, expected) &
      result(matches)
    real(dp), intent(in) :: actual, expected

    matches = .false.
    if (.not. ieee_is_finite(actual) .or. &
        .not. ieee_is_finite(expected)) return
    matches = abs(actual - expected) <= 512.0_dp * epsilon(1.0_dp) * &
      max(tiny(1.0_dp), abs(actual), abs(expected))
  end function checkpoint_real_matches

  subroutine recover_checkpoint_level_temperatures_2d( &
      species, state, temperature, geometry, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :)
    real(dp), intent(inout) :: temperature(:, :)
    type(eb_geometry_2d), intent(in) :: geometry
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: recovered_temperature, sound_speed
    logical :: local_ok
    integer :: i, j

    ok = .false.
    if (.not. geometry%is_valid() .or. &
        size(state, 1) /= reactive_nvar(size(species)) .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        any(shape(temperature) /= [geometry%nx, geometry%ny]) .or. &
        .not. all(ieee_is_finite(state)) .or. &
        .not. all(ieee_is_finite(temperature))) return
    allocate(primitive(reactive_nprim(size(species))))
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) == eb_covered_cell) cycle
        call reactive_conserved_to_primitive( &
          species, state(:, i, j), temperature(i, j), primitive, &
          recovered_temperature, sound_speed, local_ok)
        if (.not. local_ok) return
        temperature(i, j) = recovered_temperature
      end do
    end do
    ok = .true.
  end subroutine recover_checkpoint_level_temperatures_2d

  subroutine write_three_level_checkpoint_field_2d( &
      unit, state, temperature, status)
    integer, intent(in) :: unit
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    integer, intent(out) :: status

    integer :: i, j

    write(unit, '(*(i0,1x))', iostat=status) size(state, 2), size(state, 3)
    if (status /= 0) return
    do j = 1, size(state, 3)
      do i = 1, size(state, 2)
        write(unit, '(*(es27.18e3,1x))', iostat=status) &
          state(:, i, j), temperature(i, j)
        if (status /= 0) return
      end do
    end do
  end subroutine write_three_level_checkpoint_field_2d

  subroutine read_three_level_checkpoint_field_2d( &
      unit, nvar, geometry, state, temperature, status)
    integer, intent(in) :: unit, nvar
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), allocatable, intent(out) :: state(:, :, :)
    real(dp), allocatable, intent(out) :: temperature(:, :)
    integer, intent(out) :: status

    integer :: nx, ny, i, j

    read(unit, *, iostat=status) nx, ny
    if (status /= 0 .or. nx /= geometry%nx .or. ny /= geometry%ny) then
      status = 1
      return
    end if
    allocate(state(nvar, nx, ny), temperature(nx, ny))
    do j = 1, ny
      do i = 1, nx
        read(unit, *, iostat=status) state(:, i, j), temperature(i, j)
        if (status /= 0) return
      end do
    end do
  end subroutine read_three_level_checkpoint_field_2d

  subroutine write_reactive_eb_amr_2d_checkpoint( &
      path, species, config, coarse_state, coarse_temperature, &
      coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
      fine_active, time, steps, regrids, minimum_dt, base_density, ok)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_eb_amr_2d_config), intent(in) :: config
    real(dp), intent(in) :: coarse_state(:, :, :), coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    real(dp), allocatable, intent(in) :: fine_state(:, :, :)
    real(dp), allocatable, intent(in) :: fine_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    logical, intent(in) :: fine_active
    real(dp), intent(in) :: time, minimum_dt, base_density
    integer, intent(in) :: steps, regrids
    logical, intent(out) :: ok

    real(dp) :: time_tolerance
    integer :: unit, status, nvar, i, j, species_index

    ok = .false.
    nvar = reactive_nvar(size(species))
    time_tolerance = 64.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, abs(config%eb%flow%final_time))
    if (len_trim(path) == 0 .or. size(species) < 1 .or. &
        .not. supported_reactive_eb_amr_config(config) .or. &
        config%multipatch_enabled .or. config%three_level_enabled .or. &
        .not. coarse_geometry%is_valid() .or. &
        size(coarse_state, 1) /= nvar .or. &
        size(coarse_state, 2) /= coarse_geometry%nx .or. &
        size(coarse_state, 3) /= coarse_geometry%ny .or. &
        any(shape(coarse_temperature) /= &
          [coarse_geometry%nx, coarse_geometry%ny]) .or. &
        .not. all(ieee_is_finite(coarse_state)) .or. &
        .not. all(ieee_is_finite(coarse_temperature)) .or. &
        .not. ieee_is_finite(time) .or. time <= 0.0_dp .or. &
        time > config%eb%flow%final_time + time_tolerance .or. &
        steps < 1 .or. regrids < 0 .or. &
        real(regrids, dp) > real(steps, dp) + 1.0_dp .or. &
        .not. ieee_is_finite(minimum_dt) .or. minimum_dt <= 0.0_dp .or. &
        minimum_dt > time + time_tolerance .or. &
        .not. ieee_is_finite(base_density) .or. base_density <= 0.0_dp) return
    if (fine_active) then
      if (.not. allocated(fine_state) .or. &
          .not. allocated(fine_temperature)) return
      if (.not. patch%is_valid(coarse_geometry, fine_geometry) .or. &
          size(fine_state, 1) /= nvar .or. &
          size(fine_state, 2) /= fine_geometry%nx .or. &
          size(fine_state, 3) /= fine_geometry%ny .or. &
          any(shape(fine_temperature) /= &
            [fine_geometry%nx, fine_geometry%ny]) .or. &
          .not. all(ieee_is_finite(fine_state)) .or. &
          .not. all(ieee_is_finite(fine_temperature))) return
    else
      if (allocated(fine_state) .or. allocated(fine_temperature) .or. &
          fine_geometry%is_valid() .or. patch%refinement_ratio /= 0) return
    end if

    open(newunit=unit, file=trim(path), status="replace", action="write", &
      form="formatted", iostat=status)
    if (status /= 0) return
    write(unit, '(a)', iostat=status) reactive_eb_amr_checkpoint_magic
    if (status /= 0) go to 900
    write(unit, '(*(i0,1x))', iostat=status) &
      reactive_eb_amr_checkpoint_schema, size(species), nvar, &
      merge(1, 0, fine_active)
    if (status /= 0) go to 900
    do species_index = 1, size(species)
      write(unit, '(a)', iostat=status) trim(species(species_index)%name)
      if (status /= 0) go to 900
    end do
    write(unit, '(a)', iostat=status) trim(config%eb%geometry)
    if (status /= 0) go to 900
    write(unit, '(*(i0,1x))', iostat=status) &
      config%eb%flow%nx, config%eb%flow%ny, config%refinement_ratio
    if (status /= 0) go to 900
    write(unit, '(*(es27.18e3,1x))', iostat=status) &
      config%eb%flow%x_lower, config%eb%flow%x_upper, &
      config%eb%flow%y_lower, config%eb%flow%y_upper
    if (status /= 0) go to 900
    write(unit, '(*(es27.18e3,1x))', iostat=status) &
      config%eb%plane_normal_x, config%eb%plane_normal_y, &
      config%eb%plane_offset, config%eb%circle_center_x, &
      config%eb%circle_center_y, config%eb%circle_radius
    if (status /= 0) go to 900
    write(unit, '(i0)', iostat=status) &
      merge(1, 0, config%eb%circle_fluid_inside)
    if (status /= 0) go to 900
    write(unit, '(a)', iostat=status) trim(config%eb%flow%chemistry_model)
    if (status /= 0) go to 900
    write(unit, '(a)', iostat=status) trim(config%eb%flow%riemann_solver)
    if (status /= 0) go to 900
    write(unit, '(a)', iostat=status) trim(config%eb%flow%reconstruction)
    if (status /= 0) go to 900
    write(unit, '(a)', iostat=status) trim(config%eb%flow%limiter)
    if (status /= 0) go to 900
    write(unit, '(*(i0,1x))', iostat=status) &
      merge(1, 0, config%eb%flow%chemistry_enabled), &
      merge(1, 0, config%dynamic_regridding), &
      merge(1, 0, config%regrid_at_initialization), &
      merge(1, 0, config%remove_fine_patch_when_untagged), &
      config%eb%state_redist_max_order, config%regrid_interval, &
      config%regrid_buffer_cells, config%regrid_minimum_patch_cells_x, &
      config%regrid_minimum_patch_cells_y
    if (status /= 0) go to 900
    write(unit, '(*(es27.18e3,1x))', iostat=status) &
      config%eb%flow%cfl, config%eb%flow%chemistry_relative_tolerance, &
      config%eb%flow%chemistry_absolute_tolerance, &
      config%eb%state_redist_target_volume_fraction, &
      config%regrid_relative_temperature_gradient, &
      config%regrid_absolute_temperature_gradient, &
      config%regrid_temperature_scale_floor
    if (status /= 0) go to 900
    write(unit, '(*(i0,1x))', iostat=status) &
      patch%coarse_i_lower, patch%coarse_i_upper, &
      patch%coarse_j_lower, patch%coarse_j_upper, patch%refinement_ratio
    if (status /= 0) go to 900
    write(unit, '(2(es27.18e3,1x),2(i0,1x),es27.18e3)', iostat=status) &
      time, minimum_dt, steps, regrids, base_density
    if (status /= 0) go to 900
    write(unit, '(*(i0,1x))', iostat=status) &
      coarse_geometry%nx, coarse_geometry%ny
    if (status /= 0) go to 900
    do j = 1, coarse_geometry%ny
      do i = 1, coarse_geometry%nx
        write(unit, '(*(es27.18e3,1x))', iostat=status) &
          coarse_state(:, i, j), coarse_temperature(i, j)
        if (status /= 0) go to 900
      end do
    end do
    if (fine_active) then
      write(unit, '(*(i0,1x))', iostat=status) &
        fine_geometry%nx, fine_geometry%ny
      if (status /= 0) go to 900
      do j = 1, fine_geometry%ny
        do i = 1, fine_geometry%nx
          write(unit, '(*(es27.18e3,1x))', iostat=status) &
            fine_state(:, i, j), fine_temperature(i, j)
          if (status /= 0) go to 900
        end do
      end do
    else
      write(unit, '(*(i0,1x))', iostat=status) 0, 0
      if (status /= 0) go to 900
    end if
    write(unit, '(a)', iostat=status) "END_CHECKPOINT"
    if (status /= 0) go to 900
    close(unit, iostat=status)
    ok = status == 0
    return

900 continue
    close(unit)
  end subroutine write_reactive_eb_amr_2d_checkpoint

  subroutine read_reactive_eb_amr_2d_checkpoint( &
      path, species, config, coarse_state, coarse_temperature, &
      coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
      fine_active, time, steps, regrids, minimum_dt, base_density, ok)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_eb_amr_2d_config), intent(in) :: config
    real(dp), allocatable, intent(out) :: coarse_state(:, :, :)
    real(dp), allocatable, intent(out) :: coarse_temperature(:, :)
    type(eb_geometry_2d), intent(out) :: coarse_geometry
    real(dp), allocatable, intent(out) :: fine_state(:, :, :)
    real(dp), allocatable, intent(out) :: fine_temperature(:, :)
    type(eb_geometry_2d), intent(out) :: fine_geometry
    type(amr_eb_patch_2d), intent(out) :: patch
    logical, intent(out) :: fine_active
    real(dp), intent(out) :: time, minimum_dt, base_density
    integer, intent(out) :: steps, regrids
    logical, intent(out) :: ok

    type(eb_geometry_2d) :: candidate_coarse_geometry
    type(eb_geometry_2d) :: candidate_fine_geometry
    type(amr_eb_patch_2d) :: candidate_patch
    real(dp), allocatable :: candidate_coarse_state(:, :, :)
    real(dp), allocatable :: candidate_coarse_temperature(:, :)
    real(dp), allocatable :: candidate_fine_state(:, :, :)
    real(dp), allocatable :: candidate_fine_temperature(:, :)
    character(len=1024) :: magic, stored_name, stored_geometry
    character(len=1024) :: stored_chemistry_model, stored_solver
    character(len=1024) :: stored_reconstruction, stored_limiter, end_marker
    real(dp) :: stored_domain(4), stored_geometry_values(6)
    real(dp) :: stored_numerics(7), stored_time, stored_minimum_dt
    real(dp) :: stored_base_density, time_tolerance
    integer :: unit, status, schema, stored_species, stored_nvar
    integer :: stored_fine_flag, stored_nx, stored_ny, stored_ratio
    integer :: stored_circle_inside, stored_flags(9), stored_patch(5)
    integer :: stored_coarse_nx, stored_coarse_ny
    integer :: stored_fine_nx, stored_fine_ny
    integer :: stored_steps, stored_regrids, i, j, species_index
    logical :: local_ok

    magic = ""
    stored_name = ""
    stored_geometry = ""
    stored_chemistry_model = ""
    stored_solver = ""
    stored_reconstruction = ""
    stored_limiter = ""
    end_marker = ""
    stored_domain = 0.0_dp
    stored_geometry_values = 0.0_dp
    stored_numerics = 0.0_dp
    stored_time = 0.0_dp
    stored_minimum_dt = 0.0_dp
    stored_base_density = 0.0_dp
    schema = 0
    stored_species = 0
    stored_nvar = 0
    stored_fine_flag = 0
    stored_nx = 0
    stored_ny = 0
    stored_ratio = 0
    stored_circle_inside = 0
    stored_flags = 0
    stored_patch = 0
    stored_coarse_nx = 0
    stored_coarse_ny = 0
    stored_fine_nx = 0
    stored_fine_ny = 0
    stored_steps = 0
    stored_regrids = 0
    coarse_geometry = eb_geometry_2d()
    fine_geometry = eb_geometry_2d()
    patch = amr_eb_patch_2d()
    fine_active = .false.
    time = 0.0_dp
    minimum_dt = 0.0_dp
    base_density = 0.0_dp
    steps = 0
    regrids = 0
    ok = .false.
    if (len_trim(path) == 0 .or. size(species) < 1 .or. &
        .not. supported_reactive_eb_amr_config(config) .or. &
        config%multipatch_enabled .or. config%three_level_enabled) return
    open(newunit=unit, file=trim(path), status="old", action="read", &
      form="formatted", iostat=status)
    if (status /= 0) return
    read(unit, '(a)', iostat=status) magic
    if (status /= 0 .or. &
        trim(magic) /= reactive_eb_amr_checkpoint_magic) go to 900
    read(unit, *, iostat=status) &
      schema, stored_species, stored_nvar, stored_fine_flag
    if (status /= 0 .or. schema /= reactive_eb_amr_checkpoint_schema .or. &
        stored_species /= size(species) .or. &
        stored_nvar /= reactive_nvar(size(species)) .or. &
        (stored_fine_flag /= 0 .and. stored_fine_flag /= 1)) go to 900
    do species_index = 1, stored_species
      read(unit, '(a)', iostat=status) stored_name
      if (status /= 0 .or. &
          trim(stored_name) /= trim(species(species_index)%name)) go to 900
    end do
    read(unit, '(a)', iostat=status) stored_geometry
    if (status /= 0 .or. &
        trim(stored_geometry) /= trim(config%eb%geometry)) go to 900
    read(unit, *, iostat=status) stored_nx, stored_ny, stored_ratio
    if (status /= 0 .or. stored_nx /= config%eb%flow%nx .or. &
        stored_ny /= config%eb%flow%ny .or. &
        stored_ratio /= config%refinement_ratio) go to 900
    read(unit, *, iostat=status) stored_domain
    if (status /= 0 .or. .not. all(checkpoint_real_matches( &
        stored_domain, [config%eb%flow%x_lower, config%eb%flow%x_upper, &
        config%eb%flow%y_lower, config%eb%flow%y_upper]))) go to 900
    read(unit, *, iostat=status) stored_geometry_values
    if (status /= 0 .or. .not. all(checkpoint_real_matches( &
        stored_geometry_values, [config%eb%plane_normal_x, &
        config%eb%plane_normal_y, config%eb%plane_offset, &
        config%eb%circle_center_x, config%eb%circle_center_y, &
        config%eb%circle_radius]))) go to 900
    read(unit, *, iostat=status) stored_circle_inside
    if (status /= 0 .or. stored_circle_inside /= &
        merge(1, 0, config%eb%circle_fluid_inside)) go to 900
    read(unit, '(a)', iostat=status) stored_chemistry_model
    if (status /= 0 .or. trim(stored_chemistry_model) /= &
        trim(config%eb%flow%chemistry_model)) go to 900
    read(unit, '(a)', iostat=status) stored_solver
    if (status /= 0 .or. &
        trim(stored_solver) /= trim(config%eb%flow%riemann_solver)) go to 900
    read(unit, '(a)', iostat=status) stored_reconstruction
    if (status /= 0 .or. trim(stored_reconstruction) /= &
        trim(config%eb%flow%reconstruction)) go to 900
    read(unit, '(a)', iostat=status) stored_limiter
    if (status /= 0 .or. &
        trim(stored_limiter) /= trim(config%eb%flow%limiter)) go to 900
    read(unit, *, iostat=status) stored_flags
    if (status /= 0 .or. any(stored_flags /= [ &
        merge(1, 0, config%eb%flow%chemistry_enabled), &
        merge(1, 0, config%dynamic_regridding), &
        merge(1, 0, config%regrid_at_initialization), &
        merge(1, 0, config%remove_fine_patch_when_untagged), &
        config%eb%state_redist_max_order, config%regrid_interval, &
        config%regrid_buffer_cells, config%regrid_minimum_patch_cells_x, &
        config%regrid_minimum_patch_cells_y])) go to 900
    read(unit, *, iostat=status) stored_numerics
    if (status /= 0 .or. .not. all(checkpoint_real_matches( &
        stored_numerics, [config%eb%flow%cfl, &
        config%eb%flow%chemistry_relative_tolerance, &
        config%eb%flow%chemistry_absolute_tolerance, &
        config%eb%state_redist_target_volume_fraction, &
        config%regrid_relative_temperature_gradient, &
        config%regrid_absolute_temperature_gradient, &
        config%regrid_temperature_scale_floor]))) go to 900
    read(unit, *, iostat=status) stored_patch
    if (status /= 0) go to 900
    read(unit, *, iostat=status) stored_time, stored_minimum_dt, &
      stored_steps, stored_regrids, stored_base_density
    time_tolerance = 64.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, abs(config%eb%flow%final_time))
    if (status /= 0) go to 900
    if (.not. ieee_is_finite(stored_time) .or. &
        .not. ieee_is_finite(stored_minimum_dt) .or. &
        .not. ieee_is_finite(stored_base_density)) go to 900
    if (stored_time <= 0.0_dp .or. &
        stored_time > config%eb%flow%final_time + time_tolerance .or. &
        stored_minimum_dt <= 0.0_dp .or. stored_steps < 1 .or. &
        stored_steps > config%eb%flow%maximum_steps .or. &
        stored_regrids < 0 .or. &
        real(stored_regrids, dp) > real(stored_steps, dp) + 1.0_dp .or. &
        stored_minimum_dt > stored_time + time_tolerance .or. &
        stored_base_density <= 0.0_dp) go to 900

    call build_configured_eb_geometry_2d( &
      config%eb, candidate_coarse_geometry, local_ok)
    if (.not. local_ok) go to 900
    if (stored_fine_flag == 1) then
      if (stored_patch(5) /= config%refinement_ratio) go to 900
      call build_reactive_eb_amr_patch_geometry_2d( &
        config, candidate_coarse_geometry, stored_patch(1), stored_patch(2), &
        stored_patch(3), stored_patch(4), candidate_fine_geometry, &
        candidate_patch, local_ok)
      if (.not. local_ok) go to 900
    else
      if (any(stored_patch /= [0, -1, 0, -1, 0])) go to 900
    end if
    read(unit, *, iostat=status) stored_coarse_nx, stored_coarse_ny
    if (status /= 0 .or. stored_coarse_nx /= candidate_coarse_geometry%nx &
        .or. stored_coarse_ny /= candidate_coarse_geometry%ny) go to 900
    allocate(candidate_coarse_state( &
      stored_nvar, stored_coarse_nx, stored_coarse_ny))
    allocate(candidate_coarse_temperature(stored_coarse_nx, stored_coarse_ny))
    do j = 1, stored_coarse_ny
      do i = 1, stored_coarse_nx
        read(unit, *, iostat=status) &
          candidate_coarse_state(:, i, j), &
          candidate_coarse_temperature(i, j)
        if (status /= 0) go to 900
      end do
    end do
    call recover_checkpoint_level_temperatures_2d( &
      species, candidate_coarse_state, candidate_coarse_temperature, &
      candidate_coarse_geometry, local_ok)
    if (.not. local_ok) go to 900
    read(unit, *, iostat=status) stored_fine_nx, stored_fine_ny
    if (status /= 0) go to 900
    if (stored_fine_flag == 1) then
      if (stored_fine_nx /= candidate_fine_geometry%nx .or. &
          stored_fine_ny /= candidate_fine_geometry%ny) go to 900
      allocate(candidate_fine_state( &
        stored_nvar, stored_fine_nx, stored_fine_ny))
      allocate(candidate_fine_temperature(stored_fine_nx, stored_fine_ny))
      do j = 1, stored_fine_ny
        do i = 1, stored_fine_nx
          read(unit, *, iostat=status) candidate_fine_state(:, i, j), &
            candidate_fine_temperature(i, j)
          if (status /= 0) go to 900
        end do
      end do
      call recover_checkpoint_level_temperatures_2d( &
        species, candidate_fine_state, candidate_fine_temperature, &
        candidate_fine_geometry, local_ok)
      if (.not. local_ok) go to 900
    else
      if (stored_fine_nx /= 0 .or. stored_fine_ny /= 0) go to 900
    end if
    read(unit, '(a)', iostat=status) end_marker
    if (status /= 0 .or. trim(end_marker) /= "END_CHECKPOINT") go to 900
    close(unit, iostat=status)
    if (status /= 0) return

    call move_alloc(candidate_coarse_state, coarse_state)
    call move_alloc(candidate_coarse_temperature, coarse_temperature)
    coarse_geometry = candidate_coarse_geometry
    if (stored_fine_flag == 1) then
      call move_alloc(candidate_fine_state, fine_state)
      call move_alloc(candidate_fine_temperature, fine_temperature)
      fine_geometry = candidate_fine_geometry
      patch = candidate_patch
      fine_active = .true.
    end if
    time = stored_time
    minimum_dt = stored_minimum_dt
    base_density = stored_base_density
    steps = stored_steps
    regrids = stored_regrids
    ok = .true.
    return

900 continue
    close(unit)
  end subroutine read_reactive_eb_amr_2d_checkpoint

  subroutine write_reactive_eb_amr_patch_set_2d_checkpoint( &
      path, species, config, coarse_state, coarse_temperature, &
      coarse_geometry, patch_set, time, steps, regrids, minimum_dt, &
      base_density, ok)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_eb_amr_2d_config), intent(in) :: config
    real(dp), intent(in) :: coarse_state(:, :, :), coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set
    real(dp), intent(in) :: time, minimum_dt, base_density
    integer, intent(in) :: steps, regrids
    logical, intent(out) :: ok

    real(dp) :: time_tolerance
    integer :: unit, status, nvar, child, i, j, species_index

    ok = .false.
    nvar = reactive_nvar(size(species))
    time_tolerance = 64.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, abs(config%eb%flow%final_time))
    if (len_trim(path) == 0 .or. size(species) < 1 .or. &
        .not. supported_reactive_eb_amr_config(config) .or. &
        .not. config%multipatch_enabled .or. config%three_level_enabled .or. &
        .not. coarse_geometry%is_valid() .or. &
        size(coarse_state, 1) /= nvar .or. &
        size(coarse_state, 2) /= coarse_geometry%nx .or. &
        size(coarse_state, 3) /= coarse_geometry%ny .or. &
        any(shape(coarse_temperature) /= &
          [coarse_geometry%nx, coarse_geometry%ny]) .or. &
        .not. all(ieee_is_finite(coarse_state)) .or. &
        .not. all(ieee_is_finite(coarse_temperature)) .or. &
        .not. patch_set%is_valid(coarse_geometry, nvar) .or. &
        .not. ieee_is_finite(time) .or. time <= 0.0_dp .or. &
        time > config%eb%flow%final_time + time_tolerance .or. &
        steps < 1 .or. regrids < 0 .or. &
        real(regrids, dp) > real(steps, dp) + 1.0_dp .or. &
        .not. ieee_is_finite(minimum_dt) .or. minimum_dt <= 0.0_dp .or. &
        minimum_dt > time + time_tolerance .or. &
        .not. ieee_is_finite(base_density) .or. base_density <= 0.0_dp) return

    open(newunit=unit, file=trim(path), status="replace", action="write", &
      form="formatted", iostat=status)
    if (status /= 0) return
    write(unit, '(a)', iostat=status) reactive_eb_patch_set_checkpoint_magic
    if (status /= 0) go to 900
    write(unit, '(*(i0,1x))', iostat=status) &
      reactive_eb_patch_set_checkpoint_schema, size(species), nvar, &
      patch_set%patch_count()
    if (status /= 0) go to 900
    do species_index = 1, size(species)
      write(unit, '(a)', iostat=status) trim(species(species_index)%name)
      if (status /= 0) go to 900
    end do
    write(unit, '(a)', iostat=status) trim(config%eb%geometry)
    if (status /= 0) go to 900
    write(unit, '(*(i0,1x))', iostat=status) &
      config%eb%flow%nx, config%eb%flow%ny, config%refinement_ratio
    if (status /= 0) go to 900
    write(unit, '(*(es27.18e3,1x))', iostat=status) &
      config%eb%flow%x_lower, config%eb%flow%x_upper, &
      config%eb%flow%y_lower, config%eb%flow%y_upper
    if (status /= 0) go to 900
    write(unit, '(*(es27.18e3,1x))', iostat=status) &
      config%eb%plane_normal_x, config%eb%plane_normal_y, &
      config%eb%plane_offset, config%eb%circle_center_x, &
      config%eb%circle_center_y, config%eb%circle_radius
    if (status /= 0) go to 900
    write(unit, '(i0)', iostat=status) &
      merge(1, 0, config%eb%circle_fluid_inside)
    if (status /= 0) go to 900
    write(unit, '(a)', iostat=status) trim(config%eb%flow%chemistry_model)
    if (status /= 0) go to 900
    write(unit, '(a)', iostat=status) trim(config%eb%flow%riemann_solver)
    if (status /= 0) go to 900
    write(unit, '(a)', iostat=status) trim(config%eb%flow%reconstruction)
    if (status /= 0) go to 900
    write(unit, '(a)', iostat=status) trim(config%eb%flow%limiter)
    if (status /= 0) go to 900
    write(unit, '(*(i0,1x))', iostat=status) &
      merge(1, 0, config%eb%flow%chemistry_enabled), &
      merge(1, 0, config%dynamic_regridding), &
      merge(1, 0, config%regrid_at_initialization), &
      merge(1, 0, config%remove_fine_patch_when_untagged), &
      config%eb%state_redist_max_order, config%regrid_interval, &
      config%regrid_buffer_cells, config%regrid_minimum_patch_cells_x, &
      config%regrid_minimum_patch_cells_y, &
      config%regrid_maximum_patch_gap_cells, &
      merge(1, 0, config%multipatch_enabled)
    if (status /= 0) go to 900
    write(unit, '(*(es27.18e3,1x))', iostat=status) &
      config%eb%flow%cfl, config%eb%flow%chemistry_relative_tolerance, &
      config%eb%flow%chemistry_absolute_tolerance, &
      config%eb%state_redist_target_volume_fraction, &
      config%regrid_relative_temperature_gradient, &
      config%regrid_absolute_temperature_gradient, &
      config%regrid_temperature_scale_floor
    if (status /= 0) go to 900
    write(unit, '(2(es27.18e3,1x),2(i0,1x),es27.18e3)', iostat=status) &
      time, minimum_dt, steps, regrids, base_density
    if (status /= 0) go to 900
    write(unit, '(*(i0,1x))', iostat=status) &
      coarse_geometry%nx, coarse_geometry%ny
    if (status /= 0) go to 900
    do j = 1, coarse_geometry%ny
      do i = 1, coarse_geometry%nx
        write(unit, '(*(es27.18e3,1x))', iostat=status) &
          coarse_state(:, i, j), coarse_temperature(i, j)
        if (status /= 0) go to 900
      end do
    end do
    do child = 1, patch_set%patch_count()
      write(unit, '(*(i0,1x))', iostat=status) &
        patch_set%children(child)%patch%coarse_i_lower, &
        patch_set%children(child)%patch%coarse_i_upper, &
        patch_set%children(child)%patch%coarse_j_lower, &
        patch_set%children(child)%patch%coarse_j_upper, &
        patch_set%children(child)%patch%refinement_ratio
      if (status /= 0) go to 900
      write(unit, '(*(i0,1x))', iostat=status) &
        patch_set%children(child)%geometry%nx, &
        patch_set%children(child)%geometry%ny
      if (status /= 0) go to 900
      do j = 1, patch_set%children(child)%geometry%ny
        do i = 1, patch_set%children(child)%geometry%nx
          write(unit, '(*(es27.18e3,1x))', iostat=status) &
            patch_set%children(child)%state(:, i, j), &
            patch_set%children(child)%temperature(i, j)
          if (status /= 0) go to 900
        end do
      end do
    end do
    write(unit, '(a)', iostat=status) "END_CHECKPOINT"
    if (status /= 0) go to 900
    close(unit, iostat=status)
    ok = status == 0
    return

900 continue
    close(unit)
  end subroutine write_reactive_eb_amr_patch_set_2d_checkpoint

  subroutine read_reactive_eb_amr_patch_set_2d_checkpoint( &
      path, species, config, coarse_state, coarse_temperature, &
      coarse_geometry, patch_set, time, steps, regrids, minimum_dt, &
      base_density, ok)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_eb_amr_2d_config), intent(in) :: config
    real(dp), allocatable, intent(out) :: coarse_state(:, :, :)
    real(dp), allocatable, intent(out) :: coarse_temperature(:, :)
    type(eb_geometry_2d), intent(out) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(out) :: patch_set
    real(dp), intent(out) :: time, minimum_dt, base_density
    integer, intent(out) :: steps, regrids
    logical, intent(out) :: ok

    type(eb_geometry_2d) :: candidate_coarse_geometry
    type(reactive_eb_patch_set_2d) :: candidate_set
    real(dp), allocatable :: candidate_coarse_state(:, :, :)
    real(dp), allocatable :: candidate_coarse_temperature(:, :)
    character(len=1024) :: magic, stored_name, stored_geometry
    character(len=1024) :: stored_chemistry_model, stored_solver
    character(len=1024) :: stored_reconstruction, stored_limiter, end_marker
    real(dp) :: stored_domain(4), stored_geometry_values(6)
    real(dp) :: stored_numerics(7), stored_time, stored_minimum_dt
    real(dp) :: stored_base_density, time_tolerance
    integer :: unit, status, schema, stored_species, stored_nvar
    integer :: stored_patch_count, stored_nx, stored_ny, stored_ratio
    integer :: stored_circle_inside, stored_flags(11), stored_patch(5)
    integer :: stored_coarse_nx, stored_coarse_ny
    integer :: stored_fine_nx, stored_fine_ny
    integer :: stored_steps, stored_regrids, child, i, j, species_index
    logical :: local_ok

    magic = ""
    stored_name = ""
    stored_geometry = ""
    stored_chemistry_model = ""
    stored_solver = ""
    stored_reconstruction = ""
    stored_limiter = ""
    end_marker = ""
    stored_domain = 0.0_dp
    stored_geometry_values = 0.0_dp
    stored_numerics = 0.0_dp
    stored_time = 0.0_dp
    stored_minimum_dt = 0.0_dp
    stored_base_density = 0.0_dp
    schema = 0
    stored_species = 0
    stored_nvar = 0
    stored_patch_count = 0
    stored_nx = 0
    stored_ny = 0
    stored_ratio = 0
    stored_circle_inside = 0
    stored_flags = 0
    stored_patch = 0
    stored_coarse_nx = 0
    stored_coarse_ny = 0
    stored_fine_nx = 0
    stored_fine_ny = 0
    stored_steps = 0
    stored_regrids = 0
    coarse_geometry = eb_geometry_2d()
    patch_set = reactive_eb_patch_set_2d()
    time = 0.0_dp
    minimum_dt = 0.0_dp
    base_density = 0.0_dp
    steps = 0
    regrids = 0
    ok = .false.
    if (len_trim(path) == 0 .or. size(species) < 1 .or. &
        .not. supported_reactive_eb_amr_config(config) .or. &
        .not. config%multipatch_enabled .or. &
        config%three_level_enabled) return
    open(newunit=unit, file=trim(path), status="old", action="read", &
      form="formatted", iostat=status)
    if (status /= 0) return
    read(unit, '(a)', iostat=status) magic
    if (status /= 0 .or. &
        trim(magic) /= reactive_eb_patch_set_checkpoint_magic) go to 900
    read(unit, *, iostat=status) &
      schema, stored_species, stored_nvar, stored_patch_count
    if (status /= 0 .or. &
        schema /= reactive_eb_patch_set_checkpoint_schema .or. &
        stored_species /= size(species) .or. &
        stored_nvar /= reactive_nvar(size(species)) .or. &
        stored_patch_count < 0) go to 900
    do species_index = 1, stored_species
      read(unit, '(a)', iostat=status) stored_name
      if (status /= 0 .or. &
          trim(stored_name) /= trim(species(species_index)%name)) go to 900
    end do
    read(unit, '(a)', iostat=status) stored_geometry
    if (status /= 0 .or. &
        trim(stored_geometry) /= trim(config%eb%geometry)) go to 900
    read(unit, *, iostat=status) stored_nx, stored_ny, stored_ratio
    if (status /= 0 .or. stored_nx /= config%eb%flow%nx .or. &
        stored_ny /= config%eb%flow%ny .or. &
        stored_ratio /= config%refinement_ratio) go to 900
    read(unit, *, iostat=status) stored_domain
    if (status /= 0 .or. .not. all(checkpoint_real_matches( &
        stored_domain, [config%eb%flow%x_lower, config%eb%flow%x_upper, &
        config%eb%flow%y_lower, config%eb%flow%y_upper]))) go to 900
    read(unit, *, iostat=status) stored_geometry_values
    if (status /= 0 .or. .not. all(checkpoint_real_matches( &
        stored_geometry_values, [config%eb%plane_normal_x, &
        config%eb%plane_normal_y, config%eb%plane_offset, &
        config%eb%circle_center_x, config%eb%circle_center_y, &
        config%eb%circle_radius]))) go to 900
    read(unit, *, iostat=status) stored_circle_inside
    if (status /= 0 .or. stored_circle_inside /= &
        merge(1, 0, config%eb%circle_fluid_inside)) go to 900
    read(unit, '(a)', iostat=status) stored_chemistry_model
    if (status /= 0 .or. trim(stored_chemistry_model) /= &
        trim(config%eb%flow%chemistry_model)) go to 900
    read(unit, '(a)', iostat=status) stored_solver
    if (status /= 0 .or. &
        trim(stored_solver) /= trim(config%eb%flow%riemann_solver)) go to 900
    read(unit, '(a)', iostat=status) stored_reconstruction
    if (status /= 0 .or. trim(stored_reconstruction) /= &
        trim(config%eb%flow%reconstruction)) go to 900
    read(unit, '(a)', iostat=status) stored_limiter
    if (status /= 0 .or. &
        trim(stored_limiter) /= trim(config%eb%flow%limiter)) go to 900
    read(unit, *, iostat=status) stored_flags
    if (status /= 0 .or. any(stored_flags /= [ &
        merge(1, 0, config%eb%flow%chemistry_enabled), &
        merge(1, 0, config%dynamic_regridding), &
        merge(1, 0, config%regrid_at_initialization), &
        merge(1, 0, config%remove_fine_patch_when_untagged), &
        config%eb%state_redist_max_order, config%regrid_interval, &
        config%regrid_buffer_cells, config%regrid_minimum_patch_cells_x, &
        config%regrid_minimum_patch_cells_y, &
        config%regrid_maximum_patch_gap_cells, &
        merge(1, 0, config%multipatch_enabled)])) go to 900
    read(unit, *, iostat=status) stored_numerics
    if (status /= 0 .or. .not. all(checkpoint_real_matches( &
        stored_numerics, [config%eb%flow%cfl, &
        config%eb%flow%chemistry_relative_tolerance, &
        config%eb%flow%chemistry_absolute_tolerance, &
        config%eb%state_redist_target_volume_fraction, &
        config%regrid_relative_temperature_gradient, &
        config%regrid_absolute_temperature_gradient, &
        config%regrid_temperature_scale_floor]))) go to 900
    read(unit, *, iostat=status) stored_time, stored_minimum_dt, &
      stored_steps, stored_regrids, stored_base_density
    time_tolerance = 64.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, abs(config%eb%flow%final_time))
    if (status /= 0 .or. .not. ieee_is_finite(stored_time) .or. &
        .not. ieee_is_finite(stored_minimum_dt) .or. &
        .not. ieee_is_finite(stored_base_density)) go to 900
    if (stored_time <= 0.0_dp .or. &
        stored_time > config%eb%flow%final_time + time_tolerance .or. &
        stored_minimum_dt <= 0.0_dp .or. stored_steps < 1 .or. &
        stored_steps > config%eb%flow%maximum_steps .or. &
        stored_regrids < 0 .or. &
        real(stored_regrids, dp) > real(stored_steps, dp) + 1.0_dp .or. &
        stored_minimum_dt > stored_time + time_tolerance .or. &
        stored_base_density <= 0.0_dp) go to 900

    call build_configured_eb_geometry_2d( &
      config%eb, candidate_coarse_geometry, local_ok)
    if (.not. local_ok .or. stored_patch_count > &
        candidate_coarse_geometry%nx * candidate_coarse_geometry%ny) go to 900
    read(unit, *, iostat=status) stored_coarse_nx, stored_coarse_ny
    if (status /= 0 .or. stored_coarse_nx /= candidate_coarse_geometry%nx &
        .or. stored_coarse_ny /= candidate_coarse_geometry%ny) go to 900
    allocate(candidate_coarse_state( &
      stored_nvar, stored_coarse_nx, stored_coarse_ny))
    allocate(candidate_coarse_temperature(stored_coarse_nx, stored_coarse_ny))
    do j = 1, stored_coarse_ny
      do i = 1, stored_coarse_nx
        read(unit, *, iostat=status) candidate_coarse_state(:, i, j), &
          candidate_coarse_temperature(i, j)
        if (status /= 0) go to 900
      end do
    end do
    call recover_checkpoint_level_temperatures_2d( &
      species, candidate_coarse_state, candidate_coarse_temperature, &
      candidate_coarse_geometry, local_ok)
    if (.not. local_ok) go to 900

    allocate(candidate_set%children(stored_patch_count))
    do child = 1, stored_patch_count
      read(unit, *, iostat=status) stored_patch
      if (status /= 0 .or. stored_patch(5) /= config%refinement_ratio) &
        go to 900
      call build_reactive_eb_amr_patch_geometry_2d( &
        config, candidate_coarse_geometry, stored_patch(1), stored_patch(2), &
        stored_patch(3), stored_patch(4), &
        candidate_set%children(child)%geometry, &
        candidate_set%children(child)%patch, local_ok)
      if (.not. local_ok) go to 900
      read(unit, *, iostat=status) stored_fine_nx, stored_fine_ny
      if (status /= 0 .or. &
          stored_fine_nx /= candidate_set%children(child)%geometry%nx .or. &
          stored_fine_ny /= candidate_set%children(child)%geometry%ny) &
        go to 900
      allocate(candidate_set%children(child)%state( &
        stored_nvar, stored_fine_nx, stored_fine_ny))
      allocate(candidate_set%children(child)%temperature( &
        stored_fine_nx, stored_fine_ny))
      do j = 1, stored_fine_ny
        do i = 1, stored_fine_nx
          read(unit, *, iostat=status) &
            candidate_set%children(child)%state(:, i, j), &
            candidate_set%children(child)%temperature(i, j)
          if (status /= 0) go to 900
        end do
      end do
      call recover_checkpoint_level_temperatures_2d( &
        species, candidate_set%children(child)%state, &
        candidate_set%children(child)%temperature, &
        candidate_set%children(child)%geometry, local_ok)
      if (.not. local_ok) go to 900
    end do
    if (.not. candidate_set%is_valid( &
        candidate_coarse_geometry, stored_nvar)) go to 900
    read(unit, '(a)', iostat=status) end_marker
    if (status /= 0 .or. trim(end_marker) /= "END_CHECKPOINT") go to 900
    close(unit, iostat=status)
    if (status /= 0) return

    call move_alloc(candidate_coarse_state, coarse_state)
    call move_alloc(candidate_coarse_temperature, coarse_temperature)
    coarse_geometry = candidate_coarse_geometry
    patch_set = candidate_set
    time = stored_time
    minimum_dt = stored_minimum_dt
    base_density = stored_base_density
    steps = stored_steps
    regrids = stored_regrids
    ok = .true.
    return

900 continue
    close(unit)
  end subroutine read_reactive_eb_amr_patch_set_2d_checkpoint

  subroutine write_reactive_eb_amr_three_level_2d_checkpoint( &
      path, species, config, root_state, root_temperature, root_geometry, &
      level_one_state, level_one_temperature, level_one_geometry, &
      root_patch, level_two_state, level_two_temperature, &
      level_two_geometry, level_one_patch, time, steps, minimum_dt, &
      base_density, ok, regrids)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_eb_amr_2d_config), intent(in) :: config
    real(dp), intent(in) :: root_state(:, :, :), root_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: root_geometry
    real(dp), intent(in) :: level_one_state(:, :, :)
    real(dp), intent(in) :: level_one_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: level_one_geometry
    type(amr_eb_patch_2d), intent(in) :: root_patch
    real(dp), intent(in) :: level_two_state(:, :, :)
    real(dp), intent(in) :: level_two_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: level_two_geometry
    type(amr_eb_patch_2d), intent(in) :: level_one_patch
    real(dp), intent(in) :: time, minimum_dt, base_density
    integer, intent(in) :: steps
    logical, intent(out) :: ok
    integer, intent(in), optional :: regrids

    real(dp) :: time_tolerance
    integer :: unit, status, nvar, species_index, checkpoint_regrids

    ok = .false.
    checkpoint_regrids = 0
    if (present(regrids)) checkpoint_regrids = regrids
    nvar = reactive_nvar(size(species))
    time_tolerance = 64.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, abs(config%eb%flow%final_time))
    if (len_trim(path) == 0 .or. size(species) < 1 .or. &
        .not. supported_three_level_reactive_eb_amr_config(config) .or. &
        .not. root_patch%is_valid(root_geometry, level_one_geometry) .or. &
        .not. level_one_patch%is_valid( &
          level_one_geometry, level_two_geometry) .or. &
        root_patch%coarse_i_lower /= config%coarse_i_lower .or. &
        root_patch%coarse_i_upper /= config%coarse_i_upper .or. &
        root_patch%coarse_j_lower /= config%coarse_j_lower .or. &
        root_patch%coarse_j_upper /= config%coarse_j_upper .or. &
        size(root_state, 1) /= nvar .or. &
        size(root_state, 2) /= root_geometry%nx .or. &
        size(root_state, 3) /= root_geometry%ny .or. &
        size(level_one_state, 1) /= nvar .or. &
        size(level_one_state, 2) /= level_one_geometry%nx .or. &
        size(level_one_state, 3) /= level_one_geometry%ny .or. &
        size(level_two_state, 1) /= nvar .or. &
        size(level_two_state, 2) /= level_two_geometry%nx .or. &
        size(level_two_state, 3) /= level_two_geometry%ny .or. &
        any(shape(root_temperature) /= &
          [root_geometry%nx, root_geometry%ny]) .or. &
        any(shape(level_one_temperature) /= &
          [level_one_geometry%nx, level_one_geometry%ny]) .or. &
        any(shape(level_two_temperature) /= &
          [level_two_geometry%nx, level_two_geometry%ny]) .or. &
        .not. all(ieee_is_finite(root_state)) .or. &
        .not. all(ieee_is_finite(level_one_state)) .or. &
        .not. all(ieee_is_finite(level_two_state)) .or. &
        .not. all(ieee_is_finite(root_temperature)) .or. &
        .not. all(ieee_is_finite(level_one_temperature)) .or. &
        .not. all(ieee_is_finite(level_two_temperature)) .or. &
        .not. ieee_is_finite(time) .or. time <= 0.0_dp .or. &
        time > config%eb%flow%final_time + time_tolerance .or. &
        steps < 1 .or. steps > config%eb%flow%maximum_steps .or. &
        .not. ieee_is_finite(minimum_dt) .or. minimum_dt <= 0.0_dp .or. &
        minimum_dt > time + time_tolerance .or. &
        .not. ieee_is_finite(base_density) .or. base_density <= 0.0_dp) return
    if (config%dynamic_regridding) then
      if (.not. present(regrids) .or. checkpoint_regrids < 0 .or. &
          checkpoint_regrids > steps + 1 .or. &
          level_one_patch%coarse_i_lower < 3 .or. &
          level_one_patch%coarse_i_upper > level_one_geometry%nx - 2 .or. &
          level_one_patch%coarse_j_lower < 3 .or. &
          level_one_patch%coarse_j_upper > level_one_geometry%ny - 2) return
    else
      if (level_one_patch%coarse_i_lower /= config%level_two_i_lower .or. &
          level_one_patch%coarse_i_upper /= config%level_two_i_upper .or. &
          level_one_patch%coarse_j_lower /= config%level_two_j_lower .or. &
          level_one_patch%coarse_j_upper /= config%level_two_j_upper) return
    end if

    open(newunit=unit, file=trim(path), status="replace", action="write", &
      form="formatted", iostat=status)
    if (status /= 0) return
    if (config%dynamic_regridding) then
      write(unit, '(a)', iostat=status) &
        reactive_eb_dynamic_three_level_checkpoint_magic
    else
      write(unit, '(a)', iostat=status) reactive_eb_three_level_checkpoint_magic
    end if
    if (status /= 0) go to 900
    if (config%dynamic_regridding) then
      write(unit, '(*(i0,1x))', iostat=status) &
        reactive_eb_dynamic_three_level_checkpoint_schema, &
        size(species), nvar
    else
      write(unit, '(*(i0,1x))', iostat=status) &
        reactive_eb_three_level_checkpoint_schema, size(species), nvar
    end if
    if (status /= 0) go to 900
    do species_index = 1, size(species)
      write(unit, '(a)', iostat=status) trim(species(species_index)%name)
      if (status /= 0) go to 900
    end do
    write(unit, '(a)', iostat=status) trim(config%eb%geometry)
    if (status /= 0) go to 900
    write(unit, '(*(i0,1x))', iostat=status) &
      config%eb%flow%nx, config%eb%flow%ny, config%refinement_ratio
    if (status /= 0) go to 900
    write(unit, '(*(es27.18e3,1x))', iostat=status) &
      config%eb%flow%x_lower, config%eb%flow%x_upper, &
      config%eb%flow%y_lower, config%eb%flow%y_upper
    if (status /= 0) go to 900
    write(unit, '(*(es27.18e3,1x))', iostat=status) &
      config%eb%plane_normal_x, config%eb%plane_normal_y, &
      config%eb%plane_offset, config%eb%circle_center_x, &
      config%eb%circle_center_y, config%eb%circle_radius
    if (status /= 0) go to 900
    write(unit, '(i0)', iostat=status) &
      merge(1, 0, config%eb%circle_fluid_inside)
    if (status /= 0) go to 900
    write(unit, '(a)', iostat=status) trim(config%eb%flow%chemistry_model)
    if (status /= 0) go to 900
    write(unit, '(a)', iostat=status) trim(config%eb%flow%riemann_solver)
    if (status /= 0) go to 900
    write(unit, '(a)', iostat=status) trim(config%eb%flow%reconstruction)
    if (status /= 0) go to 900
    write(unit, '(a)', iostat=status) trim(config%eb%flow%limiter)
    if (status /= 0) go to 900
    write(unit, '(*(i0,1x))', iostat=status) &
      merge(1, 0, config%eb%flow%chemistry_enabled), &
      config%eb%state_redist_max_order
    if (status /= 0) go to 900
    write(unit, '(*(es27.18e3,1x))', iostat=status) &
      config%eb%flow%cfl, &
      config%eb%flow%chemistry_relative_tolerance, &
      config%eb%flow%chemistry_absolute_tolerance, &
      config%eb%state_redist_target_volume_fraction
    if (status /= 0) go to 900
    if (config%dynamic_regridding) then
      write(unit, '(*(i0,1x))', iostat=status) &
        merge(1, 0, config%regrid_at_initialization), &
        config%regrid_interval, config%regrid_buffer_cells, &
        config%regrid_minimum_patch_cells_x, &
        config%regrid_minimum_patch_cells_y, &
        config%regrid_maximum_patch_gap_cells
      if (status /= 0) go to 900
      write(unit, '(*(es27.18e3,1x))', iostat=status) &
        config%regrid_relative_temperature_gradient, &
        config%regrid_absolute_temperature_gradient, &
        config%regrid_temperature_scale_floor
      if (status /= 0) go to 900
    end if
    write(unit, '(*(i0,1x))', iostat=status) &
      root_patch%coarse_i_lower, root_patch%coarse_i_upper, &
      root_patch%coarse_j_lower, root_patch%coarse_j_upper, &
      root_patch%refinement_ratio
    if (status /= 0) go to 900
    write(unit, '(*(i0,1x))', iostat=status) &
      level_one_patch%coarse_i_lower, level_one_patch%coarse_i_upper, &
      level_one_patch%coarse_j_lower, level_one_patch%coarse_j_upper, &
      level_one_patch%refinement_ratio
    if (status /= 0) go to 900
    if (config%dynamic_regridding) then
      write(unit, '(2(es27.18e3,1x),2(i0,1x),es27.18e3)', &
        iostat=status) time, minimum_dt, steps, checkpoint_regrids, &
        base_density
    else
      write(unit, '(2(es27.18e3,1x),i0,1x,es27.18e3)', iostat=status) &
        time, minimum_dt, steps, base_density
    end if
    if (status /= 0) go to 900
    call write_three_level_checkpoint_field_2d( &
      unit, root_state, root_temperature, status)
    if (status /= 0) go to 900
    call write_three_level_checkpoint_field_2d( &
      unit, level_one_state, level_one_temperature, status)
    if (status /= 0) go to 900
    call write_three_level_checkpoint_field_2d( &
      unit, level_two_state, level_two_temperature, status)
    if (status /= 0) go to 900
    write(unit, '(a)', iostat=status) "END_CHECKPOINT"
    if (status /= 0) go to 900
    close(unit, iostat=status)
    ok = status == 0
    return

900 continue
    close(unit)
  end subroutine write_reactive_eb_amr_three_level_2d_checkpoint

  subroutine read_reactive_eb_amr_three_level_2d_checkpoint( &
      path, species, config, root_state, root_temperature, root_geometry, &
      level_one_state, level_one_temperature, level_one_geometry, &
      root_patch, level_two_state, level_two_temperature, &
      level_two_geometry, level_one_patch, time, steps, minimum_dt, &
      base_density, ok, regrids)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_eb_amr_2d_config), intent(in) :: config
    real(dp), allocatable, intent(out) :: root_state(:, :, :)
    real(dp), allocatable, intent(out) :: root_temperature(:, :)
    type(eb_geometry_2d), intent(out) :: root_geometry
    real(dp), allocatable, intent(out) :: level_one_state(:, :, :)
    real(dp), allocatable, intent(out) :: level_one_temperature(:, :)
    type(eb_geometry_2d), intent(out) :: level_one_geometry
    type(amr_eb_patch_2d), intent(out) :: root_patch
    real(dp), allocatable, intent(out) :: level_two_state(:, :, :)
    real(dp), allocatable, intent(out) :: level_two_temperature(:, :)
    type(eb_geometry_2d), intent(out) :: level_two_geometry
    type(amr_eb_patch_2d), intent(out) :: level_one_patch
    real(dp), intent(out) :: time, minimum_dt, base_density
    integer, intent(out) :: steps
    logical, intent(out) :: ok
    integer, intent(out), optional :: regrids

    type(eb_geometry_2d) :: candidate_root_geometry
    type(eb_geometry_2d) :: candidate_level_one_geometry
    type(eb_geometry_2d) :: candidate_level_two_geometry
    type(amr_eb_patch_2d) :: candidate_root_patch
    type(amr_eb_patch_2d) :: candidate_level_one_patch
    real(dp), allocatable :: candidate_root_state(:, :, :)
    real(dp), allocatable :: candidate_root_temperature(:, :)
    real(dp), allocatable :: candidate_level_one_state(:, :, :)
    real(dp), allocatable :: candidate_level_one_temperature(:, :)
    real(dp), allocatable :: candidate_level_two_state(:, :, :)
    real(dp), allocatable :: candidate_level_two_temperature(:, :)
    character(len=1024) :: magic, stored_name, stored_geometry
    character(len=1024) :: stored_chemistry_model, stored_solver
    character(len=1024) :: stored_reconstruction, stored_limiter, end_marker
    real(dp) :: stored_domain(4), stored_geometry_values(6)
    real(dp) :: stored_numerics(4), stored_time, stored_minimum_dt
    real(dp) :: stored_base_density, time_tolerance
    real(dp) :: stored_regrid_numerics(3)
    integer :: unit, status, schema, stored_species, stored_nvar
    integer :: stored_nx, stored_ny, stored_ratio, stored_circle_inside
    integer :: stored_flags(2), stored_root_patch(5), stored_level_one_patch(5)
    integer :: stored_regrid_controls(6)
    integer :: stored_steps, stored_regrids, species_index
    logical :: local_ok, dynamic_checkpoint

    root_geometry = eb_geometry_2d()
    level_one_geometry = eb_geometry_2d()
    level_two_geometry = eb_geometry_2d()
    root_patch = amr_eb_patch_2d()
    level_one_patch = amr_eb_patch_2d()
    time = 0.0_dp
    minimum_dt = 0.0_dp
    base_density = 0.0_dp
    steps = 0
    stored_regrids = 0
    if (present(regrids)) regrids = 0
    ok = .false.
    if (len_trim(path) == 0 .or. size(species) < 1 .or. &
        .not. supported_three_level_reactive_eb_amr_config(config)) return
    open(newunit=unit, file=trim(path), status="old", action="read", &
      form="formatted", iostat=status)
    if (status /= 0) return
    read(unit, '(a)', iostat=status) magic
    if (status /= 0) go to 900
    dynamic_checkpoint = &
      trim(magic) == reactive_eb_dynamic_three_level_checkpoint_magic
    if (.not. dynamic_checkpoint .and. &
        trim(magic) /= reactive_eb_three_level_checkpoint_magic) go to 900
    if (dynamic_checkpoint .neqv. config%dynamic_regridding) go to 900
    if (dynamic_checkpoint .and. .not. present(regrids)) go to 900
    read(unit, *, iostat=status) &
      schema, stored_species, stored_nvar
    if (status /= 0 .or. &
        (dynamic_checkpoint .and. &
         schema /= reactive_eb_dynamic_three_level_checkpoint_schema) .or. &
        (.not. dynamic_checkpoint .and. &
         schema /= reactive_eb_three_level_checkpoint_schema) .or. &
        stored_species /= size(species) .or. &
        stored_nvar /= reactive_nvar(size(species))) go to 900
    do species_index = 1, stored_species
      read(unit, '(a)', iostat=status) stored_name
      if (status /= 0 .or. &
          trim(stored_name) /= trim(species(species_index)%name)) go to 900
    end do
    read(unit, '(a)', iostat=status) stored_geometry
    if (status /= 0 .or. &
        trim(stored_geometry) /= trim(config%eb%geometry)) go to 900
    read(unit, *, iostat=status) stored_nx, stored_ny, stored_ratio
    if (status /= 0 .or. stored_nx /= config%eb%flow%nx .or. &
        stored_ny /= config%eb%flow%ny .or. &
        stored_ratio /= config%refinement_ratio) go to 900
    read(unit, *, iostat=status) stored_domain
    if (status /= 0 .or. .not. all(checkpoint_real_matches( &
        stored_domain, [config%eb%flow%x_lower, config%eb%flow%x_upper, &
        config%eb%flow%y_lower, config%eb%flow%y_upper]))) go to 900
    read(unit, *, iostat=status) stored_geometry_values
    if (status /= 0 .or. .not. all(checkpoint_real_matches( &
        stored_geometry_values, [config%eb%plane_normal_x, &
        config%eb%plane_normal_y, config%eb%plane_offset, &
        config%eb%circle_center_x, config%eb%circle_center_y, &
        config%eb%circle_radius]))) go to 900
    read(unit, *, iostat=status) stored_circle_inside
    if (status /= 0 .or. stored_circle_inside /= &
        merge(1, 0, config%eb%circle_fluid_inside)) go to 900
    read(unit, '(a)', iostat=status) stored_chemistry_model
    if (status /= 0 .or. trim(stored_chemistry_model) /= &
        trim(config%eb%flow%chemistry_model)) go to 900
    read(unit, '(a)', iostat=status) stored_solver
    if (status /= 0 .or. &
        trim(stored_solver) /= trim(config%eb%flow%riemann_solver)) go to 900
    read(unit, '(a)', iostat=status) stored_reconstruction
    if (status /= 0 .or. trim(stored_reconstruction) /= &
        trim(config%eb%flow%reconstruction)) go to 900
    read(unit, '(a)', iostat=status) stored_limiter
    if (status /= 0 .or. &
        trim(stored_limiter) /= trim(config%eb%flow%limiter)) go to 900
    read(unit, *, iostat=status) stored_flags
    if (status /= 0 .or. any(stored_flags /= [ &
        merge(1, 0, config%eb%flow%chemistry_enabled), &
        config%eb%state_redist_max_order])) go to 900
    read(unit, *, iostat=status) stored_numerics
    if (status /= 0 .or. .not. all(checkpoint_real_matches( &
        stored_numerics, [config%eb%flow%cfl, &
        config%eb%flow%chemistry_relative_tolerance, &
        config%eb%flow%chemistry_absolute_tolerance, &
        config%eb%state_redist_target_volume_fraction]))) go to 900
    if (dynamic_checkpoint) then
      read(unit, *, iostat=status) stored_regrid_controls
      if (status /= 0 .or. any(stored_regrid_controls /= [ &
          merge(1, 0, config%regrid_at_initialization), &
          config%regrid_interval, config%regrid_buffer_cells, &
          config%regrid_minimum_patch_cells_x, &
          config%regrid_minimum_patch_cells_y, &
          config%regrid_maximum_patch_gap_cells])) go to 900
      read(unit, *, iostat=status) stored_regrid_numerics
      if (status /= 0 .or. .not. all(checkpoint_real_matches( &
          stored_regrid_numerics, [ &
          config%regrid_relative_temperature_gradient, &
          config%regrid_absolute_temperature_gradient, &
          config%regrid_temperature_scale_floor]))) go to 900
    end if
    read(unit, *, iostat=status) stored_root_patch
    if (status /= 0 .or. any(stored_root_patch /= [ &
        config%coarse_i_lower, config%coarse_i_upper, &
        config%coarse_j_lower, config%coarse_j_upper, &
        config%refinement_ratio])) go to 900
    read(unit, *, iostat=status) stored_level_one_patch
    if (status /= 0) go to 900
    if (dynamic_checkpoint) then
      if (stored_level_one_patch(1) < 3 .or. &
          stored_level_one_patch(2) > &
            (config%coarse_i_upper - config%coarse_i_lower + 1) * &
              config%refinement_ratio - 2 .or. &
          stored_level_one_patch(3) < 3 .or. &
          stored_level_one_patch(4) > &
            (config%coarse_j_upper - config%coarse_j_lower + 1) * &
              config%refinement_ratio - 2 .or. &
          stored_level_one_patch(2) < stored_level_one_patch(1) .or. &
          stored_level_one_patch(4) < stored_level_one_patch(3) .or. &
          stored_level_one_patch(5) /= config%refinement_ratio) go to 900
      read(unit, *, iostat=status) stored_time, stored_minimum_dt, &
        stored_steps, stored_regrids, stored_base_density
    else
      if (any(stored_level_one_patch /= [ &
          config%level_two_i_lower, config%level_two_i_upper, &
          config%level_two_j_lower, config%level_two_j_upper, &
          config%refinement_ratio])) go to 900
      read(unit, *, iostat=status) stored_time, stored_minimum_dt, &
        stored_steps, stored_base_density
    end if
    time_tolerance = 64.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, abs(config%eb%flow%final_time))
    if (status /= 0 .or. .not. ieee_is_finite(stored_time) .or. &
        .not. ieee_is_finite(stored_minimum_dt) .or. &
        .not. ieee_is_finite(stored_base_density) .or. &
        stored_time <= 0.0_dp .or. &
        stored_time > config%eb%flow%final_time + time_tolerance .or. &
        stored_minimum_dt <= 0.0_dp .or. &
        stored_minimum_dt > stored_time + time_tolerance .or. &
        stored_steps < 1 .or. &
        stored_steps > config%eb%flow%maximum_steps .or. &
        stored_regrids < 0 .or. stored_regrids > stored_steps + 1 .or. &
        stored_base_density <= 0.0_dp) go to 900

    call build_configured_eb_geometry_2d( &
      config%eb, candidate_root_geometry, local_ok)
    if (.not. local_ok) go to 900
    call build_reactive_eb_amr_patch_geometry_2d( &
      config, candidate_root_geometry, config%coarse_i_lower, &
      config%coarse_i_upper, config%coarse_j_lower, &
      config%coarse_j_upper, candidate_level_one_geometry, &
      candidate_root_patch, local_ok)
    if (.not. local_ok) go to 900
    call build_reactive_eb_amr_patch_geometry_2d( &
      config, candidate_level_one_geometry, stored_level_one_patch(1), &
      stored_level_one_patch(2), stored_level_one_patch(3), &
      stored_level_one_patch(4), candidate_level_two_geometry, &
      candidate_level_one_patch, local_ok)
    if (.not. local_ok) go to 900
    call read_three_level_checkpoint_field_2d( &
      unit, stored_nvar, candidate_root_geometry, candidate_root_state, &
      candidate_root_temperature, status)
    if (status /= 0) go to 900
    call read_three_level_checkpoint_field_2d( &
      unit, stored_nvar, candidate_level_one_geometry, &
      candidate_level_one_state, candidate_level_one_temperature, status)
    if (status /= 0) go to 900
    call read_three_level_checkpoint_field_2d( &
      unit, stored_nvar, candidate_level_two_geometry, &
      candidate_level_two_state, candidate_level_two_temperature, status)
    if (status /= 0) go to 900
    call recover_checkpoint_level_temperatures_2d( &
      species, candidate_root_state, candidate_root_temperature, &
      candidate_root_geometry, local_ok)
    if (.not. local_ok) go to 900
    call recover_checkpoint_level_temperatures_2d( &
      species, candidate_level_one_state, candidate_level_one_temperature, &
      candidate_level_one_geometry, local_ok)
    if (.not. local_ok) go to 900
    call recover_checkpoint_level_temperatures_2d( &
      species, candidate_level_two_state, candidate_level_two_temperature, &
      candidate_level_two_geometry, local_ok)
    if (.not. local_ok) go to 900
    read(unit, '(a)', iostat=status) end_marker
    if (status /= 0 .or. trim(end_marker) /= "END_CHECKPOINT") go to 900
    close(unit, iostat=status)
    if (status /= 0) return

    call move_alloc(candidate_root_state, root_state)
    call move_alloc(candidate_root_temperature, root_temperature)
    call move_alloc(candidate_level_one_state, level_one_state)
    call move_alloc(candidate_level_one_temperature, level_one_temperature)
    call move_alloc(candidate_level_two_state, level_two_state)
    call move_alloc(candidate_level_two_temperature, level_two_temperature)
    root_geometry = candidate_root_geometry
    level_one_geometry = candidate_level_one_geometry
    level_two_geometry = candidate_level_two_geometry
    root_patch = candidate_root_patch
    level_one_patch = candidate_level_one_patch
    time = stored_time
    minimum_dt = stored_minimum_dt
    base_density = stored_base_density
    steps = stored_steps
    if (present(regrids)) regrids = stored_regrids
    ok = .true.
    return

900 continue
    close(unit)
  end subroutine read_reactive_eb_amr_three_level_2d_checkpoint

  subroutine simulate_reactive_eb_amr_2d( &
      species, reactions, config, coarse_state, coarse_temperature, &
      coarse_geometry, &
      fine_state, fine_temperature, fine_geometry, patch, fine_active, time, &
      steps, regrids, initial_integrals, final_integrals, minimum_dt, &
      base_density, ok, transport, minimum_transport_theta)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_eb_amr_2d_config), intent(in) :: config
    real(dp), allocatable, intent(out) :: coarse_state(:, :, :)
    real(dp), allocatable, intent(out) :: coarse_temperature(:, :)
    type(eb_geometry_2d), intent(out) :: coarse_geometry
    real(dp), allocatable, intent(out) :: fine_state(:, :, :)
    real(dp), allocatable, intent(out) :: fine_temperature(:, :)
    type(eb_geometry_2d), intent(out) :: fine_geometry
    type(amr_eb_patch_2d), intent(out) :: patch
    logical, intent(out) :: fine_active
    real(dp), intent(out) :: time, minimum_dt, base_density
    integer, intent(out) :: steps, regrids
    real(dp), allocatable, intent(out) :: initial_integrals(:)
    real(dp), allocatable, intent(out) :: final_integrals(:)
    logical, intent(out) :: ok
    type(gas_transport_species), intent(in), optional :: transport(:)
    real(dp), intent(out), optional :: minimum_transport_theta

    real(dp), allocatable :: coarse_candidate(:, :, :)
    real(dp), allocatable :: coarse_candidate_temperature(:, :)
    real(dp), allocatable :: fine_candidate(:, :, :)
    real(dp), allocatable :: fine_candidate_temperature(:, :)
    real(dp) :: coarse_dx, coarse_dy, dt, remaining, time_tolerance
    real(dp) :: coarse_transport_dt, fine_transport_dt
    real(dp) :: maximum_diffusivity, step_transport_theta
    real(dp) :: local_minimum_transport_theta
    logical :: changed, local_ok, stopped_after_checkpoint
    integer :: fine_nx, fine_ny, nvar, last_checkpoint_step
    type(reactive_boundary_set_2d) :: boundaries

    ok = .false.
    time = 0.0_dp
    steps = 0
    regrids = 0
    fine_active = .false.
    minimum_dt = 0.0_dp
    base_density = 0.0_dp
    stopped_after_checkpoint = .false.
    last_checkpoint_step = -1
    local_minimum_transport_theta = 1.0_dp
    if (present(minimum_transport_theta)) minimum_transport_theta = 1.0_dp
    if (.not. supported_reactive_eb_amr_config(config)) return
    if (config%multipatch_enabled .or. config%three_level_enabled) return
    if (config%eb%flow%chemistry_enabled .and. size(reactions) < 1) return
    if (config%eb%flow%transport_enabled .and. .not. present(transport)) return
    call build_reactive_boundary_set_2d( &
      species, config%eb%flow, boundaries, local_ok)
    if (.not. local_ok) return
    if (len_trim(config%restart_file) > 0) then
      call read_reactive_eb_amr_2d_checkpoint( &
        config%restart_file, species, config, coarse_state, &
        coarse_temperature, coarse_geometry, fine_state, fine_temperature, &
        fine_geometry, patch, fine_active, time, steps, regrids, minimum_dt, &
        base_density, local_ok)
      if (.not. local_ok) return
      nvar = size(coarse_state, 1)
    else
      call build_configured_eb_geometry_2d( &
        config%eb, coarse_geometry, local_ok)
      if (.not. local_ok) return
      call build_reactive_eb_amr_patch_geometry_2d( &
        config, coarse_geometry, config%coarse_i_lower, &
        config%coarse_i_upper, config%coarse_j_lower, &
        config%coarse_j_upper, fine_geometry, patch, local_ok)
      if (.not. local_ok) return
      fine_nx = fine_geometry%nx
      fine_ny = fine_geometry%ny

      call initialize_reactive_2d( &
        species, config%eb%flow, coarse_state, coarse_temperature, &
        coarse_dx, coarse_dy, base_density, local_ok)
      if (.not. local_ok) return
      if (abs(coarse_dx - coarse_geometry%dx) > &
          8.0_dp * epsilon(1.0_dp) * coarse_geometry%dx .or. &
          abs(coarse_dy - coarse_geometry%dy) > &
          8.0_dp * epsilon(1.0_dp) * coarse_geometry%dy) return
      nvar = size(coarse_state, 1)
      allocate(fine_state(nvar, fine_nx, fine_ny))
      allocate(fine_temperature(fine_nx, fine_ny))
      call prolong_reactive_eb_patch_pcm_2d( &
        species, coarse_state, coarse_temperature, coarse_geometry, &
        fine_geometry, patch, fine_state, fine_temperature, local_ok)
      if (.not. local_ok) return
      fine_active = .true.
      if (config%dynamic_regridding .and. config%regrid_at_initialization) then
        call regrid_reactive_eb_amr_hierarchy_2d( &
          species, config, coarse_state, coarse_temperature, &
          coarse_geometry, fine_state, fine_temperature, fine_geometry, &
          patch, fine_active, changed, local_ok)
        if (.not. local_ok) return
        if (changed) regrids = regrids + 1
      end if
      minimum_dt = huge(1.0_dp)
    end if

    allocate(initial_integrals(nvar), final_integrals(nvar))
    call compute_reactive_eb_amr_integrals_2d( &
      coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
      fine_active, initial_integrals, local_ok)
    if (.not. local_ok) return
    time_tolerance = 16.0_dp * epsilon(1.0_dp) * &
      max(tiny(1.0_dp), abs(config%eb%flow%final_time))

    do
      remaining = config%eb%flow%final_time - time
      if (remaining <= time_tolerance) exit
      if (steps >= config%eb%flow%maximum_steps) return
      if (allocated(coarse_candidate)) deallocate(coarse_candidate)
      if (allocated(coarse_candidate_temperature)) &
        deallocate(coarse_candidate_temperature)
      if (allocated(fine_candidate)) deallocate(fine_candidate)
      if (allocated(fine_candidate_temperature)) &
        deallocate(fine_candidate_temperature)
      allocate(coarse_candidate, mold=coarse_state)
      allocate(coarse_candidate_temperature, mold=coarse_temperature)
      if (fine_active) then
        allocate(fine_candidate, mold=fine_state)
        allocate(fine_candidate_temperature, mold=fine_temperature)
        call compute_reactive_eb_amr_cfl_timestep_2d( &
          species, coarse_state, coarse_temperature, coarse_geometry, &
          fine_state, fine_temperature, fine_geometry, &
          config%refinement_ratio, config%eb%flow%cfl, dt, local_ok)
      else
        call compute_reactive_eb_cfl_timestep_2d( &
          species, coarse_state, coarse_temperature, coarse_geometry, &
          config%eb%flow%cfl, dt, local_ok)
      end if
      if (.not. local_ok) return
      if (config%eb%flow%transport_enabled) then
        call reactive_eb_transport_timestep_2d( &
          species, transport, coarse_state, coarse_temperature, &
          coarse_geometry, config%eb%flow%transport_cfl, &
          config%eb%flow%viscosity_enabled, &
          config%eb%flow%thermal_conduction_enabled, &
          config%eb%flow%species_diffusion_enabled, coarse_transport_dt, &
          maximum_diffusivity, local_ok)
        if (.not. local_ok) return
        dt = min(dt, coarse_transport_dt)
        if (fine_active) then
          call reactive_eb_transport_timestep_2d( &
            species, transport, fine_state, fine_temperature, fine_geometry, &
            config%eb%flow%transport_cfl, &
            config%eb%flow%viscosity_enabled, &
            config%eb%flow%thermal_conduction_enabled, &
            config%eb%flow%species_diffusion_enabled, fine_transport_dt, &
            maximum_diffusivity, local_ok)
          if (.not. local_ok) return
          dt = min(dt, real(config%refinement_ratio, dp) * fine_transport_dt)
        end if
      end if
      dt = min(dt, remaining)
      if (fine_active) then
        if (config%eb%flow%transport_enabled) then
          call advance_two_level_reactive_eb_strang_2d( &
            species, reactions, coarse_state, coarse_temperature, &
            coarse_geometry, fine_state, fine_temperature, fine_geometry, &
            patch, config%eb%flow%riemann_solver, &
            config%eb%flow%reconstruction, config%eb%flow%limiter, &
            config%eb%state_redist_max_order, dt, &
            config%eb%flow%chemistry_enabled, &
            config%eb%flow%chemistry_relative_tolerance, &
            config%eb%flow%chemistry_absolute_tolerance, coarse_candidate, &
            coarse_candidate_temperature, fine_candidate, &
            fine_candidate_temperature, local_ok, &
            config%eb%state_redist_target_volume_fraction, transport, &
            config%eb%flow%transport_enabled, &
            config%eb%flow%viscosity_enabled, &
            config%eb%flow%thermal_conduction_enabled, &
            config%eb%flow%species_diffusion_enabled, &
            config%eb%flow%barodiffusion_enabled, step_transport_theta, &
            boundaries)
        else
          call advance_two_level_reactive_eb_strang_2d( &
            species, reactions, coarse_state, coarse_temperature, &
            coarse_geometry, fine_state, fine_temperature, fine_geometry, &
            patch, config%eb%flow%riemann_solver, &
            config%eb%flow%reconstruction, config%eb%flow%limiter, &
            config%eb%state_redist_max_order, dt, &
            config%eb%flow%chemistry_enabled, &
            config%eb%flow%chemistry_relative_tolerance, &
            config%eb%flow%chemistry_absolute_tolerance, coarse_candidate, &
            coarse_candidate_temperature, fine_candidate, &
            fine_candidate_temperature, local_ok, &
            config%eb%state_redist_target_volume_fraction)
          step_transport_theta = 1.0_dp
        end if
      else
        if (config%eb%flow%transport_enabled) then
          call advance_reactive_eb_strang_2d( &
            species, reactions, coarse_state, coarse_temperature, &
            coarse_geometry, config%eb%flow%riemann_solver, dt, &
            config%eb%flow%chemistry_enabled, &
            config%eb%flow%chemistry_relative_tolerance, &
            config%eb%flow%chemistry_absolute_tolerance, coarse_candidate, &
            coarse_candidate_temperature, local_ok, &
            config%eb%state_redist_target_volume_fraction, &
            config%eb%flow%reconstruction, config%eb%flow%limiter, &
            config%eb%state_redist_max_order, transport, &
            config%eb%flow%transport_enabled, &
            config%eb%flow%viscosity_enabled, &
            config%eb%flow%thermal_conduction_enabled, &
            config%eb%flow%species_diffusion_enabled, &
            config%eb%flow%barodiffusion_enabled, step_transport_theta, &
            boundaries)
        else
          call advance_reactive_eb_strang_2d( &
            species, reactions, coarse_state, coarse_temperature, &
            coarse_geometry, config%eb%flow%riemann_solver, dt, &
            config%eb%flow%chemistry_enabled, &
            config%eb%flow%chemistry_relative_tolerance, &
            config%eb%flow%chemistry_absolute_tolerance, coarse_candidate, &
            coarse_candidate_temperature, local_ok, &
            config%eb%state_redist_target_volume_fraction, &
            config%eb%flow%reconstruction, config%eb%flow%limiter, &
            config%eb%state_redist_max_order)
          step_transport_theta = 1.0_dp
        end if
      end if
      if (.not. local_ok) return
      coarse_state = coarse_candidate
      coarse_temperature = coarse_candidate_temperature
      if (fine_active) then
        fine_state = fine_candidate
        fine_temperature = fine_candidate_temperature
      end if
      time = time + dt
      minimum_dt = min(minimum_dt, dt)
      local_minimum_transport_theta = min( &
        local_minimum_transport_theta, step_transport_theta)
      steps = steps + 1
      if (config%dynamic_regridding .and. &
          modulo(steps, config%regrid_interval) == 0) then
        call regrid_reactive_eb_amr_hierarchy_2d( &
          species, config, coarse_state, coarse_temperature, &
          coarse_geometry, fine_state, fine_temperature, fine_geometry, &
          patch, fine_active, changed, local_ok)
        if (.not. local_ok) return
        if (changed) regrids = regrids + 1
      end if
      if (config%checkpoint_interval > 0) then
        if (modulo(steps, config%checkpoint_interval) == 0) then
          call write_reactive_eb_amr_2d_checkpoint( &
            config%checkpoint_file, species, config, coarse_state, &
            coarse_temperature, coarse_geometry, fine_state, &
            fine_temperature, fine_geometry, patch, fine_active, time, steps, &
            regrids, minimum_dt, base_density, local_ok)
          if (.not. local_ok) return
          last_checkpoint_step = steps
          if (config%checkpoint_stop_after_write) then
            stopped_after_checkpoint = .true.
            exit
          end if
        end if
      end if
    end do
    if (.not. stopped_after_checkpoint) &
      time = config%eb%flow%final_time
    if (len_trim(config%checkpoint_file) > 0 .and. &
        last_checkpoint_step /= steps) then
      call write_reactive_eb_amr_2d_checkpoint( &
        config%checkpoint_file, species, config, coarse_state, &
        coarse_temperature, coarse_geometry, fine_state, fine_temperature, &
        fine_geometry, patch, fine_active, time, steps, regrids, minimum_dt, &
        base_density, local_ok)
      if (.not. local_ok) return
    end if
    call compute_reactive_eb_amr_integrals_2d( &
      coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
      fine_active, final_integrals, local_ok)
    if (.not. local_ok) return
    if (present(minimum_transport_theta)) &
      minimum_transport_theta = local_minimum_transport_theta
    ok = steps > 0 .and. ieee_is_finite(minimum_dt) .and. minimum_dt > 0.0_dp
  end subroutine simulate_reactive_eb_amr_2d

  subroutine simulate_three_level_reactive_eb_amr_2d( &
      species, reactions, config, root_state, root_temperature, &
      root_geometry, level_one_state, level_one_temperature, &
      level_one_geometry, root_patch, level_two_state, &
      level_two_temperature, level_two_geometry, level_one_patch, time, &
      steps, regrids, initial_integrals, final_integrals, minimum_dt, &
      base_density, ok, failure_context, transport, minimum_transport_theta)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_eb_amr_2d_config), intent(in) :: config
    real(dp), allocatable, intent(out) :: root_state(:, :, :)
    real(dp), allocatable, intent(out) :: root_temperature(:, :)
    type(eb_geometry_2d), intent(out) :: root_geometry
    real(dp), allocatable, intent(out) :: level_one_state(:, :, :)
    real(dp), allocatable, intent(out) :: level_one_temperature(:, :)
    type(eb_geometry_2d), intent(out) :: level_one_geometry
    type(amr_eb_patch_2d), intent(out) :: root_patch
    real(dp), allocatable, intent(out) :: level_two_state(:, :, :)
    real(dp), allocatable, intent(out) :: level_two_temperature(:, :)
    type(eb_geometry_2d), intent(out) :: level_two_geometry
    type(amr_eb_patch_2d), intent(out) :: level_one_patch
    real(dp), intent(out) :: time, minimum_dt, base_density
    integer, intent(out) :: steps, regrids
    real(dp), allocatable, intent(out) :: initial_integrals(:)
    real(dp), allocatable, intent(out) :: final_integrals(:)
    logical, intent(out) :: ok
    character(len=*), intent(out), optional :: failure_context
    type(gas_transport_species), intent(in), optional :: transport(:)
    real(dp), intent(out), optional :: minimum_transport_theta

    real(dp), allocatable :: root_candidate(:, :, :)
    real(dp), allocatable :: root_candidate_temperature(:, :)
    real(dp), allocatable :: level_one_candidate(:, :, :)
    real(dp), allocatable :: level_one_candidate_temperature(:, :)
    real(dp), allocatable :: level_two_candidate(:, :, :)
    real(dp), allocatable :: level_two_candidate_temperature(:, :)
    real(dp) :: root_dx, root_dy, dt, remaining, time_tolerance
    real(dp) :: root_transport_dt, level_one_transport_dt
    real(dp) :: level_two_transport_dt, maximum_diffusivity
    real(dp) :: step_transport_theta, local_minimum_transport_theta
    logical :: changed, local_ok, stopped_after_checkpoint
    integer :: nvar, last_checkpoint_step
    type(reactive_boundary_set_2d) :: boundaries

    root_geometry = eb_geometry_2d()
    level_one_geometry = eb_geometry_2d()
    level_two_geometry = eb_geometry_2d()
    root_patch = amr_eb_patch_2d()
    level_one_patch = amr_eb_patch_2d()
    time = 0.0_dp
    steps = 0
    regrids = 0
    minimum_dt = 0.0_dp
    base_density = 0.0_dp
    stopped_after_checkpoint = .false.
    last_checkpoint_step = -1
    local_minimum_transport_theta = 1.0_dp
    if (present(minimum_transport_theta)) minimum_transport_theta = 1.0_dp
    ok = .false.
    if (present(failure_context)) failure_context = "validation"
    if (.not. supported_three_level_reactive_eb_amr_config(config)) return
    if (config%eb%flow%chemistry_enabled .and. size(reactions) < 1) return
    if (config%eb%flow%transport_enabled .and. .not. present(transport)) return
    call build_reactive_boundary_set_2d( &
      species, config%eb%flow, boundaries, local_ok)
    if (.not. local_ok) return

    if (len_trim(config%restart_file) > 0) then
      if (present(failure_context)) failure_context = "restart"
      call read_reactive_eb_amr_three_level_2d_checkpoint( &
        config%restart_file, species, config, root_state, root_temperature, &
        root_geometry, level_one_state, level_one_temperature, &
        level_one_geometry, root_patch, level_two_state, &
        level_two_temperature, level_two_geometry, level_one_patch, time, &
        steps, minimum_dt, base_density, local_ok, regrids)
      if (.not. local_ok) return
      nvar = size(root_state, 1)
    else
      if (present(failure_context)) failure_context = "root geometry"
      call build_configured_eb_geometry_2d(config%eb, root_geometry, local_ok)
    if (.not. local_ok) return
    if (present(failure_context)) failure_context = "middle geometry"
    call build_reactive_eb_amr_patch_geometry_2d( &
      config, root_geometry, config%coarse_i_lower, &
      config%coarse_i_upper, config%coarse_j_lower, &
      config%coarse_j_upper, level_one_geometry, root_patch, local_ok)
    if (.not. local_ok) return
    if (present(failure_context)) failure_context = "finest geometry"
    call build_reactive_eb_amr_patch_geometry_2d( &
      config, level_one_geometry, config%level_two_i_lower, &
      config%level_two_i_upper, config%level_two_j_lower, &
      config%level_two_j_upper, level_two_geometry, level_one_patch, local_ok)
    if (.not. local_ok) return

    if (present(failure_context)) failure_context = "root initialization"
    call initialize_reactive_2d( &
      species, config%eb%flow, root_state, root_temperature, root_dx, &
      root_dy, base_density, local_ok)
    if (.not. local_ok) return
    if (abs(root_dx - root_geometry%dx) > &
        8.0_dp * epsilon(1.0_dp) * root_geometry%dx .or. &
        abs(root_dy - root_geometry%dy) > &
        8.0_dp * epsilon(1.0_dp) * root_geometry%dy) return
    nvar = size(root_state, 1)
    allocate(level_one_state( &
      nvar, level_one_geometry%nx, level_one_geometry%ny))
    allocate(level_one_temperature( &
      level_one_geometry%nx, level_one_geometry%ny))
    if (present(failure_context)) failure_context = "middle prolongation"
    call prolong_reactive_eb_patch_pcm_2d( &
      species, root_state, root_temperature, root_geometry, &
      level_one_geometry, root_patch, level_one_state, &
      level_one_temperature, local_ok)
    if (.not. local_ok) return
    allocate(level_two_state( &
      nvar, level_two_geometry%nx, level_two_geometry%ny))
    allocate(level_two_temperature( &
      level_two_geometry%nx, level_two_geometry%ny))
    if (present(failure_context)) failure_context = "finest prolongation"
    call prolong_reactive_eb_patch_pcm_2d( &
      species, level_one_state, level_one_temperature, level_one_geometry, &
      level_two_geometry, level_one_patch, level_two_state, &
      level_two_temperature, local_ok)
    if (.not. local_ok) return
      if (config%dynamic_regridding .and. &
          config%regrid_at_initialization) then
        if (present(failure_context)) failure_context = "initial finest regrid"
        call regrid_three_level_reactive_eb_amr_finest_2d( &
          species, config, level_one_state, level_one_temperature, &
          level_one_geometry, level_two_state, level_two_temperature, &
          level_two_geometry, level_one_patch, changed, local_ok)
        if (.not. local_ok) return
        if (changed) regrids = regrids + 1
      end if
      minimum_dt = huge(1.0_dp)
    end if

    allocate(initial_integrals(nvar), final_integrals(nvar))
    if (present(failure_context)) failure_context = "initial composite integral"
    call composite_three_level_eb_integral_2d( &
      root_state, root_geometry, level_one_state, level_one_geometry, &
      root_patch, level_two_state, level_two_geometry, level_one_patch, &
      initial_integrals, local_ok)
    if (.not. local_ok) return
    time_tolerance = 16.0_dp * epsilon(1.0_dp) * &
      max(tiny(1.0_dp), abs(config%eb%flow%final_time))

    do
      remaining = config%eb%flow%final_time - time
      if (remaining <= time_tolerance) exit
      if (steps >= config%eb%flow%maximum_steps) return
      if (present(failure_context)) failure_context = "three-level CFL"
      call compute_three_level_reactive_eb_cfl_timestep_2d( &
        species, root_state, root_temperature, root_geometry, &
        level_one_state, level_one_temperature, level_one_geometry, &
        root_patch, level_two_state, level_two_temperature, &
        level_two_geometry, level_one_patch, config%eb%flow%cfl, dt, local_ok)
      if (.not. local_ok) return
      if (config%eb%flow%transport_enabled) then
        if (present(failure_context)) failure_context = "three-level transport CFL"
        call reactive_eb_transport_timestep_2d( &
          species, transport, root_state, root_temperature, root_geometry, &
          config%eb%flow%transport_cfl, config%eb%flow%viscosity_enabled, &
          config%eb%flow%thermal_conduction_enabled, &
          config%eb%flow%species_diffusion_enabled, root_transport_dt, &
          maximum_diffusivity, local_ok)
        if (.not. local_ok) return
        call reactive_eb_transport_timestep_2d( &
          species, transport, level_one_state, level_one_temperature, &
          level_one_geometry, config%eb%flow%transport_cfl, &
          config%eb%flow%viscosity_enabled, &
          config%eb%flow%thermal_conduction_enabled, &
          config%eb%flow%species_diffusion_enabled, level_one_transport_dt, &
          maximum_diffusivity, local_ok)
        if (.not. local_ok) return
        call reactive_eb_transport_timestep_2d( &
          species, transport, level_two_state, level_two_temperature, &
          level_two_geometry, config%eb%flow%transport_cfl, &
          config%eb%flow%viscosity_enabled, &
          config%eb%flow%thermal_conduction_enabled, &
          config%eb%flow%species_diffusion_enabled, level_two_transport_dt, &
          maximum_diffusivity, local_ok)
        if (.not. local_ok) return
        dt = min(dt, root_transport_dt)
        dt = min(dt, real(root_patch%refinement_ratio, dp) * &
          level_one_transport_dt)
        dt = min(dt, real(root_patch%refinement_ratio * &
          level_one_patch%refinement_ratio, dp) * level_two_transport_dt)
      end if
      dt = min(dt, remaining)
      if (allocated(root_candidate)) deallocate(root_candidate)
      if (allocated(root_candidate_temperature)) &
        deallocate(root_candidate_temperature)
      if (allocated(level_one_candidate)) deallocate(level_one_candidate)
      if (allocated(level_one_candidate_temperature)) &
        deallocate(level_one_candidate_temperature)
      if (allocated(level_two_candidate)) deallocate(level_two_candidate)
      if (allocated(level_two_candidate_temperature)) &
        deallocate(level_two_candidate_temperature)
      allocate(root_candidate, mold=root_state)
      allocate(root_candidate_temperature, mold=root_temperature)
      allocate(level_one_candidate, mold=level_one_state)
      allocate(level_one_candidate_temperature, mold=level_one_temperature)
      allocate(level_two_candidate, mold=level_two_state)
      allocate(level_two_candidate_temperature, mold=level_two_temperature)
      if (present(failure_context)) failure_context = "three-level advance"
      if (config%eb%flow%transport_enabled) then
        call advance_three_level_reactive_eb_strang_2d( &
          species, reactions, root_state, root_temperature, root_geometry, &
          level_one_state, level_one_temperature, level_one_geometry, &
          root_patch, level_two_state, level_two_temperature, &
          level_two_geometry, level_one_patch, &
          config%eb%flow%riemann_solver, config%eb%flow%reconstruction, &
          config%eb%flow%limiter, config%eb%state_redist_max_order, dt, &
          config%eb%flow%chemistry_enabled, &
          config%eb%flow%chemistry_relative_tolerance, &
          config%eb%flow%chemistry_absolute_tolerance, root_candidate, &
          root_candidate_temperature, level_one_candidate, &
          level_one_candidate_temperature, level_two_candidate, &
          level_two_candidate_temperature, local_ok, &
          config%eb%state_redist_target_volume_fraction, transport, &
          config%eb%flow%transport_enabled, &
          config%eb%flow%viscosity_enabled, &
          config%eb%flow%thermal_conduction_enabled, &
          config%eb%flow%species_diffusion_enabled, &
          config%eb%flow%barodiffusion_enabled, step_transport_theta, &
          boundaries)
      else
        call advance_three_level_reactive_eb_strang_2d( &
          species, reactions, root_state, root_temperature, root_geometry, &
          level_one_state, level_one_temperature, level_one_geometry, &
          root_patch, level_two_state, level_two_temperature, &
          level_two_geometry, level_one_patch, &
          config%eb%flow%riemann_solver, config%eb%flow%reconstruction, &
          config%eb%flow%limiter, config%eb%state_redist_max_order, dt, &
          config%eb%flow%chemistry_enabled, &
          config%eb%flow%chemistry_relative_tolerance, &
          config%eb%flow%chemistry_absolute_tolerance, root_candidate, &
          root_candidate_temperature, level_one_candidate, &
          level_one_candidate_temperature, level_two_candidate, &
          level_two_candidate_temperature, local_ok, &
          config%eb%state_redist_target_volume_fraction)
        step_transport_theta = 1.0_dp
      end if
      if (.not. local_ok) return
      root_state = root_candidate
      root_temperature = root_candidate_temperature
      level_one_state = level_one_candidate
      level_one_temperature = level_one_candidate_temperature
      level_two_state = level_two_candidate
      level_two_temperature = level_two_candidate_temperature
      time = time + dt
      minimum_dt = min(minimum_dt, dt)
      local_minimum_transport_theta = min( &
        local_minimum_transport_theta, step_transport_theta)
      steps = steps + 1
      if (config%dynamic_regridding .and. &
          modulo(steps, config%regrid_interval) == 0) then
        if (present(failure_context)) failure_context = "scheduled finest regrid"
        call regrid_three_level_reactive_eb_amr_finest_2d( &
          species, config, level_one_state, level_one_temperature, &
          level_one_geometry, level_two_state, level_two_temperature, &
          level_two_geometry, level_one_patch, changed, local_ok)
        if (.not. local_ok) return
        if (changed) regrids = regrids + 1
      end if
      if (config%checkpoint_interval > 0) then
        if (modulo(steps, config%checkpoint_interval) == 0) then
          if (present(failure_context)) failure_context = "scheduled checkpoint"
          call write_reactive_eb_amr_three_level_2d_checkpoint( &
            config%checkpoint_file, species, config, root_state, &
            root_temperature, root_geometry, level_one_state, &
            level_one_temperature, level_one_geometry, root_patch, &
            level_two_state, level_two_temperature, level_two_geometry, &
            level_one_patch, time, steps, minimum_dt, base_density, local_ok, &
            regrids)
          if (.not. local_ok) return
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
      if (present(failure_context)) failure_context = "final checkpoint"
      call write_reactive_eb_amr_three_level_2d_checkpoint( &
        config%checkpoint_file, species, config, root_state, &
        root_temperature, root_geometry, level_one_state, &
        level_one_temperature, level_one_geometry, root_patch, &
        level_two_state, level_two_temperature, level_two_geometry, &
        level_one_patch, time, steps, minimum_dt, base_density, local_ok, &
        regrids)
      if (.not. local_ok) return
    end if
    if (present(failure_context)) failure_context = "final composite integral"
    call composite_three_level_eb_integral_2d( &
      root_state, root_geometry, level_one_state, level_one_geometry, &
      root_patch, level_two_state, level_two_geometry, level_one_patch, &
      final_integrals, local_ok)
    if (.not. local_ok) return
    ok = steps > 0 .and. ieee_is_finite(minimum_dt) .and. minimum_dt > 0.0_dp
    if (present(minimum_transport_theta)) &
      minimum_transport_theta = local_minimum_transport_theta
    if (ok .and. present(failure_context)) failure_context = "none"
  end subroutine simulate_three_level_reactive_eb_amr_2d

  subroutine simulate_reactive_eb_amr_patch_set_2d( &
      species, reactions, config, coarse_state, coarse_temperature, &
      coarse_geometry, patch_set, time, steps, regrids, initial_integrals, &
      final_integrals, minimum_dt, base_density, ok, failure_context, &
      transport, minimum_transport_theta)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_eb_amr_2d_config), intent(in) :: config
    real(dp), allocatable, intent(out) :: coarse_state(:, :, :)
    real(dp), allocatable, intent(out) :: coarse_temperature(:, :)
    type(eb_geometry_2d), intent(out) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(out) :: patch_set
    real(dp), intent(out) :: time, minimum_dt, base_density
    integer, intent(out) :: steps, regrids
    real(dp), allocatable, intent(out) :: initial_integrals(:)
    real(dp), allocatable, intent(out) :: final_integrals(:)
    logical, intent(out) :: ok
    character(len=*), intent(out), optional :: failure_context
    type(gas_transport_species), intent(in), optional :: transport(:)
    real(dp), intent(out), optional :: minimum_transport_theta

    type(amr_eb_regrid_plan_collection_2d) :: initial_collection
    type(reactive_eb_patch_set_2d) :: candidate_set
    type(reactive_boundary_set_2d) :: boundaries
    type(eb_geometry_2d), allocatable :: fine_geometries(:)
    real(dp), allocatable :: candidate_state(:, :, :)
    real(dp), allocatable :: candidate_temperature(:, :)
    real(dp) :: coarse_dx, coarse_dy, dt, remaining, time_tolerance
    real(dp) :: coarse_transport_dt, fine_transport_dt, maximum_diffusivity
    real(dp) :: step_transport_theta, local_minimum_transport_theta
    logical :: changed, local_ok, stopped_after_checkpoint
    integer :: child, nvar, last_checkpoint_step

    ok = .false.
    time = 0.0_dp
    steps = 0
    regrids = 0
    minimum_dt = 0.0_dp
    base_density = 0.0_dp
    stopped_after_checkpoint = .false.
    last_checkpoint_step = -1
    local_minimum_transport_theta = 1.0_dp
    if (present(minimum_transport_theta)) minimum_transport_theta = 1.0_dp
    patch_set = reactive_eb_patch_set_2d()
    if (present(failure_context)) failure_context = "input validation"
    if (.not. supported_reactive_eb_amr_config(config) .or. &
        .not. config%multipatch_enabled .or. &
        config%three_level_enabled) return
    if (config%eb%flow%chemistry_enabled .and. size(reactions) < 1) return
    if (config%eb%flow%transport_enabled .and. .not. present(transport)) return
    call build_reactive_boundary_set_2d( &
      species, config%eb%flow, boundaries, local_ok)
    if (.not. local_ok) return
    if (len_trim(config%restart_file) > 0) then
      if (present(failure_context)) failure_context = "checkpoint restart"
      call read_reactive_eb_amr_patch_set_2d_checkpoint( &
        config%restart_file, species, config, coarse_state, &
        coarse_temperature, coarse_geometry, patch_set, time, steps, &
        regrids, minimum_dt, base_density, local_ok)
      if (.not. local_ok) return
      nvar = size(coarse_state, 1)
    else
      if (present(failure_context)) failure_context = "coarse geometry"
      call build_configured_eb_geometry_2d( &
        config%eb, coarse_geometry, local_ok)
      if (.not. local_ok) return
      if (present(failure_context)) failure_context = "coarse initialization"
      call initialize_reactive_2d( &
        species, config%eb%flow, coarse_state, coarse_temperature, &
        coarse_dx, coarse_dy, base_density, local_ok)
      if (.not. local_ok) return
      if (abs(coarse_dx - coarse_geometry%dx) > &
          8.0_dp * epsilon(1.0_dp) * coarse_geometry%dx .or. &
          abs(coarse_dy - coarse_geometry%dy) > &
          8.0_dp * epsilon(1.0_dp) * coarse_geometry%dy) return
      nvar = size(coarse_state, 1)
      if (present(failure_context)) failure_context = "initial patch plan"
      call build_initial_reactive_eb_patch_collection_2d( &
        config, initial_collection)
      if (.not. initial_collection%is_valid()) return
      if (present(failure_context)) failure_context = &
        "initial fine geometries"
      call build_reactive_eb_patch_set_geometries_2d( &
        config, coarse_geometry, initial_collection, fine_geometries, local_ok)
      if (.not. local_ok) return
      if (present(failure_context)) failure_context = "initial patch set"
      call initialize_reactive_eb_patch_set_2d( &
        species, coarse_state, coarse_temperature, coarse_geometry, &
        fine_geometries, initial_collection, config%refinement_ratio, &
        patch_set, local_ok)
      if (.not. local_ok) return
      if (config%regrid_at_initialization) then
        if (present(failure_context)) failure_context = "initial regrid"
        call regrid_reactive_eb_amr_patch_set_2d( &
          species, config, coarse_state, coarse_temperature, coarse_geometry, &
          patch_set, changed, local_ok)
        if (.not. local_ok) return
        if (changed) regrids = regrids + 1
      end if
      minimum_dt = huge(1.0_dp)
    end if
    allocate(initial_integrals(nvar), final_integrals(nvar))
    if (present(failure_context)) failure_context = "initial integral"
    call composite_reactive_eb_patch_set_integral_2d( &
      coarse_state, coarse_geometry, patch_set, initial_integrals, local_ok)
    if (.not. local_ok) return
    time_tolerance = 16.0_dp * epsilon(1.0_dp) * &
      max(tiny(1.0_dp), abs(config%eb%flow%final_time))

    do
      remaining = config%eb%flow%final_time - time
      if (remaining <= time_tolerance) exit
      if (steps >= config%eb%flow%maximum_steps) return
      if (present(failure_context)) failure_context = "CFL selection"
      call compute_reactive_eb_patch_set_cfl_timestep_2d( &
        species, coarse_state, coarse_temperature, coarse_geometry, &
        patch_set, config%eb%flow%cfl, dt, local_ok)
      if (.not. local_ok) return
      if (config%eb%flow%transport_enabled) then
        if (present(failure_context)) &
          failure_context = "multipatch transport CFL"
        call reactive_eb_transport_timestep_2d( &
          species, transport, coarse_state, coarse_temperature, &
          coarse_geometry, config%eb%flow%transport_cfl, &
          config%eb%flow%viscosity_enabled, &
          config%eb%flow%thermal_conduction_enabled, &
          config%eb%flow%species_diffusion_enabled, coarse_transport_dt, &
          maximum_diffusivity, local_ok)
        if (.not. local_ok) return
        dt = min(dt, coarse_transport_dt)
        do child = 1, patch_set%patch_count()
          call reactive_eb_transport_timestep_2d( &
            species, transport, patch_set%children(child)%state, &
            patch_set%children(child)%temperature, &
            patch_set%children(child)%geometry, &
            config%eb%flow%transport_cfl, &
            config%eb%flow%viscosity_enabled, &
            config%eb%flow%thermal_conduction_enabled, &
            config%eb%flow%species_diffusion_enabled, fine_transport_dt, &
            maximum_diffusivity, local_ok)
          if (.not. local_ok) return
          dt = min(dt, real( &
            patch_set%children(child)%patch%refinement_ratio, dp) * &
            fine_transport_dt)
        end do
      end if
      dt = min(dt, remaining)
      if (allocated(candidate_state)) deallocate(candidate_state)
      if (allocated(candidate_temperature)) deallocate(candidate_temperature)
      allocate(candidate_state, mold=coarse_state)
      allocate(candidate_temperature, mold=coarse_temperature)
      if (present(failure_context)) failure_context = "patch-set advance"
      if (config%eb%flow%transport_enabled) then
        call advance_reactive_eb_patch_set_strang_2d( &
          species, reactions, coarse_state, coarse_temperature, &
          coarse_geometry, patch_set, config%eb%flow%riemann_solver, &
          config%eb%flow%reconstruction, config%eb%flow%limiter, &
          config%eb%state_redist_max_order, dt, &
          config%eb%flow%chemistry_enabled, &
          config%eb%flow%chemistry_relative_tolerance, &
          config%eb%flow%chemistry_absolute_tolerance, candidate_state, &
          candidate_temperature, candidate_set, local_ok, &
          config%eb%state_redist_target_volume_fraction, failure_context, &
          transport, config%eb%flow%transport_enabled, &
          config%eb%flow%viscosity_enabled, &
          config%eb%flow%thermal_conduction_enabled, &
          config%eb%flow%species_diffusion_enabled, &
          config%eb%flow%barodiffusion_enabled, step_transport_theta, &
          boundaries)
      else
        call advance_reactive_eb_patch_set_strang_2d( &
          species, reactions, coarse_state, coarse_temperature, &
          coarse_geometry, patch_set, config%eb%flow%riemann_solver, &
          config%eb%flow%reconstruction, config%eb%flow%limiter, &
          config%eb%state_redist_max_order, dt, &
          config%eb%flow%chemistry_enabled, &
          config%eb%flow%chemistry_relative_tolerance, &
          config%eb%flow%chemistry_absolute_tolerance, candidate_state, &
          candidate_temperature, candidate_set, local_ok, &
          config%eb%state_redist_target_volume_fraction, failure_context)
        step_transport_theta = 1.0_dp
      end if
      if (.not. local_ok) return
      coarse_state = candidate_state
      coarse_temperature = candidate_temperature
      patch_set = candidate_set
      time = time + dt
      minimum_dt = min(minimum_dt, dt)
      local_minimum_transport_theta = min( &
        local_minimum_transport_theta, step_transport_theta)
      steps = steps + 1
      if (modulo(steps, config%regrid_interval) == 0) then
        if (present(failure_context)) failure_context = "periodic regrid"
        call regrid_reactive_eb_amr_patch_set_2d( &
          species, config, coarse_state, coarse_temperature, &
          coarse_geometry, patch_set, changed, local_ok)
        if (.not. local_ok) return
        if (changed) regrids = regrids + 1
      end if
      if (config%checkpoint_interval > 0) then
        if (modulo(steps, config%checkpoint_interval) == 0) then
          if (present(failure_context)) failure_context = "checkpoint write"
          call write_reactive_eb_amr_patch_set_2d_checkpoint( &
            config%checkpoint_file, species, config, coarse_state, &
            coarse_temperature, coarse_geometry, patch_set, time, steps, &
            regrids, minimum_dt, base_density, local_ok)
          if (.not. local_ok) return
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
      if (present(failure_context)) failure_context = "final checkpoint write"
      call write_reactive_eb_amr_patch_set_2d_checkpoint( &
        config%checkpoint_file, species, config, coarse_state, &
        coarse_temperature, coarse_geometry, patch_set, time, steps, regrids, &
        minimum_dt, base_density, local_ok)
      if (.not. local_ok) return
    end if
    if (present(failure_context)) failure_context = "final integral"
    call composite_reactive_eb_patch_set_integral_2d( &
      coarse_state, coarse_geometry, patch_set, final_integrals, local_ok)
    if (.not. local_ok) return
    ok = steps > 0 .and. patch_set%is_valid(coarse_geometry, nvar) .and. &
      ieee_is_finite(minimum_dt) .and. minimum_dt > 0.0_dp
    if (present(minimum_transport_theta)) &
      minimum_transport_theta = local_minimum_transport_theta
    if (ok .and. present(failure_context)) failure_context = "none"
  end subroutine simulate_reactive_eb_amr_patch_set_2d

  subroutine simulate_reactive_amr_eb_patch_tree_2d( &
      species, reactions, transport, config, solution, time, steps, regrids, &
      initial_integrals, final_integrals, minimum_dt, base_density, ok, &
      failure_context, minimum_transport_theta)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_eb_amr_2d_config), intent(in) :: config
    type(reactive_amr_eb_patch_tree_2d), intent(out) :: solution
    real(dp), intent(out) :: time, minimum_dt, base_density
    integer, intent(out) :: steps, regrids
    real(dp), allocatable, intent(out) :: initial_integrals(:)
    real(dp), allocatable, intent(out) :: final_integrals(:)
    logical, intent(out) :: ok
    character(len=*), intent(out), optional :: failure_context
    real(dp), intent(out), optional :: minimum_transport_theta

    type(amr_eb_patch_tree_level_plan_2d), allocatable :: empty_plans(:)
    type(amr_eb_patch_tree_topology_2d) :: topology
    type(amr_eb_tagging_criteria_2d) :: criteria
    type(reactive_boundary_set_2d) :: boundaries
    type(eb_geometry_2d) :: root_geometry
    real(dp), allocatable :: root_state(:, :, :), root_temperature(:, :)
    real(dp) :: dx, dy, dt, remaining, step_theta, time_tolerance
    real(dp) :: local_minimum_transport_theta
    logical :: changed, local_ok, stopped_after_checkpoint
    integer :: last_checkpoint_step, nvar, tagged_cells

    solution = reactive_amr_eb_patch_tree_2d()
    time = 0.0_dp
    minimum_dt = 0.0_dp
    base_density = 0.0_dp
    steps = 0
    regrids = 0
    ok = .false.
    local_minimum_transport_theta = 1.0_dp
    stopped_after_checkpoint = .false.
    last_checkpoint_step = -1
    if (present(failure_context)) failure_context = "input validation"
    if (present(minimum_transport_theta)) minimum_transport_theta = 1.0_dp
    if (.not. supported_reactive_eb_amr_config(config) .or. &
        config%three_level_enabled .or. config%multipatch_enabled .or. &
        config%patch_tree_maximum_levels < 1 .or. &
        config%patch_tree_maximum_levels > 64 .or. &
        (config%eb%flow%chemistry_enabled .and. size(reactions) < 1) .or. &
        (config%eb%flow%transport_enabled .and. &
         size(transport) /= size(species))) return
    call build_reactive_boundary_set_2d( &
      species, config%eb%flow, boundaries, local_ok)
    if (.not. local_ok) return

    criteria%relative_gradient_threshold = &
      config%regrid_relative_temperature_gradient
    criteria%absolute_gradient_threshold = &
      config%regrid_absolute_temperature_gradient
    criteria%scale_floor = config%regrid_temperature_scale_floor
    criteria%buffer_cells = config%regrid_buffer_cells
    criteria%minimum_patch_cells_x = config%regrid_minimum_patch_cells_x
    criteria%minimum_patch_cells_y = config%regrid_minimum_patch_cells_y
    criteria%maximum_patch_gap_cells = config%regrid_maximum_patch_gap_cells

    if (len_trim(config%restart_file) > 0) then
      if (present(failure_context)) failure_context = "checkpoint restart"
      call read_reactive_amr_eb_patch_tree_2d_checkpoint( &
        config%restart_file, species, config%patch_tree_maximum_levels, &
        solution, time, steps, regrids, minimum_dt, local_ok)
      if (.not. local_ok) return
      nvar = solution%nvar
      base_density = solution%levels(1)%patches(1)%state(irho, 1, 1)
    else
      if (present(failure_context)) failure_context = "root geometry"
      call build_configured_eb_geometry_2d(config%eb, root_geometry, local_ok)
      if (.not. local_ok) return
      if (present(failure_context)) failure_context = "root initialization"
      call initialize_reactive_2d( &
        species, config%eb%flow, root_state, root_temperature, dx, dy, &
        base_density, local_ok)
      if (.not. local_ok) return
      if (abs(dx - root_geometry%dx) > &
            8.0_dp * epsilon(1.0_dp) * root_geometry%dx .or. &
          abs(dy - root_geometry%dy) > &
            8.0_dp * epsilon(1.0_dp) * root_geometry%dy) return
      allocate(empty_plans(0))
      if (present(failure_context)) failure_context = "root topology"
      call initialize_amr_eb_patch_tree_topology_2d( &
        root_geometry, empty_plans, topology, local_ok)
      if (.not. local_ok) return
      if (present(failure_context)) failure_context = "root tree"
      call initialize_reactive_amr_eb_patch_tree_2d( &
        species, root_state, root_temperature, topology, solution, local_ok)
      if (.not. local_ok) return
      nvar = solution%nvar
      if (config%dynamic_regridding .and. &
          config%regrid_at_initialization) then
        if (present(failure_context)) failure_context = "initial regrid"
        call regrid_tagged_reactive_amr_eb_patch_tree_2d( &
          species, solution, criteria, config%patch_tree_maximum_levels, &
          config%refinement_ratio, build_patch_tree_geometry, local_ok, &
          changed, tagged_cells)
        if (.not. local_ok) return
        if (changed) regrids = regrids + 1
      end if
      minimum_dt = huge(1.0_dp)
    end if

    allocate(initial_integrals(nvar), final_integrals(nvar))
    if (present(failure_context)) failure_context = "initial integral"
    call composite_integral_reactive_amr_eb_patch_tree_2d( &
      solution, initial_integrals, local_ok)
    if (.not. local_ok) return
    time_tolerance = 16.0_dp * epsilon(1.0_dp) * &
      max(tiny(1.0_dp), abs(config%eb%flow%final_time))

    do
      remaining = config%eb%flow%final_time - time
      if (remaining <= time_tolerance) exit
      if (steps >= config%eb%flow%maximum_steps) then
        if (present(failure_context)) failure_context = "maximum steps"
        return
      end if
      if (present(failure_context)) failure_context = "timestep"
      call compute_reactive_amr_eb_patch_tree_timestep_2d( &
        species, transport, solution, config%eb%flow%cfl, &
        config%eb%flow%transport_cfl, &
        config%eb%flow%viscosity_enabled, &
        config%eb%flow%thermal_conduction_enabled, &
        config%eb%flow%species_diffusion_enabled, dt, local_ok)
      if (.not. local_ok) return
      dt = min(dt, remaining)
      if (present(failure_context)) failure_context = "full physics"
      call advance_reactive_amr_eb_patch_tree_full_physics_2d( &
        species, reactions, transport, solution, &
        config%eb%flow%riemann_solver, config%eb%flow%reconstruction, &
        config%eb%flow%limiter, config%eb%state_redist_max_order, dt, &
        config%eb%flow%chemistry_enabled, &
        config%eb%flow%chemistry_relative_tolerance, &
        config%eb%flow%chemistry_absolute_tolerance, &
        config%eb%flow%viscosity_enabled, &
        config%eb%flow%thermal_conduction_enabled, &
        config%eb%flow%species_diffusion_enabled, &
        config%eb%flow%barodiffusion_enabled, boundaries, &
        config%eb%state_redist_target_volume_fraction, step_theta, local_ok)
      if (.not. local_ok) return
      time = time + dt
      minimum_dt = min(minimum_dt, dt)
      local_minimum_transport_theta = min( &
        local_minimum_transport_theta, step_theta)
      steps = steps + 1

      if (config%dynamic_regridding .and. &
          modulo(steps, config%regrid_interval) == 0) then
        if (present(failure_context)) failure_context = "periodic regrid"
        call regrid_tagged_reactive_amr_eb_patch_tree_2d( &
          species, solution, criteria, config%patch_tree_maximum_levels, &
          config%refinement_ratio, build_patch_tree_geometry, local_ok, &
          changed, tagged_cells)
        if (.not. local_ok) return
        if (changed) regrids = regrids + 1
      end if
      if (config%checkpoint_interval > 0 .and. &
          modulo(steps, config%checkpoint_interval) == 0) then
        if (present(failure_context)) failure_context = "checkpoint write"
        call write_reactive_amr_eb_patch_tree_2d_checkpoint( &
          config%checkpoint_file, species, solution, time, steps, regrids, &
          minimum_dt, local_ok)
        if (.not. local_ok) return
        last_checkpoint_step = steps
        if (config%checkpoint_stop_after_write) then
          stopped_after_checkpoint = .true.
          exit
        end if
      end if
    end do

    if (.not. stopped_after_checkpoint) time = config%eb%flow%final_time
    if (len_trim(config%checkpoint_file) > 0 .and. &
        last_checkpoint_step /= steps) then
      if (present(failure_context)) failure_context = "final checkpoint write"
      call write_reactive_amr_eb_patch_tree_2d_checkpoint( &
        config%checkpoint_file, species, solution, time, steps, regrids, &
        minimum_dt, local_ok)
      if (.not. local_ok) return
    end if
    if (present(failure_context)) failure_context = "final integral"
    call composite_integral_reactive_amr_eb_patch_tree_2d( &
      solution, final_integrals, local_ok)
    if (.not. local_ok) return
    ok = solution%is_valid() .and. steps > 0 .and. &
      ieee_is_finite(minimum_dt) .and. minimum_dt > 0.0_dp
    if (present(minimum_transport_theta)) &
      minimum_transport_theta = local_minimum_transport_theta
    if (ok .and. present(failure_context)) failure_context = "none"

  contains

    subroutine build_patch_tree_geometry( &
        parent_geometry, coarse_i_lower, coarse_i_upper, coarse_j_lower, &
        coarse_j_upper, refinement_ratio, child_geometry, geometry_ok)
      type(eb_geometry_2d), intent(in) :: parent_geometry
      integer, intent(in) :: coarse_i_lower, coarse_i_upper
      integer, intent(in) :: coarse_j_lower, coarse_j_upper
      integer, intent(in) :: refinement_ratio
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

  end subroutine simulate_reactive_amr_eb_patch_tree_2d

end module reactive_eb_amr_2d_driver_mod
