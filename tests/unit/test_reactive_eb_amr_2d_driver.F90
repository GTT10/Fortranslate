program test_reactive_eb_amr_2d_driver
  use precision_mod, only: dp
  use state_indices_mod, only: irho, iet
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_cut_cell
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, composite_eb_integral_2d
  use amr_eb_regrid_2d_mod, only: reactive_eb_patch_set_2d
  use simulation_config_reactive_eb_amr_2d_mod, only: &
    reactive_eb_amr_2d_config
  use reactive_eb_amr_2d_driver_mod, only: &
    advance_two_level_reactive_eb_strang_2d, &
    compute_reactive_eb_amr_cfl_timestep_2d, &
    regrid_reactive_eb_amr_hierarchy_2d, &
    write_reactive_eb_amr_2d_checkpoint, &
    read_reactive_eb_amr_2d_checkpoint, &
    write_reactive_eb_amr_patch_set_2d_checkpoint, &
    read_reactive_eb_amr_patch_set_2d_checkpoint, &
    simulate_reactive_eb_amr_2d, &
    compute_reactive_eb_patch_set_cfl_timestep_2d, &
    simulate_reactive_eb_amr_patch_set_2d, &
    compute_three_level_reactive_eb_cfl_timestep_2d, &
    simulate_three_level_reactive_eb_amr_2d
  use reactive_eb_2d_driver_mod, only: reactive_eb_integrals_2d
  implicit none

  type(reactive_eb_amr_2d_config) :: config
  type(eb_geometry_2d) :: coarse_geometry, fine_geometry
  type(eb_geometry_2d) :: level_two_geometry
  type(eb_geometry_2d) :: checkpoint_coarse_geometry
  type(eb_geometry_2d) :: checkpoint_fine_geometry
  type(amr_eb_patch_2d) :: patch
  type(amr_eb_patch_2d) :: level_two_patch
  type(amr_eb_patch_2d) :: checkpoint_patch
  type(reactive_eb_patch_set_2d) :: multipatch_set
  type(reactive_eb_patch_set_2d) :: checkpoint_multipatch_set
  type(reactive_eb_patch_set_2d) :: empty_multipatch_set
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp), allocatable :: coarse_state(:, :, :), coarse_temperature(:, :)
  real(dp), allocatable :: fine_state(:, :, :), fine_temperature(:, :)
  real(dp), allocatable :: level_two_state(:, :, :)
  real(dp), allocatable :: level_two_temperature(:, :)
  real(dp), allocatable :: initial_integrals(:), final_integrals(:)
  real(dp), allocatable :: lifecycle_integrals(:)
  real(dp), allocatable :: rollback_coarse_state(:, :, :)
  real(dp), allocatable :: rollback_coarse_temperature(:, :)
  real(dp), allocatable :: rollback_fine_state(:, :, :)
  real(dp), allocatable :: rollback_fine_temperature(:, :)
  real(dp), allocatable :: reference_state(:)
  real(dp), allocatable :: checkpoint_coarse_state(:, :, :)
  real(dp), allocatable :: checkpoint_coarse_temperature(:, :)
  real(dp), allocatable :: checkpoint_fine_state(:, :, :)
  real(dp), allocatable :: checkpoint_fine_temperature(:, :)
  real(dp) :: time, minimum_dt, base_density, cfl_dt, conservation_error, scale
  real(dp) :: checkpoint_time, checkpoint_minimum_dt
  real(dp) :: checkpoint_base_density
  logical :: changed, fine_active, checkpoint_fine_active, ok
  integer :: initial_i_lower, initial_i_upper
  integer :: initial_j_lower, initial_j_upper, regrids, steps
  integer :: checkpoint_regrids, checkpoint_steps, child
  character(len=*), parameter :: checkpoint_path = &
    "reactive_eb_amr_2d_driver.chk"
  character(len=*), parameter :: patch_set_checkpoint_path = &
    "reactive_eb_amr_patch_set_2d_driver.chk"
  character(len=64) :: multipatch_failure_context

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamic database load")
  call load_h2o2_elementary_mechanism(reactions, ok)
  call require(ok, "elementary mechanism load")
  config%eb%flow%nx = 8
  config%eb%flow%ny = 8
  config%eb%flow%x_lower = 0.0_dp
  config%eb%flow%x_upper = 1.0_dp
  config%eb%flow%y_lower = 0.0_dp
  config%eb%flow%y_upper = 1.0_dp
  config%eb%flow%final_time = 1.0e-6_dp
  config%eb%flow%cfl = 0.2_dp
  config%eb%flow%maximum_steps = 20
  config%eb%flow%problem = "uniform_reactor"
  config%eb%flow%reconstruction = "characteristic_plm"
  config%eb%flow%limiter = "mc"
  config%eb%flow%riemann_solver = "hllc"
  config%eb%flow%use_transverse_correction = .false.
  config%eb%flow%chemistry_enabled = .false.
  config%eb%flow%transport_enabled = .false.
  config%eb%flow%boundary_x_lower = "outflow"
  config%eb%flow%boundary_x_upper = "outflow"
  config%eb%flow%boundary_y_lower = "outflow"
  config%eb%flow%boundary_y_upper = "outflow"
  config%eb%flow%initial_temperature = 1000.0_dp
  config%eb%flow%initial_pressure = 101325.0_dp
  config%eb%flow%initial_velocity_x = 0.0_dp
  config%eb%flow%initial_velocity_y = 0.0_dp
  config%eb%geometry = "plane"
  config%eb%plane_normal_x = 1.0_dp
  config%eb%plane_normal_y = 1.0_dp
  config%eb%plane_offset = 0.78_dp
  config%eb%state_redist_target_volume_fraction = 0.5_dp
  config%eb%state_redist_max_order = 2
  config%coarse_i_lower = 2
  config%coarse_i_upper = 6
  config%coarse_j_lower = 2
  config%coarse_j_upper = 6
  config%refinement_ratio = 2

  call simulate_reactive_eb_amr_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, patch, fine_active, time, &
    steps, regrids, initial_integrals, final_integrals, minimum_dt, &
    base_density, ok)
  call require(ok, "runnable static EB AMR simulation")
  call require(fine_active .and. steps == 1 .and. regrids == 0 .and. &
    time == config%eb%flow%final_time .and. &
    minimum_dt == config%eb%flow%final_time, "time-loop completion")
  call require(coarse_geometry%nx == 8 .and. coarse_geometry%ny == 8 .and. &
    fine_geometry%nx == 10 .and. fine_geometry%ny == 10, &
    "two-level dimensions")
  call require(patch%is_valid(coarse_geometry, fine_geometry), &
    "qualified static patch")
  call require(count(coarse_geometry%cell_type == eb_cut_cell) > 0 .and. &
    count(fine_geometry%cell_type == eb_cut_cell) > 0, &
    "two-level cut-cell coverage")
  scale = max(1.0_dp, maxval(abs(initial_integrals)))
  call require(maxval(abs(final_integrals - initial_integrals)) <= &
    3.0e-12_dp * scale, "static EB AMR composite conservation")
  allocate(reference_state, source=coarse_state(:, 8, 8))
  scale = max(1.0_dp, maxval(abs(reference_state)))
  call require(maxval(abs(coarse_state - &
    spread(spread(reference_state, 2, 8), 3, 8))) <= &
    3.0e-12_dp * scale, "coarse stationary state")
  call require(maxval(abs(fine_state - &
    spread(spread(reference_state, 2, 10), 3, 10))) <= &
    3.0e-12_dp * scale, "fine stationary state")
  call require(maxval(abs(coarse_temperature - 1000.0_dp)) <= 3.0e-8_dp .and. &
    maxval(abs(fine_temperature - 1000.0_dp)) <= 3.0e-8_dp, &
    "two-level stationary temperature")
  call compute_reactive_eb_amr_cfl_timestep_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, config%refinement_ratio, &
    config%eb%flow%cfl, cfl_dt, ok)
  call require(ok .and. cfl_dt > config%eb%flow%final_time, &
    "two-level CFL selection")

  call write_reactive_eb_amr_2d_checkpoint( &
    checkpoint_path, species, config, coarse_state, coarse_temperature, &
    coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
    fine_active, time, steps, regrids, minimum_dt, base_density, ok)
  call require(ok, "active fine checkpoint write")
  call read_reactive_eb_amr_2d_checkpoint( &
    checkpoint_path, species, config, checkpoint_coarse_state, &
    checkpoint_coarse_temperature, checkpoint_coarse_geometry, &
    checkpoint_fine_state, checkpoint_fine_temperature, &
    checkpoint_fine_geometry, checkpoint_patch, checkpoint_fine_active, &
    checkpoint_time, checkpoint_steps, checkpoint_regrids, &
    checkpoint_minimum_dt, checkpoint_base_density, ok)
  call require(ok .and. checkpoint_fine_active .and. &
    checkpoint_patch%is_valid( &
      checkpoint_coarse_geometry, checkpoint_fine_geometry) .and. &
    checkpoint_time == time .and. checkpoint_steps == steps .and. &
    checkpoint_regrids == regrids .and. &
    checkpoint_minimum_dt == minimum_dt .and. &
    checkpoint_base_density == base_density .and. &
    all(checkpoint_coarse_state == coarse_state) .and. &
    all(checkpoint_fine_state == fine_state), &
    "active fine checkpoint round trip")
  scale = max(1.0_dp, maxval(abs(coarse_temperature)), &
    maxval(abs(fine_temperature)))
  call require(maxval(abs(checkpoint_coarse_temperature - &
    coarse_temperature)) <= 3.0e-12_dp * scale .and. &
    maxval(abs(checkpoint_fine_temperature - fine_temperature)) <= &
    3.0e-12_dp * scale, "checkpoint EOS temperature recovery")

  allocate(rollback_coarse_state, mold=coarse_state)
  allocate(rollback_coarse_temperature, mold=coarse_temperature)
  allocate(rollback_fine_state, mold=fine_state)
  allocate(rollback_fine_temperature, mold=fine_temperature)
  call advance_two_level_reactive_eb_strang_2d( &
    species, reactions, coarse_state, coarse_temperature, coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, patch, "unknown", &
    config%eb%flow%reconstruction, config%eb%flow%limiter, &
    config%eb%state_redist_max_order, config%eb%flow%final_time, .true., &
    config%eb%flow%chemistry_relative_tolerance, &
    config%eb%flow%chemistry_absolute_tolerance, rollback_coarse_state, &
    rollback_coarse_temperature, rollback_fine_state, &
    rollback_fine_temperature, ok, &
    config%eb%state_redist_target_volume_fraction)
  call require(.not. ok .and. &
    all(rollback_coarse_state == coarse_state) .and. &
    all(rollback_coarse_temperature == coarse_temperature) .and. &
    all(rollback_fine_state == fine_state) .and. &
    all(rollback_fine_temperature == fine_temperature), &
    "two-level chemistry transaction rollback")

  config%eb%flow%nx = 12
  config%eb%flow%ny = 12
  config%eb%flow%final_time = 1.0e-8_dp
  config%eb%flow%maximum_steps = 10
  config%eb%flow%problem = "reactive_hotspot"
  config%eb%flow%hotspot_temperature_rise = 350.0_dp
  config%eb%flow%hotspot_center_x = 0.72_dp
  config%eb%flow%hotspot_center_y = 0.62_dp
  config%eb%flow%hotspot_width = 0.08_dp
  config%eb%plane_offset = 0.30_dp
  config%coarse_i_lower = 2
  config%coarse_i_upper = 5
  config%coarse_j_lower = 2
  config%coarse_j_upper = 5
  config%dynamic_regridding = .true.
  config%regrid_interval = 1
  config%regrid_relative_temperature_gradient = 0.02_dp
  config%regrid_absolute_temperature_gradient = 5.0_dp
  config%regrid_temperature_scale_floor = 1.0_dp
  config%regrid_buffer_cells = 1
  config%regrid_minimum_patch_cells_x = 4
  config%regrid_minimum_patch_cells_y = 4
  initial_i_lower = config%coarse_i_lower
  initial_i_upper = config%coarse_i_upper
  initial_j_lower = config%coarse_j_lower
  initial_j_upper = config%coarse_j_upper
  call simulate_reactive_eb_amr_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, patch, fine_active, time, &
    steps, regrids, initial_integrals, final_integrals, minimum_dt, &
    base_density, ok)
  call require(ok .and. fine_active .and. regrids >= 1, &
    "dynamic EB AMR simulation")
  call require(patch%is_valid(coarse_geometry, fine_geometry), &
    "dynamically selected patch")
  call require(patch%coarse_i_lower /= initial_i_lower .or. &
    patch%coarse_i_upper /= initial_i_upper .or. &
    patch%coarse_j_lower /= initial_j_lower .or. &
    patch%coarse_j_upper /= initial_j_upper, "hotspot moves static patch")
  scale = max(1.0_dp, maxval(abs(initial_integrals)))
  call require(maxval(abs(final_integrals - initial_integrals)) <= &
    2.0e-10_dp * scale, "dynamic EB AMR composite conservation")

  config%eb%flow%problem = "uniform_reactor"
  config%eb%flow%initial_temperature = 1000.0_dp
  config%eb%flow%initial_velocity_x = 0.0_dp
  config%eb%flow%initial_velocity_y = 0.0_dp
  config%eb%flow%final_time = 2.0e-5_dp
  config%regrid_at_initialization = .false.
  config%remove_fine_patch_when_untagged = .true.
  call simulate_reactive_eb_amr_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, patch, fine_active, time, &
    steps, regrids, initial_integrals, final_integrals, minimum_dt, &
    base_density, ok)
  call require(ok .and. .not. fine_active .and. steps == 2 .and. &
    regrids == 1, "time-loop fine-patch removal")
  call require(.not. allocated(fine_state) .and. &
    .not. allocated(fine_temperature) .and. &
    .not. fine_geometry%is_valid() .and. patch%refinement_ratio == 0, &
    "inactive fine storage released")
  scale = max(1.0_dp, maxval(abs(initial_integrals)))
  conservation_error = maxval(abs(final_integrals - initial_integrals)) / scale
  write(*, '(a,1x,es16.8)') &
    "Fine-patch removal conservation error:", conservation_error
  call require(conservation_error <= 3.0e-12_dp, &
    "fine-patch removal conservation")

  allocate(lifecycle_integrals(size(final_integrals)))
  coarse_temperature = 1000.0_dp
  coarse_temperature(9, 8) = 2000.0_dp
  call regrid_reactive_eb_amr_hierarchy_2d( &
    species, config, coarse_state, coarse_temperature, coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, patch, fine_active, &
    changed, ok)
  call require(ok .and. changed .and. fine_active .and. &
    allocated(fine_state) .and. allocated(fine_temperature), &
    "fine-patch re-creation from root-only state")
  call composite_eb_integral_2d( &
    coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
    lifecycle_integrals, ok)
  call require(ok .and. maxval(abs(lifecycle_integrals - &
    final_integrals)) <= 3.0e-12_dp * scale, &
    "fine-patch creation conservation")

  coarse_temperature = 1000.0_dp
  call regrid_reactive_eb_amr_hierarchy_2d( &
    species, config, coarse_state, coarse_temperature, coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, patch, fine_active, &
    changed, ok)
  call require(ok .and. changed .and. .not. fine_active, &
    "re-created fine patch collapses on empty tags")
  call reactive_eb_integrals_2d( &
    coarse_state, coarse_geometry, lifecycle_integrals, ok)
  call require(ok .and. maxval(abs(lifecycle_integrals - &
    final_integrals)) <= 3.0e-12_dp * scale, &
    "re-created patch collapse conservation")

  config%eb%flow%chemistry_enabled = .true.
  config%eb%flow%initial_temperature = 1200.0_dp
  config%eb%flow%final_time = 2.0e-7_dp
  config%eb%flow%cfl = 1.0e-3_dp
  call simulate_reactive_eb_amr_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
    fine_active, time, steps, regrids, initial_integrals, final_integrals, &
    minimum_dt, base_density, ok)
  call require(ok .and. .not. fine_active .and. steps >= 2 .and. &
    regrids == 1, "reacting fine-to-root lifecycle")
  call require(maxval(abs(coarse_temperature - 1200.0_dp)) > 1.0e-10_dp, &
    "EB AMR chemistry changes active state")
  scale = max(1.0_dp, abs(initial_integrals(irho)))
  call require(abs(final_integrals(irho) - initial_integrals(irho)) <= &
    3.0e-11_dp * scale, "EB AMR chemistry mass conservation")
  scale = max(1.0_dp, abs(initial_integrals(iet)))
  call require(abs(final_integrals(iet) - initial_integrals(iet)) <= &
    3.0e-11_dp * scale, "EB AMR chemistry energy conservation")

  call write_reactive_eb_amr_2d_checkpoint( &
    checkpoint_path, species, config, coarse_state, coarse_temperature, &
    coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
    fine_active, time, steps, regrids, minimum_dt, base_density, ok)
  call require(ok, "root-only checkpoint write")
  call read_reactive_eb_amr_2d_checkpoint( &
    checkpoint_path, species, config, checkpoint_coarse_state, &
    checkpoint_coarse_temperature, checkpoint_coarse_geometry, &
    checkpoint_fine_state, checkpoint_fine_temperature, &
    checkpoint_fine_geometry, checkpoint_patch, checkpoint_fine_active, &
    checkpoint_time, checkpoint_steps, checkpoint_regrids, &
    checkpoint_minimum_dt, checkpoint_base_density, ok)
  call require(ok .and. .not. checkpoint_fine_active .and. &
    .not. allocated(checkpoint_fine_state) .and. &
    .not. allocated(checkpoint_fine_temperature) .and. &
    .not. checkpoint_fine_geometry%is_valid() .and. &
    checkpoint_patch%refinement_ratio == 0 .and. &
    all(checkpoint_coarse_state == coarse_state), &
    "root-only checkpoint round trip")
  call write_truncated_checkpoint(checkpoint_path)
  call read_reactive_eb_amr_2d_checkpoint( &
    checkpoint_path, species, config, checkpoint_coarse_state, &
    checkpoint_coarse_temperature, checkpoint_coarse_geometry, &
    checkpoint_fine_state, checkpoint_fine_temperature, &
    checkpoint_fine_geometry, checkpoint_patch, checkpoint_fine_active, &
    checkpoint_time, checkpoint_steps, checkpoint_regrids, &
    checkpoint_minimum_dt, checkpoint_base_density, ok)
  call require(.not. ok .and. &
    .not. allocated(checkpoint_coarse_state) .and. &
    .not. allocated(checkpoint_coarse_temperature) .and. &
    .not. allocated(checkpoint_fine_state) .and. &
    .not. allocated(checkpoint_fine_temperature) .and. &
    .not. checkpoint_coarse_geometry%is_valid() .and. &
    .not. checkpoint_fine_geometry%is_valid(), &
    "truncated checkpoint transactional rejection")
  call delete_checkpoint(checkpoint_path)

  config%eb%flow%transport_enabled = .true.
  call simulate_reactive_eb_amr_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
    fine_active, time, steps, regrids, initial_integrals, final_integrals, &
    minimum_dt, base_density, ok)
  call require(.not. ok .and. steps == 0 .and. regrids == 0 .and. &
    time == 0.0_dp, "unsupported AMR transport rejection")

  config%eb%flow%transport_enabled = .false.
  config%eb%flow%nx = 14
  config%eb%flow%ny = 14
  config%eb%flow%x_lower = 0.0_dp
  config%eb%flow%x_upper = 1.0_dp
  config%eb%flow%y_lower = 0.0_dp
  config%eb%flow%y_upper = 1.0_dp
  config%eb%flow%problem = "reactive_double_hotspot"
  config%eb%flow%initial_temperature = 1200.0_dp
  config%eb%flow%initial_pressure = 135000.0_dp
  config%eb%flow%hotspot_temperature_rise = 20.0_dp
  config%eb%flow%hotspot_center_x = 0.25_dp
  config%eb%flow%hotspot_center_y = 0.55_dp
  config%eb%flow%hotspot_width = 0.03_dp
  config%eb%flow%second_hotspot_temperature_rise = 18.0_dp
  config%eb%flow%second_hotspot_center_x = 0.75_dp
  config%eb%flow%second_hotspot_center_y = 0.75_dp
  config%eb%flow%initial_velocity_x = 0.0_dp
  config%eb%flow%initial_velocity_y = 0.0_dp
  config%eb%flow%reconstruction = "pcm"
  config%eb%flow%chemistry_enabled = .true.
  config%eb%flow%final_time = 1.0e-8_dp
  config%eb%flow%cfl = 0.02_dp
  config%eb%flow%maximum_steps = 5
  config%eb%plane_offset = 0.78_dp
  config%eb%state_redist_max_order = 0
  config%coarse_i_lower = 2
  config%coarse_i_upper = 6
  config%coarse_j_lower = 2
  config%coarse_j_upper = 6
  config%multipatch_enabled = .true.
  config%dynamic_regridding = .true.
  config%regrid_at_initialization = .true.
  config%remove_fine_patch_when_untagged = .true.
  config%regrid_interval = 1
  config%regrid_relative_temperature_gradient = 1.0e-4_dp
  config%regrid_absolute_temperature_gradient = 0.1_dp
  config%regrid_temperature_scale_floor = 1.0_dp
  config%regrid_buffer_cells = 0
  config%regrid_minimum_patch_cells_x = 5
  config%regrid_minimum_patch_cells_y = 5
  config%regrid_maximum_patch_gap_cells = 0
  config%checkpoint_interval = 0
  config%checkpoint_stop_after_write = .false.
  config%checkpoint_file = ""
  config%restart_file = ""
  call simulate_reactive_eb_amr_patch_set_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, multipatch_set, time, steps, regrids, &
    initial_integrals, final_integrals, minimum_dt, base_density, ok, &
    multipatch_failure_context)
  write(*, '(a,l2,a,i0,a,i0,a,i0)') &
    "Public multipatch result: ok=", ok, ", patches=", &
    multipatch_set%patch_count(), ", steps=", steps, ", regrids=", regrids
  if (.not. ok) write(*, '(a,1x,a)') &
    "Public multipatch failure stage:", trim(multipatch_failure_context)
  call require(ok .and. multipatch_set%patch_count() == 2 .and. &
    multipatch_set%is_valid(coarse_geometry, size(coarse_state, 1)) .and. &
    steps == 1 .and. regrids >= 1, &
    "public multipatch EB AMR lifecycle")
  call require(count( &
    multipatch_set%children(1)%geometry%cell_type == eb_cut_cell) > 0, &
    "public multipatch EB cut-cell coverage")
  scale = max(1.0_dp, abs(initial_integrals(irho)))
  call require(abs(final_integrals(irho) - initial_integrals(irho)) <= &
    5.0e-10_dp * scale, "public multipatch mass conservation")
  scale = max(1.0_dp, abs(initial_integrals(iet)))
  call require(abs(final_integrals(iet) - initial_integrals(iet)) <= &
    5.0e-10_dp * scale, "public multipatch energy conservation")
  call compute_reactive_eb_patch_set_cfl_timestep_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, &
    multipatch_set, config%eb%flow%cfl, cfl_dt, ok)
  call require(ok .and. cfl_dt > 0.0_dp, &
    "public multipatch CFL selection")

  call write_reactive_eb_amr_patch_set_2d_checkpoint( &
    patch_set_checkpoint_path, species, config, coarse_state, &
    coarse_temperature, coarse_geometry, multipatch_set, time, steps, &
    regrids, minimum_dt, base_density, ok)
  call require(ok, "multipatch checkpoint write")
  call read_reactive_eb_amr_patch_set_2d_checkpoint( &
    patch_set_checkpoint_path, species, config, checkpoint_coarse_state, &
    checkpoint_coarse_temperature, checkpoint_coarse_geometry, &
    checkpoint_multipatch_set, checkpoint_time, checkpoint_steps, &
    checkpoint_regrids, checkpoint_minimum_dt, checkpoint_base_density, ok)
  call require(ok .and. checkpoint_multipatch_set%patch_count() == &
    multipatch_set%patch_count() .and. &
    checkpoint_multipatch_set%is_valid( &
      checkpoint_coarse_geometry, size(checkpoint_coarse_state, 1)) .and. &
    checkpoint_time == time .and. checkpoint_steps == steps .and. &
    checkpoint_regrids == regrids .and. &
    checkpoint_minimum_dt == minimum_dt .and. &
    checkpoint_base_density == base_density .and. &
    all(checkpoint_coarse_state == coarse_state), &
    "multipatch checkpoint root round trip")
  scale = max(1.0_dp, maxval(abs(coarse_temperature)))
  call require(maxval(abs(checkpoint_coarse_temperature - &
    coarse_temperature)) <= 3.0e-12_dp * scale, &
    "multipatch checkpoint root temperature recovery")
  do child = 1, multipatch_set%patch_count()
    call require( &
      checkpoint_multipatch_set%children(child)%patch%coarse_i_lower == &
        multipatch_set%children(child)%patch%coarse_i_lower .and. &
      checkpoint_multipatch_set%children(child)%patch%coarse_i_upper == &
        multipatch_set%children(child)%patch%coarse_i_upper .and. &
      checkpoint_multipatch_set%children(child)%patch%coarse_j_lower == &
        multipatch_set%children(child)%patch%coarse_j_lower .and. &
      checkpoint_multipatch_set%children(child)%patch%coarse_j_upper == &
        multipatch_set%children(child)%patch%coarse_j_upper .and. &
      checkpoint_multipatch_set%children(child)%patch%refinement_ratio == &
        multipatch_set%children(child)%patch%refinement_ratio .and. &
      all(checkpoint_multipatch_set%children(child)%state == &
        multipatch_set%children(child)%state), &
      "multipatch checkpoint child round trip")
    scale = max(1.0_dp, &
      maxval(abs(multipatch_set%children(child)%temperature)))
    call require(maxval(abs( &
      checkpoint_multipatch_set%children(child)%temperature - &
      multipatch_set%children(child)%temperature)) <= 3.0e-12_dp * scale, &
      "multipatch checkpoint child temperature recovery")
  end do
  allocate(empty_multipatch_set%children(0))
  call write_reactive_eb_amr_patch_set_2d_checkpoint( &
    patch_set_checkpoint_path, species, config, coarse_state, &
    coarse_temperature, coarse_geometry, empty_multipatch_set, time, steps, &
    regrids, minimum_dt, base_density, ok)
  call require(ok, "empty multipatch checkpoint write")
  call read_reactive_eb_amr_patch_set_2d_checkpoint( &
    patch_set_checkpoint_path, species, config, checkpoint_coarse_state, &
    checkpoint_coarse_temperature, checkpoint_coarse_geometry, &
    checkpoint_multipatch_set, checkpoint_time, checkpoint_steps, &
    checkpoint_regrids, checkpoint_minimum_dt, checkpoint_base_density, ok)
  call require(ok .and. checkpoint_multipatch_set%patch_count() == 0 .and. &
    checkpoint_multipatch_set%is_valid( &
      checkpoint_coarse_geometry, size(checkpoint_coarse_state, 1)) .and. &
    all(checkpoint_coarse_state == coarse_state), &
    "empty multipatch checkpoint round trip")
  call write_truncated_patch_set_checkpoint(patch_set_checkpoint_path)
  call read_reactive_eb_amr_patch_set_2d_checkpoint( &
    patch_set_checkpoint_path, species, config, checkpoint_coarse_state, &
    checkpoint_coarse_temperature, checkpoint_coarse_geometry, &
    checkpoint_multipatch_set, checkpoint_time, checkpoint_steps, &
    checkpoint_regrids, checkpoint_minimum_dt, checkpoint_base_density, ok)
  call require(.not. ok .and. &
    .not. allocated(checkpoint_coarse_state) .and. &
    .not. allocated(checkpoint_coarse_temperature) .and. &
    .not. checkpoint_coarse_geometry%is_valid() .and. &
    .not. allocated(checkpoint_multipatch_set%children), &
    "truncated multipatch checkpoint transactional rejection")
  call delete_checkpoint(patch_set_checkpoint_path)

  config = reactive_eb_amr_2d_config()
  config%eb%flow%nx = 8
  config%eb%flow%ny = 8
  config%eb%flow%x_lower = 0.0_dp
  config%eb%flow%x_upper = 1.0_dp
  config%eb%flow%y_lower = 0.0_dp
  config%eb%flow%y_upper = 1.0_dp
  config%eb%flow%final_time = 1.0e-8_dp
  config%eb%flow%cfl = 0.02_dp
  config%eb%flow%maximum_steps = 5
  config%eb%flow%problem = "uniform_reactor"
  config%eb%flow%reconstruction = "pcm"
  config%eb%flow%limiter = "mc"
  config%eb%flow%riemann_solver = "hllc"
  config%eb%flow%use_transverse_correction = .false.
  config%eb%flow%chemistry_enabled = .true.
  config%eb%flow%transport_enabled = .false.
  config%eb%flow%boundary_x_lower = "outflow"
  config%eb%flow%boundary_x_upper = "outflow"
  config%eb%flow%boundary_y_lower = "outflow"
  config%eb%flow%boundary_y_upper = "outflow"
  config%eb%flow%initial_temperature = 1350.0_dp
  config%eb%flow%initial_pressure = 135000.0_dp
  config%eb%flow%initial_velocity_x = 0.0_dp
  config%eb%flow%initial_velocity_y = 0.0_dp
  config%eb%geometry = "plane"
  config%eb%plane_normal_x = 1.0_dp
  config%eb%plane_normal_y = 1.0_dp
  config%eb%plane_offset = 0.78_dp
  config%eb%state_redist_target_volume_fraction = 0.5_dp
  config%eb%state_redist_max_order = 2
  config%coarse_i_lower = 2
  config%coarse_i_upper = 7
  config%coarse_j_lower = 2
  config%coarse_j_upper = 7
  config%refinement_ratio = 2
  config%three_level_enabled = .true.
  config%level_two_i_lower = 3
  config%level_two_i_upper = 10
  config%level_two_j_lower = 3
  config%level_two_j_upper = 10
  call simulate_three_level_reactive_eb_amr_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
    level_two_state, level_two_temperature, level_two_geometry, &
    level_two_patch, time, steps, initial_integrals, final_integrals, &
    minimum_dt, base_density, ok)
  call require(ok .and. steps == 1 .and. &
    time == config%eb%flow%final_time .and. &
    minimum_dt == config%eb%flow%final_time, &
    "public three-level time loop")
  call require(coarse_geometry%nx == 8 .and. coarse_geometry%ny == 8 .and. &
    fine_geometry%nx == 12 .and. fine_geometry%ny == 12 .and. &
    level_two_geometry%nx == 16 .and. level_two_geometry%ny == 16 .and. &
    patch%is_valid(coarse_geometry, fine_geometry) .and. &
    level_two_patch%is_valid(fine_geometry, level_two_geometry), &
    "public three-level hierarchy")
  call require(count(coarse_geometry%cell_type == eb_cut_cell) > 0 .and. &
    count(fine_geometry%cell_type == eb_cut_cell) > 0 .and. &
    count(level_two_geometry%cell_type == eb_cut_cell) > 0, &
    "public three-level EB coverage")
  scale = max(1.0_dp, maxval(abs(initial_integrals)))
  call require(abs(final_integrals(irho) - initial_integrals(irho)) <= &
    2.0e-8_dp * scale .and. &
    abs(final_integrals(iet) - initial_integrals(iet)) <= &
      2.0e-8_dp * scale, "public three-level conservation")
  call compute_three_level_reactive_eb_cfl_timestep_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, patch, &
    level_two_state, level_two_temperature, level_two_geometry, &
    level_two_patch, config%eb%flow%cfl, cfl_dt, ok)
  call require(ok .and. cfl_dt > 0.0_dp, &
    "public three-level CFL selection")

  write(*, '(a)') "test_reactive_eb_amr_2d_driver: PASS"

contains

  subroutine write_truncated_checkpoint(path)
    character(len=*), intent(in) :: path

    integer :: unit, status

    open(newunit=unit, file=trim(path), status="replace", action="write", &
      iostat=status)
    if (status /= 0) error stop "Could not create truncated checkpoint"
    write(unit, '(a)', iostat=status) &
      "PELEF_REACTIVE_EB_AMR_2D_CHECKPOINT"
    if (status /= 0) error stop "Could not write truncated checkpoint"
    close(unit, iostat=status)
    if (status /= 0) error stop "Could not close truncated checkpoint"
  end subroutine write_truncated_checkpoint

  subroutine write_truncated_patch_set_checkpoint(path)
    character(len=*), intent(in) :: path

    integer :: unit, status

    open(newunit=unit, file=trim(path), status="replace", action="write", &
      iostat=status)
    if (status /= 0) error stop &
      "Could not create truncated multipatch checkpoint"
    write(unit, '(a)', iostat=status) &
      "PELEF_REACTIVE_EB_AMR_PATCH_SET_2D_CHECKPOINT"
    if (status /= 0) error stop &
      "Could not write truncated multipatch checkpoint"
    close(unit, iostat=status)
    if (status /= 0) error stop &
      "Could not close truncated multipatch checkpoint"
  end subroutine write_truncated_patch_set_checkpoint

  subroutine delete_checkpoint(path)
    character(len=*), intent(in) :: path

    logical :: exists
    integer :: unit, status

    inquire(file=trim(path), exist=exists)
    if (.not. exists) return
    open(newunit=unit, file=trim(path), status="old", action="read", &
      iostat=status)
    if (status /= 0) error stop "Could not open checkpoint for deletion"
    close(unit, status="delete", iostat=status)
    if (status /= 0) error stop "Could not delete checkpoint"
  end subroutine delete_checkpoint

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_reactive_eb_amr_2d_driver
