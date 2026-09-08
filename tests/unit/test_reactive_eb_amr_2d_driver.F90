program test_reactive_eb_amr_2d_driver
  use precision_mod, only: dp
  use state_indices_mod, only: irho, iet
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_cut_cell
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, composite_eb_integral_2d
  use amr_eb_multilevel_2d_mod, only: &
    composite_three_level_eb_integral_2d
  use amr_eb_regrid_2d_mod, only: reactive_eb_patch_set_2d
  use amr_eb_patch_tree_reactive_2d_mod, only: &
    reactive_amr_eb_patch_tree_checkpoint_fingerprint_2d
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
    advance_reactive_eb_patch_set_strang_2d, &
    simulate_reactive_eb_amr_patch_set_2d, &
    compute_three_level_reactive_eb_cfl_timestep_2d, &
    regrid_three_level_reactive_eb_amr_parent_2d, &
    write_reactive_eb_amr_three_level_2d_checkpoint, &
    read_reactive_eb_amr_three_level_2d_checkpoint, &
    simulate_three_level_reactive_eb_amr_2d, &
    build_reactive_amr_eb_patch_tree_checkpoint_fingerprint_2d
  use reactive_eb_2d_driver_mod, only: reactive_eb_integrals_2d
  implicit none

  type(reactive_eb_amr_2d_config) :: config
  type(reactive_amr_eb_patch_tree_checkpoint_fingerprint_2d) :: &
    wall_fingerprint, changed_wall_fingerprint
  type(eb_geometry_2d) :: coarse_geometry, fine_geometry
  type(eb_geometry_2d) :: level_two_geometry
  type(eb_geometry_2d) :: checkpoint_coarse_geometry
  type(eb_geometry_2d) :: checkpoint_fine_geometry
  type(eb_geometry_2d) :: checkpoint_level_two_geometry
  type(amr_eb_patch_2d) :: patch
  type(amr_eb_patch_2d) :: level_two_patch
  type(amr_eb_patch_2d) :: checkpoint_patch
  type(amr_eb_patch_2d) :: checkpoint_level_two_patch
  type(reactive_eb_patch_set_2d) :: multipatch_set
  type(reactive_eb_patch_set_2d) :: checkpoint_multipatch_set
  type(reactive_eb_patch_set_2d) :: empty_multipatch_set
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
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
  real(dp), allocatable :: parent_rollback_root_state(:, :, :)
  real(dp), allocatable :: parent_rollback_root_temperature(:, :)
  real(dp), allocatable :: parent_rollback_level_one_state(:, :, :)
  real(dp), allocatable :: parent_rollback_level_one_temperature(:, :)
  real(dp), allocatable :: parent_rollback_level_two_state(:, :, :)
  real(dp), allocatable :: parent_rollback_level_two_temperature(:, :)
  real(dp), allocatable :: reference_state(:)
  real(dp), allocatable :: checkpoint_coarse_state(:, :, :)
  real(dp), allocatable :: checkpoint_coarse_temperature(:, :)
  real(dp), allocatable :: checkpoint_fine_state(:, :, :)
  real(dp), allocatable :: checkpoint_fine_temperature(:, :)
  real(dp), allocatable :: checkpoint_level_two_state(:, :, :)
  real(dp), allocatable :: checkpoint_level_two_temperature(:, :)
  real(dp), allocatable :: checkpoint_initial_integrals(:)
  real(dp), allocatable :: invalid_initial_integrals(:)
  real(dp), allocatable :: selected_mole_fractions(:)
  real(dp), allocatable :: changed_mole_fractions(:)
  real(dp), allocatable :: invalid_mole_fractions(:)
  real(dp) :: time, minimum_dt, base_density, cfl_dt, conservation_error, scale
  real(dp) :: minimum_transport_theta
  real(dp) :: checkpoint_time, checkpoint_minimum_dt
  real(dp) :: checkpoint_base_density
  logical :: changed, fine_active, checkpoint_fine_active, ok
  integer :: initial_i_lower, initial_i_upper
  integer :: initial_j_lower, initial_j_upper, regrids, steps
  integer :: parent_level_two_i_lower, parent_level_two_i_upper
  integer :: parent_level_two_j_lower, parent_level_two_j_upper
  integer :: checkpoint_regrids, checkpoint_steps, child
  character(len=*), parameter :: checkpoint_path = &
    "reactive_eb_amr_2d_driver.chk"
  character(len=*), parameter :: selected_checkpoint_path = &
    "selected_reactive_eb_amr_2d_driver.chk"
  character(len=*), parameter :: selected_dynamic_checkpoint_path = &
    "selected_dynamic_reactive_eb_amr_2d_driver.chk"
  character(len=*), parameter :: selected_dynamic_bad_baseline_path = &
    "selected_dynamic_reactive_eb_amr_2d_driver_bad_baseline.chk"
  character(len=*), parameter :: three_level_checkpoint_path = &
    "reactive_eb_amr_three_level_2d_driver.chk"
  character(len=*), parameter :: selected_three_level_checkpoint_path = &
    "selected_reactive_eb_amr_three_level_2d_driver.chk"
  character(len=*), parameter :: selected_three_level_bad_baseline_path = &
    "selected_reactive_eb_amr_three_level_2d_driver_bad_baseline.chk"
  character(len=*), parameter :: dynamic_three_level_checkpoint_path = &
    "dynamic_reactive_eb_amr_three_level_2d_driver.chk"
  character(len=*), parameter :: selected_dynamic_three_level_checkpoint_path = &
    "selected_dynamic_reactive_eb_amr_three_level_2d_driver.chk"
  character(len=*), parameter :: &
    selected_dynamic_three_level_corrupt_path = &
      "selected_dynamic_reactive_eb_amr_three_level_2d_driver_corrupt.chk"
  character(len=*), parameter :: selected_bad_context_path = &
    "selected_reactive_eb_amr_2d_driver_bad_context.chk"
  character(len=*), parameter :: selected_truncated_path = &
    "selected_reactive_eb_amr_2d_driver_truncated.chk"
  character(len=*), parameter :: selected_invalid_write_path = &
    "selected_reactive_eb_amr_2d_driver_invalid_write.chk"
  character(len=*), parameter :: selected_bundle_sha256 = &
    "f65e1c02e77618d188bc95f0868f3749d6345afb355fda924297521f69ce04c3"
  character(len=*), parameter :: changed_bundle_sha256 = &
    "e65e1c02e77618d188bc95f0868f3749d6345afb355fda924297521f69ce04c3"
  character(len=*), parameter :: patch_set_checkpoint_path = &
    "reactive_eb_amr_patch_set_2d_driver.chk"
  character(len=*), parameter :: selected_patch_set_checkpoint_path = &
    "selected_reactive_eb_amr_patch_set_2d_driver.chk"
  character(len=*), parameter :: selected_patch_set_corrupt_path = &
    "selected_reactive_eb_amr_patch_set_2d_driver_corrupt.chk"
  character(len=*), parameter :: selected_patch_set_invalid_write_path = &
    "selected_reactive_eb_amr_patch_set_2d_driver_invalid.chk"
  character(len=64) :: multipatch_failure_context, three_level_failure_context
  character(len=1024) :: checkpoint_failure_context
  character(len=8192) :: invalid_baseline_record

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamic database load")
  call load_h2o2_elementary_mechanism(reactions, ok)
  call require(ok, "elementary mechanism load")
  call load_h2o2_elementary_transport(transport, ok)
  call require(ok, "elementary transport load")
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

  allocate(selected_mole_fractions(size(species)), &
    changed_mole_fractions(size(species)), &
    invalid_mole_fractions(size(species)))
  selected_mole_fractions = 0.0_dp
  selected_mole_fractions(1) = 0.25_dp
  selected_mole_fractions(4) = 0.25_dp
  selected_mole_fractions(7) = 0.50_dp
  call write_reactive_eb_amr_2d_checkpoint( &
    selected_checkpoint_path, species, config, coarse_state, &
    coarse_temperature, coarse_geometry, fine_state, fine_temperature, &
    fine_geometry, patch, fine_active, time, steps, regrids, minimum_dt, &
    base_density, ok, bundle_sha256=selected_bundle_sha256, &
    chemistry_integrator="implicit", &
    base_mole_fractions=selected_mole_fractions, &
    failure_context=checkpoint_failure_context)
  call require(ok .and. len_trim(checkpoint_failure_context) == 0, &
    "selected two-level checkpoint write")
  call require_selected_checkpoint_context( &
    selected_checkpoint_path, selected_bundle_sha256, "implicit", &
    selected_mole_fractions)

  call delete_checkpoint(selected_invalid_write_path)
  call write_reactive_eb_amr_2d_checkpoint( &
    selected_invalid_write_path, species, config, coarse_state, &
    coarse_temperature, coarse_geometry, fine_state, fine_temperature, &
    fine_geometry, patch, fine_active, time, steps, regrids, minimum_dt, &
    base_density, ok, bundle_sha256=selected_bundle_sha256, &
    failure_context=checkpoint_failure_context)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "selected checkpoint context is incomplete") > 0, &
    "incomplete selected checkpoint write rejection")
  call require_file_absent( &
    selected_invalid_write_path, "incomplete selected write creates no file")

  call read_reactive_eb_amr_2d_checkpoint( &
    selected_checkpoint_path, species, config, checkpoint_coarse_state, &
    checkpoint_coarse_temperature, checkpoint_coarse_geometry, &
    checkpoint_fine_state, checkpoint_fine_temperature, &
    checkpoint_fine_geometry, checkpoint_patch, checkpoint_fine_active, &
    checkpoint_time, checkpoint_steps, checkpoint_regrids, &
    checkpoint_minimum_dt, checkpoint_base_density, ok, &
    failure_context=checkpoint_failure_context)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "not a fixed-runtime schema") > 0, &
    "fixed reader rejects selected checkpoint")
  call require_checkpoint_targets_empty( &
    "selected checkpoint fixed-reader rollback")

  call read_reactive_eb_amr_2d_checkpoint( &
    checkpoint_path, species, config, checkpoint_coarse_state, &
    checkpoint_coarse_temperature, checkpoint_coarse_geometry, &
    checkpoint_fine_state, checkpoint_fine_temperature, &
    checkpoint_fine_geometry, checkpoint_patch, checkpoint_fine_active, &
    checkpoint_time, checkpoint_steps, checkpoint_regrids, &
    checkpoint_minimum_dt, checkpoint_base_density, ok, &
    bundle_sha256=selected_bundle_sha256, &
    chemistry_integrator="implicit", &
    base_mole_fractions=selected_mole_fractions, &
    failure_context=checkpoint_failure_context)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "lacks selected mechanism context") > 0, &
    "selected reader rejects fixed checkpoint")
  call require_checkpoint_targets_empty( &
    "fixed checkpoint selected-reader rollback")

  call read_selected_checkpoint( &
    selected_checkpoint_path, selected_bundle_sha256, "implicit", &
    selected_mole_fractions)
  call require(ok .and. checkpoint_fine_active .and. &
    checkpoint_patch%is_valid( &
      checkpoint_coarse_geometry, checkpoint_fine_geometry) .and. &
    checkpoint_time == time .and. checkpoint_steps == steps .and. &
    checkpoint_regrids == regrids .and. &
    checkpoint_minimum_dt == minimum_dt .and. &
    checkpoint_base_density == base_density .and. &
    all(checkpoint_coarse_state == coarse_state) .and. &
    all(checkpoint_coarse_temperature == coarse_temperature) .and. &
    all(checkpoint_fine_state == fine_state) .and. &
    all(checkpoint_fine_temperature == fine_temperature), &
    "selected two-level checkpoint round trip")

  config%dynamic_regridding = .true.
  config%regrid_at_initialization = .false.
  call write_reactive_eb_amr_2d_checkpoint( &
    selected_dynamic_checkpoint_path, species, config, coarse_state, &
    coarse_temperature, coarse_geometry, fine_state, fine_temperature, &
    fine_geometry, patch, fine_active, time, steps, regrids, minimum_dt, &
    base_density, ok, bundle_sha256=selected_bundle_sha256, &
    chemistry_integrator="implicit", &
    base_mole_fractions=selected_mole_fractions, &
    initial_integrals=initial_integrals, &
    failure_context=checkpoint_failure_context)
  call require(ok .and. len_trim(checkpoint_failure_context) == 0, &
    "selected dynamic checkpoint write")
  call require_selected_checkpoint_context( &
    selected_dynamic_checkpoint_path, selected_bundle_sha256, "implicit", &
    selected_mole_fractions, 5, initial_integrals)
  call read_selected_checkpoint( &
    selected_dynamic_checkpoint_path, selected_bundle_sha256, "implicit", &
    selected_mole_fractions)
  call require(ok .and. checkpoint_fine_active .and. &
    checkpoint_time == time .and. checkpoint_steps == steps .and. &
    checkpoint_regrids == regrids .and. &
    allocated(checkpoint_initial_integrals) .and. &
    all(checkpoint_initial_integrals == initial_integrals) .and. &
    all(checkpoint_coarse_state == coarse_state) .and. &
    all(checkpoint_fine_state == fine_state), &
    "selected dynamic checkpoint round trip")
  call read_selected_checkpoint( &
    selected_checkpoint_path, selected_bundle_sha256, "implicit", &
    selected_mole_fractions)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "selected static checkpoint cannot restart dynamic hierarchy") > 0, &
    "selected dynamic reader rejects static schema")
  call require_checkpoint_targets_empty( &
    "selected static checkpoint dynamic-reader rollback")

  config%dynamic_regridding = .false.
  config%regrid_at_initialization = .true.
  call read_selected_checkpoint( &
    selected_dynamic_checkpoint_path, selected_bundle_sha256, "implicit", &
    selected_mole_fractions)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "selected dynamic checkpoint cannot restart static hierarchy") > 0, &
    "selected static reader rejects dynamic schema")
  call require_checkpoint_targets_empty( &
    "selected dynamic checkpoint static-reader rollback")
  config%dynamic_regridding = .true.
  config%regrid_at_initialization = .false.
  call replace_checkpoint_line( &
    selected_dynamic_checkpoint_path, selected_dynamic_bad_baseline_path, 8, &
    "BROKEN_DYNAMIC_BASELINE")
  call read_selected_checkpoint( &
    selected_dynamic_bad_baseline_path, selected_bundle_sha256, "implicit", &
    selected_mole_fractions)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "selected dynamic checkpoint baseline is invalid") > 0, &
    "selected dynamic checkpoint baseline rejection")
  call require_checkpoint_targets_empty( &
    "selected dynamic checkpoint baseline rollback")
  call delete_checkpoint(selected_dynamic_bad_baseline_path)
  config%dynamic_regridding = .false.
  config%regrid_at_initialization = .true.
  call delete_checkpoint(selected_dynamic_checkpoint_path)

  call read_reactive_eb_amr_2d_checkpoint( &
    selected_checkpoint_path, species, config, checkpoint_coarse_state, &
    checkpoint_coarse_temperature, checkpoint_coarse_geometry, &
    checkpoint_fine_state, checkpoint_fine_temperature, &
    checkpoint_fine_geometry, checkpoint_patch, checkpoint_fine_active, &
    checkpoint_time, checkpoint_steps, checkpoint_regrids, &
    checkpoint_minimum_dt, checkpoint_base_density, ok, &
    bundle_sha256=selected_bundle_sha256, &
    failure_context=checkpoint_failure_context)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "selected restart context is incomplete") > 0, &
    "incomplete selected restart context rejection")
  call require_checkpoint_targets_empty( &
    "incomplete selected restart context rollback")

  call read_selected_checkpoint( &
    selected_checkpoint_path, changed_bundle_sha256, "implicit", &
    selected_mole_fractions)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "bundle SHA-256 mismatch") > 0, &
    "selected checkpoint bundle mismatch rejection")
  call require_checkpoint_targets_empty( &
    "selected checkpoint bundle mismatch rollback")

  call read_selected_checkpoint( &
    selected_checkpoint_path, selected_bundle_sha256, "explicit", &
    selected_mole_fractions)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "chemistry integrator mismatch") > 0, &
    "selected checkpoint integrator mismatch rejection")
  call require_checkpoint_targets_empty( &
    "selected checkpoint integrator mismatch rollback")

  changed_mole_fractions = selected_mole_fractions
  changed_mole_fractions(1) = 0.24_dp
  changed_mole_fractions(4) = 0.26_dp
  call read_selected_checkpoint( &
    selected_checkpoint_path, selected_bundle_sha256, "implicit", &
    changed_mole_fractions)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "composition mismatch") > 0, &
    "selected checkpoint composition mismatch rejection")
  call require_checkpoint_targets_empty( &
    "selected checkpoint composition mismatch rollback")

  call replace_checkpoint_line( &
    selected_checkpoint_path, selected_bad_context_path, 3, &
    "BROKEN_SELECTED_CONTEXT")
  call read_selected_checkpoint( &
    selected_bad_context_path, selected_bundle_sha256, "implicit", &
    selected_mole_fractions)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "selected checkpoint context is invalid") > 0, &
    "selected checkpoint marker rejection")
  call require_checkpoint_targets_empty("selected checkpoint marker rollback")

  call make_prefix_checkpoint( &
    selected_checkpoint_path, selected_truncated_path, 5)
  call read_selected_checkpoint( &
    selected_truncated_path, selected_bundle_sha256, "implicit", &
    selected_mole_fractions)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "selected checkpoint context is invalid") > 0, &
    "selected checkpoint context truncation rejection")
  call require_checkpoint_targets_empty( &
    "selected checkpoint context truncation rollback")

  invalid_mole_fractions = selected_mole_fractions
  invalid_mole_fractions(1) = -0.25_dp
  invalid_mole_fractions(7) = 1.0_dp
  call write_reactive_eb_amr_2d_checkpoint( &
    selected_checkpoint_path, species, config, coarse_state, &
    coarse_temperature, coarse_geometry, fine_state, fine_temperature, &
    fine_geometry, patch, fine_active, time, steps, regrids, minimum_dt, &
    base_density, ok, bundle_sha256=selected_bundle_sha256, &
    chemistry_integrator="implicit", &
    base_mole_fractions=invalid_mole_fractions, &
    failure_context=checkpoint_failure_context)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "selected checkpoint context is invalid") > 0, &
    "invalid selected checkpoint write rejection")
  call read_selected_checkpoint( &
    selected_checkpoint_path, selected_bundle_sha256, "implicit", &
    selected_mole_fractions)
  call require(ok .and. all(checkpoint_coarse_state == coarse_state) .and. &
    all(checkpoint_fine_state == fine_state), &
    "invalid selected write preserves existing checkpoint")

  call delete_checkpoint(selected_checkpoint_path)
  call delete_checkpoint(selected_bad_context_path)
  call delete_checkpoint(selected_truncated_path)

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

  config%dynamic_regridding = .false.
  config%remove_fine_patch_when_untagged = .false.
  config%eb%flow%chemistry_enabled = .false.
  config%eb%flow%maximum_steps = 20
  config%eb%flow%transport_enabled = .true.
  call simulate_reactive_eb_amr_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
    fine_active, time, steps, regrids, initial_integrals, final_integrals, &
    minimum_dt, base_density, ok)
  call require(.not. ok .and. steps == 0 .and. regrids == 0 .and. &
    time == 0.0_dp, "missing AMR transport database rejection")

  config%eb%embedded_wall_thermal = "isothermal"
  config%eb%embedded_wall_temperature = 1500.0_dp
  call build_reactive_amr_eb_patch_tree_checkpoint_fingerprint_2d( &
    config, wall_fingerprint, ok)
  call require(ok .and. &
    trim(wall_fingerprint%embedded_wall_thermal) == "isothermal" .and. &
    wall_fingerprint%embedded_wall_values(1) == 1500.0_dp, &
    "patch-tree embedded-wall checkpoint fingerprint")
  changed_wall_fingerprint = wall_fingerprint
  changed_wall_fingerprint%embedded_wall_values(1) = 1510.0_dp
  call require(.not. wall_fingerprint%matches(changed_wall_fingerprint), &
    "patch-tree wall-temperature fingerprint mismatch")
  config%prolongation_method = "linear"
  call build_reactive_amr_eb_patch_tree_checkpoint_fingerprint_2d( &
    config, changed_wall_fingerprint, ok)
  call require(ok .and. &
    trim(changed_wall_fingerprint%prolongation_method) == "linear" .and. &
    .not. wall_fingerprint%matches(changed_wall_fingerprint), &
    "patch-tree prolongation fingerprint mismatch")
  call simulate_reactive_eb_amr_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
    fine_active, time, steps, regrids, initial_integrals, final_integrals, &
    minimum_dt, base_density, ok, transport, minimum_transport_theta)
  call require(ok .and. steps >= 1 .and. &
    time == config%eb%flow%final_time, &
    "runnable two-level AMR isothermal-wall transport")
  call require(final_integrals(iet) > initial_integrals(iet), &
    "two-level AMR isothermal-wall heating")
  config%eb%embedded_wall_thermal = "adiabatic"
  call simulate_reactive_eb_amr_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
    fine_active, time, steps, regrids, initial_integrals, final_integrals, &
    minimum_dt, base_density, ok, transport, minimum_transport_theta)
  conservation_error = maxval(abs(final_integrals - initial_integrals) / &
    max(1.0_dp, abs(initial_integrals)))
  call require(ok, "runnable two-level AMR transport")
  call require(steps >= 1 .and. &
    time == config%eb%flow%final_time .and. minimum_dt > 0.0_dp, &
    "two-level AMR transport time loop")
  call require(fine_active, "two-level AMR transport fine patch")
  call require(minimum_transport_theta > 0.999999999_dp, &
    "two-level AMR transport limiter")
  call require(conservation_error <= 8.0e-11_dp, &
    "two-level AMR transport conservation")

  config%eb%embedded_wall_thermal = "isothermal"
  config%checkpoint_file = checkpoint_path
  call simulate_reactive_eb_amr_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
    fine_active, time, steps, regrids, initial_integrals, final_integrals, &
    minimum_dt, base_density, ok, transport, minimum_transport_theta)
  call require(ok .and. steps >= 1 .and. &
    time == config%eb%flow%final_time, &
    "linear-prolongation checkpoint write")
  config%checkpoint_file = ""
  config%restart_file = checkpoint_path
  call simulate_reactive_eb_amr_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
    fine_active, time, steps, regrids, initial_integrals, final_integrals, &
    minimum_dt, base_density, ok, transport, minimum_transport_theta)
  call require(ok .and. steps >= 1 .and. &
    time == config%eb%flow%final_time, &
    "linear-prolongation checkpoint restart")
  config%prolongation_method = "pcm"
  call simulate_reactive_eb_amr_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
    fine_active, time, steps, regrids, initial_integrals, final_integrals, &
    minimum_dt, base_density, ok, transport, minimum_transport_theta)
  call require(.not. ok .and. steps == 0 .and. time == 0.0_dp, &
    "prolongation-method checkpoint mismatch rejection")
  config%restart_file = ""
  config%checkpoint_file = checkpoint_path
  call simulate_reactive_eb_amr_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
    fine_active, time, steps, regrids, initial_integrals, final_integrals, &
    minimum_dt, base_density, ok, transport, minimum_transport_theta)
  call require(ok .and. steps >= 1 .and. &
    time == config%eb%flow%final_time, &
    "AMR wall-transport checkpoint write")
  config%checkpoint_file = ""
  config%restart_file = checkpoint_path
  call simulate_reactive_eb_amr_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
    fine_active, time, steps, regrids, initial_integrals, final_integrals, &
    minimum_dt, base_density, ok, transport, minimum_transport_theta)
  call require(ok .and. steps >= 1 .and. &
    time == config%eb%flow%final_time, &
    "AMR wall-transport checkpoint restart")
  config%eb%embedded_wall_temperature = 1510.0_dp
  call simulate_reactive_eb_amr_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
    fine_active, time, steps, regrids, initial_integrals, final_integrals, &
    minimum_dt, base_density, ok, transport, minimum_transport_theta)
  call require(.not. ok .and. steps == 0 .and. time == 0.0_dp, &
    "AMR wall-temperature checkpoint mismatch rejection")
  config%restart_file = ""
  config%eb%embedded_wall_thermal = "adiabatic"
  config%eb%embedded_wall_temperature = 300.0_dp
  call delete_checkpoint(checkpoint_path)

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
  config%prolongation_method = "linear"
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

  allocate(checkpoint_coarse_state, mold=coarse_state)
  allocate(checkpoint_coarse_temperature, mold=coarse_temperature)
  call advance_reactive_eb_patch_set_strang_2d( &
    species, reactions, coarse_state, coarse_temperature, coarse_geometry, &
    multipatch_set, config%eb%flow%riemann_solver, &
    config%eb%flow%reconstruction, config%eb%flow%limiter, &
    config%eb%state_redist_max_order, 1.0e-12_dp, .true., &
    config%eb%flow%chemistry_relative_tolerance, &
    config%eb%flow%chemistry_absolute_tolerance, checkpoint_coarse_state, &
    checkpoint_coarse_temperature, checkpoint_multipatch_set, ok, &
    config%eb%state_redist_target_volume_fraction, &
    multipatch_failure_context, chemistry_integrator="invalid")
  call require(.not. ok .and. index(multipatch_failure_context, &
    "first coarse chemistry half-step") > 0, &
    "multipatch chemistry-integrator propagation")
  deallocate(checkpoint_coarse_state, checkpoint_coarse_temperature)
  checkpoint_multipatch_set = reactive_eb_patch_set_2d()

  config%prolongation_method = "pcm"
  call simulate_reactive_eb_amr_patch_set_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, multipatch_set, time, steps, regrids, &
    initial_integrals, final_integrals, minimum_dt, base_density, ok, &
    multipatch_failure_context)
  call require(ok .and. multipatch_set%patch_count() == 2, &
    "public multipatch PCM checkpoint baseline")
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

  call write_reactive_eb_amr_patch_set_2d_checkpoint( &
    selected_patch_set_checkpoint_path, species, config, coarse_state, &
    coarse_temperature, coarse_geometry, multipatch_set, time, steps, &
    regrids, minimum_dt, base_density, ok, &
    bundle_sha256=selected_bundle_sha256, &
    chemistry_integrator="implicit", &
    base_mole_fractions=selected_mole_fractions, &
    initial_integrals=initial_integrals, &
    failure_context=checkpoint_failure_context)
  call require(ok .and. len_trim(checkpoint_failure_context) == 0, &
    "selected multipatch checkpoint write")
  call require_selected_patch_set_checkpoint_context( &
    selected_patch_set_checkpoint_path, selected_bundle_sha256, &
    "implicit", selected_mole_fractions, initial_integrals, &
    multipatch_set%patch_count())
  call read_reactive_eb_amr_patch_set_2d_checkpoint( &
    selected_patch_set_checkpoint_path, species, config, &
    checkpoint_coarse_state, checkpoint_coarse_temperature, &
    checkpoint_coarse_geometry, checkpoint_multipatch_set, checkpoint_time, &
    checkpoint_steps, checkpoint_regrids, checkpoint_minimum_dt, &
    checkpoint_base_density, ok, bundle_sha256=selected_bundle_sha256, &
    chemistry_integrator="implicit", &
    base_mole_fractions=selected_mole_fractions, &
    initial_integrals=checkpoint_initial_integrals, &
    failure_context=checkpoint_failure_context)
  call require(ok .and. &
    checkpoint_multipatch_set%patch_count() == &
      multipatch_set%patch_count() .and. &
    all(checkpoint_coarse_state == coarse_state) .and. &
    all(checkpoint_initial_integrals == initial_integrals), &
    "selected multipatch checkpoint round trip")

  deallocate(checkpoint_initial_integrals)
  call read_reactive_eb_amr_patch_set_2d_checkpoint( &
    selected_patch_set_checkpoint_path, species, config, &
    checkpoint_coarse_state, checkpoint_coarse_temperature, &
    checkpoint_coarse_geometry, checkpoint_multipatch_set, checkpoint_time, &
    checkpoint_steps, checkpoint_regrids, checkpoint_minimum_dt, &
    checkpoint_base_density, ok, failure_context=checkpoint_failure_context)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "requires selected reader") > 0, &
    "selected multipatch checkpoint fixed-reader rejection")
  call require_patch_set_checkpoint_targets_empty( &
    "selected multipatch fixed-reader rollback")
  call read_reactive_eb_amr_patch_set_2d_checkpoint( &
    patch_set_checkpoint_path, species, config, checkpoint_coarse_state, &
    checkpoint_coarse_temperature, checkpoint_coarse_geometry, &
    checkpoint_multipatch_set, checkpoint_time, checkpoint_steps, &
    checkpoint_regrids, checkpoint_minimum_dt, checkpoint_base_density, ok, &
    bundle_sha256=selected_bundle_sha256, &
    chemistry_integrator="implicit", &
    base_mole_fractions=selected_mole_fractions, &
    initial_integrals=checkpoint_initial_integrals, &
    failure_context=checkpoint_failure_context)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "lacks selected context") > 0, &
    "fixed multipatch checkpoint selected-reader rejection")
  call require_patch_set_checkpoint_targets_empty( &
    "fixed multipatch selected-reader rollback")

  call read_reactive_eb_amr_patch_set_2d_checkpoint( &
    selected_patch_set_checkpoint_path, species, config, &
    checkpoint_coarse_state, checkpoint_coarse_temperature, &
    checkpoint_coarse_geometry, checkpoint_multipatch_set, checkpoint_time, &
    checkpoint_steps, checkpoint_regrids, checkpoint_minimum_dt, &
    checkpoint_base_density, ok, bundle_sha256=selected_bundle_sha256, &
    chemistry_integrator="implicit", &
    base_mole_fractions=changed_mole_fractions, &
    initial_integrals=checkpoint_initial_integrals, &
    failure_context=checkpoint_failure_context)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "composition mismatch") > 0, &
    "selected multipatch composition mismatch rejection")
  call require_patch_set_checkpoint_targets_empty( &
    "selected multipatch composition rollback")

  call replace_checkpoint_record( &
    selected_patch_set_checkpoint_path, selected_patch_set_corrupt_path, &
    "END_CHECKPOINT", "CORRUPT_END")
  call read_reactive_eb_amr_patch_set_2d_checkpoint( &
    selected_patch_set_corrupt_path, species, config, checkpoint_coarse_state, &
    checkpoint_coarse_temperature, checkpoint_coarse_geometry, &
    checkpoint_multipatch_set, checkpoint_time, checkpoint_steps, &
    checkpoint_regrids, checkpoint_minimum_dt, checkpoint_base_density, ok, &
    bundle_sha256=selected_bundle_sha256, &
    chemistry_integrator="implicit", &
    base_mole_fractions=selected_mole_fractions, &
    initial_integrals=checkpoint_initial_integrals, &
    failure_context=checkpoint_failure_context)
  call require(.not. ok, &
    "selected multipatch corrupt end-marker rejection")
  call require_patch_set_checkpoint_targets_empty( &
    "selected multipatch corrupt end-marker rollback")
  call append_checkpoint_record( &
    selected_patch_set_checkpoint_path, selected_patch_set_corrupt_path, &
    "TRAILING_CONTENT")
  call read_reactive_eb_amr_patch_set_2d_checkpoint( &
    selected_patch_set_corrupt_path, species, config, checkpoint_coarse_state, &
    checkpoint_coarse_temperature, checkpoint_coarse_geometry, &
    checkpoint_multipatch_set, checkpoint_time, checkpoint_steps, &
    checkpoint_regrids, checkpoint_minimum_dt, checkpoint_base_density, ok, &
    bundle_sha256=selected_bundle_sha256, &
    chemistry_integrator="implicit", &
    base_mole_fractions=selected_mole_fractions, &
    initial_integrals=checkpoint_initial_integrals, &
    failure_context=checkpoint_failure_context)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "trailing content") > 0, &
    "selected multipatch trailing-content rejection")
  call require_patch_set_checkpoint_targets_empty( &
    "selected multipatch trailing-content rollback")

  allocate(invalid_initial_integrals, source=initial_integrals)
  invalid_initial_integrals(iet + 1) = &
    invalid_initial_integrals(iet + 1) + 1.0e-3_dp * &
      invalid_initial_integrals(irho)
  call delete_checkpoint(selected_patch_set_invalid_write_path)
  call write_reactive_eb_amr_patch_set_2d_checkpoint( &
    selected_patch_set_invalid_write_path, species, config, coarse_state, &
    coarse_temperature, coarse_geometry, multipatch_set, time, steps, &
    regrids, minimum_dt, base_density, ok, &
    bundle_sha256=selected_bundle_sha256, &
    chemistry_integrator="implicit", &
    base_mole_fractions=selected_mole_fractions, &
    initial_integrals=invalid_initial_integrals, &
    failure_context=checkpoint_failure_context)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "baseline is invalid") > 0, &
    "selected multipatch invalid baseline rejection")
  call require_file_absent(selected_patch_set_invalid_write_path, &
    "selected multipatch invalid baseline creates no file")
  deallocate(invalid_initial_integrals)

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
  call delete_checkpoint(selected_patch_set_checkpoint_path)
  call delete_checkpoint(selected_patch_set_corrupt_path)

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
  config%prolongation_method = "linear"
  config%three_level_enabled = .true.
  config%level_two_i_lower = 3
  config%level_two_i_upper = 10
  config%level_two_j_lower = 3
  config%level_two_j_upper = 10
  call simulate_three_level_reactive_eb_amr_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
      level_two_state, level_two_temperature, level_two_geometry, &
      level_two_patch, time, steps, regrids, initial_integrals, &
      final_integrals, minimum_dt, base_density, ok)
  call require(ok .and. steps == 1 .and. regrids == 0 .and. &
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

  call write_reactive_eb_amr_three_level_2d_checkpoint( &
    three_level_checkpoint_path, species, config, coarse_state, &
    coarse_temperature, coarse_geometry, fine_state, fine_temperature, &
    fine_geometry, patch, level_two_state, level_two_temperature, &
    level_two_geometry, level_two_patch, time, steps, minimum_dt, &
    base_density, ok, regrids, failure_context=checkpoint_failure_context)
  call require(ok .and. len_trim(checkpoint_failure_context) == 0, &
    "fixed three-level checkpoint write")
  call write_reactive_eb_amr_three_level_2d_checkpoint( &
    selected_three_level_checkpoint_path, species, config, coarse_state, &
    coarse_temperature, coarse_geometry, fine_state, fine_temperature, &
    fine_geometry, patch, level_two_state, level_two_temperature, &
    level_two_geometry, level_two_patch, time, steps, minimum_dt, &
    base_density, ok, regrids, bundle_sha256=selected_bundle_sha256, &
    chemistry_integrator="implicit", &
    base_mole_fractions=selected_mole_fractions, &
    initial_integrals=initial_integrals, &
    failure_context=checkpoint_failure_context)
  call require(ok .and. len_trim(checkpoint_failure_context) == 0, &
    "selected three-level checkpoint write")
  call require_selected_three_level_checkpoint_context( &
    selected_three_level_checkpoint_path, selected_bundle_sha256, &
    "implicit", selected_mole_fractions, initial_integrals)
  call read_selected_three_level_checkpoint( &
    selected_three_level_checkpoint_path, selected_bundle_sha256, &
    "implicit", selected_mole_fractions)
  call require(ok .and. checkpoint_time == time .and. &
    checkpoint_steps == steps .and. checkpoint_regrids == regrids .and. &
    checkpoint_minimum_dt == minimum_dt .and. &
    checkpoint_base_density == base_density .and. &
    checkpoint_patch%is_valid( &
      checkpoint_coarse_geometry, checkpoint_fine_geometry) .and. &
    checkpoint_level_two_patch%is_valid( &
      checkpoint_fine_geometry, checkpoint_level_two_geometry) .and. &
    all(checkpoint_coarse_state == coarse_state) .and. &
    all(checkpoint_fine_state == fine_state) .and. &
    all(checkpoint_level_two_state == level_two_state) .and. &
    allocated(checkpoint_initial_integrals) .and. &
    all(checkpoint_initial_integrals == initial_integrals), &
    "selected three-level checkpoint round trip")

  allocate(invalid_initial_integrals, source=initial_integrals)
  invalid_initial_integrals(irho) = -1.0_dp
  call write_reactive_eb_amr_three_level_2d_checkpoint( &
    selected_three_level_checkpoint_path, species, config, coarse_state, &
    coarse_temperature, coarse_geometry, fine_state, fine_temperature, &
    fine_geometry, patch, level_two_state, level_two_temperature, &
    level_two_geometry, level_two_patch, time, steps, minimum_dt, &
    base_density, ok, regrids, bundle_sha256=selected_bundle_sha256, &
    chemistry_integrator="implicit", &
    base_mole_fractions=selected_mole_fractions, &
    initial_integrals=invalid_initial_integrals, &
    failure_context=checkpoint_failure_context)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "selected checkpoint baseline is invalid") > 0, &
    "invalid selected three-level checkpoint write rejection")
  call read_selected_three_level_checkpoint( &
    selected_three_level_checkpoint_path, selected_bundle_sha256, &
    "implicit", selected_mole_fractions)
  call require(ok .and. all(checkpoint_coarse_state == coarse_state) .and. &
    all(checkpoint_fine_state == fine_state) .and. &
    all(checkpoint_level_two_state == level_two_state) .and. &
    all(checkpoint_initial_integrals == initial_integrals), &
    "invalid selected three-level write preserves existing checkpoint")

  if (allocated(checkpoint_initial_integrals)) &
    deallocate(checkpoint_initial_integrals)
  call read_reactive_eb_amr_three_level_2d_checkpoint( &
    selected_three_level_checkpoint_path, species, config, &
    checkpoint_coarse_state, checkpoint_coarse_temperature, &
    checkpoint_coarse_geometry, checkpoint_fine_state, &
    checkpoint_fine_temperature, checkpoint_fine_geometry, &
    checkpoint_patch, checkpoint_level_two_state, &
    checkpoint_level_two_temperature, checkpoint_level_two_geometry, &
    checkpoint_level_two_patch, checkpoint_time, checkpoint_steps, &
    checkpoint_minimum_dt, checkpoint_base_density, ok, &
    checkpoint_regrids, failure_context=checkpoint_failure_context)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "not a fixed-runtime schema") > 0, &
    "fixed reader rejects selected three-level checkpoint")
  call require_three_level_checkpoint_targets_empty( &
    "selected three-level checkpoint fixed-reader rollback")

  call read_selected_three_level_checkpoint( &
    three_level_checkpoint_path, selected_bundle_sha256, "implicit", &
    selected_mole_fractions)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "lacks selected mechanism context") > 0, &
    "selected reader rejects fixed three-level checkpoint")
  call require_three_level_checkpoint_targets_empty( &
    "fixed three-level checkpoint selected-reader rollback")

  call read_selected_three_level_checkpoint( &
    selected_three_level_checkpoint_path, changed_bundle_sha256, &
    "implicit", selected_mole_fractions)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "bundle SHA-256 mismatch") > 0, &
    "selected three-level bundle mismatch rejection")
  call require_three_level_checkpoint_targets_empty( &
    "selected three-level bundle mismatch rollback")
  call read_selected_three_level_checkpoint( &
    selected_three_level_checkpoint_path, selected_bundle_sha256, &
    "explicit", selected_mole_fractions)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "chemistry integrator mismatch") > 0, &
    "selected three-level integrator mismatch rejection")
  call require_three_level_checkpoint_targets_empty( &
    "selected three-level integrator mismatch rollback")
  call read_selected_three_level_checkpoint( &
    selected_three_level_checkpoint_path, selected_bundle_sha256, &
    "implicit", changed_mole_fractions)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "composition mismatch") > 0, &
    "selected three-level composition mismatch rejection")
  call require_three_level_checkpoint_targets_empty( &
    "selected three-level composition mismatch rollback")

  call replace_checkpoint_line( &
    selected_three_level_checkpoint_path, &
    selected_three_level_bad_baseline_path, 8, "BROKEN_COMPOSITE_BASELINE")
  call read_selected_three_level_checkpoint( &
    selected_three_level_bad_baseline_path, selected_bundle_sha256, &
    "implicit", selected_mole_fractions)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "selected checkpoint baseline is invalid") > 0, &
    "selected three-level baseline marker rejection")
  call require_three_level_checkpoint_targets_empty( &
    "selected three-level baseline marker rollback")
  write(invalid_baseline_record, '(*(es27.18e3,1x))') &
    invalid_initial_integrals
  call replace_checkpoint_line( &
    selected_three_level_checkpoint_path, &
    selected_three_level_bad_baseline_path, 10, &
    trim(invalid_baseline_record))
  call read_selected_three_level_checkpoint( &
    selected_three_level_bad_baseline_path, selected_bundle_sha256, &
    "implicit", selected_mole_fractions)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "selected checkpoint baseline is invalid") > 0, &
    "selected three-level numeric baseline rejection")
  call require_three_level_checkpoint_targets_empty( &
    "selected three-level numeric baseline rollback")
  deallocate(invalid_initial_integrals)
  call delete_checkpoint(selected_three_level_bad_baseline_path)

  config%eb%flow%chemistry_enabled = .false.
  config%eb%flow%transport_enabled = .true.
  config%eb%flow%viscosity_enabled = .false.
  config%eb%flow%thermal_conduction_enabled = .true.
  config%eb%flow%species_diffusion_enabled = .false.
  config%eb%flow%barodiffusion_enabled = .false.
  config%eb%flow%problem = "reactive_hotspot"
  config%eb%flow%hotspot_center_x = 0.5_dp
  config%eb%flow%hotspot_center_y = 0.5_dp
  config%eb%flow%hotspot_width = 0.08_dp
  config%eb%flow%hotspot_temperature_rise = 250.0_dp
  call simulate_three_level_reactive_eb_amr_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
    level_two_state, level_two_temperature, level_two_geometry, &
    level_two_patch, time, steps, regrids, initial_integrals, &
    final_integrals, minimum_dt, base_density, ok)
  call require(.not. ok .and. steps == 0 .and. time == 0.0_dp, &
    "missing three-level transport database rejection")
  call simulate_three_level_reactive_eb_amr_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
    level_two_state, level_two_temperature, level_two_geometry, &
    level_two_patch, time, steps, regrids, initial_integrals, &
    final_integrals, minimum_dt, base_density, ok, transport=transport, &
    minimum_transport_theta=minimum_transport_theta)
  scale = max(1.0_dp, maxval(abs(initial_integrals)))
  conservation_error = maxval(abs(final_integrals - initial_integrals)) / scale
  call require(ok .and. steps >= 1 .and. &
    time == config%eb%flow%final_time .and. minimum_dt > 0.0_dp, &
    "public three-level transport time loop")
  call require(minimum_transport_theta > 0.999999999_dp, &
    "public three-level transport limiter")
  call require(conservation_error <= 2.0e-8_dp, &
    "public three-level transport conservation")

  config%eb%flow%chemistry_enabled = .true.
  config%eb%flow%transport_enabled = .true.
  config%eb%flow%problem = "reactive_hotspot"
  config%eb%flow%hotspot_center_x = 0.5_dp
  config%eb%flow%hotspot_center_y = 0.5_dp
  config%eb%flow%hotspot_width = 0.08_dp
  config%eb%flow%hotspot_temperature_rise = 250.0_dp
  config%eb%flow%final_time = 1.0e-8_dp
  config%level_two_i_upper = 6
  config%level_two_j_upper = 6
  config%dynamic_regridding = .true.
  config%regrid_at_initialization = .true.
  config%remove_fine_patch_when_untagged = .false.
  config%regrid_interval = 1
  config%regrid_relative_temperature_gradient = 1.0e-4_dp
  config%regrid_absolute_temperature_gradient = 0.1_dp
  config%regrid_temperature_scale_floor = 1.0_dp
  config%regrid_buffer_cells = 0
  config%regrid_minimum_patch_cells_x = 4
  config%regrid_minimum_patch_cells_y = 4
  call simulate_three_level_reactive_eb_amr_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
      level_two_state, level_two_temperature, level_two_geometry, &
      level_two_patch, time, steps, regrids, initial_integrals, &
      final_integrals, minimum_dt, base_density, ok, &
      failure_context=three_level_failure_context, transport=transport, &
      minimum_transport_theta=minimum_transport_theta)
  if (.not. ok .or. steps <= 0 .or. regrids <= 0 .or. &
      .not. level_two_patch%is_valid(fine_geometry, level_two_geometry) .or. &
      level_two_patch%coarse_i_lower < 3 .or. &
      level_two_patch%coarse_i_upper > fine_geometry%nx - 2 .or. &
      level_two_patch%coarse_j_lower < 3 .or. &
      level_two_patch%coarse_j_upper > fine_geometry%ny - 2 .or. &
      (level_two_patch%coarse_i_upper == 6 .and. &
       level_two_patch%coarse_j_upper == 6)) then
    write(*, *) "Dynamic three-level regrid diagnostics:", &
      trim(three_level_failure_context), ok, steps, regrids, &
      level_two_patch%coarse_i_lower, level_two_patch%coarse_i_upper, &
      level_two_patch%coarse_j_lower, level_two_patch%coarse_j_upper
  end if
  call require(ok .and. steps > 0 .and. regrids > 0 .and. &
    level_two_patch%is_valid(fine_geometry, level_two_geometry) .and. &
    level_two_patch%coarse_i_lower >= 3 .and. &
    level_two_patch%coarse_i_upper <= fine_geometry%nx - 2 .and. &
    level_two_patch%coarse_j_lower >= 3 .and. &
    level_two_patch%coarse_j_upper <= fine_geometry%ny - 2 .and. &
    (level_two_patch%coarse_i_upper /= 6 .or. &
     level_two_patch%coarse_j_upper /= 6), &
    "public dynamic three-level finest regrid")
  call require(minimum_transport_theta > 0.999999999_dp, &
    "public dynamic three-level transport limiter")

  config%eb%flow%chemistry_enabled = .false.
  config%eb%flow%transport_enabled = .false.
  config%eb%flow%hotspot_center_x = 0.72_dp
  config%eb%flow%hotspot_center_y = 0.68_dp
  config%dynamic_regridding = .false.
  config%regrid_at_initialization = .false.
  config%regrid_buffer_cells = 1
  call simulate_three_level_reactive_eb_amr_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
    level_two_state, level_two_temperature, level_two_geometry, &
    level_two_patch, time, steps, regrids, initial_integrals, &
    final_integrals, minimum_dt, base_density, ok)
  call require(ok .and. steps > 0 .and. regrids == 0, &
    "static three-level parent-regrid seed")
  initial_i_lower = patch%coarse_i_lower
  initial_i_upper = patch%coarse_i_upper
  initial_j_lower = patch%coarse_j_lower
  initial_j_upper = patch%coarse_j_upper
  call composite_three_level_eb_integral_2d( &
    coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
    level_two_state, level_two_geometry, level_two_patch, &
    initial_integrals, ok)
  call require(ok, "three-level parent-regrid initial integral")
  config%dynamic_regridding = .true.
  call regrid_three_level_reactive_eb_amr_parent_2d( &
    species, config, coarse_state, coarse_temperature, coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, patch, level_two_state, &
    level_two_temperature, level_two_geometry, level_two_patch, changed, ok)
  call require(ok .and. changed .and. &
    patch%is_valid(coarse_geometry, fine_geometry) .and. &
    level_two_patch%is_valid(fine_geometry, level_two_geometry) .and. &
    (patch%coarse_i_lower /= initial_i_lower .or. &
     patch%coarse_i_upper /= initial_i_upper .or. &
     patch%coarse_j_lower /= initial_j_lower .or. &
     patch%coarse_j_upper /= initial_j_upper) .and. &
    patch%coarse_i_upper == coarse_geometry%nx .and. &
    patch%coarse_j_upper == coarse_geometry%ny .and. &
    level_two_patch%coarse_i_lower >= 3 .and. &
    level_two_patch%coarse_i_upper <= fine_geometry%nx - 2 .and. &
    level_two_patch%coarse_j_lower >= 3 .and. &
    level_two_patch%coarse_j_upper <= fine_geometry%ny - 2, &
    "transactional three-level parent regrid")
  call composite_three_level_eb_integral_2d( &
    coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
    level_two_state, level_two_geometry, level_two_patch, &
    final_integrals, ok)
  scale = max(1.0_dp, maxval(abs(initial_integrals)))
  call require(ok .and. maxval(abs(final_integrals - &
    initial_integrals)) <= 1.0e-10_dp * scale, &
    "three-level parent-regrid conservation")

  allocate(parent_rollback_root_state, source=coarse_state)
  allocate(parent_rollback_root_temperature, source=coarse_temperature)
  allocate(parent_rollback_level_one_state, source=fine_state)
  allocate(parent_rollback_level_one_temperature, source=fine_temperature)
  allocate(parent_rollback_level_two_state, source=level_two_state)
  allocate(parent_rollback_level_two_temperature, &
    source=level_two_temperature)
  initial_i_lower = patch%coarse_i_lower
  initial_i_upper = patch%coarse_i_upper
  initial_j_lower = patch%coarse_j_lower
  initial_j_upper = patch%coarse_j_upper
  parent_level_two_i_lower = level_two_patch%coarse_i_lower
  parent_level_two_i_upper = level_two_patch%coarse_i_upper
  parent_level_two_j_lower = level_two_patch%coarse_j_lower
  parent_level_two_j_upper = level_two_patch%coarse_j_upper
  config%prolongation_method = "invalid"
  call regrid_three_level_reactive_eb_amr_parent_2d( &
    species, config, coarse_state, coarse_temperature, coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, patch, level_two_state, &
    level_two_temperature, level_two_geometry, level_two_patch, changed, ok)
  call require(.not. ok .and. .not. changed .and. &
    all(coarse_state == parent_rollback_root_state) .and. &
    all(coarse_temperature == parent_rollback_root_temperature) .and. &
    all(fine_state == parent_rollback_level_one_state) .and. &
    all(fine_temperature == parent_rollback_level_one_temperature) .and. &
    all(level_two_state == parent_rollback_level_two_state) .and. &
    all(level_two_temperature == parent_rollback_level_two_temperature) .and. &
    patch%coarse_i_lower == initial_i_lower .and. &
    patch%coarse_i_upper == initial_i_upper .and. &
    patch%coarse_j_lower == initial_j_lower .and. &
    patch%coarse_j_upper == initial_j_upper .and. &
    level_two_patch%coarse_i_lower == parent_level_two_i_lower .and. &
    level_two_patch%coarse_i_upper == parent_level_two_i_upper .and. &
    level_two_patch%coarse_j_lower == parent_level_two_j_lower .and. &
    level_two_patch%coarse_j_upper == parent_level_two_j_upper, &
    "three-level parent-regrid rollback")

  config%prolongation_method = "linear"
  config%dynamic_parent_regridding = .true.
  config%regrid_at_initialization = .true.
  call simulate_three_level_reactive_eb_amr_2d( &
    species, reactions, config, coarse_state, coarse_temperature, &
    coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
    level_two_state, level_two_temperature, level_two_geometry, &
    level_two_patch, time, steps, regrids, initial_integrals, &
    final_integrals, minimum_dt, base_density, ok)
  scale = max(1.0_dp, maxval(abs(initial_integrals)))
  call require(ok .and. steps > 0 .and. regrids > 0 .and. &
    patch%is_valid(coarse_geometry, fine_geometry) .and. &
    level_two_patch%is_valid(fine_geometry, level_two_geometry) .and. &
    (patch%coarse_i_lower /= config%coarse_i_lower .or. &
     patch%coarse_i_upper /= config%coarse_i_upper .or. &
     patch%coarse_j_lower /= config%coarse_j_lower .or. &
    patch%coarse_j_upper /= config%coarse_j_upper) .and. &
    maxval(abs(final_integrals - initial_integrals)) <= &
      2.0e-8_dp * scale, &
    "public scheduled three-level parent regrid")

  call write_reactive_eb_amr_three_level_2d_checkpoint( &
    dynamic_three_level_checkpoint_path, species, config, coarse_state, &
    coarse_temperature, coarse_geometry, fine_state, fine_temperature, &
    fine_geometry, patch, level_two_state, level_two_temperature, &
    level_two_geometry, level_two_patch, time, steps, minimum_dt, &
    base_density, ok, regrids, failure_context=checkpoint_failure_context)
  call require(ok .and. len_trim(checkpoint_failure_context) == 0, &
    "fixed dynamic three-level checkpoint write")
  call write_reactive_eb_amr_three_level_2d_checkpoint( &
    selected_dynamic_three_level_checkpoint_path, species, config, &
    coarse_state, coarse_temperature, coarse_geometry, fine_state, &
    fine_temperature, fine_geometry, patch, level_two_state, &
    level_two_temperature, level_two_geometry, level_two_patch, time, &
    steps, minimum_dt, base_density, ok, regrids, &
    bundle_sha256=selected_bundle_sha256, chemistry_integrator="implicit", &
    base_mole_fractions=selected_mole_fractions, &
    initial_integrals=initial_integrals, &
    failure_context=checkpoint_failure_context)
  call require(ok .and. len_trim(checkpoint_failure_context) == 0, &
    "selected dynamic three-level checkpoint write")
  call require_selected_three_level_checkpoint_context( &
    selected_dynamic_three_level_checkpoint_path, selected_bundle_sha256, &
    "implicit", selected_mole_fractions, initial_integrals, &
    expected_schema=5)
  call read_selected_three_level_checkpoint( &
    selected_dynamic_three_level_checkpoint_path, selected_bundle_sha256, &
    "implicit", selected_mole_fractions)
  call require(ok .and. checkpoint_time == time .and. &
    checkpoint_steps == steps .and. checkpoint_regrids == regrids .and. &
    checkpoint_minimum_dt == minimum_dt .and. &
    checkpoint_base_density == base_density .and. &
    checkpoint_patch%is_valid( &
      checkpoint_coarse_geometry, checkpoint_fine_geometry) .and. &
    checkpoint_level_two_patch%is_valid( &
      checkpoint_fine_geometry, checkpoint_level_two_geometry) .and. &
    all(checkpoint_coarse_state == coarse_state) .and. &
    all(checkpoint_fine_state == fine_state) .and. &
    all(checkpoint_level_two_state == level_two_state) .and. &
    allocated(checkpoint_initial_integrals) .and. &
    all(checkpoint_initial_integrals == initial_integrals), &
    "selected dynamic three-level checkpoint round trip")

  allocate(invalid_initial_integrals, source=initial_integrals)
  invalid_initial_integrals(iet + 1) = &
    invalid_initial_integrals(iet + 1) + 1.0e-3_dp * &
      invalid_initial_integrals(irho)
  call write_reactive_eb_amr_three_level_2d_checkpoint( &
    selected_dynamic_three_level_checkpoint_path, species, config, &
    coarse_state, coarse_temperature, coarse_geometry, fine_state, &
    fine_temperature, fine_geometry, patch, level_two_state, &
    level_two_temperature, level_two_geometry, level_two_patch, time, &
    steps, minimum_dt, base_density, ok, regrids, &
    bundle_sha256=selected_bundle_sha256, chemistry_integrator="implicit", &
    base_mole_fractions=selected_mole_fractions, &
    initial_integrals=invalid_initial_integrals, &
    failure_context=checkpoint_failure_context)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "selected checkpoint baseline is invalid") > 0, &
    "selected dynamic three-level baseline closure rejection")
  deallocate(invalid_initial_integrals)
  call read_selected_three_level_checkpoint( &
    selected_dynamic_three_level_checkpoint_path, selected_bundle_sha256, &
    "implicit", selected_mole_fractions)
  call require(ok .and. all(checkpoint_coarse_state == coarse_state) .and. &
    all(checkpoint_fine_state == fine_state) .and. &
    all(checkpoint_level_two_state == level_two_state), &
    "invalid dynamic baseline write preserves existing checkpoint")

  call read_fixed_three_level_checkpoint( &
    selected_dynamic_three_level_checkpoint_path)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "not a fixed-runtime schema") > 0, &
    "fixed reader rejects selected dynamic three-level checkpoint")
  call require_three_level_checkpoint_targets_empty( &
    "selected dynamic checkpoint fixed-reader rollback")
  call read_selected_three_level_checkpoint( &
    dynamic_three_level_checkpoint_path, selected_bundle_sha256, &
    "implicit", selected_mole_fractions)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "lacks selected mechanism context") > 0, &
    "selected reader rejects fixed dynamic three-level checkpoint")
  call require_three_level_checkpoint_targets_empty( &
    "fixed dynamic checkpoint selected-reader rollback")

  call read_selected_three_level_checkpoint( &
    selected_three_level_checkpoint_path, selected_bundle_sha256, &
    "implicit", selected_mole_fractions)
  call require(.not. ok, &
    "selected dynamic reader rejects selected static checkpoint")
  call require_three_level_checkpoint_targets_empty( &
    "selected static checkpoint dynamic-reader rollback")
  config%dynamic_regridding = .false.
  config%dynamic_parent_regridding = .false.
  call read_selected_three_level_checkpoint( &
    selected_dynamic_three_level_checkpoint_path, selected_bundle_sha256, &
    "implicit", selected_mole_fractions)
  call require(.not. ok, &
    "selected static reader rejects selected dynamic checkpoint")
  call require_three_level_checkpoint_targets_empty( &
    "selected dynamic checkpoint static-reader rollback")
  config%dynamic_regridding = .true.
  config%dynamic_parent_regridding = .true.

  call replace_checkpoint_record( &
    selected_dynamic_three_level_checkpoint_path, &
    selected_dynamic_three_level_corrupt_path, "COMPOSITE_BASELINE", &
    "BROKEN_COMPOSITE_BASELINE")
  call read_selected_three_level_checkpoint( &
    selected_dynamic_three_level_corrupt_path, selected_bundle_sha256, &
    "implicit", selected_mole_fractions)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "selected checkpoint baseline is invalid") > 0, &
    "selected dynamic baseline marker rejection")
  call require_three_level_checkpoint_targets_empty( &
    "selected dynamic baseline marker rollback")
  allocate(invalid_initial_integrals, source=initial_integrals)
  invalid_initial_integrals(irho) = -1.0_dp
  write(invalid_baseline_record, '(*(es27.18e3,1x))') &
    invalid_initial_integrals
  call replace_checkpoint_line( &
    selected_dynamic_three_level_checkpoint_path, &
    selected_dynamic_three_level_corrupt_path, 10, &
    trim(invalid_baseline_record))
  call read_selected_three_level_checkpoint( &
    selected_dynamic_three_level_corrupt_path, selected_bundle_sha256, &
    "implicit", selected_mole_fractions)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "selected checkpoint baseline is invalid") > 0, &
    "selected dynamic numeric baseline rejection")
  call require_three_level_checkpoint_targets_empty( &
    "selected dynamic numeric baseline rollback")
  deallocate(invalid_initial_integrals)

  call replace_checkpoint_line( &
    selected_dynamic_three_level_checkpoint_path, &
    selected_dynamic_three_level_corrupt_path, 33, "1 1 2 0 4 4 0")
  call read_selected_three_level_checkpoint( &
    selected_dynamic_three_level_corrupt_path, selected_bundle_sha256, &
    "implicit", selected_mole_fractions)
  call require(.not. ok, "selected dynamic control corruption rejection")
  call require_three_level_checkpoint_targets_empty( &
    "selected dynamic control corruption rollback")
  call replace_checkpoint_line( &
    selected_dynamic_three_level_checkpoint_path, &
    selected_dynamic_three_level_corrupt_path, 35, "0 1 1 1 2")
  call read_selected_three_level_checkpoint( &
    selected_dynamic_three_level_corrupt_path, selected_bundle_sha256, &
    "implicit", selected_mole_fractions)
  call require(.not. ok, "selected dynamic parent patch corruption rejection")
  call require_three_level_checkpoint_targets_empty( &
    "selected dynamic parent patch corruption rollback")
  call replace_checkpoint_line( &
    selected_dynamic_three_level_checkpoint_path, &
    selected_dynamic_three_level_corrupt_path, 39, "BROKEN_FIELD")
  call read_selected_three_level_checkpoint( &
    selected_dynamic_three_level_corrupt_path, selected_bundle_sha256, &
    "implicit", selected_mole_fractions)
  call require(.not. ok, "selected dynamic field corruption rejection")
  call require_three_level_checkpoint_targets_empty( &
    "selected dynamic field corruption rollback")
  call replace_checkpoint_record( &
    selected_dynamic_three_level_checkpoint_path, &
    selected_dynamic_three_level_corrupt_path, "END_CHECKPOINT", &
    "BROKEN_END_CHECKPOINT")
  call read_selected_three_level_checkpoint( &
    selected_dynamic_three_level_corrupt_path, selected_bundle_sha256, &
    "implicit", selected_mole_fractions)
  call require(.not. ok, "selected dynamic terminal marker rejection")
  call require_three_level_checkpoint_targets_empty( &
    "selected dynamic terminal marker rollback")
  call append_checkpoint_record( &
    selected_dynamic_three_level_checkpoint_path, &
    selected_dynamic_three_level_corrupt_path, "TRAILING_CONTENT")
  call read_selected_three_level_checkpoint( &
    selected_dynamic_three_level_corrupt_path, selected_bundle_sha256, &
    "implicit", selected_mole_fractions)
  call require(.not. ok .and. index(checkpoint_failure_context, &
    "trailing content") > 0, &
    "selected dynamic trailing content rejection")
  call require_three_level_checkpoint_targets_empty( &
    "selected dynamic trailing content rollback")

  call delete_checkpoint(three_level_checkpoint_path)
  call delete_checkpoint(selected_three_level_checkpoint_path)
  call delete_checkpoint(dynamic_three_level_checkpoint_path)
  call delete_checkpoint(selected_dynamic_three_level_checkpoint_path)
  call delete_checkpoint(selected_dynamic_three_level_corrupt_path)

  write(*, '(a)') "test_reactive_eb_amr_2d_driver: PASS"

contains

  subroutine read_selected_checkpoint( &
      path, bundle_sha256, chemistry_integrator, base_mole_fractions)
    character(len=*), intent(in) :: path, bundle_sha256
    character(len=*), intent(in) :: chemistry_integrator
    real(dp), intent(in) :: base_mole_fractions(:)

    if (config%dynamic_regridding) then
      call read_reactive_eb_amr_2d_checkpoint( &
        path, species, config, checkpoint_coarse_state, &
        checkpoint_coarse_temperature, checkpoint_coarse_geometry, &
        checkpoint_fine_state, checkpoint_fine_temperature, &
        checkpoint_fine_geometry, checkpoint_patch, checkpoint_fine_active, &
        checkpoint_time, checkpoint_steps, checkpoint_regrids, &
        checkpoint_minimum_dt, checkpoint_base_density, ok, &
        bundle_sha256=bundle_sha256, &
        chemistry_integrator=chemistry_integrator, &
        base_mole_fractions=base_mole_fractions, &
        initial_integrals=checkpoint_initial_integrals, &
        failure_context=checkpoint_failure_context)
    else
      call read_reactive_eb_amr_2d_checkpoint( &
        path, species, config, checkpoint_coarse_state, &
        checkpoint_coarse_temperature, checkpoint_coarse_geometry, &
        checkpoint_fine_state, checkpoint_fine_temperature, &
        checkpoint_fine_geometry, checkpoint_patch, checkpoint_fine_active, &
        checkpoint_time, checkpoint_steps, checkpoint_regrids, &
        checkpoint_minimum_dt, checkpoint_base_density, ok, &
        bundle_sha256=bundle_sha256, &
        chemistry_integrator=chemistry_integrator, &
        base_mole_fractions=base_mole_fractions, &
        failure_context=checkpoint_failure_context)
    end if
  end subroutine read_selected_checkpoint

  subroutine read_selected_three_level_checkpoint( &
      path, bundle_sha256, chemistry_integrator, base_mole_fractions)
    character(len=*), intent(in) :: path, bundle_sha256
    character(len=*), intent(in) :: chemistry_integrator
    real(dp), intent(in) :: base_mole_fractions(:)

    call read_reactive_eb_amr_three_level_2d_checkpoint( &
      path, species, config, checkpoint_coarse_state, &
      checkpoint_coarse_temperature, checkpoint_coarse_geometry, &
      checkpoint_fine_state, checkpoint_fine_temperature, &
      checkpoint_fine_geometry, checkpoint_patch, &
      checkpoint_level_two_state, checkpoint_level_two_temperature, &
      checkpoint_level_two_geometry, checkpoint_level_two_patch, &
      checkpoint_time, checkpoint_steps, checkpoint_minimum_dt, &
      checkpoint_base_density, ok, checkpoint_regrids, &
      bundle_sha256=bundle_sha256, &
      chemistry_integrator=chemistry_integrator, &
      base_mole_fractions=base_mole_fractions, &
      initial_integrals=checkpoint_initial_integrals, &
      failure_context=checkpoint_failure_context)
  end subroutine read_selected_three_level_checkpoint

  subroutine read_fixed_three_level_checkpoint(path)
    character(len=*), intent(in) :: path

    if (allocated(checkpoint_initial_integrals)) &
      deallocate(checkpoint_initial_integrals)
    call read_reactive_eb_amr_three_level_2d_checkpoint( &
      path, species, config, checkpoint_coarse_state, &
      checkpoint_coarse_temperature, checkpoint_coarse_geometry, &
      checkpoint_fine_state, checkpoint_fine_temperature, &
      checkpoint_fine_geometry, checkpoint_patch, &
      checkpoint_level_two_state, checkpoint_level_two_temperature, &
      checkpoint_level_two_geometry, checkpoint_level_two_patch, &
      checkpoint_time, checkpoint_steps, checkpoint_minimum_dt, &
      checkpoint_base_density, ok, checkpoint_regrids, &
      failure_context=checkpoint_failure_context)
  end subroutine read_fixed_three_level_checkpoint

  subroutine require_selected_three_level_checkpoint_context( &
      path, expected_sha256, expected_integrator, expected_composition, &
      expected_initial_integrals, expected_schema)
    character(len=*), intent(in) :: path, expected_sha256
    character(len=*), intent(in) :: expected_integrator
    real(dp), intent(in) :: expected_composition(:)
    real(dp), intent(in) :: expected_initial_integrals(:)
    integer, intent(in), optional :: expected_schema

    character(len=1024) :: magic, marker, stored_sha256, stored_integrator
    character(len=1024) :: expected_magic
    character(len=1024) :: baseline_marker
    real(dp), allocatable :: stored_composition(:), stored_baseline(:)
    integer :: unit, status, header(3), stored_size, stored_baseline_size, schema

    schema = 4
    if (present(expected_schema)) schema = expected_schema
    expected_magic = "PELEF_REACTIVE_EB_AMR_THREE_LEVEL_2D_CHECKPOINT"
    if (schema == 5) expected_magic = &
      "PELEF_REACTIVE_EB_AMR_DYNAMIC_THREE_LEVEL_2D_CHECKPOINT"

    open(newunit=unit, file=path, status="old", action="read", iostat=status)
    call require(status == 0, &
      "open selected three-level checkpoint context probe")
    read(unit, '(a)', iostat=status) magic
    call require(status == 0 .and. trim(magic) == trim(expected_magic), &
      "selected three-level checkpoint magic")
    read(unit, *, iostat=status) header
    call require(status == 0 .and. header(1) == schema .and. &
      header(2) == size(expected_composition) .and. &
      header(3) == size(expected_initial_integrals), &
      "selected three-level checkpoint schema")
    read(unit, '(a)', iostat=status) marker
    call require(status == 0 .and. trim(marker) == "SELECTED_CONTEXT", &
      "selected three-level checkpoint context marker")
    read(unit, '(a)', iostat=status) stored_sha256
    call require(status == 0 .and. &
      trim(stored_sha256) == trim(expected_sha256), &
      "selected three-level checkpoint bundle SHA")
    read(unit, '(a)', iostat=status) stored_integrator
    call require(status == 0 .and. &
      trim(stored_integrator) == trim(expected_integrator), &
      "selected three-level checkpoint integrator")
    read(unit, *, iostat=status) stored_size
    call require(status == 0 .and. &
      stored_size == size(expected_composition), &
      "selected three-level checkpoint composition size")
    allocate(stored_composition(stored_size))
    read(unit, *, iostat=status) stored_composition
    call require(status == 0 .and. &
      all(stored_composition == expected_composition), &
      "selected three-level checkpoint composition")
    read(unit, '(a)', iostat=status) baseline_marker
    call require(status == 0 .and. &
      trim(baseline_marker) == "COMPOSITE_BASELINE", &
      "selected three-level checkpoint baseline marker")
    read(unit, *, iostat=status) stored_baseline_size
    call require(status == 0 .and. &
      stored_baseline_size == size(expected_initial_integrals), &
      "selected three-level checkpoint baseline size")
    allocate(stored_baseline(stored_baseline_size))
    read(unit, *, iostat=status) stored_baseline
    call require(status == 0 .and. &
      all(stored_baseline == expected_initial_integrals), &
      "selected three-level checkpoint baseline")
    close(unit)
  end subroutine require_selected_three_level_checkpoint_context

  subroutine require_three_level_checkpoint_targets_empty(label)
    character(len=*), intent(in) :: label

    call require( &
      .not. allocated(checkpoint_coarse_state) .and. &
      .not. allocated(checkpoint_coarse_temperature) .and. &
      .not. allocated(checkpoint_fine_state) .and. &
      .not. allocated(checkpoint_fine_temperature) .and. &
      .not. allocated(checkpoint_level_two_state) .and. &
      .not. allocated(checkpoint_level_two_temperature) .and. &
      .not. allocated(checkpoint_initial_integrals) .and. &
      .not. checkpoint_coarse_geometry%is_valid() .and. &
      .not. checkpoint_fine_geometry%is_valid() .and. &
      .not. checkpoint_level_two_geometry%is_valid() .and. &
      checkpoint_time == 0.0_dp .and. checkpoint_steps == 0 .and. &
      checkpoint_regrids == 0 .and. &
      checkpoint_minimum_dt == 0.0_dp .and. &
      checkpoint_base_density == 0.0_dp, label)
  end subroutine require_three_level_checkpoint_targets_empty

  subroutine require_selected_patch_set_checkpoint_context( &
      path, expected_sha256, expected_integrator, expected_composition, &
      expected_initial_integrals, expected_patch_count)
    character(len=*), intent(in) :: path, expected_sha256
    character(len=*), intent(in) :: expected_integrator
    real(dp), intent(in) :: expected_composition(:)
    real(dp), intent(in) :: expected_initial_integrals(:)
    integer, intent(in) :: expected_patch_count

    character(len=1024) :: magic, marker, stored_sha256, stored_integrator
    real(dp), allocatable :: stored_composition(:), stored_integrals(:)
    integer :: unit, status, header(4), stored_size

    open(newunit=unit, file=path, status="old", action="read", &
      iostat=status)
    call require(status == 0, &
      "open selected multipatch checkpoint context probe")
    read(unit, '(a)', iostat=status) magic
    call require(status == 0 .and. trim(magic) == &
      "PELEF_REACTIVE_EB_AMR_PATCH_SET_2D_CHECKPOINT", &
      "selected multipatch checkpoint context magic")
    read(unit, *, iostat=status) header
    call require(status == 0 .and. header(1) == 4 .and. &
      header(2) == size(expected_composition) .and. &
      header(4) == expected_patch_count, &
      "selected multipatch checkpoint context schema")
    read(unit, '(a)', iostat=status) marker
    call require(status == 0 .and. trim(marker) == "SELECTED_CONTEXT", &
      "selected multipatch checkpoint context marker")
    read(unit, '(a)', iostat=status) stored_sha256
    call require(status == 0 .and. &
      trim(stored_sha256) == trim(expected_sha256), &
      "selected multipatch checkpoint stored bundle SHA")
    read(unit, '(a)', iostat=status) stored_integrator
    call require(status == 0 .and. &
      trim(stored_integrator) == trim(expected_integrator), &
      "selected multipatch checkpoint stored integrator")
    read(unit, *, iostat=status) stored_size
    call require(status == 0 .and. &
      stored_size == size(expected_composition), &
      "selected multipatch checkpoint stored composition size")
    allocate(stored_composition(stored_size))
    read(unit, *, iostat=status) stored_composition
    call require(status == 0 .and. &
      all(stored_composition == expected_composition), &
      "selected multipatch checkpoint stored composition")
    read(unit, '(a)', iostat=status) marker
    call require(status == 0 .and. trim(marker) == "COMPOSITE_BASELINE", &
      "selected multipatch checkpoint baseline marker")
    read(unit, *, iostat=status) stored_size
    call require(status == 0 .and. &
      stored_size == size(expected_initial_integrals), &
      "selected multipatch checkpoint baseline size")
    allocate(stored_integrals(stored_size))
    read(unit, *, iostat=status) stored_integrals
    call require(status == 0 .and. &
      all(stored_integrals == expected_initial_integrals), &
      "selected multipatch checkpoint baseline values")
    close(unit)
  end subroutine require_selected_patch_set_checkpoint_context

  subroutine require_patch_set_checkpoint_targets_empty(label)
    character(len=*), intent(in) :: label

    call require( &
      .not. allocated(checkpoint_coarse_state) .and. &
      .not. allocated(checkpoint_coarse_temperature) .and. &
      .not. allocated(checkpoint_initial_integrals) .and. &
      .not. checkpoint_coarse_geometry%is_valid() .and. &
      .not. allocated(checkpoint_multipatch_set%children) .and. &
      checkpoint_time == 0.0_dp .and. checkpoint_steps == 0 .and. &
      checkpoint_regrids == 0 .and. &
      checkpoint_minimum_dt == 0.0_dp .and. &
      checkpoint_base_density == 0.0_dp, label)
  end subroutine require_patch_set_checkpoint_targets_empty

  subroutine require_selected_checkpoint_context( &
      path, expected_sha256, expected_integrator, expected_composition, &
      expected_schema, expected_initial_integrals)
    character(len=*), intent(in) :: path, expected_sha256
    character(len=*), intent(in) :: expected_integrator
    real(dp), intent(in) :: expected_composition(:)
    integer, intent(in), optional :: expected_schema
    real(dp), intent(in), optional :: expected_initial_integrals(:)

    character(len=1024) :: magic, marker, stored_sha256, stored_integrator
    character(len=1024) :: baseline_marker
    real(dp), allocatable :: stored_composition(:), stored_initial_integrals(:)
    integer :: unit, status, header(4), stored_size, schema
    integer :: stored_baseline_size

    schema = 4
    if (present(expected_schema)) schema = expected_schema

    open(newunit=unit, file=path, status="old", action="read", &
      iostat=status)
    call require(status == 0, "open selected AMR checkpoint context probe")
    read(unit, '(a)', iostat=status) magic
    call require(status == 0 .and. &
      trim(magic) == "PELEF_REACTIVE_EB_AMR_2D_CHECKPOINT", &
      "selected AMR checkpoint context magic")
    read(unit, *, iostat=status) header
    call require(status == 0 .and. header(1) == schema .and. &
      header(2) == size(expected_composition) .and. header(4) == 1, &
      "selected AMR checkpoint context schema")
    read(unit, '(a)', iostat=status) marker
    call require(status == 0 .and. trim(marker) == "SELECTED_CONTEXT", &
      "selected AMR checkpoint context marker")
    read(unit, '(a)', iostat=status) stored_sha256
    call require(status == 0 .and. &
      trim(stored_sha256) == trim(expected_sha256), &
      "selected AMR checkpoint stored bundle SHA")
    read(unit, '(a)', iostat=status) stored_integrator
    call require(status == 0 .and. &
      trim(stored_integrator) == trim(expected_integrator), &
      "selected AMR checkpoint stored integrator")
    read(unit, *, iostat=status) stored_size
    call require(status == 0 .and. &
      stored_size == size(expected_composition), &
      "selected AMR checkpoint stored composition size")
    allocate(stored_composition(stored_size))
    read(unit, *, iostat=status) stored_composition
    call require(status == 0 .and. &
      all(stored_composition == expected_composition), &
      "selected AMR checkpoint stored composition")
    if (schema == 5) then
      call require(present(expected_initial_integrals), &
        "selected dynamic checkpoint expected baseline")
      read(unit, '(a)', iostat=status) baseline_marker
      call require(status == 0 .and. &
        trim(baseline_marker) == "DYNAMIC_BASELINE", &
        "selected dynamic checkpoint baseline marker")
      read(unit, *, iostat=status) stored_baseline_size
      call require(status == 0 .and. &
        stored_baseline_size == size(expected_initial_integrals), &
        "selected dynamic checkpoint baseline size")
      allocate(stored_initial_integrals(stored_baseline_size))
      read(unit, *, iostat=status) stored_initial_integrals
      call require(status == 0 .and. &
        all(stored_initial_integrals == expected_initial_integrals), &
        "selected dynamic checkpoint baseline values")
    else
      call require(.not. present(expected_initial_integrals), &
        "selected static checkpoint has no baseline record")
    end if
    close(unit)
  end subroutine require_selected_checkpoint_context

  subroutine require_checkpoint_targets_empty(label)
    character(len=*), intent(in) :: label

    call require( &
      .not. allocated(checkpoint_coarse_state) .and. &
      .not. allocated(checkpoint_coarse_temperature) .and. &
      .not. allocated(checkpoint_fine_state) .and. &
      .not. allocated(checkpoint_fine_temperature) .and. &
      .not. allocated(checkpoint_initial_integrals) .and. &
      .not. checkpoint_coarse_geometry%is_valid() .and. &
      .not. checkpoint_fine_geometry%is_valid() .and. &
      .not. checkpoint_fine_active .and. checkpoint_time == 0.0_dp .and. &
      checkpoint_steps == 0 .and. checkpoint_regrids == 0 .and. &
      checkpoint_minimum_dt == 0.0_dp .and. &
      checkpoint_base_density == 0.0_dp, label)
  end subroutine require_checkpoint_targets_empty

  subroutine require_file_absent(path, label)
    character(len=*), intent(in) :: path, label

    logical :: exists

    inquire(file=path, exist=exists)
    call require(.not. exists, label)
  end subroutine require_file_absent

  subroutine replace_checkpoint_line( &
      source_path, destination_path, line_number, replacement)
    character(len=*), intent(in) :: source_path, destination_path
    integer, intent(in) :: line_number
    character(len=*), intent(in) :: replacement

    character(len=8192) :: line
    integer :: input_unit, output_unit, status, current_line

    open(newunit=input_unit, file=source_path, status="old", action="read", &
      iostat=status)
    if (status /= 0) error stop "Could not open checkpoint mutation source"
    open(newunit=output_unit, file=destination_path, status="replace", &
      action="write", iostat=status)
    if (status /= 0) error stop "Could not open checkpoint mutation target"
    current_line = 0
    do
      read(input_unit, '(a)', iostat=status) line
      if (status < 0) exit
      if (status > 0) error stop "Could not read checkpoint mutation source"
      current_line = current_line + 1
      if (current_line == line_number) then
        write(output_unit, '(a)', iostat=status) trim(replacement)
      else
        write(output_unit, '(a)', iostat=status) trim(line)
      end if
      if (status /= 0) error stop "Could not write checkpoint mutation target"
    end do
    close(input_unit)
    close(output_unit)
    call require(current_line >= line_number, &
      "checkpoint mutation line is present")
  end subroutine replace_checkpoint_line

  subroutine replace_checkpoint_record( &
      source_path, destination_path, target, replacement)
    character(len=*), intent(in) :: source_path, destination_path
    character(len=*), intent(in) :: target, replacement

    character(len=8192) :: line
    integer :: input_unit, output_unit, status
    logical :: replaced

    open(newunit=input_unit, file=source_path, status="old", action="read", &
      iostat=status)
    if (status /= 0) error stop "Could not open checkpoint mutation source"
    open(newunit=output_unit, file=destination_path, status="replace", &
      action="write", iostat=status)
    if (status /= 0) error stop "Could not open checkpoint mutation target"
    replaced = .false.
    do
      read(input_unit, '(a)', iostat=status) line
      if (status < 0) exit
      if (status > 0) error stop "Could not read checkpoint mutation source"
      if (trim(line) == trim(target)) then
        write(output_unit, '(a)', iostat=status) trim(replacement)
        replaced = .true.
      else
        write(output_unit, '(a)', iostat=status) trim(line)
      end if
      if (status /= 0) error stop "Could not write checkpoint mutation target"
    end do
    close(input_unit)
    close(output_unit)
    call require(replaced, "checkpoint mutation record is present")
  end subroutine replace_checkpoint_record

  subroutine append_checkpoint_record( &
      source_path, destination_path, appended_record)
    character(len=*), intent(in) :: source_path, destination_path
    character(len=*), intent(in) :: appended_record

    character(len=8192) :: line
    integer :: input_unit, output_unit, status

    open(newunit=input_unit, file=source_path, status="old", action="read", &
      iostat=status)
    if (status /= 0) error stop "Could not open checkpoint append source"
    open(newunit=output_unit, file=destination_path, status="replace", &
      action="write", iostat=status)
    if (status /= 0) error stop "Could not open checkpoint append target"
    do
      read(input_unit, '(a)', iostat=status) line
      if (status < 0) exit
      if (status > 0) error stop "Could not read checkpoint append source"
      write(output_unit, '(a)', iostat=status) trim(line)
      if (status /= 0) error stop "Could not write checkpoint append target"
    end do
    write(output_unit, '(a)', iostat=status) trim(appended_record)
    if (status /= 0) error stop "Could not append checkpoint record"
    close(input_unit)
    close(output_unit)
  end subroutine append_checkpoint_record

  subroutine make_prefix_checkpoint( &
      source_path, destination_path, retained_lines)
    character(len=*), intent(in) :: source_path, destination_path
    integer, intent(in) :: retained_lines

    character(len=8192) :: line
    integer :: input_unit, output_unit, status, line_number

    open(newunit=input_unit, file=source_path, status="old", action="read", &
      iostat=status)
    if (status /= 0) error stop "Could not open checkpoint prefix source"
    open(newunit=output_unit, file=destination_path, status="replace", &
      action="write", iostat=status)
    if (status /= 0) error stop "Could not open checkpoint prefix target"
    do line_number = 1, retained_lines
      read(input_unit, '(a)', iostat=status) line
      if (status /= 0) error stop "Could not read checkpoint prefix source"
      write(output_unit, '(a)', iostat=status) trim(line)
      if (status /= 0) error stop "Could not write checkpoint prefix target"
    end do
    close(input_unit)
    close(output_unit)
  end subroutine make_prefix_checkpoint

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
