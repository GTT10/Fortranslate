program test_reactive_eb_3d_checkpoint
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport
  use reactive_1d_mod, only: reactive_nvar
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, build_axis_plane_eb_geometry_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: &
    initialize_reactive_eb_density_sheet_3d, reactive_eb_integrals_3d, &
    reactive_eb_element_integrals_3d
  use reactive_eb_3d_checkpoint_mod, only: &
    write_reactive_eb_3d_checkpoint, read_reactive_eb_3d_checkpoint
  implicit none

  character(len=*), parameter :: checkpoint_path = &
    "reactive_eb_3d_checkpoint_test.chk"
  character(len=*), parameter :: truncated_path = &
    "reactive_eb_3d_checkpoint_truncated.chk"
  character(len=*), parameter :: early_truncated_path = &
    "reactive_eb_3d_checkpoint_early_truncated.chk"
  character(len=*), parameter :: bad_magic_path = &
    "reactive_eb_3d_checkpoint_bad_magic.chk"
  character(len=*), parameter :: bad_schema_path = &
    "reactive_eb_3d_checkpoint_bad_schema.chk"
  character(len=*), parameter :: bad_temperature_path = &
    "reactive_eb_3d_checkpoint_bad_temperature.chk"
  character(len=*), parameter :: selected_checkpoint_path = &
    "reactive_eb_3d_checkpoint_selected.chk"
  character(len=*), parameter :: selected_bad_context_path = &
    "reactive_eb_3d_checkpoint_selected_bad_context.chk"
  character(len=*), parameter :: selected_truncated_path = &
    "reactive_eb_3d_checkpoint_selected_truncated.chk"
  character(len=*), parameter :: selected_invalid_write_path = &
    "reactive_eb_3d_checkpoint_selected_invalid_write.chk"
  character(len=*), parameter :: selected_bundle_sha256 = &
    "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef"
  character(len=*), parameter :: changed_bundle_sha256 = &
    "f123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef"
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:), changed_reactions(:)
  type(gas_transport_species), allocatable :: transport(:), changed_transport(:)
  type(reactive_eb_3d_config) :: write_config, restart_config, changed_config
  type(eb_geometry_3d) :: geometry, read_geometry, saved_geometry
  real(dp), allocatable :: state(:, :, :, :), temperature(:, :, :)
  real(dp), allocatable :: read_state(:, :, :, :), read_temperature(:, :, :)
  real(dp), allocatable :: saved_state(:, :, :, :), saved_temperature(:, :, :)
  real(dp), allocatable :: wrong_state(:, :, :, :), wrong_temperature(:, :, :)
  real(dp), allocatable :: initial_integrals(:), initial_l1_integrals(:)
  real(dp), allocatable :: read_integrals(:), read_l1_integrals(:)
  real(dp), allocatable :: saved_integrals(:), saved_l1_integrals(:)
  real(dp), allocatable :: selected_mole_fractions(:)
  real(dp), allocatable :: changed_mole_fractions(:)
  real(dp) :: initial_elements(3), read_elements(3), saved_elements(3)
  real(dp) :: time, minimum_dt, maximum_diffusivity, minimum_theta
  real(dp) :: read_time, read_minimum_dt, read_maximum_diffusivity
  real(dp) :: read_minimum_theta, saved_time, saved_minimum_dt
  real(dp) :: saved_maximum_diffusivity, saved_minimum_theta
  integer :: steps, read_steps, saved_steps, nvar
  logical :: ok
  character(len=256) :: message

  call remove_file(checkpoint_path)
  call remove_file(truncated_path)
  call remove_file(early_truncated_path)
  call remove_file(bad_magic_path)
  call remove_file(bad_schema_path)
  call remove_file(bad_temperature_path)
  call remove_file(selected_checkpoint_path)
  call remove_file(selected_bad_context_path)
  call remove_file(selected_truncated_path)
  call remove_file(selected_invalid_write_path)
  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "load checkpoint thermodynamics")
  call load_h2o2_elementary_mechanism(reactions, ok)
  call require(ok, "load checkpoint mechanism")
  call load_h2o2_elementary_transport(transport, ok)
  call require(ok, "load checkpoint transport")

  write_config = reactive_eb_3d_config()
  write_config%nx = 6
  write_config%ny = 4
  write_config%nz = 3
  write_config%x_upper = 6.0e-5_dp
  write_config%y_upper = 4.0e-5_dp
  write_config%z_upper = 3.0e-5_dp
  write_config%plane_position = 1.95e-5_dp
  write_config%final_time = 1.0e-8_dp
  write_config%maximum_steps = 10
  write_config%chemistry_enabled = .true.
  write_config%transport_enabled = .true.
  write_config%checkpoint_interval_steps = 1
  write_config%stop_after_checkpoint = .true.
  write_config%checkpoint_file = checkpoint_path
  write_config%output_file = "checkpoint_stopped.csv"
  call build_axis_plane_eb_geometry_3d( &
    write_config%nx, write_config%ny, write_config%nz, &
    write_config%x_lower, write_config%x_upper, &
    write_config%y_lower, write_config%y_upper, &
    write_config%z_lower, write_config%z_upper, &
    write_config%plane_axis, write_config%plane_position, geometry, ok)
  call require(ok, "build checkpoint geometry")

  nvar = reactive_nvar(size(species))
  allocate(state(nvar, geometry%nx, geometry%ny, geometry%nz))
  allocate(temperature(geometry%nx, geometry%ny, geometry%nz))
  allocate(read_state, mold=state)
  allocate(read_temperature, mold=temperature)
  allocate(saved_state, mold=state)
  allocate(saved_temperature, mold=temperature)
  allocate(initial_integrals(nvar), initial_l1_integrals(nvar))
  allocate(read_integrals(nvar), read_l1_integrals(nvar))
  allocate(saved_integrals(nvar), saved_l1_integrals(nvar))
  call initialize_reactive_eb_density_sheet_3d( &
    species, write_config, geometry, state, temperature, ok)
  call require(ok, "initialize checkpoint state")
  call reactive_eb_integrals_3d( &
    state, geometry, initial_integrals, ok, initial_l1_integrals)
  call require(ok, "integrate checkpoint state")
  call reactive_eb_element_integrals_3d( &
    species, state, geometry, initial_elements, ok)
  call require(ok, "integrate checkpoint elements")

  time = 2.0e-9_dp
  steps = 1
  minimum_dt = 2.0e-9_dp
  maximum_diffusivity = 9.7180419287635654e-4_dp
  minimum_theta = 0.875_dp
  call write_reactive_eb_3d_checkpoint( &
    checkpoint_path, species, reactions, transport, write_config, geometry, &
    state, temperature, time, steps, initial_integrals, &
    initial_l1_integrals, initial_elements, minimum_dt, &
    maximum_diffusivity, minimum_theta, ok, message)
  call require(ok, "write checkpoint")
  call require_checkpoint_schema(checkpoint_path, 1, "SPECIES")

  restart_config = write_config
  restart_config%maximum_steps = 20
  restart_config%final_time = 2.0e-8_dp
  restart_config%checkpoint_interval_steps = 0
  restart_config%stop_after_checkpoint = .false.
  restart_config%checkpoint_file = "different_schedule.chk"
  restart_config%restart_file = checkpoint_path
  restart_config%output_file = "checkpoint_restart.csv"
  read_geometry = geometry
  read_state = -7.0_dp
  read_temperature = -8.0_dp
  read_integrals = -9.0_dp
  read_l1_integrals = -10.0_dp
  read_elements = -11.0_dp
  read_time = -12.0_dp
  read_steps = -13
  read_minimum_dt = -14.0_dp
  read_maximum_diffusivity = -15.0_dp
  read_minimum_theta = -16.0_dp
  call read_checkpoint(restart_config, reactions, transport, checkpoint_path)
  call require(ok, "read checkpoint with changed continuation policy")
  call require(same_geometry(read_geometry, geometry), &
    "checkpoint geometry round trip")
  call require(all(read_state == state), "checkpoint state round trip")
  call require(all(read_temperature == temperature), &
    "checkpoint temperature round trip")
  call require(all(read_integrals == initial_integrals), &
    "checkpoint integral round trip")
  call require(all(read_l1_integrals == initial_l1_integrals), &
    "checkpoint L1 round trip")
  call require(all(read_elements == initial_elements), &
    "checkpoint element round trip")
  call require(read_time == time .and. read_steps == steps .and. &
    read_minimum_dt == minimum_dt .and. &
    read_maximum_diffusivity == maximum_diffusivity .and. &
    read_minimum_theta == minimum_theta, "checkpoint metadata round trip")

  changed_config = restart_config
  changed_config%cfl = 0.45_dp
  call save_targets()
  call read_checkpoint(changed_config, reactions, transport, checkpoint_path)
  call require(.not. ok, "changed checkpoint configuration rejection")
  call require_targets_unchanged("changed configuration rollback")

  allocate(changed_reactions, source=reactions)
  changed_reactions(1)%forward_rate%pre_exponential = &
    1.01_dp * changed_reactions(1)%forward_rate%pre_exponential
  call save_targets()
  call read_checkpoint( &
    restart_config, changed_reactions, transport, checkpoint_path)
  call require(.not. ok, "changed checkpoint mechanism rejection")
  call require_targets_unchanged("changed mechanism rollback")

  allocate(changed_transport, source=transport)
  changed_transport(1)%diameter = 1.01_dp * changed_transport(1)%diameter
  call save_targets()
  call read_checkpoint( &
    restart_config, reactions, changed_transport, checkpoint_path)
  call require(.not. ok, "changed checkpoint transport rejection")
  call require_targets_unchanged("changed transport rollback")

  call make_truncated_checkpoint(checkpoint_path, truncated_path)
  call save_targets()
  call read_checkpoint(restart_config, reactions, transport, truncated_path)
  call require(.not. ok, "truncated checkpoint rejection")
  call require_targets_unchanged("truncated checkpoint rollback")

  call make_prefix_checkpoint(checkpoint_path, early_truncated_path, 0)
  call save_targets()
  call read_checkpoint( &
    restart_config, reactions, transport, early_truncated_path)
  call require(.not. ok, "empty checkpoint rejection")
  call require_targets_unchanged("empty checkpoint rollback")

  call make_prefix_checkpoint(checkpoint_path, early_truncated_path, 1)
  call save_targets()
  call read_checkpoint( &
    restart_config, reactions, transport, early_truncated_path)
  call require(.not. ok, "magic-only checkpoint rejection")
  call require_targets_unchanged("magic-only checkpoint rollback")

  call make_prefix_checkpoint(checkpoint_path, early_truncated_path, 3)
  call save_targets()
  call read_checkpoint( &
    restart_config, reactions, transport, early_truncated_path)
  call require(.not. ok, "species-header-only checkpoint rejection")
  call require_targets_unchanged("species-header-only checkpoint rollback")

  call replace_checkpoint_line( &
    checkpoint_path, bad_magic_path, 1, "NOT_A_PELEF_CHECKPOINT")
  call save_targets()
  call read_checkpoint(restart_config, reactions, transport, bad_magic_path)
  call require(.not. ok, "bad checkpoint magic rejection")
  call require_targets_unchanged("bad checkpoint magic rollback")

  call replace_checkpoint_line( &
    checkpoint_path, bad_schema_path, 2, "999 0 0 0 0")
  call save_targets()
  call read_checkpoint(restart_config, reactions, transport, bad_schema_path)
  call require(.not. ok, "bad checkpoint schema rejection")
  call require_targets_unchanged("bad checkpoint schema rollback")

  call make_bad_temperature_checkpoint( &
    checkpoint_path, bad_temperature_path, nvar)
  call save_targets()
  call read_checkpoint( &
    restart_config, reactions, transport, bad_temperature_path)
  call require(.not. ok, "inconsistent checkpoint temperature rejection")
  call require_targets_unchanged("inconsistent temperature rollback")

  allocate(wrong_state(nvar, geometry%nx - 1, geometry%ny, geometry%nz))
  allocate(wrong_temperature(geometry%nx - 1, geometry%ny, geometry%nz))
  wrong_state = 17.0_dp
  wrong_temperature = 18.0_dp
  call read_reactive_eb_3d_checkpoint( &
    checkpoint_path, species, reactions, transport, restart_config, &
    read_geometry, wrong_state, wrong_temperature, read_time, read_steps, &
    read_integrals, read_l1_integrals, read_elements, read_minimum_dt, &
    read_maximum_diffusivity, read_minimum_theta, ok, message)
  call require(.not. ok .and. all(wrong_state == 17.0_dp) .and. &
    all(wrong_temperature == 18.0_dp), "wrong restart shape rollback")

  allocate(selected_mole_fractions(size(species)))
  allocate(changed_mole_fractions(size(species)))
  selected_mole_fractions = 0.0_dp
  selected_mole_fractions(1) = 0.25_dp
  selected_mole_fractions(4) = 0.25_dp
  selected_mole_fractions(7) = 0.50_dp
  call write_reactive_eb_3d_checkpoint( &
    selected_checkpoint_path, species, reactions, transport, write_config, &
    geometry, state, temperature, time, steps, initial_integrals, &
    initial_l1_integrals, initial_elements, minimum_dt, &
    maximum_diffusivity, minimum_theta, ok, message, &
    bundle_sha256=selected_bundle_sha256, &
    chemistry_integrator="explicit", &
    base_mole_fractions=selected_mole_fractions)
  call require(ok, "write selected checkpoint")
  call require_selected_checkpoint_context( &
    selected_checkpoint_path, selected_bundle_sha256, "explicit", &
    selected_mole_fractions)

  call write_reactive_eb_3d_checkpoint( &
    selected_invalid_write_path, species, reactions, transport, &
    write_config, geometry, state, temperature, time, steps, &
    initial_integrals, initial_l1_integrals, initial_elements, minimum_dt, &
    maximum_diffusivity, minimum_theta, ok, message, &
    bundle_sha256=selected_bundle_sha256)
  call require(.not. ok .and. &
    index(message, "selected checkpoint context is incomplete") > 0, &
    "incomplete selected write context rejection")
  call require_file_absent( &
    selected_invalid_write_path, "incomplete selected write creates no file")

  call save_targets()
  call read_checkpoint( &
    restart_config, reactions, transport, selected_checkpoint_path)
  call require(.not. ok .and. &
    index(message, "not a fixed-runtime schema") > 0, &
    "fixed reader rejects selected checkpoint")
  call require_targets_unchanged("selected checkpoint fixed-reader rollback")

  call save_targets()
  call read_selected_checkpoint( &
    restart_config, reactions, transport, checkpoint_path, &
    selected_bundle_sha256, "explicit", selected_mole_fractions)
  call require(.not. ok .and. &
    index(message, "lacks selected mechanism context") > 0, &
    "selected reader rejects fixed checkpoint")
  call require_targets_unchanged("fixed checkpoint selected-reader rollback")

  call read_selected_checkpoint( &
    restart_config, reactions, transport, selected_checkpoint_path, &
    selected_bundle_sha256, "explicit", selected_mole_fractions)
  call require(ok, "read selected checkpoint")
  call require(same_geometry(read_geometry, geometry), &
    "selected checkpoint geometry round trip")
  call require(all(read_state == state), &
    "selected checkpoint state round trip")
  call require(all(read_temperature == temperature), &
    "selected checkpoint temperature round trip")
  call require(all(read_integrals == initial_integrals) .and. &
    all(read_l1_integrals == initial_l1_integrals) .and. &
    all(read_elements == initial_elements), &
    "selected checkpoint diagnostics round trip")
  call require(read_time == time .and. read_steps == steps .and. &
    read_minimum_dt == minimum_dt .and. &
    read_maximum_diffusivity == maximum_diffusivity .and. &
    read_minimum_theta == minimum_theta, &
    "selected checkpoint metadata round trip")

  call save_targets()
  call read_selected_checkpoint( &
    restart_config, reactions, transport, selected_checkpoint_path, &
    changed_bundle_sha256, "explicit", selected_mole_fractions)
  call require(.not. ok .and. &
    index(message, "bundle SHA-256 mismatch") > 0, &
    "selected checkpoint bundle mismatch rejection")
  call require_targets_unchanged("selected bundle mismatch rollback")

  call save_targets()
  call read_selected_checkpoint( &
    restart_config, reactions, transport, selected_checkpoint_path, &
    selected_bundle_sha256, "implicit", selected_mole_fractions)
  call require(.not. ok .and. &
    index(message, "chemistry integrator mismatch") > 0, &
    "selected checkpoint integrator mismatch rejection")
  call require_targets_unchanged("selected integrator mismatch rollback")

  changed_mole_fractions = selected_mole_fractions
  changed_mole_fractions(1) = 0.24_dp
  changed_mole_fractions(4) = 0.26_dp
  call save_targets()
  call read_selected_checkpoint( &
    restart_config, reactions, transport, selected_checkpoint_path, &
    selected_bundle_sha256, "explicit", changed_mole_fractions)
  call require(.not. ok .and. &
    index(message, "composition mismatch") > 0, &
    "selected checkpoint composition mismatch rejection")
  call require_targets_unchanged("selected composition mismatch rollback")

  call replace_checkpoint_line( &
    selected_checkpoint_path, selected_bad_context_path, 3, &
    "BROKEN_SELECTED_CONTEXT")
  call save_targets()
  call read_selected_checkpoint( &
    restart_config, reactions, transport, selected_bad_context_path, &
    selected_bundle_sha256, "explicit", selected_mole_fractions)
  call require(.not. ok .and. &
    index(message, "selected checkpoint context is invalid") > 0, &
    "selected checkpoint marker rejection")
  call require_targets_unchanged("selected marker rollback")

  call replace_checkpoint_line( &
    selected_checkpoint_path, selected_bad_context_path, 4, "abc")
  call save_targets()
  call read_selected_checkpoint( &
    restart_config, reactions, transport, selected_bad_context_path, &
    selected_bundle_sha256, "explicit", selected_mole_fractions)
  call require(.not. ok .and. &
    index(message, "selected checkpoint context is invalid") > 0, &
    "selected checkpoint malformed SHA rejection")
  call require_targets_unchanged("selected malformed SHA rollback")

  call replace_checkpoint_line( &
    selected_checkpoint_path, selected_bad_context_path, 5, "auto")
  call save_targets()
  call read_selected_checkpoint( &
    restart_config, reactions, transport, selected_bad_context_path, &
    selected_bundle_sha256, "explicit", selected_mole_fractions)
  call require(.not. ok .and. &
    index(message, "selected checkpoint context is invalid") > 0, &
    "selected checkpoint invalid integrator rejection")
  call require_targets_unchanged("selected invalid integrator rollback")

  call replace_checkpoint_line( &
    selected_checkpoint_path, selected_bad_context_path, 6, "999")
  call save_targets()
  call read_selected_checkpoint( &
    restart_config, reactions, transport, selected_bad_context_path, &
    selected_bundle_sha256, "explicit", selected_mole_fractions)
  call require(.not. ok .and. &
    index(message, "selected checkpoint context is invalid") > 0, &
    "selected checkpoint composition size rejection")
  call require_targets_unchanged("selected composition size rollback")

  call replace_checkpoint_line( &
    selected_checkpoint_path, selected_bad_context_path, 7, &
    "-1.0 1.0 0.0 0.0 0.0 0.0 1.0")
  call save_targets()
  call read_selected_checkpoint( &
    restart_config, reactions, transport, selected_bad_context_path, &
    selected_bundle_sha256, "explicit", selected_mole_fractions)
  call require(.not. ok .and. &
    index(message, "selected checkpoint context is invalid") > 0, &
    "selected checkpoint negative composition rejection")
  call require_targets_unchanged("selected negative composition rollback")

  call replace_checkpoint_line( &
    selected_checkpoint_path, selected_bad_context_path, 7, &
    "NaN 0.0 0.0 0.0 0.0 0.0 1.0")
  call save_targets()
  call read_selected_checkpoint( &
    restart_config, reactions, transport, selected_bad_context_path, &
    selected_bundle_sha256, "explicit", selected_mole_fractions)
  call require(.not. ok .and. &
    index(message, "selected checkpoint context is invalid") > 0, &
    "selected checkpoint nonfinite composition rejection")
  call require_targets_unchanged("selected nonfinite composition rollback")

  call make_prefix_checkpoint( &
    selected_checkpoint_path, selected_truncated_path, 5)
  call save_targets()
  call read_selected_checkpoint( &
    restart_config, reactions, transport, selected_truncated_path, &
    selected_bundle_sha256, "explicit", selected_mole_fractions)
  call require(.not. ok, "selected checkpoint context truncation rejection")
  call require_targets_unchanged("selected context truncation rollback")

  call remove_file(checkpoint_path)
  call remove_file(truncated_path)
  call remove_file(early_truncated_path)
  call remove_file(bad_magic_path)
  call remove_file(bad_schema_path)
  call remove_file(bad_temperature_path)
  call remove_file(selected_checkpoint_path)
  call remove_file(selected_bad_context_path)
  call remove_file(selected_truncated_path)
  call remove_file(selected_invalid_write_path)
  write(*, '(a)') "test_reactive_eb_3d_checkpoint: PASS"

contains

  subroutine read_checkpoint(config, active_reactions, active_transport, path)
    type(reactive_eb_3d_config), intent(in) :: config
    type(elementary_reaction), intent(in) :: active_reactions(:)
    type(gas_transport_species), intent(in) :: active_transport(:)
    character(len=*), intent(in) :: path

    call read_reactive_eb_3d_checkpoint( &
      path, species, active_reactions, active_transport, config, &
      read_geometry, read_state, read_temperature, read_time, read_steps, &
      read_integrals, read_l1_integrals, read_elements, read_minimum_dt, &
      read_maximum_diffusivity, read_minimum_theta, ok, message)
  end subroutine read_checkpoint

  subroutine read_selected_checkpoint( &
      config, active_reactions, active_transport, path, bundle_sha256, &
      chemistry_integrator, base_mole_fractions)
    type(reactive_eb_3d_config), intent(in) :: config
    type(elementary_reaction), intent(in) :: active_reactions(:)
    type(gas_transport_species), intent(in) :: active_transport(:)
    character(len=*), intent(in) :: path, bundle_sha256
    character(len=*), intent(in) :: chemistry_integrator
    real(dp), intent(in) :: base_mole_fractions(:)

    call read_reactive_eb_3d_checkpoint( &
      path, species, active_reactions, active_transport, config, &
      read_geometry, read_state, read_temperature, read_time, read_steps, &
      read_integrals, read_l1_integrals, read_elements, read_minimum_dt, &
      read_maximum_diffusivity, read_minimum_theta, ok, message, &
      bundle_sha256=bundle_sha256, &
      chemistry_integrator=chemistry_integrator, &
      base_mole_fractions=base_mole_fractions)
  end subroutine read_selected_checkpoint

  subroutine require_checkpoint_schema(path, expected_schema, expected_marker)
    character(len=*), intent(in) :: path, expected_marker
    integer, intent(in) :: expected_schema
    character(len=1024) :: magic, marker
    integer :: unit, status, header(5)

    open(newunit=unit, file=path, status="old", action="read", iostat=status)
    call require(status == 0, "open checkpoint schema probe")
    read(unit, '(a)', iostat=status) magic
    call require(status == 0 .and. &
      trim(magic) == "PELEF_REACTIVE_EB_3D_CHECKPOINT", &
      "checkpoint schema probe magic")
    read(unit, *, iostat=status) header
    call require(status == 0 .and. header(1) == expected_schema, &
      "checkpoint schema probe version")
    read(unit, '(a)', iostat=status) marker
    call require(status == 0 .and. trim(marker) == trim(expected_marker), &
      "checkpoint schema probe first marker")
    close(unit)
  end subroutine require_checkpoint_schema

  subroutine require_selected_checkpoint_context( &
      path, expected_sha256, expected_integrator, expected_composition)
    character(len=*), intent(in) :: path, expected_sha256
    character(len=*), intent(in) :: expected_integrator
    real(dp), intent(in) :: expected_composition(:)
    character(len=1024) :: magic, marker, stored_sha256, stored_integrator
    real(dp), allocatable :: stored_composition(:)
    integer :: unit, status, header(5), stored_size

    open(newunit=unit, file=path, status="old", action="read", iostat=status)
    call require(status == 0, "open selected checkpoint context probe")
    read(unit, '(a)', iostat=status) magic
    call require(status == 0 .and. &
      trim(magic) == "PELEF_REACTIVE_EB_3D_CHECKPOINT", &
      "selected checkpoint context magic")
    read(unit, *, iostat=status) header
    call require(status == 0 .and. header(1) == 2, &
      "selected checkpoint context schema")
    read(unit, '(a)', iostat=status) marker
    call require(status == 0 .and. trim(marker) == "SELECTED_CONTEXT", &
      "selected checkpoint context marker")
    read(unit, '(a)', iostat=status) stored_sha256
    call require(status == 0 .and. &
      trim(stored_sha256) == trim(expected_sha256), &
      "selected checkpoint stored bundle SHA")
    read(unit, '(a)', iostat=status) stored_integrator
    call require(status == 0 .and. &
      trim(stored_integrator) == trim(expected_integrator), &
      "selected checkpoint stored integrator")
    read(unit, *, iostat=status) stored_size
    call require(status == 0 .and. &
      stored_size == size(expected_composition), &
      "selected checkpoint stored composition size")
    allocate(stored_composition(stored_size))
    read(unit, *, iostat=status) stored_composition
    call require(status == 0 .and. &
      all(stored_composition == expected_composition), &
      "selected checkpoint stored composition")
    close(unit)
  end subroutine require_selected_checkpoint_context

  subroutine require_file_absent(path, label)
    character(len=*), intent(in) :: path, label
    logical :: exists

    inquire(file=path, exist=exists)
    call require(.not. exists, label)
  end subroutine require_file_absent

  subroutine save_targets()
    saved_geometry = read_geometry
    saved_state = read_state
    saved_temperature = read_temperature
    saved_integrals = read_integrals
    saved_l1_integrals = read_l1_integrals
    saved_elements = read_elements
    saved_time = read_time
    saved_steps = read_steps
    saved_minimum_dt = read_minimum_dt
    saved_maximum_diffusivity = read_maximum_diffusivity
    saved_minimum_theta = read_minimum_theta
  end subroutine save_targets

  subroutine require_targets_unchanged(label)
    character(len=*), intent(in) :: label

    call require(same_geometry(read_geometry, saved_geometry) .and. &
      all(read_state == saved_state) .and. &
      all(read_temperature == saved_temperature) .and. &
      all(read_integrals == saved_integrals) .and. &
      all(read_l1_integrals == saved_l1_integrals) .and. &
      all(read_elements == saved_elements) .and. &
      read_time == saved_time .and. read_steps == saved_steps .and. &
      read_minimum_dt == saved_minimum_dt .and. &
      read_maximum_diffusivity == saved_maximum_diffusivity .and. &
      read_minimum_theta == saved_minimum_theta, label)
  end subroutine require_targets_unchanged

  logical function same_geometry(left, right) result(same)
    type(eb_geometry_3d), intent(in) :: left, right

    same = left%nx == right%nx .and. left%ny == right%ny .and. &
      left%nz == right%nz .and. left%x_lower == right%x_lower .and. &
      left%x_upper == right%x_upper .and. left%y_lower == right%y_lower .and. &
      left%y_upper == right%y_upper .and. left%z_lower == right%z_lower .and. &
      left%z_upper == right%z_upper .and. left%dx == right%dx .and. &
      left%dy == right%dy .and. left%dz == right%dz
    if (.not. same) return
    same = all(left%volume_fraction == right%volume_fraction) .and. &
      all(left%cell_centroid_x == right%cell_centroid_x) .and. &
      all(left%cell_centroid_y == right%cell_centroid_y) .and. &
      all(left%cell_centroid_z == right%cell_centroid_z) .and. &
      all(left%x_face_fraction == right%x_face_fraction) .and. &
      all(left%y_face_fraction == right%y_face_fraction) .and. &
      all(left%z_face_fraction == right%z_face_fraction) .and. &
      all(left%x_face_centroid_y == right%x_face_centroid_y) .and. &
      all(left%x_face_centroid_z == right%x_face_centroid_z) .and. &
      all(left%y_face_centroid_x == right%y_face_centroid_x) .and. &
      all(left%y_face_centroid_z == right%y_face_centroid_z) .and. &
      all(left%z_face_centroid_x == right%z_face_centroid_x) .and. &
      all(left%z_face_centroid_y == right%z_face_centroid_y) .and. &
      all(left%boundary_area == right%boundary_area) .and. &
      all(left%boundary_centroid_x == right%boundary_centroid_x) .and. &
      all(left%boundary_centroid_y == right%boundary_centroid_y) .and. &
      all(left%boundary_centroid_z == right%boundary_centroid_z) .and. &
      all(left%boundary_normal_x == right%boundary_normal_x) .and. &
      all(left%boundary_normal_y == right%boundary_normal_y) .and. &
      all(left%boundary_normal_z == right%boundary_normal_z) .and. &
      all(left%boundary_normal_integral_x == &
        right%boundary_normal_integral_x) .and. &
      all(left%boundary_normal_integral_y == &
        right%boundary_normal_integral_y) .and. &
      all(left%boundary_normal_integral_z == &
        right%boundary_normal_integral_z) .and. &
      all(left%cell_type == right%cell_type)
  end function same_geometry

  subroutine make_truncated_checkpoint(source, destination)
    character(len=*), intent(in) :: source, destination
    character(len=4096) :: line
    integer :: input_unit, output_unit, status

    open(newunit=input_unit, file=source, status="old", action="read")
    open(newunit=output_unit, file=destination, status="replace", &
      action="write")
    do
      read(input_unit, '(a)', iostat=status) line
      if (status /= 0) exit
      write(output_unit, '(a)') trim(line)
      if (trim(line) == "STATE") exit
    end do
    close(input_unit)
    close(output_unit)
  end subroutine make_truncated_checkpoint

  subroutine make_prefix_checkpoint(source, destination, kept_lines)
    character(len=*), intent(in) :: source, destination
    integer, intent(in) :: kept_lines
    character(len=4096) :: line
    integer :: input_unit, output_unit, status, line_number

    call require(kept_lines >= 0, "nonnegative checkpoint prefix length")
    open(newunit=input_unit, file=source, status="old", action="read")
    open(newunit=output_unit, file=destination, status="replace", &
      action="write")
    do line_number = 1, kept_lines
      read(input_unit, '(a)', iostat=status) line
      call require(status == 0, "checkpoint prefix source length")
      write(output_unit, '(a)') trim(line)
    end do
    close(input_unit)
    close(output_unit)
  end subroutine make_prefix_checkpoint

  subroutine replace_checkpoint_line( &
      source, destination, selected_line, replacement)
    character(len=*), intent(in) :: source, destination, replacement
    integer, intent(in) :: selected_line
    character(len=4096) :: line
    integer :: input_unit, output_unit, status, line_number
    logical :: changed

    line_number = 0
    changed = .false.
    open(newunit=input_unit, file=source, status="old", action="read")
    open(newunit=output_unit, file=destination, status="replace", &
      action="write")
    do
      read(input_unit, '(a)', iostat=status) line
      if (status /= 0) exit
      line_number = line_number + 1
      if (line_number == selected_line) then
        write(output_unit, '(a)') trim(replacement)
        changed = .true.
      else
        write(output_unit, '(a)') trim(line)
      end if
    end do
    close(input_unit)
    close(output_unit)
    call require(changed, "replace checkpoint line")
  end subroutine replace_checkpoint_line

  subroutine make_bad_temperature_checkpoint(source, destination, row_nvar)
    character(len=*), intent(in) :: source, destination
    integer, intent(in) :: row_nvar
    character(len=4096) :: line
    real(dp), allocatable :: row(:)
    integer :: input_unit, output_unit, status, stage, state_row
    logical :: changed

    allocate(row(row_nvar + 1))
    stage = 0
    state_row = 0
    changed = .false.
    open(newunit=input_unit, file=source, status="old", action="read")
    open(newunit=output_unit, file=destination, status="replace", &
      action="write")
    do
      read(input_unit, '(a)', iostat=status) line
      if (status /= 0) exit
      if (stage == 0 .and. trim(line) == "STATE") then
        stage = 1
        write(output_unit, '(a)') trim(line)
      else if (stage == 1) then
        stage = 2
        write(output_unit, '(a)') trim(line)
      else if (stage == 2) then
        state_row = state_row + 1
        if (state_row == 2) then
          read(line, *, iostat=status) row
          call require(status == 0, "parse checkpoint state row")
          row(row_nvar + 1) = 1.01_dp * row(row_nvar + 1)
          write(output_unit, '(*(es27.18e3,1x))') row
          changed = .true.
        else
          write(output_unit, '(a)') trim(line)
        end if
      else
        write(output_unit, '(a)') trim(line)
      end if
    end do
    close(input_unit)
    close(output_unit)
    call require(changed, "mutate checkpoint temperature")
  end subroutine make_bad_temperature_checkpoint

  subroutine remove_file(path)
    character(len=*), intent(in) :: path
    logical :: exists
    integer :: unit

    inquire(file=path, exist=exists)
    if (.not. exists) return
    open(newunit=unit, file=path, status="old")
    close(unit, status="delete")
  end subroutine remove_file

  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      write(*, '(a,1x,a)') trim(label), trim(message)
      error stop label
    end if
  end subroutine require

end program test_reactive_eb_3d_checkpoint
