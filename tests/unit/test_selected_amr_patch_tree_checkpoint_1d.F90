program test_selected_amr_patch_tree_checkpoint_1d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use simulation_config_reactive_1d_mod, only: &
    reactive_1d_config, reactive_1d_mole_fractions
  use reactive_1d_mod, only: reactive_nvar, reactive_species_component
  use amr_patch_tree_1d_mod, only: amr_patch_level_plan_1d
  use amr_patch_tree_reactive_1d_mod, only: &
    amr_patch_tree_reactive_solution_1d, &
    initialize_patch_tree_reactive_1d, &
    patch_tree_reactive_integrals_1d, &
    write_patch_tree_reactive_1d_checkpoint, &
    read_patch_tree_reactive_1d_checkpoint, &
    write_patch_tree_reactive_1d_selected_checkpoint, &
    read_patch_tree_reactive_1d_selected_checkpoint
  implicit none

  character(len=*), parameter :: selected_path = &
    "selected_amr_patch_tree_checkpoint_1d.chk"
  character(len=*), parameter :: fixed_path = &
    "fixed_amr_patch_tree_checkpoint_1d.chk"
  character(len=*), parameter :: trailing_path = &
    "selected_amr_patch_tree_checkpoint_1d_trailing.chk"
  character(len=*), parameter :: truncated_path = &
    "selected_amr_patch_tree_checkpoint_1d_truncated.chk"
  character(len=*), parameter :: bundle_sha256 = &
    "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef"

  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_level_plan_1d), allocatable :: empty_plans(:)
  type(amr_patch_tree_reactive_solution_1d) :: solution, restored
  type(amr_patch_tree_reactive_solution_1d) :: sentinel, fixed_candidate
  type(amr_patch_tree_reactive_solution_1d) :: invalid_solution
  type(reactive_1d_config) :: config
  real(dp), allocatable :: baseline(:), restored_baseline(:)
  real(dp), allocatable :: sentinel_baseline(:), invalid_baseline(:)
  real(dp), allocatable :: composition(:), mismatched_composition(:)
  real(dp) :: species_transfer
  character(len=64) :: mismatched_sha256
  logical :: ok
  integer :: unit

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "selected checkpoint thermodynamics load")
  call configure_case(config)
  allocate(composition(size(species)))
  call reactive_1d_mole_fractions(config, size(species), composition, ok)
  call require(ok, "selected checkpoint composition")
  allocate(empty_plans(0))
  call initialize_patch_tree_reactive_1d( &
    species, config, empty_plans, solution, ok, composition)
  call require(ok .and. solution%is_valid(), &
    "selected checkpoint initialization")
  allocate(baseline(reactive_nvar(size(species))))
  call patch_tree_reactive_integrals_1d(solution, baseline, ok)
  call require(ok, "selected checkpoint baseline")

  call write_patch_tree_reactive_1d_selected_checkpoint( &
    selected_path, species, solution, bundle_sha256, "implicit", &
    composition, baseline, ok)
  call require(ok, "selected checkpoint write")
  call read_patch_tree_reactive_1d_selected_checkpoint( &
    selected_path, species, config, bundle_sha256, "implicit", composition, &
    restored, restored_baseline, ok)
  call require(ok .and. restored%is_valid(), "selected checkpoint read")
  call require(same_solution_state(solution, restored), &
    "selected checkpoint state round trip")
  call require(allocated(restored_baseline) .and. &
    size(restored_baseline) == size(baseline) .and. &
    all(abs(restored_baseline - baseline) <= 0.0_dp), &
    "selected checkpoint baseline round trip")

  call read_patch_tree_reactive_1d_checkpoint( &
    selected_path, species, config, fixed_candidate, ok)
  call require(.not. ok, "fixed reader rejects selected schema")
  call write_patch_tree_reactive_1d_checkpoint( &
    fixed_path, species, solution, ok)
  call require(ok, "fixed checkpoint write")
  sentinel = solution
  sentinel_baseline = baseline
  call read_patch_tree_reactive_1d_selected_checkpoint( &
    fixed_path, species, config, bundle_sha256, "implicit", composition, &
    sentinel, sentinel_baseline, ok)
  call require_rejected_without_mutation( &
    ok, sentinel, solution, sentinel_baseline, baseline, &
    "selected reader rejects fixed schema")

  mismatched_sha256 = bundle_sha256
  mismatched_sha256(64:64) = "0"
  call read_patch_tree_reactive_1d_selected_checkpoint( &
    selected_path, species, config, mismatched_sha256, "implicit", &
    composition, sentinel, sentinel_baseline, ok)
  call require_rejected_without_mutation( &
    ok, sentinel, solution, sentinel_baseline, baseline, &
    "selected checkpoint SHA mismatch")
  call read_patch_tree_reactive_1d_selected_checkpoint( &
    selected_path, species, config, bundle_sha256, "explicit", composition, &
    sentinel, sentinel_baseline, ok)
  call require_rejected_without_mutation( &
    ok, sentinel, solution, sentinel_baseline, baseline, &
    "selected checkpoint integrator mismatch")
  mismatched_composition = composition
  mismatched_composition(1) = mismatched_composition(1) + 1.0e-6_dp
  mismatched_composition(2) = mismatched_composition(2) - 1.0e-6_dp
  call read_patch_tree_reactive_1d_selected_checkpoint( &
    selected_path, species, config, bundle_sha256, "implicit", &
    mismatched_composition, sentinel, sentinel_baseline, ok)
  call require_rejected_without_mutation( &
    ok, sentinel, solution, sentinel_baseline, baseline, &
    "selected checkpoint composition mismatch")

  call write_patch_tree_reactive_1d_selected_checkpoint( &
    trailing_path, species, solution, bundle_sha256, "implicit", &
    composition, baseline, ok)
  call require(ok, "selected trailing checkpoint write")
  open(newunit=unit, file=trailing_path, status="old", position="append", &
    action="write")
  write(unit, '(a)') "TRAILING_CONTENT"
  close(unit)
  call read_patch_tree_reactive_1d_selected_checkpoint( &
    trailing_path, species, config, bundle_sha256, "implicit", composition, &
    sentinel, sentinel_baseline, ok)
  call require_rejected_without_mutation( &
    ok, sentinel, solution, sentinel_baseline, baseline, &
    "selected checkpoint trailing content")

  open(newunit=unit, file=truncated_path, status="replace", action="write")
  write(unit, '(a)') "PELEF_PATCH_TREE_REACTIVE_1D_CHECKPOINT"
  write(unit, '(*(i0,1x))') 2, size(species), &
    reactive_nvar(size(species)), 1
  close(unit)
  call read_patch_tree_reactive_1d_selected_checkpoint( &
    truncated_path, species, config, bundle_sha256, "implicit", composition, &
    sentinel, sentinel_baseline, ok)
  call require_rejected_without_mutation( &
    ok, sentinel, solution, sentinel_baseline, baseline, &
    "selected checkpoint truncation")

  invalid_baseline = baseline
  invalid_baseline(6) = invalid_baseline(6) + 1.0_dp
  call write_patch_tree_reactive_1d_selected_checkpoint( &
    selected_path, species, solution, bundle_sha256, "implicit", &
    composition, invalid_baseline, ok)
  call require(.not. ok, "selected checkpoint invalid baseline rejected")
  call read_patch_tree_reactive_1d_selected_checkpoint( &
    selected_path, species, config, bundle_sha256, "implicit", composition, &
    restored, restored_baseline, ok)
  call require(ok .and. same_solution_state(solution, restored), &
    "invalid selected write preserves prior checkpoint")

  invalid_solution = solution
  invalid_solution%levels(1)%patches(1)%state(1, 1) = -1.0_dp
  call require_invalid_state_write_preserves_checkpoint( &
    invalid_solution, "selected checkpoint negative density rejected")
  invalid_solution = solution
  invalid_solution%levels(1)%patches(1)%state( &
    reactive_species_component(1), 1) = -1.0_dp
  call require_invalid_state_write_preserves_checkpoint( &
    invalid_solution, "selected checkpoint negative species rejected")
  invalid_solution = solution
  species_transfer = invalid_solution%levels(1)%patches(1)%state( &
    reactive_species_component(1), 1) + 4.0e-11_dp
  invalid_solution%levels(1)%patches(1)%state( &
    reactive_species_component(1), 1) = -4.0e-11_dp
  invalid_solution%levels(1)%patches(1)%state( &
    reactive_species_component(2), 1) = &
      invalid_solution%levels(1)%patches(1)%state( &
        reactive_species_component(2), 1) + species_transfer
  call require_invalid_state_write_preserves_checkpoint( &
    invalid_solution, "selected checkpoint near-floor species rejected")
  invalid_solution = solution
  invalid_solution%levels(1)%patches(1)%state( &
    reactive_species_component(1), 1) = &
      invalid_solution%levels(1)%patches(1)%state( &
        reactive_species_component(1), 1) + 1.0_dp
  call require_invalid_state_write_preserves_checkpoint( &
    invalid_solution, "selected checkpoint species closure rejected")

  write(*, '(a)') "test_selected_amr_patch_tree_checkpoint_1d: PASS"

contains

  subroutine configure_case(local_config)
    type(reactive_1d_config), intent(out) :: local_config

    local_config = reactive_1d_config()
    local_config%nx = 8
    local_config%x_lower = 0.0_dp
    local_config%x_upper = 0.008_dp
    local_config%final_time = 1.0e-9_dp
    local_config%problem = "uniform_reactor"
    local_config%chemistry_model = "elementary"
    local_config%chemistry_enabled = .false.
    local_config%transport_enabled = .false.
    local_config%initial_temperature = 1000.0_dp
    local_config%initial_pressure = 101325.0_dp
    local_config%x_h2 = 0.29570_dp
    local_config%x_h = 1.0e-5_dp
    local_config%x_o = 1.0e-5_dp
    local_config%x_o2 = 0.14784_dp
    local_config%x_oh = 1.0e-5_dp
    local_config%x_h2o = 0.0_dp
    local_config%x_n2 = 0.55643_dp
    local_config%amr_enabled = .true.
    local_config%amr_max_levels = 2
    local_config%amr_minimum_patch_cells = 2
  end subroutine configure_case

  pure logical function same_solution_state(left, right) result(same)
    type(amr_patch_tree_reactive_solution_1d), intent(in) :: left, right

    same = left%is_valid() .and. right%is_valid()
    if (.not. same) return
    same = left%level_count() == right%level_count() .and. &
      left%steps == right%steps .and. left%time == right%time
    if (.not. same) return
    same = all(abs( &
      left%levels(1)%patches(1)%state - &
      right%levels(1)%patches(1)%state) <= 0.0_dp)
  end function same_solution_state

  subroutine require_rejected_without_mutation( &
      read_ok, candidate, expected, candidate_baseline, expected_baseline, &
      label)
    logical, intent(in) :: read_ok
    type(amr_patch_tree_reactive_solution_1d), intent(in) :: candidate
    type(amr_patch_tree_reactive_solution_1d), intent(in) :: expected
    real(dp), intent(in) :: candidate_baseline(:), expected_baseline(:)
    character(len=*), intent(in) :: label

    call require(.not. read_ok .and. same_solution_state(candidate, expected) &
      .and. size(candidate_baseline) == size(expected_baseline) .and. &
      all(abs(candidate_baseline - expected_baseline) <= 0.0_dp), label)
  end subroutine require_rejected_without_mutation

  subroutine require_invalid_state_write_preserves_checkpoint( &
      invalid_state, label)
    type(amr_patch_tree_reactive_solution_1d), intent(in) :: invalid_state
    character(len=*), intent(in) :: label

    call write_patch_tree_reactive_1d_selected_checkpoint( &
      selected_path, species, invalid_state, bundle_sha256, "implicit", &
      composition, baseline, ok)
    call require(.not. ok, label)
    call read_patch_tree_reactive_1d_selected_checkpoint( &
      selected_path, species, config, bundle_sha256, "implicit", composition, &
      restored, restored_baseline, ok)
    call require(ok .and. same_solution_state(solution, restored), &
      trim(label) // " without replacement")
  end subroutine require_invalid_state_write_preserves_checkpoint

  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      write(*, '(a)') "FAILED: " // trim(label)
      error stop 1
    end if
  end subroutine require

end program test_selected_amr_patch_tree_checkpoint_1d
