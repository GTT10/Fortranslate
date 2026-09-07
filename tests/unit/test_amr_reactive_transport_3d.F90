program test_amr_reactive_transport_3d
  use, intrinsic :: ieee_arithmetic, only: &
    ieee_is_finite, ieee_quiet_nan, ieee_value
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use transport_database_mod, only: load_h2o2_elementary_transport
  use mesh_3d_mod, only: uniform_cell_centers_3d
  use reactive_1d_mod, only: reactive_nvar, reactive_species_component
  use simulation_config_reactive_3d_mod, only: reactive_3d_config
  use reactive_entropy_wave_3d_problem_mod, only: &
    initialize_reactive_problem_3d
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use reactive_transport_3d_mod, only: reactive_transport_timestep_3d
  use amr_hierarchy_3d_mod, only: &
    amr_patch_3d, initialize_amr_patch_3d, average_down_3d, &
    restrict_average_3d, composite_integrals_amr_3d
  use state_indices_mod, only: irho
  use amr_reactive_3d_mod, only: advance_amr_reactive_strang_3d
  use amr_reactive_transport_3d_mod, only: &
    compute_amr_reactive_transport_timestep_3d, &
    advance_amr_reactive_transport_3d, advance_amr_reactive_full_3d
  implicit none

  integer, parameter :: nx = 4, ny = 4, nz = 4, ratio = 2
  real(dp), parameter :: dx = 1.0e-3_dp
  real(dp), parameter :: dy = 1.0e-3_dp
  real(dp), parameter :: dz = 1.0e-3_dp
  real(dp), parameter :: interval = 1.0e-8_dp
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_3d) :: patch
  real(dp), allocatable :: coarse_state(:, :, :, :)
  real(dp), allocatable :: coarse_temperature(:, :, :)
  real(dp), allocatable :: fine_state(:, :, :, :)
  real(dp), allocatable :: fine_temperature(:, :, :)
  real(dp), allocatable :: saved_coarse(:, :, :, :)
  real(dp), allocatable :: saved_coarse_temperature(:, :, :)
  real(dp), allocatable :: saved_fine(:, :, :, :)
  real(dp), allocatable :: saved_fine_temperature(:, :, :)
  real(dp), allocatable :: reference_coarse(:, :, :, :)
  real(dp), allocatable :: reference_coarse_temperature(:, :, :)
  real(dp), allocatable :: reference_fine(:, :, :, :)
  real(dp), allocatable :: reference_fine_temperature(:, :, :)
  real(dp), allocatable :: initial_integrals(:), final_integrals(:)
  real(dp) :: dt, coarse_dt, fine_dt, maximum_diffusivity
  real(dp) :: coarse_diffusivity, fine_diffusivity
  real(dp) :: theta, reflux, reference_reflux, conservation_error
  real(dp) :: initial_temperature_range, final_temperature_range
  real(dp) :: nan_value
  logical :: ok

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "elementary thermodynamics load")
  call load_h2o2_elementary_mechanism(reactions, ok)
  call require(ok, "elementary mechanism load")
  call load_h2o2_elementary_transport(transport, ok)
  call require(ok, "elementary transport load")
  call initialize_amr_patch_3d( &
    nx, ny, nz, 2, 3, 2, 3, 2, 3, ratio, patch, ok)
  call require(ok .and. patch%is_strictly_interior(), &
    "strictly interior transport patch")

  call initialize_hierarchy(.true., coarse_state, coarse_temperature, &
    fine_state, fine_temperature, ok)
  call require(ok, "hotspot hierarchy initialization")
  call compute_amr_reactive_transport_timestep_3d( &
    species, transport, patch, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, dx, dy, dz, 0.35_dp, &
    .true., .true., .true., dt, maximum_diffusivity, ok)
  call require(ok, "AMR parabolic timestep")
  call reactive_transport_timestep_3d( &
    species, transport, coarse_state, coarse_temperature, nx, ny, nz, &
    dx, dy, dz, 0.35_dp, .true., .true., .true., &
    coarse_dt, coarse_diffusivity, ok)
  call require(ok, "coarse parabolic timestep")
  call reactive_transport_timestep_3d( &
    species, transport, fine_state, fine_temperature, &
    patch%fine_nx(), patch%fine_ny(), patch%fine_nz(), &
    dx / real(ratio, dp), dy / real(ratio, dp), dz / real(ratio, dp), &
    0.35_dp, .true., .true., .true., fine_dt, fine_diffusivity, ok)
  call require(ok, "fine parabolic timestep")
  call require(dt == min(coarse_dt, real(ratio**2, dp) * fine_dt), &
    "r-squared transport timestep scaling")
  call require(maximum_diffusivity == &
    max(coarse_diffusivity, fine_diffusivity), &
    "hierarchy maximum diffusivity")

  saved_coarse = coarse_state
  saved_coarse_temperature = coarse_temperature
  saved_fine = fine_state
  saved_fine_temperature = fine_temperature
  allocate(initial_integrals(reactive_nvar(size(species))))
  allocate(final_integrals(reactive_nvar(size(species))))
  call composite_integrals_amr_3d( &
    coarse_state, fine_state, patch, dx, dy, dz, initial_integrals, ok)
  call require(ok, "initial composite transport integral")
  initial_temperature_range = maxval(fine_temperature) - &
    minval(fine_temperature)
  call advance_amr_reactive_transport_3d( &
    species, transport, patch, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, dx, dy, dz, interval, &
    .true., .true., .true., .true., "mc", theta, reflux, ok)
  call require(ok, "transactional hierarchy transport")
  call require(theta > 0.0_dp .and. theta <= 1.0_dp, &
    "hierarchy species limiter bound")
  call require(reflux > 0.0_dp, "active transport interface reflux")
  call require(max(maxval(abs(coarse_state - saved_coarse)), &
    maxval(abs(fine_state - saved_fine))) > 1.0e-12_dp, &
    "transport changes hotspot hierarchy")
  final_temperature_range = maxval(fine_temperature) - &
    minval(fine_temperature)
  call require(final_temperature_range < initial_temperature_range, &
    "fine hotspot thermal smoothing")
  call require(levels_are_synchronized(coarse_state, fine_state), &
    "transport average-down synchronization")
  call composite_integrals_amr_3d( &
    coarse_state, fine_state, patch, dx, dy, dz, final_integrals, ok)
  call require(ok, "final composite transport integral")
  conservation_error = maxval(abs(final_integrals - initial_integrals) / &
    max(1.0_dp, abs(initial_integrals)))
  call require(conservation_error <= 5.0e-13_dp, &
    "composite transport conservation")

  saved_coarse = coarse_state
  saved_coarse_temperature = coarse_temperature
  saved_fine = fine_state
  saved_fine_temperature = fine_temperature
  call advance_amr_reactive_transport_3d( &
    species, transport, patch, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, dx, dy, dz, -interval, &
    .true., .true., .true., .true., "mc", theta, reflux, ok)
  call require(.not. ok .and. hierarchy_matches_saved(), &
    "negative interval rollback")
  call advance_amr_reactive_transport_3d( &
    species, transport, patch, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, dx, dy, dz, interval, &
    .true., .true., .false., .true., "mc", theta, reflux, ok)
  call require(.not. ok .and. hierarchy_matches_saved(), &
    "invalid barodiffusion rollback")
  nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
  call compute_amr_reactive_transport_timestep_3d( &
    species, transport, patch, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, nan_value, dy, dz, 0.35_dp, &
    .true., .true., .true., dt, maximum_diffusivity, ok)
  call require(.not. ok .and. dt == 0.0_dp .and. &
    maximum_diffusivity == 0.0_dp, "NaN spacing timestep rejection")
  call advance_amr_reactive_transport_3d( &
    species, transport, patch, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, dx, dy, dz, nan_value, &
    .true., .true., .true., .true., "mc", theta, reflux, ok)
  call require(.not. ok .and. hierarchy_matches_saved(), &
    "NaN interval rollback")
  call advance_amr_reactive_full_3d( &
    species, reactions, transport, patch, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, dx, dy, dz, interval, "rusanov", &
    .true., nan_value, 1.0e-12_dp, .true., .true., .true., .true., .true., &
    theta, reflux, ok, "characteristic_plm", "mc", "implicit")
  call require(.not. ok .and. hierarchy_matches_saved(), &
    "NaN chemistry tolerance rollback")

  call initialize_hierarchy(.false., coarse_state, coarse_temperature, &
    fine_state, fine_temperature, ok)
  call require(ok, "uniform hierarchy initialization")
  saved_coarse = coarse_state
  saved_coarse_temperature = coarse_temperature
  saved_fine = fine_state
  saved_fine_temperature = fine_temperature
  call advance_amr_reactive_transport_3d( &
    species, transport, patch, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, dx, dy, dz, interval, &
    .true., .true., .true., .true., "mc", theta, reflux, ok)
  call require(ok, "uniform hierarchy transport")
  call require(maxval(abs(coarse_state - saved_coarse)) == 0.0_dp .and. &
    maxval(abs(fine_state - saved_fine)) == 0.0_dp, &
    "uniform transport state no-op")
  call require(maxval(abs(coarse_temperature - &
      saved_coarse_temperature)) <= 5.0e-12_dp .and. &
    maxval(abs(fine_temperature - saved_fine_temperature)) <= 5.0e-12_dp, &
    "uniform transport temperature no-op")

  call initialize_hierarchy(.true., coarse_state, coarse_temperature, &
    fine_state, fine_temperature, ok)
  call require(ok, "hydro parity hierarchy initialization")
  reference_coarse = coarse_state
  reference_coarse_temperature = coarse_temperature
  reference_fine = fine_state
  reference_fine_temperature = fine_temperature
  call advance_amr_reactive_strang_3d( &
    species, reactions, patch, reference_coarse, &
    reference_coarse_temperature, reference_fine, &
    reference_fine_temperature, dx, dy, dz, interval, "rusanov", .false., &
    2.0e-7_dp, 1.0e-12_dp, reference_reflux, ok, &
    "characteristic_plm", "mc")
  call require(ok, "legacy hydro-only reference")
  call advance_amr_reactive_full_3d( &
    species, reactions, transport, patch, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, dx, dy, dz, interval, "rusanov", &
    .false., 2.0e-7_dp, 1.0e-12_dp, .false., &
    .true., .true., .true., .true., theta, reflux, ok, &
    "characteristic_plm", "mc")
  call require(ok, "transport-disabled full operator")
  call require(all(coarse_state == reference_coarse) .and. &
    all(coarse_temperature == reference_coarse_temperature) .and. &
    all(fine_state == reference_fine) .and. &
    all(fine_temperature == reference_fine_temperature) .and. &
    reflux == reference_reflux .and. theta == 1.0_dp, &
    "transport-disabled legacy bit path")

  call run_adversarial_interface_tests(species, transport, ok)
  call require(ok, "adversarial six-face interface transport tests")
  call run_individual_transport_mode_tests(species, transport, ok)
  call require(ok, "individual AMR transport process tests")
  call run_full_operator_test(species, reactions, transport, ok)
  call require(ok, "active AMR R-T-H-T-R transaction test")
  call run_ratio_three_transport_test(species, transport, ok)
  call require(ok, "ratio-three AMR transport test")

  write(*, '(a)') "test_amr_reactive_transport_3d: PASS"

contains

  subroutine run_individual_transport_mode_tests(species, transport, ok_out)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    logical, intent(out) :: ok_out

    logical, parameter :: viscosity(3) = [.true., .false., .false.]
    logical, parameter :: conduction(3) = [.false., .true., .false.]
    logical, parameter :: diffusion(3) = [.false., .false., .true.]
    character(len=12), parameter :: mode_name(3) = [ &
      character(len=12) :: "viscosity", "conduction", "diffusion"]
    real(dp), allocatable :: case_coarse(:, :, :, :)
    real(dp), allocatable :: case_coarse_temperature(:, :, :)
    real(dp), allocatable :: case_fine(:, :, :, :)
    real(dp), allocatable :: case_fine_temperature(:, :, :)
    real(dp), allocatable :: before(:), after(:)
    real(dp) :: mode_theta, mode_reflux, mode_error
    logical :: local_ok
    integer :: mode

    ok_out = .false.
    allocate(before(reactive_nvar(size(species))))
    allocate(after(reactive_nvar(size(species))))
    do mode = 1, 3
      call initialize_hierarchy( &
        .true., case_coarse, case_coarse_temperature, case_fine, &
        case_fine_temperature, local_ok)
      if (.not. local_ok) return
      call composite_integrals_amr_3d( &
        case_coarse, case_fine, patch, dx, dy, dz, before, local_ok)
      if (.not. local_ok) return
      call advance_amr_reactive_transport_3d( &
        species, transport, patch, case_coarse, case_coarse_temperature, &
        case_fine, case_fine_temperature, dx, dy, dz, interval, &
        viscosity(mode), conduction(mode), diffusion(mode), .false., &
        "mc", mode_theta, mode_reflux, local_ok)
      if (.not. local_ok) return
      call composite_integrals_amr_3d( &
        case_coarse, case_fine, patch, dx, dy, dz, after, local_ok)
      if (.not. local_ok) return
      mode_error = maxval(abs(after - before) / max(1.0_dp, abs(before)))
      call require(mode_theta >= 0.0_dp .and. mode_theta <= 1.0_dp, &
        trim(mode_name(mode)) // " limiter bound")
      call require(mode_reflux >= 0.0_dp .and. &
        ieee_is_finite(mode_reflux), trim(mode_name(mode)) // &
        " reflux bound")
      call require(mode_error <= 5.0e-13_dp, &
        trim(mode_name(mode)) // " composite conservation")
      call require(all(ieee_is_finite(case_coarse)) .and. &
        all(ieee_is_finite(case_coarse_temperature)) .and. &
        all(ieee_is_finite(case_fine)) .and. &
        all(ieee_is_finite(case_fine_temperature)), &
        trim(mode_name(mode)) // " finite hierarchy")
    end do
    ok_out = .true.
  end subroutine run_individual_transport_mode_tests

  subroutine run_full_operator_test( &
      species, reactions, transport, ok_out)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    logical, intent(out) :: ok_out

    real(dp), allocatable :: case_coarse(:, :, :, :)
    real(dp), allocatable :: case_coarse_temperature(:, :, :)
    real(dp), allocatable :: case_fine(:, :, :, :)
    real(dp), allocatable :: case_fine_temperature(:, :, :)
    real(dp), allocatable :: before_coarse(:, :, :, :)
    real(dp), allocatable :: before_coarse_temperature(:, :, :)
    real(dp), allocatable :: before_fine(:, :, :, :)
    real(dp), allocatable :: before_fine_temperature(:, :, :)
    real(dp) :: full_theta, full_reflux
    logical :: local_ok

    ok_out = .false.
    call initialize_hierarchy( &
      .true., case_coarse, case_coarse_temperature, case_fine, &
      case_fine_temperature, local_ok)
    if (.not. local_ok) return
    allocate(before_coarse, source=case_coarse)
    allocate(before_coarse_temperature, source=case_coarse_temperature)
    allocate(before_fine, source=case_fine)
    allocate(before_fine_temperature, source=case_fine_temperature)
    call advance_amr_reactive_full_3d( &
      species, reactions, transport, patch, case_coarse, &
      case_coarse_temperature, case_fine, case_fine_temperature, &
      dx, dy, dz, interval, "rusanov", .true., 2.0e-7_dp, 1.0e-12_dp, &
      .true., .true., .true., .true., .true., full_theta, full_reflux, &
      local_ok, "characteristic_plm", "mc", "implicit")
    if (.not. local_ok) return
    call require(full_theta >= 0.0_dp .and. full_theta <= 1.0_dp, &
      "full operator transport limiter bound")
    call require(full_reflux > 0.0_dp .and. ieee_is_finite(full_reflux), &
      "full operator active reflux")
    call require(levels_are_synchronized(case_coarse, case_fine), &
      "full operator average-down synchronization")
    call require(all(ieee_is_finite(case_coarse)) .and. &
      all(ieee_is_finite(case_coarse_temperature)) .and. &
      all(ieee_is_finite(case_fine)) .and. &
      all(ieee_is_finite(case_fine_temperature)), &
      "full operator finite hierarchy")
    call require(max(maxval(abs(case_coarse - before_coarse)), &
      maxval(abs(case_fine - before_fine))) > 1.0e-12_dp, &
      "full operator changes hierarchy")

    case_coarse = before_coarse
    case_coarse_temperature = before_coarse_temperature
    case_fine = before_fine
    case_fine_temperature = before_fine_temperature
    call advance_amr_reactive_full_3d( &
      species, reactions, transport, patch, case_coarse, &
      case_coarse_temperature, case_fine, case_fine_temperature, &
      dx, dy, dz, interval, "invalid", .true., 2.0e-7_dp, 1.0e-12_dp, &
      .true., .true., .true., .true., .true., full_theta, full_reflux, &
      local_ok, "characteristic_plm", "mc", "implicit")
    call require(.not. local_ok .and. full_theta == 1.0_dp .and. &
      full_reflux == 0.0_dp .and. &
      all(case_coarse == before_coarse) .and. &
      all(case_coarse_temperature == before_coarse_temperature) .and. &
      all(case_fine == before_fine) .and. &
      all(case_fine_temperature == before_fine_temperature), &
      "full operator middle-hydro rollback")
    ok_out = .true.
  end subroutine run_full_operator_test

  subroutine run_ratio_three_transport_test(species, transport, ok_out)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    logical, intent(out) :: ok_out

    type(amr_patch_3d) :: ratio_patch
    real(dp), allocatable :: case_coarse(:, :, :, :)
    real(dp), allocatable :: case_coarse_temperature(:, :, :)
    real(dp), allocatable :: case_fine(:, :, :, :)
    real(dp), allocatable :: case_fine_temperature(:, :, :)
    real(dp), allocatable :: before(:), after(:)
    real(dp), allocatable :: restricted(:, :, :, :)
    real(dp) :: ratio_theta, ratio_reflux, ratio_error
    logical :: local_ok

    ok_out = .false.
    call initialize_amr_patch_3d( &
      nx, ny, nz, 2, 3, 2, 3, 2, 3, 3, ratio_patch, local_ok)
    if (.not. local_ok) return
    call initialize_hierarchy( &
      .true., case_coarse, case_coarse_temperature, case_fine, &
      case_fine_temperature, local_ok, ratio_patch)
    if (.not. local_ok) return
    allocate(before(reactive_nvar(size(species))))
    allocate(after(reactive_nvar(size(species))))
    call composite_integrals_amr_3d( &
      case_coarse, case_fine, ratio_patch, dx, dy, dz, before, local_ok)
    if (.not. local_ok) return
    call advance_amr_reactive_transport_3d( &
      species, transport, ratio_patch, case_coarse, &
      case_coarse_temperature, case_fine, case_fine_temperature, &
      dx, dy, dz, interval, .true., .true., .true., .true., "mc", &
      ratio_theta, ratio_reflux, local_ok)
    if (.not. local_ok) return
    call composite_integrals_amr_3d( &
      case_coarse, case_fine, ratio_patch, dx, dy, dz, after, local_ok)
    if (.not. local_ok) return
    ratio_error = maxval(abs(after - before) / max(1.0_dp, abs(before)))
    allocate(restricted(size(case_coarse, 1), 2, 2, 2))
    call restrict_average_3d(case_fine, ratio_patch, restricted, local_ok)
    if (.not. local_ok) return
    call require(all(restricted == case_coarse(:, 2:3, 2:3, 2:3)), &
      "ratio-three average-down synchronization")
    call require(ratio_theta >= 0.0_dp .and. ratio_theta <= 1.0_dp, &
      "ratio-three limiter bound")
    call require(ratio_reflux > 0.0_dp .and. &
      ieee_is_finite(ratio_reflux), "ratio-three active reflux")
    call require(ratio_error <= 5.0e-13_dp, &
      "ratio-three composite conservation")
    ok_out = .true.
  end subroutine run_ratio_three_transport_test

  subroutine run_adversarial_interface_tests(species, transport, ok_out)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    logical, intent(out) :: ok_out

    character(len=8), parameter :: orientation_name(6) = [ &
      character(len=8) :: "x-lower", "x-upper", "y-lower", "y-upper", &
      "z-lower", "z-upper"]
    real(dp), parameter :: adversarial_interval = 3.0e-5_dp
    real(dp), parameter :: adversarial_low_fraction = 1.0e-12_dp
    real(dp), parameter :: adversarial_high_fraction = 0.70_dp
    real(dp), parameter :: conservation_tolerance = 2.0e-11_dp
    real(dp), allocatable :: case_coarse(:, :, :, :)
    real(dp), allocatable :: case_coarse_temperature(:, :, :)
    real(dp), allocatable :: case_fine(:, :, :, :)
    real(dp), allocatable :: case_fine_temperature(:, :, :)
    real(dp), allocatable :: case_initial_integrals(:)
    real(dp), allocatable :: case_final_integrals(:)
    real(dp), allocatable :: recovered_coarse_temperature(:, :, :)
    real(dp), allocatable :: recovered_fine_temperature(:, :, :)
    real(dp) :: case_theta, case_reflux, case_conservation_error
    real(dp) :: case_minimum_species, density_scale
    integer :: orientation, nvar, first_species, last_species
    logical :: local_ok

    ok_out = .false.
    nvar = reactive_nvar(size(species))
    first_species = reactive_species_component(1)
    last_species = reactive_species_component(size(species))
    allocate(case_initial_integrals(nvar), case_final_integrals(nvar))
    allocate(recovered_coarse_temperature(nx, ny, nz))
    allocate(recovered_fine_temperature( &
      patch%fine_nx(), patch%fine_ny(), patch%fine_nz()))

    do orientation = 1, 6
      call initialize_hierarchy( &
        .false., case_coarse, case_coarse_temperature, case_fine, &
        case_fine_temperature, local_ok)
      if (.not. local_ok) then
        write(*, '(a,i0)') "adversarial init failed orientation ", orientation
        return
      end if
      call configure_adversarial_interface_case( &
        species, orientation, case_coarse, case_coarse_temperature, &
        case_fine, case_fine_temperature, adversarial_low_fraction, &
        adversarial_high_fraction, local_ok)
      if (.not. local_ok) then
        write(*, '(a,i0)') "adversarial configure failed orientation ", &
          orientation
        return
      end if
      call require(all(ieee_is_finite(case_coarse)) .and. &
        all(ieee_is_finite(case_fine)) .and. &
        minval(case_coarse(first_species:last_species, :, :, :)) >= 0.0_dp &
        .and. minval(case_fine(first_species:last_species, :, :, :)) >= &
        0.0_dp, "adversarial initial nonnegative species")
      call composite_integrals_amr_3d( &
        case_coarse, case_fine, patch, dx, dy, dz, &
        case_initial_integrals, local_ok)
      if (.not. local_ok) then
        write(*, '(a,i0)') "adversarial initial integrals failed orientation ", &
          orientation
        return
      end if

      call advance_amr_reactive_transport_3d( &
        species, transport, patch, case_coarse, case_coarse_temperature, &
        case_fine, case_fine_temperature, dx, dy, dz, &
        adversarial_interval, .false., .false., .true., .false., "mc", &
        case_theta, case_reflux, local_ok)
      if (.not. local_ok) then
        write(*, '(a,i0,1x,a,es12.4,1x,a,es12.4)') &
          "adversarial advance failed orientation ", orientation, &
          "theta=", case_theta, "reflux=", case_reflux
        return
      end if

      call require(case_theta >= 0.0_dp .and. case_theta < 1.0_dp, &
        "active species limiter at " // trim(orientation_name(orientation)))
      call require(case_reflux > 0.0_dp, &
        "nonzero interface reflux at " // trim(orientation_name(orientation)))
      call require(all(ieee_is_finite(case_coarse)) .and. &
        all(ieee_is_finite(case_fine)) .and. &
        minval(case_coarse(first_species:last_species, :, :, :)) >= &
        -1.0e-14_dp .and. &
        minval(case_fine(first_species:last_species, :, :, :)) >= &
        -1.0e-14_dp, &
        "nonnegative species after " // trim(orientation_name(orientation)))

      call recover_reactive_temperatures_3d( &
        species, case_coarse, case_coarse_temperature, nx, ny, nz, &
        recovered_coarse_temperature, local_ok)
      call require(local_ok .and. all(ieee_is_finite( &
        recovered_coarse_temperature)) .and. &
        minval(recovered_coarse_temperature) > 0.0_dp, &
        "coarse temperature recovery at " // &
        trim(orientation_name(orientation)))
      call recover_reactive_temperatures_3d( &
        species, case_fine, case_fine_temperature, patch%fine_nx(), &
        patch%fine_ny(), patch%fine_nz(), recovered_fine_temperature, &
        local_ok)
      call require(local_ok .and. all(ieee_is_finite( &
        recovered_fine_temperature)) .and. &
        minval(recovered_fine_temperature) > 0.0_dp, &
        "fine temperature recovery at " // &
        trim(orientation_name(orientation)))

      call composite_integrals_amr_3d( &
        case_coarse, case_fine, patch, dx, dy, dz, &
        case_final_integrals, local_ok)
      if (.not. local_ok) return
      case_conservation_error = maxval(abs( &
        case_final_integrals - case_initial_integrals) / &
        max(1.0_dp, abs(case_initial_integrals)))
      call require(case_conservation_error <= conservation_tolerance, &
        "composite conservation at " // trim(orientation_name(orientation)))

      case_minimum_species = min( &
        minval(case_coarse(first_species:last_species, :, :, :)), &
        minval(case_fine(first_species:last_species, :, :, :)))
      density_scale = max(1.0_dp, maxval(abs( &
        case_coarse(irho, :, :, :))))
      call require(case_minimum_species >= -1.0e-14_dp * density_scale, &
        "species floor at " // trim(orientation_name(orientation)))
      write(*, '(a,1x,a,1x,a,es12.4,1x,a,es12.4,1x,a,es12.4,1x,a,es12.4)') &
        "adversarial", trim(orientation_name(orientation)), "theta=", &
        case_theta, "reflux=", case_reflux, "conservation=", &
        case_conservation_error, "min_species=", case_minimum_species
    end do
    ok_out = .true.
  end subroutine run_adversarial_interface_tests

  subroutine configure_adversarial_interface_case( &
      species, orientation, case_coarse, case_coarse_temperature, &
      case_fine, case_fine_temperature, low_fraction, high_fraction, &
      ok_out)
    type(nasa7_species), intent(in) :: species(:)
    integer, intent(in) :: orientation
    real(dp), intent(inout) :: case_coarse(:, :, :, :)
    real(dp), intent(inout) :: case_coarse_temperature(:, :, :)
    real(dp), intent(inout) :: case_fine(:, :, :, :)
    real(dp), intent(inout) :: case_fine_temperature(:, :, :)
    real(dp), intent(in) :: low_fraction, high_fraction
    logical, intent(out) :: ok_out

    real(dp), allocatable :: recovered(:, :, :)
    real(dp) :: fraction
    logical :: local_ok
    integer :: i, j, k, coarse_i, coarse_j, coarse_k, ratio

    ok_out = .false.
    if (orientation < 1 .or. orientation > 6 .or. &
        low_fraction <= 0.0_dp .or. high_fraction <= low_fraction .or. &
        high_fraction >= 1.0_dp) return
    if (any(shape(case_coarse) /= [reactive_nvar(size(species)), nx, ny, nz]) .or. &
        any(shape(case_coarse_temperature) /= [nx, ny, nz]) .or. &
        any(shape(case_fine) /= [reactive_nvar(size(species)), &
          patch%fine_nx(), patch%fine_ny(), patch%fine_nz()]) .or. &
        any(shape(case_fine_temperature) /= &
          [patch%fine_nx(), patch%fine_ny(), patch%fine_nz()])) return

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (adversarial_low_cell(orientation, i, j, k)) then
            fraction = low_fraction
          else
            fraction = high_fraction
          end if
          call set_adversarial_species_fraction( &
            case_coarse(:, i, j, k), size(species), fraction, local_ok)
          if (.not. local_ok) return
        end do
      end do
    end do

    ratio = patch%refinement_ratio
    do k = 1, patch%fine_nz()
      coarse_k = patch%coarse_k_lower + (k - 1) / ratio
      do j = 1, patch%fine_ny()
        coarse_j = patch%coarse_j_lower + (j - 1) / ratio
        do i = 1, patch%fine_nx()
          coarse_i = patch%coarse_i_lower + (i - 1) / ratio
          case_fine(:, i, j, k) = case_coarse(:, coarse_i, coarse_j, coarse_k)
        end do
      end do
    end do

    allocate(recovered(nx, ny, nz))
    call recover_reactive_temperatures_3d( &
      species, case_coarse, case_coarse_temperature, nx, ny, nz, &
      recovered, local_ok)
    if (.not. local_ok) return
    case_coarse_temperature = recovered
    deallocate(recovered)
    allocate(recovered(patch%fine_nx(), patch%fine_ny(), patch%fine_nz()))
    call recover_reactive_temperatures_3d( &
      species, case_fine, case_fine_temperature, patch%fine_nx(), &
      patch%fine_ny(), patch%fine_nz(), recovered, local_ok)
    if (.not. local_ok) return
    case_fine_temperature = recovered
    ok_out = all(ieee_is_finite(case_coarse_temperature)) .and. &
      all(ieee_is_finite(case_fine_temperature)) .and. &
      minval(case_coarse_temperature) > 0.0_dp .and. &
      minval(case_fine_temperature) > 0.0_dp
  end subroutine configure_adversarial_interface_case

  logical function adversarial_low_cell(orientation, i, j, k) result(low)
    integer, intent(in) :: orientation, i, j, k

    low = .false.
    select case (orientation)
    case (1)
      low = i == patch%coarse_i_lower - 1
    case (2)
      low = i == patch%coarse_i_upper + 1
    case (3)
      low = j == patch%coarse_j_lower - 1
    case (4)
      low = j == patch%coarse_j_upper + 1
    case (5)
      low = k == patch%coarse_k_lower - 1
    case (6)
      low = k == patch%coarse_k_upper + 1
    end select
  end function adversarial_low_cell

  subroutine set_adversarial_species_fraction(state, nspecies, fraction, ok_out)
    real(dp), intent(inout) :: state(:)
    integer, intent(in) :: nspecies
    real(dp), intent(in) :: fraction
    logical, intent(out) :: ok_out

    integer :: target_component, reservoir_component
    real(dp) :: rho, other_species_mass, reservoir_mass

    ok_out = .false.
    target_component = reactive_species_component(1)
    reservoir_component = reactive_species_component(nspecies)
    if (size(state) /= reactive_nvar(nspecies) .or. &
        fraction <= 0.0_dp .or. fraction >= 1.0_dp .or. &
        .not. all(ieee_is_finite(state))) return
    rho = state(irho)
    other_species_mass = rho - state(target_component) - &
      state(reservoir_component)
    reservoir_mass = rho * (1.0_dp - fraction) - other_species_mass
    if (rho <= 0.0_dp .or. reservoir_mass < 0.0_dp .or. &
        .not. ieee_is_finite(reservoir_mass)) return
    state(target_component) = rho * fraction
    state(reservoir_component) = reservoir_mass
    ok_out = all(ieee_is_finite(state)) .and. &
      minval(state(target_component:reservoir_component)) >= 0.0_dp
  end subroutine set_adversarial_species_fraction

  subroutine initialize_hierarchy( &
      hotspot, output_coarse, output_coarse_temperature, &
      output_fine, output_fine_temperature, ok_out, patch_override)
    logical, intent(in) :: hotspot
    real(dp), allocatable, intent(out) :: output_coarse(:, :, :, :)
    real(dp), allocatable, intent(out) :: output_coarse_temperature(:, :, :)
    real(dp), allocatable, intent(out) :: output_fine(:, :, :, :)
    real(dp), allocatable, intent(out) :: output_fine_temperature(:, :, :)
    logical, intent(out) :: ok_out
    type(amr_patch_3d), intent(in), optional :: patch_override

    type(amr_patch_3d) :: active_patch
    type(reactive_3d_config) :: config, fine_config
    real(dp), allocatable :: x(:), y(:), z(:), xf(:), yf(:), zf(:)
    real(dp), allocatable :: mass_fractions(:), recovered(:, :, :)
    real(dp) :: local_dx, local_dy, local_dz, density, fine_density
    logical :: local_ok
    integer :: nvar

    ok_out = .false.
    active_patch = patch
    if (present(patch_override)) active_patch = patch_override
    config = reactive_3d_config()
    config%nx = nx
    config%ny = ny
    config%nz = nz
    config%x_upper = real(nx, dp) * dx
    config%y_upper = real(ny, dp) * dy
    config%z_upper = real(nz, dp) * dz
    config%problem = "uniform_reactor"
    if (hotspot) config%problem = "reactive_hotspot"
    config%thermo_model = "elementary"
    config%initial_temperature = 1200.0_dp
    config%initial_velocity_x = 15.0_dp
    config%initial_velocity_y = -7.0_dp
    config%initial_velocity_z = 3.0_dp
    config%hotspot_temperature_rise = 200.0_dp
    config%hotspot_center_x = 0.5_dp * config%x_upper
    config%hotspot_center_y = 0.5_dp * config%y_upper
    config%hotspot_center_z = 0.5_dp * config%z_upper
    config%hotspot_width = 0.6e-3_dp
    nvar = reactive_nvar(size(species))
    allocate(output_coarse(nvar, nx, ny, nz))
    allocate(output_coarse_temperature(nx, ny, nz))
    allocate(output_fine(nvar, active_patch%fine_nx(), &
      active_patch%fine_ny(), active_patch%fine_nz()))
    allocate(output_fine_temperature( &
      active_patch%fine_nx(), active_patch%fine_ny(), &
      active_patch%fine_nz()))
    allocate(x(nx), y(ny), z(nz))
    allocate(xf(active_patch%fine_nx()), yf(active_patch%fine_ny()), &
      zf(active_patch%fine_nz()))
    allocate(mass_fractions(size(species)))
    allocate(recovered(nx, ny, nz))
    call uniform_cell_centers_3d( &
      nx, ny, nz, config%x_lower, config%x_upper, &
      config%y_lower, config%y_upper, config%z_lower, config%z_upper, &
      x, y, z, local_dx, local_dy, local_dz)
    call initialize_reactive_problem_3d( &
      species, config, x, y, z, output_coarse, output_coarse_temperature, &
      density, mass_fractions, local_ok)
    if (.not. local_ok) return
    fine_config = config
    fine_config%nx = active_patch%fine_nx()
    fine_config%ny = active_patch%fine_ny()
    fine_config%nz = active_patch%fine_nz()
    fine_config%x_lower = real(active_patch%coarse_i_lower - 1, dp) * dx
    fine_config%x_upper = real(active_patch%coarse_i_upper, dp) * dx
    fine_config%y_lower = real(active_patch%coarse_j_lower - 1, dp) * dy
    fine_config%y_upper = real(active_patch%coarse_j_upper, dp) * dy
    fine_config%z_lower = real(active_patch%coarse_k_lower - 1, dp) * dz
    fine_config%z_upper = real(active_patch%coarse_k_upper, dp) * dz
    call uniform_cell_centers_3d( &
      fine_config%nx, fine_config%ny, fine_config%nz, &
      fine_config%x_lower, fine_config%x_upper, fine_config%y_lower, &
      fine_config%y_upper, fine_config%z_lower, fine_config%z_upper, &
      xf, yf, zf, local_dx, local_dy, local_dz)
    call initialize_reactive_problem_3d( &
      species, fine_config, xf, yf, zf, output_fine, &
      output_fine_temperature, fine_density, mass_fractions, local_ok)
    if (.not. local_ok) return
    call average_down_3d( &
      output_coarse, output_fine, active_patch, local_ok)
    if (.not. local_ok) return
    call recover_reactive_temperatures_3d( &
      species, output_coarse, output_coarse_temperature, nx, ny, nz, &
      recovered, local_ok)
    if (.not. local_ok) return
    output_coarse_temperature = recovered
    ok_out = .true.
  end subroutine initialize_hierarchy

  logical function levels_are_synchronized(coarse, fine) result(match)
    real(dp), intent(in) :: coarse(:, :, :, :), fine(:, :, :, :)

    real(dp), allocatable :: restricted(:, :, :, :)
    logical :: local_ok

    allocate(restricted(size(coarse, 1), 2, 2, 2))
    call restrict_average_3d(fine, patch, restricted, local_ok)
    match = local_ok
    if (match) match = all(restricted == coarse(:, 2:3, 2:3, 2:3))
  end function levels_are_synchronized

  logical function hierarchy_matches_saved() result(match)
    match = all(coarse_state == saved_coarse) .and. &
      all(coarse_temperature == saved_coarse_temperature) .and. &
      all(fine_state == saved_fine) .and. &
      all(fine_temperature == saved_fine_temperature)
  end function hierarchy_matches_saved

  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      write(*, '(a,1x,a)') "FAILED:", trim(label)
      error stop 1
    end if
  end subroutine require

end program test_amr_reactive_transport_3d
