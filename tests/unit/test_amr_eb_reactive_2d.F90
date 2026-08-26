program test_amr_eb_reactive_2d
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use precision_mod, only: dp
  use state_indices_mod, only: irho, iet
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_mass_fraction_component, reactive_primitive_to_conserved
  use eb_geometry_2d_mod, only: &
    eb_geometry_2d, eb_covered_cell, build_eb_geometry_2d
  use eb_reactive_reconstruction_2d_mod, only: &
    reactive_eb_exterior_state_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, build_amr_eb_patch_2d, &
    average_down_reactive_eb_state_patch_2d, composite_eb_integral_2d
  use amr_eb_reactive_2d_mod, only: &
    reactive_eb_patch_exterior_context_2d, &
    prolong_reactive_eb_patch_pcm_2d, &
    extract_reactive_eb_patch_exterior_context_support_2d, &
    extract_reactive_eb_patch_exterior_context_2d, &
    build_reactive_eb_patch_exterior_from_context_2d, &
    build_reactive_eb_patch_exterior_2d, &
    advance_two_level_reactive_eb_hydro_2d
  implicit none

  integer, parameter :: coarse_nx = 8, coarse_ny = 8
  integer, parameter :: coarse_i_lower = 2, coarse_i_upper = 6
  integer, parameter :: coarse_j_lower = 2, coarse_j_upper = 6
  integer, parameter :: ratio = 2
  integer, parameter :: fine_nx = &
    (coarse_i_upper - coarse_i_lower + 1) * ratio
  integer, parameter :: fine_ny = &
    (coarse_j_upper - coarse_j_lower + 1) * ratio
  integer, parameter :: support_i_lower = max(1, coarse_i_lower - 1)
  integer, parameter :: support_i_upper = min(coarse_nx, coarse_i_upper + 1)
  integer, parameter :: support_j_lower = max(1, coarse_j_lower - 1)
  integer, parameter :: support_j_upper = min(coarse_ny, coarse_j_upper + 1)
  type(eb_geometry_2d) :: coarse_geometry, fine_geometry
  type(amr_eb_patch_2d) :: patch
  type(reactive_eb_exterior_state_2d) :: exterior, context_exterior
  type(reactive_eb_exterior_state_2d) :: support_exterior
  type(reactive_eb_patch_exterior_context_2d) :: exterior_context
  type(reactive_eb_patch_exterior_context_2d) :: support_context
  type(nasa7_species), allocatable :: species(:)
  real(dp) :: coarse_level_set(0:coarse_nx, 0:coarse_ny)
  real(dp) :: fine_level_set(0:fine_nx, 0:fine_ny)
  real(dp), allocatable :: primitive(:), state_cell(:), mass_fractions(:)
  real(dp), allocatable :: coarse_state(:, :, :), coarse_end(:, :, :)
  real(dp), allocatable :: coarse_start_support(:, :, :)
  real(dp), allocatable :: coarse_end_support(:, :, :)
  real(dp), allocatable :: fine_state(:, :, :), restricted_state(:, :, :)
  real(dp), allocatable :: new_coarse_state(:, :, :), new_fine_state(:, :, :)
  real(dp), allocatable :: coarse_temperature(:, :), coarse_end_temperature(:, :)
  real(dp), allocatable :: coarse_start_temperature_support(:, :)
  real(dp), allocatable :: coarse_end_temperature_support(:, :)
  real(dp), allocatable :: fine_temperature(:, :), restricted_temperature(:, :)
  real(dp), allocatable :: new_coarse_temperature(:, :)
  real(dp), allocatable :: new_fine_temperature(:, :)
  real(dp), allocatable :: integral_before(:), integral_after(:)
  real(dp) :: mole_fractions(7), x, y, fine_x_lower, fine_x_upper
  real(dp) :: fine_y_lower, fine_y_upper, temperature_cell, sound_speed
  real(dp) :: state_scale, integral_scale, dt, expected_scale
  logical :: ok, found_open_boundary
  integer :: i, j, k, nvar, component

  do j = 0, coarse_ny
    y = real(j, dp) / real(coarse_ny, dp)
    do i = 0, coarse_nx
      x = real(i, dp) / real(coarse_nx, dp)
      coarse_level_set(i, j) = x + y - 0.78_dp
    end do
  end do
  call build_eb_geometry_2d( &
    coarse_level_set, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    coarse_geometry, ok)
  call require(ok, "coarse EB geometry")

  fine_x_lower = real(coarse_i_lower - 1, dp) / real(coarse_nx, dp)
  fine_x_upper = real(coarse_i_upper, dp) / real(coarse_nx, dp)
  fine_y_lower = real(coarse_j_lower - 1, dp) / real(coarse_ny, dp)
  fine_y_upper = real(coarse_j_upper, dp) / real(coarse_ny, dp)
  do j = 0, fine_ny
    y = fine_y_lower + real(j, dp) * &
      (fine_y_upper - fine_y_lower) / real(fine_ny, dp)
    do i = 0, fine_nx
      x = fine_x_lower + real(i, dp) * &
        (fine_x_upper - fine_x_lower) / real(fine_nx, dp)
      fine_level_set(i, j) = x + y - 0.78_dp
    end do
  end do
  call build_eb_geometry_2d( &
    fine_level_set, fine_x_lower, fine_x_upper, &
    fine_y_lower, fine_y_upper, fine_geometry, ok)
  call require(ok, "fine EB geometry")
  call build_amr_eb_patch_2d( &
    coarse_geometry, fine_geometry, coarse_i_lower, coarse_i_upper, &
    coarse_j_lower, coarse_j_upper, ratio, patch, ok)
  call require(ok, "EB AMR patch")

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamic database load")
  nvar = reactive_nvar(size(species))
  allocate(primitive(reactive_nprim(size(species))))
  allocate(state_cell(nvar), mass_fractions(size(species)))
  mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.55643_dp]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  call require(ok, "composition conversion")
  primitive(1:5) = [0.31_dp, 0.0_dp, 0.0_dp, 0.0_dp, 135000.0_dp]
  do k = 1, size(species)
    primitive(reactive_mass_fraction_component(k)) = mass_fractions(k)
  end do
  call reactive_primitive_to_conserved( &
    species, primitive, state_cell, temperature_cell, sound_speed, ok)
  call require(ok, "reference reactive state")

  allocate(coarse_state(nvar, coarse_nx, coarse_ny))
  allocate(coarse_end(nvar, coarse_nx, coarse_ny))
  allocate(fine_state(nvar, fine_nx, fine_ny))
  allocate(restricted_state(nvar, coarse_nx, coarse_ny))
  allocate(new_coarse_state(nvar, coarse_nx, coarse_ny))
  allocate(new_fine_state(nvar, fine_nx, fine_ny))
  allocate(coarse_temperature(coarse_nx, coarse_ny))
  allocate(coarse_end_temperature(coarse_nx, coarse_ny))
  allocate(fine_temperature(fine_nx, fine_ny))
  allocate(restricted_temperature(coarse_nx, coarse_ny))
  allocate(new_coarse_temperature(coarse_nx, coarse_ny))
  allocate(new_fine_temperature(fine_nx, fine_ny))
  allocate(integral_before(nvar), integral_after(nvar))
  do j = 1, coarse_ny
    do i = 1, coarse_nx
      coarse_state(:, i, j) = state_cell
      coarse_temperature(i, j) = temperature_cell
    end do
  end do

  call prolong_reactive_eb_patch_pcm_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, &
    fine_geometry, patch, fine_state, fine_temperature, ok)
  call require(ok, "reactive PCM prolongation")
  state_scale = max(1.0_dp, maxval(abs(state_cell)))
  call require(maxval(abs(fine_state - &
    spread(spread(state_cell, 2, fine_nx), 3, fine_ny))) == 0.0_dp, &
    "PCM child-state injection")
  call require(maxval(abs(fine_temperature - temperature_cell)) <= &
    2.0e-8_dp, "PCM child temperature recovery")
  call average_down_reactive_eb_state_patch_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, fine_state, &
    fine_geometry, patch, restricted_state, restricted_temperature, ok)
  call require(ok .and. maxval(abs(restricted_state - coarse_state)) <= &
    5.0e-14_dp * state_scale, "prolong/restrict constant preservation")

  coarse_end = 1.02_dp * coarse_state
  coarse_end_temperature = coarse_temperature
  expected_scale = 1.005_dp
  call build_reactive_eb_patch_exterior_2d( &
    species, coarse_state, coarse_temperature, coarse_end, &
    coarse_end_temperature, coarse_geometry, fine_geometry, patch, &
    0.25_dp, exterior, ok)
  call require(ok .and. exterior%is_valid(fine_geometry, nvar), &
    "time-interpolated coarse exterior")
  call extract_reactive_eb_patch_exterior_context_2d( &
    coarse_state, coarse_temperature, coarse_end, coarse_end_temperature, &
    coarse_geometry, fine_geometry, patch, nvar, exterior_context, ok)
  call require(ok .and. exterior_context%is_valid(fine_geometry, nvar), &
    "compact coarse exterior context")
  allocate(coarse_start_support( &
    nvar, support_i_lower:support_i_upper, &
    support_j_lower:support_j_upper))
  allocate(coarse_end_support, mold=coarse_start_support)
  allocate(coarse_start_temperature_support( &
    support_i_lower:support_i_upper, support_j_lower:support_j_upper))
  allocate(coarse_end_temperature_support, &
    mold=coarse_start_temperature_support)
  coarse_start_support = coarse_state( &
    :, support_i_lower:support_i_upper, support_j_lower:support_j_upper)
  coarse_end_support = coarse_end( &
    :, support_i_lower:support_i_upper, support_j_lower:support_j_upper)
  coarse_start_temperature_support = coarse_temperature( &
    support_i_lower:support_i_upper, support_j_lower:support_j_upper)
  coarse_end_temperature_support = coarse_end_temperature( &
    support_i_lower:support_i_upper, support_j_lower:support_j_upper)
  call require(size(coarse_start_support) < size(coarse_state), &
    "exterior context uses compact coarse support")
  call extract_reactive_eb_patch_exterior_context_support_2d( &
    support_i_lower, support_j_lower, coarse_start_support, &
    coarse_start_temperature_support, coarse_end_support, &
    coarse_end_temperature_support, coarse_geometry, fine_geometry, patch, &
    nvar, support_context, ok)
  call require(ok .and. support_context%is_valid(fine_geometry, nvar), &
    "support coarse exterior context")
  call build_reactive_eb_patch_exterior_from_context_2d( &
    species, exterior_context, coarse_geometry, fine_geometry, patch, &
    0.25_dp, context_exterior, ok)
  call require(ok .and. &
    maxval(abs(context_exterior%x_lower_state - &
      exterior%x_lower_state)) == 0.0_dp .and. &
    maxval(abs(context_exterior%x_upper_state - &
      exterior%x_upper_state)) == 0.0_dp .and. &
    maxval(abs(context_exterior%y_lower_state - &
      exterior%y_lower_state)) == 0.0_dp .and. &
    maxval(abs(context_exterior%y_upper_state - &
      exterior%y_upper_state)) == 0.0_dp .and. &
    maxval(abs(context_exterior%x_lower_temperature - &
      exterior%x_lower_temperature)) == 0.0_dp .and. &
    maxval(abs(context_exterior%x_upper_temperature - &
      exterior%x_upper_temperature)) == 0.0_dp .and. &
    maxval(abs(context_exterior%y_lower_temperature - &
      exterior%y_lower_temperature)) == 0.0_dp .and. &
    maxval(abs(context_exterior%y_upper_temperature - &
      exterior%y_upper_temperature)) == 0.0_dp, &
    "compact exterior context parity")
  call build_reactive_eb_patch_exterior_from_context_2d( &
    species, support_context, coarse_geometry, fine_geometry, patch, &
    0.25_dp, support_exterior, ok)
  call require(ok .and. &
    maxval(abs(support_exterior%x_lower_state - &
      context_exterior%x_lower_state)) == 0.0_dp .and. &
    maxval(abs(support_exterior%x_upper_state - &
      context_exterior%x_upper_state)) == 0.0_dp .and. &
    maxval(abs(support_exterior%y_lower_state - &
      context_exterior%y_lower_state)) == 0.0_dp .and. &
    maxval(abs(support_exterior%y_upper_state - &
      context_exterior%y_upper_state)) == 0.0_dp .and. &
    maxval(abs(support_exterior%x_lower_temperature - &
      context_exterior%x_lower_temperature)) == 0.0_dp .and. &
    maxval(abs(support_exterior%x_upper_temperature - &
      context_exterior%x_upper_temperature)) == 0.0_dp .and. &
    maxval(abs(support_exterior%y_lower_temperature - &
      context_exterior%y_lower_temperature)) == 0.0_dp .and. &
    maxval(abs(support_exterior%y_upper_temperature - &
      context_exterior%y_upper_temperature)) == 0.0_dp, &
    "support/full exterior context parity")
  call extract_reactive_eb_patch_exterior_context_support_2d( &
    support_i_lower + 1, support_j_lower, &
    coarse_start_support(:, support_i_lower + 1:support_i_upper, :), &
    coarse_start_temperature_support(support_i_lower + 1:support_i_upper, :), &
    coarse_end_support(:, support_i_lower + 1:support_i_upper, :), &
    coarse_end_temperature_support(support_i_lower + 1:support_i_upper, :), &
    coarse_geometry, fine_geometry, patch, nvar, support_context, ok)
  call require(.not. ok, "incomplete exterior context support rejection")
  call extract_reactive_eb_patch_exterior_context_support_2d( &
    0, support_j_lower, coarse_start_support, &
    coarse_start_temperature_support, coarse_end_support, &
    coarse_end_temperature_support, coarse_geometry, fine_geometry, patch, &
    nvar, support_context, ok)
  call require(.not. ok, "out-of-root exterior context support rejection")
  coarse_start_support(1, support_i_lower, support_j_lower) = &
    ieee_value(0.0_dp, ieee_quiet_nan)
  call extract_reactive_eb_patch_exterior_context_support_2d( &
    support_i_lower, support_j_lower, coarse_start_support, &
    coarse_start_temperature_support, coarse_end_support, &
    coarse_end_temperature_support, coarse_geometry, fine_geometry, patch, &
    nvar, support_context, ok)
  call require(.not. ok, "nonfinite exterior context support rejection")
  coarse_start_support(1, support_i_lower, support_j_lower) = &
    coarse_state(1, support_i_lower, support_j_lower)
  found_open_boundary = .false.
  do j = 1, fine_ny
    if (fine_geometry%x_face_fraction(0, j) > 0.0_dp) then
      call require(maxval(abs(exterior%x_lower_state(:, j) - &
        expected_scale * state_cell)) <= 5.0e-14_dp * state_scale, &
        "lower-x exterior interpolation")
      found_open_boundary = .true.
    end if
    if (fine_geometry%x_face_fraction(fine_nx, j) > 0.0_dp) then
      call require(maxval(abs(exterior%x_upper_state(:, j) - &
        expected_scale * state_cell)) <= 5.0e-14_dp * state_scale, &
        "upper-x exterior interpolation")
      found_open_boundary = .true.
    end if
  end do
  do i = 1, fine_nx
    if (fine_geometry%y_face_fraction(i, 0) > 0.0_dp) then
      call require(maxval(abs(exterior%y_lower_state(:, i) - &
        expected_scale * state_cell)) <= 5.0e-14_dp * state_scale, &
        "lower-y exterior interpolation")
      found_open_boundary = .true.
    end if
    if (fine_geometry%y_face_fraction(i, fine_ny) > 0.0_dp) then
      call require(maxval(abs(exterior%y_upper_state(:, i) - &
        expected_scale * state_cell)) <= 5.0e-14_dp * state_scale, &
        "upper-y exterior interpolation")
      found_open_boundary = .true.
    end if
  end do
  call require(found_open_boundary, "open coarse/fine boundary coverage")
  call build_reactive_eb_patch_exterior_2d( &
    species, coarse_state, coarse_temperature, coarse_end, &
    coarse_end_temperature, coarse_geometry, fine_geometry, patch, &
    1.01_dp, exterior, ok)
  call require(.not. ok, "invalid exterior interpolation time rejection")

  call composite_eb_integral_2d( &
    coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
    integral_before, ok)
  call require(ok, "initial composite integral")
  dt = 0.1_dp * min(coarse_geometry%dx, coarse_geometry%dy) / sound_speed
  call advance_two_level_reactive_eb_hydro_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, patch, "hllc", &
    "characteristic_plm", "mc", 2, dt, new_coarse_state, &
    new_coarse_temperature, new_fine_state, new_fine_temperature, ok)
  call require(ok, "subcycled two-level reactive EB advance")
  call composite_eb_integral_2d( &
    new_coarse_state, coarse_geometry, new_fine_state, fine_geometry, &
    patch, integral_after, ok)
  call require(ok .and. maxval(abs(integral_after - integral_before)) <= &
    2.0e-12_dp * state_scale, "two-level composite conservation")
  call require(maxval(abs(new_coarse_state - coarse_state)) <= &
    3.0e-12_dp * state_scale .and. &
    maxval(abs(new_fine_state - fine_state)) <= &
    3.0e-12_dp * state_scale, "uniform state preservation")
  call require(maxval(abs(new_coarse_temperature - coarse_temperature)) <= &
    3.0e-8_dp .and. maxval(abs(new_fine_temperature - fine_temperature)) <= &
    3.0e-8_dp, "uniform temperature preservation")
  call assert_covered_unchanged( &
    coarse_state, new_coarse_state, coarse_geometry, &
    "coarse covered-state preservation")
  call assert_covered_unchanged( &
    fine_state, new_fine_state, fine_geometry, &
    "fine covered-state preservation")

  fine_state = 1.01_dp * &
    spread(spread(state_cell, 2, fine_nx), 3, fine_ny)
  fine_temperature = temperature_cell
  call composite_eb_integral_2d( &
    coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
    integral_before, ok)
  call require(ok, "mismatched initial composite integral")
  dt = 0.02_dp * min(coarse_geometry%dx, coarse_geometry%dy) / sound_speed
  call advance_two_level_reactive_eb_hydro_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, patch, "hllc", &
    "pcm", "mc", 2, dt, new_coarse_state, new_coarse_temperature, &
    new_fine_state, new_fine_temperature, ok, 0.5_dp)
  call require(ok, "nonmatching coarse/fine EB advance")
  call composite_eb_integral_2d( &
    new_coarse_state, coarse_geometry, new_fine_state, fine_geometry, &
    patch, integral_after, ok)
  integral_scale = max(1.0_dp, maxval(abs(integral_before)))
  call require(ok .and. &
    abs(integral_after(irho) - integral_before(irho)) <= &
      5.0e-11_dp * integral_scale .and. &
    abs(integral_after(iet) - integral_before(iet)) <= &
      5.0e-11_dp * integral_scale, &
    "nonmatching composite mass and energy conservation")
  do k = 1, size(species)
    component = reactive_species_component(k)
    call require(abs(integral_after(component) - &
      integral_before(component)) <= 5.0e-11_dp * integral_scale, &
      "nonmatching composite species conservation")
  end do
  call require(maxval(abs(new_fine_state - fine_state)) > &
    1.0e-12_dp * state_scale, "nonmatching interface evolves")

  call advance_two_level_reactive_eb_hydro_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, patch, "unknown", &
    "pcm", "mc", 0, dt, new_coarse_state, new_coarse_temperature, &
    new_fine_state, new_fine_temperature, ok)
  call require(.not. ok .and. &
    maxval(abs(new_coarse_state - coarse_state)) == 0.0_dp .and. &
    maxval(abs(new_fine_state - fine_state)) == 0.0_dp .and. &
    maxval(abs(new_coarse_temperature - coarse_temperature)) == 0.0_dp .and. &
    maxval(abs(new_fine_temperature - fine_temperature)) == 0.0_dp, &
    "two-level advance rollback")

  call advance_two_level_reactive_eb_hydro_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, patch, "hllc", &
    "pcm", "mc", 0, dt, new_coarse_state, new_coarse_temperature, &
    new_fine_state, new_fine_temperature, ok, 1.01_dp)
  call require(.not. ok .and. &
    maxval(abs(new_coarse_state - coarse_state)) == 0.0_dp .and. &
    maxval(abs(new_fine_state - fine_state)) == 0.0_dp, &
    "invalid StateRedist target rollback")

  write(*, '(a)') "test_amr_eb_reactive_2d: PASS"

contains

  subroutine assert_covered_unchanged( &
      original, candidate, geometry, message)
    real(dp), intent(in) :: original(:, :, :), candidate(:, :, :)
    type(eb_geometry_2d), intent(in) :: geometry
    character(len=*), intent(in) :: message
    integer :: cell_i, cell_j

    do cell_j = 1, geometry%ny
      do cell_i = 1, geometry%nx
        if (geometry%cell_type(cell_i, cell_j) /= eb_covered_cell) cycle
        call require(maxval(abs(candidate(:, cell_i, cell_j) - &
          original(:, cell_i, cell_j))) == 0.0_dp, message)
      end do
    end do
  end subroutine assert_covered_unchanged

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_amr_eb_reactive_2d
