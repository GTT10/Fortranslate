program test_eb_reactive_redistribution_2d
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use precision_mod, only: dp
  use state_indices_mod, only: irho
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved
  use eb_geometry_2d_mod, only: &
    eb_geometry_2d, eb_covered_cell, eb_cut_cell, build_eb_geometry_2d
  use eb_reactive_redistribution_2d_mod, only: &
    reactive_eb_flux_redistribute_2d, &
    advance_reactive_eb_redistributed_2d, &
    reactive_eb_weighted_state_redistribute_2d, &
    advance_reactive_eb_state_redistributed_2d
  implicit none

  integer, parameter :: nx = 8
  integer, parameter :: ny = 1
  integer, parameter :: cut_i = 3
  integer, parameter :: neighbor_i = 4
  integer, parameter :: overlap_nx = 4
  integer, parameter :: overlap_ny = 4
  type(nasa7_species), allocatable :: species(:)
  type(eb_geometry_2d) :: geometry, overlap_geometry
  real(dp), allocatable :: primitive(:), state_cell(:), mass_fractions(:)
  real(dp), allocatable :: conservative_rhs(:, :, :)
  real(dp), allocatable :: redistributed_rhs(:, :, :)
  real(dp), allocatable :: state(:, :, :), temperature(:, :)
  real(dp), allocatable :: new_state(:, :, :), new_temperature(:, :)
  real(dp), allocatable :: provisional_state(:, :, :)
  real(dp), allocatable :: redistributed_state(:, :, :)
  real(dp), allocatable :: scalar_state(:, :, :), scalar_redistributed(:, :, :)
  real(dp), allocatable :: scalar_zeroth_order(:, :, :)
  real(dp) :: level_set(0:nx, 0:ny), mole_fractions(7)
  real(dp) :: overlap_level_set(0:overlap_nx, 0:overlap_ny)
  real(dp) :: reference_rhs(5), temperature_cell, sound_speed
  real(dp) :: kappa, neighborhood_value, expected_cut, expected_neighbor
  real(dp) :: original_integral, redistributed_integral, tolerance
  real(dp) :: alpha, large_kappa, neighborhood_a, neighborhood_b
  real(dp) :: qhat_a, qhat_b, expected_shared
  real(dp) :: cell_coordinate_x, cell_coordinate_y, linear_error
  logical :: ok
  integer :: i, j, k, nvar

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamic database load")
  nvar = reactive_nvar(size(species))
  allocate(primitive(reactive_nprim(size(species))), state_cell(nvar))
  allocate(mass_fractions(size(species)))
  allocate(conservative_rhs(nvar, nx, ny), redistributed_rhs(nvar, nx, ny))
  allocate(state(nvar, nx, ny), new_state(nvar, nx, ny))
  allocate(temperature(nx, ny), new_temperature(nx, ny))
  allocate(provisional_state(nvar, nx, ny))
  allocate(redistributed_state(nvar, nx, ny))

  do j = 0, ny
    do i = 0, nx
      level_set(i, j) = real(i, dp) / real(nx, dp) - 0.36875_dp
    end do
  end do
  call build_eb_geometry_2d( &
    level_set, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, geometry, ok)
  call require(ok .and. geometry%is_valid(), "small-cell geometry")
  call require(geometry%cell_type(cut_i, 1) == eb_cut_cell, &
    "small cut-cell location")
  call require(all(geometry%cell_type(1:cut_i - 1, 1) == &
    eb_covered_cell), "covered cells")
  kappa = geometry%volume_fraction(cut_i, 1)
  call assert_close(kappa, 0.05_dp, 2.0e-13_dp, "small volume fraction")

  conservative_rhs = 0.0_dp
  reference_rhs = [1000.0_dp, -40.0_dp, 7.0_dp, 13.0_dp, -2800.0_dp]
  conservative_rhs(1:5, cut_i, 1) = reference_rhs
  call reactive_eb_flux_redistribute_2d( &
    geometry, conservative_rhs, redistributed_rhs, ok)
  call require(ok, "small-cell flux redistribution")
  do k = 1, 5
    neighborhood_value = kappa * reference_rhs(k) / (1.0_dp + kappa)
    expected_cut = kappa * reference_rhs(k) + &
      (1.0_dp - kappa) * neighborhood_value
    expected_neighbor = kappa * (1.0_dp - kappa) * &
      (reference_rhs(k) - neighborhood_value)
    tolerance = 2.0e-13_dp * max(1.0_dp, abs(reference_rhs(k)))
    call assert_close(redistributed_rhs(k, cut_i, 1), expected_cut, &
      tolerance, "stabilized cut-cell rhs")
    call assert_close(redistributed_rhs(k, neighbor_i, 1), &
      expected_neighbor, tolerance, "redistributed neighbor rhs")
    original_integral = sum( &
      geometry%volume_fraction * conservative_rhs(k, :, :))
    redistributed_integral = sum( &
      geometry%volume_fraction * redistributed_rhs(k, :, :))
    call assert_close(redistributed_integral, original_integral, &
      tolerance, "volume-weighted rhs conservation")
  end do
  call require(maxval(abs(redistributed_rhs(:, 1:cut_i - 1, :))) == &
    0.0_dp, "covered-cell rhs remains zero")
  call require(maxval(abs(redistributed_rhs(:, neighbor_i + 1:nx, :))) == &
    0.0_dp, "redistribution has compact support")
  call require(abs(redistributed_rhs(1, cut_i, 1)) < &
    0.11_dp * abs(conservative_rhs(1, cut_i, 1)), &
    "small-cell stiffness reduction")

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
  call require(ok, "reference state construction")
  do j = 1, ny
    do i = 1, nx
      state(:, i, j) = state_cell
      temperature(i, j) = temperature_cell
    end do
  end do

  conservative_rhs = 0.0_dp
  do i = cut_i, nx
    conservative_rhs(:, i, 1) = 0.25_dp * state_cell
  end do
  call reactive_eb_flux_redistribute_2d( &
    geometry, conservative_rhs, redistributed_rhs, ok)
  call require(ok, "uniform active rhs redistribution")
  call require(maxval(abs(redistributed_rhs(:, cut_i:nx, :) - &
    conservative_rhs(:, cut_i:nx, :))) <= &
    2.0e-13_dp * maxval(abs(conservative_rhs)), &
    "uniform active rhs preservation")

  conservative_rhs = 0.0_dp
  conservative_rhs(:, cut_i, 1) = -2.0_dp * state_cell
  call advance_reactive_eb_redistributed_2d( &
    species, state, temperature, geometry, conservative_rhs, 1.0_dp, &
    new_state, new_temperature, ok)
  call require(ok, "positive small-cell redistributed advance")
  expected_cut = 1.0_dp - 4.0_dp * kappa / (1.0_dp + kappa)
  expected_neighbor = 1.0_dp - &
    2.0_dp * kappa * (1.0_dp - kappa) / (1.0_dp + kappa)
  call require(expected_cut > 0.0_dp, "stabilized cut-cell scale")
  call assert_close(new_state(irho, cut_i, 1), &
    expected_cut * state_cell(irho), 2.0e-13_dp, &
    "positive cut-cell density")
  call assert_close(new_state(irho, neighbor_i, 1), &
    expected_neighbor * state_cell(irho), 2.0e-13_dp, &
    "positive neighbor density")
  call assert_close(new_temperature(cut_i, 1), temperature_cell, &
    2.0e-9_dp, "cut-cell temperature recovery")
  do k = 1, nvar
    original_integral = sum(geometry%volume_fraction * &
      (new_state(k, :, :) - state(k, :, :)))
    redistributed_integral = sum(geometry%volume_fraction * &
      conservative_rhs(k, :, :))
    tolerance = 4.0e-13_dp * max(1.0_dp, abs(redistributed_integral))
    call assert_close(original_integral, redistributed_integral, tolerance, &
      "advanced state conservation")
  end do

  conservative_rhs(:, cut_i, 1) = -50.0_dp * state_cell
  call advance_reactive_eb_redistributed_2d( &
    species, state, temperature, geometry, conservative_rhs, 1.0_dp, &
    new_state, new_temperature, ok)
  call require(.not. ok .and. &
    maxval(abs(new_state - state)) == 0.0_dp .and. &
    maxval(abs(new_temperature - temperature)) == 0.0_dp, &
    "nonphysical advance transaction")

  provisional_state = state
  provisional_state(:, cut_i, 1) = -state_cell
  call reactive_eb_weighted_state_redistribute_2d( &
    geometry, provisional_state, redistributed_state, ok)
  call require(ok, "weighted state redistribution")
  alpha = (0.5_dp - kappa) / &
    geometry%volume_fraction(neighbor_i, 1)
  neighborhood_value = (-kappa + 0.5_dp * alpha) / &
    (kappa + 0.5_dp * alpha)
  expected_cut = neighborhood_value
  expected_neighbor = 1.0_dp - 0.5_dp * alpha + &
    0.5_dp * alpha * neighborhood_value
  call assert_close(redistributed_state(irho, cut_i, 1), &
    expected_cut * state_cell(irho), 3.0e-13_dp, &
    "weighted cut-cell state")
  call assert_close(redistributed_state(irho, neighbor_i, 1), &
    expected_neighbor * state_cell(irho), 3.0e-13_dp, &
    "weighted receiving state")
  do k = 1, nvar
    original_integral = sum(geometry%volume_fraction * &
      provisional_state(k, :, :))
    redistributed_integral = sum(geometry%volume_fraction * &
      redistributed_state(k, :, :))
    tolerance = 5.0e-13_dp * max(1.0_dp, abs(original_integral))
    call assert_close(redistributed_integral, original_integral, tolerance, &
      "weighted state component conservation")
  end do

  conservative_rhs = 0.0_dp
  conservative_rhs(:, cut_i, 1) = -2.0_dp * state_cell
  call advance_reactive_eb_state_redistributed_2d( &
    species, state, temperature, geometry, conservative_rhs, 1.0_dp, &
    new_state, new_temperature, ok)
  call require(ok, "weighted reactive state advance")
  call assert_close(new_state(irho, cut_i, 1), &
    expected_cut * state_cell(irho), 3.0e-13_dp, &
    "weighted positive cut-cell density")
  call assert_close(new_state(irho, neighbor_i, 1), &
    expected_neighbor * state_cell(irho), 3.0e-13_dp, &
    "weighted positive receiving density")
  call assert_close(new_temperature(cut_i, 1), temperature_cell, &
    2.0e-9_dp, "weighted cut-cell temperature recovery")
  do k = 1, nvar
    original_integral = sum(geometry%volume_fraction * &
      (state(k, :, :) + conservative_rhs(k, :, :)))
    redistributed_integral = sum(geometry%volume_fraction * &
      new_state(k, :, :))
    tolerance = 5.0e-13_dp * max(1.0_dp, abs(original_integral))
    call assert_close(redistributed_integral, original_integral, tolerance, &
      "weighted reactive advance conservation")
  end do

  call advance_reactive_eb_state_redistributed_2d( &
    species, state, temperature, geometry, conservative_rhs, 1.0_dp, &
    new_state, new_temperature, ok, kappa)
  call require(.not. ok .and. &
    maxval(abs(new_state - state)) == 0.0_dp .and. &
    maxval(abs(new_temperature - temperature)) == 0.0_dp, &
    "custom weighted target forwarding")

  conservative_rhs(:, cut_i, 1) = -50.0_dp * state_cell
  call advance_reactive_eb_state_redistributed_2d( &
    species, state, temperature, geometry, conservative_rhs, 1.0_dp, &
    new_state, new_temperature, ok)
  call require(.not. ok .and. &
    maxval(abs(new_state - state)) == 0.0_dp .and. &
    maxval(abs(new_temperature - temperature)) == 0.0_dp, &
    "weighted nonphysical advance transaction")

  do j = 0, overlap_ny
    do i = 0, overlap_nx
      overlap_level_set(i, j) = &
        real(i + j, dp) / real(overlap_nx, dp) - &
        (1.0_dp - 0.25_dp * sqrt(0.1_dp))
    end do
  end do
  call build_eb_geometry_2d( &
    overlap_level_set, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    overlap_geometry, ok)
  call require(ok .and. overlap_geometry%is_valid(), &
    "overlapping neighborhood geometry")
  call assert_close(overlap_geometry%volume_fraction(1, 3), &
    0.05_dp, 3.0e-13_dp, "overlap small-cell fraction")
  call assert_close(overlap_geometry%volume_fraction(2, 2), &
    0.05_dp, 3.0e-13_dp, "second overlap small-cell fraction")
  large_kappa = overlap_geometry%volume_fraction(2, 3)
  call assert_close(large_kappa, 0.45_dp + sqrt(0.1_dp), &
    3.0e-13_dp, "overlap receiving fraction")
  allocate(scalar_state(1, overlap_nx, overlap_ny))
  allocate(scalar_redistributed(1, overlap_nx, overlap_ny))
  allocate(scalar_zeroth_order(1, overlap_nx, overlap_ny))
  scalar_state = 0.0_dp
  scalar_state(1, 2, 3) = 1.0_dp
  call reactive_eb_weighted_state_redistribute_2d( &
    overlap_geometry, scalar_state, scalar_redistributed, ok)
  call require(ok, "overlapping weighted neighborhoods")
  alpha = 0.45_dp / (2.0_dp * large_kappa + 1.0_dp)
  neighborhood_a = 0.05_dp + alpha * &
    (large_kappa / 2.0_dp + large_kappa / 3.0_dp + 0.5_dp)
  neighborhood_b = 0.05_dp + alpha * &
    (large_kappa / 3.0_dp + large_kappa / 3.0_dp + 0.5_dp)
  qhat_a = alpha * large_kappa / (3.0_dp * neighborhood_a)
  qhat_b = alpha * large_kappa / (3.0_dp * neighborhood_b)
  expected_shared = 1.0_dp - 2.0_dp * alpha / 3.0_dp + &
    alpha * (qhat_a + qhat_b) / 3.0_dp
  call assert_close(scalar_redistributed(1, 1, 3), qhat_a, &
    5.0e-13_dp, "first overlapping neighborhood average")
  call assert_close(scalar_redistributed(1, 2, 2), qhat_b, &
    5.0e-13_dp, "second overlapping neighborhood average")
  call assert_close(scalar_redistributed(1, 2, 3), expected_shared, &
    5.0e-13_dp, "shared-cell neighborhood accumulation")
  original_integral = sum(overlap_geometry%volume_fraction * &
    scalar_state(1, :, :))
  redistributed_integral = sum(overlap_geometry%volume_fraction * &
    scalar_redistributed(1, :, :))
  call assert_close(redistributed_integral, original_integral, &
    8.0e-13_dp, "overlapping neighborhood conservation")

  scalar_state = 2.5_dp
  call reactive_eb_weighted_state_redistribute_2d( &
    overlap_geometry, scalar_state, scalar_redistributed, ok)
  call require(ok, "uniform weighted state redistribution")
  call require(maxval(abs(scalar_redistributed(1, :, :) - 2.5_dp), &
    mask=overlap_geometry%cell_type /= eb_covered_cell) <= 8.0e-13_dp, &
    "uniform weighted state preservation")
  call require(maxval(abs(scalar_redistributed(1, :, :)), &
    mask=overlap_geometry%cell_type == eb_covered_cell) == 0.0_dp, &
    "weighted covered-cell state remains zero")

  do j = 1, overlap_ny
    do i = 1, overlap_nx
      cell_coordinate_x = real(i, dp) - 0.5_dp + &
        overlap_geometry%cell_centroid_x(i, j)
      cell_coordinate_y = real(j, dp) - 0.5_dp + &
        overlap_geometry%cell_centroid_y(i, j)
      scalar_state(1, i, j) = 2.0_dp + 0.17_dp * cell_coordinate_x - &
        0.11_dp * cell_coordinate_y
    end do
  end do
  call reactive_eb_weighted_state_redistribute_2d( &
    overlap_geometry, scalar_state, scalar_zeroth_order, ok, 0.5_dp, 0)
  call require(ok, "zeroth-order affine StateRedist")
  call reactive_eb_weighted_state_redistribute_2d( &
    overlap_geometry, scalar_state, scalar_redistributed, ok, 0.5_dp, 2)
  call require(ok, "second-order affine StateRedist")
  linear_error = maxval(abs(scalar_redistributed - scalar_state), &
    mask=spread(overlap_geometry%cell_type /= eb_covered_cell, 1, 1))
  call require(linear_error <= 2.0e-12_dp, &
    "second-order affine state preservation")
  call require(maxval(abs(scalar_zeroth_order - scalar_state), &
    mask=spread(overlap_geometry%cell_type /= eb_covered_cell, 1, 1)) > &
    1.0e-5_dp, "zeroth-order affine diffusion")
  original_integral = sum(overlap_geometry%volume_fraction * &
    scalar_state(1, :, :))
  redistributed_integral = sum(overlap_geometry%volume_fraction * &
    scalar_redistributed(1, :, :))
  call assert_close(redistributed_integral, original_integral, &
    2.0e-12_dp, "second-order affine conservation")

  scalar_state = 0.0_dp
  do j = 1, overlap_ny
    do i = 1, overlap_nx
      if (overlap_geometry%cell_type(i, j) /= eb_covered_cell .and. &
          i + j >= 6) scalar_state(1, i, j) = 1.0_dp
    end do
  end do
  call reactive_eb_weighted_state_redistribute_2d( &
    overlap_geometry, scalar_state, scalar_redistributed, ok, 0.5_dp, 2)
  call require(ok, "limited second-order discontinuous StateRedist")
  call require(minval(scalar_redistributed, &
    mask=spread(overlap_geometry%cell_type /= eb_covered_cell, 1, 1)) >= &
    -2.0e-12_dp .and. maxval(scalar_redistributed, &
    mask=spread(overlap_geometry%cell_type /= eb_covered_cell, 1, 1)) <= &
    1.0_dp + 2.0e-12_dp, "second-order StateRedist monotonicity")
  original_integral = sum(overlap_geometry%volume_fraction * &
    scalar_state(1, :, :))
  redistributed_integral = sum(overlap_geometry%volume_fraction * &
    scalar_redistributed(1, :, :))
  call assert_close(redistributed_integral, original_integral, &
    2.0e-12_dp, "limited second-order conservation")

  call reactive_eb_weighted_state_redistribute_2d( &
    overlap_geometry, scalar_state, scalar_redistributed, ok, 0.0_dp)
  call require(.not. ok .and. maxval(abs(scalar_redistributed)) == 0.0_dp, &
    "invalid weighted target transaction")
  call reactive_eb_weighted_state_redistribute_2d( &
    overlap_geometry, scalar_state, scalar_redistributed, ok, 0.5_dp, 1)
  call require(.not. ok .and. maxval(abs(scalar_redistributed)) == 0.0_dp, &
    "invalid StateRedist max-order transaction")

  conservative_rhs = 0.0_dp
  call advance_reactive_eb_state_redistributed_2d( &
    species, state, temperature, geometry, conservative_rhs, 1.0_dp, &
    new_state, new_temperature, ok, 0.5_dp, 1)
  call require(.not. ok .and. &
    maxval(abs(new_state - state)) == 0.0_dp .and. &
    maxval(abs(new_temperature - temperature)) == 0.0_dp, &
    "invalid max-order reactive advance rollback")

  conservative_rhs = 0.0_dp
  conservative_rhs(1, cut_i, 1) = &
    ieee_value(0.0_dp, ieee_quiet_nan)
  call reactive_eb_flux_redistribute_2d( &
    geometry, conservative_rhs, redistributed_rhs, ok)
  call require(.not. ok .and. maxval(abs(redistributed_rhs)) == 0.0_dp, &
    "nonfinite rhs transaction")
  provisional_state = state
  provisional_state(1, cut_i, 1) = &
    ieee_value(0.0_dp, ieee_quiet_nan)
  call reactive_eb_weighted_state_redistribute_2d( &
    geometry, provisional_state, redistributed_state, ok)
  call require(.not. ok .and. maxval(abs(redistributed_state)) == 0.0_dp, &
    "nonfinite weighted state transaction")

  write(*, '(a)') "test_eb_reactive_redistribution_2d: PASS"

contains

  subroutine assert_close(actual, expected, local_tolerance, message)
    real(dp), intent(in) :: actual, expected, local_tolerance
    character(len=*), intent(in) :: message

    call require(abs(actual - expected) <= local_tolerance, message)
  end subroutine assert_close

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_eb_reactive_redistribution_2d
