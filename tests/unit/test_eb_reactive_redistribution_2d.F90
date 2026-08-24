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
    advance_reactive_eb_redistributed_2d
  implicit none

  integer, parameter :: nx = 8
  integer, parameter :: ny = 1
  integer, parameter :: cut_i = 3
  integer, parameter :: neighbor_i = 4
  type(nasa7_species), allocatable :: species(:)
  type(eb_geometry_2d) :: geometry
  real(dp), allocatable :: primitive(:), state_cell(:), mass_fractions(:)
  real(dp), allocatable :: conservative_rhs(:, :, :)
  real(dp), allocatable :: redistributed_rhs(:, :, :)
  real(dp), allocatable :: state(:, :, :), temperature(:, :)
  real(dp), allocatable :: new_state(:, :, :), new_temperature(:, :)
  real(dp) :: level_set(0:nx, 0:ny), mole_fractions(7)
  real(dp) :: reference_rhs(5), temperature_cell, sound_speed
  real(dp) :: kappa, neighborhood_value, expected_cut, expected_neighbor
  real(dp) :: original_integral, redistributed_integral, tolerance
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

  conservative_rhs = 0.0_dp
  conservative_rhs(1, cut_i, 1) = &
    ieee_value(0.0_dp, ieee_quiet_nan)
  call reactive_eb_flux_redistribute_2d( &
    geometry, conservative_rhs, redistributed_rhs, ok)
  call require(.not. ok .and. maxval(abs(redistributed_rhs)) == 0.0_dp, &
    "nonfinite rhs transaction")

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
