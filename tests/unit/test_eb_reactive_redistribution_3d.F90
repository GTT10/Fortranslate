program test_eb_reactive_redistribution_3d
  use, intrinsic :: ieee_arithmetic, only: &
    ieee_value, ieee_quiet_nan
  use precision_mod, only: dp
  use state_indices_mod, only: irho
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d, eb_cut_cell_3d, &
    build_axis_plane_eb_geometry_3d
  use eb_reactive_redistribution_3d_mod, only: &
    reactive_eb_flux_redistribute_3d, &
    advance_reactive_eb_redistributed_3d, &
    reactive_eb_weighted_state_redistribute_3d, &
    advance_reactive_eb_state_redistributed_3d
  implicit none

  call check_regular_and_covered()
  call check_axis("x", 8, 1, 1, [3, 1, 1], [4, 1, 1])
  call check_axis("y", 1, 8, 1, [1, 3, 1], [1, 4, 1])
  call check_axis("z", 1, 1, 8, [1, 1, 3], [1, 1, 4])
  call check_state_axis("x", 8, 1, 1, [3, 1, 1], [4, 1, 1])
  call check_state_axis("y", 1, 8, 1, [1, 3, 1], [1, 4, 1])
  call check_state_axis("z", 1, 1, 8, [1, 1, 3], [1, 1, 4])
  call check_transverse_neighbors()
  call check_reactive_advance()

  write(*, '(a)') "test_eb_reactive_redistribution_3d: PASS"

contains

  subroutine check_regular_and_covered()
    integer, parameter :: nx = 3
    integer, parameter :: ny = 2
    integer, parameter :: nz = 2
    type(eb_geometry_3d) :: geometry
    real(dp) :: rhs(2, nx, ny, nz), redistributed(2, nx, ny, nz)
    logical :: ok
    integer :: i, j, k

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          rhs(:, i, j, k) = &
            [real(i + 3 * j + 7 * k, dp), real(2 * i - j + k, dp)]
        end do
      end do
    end do
    call build_axis_plane_eb_geometry_3d( &
      nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      "x", -0.1_dp, geometry, ok)
    call require(ok, "fully regular geometry")
    call reactive_eb_flux_redistribute_3d( &
      geometry, rhs, redistributed, ok)
    call require(ok .and. all(redistributed == rhs), &
      "fully regular redistribution identity")

    call build_axis_plane_eb_geometry_3d( &
      nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      "x", 1.1_dp, geometry, ok)
    call require(ok, "fully covered geometry")
    call reactive_eb_flux_redistribute_3d( &
      geometry, rhs, redistributed, ok)
    call require(ok .and. maxval(abs(redistributed)) == 0.0_dp, &
      "fully covered redistribution is zero")

    geometry%volume_fraction(1, 1, 1) = 2.0_dp
    call reactive_eb_flux_redistribute_3d( &
      geometry, rhs, redistributed, ok)
    call require(.not. ok .and. maxval(abs(redistributed)) == 0.0_dp, &
      "invalid geometry rejection")
  end subroutine check_regular_and_covered

  subroutine check_axis(axis, nx, ny, nz, cut, neighbor)
    character(len=*), intent(in) :: axis
    integer, intent(in) :: nx, ny, nz, cut(3), neighbor(3)

    integer, parameter :: ncomp = 5
    type(eb_geometry_3d) :: geometry
    real(dp), allocatable :: rhs(:, :, :, :), redistributed(:, :, :, :)
    real(dp) :: reference(ncomp), neighborhood_value
    real(dp) :: expected_cut, expected_neighbor
    real(dp) :: original_integral, redistributed_integral
    real(dp) :: kappa, tolerance
    logical :: ok
    integer :: component, i, j, k

    allocate(rhs(ncomp, nx, ny, nz), redistributed(ncomp, nx, ny, nz))
    call build_axis_plane_eb_geometry_3d( &
      nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      axis, 0.36875_dp, geometry, ok)
    call require(ok .and. geometry%is_valid(), &
      trim(axis)//"-axis small-cell geometry")
    call require(geometry%cell_type(cut(1), cut(2), cut(3)) == &
      eb_cut_cell_3d, trim(axis)//"-axis cut-cell location")
    kappa = geometry%volume_fraction(cut(1), cut(2), cut(3))
    call assert_close(kappa, 0.05_dp, 2.0e-13_dp, &
      trim(axis)//"-axis small volume fraction")

    reference = [1000.0_dp, -40.0_dp, 7.0_dp, 13.0_dp, -2800.0_dp]
    rhs = 0.0_dp
    rhs(:, cut(1), cut(2), cut(3)) = reference
    call reactive_eb_flux_redistribute_3d( &
      geometry, rhs, redistributed, ok)
    call require(ok, trim(axis)//"-axis flux redistribution")
    do component = 1, ncomp
      neighborhood_value = &
        kappa * reference(component) / (1.0_dp + kappa)
      expected_cut = kappa * reference(component) + &
        (1.0_dp - kappa) * neighborhood_value
      expected_neighbor = kappa * (1.0_dp - kappa) * &
        (reference(component) - neighborhood_value)
      tolerance = 3.0e-13_dp * max(1.0_dp, abs(reference(component)))
      call assert_close( &
        redistributed(component, cut(1), cut(2), cut(3)), &
        expected_cut, tolerance, trim(axis)//"-axis stabilized cut rhs")
      call assert_close( &
        redistributed(component, neighbor(1), neighbor(2), neighbor(3)), &
        expected_neighbor, tolerance, &
        trim(axis)//"-axis receiving-cell rhs")
      original_integral = sum( &
        geometry%volume_fraction * rhs(component, :, :, :))
      redistributed_integral = sum( &
        geometry%volume_fraction * redistributed(component, :, :, :))
      call assert_close(redistributed_integral, original_integral, &
        tolerance, trim(axis)//"-axis volume-weighted conservation")
    end do
    call require(abs(redistributed(1, cut(1), cut(2), cut(3))) < &
      0.11_dp * abs(rhs(1, cut(1), cut(2), cut(3))), &
      trim(axis)//"-axis small-cell stiffness reduction")
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (all([i, j, k] == cut) .or. &
              all([i, j, k] == neighbor)) cycle
          call require(maxval(abs(redistributed(:, i, j, k))) == 0.0_dp, &
            trim(axis)//"-axis compact support")
        end do
      end do
    end do

    rhs = 0.0_dp
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (geometry%cell_type(i, j, k) /= eb_covered_cell_3d) &
            rhs(:, i, j, k) = reference
        end do
      end do
    end do
    call reactive_eb_flux_redistribute_3d( &
      geometry, rhs, redistributed, ok)
    call require(ok, trim(axis)//"-axis uniform redistribution")
    call require(maxval(abs(redistributed - rhs)) <= &
      3.0e-13_dp * maxval(abs(rhs)), &
      trim(axis)//"-axis uniform active rhs preservation")

    if (axis == "x") then
      rhs = 0.0_dp
      rhs(1, cut(1), cut(2), cut(3)) = &
        ieee_value(0.0_dp, ieee_quiet_nan)
      call reactive_eb_flux_redistribute_3d( &
        geometry, rhs, redistributed, ok)
      call require(.not. ok .and. maxval(abs(redistributed)) == 0.0_dp, &
        "nonfinite rhs rejection")
    end if
  end subroutine check_axis

  subroutine check_state_axis(axis, nx, ny, nz, cut, neighbor)
    character(len=*), intent(in) :: axis
    integer, intent(in) :: nx, ny, nz, cut(3), neighbor(3)

    integer, parameter :: ncomp = 3
    type(eb_geometry_3d) :: geometry
    real(dp), allocatable :: provisional(:, :, :, :)
    real(dp), allocatable :: redistributed(:, :, :, :)
    real(dp) :: reference(ncomp), kappa, alpha
    real(dp) :: expected_cut, expected_neighbor
    real(dp) :: original_integral, redistributed_integral, tolerance
    logical :: ok
    integer :: component, i, j, k

    allocate(provisional(ncomp, nx, ny, nz))
    allocate(redistributed(ncomp, nx, ny, nz))
    call build_axis_plane_eb_geometry_3d( &
      nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      axis, 0.36875_dp, geometry, ok)
    call require(ok, trim(axis)//"-axis weighted-state geometry")
    kappa = geometry%volume_fraction(cut(1), cut(2), cut(3))
    alpha = (0.5_dp - kappa) / &
      geometry%volume_fraction(neighbor(1), neighbor(2), neighbor(3))
    expected_cut = (-kappa + 0.5_dp * alpha) / &
      (kappa + 0.5_dp * alpha)
    expected_neighbor = 1.0_dp - 0.5_dp * alpha + &
      0.5_dp * alpha * expected_cut

    reference = [0.31_dp, -1.7_dp, 2.4_dp]
    do component = 1, ncomp
      provisional(component, :, :, :) = reference(component)
    end do
    provisional(:, cut(1), cut(2), cut(3)) = -reference
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) &
            provisional(:, i, j, k) = 9.0_dp * reference
        end do
      end do
    end do
    call reactive_eb_weighted_state_redistribute_3d( &
      geometry, provisional, redistributed, ok)
    call require(ok, trim(axis)//"-axis weighted StateRedist")
    do component = 1, ncomp
      tolerance = 5.0e-13_dp * max(1.0_dp, abs(reference(component)))
      call assert_close( &
        redistributed(component, cut(1), cut(2), cut(3)), &
        expected_cut * reference(component), tolerance, &
        trim(axis)//"-axis weighted cut state")
      call assert_close( &
        redistributed(component, neighbor(1), neighbor(2), neighbor(3)), &
        expected_neighbor * reference(component), tolerance, &
        trim(axis)//"-axis weighted receiving state")
      original_integral = sum( &
        geometry%volume_fraction * provisional(component, :, :, :))
      redistributed_integral = sum( &
        geometry%volume_fraction * redistributed(component, :, :, :))
      call assert_close(redistributed_integral, original_integral, &
        tolerance, trim(axis)//"-axis weighted-state conservation")
    end do
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) then
            call require(maxval(abs(redistributed(:, i, j, k))) == &
              0.0_dp, trim(axis)//"-axis covered state is zero")
          else if (.not. all([i, j, k] == cut) .and. &
                   .not. all([i, j, k] == neighbor)) then
            call require(all(redistributed(:, i, j, k) == reference), &
              trim(axis)//"-axis weighted compact support")
          end if
        end do
      end do
    end do

    do component = 1, ncomp
      provisional(component, :, :, :) = reference(component)
    end do
    call reactive_eb_weighted_state_redistribute_3d( &
      geometry, provisional, redistributed, ok)
    call require(ok, trim(axis)//"-axis uniform weighted StateRedist")
    call require(maxval(abs(redistributed - provisional), &
      mask=spread(geometry%cell_type /= eb_covered_cell_3d, 1, ncomp)) <= &
      8.0e-13_dp, trim(axis)//"-axis uniform state preservation")
    call require(maxval(abs(redistributed), &
      mask=spread(geometry%cell_type == eb_covered_cell_3d, 1, ncomp)) == &
      0.0_dp, trim(axis)//"-axis uniform covered state is zero")

    call reactive_eb_weighted_state_redistribute_3d( &
      geometry, provisional, redistributed, ok, kappa)
    call require(ok .and. maxval(abs(redistributed - provisional), &
      mask=spread(geometry%cell_type /= eb_covered_cell_3d, 1, ncomp)) == &
      0.0_dp, trim(axis)//"-axis custom target identity")

    call reactive_eb_weighted_state_redistribute_3d( &
      geometry, provisional, redistributed, ok, 0.0_dp)
    call require(.not. ok .and. maxval(abs(redistributed)) == 0.0_dp, &
      trim(axis)//"-axis invalid weighted target transaction")
    if (axis == "x") then
      provisional(1, cut(1), cut(2), cut(3)) = &
        ieee_value(0.0_dp, ieee_quiet_nan)
      call reactive_eb_weighted_state_redistribute_3d( &
        geometry, provisional, redistributed, ok)
      call require(.not. ok .and. maxval(abs(redistributed)) == 0.0_dp, &
        "nonfinite weighted state transaction")
      call build_axis_plane_eb_geometry_3d( &
        nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
        axis, 0.99375_dp, geometry, ok)
      call require(ok, "boundary small-cell geometry")
      provisional = 1.0_dp
      call reactive_eb_weighted_state_redistribute_3d( &
        geometry, provisional, redistributed, ok)
      call require(.not. ok .and. maxval(abs(redistributed)) == 0.0_dp, &
        "missing regular receiver transaction")
    end if
  end subroutine check_state_axis

  subroutine check_transverse_neighbors()
    integer, parameter :: nx = 8
    integer, parameter :: ny = 3
    integer, parameter :: nz = 3
    type(eb_geometry_3d) :: geometry
    real(dp) :: rhs(1, nx, ny, nz), redistributed(1, nx, ny, nz)
    real(dp) :: original_integral, redistributed_integral
    logical :: ok

    call build_axis_plane_eb_geometry_3d( &
      nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      "x", 0.36875_dp, geometry, ok)
    call require(ok, "transverse-neighbor geometry")
    rhs = 0.0_dp
    rhs(1, 3, 2, 2) = 1.0_dp
    call reactive_eb_flux_redistribute_3d( &
      geometry, rhs, redistributed, ok)
    call require(ok, "transverse-neighbor redistribution")
    call require(redistributed(1, 4, 2, 2) > 0.0_dp, &
      "positive-x neighbor receives excess")
    call require(redistributed(1, 3, 1, 2) > 0.0_dp .and. &
      redistributed(1, 3, 3, 2) > 0.0_dp, &
      "both y neighbors receive excess")
    call require(redistributed(1, 3, 2, 1) > 0.0_dp .and. &
      redistributed(1, 3, 2, 3) > 0.0_dp, &
      "both z neighbors receive excess")
    original_integral = sum(geometry%volume_fraction * rhs(1, :, :, :))
    redistributed_integral = sum( &
      geometry%volume_fraction * redistributed(1, :, :, :))
    call assert_close(redistributed_integral, original_integral, &
      5.0e-14_dp, "overlapping-neighborhood conservation")
  end subroutine check_transverse_neighbors

  subroutine check_reactive_advance()
    integer, parameter :: nx = 8
    integer, parameter :: ny = 1
    integer, parameter :: nz = 1
    integer, parameter :: cut_i = 3
    integer, parameter :: neighbor_i = 4
    type(nasa7_species), allocatable :: species(:)
    type(eb_geometry_3d) :: geometry
    real(dp), allocatable :: primitive(:), state_cell(:)
    real(dp), allocatable :: mass_fractions(:)
    real(dp), allocatable :: state(:, :, :, :), rhs(:, :, :, :)
    real(dp), allocatable :: new_state(:, :, :, :)
    real(dp), allocatable :: temperature(:, :, :)
    real(dp), allocatable :: new_temperature(:, :, :)
    real(dp) :: mole_fractions(7), cell_temperature, sound_speed
    real(dp) :: kappa, alpha, expected_cut_scale, expected_neighbor_scale
    real(dp) :: original_integral, advanced_integral, tolerance
    logical :: ok
    integer :: component, nvar, species_index

    call load_h2o2_elementary_thermo(species, ok)
    call require(ok, "reactive thermodynamic database load")
    nvar = reactive_nvar(size(species))
    allocate(primitive(reactive_nprim(size(species))))
    allocate(state_cell(nvar), mass_fractions(size(species)))
    allocate(state(nvar, nx, ny, nz), rhs(nvar, nx, ny, nz))
    allocate(new_state(nvar, nx, ny, nz))
    allocate(temperature(nx, ny, nz), new_temperature(nx, ny, nz))

    mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
      1.0e-5_dp, 0.0_dp, 0.55643_dp]
    call mass_fractions_from_mole_fractions( &
      species, mole_fractions, mass_fractions, ok)
    call require(ok, "reactive composition conversion")
    primitive(1:5) = [0.31_dp, 0.0_dp, 0.0_dp, 0.0_dp, 135000.0_dp]
    do species_index = 1, size(species)
      primitive(reactive_mass_fraction_component(species_index)) = &
        mass_fractions(species_index)
    end do
    call reactive_primitive_to_conserved( &
      species, primitive, state_cell, cell_temperature, sound_speed, ok)
    call require(ok, "reactive reference-state construction")
    do component = 1, nvar
      state(component, :, :, :) = state_cell(component)
    end do
    temperature = cell_temperature

    call build_axis_plane_eb_geometry_3d( &
      nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      "x", 0.36875_dp, geometry, ok)
    call require(ok, "reactive small-cell geometry")
    kappa = geometry%volume_fraction(cut_i, 1, 1)
    rhs = 0.0_dp
    rhs(:, cut_i, 1, 1) = -2.0_dp * state_cell
    call require(state(irho, cut_i, 1, 1) + &
      rhs(irho, cut_i, 1, 1) < 0.0_dp, &
      "unredistributed update is nonphysical")
    call advance_reactive_eb_redistributed_3d( &
      species, state, temperature, geometry, rhs, 1.0_dp, &
      new_state, new_temperature, ok)
    call require(ok, "positive redistributed reactive advance")
    expected_cut_scale = 1.0_dp - 4.0_dp * kappa / (1.0_dp + kappa)
    expected_neighbor_scale = 1.0_dp - &
      2.0_dp * kappa * (1.0_dp - kappa) / (1.0_dp + kappa)
    call assert_close(new_state(irho, cut_i, 1, 1), &
      expected_cut_scale * state_cell(irho), 3.0e-13_dp, &
      "positive redistributed cut density")
    call assert_close(new_state(irho, neighbor_i, 1, 1), &
      expected_neighbor_scale * state_cell(irho), 3.0e-13_dp, &
      "positive redistributed neighbor density")
    call assert_close(new_temperature(cut_i, 1, 1), cell_temperature, &
      3.0e-9_dp, "redistributed cut temperature recovery")
    do component = 1, nvar
      original_integral = sum(geometry%volume_fraction * &
        rhs(component, :, :, :))
      advanced_integral = sum(geometry%volume_fraction * &
        (new_state(component, :, :, :) - state(component, :, :, :)))
      tolerance = 5.0e-13_dp * max(1.0_dp, abs(original_integral))
      call assert_close(advanced_integral, original_integral, tolerance, &
        "reactive advanced-state conservation")
    end do

    call advance_reactive_eb_state_redistributed_3d( &
      species, state, temperature, geometry, rhs, 1.0_dp, &
      new_state, new_temperature, ok)
    call require(ok, "positive weighted-state reactive advance")
    alpha = (0.5_dp - kappa) / &
      geometry%volume_fraction(neighbor_i, 1, 1)
    expected_cut_scale = (-kappa + 0.5_dp * alpha) / &
      (kappa + 0.5_dp * alpha)
    expected_neighbor_scale = 1.0_dp - 0.5_dp * alpha + &
      0.5_dp * alpha * expected_cut_scale
    call assert_close(new_state(irho, cut_i, 1, 1), &
      expected_cut_scale * state_cell(irho), 3.0e-13_dp, &
      "positive weighted-state cut density")
    call assert_close(new_state(irho, neighbor_i, 1, 1), &
      expected_neighbor_scale * state_cell(irho), 3.0e-13_dp, &
      "positive weighted-state neighbor density")
    call assert_close(new_temperature(cut_i, 1, 1), cell_temperature, &
      3.0e-9_dp, "weighted-state cut temperature recovery")
    do component = 1, nvar
      original_integral = sum(geometry%volume_fraction * &
        (state(component, :, :, :) + rhs(component, :, :, :)))
      advanced_integral = sum(geometry%volume_fraction * &
        new_state(component, :, :, :))
      tolerance = 5.0e-13_dp * max(1.0_dp, abs(original_integral))
      call assert_close(advanced_integral, original_integral, tolerance, &
        "weighted-state reactive advance conservation")
    end do

    call advance_reactive_eb_state_redistributed_3d( &
      species, state, temperature, geometry, rhs, 1.0_dp, &
      new_state, new_temperature, ok, kappa)
    call require(.not. ok .and. all(new_state == state) .and. &
      all(new_temperature == temperature), &
      "custom weighted target reactive rollback")

    rhs(:, cut_i, 1, 1) = -50.0_dp * state_cell
    call advance_reactive_eb_redistributed_3d( &
      species, state, temperature, geometry, rhs, 1.0_dp, &
      new_state, new_temperature, ok)
    call require(.not. ok .and. all(new_state == state) .and. &
      all(new_temperature == temperature), &
      "nonphysical redistributed advance rollback")
    call advance_reactive_eb_state_redistributed_3d( &
      species, state, temperature, geometry, rhs, 1.0_dp, &
      new_state, new_temperature, ok)
    call require(.not. ok .and. all(new_state == state) .and. &
      all(new_temperature == temperature), &
      "nonphysical weighted-state advance rollback")
  end subroutine check_reactive_advance

  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message

    if (abs(actual - expected) > tolerance) then
      write(*, '(a,2(1x,es24.16))') trim(message), actual, expected
      error stop message
    end if
  end subroutine assert_close

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_eb_reactive_redistribution_3d
