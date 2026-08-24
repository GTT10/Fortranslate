program test_amr_eb_hierarchy_2d
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved
  use eb_geometry_2d_mod, only: &
    eb_geometry_2d, eb_covered_cell, build_eb_geometry_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, build_amr_eb_patch_2d, &
    average_down_eb_state_patch_2d, &
    average_down_reactive_eb_state_patch_2d, composite_eb_integral_2d
  implicit none

  integer, parameter :: coarse_nx = 8
  integer, parameter :: coarse_ny = 8
  integer, parameter :: coarse_i_lower = 2
  integer, parameter :: coarse_i_upper = 6
  integer, parameter :: coarse_j_lower = 2
  integer, parameter :: coarse_j_upper = 6
  integer, parameter :: ratio = 2
  integer, parameter :: fine_nx = &
    (coarse_i_upper - coarse_i_lower + 1) * ratio
  integer, parameter :: fine_ny = &
    (coarse_j_upper - coarse_j_lower + 1) * ratio
  type(eb_geometry_2d) :: coarse_geometry, fine_geometry
  type(amr_eb_patch_2d) :: patch, invalid_patch
  type(nasa7_species), allocatable :: species(:)
  real(dp) :: coarse_level_set(0:coarse_nx, 0:coarse_ny)
  real(dp) :: fine_level_set(0:fine_nx, 0:fine_ny)
  real(dp), allocatable :: coarse_state(:, :, :), fine_state(:, :, :)
  real(dp), allocatable :: averaged_state(:, :, :), constant_state(:, :, :)
  real(dp), allocatable :: composite_integral(:), coarse_integral(:)
  real(dp), allocatable :: reactive_coarse(:, :, :), reactive_fine(:, :, :)
  real(dp), allocatable :: reactive_averaged(:, :, :), state_cell(:)
  real(dp), allocatable :: coarse_temperature(:, :), averaged_temperature(:, :)
  real(dp), allocatable :: primitive(:), mass_fractions(:)
  real(dp) :: mole_fractions(7), x, y, fine_x_lower, fine_x_upper
  real(dp) :: fine_y_lower, fine_y_upper, temperature_cell, sound_speed
  real(dp) :: error, scale
  logical :: ok
  integer :: i, j, nvar

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
  call require(ok .and. coarse_geometry%is_valid(), "coarse EB geometry")

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
  call require(ok .and. fine_geometry%is_valid(), "fine EB geometry")
  call build_amr_eb_patch_2d( &
    coarse_geometry, fine_geometry, coarse_i_lower, coarse_i_upper, &
    coarse_j_lower, coarse_j_upper, ratio, patch, ok)
  call require(ok .and. patch%is_valid(coarse_geometry, fine_geometry), &
    "aligned EB AMR patch")
  call require(patch%coarse_cell_count_x() == 5 .and. &
    patch%coarse_cell_count_y() == 5, "patch extents")
  call require(any(coarse_geometry%cell_type( &
    coarse_i_lower:coarse_i_upper, coarse_j_lower:coarse_j_upper) == &
    eb_covered_cell), "patch contains covered cells")

  allocate(coarse_state(1, coarse_nx, coarse_ny))
  allocate(fine_state(1, fine_nx, fine_ny))
  allocate(averaged_state(1, coarse_nx, coarse_ny))
  allocate(constant_state(1, coarse_nx, coarse_ny))
  allocate(composite_integral(1), coarse_integral(1))
  do j = 1, coarse_ny
    do i = 1, coarse_nx
      x = coarse_geometry%x_lower + &
        (real(i, dp) - 0.5_dp + coarse_geometry%cell_centroid_x(i, j)) * &
        coarse_geometry%dx
      y = coarse_geometry%y_lower + &
        (real(j, dp) - 0.5_dp + coarse_geometry%cell_centroid_y(i, j)) * &
        coarse_geometry%dy
      coarse_state(1, i, j) = 2.0_dp + 0.17_dp * x - 0.11_dp * y
    end do
  end do
  do j = 1, fine_ny
    do i = 1, fine_nx
      x = fine_geometry%x_lower + &
        (real(i, dp) - 0.5_dp + fine_geometry%cell_centroid_x(i, j)) * &
        fine_geometry%dx
      y = fine_geometry%y_lower + &
        (real(j, dp) - 0.5_dp + fine_geometry%cell_centroid_y(i, j)) * &
        fine_geometry%dy
      fine_state(1, i, j) = 2.0_dp + 0.17_dp * x - 0.11_dp * y
    end do
  end do
  call average_down_eb_state_patch_2d( &
    coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
    averaged_state, ok)
  call require(ok, "EB volume-weighted average down")
  error = maxval(abs(averaged_state - coarse_state), &
    mask=spread(coarse_geometry%cell_type /= eb_covered_cell, 1, 1))
  call require(error <= 5.0e-13_dp, "affine fluid-centroid restriction")
  call require(maxval(abs(averaged_state(:, 1, :) - &
    coarse_state(:, 1, :))) == 0.0_dp .and. &
    maxval(abs(averaged_state(:, 7:8, :) - &
      coarse_state(:, 7:8, :))) == 0.0_dp, "outside-patch state unchanged")

  fine_state = 2.5_dp
  call average_down_eb_state_patch_2d( &
    coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
    constant_state, ok)
  call require(ok, "constant EB average down")
  call require(maxval(abs(constant_state(1, &
    coarse_i_lower:coarse_i_upper, coarse_j_lower:coarse_j_upper) - &
    2.5_dp)) == 0.0_dp, "constant patch preservation")

  do j = 1, fine_ny
    do i = 1, fine_nx
      x = fine_geometry%x_lower + &
        (real(i, dp) - 0.5_dp + fine_geometry%cell_centroid_x(i, j)) * &
        fine_geometry%dx
      y = fine_geometry%y_lower + &
        (real(j, dp) - 0.5_dp + fine_geometry%cell_centroid_y(i, j)) * &
        fine_geometry%dy
      fine_state(1, i, j) = 2.0_dp + 0.17_dp * x - 0.11_dp * y
    end do
  end do
  call composite_eb_integral_2d( &
    coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
    composite_integral, ok)
  call require(ok, "EB composite integral")
  call average_down_eb_state_patch_2d( &
    coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
    averaged_state, ok)
  call require(ok, "composite restriction")
  coarse_integral(1) = sum(coarse_geometry%volume_fraction * &
    averaged_state(1, :, :)) * coarse_geometry%dx * coarse_geometry%dy
  call assert_close(coarse_integral(1), composite_integral(1), &
    8.0e-13_dp, "restriction preserves composite integral")

  fine_state(1, fine_nx, fine_ny) = &
    ieee_value(0.0_dp, ieee_quiet_nan)
  call average_down_eb_state_patch_2d( &
    coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
    averaged_state, ok)
  call require(.not. ok .and. &
    maxval(abs(averaged_state - coarse_state)) == 0.0_dp, &
    "nonfinite average-down transaction")
  fine_state(1, fine_nx, fine_ny) = coarse_state(1, coarse_i_upper, &
    coarse_j_upper)

  call build_amr_eb_patch_2d( &
    coarse_geometry, fine_geometry, 1, coarse_i_upper, &
    coarse_j_lower, coarse_j_upper, ratio, invalid_patch, ok)
  call require(.not. ok, "misaligned fine-patch rejection")

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamic database load")
  nvar = reactive_nvar(size(species))
  allocate(primitive(reactive_nprim(size(species))))
  allocate(mass_fractions(size(species)), state_cell(nvar))
  mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.55643_dp]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  call require(ok, "composition conversion")
  primitive(1:5) = [0.31_dp, 4.0_dp, -2.0_dp, 0.0_dp, 135000.0_dp]
  do i = 1, size(species)
    primitive(reactive_mass_fraction_component(i)) = mass_fractions(i)
  end do
  call reactive_primitive_to_conserved( &
    species, primitive, state_cell, temperature_cell, sound_speed, ok)
  call require(ok, "reference reactive state")
  allocate(reactive_coarse(nvar, coarse_nx, coarse_ny))
  allocate(reactive_fine(nvar, fine_nx, fine_ny))
  allocate(reactive_averaged(nvar, coarse_nx, coarse_ny))
  allocate(coarse_temperature(coarse_nx, coarse_ny))
  allocate(averaged_temperature(coarse_nx, coarse_ny))
  do j = 1, coarse_ny
    do i = 1, coarse_nx
      reactive_coarse(:, i, j) = state_cell
      coarse_temperature(i, j) = temperature_cell
    end do
  end do
  do j = 1, fine_ny
    do i = 1, fine_nx
      reactive_fine(:, i, j) = state_cell
    end do
  end do
  call average_down_reactive_eb_state_patch_2d( &
    species, reactive_coarse, coarse_temperature, coarse_geometry, &
    reactive_fine, fine_geometry, patch, reactive_averaged, &
    averaged_temperature, ok)
  call require(ok, "reactive EB average-down transaction")
  scale = max(1.0_dp, maxval(abs(state_cell)))
  call require(maxval(abs(reactive_averaged - reactive_coarse)) <= &
    5.0e-14_dp * scale, "uniform reactive state preservation")
  call require(maxval(abs(averaged_temperature - coarse_temperature)) <= &
    2.0e-9_dp, "reactive average-down temperature recovery")

  reactive_fine(:, 9:10, 9:10) = -spread(spread(state_cell, 2, 2), 3, 2)
  call average_down_reactive_eb_state_patch_2d( &
    species, reactive_coarse, coarse_temperature, coarse_geometry, &
    reactive_fine, fine_geometry, patch, reactive_averaged, &
    averaged_temperature, ok)
  call require(.not. ok .and. &
    maxval(abs(reactive_averaged - reactive_coarse)) == 0.0_dp .and. &
    maxval(abs(averaged_temperature - coarse_temperature)) == 0.0_dp, &
    "nonphysical reactive average-down rollback")

  write(*, '(a)') "test_amr_eb_hierarchy_2d: PASS"

contains

  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message

    call require(abs(actual - expected) <= &
      tolerance * max(1.0_dp, abs(expected)), message)
  end subroutine assert_close

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_amr_eb_hierarchy_2d
