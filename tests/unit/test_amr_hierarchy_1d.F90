program test_amr_hierarchy_1d
  use precision_mod, only: dp
  use amr_hierarchy_1d_mod, only: &
    amr_two_level_hierarchy_1d, amr_flux_register_1d, &
    initialize_two_level_hierarchy_1d, prolong_conservative_1d, &
    restrict_average_1d, average_down_1d, &
    level_subcycle_time_steps_1d, initialize_flux_register_1d, &
    accumulate_coarse_flux_1d, accumulate_fine_flux_1d, reflux_1d, &
    composite_integral_1d
  implicit none

  integer, parameter :: variable_count = 2
  integer, parameter :: coarse_cells = 8
  integer, parameter :: fine_cells = 8
  real(dp), parameter :: tolerance = 2.0e-13_dp
  type(amr_two_level_hierarchy_1d) :: hierarchy, invalid_hierarchy
  type(amr_flux_register_1d) :: flux_register
  real(dp) :: coarse(variable_count, coarse_cells)
  real(dp) :: averaged(variable_count, coarse_cells)
  real(dp) :: fine(variable_count, fine_cells)
  real(dp) :: restricted(variable_count, 4)
  real(dp) :: expected_fine(variable_count, fine_cells)
  real(dp) :: fine_time_steps(2)
  real(dp) :: composite_coarse(variable_count, coarse_cells)
  real(dp) :: composite_fine(variable_count, fine_cells)
  real(dp) :: composite_before(variable_count)
  real(dp) :: composite_after(variable_count)
  real(dp) :: coarse_left_flux(variable_count)
  real(dp) :: coarse_right_flux(variable_count)
  real(dp) :: fine_left_flux(variable_count)
  real(dp) :: fine_right_flux(variable_count)
  real(dp) :: fine_integral(variable_count)
  real(dp) :: offset, coarse_time_step
  logical :: ok
  integer :: cell, child, fine_cell

  call initialize_two_level_hierarchy_1d( &
    coarse_cells, 3, 6, 2, 0.0_dp, 1.0_dp, hierarchy, ok)
  call assert_true(ok, "valid hierarchy")
  call assert_true(hierarchy%coarse%cell_count() == coarse_cells, &
    "coarse cell count")
  call assert_true(hierarchy%fine%cell_count() == fine_cells, &
    "fine cell count")
  call assert_true(hierarchy%covered_coarse_cells() == 4, &
    "covered coarse count")
  call assert_close(hierarchy%coarse_dx, 0.125_dp, tolerance, &
    "coarse spacing")
  call assert_close(hierarchy%fine_dx, 0.0625_dp, tolerance, &
    "fine spacing")

  call initialize_two_level_hierarchy_1d( &
    coarse_cells, 1, 6, 2, 0.0_dp, 1.0_dp, invalid_hierarchy, ok)
  call assert_true(.not. ok, "boundary-touching patch rejection")
  call initialize_two_level_hierarchy_1d( &
    coarse_cells, 3, 6, 1, 0.0_dp, 1.0_dp, invalid_hierarchy, ok)
  call assert_true(.not. ok, "invalid refinement ratio rejection")

  do cell = 1, coarse_cells
    coarse(1, cell) = real(cell, dp)
    coarse(2, cell) = 3.0_dp
  end do
  call prolong_conservative_1d(coarse, hierarchy, fine, ok)
  call assert_true(ok, "conservative prolongation")
  expected_fine = 0.0_dp
  do cell = hierarchy%fine_coarse_lower, hierarchy%fine_coarse_upper
    do child = 1, hierarchy%refinement_ratio
      fine_cell = (cell - hierarchy%fine_coarse_lower) * &
        hierarchy%refinement_ratio + child
      offset = (real(child, dp) - 0.5_dp) / &
        real(hierarchy%refinement_ratio, dp) - 0.5_dp
      expected_fine(1, fine_cell) = real(cell, dp) + offset
      expected_fine(2, fine_cell) = 3.0_dp
    end do
  end do
  call assert_close(maxval(abs(fine - expected_fine)), 0.0_dp, &
    tolerance, "limited linear prolongation")

  call restrict_average_1d(fine, hierarchy, restricted, ok)
  call assert_true(ok, "restriction")
  call assert_close( &
    maxval(abs(restricted - coarse(:, 3:6))), 0.0_dp, tolerance, &
    "prolongation restriction identity")
  averaged = -1.0_dp
  call average_down_1d(fine, hierarchy, averaged, ok)
  call assert_true(ok, "average down")
  call assert_close( &
    maxval(abs(averaged(:, 3:6) - coarse(:, 3:6))), &
    0.0_dp, tolerance, "covered-cell synchronization")
  call assert_close(maxval(abs(averaged(:, 1:2) + 1.0_dp)), &
    0.0_dp, tolerance, "uncovered lower cells unchanged")
  call assert_close(maxval(abs(averaged(:, 7:8) + 1.0_dp)), &
    0.0_dp, tolerance, "uncovered upper cells unchanged")

  coarse_time_step = 0.1_dp
  call level_subcycle_time_steps_1d( &
    hierarchy, coarse_time_step, fine_time_steps, ok)
  call assert_true(ok, "subcycle schedule")
  call assert_close(sum(fine_time_steps), coarse_time_step, tolerance, &
    "subcycle time closure")
  call assert_close(maxval(abs(fine_time_steps - 0.05_dp)), &
    0.0_dp, tolerance, "subcycle step size")

  call initialize_flux_register_1d( &
    flux_register, variable_count, ok)
  call assert_true(ok, "flux-register allocation")
  coarse_left_flux = [1.0_dp, 2.0_dp]
  coarse_right_flux = [3.0_dp, 4.0_dp]
  fine_left_flux = [1.5_dp, 1.0_dp]
  fine_right_flux = [2.5_dp, 5.0_dp]
  call accumulate_coarse_flux_1d( &
    flux_register, coarse_left_flux, coarse_right_flux, &
    coarse_time_step, ok)
  call assert_true(ok, "coarse flux accumulation")
  call accumulate_fine_flux_1d( &
    flux_register, fine_left_flux, fine_right_flux, &
    fine_time_steps(1), ok)
  call assert_true(ok, "first fine flux accumulation")
  call accumulate_fine_flux_1d( &
    flux_register, fine_left_flux, fine_right_flux, &
    fine_time_steps(2), ok)
  call assert_true(ok, "second fine flux accumulation")

  composite_coarse = 0.0_dp
  composite_coarse(:, hierarchy%fine_coarse_lower - 1) = &
    -coarse_time_step * coarse_left_flux / hierarchy%coarse_dx
  composite_coarse(:, hierarchy%fine_coarse_upper + 1) = &
    coarse_time_step * coarse_right_flux / hierarchy%coarse_dx
  fine_integral = &
    coarse_time_step * (fine_left_flux - fine_right_flux)
  do cell = 1, fine_cells
    composite_fine(:, cell) = fine_integral / &
      (real(fine_cells, dp) * hierarchy%fine_dx)
  end do
  call composite_integral_1d( &
    composite_coarse, composite_fine, hierarchy, composite_before, ok)
  call assert_true(ok, "composite integral before reflux")
  call assert_true(maxval(abs(composite_before)) > 1.0e-3_dp, &
    "nontrivial coarse-fine flux mismatch")
  call reflux_1d(composite_coarse, hierarchy, flux_register, ok)
  call assert_true(ok, "reflux")
  call composite_integral_1d( &
    composite_coarse, composite_fine, hierarchy, composite_after, ok)
  call assert_true(ok, "composite integral after reflux")
  call assert_close(maxval(abs(composite_after)), 0.0_dp, tolerance, &
    "reflux restores composite conservation")
  call assert_close(maxval(abs(flux_register%left)), 0.0_dp, tolerance, &
    "left register reset")
  call assert_close(maxval(abs(flux_register%right)), 0.0_dp, tolerance, &
    "right register reset")

  write(*, '(a)') "test_amr_hierarchy_1d: PASS"

contains

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*, '(a,1x,a)') "FAIL:", trim(label)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected, tol, label)
    real(dp), intent(in) :: actual, expected, tol
    character(len=*), intent(in) :: label
    if (abs(actual - expected) > tol) then
      write(*, '(a,1x,a,2(1x,es24.16))') &
        "FAIL:", trim(label), actual, expected
      error stop 1
    end if
  end subroutine assert_close

end program test_amr_hierarchy_1d
