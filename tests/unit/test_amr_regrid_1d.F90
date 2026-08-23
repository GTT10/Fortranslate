program test_amr_regrid_1d
  use precision_mod, only: dp
  use amr_hierarchy_1d_mod, only: &
    amr_two_level_hierarchy_1d, initialize_two_level_hierarchy_1d, &
    prolong_conservative_1d, restrict_average_1d, composite_integral_1d
  use amr_regrid_1d_mod, only: &
    amr_tagging_criteria_1d, amr_regrid_plan_1d, &
    tag_gradient_1d, build_regrid_plan_1d, &
    plan_gradient_regrid_1d, regrid_two_level_state_1d
  implicit none

  integer, parameter :: variable_count = 2
  integer, parameter :: coarse_cells = 10
  integer, parameter :: refinement_ratio = 2
  real(dp), parameter :: tolerance = 3.0e-13_dp
  type(amr_tagging_criteria_1d) :: criteria
  type(amr_regrid_plan_1d) :: plan, inactive_plan
  type(amr_two_level_hierarchy_1d) :: old_hierarchy, new_hierarchy
  type(amr_two_level_hierarchy_1d) :: removed_hierarchy, created_hierarchy
  real(dp) :: state(variable_count, coarse_cells)
  real(dp) :: boundary_state(1, coarse_cells)
  real(dp) :: coarse(variable_count, coarse_cells)
  real(dp) :: restricted(variable_count, 3)
  real(dp) :: old_integral(variable_count), new_integral(variable_count)
  real(dp) :: uniform_integral(variable_count)
  real(dp) :: expected_overlap(variable_count, 4)
  real(dp), allocatable :: old_fine(:, :), new_fine(:, :)
  real(dp), allocatable :: removed_fine(:, :), created_fine(:, :)
  logical :: tags(coarse_cells), empty_tags(coarse_cells)
  logical :: boundary_tags(coarse_cells)
  logical :: ok
  integer :: cell

  state = 0.0_dp
  state(1, 6:coarse_cells) = 1.0_dp
  state(2, :) = 4.0_dp
  criteria%component = 1
  criteria%relative_gradient_threshold = 0.5_dp
  criteria%absolute_gradient_threshold = 0.1_dp
  criteria%scale_floor = 1.0e-12_dp
  criteria%buffer_cells = 1
  criteria%minimum_patch_cells = 4
  call plan_gradient_regrid_1d(state, criteria, tags, plan, ok)
  call assert_true(ok, "gradient regrid plan")
  call assert_true(count(tags) == 2 .and. tags(5) .and. tags(6), &
    "two cells tagged around discontinuity")
  call assert_true(plan%active, "discontinuity activates refinement")
  call assert_true(plan%tag_lower == 5 .and. plan%tag_upper == 6, &
    "tag bounds")
  call assert_true(plan%patch_lower == 4 .and. plan%patch_upper == 7, &
    "buffered patch bounds")

  criteria%component = 2
  call plan_gradient_regrid_1d(state, criteria, tags, plan, ok)
  call assert_true(ok .and. .not. any(tags), "flat component untagged")
  call assert_true(.not. plan%active, "empty tags remove refinement")
  criteria%component = 3
  call tag_gradient_1d(state, criteria, tags, ok)
  call assert_true(.not. ok, "invalid tagging component rejected")

  boundary_state = 0.0_dp
  boundary_state(1, 1) = 1.0_dp
  criteria%component = 1
  call tag_gradient_1d(boundary_state, criteria, boundary_tags, ok)
  call assert_true(ok .and. boundary_tags(1), "boundary feature tagged")
  call build_regrid_plan_1d(boundary_tags, 0, 1, plan, ok)
  call assert_true(ok .and. plan%active, "boundary patch accepted")
  call assert_true(plan%tag_lower == 1 .and. plan%patch_lower == 1, &
    "boundary plan retains physical edge")

  empty_tags = .false.
  tags = .false.
  tags(5) = .true.
  call build_regrid_plan_1d(tags, 0, 4, plan, ok)
  call assert_true(ok, "minimum-width patch plan")
  call assert_true(plan%patch_lower == 3 .and. plan%patch_upper == 6, &
    "minimum patch expansion")

  call initialize_two_level_hierarchy_1d( &
    coarse_cells, 4, 7, refinement_ratio, 0.0_dp, 1.0_dp, &
    old_hierarchy, ok)
  call assert_true(ok, "old hierarchy")
  do cell = 1, coarse_cells
    coarse(1, cell) = real(cell, dp)
    coarse(2, cell) = 10.0_dp + 2.0_dp * real(cell, dp)
  end do
  allocate(old_fine(variable_count, old_hierarchy%fine%cell_count()))
  call prolong_conservative_1d(coarse, old_hierarchy, old_fine, ok)
  call assert_true(ok, "old fine initialization")
  old_fine(1, 1:2) = old_fine(1, 1:2) + 0.5_dp
  old_fine(1, 5) = old_fine(1, 5) - 0.2_dp
  old_fine(1, 6) = old_fine(1, 6) + 0.2_dp
  old_fine(2, 7) = old_fine(2, 7) - 0.3_dp
  old_fine(2, 8) = old_fine(2, 8) + 0.3_dp
  expected_overlap = old_fine(:, 5:8)
  call composite_integral_1d( &
    coarse, old_fine, old_hierarchy, old_integral, ok)
  call assert_true(ok, "old composite integral")

  tags = .false.
  tags(6:8) = .true.
  call build_regrid_plan_1d(tags, 0, 3, plan, ok)
  call assert_true(ok, "shifted patch plan")
  call regrid_two_level_state_1d( &
    coarse, old_hierarchy, old_fine, plan, refinement_ratio, &
    0.0_dp, 1.0_dp, new_hierarchy, new_fine, ok)
  call assert_true(ok, "dynamic patch regrid")
  call assert_true(new_hierarchy%fine_coarse_lower == 6 .and. &
    new_hierarchy%fine_coarse_upper == 8, "shifted hierarchy bounds")
  call assert_close(coarse(1, 4), 4.5_dp, tolerance, &
    "departing fine cell averaged down")
  call assert_close(maxval(abs(new_fine(:, 1:4) - expected_overlap)), &
    0.0_dp, tolerance, "overlapping fine data retained")
  call restrict_average_1d(new_fine, new_hierarchy, restricted, ok)
  call assert_true(ok, "new fine restriction")
  call assert_close(maxval(abs(restricted - coarse(:, 6:8))), &
    0.0_dp, tolerance, "new patch preserves coarse averages")
  call composite_integral_1d( &
    coarse, new_fine, new_hierarchy, new_integral, ok)
  call assert_true(ok, "new composite integral")
  call assert_close(maxval(abs(new_integral - old_integral)), &
    0.0_dp, tolerance, "regrid preserves composite integral")

  call build_regrid_plan_1d(empty_tags, 1, 2, inactive_plan, ok)
  call assert_true(ok .and. .not. inactive_plan%active, &
    "inactive removal plan")
  call regrid_two_level_state_1d( &
    coarse, new_hierarchy, new_fine, inactive_plan, refinement_ratio, &
    0.0_dp, 1.0_dp, removed_hierarchy, removed_fine, ok)
  call assert_true(ok, "refinement removal")
  call assert_true(.not. removed_hierarchy%is_valid(), &
    "removed hierarchy inactive")
  call assert_true(.not. allocated(removed_fine), &
    "removed fine storage released")
  uniform_integral = sum(coarse, dim=2) / real(coarse_cells, dp)
  call assert_close(maxval(abs(uniform_integral - new_integral)), &
    0.0_dp, tolerance, "removal preserves composite integral")

  tags = .false.
  tags(3:5) = .true.
  call build_regrid_plan_1d(tags, 0, 3, plan, ok)
  call assert_true(ok, "new refinement plan")
  call regrid_two_level_state_1d( &
    coarse, removed_hierarchy, removed_fine, plan, refinement_ratio, &
    0.0_dp, 1.0_dp, created_hierarchy, created_fine, ok)
  call assert_true(ok .and. allocated(created_fine), &
    "refinement creation")
  call composite_integral_1d( &
    coarse, created_fine, created_hierarchy, new_integral, ok)
  call assert_true(ok, "created composite integral")
  call assert_close(maxval(abs(new_integral - uniform_integral)), &
    0.0_dp, tolerance, "creation preserves uniform-grid integral")

  write(*, '(a)') "test_amr_regrid_1d: PASS"

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

end program test_amr_regrid_1d
