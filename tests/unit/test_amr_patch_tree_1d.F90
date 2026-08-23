program test_amr_patch_tree_1d
  use precision_mod, only: dp
  use amr_hierarchy_1d_mod, only: &
    amr_two_level_hierarchy_1d, restrict_average_1d
  use amr_patch_tree_1d_mod, only: &
    amr_patch_level_plan_1d, amr_patch_tree_hierarchy_1d, &
    amr_patch_tree_level_fields_1d, initialize_patch_tree_1d, &
    prolong_patch_tree_1d, average_down_patch_tree_1d, &
    composite_integral_patch_tree_1d, patch_tree_child_geometry_1d
  implicit none

  integer, parameter :: variable_count = 2
  integer, parameter :: base_cells = 16
  real(dp), parameter :: tolerance = 5.0e-12_dp
  type(amr_patch_level_plan_1d), allocatable :: plans(:), invalid_plans(:)
  type(amr_patch_tree_hierarchy_1d) :: hierarchy, invalid_hierarchy
  type(amr_patch_tree_level_fields_1d), allocatable :: fields(:)
  type(amr_two_level_hierarchy_1d) :: geometry
  real(dp) :: root(variable_count, base_cells)
  real(dp) :: root_integral(variable_count)
  real(dp) :: composite_before(variable_count)
  real(dp) :: composite_after(variable_count)
  real(dp) :: x, child_lower, child_upper
  logical :: ok
  integer :: cell, last_cell

  call configure_plans(plans)
  call initialize_patch_tree_1d( &
    base_cells, 0.0_dp, 1.0_dp, plans, hierarchy, ok)
  call assert_true(ok .and. hierarchy%is_valid(), &
    "four-level patch-tree initialization")
  call assert_true(hierarchy%level_count() == 4, &
    "patch-tree level count")
  call assert_true(all([ &
    hierarchy%level_patch_count(0), hierarchy%level_patch_count(1), &
    hierarchy%level_patch_count(2), hierarchy%level_patch_count(3)] == &
    [1, 2, 3, 2]), "branched patch counts")
  call assert_true(all([ &
    hierarchy%relations(1)%parent_patch_count(), &
    hierarchy%relations(2)%parent_patch_count(), &
    hierarchy%relations(3)%parent_patch_count()] == [1, 2, 3]), &
    "relation parent counts")
  call assert_true( &
    hierarchy%relations(2)%child_index(1, 2) == 2 .and. &
    hierarchy%relations(2)%child_index(2, 1) == 3, &
    "flattened child indexing")
  call assert_true( &
    hierarchy%relations(3)%child_sets(2)%patch_count() == 0, &
    "parent without deeper children")
  call assert_close(hierarchy%level_dx(0), 1.0_dp / 16.0_dp, &
    tolerance, "root spacing")
  call assert_close(hierarchy%level_dx(3), 1.0_dp / 192.0_dp, &
    tolerance, "mixed-ratio deepest spacing")

  call patch_tree_child_geometry_1d( &
    hierarchy%relations(2), 3, geometry, ok)
  call assert_true(ok, "branched child geometry lookup")
  child_lower = geometry%x_lower + &
    real(geometry%fine_coarse_lower - 1, dp) * geometry%coarse_dx
  child_upper = geometry%x_lower + &
    real(geometry%fine_coarse_upper, dp) * geometry%coarse_dx
  call assert_close(child_lower, 0.625_dp, tolerance, &
    "branched child physical lower bound")
  call assert_close(child_upper, 0.78125_dp, tolerance, &
    "branched child physical upper bound")

  do cell = 1, base_cells
    x = (real(cell, dp) - 0.5_dp) / real(base_cells, dp)
    root(1, cell) = 2.0_dp + 0.25_dp * x
    root(2, cell) = -1.0_dp + 0.50_dp * x
  end do
  root_integral = hierarchy%level_dx(0) * sum(root, dim=2)
  call prolong_patch_tree_1d(root, hierarchy, fields, ok)
  call assert_true(ok, "recursive patch-tree prolongation")
  call composite_integral_patch_tree_1d( &
    fields, hierarchy, composite_before, ok)
  call assert_true(ok, "initial patch-tree composite integral")
  call assert_close(maxval(abs(composite_before - root_integral)), &
    0.0_dp, tolerance, "recursive prolongation conservation")
  call assert_synchronized(fields, "initial patch-tree synchronization")

  fields(4)%patches(1)%values(1, 1) = &
    fields(4)%patches(1)%values(1, 1) + 0.75_dp
  last_cell = size(fields(4)%patches(2)%values, 2)
  fields(4)%patches(2)%values(2, last_cell) = &
    fields(4)%patches(2)%values(2, last_cell) - 0.40_dp
  call composite_integral_patch_tree_1d( &
    fields, hierarchy, composite_before, ok)
  call assert_true(ok, "perturbed patch-tree composite integral")
  call average_down_patch_tree_1d(fields, hierarchy, ok)
  call assert_true(ok, "deepest-to-root patch-tree average down")
  call assert_synchronized(fields, "recursive patch-tree average down")
  call composite_integral_patch_tree_1d( &
    fields, hierarchy, composite_after, ok)
  call assert_true(ok, "averaged patch-tree composite integral")
  call assert_close(maxval(abs(composite_after - composite_before)), &
    0.0_dp, tolerance, "patch-tree average-down conservation")

  invalid_plans = plans
  invalid_plans(2)%patches(3)%parent_patch = 3
  call initialize_patch_tree_1d( &
    base_cells, 0.0_dp, 1.0_dp, invalid_plans, invalid_hierarchy, ok)
  call assert_true(.not. ok, "invalid parent ownership rejected")

  write(*, '(a)') "test_amr_patch_tree_1d: PASS"

contains

  subroutine configure_plans(local_plans)
    type(amr_patch_level_plan_1d), allocatable, intent(out) :: local_plans(:)

    allocate(local_plans(3))

    local_plans(1)%refinement_ratio = 2
    allocate(local_plans(1)%patches(2))
    local_plans(1)%patches(1)%parent_patch = 1
    local_plans(1)%patches(1)%lower = 2
    local_plans(1)%patches(1)%upper = 6
    local_plans(1)%patches(2)%parent_patch = 1
    local_plans(1)%patches(2)%lower = 10
    local_plans(1)%patches(2)%upper = 14

    local_plans(2)%refinement_ratio = 2
    allocate(local_plans(2)%patches(3))
    local_plans(2)%patches(1)%parent_patch = 1
    local_plans(2)%patches(1)%lower = 2
    local_plans(2)%patches(1)%upper = 4
    local_plans(2)%patches(2)%parent_patch = 1
    local_plans(2)%patches(2)%lower = 7
    local_plans(2)%patches(2)%upper = 9
    local_plans(2)%patches(3)%parent_patch = 2
    local_plans(2)%patches(3)%lower = 3
    local_plans(2)%patches(3)%upper = 7

    local_plans(3)%refinement_ratio = 3
    allocate(local_plans(3)%patches(2))
    local_plans(3)%patches(1)%parent_patch = 1
    local_plans(3)%patches(1)%lower = 2
    local_plans(3)%patches(1)%upper = 5
    local_plans(3)%patches(2)%parent_patch = 3
    local_plans(3)%patches(2)%lower = 3
    local_plans(3)%patches(2)%upper = 8
  end subroutine configure_plans

  subroutine assert_synchronized(local_fields, label)
    type(amr_patch_tree_level_fields_1d), intent(in) :: local_fields(:)
    character(len=*), intent(in) :: label

    real(dp), allocatable :: restricted(:, :)
    logical :: local_ok
    integer :: relation, parent, child, index, covered_cells

    do relation = 1, size(hierarchy%relations)
      do parent = 1, hierarchy%relations(relation)%parent_patch_count()
        do child = 1, &
            hierarchy%relations(relation)%child_sets(parent)%patch_count()
          geometry = &
            hierarchy%relations(relation)%child_sets(parent)%patches(child)
          index = hierarchy%relations(relation)%child_index(parent, child)
          covered_cells = geometry%covered_coarse_cells()
          allocate(restricted(variable_count, covered_cells))
          call restrict_average_1d( &
            local_fields(relation + 1)%patches(index)%values, geometry, &
            restricted, local_ok)
          call assert_true(local_ok, trim(label) // " restriction")
          call assert_close(maxval(abs(restricted - &
            local_fields(relation)%patches(parent)%values(:, &
              geometry%fine_coarse_lower:geometry%fine_coarse_upper))), &
            0.0_dp, tolerance, label)
          deallocate(restricted)
        end do
      end do
    end do
  end subroutine assert_synchronized

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

end program test_amr_patch_tree_1d
