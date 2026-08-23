program test_amr_multilevel_hierarchy_1d
  use precision_mod, only: dp
  use amr_hierarchy_1d_mod, only: &
    amr_multilevel_hierarchy_1d, amr_level_field_1d, &
    amr_flux_register_1d, initialize_multilevel_hierarchy_1d, &
    prolong_multilevel_1d, restrict_average_1d, &
    average_down_multilevel_1d, multilevel_subcycle_counts_1d, &
    multilevel_subcycle_time_steps_1d, &
    initialize_multilevel_flux_registers_1d, &
    accumulate_coarse_flux_1d, accumulate_fine_flux_1d, &
    synchronize_multilevel_1d, composite_integral_multilevel_1d
  implicit none

  integer, parameter :: variable_count = 2
  integer, parameter :: base_cells = 16
  integer, parameter :: interface_count = 3
  integer, parameter :: level_count = interface_count + 1
  real(dp), parameter :: tolerance = 5.0e-12_dp
  integer :: patch_lower(interface_count), patch_upper(interface_count)
  integer :: refinement_ratios(interface_count), expected_cells(level_count)
  integer :: subcycle_counts(level_count), bad_lower(interface_count)
  integer, allocatable :: empty(:)
  type(amr_multilevel_hierarchy_1d) :: hierarchy, invalid_hierarchy
  type(amr_multilevel_hierarchy_1d) :: base_hierarchy
  type(amr_level_field_1d), allocatable :: fields(:), averaged_fields(:)
  type(amr_level_field_1d), allocatable :: flux_fields(:)
  type(amr_flux_register_1d), allocatable :: registers(:)
  real(dp) :: root(variable_count, base_cells)
  real(dp), allocatable :: restricted(:, :)
  real(dp) :: time_steps(level_count), expected_dx(level_count)
  real(dp) :: uniform_integral(variable_count)
  real(dp) :: integral_before(variable_count), integral_after(variable_count)
  real(dp) :: coarse_left(variable_count), coarse_right(variable_count)
  real(dp) :: fine_left(variable_count), fine_right(variable_count)
  real(dp) :: fine_density(variable_count)
  real(dp) :: x_lower, x_upper, x, coarse_time_step
  logical :: ok
  integer :: relation, level, cell, covered_cells

  patch_lower = [3, 4, 5]
  patch_upper = [14, 21, 50]
  refinement_ratios = [2, 3, 2]
  expected_cells = [16, 24, 54, 92]
  expected_dx = [1.0_dp / 16.0_dp, 1.0_dp / 32.0_dp, &
    1.0_dp / 96.0_dp, 1.0_dp / 192.0_dp]

  call initialize_multilevel_hierarchy_1d( &
    base_cells, patch_lower, patch_upper, refinement_ratios, &
    0.0_dp, 1.0_dp, hierarchy, ok)
  call assert_true(ok, "four-level hierarchy initialization")
  call assert_true(hierarchy%is_valid(), "four-level hierarchy validity")
  call assert_true(hierarchy%level_count() == level_count, &
    "arbitrary level count")
  do level = 0, level_count - 1
    call assert_true(hierarchy%level_cell_count(level) == &
      expected_cells(level + 1), "level cell count")
    call assert_close(hierarchy%level_dx(level), expected_dx(level + 1), &
      tolerance, "cumulative level spacing")
    call hierarchy%level_bounds(level, x_lower, x_upper, ok)
    call assert_true(ok .and. x_upper > x_lower, "physical level bounds")
    if (level > 0) then
      call assert_true(hierarchy%interfaces(level)%coarse%level == &
        level - 1, "parent level numbering")
      call assert_true(hierarchy%interfaces(level)%fine%level == level, &
        "child level numbering")
    end if
  end do

  allocate(empty(0))
  call initialize_multilevel_hierarchy_1d( &
    base_cells, empty, empty, empty, 0.0_dp, 1.0_dp, base_hierarchy, ok)
  call assert_true(ok .and. base_hierarchy%level_count() == 1, &
    "base-only hierarchy")
  bad_lower = patch_lower
  bad_lower(2) = 1
  call initialize_multilevel_hierarchy_1d( &
    base_cells, bad_lower, patch_upper, refinement_ratios, &
    0.0_dp, 1.0_dp, invalid_hierarchy, ok)
  call assert_true(ok, "nested boundary-touching patch acceptance")
  call assert_true( &
    invalid_hierarchy%interfaces(2)%touches_left_boundary(), &
    "nested left physical boundary geometry")

  do cell = 1, base_cells
    x = (real(cell, dp) - 0.5_dp) / real(base_cells, dp)
    root(1, cell) = 2.0_dp + 0.25_dp * x
    root(2, cell) = 3.0_dp
  end do
  call prolong_multilevel_1d(root, hierarchy, fields, ok)
  call assert_true(ok, "recursive conservative prolongation")
  call assert_synchronized(fields, "prolongation restriction identity")
  uniform_integral = hierarchy%level_dx(0) * sum(root, dim=2)
  call composite_integral_multilevel_1d( &
    fields, hierarchy, integral_before, ok)
  call assert_true(ok, "multilevel composite integral")
  call assert_close(maxval(abs(integral_before - uniform_integral)), &
    0.0_dp, tolerance, "prolongation preserves root integral")

  coarse_time_step = 0.12_dp
  call multilevel_subcycle_counts_1d(hierarchy, subcycle_counts, ok)
  call assert_true(ok, "multilevel subcycle counts")
  call assert_true(all(subcycle_counts == [1, 2, 6, 12]), &
    "cumulative subcycle products")
  call multilevel_subcycle_time_steps_1d( &
    hierarchy, coarse_time_step, time_steps, ok)
  call assert_true(ok, "multilevel subcycle time steps")
  call assert_close(maxval(abs( &
    real(subcycle_counts, dp) * time_steps - coarse_time_step)), &
    0.0_dp, tolerance, "every level closes the coarse interval")

  averaged_fields = fields
  averaged_fields(level_count)%values(1, 1) = &
    averaged_fields(level_count)%values(1, 1) + 0.75_dp
  averaged_fields(level_count)%values(2, &
    expected_cells(level_count)) = &
    averaged_fields(level_count)%values(2, &
      expected_cells(level_count)) - 0.40_dp
  call composite_integral_multilevel_1d( &
    averaged_fields, hierarchy, integral_before, ok)
  call assert_true(ok, "perturbed composite integral")
  call average_down_multilevel_1d(averaged_fields, hierarchy, ok)
  call assert_true(ok, "deepest-to-root average down")
  call assert_synchronized(averaged_fields, "recursive average down")
  call composite_integral_multilevel_1d( &
    averaged_fields, hierarchy, integral_after, ok)
  call assert_true(ok, "averaged composite integral")
  call assert_close(maxval(abs(integral_after - integral_before)), &
    0.0_dp, tolerance, "average down preserves composite integral")

  flux_fields = fields
  do level = 1, level_count
    flux_fields(level)%values = 0.0_dp
  end do
  call initialize_multilevel_flux_registers_1d( &
    hierarchy, variable_count, registers, ok)
  call assert_true(ok, "multilevel flux-register allocation")
  do relation = 1, interface_count
    coarse_left = [real(relation, dp), -0.5_dp * real(relation, dp)]
    coarse_right = [2.0_dp * real(relation, dp), &
      0.25_dp * real(relation, dp)]
    fine_left = coarse_left + [0.5_dp, 0.75_dp]
    fine_right = coarse_right + [-0.25_dp, 0.50_dp]
    call accumulate_coarse_flux_1d( &
      registers(relation), coarse_left, coarse_right, &
      coarse_time_step, ok)
    call assert_true(ok, "multilevel coarse flux accumulation")
    call accumulate_fine_flux_1d( &
      registers(relation), fine_left, fine_right, coarse_time_step, ok)
    call assert_true(ok, "multilevel fine flux accumulation")

    flux_fields(relation)%values(:, &
      hierarchy%interfaces(relation)%fine_coarse_lower - 1) = &
      flux_fields(relation)%values(:, &
        hierarchy%interfaces(relation)%fine_coarse_lower - 1) - &
      coarse_time_step * coarse_left / &
        hierarchy%interfaces(relation)%coarse_dx
    flux_fields(relation)%values(:, &
      hierarchy%interfaces(relation)%fine_coarse_upper + 1) = &
      flux_fields(relation)%values(:, &
        hierarchy%interfaces(relation)%fine_coarse_upper + 1) + &
      coarse_time_step * coarse_right / &
        hierarchy%interfaces(relation)%coarse_dx
    fine_density = coarse_time_step * (fine_left - fine_right) / &
      (real(hierarchy%interfaces(relation)%fine%cell_count(), dp) * &
        hierarchy%interfaces(relation)%fine_dx)
    do level = relation + 1, level_count
      do cell = 1, hierarchy%level_cell_count(level - 1)
        flux_fields(level)%values(:, cell) = &
          flux_fields(level)%values(:, cell) + fine_density
      end do
    end do
  end do
  call composite_integral_multilevel_1d( &
    flux_fields, hierarchy, integral_before, ok)
  call assert_true(ok .and. maxval(abs(integral_before)) > 1.0e-3_dp, &
    "nontrivial multilevel interface mismatch")
  call synchronize_multilevel_1d( &
    flux_fields, hierarchy, registers, ok)
  call assert_true(ok, "deepest-to-root reflux synchronization")
  call composite_integral_multilevel_1d( &
    flux_fields, hierarchy, integral_after, ok)
  call assert_true(ok, "post-reflux composite integral")
  call assert_close(maxval(abs(integral_after)), 0.0_dp, tolerance, &
    "all-interface reflux restores conservation")
  call assert_synchronized(flux_fields, "post-reflux average down")
  do relation = 1, interface_count
    call assert_close(maxval(abs(registers(relation)%left)), 0.0_dp, &
      tolerance, "left multilevel register reset")
    call assert_close(maxval(abs(registers(relation)%right)), 0.0_dp, &
      tolerance, "right multilevel register reset")
  end do

  write(*, '(a)') "test_amr_multilevel_hierarchy_1d: PASS"

contains

  subroutine assert_synchronized(local_fields, label)
    type(amr_level_field_1d), intent(in) :: local_fields(:)
    character(len=*), intent(in) :: label

    logical :: local_ok
    integer :: local_relation

    do local_relation = 1, interface_count
      covered_cells = &
        hierarchy%interfaces(local_relation)%covered_coarse_cells()
      if (allocated(restricted)) deallocate(restricted)
      allocate(restricted(variable_count, covered_cells))
      call restrict_average_1d( &
        local_fields(local_relation + 1)%values, &
        hierarchy%interfaces(local_relation), restricted, local_ok)
      call assert_true(local_ok, trim(label) // " restriction")
      call assert_close(maxval(abs(restricted - &
        local_fields(local_relation)%values(:, &
          hierarchy%interfaces(local_relation)%fine_coarse_lower: &
          hierarchy%interfaces(local_relation)%fine_coarse_upper))), &
        0.0_dp, tolerance, label)
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

end program test_amr_multilevel_hierarchy_1d
