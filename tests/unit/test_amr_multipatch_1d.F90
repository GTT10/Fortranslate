program test_amr_multipatch_1d
  use precision_mod, only: dp
  use amr_hierarchy_1d_mod, only: &
    amr_level_field_1d, amr_flux_register_1d, restrict_average_1d, &
    accumulate_coarse_flux_1d, accumulate_fine_flux_1d
  use amr_multipatch_1d_mod, only: &
    amr_patch_set_1d, initialize_patch_set_1d, &
    prolong_patch_set_1d, average_down_patch_set_1d, &
    initialize_patch_flux_registers_1d, synchronize_patch_set_1d, &
    composite_integral_patch_set_1d
  implicit none

  integer, parameter :: variable_count = 2
  integer, parameter :: coarse_cells = 14
  integer, parameter :: patch_count = 2
  integer, parameter :: ratio = 2
  real(dp), parameter :: tolerance = 5.0e-12_dp
  integer :: patch_lower(patch_count), patch_upper(patch_count)
  integer :: adjacent_lower(patch_count), adjacent_upper(patch_count)
  integer, allocatable :: empty(:)
  type(amr_patch_set_1d) :: patch_set, invalid_set, empty_set
  type(amr_level_field_1d), allocatable :: fine_fields(:)
  type(amr_level_field_1d), allocatable :: flux_fields(:)
  type(amr_level_field_1d), allocatable :: empty_fields(:)
  type(amr_flux_register_1d), allocatable :: registers(:)
  real(dp) :: coarse(variable_count, coarse_cells)
  real(dp) :: averaged(variable_count, coarse_cells)
  real(dp) :: flux_coarse(variable_count, coarse_cells)
  real(dp) :: restricted(variable_count, 3)
  real(dp) :: root_integral(variable_count)
  real(dp) :: composite_before(variable_count)
  real(dp) :: composite_after(variable_count)
  real(dp) :: coarse_left(variable_count), coarse_right(variable_count)
  real(dp) :: fine_left(variable_count), fine_right(variable_count)
  real(dp) :: fine_density(variable_count), dt, x
  logical :: ok
  integer :: cell, child_patch

  patch_lower = [1, 9]
  patch_upper = [3, 11]
  call initialize_patch_set_1d( &
    coarse_cells, patch_lower, patch_upper, ratio, 0.0_dp, 1.0_dp, &
    patch_set, ok)
  call assert_true(ok .and. patch_set%is_valid(), &
    "two-patch hierarchy initialization")
  call assert_true(patch_set%patch_count() == patch_count, &
    "patch count")
  call assert_true(patch_set%covered_coarse_cell_count() == 6, &
    "covered parent-cell count")
  call assert_true(patch_set%fine_cell_count() == 12, &
    "aggregate fine-cell count")
  call assert_true(patch_set%patches(1)%touches_left_boundary(), &
    "first patch physical-boundary contact")
  call assert_true(.not. patch_set%patches(2)%touches_left_boundary(), &
    "second patch interior left interface")
  call assert_true(patch_set%parent_cell_is_covered(2) .and. &
    patch_set%parent_cell_is_covered(10), "coverage lookup")
  call assert_true(.not. patch_set%parent_cell_is_covered(6), &
    "uncovered lookup")

  adjacent_lower = [2, 5]
  adjacent_upper = [4, 7]
  call initialize_patch_set_1d( &
    coarse_cells, adjacent_lower, adjacent_upper, ratio, 0.0_dp, 1.0_dp, &
    invalid_set, ok)
  call assert_true(.not. ok, "adjacent patches require coalescing")

  allocate(empty(0))
  call initialize_patch_set_1d( &
    coarse_cells, empty, empty, ratio, 0.0_dp, 1.0_dp, empty_set, ok)
  call assert_true(ok .and. empty_set%patch_count() == 0, &
    "empty patch set")

  do cell = 1, coarse_cells
    x = (real(cell, dp) - 0.5_dp) / real(coarse_cells, dp)
    coarse(1, cell) = 2.0_dp + 0.25_dp * x
    coarse(2, cell) = 3.0_dp
  end do
  root_integral = patch_set%coarse_dx * sum(coarse, dim=2)
  call prolong_patch_set_1d(coarse, patch_set, fine_fields, ok)
  call assert_true(ok, "multipatch conservative prolongation")
  call composite_integral_patch_set_1d( &
    coarse, fine_fields, patch_set, composite_before, ok)
  call assert_true(ok, "multipatch composite integral")
  call assert_close(maxval(abs(composite_before - root_integral)), &
    0.0_dp, tolerance, "multipatch prolongation conservation")

  averaged = -1.0_dp
  call average_down_patch_set_1d(fine_fields, patch_set, averaged, ok)
  call assert_true(ok, "multipatch average down")
  do child_patch = 1, patch_count
    call restrict_average_1d( &
      fine_fields(child_patch)%values, patch_set%patches(child_patch), &
      restricted, ok)
    call assert_true(ok, "patch restriction")
    call assert_close(maxval(abs(restricted - averaged(:, &
      patch_set%patches(child_patch)%fine_coarse_lower: &
      patch_set%patches(child_patch)%fine_coarse_upper))), &
      0.0_dp, tolerance, "covered cells synchronized")
  end do
  call assert_close(maxval(abs(averaged(:, 4:8) + 1.0_dp)), &
    0.0_dp, tolerance, "inter-patch coarse cells unchanged")
  call assert_close(maxval(abs(averaged(:, 12:14) + 1.0_dp)), &
    0.0_dp, tolerance, "upper coarse cells unchanged")

  call prolong_patch_set_1d( &
    0.0_dp * coarse, patch_set, flux_fields, ok)
  call assert_true(ok, "zero multipatch fields")
  flux_coarse = 0.0_dp
  dt = 0.08_dp
  call initialize_patch_flux_registers_1d( &
    patch_set, variable_count, registers, ok)
  call assert_true(ok, "per-patch flux registers")

  do child_patch = 1, patch_count
    if (child_patch == 1) then
      coarse_left = 0.0_dp
      fine_left = 0.0_dp
      coarse_right = [1.5_dp, -0.5_dp]
      fine_right = [0.5_dp, 0.75_dp]
    else
      coarse_left = [0.25_dp, 1.0_dp]
      fine_left = [1.0_dp, -0.25_dp]
      coarse_right = [2.0_dp, 0.5_dp]
      fine_right = [1.25_dp, 1.5_dp]
    end if
    call accumulate_coarse_flux_1d( &
      registers(child_patch), coarse_left, coarse_right, dt, ok)
    call assert_true(ok, "patch coarse-flux accumulation")
    call accumulate_fine_flux_1d( &
      registers(child_patch), fine_left, fine_right, dt, ok)
    call assert_true(ok, "patch fine-flux accumulation")
    if (.not. patch_set%patches(child_patch)%touches_left_boundary()) then
      cell = patch_set%patches(child_patch)%fine_coarse_lower - 1
      flux_coarse(:, cell) = flux_coarse(:, cell) - &
        dt * coarse_left / patch_set%coarse_dx
    end if
    if (.not. patch_set%patches(child_patch)%touches_right_boundary()) then
      cell = patch_set%patches(child_patch)%fine_coarse_upper + 1
      flux_coarse(:, cell) = flux_coarse(:, cell) + &
        dt * coarse_right / patch_set%coarse_dx
    end if
    fine_density = dt * (fine_left - fine_right) / &
      (real(size(flux_fields(child_patch)%values, 2), dp) * &
        patch_set%fine_dx)
    do cell = 1, size(flux_fields(child_patch)%values, 2)
      flux_fields(child_patch)%values(:, cell) = fine_density
    end do
  end do

  call composite_integral_patch_set_1d( &
    flux_coarse, flux_fields, patch_set, composite_before, ok)
  call assert_true(ok .and. maxval(abs(composite_before)) > 1.0e-3_dp, &
    "nontrivial multipatch flux mismatch")
  call synchronize_patch_set_1d( &
    flux_coarse, flux_fields, patch_set, registers, ok)
  call assert_true(ok, "multipatch reflux and average down")
  call composite_integral_patch_set_1d( &
    flux_coarse, flux_fields, patch_set, composite_after, ok)
  call assert_true(ok, "post-reflux multipatch integral")
  call assert_close(maxval(abs(composite_after)), 0.0_dp, tolerance, &
    "multipatch reflux restores conservation")
  do child_patch = 1, patch_count
    call restrict_average_1d( &
      flux_fields(child_patch)%values, patch_set%patches(child_patch), &
      restricted, ok)
    call assert_true(ok, "post-reflux patch restriction")
    call assert_close(maxval(abs(restricted - flux_coarse(:, &
      patch_set%patches(child_patch)%fine_coarse_lower: &
      patch_set%patches(child_patch)%fine_coarse_upper))), &
      0.0_dp, tolerance, "post-reflux covered cells synchronized")
    call assert_close(maxval(abs(registers(child_patch)%left)), &
      0.0_dp, tolerance, "patch left register reset")
    call assert_close(maxval(abs(registers(child_patch)%right)), &
      0.0_dp, tolerance, "patch right register reset")
  end do

  call prolong_patch_set_1d(coarse, empty_set, empty_fields, ok)
  call assert_true(ok .and. size(empty_fields) == 0, &
    "empty patch prolongation")
  call composite_integral_patch_set_1d( &
    coarse, empty_fields, empty_set, composite_after, ok)
  call assert_true(ok, "empty patch composite integral")
  call assert_close(maxval(abs(composite_after - root_integral)), &
    0.0_dp, tolerance, "empty patch retains root integral")

  write(*, '(a)') "test_amr_multipatch_1d: PASS"

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

end program test_amr_multipatch_1d
