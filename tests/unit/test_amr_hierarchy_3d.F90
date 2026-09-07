program test_amr_hierarchy_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: &
    amr_patch_3d, initialize_amr_patch_3d, prolong_pcm_3d, &
    restrict_average_3d, average_down_3d, composite_integrals_amr_3d
  implicit none

  integer, parameter :: nvar = 4, nx = 6, ny = 5, nz = 4
  type(amr_patch_3d) :: patch
  real(dp) :: coarse(nvar, nx, ny, nz), averaged(nvar, nx, ny, nz)
  real(dp), allocatable :: fine(:, :, :, :), restricted(:, :, :, :)
  real(dp) :: composite(nvar), reference(nvar), error
  logical :: ok
  integer :: component, i, j, k

  call initialize_amr_patch_3d( &
    nx, ny, nz, 2, 4, 2, 4, 2, 3, 2, patch, ok)
  call require(ok, "patch initialization")
  call require(patch%is_strictly_interior(), "strictly interior patch")
  call require( &
    patch%fine_nx() == 6 .and. patch%fine_ny() == 6 .and. &
    patch%fine_nz() == 4, "fine extents")

  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        do component = 1, nvar
          coarse(component, i, j, k) = real( &
            1000 * component + 100 * k + 10 * j + i, dp)
        end do
      end do
    end do
  end do
  allocate(fine(nvar, patch%fine_nx(), patch%fine_ny(), patch%fine_nz()))
  allocate(restricted(nvar, 3, 3, 2))
  call prolong_pcm_3d(coarse, patch, fine, ok)
  call require(ok, "PCM prolongation")
  call restrict_average_3d(fine, patch, restricted, ok)
  call require(ok, "restriction")
  error = maxval(abs(restricted - coarse(:, 2:4, 2:4, 2:3)))
  call require(error == 0.0_dp, "PCM restriction identity")

  call composite_integrals_amr_3d( &
    coarse, fine, patch, 0.2_dp, 0.3_dp, 0.4_dp, composite, ok)
  call require(ok, "composite integration")
  do component = 1, nvar
    reference(component) = 0.2_dp * 0.3_dp * 0.4_dp * &
      sum(coarse(component, :, :, :))
  end do
  error = maxval(abs(composite - reference) / max(1.0_dp, abs(reference)))
  call require(error <= 5.0e-15_dp, "composite PCM integral")
  call composite_integrals_amr_3d( &
    coarse, fine, patch, ieee_value(0.0_dp, ieee_quiet_nan), 0.3_dp, &
    0.4_dp, composite, ok)
  call require(.not. ok .and. all(composite == 0.0_dp), &
    "non-finite composite spacing rejection")

  fine = fine + 0.125_dp
  averaged = coarse
  call average_down_3d(averaged, fine, patch, ok)
  call require(ok, "average down")
  error = maxval(abs(averaged(:, 2:4, 2:4, 2:3) - &
    (coarse(:, 2:4, 2:4, 2:3) + 0.125_dp)))
  call require(error <= 5.0e-13_dp, "covered coarse synchronization")
  call require(maxval(abs(averaged(:, 1, :, :) - coarse(:, 1, :, :))) == &
    0.0_dp, "uncovered coarse preservation")

  call initialize_amr_patch_3d( &
    nx, ny, nz, 0, 4, 2, 4, 2, 3, 2, patch, ok)
  call require(.not. ok, "invalid patch rejection")
  write(*, '(a)') "test_amr_hierarchy_3d: PASS"

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
      write(*, '(a)') "FAIL: " // trim(message)
      error stop 1
    end if
  end subroutine require

end program test_amr_hierarchy_3d
