program test_pelec_fourth_order_flattening
  use precision_mod, only: dp
  use state_indices_mod, only: nprim, qrho, qu, qv, qw, qp
  use reconstruction_pelec_plm_mod, only: &
    pelec_limited_slope, pelec_flattening_coefficient
  implicit none

  integer, parameter :: nx = 11
  integer, parameter :: center_cell = 6
  real(dp), parameter :: tolerance = 2.0e-14_dp
  real(dp) :: primitive(nprim, -2:nx + 3)
  real(dp) :: slope, flat
  logical :: ok
  integer :: i

  call pelec_limited_slope( &
    0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
    1.0_dp, 4, slope, ok)
  call assert_true(ok, "fourth-order linear slope")
  call assert_close(slope, 1.0_dp, tolerance, "fourth-order linear exactness")

  call pelec_limited_slope( &
    0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
    1.0_dp, 2, slope, ok)
  call assert_true(ok, "second-order PeleC slope")
  call assert_close( &
    slope, 4.0_dp / 3.0_dp, tolerance, "second-order slope formula")

  call pelec_limited_slope( &
    0.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 0.0_dp, &
    1.0_dp, 4, slope, ok)
  call assert_true(ok, "extremum slope")
  call assert_close(slope, 0.0_dp, tolerance, "extremum limiting")

  call pelec_limited_slope( &
    0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
    0.25_dp, 4, slope, ok)
  call assert_true(ok, "flattened slope")
  call assert_close(slope, 0.25_dp, tolerance, "flattening multiplication")

  call pelec_limited_slope( &
    0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
    1.0_dp, 3, slope, ok)
  call assert_true(.not. ok, "invalid order rejection")

  primitive = 0.0_dp
  primitive(qrho, :) = 1.0_dp
  primitive(qp, :) = 1.0_dp
  flat = pelec_flattening_coefficient( &
    primitive, nx, center_cell, "periodic")
  call assert_close(flat, 1.0_dp, tolerance, "smooth-region flattening")

  do i = -2, nx + 3
    primitive(qrho, i) = 1.0_dp
    primitive(qv, i) = 0.0_dp
    primitive(qw, i) = 0.0_dp
    if (i <= center_cell) then
      primitive(qp, i) = 1.0_dp
      primitive(qu, i) = 1.0_dp
    else
      primitive(qp, i) = 10.0_dp
      primitive(qu, i) = -1.0_dp
    end if
  end do
  flat = pelec_flattening_coefficient( &
    primitive, nx, center_cell, "periodic")
  call assert_close(flat, 0.0_dp, tolerance, "compressive shock flattening")

  do i = -2, nx + 3
    if (i <= center_cell) then
      primitive(qu, i) = -1.0_dp
    else
      primitive(qu, i) = 1.0_dp
    end if
  end do
  flat = pelec_flattening_coefficient( &
    primitive, nx, center_cell, "periodic")
  call assert_close(flat, 1.0_dp, tolerance, "expansion is not flattened")

  write(*, '(a)') "test_pelec_fourth_order_flattening: PASS"

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

end program test_pelec_fourth_order_flattening
