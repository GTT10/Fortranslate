program test_slope_limiter
  use precision_mod, only: dp
  use slope_limiter_mod, only: minmod2, minmod3, limited_slope
  implicit none

  real(dp), parameter :: tolerance = 1.0e-14_dp
  real(dp) :: slope
  logical :: ok

  call assert_close(minmod2(2.0_dp, 1.0_dp), 1.0_dp, tolerance, "minmod positive")
  call assert_close(minmod2(-2.0_dp, -1.0_dp), -1.0_dp, tolerance, "minmod negative")
  call assert_close(minmod2(-1.0_dp, 1.0_dp), 0.0_dp, tolerance, "minmod extremum")
  call assert_close(minmod3(3.0_dp, 2.0_dp, 1.0_dp), 1.0_dp, tolerance, "minmod3")

  call limited_slope(1.0_dp, 3.0_dp, "minmod", slope, ok)
  call assert_true(ok, "minmod selection")
  call assert_close(slope, 1.0_dp, tolerance, "minmod slope")

  call limited_slope(1.0_dp, 3.0_dp, "mc", slope, ok)
  call assert_true(ok, "MC selection")
  call assert_close(slope, 2.0_dp, tolerance, "MC slope")

  call limited_slope(1.0_dp, -1.0_dp, "mc", slope, ok)
  call assert_true(ok, "MC extremum selection")
  call assert_close(slope, 0.0_dp, tolerance, "MC extremum")

  call limited_slope(1.0_dp, 1.0_dp, "unknown", slope, ok)
  call assert_true(.not. ok, "unknown limiter rejection")

  write(*, '(a)') "test_slope_limiter: PASS"

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
      write(*, '(a,1x,a,2(1x,es24.16))') "FAIL:", trim(label), actual, expected
      error stop 1
    end if
  end subroutine assert_close

end program test_slope_limiter
