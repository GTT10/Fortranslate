program test_weno_reconstruction
  use precision_mod, only: dp
  use reconstruction_weno_mod, only: &
    weno_reconstruct_5js, weno_reconstruct_5z
  implicit none

  real(dp) :: stencil(5), left_edge, right_edge

  stencil = 3.0_dp
  call weno_reconstruct_5js(stencil, left_edge, right_edge)
  call assert_close(left_edge, 3.0_dp, 2.0e-14_dp, "WENO5-JS constant left")
  call assert_close(right_edge, 3.0_dp, 2.0e-14_dp, &
    "WENO5-JS constant right")
  call weno_reconstruct_5z(stencil, left_edge, right_edge)
  call assert_close(left_edge, 3.0_dp, 2.0e-14_dp, "WENO5-Z constant left")
  call assert_close(right_edge, 3.0_dp, 2.0e-14_dp, &
    "WENO5-Z constant right")

  stencil = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
  call weno_reconstruct_5js(stencil, left_edge, right_edge)
  call assert_close(left_edge, 2.5_dp, 3.0e-14_dp, "WENO5-JS linear left")
  call assert_close(right_edge, 3.5_dp, 3.0e-14_dp, &
    "WENO5-JS linear right")
  call weno_reconstruct_5z(stencil, left_edge, right_edge)
  call assert_close(left_edge, 2.5_dp, 3.0e-14_dp, "WENO5-Z linear left")
  call assert_close(right_edge, 3.5_dp, 3.0e-14_dp, &
    "WENO5-Z linear right")

  ! These nonsymmetric values are direct parity points for PeleC WENO.H.
  stencil = [2.0_dp, -1.0_dp, 3.0_dp, 0.0_dp, 4.0_dp]
  call weno_reconstruct_5js(stencil, left_edge, right_edge)
  call assert_close(left_edge, 2.1374741359913365_dp, 3.0e-14_dp, &
    "WENO5-JS PeleC left")
  call assert_close(right_edge, 2.5161323342134514_dp, 3.0e-14_dp, &
    "WENO5-JS PeleC right")
  call weno_reconstruct_5z(stencil, left_edge, right_edge)
  call assert_close(left_edge, 1.9427611938090070_dp, 3.0e-14_dp, &
    "WENO5-Z PeleC left")
  call assert_close(right_edge, 2.4372350525428820_dp, 3.0e-14_dp, &
    "WENO5-Z PeleC right")

  write(*, '(a)') "test_weno_reconstruction: PASS"

contains

  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label

    if (abs(actual - expected) > tolerance) then
      write(*, '(a,1x,a,2(1x,es24.16))') &
        "FAIL:", trim(label), actual, expected
      error stop 1
    end if
  end subroutine assert_close

end program test_weno_reconstruction
