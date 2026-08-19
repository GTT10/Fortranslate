program test_reactive_ppm
  use precision_mod, only: dp
  use reactive_1d_mod, only: &
    reactive_ppm_interface_value, reactive_ppm_monotone_edges
  implicit none

  real(dp) :: left_edge, right_edge, value

  value = reactive_ppm_interface_value(-1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp)
  call assert_close(value, 0.5_dp, 2.0e-15_dp, "linear interface")

  value = reactive_ppm_interface_value(0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp)
  if (value < 0.0_dp .or. value > 1.0_dp) then
    error stop "PPM interface escaped adjacent bounds"
  end if

  left_edge = 0.5_dp
  right_edge = 1.5_dp
  call reactive_ppm_monotone_edges(1.0_dp, left_edge, right_edge)
  call assert_close(left_edge, 0.5_dp, 2.0e-15_dp, "monotone left edge")
  call assert_close(right_edge, 1.5_dp, 2.0e-15_dp, "monotone right edge")

  left_edge = 0.0_dp
  right_edge = 0.1_dp
  call reactive_ppm_monotone_edges(1.0_dp, left_edge, right_edge)
  call assert_close(left_edge, 1.0_dp, 2.0e-15_dp, "extremum left edge")
  call assert_close(right_edge, 1.0_dp, 2.0e-15_dp, "extremum right edge")

  left_edge = 0.0_dp
  right_edge = 1.0_dp
  call reactive_ppm_monotone_edges(0.9_dp, left_edge, right_edge)
  call assert_close(left_edge, 0.7_dp, 2.0e-15_dp, "curvature-limited left")
  call assert_close(right_edge, 1.0_dp, 2.0e-15_dp, "curvature-limited right")

  write(*, '(a)') "test_reactive_ppm: PASS"

contains

  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label

    if (abs(actual - expected) > tolerance) then
      write(*, '(a,2(1x,es24.16))') trim(label), actual, expected
      error stop "Reactive PPM scalar mismatch"
    end if
  end subroutine assert_close

end program test_reactive_ppm
