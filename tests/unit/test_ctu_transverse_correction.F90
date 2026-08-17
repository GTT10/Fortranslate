program test_ctu_transverse_correction
  use precision_mod, only: dp
  use state_indices_mod, only: &
    ncons, nprim, irho, imx, iet, qrho, qu, qv, qw, qp
  use state_conversion_mod, only: primitive_to_conserved, state_is_physical
  use ctu_2d_mod, only: apply_transverse_flux_correction
  implicit none

  real(dp), parameter :: gamma = 1.4_dp
  real(dp), parameter :: tolerance = 5.0e-13_dp
  real(dp) :: primitive(nprim), base_state(ncons)
  real(dp) :: flux_high(ncons), flux_low(ncons), corrected(ncons), expected(ncons)
  real(dp) :: theta, scale
  logical :: ok

  primitive(qrho) = 1.0_dp
  primitive(qu) = 0.2_dp
  primitive(qv) = -0.1_dp
  primitive(qw) = 0.0_dp
  primitive(qp) = 1.0_dp
  call primitive_to_conserved(primitive, gamma, base_state, ok)
  call assert_true(ok, "base state conversion")

  flux_high = 0.0_dp
  flux_low = 0.0_dp
  call apply_transverse_flux_correction( &
    base_state, flux_high, flux_low, 0.25_dp, gamma, corrected, theta, ok)
  call assert_true(ok, "zero correction")
  call assert_close(theta, 1.0_dp, tolerance, "zero correction theta")
  call assert_close(maxval(abs(corrected - base_state)), 0.0_dp, tolerance, &
    "zero correction identity")

  flux_high = 0.0_dp
  flux_low = 0.0_dp
  flux_high(imx) = 0.3_dp
  flux_high(iet) = 0.2_dp
  scale = 0.1_dp
  expected = base_state - scale * (flux_high - flux_low)
  call apply_transverse_flux_correction( &
    base_state, flux_high, flux_low, scale, gamma, corrected, theta, ok)
  call assert_true(ok, "unlimited physical correction")
  call assert_close(theta, 1.0_dp, tolerance, "unlimited theta")
  call assert_close(maxval(abs(corrected - expected)), 0.0_dp, tolerance, &
    "unlimited correction value")

  flux_high = 0.0_dp
  flux_low = 0.0_dp
  flux_high(irho) = 20.0_dp
  flux_high(iet) = 40.0_dp
  scale = 0.1_dp
  call apply_transverse_flux_correction( &
    base_state, flux_high, flux_low, scale, gamma, corrected, theta, ok)
  call assert_true(ok, "limited correction")
  call assert_true(theta > 0.0_dp .and. theta < 1.0_dp, &
    "limited theta range")
  call assert_true(state_is_physical(corrected, gamma), &
    "limited correction physical state")
  expected = base_state - theta * scale * (flux_high - flux_low)
  call assert_close(maxval(abs(corrected - expected)), 0.0_dp, tolerance, &
    "limited correction line search")

  call apply_transverse_flux_correction( &
    base_state, flux_high, flux_low, -0.1_dp, gamma, corrected, theta, ok)
  call assert_true(.not. ok, "negative scale rejection")

  write(*, '(a)') "test_ctu_transverse_correction: PASS"

contains

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      write(*, '(a,1x,a)') "FAIL:", trim(label)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected_value, tolerance_value, label)
    real(dp), intent(in) :: actual, expected_value, tolerance_value
    character(len=*), intent(in) :: label

    if (abs(actual - expected_value) > tolerance_value) then
      write(*, '(a,1x,a,2(1x,es24.16))') &
        "FAIL:", trim(label), actual, expected_value
      error stop 1
    end if
  end subroutine assert_close

end program test_ctu_transverse_correction
