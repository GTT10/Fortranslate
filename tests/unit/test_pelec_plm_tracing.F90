program test_pelec_plm_tracing
  use precision_mod, only: dp
  use state_indices_mod, only: nprim, qrho, qu, qv, qw, qp
  use reconstruction_pelec_plm_mod, only: &
    n_characteristic_waves, wave_minus, wave_contact, wave_shear_y, &
    wave_shear_z, primitive_slope_to_characteristics, &
    characteristics_to_primitive_slope, trace_primitive_characteristics
  implicit none

  real(dp), parameter :: gamma = 1.4_dp
  real(dp), parameter :: tolerance = 5.0e-13_dp
  real(dp) :: center(nprim), slope(nprim), recovered(nprim)
  real(dp) :: characteristic(n_characteristic_waves)
  real(dp) :: recovered_characteristic(n_characteristic_waves)
  real(dp) :: left_state(nprim), right_state(nprim)
  logical :: ok

  center(qrho) = 1.3_dp
  center(qu) = 0.4_dp
  center(qv) = -0.2_dp
  center(qw) = 0.1_dp
  center(qp) = 1.7_dp

  slope(qrho) = 0.07_dp
  slope(qu) = -0.03_dp
  slope(qv) = 0.02_dp
  slope(qw) = -0.01_dp
  slope(qp) = 0.11_dp

  call primitive_slope_to_characteristics( &
    center, slope, gamma, characteristic, ok)
  call assert_true(ok, "primitive-to-characteristic projection")
  call characteristics_to_primitive_slope( &
    center, characteristic, gamma, recovered, ok)
  call assert_true(ok, "characteristic-to-primitive projection")
  call assert_close(maxval(abs(recovered - slope)), 0.0_dp, tolerance, &
    "projection round trip")

  characteristic = 0.0_dp
  characteristic(wave_minus) = 0.025_dp
  call characteristics_to_primitive_slope( &
    center, characteristic, gamma, slope, ok)
  call assert_true(ok, "left acoustic synthesis")
  call primitive_slope_to_characteristics( &
    center, slope, gamma, recovered_characteristic, ok)
  call assert_true(ok, "left acoustic projection")
  call assert_close(recovered_characteristic(wave_minus), 0.025_dp, &
    tolerance, "left acoustic amplitude")
  call assert_close(maxval(abs(recovered_characteristic(2:5))), 0.0_dp, &
    tolerance, "left acoustic isolation")

  slope(qrho) = 0.07_dp
  slope(qu) = -0.03_dp
  slope(qv) = 0.02_dp
  slope(qw) = -0.01_dp
  slope(qp) = 0.11_dp
  call trace_primitive_characteristics( &
    center, slope, gamma, 0.0_dp, left_state, right_state, ok)
  call assert_true(ok, "zero-Courant tracing")
  call assert_close(maxval(abs(left_state - (center - 0.5_dp * slope))), &
    0.0_dp, tolerance, "zero-Courant left state")
  call assert_close(maxval(abs(right_state - (center + 0.5_dp * slope))), &
    0.0_dp, tolerance, "zero-Courant right state")

  characteristic = 0.0_dp
  characteristic(wave_contact) = 0.08_dp
  characteristic(wave_shear_y) = -0.03_dp
  characteristic(wave_shear_z) = 0.02_dp
  call characteristics_to_primitive_slope( &
    center, characteristic, gamma, slope, ok)
  call assert_true(ok, "contact/shear synthesis")
  call trace_primitive_characteristics( &
    center, slope, gamma, 0.1_dp, left_state, right_state, ok)
  call assert_true(ok, "contact/shear tracing")
  call assert_close(left_state(qu), center(qu), tolerance, &
    "left contact velocity")
  call assert_close(right_state(qu), center(qu), tolerance, &
    "right contact velocity")
  call assert_close(left_state(qp), center(qp), tolerance, &
    "left contact pressure")
  call assert_close(right_state(qp), center(qp), tolerance, &
    "right contact pressure")

  call trace_primitive_characteristics( &
    center, slope, gamma, -0.1_dp, left_state, right_state, ok)
  call assert_true(.not. ok, "negative dtdx rejection")

  write(*, '(a)') "test_pelec_plm_tracing: PASS"

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

end program test_pelec_plm_tracing
