program test_reactive_ppm
  use precision_mod, only: dp
  use reactive_1d_mod, only: &
    reactive_ppm_interface_value, reactive_ppm_monotone_edges, &
    reactive_ppm_reconstruct_five, reactive_ppm_integrate_profile, &
    reactive_ppm_flattening_coefficient, &
    reactive_ppm_contact_steepening_factor, &
    reactive_ppm_apply_contact_steepening
  implicit none

  real(dp), parameter :: tolerance = 3.0e-14_dp
  real(dp) :: left_edge, right_edge, value, flattening, eta
  real(dp) :: stencil(5), integral_right(3), integral_left(3)
  real(dp) :: pressure_wide(-3:3), velocity_wide(-3:3)
  real(dp) :: density_contact(-2:2), pressure_contact(-2:2)
  logical :: ok
  integer :: i

  value = reactive_ppm_interface_value(-1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp)
  call assert_close(value, 0.5_dp, 2.0e-15_dp, "linear interface")

  value = reactive_ppm_interface_value(0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp)
  call assert_true(value >= 0.0_dp .and. value <= 1.0_dp, &
    "PPM interface adjacent bounds")

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

  stencil = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
  call reactive_ppm_reconstruct_five(stencil, 1.0_dp, left_edge, right_edge)
  call assert_close(left_edge, 2.5_dp, tolerance, &
    "PeleC five-point linear left")
  call assert_close(right_edge, 3.5_dp, tolerance, &
    "PeleC five-point linear right")

  call reactive_ppm_reconstruct_five(stencil, 0.0_dp, left_edge, right_edge)
  call assert_close(left_edge, 3.0_dp, tolerance, "fully flattened left")
  call assert_close(right_edge, 3.0_dp, tolerance, "fully flattened right")

  call reactive_ppm_integrate_profile( &
    2.0_dp, 4.0_dp, 3.0_dp, 1.0_dp, 2.0_dp, 0.1_dp, &
    integral_right, integral_left, ok)
  call assert_true(ok, "PeleC parabolic profile integration")
  call assert_close(integral_right(1), 4.0_dp, tolerance, &
    "left acoustic right integral")
  call assert_close(integral_left(1), 2.1_dp, tolerance, &
    "left acoustic left integral")
  call assert_close(integral_right(2), 3.9_dp, tolerance, &
    "contact right integral")
  call assert_close(integral_left(2), 2.0_dp, tolerance, &
    "contact left integral")
  call assert_close(integral_right(3), 3.7_dp, tolerance, &
    "right acoustic right integral")
  call assert_close(integral_left(3), 2.0_dp, tolerance, &
    "right acoustic left integral")

  call reactive_ppm_integrate_profile( &
    2.0_dp, 4.0_dp, 3.0_dp, 2.0_dp, 2.0_dp, 0.3_dp, &
    integral_right, integral_left, ok)
  call assert_true(.not. ok, "super-CFL profile rejection")

  pressure_wide = 1.0_dp
  velocity_wide = 0.0_dp
  flattening = reactive_ppm_flattening_coefficient( &
    pressure_wide, velocity_wide)
  call assert_close(flattening, 1.0_dp, tolerance, &
    "smooth reactive PPM flattening")

  do i = -3, 3
    if (i <= 0) then
      pressure_wide(i) = 1.0_dp
      velocity_wide(i) = 1.0_dp
    else
      pressure_wide(i) = 10.0_dp
      velocity_wide(i) = -1.0_dp
    end if
  end do
  flattening = reactive_ppm_flattening_coefficient( &
    pressure_wide, velocity_wide)
  call assert_close(flattening, 0.0_dp, tolerance, &
    "compressive reactive PPM flattening")

  do i = -3, 3
    if (i <= 0) then
      velocity_wide(i) = -1.0_dp
    else
      velocity_wide(i) = 1.0_dp
    end if
  end do
  flattening = reactive_ppm_flattening_coefficient( &
    pressure_wide, velocity_wide)
  call assert_close(flattening, 1.0_dp, tolerance, &
    "expansion is not flattened")

  density_contact = [1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp]
  pressure_contact = 1.0_dp
  eta = reactive_ppm_contact_steepening_factor( &
    density_contact, pressure_contact, 1.4_dp)
  call assert_close(eta, 1.0_dp, tolerance, "contact steepening detector")

  left_edge = 1.25_dp
  right_edge = 1.75_dp
  call reactive_ppm_apply_contact_steepening( &
    density_contact, eta, left_edge, right_edge)
  call assert_close(left_edge, 1.0_dp, tolerance, "steepened contact left")
  call assert_close(right_edge, 2.0_dp, tolerance, "steepened contact right")

  pressure_contact = [1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp]
  eta = reactive_ppm_contact_steepening_factor( &
    density_contact, pressure_contact, 1.4_dp)
  call assert_close(eta, 0.0_dp, tolerance, &
    "pressure jump blocks contact steepening")

  write(*, '(a)') "test_reactive_ppm: PASS"

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

end program test_reactive_ppm
