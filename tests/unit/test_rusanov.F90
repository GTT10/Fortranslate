program test_rusanov
  use precision_mod, only: dp
  use state_indices_mod, only: &
    ncons, nprim, qrho, qu, qv, qw, qp, irho, imx, imy, imz, iet
  use state_conversion_mod, only: primitive_to_conserved
  use riemann_rusanov_mod, only: euler_physical_flux_x, rusanov_flux_x
  implicit none

  real(dp), parameter :: gamma = 1.4_dp
  real(dp), parameter :: tolerance = 2.0e-13_dp
  real(dp) :: primitive(nprim), conserved(ncons), physical_flux(ncons), numerical_flux(ncons)
  logical :: ok

  primitive = 0.0_dp
  primitive(qrho) = 1.0_dp
  primitive(qu) = 0.5_dp
  primitive(qv) = 0.2_dp
  primitive(qw) = -0.1_dp
  primitive(qp) = 1.0_dp

  call primitive_to_conserved(primitive, gamma, conserved, ok)
  call assert_true(ok, "state conversion")
  call euler_physical_flux_x(conserved, gamma, physical_flux, ok)
  call assert_true(ok, "physical flux")
  call rusanov_flux_x(conserved, conserved, gamma, numerical_flux, ok)
  call assert_true(ok, "Rusanov flux")

  call assert_close(maxval(abs(numerical_flux - physical_flux)), 0.0_dp, tolerance, &
    "equal-state consistency")
  call assert_close(physical_flux(irho), 0.5_dp, tolerance, "mass flux")
  call assert_close(physical_flux(imx), 1.25_dp, tolerance, "x momentum flux")
  call assert_close(physical_flux(imy), 0.1_dp, tolerance, "y momentum flux")
  call assert_close(physical_flux(imz), -0.05_dp, tolerance, "z momentum flux")
  call assert_close(physical_flux(iet), 1.825_dp, tolerance, "energy flux")

  write(*, '(a)') "test_rusanov: PASS"

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

end program test_rusanov
