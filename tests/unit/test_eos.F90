program test_eos
  use precision_mod, only: dp
  use state_indices_mod, only: &
    ncons, nprim, qrho, qu, qv, qw, qp, irho, imx, imy, imz, iet
  use state_conversion_mod, only: primitive_to_conserved, conserved_to_primitive
  use eos_ideal_mod, only: ideal_gas_sound_speed
  implicit none

  real(dp), parameter :: gamma = 1.4_dp
  real(dp), parameter :: tolerance = 1.0e-13_dp
  real(dp) :: primitive(nprim), recovered(nprim), conserved(ncons)
  logical :: ok

  primitive = 0.0_dp
  primitive(qrho) = 1.0_dp
  primitive(qu) = 2.0_dp
  primitive(qv) = -0.5_dp
  primitive(qw) = 0.25_dp
  primitive(qp) = 1.0_dp

  call primitive_to_conserved(primitive, gamma, conserved, ok)
  call assert_true(ok, "primitive_to_conserved returned false")
  call assert_close(conserved(irho), 1.0_dp, tolerance, "density")
  call assert_close(conserved(imx), 2.0_dp, tolerance, "x momentum")
  call assert_close(conserved(imy), -0.5_dp, tolerance, "y momentum")
  call assert_close(conserved(imz), 0.25_dp, tolerance, "z momentum")
  call assert_close(conserved(iet), 4.65625_dp, tolerance, "total energy")

  call conserved_to_primitive(conserved, gamma, recovered, ok)
  call assert_true(ok, "conserved_to_primitive returned false")
  call assert_close(maxval(abs(recovered - primitive)), 0.0_dp, tolerance, "round trip")
  call assert_close( &
    ideal_gas_sound_speed(1.0_dp, 1.0_dp, gamma), sqrt(gamma), tolerance, &
    "sound speed")

  write(*, '(a)') "test_eos: PASS"

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

end program test_eos
