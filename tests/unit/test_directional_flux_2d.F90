program test_directional_flux_2d
  use precision_mod, only: dp
  use state_indices_mod, only: &
    ncons, nprim, irho, imx, imy, imz, iet, qrho, qu, qv, qw, qp
  use state_conversion_mod, only: primitive_to_conserved
  use directional_flux_mod, only: &
    rotate_conserved_y_to_x, rotate_conserved_x_to_y, &
    euler_physical_flux_y, compute_riemann_flux_y
  implicit none

  real(dp), parameter :: gamma = 1.4_dp
  real(dp), parameter :: tolerance = 3.0e-13_dp
  real(dp) :: primitive(nprim), conserved(ncons), rotated(ncons), recovered(ncons)
  real(dp) :: physical_flux(ncons), numerical_flux(ncons)
  real(dp) :: expected_energy
  logical :: ok

  primitive(qrho) = 1.2_dp
  primitive(qu) = 0.3_dp
  primitive(qv) = -0.2_dp
  primitive(qw) = 0.1_dp
  primitive(qp) = 1.1_dp

  call primitive_to_conserved(primitive, gamma, conserved, ok)
  call assert_true(ok, "primitive conversion")

  call rotate_conserved_y_to_x(conserved, rotated)
  call rotate_conserved_x_to_y(rotated, recovered)
  call assert_close(maxval(abs(recovered - conserved)), 0.0_dp, tolerance, &
    "conserved rotation round trip")

  call euler_physical_flux_y(conserved, gamma, physical_flux, ok)
  call assert_true(ok, "physical y flux")

  expected_energy = conserved(iet)
  call assert_close(physical_flux(irho), 1.2_dp * (-0.2_dp), tolerance, &
    "mass flux")
  call assert_close(physical_flux(imx), 1.2_dp * 0.3_dp * (-0.2_dp), &
    tolerance, "x momentum flux")
  call assert_close(physical_flux(imy), 1.2_dp * (-0.2_dp)**2 + 1.1_dp, &
    tolerance, "y momentum flux")
  call assert_close(physical_flux(imz), 1.2_dp * 0.1_dp * (-0.2_dp), &
    tolerance, "z momentum flux")
  call assert_close(physical_flux(iet), (expected_energy + 1.1_dp) * (-0.2_dp), &
    tolerance, "energy flux")

  call compute_riemann_flux_y( &
    conserved, conserved, gamma, "rusanov", numerical_flux, ok)
  call assert_true(ok, "Rusanov y flux")
  call assert_close(maxval(abs(numerical_flux - physical_flux)), 0.0_dp, &
    tolerance, "Rusanov equal-state consistency")

  call compute_riemann_flux_y( &
    conserved, conserved, gamma, "pelec", numerical_flux, ok)
  call assert_true(ok, "PeleC-style y flux")
  call assert_close(maxval(abs(numerical_flux - physical_flux)), 0.0_dp, &
    tolerance, "PeleC equal-state consistency")

  write(*, '(a)') "test_directional_flux_2d: PASS"

contains

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      write(*, '(a,1x,a)') "FAIL:", trim(label)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected, tolerance_value, label)
    real(dp), intent(in) :: actual, expected, tolerance_value
    character(len=*), intent(in) :: label

    if (abs(actual - expected) > tolerance_value) then
      write(*, '(a,1x,a,2(1x,es24.16))') &
        "FAIL:", trim(label), actual, expected
      error stop 1
    end if
  end subroutine assert_close

end program test_directional_flux_2d
