program test_multispecies_reconstruction
  use precision_mod, only: dp
  use state_indices_mod, only: ncons, nbase, nprim, qrho, qu, qv, qw, qp
  use state_conversion_mod, only: primitive_to_conserved
  use multispecies_state_mod, only: &
    multispecies_state_from_base, mass_fractions_from_state
  use reconstruction_multispecies_mod, only: reconstruct_multispecies_faces
  implicit none

  integer, parameter :: nx = 16
  integer, parameter :: nspecies = 2
  real(dp), parameter :: gamma = 1.4_dp
  real(dp), parameter :: dx = 1.0_dp / real(nx, dp)
  real(dp), parameter :: tolerance = 3.0e-12_dp
  real(dp) :: conserved(nbase + nspecies, 0:nx + 1)
  real(dp) :: left_faces(nbase + nspecies, 0:nx)
  real(dp) :: right_faces(nbase + nspecies, 0:nx)
  real(dp) :: primitive(nprim), base_state(ncons), y(nspecies), recovered(nspecies)
  real(dp) :: x_center, x_face, expected_y1
  logical :: ok
  integer :: i

  do i = 0, nx + 1
    x_center = (real(i, dp) - 0.5_dp) * dx
    primitive(qrho) = 1.0_dp
    primitive(qu) = 0.0_dp
    primitive(qv) = 0.0_dp
    primitive(qw) = 0.0_dp
    primitive(qp) = 1.0_dp
    call primitive_to_conserved(primitive, gamma, base_state, ok)
    call assert_true(ok, "base conversion")
    y(1) = 0.25_dp + 0.2_dp * x_center
    y(2) = 1.0_dp - y(1)
    call multispecies_state_from_base( &
      base_state, y, nspecies, gamma, conserved(:, i), ok)
    call assert_true(ok, "state assembly")
  end do

  call reconstruct_multispecies_faces( &
    conserved, nx, nspecies, gamma, "plm", "mc", "outflow", &
    left_faces, right_faces, ok)
  call assert_true(ok, "PLM reconstruction")

  do i = 2, nx - 2
    x_face = real(i, dp) * dx
    expected_y1 = 0.25_dp + 0.2_dp * x_face
    call mass_fractions_from_state(left_faces(:, i), nspecies, recovered, ok)
    call assert_true(ok, "left face extraction")
    call assert_close(recovered(1), expected_y1, tolerance, "left linear Y1")
    call assert_close(sum(recovered), 1.0_dp, tolerance, "left closure")

    call mass_fractions_from_state(right_faces(:, i), nspecies, recovered, ok)
    call assert_true(ok, "right face extraction")
    call assert_close(recovered(1), expected_y1, tolerance, "right linear Y1")
    call assert_close(sum(recovered), 1.0_dp, tolerance, "right closure")
  end do

  call reconstruct_multispecies_faces( &
    conserved, nx, nspecies, gamma, "pelec_plm", "mc", "outflow", &
    left_faces, right_faces, ok, dtdx=0.0_dp, plm_order=4, &
    use_flattening=.false.)
  call assert_true(ok, "PeleC PLM reconstruction")
  do i = 3, nx - 3
    x_face = real(i, dp) * dx
    expected_y1 = 0.25_dp + 0.2_dp * x_face
    call mass_fractions_from_state(left_faces(:, i), nspecies, recovered, ok)
    call assert_true(ok, "PeleC left extraction")
    call assert_close(recovered(1), expected_y1, tolerance, &
      "PeleC left linear Y1")
  end do

  write(*, '(a)') "test_multispecies_reconstruction: PASS"

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

end program test_multispecies_reconstruction
