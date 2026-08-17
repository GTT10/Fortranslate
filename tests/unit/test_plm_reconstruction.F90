program test_plm_reconstruction
  use precision_mod, only: dp
  use state_indices_mod, only: ncons, nprim, qrho, qu, qv, qw, qp
  use state_conversion_mod, only: primitive_to_conserved, conserved_to_primitive
  use reconstruction_plm_mod, only: reconstruct_plm_faces
  implicit none

  integer, parameter :: nx = 12
  real(dp), parameter :: gamma = 1.4_dp
  real(dp), parameter :: dx = 1.0_dp / real(nx, dp)
  real(dp), parameter :: tolerance = 2.0e-13_dp
  real(dp) :: conserved(ncons, 0:nx + 1)
  real(dp) :: left_faces(ncons, 0:nx), right_faces(ncons, 0:nx)
  real(dp) :: primitive(nprim), face_primitive(nprim), x_center, x_face
  logical :: ok
  integer :: i

  do i = 0, nx + 1
    x_center = (real(i, dp) - 0.5_dp) * dx
    primitive(qrho) = 1.0_dp + 0.1_dp * x_center
    primitive(qu) = 0.2_dp + 0.05_dp * x_center
    primitive(qv) = -0.03_dp * x_center
    primitive(qw) = 0.02_dp * x_center
    primitive(qp) = 1.0_dp + 0.2_dp * x_center
    call primitive_to_conserved(primitive, gamma, conserved(:, i), ok)
    call assert_true(ok, "linear-state conversion")
  end do

  call reconstruct_plm_faces( &
    conserved, nx, gamma, "mc", "outflow", left_faces, right_faces, ok)
  call assert_true(ok, "PLM reconstruction")

  do i = 2, nx - 2
    x_face = real(i, dp) * dx

    call conserved_to_primitive(left_faces(:, i), gamma, face_primitive, ok)
    call assert_true(ok, "left face conversion")
    call check_linear_state(face_primitive, x_face, "left face")

    call conserved_to_primitive(right_faces(:, i), gamma, face_primitive, ok)
    call assert_true(ok, "right face conversion")
    call check_linear_state(face_primitive, x_face, "right face")
  end do

  write(*, '(a)') "test_plm_reconstruction: PASS"

contains

  subroutine check_linear_state(state, coordinate, label)
    real(dp), intent(in) :: state(nprim), coordinate
    character(len=*), intent(in) :: label

    call assert_close(state(qrho), 1.0_dp + 0.1_dp * coordinate, tolerance, &
      trim(label) // " density")
    call assert_close(state(qu), 0.2_dp + 0.05_dp * coordinate, tolerance, &
      trim(label) // " velocity x")
    call assert_close(state(qv), -0.03_dp * coordinate, tolerance, &
      trim(label) // " velocity y")
    call assert_close(state(qw), 0.02_dp * coordinate, tolerance, &
      trim(label) // " velocity z")
    call assert_close(state(qp), 1.0_dp + 0.2_dp * coordinate, tolerance, &
      trim(label) // " pressure")
  end subroutine check_linear_state

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

end program test_plm_reconstruction
