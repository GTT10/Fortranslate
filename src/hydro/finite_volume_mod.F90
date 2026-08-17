module finite_volume_mod
  use precision_mod, only: dp
  use state_indices_mod, only: ncons
  use riemann_rusanov_mod, only: rusanov_flux_x
  use reconstruction_plm_mod, only: reconstruct_plm_faces
  implicit none
  private

  public :: compute_euler_rhs

contains

  subroutine compute_euler_rhs( &
      conserved, nx, dx, gamma, rhs, ok, reconstruction, limiter, boundary_condition)
    integer, intent(in) :: nx
    real(dp), intent(in) :: conserved(ncons, 0:nx + 1)
    real(dp), intent(in) :: dx, gamma
    real(dp), intent(out) :: rhs(ncons, nx)
    logical, intent(out) :: ok
    character(len=*), intent(in), optional :: reconstruction, limiter, boundary_condition

    real(dp), allocatable :: face_flux(:, :), left_faces(:, :), right_faces(:, :)
    character(len=32) :: reconstruction_name, limiter_name, boundary_name
    logical :: face_ok, reconstruction_ok
    integer :: i

    reconstruction_name = "pcm"
    limiter_name = "mc"
    boundary_name = "outflow"
    if (present(reconstruction)) reconstruction_name = trim(reconstruction)
    if (present(limiter)) limiter_name = trim(limiter)
    if (present(boundary_condition)) boundary_name = trim(boundary_condition)

    allocate(face_flux(ncons, 0:nx))
    rhs = 0.0_dp
    ok = .true.

    select case (trim(reconstruction_name))
    case ("pcm")
      do i = 0, nx
        call rusanov_flux_x( &
          conserved(:, i), conserved(:, i + 1), gamma, face_flux(:, i), face_ok)
        if (.not. face_ok) then
          ok = .false.
          return
        end if
      end do

    case ("plm")
      allocate(left_faces(ncons, 0:nx))
      allocate(right_faces(ncons, 0:nx))
      call reconstruct_plm_faces( &
        conserved, nx, gamma, limiter_name, boundary_name, &
        left_faces, right_faces, reconstruction_ok)
      if (.not. reconstruction_ok) then
        ok = .false.
        return
      end if

      do i = 0, nx
        call rusanov_flux_x( &
          left_faces(:, i), right_faces(:, i), gamma, face_flux(:, i), face_ok)
        if (.not. face_ok) then
          ok = .false.
          return
        end if
      end do

    case default
      ok = .false.
      return
    end select

    do concurrent (i = 1:nx)
      rhs(:, i) = -(face_flux(:, i) - face_flux(:, i - 1)) / dx
    end do
  end subroutine compute_euler_rhs

end module finite_volume_mod
