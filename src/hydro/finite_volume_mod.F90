module finite_volume_mod
  use precision_mod, only: dp
  use state_indices_mod, only: ncons
  use riemann_rusanov_mod, only: rusanov_flux_x
  implicit none
  private

  public :: compute_euler_rhs

contains

  subroutine compute_euler_rhs(conserved, nx, dx, gamma, rhs, ok)
    integer, intent(in) :: nx
    real(dp), intent(in) :: conserved(ncons, 0:nx + 1)
    real(dp), intent(in) :: dx, gamma
    real(dp), intent(out) :: rhs(ncons, nx)
    logical, intent(out) :: ok
    real(dp), allocatable :: face_flux(:, :)
    logical :: face_ok
    integer :: i

    allocate(face_flux(ncons, 0:nx))
    rhs = 0.0_dp
    ok = .true.

    do i = 0, nx
      call rusanov_flux_x( &
        conserved(:, i), conserved(:, i + 1), gamma, face_flux(:, i), face_ok)
      if (.not. face_ok) then
        ok = .false.
        return
      end if
    end do

    do concurrent (i = 1:nx)
      rhs(:, i) = -(face_flux(:, i) - face_flux(:, i - 1)) / dx
    end do
  end subroutine compute_euler_rhs

end module finite_volume_mod
