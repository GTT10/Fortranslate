module mesh_mod
  use precision_mod, only: dp
  implicit none
  private

  public :: uniform_cell_centers

contains

  pure subroutine uniform_cell_centers(nx, x_min, x_max, x, dx)
    integer, intent(in) :: nx
    real(dp), intent(in) :: x_min, x_max
    real(dp), intent(out) :: x(nx)
    real(dp), intent(out) :: dx
    integer :: i

    dx = (x_max - x_min) / real(nx, dp)
    do concurrent (i = 1:nx)
      x(i) = x_min + (real(i, dp) - 0.5_dp) * dx
    end do
  end subroutine uniform_cell_centers

end module mesh_mod
