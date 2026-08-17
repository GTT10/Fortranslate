module mesh_2d_mod
  use precision_mod, only: dp
  implicit none
  private

  public :: uniform_cell_centers_2d

contains

  pure subroutine uniform_cell_centers_2d( &
      nx, ny, x_min, x_max, y_min, y_max, x, y, dx, dy)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: x_min, x_max, y_min, y_max
    real(dp), intent(out) :: x(nx), y(ny)
    real(dp), intent(out) :: dx, dy
    integer :: i, j

    dx = (x_max - x_min) / real(nx, dp)
    dy = (y_max - y_min) / real(ny, dp)

    do concurrent (i = 1:nx)
      x(i) = x_min + (real(i, dp) - 0.5_dp) * dx
    end do
    do concurrent (j = 1:ny)
      y(j) = y_min + (real(j, dp) - 0.5_dp) * dy
    end do
  end subroutine uniform_cell_centers_2d

end module mesh_2d_mod
