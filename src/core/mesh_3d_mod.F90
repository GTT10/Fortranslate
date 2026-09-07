module mesh_3d_mod
  use precision_mod, only: dp
  use mesh_mod, only: uniform_cell_centers
  implicit none
  private

  public :: uniform_cell_centers_3d

contains

  pure subroutine uniform_cell_centers_3d( &
      nx, ny, nz, x_min, x_max, y_min, y_max, z_min, z_max, &
      x, y, z, dx, dy, dz)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: x_min, x_max, y_min, y_max, z_min, z_max
    real(dp), intent(out) :: x(nx), y(ny), z(nz)
    real(dp), intent(out) :: dx, dy, dz

    call uniform_cell_centers(nx, x_min, x_max, x, dx)
    call uniform_cell_centers(ny, y_min, y_max, y, dy)
    call uniform_cell_centers(nz, z_min, z_max, z, dz)
  end subroutine uniform_cell_centers_3d

end module mesh_3d_mod
