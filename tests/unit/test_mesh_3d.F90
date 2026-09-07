program test_mesh_3d
  use precision_mod, only: dp
  use mesh_3d_mod, only: uniform_cell_centers_3d
  implicit none

  integer, parameter :: nx = 4, ny = 3, nz = 2
  real(dp), parameter :: tolerance = 5.0e-15_dp
  real(dp) :: x(nx), y(ny), z(nz), dx, dy, dz

  call uniform_cell_centers_3d( &
    nx, ny, nz, -1.0_dp, 1.0_dp, 2.0_dp, 5.0_dp, -2.0_dp, 2.0_dp, &
    x, y, z, dx, dy, dz)

  call assert_close(dx, 0.5_dp, tolerance, "x spacing")
  call assert_close(dy, 1.0_dp, tolerance, "y spacing")
  call assert_close(dz, 2.0_dp, tolerance, "z spacing")
  call assert_close(x(1), -0.75_dp, tolerance, "first x center")
  call assert_close(x(nx), 0.75_dp, tolerance, "last x center")
  call assert_close(y(1), 2.5_dp, tolerance, "first y center")
  call assert_close(y(ny), 4.5_dp, tolerance, "last y center")
  call assert_close(z(1), -1.0_dp, tolerance, "first z center")
  call assert_close(z(nz), 1.0_dp, tolerance, "last z center")

  write(*, '(a)') "test_mesh_3d: PASS"

contains

  subroutine assert_close(actual, expected, tol, label)
    real(dp), intent(in) :: actual, expected, tol
    character(len=*), intent(in) :: label

    if (abs(actual - expected) > tol) then
      write(*, '(a,1x,a,2(1x,es24.16))') &
        "FAIL:", trim(label), actual, expected
      error stop 1
    end if
  end subroutine assert_close

end program test_mesh_3d
