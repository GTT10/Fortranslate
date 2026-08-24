program test_eb_geometry_2d
  use precision_mod, only: dp
  use eb_geometry_2d_mod, only: &
    eb_geometry_2d, eb_covered_cell, eb_cut_cell, eb_regular_cell, &
    build_eb_geometry_2d
  implicit none

  integer, parameter :: nx = 10
  integer, parameter :: ny = 12
  real(dp), parameter :: tolerance = 3.0e-13_dp
  type(eb_geometry_2d) :: geometry
  real(dp) :: level_set(0:nx, 0:ny)
  real(dp) :: area, coarse_error, fine_error
  integer :: i, j, coarse_cut_cells, fine_cut_cells
  logical :: ok

  level_set = 1.0_dp
  call build_eb_geometry_2d( &
    level_set, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, geometry, ok)
  call require(ok .and. geometry%is_valid(), "regular geometry validity")
  call require(all(geometry%cell_type == eb_regular_cell), &
    "regular cell classification")
  call require(maxval(abs(geometry%volume_fraction - 1.0_dp)) == 0.0_dp, &
    "regular cell volume fractions")
  call require(maxval(abs(geometry%x_face_fraction - 1.0_dp)) == 0.0_dp &
    .and. maxval(abs(geometry%y_face_fraction - 1.0_dp)) == 0.0_dp, &
    "regular face fractions")

  level_set = -1.0_dp
  call build_eb_geometry_2d( &
    level_set, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, geometry, ok)
  call require(ok .and. geometry%is_valid(), "covered geometry validity")
  call require(all(geometry%cell_type == eb_covered_cell), &
    "covered cell classification")
  call require(maxval(abs(geometry%volume_fraction)) == 0.0_dp, &
    "covered cell volume fractions")
  call require(maxval(abs(geometry%x_face_fraction)) == 0.0_dp .and. &
    maxval(abs(geometry%y_face_fraction)) == 0.0_dp, &
    "covered face fractions")

  do j = 0, ny
    do i = 0, nx
      level_set(i, j) = real(i, dp) / real(nx, dp) - 0.37_dp
    end do
  end do
  call build_eb_geometry_2d( &
    level_set, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, geometry, ok)
  call require(ok .and. geometry%is_valid(), "vertical plane validity")
  area = sum(geometry%volume_fraction) * geometry%dx * geometry%dy
  call assert_close(area, 0.63_dp, tolerance, "vertical plane fluid area")
  call require(count(geometry%cell_type == eb_cut_cell) == ny, &
    "vertical plane cut-cell count")
  call assert_close( &
    maxval(abs(geometry%volume_fraction(4, :) - 0.30_dp)), &
    0.0_dp, tolerance, "vertical plane cut fraction")
  call assert_close( &
    maxval(abs(geometry%y_face_fraction(4, :) - 0.30_dp)), &
    0.0_dp, tolerance, "vertical plane open face fraction")
  call require(all(geometry%x_face_fraction(3, :) == 0.0_dp) .and. &
    all(geometry%x_face_fraction(4, :) == 1.0_dp), &
    "vertical plane closed and open faces")

  do j = 0, ny
    do i = 0, nx
      level_set(i, j) = real(i, dp) / real(nx, dp) + &
        real(j, dp) / real(ny, dp) - 0.8_dp
    end do
  end do
  call build_eb_geometry_2d( &
    level_set, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, geometry, ok)
  call require(ok .and. geometry%is_valid(), "diagonal plane validity")
  area = sum(geometry%volume_fraction) * geometry%dx * geometry%dy
  call assert_close(area, 0.68_dp, tolerance, "diagonal plane fluid area")
  call require(count(geometry%cell_type == eb_cut_cell) > 0, &
    "diagonal plane cut cells")

  call circle_area_error(20, coarse_error, coarse_cut_cells)
  call circle_area_error(40, fine_error, fine_cut_cells)
  call require(coarse_error > 0.0_dp .and. fine_error > 0.0_dp, &
    "curved geometry has nontrivial discretization error")
  call require(fine_error < 0.45_dp * coarse_error, &
    "circle area second-order refinement")
  call require(coarse_cut_cells > 0 .and. fine_cut_cells > coarse_cut_cells, &
    "circle cut-cell refinement")

  call build_eb_geometry_2d( &
    level_set, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, geometry, ok)
  call require(.not. ok, "invalid physical bounds rejection")

  write(*, '(a)') "test_eb_geometry_2d: PASS"

contains

  subroutine circle_area_error(n, error, cut_cells)
    integer, intent(in) :: n
    real(dp), intent(out) :: error
    integer, intent(out) :: cut_cells

    real(dp), parameter :: radius = 0.27_dp
    type(eb_geometry_2d) :: circle_geometry
    real(dp) :: circle_level_set(0:n, 0:n)
    real(dp) :: x, y, measured_area, exact_area
    logical :: local_ok
    integer :: local_i, local_j

    do local_j = 0, n
      y = real(local_j, dp) / real(n, dp)
      do local_i = 0, n
        x = real(local_i, dp) / real(n, dp)
        circle_level_set(local_i, local_j) = radius - &
          sqrt((x - 0.5_dp)**2 + (y - 0.5_dp)**2)
      end do
    end do
    call build_eb_geometry_2d( &
      circle_level_set, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      circle_geometry, local_ok)
    call require(local_ok .and. circle_geometry%is_valid(), &
      "circle geometry validity")
    measured_area = sum(circle_geometry%volume_fraction) * &
      circle_geometry%dx * circle_geometry%dy
    exact_area = acos(-1.0_dp) * radius * radius
    error = abs(measured_area - exact_area)
    cut_cells = count(circle_geometry%cell_type == eb_cut_cell)
  end subroutine circle_area_error

  subroutine assert_close(actual, expected, local_tolerance, message)
    real(dp), intent(in) :: actual, expected, local_tolerance
    character(len=*), intent(in) :: message

    call require(abs(actual - expected) <= local_tolerance, message)
  end subroutine assert_close

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_eb_geometry_2d
