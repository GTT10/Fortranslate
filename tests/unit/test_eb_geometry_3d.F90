program test_eb_geometry_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use precision_mod, only: dp
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d, eb_cut_cell_3d, &
    eb_regular_cell_3d, build_axis_plane_eb_geometry_3d
  implicit none

  integer, parameter :: nx = 10
  integer, parameter :: ny = 8
  integer, parameter :: nz = 6
  real(dp), parameter :: tolerance = 4.0e-13_dp
  type(eb_geometry_3d) :: geometry
  real(dp) :: fluid_volume, nan_value
  logical :: ok

  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", -0.1_dp, geometry, ok)
  call require(ok .and. geometry%is_valid(), "regular geometry validity")
  call require(all(geometry%cell_type == eb_regular_cell_3d), &
    "regular classification")
  call require(maxval(abs(geometry%volume_fraction - 1.0_dp)) == 0.0_dp, &
    "regular volume fractions")
  call require(maxval(abs(geometry%x_face_fraction - 1.0_dp)) == 0.0_dp &
      .and. maxval(abs(geometry%y_face_fraction - 1.0_dp)) == 0.0_dp &
      .and. maxval(abs(geometry%z_face_fraction - 1.0_dp)) == 0.0_dp, &
    "regular face fractions")
  call require(maxval(abs(geometry%boundary_area)) == 0.0_dp, &
    "regular boundary area")
  call assert_close(maximum_geometric_residual(geometry), 0.0_dp, &
    tolerance, "regular geometric identity")
  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", 0.0_dp, geometry, ok)
  call require(ok .and. all(geometry%cell_type == eb_regular_cell_3d), &
    "domain-lower plane regular collapse")

  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "z", 1.1_dp, geometry, ok)
  call require(ok .and. geometry%is_valid(), "covered geometry validity")
  call require(all(geometry%cell_type == eb_covered_cell_3d), &
    "covered classification")
  call require(maxval(abs(geometry%volume_fraction)) == 0.0_dp, &
    "covered volume fractions")
  call require(maxval(abs(geometry%x_face_fraction)) == 0.0_dp .and. &
      maxval(abs(geometry%y_face_fraction)) == 0.0_dp .and. &
      maxval(abs(geometry%z_face_fraction)) == 0.0_dp, &
    "covered face fractions")
  call assert_close(maximum_geometric_residual(geometry), 0.0_dp, &
    tolerance, "covered geometric identity")
  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "z", 1.0_dp, geometry, ok)
  call require(ok .and. all(geometry%cell_type == eb_covered_cell_3d), &
    "domain-upper plane covered collapse")

  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", 0.37_dp, geometry, ok)
  call require(ok .and. geometry%is_valid(), "x-plane geometry validity")
  fluid_volume = sum(geometry%volume_fraction) * &
    geometry%dx * geometry%dy * geometry%dz
  call assert_close(fluid_volume, 0.63_dp, tolerance, &
    "x-plane fluid volume")
  call require(count(geometry%cell_type == eb_cut_cell_3d) == ny * nz, &
    "x-plane cut-cell count")
  call assert_close(maxval(abs(geometry%volume_fraction(4, :, :) - &
    0.30_dp)), 0.0_dp, tolerance, "x-plane cut fraction")
  call assert_close(maxval(abs(geometry%cell_centroid_x(4, :, :) - &
    0.35_dp)), 0.0_dp, tolerance, "x-plane cell centroid")
  call require(maxval(abs(geometry%cell_centroid_y)) == 0.0_dp .and. &
    maxval(abs(geometry%cell_centroid_z)) == 0.0_dp, &
    "x-plane transverse cell centroids")
  call require(all(geometry%x_face_fraction(3, :, :) == 0.0_dp) .and. &
    all(geometry%x_face_fraction(4, :, :) == 1.0_dp), &
    "x-plane normal face fractions")
  call assert_close(maxval(abs(geometry%y_face_fraction(4, :, :) - &
    0.30_dp)), 0.0_dp, tolerance, "x-plane y-face fraction")
  call assert_close(maxval(abs(geometry%z_face_fraction(4, :, :) - &
    0.30_dp)), 0.0_dp, tolerance, "x-plane z-face fraction")
  call assert_close(maxval(abs(geometry%y_face_centroid_x(4, :, :) - &
    0.35_dp)), 0.0_dp, tolerance, "x-plane y-face centroid")
  call assert_close(maxval(abs(geometry%z_face_centroid_x(4, :, :) - &
    0.35_dp)), 0.0_dp, tolerance, "x-plane z-face centroid")
  call assert_close(sum(geometry%boundary_area), 1.0_dp, tolerance, &
    "x-plane boundary area")
  call assert_close(maxval(abs(geometry%boundary_centroid_x(4, :, :) - &
    0.37_dp)), 0.0_dp, tolerance, "x-plane boundary centroid")
  call assert_close(maxval(abs(geometry%boundary_normal_x(4, :, :) - &
    1.0_dp)), 0.0_dp, tolerance, "x-plane boundary normal")
  call require(maxval(abs(geometry%boundary_normal_y)) == 0.0_dp .and. &
    maxval(abs(geometry%boundary_normal_z)) == 0.0_dp, &
    "x-plane transverse boundary normals")
  call assert_close(maximum_geometric_residual(geometry), 0.0_dp, &
    tolerance, "x-plane geometric identity")

  geometry%boundary_normal_integral_x(4, 2, 3) = &
    0.9_dp * geometry%boundary_normal_integral_x(4, 2, 3)
  call require(.not. geometry%is_valid(), &
    "normal-integral corruption rejection")
  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", 0.37_dp, geometry, ok)
  call require(ok, "x-plane rebuild")
  geometry%y_face_fraction(4, 2, 3) = 0.4_dp
  call require(.not. geometry%is_valid(), "face-divergence corruption rejection")
  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", 0.37_dp, geometry, ok)
  call require(ok, "x-plane second rebuild")
  geometry%y_face_centroid_x(4, 2, 3) = 0.6_dp
  call require(.not. geometry%is_valid(), "face-centroid corruption rejection")

  call verify_rotated_plane("y", nx * nz)
  call verify_rotated_plane("z", nx * ny)

  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "q", 0.37_dp, geometry, ok)
  call require(.not. ok, "invalid axis rejection")
  call build_axis_plane_eb_geometry_3d( &
    0, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", 0.37_dp, geometry, ok)
  call require(.not. ok, "invalid extent rejection")
  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", 0.37_dp, geometry, ok)
  call require(.not. ok, "invalid bound rejection")
  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", 0.50_dp, geometry, ok)
  call require(.not. ok, "face-aligned interior plane rejection")
  nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", nan_value, geometry, ok)
  call require(.not. ok, "nonfinite plane rejection")

  write(*, '(a)') "test_eb_geometry_3d: PASS"

contains

  subroutine verify_rotated_plane(axis, expected_cut_cells)
    character(len=*), intent(in) :: axis
    integer, intent(in) :: expected_cut_cells

    real(dp) :: local_volume
    logical :: local_ok

    call build_axis_plane_eb_geometry_3d( &
      nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      axis, 0.37_dp, geometry, local_ok)
    call require(local_ok .and. geometry%is_valid(), &
      trim(axis)//"-plane geometry validity")
    local_volume = sum(geometry%volume_fraction) * &
      geometry%dx * geometry%dy * geometry%dz
    call assert_close(local_volume, 0.63_dp, tolerance, &
      trim(axis)//"-plane fluid volume")
    call require(count(geometry%cell_type == eb_cut_cell_3d) == &
      expected_cut_cells, trim(axis)//"-plane cut-cell count")
    call assert_close(sum(geometry%boundary_area), 1.0_dp, tolerance, &
      trim(axis)//"-plane boundary area")
    call assert_close(maximum_geometric_residual(geometry), 0.0_dp, &
      tolerance, trim(axis)//"-plane geometric identity")
    select case (axis)
    case ("y")
      call assert_close(maxval(abs(geometry%volume_fraction(:, 3, :) - &
        0.04_dp)), 0.0_dp, tolerance, "y-plane cut fraction")
      call assert_close(maxval(abs(geometry%cell_centroid_y(:, 3, :) - &
        0.48_dp)), 0.0_dp, tolerance, "y-plane cell centroid")
      call assert_close(maxval(abs(geometry%x_face_centroid_y(:, 3, :) - &
        0.48_dp)), 0.0_dp, tolerance, "y-plane x-face centroid")
      call assert_close(maxval(abs(geometry%z_face_centroid_y(:, 3, :) - &
        0.48_dp)), 0.0_dp, tolerance, "y-plane z-face centroid")
      call assert_close(maxval(abs(geometry%boundary_normal_y(:, 3, :) - &
        1.0_dp)), 0.0_dp, tolerance, "y-plane boundary normal")
    case ("z")
      call assert_close(maxval(abs(geometry%volume_fraction(:, :, 3) - &
        0.78_dp)), 0.0_dp, tolerance, "z-plane cut fraction")
      call assert_close(maxval(abs(geometry%cell_centroid_z(:, :, 3) - &
        0.11_dp)), 0.0_dp, tolerance, "z-plane cell centroid")
      call assert_close(maxval(abs(geometry%x_face_centroid_z(:, :, 3) - &
        0.11_dp)), 0.0_dp, tolerance, "z-plane x-face centroid")
      call assert_close(maxval(abs(geometry%y_face_centroid_z(:, :, 3) - &
        0.11_dp)), 0.0_dp, tolerance, "z-plane y-face centroid")
      call assert_close(maxval(abs(geometry%boundary_normal_z(:, :, 3) - &
        1.0_dp)), 0.0_dp, tolerance, "z-plane boundary normal")
    end select
  end subroutine verify_rotated_plane

  real(dp) function maximum_geometric_residual(local_geometry) result(error)
    type(eb_geometry_3d), intent(in) :: local_geometry

    real(dp) :: residual(3)
    integer :: i, j, k

    error = 0.0_dp
    do k = 1, local_geometry%nz
      do j = 1, local_geometry%ny
        do i = 1, local_geometry%nx
          residual = [ &
            local_geometry%dy * local_geometry%dz * ( &
              local_geometry%x_face_fraction(i, j, k) - &
              local_geometry%x_face_fraction(i - 1, j, k)) - &
                local_geometry%boundary_normal_integral_x(i, j, k), &
            local_geometry%dx * local_geometry%dz * ( &
              local_geometry%y_face_fraction(i, j, k) - &
              local_geometry%y_face_fraction(i, j - 1, k)) - &
                local_geometry%boundary_normal_integral_y(i, j, k), &
            local_geometry%dx * local_geometry%dy * ( &
              local_geometry%z_face_fraction(i, j, k) - &
              local_geometry%z_face_fraction(i, j, k - 1)) - &
                local_geometry%boundary_normal_integral_z(i, j, k)]
          error = max(error, maxval(abs(residual)))
        end do
      end do
    end do
  end function maximum_geometric_residual

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

end program test_eb_geometry_3d
