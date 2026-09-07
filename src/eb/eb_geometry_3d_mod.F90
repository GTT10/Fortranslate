module eb_geometry_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  implicit none
  private

  integer, parameter, public :: eb_covered_cell_3d = 0
  integer, parameter, public :: eb_cut_cell_3d = 1
  integer, parameter, public :: eb_regular_cell_3d = 2
  real(dp), parameter :: classification_tolerance = &
    128.0_dp * epsilon(1.0_dp)

  type, public :: eb_geometry_3d
    integer :: nx = 0
    integer :: ny = 0
    integer :: nz = 0
    real(dp) :: x_lower = 0.0_dp
    real(dp) :: x_upper = 0.0_dp
    real(dp) :: y_lower = 0.0_dp
    real(dp) :: y_upper = 0.0_dp
    real(dp) :: z_lower = 0.0_dp
    real(dp) :: z_upper = 0.0_dp
    real(dp) :: dx = 0.0_dp
    real(dp) :: dy = 0.0_dp
    real(dp) :: dz = 0.0_dp
    real(dp), allocatable :: volume_fraction(:, :, :)
    ! Cell-centroid offsets relative to the Cartesian cell center,
    ! normalized by the corresponding cell width (AMReX convention).
    real(dp), allocatable :: cell_centroid_x(:, :, :)
    real(dp), allocatable :: cell_centroid_y(:, :, :)
    real(dp), allocatable :: cell_centroid_z(:, :, :)
    real(dp), allocatable :: x_face_fraction(:, :, :)
    real(dp), allocatable :: y_face_fraction(:, :, :)
    real(dp), allocatable :: z_face_fraction(:, :, :)
    ! Tangential face-centroid offsets use the corresponding cell width.
    real(dp), allocatable :: x_face_centroid_y(:, :, :)
    real(dp), allocatable :: x_face_centroid_z(:, :, :)
    real(dp), allocatable :: y_face_centroid_x(:, :, :)
    real(dp), allocatable :: y_face_centroid_z(:, :, :)
    real(dp), allocatable :: z_face_centroid_x(:, :, :)
    real(dp), allocatable :: z_face_centroid_y(:, :, :)
    real(dp), allocatable :: boundary_area(:, :, :)
    real(dp), allocatable :: boundary_centroid_x(:, :, :)
    real(dp), allocatable :: boundary_centroid_y(:, :, :)
    real(dp), allocatable :: boundary_centroid_z(:, :, :)
    ! Embedded-boundary normals point from solid to fluid.
    real(dp), allocatable :: boundary_normal_x(:, :, :)
    real(dp), allocatable :: boundary_normal_y(:, :, :)
    real(dp), allocatable :: boundary_normal_z(:, :, :)
    real(dp), allocatable :: boundary_normal_integral_x(:, :, :)
    real(dp), allocatable :: boundary_normal_integral_y(:, :, :)
    real(dp), allocatable :: boundary_normal_integral_z(:, :, :)
    integer, allocatable :: cell_type(:, :, :)
  contains
    procedure :: is_valid => eb_geometry_3d_is_valid
  end type eb_geometry_3d

  public :: build_axis_plane_eb_geometry_3d

contains

  subroutine build_axis_plane_eb_geometry_3d( &
      nx, ny, nz, x_lower, x_upper, y_lower, y_upper, z_lower, z_upper, &
      axis, plane_position, geometry, ok)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: x_lower, x_upper, y_lower, y_upper
    real(dp), intent(in) :: z_lower, z_upper, plane_position
    character(len=*), intent(in) :: axis
    type(eb_geometry_3d), intent(out) :: geometry
    logical, intent(out) :: ok

    real(dp) :: axis_lower, axis_upper, axis_spacing
    real(dp) :: fraction, centroid_offset, coordinate, alignment_tolerance
    integer :: axis_index, i, j, k

    geometry = eb_geometry_3d()
    ok = .false.
    if (nx < 1 .or. ny < 1 .or. nz < 1 .or. &
        .not. all(ieee_is_finite([ &
          x_lower, x_upper, y_lower, y_upper, z_lower, z_upper, &
          plane_position])) .or. &
        x_upper <= x_lower .or. y_upper <= y_lower .or. &
        z_upper <= z_lower) return

    select case (trim(adjustl(axis)))
    case ("x", "X")
      axis_index = 1
      axis_lower = x_lower
      axis_upper = x_upper
      axis_spacing = (x_upper - x_lower) / real(nx, dp)
    case ("y", "Y")
      axis_index = 2
      axis_lower = y_lower
      axis_upper = y_upper
      axis_spacing = (y_upper - y_lower) / real(ny, dp)
    case ("z", "Z")
      axis_index = 3
      axis_lower = z_lower
      axis_upper = z_upper
      axis_spacing = (z_upper - z_lower) / real(nz, dp)
    case default
      return
    end select

    ! A face-aligned interior EB has no positive-volume cut cell to own its
    ! wall metric.  Reject that degenerate representation explicitly.
    if (plane_position > axis_lower .and. plane_position < axis_upper) then
      coordinate = (plane_position - axis_lower) / axis_spacing
      alignment_tolerance = 512.0_dp * epsilon(1.0_dp) * &
        max(1.0_dp, abs(coordinate))
      if (abs(coordinate - anint(coordinate)) <= alignment_tolerance) return
    end if

    geometry%nx = nx
    geometry%ny = ny
    geometry%nz = nz
    geometry%x_lower = x_lower
    geometry%x_upper = x_upper
    geometry%y_lower = y_lower
    geometry%y_upper = y_upper
    geometry%z_lower = z_lower
    geometry%z_upper = z_upper
    geometry%dx = (x_upper - x_lower) / real(nx, dp)
    geometry%dy = (y_upper - y_lower) / real(ny, dp)
    geometry%dz = (z_upper - z_lower) / real(nz, dp)
    call allocate_geometry_arrays(geometry)
    call clear_geometry_arrays(geometry)

    if (plane_position <= axis_lower) then
      geometry%volume_fraction = 1.0_dp
      geometry%x_face_fraction = 1.0_dp
      geometry%y_face_fraction = 1.0_dp
      geometry%z_face_fraction = 1.0_dp
      geometry%cell_type = eb_regular_cell_3d
      ok = geometry%is_valid()
      return
    else if (plane_position >= axis_upper) then
      geometry%cell_type = eb_covered_cell_3d
      ok = geometry%is_valid()
      return
    end if

    select case (axis_index)
    case (1)
      do i = 1, nx
        call positive_interval_metrics( &
          x_lower + real(i - 1, dp) * geometry%dx, &
          x_lower + real(i, dp) * geometry%dx, plane_position, &
          fraction, centroid_offset)
        geometry%volume_fraction(i, :, :) = fraction
        geometry%cell_centroid_x(i, :, :) = centroid_offset
        geometry%y_face_fraction(i, :, :) = fraction
        geometry%z_face_fraction(i, :, :) = fraction
        geometry%y_face_centroid_x(i, :, :) = centroid_offset
        geometry%z_face_centroid_x(i, :, :) = centroid_offset
        call classify_x_slab(geometry, i, fraction, plane_position)
      end do
      do i = 0, nx
        coordinate = x_lower + real(i, dp) * geometry%dx
        if (coordinate > plane_position) &
          geometry%x_face_fraction(i, :, :) = 1.0_dp
      end do
    case (2)
      do j = 1, ny
        call positive_interval_metrics( &
          y_lower + real(j - 1, dp) * geometry%dy, &
          y_lower + real(j, dp) * geometry%dy, plane_position, &
          fraction, centroid_offset)
        geometry%volume_fraction(:, j, :) = fraction
        geometry%cell_centroid_y(:, j, :) = centroid_offset
        geometry%x_face_fraction(:, j, :) = fraction
        geometry%z_face_fraction(:, j, :) = fraction
        geometry%x_face_centroid_y(:, j, :) = centroid_offset
        geometry%z_face_centroid_y(:, j, :) = centroid_offset
        call classify_y_slab(geometry, j, fraction, plane_position)
      end do
      do j = 0, ny
        coordinate = y_lower + real(j, dp) * geometry%dy
        if (coordinate > plane_position) &
          geometry%y_face_fraction(:, j, :) = 1.0_dp
      end do
    case (3)
      do k = 1, nz
        call positive_interval_metrics( &
          z_lower + real(k - 1, dp) * geometry%dz, &
          z_lower + real(k, dp) * geometry%dz, plane_position, &
          fraction, centroid_offset)
        geometry%volume_fraction(:, :, k) = fraction
        geometry%cell_centroid_z(:, :, k) = centroid_offset
        geometry%x_face_fraction(:, :, k) = fraction
        geometry%y_face_fraction(:, :, k) = fraction
        geometry%x_face_centroid_z(:, :, k) = centroid_offset
        geometry%y_face_centroid_z(:, :, k) = centroid_offset
        call classify_z_slab(geometry, k, fraction, plane_position)
      end do
      do k = 0, nz
        coordinate = z_lower + real(k, dp) * geometry%dz
        if (coordinate > plane_position) &
          geometry%z_face_fraction(:, :, k) = 1.0_dp
      end do
    end select

    ok = geometry%is_valid()
  end subroutine build_axis_plane_eb_geometry_3d

  subroutine allocate_geometry_arrays(geometry)
    type(eb_geometry_3d), intent(inout) :: geometry

    integer :: nx, ny, nz

    nx = geometry%nx
    ny = geometry%ny
    nz = geometry%nz
    allocate(geometry%volume_fraction(nx, ny, nz))
    allocate(geometry%cell_centroid_x(nx, ny, nz))
    allocate(geometry%cell_centroid_y(nx, ny, nz))
    allocate(geometry%cell_centroid_z(nx, ny, nz))
    allocate(geometry%x_face_fraction(0:nx, ny, nz))
    allocate(geometry%y_face_fraction(nx, 0:ny, nz))
    allocate(geometry%z_face_fraction(nx, ny, 0:nz))
    allocate(geometry%x_face_centroid_y(0:nx, ny, nz))
    allocate(geometry%x_face_centroid_z(0:nx, ny, nz))
    allocate(geometry%y_face_centroid_x(nx, 0:ny, nz))
    allocate(geometry%y_face_centroid_z(nx, 0:ny, nz))
    allocate(geometry%z_face_centroid_x(nx, ny, 0:nz))
    allocate(geometry%z_face_centroid_y(nx, ny, 0:nz))
    allocate(geometry%boundary_area(nx, ny, nz))
    allocate(geometry%boundary_centroid_x(nx, ny, nz))
    allocate(geometry%boundary_centroid_y(nx, ny, nz))
    allocate(geometry%boundary_centroid_z(nx, ny, nz))
    allocate(geometry%boundary_normal_x(nx, ny, nz))
    allocate(geometry%boundary_normal_y(nx, ny, nz))
    allocate(geometry%boundary_normal_z(nx, ny, nz))
    allocate(geometry%boundary_normal_integral_x(nx, ny, nz))
    allocate(geometry%boundary_normal_integral_y(nx, ny, nz))
    allocate(geometry%boundary_normal_integral_z(nx, ny, nz))
    allocate(geometry%cell_type(nx, ny, nz))
  end subroutine allocate_geometry_arrays

  subroutine clear_geometry_arrays(geometry)
    type(eb_geometry_3d), intent(inout) :: geometry

    geometry%volume_fraction = 0.0_dp
    geometry%cell_centroid_x = 0.0_dp
    geometry%cell_centroid_y = 0.0_dp
    geometry%cell_centroid_z = 0.0_dp
    geometry%x_face_fraction = 0.0_dp
    geometry%y_face_fraction = 0.0_dp
    geometry%z_face_fraction = 0.0_dp
    geometry%x_face_centroid_y = 0.0_dp
    geometry%x_face_centroid_z = 0.0_dp
    geometry%y_face_centroid_x = 0.0_dp
    geometry%y_face_centroid_z = 0.0_dp
    geometry%z_face_centroid_x = 0.0_dp
    geometry%z_face_centroid_y = 0.0_dp
    geometry%boundary_area = 0.0_dp
    geometry%boundary_centroid_x = 0.0_dp
    geometry%boundary_centroid_y = 0.0_dp
    geometry%boundary_centroid_z = 0.0_dp
    geometry%boundary_normal_x = 0.0_dp
    geometry%boundary_normal_y = 0.0_dp
    geometry%boundary_normal_z = 0.0_dp
    geometry%boundary_normal_integral_x = 0.0_dp
    geometry%boundary_normal_integral_y = 0.0_dp
    geometry%boundary_normal_integral_z = 0.0_dp
    geometry%cell_type = eb_covered_cell_3d
  end subroutine clear_geometry_arrays

  pure subroutine positive_interval_metrics( &
      lower, upper, plane_position, fraction, centroid_offset)
    real(dp), intent(in) :: lower, upper, plane_position
    real(dp), intent(out) :: fraction, centroid_offset

    if (plane_position <= lower) then
      fraction = 1.0_dp
      centroid_offset = 0.0_dp
    else if (plane_position >= upper) then
      fraction = 0.0_dp
      centroid_offset = 0.0_dp
    else
      fraction = (upper - plane_position) / (upper - lower)
      fraction = min(1.0_dp, max(0.0_dp, fraction))
      centroid_offset = 0.5_dp * (1.0_dp - fraction)
    end if
  end subroutine positive_interval_metrics

  subroutine classify_x_slab(geometry, i, fraction, plane_position)
    type(eb_geometry_3d), intent(inout) :: geometry
    integer, intent(in) :: i
    real(dp), intent(in) :: fraction, plane_position

    integer :: j, k

    if (fraction <= classification_tolerance) then
      geometry%cell_type(i, :, :) = eb_covered_cell_3d
    else if (fraction >= 1.0_dp - classification_tolerance) then
      geometry%cell_type(i, :, :) = eb_regular_cell_3d
    else
      geometry%cell_type(i, :, :) = eb_cut_cell_3d
      do k = 1, geometry%nz
        do j = 1, geometry%ny
          geometry%boundary_area(i, j, k) = geometry%dy * geometry%dz
          geometry%boundary_centroid_x(i, j, k) = plane_position
          geometry%boundary_centroid_y(i, j, k) = geometry%y_lower + &
            (real(j, dp) - 0.5_dp) * geometry%dy
          geometry%boundary_centroid_z(i, j, k) = geometry%z_lower + &
            (real(k, dp) - 0.5_dp) * geometry%dz
          geometry%boundary_normal_x(i, j, k) = 1.0_dp
          geometry%boundary_normal_integral_x(i, j, k) = &
            geometry%boundary_area(i, j, k)
        end do
      end do
    end if
  end subroutine classify_x_slab

  subroutine classify_y_slab(geometry, j, fraction, plane_position)
    type(eb_geometry_3d), intent(inout) :: geometry
    integer, intent(in) :: j
    real(dp), intent(in) :: fraction, plane_position

    integer :: i, k

    if (fraction <= classification_tolerance) then
      geometry%cell_type(:, j, :) = eb_covered_cell_3d
    else if (fraction >= 1.0_dp - classification_tolerance) then
      geometry%cell_type(:, j, :) = eb_regular_cell_3d
    else
      geometry%cell_type(:, j, :) = eb_cut_cell_3d
      do k = 1, geometry%nz
        do i = 1, geometry%nx
          geometry%boundary_area(i, j, k) = geometry%dx * geometry%dz
          geometry%boundary_centroid_x(i, j, k) = geometry%x_lower + &
            (real(i, dp) - 0.5_dp) * geometry%dx
          geometry%boundary_centroid_y(i, j, k) = plane_position
          geometry%boundary_centroid_z(i, j, k) = geometry%z_lower + &
            (real(k, dp) - 0.5_dp) * geometry%dz
          geometry%boundary_normal_y(i, j, k) = 1.0_dp
          geometry%boundary_normal_integral_y(i, j, k) = &
            geometry%boundary_area(i, j, k)
        end do
      end do
    end if
  end subroutine classify_y_slab

  subroutine classify_z_slab(geometry, k, fraction, plane_position)
    type(eb_geometry_3d), intent(inout) :: geometry
    integer, intent(in) :: k
    real(dp), intent(in) :: fraction, plane_position

    integer :: i, j

    if (fraction <= classification_tolerance) then
      geometry%cell_type(:, :, k) = eb_covered_cell_3d
    else if (fraction >= 1.0_dp - classification_tolerance) then
      geometry%cell_type(:, :, k) = eb_regular_cell_3d
    else
      geometry%cell_type(:, :, k) = eb_cut_cell_3d
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          geometry%boundary_area(i, j, k) = geometry%dx * geometry%dy
          geometry%boundary_centroid_x(i, j, k) = geometry%x_lower + &
            (real(i, dp) - 0.5_dp) * geometry%dx
          geometry%boundary_centroid_y(i, j, k) = geometry%y_lower + &
            (real(j, dp) - 0.5_dp) * geometry%dy
          geometry%boundary_centroid_z(i, j, k) = plane_position
          geometry%boundary_normal_z(i, j, k) = 1.0_dp
          geometry%boundary_normal_integral_z(i, j, k) = &
            geometry%boundary_area(i, j, k)
        end do
      end do
    end if
  end subroutine classify_z_slab

  pure logical function eb_geometry_3d_is_valid(self) result(valid)
    class(eb_geometry_3d), intent(in) :: self

    real(dp), parameter :: tolerance = 512.0_dp * epsilon(1.0_dp)

    valid = self%nx >= 1 .and. self%ny >= 1 .and. self%nz >= 1 .and. &
      all(ieee_is_finite([ &
        self%x_lower, self%x_upper, self%y_lower, self%y_upper, &
        self%z_lower, self%z_upper, self%dx, self%dy, self%dz])) .and. &
      self%x_upper > self%x_lower .and. &
      self%y_upper > self%y_lower .and. &
      self%z_upper > self%z_lower .and. &
      self%dx > 0.0_dp .and. self%dy > 0.0_dp .and. self%dz > 0.0_dp
    if (.not. valid) return
    valid = geometry_arrays_are_valid(self)
    if (.not. valid) return
    valid = spacing_is_valid(self, tolerance)
    if (.not. valid) return
    valid = geometry_values_are_valid(self, tolerance)
    if (.not. valid) return
    valid = cell_metrics_are_valid(self, tolerance)
    if (.not. valid) return
    valid = geometric_identity_is_valid(self, tolerance)
  end function eb_geometry_3d_is_valid

  pure logical function geometry_arrays_are_valid(self) result(valid)
    class(eb_geometry_3d), intent(in) :: self

    valid = allocated(self%volume_fraction) .and. &
      allocated(self%cell_centroid_x) .and. &
      allocated(self%cell_centroid_y) .and. &
      allocated(self%cell_centroid_z) .and. &
      allocated(self%x_face_fraction) .and. &
      allocated(self%y_face_fraction) .and. &
      allocated(self%z_face_fraction) .and. &
      allocated(self%x_face_centroid_y) .and. &
      allocated(self%x_face_centroid_z) .and. &
      allocated(self%y_face_centroid_x) .and. &
      allocated(self%y_face_centroid_z) .and. &
      allocated(self%z_face_centroid_x) .and. &
      allocated(self%z_face_centroid_y) .and. &
      allocated(self%boundary_area) .and. &
      allocated(self%boundary_centroid_x) .and. &
      allocated(self%boundary_centroid_y) .and. &
      allocated(self%boundary_centroid_z) .and. &
      allocated(self%boundary_normal_x) .and. &
      allocated(self%boundary_normal_y) .and. &
      allocated(self%boundary_normal_z) .and. &
      allocated(self%boundary_normal_integral_x) .and. &
      allocated(self%boundary_normal_integral_y) .and. &
      allocated(self%boundary_normal_integral_z) .and. &
      allocated(self%cell_type)
    if (.not. valid) return
    valid = all(shape(self%volume_fraction) == &
        [self%nx, self%ny, self%nz]) .and. &
      all(shape(self%cell_centroid_x) == [self%nx, self%ny, self%nz]) .and. &
      all(shape(self%cell_centroid_y) == [self%nx, self%ny, self%nz]) .and. &
      all(shape(self%cell_centroid_z) == [self%nx, self%ny, self%nz]) .and. &
      all(shape(self%x_face_fraction) == &
        [self%nx + 1, self%ny, self%nz]) .and. &
      all(shape(self%y_face_fraction) == &
        [self%nx, self%ny + 1, self%nz]) .and. &
      all(shape(self%z_face_fraction) == &
        [self%nx, self%ny, self%nz + 1]) .and. &
      all(shape(self%x_face_centroid_y) == &
        [self%nx + 1, self%ny, self%nz]) .and. &
      all(shape(self%x_face_centroid_z) == &
        [self%nx + 1, self%ny, self%nz]) .and. &
      all(shape(self%y_face_centroid_x) == &
        [self%nx, self%ny + 1, self%nz]) .and. &
      all(shape(self%y_face_centroid_z) == &
        [self%nx, self%ny + 1, self%nz]) .and. &
      all(shape(self%z_face_centroid_x) == &
        [self%nx, self%ny, self%nz + 1]) .and. &
      all(shape(self%z_face_centroid_y) == &
        [self%nx, self%ny, self%nz + 1]) .and. &
      all(shape(self%boundary_area) == [self%nx, self%ny, self%nz]) .and. &
      all(shape(self%boundary_centroid_x) == &
        [self%nx, self%ny, self%nz]) .and. &
      all(shape(self%boundary_centroid_y) == &
        [self%nx, self%ny, self%nz]) .and. &
      all(shape(self%boundary_centroid_z) == &
        [self%nx, self%ny, self%nz]) .and. &
      all(shape(self%boundary_normal_x) == &
        [self%nx, self%ny, self%nz]) .and. &
      all(shape(self%boundary_normal_y) == &
        [self%nx, self%ny, self%nz]) .and. &
      all(shape(self%boundary_normal_z) == &
        [self%nx, self%ny, self%nz]) .and. &
      all(shape(self%boundary_normal_integral_x) == &
        [self%nx, self%ny, self%nz]) .and. &
      all(shape(self%boundary_normal_integral_y) == &
        [self%nx, self%ny, self%nz]) .and. &
      all(shape(self%boundary_normal_integral_z) == &
        [self%nx, self%ny, self%nz]) .and. &
      all(shape(self%cell_type) == [self%nx, self%ny, self%nz])
    if (.not. valid) return
    valid = all(lbound(self%volume_fraction) == [1, 1, 1]) .and. &
      all(lbound(self%cell_centroid_x) == [1, 1, 1]) .and. &
      all(lbound(self%cell_centroid_y) == [1, 1, 1]) .and. &
      all(lbound(self%cell_centroid_z) == [1, 1, 1]) .and. &
      all(lbound(self%x_face_fraction) == [0, 1, 1]) .and. &
      all(lbound(self%y_face_fraction) == [1, 0, 1]) .and. &
      all(lbound(self%z_face_fraction) == [1, 1, 0]) .and. &
      all(lbound(self%x_face_centroid_y) == [0, 1, 1]) .and. &
      all(lbound(self%x_face_centroid_z) == [0, 1, 1]) .and. &
      all(lbound(self%y_face_centroid_x) == [1, 0, 1]) .and. &
      all(lbound(self%y_face_centroid_z) == [1, 0, 1]) .and. &
      all(lbound(self%z_face_centroid_x) == [1, 1, 0]) .and. &
      all(lbound(self%z_face_centroid_y) == [1, 1, 0]) .and. &
      all(lbound(self%boundary_area) == [1, 1, 1]) .and. &
      all(lbound(self%boundary_centroid_x) == [1, 1, 1]) .and. &
      all(lbound(self%boundary_centroid_y) == [1, 1, 1]) .and. &
      all(lbound(self%boundary_centroid_z) == [1, 1, 1]) .and. &
      all(lbound(self%boundary_normal_x) == [1, 1, 1]) .and. &
      all(lbound(self%boundary_normal_y) == [1, 1, 1]) .and. &
      all(lbound(self%boundary_normal_z) == [1, 1, 1]) .and. &
      all(lbound(self%boundary_normal_integral_x) == [1, 1, 1]) .and. &
      all(lbound(self%boundary_normal_integral_y) == [1, 1, 1]) .and. &
      all(lbound(self%boundary_normal_integral_z) == [1, 1, 1]) .and. &
      all(lbound(self%cell_type) == [1, 1, 1])
  end function geometry_arrays_are_valid

  pure logical function spacing_is_valid(self, tolerance) result(valid)
    class(eb_geometry_3d), intent(in) :: self
    real(dp), intent(in) :: tolerance

    valid = abs(real(self%nx, dp) * self%dx - &
        (self%x_upper - self%x_lower)) <= tolerance * &
          max(1.0_dp, abs(self%x_lower), abs(self%x_upper)) .and. &
      abs(real(self%ny, dp) * self%dy - &
        (self%y_upper - self%y_lower)) <= tolerance * &
          max(1.0_dp, abs(self%y_lower), abs(self%y_upper)) .and. &
      abs(real(self%nz, dp) * self%dz - &
        (self%z_upper - self%z_lower)) <= tolerance * &
          max(1.0_dp, abs(self%z_lower), abs(self%z_upper))
  end function spacing_is_valid

  pure logical function geometry_values_are_valid(self, tolerance) &
      result(valid)
    class(eb_geometry_3d), intent(in) :: self
    real(dp), intent(in) :: tolerance

    valid = all(ieee_is_finite(self%volume_fraction)) .and. &
      all(ieee_is_finite(self%cell_centroid_x)) .and. &
      all(ieee_is_finite(self%cell_centroid_y)) .and. &
      all(ieee_is_finite(self%cell_centroid_z)) .and. &
      all(ieee_is_finite(self%x_face_fraction)) .and. &
      all(ieee_is_finite(self%y_face_fraction)) .and. &
      all(ieee_is_finite(self%z_face_fraction)) .and. &
      all(ieee_is_finite(self%x_face_centroid_y)) .and. &
      all(ieee_is_finite(self%x_face_centroid_z)) .and. &
      all(ieee_is_finite(self%y_face_centroid_x)) .and. &
      all(ieee_is_finite(self%y_face_centroid_z)) .and. &
      all(ieee_is_finite(self%z_face_centroid_x)) .and. &
      all(ieee_is_finite(self%z_face_centroid_y)) .and. &
      all(ieee_is_finite(self%boundary_area)) .and. &
      all(ieee_is_finite(self%boundary_centroid_x)) .and. &
      all(ieee_is_finite(self%boundary_centroid_y)) .and. &
      all(ieee_is_finite(self%boundary_centroid_z)) .and. &
      all(ieee_is_finite(self%boundary_normal_x)) .and. &
      all(ieee_is_finite(self%boundary_normal_y)) .and. &
      all(ieee_is_finite(self%boundary_normal_z)) .and. &
      all(ieee_is_finite(self%boundary_normal_integral_x)) .and. &
      all(ieee_is_finite(self%boundary_normal_integral_y)) .and. &
      all(ieee_is_finite(self%boundary_normal_integral_z)) .and. &
      minval(self%volume_fraction) >= -tolerance .and. &
      maxval(self%volume_fraction) <= 1.0_dp + tolerance .and. &
      minval(self%x_face_fraction) >= -tolerance .and. &
      maxval(self%x_face_fraction) <= 1.0_dp + tolerance .and. &
      minval(self%y_face_fraction) >= -tolerance .and. &
      maxval(self%y_face_fraction) <= 1.0_dp + tolerance .and. &
      minval(self%z_face_fraction) >= -tolerance .and. &
      maxval(self%z_face_fraction) <= 1.0_dp + tolerance .and. &
      maxval(abs(self%cell_centroid_x)) <= 0.5_dp + tolerance .and. &
      maxval(abs(self%cell_centroid_y)) <= 0.5_dp + tolerance .and. &
      maxval(abs(self%cell_centroid_z)) <= 0.5_dp + tolerance .and. &
      maxval(abs(self%x_face_centroid_y)) <= 0.5_dp + tolerance .and. &
      maxval(abs(self%x_face_centroid_z)) <= 0.5_dp + tolerance .and. &
      maxval(abs(self%y_face_centroid_x)) <= 0.5_dp + tolerance .and. &
      maxval(abs(self%y_face_centroid_z)) <= 0.5_dp + tolerance .and. &
      maxval(abs(self%z_face_centroid_x)) <= 0.5_dp + tolerance .and. &
      maxval(abs(self%z_face_centroid_y)) <= 0.5_dp + tolerance .and. &
      minval(self%boundary_area) >= 0.0_dp .and. &
      all(self%cell_type >= eb_covered_cell_3d) .and. &
      all(self%cell_type <= eb_regular_cell_3d)
  end function geometry_values_are_valid

  pure logical function cell_metrics_are_valid(self, tolerance) result(valid)
    class(eb_geometry_3d), intent(in) :: self
    real(dp), intent(in) :: tolerance

    real(dp) :: area_tolerance, coordinate_tolerance, normal_norm
    real(dp) :: x0, x1, y0, y1, z0, z1
    integer :: i, j, k

    valid = .true.
    do k = 1, self%nz
      do j = 1, self%ny
        do i = 1, self%nx
          if (((self%volume_fraction(i, j, k) <= &
                  classification_tolerance) .neqv. &
                (self%cell_type(i, j, k) == eb_covered_cell_3d)) .or. &
              ((self%volume_fraction(i, j, k) >= &
                  1.0_dp - classification_tolerance) .neqv. &
                (self%cell_type(i, j, k) == eb_regular_cell_3d))) then
            valid = .false.
            return
          end if
          if (self%cell_type(i, j, k) /= eb_cut_cell_3d) then
            if (abs(self%cell_centroid_x(i, j, k)) > 8.0_dp * tolerance &
                .or. abs(self%cell_centroid_y(i, j, k)) > &
                  8.0_dp * tolerance .or. &
                abs(self%cell_centroid_z(i, j, k)) > &
                  8.0_dp * tolerance .or. &
                any(abs([ &
                  self%boundary_area(i, j, k), &
                  self%boundary_centroid_x(i, j, k), &
                  self%boundary_centroid_y(i, j, k), &
                  self%boundary_centroid_z(i, j, k), &
                  self%boundary_normal_x(i, j, k), &
                  self%boundary_normal_y(i, j, k), &
                  self%boundary_normal_z(i, j, k), &
                  self%boundary_normal_integral_x(i, j, k), &
                  self%boundary_normal_integral_y(i, j, k), &
                  self%boundary_normal_integral_z(i, j, k)]) > &
                    8.0_dp * tolerance)) then
              valid = .false.
              return
            end if
            cycle
          end if

          if (self%boundary_area(i, j, k) <= 0.0_dp) then
            valid = .false.
            return
          end if
          normal_norm = sqrt( &
            self%boundary_normal_x(i, j, k)**2 + &
            self%boundary_normal_y(i, j, k)**2 + &
            self%boundary_normal_z(i, j, k)**2)
          area_tolerance = 32.0_dp * tolerance * max( &
            self%dx * self%dy, self%dx * self%dz, self%dy * self%dz)
          if (abs(normal_norm - 1.0_dp) > 8.0_dp * tolerance .or. &
              abs(self%boundary_normal_integral_x(i, j, k) - &
                self%boundary_area(i, j, k) * &
                  self%boundary_normal_x(i, j, k)) > area_tolerance .or. &
              abs(self%boundary_normal_integral_y(i, j, k) - &
                self%boundary_area(i, j, k) * &
                  self%boundary_normal_y(i, j, k)) > area_tolerance .or. &
              abs(self%boundary_normal_integral_z(i, j, k) - &
                self%boundary_area(i, j, k) * &
                  self%boundary_normal_z(i, j, k)) > area_tolerance) then
            valid = .false.
            return
          end if
          x0 = self%x_lower + real(i - 1, dp) * self%dx
          x1 = x0 + self%dx
          y0 = self%y_lower + real(j - 1, dp) * self%dy
          y1 = y0 + self%dy
          z0 = self%z_lower + real(k - 1, dp) * self%dz
          z1 = z0 + self%dz
          coordinate_tolerance = 16.0_dp * tolerance * max( &
            1.0_dp, abs(x0), abs(x1), abs(y0), abs(y1), abs(z0), abs(z1))
          if (self%boundary_centroid_x(i, j, k) < &
                x0 - coordinate_tolerance .or. &
              self%boundary_centroid_x(i, j, k) > &
                x1 + coordinate_tolerance .or. &
              self%boundary_centroid_y(i, j, k) < &
                y0 - coordinate_tolerance .or. &
              self%boundary_centroid_y(i, j, k) > &
                y1 + coordinate_tolerance .or. &
              self%boundary_centroid_z(i, j, k) < &
                z0 - coordinate_tolerance .or. &
              self%boundary_centroid_z(i, j, k) > &
                z1 + coordinate_tolerance) then
            valid = .false.
            return
          end if
        end do
      end do
    end do
  end function cell_metrics_are_valid

  pure logical function geometric_identity_is_valid(self, tolerance) &
      result(valid)
    class(eb_geometry_3d), intent(in) :: self
    real(dp), intent(in) :: tolerance

    real(dp) :: metric_tolerance, normal_integral(3)
    integer :: i, j, k

    metric_tolerance = 64.0_dp * tolerance * max( &
      self%dx * self%dy, self%dx * self%dz, self%dy * self%dz)
    valid = .true.
    do k = 1, self%nz
      do j = 1, self%ny
        do i = 1, self%nx
          normal_integral = [ &
            self%dy * self%dz * ( &
              self%x_face_fraction(i, j, k) - &
              self%x_face_fraction(i - 1, j, k)), &
            self%dx * self%dz * ( &
              self%y_face_fraction(i, j, k) - &
              self%y_face_fraction(i, j - 1, k)), &
            self%dx * self%dy * ( &
              self%z_face_fraction(i, j, k) - &
              self%z_face_fraction(i, j, k - 1))]
          if (any(abs(normal_integral - [ &
              self%boundary_normal_integral_x(i, j, k), &
              self%boundary_normal_integral_y(i, j, k), &
              self%boundary_normal_integral_z(i, j, k)]) > &
                metric_tolerance)) then
            valid = .false.
            return
          end if
        end do
      end do
    end do
  end function geometric_identity_is_valid

end module eb_geometry_3d_mod
