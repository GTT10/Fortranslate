module eb_geometry_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  implicit none
  private

  integer, parameter, public :: eb_covered_cell = 0
  integer, parameter, public :: eb_cut_cell = 1
  integer, parameter, public :: eb_regular_cell = 2
  real(dp), parameter :: eb_classification_tolerance = &
    128.0_dp * epsilon(1.0_dp)

  type, public :: eb_geometry_2d
    integer :: nx = 0
    integer :: ny = 0
    real(dp) :: x_lower = 0.0_dp
    real(dp) :: x_upper = 0.0_dp
    real(dp) :: y_lower = 0.0_dp
    real(dp) :: y_upper = 0.0_dp
    real(dp) :: dx = 0.0_dp
    real(dp) :: dy = 0.0_dp
    real(dp), allocatable :: volume_fraction(:, :)
    ! Fluid-volume centroid offsets relative to the Cartesian cell center,
    ! normalized by dx and dy respectively (AMReX convention).
    real(dp), allocatable :: cell_centroid_x(:, :)
    real(dp), allocatable :: cell_centroid_y(:, :)
    real(dp), allocatable :: x_face_fraction(:, :)
    real(dp), allocatable :: y_face_fraction(:, :)
    ! Tangential face-centroid offsets relative to the Cartesian face center,
    ! normalized by dy on x-faces and dx on y-faces (AMReX convention).
    real(dp), allocatable :: x_face_centroid_y(:, :)
    real(dp), allocatable :: y_face_centroid_x(:, :)
    real(dp), allocatable :: boundary_length(:, :)
    real(dp), allocatable :: boundary_centroid_x(:, :)
    real(dp), allocatable :: boundary_centroid_y(:, :)
    real(dp), allocatable :: boundary_normal_x(:, :)
    real(dp), allocatable :: boundary_normal_y(:, :)
    real(dp), allocatable :: boundary_normal_integral_x(:, :)
    real(dp), allocatable :: boundary_normal_integral_y(:, :)
    integer, allocatable :: cell_type(:, :)
  contains
    procedure :: is_valid => eb_geometry_is_valid
  end type eb_geometry_2d

  public :: build_eb_geometry_2d

contains

  subroutine build_eb_geometry_2d( &
      node_level_set, x_lower, x_upper, y_lower, y_upper, geometry, ok)
    real(dp), intent(in) :: node_level_set(0:, 0:)
    real(dp), intent(in) :: x_lower, x_upper, y_lower, y_upper
    type(eb_geometry_2d), intent(out) :: geometry
    logical, intent(out) :: ok

    real(dp) :: fraction, centroid_x, centroid_y, x0, x1, y0, y1
    logical :: metrics_ok
    integer :: i, j

    geometry = eb_geometry_2d()
    ok = .false.
    if (size(node_level_set, 1) < 2 .or. &
        size(node_level_set, 2) < 2 .or. &
        x_upper <= x_lower .or. y_upper <= y_lower .or. &
        .not. all(ieee_is_finite(node_level_set))) return

    geometry%nx = size(node_level_set, 1) - 1
    geometry%ny = size(node_level_set, 2) - 1
    geometry%x_lower = x_lower
    geometry%x_upper = x_upper
    geometry%y_lower = y_lower
    geometry%y_upper = y_upper
    geometry%dx = (x_upper - x_lower) / real(geometry%nx, dp)
    geometry%dy = (y_upper - y_lower) / real(geometry%ny, dp)
    allocate(geometry%volume_fraction(1:geometry%nx, 1:geometry%ny))
    allocate(geometry%cell_centroid_x(1:geometry%nx, 1:geometry%ny))
    allocate(geometry%cell_centroid_y(1:geometry%nx, 1:geometry%ny))
    allocate(geometry%cell_type(1:geometry%nx, 1:geometry%ny))
    allocate(geometry%x_face_fraction(0:geometry%nx, 1:geometry%ny))
    allocate(geometry%y_face_fraction(1:geometry%nx, 0:geometry%ny))
    allocate(geometry%x_face_centroid_y(0:geometry%nx, 1:geometry%ny))
    allocate(geometry%y_face_centroid_x(1:geometry%nx, 0:geometry%ny))
    allocate(geometry%boundary_length(1:geometry%nx, 1:geometry%ny))
    allocate(geometry%boundary_centroid_x(1:geometry%nx, 1:geometry%ny))
    allocate(geometry%boundary_centroid_y(1:geometry%nx, 1:geometry%ny))
    allocate(geometry%boundary_normal_x(1:geometry%nx, 1:geometry%ny))
    allocate(geometry%boundary_normal_y(1:geometry%nx, 1:geometry%ny))
    allocate(geometry%boundary_normal_integral_x( &
      1:geometry%nx, 1:geometry%ny))
    allocate(geometry%boundary_normal_integral_y( &
      1:geometry%nx, 1:geometry%ny))
    geometry%boundary_length = 0.0_dp
    geometry%boundary_centroid_x = 0.0_dp
    geometry%boundary_centroid_y = 0.0_dp
    geometry%boundary_normal_x = 0.0_dp
    geometry%boundary_normal_y = 0.0_dp
    geometry%boundary_normal_integral_x = 0.0_dp
    geometry%boundary_normal_integral_y = 0.0_dp

    do j = 1, geometry%ny
      do i = 0, geometry%nx
        call positive_segment_metrics( &
          node_level_set(i, j - 1), node_level_set(i, j), &
          geometry%x_face_fraction(i, j), &
          geometry%x_face_centroid_y(i, j))
      end do
    end do
    do j = 0, geometry%ny
      do i = 1, geometry%nx
        call positive_segment_metrics( &
          node_level_set(i - 1, j), node_level_set(i, j), &
          geometry%y_face_fraction(i, j), &
          geometry%y_face_centroid_x(i, j))
      end do
    end do

    do j = 1, geometry%ny
      do i = 1, geometry%nx
        call cell_positive_metrics( &
          node_level_set(i - 1, j - 1), &
          node_level_set(i, j - 1), &
          node_level_set(i, j), &
          node_level_set(i - 1, j), fraction, centroid_x, centroid_y)
        geometry%volume_fraction(i, j) = fraction
        geometry%cell_centroid_x(i, j) = centroid_x
        geometry%cell_centroid_y(i, j) = centroid_y
        if (fraction <= eb_classification_tolerance) then
          geometry%cell_type(i, j) = eb_covered_cell
        else if (fraction >= 1.0_dp - eb_classification_tolerance) then
          geometry%cell_type(i, j) = eb_regular_cell
        else
          geometry%cell_type(i, j) = eb_cut_cell
          x0 = x_lower + real(i - 1, dp) * geometry%dx
          x1 = x0 + geometry%dx
          y0 = y_lower + real(j - 1, dp) * geometry%dy
          y1 = y0 + geometry%dy
          call cell_interface_metrics( &
            [x0, x1, x1, x0], [y0, y0, y1, y1], &
            [node_level_set(i - 1, j - 1), &
             node_level_set(i, j - 1), node_level_set(i, j), &
             node_level_set(i - 1, j)], &
            geometry%boundary_length(i, j), &
            geometry%boundary_centroid_x(i, j), &
            geometry%boundary_centroid_y(i, j), &
            geometry%boundary_normal_x(i, j), &
            geometry%boundary_normal_y(i, j), &
            geometry%boundary_normal_integral_x(i, j), &
            geometry%boundary_normal_integral_y(i, j), metrics_ok)
          if (.not. metrics_ok) return
        end if
      end do
    end do

    ok = geometry%is_valid()
  end subroutine build_eb_geometry_2d

  pure logical function eb_geometry_is_valid(self) result(valid)
    class(eb_geometry_2d), intent(in) :: self

    real(dp), parameter :: tolerance = 256.0_dp * epsilon(1.0_dp)

    valid = self%nx >= 1 .and. self%ny >= 1 .and. &
      self%x_upper > self%x_lower .and. self%y_upper > self%y_lower .and. &
      self%dx > 0.0_dp .and. self%dy > 0.0_dp .and. &
      allocated(self%volume_fraction) .and. &
      allocated(self%cell_centroid_x) .and. &
      allocated(self%cell_centroid_y) .and. &
      allocated(self%x_face_fraction) .and. &
      allocated(self%y_face_fraction) .and. &
      allocated(self%x_face_centroid_y) .and. &
      allocated(self%y_face_centroid_x) .and. &
      allocated(self%boundary_length) .and. &
      allocated(self%boundary_centroid_x) .and. &
      allocated(self%boundary_centroid_y) .and. &
      allocated(self%boundary_normal_x) .and. &
      allocated(self%boundary_normal_y) .and. &
      allocated(self%boundary_normal_integral_x) .and. &
      allocated(self%boundary_normal_integral_y) .and. &
      allocated(self%cell_type)
    if (.not. valid) return
    valid = abs(real(self%nx, dp) * self%dx - &
      (self%x_upper - self%x_lower)) <= tolerance * &
        max(1.0_dp, abs(self%x_lower), abs(self%x_upper)) .and. &
      abs(real(self%ny, dp) * self%dy - &
        (self%y_upper - self%y_lower)) <= tolerance * &
          max(1.0_dp, abs(self%y_lower), abs(self%y_upper))
    if (.not. valid) return
    valid = all(shape(self%volume_fraction) == [self%nx, self%ny]) .and. &
      all(shape(self%cell_centroid_x) == [self%nx, self%ny]) .and. &
      all(shape(self%cell_centroid_y) == [self%nx, self%ny]) .and. &
      all(shape(self%x_face_fraction) == [self%nx + 1, self%ny]) .and. &
      all(shape(self%y_face_fraction) == [self%nx, self%ny + 1]) .and. &
      all(shape(self%x_face_centroid_y) == [self%nx + 1, self%ny]) .and. &
      all(shape(self%y_face_centroid_x) == [self%nx, self%ny + 1]) .and. &
      all(shape(self%boundary_length) == [self%nx, self%ny]) .and. &
      all(shape(self%boundary_centroid_x) == [self%nx, self%ny]) .and. &
      all(shape(self%boundary_centroid_y) == [self%nx, self%ny]) .and. &
      all(shape(self%boundary_normal_x) == [self%nx, self%ny]) .and. &
      all(shape(self%boundary_normal_y) == [self%nx, self%ny]) .and. &
      all(shape(self%boundary_normal_integral_x) == &
        [self%nx, self%ny]) .and. &
      all(shape(self%boundary_normal_integral_y) == &
        [self%nx, self%ny]) .and. &
      all(shape(self%cell_type) == [self%nx, self%ny]) .and. &
      all(lbound(self%volume_fraction) == [1, 1]) .and. &
      all(lbound(self%cell_centroid_x) == [1, 1]) .and. &
      all(lbound(self%cell_centroid_y) == [1, 1]) .and. &
      all(lbound(self%boundary_length) == [1, 1]) .and. &
      all(lbound(self%boundary_centroid_x) == [1, 1]) .and. &
      all(lbound(self%boundary_centroid_y) == [1, 1]) .and. &
      all(lbound(self%boundary_normal_x) == [1, 1]) .and. &
      all(lbound(self%boundary_normal_y) == [1, 1]) .and. &
      all(lbound(self%boundary_normal_integral_x) == [1, 1]) .and. &
      all(lbound(self%boundary_normal_integral_y) == [1, 1]) .and. &
      all(lbound(self%cell_type) == [1, 1]) .and. &
      all(lbound(self%x_face_fraction) == [0, 1]) .and. &
      all(lbound(self%y_face_fraction) == [1, 0]) .and. &
      all(lbound(self%x_face_centroid_y) == [0, 1]) .and. &
      all(lbound(self%y_face_centroid_x) == [1, 0])
    if (.not. valid) return
    valid = all(ieee_is_finite(self%volume_fraction)) .and. &
      all(ieee_is_finite(self%cell_centroid_x)) .and. &
      all(ieee_is_finite(self%cell_centroid_y)) .and. &
      all(ieee_is_finite(self%x_face_fraction)) .and. &
      all(ieee_is_finite(self%y_face_fraction)) .and. &
      all(ieee_is_finite(self%x_face_centroid_y)) .and. &
      all(ieee_is_finite(self%y_face_centroid_x)) .and. &
      all(ieee_is_finite(self%boundary_length)) .and. &
      all(ieee_is_finite(self%boundary_centroid_x)) .and. &
      all(ieee_is_finite(self%boundary_centroid_y)) .and. &
      all(ieee_is_finite(self%boundary_normal_x)) .and. &
      all(ieee_is_finite(self%boundary_normal_y)) .and. &
      all(ieee_is_finite(self%boundary_normal_integral_x)) .and. &
      all(ieee_is_finite(self%boundary_normal_integral_y)) .and. &
      minval(self%volume_fraction) >= -tolerance .and. &
      maxval(self%volume_fraction) <= 1.0_dp + tolerance .and. &
      maxval(abs(self%cell_centroid_x)) <= 0.5_dp + tolerance .and. &
      maxval(abs(self%cell_centroid_y)) <= 0.5_dp + tolerance .and. &
      minval(self%x_face_fraction) >= -tolerance .and. &
      maxval(self%x_face_fraction) <= 1.0_dp + tolerance .and. &
      minval(self%y_face_fraction) >= -tolerance .and. &
      maxval(self%y_face_fraction) <= 1.0_dp + tolerance .and. &
      maxval(abs(self%x_face_centroid_y)) <= 0.5_dp + tolerance .and. &
      maxval(abs(self%y_face_centroid_x)) <= 0.5_dp + tolerance .and. &
      minval(self%boundary_length) >= 0.0_dp .and. &
      all(self%cell_type >= eb_covered_cell) .and. &
      all(self%cell_type <= eb_regular_cell)
    if (.not. valid) return
    valid = all(((self%cell_type == eb_cut_cell) .or. &
      (abs(self%cell_centroid_x) <= 8.0_dp * tolerance .and. &
       abs(self%cell_centroid_y) <= 8.0_dp * tolerance)))
    if (.not. valid) return
    valid = all((self%volume_fraction >= &
      1.0_dp - eb_classification_tolerance) .eqv. &
      (self%cell_type == eb_regular_cell)) .and. &
      all((self%volume_fraction <= eb_classification_tolerance) .eqv. &
        (self%cell_type == eb_covered_cell))
    if (.not. valid) return
    valid = validate_face_centroids(self, tolerance)
    if (.not. valid) return
    valid = validate_interface_metrics(self, tolerance)
  end function eb_geometry_is_valid

  pure logical function validate_face_centroids(self, tolerance) result(valid)
    class(eb_geometry_2d), intent(in) :: self
    real(dp), intent(in) :: tolerance

    real(dp) :: expected
    integer :: i, j

    valid = .true.
    do j = 1, self%ny
      do i = 0, self%nx
        expected = 0.5_dp * (1.0_dp - self%x_face_fraction(i, j))
        if (abs(abs(self%x_face_centroid_y(i, j)) - expected) > &
            8.0_dp * tolerance .and. &
            self%x_face_fraction(i, j) > 0.0_dp .and. &
            self%x_face_fraction(i, j) < 1.0_dp) then
          valid = .false.
          return
        end if
        if ((self%x_face_fraction(i, j) == 0.0_dp .or. &
             self%x_face_fraction(i, j) == 1.0_dp) .and. &
            abs(self%x_face_centroid_y(i, j)) > 8.0_dp * tolerance) then
          valid = .false.
          return
        end if
      end do
    end do
    do j = 0, self%ny
      do i = 1, self%nx
        expected = 0.5_dp * (1.0_dp - self%y_face_fraction(i, j))
        if (abs(abs(self%y_face_centroid_x(i, j)) - expected) > &
            8.0_dp * tolerance .and. &
            self%y_face_fraction(i, j) > 0.0_dp .and. &
            self%y_face_fraction(i, j) < 1.0_dp) then
          valid = .false.
          return
        end if
        if ((self%y_face_fraction(i, j) == 0.0_dp .or. &
             self%y_face_fraction(i, j) == 1.0_dp) .and. &
            abs(self%y_face_centroid_x(i, j)) > 8.0_dp * tolerance) then
          valid = .false.
          return
        end if
      end do
    end do
  end function validate_face_centroids

  pure logical function validate_interface_metrics( &
      self, tolerance) result(valid)
    class(eb_geometry_2d), intent(in) :: self
    real(dp), intent(in) :: tolerance

    real(dp) :: x0, x1, y0, y1, normal_norm, normal_integral_norm
    real(dp) :: coordinate_tolerance, metric_tolerance
    integer :: i, j

    valid = .true.
    do j = 1, self%ny
      do i = 1, self%nx
        if (self%cell_type(i, j) == eb_cut_cell) then
          if (self%boundary_length(i, j) <= 0.0_dp) then
            valid = .false.
            return
          end if
          normal_norm = sqrt(self%boundary_normal_x(i, j)**2 + &
            self%boundary_normal_y(i, j)**2)
          if (abs(normal_norm - 1.0_dp) > 8.0_dp * tolerance) then
            valid = .false.
            return
          end if
          normal_integral_norm = sqrt( &
            self%boundary_normal_integral_x(i, j)**2 + &
            self%boundary_normal_integral_y(i, j)**2)
          metric_tolerance = 32.0_dp * tolerance * max(self%dx, self%dy)
          if (normal_integral_norm <= 0.0_dp .or. &
              normal_integral_norm > &
                self%boundary_length(i, j) + metric_tolerance .or. &
              abs(self%boundary_normal_integral_x(i, j) - &
                normal_integral_norm * self%boundary_normal_x(i, j)) > &
                  metric_tolerance .or. &
              abs(self%boundary_normal_integral_y(i, j) - &
                normal_integral_norm * self%boundary_normal_y(i, j)) > &
                  metric_tolerance) then
            valid = .false.
            return
          end if
          x0 = self%x_lower + real(i - 1, dp) * self%dx
          x1 = x0 + self%dx
          y0 = self%y_lower + real(j - 1, dp) * self%dy
          y1 = y0 + self%dy
          coordinate_tolerance = 8.0_dp * tolerance * &
            max(self%dx, self%dy, abs(x0), abs(x1), abs(y0), abs(y1))
          if (self%boundary_centroid_x(i, j) < x0 - coordinate_tolerance &
              .or. self%boundary_centroid_x(i, j) > &
                x1 + coordinate_tolerance .or. &
              self%boundary_centroid_y(i, j) < y0 - coordinate_tolerance &
              .or. self%boundary_centroid_y(i, j) > &
                y1 + coordinate_tolerance) then
            valid = .false.
            return
          end if
        else if (self%boundary_length(i, j) /= 0.0_dp .or. &
            self%boundary_centroid_x(i, j) /= 0.0_dp .or. &
            self%boundary_centroid_y(i, j) /= 0.0_dp .or. &
            self%boundary_normal_x(i, j) /= 0.0_dp .or. &
            self%boundary_normal_y(i, j) /= 0.0_dp .or. &
            self%boundary_normal_integral_x(i, j) /= 0.0_dp .or. &
            self%boundary_normal_integral_y(i, j) /= 0.0_dp) then
          valid = .false.
          return
        end if
      end do
    end do
  end function validate_interface_metrics

  pure subroutine cell_interface_metrics( &
      vertex_x, vertex_y, vertex_level_set, length, centroid_x, centroid_y, &
      normal_x, normal_y, normal_integral_x, normal_integral_y, ok)
    real(dp), intent(in) :: vertex_x(4), vertex_y(4), vertex_level_set(4)
    real(dp), intent(out) :: length, centroid_x, centroid_y
    real(dp), intent(out) :: normal_x, normal_y
    real(dp), intent(out) :: normal_integral_x, normal_integral_y
    logical, intent(out) :: ok

    real(dp) :: segment_length(2), segment_centroid_x(2)
    real(dp) :: segment_centroid_y(2), segment_normal_x(2)
    real(dp) :: segment_normal_y(2), scale, duplicate_tolerance, norm
    logical :: duplicate_interface, segment_ok

    call triangle_interface_metrics( &
      vertex_x([1, 2, 3]), vertex_y([1, 2, 3]), &
      vertex_level_set([1, 2, 3]), segment_length(1), &
      segment_centroid_x(1), segment_centroid_y(1), &
      segment_normal_x(1), segment_normal_y(1), segment_ok)
    if (.not. segment_ok) then
      ok = .false.
      return
    end if
    call triangle_interface_metrics( &
      vertex_x([1, 3, 4]), vertex_y([1, 3, 4]), &
      vertex_level_set([1, 3, 4]), segment_length(2), &
      segment_centroid_x(2), segment_centroid_y(2), &
      segment_normal_x(2), segment_normal_y(2), segment_ok)
    if (.not. segment_ok) then
      ok = .false.
      return
    end if

    scale = max(abs(vertex_x(2) - vertex_x(1)), &
      abs(vertex_y(4) - vertex_y(1)), tiny(1.0_dp))
    duplicate_tolerance = max( &
      256.0_dp * epsilon(1.0_dp) * scale, &
      8.0_dp * spacing(maxval(abs(vertex_x))), &
      8.0_dp * spacing(maxval(abs(vertex_y))))
    duplicate_interface = &
      segment_length(1) > 0.0_dp .and. segment_length(2) > 0.0_dp .and. &
      abs(segment_length(1) - segment_length(2)) <= &
        duplicate_tolerance .and. &
      sqrt((segment_centroid_x(1) - segment_centroid_x(2))**2 + &
        (segment_centroid_y(1) - segment_centroid_y(2))**2) <= &
          duplicate_tolerance
    if (duplicate_interface) then
      length = segment_length(1)
      centroid_x = 0.5_dp * &
        (segment_centroid_x(1) + segment_centroid_x(2))
      centroid_y = 0.5_dp * &
        (segment_centroid_y(1) + segment_centroid_y(2))
      normal_x = segment_normal_x(1) + segment_normal_x(2)
      normal_y = segment_normal_y(1) + segment_normal_y(2)
    else
      length = sum(segment_length)
      if (length <= 0.0_dp) then
        ok = .false.
        return
      end if
      centroid_x = sum(segment_length * segment_centroid_x) / length
      centroid_y = sum(segment_length * segment_centroid_y) / length
      normal_x = sum(segment_length * segment_normal_x)
      normal_y = sum(segment_length * segment_normal_y)
    end if
    norm = sqrt(normal_x**2 + normal_y**2)
    if (norm <= 0.0_dp) then
      ok = .false.
      return
    end if
    if (duplicate_interface) then
      normal_integral_x = length * normal_x / norm
      normal_integral_y = length * normal_y / norm
    else
      normal_integral_x = normal_x
      normal_integral_y = normal_y
    end if
    normal_x = normal_integral_x / &
      sqrt(normal_integral_x**2 + normal_integral_y**2)
    normal_y = normal_integral_y / &
      sqrt(normal_integral_x**2 + normal_integral_y**2)
    ok = .true.
  end subroutine cell_interface_metrics

  pure subroutine triangle_interface_metrics( &
      vertex_x, vertex_y, vertex_level_set, length, centroid_x, centroid_y, &
      normal_x, normal_y, ok)
    real(dp), intent(in) :: vertex_x(3), vertex_y(3), vertex_level_set(3)
    real(dp), intent(out) :: length, centroid_x, centroid_y
    real(dp), intent(out) :: normal_x, normal_y
    logical, intent(out) :: ok

    real(dp) :: point_x(3), point_y(3), crossing, determinant, norm
    real(dp) :: candidate_x, candidate_y, point_tolerance, scale
    integer :: edge, next_edge, point_count

    length = 0.0_dp
    centroid_x = 0.0_dp
    centroid_y = 0.0_dp
    normal_x = 0.0_dp
    normal_y = 0.0_dp
    ok = .true.
    if (all(vertex_level_set == 0.0_dp)) return

    scale = max(maxval(vertex_x) - minval(vertex_x), &
      maxval(vertex_y) - minval(vertex_y), tiny(1.0_dp))
    point_tolerance = max(256.0_dp * epsilon(1.0_dp) * scale, &
      8.0_dp * spacing(maxval(abs(vertex_x))), &
      8.0_dp * spacing(maxval(abs(vertex_y))))
    point_count = 0
    do edge = 1, 3
      next_edge = merge(edge + 1, 1, edge < 3)
      if (vertex_level_set(edge) == 0.0_dp) then
        call append_unique_point(vertex_x(edge), vertex_y(edge), &
          point_x, point_y, point_count, point_tolerance)
      end if
      if ((vertex_level_set(edge) > 0.0_dp .and. &
           vertex_level_set(next_edge) < 0.0_dp) .or. &
          (vertex_level_set(edge) < 0.0_dp .and. &
           vertex_level_set(next_edge) > 0.0_dp)) then
        crossing = vertex_level_set(edge) / &
          (vertex_level_set(edge) - vertex_level_set(next_edge))
        candidate_x = vertex_x(edge) + crossing * &
          (vertex_x(next_edge) - vertex_x(edge))
        candidate_y = vertex_y(edge) + crossing * &
          (vertex_y(next_edge) - vertex_y(edge))
        call append_unique_point(candidate_x, candidate_y, &
          point_x, point_y, point_count, point_tolerance)
      end if
    end do
    if (point_count < 2) return
    if (point_count > 2) then
      ok = .false.
      return
    end if

    length = sqrt((point_x(2) - point_x(1))**2 + &
      (point_y(2) - point_y(1))**2)
    if (length <= point_tolerance) then
      length = 0.0_dp
      return
    end if
    centroid_x = 0.5_dp * (point_x(1) + point_x(2))
    centroid_y = 0.5_dp * (point_y(1) + point_y(2))
    determinant = (vertex_x(2) - vertex_x(1)) * &
      (vertex_y(3) - vertex_y(1)) - &
      (vertex_x(3) - vertex_x(1)) * &
        (vertex_y(2) - vertex_y(1))
    if (abs(determinant) <= tiny(1.0_dp)) then
      ok = .false.
      return
    end if
    normal_x = ( &
      (vertex_level_set(2) - vertex_level_set(1)) * &
        (vertex_y(3) - vertex_y(1)) - &
      (vertex_level_set(3) - vertex_level_set(1)) * &
        (vertex_y(2) - vertex_y(1))) / determinant
    normal_y = ( &
      (vertex_x(2) - vertex_x(1)) * &
        (vertex_level_set(3) - vertex_level_set(1)) - &
      (vertex_x(3) - vertex_x(1)) * &
        (vertex_level_set(2) - vertex_level_set(1))) / determinant
    norm = sqrt(normal_x**2 + normal_y**2)
    if (norm <= 0.0_dp) then
      ok = .false.
      return
    end if
    normal_x = normal_x / norm
    normal_y = normal_y / norm
  end subroutine triangle_interface_metrics

  pure subroutine append_unique_point( &
      x, y, point_x, point_y, point_count, tolerance)
    real(dp), intent(in) :: x, y, tolerance
    real(dp), intent(inout) :: point_x(3), point_y(3)
    integer, intent(inout) :: point_count

    integer :: point

    do point = 1, point_count
      if (sqrt((x - point_x(point))**2 + &
          (y - point_y(point))**2) <= tolerance) return
    end do
    if (point_count >= size(point_x)) return
    point_count = point_count + 1
    point_x(point_count) = x
    point_y(point_count) = y
  end subroutine append_unique_point

  pure subroutine positive_segment_metrics( &
      first_value, second_value, fraction, centroid_offset)
    real(dp), intent(in) :: first_value, second_value
    real(dp), intent(out) :: fraction, centroid_offset

    real(dp) :: crossing

    if (first_value > 0.0_dp .and. second_value > 0.0_dp) then
      fraction = 1.0_dp
      centroid_offset = 0.0_dp
    else if (first_value <= 0.0_dp .and. second_value <= 0.0_dp) then
      fraction = 0.0_dp
      centroid_offset = 0.0_dp
    else
      crossing = first_value / (first_value - second_value)
      if (first_value > 0.0_dp) then
        fraction = crossing
        centroid_offset = -0.5_dp * (1.0_dp - fraction)
      else
        fraction = 1.0_dp - crossing
        centroid_offset = 0.5_dp * (1.0_dp - fraction)
      end if
      fraction = min(1.0_dp, max(0.0_dp, fraction))
      if (fraction <= 0.0_dp .or. fraction >= 1.0_dp) &
        centroid_offset = 0.0_dp
    end if
  end subroutine positive_segment_metrics

  pure subroutine cell_positive_metrics( &
      lower_left, lower_right, upper_right, upper_left, fraction, &
      centroid_offset_x, centroid_offset_y)
    real(dp), intent(in) :: lower_left, lower_right, upper_right, upper_left
    real(dp), intent(out) :: fraction, centroid_offset_x, centroid_offset_y

    real(dp) :: first_area, second_area
    real(dp) :: first_centroid_x, first_centroid_y
    real(dp) :: second_centroid_x, second_centroid_y

    call positive_triangle_metrics( &
      [0.0_dp, 1.0_dp, 1.0_dp], [0.0_dp, 0.0_dp, 1.0_dp], &
      [lower_left, lower_right, upper_right], first_area, &
      first_centroid_x, first_centroid_y)
    call positive_triangle_metrics( &
      [0.0_dp, 1.0_dp, 0.0_dp], [0.0_dp, 1.0_dp, 1.0_dp], &
      [lower_left, upper_right, upper_left], second_area, &
      second_centroid_x, second_centroid_y)
    fraction = first_area + second_area
    fraction = min(1.0_dp, max(0.0_dp, fraction))
    if (fraction > tiny(1.0_dp)) then
      centroid_offset_x = (first_area * first_centroid_x + &
        second_area * second_centroid_x) / fraction - 0.5_dp
      centroid_offset_y = (first_area * first_centroid_y + &
        second_area * second_centroid_y) / fraction - 0.5_dp
      centroid_offset_x = min(0.5_dp, max(-0.5_dp, centroid_offset_x))
      centroid_offset_y = min(0.5_dp, max(-0.5_dp, centroid_offset_y))
    else
      centroid_offset_x = 0.0_dp
      centroid_offset_y = 0.0_dp
    end if
    if (fraction >= 1.0_dp - eb_classification_tolerance) then
      centroid_offset_x = 0.0_dp
      centroid_offset_y = 0.0_dp
    end if
  end subroutine cell_positive_metrics

  pure subroutine positive_triangle_metrics( &
      vertex_x, vertex_y, vertex_level_set, area, centroid_x, centroid_y)
    real(dp), intent(in) :: vertex_x(3), vertex_y(3)
    real(dp), intent(in) :: vertex_level_set(3)
    real(dp), intent(out) :: area, centroid_x, centroid_y

    real(dp) :: clipped_x(4), clipped_y(4)
    real(dp) :: first_x, first_y, second_x, second_y
    real(dp) :: first_value, second_value, crossing, cross
    real(dp) :: signed_area_twice, centroid_numerator_x
    real(dp) :: centroid_numerator_y
    integer :: edge, next_edge, count, vertex

    area = 0.0_dp
    centroid_x = 0.0_dp
    centroid_y = 0.0_dp
    if (maxval(vertex_level_set) <= 0.0_dp) then
      return
    else if (minval(vertex_level_set) > 0.0_dp) then
      signed_area_twice = &
        vertex_x(1) * vertex_y(2) - vertex_y(1) * vertex_x(2) + &
        vertex_x(2) * vertex_y(3) - vertex_y(2) * vertex_x(3) + &
        vertex_x(3) * vertex_y(1) - vertex_y(3) * vertex_x(1)
      area = 0.5_dp * abs(signed_area_twice)
      centroid_x = sum(vertex_x) / 3.0_dp
      centroid_y = sum(vertex_y) / 3.0_dp
      return
    end if

    count = 0
    do edge = 1, 3
      next_edge = merge(edge + 1, 1, edge < 3)
      first_x = vertex_x(edge)
      first_y = vertex_y(edge)
      first_value = vertex_level_set(edge)
      second_x = vertex_x(next_edge)
      second_y = vertex_y(next_edge)
      second_value = vertex_level_set(next_edge)
      if (first_value > 0.0_dp .and. second_value > 0.0_dp) then
        count = count + 1
        clipped_x(count) = second_x
        clipped_y(count) = second_y
      else if (first_value > 0.0_dp .and. second_value <= 0.0_dp) then
        crossing = first_value / (first_value - second_value)
        count = count + 1
        clipped_x(count) = first_x + crossing * (second_x - first_x)
        clipped_y(count) = first_y + crossing * (second_y - first_y)
      else if (first_value <= 0.0_dp .and. second_value > 0.0_dp) then
        crossing = first_value / (first_value - second_value)
        count = count + 1
        clipped_x(count) = first_x + crossing * (second_x - first_x)
        clipped_y(count) = first_y + crossing * (second_y - first_y)
        count = count + 1
        clipped_x(count) = second_x
        clipped_y(count) = second_y
      end if
    end do

    if (count < 3) return
    signed_area_twice = 0.0_dp
    centroid_numerator_x = 0.0_dp
    centroid_numerator_y = 0.0_dp
    do vertex = 1, count
      next_edge = merge(vertex + 1, 1, vertex < count)
      cross = clipped_x(vertex) * clipped_y(next_edge) - &
        clipped_y(vertex) * clipped_x(next_edge)
      signed_area_twice = signed_area_twice + cross
      centroid_numerator_x = centroid_numerator_x + &
        (clipped_x(vertex) + clipped_x(next_edge)) * cross
      centroid_numerator_y = centroid_numerator_y + &
        (clipped_y(vertex) + clipped_y(next_edge)) * cross
    end do
    if (abs(signed_area_twice) <= tiny(1.0_dp)) return
    area = min(0.5_dp, max(0.0_dp, 0.5_dp * abs(signed_area_twice)))
    centroid_x = centroid_numerator_x / (3.0_dp * signed_area_twice)
    centroid_y = centroid_numerator_y / (3.0_dp * signed_area_twice)
  end subroutine positive_triangle_metrics

end module eb_geometry_2d_mod
