module eb_geometry_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  implicit none
  private

  integer, parameter, public :: eb_covered_cell = 0
  integer, parameter, public :: eb_cut_cell = 1
  integer, parameter, public :: eb_regular_cell = 2

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
    real(dp), allocatable :: x_face_fraction(:, :)
    real(dp), allocatable :: y_face_fraction(:, :)
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

    real(dp), parameter :: classification_tolerance = &
      128.0_dp * epsilon(1.0_dp)
    real(dp) :: fraction
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
    allocate(geometry%cell_type(1:geometry%nx, 1:geometry%ny))
    allocate(geometry%x_face_fraction(0:geometry%nx, 1:geometry%ny))
    allocate(geometry%y_face_fraction(1:geometry%nx, 0:geometry%ny))

    do j = 1, geometry%ny
      do i = 0, geometry%nx
        geometry%x_face_fraction(i, j) = positive_segment_fraction( &
          node_level_set(i, j - 1), node_level_set(i, j))
      end do
    end do
    do j = 0, geometry%ny
      do i = 1, geometry%nx
        geometry%y_face_fraction(i, j) = positive_segment_fraction( &
          node_level_set(i - 1, j), node_level_set(i, j))
      end do
    end do

    do j = 1, geometry%ny
      do i = 1, geometry%nx
        fraction = cell_positive_fraction( &
          node_level_set(i - 1, j - 1), &
          node_level_set(i, j - 1), &
          node_level_set(i, j), &
          node_level_set(i - 1, j))
        geometry%volume_fraction(i, j) = fraction
        if (fraction <= classification_tolerance) then
          geometry%cell_type(i, j) = eb_covered_cell
        else if (fraction >= 1.0_dp - classification_tolerance) then
          geometry%cell_type(i, j) = eb_regular_cell
        else
          geometry%cell_type(i, j) = eb_cut_cell
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
      allocated(self%x_face_fraction) .and. &
      allocated(self%y_face_fraction) .and. allocated(self%cell_type)
    if (.not. valid) return
    valid = all(shape(self%volume_fraction) == [self%nx, self%ny]) .and. &
      all(shape(self%x_face_fraction) == [self%nx + 1, self%ny]) .and. &
      all(shape(self%y_face_fraction) == [self%nx, self%ny + 1]) .and. &
      all(shape(self%cell_type) == [self%nx, self%ny]) .and. &
      all(lbound(self%x_face_fraction) == [0, 1]) .and. &
      all(lbound(self%y_face_fraction) == [1, 0])
    if (.not. valid) return
    valid = all(ieee_is_finite(self%volume_fraction)) .and. &
      all(ieee_is_finite(self%x_face_fraction)) .and. &
      all(ieee_is_finite(self%y_face_fraction)) .and. &
      minval(self%volume_fraction) >= -tolerance .and. &
      maxval(self%volume_fraction) <= 1.0_dp + tolerance .and. &
      minval(self%x_face_fraction) >= -tolerance .and. &
      maxval(self%x_face_fraction) <= 1.0_dp + tolerance .and. &
      minval(self%y_face_fraction) >= -tolerance .and. &
      maxval(self%y_face_fraction) <= 1.0_dp + tolerance .and. &
      all(self%cell_type >= eb_covered_cell) .and. &
      all(self%cell_type <= eb_regular_cell)
    if (.not. valid) return
    valid = all((self%volume_fraction >= 1.0_dp - tolerance) .eqv. &
      (self%cell_type == eb_regular_cell)) .and. &
      all((self%volume_fraction <= tolerance) .eqv. &
        (self%cell_type == eb_covered_cell))
  end function eb_geometry_is_valid

  pure real(dp) function positive_segment_fraction( &
      first_value, second_value) result(fraction)
    real(dp), intent(in) :: first_value, second_value

    real(dp) :: crossing

    if (first_value > 0.0_dp .and. second_value > 0.0_dp) then
      fraction = 1.0_dp
    else if (first_value <= 0.0_dp .and. second_value <= 0.0_dp) then
      fraction = 0.0_dp
    else
      crossing = first_value / (first_value - second_value)
      if (first_value > 0.0_dp) then
        fraction = crossing
      else
        fraction = 1.0_dp - crossing
      end if
      fraction = min(1.0_dp, max(0.0_dp, fraction))
    end if
  end function positive_segment_fraction

  pure real(dp) function cell_positive_fraction( &
      lower_left, lower_right, upper_right, upper_left) result(fraction)
    real(dp), intent(in) :: lower_left, lower_right, upper_right, upper_left

    fraction = positive_triangle_area( &
      [0.0_dp, 1.0_dp, 1.0_dp], [0.0_dp, 0.0_dp, 1.0_dp], &
      [lower_left, lower_right, upper_right]) + &
      positive_triangle_area( &
        [0.0_dp, 1.0_dp, 0.0_dp], [0.0_dp, 1.0_dp, 1.0_dp], &
        [lower_left, upper_right, upper_left])
    fraction = min(1.0_dp, max(0.0_dp, fraction))
  end function cell_positive_fraction

  pure real(dp) function positive_triangle_area( &
      vertex_x, vertex_y, vertex_level_set) result(area)
    real(dp), intent(in) :: vertex_x(3), vertex_y(3)
    real(dp), intent(in) :: vertex_level_set(3)

    real(dp) :: clipped_x(4), clipped_y(4)
    real(dp) :: first_x, first_y, second_x, second_y
    real(dp) :: first_value, second_value, crossing
    integer :: edge, next_edge, count, vertex

    if (maxval(vertex_level_set) <= 0.0_dp) then
      area = 0.0_dp
      return
    else if (minval(vertex_level_set) > 0.0_dp) then
      area = 0.5_dp
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

    area = 0.0_dp
    if (count < 3) return
    do vertex = 1, count
      next_edge = merge(vertex + 1, 1, vertex < count)
      area = area + clipped_x(vertex) * clipped_y(next_edge) - &
        clipped_y(vertex) * clipped_x(next_edge)
    end do
    area = min(0.5_dp, max(0.0_dp, 0.5_dp * abs(area)))
  end function positive_triangle_area

end module eb_geometry_2d_mod
