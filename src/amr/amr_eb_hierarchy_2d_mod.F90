module amr_eb_hierarchy_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_covered_cell
  implicit none
  private

  real(dp), parameter :: geometry_consistency_tolerance = 1.0e-11_dp

  type, public :: amr_eb_patch_2d
    integer :: coarse_i_lower = 0
    integer :: coarse_i_upper = -1
    integer :: coarse_j_lower = 0
    integer :: coarse_j_upper = -1
    integer :: refinement_ratio = 0
  contains
    procedure :: is_valid => amr_eb_patch_is_valid
    procedure :: coarse_cell_count_x => amr_eb_patch_coarse_cell_count_x
    procedure :: coarse_cell_count_y => amr_eb_patch_coarse_cell_count_y
  end type amr_eb_patch_2d

  public :: build_amr_eb_patch_2d
  public :: average_down_eb_state_patch_2d
  public :: average_down_reactive_eb_state_patch_2d
  public :: composite_eb_integral_2d

contains

  pure integer function amr_eb_patch_coarse_cell_count_x(self) result(count)
    class(amr_eb_patch_2d), intent(in) :: self

    count = max(0, self%coarse_i_upper - self%coarse_i_lower + 1)
  end function amr_eb_patch_coarse_cell_count_x

  pure integer function amr_eb_patch_coarse_cell_count_y(self) result(count)
    class(amr_eb_patch_2d), intent(in) :: self

    count = max(0, self%coarse_j_upper - self%coarse_j_lower + 1)
  end function amr_eb_patch_coarse_cell_count_y

  pure logical function amr_eb_patch_is_valid( &
      self, coarse_geometry, fine_geometry) result(valid)
    class(amr_eb_patch_2d), intent(in) :: self
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry

    real(dp) :: expected_x_lower, expected_x_upper
    real(dp) :: expected_y_lower, expected_y_upper, scale, tolerance
    real(dp) :: restricted_volume_fraction
    integer :: coarse_i, coarse_j, fine_i_lower, fine_i_upper
    integer :: fine_j_lower, fine_j_upper, ratio

    valid = coarse_geometry%is_valid() .and. fine_geometry%is_valid()
    if (.not. valid) return
    ratio = self%refinement_ratio
    valid = ratio >= 2 .and. &
      self%coarse_i_lower >= 1 .and. &
      self%coarse_i_upper <= coarse_geometry%nx .and. &
      self%coarse_i_upper >= self%coarse_i_lower .and. &
      self%coarse_j_lower >= 1 .and. &
      self%coarse_j_upper <= coarse_geometry%ny .and. &
      self%coarse_j_upper >= self%coarse_j_lower
    if (.not. valid) return
    valid = fine_geometry%nx == self%coarse_cell_count_x() * ratio .and. &
      fine_geometry%ny == self%coarse_cell_count_y() * ratio
    if (.not. valid) return

    expected_x_lower = coarse_geometry%x_lower + &
      real(self%coarse_i_lower - 1, dp) * coarse_geometry%dx
    expected_x_upper = coarse_geometry%x_lower + &
      real(self%coarse_i_upper, dp) * coarse_geometry%dx
    expected_y_lower = coarse_geometry%y_lower + &
      real(self%coarse_j_lower - 1, dp) * coarse_geometry%dy
    expected_y_upper = coarse_geometry%y_lower + &
      real(self%coarse_j_upper, dp) * coarse_geometry%dy
    scale = max(1.0_dp, abs(expected_x_lower), abs(expected_x_upper), &
      abs(expected_y_lower), abs(expected_y_upper))
    tolerance = 1024.0_dp * epsilon(1.0_dp) * scale
    valid = abs(fine_geometry%x_lower - expected_x_lower) <= tolerance .and. &
      abs(fine_geometry%x_upper - expected_x_upper) <= tolerance .and. &
      abs(fine_geometry%y_lower - expected_y_lower) <= tolerance .and. &
      abs(fine_geometry%y_upper - expected_y_upper) <= tolerance .and. &
      abs(real(ratio, dp) * fine_geometry%dx - coarse_geometry%dx) <= &
        tolerance .and. &
      abs(real(ratio, dp) * fine_geometry%dy - coarse_geometry%dy) <= &
        tolerance
    if (.not. valid) return

    ! AMReX EB average-down assumes level geometry has a compatible volume
    ! measure. Reject a hierarchy whose parent and child measures disagree.
    do coarse_j = self%coarse_j_lower, self%coarse_j_upper
      fine_j_lower = (coarse_j - self%coarse_j_lower) * ratio + 1
      fine_j_upper = fine_j_lower + ratio - 1
      do coarse_i = self%coarse_i_lower, self%coarse_i_upper
        fine_i_lower = (coarse_i - self%coarse_i_lower) * ratio + 1
        fine_i_upper = fine_i_lower + ratio - 1
        restricted_volume_fraction = sum(fine_geometry%volume_fraction( &
          fine_i_lower:fine_i_upper, fine_j_lower:fine_j_upper)) / &
          real(ratio * ratio, dp)
        if (abs(restricted_volume_fraction - &
            coarse_geometry%volume_fraction(coarse_i, coarse_j)) > &
            geometry_consistency_tolerance) then
          valid = .false.
          return
        end if
      end do
    end do
  end function amr_eb_patch_is_valid

  subroutine build_amr_eb_patch_2d( &
      coarse_geometry, fine_geometry, coarse_i_lower, coarse_i_upper, &
      coarse_j_lower, coarse_j_upper, refinement_ratio, patch, ok)
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    integer, intent(in) :: coarse_i_lower, coarse_i_upper
    integer, intent(in) :: coarse_j_lower, coarse_j_upper, refinement_ratio
    type(amr_eb_patch_2d), intent(out) :: patch
    logical, intent(out) :: ok

    patch = amr_eb_patch_2d( &
      coarse_i_lower=coarse_i_lower, coarse_i_upper=coarse_i_upper, &
      coarse_j_lower=coarse_j_lower, coarse_j_upper=coarse_j_upper, &
      refinement_ratio=refinement_ratio)
    ok = patch%is_valid(coarse_geometry, fine_geometry)
  end subroutine build_amr_eb_patch_2d

  subroutine average_down_eb_state_patch_2d( &
      coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
      averaged_state, ok)
    real(dp), intent(in) :: coarse_state(:, :, :), fine_state(:, :, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    real(dp), intent(out) :: averaged_state(:, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate(:, :, :)
    real(dp) :: fine_volume
    integer :: coarse_i, coarse_j, fine_i_lower, fine_i_upper
    integer :: fine_j_lower, fine_j_upper, component, ratio

    averaged_state = 0.0_dp
    ok = .false.
    if (size(coarse_state, 1) < 1 .or. &
        size(coarse_state, 2) /= coarse_geometry%nx .or. &
        size(coarse_state, 3) /= coarse_geometry%ny .or. &
        size(fine_state, 1) /= size(coarse_state, 1) .or. &
        size(fine_state, 2) /= fine_geometry%nx .or. &
        size(fine_state, 3) /= fine_geometry%ny .or. &
        any(shape(averaged_state) /= shape(coarse_state))) return
    averaged_state = coarse_state
    if (.not. patch%is_valid(coarse_geometry, fine_geometry) .or. &
        any(.not. ieee_is_finite(coarse_state)) .or. &
        any(.not. ieee_is_finite(fine_state))) return

    allocate(candidate, source=coarse_state)
    ratio = patch%refinement_ratio
    do coarse_j = patch%coarse_j_lower, patch%coarse_j_upper
      fine_j_lower = (coarse_j - patch%coarse_j_lower) * ratio + 1
      fine_j_upper = fine_j_lower + ratio - 1
      do coarse_i = patch%coarse_i_lower, patch%coarse_i_upper
        fine_i_lower = (coarse_i - patch%coarse_i_lower) * ratio + 1
        fine_i_upper = fine_i_lower + ratio - 1
        fine_volume = sum(fine_geometry%volume_fraction( &
          fine_i_lower:fine_i_upper, fine_j_lower:fine_j_upper))
        if (fine_volume > tiny(1.0_dp)) then
          do component = 1, size(coarse_state, 1)
            candidate(component, coarse_i, coarse_j) = sum( &
              fine_geometry%volume_fraction( &
                fine_i_lower:fine_i_upper, fine_j_lower:fine_j_upper) * &
              fine_state(component, fine_i_lower:fine_i_upper, &
                fine_j_lower:fine_j_upper)) / fine_volume
          end do
        else
          candidate(:, coarse_i, coarse_j) = &
            fine_state(:, fine_i_lower, fine_j_lower)
        end if
      end do
    end do
    if (any(.not. ieee_is_finite(candidate))) return
    averaged_state = candidate
    ok = .true.
  end subroutine average_down_eb_state_patch_2d

  subroutine average_down_reactive_eb_state_patch_2d( &
      species, coarse_state, coarse_temperature, coarse_geometry, &
      fine_state, fine_geometry, patch, averaged_state, &
      averaged_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: coarse_state(:, :, :), coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    real(dp), intent(in) :: fine_state(:, :, :)
    type(amr_eb_patch_2d), intent(in) :: patch
    real(dp), intent(out) :: averaged_state(:, :, :)
    real(dp), intent(out) :: averaged_temperature(:, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate_state(:, :, :)
    real(dp), allocatable :: candidate_temperature(:, :), primitive(:)
    real(dp) :: recovered_temperature, sound_speed
    logical :: local_ok
    integer :: coarse_i, coarse_j, nvar

    averaged_state = 0.0_dp
    averaged_temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar < 1 .or. size(coarse_state, 1) /= nvar .or. &
        any(shape(averaged_state) /= shape(coarse_state)) .or. &
        any(shape(coarse_temperature) /= &
          [coarse_geometry%nx, coarse_geometry%ny]) .or. &
        any(shape(averaged_temperature) /= shape(coarse_temperature))) return
    averaged_state = coarse_state
    averaged_temperature = coarse_temperature
    if (any(.not. ieee_is_finite(coarse_temperature))) return

    allocate(candidate_state, mold=coarse_state)
    call average_down_eb_state_patch_2d( &
      coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
      candidate_state, local_ok)
    if (.not. local_ok) return
    allocate(candidate_temperature, source=coarse_temperature)
    allocate(primitive(reactive_nprim(size(species))))
    do coarse_j = patch%coarse_j_lower, patch%coarse_j_upper
      do coarse_i = patch%coarse_i_lower, patch%coarse_i_upper
        if (coarse_geometry%cell_type(coarse_i, coarse_j) == &
            eb_covered_cell) then
          candidate_state(:, coarse_i, coarse_j) = &
            coarse_state(:, coarse_i, coarse_j)
          cycle
        end if
        if (coarse_temperature(coarse_i, coarse_j) <= 0.0_dp) return
        call reactive_conserved_to_primitive( &
          species, candidate_state(:, coarse_i, coarse_j), &
          coarse_temperature(coarse_i, coarse_j), primitive, &
          recovered_temperature, sound_speed, local_ok)
        if (.not. local_ok) return
        candidate_temperature(coarse_i, coarse_j) = recovered_temperature
      end do
    end do
    averaged_state = candidate_state
    averaged_temperature = candidate_temperature
    ok = .true.
  end subroutine average_down_reactive_eb_state_patch_2d

  subroutine composite_eb_integral_2d( &
      coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
      integral, ok)
    real(dp), intent(in) :: coarse_state(:, :, :), fine_state(:, :, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    real(dp), intent(out) :: integral(:)
    logical, intent(out) :: ok

    integer :: i, j

    integral = 0.0_dp
    ok = .false.
    if (size(integral) /= size(coarse_state, 1) .or. &
        size(coarse_state, 2) /= coarse_geometry%nx .or. &
        size(coarse_state, 3) /= coarse_geometry%ny .or. &
        size(fine_state, 1) /= size(coarse_state, 1) .or. &
        size(fine_state, 2) /= fine_geometry%nx .or. &
        size(fine_state, 3) /= fine_geometry%ny .or. &
        .not. patch%is_valid(coarse_geometry, fine_geometry) .or. &
        any(.not. ieee_is_finite(coarse_state)) .or. &
        any(.not. ieee_is_finite(fine_state))) return

    do j = 1, coarse_geometry%ny
      do i = 1, coarse_geometry%nx
        if (i >= patch%coarse_i_lower .and. &
            i <= patch%coarse_i_upper .and. &
            j >= patch%coarse_j_lower .and. &
            j <= patch%coarse_j_upper) cycle
        integral = integral + coarse_geometry%volume_fraction(i, j) * &
          coarse_state(:, i, j) * coarse_geometry%dx * coarse_geometry%dy
      end do
    end do
    do j = 1, fine_geometry%ny
      do i = 1, fine_geometry%nx
        integral = integral + fine_geometry%volume_fraction(i, j) * &
          fine_state(:, i, j) * fine_geometry%dx * fine_geometry%dy
      end do
    end do
    if (any(.not. ieee_is_finite(integral))) then
      integral = 0.0_dp
      return
    end if
    ok = .true.
  end subroutine composite_eb_integral_2d

end module amr_eb_hierarchy_2d_mod
