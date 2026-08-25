module amr_eb_multilevel_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use eb_geometry_2d_mod, only: eb_geometry_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, average_down_eb_state_patch_2d, &
    average_down_reactive_eb_state_patch_2d
  implicit none
  private

  public :: average_down_three_level_eb_state_2d
  public :: average_down_three_level_reactive_eb_state_2d
  public :: composite_three_level_eb_integral_2d

contains

  subroutine average_down_three_level_eb_state_2d( &
      root_state, root_geometry, level_one_state, level_one_geometry, &
      root_patch, level_two_state, level_two_geometry, level_one_patch, &
      synchronized_root_state, synchronized_level_one_state, ok)
    real(dp), intent(in) :: root_state(:, :, :)
    real(dp), intent(in) :: level_one_state(:, :, :)
    real(dp), intent(in) :: level_two_state(:, :, :)
    type(eb_geometry_2d), intent(in) :: root_geometry
    type(eb_geometry_2d), intent(in) :: level_one_geometry
    type(eb_geometry_2d), intent(in) :: level_two_geometry
    type(amr_eb_patch_2d), intent(in) :: root_patch, level_one_patch
    real(dp), intent(out) :: synchronized_root_state(:, :, :)
    real(dp), intent(out) :: synchronized_level_one_state(:, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: root_candidate(:, :, :)
    real(dp), allocatable :: level_one_candidate(:, :, :)
    logical :: local_ok
    integer :: nvar

    synchronized_root_state = 0.0_dp
    synchronized_level_one_state = 0.0_dp
    ok = .false.
    nvar = size(root_state, 1)
    if (nvar < 1 .or. &
        any(shape(root_state) /= &
          [nvar, root_geometry%nx, root_geometry%ny]) .or. &
        any(shape(level_one_state) /= &
          [nvar, level_one_geometry%nx, level_one_geometry%ny]) .or. &
        any(shape(level_two_state) /= &
          [nvar, level_two_geometry%nx, level_two_geometry%ny]) .or. &
        any(shape(synchronized_root_state) /= shape(root_state)) .or. &
        any(shape(synchronized_level_one_state) /= &
          shape(level_one_state))) return
    synchronized_root_state = root_state
    synchronized_level_one_state = level_one_state
    if (.not. root_patch%is_valid(root_geometry, level_one_geometry) .or. &
        .not. level_one_patch%is_valid( &
          level_one_geometry, level_two_geometry) .or. &
        any(.not. ieee_is_finite(root_state)) .or. &
        any(.not. ieee_is_finite(level_one_state)) .or. &
        any(.not. ieee_is_finite(level_two_state))) return

    allocate(level_one_candidate, mold=level_one_state)
    call average_down_eb_state_patch_2d( &
      level_one_state, level_one_geometry, level_two_state, &
      level_two_geometry, level_one_patch, level_one_candidate, local_ok)
    if (.not. local_ok) return
    allocate(root_candidate, mold=root_state)
    call average_down_eb_state_patch_2d( &
      root_state, root_geometry, level_one_candidate, level_one_geometry, &
      root_patch, root_candidate, local_ok)
    if (.not. local_ok) return
    if (any(.not. ieee_is_finite(root_candidate)) .or. &
        any(.not. ieee_is_finite(level_one_candidate))) return
    synchronized_root_state = root_candidate
    synchronized_level_one_state = level_one_candidate
    ok = .true.
  end subroutine average_down_three_level_eb_state_2d

  subroutine average_down_three_level_reactive_eb_state_2d( &
      species, root_state, root_temperature, root_geometry, &
      level_one_state, level_one_temperature, level_one_geometry, &
      root_patch, level_two_state, level_two_temperature, &
      level_two_geometry, level_one_patch, synchronized_root_state, &
      synchronized_root_temperature, synchronized_level_one_state, &
      synchronized_level_one_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: root_state(:, :, :), root_temperature(:, :)
    real(dp), intent(in) :: level_one_state(:, :, :)
    real(dp), intent(in) :: level_one_temperature(:, :)
    real(dp), intent(in) :: level_two_state(:, :, :)
    real(dp), intent(in) :: level_two_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: root_geometry
    type(eb_geometry_2d), intent(in) :: level_one_geometry
    type(eb_geometry_2d), intent(in) :: level_two_geometry
    type(amr_eb_patch_2d), intent(in) :: root_patch, level_one_patch
    real(dp), intent(out) :: synchronized_root_state(:, :, :)
    real(dp), intent(out) :: synchronized_root_temperature(:, :)
    real(dp), intent(out) :: synchronized_level_one_state(:, :, :)
    real(dp), intent(out) :: synchronized_level_one_temperature(:, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: root_candidate(:, :, :)
    real(dp), allocatable :: root_temperature_candidate(:, :)
    real(dp), allocatable :: level_one_candidate(:, :, :)
    real(dp), allocatable :: level_one_temperature_candidate(:, :)
    logical :: local_ok
    integer :: nvar

    synchronized_root_state = 0.0_dp
    synchronized_root_temperature = 0.0_dp
    synchronized_level_one_state = 0.0_dp
    synchronized_level_one_temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar < 1 .or. &
        any(shape(root_state) /= &
          [nvar, root_geometry%nx, root_geometry%ny]) .or. &
        any(shape(root_temperature) /= &
          [root_geometry%nx, root_geometry%ny]) .or. &
        any(shape(level_one_state) /= &
          [nvar, level_one_geometry%nx, level_one_geometry%ny]) .or. &
        any(shape(level_one_temperature) /= &
          [level_one_geometry%nx, level_one_geometry%ny]) .or. &
        any(shape(level_two_state) /= &
          [nvar, level_two_geometry%nx, level_two_geometry%ny]) .or. &
        any(shape(level_two_temperature) /= &
          [level_two_geometry%nx, level_two_geometry%ny]) .or. &
        any(shape(synchronized_root_state) /= shape(root_state)) .or. &
        any(shape(synchronized_root_temperature) /= &
          shape(root_temperature)) .or. &
        any(shape(synchronized_level_one_state) /= &
          shape(level_one_state)) .or. &
        any(shape(synchronized_level_one_temperature) /= &
          shape(level_one_temperature))) return
    synchronized_root_state = root_state
    synchronized_root_temperature = root_temperature
    synchronized_level_one_state = level_one_state
    synchronized_level_one_temperature = level_one_temperature
    if (.not. root_patch%is_valid(root_geometry, level_one_geometry) .or. &
        .not. level_one_patch%is_valid( &
          level_one_geometry, level_two_geometry) .or. &
        any(.not. ieee_is_finite(root_state)) .or. &
        any(.not. ieee_is_finite(root_temperature)) .or. &
        any(.not. ieee_is_finite(level_one_state)) .or. &
        any(.not. ieee_is_finite(level_one_temperature)) .or. &
        any(.not. ieee_is_finite(level_two_state)) .or. &
        any(.not. ieee_is_finite(level_two_temperature)) .or. &
        any(root_temperature <= 0.0_dp) .or. &
        any(level_one_temperature <= 0.0_dp) .or. &
        any(level_two_temperature <= 0.0_dp)) return

    allocate(level_one_candidate, mold=level_one_state)
    allocate(level_one_temperature_candidate, mold=level_one_temperature)
    call average_down_reactive_eb_state_patch_2d( &
      species, level_one_state, level_one_temperature, level_one_geometry, &
      level_two_state, level_two_geometry, level_one_patch, &
      level_one_candidate, level_one_temperature_candidate, local_ok)
    if (.not. local_ok) return
    allocate(root_candidate, mold=root_state)
    allocate(root_temperature_candidate, mold=root_temperature)
    call average_down_reactive_eb_state_patch_2d( &
      species, root_state, root_temperature, root_geometry, &
      level_one_candidate, level_one_geometry, root_patch, root_candidate, &
      root_temperature_candidate, local_ok)
    if (.not. local_ok) return
    if (any(.not. ieee_is_finite(root_candidate)) .or. &
        any(.not. ieee_is_finite(root_temperature_candidate)) .or. &
        any(.not. ieee_is_finite(level_one_candidate)) .or. &
        any(.not. ieee_is_finite(level_one_temperature_candidate)) .or. &
        any(root_temperature_candidate <= 0.0_dp) .or. &
        any(level_one_temperature_candidate <= 0.0_dp)) return
    synchronized_root_state = root_candidate
    synchronized_root_temperature = root_temperature_candidate
    synchronized_level_one_state = level_one_candidate
    synchronized_level_one_temperature = level_one_temperature_candidate
    ok = .true.
  end subroutine average_down_three_level_reactive_eb_state_2d

  subroutine composite_three_level_eb_integral_2d( &
      root_state, root_geometry, level_one_state, level_one_geometry, &
      root_patch, level_two_state, level_two_geometry, level_one_patch, &
      integral, ok)
    real(dp), intent(in) :: root_state(:, :, :)
    real(dp), intent(in) :: level_one_state(:, :, :)
    real(dp), intent(in) :: level_two_state(:, :, :)
    type(eb_geometry_2d), intent(in) :: root_geometry
    type(eb_geometry_2d), intent(in) :: level_one_geometry
    type(eb_geometry_2d), intent(in) :: level_two_geometry
    type(amr_eb_patch_2d), intent(in) :: root_patch, level_one_patch
    real(dp), intent(out) :: integral(:)
    logical, intent(out) :: ok

    integer :: i, j, nvar

    integral = 0.0_dp
    ok = .false.
    nvar = size(root_state, 1)
    if (nvar < 1 .or. size(integral) /= nvar .or. &
        any(shape(root_state) /= &
          [nvar, root_geometry%nx, root_geometry%ny]) .or. &
        any(shape(level_one_state) /= &
          [nvar, level_one_geometry%nx, level_one_geometry%ny]) .or. &
        any(shape(level_two_state) /= &
          [nvar, level_two_geometry%nx, level_two_geometry%ny]) .or. &
        .not. root_patch%is_valid(root_geometry, level_one_geometry) .or. &
        .not. level_one_patch%is_valid( &
          level_one_geometry, level_two_geometry) .or. &
        any(.not. ieee_is_finite(root_state)) .or. &
        any(.not. ieee_is_finite(level_one_state)) .or. &
        any(.not. ieee_is_finite(level_two_state))) return

    do j = 1, root_geometry%ny
      do i = 1, root_geometry%nx
        if (cell_is_in_patch(root_patch, i, j)) cycle
        integral = integral + root_geometry%volume_fraction(i, j) * &
          root_state(:, i, j) * root_geometry%dx * root_geometry%dy
      end do
    end do
    do j = 1, level_one_geometry%ny
      do i = 1, level_one_geometry%nx
        if (cell_is_in_patch(level_one_patch, i, j)) cycle
        integral = integral + level_one_geometry%volume_fraction(i, j) * &
          level_one_state(:, i, j) * level_one_geometry%dx * &
          level_one_geometry%dy
      end do
    end do
    do j = 1, level_two_geometry%ny
      do i = 1, level_two_geometry%nx
        integral = integral + level_two_geometry%volume_fraction(i, j) * &
          level_two_state(:, i, j) * level_two_geometry%dx * &
          level_two_geometry%dy
      end do
    end do
    if (any(.not. ieee_is_finite(integral))) then
      integral = 0.0_dp
      return
    end if
    ok = .true.
  end subroutine composite_three_level_eb_integral_2d

  pure logical function cell_is_in_patch(patch, i, j) result(inside)
    type(amr_eb_patch_2d), intent(in) :: patch
    integer, intent(in) :: i, j

    inside = i >= patch%coarse_i_lower .and. &
      i <= patch%coarse_i_upper .and. &
      j >= patch%coarse_j_lower .and. j <= patch%coarse_j_upper
  end function cell_is_in_patch

end module amr_eb_multilevel_2d_mod
