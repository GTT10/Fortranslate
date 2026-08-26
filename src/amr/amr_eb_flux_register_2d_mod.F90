module amr_eb_flux_register_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use eb_geometry_2d_mod, only: &
    eb_geometry_2d, eb_covered_cell, eb_cut_cell, eb_regular_cell
  use amr_eb_hierarchy_2d_mod, only: amr_eb_patch_2d
  implicit none
  private

  type, public :: amr_eb_flux_register_2d
    integer :: component_count = 0
    integer :: coarse_nx = 0
    integer :: coarse_ny = 0
    integer :: correction_i_lower = 0
    integer :: correction_i_upper = -1
    integer :: correction_j_lower = 0
    integer :: correction_j_upper = -1
    type(amr_eb_patch_2d) :: patch
    real(dp), allocatable :: correction(:, :, :)
  contains
    procedure :: is_valid => amr_eb_flux_register_is_valid
    procedure :: reset => reset_amr_eb_flux_register_2d
  end type amr_eb_flux_register_2d

  public :: initialize_amr_eb_flux_register_2d
  public :: accumulate_coarse_eb_fluxes_2d
  public :: accumulate_fine_eb_fluxes_2d
  public :: reflux_eb_state_patch_support_2d
  public :: reflux_eb_state_patch_2d
  public :: reflux_reactive_eb_state_patch_support_2d
  public :: reflux_reactive_eb_state_patch_2d

contains

  pure logical function patches_match(first, second) result(matches)
    type(amr_eb_patch_2d), intent(in) :: first, second

    matches = first%coarse_i_lower == second%coarse_i_lower .and. &
      first%coarse_i_upper == second%coarse_i_upper .and. &
      first%coarse_j_lower == second%coarse_j_lower .and. &
      first%coarse_j_upper == second%coarse_j_upper .and. &
      first%refinement_ratio == second%refinement_ratio
  end function patches_match

  pure logical function amr_eb_flux_register_is_valid( &
      self, coarse_geometry, fine_geometry, patch) result(valid)
    class(amr_eb_flux_register_2d), intent(in) :: self
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch

    valid = self%component_count >= 1 .and. &
      self%coarse_nx == coarse_geometry%nx .and. &
      self%coarse_ny == coarse_geometry%ny .and. &
      self%correction_i_lower == max(1, patch%coarse_i_lower - 1) .and. &
      self%correction_i_upper == &
        min(coarse_geometry%nx, patch%coarse_i_upper + 1) .and. &
      self%correction_j_lower == max(1, patch%coarse_j_lower - 1) .and. &
      self%correction_j_upper == &
        min(coarse_geometry%ny, patch%coarse_j_upper + 1) .and. &
      patches_match(self%patch, patch) .and. &
      patch%is_valid(coarse_geometry, fine_geometry) .and. &
      allocated(self%correction)
    if (.not. valid) return
    valid = all(shape(self%correction) == &
      [self%component_count, &
       self%correction_i_upper - self%correction_i_lower + 1, &
       self%correction_j_upper - self%correction_j_lower + 1]) .and. &
      lbound(self%correction, 2) == self%correction_i_lower .and. &
      ubound(self%correction, 2) == self%correction_i_upper .and. &
      lbound(self%correction, 3) == self%correction_j_lower .and. &
      ubound(self%correction, 3) == self%correction_j_upper .and. &
      all(ieee_is_finite(self%correction))
  end function amr_eb_flux_register_is_valid

  subroutine initialize_amr_eb_flux_register_2d( &
      coarse_geometry, fine_geometry, patch, component_count, register, ok)
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    integer, intent(in) :: component_count
    type(amr_eb_flux_register_2d), intent(out) :: register
    logical, intent(out) :: ok

    ok = .false.
    if (component_count < 1 .or. &
        .not. patch%is_valid(coarse_geometry, fine_geometry)) return
    register%component_count = component_count
    register%coarse_nx = coarse_geometry%nx
    register%coarse_ny = coarse_geometry%ny
    register%correction_i_lower = max(1, patch%coarse_i_lower - 1)
    register%correction_i_upper = &
      min(coarse_geometry%nx, patch%coarse_i_upper + 1)
    register%correction_j_lower = max(1, patch%coarse_j_lower - 1)
    register%correction_j_upper = &
      min(coarse_geometry%ny, patch%coarse_j_upper + 1)
    register%patch = patch
    allocate(register%correction( &
      component_count, &
      register%correction_i_lower:register%correction_i_upper, &
      register%correction_j_lower:register%correction_j_upper))
    register%correction = 0.0_dp
    ok = register%is_valid(coarse_geometry, fine_geometry, patch)
  end subroutine initialize_amr_eb_flux_register_2d

  subroutine reset_amr_eb_flux_register_2d(self)
    class(amr_eb_flux_register_2d), intent(inout) :: self

    if (allocated(self%correction)) self%correction = 0.0_dp
  end subroutine reset_amr_eb_flux_register_2d

  subroutine accumulate_coarse_eb_fluxes_2d( &
      register, coarse_geometry, fine_geometry, patch, x_flux, y_flux, dt, ok)
    type(amr_eb_flux_register_2d), intent(inout) :: register
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    real(dp), intent(in) :: x_flux(:, 0:, :), y_flux(:, :, 0:), dt
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate(:, :, :)
    real(dp) :: kappa, scale
    integer :: i, j

    ok = .false.
    if (.not. register%is_valid(coarse_geometry, fine_geometry, patch) .or. &
        .not. ieee_is_finite(dt) .or. dt < 0.0_dp .or. &
        size(x_flux, 1) /= register%component_count .or. &
        size(x_flux, 2) /= coarse_geometry%nx + 1 .or. &
        size(x_flux, 3) /= coarse_geometry%ny .or. &
        size(y_flux, 1) /= register%component_count .or. &
        size(y_flux, 2) /= coarse_geometry%nx .or. &
        size(y_flux, 3) /= coarse_geometry%ny + 1 .or. &
        any(.not. ieee_is_finite(x_flux)) .or. &
        any(.not. ieee_is_finite(y_flux))) return
    allocate(candidate( &
      register%component_count, &
      register%correction_i_lower:register%correction_i_upper, &
      register%correction_j_lower:register%correction_j_upper))
    candidate = register%correction

    if (patch%coarse_i_lower > 1) then
      i = patch%coarse_i_lower - 1
      do j = patch%coarse_j_lower, patch%coarse_j_upper
        kappa = coarse_geometry%volume_fraction(i, j)
        if (kappa <= tiny(1.0_dp)) cycle
        scale = dt * coarse_geometry%x_face_fraction(i, j) / &
          (kappa * coarse_geometry%dx)
        candidate(:, i, j) = candidate(:, i, j) + scale * x_flux(:, i, j)
      end do
    end if
    if (patch%coarse_i_upper < coarse_geometry%nx) then
      i = patch%coarse_i_upper + 1
      do j = patch%coarse_j_lower, patch%coarse_j_upper
        kappa = coarse_geometry%volume_fraction(i, j)
        if (kappa <= tiny(1.0_dp)) cycle
        scale = dt * coarse_geometry%x_face_fraction(i - 1, j) / &
          (kappa * coarse_geometry%dx)
        candidate(:, i, j) = candidate(:, i, j) - scale * &
          x_flux(:, i - 1, j)
      end do
    end if
    if (patch%coarse_j_lower > 1) then
      j = patch%coarse_j_lower - 1
      do i = patch%coarse_i_lower, patch%coarse_i_upper
        kappa = coarse_geometry%volume_fraction(i, j)
        if (kappa <= tiny(1.0_dp)) cycle
        scale = dt * coarse_geometry%y_face_fraction(i, j) / &
          (kappa * coarse_geometry%dy)
        candidate(:, i, j) = candidate(:, i, j) + scale * y_flux(:, i, j)
      end do
    end if
    if (patch%coarse_j_upper < coarse_geometry%ny) then
      j = patch%coarse_j_upper + 1
      do i = patch%coarse_i_lower, patch%coarse_i_upper
        kappa = coarse_geometry%volume_fraction(i, j)
        if (kappa <= tiny(1.0_dp)) cycle
        scale = dt * coarse_geometry%y_face_fraction(i, j - 1) / &
          (kappa * coarse_geometry%dy)
        candidate(:, i, j) = candidate(:, i, j) - scale * &
          y_flux(:, i, j - 1)
      end do
    end if
    if (any(.not. ieee_is_finite(candidate))) return
    register%correction = candidate
    ok = .true.
  end subroutine accumulate_coarse_eb_fluxes_2d

  subroutine accumulate_fine_eb_fluxes_2d( &
      register, coarse_geometry, fine_geometry, patch, x_flux, y_flux, dt, ok)
    type(amr_eb_flux_register_2d), intent(inout) :: register
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    real(dp), intent(in) :: x_flux(:, 0:, :), y_flux(:, :, 0:), dt
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate(:, :, :), integrated_flux(:)
    real(dp) :: kappa, scale
    integer :: coarse_i, coarse_j, fine_i, fine_j, lower, upper, ratio

    ok = .false.
    if (.not. register%is_valid(coarse_geometry, fine_geometry, patch) .or. &
        .not. ieee_is_finite(dt) .or. dt < 0.0_dp .or. &
        size(x_flux, 1) /= register%component_count .or. &
        size(x_flux, 2) /= fine_geometry%nx + 1 .or. &
        size(x_flux, 3) /= fine_geometry%ny .or. &
        size(y_flux, 1) /= register%component_count .or. &
        size(y_flux, 2) /= fine_geometry%nx .or. &
        size(y_flux, 3) /= fine_geometry%ny + 1 .or. &
        any(.not. ieee_is_finite(x_flux)) .or. &
        any(.not. ieee_is_finite(y_flux))) return
    allocate(candidate( &
      register%component_count, &
      register%correction_i_lower:register%correction_i_upper, &
      register%correction_j_lower:register%correction_j_upper))
    candidate = register%correction
    allocate(integrated_flux(register%component_count))
    ratio = patch%refinement_ratio

    if (patch%coarse_i_lower > 1) then
      coarse_i = patch%coarse_i_lower - 1
      fine_i = 0
      do coarse_j = patch%coarse_j_lower, patch%coarse_j_upper
        kappa = coarse_geometry%volume_fraction(coarse_i, coarse_j)
        if (kappa <= tiny(1.0_dp)) cycle
        lower = (coarse_j - patch%coarse_j_lower) * ratio + 1
        upper = lower + ratio - 1
        integrated_flux = 0.0_dp
        do fine_j = lower, upper
          integrated_flux = integrated_flux + &
            fine_geometry%x_face_fraction(fine_i, fine_j) * &
            x_flux(:, fine_i, fine_j) * fine_geometry%dy
        end do
        scale = dt / (kappa * coarse_geometry%dx * coarse_geometry%dy)
        candidate(:, coarse_i, coarse_j) = &
          candidate(:, coarse_i, coarse_j) - scale * integrated_flux
      end do
    end if
    if (patch%coarse_i_upper < coarse_geometry%nx) then
      coarse_i = patch%coarse_i_upper + 1
      fine_i = fine_geometry%nx
      do coarse_j = patch%coarse_j_lower, patch%coarse_j_upper
        kappa = coarse_geometry%volume_fraction(coarse_i, coarse_j)
        if (kappa <= tiny(1.0_dp)) cycle
        lower = (coarse_j - patch%coarse_j_lower) * ratio + 1
        upper = lower + ratio - 1
        integrated_flux = 0.0_dp
        do fine_j = lower, upper
          integrated_flux = integrated_flux + &
            fine_geometry%x_face_fraction(fine_i, fine_j) * &
            x_flux(:, fine_i, fine_j) * fine_geometry%dy
        end do
        scale = dt / (kappa * coarse_geometry%dx * coarse_geometry%dy)
        candidate(:, coarse_i, coarse_j) = &
          candidate(:, coarse_i, coarse_j) + scale * integrated_flux
      end do
    end if
    if (patch%coarse_j_lower > 1) then
      coarse_j = patch%coarse_j_lower - 1
      fine_j = 0
      do coarse_i = patch%coarse_i_lower, patch%coarse_i_upper
        kappa = coarse_geometry%volume_fraction(coarse_i, coarse_j)
        if (kappa <= tiny(1.0_dp)) cycle
        lower = (coarse_i - patch%coarse_i_lower) * ratio + 1
        upper = lower + ratio - 1
        integrated_flux = 0.0_dp
        do fine_i = lower, upper
          integrated_flux = integrated_flux + &
            fine_geometry%y_face_fraction(fine_i, fine_j) * &
            y_flux(:, fine_i, fine_j) * fine_geometry%dx
        end do
        scale = dt / (kappa * coarse_geometry%dx * coarse_geometry%dy)
        candidate(:, coarse_i, coarse_j) = &
          candidate(:, coarse_i, coarse_j) - scale * integrated_flux
      end do
    end if
    if (patch%coarse_j_upper < coarse_geometry%ny) then
      coarse_j = patch%coarse_j_upper + 1
      fine_j = fine_geometry%ny
      do coarse_i = patch%coarse_i_lower, patch%coarse_i_upper
        kappa = coarse_geometry%volume_fraction(coarse_i, coarse_j)
        if (kappa <= tiny(1.0_dp)) cycle
        lower = (coarse_i - patch%coarse_i_lower) * ratio + 1
        upper = lower + ratio - 1
        integrated_flux = 0.0_dp
        do fine_i = lower, upper
          integrated_flux = integrated_flux + &
            fine_geometry%y_face_fraction(fine_i, fine_j) * &
            y_flux(:, fine_i, fine_j) * fine_geometry%dx
        end do
        scale = dt / (kappa * coarse_geometry%dx * coarse_geometry%dy)
        candidate(:, coarse_i, coarse_j) = &
          candidate(:, coarse_i, coarse_j) + scale * integrated_flux
      end do
    end if
    if (any(.not. ieee_is_finite(candidate))) return
    register%correction = candidate
    ok = .true.
  end subroutine accumulate_fine_eb_fluxes_2d

  pure logical function cell_is_inside_patch(patch, i, j) result(inside)
    type(amr_eb_patch_2d), intent(in) :: patch
    integer, intent(in) :: i, j

    inside = i >= patch%coarse_i_lower .and. &
      i <= patch%coarse_i_upper .and. &
      j >= patch%coarse_j_lower .and. j <= patch%coarse_j_upper
  end function cell_is_inside_patch

  pure logical function cardinal_cells_connected( &
      geometry, first_i, first_j, second_i, second_j) result(connected)
    type(eb_geometry_2d), intent(in) :: geometry
    integer, intent(in) :: first_i, first_j, second_i, second_j

    connected = .false.
    if (first_i < 1 .or. first_i > geometry%nx .or. &
        second_i < 1 .or. second_i > geometry%nx .or. &
        first_j < 1 .or. first_j > geometry%ny .or. &
        second_j < 1 .or. second_j > geometry%ny) return
    if (geometry%cell_type(first_i, first_j) == eb_covered_cell .or. &
        geometry%cell_type(second_i, second_j) == eb_covered_cell) return
    if (second_i == first_i + 1 .and. second_j == first_j) then
      connected = geometry%x_face_fraction(first_i, first_j) > 0.0_dp
    else if (second_i == first_i - 1 .and. second_j == first_j) then
      connected = geometry%x_face_fraction(second_i, first_j) > 0.0_dp
    else if (second_j == first_j + 1 .and. second_i == first_i) then
      connected = geometry%y_face_fraction(first_i, first_j) > 0.0_dp
    else if (second_j == first_j - 1 .and. second_i == first_i) then
      connected = geometry%y_face_fraction(first_i, second_j) > 0.0_dp
    end if
  end function cardinal_cells_connected

  pure logical function cells_connected( &
      geometry, i, j, offset_i, offset_j) result(connected)
    type(eb_geometry_2d), intent(in) :: geometry
    integer, intent(in) :: i, j, offset_i, offset_j
    integer :: neighbor_i, neighbor_j
    logical :: horizontal_first, vertical_first

    connected = .false.
    if (abs(offset_i) > 1 .or. abs(offset_j) > 1 .or. &
        (offset_i == 0 .and. offset_j == 0)) return
    neighbor_i = i + offset_i
    neighbor_j = j + offset_j
    if (neighbor_i < 1 .or. neighbor_i > geometry%nx .or. &
        neighbor_j < 1 .or. neighbor_j > geometry%ny) return
    if (offset_i == 0 .or. offset_j == 0) then
      connected = cardinal_cells_connected( &
        geometry, i, j, neighbor_i, neighbor_j)
    else
      horizontal_first = cardinal_cells_connected( &
        geometry, i, j, neighbor_i, j)
      if (horizontal_first) horizontal_first = cardinal_cells_connected( &
        geometry, neighbor_i, j, neighbor_i, neighbor_j)
      vertical_first = cardinal_cells_connected( &
        geometry, i, j, i, neighbor_j)
      if (vertical_first) vertical_first = cardinal_cells_connected( &
        geometry, i, neighbor_j, neighbor_i, neighbor_j)
      connected = horizontal_first .or. vertical_first
    end if
  end function cells_connected

  subroutine reflux_eb_state_patch_support_2d( &
      coarse_i_lower, coarse_j_lower, coarse_state, coarse_geometry, &
      fine_state, fine_geometry, patch, register, refluxed_coarse_state, &
      refluxed_fine_state, ok)
    integer, intent(in) :: coarse_i_lower, coarse_j_lower
    real(dp), intent(in) :: coarse_state(:, coarse_i_lower:, coarse_j_lower:)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    real(dp), intent(in) :: fine_state(:, :, :)
    type(eb_geometry_2d), intent(in) :: fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    type(amr_eb_flux_register_2d), intent(inout) :: register
    real(dp), intent(out) :: &
      refluxed_coarse_state(:, coarse_i_lower:, coarse_j_lower:)
    real(dp), intent(out) :: refluxed_fine_state(:, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: coarse_increment(:, :, :), fine_increment(:, :, :)
    real(dp), allocatable :: correction(:), extensive_correction(:)
    real(dp), allocatable :: recipient_increment(:)
    real(dp) :: kappa, neighbor_volume
    integer :: coarse_i_upper, coarse_j_upper
    integer :: expected_i_lower, expected_i_upper
    integer :: expected_j_lower, expected_j_upper
    integer :: i, j, offset_i, offset_j, neighbor_i, neighbor_j
    integer :: fine_i_lower, fine_i_upper, fine_j_lower, fine_j_upper, ratio

    refluxed_coarse_state = 0.0_dp
    refluxed_fine_state = 0.0_dp
    ok = .false.
    coarse_i_upper = coarse_i_lower + size(coarse_state, 2) - 1
    coarse_j_upper = coarse_j_lower + size(coarse_state, 3) - 1
    expected_i_lower = max(1, patch%coarse_i_lower - 2)
    expected_i_upper = min(coarse_geometry%nx, patch%coarse_i_upper + 2)
    expected_j_lower = max(1, patch%coarse_j_lower - 2)
    expected_j_upper = min(coarse_geometry%ny, patch%coarse_j_upper + 2)
    if (size(coarse_state, 1) /= register%component_count .or. &
        coarse_i_lower < 1 .or. coarse_i_upper > coarse_geometry%nx .or. &
        coarse_j_lower < 1 .or. coarse_j_upper > coarse_geometry%ny .or. &
        coarse_i_lower > expected_i_lower .or. &
        coarse_i_upper < expected_i_upper .or. &
        coarse_j_lower > expected_j_lower .or. &
        coarse_j_upper < expected_j_upper .or. &
        size(fine_state, 1) /= register%component_count .or. &
        size(fine_state, 2) /= fine_geometry%nx .or. &
        size(fine_state, 3) /= fine_geometry%ny .or. &
        any(shape(refluxed_coarse_state) /= shape(coarse_state)) .or. &
        any(shape(refluxed_fine_state) /= shape(fine_state))) return
    refluxed_coarse_state = coarse_state
    refluxed_fine_state = fine_state
    if (.not. register%is_valid(coarse_geometry, fine_geometry, patch) .or. &
        any(.not. ieee_is_finite(coarse_state)) .or. &
        any(.not. ieee_is_finite(fine_state))) return

    allocate(coarse_increment( &
      register%component_count, coarse_i_lower:coarse_i_upper, &
      coarse_j_lower:coarse_j_upper))
    allocate(fine_increment, mold=fine_state)
    allocate(correction(register%component_count))
    allocate(extensive_correction(register%component_count))
    allocate(recipient_increment(register%component_count))
    coarse_increment = 0.0_dp
    fine_increment = 0.0_dp
    ratio = patch%refinement_ratio

    do j = register%correction_j_lower, register%correction_j_upper
      do i = register%correction_i_lower, register%correction_i_upper
        correction = register%correction(:, i, j)
        if (all(correction == 0.0_dp)) cycle
        select case (coarse_geometry%cell_type(i, j))
        case (eb_regular_cell)
          coarse_increment(:, i, j) = &
            coarse_increment(:, i, j) + correction
        case (eb_cut_cell)
          kappa = coarse_geometry%volume_fraction(i, j)
          extensive_correction = kappa * correction
          coarse_increment(:, i, j) = &
            coarse_increment(:, i, j) + extensive_correction
          neighbor_volume = 0.0_dp
          do offset_j = -1, 1
            do offset_i = -1, 1
              if (.not. cells_connected( &
                  coarse_geometry, i, j, offset_i, offset_j)) cycle
              neighbor_volume = neighbor_volume + &
                coarse_geometry%volume_fraction(i + offset_i, j + offset_j)
            end do
          end do
          if (neighbor_volume <= tiny(1.0_dp)) return
          recipient_increment = extensive_correction * &
            (1.0_dp - kappa) / neighbor_volume
          do offset_j = -1, 1
            do offset_i = -1, 1
              if (.not. cells_connected( &
                  coarse_geometry, i, j, offset_i, offset_j)) cycle
              neighbor_i = i + offset_i
              neighbor_j = j + offset_j
              if (cell_is_inside_patch(patch, neighbor_i, neighbor_j)) then
                fine_i_lower = &
                  (neighbor_i - patch%coarse_i_lower) * ratio + 1
                fine_i_upper = fine_i_lower + ratio - 1
                fine_j_lower = &
                  (neighbor_j - patch%coarse_j_lower) * ratio + 1
                fine_j_upper = fine_j_lower + ratio - 1
                fine_increment(:, fine_i_lower:fine_i_upper, &
                  fine_j_lower:fine_j_upper) = &
                  fine_increment(:, fine_i_lower:fine_i_upper, &
                    fine_j_lower:fine_j_upper) + &
                  spread(spread(recipient_increment, 2, ratio), 3, ratio)
              else
                coarse_increment(:, neighbor_i, neighbor_j) = &
                  coarse_increment(:, neighbor_i, neighbor_j) + &
                  recipient_increment
              end if
            end do
          end do
        case (eb_covered_cell)
          return
        case default
          return
        end select
      end do
    end do

    refluxed_coarse_state = coarse_state + coarse_increment
    refluxed_fine_state = fine_state + fine_increment
    if (any(.not. ieee_is_finite(refluxed_coarse_state)) .or. &
        any(.not. ieee_is_finite(refluxed_fine_state))) then
      refluxed_coarse_state = coarse_state
      refluxed_fine_state = fine_state
      return
    end if
    call register%reset()
    ok = .true.
  end subroutine reflux_eb_state_patch_support_2d

  subroutine reflux_eb_state_patch_2d( &
      coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
      register, refluxed_coarse_state, refluxed_fine_state, ok)
    real(dp), intent(in) :: coarse_state(:, :, :), fine_state(:, :, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    type(amr_eb_flux_register_2d), intent(inout) :: register
    real(dp), intent(out) :: refluxed_coarse_state(:, :, :)
    real(dp), intent(out) :: refluxed_fine_state(:, :, :)
    logical, intent(out) :: ok

    call reflux_eb_state_patch_support_2d( &
      1, 1, coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
      register, refluxed_coarse_state, refluxed_fine_state, ok)
  end subroutine reflux_eb_state_patch_2d

  subroutine reflux_reactive_eb_state_patch_support_2d( &
      species, coarse_i_lower, coarse_j_lower, coarse_state, &
      coarse_temperature, coarse_geometry, fine_state, fine_temperature, &
      fine_geometry, patch, register, refluxed_coarse_state, &
      refluxed_coarse_temperature, refluxed_fine_state, &
      refluxed_fine_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    integer, intent(in) :: coarse_i_lower, coarse_j_lower
    real(dp), intent(in) :: coarse_state(:, coarse_i_lower:, coarse_j_lower:)
    real(dp), intent(in) :: &
      coarse_temperature(coarse_i_lower:, coarse_j_lower:)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    real(dp), intent(in) :: fine_state(:, :, :), fine_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    type(amr_eb_flux_register_2d), intent(inout) :: register
    real(dp), intent(out) :: &
      refluxed_coarse_state(:, coarse_i_lower:, coarse_j_lower:)
    real(dp), intent(out) :: &
      refluxed_coarse_temperature(coarse_i_lower:, coarse_j_lower:)
    real(dp), intent(out) :: refluxed_fine_state(:, :, :)
    real(dp), intent(out) :: refluxed_fine_temperature(:, :)
    logical, intent(out) :: ok

    type(amr_eb_flux_register_2d) :: candidate_register
    real(dp), allocatable :: candidate_coarse_state(:, :, :)
    real(dp), allocatable :: candidate_fine_state(:, :, :)
    real(dp), allocatable :: candidate_coarse_temperature(:, :)
    real(dp), allocatable :: candidate_fine_temperature(:, :), primitive(:)
    real(dp) :: recovered_temperature, sound_speed
    logical :: local_ok
    integer :: coarse_i_upper, coarse_j_upper, i, j, nvar

    refluxed_coarse_state = 0.0_dp
    refluxed_coarse_temperature = 0.0_dp
    refluxed_fine_state = 0.0_dp
    refluxed_fine_temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    coarse_i_upper = coarse_i_lower + size(coarse_state, 2) - 1
    coarse_j_upper = coarse_j_lower + size(coarse_state, 3) - 1
    if (nvar < 1 .or. size(coarse_state, 1) /= nvar .or. &
        size(fine_state, 1) /= nvar .or. &
        any(shape(coarse_temperature) /= &
          [size(coarse_state, 2), size(coarse_state, 3)]) .or. &
        any(shape(fine_temperature) /= &
          [fine_geometry%nx, fine_geometry%ny]) .or. &
        any(shape(refluxed_coarse_state) /= shape(coarse_state)) .or. &
        any(shape(refluxed_coarse_temperature) /= &
          shape(coarse_temperature)) .or. &
        any(shape(refluxed_fine_state) /= shape(fine_state)) .or. &
        any(shape(refluxed_fine_temperature) /= &
          shape(fine_temperature))) return
    refluxed_coarse_state = coarse_state
    refluxed_coarse_temperature = coarse_temperature
    refluxed_fine_state = fine_state
    refluxed_fine_temperature = fine_temperature
    if (any(.not. ieee_is_finite(coarse_temperature)) .or. &
        any(.not. ieee_is_finite(fine_temperature))) return

    candidate_register = register
    allocate(candidate_coarse_state( &
      nvar, coarse_i_lower:coarse_i_upper, coarse_j_lower:coarse_j_upper))
    allocate(candidate_fine_state, mold=fine_state)
    call reflux_eb_state_patch_support_2d( &
      coarse_i_lower, coarse_j_lower, coarse_state, coarse_geometry, &
      fine_state, fine_geometry, patch, candidate_register, &
      candidate_coarse_state, candidate_fine_state, local_ok)
    if (.not. local_ok) return
    allocate(candidate_coarse_temperature( &
      coarse_i_lower:coarse_i_upper, coarse_j_lower:coarse_j_upper))
    candidate_coarse_temperature = coarse_temperature
    allocate(candidate_fine_temperature, source=fine_temperature)
    allocate(primitive(reactive_nprim(size(species))))

    do j = coarse_j_lower, coarse_j_upper
      do i = coarse_i_lower, coarse_i_upper
        if (coarse_geometry%cell_type(i, j) == eb_covered_cell) then
          candidate_coarse_state(:, i, j) = coarse_state(:, i, j)
          cycle
        end if
        if (coarse_temperature(i, j) <= 0.0_dp) return
        call reactive_conserved_to_primitive( &
          species, candidate_coarse_state(:, i, j), &
          coarse_temperature(i, j), primitive, recovered_temperature, &
          sound_speed, local_ok)
        if (.not. local_ok) return
        candidate_coarse_temperature(i, j) = recovered_temperature
      end do
    end do
    do j = 1, fine_geometry%ny
      do i = 1, fine_geometry%nx
        if (fine_geometry%cell_type(i, j) == eb_covered_cell) then
          candidate_fine_state(:, i, j) = fine_state(:, i, j)
          cycle
        end if
        if (fine_temperature(i, j) <= 0.0_dp) return
        call reactive_conserved_to_primitive( &
          species, candidate_fine_state(:, i, j), fine_temperature(i, j), &
          primitive, recovered_temperature, sound_speed, local_ok)
        if (.not. local_ok) return
        candidate_fine_temperature(i, j) = recovered_temperature
      end do
    end do

    refluxed_coarse_state = candidate_coarse_state
    refluxed_coarse_temperature = candidate_coarse_temperature
    refluxed_fine_state = candidate_fine_state
    refluxed_fine_temperature = candidate_fine_temperature
    register = candidate_register
    ok = .true.
  end subroutine reflux_reactive_eb_state_patch_support_2d

  subroutine reflux_reactive_eb_state_patch_2d( &
      species, coarse_state, coarse_temperature, coarse_geometry, &
      fine_state, fine_temperature, fine_geometry, patch, register, &
      refluxed_coarse_state, refluxed_coarse_temperature, &
      refluxed_fine_state, refluxed_fine_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: coarse_state(:, :, :), coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    real(dp), intent(in) :: fine_state(:, :, :), fine_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    type(amr_eb_flux_register_2d), intent(inout) :: register
    real(dp), intent(out) :: refluxed_coarse_state(:, :, :)
    real(dp), intent(out) :: refluxed_coarse_temperature(:, :)
    real(dp), intent(out) :: refluxed_fine_state(:, :, :)
    real(dp), intent(out) :: refluxed_fine_temperature(:, :)
    logical, intent(out) :: ok

    call reflux_reactive_eb_state_patch_support_2d( &
      species, 1, 1, coarse_state, coarse_temperature, coarse_geometry, &
      fine_state, fine_temperature, fine_geometry, patch, register, &
      refluxed_coarse_state, refluxed_coarse_temperature, &
      refluxed_fine_state, refluxed_fine_temperature, ok)
  end subroutine reflux_reactive_eb_state_patch_2d

end module amr_eb_flux_register_2d_mod
