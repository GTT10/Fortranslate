module amr_hierarchy_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  implicit none
  private

  type, public :: amr_patch_3d
    integer :: coarse_nx = 0
    integer :: coarse_ny = 0
    integer :: coarse_nz = 0
    integer :: coarse_i_lower = 1
    integer :: coarse_i_upper = 0
    integer :: coarse_j_lower = 1
    integer :: coarse_j_upper = 0
    integer :: coarse_k_lower = 1
    integer :: coarse_k_upper = 0
    integer :: refinement_ratio = 0
  contains
    procedure :: fine_nx => patch_fine_nx
    procedure :: fine_ny => patch_fine_ny
    procedure :: fine_nz => patch_fine_nz
    procedure :: is_valid => patch_is_valid
    procedure :: is_strictly_interior => patch_is_strictly_interior
  end type amr_patch_3d

  public :: initialize_amr_patch_3d
  public :: prolong_pcm_3d
  public :: restrict_average_3d
  public :: average_down_3d
  public :: composite_integrals_amr_3d

contains

  pure integer function patch_fine_nx(self) result(count)
    class(amr_patch_3d), intent(in) :: self
    count = max(0, self%coarse_i_upper - self%coarse_i_lower + 1) * &
      max(0, self%refinement_ratio)
  end function patch_fine_nx

  pure integer function patch_fine_ny(self) result(count)
    class(amr_patch_3d), intent(in) :: self
    count = max(0, self%coarse_j_upper - self%coarse_j_lower + 1) * &
      max(0, self%refinement_ratio)
  end function patch_fine_ny

  pure integer function patch_fine_nz(self) result(count)
    class(amr_patch_3d), intent(in) :: self
    count = max(0, self%coarse_k_upper - self%coarse_k_lower + 1) * &
      max(0, self%refinement_ratio)
  end function patch_fine_nz

  pure logical function patch_is_valid(self) result(valid)
    class(amr_patch_3d), intent(in) :: self

    valid = self%coarse_nx >= 2 .and. self%coarse_ny >= 2 .and. &
      self%coarse_nz >= 2 .and. self%refinement_ratio >= 2 .and. &
      self%coarse_i_lower >= 1 .and. &
      self%coarse_i_upper <= self%coarse_nx .and. &
      self%coarse_i_upper >= self%coarse_i_lower .and. &
      self%coarse_j_lower >= 1 .and. &
      self%coarse_j_upper <= self%coarse_ny .and. &
      self%coarse_j_upper >= self%coarse_j_lower .and. &
      self%coarse_k_lower >= 1 .and. &
      self%coarse_k_upper <= self%coarse_nz .and. &
      self%coarse_k_upper >= self%coarse_k_lower
  end function patch_is_valid

  pure logical function patch_is_strictly_interior(self) result(interior)
    class(amr_patch_3d), intent(in) :: self

    interior = self%is_valid() .and. self%coarse_i_lower > 1 .and. &
      self%coarse_i_upper < self%coarse_nx .and. &
      self%coarse_j_lower > 1 .and. &
      self%coarse_j_upper < self%coarse_ny .and. &
      self%coarse_k_lower > 1 .and. &
      self%coarse_k_upper < self%coarse_nz
  end function patch_is_strictly_interior

  pure subroutine initialize_amr_patch_3d( &
      coarse_nx, coarse_ny, coarse_nz, coarse_i_lower, coarse_i_upper, &
      coarse_j_lower, coarse_j_upper, coarse_k_lower, coarse_k_upper, &
      refinement_ratio, patch, ok)
    integer, intent(in) :: coarse_nx, coarse_ny, coarse_nz
    integer, intent(in) :: coarse_i_lower, coarse_i_upper
    integer, intent(in) :: coarse_j_lower, coarse_j_upper
    integer, intent(in) :: coarse_k_lower, coarse_k_upper
    integer, intent(in) :: refinement_ratio
    type(amr_patch_3d), intent(out) :: patch
    logical, intent(out) :: ok

    patch%coarse_nx = coarse_nx
    patch%coarse_ny = coarse_ny
    patch%coarse_nz = coarse_nz
    patch%coarse_i_lower = coarse_i_lower
    patch%coarse_i_upper = coarse_i_upper
    patch%coarse_j_lower = coarse_j_lower
    patch%coarse_j_upper = coarse_j_upper
    patch%coarse_k_lower = coarse_k_lower
    patch%coarse_k_upper = coarse_k_upper
    patch%refinement_ratio = refinement_ratio
    ok = patch%is_valid()
  end subroutine initialize_amr_patch_3d

  subroutine prolong_pcm_3d(coarse, patch, fine, ok)
    real(dp), intent(in) :: coarse(:, :, :, :)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(out) :: fine(:, :, :, :)
    logical, intent(out) :: ok

    integer :: i, j, k, coarse_i, coarse_j, coarse_k, ratio

    fine = 0.0_dp
    ok = patch%is_valid() .and. size(coarse, 1) == size(fine, 1) .and. &
      size(coarse, 2) == patch%coarse_nx .and. &
      size(coarse, 3) == patch%coarse_ny .and. &
      size(coarse, 4) == patch%coarse_nz .and. &
      size(fine, 2) == patch%fine_nx() .and. &
      size(fine, 3) == patch%fine_ny() .and. &
      size(fine, 4) == patch%fine_nz()
    if (.not. ok) return
    ratio = patch%refinement_ratio
    do k = 1, patch%fine_nz()
      coarse_k = patch%coarse_k_lower + (k - 1) / ratio
      do j = 1, patch%fine_ny()
        coarse_j = patch%coarse_j_lower + (j - 1) / ratio
        do i = 1, patch%fine_nx()
          coarse_i = patch%coarse_i_lower + (i - 1) / ratio
          fine(:, i, j, k) = coarse(:, coarse_i, coarse_j, coarse_k)
        end do
      end do
    end do
    ok = all(ieee_is_finite(fine))
  end subroutine prolong_pcm_3d

  subroutine restrict_average_3d(fine, patch, restricted, ok)
    real(dp), intent(in) :: fine(:, :, :, :)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(out) :: restricted(:, :, :, :)
    logical, intent(out) :: ok

    integer :: coarse_i, coarse_j, coarse_k, child_i, child_j, child_k
    integer :: i_lower, j_lower, k_lower, ratio
    real(dp) :: inverse_children

    restricted = 0.0_dp
    ok = patch%is_valid() .and. size(fine, 1) == size(restricted, 1) .and. &
      size(fine, 2) == patch%fine_nx() .and. &
      size(fine, 3) == patch%fine_ny() .and. &
      size(fine, 4) == patch%fine_nz() .and. &
      size(restricted, 2) == &
        patch%coarse_i_upper - patch%coarse_i_lower + 1 .and. &
      size(restricted, 3) == &
        patch%coarse_j_upper - patch%coarse_j_lower + 1 .and. &
      size(restricted, 4) == &
        patch%coarse_k_upper - patch%coarse_k_lower + 1
    if (.not. ok) return
    ratio = patch%refinement_ratio
    inverse_children = 1.0_dp / real(ratio**3, dp)
    do coarse_k = 1, size(restricted, 4)
      k_lower = (coarse_k - 1) * ratio + 1
      do coarse_j = 1, size(restricted, 3)
        j_lower = (coarse_j - 1) * ratio + 1
        do coarse_i = 1, size(restricted, 2)
          i_lower = (coarse_i - 1) * ratio + 1
          do child_k = k_lower, k_lower + ratio - 1
            do child_j = j_lower, j_lower + ratio - 1
              do child_i = i_lower, i_lower + ratio - 1
                restricted(:, coarse_i, coarse_j, coarse_k) = &
                  restricted(:, coarse_i, coarse_j, coarse_k) + &
                  inverse_children * fine(:, child_i, child_j, child_k)
              end do
            end do
          end do
        end do
      end do
    end do
    ok = all(ieee_is_finite(restricted))
  end subroutine restrict_average_3d

  subroutine average_down_3d(coarse, fine, patch, ok)
    real(dp), intent(inout) :: coarse(:, :, :, :)
    real(dp), intent(in) :: fine(:, :, :, :)
    type(amr_patch_3d), intent(in) :: patch
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate(:, :, :, :)
    real(dp), allocatable :: restricted(:, :, :, :)
    integer :: covered_nx, covered_ny, covered_nz

    ok = .false.
    if (.not. patch%is_valid() .or. &
        size(coarse, 2) /= patch%coarse_nx .or. &
        size(coarse, 3) /= patch%coarse_ny .or. &
        size(coarse, 4) /= patch%coarse_nz) return
    covered_nx = patch%coarse_i_upper - patch%coarse_i_lower + 1
    covered_ny = patch%coarse_j_upper - patch%coarse_j_lower + 1
    covered_nz = patch%coarse_k_upper - patch%coarse_k_lower + 1
    allocate(restricted(size(coarse, 1), covered_nx, covered_ny, covered_nz))
    call restrict_average_3d(fine, patch, restricted, ok)
    if (.not. ok) return
    allocate(candidate, source=coarse)
    candidate(:, &
      patch%coarse_i_lower:patch%coarse_i_upper, &
      patch%coarse_j_lower:patch%coarse_j_upper, &
      patch%coarse_k_lower:patch%coarse_k_upper) = restricted
    coarse = candidate
    ok = .true.
  end subroutine average_down_3d

  subroutine composite_integrals_amr_3d( &
      coarse, fine, patch, dx, dy, dz, integrals, ok)
    real(dp), intent(in) :: coarse(:, :, :, :), fine(:, :, :, :)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: dx, dy, dz
    real(dp), intent(out) :: integrals(:)
    logical, intent(out) :: ok

    real(dp) :: coarse_volume, fine_volume
    integer :: component

    integrals = 0.0_dp
    ok = .false.
    if (.not. patch%is_valid()) return
    if (size(coarse, 1) /= size(fine, 1) .or. &
        size(integrals) /= size(coarse, 1) .or. &
        size(coarse, 2) /= patch%coarse_nx .or. &
        size(coarse, 3) /= patch%coarse_ny .or. &
        size(coarse, 4) /= patch%coarse_nz .or. &
        size(fine, 2) /= patch%fine_nx() .or. &
        size(fine, 3) /= patch%fine_ny() .or. &
        size(fine, 4) /= patch%fine_nz()) return
    if (.not. all(ieee_is_finite([dx, dy, dz]))) return
    if (dx <= 0.0_dp .or. dy <= 0.0_dp .or. dz <= 0.0_dp) return
    if (.not. all(ieee_is_finite(coarse))) return
    if (.not. all(ieee_is_finite(fine))) return
    coarse_volume = dx * dy * dz
    fine_volume = coarse_volume / real(patch%refinement_ratio**3, dp)
    do component = 1, size(integrals)
      integrals(component) = coarse_volume * ( &
        sum(coarse(component, :, :, :)) - sum(coarse(component, &
          patch%coarse_i_lower:patch%coarse_i_upper, &
          patch%coarse_j_lower:patch%coarse_j_upper, &
          patch%coarse_k_lower:patch%coarse_k_upper))) + &
        fine_volume * sum(fine(component, :, :, :))
    end do
    ok = all(ieee_is_finite(integrals))
  end subroutine composite_integrals_amr_3d

end module amr_hierarchy_3d_mod
