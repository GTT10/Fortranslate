module amr_multipatch_1d_mod
  use precision_mod, only: dp
  use amr_hierarchy_1d_mod, only: &
    amr_two_level_hierarchy_1d, amr_level_field_1d, &
    amr_flux_register_1d, initialize_two_level_hierarchy_1d, &
    prolong_conservative_1d, average_down_1d, &
    initialize_flux_register_1d, reflux_1d
  implicit none
  private

  type, public :: amr_patch_set_1d
    integer :: coarse_level = -1
    integer :: coarse_cells = 0
    integer :: refinement_ratio = 0
    real(dp) :: x_lower = 0.0_dp
    real(dp) :: x_upper = 0.0_dp
    real(dp) :: coarse_dx = 0.0_dp
    real(dp) :: fine_dx = 0.0_dp
    type(amr_two_level_hierarchy_1d), allocatable :: patches(:)
  contains
    procedure :: patch_count => amr_patch_set_patch_count
    procedure :: covered_coarse_cell_count => &
      amr_patch_set_covered_coarse_cell_count
    procedure :: fine_cell_count => amr_patch_set_fine_cell_count
    procedure :: parent_cell_is_covered => amr_patch_set_parent_cell_is_covered
    procedure :: is_valid => amr_patch_set_is_valid
  end type amr_patch_set_1d

  public :: initialize_patch_set_1d
  public :: prolong_patch_set_1d
  public :: average_down_patch_set_1d
  public :: initialize_patch_flux_registers_1d
  public :: reflux_patch_set_1d
  public :: synchronize_patch_set_1d
  public :: composite_integral_patch_set_1d

contains

  pure integer function amr_patch_set_patch_count(self) result(count)
    class(amr_patch_set_1d), intent(in) :: self

    count = 0
    if (allocated(self%patches)) count = size(self%patches)
  end function amr_patch_set_patch_count

  pure integer function amr_patch_set_covered_coarse_cell_count(self) &
      result(count)
    class(amr_patch_set_1d), intent(in) :: self

    integer :: patch

    count = 0
    if (.not. allocated(self%patches)) return
    do patch = 1, size(self%patches)
      count = count + self%patches(patch)%covered_coarse_cells()
    end do
  end function amr_patch_set_covered_coarse_cell_count

  pure integer function amr_patch_set_fine_cell_count(self) result(count)
    class(amr_patch_set_1d), intent(in) :: self

    integer :: patch

    count = 0
    if (.not. allocated(self%patches)) return
    do patch = 1, size(self%patches)
      count = count + self%patches(patch)%fine%cell_count()
    end do
  end function amr_patch_set_fine_cell_count

  pure logical function amr_patch_set_parent_cell_is_covered( &
      self, cell) result(covered)
    class(amr_patch_set_1d), intent(in) :: self
    integer, intent(in) :: cell

    integer :: patch

    covered = .false.
    if (cell < 1 .or. cell > self%coarse_cells) return
    if (.not. allocated(self%patches)) return
    do patch = 1, size(self%patches)
      if (cell >= self%patches(patch)%fine_coarse_lower .and. &
          cell <= self%patches(patch)%fine_coarse_upper) then
        covered = .true.
        return
      end if
    end do
  end function amr_patch_set_parent_cell_is_covered

  pure logical function amr_patch_set_is_valid(self) result(valid)
    class(amr_patch_set_1d), intent(in) :: self

    real(dp) :: spacing_tolerance
    integer :: patch, previous_upper

    valid = self%coarse_level >= 0 .and. self%coarse_cells >= 1 .and. &
      self%refinement_ratio >= 2 .and. self%x_upper > self%x_lower .and. &
      self%coarse_dx > 0.0_dp .and. self%fine_dx > 0.0_dp .and. &
      allocated(self%patches)
    if (.not. valid) return
    spacing_tolerance = 32.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, abs(self%coarse_dx), abs(self%fine_dx))
    valid = abs(self%fine_dx * real(self%refinement_ratio, dp) - &
      self%coarse_dx) <= spacing_tolerance
    if (.not. valid) return

    previous_upper = -1
    do patch = 1, size(self%patches)
      valid = self%patches(patch)%is_valid() .and. &
        self%patches(patch)%coarse%level == self%coarse_level .and. &
        self%patches(patch)%coarse%cell_count() == self%coarse_cells .and. &
        self%patches(patch)%refinement_ratio == self%refinement_ratio .and. &
        abs(self%patches(patch)%x_lower - self%x_lower) <= &
          spacing_tolerance .and. &
        abs(self%patches(patch)%x_upper - self%x_upper) <= &
          spacing_tolerance
      if (.not. valid) return
      if (patch > 1) then
        valid = self%patches(patch)%fine_coarse_lower > previous_upper + 1
        if (.not. valid) return
      end if
      previous_upper = self%patches(patch)%fine_coarse_upper
    end do
  end function amr_patch_set_is_valid

  subroutine initialize_patch_set_1d( &
      coarse_cells, patch_parent_lower, patch_parent_upper, &
      refinement_ratio, x_lower, x_upper, patch_set, ok, coarse_level)
    integer, intent(in) :: coarse_cells
    integer, intent(in) :: patch_parent_lower(:), patch_parent_upper(:)
    integer, intent(in) :: refinement_ratio
    real(dp), intent(in) :: x_lower, x_upper
    type(amr_patch_set_1d), intent(out) :: patch_set
    logical, intent(out) :: ok
    integer, intent(in), optional :: coarse_level

    logical :: local_ok
    integer :: parent_level, patch

    parent_level = 0
    if (present(coarse_level)) parent_level = coarse_level
    patch_set%coarse_level = parent_level
    patch_set%coarse_cells = coarse_cells
    patch_set%refinement_ratio = refinement_ratio
    patch_set%x_lower = x_lower
    patch_set%x_upper = x_upper
    if (coarse_cells > 0 .and. refinement_ratio > 0 .and. &
        x_upper > x_lower) then
      patch_set%coarse_dx = &
        (x_upper - x_lower) / real(coarse_cells, dp)
      patch_set%fine_dx = &
        patch_set%coarse_dx / real(refinement_ratio, dp)
    end if
    ok = size(patch_parent_lower) == size(patch_parent_upper)
    if (.not. ok) return
    allocate(patch_set%patches(size(patch_parent_lower)))
    do patch = 1, size(patch_set%patches)
      call initialize_two_level_hierarchy_1d( &
        coarse_cells, patch_parent_lower(patch), &
        patch_parent_upper(patch), refinement_ratio, x_lower, x_upper, &
        patch_set%patches(patch), local_ok, parent_level)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
    end do
    ok = patch_set%is_valid()
  end subroutine initialize_patch_set_1d

  subroutine prolong_patch_set_1d(coarse, patch_set, fine_fields, ok)
    real(dp), intent(in) :: coarse(:, :)
    type(amr_patch_set_1d), intent(in) :: patch_set
    type(amr_level_field_1d), allocatable, intent(out) :: fine_fields(:)
    logical, intent(out) :: ok

    logical :: local_ok
    integer :: patch, variable_count

    ok = patch_set%is_valid() .and. size(coarse, 1) >= 1 .and. &
      size(coarse, 2) == patch_set%coarse_cells
    if (.not. ok) return
    variable_count = size(coarse, 1)
    allocate(fine_fields(patch_set%patch_count()))
    do patch = 1, size(fine_fields)
      allocate(fine_fields(patch)%values( &
        variable_count, patch_set%patches(patch)%fine%cell_count()))
      call prolong_conservative_1d( &
        coarse, patch_set%patches(patch), &
        fine_fields(patch)%values, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
    end do
    ok = valid_patch_fields(fine_fields, patch_set, variable_count)
  end subroutine prolong_patch_set_1d

  subroutine average_down_patch_set_1d( &
      fine_fields, patch_set, coarse, ok)
    type(amr_level_field_1d), intent(in) :: fine_fields(:)
    type(amr_patch_set_1d), intent(in) :: patch_set
    real(dp), intent(inout) :: coarse(:, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: coarse_backup(:, :)
    logical :: local_ok
    integer :: patch

    ok = size(coarse, 1) >= 1 .and. &
      size(coarse, 2) == patch_set%coarse_cells .and. &
      valid_patch_fields(fine_fields, patch_set, size(coarse, 1))
    if (.not. ok) return
    coarse_backup = coarse
    do patch = 1, size(fine_fields)
      call average_down_1d( &
        fine_fields(patch)%values, patch_set%patches(patch), &
        coarse, local_ok)
      if (.not. local_ok) then
        coarse = coarse_backup
        ok = .false.
        return
      end if
    end do
    ok = .true.
  end subroutine average_down_patch_set_1d

  subroutine initialize_patch_flux_registers_1d( &
      patch_set, variable_count, registers, ok)
    type(amr_patch_set_1d), intent(in) :: patch_set
    integer, intent(in) :: variable_count
    type(amr_flux_register_1d), allocatable, intent(out) :: registers(:)
    logical, intent(out) :: ok

    logical :: local_ok
    integer :: patch

    ok = patch_set%is_valid() .and. variable_count >= 1
    if (.not. ok) return
    allocate(registers(patch_set%patch_count()))
    do patch = 1, size(registers)
      call initialize_flux_register_1d( &
        registers(patch), variable_count, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
    end do
    ok = .true.
  end subroutine initialize_patch_flux_registers_1d

  subroutine reflux_patch_set_1d(coarse, patch_set, registers, ok)
    real(dp), intent(inout) :: coarse(:, :)
    type(amr_patch_set_1d), intent(in) :: patch_set
    type(amr_flux_register_1d), intent(inout) :: registers(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: coarse_backup(:, :)
    type(amr_flux_register_1d), allocatable :: register_backup(:)
    logical :: local_ok
    integer :: patch

    ok = patch_set%is_valid() .and. size(coarse, 1) >= 1 .and. &
      size(coarse, 2) == patch_set%coarse_cells .and. &
      size(registers) == patch_set%patch_count()
    if (.not. ok) return
    coarse_backup = coarse
    register_backup = registers
    do patch = 1, size(registers)
      call reflux_1d( &
        coarse, patch_set%patches(patch), registers(patch), local_ok)
      if (.not. local_ok) then
        coarse = coarse_backup
        registers = register_backup
        ok = .false.
        return
      end if
    end do
    ok = .true.
  end subroutine reflux_patch_set_1d

  subroutine synchronize_patch_set_1d( &
      coarse, fine_fields, patch_set, registers, ok)
    real(dp), intent(inout) :: coarse(:, :)
    type(amr_level_field_1d), intent(in) :: fine_fields(:)
    type(amr_patch_set_1d), intent(in) :: patch_set
    type(amr_flux_register_1d), intent(inout) :: registers(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: coarse_backup(:, :)
    type(amr_flux_register_1d), allocatable :: register_backup(:)

    ok = valid_patch_fields(fine_fields, patch_set, size(coarse, 1))
    if (.not. ok) return
    coarse_backup = coarse
    register_backup = registers
    call reflux_patch_set_1d(coarse, patch_set, registers, ok)
    if (.not. ok) then
      coarse = coarse_backup
      registers = register_backup
      return
    end if
    call average_down_patch_set_1d(fine_fields, patch_set, coarse, ok)
    if (.not. ok) then
      coarse = coarse_backup
      registers = register_backup
    end if
  end subroutine synchronize_patch_set_1d

  pure subroutine composite_integral_patch_set_1d( &
      coarse, fine_fields, patch_set, integral, ok)
    real(dp), intent(in) :: coarse(:, :)
    type(amr_level_field_1d), intent(in) :: fine_fields(:)
    type(amr_patch_set_1d), intent(in) :: patch_set
    real(dp), intent(out) :: integral(:)
    logical, intent(out) :: ok

    integer :: component, patch, lower, upper

    integral = 0.0_dp
    ok = size(coarse, 1) >= 1 .and. &
      size(coarse, 2) == patch_set%coarse_cells .and. &
      size(integral) == size(coarse, 1) .and. &
      valid_patch_fields(fine_fields, patch_set, size(coarse, 1))
    if (.not. ok) return

    do component = 1, size(integral)
      integral(component) = patch_set%coarse_dx * &
        sum(coarse(component, :))
      do patch = 1, size(fine_fields)
        lower = patch_set%patches(patch)%fine_coarse_lower
        upper = patch_set%patches(patch)%fine_coarse_upper
        integral(component) = integral(component) - &
          patch_set%coarse_dx * sum(coarse(component, lower:upper)) + &
          patch_set%fine_dx * sum(fine_fields(patch)%values(component, :))
      end do
    end do
    ok = .true.
  end subroutine composite_integral_patch_set_1d

  pure logical function valid_patch_fields( &
      fine_fields, patch_set, variable_count) result(valid)
    type(amr_level_field_1d), intent(in) :: fine_fields(:)
    type(amr_patch_set_1d), intent(in) :: patch_set
    integer, intent(in) :: variable_count

    integer :: patch

    valid = patch_set%is_valid() .and. variable_count >= 1 .and. &
      size(fine_fields) == patch_set%patch_count()
    if (.not. valid) return
    do patch = 1, size(fine_fields)
      valid = allocated(fine_fields(patch)%values)
      if (.not. valid) return
      valid = size(fine_fields(patch)%values, 1) == variable_count .and. &
        size(fine_fields(patch)%values, 2) == &
          patch_set%patches(patch)%fine%cell_count()
      if (.not. valid) return
    end do
  end function valid_patch_fields

end module amr_multipatch_1d_mod
