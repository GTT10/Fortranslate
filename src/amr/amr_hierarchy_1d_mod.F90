module amr_hierarchy_1d_mod
  use precision_mod, only: dp
  use slope_limiter_mod, only: limited_slope
  implicit none
  private

  type, public :: amr_box_1d
    integer :: level = -1
    integer :: lower = 1
    integer :: upper = 0
  contains
    procedure :: cell_count => amr_box_cell_count
    procedure :: is_valid => amr_box_is_valid
  end type amr_box_1d

  type, public :: amr_two_level_hierarchy_1d
    type(amr_box_1d) :: coarse
    type(amr_box_1d) :: fine
    integer :: refinement_ratio = 0
    integer :: fine_coarse_lower = 1
    integer :: fine_coarse_upper = 0
    real(dp) :: x_lower = 0.0_dp
    real(dp) :: x_upper = 0.0_dp
    real(dp) :: coarse_dx = 0.0_dp
    real(dp) :: fine_dx = 0.0_dp
  contains
    procedure :: covered_coarse_cells => amr_covered_coarse_cells
    procedure :: is_valid => amr_hierarchy_is_valid
  end type amr_two_level_hierarchy_1d

  type, public :: amr_level_field_1d
    real(dp), allocatable :: values(:, :)
  end type amr_level_field_1d

  type, public :: amr_multilevel_hierarchy_1d
    integer :: base_cells = 0
    real(dp) :: x_lower = 0.0_dp
    real(dp) :: x_upper = 0.0_dp
    type(amr_two_level_hierarchy_1d), allocatable :: interfaces(:)
  contains
    procedure :: level_count => amr_multilevel_level_count
    procedure :: level_cell_count => amr_multilevel_level_cell_count
    procedure :: level_dx => amr_multilevel_level_dx
    procedure :: level_bounds => amr_multilevel_level_bounds
    procedure :: is_valid => amr_multilevel_is_valid
  end type amr_multilevel_hierarchy_1d

  type, public :: amr_flux_register_1d
    real(dp), allocatable :: left(:)
    real(dp), allocatable :: right(:)
  contains
    procedure :: reset => reset_flux_register_1d
  end type amr_flux_register_1d

  public :: initialize_two_level_hierarchy_1d
  public :: initialize_multilevel_hierarchy_1d
  public :: prolong_conservative_1d
  public :: prolong_multilevel_1d
  public :: restrict_average_1d
  public :: average_down_1d
  public :: average_down_multilevel_1d
  public :: level_subcycle_time_steps_1d
  public :: multilevel_subcycle_counts_1d
  public :: multilevel_subcycle_time_steps_1d
  public :: initialize_flux_register_1d
  public :: initialize_multilevel_flux_registers_1d
  public :: accumulate_coarse_flux_1d
  public :: accumulate_fine_flux_1d
  public :: reflux_1d
  public :: composite_integral_1d
  public :: synchronize_multilevel_1d
  public :: composite_integral_multilevel_1d

contains

  pure integer function amr_box_cell_count(self) result(count)
    class(amr_box_1d), intent(in) :: self

    count = max(0, self%upper - self%lower + 1)
  end function amr_box_cell_count

  pure logical function amr_box_is_valid(self) result(valid)
    class(amr_box_1d), intent(in) :: self

    valid = self%level >= 0 .and. self%lower >= 1 .and. &
      self%upper >= self%lower
  end function amr_box_is_valid

  pure integer function amr_multilevel_level_count(self) result(count)
    class(amr_multilevel_hierarchy_1d), intent(in) :: self

    count = 0
    if (allocated(self%interfaces)) count = size(self%interfaces) + 1
  end function amr_multilevel_level_count

  pure integer function amr_multilevel_level_cell_count(self, level) &
      result(count)
    class(amr_multilevel_hierarchy_1d), intent(in) :: self
    integer, intent(in) :: level

    count = 0
    if (level < 0 .or. level >= self%level_count()) return
    if (level == 0) then
      count = self%base_cells
    else
      count = self%interfaces(level)%fine%cell_count()
    end if
  end function amr_multilevel_level_cell_count

  pure real(dp) function amr_multilevel_level_dx(self, level) result(dx)
    class(amr_multilevel_hierarchy_1d), intent(in) :: self
    integer, intent(in) :: level

    dx = 0.0_dp
    if (level < 0 .or. level >= self%level_count()) return
    if (level == 0) then
      if (self%base_cells > 0) &
        dx = (self%x_upper - self%x_lower) / real(self%base_cells, dp)
    else
      dx = self%interfaces(level)%fine_dx
    end if
  end function amr_multilevel_level_dx

  pure subroutine amr_multilevel_level_bounds( &
      self, level, x_lower, x_upper, ok)
    class(amr_multilevel_hierarchy_1d), intent(in) :: self
    integer, intent(in) :: level
    real(dp), intent(out) :: x_lower, x_upper
    logical, intent(out) :: ok

    type(amr_two_level_hierarchy_1d) :: relation

    x_lower = 0.0_dp
    x_upper = 0.0_dp
    ok = level >= 0 .and. level < self%level_count()
    if (.not. ok) return
    if (level == 0) then
      x_lower = self%x_lower
      x_upper = self%x_upper
    else
      relation = self%interfaces(level)
      x_lower = relation%x_lower + &
        real(relation%fine_coarse_lower - 1, dp) * relation%coarse_dx
      x_upper = relation%x_lower + &
        real(relation%fine_coarse_upper, dp) * relation%coarse_dx
    end if
    ok = x_upper > x_lower
  end subroutine amr_multilevel_level_bounds

  pure logical function amr_multilevel_is_valid(self) result(valid)
    class(amr_multilevel_hierarchy_1d), intent(in) :: self

    real(dp) :: expected_lower, expected_upper, child_lower, child_upper
    real(dp) :: coordinate_tolerance
    logical :: bounds_ok
    integer :: relation, expected_cells

    valid = allocated(self%interfaces)
    if (.not. valid) return
    valid = self%base_cells >= 1 .and. self%x_upper > self%x_lower
    if (.not. valid) return
    if (size(self%interfaces) == 0) return

    expected_cells = self%base_cells
    expected_lower = self%x_lower
    expected_upper = self%x_upper
    do relation = 1, size(self%interfaces)
      coordinate_tolerance = 64.0_dp * epsilon(1.0_dp) * &
        max(1.0_dp, abs(expected_lower), abs(expected_upper))
      valid = self%interfaces(relation)%is_valid() .and. &
        self%interfaces(relation)%coarse%level == relation - 1 .and. &
        self%interfaces(relation)%fine%level == relation .and. &
        self%interfaces(relation)%coarse%cell_count() == expected_cells .and. &
        abs(self%interfaces(relation)%x_lower - expected_lower) <= &
          coordinate_tolerance .and. &
        abs(self%interfaces(relation)%x_upper - expected_upper) <= &
          coordinate_tolerance
      if (.not. valid) return
      call self%level_bounds( &
        relation, child_lower, child_upper, bounds_ok)
      if (.not. bounds_ok) then
        valid = .false.
        return
      end if
      expected_cells = self%interfaces(relation)%fine%cell_count()
      expected_lower = child_lower
      expected_upper = child_upper
    end do
  end function amr_multilevel_is_valid

  pure integer function amr_covered_coarse_cells(self) result(count)
    class(amr_two_level_hierarchy_1d), intent(in) :: self

    count = max(0, self%fine_coarse_upper - self%fine_coarse_lower + 1)
  end function amr_covered_coarse_cells

  pure logical function amr_hierarchy_is_valid(self) result(valid)
    class(amr_two_level_hierarchy_1d), intent(in) :: self
    real(dp) :: spacing_tolerance

    spacing_tolerance = 32.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, abs(self%coarse_dx), abs(self%fine_dx))
    valid = self%coarse%is_valid() .and. self%fine%is_valid() .and. &
      self%fine%level == self%coarse%level + 1 .and. &
      self%refinement_ratio >= 2 .and. self%x_upper > self%x_lower .and. &
      self%coarse_dx > 0.0_dp .and. self%fine_dx > 0.0_dp .and. &
      self%fine_coarse_lower > self%coarse%lower .and. &
      self%fine_coarse_upper < self%coarse%upper .and. &
      self%fine_coarse_upper >= self%fine_coarse_lower .and. &
      self%fine%lower == &
        (self%fine_coarse_lower - 1) * self%refinement_ratio + 1 .and. &
      self%fine%upper == &
        self%fine_coarse_upper * self%refinement_ratio .and. &
      self%fine%cell_count() == &
        self%covered_coarse_cells() * self%refinement_ratio .and. &
      abs(self%fine_dx * real(self%refinement_ratio, dp) - &
        self%coarse_dx) <= spacing_tolerance
  end function amr_hierarchy_is_valid

  pure subroutine initialize_two_level_hierarchy_1d( &
      coarse_cells, fine_coarse_lower, fine_coarse_upper, &
      refinement_ratio, x_lower, x_upper, hierarchy, ok, coarse_level)
    integer, intent(in) :: coarse_cells
    integer, intent(in) :: fine_coarse_lower, fine_coarse_upper
    integer, intent(in) :: refinement_ratio
    real(dp), intent(in) :: x_lower, x_upper
    type(amr_two_level_hierarchy_1d), intent(out) :: hierarchy
    logical, intent(out) :: ok
    integer, intent(in), optional :: coarse_level

    integer :: parent_level

    parent_level = 0
    if (present(coarse_level)) parent_level = coarse_level
    hierarchy%coarse%level = parent_level
    hierarchy%coarse%lower = 1
    hierarchy%coarse%upper = coarse_cells
    hierarchy%fine%level = parent_level + 1
    hierarchy%fine%lower = &
      (fine_coarse_lower - 1) * refinement_ratio + 1
    hierarchy%fine%upper = fine_coarse_upper * refinement_ratio
    hierarchy%refinement_ratio = refinement_ratio
    hierarchy%fine_coarse_lower = fine_coarse_lower
    hierarchy%fine_coarse_upper = fine_coarse_upper
    hierarchy%x_lower = x_lower
    hierarchy%x_upper = x_upper
    hierarchy%coarse_dx = 0.0_dp
    hierarchy%fine_dx = 0.0_dp
    if (coarse_cells > 0 .and. refinement_ratio > 0 .and. &
        x_upper > x_lower) then
      hierarchy%coarse_dx = &
        (x_upper - x_lower) / real(coarse_cells, dp)
      hierarchy%fine_dx = &
        hierarchy%coarse_dx / real(refinement_ratio, dp)
    end if
    ok = hierarchy%is_valid()
  end subroutine initialize_two_level_hierarchy_1d

  subroutine initialize_multilevel_hierarchy_1d( &
      base_cells, patch_parent_lower, patch_parent_upper, &
      refinement_ratios, x_lower, x_upper, hierarchy, ok)
    integer, intent(in) :: base_cells
    integer, intent(in) :: patch_parent_lower(:), patch_parent_upper(:)
    integer, intent(in) :: refinement_ratios(:)
    real(dp), intent(in) :: x_lower, x_upper
    type(amr_multilevel_hierarchy_1d), intent(out) :: hierarchy
    logical, intent(out) :: ok

    real(dp) :: parent_lower, parent_upper, parent_dx
    logical :: local_ok
    integer :: parent_cells, relation

    hierarchy%base_cells = base_cells
    hierarchy%x_lower = x_lower
    hierarchy%x_upper = x_upper
    ok = base_cells >= 1 .and. x_upper > x_lower .and. &
      size(patch_parent_lower) == size(refinement_ratios) .and. &
      size(patch_parent_upper) == size(refinement_ratios)
    if (.not. ok) return
    allocate(hierarchy%interfaces(size(refinement_ratios)))
    if (size(refinement_ratios) == 0) then
      ok = hierarchy%is_valid()
      return
    end if

    parent_cells = base_cells
    parent_lower = x_lower
    parent_upper = x_upper
    do relation = 1, size(refinement_ratios)
      call initialize_two_level_hierarchy_1d( &
        parent_cells, patch_parent_lower(relation), &
        patch_parent_upper(relation), refinement_ratios(relation), &
        parent_lower, parent_upper, hierarchy%interfaces(relation), &
        local_ok, relation - 1)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      parent_dx = hierarchy%interfaces(relation)%coarse_dx
      parent_upper = parent_lower + &
        real(patch_parent_upper(relation), dp) * parent_dx
      parent_lower = parent_lower + &
        real(patch_parent_lower(relation) - 1, dp) * parent_dx
      parent_cells = hierarchy%interfaces(relation)%fine%cell_count()
    end do
    ok = hierarchy%is_valid()
  end subroutine initialize_multilevel_hierarchy_1d

  subroutine prolong_conservative_1d(coarse, hierarchy, fine, ok)
    real(dp), intent(in) :: coarse(:, :)
    type(amr_two_level_hierarchy_1d), intent(in) :: hierarchy
    real(dp), intent(out) :: fine(:, :)
    logical, intent(out) :: ok

    real(dp) :: delta_minus, delta_plus, slope, offset
    logical :: slope_ok
    integer :: component, coarse_cell, child, fine_cell

    fine = 0.0_dp
    ok = hierarchy%is_valid() .and. size(coarse, 1) == size(fine, 1) .and. &
      size(coarse, 1) >= 1 .and. &
      size(coarse, 2) == hierarchy%coarse%cell_count() .and. &
      size(fine, 2) == hierarchy%fine%cell_count()
    if (.not. ok) return

    do component = 1, size(coarse, 1)
      do coarse_cell = hierarchy%fine_coarse_lower, &
          hierarchy%fine_coarse_upper
        delta_minus = coarse(component, coarse_cell) - &
          coarse(component, coarse_cell - 1)
        delta_plus = coarse(component, coarse_cell + 1) - &
          coarse(component, coarse_cell)
        call limited_slope( &
          delta_minus, delta_plus, "mc", slope, slope_ok)
        if (.not. slope_ok) then
          ok = .false.
          return
        end if
        do child = 1, hierarchy%refinement_ratio
          fine_cell = &
            (coarse_cell - hierarchy%fine_coarse_lower) * &
            hierarchy%refinement_ratio + child
          offset = (real(child, dp) - 0.5_dp) / &
            real(hierarchy%refinement_ratio, dp) - 0.5_dp
          fine(component, fine_cell) = &
            coarse(component, coarse_cell) + slope * offset
        end do
      end do
    end do
  end subroutine prolong_conservative_1d

  subroutine prolong_multilevel_1d(root, hierarchy, fields, ok)
    real(dp), intent(in) :: root(:, :)
    type(amr_multilevel_hierarchy_1d), intent(in) :: hierarchy
    type(amr_level_field_1d), allocatable, intent(out) :: fields(:)
    logical, intent(out) :: ok

    logical :: local_ok
    integer :: relation, variable_count

    ok = hierarchy%is_valid() .and. size(root, 1) >= 1 .and. &
      size(root, 2) == hierarchy%base_cells
    if (.not. ok) return
    variable_count = size(root, 1)
    allocate(fields(hierarchy%level_count()))
    allocate(fields(1)%values(variable_count, hierarchy%base_cells))
    fields(1)%values = root
    do relation = 1, size(hierarchy%interfaces)
      allocate(fields(relation + 1)%values( &
        variable_count, hierarchy%interfaces(relation)%fine%cell_count()))
      call prolong_conservative_1d( &
        fields(relation)%values, hierarchy%interfaces(relation), &
        fields(relation + 1)%values, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
    end do
    ok = valid_multilevel_fields(fields, hierarchy)
  end subroutine prolong_multilevel_1d

  pure subroutine restrict_average_1d(fine, hierarchy, restricted, ok)
    real(dp), intent(in) :: fine(:, :)
    type(amr_two_level_hierarchy_1d), intent(in) :: hierarchy
    real(dp), intent(out) :: restricted(:, :)
    logical, intent(out) :: ok

    integer :: component, covered_cell, first_child, last_child

    restricted = 0.0_dp
    ok = hierarchy%is_valid() .and. &
      size(fine, 1) == size(restricted, 1) .and. &
      size(fine, 1) >= 1 .and. &
      size(fine, 2) == hierarchy%fine%cell_count() .and. &
      size(restricted, 2) == hierarchy%covered_coarse_cells()
    if (.not. ok) return

    do covered_cell = 1, hierarchy%covered_coarse_cells()
      first_child = &
        (covered_cell - 1) * hierarchy%refinement_ratio + 1
      last_child = first_child + hierarchy%refinement_ratio - 1
      do component = 1, size(fine, 1)
        restricted(component, covered_cell) = &
          sum(fine(component, first_child:last_child)) / &
          real(hierarchy%refinement_ratio, dp)
      end do
    end do
  end subroutine restrict_average_1d

  subroutine average_down_1d(fine, hierarchy, coarse, ok)
    real(dp), intent(in) :: fine(:, :)
    type(amr_two_level_hierarchy_1d), intent(in) :: hierarchy
    real(dp), intent(inout) :: coarse(:, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: restricted(:, :)

    ok = hierarchy%is_valid() .and. size(fine, 1) == size(coarse, 1) .and. &
      size(fine, 1) >= 1 .and. &
      size(fine, 2) == hierarchy%fine%cell_count() .and. &
      size(coarse, 2) == hierarchy%coarse%cell_count()
    if (.not. ok) return
    allocate(restricted(size(coarse, 1), hierarchy%covered_coarse_cells()))
    call restrict_average_1d(fine, hierarchy, restricted, ok)
    if (.not. ok) return
    coarse(:, hierarchy%fine_coarse_lower:hierarchy%fine_coarse_upper) = &
      restricted
  end subroutine average_down_1d

  subroutine average_down_multilevel_1d(fields, hierarchy, ok)
    type(amr_level_field_1d), intent(inout) :: fields(:)
    type(amr_multilevel_hierarchy_1d), intent(in) :: hierarchy
    logical, intent(out) :: ok

    logical :: local_ok
    integer :: relation

    ok = valid_multilevel_fields(fields, hierarchy)
    if (.not. ok) return
    do relation = size(hierarchy%interfaces), 1, -1
      call average_down_1d( &
        fields(relation + 1)%values, hierarchy%interfaces(relation), &
        fields(relation)%values, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
    end do
    ok = .true.
  end subroutine average_down_multilevel_1d

  pure subroutine level_subcycle_time_steps_1d( &
      hierarchy, coarse_time_step, fine_time_steps, ok)
    type(amr_two_level_hierarchy_1d), intent(in) :: hierarchy
    real(dp), intent(in) :: coarse_time_step
    real(dp), intent(out) :: fine_time_steps(:)
    logical, intent(out) :: ok

    fine_time_steps = 0.0_dp
    ok = hierarchy%is_valid() .and. coarse_time_step > 0.0_dp .and. &
      size(fine_time_steps) == hierarchy%refinement_ratio
    if (.not. ok) return
    fine_time_steps = coarse_time_step / real(hierarchy%refinement_ratio, dp)
  end subroutine level_subcycle_time_steps_1d

  pure subroutine multilevel_subcycle_counts_1d(hierarchy, counts, ok)
    type(amr_multilevel_hierarchy_1d), intent(in) :: hierarchy
    integer, intent(out) :: counts(:)
    logical, intent(out) :: ok

    integer :: relation, ratio

    counts = 0
    ok = hierarchy%is_valid() .and. size(counts) == hierarchy%level_count()
    if (.not. ok) return
    counts(1) = 1
    do relation = 1, size(hierarchy%interfaces)
      ratio = hierarchy%interfaces(relation)%refinement_ratio
      if (counts(relation) > huge(counts(1)) / ratio) then
        ok = .false.
        return
      end if
      counts(relation + 1) = counts(relation) * ratio
    end do
  end subroutine multilevel_subcycle_counts_1d

  pure subroutine multilevel_subcycle_time_steps_1d( &
      hierarchy, coarse_time_step, time_steps, ok)
    type(amr_multilevel_hierarchy_1d), intent(in) :: hierarchy
    real(dp), intent(in) :: coarse_time_step
    real(dp), intent(out) :: time_steps(:)
    logical, intent(out) :: ok

    integer :: relation

    time_steps = 0.0_dp
    ok = hierarchy%is_valid() .and. coarse_time_step > 0.0_dp .and. &
      size(time_steps) == hierarchy%level_count()
    if (.not. ok) return
    time_steps(1) = coarse_time_step
    do relation = 1, size(hierarchy%interfaces)
      time_steps(relation + 1) = time_steps(relation) / &
        real(hierarchy%interfaces(relation)%refinement_ratio, dp)
    end do
  end subroutine multilevel_subcycle_time_steps_1d

  subroutine initialize_flux_register_1d(register, variable_count, ok)
    type(amr_flux_register_1d), intent(out) :: register
    integer, intent(in) :: variable_count
    logical, intent(out) :: ok

    ok = variable_count >= 1
    if (.not. ok) return
    allocate(register%left(variable_count), register%right(variable_count))
    call register%reset()
  end subroutine initialize_flux_register_1d

  subroutine initialize_multilevel_flux_registers_1d( &
      hierarchy, variable_count, registers, ok)
    type(amr_multilevel_hierarchy_1d), intent(in) :: hierarchy
    integer, intent(in) :: variable_count
    type(amr_flux_register_1d), allocatable, intent(out) :: registers(:)
    logical, intent(out) :: ok

    logical :: local_ok
    integer :: relation

    ok = hierarchy%is_valid() .and. variable_count >= 1
    if (.not. ok) return
    allocate(registers(size(hierarchy%interfaces)))
    do relation = 1, size(registers)
      call initialize_flux_register_1d( &
        registers(relation), variable_count, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
    end do
    ok = .true.
  end subroutine initialize_multilevel_flux_registers_1d

  subroutine reset_flux_register_1d(self)
    class(amr_flux_register_1d), intent(inout) :: self

    if (allocated(self%left)) self%left = 0.0_dp
    if (allocated(self%right)) self%right = 0.0_dp
  end subroutine reset_flux_register_1d

  subroutine accumulate_coarse_flux_1d( &
      register, left_flux, right_flux, time_step, ok)
    type(amr_flux_register_1d), intent(inout) :: register
    real(dp), intent(in) :: left_flux(:), right_flux(:), time_step
    logical, intent(out) :: ok

    ok = valid_flux_register(register, size(left_flux)) .and. &
      size(right_flux) == size(left_flux) .and. time_step >= 0.0_dp
    if (.not. ok) return
    register%left = register%left - time_step * left_flux
    register%right = register%right - time_step * right_flux
  end subroutine accumulate_coarse_flux_1d

  subroutine accumulate_fine_flux_1d( &
      register, left_flux, right_flux, time_step, ok)
    type(amr_flux_register_1d), intent(inout) :: register
    real(dp), intent(in) :: left_flux(:), right_flux(:), time_step
    logical, intent(out) :: ok

    ok = valid_flux_register(register, size(left_flux)) .and. &
      size(right_flux) == size(left_flux) .and. time_step >= 0.0_dp
    if (.not. ok) return
    register%left = register%left + time_step * left_flux
    register%right = register%right + time_step * right_flux
  end subroutine accumulate_fine_flux_1d

  subroutine reflux_1d(coarse, hierarchy, register, ok)
    real(dp), intent(inout) :: coarse(:, :)
    type(amr_two_level_hierarchy_1d), intent(in) :: hierarchy
    type(amr_flux_register_1d), intent(inout) :: register
    logical, intent(out) :: ok

    ok = hierarchy%is_valid() .and. &
      size(coarse, 2) == hierarchy%coarse%cell_count() .and. &
      valid_flux_register(register, size(coarse, 1))
    if (.not. ok) return
    coarse(:, hierarchy%fine_coarse_lower - 1) = &
      coarse(:, hierarchy%fine_coarse_lower - 1) - &
      register%left / hierarchy%coarse_dx
    coarse(:, hierarchy%fine_coarse_upper + 1) = &
      coarse(:, hierarchy%fine_coarse_upper + 1) + &
      register%right / hierarchy%coarse_dx
    call register%reset()
  end subroutine reflux_1d

  pure subroutine composite_integral_1d( &
      coarse, fine, hierarchy, integral, ok)
    real(dp), intent(in) :: coarse(:, :), fine(:, :)
    type(amr_two_level_hierarchy_1d), intent(in) :: hierarchy
    real(dp), intent(out) :: integral(:)
    logical, intent(out) :: ok

    integer :: component

    integral = 0.0_dp
    ok = hierarchy%is_valid() .and. size(coarse, 1) == size(fine, 1) .and. &
      size(integral) == size(coarse, 1) .and. &
      size(coarse, 2) == hierarchy%coarse%cell_count() .and. &
      size(fine, 2) == hierarchy%fine%cell_count()
    if (.not. ok) return

    do component = 1, size(integral)
      integral(component) = hierarchy%coarse_dx * ( &
        sum(coarse(component, 1:hierarchy%fine_coarse_lower - 1)) + &
        sum(coarse(component, hierarchy%fine_coarse_upper + 1: &
          hierarchy%coarse%upper))) + &
        hierarchy%fine_dx * sum(fine(component, :))
    end do
  end subroutine composite_integral_1d

  subroutine synchronize_multilevel_1d(fields, hierarchy, registers, ok)
    type(amr_level_field_1d), intent(inout) :: fields(:)
    type(amr_multilevel_hierarchy_1d), intent(in) :: hierarchy
    type(amr_flux_register_1d), intent(inout) :: registers(:)
    logical, intent(out) :: ok

    type(amr_level_field_1d), allocatable :: field_backup(:)
    type(amr_flux_register_1d), allocatable :: register_backup(:)
    logical :: local_ok
    integer :: relation, variable_count

    ok = valid_multilevel_fields(fields, hierarchy)
    if (.not. ok) return
    ok = size(registers) == size(hierarchy%interfaces)
    if (.not. ok) return
    variable_count = size(fields(1)%values, 1)
    do relation = 1, size(registers)
      if (.not. valid_flux_register(registers(relation), variable_count)) then
        ok = .false.
        return
      end if
    end do

    field_backup = fields
    register_backup = registers
    do relation = size(hierarchy%interfaces), 1, -1
      call reflux_1d( &
        fields(relation)%values, hierarchy%interfaces(relation), &
        registers(relation), local_ok)
      if (.not. local_ok) then
        fields = field_backup
        registers = register_backup
        ok = .false.
        return
      end if
      call average_down_1d( &
        fields(relation + 1)%values, hierarchy%interfaces(relation), &
        fields(relation)%values, local_ok)
      if (.not. local_ok) then
        fields = field_backup
        registers = register_backup
        ok = .false.
        return
      end if
    end do
    ok = .true.
  end subroutine synchronize_multilevel_1d

  pure subroutine composite_integral_multilevel_1d( &
      fields, hierarchy, integral, ok)
    type(amr_level_field_1d), intent(in) :: fields(:)
    type(amr_multilevel_hierarchy_1d), intent(in) :: hierarchy
    real(dp), intent(out) :: integral(:)
    logical, intent(out) :: ok

    integer :: relation, deepest, component

    integral = 0.0_dp
    ok = valid_multilevel_fields(fields, hierarchy)
    if (.not. ok) return
    ok = size(integral) == size(fields(1)%values, 1)
    if (.not. ok) return

    do relation = 1, size(hierarchy%interfaces)
      do component = 1, size(integral)
        integral(component) = integral(component) + &
          hierarchy%interfaces(relation)%coarse_dx * ( &
            sum(fields(relation)%values(component, &
              1:hierarchy%interfaces(relation)%fine_coarse_lower - 1)) + &
            sum(fields(relation)%values(component, &
              hierarchy%interfaces(relation)%fine_coarse_upper + 1: &
              hierarchy%interfaces(relation)%coarse%upper)))
      end do
    end do
    deepest = hierarchy%level_count()
    integral = integral + hierarchy%level_dx(deepest - 1) * &
      sum(fields(deepest)%values, dim=2)
    ok = .true.
  end subroutine composite_integral_multilevel_1d

  pure logical function valid_multilevel_fields(fields, hierarchy) &
      result(valid)
    type(amr_level_field_1d), intent(in) :: fields(:)
    type(amr_multilevel_hierarchy_1d), intent(in) :: hierarchy

    integer :: level, variable_count

    valid = hierarchy%is_valid() .and. &
      size(fields) == hierarchy%level_count()
    if (.not. valid) return
    valid = allocated(fields(1)%values)
    if (.not. valid) return
    variable_count = size(fields(1)%values, 1)
    valid = variable_count >= 1
    if (.not. valid) return
    do level = 1, size(fields)
      valid = allocated(fields(level)%values)
      if (.not. valid) return
      valid = size(fields(level)%values, 1) == variable_count .and. &
        size(fields(level)%values, 2) == &
          hierarchy%level_cell_count(level - 1)
      if (.not. valid) return
    end do
  end function valid_multilevel_fields

  pure logical function valid_flux_register(register, variable_count) &
      result(valid)
    type(amr_flux_register_1d), intent(in) :: register
    integer, intent(in) :: variable_count

    valid = allocated(register%left) .and. allocated(register%right)
    if (.not. valid) return
    valid = variable_count >= 1 .and. &
      size(register%left) == variable_count .and. &
      size(register%right) == variable_count
  end function valid_flux_register

end module amr_hierarchy_1d_mod
