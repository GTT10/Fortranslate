module amr_patch_tree_1d_mod
  use precision_mod, only: dp
  use amr_hierarchy_1d_mod, only: &
    amr_two_level_hierarchy_1d, amr_level_field_1d, amr_flux_register_1d
  use amr_multipatch_1d_mod, only: &
    amr_patch_set_1d, initialize_patch_set_1d, prolong_patch_set_1d, &
    average_down_patch_set_1d, initialize_patch_flux_registers_1d, &
    synchronize_patch_set_1d
  implicit none
  private

  type, public :: amr_child_patch_plan_1d
    integer :: parent_patch = 0
    integer :: lower = 1
    integer :: upper = 0
  end type amr_child_patch_plan_1d

  type, public :: amr_patch_level_plan_1d
    integer :: refinement_ratio = 0
    type(amr_child_patch_plan_1d), allocatable :: patches(:)
  contains
    procedure :: patch_count => patch_level_plan_patch_count
  end type amr_patch_level_plan_1d

  type, public :: amr_patch_tree_relation_1d
    integer :: level = 0
    integer :: refinement_ratio = 0
    type(amr_patch_set_1d), allocatable :: child_sets(:)
    integer, allocatable :: child_offsets(:)
  contains
    procedure :: parent_patch_count => patch_tree_parent_patch_count
    procedure :: child_patch_count => patch_tree_child_patch_count
    procedure :: child_index => patch_tree_child_index
    procedure :: is_valid => patch_tree_relation_is_valid
  end type amr_patch_tree_relation_1d

  type, public :: amr_patch_tree_hierarchy_1d
    integer :: base_cells = 0
    real(dp) :: x_lower = 0.0_dp
    real(dp) :: x_upper = 0.0_dp
    type(amr_patch_tree_relation_1d), allocatable :: relations(:)
  contains
    procedure :: level_count => patch_tree_level_count
    procedure :: level_patch_count => patch_tree_level_patch_count
    procedure :: level_dx => patch_tree_level_dx
    procedure :: is_valid => patch_tree_hierarchy_is_valid
  end type amr_patch_tree_hierarchy_1d

  type, public :: amr_patch_tree_level_fields_1d
    type(amr_level_field_1d), allocatable :: patches(:)
  end type amr_patch_tree_level_fields_1d

  type, public :: amr_patch_tree_parent_flux_registers_1d
    type(amr_flux_register_1d), allocatable :: children(:)
  end type amr_patch_tree_parent_flux_registers_1d

  type, public :: amr_patch_tree_relation_flux_registers_1d
    type(amr_patch_tree_parent_flux_registers_1d), allocatable :: parents(:)
  end type amr_patch_tree_relation_flux_registers_1d

  public :: initialize_patch_tree_1d
  public :: prolong_patch_tree_1d
  public :: average_down_patch_tree_1d
  public :: initialize_patch_tree_flux_registers_1d
  public :: synchronize_patch_tree_1d
  public :: composite_integral_patch_tree_1d
  public :: patch_tree_fields_are_valid_1d
  public :: patch_tree_flux_registers_are_valid_1d
  public :: patch_tree_child_geometry_1d

contains

  pure integer function patch_level_plan_patch_count(self) result(count)
    class(amr_patch_level_plan_1d), intent(in) :: self

    count = 0
    if (allocated(self%patches)) count = size(self%patches)
  end function patch_level_plan_patch_count

  pure integer function patch_tree_parent_patch_count(self) result(count)
    class(amr_patch_tree_relation_1d), intent(in) :: self

    count = 0
    if (allocated(self%child_sets)) count = size(self%child_sets)
  end function patch_tree_parent_patch_count

  pure integer function patch_tree_child_patch_count(self) result(count)
    class(amr_patch_tree_relation_1d), intent(in) :: self

    count = 0
    if (.not. allocated(self%child_offsets)) return
    if (size(self%child_offsets) < 1) return
    count = self%child_offsets(size(self%child_offsets))
  end function patch_tree_child_patch_count

  pure integer function patch_tree_child_index( &
      self, parent_patch, local_child) result(index)
    class(amr_patch_tree_relation_1d), intent(in) :: self
    integer, intent(in) :: parent_patch, local_child

    index = 0
    if (.not. self%is_valid()) return
    if (parent_patch < 1 .or. &
        parent_patch > self%parent_patch_count()) return
    if (local_child < 1 .or. &
        local_child > self%child_sets(parent_patch)%patch_count()) return
    index = self%child_offsets(parent_patch) + local_child
  end function patch_tree_child_index

  pure logical function patch_tree_relation_is_valid(self) result(valid)
    class(amr_patch_tree_relation_1d), intent(in) :: self

    integer :: parent

    valid = self%level >= 1 .and. self%refinement_ratio >= 2 .and. &
      allocated(self%child_sets) .and. allocated(self%child_offsets)
    if (.not. valid) return
    valid = size(self%child_sets) >= 1 .and. &
      size(self%child_offsets) == size(self%child_sets) + 1 .and. &
      self%child_offsets(1) == 0
    if (.not. valid) return
    do parent = 1, size(self%child_sets)
      valid = self%child_sets(parent)%is_valid() .and. &
        self%child_sets(parent)%coarse_level == self%level - 1 .and. &
        self%child_sets(parent)%refinement_ratio == self%refinement_ratio .and. &
        self%child_offsets(parent + 1) == &
          self%child_offsets(parent) + &
          self%child_sets(parent)%patch_count()
      if (.not. valid) return
    end do
    valid = self%child_patch_count() >= 1
  end function patch_tree_relation_is_valid

  pure integer function patch_tree_level_count(self) result(count)
    class(amr_patch_tree_hierarchy_1d), intent(in) :: self

    count = 0
    if (allocated(self%relations)) count = size(self%relations) + 1
  end function patch_tree_level_count

  pure integer function patch_tree_level_patch_count( &
      self, level) result(count)
    class(amr_patch_tree_hierarchy_1d), intent(in) :: self
    integer, intent(in) :: level

    count = 0
    if (level < 0 .or. level >= self%level_count()) return
    if (level == 0) then
      count = 1
    else
      count = self%relations(level)%child_patch_count()
    end if
  end function patch_tree_level_patch_count

  pure real(dp) function patch_tree_level_dx(self, level) result(dx)
    class(amr_patch_tree_hierarchy_1d), intent(in) :: self
    integer, intent(in) :: level

    integer :: relation

    dx = 0.0_dp
    if (level < 0 .or. level >= self%level_count()) return
    if (self%base_cells < 1 .or. self%x_upper <= self%x_lower) return
    dx = (self%x_upper - self%x_lower) / real(self%base_cells, dp)
    do relation = 1, level
      dx = dx / real(self%relations(relation)%refinement_ratio, dp)
    end do
  end function patch_tree_level_dx

  pure logical function patch_tree_hierarchy_is_valid(self) result(valid)
    class(amr_patch_tree_hierarchy_1d), intent(in) :: self

    real(dp) :: expected_lower, expected_upper, tolerance
    logical :: local_ok
    integer :: relation, parent, expected_cells, expected_parent_count

    valid = self%base_cells >= 1 .and. self%x_upper > self%x_lower .and. &
      allocated(self%relations)
    if (.not. valid) return
    expected_parent_count = 1
    do relation = 1, size(self%relations)
      valid = self%relations(relation)%is_valid() .and. &
        self%relations(relation)%level == relation .and. &
        self%relations(relation)%parent_patch_count() == &
          expected_parent_count
      if (.not. valid) return
      do parent = 1, expected_parent_count
        call patch_tree_parent_description( &
          self, relation, parent, expected_cells, expected_lower, &
          expected_upper, local_ok)
        if (.not. local_ok) then
          valid = .false.
          return
        end if
        tolerance = 64.0_dp * epsilon(1.0_dp) * &
          max(1.0_dp, abs(expected_lower), abs(expected_upper))
        valid = self%relations(relation)%child_sets(parent)%coarse_cells == &
          expected_cells .and. &
          abs(self%relations(relation)%child_sets(parent)%x_lower - &
            expected_lower) <= tolerance .and. &
          abs(self%relations(relation)%child_sets(parent)%x_upper - &
            expected_upper) <= tolerance
        if (.not. valid) return
      end do
      expected_parent_count = &
        self%relations(relation)%child_patch_count()
    end do
  end function patch_tree_hierarchy_is_valid

  subroutine initialize_patch_tree_1d( &
      base_cells, x_lower, x_upper, plans, hierarchy, ok)
    integer, intent(in) :: base_cells
    real(dp), intent(in) :: x_lower, x_upper
    type(amr_patch_level_plan_1d), intent(in) :: plans(:)
    type(amr_patch_tree_hierarchy_1d), intent(out) :: hierarchy
    logical, intent(out) :: ok

    integer, allocatable :: lower(:), upper(:)
    real(dp) :: parent_lower, parent_upper
    logical :: local_ok
    integer :: relation, parent, entry, child, child_count
    integer :: parent_count, parent_cells

    hierarchy%base_cells = base_cells
    hierarchy%x_lower = x_lower
    hierarchy%x_upper = x_upper
    ok = base_cells >= 1 .and. x_upper > x_lower
    if (.not. ok) return
    allocate(hierarchy%relations(size(plans)))
    if (size(plans) == 0) then
      ok = hierarchy%is_valid()
      return
    end if

    parent_count = 1
    do relation = 1, size(plans)
      ok = plans(relation)%refinement_ratio >= 2 .and. &
        allocated(plans(relation)%patches) .and. &
        plans(relation)%patch_count() >= 1 .and. parent_count >= 1
      if (.not. ok) return
      hierarchy%relations(relation)%level = relation
      hierarchy%relations(relation)%refinement_ratio = &
        plans(relation)%refinement_ratio
      allocate(hierarchy%relations(relation)%child_sets(parent_count))
      allocate(hierarchy%relations(relation)%child_offsets(parent_count + 1))
      hierarchy%relations(relation)%child_offsets(1) = 0

      do parent = 1, parent_count
        child_count = 0
        do entry = 1, plans(relation)%patch_count()
          if (plans(relation)%patches(entry)%parent_patch == parent) &
            child_count = child_count + 1
        end do
        allocate(lower(child_count), upper(child_count))
        child = 0
        do entry = 1, plans(relation)%patch_count()
          if (plans(relation)%patches(entry)%parent_patch /= parent) cycle
          child = child + 1
          lower(child) = plans(relation)%patches(entry)%lower
          upper(child) = plans(relation)%patches(entry)%upper
        end do
        call patch_tree_parent_description( &
          hierarchy, relation, parent, parent_cells, parent_lower, &
          parent_upper, local_ok)
        if (.not. local_ok) then
          ok = .false.
          return
        end if
        call initialize_patch_set_1d( &
          parent_cells, lower, upper, plans(relation)%refinement_ratio, &
          parent_lower, parent_upper, &
          hierarchy%relations(relation)%child_sets(parent), local_ok, &
          relation - 1)
        deallocate(lower, upper)
        if (.not. local_ok) then
          ok = .false.
          return
        end if
        hierarchy%relations(relation)%child_offsets(parent + 1) = &
          hierarchy%relations(relation)%child_offsets(parent) + child_count
      end do

      ok = all(plans(relation)%patches%parent_patch >= 1) .and. &
        all(plans(relation)%patches%parent_patch <= parent_count) .and. &
        hierarchy%relations(relation)%is_valid()
      if (.not. ok) return
      parent_count = hierarchy%relations(relation)%child_patch_count()
    end do
    ok = hierarchy%is_valid()
  end subroutine initialize_patch_tree_1d

  subroutine prolong_patch_tree_1d(root, hierarchy, fields, ok)
    real(dp), intent(in) :: root(:, :)
    type(amr_patch_tree_hierarchy_1d), intent(in) :: hierarchy
    type(amr_patch_tree_level_fields_1d), allocatable, intent(out) :: fields(:)
    logical, intent(out) :: ok

    type(amr_level_field_1d), allocatable :: children(:)
    logical :: local_ok
    integer :: relation, parent, child, index, variable_count

    ok = hierarchy%is_valid() .and. size(root, 1) >= 1 .and. &
      size(root, 2) == hierarchy%base_cells
    if (.not. ok) return
    variable_count = size(root, 1)
    allocate(fields(hierarchy%level_count()))
    allocate(fields(1)%patches(1))
    allocate(fields(1)%patches(1)%values(variable_count, hierarchy%base_cells))
    fields(1)%patches(1)%values = root

    do relation = 1, size(hierarchy%relations)
      allocate(fields(relation + 1)%patches( &
        hierarchy%relations(relation)%child_patch_count()))
      do parent = 1, &
          hierarchy%relations(relation)%parent_patch_count()
        call prolong_patch_set_1d( &
          fields(relation)%patches(parent)%values, &
          hierarchy%relations(relation)%child_sets(parent), &
          children, local_ok)
        if (.not. local_ok) then
          ok = .false.
          return
        end if
        do child = 1, size(children)
          index = hierarchy%relations(relation)%child_index(parent, child)
          fields(relation + 1)%patches(index)%values = &
            children(child)%values
        end do
      end do
    end do
    ok = patch_tree_fields_are_valid_1d(fields, hierarchy)
  end subroutine prolong_patch_tree_1d

  subroutine average_down_patch_tree_1d(fields, hierarchy, ok)
    type(amr_patch_tree_level_fields_1d), intent(inout) :: fields(:)
    type(amr_patch_tree_hierarchy_1d), intent(in) :: hierarchy
    logical, intent(out) :: ok

    type(amr_patch_tree_level_fields_1d), allocatable :: backup(:)
    logical :: local_ok
    integer :: relation, parent, first_child, last_child

    ok = patch_tree_fields_are_valid_1d(fields, hierarchy)
    if (.not. ok) return
    backup = fields
    do relation = size(hierarchy%relations), 1, -1
      do parent = 1, &
          hierarchy%relations(relation)%parent_patch_count()
        first_child = &
          hierarchy%relations(relation)%child_offsets(parent) + 1
        last_child = &
          hierarchy%relations(relation)%child_offsets(parent + 1)
        if (last_child < first_child) cycle
        call average_down_patch_set_1d( &
          fields(relation + 1)%patches(first_child:last_child), &
          hierarchy%relations(relation)%child_sets(parent), &
          fields(relation)%patches(parent)%values, local_ok)
        if (.not. local_ok) then
          fields = backup
          ok = .false.
          return
        end if
      end do
    end do
    ok = patch_tree_fields_are_valid_1d(fields, hierarchy)
  end subroutine average_down_patch_tree_1d

  subroutine initialize_patch_tree_flux_registers_1d( &
      hierarchy, variable_count, registers, ok)
    type(amr_patch_tree_hierarchy_1d), intent(in) :: hierarchy
    integer, intent(in) :: variable_count
    type(amr_patch_tree_relation_flux_registers_1d), allocatable, &
      intent(out) :: registers(:)
    logical, intent(out) :: ok

    logical :: local_ok
    integer :: relation, parent

    ok = hierarchy%is_valid() .and. variable_count >= 1
    if (.not. ok) return
    allocate(registers(size(hierarchy%relations)))
    do relation = 1, size(hierarchy%relations)
      allocate(registers(relation)%parents( &
        hierarchy%relations(relation)%parent_patch_count()))
      do parent = 1, hierarchy%relations(relation)%parent_patch_count()
        call initialize_patch_flux_registers_1d( &
          hierarchy%relations(relation)%child_sets(parent), &
          variable_count, registers(relation)%parents(parent)%children, &
          local_ok)
        if (.not. local_ok) then
          ok = .false.
          return
        end if
      end do
    end do
    ok = patch_tree_flux_registers_are_valid_1d( &
      registers, hierarchy, variable_count)
  end subroutine initialize_patch_tree_flux_registers_1d

  subroutine synchronize_patch_tree_1d( &
      fields, hierarchy, registers, ok)
    type(amr_patch_tree_level_fields_1d), intent(inout) :: fields(:)
    type(amr_patch_tree_hierarchy_1d), intent(in) :: hierarchy
    type(amr_patch_tree_relation_flux_registers_1d), &
      intent(inout) :: registers(:)
    logical, intent(out) :: ok

    type(amr_patch_tree_level_fields_1d), allocatable :: field_backup(:)
    type(amr_patch_tree_relation_flux_registers_1d), allocatable :: &
      register_backup(:)
    logical :: local_ok
    integer :: relation, parent, first_child, last_child, variable_count

    ok = patch_tree_fields_are_valid_1d(fields, hierarchy)
    if (.not. ok) return
    variable_count = size(fields(1)%patches(1)%values, 1)
    ok = patch_tree_flux_registers_are_valid_1d( &
      registers, hierarchy, variable_count)
    if (.not. ok) return
    field_backup = fields
    register_backup = registers
    do relation = size(hierarchy%relations), 1, -1
      do parent = 1, hierarchy%relations(relation)%parent_patch_count()
        first_child = &
          hierarchy%relations(relation)%child_offsets(parent) + 1
        last_child = &
          hierarchy%relations(relation)%child_offsets(parent + 1)
        if (last_child < first_child) cycle
        call synchronize_patch_set_1d( &
          fields(relation)%patches(parent)%values, &
          fields(relation + 1)%patches(first_child:last_child), &
          hierarchy%relations(relation)%child_sets(parent), &
          registers(relation)%parents(parent)%children, local_ok)
        if (.not. local_ok) then
          fields = field_backup
          registers = register_backup
          ok = .false.
          return
        end if
      end do
    end do
    ok = patch_tree_fields_are_valid_1d(fields, hierarchy) .and. &
      patch_tree_flux_registers_are_valid_1d( &
        registers, hierarchy, variable_count)
  end subroutine synchronize_patch_tree_1d

  pure subroutine composite_integral_patch_tree_1d( &
      fields, hierarchy, integral, ok)
    type(amr_patch_tree_level_fields_1d), intent(in) :: fields(:)
    type(amr_patch_tree_hierarchy_1d), intent(in) :: hierarchy
    real(dp), intent(out) :: integral(:)
    logical, intent(out) :: ok

    type(amr_two_level_hierarchy_1d) :: geometry
    logical :: local_ok
    integer :: relation, parent, child, index, component, lower, upper

    integral = 0.0_dp
    ok = patch_tree_fields_are_valid_1d(fields, hierarchy) .and. &
      size(integral) == size(fields(1)%patches(1)%values, 1)
    if (.not. ok) return
    integral = hierarchy%level_dx(0) * &
      sum(fields(1)%patches(1)%values, dim=2)
    do relation = 1, size(hierarchy%relations)
      do parent = 1, &
          hierarchy%relations(relation)%parent_patch_count()
        do child = 1, hierarchy%relations(relation)% &
            child_sets(parent)%patch_count()
          index = hierarchy%relations(relation)%child_index(parent, child)
          call patch_tree_child_geometry_1d( &
            hierarchy%relations(relation), index, geometry, local_ok)
          if (.not. local_ok) then
            ok = .false.
            return
          end if
          lower = geometry%fine_coarse_lower
          upper = geometry%fine_coarse_upper
          do component = 1, size(integral)
            integral(component) = integral(component) - &
              geometry%coarse_dx * sum( &
                fields(relation)%patches(parent)%values( &
                  component, lower:upper)) + &
              geometry%fine_dx * sum( &
                fields(relation + 1)%patches(index)%values(component, :))
          end do
        end do
      end do
    end do
    ok = .true.
  end subroutine composite_integral_patch_tree_1d

  pure logical function patch_tree_fields_are_valid_1d( &
      fields, hierarchy) result(valid)
    type(amr_patch_tree_level_fields_1d), intent(in) :: fields(:)
    type(amr_patch_tree_hierarchy_1d), intent(in) :: hierarchy

    type(amr_two_level_hierarchy_1d) :: geometry
    logical :: local_ok
    integer :: relation, patch, variable_count

    valid = hierarchy%is_valid() .and. &
      size(fields) == hierarchy%level_count()
    if (.not. valid) return
    valid = allocated(fields(1)%patches) .and. &
      size(fields(1)%patches) == 1 .and. &
      allocated(fields(1)%patches(1)%values)
    if (.not. valid) return
    variable_count = size(fields(1)%patches(1)%values, 1)
    valid = variable_count >= 1 .and. &
      size(fields(1)%patches(1)%values, 2) == hierarchy%base_cells
    if (.not. valid) return

    do relation = 1, size(hierarchy%relations)
      valid = allocated(fields(relation + 1)%patches) .and. &
        size(fields(relation + 1)%patches) == &
          hierarchy%relations(relation)%child_patch_count()
      if (.not. valid) return
      do patch = 1, size(fields(relation + 1)%patches)
        call patch_tree_child_geometry_1d( &
          hierarchy%relations(relation), patch, geometry, local_ok)
        if (.not. local_ok) then
          valid = .false.
          return
        end if
        valid = allocated(fields(relation + 1)%patches(patch)%values)
        if (.not. valid) return
        valid = size(fields(relation + 1)%patches(patch)%values, 1) == &
          variable_count .and. &
          size(fields(relation + 1)%patches(patch)%values, 2) == &
            geometry%fine%cell_count()
        if (.not. valid) return
      end do
    end do
  end function patch_tree_fields_are_valid_1d

  pure logical function patch_tree_flux_registers_are_valid_1d( &
      registers, hierarchy, variable_count) result(valid)
    type(amr_patch_tree_relation_flux_registers_1d), &
      intent(in) :: registers(:)
    type(amr_patch_tree_hierarchy_1d), intent(in) :: hierarchy
    integer, intent(in) :: variable_count

    integer :: relation, parent, child

    valid = hierarchy%is_valid() .and. variable_count >= 1 .and. &
      size(registers) == size(hierarchy%relations)
    if (.not. valid) return
    do relation = 1, size(registers)
      valid = allocated(registers(relation)%parents) .and. &
        size(registers(relation)%parents) == &
          hierarchy%relations(relation)%parent_patch_count()
      if (.not. valid) return
      do parent = 1, size(registers(relation)%parents)
        valid = allocated(registers(relation)%parents(parent)%children) .and. &
          size(registers(relation)%parents(parent)%children) == &
            hierarchy%relations(relation)%child_sets(parent)%patch_count()
        if (.not. valid) return
        do child = 1, &
            size(registers(relation)%parents(parent)%children)
          valid = allocated( &
            registers(relation)%parents(parent)%children(child)%left) .and. &
            allocated( &
            registers(relation)%parents(parent)%children(child)%right)
          if (.not. valid) return
          valid = size( &
            registers(relation)%parents(parent)%children(child)%left) == &
              variable_count .and. &
            size( &
            registers(relation)%parents(parent)%children(child)%right) == &
              variable_count
          if (.not. valid) return
        end do
      end do
    end do
  end function patch_tree_flux_registers_are_valid_1d

  pure subroutine patch_tree_child_geometry_1d( &
      relation, child_index, geometry, ok)
    type(amr_patch_tree_relation_1d), intent(in) :: relation
    integer, intent(in) :: child_index
    type(amr_two_level_hierarchy_1d), intent(out) :: geometry
    logical, intent(out) :: ok

    integer :: parent, local_child

    ok = relation%is_valid() .and. child_index >= 1 .and. &
      child_index <= relation%child_patch_count()
    if (.not. ok) return
    do parent = 1, relation%parent_patch_count()
      if (child_index <= relation%child_offsets(parent) .or. &
          child_index > relation%child_offsets(parent + 1)) cycle
      local_child = child_index - relation%child_offsets(parent)
      geometry = relation%child_sets(parent)%patches(local_child)
      ok = geometry%is_valid()
      return
    end do
    ok = .false.
  end subroutine patch_tree_child_geometry_1d

  pure subroutine patch_tree_parent_description( &
      hierarchy, relation, parent_patch, cells, x_lower, x_upper, ok)
    type(amr_patch_tree_hierarchy_1d), intent(in) :: hierarchy
    integer, intent(in) :: relation, parent_patch
    integer, intent(out) :: cells
    real(dp), intent(out) :: x_lower, x_upper
    logical, intent(out) :: ok

    type(amr_two_level_hierarchy_1d) :: geometry

    cells = 0
    x_lower = 0.0_dp
    x_upper = 0.0_dp
    ok = relation >= 1 .and. relation <= size(hierarchy%relations)
    if (.not. ok) return
    if (relation == 1) then
      ok = parent_patch == 1 .and. hierarchy%base_cells >= 1 .and. &
        hierarchy%x_upper > hierarchy%x_lower
      if (.not. ok) return
      cells = hierarchy%base_cells
      x_lower = hierarchy%x_lower
      x_upper = hierarchy%x_upper
      return
    end if
    call patch_tree_child_geometry_1d( &
      hierarchy%relations(relation - 1), parent_patch, geometry, ok)
    if (.not. ok) return
    cells = geometry%fine%cell_count()
    x_lower = geometry%x_lower + &
      real(geometry%fine_coarse_lower - 1, dp) * geometry%coarse_dx
    x_upper = geometry%x_lower + &
      real(geometry%fine_coarse_upper, dp) * geometry%coarse_dx
    ok = cells >= 1 .and. x_upper > x_lower
  end subroutine patch_tree_parent_description

end module amr_patch_tree_1d_mod
