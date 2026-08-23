module amr_regrid_1d_mod
  use precision_mod, only: dp
  use amr_hierarchy_1d_mod, only: &
    amr_two_level_hierarchy_1d, amr_level_field_1d, &
    initialize_two_level_hierarchy_1d, prolong_conservative_1d, &
    average_down_1d
  use amr_multipatch_1d_mod, only: &
    amr_patch_set_1d, initialize_patch_set_1d, &
    prolong_patch_set_1d, average_down_patch_set_1d, &
    patch_fields_are_valid_1d
  implicit none
  private

  type, public :: amr_tagging_criteria_1d
    integer :: component = 1
    real(dp) :: relative_gradient_threshold = 0.1_dp
    real(dp) :: absolute_gradient_threshold = 0.0_dp
    real(dp) :: scale_floor = 1.0e-12_dp
    integer :: buffer_cells = 1
    integer :: minimum_patch_cells = 1
    integer :: maximum_patch_gap_cells = 0
  contains
    procedure :: is_valid => tagging_criteria_is_valid
  end type amr_tagging_criteria_1d

  type, public :: amr_regrid_plan_1d
    logical :: active = .false.
    integer :: coarse_cells = 0
    integer :: tagged_cell_count = 0
    integer :: tag_lower = 1
    integer :: tag_upper = 0
    integer :: patch_lower = 1
    integer :: patch_upper = 0
  contains
    procedure :: is_valid => regrid_plan_is_valid
  end type amr_regrid_plan_1d

  type, public :: amr_regrid_plan_collection_1d
    integer :: coarse_cells = 0
    integer :: tagged_cell_count = 0
    type(amr_regrid_plan_1d), allocatable :: plans(:)
  contains
    procedure :: patch_count => regrid_plan_collection_patch_count
    procedure :: is_valid => regrid_plan_collection_is_valid
  end type amr_regrid_plan_collection_1d

  public :: tag_gradient_1d
  public :: build_regrid_plan_1d
  public :: build_regrid_plan_collection_1d
  public :: plan_gradient_regrid_1d
  public :: plan_gradient_regrid_collection_1d
  public :: regrid_two_level_state_1d
  public :: regrid_patch_set_state_1d

contains

  pure logical function tagging_criteria_is_valid(self, variable_count) &
      result(valid)
    class(amr_tagging_criteria_1d), intent(in) :: self
    integer, intent(in), optional :: variable_count

    valid = self%component >= 1 .and. &
      self%relative_gradient_threshold >= 0.0_dp .and. &
      self%absolute_gradient_threshold >= 0.0_dp .and. &
      self%scale_floor > 0.0_dp .and. self%buffer_cells >= 0 .and. &
      self%minimum_patch_cells >= 1 .and. &
      self%maximum_patch_gap_cells >= 0
    if (present(variable_count)) then
      valid = valid .and. self%component <= variable_count
    end if
  end function tagging_criteria_is_valid

  pure logical function regrid_plan_is_valid(self) result(valid)
    class(amr_regrid_plan_1d), intent(in) :: self

    valid = self%coarse_cells >= 3 .and. self%tagged_cell_count >= 0
    if (.not. valid) return
    if (.not. self%active) then
      valid = self%tagged_cell_count == 0 .and. &
        self%tag_lower > self%tag_upper .and. &
        self%patch_lower > self%patch_upper
      return
    end if
    valid = self%tagged_cell_count >= 1 .and. &
      self%tag_lower >= 1 .and. &
      self%tag_upper <= self%coarse_cells .and. &
      self%tag_lower <= self%tag_upper .and. &
      self%patch_lower >= 1 .and. &
      self%patch_upper <= self%coarse_cells .and. &
      self%patch_lower <= self%tag_lower .and. &
      self%patch_upper >= self%tag_upper
  end function regrid_plan_is_valid

  pure integer function regrid_plan_collection_patch_count(self) &
      result(count)
    class(amr_regrid_plan_collection_1d), intent(in) :: self

    count = 0
    if (allocated(self%plans)) count = size(self%plans)
  end function regrid_plan_collection_patch_count

  pure logical function regrid_plan_collection_is_valid(self) result(valid)
    class(amr_regrid_plan_collection_1d), intent(in) :: self

    integer :: patch, previous_upper, tagged_count

    valid = self%coarse_cells >= 3 .and. &
      self%tagged_cell_count >= 0 .and. allocated(self%plans)
    if (.not. valid) return
    if (size(self%plans) == 0) then
      valid = self%tagged_cell_count == 0
      return
    end if

    previous_upper = -1
    tagged_count = 0
    do patch = 1, size(self%plans)
      valid = self%plans(patch)%is_valid() .and. &
        self%plans(patch)%active .and. &
        self%plans(patch)%coarse_cells == self%coarse_cells
      if (.not. valid) return
      if (patch > 1) then
        valid = self%plans(patch)%patch_lower > previous_upper + 1
        if (.not. valid) return
      end if
      previous_upper = self%plans(patch)%patch_upper
      tagged_count = tagged_count + self%plans(patch)%tagged_cell_count
    end do
    valid = tagged_count == self%tagged_cell_count
  end function regrid_plan_collection_is_valid

  pure subroutine tag_gradient_1d(state, criteria, tags, ok)
    real(dp), intent(in) :: state(:, :)
    type(amr_tagging_criteria_1d), intent(in) :: criteria
    logical, intent(out) :: tags(:)
    logical, intent(out) :: ok

    real(dp) :: jump, local_scale
    integer :: cell, left_cell, right_cell

    tags = .false.
    ok = size(state, 1) >= 1 .and. size(state, 2) >= 3 .and. &
      size(tags) == size(state, 2) .and. &
      criteria%is_valid(size(state, 1))
    if (.not. ok) return

    do cell = 1, size(state, 2)
      left_cell = max(1, cell - 1)
      right_cell = min(size(state, 2), cell + 1)
      jump = max( &
        abs(state(criteria%component, cell) - &
          state(criteria%component, left_cell)), &
        abs(state(criteria%component, right_cell) - &
          state(criteria%component, cell)))
      local_scale = max( &
        criteria%scale_floor, &
        abs(state(criteria%component, left_cell)), &
        abs(state(criteria%component, cell)), &
        abs(state(criteria%component, right_cell)))
      tags(cell) = jump > criteria%absolute_gradient_threshold .and. &
        jump / local_scale >= criteria%relative_gradient_threshold
    end do
  end subroutine tag_gradient_1d

  pure subroutine build_regrid_plan_1d( &
      tags, buffer_cells, minimum_patch_cells, plan, ok)
    logical, intent(in) :: tags(:)
    integer, intent(in) :: buffer_cells, minimum_patch_cells
    type(amr_regrid_plan_1d), intent(out) :: plan
    logical, intent(out) :: ok

    integer :: cell, first_tag, last_tag

    plan = amr_regrid_plan_1d()
    plan%coarse_cells = size(tags)
    ok = size(tags) >= 3 .and. buffer_cells >= 0 .and. &
      minimum_patch_cells >= 1 .and. &
      minimum_patch_cells <= size(tags)
    if (.not. ok) return

    plan%tagged_cell_count = count(tags)
    if (plan%tagged_cell_count == 0) then
      ok = plan%is_valid()
      return
    end if

    first_tag = 0
    last_tag = 0
    do cell = 1, size(tags)
      if (tags(cell)) then
        if (first_tag == 0) first_tag = cell
        last_tag = cell
      end if
    end do

    plan%active = .true.
    plan%tag_lower = first_tag
    plan%tag_upper = last_tag
    plan%patch_lower = max(1, first_tag - buffer_cells)
    plan%patch_upper = min(size(tags), last_tag + buffer_cells)
    do while (plan%patch_upper - plan%patch_lower + 1 < &
        minimum_patch_cells)
      if (plan%patch_lower > 1) plan%patch_lower = plan%patch_lower - 1
      if (plan%patch_upper - plan%patch_lower + 1 >= &
          minimum_patch_cells) exit
      if (plan%patch_upper < size(tags)) then
        plan%patch_upper = plan%patch_upper + 1
      end if
    end do
    ok = plan%is_valid() .and. &
      plan%patch_upper - plan%patch_lower + 1 >= minimum_patch_cells
  end subroutine build_regrid_plan_1d

  pure subroutine build_regrid_plan_collection_1d( &
      tags, buffer_cells, minimum_patch_cells, maximum_gap_cells, &
      collection, ok)
    logical, intent(in) :: tags(:)
    integer, intent(in) :: buffer_cells, minimum_patch_cells
    integer, intent(in) :: maximum_gap_cells
    type(amr_regrid_plan_collection_1d), intent(out) :: collection
    logical, intent(out) :: ok

    type(amr_regrid_plan_1d), allocatable :: candidates(:), merged(:)
    integer :: candidate_count, merged_count
    integer :: cell, first_tag, last_tag, cluster_tag_count

    collection = amr_regrid_plan_collection_1d()
    collection%coarse_cells = size(tags)
    collection%tagged_cell_count = count(tags)
    ok = size(tags) >= 3 .and. buffer_cells >= 0 .and. &
      minimum_patch_cells >= 1 .and. &
      minimum_patch_cells <= size(tags) .and. maximum_gap_cells >= 0
    if (.not. ok) return
    if (collection%tagged_cell_count == 0) then
      allocate(collection%plans(0))
      ok = collection%is_valid()
      return
    end if

    allocate(candidates(collection%tagged_cell_count))
    candidate_count = 0
    first_tag = 0
    last_tag = 0
    cluster_tag_count = 0
    do cell = 1, size(tags)
      if (.not. tags(cell)) cycle
      if (first_tag == 0) then
        first_tag = cell
        last_tag = cell
        cluster_tag_count = 1
      else if (cell - last_tag - 1 <= maximum_gap_cells) then
        last_tag = cell
        cluster_tag_count = cluster_tag_count + 1
      else
        candidate_count = candidate_count + 1
        call initialize_cluster_plan( &
          size(tags), first_tag, last_tag, cluster_tag_count, &
          buffer_cells, minimum_patch_cells, &
          candidates(candidate_count), ok)
        if (.not. ok) return
        first_tag = cell
        last_tag = cell
        cluster_tag_count = 1
      end if
    end do
    candidate_count = candidate_count + 1
    call initialize_cluster_plan( &
      size(tags), first_tag, last_tag, cluster_tag_count, &
      buffer_cells, minimum_patch_cells, candidates(candidate_count), ok)
    if (.not. ok) return

    allocate(merged(candidate_count))
    merged_count = 0
    do cell = 1, candidate_count
      if (merged_count == 0) then
        merged_count = merged_count + 1
        merged(merged_count) = candidates(cell)
      else if (candidates(cell)%patch_lower > &
          merged(merged_count)%patch_upper + 1) then
        merged_count = merged_count + 1
        merged(merged_count) = candidates(cell)
      else
        merged(merged_count)%tagged_cell_count = &
          merged(merged_count)%tagged_cell_count + &
          candidates(cell)%tagged_cell_count
        merged(merged_count)%tag_lower = min( &
          merged(merged_count)%tag_lower, candidates(cell)%tag_lower)
        merged(merged_count)%tag_upper = max( &
          merged(merged_count)%tag_upper, candidates(cell)%tag_upper)
        merged(merged_count)%patch_lower = min( &
          merged(merged_count)%patch_lower, candidates(cell)%patch_lower)
        merged(merged_count)%patch_upper = max( &
          merged(merged_count)%patch_upper, candidates(cell)%patch_upper)
      end if
    end do
    allocate(collection%plans(merged_count))
    collection%plans = merged(1:merged_count)
    ok = collection%is_valid()
  end subroutine build_regrid_plan_collection_1d

  pure subroutine plan_gradient_regrid_1d( &
      state, criteria, tags, plan, ok)
    real(dp), intent(in) :: state(:, :)
    type(amr_tagging_criteria_1d), intent(in) :: criteria
    logical, intent(out) :: tags(:)
    type(amr_regrid_plan_1d), intent(out) :: plan
    logical, intent(out) :: ok

    call tag_gradient_1d(state, criteria, tags, ok)
    if (.not. ok) return
    call build_regrid_plan_1d( &
      tags, criteria%buffer_cells, criteria%minimum_patch_cells, plan, ok)
  end subroutine plan_gradient_regrid_1d

  pure subroutine plan_gradient_regrid_collection_1d( &
      state, criteria, tags, collection, ok)
    real(dp), intent(in) :: state(:, :)
    type(amr_tagging_criteria_1d), intent(in) :: criteria
    logical, intent(out) :: tags(:)
    type(amr_regrid_plan_collection_1d), intent(out) :: collection
    logical, intent(out) :: ok

    call tag_gradient_1d(state, criteria, tags, ok)
    if (.not. ok) return
    call build_regrid_plan_collection_1d( &
      tags, criteria%buffer_cells, criteria%minimum_patch_cells, &
      criteria%maximum_patch_gap_cells, collection, ok)
  end subroutine plan_gradient_regrid_collection_1d

  pure subroutine initialize_cluster_plan( &
      coarse_cells, first_tag, last_tag, tagged_cell_count, &
      buffer_cells, minimum_patch_cells, plan, ok)
    integer, intent(in) :: coarse_cells, first_tag, last_tag
    integer, intent(in) :: tagged_cell_count
    integer, intent(in) :: buffer_cells, minimum_patch_cells
    type(amr_regrid_plan_1d), intent(out) :: plan
    logical, intent(out) :: ok

    plan = amr_regrid_plan_1d()
    plan%active = .true.
    plan%coarse_cells = coarse_cells
    plan%tagged_cell_count = tagged_cell_count
    plan%tag_lower = first_tag
    plan%tag_upper = last_tag
    plan%patch_lower = max(1, first_tag - buffer_cells)
    plan%patch_upper = min(coarse_cells, last_tag + buffer_cells)
    do while (plan%patch_upper - plan%patch_lower + 1 < &
        minimum_patch_cells)
      if (plan%patch_lower > 1) plan%patch_lower = plan%patch_lower - 1
      if (plan%patch_upper - plan%patch_lower + 1 >= &
          minimum_patch_cells) exit
      if (plan%patch_upper < coarse_cells) &
        plan%patch_upper = plan%patch_upper + 1
    end do
    ok = plan%is_valid() .and. plan%patch_upper - plan%patch_lower + 1 >= &
      minimum_patch_cells
  end subroutine initialize_cluster_plan

  subroutine regrid_patch_set_state_1d( &
      coarse, old_patch_set, old_fine_fields, collection, &
      refinement_ratio, x_lower, x_upper, &
      new_patch_set, new_fine_fields, ok, coarse_level)
    real(dp), intent(inout) :: coarse(:, :)
    type(amr_patch_set_1d), intent(in) :: old_patch_set
    type(amr_level_field_1d), intent(in) :: old_fine_fields(:)
    type(amr_regrid_plan_collection_1d), intent(in) :: collection
    integer, intent(in) :: refinement_ratio
    real(dp), intent(in) :: x_lower, x_upper
    type(amr_patch_set_1d), intent(out) :: new_patch_set
    type(amr_level_field_1d), allocatable, intent(out) :: new_fine_fields(:)
    logical, intent(out) :: ok
    integer, intent(in), optional :: coarse_level

    real(dp), allocatable :: coarse_backup(:, :)
    real(dp) :: domain_tolerance
    integer, allocatable :: patch_lower(:), patch_upper(:)
    integer :: parent_level, old_patch, new_patch
    integer :: overlap_lower, overlap_upper
    integer :: old_first, old_last, new_first, new_last
    logical :: retain_overlap

    parent_level = old_patch_set%coarse_level
    if (present(coarse_level)) parent_level = coarse_level
    ok = collection%is_valid() .and. old_patch_set%is_valid() .and. &
      collection%coarse_cells == size(coarse, 2) .and. &
      old_patch_set%coarse_cells == size(coarse, 2) .and. &
      size(coarse, 1) >= 1 .and. &
      patch_fields_are_valid_1d( &
        old_fine_fields, old_patch_set, size(coarse, 1))
    if (.not. ok) return

    coarse_backup = coarse
    call average_down_patch_set_1d( &
      old_fine_fields, old_patch_set, coarse, ok)
    if (.not. ok) then
      coarse = coarse_backup
      return
    end if

    allocate(patch_lower(collection%patch_count()))
    allocate(patch_upper(collection%patch_count()))
    do new_patch = 1, collection%patch_count()
      patch_lower(new_patch) = collection%plans(new_patch)%patch_lower
      patch_upper(new_patch) = collection%plans(new_patch)%patch_upper
    end do
    call initialize_patch_set_1d( &
      size(coarse, 2), patch_lower, patch_upper, refinement_ratio, &
      x_lower, x_upper, new_patch_set, ok, parent_level)
    if (.not. ok) then
      coarse = coarse_backup
      return
    end if
    call prolong_patch_set_1d( &
      coarse, new_patch_set, new_fine_fields, ok)
    if (.not. ok) then
      coarse = coarse_backup
      return
    end if

    domain_tolerance = 64.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, abs(x_lower), abs(x_upper))
    retain_overlap = &
      old_patch_set%refinement_ratio == refinement_ratio .and. &
      old_patch_set%coarse_level == parent_level .and. &
      abs(old_patch_set%x_lower - x_lower) <= domain_tolerance .and. &
      abs(old_patch_set%x_upper - x_upper) <= domain_tolerance
    if (retain_overlap) then
      do new_patch = 1, new_patch_set%patch_count()
        do old_patch = 1, old_patch_set%patch_count()
          overlap_lower = max( &
            new_patch_set%patches(new_patch)%fine_coarse_lower, &
            old_patch_set%patches(old_patch)%fine_coarse_lower)
          overlap_upper = min( &
            new_patch_set%patches(new_patch)%fine_coarse_upper, &
            old_patch_set%patches(old_patch)%fine_coarse_upper)
          if (overlap_lower > overlap_upper) cycle
          old_first = (overlap_lower - &
            old_patch_set%patches(old_patch)%fine_coarse_lower) * &
            refinement_ratio + 1
          old_last = (overlap_upper - &
            old_patch_set%patches(old_patch)%fine_coarse_lower + 1) * &
            refinement_ratio
          new_first = (overlap_lower - &
            new_patch_set%patches(new_patch)%fine_coarse_lower) * &
            refinement_ratio + 1
          new_last = (overlap_upper - &
            new_patch_set%patches(new_patch)%fine_coarse_lower + 1) * &
            refinement_ratio
          new_fine_fields(new_patch)%values(:, new_first:new_last) = &
            old_fine_fields(old_patch)%values(:, old_first:old_last)
        end do
      end do
    end if
    ok = patch_fields_are_valid_1d( &
      new_fine_fields, new_patch_set, size(coarse, 1))
    if (.not. ok) coarse = coarse_backup
  end subroutine regrid_patch_set_state_1d

  subroutine regrid_two_level_state_1d( &
      coarse, old_hierarchy, old_fine, plan, refinement_ratio, &
      x_lower, x_upper, new_hierarchy, new_fine, ok)
    real(dp), intent(inout) :: coarse(:, :)
    type(amr_two_level_hierarchy_1d), intent(in) :: old_hierarchy
    real(dp), allocatable, intent(in) :: old_fine(:, :)
    type(amr_regrid_plan_1d), intent(in) :: plan
    integer, intent(in) :: refinement_ratio
    real(dp), intent(in) :: x_lower, x_upper
    type(amr_two_level_hierarchy_1d), intent(out) :: new_hierarchy
    real(dp), allocatable, intent(out) :: new_fine(:, :)
    logical, intent(out) :: ok

    real(dp) :: domain_tolerance
    logical :: old_active
    integer :: overlap_lower, overlap_upper
    integer :: old_first, old_last, new_first, new_last

    new_hierarchy = amr_two_level_hierarchy_1d()
    old_active = old_hierarchy%is_valid()
    ok = size(coarse, 1) >= 1 .and. size(coarse, 2) >= 3 .and. &
      plan%is_valid() .and. plan%coarse_cells == size(coarse, 2) .and. &
      refinement_ratio >= 2 .and. x_upper > x_lower
    if (.not. ok) return

    if (old_active) then
      domain_tolerance = 32.0_dp * epsilon(1.0_dp) * &
        max(1.0_dp, abs(x_lower), abs(x_upper))
      ok = allocated(old_fine) .and. &
        old_hierarchy%coarse%cell_count() == size(coarse, 2) .and. &
        abs(old_hierarchy%x_lower - x_lower) <= domain_tolerance .and. &
        abs(old_hierarchy%x_upper - x_upper) <= domain_tolerance
      if (ok) then
        ok = size(old_fine, 1) == size(coarse, 1) .and. &
          size(old_fine, 2) == old_hierarchy%fine%cell_count()
      end if
    else
      ok = .not. allocated(old_fine)
    end if
    if (.not. ok) return

    if (old_active) then
      call average_down_1d(old_fine, old_hierarchy, coarse, ok)
      if (.not. ok) return
    end if
    if (.not. plan%active) then
      ok = .true.
      return
    end if

    call initialize_two_level_hierarchy_1d( &
      size(coarse, 2), plan%patch_lower, plan%patch_upper, &
      refinement_ratio, x_lower, x_upper, new_hierarchy, ok)
    if (.not. ok) return
    allocate(new_fine(size(coarse, 1), new_hierarchy%fine%cell_count()))
    call prolong_conservative_1d(coarse, new_hierarchy, new_fine, ok)
    if (.not. ok) return

    if (old_active .and. &
        old_hierarchy%refinement_ratio == refinement_ratio) then
      overlap_lower = max( &
        old_hierarchy%fine_coarse_lower, &
        new_hierarchy%fine_coarse_lower)
      overlap_upper = min( &
        old_hierarchy%fine_coarse_upper, &
        new_hierarchy%fine_coarse_upper)
      if (overlap_lower <= overlap_upper) then
        old_first = &
          (overlap_lower - old_hierarchy%fine_coarse_lower) * &
          refinement_ratio + 1
        old_last = &
          (overlap_upper - old_hierarchy%fine_coarse_lower + 1) * &
          refinement_ratio
        new_first = &
          (overlap_lower - new_hierarchy%fine_coarse_lower) * &
          refinement_ratio + 1
        new_last = &
          (overlap_upper - new_hierarchy%fine_coarse_lower + 1) * &
          refinement_ratio
        new_fine(:, new_first:new_last) = old_fine(:, old_first:old_last)
      end if
    end if
  end subroutine regrid_two_level_state_1d

end module amr_regrid_1d_mod
