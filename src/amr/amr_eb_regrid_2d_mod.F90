module amr_eb_regrid_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_covered_cell
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, average_down_reactive_eb_state_patch_2d
  use amr_eb_reactive_2d_mod, only: prolong_reactive_eb_patch_pcm_2d
  implicit none
  private

  type, public :: amr_eb_tagging_criteria_2d
    real(dp) :: relative_gradient_threshold = 0.1_dp
    real(dp) :: absolute_gradient_threshold = 0.0_dp
    real(dp) :: scale_floor = 1.0e-12_dp
    integer :: buffer_cells = 1
    integer :: minimum_patch_cells_x = 2
    integer :: minimum_patch_cells_y = 2
  contains
    procedure :: is_valid => amr_eb_tagging_criteria_is_valid
  end type amr_eb_tagging_criteria_2d

  type, public :: amr_eb_regrid_plan_2d
    logical :: active = .false.
    integer :: coarse_nx = 0
    integer :: coarse_ny = 0
    integer :: tagged_cell_count = 0
    integer :: tag_i_lower = 1
    integer :: tag_i_upper = 0
    integer :: tag_j_lower = 1
    integer :: tag_j_upper = 0
    integer :: coarse_i_lower = 1
    integer :: coarse_i_upper = 0
    integer :: coarse_j_lower = 1
    integer :: coarse_j_upper = 0
  contains
    procedure :: is_valid => amr_eb_regrid_plan_is_valid
  end type amr_eb_regrid_plan_2d

  public :: tag_reactive_eb_temperature_gradient_2d
  public :: build_amr_eb_regrid_plan_2d
  public :: plan_reactive_eb_temperature_regrid_2d
  public :: regrid_two_level_reactive_eb_patch_2d

contains

  pure logical function amr_eb_tagging_criteria_is_valid( &
      self, coarse_nx, coarse_ny) result(valid)
    class(amr_eb_tagging_criteria_2d), intent(in) :: self
    integer, intent(in), optional :: coarse_nx, coarse_ny

    valid = ieee_is_finite(self%relative_gradient_threshold) .and. &
      ieee_is_finite(self%absolute_gradient_threshold) .and. &
      ieee_is_finite(self%scale_floor) .and. &
      self%relative_gradient_threshold >= 0.0_dp .and. &
      self%absolute_gradient_threshold >= 0.0_dp .and. &
      self%scale_floor > 0.0_dp .and. self%buffer_cells >= 0 .and. &
      self%minimum_patch_cells_x >= 1 .and. &
      self%minimum_patch_cells_y >= 1
    if (.not. valid) return
    if (present(coarse_nx)) then
      valid = coarse_nx >= 4 .and. &
        self%minimum_patch_cells_x <= coarse_nx - 2
    end if
    if (.not. valid) return
    if (present(coarse_ny)) then
      valid = coarse_ny >= 4 .and. &
        self%minimum_patch_cells_y <= coarse_ny - 2
    end if
  end function amr_eb_tagging_criteria_is_valid

  pure logical function amr_eb_regrid_plan_is_valid(self) result(valid)
    class(amr_eb_regrid_plan_2d), intent(in) :: self

    valid = self%coarse_nx >= 4 .and. self%coarse_ny >= 4 .and. &
      self%tagged_cell_count >= 0
    if (.not. valid) return
    if (.not. self%active) then
      valid = self%tagged_cell_count == 0 .and. &
        self%tag_i_lower > self%tag_i_upper .and. &
        self%tag_j_lower > self%tag_j_upper .and. &
        self%coarse_i_lower > self%coarse_i_upper .and. &
        self%coarse_j_lower > self%coarse_j_upper
      return
    end if
    valid = self%tagged_cell_count >= 1 .and. &
      self%tag_i_lower >= 2 .and. &
      self%tag_i_upper <= self%coarse_nx - 1 .and. &
      self%tag_j_lower >= 2 .and. &
      self%tag_j_upper <= self%coarse_ny - 1 .and. &
      self%tag_i_lower <= self%tag_i_upper .and. &
      self%tag_j_lower <= self%tag_j_upper .and. &
      self%coarse_i_lower >= 2 .and. &
      self%coarse_i_upper <= self%coarse_nx - 1 .and. &
      self%coarse_j_lower >= 2 .and. &
      self%coarse_j_upper <= self%coarse_ny - 1 .and. &
      self%coarse_i_lower <= self%tag_i_lower .and. &
      self%coarse_i_upper >= self%tag_i_upper .and. &
      self%coarse_j_lower <= self%tag_j_lower .and. &
      self%coarse_j_upper >= self%tag_j_upper
  end function amr_eb_regrid_plan_is_valid

  subroutine tag_reactive_eb_temperature_gradient_2d( &
      temperature, geometry, criteria, tags, ok)
    real(dp), intent(in) :: temperature(:, :)
    type(eb_geometry_2d), intent(in) :: geometry
    type(amr_eb_tagging_criteria_2d), intent(in) :: criteria
    logical, intent(out) :: tags(:, :)
    logical, intent(out) :: ok

    real(dp) :: jump, local_scale
    integer :: i, j

    tags = .false.
    ok = geometry%is_valid() .and. &
      all(shape(temperature) == [geometry%nx, geometry%ny]) .and. &
      all(shape(tags) == [geometry%nx, geometry%ny]) .and. &
      criteria%is_valid(geometry%nx, geometry%ny) .and. &
      all(ieee_is_finite(temperature)) .and. all(temperature > 0.0_dp)
    if (.not. ok) return

    do j = 2, geometry%ny - 1
      do i = 2, geometry%nx - 1
        if (geometry%cell_type(i, j) == eb_covered_cell) cycle
        jump = 0.0_dp
        local_scale = max(criteria%scale_floor, abs(temperature(i, j)))
        call accumulate_neighbor(i - 1, j)
        call accumulate_neighbor(i + 1, j)
        call accumulate_neighbor(i, j - 1)
        call accumulate_neighbor(i, j + 1)
        tags(i, j) = jump > criteria%absolute_gradient_threshold .and. &
          jump / local_scale >= criteria%relative_gradient_threshold
      end do
    end do
    ok = .true.

  contains

    subroutine accumulate_neighbor(neighbor_i, neighbor_j)
      integer, intent(in) :: neighbor_i, neighbor_j

      if (geometry%cell_type(neighbor_i, neighbor_j) == &
          eb_covered_cell) return
      jump = max(jump, abs(temperature(neighbor_i, neighbor_j) - &
        temperature(i, j)))
      local_scale = max(local_scale, &
        abs(temperature(neighbor_i, neighbor_j)))
    end subroutine accumulate_neighbor

  end subroutine tag_reactive_eb_temperature_gradient_2d

  pure subroutine build_amr_eb_regrid_plan_2d( &
      tags, criteria, plan, ok)
    logical, intent(in) :: tags(:, :)
    type(amr_eb_tagging_criteria_2d), intent(in) :: criteria
    type(amr_eb_regrid_plan_2d), intent(out) :: plan
    logical, intent(out) :: ok

    integer :: i, j

    plan = amr_eb_regrid_plan_2d()
    plan%coarse_nx = size(tags, 1)
    plan%coarse_ny = size(tags, 2)
    ok = criteria%is_valid(plan%coarse_nx, plan%coarse_ny) .and. &
      .not. any(tags(1, :)) .and. .not. any(tags(plan%coarse_nx, :)) .and. &
      .not. any(tags(:, 1)) .and. .not. any(tags(:, plan%coarse_ny))
    if (.not. ok) return
    plan%tagged_cell_count = count(tags)
    if (plan%tagged_cell_count == 0) then
      ok = plan%is_valid()
      return
    end if

    plan%active = .true.
    plan%tag_i_lower = plan%coarse_nx
    plan%tag_i_upper = 1
    plan%tag_j_lower = plan%coarse_ny
    plan%tag_j_upper = 1
    do j = 2, plan%coarse_ny - 1
      do i = 2, plan%coarse_nx - 1
        if (.not. tags(i, j)) cycle
        plan%tag_i_lower = min(plan%tag_i_lower, i)
        plan%tag_i_upper = max(plan%tag_i_upper, i)
        plan%tag_j_lower = min(plan%tag_j_lower, j)
        plan%tag_j_upper = max(plan%tag_j_upper, j)
      end do
    end do
    plan%coarse_i_lower = max(2, &
      plan%tag_i_lower - criteria%buffer_cells)
    plan%coarse_i_upper = min(plan%coarse_nx - 1, &
      plan%tag_i_upper + criteria%buffer_cells)
    plan%coarse_j_lower = max(2, &
      plan%tag_j_lower - criteria%buffer_cells)
    plan%coarse_j_upper = min(plan%coarse_ny - 1, &
      plan%tag_j_upper + criteria%buffer_cells)
    call grow_interval( &
      plan%coarse_i_lower, plan%coarse_i_upper, plan%coarse_nx, &
      criteria%minimum_patch_cells_x)
    call grow_interval( &
      plan%coarse_j_lower, plan%coarse_j_upper, plan%coarse_ny, &
      criteria%minimum_patch_cells_y)
    ok = plan%is_valid() .and. &
      plan%coarse_i_upper - plan%coarse_i_lower + 1 >= &
        criteria%minimum_patch_cells_x .and. &
      plan%coarse_j_upper - plan%coarse_j_lower + 1 >= &
        criteria%minimum_patch_cells_y

  contains

    pure subroutine grow_interval(lower, upper, cell_count, minimum_cells)
      integer, intent(inout) :: lower, upper
      integer, intent(in) :: cell_count, minimum_cells

      do while (upper - lower + 1 < minimum_cells)
        if (lower > 2) lower = lower - 1
        if (upper - lower + 1 >= minimum_cells) exit
        if (upper < cell_count - 1) upper = upper + 1
      end do
    end subroutine grow_interval

  end subroutine build_amr_eb_regrid_plan_2d

  subroutine plan_reactive_eb_temperature_regrid_2d( &
      temperature, geometry, criteria, tags, plan, ok)
    real(dp), intent(in) :: temperature(:, :)
    type(eb_geometry_2d), intent(in) :: geometry
    type(amr_eb_tagging_criteria_2d), intent(in) :: criteria
    logical, intent(out) :: tags(:, :)
    type(amr_eb_regrid_plan_2d), intent(out) :: plan
    logical, intent(out) :: ok

    call tag_reactive_eb_temperature_gradient_2d( &
      temperature, geometry, criteria, tags, ok)
    if (.not. ok) return
    call build_amr_eb_regrid_plan_2d(tags, criteria, plan, ok)
  end subroutine plan_reactive_eb_temperature_regrid_2d

  subroutine regrid_two_level_reactive_eb_patch_2d( &
      species, coarse_state, coarse_temperature, coarse_geometry, &
      old_fine_state, old_fine_temperature, old_fine_geometry, old_patch, &
      new_fine_geometry, new_patch, new_coarse_state, &
      new_coarse_temperature, new_fine_state, new_fine_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: coarse_state(:, :, :), coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    real(dp), intent(in) :: old_fine_state(:, :, :)
    real(dp), intent(in) :: old_fine_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: old_fine_geometry
    type(amr_eb_patch_2d), intent(in) :: old_patch
    type(eb_geometry_2d), intent(in) :: new_fine_geometry
    type(amr_eb_patch_2d), intent(in) :: new_patch
    real(dp), intent(out) :: new_coarse_state(:, :, :)
    real(dp), intent(out) :: new_coarse_temperature(:, :)
    real(dp), intent(out) :: new_fine_state(:, :, :)
    real(dp), intent(out) :: new_fine_temperature(:, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate_coarse(:, :, :)
    real(dp), allocatable :: candidate_coarse_temperature(:, :)
    real(dp), allocatable :: candidate_fine(:, :, :)
    real(dp), allocatable :: candidate_fine_temperature(:, :)
    real(dp), allocatable :: primitive(:)
    real(dp) :: recovered_temperature, sound_speed
    logical :: local_ok
    integer :: global_i, global_j, i, j, old_i, old_j, nvar, ratio

    new_coarse_state = coarse_state
    new_coarse_temperature = coarse_temperature
    new_fine_state = 0.0_dp
    new_fine_temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar < 1 .or. old_patch%refinement_ratio /= &
        new_patch%refinement_ratio .or. &
        .not. old_patch%is_valid(coarse_geometry, old_fine_geometry) .or. &
        .not. new_patch%is_valid(coarse_geometry, new_fine_geometry) .or. &
        any(shape(new_coarse_state) /= shape(coarse_state)) .or. &
        any(shape(new_coarse_temperature) /= shape(coarse_temperature)) .or. &
        any(shape(new_fine_state) /= &
          [nvar, new_fine_geometry%nx, new_fine_geometry%ny]) .or. &
        any(shape(new_fine_temperature) /= &
          [new_fine_geometry%nx, new_fine_geometry%ny])) return

    allocate(candidate_coarse, mold=coarse_state)
    allocate(candidate_coarse_temperature, mold=coarse_temperature)
    call average_down_reactive_eb_state_patch_2d( &
      species, coarse_state, coarse_temperature, coarse_geometry, &
      old_fine_state, old_fine_geometry, old_patch, candidate_coarse, &
      candidate_coarse_temperature, local_ok)
    if (.not. local_ok) return
    allocate(candidate_fine, mold=new_fine_state)
    allocate(candidate_fine_temperature, mold=new_fine_temperature)
    call prolong_reactive_eb_patch_pcm_2d( &
      species, candidate_coarse, candidate_coarse_temperature, &
      coarse_geometry, new_fine_geometry, new_patch, candidate_fine, &
      candidate_fine_temperature, local_ok)
    if (.not. local_ok) return

    ratio = new_patch%refinement_ratio
    do j = 1, new_fine_geometry%ny
      global_j = (new_patch%coarse_j_lower - 1) * ratio + j
      old_j = global_j - (old_patch%coarse_j_lower - 1) * ratio
      if (old_j < 1 .or. old_j > old_fine_geometry%ny) cycle
      do i = 1, new_fine_geometry%nx
        global_i = (new_patch%coarse_i_lower - 1) * ratio + i
        old_i = global_i - (old_patch%coarse_i_lower - 1) * ratio
        if (old_i < 1 .or. old_i > old_fine_geometry%nx) cycle
        candidate_fine(:, i, j) = old_fine_state(:, old_i, old_j)
        candidate_fine_temperature(i, j) = &
          old_fine_temperature(old_i, old_j)
      end do
    end do
    if (any(.not. ieee_is_finite(candidate_coarse)) .or. &
        any(.not. ieee_is_finite(candidate_coarse_temperature)) .or. &
        any(.not. ieee_is_finite(candidate_fine)) .or. &
        any(.not. ieee_is_finite(candidate_fine_temperature))) return
    allocate(primitive(reactive_nprim(size(species))))
    do j = 1, new_fine_geometry%ny
      do i = 1, new_fine_geometry%nx
        if (new_fine_geometry%cell_type(i, j) == eb_covered_cell) cycle
        if (candidate_fine_temperature(i, j) <= 0.0_dp) return
        call reactive_conserved_to_primitive( &
          species, candidate_fine(:, i, j), &
          candidate_fine_temperature(i, j), primitive, &
          recovered_temperature, sound_speed, local_ok)
        if (.not. local_ok) return
        candidate_fine_temperature(i, j) = recovered_temperature
      end do
    end do

    new_coarse_state = candidate_coarse
    new_coarse_temperature = candidate_coarse_temperature
    new_fine_state = candidate_fine
    new_fine_temperature = candidate_fine_temperature
    ok = .true.
  end subroutine regrid_two_level_reactive_eb_patch_2d

end module amr_eb_regrid_2d_mod
