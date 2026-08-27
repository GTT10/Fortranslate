module amr_eb_patch_tree_2d_mod
  use precision_mod, only: dp
  use eb_geometry_2d_mod, only: eb_geometry_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, build_amr_eb_patch_2d
  implicit none
  private

  integer, parameter :: patch_tree_separation_cells = 2

  type, public :: amr_eb_patch_tree_child_plan_2d
    integer :: parent_patch = 0
    integer :: coarse_i_lower = 1
    integer :: coarse_i_upper = 0
    integer :: coarse_j_lower = 1
    integer :: coarse_j_upper = 0
    type(eb_geometry_2d) :: geometry
  end type amr_eb_patch_tree_child_plan_2d

  type, public :: amr_eb_patch_tree_level_plan_2d
    integer :: refinement_ratio = 0
    type(amr_eb_patch_tree_child_plan_2d), allocatable :: children(:)
  contains
    procedure :: patch_count => amr_eb_patch_tree_plan_patch_count
  end type amr_eb_patch_tree_level_plan_2d

  type, public :: amr_eb_patch_tree_node_2d
    integer :: parent_patch = 0
    type(eb_geometry_2d) :: geometry
    type(amr_eb_patch_2d) :: patch
  end type amr_eb_patch_tree_node_2d

  type, public :: amr_eb_patch_tree_relation_2d
    integer :: level = 0
    integer :: refinement_ratio = 0
    type(amr_eb_patch_tree_node_2d), allocatable :: children(:)
    integer, allocatable :: child_offsets(:)
  contains
    procedure :: parent_patch_count => amr_eb_patch_tree_parent_patch_count
    procedure :: child_patch_count => amr_eb_patch_tree_child_patch_count
    procedure :: child_index => amr_eb_patch_tree_child_index
    procedure :: is_valid => amr_eb_patch_tree_relation_is_valid
  end type amr_eb_patch_tree_relation_2d

  type, public :: amr_eb_patch_tree_topology_2d
    type(eb_geometry_2d) :: root_geometry
    type(amr_eb_patch_tree_relation_2d), allocatable :: relations(:)
  contains
    procedure :: level_count => amr_eb_patch_tree_level_count
    procedure :: level_patch_count => amr_eb_patch_tree_level_patch_count
    procedure :: is_valid => amr_eb_patch_tree_topology_is_valid
  end type amr_eb_patch_tree_topology_2d

  public :: initialize_amr_eb_patch_tree_topology_2d
  public :: rebuild_amr_eb_patch_tree_topology_2d
  public :: patch_tree_topologies_match_2d
  public :: eb_geometries_match_2d

contains

  pure integer function amr_eb_patch_tree_plan_patch_count(self) &
      result(count)
    class(amr_eb_patch_tree_level_plan_2d), intent(in) :: self

    count = 0
    if (allocated(self%children)) count = size(self%children)
  end function amr_eb_patch_tree_plan_patch_count

  pure integer function amr_eb_patch_tree_parent_patch_count(self) &
      result(count)
    class(amr_eb_patch_tree_relation_2d), intent(in) :: self

    count = 0
    if (.not. allocated(self%child_offsets)) return
    if (size(self%child_offsets) < 1) return
    count = size(self%child_offsets) - 1
  end function amr_eb_patch_tree_parent_patch_count

  pure integer function amr_eb_patch_tree_child_patch_count(self) &
      result(count)
    class(amr_eb_patch_tree_relation_2d), intent(in) :: self

    count = 0
    if (allocated(self%children)) count = size(self%children)
  end function amr_eb_patch_tree_child_patch_count

  pure integer function amr_eb_patch_tree_child_index( &
      self, parent_patch, local_child) result(index)
    class(amr_eb_patch_tree_relation_2d), intent(in) :: self
    integer, intent(in) :: parent_patch, local_child

    integer :: first_child, last_child

    index = 0
    if (.not. allocated(self%child_offsets)) return
    if (parent_patch < 1 .or. &
        parent_patch > self%parent_patch_count()) return
    first_child = self%child_offsets(parent_patch) + 1
    last_child = self%child_offsets(parent_patch + 1)
    if (local_child < 1 .or. &
        local_child > last_child - first_child + 1) return
    index = first_child + local_child - 1
  end function amr_eb_patch_tree_child_index

  pure logical function amr_eb_patch_tree_relation_is_valid( &
      self, parent_geometries) result(valid)
    class(amr_eb_patch_tree_relation_2d), intent(in) :: self
    type(eb_geometry_2d), intent(in) :: parent_geometries(:)

    integer :: child, first_child, local_child, parent, second_child

    valid = self%level >= 1 .and. self%refinement_ratio >= 2 .and. &
      size(parent_geometries) >= 1 .and. allocated(self%children) .and. &
      allocated(self%child_offsets)
    if (.not. valid) return
    valid = size(self%children) >= 1 .and. &
      size(self%child_offsets) == size(parent_geometries) + 1 .and. &
      self%child_offsets(1) == 0 .and. &
      self%child_offsets(size(self%child_offsets)) == size(self%children)
    if (.not. valid) return

    do parent = 1, size(parent_geometries)
      valid = parent_geometries(parent)%is_valid() .and. &
        self%child_offsets(parent + 1) >= self%child_offsets(parent)
      if (.not. valid) return
      first_child = self%child_offsets(parent) + 1
      do child = first_child, self%child_offsets(parent + 1)
        local_child = child - first_child + 1
        valid = self%children(child)%parent_patch == parent .and. &
          self%child_index(parent, local_child) == child .and. &
          self%children(child)%patch%refinement_ratio == &
            self%refinement_ratio .and. &
          self%children(child)%patch%is_valid( &
            parent_geometries(parent), self%children(child)%geometry)
        if (.not. valid) return
        do second_child = child + 1, self%child_offsets(parent + 1)
          valid = patch_tree_siblings_are_separated_2d( &
            self%children(child)%patch, &
            self%children(second_child)%patch)
          if (.not. valid) return
        end do
      end do
    end do
  end function amr_eb_patch_tree_relation_is_valid

  pure integer function amr_eb_patch_tree_level_count(self) result(count)
    class(amr_eb_patch_tree_topology_2d), intent(in) :: self

    count = 0
    if (allocated(self%relations)) count = size(self%relations) + 1
  end function amr_eb_patch_tree_level_count

  pure integer function amr_eb_patch_tree_level_patch_count( &
      self, level) result(count)
    class(amr_eb_patch_tree_topology_2d), intent(in) :: self
    integer, intent(in) :: level

    count = 0
    if (level < 0 .or. level >= self%level_count()) return
    if (level == 0) then
      count = 1
    else
      count = self%relations(level)%child_patch_count()
    end if
  end function amr_eb_patch_tree_level_patch_count

  logical function amr_eb_patch_tree_topology_is_valid(self) result(valid)
    class(amr_eb_patch_tree_topology_2d), intent(in) :: self

    type(eb_geometry_2d), allocatable :: parent_geometries(:)
    integer :: child, relation

    valid = self%root_geometry%is_valid() .and. allocated(self%relations)
    if (.not. valid) return
    allocate(parent_geometries(1))
    parent_geometries(1) = self%root_geometry
    do relation = 1, size(self%relations)
      valid = self%relations(relation)%level == relation .and. &
        self%relations(relation)%is_valid(parent_geometries)
      if (.not. valid) return
      deallocate(parent_geometries)
      allocate(parent_geometries( &
        self%relations(relation)%child_patch_count()))
      do child = 1, size(parent_geometries)
        parent_geometries(child) = &
          self%relations(relation)%children(child)%geometry
      end do
    end do
  end function amr_eb_patch_tree_topology_is_valid

  subroutine initialize_amr_eb_patch_tree_topology_2d( &
      root_geometry, plans, topology, ok)
    type(eb_geometry_2d), intent(in) :: root_geometry
    type(amr_eb_patch_tree_level_plan_2d), intent(in) :: plans(:)
    type(amr_eb_patch_tree_topology_2d), intent(out) :: topology
    logical, intent(out) :: ok

    type(amr_eb_patch_tree_topology_2d) :: candidate
    type(eb_geometry_2d), allocatable :: parent_geometries(:)
    logical :: local_ok
    integer :: child, parent, relation

    topology = amr_eb_patch_tree_topology_2d()
    ok = root_geometry%is_valid()
    if (.not. ok) return
    candidate%root_geometry = root_geometry
    allocate(candidate%relations(size(plans)))
    if (size(plans) == 0) then
      topology = candidate
      ok = topology%is_valid()
      return
    end if

    allocate(parent_geometries(1))
    parent_geometries(1) = root_geometry
    do relation = 1, size(plans)
      local_ok = plans(relation)%refinement_ratio >= 2 .and. &
        allocated(plans(relation)%children)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      local_ok = size(plans(relation)%children) >= 1
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      candidate%relations(relation)%level = relation
      candidate%relations(relation)%refinement_ratio = &
        plans(relation)%refinement_ratio
      allocate(candidate%relations(relation)%children( &
        plans(relation)%patch_count()))
      allocate(candidate%relations(relation)%child_offsets( &
        size(parent_geometries) + 1))
      candidate%relations(relation)%child_offsets(1) = 0
      child = 1
      do parent = 1, size(parent_geometries)
        do while (child <= plans(relation)%patch_count())
          if (plans(relation)%children(child)%parent_patch /= parent) exit
          candidate%relations(relation)%children(child)%parent_patch = parent
          candidate%relations(relation)%children(child)%geometry = &
            plans(relation)%children(child)%geometry
          call build_amr_eb_patch_2d( &
            parent_geometries(parent), &
            plans(relation)%children(child)%geometry, &
            plans(relation)%children(child)%coarse_i_lower, &
            plans(relation)%children(child)%coarse_i_upper, &
            plans(relation)%children(child)%coarse_j_lower, &
            plans(relation)%children(child)%coarse_j_upper, &
            plans(relation)%refinement_ratio, &
            candidate%relations(relation)%children(child)%patch, local_ok)
          if (.not. local_ok) then
            ok = .false.
            return
          end if
          child = child + 1
        end do
        candidate%relations(relation)%child_offsets(parent + 1) = child - 1
      end do
      if (child <= plans(relation)%patch_count()) then
        ok = .false.
        return
      end if
      local_ok = candidate%relations(relation)%is_valid(parent_geometries)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      deallocate(parent_geometries)
      allocate(parent_geometries( &
        candidate%relations(relation)%child_patch_count()))
      do child = 1, size(parent_geometries)
        parent_geometries(child) = &
          candidate%relations(relation)%children(child)%geometry
      end do
    end do

    ok = candidate%is_valid()
    if (.not. ok) return
    topology = candidate
  end subroutine initialize_amr_eb_patch_tree_topology_2d

  subroutine rebuild_amr_eb_patch_tree_topology_2d( &
      topology, plans, ok, changed)
    type(amr_eb_patch_tree_topology_2d), intent(inout) :: topology
    type(amr_eb_patch_tree_level_plan_2d), intent(in) :: plans(:)
    logical, intent(out) :: ok, changed

    type(amr_eb_patch_tree_topology_2d) :: candidate

    ok = .false.
    changed = .false.
    if (.not. topology%is_valid()) return
    call initialize_amr_eb_patch_tree_topology_2d( &
      topology%root_geometry, plans, candidate, ok)
    if (.not. ok) return
    changed = .not. patch_tree_topologies_match_2d(topology, candidate)
    if (changed) topology = candidate
  end subroutine rebuild_amr_eb_patch_tree_topology_2d

  pure logical function patch_tree_siblings_are_separated_2d( &
      first, second) result(separated)
    type(amr_eb_patch_2d), intent(in) :: first, second

    separated = &
      first%coarse_i_upper + patch_tree_separation_cells < &
        second%coarse_i_lower .or. &
      second%coarse_i_upper + patch_tree_separation_cells < &
        first%coarse_i_lower .or. &
      first%coarse_j_upper + patch_tree_separation_cells < &
        second%coarse_j_lower .or. &
      second%coarse_j_upper + patch_tree_separation_cells < &
        first%coarse_j_lower
  end function patch_tree_siblings_are_separated_2d

  logical function patch_tree_topologies_match_2d(first, second) &
      result(matches)
    type(amr_eb_patch_tree_topology_2d), intent(in) :: first, second

    integer :: child, relation

    matches = first%is_valid() .and. second%is_valid()
    if (.not. matches) return
    matches = eb_geometries_match_2d( &
      first%root_geometry, second%root_geometry) .and. &
      size(first%relations) == size(second%relations)
    if (.not. matches) return
    do relation = 1, size(first%relations)
      matches = first%relations(relation)%level == &
          second%relations(relation)%level .and. &
        first%relations(relation)%refinement_ratio == &
          second%relations(relation)%refinement_ratio .and. &
        size(first%relations(relation)%children) == &
          size(second%relations(relation)%children) .and. &
        all(first%relations(relation)%child_offsets == &
          second%relations(relation)%child_offsets)
      if (.not. matches) return
      do child = 1, first%relations(relation)%child_patch_count()
        matches = first%relations(relation)%children(child)%parent_patch == &
            second%relations(relation)%children(child)%parent_patch .and. &
          patches_match_2d( &
            first%relations(relation)%children(child)%patch, &
            second%relations(relation)%children(child)%patch) .and. &
          eb_geometries_match_2d( &
            first%relations(relation)%children(child)%geometry, &
            second%relations(relation)%children(child)%geometry)
        if (.not. matches) return
      end do
    end do
  end function patch_tree_topologies_match_2d

  pure logical function patches_match_2d(first, second) result(matches)
    type(amr_eb_patch_2d), intent(in) :: first, second

    matches = all([ &
      first%coarse_i_lower, first%coarse_i_upper, &
      first%coarse_j_lower, first%coarse_j_upper, &
      first%refinement_ratio] == [ &
      second%coarse_i_lower, second%coarse_i_upper, &
      second%coarse_j_lower, second%coarse_j_upper, &
      second%refinement_ratio])
  end function patches_match_2d

  pure logical function eb_geometries_match_2d(first, second) result(matches)
    type(eb_geometry_2d), intent(in) :: first, second

    real(dp), parameter :: tolerance = 5.0e3_dp * epsilon(1.0_dp)

    matches = first%is_valid() .and. second%is_valid()
    if (.not. matches) return
    matches = first%nx == second%nx .and. first%ny == second%ny
    if (.not. matches) return
    matches = all(abs([ &
      first%x_lower, first%x_upper, first%y_lower, first%y_upper, &
      first%dx, first%dy] - [ &
      second%x_lower, second%x_upper, second%y_lower, second%y_upper, &
      second%dx, second%dy]) <= tolerance) .and. &
      all(first%cell_type == second%cell_type) .and. &
      all(abs(first%volume_fraction - second%volume_fraction) <= tolerance) &
      .and. all(abs(first%cell_centroid_x - &
        second%cell_centroid_x) <= tolerance) .and. &
      all(abs(first%cell_centroid_y - &
        second%cell_centroid_y) <= tolerance) .and. &
      all(abs(first%x_face_fraction - &
        second%x_face_fraction) <= tolerance) .and. &
      all(abs(first%y_face_fraction - &
        second%y_face_fraction) <= tolerance) .and. &
      all(abs(first%x_face_centroid_y - &
        second%x_face_centroid_y) <= tolerance) .and. &
      all(abs(first%y_face_centroid_x - &
        second%y_face_centroid_x) <= tolerance) .and. &
      all(abs(first%boundary_length - &
        second%boundary_length) <= tolerance) .and. &
      all(abs(first%boundary_centroid_x - &
        second%boundary_centroid_x) <= tolerance) .and. &
      all(abs(first%boundary_centroid_y - &
        second%boundary_centroid_y) <= tolerance) .and. &
      all(abs(first%boundary_normal_x - &
        second%boundary_normal_x) <= tolerance) .and. &
      all(abs(first%boundary_normal_y - &
        second%boundary_normal_y) <= tolerance) .and. &
      all(abs(first%boundary_normal_integral_x - &
        second%boundary_normal_integral_x) <= tolerance) .and. &
      all(abs(first%boundary_normal_integral_y - &
        second%boundary_normal_integral_y) <= tolerance)
  end function eb_geometries_match_2d

end module amr_eb_patch_tree_2d_mod
