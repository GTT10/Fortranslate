module amr_eb_patch_tree_reactive_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_covered_cell
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, average_down_reactive_eb_state_patch_2d
  use amr_eb_reactive_2d_mod, only: prolong_reactive_eb_patch_pcm_2d
  use amr_eb_patch_tree_2d_mod, only: &
    amr_eb_patch_tree_level_plan_2d, amr_eb_patch_tree_topology_2d, &
    rebuild_amr_eb_patch_tree_topology_2d
  implicit none
  private

  real(dp), parameter :: geometry_tolerance = &
    5.0e3_dp * epsilon(1.0_dp)
  real(dp), parameter :: conservation_tolerance = &
    5.0e4_dp * epsilon(1.0_dp)

  type, public :: reactive_amr_eb_patch_tree_node_2d
    real(dp), allocatable :: state(:, :, :)
    real(dp), allocatable :: temperature(:, :)
  end type reactive_amr_eb_patch_tree_node_2d

  type, public :: reactive_amr_eb_patch_tree_level_2d
    type(reactive_amr_eb_patch_tree_node_2d), allocatable :: patches(:)
  contains
    procedure :: patch_count => reactive_amr_eb_patch_tree_level_patch_count
  end type reactive_amr_eb_patch_tree_level_2d

  type, public :: reactive_amr_eb_patch_tree_2d
    integer :: nvar = 0
    type(amr_eb_patch_tree_topology_2d) :: topology
    type(reactive_amr_eb_patch_tree_level_2d), allocatable :: levels(:)
  contains
    procedure :: level_count => reactive_amr_eb_patch_tree_level_count
    procedure :: level_patch_count => &
      reactive_amr_eb_patch_tree_level_patch_count_at
    procedure :: is_valid => reactive_amr_eb_patch_tree_is_valid
  end type reactive_amr_eb_patch_tree_2d

  public :: initialize_reactive_amr_eb_patch_tree_2d
  public :: synchronize_reactive_amr_eb_patch_tree_2d
  public :: rebuild_reactive_amr_eb_patch_tree_2d
  public :: composite_integral_reactive_amr_eb_patch_tree_2d

contains

  pure integer function reactive_amr_eb_patch_tree_level_patch_count(self) &
      result(count)
    class(reactive_amr_eb_patch_tree_level_2d), intent(in) :: self

    count = 0
    if (allocated(self%patches)) count = size(self%patches)
  end function reactive_amr_eb_patch_tree_level_patch_count

  pure integer function reactive_amr_eb_patch_tree_level_count(self) &
      result(count)
    class(reactive_amr_eb_patch_tree_2d), intent(in) :: self

    count = 0
    if (allocated(self%levels)) count = size(self%levels)
  end function reactive_amr_eb_patch_tree_level_count

  pure integer function reactive_amr_eb_patch_tree_level_patch_count_at( &
      self, level) result(count)
    class(reactive_amr_eb_patch_tree_2d), intent(in) :: self
    integer, intent(in) :: level

    count = 0
    if (.not. allocated(self%levels)) return
    if (level < 0 .or. level >= size(self%levels)) return
    count = self%levels(level + 1)%patch_count()
  end function reactive_amr_eb_patch_tree_level_patch_count_at

  logical function reactive_amr_eb_patch_tree_is_valid(self) result(valid)
    class(reactive_amr_eb_patch_tree_2d), intent(in) :: self

    type(eb_geometry_2d) :: geometry
    integer :: level, patch

    valid = self%nvar >= 1 .and. self%topology%is_valid() .and. &
      allocated(self%levels)
    if (.not. valid) return
    valid = size(self%levels) == self%topology%level_count()
    if (.not. valid) return

    do level = 1, size(self%levels)
      valid = allocated(self%levels(level)%patches) .and. &
        self%levels(level)%patch_count() == &
          self%topology%level_patch_count(level - 1)
      if (.not. valid) return
      do patch = 1, self%levels(level)%patch_count()
        call patch_geometry_at(self%topology, level, patch, geometry, valid)
        if (.not. valid) return
        valid = allocated(self%levels(level)%patches(patch)%state) .and. &
          allocated(self%levels(level)%patches(patch)%temperature)
        if (.not. valid) return
        valid = all(shape(self%levels(level)%patches(patch)%state) == &
            [self%nvar, geometry%nx, geometry%ny]) .and. &
          all(shape(self%levels(level)%patches(patch)%temperature) == &
            [geometry%nx, geometry%ny])
        if (.not. valid) return
        valid = all(ieee_is_finite( &
            self%levels(level)%patches(patch)%state)) .and. &
          all(ieee_is_finite( &
            self%levels(level)%patches(patch)%temperature)) .and. &
          all(self%levels(level)%patches(patch)%temperature > 0.0_dp)
        if (.not. valid) return
      end do
    end do
  end function reactive_amr_eb_patch_tree_is_valid

  subroutine initialize_reactive_amr_eb_patch_tree_2d( &
      species, root_state, root_temperature, topology, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: root_state(:, :, :), root_temperature(:, :)
    type(amr_eb_patch_tree_topology_2d), intent(in) :: topology
    type(reactive_amr_eb_patch_tree_2d), intent(out) :: solution
    logical, intent(out) :: ok

    type(reactive_amr_eb_patch_tree_2d) :: candidate
    type(eb_geometry_2d) :: parent_geometry
    logical :: local_ok
    integer :: child, nvar, parent, relation

    solution = reactive_amr_eb_patch_tree_2d()
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar < 1 .or. .not. topology%is_valid()) return
    if (any(shape(root_state) /= &
          [nvar, topology%root_geometry%nx, topology%root_geometry%ny]) &
        .or. any(shape(root_temperature) /= &
          [topology%root_geometry%nx, topology%root_geometry%ny]) .or. &
        any(.not. ieee_is_finite(root_state)) .or. &
        any(.not. ieee_is_finite(root_temperature)) .or. &
        any(root_temperature <= 0.0_dp)) return

    candidate%nvar = nvar
    candidate%topology = topology
    allocate(candidate%levels(topology%level_count()))
    allocate(candidate%levels(1)%patches(1))
    allocate(candidate%levels(1)%patches(1)%state, source=root_state)
    allocate(candidate%levels(1)%patches(1)%temperature, &
      source=root_temperature)
    call recover_patch_temperature( &
      species, candidate%levels(1)%patches(1), topology%root_geometry, &
      local_ok)
    if (.not. local_ok) return

    do relation = 1, size(topology%relations)
      allocate(candidate%levels(relation + 1)%patches( &
        topology%relations(relation)%child_patch_count()))
      do child = 1, topology%relations(relation)%child_patch_count()
        parent = topology%relations(relation)%children(child)%parent_patch
        call patch_geometry_at( &
          topology, relation, parent, parent_geometry, local_ok)
        if (.not. local_ok) return
        allocate(candidate%levels(relation + 1)%patches(child)%state( &
          nvar, topology%relations(relation)%children(child)%geometry%nx, &
          topology%relations(relation)%children(child)%geometry%ny))
        allocate(candidate%levels(relation + 1)%patches(child)%temperature( &
          topology%relations(relation)%children(child)%geometry%nx, &
          topology%relations(relation)%children(child)%geometry%ny))
        call prolong_reactive_eb_patch_pcm_2d( &
          species, candidate%levels(relation)%patches(parent)%state, &
          candidate%levels(relation)%patches(parent)%temperature, &
          parent_geometry, &
          topology%relations(relation)%children(child)%geometry, &
          topology%relations(relation)%children(child)%patch, &
          candidate%levels(relation + 1)%patches(child)%state, &
          candidate%levels(relation + 1)%patches(child)%temperature, &
          local_ok)
        if (.not. local_ok) return
      end do
    end do

    if (.not. candidate%is_valid()) return
    solution = candidate
    ok = .true.
  end subroutine initialize_reactive_amr_eb_patch_tree_2d

  subroutine synchronize_reactive_amr_eb_patch_tree_2d( &
      species, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_amr_eb_patch_tree_2d), intent(inout) :: solution
    logical, intent(out) :: ok

    type(reactive_amr_eb_patch_tree_2d) :: candidate

    ok = .false.
    if (.not. solution%is_valid()) return
    if (solution%nvar /= reactive_nvar(size(species))) return
    candidate = solution
    call synchronize_candidate(species, candidate, ok)
    if (.not. ok) return
    if (.not. candidate%is_valid()) then
      ok = .false.
      return
    end if
    solution = candidate
  end subroutine synchronize_reactive_amr_eb_patch_tree_2d

  subroutine rebuild_reactive_amr_eb_patch_tree_2d( &
      species, solution, plans, ok, changed, failure_context)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_amr_eb_patch_tree_2d), intent(inout) :: solution
    type(amr_eb_patch_tree_level_plan_2d), intent(in) :: plans(:)
    logical, intent(out) :: ok, changed
    character(len=*), intent(out), optional :: failure_context

    type(reactive_amr_eb_patch_tree_2d) :: collapsed, candidate
    type(amr_eb_patch_tree_topology_2d) :: new_topology
    type(eb_geometry_2d) :: parent_geometry
    real(dp), allocatable :: old_integral(:), new_integral(:)
    real(dp) :: integral_scale
    logical :: local_ok, topology_changed
    integer :: child, level, old_patch, parent, relation
    logical, allocatable :: copied(:, :)

    ok = .false.
    changed = .false.
    if (present(failure_context)) failure_context = "source validation"
    if (.not. solution%is_valid()) return
    if (solution%nvar /= reactive_nvar(size(species))) return

    if (present(failure_context)) failure_context = "topology rebuild"
    new_topology = solution%topology
    call rebuild_amr_eb_patch_tree_topology_2d( &
      new_topology, plans, local_ok, topology_changed)
    if (.not. local_ok) return
    if (.not. topology_changed) then
      ok = .true.
      if (present(failure_context)) failure_context = "none"
      return
    end if

    if (present(failure_context)) failure_context = "source integral"
    allocate(old_integral(solution%nvar), new_integral(solution%nvar))
    call composite_integral_reactive_amr_eb_patch_tree_2d( &
      solution, old_integral, local_ok)
    if (.not. local_ok) return
    if (present(failure_context)) failure_context = "source synchronization"
    collapsed = solution
    call synchronize_candidate(species, collapsed, local_ok)
    if (.not. local_ok) return
    if (present(failure_context)) failure_context = "candidate initialization"
    call initialize_reactive_amr_eb_patch_tree_2d( &
      species, collapsed%levels(1)%patches(1)%state, &
      collapsed%levels(1)%patches(1)%temperature, new_topology, &
      candidate, local_ok)
    if (.not. local_ok) return

    do level = 2, candidate%level_count()
      relation = level - 1
      do child = 1, candidate%levels(level)%patch_count()
        parent = candidate%topology%relations(relation)% &
          children(child)%parent_patch
        call patch_geometry_at( &
          candidate%topology, level - 1, parent, parent_geometry, local_ok)
        if (.not. local_ok) return
        if (present(failure_context)) &
          write(failure_context, '(a,i0,a,i0)') &
            "PCM initialization level ", level - 1, " patch ", child
        call prolong_reactive_eb_patch_pcm_2d( &
          species, candidate%levels(level - 1)%patches(parent)%state, &
          candidate%levels(level - 1)%patches(parent)%temperature, &
          parent_geometry, &
          candidate%topology%relations(relation)%children(child)%geometry, &
          candidate%topology%relations(relation)%children(child)%patch, &
          candidate%levels(level)%patches(child)%state, &
          candidate%levels(level)%patches(child)%temperature, local_ok)
        if (.not. local_ok) return
        if (present(failure_context)) &
          write(failure_context, '(a,i0,a,i0)') &
            "overlap retention level ", level - 1, " patch ", child
        allocate(copied( &
          size(candidate%levels(level)%patches(child)%temperature, 1), &
          size(candidate%levels(level)%patches(child)%temperature, 2)))
        copied = .false.
        if (level <= collapsed%level_count()) then
          do old_patch = 1, collapsed%levels(level)%patch_count()
            call retain_same_resolution_overlap( &
              candidate%levels(level)%patches(child), &
              candidate%topology%relations(relation)%children(child)%geometry, &
              collapsed%levels(level)%patches(old_patch), &
              collapsed%topology%relations(relation)%children(old_patch)% &
                geometry, copied, local_ok)
            if (.not. local_ok) return
          end do
        end if
        deallocate(copied)
        if (present(failure_context)) &
          write(failure_context, '(a,i0,a,i0)') &
            "temperature recovery level ", level - 1, " patch ", child
        call recover_patch_temperature( &
          species, candidate%levels(level)%patches(child), &
          candidate%topology%relations(relation)%children(child)%geometry, &
          local_ok)
        if (.not. local_ok) return
      end do
    end do

    if (present(failure_context)) failure_context = "candidate synchronization"
    call synchronize_candidate(species, candidate, local_ok)
    if (.not. local_ok) return
    if (.not. candidate%is_valid()) return
    if (present(failure_context)) failure_context = "candidate integral"
    call composite_integral_reactive_amr_eb_patch_tree_2d( &
      candidate, new_integral, local_ok)
    if (.not. local_ok) return
    integral_scale = max(1.0_dp, maxval(abs(old_integral)))
    if (present(failure_context)) failure_context = "conservation check"
    if (maxval(abs(new_integral - old_integral)) > &
        conservation_tolerance * integral_scale) return

    solution = candidate
    ok = .true.
    changed = .true.
    if (present(failure_context)) failure_context = "none"
  end subroutine rebuild_reactive_amr_eb_patch_tree_2d

  subroutine composite_integral_reactive_amr_eb_patch_tree_2d( &
      solution, integral, ok)
    type(reactive_amr_eb_patch_tree_2d), intent(in) :: solution
    real(dp), intent(out) :: integral(:)
    logical, intent(out) :: ok

    type(eb_geometry_2d) :: parent_geometry
    type(amr_eb_patch_2d) :: patch
    integer :: child, component, i, j, parent, relation

    integral = 0.0_dp
    ok = .false.
    if (.not. solution%is_valid()) return
    if (size(integral) /= solution%nvar) return

    do component = 1, solution%nvar
      integral(component) = sum( &
        solution%topology%root_geometry%volume_fraction * &
        solution%levels(1)%patches(1)%state(component, :, :)) * &
        solution%topology%root_geometry%dx * &
        solution%topology%root_geometry%dy
    end do

    do relation = 1, size(solution%topology%relations)
      do child = 1, &
          solution%topology%relations(relation)%child_patch_count()
        parent = solution%topology%relations(relation)% &
          children(child)%parent_patch
        call patch_geometry_at( &
          solution%topology, relation, parent, parent_geometry, ok)
        if (.not. ok) return
        patch = solution%topology%relations(relation)%children(child)%patch
        do component = 1, solution%nvar
          integral(component) = integral(component) + sum( &
            solution%topology%relations(relation)%children(child)% &
              geometry%volume_fraction * &
            solution%levels(relation + 1)%patches(child)% &
              state(component, :, :)) * &
            solution%topology%relations(relation)%children(child)% &
              geometry%dx * &
            solution%topology%relations(relation)%children(child)% &
              geometry%dy
          do j = patch%coarse_j_lower, patch%coarse_j_upper
            do i = patch%coarse_i_lower, patch%coarse_i_upper
              integral(component) = integral(component) - &
                parent_geometry%volume_fraction(i, j) * &
                solution%levels(relation)%patches(parent)% &
                  state(component, i, j) * &
                parent_geometry%dx * parent_geometry%dy
            end do
          end do
        end do
      end do
    end do
    ok = all(ieee_is_finite(integral))
  end subroutine composite_integral_reactive_amr_eb_patch_tree_2d

  subroutine synchronize_candidate(species, candidate, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_amr_eb_patch_tree_2d), intent(inout) :: candidate
    logical, intent(out) :: ok

    type(eb_geometry_2d) :: geometry, parent_geometry
    real(dp), allocatable :: state_work(:, :, :), temperature_work(:, :)
    logical :: local_ok
    integer :: child, level, parent, patch, relation

    ok = .false.
    if (.not. candidate%is_valid()) return
    do level = 1, candidate%level_count()
      do patch = 1, candidate%levels(level)%patch_count()
        call patch_geometry_at( &
          candidate%topology, level, patch, geometry, local_ok)
        if (.not. local_ok) return
        call recover_patch_temperature( &
          species, candidate%levels(level)%patches(patch), geometry, local_ok)
        if (.not. local_ok) return
      end do
    end do
    do relation = size(candidate%topology%relations), 1, -1
      do child = 1, &
          candidate%topology%relations(relation)%child_patch_count()
        parent = candidate%topology%relations(relation)% &
          children(child)%parent_patch
        call patch_geometry_at( &
          candidate%topology, relation, parent, parent_geometry, local_ok)
        if (.not. local_ok) return
        allocate(state_work, mold= &
          candidate%levels(relation)%patches(parent)%state)
        allocate(temperature_work, mold= &
          candidate%levels(relation)%patches(parent)%temperature)
        call average_down_reactive_eb_state_patch_2d( &
          species, candidate%levels(relation)%patches(parent)%state, &
          candidate%levels(relation)%patches(parent)%temperature, &
          parent_geometry, &
          candidate%levels(relation + 1)%patches(child)%state, &
          candidate%topology%relations(relation)%children(child)%geometry, &
          candidate%topology%relations(relation)%children(child)%patch, &
          state_work, temperature_work, local_ok)
        if (.not. local_ok) return
        candidate%levels(relation)%patches(parent)%state = state_work
        candidate%levels(relation)%patches(parent)%temperature = &
          temperature_work
        deallocate(state_work, temperature_work)
      end do
    end do
    ok = candidate%is_valid()
  end subroutine synchronize_candidate

  subroutine retain_same_resolution_overlap( &
      new_node, new_geometry, old_node, old_geometry, copied, ok)
    type(reactive_amr_eb_patch_tree_node_2d), intent(inout) :: new_node
    type(eb_geometry_2d), intent(in) :: new_geometry
    type(reactive_amr_eb_patch_tree_node_2d), intent(in) :: old_node
    type(eb_geometry_2d), intent(in) :: old_geometry
    logical, intent(inout) :: copied(:, :)
    logical, intent(out) :: ok

    real(dp) :: offset_i_real, offset_j_real, spacing_scale
    integer :: i, j, offset_i, offset_j, old_i, old_j

    ok = .true.
    spacing_scale = max(1.0_dp, abs(new_geometry%dx), &
      abs(new_geometry%dy), abs(old_geometry%dx), abs(old_geometry%dy))
    if (abs(new_geometry%dx - old_geometry%dx) > &
          geometry_tolerance * spacing_scale .or. &
        abs(new_geometry%dy - old_geometry%dy) > &
          geometry_tolerance * spacing_scale) return
    offset_i_real = (new_geometry%x_lower - old_geometry%x_lower) / &
      old_geometry%dx
    offset_j_real = (new_geometry%y_lower - old_geometry%y_lower) / &
      old_geometry%dy
    offset_i = nint(offset_i_real)
    offset_j = nint(offset_j_real)
    if (abs(offset_i_real - real(offset_i, dp)) > geometry_tolerance .or. &
        abs(offset_j_real - real(offset_j, dp)) > geometry_tolerance) return

    do j = 1, new_geometry%ny
      old_j = j + offset_j
      if (old_j < 1 .or. old_j > old_geometry%ny) cycle
      do i = 1, new_geometry%nx
        if (copied(i, j)) cycle
        old_i = i + offset_i
        if (old_i < 1 .or. old_i > old_geometry%nx) cycle
        if (.not. overlap_cell_geometry_matches( &
            new_geometry, i, j, old_geometry, old_i, old_j)) then
          ok = .false.
          return
        end if
        new_node%state(:, i, j) = old_node%state(:, old_i, old_j)
        new_node%temperature(i, j) = old_node%temperature(old_i, old_j)
        copied(i, j) = .true.
      end do
    end do
  end subroutine retain_same_resolution_overlap

  pure logical function overlap_cell_geometry_matches( &
      first, first_i, first_j, second, second_i, second_j) result(matches)
    type(eb_geometry_2d), intent(in) :: first, second
    integer, intent(in) :: first_i, first_j, second_i, second_j

    matches = first%cell_type(first_i, first_j) == &
        second%cell_type(second_i, second_j) .and. &
      abs(first%volume_fraction(first_i, first_j) - &
        second%volume_fraction(second_i, second_j)) <= &
        geometry_tolerance .and. &
      abs(first%cell_centroid_x(first_i, first_j) - &
        second%cell_centroid_x(second_i, second_j)) <= &
        geometry_tolerance .and. &
      abs(first%cell_centroid_y(first_i, first_j) - &
        second%cell_centroid_y(second_i, second_j)) <= &
        geometry_tolerance .and. &
      abs(first%boundary_length(first_i, first_j) - &
        second%boundary_length(second_i, second_j)) <= &
        geometry_tolerance .and. &
      abs(first%boundary_centroid_x(first_i, first_j) - &
        second%boundary_centroid_x(second_i, second_j)) <= &
        geometry_tolerance .and. &
      abs(first%boundary_centroid_y(first_i, first_j) - &
        second%boundary_centroid_y(second_i, second_j)) <= &
        geometry_tolerance .and. &
      abs(first%boundary_normal_x(first_i, first_j) - &
        second%boundary_normal_x(second_i, second_j)) <= &
        geometry_tolerance .and. &
      abs(first%boundary_normal_y(first_i, first_j) - &
        second%boundary_normal_y(second_i, second_j)) <= &
        geometry_tolerance .and. &
      abs(first%boundary_normal_integral_x(first_i, first_j) - &
        second%boundary_normal_integral_x(second_i, second_j)) <= &
        geometry_tolerance .and. &
      abs(first%boundary_normal_integral_y(first_i, first_j) - &
        second%boundary_normal_integral_y(second_i, second_j)) <= &
        geometry_tolerance .and. &
      all(abs(first%x_face_fraction(first_i - 1:first_i, first_j) - &
        second%x_face_fraction(second_i - 1:second_i, second_j)) <= &
        geometry_tolerance) .and. &
      all(abs(first%x_face_centroid_y(first_i - 1:first_i, first_j) - &
        second%x_face_centroid_y(second_i - 1:second_i, second_j)) <= &
        geometry_tolerance) .and. &
      all(abs(first%y_face_fraction(first_i, first_j - 1:first_j) - &
        second%y_face_fraction(second_i, second_j - 1:second_j)) <= &
        geometry_tolerance) .and. &
      all(abs(first%y_face_centroid_x(first_i, first_j - 1:first_j) - &
        second%y_face_centroid_x(second_i, second_j - 1:second_j)) <= &
        geometry_tolerance)
  end function overlap_cell_geometry_matches

  subroutine recover_patch_temperature(species, node, geometry, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_amr_eb_patch_tree_node_2d), intent(inout) :: node
    type(eb_geometry_2d), intent(in) :: geometry
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: recovered_temperature, sound_speed
    logical :: local_ok
    integer :: i, j

    ok = .false.
    if (any(.not. ieee_is_finite(node%state)) .or. &
        any(.not. ieee_is_finite(node%temperature)) .or. &
        any(node%temperature <= 0.0_dp)) return
    allocate(primitive(reactive_nprim(size(species))))
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) == eb_covered_cell) cycle
        call reactive_conserved_to_primitive( &
          species, node%state(:, i, j), node%temperature(i, j), &
          primitive, recovered_temperature, sound_speed, local_ok)
        if (.not. local_ok) return
        node%temperature(i, j) = recovered_temperature
      end do
    end do
    ok = .true.
  end subroutine recover_patch_temperature

  subroutine patch_geometry_at( &
      topology, level_index, patch_index, geometry, ok)
    type(amr_eb_patch_tree_topology_2d), intent(in) :: topology
    integer, intent(in) :: level_index, patch_index
    type(eb_geometry_2d), intent(out) :: geometry
    logical, intent(out) :: ok

    ok = .false.
    if (level_index < 1 .or. &
        level_index > topology%level_count()) return
    if (patch_index < 1 .or. &
        patch_index > topology%level_patch_count(level_index - 1)) return
    if (level_index == 1) then
      geometry = topology%root_geometry
    else
      geometry = topology%relations(level_index - 1)% &
        children(patch_index)%geometry
    end if
    ok = geometry%is_valid()
  end subroutine patch_geometry_at

end module amr_eb_patch_tree_reactive_2d_mod
