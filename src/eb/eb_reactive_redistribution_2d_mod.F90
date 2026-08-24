module eb_reactive_redistribution_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use eb_geometry_2d_mod, only: &
    eb_geometry_2d, eb_covered_cell, eb_cut_cell, eb_regular_cell
  implicit none
  private

  public :: reactive_eb_flux_redistribute_2d
  public :: advance_reactive_eb_redistributed_2d
  public :: reactive_eb_weighted_state_redistribute_2d
  public :: advance_reactive_eb_state_redistributed_2d

contains

  subroutine reactive_eb_flux_redistribute_2d( &
      geometry, conservative_rhs, redistributed_rhs, ok)
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), intent(in) :: conservative_rhs(:, :, :)
    real(dp), intent(out) :: redistributed_rhs(:, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate(:, :, :)
    real(dp), allocatable :: neighborhood_rhs(:, :, :)
    real(dp), allocatable :: neighbor_volume_fraction(:, :)
    real(dp), allocatable :: excess(:)
    real(dp) :: kappa, total_volume_fraction
    integer :: i, j, ncomp

    redistributed_rhs = 0.0_dp
    ok = .false.
    if (.not. geometry%is_valid()) return
    ncomp = size(conservative_rhs, 1)
    if (ncomp < 1 .or. &
        size(conservative_rhs, 2) /= geometry%nx .or. &
        size(conservative_rhs, 3) /= geometry%ny .or. &
        any(shape(redistributed_rhs) /= shape(conservative_rhs)) .or. &
        any(.not. ieee_is_finite(conservative_rhs))) return

    allocate(candidate(ncomp, geometry%nx, geometry%ny))
    allocate(neighborhood_rhs(ncomp, geometry%nx, geometry%ny))
    allocate(neighbor_volume_fraction(geometry%nx, geometry%ny))
    allocate(excess(ncomp))
    candidate = 0.0_dp
    neighborhood_rhs = 0.0_dp
    neighbor_volume_fraction = 0.0_dp

    do j = 1, geometry%ny
      do i = 1, geometry%nx
        select case (geometry%cell_type(i, j))
        case (eb_covered_cell)
          cycle
        case (eb_regular_cell)
          candidate(:, i, j) = conservative_rhs(:, i, j)
        case (eb_cut_cell)
          kappa = geometry%volume_fraction(i, j)
          total_volume_fraction = kappa
          neighborhood_rhs(:, i, j) = &
            kappa * conservative_rhs(:, i, j)
          call accumulate_neighbors( &
            geometry, conservative_rhs, i, j, neighborhood_rhs(:, i, j), &
            total_volume_fraction, neighbor_volume_fraction(i, j))
          if (neighbor_volume_fraction(i, j) <= 0.0_dp .or. &
              total_volume_fraction <= kappa) return
          neighborhood_rhs(:, i, j) = &
            neighborhood_rhs(:, i, j) / total_volume_fraction
          candidate(:, i, j) = &
            kappa * conservative_rhs(:, i, j) + &
            (1.0_dp - kappa) * neighborhood_rhs(:, i, j)
        case default
          return
        end select
      end do
    end do

    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) /= eb_cut_cell) cycle
        kappa = geometry%volume_fraction(i, j)
        excess = kappa * (1.0_dp - kappa) * &
          (conservative_rhs(:, i, j) - neighborhood_rhs(:, i, j))
        call distribute_to_neighbors( &
          geometry, i, j, excess / neighbor_volume_fraction(i, j), candidate)
      end do
    end do
    if (any(.not. ieee_is_finite(candidate))) return

    redistributed_rhs = candidate
    ok = .true.
  end subroutine reactive_eb_flux_redistribute_2d

  subroutine accumulate_neighbors( &
      geometry, rhs, i, j, weighted_rhs, total_volume, neighbor_volume)
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), intent(in) :: rhs(:, :, :)
    integer, intent(in) :: i, j
    real(dp), intent(inout) :: weighted_rhs(:), total_volume
    real(dp), intent(out) :: neighbor_volume

    neighbor_volume = 0.0_dp
    if (i > 1 .and. geometry%x_face_fraction(i - 1, j) > 0.0_dp) &
      call accumulate_cell(i - 1, j)
    if (i < geometry%nx .and. geometry%x_face_fraction(i, j) > 0.0_dp) &
      call accumulate_cell(i + 1, j)
    if (j > 1 .and. geometry%y_face_fraction(i, j - 1) > 0.0_dp) &
      call accumulate_cell(i, j - 1)
    if (j < geometry%ny .and. geometry%y_face_fraction(i, j) > 0.0_dp) &
      call accumulate_cell(i, j + 1)

  contains

    subroutine accumulate_cell(neighbor_i, neighbor_j)
      integer, intent(in) :: neighbor_i, neighbor_j
      real(dp) :: neighbor_kappa

      if (geometry%cell_type(neighbor_i, neighbor_j) == &
          eb_covered_cell) return
      neighbor_kappa = geometry%volume_fraction(neighbor_i, neighbor_j)
      weighted_rhs = weighted_rhs + &
        neighbor_kappa * rhs(:, neighbor_i, neighbor_j)
      total_volume = total_volume + neighbor_kappa
      neighbor_volume = neighbor_volume + neighbor_kappa
    end subroutine accumulate_cell

  end subroutine accumulate_neighbors

  subroutine distribute_to_neighbors( &
      geometry, i, j, increment, redistributed_rhs)
    type(eb_geometry_2d), intent(in) :: geometry
    integer, intent(in) :: i, j
    real(dp), intent(in) :: increment(:)
    real(dp), intent(inout) :: redistributed_rhs(:, :, :)

    if (i > 1 .and. geometry%x_face_fraction(i - 1, j) > 0.0_dp) &
      call add_to_cell(i - 1, j)
    if (i < geometry%nx .and. geometry%x_face_fraction(i, j) > 0.0_dp) &
      call add_to_cell(i + 1, j)
    if (j > 1 .and. geometry%y_face_fraction(i, j - 1) > 0.0_dp) &
      call add_to_cell(i, j - 1)
    if (j < geometry%ny .and. geometry%y_face_fraction(i, j) > 0.0_dp) &
      call add_to_cell(i, j + 1)

  contains

    subroutine add_to_cell(neighbor_i, neighbor_j)
      integer, intent(in) :: neighbor_i, neighbor_j

      if (geometry%cell_type(neighbor_i, neighbor_j) == &
          eb_covered_cell) return
      redistributed_rhs(:, neighbor_i, neighbor_j) = &
        redistributed_rhs(:, neighbor_i, neighbor_j) + increment
    end subroutine add_to_cell

  end subroutine distribute_to_neighbors

  subroutine advance_reactive_eb_redistributed_2d( &
      species, state, temperature, geometry, conservative_rhs, dt, &
      new_state, new_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), intent(in) :: conservative_rhs(:, :, :), dt
    real(dp), intent(out) :: new_state(:, :, :), new_temperature(:, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate_state(:, :, :)
    real(dp), allocatable :: candidate_temperature(:, :)
    real(dp), allocatable :: redistributed_rhs(:, :, :), primitive(:)
    real(dp) :: recovered_temperature, sound_speed
    logical :: local_ok
    integer :: i, j, nvar

    new_state = 0.0_dp
    new_temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        any(shape(temperature) /= [geometry%nx, geometry%ny]) .or. &
        any(shape(conservative_rhs) /= shape(state)) .or. &
        any(shape(new_state) /= shape(state)) .or. &
        any(shape(new_temperature) /= shape(temperature))) return
    new_state = state
    new_temperature = temperature
    if (.not. ieee_is_finite(dt) .or. dt < 0.0_dp) return

    allocate(redistributed_rhs(nvar, geometry%nx, geometry%ny))
    call reactive_eb_flux_redistribute_2d( &
      geometry, conservative_rhs, redistributed_rhs, local_ok)
    if (.not. local_ok) return
    allocate(candidate_state(nvar, geometry%nx, geometry%ny))
    allocate(candidate_temperature(geometry%nx, geometry%ny))
    allocate(primitive(reactive_nprim(size(species))))
    candidate_state = state
    candidate_temperature = temperature
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) == eb_covered_cell) cycle
        if (.not. ieee_is_finite(temperature(i, j)) .or. &
            temperature(i, j) <= 0.0_dp) return
        candidate_state(:, i, j) = state(:, i, j) + &
          dt * redistributed_rhs(:, i, j)
        call reactive_conserved_to_primitive( &
          species, candidate_state(:, i, j), temperature(i, j), primitive, &
          recovered_temperature, sound_speed, local_ok)
        if (.not. local_ok) return
        candidate_temperature(i, j) = recovered_temperature
      end do
    end do

    new_state = candidate_state
    new_temperature = candidate_temperature
    ok = .true.
  end subroutine advance_reactive_eb_redistributed_2d

  subroutine reactive_eb_weighted_state_redistribute_2d( &
      geometry, provisional_state, redistributed_state, ok, &
      target_volume_fraction)
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), intent(in) :: provisional_state(:, :, :)
    real(dp), intent(out) :: redistributed_state(:, :, :)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: target_volume_fraction

    integer, allocatable :: neighbor_count(:, :)
    integer, allocatable :: neighbor_offset_i(:, :, :)
    integer, allocatable :: neighbor_offset_j(:, :, :)
    integer, allocatable :: neighborhood_count(:, :)
    real(dp), allocatable :: alpha_self(:, :), alpha_neighbor(:, :)
    real(dp), allocatable :: neighborhood_volume(:, :)
    real(dp), allocatable :: neighborhood_state(:, :, :)
    real(dp), allocatable :: candidate(:, :, :)
    real(dp) :: target
    integer :: i, j, neighbor, neighbor_i, neighbor_j, ncomp

    redistributed_state = 0.0_dp
    ok = .false.
    if (.not. geometry%is_valid()) return
    ncomp = size(provisional_state, 1)
    if (ncomp < 1 .or. &
        size(provisional_state, 2) /= geometry%nx .or. &
        size(provisional_state, 3) /= geometry%ny .or. &
        any(shape(redistributed_state) /= shape(provisional_state)) .or. &
        any(.not. ieee_is_finite(provisional_state))) return

    target = 0.5_dp
    if (present(target_volume_fraction)) target = target_volume_fraction
    if (.not. ieee_is_finite(target) .or. &
        target <= 0.0_dp .or. target > 1.0_dp) return

    allocate(neighbor_count(geometry%nx, geometry%ny))
    allocate(neighbor_offset_i(3, geometry%nx, geometry%ny))
    allocate(neighbor_offset_j(3, geometry%nx, geometry%ny))
    allocate(neighborhood_count(geometry%nx, geometry%ny))
    allocate(alpha_self(geometry%nx, geometry%ny))
    allocate(alpha_neighbor(geometry%nx, geometry%ny))
    allocate(neighborhood_volume(geometry%nx, geometry%ny))
    ! Match AMReX StateRedist with max_order=0: build the overlapping
    ! neighborhoods and their partition before averaging the provisional state.
    call build_state_redistribution_neighborhoods( &
      geometry, target, neighbor_count, neighbor_offset_i, &
      neighbor_offset_j, neighborhood_count, alpha_self, alpha_neighbor, &
      neighborhood_volume, ok)
    if (.not. ok) return

    allocate(neighborhood_state(ncomp, geometry%nx, geometry%ny))
    allocate(candidate(ncomp, geometry%nx, geometry%ny))
    neighborhood_state = 0.0_dp
    candidate = 0.0_dp
    ! Form each zeroth-order neighborhood state Qhat.
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) == eb_covered_cell) cycle
        neighborhood_state(:, i, j) = alpha_self(i, j) * &
          geometry%volume_fraction(i, j) * provisional_state(:, i, j)
        do neighbor = 1, neighbor_count(i, j)
          neighbor_i = i + neighbor_offset_i(neighbor, i, j)
          neighbor_j = j + neighbor_offset_j(neighbor, i, j)
          neighborhood_state(:, i, j) = neighborhood_state(:, i, j) + &
            alpha_neighbor(i, j) * &
            geometry%volume_fraction(neighbor_i, neighbor_j) * &
            provisional_state(:, neighbor_i, neighbor_j) / &
            real(neighborhood_count(neighbor_i, neighbor_j), dp)
        end do
        neighborhood_state(:, i, j) = neighborhood_state(:, i, j) / &
          neighborhood_volume(i, j)
      end do
    end do

    ! Scatter every Qhat through the same self/neighbor partition.
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) == eb_covered_cell) cycle
        candidate(:, i, j) = candidate(:, i, j) + &
          alpha_self(i, j) * real(neighborhood_count(i, j), dp) * &
          neighborhood_state(:, i, j)
        do neighbor = 1, neighbor_count(i, j)
          neighbor_i = i + neighbor_offset_i(neighbor, i, j)
          neighbor_j = j + neighbor_offset_j(neighbor, i, j)
          candidate(:, neighbor_i, neighbor_j) = &
            candidate(:, neighbor_i, neighbor_j) + &
            alpha_neighbor(i, j) * neighborhood_state(:, i, j)
        end do
      end do
    end do
    ! Complete the overlapping-neighborhood average at every recipient.
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) == eb_covered_cell) cycle
        candidate(:, i, j) = candidate(:, i, j) / &
          real(neighborhood_count(i, j), dp)
      end do
    end do
    if (any(.not. ieee_is_finite(candidate))) then
      ok = .false.
      return
    end if

    redistributed_state = candidate
    ok = .true.
  end subroutine reactive_eb_weighted_state_redistribute_2d

  subroutine build_state_redistribution_neighborhoods( &
      geometry, target, neighbor_count, neighbor_offset_i, &
      neighbor_offset_j, neighborhood_count, alpha_self, alpha_neighbor, &
      neighborhood_volume, ok)
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), intent(in) :: target
    integer, intent(out) :: neighbor_count(:, :)
    integer, intent(out) :: neighbor_offset_i(:, :, :)
    integer, intent(out) :: neighbor_offset_j(:, :, :)
    integer, intent(out) :: neighborhood_count(:, :)
    real(dp), intent(out) :: alpha_self(:, :), alpha_neighbor(:, :)
    real(dp), intent(out) :: neighborhood_volume(:, :)
    logical, intent(out) :: ok

    real(dp), parameter :: equal_normal_tolerance = 1.0e-8_dp
    real(dp), parameter :: weight_tolerance = &
      1024.0_dp * epsilon(1.0_dp)
    real(dp) :: normal_x, normal_y, normal_norm, sum_volume
    real(dp) :: neighbor_volume
    logical :: equal_normal_components, local_ok
    integer :: i, j, neighbor, offset_i, offset_j

    neighbor_count = 0
    neighbor_offset_i = 0
    neighbor_offset_j = 0
    neighborhood_count = 1
    alpha_self = 0.0_dp
    alpha_neighbor = 0.0_dp
    neighborhood_volume = 0.0_dp
    ok = .false.

    ! Build each small-cell merge list from the aperture-difference normal.
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) == eb_covered_cell) cycle
        alpha_self(i, j) = 1.0_dp
        alpha_neighbor(i, j) = 1.0_dp
        if (geometry%volume_fraction(i, j) >= target) cycle

        normal_x = geometry%x_face_fraction(i, j) - &
          geometry%x_face_fraction(i - 1, j)
        normal_y = geometry%y_face_fraction(i, j) - &
          geometry%y_face_fraction(i, j - 1)
        normal_norm = sqrt(normal_x**2 + normal_y**2)
        if (.not. ieee_is_finite(normal_norm) .or. &
            normal_norm <= tiny(1.0_dp)) return
        normal_x = normal_x / normal_norm
        normal_y = normal_y / normal_norm
        equal_normal_components = &
          abs(normal_x - normal_y) < equal_normal_tolerance .or. &
          abs(normal_x + normal_y) < equal_normal_tolerance

        if (abs(normal_x) > abs(normal_y)) then
          offset_i = merge(1, -1, normal_x > 0.0_dp)
          offset_j = 0
        else
          offset_i = 0
          offset_j = merge(1, -1, normal_y > 0.0_dp)
        end if
        if (i + offset_i < 1 .or. i + offset_i > geometry%nx) then
          offset_i = 0
          offset_j = merge(1, -1, normal_y > 0.0_dp)
        end if
        if (j + offset_j < 1 .or. j + offset_j > geometry%ny) then
          offset_i = merge(1, -1, normal_x > 0.0_dp)
          offset_j = 0
        end if
        sum_volume = geometry%volume_fraction(i, j)
        call append_state_redistribution_neighbor( &
          geometry, i, j, offset_i, offset_j, neighbor_count, &
          neighbor_offset_i, neighbor_offset_j, sum_volume, local_ok)
        if (.not. local_ok) return

        if (sum_volume < target .or. equal_normal_components) then
          if (offset_i == 0) then
            if (normal_x >= 0.0_dp .and. i < geometry%nx) then
              call append_state_redistribution_neighbor( &
                geometry, i, j, 1, 0, neighbor_count, &
                neighbor_offset_i, neighbor_offset_j, sum_volume, local_ok)
            else if (normal_x <= 0.0_dp .and. i > 1) then
              call append_state_redistribution_neighbor( &
                geometry, i, j, -1, 0, neighbor_count, &
                neighbor_offset_i, neighbor_offset_j, sum_volume, local_ok)
            else
              local_ok = .true.
            end if
          else
            if (normal_y >= 0.0_dp .and. j < geometry%ny) then
              call append_state_redistribution_neighbor( &
                geometry, i, j, 0, 1, neighbor_count, &
                neighbor_offset_i, neighbor_offset_j, sum_volume, local_ok)
            else if (normal_y <= 0.0_dp .and. j > 1) then
              call append_state_redistribution_neighbor( &
                geometry, i, j, 0, -1, neighbor_count, &
                neighbor_offset_i, neighbor_offset_j, sum_volume, local_ok)
            else
              local_ok = .true.
            end if
          end if
          if (.not. local_ok) return
        end if

        if (neighbor_count(i, j) == 2) then
          offset_i = neighbor_offset_i(1, i, j) + &
            neighbor_offset_i(2, i, j)
          offset_j = neighbor_offset_j(1, i, j) + &
            neighbor_offset_j(2, i, j)
          call append_state_redistribution_neighbor( &
            geometry, i, j, offset_i, offset_j, neighbor_count, &
            neighbor_offset_i, neighbor_offset_j, sum_volume, local_ok)
          if (.not. local_ok) return
        end if
        if (sum_volume + weight_tolerance < target) return
      end do
    end do

    ! nrs: every active cell belongs to its own neighborhood, plus each
    ! small-cell merge neighborhood that names it.
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        do neighbor = 1, neighbor_count(i, j)
          offset_i = i + neighbor_offset_i(neighbor, i, j)
          offset_j = j + neighbor_offset_j(neighbor, i, j)
          neighborhood_count(offset_i, offset_j) = &
            neighborhood_count(offset_i, offset_j) + 1
        end do
      end do
    end do

    ! The neighbor weight raises each small neighborhood to the target volume.
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) == eb_covered_cell) cycle
        if (neighbor_count(i, j) == 0) cycle
        neighbor_volume = 0.0_dp
        do neighbor = 1, neighbor_count(i, j)
          offset_i = i + neighbor_offset_i(neighbor, i, j)
          offset_j = j + neighbor_offset_j(neighbor, i, j)
          neighbor_volume = neighbor_volume + &
            geometry%volume_fraction(offset_i, offset_j)
        end do
        if (neighbor_volume <= 0.0_dp) return
        alpha_neighbor(i, j) = &
          (target - geometry%volume_fraction(i, j)) / neighbor_volume
      end do
    end do

    ! Subtract the shared neighbor partitions from each cell's self partition.
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        do neighbor = 1, neighbor_count(i, j)
          offset_i = i + neighbor_offset_i(neighbor, i, j)
          offset_j = j + neighbor_offset_j(neighbor, i, j)
          alpha_self(offset_i, offset_j) = &
            alpha_self(offset_i, offset_j) - alpha_neighbor(i, j) / &
            real(neighborhood_count(offset_i, offset_j), dp)
        end do
      end do
    end do

    ! The same partition defines the neighborhood volume used by Qhat.
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) == eb_covered_cell) cycle
        if (.not. ieee_is_finite(alpha_self(i, j)) .or. &
            alpha_self(i, j) < -weight_tolerance .or. &
            alpha_neighbor(i, j) < 0.0_dp) return
        neighborhood_volume(i, j) = alpha_self(i, j) * &
          geometry%volume_fraction(i, j)
        do neighbor = 1, neighbor_count(i, j)
          offset_i = i + neighbor_offset_i(neighbor, i, j)
          offset_j = j + neighbor_offset_j(neighbor, i, j)
          neighborhood_volume(i, j) = neighborhood_volume(i, j) + &
            alpha_neighbor(i, j) * &
            geometry%volume_fraction(offset_i, offset_j) / &
            real(neighborhood_count(offset_i, offset_j), dp)
        end do
        if (.not. ieee_is_finite(neighborhood_volume(i, j)) .or. &
            neighborhood_volume(i, j) <= 0.0_dp) return
      end do
    end do
    ok = .true.
  end subroutine build_state_redistribution_neighborhoods

  subroutine append_state_redistribution_neighbor( &
      geometry, i, j, offset_i, offset_j, neighbor_count, &
      neighbor_offset_i, neighbor_offset_j, sum_volume, ok)
    type(eb_geometry_2d), intent(in) :: geometry
    integer, intent(in) :: i, j, offset_i, offset_j
    integer, intent(inout) :: neighbor_count(:, :)
    integer, intent(inout) :: neighbor_offset_i(:, :, :)
    integer, intent(inout) :: neighbor_offset_j(:, :, :)
    real(dp), intent(inout) :: sum_volume
    logical, intent(out) :: ok

    integer :: count, neighbor_i, neighbor_j

    ok = .false.
    neighbor_i = i + offset_i
    neighbor_j = j + offset_j
    if ((offset_i == 0 .and. offset_j == 0) .or. &
        neighbor_i < 1 .or. neighbor_i > geometry%nx .or. &
        neighbor_j < 1 .or. neighbor_j > geometry%ny .or. &
        geometry%cell_type(neighbor_i, neighbor_j) == eb_covered_cell) return
    count = neighbor_count(i, j) + 1
    if (count > size(neighbor_offset_i, 1)) return
    neighbor_count(i, j) = count
    neighbor_offset_i(count, i, j) = offset_i
    neighbor_offset_j(count, i, j) = offset_j
    sum_volume = sum_volume + &
      geometry%volume_fraction(neighbor_i, neighbor_j)
    ok = .true.
  end subroutine append_state_redistribution_neighbor

  subroutine advance_reactive_eb_state_redistributed_2d( &
      species, state, temperature, geometry, conservative_rhs, dt, &
      new_state, new_temperature, ok, target_volume_fraction)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), intent(in) :: conservative_rhs(:, :, :), dt
    real(dp), intent(out) :: new_state(:, :, :), new_temperature(:, :)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: target_volume_fraction

    real(dp), allocatable :: provisional_state(:, :, :)
    real(dp), allocatable :: redistributed_state(:, :, :)
    real(dp), allocatable :: candidate_state(:, :, :)
    real(dp), allocatable :: candidate_temperature(:, :), primitive(:)
    real(dp) :: recovered_temperature, sound_speed
    logical :: local_ok
    integer :: i, j, nvar

    new_state = 0.0_dp
    new_temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        any(shape(temperature) /= [geometry%nx, geometry%ny]) .or. &
        any(shape(conservative_rhs) /= shape(state)) .or. &
        any(shape(new_state) /= shape(state)) .or. &
        any(shape(new_temperature) /= shape(temperature))) return
    new_state = state
    new_temperature = temperature
    if (.not. ieee_is_finite(dt) .or. dt < 0.0_dp .or. &
        any(.not. ieee_is_finite(state)) .or. &
        any(.not. ieee_is_finite(temperature)) .or. &
        any(.not. ieee_is_finite(conservative_rhs))) return

    allocate(provisional_state(nvar, geometry%nx, geometry%ny))
    allocate(redistributed_state(nvar, geometry%nx, geometry%ny))
    provisional_state = state
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) == eb_covered_cell) cycle
        provisional_state(:, i, j) = state(:, i, j) + &
          dt * conservative_rhs(:, i, j)
      end do
    end do
    if (present(target_volume_fraction)) then
      call reactive_eb_weighted_state_redistribute_2d( &
        geometry, provisional_state, redistributed_state, local_ok, &
        target_volume_fraction)
    else
      call reactive_eb_weighted_state_redistribute_2d( &
        geometry, provisional_state, redistributed_state, local_ok)
    end if
    if (.not. local_ok) return

    allocate(candidate_state(nvar, geometry%nx, geometry%ny))
    allocate(candidate_temperature(geometry%nx, geometry%ny))
    allocate(primitive(reactive_nprim(size(species))))
    candidate_state = state
    candidate_temperature = temperature
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) == eb_covered_cell) cycle
        if (temperature(i, j) <= 0.0_dp) return
        candidate_state(:, i, j) = redistributed_state(:, i, j)
        call reactive_conserved_to_primitive( &
          species, candidate_state(:, i, j), temperature(i, j), primitive, &
          recovered_temperature, sound_speed, local_ok)
        if (.not. local_ok) return
        candidate_temperature(i, j) = recovered_temperature
      end do
    end do

    new_state = candidate_state
    new_temperature = candidate_temperature
    ok = .true.
  end subroutine advance_reactive_eb_state_redistributed_2d

end module eb_reactive_redistribution_2d_mod
