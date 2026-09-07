module eb_reactive_redistribution_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d, eb_cut_cell_3d, &
    eb_regular_cell_3d
  implicit none
  private

  public :: reactive_eb_flux_redistribute_3d
  public :: advance_reactive_eb_redistributed_3d
  public :: reactive_eb_weighted_state_redistribute_3d
  public :: advance_reactive_eb_state_redistributed_3d

contains

  subroutine reactive_eb_flux_redistribute_3d( &
      geometry, conservative_rhs, redistributed_rhs, ok)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: conservative_rhs(:, :, :, :)
    real(dp), intent(out) :: redistributed_rhs(:, :, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate(:, :, :, :)
    real(dp), allocatable :: neighborhood_rhs(:, :, :, :)
    real(dp), allocatable :: neighbor_volume_fraction(:, :, :)
    real(dp), allocatable :: excess(:)
    real(dp) :: kappa, total_volume_fraction
    integer :: i, j, k, ncomp

    redistributed_rhs = 0.0_dp
    ok = .false.
    if (.not. geometry%is_valid()) return
    ncomp = size(conservative_rhs, 1)
    if (ncomp < 1 .or. &
        size(conservative_rhs, 2) /= geometry%nx .or. &
        size(conservative_rhs, 3) /= geometry%ny .or. &
        size(conservative_rhs, 4) /= geometry%nz .or. &
        any(shape(redistributed_rhs) /= shape(conservative_rhs)) .or. &
        any(.not. ieee_is_finite(conservative_rhs))) return

    allocate(candidate(ncomp, geometry%nx, geometry%ny, geometry%nz))
    allocate(neighborhood_rhs( &
      ncomp, geometry%nx, geometry%ny, geometry%nz))
    allocate(neighbor_volume_fraction( &
      geometry%nx, geometry%ny, geometry%nz))
    allocate(excess(ncomp))
    candidate = 0.0_dp
    neighborhood_rhs = 0.0_dp
    neighbor_volume_fraction = 0.0_dp

    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          select case (geometry%cell_type(i, j, k))
          case (eb_covered_cell_3d)
            cycle
          case (eb_regular_cell_3d)
            candidate(:, i, j, k) = conservative_rhs(:, i, j, k)
          case (eb_cut_cell_3d)
            kappa = geometry%volume_fraction(i, j, k)
            total_volume_fraction = kappa
            neighborhood_rhs(:, i, j, k) = &
              kappa * conservative_rhs(:, i, j, k)
            call accumulate_neighbors( &
              geometry, conservative_rhs, i, j, k, &
              neighborhood_rhs(:, i, j, k), total_volume_fraction, &
              neighbor_volume_fraction(i, j, k))
            if (neighbor_volume_fraction(i, j, k) <= 0.0_dp .or. &
                total_volume_fraction <= kappa) return
            neighborhood_rhs(:, i, j, k) = &
              neighborhood_rhs(:, i, j, k) / total_volume_fraction
            candidate(:, i, j, k) = &
              kappa * conservative_rhs(:, i, j, k) + &
              (1.0_dp - kappa) * neighborhood_rhs(:, i, j, k)
          case default
            return
          end select
        end do
      end do
    end do

    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) /= eb_cut_cell_3d) cycle
          kappa = geometry%volume_fraction(i, j, k)
          excess = kappa * (1.0_dp - kappa) * &
            (conservative_rhs(:, i, j, k) - &
              neighborhood_rhs(:, i, j, k))
          call distribute_to_neighbors( &
            geometry, i, j, k, &
            excess / neighbor_volume_fraction(i, j, k), candidate)
        end do
      end do
    end do
    if (any(.not. ieee_is_finite(candidate))) return

    redistributed_rhs = candidate
    ok = .true.
  end subroutine reactive_eb_flux_redistribute_3d

  subroutine accumulate_neighbors( &
      geometry, rhs, i, j, k, weighted_rhs, total_volume, &
      neighbor_volume)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: rhs(:, :, :, :)
    integer, intent(in) :: i, j, k
    real(dp), intent(inout) :: weighted_rhs(:), total_volume
    real(dp), intent(out) :: neighbor_volume

    neighbor_volume = 0.0_dp
    if (i > 1) then
      if (geometry%x_face_fraction(i - 1, j, k) > 0.0_dp) &
        call accumulate_cell(i - 1, j, k)
    end if
    if (i < geometry%nx) then
      if (geometry%x_face_fraction(i, j, k) > 0.0_dp) &
        call accumulate_cell(i + 1, j, k)
    end if
    if (j > 1) then
      if (geometry%y_face_fraction(i, j - 1, k) > 0.0_dp) &
        call accumulate_cell(i, j - 1, k)
    end if
    if (j < geometry%ny) then
      if (geometry%y_face_fraction(i, j, k) > 0.0_dp) &
        call accumulate_cell(i, j + 1, k)
    end if
    if (k > 1) then
      if (geometry%z_face_fraction(i, j, k - 1) > 0.0_dp) &
        call accumulate_cell(i, j, k - 1)
    end if
    if (k < geometry%nz) then
      if (geometry%z_face_fraction(i, j, k) > 0.0_dp) &
        call accumulate_cell(i, j, k + 1)
    end if

  contains

    subroutine accumulate_cell(neighbor_i, neighbor_j, neighbor_k)
      integer, intent(in) :: neighbor_i, neighbor_j, neighbor_k
      real(dp) :: neighbor_kappa

      if (geometry%cell_type(neighbor_i, neighbor_j, neighbor_k) == &
          eb_covered_cell_3d) return
      neighbor_kappa = &
        geometry%volume_fraction(neighbor_i, neighbor_j, neighbor_k)
      weighted_rhs = weighted_rhs + neighbor_kappa * &
        rhs(:, neighbor_i, neighbor_j, neighbor_k)
      total_volume = total_volume + neighbor_kappa
      neighbor_volume = neighbor_volume + neighbor_kappa
    end subroutine accumulate_cell

  end subroutine accumulate_neighbors

  subroutine distribute_to_neighbors( &
      geometry, i, j, k, increment, redistributed_rhs)
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: i, j, k
    real(dp), intent(in) :: increment(:)
    real(dp), intent(inout) :: redistributed_rhs(:, :, :, :)

    if (i > 1) then
      if (geometry%x_face_fraction(i - 1, j, k) > 0.0_dp) &
        call add_to_cell(i - 1, j, k)
    end if
    if (i < geometry%nx) then
      if (geometry%x_face_fraction(i, j, k) > 0.0_dp) &
        call add_to_cell(i + 1, j, k)
    end if
    if (j > 1) then
      if (geometry%y_face_fraction(i, j - 1, k) > 0.0_dp) &
        call add_to_cell(i, j - 1, k)
    end if
    if (j < geometry%ny) then
      if (geometry%y_face_fraction(i, j, k) > 0.0_dp) &
        call add_to_cell(i, j + 1, k)
    end if
    if (k > 1) then
      if (geometry%z_face_fraction(i, j, k - 1) > 0.0_dp) &
        call add_to_cell(i, j, k - 1)
    end if
    if (k < geometry%nz) then
      if (geometry%z_face_fraction(i, j, k) > 0.0_dp) &
        call add_to_cell(i, j, k + 1)
    end if

  contains

    subroutine add_to_cell(neighbor_i, neighbor_j, neighbor_k)
      integer, intent(in) :: neighbor_i, neighbor_j, neighbor_k

      if (geometry%cell_type(neighbor_i, neighbor_j, neighbor_k) == &
          eb_covered_cell_3d) return
      redistributed_rhs(:, neighbor_i, neighbor_j, neighbor_k) = &
        redistributed_rhs(:, neighbor_i, neighbor_j, neighbor_k) + &
        increment
    end subroutine add_to_cell

  end subroutine distribute_to_neighbors

  subroutine advance_reactive_eb_redistributed_3d( &
      species, state, temperature, geometry, conservative_rhs, dt, &
      new_state, new_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: conservative_rhs(:, :, :, :), dt
    real(dp), intent(out) :: new_state(:, :, :, :)
    real(dp), intent(out) :: new_temperature(:, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    real(dp), allocatable :: redistributed_rhs(:, :, :, :)
    real(dp), allocatable :: primitive(:)
    real(dp) :: recovered_temperature, sound_speed
    logical :: local_ok
    integer :: i, j, k, nvar

    new_state = 0.0_dp
    new_temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        size(state, 4) /= geometry%nz .or. &
        any(shape(temperature) /= &
          [geometry%nx, geometry%ny, geometry%nz]) .or. &
        any(shape(conservative_rhs) /= shape(state)) .or. &
        any(shape(new_state) /= shape(state)) .or. &
        any(shape(new_temperature) /= shape(temperature))) return
    new_state = state
    new_temperature = temperature
    if (.not. ieee_is_finite(dt) .or. dt < 0.0_dp .or. &
        any(.not. ieee_is_finite(state)) .or. &
        any(.not. ieee_is_finite(temperature)) .or. &
        any(.not. ieee_is_finite(conservative_rhs))) return

    allocate(redistributed_rhs(nvar, geometry%nx, geometry%ny, geometry%nz))
    call reactive_eb_flux_redistribute_3d( &
      geometry, conservative_rhs, redistributed_rhs, local_ok)
    if (.not. local_ok) return

    allocate(candidate_state, source=state)
    allocate(candidate_temperature, source=temperature)
    allocate(primitive(reactive_nprim(size(species))))
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          if (temperature(i, j, k) <= 0.0_dp) return
          candidate_state(:, i, j, k) = state(:, i, j, k) + &
            dt * redistributed_rhs(:, i, j, k)
          call reactive_conserved_to_primitive( &
            species, candidate_state(:, i, j, k), temperature(i, j, k), &
            primitive, recovered_temperature, sound_speed, local_ok)
          if (.not. local_ok) return
          candidate_temperature(i, j, k) = recovered_temperature
        end do
      end do
    end do
    if (any(.not. ieee_is_finite(candidate_state)) .or. &
        any(.not. ieee_is_finite(candidate_temperature))) return

    new_state = candidate_state
    new_temperature = candidate_temperature
    ok = .true.
  end subroutine advance_reactive_eb_redistributed_3d

  subroutine reactive_eb_weighted_state_redistribute_3d( &
      geometry, provisional_state, redistributed_state, ok, &
      target_volume_fraction)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: provisional_state(:, :, :, :)
    real(dp), intent(out) :: redistributed_state(:, :, :, :)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: target_volume_fraction

    logical, allocatable :: has_neighbor(:, :, :)
    integer, allocatable :: neighbor_i(:, :, :)
    integer, allocatable :: neighbor_j(:, :, :)
    integer, allocatable :: neighbor_k(:, :, :)
    integer, allocatable :: neighborhood_count(:, :, :)
    real(dp), allocatable :: alpha_self(:, :, :)
    real(dp), allocatable :: alpha_neighbor(:, :, :)
    real(dp), allocatable :: neighborhood_volume(:, :, :)
    real(dp), allocatable :: neighborhood_state(:, :, :, :)
    real(dp), allocatable :: candidate(:, :, :, :)
    real(dp), parameter :: weight_tolerance = &
      1024.0_dp * epsilon(1.0_dp)
    real(dp) :: target, kappa, receiving_kappa
    logical :: local_ok
    integer :: i, j, k, ni, nj, nk, ncomp

    redistributed_state = 0.0_dp
    ok = .false.
    if (.not. geometry%is_valid()) return
    ncomp = size(provisional_state, 1)
    if (ncomp < 1 .or. &
        size(provisional_state, 2) /= geometry%nx .or. &
        size(provisional_state, 3) /= geometry%ny .or. &
        size(provisional_state, 4) /= geometry%nz .or. &
        any(shape(redistributed_state) /= shape(provisional_state)) .or. &
        any(.not. ieee_is_finite(provisional_state))) return

    target = 0.5_dp
    if (present(target_volume_fraction)) target = target_volume_fraction
    if (.not. ieee_is_finite(target) .or. &
        target <= 0.0_dp .or. target > 1.0_dp) return

    allocate(has_neighbor(geometry%nx, geometry%ny, geometry%nz))
    allocate(neighbor_i(geometry%nx, geometry%ny, geometry%nz))
    allocate(neighbor_j(geometry%nx, geometry%ny, geometry%nz))
    allocate(neighbor_k(geometry%nx, geometry%ny, geometry%nz))
    allocate(neighborhood_count(geometry%nx, geometry%ny, geometry%nz))
    allocate(alpha_self(geometry%nx, geometry%ny, geometry%nz))
    allocate(alpha_neighbor(geometry%nx, geometry%ny, geometry%nz))
    allocate(neighborhood_volume(geometry%nx, geometry%ny, geometry%nz))
    has_neighbor = .false.
    neighbor_i = 0
    neighbor_j = 0
    neighbor_k = 0
    neighborhood_count = 1
    alpha_self = 0.0_dp
    alpha_neighbor = 0.0_dp
    neighborhood_volume = 0.0_dp

    ! The qualified geometry has one axis-normal regular receiver for every
    ! small cut cell. Keep the weighting general enough for shared receivers.
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          alpha_self(i, j, k) = 1.0_dp
          kappa = geometry%volume_fraction(i, j, k)
          if (kappa >= target) cycle
          call select_axis_plane_state_receiver( &
            geometry, i, j, k, ni, nj, nk, local_ok)
          if (.not. local_ok) return
          receiving_kappa = geometry%volume_fraction(ni, nj, nk)
          if (kappa + receiving_kappa + weight_tolerance < target) return
          has_neighbor(i, j, k) = .true.
          neighbor_i(i, j, k) = ni
          neighbor_j(i, j, k) = nj
          neighbor_k(i, j, k) = nk
          alpha_neighbor(i, j, k) = &
            (target - kappa) / receiving_kappa
          neighborhood_count(ni, nj, nk) = &
            neighborhood_count(ni, nj, nk) + 1
        end do
      end do
    end do

    ! Each receiving cell reserves an equal share for every neighborhood to
    ! which it belongs, matching the zeroth-order weighted StateRedist rule.
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (.not. has_neighbor(i, j, k)) cycle
          ni = neighbor_i(i, j, k)
          nj = neighbor_j(i, j, k)
          nk = neighbor_k(i, j, k)
          alpha_self(ni, nj, nk) = alpha_self(ni, nj, nk) - &
            alpha_neighbor(i, j, k) / &
              real(neighborhood_count(ni, nj, nk), dp)
        end do
      end do
    end do

    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          if (.not. ieee_is_finite(alpha_self(i, j, k)) .or. &
              alpha_self(i, j, k) < -weight_tolerance .or. &
              alpha_neighbor(i, j, k) < 0.0_dp) return
          neighborhood_volume(i, j, k) = alpha_self(i, j, k) * &
            geometry%volume_fraction(i, j, k)
          if (has_neighbor(i, j, k)) then
            ni = neighbor_i(i, j, k)
            nj = neighbor_j(i, j, k)
            nk = neighbor_k(i, j, k)
            neighborhood_volume(i, j, k) = &
              neighborhood_volume(i, j, k) + &
              alpha_neighbor(i, j, k) * &
                geometry%volume_fraction(ni, nj, nk) / &
                real(neighborhood_count(ni, nj, nk), dp)
          end if
          if (.not. ieee_is_finite(neighborhood_volume(i, j, k)) .or. &
              neighborhood_volume(i, j, k) <= 0.0_dp) return
        end do
      end do
    end do

    allocate(neighborhood_state( &
      ncomp, geometry%nx, geometry%ny, geometry%nz))
    allocate(candidate(ncomp, geometry%nx, geometry%ny, geometry%nz))
    neighborhood_state = 0.0_dp
    candidate = 0.0_dp
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          neighborhood_state(:, i, j, k) = alpha_self(i, j, k) * &
            geometry%volume_fraction(i, j, k) * &
            provisional_state(:, i, j, k)
          if (has_neighbor(i, j, k)) then
            ni = neighbor_i(i, j, k)
            nj = neighbor_j(i, j, k)
            nk = neighbor_k(i, j, k)
            neighborhood_state(:, i, j, k) = &
              neighborhood_state(:, i, j, k) + &
              alpha_neighbor(i, j, k) * &
                geometry%volume_fraction(ni, nj, nk) * &
                provisional_state(:, ni, nj, nk) / &
                real(neighborhood_count(ni, nj, nk), dp)
          end if
          neighborhood_state(:, i, j, k) = &
            neighborhood_state(:, i, j, k) / &
              neighborhood_volume(i, j, k)
        end do
      end do
    end do

    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          candidate(:, i, j, k) = candidate(:, i, j, k) + &
            alpha_self(i, j, k) * &
              real(neighborhood_count(i, j, k), dp) * &
              neighborhood_state(:, i, j, k)
          if (has_neighbor(i, j, k)) then
            ni = neighbor_i(i, j, k)
            nj = neighbor_j(i, j, k)
            nk = neighbor_k(i, j, k)
            candidate(:, ni, nj, nk) = candidate(:, ni, nj, nk) + &
              alpha_neighbor(i, j, k) * &
                neighborhood_state(:, i, j, k)
          end if
        end do
      end do
    end do
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          candidate(:, i, j, k) = candidate(:, i, j, k) / &
            real(neighborhood_count(i, j, k), dp)
        end do
      end do
    end do
    if (any(.not. ieee_is_finite(candidate))) return

    redistributed_state = candidate
    ok = .true.
  end subroutine reactive_eb_weighted_state_redistribute_3d

  subroutine select_axis_plane_state_receiver( &
      geometry, i, j, k, neighbor_i, neighbor_j, neighbor_k, ok)
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: i, j, k
    integer, intent(out) :: neighbor_i, neighbor_j, neighbor_k
    logical, intent(out) :: ok

    real(dp), parameter :: axis_tolerance = &
      512.0_dp * epsilon(1.0_dp)
    real(dp) :: normal(3), normal_norm
    integer :: axis, offset

    neighbor_i = i
    neighbor_j = j
    neighbor_k = k
    ok = .false.
    normal = [ &
      geometry%x_face_fraction(i, j, k) - &
        geometry%x_face_fraction(i - 1, j, k), &
      geometry%y_face_fraction(i, j, k) - &
        geometry%y_face_fraction(i, j - 1, k), &
      geometry%z_face_fraction(i, j, k) - &
        geometry%z_face_fraction(i, j, k - 1)]
    normal_norm = sqrt(sum(normal**2))
    if (.not. ieee_is_finite(normal_norm) .or. &
        normal_norm <= tiny(1.0_dp)) return
    if (count(abs(normal) > axis_tolerance * normal_norm) /= 1) return
    axis = maxloc(abs(normal), dim=1)
    offset = merge(1, -1, normal(axis) > 0.0_dp)
    select case (axis)
    case (1)
      neighbor_i = i + offset
      if (neighbor_i < 1 .or. neighbor_i > geometry%nx) return
      if (offset > 0) then
        if (geometry%x_face_fraction(i, j, k) <= 0.0_dp) return
      else
        if (geometry%x_face_fraction(i - 1, j, k) <= 0.0_dp) return
      end if
    case (2)
      neighbor_j = j + offset
      if (neighbor_j < 1 .or. neighbor_j > geometry%ny) return
      if (offset > 0) then
        if (geometry%y_face_fraction(i, j, k) <= 0.0_dp) return
      else
        if (geometry%y_face_fraction(i, j - 1, k) <= 0.0_dp) return
      end if
    case (3)
      neighbor_k = k + offset
      if (neighbor_k < 1 .or. neighbor_k > geometry%nz) return
      if (offset > 0) then
        if (geometry%z_face_fraction(i, j, k) <= 0.0_dp) return
      else
        if (geometry%z_face_fraction(i, j, k - 1) <= 0.0_dp) return
      end if
    case default
      return
    end select
    if (geometry%cell_type(neighbor_i, neighbor_j, neighbor_k) /= &
        eb_regular_cell_3d) return
    ok = .true.
  end subroutine select_axis_plane_state_receiver

  subroutine advance_reactive_eb_state_redistributed_3d( &
      species, state, temperature, geometry, conservative_rhs, dt, &
      new_state, new_temperature, ok, target_volume_fraction)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: conservative_rhs(:, :, :, :), dt
    real(dp), intent(out) :: new_state(:, :, :, :)
    real(dp), intent(out) :: new_temperature(:, :, :)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: target_volume_fraction

    real(dp), allocatable :: provisional_state(:, :, :, :)
    real(dp), allocatable :: redistributed_state(:, :, :, :)
    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    real(dp), allocatable :: primitive(:)
    real(dp) :: recovered_temperature, sound_speed, target
    logical :: local_ok
    integer :: i, j, k, nvar

    new_state = 0.0_dp
    new_temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        size(state, 4) /= geometry%nz .or. &
        any(shape(temperature) /= &
          [geometry%nx, geometry%ny, geometry%nz]) .or. &
        any(shape(conservative_rhs) /= shape(state)) .or. &
        any(shape(new_state) /= shape(state)) .or. &
        any(shape(new_temperature) /= shape(temperature))) return
    new_state = state
    new_temperature = temperature
    if (.not. ieee_is_finite(dt) .or. dt < 0.0_dp .or. &
        any(.not. ieee_is_finite(state)) .or. &
        any(.not. ieee_is_finite(temperature)) .or. &
        any(.not. ieee_is_finite(conservative_rhs))) return

    allocate(provisional_state, source=state)
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          provisional_state(:, i, j, k) = state(:, i, j, k) + &
            dt * conservative_rhs(:, i, j, k)
        end do
      end do
    end do
    target = 0.5_dp
    if (present(target_volume_fraction)) target = target_volume_fraction
    allocate(redistributed_state(nvar, geometry%nx, geometry%ny, geometry%nz))
    call reactive_eb_weighted_state_redistribute_3d( &
      geometry, provisional_state, redistributed_state, local_ok, target)
    if (.not. local_ok) return

    allocate(candidate_state, source=state)
    allocate(candidate_temperature, source=temperature)
    allocate(primitive(reactive_nprim(size(species))))
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          if (temperature(i, j, k) <= 0.0_dp) return
          candidate_state(:, i, j, k) = &
            redistributed_state(:, i, j, k)
          call reactive_conserved_to_primitive( &
            species, candidate_state(:, i, j, k), temperature(i, j, k), &
            primitive, recovered_temperature, sound_speed, local_ok)
          if (.not. local_ok) return
          candidate_temperature(i, j, k) = recovered_temperature
        end do
      end do
    end do
    if (any(.not. ieee_is_finite(candidate_state)) .or. &
        any(.not. ieee_is_finite(candidate_temperature))) return

    new_state = candidate_state
    new_temperature = candidate_temperature
    ok = .true.
  end subroutine advance_reactive_eb_state_redistributed_3d

end module eb_reactive_redistribution_3d_mod
