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

end module eb_reactive_redistribution_2d_mod
