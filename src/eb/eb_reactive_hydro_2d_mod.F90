module eb_reactive_hydro_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar, reactive_riemann_flux_x
  use reactive_2d_mod, only: reactive_riemann_flux_y
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_covered_cell
  use eb_reactive_wall_flux_2d_mod, only: reactive_eb_flux_divergence_2d
  use eb_reactive_redistribution_2d_mod, only: &
    advance_reactive_eb_state_redistributed_2d
  implicit none
  private

  public :: reactive_eb_outflow_riemann_fluxes_2d
  public :: advance_reactive_eb_hydro_2d

contains

  subroutine reactive_eb_outflow_riemann_fluxes_2d( &
      species, state, temperature, geometry, solver, x_flux, y_flux, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    type(eb_geometry_2d), intent(in) :: geometry
    character(len=*), intent(in) :: solver
    real(dp), intent(out) :: x_flux(:, 0:, :), y_flux(:, :, 0:)
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate_x(:, :, :), candidate_y(:, :, :)
    real(dp), allocatable :: face_flux(:)
    logical :: local_ok
    integer :: i, j, left_i, right_i, lower_j, upper_j, nvar

    x_flux = 0.0_dp
    y_flux = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        any(shape(temperature) /= [geometry%nx, geometry%ny]) .or. &
        size(x_flux, 1) /= nvar .or. &
        size(x_flux, 2) /= geometry%nx + 1 .or. &
        size(x_flux, 3) /= geometry%ny .or. &
        size(y_flux, 1) /= nvar .or. &
        size(y_flux, 2) /= geometry%nx .or. &
        size(y_flux, 3) /= geometry%ny + 1 .or. &
        len_trim(solver) == 0 .or. &
        any(.not. ieee_is_finite(state)) .or. &
        any(.not. ieee_is_finite(temperature))) return

    allocate(candidate_x(nvar, 0:geometry%nx, geometry%ny))
    allocate(candidate_y(nvar, geometry%nx, 0:geometry%ny))
    allocate(face_flux(nvar))
    candidate_x = 0.0_dp
    candidate_y = 0.0_dp

    do j = 1, geometry%ny
      do i = 0, geometry%nx
        if (geometry%x_face_fraction(i, j) <= 0.0_dp) cycle
        left_i = max(1, i)
        right_i = min(geometry%nx, i + 1)
        if (geometry%cell_type(left_i, j) == eb_covered_cell .or. &
            geometry%cell_type(right_i, j) == eb_covered_cell .or. &
            temperature(left_i, j) <= 0.0_dp .or. &
            temperature(right_i, j) <= 0.0_dp) return
        call reactive_riemann_flux_x( &
          species, state(:, left_i, j), state(:, right_i, j), &
          temperature(left_i, j), temperature(right_i, j), solver, &
          face_flux, local_ok)
        if (.not. local_ok) return
        candidate_x(:, i, j) = face_flux
      end do
    end do

    do j = 0, geometry%ny
      do i = 1, geometry%nx
        if (geometry%y_face_fraction(i, j) <= 0.0_dp) cycle
        lower_j = max(1, j)
        upper_j = min(geometry%ny, j + 1)
        if (geometry%cell_type(i, lower_j) == eb_covered_cell .or. &
            geometry%cell_type(i, upper_j) == eb_covered_cell .or. &
            temperature(i, lower_j) <= 0.0_dp .or. &
            temperature(i, upper_j) <= 0.0_dp) return
        call reactive_riemann_flux_y( &
          species, state(:, i, lower_j), state(:, i, upper_j), &
          temperature(i, lower_j), temperature(i, upper_j), solver, &
          face_flux, local_ok)
        if (.not. local_ok) return
        candidate_y(:, i, j) = face_flux
      end do
    end do
    if (any(.not. ieee_is_finite(candidate_x)) .or. &
        any(.not. ieee_is_finite(candidate_y))) return

    x_flux = candidate_x
    y_flux = candidate_y
    ok = .true.
  end subroutine reactive_eb_outflow_riemann_fluxes_2d

  subroutine advance_reactive_eb_hydro_2d( &
      species, state, temperature, geometry, solver, dt, &
      new_state, new_temperature, ok, target_volume_fraction)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    type(eb_geometry_2d), intent(in) :: geometry
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    real(dp), intent(out) :: new_state(:, :, :), new_temperature(:, :)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: target_volume_fraction

    real(dp), allocatable :: x_flux(:, :, :), y_flux(:, :, :)
    real(dp), allocatable :: conservative_rhs(:, :, :)
    logical :: local_ok
    integer :: nvar

    new_state = 0.0_dp
    new_temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        any(shape(temperature) /= [geometry%nx, geometry%ny]) .or. &
        any(shape(new_state) /= shape(state)) .or. &
        any(shape(new_temperature) /= shape(temperature))) return
    new_state = state
    new_temperature = temperature
    if (.not. ieee_is_finite(dt) .or. dt < 0.0_dp) return

    allocate(x_flux(nvar, 0:geometry%nx, geometry%ny))
    allocate(y_flux(nvar, geometry%nx, 0:geometry%ny))
    call reactive_eb_outflow_riemann_fluxes_2d( &
      species, state, temperature, geometry, solver, x_flux, y_flux, local_ok)
    if (.not. local_ok) return

    allocate(conservative_rhs(nvar, geometry%nx, geometry%ny))
    call reactive_eb_flux_divergence_2d( &
      species, state, temperature, geometry, x_flux, y_flux, &
      conservative_rhs, local_ok)
    if (.not. local_ok) return

    if (present(target_volume_fraction)) then
      call advance_reactive_eb_state_redistributed_2d( &
        species, state, temperature, geometry, conservative_rhs, dt, &
        new_state, new_temperature, local_ok, target_volume_fraction)
    else
      call advance_reactive_eb_state_redistributed_2d( &
        species, state, temperature, geometry, conservative_rhs, dt, &
        new_state, new_temperature, local_ok)
    end if
    if (.not. local_ok) return
    ok = .true.
  end subroutine advance_reactive_eb_hydro_2d

end module eb_reactive_hydro_2d_mod
