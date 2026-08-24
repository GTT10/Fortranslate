module eb_reactive_hydro_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use eb_geometry_2d_mod, only: eb_geometry_2d
  use eb_reactive_reconstruction_2d_mod, only: &
    build_reactive_eb_face_center_fluxes_2d, &
    interpolate_reactive_eb_face_centroid_fluxes_2d
  use eb_reactive_wall_flux_2d_mod, only: reactive_eb_flux_divergence_2d
  use eb_reactive_redistribution_2d_mod, only: &
    advance_reactive_eb_state_redistributed_2d
  implicit none
  private

  public :: reactive_eb_outflow_riemann_fluxes_2d
  public :: advance_reactive_eb_hydro_2d

contains

  subroutine reactive_eb_outflow_riemann_fluxes_2d( &
      species, state, temperature, geometry, solver, x_flux, y_flux, ok, &
      reconstruction, limiter, dt)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    type(eb_geometry_2d), intent(in) :: geometry
    character(len=*), intent(in) :: solver
    real(dp), intent(out) :: x_flux(:, 0:, :), y_flux(:, :, 0:)
    logical, intent(out) :: ok
    character(len=*), intent(in), optional :: reconstruction, limiter
    real(dp), intent(in), optional :: dt

    real(dp), allocatable :: center_x(:, :, :), center_y(:, :, :)
    character(len=32) :: selected_reconstruction, selected_limiter
    real(dp) :: selected_dt
    logical :: local_ok
    integer :: nvar

    x_flux = 0.0_dp
    y_flux = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    selected_reconstruction = "pcm"
    selected_limiter = "mc"
    selected_dt = 0.0_dp
    if (present(reconstruction)) selected_reconstruction = trim(reconstruction)
    if (present(limiter)) selected_limiter = trim(limiter)
    if (present(dt)) selected_dt = dt
    allocate(center_x(nvar, 0:geometry%nx, geometry%ny))
    allocate(center_y(nvar, geometry%nx, 0:geometry%ny))
    call build_reactive_eb_face_center_fluxes_2d( &
      species, state, temperature, geometry, solver, &
      selected_reconstruction, selected_limiter, selected_dt, center_x, &
      center_y, local_ok)
    if (.not. local_ok) return
    call interpolate_reactive_eb_face_centroid_fluxes_2d( &
      geometry, center_x, center_y, x_flux, y_flux, local_ok)
    if (.not. local_ok) return
    ok = .true.
  end subroutine reactive_eb_outflow_riemann_fluxes_2d

  subroutine advance_reactive_eb_hydro_2d( &
      species, state, temperature, geometry, solver, dt, &
      new_state, new_temperature, ok, target_volume_fraction, &
      reconstruction, limiter, state_redist_max_order)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    type(eb_geometry_2d), intent(in) :: geometry
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    real(dp), intent(out) :: new_state(:, :, :), new_temperature(:, :)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: target_volume_fraction
    character(len=*), intent(in), optional :: reconstruction, limiter
    integer, intent(in), optional :: state_redist_max_order

    real(dp), allocatable :: x_flux(:, :, :), y_flux(:, :, :)
    real(dp), allocatable :: conservative_rhs(:, :, :)
    logical :: local_ok
    integer :: nvar, selected_max_order
    real(dp) :: selected_target

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
      species, state, temperature, geometry, solver, x_flux, y_flux, &
      local_ok, reconstruction, limiter, dt)
    if (.not. local_ok) return

    allocate(conservative_rhs(nvar, geometry%nx, geometry%ny))
    call reactive_eb_flux_divergence_2d( &
      species, state, temperature, geometry, x_flux, y_flux, &
      conservative_rhs, local_ok)
    if (.not. local_ok) return

    selected_target = 0.5_dp
    if (present(target_volume_fraction)) selected_target = target_volume_fraction
    selected_max_order = 0
    if (present(state_redist_max_order)) &
      selected_max_order = state_redist_max_order
    call advance_reactive_eb_state_redistributed_2d( &
      species, state, temperature, geometry, conservative_rhs, dt, &
      new_state, new_temperature, local_ok, selected_target, &
      selected_max_order)
    if (.not. local_ok) return
    ok = .true.
  end subroutine advance_reactive_eb_hydro_2d

end module eb_reactive_hydro_2d_mod
