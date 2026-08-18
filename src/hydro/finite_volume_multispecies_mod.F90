module finite_volume_multispecies_mod
  use precision_mod, only: dp
  use multispecies_state_mod, only: multispecies_nvar
  use reconstruction_multispecies_mod, only: reconstruct_multispecies_faces
  use multispecies_flux_mod, only: compute_multispecies_flux_x
  implicit none
  private

  public :: compute_multispecies_rhs

contains

  subroutine compute_multispecies_rhs( &
      conserved, nx, nspecies, dx, gamma, rhs, ok, reconstruction, limiter, &
      boundary_condition, riemann_solver, dt, plm_order, use_flattening)
    integer, intent(in) :: nx, nspecies
    real(dp), intent(in) :: conserved(:, 0:), dx, gamma
    real(dp), intent(out) :: rhs(:, 1:)
    logical, intent(out) :: ok
    character(len=*), intent(in), optional :: reconstruction, limiter
    character(len=*), intent(in), optional :: boundary_condition
    character(len=*), intent(in), optional :: riemann_solver
    real(dp), intent(in), optional :: dt
    integer, intent(in), optional :: plm_order
    logical, intent(in), optional :: use_flattening

    real(dp), allocatable :: left_faces(:, :), right_faces(:, :)
    real(dp), allocatable :: face_flux(:, :)
    character(len=32) :: reconstruction_name, limiter_name, boundary_name
    character(len=32) :: riemann_name
    real(dp) :: local_dtdx
    integer :: slope_order, nvar, i
    logical :: flattening_enabled, reconstruction_ok, flux_ok

    ok = .false.
    nvar = multispecies_nvar(nspecies)
    if (nvar == 0 .or. nx < 4 .or. dx <= 0.0_dp) return
    if (size(conserved, 1) /= nvar .or. ubound(conserved, 2) < nx + 1) return
    if (size(rhs, 1) /= nvar .or. ubound(rhs, 2) < nx) return

    reconstruction_name = "pcm"
    limiter_name = "mc"
    boundary_name = "outflow"
    riemann_name = "rusanov"
    slope_order = 2
    flattening_enabled = .false.
    local_dtdx = 0.0_dp
    if (present(reconstruction)) reconstruction_name = trim(reconstruction)
    if (present(limiter)) limiter_name = trim(limiter)
    if (present(boundary_condition)) boundary_name = trim(boundary_condition)
    if (present(riemann_solver)) riemann_name = trim(riemann_solver)
    if (present(plm_order)) slope_order = plm_order
    if (present(use_flattening)) flattening_enabled = use_flattening
    if (present(dt)) local_dtdx = dt / dx

    allocate(left_faces(nvar, 0:nx), right_faces(nvar, 0:nx))
    allocate(face_flux(nvar, 0:nx))
    rhs = 0.0_dp

    call reconstruct_multispecies_faces( &
      conserved, nx, nspecies, gamma, reconstruction_name, limiter_name, &
      boundary_name, left_faces, right_faces, reconstruction_ok, &
      dtdx=local_dtdx, plm_order=slope_order, &
      use_flattening=flattening_enabled)
    if (.not. reconstruction_ok) return

    do i = 0, nx
      call compute_multispecies_flux_x( &
        left_faces(:, i), right_faces(:, i), nspecies, gamma, riemann_name, &
        face_flux(:, i), flux_ok)
      if (.not. flux_ok) return
    end do

    do i = 1, nx
      rhs(:, i) = -(face_flux(:, i) - face_flux(:, i - 1)) / dx
    end do
    ok = .true.
  end subroutine compute_multispecies_rhs

end module finite_volume_multispecies_mod
