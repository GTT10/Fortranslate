module directional_flux_mod
  use precision_mod, only: dp
  use state_indices_mod, only: &
    irho, imx, imy, imz, iet, ncons, qrho, qu, qv, qw, qp, nprim
  use riemann_rusanov_mod, only: euler_physical_flux_x
  use riemann_flux_mod, only: compute_riemann_flux_x
  implicit none
  private

  public :: rotate_conserved_y_to_x
  public :: rotate_conserved_x_to_y
  public :: rotate_primitive_y_to_x
  public :: rotate_primitive_x_to_y
  public :: rotate_flux_x_to_y
  public :: euler_physical_flux_y
  public :: compute_riemann_flux_y

contains

  pure subroutine rotate_conserved_y_to_x(state_y, state_x)
    real(dp), intent(in) :: state_y(ncons)
    real(dp), intent(out) :: state_x(ncons)

    state_x(irho) = state_y(irho)
    state_x(imx) = state_y(imy)
    state_x(imy) = state_y(imx)
    state_x(imz) = state_y(imz)
    state_x(iet) = state_y(iet)
  end subroutine rotate_conserved_y_to_x

  pure subroutine rotate_conserved_x_to_y(state_x, state_y)
    real(dp), intent(in) :: state_x(ncons)
    real(dp), intent(out) :: state_y(ncons)

    call rotate_conserved_y_to_x(state_x, state_y)
  end subroutine rotate_conserved_x_to_y

  pure subroutine rotate_primitive_y_to_x(primitive_y, primitive_x)
    real(dp), intent(in) :: primitive_y(nprim)
    real(dp), intent(out) :: primitive_x(nprim)

    primitive_x(qrho) = primitive_y(qrho)
    primitive_x(qu) = primitive_y(qv)
    primitive_x(qv) = primitive_y(qu)
    primitive_x(qw) = primitive_y(qw)
    primitive_x(qp) = primitive_y(qp)
  end subroutine rotate_primitive_y_to_x

  pure subroutine rotate_primitive_x_to_y(primitive_x, primitive_y)
    real(dp), intent(in) :: primitive_x(nprim)
    real(dp), intent(out) :: primitive_y(nprim)

    call rotate_primitive_y_to_x(primitive_x, primitive_y)
  end subroutine rotate_primitive_x_to_y

  pure subroutine rotate_flux_x_to_y(flux_x, flux_y)
    real(dp), intent(in) :: flux_x(ncons)
    real(dp), intent(out) :: flux_y(ncons)

    flux_y(irho) = flux_x(irho)
    flux_y(imx) = flux_x(imy)
    flux_y(imy) = flux_x(imx)
    flux_y(imz) = flux_x(imz)
    flux_y(iet) = flux_x(iet)
  end subroutine rotate_flux_x_to_y

  pure subroutine euler_physical_flux_y(conserved, gamma, flux, ok)
    real(dp), intent(in) :: conserved(ncons)
    real(dp), intent(in) :: gamma
    real(dp), intent(out) :: flux(ncons)
    logical, intent(out) :: ok

    real(dp) :: rotated_state(ncons), rotated_flux(ncons)

    call rotate_conserved_y_to_x(conserved, rotated_state)
    call euler_physical_flux_x(rotated_state, gamma, rotated_flux, ok)
    if (.not. ok) then
      flux = 0.0_dp
      return
    end if
    call rotate_flux_x_to_y(rotated_flux, flux)
  end subroutine euler_physical_flux_y

  pure subroutine compute_riemann_flux_y( &
      lower_state, upper_state, gamma, solver, flux, ok)
    real(dp), intent(in) :: lower_state(ncons), upper_state(ncons)
    real(dp), intent(in) :: gamma
    character(len=*), intent(in) :: solver
    real(dp), intent(out) :: flux(ncons)
    logical, intent(out) :: ok

    real(dp) :: lower_rotated(ncons), upper_rotated(ncons)
    real(dp) :: rotated_flux(ncons)

    call rotate_conserved_y_to_x(lower_state, lower_rotated)
    call rotate_conserved_y_to_x(upper_state, upper_rotated)
    call compute_riemann_flux_x( &
      lower_rotated, upper_rotated, gamma, solver, rotated_flux, ok)
    if (.not. ok) then
      flux = 0.0_dp
      return
    end if
    call rotate_flux_x_to_y(rotated_flux, flux)
  end subroutine compute_riemann_flux_y

end module directional_flux_mod
