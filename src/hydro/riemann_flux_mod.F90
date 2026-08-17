module riemann_flux_mod
  use precision_mod, only: dp
  use state_indices_mod, only: ncons
  use riemann_rusanov_mod, only: rusanov_flux_x
  use riemann_pelec_mod, only: pelec_riemann_flux_x
  implicit none
  private

  public :: compute_riemann_flux_x

contains

  pure subroutine compute_riemann_flux_x( &
      left_state, right_state, gamma, solver, flux, ok)
    real(dp), intent(in) :: left_state(ncons), right_state(ncons)
    real(dp), intent(in) :: gamma
    character(len=*), intent(in) :: solver
    real(dp), intent(out) :: flux(ncons)
    logical, intent(out) :: ok

    select case (trim(solver))
    case ("rusanov")
      call rusanov_flux_x(left_state, right_state, gamma, flux, ok)
    case ("pelec")
      call pelec_riemann_flux_x(left_state, right_state, gamma, flux, ok)
    case default
      flux = 0.0_dp
      ok = .false.
    end select
  end subroutine compute_riemann_flux_x

end module riemann_flux_mod
