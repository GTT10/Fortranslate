module diagnostics_mod
  use precision_mod, only: dp
  use state_indices_mod, only: ncons, nprim, qrho, qp
  use state_conversion_mod, only: conserved_to_primitive
  implicit none
  private

  public :: integrated_conserved_quantities
  public :: primitive_extrema

contains

  pure subroutine integrated_conserved_quantities(conserved, nx, dx, totals)
    integer, intent(in) :: nx
    real(dp), intent(in) :: conserved(ncons, 0:nx + 1)
    real(dp), intent(in) :: dx
    real(dp), intent(out) :: totals(ncons)

    totals = sum(conserved(:, 1:nx), dim=2) * dx
  end subroutine integrated_conserved_quantities

  subroutine primitive_extrema( &
      conserved, nx, gamma, minimum_density, minimum_pressure, ok)
    integer, intent(in) :: nx
    real(dp), intent(in) :: conserved(ncons, 0:nx + 1)
    real(dp), intent(in) :: gamma
    real(dp), intent(out) :: minimum_density, minimum_pressure
    logical, intent(out) :: ok
    real(dp) :: primitive(nprim)
    logical :: cell_ok
    integer :: i

    minimum_density = huge(1.0_dp)
    minimum_pressure = huge(1.0_dp)
    ok = .true.

    do i = 1, nx
      call conserved_to_primitive(conserved(:, i), gamma, primitive, cell_ok)
      if (.not. cell_ok) then
        ok = .false.
        return
      end if
      minimum_density = min(minimum_density, primitive(qrho))
      minimum_pressure = min(minimum_pressure, primitive(qp))
    end do
  end subroutine primitive_extrema

end module diagnostics_mod
