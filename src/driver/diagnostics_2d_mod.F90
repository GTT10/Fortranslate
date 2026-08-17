module diagnostics_2d_mod
  use precision_mod, only: dp
  use state_indices_mod, only: ncons, nprim, qrho, qp
  use state_conversion_mod, only: conserved_to_primitive
  implicit none
  private

  public :: integrated_conserved_quantities_2d
  public :: primitive_extrema_2d

contains

  pure subroutine integrated_conserved_quantities_2d( &
      conserved, nx, ny, dx, dy, totals)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: conserved(ncons, nx, ny)
    real(dp), intent(in) :: dx, dy
    real(dp), intent(out) :: totals(ncons)

    totals = sum(sum(conserved, dim=3), dim=2) * dx * dy
  end subroutine integrated_conserved_quantities_2d

  subroutine primitive_extrema_2d( &
      conserved, nx, ny, gamma, minimum_density, maximum_density, &
      minimum_pressure, maximum_pressure, ok)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: conserved(ncons, nx, ny)
    real(dp), intent(in) :: gamma
    real(dp), intent(out) :: minimum_density, maximum_density
    real(dp), intent(out) :: minimum_pressure, maximum_pressure
    logical, intent(out) :: ok

    real(dp) :: primitive(nprim)
    logical :: cell_ok
    integer :: i, j

    minimum_density = huge(1.0_dp)
    maximum_density = -huge(1.0_dp)
    minimum_pressure = huge(1.0_dp)
    maximum_pressure = -huge(1.0_dp)
    ok = .true.

    do j = 1, ny
      do i = 1, nx
        call conserved_to_primitive( &
          conserved(:, i, j), gamma, primitive, cell_ok)
        if (.not. cell_ok) then
          ok = .false.
          return
        end if
        minimum_density = min(minimum_density, primitive(qrho))
        maximum_density = max(maximum_density, primitive(qrho))
        minimum_pressure = min(minimum_pressure, primitive(qp))
        maximum_pressure = max(maximum_pressure, primitive(qp))
      end do
    end do
  end subroutine primitive_extrema_2d

end module diagnostics_2d_mod
