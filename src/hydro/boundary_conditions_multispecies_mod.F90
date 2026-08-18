module boundary_conditions_multispecies_mod
  use precision_mod, only: dp
  use state_indices_mod, only: ncons
  implicit none
  private

  public :: apply_multispecies_boundary_conditions

contains

  pure subroutine apply_multispecies_boundary_conditions( &
      conserved, nx, boundary_condition, ok)
    integer, intent(in) :: nx
    real(dp), intent(inout) :: conserved(:, 0:)
    character(len=*), intent(in) :: boundary_condition
    logical, intent(out) :: ok

    ok = .false.
    if (size(conserved, 1) <= ncons) return
    if (ubound(conserved, 2) < nx + 1) return

    select case (trim(boundary_condition))
    case ("outflow")
      conserved(:, 0) = conserved(:, 1)
      conserved(:, nx + 1) = conserved(:, nx)
      ok = .true.
    case ("periodic")
      conserved(:, 0) = conserved(:, nx)
      conserved(:, nx + 1) = conserved(:, 1)
      ok = .true.
    case default
      continue
    end select
  end subroutine apply_multispecies_boundary_conditions

end module boundary_conditions_multispecies_mod
