module boundary_conditions_mod
  use precision_mod, only: dp
  use state_indices_mod, only: ncons
  implicit none
  private

  public :: apply_outflow_boundaries

contains

  pure subroutine apply_outflow_boundaries(conserved, nx)
    integer, intent(in) :: nx
    real(dp), intent(inout) :: conserved(ncons, 0:nx + 1)

    conserved(:, 0) = conserved(:, 1)
    conserved(:, nx + 1) = conserved(:, nx)
  end subroutine apply_outflow_boundaries

end module boundary_conditions_mod
