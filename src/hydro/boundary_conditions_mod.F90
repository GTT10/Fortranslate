module boundary_conditions_mod
  use precision_mod, only: dp
  use state_indices_mod, only: ncons
  implicit none
  private

  public :: apply_boundary_conditions
  public :: apply_outflow_boundaries
  public :: apply_periodic_boundaries

contains

  pure subroutine apply_boundary_conditions(conserved, nx, boundary_condition, ok)
    integer, intent(in) :: nx
    real(dp), intent(inout) :: conserved(ncons, 0:nx + 1)
    character(len=*), intent(in) :: boundary_condition
    logical, intent(out) :: ok

    select case (trim(boundary_condition))
    case ("outflow")
      call apply_outflow_boundaries(conserved, nx)
      ok = .true.
    case ("periodic")
      call apply_periodic_boundaries(conserved, nx)
      ok = .true.
    case default
      ok = .false.
    end select
  end subroutine apply_boundary_conditions

  pure subroutine apply_outflow_boundaries(conserved, nx)
    integer, intent(in) :: nx
    real(dp), intent(inout) :: conserved(ncons, 0:nx + 1)

    conserved(:, 0) = conserved(:, 1)
    conserved(:, nx + 1) = conserved(:, nx)
  end subroutine apply_outflow_boundaries

  pure subroutine apply_periodic_boundaries(conserved, nx)
    integer, intent(in) :: nx
    real(dp), intent(inout) :: conserved(ncons, 0:nx + 1)

    conserved(:, 0) = conserved(:, nx)
    conserved(:, nx + 1) = conserved(:, 1)
  end subroutine apply_periodic_boundaries

end module boundary_conditions_mod
