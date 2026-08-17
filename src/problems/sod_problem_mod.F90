module sod_problem_mod
  use precision_mod, only: dp
  use state_indices_mod, only: ncons, nprim, qrho, qu, qv, qw, qp
  use state_conversion_mod, only: primitive_to_conserved
  use simulation_config_mod, only: sod_config
  implicit none
  private

  public :: initialize_sod_problem

contains

  subroutine initialize_sod_problem(x, nx, gamma, sod, conserved, ok)
    integer, intent(in) :: nx
    real(dp), intent(in) :: x(nx), gamma
    type(sod_config), intent(in) :: sod
    real(dp), intent(out) :: conserved(ncons, 0:nx + 1)
    logical, intent(out) :: ok
    real(dp) :: primitive(nprim)
    logical :: cell_ok
    integer :: i

    conserved = 0.0_dp
    ok = .true.

    do i = 1, nx
      if (x(i) < sod%discontinuity) then
        primitive(qrho) = sod%rho_left
        primitive(qu) = sod%velocity_left
        primitive(qp) = sod%pressure_left
      else
        primitive(qrho) = sod%rho_right
        primitive(qu) = sod%velocity_right
        primitive(qp) = sod%pressure_right
      end if
      primitive(qv) = 0.0_dp
      primitive(qw) = 0.0_dp

      call primitive_to_conserved(primitive, gamma, conserved(:, i), cell_ok)
      if (.not. cell_ok) then
        ok = .false.
        return
      end if
    end do

    conserved(:, 0) = conserved(:, 1)
    conserved(:, nx + 1) = conserved(:, nx)
  end subroutine initialize_sod_problem

end module sod_problem_mod
