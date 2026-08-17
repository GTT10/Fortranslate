module sedov_problem_mod
  use precision_mod, only: dp
  use state_indices_mod, only: ncons, nprim, qrho, qu, qv, qw, qp
  use state_conversion_mod, only: primitive_to_conserved
  use simulation_config_mod, only: sedov_config
  implicit none
  private

  public :: initialize_sedov_problem

contains

  subroutine initialize_sedov_problem( &
      x, nx, gamma, sedov, conserved, ok)
    integer, intent(in) :: nx
    real(dp), intent(in) :: x(nx), gamma
    type(sedov_config), intent(in) :: sedov
    real(dp), intent(out) :: conserved(ncons, 0:nx + 1)
    logical, intent(out) :: ok

    real(dp) :: primitive(nprim)
    logical :: cell_ok
    integer :: i

    conserved = 0.0_dp
    ok = .true.

    do i = 1, nx
      primitive(qrho) = sedov%ambient_density
      primitive(qu) = sedov%initial_velocity
      primitive(qv) = 0.0_dp
      primitive(qw) = 0.0_dp
      if (abs(x(i) - sedov%blast_center) <= sedov%blast_radius) then
        primitive(qp) = sedov%blast_pressure
      else
        primitive(qp) = sedov%ambient_pressure
      end if

      call primitive_to_conserved(primitive, gamma, conserved(:, i), cell_ok)
      if (.not. cell_ok) then
        ok = .false.
        return
      end if
    end do

    conserved(:, 0) = conserved(:, 1)
    conserved(:, nx + 1) = conserved(:, nx)
  end subroutine initialize_sedov_problem

end module sedov_problem_mod
