module shu_osher_problem_mod
  use precision_mod, only: dp
  use state_indices_mod, only: ncons, nprim, qrho, qu, qv, qw, qp
  use state_conversion_mod, only: primitive_to_conserved
  use simulation_config_mod, only: shu_osher_config
  implicit none
  private

  public :: initialize_shu_osher_problem

contains

  subroutine initialize_shu_osher_problem( &
      x, nx, gamma, shu_osher, conserved, ok)
    integer, intent(in) :: nx
    real(dp), intent(in) :: x(nx), gamma
    type(shu_osher_config), intent(in) :: shu_osher
    real(dp), intent(out) :: conserved(ncons, 0:nx + 1)
    logical, intent(out) :: ok

    real(dp) :: primitive(nprim)
    logical :: cell_ok
    integer :: i

    conserved = 0.0_dp
    ok = .true.

    do i = 1, nx
      if (x(i) < shu_osher%shock_location) then
        primitive(qrho) = shu_osher%left_density
        primitive(qu) = shu_osher%left_velocity
        primitive(qp) = shu_osher%left_pressure
      else
        primitive(qrho) = shu_osher%density_base + &
          shu_osher%density_amplitude * &
          sin(shu_osher%density_wavenumber * x(i))
        primitive(qu) = shu_osher%right_velocity
        primitive(qp) = shu_osher%right_pressure
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
  end subroutine initialize_shu_osher_problem

end module shu_osher_problem_mod
