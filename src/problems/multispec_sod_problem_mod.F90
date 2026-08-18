module multispec_sod_problem_mod
  use precision_mod, only: dp
  use state_indices_mod, only: ncons, nprim, qrho, qu, qv, qw, qp
  use state_conversion_mod, only: primitive_to_conserved
  use multispecies_state_mod, only: &
    max_supported_species, multispecies_nvar, multispecies_state_from_base
  use simulation_config_multispecies_mod, only: multispec_sod_config
  implicit none
  private

  public :: initialize_multispec_sod_problem

contains

  subroutine initialize_multispec_sod_problem( &
      x, nx, nspecies, gamma, problem, conserved, ok)
    integer, intent(in) :: nx, nspecies
    real(dp), intent(in) :: x(nx), gamma
    type(multispec_sod_config), intent(in) :: problem
    real(dp), intent(out) :: conserved(:, 0:)
    logical, intent(out) :: ok

    real(dp) :: primitive(nprim), base_state(ncons)
    real(dp) :: mass_fractions(max_supported_species)
    logical :: cell_ok
    integer :: i, nvar

    conserved = 0.0_dp
    ok = .false.
    nvar = multispecies_nvar(nspecies)
    if (nvar == 0 .or. size(conserved, 1) /= nvar) return
    if (ubound(conserved, 2) < nx + 1) return

    do i = 1, nx
      mass_fractions = 0.0_dp
      if (x(i) < problem%discontinuity) then
        primitive(qrho) = problem%rho_left
        primitive(qu) = problem%velocity_left
        primitive(qp) = problem%pressure_left
        mass_fractions(1:nspecies) = &
          problem%mass_fractions_left(1:nspecies)
      else
        primitive(qrho) = problem%rho_right
        primitive(qu) = problem%velocity_right
        primitive(qp) = problem%pressure_right
        mass_fractions(1:nspecies) = &
          problem%mass_fractions_right(1:nspecies)
      end if
      primitive(qv) = 0.0_dp
      primitive(qw) = 0.0_dp

      call primitive_to_conserved(primitive, gamma, base_state, cell_ok)
      if (.not. cell_ok) return
      call multispecies_state_from_base( &
        base_state, mass_fractions, nspecies, gamma, conserved(:, i), cell_ok)
      if (.not. cell_ok) return
    end do

    conserved(:, 0) = conserved(:, 1)
    conserved(:, nx + 1) = conserved(:, nx)
    ok = .true.
  end subroutine initialize_multispec_sod_problem

end module multispec_sod_problem_mod
