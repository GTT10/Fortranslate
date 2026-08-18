module multispecies_flux_mod
  use precision_mod, only: dp
  use state_indices_mod, only: irho, ncons, nbase
  use multispecies_state_mod, only: &
    max_supported_species, multispecies_nvar, species_component, &
    mass_fractions_from_state
  use riemann_flux_mod, only: compute_riemann_flux_x
  use directional_flux_mod, only: &
    rotate_conserved_y_to_x, rotate_flux_x_to_y
  implicit none
  private

  public :: compute_multispecies_flux_x
  public :: compute_multispecies_flux_y

contains

  pure subroutine compute_multispecies_flux_x( &
      left_state, right_state, nspecies, gamma, solver, flux, ok)
    integer, intent(in) :: nspecies
    real(dp), intent(in) :: left_state(:), right_state(:), gamma
    character(len=*), intent(in) :: solver
    real(dp), intent(out) :: flux(:)
    logical, intent(out) :: ok

    real(dp) :: base_flux(ncons)
    real(dp) :: left_mass_fractions(max_supported_species)
    real(dp) :: right_mass_fractions(max_supported_species)
    real(dp) :: donor_mass_fractions(max_supported_species)
    real(dp) :: mass_flux, zero_threshold
    logical :: left_ok, right_ok, base_ok
    integer :: species, nvar

    flux = 0.0_dp
    ok = .false.
    nvar = multispecies_nvar(nspecies)
    if (nvar == 0) return
    if (size(left_state) /= nvar .or. size(right_state) /= nvar) return
    if (size(flux) /= nvar) return

    call mass_fractions_from_state( &
      left_state, nspecies, left_mass_fractions, left_ok)
    call mass_fractions_from_state( &
      right_state, nspecies, right_mass_fractions, right_ok)
    if (.not. (left_ok .and. right_ok)) return

    call compute_riemann_flux_x( &
      left_state(1:ncons), right_state(1:ncons), gamma, solver, &
      base_flux, base_ok)
    if (.not. base_ok) return

    flux(1:ncons) = base_flux
    mass_flux = base_flux(irho)
    zero_threshold = sqrt(epsilon(1.0_dp)) * max(1.0_dp, abs(mass_flux))
    if (mass_flux > zero_threshold) then
      donor_mass_fractions(1:nspecies) = left_mass_fractions(1:nspecies)
    else if (mass_flux < -zero_threshold) then
      donor_mass_fractions(1:nspecies) = right_mass_fractions(1:nspecies)
    else
      donor_mass_fractions(1:nspecies) = 0.5_dp * &
        (left_mass_fractions(1:nspecies) + &
         right_mass_fractions(1:nspecies))
    end if

    if (nspecies == 1) then
      flux(species_component(1)) = mass_flux
    else
      do species = 1, nspecies - 1
        flux(species_component(species)) = &
          mass_flux * donor_mass_fractions(species)
      end do
      flux(species_component(nspecies)) = mass_flux - &
        sum(flux(species_component(1):species_component(nspecies - 1)))
    end if
    ok = .true.
  end subroutine compute_multispecies_flux_x

  pure subroutine compute_multispecies_flux_y( &
      lower_state, upper_state, nspecies, gamma, solver, flux, ok)
    integer, intent(in) :: nspecies
    real(dp), intent(in) :: lower_state(:), upper_state(:), gamma
    character(len=*), intent(in) :: solver
    real(dp), intent(out) :: flux(:)
    logical, intent(out) :: ok

    real(dp) :: lower_rotated(ncons), upper_rotated(ncons)
    real(dp) :: lower_extended(nbase + max_supported_species)
    real(dp) :: upper_extended(nbase + max_supported_species)
    real(dp) :: rotated_flux(nbase + max_supported_species)
    real(dp) :: base_flux_y(ncons)
    integer :: nvar

    flux = 0.0_dp
    ok = .false.
    nvar = multispecies_nvar(nspecies)
    if (nvar == 0) return
    if (size(lower_state) /= nvar .or. size(upper_state) /= nvar) return
    if (size(flux) /= nvar) return

    call rotate_conserved_y_to_x(lower_state(1:ncons), lower_rotated)
    call rotate_conserved_y_to_x(upper_state(1:ncons), upper_rotated)
    lower_extended = 0.0_dp
    upper_extended = 0.0_dp
    lower_extended(1:ncons) = lower_rotated
    upper_extended(1:ncons) = upper_rotated
    lower_extended(ncons + 1:nvar) = lower_state(ncons + 1:nvar)
    upper_extended(ncons + 1:nvar) = upper_state(ncons + 1:nvar)

    call compute_multispecies_flux_x( &
      lower_extended(1:nvar), upper_extended(1:nvar), nspecies, gamma, &
      solver, rotated_flux(1:nvar), ok)
    if (.not. ok) return

    call rotate_flux_x_to_y(rotated_flux(1:ncons), base_flux_y)
    flux(1:ncons) = base_flux_y
    flux(ncons + 1:nvar) = rotated_flux(ncons + 1:nvar)
  end subroutine compute_multispecies_flux_y

end module multispecies_flux_mod
