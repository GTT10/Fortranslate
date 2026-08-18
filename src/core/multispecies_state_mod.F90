module multispecies_state_mod
  use precision_mod, only: dp
  use constants_mod, only: density_floor
  use state_indices_mod, only: &
    irho, imx, imy, imz, iet, iei, item, ncons, nbase, nprim, qrho, qp
  use state_conversion_mod, only: state_is_physical, conserved_to_primitive
  implicit none
  private

  integer, parameter, public :: max_supported_species = 32
  real(dp), parameter, public :: species_negative_tolerance = 1.0e-12_dp
  real(dp), parameter, public :: species_closure_tolerance = 5.0e-11_dp
  real(dp), parameter, public :: thermodynamic_layout_tolerance = 5.0e-11_dp

  public :: multispecies_nvar
  public :: species_component
  public :: normalize_mass_fractions
  public :: multispecies_state_from_base
  public :: mass_fractions_from_state
  public :: synchronize_multispecies_thermodynamics
  public :: multispecies_state_is_physical
  public :: species_closure_error
  public :: thermodynamic_layout_error

contains

  pure integer function multispecies_nvar(nspecies) result(nvar)
    integer, intent(in) :: nspecies

    if (nspecies < 1 .or. nspecies > max_supported_species) then
      nvar = 0
    else
      nvar = nbase + nspecies
    end if
  end function multispecies_nvar

  pure integer function species_component(species) result(component)
    integer, intent(in) :: species

    if (species < 1 .or. species > max_supported_species) then
      component = 0
    else
      component = nbase + species
    end if
  end function species_component

  pure subroutine normalize_mass_fractions(mass_fractions, nspecies, ok)
    integer, intent(in) :: nspecies
    real(dp), intent(inout) :: mass_fractions(:)
    logical, intent(out) :: ok

    real(dp) :: total

    ok = .false.
    if (multispecies_nvar(nspecies) == 0) return
    if (size(mass_fractions) < nspecies) return
    if (any(mass_fractions(1:nspecies) < -species_negative_tolerance)) return
    if (any(mass_fractions(1:nspecies) > 1.0_dp + species_negative_tolerance)) return

    mass_fractions(1:nspecies) = max(0.0_dp, mass_fractions(1:nspecies))
    total = sum(mass_fractions(1:nspecies))
    if (total <= tiny(1.0_dp)) return

    mass_fractions(1:nspecies) = mass_fractions(1:nspecies) / total
    ok = .true.
  end subroutine normalize_mass_fractions

  pure subroutine multispecies_state_from_base( &
      base_state, mass_fractions, nspecies, gamma, state, ok)
    integer, intent(in) :: nspecies
    real(dp), intent(in) :: base_state(ncons), mass_fractions(:), gamma
    real(dp), intent(out) :: state(:)
    logical, intent(out) :: ok

    real(dp) :: normalized(max_supported_species)
    integer :: species, nvar

    state = 0.0_dp
    ok = .false.
    nvar = multispecies_nvar(nspecies)
    if (nvar == 0 .or. size(state) /= nvar) return
    if (size(mass_fractions) < nspecies) return
    if (base_state(irho) <= density_floor) return

    normalized = 0.0_dp
    normalized(1:nspecies) = mass_fractions(1:nspecies)
    call normalize_mass_fractions(normalized, nspecies, ok)
    if (.not. ok) return

    state(1:ncons) = base_state
    call synchronize_multispecies_thermodynamics(state, gamma, ok)
    if (.not. ok) return
    do species = 1, nspecies
      state(species_component(species)) = &
        base_state(irho) * normalized(species)
    end do
    ok = .true.
  end subroutine multispecies_state_from_base

  pure subroutine synchronize_multispecies_thermodynamics(state, gamma, ok)
    real(dp), intent(inout) :: state(:)
    real(dp), intent(in) :: gamma
    logical, intent(out) :: ok

    real(dp) :: primitive(nprim), rho, kinetic_energy_density

    ok = .false.
    if (size(state) < nbase) return
    call conserved_to_primitive(state(1:ncons), gamma, primitive, ok)
    if (.not. ok) return

    rho = state(irho)
    kinetic_energy_density = 0.5_dp * &
      (state(imx)**2 + state(imy)**2 + state(imz)**2) / rho
    state(iei) = state(iet) - kinetic_energy_density
    state(item) = primitive(qp) / primitive(qrho)
    ok = state(iei) > 0.0_dp .and. state(item) > 0.0_dp
  end subroutine synchronize_multispecies_thermodynamics

  pure subroutine mass_fractions_from_state( &
      state, nspecies, mass_fractions, ok, closure_error)
    integer, intent(in) :: nspecies
    real(dp), intent(in) :: state(:)
    real(dp), intent(out) :: mass_fractions(:)
    logical, intent(out) :: ok
    real(dp), intent(out), optional :: closure_error

    real(dp) :: rho, species_total, local_error
    integer :: species, nvar

    mass_fractions = 0.0_dp
    ok = .false.
    local_error = huge(1.0_dp)
    nvar = multispecies_nvar(nspecies)
    if (nvar == 0 .or. size(state) /= nvar) then
      if (present(closure_error)) closure_error = local_error
      return
    end if
    if (size(mass_fractions) < nspecies) then
      if (present(closure_error)) closure_error = local_error
      return
    end if

    rho = state(irho)
    if (rho <= density_floor) then
      if (present(closure_error)) closure_error = local_error
      return
    end if

    species_total = 0.0_dp
    do species = 1, nspecies
      if (state(species_component(species)) < &
          -species_negative_tolerance * max(1.0_dp, rho)) then
        if (present(closure_error)) closure_error = local_error
        return
      end if
      mass_fractions(species) = &
        max(0.0_dp, state(species_component(species))) / rho
      species_total = species_total + state(species_component(species))
    end do

    local_error = abs(species_total - rho) / max(rho, density_floor)
    if (present(closure_error)) closure_error = local_error
    if (local_error > species_closure_tolerance) return

    call normalize_mass_fractions(mass_fractions, nspecies, ok)
  end subroutine mass_fractions_from_state

  pure logical function multispecies_state_is_physical( &
      state, gamma, nspecies) result(is_physical)
    integer, intent(in) :: nspecies
    real(dp), intent(in) :: state(:), gamma

    real(dp) :: mass_fractions(max_supported_species)
    logical :: species_ok

    is_physical = .false.
    if (size(state) /= multispecies_nvar(nspecies)) return
    if (.not. state_is_physical(state(1:ncons), gamma)) return
    if (thermodynamic_layout_error(state, gamma) > &
        thermodynamic_layout_tolerance) return

    call mass_fractions_from_state( &
      state, nspecies, mass_fractions, species_ok)
    is_physical = species_ok
  end function multispecies_state_is_physical

  pure real(dp) function species_closure_error(state, nspecies) result(error)
    integer, intent(in) :: nspecies
    real(dp), intent(in) :: state(:)

    real(dp) :: rho, species_total
    integer :: species

    error = huge(1.0_dp)
    if (size(state) /= multispecies_nvar(nspecies)) return
    rho = state(irho)
    if (rho <= density_floor) return

    species_total = 0.0_dp
    do species = 1, nspecies
      species_total = species_total + state(species_component(species))
    end do
    error = abs(species_total - rho) / max(rho, density_floor)
  end function species_closure_error

  pure real(dp) function thermodynamic_layout_error(state, gamma) result(error)
    real(dp), intent(in) :: state(:), gamma

    real(dp) :: primitive(nprim), internal_energy_density, temperature
    real(dp) :: kinetic_energy_density, scale
    logical :: ok

    error = huge(1.0_dp)
    if (size(state) < nbase) return
    call conserved_to_primitive(state(1:ncons), gamma, primitive, ok)
    if (.not. ok) return

    kinetic_energy_density = 0.5_dp * &
      (state(imx)**2 + state(imy)**2 + state(imz)**2) / state(irho)
    internal_energy_density = state(iet) - kinetic_energy_density
    temperature = primitive(qp) / primitive(qrho)
    scale = max(1.0_dp, abs(internal_energy_density), abs(temperature))
    error = max(abs(state(iei) - internal_energy_density), &
      abs(state(item) - temperature)) / scale
  end function thermodynamic_layout_error

end module multispecies_state_mod
