module full_reactor_rhs_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use mixture_thermo_mod, only: temperature_from_internal_energy
  use pressure_dependent_kinetics_mod, only: &
    pressure_dependent_reaction, pressure_dependent_mass_fraction_rhs, &
    pressure_dependent_production_rates
  implicit none
  private

  public :: evaluate_full_reactor_rhs
  public :: normalize_reactor_composition

contains

  subroutine evaluate_full_reactor_rhs( &
      species, reactions, density, target_internal_energy, mass_fractions, &
      temperature_guess, derivative, temperature, production_rates, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(pressure_dependent_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: density, target_internal_energy
    real(dp), intent(in) :: mass_fractions(:), temperature_guess
    real(dp), intent(out) :: derivative(:), temperature
    real(dp), intent(out), optional :: production_rates(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: normalized(:), local_rates(:)
    integer :: iterations

    derivative = 0.0_dp
    temperature = temperature_guess
    ok = density > 0.0_dp
    ok = ok .and. size(species) == size(mass_fractions)
    ok = ok .and. size(derivative) == size(species)
    if (present(production_rates)) then
      ok = ok .and. size(production_rates) == size(species)
      production_rates = 0.0_dp
    end if
    if (.not. ok) return

    allocate(normalized(size(species)), local_rates(size(species)))
    call normalize_reactor_composition(mass_fractions, normalized, ok)
    if (.not. ok) return
    call temperature_from_internal_energy( &
      species, normalized, target_internal_energy, temperature_guess, &
      temperature, ok, iterations)
    if (.not. ok) return
    call pressure_dependent_mass_fraction_rhs( &
      species, reactions, density, temperature, normalized, derivative, ok)
    if (.not. ok) return
    if (present(production_rates)) then
      call pressure_dependent_production_rates( &
        species, reactions, density, temperature, normalized, local_rates, ok)
      if (.not. ok) return
      production_rates = local_rates
    end if
    ok = ieee_is_finite(temperature) .and. all(ieee_is_finite(derivative))
  end subroutine evaluate_full_reactor_rhs

  subroutine normalize_reactor_composition(input, output, ok)
    real(dp), intent(in) :: input(:)
    real(dp), intent(out) :: output(:)
    logical, intent(out) :: ok
    real(dp) :: total

    output = 0.0_dp
    ok = size(input) == size(output) .and. size(input) > 0
    if (ok) ok = all(ieee_is_finite(input))
    if (ok) ok = minval(input) >= -1.0e-11_dp
    if (.not. ok) return
    output = max(input, 0.0_dp)
    total = sum(output)
    ok = total > 0.0_dp .and. ieee_is_finite(total)
    if (.not. ok) return
    output = output / total
  end subroutine normalize_reactor_composition

end module full_reactor_rhs_mod
