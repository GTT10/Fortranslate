module pressure_dependent_kinetics_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species, nasa7_mass_properties
  implicit none
  private

  integer, parameter, public :: reaction_elementary = 1
  integer, parameter, public :: reaction_three_body = 2
  integer, parameter, public :: reaction_falloff_troe = 3

  real(dp), parameter :: universal_gas_constant = 8314.46261815324_dp
  real(dp), parameter :: standard_pressure = 101325.0_dp
  real(dp), parameter :: tiny_positive = 1.0e-300_dp

  type, public :: arrhenius_rate
    real(dp) :: pre_exponential = 0.0_dp
    real(dp) :: temperature_exponent = 0.0_dp
    real(dp) :: activation_temperature = 0.0_dp
  end type arrhenius_rate

  type, public :: pressure_dependent_reaction
    integer :: kind = reaction_elementary
    logical :: reversible = .true.
    real(dp), allocatable :: reactant_stoich(:)
    real(dp), allocatable :: product_stoich(:)
    real(dp), allocatable :: efficiencies(:)
    type(arrhenius_rate) :: high_rate
    type(arrhenius_rate) :: low_rate
    real(dp) :: troe_a = 0.0_dp
    real(dp) :: troe_t3 = 0.0_dp
    real(dp) :: troe_t1 = 0.0_dp
    real(dp) :: troe_t2 = 0.0_dp
  end type pressure_dependent_reaction

  public :: valid_pressure_dependent_reaction
  public :: arrhenius_rate_constant
  public :: effective_third_body_concentration
  public :: troe_broadening_factor
  public :: reaction_rate_constants
  public :: pressure_dependent_production_rates
  public :: pressure_dependent_mass_fraction_rhs

contains

  logical function valid_pressure_dependent_reaction(reaction, nspecies) result(valid)
    type(pressure_dependent_reaction), intent(in) :: reaction
    integer, intent(in) :: nspecies

    valid = nspecies > 0
    valid = valid .and. allocated(reaction%reactant_stoich)
    valid = valid .and. allocated(reaction%product_stoich)
    valid = valid .and. allocated(reaction%efficiencies)
    if (.not. valid) return
    valid = size(reaction%reactant_stoich) == nspecies
    valid = valid .and. size(reaction%product_stoich) == nspecies
    valid = valid .and. size(reaction%efficiencies) == nspecies
    valid = valid .and. all(reaction%reactant_stoich >= 0.0_dp)
    valid = valid .and. all(reaction%product_stoich >= 0.0_dp)
    valid = valid .and. all(reaction%efficiencies >= 0.0_dp)
    valid = valid .and. sum(reaction%reactant_stoich) > 0.0_dp
    valid = valid .and. sum(reaction%product_stoich) > 0.0_dp
    valid = valid .and. reaction%high_rate%pre_exponential > 0.0_dp
    valid = valid .and. reaction%kind >= reaction_elementary
    valid = valid .and. reaction%kind <= reaction_falloff_troe
    if (reaction%kind == reaction_falloff_troe) then
      valid = valid .and. reaction%low_rate%pre_exponential > 0.0_dp
      valid = valid .and. reaction%troe_a > 0.0_dp
      valid = valid .and. reaction%troe_a < 1.0_dp
      valid = valid .and. reaction%troe_t3 > 0.0_dp
      valid = valid .and. reaction%troe_t1 > 0.0_dp
      valid = valid .and. reaction%troe_t2 >= 0.0_dp
    end if
  end function valid_pressure_dependent_reaction

  real(dp) function arrhenius_rate_constant(rate, temperature, ok) result(value)
    type(arrhenius_rate), intent(in) :: rate
    real(dp), intent(in) :: temperature
    logical, intent(out) :: ok
    real(dp) :: logarithm

    ok = temperature > 0.0_dp .and. rate%pre_exponential > 0.0_dp
    if (.not. ok) then
      value = 0.0_dp
      return
    end if
    logarithm = log(rate%pre_exponential) + &
      rate%temperature_exponent * log(temperature) - &
      rate%activation_temperature / temperature
    if (logarithm > log(huge(1.0_dp)) - 4.0_dp) then
      ok = .false.
      value = 0.0_dp
      return
    end if
    value = exp(max(logarithm, log(tiny_positive)))
    ok = ieee_is_finite(value) .and. value >= 0.0_dp
  end function arrhenius_rate_constant

  real(dp) function effective_third_body_concentration( &
      reaction, concentrations, ok) result(concentration)
    type(pressure_dependent_reaction), intent(in) :: reaction
    real(dp), intent(in) :: concentrations(:)
    logical, intent(out) :: ok

    ok = allocated(reaction%efficiencies)
    if (ok) ok = size(reaction%efficiencies) == size(concentrations)
    if (ok) ok = all(concentrations >= 0.0_dp)
    if (.not. ok) then
      concentration = 0.0_dp
      return
    end if
    concentration = dot_product(reaction%efficiencies, concentrations)
    ok = ieee_is_finite(concentration) .and. concentration >= 0.0_dp
  end function effective_third_body_concentration

  real(dp) function troe_broadening_factor( &
      reaction, temperature, reduced_pressure, ok) result(factor)
    type(pressure_dependent_reaction), intent(in) :: reaction
    real(dp), intent(in) :: temperature, reduced_pressure
    logical, intent(out) :: ok
    real(dp) :: fcent, log_fcent, log_pr, c_value, n_value, denominator, f_value

    ok = reaction%kind == reaction_falloff_troe
    ok = ok .and. temperature > 0.0_dp .and. reduced_pressure >= 0.0_dp
    ok = ok .and. reaction%troe_t3 > 0.0_dp .and. reaction%troe_t1 > 0.0_dp
    if (.not. ok) then
      factor = 0.0_dp
      return
    end if

    fcent = (1.0_dp - reaction%troe_a) * &
      exp(-temperature / reaction%troe_t3) + &
      reaction%troe_a * exp(-temperature / reaction%troe_t1)
    if (reaction%troe_t2 > 0.0_dp) then
      fcent = fcent + exp(-reaction%troe_t2 / temperature)
    end if
    fcent = max(fcent, tiny_positive)
    log_fcent = log10(fcent)
    log_pr = log10(max(reduced_pressure, tiny_positive))
    c_value = -0.4_dp - 0.67_dp * log_fcent
    n_value = 0.75_dp - 1.27_dp * log_fcent
    denominator = n_value - 0.14_dp * (log_pr + c_value)
    if (abs(denominator) < 100.0_dp * epsilon(1.0_dp)) then
      ok = .false.
      factor = 0.0_dp
      return
    end if
    f_value = (log_pr + c_value) / denominator
    factor = 10.0_dp ** (log_fcent / (1.0_dp + f_value * f_value))
    ok = ieee_is_finite(factor) .and. factor > 0.0_dp
  end function troe_broadening_factor

  subroutine reaction_rate_constants( &
      species, reaction, temperature, concentrations, forward_rate, &
      reverse_rate, collider, reduced_pressure, broadening, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(pressure_dependent_reaction), intent(in) :: reaction
    real(dp), intent(in) :: temperature, concentrations(:)
    real(dp), intent(out) :: forward_rate, reverse_rate
    real(dp), intent(out) :: collider, reduced_pressure, broadening
    logical, intent(out) :: ok
    real(dp) :: high_rate, low_rate, equilibrium_constant
    logical :: local_ok

    forward_rate = 0.0_dp
    reverse_rate = 0.0_dp
    collider = 1.0_dp
    reduced_pressure = 0.0_dp
    broadening = 1.0_dp
    ok = size(species) == size(concentrations)
    if (ok) ok = valid_pressure_dependent_reaction(reaction, size(species))
    if (.not. ok) return

    high_rate = arrhenius_rate_constant(reaction%high_rate, temperature, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if

    select case (reaction%kind)
    case (reaction_elementary)
      forward_rate = high_rate
    case (reaction_three_body)
      collider = effective_third_body_concentration(reaction, concentrations, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      forward_rate = high_rate * collider
    case (reaction_falloff_troe)
      collider = effective_third_body_concentration(reaction, concentrations, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      low_rate = arrhenius_rate_constant(reaction%low_rate, temperature, local_ok)
      if (.not. local_ok .or. high_rate <= 0.0_dp) then
        ok = .false.
        return
      end if
      reduced_pressure = low_rate * collider / high_rate
      broadening = troe_broadening_factor( &
        reaction, temperature, reduced_pressure, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      forward_rate = high_rate * reduced_pressure / &
        (1.0_dp + reduced_pressure) * broadening
    case default
      ok = .false.
      return
    end select

    if (reaction%reversible) then
      equilibrium_constant = concentration_equilibrium_constant( &
        species, reaction, temperature, local_ok)
      if (.not. local_ok .or. equilibrium_constant <= 0.0_dp) then
        ok = .false.
        return
      end if
      reverse_rate = forward_rate / equilibrium_constant
    end if
    ok = ieee_is_finite(forward_rate) .and. ieee_is_finite(reverse_rate)
    ok = ok .and. forward_rate >= 0.0_dp .and. reverse_rate >= 0.0_dp
  end subroutine reaction_rate_constants

  subroutine pressure_dependent_production_rates( &
      species, reactions, density, temperature, mass_fractions, &
      production_rates, ok, forward_progress, reverse_progress)
    type(nasa7_species), intent(in) :: species(:)
    type(pressure_dependent_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: density, temperature, mass_fractions(:)
    real(dp), intent(out) :: production_rates(:)
    logical, intent(out) :: ok
    real(dp), intent(out), optional :: forward_progress(:), reverse_progress(:)
    real(dp), allocatable :: concentrations(:)
    real(dp) :: kf, kr, collider, reduced_pressure, broadening, qf, qr
    logical :: local_ok
    integer :: reaction_index, species_index

    production_rates = 0.0_dp
    ok = density > 0.0_dp .and. temperature > 0.0_dp
    ok = ok .and. size(species) == size(mass_fractions)
    ok = ok .and. size(production_rates) == size(species)
    ok = ok .and. all(mass_fractions >= -1.0e-13_dp)
    ok = ok .and. abs(sum(mass_fractions) - 1.0_dp) <= 1.0e-8_dp
    if (present(forward_progress)) ok = ok .and. &
      size(forward_progress) == size(reactions)
    if (present(reverse_progress)) ok = ok .and. &
      size(reverse_progress) == size(reactions)
    if (.not. ok) return

    allocate(concentrations(size(species)))
    do species_index = 1, size(species)
      if (species(species_index)%molecular_weight <= 0.0_dp) then
        ok = .false.
        return
      end if
      concentrations(species_index) = density * &
        max(0.0_dp, mass_fractions(species_index)) / &
        species(species_index)%molecular_weight
    end do

    do reaction_index = 1, size(reactions)
      call reaction_rate_constants( &
        species, reactions(reaction_index), temperature, concentrations, &
        kf, kr, collider, reduced_pressure, broadening, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      qf = kf * concentration_product( &
        concentrations, reactions(reaction_index)%reactant_stoich)
      qr = kr * concentration_product( &
        concentrations, reactions(reaction_index)%product_stoich)
      if (present(forward_progress)) forward_progress(reaction_index) = qf
      if (present(reverse_progress)) reverse_progress(reaction_index) = qr
      do species_index = 1, size(species)
        production_rates(species_index) = production_rates(species_index) + &
          (reactions(reaction_index)%product_stoich(species_index) - &
           reactions(reaction_index)%reactant_stoich(species_index)) * (qf - qr)
      end do
    end do
    ok = all(ieee_is_finite(production_rates))
  end subroutine pressure_dependent_production_rates

  subroutine pressure_dependent_mass_fraction_rhs( &
      species, reactions, density, temperature, mass_fractions, derivative, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(pressure_dependent_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: density, temperature, mass_fractions(:)
    real(dp), intent(out) :: derivative(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: production_rates(:)
    integer :: species_index

    derivative = 0.0_dp
    ok = size(derivative) == size(species)
    if (.not. ok) return
    allocate(production_rates(size(species)))
    call pressure_dependent_production_rates( &
      species, reactions, density, temperature, mass_fractions, &
      production_rates, ok)
    if (.not. ok) return
    do species_index = 1, size(species)
      derivative(species_index) = species(species_index)%molecular_weight * &
        production_rates(species_index) / density
    end do
    ok = all(ieee_is_finite(derivative))
  end subroutine pressure_dependent_mass_fraction_rhs

  real(dp) function concentration_product(concentrations, orders) result(product)
    real(dp), intent(in) :: concentrations(:), orders(:)
    integer :: index

    product = 1.0_dp
    do index = 1, size(concentrations)
      if (orders(index) <= 0.0_dp) cycle
      if (concentrations(index) <= 0.0_dp) then
        product = 0.0_dp
        return
      end if
      product = product * concentrations(index) ** orders(index)
    end do
  end function concentration_product

  real(dp) function concentration_equilibrium_constant( &
      species, reaction, temperature, ok) result(constant)
    type(nasa7_species), intent(in) :: species(:)
    type(pressure_dependent_reaction), intent(in) :: reaction
    real(dp), intent(in) :: temperature
    logical, intent(out) :: ok
    real(dp) :: cp, cv, enthalpy, internal_energy, entropy
    real(dp) :: delta_gibbs, delta_nu, logarithm
    logical :: local_ok
    integer :: index

    constant = 0.0_dp
    ok = temperature > 0.0_dp .and. size(species) == &
      size(reaction%reactant_stoich)
    if (.not. ok) return
    delta_gibbs = 0.0_dp
    do index = 1, size(species)
      call nasa7_mass_properties( &
        species(index), temperature, cp, cv, enthalpy, internal_energy, &
        entropy, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      delta_gibbs = delta_gibbs + &
        (reaction%product_stoich(index) - reaction%reactant_stoich(index)) * &
        (enthalpy - temperature * entropy) * species(index)%molecular_weight
    end do
    delta_nu = sum(reaction%product_stoich - reaction%reactant_stoich)
    logarithm = -delta_gibbs / (universal_gas_constant * temperature) + &
      delta_nu * log(standard_pressure / &
      (universal_gas_constant * temperature))
    if (logarithm > 700.0_dp .or. logarithm < -700.0_dp) then
      constant = exp(max(-700.0_dp, min(700.0_dp, logarithm)))
    else
      constant = exp(logarithm)
    end if
    ok = ieee_is_finite(constant) .and. constant > 0.0_dp
  end function concentration_equilibrium_constant

end module pressure_dependent_kinetics_mod
