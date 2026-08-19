module elementary_kinetics_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: &
    nasa7_species, nasa7_dimensionless_properties, universal_gas_constant
  use mixture_thermo_mod, only: valid_mixture_composition
  implicit none
  private

  real(dp), parameter, public :: kinetics_standard_pressure = 101325.0_dp
  real(dp), parameter :: logarithm_limit = 700.0_dp
  real(dp), parameter :: minimum_reduced_pressure = 1.0e-300_dp

  integer, parameter, public :: reaction_kind_elementary = 1
  integer, parameter, public :: reaction_kind_three_body = 2
  integer, parameter, public :: reaction_kind_falloff = 3

  type, public :: arrhenius_rate
    real(dp) :: pre_exponential = 0.0_dp
    real(dp) :: temperature_exponent = 0.0_dp
    real(dp) :: activation_energy = 0.0_dp ! J / kmol
  end type arrhenius_rate

  type, public :: troe_parameters
    logical :: enabled = .false.
    real(dp) :: alpha = 0.0_dp
    real(dp) :: temperature_3 = 0.0_dp
    real(dp) :: temperature_1 = 0.0_dp
    real(dp) :: temperature_2 = 0.0_dp
  end type troe_parameters

  type, public :: elementary_reaction
    character(len=128) :: equation = ""
    integer :: kind = reaction_kind_elementary
    real(dp), allocatable :: reactant_stoich(:)
    real(dp), allocatable :: product_stoich(:)
    type(arrhenius_rate) :: forward_rate
    type(arrhenius_rate) :: low_pressure_rate
    type(arrhenius_rate) :: high_pressure_rate
    real(dp), allocatable :: third_body_efficiencies(:)
    type(troe_parameters) :: troe
    logical :: reversible = .true.
  end type elementary_reaction

  public :: valid_arrhenius_rate
  public :: valid_troe_parameters
  public :: valid_elementary_reaction
  public :: arrhenius_rate_constant
  public :: reaction_equilibrium_constant
  public :: effective_third_body_concentration
  public :: troe_falloff_factor
  public :: reaction_effective_forward_rate
  public :: reaction_progress_rate
  public :: molar_concentrations_from_mass_fractions
  public :: elementary_production_rates
  public :: elementary_production_rates_from_concentrations
  public :: elementary_production_jacobian
  public :: elementary_mass_fraction_rhs
  public :: elementary_mass_fraction_jacobian

contains

  logical function valid_arrhenius_rate(rate) result(valid)
    type(arrhenius_rate), intent(in) :: rate

    valid = rate%pre_exponential >= 0.0_dp .and. &
      all(ieee_is_finite([rate%pre_exponential, &
        rate%temperature_exponent, rate%activation_energy]))
  end function valid_arrhenius_rate

  logical function valid_troe_parameters(parameters) result(valid)
    type(troe_parameters), intent(in) :: parameters

    if (.not. parameters%enabled) then
      valid = .true.
      return
    end if
    valid = parameters%alpha > 0.0_dp .and. parameters%alpha < 1.0_dp .and. &
      parameters%temperature_3 > 0.0_dp .and. &
      parameters%temperature_1 > 0.0_dp .and. &
      parameters%temperature_2 >= 0.0_dp .and. &
      all(ieee_is_finite([parameters%alpha, parameters%temperature_3, &
        parameters%temperature_1, parameters%temperature_2]))
  end function valid_troe_parameters

  logical function valid_elementary_reaction(reaction, nspecies) result(valid)
    type(elementary_reaction), intent(in) :: reaction
    integer, intent(in) :: nspecies

    valid = .false.
    if (nspecies < 1) return
    if (.not. allocated(reaction%reactant_stoich)) return
    if (.not. allocated(reaction%product_stoich)) return
    if (size(reaction%reactant_stoich) /= nspecies) return
    if (size(reaction%product_stoich) /= nspecies) return
    if (any(.not. ieee_is_finite(reaction%reactant_stoich))) return
    if (any(.not. ieee_is_finite(reaction%product_stoich))) return
    if (any(reaction%reactant_stoich < 0.0_dp)) return
    if (any(reaction%product_stoich < 0.0_dp)) return
    if (sum(reaction%reactant_stoich) <= 0.0_dp) return
    if (sum(reaction%product_stoich) <= 0.0_dp) return
    if (maxval(abs(reaction%product_stoich - &
        reaction%reactant_stoich)) <= 0.0_dp) return

    select case (reaction%kind)
    case (reaction_kind_elementary)
      if (.not. valid_arrhenius_rate(reaction%forward_rate)) return
    case (reaction_kind_three_body)
      if (.not. valid_arrhenius_rate(reaction%forward_rate)) return
      if (.not. valid_efficiencies(reaction, nspecies)) return
    case (reaction_kind_falloff)
      if (.not. valid_arrhenius_rate(reaction%low_pressure_rate)) return
      if (.not. valid_arrhenius_rate(reaction%high_pressure_rate)) return
      if (.not. valid_efficiencies(reaction, nspecies)) return
      if (.not. valid_troe_parameters(reaction%troe)) return
    case default
      return
    end select
    valid = .true.
  end function valid_elementary_reaction

  logical function valid_efficiencies(reaction, nspecies) result(valid)
    type(elementary_reaction), intent(in) :: reaction
    integer, intent(in) :: nspecies

    valid = allocated(reaction%third_body_efficiencies)
    if (.not. valid) return
    valid = size(reaction%third_body_efficiencies) == nspecies
    if (.not. valid) return
    valid = all(ieee_is_finite(reaction%third_body_efficiencies)) .and. &
      all(reaction%third_body_efficiencies >= 0.0_dp)
  end function valid_efficiencies

  real(dp) function arrhenius_rate_constant(rate, temperature, ok) &
      result(rate_constant)
    type(arrhenius_rate), intent(in) :: rate
    real(dp), intent(in) :: temperature
    logical, intent(out) :: ok

    real(dp) :: logarithm

    rate_constant = 0.0_dp
    ok = valid_arrhenius_rate(rate) .and. temperature > 0.0_dp .and. &
      ieee_is_finite(temperature)
    if (.not. ok) return
    if (rate%pre_exponential <= tiny(1.0_dp)) return

    logarithm = log(rate%pre_exponential) + &
      rate%temperature_exponent * log(temperature) - &
      rate%activation_energy / (universal_gas_constant * temperature)
    if (.not. ieee_is_finite(logarithm) .or. &
        logarithm > logarithm_limit) then
      ok = .false.
      return
    end if
    if (logarithm < -logarithm_limit) then
      rate_constant = 0.0_dp
    else
      rate_constant = exp(logarithm)
    end if
    ok = ieee_is_finite(rate_constant) .and. rate_constant >= 0.0_dp
  end function arrhenius_rate_constant

  real(dp) function reaction_equilibrium_constant( &
      species, reaction, temperature, ok) result(equilibrium_constant)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reaction
    real(dp), intent(in) :: temperature
    logical, intent(out) :: ok

    real(dp) :: cp_over_r, h_over_rt, s_over_r
    real(dp) :: delta_g_over_rt, delta_nu, log_kc
    real(dp) :: reference_concentration, net_stoich
    logical :: species_ok
    integer :: i

    equilibrium_constant = 0.0_dp
    ok = valid_elementary_reaction(reaction, size(species)) .and. &
      temperature > 0.0_dp
    if (.not. ok) return

    delta_g_over_rt = 0.0_dp
    delta_nu = 0.0_dp
    do i = 1, size(species)
      call nasa7_dimensionless_properties( &
        species(i), temperature, cp_over_r, h_over_rt, s_over_r, species_ok)
      if (.not. species_ok) then
        ok = .false.
        return
      end if
      net_stoich = reaction%product_stoich(i) - &
        reaction%reactant_stoich(i)
      delta_g_over_rt = delta_g_over_rt + &
        net_stoich * (h_over_rt - s_over_r)
      delta_nu = delta_nu + net_stoich
    end do

    reference_concentration = kinetics_standard_pressure / &
      (universal_gas_constant * temperature)
    log_kc = -delta_g_over_rt + delta_nu * log(reference_concentration)
    if (.not. ieee_is_finite(log_kc) .or. abs(log_kc) > logarithm_limit) then
      ok = .false.
      return
    end if
    equilibrium_constant = exp(log_kc)
    ok = ieee_is_finite(equilibrium_constant) .and. &
      equilibrium_constant > 0.0_dp
  end function reaction_equilibrium_constant

  real(dp) function effective_third_body_concentration( &
      reaction, concentrations, ok) result(effective_concentration)
    type(elementary_reaction), intent(in) :: reaction
    real(dp), intent(in) :: concentrations(:)
    logical, intent(out) :: ok

    effective_concentration = 0.0_dp
    ok = allocated(reaction%third_body_efficiencies)
    if (.not. ok) return
    ok = size(reaction%third_body_efficiencies) == size(concentrations) .and. &
      all(ieee_is_finite(concentrations)) .and. &
      all(concentrations >= 0.0_dp)
    if (.not. ok) return
    effective_concentration = sum( &
      reaction%third_body_efficiencies * concentrations)
    ok = ieee_is_finite(effective_concentration) .and. &
      effective_concentration >= 0.0_dp
  end function effective_third_body_concentration

  real(dp) function troe_falloff_factor( &
      parameters, temperature, reduced_pressure, ok, &
      logarithmic_derivative) result(falloff_factor)
    type(troe_parameters), intent(in) :: parameters
    real(dp), intent(in) :: temperature, reduced_pressure
    logical, intent(out) :: ok
    real(dp), intent(out), optional :: logarithmic_derivative

    real(dp) :: f_center, log_f_center, log_pressure
    real(dp) :: c_parameter, n_parameter, numerator, denominator
    real(dp) :: broadening_argument, log_factor
    real(dp) :: derivative_argument, derivative_log_factor

    falloff_factor = 1.0_dp
    if (present(logarithmic_derivative)) logarithmic_derivative = 0.0_dp
    ok = temperature > 0.0_dp .and. reduced_pressure >= 0.0_dp .and. &
      ieee_is_finite(temperature) .and. ieee_is_finite(reduced_pressure) .and. &
      valid_troe_parameters(parameters)
    if (.not. ok) return
    if (.not. parameters%enabled) return

    f_center = (1.0_dp - parameters%alpha) * &
      exp(-temperature / parameters%temperature_3) + &
      parameters%alpha * exp(-temperature / parameters%temperature_1)
    if (parameters%temperature_2 > 0.0_dp) then
      f_center = f_center + exp(-parameters%temperature_2 / temperature)
    end if
    if (f_center <= 0.0_dp .or. .not. ieee_is_finite(f_center)) then
      ok = .false.
      return
    end if

    log_f_center = log10(f_center)
    log_pressure = log10(max(reduced_pressure, minimum_reduced_pressure))
    c_parameter = -0.4_dp - 0.67_dp * log_f_center
    n_parameter = 0.75_dp - 1.27_dp * log_f_center
    numerator = log_pressure + c_parameter
    denominator = n_parameter - 0.14_dp * numerator
    if (abs(denominator) <= tiny(1.0_dp)) then
      ok = .false.
      return
    end if
    broadening_argument = numerator / denominator
    log_factor = log_f_center / &
      (1.0_dp + broadening_argument * broadening_argument)
    falloff_factor = 10.0_dp**log_factor

    if (present(logarithmic_derivative)) then
      derivative_argument = n_parameter / (denominator * denominator)
      derivative_log_factor = -2.0_dp * log_f_center * &
        broadening_argument * derivative_argument / &
        (1.0_dp + broadening_argument * broadening_argument)**2
      logarithmic_derivative = derivative_log_factor
    end if
    ok = ieee_is_finite(falloff_factor) .and. falloff_factor > 0.0_dp
  end function troe_falloff_factor

  subroutine reaction_effective_forward_rate( &
      reaction, temperature, concentrations, rate_constant, &
      concentration_derivative, ok)
    type(elementary_reaction), intent(in) :: reaction
    real(dp), intent(in) :: temperature, concentrations(:)
    real(dp), intent(out) :: rate_constant
    real(dp), intent(out), optional :: concentration_derivative(:)
    logical, intent(out) :: ok

    real(dp) :: base_rate, low_rate, high_rate
    real(dp) :: third_body, reduced_pressure, falloff_factor
    real(dp) :: log_factor_derivative, derivative_wrt_pressure
    real(dp) :: derivative_wrt_third_body
    integer :: i

    rate_constant = 0.0_dp
    if (present(concentration_derivative)) concentration_derivative = 0.0_dp
    ok = size(concentrations) > 0 .and. &
      all(ieee_is_finite(concentrations)) .and. all(concentrations >= 0.0_dp)
    if (present(concentration_derivative)) then
      ok = ok .and. size(concentration_derivative) == size(concentrations)
    end if
    if (.not. ok) return

    select case (reaction%kind)
    case (reaction_kind_elementary)
      rate_constant = arrhenius_rate_constant( &
        reaction%forward_rate, temperature, ok)
    case (reaction_kind_three_body)
      if (.not. valid_efficiencies(reaction, size(concentrations))) then
        ok = .false.
        return
      end if
      base_rate = arrhenius_rate_constant( &
        reaction%forward_rate, temperature, ok)
      if (.not. ok) return
      third_body = effective_third_body_concentration( &
        reaction, concentrations, ok)
      if (.not. ok) return
      rate_constant = base_rate * third_body
      if (present(concentration_derivative)) then
        concentration_derivative = base_rate * &
          reaction%third_body_efficiencies
      end if
    case (reaction_kind_falloff)
      if (.not. valid_efficiencies(reaction, size(concentrations))) then
        ok = .false.
        return
      end if
      low_rate = arrhenius_rate_constant( &
        reaction%low_pressure_rate, temperature, ok)
      if (.not. ok) return
      high_rate = arrhenius_rate_constant( &
        reaction%high_pressure_rate, temperature, ok)
      if (.not. ok .or. high_rate <= 0.0_dp) then
        ok = .false.
        return
      end if
      third_body = effective_third_body_concentration( &
        reaction, concentrations, ok)
      if (.not. ok) return
      reduced_pressure = low_rate * third_body / high_rate
      falloff_factor = troe_falloff_factor( &
        reaction%troe, temperature, reduced_pressure, ok, &
        log_factor_derivative)
      if (.not. ok) return
      rate_constant = high_rate * reduced_pressure / &
        (1.0_dp + reduced_pressure) * falloff_factor
      if (present(concentration_derivative)) then
        if (third_body <= tiny(1.0_dp)) then
          derivative_wrt_third_body = low_rate * falloff_factor
        else
          derivative_wrt_pressure = high_rate * falloff_factor * &
            (1.0_dp / (1.0_dp + reduced_pressure)**2 + &
             log_factor_derivative / (1.0_dp + reduced_pressure))
          derivative_wrt_third_body = derivative_wrt_pressure * &
            low_rate / high_rate
        end if
        do i = 1, size(concentrations)
          concentration_derivative(i) = derivative_wrt_third_body * &
            reaction%third_body_efficiencies(i)
        end do
      end if
    case default
      ok = .false.
    end select
    ok = ok .and. ieee_is_finite(rate_constant) .and. rate_constant >= 0.0_dp
    if (present(concentration_derivative)) then
      ok = ok .and. all(ieee_is_finite(concentration_derivative))
    end if
  end subroutine reaction_effective_forward_rate

  subroutine reaction_progress_rate( &
      species, reaction, temperature, concentrations, forward_progress, &
      reverse_progress, net_progress, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reaction
    real(dp), intent(in) :: temperature, concentrations(:)
    real(dp), intent(out) :: forward_progress, reverse_progress, net_progress
    logical, intent(out) :: ok

    real(dp) :: forward_constant, reverse_constant, equilibrium_constant
    real(dp) :: reactant_product, product_product
    real(dp), allocatable :: unused_derivative(:)

    forward_progress = 0.0_dp
    reverse_progress = 0.0_dp
    net_progress = 0.0_dp
    ok = valid_elementary_reaction(reaction, size(species)) .and. &
      size(concentrations) == size(species) .and. &
      all(ieee_is_finite(concentrations)) .and. &
      all(concentrations >= 0.0_dp)
    if (.not. ok) return

    allocate(unused_derivative(size(species)))
    call reaction_effective_forward_rate( &
      reaction, temperature, concentrations, forward_constant, &
      unused_derivative, ok)
    if (.not. ok) return
    call concentration_product( &
      concentrations, reaction%reactant_stoich, reactant_product, ok)
    if (.not. ok) return
    forward_progress = forward_constant * reactant_product

    if (reaction%reversible) then
      equilibrium_constant = reaction_equilibrium_constant( &
        species, reaction, temperature, ok)
      if (.not. ok) return
      reverse_constant = forward_constant / equilibrium_constant
      call concentration_product( &
        concentrations, reaction%product_stoich, product_product, ok)
      if (.not. ok) return
      reverse_progress = reverse_constant * product_product
    end if

    net_progress = forward_progress - reverse_progress
    ok = all(ieee_is_finite( &
      [forward_progress, reverse_progress, net_progress]))
  end subroutine reaction_progress_rate

  subroutine concentration_product(concentrations, stoich, product, ok)
    real(dp), intent(in) :: concentrations(:), stoich(:)
    real(dp), intent(out) :: product
    logical, intent(out) :: ok
    integer :: i

    product = 1.0_dp
    ok = size(concentrations) == size(stoich) .and. &
      all(ieee_is_finite(concentrations)) .and. &
      all(ieee_is_finite(stoich)) .and. &
      all(concentrations >= 0.0_dp) .and. all(stoich >= 0.0_dp)
    if (.not. ok) return
    do i = 1, size(concentrations)
      if (stoich(i) > 0.0_dp) product = product * &
        concentrations(i)**stoich(i)
    end do
    ok = ieee_is_finite(product) .and. product >= 0.0_dp
  end subroutine concentration_product

  subroutine concentration_product_gradient( &
      concentrations, stoich, product, gradient, ok)
    real(dp), intent(in) :: concentrations(:), stoich(:)
    real(dp), intent(out) :: product, gradient(:)
    logical, intent(out) :: ok
    real(dp) :: factor
    integer :: i, j

    gradient = 0.0_dp
    call concentration_product(concentrations, stoich, product, ok)
    if (.not. ok .or. size(gradient) /= size(concentrations)) then
      ok = .false.
      return
    end if

    do j = 1, size(concentrations)
      if (stoich(j) <= 0.0_dp) cycle
      if (concentrations(j) > 0.0_dp) then
        gradient(j) = product * stoich(j) / concentrations(j)
      else
        if (stoich(j) < 1.0_dp) then
          ok = .false.
          return
        else if (abs(stoich(j) - 1.0_dp) <= epsilon(1.0_dp)) then
          factor = 1.0_dp
          do i = 1, size(concentrations)
            if (i == j .or. stoich(i) <= 0.0_dp) cycle
            factor = factor * concentrations(i)**stoich(i)
          end do
          gradient(j) = factor
        else
          gradient(j) = 0.0_dp
        end if
      end if
    end do
    ok = all(ieee_is_finite(gradient))
  end subroutine concentration_product_gradient

  subroutine molar_concentrations_from_mass_fractions( &
      species, density, mass_fractions, concentrations, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: density, mass_fractions(:)
    real(dp), intent(out) :: concentrations(:)
    logical, intent(out) :: ok

    integer :: i

    concentrations = 0.0_dp
    ok = density > 0.0_dp .and. size(concentrations) == size(species)
    if (.not. ok) return
    ok = valid_mixture_composition(species, mass_fractions)
    if (.not. ok) return
    do i = 1, size(species)
      concentrations(i) = density * max(0.0_dp, mass_fractions(i)) / &
        species(i)%molecular_weight
    end do
    ok = all(ieee_is_finite(concentrations)) .and. &
      all(concentrations >= 0.0_dp)
  end subroutine molar_concentrations_from_mass_fractions

  subroutine elementary_production_rates( &
      species, reactions, temperature, density, mass_fractions, &
      molar_production_rates, ok, forward_progress, reverse_progress)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: temperature, density, mass_fractions(:)
    real(dp), intent(out) :: molar_production_rates(:)
    logical, intent(out) :: ok
    real(dp), intent(out), optional :: forward_progress(:), reverse_progress(:)

    real(dp), allocatable :: concentrations(:)

    molar_production_rates = 0.0_dp
    ok = size(molar_production_rates) == size(species) .and. &
      size(reactions) >= 1
    if (.not. ok) return
    allocate(concentrations(size(species)))
    call molar_concentrations_from_mass_fractions( &
      species, density, mass_fractions, concentrations, ok)
    if (.not. ok) return
    call elementary_production_rates_from_concentrations( &
      species, reactions, temperature, concentrations, &
      molar_production_rates, ok, forward_progress, reverse_progress)
  end subroutine elementary_production_rates

  subroutine elementary_production_rates_from_concentrations( &
      species, reactions, temperature, concentrations, &
      molar_production_rates, ok, forward_progress, reverse_progress)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: temperature, concentrations(:)
    real(dp), intent(out) :: molar_production_rates(:)
    logical, intent(out) :: ok
    real(dp), intent(out), optional :: forward_progress(:), reverse_progress(:)

    real(dp) :: q_forward, q_reverse, q_net
    logical :: reaction_ok
    integer :: i, reaction_index

    molar_production_rates = 0.0_dp
    ok = size(molar_production_rates) == size(species) .and. &
      size(concentrations) == size(species) .and. size(reactions) >= 1 .and. &
      all(ieee_is_finite(concentrations)) .and. &
      all(concentrations >= 0.0_dp)
    if (.not. ok) return
    if (present(forward_progress)) then
      if (size(forward_progress) /= size(reactions)) then
        ok = .false.
        return
      end if
      forward_progress = 0.0_dp
    end if
    if (present(reverse_progress)) then
      if (size(reverse_progress) /= size(reactions)) then
        ok = .false.
        return
      end if
      reverse_progress = 0.0_dp
    end if

    do reaction_index = 1, size(reactions)
      call reaction_progress_rate( &
        species, reactions(reaction_index), temperature, concentrations, &
        q_forward, q_reverse, q_net, reaction_ok)
      if (.not. reaction_ok) then
        ok = .false.
        return
      end if
      if (present(forward_progress)) &
        forward_progress(reaction_index) = q_forward
      if (present(reverse_progress)) &
        reverse_progress(reaction_index) = q_reverse
      do i = 1, size(species)
        molar_production_rates(i) = molar_production_rates(i) + &
          (reactions(reaction_index)%product_stoich(i) - &
           reactions(reaction_index)%reactant_stoich(i)) * q_net
      end do
    end do
    ok = all(ieee_is_finite(molar_production_rates))
  end subroutine elementary_production_rates_from_concentrations

  subroutine elementary_production_jacobian( &
      species, reactions, temperature, concentrations, jacobian, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: temperature, concentrations(:)
    real(dp), intent(out) :: jacobian(:, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: rate_derivative(:)
    real(dp), allocatable :: reactant_gradient(:), product_gradient(:)
    real(dp), allocatable :: progress_derivative(:)
    real(dp) :: forward_constant, reverse_constant, equilibrium_constant
    real(dp) :: reactant_product, product_product, net_stoich
    logical :: reaction_ok
    integer :: i, j, reaction_index, nspecies

    jacobian = 0.0_dp
    nspecies = size(species)
    ok = size(concentrations) == nspecies .and. &
      size(jacobian, 1) == nspecies .and. size(jacobian, 2) == nspecies .and. &
      all(ieee_is_finite(concentrations)) .and. all(concentrations >= 0.0_dp)
    if (.not. ok) return
    allocate(rate_derivative(nspecies), reactant_gradient(nspecies))
    allocate(product_gradient(nspecies), progress_derivative(nspecies))

    do reaction_index = 1, size(reactions)
      if (.not. valid_elementary_reaction( &
          reactions(reaction_index), nspecies)) then
        ok = .false.
        return
      end if
      call reaction_effective_forward_rate( &
        reactions(reaction_index), temperature, concentrations, &
        forward_constant, rate_derivative, reaction_ok)
      if (.not. reaction_ok) then
        ok = .false.
        return
      end if
      call concentration_product_gradient( &
        concentrations, reactions(reaction_index)%reactant_stoich, &
        reactant_product, reactant_gradient, reaction_ok)
      if (.not. reaction_ok) then
        ok = .false.
        return
      end if
      progress_derivative = rate_derivative * reactant_product + &
        forward_constant * reactant_gradient

      if (reactions(reaction_index)%reversible) then
        equilibrium_constant = reaction_equilibrium_constant( &
          species, reactions(reaction_index), temperature, reaction_ok)
        if (.not. reaction_ok) then
          ok = .false.
          return
        end if
        reverse_constant = forward_constant / equilibrium_constant
        call concentration_product_gradient( &
          concentrations, reactions(reaction_index)%product_stoich, &
          product_product, product_gradient, reaction_ok)
        if (.not. reaction_ok) then
          ok = .false.
          return
        end if
        progress_derivative = progress_derivative - &
          rate_derivative * product_product / equilibrium_constant - &
          reverse_constant * product_gradient
      end if

      do i = 1, nspecies
        net_stoich = reactions(reaction_index)%product_stoich(i) - &
          reactions(reaction_index)%reactant_stoich(i)
        if (abs(net_stoich) <= tiny(1.0_dp)) cycle
        do j = 1, nspecies
          jacobian(i, j) = jacobian(i, j) + &
            net_stoich * progress_derivative(j)
        end do
      end do
    end do
    ok = all(ieee_is_finite(jacobian))
  end subroutine elementary_production_jacobian

  subroutine elementary_mass_fraction_rhs( &
      species, reactions, temperature, density, mass_fractions, derivative, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: temperature, density, mass_fractions(:)
    real(dp), intent(out) :: derivative(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: molar_production_rates(:)
    integer :: i

    derivative = 0.0_dp
    ok = size(derivative) == size(species) .and. density > 0.0_dp
    if (.not. ok) return
    allocate(molar_production_rates(size(species)))
    call elementary_production_rates( &
      species, reactions, temperature, density, mass_fractions, &
      molar_production_rates, ok)
    if (.not. ok) return
    do i = 1, size(species)
      derivative(i) = species(i)%molecular_weight * &
        molar_production_rates(i) / density
    end do
    ok = all(ieee_is_finite(derivative))
  end subroutine elementary_mass_fraction_rhs

  subroutine elementary_mass_fraction_jacobian( &
      species, reactions, temperature, density, mass_fractions, jacobian, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: temperature, density, mass_fractions(:)
    real(dp), intent(out) :: jacobian(:, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: concentrations(:), molar_jacobian(:, :)
    integer :: i, j, nspecies

    jacobian = 0.0_dp
    nspecies = size(species)
    ok = density > 0.0_dp .and. size(jacobian, 1) == nspecies .and. &
      size(jacobian, 2) == nspecies
    if (.not. ok) return
    allocate(concentrations(nspecies), molar_jacobian(nspecies, nspecies))
    call molar_concentrations_from_mass_fractions( &
      species, density, mass_fractions, concentrations, ok)
    if (.not. ok) return
    call elementary_production_jacobian( &
      species, reactions, temperature, concentrations, molar_jacobian, ok)
    if (.not. ok) return
    do i = 1, nspecies
      do j = 1, nspecies
        jacobian(i, j) = species(i)%molecular_weight / &
          species(j)%molecular_weight * molar_jacobian(i, j)
      end do
    end do
    ok = all(ieee_is_finite(jacobian))
  end subroutine elementary_mass_fraction_jacobian

end module elementary_kinetics_mod
