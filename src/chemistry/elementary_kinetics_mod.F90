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

  type, public :: arrhenius_rate
    real(dp) :: pre_exponential = 0.0_dp
    real(dp) :: temperature_exponent = 0.0_dp
    real(dp) :: activation_energy = 0.0_dp ! J / kmol
  end type arrhenius_rate

  type, public :: elementary_reaction
    character(len=96) :: equation = ""
    real(dp), allocatable :: reactant_stoich(:)
    real(dp), allocatable :: product_stoich(:)
    type(arrhenius_rate) :: forward_rate
    logical :: reversible = .true.
  end type elementary_reaction

  public :: valid_arrhenius_rate
  public :: valid_elementary_reaction
  public :: arrhenius_rate_constant
  public :: reaction_equilibrium_constant
  public :: reaction_progress_rate
  public :: molar_concentrations_from_mass_fractions
  public :: elementary_production_rates
  public :: elementary_mass_fraction_rhs

contains

  logical function valid_arrhenius_rate(rate) result(valid)
    type(arrhenius_rate), intent(in) :: rate

    valid = rate%pre_exponential >= 0.0_dp .and. &
      all(ieee_is_finite([rate%pre_exponential, &
        rate%temperature_exponent, rate%activation_energy]))
  end function valid_arrhenius_rate

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
    if (.not. valid_arrhenius_rate(reaction%forward_rate)) return
    valid = .true.
  end function valid_elementary_reaction

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
    real(dp) :: reference_concentration
    real(dp) :: net_stoich
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

  subroutine reaction_progress_rate( &
      species, reaction, temperature, concentrations, forward_progress, &
      reverse_progress, net_progress, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reaction
    real(dp), intent(in) :: temperature, concentrations(:)
    real(dp), intent(out) :: forward_progress, reverse_progress, net_progress
    logical, intent(out) :: ok

    real(dp) :: forward_constant, reverse_constant, equilibrium_constant
    integer :: i

    forward_progress = 0.0_dp
    reverse_progress = 0.0_dp
    net_progress = 0.0_dp
    ok = valid_elementary_reaction(reaction, size(species)) .and. &
      size(concentrations) == size(species) .and. &
      all(ieee_is_finite(concentrations)) .and. &
      all(concentrations >= 0.0_dp)
    if (.not. ok) return

    forward_constant = arrhenius_rate_constant( &
      reaction%forward_rate, temperature, ok)
    if (.not. ok) return
    forward_progress = forward_constant
    do i = 1, size(species)
      if (reaction%reactant_stoich(i) > 0.0_dp) then
        forward_progress = forward_progress * &
          concentrations(i)**reaction%reactant_stoich(i)
      end if
    end do

    if (reaction%reversible) then
      equilibrium_constant = reaction_equilibrium_constant( &
        species, reaction, temperature, ok)
      if (.not. ok) return
      reverse_constant = forward_constant / equilibrium_constant
      reverse_progress = reverse_constant
      do i = 1, size(species)
        if (reaction%product_stoich(i) > 0.0_dp) then
          reverse_progress = reverse_progress * &
            concentrations(i)**reaction%product_stoich(i)
        end if
      end do
    end if

    net_progress = forward_progress - reverse_progress
    ok = all(ieee_is_finite( &
      [forward_progress, reverse_progress, net_progress]))
  end subroutine reaction_progress_rate

  subroutine molar_concentrations_from_mass_fractions( &
      species, density, mass_fractions, concentrations, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: density, mass_fractions(:)
    real(dp), intent(out) :: concentrations(:)
    logical, intent(out) :: ok

    integer :: i

    concentrations = 0.0_dp
    ok = density > 0.0_dp .and. size(concentrations) == size(species) .and. &
      valid_mixture_composition(species, mass_fractions)
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
    real(dp) :: q_forward, q_reverse, q_net
    logical :: reaction_ok
    integer :: i, reaction_index

    molar_production_rates = 0.0_dp
    ok = size(molar_production_rates) == size(species) .and. &
      size(reactions) >= 1
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

    allocate(concentrations(size(species)))
    call molar_concentrations_from_mass_fractions( &
      species, density, mass_fractions, concentrations, ok)
    if (.not. ok) return

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
  end subroutine elementary_production_rates

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

end module elementary_kinetics_mod
