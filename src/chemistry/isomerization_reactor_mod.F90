module isomerization_reactor_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use mixture_thermo_mod, only: &
    valid_mixture_composition, mixture_mass_properties, &
    temperature_from_internal_energy
  implicit none
  private

  type, public :: isomerization_reaction
    integer :: reactant = 1
    integer :: product = 2
    real(dp) :: pre_exponential = 0.0_dp
    real(dp) :: temperature_exponent = 0.0_dp
    real(dp) :: activation_temperature = 0.0_dp ! Ea/R in kelvin
  end type isomerization_reaction

  public :: valid_isomerization_reaction
  public :: isomerization_rate_constant
  public :: advance_isomerization_rk4
  public :: reactor_internal_energy

contains

  logical function valid_isomerization_reaction( &
      reaction, nspecies) result(valid)
    type(isomerization_reaction), intent(in) :: reaction
    integer, intent(in) :: nspecies

    valid = nspecies >= 2 .and. reaction%reactant >= 1 .and. &
      reaction%reactant <= nspecies .and. reaction%product >= 1 .and. &
      reaction%product <= nspecies .and. &
      reaction%reactant /= reaction%product .and. &
      reaction%pre_exponential >= 0.0_dp .and. &
      reaction%activation_temperature >= 0.0_dp
  end function valid_isomerization_reaction

  real(dp) function isomerization_rate_constant( &
      reaction, temperature, ok) result(rate_constant)
    type(isomerization_reaction), intent(in) :: reaction
    real(dp), intent(in) :: temperature
    logical, intent(out) :: ok

    rate_constant = 0.0_dp
    ok = temperature > 0.0_dp .and. reaction%pre_exponential >= 0.0_dp .and. &
      reaction%activation_temperature >= 0.0_dp
    if (.not. ok) return

    rate_constant = reaction%pre_exponential * &
      temperature**reaction%temperature_exponent * &
      exp(-reaction%activation_temperature / temperature)
    ok = ieee_is_finite(rate_constant) .and. rate_constant >= 0.0_dp
  end function isomerization_rate_constant

  real(dp) function reactor_internal_energy( &
      species, mass_fractions, temperature, ok) result(internal_energy)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: mass_fractions(:), temperature
    logical, intent(out) :: ok

    real(dp) :: molecular_weight, gas_constant, cp, cv, gamma
    real(dp) :: enthalpy, entropy

    internal_energy = 0.0_dp
    call mixture_mass_properties( &
      species, mass_fractions, temperature, molecular_weight, gas_constant, &
      cp, cv, gamma, enthalpy, internal_energy, entropy, ok)
  end function reactor_internal_energy

  subroutine advance_isomerization_rk4( &
      species, reaction, time_step, adiabatic, target_internal_energy, &
      mass_fractions, temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(isomerization_reaction), intent(in) :: reaction
    real(dp), intent(in) :: time_step, target_internal_energy
    logical, intent(in) :: adiabatic
    real(dp), intent(inout) :: mass_fractions(:), temperature
    logical, intent(out) :: ok

    real(dp), allocatable :: k1(:), k2(:), k3(:), k4(:), stage(:), updated(:)
    real(dp) :: fixed_temperature
    real(dp) :: stage_temperature_1, stage_temperature_2
    real(dp) :: stage_temperature_3, stage_temperature_4
    logical :: stage_ok

    ok = .false.
    if (time_step <= 0.0_dp) return
    if (size(species) /= size(mass_fractions)) return
    if (.not. valid_mixture_composition(species, mass_fractions)) return
    if (.not. valid_isomerization_reaction(reaction, size(species))) return

    allocate(k1(size(species)), k2(size(species)), k3(size(species)))
    allocate(k4(size(species)), stage(size(species)), updated(size(species)))
    fixed_temperature = temperature

    call evaluate_composition_rhs( &
      mass_fractions, temperature, k1, stage_temperature_1, stage_ok)
    if (.not. stage_ok) return

    stage = mass_fractions + 0.5_dp * time_step * k1
    call evaluate_composition_rhs( &
      stage, stage_temperature_1, k2, stage_temperature_2, stage_ok)
    if (.not. stage_ok) return

    stage = mass_fractions + 0.5_dp * time_step * k2
    call evaluate_composition_rhs( &
      stage, stage_temperature_2, k3, stage_temperature_3, stage_ok)
    if (.not. stage_ok) return

    stage = mass_fractions + time_step * k3
    call evaluate_composition_rhs( &
      stage, stage_temperature_3, k4, stage_temperature_4, stage_ok)
    if (.not. stage_ok) return

    updated = mass_fractions + time_step * &
      (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4) / 6.0_dp
    call enforce_roundoff_composition(updated, stage_ok)
    if (.not. stage_ok) return

    if (adiabatic) then
      call temperature_from_internal_energy( &
        species, updated, target_internal_energy, stage_temperature_4, &
        temperature, stage_ok)
      if (.not. stage_ok) return
    else
      temperature = fixed_temperature
    end if

    mass_fractions = updated
    ok = .true.

  contains

    subroutine evaluate_composition_rhs( &
        composition, temperature_guess, derivative, evaluated_temperature, &
        rhs_ok)
      real(dp), intent(in) :: composition(:), temperature_guess
      real(dp), intent(out) :: derivative(:), evaluated_temperature
      logical, intent(out) :: rhs_ok

      real(dp) :: rate_constant, progress_rate

      derivative = 0.0_dp
      rhs_ok = valid_mixture_composition(species, composition)
      if (.not. rhs_ok) return

      if (adiabatic) then
        call temperature_from_internal_energy( &
          species, composition, target_internal_energy, temperature_guess, &
          evaluated_temperature, rhs_ok)
        if (.not. rhs_ok) return
      else
        evaluated_temperature = fixed_temperature
      end if

      rate_constant = isomerization_rate_constant( &
        reaction, evaluated_temperature, rhs_ok)
      if (.not. rhs_ok) return
      progress_rate = rate_constant * composition(reaction%reactant)
      derivative(reaction%reactant) = -progress_rate
      derivative(reaction%product) = progress_rate
    end subroutine evaluate_composition_rhs

  end subroutine advance_isomerization_rk4

  subroutine enforce_roundoff_composition(mass_fractions, ok)
    real(dp), intent(inout) :: mass_fractions(:)
    logical, intent(out) :: ok

    real(dp), parameter :: tolerance = 1.0e-11_dp
    real(dp) :: total

    ok = .false.
    if (any(mass_fractions < -tolerance)) return
    mass_fractions = max(0.0_dp, mass_fractions)
    total = sum(mass_fractions)
    if (abs(total - 1.0_dp) > tolerance .or. total <= 0.0_dp) return
    mass_fractions = mass_fractions / total
    ok = .true.
  end subroutine enforce_roundoff_composition

end module isomerization_reactor_mod
