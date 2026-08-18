module mixture_thermo_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: &
    nasa7_species, valid_nasa7_species, nasa7_mass_properties, &
    universal_gas_constant
  implicit none
  private

  real(dp), parameter, public :: mixture_composition_tolerance = 5.0e-12_dp
  integer, parameter, public :: temperature_inversion_max_iterations = 100

  public :: valid_mixture_composition
  public :: mixture_temperature_bounds
  public :: mixture_mass_properties
  public :: mixture_molecular_weight
  public :: mixture_specific_gas_constant
  public :: mixture_pressure
  public :: mixture_density
  public :: mixture_sound_speed
  public :: temperature_from_internal_energy

contains

  logical function valid_mixture_composition(species, mass_fractions) &
      result(valid)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: mass_fractions(:)

    integer :: i

    valid = .false.
    if (size(species) < 1 .or. size(mass_fractions) /= size(species)) return
    do i = 1, size(species)
      if (.not. valid_nasa7_species(species(i))) return
    end do
    if (any(.not. ieee_is_finite(mass_fractions))) return
    if (any(mass_fractions < -mixture_composition_tolerance)) return
    if (any(mass_fractions > 1.0_dp + mixture_composition_tolerance)) return
    if (abs(sum(mass_fractions) - 1.0_dp) > &
        mixture_composition_tolerance) return
    valid = .true.
  end function valid_mixture_composition

  subroutine mixture_temperature_bounds(species, lower, upper, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(out) :: lower, upper
    logical, intent(out) :: ok

    integer :: i

    lower = 0.0_dp
    upper = 0.0_dp
    ok = .false.
    if (size(species) < 1) return
    if (.not. valid_nasa7_species(species(1))) return

    lower = species(1)%temperature_min
    upper = species(1)%temperature_max
    do i = 2, size(species)
      if (.not. valid_nasa7_species(species(i))) return
      lower = max(lower, species(i)%temperature_min)
      upper = min(upper, species(i)%temperature_max)
    end do
    ok = upper > lower
  end subroutine mixture_temperature_bounds

  real(dp) function mixture_molecular_weight( &
      species, mass_fractions, ok) result(molecular_weight)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: mass_fractions(:)
    logical, intent(out) :: ok

    real(dp) :: inverse_weight
    integer :: i

    molecular_weight = 0.0_dp
    ok = valid_mixture_composition(species, mass_fractions)
    if (.not. ok) return

    inverse_weight = 0.0_dp
    do i = 1, size(species)
      inverse_weight = inverse_weight + &
        max(0.0_dp, mass_fractions(i)) / species(i)%molecular_weight
    end do
    if (inverse_weight <= 0.0_dp) then
      ok = .false.
      return
    end if
    molecular_weight = 1.0_dp / inverse_weight
  end function mixture_molecular_weight

  real(dp) function mixture_specific_gas_constant( &
      species, mass_fractions, ok) result(gas_constant)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: mass_fractions(:)
    logical, intent(out) :: ok

    real(dp) :: molecular_weight

    molecular_weight = mixture_molecular_weight(species, mass_fractions, ok)
    if (.not. ok) then
      gas_constant = 0.0_dp
    else
      gas_constant = universal_gas_constant / molecular_weight
    end if
  end function mixture_specific_gas_constant

  subroutine mixture_mass_properties( &
      species, mass_fractions, temperature, molecular_weight, gas_constant, &
      cp, cv, gamma, enthalpy, internal_energy, entropy, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: mass_fractions(:), temperature
    real(dp), intent(out) :: molecular_weight, gas_constant
    real(dp), intent(out) :: cp, cv, gamma
    real(dp), intent(out) :: enthalpy, internal_energy, entropy
    logical, intent(out) :: ok

    real(dp) :: species_cp, species_cv, species_h, species_u, species_s
    logical :: species_ok
    integer :: i

    molecular_weight = 0.0_dp
    gas_constant = 0.0_dp
    cp = 0.0_dp
    cv = 0.0_dp
    gamma = 0.0_dp
    enthalpy = 0.0_dp
    internal_energy = 0.0_dp
    entropy = 0.0_dp

    molecular_weight = mixture_molecular_weight(species, mass_fractions, ok)
    if (.not. ok) return
    gas_constant = universal_gas_constant / molecular_weight

    do i = 1, size(species)
      call nasa7_mass_properties( &
        species(i), temperature, species_cp, species_cv, species_h, &
        species_u, species_s, species_ok)
      if (.not. species_ok) then
        ok = .false.
        return
      end if
      cp = cp + max(0.0_dp, mass_fractions(i)) * species_cp
      cv = cv + max(0.0_dp, mass_fractions(i)) * species_cv
      enthalpy = enthalpy + max(0.0_dp, mass_fractions(i)) * species_h
      internal_energy = internal_energy + &
        max(0.0_dp, mass_fractions(i)) * species_u
      entropy = entropy + max(0.0_dp, mass_fractions(i)) * species_s
    end do

    if (cv <= 0.0_dp .or. cp <= cv) then
      ok = .false.
      return
    end if
    gamma = cp / cv
    ok = all(ieee_is_finite([molecular_weight, gas_constant, cp, cv, gamma, &
      enthalpy, internal_energy, entropy]))
  end subroutine mixture_mass_properties

  real(dp) function mixture_pressure( &
      species, mass_fractions, density, temperature, ok) result(pressure)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: mass_fractions(:), density, temperature
    logical, intent(out) :: ok

    real(dp) :: gas_constant

    pressure = 0.0_dp
    if (density <= 0.0_dp .or. temperature <= 0.0_dp) then
      ok = .false.
      return
    end if
    gas_constant = mixture_specific_gas_constant(species, mass_fractions, ok)
    if (.not. ok) return
    pressure = density * gas_constant * temperature
    ok = ieee_is_finite(pressure) .and. pressure > 0.0_dp
  end function mixture_pressure

  real(dp) function mixture_density( &
      species, mass_fractions, pressure, temperature, ok) result(density)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: mass_fractions(:), pressure, temperature
    logical, intent(out) :: ok

    real(dp) :: gas_constant

    density = 0.0_dp
    if (pressure <= 0.0_dp .or. temperature <= 0.0_dp) then
      ok = .false.
      return
    end if
    gas_constant = mixture_specific_gas_constant(species, mass_fractions, ok)
    if (.not. ok) return
    density = pressure / (gas_constant * temperature)
    ok = ieee_is_finite(density) .and. density > 0.0_dp
  end function mixture_density

  real(dp) function mixture_sound_speed( &
      species, mass_fractions, temperature, ok) result(sound_speed)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: mass_fractions(:), temperature
    logical, intent(out) :: ok

    real(dp) :: molecular_weight, gas_constant, cp, cv, gamma
    real(dp) :: enthalpy, internal_energy, entropy

    sound_speed = 0.0_dp
    call mixture_mass_properties( &
      species, mass_fractions, temperature, molecular_weight, gas_constant, &
      cp, cv, gamma, enthalpy, internal_energy, entropy, ok)
    if (.not. ok) return
    sound_speed = sqrt(gamma * gas_constant * temperature)
    ok = ieee_is_finite(sound_speed) .and. sound_speed > 0.0_dp
  end function mixture_sound_speed

  subroutine temperature_from_internal_energy( &
      species, mass_fractions, target_internal_energy, initial_guess, &
      temperature, ok, iterations)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: mass_fractions(:)
    real(dp), intent(in) :: target_internal_energy, initial_guess
    real(dp), intent(out) :: temperature
    logical, intent(out) :: ok
    integer, intent(out), optional :: iterations

    real(dp) :: lower, upper, energy_lower, energy_upper
    real(dp) :: molecular_weight, gas_constant, cp, cv, gamma
    real(dp) :: enthalpy, internal_energy, entropy
    real(dp) :: residual, candidate, tolerance
    logical :: properties_ok
    integer :: iteration

    temperature = 0.0_dp
    ok = .false.
    if (present(iterations)) iterations = 0
    if (.not. valid_mixture_composition(species, mass_fractions)) return
    if (.not. ieee_is_finite(target_internal_energy)) return

    call mixture_temperature_bounds(species, lower, upper, properties_ok)
    if (.not. properties_ok) return

    call mixture_mass_properties( &
      species, mass_fractions, lower, molecular_weight, gas_constant, cp, cv, &
      gamma, enthalpy, energy_lower, entropy, properties_ok)
    if (.not. properties_ok) return
    call mixture_mass_properties( &
      species, mass_fractions, upper, molecular_weight, gas_constant, cp, cv, &
      gamma, enthalpy, energy_upper, entropy, properties_ok)
    if (.not. properties_ok) return

    tolerance = 1.0e-11_dp * max(1.0_dp, abs(target_internal_energy))
    if (target_internal_energy < energy_lower - tolerance .or. &
        target_internal_energy > energy_upper + tolerance) return

    temperature = min(upper, max(lower, initial_guess))
    do iteration = 1, temperature_inversion_max_iterations
      call mixture_mass_properties( &
        species, mass_fractions, temperature, molecular_weight, gas_constant, &
        cp, cv, gamma, enthalpy, internal_energy, entropy, properties_ok)
      if (.not. properties_ok) return

      residual = internal_energy - target_internal_energy
      if (abs(residual) <= tolerance) then
        ok = .true.
        if (present(iterations)) iterations = iteration
        return
      end if

      if (residual > 0.0_dp) then
        upper = temperature
      else
        lower = temperature
      end if

      candidate = temperature - residual / cv
      if (.not. ieee_is_finite(candidate) .or. candidate <= lower .or. &
          candidate >= upper) candidate = 0.5_dp * (lower + upper)
      temperature = candidate
    end do

    if (present(iterations)) iterations = &
      temperature_inversion_max_iterations
  end subroutine temperature_from_internal_energy

end module mixture_thermo_mod
