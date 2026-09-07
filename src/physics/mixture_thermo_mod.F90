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
  real(dp), parameter :: finite_guard_limit = huge(1.0_dp) / 2.0_dp

  public :: valid_mixture_composition
  public :: mixture_temperature_bounds
  public :: mixture_mass_properties
  public :: mixture_molecular_weight
  public :: mixture_specific_gas_constant
  public :: mixture_pressure
  public :: mixture_density
  public :: mixture_sound_speed
  public :: temperature_from_internal_energy
  public :: mass_fractions_from_mole_fractions
  public :: mole_fractions_from_mass_fractions

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
    ! Every term is finite and bounded near [0, 1], so this reduction and the
    ! subtraction cannot overflow for any representable array extent.
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
    if (.not. all(ieee_is_finite([lower, upper]))) then
      lower = 0.0_dp
      upper = 0.0_dp
      return
    end if
    do i = 2, size(species)
      if (.not. valid_nasa7_species(species(i))) then
        lower = 0.0_dp
        upper = 0.0_dp
        return
      end if
      lower = max(lower, species(i)%temperature_min)
      upper = min(upper, species(i)%temperature_max)
      if (.not. all(ieee_is_finite([lower, upper]))) then
        lower = 0.0_dp
        upper = 0.0_dp
        return
      end if
    end do
    if (upper <= lower) then
      lower = 0.0_dp
      upper = 0.0_dp
      return
    end if
    ok = .true.
  end subroutine mixture_temperature_bounds

  real(dp) function mixture_molecular_weight( &
      species, mass_fractions, ok) result(molecular_weight)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: mass_fractions(:)
    logical, intent(out) :: ok

    real(dp) :: inverse_weight
    real(dp) :: numerator, contribution
    integer :: i

    molecular_weight = 0.0_dp
    ok = .false.
    if (.not. valid_mixture_composition(species, mass_fractions)) return

    inverse_weight = 0.0_dp
    do i = 1, size(species)
      numerator = max(0.0_dp, mass_fractions(i))
      if (species(i)%molecular_weight < 1.0_dp) then
        if (numerator > &
            finite_guard_limit * species(i)%molecular_weight) return
      end if
      contribution = numerator / species(i)%molecular_weight
      if (contribution > finite_guard_limit - inverse_weight) return
      inverse_weight = inverse_weight + contribution
    end do
    if (inverse_weight <= 0.0_dp) return
    if (inverse_weight < 1.0_dp / finite_guard_limit) return
    molecular_weight = 1.0_dp / inverse_weight
    if (.not. ieee_is_finite(molecular_weight)) then
      molecular_weight = 0.0_dp
      return
    end if
    if (molecular_weight <= 0.0_dp) then
      molecular_weight = 0.0_dp
      return
    end if
    ok = .true.
  end function mixture_molecular_weight

  real(dp) function mixture_specific_gas_constant( &
      species, mass_fractions, ok) result(gas_constant)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: mass_fractions(:)
    logical, intent(out) :: ok

    real(dp) :: molecular_weight

    gas_constant = 0.0_dp
    molecular_weight = mixture_molecular_weight(species, mass_fractions, ok)
    if (.not. ok) return
    if (molecular_weight < &
        universal_gas_constant / finite_guard_limit) then
      ok = .false.
      return
    end if
    gas_constant = universal_gas_constant / molecular_weight
    if (.not. ieee_is_finite(gas_constant)) then
      gas_constant = 0.0_dp
      ok = .false.
      return
    end if
    if (gas_constant <= 0.0_dp) then
      gas_constant = 0.0_dp
      ok = .false.
      return
    end if
    ok = .true.
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
    real(dp) :: candidate_molecular_weight, candidate_gas_constant
    real(dp) :: candidate_cp, candidate_cv, candidate_gamma
    real(dp) :: candidate_enthalpy, candidate_internal_energy, candidate_entropy
    real(dp) :: weight, product, overflow_limit
    real(dp) :: species_values(5), candidate_values(5)
    logical :: species_ok, local_ok
    integer :: i, property

    molecular_weight = 0.0_dp
    gas_constant = 0.0_dp
    cp = 0.0_dp
    cv = 0.0_dp
    gamma = 0.0_dp
    enthalpy = 0.0_dp
    internal_energy = 0.0_dp
    entropy = 0.0_dp
    ok = .false.
    if (.not. ieee_is_finite(temperature)) return
    if (temperature <= 0.0_dp) return

    candidate_molecular_weight = &
      mixture_molecular_weight(species, mass_fractions, local_ok)
    if (.not. local_ok) return
    if (candidate_molecular_weight < &
        universal_gas_constant / finite_guard_limit) return
    candidate_gas_constant = &
      universal_gas_constant / candidate_molecular_weight
    if (.not. ieee_is_finite(candidate_gas_constant)) return
    if (candidate_gas_constant <= 0.0_dp) return

    candidate_values = 0.0_dp

    do i = 1, size(species)
      call nasa7_mass_properties( &
        species(i), temperature, species_cp, species_cv, species_h, &
        species_u, species_s, species_ok)
      if (.not. species_ok) return
      if (.not. all(ieee_is_finite([ &
          species_cp, species_cv, species_h, species_u, species_s]))) return
      weight = max(0.0_dp, mass_fractions(i))
      species_values = [species_cp, species_cv, species_h, species_u, species_s]
      if (weight > 1.0_dp) then
        overflow_limit = finite_guard_limit / weight
        if (any(abs(species_values) > overflow_limit)) return
      end if
      do property = 1, size(species_values)
        product = weight * species_values(property)
        if (product > 0.0_dp) then
          if (candidate_values(property) > &
              finite_guard_limit - product) return
        else if (product < 0.0_dp) then
          if (candidate_values(property) < &
              -finite_guard_limit - product) return
        end if
        candidate_values(property) = candidate_values(property) + product
      end do
    end do

    candidate_cp = candidate_values(1)
    candidate_cv = candidate_values(2)
    candidate_enthalpy = candidate_values(3)
    candidate_internal_energy = candidate_values(4)
    candidate_entropy = candidate_values(5)

    if (.not. all(ieee_is_finite([ &
        candidate_cp, candidate_cv, candidate_enthalpy, &
        candidate_internal_energy, candidate_entropy]))) return
    if (candidate_cv <= 0.0_dp .or. candidate_cp <= candidate_cv) return
    if (candidate_cv < candidate_cp / finite_guard_limit) return
    candidate_gamma = candidate_cp / candidate_cv
    if (.not. ieee_is_finite(candidate_gamma)) return
    if (candidate_gamma <= 0.0_dp) return
    if (.not. all(ieee_is_finite([ &
        candidate_molecular_weight, candidate_gas_constant, candidate_cp, &
        candidate_cv, candidate_gamma, candidate_enthalpy, &
        candidate_internal_energy, candidate_entropy]))) return

    molecular_weight = candidate_molecular_weight
    gas_constant = candidate_gas_constant
    cp = candidate_cp
    cv = candidate_cv
    gamma = candidate_gamma
    enthalpy = candidate_enthalpy
    internal_energy = candidate_internal_energy
    entropy = candidate_entropy
    ok = .true.
  end subroutine mixture_mass_properties

  real(dp) function mixture_pressure( &
      species, mass_fractions, density, temperature, ok) result(pressure)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: mass_fractions(:), density, temperature
    logical, intent(out) :: ok

    real(dp) :: gas_constant, density_gas_constant
    logical :: local_ok

    pressure = 0.0_dp
    ok = .false.
    if (.not. all(ieee_is_finite([density, temperature]))) return
    if (density <= 0.0_dp .or. temperature <= 0.0_dp) return
    gas_constant = mixture_specific_gas_constant(species, mass_fractions, ok)
    if (.not. ok) return
    call multiply_finite( &
      density, gas_constant, density_gas_constant, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call multiply_finite( &
      density_gas_constant, temperature, pressure, local_ok)
    if (.not. local_ok .or. pressure <= 0.0_dp) then
      pressure = 0.0_dp
      ok = .false.
      return
    end if
    ok = .true.
  end function mixture_pressure

  real(dp) function mixture_density( &
      species, mass_fractions, pressure, temperature, ok) result(density)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: mass_fractions(:), pressure, temperature
    logical, intent(out) :: ok

    real(dp) :: gas_constant, denominator
    logical :: local_ok

    density = 0.0_dp
    ok = .false.
    if (.not. all(ieee_is_finite([pressure, temperature]))) return
    if (pressure <= 0.0_dp .or. temperature <= 0.0_dp) return
    gas_constant = mixture_specific_gas_constant(species, mass_fractions, ok)
    if (.not. ok) return
    call multiply_finite(gas_constant, temperature, denominator, local_ok)
    if (.not. local_ok .or. denominator <= 0.0_dp) then
      density = 0.0_dp
      ok = .false.
      return
    end if
    call divide_finite(pressure, denominator, density, local_ok)
    if (.not. local_ok .or. density <= 0.0_dp) then
      density = 0.0_dp
      ok = .false.
      return
    end if
    ok = .true.
  end function mixture_density

  real(dp) function mixture_sound_speed( &
      species, mass_fractions, temperature, ok) result(sound_speed)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: mass_fractions(:), temperature
    logical, intent(out) :: ok

    real(dp) :: molecular_weight, gas_constant, cp, cv, gamma
    real(dp) :: enthalpy, internal_energy, entropy, gamma_gas_constant
    real(dp) :: sound_speed_squared
    logical :: local_ok

    sound_speed = 0.0_dp
    ok = .false.
    call mixture_mass_properties( &
      species, mass_fractions, temperature, molecular_weight, gas_constant, &
      cp, cv, gamma, enthalpy, internal_energy, entropy, ok)
    if (.not. ok) return
    if (.not. all(ieee_is_finite([gamma, gas_constant, temperature]))) then
      ok = .false.
      return
    end if
    if (gamma <= 0.0_dp .or. gas_constant <= 0.0_dp .or. &
        temperature <= 0.0_dp) then
      ok = .false.
      return
    end if
    call multiply_finite(gamma, gas_constant, gamma_gas_constant, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call multiply_finite( &
      gamma_gas_constant, temperature, sound_speed_squared, local_ok)
    if (.not. local_ok .or. sound_speed_squared <= 0.0_dp) then
      sound_speed = 0.0_dp
      ok = .false.
      return
    end if
    sound_speed = sqrt(sound_speed_squared)
    if (.not. ieee_is_finite(sound_speed)) then
      sound_speed = 0.0_dp
      ok = .false.
      return
    end if
    if (sound_speed <= 0.0_dp) then
      sound_speed = 0.0_dp
      ok = .false.
      return
    end if
    ok = .true.
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
    real(dp) :: residual, candidate, tolerance, step
    real(dp) :: lower_limit, upper_limit
    real(dp) :: working_temperature
    logical :: properties_ok, local_ok, candidate_ok
    integer :: iteration

    temperature = 0.0_dp
    ok = .false.
    if (present(iterations)) iterations = 0
    if (.not. valid_mixture_composition(species, mass_fractions)) return
    if (.not. all(ieee_is_finite([target_internal_energy, initial_guess]))) return

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
    if (.not. all(ieee_is_finite([energy_lower, energy_upper]))) return

    tolerance = 1.0e-11_dp * max(1.0_dp, abs(target_internal_energy))
    if (.not. ieee_is_finite(tolerance)) return
    if (tolerance <= 0.0_dp) return
    call subtract_finite(energy_lower, tolerance, lower_limit, local_ok)
    if (.not. local_ok) return
    upper_limit = energy_upper
    call add_finite(upper_limit, tolerance, local_ok)
    if (.not. local_ok) return
    if (target_internal_energy < lower_limit .or. &
        target_internal_energy > upper_limit) return

    working_temperature = min(upper, max(lower, initial_guess))
    if (.not. ieee_is_finite(working_temperature)) return
    do iteration = 1, temperature_inversion_max_iterations
      call mixture_mass_properties( &
        species, mass_fractions, working_temperature, molecular_weight, gas_constant, &
        cp, cv, gamma, enthalpy, internal_energy, entropy, properties_ok)
      if (.not. properties_ok) return

      call subtract_finite( &
        internal_energy, target_internal_energy, residual, local_ok)
      if (.not. local_ok) return
      if (abs(residual) <= tolerance) then
        temperature = working_temperature
        ok = .true.
        if (present(iterations)) iterations = iteration
        return
      end if

      if (residual > 0.0_dp) then
        upper = working_temperature
      else
        lower = working_temperature
      end if

      call divide_finite(residual, cv, step, local_ok)
      if (.not. local_ok) return
      call subtract_finite(working_temperature, step, candidate, candidate_ok)
      if (.not. candidate_ok) then
        call midpoint_finite(lower, upper, candidate, candidate_ok)
      else
        if (.not. ieee_is_finite(candidate)) return
        if (candidate <= lower .or. candidate >= upper) then
          call midpoint_finite(lower, upper, candidate, candidate_ok)
        end if
      end if
      if (.not. candidate_ok .or. .not. ieee_is_finite(candidate)) return
      working_temperature = candidate
    end do

    if (present(iterations)) iterations = &
      temperature_inversion_max_iterations
  end subroutine temperature_from_internal_energy

  subroutine mass_fractions_from_mole_fractions( &
      species, mole_fractions, mass_fractions, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: mole_fractions(:)
    real(dp), intent(out) :: mass_fractions(:)
    logical, intent(out) :: ok

    real(dp) :: denominator, numerator, total, sum_difference
    logical :: local_ok
    integer :: i

    mass_fractions = 0.0_dp
    ok = .false.
    if (size(species) < 1 .or. size(mole_fractions) /= size(species) .or. &
        size(mass_fractions) /= size(species)) return
    if (any(.not. ieee_is_finite(mole_fractions))) then
      return
    end if
    if (any(mole_fractions < -mixture_composition_tolerance)) return
    if (any(mole_fractions > 1.0_dp + mixture_composition_tolerance)) return
    call finite_sum(mole_fractions, total, local_ok)
    if (.not. local_ok) return
    call subtract_finite(total, 1.0_dp, sum_difference, local_ok)
    if (.not. local_ok .or. abs(sum_difference) > &
        mixture_composition_tolerance) return

    denominator = 0.0_dp
    do i = 1, size(species)
      if (.not. valid_nasa7_species(species(i))) return
      call multiply_finite( &
        max(0.0_dp, mole_fractions(i)), species(i)%molecular_weight, &
        numerator, local_ok)
      if (.not. local_ok) return
      call add_finite(denominator, numerator, local_ok)
      if (.not. local_ok) return
    end do
    if (denominator <= 0.0_dp) return
    do i = 1, size(species)
      call multiply_finite( &
        max(0.0_dp, mole_fractions(i)), species(i)%molecular_weight, &
        numerator, local_ok)
      if (.not. local_ok) then
        mass_fractions = 0.0_dp
        return
      end if
      call divide_finite(numerator, denominator, mass_fractions(i), local_ok)
      if (.not. local_ok .or. mass_fractions(i) < 0.0_dp) then
        mass_fractions = 0.0_dp
        return
      end if
    end do
    ok = valid_mixture_composition(species, mass_fractions)
    if (.not. ok) mass_fractions = 0.0_dp
  end subroutine mass_fractions_from_mole_fractions

  subroutine mole_fractions_from_mass_fractions( &
      species, mass_fractions, mole_fractions, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: mass_fractions(:)
    real(dp), intent(out) :: mole_fractions(:)
    logical, intent(out) :: ok

    real(dp) :: denominator, numerator, total, sum_difference
    logical :: local_ok
    integer :: i

    mole_fractions = 0.0_dp
    ok = .false.
    if (size(mole_fractions) /= size(species)) return
    if (.not. valid_mixture_composition(species, mass_fractions)) return
    denominator = 0.0_dp
    do i = 1, size(species)
      call divide_finite( &
        max(0.0_dp, mass_fractions(i)), species(i)%molecular_weight, &
        numerator, local_ok)
      if (.not. local_ok) return
      call add_finite(denominator, numerator, local_ok)
      if (.not. local_ok) return
    end do
    if (denominator <= 0.0_dp) return
    do i = 1, size(species)
      call divide_finite( &
        max(0.0_dp, mass_fractions(i)), species(i)%molecular_weight, &
        numerator, local_ok)
      if (.not. local_ok) then
        mole_fractions = 0.0_dp
        return
      end if
      call divide_finite(numerator, denominator, mole_fractions(i), local_ok)
      if (.not. local_ok .or. mole_fractions(i) < 0.0_dp) then
        mole_fractions = 0.0_dp
        return
      end if
    end do
    if (.not. all(ieee_is_finite(mole_fractions))) then
      mole_fractions = 0.0_dp
      return
    end if
    call finite_sum(mole_fractions, total, local_ok)
    if (.not. local_ok) then
      mole_fractions = 0.0_dp
      return
    end if
    call subtract_finite(total, 1.0_dp, sum_difference, local_ok)
    if (.not. local_ok .or. abs(sum_difference) > &
        mixture_composition_tolerance) then
      mole_fractions = 0.0_dp
      return
    end if
    ok = .true.
  end subroutine mole_fractions_from_mass_fractions

  subroutine finite_sum(values, total, ok)
    real(dp), intent(in) :: values(:)
    real(dp), intent(out) :: total
    logical, intent(out) :: ok

    integer :: i

    total = 0.0_dp
    ok = .true.
    do i = 1, size(values)
      call add_finite(total, values(i), ok)
      if (.not. ok) then
        total = 0.0_dp
        return
      end if
    end do
  end subroutine finite_sum

  subroutine add_finite(accumulator, increment, ok)
    real(dp), intent(inout) :: accumulator
    real(dp), intent(in) :: increment
    logical, intent(out) :: ok

    ok = .false.
    if (.not. all(ieee_is_finite([accumulator, increment]))) return
    if (increment > 0.0_dp) then
      if (accumulator > finite_guard_limit - increment) return
    else if (increment < 0.0_dp) then
      if (accumulator < -finite_guard_limit - increment) return
    end if
    accumulator = accumulator + increment
    ok = ieee_is_finite(accumulator)
  end subroutine add_finite

  subroutine multiply_finite(left, right, product, ok)
    real(dp), intent(in) :: left, right
    real(dp), intent(out) :: product
    logical, intent(out) :: ok

    product = 0.0_dp
    ok = .false.
    if (.not. all(ieee_is_finite([left, right]))) return
    if ((left <= 0.0_dp .and. left >= 0.0_dp) .or. &
        (right <= 0.0_dp .and. right >= 0.0_dp)) then
      ok = .true.
      return
    end if
    if (abs(left) > 1.0_dp .and. abs(right) > 1.0_dp) then
      if (abs(left) > finite_guard_limit / abs(right)) return
    end if
    product = left * right
    ok = ieee_is_finite(product)
    if (.not. ok) product = 0.0_dp
  end subroutine multiply_finite

  subroutine divide_finite(numerator, denominator, quotient, ok)
    real(dp), intent(in) :: numerator, denominator
    real(dp), intent(out) :: quotient
    logical, intent(out) :: ok

    quotient = 0.0_dp
    ok = .false.
    if (.not. all(ieee_is_finite([numerator, denominator]))) return
    if (denominator <= 0.0_dp .and. denominator >= 0.0_dp) return
    if (numerator <= 0.0_dp .and. numerator >= 0.0_dp) then
      ok = .true.
      return
    end if
    if (abs(denominator) < abs(numerator) / finite_guard_limit) return
    quotient = numerator / denominator
    ok = ieee_is_finite(quotient)
    if (.not. ok) quotient = 0.0_dp
  end subroutine divide_finite

  subroutine subtract_finite(left, right, difference, ok)
    real(dp), intent(in) :: left, right
    real(dp), intent(out) :: difference
    logical, intent(out) :: ok

    difference = 0.0_dp
    ok = .false.
    if (.not. all(ieee_is_finite([left, right]))) return
    difference = left
    call add_finite(difference, -right, ok)
    if (.not. ok) difference = 0.0_dp
  end subroutine subtract_finite

  subroutine midpoint_finite(lower, upper, midpoint, ok)
    real(dp), intent(in) :: lower, upper
    real(dp), intent(out) :: midpoint
    logical, intent(out) :: ok

    real(dp) :: sum

    midpoint = 0.0_dp
    ok = .false.
    if (.not. all(ieee_is_finite([lower, upper]))) return
    sum = lower
    call add_finite(sum, upper, ok)
    if (ok) then
      midpoint = 0.5_dp * sum
    else
      midpoint = 0.5_dp * lower + 0.5_dp * upper
      ok = ieee_is_finite(midpoint)
    end if
    if (.not. ok) midpoint = 0.0_dp
  end subroutine midpoint_finite

end module mixture_thermo_mod
