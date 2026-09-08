module nasa7_thermo_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  implicit none
  private

  real(dp), parameter, public :: universal_gas_constant = &
    8.31446261815324e3_dp ! J / (kmol K)
  real(dp), parameter :: finite_guard_limit = huge(1.0_dp) / 2.0_dp

  type, public :: nasa7_species
    character(len=24) :: name = ""
    real(dp) :: molecular_weight = 0.0_dp ! kg / kmol
    real(dp) :: temperature_min = 0.0_dp
    real(dp) :: temperature_mid = 0.0_dp
    real(dp) :: temperature_max = 0.0_dp
    real(dp) :: low_coefficients(7) = 0.0_dp
    real(dp) :: high_coefficients(7) = 0.0_dp
  end type nasa7_species

  public :: valid_nasa7_species
  public :: nasa7_dimensionless_properties
  public :: nasa7_molar_properties
  public :: nasa7_mass_properties
  public :: nasa7_specific_gas_constant

contains

  logical function valid_nasa7_species(species) result(valid)
    type(nasa7_species), intent(in) :: species

    valid = .false.
    if (.not. all(ieee_is_finite([species%molecular_weight, &
        species%temperature_min, species%temperature_mid, &
        species%temperature_max]))) return
    if (.not. all(ieee_is_finite(species%low_coefficients))) return
    if (.not. all(ieee_is_finite(species%high_coefficients))) return
    if (species%molecular_weight <= 0.0_dp) return
    if (species%temperature_min <= 0.0_dp) return
    if (species%temperature_mid <= species%temperature_min) return
    if (species%temperature_max <= species%temperature_mid) return
    valid = .true.
  end function valid_nasa7_species

  real(dp) function nasa7_specific_gas_constant(species) result(gas_constant)
    type(nasa7_species), intent(in) :: species

    if (.not. valid_nasa7_species(species)) then
      gas_constant = -huge(1.0_dp)
    else if (species%molecular_weight < &
        universal_gas_constant / finite_guard_limit) then
      gas_constant = -huge(1.0_dp)
    else
      gas_constant = universal_gas_constant / species%molecular_weight
    end if
  end function nasa7_specific_gas_constant

  subroutine nasa7_dimensionless_properties( &
      species, temperature, cp_over_r, h_over_rt, s_over_r, ok)
    type(nasa7_species), intent(in) :: species
    real(dp), intent(in) :: temperature
    real(dp), intent(out) :: cp_over_r, h_over_rt, s_over_r
    logical, intent(out) :: ok

    real(dp) :: coefficients(7), temperature_squared
    real(dp) :: temperature_cubed, temperature_fourth, log_temperature
    real(dp) :: inverse_temperature, magnitude_scale, coefficient_limit

    cp_over_r = 0.0_dp
    h_over_rt = 0.0_dp
    s_over_r = 0.0_dp
    ok = .false.

    if (.not. valid_nasa7_species(species)) return
    if (.not. ieee_is_finite(temperature)) return
    if (temperature < species%temperature_min .or. &
        temperature > species%temperature_max) return

    if (temperature <= species%temperature_mid) then
      coefficients = species%low_coefficients
    else
      coefficients = species%high_coefficients
    end if

    if (temperature > 1.0_dp) then
      if (temperature > sqrt(sqrt(finite_guard_limit / 8.0_dp))) return
    end if
    if (temperature < 1.0_dp / finite_guard_limit) return
    temperature_squared = temperature * temperature
    temperature_cubed = temperature_squared * temperature
    temperature_fourth = temperature_squared * temperature_squared
    inverse_temperature = 1.0_dp / temperature
    log_temperature = log(temperature)
    if (.not. ieee_is_finite(log_temperature)) return
    magnitude_scale = max(1.0_dp, temperature, temperature_squared, &
      temperature_cubed, temperature_fourth, inverse_temperature, &
      abs(log_temperature))
    coefficient_limit = (finite_guard_limit / 8.0_dp) / magnitude_scale
    if (any(abs(coefficients) > coefficient_limit)) return

    cp_over_r = coefficients(1) + &
      coefficients(2) * temperature + &
      coefficients(3) * temperature_squared + &
      coefficients(4) * temperature_cubed + &
      coefficients(5) * temperature_fourth

    h_over_rt = coefficients(1) + &
      0.5_dp * coefficients(2) * temperature + &
      coefficients(3) * temperature_squared / 3.0_dp + &
      0.25_dp * coefficients(4) * temperature_cubed + &
      coefficients(5) * temperature_fourth / 5.0_dp + &
      coefficients(6) / temperature

    s_over_r = coefficients(1) * log_temperature + &
      coefficients(2) * temperature + &
      0.5_dp * coefficients(3) * temperature_squared + &
      coefficients(4) * temperature_cubed / 3.0_dp + &
      0.25_dp * coefficients(5) * temperature_fourth + &
      coefficients(7)

    if (.not. all(ieee_is_finite([cp_over_r, h_over_rt, s_over_r]))) then
      cp_over_r = 0.0_dp
      h_over_rt = 0.0_dp
      s_over_r = 0.0_dp
      return
    end if
    if (cp_over_r <= 1.0_dp) then
      cp_over_r = 0.0_dp
      h_over_rt = 0.0_dp
      s_over_r = 0.0_dp
      return
    end if
    ok = .true.
  end subroutine nasa7_dimensionless_properties

  subroutine nasa7_molar_properties( &
      species, temperature, cp, cv, enthalpy, internal_energy, entropy, ok)
    type(nasa7_species), intent(in) :: species
    real(dp), intent(in) :: temperature
    real(dp), intent(out) :: cp, cv, enthalpy, internal_energy, entropy
    logical, intent(out) :: ok

    real(dp) :: cp_over_r, h_over_rt, s_over_r, gas_temperature

    cp = 0.0_dp
    cv = 0.0_dp
    enthalpy = 0.0_dp
    internal_energy = 0.0_dp
    entropy = 0.0_dp

    call nasa7_dimensionless_properties( &
      species, temperature, cp_over_r, h_over_rt, s_over_r, ok)
    if (.not. ok) return
    ! Polynomial success does not imply that dimensional conversion is safe.
    ! Every guard below must leave failure until all outputs are committed.
    ok = .false.

    if (cp_over_r > finite_guard_limit / universal_gas_constant) return
    if (temperature > finite_guard_limit / universal_gas_constant) return
    gas_temperature = universal_gas_constant * temperature
    if (.not. ieee_is_finite(gas_temperature)) return
    if (gas_temperature >= 1.0_dp .and. abs(h_over_rt) >= 1.0_dp) then
      if (abs(h_over_rt) > finite_guard_limit / gas_temperature) return
    end if
    if (abs(s_over_r) > &
        finite_guard_limit / universal_gas_constant) return

    cp = universal_gas_constant * cp_over_r
    cv = cp - universal_gas_constant
    enthalpy = gas_temperature * h_over_rt
    if (.not. ieee_is_finite(enthalpy)) then
      cp = 0.0_dp
      cv = 0.0_dp
      enthalpy = 0.0_dp
      ok = .false.
      return
    end if
    if (abs(enthalpy) > finite_guard_limit) then
      cp = 0.0_dp
      cv = 0.0_dp
      enthalpy = 0.0_dp
      ok = .false.
      return
    end if
    if (enthalpy < 0.0_dp) then
      if (enthalpy < -finite_guard_limit + gas_temperature) then
        cp = 0.0_dp
        cv = 0.0_dp
        enthalpy = 0.0_dp
        ok = .false.
        return
      end if
    end if
    internal_energy = enthalpy - gas_temperature
    entropy = universal_gas_constant * s_over_r
    if (.not. all(ieee_is_finite( &
        [cp, cv, enthalpy, internal_energy, entropy]))) then
      cp = 0.0_dp
      cv = 0.0_dp
      enthalpy = 0.0_dp
      internal_energy = 0.0_dp
      entropy = 0.0_dp
      ok = .false.
      return
    end if
    if (cv <= 0.0_dp) then
      cp = 0.0_dp
      cv = 0.0_dp
      enthalpy = 0.0_dp
      internal_energy = 0.0_dp
      entropy = 0.0_dp
      ok = .false.
      return
    end if
    ok = .true.
  end subroutine nasa7_molar_properties

  subroutine nasa7_mass_properties( &
      species, temperature, cp, cv, enthalpy, internal_energy, entropy, ok)
    type(nasa7_species), intent(in) :: species
    real(dp), intent(in) :: temperature
    real(dp), intent(out) :: cp, cv, enthalpy, internal_energy, entropy
    logical, intent(out) :: ok

    real(dp) :: molar_cp, molar_cv, molar_enthalpy
    real(dp) :: molar_internal_energy, molar_entropy

    cp = 0.0_dp
    cv = 0.0_dp
    enthalpy = 0.0_dp
    internal_energy = 0.0_dp
    entropy = 0.0_dp

    call nasa7_molar_properties( &
      species, temperature, molar_cp, molar_cv, molar_enthalpy, &
      molar_internal_energy, molar_entropy, ok)
    if (.not. ok) return

    if (species%molecular_weight < maxval(abs([ &
        molar_cp, molar_cv, molar_enthalpy, molar_internal_energy, &
        molar_entropy])) / finite_guard_limit) then
      ok = .false.
      return
    end if

    cp = molar_cp / species%molecular_weight
    cv = molar_cv / species%molecular_weight
    enthalpy = molar_enthalpy / species%molecular_weight
    internal_energy = molar_internal_energy / species%molecular_weight
    entropy = molar_entropy / species%molecular_weight
    if (.not. all(ieee_is_finite( &
        [cp, cv, enthalpy, internal_energy, entropy]))) then
      cp = 0.0_dp
      cv = 0.0_dp
      enthalpy = 0.0_dp
      internal_energy = 0.0_dp
      entropy = 0.0_dp
      ok = .false.
      return
    end if
    ok = .true.
  end subroutine nasa7_mass_properties

end module nasa7_thermo_mod
