module nasa7_thermo_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  implicit none
  private

  real(dp), parameter, public :: universal_gas_constant = &
    8.31446261815324e3_dp ! J / (kmol K)

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

    valid = species%molecular_weight > 0.0_dp .and. &
      species%temperature_min > 0.0_dp .and. &
      species%temperature_mid > species%temperature_min .and. &
      species%temperature_max > species%temperature_mid .and. &
      all(ieee_is_finite(species%low_coefficients)) .and. &
      all(ieee_is_finite(species%high_coefficients))
  end function valid_nasa7_species

  real(dp) function nasa7_specific_gas_constant(species) result(gas_constant)
    type(nasa7_species), intent(in) :: species

    if (.not. valid_nasa7_species(species)) then
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
    real(dp) :: temperature_cubed, temperature_fourth

    cp_over_r = 0.0_dp
    h_over_rt = 0.0_dp
    s_over_r = 0.0_dp
    ok = .false.

    if (.not. valid_nasa7_species(species)) return
    if (temperature < species%temperature_min .or. &
        temperature > species%temperature_max) return

    if (temperature <= species%temperature_mid) then
      coefficients = species%low_coefficients
    else
      coefficients = species%high_coefficients
    end if

    temperature_squared = temperature * temperature
    temperature_cubed = temperature_squared * temperature
    temperature_fourth = temperature_squared * temperature_squared

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

    s_over_r = coefficients(1) * log(temperature) + &
      coefficients(2) * temperature + &
      0.5_dp * coefficients(3) * temperature_squared + &
      coefficients(4) * temperature_cubed / 3.0_dp + &
      0.25_dp * coefficients(5) * temperature_fourth + &
      coefficients(7)

    ok = ieee_is_finite(cp_over_r) .and. ieee_is_finite(h_over_rt) .and. &
      ieee_is_finite(s_over_r) .and. cp_over_r > 1.0_dp
  end subroutine nasa7_dimensionless_properties

  subroutine nasa7_molar_properties( &
      species, temperature, cp, cv, enthalpy, internal_energy, entropy, ok)
    type(nasa7_species), intent(in) :: species
    real(dp), intent(in) :: temperature
    real(dp), intent(out) :: cp, cv, enthalpy, internal_energy, entropy
    logical, intent(out) :: ok

    real(dp) :: cp_over_r, h_over_rt, s_over_r

    cp = 0.0_dp
    cv = 0.0_dp
    enthalpy = 0.0_dp
    internal_energy = 0.0_dp
    entropy = 0.0_dp

    call nasa7_dimensionless_properties( &
      species, temperature, cp_over_r, h_over_rt, s_over_r, ok)
    if (.not. ok) return

    cp = universal_gas_constant * cp_over_r
    cv = cp - universal_gas_constant
    enthalpy = universal_gas_constant * temperature * h_over_rt
    internal_energy = enthalpy - universal_gas_constant * temperature
    entropy = universal_gas_constant * s_over_r
    ok = cv > 0.0_dp .and. all(ieee_is_finite( &
      [cp, cv, enthalpy, internal_energy, entropy]))
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

    cp = molar_cp / species%molecular_weight
    cv = molar_cv / species%molecular_weight
    enthalpy = molar_enthalpy / species%molecular_weight
    internal_energy = molar_internal_energy / species%molecular_weight
    entropy = molar_entropy / species%molecular_weight
  end subroutine nasa7_mass_properties

end module nasa7_thermo_mod
