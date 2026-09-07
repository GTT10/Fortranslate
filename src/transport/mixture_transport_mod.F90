module mixture_transport_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: &
    nasa7_species, valid_nasa7_species, nasa7_mass_properties, &
    nasa7_specific_gas_constant
  use mixture_thermo_mod, only: &
    valid_mixture_composition, mole_fractions_from_mass_fractions
  use gas_transport_mod, only: &
    gas_transport_species, compatible_transport_database
  implicit none
  private

  real(dp), parameter, public :: standard_atmosphere = 101325.0_dp
  real(dp), parameter :: trace_fraction = 1.0e-15_dp
  real(dp), parameter :: viscosity_prefactor = 2.6693e-6_dp
  real(dp), parameter :: diffusion_prefactor = 1.8580e-7_dp

  public :: collision_integral_viscosity
  public :: collision_integral_diffusion
  public :: pure_species_viscosities
  public :: pure_species_thermal_conductivities
  public :: mixture_viscosity_wilke
  public :: mixture_thermal_conductivity_mathur
  public :: binary_diffusion_coefficients
  public :: mixture_averaged_diffusion_coefficients
  public :: mixture_transport_coefficients

contains

  pure logical function finite_product_fits(left, right) result(fits)
    real(dp), intent(in) :: left, right

    fits = .false.
    if (.not. ieee_is_finite(left)) return
    if (.not. ieee_is_finite(right)) return
    if (abs(left) <= 1.0_dp .or. abs(right) <= 1.0_dp) then
      fits = .true.
      return
    end if
    fits = abs(left) <= huge(1.0_dp) / abs(right)
  end function finite_product_fits

  pure logical function finite_quotient_fits(numerator, denominator) &
      result(fits)
    real(dp), intent(in) :: numerator, denominator

    fits = .false.
    if (.not. ieee_is_finite(numerator)) return
    if (.not. ieee_is_finite(denominator)) return
    if (abs(denominator) <= 0.0_dp) return
    if (abs(denominator) >= 1.0_dp) then
      fits = .true.
      return
    end if
    fits = abs(numerator) <= huge(1.0_dp) * abs(denominator)
  end function finite_quotient_fits

  pure logical function finite_sum_fits(left, right) result(fits)
    real(dp), intent(in) :: left, right

    fits = .false.
    if (.not. ieee_is_finite(left)) return
    if (.not. ieee_is_finite(right)) return
    if (right > 0.0_dp) then
      if (left > huge(1.0_dp) - right) return
    else if (right < 0.0_dp) then
      if (left < -huge(1.0_dp) - right) return
    end if
    fits = .true.
  end function finite_sum_fits

  pure subroutine inverse_exponential_term( &
      coefficient, exponent, argument, value, ok)
    real(dp), intent(in) :: coefficient, exponent, argument
    real(dp), intent(out) :: value
    logical, intent(out) :: ok

    real(dp) :: exponential_argument

    value = 0.0_dp
    ok = .false.
    if (.not. all(ieee_is_finite([coefficient, exponent, argument]))) return
    if (coefficient < 0.0_dp .or. exponent <= 0.0_dp .or. &
        argument <= 0.0_dp) return
    if (.not. finite_product_fits(exponent, argument)) then
      ok = .true.
      return
    end if
    exponential_argument = exponent * argument
    value = coefficient * exp(-exponential_argument)
    if (.not. ieee_is_finite(value)) then
      value = 0.0_dp
      return
    end if
    if (value < 0.0_dp) then
      value = 0.0_dp
      return
    end if
    ok = .true.
  end subroutine inverse_exponential_term

  logical function valid_transport_inputs(species, transport) result(valid)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    integer :: k

    valid = .false.
    if (.not. compatible_transport_database(species, transport)) return
    do k = 1, size(species)
      if (.not. valid_nasa7_species(species(k))) return
    end do
    valid = .true.
  end function valid_transport_inputs

  pure real(dp) function collision_integral_viscosity(reduced_temperature) &
      result(omega)
    real(dp), intent(in) :: reduced_temperature

    real(dp) :: exponential_term_one, exponential_term_two
    logical :: term_ok

    if (.not. ieee_is_finite(reduced_temperature)) then
      omega = -huge(1.0_dp)
      return
    end if
    if (reduced_temperature <= 0.0_dp) then
      omega = -huge(1.0_dp)
      return
    end if
    if (reduced_temperature <= &
        (log(huge(1.0_dp)) - 1.0_dp) / 2.43787_dp) then
      omega = 1.16145_dp / reduced_temperature**0.14874_dp + &
        0.52487_dp / exp(0.77320_dp * reduced_temperature) + &
        2.16178_dp / exp(2.43787_dp * reduced_temperature)
    else
      call inverse_exponential_term( &
        0.52487_dp, 0.77320_dp, reduced_temperature, &
        exponential_term_one, term_ok)
      if (.not. term_ok) then
        omega = -huge(1.0_dp)
        return
      end if
      call inverse_exponential_term( &
        2.16178_dp, 2.43787_dp, reduced_temperature, &
        exponential_term_two, term_ok)
      if (.not. term_ok) then
        omega = -huge(1.0_dp)
        return
      end if
      omega = 1.16145_dp / reduced_temperature**0.14874_dp + &
        exponential_term_one + exponential_term_two
    end if
    if (.not. ieee_is_finite(omega)) then
      omega = -huge(1.0_dp)
      return
    end if
    if (omega <= 0.0_dp) omega = -huge(1.0_dp)
  end function collision_integral_viscosity

  pure real(dp) function collision_integral_diffusion(reduced_temperature) &
      result(omega)
    real(dp), intent(in) :: reduced_temperature

    real(dp) :: exponential_term_one, exponential_term_two
    real(dp) :: exponential_term_three
    logical :: term_ok

    if (.not. ieee_is_finite(reduced_temperature)) then
      omega = -huge(1.0_dp)
      return
    end if
    if (reduced_temperature <= 0.0_dp) then
      omega = -huge(1.0_dp)
      return
    end if
    if (reduced_temperature <= &
        (log(huge(1.0_dp)) - 1.0_dp) / 3.89411_dp) then
      omega = 1.06036_dp / reduced_temperature**0.15610_dp + &
        0.19300_dp / exp(0.47635_dp * reduced_temperature) + &
        1.03587_dp / exp(1.52996_dp * reduced_temperature) + &
        1.76474_dp / exp(3.89411_dp * reduced_temperature)
    else
      call inverse_exponential_term( &
        0.19300_dp, 0.47635_dp, reduced_temperature, &
        exponential_term_one, term_ok)
      if (.not. term_ok) then
        omega = -huge(1.0_dp)
        return
      end if
      call inverse_exponential_term( &
        1.03587_dp, 1.52996_dp, reduced_temperature, &
        exponential_term_two, term_ok)
      if (.not. term_ok) then
        omega = -huge(1.0_dp)
        return
      end if
      call inverse_exponential_term( &
        1.76474_dp, 3.89411_dp, reduced_temperature, &
        exponential_term_three, term_ok)
      if (.not. term_ok) then
        omega = -huge(1.0_dp)
        return
      end if
      omega = 1.06036_dp / reduced_temperature**0.15610_dp + &
        exponential_term_one + exponential_term_two + exponential_term_three
    end if
    if (.not. ieee_is_finite(omega)) then
      omega = -huge(1.0_dp)
      return
    end if
    if (omega <= 0.0_dp) omega = -huge(1.0_dp)
  end function collision_integral_diffusion

  subroutine pure_species_viscosities( &
      species, transport, temperature, viscosities, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: temperature
    real(dp), intent(out) :: viscosities(:)
    logical, intent(out) :: ok

    real(dp) :: reduced_temperature, omega
    real(dp) :: molecular_temperature, diameter_squared, denominator
    real(dp) :: numerator, sqrt_molecular_temperature
    real(dp), allocatable :: candidate(:)
    integer :: k

    viscosities = 0.0_dp
    ok = .false.
    if (.not. valid_transport_inputs(species, transport)) return
    if (size(viscosities) /= size(species)) return
    if (.not. ieee_is_finite(temperature)) return
    if (temperature <= 0.0_dp) return
    allocate(candidate(size(species)))
    candidate = 0.0_dp
    do k = 1, size(species)
      if (.not. finite_product_fits( &
          temperature, species(k)%molecular_weight)) return
      molecular_temperature = species(k)%molecular_weight * temperature
      if (.not. ieee_is_finite(molecular_temperature)) return
      if (molecular_temperature <= 0.0_dp) return
      if (.not. finite_product_fits( &
          transport(k)%diameter, transport(k)%diameter)) return
      diameter_squared = transport(k)%diameter**2
      if (.not. ieee_is_finite(diameter_squared)) return
      if (diameter_squared <= 0.0_dp) return
      if (.not. finite_quotient_fits( &
          temperature, transport(k)%well_depth)) return
      reduced_temperature = temperature / transport(k)%well_depth
      if (.not. ieee_is_finite(reduced_temperature)) return
      if (reduced_temperature <= 0.0_dp) return
      omega = collision_integral_viscosity(reduced_temperature)
      if (.not. ieee_is_finite(omega)) return
      if (omega <= 0.0_dp) return
      ! Chapman-Enskog dilute-gas expression. Molecular weight is numerically
      ! identical in kg/kmol and g/mol; sigma is in angstrom.
      if (.not. finite_product_fits(diameter_squared, omega)) return
      denominator = diameter_squared * omega
      if (.not. ieee_is_finite(denominator)) return
      if (denominator <= 0.0_dp) return
      sqrt_molecular_temperature = sqrt(molecular_temperature)
      if (.not. ieee_is_finite(sqrt_molecular_temperature)) return
      numerator = viscosity_prefactor * sqrt_molecular_temperature
      if (.not. ieee_is_finite(numerator)) return
      if (numerator <= 0.0_dp) return
      if (.not. finite_quotient_fits(numerator, denominator)) return
      candidate(k) = viscosity_prefactor * &
        sqrt(species(k)%molecular_weight * temperature) / &
        (transport(k)%diameter**2 * omega)
      if (.not. ieee_is_finite(candidate(k))) return
      if (candidate(k) <= 0.0_dp) return
    end do
    viscosities = candidate
    ok = .true.
  end subroutine pure_species_viscosities

  subroutine pure_species_thermal_conductivities( &
      species, transport, temperature, conductivities, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: temperature
    real(dp), intent(out) :: conductivities(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: viscosities(:), candidate(:)
    real(dp) :: cp, cv, enthalpy, internal_energy, entropy, r_species
    real(dp) :: conductivity_factor, rotational_term
    logical :: local_ok
    integer :: k

    conductivities = 0.0_dp
    ok = .false.
    if (.not. valid_transport_inputs(species, transport)) return
    if (size(conductivities) /= size(species)) return
    if (.not. ieee_is_finite(temperature)) return
    if (temperature <= 0.0_dp) return
    allocate(viscosities(size(species)))
    allocate(candidate(size(species)))
    candidate = 0.0_dp
    call pure_species_viscosities( &
      species, transport, temperature, viscosities, ok)
    if (.not. ok) return
    do k = 1, size(species)
      call nasa7_mass_properties( &
        species(k), temperature, cp, cv, enthalpy, internal_energy, entropy, &
        local_ok)
      if (.not. local_ok) return
      if (.not. all(ieee_is_finite([cp, cv, enthalpy, internal_energy, entropy]))) &
        return
      r_species = nasa7_specific_gas_constant(species(k))
      if (.not. ieee_is_finite(r_species)) return
      if (r_species <= 0.0_dp) return
      ! Modified Eucken relation. This is a qualified dilute-gas subset; it
      ! does not include the full rotational/vibrational transport model used
      ! by Cantera or PelePhysics polynomial fits.
      if (.not. finite_product_fits(1.25_dp, r_species)) return
      rotational_term = 1.25_dp * r_species
      if (.not. ieee_is_finite(rotational_term)) return
      if (.not. finite_sum_fits(cp, rotational_term)) return
      conductivity_factor = cp + rotational_term
      if (.not. ieee_is_finite(conductivity_factor)) return
      if (conductivity_factor <= 0.0_dp) return
      if (.not. finite_product_fits( &
          viscosities(k), conductivity_factor)) return
      candidate(k) = viscosities(k) * (cp + 1.25_dp * r_species)
      if (.not. ieee_is_finite(candidate(k))) return
      if (candidate(k) <= 0.0_dp) return
    end do
    conductivities = candidate
    ok = .true.
  end subroutine pure_species_thermal_conductivities

  subroutine mixture_viscosity_wilke( &
      species, transport, mass_fractions, temperature, viscosity, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: mass_fractions(:), temperature
    real(dp), intent(out) :: viscosity
    logical, intent(out) :: ok

    real(dp), allocatable :: mole_fractions(:), pure_mu(:)
    real(dp) :: denominator, phi, denominator_term
    real(dp) :: molecular_viscosity_ratio, molecular_weight_ratio
    real(dp) :: sqrt_viscosity_ratio, fourth_root_weight_ratio
    real(dp) :: phi_numerator, phi_denominator, numerator, contribution
    real(dp) :: phi_product, phi_base, weight_base, weight_argument
    real(dp) :: candidate_viscosity, guard_denominator, guard_viscosity
    integer :: i, j

    viscosity = 0.0_dp
    ok = .false.
    if (.not. valid_transport_inputs(species, transport)) return
    if (.not. valid_mixture_composition(species, mass_fractions)) return
    if (.not. ieee_is_finite(temperature)) return
    if (temperature <= 0.0_dp) return
    allocate(mole_fractions(size(species)), pure_mu(size(species)))
    call mole_fractions_from_mass_fractions( &
      species, mass_fractions, mole_fractions, ok)
    if (.not. ok) return
    if (.not. all(ieee_is_finite(mole_fractions))) then
      ok = .false.
      return
    end if
    if (any(mole_fractions < 0.0_dp)) then
      ok = .false.
      return
    end if
    call pure_species_viscosities( &
      species, transport, temperature, pure_mu, ok)
    if (.not. ok) return
    if (.not. all(ieee_is_finite(pure_mu))) return
    if (any(pure_mu <= 0.0_dp)) return

    candidate_viscosity = 0.0_dp
    guard_viscosity = 0.0_dp
    do i = 1, size(species)
      denominator = 0.0_dp
      guard_denominator = 0.0_dp
      do j = 1, size(species)
        if (.not. finite_quotient_fits(pure_mu(i), pure_mu(j))) return
        molecular_viscosity_ratio = pure_mu(i) / pure_mu(j)
        if (.not. ieee_is_finite(molecular_viscosity_ratio)) return
        if (molecular_viscosity_ratio <= 0.0_dp) return
        sqrt_viscosity_ratio = sqrt(molecular_viscosity_ratio)
        if (.not. ieee_is_finite(sqrt_viscosity_ratio)) return
        if (.not. finite_quotient_fits( &
            species(j)%molecular_weight, &
            species(i)%molecular_weight)) return
        molecular_weight_ratio = species(j)%molecular_weight / &
          species(i)%molecular_weight
        if (.not. ieee_is_finite(molecular_weight_ratio)) return
        if (molecular_weight_ratio <= 0.0_dp) return
        fourth_root_weight_ratio = molecular_weight_ratio**0.25_dp
        if (.not. ieee_is_finite(fourth_root_weight_ratio)) return
        if (.not. finite_product_fits( &
            sqrt_viscosity_ratio, fourth_root_weight_ratio)) return
        phi_product = sqrt_viscosity_ratio * fourth_root_weight_ratio
        if (.not. ieee_is_finite(phi_product)) return
        if (phi_product < 0.0_dp) return
        if (.not. finite_sum_fits(1.0_dp, phi_product)) return
        phi_base = 1.0_dp + phi_product
        if (.not. ieee_is_finite(phi_base)) return
        if (phi_base <= 0.0_dp) return
        if (phi_base > sqrt(huge(1.0_dp))) return
        phi_numerator = phi_base**2
        if (.not. ieee_is_finite(phi_numerator)) return
        if (phi_numerator <= 0.0_dp) return
        if (.not. finite_quotient_fits( &
            species(i)%molecular_weight, &
            species(j)%molecular_weight)) return
        molecular_weight_ratio = species(i)%molecular_weight / &
          species(j)%molecular_weight
        if (.not. ieee_is_finite(molecular_weight_ratio)) return
        if (molecular_weight_ratio <= 0.0_dp) return
        if (.not. finite_sum_fits(1.0_dp, molecular_weight_ratio)) return
        weight_base = 1.0_dp + molecular_weight_ratio
        if (.not. ieee_is_finite(weight_base)) return
        if (weight_base <= 0.0_dp) return
        if (.not. finite_product_fits(8.0_dp, weight_base)) return
        weight_argument = 8.0_dp * weight_base
        if (.not. ieee_is_finite(weight_argument)) return
        if (weight_argument <= 0.0_dp) return
        phi_denominator = sqrt(weight_argument)
        if (.not. ieee_is_finite(phi_denominator)) return
        if (phi_denominator <= 0.0_dp) return
        if (.not. finite_quotient_fits( &
            phi_numerator, phi_denominator)) return
        phi = phi_numerator / phi_denominator
        if (.not. ieee_is_finite(phi)) return
        if (phi <= 0.0_dp) return
        if (mole_fractions(j) > 0.0_dp) then
          if (.not. finite_product_fits(mole_fractions(j), phi)) return
        end if
        denominator_term = mole_fractions(j) * phi
        if (.not. ieee_is_finite(denominator_term)) return
        if (denominator_term < 0.0_dp) return
        if (.not. finite_sum_fits(guard_denominator, denominator_term)) return
        guard_denominator = guard_denominator + denominator_term
        if (.not. ieee_is_finite(guard_denominator)) return

        phi = (1.0_dp + sqrt(pure_mu(i) / pure_mu(j)) * &
          (species(j)%molecular_weight / species(i)%molecular_weight)**0.25_dp)**2 / &
          sqrt(8.0_dp * (1.0_dp + species(i)%molecular_weight / &
            species(j)%molecular_weight))
        if (.not. ieee_is_finite(phi)) return
        if (phi <= 0.0_dp) return
        if (.not. finite_product_fits(mole_fractions(j), phi)) return
        denominator_term = mole_fractions(j) * phi
        if (.not. ieee_is_finite(denominator_term)) return
        if (.not. finite_sum_fits(denominator, denominator_term)) return
        denominator = denominator + mole_fractions(j) * phi
        if (.not. ieee_is_finite(denominator)) return
      end do
      if (denominator <= 0.0_dp) return
      if (mole_fractions(i) > 0.0_dp) then
        if (.not. finite_product_fits( &
            mole_fractions(i), pure_mu(i))) return
      end if
      numerator = mole_fractions(i) * pure_mu(i)
      if (.not. ieee_is_finite(numerator)) return
      if (numerator < 0.0_dp) return
      if (.not. finite_quotient_fits(numerator, denominator)) return
      contribution = numerator / denominator
      if (.not. ieee_is_finite(contribution)) return
      if (contribution < 0.0_dp) return
      if (.not. finite_sum_fits(guard_viscosity, contribution)) return
      guard_viscosity = guard_viscosity + contribution
      if (.not. ieee_is_finite(guard_viscosity)) return
      candidate_viscosity = candidate_viscosity + &
        mole_fractions(i) * pure_mu(i) / denominator
      if (.not. ieee_is_finite(candidate_viscosity)) return
    end do
    if (candidate_viscosity <= 0.0_dp) return

    ! The loop above is a fail-closed preflight.  Re-evaluate the accepted
    ! state with the established Wilke operation order so optimized builds
    ! retain the frozen checkpoint rounding contract.
    viscosity = 0.0_dp
    do i = 1, size(species)
      denominator = 0.0_dp
      do j = 1, size(species)
        phi = (1.0_dp + sqrt(pure_mu(i) / pure_mu(j)) * &
          (species(j)%molecular_weight / species(i)%molecular_weight)**0.25_dp)**2 / &
          sqrt(8.0_dp * (1.0_dp + species(i)%molecular_weight / &
            species(j)%molecular_weight))
        denominator = denominator + mole_fractions(j) * phi
      end do
      if (denominator <= 0.0_dp) then
        viscosity = 0.0_dp
        return
      end if
      viscosity = viscosity + mole_fractions(i) * pure_mu(i) / denominator
    end do
    if (.not. ieee_is_finite(viscosity) .or. viscosity <= 0.0_dp) then
      viscosity = 0.0_dp
      return
    end if
    ok = .true.
  end subroutine mixture_viscosity_wilke

  subroutine mixture_thermal_conductivity_mathur( &
      species, transport, mass_fractions, temperature, conductivity, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: mass_fractions(:), temperature
    real(dp), intent(out) :: conductivity
    logical, intent(out) :: ok

    real(dp), allocatable :: mole_fractions(:), pure_lambda(:)
    real(dp) :: arithmetic_mean, harmonic_denominator
    real(dp) :: arithmetic_term, harmonic_term, inverse_harmonic
    real(dp) :: guard_arithmetic_mean, guard_harmonic_denominator
    integer :: k

    conductivity = 0.0_dp
    ok = .false.
    if (.not. valid_transport_inputs(species, transport)) return
    if (.not. valid_mixture_composition(species, mass_fractions)) return
    if (.not. ieee_is_finite(temperature)) return
    if (temperature <= 0.0_dp) return
    allocate(mole_fractions(size(species)), pure_lambda(size(species)))
    call mole_fractions_from_mass_fractions( &
      species, mass_fractions, mole_fractions, ok)
    if (.not. ok) return
    if (.not. all(ieee_is_finite(mole_fractions))) then
      ok = .false.
      return
    end if
    if (any(mole_fractions < 0.0_dp)) then
      ok = .false.
      return
    end if
    call pure_species_thermal_conductivities( &
      species, transport, temperature, pure_lambda, ok)
    if (.not. ok) return
    if (.not. all(ieee_is_finite(pure_lambda))) return
    if (any(pure_lambda <= 0.0_dp)) return

    guard_arithmetic_mean = 0.0_dp
    guard_harmonic_denominator = 0.0_dp
    do k = 1, size(species)
      if (mole_fractions(k) > 0.0_dp) then
        if (.not. finite_product_fits( &
            mole_fractions(k), pure_lambda(k))) return
      end if
      arithmetic_term = mole_fractions(k) * pure_lambda(k)
      if (.not. finite_quotient_fits( &
          mole_fractions(k), pure_lambda(k))) return
      harmonic_term = mole_fractions(k) / pure_lambda(k)
      if (.not. ieee_is_finite(arithmetic_term)) return
      if (.not. ieee_is_finite(harmonic_term)) return
      if (arithmetic_term < 0.0_dp) return
      if (harmonic_term < 0.0_dp) return
      if (.not. finite_sum_fits( &
          guard_arithmetic_mean, arithmetic_term)) return
      if (.not. finite_sum_fits( &
          guard_harmonic_denominator, harmonic_term)) return
      guard_arithmetic_mean = guard_arithmetic_mean + arithmetic_term
      guard_harmonic_denominator = guard_harmonic_denominator + harmonic_term
      if (.not. ieee_is_finite(guard_arithmetic_mean)) return
      if (.not. ieee_is_finite(guard_harmonic_denominator)) return
    end do
    if (guard_arithmetic_mean <= 0.0_dp) return
    if (guard_harmonic_denominator <= 0.0_dp) return
    arithmetic_mean = sum(mole_fractions * pure_lambda)
    harmonic_denominator = sum(mole_fractions / pure_lambda)
    if (.not. ieee_is_finite(arithmetic_mean)) return
    if (.not. ieee_is_finite(harmonic_denominator)) return
    if (arithmetic_mean <= 0.0_dp) return
    if (harmonic_denominator <= 0.0_dp) return
    if (.not. finite_quotient_fits(1.0_dp, harmonic_denominator)) return
    inverse_harmonic = 1.0_dp / harmonic_denominator
    if (.not. ieee_is_finite(inverse_harmonic)) return
    if (inverse_harmonic <= 0.0_dp) return
    if (inverse_harmonic > huge(1.0_dp) - arithmetic_mean) return
    if (.not. finite_sum_fits(arithmetic_mean, inverse_harmonic)) return
    conductivity = 0.5_dp * &
      (arithmetic_mean + 1.0_dp / harmonic_denominator)
    if (.not. ieee_is_finite(conductivity)) then
      conductivity = 0.0_dp
      return
    end if
    if (conductivity <= 0.0_dp) then
      conductivity = 0.0_dp
      return
    end if
    ok = .true.
  end subroutine mixture_thermal_conductivity_mathur

  subroutine binary_diffusion_coefficients( &
      species, transport, temperature, pressure, binary_diffusion, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: temperature, pressure
    real(dp), intent(out) :: binary_diffusion(:, :)
    logical, intent(out) :: ok

    real(dp) :: sigma_ij, epsilon_ij, reduced_temperature, omega
    real(dp) :: pressure_atmospheres, sigma_squared, epsilon_product
    real(dp) :: inverse_weight_i, inverse_weight_j, inverse_weight_sum
    real(dp) :: temperature_power, sqrt_weight_sum
    real(dp) :: numerator_factor, numerator, denominator, coefficient
    real(dp), allocatable :: candidate(:, :)
    integer :: i, j

    binary_diffusion = 0.0_dp
    ok = .false.
    if (.not. valid_transport_inputs(species, transport)) return
    if (size(binary_diffusion, 1) /= size(species)) return
    if (size(binary_diffusion, 2) /= size(species)) return
    if (.not. ieee_is_finite(temperature)) return
    if (temperature <= 0.0_dp) return
    if (.not. ieee_is_finite(pressure)) return
    if (pressure <= 0.0_dp) return
    allocate(candidate(size(species), size(species)))
    candidate = 0.0_dp
    pressure_atmospheres = pressure / standard_atmosphere
    if (.not. ieee_is_finite(pressure_atmospheres)) return
    if (pressure_atmospheres <= 0.0_dp) return
    if (temperature > &
        (0.5_dp * huge(1.0_dp))**(2.0_dp / 3.0_dp)) return
    temperature_power = temperature**1.5_dp
    if (.not. ieee_is_finite(temperature_power)) return
    if (temperature_power <= 0.0_dp) return
    do i = 1, size(species)
      do j = i + 1, size(species)
        if (.not. finite_product_fits( &
            transport(i)%diameter, transport(i)%diameter)) return
        if (.not. finite_product_fits( &
            transport(j)%diameter, transport(j)%diameter)) return
        if (.not. finite_sum_fits( &
            transport(i)%diameter, transport(j)%diameter)) return
        sigma_ij = 0.5_dp * &
          (transport(i)%diameter + transport(j)%diameter)
        if (.not. ieee_is_finite(sigma_ij)) return
        if (sigma_ij <= 0.0_dp) return
        if (.not. finite_product_fits(sigma_ij, sigma_ij)) return
        sigma_squared = sigma_ij**2
        if (.not. ieee_is_finite(sigma_squared)) return
        if (sigma_squared <= 0.0_dp) return
        if (.not. finite_product_fits( &
            transport(i)%well_depth, transport(j)%well_depth)) return
        epsilon_product = transport(i)%well_depth * transport(j)%well_depth
        if (.not. ieee_is_finite(epsilon_product)) return
        if (epsilon_product <= 0.0_dp) return
        epsilon_ij = sqrt(epsilon_product)
        if (.not. ieee_is_finite(epsilon_ij)) return
        if (epsilon_ij <= 0.0_dp) return
        if (.not. finite_quotient_fits(temperature, epsilon_ij)) return
        reduced_temperature = temperature / epsilon_ij
        if (.not. ieee_is_finite(reduced_temperature)) return
        if (reduced_temperature <= 0.0_dp) return
        omega = collision_integral_diffusion(reduced_temperature)
        if (.not. ieee_is_finite(omega)) return
        if (omega <= 0.0_dp) return
        if (.not. finite_quotient_fits( &
            1.0_dp, species(i)%molecular_weight)) return
        inverse_weight_i = 1.0_dp / species(i)%molecular_weight
        if (.not. ieee_is_finite(inverse_weight_i)) return
        if (inverse_weight_i <= 0.0_dp) return
        if (.not. finite_quotient_fits( &
            1.0_dp, species(j)%molecular_weight)) return
        inverse_weight_j = 1.0_dp / species(j)%molecular_weight
        if (.not. ieee_is_finite(inverse_weight_j)) return
        if (inverse_weight_j <= 0.0_dp) return
        if (inverse_weight_j > huge(1.0_dp) - inverse_weight_i) return
        inverse_weight_sum = inverse_weight_i + inverse_weight_j
        if (.not. ieee_is_finite(inverse_weight_sum)) return
        if (inverse_weight_sum <= 0.0_dp) return
        sqrt_weight_sum = sqrt(inverse_weight_sum)
        if (.not. ieee_is_finite(sqrt_weight_sum)) return
        if (sqrt_weight_sum <= 0.0_dp) return
        ! Chapman-Enskog binary coefficient. The conventional 0.001858
        ! expression returns cm^2/s; diffusion_prefactor includes 1e-4 to SI.
        numerator_factor = diffusion_prefactor * temperature_power
        if (.not. ieee_is_finite(numerator_factor)) return
        if (.not. finite_product_fits( &
            numerator_factor, sqrt_weight_sum)) return
        numerator = numerator_factor * sqrt_weight_sum
        if (.not. ieee_is_finite(numerator)) return
        if (numerator <= 0.0_dp) return
        if (.not. finite_product_fits( &
            pressure_atmospheres, sigma_squared)) return
        denominator = pressure_atmospheres * sigma_squared
        if (.not. ieee_is_finite(denominator)) return
        if (denominator <= 0.0_dp) return
        if (.not. finite_product_fits(denominator, omega)) return
        denominator = denominator * omega
        if (.not. ieee_is_finite(denominator)) return
        if (denominator <= 0.0_dp) return
        if (.not. finite_quotient_fits(numerator, denominator)) return
        coefficient = numerator / denominator
        if (.not. ieee_is_finite(coefficient)) return
        if (coefficient <= 0.0_dp) return
        candidate(i, j) = diffusion_prefactor * temperature**1.5_dp * &
          sqrt(1.0_dp / species(i)%molecular_weight + &
            1.0_dp / species(j)%molecular_weight) / &
          (pressure_atmospheres * sigma_ij**2 * omega)
        if (.not. ieee_is_finite(candidate(i, j))) return
        if (candidate(i, j) <= 0.0_dp) return
        candidate(j, i) = candidate(i, j)
      end do
    end do
    binary_diffusion = candidate
    ok = .true.
  end subroutine binary_diffusion_coefficients

  subroutine mixture_averaged_diffusion_coefficients( &
      species, transport, mass_fractions, temperature, pressure, &
      diffusion_coefficients, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: mass_fractions(:), temperature, pressure
    real(dp), intent(out) :: diffusion_coefficients(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: y_modified(:), mole_fractions(:)
    real(dp), allocatable :: binary_diffusion(:, :), candidate(:)
    real(dp) :: denominator, total, modified_total, mass_fraction_mean
    real(dp) :: mass_difference, correction, numerator, coefficient, term
    integer :: i, j, nspecies

    diffusion_coefficients = 0.0_dp
    nspecies = size(species)
    ok = .false.
    if (.not. valid_transport_inputs(species, transport)) return
    if (size(diffusion_coefficients) /= nspecies) return
    if (.not. valid_mixture_composition(species, mass_fractions)) return
    if (.not. ieee_is_finite(temperature)) return
    if (temperature <= 0.0_dp) return
    if (.not. ieee_is_finite(pressure)) return
    if (pressure <= 0.0_dp) return
    allocate(y_modified(nspecies), mole_fractions(nspecies))
    allocate(binary_diffusion(nspecies, nspecies))
    allocate(candidate(nspecies))
    candidate = 0.0_dp

    total = sum(mass_fractions)
    if (.not. ieee_is_finite(total)) return
    if (total <= 0.0_dp) return
    mass_fraction_mean = total / real(nspecies, dp)
    if (.not. ieee_is_finite(mass_fraction_mean)) return
    if (mass_fraction_mean <= 0.0_dp) return
    do i = 1, nspecies
      mass_difference = mass_fraction_mean - mass_fractions(i)
      if (.not. ieee_is_finite(mass_difference)) return
      correction = trace_fraction * mass_difference
      if (.not. ieee_is_finite(correction)) return
      y_modified(i) = mass_fractions(i) + correction
      if (.not. ieee_is_finite(y_modified(i))) return
    end do
    modified_total = sum(y_modified)
    if (.not. ieee_is_finite(modified_total)) return
    if (modified_total <= 0.0_dp) return
    do i = 1, nspecies
      y_modified(i) = y_modified(i) / modified_total
      if (.not. ieee_is_finite(y_modified(i))) return
    end do
    call mole_fractions_from_mass_fractions( &
      species, y_modified, mole_fractions, ok)
    if (.not. ok) return
    if (.not. all(ieee_is_finite(mole_fractions))) return
    if (any(mole_fractions < 0.0_dp)) return
    call binary_diffusion_coefficients( &
      species, transport, temperature, pressure, binary_diffusion, ok)
    if (.not. ok) return
    if (.not. all(ieee_is_finite(binary_diffusion))) return

    do i = 1, nspecies
      denominator = 0.0_dp
      do j = 1, nspecies
        if (j /= i) then
          if (.not. ieee_is_finite(binary_diffusion(i, j))) return
          if (binary_diffusion(i, j) <= 0.0_dp) return
          if (.not. finite_quotient_fits( &
              mole_fractions(j), binary_diffusion(i, j))) return
          term = mole_fractions(j) / binary_diffusion(i, j)
          if (.not. ieee_is_finite(term)) return
          if (term < 0.0_dp) return
          if (term > huge(1.0_dp) - denominator) return
          denominator = denominator + term
          if (.not. ieee_is_finite(denominator)) return
        end if
      end do
      if (denominator <= 0.0_dp) return
      numerator = 1.0_dp - y_modified(i)
      if (.not. ieee_is_finite(numerator)) return
      if (numerator <= 0.0_dp) return
      if (.not. finite_quotient_fits(numerator, denominator)) return
      coefficient = numerator / denominator
      if (.not. ieee_is_finite(coefficient)) return
      if (coefficient <= 0.0_dp) return
      candidate(i) = coefficient
    end do
    diffusion_coefficients = candidate
    ok = .true.
  end subroutine mixture_averaged_diffusion_coefficients

  subroutine mixture_transport_coefficients( &
      species, transport, mass_fractions, temperature, pressure, viscosity, &
      conductivity, diffusion_coefficients, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: mass_fractions(:), temperature, pressure
    real(dp), intent(out) :: viscosity, conductivity
    real(dp), intent(out) :: diffusion_coefficients(:)
    logical, intent(out) :: ok

    logical :: local_ok
    real(dp) :: candidate_viscosity, candidate_conductivity
    real(dp), allocatable :: candidate_diffusion(:)

    viscosity = 0.0_dp
    conductivity = 0.0_dp
    diffusion_coefficients = 0.0_dp
    ok = .false.
    if (.not. valid_transport_inputs(species, transport)) return
    if (size(diffusion_coefficients) /= size(species)) return
    if (.not. valid_mixture_composition(species, mass_fractions)) return
    if (.not. ieee_is_finite(temperature)) return
    if (temperature <= 0.0_dp) return
    if (.not. ieee_is_finite(pressure)) return
    if (pressure <= 0.0_dp) return
    allocate(candidate_diffusion(size(species)))
    candidate_viscosity = 0.0_dp
    candidate_conductivity = 0.0_dp
    candidate_diffusion = 0.0_dp
    call mixture_viscosity_wilke( &
      species, transport, mass_fractions, temperature, &
      candidate_viscosity, local_ok)
    if (.not. local_ok) then
      return
    end if
    call mixture_thermal_conductivity_mathur( &
      species, transport, mass_fractions, temperature, &
      candidate_conductivity, local_ok)
    if (.not. local_ok) then
      return
    end if
    call mixture_averaged_diffusion_coefficients( &
      species, transport, mass_fractions, temperature, pressure, &
      candidate_diffusion, local_ok)
    if (.not. local_ok) return
    viscosity = candidate_viscosity
    conductivity = candidate_conductivity
    diffusion_coefficients = candidate_diffusion
    ok = .true.
  end subroutine mixture_transport_coefficients

end module mixture_transport_mod
