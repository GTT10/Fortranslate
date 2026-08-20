module mixture_transport_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: &
    nasa7_species, nasa7_mass_properties, nasa7_specific_gas_constant
  use mixture_thermo_mod, only: &
    valid_mixture_composition, mole_fractions_from_mass_fractions
  use transport_database_mod, only: &
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

  pure real(dp) function collision_integral_viscosity(reduced_temperature) &
      result(omega)
    real(dp), intent(in) :: reduced_temperature

    if (reduced_temperature <= 0.0_dp) then
      omega = -huge(1.0_dp)
      return
    end if
    omega = 1.16145_dp / reduced_temperature**0.14874_dp + &
      0.52487_dp / exp(0.77320_dp * reduced_temperature) + &
      2.16178_dp / exp(2.43787_dp * reduced_temperature)
  end function collision_integral_viscosity

  pure real(dp) function collision_integral_diffusion(reduced_temperature) &
      result(omega)
    real(dp), intent(in) :: reduced_temperature

    if (reduced_temperature <= 0.0_dp) then
      omega = -huge(1.0_dp)
      return
    end if
    omega = 1.06036_dp / reduced_temperature**0.15610_dp + &
      0.19300_dp / exp(0.47635_dp * reduced_temperature) + &
      1.03587_dp / exp(1.52996_dp * reduced_temperature) + &
      1.76474_dp / exp(3.89411_dp * reduced_temperature)
  end function collision_integral_diffusion

  subroutine pure_species_viscosities( &
      species, transport, temperature, viscosities, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: temperature
    real(dp), intent(out) :: viscosities(:)
    logical, intent(out) :: ok

    real(dp) :: reduced_temperature, omega
    integer :: k

    viscosities = 0.0_dp
    ok = compatible_transport_database(species, transport) .and. &
      size(viscosities) == size(species) .and. temperature > 0.0_dp
    if (.not. ok) return
    do k = 1, size(species)
      reduced_temperature = temperature / transport(k)%well_depth
      omega = collision_integral_viscosity(reduced_temperature)
      if (.not. ieee_is_finite(omega) .or. omega <= 0.0_dp) then
        ok = .false.
        return
      end if
      ! Chapman-Enskog dilute-gas expression. Molecular weight is numerically
      ! identical in kg/kmol and g/mol; sigma is in angstrom.
      viscosities(k) = viscosity_prefactor * &
        sqrt(species(k)%molecular_weight * temperature) / &
        (transport(k)%diameter**2 * omega)
    end do
    ok = all(ieee_is_finite(viscosities)) .and. all(viscosities > 0.0_dp)
  end subroutine pure_species_viscosities

  subroutine pure_species_thermal_conductivities( &
      species, transport, temperature, conductivities, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: temperature
    real(dp), intent(out) :: conductivities(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: viscosities(:)
    real(dp) :: cp, cv, enthalpy, internal_energy, entropy, r_species
    logical :: local_ok
    integer :: k

    conductivities = 0.0_dp
    allocate(viscosities(size(species)))
    call pure_species_viscosities( &
      species, transport, temperature, viscosities, ok)
    if (.not. ok .or. size(conductivities) /= size(species)) then
      ok = .false.
      return
    end if
    do k = 1, size(species)
      call nasa7_mass_properties( &
        species(k), temperature, cp, cv, enthalpy, internal_energy, entropy, &
        local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      r_species = nasa7_specific_gas_constant(species(k))
      if (.not. ieee_is_finite(r_species) .or. r_species <= 0.0_dp) then
        ok = .false.
        return
      end if
      ! Modified Eucken relation. This is a qualified dilute-gas subset; it
      ! does not include the full rotational/vibrational transport model used
      ! by Cantera or PelePhysics polynomial fits.
      conductivities(k) = viscosities(k) * (cp + 1.25_dp * r_species)
    end do
    ok = all(ieee_is_finite(conductivities)) .and. &
      all(conductivities > 0.0_dp)
  end subroutine pure_species_thermal_conductivities

  subroutine mixture_viscosity_wilke( &
      species, transport, mass_fractions, temperature, viscosity, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: mass_fractions(:), temperature
    real(dp), intent(out) :: viscosity
    logical, intent(out) :: ok

    real(dp), allocatable :: mole_fractions(:), pure_mu(:)
    real(dp) :: denominator, phi
    integer :: i, j

    viscosity = 0.0_dp
    ok = compatible_transport_database(species, transport) .and. &
      valid_mixture_composition(species, mass_fractions)
    if (.not. ok) return
    allocate(mole_fractions(size(species)), pure_mu(size(species)))
    call mole_fractions_from_mass_fractions( &
      species, mass_fractions, mole_fractions, ok)
    if (.not. ok) return
    call pure_species_viscosities( &
      species, transport, temperature, pure_mu, ok)
    if (.not. ok) return

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
        ok = .false.
        return
      end if
      viscosity = viscosity + mole_fractions(i) * pure_mu(i) / denominator
    end do
    ok = ieee_is_finite(viscosity) .and. viscosity > 0.0_dp
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

    conductivity = 0.0_dp
    ok = compatible_transport_database(species, transport) .and. &
      valid_mixture_composition(species, mass_fractions)
    if (.not. ok) return
    allocate(mole_fractions(size(species)), pure_lambda(size(species)))
    call mole_fractions_from_mass_fractions( &
      species, mass_fractions, mole_fractions, ok)
    if (.not. ok) return
    call pure_species_thermal_conductivities( &
      species, transport, temperature, pure_lambda, ok)
    if (.not. ok) return

    arithmetic_mean = sum(mole_fractions * pure_lambda)
    harmonic_denominator = sum(mole_fractions / pure_lambda)
    if (arithmetic_mean <= 0.0_dp .or. harmonic_denominator <= 0.0_dp) then
      ok = .false.
      return
    end if
    conductivity = 0.5_dp * &
      (arithmetic_mean + 1.0_dp / harmonic_denominator)
    ok = ieee_is_finite(conductivity) .and. conductivity > 0.0_dp
  end subroutine mixture_thermal_conductivity_mathur

  subroutine binary_diffusion_coefficients( &
      species, transport, temperature, pressure, binary_diffusion, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: temperature, pressure
    real(dp), intent(out) :: binary_diffusion(:, :)
    logical, intent(out) :: ok

    real(dp) :: sigma_ij, epsilon_ij, reduced_temperature, omega
    real(dp) :: pressure_atmospheres
    integer :: i, j

    binary_diffusion = 0.0_dp
    ok = compatible_transport_database(species, transport) .and. &
      size(binary_diffusion, 1) == size(species) .and. &
      size(binary_diffusion, 2) == size(species) .and. &
      temperature > 0.0_dp .and. pressure > 0.0_dp
    if (.not. ok) return
    pressure_atmospheres = pressure / standard_atmosphere
    do i = 1, size(species)
      do j = i + 1, size(species)
        sigma_ij = 0.5_dp * &
          (transport(i)%diameter + transport(j)%diameter)
        epsilon_ij = sqrt(transport(i)%well_depth * transport(j)%well_depth)
        reduced_temperature = temperature / epsilon_ij
        omega = collision_integral_diffusion(reduced_temperature)
        if (.not. ieee_is_finite(omega) .or. omega <= 0.0_dp) then
          ok = .false.
          return
        end if
        ! Chapman-Enskog binary coefficient. The conventional 0.001858
        ! expression returns cm^2/s; diffusion_prefactor includes 1e-4 to SI.
        binary_diffusion(i, j) = diffusion_prefactor * temperature**1.5_dp * &
          sqrt(1.0_dp / species(i)%molecular_weight + &
            1.0_dp / species(j)%molecular_weight) / &
          (pressure_atmospheres * sigma_ij**2 * omega)
        binary_diffusion(j, i) = binary_diffusion(i, j)
      end do
    end do
    ok = all(ieee_is_finite(binary_diffusion)) .and. &
      all(binary_diffusion >= 0.0_dp)
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
    real(dp), allocatable :: binary_diffusion(:, :)
    real(dp) :: denominator, total
    integer :: i, j, nspecies

    diffusion_coefficients = 0.0_dp
    nspecies = size(species)
    ok = compatible_transport_database(species, transport) .and. &
      size(diffusion_coefficients) == nspecies .and. &
      valid_mixture_composition(species, mass_fractions)
    if (.not. ok) return
    allocate(y_modified(nspecies), mole_fractions(nspecies))
    allocate(binary_diffusion(nspecies, nspecies))

    total = sum(mass_fractions)
    do i = 1, nspecies
      y_modified(i) = mass_fractions(i) + trace_fraction * &
        (total / real(nspecies, dp) - mass_fractions(i))
    end do
    y_modified = y_modified / sum(y_modified)
    call mole_fractions_from_mass_fractions( &
      species, y_modified, mole_fractions, ok)
    if (.not. ok) return
    call binary_diffusion_coefficients( &
      species, transport, temperature, pressure, binary_diffusion, ok)
    if (.not. ok) return

    do i = 1, nspecies
      denominator = 0.0_dp
      do j = 1, nspecies
        if (j /= i) denominator = denominator + &
          mole_fractions(j) / binary_diffusion(i, j)
      end do
      if (denominator <= 0.0_dp) then
        ok = .false.
        return
      end if
      diffusion_coefficients(i) = (1.0_dp - y_modified(i)) / denominator
    end do
    ok = all(ieee_is_finite(diffusion_coefficients)) .and. &
      all(diffusion_coefficients > 0.0_dp)
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

    viscosity = 0.0_dp
    conductivity = 0.0_dp
    diffusion_coefficients = 0.0_dp
    call mixture_viscosity_wilke( &
      species, transport, mass_fractions, temperature, viscosity, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call mixture_thermal_conductivity_mathur( &
      species, transport, mass_fractions, temperature, conductivity, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call mixture_averaged_diffusion_coefficients( &
      species, transport, mass_fractions, temperature, pressure, &
      diffusion_coefficients, local_ok)
    ok = local_ok
  end subroutine mixture_transport_coefficients

end module mixture_transport_mod
