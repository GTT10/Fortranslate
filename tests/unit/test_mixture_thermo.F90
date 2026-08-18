program test_mixture_thermo
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: &
    load_gri30_thermo_subset, gri_h2_index, gri_o2_index, &
    gri_h2o_index, gri_n2_index
  use mixture_thermo_mod, only: &
    valid_mixture_composition, mixture_mass_properties, mixture_pressure, &
    mixture_density, mixture_sound_speed, temperature_from_internal_energy
  implicit none

  type(nasa7_species), allocatable :: species(:)
  real(dp) :: mass_fractions(4), invalid_mass_fractions(4)
  real(dp) :: molecular_weight, gas_constant, cp, cv, gamma
  real(dp) :: enthalpy, internal_energy, entropy, pressure, density
  real(dp) :: sound_speed, recovered_temperature, target_temperature
  logical :: ok
  integer :: iterations, i
  real(dp), parameter :: temperatures(3) = [300.0_dp, 1200.0_dp, 2500.0_dp]

  call load_gri30_thermo_subset(species, ok)
  if (.not. ok) error stop "Failed to load mixture species"

  mass_fractions = 0.0_dp
  mass_fractions(gri_o2_index) = 0.23291751145757963_dp
  mass_fractions(gri_n2_index) = 0.7670824885424203_dp
  if (.not. valid_mixture_composition(species, mass_fractions)) then
    error stop "Valid air composition was rejected"
  end if

  call mixture_mass_properties( &
    species, mass_fractions, 1200.0_dp, molecular_weight, gas_constant, &
    cp, cv, gamma, enthalpy, internal_energy, entropy, ok)
  if (.not. ok) error stop "Mixture property evaluation failed"

  call assert_close(molecular_weight, 28.850334_dp, 2.0e-12_dp, "Wmix")
  call assert_close( &
    gas_constant, 288.1929414804432_dp, 2.0e-12_dp, "Rmix")
  call assert_close(cp, 1182.5018900939492_dp, 2.0e-12_dp, "cp")
  call assert_close(cv, 894.308948613506_dp, 2.0e-12_dp, "cv")
  call assert_close(gamma, 1.3222521053012428_dp, 2.0e-12_dp, "gamma")
  call assert_close( &
    enthalpy, 986641.1625680961_dp, 2.0e-12_dp, "enthalpy")
  call assert_close( &
    internal_energy, 640809.6327915643_dp, 2.0e-12_dp, "internal energy")

  pressure = mixture_pressure( &
    species, mass_fractions, 1.2_dp, 1200.0_dp, ok)
  if (.not. ok) error stop "Mixture pressure evaluation failed"
  call assert_close( &
    pressure, 414997.83573183825_dp, 2.0e-12_dp, "pressure")

  density = mixture_density( &
    species, mass_fractions, pressure, 1200.0_dp, ok)
  if (.not. ok) error stop "Mixture density evaluation failed"
  call assert_close(density, 1.2_dp, 2.0e-13_dp, "density round trip")

  sound_speed = mixture_sound_speed( &
    species, mass_fractions, 1200.0_dp, ok)
  if (.not. ok) error stop "Mixture sound-speed evaluation failed"
  call assert_close( &
    sound_speed, 676.2222033670358_dp, 2.0e-12_dp, "sound speed")

  do i = 1, size(temperatures)
    target_temperature = temperatures(i)
    call mixture_mass_properties( &
      species, mass_fractions, target_temperature, molecular_weight, &
      gas_constant, cp, cv, gamma, enthalpy, internal_energy, entropy, ok)
    if (.not. ok) error stop "Mixture energy evaluation failed"
    call temperature_from_internal_energy( &
      species, mass_fractions, internal_energy, 900.0_dp, &
      recovered_temperature, ok, iterations)
    if (.not. ok) error stop "e-to-T inversion failed"
    call assert_close( &
      recovered_temperature, target_temperature, 2.0e-11_dp, &
      "e-to-T round trip")
    if (iterations < 1 .or. iterations > 100) then
      error stop "Invalid temperature inversion iteration count"
    end if
  end do

  call mixture_mass_properties( &
    species, mass_fractions, 300.0_dp, molecular_weight, gas_constant, &
    cp, cv, gamma, enthalpy, internal_energy, entropy, ok)
  call temperature_from_internal_energy( &
    species, mass_fractions, internal_energy - 1.0e6_dp, 300.0_dp, &
    recovered_temperature, ok)
  if (ok) error stop "Out-of-range internal energy was accepted"

  invalid_mass_fractions = mass_fractions
  invalid_mass_fractions(gri_h2_index) = 0.01_dp
  if (valid_mixture_composition(species, invalid_mass_fractions)) then
    error stop "Non-unit mixture composition was accepted"
  end if

  invalid_mass_fractions = mass_fractions
  invalid_mass_fractions(gri_h2o_index) = -1.0e-4_dp
  invalid_mass_fractions(gri_n2_index) = &
    invalid_mass_fractions(gri_n2_index) + 1.0e-4_dp
  if (valid_mixture_composition(species, invalid_mass_fractions)) then
    error stop "Negative mixture composition was accepted"
  end if

  write(*, '(a)') "test_mixture_thermo: PASS"

contains

  subroutine assert_close(actual, expected, relative_tolerance, label)
    real(dp), intent(in) :: actual, expected, relative_tolerance
    character(len=*), intent(in) :: label
    real(dp) :: error

    error = abs(actual - expected) / max(1.0_dp, abs(expected))
    if (error > relative_tolerance) then
      write(*, '(a,2(1x,es24.16),1x,es12.4)') &
        trim(label), actual, expected, error
      error stop "Mixture thermodynamics reference mismatch"
    end if
  end subroutine assert_close

end program test_mixture_thermo
