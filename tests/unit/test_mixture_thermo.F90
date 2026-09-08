program test_mixture_thermo
  use, intrinsic :: ieee_arithmetic, only: &
    ieee_is_finite, ieee_positive_inf, ieee_quiet_nan, ieee_value
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: &
    load_gri30_thermo_subset, gri_h2_index, gri_o2_index, &
    gri_h2o_index, gri_n2_index
  use mixture_thermo_mod, only: &
    valid_mixture_composition, mixture_mass_properties, mixture_pressure, &
    mixture_density, mixture_sound_speed, temperature_from_internal_energy, &
    mixture_temperature_bounds, mixture_molecular_weight, &
    mixture_specific_gas_constant, mass_fractions_from_mole_fractions, &
    mole_fractions_from_mass_fractions
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

  call assert_close( &
    molecular_weight, 28.85067068107281_dp, 2.0e-12_dp, "Wmix")
  call assert_close( &
    gas_constant, 288.18957833128849_dp, 2.0e-12_dp, "Rmix")
  call assert_close(cp, 1182.4886092588142_dp, 2.0e-12_dp, "cp")
  call assert_close(cv, 894.29903092752579_dp, 2.0e-12_dp, "cv")
  call assert_close(gamma, 1.3222519183906434_dp, 2.0e-12_dp, "gamma")
  call assert_close( &
    enthalpy, 986630.08716612426_dp, 2.0e-12_dp, "enthalpy")
  call assert_close( &
    internal_energy, 640802.59316857834_dp, 2.0e-12_dp, &
    "internal energy")

  pressure = mixture_pressure( &
    species, mass_fractions, 1.2_dp, 1200.0_dp, ok)
  if (.not. ok) error stop "Mixture pressure evaluation failed"
  call assert_close( &
    pressure, 414992.99279705540_dp, 2.0e-12_dp, "pressure")

  density = mixture_density( &
    species, mass_fractions, pressure, 1200.0_dp, ok)
  if (.not. ok) error stop "Mixture density evaluation failed"
  call assert_close(density, 1.2_dp, 2.0e-13_dp, "density round trip")

  sound_speed = mixture_sound_speed( &
    species, mass_fractions, 1200.0_dp, ok)
  if (.not. ok) error stop "Mixture sound-speed evaluation failed"
  call assert_close( &
    sound_speed, 676.21820987790932_dp, 2.0e-12_dp, "sound speed")

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

  call check_nonfinite_contract(species, mass_fractions, ok)
  if (.not. ok) error stop "Mixture thermodynamics finite contract failed"

  write(*, '(a)') "test_mixture_thermo: PASS"

contains

  subroutine check_nonfinite_contract(active_species, valid_mass_fractions, ok_out)
    type(nasa7_species), intent(in) :: active_species(:)
    real(dp), intent(in) :: valid_mass_fractions(:)
    logical, intent(out) :: ok_out

    type(nasa7_species), allocatable :: bad_species(:)
    real(dp) :: bad_values(2), bad_value
    real(dp) :: bad_mass_fractions(4), bad_mole_fractions(4)
    real(dp) :: output_mass_fractions(4), output_mole_fractions(4)
    real(dp) :: lower, upper, molecular_weight, gas_constant
    real(dp) :: cp, cv, gamma, enthalpy, internal_energy, entropy
    real(dp) :: pressure, density, sound_speed, recovered_temperature
    logical :: local_ok
    integer :: bad_index

    ok_out = .false.
    if (size(active_species) /= 4 .or. size(valid_mass_fractions) /= 4) return
    allocate(bad_species, source=active_species)
    bad_values(1) = ieee_value(0.0_dp, ieee_quiet_nan)
    bad_values(2) = ieee_value(0.0_dp, ieee_positive_inf)

    do bad_index = 1, 2
      bad_value = bad_values(bad_index)

      bad_mass_fractions = valid_mass_fractions
      bad_mass_fractions(1) = bad_value
      if (valid_mixture_composition(active_species, bad_mass_fractions)) return

      bad_species = active_species
      bad_species(1)%molecular_weight = bad_value
      call mixture_temperature_bounds(bad_species, lower, upper, local_ok)
      if (local_ok .or. lower /= 0.0_dp .or. upper /= 0.0_dp) return
      molecular_weight = mixture_molecular_weight( &
        bad_species, valid_mass_fractions, local_ok)
      if (local_ok .or. molecular_weight /= 0.0_dp) return
      gas_constant = mixture_specific_gas_constant( &
        bad_species, valid_mass_fractions, local_ok)
      if (local_ok .or. gas_constant /= 0.0_dp) return

      call mixture_mass_properties( &
        active_species, valid_mass_fractions, bad_value, molecular_weight, &
        gas_constant, cp, cv, gamma, enthalpy, internal_energy, entropy, local_ok)
      if (local_ok .or. .not. all([molecular_weight, gas_constant, cp, cv, &
          gamma, enthalpy, internal_energy, entropy] == 0.0_dp)) return

      pressure = mixture_pressure( &
        active_species, valid_mass_fractions, bad_value, 1200.0_dp, local_ok)
      if (local_ok .or. pressure /= 0.0_dp) return
      pressure = mixture_pressure( &
        active_species, valid_mass_fractions, 1.2_dp, bad_value, local_ok)
      if (local_ok .or. pressure /= 0.0_dp) return
      pressure = mixture_pressure( &
        active_species, valid_mass_fractions, huge(1.0_dp), &
        huge(1.0_dp), local_ok)
      if (local_ok .or. pressure /= 0.0_dp) return

      density = mixture_density( &
        active_species, valid_mass_fractions, bad_value, 1200.0_dp, local_ok)
      if (local_ok .or. density /= 0.0_dp) return
      density = mixture_density( &
        active_species, valid_mass_fractions, 101325.0_dp, bad_value, local_ok)
      if (local_ok .or. density /= 0.0_dp) return
      density = mixture_density( &
        active_species, valid_mass_fractions, huge(1.0_dp), &
        tiny(1.0_dp), local_ok)
      if (local_ok .or. density /= 0.0_dp) return

      sound_speed = mixture_sound_speed( &
        active_species, valid_mass_fractions, bad_value, local_ok)
      if (local_ok .or. sound_speed /= 0.0_dp) return

      call mixture_mass_properties( &
        active_species, valid_mass_fractions, 1200.0_dp, molecular_weight, &
        gas_constant, cp, cv, gamma, enthalpy, internal_energy, entropy, local_ok)
      if (.not. local_ok .or. .not. ieee_is_finite(internal_energy)) return
      recovered_temperature = 77.0_dp
      call temperature_from_internal_energy( &
        active_species, valid_mass_fractions, bad_value, 900.0_dp, &
        recovered_temperature, local_ok)
      if (local_ok .or. recovered_temperature /= 0.0_dp) return
      recovered_temperature = 77.0_dp
      call temperature_from_internal_energy( &
        active_species, valid_mass_fractions, internal_energy, bad_value, &
        recovered_temperature, local_ok)
      if (local_ok .or. recovered_temperature /= 0.0_dp) return

      bad_mole_fractions = valid_mass_fractions
      bad_mole_fractions(1) = bad_value
      output_mass_fractions = 1.0_dp
      call mass_fractions_from_mole_fractions( &
        active_species, bad_mole_fractions, output_mass_fractions, local_ok)
      if (local_ok .or. .not. all(output_mass_fractions == 0.0_dp)) return
      output_mole_fractions = 1.0_dp
      call mole_fractions_from_mass_fractions( &
        active_species, bad_mass_fractions, output_mole_fractions, local_ok)
      if (local_ok .or. .not. all(output_mole_fractions == 0.0_dp)) return
    end do

    bad_species = active_species
    bad_species(1)%molecular_weight = 0.125_dp * tiny(1.0_dp)
    bad_mass_fractions = 0.0_dp
    bad_mass_fractions(1) = 1.0_dp
    molecular_weight = mixture_molecular_weight( &
      bad_species, bad_mass_fractions, local_ok)
    if (local_ok .or. molecular_weight /= 0.0_dp) return
    gas_constant = mixture_specific_gas_constant( &
      bad_species, bad_mass_fractions, local_ok)
    if (local_ok .or. gas_constant /= 0.0_dp) return
    output_mole_fractions = 1.0_dp
    call mole_fractions_from_mass_fractions( &
      bad_species, bad_mass_fractions, output_mole_fractions, local_ok)
    if (local_ok .or. .not. all(output_mole_fractions == 0.0_dp)) return

    bad_species = active_species
    bad_species(1)%molecular_weight = huge(1.0_dp)
    molecular_weight = mixture_molecular_weight( &
      bad_species, bad_mass_fractions, local_ok)
    if (local_ok .or. molecular_weight /= 0.0_dp) return
    gas_constant = mixture_specific_gas_constant( &
      bad_species, bad_mass_fractions, local_ok)
    if (local_ok .or. gas_constant /= 0.0_dp) return

    ok_out = .true.
  end subroutine check_nonfinite_contract

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
