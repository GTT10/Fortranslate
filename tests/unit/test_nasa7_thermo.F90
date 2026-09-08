program test_nasa7_thermo
  use, intrinsic :: ieee_arithmetic, only: &
    ieee_positive_inf, ieee_quiet_nan, ieee_value
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: &
    nasa7_species, valid_nasa7_species, nasa7_mass_properties, &
    nasa7_specific_gas_constant, nasa7_molar_properties, &
    nasa7_dimensionless_properties, universal_gas_constant
  use thermo_database_mod, only: &
    load_gri30_thermo_subset, gri_h2_index, gri_o2_index
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(nasa7_species) :: invalid_species
  real(dp) :: cp, cv, enthalpy, internal_energy, entropy, gas_constant
  real(dp) :: nan_value, positive_inf
  logical :: ok

  call load_gri30_thermo_subset(species, ok)
  if (.not. ok) error stop "Failed to load NASA7 database subset"

  call nasa7_mass_properties( &
    species(gri_h2_index), 300.0_dp, cp, cv, enthalpy, &
    internal_energy, entropy, ok)
  if (.not. ok) error stop "H2 NASA7 evaluation failed"
  call assert_close(cp, 14310.905255369766_dp, 2.0e-12_dp, "H2 cp")
  call assert_close(cv, 10186.667845571532_dp, 2.0e-12_dp, "H2 cv")
  call assert_close(enthalpy, 26468.504562941045_dp, 2.0e-12_dp, "H2 h")
  call assert_close( &
    internal_energy, -1210802.7183765292_dp, 2.0e-12_dp, "H2 u")
  gas_constant = nasa7_specific_gas_constant(species(gri_h2_index))
  call assert_close(cp - cv, gas_constant, 2.0e-13_dp, "H2 cp-cv")
  call assert_close( &
    enthalpy - internal_energy, gas_constant * 300.0_dp, &
    2.0e-13_dp, "H2 h-u")

  call nasa7_mass_properties( &
    species(gri_o2_index), 1500.0_dp, cp, cv, enthalpy, &
    internal_energy, entropy, ok)
  if (.not. ok) error stop "O2 NASA7 evaluation failed"
  call assert_close(cp, 1143.0486346493860_dp, 2.0e-12_dp, "O2 cp")
  call assert_close(cv, 883.20543763228386_dp, 2.0e-12_dp, "O2 cv")
  call assert_close(enthalpy, 1268894.1486499966_dp, 2.0e-12_dp, "O2 h")
  call assert_close( &
    internal_energy, 879129.35312434332_dp, 2.0e-12_dp, "O2 u")
  gas_constant = nasa7_specific_gas_constant(species(gri_o2_index))
  call assert_close(cp - cv, gas_constant, 2.0e-13_dp, "O2 cp-cv")
  call assert_close( &
    enthalpy - internal_energy, gas_constant * 1500.0_dp, &
    2.0e-13_dp, "O2 h-u")

  call nasa7_mass_properties( &
    species(gri_h2_index), 199.0_dp, cp, cv, enthalpy, &
    internal_energy, entropy, ok)
  if (ok) error stop "NASA7 accepted a temperature below its valid range"

  nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
  positive_inf = ieee_value(0.0_dp, ieee_positive_inf)

  call nasa7_mass_properties( &
    species(gri_h2_index), nan_value, cp, cv, enthalpy, &
    internal_energy, entropy, ok)
  if (ok .or. any([cp, cv, enthalpy, internal_energy, entropy] /= 0.0_dp)) &
    error stop "NASA7 accepted a NaN temperature"

  call nasa7_mass_properties( &
    species(gri_h2_index), positive_inf, cp, cv, enthalpy, &
    internal_energy, entropy, ok)
  if (ok .or. any([cp, cv, enthalpy, internal_energy, entropy] /= 0.0_dp)) &
    error stop "NASA7 accepted an infinite temperature"

  invalid_species = species(gri_h2_index)
  invalid_species%temperature_mid = nan_value
  if (valid_nasa7_species(invalid_species)) &
    error stop "NASA7 accepted a NaN temperature bound"

  invalid_species = species(gri_h2_index)
  invalid_species%low_coefficients(1) = positive_inf
  if (valid_nasa7_species(invalid_species)) &
    error stop "NASA7 accepted an infinite coefficient"

  invalid_species = species(gri_h2_index)
  invalid_species%low_coefficients(1) = huge(1.0_dp)
  call nasa7_mass_properties( &
    invalid_species, 300.0_dp, cp, cv, enthalpy, internal_energy, entropy, ok)
  if (ok .or. any([cp, cv, enthalpy, internal_energy, entropy] /= 0.0_dp)) &
    error stop "NASA7 failed to reject an overflowing polynomial"

  invalid_species = species(gri_h2_index)
  invalid_species%temperature_min = 1.0_dp / huge(1.0_dp)
  invalid_species%temperature_mid = 1.0_dp
  invalid_species%temperature_max = 2.0_dp
  call nasa7_mass_properties( &
    invalid_species, invalid_species%temperature_min, cp, cv, enthalpy, &
    internal_energy, entropy, ok)
  if (ok .or. any([cp, cv, enthalpy, internal_energy, entropy] /= 0.0_dp)) &
    error stop "NASA7 failed to reject an overflowing reciprocal term"

  invalid_species = species(gri_h2_index)
  invalid_species%molecular_weight = tiny(1.0_dp)
  gas_constant = nasa7_specific_gas_constant(invalid_species)
  if (gas_constant /= -huge(1.0_dp)) &
    error stop "NASA7 failed to reject an overflowing gas constant"
  call nasa7_mass_properties( &
    invalid_species, 300.0_dp, cp, cv, enthalpy, internal_energy, entropy, ok)
  if (ok .or. any([cp, cv, enthalpy, internal_energy, entropy] /= 0.0_dp)) &
    error stop "NASA7 failed to reject overflowing mass properties"

  call check_molar_conversion_rejection()

  write(*, '(a)') "test_nasa7_thermo: PASS"

contains

  subroutine check_molar_conversion_rejection()
    type(nasa7_species) :: synthetic_species
    real(dp) :: cp_over_r, h_over_rt, s_over_r, oversized_coefficient
    real(dp) :: local_cp, local_cv, local_h, local_u, local_s
    logical :: local_ok
    integer :: coefficient_index, test_case
    integer, parameter :: conversion_coefficients(3) = [1, 6, 7]

    synthetic_species = species(gri_h2_index)
    synthetic_species%temperature_min = 0.5_dp
    synthetic_species%temperature_mid = 1.0_dp
    synthetic_species%temperature_max = 2.0_dp
    ! These finite dimensionless values pass polynomial evaluation, but their
    ! dimensional conversion would overflow. Test cp, enthalpy, and entropy.
    oversized_coefficient = 2.0_dp * (huge(1.0_dp) / universal_gas_constant)
    do test_case = 1, size(conversion_coefficients)
      synthetic_species%low_coefficients = 0.0_dp
      synthetic_species%low_coefficients(1) = 3.5_dp
      coefficient_index = conversion_coefficients(test_case)
      synthetic_species%low_coefficients(coefficient_index) = &
        oversized_coefficient
      call nasa7_dimensionless_properties( &
        synthetic_species, 1.0_dp, cp_over_r, h_over_rt, s_over_r, local_ok)
      if (.not. local_ok) &
        error stop "Molar overflow fixture failed before dimensional conversion"

      call nasa7_molar_properties( &
        synthetic_species, 1.0_dp, local_cp, local_cv, local_h, local_u, &
        local_s, local_ok)
      if (local_ok .or. any([local_cp, local_cv, local_h, local_u, local_s] /= 0.0_dp)) &
        error stop "NASA7 molar conversion reported success after overflow rejection"
      call nasa7_mass_properties( &
        synthetic_species, 1.0_dp, local_cp, local_cv, local_h, local_u, &
        local_s, local_ok)
      if (local_ok .or. any([local_cp, local_cv, local_h, local_u, local_s] /= 0.0_dp)) &
        error stop "NASA7 mass conversion did not propagate molar overflow failure"
    end do
  end subroutine check_molar_conversion_rejection

  subroutine assert_close(actual, expected, relative_tolerance, label)
    real(dp), intent(in) :: actual, expected, relative_tolerance
    character(len=*), intent(in) :: label
    real(dp) :: error

    error = abs(actual - expected) / max(1.0_dp, abs(expected))
    if (error > relative_tolerance) then
      write(*, '(a,2(1x,es24.16),1x,es12.4)') &
        trim(label), actual, expected, error
      error stop "NASA7 reference mismatch"
    end if
  end subroutine assert_close

end program test_nasa7_thermo
