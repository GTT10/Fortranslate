program test_nasa7_thermo
  use, intrinsic :: ieee_arithmetic, only: &
    ieee_positive_inf, ieee_quiet_nan, ieee_value
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: &
    nasa7_species, valid_nasa7_species, nasa7_mass_properties, &
    nasa7_specific_gas_constant
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

  write(*, '(a)') "test_nasa7_thermo: PASS"

contains

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
