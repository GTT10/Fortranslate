program test_nasa7_thermo
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: &
    nasa7_species, nasa7_mass_properties, nasa7_specific_gas_constant
  use thermo_database_mod, only: &
    load_gri30_thermo_subset, gri_h2_index, gri_o2_index
  implicit none

  type(nasa7_species), allocatable :: species(:)
  real(dp) :: cp, cv, enthalpy, internal_energy, entropy, gas_constant
  logical :: ok

  call load_gri30_thermo_subset(species, ok)
  if (.not. ok) error stop "Failed to load NASA7 database subset"

  call nasa7_mass_properties( &
    species(gri_h2_index), 300.0_dp, cp, cv, enthalpy, &
    internal_energy, entropy, ok)
  if (.not. ok) error stop "H2 NASA7 evaluation failed"
  call assert_close(cp, 14311.757145676056_dp, 2.0e-12_dp, "H2 cp")
  call assert_close(cv, 10187.274230942421_dp, 2.0e-12_dp, "H2 cv")
  call assert_close(enthalpy, 26470.080162950748_dp, 2.0e-12_dp, "H2 h")
  call assert_close( &
    internal_energy, -1210874.7942571398_dp, 2.0e-12_dp, "H2 u")
  gas_constant = nasa7_specific_gas_constant(species(gri_h2_index))
  call assert_close(cp - cv, gas_constant, 2.0e-13_dp, "H2 cp-cv")
  call assert_close( &
    enthalpy - internal_energy, gas_constant * 300.0_dp, &
    2.0e-13_dp, "H2 h-u")

  call nasa7_mass_properties( &
    species(gri_o2_index), 1500.0_dp, cp, cv, enthalpy, &
    internal_energy, entropy, ok)
  if (.not. ok) error stop "O2 NASA7 evaluation failed"
  call assert_close(cp, 1143.0200573618715_dp, 2.0e-12_dp, "O2 cp")
  call assert_close(cv, 883.1833566683067_dp, 2.0e-12_dp, "O2 cv")
  call assert_close(enthalpy, 1268862.4251066474_dp, 2.0e-12_dp, "O2 h")
  call assert_close( &
    internal_energy, 879107.3740663002_dp, 2.0e-12_dp, "O2 u")
  gas_constant = nasa7_specific_gas_constant(species(gri_o2_index))
  call assert_close(cp - cv, gas_constant, 2.0e-13_dp, "O2 cp-cv")
  call assert_close( &
    enthalpy - internal_energy, gas_constant * 1500.0_dp, &
    2.0e-13_dp, "O2 h-u")

  call nasa7_mass_properties( &
    species(gri_h2_index), 199.0_dp, cp, cv, enthalpy, &
    internal_energy, entropy, ok)
  if (ok) error stop "NASA7 accepted a temperature below its valid range"

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
