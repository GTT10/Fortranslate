program test_isomerization_reactor
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_toy_isomerization_thermo
  use isomerization_reactor_mod, only: &
    isomerization_reaction, advance_isomerization_rk4, reactor_internal_energy
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(isomerization_reaction) :: reaction
  real(dp) :: mass_fractions(2), temperature, target_energy, final_energy
  real(dp) :: expected_reactant, relative_energy_error
  logical :: ok
  integer :: step

  call load_toy_isomerization_thermo(species, ok)
  if (.not. ok) error stop "Failed to load toy reactor thermodynamics"

  reaction%reactant = 1
  reaction%product = 2
  reaction%pre_exponential = 3.0_dp
  reaction%temperature_exponent = 0.0_dp
  reaction%activation_temperature = 0.0_dp

  mass_fractions = [1.0_dp, 0.0_dp]
  temperature = 800.0_dp
  target_energy = reactor_internal_energy( &
    species, mass_fractions, temperature, ok)
  if (.not. ok) error stop "Initial isothermal energy evaluation failed"

  do step = 1, 1000
    call advance_isomerization_rk4( &
      species, reaction, 1.0e-3_dp, .false., target_energy, &
      mass_fractions, temperature, ok)
    if (.not. ok) error stop "Isothermal reactor integration failed"
  end do
  expected_reactant = exp(-3.0_dp)
  call assert_close( &
    mass_fractions(1), expected_reactant, 2.0e-11_dp, &
    "isothermal analytical solution")
  call assert_close(temperature, 800.0_dp, 1.0e-14_dp, &
    "isothermal temperature")
  call assert_close(sum(mass_fractions), 1.0_dp, 1.0e-14_dp, &
    "isothermal closure")

  reaction%pre_exponential = 2.0e4_dp
  reaction%activation_temperature = 5000.0_dp
  mass_fractions = [0.95_dp, 0.05_dp]
  temperature = 700.0_dp
  target_energy = reactor_internal_energy( &
    species, mass_fractions, temperature, ok)
  if (.not. ok) error stop "Initial adiabatic energy evaluation failed"

  do step = 1, 4000
    call advance_isomerization_rk4( &
      species, reaction, 5.0e-5_dp, .true., target_energy, &
      mass_fractions, temperature, ok)
    if (.not. ok) error stop "Adiabatic reactor integration failed"
  end do

  final_energy = reactor_internal_energy( &
    species, mass_fractions, temperature, ok)
  if (.not. ok) error stop "Final adiabatic energy evaluation failed"
  relative_energy_error = abs(final_energy - target_energy) / &
    max(1.0_dp, abs(target_energy))

  if (mass_fractions(1) > 1.0e-8_dp) then
    error stop "Adiabatic reaction did not approach completion"
  end if
  call assert_close(sum(mass_fractions), 1.0_dp, 2.0e-13_dp, &
    "adiabatic closure")
  call assert_close(temperature, 1840.0_dp, 2.0e-11_dp, &
    "adiabatic final temperature")
  if (relative_energy_error > 1.0e-11_dp) then
    write(*, '(a,es24.16)') "relative energy error=", relative_energy_error
    error stop "Adiabatic internal energy drift"
  end if

  write(*, '(a)') "test_isomerization_reactor: PASS"

contains

  subroutine assert_close(actual, expected, relative_tolerance, label)
    real(dp), intent(in) :: actual, expected, relative_tolerance
    character(len=*), intent(in) :: label
    real(dp) :: error

    error = abs(actual - expected) / max(1.0_dp, abs(expected))
    if (error > relative_tolerance) then
      write(*, '(a,2(1x,es24.16),1x,es12.4)') &
        trim(label), actual, expected, error
      error stop "Reactor reference mismatch"
    end if
  end subroutine assert_close

end program test_isomerization_reactor
