program test_multispecies_state
  use precision_mod, only: dp
  use state_indices_mod, only: &
    iei, item, ncons, nbase, nprim, qrho, qu, qv, qw, qp
  use state_conversion_mod, only: primitive_to_conserved
  use multispecies_state_mod, only: &
    multispecies_nvar, species_component, normalize_mass_fractions, &
    multispecies_state_from_base, mass_fractions_from_state, &
    multispecies_state_is_physical, species_closure_error, &
    thermodynamic_layout_error
  implicit none

  integer, parameter :: nspecies = 3
  real(dp), parameter :: gamma = 1.4_dp
  real(dp), parameter :: tolerance = 2.0e-13_dp
  real(dp) :: primitive(nprim), base_state(ncons)
  real(dp) :: state(nbase + nspecies), altered(nbase + nspecies)
  real(dp) :: mass_fractions(nspecies), recovered(nspecies)
  logical :: ok

  primitive(qrho) = 1.7_dp
  primitive(qu) = 0.4_dp
  primitive(qv) = -0.2_dp
  primitive(qw) = 0.1_dp
  primitive(qp) = 2.3_dp
  call primitive_to_conserved(primitive, gamma, base_state, ok)
  call assert_true(ok, "base conversion")

  mass_fractions = [0.2_dp, 0.3_dp, 0.5_dp]
  call multispecies_state_from_base( &
    base_state, mass_fractions, nspecies, gamma, state, ok)
  call assert_true(ok, "state assembly")
  call assert_true(multispecies_nvar(nspecies) == nbase + nspecies, &
    "state size")
  call assert_true(species_component(2) == nbase + 2, "species index")

  call mass_fractions_from_state(state, nspecies, recovered, ok)
  call assert_true(ok, "mass-fraction extraction")
  call assert_close(maxval(abs(recovered - mass_fractions)), &
    0.0_dp, tolerance, "mass-fraction round trip")
  call assert_true(multispecies_state_is_physical(state, gamma, nspecies), &
    "physical state")
  call assert_close(species_closure_error(state, nspecies), &
    0.0_dp, tolerance, "species closure")
  call assert_close(state(iei), primitive(qp) / (gamma - 1.0_dp), &
    tolerance, "internal-energy density slot")
  call assert_close(state(item), primitive(qp) / primitive(qrho), &
    tolerance, "temperature slot")
  call assert_close(thermodynamic_layout_error(state, gamma), &
    0.0_dp, tolerance, "thermodynamic layout")

  altered = state
  altered(item) = 1.1_dp * altered(item)
  call assert_true(.not. multispecies_state_is_physical( &
    altered, gamma, nspecies), "thermodynamic-layout rejection")

  altered = state
  altered(species_component(1)) = altered(species_component(1)) + 1.0e-5_dp
  call assert_true(.not. multispecies_state_is_physical( &
    altered, gamma, nspecies), "closure rejection")

  recovered = [0.2_dp, 0.2_dp, 0.2_dp]
  call normalize_mass_fractions(recovered, nspecies, ok)
  call assert_true(ok, "normalization")
  call assert_close(sum(recovered), 1.0_dp, tolerance, "normalized sum")

  recovered = [0.5_dp, -1.0e-4_dp, 0.5_dp]
  call normalize_mass_fractions(recovered, nspecies, ok)
  call assert_true(.not. ok, "negative mass-fraction rejection")

  write(*, '(a)') "test_multispecies_state: PASS"

contains

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*, '(a,1x,a)') "FAIL:", trim(label)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected, tol, label)
    real(dp), intent(in) :: actual, expected, tol
    character(len=*), intent(in) :: label
    if (abs(actual - expected) > tol) then
      write(*, '(a,1x,a,2(1x,es24.16))') &
        "FAIL:", trim(label), actual, expected
      error stop 1
    end if
  end subroutine assert_close

end program test_multispecies_state
