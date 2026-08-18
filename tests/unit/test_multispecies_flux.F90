program test_multispecies_flux
  use precision_mod, only: dp
  use state_indices_mod, only: irho, ncons, nbase, nprim, qrho, qu, qv, qw, qp
  use state_conversion_mod, only: primitive_to_conserved
  use multispecies_state_mod, only: &
    species_component, multispecies_state_from_base
  use multispecies_flux_mod, only: &
    compute_multispecies_flux_x, compute_multispecies_flux_y
  implicit none

  integer, parameter :: nspecies = 2
  real(dp), parameter :: gamma = 1.4_dp
  real(dp), parameter :: tolerance = 3.0e-12_dp
  real(dp) :: primitive(nprim), base_state(ncons)
  real(dp) :: left_state(nbase + nspecies), right_state(nbase + nspecies)
  real(dp) :: flux(nbase + nspecies), y(nspecies)
  logical :: ok

  primitive(qrho) = 1.0_dp
  primitive(qu) = 2.0_dp
  primitive(qv) = -0.5_dp
  primitive(qw) = 0.1_dp
  primitive(qp) = 1.0_dp
  call primitive_to_conserved(primitive, gamma, base_state, ok)
  call assert_true(ok, "base state")
  y = [0.25_dp, 0.75_dp]
  call multispecies_state_from_base(base_state, y, nspecies, gamma, left_state, ok)
  call assert_true(ok, "left assembly")
  right_state = left_state

  call compute_multispecies_flux_x( &
    left_state, right_state, nspecies, gamma, "pelec", flux, ok)
  call assert_true(ok, "x flux")
  call assert_close(flux(irho), 2.0_dp, tolerance, "x mass flux")
  call assert_close(flux(species_component(1)), 0.5_dp, tolerance, &
    "x species 1 flux")
  call assert_close(flux(species_component(2)), 1.5_dp, tolerance, &
    "x species 2 flux")
  call assert_close(sum(flux(nbase + 1:nbase + nspecies)), &
    flux(irho), tolerance, "x species-flux closure")

  call compute_multispecies_flux_y( &
    left_state, right_state, nspecies, gamma, "rusanov", flux, ok)
  call assert_true(ok, "y flux")
  call assert_close(flux(irho), -0.5_dp, tolerance, "y mass flux")
  call assert_close(flux(species_component(1)), -0.125_dp, tolerance, &
    "y species 1 flux")
  call assert_close(sum(flux(nbase + 1:nbase + nspecies)), &
    flux(irho), tolerance, "y species-flux closure")

  primitive(qrho) = 0.5_dp
  call primitive_to_conserved(primitive, gamma, base_state, ok)
  y = [0.0_dp, 1.0_dp]
  call multispecies_state_from_base(base_state, y, nspecies, gamma, right_state, ok)
  call assert_true(ok, "contact right state")
  call compute_multispecies_flux_x( &
    left_state, right_state, nspecies, gamma, "pelec", flux, ok)
  call assert_true(ok, "contact flux")
  call assert_true(flux(irho) > 0.0_dp, "positive contact mass flux")
  call assert_close( &
    flux(species_component(1)), 0.25_dp * flux(irho), tolerance, &
    "upwind contact species")

  call compute_multispecies_flux_x( &
    left_state, right_state, nspecies, gamma, "unknown", flux, ok)
  call assert_true(.not. ok, "unknown solver rejection")

  write(*, '(a)') "test_multispecies_flux: PASS"

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

end program test_multispecies_flux
