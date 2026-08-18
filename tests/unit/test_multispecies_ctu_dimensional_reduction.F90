program test_multispecies_ctu_dimensional_reduction
  use precision_mod, only: dp
  use state_indices_mod, only: ncons, nprim, qrho, qu, qv, qw, qp
  use state_conversion_mod, only: primitive_to_conserved
  use multispecies_state_mod, only: &
    multispecies_nvar, multispecies_state_from_base
  use boundary_conditions_multispecies_mod, only: &
    apply_multispecies_boundary_conditions
  use time_integrator_multispecies_mod, only: advance_multispecies_hydro_step
  use ctu_multispecies_2d_mod, only: advance_ctu_multispecies_2d
  implicit none

  integer, parameter :: nx = 32, ny = 5, nspecies = 2
  real(dp), parameter :: gamma = 1.4_dp
  real(dp), parameter :: dx = 1.0_dp / real(nx, dp)
  real(dp), parameter :: dy = 1.0_dp / real(ny, dp)
  real(dp), parameter :: dt = 0.12_dp * dx
  real(dp), parameter :: tolerance = 5.0e-12_dp
  real(dp), allocatable :: one_d(:, :), two_d(:, :, :)
  real(dp) :: primitive(nprim), base_state(ncons), mass_fractions(nspecies)
  real(dp) :: x, difference, minimum_theta
  logical :: ok
  integer :: i, j, nvar, fallback_count

  nvar = multispecies_nvar(nspecies)
  allocate(one_d(nvar, 0:nx + 1), two_d(nvar, nx, ny))

  do i = 1, nx
    x = (real(i, dp) - 0.5_dp) * dx
    primitive(qrho) = 1.0_dp + 0.1_dp * sin(2.0_dp * acos(-1.0_dp) * x)
    primitive(qu) = 0.7_dp
    primitive(qv) = 0.0_dp
    primitive(qw) = 0.0_dp
    primitive(qp) = 1.0_dp
    call primitive_to_conserved(primitive, gamma, base_state, ok)
    call assert_true(ok, "base conversion")

    mass_fractions(1) = 0.5_dp + 0.2_dp * &
      cos(2.0_dp * acos(-1.0_dp) * x)
    mass_fractions(2) = 1.0_dp - mass_fractions(1)
    call multispecies_state_from_base( &
      base_state, mass_fractions, nspecies, gamma, one_d(:, i), ok)
    call assert_true(ok, "multispecies assembly")
    do j = 1, ny
      two_d(:, i, j) = one_d(:, i)
    end do
  end do

  call apply_multispecies_boundary_conditions(one_d, nx, "periodic", ok)
  call assert_true(ok, "periodic fill")

  call advance_multispecies_hydro_step( &
    one_d, nx, nspecies, dx, dt, gamma, ok, &
    reconstruction="pelec_plm", limiter="mc", &
    boundary_condition="periodic", riemann_solver="pelec", &
    plm_order=2, use_flattening=.false.)
  call assert_true(ok, "one-dimensional step")

  call advance_ctu_multispecies_2d( &
    two_d, nx, ny, nspecies, dx, dy, dt, gamma, "mc", "pelec", &
    .true., ok, minimum_theta, fallback_count)
  call assert_true(ok, "two-dimensional step")
  call assert_true(fallback_count == 0, "no species face fallback")
  call assert_close(minimum_theta, 1.0_dp, tolerance, "transverse theta")

  difference = 0.0_dp
  do j = 1, ny
    difference = max(difference, maxval(abs(two_d(:, :, j) - one_d(:, 1:nx))))
  end do
  call assert_close(difference, 0.0_dp, tolerance, "dimensional reduction")

  write(*, '(a)') "test_multispecies_ctu_dimensional_reduction: PASS"

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

end program test_multispecies_ctu_dimensional_reduction
