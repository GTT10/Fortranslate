program test_ctu_dimensional_reduction
  use precision_mod, only: dp
  use state_indices_mod, only: &
    ncons, nprim, qrho, qu, qv, qw, qp
  use state_conversion_mod, only: primitive_to_conserved
  use boundary_conditions_mod, only: apply_periodic_boundaries
  use time_integrator_mod, only: compute_cfl_timestep, advance_hydro_step
  use ctu_2d_mod, only: advance_ctu_2d
  implicit none

  integer, parameter :: nx = 32, ny = 4
  real(dp), parameter :: gamma = 1.4_dp
  real(dp), parameter :: dx = 1.0_dp / real(nx, dp)
  real(dp), parameter :: dy = 0.25_dp
  real(dp), parameter :: tolerance = 2.0e-12_dp
  real(dp) :: state_1d(ncons, 0:nx + 1)
  real(dp) :: state_2d(ncons, nx, ny)
  real(dp) :: primitive(nprim), x, dt, cfl_dt, minimum_theta
  real(dp) :: maximum_difference
  logical :: ok
  integer :: i, j

  state_1d = 0.0_dp
  do i = 1, nx
    x = (real(i, dp) - 0.5_dp) * dx
    primitive(qrho) = 1.0_dp + 0.1_dp * sin(2.0_dp * acos(-1.0_dp) * x)
    primitive(qu) = 0.7_dp
    primitive(qv) = 0.2_dp
    primitive(qw) = -0.05_dp
    primitive(qp) = 1.0_dp
    call primitive_to_conserved(primitive, gamma, state_1d(:, i), ok)
    call assert_true(ok, "initial state conversion")
  end do
  call apply_periodic_boundaries(state_1d, nx)

  do j = 1, ny
    state_2d(:, :, j) = state_1d(:, 1:nx)
  end do

  call compute_cfl_timestep( &
    state_1d, nx, dx, gamma, 0.4_dp, cfl_dt, ok)
  call assert_true(ok, "1D CFL timestep")
  dt = 0.5_dp * cfl_dt

  call advance_hydro_step( &
    state_1d, nx, dx, dt, gamma, ok, reconstruction="pelec_plm", &
    limiter="mc", boundary_condition="periodic", riemann_solver="pelec", &
    plm_order=2, use_flattening=.false.)
  call assert_true(ok, "1D characteristic Godunov step")

  call advance_ctu_2d( &
    state_2d, nx, ny, dx, dy, dt, gamma, "mc", "pelec", .true., ok, &
    minimum_theta)
  call assert_true(ok, "2D CTU step")
  call assert_true(minimum_theta > 1.0_dp - tolerance, &
    "dimensionally reduced transverse limiter inactive")

  maximum_difference = 0.0_dp
  do j = 1, ny
    maximum_difference = max(maximum_difference, &
      maxval(abs(state_2d(:, :, j) - state_1d(:, 1:nx))))
  end do
  call assert_close(maximum_difference, 0.0_dp, tolerance, &
    "2D-to-1D dimensional reduction")

  write(*, '(a)') "test_ctu_dimensional_reduction: PASS"

contains

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      write(*, '(a,1x,a)') "FAIL:", trim(label)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected, tolerance_value, label)
    real(dp), intent(in) :: actual, expected, tolerance_value
    character(len=*), intent(in) :: label

    if (abs(actual - expected) > tolerance_value) then
      write(*, '(a,1x,a,2(1x,es24.16))') &
        "FAIL:", trim(label), actual, expected
      error stop 1
    end if
  end subroutine assert_close

end program test_ctu_dimensional_reduction
