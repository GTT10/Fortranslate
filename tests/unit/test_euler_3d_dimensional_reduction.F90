program test_euler_3d_dimensional_reduction
  use precision_mod, only: dp
  use state_indices_mod, only: &
    ncons, nprim, qrho, qu, qv, qw, qp
  use state_conversion_mod, only: primitive_to_conserved
  use boundary_conditions_mod, only: apply_periodic_boundaries
  use time_integrator_mod, only: advance_ssprk2
  use directional_flux_mod, only: &
    rotate_conserved_x_to_y, rotate_conserved_y_to_x, &
    rotate_conserved_x_to_z, rotate_conserved_z_to_x
  use finite_volume_3d_mod, only: &
    advance_euler_ssprk2_3d, compute_euler_cfl_timestep_3d
  implicit none

  integer, parameter :: n = 12
  real(dp), parameter :: gamma = 1.4_dp
  real(dp), parameter :: tolerance = 4.0e-13_dp
  real(dp) :: one_d(ncons, 0:n + 1)
  real(dp) :: three_d(ncons, n, n, n)
  real(dp) :: initial_one_d(ncons, 0:n + 1)
  real(dp) :: snapshot(ncons, n, n, n)
  real(dp) :: primitive(nprim), rotated(ncons)
  real(dp) :: cell_spacing, dx, dy, dz, dt, cfl_dt, x, difference
  character(len=8), parameter :: solvers(2) = &
    [character(len=8) :: "rusanov", "pelec"]
  logical :: ok
  integer :: i, j, k, direction, solver_index

  cell_spacing = 1.0_dp / real(n, dp)
  one_d = 0.0_dp
  do i = 1, n
    x = (real(i, dp) - 0.5_dp) * cell_spacing
    primitive(qrho) = 1.0_dp + 0.15_dp * sin(2.0_dp * acos(-1.0_dp) * x)
    primitive(qu) = 0.4_dp
    primitive(qv) = -0.2_dp
    primitive(qw) = 0.1_dp
    primitive(qp) = 1.0_dp
    call primitive_to_conserved(primitive, gamma, one_d(:, i), ok)
    call require(ok, "one-dimensional initial state")
  end do
  call apply_periodic_boundaries(one_d, n)
  initial_one_d = one_d
  do k = 1, n
    do j = 1, n
      three_d(:, :, j, k) = initial_one_d(:, 1:n)
    end do
  end do

  call compute_euler_cfl_timestep_3d( &
    three_d, n, n, n, cell_spacing, 1.0_dp, 1.0_dp, &
    gamma, 0.4_dp, cfl_dt, ok)
  call require(ok .and. cfl_dt > 0.0_dp, "three-dimensional CFL")
  dt = min(0.05_dp * cell_spacing, 0.5_dp * cfl_dt)
  difference = 0.0_dp
  do solver_index = 1, size(solvers)
    do direction = 1, 3
      one_d = initial_one_d
      call advance_ssprk2( &
        one_d, n, cell_spacing, dt, gamma, ok, reconstruction="pcm", &
        boundary_condition="periodic", &
        riemann_solver=trim(solvers(solver_index)))
      call require(ok, "one-dimensional reference advance")

      do k = 1, n
        do j = 1, n
          do i = 1, n
            select case (direction)
            case (1)
              three_d(:, i, j, k) = initial_one_d(:, i)
            case (2)
              call rotate_conserved_x_to_y( &
                initial_one_d(:, j), three_d(:, i, j, k))
            case (3)
              call rotate_conserved_x_to_z( &
                initial_one_d(:, k), three_d(:, i, j, k))
            end select
          end do
        end do
      end do
      dx = 1.0_dp
      dy = 1.0_dp
      dz = 1.0_dp
      select case (direction)
      case (1)
        dx = cell_spacing
      case (2)
        dy = cell_spacing
      case (3)
        dz = cell_spacing
      end select
      call advance_euler_ssprk2_3d( &
        three_d, n, n, n, dx, dy, dz, dt, gamma, &
        trim(solvers(solver_index)), ok)
      call require(ok, "three-dimensional advance")

      do k = 1, n
        do j = 1, n
          do i = 1, n
            select case (direction)
            case (1)
              rotated = three_d(:, i, j, k)
              difference = max(difference, &
                maxval(abs(rotated - one_d(:, i))))
            case (2)
              call rotate_conserved_y_to_x( &
                three_d(:, i, j, k), rotated)
              difference = max(difference, &
                maxval(abs(rotated - one_d(:, j))))
            case (3)
              call rotate_conserved_z_to_x( &
                three_d(:, i, j, k), rotated)
              difference = max(difference, &
                maxval(abs(rotated - one_d(:, k))))
            end select
          end do
        end do
      end do
    end do
  end do
  call require(difference <= tolerance, "x/y/z dimensional reductions")

  snapshot = three_d
  call advance_euler_ssprk2_3d( &
    three_d, n, n, n, dx, dy, dz, dt, gamma, "unknown", ok)
  call require(.not. ok, "unknown solver rejection")
  call require(maxval(abs(three_d - snapshot)) == 0.0_dp, &
    "rejected update is transactional")

  write(*, '(a,1x,es16.8)') "3D dimensional-reduction error:", difference
  write(*, '(a)') "test_euler_3d_dimensional_reduction: PASS"

contains

  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      write(*, '(a,1x,a)') "FAIL:", trim(label)
      error stop 1
    end if
  end subroutine require

end program test_euler_3d_dimensional_reduction
