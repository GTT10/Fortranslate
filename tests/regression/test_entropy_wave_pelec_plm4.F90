program test_entropy_wave_pelec_plm4
  use precision_mod, only: dp
  use state_indices_mod, only: ncons, nprim, irho, qrho, qu, qv, qw, qp
  use mesh_mod, only: uniform_cell_centers
  use state_conversion_mod, only: primitive_to_conserved
  use boundary_conditions_mod, only: apply_periodic_boundaries
  use time_integrator_mod, only: compute_cfl_timestep, advance_hydro_step
  implicit none

  integer, parameter :: number_of_resolutions = 3
  integer, parameter :: resolutions(number_of_resolutions) = [40, 80, 160]
  real(dp), parameter :: minimum_order = 1.9_dp
  real(dp) :: errors(number_of_resolutions), order
  logical :: ok
  integer :: case_index

  do case_index = 1, number_of_resolutions
    call run_entropy_wave(resolutions(case_index), errors(case_index), ok)
    if (.not. ok) error stop "PeleC fourth-order-slope entropy wave failed"
    write(*, '(a,i0,a,es24.16)') &
      "nx=", resolutions(case_index), ", L1=", errors(case_index)
  end do

  do case_index = 1, number_of_resolutions - 1
    order = log(errors(case_index) / errors(case_index + 1)) / log(2.0_dp)
    write(*, '(a,i0,a,es24.16)') &
      "refinement pair ", case_index, ", order=", order
    if (order < minimum_order) error stop "Observed order below threshold"
  end do

  write(*, '(a)') "test_entropy_wave_pelec_plm4: PASS"

contains

  subroutine run_entropy_wave(nx, density_l1_error, ok)
    integer, intent(in) :: nx
    real(dp), intent(out) :: density_l1_error
    logical, intent(out) :: ok

    real(dp), parameter :: x_min = 0.0_dp, x_max = 1.0_dp
    real(dp), parameter :: final_time = 0.1_dp
    real(dp), parameter :: gamma = 1.4_dp, cfl = 0.45_dp
    real(dp), parameter :: advection_velocity = 1.0_dp
    real(dp), parameter :: density_amplitude = 0.2_dp
    real(dp), parameter :: background_pressure = 1.0_dp
    real(dp), allocatable :: x(:), conserved(:, :)
    real(dp) :: primitive(nprim), dx, time, dt, phase, exact_density
    logical :: cell_ok, step_ok
    integer :: i

    allocate(x(nx))
    allocate(conserved(ncons, 0:nx + 1))
    call uniform_cell_centers(nx, x_min, x_max, x, dx)

    conserved = 0.0_dp
    do i = 1, nx
      primitive(qrho) = 1.0_dp + density_amplitude * &
        sin(2.0_dp * acos(-1.0_dp) * x(i))
      primitive(qu) = advection_velocity
      primitive(qv) = 0.0_dp
      primitive(qw) = 0.0_dp
      primitive(qp) = background_pressure
      call primitive_to_conserved( &
        primitive, gamma, conserved(:, i), cell_ok)
      if (.not. cell_ok) then
        ok = .false.
        density_l1_error = huge(1.0_dp)
        return
      end if
    end do
    call apply_periodic_boundaries(conserved, nx)

    time = 0.0_dp
    do while (time < final_time)
      call compute_cfl_timestep( &
        conserved, nx, dx, gamma, cfl, dt, step_ok)
      if (.not. step_ok) then
        ok = .false.
        density_l1_error = huge(1.0_dp)
        return
      end if

      dt = min(dt, final_time - time)
      call advance_hydro_step( &
        conserved, nx, dx, dt, gamma, step_ok, &
        reconstruction="pelec_plm", limiter="mc", &
        boundary_condition="periodic", riemann_solver="pelec", &
        plm_order=4, use_flattening=.false.)
      if (.not. step_ok) then
        ok = .false.
        density_l1_error = huge(1.0_dp)
        return
      end if
      time = time + dt
    end do

    density_l1_error = 0.0_dp
    do i = 1, nx
      phase = x_min + modulo( &
        x(i) - x_min - advection_velocity * final_time, x_max - x_min)
      exact_density = 1.0_dp + density_amplitude * &
        sin(2.0_dp * acos(-1.0_dp) * phase)
      density_l1_error = density_l1_error + &
        abs(conserved(irho, i) - exact_density)
    end do
    density_l1_error = density_l1_error / real(nx, dp)
    ok = density_l1_error > 0.0_dp .and. &
      density_l1_error < huge(1.0_dp)
  end subroutine run_entropy_wave

end program test_entropy_wave_pelec_plm4
