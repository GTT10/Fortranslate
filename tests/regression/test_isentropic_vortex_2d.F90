program test_isentropic_vortex_2d
  use precision_mod, only: dp
  use state_indices_mod, only: ncons, nprim, qrho, qp
  use mesh_2d_mod, only: uniform_cell_centers_2d
  use simulation_config_2d_mod, only: isentropic_vortex_config
  use isentropic_vortex_problem_mod, only: &
    initialize_isentropic_vortex, isentropic_vortex_primitive
  use state_conversion_mod, only: conserved_to_primitive
  use ctu_2d_mod, only: compute_cfl_timestep_2d, advance_ctu_2d
  use diagnostics_2d_mod, only: integrated_conserved_quantities_2d
  implicit none

  integer, parameter :: number_of_resolutions = 3
  integer, parameter :: resolutions(number_of_resolutions) = [24, 48, 96]
  real(dp), parameter :: minimum_order = 1.8_dp
  real(dp), parameter :: maximum_conservation_error = 2.0e-10_dp
  real(dp) :: density_errors(number_of_resolutions)
  real(dp) :: pressure_errors(number_of_resolutions)
  real(dp) :: conservation_errors(number_of_resolutions)
  real(dp) :: minimum_thetas(number_of_resolutions)
  real(dp) :: uncorrected_error, dummy_pressure, dummy_conservation, dummy_theta
  real(dp) :: order
  logical :: ok
  integer :: case_index

  do case_index = 1, number_of_resolutions
    call run_vortex( &
      resolutions(case_index), .true., density_errors(case_index), &
      pressure_errors(case_index), conservation_errors(case_index), &
      minimum_thetas(case_index), ok)
    if (.not. ok) error stop "2D isentropic-vortex run failed"

    write(*, '(a,i0,4(a,es24.16))') &
      "nx=", resolutions(case_index), &
      ", density L1=", density_errors(case_index), &
      ", pressure L1=", pressure_errors(case_index), &
      ", conservation=", conservation_errors(case_index), &
      ", min theta=", minimum_thetas(case_index)

    if (conservation_errors(case_index) > maximum_conservation_error) then
      error stop "2D periodic conservation error exceeds threshold"
    end if
  end do

  do case_index = 1, number_of_resolutions - 1
    order = log(density_errors(case_index) / &
      density_errors(case_index + 1)) / log(2.0_dp)
    write(*, '(a,i0,a,es24.16)') &
      "refinement pair ", case_index, ", order=", order
    if (order < minimum_order) then
      error stop "2D vortex observed order is below threshold"
    end if
  end do

  call run_vortex( &
    48, .false., uncorrected_error, dummy_pressure, dummy_conservation, &
    dummy_theta, ok)
  if (.not. ok) error stop "2D no-transverse comparison run failed"
  write(*, '(a,es24.16)') &
    "nx=48 no-transverse density L1=", uncorrected_error

  if (density_errors(2) >= 0.65_dp * uncorrected_error) then
    error stop "transverse correction did not materially improve vortex error"
  end if

  write(*, '(a)') "test_isentropic_vortex_2d: PASS"

contains

  subroutine run_vortex( &
      nx, use_transverse, density_l1_error, pressure_l1_error, &
      conservation_error, minimum_theta, ok)
    integer, intent(in) :: nx
    logical, intent(in) :: use_transverse
    real(dp), intent(out) :: density_l1_error, pressure_l1_error
    real(dp), intent(out) :: conservation_error, minimum_theta
    logical, intent(out) :: ok

    integer :: ny, i, j
    real(dp), parameter :: x_min = 0.0_dp, x_max = 10.0_dp
    real(dp), parameter :: y_min = 0.0_dp, y_max = 10.0_dp
    real(dp), parameter :: final_time = 1.0_dp
    real(dp), parameter :: gamma = 1.4_dp, cfl = 0.4_dp
    real(dp), allocatable :: x(:), y(:), conserved(:, :, :)
    real(dp) :: dx, dy, time, dt, step_theta
    real(dp) :: primitive(nprim), exact(nprim)
    real(dp) :: initial_totals(ncons), final_totals(ncons)
    type(isentropic_vortex_config) :: vortex
    logical :: cell_ok, step_ok

    ny = nx
    vortex = isentropic_vortex_config()
    allocate(x(nx), y(ny), conserved(ncons, nx, ny))

    call uniform_cell_centers_2d( &
      nx, ny, x_min, x_max, y_min, y_max, x, y, dx, dy)
    call initialize_isentropic_vortex( &
      x, y, nx, ny, x_min, x_max, y_min, y_max, gamma, vortex, &
      conserved, ok)
    if (.not. ok) return

    call integrated_conserved_quantities_2d( &
      conserved, nx, ny, dx, dy, initial_totals)

    time = 0.0_dp
    minimum_theta = 1.0_dp
    do while (time < final_time)
      call compute_cfl_timestep_2d( &
        conserved, nx, ny, dx, dy, gamma, cfl, dt, step_ok)
      if (.not. step_ok) then
        ok = .false.
        return
      end if

      dt = min(dt, final_time - time)
      call advance_ctu_2d( &
        conserved, nx, ny, dx, dy, dt, gamma, "mc", "pelec", &
        use_transverse, step_ok, step_theta)
      if (.not. step_ok) then
        ok = .false.
        return
      end if
      minimum_theta = min(minimum_theta, step_theta)
      time = time + dt
    end do

    call integrated_conserved_quantities_2d( &
      conserved, nx, ny, dx, dy, final_totals)
    conservation_error = maxval(abs(final_totals - initial_totals))

    density_l1_error = 0.0_dp
    pressure_l1_error = 0.0_dp
    do j = 1, ny
      do i = 1, nx
        call conserved_to_primitive( &
          conserved(:, i, j), gamma, primitive, cell_ok)
        if (.not. cell_ok) then
          ok = .false.
          return
        end if
        call isentropic_vortex_primitive( &
          x(i), y(j), final_time, x_min, x_max, y_min, y_max, gamma, &
          vortex, exact, cell_ok)
        if (.not. cell_ok) then
          ok = .false.
          return
        end if
        density_l1_error = density_l1_error + &
          abs(primitive(qrho) - exact(qrho))
        pressure_l1_error = pressure_l1_error + &
          abs(primitive(qp) - exact(qp))
      end do
    end do

    density_l1_error = density_l1_error / real(nx * ny, dp)
    pressure_l1_error = pressure_l1_error / real(nx * ny, dp)
    ok = density_l1_error > 0.0_dp .and. pressure_l1_error > 0.0_dp
  end subroutine run_vortex

end program test_isentropic_vortex_2d
