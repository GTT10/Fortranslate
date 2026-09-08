program test_entropy_wave_3d
  use precision_mod, only: dp
  use state_indices_mod, only: irho, ncons, nprim, qrho
  use mesh_3d_mod, only: uniform_cell_centers_3d
  use simulation_config_3d_mod, only: entropy_wave_3d_config
  use entropy_wave_3d_problem_mod, only: &
    initialize_entropy_wave_3d, entropy_wave_3d_primitive
  use finite_volume_3d_mod, only: &
    compute_euler_cfl_timestep_3d, advance_euler_ssprk2_3d
  use diagnostics_3d_mod, only: integrated_conserved_quantities_3d
  implicit none

  integer, parameter :: case_count = 3
  integer, parameter :: resolutions(case_count) = [16, 32, 64]
  real(dp), parameter :: minimum_order = 0.85_dp
  real(dp), parameter :: conservation_tolerance = 2.0e-11_dp
  real(dp) :: errors(case_count), conservation_errors(case_count), order
  logical :: ok
  integer :: case_index

  do case_index = 1, case_count
    call run_case( &
      resolutions(case_index), errors(case_index), &
      conservation_errors(case_index), ok)
    if (.not. ok) error stop "3D entropy-wave run failed"
    write(*, '(a,i0,2(a,es24.16))') &
      "n=", resolutions(case_index), ", density L1=", errors(case_index), &
      ", conservation=", conservation_errors(case_index)
    if (conservation_errors(case_index) > conservation_tolerance) then
      error stop "3D periodic conservation error exceeds threshold"
    end if
  end do

  do case_index = 1, case_count - 1
    order = log(errors(case_index) / errors(case_index + 1)) / log(2.0_dp)
    write(*, '(a,i0,a,es24.16)') &
      "refinement pair ", case_index, ", order=", order
    if (order < minimum_order) then
      error stop "3D entropy-wave observed order is below threshold"
    end if
  end do

  write(*, '(a)') "test_entropy_wave_3d: PASS"

contains

  subroutine run_case(n, density_l1_error, conservation_error, ok)
    integer, intent(in) :: n
    real(dp), intent(out) :: density_l1_error, conservation_error
    logical, intent(out) :: ok

    real(dp), parameter :: x_min = 0.0_dp, x_max = 1.0_dp
    real(dp), parameter :: y_min = 0.0_dp, y_max = 1.0_dp
    real(dp), parameter :: z_min = 0.0_dp, z_max = 1.0_dp
    real(dp), parameter :: gamma = 1.4_dp, cfl = 0.4_dp
    real(dp), parameter :: final_time = 0.05_dp
    real(dp), allocatable :: x(:), y(:), z(:), conserved(:, :, :, :)
    real(dp) :: dx, dy, dz, time, dt
    real(dp) :: initial_totals(ncons), final_totals(ncons)
    real(dp) :: exact_primitive(nprim)
    type(entropy_wave_3d_config) :: wave
    logical :: step_ok, exact_ok
    integer :: i, j, k, step

    wave = entropy_wave_3d_config()
    allocate(x(n), y(n), z(n), conserved(ncons, n, n, n))
    call uniform_cell_centers_3d( &
      n, n, n, x_min, x_max, y_min, y_max, z_min, z_max, &
      x, y, z, dx, dy, dz)
    call initialize_entropy_wave_3d( &
      x, y, z, n, n, n, x_min, x_max, y_min, y_max, z_min, z_max, &
      gamma, wave, conserved, ok)
    if (.not. ok) return
    call integrated_conserved_quantities_3d( &
      conserved, n, n, n, dx, dy, dz, initial_totals)

    time = 0.0_dp
    step = 0
    do while (time < final_time)
      call compute_euler_cfl_timestep_3d( &
        conserved, n, n, n, dx, dy, dz, gamma, cfl, dt, step_ok)
      if (.not. step_ok) then
        ok = .false.
        return
      end if
      dt = min(dt, final_time - time)
      call advance_euler_ssprk2_3d( &
        conserved, n, n, n, dx, dy, dz, dt, gamma, "rusanov", step_ok)
      if (.not. step_ok) then
        ok = .false.
        return
      end if
      time = time + dt
      step = step + 1
      if (step > 10000) then
        ok = .false.
        return
      end if
    end do

    call integrated_conserved_quantities_3d( &
      conserved, n, n, n, dx, dy, dz, final_totals)
    conservation_error = maxval(abs(final_totals - initial_totals) / &
      max(1.0_dp, abs(initial_totals)))

    density_l1_error = 0.0_dp
    do k = 1, n
      do j = 1, n
        do i = 1, n
          call entropy_wave_3d_primitive( &
            x(i), y(j), z(k), final_time, &
            x_min, x_max, y_min, y_max, z_min, z_max, &
            wave, exact_primitive, exact_ok)
          if (.not. exact_ok) then
            ok = .false.
            return
          end if
          density_l1_error = density_l1_error + &
            abs(conserved(irho, i, j, k) - exact_primitive(qrho))
        end do
      end do
    end do
    density_l1_error = density_l1_error / real(n * n * n, dp)
    ok = density_l1_error > 0.0_dp .and. &
      density_l1_error < huge(1.0_dp)
  end subroutine run_case

end program test_entropy_wave_3d
