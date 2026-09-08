program test_reactive_entropy_wave_plm_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mesh_3d_mod, only: uniform_cell_centers_3d
  use reactive_1d_mod, only: reactive_nvar
  use simulation_config_reactive_3d_mod, only: reactive_3d_config
  use reactive_entropy_wave_3d_problem_mod, only: &
    initialize_reactive_entropy_wave_3d, reactive_entropy_wave_density_3d
  use reactive_3d_mod, only: &
    compute_reactive_cfl_timestep_3d, &
    advance_reactive_euler_ssprk2_plm_3d, reactive_integrals_3d, &
    reactive_extrema_3d
  implicit none

  integer, parameter :: case_count = 3
  integer, parameter :: resolutions(case_count) = [6, 12, 24]
  real(dp), parameter :: minimum_order = 1.65_dp
  real(dp), parameter :: conservation_tolerance = 5.0e-11_dp
  real(dp), parameter :: closure_tolerance = 5.0e-11_dp
  type(nasa7_species), allocatable :: species(:)
  real(dp) :: errors(case_count), conservation_errors(case_count)
  real(dp) :: closure_errors(case_count), order
  logical :: ok
  integer :: case_index

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "Failed to load elementary thermodynamics"
  do case_index = 1, case_count
    call run_case( &
      species, resolutions(case_index), errors(case_index), &
      conservation_errors(case_index), closure_errors(case_index), ok)
    if (.not. ok) error stop "Reactive PLM 3D entropy-wave run failed"
    write(*, '(a,i0,3(a,es24.16))') &
      "n=", resolutions(case_index), ", density L1=", errors(case_index), &
      ", conservation=", conservation_errors(case_index), &
      ", closure=", closure_errors(case_index)
    if (conservation_errors(case_index) > conservation_tolerance) then
      error stop "Reactive PLM 3D conservation error exceeds threshold"
    end if
    if (closure_errors(case_index) > closure_tolerance) then
      error stop "Reactive PLM 3D species closure exceeds threshold"
    end if
  end do

  do case_index = 1, case_count - 1
    order = log(errors(case_index) / errors(case_index + 1)) / log(2.0_dp)
    write(*, '(a,i0,a,es24.16)') &
      "PLM refinement pair ", case_index, ", order=", order
    if (order < minimum_order) then
      error stop "Reactive PLM 3D entropy-wave order is below threshold"
    end if
  end do
  write(*, '(a)') "test_reactive_entropy_wave_plm_3d: PASS"

contains

  subroutine run_case( &
      species, n, density_l1_error, conservation_error, closure_error, ok)
    type(nasa7_species), intent(in) :: species(:)
    integer, intent(in) :: n
    real(dp), intent(out) :: density_l1_error, conservation_error
    real(dp), intent(out) :: closure_error
    logical, intent(out) :: ok

    type(reactive_3d_config) :: config
    real(dp), allocatable :: x(:), y(:), z(:)
    real(dp), allocatable :: state(:, :, :, :), temperature(:, :, :)
    real(dp), allocatable :: initial_totals(:), final_totals(:)
    real(dp), allocatable :: mass_fractions(:)
    real(dp) :: dx, dy, dz, time, dt, base_density, exact_density
    real(dp) :: minimum_density, maximum_density
    real(dp) :: minimum_pressure, maximum_pressure
    real(dp) :: minimum_temperature, maximum_temperature, maximum_speed
    logical :: step_ok
    integer :: i, j, k, step, nvar

    config = reactive_3d_config()
    config%nx = n
    config%ny = n
    config%nz = n
    config%thermo_model = "elementary"
    config%final_time = 5.0e-7_dp
    config%cfl = 0.30_dp
    nvar = reactive_nvar(size(species))
    allocate(x(n), y(n), z(n))
    allocate(state(nvar, n, n, n), temperature(n, n, n))
    allocate(initial_totals(nvar), final_totals(nvar))
    allocate(mass_fractions(size(species)))
    call uniform_cell_centers_3d( &
      n, n, n, config%x_lower, config%x_upper, &
      config%y_lower, config%y_upper, config%z_lower, config%z_upper, &
      x, y, z, dx, dy, dz)
    call initialize_reactive_entropy_wave_3d( &
      species, config, x, y, z, state, temperature, base_density, &
      mass_fractions, ok)
    if (.not. ok) return
    call reactive_integrals_3d( &
      state, n, n, n, dx, dy, dz, initial_totals, ok)
    if (.not. ok) return

    time = 0.0_dp
    step = 0
    do while (time < config%final_time)
      call compute_reactive_cfl_timestep_3d( &
        species, state, temperature, n, n, n, dx, dy, dz, &
        config%cfl, dt, step_ok)
      if (.not. step_ok) then
        ok = .false.
        return
      end if
      dt = min(dt, config%final_time - time)
      call advance_reactive_euler_ssprk2_plm_3d( &
        species, state, temperature, n, n, n, dx, dy, dz, dt, &
        "mc", "rusanov", step_ok)
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

    call reactive_integrals_3d( &
      state, n, n, n, dx, dy, dz, final_totals, ok)
    if (.not. ok) return
    conservation_error = maxval( &
      abs(final_totals - initial_totals) / max(1.0_dp, abs(initial_totals)))
    call reactive_extrema_3d( &
      species, state, temperature, n, n, n, &
      minimum_density, maximum_density, minimum_pressure, maximum_pressure, &
      minimum_temperature, maximum_temperature, maximum_speed, &
      closure_error, ok)
    if (.not. ok) return
    if (minimum_density <= 0.0_dp .or. minimum_pressure <= 0.0_dp .or. &
        minimum_temperature <= 0.0_dp .or. maximum_density <= 0.0_dp .or. &
        maximum_pressure <= 0.0_dp .or. maximum_temperature <= 0.0_dp .or. &
        maximum_speed <= 0.0_dp) then
      ok = .false.
      return
    end if

    density_l1_error = 0.0_dp
    do k = 1, n
      do j = 1, n
        do i = 1, n
          exact_density = reactive_entropy_wave_density_3d( &
            x(i), y(j), z(k), time, config, base_density)
          density_l1_error = density_l1_error + &
            abs(state(1, i, j, k) - exact_density)
        end do
      end do
    end do
    density_l1_error = density_l1_error / real(n * n * n, dp)
    ok = density_l1_error > 0.0_dp .and. &
      density_l1_error < huge(1.0_dp)
  end subroutine run_case

end program test_reactive_entropy_wave_plm_3d
