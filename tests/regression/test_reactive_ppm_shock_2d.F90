program test_reactive_ppm_shock_2d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved, reactive_conserved_to_primitive
  use reactive_2d_mod, only: &
    compute_reactive_cfl_timestep_2d, advance_reactive_hydro_2d, &
    reactive_integrals_2d
  implicit none

  integer, parameter :: n = 32
  type(nasa7_species), allocatable :: species(:)
  real(dp), allocatable, save :: reference_state(:, :, :)
  real(dp) :: unflattened_tv, flattened_tv
  real(dp) :: unflattened_overshoot, flattened_overshoot
  real(dp) :: solution_difference
  logical :: ok

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "thermodynamic database load failed"
  call run_case(.false., unflattened_tv, unflattened_overshoot, &
    solution_difference, .true.)
  call run_case(.true., flattened_tv, flattened_overshoot, &
    solution_difference, .false.)

  write(*, '(a,1x,es16.8)') "unflattened pressure TV:", unflattened_tv
  write(*, '(a,1x,es16.8)') "flattened pressure TV:", flattened_tv
  write(*, '(a,1x,es16.8)') "flattened/unflattened state difference:", &
    solution_difference
  if (max(unflattened_overshoot, flattened_overshoot) > 3.0e-3_dp) &
    error stop "2D characteristic PPM shock produced pressure overshoot"
  if (flattened_tv > 1.01_dp * unflattened_tv) &
    error stop "2D shock flattening materially increased pressure variation"
  if (solution_difference <= 1.0e-3_dp) &
    error stop "2D shock flattening was not activated"

contains

  subroutine run_case(use_flattening, pressure_tv, pressure_overshoot, &
      state_difference, store_reference)
    logical, intent(in) :: use_flattening, store_reference
    real(dp), intent(out) :: pressure_tv, pressure_overshoot, state_difference
    real(dp), allocatable :: state(:, :, :), temperature(:, :), primitive(:)
    real(dp) :: dx, dy, time, dt, theta, x, y, phase
    real(dp) :: sound_speed, local_temperature, pressure_here
    real(dp) :: initial_integrals(5), final_integrals(5), conservation
    logical :: local_ok
    integer :: i, j, k, im, jm, steps, nvar

    nvar = reactive_nvar(size(species))
    allocate(state(nvar, n, n), temperature(n, n))
    allocate(primitive(reactive_nprim(size(species))))
    dx = 1.0_dp / real(n, dp)
    dy = dx
    do j = 1, n
      y = (real(j, dp) - 0.5_dp) * dy
      do i = 1, n
        x = (real(i, dp) - 0.5_dp) * dx
        phase = modulo(x + y, 1.0_dp)
        primitive(1:5) = [merge(1.0_dp, 0.125_dp, phase < 0.5_dp), &
          0.0_dp, 0.0_dp, 0.0_dp, &
          merge(3.0e5_dp, 1.0e5_dp, phase < 0.5_dp)]
        do k = 1, size(species)
          primitive(reactive_mass_fraction_component(k)) = 0.0_dp
        end do
        primitive(reactive_mass_fraction_component(7)) = 1.0_dp
        call reactive_primitive_to_conserved( &
          species, primitive, state(:, i, j), temperature(i, j), &
          sound_speed, local_ok)
        if (.not. local_ok) error stop "2D shock state construction failed"
      end do
    end do
    call reactive_integrals_2d(state, n, n, dx, dy, initial_integrals)
    time = 0.0_dp
    steps = 0
    do while (time < 2.0e-5_dp - 1.0e-16_dp)
      call compute_reactive_cfl_timestep_2d( &
        species, state, temperature, n, n, dx, dy, 0.10_dp, dt, local_ok)
      if (.not. local_ok) error stop "2D shock timestep failed"
      dt = min(dt, 2.0e-5_dp - time)
      call advance_reactive_hydro_2d( &
        species, state, temperature, n, n, dx, dy, dt, &
        "characteristic_ppm", "mc", "hllc", .true., local_ok, theta, &
        .false., use_flattening)
      if (.not. local_ok) error stop "2D shock hydro update failed"
      time = time + dt
      steps = steps + 1
    end do

    pressure_tv = 0.0_dp
    pressure_overshoot = 0.0_dp
    do j = 1, n
      jm = 1 + modulo(j - 2, n)
      do i = 1, n
        im = 1 + modulo(i - 2, n)
        call reactive_conserved_to_primitive( &
          species, state(:, i, j), temperature(i, j), primitive, &
          local_temperature, sound_speed, local_ok)
        if (.not. local_ok) error stop "invalid final 2D shock state"
        if (primitive(1) <= 0.0_dp .or. primitive(5) <= 0.0_dp .or. &
            local_temperature <= 0.0_dp) error stop "2D shock lost positivity"
        pressure_here = primitive(5)
        pressure_overshoot = max(pressure_overshoot, &
          max(0.0_dp, pressure_here - 3.0e5_dp), &
          max(0.0_dp, 1.0e5_dp - pressure_here))
        call reactive_conserved_to_primitive( &
          species, state(:, im, j), temperature(im, j), primitive, &
          local_temperature, sound_speed, local_ok)
        if (.not. local_ok) error stop "invalid x-neighbor shock state"
        pressure_tv = pressure_tv + abs(pressure_here - primitive(5))
        call reactive_conserved_to_primitive( &
          species, state(:, i, jm), temperature(i, jm), primitive, &
          local_temperature, sound_speed, local_ok)
        if (.not. local_ok) error stop "invalid y-neighbor shock state"
        pressure_tv = pressure_tv + abs(pressure_here - primitive(5))
      end do
    end do
    call reactive_integrals_2d(state, n, n, dx, dy, final_integrals)
    conservation = maxval(abs(final_integrals - initial_integrals) / &
      max(1.0_dp, abs(initial_integrals)))
    if (conservation > 3.0e-11_dp) error stop "2D shock lost conservation"

    if (store_reference) then
      allocate(reference_state(nvar, n, n))
      reference_state = state
      state_difference = 0.0_dp
    else
      state_difference = sum(abs(state - reference_state)) / real(size(state), dp)
    end if
  end subroutine run_case
end program test_reactive_ppm_shock_2d
