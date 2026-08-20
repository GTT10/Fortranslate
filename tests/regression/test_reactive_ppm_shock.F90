program test_reactive_ppm_shock
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved, reactive_conserved_to_primitive, &
    advance_reactive_hydro, reactive_integrals
  implicit none

  type(nasa7_species), allocatable :: species(:)
  real(dp) :: unflattened_tv, flattened_tv, unflattened_overshoot
  real(dp) :: flattened_overshoot, solution_difference
  logical :: ok

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "Failed to load reactive shock thermodynamics"

  call run_shock(.false., unflattened_tv, unflattened_overshoot, &
    solution_difference, .true.)
  call run_shock(.true., flattened_tv, flattened_overshoot, &
    solution_difference, .false.)

  write(*, '(a,1x,es16.8)') "unflattened pressure TV:", unflattened_tv
  write(*, '(a,1x,es16.8)') "flattened pressure TV:", flattened_tv
  write(*, '(a,1x,es16.8)') "unflattened pressure overshoot:", &
    unflattened_overshoot
  write(*, '(a,1x,es16.8)') "flattened pressure overshoot:", &
    flattened_overshoot
  write(*, '(a,1x,es16.8)') "flattened/unflattened state difference:", &
    solution_difference

  if (max(unflattened_overshoot, flattened_overshoot) > 1.0e-8_dp * 3.0e5_dp) &
    error stop "Reactive PPM shock produced a pressure overshoot"
  if (flattened_tv > 1.05_dp * unflattened_tv) then
    error stop "Reactive PPM shock flattening increased pressure variation"
  end if
  if (solution_difference <= 1.0e-6_dp) then
    error stop "Reactive PPM shock flattening was not activated"
  end if

  write(*, '(a)') "test_reactive_ppm_shock: PASS"

contains

  subroutine run_shock(use_flattening, pressure_tv, pressure_overshoot, &
      state_difference, store_reference)
    logical, intent(in) :: use_flattening, store_reference
    real(dp), intent(out) :: pressure_tv, pressure_overshoot, state_difference
    integer, parameter :: nx = 200
    real(dp), parameter :: x_lower = 0.0_dp
    real(dp), parameter :: x_upper = 1.0_dp
    real(dp), parameter :: final_time = 5.0e-5_dp
    real(dp), parameter :: cfl = 0.10_dp
    real(dp), parameter :: left_density = 1.0_dp
    real(dp), parameter :: right_density = 0.125_dp
    real(dp), parameter :: left_pressure = 3.0e5_dp
    real(dp), parameter :: right_pressure = 1.0e5_dp
    real(dp), save, allocatable :: reference_state(:, :)
    real(dp), allocatable :: state(:, :), temperature(:), primitive(:), q(:)
    real(dp) :: y(7), dx, time, dt, max_speed, local_t, sound_speed, x
    real(dp) :: initial_integrals(5), final_integrals(5), conservation_error
    logical :: local_ok
    integer :: cell, step

    allocate(state(reactive_nvar(7), 0:nx + 1), temperature(0:nx + 1))
    allocate(primitive(reactive_nprim(7)), q(reactive_nprim(7)))
    y = 0.0_dp
    y(7) = 1.0_dp
    dx = (x_upper - x_lower) / real(nx, dp)
    do cell = 1, nx
      x = x_lower + (real(cell, dp) - 0.5_dp) * dx
      if (x < 0.5_dp) then
        call make_state(left_density, left_pressure, y, q, state(:, cell), &
          temperature(cell), sound_speed)
      else
        call make_state(right_density, right_pressure, y, q, state(:, cell), &
          temperature(cell), sound_speed)
      end if
    end do
    state(:, 0) = state(:, 1)
    state(:, nx + 1) = state(:, nx)
    temperature(0) = temperature(1)
    temperature(nx + 1) = temperature(nx)
    call reactive_integrals(state, nx, dx, initial_integrals)

    time = 0.0_dp
    step = 0
    do while (time < final_time - 1.0e-16_dp)
      max_speed = 0.0_dp
      do cell = 1, nx
        call reactive_conserved_to_primitive( &
          species, state(:, cell), temperature(cell), primitive, local_t, &
          sound_speed, local_ok)
        if (.not. local_ok) error stop "Invalid reactive shock state"
        max_speed = max(max_speed, abs(primitive(2)) + sound_speed)
      end do
      dt = min(cfl * dx / max_speed, final_time - time)
      call advance_reactive_hydro( &
        species, state, temperature, nx, dx, dt, "characteristic_ppm", &
        "mc", "periodic", local_ok, "hllc", .false., use_flattening)
      if (.not. local_ok) error stop "Reactive shock hydro step failed"
      time = time + dt
      step = step + 1
      if (step > 10000) error stop "Reactive shock step limit reached"
    end do

    pressure_tv = 0.0_dp
    pressure_overshoot = 0.0_dp
    do cell = 1, nx
      call reactive_conserved_to_primitive( &
        species, state(:, cell), temperature(cell), primitive, local_t, &
        sound_speed, local_ok)
      if (.not. local_ok) error stop "Invalid final reactive shock state"
      if (primitive(1) <= 0.0_dp .or. primitive(5) <= 0.0_dp .or. &
          local_t <= 0.0_dp) error stop "Reactive shock lost positivity"
      pressure_overshoot = max(pressure_overshoot, &
        max(0.0_dp, primitive(5) - left_pressure), &
        max(0.0_dp, right_pressure - primitive(5)))
      if (cell > 1) then
        q = primitive
        call reactive_conserved_to_primitive( &
          species, state(:, cell - 1), temperature(cell - 1), primitive, &
          local_t, sound_speed, local_ok)
        if (.not. local_ok) error stop "Invalid neighboring shock state"
        pressure_tv = pressure_tv + abs(q(5) - primitive(5))
      end if
    end do

    call reactive_integrals(state, nx, dx, final_integrals)
    conservation_error = maxval(abs(final_integrals - initial_integrals) / &
      max(1.0_dp, abs(initial_integrals)))
    if (conservation_error > 2.0e-11_dp) then
      error stop "Reactive shock lost conservation"
    end if

    if (store_reference) then
      if (allocated(reference_state)) deallocate(reference_state)
      allocate(reference_state(size(state, 1), 1:nx))
      reference_state = state(:, 1:nx)
      state_difference = 0.0_dp
    else
      if (.not. allocated(reference_state)) error stop "Missing shock reference"
      state_difference = sum(abs(state(:, 1:nx) - reference_state)) / &
        real(size(reference_state), dp)
      if (state_difference <= 1.0e-12_dp) then
        error stop "Shock flattening did not change the solution"
      end if
    end if
  end subroutine run_shock

  subroutine make_state(rho, pressure, y, primitive, conserved, temperature, &
      sound_speed)
    real(dp), intent(in) :: rho, pressure, y(:)
    real(dp), intent(out) :: primitive(:), conserved(:), temperature, sound_speed
    logical :: local_ok
    integer :: k

    primitive(1:5) = [rho, 0.0_dp, 0.0_dp, 0.0_dp, pressure]
    do k = 1, size(y)
      primitive(reactive_mass_fraction_component(k)) = y(k)
    end do
    call reactive_primitive_to_conserved( &
      species, primitive, conserved, temperature, sound_speed, local_ok)
    if (.not. local_ok) error stop "Failed to construct shock state"
  end subroutine make_state

end program test_reactive_ppm_shock
