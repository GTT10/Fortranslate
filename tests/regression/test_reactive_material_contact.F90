program test_reactive_material_contact
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved, reactive_conserved_to_primitive, &
    advance_reactive_hydro
  implicit none

  type(nasa7_species), allocatable :: species(:)
  real(dp) :: hllc_error, ppm_error, rusanov_error, improvement
  logical :: ok

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "Failed to load material-contact thermodynamics"
  call run_contact("characteristic_plm", "hllc", hllc_error)
  call run_contact("ppm", "hllc", ppm_error)
  call run_contact("characteristic_plm", "rusanov", rusanov_error)
  improvement = rusanov_error / hllc_error
  write(*, '(a,1x,es16.8)') "HLLC material-contact Y_H2 L1:", hllc_error
  write(*, '(a,1x,es16.8)') "PPM HLLC material-contact Y_H2 L1:", ppm_error
  write(*, '(a,1x,es16.8)') "Rusanov material-contact Y_H2 L1:", rusanov_error
  write(*, '(a,1x,f10.4)') "HLLC material-contact improvement:", improvement
  if (improvement < 1.20_dp) then
    error stop "HLLC did not improve material-contact resolution"
  end if
  if (ppm_error >= 0.90_dp * hllc_error) then
    error stop "Reactive PPM did not sharpen the material contact"
  end if
  write(*, '(a)') "test_reactive_material_contact: PASS"

contains

  subroutine run_contact(reconstruction, solver, l1_error)
    character(len=*), intent(in) :: reconstruction, solver
    real(dp), intent(out) :: l1_error
    integer, parameter :: nx = 200
    real(dp), parameter :: x_lower = 0.0_dp
    real(dp), parameter :: x_upper = 0.012_dp
    real(dp), parameter :: contact_initial = 0.004_dp
    real(dp), parameter :: final_time = 2.0e-5_dp
    real(dp), parameter :: pressure = 101325.0_dp
    real(dp), parameter :: initial_temperature = 1000.0_dp
    real(dp), parameter :: velocity = 150.0_dp
    real(dp), parameter :: cfl = 0.4_dp
    real(dp) :: xl(7), xr(7), yl(7), yr(7)
    real(dp), allocatable :: state(:, :), temperature(:), q(:), primitive(:)
    real(dp) :: rho_l, rho_r, tl, tr, cl, cr, dx, dt, time, x, exact_y
    real(dp) :: local_temperature, sound_speed
    logical :: local_ok
    integer :: cell, step

    xl = [0.35570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
      1.0e-5_dp, 0.0_dp, 0.49643_dp]
    xr = [0.23570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
      1.0e-5_dp, 0.0_dp, 0.61643_dp]
    call mass_fractions_from_mole_fractions(species, xl, yl, local_ok)
    if (.not. local_ok) error stop "Failed to build left material composition"
    call mass_fractions_from_mole_fractions(species, xr, yr, local_ok)
    if (.not. local_ok) error stop "Failed to build right material composition"
    rho_l = mixture_density(species, yl, pressure, initial_temperature, local_ok)
    if (.not. local_ok) error stop "Failed to build left material density"
    rho_r = mixture_density(species, yr, pressure, initial_temperature, local_ok)
    if (.not. local_ok) error stop "Failed to build right material density"

    allocate(state(reactive_nvar(7), 0:nx + 1), temperature(0:nx + 1))
    allocate(q(reactive_nprim(7)), primitive(reactive_nprim(7)))
    dx = (x_upper - x_lower) / real(nx, dp)
    do cell = 1, nx
      x = x_lower + (real(cell, dp) - 0.5_dp) * dx
      if (x < contact_initial) then
        call make_state(rho_l, velocity, pressure, yl, q, state(:, cell), tl, cl)
        temperature(cell) = tl
      else
        call make_state(rho_r, velocity, pressure, yr, q, state(:, cell), tr, cr)
        temperature(cell) = tr
      end if
    end do
    state(:, 0) = state(:, 1)
    state(:, nx + 1) = state(:, nx)
    temperature(0) = temperature(1)
    temperature(nx + 1) = temperature(nx)

    time = 0.0_dp
    step = 0
    do while (time < final_time - 1.0e-16_dp)
      dt = cfl * dx / max(abs(velocity) + cl, abs(velocity) + cr)
      dt = min(dt, final_time - time)
      call advance_reactive_hydro( &
        species, state, temperature, nx, dx, dt, reconstruction, &
        "mc", "outflow", local_ok, solver)
      if (.not. local_ok) error stop "Material-contact hydro step failed"
      time = time + dt
      step = step + 1
      if (step > 10000) error stop "Material-contact step limit reached"
    end do

    l1_error = 0.0_dp
    do cell = 1, nx
      x = x_lower + (real(cell, dp) - 0.5_dp) * dx
      if (x < contact_initial + velocity * final_time) then
        exact_y = yl(1)
      else
        exact_y = yr(1)
      end if
      call reactive_conserved_to_primitive( &
        species, state(:, cell), temperature(cell), primitive, &
        local_temperature, sound_speed, local_ok)
      if (.not. local_ok) error stop "Invalid material-contact state"
      l1_error = l1_error + abs( &
        primitive(reactive_mass_fraction_component(1)) - exact_y)
    end do
    l1_error = l1_error / real(nx, dp)
  end subroutine run_contact

  subroutine make_state(rho, velocity_value, pressure_value, y, q, conserved, temperature, sound_speed)
    real(dp), intent(in) :: rho, velocity_value, pressure_value, y(:)
    real(dp), intent(out) :: q(:), conserved(:), temperature, sound_speed
    logical :: local_ok
    integer :: k

    q(1:5) = [rho, velocity_value, 0.0_dp, 0.0_dp, pressure_value]
    do k = 1, size(y)
      q(reactive_mass_fraction_component(k)) = y(k)
    end do
    call reactive_primitive_to_conserved( &
      species, q, conserved, temperature, sound_speed, local_ok)
    if (.not. local_ok) error stop "Failed to construct material-contact state"
  end subroutine make_state

end program test_reactive_material_contact
