program test_reactive_material_contact_2d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved, reactive_conserved_to_primitive
  use reactive_2d_mod, only: advance_reactive_hydro_2d
  implicit none

  type(nasa7_species), allocatable :: species(:)
  real(dp) :: plain_error, steepened_error
  logical :: ok

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "2D contact thermodynamic database load failed"
  call run_contact(.false., plain_error)
  call run_contact(.true., steepened_error)

  write(*, '(a,1x,es16.8)') &
    "2D characteristic-PPM contact Y_H2 L1:", plain_error
  write(*, '(a,1x,es16.8)') &
    "2D steepened characteristic-PPM Y_H2 L1:", steepened_error
  if (steepened_error >= 0.75_dp * plain_error) &
    error stop "2D contact steepening did not sharpen the material interface"

contains

  subroutine run_contact(contact_steepening, l1_error)
    logical, intent(in) :: contact_steepening
    real(dp), intent(out) :: l1_error
    integer, parameter :: nx = 100, ny = 4
    real(dp), parameter :: length = 0.012_dp
    real(dp), parameter :: final_time = 1.0e-5_dp
    real(dp), parameter :: pressure = 101325.0_dp
    real(dp), parameter :: initial_temperature = 1000.0_dp
    real(dp), parameter :: velocity = 150.0_dp
    real(dp), parameter :: cfl = 0.35_dp
    real(dp) :: high_x(7), low_x(7), high_y(7), low_y(7)
    real(dp), allocatable :: state(:, :, :), temperature(:, :), q(:), primitive(:)
    real(dp) :: rho_high, rho_low, dx, dy, dt, time, x, source_x, exact_y
    real(dp) :: local_temperature, sound_speed, c_high, c_low, theta
    logical :: local_ok
    integer :: i, j, k, step

    high_x = [0.35570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
      1.0e-5_dp, 0.0_dp, 0.49643_dp]
    low_x = [0.23570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
      1.0e-5_dp, 0.0_dp, 0.61643_dp]
    call mass_fractions_from_mole_fractions(species, high_x, high_y, local_ok)
    if (.not. local_ok) error stop "2D high composition conversion failed"
    call mass_fractions_from_mole_fractions(species, low_x, low_y, local_ok)
    if (.not. local_ok) error stop "2D low composition conversion failed"
    rho_high = mixture_density( &
      species, high_y, pressure, initial_temperature, local_ok)
    if (.not. local_ok) error stop "2D high density construction failed"
    rho_low = mixture_density( &
      species, low_y, pressure, initial_temperature, local_ok)
    if (.not. local_ok) error stop "2D low density construction failed"

    allocate(state(reactive_nvar(7), nx, ny), temperature(nx, ny))
    allocate(q(reactive_nprim(7)), primitive(reactive_nprim(7)))
    dx = length / real(nx, dp)
    dy = 0.002_dp / real(ny, dp)
    do i = 1, nx
      x = (real(i, dp) - 0.5_dp) * dx
      if (x < 0.5_dp * length) then
        call make_state(rho_high, velocity, pressure, high_y, q, state(:, i, 1), &
          local_temperature, c_high)
      else
        call make_state(rho_low, velocity, pressure, low_y, q, state(:, i, 1), &
          local_temperature, c_low)
      end if
      do j = 1, ny
        state(:, i, j) = state(:, i, 1)
        temperature(i, j) = local_temperature
      end do
    end do

    time = 0.0_dp
    step = 0
    do while (time < final_time - 1.0e-16_dp)
      dt = cfl * dx / max(abs(velocity) + c_high, abs(velocity) + c_low)
      dt = min(dt, final_time - time)
      call advance_reactive_hydro_2d( &
        species, state, temperature, nx, ny, dx, dy, dt, &
        "characteristic_ppm", "mc", "hllc", .true., local_ok, theta, &
        contact_steepening, .false.)
      if (.not. local_ok) error stop "2D material-contact update failed"
      if (theta < 0.999999999_dp) &
        error stop "unexpected 2D contact transverse limiter activation"
      time = time + dt
      step = step + 1
      if (step > 10000) error stop "2D material-contact step limit reached"
    end do

    l1_error = 0.0_dp
    do i = 1, nx
      x = (real(i, dp) - 0.5_dp) * dx
      source_x = modulo(x - velocity * final_time, length)
      exact_y = merge(high_y(1), low_y(1), source_x < 0.5_dp * length)
      call reactive_conserved_to_primitive( &
        species, state(:, i, 1), temperature(i, 1), primitive, &
        local_temperature, sound_speed, local_ok)
      if (.not. local_ok) error stop "invalid 2D material-contact state"
      l1_error = l1_error + abs( &
        primitive(reactive_mass_fraction_component(1)) - exact_y)
    end do
    l1_error = l1_error / real(nx, dp)
  end subroutine run_contact

  subroutine make_state( &
      rho, velocity_value, pressure_value, mass_fractions, q, conserved, &
      temperature_value, c)
    real(dp), intent(in) :: rho, velocity_value, pressure_value
    real(dp), intent(in) :: mass_fractions(:)
    real(dp), intent(out) :: q(:), conserved(:), temperature_value, c
    logical :: local_ok
    integer :: k

    q(1:5) = [rho, velocity_value, 0.0_dp, 0.0_dp, pressure_value]
    do k = 1, size(mass_fractions)
      q(reactive_mass_fraction_component(k)) = mass_fractions(k)
    end do
    call reactive_primitive_to_conserved( &
      species, q, conserved, temperature_value, c, local_ok)
    if (.not. local_ok) error stop "2D material-contact state construction failed"
  end subroutine make_state
end program test_reactive_material_contact_2d
