module riemann_pelec_mod
  use precision_mod, only: dp
  use constants_mod, only: density_floor, pressure_floor, tiny_speed
  use state_indices_mod, only: &
    irho, imx, imy, imz, iet, ncons, nprim, qrho, qu, qv, qw, qp
  use state_conversion_mod, only: conserved_to_primitive
  use eos_ideal_mod, only: ideal_gas_sound_speed
  implicit none
  private

  public :: pelec_riemann_flux_x

contains

  pure subroutine pelec_riemann_flux_x( &
      left_state, right_state, gamma, flux, ok, interface_density, &
      interface_velocity, interface_pressure)
    real(dp), intent(in) :: left_state(ncons), right_state(ncons)
    real(dp), intent(in) :: gamma
    real(dp), intent(out) :: flux(ncons)
    logical, intent(out) :: ok
    real(dp), intent(out), optional :: interface_density
    real(dp), intent(out), optional :: interface_velocity
    real(dp), intent(out), optional :: interface_pressure

    real(dp) :: left_primitive(nprim), right_primitive(nprim)
    real(dp) :: rho_left, rho_right, velocity_left, velocity_right
    real(dp) :: transverse_y_left, transverse_y_right
    real(dp) :: transverse_z_left, transverse_z_right
    real(dp) :: pressure_left, pressure_right
    real(dp) :: sound_left, sound_right, impedance_left, impedance_right
    real(dp) :: pressure_star, velocity_star
    real(dp) :: rho_origin, velocity_origin, pressure_origin, sound_origin
    real(dp) :: rho_star, sound_star, density_increment
    real(dp) :: transverse_y, transverse_z
    real(dp) :: wave_out, wave_in, shock_speed, speed_difference
    real(dp) :: interpolation_fraction, average_sound_speed
    real(dp) :: rho_interface, velocity_interface, pressure_interface
    real(dp) :: internal_energy_density, total_energy_density
    real(dp) :: stationary_threshold, regularization
    logical :: left_ok, right_ok, stationary_interface

    flux = 0.0_dp
    ok = .false.

    call conserved_to_primitive(left_state, gamma, left_primitive, left_ok)
    call conserved_to_primitive(right_state, gamma, right_primitive, right_ok)
    if (.not. (left_ok .and. right_ok)) return

    rho_left = left_primitive(qrho)
    velocity_left = left_primitive(qu)
    transverse_y_left = left_primitive(qv)
    transverse_z_left = left_primitive(qw)
    pressure_left = left_primitive(qp)

    rho_right = right_primitive(qrho)
    velocity_right = right_primitive(qu)
    transverse_y_right = right_primitive(qv)
    transverse_z_right = right_primitive(qw)
    pressure_right = right_primitive(qp)

    sound_left = ideal_gas_sound_speed(rho_left, pressure_left, gamma)
    sound_right = ideal_gas_sound_speed(rho_right, pressure_right, gamma)
    if (sound_left <= 0.0_dp .or. sound_right <= 0.0_dp) return

    impedance_left = max(tiny(1.0_dp), sound_left * rho_left)
    impedance_right = max(tiny(1.0_dp), sound_right * rho_right)

    pressure_star = max(pressure_floor, &
      ((impedance_right * pressure_left + impedance_left * pressure_right) + &
       impedance_left * impedance_right * (velocity_left - velocity_right)) / &
      (impedance_left + impedance_right))
    velocity_star = &
      ((impedance_left * velocity_left + impedance_right * velocity_right) + &
       (pressure_left - pressure_right)) / &
      (impedance_left + impedance_right)

    if (velocity_star > 0.0_dp) then
      rho_origin = rho_left
      velocity_origin = velocity_left
      pressure_origin = pressure_left
      transverse_y = transverse_y_left
      transverse_z = transverse_z_left
    else if (velocity_star < 0.0_dp) then
      rho_origin = rho_right
      velocity_origin = velocity_right
      pressure_origin = pressure_right
      transverse_y = transverse_y_right
      transverse_z = transverse_z_right
    else
      rho_origin = 0.5_dp * (rho_left + rho_right)
      velocity_origin = 0.5_dp * (velocity_left + velocity_right)
      pressure_origin = 0.5_dp * (pressure_left + pressure_right)
      transverse_y = 0.5_dp * (transverse_y_left + transverse_y_right)
      transverse_z = 0.5_dp * (transverse_z_left + transverse_z_right)
    end if

    stationary_threshold = tiny_speed * &
      0.5_dp * (abs(velocity_left) + abs(velocity_right))
    stationary_interface = abs(velocity_star) <= stationary_threshold
    if (stationary_interface) then
      velocity_star = 0.0_dp
      rho_origin = 0.5_dp * (rho_left + rho_right)
      velocity_origin = 0.5_dp * (velocity_left + velocity_right)
      pressure_origin = 0.5_dp * (pressure_left + pressure_right)
      transverse_y = 0.5_dp * (transverse_y_left + transverse_y_right)
      transverse_z = 0.5_dp * (transverse_z_left + transverse_z_right)
    end if

    sound_origin = ideal_gas_sound_speed( &
      rho_origin, pressure_origin, gamma)
    if (sound_origin <= 0.0_dp) return

    density_increment = &
      (pressure_star - pressure_origin) / (sound_origin * sound_origin)
    rho_star = rho_origin + density_increment
    if (rho_star <= density_floor) return

    sound_star = ideal_gas_sound_speed(rho_star, pressure_star, gamma)
    if (sound_star <= 0.0_dp) return

    wave_out = sound_origin - sign(1.0_dp, velocity_star) * velocity_origin
    wave_in = sound_star - sign(1.0_dp, velocity_star) * velocity_star
    shock_speed = 0.5_dp * (wave_in + wave_out)

    if (pressure_star >= pressure_origin) then
      wave_out = shock_speed
      wave_in = shock_speed
    end if

    average_sound_speed = 0.5_dp * (sound_left + sound_right)
    regularization = sqrt(epsilon(1.0_dp)) * max(1.0_dp, average_sound_speed)
    if (abs(wave_out - wave_in) < regularization) then
      speed_difference = regularization
    else
      speed_difference = wave_out - wave_in
    end if

    interpolation_fraction = max(0.0_dp, min(1.0_dp, &
      0.5_dp * (1.0_dp + (wave_out + wave_in) / speed_difference)))

    rho_interface = interpolation_fraction * rho_star + &
      (1.0_dp - interpolation_fraction) * rho_origin
    velocity_interface = interpolation_fraction * velocity_star + &
      (1.0_dp - interpolation_fraction) * velocity_origin
    pressure_interface = interpolation_fraction * pressure_star + &
      (1.0_dp - interpolation_fraction) * pressure_origin

    if (wave_out < 0.0_dp) then
      rho_interface = rho_origin
      velocity_interface = velocity_origin
      pressure_interface = pressure_origin
    end if

    if (wave_in >= 0.0_dp) then
      rho_interface = rho_star
      velocity_interface = velocity_star
      pressure_interface = pressure_star
    end if

    if (rho_interface <= density_floor .or. &
        pressure_interface <= pressure_floor) return

    internal_energy_density = pressure_interface / (gamma - 1.0_dp)
    total_energy_density = internal_energy_density + &
      0.5_dp * rho_interface * &
      (velocity_interface**2 + transverse_y**2 + transverse_z**2)

    flux(irho) = rho_interface * velocity_interface
    flux(imx) = flux(irho) * velocity_interface + pressure_interface
    flux(imy) = flux(irho) * transverse_y
    flux(imz) = flux(irho) * transverse_z
    flux(iet) = velocity_interface * &
      (total_energy_density + pressure_interface)

    if (present(interface_density)) interface_density = rho_interface
    if (present(interface_velocity)) interface_velocity = velocity_interface
    if (present(interface_pressure)) interface_pressure = pressure_interface
    ok = .true.
  end subroutine pelec_riemann_flux_x

end module riemann_pelec_mod
