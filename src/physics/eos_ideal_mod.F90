module eos_ideal_mod
  use precision_mod, only: dp
  use constants_mod, only: density_floor, pressure_floor
  use state_indices_mod, only: irho, imx, imy, imz, iet, ncons
  implicit none
  private

  public :: ideal_gas_pressure
  public :: ideal_gas_sound_speed
  public :: ideal_gas_specific_internal_energy
  public :: ideal_gas_total_energy_density
  public :: valid_ideal_gas_gamma

contains

  pure logical function valid_ideal_gas_gamma(gamma) result(valid)
    real(dp), intent(in) :: gamma
    valid = gamma > 1.0_dp
  end function valid_ideal_gas_gamma

  pure real(dp) function ideal_gas_pressure(conserved, gamma) result(pressure)
    real(dp), intent(in) :: conserved(ncons)
    real(dp), intent(in) :: gamma
    real(dp) :: rho, kinetic_energy_density

    rho = conserved(irho)
    if (rho <= density_floor .or. .not. valid_ideal_gas_gamma(gamma)) then
      pressure = -huge(1.0_dp)
      return
    end if

    kinetic_energy_density = 0.5_dp * &
      (conserved(imx)**2 + conserved(imy)**2 + conserved(imz)**2) / rho
    pressure = (gamma - 1.0_dp) * (conserved(iet) - kinetic_energy_density)
  end function ideal_gas_pressure

  pure real(dp) function ideal_gas_sound_speed(rho, pressure, gamma) result(sound_speed)
    real(dp), intent(in) :: rho, pressure, gamma

    if (rho <= density_floor .or. pressure <= pressure_floor .or. &
        .not. valid_ideal_gas_gamma(gamma)) then
      sound_speed = -huge(1.0_dp)
      return
    end if

    sound_speed = sqrt(gamma * pressure / rho)
  end function ideal_gas_sound_speed

  pure real(dp) function ideal_gas_specific_internal_energy(rho, pressure, gamma) result(e_int)
    real(dp), intent(in) :: rho, pressure, gamma

    if (rho <= density_floor .or. pressure <= pressure_floor .or. &
        .not. valid_ideal_gas_gamma(gamma)) then
      e_int = -huge(1.0_dp)
      return
    end if

    e_int = pressure / ((gamma - 1.0_dp) * rho)
  end function ideal_gas_specific_internal_energy

  pure real(dp) function ideal_gas_total_energy_density( &
      rho, velocity_x, velocity_y, velocity_z, pressure, gamma) result(total_energy_density)
    real(dp), intent(in) :: rho, velocity_x, velocity_y, velocity_z
    real(dp), intent(in) :: pressure, gamma
    real(dp) :: kinetic_energy_density

    if (rho <= density_floor .or. pressure <= pressure_floor .or. &
        .not. valid_ideal_gas_gamma(gamma)) then
      total_energy_density = -huge(1.0_dp)
      return
    end if

    kinetic_energy_density = 0.5_dp * rho * &
      (velocity_x**2 + velocity_y**2 + velocity_z**2)
    total_energy_density = pressure / (gamma - 1.0_dp) + kinetic_energy_density
  end function ideal_gas_total_energy_density

end module eos_ideal_mod
