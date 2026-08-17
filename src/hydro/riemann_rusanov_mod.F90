module riemann_rusanov_mod
  use precision_mod, only: dp
  use state_indices_mod, only: &
    irho, imx, imy, imz, iet, ncons, qrho, qu, qv, qw, qp, nprim
  use state_conversion_mod, only: conserved_to_primitive
  use eos_ideal_mod, only: ideal_gas_sound_speed
  implicit none
  private

  public :: euler_physical_flux_x
  public :: rusanov_flux_x

contains

  pure subroutine euler_physical_flux_x(conserved, gamma, flux, ok)
    real(dp), intent(in) :: conserved(ncons)
    real(dp), intent(in) :: gamma
    real(dp), intent(out) :: flux(ncons)
    logical, intent(out) :: ok
    real(dp) :: primitive(nprim)
    real(dp) :: rho, velocity_x, velocity_y, velocity_z, pressure

    call conserved_to_primitive(conserved, gamma, primitive, ok)
    if (.not. ok) then
      flux = 0.0_dp
      return
    end if

    rho = primitive(qrho)
    velocity_x = primitive(qu)
    velocity_y = primitive(qv)
    velocity_z = primitive(qw)
    pressure = primitive(qp)

    flux(irho) = rho * velocity_x
    flux(imx) = rho * velocity_x**2 + pressure
    flux(imy) = rho * velocity_x * velocity_y
    flux(imz) = rho * velocity_x * velocity_z
    flux(iet) = (conserved(iet) + pressure) * velocity_x
  end subroutine euler_physical_flux_x

  pure subroutine rusanov_flux_x(left_state, right_state, gamma, flux, ok)
    real(dp), intent(in) :: left_state(ncons), right_state(ncons)
    real(dp), intent(in) :: gamma
    real(dp), intent(out) :: flux(ncons)
    logical, intent(out) :: ok
    real(dp) :: left_flux(ncons), right_flux(ncons)
    real(dp) :: left_primitive(nprim), right_primitive(nprim)
    real(dp) :: left_sound_speed, right_sound_speed, max_wave_speed
    logical :: left_ok, right_ok, left_flux_ok, right_flux_ok

    call conserved_to_primitive(left_state, gamma, left_primitive, left_ok)
    call conserved_to_primitive(right_state, gamma, right_primitive, right_ok)
    call euler_physical_flux_x(left_state, gamma, left_flux, left_flux_ok)
    call euler_physical_flux_x(right_state, gamma, right_flux, right_flux_ok)

    ok = left_ok .and. right_ok .and. left_flux_ok .and. right_flux_ok
    if (.not. ok) then
      flux = 0.0_dp
      return
    end if

    left_sound_speed = ideal_gas_sound_speed( &
      left_primitive(qrho), left_primitive(qp), gamma)
    right_sound_speed = ideal_gas_sound_speed( &
      right_primitive(qrho), right_primitive(qp), gamma)
    max_wave_speed = max( &
      abs(left_primitive(qu)) + left_sound_speed, &
      abs(right_primitive(qu)) + right_sound_speed)

    flux = 0.5_dp * (left_flux + right_flux) - &
      0.5_dp * max_wave_speed * (right_state - left_state)
  end subroutine rusanov_flux_x

end module riemann_rusanov_mod
