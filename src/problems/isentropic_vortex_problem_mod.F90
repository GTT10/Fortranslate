module isentropic_vortex_problem_mod
  use precision_mod, only: dp
  use state_indices_mod, only: &
    ncons, nprim, qrho, qu, qv, qw, qp
  use state_conversion_mod, only: primitive_to_conserved
  use simulation_config_2d_mod, only: isentropic_vortex_config
  implicit none
  private

  public :: initialize_isentropic_vortex
  public :: isentropic_vortex_primitive

contains

  subroutine initialize_isentropic_vortex( &
      x, y, nx, ny, x_min, x_max, y_min, y_max, gamma, vortex, &
      conserved, ok)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: x(nx), y(ny)
    real(dp), intent(in) :: x_min, x_max, y_min, y_max, gamma
    type(isentropic_vortex_config), intent(in) :: vortex
    real(dp), intent(out) :: conserved(ncons, nx, ny)
    logical, intent(out) :: ok

    real(dp) :: primitive(nprim)
    logical :: cell_ok
    integer :: i, j

    conserved = 0.0_dp
    ok = .true.

    do j = 1, ny
      do i = 1, nx
        call isentropic_vortex_primitive( &
          x(i), y(j), 0.0_dp, x_min, x_max, y_min, y_max, gamma, &
          vortex, primitive, cell_ok)
        if (.not. cell_ok) then
          ok = .false.
          return
        end if
        call primitive_to_conserved( &
          primitive, gamma, conserved(:, i, j), cell_ok)
        if (.not. cell_ok) then
          ok = .false.
          return
        end if
      end do
    end do
  end subroutine initialize_isentropic_vortex

  pure subroutine isentropic_vortex_primitive( &
      x, y, time, x_min, x_max, y_min, y_max, gamma, vortex, primitive, ok)
    real(dp), intent(in) :: x, y, time
    real(dp), intent(in) :: x_min, x_max, y_min, y_max, gamma
    type(isentropic_vortex_config), intent(in) :: vortex
    real(dp), intent(out) :: primitive(nprim)
    logical, intent(out) :: ok

    real(dp) :: domain_x, domain_y, center_x, center_y
    real(dp) :: displacement_x, displacement_y, radius_squared
    real(dp) :: pi, velocity_factor, temperature_perturbation
    real(dp) :: base_temperature, temperature, temperature_ratio

    primitive = 0.0_dp
    ok = .false.

    domain_x = x_max - x_min
    domain_y = y_max - y_min
    if (domain_x <= 0.0_dp .or. domain_y <= 0.0_dp .or. gamma <= 1.0_dp) return
    if (vortex%base_density <= 0.0_dp .or. vortex%base_pressure <= 0.0_dp) return

    center_x = x_min + modulo( &
      vortex%center_x + vortex%base_velocity_x * time - x_min, domain_x)
    center_y = y_min + modulo( &
      vortex%center_y + vortex%base_velocity_y * time - y_min, domain_y)

    displacement_x = periodic_displacement(x - center_x, domain_x)
    displacement_y = periodic_displacement(y - center_y, domain_y)
    radius_squared = displacement_x**2 + displacement_y**2

    pi = acos(-1.0_dp)
    velocity_factor = vortex%strength / (2.0_dp * pi) * &
      exp(0.5_dp * (1.0_dp - radius_squared))
    temperature_perturbation = -(gamma - 1.0_dp) * vortex%strength**2 / &
      (8.0_dp * gamma * pi**2) * exp(1.0_dp - radius_squared)

    base_temperature = vortex%base_pressure / vortex%base_density
    temperature = base_temperature + temperature_perturbation
    if (temperature <= 0.0_dp) return
    temperature_ratio = temperature / base_temperature

    primitive(qrho) = vortex%base_density * &
      temperature_ratio**(1.0_dp / (gamma - 1.0_dp))
    primitive(qu) = vortex%base_velocity_x - &
      velocity_factor * displacement_y
    primitive(qv) = vortex%base_velocity_y + &
      velocity_factor * displacement_x
    primitive(qw) = 0.0_dp
    primitive(qp) = vortex%base_pressure * &
      temperature_ratio**(gamma / (gamma - 1.0_dp))

    ok = primitive(qrho) > 0.0_dp .and. primitive(qp) > 0.0_dp
  end subroutine isentropic_vortex_primitive

  pure real(dp) function periodic_displacement(delta, period) result(wrapped)
    real(dp), intent(in) :: delta, period

    wrapped = modulo(delta + 0.5_dp * period, period) - 0.5_dp * period
  end function periodic_displacement

end module isentropic_vortex_problem_mod
