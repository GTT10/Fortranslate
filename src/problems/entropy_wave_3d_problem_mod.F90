module entropy_wave_3d_problem_mod
  use precision_mod, only: dp
  use state_indices_mod, only: &
    ncons, nprim, qrho, qu, qv, qw, qp
  use state_conversion_mod, only: primitive_to_conserved
  use simulation_config_3d_mod, only: entropy_wave_3d_config
  implicit none
  private

  public :: initialize_entropy_wave_3d
  public :: entropy_wave_3d_primitive

contains

  subroutine initialize_entropy_wave_3d( &
      x, y, z, nx, ny, nz, x_min, x_max, y_min, y_max, z_min, z_max, &
      gamma, wave, conserved, ok)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: x(nx), y(ny), z(nz)
    real(dp), intent(in) :: x_min, x_max, y_min, y_max, z_min, z_max
    real(dp), intent(in) :: gamma
    type(entropy_wave_3d_config), intent(in) :: wave
    real(dp), intent(out) :: conserved(ncons, nx, ny, nz)
    logical, intent(out) :: ok

    real(dp) :: primitive(nprim)
    logical :: cell_ok
    integer :: i, j, k

    conserved = 0.0_dp
    ok = .true.
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          call entropy_wave_3d_primitive( &
            x(i), y(j), z(k), 0.0_dp, &
            x_min, x_max, y_min, y_max, z_min, z_max, &
            wave, primitive, cell_ok)
          if (.not. cell_ok) then
            ok = .false.
            return
          end if
          call primitive_to_conserved( &
            primitive, gamma, conserved(:, i, j, k), cell_ok)
          if (.not. cell_ok) then
            ok = .false.
            return
          end if
        end do
      end do
    end do
  end subroutine initialize_entropy_wave_3d

  pure subroutine entropy_wave_3d_primitive( &
      x, y, z, time, x_min, x_max, y_min, y_max, z_min, z_max, &
      wave, primitive, ok)
    real(dp), intent(in) :: x, y, z, time
    real(dp), intent(in) :: x_min, x_max, y_min, y_max, z_min, z_max
    type(entropy_wave_3d_config), intent(in) :: wave
    real(dp), intent(out) :: primitive(nprim)
    logical, intent(out) :: ok

    real(dp) :: length_x, length_y, length_z, phase, phase_speed, pi

    primitive = 0.0_dp
    ok = .false.
    length_x = x_max - x_min
    length_y = y_max - y_min
    length_z = z_max - z_min
    if (length_x <= 0.0_dp .or. length_y <= 0.0_dp .or. &
        length_z <= 0.0_dp) return
    if (wave%base_density <= 0.0_dp .or. wave%base_pressure <= 0.0_dp) return
    if (abs(wave%density_amplitude) >= wave%base_density) return

    phase_speed = real(wave%wave_number_x, dp) * wave%velocity_x / length_x + &
      real(wave%wave_number_y, dp) * wave%velocity_y / length_y + &
      real(wave%wave_number_z, dp) * wave%velocity_z / length_z
    pi = acos(-1.0_dp)
    phase = 2.0_dp * pi * ( &
      real(wave%wave_number_x, dp) * (x - x_min) / length_x + &
      real(wave%wave_number_y, dp) * (y - y_min) / length_y + &
      real(wave%wave_number_z, dp) * (z - z_min) / length_z - &
      time * phase_speed)

    primitive(qrho) = wave%base_density + &
      wave%density_amplitude * sin(phase)
    primitive(qu) = wave%velocity_x
    primitive(qv) = wave%velocity_y
    primitive(qw) = wave%velocity_z
    primitive(qp) = wave%base_pressure
    ok = primitive(qrho) > 0.0_dp
  end subroutine entropy_wave_3d_primitive

end module entropy_wave_3d_problem_mod
