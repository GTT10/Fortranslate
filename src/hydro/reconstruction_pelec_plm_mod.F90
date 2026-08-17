module reconstruction_pelec_plm_mod
  use precision_mod, only: dp
  use constants_mod, only: density_floor, pressure_floor
  use state_indices_mod, only: ncons, nprim, qrho, qu, qv, qw, qp
  use state_conversion_mod, only: conserved_to_primitive, primitive_to_conserved
  use eos_ideal_mod, only: ideal_gas_sound_speed
  use slope_limiter_mod, only: limited_slope
  implicit none
  private

  integer, parameter, public :: n_characteristic_waves = 5
  integer, parameter, public :: wave_minus = 1
  integer, parameter, public :: wave_contact = 2
  integer, parameter, public :: wave_shear_y = 3
  integer, parameter, public :: wave_shear_z = 4
  integer, parameter, public :: wave_plus = 5

  public :: reconstruct_pelec_plm_faces
  public :: primitive_slope_to_characteristics
  public :: characteristics_to_primitive_slope
  public :: trace_primitive_characteristics
  public :: pelec_limited_slope
  public :: pelec_flattening_coefficient

contains

  pure subroutine primitive_slope_to_characteristics( &
      center, slope, gamma, characteristic, ok)
    real(dp), intent(in) :: center(nprim), slope(nprim)
    real(dp), intent(in) :: gamma
    real(dp), intent(out) :: characteristic(n_characteristic_waves)
    logical, intent(out) :: ok

    real(dp) :: rho, pressure, sound_speed, sound_speed_squared

    characteristic = 0.0_dp
    ok = .false.

    rho = center(qrho)
    pressure = center(qp)
    if (rho <= density_floor .or. pressure <= pressure_floor) return

    sound_speed = ideal_gas_sound_speed(rho, pressure, gamma)
    if (sound_speed <= 0.0_dp) return
    sound_speed_squared = sound_speed * sound_speed

    characteristic(wave_minus) = 0.5_dp * &
      (slope(qp) / (rho * sound_speed) - slope(qu)) * rho / sound_speed
    characteristic(wave_plus) = 0.5_dp * &
      (slope(qp) / (rho * sound_speed) + slope(qu)) * rho / sound_speed
    characteristic(wave_contact) = slope(qrho) - slope(qp) / sound_speed_squared
    characteristic(wave_shear_y) = slope(qv)
    characteristic(wave_shear_z) = slope(qw)
    ok = .true.
  end subroutine primitive_slope_to_characteristics

  pure subroutine characteristics_to_primitive_slope( &
      center, characteristic, gamma, slope, ok)
    real(dp), intent(in) :: center(nprim)
    real(dp), intent(in) :: characteristic(n_characteristic_waves)
    real(dp), intent(in) :: gamma
    real(dp), intent(out) :: slope(nprim)
    logical, intent(out) :: ok

    real(dp) :: rho, pressure, sound_speed, sound_speed_squared

    slope = 0.0_dp
    ok = .false.

    rho = center(qrho)
    pressure = center(qp)
    if (rho <= density_floor .or. pressure <= pressure_floor) return

    sound_speed = ideal_gas_sound_speed(rho, pressure, gamma)
    if (sound_speed <= 0.0_dp) return
    sound_speed_squared = sound_speed * sound_speed

    slope(qrho) = characteristic(wave_minus) + &
      characteristic(wave_contact) + characteristic(wave_plus)
    slope(qu) = (characteristic(wave_plus) - characteristic(wave_minus)) * &
      sound_speed / rho
    slope(qv) = characteristic(wave_shear_y)
    slope(qw) = characteristic(wave_shear_z)
    slope(qp) = (characteristic(wave_minus) + characteristic(wave_plus)) * &
      sound_speed_squared
    ok = .true.
  end subroutine characteristics_to_primitive_slope

  pure subroutine trace_primitive_characteristics( &
      center, slope, gamma, dtdx, left_state, right_state, ok)
    real(dp), intent(in) :: center(nprim), slope(nprim)
    real(dp), intent(in) :: gamma, dtdx
    real(dp), intent(out) :: left_state(nprim), right_state(nprim)
    logical, intent(out) :: ok

    real(dp) :: characteristic(n_characteristic_waves)
    real(dp) :: rho, velocity, pressure, sound_speed, sound_speed_squared
    real(dp) :: wave_speed_minus, wave_speed_contact, wave_speed_plus
    real(dp) :: left_factor, right_factor
    real(dp) :: acoustic_plus, acoustic_minus
    real(dp) :: contact_left, contact_right
    real(dp) :: shear_y_left, shear_y_right
    real(dp) :: shear_z_left, shear_z_right
    logical :: projection_ok

    left_state = 0.0_dp
    right_state = 0.0_dp
    ok = .false.
    if (dtdx < 0.0_dp) return

    call primitive_slope_to_characteristics( &
      center, slope, gamma, characteristic, projection_ok)
    if (.not. projection_ok) return

    rho = center(qrho)
    velocity = center(qu)
    pressure = center(qp)
    sound_speed = ideal_gas_sound_speed(rho, pressure, gamma)
    if (sound_speed <= 0.0_dp) return
    sound_speed_squared = sound_speed * sound_speed

    wave_speed_minus = velocity - sound_speed
    wave_speed_contact = velocity
    wave_speed_plus = velocity + sound_speed

    left_factor = 0.5_dp * &
      (1.0_dp + dtdx * min(wave_speed_minus, 0.0_dp))
    left_state = center - left_factor * slope

    acoustic_plus = 0.25_dp * dtdx * &
      (wave_speed_minus - wave_speed_plus) * &
      (1.0_dp - sign(1.0_dp, wave_speed_plus)) * &
      characteristic(wave_plus)
    contact_left = 0.25_dp * dtdx * &
      (wave_speed_minus - wave_speed_contact) * &
      (1.0_dp - sign(1.0_dp, wave_speed_contact)) * &
      characteristic(wave_contact)
    shear_y_left = 0.25_dp * dtdx * &
      (wave_speed_minus - wave_speed_contact) * &
      (1.0_dp - sign(1.0_dp, wave_speed_contact)) * &
      characteristic(wave_shear_y)
    shear_z_left = 0.25_dp * dtdx * &
      (wave_speed_minus - wave_speed_contact) * &
      (1.0_dp - sign(1.0_dp, wave_speed_contact)) * &
      characteristic(wave_shear_z)

    left_state(qrho) = left_state(qrho) + acoustic_plus + contact_left
    left_state(qu) = left_state(qu) + acoustic_plus * sound_speed / rho
    left_state(qv) = left_state(qv) + shear_y_left
    left_state(qw) = left_state(qw) + shear_z_left
    left_state(qp) = left_state(qp) + acoustic_plus * sound_speed_squared

    right_factor = 0.5_dp * &
      (1.0_dp - dtdx * max(wave_speed_plus, 0.0_dp))
    right_state = center + right_factor * slope

    acoustic_minus = 0.25_dp * dtdx * &
      (wave_speed_plus - wave_speed_minus) * &
      (1.0_dp + sign(1.0_dp, wave_speed_minus)) * &
      characteristic(wave_minus)
    contact_right = 0.25_dp * dtdx * &
      (wave_speed_plus - wave_speed_contact) * &
      (1.0_dp + sign(1.0_dp, wave_speed_contact)) * &
      characteristic(wave_contact)
    shear_y_right = 0.25_dp * dtdx * &
      (wave_speed_plus - wave_speed_contact) * &
      (1.0_dp + sign(1.0_dp, wave_speed_contact)) * &
      characteristic(wave_shear_y)
    shear_z_right = 0.25_dp * dtdx * &
      (wave_speed_plus - wave_speed_contact) * &
      (1.0_dp + sign(1.0_dp, wave_speed_contact)) * &
      characteristic(wave_shear_z)

    right_state(qrho) = right_state(qrho) + acoustic_minus + contact_right
    right_state(qu) = right_state(qu) - acoustic_minus * sound_speed / rho
    right_state(qv) = right_state(qv) + shear_y_right
    right_state(qw) = right_state(qw) + shear_z_right
    right_state(qp) = right_state(qp) + acoustic_minus * sound_speed_squared

    if (left_state(qrho) <= density_floor .or. &
        right_state(qrho) <= density_floor) return
    if (left_state(qp) <= pressure_floor .or. &
        right_state(qp) <= pressure_floor) return

    ok = .true.
  end subroutine trace_primitive_characteristics

  pure subroutine pelec_limited_slope( &
      qm2, qm, qc, qp, qp2, flat, order, slope, ok)
    real(dp), intent(in) :: qm2, qm, qc, qp, qp2, flat
    integer, intent(in) :: order
    real(dp), intent(out) :: slope
    logical, intent(out) :: ok

    real(dp) :: dlft, drgt, dcen, dfm, dfp
    real(dp) :: dlim, dsgn, dtemp

    slope = 0.0_dp
    ok = .false.
    if (order /= 2 .and. order /= 4) return
    if (flat < 0.0_dp .or. flat > 1.0_dp) return

    dfm = 0.0_dp
    dfp = 0.0_dp
    if (order == 4) then
      dlft = qm - qm2
      drgt = qc - qm
      dcen = 0.5_dp * (dlft + drgt)
      dsgn = sign(1.0_dp, dcen)
      if (dlft * drgt >= 0.0_dp) then
        dlim = 2.0_dp * min(abs(dlft), abs(drgt))
      else
        dlim = 0.0_dp
      end if
      dfm = dsgn * min(dlim, abs(dcen))

      dlft = qp - qc
      drgt = qp2 - qp
      dcen = 0.5_dp * (dlft + drgt)
      dsgn = sign(1.0_dp, dcen)
      if (dlft * drgt >= 0.0_dp) then
        dlim = 2.0_dp * min(abs(dlft), abs(drgt))
      else
        dlim = 0.0_dp
      end if
      dfp = dsgn * min(dlim, abs(dcen))
    end if

    dlft = qc - qm
    drgt = qp - qc
    dcen = 0.5_dp * (dlft + drgt)
    dsgn = sign(1.0_dp, dcen)
    if (dlft * drgt >= 0.0_dp) then
      dlim = 2.0_dp * min(abs(dlft), abs(drgt))
    else
      dlim = 0.0_dp
    end if

    dtemp = (4.0_dp / 3.0_dp) * dcen - (dfp + dfm) / 6.0_dp
    slope = flat * dsgn * min(dlim, abs(dtemp))
    ok = .true.
  end subroutine pelec_limited_slope

  pure real(dp) function pelec_flattening_coefficient( &
      primitive, nx, cell, boundary_condition) result(flat)
    integer, intent(in) :: nx, cell
    real(dp), intent(in) :: primitive(:, -2:)
    character(len=*), intent(in) :: boundary_condition

    real(dp), parameter :: small_pressure = 1.0e-200_dp
    real(dp), parameter :: shock_threshold = 0.33_dp
    real(dp), parameter :: zcut1 = 0.75_dp
    real(dp), parameter :: zcut2 = 0.85_dp
    real(dp), parameter :: inverse_zcut_width = 1.0_dp / (zcut2 - zcut1)

    real(dp) :: pressure_minus, pressure_plus
    real(dp) :: pressure_minus2, pressure_plus2
    real(dp) :: velocity_minus, velocity_plus
    real(dp) :: pressure_jump, denominator, zeta, z
    real(dp) :: compression_test, minimum_pressure, chi
    real(dp) :: shifted_pressure_minus, shifted_pressure_plus
    real(dp) :: shifted_pressure_minus2, shifted_pressure_plus2
    real(dp) :: shifted_velocity_minus, shifted_velocity_plus
    real(dp) :: zeta_shifted, z_shifted, chi_shifted
    integer :: shift_direction

    flat = 1.0_dp
    if (cell < 1 .or. cell > nx) return

    select case (trim(boundary_condition))
    case ("outflow")
      if (cell < 3 .or. cell > nx - 2) return
    case ("periodic")
      continue
    case default
      return
    end select

    pressure_minus = primitive(qp, cell - 1)
    pressure_plus = primitive(qp, cell + 1)
    velocity_minus = primitive(qu, cell - 1)
    velocity_plus = primitive(qu, cell + 1)
    pressure_minus2 = primitive(qp, cell - 2)
    pressure_plus2 = primitive(qp, cell + 2)

    pressure_jump = pressure_plus - pressure_minus
    if (pressure_jump > 0.0_dp) then
      shift_direction = 1
    else
      shift_direction = -1
    end if

    denominator = max(small_pressure, abs(pressure_plus2 - pressure_minus2))
    zeta = abs(pressure_jump) / denominator
    z = min(1.0_dp, max(0.0_dp, inverse_zcut_width * (zeta - zcut1)))

    if (velocity_minus - velocity_plus >= 0.0_dp) then
      compression_test = 1.0_dp
    else
      compression_test = 0.0_dp
    end if
    minimum_pressure = max( &
      small_pressure, min(pressure_plus, pressure_minus))
    if (abs(pressure_jump) / minimum_pressure > shock_threshold) then
      chi = compression_test
    else
      chi = 0.0_dp
    end if

    shifted_pressure_plus = primitive(qp, cell + 1 - shift_direction)
    shifted_pressure_minus = primitive(qp, cell - 1 - shift_direction)
    shifted_velocity_plus = primitive(qu, cell + 1 - shift_direction)
    shifted_velocity_minus = primitive(qu, cell - 1 - shift_direction)
    shifted_pressure_plus2 = primitive(qp, cell + 2 - shift_direction)
    shifted_pressure_minus2 = primitive(qp, cell - 2 - shift_direction)

    pressure_jump = shifted_pressure_plus - shifted_pressure_minus
    denominator = max( &
      small_pressure, abs(shifted_pressure_plus2 - shifted_pressure_minus2))
    zeta_shifted = abs(pressure_jump) / denominator
    z_shifted = min(1.0_dp, max(0.0_dp, &
      inverse_zcut_width * (zeta_shifted - zcut1)))

    if (shifted_velocity_minus - shifted_velocity_plus >= 0.0_dp) then
      compression_test = 1.0_dp
    else
      compression_test = 0.0_dp
    end if
    minimum_pressure = max( &
      small_pressure, min(shifted_pressure_plus, shifted_pressure_minus))
    if (abs(pressure_jump) / minimum_pressure > shock_threshold) then
      chi_shifted = compression_test
    else
      chi_shifted = 0.0_dp
    end if

    flat = 1.0_dp - max(chi_shifted * z_shifted, chi * z)
    flat = min(1.0_dp, max(0.0_dp, flat))
  end function pelec_flattening_coefficient

  subroutine reconstruct_pelec_plm_faces( &
      conserved, nx, gamma, limiter, boundary_condition, dtdx, &
      left_faces, right_faces, ok, slope_order, use_flattening)
    integer, intent(in) :: nx
    real(dp), intent(in) :: conserved(ncons, 0:nx + 1)
    real(dp), intent(in) :: gamma, dtdx
    character(len=*), intent(in) :: limiter, boundary_condition
    real(dp), intent(out) :: left_faces(ncons, 0:nx)
    real(dp), intent(out) :: right_faces(ncons, 0:nx)
    logical, intent(out) :: ok
    integer, intent(in), optional :: slope_order
    logical, intent(in), optional :: use_flattening

    real(dp), allocatable :: primitive(:, :), slopes(:, :)
    real(dp), allocatable :: cell_left(:, :), cell_right(:, :)
    real(dp) :: theta, flat
    integer :: order
    logical :: flattening_enabled
    logical :: cell_ok, limiter_ok, trace_ok, left_ok, right_ok, slope_ok
    integer :: i, component

    allocate(primitive(nprim, -2:nx + 3))
    allocate(slopes(nprim, -2:nx + 3))
    allocate(cell_left(nprim, 0:nx + 1))
    allocate(cell_right(nprim, 0:nx + 1))

    left_faces = 0.0_dp
    right_faces = 0.0_dp
    primitive = 0.0_dp
    slopes = 0.0_dp
    ok = .false.
    if (dtdx < 0.0_dp) return

    order = 2
    flattening_enabled = .false.
    if (present(slope_order)) order = slope_order
    if (present(use_flattening)) flattening_enabled = use_flattening
    if (order /= 2 .and. order /= 4) return

    do i = 1, nx
      call conserved_to_primitive(conserved(:, i), gamma, primitive(:, i), cell_ok)
      if (.not. cell_ok) return
    end do
    call fill_primitive_ghosts(primitive, nx, boundary_condition, cell_ok)
    if (.not. cell_ok) return

    do i = 1, nx
      flat = 1.0_dp
      if (flattening_enabled) then
        flat = pelec_flattening_coefficient( &
          primitive, nx, i, boundary_condition)
      end if

      do component = 1, nprim
        if (order == 4) then
          call pelec_limited_slope( &
            primitive(component, i - 2), primitive(component, i - 1), &
            primitive(component, i), primitive(component, i + 1), &
            primitive(component, i + 2), flat, order, &
            slopes(component, i), slope_ok)
          if (.not. slope_ok) return
        else
          call limited_slope( &
            primitive(component, i) - primitive(component, i - 1), &
            primitive(component, i + 1) - primitive(component, i), &
            limiter, slopes(component, i), limiter_ok)
          if (.not. limiter_ok) return
          slopes(component, i) = flat * slopes(component, i)
        end if
      end do

      theta = min( &
        positivity_scale(primitive(qrho, i), slopes(qrho, i), density_floor), &
        positivity_scale(primitive(qp, i), slopes(qp, i), pressure_floor))
      slopes(:, i) = theta * slopes(:, i)
    end do

    call fill_slope_ghosts(slopes, nx, boundary_condition, cell_ok)
    if (.not. cell_ok) return

    do i = 0, nx + 1
      call trace_primitive_characteristics( &
        primitive(:, i), slopes(:, i), gamma, dtdx, &
        cell_left(:, i), cell_right(:, i), trace_ok)
      if (.not. trace_ok) then
        cell_left(:, i) = primitive(:, i)
        cell_right(:, i) = primitive(:, i)
      end if
    end do

    do i = 0, nx
      call primitive_to_conserved( &
        cell_right(:, i), gamma, left_faces(:, i), left_ok)
      if (.not. left_ok) left_faces(:, i) = conserved(:, i)

      call primitive_to_conserved( &
        cell_left(:, i + 1), gamma, right_faces(:, i), right_ok)
      if (.not. right_ok) right_faces(:, i) = conserved(:, i + 1)
    end do

    ok = .true.
  end subroutine reconstruct_pelec_plm_faces

  pure subroutine fill_primitive_ghosts( &
      primitive, nx, boundary_condition, ok)
    integer, intent(in) :: nx
    real(dp), intent(inout) :: primitive(nprim, -2:nx + 3)
    character(len=*), intent(in) :: boundary_condition
    logical, intent(out) :: ok
    integer :: i, wrapped

    select case (trim(boundary_condition))
    case ("outflow")
      do i = -2, 0
        primitive(:, i) = primitive(:, 1)
      end do
      do i = nx + 1, nx + 3
        primitive(:, i) = primitive(:, nx)
      end do
      ok = .true.
    case ("periodic")
      do i = -2, 0
        wrapped = 1 + modulo(i - 1, nx)
        primitive(:, i) = primitive(:, wrapped)
      end do
      do i = nx + 1, nx + 3
        wrapped = 1 + modulo(i - 1, nx)
        primitive(:, i) = primitive(:, wrapped)
      end do
      ok = .true.
    case default
      ok = .false.
    end select
  end subroutine fill_primitive_ghosts

  pure subroutine fill_slope_ghosts(slopes, nx, boundary_condition, ok)
    integer, intent(in) :: nx
    real(dp), intent(inout) :: slopes(nprim, -2:nx + 3)
    character(len=*), intent(in) :: boundary_condition
    logical, intent(out) :: ok
    integer :: i, wrapped

    select case (trim(boundary_condition))
    case ("outflow")
      slopes(:, -2:0) = 0.0_dp
      slopes(:, nx + 1:nx + 3) = 0.0_dp
      slopes(:, 1) = 0.0_dp
      slopes(:, nx) = 0.0_dp
      ok = .true.
    case ("periodic")
      do i = -2, 0
        wrapped = 1 + modulo(i - 1, nx)
        slopes(:, i) = slopes(:, wrapped)
      end do
      do i = nx + 1, nx + 3
        wrapped = 1 + modulo(i - 1, nx)
        slopes(:, i) = slopes(:, wrapped)
      end do
      ok = .true.
    case default
      ok = .false.
    end select
  end subroutine fill_slope_ghosts

  pure real(dp) function positivity_scale(center, slope, lower_bound) result(theta)
    real(dp), intent(in) :: center, slope, lower_bound
    real(dp) :: slope_magnitude

    slope_magnitude = abs(slope)
    if (slope_magnitude <= tiny(1.0_dp)) then
      theta = 1.0_dp
    else if (center - 0.5_dp * slope_magnitude > lower_bound) then
      theta = 1.0_dp
    else
      theta = max(0.0_dp, min(1.0_dp, &
        2.0_dp * (center - lower_bound) / slope_magnitude))
    end if
  end function positivity_scale

end module reconstruction_pelec_plm_mod
