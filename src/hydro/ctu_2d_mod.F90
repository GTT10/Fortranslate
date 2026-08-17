module ctu_2d_mod
  use precision_mod, only: dp
  use constants_mod, only: density_floor, pressure_floor
  use state_indices_mod, only: &
    ncons, nprim, qrho, qu, qv, qp
  use state_conversion_mod, only: &
    conserved_to_primitive, primitive_to_conserved, state_is_physical
  use eos_ideal_mod, only: ideal_gas_sound_speed
  use slope_limiter_mod, only: limited_slope
  use reconstruction_pelec_plm_mod, only: trace_primitive_characteristics
  use riemann_flux_mod, only: compute_riemann_flux_x
  use directional_flux_mod, only: &
    rotate_primitive_y_to_x, rotate_primitive_x_to_y, compute_riemann_flux_y
  implicit none
  private

  public :: compute_cfl_timestep_2d
  public :: advance_ctu_2d
  public :: all_cells_physical_2d
  public :: apply_transverse_flux_correction

contains

  pure integer function periodic_index(index, extent) result(wrapped)
    integer, intent(in) :: index, extent

    wrapped = 1 + modulo(index - 1, extent)
  end function periodic_index

  subroutine compute_cfl_timestep_2d( &
      conserved, nx, ny, dx, dy, gamma, cfl, dt, ok)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: conserved(ncons, nx, ny)
    real(dp), intent(in) :: dx, dy, gamma, cfl
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok

    real(dp) :: primitive(nprim), sound_speed, local_rate, maximum_rate
    logical :: cell_ok
    integer :: i, j

    maximum_rate = 0.0_dp
    ok = nx > 1 .and. ny > 1 .and. dx > 0.0_dp .and. dy > 0.0_dp .and. &
      cfl > 0.0_dp .and. cfl <= 1.0_dp
    if (.not. ok) then
      dt = 0.0_dp
      return
    end if

    do j = 1, ny
      do i = 1, nx
        call conserved_to_primitive( &
          conserved(:, i, j), gamma, primitive, cell_ok)
        if (.not. cell_ok) then
          ok = .false.
          dt = 0.0_dp
          return
        end if
        sound_speed = ideal_gas_sound_speed( &
          primitive(qrho), primitive(qp), gamma)
        if (sound_speed <= 0.0_dp) then
          ok = .false.
          dt = 0.0_dp
          return
        end if
        local_rate = (abs(primitive(qu)) + sound_speed) / dx + &
          (abs(primitive(qv)) + sound_speed) / dy
        maximum_rate = max(maximum_rate, local_rate)
      end do
    end do

    if (maximum_rate <= 0.0_dp) then
      ok = .false.
      dt = 0.0_dp
      return
    end if

    dt = cfl / maximum_rate
  end subroutine compute_cfl_timestep_2d

  subroutine advance_ctu_2d( &
      conserved, nx, ny, dx, dy, dt, gamma, limiter, riemann_solver, &
      use_transverse_correction, ok, minimum_transverse_theta)
    integer, intent(in) :: nx, ny
    real(dp), intent(inout) :: conserved(ncons, nx, ny)
    real(dp), intent(in) :: dx, dy, dt, gamma
    character(len=*), intent(in) :: limiter, riemann_solver
    logical, intent(in) :: use_transverse_correction
    logical, intent(out) :: ok
    real(dp), intent(out), optional :: minimum_transverse_theta

    real(dp), allocatable :: primitive(:, :, :)
    real(dp), allocatable :: slope_x(:, :, :), slope_y(:, :, :)
    real(dp), allocatable :: x_minus(:, :, :), x_plus(:, :, :)
    real(dp), allocatable :: y_minus(:, :, :), y_plus(:, :, :)
    real(dp), allocatable :: x_left_base(:, :, :), x_right_base(:, :, :)
    real(dp), allocatable :: y_lower_base(:, :, :), y_upper_base(:, :, :)
    real(dp), allocatable :: provisional_x_flux(:, :, :)
    real(dp), allocatable :: provisional_y_flux(:, :, :)
    real(dp), allocatable :: final_x_flux(:, :, :), final_y_flux(:, :, :)
    real(dp), allocatable :: new_state(:, :, :)

    real(dp) :: rotated_center(nprim), rotated_slope(nprim)
    real(dp) :: rotated_minus(nprim), rotated_plus(nprim)
    real(dp) :: x_left(ncons), x_right(ncons)
    real(dp) :: y_lower(ncons), y_upper(ncons)
    real(dp) :: theta_x_left, theta_x_right
    real(dp) :: theta_y_lower, theta_y_upper
    real(dp) :: local_minimum_theta, theta_x, theta_y
    logical :: cell_ok, limiter_ok, trace_ok, face_ok, correction_ok
    integer :: i, j, component, im, ip, jm, jp

    ok = .false.
    if (nx < 4 .or. ny < 4 .or. dx <= 0.0_dp .or. dy <= 0.0_dp .or. &
        dt <= 0.0_dp) then
      if (present(minimum_transverse_theta)) minimum_transverse_theta = 0.0_dp
      return
    end if

    allocate(primitive(nprim, nx, ny))
    allocate(slope_x(nprim, nx, ny), slope_y(nprim, nx, ny))
    allocate(x_minus(nprim, nx, ny), x_plus(nprim, nx, ny))
    allocate(y_minus(nprim, nx, ny), y_plus(nprim, nx, ny))
    allocate(x_left_base(ncons, nx, ny), x_right_base(ncons, nx, ny))
    allocate(y_lower_base(ncons, nx, ny), y_upper_base(ncons, nx, ny))
    allocate(provisional_x_flux(ncons, nx, ny))
    allocate(provisional_y_flux(ncons, nx, ny))
    allocate(final_x_flux(ncons, nx, ny), final_y_flux(ncons, nx, ny))
    allocate(new_state(ncons, nx, ny))

    slope_x = 0.0_dp
    slope_y = 0.0_dp

    do j = 1, ny
      do i = 1, nx
        call conserved_to_primitive( &
          conserved(:, i, j), gamma, primitive(:, i, j), cell_ok)
        if (.not. cell_ok) then
          if (present(minimum_transverse_theta)) minimum_transverse_theta = 0.0_dp
          return
        end if
      end do
    end do

    do j = 1, ny
      jm = periodic_index(j - 1, ny)
      jp = periodic_index(j + 1, ny)
      do i = 1, nx
        im = periodic_index(i - 1, nx)
        ip = periodic_index(i + 1, nx)

        do component = 1, nprim
          call limited_slope( &
            primitive(component, i, j) - primitive(component, im, j), &
            primitive(component, ip, j) - primitive(component, i, j), &
            limiter, slope_x(component, i, j), limiter_ok)
          if (.not. limiter_ok) then
            if (present(minimum_transverse_theta)) minimum_transverse_theta = 0.0_dp
            return
          end if

          call limited_slope( &
            primitive(component, i, j) - primitive(component, i, jm), &
            primitive(component, i, jp) - primitive(component, i, j), &
            limiter, slope_y(component, i, j), limiter_ok)
          if (.not. limiter_ok) then
            if (present(minimum_transverse_theta)) minimum_transverse_theta = 0.0_dp
            return
          end if
        end do

        theta_x = min( &
          positivity_scale_2d( &
            primitive(qrho, i, j), slope_x(qrho, i, j), density_floor), &
          positivity_scale_2d( &
            primitive(qp, i, j), slope_x(qp, i, j), pressure_floor))
        slope_x(:, i, j) = theta_x * slope_x(:, i, j)

        theta_y = min( &
          positivity_scale_2d( &
            primitive(qrho, i, j), slope_y(qrho, i, j), density_floor), &
          positivity_scale_2d( &
            primitive(qp, i, j), slope_y(qp, i, j), pressure_floor))
        slope_y(:, i, j) = theta_y * slope_y(:, i, j)
      end do
    end do

    do j = 1, ny
      do i = 1, nx
        call trace_primitive_characteristics( &
          primitive(:, i, j), slope_x(:, i, j), gamma, dt / dx, &
          x_minus(:, i, j), x_plus(:, i, j), trace_ok)
        if (.not. trace_ok) then
          x_minus(:, i, j) = primitive(:, i, j)
          x_plus(:, i, j) = primitive(:, i, j)
        end if

        call rotate_primitive_y_to_x(primitive(:, i, j), rotated_center)
        call rotate_primitive_y_to_x(slope_y(:, i, j), rotated_slope)
        call trace_primitive_characteristics( &
          rotated_center, rotated_slope, gamma, dt / dy, &
          rotated_minus, rotated_plus, trace_ok)
        if (.not. trace_ok) then
          rotated_minus = rotated_center
          rotated_plus = rotated_center
        end if
        call rotate_primitive_x_to_y(rotated_minus, y_minus(:, i, j))
        call rotate_primitive_x_to_y(rotated_plus, y_plus(:, i, j))
      end do
    end do

    do j = 1, ny
      jp = periodic_index(j + 1, ny)
      do i = 1, nx
        ip = periodic_index(i + 1, nx)

        call primitive_to_conserved( &
          x_plus(:, i, j), gamma, x_left_base(:, i, j), cell_ok)
        if (.not. cell_ok) x_left_base(:, i, j) = conserved(:, i, j)
        call primitive_to_conserved( &
          x_minus(:, ip, j), gamma, x_right_base(:, i, j), cell_ok)
        if (.not. cell_ok) x_right_base(:, i, j) = conserved(:, ip, j)

        call primitive_to_conserved( &
          y_plus(:, i, j), gamma, y_lower_base(:, i, j), cell_ok)
        if (.not. cell_ok) y_lower_base(:, i, j) = conserved(:, i, j)
        call primitive_to_conserved( &
          y_minus(:, i, jp), gamma, y_upper_base(:, i, j), cell_ok)
        if (.not. cell_ok) y_upper_base(:, i, j) = conserved(:, i, jp)

        call compute_riemann_flux_x( &
          x_left_base(:, i, j), x_right_base(:, i, j), gamma, &
          riemann_solver, provisional_x_flux(:, i, j), face_ok)
        if (.not. face_ok) then
          if (present(minimum_transverse_theta)) minimum_transverse_theta = 0.0_dp
          return
        end if

        call compute_riemann_flux_y( &
          y_lower_base(:, i, j), y_upper_base(:, i, j), gamma, &
          riemann_solver, provisional_y_flux(:, i, j), face_ok)
        if (.not. face_ok) then
          if (present(minimum_transverse_theta)) minimum_transverse_theta = 0.0_dp
          return
        end if
      end do
    end do

    local_minimum_theta = 1.0_dp
    do j = 1, ny
      jm = periodic_index(j - 1, ny)
      jp = periodic_index(j + 1, ny)
      do i = 1, nx
        im = periodic_index(i - 1, nx)
        ip = periodic_index(i + 1, nx)

        if (use_transverse_correction) then
          call apply_transverse_flux_correction( &
            x_left_base(:, i, j), provisional_y_flux(:, i, j), &
            provisional_y_flux(:, i, jm), 0.5_dp * dt / dy, gamma, &
            x_left, theta_x_left, correction_ok)
          if (.not. correction_ok) then
            if (present(minimum_transverse_theta)) minimum_transverse_theta = 0.0_dp
            return
          end if

          call apply_transverse_flux_correction( &
            x_right_base(:, i, j), provisional_y_flux(:, ip, j), &
            provisional_y_flux(:, ip, jm), 0.5_dp * dt / dy, gamma, &
            x_right, theta_x_right, correction_ok)
          if (.not. correction_ok) then
            if (present(minimum_transverse_theta)) minimum_transverse_theta = 0.0_dp
            return
          end if

          call apply_transverse_flux_correction( &
            y_lower_base(:, i, j), provisional_x_flux(:, i, j), &
            provisional_x_flux(:, im, j), 0.5_dp * dt / dx, gamma, &
            y_lower, theta_y_lower, correction_ok)
          if (.not. correction_ok) then
            if (present(minimum_transverse_theta)) minimum_transverse_theta = 0.0_dp
            return
          end if

          call apply_transverse_flux_correction( &
            y_upper_base(:, i, j), provisional_x_flux(:, i, jp), &
            provisional_x_flux(:, im, jp), 0.5_dp * dt / dx, gamma, &
            y_upper, theta_y_upper, correction_ok)
          if (.not. correction_ok) then
            if (present(minimum_transverse_theta)) minimum_transverse_theta = 0.0_dp
            return
          end if
        else
          x_left = x_left_base(:, i, j)
          x_right = x_right_base(:, i, j)
          y_lower = y_lower_base(:, i, j)
          y_upper = y_upper_base(:, i, j)
          theta_x_left = 1.0_dp
          theta_x_right = 1.0_dp
          theta_y_lower = 1.0_dp
          theta_y_upper = 1.0_dp
        end if

        local_minimum_theta = min(local_minimum_theta, theta_x_left, &
          theta_x_right, theta_y_lower, theta_y_upper)

        call compute_riemann_flux_x( &
          x_left, x_right, gamma, riemann_solver, final_x_flux(:, i, j), &
          face_ok)
        if (.not. face_ok) then
          if (present(minimum_transverse_theta)) minimum_transverse_theta = 0.0_dp
          return
        end if

        call compute_riemann_flux_y( &
          y_lower, y_upper, gamma, riemann_solver, final_y_flux(:, i, j), &
          face_ok)
        if (.not. face_ok) then
          if (present(minimum_transverse_theta)) minimum_transverse_theta = 0.0_dp
          return
        end if
      end do
    end do

    do j = 1, ny
      jm = periodic_index(j - 1, ny)
      do i = 1, nx
        im = periodic_index(i - 1, nx)
        new_state(:, i, j) = conserved(:, i, j) - &
          dt / dx * (final_x_flux(:, i, j) - final_x_flux(:, im, j)) - &
          dt / dy * (final_y_flux(:, i, j) - final_y_flux(:, i, jm))
      end do
    end do

    if (.not. all_cells_physical_2d(new_state, nx, ny, gamma)) then
      if (present(minimum_transverse_theta)) minimum_transverse_theta = 0.0_dp
      return
    end if

    conserved = new_state
    if (present(minimum_transverse_theta)) then
      minimum_transverse_theta = local_minimum_theta
    end if
    ok = .true.
  end subroutine advance_ctu_2d

  pure subroutine apply_transverse_flux_correction( &
      base_state, flux_high, flux_low, scale, gamma, corrected_state, theta, ok)
    real(dp), intent(in) :: base_state(ncons), flux_high(ncons), flux_low(ncons)
    real(dp), intent(in) :: scale, gamma
    real(dp), intent(out) :: corrected_state(ncons)
    real(dp), intent(out) :: theta
    logical, intent(out) :: ok

    real(dp) :: correction(ncons), trial_state(ncons)
    real(dp) :: lower_theta, upper_theta, midpoint
    integer :: iteration

    corrected_state = base_state
    theta = 0.0_dp
    ok = .false.

    if (scale < 0.0_dp .or. .not. state_is_physical(base_state, gamma)) return

    correction = -scale * (flux_high - flux_low)
    trial_state = base_state + correction
    if (state_is_physical(trial_state, gamma)) then
      corrected_state = trial_state
      theta = 1.0_dp
      ok = .true.
      return
    end if

    lower_theta = 0.0_dp
    upper_theta = 1.0_dp
    do iteration = 1, 60
      midpoint = 0.5_dp * (lower_theta + upper_theta)
      trial_state = base_state + midpoint * correction
      if (state_is_physical(trial_state, gamma)) then
        lower_theta = midpoint
      else
        upper_theta = midpoint
      end if
    end do

    theta = lower_theta
    corrected_state = base_state + theta * correction
    ok = state_is_physical(corrected_state, gamma)
  end subroutine apply_transverse_flux_correction

  pure logical function all_cells_physical_2d( &
      conserved, nx, ny, gamma) result(all_physical)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: conserved(ncons, nx, ny)
    real(dp), intent(in) :: gamma
    integer :: i, j

    all_physical = .true.
    do j = 1, ny
      do i = 1, nx
        if (.not. state_is_physical(conserved(:, i, j), gamma)) then
          all_physical = .false.
          return
        end if
      end do
    end do
  end function all_cells_physical_2d

  pure real(dp) function positivity_scale_2d( &
      center, slope, lower_bound) result(theta)
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
  end function positivity_scale_2d

end module ctu_2d_mod
