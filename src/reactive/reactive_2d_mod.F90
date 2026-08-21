module reactive_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use constants_mod, only: density_floor, pressure_floor
  use state_indices_mod, only: irho, imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: gas_transport_species
  use reactive_transport_2d_mod, only: &
    reactive_transport_timestep_2d, advance_reactive_transport_2d
  use constant_volume_reactor_mod, only: advance_constant_volume_adaptive
  use simulation_config_reactive_2d_mod, only: reactive_2d_config
  use reactive_boundary_2d_mod, only: &
    reactive_boundary_set_2d, initialize_periodic_boundary_set_2d, &
    build_reactive_boundary_set_2d, sample_reactive_primitive_2d, &
    reactive_boundary_is_periodic, reactive_boundary_is_wall
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_mass_fraction_component, reactive_primitive_to_conserved, &
    reactive_conserved_to_primitive, reactive_riemann_flux_x, &
    characteristic_limited_slope, trace_reactive_characteristics, &
    reactive_ppm_reconstruct_five, reactive_ppm_integrate_profile, &
    reactive_ppm_flattening_coefficient, &
    reactive_ppm_contact_steepening_factor, &
    reactive_ppm_apply_contact_steepening, &
    build_characteristic_ppm_states
  implicit none
  private

  integer, parameter :: max_chemistry_substeps = 100000
  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter :: contact_steepening_cap_2d = 0.5_dp

  public :: reactive_riemann_flux_y
  public :: compute_reactive_cfl_timestep_2d
  public :: apply_reactive_transverse_correction_2d
  public :: reconstruct_reactive_characteristic_ppm_cell_2d
  public :: advance_reactive_hydro_2d
  public :: advance_reactive_chemistry_2d
  public :: advance_reactive_strang_2d
  public :: initialize_reactive_2d
  public :: simulate_reactive_2d
  public :: reactive_integrals_2d
  public :: reactive_extrema_2d
  public :: reactive_diagonal_wave_density
  public :: reactive_diagonal_composition_wave_exact
  public :: write_reactive_2d_csv

contains

  pure integer function periodic_index(index, extent) result(wrapped)
    integer, intent(in) :: index, extent
    wrapped = 1 + modulo(index - 1, extent)
  end function periodic_index

  pure subroutine rotate_primitive_y_to_x(input, output)
    real(dp), intent(in) :: input(:)
    real(dp), intent(out) :: output(:)
    output = input
    if (size(input) /= size(output) .or. size(input) < 5) return
    output(2) = input(3)
    output(3) = input(2)
  end subroutine rotate_primitive_y_to_x

  pure subroutine rotate_primitive_x_to_y(input, output)
    real(dp), intent(in) :: input(:)
    real(dp), intent(out) :: output(:)
    call rotate_primitive_y_to_x(input, output)
  end subroutine rotate_primitive_x_to_y

  pure subroutine rotate_conserved_y_to_x(input, output)
    real(dp), intent(in) :: input(:)
    real(dp), intent(out) :: output(:)
    output = input
    if (size(input) /= size(output) .or. size(input) < iet) return
    output(imx) = input(imy)
    output(imy) = input(imx)
  end subroutine rotate_conserved_y_to_x

  pure subroutine rotate_flux_x_to_y(input, output)
    real(dp), intent(in) :: input(:)
    real(dp), intent(out) :: output(:)
    call rotate_conserved_y_to_x(input, output)
  end subroutine rotate_flux_x_to_y

  subroutine reactive_riemann_flux_y( &
      species, lower_state, upper_state, lower_temperature_guess, &
      upper_temperature_guess, solver, flux, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: lower_state(:), upper_state(:)
    real(dp), intent(in) :: lower_temperature_guess, upper_temperature_guess
    character(len=*), intent(in) :: solver
    real(dp), intent(out) :: flux(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: lower_rotated(:), upper_rotated(:), flux_rotated(:)
    integer :: nvar

    flux = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (size(lower_state) /= nvar .or. size(upper_state) /= nvar .or. &
        size(flux) /= nvar) return
    allocate(lower_rotated(nvar), upper_rotated(nvar), flux_rotated(nvar))
    call rotate_conserved_y_to_x(lower_state, lower_rotated)
    call rotate_conserved_y_to_x(upper_state, upper_rotated)
    call reactive_riemann_flux_x( &
      species, lower_rotated, upper_rotated, lower_temperature_guess, &
      upper_temperature_guess, solver, flux_rotated, ok)
    if (.not. ok) return
    call rotate_flux_x_to_y(flux_rotated, flux)
    ok = all(ieee_is_finite(flux))
  end subroutine reactive_riemann_flux_y

  pure real(dp) function lower_scale(center, slope, lower) result(theta)
    real(dp), intent(in) :: center, slope, lower
    real(dp) :: magnitude
    magnitude = abs(slope)
    if (magnitude <= tiny(1.0_dp) .or. center - magnitude > lower) then
      theta = 1.0_dp
    else
      theta = max(0.0_dp, min(1.0_dp, (center - lower) / magnitude))
    end if
  end function lower_scale

  pure real(dp) function upper_scale(center, slope, upper) result(theta)
    real(dp), intent(in) :: center, slope, upper
    real(dp) :: magnitude
    magnitude = abs(slope)
    if (magnitude <= tiny(1.0_dp) .or. center + magnitude < upper) then
      theta = 1.0_dp
    else
      theta = max(0.0_dp, min(1.0_dp, (upper - center) / magnitude))
    end if
  end function upper_scale

  pure real(dp) function primitive_slope_scale(center, slope, nspecies) &
      result(theta)
    real(dp), intent(in) :: center(:), slope(:)
    integer, intent(in) :: nspecies
    integer :: k, component

    theta = 1.0_dp
    if (size(center) /= size(slope)) then
      theta = 0.0_dp
      return
    end if
    theta = min(theta, lower_scale(center(1), slope(1), density_floor))
    theta = min(theta, lower_scale(center(5), slope(5), pressure_floor))
    do k = 1, nspecies
      component = reactive_mass_fraction_component(k)
      theta = min(theta, lower_scale(center(component), slope(component), 0.0_dp))
      theta = min(theta, upper_scale(center(component), slope(component), 1.0_dp))
    end do
  end function primitive_slope_scale

  pure subroutine sanitize_primitive(q, fallback, nspecies)
    real(dp), intent(inout) :: q(:)
    real(dp), intent(in) :: fallback(:)
    integer, intent(in) :: nspecies
    real(dp) :: total
    integer :: k, component

    if (size(q) /= size(fallback) .or. q(1) <= density_floor .or. &
        q(5) <= pressure_floor .or. any(.not. ieee_is_finite(q))) then
      q = fallback
      return
    end if
    total = 0.0_dp
    do k = 1, nspecies
      component = reactive_mass_fraction_component(k)
      if (q(component) < -1.0e-12_dp) then
        q = fallback
        return
      end if
      q(component) = max(0.0_dp, q(component))
      total = total + q(component)
    end do
    if (total <= tiny(1.0_dp)) then
      q = fallback
      return
    end if
    do k = 1, nspecies
      component = reactive_mass_fraction_component(k)
      q(component) = q(component) / total
    end do
  end subroutine sanitize_primitive

  subroutine primitive_face_to_state( &
      species, face_primitive, fallback_primitive, state, temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: face_primitive(:), fallback_primitive(:)
    real(dp), intent(out) :: state(:), temperature
    logical, intent(out) :: ok

    real(dp), allocatable :: work(:)
    real(dp) :: sound_speed

    allocate(work(size(face_primitive)))
    work = face_primitive
    call sanitize_primitive(work, fallback_primitive, size(species))
    call reactive_primitive_to_conserved( &
      species, work, state, temperature, sound_speed, ok)
    if (ok) return
    call reactive_primitive_to_conserved( &
      species, fallback_primitive, state, temperature, sound_speed, ok)
  end subroutine primitive_face_to_state

  subroutine reconstruct_reactive_characteristic_ppm_cell_2d( &
      species, stencil, sound_speed, dtdn, use_contact_steepening, &
      use_shock_flattening, minus_state, plus_state, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: stencil(:, -3:)
    real(dp), intent(in) :: sound_speed, dtdn
    logical, intent(in) :: use_contact_steepening, use_shock_flattening
    real(dp), intent(out) :: minus_state(:), plus_state(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: edge_minus(:), edge_plus(:)
    real(dp), allocatable :: integral_right(:, :), integral_left(:, :)
    real(dp) :: five_point(5), contact_stencil(-2:2)
    real(dp) :: density_stencil(-2:2), pressure_stencil(-2:2)
    real(dp) :: pressure_wide(-3:3), velocity_wide(-3:3)
    real(dp) :: profile_right(3), profile_left(3)
    real(dp) :: flattening, eta, gamma_effective
    logical :: local_ok
    integer :: nspecies, nprimitive, component, offset, k

    minus_state = 0.0_dp
    plus_state = 0.0_dp
    ok = .false.
    nspecies = size(species)
    nprimitive = reactive_nprim(nspecies)
    if (size(stencil, 1) /= nprimitive .or. size(stencil, 2) < 7 .or. &
        size(minus_state) /= nprimitive .or. size(plus_state) /= nprimitive .or. &
        sound_speed <= 0.0_dp .or. dtdn < 0.0_dp) return

    allocate(edge_minus(nprimitive), edge_plus(nprimitive))
    allocate(integral_right(nprimitive, 3), integral_left(nprimitive, 3))

    flattening = 1.0_dp
    if (use_shock_flattening) then
      do offset = -3, 3
        pressure_wide(offset) = stencil(5, offset)
        velocity_wide(offset) = stencil(2, offset)
      end do
      flattening = reactive_ppm_flattening_coefficient( &
        pressure_wide, velocity_wide)
    end if

    do component = 1, nprimitive
      five_point = stencil(component, -2:2)
      call reactive_ppm_reconstruct_five( &
        five_point, flattening, edge_minus(component), edge_plus(component))
    end do

    if (use_contact_steepening .and. flattening > 0.999_dp) then
      do offset = -2, 2
        density_stencil(offset) = stencil(1, offset)
        pressure_stencil(offset) = stencil(5, offset)
      end do
      gamma_effective = sound_speed**2 * stencil(1, 0) / stencil(5, 0)
      eta = min(contact_steepening_cap_2d, &
        reactive_ppm_contact_steepening_factor( &
          density_stencil, pressure_stencil, gamma_effective))
      if (eta > 0.0_dp) then
        call reactive_ppm_apply_contact_steepening( &
          density_stencil, eta, edge_minus(1), edge_plus(1))
        do k = 1, nspecies
          component = reactive_mass_fraction_component(k)
          do offset = -2, 2
            contact_stencil(offset) = stencil(component, offset)
          end do
          call reactive_ppm_apply_contact_steepening( &
            contact_stencil, eta, edge_minus(component), edge_plus(component))
        end do
      end if
    end if

    call sanitize_primitive(edge_minus, stencil(:, 0), nspecies)
    call sanitize_primitive(edge_plus, stencil(:, 0), nspecies)
    do component = 1, nprimitive
      call reactive_ppm_integrate_profile( &
        edge_minus(component), edge_plus(component), stencil(component, 0), &
        stencil(2, 0), sound_speed, dtdn, profile_right, profile_left, local_ok)
      if (.not. local_ok) return
      integral_right(component, :) = profile_right
      integral_left(component, :) = profile_left
    end do
    call build_characteristic_ppm_states( &
      species, stencil(:, 0), sound_speed, integral_right, integral_left, &
      minus_state, plus_state, local_ok)
    if (.not. local_ok) return
    call sanitize_primitive(minus_state, stencil(:, 0), nspecies)
    call sanitize_primitive(plus_state, stencil(:, 0), nspecies)
    ok = all(ieee_is_finite(minus_state)) .and. &
      all(ieee_is_finite(plus_state))
  end subroutine reconstruct_reactive_characteristic_ppm_cell_2d

  subroutine recover_state_temperature( &
      species, state, temperature_guess, temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:), temperature_guess
    real(dp), intent(out) :: temperature
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: sound_speed

    allocate(primitive(reactive_nprim(size(species))))
    call reactive_conserved_to_primitive( &
      species, state, temperature_guess, primitive, temperature, sound_speed, ok)
  end subroutine recover_state_temperature

  subroutine compute_reactive_cfl_timestep_2d( &
      species, state, temperature, nx, ny, dx, dy, cfl, dt, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: dx, dy, cfl
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: local_temperature, sound_speed, rate, maximum_rate
    logical :: local_ok
    integer :: i, j

    dt = 0.0_dp
    ok = .false.
    if (size(state, 2) /= nx .or. size(state, 3) /= ny .or. &
        size(temperature, 1) /= nx .or. size(temperature, 2) /= ny .or. &
        nx < 2 .or. ny < 2 .or. dx <= 0.0_dp .or. dy <= 0.0_dp .or. &
        cfl <= 0.0_dp .or. cfl > 1.0_dp) return
    allocate(primitive(reactive_nprim(size(species))))
    maximum_rate = 0.0_dp
    do j = 1, ny
      do i = 1, nx
        call reactive_conserved_to_primitive( &
          species, state(:, i, j), temperature(i, j), primitive, &
          local_temperature, sound_speed, local_ok)
        if (.not. local_ok) return
        rate = (abs(primitive(2)) + sound_speed) / dx + &
          (abs(primitive(3)) + sound_speed) / dy
        maximum_rate = max(maximum_rate, rate)
      end do
    end do
    if (maximum_rate <= 0.0_dp) return
    dt = cfl / maximum_rate
    ok = ieee_is_finite(dt) .and. dt > 0.0_dp
  end subroutine compute_reactive_cfl_timestep_2d

  subroutine apply_reactive_transverse_correction_2d( &
      species, base_state, base_temperature, flux_high, flux_low, scale, &
      corrected_state, corrected_temperature, theta, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: base_state(:), base_temperature
    real(dp), intent(in) :: flux_high(:), flux_low(:), scale
    real(dp), intent(out) :: corrected_state(:), corrected_temperature, theta
    logical, intent(out) :: ok

    real(dp), allocatable :: correction(:), trial(:)
    real(dp) :: lower_theta, upper_theta, midpoint, trial_temperature
    logical :: local_ok
    integer :: iteration

    corrected_state = base_state
    corrected_temperature = base_temperature
    theta = 0.0_dp
    ok = .false.
    if (size(base_state) /= size(flux_high) .or. &
        size(base_state) /= size(flux_low) .or. &
        size(base_state) /= size(corrected_state) .or. scale < 0.0_dp) return
    call recover_state_temperature( &
      species, base_state, base_temperature, trial_temperature, local_ok)
    if (.not. local_ok) return
    allocate(correction(size(base_state)), trial(size(base_state)))
    correction = -scale * (flux_high - flux_low)
    trial = base_state + correction
    call recover_state_temperature( &
      species, trial, base_temperature, trial_temperature, local_ok)
    if (local_ok) then
      corrected_state = trial
      corrected_temperature = trial_temperature
      theta = 1.0_dp
      ok = .true.
      return
    end if

    lower_theta = 0.0_dp
    upper_theta = 1.0_dp
    do iteration = 1, 60
      midpoint = 0.5_dp * (lower_theta + upper_theta)
      trial = base_state + midpoint * correction
      call recover_state_temperature( &
        species, trial, base_temperature, trial_temperature, local_ok)
      if (local_ok) then
        lower_theta = midpoint
      else
        upper_theta = midpoint
      end if
    end do
    theta = lower_theta
    corrected_state = base_state + theta * correction
    call recover_state_temperature( &
      species, corrected_state, base_temperature, corrected_temperature, ok)
  end subroutine apply_reactive_transverse_correction_2d

  subroutine advance_reactive_hydro_periodic_2d( &
      species, state, temperature, nx, ny, dx, dy, dt, reconstruction, &
      limiter, riemann_solver, use_transverse_correction, ok, &
      minimum_transverse_theta, ppm_contact_steepening, ppm_shock_flattening)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:, :, :), temperature(:, :)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: dx, dy, dt
    character(len=*), intent(in) :: reconstruction, limiter, riemann_solver
    logical, intent(in) :: use_transverse_correction
    logical, intent(out) :: ok
    real(dp), intent(out), optional :: minimum_transverse_theta
    logical, intent(in), optional :: ppm_contact_steepening
    logical, intent(in), optional :: ppm_shock_flattening

    real(dp), allocatable :: primitive(:, :, :), sound_speed(:, :)
    real(dp), allocatable :: slope_x(:, :, :), slope_y(:, :, :)
    real(dp), allocatable :: x_minus(:, :, :), x_plus(:, :, :)
    real(dp), allocatable :: y_minus(:, :, :), y_plus(:, :, :)
    real(dp), allocatable :: x_left_base(:, :, :), x_right_base(:, :, :)
    real(dp), allocatable :: y_lower_base(:, :, :), y_upper_base(:, :, :)
    real(dp), allocatable :: provisional_x_flux(:, :, :)
    real(dp), allocatable :: provisional_y_flux(:, :, :)
    real(dp), allocatable :: final_x_flux(:, :, :), final_y_flux(:, :, :)
    real(dp), allocatable :: x_left_temperature(:, :), x_right_temperature(:, :)
    real(dp), allocatable :: y_lower_temperature(:, :), y_upper_temperature(:, :)
    real(dp), allocatable :: new_state(:, :, :), new_temperature(:, :)
    real(dp), allocatable :: dl(:), dr(:), rotated_center(:), rotated_dl(:)
    real(dp), allocatable :: rotated_dr(:), rotated_slope(:)
    real(dp), allocatable :: rotated_minus(:), rotated_plus(:)
    real(dp), allocatable :: directional_stencil(:, :)
    real(dp), allocatable :: x_left(:), x_right(:), y_lower(:), y_upper(:)
    real(dp) :: local_temperature, theta, theta_x_left, theta_x_right
    real(dp) :: theta_y_lower, theta_y_upper, local_minimum_theta
    real(dp) :: xlt, xrt, ylt, yut
    logical :: local_ok, face_ok, correction_ok
    logical :: use_contact_steepening, use_shock_flattening
    integer :: nvar, nprimitive, i, j, im, ip, jm, jp, offset

    ok = .false.
    use_contact_steepening = .false.
    use_shock_flattening = .false.
    if (present(ppm_contact_steepening)) &
      use_contact_steepening = ppm_contact_steepening
    if (present(ppm_shock_flattening)) &
      use_shock_flattening = ppm_shock_flattening
    if (present(minimum_transverse_theta)) minimum_transverse_theta = 0.0_dp
    nvar = reactive_nvar(size(species))
    nprimitive = reactive_nprim(size(species))
    if (size(state, 1) /= nvar .or. size(state, 2) /= nx .or. &
        size(state, 3) /= ny .or. size(temperature, 1) /= nx .or. &
        size(temperature, 2) /= ny .or. nx < 4 .or. ny < 4 .or. &
        dx <= 0.0_dp .or. dy <= 0.0_dp .or. dt <= 0.0_dp) return
    if (trim(reconstruction) /= "pcm" .and. &
        trim(reconstruction) /= "characteristic_plm" .and. &
        trim(reconstruction) /= "characteristic_ppm") return
    if ((use_contact_steepening .or. use_shock_flattening) .and. &
        trim(reconstruction) /= "characteristic_ppm") return

    allocate(primitive(nprimitive, nx, ny), sound_speed(nx, ny))
    allocate(slope_x(nprimitive, nx, ny), slope_y(nprimitive, nx, ny))
    allocate(x_minus(nprimitive, nx, ny), x_plus(nprimitive, nx, ny))
    allocate(y_minus(nprimitive, nx, ny), y_plus(nprimitive, nx, ny))
    allocate(x_left_base(nvar, nx, ny), x_right_base(nvar, nx, ny))
    allocate(y_lower_base(nvar, nx, ny), y_upper_base(nvar, nx, ny))
    allocate(provisional_x_flux(nvar, nx, ny), provisional_y_flux(nvar, nx, ny))
    allocate(final_x_flux(nvar, nx, ny), final_y_flux(nvar, nx, ny))
    allocate(x_left_temperature(nx, ny), x_right_temperature(nx, ny))
    allocate(y_lower_temperature(nx, ny), y_upper_temperature(nx, ny))
    allocate(new_state(nvar, nx, ny), new_temperature(nx, ny))
    allocate(dl(nprimitive), dr(nprimitive), rotated_center(nprimitive))
    allocate(rotated_dl(nprimitive), rotated_dr(nprimitive))
    allocate(rotated_slope(nprimitive), rotated_minus(nprimitive))
    allocate(rotated_plus(nprimitive))
    allocate(directional_stencil(nprimitive, -3:3))
    allocate(x_left(nvar), x_right(nvar), y_lower(nvar), y_upper(nvar))

    slope_x = 0.0_dp
    slope_y = 0.0_dp
    do j = 1, ny
      do i = 1, nx
        call reactive_conserved_to_primitive( &
          species, state(:, i, j), temperature(i, j), primitive(:, i, j), &
          local_temperature, sound_speed(i, j), local_ok)
        if (.not. local_ok) return
      end do
    end do

    if (trim(reconstruction) == "characteristic_plm") then
      do j = 1, ny
        jm = periodic_index(j - 1, ny)
        jp = periodic_index(j + 1, ny)
        do i = 1, nx
          im = periodic_index(i - 1, nx)
          ip = periodic_index(i + 1, nx)
          dl = primitive(:, i, j) - primitive(:, im, j)
          dr = primitive(:, ip, j) - primitive(:, i, j)
          call characteristic_limited_slope( &
            primitive(:, i, j), dl, dr, sound_speed(i, j), limiter, &
            slope_x(:, i, j), local_ok)
          if (.not. local_ok) return
          theta = primitive_slope_scale( &
            primitive(:, i, j), slope_x(:, i, j), size(species))
          slope_x(:, i, j) = theta * slope_x(:, i, j)

          call rotate_primitive_y_to_x(primitive(:, i, j), rotated_center)
          call rotate_primitive_y_to_x( &
            primitive(:, i, j) - primitive(:, i, jm), rotated_dl)
          call rotate_primitive_y_to_x( &
            primitive(:, i, jp) - primitive(:, i, j), rotated_dr)
          call characteristic_limited_slope( &
            rotated_center, rotated_dl, rotated_dr, sound_speed(i, j), &
            limiter, rotated_slope, local_ok)
          if (.not. local_ok) return
          theta = primitive_slope_scale( &
            rotated_center, rotated_slope, size(species))
          rotated_slope = theta * rotated_slope
          call rotate_primitive_x_to_y(rotated_slope, slope_y(:, i, j))
        end do
      end do
    end if

    do j = 1, ny
      do i = 1, nx
        select case (trim(reconstruction))
        case ("pcm")
          x_minus(:, i, j) = primitive(:, i, j)
          x_plus(:, i, j) = primitive(:, i, j)
          y_minus(:, i, j) = primitive(:, i, j)
          y_plus(:, i, j) = primitive(:, i, j)

        case ("characteristic_plm")
          call trace_reactive_characteristics( &
            primitive(:, i, j), slope_x(:, i, j), sound_speed(i, j), &
            dt / dx, x_minus(:, i, j), x_plus(:, i, j), local_ok)
          if (.not. local_ok) return
          call sanitize_primitive( &
            x_minus(:, i, j), primitive(:, i, j), size(species))
          call sanitize_primitive( &
            x_plus(:, i, j), primitive(:, i, j), size(species))

          call rotate_primitive_y_to_x(primitive(:, i, j), rotated_center)
          call rotate_primitive_y_to_x(slope_y(:, i, j), rotated_slope)
          call trace_reactive_characteristics( &
            rotated_center, rotated_slope, sound_speed(i, j), dt / dy, &
            rotated_minus, rotated_plus, local_ok)
          if (.not. local_ok) return
          call rotate_primitive_x_to_y(rotated_minus, y_minus(:, i, j))
          call rotate_primitive_x_to_y(rotated_plus, y_plus(:, i, j))
          call sanitize_primitive( &
            y_minus(:, i, j), primitive(:, i, j), size(species))
          call sanitize_primitive( &
            y_plus(:, i, j), primitive(:, i, j), size(species))

        case ("characteristic_ppm")
          do offset = -3, 3
            directional_stencil(:, offset) = &
              primitive(:, periodic_index(i + offset, nx), j)
          end do
          call reconstruct_reactive_characteristic_ppm_cell_2d( &
            species, directional_stencil, sound_speed(i, j), dt / dx, &
            use_contact_steepening, use_shock_flattening, &
            x_minus(:, i, j), x_plus(:, i, j), local_ok)
          if (.not. local_ok) return

          do offset = -3, 3
            call rotate_primitive_y_to_x( &
              primitive(:, i, periodic_index(j + offset, ny)), &
              directional_stencil(:, offset))
          end do
          call reconstruct_reactive_characteristic_ppm_cell_2d( &
            species, directional_stencil, sound_speed(i, j), dt / dy, &
            use_contact_steepening, use_shock_flattening, &
            rotated_minus, rotated_plus, local_ok)
          if (.not. local_ok) return
          call rotate_primitive_x_to_y(rotated_minus, y_minus(:, i, j))
          call rotate_primitive_x_to_y(rotated_plus, y_plus(:, i, j))
          call sanitize_primitive( &
            y_minus(:, i, j), primitive(:, i, j), size(species))
          call sanitize_primitive( &
            y_plus(:, i, j), primitive(:, i, j), size(species))

        case default
          return
        end select
      end do
    end do

    do j = 1, ny
      jp = periodic_index(j + 1, ny)
      do i = 1, nx
        ip = periodic_index(i + 1, nx)
        call primitive_face_to_state( &
          species, x_plus(:, i, j), primitive(:, i, j), &
          x_left_base(:, i, j), x_left_temperature(i, j), local_ok)
        if (.not. local_ok) return
        call primitive_face_to_state( &
          species, x_minus(:, ip, j), primitive(:, ip, j), &
          x_right_base(:, i, j), x_right_temperature(i, j), local_ok)
        if (.not. local_ok) return
        call primitive_face_to_state( &
          species, y_plus(:, i, j), primitive(:, i, j), &
          y_lower_base(:, i, j), y_lower_temperature(i, j), local_ok)
        if (.not. local_ok) return
        call primitive_face_to_state( &
          species, y_minus(:, i, jp), primitive(:, i, jp), &
          y_upper_base(:, i, j), y_upper_temperature(i, j), local_ok)
        if (.not. local_ok) return

        call reactive_riemann_flux_x( &
          species, x_left_base(:, i, j), x_right_base(:, i, j), &
          x_left_temperature(i, j), x_right_temperature(i, j), &
          riemann_solver, provisional_x_flux(:, i, j), face_ok)
        if (.not. face_ok) return
        call reactive_riemann_flux_y( &
          species, y_lower_base(:, i, j), y_upper_base(:, i, j), &
          y_lower_temperature(i, j), y_upper_temperature(i, j), &
          riemann_solver, provisional_y_flux(:, i, j), face_ok)
        if (.not. face_ok) return
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
          call apply_reactive_transverse_correction_2d( &
            species, x_left_base(:, i, j), x_left_temperature(i, j), &
            provisional_y_flux(:, i, j), provisional_y_flux(:, i, jm), &
            0.5_dp * dt / dy, x_left, xlt, theta_x_left, correction_ok)
          if (.not. correction_ok) return
          call apply_reactive_transverse_correction_2d( &
            species, x_right_base(:, i, j), x_right_temperature(i, j), &
            provisional_y_flux(:, ip, j), provisional_y_flux(:, ip, jm), &
            0.5_dp * dt / dy, x_right, xrt, theta_x_right, correction_ok)
          if (.not. correction_ok) return
          call apply_reactive_transverse_correction_2d( &
            species, y_lower_base(:, i, j), y_lower_temperature(i, j), &
            provisional_x_flux(:, i, j), provisional_x_flux(:, im, j), &
            0.5_dp * dt / dx, y_lower, ylt, theta_y_lower, correction_ok)
          if (.not. correction_ok) return
          call apply_reactive_transverse_correction_2d( &
            species, y_upper_base(:, i, j), y_upper_temperature(i, j), &
            provisional_x_flux(:, i, jp), provisional_x_flux(:, im, jp), &
            0.5_dp * dt / dx, y_upper, yut, theta_y_upper, correction_ok)
          if (.not. correction_ok) return
        else
          x_left = x_left_base(:, i, j)
          x_right = x_right_base(:, i, j)
          y_lower = y_lower_base(:, i, j)
          y_upper = y_upper_base(:, i, j)
          xlt = x_left_temperature(i, j)
          xrt = x_right_temperature(i, j)
          ylt = y_lower_temperature(i, j)
          yut = y_upper_temperature(i, j)
          theta_x_left = 1.0_dp
          theta_x_right = 1.0_dp
          theta_y_lower = 1.0_dp
          theta_y_upper = 1.0_dp
        end if
        local_minimum_theta = min(local_minimum_theta, theta_x_left, &
          theta_x_right, theta_y_lower, theta_y_upper)
        call reactive_riemann_flux_x( &
          species, x_left, x_right, xlt, xrt, riemann_solver, &
          final_x_flux(:, i, j), face_ok)
        if (.not. face_ok) return
        call reactive_riemann_flux_y( &
          species, y_lower, y_upper, ylt, yut, riemann_solver, &
          final_y_flux(:, i, j), face_ok)
        if (.not. face_ok) return
      end do
    end do

    do j = 1, ny
      jm = periodic_index(j - 1, ny)
      do i = 1, nx
        im = periodic_index(i - 1, nx)
        new_state(:, i, j) = state(:, i, j) - &
          dt / dx * (final_x_flux(:, i, j) - final_x_flux(:, im, j)) - &
          dt / dy * (final_y_flux(:, i, j) - final_y_flux(:, i, jm))
        call recover_state_temperature( &
          species, new_state(:, i, j), temperature(i, j), &
          new_temperature(i, j), local_ok)
        if (.not. local_ok) return
      end do
    end do
    state = new_state
    temperature = new_temperature
    if (present(minimum_transverse_theta)) then
      minimum_transverse_theta = local_minimum_theta
    end if
    ok = .true.
  end subroutine advance_reactive_hydro_periodic_2d

  pure subroutine reactive_wall_flux_x(pressure, flux)
    real(dp), intent(in) :: pressure
    real(dp), intent(out) :: flux(:)
    flux = 0.0_dp
    if (size(flux) >= imx) flux(imx) = pressure
  end subroutine reactive_wall_flux_x

  pure subroutine reactive_wall_flux_y(pressure, flux)
    real(dp), intent(in) :: pressure
    real(dp), intent(out) :: flux(:)
    flux = 0.0_dp
    if (size(flux) >= imy) flux(imy) = pressure
  end subroutine reactive_wall_flux_y

  subroutine advance_reactive_hydro_physical_2d( &
      species, state, temperature, nx, ny, dx, dy, dt, reconstruction, &
      limiter, riemann_solver, use_transverse_correction, boundaries, ok, &
      minimum_transverse_theta, ppm_contact_steepening, ppm_shock_flattening)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:, :, :), temperature(:, :)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: dx, dy, dt
    character(len=*), intent(in) :: reconstruction, limiter, riemann_solver
    logical, intent(in) :: use_transverse_correction
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    logical, intent(out) :: ok
    real(dp), intent(out) :: minimum_transverse_theta
    logical, intent(in) :: ppm_contact_steepening, ppm_shock_flattening

    real(dp), allocatable :: primitive(:, :, :), sound_speed(:, :)
    real(dp), allocatable :: slope_x(:, :, :), slope_y(:, :, :)
    real(dp), allocatable :: x_minus(:, :, :), x_plus(:, :, :)
    real(dp), allocatable :: y_minus(:, :, :), y_plus(:, :, :)
    real(dp), allocatable :: x_left_base(:, :, :), x_right_base(:, :, :)
    real(dp), allocatable :: y_lower_base(:, :, :), y_upper_base(:, :, :)
    real(dp), allocatable :: provisional_x_flux(:, :, :)
    real(dp), allocatable :: provisional_y_flux(:, :, :)
    real(dp), allocatable :: final_x_flux(:, :, :), final_y_flux(:, :, :)
    real(dp), allocatable :: x_left_temperature(:, :), x_right_temperature(:, :)
    real(dp), allocatable :: y_lower_temperature(:, :), y_upper_temperature(:, :)
    real(dp), allocatable :: new_state(:, :, :), new_temperature(:, :)
    real(dp), allocatable :: dl(:), dr(:), qminus(:), qplus(:)
    real(dp), allocatable :: rotated_center(:), rotated_dl(:), rotated_dr(:)
    real(dp), allocatable :: rotated_slope(:), rotated_minus(:), rotated_plus(:)
    real(dp), allocatable :: directional_stencil(:, :), sample_q(:)
    real(dp), allocatable :: x_left(:), x_right(:), y_lower(:), y_upper(:)
    real(dp) :: sample_t, local_temperature, theta
    real(dp) :: theta_x_left, theta_x_right, theta_y_lower, theta_y_upper
    real(dp) :: xlt, xrt, ylt, yut, wall_pressure
    logical :: local_ok, face_ok, correction_ok, physical_face
    logical :: periodic_x, periodic_y
    integer :: nvar, nprimitive, i, j, face_i, face_j, li, ri, lj, uj, offset

    ok = .false.
    minimum_transverse_theta = 1.0_dp
    nvar = reactive_nvar(size(species))
    nprimitive = reactive_nprim(size(species))
    if (size(state, 1) /= nvar .or. size(state, 2) /= nx .or. &
        size(state, 3) /= ny .or. size(temperature, 1) /= nx .or. &
        size(temperature, 2) /= ny .or. nx < 4 .or. ny < 4 .or. &
        dx <= 0.0_dp .or. dy <= 0.0_dp .or. dt <= 0.0_dp) return
    if (trim(reconstruction) /= 'pcm' .and. &
        trim(reconstruction) /= 'characteristic_plm' .and. &
        trim(reconstruction) /= 'characteristic_ppm') return
    if ((ppm_contact_steepening .or. ppm_shock_flattening) .and. &
        trim(reconstruction) /= 'characteristic_ppm') return
    periodic_x = reactive_boundary_is_periodic(boundaries%face(1))
    periodic_y = reactive_boundary_is_periodic(boundaries%face(3))

    allocate(primitive(nprimitive, nx, ny), sound_speed(nx, ny))
    allocate(slope_x(nprimitive, nx, ny), slope_y(nprimitive, nx, ny))
    allocate(x_minus(nprimitive, nx, ny), x_plus(nprimitive, nx, ny))
    allocate(y_minus(nprimitive, nx, ny), y_plus(nprimitive, nx, ny))
    allocate(x_left_base(nvar, 0:nx, ny), x_right_base(nvar, 0:nx, ny))
    allocate(y_lower_base(nvar, nx, 0:ny), y_upper_base(nvar, nx, 0:ny))
    allocate(provisional_x_flux(nvar, 0:nx, ny))
    allocate(provisional_y_flux(nvar, nx, 0:ny))
    allocate(final_x_flux(nvar, 0:nx, ny), final_y_flux(nvar, nx, 0:ny))
    allocate(x_left_temperature(0:nx, ny), x_right_temperature(0:nx, ny))
    allocate(y_lower_temperature(nx, 0:ny), y_upper_temperature(nx, 0:ny))
    allocate(new_state(nvar, nx, ny), new_temperature(nx, ny))
    allocate(dl(nprimitive), dr(nprimitive), qminus(nprimitive), qplus(nprimitive))
    allocate(rotated_center(nprimitive), rotated_dl(nprimitive))
    allocate(rotated_dr(nprimitive), rotated_slope(nprimitive))
    allocate(rotated_minus(nprimitive), rotated_plus(nprimitive))
    allocate(directional_stencil(nprimitive, -3:3), sample_q(nprimitive))
    allocate(x_left(nvar), x_right(nvar), y_lower(nvar), y_upper(nvar))

    slope_x = 0.0_dp
    slope_y = 0.0_dp
    do j = 1, ny
      do i = 1, nx
        call reactive_conserved_to_primitive( &
          species, state(:, i, j), temperature(i, j), primitive(:, i, j), &
          local_temperature, sound_speed(i, j), local_ok)
        if (.not. local_ok) return
      end do
    end do

    if (trim(reconstruction) == 'characteristic_plm') then
      do j = 1, ny
        do i = 1, nx
          call sample_reactive_primitive_2d( &
            primitive, temperature, nx, ny, i - 1, j, boundaries, qminus, &
            sample_t, local_ok)
          if (.not. local_ok) return
          call sample_reactive_primitive_2d( &
            primitive, temperature, nx, ny, i + 1, j, boundaries, qplus, &
            sample_t, local_ok)
          if (.not. local_ok) return
          dl = primitive(:, i, j) - qminus
          dr = qplus - primitive(:, i, j)
          call characteristic_limited_slope( &
            primitive(:, i, j), dl, dr, sound_speed(i, j), limiter, &
            slope_x(:, i, j), local_ok)
          if (.not. local_ok) return
          theta = primitive_slope_scale( &
            primitive(:, i, j), slope_x(:, i, j), size(species))
          slope_x(:, i, j) = theta * slope_x(:, i, j)

          call sample_reactive_primitive_2d( &
            primitive, temperature, nx, ny, i, j - 1, boundaries, qminus, &
            sample_t, local_ok)
          if (.not. local_ok) return
          call sample_reactive_primitive_2d( &
            primitive, temperature, nx, ny, i, j + 1, boundaries, qplus, &
            sample_t, local_ok)
          if (.not. local_ok) return
          call rotate_primitive_y_to_x(primitive(:, i, j), rotated_center)
          call rotate_primitive_y_to_x(primitive(:, i, j) - qminus, rotated_dl)
          call rotate_primitive_y_to_x(qplus - primitive(:, i, j), rotated_dr)
          call characteristic_limited_slope( &
            rotated_center, rotated_dl, rotated_dr, sound_speed(i, j), &
            limiter, rotated_slope, local_ok)
          if (.not. local_ok) return
          theta = primitive_slope_scale( &
            rotated_center, rotated_slope, size(species))
          rotated_slope = theta * rotated_slope
          call rotate_primitive_x_to_y(rotated_slope, slope_y(:, i, j))
        end do
      end do
    end if

    do j = 1, ny
      do i = 1, nx
        select case (trim(reconstruction))
        case ('pcm')
          x_minus(:, i, j) = primitive(:, i, j)
          x_plus(:, i, j) = primitive(:, i, j)
          y_minus(:, i, j) = primitive(:, i, j)
          y_plus(:, i, j) = primitive(:, i, j)
        case ('characteristic_plm')
          call trace_reactive_characteristics( &
            primitive(:, i, j), slope_x(:, i, j), sound_speed(i, j), &
            dt / dx, x_minus(:, i, j), x_plus(:, i, j), local_ok)
          if (.not. local_ok) return
          call sanitize_primitive(x_minus(:, i, j), primitive(:, i, j), size(species))
          call sanitize_primitive(x_plus(:, i, j), primitive(:, i, j), size(species))
          call rotate_primitive_y_to_x(primitive(:, i, j), rotated_center)
          call rotate_primitive_y_to_x(slope_y(:, i, j), rotated_slope)
          call trace_reactive_characteristics( &
            rotated_center, rotated_slope, sound_speed(i, j), dt / dy, &
            rotated_minus, rotated_plus, local_ok)
          if (.not. local_ok) return
          call rotate_primitive_x_to_y(rotated_minus, y_minus(:, i, j))
          call rotate_primitive_x_to_y(rotated_plus, y_plus(:, i, j))
          call sanitize_primitive(y_minus(:, i, j), primitive(:, i, j), size(species))
          call sanitize_primitive(y_plus(:, i, j), primitive(:, i, j), size(species))
        case ('characteristic_ppm')
          do offset = -3, 3
            call sample_reactive_primitive_2d( &
              primitive, temperature, nx, ny, i + offset, j, boundaries, &
              directional_stencil(:, offset), sample_t, local_ok)
            if (.not. local_ok) return
          end do
          call reconstruct_reactive_characteristic_ppm_cell_2d( &
            species, directional_stencil, sound_speed(i, j), dt / dx, &
            ppm_contact_steepening, ppm_shock_flattening, &
            x_minus(:, i, j), x_plus(:, i, j), local_ok)
          if (.not. local_ok) return
          do offset = -3, 3
            call sample_reactive_primitive_2d( &
              primitive, temperature, nx, ny, i, j + offset, boundaries, &
              sample_q, sample_t, local_ok)
            if (.not. local_ok) return
            call rotate_primitive_y_to_x(sample_q, directional_stencil(:, offset))
          end do
          call reconstruct_reactive_characteristic_ppm_cell_2d( &
            species, directional_stencil, sound_speed(i, j), dt / dy, &
            ppm_contact_steepening, ppm_shock_flattening, &
            rotated_minus, rotated_plus, local_ok)
          if (.not. local_ok) return
          call rotate_primitive_x_to_y(rotated_minus, y_minus(:, i, j))
          call rotate_primitive_x_to_y(rotated_plus, y_plus(:, i, j))
          call sanitize_primitive(y_minus(:, i, j), primitive(:, i, j), size(species))
          call sanitize_primitive(y_plus(:, i, j), primitive(:, i, j), size(species))
        end select
      end do
    end do

    ! Provisional x-normal fluxes on explicit faces.
    do j = 1, ny
      do face_i = 0, nx
        physical_face = (face_i == 0 .and. .not. periodic_x) .or. &
          (face_i == nx .and. .not. periodic_x)
        if (physical_face .and. &
            reactive_boundary_is_wall(boundaries%face(merge(1, 2, face_i == 0)))) then
          i = merge(1, nx, face_i == 0)
          wall_pressure = primitive(5, i, j)
          call reactive_wall_flux_x(wall_pressure, provisional_x_flux(:, face_i, j))
          x_left_base(:, face_i, j) = state(:, i, j)
          x_right_base(:, face_i, j) = state(:, i, j)
          x_left_temperature(face_i, j) = temperature(i, j)
          x_right_temperature(face_i, j) = temperature(i, j)
        else
          if (physical_face) then
            if (face_i == 0) then
              call sample_reactive_primitive_2d( &
                primitive, temperature, nx, ny, 0, j, boundaries, qminus, &
                sample_t, local_ok)
              if (.not. local_ok) return
              call primitive_face_to_state( &
                species, qminus, qminus, x_left_base(:, face_i, j), &
                x_left_temperature(face_i, j), local_ok)
              if (.not. local_ok) return
              call primitive_face_to_state( &
                species, x_minus(:, 1, j), primitive(:, 1, j), &
                x_right_base(:, face_i, j), x_right_temperature(face_i, j), local_ok)
              if (.not. local_ok) return
            else
              call primitive_face_to_state( &
                species, x_plus(:, nx, j), primitive(:, nx, j), &
                x_left_base(:, face_i, j), x_left_temperature(face_i, j), local_ok)
              if (.not. local_ok) return
              call sample_reactive_primitive_2d( &
                primitive, temperature, nx, ny, nx + 1, j, boundaries, qplus, &
                sample_t, local_ok)
              if (.not. local_ok) return
              call primitive_face_to_state( &
                species, qplus, qplus, x_right_base(:, face_i, j), &
                x_right_temperature(face_i, j), local_ok)
              if (.not. local_ok) return
            end if
          else
            li = face_i
            ri = face_i + 1
            if (face_i == 0) li = nx
            if (face_i == nx) ri = 1
            call primitive_face_to_state( &
              species, x_plus(:, li, j), primitive(:, li, j), &
              x_left_base(:, face_i, j), x_left_temperature(face_i, j), local_ok)
            if (.not. local_ok) return
            call primitive_face_to_state( &
              species, x_minus(:, ri, j), primitive(:, ri, j), &
              x_right_base(:, face_i, j), x_right_temperature(face_i, j), local_ok)
            if (.not. local_ok) return
          end if
          call reactive_riemann_flux_x( &
            species, x_left_base(:, face_i, j), x_right_base(:, face_i, j), &
            x_left_temperature(face_i, j), x_right_temperature(face_i, j), &
            riemann_solver, provisional_x_flux(:, face_i, j), face_ok)
          if (.not. face_ok) return
        end if
      end do
    end do

    ! Provisional y-normal fluxes on explicit faces.
    do face_j = 0, ny
      do i = 1, nx
        physical_face = (face_j == 0 .and. .not. periodic_y) .or. &
          (face_j == ny .and. .not. periodic_y)
        if (physical_face .and. &
            reactive_boundary_is_wall(boundaries%face(merge(3, 4, face_j == 0)))) then
          j = merge(1, ny, face_j == 0)
          wall_pressure = primitive(5, i, j)
          call reactive_wall_flux_y(wall_pressure, provisional_y_flux(:, i, face_j))
          y_lower_base(:, i, face_j) = state(:, i, j)
          y_upper_base(:, i, face_j) = state(:, i, j)
          y_lower_temperature(i, face_j) = temperature(i, j)
          y_upper_temperature(i, face_j) = temperature(i, j)
        else
          if (physical_face) then
            if (face_j == 0) then
              call sample_reactive_primitive_2d( &
                primitive, temperature, nx, ny, i, 0, boundaries, qminus, &
                sample_t, local_ok)
              if (.not. local_ok) return
              call primitive_face_to_state( &
                species, qminus, qminus, y_lower_base(:, i, face_j), &
                y_lower_temperature(i, face_j), local_ok)
              if (.not. local_ok) return
              call primitive_face_to_state( &
                species, y_minus(:, i, 1), primitive(:, i, 1), &
                y_upper_base(:, i, face_j), y_upper_temperature(i, face_j), local_ok)
              if (.not. local_ok) return
            else
              call primitive_face_to_state( &
                species, y_plus(:, i, ny), primitive(:, i, ny), &
                y_lower_base(:, i, face_j), y_lower_temperature(i, face_j), local_ok)
              if (.not. local_ok) return
              call sample_reactive_primitive_2d( &
                primitive, temperature, nx, ny, i, ny + 1, boundaries, qplus, &
                sample_t, local_ok)
              if (.not. local_ok) return
              call primitive_face_to_state( &
                species, qplus, qplus, y_upper_base(:, i, face_j), &
                y_upper_temperature(i, face_j), local_ok)
              if (.not. local_ok) return
            end if
          else
            lj = face_j
            uj = face_j + 1
            if (face_j == 0) lj = ny
            if (face_j == ny) uj = 1
            call primitive_face_to_state( &
              species, y_plus(:, i, lj), primitive(:, i, lj), &
              y_lower_base(:, i, face_j), y_lower_temperature(i, face_j), local_ok)
            if (.not. local_ok) return
            call primitive_face_to_state( &
              species, y_minus(:, i, uj), primitive(:, i, uj), &
              y_upper_base(:, i, face_j), y_upper_temperature(i, face_j), local_ok)
            if (.not. local_ok) return
          end if
          call reactive_riemann_flux_y( &
            species, y_lower_base(:, i, face_j), y_upper_base(:, i, face_j), &
            y_lower_temperature(i, face_j), y_upper_temperature(i, face_j), &
            riemann_solver, provisional_y_flux(:, i, face_j), face_ok)
          if (.not. face_ok) return
        end if
      end do
    end do

    ! Correct interior/periodic faces. Exact wall fluxes remain untouched.
    do j = 1, ny
      do face_i = 0, nx
        physical_face = (face_i == 0 .and. .not. periodic_x) .or. &
          (face_i == nx .and. .not. periodic_x)
        if (physical_face .and. &
            reactive_boundary_is_wall(boundaries%face(merge(1, 2, face_i == 0)))) then
          final_x_flux(:, face_i, j) = provisional_x_flux(:, face_i, j)
          cycle
        end if
        li = face_i
        ri = face_i + 1
        if (face_i == 0) li = nx
        if (face_i == nx) ri = 1
        if (use_transverse_correction .and. .not. physical_face) then
          call apply_reactive_transverse_correction_2d( &
            species, x_left_base(:, face_i, j), x_left_temperature(face_i, j), &
            provisional_y_flux(:, li, j), provisional_y_flux(:, li, j - 1), &
            0.5_dp * dt / dy, x_left, xlt, theta_x_left, correction_ok)
          if (.not. correction_ok) return
          call apply_reactive_transverse_correction_2d( &
            species, x_right_base(:, face_i, j), x_right_temperature(face_i, j), &
            provisional_y_flux(:, ri, j), provisional_y_flux(:, ri, j - 1), &
            0.5_dp * dt / dy, x_right, xrt, theta_x_right, correction_ok)
          if (.not. correction_ok) return
        else
          x_left = x_left_base(:, face_i, j)
          x_right = x_right_base(:, face_i, j)
          xlt = x_left_temperature(face_i, j)
          xrt = x_right_temperature(face_i, j)
          theta_x_left = 1.0_dp
          theta_x_right = 1.0_dp
        end if
        minimum_transverse_theta = min( &
          minimum_transverse_theta, theta_x_left, theta_x_right)
        call reactive_riemann_flux_x( &
          species, x_left, x_right, xlt, xrt, riemann_solver, &
          final_x_flux(:, face_i, j), face_ok)
        if (.not. face_ok) return
      end do
    end do

    do face_j = 0, ny
      do i = 1, nx
        physical_face = (face_j == 0 .and. .not. periodic_y) .or. &
          (face_j == ny .and. .not. periodic_y)
        if (physical_face .and. &
            reactive_boundary_is_wall(boundaries%face(merge(3, 4, face_j == 0)))) then
          final_y_flux(:, i, face_j) = provisional_y_flux(:, i, face_j)
          cycle
        end if
        lj = face_j
        uj = face_j + 1
        if (face_j == 0) lj = ny
        if (face_j == ny) uj = 1
        if (use_transverse_correction .and. .not. physical_face) then
          call apply_reactive_transverse_correction_2d( &
            species, y_lower_base(:, i, face_j), y_lower_temperature(i, face_j), &
            provisional_x_flux(:, i, lj), provisional_x_flux(:, i - 1, lj), &
            0.5_dp * dt / dx, y_lower, ylt, theta_y_lower, correction_ok)
          if (.not. correction_ok) return
          call apply_reactive_transverse_correction_2d( &
            species, y_upper_base(:, i, face_j), y_upper_temperature(i, face_j), &
            provisional_x_flux(:, i, uj), provisional_x_flux(:, i - 1, uj), &
            0.5_dp * dt / dx, y_upper, yut, theta_y_upper, correction_ok)
          if (.not. correction_ok) return
        else
          y_lower = y_lower_base(:, i, face_j)
          y_upper = y_upper_base(:, i, face_j)
          ylt = y_lower_temperature(i, face_j)
          yut = y_upper_temperature(i, face_j)
          theta_y_lower = 1.0_dp
          theta_y_upper = 1.0_dp
        end if
        minimum_transverse_theta = min( &
          minimum_transverse_theta, theta_y_lower, theta_y_upper)
        call reactive_riemann_flux_y( &
          species, y_lower, y_upper, ylt, yut, riemann_solver, &
          final_y_flux(:, i, face_j), face_ok)
        if (.not. face_ok) return
      end do
    end do

    do j = 1, ny
      do i = 1, nx
        new_state(:, i, j) = state(:, i, j) - &
          dt / dx * (final_x_flux(:, i, j) - final_x_flux(:, i - 1, j)) - &
          dt / dy * (final_y_flux(:, i, j) - final_y_flux(:, i, j - 1))
        call recover_state_temperature( &
          species, new_state(:, i, j), temperature(i, j), &
          new_temperature(i, j), local_ok)
        if (.not. local_ok) return
      end do
    end do
    state = new_state
    temperature = new_temperature
    ok = .true.
  end subroutine advance_reactive_hydro_physical_2d

  subroutine advance_reactive_hydro_2d( &
      species, state, temperature, nx, ny, dx, dy, dt, reconstruction, &
      limiter, riemann_solver, use_transverse_correction, ok, &
      minimum_transverse_theta, ppm_contact_steepening, ppm_shock_flattening, &
      boundaries)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:, :, :), temperature(:, :)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: dx, dy, dt
    character(len=*), intent(in) :: reconstruction, limiter, riemann_solver
    logical, intent(in) :: use_transverse_correction
    logical, intent(out) :: ok
    real(dp), intent(out), optional :: minimum_transverse_theta
    logical, intent(in), optional :: ppm_contact_steepening
    logical, intent(in), optional :: ppm_shock_flattening
    type(reactive_boundary_set_2d), intent(in), optional :: boundaries

    logical :: use_contact, use_flattening, all_periodic
    real(dp) :: theta

    use_contact = .false.
    use_flattening = .false.
    if (present(ppm_contact_steepening)) use_contact = ppm_contact_steepening
    if (present(ppm_shock_flattening)) use_flattening = ppm_shock_flattening
    all_periodic = .true.
    if (present(boundaries)) then
      all_periodic = all([ &
        reactive_boundary_is_periodic(boundaries%face(1)), &
        reactive_boundary_is_periodic(boundaries%face(2)), &
        reactive_boundary_is_periodic(boundaries%face(3)), &
        reactive_boundary_is_periodic(boundaries%face(4))])
    end if
    if (all_periodic) then
      call advance_reactive_hydro_periodic_2d( &
        species, state, temperature, nx, ny, dx, dy, dt, reconstruction, &
        limiter, riemann_solver, use_transverse_correction, ok, theta, &
        use_contact, use_flattening)
    else
      call advance_reactive_hydro_physical_2d( &
        species, state, temperature, nx, ny, dx, dy, dt, reconstruction, &
        limiter, riemann_solver, use_transverse_correction, boundaries, ok, &
        theta, use_contact, use_flattening)
    end if
    if (present(minimum_transverse_theta)) minimum_transverse_theta = theta
  end subroutine advance_reactive_hydro_2d


  subroutine advance_reactive_chemistry_2d( &
      species, reactions, state, temperature, nx, ny, interval, rtol, atol, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(inout) :: state(:, :, :), temperature(:, :)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: interval, rtol, atol
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:), mass_fractions(:)
    real(dp) :: rho, kinetic_density, target_energy
    real(dp) :: elapsed, request, accepted, next_step, tolerance
    real(dp) :: checked_temperature, sound_speed
    logical :: local_ok
    integer :: i, j, k, substeps, nspecies

    ok = .false.
    if (interval < 0.0_dp) return
    if (interval <= tiny(1.0_dp)) then
      ok = .true.
      return
    end if
    nspecies = size(species)
    allocate(primitive(reactive_nprim(nspecies)), mass_fractions(nspecies))
    tolerance = 50.0_dp * epsilon(1.0_dp) * max(1.0_dp, interval)
    do j = 1, ny
      do i = 1, nx
        call reactive_conserved_to_primitive( &
          species, state(:, i, j), temperature(i, j), primitive, &
          checked_temperature, sound_speed, local_ok)
        if (.not. local_ok) return
        rho = state(irho, i, j)
        do k = 1, nspecies
          mass_fractions(k) = primitive(reactive_mass_fraction_component(k))
        end do
        kinetic_density = 0.5_dp * &
          (state(imx, i, j)**2 + state(imy, i, j)**2 + &
           state(imz, i, j)**2) / rho
        target_energy = (state(iet, i, j) - kinetic_density) / rho
        elapsed = 0.0_dp
        request = interval
        substeps = 0
        do while (elapsed < interval - tolerance)
          if (substeps >= max_chemistry_substeps) return
          request = min(request, interval - elapsed)
          call advance_constant_volume_adaptive( &
            species, reactions, rho, target_energy, request, rtol, atol, &
            mass_fractions, temperature(i, j), accepted, next_step, local_ok)
          if (.not. local_ok .or. accepted <= 0.0_dp) return
          elapsed = elapsed + accepted
          request = min(next_step, interval - elapsed)
          substeps = substeps + 1
        end do
        do k = 1, nspecies - 1
          state(reactive_species_component(k), i, j) = &
            rho * mass_fractions(k)
        end do
        state(reactive_species_component(nspecies), i, j) = rho - &
          sum(state(reactive_species_component(1): &
            reactive_species_component(nspecies - 1), i, j))
        call reactive_conserved_to_primitive( &
          species, state(:, i, j), temperature(i, j), primitive, &
          checked_temperature, sound_speed, local_ok)
        if (.not. local_ok) return
        temperature(i, j) = checked_temperature
      end do
    end do
    ok = .true.
  end subroutine advance_reactive_chemistry_2d

  subroutine advance_reactive_strang_2d( &
      species, reactions, state, temperature, nx, ny, dx, dy, dt, &
      reconstruction, limiter, riemann_solver, use_transverse_correction, &
      chemistry_enabled, rtol, atol, ok, minimum_transverse_theta, &
      ppm_contact_steepening, ppm_shock_flattening, transport, &
      transport_enabled, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, &
      minimum_transport_theta, boundaries)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(inout) :: state(:, :, :), temperature(:, :)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: dx, dy, dt, rtol, atol
    character(len=*), intent(in) :: reconstruction, limiter, riemann_solver
    logical, intent(in) :: use_transverse_correction, chemistry_enabled
    logical, intent(out) :: ok
    real(dp), intent(out), optional :: minimum_transverse_theta
    logical, intent(in), optional :: ppm_contact_steepening
    logical, intent(in), optional :: ppm_shock_flattening
    type(gas_transport_species), intent(in), optional :: transport(:)
    logical, intent(in), optional :: transport_enabled, viscosity_enabled
    logical, intent(in), optional :: thermal_conduction_enabled
    logical, intent(in), optional :: species_diffusion_enabled
    logical, intent(in), optional :: barodiffusion_enabled
    real(dp), intent(out), optional :: minimum_transport_theta
    type(reactive_boundary_set_2d), intent(in), optional :: boundaries

    logical :: local_ok, use_contact_steepening, use_shock_flattening
    logical :: use_transport, use_viscosity, use_conduction, use_diffusion
    logical :: use_barodiffusion
    real(dp) :: theta, transport_theta, stage_transport_theta

    ok = .false.
    theta = 1.0_dp
    transport_theta = 1.0_dp
    use_contact_steepening = .false.
    use_shock_flattening = .false.
    use_transport = .false.
    use_viscosity = .true.
    use_conduction = .true.
    use_diffusion = .true.
    use_barodiffusion = .true.
    if (present(ppm_contact_steepening)) &
      use_contact_steepening = ppm_contact_steepening
    if (present(ppm_shock_flattening)) &
      use_shock_flattening = ppm_shock_flattening
    if (present(transport_enabled)) use_transport = transport_enabled
    if (present(viscosity_enabled)) use_viscosity = viscosity_enabled
    if (present(thermal_conduction_enabled)) &
      use_conduction = thermal_conduction_enabled
    if (present(species_diffusion_enabled)) &
      use_diffusion = species_diffusion_enabled
    if (present(barodiffusion_enabled)) &
      use_barodiffusion = barodiffusion_enabled
    if (use_transport .and. .not. present(transport)) return
    if (chemistry_enabled) then
      call advance_reactive_chemistry_2d( &
        species, reactions, state, temperature, nx, ny, 0.5_dp * dt, &
        rtol, atol, local_ok)
      if (.not. local_ok) return
    end if
    if (use_transport) then
      if (present(boundaries)) then
        call advance_reactive_transport_2d( &
          species, transport, state, temperature, nx, ny, dx, dy, &
          0.5_dp * dt, use_viscosity, use_conduction, use_diffusion, &
          use_barodiffusion, stage_transport_theta, local_ok, boundaries)
      else
        call advance_reactive_transport_2d( &
          species, transport, state, temperature, nx, ny, dx, dy, &
          0.5_dp * dt, use_viscosity, use_conduction, use_diffusion, &
          use_barodiffusion, stage_transport_theta, local_ok)
      end if
      if (.not. local_ok) return
      transport_theta = min(transport_theta, stage_transport_theta)
    end if
    if (present(boundaries)) then
      call advance_reactive_hydro_2d( &
        species, state, temperature, nx, ny, dx, dy, dt, reconstruction, &
        limiter, riemann_solver, use_transverse_correction, local_ok, theta, &
        use_contact_steepening, use_shock_flattening, boundaries)
    else
      call advance_reactive_hydro_2d( &
        species, state, temperature, nx, ny, dx, dy, dt, reconstruction, &
        limiter, riemann_solver, use_transverse_correction, local_ok, theta, &
        use_contact_steepening, use_shock_flattening)
    end if
    if (.not. local_ok) return
    if (use_transport) then
      if (present(boundaries)) then
        call advance_reactive_transport_2d( &
          species, transport, state, temperature, nx, ny, dx, dy, &
          0.5_dp * dt, use_viscosity, use_conduction, use_diffusion, &
          use_barodiffusion, stage_transport_theta, local_ok, boundaries)
      else
        call advance_reactive_transport_2d( &
          species, transport, state, temperature, nx, ny, dx, dy, &
          0.5_dp * dt, use_viscosity, use_conduction, use_diffusion, &
          use_barodiffusion, stage_transport_theta, local_ok)
      end if
      if (.not. local_ok) return
      transport_theta = min(transport_theta, stage_transport_theta)
    end if
    if (chemistry_enabled) then
      call advance_reactive_chemistry_2d( &
        species, reactions, state, temperature, nx, ny, 0.5_dp * dt, &
        rtol, atol, local_ok)
      if (.not. local_ok) return
    end if
    if (present(minimum_transverse_theta)) minimum_transverse_theta = theta
    if (present(minimum_transport_theta)) &
      minimum_transport_theta = transport_theta
    ok = .true.
  end subroutine advance_reactive_strang_2d

  pure real(dp) function periodic_displacement(value, center, length) result(delta)
    real(dp), intent(in) :: value, center, length
    delta = value - center
    delta = delta - length * real(nint(delta / length), dp)
  end function periodic_displacement

  pure real(dp) function reactive_diagonal_wave_density( &
      x, y, time, config, base_density) result(density)
    real(dp), intent(in) :: x, y, time, base_density
    type(reactive_2d_config), intent(in) :: config
    real(dp) :: phase, lx, ly

    lx = config%x_upper - config%x_lower
    ly = config%y_upper - config%y_lower
    phase = 2.0_dp * pi * ( &
      (x - config%x_lower - config%initial_velocity_x * time) / lx + &
      (y - config%y_lower - config%initial_velocity_y * time) / ly)
    density = base_density * &
      (1.0_dp + config%density_wave_amplitude * sin(phase))
  end function reactive_diagonal_wave_density

  subroutine reactive_diagonal_composition_wave_exact( &
      species, x, y, time, config, density, mass_fractions, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: x, y, time
    type(reactive_2d_config), intent(in) :: config
    real(dp), intent(out) :: density, mass_fractions(:)
    logical, intent(out) :: ok

    real(dp) :: mole_fractions(7), phase, lx, ly

    density = 0.0_dp
    mass_fractions = 0.0_dp
    ok = .false.
    if (size(species) /= 7 .or. size(mass_fractions) /= 7) return
    lx = config%x_upper - config%x_lower
    ly = config%y_upper - config%y_lower
    if (lx <= 0.0_dp .or. ly <= 0.0_dp) return
    phase = sin(2.0_dp * pi * ( &
      (x - config%x_lower - config%initial_velocity_x * time) / lx + &
      (y - config%y_lower - config%initial_velocity_y * time) / ly))
    mole_fractions = [config%x_h2, config%x_h, config%x_o, config%x_o2, &
      config%x_oh, config%x_h2o, config%x_n2]
    mole_fractions(1) = mole_fractions(1) + &
      config%composition_wave_amplitude * phase
    mole_fractions(7) = mole_fractions(7) - &
      config%composition_wave_amplitude * phase
    call mass_fractions_from_mole_fractions( &
      species, mole_fractions, mass_fractions, ok)
    if (.not. ok) return
    density = mixture_density( &
      species, mass_fractions, config%initial_pressure, &
      config%initial_temperature, ok)
  end subroutine reactive_diagonal_composition_wave_exact

  subroutine initialize_reactive_2d( &
      species, config, state, temperature, dx, dy, base_density, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_2d_config), intent(in) :: config
    real(dp), allocatable, intent(out) :: state(:, :, :), temperature(:, :)
    real(dp), intent(out) :: dx, dy, base_density
    logical, intent(out) :: ok

    real(dp), allocatable :: mole_fractions(:), mass_fractions(:), primitive(:)
    real(dp), allocatable :: local_mass_fractions(:)
    real(dp) :: x, y, rho, pressure, u, v, w, local_temperature
    real(dp) :: sound_speed, rx, ry, radius_squared, factor
    logical :: local_ok
    integer :: i, j, k, nvar, nprimitive

    ok = .false.
    dx = (config%x_upper - config%x_lower) / real(config%nx, dp)
    dy = (config%y_upper - config%y_lower) / real(config%ny, dp)
    if (dx <= 0.0_dp .or. dy <= 0.0_dp .or. size(species) /= 7) return
    nvar = reactive_nvar(size(species))
    nprimitive = reactive_nprim(size(species))
    allocate(state(nvar, config%nx, config%ny))
    allocate(temperature(config%nx, config%ny))
    allocate(mole_fractions(size(species)), mass_fractions(size(species)))
    allocate(local_mass_fractions(size(species)))
    allocate(primitive(nprimitive))
    mole_fractions = [config%x_h2, config%x_h, config%x_o, config%x_o2, &
      config%x_oh, config%x_h2o, config%x_n2]
    call mass_fractions_from_mole_fractions( &
      species, mole_fractions, mass_fractions, local_ok)
    if (.not. local_ok) return
    base_density = mixture_density( &
      species, mass_fractions, config%initial_pressure, &
      config%initial_temperature, local_ok)
    if (.not. local_ok) return

    do j = 1, config%ny
      y = config%y_lower + (real(j, dp) - 0.5_dp) * dy
      do i = 1, config%nx
        x = config%x_lower + (real(i, dp) - 0.5_dp) * dx
        local_mass_fractions = mass_fractions
        pressure = config%initial_pressure
        u = config%initial_velocity_x
        v = config%initial_velocity_y
        w = 0.0_dp
        select case (trim(config%problem))
        case ("diagonal_wave")
          rho = reactive_diagonal_wave_density(x, y, 0.0_dp, config, base_density)
        case ("diagonal_composition_wave")
          call reactive_diagonal_composition_wave_exact( &
            species, x, y, 0.0_dp, config, rho, local_mass_fractions, local_ok)
          if (.not. local_ok) return
        case ("reactive_vortex")
          rx = periodic_displacement( &
            x, config%vortex_center_x, config%x_upper - config%x_lower)
          ry = periodic_displacement( &
            y, config%vortex_center_y, config%y_upper - config%y_lower)
          radius_squared = (rx * rx + ry * ry) / &
            (config%vortex_radius * config%vortex_radius)
          factor = exp(-0.5_dp * radius_squared)
          u = u - config%vortex_strength * ry / config%vortex_radius * factor
          v = v + config%vortex_strength * rx / config%vortex_radius * factor
          rho = base_density
        case ("reactive_hotspot")
          rx = periodic_displacement( &
            x, config%hotspot_center_x, config%x_upper - config%x_lower)
          ry = periodic_displacement( &
            y, config%hotspot_center_y, config%y_upper - config%y_lower)
          local_temperature = config%initial_temperature + &
            config%hotspot_temperature_rise * exp(-0.5_dp * &
            (rx * rx + ry * ry) / (config%hotspot_width**2))
          rho = mixture_density( &
            species, local_mass_fractions, pressure, local_temperature, local_ok)
          if (.not. local_ok) return
        case ("uniform_reactor")
          rho = base_density
        case ("couette_channel")
          u = config%wall_velocity_y_lower(1) + &
            (config%wall_velocity_y_upper(1) - &
             config%wall_velocity_y_lower(1)) * &
            (y - config%y_lower) / (config%y_upper - config%y_lower)
          v = 0.0_dp
          rho = base_density
        case ("thermal_channel")
          local_temperature = config%wall_temperature_y_lower + &
            (config%wall_temperature_y_upper - &
             config%wall_temperature_y_lower) * &
            (y - config%y_lower) / (config%y_upper - config%y_lower)
          rho = mixture_density( &
            species, local_mass_fractions, pressure, local_temperature, local_ok)
          if (.not. local_ok) return
        case ("inflow_outflow")
          rho = base_density
        case default
          return
        end select
        primitive(1:5) = [rho, u, v, w, pressure]
        do k = 1, size(species)
          primitive(reactive_mass_fraction_component(k)) = &
            local_mass_fractions(k)
        end do
        call reactive_primitive_to_conserved( &
          species, primitive, state(:, i, j), temperature(i, j), &
          sound_speed, local_ok)
        if (.not. local_ok) return
      end do
    end do
    ok = .true.
  end subroutine initialize_reactive_2d

  subroutine reactive_integrals_2d(state, nx, ny, dx, dy, integrals)
    real(dp), intent(in) :: state(:, :, :)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: dx, dy
    real(dp), intent(out) :: integrals(5)

    integrals(1) = sum(state(irho, 1:nx, 1:ny)) * dx * dy
    integrals(2) = sum(state(imx, 1:nx, 1:ny)) * dx * dy
    integrals(3) = sum(state(imy, 1:nx, 1:ny)) * dx * dy
    integrals(4) = sum(state(imz, 1:nx, 1:ny)) * dx * dy
    integrals(5) = sum(state(iet, 1:nx, 1:ny)) * dx * dy
  end subroutine reactive_integrals_2d

  subroutine reactive_extrema_2d( &
      species, state, temperature, nx, ny, minimum_density, maximum_density, &
      minimum_pressure, maximum_pressure, minimum_temperature, &
      maximum_temperature, maximum_speed, maximum_closure_error, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    integer, intent(in) :: nx, ny
    real(dp), intent(out) :: minimum_density, maximum_density
    real(dp), intent(out) :: minimum_pressure, maximum_pressure
    real(dp), intent(out) :: minimum_temperature, maximum_temperature
    real(dp), intent(out) :: maximum_speed, maximum_closure_error
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: local_temperature, sound_speed, closure, speed
    logical :: local_ok
    integer :: i, j, k

    minimum_density = huge(1.0_dp)
    maximum_density = -huge(1.0_dp)
    minimum_pressure = huge(1.0_dp)
    maximum_pressure = -huge(1.0_dp)
    minimum_temperature = huge(1.0_dp)
    maximum_temperature = -huge(1.0_dp)
    maximum_speed = 0.0_dp
    maximum_closure_error = 0.0_dp
    ok = .false.
    allocate(primitive(reactive_nprim(size(species))))
    do j = 1, ny
      do i = 1, nx
        call reactive_conserved_to_primitive( &
          species, state(:, i, j), temperature(i, j), primitive, &
          local_temperature, sound_speed, local_ok)
        if (.not. local_ok) return
        minimum_density = min(minimum_density, primitive(1))
        maximum_density = max(maximum_density, primitive(1))
        minimum_pressure = min(minimum_pressure, primitive(5))
        maximum_pressure = max(maximum_pressure, primitive(5))
        minimum_temperature = min(minimum_temperature, local_temperature)
        maximum_temperature = max(maximum_temperature, local_temperature)
        speed = sqrt(primitive(2)**2 + primitive(3)**2 + primitive(4)**2)
        maximum_speed = max(maximum_speed, speed)
        closure = 0.0_dp
        do k = 1, size(species)
          closure = closure + &
            primitive(reactive_mass_fraction_component(k))
        end do
        maximum_closure_error = max(maximum_closure_error, abs(closure - 1.0_dp))
      end do
    end do
    ok = .true.
  end subroutine reactive_extrema_2d

  subroutine simulate_reactive_2d( &
      species, reactions, config, state, temperature, dx, dy, time, steps, &
      initial_integrals, final_integrals, minimum_transverse_theta, &
      base_density, ok, transport, minimum_transport_theta)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_2d_config), intent(in) :: config
    real(dp), allocatable, intent(out) :: state(:, :, :), temperature(:, :)
    real(dp), intent(out) :: dx, dy, time
    integer, intent(out) :: steps
    real(dp), intent(out) :: initial_integrals(5), final_integrals(5)
    real(dp), intent(out) :: minimum_transverse_theta, base_density
    logical, intent(out) :: ok
    type(gas_transport_species), intent(in), optional :: transport(:)
    real(dp), intent(out), optional :: minimum_transport_theta

    real(dp) :: dt, step_theta, transport_dt, maximum_diffusivity
    real(dp) :: step_transport_theta, local_minimum_transport_theta
    logical :: local_ok
    type(reactive_boundary_set_2d) :: boundaries

    call initialize_reactive_2d( &
      species, config, state, temperature, dx, dy, base_density, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call build_reactive_boundary_set_2d(species, config, boundaries, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call reactive_integrals_2d( &
      state, config%nx, config%ny, dx, dy, initial_integrals)
    time = 0.0_dp
    steps = 0
    minimum_transverse_theta = 1.0_dp
    local_minimum_transport_theta = 1.0_dp
    if (config%transport_enabled .and. .not. present(transport)) then
      ok = .false.
      return
    end if
    do while (time < config%final_time)
      if (steps >= config%maximum_steps) then
        ok = .false.
        return
      end if
      call compute_reactive_cfl_timestep_2d( &
        species, state, temperature, config%nx, config%ny, dx, dy, &
        config%cfl, dt, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      if (config%transport_enabled) then
        call reactive_transport_timestep_2d( &
          species, transport, state, temperature, config%nx, config%ny, &
          dx, dy, config%transport_cfl, config%viscosity_enabled, &
          config%thermal_conduction_enabled, config%species_diffusion_enabled, &
          transport_dt, maximum_diffusivity, local_ok)
        if (.not. local_ok) then
          ok = .false.
          return
        end if
        dt = min(dt, transport_dt)
      end if
      dt = min(dt, config%final_time - time)
      if (config%transport_enabled) then
        call advance_reactive_strang_2d( &
          species, reactions, state, temperature, config%nx, config%ny, &
          dx, dy, dt, config%reconstruction, config%limiter, &
          config%riemann_solver, config%use_transverse_correction, &
          config%chemistry_enabled, config%chemistry_relative_tolerance, &
          config%chemistry_absolute_tolerance, local_ok, step_theta, &
          config%ppm_contact_steepening, config%ppm_shock_flattening, &
          transport, config%transport_enabled, config%viscosity_enabled, &
          config%thermal_conduction_enabled, config%species_diffusion_enabled, &
          config%barodiffusion_enabled, step_transport_theta, boundaries)
      else
        call advance_reactive_strang_2d( &
        species, reactions, state, temperature, config%nx, config%ny, &
        dx, dy, dt, config%reconstruction, config%limiter, &
        config%riemann_solver, config%use_transverse_correction, &
        config%chemistry_enabled, config%chemistry_relative_tolerance, &
        config%chemistry_absolute_tolerance, local_ok, step_theta, &
        config%ppm_contact_steepening, config%ppm_shock_flattening, &
        boundaries=boundaries)
        step_transport_theta = 1.0_dp
      end if
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      minimum_transverse_theta = min(minimum_transverse_theta, step_theta)
      local_minimum_transport_theta = min( &
        local_minimum_transport_theta, step_transport_theta)
      time = time + dt
      steps = steps + 1
    end do
    call reactive_integrals_2d( &
      state, config%nx, config%ny, dx, dy, final_integrals)
    if (present(minimum_transport_theta)) &
      minimum_transport_theta = local_minimum_transport_theta
    ok = .true.
  end subroutine simulate_reactive_2d

  subroutine write_reactive_2d_csv( &
      path, species, config, state, temperature, dx, dy, time, ok)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_2d_config), intent(in) :: config
    real(dp), intent(in) :: state(:, :, :), temperature(:, :), dx, dy, time
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: x, y, local_temperature, sound_speed
    logical :: local_ok
    integer :: unit, status, i, j, k

    ok = .false.
    if (size(species) /= 7) return
    allocate(primitive(reactive_nprim(size(species))))
    open(newunit=unit, file=trim(path), status="replace", action="write", &
      iostat=status)
    if (status /= 0) return
    write(unit, '(a)') &
      "time,x,y,rho,u,v,w,pressure,temperature,rhoE," // &
      "Y_H2,Y_H,Y_O,Y_O2,Y_OH,Y_H2O,Y_N2"
    do j = 1, config%ny
      y = config%y_lower + (real(j, dp) - 0.5_dp) * dy
      do i = 1, config%nx
        x = config%x_lower + (real(i, dp) - 0.5_dp) * dx
        call reactive_conserved_to_primitive( &
          species, state(:, i, j), temperature(i, j), primitive, &
          local_temperature, sound_speed, local_ok)
        if (.not. local_ok) then
          close(unit)
          return
        end if
        write(unit, '(*(es25.16e3,:,","))') time, x, y, state(irho, i, j), &
          primitive(2), primitive(3), primitive(4), primitive(5), &
          local_temperature, state(iet, i, j), &
          (primitive(reactive_mass_fraction_component(k)), &
            k = 1, size(species))
      end do
    end do
    close(unit)
    ok = .true.
  end subroutine write_reactive_2d_csv

end module reactive_2d_mod
