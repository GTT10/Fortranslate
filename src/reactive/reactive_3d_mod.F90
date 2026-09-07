module reactive_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use constants_mod, only: density_floor, pressure_floor
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved, reactive_conserved_to_primitive, &
    reactive_riemann_flux_x, characteristic_limited_slope
  use reactive_2d_mod, only: &
    reactive_riemann_flux_y, advance_reactive_chemistry_2d
  use reactive_directional_flux_3d_mod, only: reactive_riemann_flux_z
  use reactive_transport_3d_mod, only: advance_reactive_transport_3d
  implicit none
  private

  public :: compute_reactive_cfl_timestep_3d
  public :: compute_reactive_euler_rhs_3d
  public :: advance_reactive_euler_ssprk2_3d
  public :: compute_reactive_face_fluxes_3d
  public :: compute_reactive_plm_face_fluxes_3d
  public :: compute_reactive_plm_slab_face_fluxes_3d
  public :: advance_reactive_euler_ssprk2_with_fluxes_3d
  public :: advance_reactive_euler_ssprk2_plm_3d
  public :: advance_reactive_euler_ssprk2_plm_with_fluxes_3d
  public :: advance_reactive_chemistry_3d
  public :: advance_reactive_strang_3d
  public :: advance_reactive_full_3d
  public :: recover_reactive_temperatures_3d
  public :: reactive_integrals_3d
  public :: reactive_extrema_3d

contains

  subroutine compute_reactive_cfl_timestep_3d( &
      species, state, temperature, nx, ny, nz, dx, dy, dz, cfl, dt, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz, cfl
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: local_temperature, sound_speed, rate, maximum_rate
    logical :: cell_ok
    integer :: i, j, k, nvar

    dt = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return
    if (.not. all(ieee_is_finite([dx, dy, dz, cfl]))) return
    if (nx < 2 .or. ny < 2 .or. nz < 2 .or. &
        dx <= 0.0_dp .or. dy <= 0.0_dp .or. dz <= 0.0_dp .or. &
        cfl <= 0.0_dp .or. cfl > 1.0_dp) return

    allocate(primitive(reactive_nprim(size(species))))
    maximum_rate = 0.0_dp
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), primitive, &
            local_temperature, sound_speed, cell_ok)
          if (.not. cell_ok) return
          rate = (abs(primitive(2)) + sound_speed) / dx + &
            (abs(primitive(3)) + sound_speed) / dy + &
            (abs(primitive(4)) + sound_speed) / dz
          maximum_rate = max(maximum_rate, rate)
        end do
      end do
    end do
    if (.not. ieee_is_finite(maximum_rate)) return
    if (maximum_rate <= 0.0_dp) return
    dt = cfl / maximum_rate
    ok = ieee_is_finite(dt)
    if (ok) ok = dt > 0.0_dp
  end subroutine compute_reactive_cfl_timestep_3d

  subroutine compute_reactive_euler_rhs_3d( &
      species, state, temperature, nx, ny, nz, dx, dy, dz, &
      riemann_solver, rhs, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz
    character(len=*), intent(in) :: riemann_solver
    real(dp), intent(out) :: rhs(:, :, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: flux_x(:, :, :, :)
    real(dp), allocatable :: flux_y(:, :, :, :)
    real(dp), allocatable :: flux_z(:, :, :, :)
    logical :: face_ok
    integer :: i, j, k, next_i, next_j, next_k
    integer :: previous_i, previous_j, previous_k, nvar

    rhs = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return
    if (size(rhs, 1) /= nvar .or. size(rhs, 2) /= nx .or. &
        size(rhs, 3) /= ny .or. size(rhs, 4) /= nz) return
    if (.not. all(ieee_is_finite([dx, dy, dz]))) return
    if (nx < 2 .or. ny < 2 .or. nz < 2 .or. &
        dx <= 0.0_dp .or. dy <= 0.0_dp .or. dz <= 0.0_dp) return

    allocate(flux_x(nvar, nx, ny, nz))
    allocate(flux_y(nvar, nx, ny, nz))
    allocate(flux_z(nvar, nx, ny, nz))

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          next_i = modulo(i, nx) + 1
          call reactive_riemann_flux_x( &
            species, state(:, i, j, k), state(:, next_i, j, k), &
            temperature(i, j, k), temperature(next_i, j, k), &
            riemann_solver, flux_x(:, i, j, k), face_ok)
          if (.not. face_ok) return
        end do
      end do
    end do

    do k = 1, nz
      do j = 1, ny
        next_j = modulo(j, ny) + 1
        do i = 1, nx
          call reactive_riemann_flux_y( &
            species, state(:, i, j, k), state(:, i, next_j, k), &
            temperature(i, j, k), temperature(i, next_j, k), &
            riemann_solver, flux_y(:, i, j, k), face_ok)
          if (.not. face_ok) return
        end do
      end do
    end do

    do k = 1, nz
      next_k = modulo(k, nz) + 1
      do j = 1, ny
        do i = 1, nx
          call reactive_riemann_flux_z( &
            species, state(:, i, j, k), state(:, i, j, next_k), &
            temperature(i, j, k), temperature(i, j, next_k), &
            riemann_solver, flux_z(:, i, j, k), face_ok)
          if (.not. face_ok) return
        end do
      end do
    end do

    do k = 1, nz
      previous_k = modulo(k - 2, nz) + 1
      do j = 1, ny
        previous_j = modulo(j - 2, ny) + 1
        do i = 1, nx
          previous_i = modulo(i - 2, nx) + 1
          rhs(:, i, j, k) = &
            -(flux_x(:, i, j, k) - &
              flux_x(:, previous_i, j, k)) / dx &
            -(flux_y(:, i, j, k) - &
              flux_y(:, i, previous_j, k)) / dy &
            -(flux_z(:, i, j, k) - &
              flux_z(:, i, j, previous_k)) / dz
        end do
      end do
    end do
    ok = all(ieee_is_finite(rhs))
  end subroutine compute_reactive_euler_rhs_3d

  subroutine advance_reactive_euler_ssprk2_3d( &
      species, state, temperature, nx, ny, nz, dx, dy, dz, dt, &
      riemann_solver, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz, dt
    character(len=*), intent(in) :: riemann_solver
    logical, intent(out) :: ok

    real(dp), allocatable :: old_state(:, :, :, :)
    real(dp), allocatable :: stage_state(:, :, :, :)
    real(dp), allocatable :: updated_state(:, :, :, :)
    real(dp), allocatable :: rhs(:, :, :, :)
    real(dp), allocatable :: old_temperature(:, :, :)
    real(dp), allocatable :: stage_temperature(:, :, :)
    real(dp), allocatable :: updated_temperature(:, :, :)
    logical :: local_ok
    integer :: nvar

    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return
    if (.not. ieee_is_finite(dt)) return
    if (nx < 2 .or. ny < 2 .or. nz < 2 .or. dt <= 0.0_dp) return

    allocate(old_state(nvar, nx, ny, nz))
    allocate(stage_state(nvar, nx, ny, nz))
    allocate(updated_state(nvar, nx, ny, nz))
    allocate(rhs(nvar, nx, ny, nz))
    allocate(old_temperature(nx, ny, nz))
    allocate(stage_temperature(nx, ny, nz))
    allocate(updated_temperature(nx, ny, nz))
    old_state = state
    call recover_reactive_temperatures_3d( &
      species, old_state, temperature, nx, ny, nz, old_temperature, local_ok)
    if (.not. local_ok) return

    call compute_reactive_euler_rhs_3d( &
      species, old_state, old_temperature, nx, ny, nz, dx, dy, dz, &
      riemann_solver, rhs, local_ok)
    if (.not. local_ok) return
    stage_state = old_state + dt * rhs
    call recover_reactive_temperatures_3d( &
      species, stage_state, old_temperature, nx, ny, nz, &
      stage_temperature, local_ok)
    if (.not. local_ok) return

    call compute_reactive_euler_rhs_3d( &
      species, stage_state, stage_temperature, nx, ny, nz, dx, dy, dz, &
      riemann_solver, rhs, local_ok)
    if (.not. local_ok) return
    updated_state = 0.5_dp * old_state + &
      0.5_dp * (stage_state + dt * rhs)
    call recover_reactive_temperatures_3d( &
      species, updated_state, stage_temperature, nx, ny, nz, &
      updated_temperature, local_ok)
    if (.not. local_ok) return

    state = updated_state
    temperature = updated_temperature
    ok = .true.
  end subroutine advance_reactive_euler_ssprk2_3d

  subroutine compute_reactive_face_fluxes_3d( &
      species, state, temperature, nx, ny, nz, riemann_solver, &
      flux_x, flux_y, flux_z, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    character(len=*), intent(in) :: riemann_solver
    real(dp), intent(out) :: flux_x(:, :, :, :)
    real(dp), intent(out) :: flux_y(:, :, :, :)
    real(dp), intent(out) :: flux_z(:, :, :, :)
    logical, intent(out) :: ok

    logical :: face_ok
    integer :: i, j, k, next_i, next_j, next_k, nvar

    flux_x = 0.0_dp
    flux_y = 0.0_dp
    flux_z = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return
    if (nx < 2 .or. ny < 2 .or. nz < 2 .or. &
        any(shape(flux_x) /= shape(state)) .or. &
        any(shape(flux_y) /= shape(state)) .or. &
        any(shape(flux_z) /= shape(state))) return

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          next_i = modulo(i, nx) + 1
          call reactive_riemann_flux_x( &
            species, state(:, i, j, k), state(:, next_i, j, k), &
            temperature(i, j, k), temperature(next_i, j, k), &
            riemann_solver, flux_x(:, i, j, k), face_ok)
          if (.not. face_ok) return
          next_j = modulo(j, ny) + 1
          call reactive_riemann_flux_y( &
            species, state(:, i, j, k), state(:, i, next_j, k), &
            temperature(i, j, k), temperature(i, next_j, k), &
            riemann_solver, flux_y(:, i, j, k), face_ok)
          if (.not. face_ok) return
          next_k = modulo(k, nz) + 1
          call reactive_riemann_flux_z( &
            species, state(:, i, j, k), state(:, i, j, next_k), &
            temperature(i, j, k), temperature(i, j, next_k), &
            riemann_solver, flux_z(:, i, j, k), face_ok)
          if (.not. face_ok) return
        end do
      end do
    end do
    ok = all(ieee_is_finite(flux_x)) .and. &
      all(ieee_is_finite(flux_y)) .and. all(ieee_is_finite(flux_z))
  end subroutine compute_reactive_face_fluxes_3d

  subroutine compute_reactive_plm_face_fluxes_3d( &
      species, state, temperature, nx, ny, nz, limiter, riemann_solver, &
      flux_x, flux_y, flux_z, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    character(len=*), intent(in) :: limiter, riemann_solver
    real(dp), intent(out) :: flux_x(:, :, :, :)
    real(dp), intent(out) :: flux_y(:, :, :, :)
    real(dp), intent(out) :: flux_z(:, :, :, :)
    logical, intent(out) :: ok

    call compute_reactive_plm_slab_face_fluxes_3d( &
      species, state, temperature, nx, ny, nz, 1, nx, limiter, &
      riemann_solver, flux_x, flux_y, flux_z, ok)
  end subroutine compute_reactive_plm_face_fluxes_3d

  subroutine compute_reactive_plm_slab_face_fluxes_3d( &
      species, state, temperature, nx, ny, nz, first_i, last_i, &
      limiter, riemann_solver, flux_x, flux_y, flux_z, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz, first_i, last_i
    character(len=*), intent(in) :: limiter, riemann_solver
    real(dp), intent(out) :: flux_x(:, :, :, :)
    real(dp), intent(out) :: flux_y(:, :, :, :)
    real(dp), intent(out) :: flux_z(:, :, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:, :, :, :), sound_speed(:, :, :)
    real(dp), allocatable :: x_minus(:, :, :, :), x_plus(:, :, :, :)
    real(dp), allocatable :: y_minus(:, :, :, :), y_plus(:, :, :, :)
    real(dp), allocatable :: z_minus(:, :, :, :), z_plus(:, :, :, :)
    real(dp), allocatable :: candidate_x(:, :, :, :)
    real(dp), allocatable :: candidate_y(:, :, :, :)
    real(dp), allocatable :: candidate_z(:, :, :, :)
    real(dp), allocatable :: center(:), dl(:), dr(:), slope(:)
    real(dp), allocatable :: rotated_center(:), rotated_dl(:), rotated_dr(:)
    real(dp), allocatable :: rotated_slope(:), rotated_minus(:)
    real(dp), allocatable :: rotated_plus(:)
    real(dp), allocatable :: left_state(:), right_state(:)
    logical, allocatable :: primitive_needed(:), reconstruction_needed(:)
    real(dp) :: recovered_temperature, theta
    real(dp) :: left_temperature, right_temperature
    logical :: local_ok
    integer :: i, j, k, im, ip, jm, jp, km, kp
    integer :: nvar, nprimitive

    flux_x = 0.0_dp
    flux_y = 0.0_dp
    flux_z = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    nprimitive = reactive_nprim(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return
    if (nx < 3 .or. ny < 3 .or. nz < 3 .or. &
        any(shape(flux_x) /= shape(state)) .or. &
        any(shape(flux_y) /= shape(state)) .or. &
        any(shape(flux_z) /= shape(state))) return
    if (first_i < 1 .or. last_i > nx .or. last_i < first_i) return
    if (trim(limiter) /= "minmod" .and. trim(limiter) /= "mc") return

    allocate(primitive(nprimitive, nx, ny, nz), sound_speed(nx, ny, nz))
    allocate(x_minus(nprimitive, nx, ny, nz))
    allocate(x_plus(nprimitive, nx, ny, nz))
    allocate(y_minus(nprimitive, nx, ny, nz))
    allocate(y_plus(nprimitive, nx, ny, nz))
    allocate(z_minus(nprimitive, nx, ny, nz))
    allocate(z_plus(nprimitive, nx, ny, nz))
    allocate(candidate_x(nvar, nx, ny, nz))
    allocate(candidate_y(nvar, nx, ny, nz))
    allocate(candidate_z(nvar, nx, ny, nz))
    allocate(center(nprimitive), dl(nprimitive), dr(nprimitive))
    allocate(slope(nprimitive), rotated_center(nprimitive))
    allocate(rotated_dl(nprimitive), rotated_dr(nprimitive))
    allocate(rotated_slope(nprimitive), rotated_minus(nprimitive))
    allocate(rotated_plus(nprimitive))
    allocate(left_state(nvar), right_state(nvar))
    allocate(primitive_needed(nx), reconstruction_needed(nx))

    primitive_needed = .false.
    reconstruction_needed = .false.
    do i = first_i, last_i
      reconstruction_needed(i) = .true.
      reconstruction_needed(periodic_index_3d(i + 1, nx)) = .true.
    end do
    do i = 1, nx
      if (.not. reconstruction_needed(i)) cycle
      primitive_needed(i) = .true.
      primitive_needed(periodic_index_3d(i - 1, nx)) = .true.
      primitive_needed(periodic_index_3d(i + 1, nx)) = .true.
    end do

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (.not. primitive_needed(i)) cycle
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), &
            primitive(:, i, j, k), recovered_temperature, &
            sound_speed(i, j, k), local_ok)
          if (.not. local_ok) return
        end do
      end do
    end do

    do k = 1, nz
      km = periodic_index_3d(k - 1, nz)
      kp = periodic_index_3d(k + 1, nz)
      do j = 1, ny
        jm = periodic_index_3d(j - 1, ny)
        jp = periodic_index_3d(j + 1, ny)
        do i = 1, nx
          if (.not. reconstruction_needed(i)) cycle
          im = periodic_index_3d(i - 1, nx)
          ip = periodic_index_3d(i + 1, nx)
          center = primitive(:, i, j, k)
          dl = center - primitive(:, im, j, k)
          dr = primitive(:, ip, j, k) - center
          call characteristic_limited_slope( &
            center, dl, dr, sound_speed(i, j, k), limiter, slope, local_ok)
          if (.not. local_ok) return
          theta = primitive_slope_scale_3d(center, slope, size(species))
          slope = theta * slope
          x_minus(:, i, j, k) = center - 0.5_dp * slope
          x_plus(:, i, j, k) = center + 0.5_dp * slope
          call sanitize_primitive_3d( &
            x_minus(:, i, j, k), center, size(species))
          call sanitize_primitive_3d( &
            x_plus(:, i, j, k), center, size(species))

          call rotate_primitive_y_to_x_3d(center, rotated_center)
          call rotate_primitive_y_to_x_3d( &
            center - primitive(:, i, jm, k), rotated_dl)
          call rotate_primitive_y_to_x_3d( &
            primitive(:, i, jp, k) - center, rotated_dr)
          call characteristic_limited_slope( &
            rotated_center, rotated_dl, rotated_dr, sound_speed(i, j, k), &
            limiter, rotated_slope, local_ok)
          if (.not. local_ok) return
          theta = primitive_slope_scale_3d( &
            rotated_center, rotated_slope, size(species))
          rotated_slope = theta * rotated_slope
          rotated_minus = rotated_center - 0.5_dp * rotated_slope
          rotated_plus = rotated_center + 0.5_dp * rotated_slope
          call rotate_primitive_x_to_y_3d( &
            rotated_minus, y_minus(:, i, j, k))
          call rotate_primitive_x_to_y_3d( &
            rotated_plus, y_plus(:, i, j, k))
          call sanitize_primitive_3d( &
            y_minus(:, i, j, k), center, size(species))
          call sanitize_primitive_3d( &
            y_plus(:, i, j, k), center, size(species))

          call rotate_primitive_z_to_x_3d(center, rotated_center)
          call rotate_primitive_z_to_x_3d( &
            center - primitive(:, i, j, km), rotated_dl)
          call rotate_primitive_z_to_x_3d( &
            primitive(:, i, j, kp) - center, rotated_dr)
          call characteristic_limited_slope( &
            rotated_center, rotated_dl, rotated_dr, sound_speed(i, j, k), &
            limiter, rotated_slope, local_ok)
          if (.not. local_ok) return
          theta = primitive_slope_scale_3d( &
            rotated_center, rotated_slope, size(species))
          rotated_slope = theta * rotated_slope
          rotated_minus = rotated_center - 0.5_dp * rotated_slope
          rotated_plus = rotated_center + 0.5_dp * rotated_slope
          call rotate_primitive_x_to_z_3d( &
            rotated_minus, z_minus(:, i, j, k))
          call rotate_primitive_x_to_z_3d( &
            rotated_plus, z_plus(:, i, j, k))
          call sanitize_primitive_3d( &
            z_minus(:, i, j, k), center, size(species))
          call sanitize_primitive_3d( &
            z_plus(:, i, j, k), center, size(species))
        end do
      end do
    end do

    candidate_x = 0.0_dp
    candidate_y = 0.0_dp
    candidate_z = 0.0_dp
    do k = 1, nz
      kp = periodic_index_3d(k + 1, nz)
      do j = 1, ny
        jp = periodic_index_3d(j + 1, ny)
        do i = first_i, last_i
          ip = periodic_index_3d(i + 1, nx)
          call primitive_face_to_state_3d( &
            species, x_plus(:, i, j, k), primitive(:, i, j, k), &
            left_state, left_temperature, local_ok)
          if (.not. local_ok) return
          call primitive_face_to_state_3d( &
            species, x_minus(:, ip, j, k), primitive(:, ip, j, k), &
            right_state, right_temperature, local_ok)
          if (.not. local_ok) return
          call reactive_riemann_flux_x( &
            species, left_state, right_state, left_temperature, &
            right_temperature, riemann_solver, candidate_x(:, i, j, k), &
            local_ok)
          if (.not. local_ok) return

          call primitive_face_to_state_3d( &
            species, y_plus(:, i, j, k), primitive(:, i, j, k), &
            left_state, left_temperature, local_ok)
          if (.not. local_ok) return
          call primitive_face_to_state_3d( &
            species, y_minus(:, i, jp, k), primitive(:, i, jp, k), &
            right_state, right_temperature, local_ok)
          if (.not. local_ok) return
          call reactive_riemann_flux_y( &
            species, left_state, right_state, left_temperature, &
            right_temperature, riemann_solver, candidate_y(:, i, j, k), &
            local_ok)
          if (.not. local_ok) return

          call primitive_face_to_state_3d( &
            species, z_plus(:, i, j, k), primitive(:, i, j, k), &
            left_state, left_temperature, local_ok)
          if (.not. local_ok) return
          call primitive_face_to_state_3d( &
            species, z_minus(:, i, j, kp), primitive(:, i, j, kp), &
            right_state, right_temperature, local_ok)
          if (.not. local_ok) return
          call reactive_riemann_flux_z( &
            species, left_state, right_state, left_temperature, &
            right_temperature, riemann_solver, candidate_z(:, i, j, k), &
            local_ok)
          if (.not. local_ok) return
        end do
      end do
    end do
    if (.not. all(ieee_is_finite(candidate_x)) .or. &
        .not. all(ieee_is_finite(candidate_y)) .or. &
        .not. all(ieee_is_finite(candidate_z))) return
    flux_x = candidate_x
    flux_y = candidate_y
    flux_z = candidate_z
    ok = .true.
  end subroutine compute_reactive_plm_slab_face_fluxes_3d

  pure integer function periodic_index_3d(index, extent) result(wrapped)
    integer, intent(in) :: index, extent

    wrapped = 1 + modulo(index - 1, extent)
  end function periodic_index_3d

  pure subroutine rotate_primitive_y_to_x_3d(input, output)
    real(dp), intent(in) :: input(:)
    real(dp), intent(out) :: output(:)

    output = input
    if (size(input) /= size(output) .or. size(input) < 5) return
    output(2) = input(3)
    output(3) = input(2)
  end subroutine rotate_primitive_y_to_x_3d

  pure subroutine rotate_primitive_x_to_y_3d(input, output)
    real(dp), intent(in) :: input(:)
    real(dp), intent(out) :: output(:)

    call rotate_primitive_y_to_x_3d(input, output)
  end subroutine rotate_primitive_x_to_y_3d

  pure subroutine rotate_primitive_z_to_x_3d(input, output)
    real(dp), intent(in) :: input(:)
    real(dp), intent(out) :: output(:)

    output = input
    if (size(input) /= size(output) .or. size(input) < 5) return
    output(2) = input(4)
    output(4) = input(2)
  end subroutine rotate_primitive_z_to_x_3d

  pure subroutine rotate_primitive_x_to_z_3d(input, output)
    real(dp), intent(in) :: input(:)
    real(dp), intent(out) :: output(:)

    call rotate_primitive_z_to_x_3d(input, output)
  end subroutine rotate_primitive_x_to_z_3d

  pure real(dp) function primitive_lower_scale_3d( &
      center, slope, lower) result(theta)
    real(dp), intent(in) :: center, slope, lower
    real(dp) :: magnitude

    magnitude = abs(slope)
    if (magnitude <= tiny(1.0_dp) .or. center - magnitude > lower) then
      theta = 1.0_dp
    else
      theta = max(0.0_dp, min(1.0_dp, (center - lower) / magnitude))
    end if
  end function primitive_lower_scale_3d

  pure real(dp) function primitive_upper_scale_3d( &
      center, slope, upper) result(theta)
    real(dp), intent(in) :: center, slope, upper
    real(dp) :: magnitude

    magnitude = abs(slope)
    if (magnitude <= tiny(1.0_dp) .or. center + magnitude < upper) then
      theta = 1.0_dp
    else
      theta = max(0.0_dp, min(1.0_dp, (upper - center) / magnitude))
    end if
  end function primitive_upper_scale_3d

  pure real(dp) function primitive_slope_scale_3d( &
      center, slope, nspecies) result(theta)
    real(dp), intent(in) :: center(:), slope(:)
    integer, intent(in) :: nspecies
    integer :: species_index, component

    theta = 1.0_dp
    if (size(center) /= size(slope)) then
      theta = 0.0_dp
      return
    end if
    theta = min(theta, &
      primitive_lower_scale_3d(center(1), slope(1), density_floor))
    theta = min(theta, &
      primitive_lower_scale_3d(center(5), slope(5), pressure_floor))
    do species_index = 1, nspecies
      component = reactive_mass_fraction_component(species_index)
      theta = min(theta, &
        primitive_lower_scale_3d(center(component), slope(component), 0.0_dp))
      theta = min(theta, &
        primitive_upper_scale_3d(center(component), slope(component), 1.0_dp))
    end do
  end function primitive_slope_scale_3d

  pure subroutine sanitize_primitive_3d(q, fallback, nspecies)
    real(dp), intent(inout) :: q(:)
    real(dp), intent(in) :: fallback(:)
    integer, intent(in) :: nspecies

    real(dp) :: total
    integer :: species_index, component

    if (size(q) /= size(fallback) .or. size(q) < 5) return
    if (.not. all(ieee_is_finite(q))) then
      q = fallback
      return
    end if
    if (q(1) <= density_floor .or. q(5) <= pressure_floor) then
      q = fallback
      return
    end if
    total = 0.0_dp
    do species_index = 1, nspecies
      component = reactive_mass_fraction_component(species_index)
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
    do species_index = 1, nspecies
      component = reactive_mass_fraction_component(species_index)
      q(component) = q(component) / total
    end do
  end subroutine sanitize_primitive_3d

  subroutine primitive_face_to_state_3d( &
      species, face_primitive, fallback_primitive, state, temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: face_primitive(:), fallback_primitive(:)
    real(dp), intent(out) :: state(:), temperature
    logical, intent(out) :: ok

    real(dp) :: work(size(face_primitive))
    real(dp) :: sound_speed

    work = face_primitive
    call sanitize_primitive_3d(work, fallback_primitive, size(species))
    call reactive_primitive_to_conserved( &
      species, work, state, temperature, sound_speed, ok)
    if (ok) return
    call reactive_primitive_to_conserved( &
      species, fallback_primitive, state, temperature, sound_speed, ok)
  end subroutine primitive_face_to_state_3d

  subroutine reactive_flux_divergence_3d( &
      flux_x, flux_y, flux_z, nx, ny, nz, dx, dy, dz, rhs, ok)
    real(dp), intent(in) :: flux_x(:, :, :, :)
    real(dp), intent(in) :: flux_y(:, :, :, :)
    real(dp), intent(in) :: flux_z(:, :, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz
    real(dp), intent(out) :: rhs(:, :, :, :)
    logical, intent(out) :: ok

    integer :: i, j, k, previous_i, previous_j, previous_k

    rhs = 0.0_dp
    ok = .false.
    if (.not. all(ieee_is_finite([dx, dy, dz]))) return
    ok = nx >= 2 .and. ny >= 2 .and. nz >= 2 .and. &
      dx > 0.0_dp .and. dy > 0.0_dp .and. dz > 0.0_dp
    if (.not. ok) return
    ok = all(shape(flux_x) == shape(rhs)) .and. &
      all(shape(flux_y) == shape(rhs)) .and. &
      all(shape(flux_z) == shape(rhs)) .and. &
      size(rhs, 2) == nx .and. size(rhs, 3) == ny .and. &
      size(rhs, 4) == nz
    if (.not. ok) return
    do k = 1, nz
      previous_k = modulo(k - 2, nz) + 1
      do j = 1, ny
        previous_j = modulo(j - 2, ny) + 1
        do i = 1, nx
          previous_i = modulo(i - 2, nx) + 1
          rhs(:, i, j, k) = &
            -(flux_x(:, i, j, k) - &
              flux_x(:, previous_i, j, k)) / dx &
            -(flux_y(:, i, j, k) - &
              flux_y(:, i, previous_j, k)) / dy &
            -(flux_z(:, i, j, k) - &
              flux_z(:, i, j, previous_k)) / dz
        end do
      end do
    end do
    ok = all(ieee_is_finite(rhs))
  end subroutine reactive_flux_divergence_3d

  subroutine advance_reactive_euler_ssprk2_with_fluxes_3d( &
      species, state, temperature, nx, ny, nz, dx, dy, dz, dt, &
      riemann_solver, flux_x, flux_y, flux_z, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz, dt
    character(len=*), intent(in) :: riemann_solver
    real(dp), intent(out) :: flux_x(:, :, :, :)
    real(dp), intent(out) :: flux_y(:, :, :, :)
    real(dp), intent(out) :: flux_z(:, :, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: old_state(:, :, :, :)
    real(dp), allocatable :: stage_state(:, :, :, :)
    real(dp), allocatable :: updated_state(:, :, :, :)
    real(dp), allocatable :: rhs(:, :, :, :)
    real(dp), allocatable :: first_flux_x(:, :, :, :)
    real(dp), allocatable :: first_flux_y(:, :, :, :)
    real(dp), allocatable :: first_flux_z(:, :, :, :)
    real(dp), allocatable :: second_flux_x(:, :, :, :)
    real(dp), allocatable :: second_flux_y(:, :, :, :)
    real(dp), allocatable :: second_flux_z(:, :, :, :)
    real(dp), allocatable :: old_temperature(:, :, :)
    real(dp), allocatable :: stage_temperature(:, :, :)
    real(dp), allocatable :: updated_temperature(:, :, :)
    logical :: local_ok
    integer :: nvar

    flux_x = 0.0_dp
    flux_y = 0.0_dp
    flux_z = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return
    if (.not. ieee_is_finite(dt)) return
    if (nx < 2 .or. ny < 2 .or. nz < 2 .or. dt <= 0.0_dp) return
    if (any(shape(flux_x) /= shape(state)) .or. &
        any(shape(flux_y) /= shape(state)) .or. &
        any(shape(flux_z) /= shape(state))) return

    allocate(old_state(nvar, nx, ny, nz))
    allocate(stage_state(nvar, nx, ny, nz))
    allocate(updated_state(nvar, nx, ny, nz))
    allocate(rhs(nvar, nx, ny, nz))
    allocate(first_flux_x(nvar, nx, ny, nz))
    allocate(first_flux_y(nvar, nx, ny, nz))
    allocate(first_flux_z(nvar, nx, ny, nz))
    allocate(second_flux_x(nvar, nx, ny, nz))
    allocate(second_flux_y(nvar, nx, ny, nz))
    allocate(second_flux_z(nvar, nx, ny, nz))
    allocate(old_temperature(nx, ny, nz))
    allocate(stage_temperature(nx, ny, nz))
    allocate(updated_temperature(nx, ny, nz))
    old_state = state
    call recover_reactive_temperatures_3d( &
      species, old_state, temperature, nx, ny, nz, old_temperature, local_ok)
    if (.not. local_ok) return

    call compute_reactive_face_fluxes_3d( &
      species, old_state, old_temperature, nx, ny, nz, riemann_solver, &
      first_flux_x, first_flux_y, first_flux_z, local_ok)
    if (.not. local_ok) return
    call reactive_flux_divergence_3d( &
      first_flux_x, first_flux_y, first_flux_z, nx, ny, nz, &
      dx, dy, dz, rhs, local_ok)
    if (.not. local_ok) return
    stage_state = old_state + dt * rhs
    call recover_reactive_temperatures_3d( &
      species, stage_state, old_temperature, nx, ny, nz, &
      stage_temperature, local_ok)
    if (.not. local_ok) return

    call compute_reactive_face_fluxes_3d( &
      species, stage_state, stage_temperature, nx, ny, nz, riemann_solver, &
      second_flux_x, second_flux_y, second_flux_z, local_ok)
    if (.not. local_ok) return
    call reactive_flux_divergence_3d( &
      second_flux_x, second_flux_y, second_flux_z, nx, ny, nz, &
      dx, dy, dz, rhs, local_ok)
    if (.not. local_ok) return
    updated_state = 0.5_dp * old_state + &
      0.5_dp * (stage_state + dt * rhs)
    call recover_reactive_temperatures_3d( &
      species, updated_state, stage_temperature, nx, ny, nz, &
      updated_temperature, local_ok)
    if (.not. local_ok) return

    flux_x = 0.5_dp * (first_flux_x + second_flux_x)
    flux_y = 0.5_dp * (first_flux_y + second_flux_y)
    flux_z = 0.5_dp * (first_flux_z + second_flux_z)
    state = updated_state
    temperature = updated_temperature
    ok = .true.
  end subroutine advance_reactive_euler_ssprk2_with_fluxes_3d

  subroutine advance_reactive_euler_ssprk2_plm_3d( &
      species, state, temperature, nx, ny, nz, dx, dy, dz, dt, &
      limiter, riemann_solver, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz, dt
    character(len=*), intent(in) :: limiter, riemann_solver
    logical, intent(out) :: ok

    real(dp), allocatable :: flux_x(:, :, :, :)
    real(dp), allocatable :: flux_y(:, :, :, :)
    real(dp), allocatable :: flux_z(:, :, :, :)
    integer :: nvar

    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return
    if (nx < 3 .or. ny < 3 .or. nz < 3) return
    allocate(flux_x(nvar, nx, ny, nz))
    allocate(flux_y(nvar, nx, ny, nz))
    allocate(flux_z(nvar, nx, ny, nz))
    call advance_reactive_euler_ssprk2_plm_with_fluxes_3d( &
      species, state, temperature, nx, ny, nz, dx, dy, dz, dt, &
      limiter, riemann_solver, flux_x, flux_y, flux_z, ok)
  end subroutine advance_reactive_euler_ssprk2_plm_3d

  subroutine advance_reactive_euler_ssprk2_plm_with_fluxes_3d( &
      species, state, temperature, nx, ny, nz, dx, dy, dz, dt, &
      limiter, riemann_solver, flux_x, flux_y, flux_z, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz, dt
    character(len=*), intent(in) :: limiter, riemann_solver
    real(dp), intent(out) :: flux_x(:, :, :, :)
    real(dp), intent(out) :: flux_y(:, :, :, :)
    real(dp), intent(out) :: flux_z(:, :, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: old_state(:, :, :, :)
    real(dp), allocatable :: stage_state(:, :, :, :)
    real(dp), allocatable :: updated_state(:, :, :, :)
    real(dp), allocatable :: rhs(:, :, :, :)
    real(dp), allocatable :: first_flux_x(:, :, :, :)
    real(dp), allocatable :: first_flux_y(:, :, :, :)
    real(dp), allocatable :: first_flux_z(:, :, :, :)
    real(dp), allocatable :: second_flux_x(:, :, :, :)
    real(dp), allocatable :: second_flux_y(:, :, :, :)
    real(dp), allocatable :: second_flux_z(:, :, :, :)
    real(dp), allocatable :: old_temperature(:, :, :)
    real(dp), allocatable :: stage_temperature(:, :, :)
    real(dp), allocatable :: updated_temperature(:, :, :)
    logical :: local_ok
    integer :: nvar

    flux_x = 0.0_dp
    flux_y = 0.0_dp
    flux_z = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return
    if (.not. ieee_is_finite(dt)) return
    if (nx < 3 .or. ny < 3 .or. nz < 3 .or. dt <= 0.0_dp) return
    if (any(shape(flux_x) /= shape(state)) .or. &
        any(shape(flux_y) /= shape(state)) .or. &
        any(shape(flux_z) /= shape(state))) return

    allocate(old_state, source=state)
    allocate(stage_state(nvar, nx, ny, nz))
    allocate(updated_state(nvar, nx, ny, nz))
    allocate(rhs(nvar, nx, ny, nz))
    allocate(first_flux_x(nvar, nx, ny, nz))
    allocate(first_flux_y(nvar, nx, ny, nz))
    allocate(first_flux_z(nvar, nx, ny, nz))
    allocate(second_flux_x(nvar, nx, ny, nz))
    allocate(second_flux_y(nvar, nx, ny, nz))
    allocate(second_flux_z(nvar, nx, ny, nz))
    allocate(old_temperature(nx, ny, nz))
    allocate(stage_temperature(nx, ny, nz))
    allocate(updated_temperature(nx, ny, nz))
    call recover_reactive_temperatures_3d( &
      species, old_state, temperature, nx, ny, nz, old_temperature, local_ok)
    if (.not. local_ok) return

    call compute_reactive_plm_face_fluxes_3d( &
      species, old_state, old_temperature, nx, ny, nz, limiter, &
      riemann_solver, first_flux_x, first_flux_y, first_flux_z, local_ok)
    if (.not. local_ok) return
    call reactive_flux_divergence_3d( &
      first_flux_x, first_flux_y, first_flux_z, nx, ny, nz, &
      dx, dy, dz, rhs, local_ok)
    if (.not. local_ok) return
    stage_state = old_state + dt * rhs
    call recover_reactive_temperatures_3d( &
      species, stage_state, old_temperature, nx, ny, nz, &
      stage_temperature, local_ok)
    if (.not. local_ok) return

    call compute_reactive_plm_face_fluxes_3d( &
      species, stage_state, stage_temperature, nx, ny, nz, limiter, &
      riemann_solver, second_flux_x, second_flux_y, second_flux_z, local_ok)
    if (.not. local_ok) return
    call reactive_flux_divergence_3d( &
      second_flux_x, second_flux_y, second_flux_z, nx, ny, nz, &
      dx, dy, dz, rhs, local_ok)
    if (.not. local_ok) return
    updated_state = 0.5_dp * old_state + &
      0.5_dp * (stage_state + dt * rhs)
    call recover_reactive_temperatures_3d( &
      species, updated_state, stage_temperature, nx, ny, nz, &
      updated_temperature, local_ok)
    if (.not. local_ok) return

    flux_x = 0.5_dp * (first_flux_x + second_flux_x)
    flux_y = 0.5_dp * (first_flux_y + second_flux_y)
    flux_z = 0.5_dp * (first_flux_z + second_flux_z)
    state = updated_state
    temperature = updated_temperature
    ok = .true.
  end subroutine advance_reactive_euler_ssprk2_plm_with_fluxes_3d

  subroutine advance_reactive_chemistry_3d( &
      species, reactions, state, temperature, nx, ny, nz, interval, &
      rtol, atol, ok, active_mask, chemistry_integrator)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: interval, rtol, atol
    logical, intent(out) :: ok
    logical, intent(in), optional :: active_mask(:, :, :)
    character(len=*), intent(in), optional :: chemistry_integrator

    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    logical :: plane_ok
    integer :: plane, nvar

    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return
    if (present(active_mask)) then
      if (any(shape(active_mask) /= [nx, ny, nz])) return
    end if
    if (.not. all(ieee_is_finite([interval, rtol, atol]))) return
    if (nx < 1 .or. ny < 1 .or. nz < 1 .or. interval < 0.0_dp .or. &
        rtol <= 0.0_dp .or. atol <= 0.0_dp) return
    if (interval <= tiny(1.0_dp)) then
      ok = .true.
      return
    end if

    allocate(candidate_state, source=state)
    allocate(candidate_temperature, source=temperature)
    do plane = 1, nz
      if (present(active_mask)) then
        call advance_reactive_chemistry_2d( &
          species, reactions, candidate_state(:, :, :, plane), &
          candidate_temperature(:, :, plane), nx, ny, interval, rtol, atol, &
          plane_ok, active_mask(:, :, plane), &
          chemistry_integrator=chemistry_integrator)
      else
        call advance_reactive_chemistry_2d( &
          species, reactions, candidate_state(:, :, :, plane), &
          candidate_temperature(:, :, plane), nx, ny, interval, rtol, atol, &
          plane_ok, chemistry_integrator=chemistry_integrator)
      end if
      if (.not. plane_ok) return
    end do
    state = candidate_state
    temperature = candidate_temperature
    ok = .true.
  end subroutine advance_reactive_chemistry_3d

  subroutine advance_reactive_strang_3d( &
      species, reactions, state, temperature, nx, ny, nz, dx, dy, dz, dt, &
      riemann_solver, chemistry_enabled, rtol, atol, ok, chemistry_integrator)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz, dt, rtol, atol
    character(len=*), intent(in) :: riemann_solver
    logical, intent(in) :: chemistry_enabled
    logical, intent(out) :: ok
    character(len=*), intent(in), optional :: chemistry_integrator

    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    logical :: local_ok
    integer :: nvar

    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return
    if (.not. ieee_is_finite(dt)) return
    if (dt <= 0.0_dp) return
    if (chemistry_enabled) then
      if (.not. all(ieee_is_finite([rtol, atol]))) return
      if (rtol <= 0.0_dp .or. atol <= 0.0_dp) return
    end if

    allocate(candidate_state, source=state)
    allocate(candidate_temperature, source=temperature)
    if (chemistry_enabled) then
      call advance_reactive_chemistry_3d( &
        species, reactions, candidate_state, candidate_temperature, &
        nx, ny, nz, 0.5_dp * dt, rtol, atol, local_ok, &
        chemistry_integrator=chemistry_integrator)
      if (.not. local_ok) return
    end if
    call advance_reactive_euler_ssprk2_3d( &
      species, candidate_state, candidate_temperature, nx, ny, nz, &
      dx, dy, dz, dt, riemann_solver, local_ok)
    if (.not. local_ok) return
    if (chemistry_enabled) then
      call advance_reactive_chemistry_3d( &
        species, reactions, candidate_state, candidate_temperature, &
        nx, ny, nz, 0.5_dp * dt, rtol, atol, local_ok, &
        chemistry_integrator=chemistry_integrator)
      if (.not. local_ok) return
    end if
    state = candidate_state
    temperature = candidate_temperature
    ok = .true.
  end subroutine advance_reactive_strang_3d

  subroutine advance_reactive_full_3d( &
      species, reactions, transport, state, temperature, nx, ny, nz, &
      dx, dy, dz, dt, riemann_solver, chemistry_enabled, rtol, atol, &
      transport_enabled, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, &
      minimum_transport_theta, ok, reconstruction, limiter, chemistry_integrator)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz, dt, rtol, atol
    character(len=*), intent(in) :: riemann_solver
    logical, intent(in) :: chemistry_enabled, transport_enabled
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: minimum_transport_theta
    logical, intent(out) :: ok
    character(len=*), intent(in), optional :: reconstruction, limiter
    character(len=*), intent(in), optional :: chemistry_integrator

    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    real(dp) :: stage_transport_theta
    logical :: local_ok
    integer :: nvar
    character(len=32) :: selected_reconstruction, selected_limiter

    ok = .false.
    minimum_transport_theta = 1.0_dp
    selected_reconstruction = "pcm"
    selected_limiter = "mc"
    if (present(reconstruction)) selected_reconstruction = &
      trim(reconstruction)
    if (present(limiter)) selected_limiter = trim(limiter)
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return
    if (.not. ieee_is_finite(dt)) return
    if (dt <= 0.0_dp) return
    if (chemistry_enabled) then
      if (.not. all(ieee_is_finite([rtol, atol]))) return
      if (rtol <= 0.0_dp .or. atol <= 0.0_dp) return
    end if
    if (transport_enabled .and. size(transport) /= size(species)) return
    if (trim(selected_reconstruction) /= "pcm" .and. &
        trim(selected_reconstruction) /= "characteristic_plm") return
    if (trim(selected_limiter) /= "minmod" .and. &
        trim(selected_limiter) /= "mc") return

    allocate(candidate_state, source=state)
    allocate(candidate_temperature, source=temperature)
    if (chemistry_enabled) then
      call advance_reactive_chemistry_3d( &
        species, reactions, candidate_state, candidate_temperature, &
        nx, ny, nz, 0.5_dp * dt, rtol, atol, local_ok, &
        chemistry_integrator=chemistry_integrator)
      if (.not. local_ok) return
    end if
    if (transport_enabled) then
      call advance_reactive_transport_3d( &
        species, transport, candidate_state, candidate_temperature, &
        nx, ny, nz, dx, dy, dz, 0.5_dp * dt, viscosity_enabled, &
        thermal_conduction_enabled, species_diffusion_enabled, &
        barodiffusion_enabled, stage_transport_theta, local_ok)
      if (.not. local_ok) return
      minimum_transport_theta = min( &
        minimum_transport_theta, stage_transport_theta)
    end if
    select case (trim(selected_reconstruction))
    case ("pcm")
      call advance_reactive_euler_ssprk2_3d( &
        species, candidate_state, candidate_temperature, nx, ny, nz, &
        dx, dy, dz, dt, riemann_solver, local_ok)
    case ("characteristic_plm")
      call advance_reactive_euler_ssprk2_plm_3d( &
        species, candidate_state, candidate_temperature, nx, ny, nz, &
        dx, dy, dz, dt, selected_limiter, riemann_solver, local_ok)
    end select
    if (.not. local_ok) return
    if (transport_enabled) then
      call advance_reactive_transport_3d( &
        species, transport, candidate_state, candidate_temperature, &
        nx, ny, nz, dx, dy, dz, 0.5_dp * dt, viscosity_enabled, &
        thermal_conduction_enabled, species_diffusion_enabled, &
        barodiffusion_enabled, stage_transport_theta, local_ok)
      if (.not. local_ok) return
      minimum_transport_theta = min( &
        minimum_transport_theta, stage_transport_theta)
    end if
    if (chemistry_enabled) then
      call advance_reactive_chemistry_3d( &
        species, reactions, candidate_state, candidate_temperature, &
        nx, ny, nz, 0.5_dp * dt, rtol, atol, local_ok, &
        chemistry_integrator=chemistry_integrator)
      if (.not. local_ok) return
    end if
    state = candidate_state
    temperature = candidate_temperature
    ok = .true.
  end subroutine advance_reactive_full_3d

  subroutine recover_reactive_temperatures_3d( &
      species, state, temperature_guess, nx, ny, nz, temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature_guess(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(out) :: temperature(:, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: sound_speed
    logical :: cell_ok
    integer :: i, j, k, nvar

    temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature_guess, nvar, nx, ny, nz)) return
    if (size(temperature, 1) /= nx .or. size(temperature, 2) /= ny .or. &
        size(temperature, 3) /= nz) return

    allocate(primitive(reactive_nprim(size(species))))
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature_guess(i, j, k), &
            primitive, temperature(i, j, k), sound_speed, cell_ok)
          if (.not. cell_ok) then
            temperature = 0.0_dp
            return
          end if
        end do
      end do
    end do
    ok = all(ieee_is_finite(temperature))
    if (.not. ok) then
      temperature = 0.0_dp
      return
    end if
    ok = minval(temperature) > 0.0_dp
    if (.not. ok) temperature = 0.0_dp
  end subroutine recover_reactive_temperatures_3d

  subroutine reactive_integrals_3d( &
      state, nx, ny, nz, dx, dy, dz, integrals, ok)
    real(dp), intent(in) :: state(:, :, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz
    real(dp), intent(out) :: integrals(:)
    logical, intent(out) :: ok

    integer :: component

    integrals = 0.0_dp
    ok = .false.
    if (.not. all(ieee_is_finite([dx, dy, dz]))) return
    if (.not. all(ieee_is_finite(state))) return
    ok = size(state, 1) == size(integrals) .and. &
      size(state, 2) == nx .and. size(state, 3) == ny .and. &
      size(state, 4) == nz .and. nx > 0 .and. ny > 0 .and. nz > 0 .and. &
      dx > 0.0_dp .and. dy > 0.0_dp .and. dz > 0.0_dp
    if (.not. ok) return
    do component = 1, size(integrals)
      integrals(component) = sum(state(component, :, :, :)) * dx * dy * dz
    end do
    ok = all(ieee_is_finite(integrals))
    if (.not. ok) integrals = 0.0_dp
  end subroutine reactive_integrals_3d

  subroutine reactive_extrema_3d( &
      species, state, temperature, nx, ny, nz, &
      minimum_density, maximum_density, minimum_pressure, maximum_pressure, &
      minimum_temperature, maximum_temperature, maximum_speed, &
      maximum_closure_error, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(out) :: minimum_density, maximum_density
    real(dp), intent(out) :: minimum_pressure, maximum_pressure
    real(dp), intent(out) :: minimum_temperature, maximum_temperature
    real(dp), intent(out) :: maximum_speed, maximum_closure_error
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: local_temperature, sound_speed, closure, speed
    logical :: cell_ok
    integer :: i, j, k, species_index, nvar

    minimum_density = huge(1.0_dp)
    maximum_density = -huge(1.0_dp)
    minimum_pressure = huge(1.0_dp)
    maximum_pressure = -huge(1.0_dp)
    minimum_temperature = huge(1.0_dp)
    maximum_temperature = -huge(1.0_dp)
    maximum_speed = 0.0_dp
    maximum_closure_error = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return

    allocate(primitive(reactive_nprim(size(species))))
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), primitive, &
            local_temperature, sound_speed, cell_ok)
          if (.not. cell_ok) return
          minimum_density = min(minimum_density, primitive(1))
          maximum_density = max(maximum_density, primitive(1))
          minimum_pressure = min(minimum_pressure, primitive(5))
          maximum_pressure = max(maximum_pressure, primitive(5))
          minimum_temperature = min(minimum_temperature, local_temperature)
          maximum_temperature = max(maximum_temperature, local_temperature)
          speed = sqrt(sum(primitive(2:4)**2))
          maximum_speed = max(maximum_speed, speed)
          closure = 0.0_dp
          do species_index = 1, size(species)
            closure = closure + primitive( &
              reactive_mass_fraction_component(species_index))
          end do
          maximum_closure_error = max( &
            maximum_closure_error, abs(closure - 1.0_dp))
        end do
      end do
    end do
    ok = all(ieee_is_finite([ &
      minimum_density, maximum_density, minimum_pressure, maximum_pressure, &
      minimum_temperature, maximum_temperature, maximum_speed, &
      maximum_closure_error]))
  end subroutine reactive_extrema_3d

  pure logical function valid_reactive_3d_shapes( &
      state, temperature, nvar, nx, ny, nz) result(valid)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nvar, nx, ny, nz

    valid = nvar > 0 .and. size(state, 1) == nvar .and. &
      size(state, 2) == nx .and. size(state, 3) == ny .and. &
      size(state, 4) == nz .and. size(temperature, 1) == nx .and. &
      size(temperature, 2) == ny .and. size(temperature, 3) == nz
  end function valid_reactive_3d_shapes

end module reactive_3d_mod
