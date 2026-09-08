program test_reactive_transport_3d
  use, intrinsic :: ieee_arithmetic, only: &
    ieee_is_finite, ieee_positive_inf, ieee_quiet_nan, ieee_value
  use precision_mod, only: dp
  use state_indices_mod, only: imx, imz
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport, &
    load_h2o2_full_transport
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_species_component, reactive_conserved_to_primitive, &
    reactive_primitive_to_conserved
  use reactive_transport_2d_mod, only: advance_reactive_transport_2d
  use reactive_transport_3d_mod, only: &
    reactive_transport_face_flux_3d, reactive_transport_fluxes_3d, &
    reactive_transport_interface_theta_3d, reactive_transport_timestep_3d, &
    reactive_transport_ghosted_fluxes_3d, reactive_transport_euler_update_3d, &
    advance_reactive_transport_3d
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  logical :: ok
  integer :: direction

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "elementary thermodynamics load")
  call load_h2o2_elementary_transport(transport, ok)
  call require(ok, "elementary transport load")
  do direction = 1, 3
    call check_reduction(species, transport, direction, ok)
    call require(ok, "elementary x/y/z transport reduction")
  end do
  call check_ghosted_external_theta_contract(species, transport, ok)
  call require(ok, "elementary ghosted external theta contract")
  call check_nonfinite_contract(species, transport, ok)
  call require(ok, "elementary nonfinite transport contract")
  call check_extreme_finite_contract(species, transport, ok)
  call require(ok, "elementary extreme finite transport contract")

  call load_h2o2_full_thermo(species, ok)
  call require(ok, "full thermodynamics load")
  call load_h2o2_full_transport(transport, ok)
  call require(ok, "full transport load")
  do direction = 1, 3
    call check_reduction(species, transport, direction, ok)
    call require(ok, "full x/y/z transport reduction")
  end do
  call check_ghosted_external_theta_contract(species, transport, ok)
  call require(ok, "full ghosted external theta contract")
  call check_nonfinite_contract(species, transport, ok)
  call require(ok, "full nonfinite transport contract")
  call check_extreme_finite_contract(species, transport, ok)
  call require(ok, "full extreme finite transport contract")

  write(*, '(a)') "test_reactive_transport_3d: PASS"

contains

  subroutine check_reduction(active_species, active_transport, direction, ok_out)
    type(nasa7_species), intent(in) :: active_species(:)
    type(gas_transport_species), intent(in) :: active_transport(:)
    integer, intent(in) :: direction
    logical, intent(out) :: ok_out

    integer, parameter :: nline = 6, ntransverse = 4
    real(dp), parameter :: length = 0.012_dp, interval = 1.0e-8_dp
    real(dp), allocatable :: line_state(:, :), line_temperature(:)
    real(dp), allocatable :: reference_state(:, :, :)
    real(dp), allocatable :: reference_temperature(:, :)
    real(dp), allocatable :: state(:, :, :, :), temperature(:, :, :)
    real(dp), allocatable :: saved_state(:, :, :, :)
    real(dp), allocatable :: saved_temperature(:, :, :)
    real(dp) :: dx, dy, dz, theta_reference, theta
    real(dp) :: dt, maximum_diffusivity, expected_dt
    real(dp) :: state_error, temperature_error
    logical :: local_ok
    integer :: nx2, ny2, nx3, ny3, nz3, nvar
    integer :: i, j, k, line_index

    ok_out = .false.
    nvar = reactive_nvar(size(active_species))
    select case (direction)
    case (1)
      nx2 = nline
      ny2 = ntransverse
      nx3 = nline
      ny3 = ntransverse
      nz3 = ntransverse
      dx = length / real(nline, dp)
      dy = length / real(ntransverse, dp)
      dz = dy
    case (2)
      nx2 = ntransverse
      ny2 = nline
      nx3 = ntransverse
      ny3 = nline
      nz3 = ntransverse
      dx = length / real(ntransverse, dp)
      dy = length / real(nline, dp)
      dz = dx
    case (3)
      nx2 = nline
      ny2 = ntransverse
      nx3 = ntransverse
      ny3 = ntransverse
      nz3 = nline
      dx = length / real(ntransverse, dp)
      dy = dx
      dz = length / real(nline, dp)
    case default
      return
    end select

    allocate(line_state(nvar, nline), line_temperature(nline))
    call build_line(active_species, line_state, line_temperature, local_ok)
    if (.not. local_ok) return
    call check_face_barodiffusion( &
      active_species, active_transport, line_state(:, 1), line_temperature(1), &
      local_ok)
    if (.not. local_ok) return
    allocate(reference_state(nvar, nx2, ny2))
    allocate(reference_temperature(nx2, ny2))
    allocate(state(nvar, nx3, ny3, nz3), temperature(nx3, ny3, nz3))
    call embed_line( &
      direction, line_state, line_temperature, reference_state, &
      reference_temperature, state, temperature)

    call reactive_transport_timestep_3d( &
      active_species, active_transport, state, temperature, nx3, ny3, nz3, &
      dx, dy, dz, 0.35_dp, .true., .true., .true., &
      dt, maximum_diffusivity, local_ok)
    if (.not. local_ok) return
    expected_dt = 0.35_dp / (maximum_diffusivity * &
      (1.0_dp / dx**2 + 1.0_dp / dy**2 + 1.0_dp / dz**2))
    if (abs(dt - expected_dt) > &
        5.0e-14_dp * max(1.0_dp, expected_dt)) return

    if (direction == 3) then
      call advance_reactive_transport_2d( &
        active_species, active_transport, reference_state, &
        reference_temperature, nx2, ny2, dz, dy, interval, &
        .true., .true., .true., .true., theta_reference, local_ok)
    else
      call advance_reactive_transport_2d( &
        active_species, active_transport, reference_state, &
        reference_temperature, nx2, ny2, dx, dy, interval, &
        .true., .true., .true., .true., theta_reference, local_ok)
    end if
    if (.not. local_ok) return
    call advance_reactive_transport_3d( &
      active_species, active_transport, state, temperature, nx3, ny3, nz3, &
      dx, dy, dz, interval, .true., .true., .true., .true., theta, local_ok)
    if (.not. local_ok) return

    state_error = 0.0_dp
    temperature_error = 0.0_dp
    do k = 1, nz3
      do j = 1, ny3
        do i = 1, nx3
          select case (direction)
          case (1)
            line_index = i
            call accumulate_state_error( &
              state(:, i, j, k), reference_state(:, line_index, 1), &
              .false., state_error)
            temperature_error = max(temperature_error, abs( &
              temperature(i, j, k) - reference_temperature(line_index, 1)))
          case (2)
            line_index = j
            call accumulate_state_error( &
              state(:, i, j, k), reference_state(:, 1, line_index), &
              .false., state_error)
            temperature_error = max(temperature_error, abs( &
              temperature(i, j, k) - reference_temperature(1, line_index)))
          case (3)
            line_index = k
            call accumulate_state_error( &
              state(:, i, j, k), reference_state(:, line_index, 1), &
              .true., state_error)
            temperature_error = max(temperature_error, abs( &
              temperature(i, j, k) - reference_temperature(line_index, 1)))
          end select
        end do
      end do
    end do
    write(*, '(a,i0,a,i0,4(a,es24.16))') &
      "nspecies=", size(active_species), ", direction=", direction, &
      ", state error=", state_error, &
      ", temperature error=", temperature_error, &
      ", theta=", theta, ", dt=", dt
    if (state_error > 5.0e-13_dp .or. &
        temperature_error > 5.0e-9_dp .or. &
        abs(theta - theta_reference) > 5.0e-13_dp) return

    saved_state = state
    saved_temperature = temperature
    call advance_reactive_transport_3d( &
      active_species, active_transport(1:size(active_transport) - 1), &
      state, temperature, nx3, ny3, nz3, dx, dy, dz, interval, &
      .true., .true., .true., .true., theta, local_ok)
    ok_out = .not. local_ok .and. &
      maxval(abs(state - saved_state)) == 0.0_dp .and. &
      maxval(abs(temperature - saved_temperature)) == 0.0_dp
  end subroutine check_reduction

  subroutine check_ghosted_external_theta_contract( &
      active_species, active_transport, ok_out)
    type(nasa7_species), intent(in) :: active_species(:)
    type(gas_transport_species), intent(in) :: active_transport(:)
    logical, intent(out) :: ok_out

    integer, parameter :: nx = 2, ny = 2, nz = 2
    real(dp), parameter :: dx = 1.0e-3_dp, dy = 1.0e-3_dp
    real(dp), parameter :: dz = 1.0e-3_dp, interval = 1.0e-8_dp
    real(dp), allocatable :: line_state(:, :), line_temperature(:)
    real(dp), allocatable :: state(:, :, :, :), temperature(:, :, :)
    real(dp), allocatable :: saved_state(:, :, :, :), saved_temperature(:, :, :)
    real(dp), allocatable :: flux_x(:, :, :, :), flux_y(:, :, :, :)
    real(dp), allocatable :: flux_z(:, :, :, :), theta_cells(:, :, :)
    real(dp), allocatable :: baseline_flux_x(:, :, :, :)
    real(dp), allocatable :: baseline_flux_y(:, :, :, :)
    real(dp), allocatable :: baseline_flux_z(:, :, :, :)
    real(dp), allocatable :: baseline_theta_cells(:, :, :)
    real(dp), allocatable :: lower_theta(:, :), upper_theta(:, :)
    real(dp) :: minimum_theta, baseline_theta, bad_theta
    logical :: local_ok
    integer :: i, j, k, bad_kind, nvar
    integer :: species_first, species_last

    ok_out = .false.
    nvar = reactive_nvar(size(active_species))
    allocate(line_state(nvar, 2), line_temperature(2))
    call build_line(active_species, line_state, line_temperature, local_ok)
    if (.not. local_ok) return

    allocate(state(nvar, 0:nx + 1, 0:ny + 1, 0:nz + 1))
    allocate(temperature(0:nx + 1, 0:ny + 1, 0:nz + 1))
    do k = 0, nz + 1
      do j = 0, ny + 1
        do i = 0, nx + 1
          state(:, i, j, k) = line_state(:, 1)
          temperature(i, j, k) = line_temperature(1)
        end do
      end do
    end do
    do k = 0, nz + 1
      do j = 0, ny + 1
        state(:, 0, j, k) = line_state(:, 2)
        temperature(0, j, k) = line_temperature(2)
        state(:, 1, j, k) = line_state(:, 1)
        temperature(1, j, k) = line_temperature(1)
        state(:, 2, j, k) = line_state(:, 2)
        temperature(2, j, k) = line_temperature(2)
        state(:, 3, j, k) = line_state(:, 1)
        temperature(3, j, k) = line_temperature(1)
      end do
    end do

    allocate(saved_state, source=state)
    allocate(saved_temperature, source=temperature)
    allocate(flux_x(nvar, 0:nx, ny, nz))
    allocate(flux_y(nvar, nx, 0:ny, nz))
    allocate(flux_z(nvar, nx, ny, 0:nz))
    allocate(theta_cells(nx, ny, nz))
    allocate(baseline_flux_x, source=flux_x)
    allocate(baseline_flux_y, source=flux_y)
    allocate(baseline_flux_z, source=flux_z)
    allocate(baseline_theta_cells, source=theta_cells)
    allocate(lower_theta(ny, nz), upper_theta(ny, nz))

    lower_theta = 1.0_dp
    upper_theta = 1.0_dp
    call reactive_transport_ghosted_fluxes_3d( &
      active_species, active_transport, state, temperature, nx, ny, nz, &
      dx, dy, dz, interval, .false., .false., .true., .false., &
      flux_x, flux_y, flux_z, minimum_theta, local_ok, &
      theta_cell_output=theta_cells, &
      exterior_theta_x_lower=lower_theta, &
      exterior_theta_x_upper=upper_theta)
    if (.not. local_ok .or. &
        .not. all(ieee_is_finite(flux_x)) .or. &
        .not. all(ieee_is_finite(flux_y)) .or. &
        .not. all(ieee_is_finite(flux_z)) .or. &
        .not. all(ieee_is_finite(theta_cells)) .or. &
        .not. ieee_is_finite(minimum_theta) .or. &
        minimum_theta < 0.0_dp .or. minimum_theta > 1.0_dp) return
    species_first = reactive_species_component(1)
    species_last = reactive_species_component(size(active_species))
    if (maxval(abs(flux_x(species_first:species_last, 0, :, :))) <= 0.0_dp .or. &
        maxval(abs(flux_x(species_first:species_last, nx, :, :))) <= 0.0_dp) return
    if (maxval(abs(state - saved_state)) /= 0.0_dp .or. &
        maxval(abs(temperature - saved_temperature)) /= 0.0_dp) return
    baseline_flux_x = flux_x
    baseline_flux_y = flux_y
    baseline_flux_z = flux_z
    baseline_theta = minimum_theta
    baseline_theta_cells = theta_cells

    lower_theta = 0.0_dp
    upper_theta = 1.0_dp
    call reactive_transport_ghosted_fluxes_3d( &
      active_species, active_transport, state, temperature, nx, ny, nz, &
      dx, dy, dz, interval, .false., .false., .true., .false., &
      flux_x, flux_y, flux_z, minimum_theta, local_ok, &
      theta_cell_output=theta_cells, &
      exterior_theta_x_lower=lower_theta, &
      exterior_theta_x_upper=upper_theta)
    if (.not. local_ok .or. minimum_theta /= 0.0_dp .or. &
        maxval(abs(flux_x(species_first:species_last, 0, :, :))) /= 0.0_dp .or. &
        maxval(abs(flux_x(species_first:species_last, nx, :, :))) <= 0.0_dp .or. &
        maxval(abs(state - saved_state)) /= 0.0_dp .or. &
        maxval(abs(temperature - saved_temperature)) /= 0.0_dp) return

    lower_theta = 1.0_dp
    upper_theta = 0.0_dp
    call reactive_transport_ghosted_fluxes_3d( &
      active_species, active_transport, state, temperature, nx, ny, nz, &
      dx, dy, dz, interval, .false., .false., .true., .false., &
      flux_x, flux_y, flux_z, minimum_theta, local_ok, &
      theta_cell_output=theta_cells, &
      exterior_theta_x_lower=lower_theta, &
      exterior_theta_x_upper=upper_theta)
    if (.not. local_ok .or. minimum_theta /= 0.0_dp .or. &
        maxval(abs(flux_x(species_first:species_last, 0, :, :))) <= 0.0_dp .or. &
        maxval(abs(flux_x(species_first:species_last, nx, :, :))) /= 0.0_dp .or. &
        maxval(abs(state - saved_state)) /= 0.0_dp .or. &
        maxval(abs(temperature - saved_temperature)) /= 0.0_dp) return

    lower_theta = 1.0_dp
    upper_theta = 1.0_dp
    call reactive_transport_ghosted_fluxes_3d( &
      active_species, active_transport, state, temperature, nx, ny, nz, &
      dx, dy, dz, interval, .false., .false., .true., .false., &
      flux_x, flux_y, flux_z, minimum_theta, local_ok, &
      theta_cell_output=theta_cells, &
      exterior_theta_x_lower=lower_theta, &
      exterior_theta_x_upper=upper_theta)
    if (.not. local_ok .or. minimum_theta /= baseline_theta .or. &
        .not. all(flux_x == baseline_flux_x) .or. &
        .not. all(flux_y == baseline_flux_y) .or. &
        .not. all(flux_z == baseline_flux_z) .or. &
        .not. all(theta_cells == baseline_theta_cells) .or. &
        maxval(abs(state - saved_state)) /= 0.0_dp .or. &
        maxval(abs(temperature - saved_temperature)) /= 0.0_dp) return

    do bad_kind = 1, 3
      select case (bad_kind)
      case (1)
        bad_theta = -1.0e-6_dp
      case (2)
        bad_theta = 1.0_dp + 1.0e-6_dp
      case (3)
        bad_theta = ieee_value(0.0_dp, ieee_quiet_nan)
      end select
      lower_theta = bad_theta
      upper_theta = 1.0_dp
      call check_ghosted_invalid_external_theta( &
        active_species, active_transport, state, temperature, saved_state, &
        saved_temperature, flux_x, flux_y, flux_z, theta_cells, lower_theta, &
        upper_theta, dx, dy, dz, interval, local_ok)
      if (.not. local_ok) return
      lower_theta = 1.0_dp
      upper_theta = bad_theta
      call check_ghosted_invalid_external_theta( &
        active_species, active_transport, state, temperature, saved_state, &
        saved_temperature, flux_x, flux_y, flux_z, theta_cells, lower_theta, &
        upper_theta, dx, dy, dz, interval, local_ok)
      if (.not. local_ok) return
    end do

    write(*, '(a,i0,4(a,es12.4))') &
      "ghosted theta nspecies=", size(active_species), &
      ", accepted lower0 theta=", 0.0_dp, ", upper0 theta=", 0.0_dp, &
      ", accepted one theta=", baseline_theta, ", x-face flux=", &
      maxval(abs(baseline_flux_x(species_first:species_last, :, :, :)))
    ok_out = .true.
  end subroutine check_ghosted_external_theta_contract

  subroutine check_ghosted_invalid_external_theta( &
      active_species, active_transport, state, temperature, saved_state, &
      saved_temperature, flux_x, flux_y, flux_z, theta_cells, lower_theta, &
      upper_theta, dx, dy, dz, interval, ok_out)
    type(nasa7_species), intent(in) :: active_species(:)
    type(gas_transport_species), intent(in) :: active_transport(:)
    real(dp), intent(in) :: state(:, 0:, 0:, 0:)
    real(dp), intent(in) :: temperature(0:, 0:, 0:)
    real(dp), intent(in) :: saved_state(:, 0:, 0:, 0:)
    real(dp), intent(in) :: saved_temperature(0:, 0:, 0:)
    real(dp), intent(inout) :: flux_x(:, 0:, :, :)
    real(dp), intent(inout) :: flux_y(:, :, 0:, :)
    real(dp), intent(inout) :: flux_z(:, :, :, 0:)
    real(dp), intent(inout) :: theta_cells(:, :, :)
    real(dp), intent(in) :: lower_theta(:, :), upper_theta(:, :)
    real(dp), intent(in) :: dx, dy, dz, interval
    logical, intent(out) :: ok_out

    real(dp) :: minimum_theta
    logical :: local_ok
    integer :: nx, ny, nz

    ok_out = .false.
    nx = size(theta_cells, 1)
    ny = size(theta_cells, 2)
    nz = size(theta_cells, 3)
    flux_x = 17.0_dp
    flux_y = 19.0_dp
    flux_z = 23.0_dp
    theta_cells = 29.0_dp
    minimum_theta = 31.0_dp
    call reactive_transport_ghosted_fluxes_3d( &
      active_species, active_transport, state, temperature, nx, ny, nz, &
      dx, dy, dz, interval, .false., .false., .true., .false., &
      flux_x, flux_y, flux_z, minimum_theta, local_ok, &
      theta_cell_output=theta_cells, &
      exterior_theta_x_lower=lower_theta, &
      exterior_theta_x_upper=upper_theta)
    ok_out = .not. local_ok .and. all(state == saved_state) .and. &
      all(temperature == saved_temperature) .and. all(flux_x == 0.0_dp) .and. &
      all(flux_y == 0.0_dp) .and. all(flux_z == 0.0_dp) .and. &
      minimum_theta == 1.0_dp .and. all(theta_cells == 1.0_dp)
  end subroutine check_ghosted_invalid_external_theta

  subroutine check_nonfinite_contract( &
      active_species, active_transport, ok_out)
    type(nasa7_species), intent(in) :: active_species(:)
    type(gas_transport_species), intent(in) :: active_transport(:)
    logical, intent(out) :: ok_out

    integer, parameter :: nx = 2, ny = 2, nz = 2
    real(dp), parameter :: dx0 = 1.0e-3_dp, dy0 = 1.1e-3_dp
    real(dp), parameter :: dz0 = 1.2e-3_dp, dt0 = 1.0e-8_dp
    real(dp), allocatable :: line_state(:, :), line_temperature(:)
    real(dp), allocatable :: state(:, :, :, :), temperature(:, :, :)
    real(dp), allocatable :: bad_state(:, :, :, :), bad_temperature(:, :, :)
    real(dp), allocatable :: ghost_state(:, :, :, :)
    real(dp), allocatable :: ghost_temperature(:, :, :)
    real(dp), allocatable :: flux_x(:, :, :, :), flux_y(:, :, :, :)
    real(dp), allocatable :: flux_z(:, :, :, :), output_state(:, :, :, :)
    real(dp), allocatable :: output_temperature(:, :, :)
    real(dp), allocatable :: ghost_flux_x(:, :, :, :)
    real(dp), allocatable :: ghost_flux_y(:, :, :, :)
    real(dp), allocatable :: ghost_flux_z(:, :, :, :)
    real(dp), allocatable :: theta_cells(:, :, :)
    real(dp), allocatable :: left_primitive(:), right_primitive(:)
    real(dp), allocatable :: face_flux(:), base_state(:)
    real(dp), allocatable :: zero_state(:), zero_flux(:)
    real(dp) :: bad, bdx, bdy, bdz, bdt, checked_temperature, sound_speed
    real(dp) :: species_energy, minimum_theta, theta, dt, maximum_diffusivity
    real(dp) :: gradient(3, 3)
    logical :: local_ok
    integer :: bad_kind, bad_input, i, j, k, nvar, nprim

    ok_out = .false.
    nvar = reactive_nvar(size(active_species))
    nprim = reactive_nprim(size(active_species))
    allocate(line_state(nvar, 2), line_temperature(2))
    call build_line(active_species, line_state, line_temperature, local_ok)
    if (.not. local_ok) return

    allocate(state(nvar, nx, ny, nz), temperature(nx, ny, nz))
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          state(:, i, j, k) = line_state(:, 1)
          temperature(i, j, k) = line_temperature(1)
        end do
      end do
    end do
    allocate(bad_state, source=state)
    allocate(bad_temperature, source=temperature)
    allocate(flux_x(nvar, nx, ny, nz), flux_y(nvar, nx, ny, nz))
    allocate(flux_z(nvar, nx, ny, nz))
    allocate(output_state(nvar, nx, ny, nz))
    allocate(output_temperature(nx, ny, nz))
    allocate(left_primitive(nprim), right_primitive(nprim))
    allocate(face_flux(nvar), base_state(nvar))
    allocate(ghost_state(nvar, 0:nx + 1, 0:ny + 1, 0:nz + 1))
    allocate(ghost_temperature(0:nx + 1, 0:ny + 1, 0:nz + 1))
    allocate(ghost_flux_x(nvar, 0:nx, ny, nz))
    allocate(ghost_flux_y(nvar, nx, 0:ny, nz))
    allocate(ghost_flux_z(nvar, nx, ny, 0:nz))
    allocate(theta_cells(nx, ny, nz))
    do k = 0, nz + 1
      do j = 0, ny + 1
        do i = 0, nx + 1
          ghost_state(:, i, j, k) = line_state(:, 1)
          ghost_temperature(i, j, k) = line_temperature(1)
        end do
      end do
    end do
    base_state = line_state(:, 1)
    call reactive_conserved_to_primitive( &
      active_species, base_state, line_temperature(1), left_primitive, &
      checked_temperature, sound_speed, local_ok)
    if (.not. local_ok) return
    right_primitive = left_primitive
    gradient = 0.0_dp

    allocate(zero_state(0), zero_flux(0))
    theta = 31.0_dp
    call reactive_transport_interface_theta_3d( &
      33, zero_state, zero_flux, 1.0_dp, theta, local_ok)
    if (local_ok .or. theta /= 0.0_dp) return

    do bad_kind = 1, 2
      do k = 0, nz + 1
        do j = 0, ny + 1
          do i = 0, nx + 1
            ghost_state(:, i, j, k) = line_state(:, 1)
            ghost_temperature(i, j, k) = line_temperature(1)
          end do
        end do
      end do
      if (bad_kind == 1) then
        bad = ieee_value(0.0_dp, ieee_quiet_nan)
      else
        bad = ieee_value(0.0_dp, ieee_positive_inf)
      end if

      do bad_input = 1, 3
        face_flux = 17.0_dp
        species_energy = 19.0_dp
        if (bad_input == 1) then
          call reactive_transport_face_flux_3d( &
            active_species, active_transport, left_primitive, right_primitive, &
            bad, checked_temperature, dx0, 1, gradient, .true., .true., &
            .true., .true., face_flux, species_energy, local_ok)
        else if (bad_input == 2) then
          call reactive_transport_face_flux_3d( &
            active_species, active_transport, left_primitive, right_primitive, &
            checked_temperature, bad, dx0, 1, gradient, .true., .true., &
            .true., .true., face_flux, species_energy, local_ok)
        else
          call reactive_transport_face_flux_3d( &
            active_species, active_transport, left_primitive, right_primitive, &
            checked_temperature, checked_temperature, bad, 1, gradient, &
            .true., .true., .true., .true., face_flux, species_energy, local_ok)
        end if
        if (.not. ieee_is_finite(species_energy) .or. &
            local_ok .or. .not. all(face_flux == 0.0_dp) .or. &
            species_energy /= 0.0_dp) return
      end do

      face_flux = 0.0_dp
      theta = 31.0_dp
      call reactive_transport_interface_theta_3d( &
        size(active_species), base_state, face_flux, bad, theta, local_ok)
      if (.not. ieee_is_finite(theta) .or. local_ok .or. theta /= 0.0_dp) return
      face_flux(reactive_species_component(1)) = bad
      theta = 31.0_dp
      call reactive_transport_interface_theta_3d( &
        size(active_species), base_state, face_flux, 1.0_dp, theta, local_ok)
      if (.not. ieee_is_finite(theta) .or. local_ok .or. theta /= 0.0_dp) return
      face_flux = 0.0_dp
      bad_state = state
      bad_state(1, 1, 1, 1) = bad
      theta = 31.0_dp
      call reactive_transport_interface_theta_3d( &
        size(active_species), bad_state(:, 1, 1, 1), face_flux, 1.0_dp, &
        theta, local_ok)
      if (.not. ieee_is_finite(theta) .or. local_ok .or. theta /= 0.0_dp) return

      do bad_input = 1, 4
        bdx = dx0
        bdy = dy0
        bdz = dz0
        bdt = dt0
        select case (bad_input)
        case (1)
          bdx = bad
        case (2)
          bdy = bad
        case (3)
          bdz = bad
        case (4)
          bdt = bad
        end select
        flux_x = 17.0_dp
        flux_y = 19.0_dp
        flux_z = 23.0_dp
        minimum_theta = 31.0_dp
        call reactive_transport_fluxes_3d( &
          active_species, active_transport, state, temperature, nx, ny, nz, &
          bdx, bdy, bdz, bdt, .true., .true., .true., .true., flux_x, &
          flux_y, flux_z, minimum_theta, local_ok)
        if (.not. ieee_is_finite(minimum_theta) .or. &
            local_ok .or. .not. all(flux_x == 0.0_dp) .or. &
            .not. all(flux_y == 0.0_dp) .or. &
            .not. all(flux_z == 0.0_dp) .or. minimum_theta /= 1.0_dp) return

        dt = 17.0_dp
        maximum_diffusivity = 19.0_dp
        if (bad_input <= 3) then
          call reactive_transport_timestep_3d( &
            active_species, active_transport, state, temperature, nx, ny, nz, &
            bdx, bdy, bdz, 0.35_dp, .true., .true., .true., dt, &
            maximum_diffusivity, local_ok)
          if (.not. ieee_is_finite(dt) .or. &
              .not. ieee_is_finite(maximum_diffusivity) .or. &
              local_ok .or. dt /= 0.0_dp .or. &
              maximum_diffusivity /= 0.0_dp) return
        end if

        output_state = 17.0_dp
        output_temperature = 19.0_dp
        minimum_theta = 23.0_dp
        call reactive_transport_euler_update_3d( &
          active_species, active_transport, state, temperature, nx, ny, nz, &
          bdx, bdy, bdz, bdt, .true., .true., .true., .true., output_state, &
          output_temperature, minimum_theta, local_ok)
        if (.not. ieee_is_finite(minimum_theta) .or. &
            local_ok .or. .not. all(output_state == 0.0_dp) .or. &
            .not. all(output_temperature == 0.0_dp) .or. &
            minimum_theta /= 1.0_dp) return

        bad_state = state
        bad_temperature = temperature
        output_state = 17.0_dp
        output_temperature = 19.0_dp
        minimum_theta = 23.0_dp
        if (bad_input == 1) then
          bad_state(1, 1, 1, 1) = bad
        else if (bad_input == 2) then
          bad_temperature(1, 1, 1) = bad
        end if
        if (bad_input <= 2) then
          call reactive_transport_fluxes_3d( &
            active_species, active_transport, bad_state, bad_temperature, &
            nx, ny, nz, dx0, dy0, dz0, dt0, .true., .true., .true., .true., &
            flux_x, flux_y, flux_z, minimum_theta, local_ok)
          if (.not. ieee_is_finite(minimum_theta) .or. &
              local_ok .or. .not. all(flux_x == 0.0_dp) .or. &
              .not. all(flux_y == 0.0_dp) .or. &
              .not. all(flux_z == 0.0_dp) .or. minimum_theta /= 1.0_dp) return
          dt = 17.0_dp
          maximum_diffusivity = 19.0_dp
          call reactive_transport_timestep_3d( &
            active_species, active_transport, bad_state, bad_temperature, &
            nx, ny, nz, dx0, dy0, dz0, 0.35_dp, .true., .true., .true., &
            dt, maximum_diffusivity, local_ok)
          if (.not. ieee_is_finite(dt) .or. &
              .not. ieee_is_finite(maximum_diffusivity) .or. local_ok .or. &
              dt /= 0.0_dp .or. maximum_diffusivity /= 0.0_dp) return
        end if
      end do

      dt = 17.0_dp
      maximum_diffusivity = 19.0_dp
      call reactive_transport_timestep_3d( &
        active_species, active_transport, state, temperature, nx, ny, nz, &
        dx0, dy0, dz0, bad, .true., .true., .true., dt, maximum_diffusivity, &
        local_ok)
      if (.not. ieee_is_finite(dt) .or. &
          .not. ieee_is_finite(maximum_diffusivity) .or. local_ok .or. &
          dt /= 0.0_dp .or. maximum_diffusivity /= 0.0_dp) return

      do bad_input = 1, 4
        bdx = dx0
        bdy = dy0
        bdz = dz0
        bdt = dt0
        select case (bad_input)
        case (1)
          bdx = bad
        case (2)
          bdy = bad
        case (3)
          bdz = bad
        case (4)
          bdt = bad
        end select
        ghost_flux_x = 17.0_dp
        ghost_flux_y = 19.0_dp
        ghost_flux_z = 23.0_dp
        theta_cells = 29.0_dp
        minimum_theta = 31.0_dp
        call reactive_transport_ghosted_fluxes_3d( &
          active_species, active_transport, ghost_state, ghost_temperature, &
          nx, ny, nz, bdx, bdy, bdz, bdt, .false., .false., .true., .false., &
          ghost_flux_x, ghost_flux_y, ghost_flux_z, minimum_theta, local_ok, &
          theta_cell_output=theta_cells)
        if (.not. ieee_is_finite(minimum_theta) .or. &
            local_ok .or. .not. all(ghost_flux_x == 0.0_dp) .or. &
            .not. all(ghost_flux_y == 0.0_dp) .or. &
            .not. all(ghost_flux_z == 0.0_dp) .or. &
            .not. all(theta_cells == 1.0_dp) .or. minimum_theta /= 1.0_dp) return
      end do
      ghost_state(1, 1, 1, 1) = bad
      ghost_flux_x = 17.0_dp
      ghost_flux_y = 19.0_dp
      ghost_flux_z = 23.0_dp
      theta_cells = 29.0_dp
      minimum_theta = 31.0_dp
      call reactive_transport_ghosted_fluxes_3d( &
        active_species, active_transport, ghost_state, ghost_temperature, &
        nx, ny, nz, dx0, dy0, dz0, dt0, .false., .false., .true., .false., &
        ghost_flux_x, ghost_flux_y, ghost_flux_z, minimum_theta, local_ok, &
        theta_cell_output=theta_cells)
      if (.not. ieee_is_finite(minimum_theta) .or. &
          local_ok .or. .not. all(ghost_flux_x == 0.0_dp) .or. &
          .not. all(ghost_flux_y == 0.0_dp) .or. &
          .not. all(ghost_flux_z == 0.0_dp) .or. &
          .not. all(theta_cells == 1.0_dp) .or. minimum_theta /= 1.0_dp) return

      do bad_input = 1, 4
        bdx = dx0
        bdy = dy0
        bdz = dz0
        bdt = dt0
        select case (bad_input)
        case (1)
          bdx = bad
        case (2)
          bdy = bad
        case (3)
          bdz = bad
        case (4)
          bdt = bad
        end select
        bad_state = state
        bad_temperature = temperature
        call advance_reactive_transport_3d( &
          active_species, active_transport, bad_state, bad_temperature, nx, ny, &
          nz, bdx, bdy, bdz, bdt, .true., .true., .true., .true., &
          minimum_theta, local_ok)
        if (.not. ieee_is_finite(minimum_theta) .or. &
            local_ok .or. .not. all(bad_state == state) .or. &
            .not. all(bad_temperature == temperature) .or. &
            minimum_theta /= 1.0_dp) return
      end do
    end do

    bad_state = state
    bad_temperature = temperature
    call advance_reactive_transport_3d( &
      active_species, active_transport, bad_state, bad_temperature, nx, ny, &
      nz, 0.0_dp, dy0, dz0, 0.0_dp, .false., .false., .false., .false., &
      minimum_theta, local_ok)
    if (local_ok .or. .not. all(bad_state == state) .or. &
        .not. all(bad_temperature == temperature) .or. &
        minimum_theta /= 1.0_dp) return

    face_flux = 0.0_dp
    face_flux(reactive_species_component(1)) = huge(1.0_dp)
    theta = 31.0_dp
    call reactive_transport_interface_theta_3d( &
      size(active_species), base_state, face_flux, huge(1.0_dp), &
      theta, local_ok)
    if (local_ok .or. theta /= 0.0_dp) return

    right_primitive = left_primitive
    left_primitive(2) = huge(1.0_dp)
    right_primitive(2) = huge(1.0_dp)
    face_flux = 17.0_dp
    species_energy = 19.0_dp
    gradient = 0.0_dp
    call reactive_transport_face_flux_3d( &
      active_species, active_transport, left_primitive, right_primitive, &
      checked_temperature, checked_temperature, dx0, 1, gradient, &
      .true., .true., .true., .true., face_flux, species_energy, local_ok)
    if (local_ok .or. .not. all(face_flux == 0.0_dp) .or. &
        species_energy /= 0.0_dp) return

    left_primitive = right_primitive
    left_primitive(2) = 0.0_dp
    right_primitive(2) = 0.0_dp
    gradient = 0.0_dp
    gradient(1, 1) = huge(1.0_dp)
    gradient(2, 2) = huge(1.0_dp)
    face_flux = 17.0_dp
    species_energy = 19.0_dp
    call reactive_transport_face_flux_3d( &
      active_species, active_transport, left_primitive, right_primitive, &
      checked_temperature, checked_temperature, dx0, 1, gradient, &
      .true., .false., .false., .false., face_flux, species_energy, local_ok)
    if (local_ok .or. .not. all(face_flux == 0.0_dp) .or. &
        species_energy /= 0.0_dp) return

    dt = 17.0_dp
    maximum_diffusivity = 19.0_dp
    call reactive_transport_timestep_3d( &
      active_species, active_transport, state, temperature, nx, ny, nz, &
      tiny(1.0_dp), dy0, dz0, 0.35_dp, .true., .true., .true., dt, &
      maximum_diffusivity, local_ok)
    if (local_ok .or. dt /= 0.0_dp .or. maximum_diffusivity /= 0.0_dp) return

    dt = 17.0_dp
    maximum_diffusivity = 19.0_dp
    call reactive_transport_timestep_3d( &
      active_species, active_transport, state, temperature, nx, ny, nz, &
      huge(1.0_dp), dy0, dz0, 0.35_dp, .true., .true., .true., dt, &
      maximum_diffusivity, local_ok)
    if (local_ok .or. dt /= 0.0_dp .or. maximum_diffusivity /= 0.0_dp) return

    ok_out = .true.
  end subroutine check_nonfinite_contract

  subroutine check_extreme_finite_contract( &
      active_species, active_transport, ok_out)
    type(nasa7_species), intent(in) :: active_species(:)
    type(gas_transport_species), intent(in) :: active_transport(:)
    logical, intent(out) :: ok_out

    integer, parameter :: nx = 3, ny = 3, nz = 3
    real(dp), parameter :: dx0 = 1.0e-3_dp, dy0 = 1.1e-3_dp
    real(dp), parameter :: dz0 = 1.2e-3_dp, dt0 = 1.0e-8_dp
    real(dp), allocatable :: line_state(:, :), line_temperature(:)
    real(dp), allocatable :: state(:, :, :, :), temperature(:, :, :)
    real(dp), allocatable :: saved_state(:, :, :, :), saved_temperature(:, :, :)
    real(dp), allocatable :: ghost_state(:, :, :, :), ghost_temperature(:, :, :)
    real(dp), allocatable :: flux_x(:, :, :, :), flux_y(:, :, :, :)
    real(dp), allocatable :: flux_z(:, :, :, :)
    real(dp), allocatable :: ghost_flux_x(:, :, :, :), ghost_flux_y(:, :, :, :)
    real(dp), allocatable :: ghost_flux_z(:, :, :, :)
    real(dp), allocatable :: output_state(:, :, :, :)
    real(dp), allocatable :: output_temperature(:, :, :)
    real(dp) :: extreme_spacing, minimum_theta
    real(dp) :: extreme_dt
    logical :: local_ok
    integer :: nvar, i, j, k, extreme_case

    ok_out = .false.
    nvar = reactive_nvar(size(active_species))
    allocate(line_state(nvar, nx), line_temperature(nx))
    call build_line(active_species, line_state, line_temperature, local_ok)
    if (.not. local_ok) return

    allocate(state(nvar, nx, ny, nz), temperature(nx, ny, nz))
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          state(:, i, j, k) = line_state(:, i)
          temperature(i, j, k) = line_temperature(i)
        end do
      end do
    end do
    allocate(saved_state, source=state)
    allocate(saved_temperature, source=temperature)
    allocate(flux_x(nvar, nx, ny, nz), flux_y(nvar, nx, ny, nz))
    allocate(flux_z(nvar, nx, ny, nz))
    allocate(output_state(nvar, nx, ny, nz))
    allocate(output_temperature(nx, ny, nz))

    do extreme_case = 1, 2
      if (extreme_case == 1) then
        extreme_spacing = 0.01_dp * tiny(1.0_dp)
      else
        extreme_spacing = huge(1.0_dp)
      end if
      flux_x = 17.0_dp
      flux_y = 19.0_dp
      flux_z = 23.0_dp
      minimum_theta = 31.0_dp
      call reactive_transport_fluxes_3d( &
        active_species, active_transport, state, temperature, nx, ny, nz, &
        extreme_spacing, dy0, dz0, dt0, .true., .false., .false., .false., &
        flux_x, flux_y, flux_z, minimum_theta, local_ok)
      if (local_ok .or. .not. ieee_is_finite(minimum_theta) .or. &
          .not. all(flux_x == 0.0_dp) .or. &
          .not. all(flux_y == 0.0_dp) .or. &
          .not. all(flux_z == 0.0_dp) .or. minimum_theta /= 1.0_dp) return
    end do

    allocate(ghost_state(nvar, 0:nx + 1, 0:ny + 1, 0:nz + 1))
    allocate(ghost_temperature(0:nx + 1, 0:ny + 1, 0:nz + 1))
    do k = 0, nz + 1
      do j = 0, ny + 1
        do i = 0, nx + 1
          ghost_state(:, i, j, k) = line_state(:, 1)
          ghost_temperature(i, j, k) = line_temperature(1)
        end do
      end do
    end do
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ghost_state(:, i, j, k) = line_state(:, i)
          ghost_temperature(i, j, k) = line_temperature(i)
        end do
      end do
    end do
    allocate(ghost_flux_x(nvar, 0:nx, ny, nz))
    allocate(ghost_flux_y(nvar, nx, 0:ny, nz))
    allocate(ghost_flux_z(nvar, nx, ny, 0:nz))

    do extreme_case = 1, 2
      if (extreme_case == 1) then
        extreme_spacing = 0.01_dp * tiny(1.0_dp)
      else
        extreme_spacing = huge(1.0_dp)
      end if
      ghost_flux_x = 17.0_dp
      ghost_flux_y = 19.0_dp
      ghost_flux_z = 23.0_dp
      minimum_theta = 31.0_dp
      call reactive_transport_ghosted_fluxes_3d( &
        active_species, active_transport, ghost_state, ghost_temperature, &
        nx, ny, nz, extreme_spacing, dy0, dz0, dt0, .true., .false., &
        .false., .false., ghost_flux_x, ghost_flux_y, ghost_flux_z, &
        minimum_theta, local_ok)
      if (local_ok .or. .not. ieee_is_finite(minimum_theta) .or. &
          .not. all(ghost_flux_x == 0.0_dp) .or. &
          .not. all(ghost_flux_y == 0.0_dp) .or. &
          .not. all(ghost_flux_z == 0.0_dp) .or. &
          minimum_theta /= 1.0_dp) return
    end do

    extreme_dt = huge(1.0_dp)
    flux_x = 17.0_dp
    flux_y = 19.0_dp
    flux_z = 23.0_dp
    minimum_theta = 31.0_dp
    call reactive_transport_fluxes_3d( &
      active_species, active_transport, state, temperature, nx, ny, nz, &
      dx0, dy0, dz0, extreme_dt, .false., .false., .true., .false., &
      flux_x, flux_y, flux_z, minimum_theta, local_ok)
    if (local_ok .or. .not. ieee_is_finite(minimum_theta) .or. &
        .not. all(flux_x == 0.0_dp) .or. &
        .not. all(flux_y == 0.0_dp) .or. &
        .not. all(flux_z == 0.0_dp) .or. minimum_theta /= 1.0_dp) return

    ghost_flux_x = 17.0_dp
    ghost_flux_y = 19.0_dp
    ghost_flux_z = 23.0_dp
    minimum_theta = 31.0_dp
    call reactive_transport_ghosted_fluxes_3d( &
      active_species, active_transport, ghost_state, ghost_temperature, &
      nx, ny, nz, dx0, dy0, dz0, extreme_dt, .false., .false., .true., &
      .false., ghost_flux_x, ghost_flux_y, ghost_flux_z, minimum_theta, local_ok)
    if (local_ok .or. .not. ieee_is_finite(minimum_theta) .or. &
        .not. all(ghost_flux_x == 0.0_dp) .or. &
        .not. all(ghost_flux_y == 0.0_dp) .or. &
        .not. all(ghost_flux_z == 0.0_dp) .or. minimum_theta /= 1.0_dp) return

    output_state = 17.0_dp
    output_temperature = 19.0_dp
    minimum_theta = 23.0_dp
    call reactive_transport_euler_update_3d( &
      active_species, active_transport, state, temperature, nx, ny, nz, &
      dx0, dy0, dz0, extreme_dt, .true., .true., .false., .false., &
      output_state, output_temperature, minimum_theta, local_ok)
    if (local_ok .or. .not. ieee_is_finite(minimum_theta) .or. &
        .not. all(output_state == 0.0_dp) .or. &
        .not. all(output_temperature == 0.0_dp) .or. minimum_theta /= 1.0_dp) return

    state = saved_state
    temperature = saved_temperature
    call advance_reactive_transport_3d( &
      active_species, active_transport, state, temperature, nx, ny, nz, &
      dx0, dy0, dz0, extreme_dt, .true., .true., .false., .false., &
      minimum_theta, local_ok)
    if (local_ok .or. .not. ieee_is_finite(minimum_theta) .or. &
        .not. all(state == saved_state) .or. &
        .not. all(temperature == saved_temperature) .or. &
        minimum_theta /= 1.0_dp) return

    ok_out = .true.
  end subroutine check_extreme_finite_contract

  subroutine check_face_barodiffusion( &
      active_species, active_transport, conserved_state, state_temperature, &
      ok_out)
    type(nasa7_species), intent(in) :: active_species(:)
    type(gas_transport_species), intent(in) :: active_transport(:)
    real(dp), intent(in) :: conserved_state(:), state_temperature
    logical, intent(out) :: ok_out

    real(dp), allocatable :: left_primitive(:), right_primitive(:), flux(:)
    real(dp) :: checked_temperature, sound_speed, species_energy
    real(dp) :: gradient(3, 3)
    logical :: local_ok

    ok_out = .false.
    allocate(left_primitive(reactive_nprim(size(active_species))))
    allocate(right_primitive(reactive_nprim(size(active_species))))
    allocate(flux(reactive_nvar(size(active_species))))
    call reactive_conserved_to_primitive( &
      active_species, conserved_state, state_temperature, left_primitive, &
      checked_temperature, sound_speed, local_ok)
    if (.not. local_ok) return
    right_primitive = left_primitive
    gradient = 0.0_dp
    call reactive_transport_face_flux_3d( &
      active_species, active_transport, left_primitive, right_primitive, &
      checked_temperature, checked_temperature, 1.0e-3_dp, 1, gradient, &
      .true., .false., .false., .true., flux, species_energy, local_ok)
    ok_out = .not. local_ok
  end subroutine check_face_barodiffusion

  subroutine embed_line( &
      direction, line_state, line_temperature, state_2d, temperature_2d, &
      state_3d, temperature_3d)
    integer, intent(in) :: direction
    real(dp), intent(in) :: line_state(:, :), line_temperature(:)
    real(dp), intent(out) :: state_2d(:, :, :), temperature_2d(:, :)
    real(dp), intent(out) :: state_3d(:, :, :, :), temperature_3d(:, :, :)

    integer :: i, j, k, line_index

    do j = 1, size(state_2d, 3)
      do i = 1, size(state_2d, 2)
        if (direction == 2) then
          line_index = j
        else
          line_index = i
        end if
        state_2d(:, i, j) = line_state(:, line_index)
        temperature_2d(i, j) = line_temperature(line_index)
      end do
    end do
    do k = 1, size(state_3d, 4)
      do j = 1, size(state_3d, 3)
        do i = 1, size(state_3d, 2)
          select case (direction)
          case (1)
            line_index = i
          case (2)
            line_index = j
          case (3)
            line_index = k
          end select
          state_3d(:, i, j, k) = line_state(:, line_index)
          if (direction == 3) then
            state_3d(imx, i, j, k) = line_state(imz, line_index)
            state_3d(imz, i, j, k) = line_state(imx, line_index)
          end if
          temperature_3d(i, j, k) = line_temperature(line_index)
        end do
      end do
    end do
  end subroutine embed_line

  subroutine build_line(active_species, state, temperature, ok_out)
    type(nasa7_species), intent(in) :: active_species(:)
    real(dp), intent(out) :: state(:, :), temperature(:)
    logical, intent(out) :: ok_out

    real(dp), allocatable :: mole_fractions(:), mass_fractions(:)
    real(dp), allocatable :: primitive(:)
    real(dp) :: phase, pressure, density, sound_speed
    logical :: local_ok
    integer :: cell, species_index, nspecies

    ok_out = .false.
    nspecies = size(active_species)
    allocate(mole_fractions(nspecies), mass_fractions(nspecies))
    allocate(primitive(reactive_nprim(nspecies)))
    do cell = 1, size(state, 2)
      call base_mole_fractions(nspecies, mole_fractions, local_ok)
      if (.not. local_ok) return
      phase = 2.0_dp * acos(-1.0_dp) * &
        (real(cell, dp) - 0.5_dp) / real(size(state, 2), dp)
      mole_fractions(1) = mole_fractions(1) + 0.005_dp * sin(phase)
      mole_fractions(nspecies) = &
        mole_fractions(nspecies) - 0.005_dp * sin(phase)
      call mass_fractions_from_mole_fractions( &
        active_species, mole_fractions, mass_fractions, local_ok)
      if (.not. local_ok) return
      temperature(cell) = 1000.0_dp + 10.0_dp * cos(phase)
      pressure = 101325.0_dp * (1.0_dp + 0.01_dp * sin(phase))
      density = mixture_density( &
        active_species, mass_fractions, pressure, temperature(cell), local_ok)
      if (.not. local_ok) return
      primitive(1:5) = [ &
        density, 0.20_dp * sin(phase), -0.15_dp * cos(phase), &
        0.10_dp * sin(2.0_dp * phase), pressure]
      do species_index = 1, nspecies
        primitive(reactive_mass_fraction_component(species_index)) = &
          mass_fractions(species_index)
      end do
      call reactive_primitive_to_conserved( &
        active_species, primitive, state(:, cell), temperature(cell), &
        sound_speed, local_ok)
      if (.not. local_ok) return
    end do
    ok_out = .true.
  end subroutine build_line

  subroutine base_mole_fractions(nspecies, mole_fractions, ok_out)
    integer, intent(in) :: nspecies
    real(dp), intent(out) :: mole_fractions(:)
    logical, intent(out) :: ok_out

    ok_out = size(mole_fractions) == nspecies
    if (.not. ok_out) return
    select case (nspecies)
    case (7)
      mole_fractions = [ &
        0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
        1.0e-5_dp, 0.0_dp, 0.55643_dp]
    case (10)
      mole_fractions = [ &
        0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
        1.0e-5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.55643_dp]
    case default
      ok_out = .false.
    end select
  end subroutine base_mole_fractions

  subroutine accumulate_state_error(actual, reference, rotate_z, error)
    real(dp), intent(in) :: actual(:), reference(:)
    logical, intent(in) :: rotate_z
    real(dp), intent(inout) :: error

    real(dp), allocatable :: canonical(:)

    allocate(canonical, source=actual)
    if (rotate_z) then
      canonical(imx) = actual(imz)
      canonical(imz) = actual(imx)
    end if
    error = max(error, maxval( &
      abs(canonical - reference) / max(1.0_dp, abs(reference))))
  end subroutine accumulate_state_error

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') "FAILED: " // trim(message)
      error stop 1
    end if
  end subroutine require

end program test_reactive_transport_3d
