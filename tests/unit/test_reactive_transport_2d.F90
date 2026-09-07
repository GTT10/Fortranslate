program test_reactive_transport_2d
  use, intrinsic :: ieee_arithmetic, only: &
    ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imx, imy
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_elementary_transport
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved, advance_reactive_transport
  use reactive_transport_2d_mod, only: &
    reactive_transport_fluxes_2d, reactive_transport_fluxes_2d_faces, &
    reactive_transport_timestep_2d, reactive_transport_euler_update_2d, &
    advance_reactive_transport_2d, face_transport_data, species_face_flux
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  logical :: ok

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamics load")
  call load_h2o2_elementary_transport(transport, ok)
  call require(ok, "transport database load")
  call test_uniform_zero_flux()
  call test_x_dimensional_reduction()
  call test_y_dimensional_reduction()
  call test_trace_species_limiter()
  call test_face_transport_fail_closed()
  call test_species_face_flux_fail_closed()
  call test_public_validation_and_rollback()

contains

  subroutine test_uniform_zero_flux()
    integer, parameter :: nx = 6, ny = 5
    real(dp), allocatable :: state(:, :, :), temperature(:, :)
    real(dp), allocatable :: flux_x(:, :, :), flux_y(:, :, :)
    real(dp), allocatable :: primitive(:), mass_fractions(:)
    real(dp) :: mole_fractions(7), density, sound_speed, theta
    logical :: local_ok
    integer :: i, j, k, nvar

    nvar = reactive_nvar(size(species))
    allocate(state(nvar, nx, ny), temperature(nx, ny))
    allocate(flux_x(nvar, nx, ny), flux_y(nvar, nx, ny))
    allocate(primitive(reactive_nprim(size(species))))
    allocate(mass_fractions(size(species)))
    mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
      1.0e-5_dp, 0.0_dp, 0.55643_dp]
    call mass_fractions_from_mole_fractions( &
      species, mole_fractions, mass_fractions, local_ok)
    call require(local_ok, "uniform composition")
    density = mixture_density( &
      species, mass_fractions, 101325.0_dp, 1000.0_dp, local_ok)
    call require(local_ok, "uniform density")
    primitive(1:5) = [density, 12.0_dp, -3.0_dp, 2.0_dp, 101325.0_dp]
    do k = 1, size(species)
      primitive(reactive_mass_fraction_component(k)) = mass_fractions(k)
    end do
    do j = 1, ny
      do i = 1, nx
        call reactive_primitive_to_conserved( &
          species, primitive, state(:, i, j), temperature(i, j), &
          sound_speed, local_ok)
        call require(local_ok, "uniform state construction")
      end do
    end do
    call reactive_transport_fluxes_2d( &
      species, transport, state, temperature, nx, ny, 1.0e-3_dp, &
      1.2e-3_dp, 1.0e-7_dp, .true., .true., .true., .true., &
      flux_x, flux_y, theta, local_ok)
    call require(local_ok, "uniform transport flux")
    call require(maxval(abs(flux_x)) < 1.0e-12_dp, &
      "uniform x transport flux is zero")
    call require(maxval(abs(flux_y)) < 1.0e-12_dp, &
      "uniform y transport flux is zero")
    call require(theta > 0.999999999999_dp, &
      "uniform transport does not activate positivity limiter")
  end subroutine test_uniform_zero_flux

  subroutine test_x_dimensional_reduction()
    integer, parameter :: nx = 18, ny = 4
    real(dp), allocatable :: state_1d(:, :), temperature_1d(:)
    real(dp), allocatable :: state_2d(:, :, :), temperature_2d(:, :)
    real(dp), allocatable :: primitive(:), mass_fractions(:)
    real(dp) :: mole_fractions(7), density, sound_speed, x, dx, dy, dt
    real(dp) :: theta, difference, scale
    logical :: local_ok
    integer :: i, j, k, nvar

    nvar = reactive_nvar(size(species))
    allocate(state_1d(nvar, 0:nx + 1), temperature_1d(0:nx + 1))
    allocate(state_2d(nvar, nx, ny), temperature_2d(nx, ny))
    allocate(primitive(reactive_nprim(size(species))))
    allocate(mass_fractions(size(species)))
    mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
      1.0e-5_dp, 0.0_dp, 0.55643_dp]
    call mass_fractions_from_mole_fractions( &
      species, mole_fractions, mass_fractions, local_ok)
    call require(local_ok, "reduction composition")
    density = mixture_density( &
      species, mass_fractions, 101325.0_dp, 1000.0_dp, local_ok)
    call require(local_ok, "reduction density")
    dx = 0.01_dp / real(nx, dp)
    dy = 0.002_dp / real(ny, dp)
    do i = 1, nx
      x = (real(i, dp) - 0.5_dp) * dx
      primitive(1:5) = [density, 0.0_dp, &
        0.02_dp * sin(2.0_dp * acos(-1.0_dp) * x / 0.01_dp), &
        0.0_dp, 101325.0_dp]
      do k = 1, size(species)
        primitive(reactive_mass_fraction_component(k)) = mass_fractions(k)
      end do
      call reactive_primitive_to_conserved( &
        species, primitive, state_1d(:, i), temperature_1d(i), sound_speed, &
        local_ok)
      call require(local_ok, "reduction state")
      do j = 1, ny
        state_2d(:, i, j) = state_1d(:, i)
        temperature_2d(i, j) = temperature_1d(i)
      end do
    end do
    state_1d(:, 0) = state_1d(:, nx)
    state_1d(:, nx + 1) = state_1d(:, 1)
    temperature_1d(0) = temperature_1d(nx)
    temperature_1d(nx + 1) = temperature_1d(1)
    dt = 2.0e-6_dp
    call advance_reactive_transport( &
      species, transport, state_1d, temperature_1d, nx, dx, dt, &
      "periodic", .true., .true., .true., .true., local_ok)
    call require(local_ok, "1D reduction reference")
    call advance_reactive_transport_2d( &
      species, transport, state_2d, temperature_2d, nx, ny, dx, dy, dt, &
      .true., .true., .true., .true., theta, local_ok)
    call require(local_ok, "2D reduction update")
    call require(theta > 0.999999999_dp, &
      "smooth reduction does not activate limiter")
    difference = 0.0_dp
    scale = 1.0_dp
    do j = 1, ny
      do i = 1, nx
        difference = max(difference, &
          maxval(abs(state_2d(:, i, j) - state_1d(:, i))))
        scale = max(scale, maxval(abs(state_1d(:, i))))
        difference = max(difference, &
          abs(temperature_2d(i, j) - temperature_1d(i)))
        scale = max(scale, abs(temperature_1d(i)))
      end do
    end do
    call require(difference / scale < 5.0e-13_dp, &
      "2D molecular transport reduces to the 1D operator")
    call require(abs(sum(state_2d(imx, :, :))) < 1.0e-12_dp, &
      "reduction keeps zero normal momentum")
    call require(ieee_safe(sum(state_2d(imy, :, :))), &
      "reduction transverse momentum remains finite")
  end subroutine test_x_dimensional_reduction



  subroutine test_y_dimensional_reduction()
    integer, parameter :: nx = 4, ny = 18
    real(dp), allocatable :: state_1d(:, :), temperature_1d(:)
    real(dp), allocatable :: state_2d(:, :, :), temperature_2d(:, :)
    real(dp), allocatable :: primitive(:), mass_fractions(:), rotated(:)
    real(dp) :: mole_fractions(7), density, sound_speed, y, dx, dy, dt
    real(dp) :: theta, difference, scale
    logical :: local_ok
    integer :: i, j, k, nvar

    nvar = reactive_nvar(size(species))
    allocate(state_1d(nvar, 0:ny + 1), temperature_1d(0:ny + 1))
    allocate(state_2d(nvar, nx, ny), temperature_2d(nx, ny))
    allocate(primitive(reactive_nprim(size(species))))
    allocate(mass_fractions(size(species)), rotated(nvar))
    mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
      1.0e-5_dp, 0.0_dp, 0.55643_dp]
    call mass_fractions_from_mole_fractions( &
      species, mole_fractions, mass_fractions, local_ok)
    call require(local_ok, "y reduction composition")
    density = mixture_density( &
      species, mass_fractions, 101325.0_dp, 1000.0_dp, local_ok)
    call require(local_ok, "y reduction density")
    dx = 0.002_dp / real(nx, dp)
    dy = 0.01_dp / real(ny, dp)
    do j = 1, ny
      y = (real(j, dp) - 0.5_dp) * dy
      ! 1D x velocity is the physical y-normal velocity.  Its first
      ! transverse velocity becomes the physical x velocity after rotation.
      primitive(1:5) = [density, 0.0_dp, &
        0.02_dp * sin(2.0_dp * acos(-1.0_dp) * y / 0.01_dp), &
        0.0_dp, 101325.0_dp]
      do k = 1, size(species)
        primitive(reactive_mass_fraction_component(k)) = mass_fractions(k)
      end do
      call reactive_primitive_to_conserved( &
        species, primitive, state_1d(:, j), temperature_1d(j), sound_speed, &
        local_ok)
      call require(local_ok, "y reduction state")
      rotated = state_1d(:, j)
      rotated(imx) = state_1d(imy, j)
      rotated(imy) = state_1d(imx, j)
      do i = 1, nx
        state_2d(:, i, j) = rotated
        temperature_2d(i, j) = temperature_1d(j)
      end do
    end do
    state_1d(:, 0) = state_1d(:, ny)
    state_1d(:, ny + 1) = state_1d(:, 1)
    temperature_1d(0) = temperature_1d(ny)
    temperature_1d(ny + 1) = temperature_1d(1)
    dt = 2.0e-6_dp
    call advance_reactive_transport( &
      species, transport, state_1d, temperature_1d, ny, dy, dt, &
      "periodic", .true., .true., .true., .true., local_ok)
    call require(local_ok, "1D y reduction reference")
    call advance_reactive_transport_2d( &
      species, transport, state_2d, temperature_2d, nx, ny, dx, dy, dt, &
      .true., .true., .true., .true., theta, local_ok)
    call require(local_ok, "2D y reduction update")
    call require(theta > 0.999999999_dp, &
      "smooth y reduction does not activate limiter")
    difference = 0.0_dp
    scale = 1.0_dp
    do j = 1, ny
      rotated = state_1d(:, j)
      rotated(imx) = state_1d(imy, j)
      rotated(imy) = state_1d(imx, j)
      do i = 1, nx
        difference = max(difference, &
          maxval(abs(state_2d(:, i, j) - rotated)))
        scale = max(scale, maxval(abs(rotated)))
        difference = max(difference, &
          abs(temperature_2d(i, j) - temperature_1d(j)))
        scale = max(scale, abs(temperature_1d(j)))
      end do
    end do
    call require(difference / scale < 5.0e-13_dp, &
      "rotated 2D molecular transport reduces to the 1D operator")
  end subroutine test_y_dimensional_reduction

  subroutine test_trace_species_limiter()
    integer, parameter :: nx = 6, ny = 4
    real(dp), allocatable :: state(:, :, :), temperature(:, :)
    real(dp), allocatable :: flux_x(:, :, :), flux_y(:, :, :)
    real(dp), allocatable :: primitive(:), mass_fractions(:)
    real(dp) :: mole_fractions(7), density, sound_speed, theta
    real(dp) :: closure_scale
    logical :: local_ok
    integer :: i, j, k, nvar

    nvar = reactive_nvar(size(species))
    allocate(state(nvar, nx, ny), temperature(nx, ny))
    allocate(flux_x(nvar, nx, ny), flux_y(nvar, nx, ny))
    allocate(primitive(reactive_nprim(size(species))))
    allocate(mass_fractions(size(species)))
    do j = 1, ny
      do i = 1, nx
        mole_fractions = [0.29570_dp, 1.0e-12_dp, 1.0e-5_dp, &
          0.14784_dp, 1.0e-5_dp, 0.0_dp, 0.556439999999_dp]
        if (mod(i, 2) == 0) then
          mole_fractions(2) = 0.02_dp
          mole_fractions(7) = mole_fractions(7) - 0.02_dp
        end if
        mole_fractions = mole_fractions / sum(mole_fractions)
        call mass_fractions_from_mole_fractions( &
          species, mole_fractions, mass_fractions, local_ok)
        call require(local_ok, "trace limiter composition")
        density = mixture_density( &
          species, mass_fractions, 101325.0_dp, 1000.0_dp, local_ok)
        call require(local_ok, "trace limiter density")
        primitive(1:5) = [density, 0.0_dp, 0.0_dp, 0.0_dp, 101325.0_dp]
        do k = 1, size(species)
          primitive(reactive_mass_fraction_component(k)) = mass_fractions(k)
        end do
        call reactive_primitive_to_conserved( &
          species, primitive, state(:, i, j), temperature(i, j), &
          sound_speed, local_ok)
        call require(local_ok, "trace limiter state")
      end do
    end do
    call reactive_transport_fluxes_2d( &
      species, transport, state, temperature, nx, ny, 1.0e-3_dp, &
      1.0e-3_dp, 2.0e-2_dp, .false., .false., .true., .false., &
      flux_x, flux_y, theta, local_ok)
    call require(local_ok, "trace limiter flux construction")
    call require(theta >= 0.0_dp .and. theta < 0.999_dp, &
      "trace-species limiter activates for an oversized explicit interval")
    closure_scale = max(1.0_dp, maxval(abs(flux_x)), maxval(abs(flux_y)))
    do j = 1, ny
      do i = 1, nx
        call require(abs(sum(flux_x(6:, i, j))) < &
          5.0e-12_dp * closure_scale, "limited x species flux closes")
        call require(abs(sum(flux_y(6:, i, j))) < &
          5.0e-12_dp * closure_scale, "limited y species flux closes")
      end do
    end do
  end subroutine test_trace_species_limiter

  subroutine test_face_transport_fail_closed()
    real(dp), allocatable :: left(:), right(:), mass_fractions(:)
    real(dp), allocatable :: diffusion(:), yleft(:), yright(:), yface(:)
    real(dp), allocatable :: xleft(:), xright(:), hface(:)
    type(gas_transport_species), allocatable :: bad_transport(:)
    real(dp) :: mole_fractions(7), density, sound_speed
    real(dp) :: viscosity, conductivity, left_temperature, right_temperature
    logical :: local_ok
    integer :: k, nspecies, nprim

    nspecies = size(species)
    nprim = reactive_nprim(nspecies)
    allocate(left(nprim), right(nprim), mass_fractions(nspecies))
    allocate(diffusion(nspecies), yleft(nspecies), yright(nspecies))
    allocate(yface(nspecies), xleft(nspecies), xright(nspecies))
    allocate(hface(nspecies), bad_transport(nspecies))
    mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
      1.0e-5_dp, 0.0_dp, 0.55643_dp]
    call mass_fractions_from_mole_fractions( &
      species, mole_fractions, mass_fractions, local_ok)
    call require(local_ok, "face failure reference composition")
    density = mixture_density( &
      species, mass_fractions, 101325.0_dp, 1000.0_dp, local_ok)
    call require(local_ok, "face failure reference density")
    left = 0.0_dp
    right = 0.0_dp
    left(1) = density
    right(1) = density
    left(5) = 101325.0_dp
    right(5) = 101325.0_dp
    do k = 1, nspecies
      left(reactive_mass_fraction_component(k)) = mass_fractions(k)
      right(reactive_mass_fraction_component(k)) = mass_fractions(k)
    end do
    left_temperature = 1000.0_dp
    right_temperature = 1000.0_dp
    call face_transport_data( &
      species, transport, left, right, left_temperature, right_temperature, &
      viscosity, conductivity, diffusion, yleft, yright, yface, xleft, &
      xright, hface, local_ok)
    call require(local_ok, "face transport reference data")
    call require(ieee_is_finite(viscosity) .and. &
      ieee_is_finite(conductivity) .and. all(ieee_is_finite(diffusion)) .and. &
      all(ieee_is_finite(hface)), "face transport reference data finite")

    left_temperature = ieee_value(0.0_dp, ieee_quiet_nan)
    call face_transport_data( &
      species, transport, left, right, left_temperature, right_temperature, &
      viscosity, conductivity, diffusion, yleft, yright, yface, xleft, &
      xright, hface, local_ok)
    call require(.not. local_ok, "face transport rejects NaN temperature")
    call require(viscosity == 0.0_dp .and. conductivity == 0.0_dp .and. &
      all(diffusion == 0.0_dp) .and. all(yleft == 0.0_dp) .and. &
      all(yright == 0.0_dp) .and. all(yface == 0.0_dp) .and. &
      all(xleft == 0.0_dp) .and. all(xright == 0.0_dp) .and. &
      all(hface == 0.0_dp), &
      "NaN face failure rolls back every output")

    left_temperature = 1000.0_dp
    right_temperature = ieee_value(0.0_dp, ieee_positive_inf)
    call face_transport_data( &
      species, transport, left, right, left_temperature, right_temperature, &
      viscosity, conductivity, diffusion, yleft, yright, yface, xleft, &
      xright, hface, local_ok)
    call require(.not. local_ok, "face transport rejects infinite temperature")
    call require(viscosity == 0.0_dp .and. conductivity == 0.0_dp .and. &
      all(diffusion == 0.0_dp) .and. all(yleft == 0.0_dp) .and. &
      all(yright == 0.0_dp) .and. all(yface == 0.0_dp) .and. &
      all(xleft == 0.0_dp) .and. all(xright == 0.0_dp) .and. &
      all(hface == 0.0_dp), &
      "infinite face failure rolls back every output")

    right_temperature = 1000.0_dp
    bad_transport = transport
    bad_transport(1)%diameter = ieee_value(0.0_dp, ieee_quiet_nan)
    call face_transport_data( &
      species, bad_transport, left, right, left_temperature, right_temperature, &
      viscosity, conductivity, diffusion, yleft, yright, yface, xleft, &
      xright, hface, local_ok)
    call require(.not. local_ok, "face transport rejects invalid coefficient")
    call require(viscosity == 0.0_dp .and. conductivity == 0.0_dp .and. &
      all(diffusion == 0.0_dp) .and. all(yleft == 0.0_dp) .and. &
      all(yright == 0.0_dp) .and. all(yface == 0.0_dp) .and. &
      all(xleft == 0.0_dp) .and. all(xright == 0.0_dp) .and. &
      all(hface == 0.0_dp), &
      "invalid face coefficient rolls back every output")
  end subroutine test_face_transport_fail_closed

  subroutine test_species_face_flux_fail_closed()
    real(dp), allocatable :: diffusion(:), yleft(:), yright(:), yface(:)
    real(dp), allocatable :: xleft(:), xright(:), hface(:), species_flux(:)
    real(dp) :: enthalpy_flux, density_face, pressure_left, pressure_right
    real(dp) :: pressure_face, spacing
    logical :: local_ok
    integer :: nspecies

    nspecies = size(species)
    allocate(diffusion(nspecies), yleft(nspecies), yright(nspecies))
    allocate(yface(nspecies), xleft(nspecies), xright(nspecies))
    allocate(hface(nspecies), species_flux(nspecies))
    diffusion = 1.0_dp
    yleft = 0.0_dp
    yright = 0.0_dp
    yface = 0.0_dp
    yface(1:2) = 0.5_dp
    xleft = 0.0_dp
    xright = 0.0_dp
    xright(1) = 1.0_dp
    hface = 1.0_dp
    density_face = 4.0_dp
    pressure_left = 101325.0_dp
    pressure_right = 101325.0_dp
    pressure_face = 101325.0_dp
    spacing = 1.0_dp
    call species_face_flux( &
      species, diffusion, yleft, yright, yface, xleft, xright, hface, &
      density_face, pressure_left, pressure_right, pressure_face, spacing, &
      .false., species_flux, enthalpy_flux, local_ok)
    call require(local_ok, "species face flux reference data")
    call require(all(ieee_is_finite(species_flux)) .and. &
      ieee_is_finite(enthalpy_flux), "species face flux reference finite")
    call require(abs(sum(species_flux)) < 1.0e-12_dp, &
      "species face flux reference closes")

    hface = 1.0_dp
    hface(1) = huge(1.0_dp)
    call species_face_flux( &
      species, diffusion, yleft, yright, yface, xleft, xright, hface, &
      density_face, pressure_left, pressure_right, pressure_face, spacing, &
      .false., species_flux, enthalpy_flux, local_ok)
    call require(.not. local_ok, "species face flux rejects extreme enthalpy")
    call require(all(species_flux == 0.0_dp) .and. enthalpy_flux == 0.0_dp, &
      "extreme enthalpy failure rolls back species flux")

    hface = 1.0_dp
    diffusion(1) = huge(1.0_dp)
    call species_face_flux( &
      species, diffusion, yleft, yright, yface, xleft, xright, hface, &
      density_face, pressure_left, pressure_right, pressure_face, spacing, &
      .false., species_flux, enthalpy_flux, local_ok)
    call require(.not. local_ok, "species face flux rejects extreme diffusion")
    call require(all(species_flux == 0.0_dp) .and. enthalpy_flux == 0.0_dp, &
      "extreme diffusion failure rolls back species flux")

    diffusion = 1.0_dp
    diffusion(1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call species_face_flux( &
      species, diffusion, yleft, yright, yface, xleft, xright, hface, &
      density_face, pressure_left, pressure_right, pressure_face, spacing, &
      .false., species_flux, enthalpy_flux, local_ok)
    call require(.not. local_ok, "species face flux rejects NaN diffusion")
    call require(all(species_flux == 0.0_dp) .and. enthalpy_flux == 0.0_dp, &
      "NaN diffusion failure rolls back species flux")

    diffusion(1) = ieee_value(0.0_dp, ieee_positive_inf)
    call species_face_flux( &
      species, diffusion, yleft, yright, yface, xleft, xright, hface, &
      density_face, pressure_left, pressure_right, pressure_face, spacing, &
      .false., species_flux, enthalpy_flux, local_ok)
    call require(.not. local_ok, "species face flux rejects infinite diffusion")
    call require(all(species_flux == 0.0_dp) .and. enthalpy_flux == 0.0_dp, &
      "infinite diffusion failure rolls back species flux")

    diffusion = 1.0_dp
    hface(1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call species_face_flux( &
      species, diffusion, yleft, yright, yface, xleft, xright, hface, &
      density_face, pressure_left, pressure_right, pressure_face, spacing, &
      .false., species_flux, enthalpy_flux, local_ok)
    call require(.not. local_ok, "species face flux rejects NaN enthalpy")
    call require(all(species_flux == 0.0_dp) .and. enthalpy_flux == 0.0_dp, &
      "NaN enthalpy failure rolls back species flux")

    hface(1) = ieee_value(0.0_dp, ieee_positive_inf)
    call species_face_flux( &
      species, diffusion, yleft, yright, yface, xleft, xright, hface, &
      density_face, pressure_left, pressure_right, pressure_face, spacing, &
      .false., species_flux, enthalpy_flux, local_ok)
    call require(.not. local_ok, "species face flux rejects infinite enthalpy")
    call require(all(species_flux == 0.0_dp) .and. enthalpy_flux == 0.0_dp, &
      "infinite enthalpy failure rolls back species flux")
  end subroutine test_species_face_flux_fail_closed

  subroutine test_public_validation_and_rollback()
    integer, parameter :: nx = 3, ny = 3
    integer :: nvar
    real(dp), allocatable :: state(:, :, :), temperature(:, :)
    real(dp), allocatable :: saved_state(:, :, :), saved_temperature(:, :)
    real(dp), allocatable :: bad_state(:, :, :), bad_temperature(:, :)
    real(dp), allocatable :: velocity_x(:, :), temperature_field(:, :)
    real(dp), allocatable :: flux_x(:, :, :), flux_y(:, :, :)
    real(dp), allocatable :: face_flux_x(:, :, :), face_flux_y(:, :, :)
    real(dp), allocatable :: output_state(:, :, :), output_temperature(:, :)
    real(dp) :: dt_value, maximum_diffusivity, theta
    real(dp) :: nan_value, inf_value
    logical :: local_ok

    nvar = reactive_nvar(size(species))
    allocate(state(nvar, nx, ny), temperature(nx, ny))
    allocate(saved_state(nvar, nx, ny), saved_temperature(nx, ny))
    allocate(bad_state(nvar, nx, ny), bad_temperature(nx, ny))
    allocate(velocity_x(nx, ny), temperature_field(nx, ny))
    allocate(flux_x(nvar, nx, ny), flux_y(nvar, nx, ny))
    allocate(face_flux_x(nvar, 0:nx, ny), face_flux_y(nvar, nx, 0:ny))
    allocate(output_state(nvar, nx, ny), output_temperature(nx, ny))

    velocity_x = 0.0_dp
    temperature_field = 1000.0_dp
    call build_test_state_2d( &
      state, temperature, velocity_x, temperature_field, local_ok)
    call require(local_ok, "public validation reference state")
    saved_state = state
    saved_temperature = temperature
    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    inf_value = ieee_value(0.0_dp, ieee_positive_inf)

    flux_x = 7.0_dp
    flux_y = 7.0_dp
    call reactive_transport_fluxes_2d( &
      species, transport, state, temperature, nx, ny, nan_value, 1.0e-3_dp, &
      1.0e-7_dp, .false., .false., .false., .false., flux_x, flux_y, theta, &
      local_ok)
    call require(.not. local_ok .and. all(flux_x == 0.0_dp) .and. &
      all(flux_y == 0.0_dp) .and. theta == 1.0_dp, &
      "NaN geometry rejects public fluxes")

    flux_x = 7.0_dp
    flux_y = 7.0_dp
    call reactive_transport_fluxes_2d( &
      species, transport, state, temperature, nx, ny, 1.0e-3_dp, inf_value, &
      1.0e-7_dp, .false., .false., .false., .false., flux_x, flux_y, theta, &
      local_ok)
    call require(.not. local_ok .and. all(flux_x == 0.0_dp) .and. &
      all(flux_y == 0.0_dp) .and. theta == 1.0_dp, &
      "infinite geometry rejects public fluxes")

    dt_value = 7.0_dp
    maximum_diffusivity = 8.0_dp
    call reactive_transport_timestep_2d( &
      species, transport, state, temperature, nx, ny, nan_value, 1.0e-3_dp, &
      0.35_dp, .false., .false., .false., dt_value, maximum_diffusivity, &
      local_ok)
    call require(.not. local_ok .and. dt_value == 0.0_dp .and. &
      maximum_diffusivity == 0.0_dp, &
      "NaN geometry rejects no-op timestep")

    output_state = 7.0_dp
    output_temperature = 8.0_dp
    call reactive_transport_euler_update_2d( &
      species, transport, state, temperature, nx, ny, nan_value, 1.0e-3_dp, &
      1.0e-7_dp, .false., .false., .false., .false., output_state, &
      output_temperature, theta, local_ok)
    call require(.not. local_ok .and. all(output_state == 0.0_dp) .and. &
      all(output_temperature == 0.0_dp) .and. theta == 1.0_dp, &
      "NaN geometry rolls back Euler outputs")

    state = saved_state
    temperature = saved_temperature
    call advance_reactive_transport_2d( &
      species, transport, state, temperature, nx, ny, nan_value, 1.0e-3_dp, &
      1.0e-7_dp, .false., .false., .false., .false., theta, local_ok)
    call require(.not. local_ok .and. all(state == saved_state) .and. &
      all(temperature == saved_temperature) .and. theta == 1.0_dp, &
      "NaN geometry rolls back advance state")

    bad_state = saved_state
    bad_temperature = saved_temperature
    bad_state(irho, 2, 1) = -1.0_dp
    dt_value = 7.0_dp
    maximum_diffusivity = 8.0_dp
    call reactive_transport_timestep_2d( &
      species, transport, bad_state, bad_temperature, nx, ny, 1.0e-3_dp, &
      1.0e-3_dp, 0.35_dp, .true., .true., .true., dt_value, &
      maximum_diffusivity, local_ok)
    call require(.not. local_ok .and. dt_value == 0.0_dp .and. &
      maximum_diffusivity == 0.0_dp, &
      "late timestep failure rolls back maximum diffusivity")

    output_state = 7.0_dp
    output_temperature = 8.0_dp
    call reactive_transport_euler_update_2d( &
      species, transport, bad_state, bad_temperature, nx, ny, 1.0e-3_dp, &
      1.0e-3_dp, 1.0e-7_dp, .true., .true., .true., .true., output_state, &
      output_temperature, theta, local_ok)
    call require(.not. local_ok .and. all(output_state == 0.0_dp) .and. &
      all(output_temperature == 0.0_dp) .and. theta == 1.0_dp, &
      "late Euler failure rolls back outputs")

    state = bad_state
    temperature = bad_temperature
    saved_state = state
    saved_temperature = temperature
    call advance_reactive_transport_2d( &
      species, transport, state, temperature, nx, ny, 1.0e-3_dp, 1.0e-3_dp, &
      1.0e-7_dp, .true., .true., .true., .true., theta, local_ok)
    call require(.not. local_ok .and. all(state == saved_state) .and. &
      all(temperature == saved_temperature) .and. theta == 1.0_dp, &
      "late advance failure rolls back state")

    velocity_x = 0.0_dp
    velocity_x(:, 2) = 1000.0_dp
    velocity_x(2, 1) = 1000.0_dp
    temperature_field = 1000.0_dp
    call build_test_state_2d( &
      state, temperature, velocity_x, temperature_field, local_ok)
    call require(local_ok, "stress guard reference state")
    face_flux_x = 7.0_dp
    face_flux_y = 7.0_dp
    call reactive_transport_fluxes_2d_faces( &
      species, transport, state, temperature, nx, ny, 1.0e-310_dp, 1.0e-2_dp, &
      1.0e-7_dp, .true., .false., .false., .false., face_flux_x, face_flux_y, &
      theta, local_ok)
    call require(.not. local_ok .and. all(face_flux_x == 0.0_dp) .and. &
      all(face_flux_y == 0.0_dp) .and. theta == 1.0_dp, &
      "overflowing x stress rolls back face fluxes")

    velocity_x = 0.0_dp
    temperature_field = 1000.0_dp
    temperature_field(2, 1) = 2000.0_dp
    call build_test_state_2d( &
      state, temperature, velocity_x, temperature_field, local_ok)
    call require(local_ok, "thermal guard reference state")
    face_flux_x = 7.0_dp
    face_flux_y = 7.0_dp
    call reactive_transport_fluxes_2d_faces( &
      species, transport, state, temperature, nx, ny, 1.0e-310_dp, 1.0e-2_dp, &
      1.0e-7_dp, .false., .true., .false., .false., face_flux_x, face_flux_y, &
      theta, local_ok)
    call require(.not. local_ok .and. all(face_flux_x == 0.0_dp) .and. &
      all(face_flux_y == 0.0_dp) .and. theta == 1.0_dp, &
      "overflowing x thermal gradient rolls back face fluxes")
  end subroutine test_public_validation_and_rollback

  subroutine build_test_state_2d( &
      state, temperature, velocity_x, temperature_field, ok)
    real(dp), intent(out) :: state(:, :, :), temperature(:, :)
    real(dp), intent(in) :: velocity_x(:, :), temperature_field(:, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:), mass_fractions(:)
    real(dp) :: mole_fractions(7), density, sound_speed
    logical :: local_ok
    integer :: i, j, k, nvar, nspecies

    ok = .false.
    nspecies = size(species)
    nvar = reactive_nvar(nspecies)
    state = 0.0_dp
    temperature = 0.0_dp
    if (size(state, 1) /= nvar .or. size(state, 2) /= size(velocity_x, 1) .or. &
        size(state, 3) /= size(velocity_x, 2) .or. &
        size(temperature, 1) /= size(velocity_x, 1) .or. &
        size(temperature, 2) /= size(velocity_x, 2) .or. &
        size(temperature_field, 1) /= size(velocity_x, 1) .or. &
        size(temperature_field, 2) /= size(velocity_x, 2)) return
    allocate(primitive(reactive_nprim(nspecies)), mass_fractions(nspecies))
    mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
      1.0e-5_dp, 0.0_dp, 0.55643_dp]
    call mass_fractions_from_mole_fractions( &
      species, mole_fractions, mass_fractions, local_ok)
    if (.not. local_ok) return
    do j = 1, size(velocity_x, 2)
      do i = 1, size(velocity_x, 1)
        density = mixture_density( &
          species, mass_fractions, 101325.0_dp, temperature_field(i, j), &
          local_ok)
        if (.not. local_ok) return
        primitive = 0.0_dp
        primitive(1:5) = [density, velocity_x(i, j), 0.0_dp, 0.0_dp, 101325.0_dp]
        do k = 1, nspecies
          primitive(reactive_mass_fraction_component(k)) = mass_fractions(k)
        end do
        call reactive_primitive_to_conserved( &
          species, primitive, state(:, i, j), temperature(i, j), sound_speed, &
          local_ok)
        if (.not. local_ok) return
      end do
    end do
    ok = .true.
  end subroutine build_test_state_2d

  logical function ieee_safe(value)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(dp), intent(in) :: value
    ieee_safe = ieee_is_finite(value)
  end function ieee_safe

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') "FAIL: " // trim(message)
      error stop
    end if
  end subroutine require
end program test_reactive_transport_2d
