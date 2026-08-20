program test_reactive_transport_2d
  use precision_mod, only: dp
  use state_indices_mod, only: imx, imy
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved, advance_reactive_transport
  use reactive_transport_2d_mod, only: &
    reactive_transport_fluxes_2d, advance_reactive_transport_2d
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
