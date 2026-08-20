program test_reactive_ctu_dimensional_reduction
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved, advance_reactive_hydro
  use reactive_2d_mod, only: advance_reactive_hydro_2d
  implicit none

  type(nasa7_species), allocatable :: species(:)
  logical :: ok

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "thermodynamic database load failed"
  call run_mode("characteristic_plm")
  call run_mode("characteristic_ppm")

contains

  subroutine run_mode(reconstruction)
    character(len=*), intent(in) :: reconstruction
    integer, parameter :: nx = 24, ny = 4
    real(dp), allocatable :: state_1d(:, :), state_2d(:, :, :)
    real(dp), allocatable :: temperature_1d(:), temperature_2d(:, :)
    real(dp), allocatable :: primitive(:), mass_fractions(:)
    real(dp) :: mole_fractions(7), dx, dy, dt, x, rho0, rho
    real(dp) :: local_temperature, sound_speed, theta, scale, difference
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
    if (.not. local_ok) error stop "composition conversion failed"
    rho0 = mixture_density( &
      species, mass_fractions, 101325.0_dp, 1000.0_dp, local_ok)
    if (.not. local_ok) error stop "base density failed"
    dx = 0.01_dp / real(nx, dp)
    dy = 0.002_dp / real(ny, dp)
    do i = 1, nx
      x = (real(i, dp) - 0.5_dp) * dx
      rho = rho0 * &
        (1.0_dp + 0.06_dp * sin(2.0_dp * acos(-1.0_dp) * x / 0.01_dp))
      primitive(1:5) = [rho, 250.0_dp, 0.0_dp, 0.0_dp, 101325.0_dp]
      do k = 1, size(species)
        primitive(reactive_mass_fraction_component(k)) = mass_fractions(k)
      end do
      call reactive_primitive_to_conserved( &
        species, primitive, state_1d(:, i), local_temperature, sound_speed, &
        local_ok)
      if (.not. local_ok) error stop "state construction failed"
      temperature_1d(i) = local_temperature
      do j = 1, ny
        state_2d(:, i, j) = state_1d(:, i)
        temperature_2d(i, j) = temperature_1d(i)
      end do
    end do
    state_1d(:, 0) = state_1d(:, nx)
    state_1d(:, nx + 1) = state_1d(:, 1)
    temperature_1d(0) = temperature_1d(nx)
    temperature_1d(nx + 1) = temperature_1d(1)
    dt = 0.18_dp * dx / (250.0_dp + sound_speed)

    call advance_reactive_hydro( &
      species, state_1d, temperature_1d, nx, dx, dt, &
      reconstruction, "mc", "periodic", local_ok, "hllc")
    if (.not. local_ok) error stop "1D reference update failed"
    call advance_reactive_hydro_2d( &
      species, state_2d, temperature_2d, nx, ny, dx, dy, dt, &
      reconstruction, "mc", "hllc", .true., local_ok, theta)
    if (.not. local_ok) error stop "2D CTU update failed"
    if (theta < 0.999999999_dp) &
      error stop "unexpected transverse limiter activation"

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
    if (difference / scale > 3.0e-12_dp) then
      write(*, '(a,1x,a,1x,es16.8)') &
        "2D-to-1D reduction mismatch", trim(reconstruction), difference / scale
      error stop
    end if
    call run_y_mode(reconstruction, mass_fractions, rho0, sound_speed)
  end subroutine run_mode

  subroutine run_y_mode(reconstruction, mass_fractions, rho0, sound_speed_guess)
    character(len=*), intent(in) :: reconstruction
    real(dp), intent(in) :: mass_fractions(:), rho0, sound_speed_guess
    integer, parameter :: nx = 4, ny = 24
    real(dp), allocatable :: state_1d(:, :), state_2d(:, :, :)
    real(dp), allocatable :: temperature_1d(:), temperature_2d(:, :)
    real(dp), allocatable :: primitive(:), rotated(:)
    real(dp) :: dy, dx, dt, y, rho, local_temperature, sound_speed
    real(dp) :: theta, scale, difference
    logical :: local_ok
    integer :: i, j, k, nvar

    nvar = reactive_nvar(size(species))
    allocate(state_1d(nvar, 0:ny + 1), temperature_1d(0:ny + 1))
    allocate(state_2d(nvar, nx, ny), temperature_2d(nx, ny))
    allocate(primitive(reactive_nprim(size(species))), rotated(nvar))
    dx = 0.002_dp / real(nx, dp)
    dy = 0.01_dp / real(ny, dp)
    do j = 1, ny
      y = (real(j, dp) - 0.5_dp) * dy
      rho = rho0 * &
        (1.0_dp + 0.06_dp * sin(2.0_dp * acos(-1.0_dp) * y / 0.01_dp))
      ! The 1D reference uses its x velocity as the normal component.  The
      ! physical 2D state rotates that normal velocity into v and keeps u as
      ! the first transverse component.
      primitive(1:5) = [rho, 180.0_dp, 40.0_dp, 0.0_dp, 101325.0_dp]
      do k = 1, size(species)
        primitive(reactive_mass_fraction_component(k)) = mass_fractions(k)
      end do
      call reactive_primitive_to_conserved( &
        species, primitive, state_1d(:, j), local_temperature, sound_speed, &
        local_ok)
      if (.not. local_ok) error stop "y-reference state construction failed"
      temperature_1d(j) = local_temperature
      rotated = state_1d(:, j)
      rotated(2) = state_1d(3, j)
      rotated(3) = state_1d(2, j)
      do i = 1, nx
        state_2d(:, i, j) = rotated
        temperature_2d(i, j) = local_temperature
      end do
    end do
    state_1d(:, 0) = state_1d(:, ny)
    state_1d(:, ny + 1) = state_1d(:, 1)
    temperature_1d(0) = temperature_1d(ny)
    temperature_1d(ny + 1) = temperature_1d(1)
    dt = 0.18_dp * dy / (180.0_dp + max(sound_speed, sound_speed_guess))

    call advance_reactive_hydro( &
      species, state_1d, temperature_1d, ny, dy, dt, reconstruction, &
      "mc", "periodic", local_ok, "hllc")
    if (.not. local_ok) error stop "1D y-reference update failed"
    call advance_reactive_hydro_2d( &
      species, state_2d, temperature_2d, nx, ny, dx, dy, dt, &
      reconstruction, "mc", "hllc", .true., local_ok, theta)
    if (.not. local_ok) error stop "2D y-direction CTU update failed"
    if (theta < 0.999999999_dp) &
      error stop "unexpected y-direction transverse limiter activation"

    difference = 0.0_dp
    scale = 1.0_dp
    do j = 1, ny
      rotated = state_1d(:, j)
      rotated(2) = state_1d(3, j)
      rotated(3) = state_1d(2, j)
      do i = 1, nx
        difference = max(difference, &
          maxval(abs(state_2d(:, i, j) - rotated)))
        scale = max(scale, maxval(abs(rotated)))
        difference = max(difference, &
          abs(temperature_2d(i, j) - temperature_1d(j)))
        scale = max(scale, abs(temperature_1d(j)))
      end do
    end do
    if (difference / scale > 3.0e-12_dp) then
      write(*, '(a,1x,a,1x,es16.8)') &
        "2D y-to-1D reduction mismatch", trim(reconstruction), difference / scale
      error stop
    end if
  end subroutine run_y_mode
end program test_reactive_ctu_dimensional_reduction
