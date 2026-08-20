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

  integer, parameter :: nx = 24, ny = 4
  type(nasa7_species), allocatable :: species(:)
  real(dp), allocatable :: state_1d(:, :), state_2d(:, :, :)
  real(dp), allocatable :: temperature_1d(:), temperature_2d(:, :)
  real(dp), allocatable :: primitive(:), mass_fractions(:)
  real(dp) :: mole_fractions(7), dx, dy, dt, x, rho0, rho
  real(dp) :: local_temperature, sound_speed, theta, scale, difference
  logical :: ok
  integer :: i, j, k, nvar

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "thermodynamic database load failed"
  nvar = reactive_nvar(size(species))
  allocate(state_1d(nvar, 0:nx + 1), temperature_1d(0:nx + 1))
  allocate(state_2d(nvar, nx, ny), temperature_2d(nx, ny))
  allocate(primitive(reactive_nprim(size(species))))
  allocate(mass_fractions(size(species)))
  mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.55643_dp]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  if (.not. ok) error stop "composition conversion failed"
  rho0 = mixture_density(species, mass_fractions, 101325.0_dp, 1000.0_dp, ok)
  if (.not. ok) error stop "base density failed"
  dx = 0.01_dp / real(nx, dp)
  dy = 0.002_dp / real(ny, dp)
  do i = 1, nx
    x = (real(i, dp) - 0.5_dp) * dx
    rho = rho0 * (1.0_dp + 0.06_dp * sin(2.0_dp * acos(-1.0_dp) * x / 0.01_dp))
    primitive(1:5) = [rho, 250.0_dp, 0.0_dp, 0.0_dp, 101325.0_dp]
    do k = 1, size(species)
      primitive(reactive_mass_fraction_component(k)) = mass_fractions(k)
    end do
    call reactive_primitive_to_conserved( &
      species, primitive, state_1d(:, i), local_temperature, sound_speed, ok)
    if (.not. ok) error stop "state construction failed"
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
    "characteristic_plm", "mc", "periodic", ok, "hllc")
  if (.not. ok) error stop "1D reference update failed"
  call advance_reactive_hydro_2d( &
    species, state_2d, temperature_2d, nx, ny, dx, dy, dt, &
    "characteristic_plm", "mc", "hllc", .true., ok, theta)
  if (.not. ok) error stop "2D CTU update failed"
  if (theta < 0.999999999_dp) error stop "unexpected transverse limiter activation"

  difference = 0.0_dp
  scale = 1.0_dp
  do j = 1, ny
    do i = 1, nx
      difference = max(difference, maxval(abs(state_2d(:, i, j) - state_1d(:, i))))
      scale = max(scale, maxval(abs(state_1d(:, i))))
      difference = max(difference, abs(temperature_2d(i, j) - temperature_1d(i)))
      scale = max(scale, abs(temperature_1d(i)))
    end do
  end do
  if (difference / scale > 3.0e-12_dp) error stop "2D-to-1D reduction mismatch"
end program test_reactive_ctu_dimensional_reduction
