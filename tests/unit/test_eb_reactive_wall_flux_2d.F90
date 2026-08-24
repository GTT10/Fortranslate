program test_eb_reactive_wall_flux_2d
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_species_component, reactive_primitive_to_conserved
  use eb_geometry_2d_mod, only: &
    eb_geometry_2d, eb_cut_cell, build_eb_geometry_2d
  use eb_reactive_wall_flux_2d_mod, only: &
    reactive_eb_slip_wall_flux_2d, reactive_eb_slip_wall_source_2d
  implicit none

  integer, parameter :: nx = 10
  integer, parameter :: ny = 10
  real(dp), parameter :: pressure = 135000.0_dp
  type(nasa7_species), allocatable :: species(:)
  type(eb_geometry_2d) :: geometry
  real(dp), allocatable :: primitive(:), rotated_primitive(:)
  real(dp), allocatable :: state_cell(:), rotated_state(:), changed_state(:)
  real(dp), allocatable :: flux(:), rotated_flux(:), changed_flux(:)
  real(dp), allocatable :: state(:, :, :), source(:, :, :)
  real(dp), allocatable :: temperature_field(:, :), mass_fractions(:)
  real(dp) :: level_set(0:nx, 0:ny)
  real(dp) :: mole_fractions(7), normal(2), rotated_normal(2)
  real(dp) :: temperature, rotated_temperature, changed_temperature
  real(dp) :: sound_speed, rotated_sound_speed, changed_sound_speed
  real(dp) :: wall_pressure, rotated_pressure, changed_pressure
  real(dp) :: integrated_force_x, integrated_force_y, tolerance
  logical :: ok
  integer :: i, j, k, nvar, species_first

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamic database load")
  nvar = reactive_nvar(size(species))
  allocate(primitive(reactive_nprim(size(species))))
  allocate(rotated_primitive(size(primitive)))
  allocate(state_cell(nvar), rotated_state(nvar), changed_state(nvar))
  allocate(flux(nvar), rotated_flux(nvar), changed_flux(nvar))
  allocate(state(nvar, nx, ny), source(nvar, nx, ny))
  allocate(temperature_field(nx, ny), mass_fractions(size(species)))

  mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.55643_dp]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  call require(ok, "composition conversion")
  primitive(1:5) = [0.31_dp, 42.0_dp, -27.0_dp, 6.0_dp, pressure]
  do k = 1, size(species)
    primitive(reactive_mass_fraction_component(k)) = mass_fractions(k)
  end do
  call reactive_primitive_to_conserved( &
    species, primitive, state_cell, temperature, sound_speed, ok)
  call require(ok, "reference state construction")

  normal = [0.6_dp, 0.8_dp]
  call reactive_eb_slip_wall_flux_2d( &
    species, state_cell, temperature, normal, flux, wall_pressure, ok)
  call require(ok, "oblique wall flux")
  tolerance = 2.0e-9_dp * pressure
  call assert_close(wall_pressure, pressure, tolerance, &
    "wall pressure recovery")
  call assert_close(flux(imx), -0.6_dp * pressure, tolerance, &
    "wall x-momentum flux")
  call assert_close(flux(imy), -0.8_dp * pressure, tolerance, &
    "wall y-momentum flux")
  call require(flux(irho) == 0.0_dp .and. flux(imz) == 0.0_dp .and. &
    flux(iet) == 0.0_dp, "wall base-state impermeability")
  species_first = reactive_species_component(1)
  call require(maxval(abs(flux(species_first:nvar))) == 0.0_dp, &
    "wall species impermeability")

  rotated_primitive = primitive
  rotated_primitive(2) = -primitive(3)
  rotated_primitive(3) = primitive(2)
  call reactive_primitive_to_conserved( &
    species, rotated_primitive, rotated_state, rotated_temperature, &
    rotated_sound_speed, ok)
  call require(ok, "rotated state construction")
  rotated_normal = [-normal(2), normal(1)]
  call reactive_eb_slip_wall_flux_2d( &
    species, rotated_state, rotated_temperature, rotated_normal, &
    rotated_flux, rotated_pressure, ok)
  call require(ok, "rotated wall flux")
  call assert_close(rotated_flux(imx), -flux(imy), tolerance, &
    "rotated x-momentum flux")
  call assert_close(rotated_flux(imy), flux(imx), tolerance, &
    "rotated y-momentum flux")
  call assert_close(rotated_pressure, wall_pressure, tolerance, &
    "rotated wall pressure")

  rotated_primitive(2:4) = [-113.0_dp, 87.0_dp, -19.0_dp]
  call reactive_primitive_to_conserved( &
    species, rotated_primitive, changed_state, changed_temperature, &
    changed_sound_speed, ok)
  call require(ok, "changed-velocity state construction")
  call reactive_eb_slip_wall_flux_2d( &
    species, changed_state, changed_temperature, normal, changed_flux, &
    changed_pressure, ok)
  call require(ok, "changed-velocity wall flux")
  call require(maxval(abs(changed_flux - flux)) <= tolerance, &
    "slip wall flux velocity independence")

  do j = 0, ny
    do i = 0, nx
      level_set(i, j) = real(i, dp) / real(nx, dp) - 0.37_dp
    end do
  end do
  call build_eb_geometry_2d( &
    level_set, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, geometry, ok)
  call require(ok .and. geometry%is_valid(), "vertical EB geometry")
  do j = 1, ny
    do i = 1, nx
      state(:, i, j) = state_cell
      temperature_field(i, j) = temperature
    end do
  end do
  call reactive_eb_slip_wall_source_2d( &
    species, state, temperature_field, geometry, source, ok)
  call require(ok, "vertical EB wall source")
  call require(count(geometry%cell_type == eb_cut_cell) == ny, &
    "vertical cut-cell count")
  call assert_close(maxval(abs(source(imx, 4, :) - &
    pressure / (0.30_dp * geometry%dx))), 0.0_dp, &
    8.0e-9_dp * pressure / geometry%dx, "cut-cell source density")
  call require(maxval(abs(source(:, 1:3, :))) == 0.0_dp .and. &
    maxval(abs(source(:, 5:nx, :))) == 0.0_dp, &
    "wall source restricted to cut cells")
  call require(maxval(abs(source(irho, :, :))) == 0.0_dp .and. &
    maxval(abs(source(imz, :, :))) == 0.0_dp .and. &
    maxval(abs(source(iet, :, :))) == 0.0_dp .and. &
    maxval(abs(source(species_first:nvar, :, :))) == 0.0_dp, &
    "wall source has no mass energy or species leakage")
  integrated_force_x = sum(source(imx, :, :) * &
    geometry%volume_fraction) * geometry%dx * geometry%dy
  integrated_force_y = sum(source(imy, :, :) * &
    geometry%volume_fraction) * geometry%dx * geometry%dy
  call assert_close(integrated_force_x, pressure, tolerance, &
    "vertical integrated pressure force x")
  call assert_close(integrated_force_y, 0.0_dp, tolerance, &
    "vertical integrated pressure force y")

  do j = 0, ny
    do i = 0, nx
      level_set(i, j) = real(i, dp) / real(nx, dp) + &
        real(j, dp) / real(ny, dp) - 0.8_dp
    end do
  end do
  call build_eb_geometry_2d( &
    level_set, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, geometry, ok)
  call require(ok .and. geometry%is_valid(), "diagonal EB geometry")
  call reactive_eb_slip_wall_source_2d( &
    species, state, temperature_field, geometry, source, ok)
  call require(ok, "diagonal EB wall source")
  integrated_force_x = sum(source(imx, :, :) * &
    geometry%volume_fraction) * geometry%dx * geometry%dy
  integrated_force_y = sum(source(imy, :, :) * &
    geometry%volume_fraction) * geometry%dx * geometry%dy
  call assert_close(integrated_force_x, 0.8_dp * pressure, tolerance, &
    "diagonal integrated pressure force x")
  call assert_close(integrated_force_y, 0.8_dp * pressure, tolerance, &
    "diagonal integrated pressure force y")

  call reactive_eb_slip_wall_flux_2d( &
    species, state_cell, temperature, [2.0_dp, 0.0_dp], flux, &
    wall_pressure, ok)
  call require(.not. ok .and. maxval(abs(flux)) == 0.0_dp, &
    "non-unit normal rejection")

  state(:, 4, 5) = state_cell
  state(irho, 4, 5) = -1.0_dp
  do j = 0, ny
    do i = 0, nx
      level_set(i, j) = real(i, dp) / real(nx, dp) - 0.37_dp
    end do
  end do
  call build_eb_geometry_2d( &
    level_set, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, geometry, ok)
  call require(ok, "transaction geometry")
  call reactive_eb_slip_wall_source_2d( &
    species, state, temperature_field, geometry, source, ok)
  call require(.not. ok .and. maxval(abs(source)) == 0.0_dp, &
    "failed wall source is transactional")

  write(*, '(a)') "test_eb_reactive_wall_flux_2d: PASS"

contains

  subroutine assert_close(actual, expected, local_tolerance, message)
    real(dp), intent(in) :: actual, expected, local_tolerance
    character(len=*), intent(in) :: message

    call require(abs(actual - expected) <= local_tolerance, message)
  end subroutine assert_close

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_eb_reactive_wall_flux_2d
