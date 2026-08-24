program test_eb_reactive_hydro_2d
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use precision_mod, only: dp
  use state_indices_mod, only: imx, imy
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved, reactive_riemann_flux_x
  use eb_geometry_2d_mod, only: &
    eb_geometry_2d, eb_covered_cell, build_eb_geometry_2d
  use eb_reactive_hydro_2d_mod, only: &
    reactive_eb_outflow_riemann_fluxes_2d, &
    advance_reactive_eb_hydro_2d
  implicit none

  integer, parameter :: nx = 20
  integer, parameter :: ny = 20
  real(dp), parameter :: pressure = 135000.0_dp
  type(nasa7_species), allocatable :: species(:)
  type(eb_geometry_2d) :: geometry
  real(dp), allocatable :: primitive(:), state_cell(:), second_state(:)
  real(dp), allocatable :: mass_fractions(:), face_flux(:)
  real(dp), allocatable :: state(:, :, :), temperature_field(:, :)
  real(dp), allocatable :: new_state(:, :, :), new_temperature(:, :)
  real(dp), allocatable :: x_flux(:, :, :), y_flux(:, :, :)
  real(dp) :: level_set(0:nx, 0:ny), mole_fractions(7)
  real(dp) :: temperature, second_temperature, sound_speed, second_sound_speed
  real(dp) :: x, y, dt, tolerance
  logical :: ok
  integer :: i, j, k, nvar

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamic database load")
  nvar = reactive_nvar(size(species))
  allocate(primitive(reactive_nprim(size(species))))
  allocate(state_cell(nvar), second_state(nvar), face_flux(nvar))
  allocate(mass_fractions(size(species)))
  allocate(state(nvar, nx, ny), new_state(nvar, nx, ny))
  allocate(temperature_field(nx, ny), new_temperature(nx, ny))
  allocate(x_flux(nvar, 0:nx, ny), y_flux(nvar, nx, 0:ny))

  mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.55643_dp]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  call require(ok, "composition conversion")
  primitive(1:5) = [0.31_dp, 0.0_dp, 0.0_dp, 0.0_dp, pressure]
  do k = 1, size(species)
    primitive(reactive_mass_fraction_component(k)) = mass_fractions(k)
  end do
  call reactive_primitive_to_conserved( &
    species, primitive, state_cell, temperature, sound_speed, ok)
  call require(ok, "stationary state construction")
  do j = 1, ny
    do i = 1, nx
      state(:, i, j) = state_cell
      temperature_field(i, j) = temperature
    end do
  end do
  dt = 0.1_dp / sound_speed / real(nx, dp)

  level_set = 1.0_dp
  call build_eb_geometry_2d( &
    level_set, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, geometry, ok)
  call require(ok .and. geometry%is_valid(), "regular geometry")
  primitive(1:5) = [0.35_dp, 60.0_dp, -20.0_dp, 5.0_dp, 150000.0_dp]
  call reactive_primitive_to_conserved( &
    species, primitive, second_state, second_temperature, &
    second_sound_speed, ok)
  call require(ok, "second state construction")
  state(:, 2, 1) = second_state
  temperature_field(2, 1) = second_temperature
  call reactive_eb_outflow_riemann_fluxes_2d( &
    species, state, temperature_field, geometry, "hllc", x_flux, y_flux, ok)
  call require(ok, "regular Riemann flux construction")
  call reactive_riemann_flux_x( &
    species, state_cell, second_state, temperature, second_temperature, &
    "hllc", face_flux, ok)
  call require(ok, "reference interior Riemann flux")
  tolerance = 3.0e-12_dp * max(1.0_dp, maxval(abs(face_flux)))
  call require(maxval(abs(x_flux(:, 1, 1) - face_flux)) <= tolerance, &
    "interior face Riemann parity")
  call reactive_riemann_flux_x( &
    species, state_cell, state_cell, temperature, temperature, &
    "hllc", face_flux, ok)
  call require(ok, "reference boundary Riemann flux")
  tolerance = 3.0e-12_dp * max(1.0_dp, maxval(abs(face_flux)))
  call require(maxval(abs(x_flux(:, 0, 1) - face_flux)) <= tolerance, &
    "zero-gradient boundary face")
  state(:, 2, 1) = state_cell
  temperature_field(2, 1) = temperature

  level_set = 1.0_dp
  call verify_uniform_hydro_step(level_set, "regular domain")

  do j = 0, ny
    do i = 0, nx
      level_set(i, j) = real(i, dp) / real(nx, dp) - 0.37_dp
    end do
  end do
  call verify_uniform_hydro_step(level_set, "vertical wall")

  do j = 0, ny
    do i = 0, nx
      level_set(i, j) = real(i, dp) / real(nx, dp) + &
        real(j, dp) / real(ny, dp) - 0.8_dp
    end do
  end do
  call verify_uniform_hydro_step(level_set, "diagonal wall")

  do j = 0, ny
    y = real(j, dp) / real(ny, dp)
    do i = 0, nx
      x = real(i, dp) / real(nx, dp)
      level_set(i, j) = 0.27_dp - &
        sqrt((x - 0.5_dp)**2 + (y - 0.5_dp)**2)
    end do
  end do
  call verify_uniform_hydro_step(level_set, "circular wall")

  call reactive_eb_outflow_riemann_fluxes_2d( &
    species, state, temperature_field, geometry, "unknown", &
    x_flux, y_flux, ok)
  call require(.not. ok .and. maxval(abs(x_flux)) == 0.0_dp .and. &
    maxval(abs(y_flux)) == 0.0_dp, "unknown solver flux transaction")
  call advance_reactive_eb_hydro_2d( &
    species, state, temperature_field, geometry, "unknown", dt, &
    new_state, new_temperature, ok)
  call require(.not. ok .and. &
    maxval(abs(new_state - state)) == 0.0_dp .and. &
    maxval(abs(new_temperature - temperature_field)) == 0.0_dp, &
    "unknown solver advance transaction")

  call require(geometry%cell_type(1, 1) == eb_covered_cell, &
    "circular corner is covered")
  x = geometry%x_face_fraction(0, 1)
  geometry%x_face_fraction(0, 1) = 1.0_dp
  call require(geometry%is_valid(), "inconsistent face test geometry")
  call reactive_eb_outflow_riemann_fluxes_2d( &
    species, state, temperature_field, geometry, "hllc", x_flux, y_flux, ok)
  call require(.not. ok .and. maxval(abs(x_flux)) == 0.0_dp .and. &
    maxval(abs(y_flux)) == 0.0_dp, "open face beside covered cell")
  geometry%x_face_fraction(0, 1) = x

  state(1, 1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
  call reactive_eb_outflow_riemann_fluxes_2d( &
    species, state, temperature_field, geometry, "hllc", x_flux, y_flux, ok)
  call require(.not. ok .and. maxval(abs(x_flux)) == 0.0_dp .and. &
    maxval(abs(y_flux)) == 0.0_dp, "nonfinite state transaction")

  write(*, '(a)') "test_eb_reactive_hydro_2d: PASS"

contains

  subroutine verify_uniform_hydro_step(node_level_set, label)
    real(dp), intent(in) :: node_level_set(0:nx, 0:ny)
    character(len=*), intent(in) :: label

    real(dp) :: extensive_change, change_tolerance, flux_tolerance
    integer :: local_i, local_j

    call build_eb_geometry_2d( &
      node_level_set, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, geometry, ok)
    call require(ok .and. geometry%is_valid(), trim(label)//" geometry")
    call reactive_eb_outflow_riemann_fluxes_2d( &
      species, state, temperature_field, geometry, "hllc", &
      x_flux, y_flux, ok)
    call require(ok, trim(label)//" Riemann fluxes")
    flux_tolerance = 2.0e-10_dp * pressure
    face_flux = 0.0_dp
    face_flux(imx) = pressure
    do local_j = 1, ny
      do local_i = 0, nx
        if (geometry%x_face_fraction(local_i, local_j) > 0.0_dp) then
          call require(maxval(abs(x_flux(:, local_i, local_j) - &
            face_flux)) <= &
            flux_tolerance, trim(label)//" open x-face pressure flux")
        else
          call require(maxval(abs(x_flux(:, local_i, local_j))) == 0.0_dp, &
            trim(label)//" closed x-face flux")
        end if
      end do
    end do
    face_flux = 0.0_dp
    face_flux(imy) = pressure
    do local_j = 0, ny
      do local_i = 1, nx
        if (geometry%y_face_fraction(local_i, local_j) > 0.0_dp) then
          call require(maxval(abs(y_flux(:, local_i, local_j) - &
            face_flux)) <= &
            flux_tolerance, trim(label)//" open y-face pressure flux")
        else
          call require(maxval(abs(y_flux(:, local_i, local_j))) == 0.0_dp, &
            trim(label)//" closed y-face flux")
        end if
      end do
    end do

    call advance_reactive_eb_hydro_2d( &
      species, state, temperature_field, geometry, "hllc", dt, &
      new_state, new_temperature, ok)
    call require(ok, trim(label)//" complete hydro step")
    change_tolerance = dt * 2.0e-10_dp * pressure * &
      max(geometry%dx, geometry%dy)
    do local_j = 1, ny
      do local_i = 1, nx
        if (geometry%cell_type(local_i, local_j) == eb_covered_cell) then
          call require(maxval(abs(new_state(:, local_i, local_j) - &
            state(:, local_i, local_j))) == 0.0_dp, &
            trim(label)//" covered state")
          call require(new_temperature(local_i, local_j) == &
            temperature_field(local_i, local_j), &
            trim(label)//" covered temperature")
          cycle
        end if
        extensive_change = maxval(abs(new_state(:, local_i, local_j) - &
          state(:, local_i, local_j))) * &
          geometry%volume_fraction(local_i, local_j) * &
          geometry%dx * geometry%dy
        call require(extensive_change <= change_tolerance, &
          trim(label)//" stationary extensive update")
        call require(abs(new_temperature(local_i, local_j) - temperature) <= &
          2.0e-8_dp, trim(label)//" recovered temperature")
      end do
    end do
  end subroutine verify_uniform_hydro_step

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_eb_reactive_hydro_2d
