program test_eb_reactive_flux_divergence_2d
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use precision_mod, only: dp
  use state_indices_mod, only: imx, imy
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved
  use eb_geometry_2d_mod, only: &
    eb_geometry_2d, eb_covered_cell, build_eb_geometry_2d
  use eb_reactive_wall_flux_2d_mod, only: &
    reactive_eb_flux_divergence_2d
  implicit none

  integer, parameter :: nx = 20
  integer, parameter :: ny = 20
  real(dp), parameter :: pressure = 135000.0_dp
  type(nasa7_species), allocatable :: species(:)
  type(eb_geometry_2d) :: geometry
  real(dp), allocatable :: primitive(:), state_cell(:), mass_fractions(:)
  real(dp), allocatable :: state(:, :, :), temperature_field(:, :)
  real(dp), allocatable :: x_flux(:, :, :), y_flux(:, :, :), rhs(:, :, :)
  real(dp), allocatable :: short_x_flux(:, :, :)
  real(dp) :: level_set(0:nx, 0:ny), mole_fractions(7)
  real(dp) :: temperature, sound_speed, x, y
  logical :: ok
  integer :: i, j, k, nvar

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamic database load")
  nvar = reactive_nvar(size(species))
  allocate(primitive(reactive_nprim(size(species))))
  allocate(state_cell(nvar), mass_fractions(size(species)))
  allocate(state(nvar, nx, ny), temperature_field(nx, ny))
  allocate(x_flux(nvar, 0:nx, ny), y_flux(nvar, nx, 0:ny))
  allocate(rhs(nvar, nx, ny))

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

  level_set = 1.0_dp
  call verify_uniform_pressure_balance(level_set, "regular domain")

  do j = 0, ny
    do i = 0, nx
      level_set(i, j) = real(i, dp) / real(nx, dp) - 0.37_dp
    end do
  end do
  call verify_uniform_pressure_balance(level_set, "vertical wall")

  do j = 0, ny
    do i = 0, nx
      level_set(i, j) = real(i, dp) / real(nx, dp) + &
        real(j, dp) / real(ny, dp) - 0.8_dp
    end do
  end do
  call verify_uniform_pressure_balance(level_set, "diagonal wall")

  do j = 0, ny
    y = real(j, dp) / real(ny, dp)
    do i = 0, nx
      x = real(i, dp) / real(nx, dp)
      level_set(i, j) = 0.27_dp - &
        sqrt((x - 0.5_dp)**2 + (y - 0.5_dp)**2)
    end do
  end do
  call verify_uniform_pressure_balance(level_set, "circular wall")

  x_flux(imx, 0, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
  call reactive_eb_flux_divergence_2d( &
    species, state, temperature_field, geometry, x_flux, y_flux, rhs, ok)
  call require(.not. ok .and. maxval(abs(rhs)) == 0.0_dp, &
    "nonfinite face flux transaction")

  allocate(short_x_flux(nvar, 0:nx - 1, ny))
  short_x_flux = 0.0_dp
  call reactive_eb_flux_divergence_2d( &
    species, state, temperature_field, geometry, short_x_flux, y_flux, rhs, ok)
  call require(.not. ok .and. maxval(abs(rhs)) == 0.0_dp, &
    "face-array extent transaction")

  write(*, '(a)') "test_eb_reactive_flux_divergence_2d: PASS"

contains

  subroutine verify_uniform_pressure_balance(node_level_set, label)
    real(dp), intent(in) :: node_level_set(0:nx, 0:ny)
    character(len=*), intent(in) :: label

    real(dp) :: integrated_residual, maximum_integrated_residual
    real(dp) :: residual_tolerance
    integer :: local_i, local_j

    call build_eb_geometry_2d( &
      node_level_set, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, geometry, ok)
    call require(ok .and. geometry%is_valid(), trim(label)//" geometry")
    x_flux = 0.0_dp
    y_flux = 0.0_dp
    x_flux(imx, :, :) = pressure
    y_flux(imy, :, :) = pressure
    call reactive_eb_flux_divergence_2d( &
      species, state, temperature_field, geometry, x_flux, y_flux, rhs, ok)
    call require(ok, trim(label)//" divergence")

    maximum_integrated_residual = 0.0_dp
    do local_j = 1, ny
      do local_i = 1, nx
        if (geometry%cell_type(local_i, local_j) == eb_covered_cell) then
          call require(maxval(abs(rhs(:, local_i, local_j))) == 0.0_dp, &
            trim(label)//" covered-cell residual")
          cycle
        end if
        integrated_residual = maxval(abs(rhs(:, local_i, local_j))) * &
          geometry%volume_fraction(local_i, local_j) * &
          geometry%dx * geometry%dy
        maximum_integrated_residual = max( &
          maximum_integrated_residual, integrated_residual)
      end do
    end do
    residual_tolerance = 2.0e-10_dp * pressure * &
      max(geometry%dx, geometry%dy)
    call require(maximum_integrated_residual <= residual_tolerance, &
      trim(label)//" uniform-pressure balance")
  end subroutine verify_uniform_pressure_balance

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_eb_reactive_flux_divergence_2d
