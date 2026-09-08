program test_reactive_eb_cfl_3d
  use precision_mod, only: dp
  use state_indices_mod, only: irho
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, build_axis_plane_eb_geometry_3d
  use reactive_eb_cfl_3d_mod, only: &
    compute_reactive_eb_cfl_timestep_3d
  implicit none

  integer, parameter :: nx = 8
  integer, parameter :: ny = 5
  integer, parameter :: nz = 4
  type(nasa7_species), allocatable :: species(:)
  type(eb_geometry_3d) :: geometry
  real(dp), allocatable :: primitive(:), conserved(:), mass_fractions(:)
  real(dp), allocatable :: state(:, :, :, :), temperature(:, :, :)
  real(dp) :: mole_fractions(7), recovered_temperature, sound_speed
  real(dp) :: expected_rate, expected_dt, dt
  logical :: ok
  integer :: component, species_index, nvar

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "CFL thermodynamic database load")
  nvar = reactive_nvar(size(species))
  allocate(primitive(reactive_nprim(size(species))))
  allocate(conserved(nvar), mass_fractions(size(species)))
  allocate(state(nvar, nx, ny, nz), temperature(nx, ny, nz))
  mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.55643_dp]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  call require(ok, "CFL composition conversion")
  primitive(1:5) = [0.31_dp, 11.0_dp, -7.0_dp, 3.0_dp, 135000.0_dp]
  do species_index = 1, size(species)
    primitive(reactive_mass_fraction_component(species_index)) = &
      mass_fractions(species_index)
  end do
  call reactive_primitive_to_conserved( &
    species, primitive, conserved, recovered_temperature, sound_speed, ok)
  call require(ok, "CFL reference state")
  do component = 1, nvar
    state(component, :, :, :) = conserved(component)
  end do
  temperature = recovered_temperature

  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 2.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.5_dp, &
    "x", 0.7375_dp, geometry, ok)
  call require(ok, "CFL axis-plane geometry")
  call compute_reactive_eb_cfl_timestep_3d( &
    species, state, temperature, geometry, 0.5_dp, dt, ok)
  call require(ok, "CFL timestep")
  expected_rate = (11.0_dp + sound_speed) / geometry%dx + &
    (7.0_dp + sound_speed) / geometry%dy + &
    (3.0_dp + sound_speed) / geometry%dz
  expected_dt = 0.5_dp / expected_rate
  call assert_close(dt, expected_dt, 2.0e-15_dp * expected_dt, &
    "full-grid CFL analytical value")

  state(irho, 1, 1, 1) = -1.0_dp
  call compute_reactive_eb_cfl_timestep_3d( &
    species, state, temperature, geometry, 0.5_dp, dt, ok)
  call require(ok .and. dt == expected_dt, &
    "covered state excluded from CFL")
  state(:, 4, 1, 1) = conserved
  state(irho, 4, 1, 1) = -1.0_dp
  call compute_reactive_eb_cfl_timestep_3d( &
    species, state, temperature, geometry, 0.5_dp, dt, ok)
  call require(.not. ok .and. dt == 0.0_dp, &
    "invalid active state CFL transaction")

  state(:, 4, 1, 1) = conserved
  call compute_reactive_eb_cfl_timestep_3d( &
    species, state, temperature, geometry, 0.0_dp, dt, ok)
  call require(.not. ok .and. dt == 0.0_dp, &
    "invalid CFL number transaction")
  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 2.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.5_dp, &
    "x", 2.1_dp, geometry, ok)
  call require(ok, "fully covered CFL geometry")
  call compute_reactive_eb_cfl_timestep_3d( &
    species, state, temperature, geometry, 0.5_dp, dt, ok)
  call require(.not. ok .and. dt == 0.0_dp, &
    "fully covered CFL transaction")

  write(*, '(a)') "test_reactive_eb_cfl_3d: PASS"

contains

  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message

    if (abs(actual - expected) > tolerance) error stop message
  end subroutine assert_close

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_reactive_eb_cfl_3d
