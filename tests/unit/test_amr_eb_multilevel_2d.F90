program test_amr_eb_multilevel_2d
  use, intrinsic :: ieee_arithmetic, only: &
    ieee_is_finite, ieee_value, ieee_quiet_nan
  use precision_mod, only: dp
  use state_indices_mod, only: irho, iet
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_mass_fraction_component, reactive_primitive_to_conserved
  use eb_geometry_2d_mod, only: eb_geometry_2d, build_eb_geometry_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, build_amr_eb_patch_2d
  use amr_eb_multilevel_2d_mod, only: &
    average_down_three_level_eb_state_2d, &
    average_down_three_level_reactive_eb_state_2d, &
    composite_three_level_eb_integral_2d
  use amr_eb_multilevel_reactive_2d_mod, only: &
    advance_three_level_reactive_eb_hydro_2d
  implicit none

  integer, parameter :: root_nx = 8, root_ny = 8, ratio = 2
  integer, parameter :: root_i_lower = 2, root_i_upper = 7
  integer, parameter :: root_j_lower = 2, root_j_upper = 7
  integer, parameter :: level_one_nx = &
    (root_i_upper - root_i_lower + 1) * ratio
  integer, parameter :: level_one_ny = &
    (root_j_upper - root_j_lower + 1) * ratio
  integer, parameter :: level_one_i_lower = 6, level_one_i_upper = 10
  integer, parameter :: level_one_j_lower = 6, level_one_j_upper = 10
  integer, parameter :: level_two_nx = &
    (level_one_i_upper - level_one_i_lower + 1) * ratio
  integer, parameter :: level_two_ny = &
    (level_one_j_upper - level_one_j_lower + 1) * ratio
  type(eb_geometry_2d) :: root_geometry
  type(eb_geometry_2d) :: level_one_geometry, level_two_geometry
  type(amr_eb_patch_2d) :: root_patch, level_one_patch
  type(nasa7_species), allocatable :: species(:)
  real(dp) :: root_level_set(0:root_nx, 0:root_ny)
  real(dp), allocatable :: root_state(:, :, :)
  real(dp), allocatable :: level_one_state(:, :, :)
  real(dp), allocatable :: level_two_state(:, :, :)
  real(dp), allocatable :: synchronized_root(:, :, :)
  real(dp), allocatable :: synchronized_level_one(:, :, :)
  real(dp), allocatable :: integral_before(:), integral_after(:)
  real(dp), allocatable :: primitive(:), mass_fractions(:), state_cell(:)
  real(dp), allocatable :: reactive_root(:, :, :)
  real(dp), allocatable :: reactive_level_one(:, :, :)
  real(dp), allocatable :: reactive_level_two(:, :, :)
  real(dp), allocatable :: reactive_root_sync(:, :, :)
  real(dp), allocatable :: reactive_level_one_sync(:, :, :)
  real(dp), allocatable :: reactive_level_two_sync(:, :, :)
  real(dp), allocatable :: root_temperature(:, :)
  real(dp), allocatable :: level_one_temperature(:, :)
  real(dp), allocatable :: level_two_temperature(:, :)
  real(dp), allocatable :: root_temperature_sync(:, :)
  real(dp), allocatable :: level_one_temperature_sync(:, :)
  real(dp), allocatable :: level_two_temperature_sync(:, :)
  real(dp) :: mole_fractions(7), x, y, temperature_cell, sound_speed
  real(dp) :: scale, dt
  logical :: ok
  integer :: i, j, k, nvar

  do j = 0, root_ny
    y = real(j, dp) / real(root_ny, dp)
    do i = 0, root_nx
      x = real(i, dp) / real(root_nx, dp)
      root_level_set(i, j) = x + y - 0.78_dp
    end do
  end do
  call build_eb_geometry_2d( &
    root_level_set, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, root_geometry, ok)
  call require(ok, "three-level root geometry")
  call build_patch_geometry( &
    root_geometry, root_i_lower, root_i_upper, root_j_lower, root_j_upper, &
    ratio, level_one_geometry, root_patch, ok)
  call require(ok, "three-level middle geometry")
  call build_patch_geometry( &
    level_one_geometry, level_one_i_lower, level_one_i_upper, &
    level_one_j_lower, level_one_j_upper, ratio, level_two_geometry, &
    level_one_patch, ok)
  call require(ok, "three-level finest geometry")

  allocate(root_state(1, root_nx, root_ny), source=1.0_dp)
  allocate(level_one_state(1, level_one_nx, level_one_ny), source=2.0_dp)
  allocate(level_two_state(1, level_two_nx, level_two_ny), source=3.0_dp)
  allocate(synchronized_root(1, root_nx, root_ny))
  allocate(synchronized_level_one(1, level_one_nx, level_one_ny))
  allocate(integral_before(1), integral_after(1))
  call composite_three_level_eb_integral_2d( &
    root_state, root_geometry, level_one_state, level_one_geometry, &
    root_patch, level_two_state, level_two_geometry, level_one_patch, &
    integral_before, ok)
  call require(ok, "three-level composite integral")
  call average_down_three_level_eb_state_2d( &
    root_state, root_geometry, level_one_state, level_one_geometry, &
    root_patch, level_two_state, level_two_geometry, level_one_patch, &
    synchronized_root, synchronized_level_one, ok)
  call require(ok, "deepest-to-root EB synchronization")
  integral_after(1) = sum(root_geometry%volume_fraction * &
    synchronized_root(1, :, :)) * root_geometry%dx * root_geometry%dy
  call assert_close(integral_after(1), integral_before(1), 8.0e-13_dp, &
    "three-level synchronization conservation")
  call require(maxval(abs(synchronized_root(:, 1, :) - &
    root_state(:, 1, :))) == 0.0_dp .and. &
    maxval(abs(synchronized_root(:, 8, :) - &
      root_state(:, 8, :))) == 0.0_dp, &
    "root cells outside middle patch unchanged")
  call require(maxval(abs(synchronized_level_one(:, &
    level_one_i_lower:level_one_i_upper, &
    level_one_j_lower:level_one_j_upper) - 3.0_dp)) == 0.0_dp, &
    "finest constant restricted into middle")

  level_two_state(1, level_two_nx, level_two_ny) = &
    ieee_value(0.0_dp, ieee_quiet_nan)
  call average_down_three_level_eb_state_2d( &
    root_state, root_geometry, level_one_state, level_one_geometry, &
    root_patch, level_two_state, level_two_geometry, level_one_patch, &
    synchronized_root, synchronized_level_one, ok)
  call require(.not. ok .and. all(synchronized_root == root_state) .and. &
    all(synchronized_level_one == level_one_state), &
    "three-level nonfinite rollback")
  level_two_state(1, level_two_nx, level_two_ny) = 3.0_dp

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "three-level thermodynamic database")
  nvar = reactive_nvar(size(species))
  allocate(primitive(reactive_nprim(size(species))))
  allocate(mass_fractions(size(species)), state_cell(nvar))
  mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.55643_dp]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  call require(ok, "three-level composition conversion")
  primitive(1:5) = [0.31_dp, 2.0_dp, -1.0_dp, 0.0_dp, 135000.0_dp]
  do i = 1, size(species)
    primitive(reactive_mass_fraction_component(i)) = mass_fractions(i)
  end do
  call reactive_primitive_to_conserved( &
    species, primitive, state_cell, temperature_cell, sound_speed, ok)
  call require(ok, "three-level reference reactive state")
  allocate(reactive_root(nvar, root_nx, root_ny))
  allocate(reactive_level_one(nvar, level_one_nx, level_one_ny))
  allocate(reactive_level_two(nvar, level_two_nx, level_two_ny))
  allocate(reactive_root_sync(nvar, root_nx, root_ny))
  allocate(reactive_level_one_sync(nvar, level_one_nx, level_one_ny))
  allocate(reactive_level_two_sync(nvar, level_two_nx, level_two_ny))
  allocate(root_temperature(root_nx, root_ny), source=temperature_cell)
  allocate(level_one_temperature( &
    level_one_nx, level_one_ny), source=temperature_cell)
  allocate(level_two_temperature( &
    level_two_nx, level_two_ny), source=temperature_cell)
  allocate(root_temperature_sync(root_nx, root_ny))
  allocate(level_one_temperature_sync(level_one_nx, level_one_ny))
  allocate(level_two_temperature_sync(level_two_nx, level_two_ny))
  reactive_root = spread(spread(state_cell, 2, root_nx), 3, root_ny)
  reactive_level_one = 1.01_dp * &
    spread(spread(state_cell, 2, level_one_nx), 3, level_one_ny)
  reactive_level_two = 0.99_dp * &
    spread(spread(state_cell, 2, level_two_nx), 3, level_two_ny)
  deallocate(integral_before, integral_after)
  allocate(integral_before(nvar), integral_after(nvar))
  call composite_three_level_eb_integral_2d( &
    reactive_root, root_geometry, reactive_level_one, level_one_geometry, &
    root_patch, reactive_level_two, level_two_geometry, level_one_patch, &
    integral_before, ok)
  call require(ok, "three-level reactive composite integral")
  call average_down_three_level_reactive_eb_state_2d( &
    species, reactive_root, root_temperature, root_geometry, &
    reactive_level_one, level_one_temperature, level_one_geometry, &
    root_patch, reactive_level_two, level_two_temperature, &
    level_two_geometry, level_one_patch, reactive_root_sync, &
    root_temperature_sync, reactive_level_one_sync, &
    level_one_temperature_sync, ok)
  call require(ok, "three-level reactive synchronization")
  do i = 1, nvar
    integral_after(i) = sum(root_geometry%volume_fraction * &
      reactive_root_sync(i, :, :)) * root_geometry%dx * root_geometry%dy
  end do
  scale = max(1.0_dp, maxval(abs(integral_before)))
  call require(maxval(abs(integral_after - integral_before)) <= &
    8.0e-12_dp * scale, "three-level reactive conservation")
  call require(all(ieee_is_finite(root_temperature_sync)) .and. &
    all(ieee_is_finite(level_one_temperature_sync)) .and. &
    all(root_temperature_sync > 0.0_dp) .and. &
    all(level_one_temperature_sync > 0.0_dp), &
    "three-level reactive thermodynamics")
  call require(all(reactive_root_sync(:, 1, :) == reactive_root(:, 1, :)) &
    .and. all(reactive_level_one_sync(:, 1:2, :) == &
      reactive_level_one(:, 1:2, :)), &
    "three-level reactive unrefined ownership")

  level_two_temperature(1, 1) = -1.0_dp
  call average_down_three_level_reactive_eb_state_2d( &
    species, reactive_root, root_temperature, root_geometry, &
    reactive_level_one, level_one_temperature, level_one_geometry, &
    root_patch, reactive_level_two, level_two_temperature, &
    level_two_geometry, level_one_patch, reactive_root_sync, &
    root_temperature_sync, reactive_level_one_sync, &
    level_one_temperature_sync, ok)
  call require(.not. ok .and. all(reactive_root_sync == reactive_root) .and. &
    all(root_temperature_sync == root_temperature) .and. &
    all(reactive_level_one_sync == reactive_level_one) .and. &
    all(level_one_temperature_sync == level_one_temperature), &
    "three-level reactive rollback")

  level_two_temperature(1, 1) = temperature_cell
  primitive(2:4) = 0.0_dp
  call reactive_primitive_to_conserved( &
    species, primitive, state_cell, temperature_cell, sound_speed, ok)
  call require(ok, "stationary three-level hydro state")
  call report_hydro_residual(1.0_dp, 1.0_dp, "uniform")
  call report_hydro_residual(1.0_dp, 0.99_dp, "inner mismatch")
  call report_hydro_residual(1.01_dp, 1.01_dp, "outer mismatch")
  reactive_root = spread(spread(state_cell, 2, root_nx), 3, root_ny)
  reactive_level_one = 1.01_dp * &
    spread(spread(state_cell, 2, level_one_nx), 3, level_one_ny)
  reactive_level_two = 0.99_dp * &
    spread(spread(state_cell, 2, level_two_nx), 3, level_two_ny)
  root_temperature = temperature_cell
  level_one_temperature = temperature_cell
  level_two_temperature = temperature_cell
  call composite_three_level_eb_integral_2d( &
    reactive_root, root_geometry, reactive_level_one, level_one_geometry, &
    root_patch, reactive_level_two, level_two_geometry, level_one_patch, &
    integral_before, ok)
  call require(ok, "three-level hydro initial composite integral")
  dt = 0.015_dp * min(root_geometry%dx, root_geometry%dy) / sound_speed
  call advance_three_level_reactive_eb_hydro_2d( &
    species, reactive_root, root_temperature, root_geometry, &
    reactive_level_one, level_one_temperature, level_one_geometry, &
    root_patch, reactive_level_two, level_two_temperature, &
    level_two_geometry, level_one_patch, "hllc", "pcm", "mc", 2, dt, &
    reactive_root_sync, root_temperature_sync, reactive_level_one_sync, &
    level_one_temperature_sync, reactive_level_two_sync, &
    level_two_temperature_sync, ok)
  call require(ok, "recursive three-level reactive EB hydro")
  call composite_three_level_eb_integral_2d( &
    reactive_root_sync, root_geometry, reactive_level_one_sync, &
    level_one_geometry, root_patch, reactive_level_two_sync, &
    level_two_geometry, level_one_patch, integral_after, ok)
  scale = max(1.0_dp, maxval(abs(integral_before)))
  write(*, '(a,5(es24.16,1x))') &
    "three-level hydro mass diagnostic: ", integral_before(irho), &
    integral_after(irho), integral_after(irho) - integral_before(irho), &
    integral_before(iet), integral_after(iet) - integral_before(iet)
  call require(ok .and. &
    abs(integral_after(irho) - integral_before(irho)) <= &
      8.0e-10_dp * scale .and. &
    abs(integral_after(iet) - integral_before(iet)) <= &
      8.0e-10_dp * scale, &
    "three-level hydro mass and energy conservation")
  do k = 1, size(species)
    call require(abs(integral_after(reactive_species_component(k)) - &
      integral_before(reactive_species_component(k))) <= &
      8.0e-10_dp * scale, "three-level hydro species conservation")
  end do
  call require(all(ieee_is_finite(root_temperature_sync)) .and. &
    all(ieee_is_finite(level_one_temperature_sync)) .and. &
    all(ieee_is_finite(level_two_temperature_sync)) .and. &
    all(root_temperature_sync > 0.0_dp) .and. &
    all(level_one_temperature_sync > 0.0_dp) .and. &
    all(level_two_temperature_sync > 0.0_dp), &
    "three-level hydro thermodynamics")

  call average_down_three_level_reactive_eb_state_2d( &
    species, reactive_root_sync, root_temperature_sync, root_geometry, &
    reactive_level_one_sync, level_one_temperature_sync, &
    level_one_geometry, root_patch, reactive_level_two_sync, &
    level_two_temperature_sync, level_two_geometry, level_one_patch, &
    reactive_root, root_temperature, reactive_level_one, &
    level_one_temperature, ok)
  call require(ok .and. maxval(abs(reactive_root - &
    reactive_root_sync)) <= 2.0e-12_dp * scale .and. &
    maxval(abs(reactive_level_one - reactive_level_one_sync)) <= &
      2.0e-12_dp * scale, "three-level hydro final synchronization")

  reactive_root = spread(spread(state_cell, 2, root_nx), 3, root_ny)
  reactive_level_one = 1.01_dp * &
    spread(spread(state_cell, 2, level_one_nx), 3, level_one_ny)
  reactive_level_two = 0.99_dp * &
    spread(spread(state_cell, 2, level_two_nx), 3, level_two_ny)
  root_temperature = temperature_cell
  level_one_temperature = temperature_cell
  level_two_temperature = temperature_cell
  call advance_three_level_reactive_eb_hydro_2d( &
    species, reactive_root, root_temperature, root_geometry, &
    reactive_level_one, level_one_temperature, level_one_geometry, &
    root_patch, reactive_level_two, level_two_temperature, &
    level_two_geometry, level_one_patch, "unknown", "pcm", "mc", 2, dt, &
    reactive_root_sync, root_temperature_sync, reactive_level_one_sync, &
    level_one_temperature_sync, reactive_level_two_sync, &
    level_two_temperature_sync, ok)
  call require(.not. ok .and. all(reactive_root_sync == reactive_root) .and. &
    all(root_temperature_sync == root_temperature) .and. &
    all(reactive_level_one_sync == reactive_level_one) .and. &
    all(level_one_temperature_sync == level_one_temperature) .and. &
    all(reactive_level_two_sync == reactive_level_two) .and. &
    all(level_two_temperature_sync == level_two_temperature), &
    "three-level hydro rollback")

  write(*, '(a)') "test_amr_eb_multilevel_2d: PASS"

contains

  subroutine report_hydro_residual( &
      level_one_scale, level_two_scale, label)
    real(dp), intent(in) :: level_one_scale, level_two_scale
    character(len=*), intent(in) :: label

    reactive_root = spread(spread(state_cell, 2, root_nx), 3, root_ny)
    reactive_level_one = level_one_scale * &
      spread(spread(state_cell, 2, level_one_nx), 3, level_one_ny)
    reactive_level_two = level_two_scale * &
      spread(spread(state_cell, 2, level_two_nx), 3, level_two_ny)
    root_temperature = temperature_cell
    level_one_temperature = temperature_cell
    level_two_temperature = temperature_cell
    call composite_three_level_eb_integral_2d( &
      reactive_root, root_geometry, reactive_level_one, level_one_geometry, &
      root_patch, reactive_level_two, level_two_geometry, level_one_patch, &
      integral_before, ok)
    call require(ok, "diagnostic initial integral")
    dt = 0.015_dp * min(root_geometry%dx, root_geometry%dy) / sound_speed
    call advance_three_level_reactive_eb_hydro_2d( &
      species, reactive_root, root_temperature, root_geometry, &
      reactive_level_one, level_one_temperature, level_one_geometry, &
      root_patch, reactive_level_two, level_two_temperature, &
      level_two_geometry, level_one_patch, "hllc", "pcm", "mc", 2, dt, &
      reactive_root_sync, root_temperature_sync, reactive_level_one_sync, &
      level_one_temperature_sync, reactive_level_two_sync, &
      level_two_temperature_sync, ok)
    call require(ok, "diagnostic three-level hydro")
    call composite_three_level_eb_integral_2d( &
      reactive_root_sync, root_geometry, reactive_level_one_sync, &
      level_one_geometry, root_patch, reactive_level_two_sync, &
      level_two_geometry, level_one_patch, integral_after, ok)
    call require(ok, "diagnostic final integral")
    write(*, '(a,1x,a,2(es24.16,1x))') &
      "three-level hydro residual", trim(label), &
      integral_after(irho) - integral_before(irho), &
      integral_after(iet) - integral_before(iet)
  end subroutine report_hydro_residual

  subroutine build_patch_geometry( &
      parent_geometry, i_lower, i_upper, j_lower, j_upper, &
      refinement_ratio, child_geometry, patch, valid)
    type(eb_geometry_2d), intent(in) :: parent_geometry
    integer, intent(in) :: i_lower, i_upper, j_lower, j_upper
    integer, intent(in) :: refinement_ratio
    type(eb_geometry_2d), intent(out) :: child_geometry
    type(amr_eb_patch_2d), intent(out) :: patch
    logical, intent(out) :: valid

    real(dp), allocatable :: level_set(:, :)
    real(dp) :: x_lower, x_upper, y_lower, y_upper, local_x, local_y
    integer :: nx, ny, local_i, local_j

    nx = (i_upper - i_lower + 1) * refinement_ratio
    ny = (j_upper - j_lower + 1) * refinement_ratio
    x_lower = parent_geometry%x_lower + real(i_lower - 1, dp) * &
      parent_geometry%dx
    x_upper = parent_geometry%x_lower + real(i_upper, dp) * &
      parent_geometry%dx
    y_lower = parent_geometry%y_lower + real(j_lower - 1, dp) * &
      parent_geometry%dy
    y_upper = parent_geometry%y_lower + real(j_upper, dp) * &
      parent_geometry%dy
    allocate(level_set(0:nx, 0:ny))
    do local_j = 0, ny
      local_y = y_lower + real(local_j, dp) * &
        (y_upper - y_lower) / real(ny, dp)
      do local_i = 0, nx
        local_x = x_lower + real(local_i, dp) * &
          (x_upper - x_lower) / real(nx, dp)
        level_set(local_i, local_j) = local_x + local_y - 0.78_dp
      end do
    end do
    call build_eb_geometry_2d( &
      level_set, x_lower, x_upper, y_lower, y_upper, child_geometry, valid)
    if (.not. valid) return
    call build_amr_eb_patch_2d( &
      parent_geometry, child_geometry, i_lower, i_upper, j_lower, j_upper, &
      refinement_ratio, patch, valid)
  end subroutine build_patch_geometry

  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message

    call require(abs(actual - expected) <= &
      tolerance * max(1.0_dp, abs(expected)), message)
  end subroutine assert_close

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_amr_eb_multilevel_2d
