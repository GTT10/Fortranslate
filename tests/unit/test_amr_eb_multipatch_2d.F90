program test_amr_eb_multipatch_2d
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved
  use eb_geometry_2d_mod, only: eb_geometry_2d, build_eb_geometry_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, build_amr_eb_patch_2d
  use amr_eb_regrid_2d_mod, only: &
    amr_eb_tagging_criteria_2d, amr_eb_regrid_plan_collection_2d, &
    reactive_eb_patch_set_2d, build_amr_eb_regrid_plan_collection_2d, &
    initialize_reactive_eb_patch_set_2d, &
    average_down_reactive_eb_patch_set_2d, &
    composite_reactive_eb_patch_set_integral_2d, &
    regrid_reactive_eb_patch_set_2d, &
    advance_reactive_eb_patch_set_hydro_2d
  implicit none

  integer, parameter :: coarse_nx = 10, coarse_ny = 10, ratio = 2
  type(eb_geometry_2d) :: coarse_geometry
  type(eb_geometry_2d), allocatable :: old_geometries(:)
  type(eb_geometry_2d), allocatable :: new_geometries(:)
  type(amr_eb_patch_2d) :: geometry_patch
  type(amr_eb_tagging_criteria_2d) :: criteria
  type(amr_eb_regrid_plan_collection_2d) :: old_collection
  type(amr_eb_regrid_plan_collection_2d) :: new_collection
  type(reactive_eb_patch_set_2d) :: old_set, new_set, removed_set
  type(reactive_eb_patch_set_2d) :: hydro_set, failed_set
  type(nasa7_species), allocatable :: species(:)
  real(dp) :: coarse_level_set(0:coarse_nx, 0:coarse_ny)
  real(dp), allocatable :: primitive(:), mass_fractions(:), state_cell(:)
  real(dp), allocatable :: coarse_state(:, :, :), coarse_temperature(:, :)
  real(dp), allocatable :: averaged_state(:, :, :)
  real(dp), allocatable :: averaged_temperature(:, :)
  real(dp), allocatable :: new_coarse_state(:, :, :)
  real(dp), allocatable :: new_coarse_temperature(:, :)
  real(dp), allocatable :: integral_before(:), integral_after(:)
  real(dp) :: mole_fractions(7), x, y, temperature_cell, sound_speed
  real(dp) :: dt, factor, integral_scale, state_scale
  character(len=64) :: hydro_failure_context
  logical :: old_tags(coarse_nx, coarse_ny)
  logical :: new_tags(coarse_nx, coarse_ny), ok
  integer :: child, i, j, nvar

  do j = 0, coarse_ny
    y = real(j, dp) / real(coarse_ny, dp)
    do i = 0, coarse_nx
      x = real(i, dp) / real(coarse_nx, dp)
      coarse_level_set(i, j) = x + y - 0.78_dp
    end do
  end do
  call build_eb_geometry_2d( &
    coarse_level_set, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    coarse_geometry, ok)
  call require(ok, "multipatch coarse EB geometry")

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "multipatch thermodynamic database")
  nvar = reactive_nvar(size(species))
  allocate(primitive(reactive_nprim(size(species))))
  allocate(mass_fractions(size(species)), state_cell(nvar))
  mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.55643_dp]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  call require(ok, "multipatch composition conversion")
  primitive(1:5) = [0.31_dp, 2.0_dp, -1.0_dp, 0.0_dp, 135000.0_dp]
  do i = 1, size(species)
    primitive(reactive_mass_fraction_component(i)) = mass_fractions(i)
  end do
  call reactive_primitive_to_conserved( &
    species, primitive, state_cell, temperature_cell, sound_speed, ok)
  call require(ok, "multipatch reference state")

  allocate(coarse_state(nvar, coarse_nx, coarse_ny))
  allocate(coarse_temperature(coarse_nx, coarse_ny))
  allocate(averaged_state(nvar, coarse_nx, coarse_ny))
  allocate(averaged_temperature(coarse_nx, coarse_ny))
  allocate(new_coarse_state(nvar, coarse_nx, coarse_ny))
  allocate(new_coarse_temperature(coarse_nx, coarse_ny))
  allocate(integral_before(nvar), integral_after(nvar))
  coarse_state = spread(spread(state_cell, 2, coarse_nx), 3, coarse_ny)
  coarse_temperature = temperature_cell

  criteria%buffer_cells = 0
  criteria%minimum_patch_cells_x = 2
  criteria%minimum_patch_cells_y = 2
  criteria%maximum_patch_gap_cells = 0
  old_tags = .false.
  old_tags(2:3, 4:5) = .true.
  old_tags(7:8, 7:8) = .true.
  call build_amr_eb_regrid_plan_collection_2d( &
    old_tags, criteria, old_collection, ok)
  call require(ok .and. old_collection%patch_count() == 2, &
    "old two-patch plan")
  allocate(old_geometries(old_collection%patch_count()))
  do child = 1, old_collection%patch_count()
    call build_patch_geometry( &
      coarse_geometry, old_collection%plans(child)%coarse_i_lower, &
      old_collection%plans(child)%coarse_i_upper, &
      old_collection%plans(child)%coarse_j_lower, &
      old_collection%plans(child)%coarse_j_upper, ratio, &
      old_geometries(child), geometry_patch, ok)
    call require(ok, "old child EB geometry")
  end do
  call initialize_reactive_eb_patch_set_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, &
    old_geometries, old_collection, ratio, old_set, ok)
  call require(ok .and. old_set%is_valid(coarse_geometry, nvar) .and. &
    old_set%patch_count() == 2, "initialize reactive EB patch set")

  call composite_reactive_eb_patch_set_integral_2d( &
    coarse_state, coarse_geometry, old_set, integral_before, ok)
  call require(ok, "initial multipatch hydro integral")
  dt = 0.02_dp * min(coarse_geometry%dx, coarse_geometry%dy) / sound_speed
  call advance_reactive_eb_patch_set_hydro_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, old_set, &
    "hllc", "pcm", "mc", 2, dt, new_coarse_state, &
    new_coarse_temperature, hydro_set, ok, &
    failure_context=hydro_failure_context)
  if (.not. ok) write(*, '(a)') trim(hydro_failure_context)
  call require(ok .and. hydro_set%is_valid(coarse_geometry, nvar), &
    "subcycled multipatch EB hydro")
  call composite_reactive_eb_patch_set_integral_2d( &
    new_coarse_state, coarse_geometry, hydro_set, integral_after, ok)
  integral_scale = max(1.0_dp, maxval(abs(integral_before)))
  call require(ok .and. maxval(abs(integral_after - integral_before)) <= &
    3.0e-12_dp * integral_scale, "multipatch hydro conservation")
  state_scale = max(1.0_dp, maxval(abs(state_cell)))
  call require(maxval(abs(new_coarse_state - coarse_state)) <= &
    4.0e-12_dp * state_scale, "uniform multipatch coarse preservation")
  do child = 1, old_set%patch_count()
    call require(maxval(abs(hydro_set%children(child)%state - &
      old_set%children(child)%state)) <= 4.0e-12_dp * state_scale, &
      "uniform multipatch fine preservation")
  end do
  call average_down_reactive_eb_patch_set_2d( &
    species, new_coarse_state, new_coarse_temperature, coarse_geometry, &
    hydro_set, averaged_state, averaged_temperature, ok)
  call require(ok .and. maxval(abs(averaged_state - &
    new_coarse_state)) <= 4.0e-13_dp * state_scale, &
    "multipatch hydro synchronized hierarchy")

  call advance_reactive_eb_patch_set_hydro_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, old_set, &
    "invalid", "pcm", "mc", 2, dt, new_coarse_state, &
    new_coarse_temperature, failed_set, ok)
  call require(.not. ok .and. &
    maxval(abs(new_coarse_state - coarse_state)) == 0.0_dp .and. &
    maxval(abs(new_coarse_temperature - coarse_temperature)) == 0.0_dp .and. &
    failed_set%patch_count() == old_set%patch_count(), &
    "multipatch hydro coarse rollback")
  do child = 1, old_set%patch_count()
    call require(maxval(abs(failed_set%children(child)%state - &
      old_set%children(child)%state)) == 0.0_dp .and. &
      maxval(abs(failed_set%children(child)%temperature - &
        old_set%children(child)%temperature)) == 0.0_dp, &
      "multipatch hydro fine rollback")
  end do

  do child = 1, old_set%patch_count()
    do j = 1, old_set%children(child)%geometry%ny
      do i = 1, old_set%children(child)%geometry%nx
        factor = 1.0_dp + 1.0e-3_dp * &
          real(10 * child + i + 2 * j, dp)
        old_set%children(child)%state(:, i, j) = factor * &
          old_set%children(child)%state(:, i, j)
      end do
    end do
  end do
  call require(old_set%is_valid(coarse_geometry, nvar), &
    "perturbed reactive EB patch set")
  call composite_reactive_eb_patch_set_integral_2d( &
    coarse_state, coarse_geometry, old_set, integral_before, ok)
  call require(ok, "old multipatch composite integral")
  call average_down_reactive_eb_patch_set_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, old_set, &
    averaged_state, averaged_temperature, ok)
  call require(ok, "multipatch reactive average down")
  integral_after = 0.0_dp
  do i = 1, nvar
    integral_after(i) = sum(coarse_geometry%volume_fraction * &
      averaged_state(i, :, :)) * coarse_geometry%dx * coarse_geometry%dy
  end do
  integral_scale = max(1.0_dp, maxval(abs(integral_before)))
  call require(maxval(abs(integral_after - integral_before)) <= &
    8.0e-12_dp * integral_scale, "multipatch average-down conservation")

  new_tags = .false.
  new_tags(3:4, 4:5) = .true.
  new_tags(7:8, 7:8) = .true.
  call build_amr_eb_regrid_plan_collection_2d( &
    new_tags, criteria, new_collection, ok)
  call require(ok .and. new_collection%patch_count() == 2, &
    "moved two-patch plan")
  allocate(new_geometries(new_collection%patch_count()))
  do child = 1, new_collection%patch_count()
    call build_patch_geometry( &
      coarse_geometry, new_collection%plans(child)%coarse_i_lower, &
      new_collection%plans(child)%coarse_i_upper, &
      new_collection%plans(child)%coarse_j_lower, &
      new_collection%plans(child)%coarse_j_upper, ratio, &
      new_geometries(child), geometry_patch, ok)
    call require(ok, "new child EB geometry")
  end do
  call regrid_reactive_eb_patch_set_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, old_set, &
    new_geometries, new_collection, ratio, new_coarse_state, &
    new_coarse_temperature, new_set, ok)
  call require(ok .and. new_set%is_valid(coarse_geometry, nvar) .and. &
    new_set%patch_count() == 2, "transactional multipatch regrid")
  call composite_reactive_eb_patch_set_integral_2d( &
    new_coarse_state, coarse_geometry, new_set, integral_after, ok)
  call require(ok .and. maxval(abs(integral_after - integral_before)) <= &
    8.0e-12_dp * integral_scale, "multipatch regrid conservation")
  call require(maxval(abs(new_set%children(1)%state(:, 1:2, :) - &
    old_set%children(1)%state(:, 3:4, :))) == 0.0_dp, &
    "moved patch exact fine overlap")
  call require(maxval(abs(new_set%children(2)%state - &
    old_set%children(2)%state)) == 0.0_dp, &
    "unchanged patch exact state retention")
  do j = 1, new_set%children(1)%geometry%ny
    call require(maxval(abs(new_set%children(1)%state(:, 3:4, j) - &
      spread(new_coarse_state(:, 4, 4 + (j - 1) / ratio), 2, 2))) <= &
      5.0e-14_dp * state_scale, &
      "new patch cells use synchronized coarse PCM")
  end do

  old_set%children(1)%state(:, 1, 1) = &
    ieee_value(0.0_dp, ieee_quiet_nan)
  call regrid_reactive_eb_patch_set_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, old_set, &
    new_geometries, new_collection, ratio, new_coarse_state, &
    new_coarse_temperature, removed_set, ok)
  call require(.not. ok .and. .not. allocated(removed_set%children) .and. &
    maxval(abs(new_coarse_state - coarse_state)) == 0.0_dp .and. &
    maxval(abs(new_coarse_temperature - coarse_temperature)) == 0.0_dp, &
    "invalid multipatch regrid rollback")
  old_set%children(1)%state(:, 1, 1) = &
    (1.0_dp + 1.0e-3_dp * 13.0_dp) * state_cell

  new_tags = .false.
  call build_amr_eb_regrid_plan_collection_2d( &
    new_tags, criteria, new_collection, ok)
  call require(ok .and. new_collection%patch_count() == 0, &
    "empty multipatch removal plan")
  deallocate(new_geometries)
  allocate(new_geometries(0))
  call regrid_reactive_eb_patch_set_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, old_set, &
    new_geometries, new_collection, ratio, new_coarse_state, &
    new_coarse_temperature, removed_set, ok)
  call require(ok .and. removed_set%is_valid(coarse_geometry, nvar) .and. &
    removed_set%patch_count() == 0, "remove complete EB patch set")
  integral_after = 0.0_dp
  do i = 1, nvar
    integral_after(i) = sum(coarse_geometry%volume_fraction * &
      new_coarse_state(i, :, :)) * coarse_geometry%dx * coarse_geometry%dy
  end do
  call require(maxval(abs(integral_after - integral_before)) <= &
    8.0e-12_dp * integral_scale, "multipatch removal conservation")

  write(*, '(a)') "test_amr_eb_multipatch_2d: PASS"

contains

  subroutine build_patch_geometry( &
      root_geometry, i_lower, i_upper, j_lower, j_upper, refinement_ratio, &
      fine_geometry, patch, valid)
    type(eb_geometry_2d), intent(in) :: root_geometry
    integer, intent(in) :: i_lower, i_upper, j_lower, j_upper
    integer, intent(in) :: refinement_ratio
    type(eb_geometry_2d), intent(out) :: fine_geometry
    type(amr_eb_patch_2d), intent(out) :: patch
    logical, intent(out) :: valid

    real(dp), allocatable :: level_set(:, :)
    real(dp) :: x_lower, x_upper, y_lower, y_upper, local_x, local_y
    integer :: fine_nx, fine_ny, local_i, local_j

    fine_nx = (i_upper - i_lower + 1) * refinement_ratio
    fine_ny = (j_upper - j_lower + 1) * refinement_ratio
    x_lower = root_geometry%x_lower + real(i_lower - 1, dp) * &
      root_geometry%dx
    x_upper = root_geometry%x_lower + real(i_upper, dp) * root_geometry%dx
    y_lower = root_geometry%y_lower + real(j_lower - 1, dp) * &
      root_geometry%dy
    y_upper = root_geometry%y_lower + real(j_upper, dp) * root_geometry%dy
    allocate(level_set(0:fine_nx, 0:fine_ny))
    do local_j = 0, fine_ny
      local_y = y_lower + real(local_j, dp) * &
        (y_upper - y_lower) / real(fine_ny, dp)
      do local_i = 0, fine_nx
        local_x = x_lower + real(local_i, dp) * &
          (x_upper - x_lower) / real(fine_nx, dp)
        level_set(local_i, local_j) = local_x + local_y - 0.78_dp
      end do
    end do
    call build_eb_geometry_2d( &
      level_set, x_lower, x_upper, y_lower, y_upper, fine_geometry, valid)
    if (.not. valid) return
    call build_amr_eb_patch_2d( &
      root_geometry, fine_geometry, i_lower, i_upper, j_lower, j_upper, &
      refinement_ratio, patch, valid)
  end subroutine build_patch_geometry

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_amr_eb_multipatch_2d
