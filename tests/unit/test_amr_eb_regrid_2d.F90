program test_amr_eb_regrid_2d
  use, intrinsic :: ieee_arithmetic, only: &
    ieee_is_finite, ieee_value, ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved
  use eb_geometry_2d_mod, only: &
    eb_geometry_2d, eb_covered_cell, build_eb_geometry_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, build_amr_eb_patch_2d, composite_eb_integral_2d
  use amr_eb_reactive_2d_mod, only: prolong_reactive_eb_patch_pcm_2d
  use amr_eb_regrid_2d_mod, only: &
    amr_eb_tagging_criteria_2d, amr_eb_regrid_plan_2d, &
    amr_eb_regrid_plan_collection_2d, &
    build_amr_eb_regrid_plan_collection_2d, &
    plan_reactive_eb_temperature_regrid_2d, &
    collapse_two_level_reactive_eb_patch_2d, &
    regrid_two_level_reactive_eb_patch_2d
  implicit none

  integer, parameter :: coarse_nx = 10, coarse_ny = 10, ratio = 2
  integer, parameter :: old_i_lower = 3, old_i_upper = 7
  integer, parameter :: old_j_lower = 3, old_j_upper = 7
  integer, parameter :: old_fine_nx = &
    (old_i_upper - old_i_lower + 1) * ratio
  integer, parameter :: old_fine_ny = &
    (old_j_upper - old_j_lower + 1) * ratio
  type(eb_geometry_2d) :: coarse_geometry
  type(eb_geometry_2d) :: old_fine_geometry, new_fine_geometry
  type(amr_eb_patch_2d) :: old_patch, new_patch
  type(amr_eb_tagging_criteria_2d) :: criteria
  type(amr_eb_tagging_criteria_2d) :: multipatch_criteria
  type(amr_eb_regrid_plan_2d) :: plan
  type(amr_eb_regrid_plan_collection_2d) :: plan_collection
  type(nasa7_species), allocatable :: species(:)
  real(dp) :: coarse_level_set(0:coarse_nx, 0:coarse_ny)
  real(dp), allocatable :: primitive(:), mass_fractions(:), state_cell(:)
  real(dp), allocatable :: coarse_state(:, :, :), old_fine_state(:, :, :)
  real(dp), allocatable :: new_coarse_state(:, :, :), new_fine_state(:, :, :)
  real(dp), allocatable :: collapsed_state(:, :, :)
  real(dp), allocatable :: coarse_temperature(:, :)
  real(dp), allocatable :: old_fine_temperature(:, :)
  real(dp), allocatable :: invalid_old_fine_temperature(:, :)
  real(dp), allocatable :: new_coarse_temperature(:, :)
  real(dp), allocatable :: new_fine_temperature(:, :)
  real(dp), allocatable :: collapsed_temperature(:, :)
  real(dp), allocatable :: integral_before(:), integral_after(:)
  real(dp) :: mole_fractions(7), x, y, temperature_cell, sound_speed
  real(dp) :: state_scale, integral_scale
  logical :: tags(coarse_nx, coarse_ny), ok, found_new_cell
  logical :: multipatch_tags(18, 12)
  integer :: component, i, j, k, nvar, new_fine_nx, new_fine_ny
  integer :: global_i, global_j, old_i, old_j, parent_i, parent_j

  do j = 0, coarse_ny
    y = real(j, dp) / real(coarse_ny, dp)
    do i = 0, coarse_nx
      x = real(i, dp) / real(coarse_nx, dp)
      coarse_level_set(i, j) = x + y - 0.30_dp
    end do
  end do
  call build_eb_geometry_2d( &
    coarse_level_set, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    coarse_geometry, ok)
  call require(ok, "coarse EB geometry")

  allocate(coarse_temperature(coarse_nx, coarse_ny))
  coarse_temperature = 1000.0_dp
  where (coarse_geometry%cell_type == eb_covered_cell)
    coarse_temperature = 3000.0_dp
  end where
  coarse_temperature(7, 6) = 2000.0_dp
  criteria%relative_gradient_threshold = 0.20_dp
  criteria%absolute_gradient_threshold = 100.0_dp
  criteria%scale_floor = 1.0_dp
  criteria%buffer_cells = 1
  criteria%minimum_patch_cells_x = 5
  criteria%minimum_patch_cells_y = 5
  call plan_reactive_eb_temperature_regrid_2d( &
    coarse_temperature, coarse_geometry, criteria, tags, plan, ok)
  call require(ok .and. plan%is_valid(), "temperature-gradient regrid plan")
  call require(plan%tagged_cell_count == 5, "five-point hotspot tagging")
  call require(.not. any(tags(1, :)) .and. &
    .not. any(tags(coarse_nx, :)) .and. .not. any(tags(:, 1)) .and. &
    .not. any(tags(:, coarse_ny)), "root-boundary cells remain untagged")
  call require(.not. any(tags .and. &
    coarse_geometry%cell_type == eb_covered_cell), &
    "covered cells remain untagged")
  call require(plan%tag_i_lower == 6 .and. plan%tag_i_upper == 8 .and. &
    plan%tag_j_lower == 5 .and. plan%tag_j_upper == 7, &
    "tag bounding box")
  call require(plan%coarse_i_lower == 5 .and. &
    plan%coarse_i_upper == 9 .and. plan%coarse_j_lower == 4 .and. &
    plan%coarse_j_upper == 8, "buffered internal patch plan")

  coarse_temperature = 1000.0_dp
  call plan_reactive_eb_temperature_regrid_2d( &
    coarse_temperature, coarse_geometry, criteria, tags, plan, ok)
  call require(ok .and. .not. plan%active .and. &
    plan%tagged_cell_count == 0, "empty-tag plan")
  coarse_temperature(7, 6) = 2000.0_dp
  call plan_reactive_eb_temperature_regrid_2d( &
    coarse_temperature, coarse_geometry, criteria, tags, plan, ok)
  call require(ok .and. plan%active, "restore active regrid plan")

  coarse_temperature = 1000.0_dp
  where (coarse_geometry%cell_type == eb_covered_cell)
    coarse_temperature = 3000.0_dp
  end where
  coarse_temperature(coarse_nx, 6) = 2000.0_dp
  call plan_reactive_eb_temperature_regrid_2d( &
    coarse_temperature, coarse_geometry, criteria, tags, plan, ok)
  call require(ok .and. plan%active .and. tags(coarse_nx, 6) .and. &
    plan%tag_i_upper == coarse_nx .and. &
    plan%coarse_i_upper == coarse_nx, &
    "one-sided physical-boundary temperature plan")
  coarse_temperature = 1000.0_dp
  where (coarse_geometry%cell_type == eb_covered_cell)
    coarse_temperature = 3000.0_dp
  end where
  coarse_temperature(7, 6) = 2000.0_dp
  call plan_reactive_eb_temperature_regrid_2d( &
    coarse_temperature, coarse_geometry, criteria, tags, plan, ok)
  call require(ok .and. plan%active .and. plan%coarse_i_upper == 9, &
    "restore internal topology plan")

  multipatch_criteria%buffer_cells = 1
  multipatch_criteria%minimum_patch_cells_x = 2
  multipatch_criteria%minimum_patch_cells_y = 2
  multipatch_criteria%maximum_patch_gap_cells = 0
  multipatch_tags = .false.
  multipatch_tags(3:4, 3:4) = .true.
  multipatch_tags(14:15, 8:9) = .true.
  call build_amr_eb_regrid_plan_collection_2d( &
    multipatch_tags, multipatch_criteria, plan_collection, ok)
  call require(ok .and. plan_collection%is_valid() .and. &
    plan_collection%patch_count() == 2 .and. &
    plan_collection%tagged_cell_count == 8, &
    "two disconnected EB tag clusters")
  call require( &
    plan_collection%plans(1)%coarse_i_lower == 2 .and. &
    plan_collection%plans(1)%coarse_i_upper == 5 .and. &
    plan_collection%plans(1)%coarse_j_lower == 2 .and. &
    plan_collection%plans(1)%coarse_j_upper == 5 .and. &
    plan_collection%plans(2)%coarse_i_lower == 13 .and. &
    plan_collection%plans(2)%coarse_i_upper == 16 .and. &
    plan_collection%plans(2)%coarse_j_lower == 7 .and. &
    plan_collection%plans(2)%coarse_j_upper == 10, &
    "deterministic buffered EB patch collection")

  multipatch_tags = .false.
  call build_amr_eb_regrid_plan_collection_2d( &
    multipatch_tags, multipatch_criteria, plan_collection, ok)
  call require(ok .and. plan_collection%is_valid() .and. &
    plan_collection%patch_count() == 0, "empty EB patch collection")

  multipatch_criteria%buffer_cells = 1
  multipatch_criteria%minimum_patch_cells_x = 1
  multipatch_criteria%minimum_patch_cells_y = 1
  multipatch_tags = .false.
  multipatch_tags(4, 4) = .true.
  multipatch_tags(8, 4) = .true.
  call build_amr_eb_regrid_plan_collection_2d( &
    multipatch_tags, multipatch_criteria, plan_collection, ok)
  call require(ok .and. plan_collection%patch_count() == 1 .and. &
    plan_collection%plans(1)%tagged_cell_count == 2 .and. &
    plan_collection%plans(1)%coarse_i_lower == 3 .and. &
    plan_collection%plans(1)%coarse_i_upper == 9 .and. &
    plan_collection%plans(1)%coarse_j_lower == 3 .and. &
    plan_collection%plans(1)%coarse_j_upper == 5, &
    "redistribution-adjacent EB plans coalesce")

  multipatch_criteria%buffer_cells = 0
  multipatch_criteria%maximum_patch_gap_cells = 1
  multipatch_tags = .false.
  multipatch_tags(4, 4) = .true.
  multipatch_tags(6, 4) = .true.
  call build_amr_eb_regrid_plan_collection_2d( &
    multipatch_tags, multipatch_criteria, plan_collection, ok)
  call require(ok .and. plan_collection%patch_count() == 1 .and. &
    plan_collection%plans(1)%tag_i_lower == 4 .and. &
    plan_collection%plans(1)%tag_i_upper == 6, &
    "configured EB tag gap joins one component")

  multipatch_tags = .false.
  multipatch_tags(1, 4) = .true.
  call build_amr_eb_regrid_plan_collection_2d( &
    multipatch_tags, multipatch_criteria, plan_collection, ok)
  call require(ok .and. plan_collection%patch_count() == 1 .and. &
    plan_collection%plans(1)%tag_i_lower == 1 .and. &
    plan_collection%plans(1)%coarse_i_lower == 1, &
    "physical-boundary EB tag collection")

  call build_patch_geometry( &
    coarse_geometry, old_i_lower, old_i_upper, old_j_lower, old_j_upper, &
    ratio, old_fine_geometry, old_patch, ok)
  call require(ok, "old fine geometry and patch")
  call build_patch_geometry( &
    coarse_geometry, plan%coarse_i_lower, plan%coarse_i_upper, &
    plan%coarse_j_lower, plan%coarse_j_upper, ratio, new_fine_geometry, &
    new_patch, ok)
  call require(ok, "new fine geometry and patch")
  new_fine_nx = new_fine_geometry%nx
  new_fine_ny = new_fine_geometry%ny

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamic database load")
  nvar = reactive_nvar(size(species))
  allocate(primitive(reactive_nprim(size(species))))
  allocate(mass_fractions(size(species)), state_cell(nvar))
  mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.55643_dp]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  call require(ok, "composition conversion")
  primitive(1:5) = [0.31_dp, 2.0_dp, -1.0_dp, 0.0_dp, 135000.0_dp]
  do k = 1, size(species)
    primitive(reactive_mass_fraction_component(k)) = mass_fractions(k)
  end do
  call reactive_primitive_to_conserved( &
    species, primitive, state_cell, temperature_cell, sound_speed, ok)
  call require(ok, "reference reactive state")

  allocate(coarse_state(nvar, coarse_nx, coarse_ny))
  allocate(old_fine_state(nvar, old_fine_nx, old_fine_ny))
  allocate(old_fine_temperature(old_fine_nx, old_fine_ny))
  allocate(new_coarse_state(nvar, coarse_nx, coarse_ny))
  allocate(new_coarse_temperature(coarse_nx, coarse_ny))
  allocate(new_fine_state(nvar, new_fine_nx, new_fine_ny))
  allocate(new_fine_temperature(new_fine_nx, new_fine_ny))
  allocate(integral_before(nvar), integral_after(nvar))
  coarse_state = spread(spread(state_cell, 2, coarse_nx), 3, coarse_ny)
  coarse_temperature = temperature_cell
  call prolong_reactive_eb_patch_pcm_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, &
    old_fine_geometry, old_patch, old_fine_state, old_fine_temperature, ok)
  call require(ok, "old-patch PCM initialization")
  old_fine_state = 1.01_dp * old_fine_state

  call composite_eb_integral_2d( &
    coarse_state, coarse_geometry, old_fine_state, old_fine_geometry, &
    old_patch, integral_before, ok)
  call require(ok, "old hierarchy composite integral")
  call regrid_two_level_reactive_eb_patch_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, &
    old_fine_state, old_fine_temperature, old_fine_geometry, old_patch, &
    new_fine_geometry, new_patch, new_coarse_state, &
    new_coarse_temperature, new_fine_state, new_fine_temperature, ok)
  call require(ok, "transactional EB hierarchy regrid")
  call composite_eb_integral_2d( &
    new_coarse_state, coarse_geometry, new_fine_state, &
    new_fine_geometry, new_patch, integral_after, ok)
  integral_scale = max(1.0_dp, maxval(abs(integral_before)))
  call require(ok .and. maxval(abs(integral_after - integral_before)) <= &
    5.0e-12_dp * integral_scale, "regrid composite conservation")

  state_scale = max(1.0_dp, maxval(abs(state_cell)))
  found_new_cell = .false.
  do j = 1, new_fine_ny
    global_j = (new_patch%coarse_j_lower - 1) * ratio + j
    old_j = global_j - (old_patch%coarse_j_lower - 1) * ratio
    parent_j = new_patch%coarse_j_lower + (j - 1) / ratio
    do i = 1, new_fine_nx
      global_i = (new_patch%coarse_i_lower - 1) * ratio + i
      old_i = global_i - (old_patch%coarse_i_lower - 1) * ratio
      parent_i = new_patch%coarse_i_lower + (i - 1) / ratio
      if (old_i >= 1 .and. old_i <= old_fine_nx .and. &
          old_j >= 1 .and. old_j <= old_fine_ny) then
        call require(maxval(abs(new_fine_state(:, i, j) - &
          old_fine_state(:, old_i, old_j))) == 0.0_dp, &
          "overlapping fine state retained")
      else if (new_fine_geometry%cell_type(i, j) /= eb_covered_cell) then
        found_new_cell = .true.
        call require(maxval(abs(new_fine_state(:, i, j) - &
          new_coarse_state(:, parent_i, parent_j))) <= &
          5.0e-14_dp * state_scale, "new fine cell uses coarse PCM")
      end if
    end do
  end do
  call require(found_new_cell, "newly refined active cells checked")
  call require(maxval(abs(new_coarse_state(:, 3:4, 3:7) - &
    1.01_dp * spread(spread(state_cell, 2, 2), 3, 5))) <= &
    5.0e-14_dp * state_scale, "retired fine region averaged down")

  allocate(collapsed_state, mold=new_coarse_state)
  allocate(collapsed_temperature, mold=new_coarse_temperature)
  call collapse_two_level_reactive_eb_patch_2d( &
    species, new_coarse_state, new_coarse_temperature, coarse_geometry, &
    new_fine_state, new_fine_geometry, new_patch, collapsed_state, &
    collapsed_temperature, ok)
  call require(ok, "reactive EB fine-patch collapse")
  do component = 1, nvar
    integral_after(component) = sum(coarse_geometry%volume_fraction * &
      collapsed_state(component, :, :)) * coarse_geometry%dx * &
      coarse_geometry%dy
  end do
  call require(maxval(abs(integral_after - integral_before)) <= &
    5.0e-12_dp * integral_scale, "fine-patch collapse conservation")
  call require(all(collapsed_temperature > 0.0_dp) .and. &
    all(ieee_is_finite(collapsed_temperature)), &
    "fine-patch collapse temperature recovery")

  old_fine_state(:, 5, 5) = &
    ieee_value(0.0_dp, ieee_quiet_nan)
  call regrid_two_level_reactive_eb_patch_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, &
    old_fine_state, old_fine_temperature, old_fine_geometry, old_patch, &
    new_fine_geometry, new_patch, new_coarse_state, &
    new_coarse_temperature, new_fine_state, new_fine_temperature, ok)
  call require(.not. ok .and. &
    maxval(abs(new_coarse_state - coarse_state)) == 0.0_dp .and. &
    maxval(abs(new_coarse_temperature - coarse_temperature)) == 0.0_dp .and. &
    maxval(abs(new_fine_state)) == 0.0_dp .and. &
    maxval(abs(new_fine_temperature)) == 0.0_dp, &
    "invalid old hierarchy regrid rollback")

  old_fine_state(:, 5, 5) = 1.01_dp * state_cell
  allocate(invalid_old_fine_temperature(1, 1), source=temperature_cell)
  call regrid_two_level_reactive_eb_patch_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, &
    old_fine_state, invalid_old_fine_temperature, old_fine_geometry, &
    old_patch, new_fine_geometry, new_patch, new_coarse_state, &
    new_coarse_temperature, new_fine_state, new_fine_temperature, ok)
  call require(.not. ok .and. &
    maxval(abs(new_coarse_state - coarse_state)) == 0.0_dp .and. &
    maxval(abs(new_coarse_temperature - coarse_temperature)) == 0.0_dp .and. &
    maxval(abs(new_fine_state)) == 0.0_dp .and. &
    maxval(abs(new_fine_temperature)) == 0.0_dp, &
    "invalid old-temperature shape rollback")

  write(*, '(a)') "test_amr_eb_regrid_2d: PASS"

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
    integer :: fine_cell_count_x, fine_cell_count_y, local_i, local_j

    fine_cell_count_x = (i_upper - i_lower + 1) * refinement_ratio
    fine_cell_count_y = (j_upper - j_lower + 1) * refinement_ratio
    x_lower = root_geometry%x_lower + real(i_lower - 1, dp) * &
      root_geometry%dx
    x_upper = root_geometry%x_lower + real(i_upper, dp) * root_geometry%dx
    y_lower = root_geometry%y_lower + real(j_lower - 1, dp) * &
      root_geometry%dy
    y_upper = root_geometry%y_lower + real(j_upper, dp) * root_geometry%dy
    allocate(level_set(0:fine_cell_count_x, 0:fine_cell_count_y))
    do local_j = 0, fine_cell_count_y
      local_y = y_lower + real(local_j, dp) * &
        (y_upper - y_lower) / real(fine_cell_count_y, dp)
      do local_i = 0, fine_cell_count_x
        local_x = x_lower + real(local_i, dp) * &
          (x_upper - x_lower) / real(fine_cell_count_x, dp)
        level_set(local_i, local_j) = local_x + local_y - 0.30_dp
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

end program test_amr_eb_regrid_2d
