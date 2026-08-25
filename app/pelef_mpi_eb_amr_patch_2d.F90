program pelef_mpi_eb_amr_patch_2d
  use, intrinsic :: iso_fortran_env, only: int64
  use mpi_f08
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
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
    initialize_reactive_eb_patch_set_2d
  use mpi_amr_eb_patch_2d_mod, only: &
    mpi_amr_eb_patch_distribution_2d, &
    initialize_mpi_amr_eb_patch_distribution_2d, &
    mpi_amr_eb_distribution_matches_patch_set_2d, &
    synchronize_owned_reactive_eb_patch_set_2d
  implicit none

  integer, parameter :: coarse_nx = 14, coarse_ny = 14, ratio = 2
  type(eb_geometry_2d) :: coarse_geometry
  type(eb_geometry_2d), allocatable :: fine_geometries(:)
  type(amr_eb_patch_2d) :: geometry_patch
  type(amr_eb_tagging_criteria_2d) :: criteria
  type(amr_eb_regrid_plan_collection_2d) :: collection
  type(reactive_eb_patch_set_2d) :: patch_set, local_patch_set
  type(reactive_eb_patch_set_2d) :: synchronized_patch_set
  type(mpi_amr_eb_patch_distribution_2d) :: distribution
  type(mpi_amr_eb_patch_distribution_2d) :: invalid_distribution
  type(mpi_amr_eb_patch_distribution_2d) :: rejected_distribution
  type(nasa7_species), allocatable :: species(:)
  real(dp) :: coarse_level_set(0:coarse_nx, 0:coarse_ny)
  real(dp), allocatable :: primitive(:), mass_fractions(:), state_cell(:)
  real(dp), allocatable :: coarse_state(:, :, :), coarse_temperature(:, :)
  real(dp), allocatable :: local_coarse_state(:, :, :)
  real(dp), allocatable :: local_coarse_temperature(:, :)
  real(dp), allocatable :: synchronized_coarse_state(:, :, :)
  real(dp), allocatable :: synchronized_coarse_temperature(:, :)
  real(dp) :: mole_fractions(7), x, y, temperature, sound_speed, factor
  logical :: tags(coarse_nx, coarse_ny), ok
  integer :: child, i, ierr, j, nvar, rank, nranks, tile
  integer :: inconsistent_exponent

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Init failed"
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Comm_rank failed"
  call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Comm_size failed"

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
  call assert_all(ok, "MPI EB AMR coarse geometry", rank)

  call load_h2o2_elementary_thermo(species, ok)
  call assert_all(ok, "MPI EB AMR thermodynamic database", rank)
  nvar = reactive_nvar(size(species))
  allocate(primitive(reactive_nprim(size(species))))
  allocate(mass_fractions(size(species)), state_cell(nvar))
  mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.55643_dp]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  call assert_all(ok, "MPI EB AMR composition", rank)
  primitive(1:5) = [0.31_dp, 0.0_dp, 0.0_dp, 0.0_dp, 135000.0_dp]
  do i = 1, size(species)
    primitive(reactive_mass_fraction_component(i)) = mass_fractions(i)
  end do
  call reactive_primitive_to_conserved( &
    species, primitive, state_cell, temperature, sound_speed, ok)
  call assert_all(ok, "MPI EB AMR reference state", rank)
  allocate(coarse_state(nvar, coarse_nx, coarse_ny))
  allocate(coarse_temperature(coarse_nx, coarse_ny))
  coarse_state = spread(spread(state_cell, 2, coarse_nx), 3, coarse_ny)
  coarse_temperature = temperature

  criteria%buffer_cells = 0
  criteria%minimum_patch_cells_x = 5
  criteria%minimum_patch_cells_y = 5
  criteria%maximum_patch_gap_cells = 0
  tags = .false.
  tags(1:5, 6:10) = .true.
  tags(9:13, 9:13) = .true.
  call build_amr_eb_regrid_plan_collection_2d( &
    tags, criteria, collection, ok)
  call assert_all(ok .and. collection%patch_count() == 2, &
    "MPI EB AMR sibling plan", rank)
  allocate(fine_geometries(collection%patch_count()))
  do child = 1, collection%patch_count()
    call build_patch_geometry( &
      coarse_geometry, collection%plans(child)%coarse_i_lower, &
      collection%plans(child)%coarse_i_upper, &
      collection%plans(child)%coarse_j_lower, &
      collection%plans(child)%coarse_j_upper, ratio, &
      fine_geometries(child), geometry_patch, ok)
    call assert_all(ok, "MPI EB AMR child geometry", rank)
  end do
  call initialize_reactive_eb_patch_set_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, &
    fine_geometries, collection, ratio, patch_set, ok)
  call assert_all(ok .and. patch_set%patch_count() == 2, &
    "MPI EB AMR patch set", rank)

  call initialize_mpi_amr_eb_patch_distribution_2d( &
    coarse_geometry, patch_set, MPI_COMM_WORLD, distribution, ok, 2)
  call assert_all(ok .and. &
    mpi_amr_eb_distribution_matches_patch_set_2d( &
      distribution, coarse_geometry, patch_set), &
    "MPI EB AMR deterministic distribution", rank)
  call assert_all(distribution%root_tile_count() == min(coarse_ny, nranks), &
    "MPI EB AMR root tiling", rank)
  call assert_all(distribution%child_count() == 2 .and. &
    sum(distribution%rank_entity_counts) == &
      distribution%root_tile_count() + 2, &
    "MPI EB AMR entity ownership", rank)
  call assert_all(sum(distribution%rank_cell_counts) == &
    coarse_nx * coarse_ny + 2 * 10 * 10, &
    "MPI EB AMR owned cell accounting", rank)
  call assert_all(sum(distribution%rank_work_counts) == 996_int64, &
    "MPI EB AMR parabolic work accounting", rank)
  call assert_all(minval(distribution%rank_entity_counts) >= 1, &
    "MPI EB AMR every rank owns a root tile", rank)

  allocate(local_coarse_state, source=coarse_state)
  allocate(local_coarse_temperature, source=coarse_temperature)
  allocate(synchronized_coarse_state, mold=coarse_state)
  allocate(synchronized_coarse_temperature, mold=coarse_temperature)
  local_patch_set = patch_set
  do tile = 1, distribution%root_tile_count()
    factor = root_tile_factor(tile)
    if (distribution%root_tile_is_local(tile)) then
      local_coarse_state(:, :, &
        distribution%root_tiles(tile)%j_lower: &
        distribution%root_tiles(tile)%j_upper) = factor * &
        coarse_state(:, :, distribution%root_tiles(tile)%j_lower: &
          distribution%root_tiles(tile)%j_upper)
      local_coarse_temperature(:, &
        distribution%root_tiles(tile)%j_lower: &
        distribution%root_tiles(tile)%j_upper) = coarse_temperature(:, &
          distribution%root_tiles(tile)%j_lower: &
          distribution%root_tiles(tile)%j_upper)
    else
      local_coarse_state(:, :, &
        distribution%root_tiles(tile)%j_lower: &
        distribution%root_tiles(tile)%j_upper) = 0.75_dp * &
        coarse_state(:, :, distribution%root_tiles(tile)%j_lower: &
          distribution%root_tiles(tile)%j_upper)
      local_coarse_temperature(:, &
        distribution%root_tiles(tile)%j_lower: &
        distribution%root_tiles(tile)%j_upper) = 0.9_dp * &
          coarse_temperature(:, distribution%root_tiles(tile)%j_lower: &
            distribution%root_tiles(tile)%j_upper)
    end if
  end do
  do child = 1, distribution%child_count()
    factor = child_factor(child)
    if (distribution%child_is_local(child)) then
      local_patch_set%children(child)%state = &
        factor * patch_set%children(child)%state
      local_patch_set%children(child)%temperature = &
        patch_set%children(child)%temperature
    else
      local_patch_set%children(child)%state = &
        0.75_dp * patch_set%children(child)%state
      local_patch_set%children(child)%temperature = &
        0.9_dp * patch_set%children(child)%temperature
    end if
  end do
  call synchronize_owned_reactive_eb_patch_set_2d( &
    distribution, size(species), local_coarse_state, &
    local_coarse_temperature, coarse_geometry, local_patch_set, &
    synchronized_coarse_state, synchronized_coarse_temperature, &
    synchronized_patch_set, ok)
  call assert_all(ok, "MPI EB AMR owner-authoritative synchronization", rank)
  do tile = 1, distribution%root_tile_count()
    factor = root_tile_factor(tile)
    call assert_all(all(synchronized_coarse_state(:, :, &
      distribution%root_tiles(tile)%j_lower: &
      distribution%root_tiles(tile)%j_upper) == factor * &
      coarse_state(:, :, distribution%root_tiles(tile)%j_lower: &
        distribution%root_tiles(tile)%j_upper)), &
      "MPI EB AMR root owner payload", rank)
    call assert_all(all(synchronized_coarse_temperature(:, &
      distribution%root_tiles(tile)%j_lower: &
      distribution%root_tiles(tile)%j_upper) == coarse_temperature(:, &
      distribution%root_tiles(tile)%j_lower: &
        distribution%root_tiles(tile)%j_upper)), &
      "MPI EB AMR root temperature payload", rank)
  end do
  do child = 1, distribution%child_count()
    factor = child_factor(child)
    call assert_all(all(synchronized_patch_set%children(child)%state == &
      factor * patch_set%children(child)%state), &
      "MPI EB AMR child owner payload", rank)
    call assert_all(all( &
      synchronized_patch_set%children(child)%temperature == &
      patch_set%children(child)%temperature), &
      "MPI EB AMR child temperature payload", rank)
  end do

  invalid_distribution = distribution
  invalid_distribution%root_tiles(1)%owner = nranks
  call synchronize_owned_reactive_eb_patch_set_2d( &
    invalid_distribution, size(species), local_coarse_state, &
    local_coarse_temperature, coarse_geometry, local_patch_set, &
    synchronized_coarse_state, synchronized_coarse_temperature, &
    synchronized_patch_set, ok)
  call assert_all(.not. ok .and. &
    all(synchronized_coarse_state == local_coarse_state) .and. &
    all(synchronized_coarse_temperature == local_coarse_temperature), &
    "MPI EB AMR invalid ownership rollback", rank)

  inconsistent_exponent = 0
  if (nranks > 1) inconsistent_exponent = modulo(rank, 2)
  call initialize_mpi_amr_eb_patch_distribution_2d( &
    coarse_geometry, patch_set, MPI_COMM_WORLD, rejected_distribution, ok, &
    inconsistent_exponent)
  call assert_all(ok .eqv. (nranks == 1), &
    "MPI EB AMR inconsistent work model rejection", rank)
  call initialize_mpi_amr_eb_patch_distribution_2d( &
    coarse_geometry, patch_set, MPI_COMM_WORLD, rejected_distribution, ok, 3)
  call assert_all(.not. ok, "MPI EB AMR invalid work model rejection", rank)

  if (rank == 0) then
    write(*, '(a)') "PeleF " // pelef_version // &
      " MPI EB AMR ownership 2D: PASS"
    write(*, '(a,i0)') "MPI ranks: ", nranks
    write(*, '(a,i0)') "Root tiles: ", distribution%root_tile_count()
    write(*, '(a,i0)') "Fine sibling patches: ", distribution%child_count()
  end if
  call MPI_Finalize(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Finalize failed"

contains

  pure real(dp) function root_tile_factor(local_tile) result(value)
    integer, intent(in) :: local_tile

    value = 1.0_dp + 1.0e-3_dp * real(local_tile, dp)
  end function root_tile_factor

  pure real(dp) function child_factor(local_child) result(value)
    integer, intent(in) :: local_child

    value = 1.0_dp + 1.0e-2_dp * real(local_child, dp)
  end function child_factor

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

  subroutine assert_all(condition, message, local_rank)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    integer, intent(in) :: local_rank

    logical :: global_condition
    integer :: local_ierr

    call MPI_Allreduce(condition, global_condition, 1, MPI_LOGICAL, MPI_LAND, &
      MPI_COMM_WORLD, local_ierr)
    if (local_ierr == MPI_SUCCESS .and. global_condition) return
    if (local_rank == 0) write(*, '(a)') "FAIL: " // trim(message)
    call MPI_Abort(MPI_COMM_WORLD, 1, local_ierr)
  end subroutine assert_all

end program pelef_mpi_eb_amr_patch_2d
