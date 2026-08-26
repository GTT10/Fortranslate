program pelef_mpi_eb_amr_patch_2d
  use, intrinsic :: iso_fortran_env, only: int64
  use mpi_f08
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
  use state_indices_mod, only: irho
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport
  use reactive_boundary_2d_mod, only: &
    reactive_boundary_set_2d, initialize_periodic_boundary_set_2d
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_mass_fraction_component, reactive_primitive_to_conserved
  use reactive_2d_mod, only: advance_reactive_chemistry_2d
  use eb_geometry_2d_mod, only: &
    eb_geometry_2d, eb_covered_cell, build_eb_geometry_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, build_amr_eb_patch_2d
  use amr_eb_regrid_2d_mod, only: &
    amr_eb_tagging_criteria_2d, amr_eb_regrid_plan_collection_2d, &
    reactive_eb_patch_set_2d, build_amr_eb_regrid_plan_collection_2d, &
    initialize_reactive_eb_patch_set_2d, &
    average_down_reactive_eb_patch_set_2d, &
    advance_reactive_eb_patch_set_hydro_2d
  use amr_eb_multipatch_transport_2d_mod, only: &
    advance_reactive_eb_patch_set_transport_2d
  use eb_reactive_transport_2d_mod, only: &
    reactive_eb_transport_timestep_2d
  use reactive_eb_amr_2d_driver_mod, only: &
    advance_reactive_eb_patch_set_strang_2d, &
    compute_reactive_eb_patch_set_cfl_timestep_2d
  use mpi_amr_eb_patch_2d_mod, only: &
    mpi_amr_eb_patch_distribution_2d, &
    mpi_amr_eb_sparse_patch_set_2d, &
    initialize_mpi_amr_eb_patch_distribution_2d, &
    mpi_amr_eb_distribution_matches_patch_set_2d, &
    synchronize_owned_reactive_eb_patch_set_2d, &
    advance_owned_reactive_eb_patch_set_chemistry_2d, &
    advance_owned_reactive_eb_patch_set_hydro_2d, &
    advance_owned_reactive_eb_patch_set_transport_2d, &
    advance_owned_reactive_eb_patch_set_strang_2d, &
    scatter_owned_reactive_eb_patch_set_2d, &
    materialize_owned_reactive_eb_patch_set_2d, &
    average_down_sparse_owned_reactive_eb_patch_set_2d, &
    advance_sparse_owned_reactive_eb_patch_set_chemistry_2d, &
    compute_sparse_owned_reactive_eb_patch_set_timestep_2d, &
    advance_sparse_owned_reactive_eb_patch_set_hydro_2d, &
    advance_sparse_owned_reactive_eb_patch_set_transport_2d, &
    advance_sparse_owned_reactive_eb_patch_set_strang_2d
  implicit none

  integer, parameter :: coarse_nx = 14, coarse_ny = 14, ratio = 2
  type(eb_geometry_2d) :: coarse_geometry
  type(eb_geometry_2d), allocatable :: fine_geometries(:)
  type(amr_eb_patch_2d) :: geometry_patch
  type(amr_eb_tagging_criteria_2d) :: criteria
  type(amr_eb_regrid_plan_collection_2d) :: collection
  type(reactive_eb_patch_set_2d) :: patch_set, local_patch_set
  type(reactive_eb_patch_set_2d) :: synchronized_patch_set
  type(reactive_eb_patch_set_2d) :: materialized_patch_set
  type(reactive_eb_patch_set_2d) :: rejected_materialized_set
  type(mpi_amr_eb_patch_distribution_2d) :: distribution
  type(mpi_amr_eb_patch_distribution_2d) :: invalid_distribution
  type(mpi_amr_eb_patch_distribution_2d) :: rejected_distribution
  type(mpi_amr_eb_sparse_patch_set_2d) :: sparse_patch_set
  type(mpi_amr_eb_sparse_patch_set_2d) :: invalid_sparse_patch_set
  type(mpi_amr_eb_sparse_patch_set_2d) :: sparse_chemistry_set
  type(mpi_amr_eb_sparse_patch_set_2d) :: sparse_failed_set
  type(mpi_amr_eb_sparse_patch_set_2d) :: sparse_failed_backup_set
  type(mpi_amr_eb_sparse_patch_set_2d) :: sparse_hydro_set
  type(mpi_amr_eb_sparse_patch_set_2d) :: sparse_hydro_failed_set
  type(mpi_amr_eb_sparse_patch_set_2d) :: sparse_hydro_failed_backup_set
  type(mpi_amr_eb_sparse_patch_set_2d) :: sparse_transport_set
  type(mpi_amr_eb_sparse_patch_set_2d) :: sparse_timestep_set
  type(mpi_amr_eb_sparse_patch_set_2d) :: sparse_timestep_failed_set
  type(mpi_amr_eb_sparse_patch_set_2d) :: sparse_timestep_failed_backup_set
  type(mpi_amr_eb_sparse_patch_set_2d) :: sparse_transport_failed_set
  type(mpi_amr_eb_sparse_patch_set_2d) :: sparse_transport_failed_backup_set
  type(mpi_amr_eb_sparse_patch_set_2d) :: sparse_full_set
  type(mpi_amr_eb_sparse_patch_set_2d) :: sparse_full_failed_set
  type(mpi_amr_eb_sparse_patch_set_2d) :: sparse_full_failed_backup_set
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(reactive_boundary_set_2d) :: boundaries
  type(reactive_eb_patch_set_2d) :: chemistry_set, reference_set
  type(reactive_eb_patch_set_2d) :: failed_set, failed_backup_set
  type(reactive_eb_patch_set_2d) :: hydro_start_set, hydro_reference_set
  type(reactive_eb_patch_set_2d) :: hydro_mpi_set, hydro_failed_set
  type(reactive_eb_patch_set_2d) :: hydro_failed_backup_set
  type(reactive_eb_patch_set_2d) :: transport_start_set
  type(reactive_eb_patch_set_2d) :: transport_reference_set
  type(reactive_eb_patch_set_2d) :: transport_mpi_set
  type(reactive_eb_patch_set_2d) :: transport_failed_set
  type(reactive_eb_patch_set_2d) :: transport_failed_backup_set
  type(reactive_eb_patch_set_2d) :: full_reference_set, full_mpi_set
  type(reactive_eb_patch_set_2d) :: full_failed_set, full_failed_backup_set
  real(dp) :: coarse_level_set(0:coarse_nx, 0:coarse_ny)
  real(dp), allocatable :: primitive(:), mass_fractions(:), state_cell(:)
  real(dp), allocatable :: coarse_state(:, :, :), coarse_temperature(:, :)
  real(dp), allocatable :: local_coarse_state(:, :, :)
  real(dp), allocatable :: local_coarse_temperature(:, :)
  real(dp), allocatable :: synchronized_coarse_state(:, :, :)
  real(dp), allocatable :: synchronized_coarse_temperature(:, :)
  real(dp), allocatable :: materialized_coarse_state(:, :, :)
  real(dp), allocatable :: materialized_coarse_temperature(:, :)
  real(dp), allocatable :: rejected_materialized_state(:, :, :)
  real(dp), allocatable :: rejected_materialized_temperature(:, :)
  real(dp), allocatable :: chemistry_coarse_state(:, :, :)
  real(dp), allocatable :: chemistry_coarse_temperature(:, :)
  real(dp), allocatable :: reference_coarse_state(:, :, :)
  real(dp), allocatable :: reference_coarse_temperature(:, :)
  real(dp), allocatable :: averaged_coarse_state(:, :, :)
  real(dp), allocatable :: averaged_coarse_temperature(:, :)
  real(dp), allocatable :: failed_coarse_state(:, :, :)
  real(dp), allocatable :: failed_coarse_temperature(:, :)
  real(dp), allocatable :: failed_backup_state(:, :, :)
  real(dp), allocatable :: failed_backup_temperature(:, :)
  real(dp), allocatable :: hydro_reference_state(:, :, :)
  real(dp), allocatable :: hydro_reference_temperature(:, :)
  real(dp), allocatable :: hydro_mpi_state(:, :, :)
  real(dp), allocatable :: hydro_mpi_temperature(:, :)
  real(dp), allocatable :: hydro_failed_state(:, :, :)
  real(dp), allocatable :: hydro_failed_temperature(:, :)
  real(dp), allocatable :: hydro_failed_backup_state(:, :, :)
  real(dp), allocatable :: hydro_failed_backup_temperature(:, :)
  real(dp), allocatable :: transport_reference_state(:, :, :)
  real(dp), allocatable :: transport_reference_temperature(:, :)
  real(dp), allocatable :: transport_mpi_state(:, :, :)
  real(dp), allocatable :: transport_mpi_temperature(:, :)
  real(dp), allocatable :: transport_failed_state(:, :, :)
  real(dp), allocatable :: transport_failed_temperature(:, :)
  real(dp), allocatable :: transport_failed_backup_state(:, :, :)
  real(dp), allocatable :: transport_failed_backup_temperature(:, :)
  real(dp), allocatable :: full_reference_state(:, :, :)
  real(dp), allocatable :: full_reference_temperature(:, :)
  real(dp), allocatable :: full_mpi_state(:, :, :)
  real(dp), allocatable :: full_mpi_temperature(:, :)
  real(dp), allocatable :: full_failed_state(:, :, :)
  real(dp), allocatable :: full_failed_temperature(:, :)
  real(dp), allocatable :: full_failed_backup_state(:, :, :)
  real(dp), allocatable :: full_failed_backup_temperature(:, :)
  logical, allocatable :: active_mask(:, :)
  logical, allocatable :: restriction_recipients(:)
  real(dp) :: mole_fractions(7), x, y, temperature, sound_speed, factor
  real(dp) :: chemistry_interval, chemistry_change, state_scale
  real(dp) :: hydro_dt, hydro_change, hydro_scale
  real(dp) :: transport_dt, transport_change, transport_scale
  real(dp) :: transport_theta, transport_reference_theta
  real(dp) :: timestep_dt, timestep_reference_dt, timestep_entity_dt
  real(dp) :: timestep_maximum_diffusivity
  real(dp) :: child_sound_speed, child_temperature
  real(dp) :: full_reference_theta, full_theta, full_scale, full_change
  logical :: tags(coarse_nx, coarse_ny), ok
  integer :: child, component, global_advances, i, ierr, j, owner
  integer :: local_advances, nvar, rank, nranks, tile
  integer :: expected_global_advances, expected_local_advances
  integer :: local_chemistry_advances, local_hydro_advances
  integer :: local_transport_advances
  integer :: global_chemistry_advances, global_hydro_advances
  integer :: global_transport_advances
  integer :: expected_local_chemistry, expected_local_hydro
  integer :: expected_local_transport, expected_global_chemistry
  integer :: expected_global_hydro, expected_global_transport
  integer :: inconsistent_exponent
  integer :: sparse_local_values, sparse_global_values
  integer :: sparse_expected_local_values, sparse_expected_global_values
  integer :: local_restriction_transfers, global_restriction_transfers
  integer :: expected_local_restriction_transfers
  integer :: expected_global_restriction_transfers
  integer :: local_root_transfers, global_root_transfers
  integer :: expected_local_root_transfers
  integer :: expected_global_root_transfers, root_owner
  character(len=160) :: full_failure_context

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
  call load_h2o2_elementary_mechanism(reactions, ok)
  call assert_all(ok, "MPI EB AMR chemistry mechanism", rank)
  call load_h2o2_elementary_transport(transport, ok)
  call assert_all(ok, "MPI EB AMR transport database", rank)
  nvar = reactive_nvar(size(species))
  call initialize_periodic_boundary_set_2d( &
    reactive_nprim(size(species)), boundaries)
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

  call scatter_owned_reactive_eb_patch_set_2d( &
    distribution, size(species), local_coarse_state, &
    local_coarse_temperature, coarse_geometry, local_patch_set, &
    sparse_patch_set, ok)
  call assert_all(ok .and. sparse_patch_set%is_valid( &
    distribution, coarse_geometry, patch_set), &
    "MPI EB AMR sparse owner storage", rank)
  sparse_expected_local_values = 0
  do tile = 1, distribution%root_tile_count()
    call assert_all( &
      (allocated(sparse_patch_set%root_tiles(tile)%state) .eqv. &
        distribution%root_tile_is_local(tile)) .and. &
      (allocated(sparse_patch_set%root_tiles(tile)%temperature) .eqv. &
        distribution%root_tile_is_local(tile)), &
      "MPI EB AMR sparse root allocation", rank)
    if (distribution%root_tile_is_local(tile)) &
      sparse_expected_local_values = sparse_expected_local_values + &
        (nvar + 1) * distribution%root_tiles(tile)%cell_count
  end do
  do child = 1, distribution%child_count()
    call assert_all( &
      (allocated(sparse_patch_set%children(child)%state) .eqv. &
        distribution%child_is_local(child)) .and. &
      (allocated(sparse_patch_set%children(child)%temperature) .eqv. &
        distribution%child_is_local(child)), &
      "MPI EB AMR sparse child allocation", rank)
    if (distribution%child_is_local(child)) &
      sparse_expected_local_values = sparse_expected_local_values + &
        (nvar + 1) * distribution%child_cell_counts(child)
  end do
  sparse_local_values = int(sparse_patch_set%local_value_count())
  call MPI_Allreduce( &
    sparse_local_values, sparse_global_values, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  sparse_expected_global_values = (nvar + 1) * &
    (coarse_nx * coarse_ny + sum(distribution%child_cell_counts))
  call assert_all(ierr == MPI_SUCCESS .and. &
    sparse_local_values == sparse_expected_local_values .and. &
    sparse_global_values == sparse_expected_global_values, &
    "MPI EB AMR sparse value accounting", rank)

  allocate(materialized_coarse_state, mold=coarse_state)
  allocate(materialized_coarse_temperature, mold=coarse_temperature)
  call materialize_owned_reactive_eb_patch_set_2d( &
    distribution, sparse_patch_set, coarse_state, coarse_temperature, &
    coarse_geometry, patch_set, materialized_coarse_state, &
    materialized_coarse_temperature, materialized_patch_set, ok)
  call assert_all(ok .and. &
    all(materialized_coarse_state == synchronized_coarse_state) .and. &
    all(materialized_coarse_temperature == &
      synchronized_coarse_temperature), &
    "MPI EB AMR sparse root materialization", rank)
  do child = 1, distribution%child_count()
    call assert_all( &
      all(materialized_patch_set%children(child)%state == &
        synchronized_patch_set%children(child)%state) .and. &
      all(materialized_patch_set%children(child)%temperature == &
        synchronized_patch_set%children(child)%temperature), &
      "MPI EB AMR sparse child materialization", rank)
  end do

  invalid_sparse_patch_set = sparse_patch_set
  if (rank == 0) then
    do tile = 1, distribution%root_tile_count()
      if (.not. distribution%root_tile_is_local(tile)) cycle
      deallocate(invalid_sparse_patch_set%root_tiles(tile)%state)
      exit
    end do
  end if
  allocate(rejected_materialized_state, mold=coarse_state)
  allocate(rejected_materialized_temperature, mold=coarse_temperature)
  call materialize_owned_reactive_eb_patch_set_2d( &
    distribution, invalid_sparse_patch_set, coarse_state, &
    coarse_temperature, coarse_geometry, patch_set, &
    rejected_materialized_state, rejected_materialized_temperature, &
    rejected_materialized_set, ok)
  call assert_all(.not. ok .and. &
    all(rejected_materialized_state == coarse_state) .and. &
    all(rejected_materialized_temperature == coarse_temperature), &
    "MPI EB AMR invalid sparse root rollback", rank)
  do child = 1, distribution%child_count()
    call assert_all( &
      all(rejected_materialized_set%children(child)%state == &
        patch_set%children(child)%state) .and. &
      all(rejected_materialized_set%children(child)%temperature == &
        patch_set%children(child)%temperature), &
      "MPI EB AMR invalid sparse child rollback", rank)
  end do

  allocate(chemistry_coarse_state, source=coarse_state)
  allocate(chemistry_coarse_temperature, source=coarse_temperature)
  allocate(reference_coarse_state, source=coarse_state)
  allocate(reference_coarse_temperature, source=coarse_temperature)
  allocate(averaged_coarse_state, mold=coarse_state)
  allocate(averaged_coarse_temperature, mold=coarse_temperature)
  chemistry_set = patch_set
  reference_set = patch_set
  chemistry_interval = &
    0.02_dp * min(coarse_geometry%dx, coarse_geometry%dy) / sound_speed
  active_mask = coarse_geometry%cell_type /= eb_covered_cell
  call advance_reactive_chemistry_2d( &
    species, reactions, reference_coarse_state, &
    reference_coarse_temperature, coarse_geometry%nx, coarse_geometry%ny, &
    chemistry_interval, 1.0e-8_dp, 1.0e-14_dp, ok, active_mask)
  call assert_all(ok, "serial EB AMR root chemistry reference", rank)
  deallocate(active_mask)
  do child = 1, reference_set%patch_count()
    active_mask = reference_set%children(child)%geometry%cell_type /= &
      eb_covered_cell
    call advance_reactive_chemistry_2d( &
      species, reactions, reference_set%children(child)%state, &
      reference_set%children(child)%temperature, &
      reference_set%children(child)%geometry%nx, &
      reference_set%children(child)%geometry%ny, chemistry_interval, &
      1.0e-8_dp, 1.0e-14_dp, ok, active_mask)
    call assert_all(ok, "serial EB AMR child chemistry reference", rank)
    deallocate(active_mask)
  end do
  call average_down_reactive_eb_patch_set_2d( &
    species, reference_coarse_state, reference_coarse_temperature, &
    coarse_geometry, reference_set, averaged_coarse_state, &
    averaged_coarse_temperature, ok)
  call assert_all(ok, "serial EB AMR chemistry average down", rank)
  reference_coarse_state = averaged_coarse_state
  reference_coarse_temperature = averaged_coarse_temperature

  call scatter_owned_reactive_eb_patch_set_2d( &
    distribution, size(species), coarse_state, coarse_temperature, &
    coarse_geometry, patch_set, sparse_chemistry_set, ok)
  call assert_all(ok, "MPI EB AMR sparse chemistry scatter", rank)
  allocate(restriction_recipients(nranks))
  expected_local_restriction_transfers = 0
  expected_global_restriction_transfers = 0
  do child = 1, distribution%child_count()
    restriction_recipients = .false.
    do tile = 1, distribution%root_tile_count()
      if (distribution%root_tiles(tile)%j_upper < &
          patch_set%children(child)%patch%coarse_j_lower .or. &
          distribution%root_tiles(tile)%j_lower > &
          patch_set%children(child)%patch%coarse_j_upper) cycle
      owner = distribution%root_tiles(tile)%owner
      restriction_recipients(owner + 1) = .true.
    end do
    owner = distribution%child_owner(child)
    restriction_recipients(owner + 1) = .false.
    expected_global_restriction_transfers = &
      expected_global_restriction_transfers + &
      count(restriction_recipients)
    if (rank == owner) expected_local_restriction_transfers = &
      expected_local_restriction_transfers + &
      count(restriction_recipients)
  end do
  call advance_sparse_owned_reactive_eb_patch_set_chemistry_2d( &
    species, reactions, chemistry_interval, 1.0e-8_dp, 1.0e-14_dp, &
    distribution, sparse_chemistry_set, coarse_geometry, patch_set, ok, &
    local_advances, local_restriction_transfers)
  call assert_all(ok .and. sparse_chemistry_set%is_valid( &
    distribution, coarse_geometry, patch_set) .and. &
    int(sparse_chemistry_set%local_value_count()) == &
      sparse_expected_local_values, &
    "MPI EB AMR direct sparse chemistry", rank)
  call MPI_Allreduce( &
    local_advances, global_advances, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_advances == distribution%rank_entity_counts(rank + 1) .and. &
    global_advances == distribution%root_tile_count() + &
      distribution%child_count(), &
    "MPI EB AMR sparse chemistry owner accounting", rank)
  call MPI_Allreduce( &
    local_restriction_transfers, global_restriction_transfers, 1, &
    MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_restriction_transfers == &
      expected_local_restriction_transfers .and. &
    global_restriction_transfers == &
      expected_global_restriction_transfers, &
    "MPI EB AMR targeted sparse restriction accounting", rank)
  call materialize_owned_reactive_eb_patch_set_2d( &
    distribution, sparse_chemistry_set, coarse_state, coarse_temperature, &
    coarse_geometry, patch_set, materialized_coarse_state, &
    materialized_coarse_temperature, materialized_patch_set, ok)
  call assert_all(ok .and. &
    all(materialized_coarse_state == reference_coarse_state) .and. &
    all(materialized_coarse_temperature == reference_coarse_temperature), &
    "MPI EB AMR sparse chemistry serial root parity", rank)
  do child = 1, distribution%child_count()
    call assert_all( &
      all(materialized_patch_set%children(child)%state == &
        reference_set%children(child)%state) .and. &
      all(materialized_patch_set%children(child)%temperature == &
        reference_set%children(child)%temperature), &
      "MPI EB AMR sparse chemistry serial child parity", rank)
  end do

  call advance_owned_reactive_eb_patch_set_chemistry_2d( &
    species, reactions, chemistry_interval, 1.0e-8_dp, 1.0e-14_dp, &
    distribution, chemistry_coarse_state, chemistry_coarse_temperature, &
    coarse_geometry, chemistry_set, ok, local_advances)
  call assert_all(ok, "MPI owner-only EB AMR chemistry", rank)
  call MPI_Allreduce( &
    local_advances, global_advances, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_advances == distribution%rank_entity_counts(rank + 1) .and. &
    global_advances == distribution%root_tile_count() + &
      distribution%child_count(), &
    "MPI EB AMR one chemistry call per owner entity", rank)
  state_scale = max(1.0_dp, maxval(abs(reference_coarse_state)))
  call assert_all(maxval(abs(chemistry_coarse_state - &
    reference_coarse_state)) <= 5.0e-13_dp * state_scale .and. &
    maxval(abs(chemistry_coarse_temperature - &
      reference_coarse_temperature)) <= 5.0e-13_dp * &
        max(1.0_dp, maxval(abs(reference_coarse_temperature))), &
    "MPI EB AMR chemistry serial root parity", rank)
  do child = 1, chemistry_set%patch_count()
    call assert_all(maxval(abs(chemistry_set%children(child)%state - &
      reference_set%children(child)%state)) <= 5.0e-13_dp * state_scale .and. &
      maxval(abs(chemistry_set%children(child)%temperature - &
        reference_set%children(child)%temperature)) <= 5.0e-13_dp * &
          max(1.0_dp, maxval(abs( &
            reference_set%children(child)%temperature))), &
      "MPI EB AMR chemistry serial child parity", rank)
  end do
  chemistry_change = 0.0_dp
  do child = 1, chemistry_set%patch_count()
    do i = 1, size(species)
      component = reactive_species_component(i)
      chemistry_change = max(chemistry_change, maxval(abs( &
        chemistry_set%children(child)%state(component, :, :) - &
        patch_set%children(child)%state(component, :, :))))
    end do
  end do
  call assert_all(chemistry_change > 1.0e-14_dp * state_scale, &
    "MPI EB AMR chemistry changes species", rank)

  allocate(failed_coarse_state, source=coarse_state)
  allocate(failed_coarse_temperature, source=coarse_temperature)
  failed_set = patch_set
  tile = distribution%root_tile_count()
  if (distribution%root_tile_is_local(tile)) then
    failed_coarse_state(irho, :, &
      distribution%root_tiles(tile)%j_lower: &
        distribution%root_tiles(tile)%j_upper) = -1.0_dp
  end if
  allocate(failed_backup_state, source=failed_coarse_state)
  allocate(failed_backup_temperature, source=failed_coarse_temperature)
  failed_backup_set = failed_set
  call advance_owned_reactive_eb_patch_set_chemistry_2d( &
    species, reactions, chemistry_interval, 1.0e-8_dp, 1.0e-14_dp, &
    distribution, failed_coarse_state, failed_coarse_temperature, &
    coarse_geometry, failed_set, ok, local_advances)
  call assert_all(.not. ok .and. local_advances == 0 .and. &
    all(failed_coarse_state == failed_backup_state) .and. &
    all(failed_coarse_temperature == failed_backup_temperature), &
    "MPI EB AMR owner chemistry failure rollback", rank)
  do child = 1, failed_set%patch_count()
    call assert_all(all(failed_set%children(child)%state == &
      failed_backup_set%children(child)%state) .and. &
      all(failed_set%children(child)%temperature == &
        failed_backup_set%children(child)%temperature), &
      "MPI EB AMR child chemistry rollback", rank)
  end do

  call scatter_owned_reactive_eb_patch_set_2d( &
    distribution, size(species), coarse_state, coarse_temperature, &
    coarse_geometry, patch_set, sparse_failed_set, ok)
  call assert_all(ok, "MPI EB AMR sparse chemistry failure scatter", rank)
  child = distribution%child_count()
  if (distribution%child_is_local(child)) &
    sparse_failed_set%children(child)%state(irho, :, :) = -1.0_dp
  sparse_failed_backup_set = sparse_failed_set
  call advance_sparse_owned_reactive_eb_patch_set_chemistry_2d( &
    species, reactions, chemistry_interval, 1.0e-8_dp, 1.0e-14_dp, &
    distribution, sparse_failed_set, coarse_geometry, patch_set, ok, &
    local_advances, local_restriction_transfers)
  call assert_all(.not. ok .and. local_advances == 0 .and. &
    local_restriction_transfers == 0 .and. &
    sparse_failed_set%local_value_count() == &
      sparse_failed_backup_set%local_value_count(), &
    "MPI EB AMR late sparse chemistry rollback", rank)
  ok = .true.
  do tile = 1, distribution%root_tile_count()
    if (.not. distribution%root_tile_is_local(tile)) cycle
    ok = ok .and. all(sparse_failed_set%root_tiles(tile)%state == &
        sparse_failed_backup_set%root_tiles(tile)%state) .and. &
      all(sparse_failed_set%root_tiles(tile)%temperature == &
        sparse_failed_backup_set%root_tiles(tile)%temperature)
  end do
  call assert_all(ok, "MPI EB AMR sparse chemistry root rollback", rank)
  ok = .true.
  do child = 1, distribution%child_count()
    if (.not. distribution%child_is_local(child)) cycle
    ok = ok .and. all(sparse_failed_set%children(child)%state == &
        sparse_failed_backup_set%children(child)%state) .and. &
      all(sparse_failed_set%children(child)%temperature == &
        sparse_failed_backup_set%children(child)%temperature)
  end do
  call assert_all(ok, "MPI EB AMR sparse chemistry child rollback", rank)

  call average_down_sparse_owned_reactive_eb_patch_set_2d( &
    species, distribution, sparse_failed_set, coarse_geometry, patch_set, &
    ok, local_restriction_transfers)
  call assert_all(.not. ok .and. local_restriction_transfers == 0 .and. &
    sparse_failed_set%local_value_count() == &
      sparse_failed_backup_set%local_value_count(), &
    "MPI EB AMR direct sparse average-down rejection", rank)
  ok = .true.
  do tile = 1, distribution%root_tile_count()
    if (.not. distribution%root_tile_is_local(tile)) cycle
    ok = ok .and. all(sparse_failed_set%root_tiles(tile)%state == &
        sparse_failed_backup_set%root_tiles(tile)%state) .and. &
      all(sparse_failed_set%root_tiles(tile)%temperature == &
        sparse_failed_backup_set%root_tiles(tile)%temperature)
  end do
  call assert_all(ok, "MPI EB AMR sparse average-down root rollback", rank)
  ok = .true.
  do child = 1, distribution%child_count()
    if (.not. distribution%child_is_local(child)) cycle
    ok = ok .and. all(sparse_failed_set%children(child)%state == &
        sparse_failed_backup_set%children(child)%state) .and. &
      all(sparse_failed_set%children(child)%temperature == &
        sparse_failed_backup_set%children(child)%temperature)
  end do
  call assert_all(ok, "MPI EB AMR sparse average-down child rollback", rank)

  hydro_start_set = patch_set
  do child = 1, hydro_start_set%patch_count()
    factor = 1.0_dp + 0.01_dp * real(3 - 2 * child, dp)
    hydro_start_set%children(child)%state = factor * &
      hydro_start_set%children(child)%state
  end do
  call assert_all(hydro_start_set%is_valid(coarse_geometry, nvar), &
    "MPI EB AMR hydro start hierarchy", rank)
  hydro_dt = chemistry_interval
  allocate(hydro_reference_state, mold=coarse_state)
  allocate(hydro_reference_temperature, mold=coarse_temperature)
  call advance_reactive_eb_patch_set_hydro_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, &
    hydro_start_set, "hllc", "pcm", "mc", 2, hydro_dt, &
    hydro_reference_state, hydro_reference_temperature, &
    hydro_reference_set, ok, &
    state_redist_target_volume_fraction=0.5_dp)
  call assert_all(ok, "serial EB AMR hydro reference", rank)

  allocate(hydro_mpi_state, source=coarse_state)
  allocate(hydro_mpi_temperature, source=coarse_temperature)
  hydro_mpi_set = hydro_start_set
  call advance_owned_reactive_eb_patch_set_hydro_2d( &
    species, distribution, hydro_mpi_state, hydro_mpi_temperature, &
    coarse_geometry, hydro_mpi_set, "hllc", "pcm", "mc", 2, hydro_dt, &
    ok, local_advances, 0.5_dp)
  call assert_all(ok, "MPI owner-only EB AMR hydro", rank)
  expected_local_advances = 0
  if (rank == distribution%root_level_owner()) &
    expected_local_advances = expected_local_advances + 1
  expected_global_advances = 1
  do child = 1, distribution%child_count()
    expected_global_advances = expected_global_advances + &
      hydro_start_set%children(child)%patch%refinement_ratio
    if (distribution%child_is_local(child)) &
      expected_local_advances = expected_local_advances + &
        hydro_start_set%children(child)%patch%refinement_ratio
  end do
  call MPI_Allreduce( &
    local_advances, global_advances, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_advances == expected_local_advances .and. &
    global_advances == expected_global_advances, &
    "MPI EB AMR one hydro advance per owner interval", rank)
  hydro_scale = max(1.0_dp, maxval(abs(hydro_reference_state)))
  call assert_all(maxval(abs(hydro_mpi_state - hydro_reference_state)) <= &
    8.0e-12_dp * hydro_scale .and. &
    maxval(abs(hydro_mpi_temperature - hydro_reference_temperature)) <= &
      8.0e-12_dp * max(1.0_dp, &
        maxval(abs(hydro_reference_temperature))), &
    "MPI EB AMR hydro serial root parity", rank)
  do child = 1, hydro_mpi_set%patch_count()
    call assert_all(maxval(abs(hydro_mpi_set%children(child)%state - &
      hydro_reference_set%children(child)%state)) <= &
        8.0e-12_dp * hydro_scale .and. &
      maxval(abs(hydro_mpi_set%children(child)%temperature - &
        hydro_reference_set%children(child)%temperature)) <= &
          8.0e-12_dp * max(1.0_dp, maxval(abs( &
            hydro_reference_set%children(child)%temperature))), &
      "MPI EB AMR hydro serial child parity", rank)
  end do
  hydro_change = maxval(abs(hydro_reference_state - coarse_state))
  do child = 1, hydro_reference_set%patch_count()
    hydro_change = max(hydro_change, maxval(abs( &
      hydro_reference_set%children(child)%state - &
      hydro_start_set%children(child)%state)))
  end do
  call assert_all(hydro_change > 1.0e-14_dp * hydro_scale, &
    "MPI EB AMR hydro changes nonuniform hierarchy", rank)

  call scatter_owned_reactive_eb_patch_set_2d( &
    distribution, size(species), coarse_state, coarse_temperature, &
    coarse_geometry, hydro_start_set, sparse_hydro_set, ok)
  call assert_all(ok, "MPI EB AMR sparse hydro scatter", rank)
  root_owner = distribution%root_level_owner()
  expected_local_root_transfers = 0
  expected_global_root_transfers = 0
  restriction_recipients = .false.
  do tile = 1, distribution%root_tile_count()
    owner = distribution%root_tiles(tile)%owner
    if (owner == root_owner) cycle
    expected_global_root_transfers = expected_global_root_transfers + 2
    if (rank == owner) expected_local_root_transfers = &
      expected_local_root_transfers + 1
    if (rank == root_owner) expected_local_root_transfers = &
      expected_local_root_transfers + 1
  end do
  do child = 1, distribution%child_count()
    owner = distribution%child_owner(child)
    if (owner == root_owner) cycle
    restriction_recipients(owner + 1) = .true.
    expected_global_root_transfers = expected_global_root_transfers + 2
    if (rank == root_owner) expected_local_root_transfers = &
      expected_local_root_transfers + 1
    if (rank == owner) expected_local_root_transfers = &
      expected_local_root_transfers + 1
  end do
  restriction_recipients(root_owner + 1) = .false.
  expected_global_root_transfers = expected_global_root_transfers + &
    count(restriction_recipients)
  if (rank == root_owner) expected_local_root_transfers = &
    expected_local_root_transfers + count(restriction_recipients)
  call advance_sparse_owned_reactive_eb_patch_set_hydro_2d( &
    species, distribution, sparse_hydro_set, coarse_geometry, &
    hydro_start_set, "hllc", "pcm", "mc", 2, hydro_dt, ok, &
    local_advances, 0.5_dp, local_root_transfers)
  call assert_all(ok .and. sparse_hydro_set%is_valid( &
    distribution, coarse_geometry, hydro_start_set) .and. &
    int(sparse_hydro_set%local_value_count()) == sparse_expected_local_values, &
    "MPI EB AMR direct sparse hydro", rank)
  call MPI_Allreduce( &
    local_advances, global_advances, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_advances == expected_local_advances .and. &
    global_advances == expected_global_advances, &
    "MPI EB AMR sparse hydro owner accounting", rank)
  call MPI_Allreduce( &
    local_root_transfers, global_root_transfers, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_root_transfers == expected_local_root_transfers .and. &
    global_root_transfers == expected_global_root_transfers, &
    "MPI EB AMR targeted sparse hydro root traffic", rank)
  call materialize_owned_reactive_eb_patch_set_2d( &
    distribution, sparse_hydro_set, coarse_state, coarse_temperature, &
    coarse_geometry, hydro_start_set, materialized_coarse_state, &
    materialized_coarse_temperature, materialized_patch_set, ok)
  call assert_all(ok .and. &
    maxval(abs(materialized_coarse_state - hydro_reference_state)) <= &
      8.0e-12_dp * hydro_scale .and. &
    maxval(abs(materialized_coarse_temperature - &
      hydro_reference_temperature)) <= &
        8.0e-12_dp * max(1.0_dp, &
          maxval(abs(hydro_reference_temperature))), &
    "MPI EB AMR sparse hydro serial root parity", rank)
  do child = 1, materialized_patch_set%patch_count()
    call assert_all( &
      maxval(abs(materialized_patch_set%children(child)%state - &
        hydro_reference_set%children(child)%state)) <= &
          8.0e-12_dp * hydro_scale .and. &
      maxval(abs(materialized_patch_set%children(child)%temperature - &
        hydro_reference_set%children(child)%temperature)) <= &
          8.0e-12_dp * max(1.0_dp, maxval(abs( &
            hydro_reference_set%children(child)%temperature))), &
      "MPI EB AMR sparse hydro serial child parity", rank)
  end do

  allocate(hydro_failed_state, source=coarse_state)
  allocate(hydro_failed_temperature, source=coarse_temperature)
  hydro_failed_set = hydro_start_set
  child = hydro_failed_set%patch_count()
  if (distribution%child_is_local(child)) &
    hydro_failed_set%children(child)%state(irho, :, :) = -1.0_dp
  allocate(hydro_failed_backup_state, source=hydro_failed_state)
  allocate(hydro_failed_backup_temperature, source=hydro_failed_temperature)
  hydro_failed_backup_set = hydro_failed_set
  call advance_owned_reactive_eb_patch_set_hydro_2d( &
    species, distribution, hydro_failed_state, hydro_failed_temperature, &
    coarse_geometry, hydro_failed_set, "hllc", "pcm", "mc", 2, hydro_dt, &
    ok, local_advances, 0.5_dp)
  call assert_all(.not. ok .and. local_advances == 0 .and. &
    all(hydro_failed_state == hydro_failed_backup_state) .and. &
    all(hydro_failed_temperature == hydro_failed_backup_temperature), &
    "MPI EB AMR late hydro failure root rollback", rank)
  do child = 1, hydro_failed_set%patch_count()
    call assert_all(all(hydro_failed_set%children(child)%state == &
      hydro_failed_backup_set%children(child)%state) .and. &
      all(hydro_failed_set%children(child)%temperature == &
        hydro_failed_backup_set%children(child)%temperature), &
      "MPI EB AMR late hydro failure child rollback", rank)
  end do

  call scatter_owned_reactive_eb_patch_set_2d( &
    distribution, size(species), coarse_state, coarse_temperature, &
    coarse_geometry, hydro_start_set, sparse_hydro_failed_set, ok)
  call assert_all(ok, "MPI EB AMR sparse hydro failure scatter", rank)
  child = distribution%child_count()
  if (distribution%child_is_local(child)) &
    sparse_hydro_failed_set%children(child)%state(irho, :, :) = -1.0_dp
  sparse_hydro_failed_backup_set = sparse_hydro_failed_set
  call advance_sparse_owned_reactive_eb_patch_set_hydro_2d( &
    species, distribution, sparse_hydro_failed_set, coarse_geometry, &
    hydro_start_set, "hllc", "pcm", "mc", 2, hydro_dt, ok, &
    local_advances, 0.5_dp, local_root_transfers)
  call assert_all(.not. ok .and. local_advances == 0 .and. &
    local_root_transfers == 0 .and. &
    sparse_hydro_failed_set%local_value_count() == &
      sparse_hydro_failed_backup_set%local_value_count(), &
    "MPI EB AMR late sparse hydro rollback", rank)
  ok = .true.
  do tile = 1, distribution%root_tile_count()
    if (.not. distribution%root_tile_is_local(tile)) cycle
    ok = ok .and. all(sparse_hydro_failed_set%root_tiles(tile)%state == &
        sparse_hydro_failed_backup_set%root_tiles(tile)%state) .and. &
      all(sparse_hydro_failed_set%root_tiles(tile)%temperature == &
        sparse_hydro_failed_backup_set%root_tiles(tile)%temperature)
  end do
  call assert_all(ok, "MPI EB AMR sparse hydro root rollback", rank)
  ok = .true.
  do child = 1, distribution%child_count()
    if (.not. distribution%child_is_local(child)) cycle
    ok = ok .and. all(sparse_hydro_failed_set%children(child)%state == &
        sparse_hydro_failed_backup_set%children(child)%state) .and. &
      all(sparse_hydro_failed_set%children(child)%temperature == &
        sparse_hydro_failed_backup_set%children(child)%temperature)
  end do
  call assert_all(ok, "MPI EB AMR sparse hydro child rollback", rank)

  transport_start_set = patch_set
  do child = 1, transport_start_set%patch_count()
    primitive(1) = 0.31_dp + 0.015_dp * real(3 - 2 * child, dp)
    primitive(2) = 0.035_dp * real(3 - 2 * child, dp)
    call reactive_primitive_to_conserved( &
      species, primitive, state_cell, child_temperature, child_sound_speed, ok)
    call assert_all(ok, "MPI EB AMR transport child state", rank)
    transport_start_set%children(child)%state = spread( &
      spread(state_cell, 2, &
        transport_start_set%children(child)%geometry%nx), 3, &
      transport_start_set%children(child)%geometry%ny)
    transport_start_set%children(child)%temperature = child_temperature
  end do
  call assert_all(transport_start_set%is_valid(coarse_geometry, nvar), &
    "MPI EB AMR transport start hierarchy", rank)
  transport_dt = 0.25_dp * hydro_dt
  allocate(transport_reference_state, mold=coarse_state)
  allocate(transport_reference_temperature, mold=coarse_temperature)
  call advance_reactive_eb_patch_set_transport_2d( &
    species, transport, coarse_state, coarse_temperature, coarse_geometry, &
    transport_start_set, transport_dt, .true., .true., .true., .true., &
    boundaries, 0.5_dp, 2, transport_reference_state, &
    transport_reference_temperature, transport_reference_set, &
    transport_reference_theta, ok)
  call assert_all(ok, "serial EB AMR transport reference", rank)

  call compute_reactive_eb_patch_set_cfl_timestep_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, &
    transport_start_set, 0.35_dp, timestep_reference_dt, ok)
  call assert_all(ok, "serial EB AMR hydro timestep reference", rank)
  call reactive_eb_transport_timestep_2d( &
    species, transport, coarse_state, coarse_temperature, coarse_geometry, &
    0.20_dp, .true., .true., .true., timestep_entity_dt, &
    timestep_maximum_diffusivity, ok)
  call assert_all(ok, "serial EB AMR root transport timestep", rank)
  timestep_reference_dt = min(timestep_reference_dt, timestep_entity_dt)
  do child = 1, transport_start_set%patch_count()
    call reactive_eb_transport_timestep_2d( &
      species, transport, transport_start_set%children(child)%state, &
      transport_start_set%children(child)%temperature, &
      transport_start_set%children(child)%geometry, 0.20_dp, .true., &
      .true., .true., timestep_entity_dt, timestep_maximum_diffusivity, ok)
    call assert_all(ok, "serial EB AMR child transport timestep", rank)
    timestep_reference_dt = min(timestep_reference_dt, &
      real(transport_start_set%children(child)%patch%refinement_ratio, dp) * &
        timestep_entity_dt)
  end do
  call scatter_owned_reactive_eb_patch_set_2d( &
    distribution, size(species), coarse_state, coarse_temperature, &
    coarse_geometry, transport_start_set, sparse_timestep_set, ok)
  call assert_all(ok, "MPI EB AMR sparse timestep scatter", rank)
  expected_local_root_transfers = 0
  expected_global_root_transfers = 0
  root_owner = distribution%root_level_owner()
  do tile = 1, distribution%root_tile_count()
    owner = distribution%root_tiles(tile)%owner
    if (owner == root_owner) cycle
    expected_global_root_transfers = expected_global_root_transfers + 1
    if (rank == owner) expected_local_root_transfers = &
      expected_local_root_transfers + 1
  end do
  call compute_sparse_owned_reactive_eb_patch_set_timestep_2d( &
    species, transport, distribution, sparse_timestep_set, coarse_geometry, &
    transport_start_set, 0.35_dp, 0.20_dp, .true., .true., .true., .true., &
    boundaries, timestep_dt, ok, local_root_transfers)
  call assert_all(ok .and. abs(timestep_dt - timestep_reference_dt) <= &
      64.0_dp * epsilon(1.0_dp) * timestep_reference_dt, &
    "MPI EB AMR sparse timestep serial parity", rank)
  call MPI_Allreduce( &
    local_root_transfers, global_root_transfers, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_root_transfers == expected_local_root_transfers .and. &
    global_root_transfers == expected_global_root_transfers, &
    "MPI EB AMR targeted sparse timestep root gather", rank)

  sparse_timestep_failed_set = sparse_timestep_set
  child = distribution%child_count()
  if (distribution%child_is_local(child)) &
    sparse_timestep_failed_set%children(child)%state(irho, :, :) = -1.0_dp
  sparse_timestep_failed_backup_set = sparse_timestep_failed_set
  call compute_sparse_owned_reactive_eb_patch_set_timestep_2d( &
    species, transport, distribution, sparse_timestep_failed_set, &
    coarse_geometry, transport_start_set, 0.35_dp, 0.20_dp, .true., .true., &
    .true., .true., boundaries, timestep_dt, ok, local_root_transfers)
  call assert_all(.not. ok .and. timestep_dt == 0.0_dp .and. &
    local_root_transfers == 0 .and. &
    sparse_timestep_failed_set%local_value_count() == &
      sparse_timestep_failed_backup_set%local_value_count(), &
    "MPI EB AMR sparse timestep rejection", rank)
  ok = .true.
  do tile = 1, distribution%root_tile_count()
    if (.not. distribution%root_tile_is_local(tile)) cycle
    ok = ok .and. &
      all(sparse_timestep_failed_set%root_tiles(tile)%state == &
        sparse_timestep_failed_backup_set%root_tiles(tile)%state) .and. &
      all(sparse_timestep_failed_set%root_tiles(tile)%temperature == &
        sparse_timestep_failed_backup_set%root_tiles(tile)%temperature)
  end do
  do child = 1, distribution%child_count()
    if (.not. distribution%child_is_local(child)) cycle
    ok = ok .and. all(sparse_timestep_failed_set%children(child)%state == &
        sparse_timestep_failed_backup_set%children(child)%state) .and. &
      all(sparse_timestep_failed_set%children(child)%temperature == &
        sparse_timestep_failed_backup_set%children(child)%temperature)
  end do
  call assert_all(ok, "MPI EB AMR sparse timestep state preservation", rank)

  allocate(transport_mpi_state, source=coarse_state)
  allocate(transport_mpi_temperature, source=coarse_temperature)
  transport_mpi_set = transport_start_set
  do tile = 1, distribution%root_tile_count()
    if (.not. distribution%root_tile_is_local(tile)) then
      transport_mpi_state(:, :, &
        distribution%root_tiles(tile)%j_lower: &
          distribution%root_tiles(tile)%j_upper) = 0.75_dp * &
        transport_mpi_state(:, :, &
          distribution%root_tiles(tile)%j_lower: &
            distribution%root_tiles(tile)%j_upper)
      transport_mpi_temperature(:, &
        distribution%root_tiles(tile)%j_lower: &
          distribution%root_tiles(tile)%j_upper) = 0.9_dp * &
        transport_mpi_temperature(:, &
          distribution%root_tiles(tile)%j_lower: &
            distribution%root_tiles(tile)%j_upper)
    end if
  end do
  do child = 1, distribution%child_count()
    if (.not. distribution%child_is_local(child)) then
      transport_mpi_set%children(child)%state = 0.75_dp * &
        transport_mpi_set%children(child)%state
      transport_mpi_set%children(child)%temperature = 0.9_dp * &
        transport_mpi_set%children(child)%temperature
    end if
  end do
  call advance_owned_reactive_eb_patch_set_transport_2d( &
    species, transport, distribution, transport_mpi_state, &
    transport_mpi_temperature, coarse_geometry, transport_mpi_set, &
    transport_dt, .true., .true., .true., .true., boundaries, 2, ok, &
    local_advances, transport_theta, 0.5_dp)
  call assert_all(ok, "MPI owner-only EB AMR transport", rank)
  expected_local_advances = 0
  if (rank == distribution%root_level_owner()) &
    expected_local_advances = expected_local_advances + 2
  expected_global_advances = 2
  do child = 1, distribution%child_count()
    expected_global_advances = expected_global_advances + 2 * &
      transport_start_set%children(child)%patch%refinement_ratio
    if (distribution%child_is_local(child)) &
      expected_local_advances = expected_local_advances + 2 * &
        transport_start_set%children(child)%patch%refinement_ratio
  end do
  call MPI_Allreduce( &
    local_advances, global_advances, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_advances == expected_local_advances .and. &
    global_advances == expected_global_advances, &
    "MPI EB AMR one transport Euler advance per owner interval", rank)
  transport_scale = max(1.0_dp, maxval(abs(transport_reference_state)))
  call assert_all(maxval(abs(transport_mpi_state - &
    transport_reference_state)) <= 2.0e-11_dp * transport_scale .and. &
    maxval(abs(transport_mpi_temperature - &
      transport_reference_temperature)) <= 2.0e-11_dp * &
        max(1.0_dp, maxval(abs(transport_reference_temperature))) .and. &
    abs(transport_theta - transport_reference_theta) <= &
      2.0e-13_dp * max(1.0_dp, abs(transport_reference_theta)), &
    "MPI EB AMR transport serial root parity", rank)
  do child = 1, transport_mpi_set%patch_count()
    call assert_all(maxval(abs(transport_mpi_set%children(child)%state - &
      transport_reference_set%children(child)%state)) <= &
        2.0e-11_dp * transport_scale .and. &
      maxval(abs(transport_mpi_set%children(child)%temperature - &
        transport_reference_set%children(child)%temperature)) <= &
          2.0e-11_dp * max(1.0_dp, maxval(abs( &
            transport_reference_set%children(child)%temperature))), &
      "MPI EB AMR transport serial child parity", rank)
  end do
  transport_change = maxval(abs( &
    transport_reference_state - coarse_state))
  do child = 1, transport_reference_set%patch_count()
    transport_change = max(transport_change, maxval(abs( &
      transport_reference_set%children(child)%state - &
      transport_start_set%children(child)%state)))
  end do
  call assert_all(transport_change > 1.0e-14_dp * transport_scale, &
    "MPI EB AMR transport changes nonuniform hierarchy", rank)

  call scatter_owned_reactive_eb_patch_set_2d( &
    distribution, size(species), coarse_state, coarse_temperature, &
    coarse_geometry, transport_start_set, sparse_transport_set, ok)
  call assert_all(ok, "MPI EB AMR sparse transport scatter", rank)
  root_owner = distribution%root_level_owner()
  expected_local_root_transfers = 0
  expected_global_root_transfers = 0
  restriction_recipients = .false.
  do tile = 1, distribution%root_tile_count()
    owner = distribution%root_tiles(tile)%owner
    if (owner == root_owner) cycle
    expected_global_root_transfers = expected_global_root_transfers + 7
    if (rank == owner) expected_local_root_transfers = &
      expected_local_root_transfers + 4
    if (rank == root_owner) expected_local_root_transfers = &
      expected_local_root_transfers + 3
  end do
  do child = 1, distribution%child_count()
    owner = distribution%child_owner(child)
    if (owner == root_owner) cycle
    restriction_recipients(owner + 1) = .true.
    expected_global_root_transfers = expected_global_root_transfers + 4
    if (rank == root_owner) expected_local_root_transfers = &
      expected_local_root_transfers + 2
    if (rank == owner) expected_local_root_transfers = &
      expected_local_root_transfers + 2
  end do
  restriction_recipients(root_owner + 1) = .false.
  expected_global_root_transfers = expected_global_root_transfers + &
    2 * count(restriction_recipients)
  if (rank == root_owner) expected_local_root_transfers = &
    expected_local_root_transfers + 2 * count(restriction_recipients)
  call advance_sparse_owned_reactive_eb_patch_set_transport_2d( &
    species, transport, distribution, sparse_transport_set, &
    coarse_geometry, transport_start_set, transport_dt, .true., .true., &
    .true., .true., boundaries, 2, ok, local_advances, transport_theta, &
    0.5_dp, local_root_transfers)
  call assert_all(ok .and. sparse_transport_set%is_valid( &
    distribution, coarse_geometry, transport_start_set) .and. &
    int(sparse_transport_set%local_value_count()) == &
      sparse_expected_local_values, &
    "MPI EB AMR direct sparse transport", rank)
  call MPI_Allreduce( &
    local_advances, global_advances, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_advances == expected_local_advances .and. &
    global_advances == expected_global_advances, &
    "MPI EB AMR sparse transport owner accounting", rank)
  call MPI_Allreduce( &
    local_root_transfers, global_root_transfers, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_root_transfers == expected_local_root_transfers .and. &
    global_root_transfers == expected_global_root_transfers, &
    "MPI EB AMR targeted sparse transport root traffic", rank)
  call materialize_owned_reactive_eb_patch_set_2d( &
    distribution, sparse_transport_set, coarse_state, coarse_temperature, &
    coarse_geometry, transport_start_set, materialized_coarse_state, &
    materialized_coarse_temperature, materialized_patch_set, ok)
  call assert_all(ok .and. &
    maxval(abs(materialized_coarse_state - transport_reference_state)) <= &
      2.0e-11_dp * transport_scale .and. &
    maxval(abs(materialized_coarse_temperature - &
      transport_reference_temperature)) <= &
        2.0e-11_dp * max(1.0_dp, &
          maxval(abs(transport_reference_temperature))) .and. &
    abs(transport_theta - transport_reference_theta) <= &
      2.0e-13_dp * max(1.0_dp, abs(transport_reference_theta)), &
    "MPI EB AMR sparse transport serial root parity", rank)
  do child = 1, materialized_patch_set%patch_count()
    call assert_all( &
      maxval(abs(materialized_patch_set%children(child)%state - &
        transport_reference_set%children(child)%state)) <= &
          2.0e-11_dp * transport_scale .and. &
      maxval(abs(materialized_patch_set%children(child)%temperature - &
        transport_reference_set%children(child)%temperature)) <= &
          2.0e-11_dp * max(1.0_dp, maxval(abs( &
            transport_reference_set%children(child)%temperature))), &
      "MPI EB AMR sparse transport serial child parity", rank)
  end do

  allocate(transport_failed_state, source=coarse_state)
  allocate(transport_failed_temperature, source=coarse_temperature)
  transport_failed_set = transport_start_set
  child = transport_failed_set%patch_count()
  if (distribution%child_is_local(child)) &
    transport_failed_set%children(child)%state(irho, :, :) = -1.0_dp
  allocate(transport_failed_backup_state, source=transport_failed_state)
  allocate(transport_failed_backup_temperature, &
    source=transport_failed_temperature)
  transport_failed_backup_set = transport_failed_set
  call advance_owned_reactive_eb_patch_set_transport_2d( &
    species, transport, distribution, transport_failed_state, &
    transport_failed_temperature, coarse_geometry, transport_failed_set, &
    transport_dt, .true., .true., .true., .true., boundaries, 2, ok, &
    local_advances, transport_theta, 0.5_dp)
  call assert_all(.not. ok .and. local_advances == 0 .and. &
    all(transport_failed_state == transport_failed_backup_state) .and. &
    all(transport_failed_temperature == &
      transport_failed_backup_temperature), &
    "MPI EB AMR late transport failure root rollback", rank)
  do child = 1, transport_failed_set%patch_count()
    call assert_all(all(transport_failed_set%children(child)%state == &
      transport_failed_backup_set%children(child)%state) .and. &
      all(transport_failed_set%children(child)%temperature == &
        transport_failed_backup_set%children(child)%temperature), &
      "MPI EB AMR late transport failure child rollback", rank)
  end do

  call scatter_owned_reactive_eb_patch_set_2d( &
    distribution, size(species), coarse_state, coarse_temperature, &
    coarse_geometry, transport_start_set, sparse_transport_failed_set, ok)
  call assert_all(ok, "MPI EB AMR sparse transport failure scatter", rank)
  child = distribution%child_count()
  if (distribution%child_is_local(child)) &
    sparse_transport_failed_set%children(child)%state(irho, :, :) = -1.0_dp
  sparse_transport_failed_backup_set = sparse_transport_failed_set
  call advance_sparse_owned_reactive_eb_patch_set_transport_2d( &
    species, transport, distribution, sparse_transport_failed_set, &
    coarse_geometry, transport_start_set, transport_dt, .true., .true., &
    .true., .true., boundaries, 2, ok, local_advances, transport_theta, &
    0.5_dp, local_root_transfers)
  call assert_all(.not. ok .and. local_advances == 0 .and. &
    transport_theta == 1.0_dp .and. &
    local_root_transfers == 0 .and. &
    sparse_transport_failed_set%local_value_count() == &
      sparse_transport_failed_backup_set%local_value_count(), &
    "MPI EB AMR late sparse transport rollback", rank)
  ok = .true.
  do tile = 1, distribution%root_tile_count()
    if (.not. distribution%root_tile_is_local(tile)) cycle
    ok = ok .and. all(sparse_transport_failed_set%root_tiles(tile)%state == &
        sparse_transport_failed_backup_set%root_tiles(tile)%state) .and. &
      all(sparse_transport_failed_set%root_tiles(tile)%temperature == &
        sparse_transport_failed_backup_set%root_tiles(tile)%temperature)
  end do
  call assert_all(ok, "MPI EB AMR sparse transport root rollback", rank)
  ok = .true.
  do child = 1, distribution%child_count()
    if (.not. distribution%child_is_local(child)) cycle
    ok = ok .and. all(sparse_transport_failed_set%children(child)%state == &
        sparse_transport_failed_backup_set%children(child)%state) .and. &
      all(sparse_transport_failed_set%children(child)%temperature == &
        sparse_transport_failed_backup_set%children(child)%temperature)
  end do
  call assert_all(ok, "MPI EB AMR sparse transport child rollback", rank)

  allocate(full_reference_state, mold=coarse_state)
  allocate(full_reference_temperature, mold=coarse_temperature)
  call advance_reactive_eb_patch_set_strang_2d( &
    species, reactions, coarse_state, coarse_temperature, coarse_geometry, &
    patch_set, "hllc", "pcm", "mc", 2, transport_dt, .true., &
    1.0e-8_dp, 1.0e-14_dp, full_reference_state, &
    full_reference_temperature, full_reference_set, ok, &
    target_volume_fraction=0.5_dp, failure_context=full_failure_context, &
    transport=transport, &
    transport_enabled=.true., viscosity_enabled=.true., &
    thermal_conduction_enabled=.true., species_diffusion_enabled=.true., &
    barodiffusion_enabled=.true., &
    minimum_transport_theta=full_reference_theta, boundaries=boundaries)
  call assert_all(ok, "serial EB AMR full-physics reference: " // &
    trim(full_failure_context), rank)

  allocate(full_mpi_state, source=coarse_state)
  allocate(full_mpi_temperature, source=coarse_temperature)
  full_mpi_set = patch_set
  call advance_owned_reactive_eb_patch_set_strang_2d( &
    species, reactions, transport, distribution, full_mpi_state, &
    full_mpi_temperature, coarse_geometry, full_mpi_set, "hllc", "pcm", &
    "mc", 2, transport_dt, 1.0e-8_dp, 1.0e-14_dp, .true., .true., &
    .true., .true., boundaries, ok, local_chemistry_advances, &
    local_hydro_advances, local_transport_advances, full_theta, 0.5_dp)
  call assert_all(ok, "MPI owner-only EB AMR full physics", rank)
  expected_local_chemistry = 2 * &
    distribution%rank_entity_counts(rank + 1)
  expected_global_chemistry = 2 * &
    (distribution%root_tile_count() + distribution%child_count())
  expected_local_hydro = 0
  expected_local_transport = 0
  if (rank == distribution%root_level_owner()) then
    expected_local_hydro = 1
    expected_local_transport = 4
  end if
  expected_global_hydro = 1
  expected_global_transport = 4
  do child = 1, distribution%child_count()
    expected_global_hydro = expected_global_hydro + &
      patch_set%children(child)%patch%refinement_ratio
    expected_global_transport = expected_global_transport + 4 * &
      patch_set%children(child)%patch%refinement_ratio
    if (distribution%child_is_local(child)) then
      expected_local_hydro = expected_local_hydro + &
        patch_set%children(child)%patch%refinement_ratio
      expected_local_transport = expected_local_transport + 4 * &
      patch_set%children(child)%patch%refinement_ratio
    end if
  end do
  call MPI_Allreduce( &
    local_chemistry_advances, global_chemistry_advances, 1, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS, &
    "MPI EB AMR full chemistry count reduction", rank)
  call MPI_Allreduce( &
    local_hydro_advances, global_hydro_advances, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS, &
    "MPI EB AMR full hydro count reduction", rank)
  call MPI_Allreduce( &
    local_transport_advances, global_transport_advances, 1, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_chemistry_advances == expected_local_chemistry .and. &
    global_chemistry_advances == expected_global_chemistry .and. &
    local_hydro_advances == expected_local_hydro .and. &
    global_hydro_advances == expected_global_hydro .and. &
    local_transport_advances == expected_local_transport .and. &
    global_transport_advances == expected_global_transport, &
    "MPI EB AMR full-physics owner accounting", rank)
  full_scale = max(1.0_dp, maxval(abs(full_reference_state)))
  call assert_all(maxval(abs(full_mpi_state - full_reference_state)) <= &
    5.0e-11_dp * full_scale .and. &
    maxval(abs(full_mpi_temperature - full_reference_temperature)) <= &
      5.0e-11_dp * max(1.0_dp, &
        maxval(abs(full_reference_temperature))) .and. &
    abs(full_theta - full_reference_theta) <= &
      3.0e-13_dp * max(1.0_dp, abs(full_reference_theta)), &
    "MPI EB AMR full-physics serial root parity", rank)
  do child = 1, full_mpi_set%patch_count()
    call assert_all(maxval(abs(full_mpi_set%children(child)%state - &
      full_reference_set%children(child)%state)) <= &
        5.0e-11_dp * full_scale .and. &
      maxval(abs(full_mpi_set%children(child)%temperature - &
        full_reference_set%children(child)%temperature)) <= &
          5.0e-11_dp * max(1.0_dp, maxval(abs( &
            full_reference_set%children(child)%temperature))), &
      "MPI EB AMR full-physics serial child parity", rank)
  end do
  full_change = maxval(abs(full_reference_state - coarse_state))
  do child = 1, full_reference_set%patch_count()
    full_change = max(full_change, maxval(abs( &
      full_reference_set%children(child)%state - &
      patch_set%children(child)%state)))
  end do
  call assert_all(full_change > 1.0e-14_dp * full_scale, &
    "MPI EB AMR full physics changes hierarchy", rank)

  call scatter_owned_reactive_eb_patch_set_2d( &
    distribution, size(species), coarse_state, coarse_temperature, &
    coarse_geometry, patch_set, sparse_full_set, ok)
  call assert_all(ok, "MPI EB AMR sparse full-physics scatter", rank)
  call advance_sparse_owned_reactive_eb_patch_set_strang_2d( &
    species, reactions, transport, distribution, sparse_full_set, &
    coarse_geometry, patch_set, "hllc", "pcm", "mc", 2, transport_dt, &
    1.0e-8_dp, 1.0e-14_dp, .true., .true., .true., .true., boundaries, &
    ok, local_chemistry_advances, local_hydro_advances, &
    local_transport_advances, full_theta, 0.5_dp)
  call assert_all(ok .and. sparse_full_set%is_valid( &
    distribution, coarse_geometry, patch_set) .and. &
    int(sparse_full_set%local_value_count()) == sparse_expected_local_values, &
    "MPI EB AMR sparse full-physics transaction", rank)
  call MPI_Allreduce( &
    local_chemistry_advances, global_chemistry_advances, 1, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS, &
    "MPI EB AMR sparse full chemistry count reduction", rank)
  call MPI_Allreduce( &
    local_hydro_advances, global_hydro_advances, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS, &
    "MPI EB AMR sparse full hydro count reduction", rank)
  call MPI_Allreduce( &
    local_transport_advances, global_transport_advances, 1, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_chemistry_advances == expected_local_chemistry .and. &
    global_chemistry_advances == expected_global_chemistry .and. &
    local_hydro_advances == expected_local_hydro .and. &
    global_hydro_advances == expected_global_hydro .and. &
    local_transport_advances == expected_local_transport .and. &
    global_transport_advances == expected_global_transport, &
    "MPI EB AMR sparse full-physics owner accounting", rank)
  call materialize_owned_reactive_eb_patch_set_2d( &
    distribution, sparse_full_set, coarse_state, coarse_temperature, &
    coarse_geometry, patch_set, materialized_coarse_state, &
    materialized_coarse_temperature, materialized_patch_set, ok)
  call assert_all(ok .and. &
    maxval(abs(materialized_coarse_state - full_reference_state)) <= &
      5.0e-11_dp * full_scale .and. &
    maxval(abs(materialized_coarse_temperature - &
      full_reference_temperature)) <= &
        5.0e-11_dp * max(1.0_dp, &
          maxval(abs(full_reference_temperature))) .and. &
    abs(full_theta - full_reference_theta) <= &
      3.0e-13_dp * max(1.0_dp, abs(full_reference_theta)), &
    "MPI EB AMR sparse full-physics serial root parity", rank)
  do child = 1, materialized_patch_set%patch_count()
    call assert_all( &
      maxval(abs(materialized_patch_set%children(child)%state - &
        full_reference_set%children(child)%state)) <= &
          5.0e-11_dp * full_scale .and. &
      maxval(abs(materialized_patch_set%children(child)%temperature - &
        full_reference_set%children(child)%temperature)) <= &
          5.0e-11_dp * max(1.0_dp, maxval(abs( &
            full_reference_set%children(child)%temperature))), &
      "MPI EB AMR sparse full-physics serial child parity", rank)
  end do

  allocate(full_failed_state, source=coarse_state)
  allocate(full_failed_temperature, source=coarse_temperature)
  full_failed_set = patch_set
  allocate(full_failed_backup_state, source=full_failed_state)
  allocate(full_failed_backup_temperature, source=full_failed_temperature)
  full_failed_backup_set = full_failed_set
  call advance_owned_reactive_eb_patch_set_strang_2d( &
    species, reactions, transport, distribution, full_failed_state, &
    full_failed_temperature, coarse_geometry, full_failed_set, &
    "missing_solver", "pcm", "mc", 2, transport_dt, 1.0e-8_dp, &
    1.0e-14_dp, .true., .true., .true., .true., boundaries, ok, &
    local_chemistry_advances, local_hydro_advances, &
    local_transport_advances, full_theta, 0.5_dp)
  call assert_all(.not. ok .and. local_chemistry_advances == 0 .and. &
    local_hydro_advances == 0 .and. local_transport_advances == 0 .and. &
    all(full_failed_state == full_failed_backup_state) .and. &
    all(full_failed_temperature == full_failed_backup_temperature), &
    "MPI EB AMR late full-physics root rollback", rank)
  do child = 1, full_failed_set%patch_count()
    call assert_all(all(full_failed_set%children(child)%state == &
      full_failed_backup_set%children(child)%state) .and. &
      all(full_failed_set%children(child)%temperature == &
        full_failed_backup_set%children(child)%temperature), &
      "MPI EB AMR late full-physics child rollback", rank)
  end do

  call scatter_owned_reactive_eb_patch_set_2d( &
    distribution, size(species), coarse_state, coarse_temperature, &
    coarse_geometry, patch_set, sparse_full_failed_set, ok)
  call assert_all(ok, "MPI EB AMR sparse full failure scatter", rank)
  sparse_full_failed_backup_set = sparse_full_failed_set
  call advance_sparse_owned_reactive_eb_patch_set_strang_2d( &
    species, reactions, transport, distribution, sparse_full_failed_set, &
    coarse_geometry, patch_set, "missing_solver", "pcm", "mc", 2, &
    transport_dt, 1.0e-8_dp, 1.0e-14_dp, .true., .true., .true., .true., &
    boundaries, ok, local_chemistry_advances, local_hydro_advances, &
    local_transport_advances, full_theta, 0.5_dp)
  call assert_all(.not. ok .and. local_chemistry_advances == 0 .and. &
    local_hydro_advances == 0 .and. local_transport_advances == 0 .and. &
    full_theta == 1.0_dp .and. &
    sparse_full_failed_set%local_value_count() == &
      sparse_full_failed_backup_set%local_value_count(), &
    "MPI EB AMR late sparse full-physics rollback", rank)
  ok = .true.
  do tile = 1, distribution%root_tile_count()
    if (.not. distribution%root_tile_is_local(tile)) cycle
    ok = ok .and. all(sparse_full_failed_set%root_tiles(tile)%state == &
        sparse_full_failed_backup_set%root_tiles(tile)%state) .and. &
      all(sparse_full_failed_set%root_tiles(tile)%temperature == &
        sparse_full_failed_backup_set%root_tiles(tile)%temperature)
  end do
  call assert_all(ok, "MPI EB AMR sparse full root rollback", rank)
  ok = .true.
  do child = 1, distribution%child_count()
    if (.not. distribution%child_is_local(child)) cycle
    ok = ok .and. all(sparse_full_failed_set%children(child)%state == &
        sparse_full_failed_backup_set%children(child)%state) .and. &
      all(sparse_full_failed_set%children(child)%temperature == &
        sparse_full_failed_backup_set%children(child)%temperature)
  end do
  call assert_all(ok, "MPI EB AMR sparse full child rollback", rank)

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
      " MPI EB AMR owner chemistry/hydro/transport 2D: PASS"
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
