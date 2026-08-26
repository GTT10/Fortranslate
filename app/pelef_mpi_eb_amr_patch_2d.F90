program pelef_mpi_eb_amr_patch_2d
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
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
    eb_geometry_2d, eb_covered_cell, eb_cut_cell, build_eb_geometry_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, build_amr_eb_patch_2d
  use amr_eb_regrid_2d_mod, only: &
    amr_eb_tagging_criteria_2d, amr_eb_regrid_plan_collection_2d, &
    reactive_eb_patch_set_2d, reactive_eb_patch_topology_2d, &
    build_amr_eb_regrid_plan_collection_2d, &
    plan_reactive_eb_temperature_regrid_collection_2d, &
    initialize_reactive_eb_patch_set_2d, &
    initialize_reactive_eb_patch_topology_2d, &
    average_down_reactive_eb_patch_set_2d, &
    advance_reactive_eb_patch_set_hydro_2d, &
    regrid_reactive_eb_patch_set_2d
  use amr_eb_multipatch_transport_2d_mod, only: &
    advance_reactive_eb_patch_set_transport_2d
  use eb_reactive_transport_2d_mod, only: &
    reactive_eb_transport_timestep_2d
  use simulation_config_reactive_eb_amr_2d_mod, only: &
    reactive_eb_amr_2d_config
  use reactive_eb_amr_2d_driver_mod, only: &
    advance_reactive_eb_patch_set_strang_2d, &
    compute_reactive_eb_patch_set_cfl_timestep_2d, &
    read_reactive_eb_amr_patch_set_2d_checkpoint
  use mpi_amr_eb_patch_2d_mod, only: &
    mpi_amr_eb_root_tile_hydro_halo_cells, &
    mpi_amr_eb_root_tile_transport_halo_cells, &
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
    gather_sparse_owned_reactive_eb_patch_set_to_root_2d, &
    scatter_root_reactive_eb_topology_to_sparse_2d, &
    regrid_sparse_owned_reactive_eb_patch_set_2d, &
    average_down_sparse_owned_reactive_eb_patch_set_2d, &
    advance_sparse_owned_reactive_eb_patch_set_chemistry_2d, &
    compute_sparse_owned_reactive_eb_patch_set_timestep_2d, &
    advance_sparse_owned_reactive_eb_patch_set_hydro_2d, &
    advance_sparse_owned_reactive_eb_patch_set_transport_2d, &
    advance_sparse_owned_reactive_eb_patch_set_strang_2d, &
    advance_sparse_owned_reactive_eb_patch_set_to_time_2d
  use mpi_amr_eb_io_2d_mod, only: &
    write_sparse_owned_reactive_eb_patch_set_2d_checkpoint, &
    read_sparse_owned_reactive_eb_topology_2d_checkpoint, &
    write_sparse_owned_reactive_eb_patch_set_2d_csv
  implicit none

  integer, parameter :: coarse_nx = 14, coarse_ny = 14, ratio = 2
  type(eb_geometry_2d) :: coarse_geometry
  type(eb_geometry_2d), allocatable :: fine_geometries(:)
  type(amr_eb_patch_2d) :: geometry_patch
  type(amr_eb_tagging_criteria_2d) :: criteria
  type(amr_eb_tagging_criteria_2d) :: scheduled_criteria
  type(amr_eb_regrid_plan_collection_2d) :: collection
  type(amr_eb_regrid_plan_collection_2d) :: regrid_collection
  type(amr_eb_regrid_plan_collection_2d) :: scheduled_collection
  type(reactive_eb_patch_set_2d) :: patch_set, local_patch_set
  type(reactive_eb_patch_set_2d) :: synchronized_patch_set
  type(reactive_eb_patch_set_2d) :: materialized_patch_set
  type(reactive_eb_patch_set_2d) :: root_materialized_patch_set
  type(reactive_eb_patch_set_2d) :: checkpoint_patch_set
  type(reactive_eb_patch_set_2d) :: rejected_materialized_set
  type(reactive_eb_patch_topology_2d) :: restart_topology
  type(reactive_eb_patch_topology_2d) :: invalid_restart_topology
  type(mpi_amr_eb_patch_distribution_2d) :: distribution
  type(mpi_amr_eb_patch_distribution_2d) :: topology_distribution
  type(mpi_amr_eb_patch_distribution_2d) :: invalid_distribution
  type(mpi_amr_eb_patch_distribution_2d) :: rejected_distribution
  type(mpi_amr_eb_patch_distribution_2d) :: scheduled_distribution
  type(mpi_amr_eb_patch_distribution_2d) :: scheduled_failed_distribution
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
  type(mpi_amr_eb_sparse_patch_set_2d) :: sparse_time_loop_set
  type(mpi_amr_eb_sparse_patch_set_2d) :: sparse_limited_time_loop_set
  type(mpi_amr_eb_sparse_patch_set_2d) :: sparse_regrid_set
  type(mpi_amr_eb_sparse_patch_set_2d) :: sparse_scheduled_set
  type(mpi_amr_eb_sparse_patch_set_2d) :: sparse_scheduled_failed_set
  type(mpi_amr_eb_sparse_patch_set_2d) :: sparse_restart_set
  type(reactive_eb_amr_2d_config) :: io_config
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
  type(reactive_eb_patch_set_2d) :: time_loop_reference_set
  type(reactive_eb_patch_set_2d) :: regrid_start_set
  type(reactive_eb_patch_set_2d) :: regrid_reference_set
  type(reactive_eb_patch_set_2d) :: sparse_regrid_template
  type(reactive_eb_patch_set_2d) :: scheduled_start_set
  type(reactive_eb_patch_set_2d) :: scheduled_reference_set
  type(reactive_eb_patch_set_2d) :: scheduled_template
  type(reactive_eb_patch_set_2d) :: scheduled_failed_template
  real(dp) :: coarse_level_set(0:coarse_nx, 0:coarse_ny)
  real(dp), allocatable :: primitive(:), mass_fractions(:), state_cell(:)
  real(dp), allocatable :: coarse_state(:, :, :), coarse_temperature(:, :)
  real(dp), allocatable :: local_coarse_state(:, :, :)
  real(dp), allocatable :: local_coarse_temperature(:, :)
  real(dp), allocatable :: synchronized_coarse_state(:, :, :)
  real(dp), allocatable :: synchronized_coarse_temperature(:, :)
  real(dp), allocatable :: materialized_coarse_state(:, :, :)
  real(dp), allocatable :: materialized_coarse_temperature(:, :)
  real(dp), allocatable :: root_materialized_state(:, :, :)
  real(dp), allocatable :: root_materialized_temperature(:, :)
  real(dp), allocatable :: checkpoint_state(:, :, :)
  real(dp), allocatable :: checkpoint_temperature(:, :)
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
  real(dp), allocatable :: time_loop_reference_state(:, :, :)
  real(dp), allocatable :: time_loop_reference_temperature(:, :)
  real(dp), allocatable :: regrid_start_state(:, :, :)
  real(dp), allocatable :: regrid_start_temperature(:, :)
  real(dp), allocatable :: regrid_reference_state(:, :, :)
  real(dp), allocatable :: regrid_reference_temperature(:, :)
  real(dp), allocatable :: scheduled_start_state(:, :, :)
  real(dp), allocatable :: scheduled_start_temperature(:, :)
  real(dp), allocatable :: scheduled_reference_state(:, :, :)
  real(dp), allocatable :: scheduled_reference_temperature(:, :)
  type(eb_geometry_2d), allocatable :: regrid_fine_geometries(:)
  type(eb_geometry_2d), allocatable :: rejected_regrid_fine_geometries(:)
  type(eb_geometry_2d) :: checkpoint_geometry
  logical, allocatable :: active_mask(:, :)
  logical, allocatable :: regrid_recipients(:)
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
  real(dp) :: time_loop_initial_dt, time_loop_final_time, time_loop_time
  real(dp) :: time_loop_minimum_dt, time_loop_reference_minimum_dt
  real(dp) :: time_loop_reference_theta
  real(dp) :: scheduled_initial_dt, scheduled_final_time
  real(dp) :: scheduled_time, scheduled_minimum_dt
  real(dp) :: scheduled_reference_minimum_dt, scheduled_reference_theta
  real(dp) :: scheduled_theta
  real(dp) :: checkpoint_time, checkpoint_minimum_dt
  real(dp) :: checkpoint_base_density
  logical :: tags(coarse_nx, coarse_ny), regrid_tags(coarse_nx, coarse_ny)
  logical :: geometry_perturbed, ok, topology_changed
  logical :: io_files_ok
  integer :: child, component, global_advances, global_i, global_j, i, ierr
  integer :: j, new_child, old_child, old_i, old_j, owner, source
  integer :: local_advances, nvar, rank, nranks, tile
  integer :: expected_global_advances, expected_local_advances
  integer :: local_root_hydro_cells, global_root_hydro_cells
  integer :: expected_local_root_hydro_cells
  integer :: expected_global_root_hydro_cells
  integer :: local_root_transport_cells, global_root_transport_cells
  integer :: expected_local_root_transport_cells
  integer :: expected_global_root_transport_cells
  integer :: local_chemistry_advances, local_hydro_advances
  integer :: local_transport_advances
  integer :: global_chemistry_advances, global_hydro_advances
  integer :: global_transport_advances
  integer :: expected_local_chemistry, expected_local_hydro
  integer :: expected_local_transport, expected_global_chemistry
  integer :: expected_global_hydro, expected_global_transport
  integer :: time_loop_steps, time_loop_reference_steps
  integer :: time_loop_advanced_steps
  integer :: scheduled_steps, scheduled_reference_steps
  integer :: scheduled_advanced_steps, scheduled_regrid_evaluations
  integer :: scheduled_reference_evaluations, scheduled_regrids
  integer :: scheduled_reference_regrids, scheduled_timestep_transfers
  integer :: scheduled_regrid_transfers, global_scheduled_regrid_transfers
  integer :: failing_geometry_calls
  integer :: inconsistent_exponent
  integer :: sparse_local_values, sparse_global_values
  integer :: sparse_expected_local_values, sparse_expected_global_values
  integer :: local_restriction_transfers, global_restriction_transfers
  integer :: expected_local_restriction_transfers
  integer :: expected_global_restriction_transfers
  integer :: local_regrid_restriction_transfers
  integer :: global_regrid_restriction_transfers
  integer :: expected_local_regrid_restriction_transfers
  integer :: expected_global_regrid_restriction_transfers
  integer :: local_regrid_prolongation_transfers
  integer :: global_regrid_prolongation_transfers
  integer :: expected_local_regrid_prolongation_transfers
  integer :: expected_global_regrid_prolongation_transfers
  integer :: local_regrid_overlap_transfers
  integer :: global_regrid_overlap_transfers
  integer :: expected_local_regrid_overlap_transfers
  integer :: expected_global_regrid_overlap_transfers
  integer :: local_root_transfers, global_root_transfers
  integer :: expected_local_root_transfers
  integer :: expected_global_root_transfers, root_owner
  integer :: local_root_materialization_transfers
  integer :: global_root_materialization_transfers
  integer :: expected_local_root_materialization_transfers
  integer :: expected_global_root_materialization_transfers
  integer :: invalid_materialization_root
  integer :: checkpoint_steps, checkpoint_regrids
  integer :: local_restart_transfers, global_restart_transfers
  integer :: expected_local_restart_transfers
  character(len=160) :: full_failure_context
  character(len=*), parameter :: sparse_checkpoint_path = &
    "pelef_mpi_eb_amr_sparse_io.chk"
  character(len=*), parameter :: sparse_root_output_path = &
    "pelef_mpi_eb_amr_sparse_root.csv"
  character(len=*), parameter :: sparse_fine_output_path = &
    "pelef_mpi_eb_amr_sparse_fine.csv"
  character(len=*), parameter :: sparse_fine_output_path_1 = &
    "pelef_mpi_eb_amr_sparse_fine_patch0001.csv"
  character(len=*), parameter :: sparse_fine_output_path_2 = &
    "pelef_mpi_eb_amr_sparse_fine_patch0002.csv"
  character(len=*), parameter :: missing_checkpoint_path = &
    "pelef_missing_output_directory/sparse_io.chk"
  character(len=*), parameter :: missing_root_output_path = &
    "pelef_missing_output_directory/sparse_root.csv"

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
  call initialize_reactive_eb_patch_topology_2d( &
    coarse_geometry, fine_geometries, collection, ratio, restart_topology, ok)
  call assert_all(ok .and. restart_topology%is_valid(coarse_geometry) .and. &
    restart_topology%patch_count() == collection%patch_count(), &
    "MPI EB AMR geometry-only restart topology", rank)
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

  root_owner = distribution%root_level_owner()
  expected_local_root_materialization_transfers = 0
  expected_global_root_materialization_transfers = 0
  do tile = 1, distribution%root_tile_count()
    owner = distribution%root_tiles(tile)%owner
    if (owner == root_owner) cycle
    expected_global_root_materialization_transfers = &
      expected_global_root_materialization_transfers + 1
    if (rank == owner) expected_local_root_materialization_transfers = &
      expected_local_root_materialization_transfers + 1
  end do
  do child = 1, distribution%child_count()
    owner = distribution%child_owner(child)
    if (owner == root_owner) cycle
    expected_global_root_materialization_transfers = &
      expected_global_root_materialization_transfers + 1
    if (rank == owner) expected_local_root_materialization_transfers = &
      expected_local_root_materialization_transfers + 1
  end do
  call gather_sparse_owned_reactive_eb_patch_set_to_root_2d( &
    distribution, sparse_patch_set, coarse_geometry, patch_set, root_owner, &
    root_materialized_state, root_materialized_temperature, &
    root_materialized_patch_set, ok, local_root_materialization_transfers)
  call assert_all(ok .and. &
    merge(allocated(root_materialized_state), &
      .not. allocated(root_materialized_state), rank == root_owner) .and. &
    merge(allocated(root_materialized_temperature), &
      .not. allocated(root_materialized_temperature), rank == root_owner) &
      .and. merge(root_materialized_patch_set%patch_count() == &
        patch_set%patch_count(), &
        root_materialized_patch_set%patch_count() == 0, &
        rank == root_owner), &
    "MPI EB AMR root-only materialization allocation", rank)
  if (rank == root_owner) then
    call assert_all( &
      all(root_materialized_state == synchronized_coarse_state) .and. &
      all(root_materialized_temperature == &
        synchronized_coarse_temperature), &
      "MPI EB AMR root-only root materialization", rank)
    do child = 1, distribution%child_count()
      call assert_all( &
        all(root_materialized_patch_set%children(child)%state == &
          synchronized_patch_set%children(child)%state) .and. &
        all(root_materialized_patch_set%children(child)%temperature == &
          synchronized_patch_set%children(child)%temperature), &
        "MPI EB AMR root-only child materialization", rank)
    end do
  else
    call assert_all(.true., "MPI EB AMR root-only root materialization", rank)
    do child = 1, distribution%child_count()
      call assert_all( &
        .true., "MPI EB AMR root-only child materialization", rank)
    end do
  end if
  call MPI_Allreduce( &
    local_root_materialization_transfers, &
    global_root_materialization_transfers, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_root_materialization_transfers == &
      expected_local_root_materialization_transfers .and. &
    global_root_materialization_transfers == &
      expected_global_root_materialization_transfers, &
    "MPI EB AMR root-only materialization traffic", rank)

  expected_local_restart_transfers = 0
  if (rank == root_owner) expected_local_restart_transfers = &
    expected_global_root_materialization_transfers
  call scatter_root_reactive_eb_topology_to_sparse_2d( &
    distribution, size(species), root_materialized_state, &
    root_materialized_temperature, coarse_geometry, &
    root_materialized_patch_set, restart_topology, root_owner, &
    sparse_restart_set, &
    ok, local_restart_transfers)
  call MPI_Allreduce( &
    local_restart_transfers, global_restart_transfers, 1, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ok .and. ierr == MPI_SUCCESS .and. &
    local_restart_transfers == expected_local_restart_transfers .and. &
    global_restart_transfers == &
      expected_global_root_materialization_transfers .and. &
    sparse_restart_set%local_value_count() == &
      sparse_patch_set%local_value_count(), &
    "MPI EB AMR direct root restart scatter accounting", rank)
  ok = .true.
  do tile = 1, distribution%root_tile_count()
    if (.not. distribution%root_tile_is_local(tile)) cycle
    ok = ok .and. all(sparse_restart_set%root_tiles(tile)%state == &
        sparse_patch_set%root_tiles(tile)%state) .and. &
      all(sparse_restart_set%root_tiles(tile)%temperature == &
        sparse_patch_set%root_tiles(tile)%temperature)
  end do
  do child = 1, distribution%child_count()
    if (.not. distribution%child_is_local(child)) cycle
    ok = ok .and. all(sparse_restart_set%children(child)%state == &
        sparse_patch_set%children(child)%state) .and. &
      all(sparse_restart_set%children(child)%temperature == &
        sparse_patch_set%children(child)%temperature)
  end do
  call assert_all(ok, "MPI EB AMR direct root restart scatter parity", rank)

  invalid_restart_topology = restart_topology
  invalid_restart_topology%children(1)%geometry%volume_fraction(1, 1) = &
    2.0_dp
  call scatter_root_reactive_eb_topology_to_sparse_2d( &
    distribution, size(species), root_materialized_state, &
    root_materialized_temperature, coarse_geometry, &
    root_materialized_patch_set, invalid_restart_topology, root_owner, &
    sparse_restart_set, ok, local_restart_transfers)
  call assert_all(.not. ok .and. local_restart_transfers == 0 .and. &
    sparse_restart_set%local_value_count() == 0, &
    "MPI EB AMR invalid geometry-only restart topology rollback", rank)

  io_config%eb%flow%nx = coarse_nx
  io_config%eb%flow%ny = coarse_ny
  io_config%eb%flow%x_lower = coarse_geometry%x_lower
  io_config%eb%flow%x_upper = coarse_geometry%x_upper
  io_config%eb%flow%y_lower = coarse_geometry%y_lower
  io_config%eb%flow%y_upper = coarse_geometry%y_upper
  io_config%eb%flow%final_time = 2.0e-6_dp
  io_config%eb%flow%reconstruction = "pcm"
  io_config%eb%flow%use_transverse_correction = .false.
  io_config%eb%flow%boundary_x_lower = "outflow"
  io_config%eb%flow%boundary_x_upper = "outflow"
  io_config%eb%flow%boundary_y_lower = "outflow"
  io_config%eb%flow%boundary_y_upper = "outflow"
  io_config%eb%flow%output_file = sparse_root_output_path
  io_config%eb%plane_normal_x = 1.0_dp
  io_config%eb%plane_normal_y = 1.0_dp
  io_config%eb%plane_offset = 0.78_dp
  io_config%coarse_i_lower = collection%plans(1)%coarse_i_lower
  io_config%coarse_i_upper = collection%plans(1)%coarse_i_upper
  io_config%coarse_j_lower = collection%plans(1)%coarse_j_lower
  io_config%coarse_j_upper = collection%plans(1)%coarse_j_upper
  io_config%refinement_ratio = ratio
  io_config%multipatch_enabled = .true.
  io_config%dynamic_regridding = .true.
  io_config%fine_output_file = sparse_fine_output_path
  call write_sparse_owned_reactive_eb_patch_set_2d_checkpoint( &
    sparse_checkpoint_path, species, io_config, distribution, &
    sparse_patch_set, coarse_geometry, patch_set, root_owner, 1.0e-6_dp, &
    1, 0, 1.0e-6_dp, 0.31_dp, ok, &
    local_root_materialization_transfers)
  call MPI_Allreduce( &
    local_root_materialization_transfers, &
    global_root_materialization_transfers, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ok .and. ierr == MPI_SUCCESS .and. &
    local_root_materialization_transfers == &
      expected_local_root_materialization_transfers .and. &
    global_root_materialization_transfers == &
      expected_global_root_materialization_transfers, &
    "MPI EB AMR sparse checkpoint transfer accounting", rank)
  io_files_ok = .true.
  if (rank == root_owner) then
    call read_reactive_eb_amr_patch_set_2d_checkpoint( &
      sparse_checkpoint_path, species, io_config, checkpoint_state, &
      checkpoint_temperature, checkpoint_geometry, checkpoint_patch_set, &
      checkpoint_time, checkpoint_steps, checkpoint_regrids, &
      checkpoint_minimum_dt, checkpoint_base_density, io_files_ok)
    if (io_files_ok) then
      io_files_ok = all(checkpoint_state == synchronized_coarse_state) .and. &
        maxval(abs(checkpoint_temperature - &
          synchronized_coarse_temperature)) <= &
          3.0e-12_dp * max(1.0_dp, &
            maxval(abs(synchronized_coarse_temperature))) .and. &
        checkpoint_geometry%is_valid() .and. &
        checkpoint_patch_set%patch_count() == patch_set%patch_count() .and. &
        checkpoint_time == 1.0e-6_dp .and. checkpoint_steps == 1 .and. &
        checkpoint_regrids == 0 .and. &
        checkpoint_minimum_dt == 1.0e-6_dp .and. &
        checkpoint_base_density == 0.31_dp
    end if
    if (io_files_ok) then
      do child = 1, patch_set%patch_count()
        io_files_ok = io_files_ok .and. &
          all(checkpoint_patch_set%children(child)%state == &
            synchronized_patch_set%children(child)%state) .and. &
          maxval(abs(checkpoint_patch_set%children(child)%temperature - &
            synchronized_patch_set%children(child)%temperature)) <= &
            3.0e-12_dp * max(1.0_dp, maxval(abs( &
              synchronized_patch_set%children(child)%temperature)))
      end do
    end if
  end if
  call assert_all(io_files_ok, &
    "MPI EB AMR sparse checkpoint round-trip parity", rank)

  call read_sparse_owned_reactive_eb_topology_2d_checkpoint( &
    sparse_checkpoint_path, species, io_config, distribution, &
    coarse_geometry, restart_topology, root_owner, sparse_restart_set, &
    checkpoint_time, checkpoint_steps, checkpoint_regrids, &
    checkpoint_minimum_dt, checkpoint_base_density, ok, &
    local_restart_transfers)
  call MPI_Allreduce( &
    local_restart_transfers, global_restart_transfers, 1, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ok .and. ierr == MPI_SUCCESS .and. &
    checkpoint_time == 1.0e-6_dp .and. checkpoint_steps == 1 .and. &
    checkpoint_regrids == 0 .and. &
    checkpoint_minimum_dt == 1.0e-6_dp .and. &
    checkpoint_base_density == 0.31_dp .and. &
    local_restart_transfers == expected_local_restart_transfers .and. &
    global_restart_transfers == &
      expected_global_root_materialization_transfers .and. &
    sparse_restart_set%local_value_count() == &
      sparse_patch_set%local_value_count(), &
    "MPI EB AMR sparse checkpoint restart accounting", rank)
  ok = .true.
  do tile = 1, distribution%root_tile_count()
    if (.not. distribution%root_tile_is_local(tile)) cycle
    ok = ok .and. all(sparse_restart_set%root_tiles(tile)%state == &
        sparse_patch_set%root_tiles(tile)%state) .and. &
      all(sparse_restart_set%root_tiles(tile)%temperature == &
        sparse_patch_set%root_tiles(tile)%temperature)
  end do
  do child = 1, distribution%child_count()
    if (.not. distribution%child_is_local(child)) cycle
    ok = ok .and. all(sparse_restart_set%children(child)%state == &
        sparse_patch_set%children(child)%state) .and. &
      all(sparse_restart_set%children(child)%temperature == &
        sparse_patch_set%children(child)%temperature)
  end do
  call assert_all(ok, "MPI EB AMR sparse checkpoint restart parity", rank)
  io_files_ok = .true.
  if (rank == root_owner) then
    call remove_nonempty_file(sparse_checkpoint_path, io_files_ok)
  end if
  call assert_all(io_files_ok, &
    "MPI EB AMR sparse checkpoint restart cleanup", rank)

  call read_sparse_owned_reactive_eb_topology_2d_checkpoint( &
    missing_checkpoint_path, species, io_config, distribution, &
    coarse_geometry, restart_topology, root_owner, sparse_restart_set, &
    checkpoint_time, checkpoint_steps, checkpoint_regrids, &
    checkpoint_minimum_dt, checkpoint_base_density, ok, &
    local_restart_transfers)
  call assert_all(.not. ok .and. local_restart_transfers == 0 .and. &
    sparse_restart_set%local_value_count() == 0 .and. &
    checkpoint_time == 0.0_dp .and. checkpoint_steps == 0 .and. &
    checkpoint_regrids == 0 .and. checkpoint_minimum_dt == 0.0_dp .and. &
    checkpoint_base_density == 0.0_dp, &
    "MPI EB AMR sparse checkpoint read failure propagation", rank)

  call write_sparse_owned_reactive_eb_patch_set_2d_checkpoint( &
    missing_checkpoint_path, species, io_config, distribution, &
    sparse_patch_set, coarse_geometry, patch_set, root_owner, 1.0e-6_dp, &
    1, 0, 1.0e-6_dp, 0.31_dp, ok, &
    local_root_materialization_transfers)
  call assert_all(.not. ok .and. &
    local_root_materialization_transfers == 0, &
    "MPI EB AMR sparse checkpoint failure propagation", rank)

  call write_sparse_owned_reactive_eb_patch_set_2d_csv( &
    species, io_config, distribution, sparse_patch_set, coarse_geometry, &
    patch_set, root_owner, 1.0e-6_dp, ok, &
    local_root_materialization_transfers)
  call MPI_Allreduce( &
    local_root_materialization_transfers, &
    global_root_materialization_transfers, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ok .and. ierr == MPI_SUCCESS .and. &
    local_root_materialization_transfers == &
      expected_local_root_materialization_transfers .and. &
    global_root_materialization_transfers == &
      expected_global_root_materialization_transfers, &
    "MPI EB AMR sparse CSV transfer accounting", rank)
  io_files_ok = .true.
  if (rank == root_owner) then
    call remove_nonempty_file(sparse_root_output_path, io_files_ok)
    call remove_nonempty_file(sparse_fine_output_path_1, io_files_ok)
    call remove_nonempty_file(sparse_fine_output_path_2, io_files_ok)
  end if
  call assert_all(io_files_ok, "MPI EB AMR sparse CSV publication", rank)

  io_config%eb%flow%output_file = missing_root_output_path
  call write_sparse_owned_reactive_eb_patch_set_2d_csv( &
    species, io_config, distribution, sparse_patch_set, coarse_geometry, &
    patch_set, root_owner, 1.0e-6_dp, ok, &
    local_root_materialization_transfers)
  call assert_all(.not. ok .and. &
    local_root_materialization_transfers == 0, &
    "MPI EB AMR sparse CSV failure propagation", rank)
  io_config%eb%flow%output_file = sparse_root_output_path

  invalid_materialization_root = root_owner
  if (nranks == 1) then
    invalid_materialization_root = nranks
  else if (rank == nranks - 1) then
    invalid_materialization_root = modulo(root_owner + 1, nranks)
  end if
  call gather_sparse_owned_reactive_eb_patch_set_to_root_2d( &
    distribution, sparse_patch_set, coarse_geometry, patch_set, &
    invalid_materialization_root, root_materialized_state, &
    root_materialized_temperature, root_materialized_patch_set, ok, &
    local_root_materialization_transfers)
  call assert_all(.not. ok .and. &
    .not. allocated(root_materialized_state) .and. &
    .not. allocated(root_materialized_temperature) .and. &
    root_materialized_patch_set%patch_count() == 0 .and. &
    local_root_materialization_transfers == 0, &
    "MPI EB AMR inconsistent root-only materialization rollback", rank)

  invalid_sparse_patch_set = sparse_patch_set
  if (rank == 0) then
    do tile = 1, distribution%root_tile_count()
      if (.not. distribution%root_tile_is_local(tile)) cycle
      deallocate(invalid_sparse_patch_set%root_tiles(tile)%state)
      exit
    end do
  end if
  call gather_sparse_owned_reactive_eb_patch_set_to_root_2d( &
    distribution, invalid_sparse_patch_set, coarse_geometry, patch_set, &
    root_owner, root_materialized_state, root_materialized_temperature, &
    root_materialized_patch_set, ok, local_root_materialization_transfers)
  call assert_all(.not. ok .and. &
    .not. allocated(root_materialized_state) .and. &
    .not. allocated(root_materialized_temperature) .and. &
    root_materialized_patch_set%patch_count() == 0 .and. &
    local_root_materialization_transfers == 0, &
    "MPI EB AMR invalid root-only materialization rollback", rank)
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

  allocate(regrid_start_state, source=coarse_state)
  allocate(regrid_start_temperature, source=coarse_temperature)
  regrid_start_set = patch_set
  do child = 1, regrid_start_set%patch_count()
    regrid_start_set%children(child)%state = child_factor(child) * &
      regrid_start_set%children(child)%state
  end do
  topology_distribution = distribution
  sparse_regrid_template = regrid_start_set
  call scatter_owned_reactive_eb_patch_set_2d( &
    topology_distribution, size(species), regrid_start_state, &
    regrid_start_temperature, coarse_geometry, sparse_regrid_template, &
    sparse_regrid_set, ok)
  call assert_all(ok, "MPI EB AMR sparse regrid scatter", rank)

  regrid_tags = .false.
  regrid_tags(2:7, 5:11) = .true.
  call build_amr_eb_regrid_plan_collection_2d( &
    regrid_tags, criteria, regrid_collection, ok)
  call assert_all(ok .and. regrid_collection%patch_count() == 1, &
    "MPI EB AMR sparse regrid plan", rank)
  allocate(regrid_fine_geometries(regrid_collection%patch_count()))
  do child = 1, regrid_collection%patch_count()
    call build_patch_geometry( &
      coarse_geometry, regrid_collection%plans(child)%coarse_i_lower, &
      regrid_collection%plans(child)%coarse_i_upper, &
      regrid_collection%plans(child)%coarse_j_lower, &
      regrid_collection%plans(child)%coarse_j_upper, ratio, &
      regrid_fine_geometries(child), geometry_patch, ok)
    call assert_all(ok, "MPI EB AMR sparse regrid geometry", rank)
  end do
  allocate(regrid_reference_state, mold=coarse_state)
  allocate(regrid_reference_temperature, mold=coarse_temperature)
  call regrid_reactive_eb_patch_set_2d( &
    species, regrid_start_state, regrid_start_temperature, coarse_geometry, &
    regrid_start_set, regrid_fine_geometries, regrid_collection, ratio, &
    regrid_reference_state, regrid_reference_temperature, &
    regrid_reference_set, ok)
  call assert_all(ok .and. regrid_reference_set%patch_count() == 1, &
    "serial EB AMR sparse regrid reference", rank)

  call regrid_sparse_owned_reactive_eb_patch_set_2d( &
    species, topology_distribution, sparse_regrid_set, coarse_geometry, &
    sparse_regrid_template, regrid_fine_geometries, regrid_collection, 1, &
    ok, topology_changed, local_regrid_restriction_transfers, &
    local_regrid_prolongation_transfers, local_regrid_overlap_transfers)
  call assert_all(.not. ok .and. .not. topology_changed .and. &
    local_regrid_restriction_transfers == 0 .and. &
    local_regrid_prolongation_transfers == 0 .and. &
    local_regrid_overlap_transfers == 0 .and. &
    topology_distribution%child_count() == distribution%child_count() .and. &
    all(topology_distribution%child_owners == &
      distribution%child_owners) .and. &
    sparse_regrid_template%patch_count() == patch_set%patch_count() .and. &
    sparse_regrid_set%is_valid( &
      topology_distribution, coarse_geometry, sparse_regrid_template), &
    "MPI EB AMR sparse regrid control rollback", rank)
  call materialize_owned_reactive_eb_patch_set_2d( &
    topology_distribution, sparse_regrid_set, coarse_state, &
    coarse_temperature, coarse_geometry, sparse_regrid_template, &
    materialized_coarse_state, materialized_coarse_temperature, &
    materialized_patch_set, ok)
  call assert_all(ok .and. &
    all(materialized_coarse_state == regrid_start_state) .and. &
    all(materialized_coarse_temperature == regrid_start_temperature), &
    "MPI EB AMR sparse regrid root rollback", rank)
  do child = 1, regrid_start_set%patch_count()
    call assert_all( &
      all(materialized_patch_set%children(child)%state == &
        regrid_start_set%children(child)%state) .and. &
      all(materialized_patch_set%children(child)%temperature == &
        regrid_start_set%children(child)%temperature), &
      "MPI EB AMR sparse regrid child rollback", rank)
  end do

  allocate(rejected_regrid_fine_geometries, source=regrid_fine_geometries)
  geometry_perturbed = .false.
  do new_child = 1, size(rejected_regrid_fine_geometries)
    do j = 1, rejected_regrid_fine_geometries(new_child)%ny
      global_j = (regrid_collection%plans(new_child)%coarse_j_lower - 1) * &
        ratio + j
      do i = 1, rejected_regrid_fine_geometries(new_child)%nx
        if (rejected_regrid_fine_geometries(new_child)%cell_type(i, j) /= &
            eb_cut_cell) cycle
        global_i = (regrid_collection%plans(new_child)%coarse_i_lower - 1) * &
          ratio + i
        do old_child = 1, regrid_start_set%patch_count()
          old_i = global_i - &
            (regrid_start_set%children(old_child)%patch%coarse_i_lower - 1) * &
              ratio
          old_j = global_j - &
            (regrid_start_set%children(old_child)%patch%coarse_j_lower - 1) * &
              ratio
          if (old_i < 1 .or. &
              old_i > regrid_start_set%children(old_child)%geometry%nx .or. &
              old_j < 1 .or. &
              old_j > regrid_start_set%children(old_child)%geometry%ny) cycle
          if (rejected_regrid_fine_geometries(new_child)% &
              volume_fraction(i, j) <= 0.5_dp) then
            rejected_regrid_fine_geometries(new_child)% &
              volume_fraction(i, j) = &
                rejected_regrid_fine_geometries(new_child)% &
                  volume_fraction(i, j) + 1.0e-8_dp
          else
            rejected_regrid_fine_geometries(new_child)% &
              volume_fraction(i, j) = &
                rejected_regrid_fine_geometries(new_child)% &
                  volume_fraction(i, j) - 1.0e-8_dp
          end if
          geometry_perturbed = .true.
          exit
        end do
        if (geometry_perturbed) exit
      end do
      if (geometry_perturbed) exit
    end do
    if (geometry_perturbed) exit
  end do
  call assert_all(geometry_perturbed .and. &
    rejected_regrid_fine_geometries(1)%is_valid(), &
    "MPI EB AMR sparse regrid mismatched geometry setup", rank)
  call regrid_sparse_owned_reactive_eb_patch_set_2d( &
    species, topology_distribution, sparse_regrid_set, coarse_geometry, &
    sparse_regrid_template, rejected_regrid_fine_geometries, &
    regrid_collection, ratio, ok, topology_changed, &
    local_regrid_restriction_transfers, &
    local_regrid_prolongation_transfers, local_regrid_overlap_transfers)
  call assert_all(.not. ok .and. .not. topology_changed .and. &
    local_regrid_restriction_transfers == 0 .and. &
    local_regrid_prolongation_transfers == 0 .and. &
    local_regrid_overlap_transfers == 0 .and. &
    topology_distribution%child_count() == distribution%child_count() .and. &
    all(topology_distribution%child_owners == distribution%child_owners) .and. &
    same_patch_topology(sparse_regrid_template, regrid_start_set), &
    "MPI EB AMR sparse regrid geometry rollback", rank)
  call materialize_owned_reactive_eb_patch_set_2d( &
    topology_distribution, sparse_regrid_set, coarse_state, &
    coarse_temperature, coarse_geometry, sparse_regrid_template, &
    materialized_coarse_state, materialized_coarse_temperature, &
    materialized_patch_set, ok)
  call assert_all(ok .and. &
    all(materialized_coarse_state == regrid_start_state) .and. &
    all(materialized_coarse_temperature == regrid_start_temperature), &
    "MPI EB AMR sparse regrid geometry root rollback", rank)
  do child = 1, regrid_start_set%patch_count()
    call assert_all( &
      all(materialized_patch_set%children(child)%state == &
        regrid_start_set%children(child)%state) .and. &
      all(materialized_patch_set%children(child)%temperature == &
        regrid_start_set%children(child)%temperature), &
      "MPI EB AMR sparse regrid geometry child rollback", rank)
  end do

  call regrid_sparse_owned_reactive_eb_patch_set_2d( &
    species, topology_distribution, sparse_regrid_set, coarse_geometry, &
    sparse_regrid_template, regrid_fine_geometries, regrid_collection, &
    ratio, ok, topology_changed, local_regrid_restriction_transfers, &
    local_regrid_prolongation_transfers, local_regrid_overlap_transfers)
  call assert_all(ok .and. topology_changed .and. &
    topology_distribution%child_count() == 1 .and. &
    sparse_regrid_template%patch_count() == 1 .and. &
    topology_distribution%subcycle_exponent == &
      distribution%subcycle_exponent .and. &
    sum(topology_distribution%rank_entity_counts) == &
      topology_distribution%root_tile_count() + 1 .and. &
    sum(topology_distribution%rank_cell_counts) == coarse_nx * coarse_ny + &
      regrid_reference_set%children(1)%geometry%nx * &
        regrid_reference_set%children(1)%geometry%ny .and. &
    sparse_regrid_set%is_valid( &
      topology_distribution, coarse_geometry, sparse_regrid_template), &
    "MPI EB AMR sparse regrid ownership rebuild", rank)
  sparse_local_values = int(sparse_regrid_set%local_value_count())
  call MPI_Allreduce( &
    sparse_local_values, sparse_global_values, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    sparse_local_values == (nvar + 1) * &
      topology_distribution%rank_cell_counts(rank + 1) .and. &
    sparse_global_values == (nvar + 1) * &
      sum(topology_distribution%rank_cell_counts), &
    "MPI EB AMR sparse regrid one-copy storage", rank)

  allocate(regrid_recipients(nranks))
  expected_local_regrid_restriction_transfers = 0
  expected_global_regrid_restriction_transfers = 0
  do child = 1, distribution%child_count()
    regrid_recipients = .false.
    do tile = 1, distribution%root_tile_count()
      if (distribution%root_tiles(tile)%j_upper < &
          regrid_start_set%children(child)%patch%coarse_j_lower .or. &
          distribution%root_tiles(tile)%j_lower > &
          regrid_start_set%children(child)%patch%coarse_j_upper) cycle
      owner = distribution%root_tiles(tile)%owner
      regrid_recipients(owner + 1) = .true.
    end do
    owner = distribution%child_owner(child)
    regrid_recipients(owner + 1) = .false.
    expected_global_regrid_restriction_transfers = &
      expected_global_regrid_restriction_transfers + &
      count(regrid_recipients)
    if (rank == owner) expected_local_regrid_restriction_transfers = &
      expected_local_regrid_restriction_transfers + &
      count(regrid_recipients)
  end do
  regrid_recipients = .false.
  do child = 1, topology_distribution%child_count()
    owner = topology_distribution%child_owner(child)
    regrid_recipients(owner + 1) = .true.
  end do
  expected_local_regrid_prolongation_transfers = 0
  expected_global_regrid_prolongation_transfers = 0
  do owner = 0, nranks - 1
    if (.not. regrid_recipients(owner + 1)) cycle
    do tile = 1, distribution%root_tile_count()
      if (distribution%root_tiles(tile)%owner == owner) cycle
      expected_global_regrid_prolongation_transfers = &
        expected_global_regrid_prolongation_transfers + 1
      if (rank == distribution%root_tiles(tile)%owner) &
        expected_local_regrid_prolongation_transfers = &
          expected_local_regrid_prolongation_transfers + 1
    end do
  end do
  expected_local_regrid_overlap_transfers = 0
  expected_global_regrid_overlap_transfers = 0
  do new_child = 1, sparse_regrid_template%patch_count()
    do old_child = 1, regrid_start_set%patch_count()
      if (regrid_start_set%children(old_child)%patch%refinement_ratio /= &
          sparse_regrid_template%children(new_child)%patch% &
            refinement_ratio) cycle
      if (max( &
          regrid_start_set%children(old_child)%patch%coarse_i_lower, &
          sparse_regrid_template%children(new_child)%patch%coarse_i_lower) > &
          min( &
            regrid_start_set%children(old_child)%patch%coarse_i_upper, &
            sparse_regrid_template%children(new_child)%patch%coarse_i_upper) &
          .or. max( &
            regrid_start_set%children(old_child)%patch%coarse_j_lower, &
            sparse_regrid_template%children(new_child)%patch%coarse_j_lower) > &
          min( &
            regrid_start_set%children(old_child)%patch%coarse_j_upper, &
            sparse_regrid_template%children(new_child)%patch%coarse_j_upper)) &
        cycle
      if (distribution%child_owner(old_child) == &
          topology_distribution%child_owner(new_child)) cycle
      expected_global_regrid_overlap_transfers = &
        expected_global_regrid_overlap_transfers + 1
      if (rank == distribution%child_owner(old_child)) &
        expected_local_regrid_overlap_transfers = &
          expected_local_regrid_overlap_transfers + 1
    end do
  end do
  call MPI_Allreduce( &
    local_regrid_restriction_transfers, &
    global_regrid_restriction_transfers, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_regrid_restriction_transfers == &
      expected_local_regrid_restriction_transfers .and. &
    global_regrid_restriction_transfers == &
      expected_global_regrid_restriction_transfers, &
    "MPI EB AMR sparse regrid restriction traffic", rank)
  call MPI_Allreduce( &
    local_regrid_prolongation_transfers, &
    global_regrid_prolongation_transfers, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_regrid_prolongation_transfers == &
      expected_local_regrid_prolongation_transfers .and. &
    global_regrid_prolongation_transfers == &
      expected_global_regrid_prolongation_transfers, &
    "MPI EB AMR sparse regrid prolongation traffic", rank)
  call MPI_Allreduce( &
    local_regrid_overlap_transfers, global_regrid_overlap_transfers, 1, &
    MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_regrid_overlap_transfers == &
      expected_local_regrid_overlap_transfers .and. &
    global_regrid_overlap_transfers == &
      expected_global_regrid_overlap_transfers, &
    "MPI EB AMR sparse regrid overlap traffic", rank)
  call materialize_owned_reactive_eb_patch_set_2d( &
    topology_distribution, sparse_regrid_set, coarse_state, &
    coarse_temperature, coarse_geometry, sparse_regrid_template, &
    materialized_coarse_state, materialized_coarse_temperature, &
    materialized_patch_set, ok)
  call assert_all(ok .and. &
    all(materialized_coarse_state == regrid_reference_state) .and. &
    all(materialized_coarse_temperature == regrid_reference_temperature), &
    "MPI EB AMR sparse regrid serial root parity", rank)
  call assert_all( &
    all(materialized_patch_set%children(1)%state == &
      regrid_reference_set%children(1)%state) .and. &
    all(materialized_patch_set%children(1)%temperature == &
      regrid_reference_set%children(1)%temperature), &
    "MPI EB AMR sparse regrid serial child parity", rank)

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
  expected_local_root_hydro_cells = 0
  do tile = 1, distribution%root_tile_count()
    if (.not. distribution%root_tile_is_local(tile)) cycle
    expected_local_root_hydro_cells = expected_local_root_hydro_cells + &
      coarse_nx * ( &
        min(coarse_ny, distribution%root_tiles(tile)%j_upper + &
          mpi_amr_eb_root_tile_hydro_halo_cells) - &
        max(1, distribution%root_tiles(tile)%j_lower - &
          mpi_amr_eb_root_tile_hydro_halo_cells) + 1)
  end do
  call advance_owned_reactive_eb_patch_set_hydro_2d( &
    species, distribution, hydro_mpi_state, hydro_mpi_temperature, &
    coarse_geometry, hydro_mpi_set, "hllc", "pcm", "mc", 2, hydro_dt, &
    ok, local_advances, 0.5_dp, local_root_hydro_cells)
  call assert_all(ok, "MPI owner-only EB AMR hydro", rank)
  expected_local_advances = 0
  do tile = 1, distribution%root_tile_count()
    if (distribution%root_tile_is_local(tile)) &
      expected_local_advances = expected_local_advances + 1
  end do
  expected_global_advances = distribution%root_tile_count()
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
  call MPI_Allreduce( &
    local_root_hydro_cells, global_root_hydro_cells, 1, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce( &
    expected_local_root_hydro_cells, expected_global_root_hydro_cells, 1, &
    MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_root_hydro_cells == expected_local_root_hydro_cells .and. &
    global_root_hydro_cells == expected_global_root_hydro_cells .and. &
    (distribution%nranks == 1 .or. global_root_hydro_cells < &
      distribution%nranks * coarse_nx * coarse_ny), &
    "MPI EB AMR bounded root-tile hydro work", rank)
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
  expected_local_root_hydro_cells = 0
  restriction_recipients = .false.
  do tile = 1, distribution%root_tile_count()
    owner = distribution%root_tiles(tile)%owner
    if (rank == owner) expected_local_root_hydro_cells = &
      expected_local_root_hydro_cells + coarse_nx * ( &
        min(coarse_ny, distribution%root_tiles(tile)%j_upper + &
          mpi_amr_eb_root_tile_hydro_halo_cells) - &
        max(1, distribution%root_tiles(tile)%j_lower - &
          mpi_amr_eb_root_tile_hydro_halo_cells) + 1)
    do source = 1, distribution%root_tile_count()
      if (distribution%root_tiles(source)%owner == owner) cycle
      if (min( &
          coarse_ny, distribution%root_tiles(tile)%j_upper + &
            mpi_amr_eb_root_tile_hydro_halo_cells) < &
          distribution%root_tiles(source)%j_lower .or. &
          max(1, distribution%root_tiles(tile)%j_lower - &
            mpi_amr_eb_root_tile_hydro_halo_cells) > &
          distribution%root_tiles(source)%j_upper) cycle
      expected_global_root_transfers = expected_global_root_transfers + 1
      if (rank == distribution%root_tiles(source)%owner) &
        expected_local_root_transfers = expected_local_root_transfers + 1
    end do
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
  expected_local_advances = 0
  do tile = 1, distribution%root_tile_count()
    if (distribution%root_tile_is_local(tile)) &
      expected_local_advances = expected_local_advances + 1
  end do
  expected_global_advances = distribution%root_tile_count()
  do child = 1, distribution%child_count()
    expected_global_advances = expected_global_advances + &
      hydro_start_set%children(child)%patch%refinement_ratio
    if (distribution%child_is_local(child)) &
      expected_local_advances = expected_local_advances + &
        hydro_start_set%children(child)%patch%refinement_ratio
  end do
  call advance_sparse_owned_reactive_eb_patch_set_hydro_2d( &
    species, distribution, sparse_hydro_set, coarse_geometry, &
    hydro_start_set, "hllc", "pcm", "mc", 2, hydro_dt, ok, &
    local_advances, 0.5_dp, local_root_transfers, local_root_hydro_cells)
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
  call MPI_Allreduce( &
    local_root_hydro_cells, global_root_hydro_cells, 1, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce( &
    expected_local_root_hydro_cells, expected_global_root_hydro_cells, 1, &
    MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_root_hydro_cells == expected_local_root_hydro_cells .and. &
    global_root_hydro_cells == expected_global_root_hydro_cells .and. &
    (distribution%nranks == 1 .or. global_root_hydro_cells < &
      distribution%nranks * coarse_nx * coarse_ny), &
    "MPI EB AMR bounded sparse root-tile hydro work", rank)
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
    ok, local_advances, 0.5_dp, local_root_hydro_cells)
  call assert_all(.not. ok .and. local_advances == 0 .and. &
    local_root_hydro_cells == 0 .and. &
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
    local_advances, 0.5_dp, local_root_transfers, local_root_hydro_cells)
  call assert_all(.not. ok .and. local_advances == 0 .and. &
    local_root_transfers == 0 .and. local_root_hydro_cells == 0 .and. &
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
    "MPI EB AMR owner-local sparse timestep zero traffic", rank)

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
  expected_local_root_transport_cells = 0
  restriction_recipients = .false.
  do tile = 1, distribution%root_tile_count()
    owner = distribution%root_tiles(tile)%owner
    if (rank == owner) expected_local_root_transport_cells = &
      expected_local_root_transport_cells + 2 * coarse_nx * ( &
        min(coarse_ny, distribution%root_tiles(tile)%j_upper + &
          mpi_amr_eb_root_tile_transport_halo_cells) - &
        max(1, distribution%root_tiles(tile)%j_lower - &
          mpi_amr_eb_root_tile_transport_halo_cells) + 1)
    do source = 1, distribution%root_tile_count()
      if (distribution%root_tiles(source)%owner == owner) cycle
      if (min( &
          coarse_ny, distribution%root_tiles(tile)%j_upper + &
            mpi_amr_eb_root_tile_transport_halo_cells) < &
          distribution%root_tiles(source)%j_lower .or. &
          max(1, distribution%root_tiles(tile)%j_lower - &
            mpi_amr_eb_root_tile_transport_halo_cells) > &
          distribution%root_tiles(source)%j_upper) cycle
      expected_global_root_transfers = expected_global_root_transfers + 2
      if (rank == distribution%root_tiles(source)%owner) &
        expected_local_root_transfers = expected_local_root_transfers + 2
    end do
    if (owner == root_owner) cycle
    expected_global_root_transfers = expected_global_root_transfers + 4
    if (rank == owner) expected_local_root_transfers = &
      expected_local_root_transfers + 2
    if (rank == root_owner) expected_local_root_transfers = &
      expected_local_root_transfers + 2
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
  expected_local_advances = 0
  do tile = 1, distribution%root_tile_count()
    if (distribution%root_tile_is_local(tile)) &
      expected_local_advances = expected_local_advances + 2
  end do
  expected_global_advances = 2 * distribution%root_tile_count()
  do child = 1, distribution%child_count()
    expected_global_advances = expected_global_advances + 2 * &
      transport_start_set%children(child)%patch%refinement_ratio
    if (distribution%child_is_local(child)) &
      expected_local_advances = expected_local_advances + 2 * &
        transport_start_set%children(child)%patch%refinement_ratio
  end do
  call advance_sparse_owned_reactive_eb_patch_set_transport_2d( &
    species, transport, distribution, sparse_transport_set, &
    coarse_geometry, transport_start_set, transport_dt, .true., .true., &
    .true., .true., boundaries, 2, ok, local_advances, transport_theta, &
    0.5_dp, local_root_transfers, local_root_transport_cells)
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
  call MPI_Allreduce( &
    local_root_transport_cells, global_root_transport_cells, 1, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  expected_global_root_transport_cells = 0
  call MPI_Allreduce( &
    expected_local_root_transport_cells, &
    expected_global_root_transport_cells, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_root_transport_cells == expected_local_root_transport_cells .and. &
    global_root_transport_cells == expected_global_root_transport_cells .and. &
    (distribution%nranks <= 2 .or. global_root_transport_cells < &
      2 * distribution%nranks * coarse_nx * coarse_ny), &
    "MPI EB AMR sparse root transport bounded work", rank)
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
    0.5_dp, local_root_transfers, local_root_transport_cells)
  call assert_all(.not. ok .and. local_advances == 0 .and. &
    transport_theta == 1.0_dp .and. &
    local_root_transfers == 0 .and. local_root_transport_cells == 0 .and. &
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
  do tile = 1, distribution%root_tile_count()
    if (distribution%root_tile_is_local(tile)) &
      expected_local_hydro = expected_local_hydro + 1
  end do
  if (rank == distribution%root_level_owner()) then
    expected_local_transport = 4
  end if
  expected_global_hydro = distribution%root_tile_count()
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

  expected_local_hydro = 0
  expected_local_transport = 0
  do tile = 1, distribution%root_tile_count()
    if (distribution%root_tile_is_local(tile)) then
      expected_local_hydro = expected_local_hydro + 1
      expected_local_transport = expected_local_transport + 4
    end if
  end do
  expected_global_hydro = distribution%root_tile_count()
  expected_global_transport = 4 * distribution%root_tile_count()
  do child = 1, distribution%child_count()
    expected_global_hydro = expected_global_hydro + &
      patch_set%children(child)%patch%refinement_ratio
    if (distribution%child_is_local(child)) &
      expected_local_hydro = expected_local_hydro + &
        patch_set%children(child)%patch%refinement_ratio
    expected_global_transport = expected_global_transport + 4 * &
      patch_set%children(child)%patch%refinement_ratio
    if (distribution%child_is_local(child)) &
      expected_local_transport = expected_local_transport + 4 * &
        patch_set%children(child)%patch%refinement_ratio
  end do
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

  call compute_serial_full_physics_timestep( &
    coarse_state, coarse_temperature, patch_set, time_loop_initial_dt, ok)
  call assert_all(ok, "serial EB AMR time-loop initial timestep", rank)
  time_loop_final_time = 1.25_dp * time_loop_initial_dt
  allocate(time_loop_reference_state, mold=coarse_state)
  allocate(time_loop_reference_temperature, mold=coarse_temperature)
  call advance_serial_full_physics_to_time( &
    coarse_state, coarse_temperature, patch_set, time_loop_final_time, &
    time_loop_reference_state, time_loop_reference_temperature, &
    time_loop_reference_set, time_loop_reference_steps, &
    time_loop_reference_minimum_dt, time_loop_reference_theta, ok)
  call assert_all(ok .and. time_loop_reference_steps >= 2, &
    "serial EB AMR clipped multi-step reference", rank)

  call scatter_owned_reactive_eb_patch_set_2d( &
    distribution, size(species), coarse_state, coarse_temperature, &
    coarse_geometry, patch_set, sparse_time_loop_set, ok)
  call assert_all(ok, "MPI EB AMR sparse time-loop scatter", rank)
  time_loop_time = 0.0_dp
  time_loop_steps = 0
  call advance_sparse_owned_reactive_eb_patch_set_to_time_2d( &
    species, reactions, transport, distribution, sparse_time_loop_set, &
    coarse_geometry, patch_set, "hllc", "pcm", "mc", 2, time_loop_time, &
    time_loop_final_time, time_loop_steps, 16, 0.35_dp, 0.20_dp, &
    1.0e-8_dp, 1.0e-14_dp, .true., .true., .true., .true., boundaries, ok, &
    time_loop_minimum_dt, time_loop_advanced_steps, &
    local_chemistry_advances, local_hydro_advances, &
    local_transport_advances, full_theta, 0.5_dp, local_root_transfers)
  call assert_all(ok .and. time_loop_time == time_loop_final_time .and. &
    time_loop_steps == time_loop_reference_steps .and. &
    time_loop_advanced_steps == time_loop_reference_steps .and. &
    abs(time_loop_minimum_dt - time_loop_reference_minimum_dt) <= &
      128.0_dp * epsilon(1.0_dp) * time_loop_reference_minimum_dt .and. &
    abs(full_theta - time_loop_reference_theta) <= &
      5.0e-13_dp * max(1.0_dp, abs(time_loop_reference_theta)), &
    "MPI EB AMR sparse public time-loop control", rank)
  call assert_all( &
    local_chemistry_advances == &
      time_loop_steps * expected_local_chemistry .and. &
    local_hydro_advances == time_loop_steps * expected_local_hydro .and. &
    local_transport_advances == &
      time_loop_steps * expected_local_transport, &
    "MPI EB AMR sparse time-loop local accounting", rank)
  call MPI_Allreduce( &
    local_chemistry_advances, global_chemistry_advances, 1, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS, &
    "MPI EB AMR sparse time-loop chemistry reduction", rank)
  call MPI_Allreduce( &
    local_hydro_advances, global_hydro_advances, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS, &
    "MPI EB AMR sparse time-loop hydro reduction", rank)
  call MPI_Allreduce( &
    local_transport_advances, global_transport_advances, 1, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    global_chemistry_advances == &
      time_loop_steps * expected_global_chemistry .and. &
    global_hydro_advances == time_loop_steps * expected_global_hydro .and. &
    global_transport_advances == &
      time_loop_steps * expected_global_transport, &
    "MPI EB AMR sparse time-loop global accounting", rank)
  expected_local_root_transfers = 0
  expected_global_root_transfers = 0
  call MPI_Allreduce( &
    local_root_transfers, global_root_transfers, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_root_transfers == &
      time_loop_steps * expected_local_root_transfers .and. &
    global_root_transfers == &
      time_loop_steps * expected_global_root_transfers, &
    "MPI EB AMR sparse time-loop timestep traffic", rank)
  call materialize_owned_reactive_eb_patch_set_2d( &
    distribution, sparse_time_loop_set, coarse_state, coarse_temperature, &
    coarse_geometry, patch_set, materialized_coarse_state, &
    materialized_coarse_temperature, materialized_patch_set, ok)
  full_scale = max(1.0_dp, maxval(abs(time_loop_reference_state)))
  call assert_all(ok .and. &
    maxval(abs(materialized_coarse_state - time_loop_reference_state)) <= &
      2.0e-10_dp * full_scale .and. &
    maxval(abs(materialized_coarse_temperature - &
      time_loop_reference_temperature)) <= 2.0e-10_dp * &
        max(1.0_dp, maxval(abs(time_loop_reference_temperature))), &
    "MPI EB AMR sparse time-loop serial root parity", rank)
  do child = 1, materialized_patch_set%patch_count()
    call assert_all( &
      maxval(abs(materialized_patch_set%children(child)%state - &
        time_loop_reference_set%children(child)%state)) <= &
          2.0e-10_dp * full_scale .and. &
      maxval(abs(materialized_patch_set%children(child)%temperature - &
        time_loop_reference_set%children(child)%temperature)) <= &
          2.0e-10_dp * max(1.0_dp, maxval(abs( &
            time_loop_reference_set%children(child)%temperature))), &
      "MPI EB AMR sparse time-loop serial child parity", rank)
  end do

  allocate(scheduled_start_state, mold=coarse_state)
  allocate(scheduled_start_temperature, mold=coarse_temperature)
  call initialize_scheduled_regrid_state( &
    scheduled_start_state, scheduled_start_temperature, ok)
  call assert_all(ok, "serial EB AMR scheduled-regrid state", rank)
  call initialize_reactive_eb_patch_set_2d( &
    species, scheduled_start_state, scheduled_start_temperature, &
    coarse_geometry, fine_geometries, collection, ratio, &
    scheduled_start_set, ok)
  call assert_all(ok .and. scheduled_start_set%patch_count() == 2, &
    "serial EB AMR scheduled-regrid hierarchy", rank)
  scheduled_criteria%relative_gradient_threshold = 1.0e-4_dp
  scheduled_criteria%absolute_gradient_threshold = 1.0e-3_dp
  scheduled_criteria%scale_floor = 1.0_dp
  scheduled_criteria%buffer_cells = 0
  scheduled_criteria%minimum_patch_cells_x = 5
  scheduled_criteria%minimum_patch_cells_y = 5
  scheduled_criteria%maximum_patch_gap_cells = 0
  call plan_reactive_eb_temperature_regrid_collection_2d( &
    scheduled_start_temperature, coarse_geometry, scheduled_criteria, &
    regrid_tags, scheduled_collection, ok)
  call assert_all(ok .and. scheduled_collection%patch_count() >= 1, &
    "serial EB AMR scheduled-regrid temperature tags", rank)
  call compute_serial_full_physics_timestep( &
    scheduled_start_state, scheduled_start_temperature, scheduled_start_set, &
    scheduled_initial_dt, ok)
  call assert_all(ok, "serial EB AMR scheduled-regrid timestep", rank)
  scheduled_final_time = 1.01_dp * scheduled_initial_dt
  allocate(scheduled_reference_state, mold=coarse_state)
  allocate(scheduled_reference_temperature, mold=coarse_temperature)
  call advance_serial_full_physics_to_time_with_regrid( &
    scheduled_start_state, scheduled_start_temperature, &
    scheduled_start_set, scheduled_final_time, 2, scheduled_criteria, &
    scheduled_reference_state, scheduled_reference_temperature, &
    scheduled_reference_set, scheduled_reference_steps, &
    scheduled_reference_minimum_dt, scheduled_reference_theta, &
    scheduled_reference_evaluations, scheduled_reference_regrids, ok)
  call assert_all(ok .and. scheduled_reference_steps == 2 .and. &
    scheduled_reference_evaluations == 1 .and. &
    scheduled_reference_regrids == 1, &
    "serial EB AMR scheduled-regrid reference", rank)

  scheduled_distribution = distribution
  scheduled_template = scheduled_start_set
  call scatter_owned_reactive_eb_patch_set_2d( &
    scheduled_distribution, size(species), scheduled_start_state, &
    scheduled_start_temperature, coarse_geometry, scheduled_template, &
    sparse_scheduled_set, ok)
  call assert_all(ok, "MPI EB AMR scheduled-regrid scatter", rank)
  scheduled_time = 0.0_dp
  scheduled_steps = 0
  scheduled_regrid_evaluations = 0
  scheduled_regrids = 0
  call advance_sparse_owned_reactive_eb_patch_set_to_time_2d( &
    species, reactions, transport, scheduled_distribution, &
    sparse_scheduled_set, coarse_geometry, scheduled_template, "hllc", &
    "pcm", "mc", 2, scheduled_time, scheduled_final_time, scheduled_steps, &
    16, 0.35_dp, 0.20_dp, 1.0e-8_dp, 1.0e-14_dp, .true., .true., .true., &
    .true., boundaries, ok, scheduled_minimum_dt, &
    advanced_steps=scheduled_advanced_steps, &
    minimum_transport_theta=scheduled_theta, &
    local_timestep_root_transfers=scheduled_timestep_transfers, &
    regrid_evaluations=scheduled_regrid_evaluations, &
    regrids=scheduled_regrids, regrid_interval=2, &
    regrid_criteria=scheduled_criteria, refinement_ratio=ratio, &
    geometry_builder=build_scheduled_patch_geometry, &
    local_regrid_root_transfers=scheduled_regrid_transfers, &
    local_regrid_restriction_transfers= &
      local_regrid_restriction_transfers, &
    local_regrid_prolongation_transfers= &
      local_regrid_prolongation_transfers, &
    local_regrid_overlap_transfers=local_regrid_overlap_transfers)
  call assert_all(ok .and. scheduled_time == scheduled_final_time .and. &
    scheduled_steps == scheduled_reference_steps .and. &
    scheduled_advanced_steps == scheduled_reference_steps .and. &
    scheduled_regrid_evaluations == scheduled_reference_evaluations .and. &
    scheduled_regrids == scheduled_reference_regrids .and. &
    abs(scheduled_minimum_dt - scheduled_reference_minimum_dt) <= &
      128.0_dp * epsilon(1.0_dp) * scheduled_reference_minimum_dt .and. &
    abs(scheduled_theta - scheduled_reference_theta) <= &
      5.0e-13_dp * max(1.0_dp, abs(scheduled_reference_theta)) .and. &
    same_patch_topology(scheduled_template, scheduled_reference_set), &
    "MPI EB AMR scheduled-regrid public control", rank)
  expected_local_regrid_restriction_transfers = 0
  expected_global_regrid_restriction_transfers = 0
  do child = 1, distribution%child_count()
    regrid_recipients = .false.
    do tile = 1, distribution%root_tile_count()
      if (distribution%root_tiles(tile)%j_upper < &
          scheduled_start_set%children(child)%patch%coarse_j_lower .or. &
          distribution%root_tiles(tile)%j_lower > &
          scheduled_start_set%children(child)%patch%coarse_j_upper) cycle
      owner = distribution%root_tiles(tile)%owner
      regrid_recipients(owner + 1) = .true.
    end do
    owner = distribution%child_owner(child)
    regrid_recipients(owner + 1) = .false.
    expected_global_regrid_restriction_transfers = &
      expected_global_regrid_restriction_transfers + &
      count(regrid_recipients)
    if (rank == owner) expected_local_regrid_restriction_transfers = &
      expected_local_regrid_restriction_transfers + &
      count(regrid_recipients)
  end do
  regrid_recipients = .false.
  do child = 1, scheduled_distribution%child_count()
    owner = scheduled_distribution%child_owner(child)
    regrid_recipients(owner + 1) = .true.
  end do
  expected_local_regrid_prolongation_transfers = 0
  expected_global_regrid_prolongation_transfers = 0
  do owner = 0, nranks - 1
    if (.not. regrid_recipients(owner + 1)) cycle
    do tile = 1, distribution%root_tile_count()
      if (distribution%root_tiles(tile)%owner == owner) cycle
      expected_global_regrid_prolongation_transfers = &
        expected_global_regrid_prolongation_transfers + 1
      if (rank == distribution%root_tiles(tile)%owner) &
        expected_local_regrid_prolongation_transfers = &
          expected_local_regrid_prolongation_transfers + 1
    end do
  end do
  expected_local_regrid_overlap_transfers = 0
  expected_global_regrid_overlap_transfers = 0
  do new_child = 1, scheduled_template%patch_count()
    do old_child = 1, scheduled_start_set%patch_count()
      if (scheduled_start_set%children(old_child)%patch%refinement_ratio /= &
          scheduled_template%children(new_child)%patch%refinement_ratio) cycle
      if (max( &
          scheduled_start_set%children(old_child)%patch%coarse_i_lower, &
          scheduled_template%children(new_child)%patch%coarse_i_lower) > &
          min( &
            scheduled_start_set%children(old_child)%patch%coarse_i_upper, &
            scheduled_template%children(new_child)%patch%coarse_i_upper) &
          .or. max( &
            scheduled_start_set%children(old_child)%patch%coarse_j_lower, &
            scheduled_template%children(new_child)%patch%coarse_j_lower) > &
          min( &
            scheduled_start_set%children(old_child)%patch%coarse_j_upper, &
            scheduled_template%children(new_child)%patch%coarse_j_upper)) &
        cycle
      if (distribution%child_owner(old_child) == &
          scheduled_distribution%child_owner(new_child)) cycle
      expected_global_regrid_overlap_transfers = &
        expected_global_regrid_overlap_transfers + 1
      if (rank == distribution%child_owner(old_child)) &
        expected_local_regrid_overlap_transfers = &
          expected_local_regrid_overlap_transfers + 1
    end do
  end do
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
  call MPI_Allreduce( &
    scheduled_regrid_transfers, global_scheduled_regrid_transfers, 1, &
    MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    scheduled_timestep_transfers == 0 .and. &
    scheduled_regrid_transfers == expected_local_root_transfers .and. &
    global_scheduled_regrid_transfers == expected_global_root_transfers, &
    "MPI EB AMR scheduled-regrid root traffic", rank)
  call MPI_Allreduce( &
    local_regrid_restriction_transfers, &
    global_regrid_restriction_transfers, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_regrid_restriction_transfers == &
      expected_local_regrid_restriction_transfers .and. &
    global_regrid_restriction_transfers == &
      expected_global_regrid_restriction_transfers, &
    "MPI EB AMR scheduled-regrid restriction traffic", rank)
  call MPI_Allreduce( &
    local_regrid_prolongation_transfers, &
    global_regrid_prolongation_transfers, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_regrid_prolongation_transfers == &
      expected_local_regrid_prolongation_transfers .and. &
    global_regrid_prolongation_transfers == &
      expected_global_regrid_prolongation_transfers, &
    "MPI EB AMR scheduled-regrid prolongation traffic", rank)
  call MPI_Allreduce( &
    local_regrid_overlap_transfers, global_regrid_overlap_transfers, 1, &
    MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    local_regrid_overlap_transfers == &
      expected_local_regrid_overlap_transfers .and. &
    global_regrid_overlap_transfers == &
      expected_global_regrid_overlap_transfers, &
    "MPI EB AMR scheduled-regrid overlap traffic", rank)
  call materialize_owned_reactive_eb_patch_set_2d( &
    scheduled_distribution, sparse_scheduled_set, scheduled_start_state, &
    scheduled_start_temperature, coarse_geometry, scheduled_template, &
    materialized_coarse_state, materialized_coarse_temperature, &
    materialized_patch_set, ok)
  full_scale = max(1.0_dp, maxval(abs(scheduled_reference_state)))
  call assert_all(ok .and. &
    maxval(abs(materialized_coarse_state - scheduled_reference_state)) <= &
      2.0e-10_dp * full_scale .and. &
    maxval(abs(materialized_coarse_temperature - &
      scheduled_reference_temperature)) <= 2.0e-10_dp * &
        max(1.0_dp, maxval(abs(scheduled_reference_temperature))), &
    "MPI EB AMR scheduled-regrid serial root parity", rank)
  do child = 1, materialized_patch_set%patch_count()
    call assert_all( &
      maxval(abs(materialized_patch_set%children(child)%state - &
        scheduled_reference_set%children(child)%state)) <= &
          2.0e-10_dp * full_scale .and. &
      maxval(abs(materialized_patch_set%children(child)%temperature - &
        scheduled_reference_set%children(child)%temperature)) <= &
          2.0e-10_dp * max(1.0_dp, maxval(abs( &
            scheduled_reference_set%children(child)%temperature))), &
      "MPI EB AMR scheduled-regrid serial child parity", rank)
  end do

  scheduled_failed_distribution = distribution
  scheduled_failed_template = scheduled_start_set
  call scatter_owned_reactive_eb_patch_set_2d( &
    scheduled_failed_distribution, size(species), scheduled_start_state, &
    scheduled_start_temperature, coarse_geometry, scheduled_failed_template, &
    sparse_scheduled_failed_set, ok)
  call assert_all(ok, "MPI EB AMR scheduled-regrid rollback scatter", rank)
  scheduled_time = 0.0_dp
  scheduled_steps = 0
  scheduled_regrid_evaluations = 0
  scheduled_regrids = 0
  failing_geometry_calls = 0
  call advance_sparse_owned_reactive_eb_patch_set_to_time_2d( &
    species, reactions, transport, scheduled_failed_distribution, &
    sparse_scheduled_failed_set, coarse_geometry, scheduled_failed_template, &
    "hllc", "pcm", "mc", 2, scheduled_time, &
    0.25_dp * scheduled_initial_dt, scheduled_steps, 1, 0.35_dp, 0.20_dp, &
    1.0e-8_dp, 1.0e-14_dp, .true., .true., .true., .true., boundaries, ok, &
    scheduled_minimum_dt, advanced_steps=scheduled_advanced_steps, &
    local_timestep_root_transfers=scheduled_timestep_transfers, &
    regrid_evaluations=scheduled_regrid_evaluations, &
    regrids=scheduled_regrids, regrid_interval=1, &
    regrid_criteria=scheduled_criteria, refinement_ratio=ratio, &
    geometry_builder=reject_scheduled_patch_geometry, &
    local_regrid_root_transfers=scheduled_regrid_transfers, &
    local_regrid_restriction_transfers= &
      local_regrid_restriction_transfers, &
    local_regrid_prolongation_transfers= &
      local_regrid_prolongation_transfers, &
    local_regrid_overlap_transfers=local_regrid_overlap_transfers)
  call assert_all(.not. ok .and. failing_geometry_calls >= 1 .and. &
    scheduled_time == 0.0_dp .and. scheduled_steps == 0 .and. &
    scheduled_advanced_steps == 0 .and. scheduled_minimum_dt == 0.0_dp .and. &
    scheduled_timestep_transfers == 0 .and. &
    scheduled_regrid_evaluations == 0 .and. scheduled_regrids == 0 .and. &
    scheduled_regrid_transfers == 0 .and. &
    local_regrid_restriction_transfers == 0 .and. &
    local_regrid_prolongation_transfers == 0 .and. &
    local_regrid_overlap_transfers == 0 .and. &
    scheduled_failed_distribution%child_count() == &
      distribution%child_count() .and. &
    all(scheduled_failed_distribution%child_owners == &
      distribution%child_owners) .and. &
    same_patch_topology(scheduled_failed_template, scheduled_start_set), &
    "MPI EB AMR scheduled-regrid atomic failure", rank)
  call materialize_owned_reactive_eb_patch_set_2d( &
    scheduled_failed_distribution, sparse_scheduled_failed_set, &
    scheduled_start_state, scheduled_start_temperature, coarse_geometry, &
    scheduled_failed_template, materialized_coarse_state, &
    materialized_coarse_temperature, materialized_patch_set, ok)
  call assert_all(ok .and. &
    all(materialized_coarse_state == scheduled_start_state) .and. &
    all(materialized_coarse_temperature == scheduled_start_temperature), &
    "MPI EB AMR scheduled-regrid failed root rollback", rank)
  do child = 1, materialized_patch_set%patch_count()
    call assert_all( &
      all(materialized_patch_set%children(child)%state == &
        scheduled_start_set%children(child)%state) .and. &
      all(materialized_patch_set%children(child)%temperature == &
        scheduled_start_set%children(child)%temperature), &
      "MPI EB AMR scheduled-regrid failed child rollback", rank)
  end do

  call scatter_owned_reactive_eb_patch_set_2d( &
    distribution, size(species), coarse_state, coarse_temperature, &
    coarse_geometry, patch_set, sparse_limited_time_loop_set, ok)
  call assert_all(ok, "MPI EB AMR sparse limited time-loop scatter", rank)
  time_loop_time = 0.0_dp
  time_loop_steps = 0
  call advance_sparse_owned_reactive_eb_patch_set_to_time_2d( &
    species, reactions, transport, distribution, &
    sparse_limited_time_loop_set, coarse_geometry, patch_set, "hllc", &
    "pcm", "mc", 2, time_loop_time, time_loop_final_time, time_loop_steps, &
    1, 0.35_dp, 0.20_dp, 1.0e-8_dp, 1.0e-14_dp, .true., .true., .true., &
    .true., boundaries, ok, time_loop_minimum_dt, &
    time_loop_advanced_steps, local_chemistry_advances, &
    local_hydro_advances, local_transport_advances, full_theta, 0.5_dp, &
    local_root_transfers)
  call assert_all(.not. ok .and. time_loop_steps == 1 .and. &
    time_loop_advanced_steps == 1 .and. time_loop_time > 0.0_dp .and. &
    time_loop_time < time_loop_final_time .and. &
    abs(time_loop_time - time_loop_initial_dt) <= &
      128.0_dp * epsilon(1.0_dp) * time_loop_initial_dt .and. &
    time_loop_minimum_dt == time_loop_time .and. &
    local_chemistry_advances == expected_local_chemistry .and. &
    local_hydro_advances == expected_local_hydro .and. &
    local_transport_advances == expected_local_transport .and. &
    local_root_transfers == 0 .and. &
    sparse_limited_time_loop_set%is_valid( &
      distribution, coarse_geometry, patch_set), &
    "MPI EB AMR sparse time-loop committed step limit", rank)

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

  subroutine remove_nonempty_file(path, valid)
    character(len=*), intent(in) :: path
    logical, intent(inout) :: valid

    logical :: exists
    integer :: file_size, status, unit

    exists = .false.
    file_size = 0
    inquire(file=trim(path), exist=exists, size=file_size, iostat=status)
    if (status /= 0 .or. .not. exists .or. file_size < 1) then
      valid = .false.
      return
    end if
    open(newunit=unit, file=trim(path), status="old", action="read", &
      iostat=status)
    if (status /= 0) then
      valid = .false.
      return
    end if
    close(unit, status="delete", iostat=status)
    valid = valid .and. status == 0
  end subroutine remove_nonempty_file

  pure real(dp) function root_tile_factor(local_tile) result(value)
    integer, intent(in) :: local_tile

    value = 1.0_dp + 1.0e-3_dp * real(local_tile, dp)
  end function root_tile_factor

  pure real(dp) function child_factor(local_child) result(value)
    integer, intent(in) :: local_child

    value = 1.0_dp + 1.0e-2_dp * real(local_child, dp)
  end function child_factor

  subroutine compute_serial_full_physics_timestep( &
      state, state_temperature, state_set, dt, valid)
    real(dp), intent(in) :: state(:, :, :), state_temperature(:, :)
    type(reactive_eb_patch_set_2d), intent(in) :: state_set
    real(dp), intent(out) :: dt
    logical, intent(out) :: valid

    real(dp) :: child_dt, maximum_diffusivity
    integer :: local_child, local_ratio

    dt = 0.0_dp
    valid = .false.
    call compute_reactive_eb_patch_set_cfl_timestep_2d( &
      species, state, state_temperature, coarse_geometry, state_set, &
      0.35_dp, dt, valid)
    if (.not. valid) return
    call reactive_eb_transport_timestep_2d( &
      species, transport, state, state_temperature, coarse_geometry, &
      0.20_dp, .true., .true., .true., child_dt, maximum_diffusivity, valid)
    if (.not. valid) return
    dt = min(dt, child_dt)
    do local_child = 1, state_set%patch_count()
      call reactive_eb_transport_timestep_2d( &
        species, transport, state_set%children(local_child)%state, &
        state_set%children(local_child)%temperature, &
        state_set%children(local_child)%geometry, 0.20_dp, .true., .true., &
        .true., child_dt, maximum_diffusivity, valid)
      if (.not. valid) return
      local_ratio = &
        state_set%children(local_child)%patch%refinement_ratio
      dt = min(dt, real(local_ratio, dp) * child_dt)
    end do
    valid = ieee_is_finite(dt) .and. dt > 0.0_dp
  end subroutine compute_serial_full_physics_timestep

  subroutine advance_serial_full_physics_to_time( &
      start_state, start_temperature, start_set, target_time, final_state, &
      final_temperature, final_set, completed_steps, minimum_dt, &
      minimum_theta, valid)
    real(dp), intent(in) :: start_state(:, :, :), start_temperature(:, :)
    type(reactive_eb_patch_set_2d), intent(in) :: start_set
    real(dp), intent(in) :: target_time
    real(dp), intent(out) :: final_state(:, :, :), final_temperature(:, :)
    type(reactive_eb_patch_set_2d), intent(out) :: final_set
    integer, intent(out) :: completed_steps
    real(dp), intent(out) :: minimum_dt, minimum_theta
    logical, intent(out) :: valid

    type(reactive_eb_patch_set_2d) :: candidate_set
    real(dp), allocatable :: candidate_state(:, :, :)
    real(dp), allocatable :: candidate_temperature(:, :)
    real(dp) :: dt, local_theta, remaining, time, tolerance
    character(len=160) :: failure_context

    final_state = start_state
    final_temperature = start_temperature
    final_set = start_set
    completed_steps = 0
    minimum_dt = 0.0_dp
    minimum_theta = 1.0_dp
    valid = .false.
    time = 0.0_dp
    tolerance = 16.0_dp * epsilon(1.0_dp) * &
      max(tiny(1.0_dp), abs(target_time))
    allocate(candidate_state, mold=start_state)
    allocate(candidate_temperature, mold=start_temperature)
    do
      remaining = target_time - time
      if (remaining <= tolerance) exit
      if (completed_steps >= 16) return
      call compute_serial_full_physics_timestep( &
        final_state, final_temperature, final_set, dt, valid)
      if (.not. valid) return
      dt = min(dt, remaining)
      call advance_reactive_eb_patch_set_strang_2d( &
        species, reactions, final_state, final_temperature, coarse_geometry, &
        final_set, "hllc", "pcm", "mc", 2, dt, .true., 1.0e-8_dp, &
        1.0e-14_dp, candidate_state, candidate_temperature, candidate_set, &
        valid, target_volume_fraction=0.5_dp, &
        failure_context=failure_context, transport=transport, &
        transport_enabled=.true., viscosity_enabled=.true., &
        thermal_conduction_enabled=.true., species_diffusion_enabled=.true., &
        barodiffusion_enabled=.true., minimum_transport_theta=local_theta, &
        boundaries=boundaries)
      if (.not. valid) return
      final_state = candidate_state
      final_temperature = candidate_temperature
      final_set = candidate_set
      time = time + dt
      completed_steps = completed_steps + 1
      if (minimum_dt == 0.0_dp) then
        minimum_dt = dt
      else
        minimum_dt = min(minimum_dt, dt)
      end if
      minimum_theta = min(minimum_theta, local_theta)
    end do
    valid = completed_steps >= 1 .and. &
      abs(time - target_time) <= tolerance .and. &
      ieee_is_finite(minimum_dt) .and. minimum_dt > 0.0_dp
  end subroutine advance_serial_full_physics_to_time

  subroutine initialize_scheduled_regrid_state( &
      state, state_temperature, valid)
    real(dp), intent(out) :: state(:, :, :), state_temperature(:, :)
    logical, intent(out) :: valid

    real(dp), allocatable :: cell_primitive(:), cell_state(:)
    real(dp) :: cell_sound_speed, cell_temperature, pressure
    integer :: local_i, local_j, local_species

    valid = .false.
    if (any(shape(state) /= [nvar, coarse_nx, coarse_ny]) .or. &
        any(shape(state_temperature) /= [coarse_nx, coarse_ny])) return
    allocate(cell_primitive(reactive_nprim(size(species))))
    allocate(cell_state(nvar))
    do local_j = 1, coarse_ny
      do local_i = 1, coarse_nx
        pressure = 135000.0_dp
        if (local_i >= 8 .and. local_i <= 11 .and. &
            local_j >= 8 .and. local_j <= 11) pressure = 1.04_dp * pressure
        cell_primitive(1:5) = [0.31_dp, 0.0_dp, 0.0_dp, 0.0_dp, pressure]
        do local_species = 1, size(species)
          cell_primitive(reactive_mass_fraction_component(local_species)) = &
            mass_fractions(local_species)
        end do
        call reactive_primitive_to_conserved( &
          species, cell_primitive, cell_state, cell_temperature, &
          cell_sound_speed, valid)
        if (.not. valid) return
        state(:, local_i, local_j) = cell_state
        state_temperature(local_i, local_j) = cell_temperature
      end do
    end do
    valid = all(ieee_is_finite(state)) .and. &
      all(ieee_is_finite(state_temperature)) .and. &
      all(state_temperature > 0.0_dp)
  end subroutine initialize_scheduled_regrid_state

  subroutine advance_serial_full_physics_to_time_with_regrid( &
      start_state, start_temperature, start_set, target_time, &
      regrid_interval, regrid_criteria, final_state, final_temperature, &
      final_set, completed_steps, minimum_dt, minimum_theta, &
      regrid_evaluations, regrids, valid)
    real(dp), intent(in) :: start_state(:, :, :), start_temperature(:, :)
    type(reactive_eb_patch_set_2d), intent(in) :: start_set
    real(dp), intent(in) :: target_time
    integer, intent(in) :: regrid_interval
    type(amr_eb_tagging_criteria_2d), intent(in) :: regrid_criteria
    real(dp), intent(out) :: final_state(:, :, :), final_temperature(:, :)
    type(reactive_eb_patch_set_2d), intent(out) :: final_set
    integer, intent(out) :: completed_steps, regrid_evaluations, regrids
    real(dp), intent(out) :: minimum_dt, minimum_theta
    logical, intent(out) :: valid

    type(amr_eb_regrid_plan_collection_2d) :: local_collection
    type(reactive_eb_patch_set_2d) :: candidate_set, regridded_set
    type(eb_geometry_2d), allocatable :: planned_geometries(:)
    real(dp), allocatable :: candidate_state(:, :, :)
    real(dp), allocatable :: candidate_temperature(:, :)
    real(dp), allocatable :: regridded_state(:, :, :)
    real(dp), allocatable :: regridded_temperature(:, :)
    logical, allocatable :: local_tags(:, :)
    real(dp) :: dt, local_theta, remaining, time, tolerance
    logical :: changed
    integer :: local_child
    character(len=160) :: failure_context

    final_state = start_state
    final_temperature = start_temperature
    final_set = start_set
    completed_steps = 0
    regrid_evaluations = 0
    regrids = 0
    minimum_dt = 0.0_dp
    minimum_theta = 1.0_dp
    valid = .false.
    if (regrid_interval < 1 .or. &
        .not. regrid_criteria%is_valid(coarse_nx, coarse_ny)) return
    time = 0.0_dp
    tolerance = 16.0_dp * epsilon(1.0_dp) * &
      max(tiny(1.0_dp), abs(target_time))
    allocate(candidate_state, mold=start_state)
    allocate(candidate_temperature, mold=start_temperature)
    allocate(regridded_state, mold=start_state)
    allocate(regridded_temperature, mold=start_temperature)
    allocate(local_tags(coarse_nx, coarse_ny))
    do
      remaining = target_time - time
      if (remaining <= tolerance) exit
      if (completed_steps >= 16) return
      call compute_serial_full_physics_timestep( &
        final_state, final_temperature, final_set, dt, valid)
      if (.not. valid) return
      dt = min(dt, remaining)
      call advance_reactive_eb_patch_set_strang_2d( &
        species, reactions, final_state, final_temperature, coarse_geometry, &
        final_set, "hllc", "pcm", "mc", 2, dt, .true., 1.0e-8_dp, &
        1.0e-14_dp, candidate_state, candidate_temperature, candidate_set, &
        valid, target_volume_fraction=0.5_dp, &
        failure_context=failure_context, transport=transport, &
        transport_enabled=.true., viscosity_enabled=.true., &
        thermal_conduction_enabled=.true., species_diffusion_enabled=.true., &
        barodiffusion_enabled=.true., minimum_transport_theta=local_theta, &
        boundaries=boundaries)
      if (.not. valid) return
      if (modulo(completed_steps + 1, regrid_interval) == 0) then
        call plan_reactive_eb_temperature_regrid_collection_2d( &
          candidate_temperature, coarse_geometry, regrid_criteria, &
          local_tags, local_collection, valid)
        if (.not. valid) return
        if (allocated(planned_geometries)) deallocate(planned_geometries)
        allocate(planned_geometries(local_collection%patch_count()))
        do local_child = 1, local_collection%patch_count()
          call build_scheduled_patch_geometry( &
            coarse_geometry, &
            local_collection%plans(local_child)%coarse_i_lower, &
            local_collection%plans(local_child)%coarse_i_upper, &
            local_collection%plans(local_child)%coarse_j_lower, &
            local_collection%plans(local_child)%coarse_j_upper, ratio, &
            planned_geometries(local_child), valid)
          if (.not. valid) return
        end do
        call regrid_reactive_eb_patch_set_2d( &
          species, candidate_state, candidate_temperature, coarse_geometry, &
          candidate_set, planned_geometries, local_collection, ratio, &
          regridded_state, regridded_temperature, regridded_set, valid)
        if (.not. valid) return
        changed = .not. same_patch_topology(candidate_set, regridded_set)
        candidate_state = regridded_state
        candidate_temperature = regridded_temperature
        candidate_set = regridded_set
        regrid_evaluations = regrid_evaluations + 1
        if (changed) regrids = regrids + 1
      end if
      final_state = candidate_state
      final_temperature = candidate_temperature
      final_set = candidate_set
      time = time + dt
      completed_steps = completed_steps + 1
      if (minimum_dt == 0.0_dp) then
        minimum_dt = dt
      else
        minimum_dt = min(minimum_dt, dt)
      end if
      minimum_theta = min(minimum_theta, local_theta)
    end do
    valid = completed_steps >= 1 .and. &
      abs(time - target_time) <= tolerance .and. &
      ieee_is_finite(minimum_dt) .and. minimum_dt > 0.0_dp
  end subroutine advance_serial_full_physics_to_time_with_regrid

  pure logical function same_patch_topology(first, second) result(same)
    type(reactive_eb_patch_set_2d), intent(in) :: first, second

    integer :: local_child

    same = first%patch_count() == second%patch_count()
    if (.not. same) return
    do local_child = 1, first%patch_count()
      same = all([ &
        first%children(local_child)%patch%coarse_i_lower, &
        first%children(local_child)%patch%coarse_i_upper, &
        first%children(local_child)%patch%coarse_j_lower, &
        first%children(local_child)%patch%coarse_j_upper, &
        first%children(local_child)%patch%refinement_ratio] == [ &
        second%children(local_child)%patch%coarse_i_lower, &
        second%children(local_child)%patch%coarse_i_upper, &
        second%children(local_child)%patch%coarse_j_lower, &
        second%children(local_child)%patch%coarse_j_upper, &
        second%children(local_child)%patch%refinement_ratio])
      if (.not. same) return
    end do
  end function same_patch_topology

  subroutine build_scheduled_patch_geometry( &
      root_geometry, i_lower, i_upper, j_lower, j_upper, refinement_ratio, &
      fine_geometry, valid)
    type(eb_geometry_2d), intent(in) :: root_geometry
    integer, intent(in) :: i_lower, i_upper, j_lower, j_upper
    integer, intent(in) :: refinement_ratio
    type(eb_geometry_2d), intent(out) :: fine_geometry
    logical, intent(out) :: valid

    type(amr_eb_patch_2d) :: local_patch

    call build_patch_geometry( &
      root_geometry, i_lower, i_upper, j_lower, j_upper, refinement_ratio, &
      fine_geometry, local_patch, valid)
  end subroutine build_scheduled_patch_geometry

  subroutine reject_scheduled_patch_geometry( &
      root_geometry, i_lower, i_upper, j_lower, j_upper, refinement_ratio, &
      fine_geometry, valid)
    type(eb_geometry_2d), intent(in) :: root_geometry
    integer, intent(in) :: i_lower, i_upper, j_lower, j_upper
    integer, intent(in) :: refinement_ratio
    type(eb_geometry_2d), intent(out) :: fine_geometry
    logical, intent(out) :: valid

    failing_geometry_calls = failing_geometry_calls + 1
    fine_geometry = eb_geometry_2d()
    valid = root_geometry%is_valid() .and. i_lower <= i_upper .and. &
      j_lower <= j_upper .and. refinement_ratio >= 2 .and. .false.
  end subroutine reject_scheduled_patch_geometry

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
