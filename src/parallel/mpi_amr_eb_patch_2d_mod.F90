module mpi_amr_eb_patch_2d_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mpi_f08
  use precision_mod, only: dp
  use state_indices_mod, only: irho
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_conserved_to_primitive
  use reactive_2d_mod, only: advance_reactive_chemistry_2d
  use reactive_boundary_2d_mod, only: &
    reactive_boundary_set_2d, validate_reactive_boundary_set_2d, &
    boundary_y_lower, reactive_boundary_is_periodic
  use transport_database_mod, only: &
    gas_transport_species, compatible_transport_database
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_covered_cell
  use eb_reactive_reconstruction_2d_mod, only: &
    reactive_eb_exterior_state_2d
  use amr_eb_reactive_2d_mod, only: &
    reactive_eb_patch_exterior_context_2d, &
    prolong_reactive_eb_patch_pcm_2d, &
    extract_reactive_eb_patch_exterior_context_support_2d, &
    build_reactive_eb_patch_exterior_from_context_2d, &
    build_reactive_eb_patch_exterior_2d, &
    advance_reactive_eb_level_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, build_amr_eb_patch_2d
  use eb_reactive_redistribution_2d_mod, only: &
    advance_reactive_eb_state_redistributed_2d
  use eb_reactive_transport_2d_mod, only: &
    reactive_eb_transport_fluxes_rhs_2d, reactive_eb_transport_timestep_2d
  use amr_eb_flux_register_2d_mod, only: &
    amr_eb_flux_register_2d, initialize_amr_eb_flux_register_2d, &
    accumulate_coarse_eb_fluxes_patch_support_2d, &
    accumulate_coarse_eb_fluxes_2d, accumulate_fine_eb_fluxes_2d, &
    reflux_reactive_eb_state_patch_support_2d, &
    reflux_reactive_eb_state_patch_2d
  use amr_eb_regrid_2d_mod, only: &
    amr_eb_tagging_criteria_2d, amr_eb_regrid_plan_collection_2d, &
    reactive_eb_patch_set_2d, reactive_eb_patch_topology_2d, &
    extract_reactive_eb_patch_topology_2d, &
    average_down_reactive_eb_patch_set_2d, &
    composite_reactive_eb_patch_set_integral_2d, &
    plan_reactive_eb_temperature_regrid_collection_2d
  use amr_eb_transport_2d_mod, only: recover_transport_temperature_2d
  use amr_eb_multilevel_2d_mod, only: &
    mark_local_coarse_fine_interface_recipients_2d
  use amr_eb_multilevel_reactive_2d_mod, only: &
    level_two_interface_is_regular
  use amr_eb_multipatch_transport_2d_mod, only: &
    close_cut_patch_set_conservation_2d
  use reactive_eb_2d_driver_mod, only: &
    compute_reactive_eb_cfl_timestep_2d
  implicit none
  private

  integer, parameter :: sparse_restriction_tag = 2701
  integer, parameter :: sparse_root_gather_tag = 2702
  integer, parameter :: sparse_regrid_prolongation_tag = 2707
  integer, parameter :: sparse_regrid_overlap_tag = 2708
  integer, parameter :: sparse_root_materialization_tag = 2709
  integer, parameter :: sparse_root_restart_scatter_tag = 2710
  integer, parameter :: sparse_root_halo_tag = 2711
  integer, parameter :: sparse_child_coarse_flux_support_tag = 2715
  integer, parameter :: sparse_child_state_support_tag = 2716
  integer, parameter :: sparse_child_state_correction_tag = 2717
  integer, parameter, public :: mpi_amr_eb_root_tile_hydro_halo_cells = 6
  integer, parameter, public :: mpi_amr_eb_root_tile_transport_halo_cells = 6

  type, public :: mpi_amr_eb_root_tile_2d
    integer :: owner = -1
    integer :: i_lower = 1
    integer :: i_upper = 0
    integer :: j_lower = 1
    integer :: j_upper = 0
    integer :: cell_count = 0
    integer(int64) :: work_count = 0_int64
  contains
    procedure :: is_valid => mpi_amr_eb_root_tile_is_valid
  end type mpi_amr_eb_root_tile_2d

  type :: mpi_amr_eb_root_tile_transport_flux_2d
    real(dp), allocatable :: x_flux(:, :, :)
    real(dp), allocatable :: y_flux(:, :, :)
  end type mpi_amr_eb_root_tile_transport_flux_2d

  type :: mpi_amr_eb_root_tile_transport_state_2d
    real(dp), allocatable :: start_state(:, :, :)
    real(dp), allocatable :: start_temperature(:, :)
    real(dp), allocatable :: end_state(:, :, :)
    real(dp), allocatable :: end_temperature(:, :)
    real(dp), allocatable :: corrected_state(:, :, :)
    real(dp), allocatable :: corrected_temperature(:, :)
  end type mpi_amr_eb_root_tile_transport_state_2d

  type, public :: mpi_amr_eb_patch_distribution_2d
    type(MPI_Comm) :: comm = MPI_COMM_NULL
    integer :: rank = -1
    integer :: nranks = 0
    integer :: subcycle_exponent = 0
    type(mpi_amr_eb_root_tile_2d), allocatable :: root_tiles(:)
    integer, allocatable :: child_owners(:)
    integer, allocatable :: child_cell_counts(:)
    integer(int64), allocatable :: child_work_counts(:)
    integer, allocatable :: rank_cell_counts(:)
    integer, allocatable :: rank_entity_counts(:)
    integer(int64), allocatable :: rank_work_counts(:)
  contains
    procedure :: root_tile_count => mpi_amr_eb_root_tile_count
    procedure :: child_count => mpi_amr_eb_child_count
    procedure :: child_owner => mpi_amr_eb_child_owner
    procedure :: root_level_owner => mpi_amr_eb_root_level_owner
    procedure :: root_tile_is_local => mpi_amr_eb_root_tile_is_local
    procedure :: child_is_local => mpi_amr_eb_child_is_local
    procedure :: is_valid => mpi_amr_eb_distribution_is_valid
  end type mpi_amr_eb_patch_distribution_2d

  type, public :: mpi_amr_eb_sparse_field_2d
    real(dp), allocatable :: state(:, :, :)
    real(dp), allocatable :: temperature(:, :)
  end type mpi_amr_eb_sparse_field_2d

  type, public :: mpi_amr_eb_sparse_patch_set_2d
    type(mpi_amr_eb_sparse_field_2d), allocatable :: root_tiles(:)
    type(mpi_amr_eb_sparse_field_2d), allocatable :: children(:)
    integer :: rank = -1
    integer :: nranks = 0
    integer :: nvar = 0
  contains
    procedure :: is_valid => mpi_amr_eb_sparse_patch_set_is_valid
    procedure :: local_value_count => &
      mpi_amr_eb_sparse_patch_set_local_value_count
  end type mpi_amr_eb_sparse_patch_set_2d

  abstract interface
    subroutine sparse_eb_geometry_builder_2d( &
        coarse_geometry, coarse_i_lower, coarse_i_upper, coarse_j_lower, &
        coarse_j_upper, refinement_ratio, fine_geometry, ok)
      import :: eb_geometry_2d
      type(eb_geometry_2d), intent(in) :: coarse_geometry
      integer, intent(in) :: coarse_i_lower, coarse_i_upper
      integer, intent(in) :: coarse_j_lower, coarse_j_upper
      integer, intent(in) :: refinement_ratio
      type(eb_geometry_2d), intent(out) :: fine_geometry
      logical, intent(out) :: ok
    end subroutine sparse_eb_geometry_builder_2d
  end interface

  public :: initialize_mpi_amr_eb_patch_distribution_2d
  public :: mpi_amr_eb_child_transport_context_value_count_2d
  public :: mpi_amr_eb_child_transport_reflux_context_value_count_2d
  public :: mpi_amr_eb_child_transport_state_context_value_count_2d
  public :: mpi_amr_eb_child_transport_tile_state_support_value_count_2d
  public :: mpi_amr_eb_child_coarse_flux_support_value_count_2d
  public :: mpi_amr_eb_distribution_matches_patch_set_2d
  public :: synchronize_owned_reactive_eb_patch_set_2d
  public :: advance_owned_reactive_eb_patch_set_chemistry_2d
  public :: advance_owned_reactive_eb_patch_set_hydro_2d
  public :: advance_owned_reactive_eb_patch_set_transport_2d
  public :: advance_owned_reactive_eb_patch_set_strang_2d
  public :: scatter_owned_reactive_eb_patch_set_2d
  public :: materialize_owned_reactive_eb_patch_set_2d
  public :: gather_sparse_owned_reactive_eb_patch_set_to_root_2d
  public :: scatter_root_reactive_eb_patch_set_to_sparse_2d
  public :: scatter_root_reactive_eb_topology_to_sparse_2d
  public :: mpi_amr_eb_distribution_matches_topology_2d
  public :: mpi_amr_eb_sparse_patch_set_matches_topology_2d
  public :: regrid_sparse_owned_reactive_eb_patch_set_2d
  public :: regrid_tagged_sparse_owned_reactive_eb_patch_set_2d
  public :: average_down_sparse_owned_reactive_eb_patch_set_2d
  public :: advance_sparse_owned_reactive_eb_patch_set_chemistry_2d
  public :: compute_sparse_owned_reactive_eb_patch_set_timestep_2d
  public :: advance_sparse_owned_reactive_eb_patch_set_hydro_2d
  public :: advance_sparse_owned_reactive_eb_patch_set_transport_2d
  public :: advance_sparse_owned_reactive_eb_patch_set_strang_2d
  public :: advance_sparse_owned_reactive_eb_patch_set_to_time_2d

contains

  pure integer(int64) function &
      mpi_amr_eb_child_transport_context_value_count_2d( &
        component_count, coarse_geometry, fine_geometry, patch) &
      result(value_count)
    integer, intent(in) :: component_count
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch

    integer(int64) :: correction_cells, edge_values

    value_count = -1_int64
    if (component_count < 1 .or. &
        .not. patch%is_valid(coarse_geometry, fine_geometry)) return
    edge_values = 4_int64 * int(component_count + 1, int64) * &
      int(fine_geometry%nx + fine_geometry%ny, int64)
    correction_cells = &
      int(min(coarse_geometry%nx, patch%coarse_i_upper + 1) - &
        max(1, patch%coarse_i_lower - 1) + 1, int64) * &
      int(min(coarse_geometry%ny, patch%coarse_j_upper + 1) - &
        max(1, patch%coarse_j_lower - 1) + 1, int64)
    value_count = edge_values + &
      int(component_count, int64) * correction_cells
  end function mpi_amr_eb_child_transport_context_value_count_2d

  pure integer(int64) function &
      mpi_amr_eb_child_transport_reflux_context_value_count_2d( &
        component_count, coarse_geometry, fine_geometry, patch) &
      result(value_count)
    integer, intent(in) :: component_count
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch

    integer(int64) :: boundary_context_values, support_cells

    value_count = -1_int64
    boundary_context_values = &
      mpi_amr_eb_child_transport_context_value_count_2d( &
        component_count, coarse_geometry, fine_geometry, patch)
    if (boundary_context_values < 1_int64) return
    support_cells = &
      int(min(coarse_geometry%nx, patch%coarse_i_upper + 2) - &
        max(1, patch%coarse_i_lower - 2) + 1, int64) * &
      int(min(coarse_geometry%ny, patch%coarse_j_upper + 2) - &
        max(1, patch%coarse_j_lower - 2) + 1, int64)
    value_count = boundary_context_values + &
      int(component_count + 1, int64) * support_cells
  end function mpi_amr_eb_child_transport_reflux_context_value_count_2d

  pure integer(int64) function &
      mpi_amr_eb_child_transport_state_context_value_count_2d( &
        component_count, coarse_geometry, fine_geometry, patch) &
      result(value_count)
    integer, intent(in) :: component_count
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch

    integer(int64) :: edge_values, support_cells

    value_count = -1_int64
    if (component_count < 1 .or. &
        .not. patch%is_valid(coarse_geometry, fine_geometry)) return
    edge_values = 4_int64 * int(component_count + 1, int64) * &
      int(fine_geometry%nx + fine_geometry%ny, int64)
    support_cells = &
      int(min(coarse_geometry%nx, patch%coarse_i_upper + 2) - &
        max(1, patch%coarse_i_lower - 2) + 1, int64) * &
      int(min(coarse_geometry%ny, patch%coarse_j_upper + 2) - &
        max(1, patch%coarse_j_lower - 2) + 1, int64)
    value_count = edge_values + &
      int(component_count + 1, int64) * support_cells
  end function mpi_amr_eb_child_transport_state_context_value_count_2d

  pure integer(int64) function &
      mpi_amr_eb_child_transport_tile_state_support_value_count_2d( &
        component_count, coarse_geometry, fine_geometry, patch) &
      result(value_count)
    integer, intent(in) :: component_count
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch

    integer(int64) :: support_cells

    value_count = -1_int64
    if (component_count < 1 .or. &
        .not. patch%is_valid(coarse_geometry, fine_geometry)) return
    support_cells = &
      int(min(coarse_geometry%nx, patch%coarse_i_upper + 2) - &
        max(1, patch%coarse_i_lower - 2) + 1, int64) * &
      int(min(coarse_geometry%ny, patch%coarse_j_upper + 2) - &
        max(1, patch%coarse_j_lower - 2) + 1, int64)
    value_count = 3_int64 * int(component_count + 1, int64) * support_cells
  end function mpi_amr_eb_child_transport_tile_state_support_value_count_2d

  pure integer(int64) function &
      mpi_amr_eb_child_coarse_flux_support_value_count_2d( &
        component_count, coarse_geometry, fine_geometry, patch) &
      result(value_count)
    integer, intent(in) :: component_count
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch

    integer(int64) :: patch_height, patch_width

    value_count = -1_int64
    if (component_count < 1 .or. &
        .not. patch%is_valid(coarse_geometry, fine_geometry)) return
    patch_width = int( &
      patch%coarse_i_upper - patch%coarse_i_lower + 1, int64)
    patch_height = int( &
      patch%coarse_j_upper - patch%coarse_j_lower + 1, int64)
    value_count = int(component_count, int64) * ( &
      (patch_width + 1_int64) * patch_height + &
      patch_width * (patch_height + 1_int64))
  end function mpi_amr_eb_child_coarse_flux_support_value_count_2d

  pure logical function mpi_amr_eb_root_tile_is_valid( &
      self, nx, ny, nranks) result(valid)
    class(mpi_amr_eb_root_tile_2d), intent(in) :: self
    integer, intent(in) :: nx, ny, nranks
    integer(int64) :: expected_cells

    valid = self%owner >= 0 .and. self%owner < nranks .and. &
      self%i_lower == 1 .and. self%i_upper == nx .and. &
      self%j_lower >= 1 .and. self%j_upper <= ny .and. &
      self%j_upper >= self%j_lower
    if (.not. valid) return
    expected_cells = int(nx, int64) * &
      int(self%j_upper - self%j_lower + 1, int64)
    valid = expected_cells <= int(huge(self%cell_count), int64) .and. &
      int(self%cell_count, int64) == expected_cells .and. &
      self%work_count == expected_cells
  end function mpi_amr_eb_root_tile_is_valid

  pure integer function mpi_amr_eb_root_tile_count(self) result(count)
    class(mpi_amr_eb_patch_distribution_2d), intent(in) :: self

    count = 0
    if (allocated(self%root_tiles)) count = size(self%root_tiles)
  end function mpi_amr_eb_root_tile_count

  pure integer function mpi_amr_eb_child_count(self) result(count)
    class(mpi_amr_eb_patch_distribution_2d), intent(in) :: self

    count = 0
    if (allocated(self%child_owners)) count = size(self%child_owners)
  end function mpi_amr_eb_child_count

  pure integer function mpi_amr_eb_child_owner(self, child) result(owner)
    class(mpi_amr_eb_patch_distribution_2d), intent(in) :: self
    integer, intent(in) :: child

    owner = -1
    if (.not. allocated(self%child_owners)) return
    if (child < 1 .or. child > size(self%child_owners)) return
    owner = self%child_owners(child)
  end function mpi_amr_eb_child_owner

  pure integer function mpi_amr_eb_root_level_owner(self) result(owner)
    class(mpi_amr_eb_patch_distribution_2d), intent(in) :: self

    owner = -1
    if (.not. allocated(self%root_tiles)) return
    if (size(self%root_tiles) < 1) return
    owner = self%root_tiles(1)%owner
  end function mpi_amr_eb_root_level_owner

  pure logical function mpi_amr_eb_root_tile_is_local( &
      self, tile) result(local)
    class(mpi_amr_eb_patch_distribution_2d), intent(in) :: self
    integer, intent(in) :: tile

    local = .false.
    if (.not. allocated(self%root_tiles)) return
    if (tile < 1 .or. tile > size(self%root_tiles)) return
    local = self%root_tiles(tile)%owner == self%rank
  end function mpi_amr_eb_root_tile_is_local

  pure logical function mpi_amr_eb_child_is_local( &
      self, child) result(local)
    class(mpi_amr_eb_patch_distribution_2d), intent(in) :: self
    integer, intent(in) :: child

    local = self%child_owner(child) == self%rank
  end function mpi_amr_eb_child_is_local

  pure logical function mpi_amr_eb_distribution_is_valid( &
      self, coarse_geometry, patch_set) result(valid)
    class(mpi_amr_eb_patch_distribution_2d), intent(in) :: self
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set

    integer, allocatable :: cells(:), entities(:)
    integer(int64), allocatable :: work(:)
    integer(int64) :: expected_cells, expected_work, level_scale
    integer :: child, exponent, owner, ratio, tile
    integer :: expected_j_lower

    valid = self%rank >= 0 .and. self%nranks >= 1 .and. &
      self%rank < self%nranks .and. &
      self%subcycle_exponent >= 0 .and. self%subcycle_exponent <= 2 .and. &
      coarse_geometry%is_valid() .and. allocated(patch_set%children) .and. &
      allocated(self%root_tiles) .and. allocated(self%child_owners) .and. &
      allocated(self%child_cell_counts) .and. &
      allocated(self%child_work_counts) .and. &
      allocated(self%rank_cell_counts) .and. &
      allocated(self%rank_entity_counts) .and. &
      allocated(self%rank_work_counts)
    if (.not. valid) return
    valid = size(self%root_tiles) >= 1 .and. &
      size(self%root_tiles) <= min(coarse_geometry%ny, self%nranks) .and. &
      size(self%child_owners) == patch_set%patch_count() .and. &
      size(self%child_cell_counts) == patch_set%patch_count() .and. &
      size(self%child_work_counts) == patch_set%patch_count() .and. &
      size(self%rank_cell_counts) == self%nranks .and. &
      size(self%rank_entity_counts) == self%nranks .and. &
      size(self%rank_work_counts) == self%nranks
    if (.not. valid) return

    allocate(cells(self%nranks), entities(self%nranks), work(self%nranks))
    cells = 0
    entities = 0
    work = 0_int64
    expected_j_lower = 1
    do tile = 1, size(self%root_tiles)
      valid = self%root_tiles(tile)%is_valid( &
        coarse_geometry%nx, coarse_geometry%ny, self%nranks) .and. &
        self%root_tiles(tile)%j_lower == expected_j_lower
      if (.not. valid) return
      expected_j_lower = self%root_tiles(tile)%j_upper + 1
      owner = self%root_tiles(tile)%owner + 1
      cells(owner) = cells(owner) + self%root_tiles(tile)%cell_count
      entities(owner) = entities(owner) + 1
      work(owner) = work(owner) + self%root_tiles(tile)%work_count
    end do
    if (expected_j_lower /= coarse_geometry%ny + 1) then
      valid = .false.
      return
    end if

    do child = 1, patch_set%patch_count()
      valid = patch_set%children(child)%geometry%is_valid() .and. &
        patch_set%children(child)%patch%is_valid( &
          coarse_geometry, patch_set%children(child)%geometry) .and. &
        self%child_owners(child) >= 0 .and. &
        self%child_owners(child) < self%nranks
      if (.not. valid) return
      expected_cells = &
        int(patch_set%children(child)%geometry%nx, int64) * &
        int(patch_set%children(child)%geometry%ny, int64)
      if (expected_cells > int(huge(self%child_cell_counts(child)), int64)) then
        valid = .false.
        return
      end if
      ratio = patch_set%children(child)%patch%refinement_ratio
      level_scale = 1_int64
      do exponent = 1, self%subcycle_exponent
        if (ratio < 1 .or. &
            level_scale > huge(level_scale) / int(ratio, int64)) then
          valid = .false.
          return
        end if
        level_scale = level_scale * int(ratio, int64)
      end do
      if (expected_cells > huge(expected_work) / level_scale) then
        valid = .false.
        return
      end if
      expected_work = expected_cells * level_scale
      valid = int(self%child_cell_counts(child), int64) == expected_cells .and. &
        self%child_work_counts(child) == expected_work
      if (.not. valid) return
      owner = self%child_owners(child) + 1
      cells(owner) = cells(owner) + self%child_cell_counts(child)
      entities(owner) = entities(owner) + 1
      work(owner) = work(owner) + self%child_work_counts(child)
    end do
    valid = all(cells == self%rank_cell_counts) .and. &
      all(entities == self%rank_entity_counts) .and. &
      all(work == self%rank_work_counts) .and. &
      sum(self%rank_entity_counts) == &
        size(self%root_tiles) + patch_set%patch_count()
  end function mpi_amr_eb_distribution_is_valid

  pure logical function mpi_amr_eb_distribution_matches_topology_2d( &
      distribution, coarse_geometry, topology) result(matches)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_topology_2d), intent(in) :: topology

    integer, allocatable :: cells(:), entities(:)
    integer(int64), allocatable :: work(:)
    integer(int64) :: expected_cells, expected_work, level_scale
    integer :: child, exponent, expected_j_lower, owner, ratio, tile

    matches = distribution%rank >= 0 .and. distribution%nranks >= 1 .and. &
      distribution%rank < distribution%nranks .and. &
      distribution%subcycle_exponent >= 0 .and. &
      distribution%subcycle_exponent <= 2 .and. &
      topology%is_valid(coarse_geometry) .and. &
      allocated(distribution%root_tiles) .and. &
      allocated(distribution%child_owners) .and. &
      allocated(distribution%child_cell_counts) .and. &
      allocated(distribution%child_work_counts) .and. &
      allocated(distribution%rank_cell_counts) .and. &
      allocated(distribution%rank_entity_counts) .and. &
      allocated(distribution%rank_work_counts)
    if (.not. matches) return
    matches = size(distribution%root_tiles) >= 1 .and. &
      size(distribution%root_tiles) <= &
        min(coarse_geometry%ny, distribution%nranks) .and. &
      size(distribution%child_owners) == topology%patch_count() .and. &
      size(distribution%child_cell_counts) == topology%patch_count() .and. &
      size(distribution%child_work_counts) == topology%patch_count() .and. &
      size(distribution%rank_cell_counts) == distribution%nranks .and. &
      size(distribution%rank_entity_counts) == distribution%nranks .and. &
      size(distribution%rank_work_counts) == distribution%nranks
    if (.not. matches) return

    allocate(cells(distribution%nranks), entities(distribution%nranks))
    allocate(work(distribution%nranks))
    cells = 0
    entities = 0
    work = 0_int64
    expected_j_lower = 1
    do tile = 1, size(distribution%root_tiles)
      matches = distribution%root_tiles(tile)%is_valid( &
        coarse_geometry%nx, coarse_geometry%ny, distribution%nranks) .and. &
        distribution%root_tiles(tile)%j_lower == expected_j_lower
      if (.not. matches) return
      expected_j_lower = distribution%root_tiles(tile)%j_upper + 1
      owner = distribution%root_tiles(tile)%owner + 1
      cells(owner) = cells(owner) + &
        distribution%root_tiles(tile)%cell_count
      entities(owner) = entities(owner) + 1
      work(owner) = work(owner) + distribution%root_tiles(tile)%work_count
    end do
    if (expected_j_lower /= coarse_geometry%ny + 1) then
      matches = .false.
      return
    end if

    do child = 1, topology%patch_count()
      matches = distribution%child_owners(child) >= 0 .and. &
        distribution%child_owners(child) < distribution%nranks
      if (.not. matches) return
      expected_cells = int(topology%children(child)%geometry%nx, int64) * &
        int(topology%children(child)%geometry%ny, int64)
      if (expected_cells > &
          int(huge(distribution%child_cell_counts(child)), int64)) then
        matches = .false.
        return
      end if
      ratio = topology%children(child)%patch%refinement_ratio
      level_scale = 1_int64
      do exponent = 1, distribution%subcycle_exponent
        if (ratio < 1 .or. &
            level_scale > huge(level_scale) / int(ratio, int64)) then
          matches = .false.
          return
        end if
        level_scale = level_scale * int(ratio, int64)
      end do
      if (expected_cells > huge(expected_work) / level_scale) then
        matches = .false.
        return
      end if
      expected_work = expected_cells * level_scale
      matches = int(distribution%child_cell_counts(child), int64) == &
          expected_cells .and. &
        distribution%child_work_counts(child) == expected_work
      if (.not. matches) return
      owner = distribution%child_owners(child) + 1
      cells(owner) = cells(owner) + distribution%child_cell_counts(child)
      entities(owner) = entities(owner) + 1
      work(owner) = work(owner) + distribution%child_work_counts(child)
    end do
    matches = all(cells == distribution%rank_cell_counts) .and. &
      all(entities == distribution%rank_entity_counts) .and. &
      all(work == distribution%rank_work_counts) .and. &
      sum(distribution%rank_entity_counts) == &
        size(distribution%root_tiles) + topology%patch_count()
  end function mpi_amr_eb_distribution_matches_topology_2d

  logical function mpi_amr_eb_sparse_patch_set_is_valid( &
      self, distribution, coarse_geometry, patch_set) result(valid)
    class(mpi_amr_eb_sparse_patch_set_2d), intent(in) :: self
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set

    logical :: local
    integer :: child, height, tile

    valid = self%rank == distribution%rank .and. &
      self%nranks == distribution%nranks .and. self%nvar >= 1 .and. &
      allocated(self%root_tiles) .and. allocated(self%children) .and. &
      size(self%root_tiles) == distribution%root_tile_count() .and. &
      size(self%children) == distribution%child_count() .and. &
      distribution%is_valid(coarse_geometry, patch_set) .and. &
      patch_set%is_valid(coarse_geometry, self%nvar)
    if (.not. valid) return
    do tile = 1, size(self%root_tiles)
      local = distribution%root_tile_is_local(tile)
      valid = (allocated(self%root_tiles(tile)%state) .eqv. local) .and. &
        (allocated(self%root_tiles(tile)%temperature) .eqv. local)
      if (.not. valid) return
      if (local) then
        height = distribution%root_tiles(tile)%j_upper - &
          distribution%root_tiles(tile)%j_lower + 1
        valid = all(shape(self%root_tiles(tile)%state) == &
          [self%nvar, coarse_geometry%nx, height]) .and. &
          all(shape(self%root_tiles(tile)%temperature) == &
            [coarse_geometry%nx, height]) .and. &
          all(ieee_is_finite(self%root_tiles(tile)%state)) .and. &
          all(ieee_is_finite(self%root_tiles(tile)%temperature))
        if (.not. valid) return
      end if
    end do
    do child = 1, size(self%children)
      local = distribution%child_is_local(child)
      valid = (allocated(self%children(child)%state) .eqv. local) .and. &
        (allocated(self%children(child)%temperature) .eqv. local)
      if (.not. valid) return
      if (local) then
        valid = all(shape(self%children(child)%state) == shape( &
          patch_set%children(child)%state)) .and. &
          all(shape(self%children(child)%temperature) == shape( &
            patch_set%children(child)%temperature)) .and. &
          all(ieee_is_finite(self%children(child)%state)) .and. &
          all(ieee_is_finite(self%children(child)%temperature))
        if (.not. valid) return
      end if
    end do
  end function mpi_amr_eb_sparse_patch_set_is_valid

  logical function mpi_amr_eb_sparse_patch_set_matches_topology_2d( &
      sparse_patch_set, distribution, coarse_geometry, topology) &
      result(matches)
    type(mpi_amr_eb_sparse_patch_set_2d), intent(in) :: sparse_patch_set
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_topology_2d), intent(in) :: topology

    logical :: local
    integer :: child, height, tile

    matches = sparse_patch_set%rank == distribution%rank .and. &
      sparse_patch_set%nranks == distribution%nranks .and. &
      sparse_patch_set%nvar >= 1 .and. &
      allocated(sparse_patch_set%root_tiles) .and. &
      allocated(sparse_patch_set%children) .and. &
      mpi_amr_eb_distribution_matches_topology_2d( &
        distribution, coarse_geometry, topology)
    if (.not. matches) return
    matches = size(sparse_patch_set%root_tiles) == &
        distribution%root_tile_count() .and. &
      size(sparse_patch_set%children) == distribution%child_count()
    if (.not. matches) return
    do tile = 1, size(sparse_patch_set%root_tiles)
      local = distribution%root_tile_is_local(tile)
      matches = &
        (allocated(sparse_patch_set%root_tiles(tile)%state) .eqv. local) &
        .and. (allocated( &
          sparse_patch_set%root_tiles(tile)%temperature) .eqv. local)
      if (.not. matches) return
      if (local) then
        height = distribution%root_tiles(tile)%j_upper - &
          distribution%root_tiles(tile)%j_lower + 1
        matches = all(shape(sparse_patch_set%root_tiles(tile)%state) == &
            [sparse_patch_set%nvar, coarse_geometry%nx, height]) .and. &
          all(shape(sparse_patch_set%root_tiles(tile)%temperature) == &
            [coarse_geometry%nx, height]) .and. &
          all(ieee_is_finite( &
            sparse_patch_set%root_tiles(tile)%state)) .and. &
          all(ieee_is_finite( &
            sparse_patch_set%root_tiles(tile)%temperature))
        if (.not. matches) return
      end if
    end do
    do child = 1, size(sparse_patch_set%children)
      local = distribution%child_is_local(child)
      matches = (allocated( &
          sparse_patch_set%children(child)%state) .eqv. local) .and. &
        (allocated(sparse_patch_set%children(child)%temperature) .eqv. local)
      if (.not. matches) return
      if (local) then
        matches = all(shape(sparse_patch_set%children(child)%state) == &
            [sparse_patch_set%nvar, &
              topology%children(child)%geometry%nx, &
              topology%children(child)%geometry%ny]) .and. &
          all(shape(sparse_patch_set%children(child)%temperature) == &
            [topology%children(child)%geometry%nx, &
              topology%children(child)%geometry%ny]) .and. &
          all(ieee_is_finite(sparse_patch_set%children(child)%state)) .and. &
          all(ieee_is_finite( &
            sparse_patch_set%children(child)%temperature))
        if (.not. matches) return
      end if
    end do
  end function mpi_amr_eb_sparse_patch_set_matches_topology_2d

  pure integer(int64) function &
      mpi_amr_eb_sparse_patch_set_local_value_count(self) result(count)
    class(mpi_amr_eb_sparse_patch_set_2d), intent(in) :: self
    integer :: child, tile

    count = 0_int64
    if (allocated(self%root_tiles)) then
      do tile = 1, size(self%root_tiles)
        if (allocated(self%root_tiles(tile)%state)) count = count + &
          int(size(self%root_tiles(tile)%state), int64)
        if (allocated(self%root_tiles(tile)%temperature)) count = count + &
          int(size(self%root_tiles(tile)%temperature), int64)
      end do
    end if
    if (allocated(self%children)) then
      do child = 1, size(self%children)
        if (allocated(self%children(child)%state)) count = count + &
          int(size(self%children(child)%state), int64)
        if (allocated(self%children(child)%temperature)) count = count + &
          int(size(self%children(child)%temperature), int64)
      end do
    end if
  end function mpi_amr_eb_sparse_patch_set_local_value_count

  subroutine scatter_owned_reactive_eb_patch_set_2d( &
      distribution, nspecies, coarse_state, coarse_temperature, &
      coarse_geometry, patch_set, sparse_patch_set, ok)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    integer, intent(in) :: nspecies
    real(dp), intent(in) :: coarse_state(:, :, :)
    real(dp), intent(in) :: coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set
    type(mpi_amr_eb_sparse_patch_set_2d), intent(out) :: sparse_patch_set
    logical, intent(out) :: ok

    type(mpi_amr_eb_sparse_patch_set_2d) :: candidate
    logical :: accepted, global_ok, local_ok
    integer :: child, j_lower, j_upper, nvar, tile

    sparse_patch_set = mpi_amr_eb_sparse_patch_set_2d()
    ok = .false.
    nvar = reactive_nvar(nspecies)
    local_ok = nspecies >= 1 .and. &
      all(shape(coarse_state) == &
        [nvar, coarse_geometry%nx, coarse_geometry%ny]) .and. &
      all(shape(coarse_temperature) == &
        [coarse_geometry%nx, coarse_geometry%ny]) .and. &
      distribution%is_valid(coarse_geometry, patch_set) .and. &
      patch_set%is_valid(coarse_geometry, nvar)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    candidate%rank = distribution%rank
    candidate%nranks = distribution%nranks
    candidate%nvar = nvar
    allocate(candidate%root_tiles(distribution%root_tile_count()))
    allocate(candidate%children(distribution%child_count()))
    do tile = 1, distribution%root_tile_count()
      if (.not. distribution%root_tile_is_local(tile)) cycle
      j_lower = distribution%root_tiles(tile)%j_lower
      j_upper = distribution%root_tiles(tile)%j_upper
      allocate(candidate%root_tiles(tile)%state, &
        source=coarse_state(:, :, j_lower:j_upper))
      allocate(candidate%root_tiles(tile)%temperature, &
        source=coarse_temperature(:, j_lower:j_upper))
    end do
    do child = 1, distribution%child_count()
      if (.not. distribution%child_is_local(child)) cycle
      allocate(candidate%children(child)%state, &
        source=patch_set%children(child)%state)
      allocate(candidate%children(child)%temperature, &
        source=patch_set%children(child)%temperature)
    end do
    local_ok = candidate%is_valid(distribution, coarse_geometry, patch_set)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    sparse_patch_set = candidate
    ok = .true.
  end subroutine scatter_owned_reactive_eb_patch_set_2d

  subroutine materialize_owned_reactive_eb_patch_set_2d( &
      distribution, sparse_patch_set, fallback_coarse_state, &
      fallback_coarse_temperature, coarse_geometry, patch_set_template, &
      coarse_state, coarse_temperature, patch_set, ok)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(in) :: sparse_patch_set
    real(dp), intent(in) :: fallback_coarse_state(:, :, :)
    real(dp), intent(in) :: fallback_coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set_template
    real(dp), intent(out) :: coarse_state(:, :, :)
    real(dp), intent(out) :: coarse_temperature(:, :)
    type(reactive_eb_patch_set_2d), intent(out) :: patch_set
    logical, intent(out) :: ok

    type(reactive_eb_patch_set_2d) :: candidate_set
    real(dp), allocatable :: candidate_state(:, :, :)
    real(dp), allocatable :: candidate_temperature(:, :)
    logical :: accepted, global_ok, local_ok
    integer :: child, ierr, j_lower, j_upper, owner, tile

    coarse_state = fallback_coarse_state
    coarse_temperature = fallback_coarse_temperature
    patch_set = patch_set_template
    ok = .false.
    local_ok = all(shape(coarse_state) == shape(fallback_coarse_state)) .and. &
      all(shape(coarse_temperature) == &
        shape(fallback_coarse_temperature)) .and. &
      sparse_patch_set%is_valid( &
        distribution, coarse_geometry, patch_set_template)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    allocate(candidate_state, source=fallback_coarse_state)
    allocate(candidate_temperature, source=fallback_coarse_temperature)
    candidate_set = patch_set_template
    do tile = 1, distribution%root_tile_count()
      j_lower = distribution%root_tiles(tile)%j_lower
      j_upper = distribution%root_tiles(tile)%j_upper
      owner = distribution%root_tiles(tile)%owner
      if (distribution%root_tile_is_local(tile)) then
        candidate_state(:, :, j_lower:j_upper) = &
          sparse_patch_set%root_tiles(tile)%state
        candidate_temperature(:, j_lower:j_upper) = &
          sparse_patch_set%root_tiles(tile)%temperature
      end if
      call MPI_Bcast( &
        candidate_state(:, :, j_lower:j_upper), &
        size(candidate_state(:, :, j_lower:j_upper)), MPI_DOUBLE_PRECISION, &
        owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
      call MPI_Bcast( &
        candidate_temperature(:, j_lower:j_upper), &
        size(candidate_temperature(:, j_lower:j_upper)), &
        MPI_DOUBLE_PRECISION, owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
    end do
    do child = 1, distribution%child_count()
      owner = distribution%child_owner(child)
      if (distribution%child_is_local(child)) then
        candidate_set%children(child)%state = &
          sparse_patch_set%children(child)%state
        candidate_set%children(child)%temperature = &
          sparse_patch_set%children(child)%temperature
      end if
      call MPI_Bcast( &
        candidate_set%children(child)%state, &
        size(candidate_set%children(child)%state), MPI_DOUBLE_PRECISION, &
        owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
      call MPI_Bcast( &
        candidate_set%children(child)%temperature, &
        size(candidate_set%children(child)%temperature), &
        MPI_DOUBLE_PRECISION, owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
    end do
    local_ok = candidate_set%is_valid( &
      coarse_geometry, sparse_patch_set%nvar) .and. &
      all(ieee_is_finite(candidate_state)) .and. &
      all(ieee_is_finite(candidate_temperature))
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    coarse_state = candidate_state
    coarse_temperature = candidate_temperature
    patch_set = candidate_set
    ok = .true.
  end subroutine materialize_owned_reactive_eb_patch_set_2d

  subroutine gather_sparse_owned_reactive_eb_patch_set_to_root_2d( &
      distribution, sparse_patch_set, coarse_geometry, patch_set_template, &
      root, coarse_state, coarse_temperature, patch_set, ok, &
      local_transfers)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(in) :: sparse_patch_set
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set_template
    integer, intent(in) :: root
    real(dp), allocatable, intent(out) :: coarse_state(:, :, :)
    real(dp), allocatable, intent(out) :: coarse_temperature(:, :)
    type(reactive_eb_patch_set_2d), intent(out) :: patch_set
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_transfers

    type(reactive_eb_patch_set_2d) :: candidate_patch_set
    type(MPI_Status) :: status
    real(dp), allocatable :: candidate_state(:, :, :)
    real(dp), allocatable :: candidate_temperature(:, :), payload(:)
    logical :: accepted, entity_ok, global_ok, local_ok
    integer :: cell_count, child, ierr, j_lower, j_upper, owner
    integer :: root_maximum, root_minimum, state_count, tile, transfers
    integer :: value_count

    patch_set = reactive_eb_patch_set_2d()
    ok = .false.
    transfers = 0
    if (present(local_transfers)) local_transfers = 0
    local_ok = root >= 0 .and. root < distribution%nranks .and. &
      sparse_patch_set%is_valid( &
        distribution, coarse_geometry, patch_set_template)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Allreduce( &
      root, root_minimum, 1, MPI_INTEGER, MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      root, root_maximum, 1, MPI_INTEGER, MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. root_minimum /= root_maximum) return

    if (distribution%rank == root) then
      allocate(candidate_state( &
        sparse_patch_set%nvar, coarse_geometry%nx, coarse_geometry%ny))
      allocate(candidate_temperature( &
        coarse_geometry%nx, coarse_geometry%ny))
      candidate_patch_set = patch_set_template
    end if
    do tile = 1, distribution%root_tile_count()
      owner = distribution%root_tiles(tile)%owner
      j_lower = distribution%root_tiles(tile)%j_lower
      j_upper = distribution%root_tiles(tile)%j_upper
      cell_count = distribution%root_tiles(tile)%cell_count
      state_count = sparse_patch_set%nvar * cell_count
      value_count = state_count + cell_count
      if (owner == root) then
        if (distribution%rank == root) then
          candidate_state(:, :, j_lower:j_upper) = &
            sparse_patch_set%root_tiles(tile)%state
          candidate_temperature(:, j_lower:j_upper) = &
            sparse_patch_set%root_tiles(tile)%temperature
        end if
      else if (distribution%rank == owner) then
        allocate(payload(value_count))
        payload(1:state_count) = reshape( &
          sparse_patch_set%root_tiles(tile)%state, [state_count])
        payload(state_count + 1:value_count) = reshape( &
          sparse_patch_set%root_tiles(tile)%temperature, [cell_count])
        call MPI_Send( &
          payload, value_count, MPI_DOUBLE_PRECISION, root, &
          sparse_root_materialization_tag, distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS) return
        transfers = transfers + 1
      else if (distribution%rank == root) then
        allocate(payload(value_count))
        call MPI_Recv( &
          payload, value_count, MPI_DOUBLE_PRECISION, owner, &
          sparse_root_materialization_tag, distribution%comm, status, ierr)
        if (ierr /= MPI_SUCCESS) return
        candidate_state(:, :, j_lower:j_upper) = reshape( &
          payload(1:state_count), [sparse_patch_set%nvar, &
            coarse_geometry%nx, j_upper - j_lower + 1])
        candidate_temperature(:, j_lower:j_upper) = reshape( &
          payload(state_count + 1:value_count), &
          [coarse_geometry%nx, j_upper - j_lower + 1])
      end if
      if (allocated(payload)) deallocate(payload)
    end do

    do child = 1, distribution%child_count()
      owner = distribution%child_owner(child)
      cell_count = distribution%child_cell_counts(child)
      state_count = sparse_patch_set%nvar * cell_count
      value_count = state_count + cell_count
      if (owner == root) then
        if (distribution%rank == root) then
          candidate_patch_set%children(child)%state = &
            sparse_patch_set%children(child)%state
          candidate_patch_set%children(child)%temperature = &
            sparse_patch_set%children(child)%temperature
        end if
      else if (distribution%rank == owner) then
        allocate(payload(value_count))
        payload(1:state_count) = reshape( &
          sparse_patch_set%children(child)%state, [state_count])
        payload(state_count + 1:value_count) = reshape( &
          sparse_patch_set%children(child)%temperature, [cell_count])
        call MPI_Send( &
          payload, value_count, MPI_DOUBLE_PRECISION, root, &
          sparse_root_materialization_tag, distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS) return
        transfers = transfers + 1
      else if (distribution%rank == root) then
        allocate(payload(value_count))
        call MPI_Recv( &
          payload, value_count, MPI_DOUBLE_PRECISION, owner, &
          sparse_root_materialization_tag, distribution%comm, status, ierr)
        if (ierr /= MPI_SUCCESS) return
        candidate_patch_set%children(child)%state = reshape( &
          payload(1:state_count), shape( &
            candidate_patch_set%children(child)%state))
        candidate_patch_set%children(child)%temperature = reshape( &
          payload(state_count + 1:value_count), shape( &
            candidate_patch_set%children(child)%temperature))
      end if
      if (allocated(payload)) deallocate(payload)
    end do

    entity_ok = .true.
    if (distribution%rank == root) entity_ok = &
      all(ieee_is_finite(candidate_state)) .and. &
      all(ieee_is_finite(candidate_temperature)) .and. &
      candidate_patch_set%is_valid(coarse_geometry, sparse_patch_set%nvar)
    call all_ranks_accept_eb_2d( &
      distribution, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    if (distribution%rank == root) then
      call move_alloc(candidate_state, coarse_state)
      call move_alloc(candidate_temperature, coarse_temperature)
      patch_set = candidate_patch_set
    end if
    ok = .true.
    if (present(local_transfers)) local_transfers = transfers
  end subroutine gather_sparse_owned_reactive_eb_patch_set_to_root_2d

  subroutine scatter_root_reactive_eb_patch_set_to_sparse_2d( &
      distribution, nspecies, root_coarse_state, root_coarse_temperature, &
      coarse_geometry, root_patch_set, patch_set_template, root, &
      sparse_patch_set, ok, local_transfers)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    integer, intent(in) :: nspecies
    real(dp), allocatable, intent(in) :: root_coarse_state(:, :, :)
    real(dp), allocatable, intent(in) :: root_coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: root_patch_set
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set_template
    integer, intent(in) :: root
    type(mpi_amr_eb_sparse_patch_set_2d), intent(out) :: sparse_patch_set
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_transfers

    type(reactive_eb_patch_topology_2d) :: topology
    logical :: extracted
    integer :: transfers

    sparse_patch_set = mpi_amr_eb_sparse_patch_set_2d()
    ok = .false.
    transfers = 0
    if (present(local_transfers)) local_transfers = 0
    extracted = .false.
    if (nspecies >= 1) then
      call extract_reactive_eb_patch_topology_2d( &
        coarse_geometry, reactive_nvar(nspecies), patch_set_template, &
        topology, extracted)
    end if
    if (.not. extracted) topology = reactive_eb_patch_topology_2d()
    call scatter_root_reactive_eb_topology_to_sparse_2d( &
      distribution, nspecies, root_coarse_state, root_coarse_temperature, &
      coarse_geometry, root_patch_set, topology, root, sparse_patch_set, &
      ok, transfers)
    if (ok .and. present(local_transfers)) local_transfers = transfers
  end subroutine scatter_root_reactive_eb_patch_set_to_sparse_2d

  subroutine scatter_root_reactive_eb_topology_to_sparse_2d( &
      distribution, nspecies, root_coarse_state, root_coarse_temperature, &
      coarse_geometry, root_patch_set, topology, root, &
      sparse_patch_set, ok, local_transfers)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    integer, intent(in) :: nspecies
    real(dp), allocatable, intent(in) :: root_coarse_state(:, :, :)
    real(dp), allocatable, intent(in) :: root_coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: root_patch_set
    type(reactive_eb_patch_topology_2d), intent(in) :: topology
    integer, intent(in) :: root
    type(mpi_amr_eb_sparse_patch_set_2d), intent(out) :: sparse_patch_set
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_transfers

    type(mpi_amr_eb_sparse_patch_set_2d) :: candidate
    type(MPI_Status) :: status
    real(dp), allocatable :: payload(:)
    logical :: accepted, global_ok, local_ok
    integer :: cell_count, child, height, ierr, j_lower, j_upper, nvar
    integer :: owner, root_maximum, root_minimum, state_count, tile
    integer :: transfers, value_count

    sparse_patch_set = mpi_amr_eb_sparse_patch_set_2d()
    ok = .false.
    transfers = 0
    if (present(local_transfers)) local_transfers = 0
    nvar = 0
    if (nspecies >= 1) nvar = reactive_nvar(nspecies)
    local_ok = nspecies >= 1 .and. root >= 0 .and. &
      root < distribution%nranks .and. &
      mpi_amr_eb_distribution_matches_topology_2d( &
        distribution, coarse_geometry, topology) .and. &
      topology%is_valid(coarse_geometry)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Allreduce( &
      root, root_minimum, 1, MPI_INTEGER, MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      root, root_maximum, 1, MPI_INTEGER, MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. root_minimum /= root_maximum) return

    if (distribution%rank == root) then
      local_ok = allocated(root_coarse_state) .and. &
        allocated(root_coarse_temperature)
      if (local_ok) then
        local_ok = all(shape(root_coarse_state) == &
          [nvar, coarse_geometry%nx, coarse_geometry%ny]) .and. &
          all(shape(root_coarse_temperature) == &
            [coarse_geometry%nx, coarse_geometry%ny]) .and. &
          all(ieee_is_finite(root_coarse_state)) .and. &
          all(ieee_is_finite(root_coarse_temperature))
      end if
      if (local_ok) local_ok = &
        root_patch_set%is_valid(coarse_geometry, nvar)
      if (local_ok) local_ok = &
        reactive_eb_patch_set_matches_topology_2d(root_patch_set, topology)
    else
      local_ok = .not. allocated(root_coarse_state) .and. &
        .not. allocated(root_coarse_temperature) .and. &
        root_patch_set%patch_count() == 0
    end if
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    candidate%rank = distribution%rank
    candidate%nranks = distribution%nranks
    candidate%nvar = nvar
    allocate(candidate%root_tiles(distribution%root_tile_count()))
    allocate(candidate%children(distribution%child_count()))
    do tile = 1, distribution%root_tile_count()
      owner = distribution%root_tiles(tile)%owner
      j_lower = distribution%root_tiles(tile)%j_lower
      j_upper = distribution%root_tiles(tile)%j_upper
      height = j_upper - j_lower + 1
      cell_count = distribution%root_tiles(tile)%cell_count
      state_count = nvar * cell_count
      value_count = state_count + cell_count
      if (owner == root) then
        if (distribution%rank == root) then
          allocate(candidate%root_tiles(tile)%state, &
            source=root_coarse_state(:, :, j_lower:j_upper))
          allocate(candidate%root_tiles(tile)%temperature, &
            source=root_coarse_temperature(:, j_lower:j_upper))
        end if
      else if (distribution%rank == root) then
        allocate(payload(value_count))
        payload(1:state_count) = reshape( &
          root_coarse_state(:, :, j_lower:j_upper), [state_count])
        payload(state_count + 1:value_count) = reshape( &
          root_coarse_temperature(:, j_lower:j_upper), [cell_count])
        call MPI_Send( &
          payload, value_count, MPI_DOUBLE_PRECISION, owner, &
          sparse_root_restart_scatter_tag, distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS) return
        transfers = transfers + 1
      else if (distribution%rank == owner) then
        allocate(payload(value_count))
        call MPI_Recv( &
          payload, value_count, MPI_DOUBLE_PRECISION, root, &
          sparse_root_restart_scatter_tag, distribution%comm, status, ierr)
        if (ierr /= MPI_SUCCESS) return
        allocate(candidate%root_tiles(tile)%state(nvar, coarse_geometry%nx, &
          height))
        allocate(candidate%root_tiles(tile)%temperature( &
          coarse_geometry%nx, height))
        candidate%root_tiles(tile)%state = reshape( &
          payload(1:state_count), [nvar, coarse_geometry%nx, height])
        candidate%root_tiles(tile)%temperature = reshape( &
          payload(state_count + 1:value_count), &
          [coarse_geometry%nx, height])
      end if
      if (allocated(payload)) deallocate(payload)
    end do

    do child = 1, distribution%child_count()
      owner = distribution%child_owner(child)
      cell_count = distribution%child_cell_counts(child)
      state_count = nvar * cell_count
      value_count = state_count + cell_count
      if (owner == root) then
        if (distribution%rank == root) then
          allocate(candidate%children(child)%state, &
            source=root_patch_set%children(child)%state)
          allocate(candidate%children(child)%temperature, &
            source=root_patch_set%children(child)%temperature)
        end if
      else if (distribution%rank == root) then
        allocate(payload(value_count))
        payload(1:state_count) = reshape( &
          root_patch_set%children(child)%state, [state_count])
        payload(state_count + 1:value_count) = reshape( &
          root_patch_set%children(child)%temperature, [cell_count])
        call MPI_Send( &
          payload, value_count, MPI_DOUBLE_PRECISION, owner, &
          sparse_root_restart_scatter_tag, distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS) return
        transfers = transfers + 1
      else if (distribution%rank == owner) then
        allocate(payload(value_count))
        call MPI_Recv( &
          payload, value_count, MPI_DOUBLE_PRECISION, root, &
          sparse_root_restart_scatter_tag, distribution%comm, status, ierr)
        if (ierr /= MPI_SUCCESS) return
        allocate(candidate%children(child)%state(nvar, &
          topology%children(child)%geometry%nx, &
          topology%children(child)%geometry%ny))
        allocate(candidate%children(child)%temperature( &
          topology%children(child)%geometry%nx, &
          topology%children(child)%geometry%ny))
        candidate%children(child)%state = reshape( &
          payload(1:state_count), shape( &
            candidate%children(child)%state))
        candidate%children(child)%temperature = reshape( &
          payload(state_count + 1:value_count), shape( &
            candidate%children(child)%temperature))
      end if
      if (allocated(payload)) deallocate(payload)
    end do

    local_ok = mpi_amr_eb_sparse_patch_set_matches_topology_2d( &
      candidate, distribution, coarse_geometry, topology)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    sparse_patch_set = candidate
    ok = .true.
    if (present(local_transfers)) local_transfers = transfers
  end subroutine scatter_root_reactive_eb_topology_to_sparse_2d

  pure logical function reactive_eb_patch_set_matches_topology_2d( &
      patch_set, topology) result(matches)
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set
    type(reactive_eb_patch_topology_2d), intent(in) :: topology

    real(dp), parameter :: tolerance = 5.0e3_dp * epsilon(1.0_dp)
    integer :: child

    matches = patch_set%patch_count() == topology%patch_count()
    if (.not. matches) return
    do child = 1, patch_set%patch_count()
      matches = all([ &
        patch_set%children(child)%patch%coarse_i_lower, &
        patch_set%children(child)%patch%coarse_i_upper, &
        patch_set%children(child)%patch%coarse_j_lower, &
        patch_set%children(child)%patch%coarse_j_upper, &
        patch_set%children(child)%patch%refinement_ratio, &
        patch_set%children(child)%geometry%nx, &
        patch_set%children(child)%geometry%ny] == [ &
        topology%children(child)%patch%coarse_i_lower, &
        topology%children(child)%patch%coarse_i_upper, &
        topology%children(child)%patch%coarse_j_lower, &
        topology%children(child)%patch%coarse_j_upper, &
        topology%children(child)%patch%refinement_ratio, &
        topology%children(child)%geometry%nx, &
        topology%children(child)%geometry%ny])
      if (.not. matches) return
      matches = all(abs([ &
        patch_set%children(child)%geometry%x_lower, &
        patch_set%children(child)%geometry%x_upper, &
        patch_set%children(child)%geometry%y_lower, &
        patch_set%children(child)%geometry%y_upper, &
        patch_set%children(child)%geometry%dx, &
        patch_set%children(child)%geometry%dy] - [ &
        topology%children(child)%geometry%x_lower, &
        topology%children(child)%geometry%x_upper, &
        topology%children(child)%geometry%y_lower, &
        topology%children(child)%geometry%y_upper, &
        topology%children(child)%geometry%dx, &
        topology%children(child)%geometry%dy]) <= tolerance) .and. &
        all(patch_set%children(child)%geometry%cell_type == &
          topology%children(child)%geometry%cell_type) .and. &
        all(abs(patch_set%children(child)%geometry%volume_fraction - &
          topology%children(child)%geometry%volume_fraction) <= &
            tolerance) .and. &
        all(abs(patch_set%children(child)%geometry%cell_centroid_x - &
          topology%children(child)%geometry%cell_centroid_x) <= tolerance) &
        .and. all(abs( &
          patch_set%children(child)%geometry%cell_centroid_y - &
          topology%children(child)%geometry%cell_centroid_y) <= tolerance) &
        .and. all(abs( &
          patch_set%children(child)%geometry%x_face_fraction - &
          topology%children(child)%geometry%x_face_fraction) <= tolerance) &
        .and. all(abs( &
          patch_set%children(child)%geometry%y_face_fraction - &
          topology%children(child)%geometry%y_face_fraction) <= tolerance) &
        .and. all(abs( &
          patch_set%children(child)%geometry%x_face_centroid_y - &
          topology%children(child)%geometry%x_face_centroid_y) <= tolerance) &
        .and. all(abs( &
          patch_set%children(child)%geometry%y_face_centroid_x - &
          topology%children(child)%geometry%y_face_centroid_x) <= tolerance) &
        .and. all(abs( &
          patch_set%children(child)%geometry%boundary_length - &
          topology%children(child)%geometry%boundary_length) <= tolerance) &
        .and. all(abs( &
          patch_set%children(child)%geometry%boundary_centroid_x - &
          topology%children(child)%geometry%boundary_centroid_x) <= &
            tolerance) .and. all(abs( &
          patch_set%children(child)%geometry%boundary_centroid_y - &
          topology%children(child)%geometry%boundary_centroid_y) <= &
            tolerance) .and. all(abs( &
          patch_set%children(child)%geometry%boundary_normal_x - &
          topology%children(child)%geometry%boundary_normal_x) <= tolerance) &
        .and. all(abs( &
          patch_set%children(child)%geometry%boundary_normal_y - &
          topology%children(child)%geometry%boundary_normal_y) <= tolerance) &
        .and. all(abs( &
          patch_set%children(child)%geometry%boundary_normal_integral_x - &
          topology%children(child)%geometry%boundary_normal_integral_x) <= &
            tolerance) .and. all(abs( &
          patch_set%children(child)%geometry%boundary_normal_integral_y - &
          topology%children(child)%geometry%boundary_normal_integral_y) <= &
            tolerance)
      if (.not. matches) return
    end do
  end function reactive_eb_patch_set_matches_topology_2d

  subroutine regrid_sparse_owned_reactive_eb_patch_set_2d( &
      species, distribution, sparse_patch_set, coarse_geometry, &
      patch_set_template, new_fine_geometries, new_collection, &
      refinement_ratio, ok, changed, local_restriction_transfers, &
      local_prolongation_transfers, local_overlap_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_distribution_2d), intent(inout) :: distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(inout) :: sparse_patch_set
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(inout) :: patch_set_template
    type(eb_geometry_2d), intent(in) :: new_fine_geometries(:)
    type(amr_eb_regrid_plan_collection_2d), intent(in) :: new_collection
    integer, intent(in) :: refinement_ratio
    logical, intent(out) :: ok, changed
    integer, intent(out), optional :: local_restriction_transfers
    integer, intent(out), optional :: local_prolongation_transfers
    integer, intent(out), optional :: local_overlap_transfers

    type(mpi_amr_eb_patch_distribution_2d) :: candidate_distribution
    type(mpi_amr_eb_sparse_patch_set_2d) :: averaged_sparse_patch_set
    type(mpi_amr_eb_sparse_patch_set_2d) :: candidate_sparse_patch_set
    type(reactive_eb_patch_set_2d) :: candidate_patch_set
    logical :: accepted, global_ok, local_changed, local_ok
    integer :: controls(5), controls_maximum(5), controls_minimum(5)
    integer :: ierr, nvar, overlap_transfers, prolongation_transfers
    integer :: restriction_transfers

    ok = .false.
    changed = .false.
    restriction_transfers = 0
    prolongation_transfers = 0
    overlap_transfers = 0
    if (present(local_restriction_transfers)) &
      local_restriction_transfers = 0
    if (present(local_prolongation_transfers)) &
      local_prolongation_transfers = 0
    if (present(local_overlap_transfers)) local_overlap_transfers = 0
    nvar = reactive_nvar(size(species))
    local_ok = size(species) >= 1 .and. refinement_ratio >= 2 .and. &
      new_collection%is_valid() .and. &
      new_collection%coarse_nx == coarse_geometry%nx .and. &
      new_collection%coarse_ny == coarse_geometry%ny .and. &
      size(new_fine_geometries) == new_collection%patch_count() .and. &
      sparse_patch_set%nvar == nvar .and. &
      sparse_patch_set%is_valid( &
        distribution, coarse_geometry, patch_set_template)
    if (local_ok) local_ok = all_sparse_regrid_geometries_valid_2d( &
      new_fine_geometries)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    controls = [ &
      refinement_ratio, new_collection%coarse_nx, &
      new_collection%coarse_ny, new_collection%tagged_cell_count, &
      new_collection%patch_count()]
    call MPI_Allreduce( &
      controls, controls_minimum, 5, MPI_INTEGER, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      controls, controls_maximum, 5, MPI_INTEGER, MPI_MAX, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. &
        any(controls_minimum /= controls_maximum)) return

    call build_sparse_regrid_template_2d( &
      nvar, coarse_geometry, new_fine_geometries, new_collection, &
      refinement_ratio, candidate_patch_set, local_ok)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call initialize_mpi_amr_eb_patch_distribution_2d( &
      coarse_geometry, candidate_patch_set, distribution%comm, &
      candidate_distribution, local_ok, distribution%subcycle_exponent)
    if (.not. local_ok) return

    averaged_sparse_patch_set = sparse_patch_set
    call average_down_sparse_owned_reactive_eb_patch_set_2d( &
      species, distribution, averaged_sparse_patch_set, coarse_geometry, &
      patch_set_template, local_ok, restriction_transfers)
    if (.not. local_ok) return

    call initialize_direct_sparse_regrid_2d( &
      species, distribution, averaged_sparse_patch_set, coarse_geometry, &
      candidate_distribution, candidate_patch_set, &
      candidate_sparse_patch_set, prolongation_transfers, local_ok)
    if (.not. local_ok) return
    call transfer_direct_sparse_regrid_overlaps_2d( &
      species, distribution, averaged_sparse_patch_set, patch_set_template, &
      candidate_distribution, candidate_sparse_patch_set, &
      candidate_patch_set, coarse_geometry, overlap_transfers, local_ok)
    if (.not. local_ok) return
    local_changed = sparse_regrid_topology_changed_2d( &
      patch_set_template, candidate_patch_set)

    local_ok = candidate_sparse_patch_set%is_valid( &
      candidate_distribution, coarse_geometry, candidate_patch_set)
    call all_ranks_accept_eb_2d( &
      candidate_distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    distribution = candidate_distribution
    sparse_patch_set = candidate_sparse_patch_set
    patch_set_template = candidate_patch_set
    changed = local_changed
    ok = .true.
    if (present(local_restriction_transfers)) &
      local_restriction_transfers = restriction_transfers
    if (present(local_prolongation_transfers)) &
      local_prolongation_transfers = prolongation_transfers
    if (present(local_overlap_transfers)) &
      local_overlap_transfers = overlap_transfers
  end subroutine regrid_sparse_owned_reactive_eb_patch_set_2d

  pure logical function all_sparse_regrid_geometries_valid_2d(geometries) &
      result(valid)
    type(eb_geometry_2d), intent(in) :: geometries(:)

    integer :: child

    valid = .true.
    do child = 1, size(geometries)
      valid = geometries(child)%is_valid()
      if (.not. valid) return
    end do
  end function all_sparse_regrid_geometries_valid_2d

  subroutine build_sparse_regrid_template_2d( &
      nvar, coarse_geometry, fine_geometries, collection, refinement_ratio, &
      patch_set, ok)
    integer, intent(in) :: nvar, refinement_ratio
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(eb_geometry_2d), intent(in) :: fine_geometries(:)
    type(amr_eb_regrid_plan_collection_2d), intent(in) :: collection
    type(reactive_eb_patch_set_2d), intent(out) :: patch_set
    logical, intent(out) :: ok

    type(reactive_eb_patch_set_2d) :: candidate
    logical :: local_ok
    integer :: child

    patch_set = reactive_eb_patch_set_2d()
    ok = .false.
    if (nvar < 1 .or. refinement_ratio < 2 .or. &
        .not. coarse_geometry%is_valid() .or. &
        .not. collection%is_valid() .or. &
        collection%coarse_nx /= coarse_geometry%nx .or. &
        collection%coarse_ny /= coarse_geometry%ny .or. &
        size(fine_geometries) /= collection%patch_count() .or. &
        .not. all_sparse_regrid_geometries_valid_2d(fine_geometries)) return

    allocate(candidate%children(collection%patch_count()))
    do child = 1, collection%patch_count()
      candidate%children(child)%geometry = fine_geometries(child)
      call build_amr_eb_patch_2d( &
        coarse_geometry, fine_geometries(child), &
        collection%plans(child)%coarse_i_lower, &
        collection%plans(child)%coarse_i_upper, &
        collection%plans(child)%coarse_j_lower, &
        collection%plans(child)%coarse_j_upper, refinement_ratio, &
        candidate%children(child)%patch, local_ok)
      if (.not. local_ok) return
      allocate(candidate%children(child)%state( &
        nvar, fine_geometries(child)%nx, fine_geometries(child)%ny), &
        source=0.0_dp)
      allocate(candidate%children(child)%temperature( &
        fine_geometries(child)%nx, fine_geometries(child)%ny), source=1.0_dp)
    end do
    if (.not. candidate%is_valid(coarse_geometry, nvar)) return
    patch_set = candidate
    ok = .true.
  end subroutine build_sparse_regrid_template_2d

  pure logical function sparse_regrid_topology_changed_2d( &
      old_patch_set, new_patch_set) result(changed)
    type(reactive_eb_patch_set_2d), intent(in) :: old_patch_set
    type(reactive_eb_patch_set_2d), intent(in) :: new_patch_set

    integer :: child

    changed = old_patch_set%patch_count() /= new_patch_set%patch_count()
    if (changed) return
    do child = 1, new_patch_set%patch_count()
      changed = any([ &
        old_patch_set%children(child)%patch%coarse_i_lower, &
        old_patch_set%children(child)%patch%coarse_i_upper, &
        old_patch_set%children(child)%patch%coarse_j_lower, &
        old_patch_set%children(child)%patch%coarse_j_upper, &
        old_patch_set%children(child)%patch%refinement_ratio] /= [ &
        new_patch_set%children(child)%patch%coarse_i_lower, &
        new_patch_set%children(child)%patch%coarse_i_upper, &
        new_patch_set%children(child)%patch%coarse_j_lower, &
        new_patch_set%children(child)%patch%coarse_j_upper, &
        new_patch_set%children(child)%patch%refinement_ratio])
      if (changed) return
    end do
  end function sparse_regrid_topology_changed_2d

  pure logical function sparse_regrid_root_layouts_match_2d( &
      old_distribution, new_distribution) result(matches)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: old_distribution
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: new_distribution

    integer :: tile

    matches = old_distribution%rank == new_distribution%rank .and. &
      old_distribution%nranks == new_distribution%nranks .and. &
      old_distribution%root_tile_count() == &
        new_distribution%root_tile_count()
    if (.not. matches) return
    do tile = 1, old_distribution%root_tile_count()
      matches = all([ &
        old_distribution%root_tiles(tile)%owner, &
        old_distribution%root_tiles(tile)%i_lower, &
        old_distribution%root_tiles(tile)%i_upper, &
        old_distribution%root_tiles(tile)%j_lower, &
        old_distribution%root_tiles(tile)%j_upper] == [ &
        new_distribution%root_tiles(tile)%owner, &
        new_distribution%root_tiles(tile)%i_lower, &
        new_distribution%root_tiles(tile)%i_upper, &
        new_distribution%root_tiles(tile)%j_lower, &
        new_distribution%root_tiles(tile)%j_upper])
      if (.not. matches) return
    end do
  end function sparse_regrid_root_layouts_match_2d

  subroutine initialize_direct_sparse_regrid_2d( &
      species, old_distribution, averaged_patch_set, coarse_geometry, &
      new_distribution, new_patch_set, new_sparse_patch_set, &
      local_transfers, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: old_distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(in) :: averaged_patch_set
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: new_distribution
    type(reactive_eb_patch_set_2d), intent(in) :: new_patch_set
    type(mpi_amr_eb_sparse_patch_set_2d), intent(out) :: new_sparse_patch_set
    integer, intent(out) :: local_transfers
    logical, intent(out) :: ok

    type(mpi_amr_eb_sparse_patch_set_2d) :: candidate
    type(MPI_Status) :: status
    real(dp), allocatable :: payload(:), root_state(:, :, :)
    real(dp), allocatable :: root_temperature(:, :)
    logical, allocatable :: recipients(:)
    logical :: accepted, entity_ok, global_ok, local_ok
    integer :: cell_count, child, ierr, j_lower, j_upper, nvar, owner
    integer :: recipient, state_count, tile, value_count

    new_sparse_patch_set = mpi_amr_eb_sparse_patch_set_2d()
    local_transfers = 0
    ok = .false.
    nvar = reactive_nvar(size(species))
    local_ok = size(species) >= 1 .and. averaged_patch_set%nvar == nvar .and. &
      sparse_regrid_root_layouts_match_2d( &
        old_distribution, new_distribution)
    call all_ranks_accept_eb_2d( &
      old_distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    candidate%rank = new_distribution%rank
    candidate%nranks = new_distribution%nranks
    candidate%nvar = nvar
    allocate(candidate%root_tiles(new_distribution%root_tile_count()))
    allocate(candidate%children(new_distribution%child_count()))
    do tile = 1, new_distribution%root_tile_count()
      if (.not. new_distribution%root_tile_is_local(tile)) cycle
      allocate(candidate%root_tiles(tile)%state, &
        source=averaged_patch_set%root_tiles(tile)%state)
      allocate(candidate%root_tiles(tile)%temperature, &
        source=averaged_patch_set%root_tiles(tile)%temperature)
    end do

    allocate(recipients(new_distribution%nranks), source=.false.)
    do child = 1, new_distribution%child_count()
      recipient = new_distribution%child_owner(child)
      recipients(recipient + 1) = .true.
    end do
    if (recipients(new_distribution%rank + 1)) then
      allocate(root_state(nvar, coarse_geometry%nx, coarse_geometry%ny))
      allocate(root_temperature(coarse_geometry%nx, coarse_geometry%ny))
    end if

    do recipient = 0, new_distribution%nranks - 1
      if (.not. recipients(recipient + 1)) cycle
      do tile = 1, old_distribution%root_tile_count()
        owner = old_distribution%root_tiles(tile)%owner
        j_lower = old_distribution%root_tiles(tile)%j_lower
        j_upper = old_distribution%root_tiles(tile)%j_upper
        cell_count = old_distribution%root_tiles(tile)%cell_count
        state_count = nvar * cell_count
        value_count = state_count + cell_count
        if (owner == recipient) then
          if (old_distribution%rank == recipient) then
            root_state(:, :, j_lower:j_upper) = &
              averaged_patch_set%root_tiles(tile)%state
            root_temperature(:, j_lower:j_upper) = &
              averaged_patch_set%root_tiles(tile)%temperature
          end if
        else if (old_distribution%rank == owner) then
          allocate(payload(value_count))
          payload(1:state_count) = reshape( &
            averaged_patch_set%root_tiles(tile)%state, [state_count])
          payload(state_count + 1:value_count) = reshape( &
            averaged_patch_set%root_tiles(tile)%temperature, [cell_count])
          call MPI_Send( &
            payload, value_count, MPI_DOUBLE_PRECISION, recipient, &
            sparse_regrid_prolongation_tag, old_distribution%comm, ierr)
          deallocate(payload)
          if (ierr /= MPI_SUCCESS) return
          local_transfers = local_transfers + 1
        else if (old_distribution%rank == recipient) then
          allocate(payload(value_count))
          call MPI_Recv( &
            payload, value_count, MPI_DOUBLE_PRECISION, owner, &
            sparse_regrid_prolongation_tag, old_distribution%comm, status, &
            ierr)
          if (ierr /= MPI_SUCCESS) return
          root_state(:, :, j_lower:j_upper) = reshape( &
            payload(1:state_count), [nvar, coarse_geometry%nx, &
              j_upper - j_lower + 1])
          root_temperature(:, j_lower:j_upper) = reshape( &
            payload(state_count + 1:value_count), &
            [coarse_geometry%nx, j_upper - j_lower + 1])
          deallocate(payload)
        end if
      end do
    end do

    entity_ok = .true.
    if (recipients(new_distribution%rank + 1)) entity_ok = &
      all(ieee_is_finite(root_state)) .and. &
      all(ieee_is_finite(root_temperature))
    if (entity_ok) then
      do child = 1, new_distribution%child_count()
        if (.not. new_distribution%child_is_local(child)) cycle
        allocate(candidate%children(child)%state( &
          nvar, new_patch_set%children(child)%geometry%nx, &
          new_patch_set%children(child)%geometry%ny))
        allocate(candidate%children(child)%temperature( &
          new_patch_set%children(child)%geometry%nx, &
          new_patch_set%children(child)%geometry%ny))
        call prolong_reactive_eb_patch_pcm_2d( &
          species, root_state, root_temperature, coarse_geometry, &
          new_patch_set%children(child)%geometry, &
          new_patch_set%children(child)%patch, &
          candidate%children(child)%state, &
          candidate%children(child)%temperature, local_ok)
        entity_ok = entity_ok .and. local_ok
      end do
    end if
    call all_ranks_accept_eb_2d( &
      new_distribution, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    local_ok = candidate%is_valid( &
      new_distribution, coarse_geometry, new_patch_set)
    call all_ranks_accept_eb_2d( &
      new_distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    new_sparse_patch_set = candidate
    ok = .true.
  end subroutine initialize_direct_sparse_regrid_2d

  pure logical function sparse_regrid_overlap_geometry_matches_2d( &
      old_patch_set, new_patch_set, old_child, new_child, &
      coarse_i_lower, coarse_i_upper, coarse_j_lower, coarse_j_upper, &
      ratio) result(matches)
    type(reactive_eb_patch_set_2d), intent(in) :: old_patch_set
    type(reactive_eb_patch_set_2d), intent(in) :: new_patch_set
    integer, intent(in) :: old_child, new_child, coarse_i_lower
    integer, intent(in) :: coarse_i_upper, coarse_j_lower, coarse_j_upper
    integer, intent(in) :: ratio

    real(dp), parameter :: geometry_tolerance = &
      5.0e3_dp * epsilon(1.0_dp)
    integer :: old_i_lower, old_i_upper, old_j_lower, old_j_upper
    integer :: new_i_lower, new_i_upper, new_j_lower, new_j_upper

    old_i_lower = (coarse_i_lower - &
      old_patch_set%children(old_child)%patch%coarse_i_lower) * ratio + 1
    old_i_upper = old_i_lower + &
      (coarse_i_upper - coarse_i_lower + 1) * ratio - 1
    old_j_lower = (coarse_j_lower - &
      old_patch_set%children(old_child)%patch%coarse_j_lower) * ratio + 1
    old_j_upper = old_j_lower + &
      (coarse_j_upper - coarse_j_lower + 1) * ratio - 1
    new_i_lower = (coarse_i_lower - &
      new_patch_set%children(new_child)%patch%coarse_i_lower) * ratio + 1
    new_i_upper = new_i_lower + &
      (coarse_i_upper - coarse_i_lower + 1) * ratio - 1
    new_j_lower = (coarse_j_lower - &
      new_patch_set%children(new_child)%patch%coarse_j_lower) * ratio + 1
    new_j_upper = new_j_lower + &
      (coarse_j_upper - coarse_j_lower + 1) * ratio - 1
    matches = all( &
      old_patch_set%children(old_child)%geometry%cell_type( &
        old_i_lower:old_i_upper, old_j_lower:old_j_upper) == &
      new_patch_set%children(new_child)%geometry%cell_type( &
        new_i_lower:new_i_upper, new_j_lower:new_j_upper)) .and. &
      all(abs( &
        old_patch_set%children(old_child)%geometry%volume_fraction( &
          old_i_lower:old_i_upper, old_j_lower:old_j_upper) - &
        new_patch_set%children(new_child)%geometry%volume_fraction( &
          new_i_lower:new_i_upper, new_j_lower:new_j_upper)) <= &
        geometry_tolerance)
  end function sparse_regrid_overlap_geometry_matches_2d

  subroutine transfer_sparse_regrid_overlap_rectangle_2d( &
      old_distribution, new_distribution, old_sparse_patch_set, &
      new_sparse_patch_set, old_patch_set, new_patch_set, old_child, &
      new_child, coarse_i_lower, coarse_i_upper, coarse_j_lower, &
      coarse_j_upper, ratio, local_transfers, ok)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: old_distribution
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: new_distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(in) :: old_sparse_patch_set
    type(mpi_amr_eb_sparse_patch_set_2d), intent(inout) :: &
      new_sparse_patch_set
    type(reactive_eb_patch_set_2d), intent(in) :: old_patch_set
    type(reactive_eb_patch_set_2d), intent(in) :: new_patch_set
    integer, intent(in) :: old_child, new_child, coarse_i_lower
    integer, intent(in) :: coarse_i_upper, coarse_j_lower, coarse_j_upper
    integer, intent(in) :: ratio
    integer, intent(inout) :: local_transfers
    logical, intent(out) :: ok

    type(MPI_Status) :: status
    real(dp), allocatable :: payload(:)
    integer :: cell_count, fine_height, fine_width, ierr, new_i_lower
    integer :: new_i_upper, new_j_lower, new_j_upper, new_owner, old_i_lower
    integer :: old_i_upper, old_j_lower, old_j_upper, old_owner, state_count
    integer :: value_count

    ok = .false.
    fine_width = (coarse_i_upper - coarse_i_lower + 1) * ratio
    fine_height = (coarse_j_upper - coarse_j_lower + 1) * ratio
    if (fine_width < 1 .or. fine_height < 1) return
    old_i_lower = (coarse_i_lower - &
      old_patch_set%children(old_child)%patch%coarse_i_lower) * ratio + 1
    old_i_upper = old_i_lower + fine_width - 1
    old_j_lower = (coarse_j_lower - &
      old_patch_set%children(old_child)%patch%coarse_j_lower) * ratio + 1
    old_j_upper = old_j_lower + fine_height - 1
    new_i_lower = (coarse_i_lower - &
      new_patch_set%children(new_child)%patch%coarse_i_lower) * ratio + 1
    new_i_upper = new_i_lower + fine_width - 1
    new_j_lower = (coarse_j_lower - &
      new_patch_set%children(new_child)%patch%coarse_j_lower) * ratio + 1
    new_j_upper = new_j_lower + fine_height - 1
    old_owner = old_distribution%child_owner(old_child)
    new_owner = new_distribution%child_owner(new_child)

    if (old_owner == new_owner) then
      if (old_distribution%rank == old_owner) then
        new_sparse_patch_set%children(new_child)%state( &
          :, new_i_lower:new_i_upper, new_j_lower:new_j_upper) = &
          old_sparse_patch_set%children(old_child)%state( &
            :, old_i_lower:old_i_upper, old_j_lower:old_j_upper)
        new_sparse_patch_set%children(new_child)%temperature( &
          new_i_lower:new_i_upper, new_j_lower:new_j_upper) = &
          old_sparse_patch_set%children(old_child)%temperature( &
            old_i_lower:old_i_upper, old_j_lower:old_j_upper)
      end if
      ok = .true.
      return
    end if

    cell_count = fine_width * fine_height
    state_count = old_sparse_patch_set%nvar * cell_count
    value_count = state_count + cell_count
    if (old_distribution%rank == old_owner) then
      allocate(payload(value_count))
      payload(1:state_count) = reshape( &
        old_sparse_patch_set%children(old_child)%state( &
          :, old_i_lower:old_i_upper, old_j_lower:old_j_upper), &
        [state_count])
      payload(state_count + 1:value_count) = reshape( &
        old_sparse_patch_set%children(old_child)%temperature( &
          old_i_lower:old_i_upper, old_j_lower:old_j_upper), [cell_count])
      call MPI_Send( &
        payload, value_count, MPI_DOUBLE_PRECISION, new_owner, &
        sparse_regrid_overlap_tag, old_distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
      local_transfers = local_transfers + 1
    else if (old_distribution%rank == new_owner) then
      allocate(payload(value_count))
      call MPI_Recv( &
        payload, value_count, MPI_DOUBLE_PRECISION, old_owner, &
        sparse_regrid_overlap_tag, old_distribution%comm, status, ierr)
      if (ierr /= MPI_SUCCESS) return
      new_sparse_patch_set%children(new_child)%state( &
        :, new_i_lower:new_i_upper, new_j_lower:new_j_upper) = reshape( &
          payload(1:state_count), &
          [old_sparse_patch_set%nvar, fine_width, fine_height])
      new_sparse_patch_set%children(new_child)%temperature( &
        new_i_lower:new_i_upper, new_j_lower:new_j_upper) = reshape( &
          payload(state_count + 1:value_count), [fine_width, fine_height])
    end if
    ok = .true.
  end subroutine transfer_sparse_regrid_overlap_rectangle_2d

  subroutine transfer_direct_sparse_regrid_overlaps_2d( &
      species, old_distribution, old_sparse_patch_set, old_patch_set, &
      new_distribution, new_sparse_patch_set, new_patch_set, &
      coarse_geometry, local_transfers, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: old_distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(in) :: old_sparse_patch_set
    type(reactive_eb_patch_set_2d), intent(in) :: old_patch_set
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: new_distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(inout) :: &
      new_sparse_patch_set
    type(reactive_eb_patch_set_2d), intent(in) :: new_patch_set
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    integer, intent(out) :: local_transfers
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: recovered_temperature, sound_speed
    logical :: accepted, global_ok, local_ok
    integer :: coarse_i_lower, coarse_i_upper, coarse_j_lower
    integer :: coarse_j_upper, i, j, new_child, old_child, ratio

    local_transfers = 0
    ok = .false.
    local_ok = old_sparse_patch_set%is_valid( &
      old_distribution, coarse_geometry, old_patch_set) .and. &
      new_sparse_patch_set%is_valid( &
        new_distribution, coarse_geometry, new_patch_set)
    call all_ranks_accept_eb_2d( &
      old_distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    local_ok = .true.
    do new_child = 1, new_patch_set%patch_count()
      ratio = new_patch_set%children(new_child)%patch%refinement_ratio
      do old_child = 1, old_patch_set%patch_count()
        if (old_patch_set%children(old_child)%patch%refinement_ratio /= &
            ratio) cycle
        coarse_i_lower = max( &
          old_patch_set%children(old_child)%patch%coarse_i_lower, &
          new_patch_set%children(new_child)%patch%coarse_i_lower)
        coarse_i_upper = min( &
          old_patch_set%children(old_child)%patch%coarse_i_upper, &
          new_patch_set%children(new_child)%patch%coarse_i_upper)
        coarse_j_lower = max( &
          old_patch_set%children(old_child)%patch%coarse_j_lower, &
          new_patch_set%children(new_child)%patch%coarse_j_lower)
        coarse_j_upper = min( &
          old_patch_set%children(old_child)%patch%coarse_j_upper, &
          new_patch_set%children(new_child)%patch%coarse_j_upper)
        if (coarse_i_lower > coarse_i_upper .or. &
            coarse_j_lower > coarse_j_upper) cycle
        local_ok = local_ok .and. &
          sparse_regrid_overlap_geometry_matches_2d( &
            old_patch_set, new_patch_set, old_child, new_child, &
            coarse_i_lower, coarse_i_upper, coarse_j_lower, &
            coarse_j_upper, ratio)
      end do
    end do
    call all_ranks_accept_eb_2d( &
      old_distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    do new_child = 1, new_patch_set%patch_count()
      ratio = new_patch_set%children(new_child)%patch%refinement_ratio
      do old_child = 1, old_patch_set%patch_count()
        if (old_patch_set%children(old_child)%patch%refinement_ratio /= &
            ratio) cycle
        coarse_i_lower = max( &
          old_patch_set%children(old_child)%patch%coarse_i_lower, &
          new_patch_set%children(new_child)%patch%coarse_i_lower)
        coarse_i_upper = min( &
          old_patch_set%children(old_child)%patch%coarse_i_upper, &
          new_patch_set%children(new_child)%patch%coarse_i_upper)
        coarse_j_lower = max( &
          old_patch_set%children(old_child)%patch%coarse_j_lower, &
          new_patch_set%children(new_child)%patch%coarse_j_lower)
        coarse_j_upper = min( &
          old_patch_set%children(old_child)%patch%coarse_j_upper, &
          new_patch_set%children(new_child)%patch%coarse_j_upper)
        if (coarse_i_lower > coarse_i_upper .or. &
            coarse_j_lower > coarse_j_upper) cycle
        call transfer_sparse_regrid_overlap_rectangle_2d( &
          old_distribution, new_distribution, old_sparse_patch_set, &
          new_sparse_patch_set, old_patch_set, new_patch_set, old_child, &
          new_child, coarse_i_lower, coarse_i_upper, coarse_j_lower, &
          coarse_j_upper, ratio, local_transfers, local_ok)
        if (.not. local_ok) return
      end do
    end do

    allocate(primitive(reactive_nprim(size(species))))
    local_ok = .true.
    do new_child = 1, new_distribution%child_count()
      if (.not. new_distribution%child_is_local(new_child)) cycle
      do j = 1, new_patch_set%children(new_child)%geometry%ny
        do i = 1, new_patch_set%children(new_child)%geometry%nx
          if (new_patch_set%children(new_child)%geometry%cell_type(i, j) == &
              eb_covered_cell) cycle
          if (new_sparse_patch_set%children(new_child)%temperature(i, j) &
              <= 0.0_dp) then
            local_ok = .false.
            cycle
          end if
          call reactive_conserved_to_primitive( &
            species, new_sparse_patch_set%children(new_child)%state(:, i, j), &
            new_sparse_patch_set%children(new_child)%temperature(i, j), &
            primitive, recovered_temperature, sound_speed, accepted)
          if (.not. accepted) then
            local_ok = .false.
            cycle
          end if
          new_sparse_patch_set%children(new_child)%temperature(i, j) = &
            recovered_temperature
        end do
      end do
    end do
    call all_ranks_accept_eb_2d( &
      new_distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    local_ok = new_sparse_patch_set%is_valid( &
      new_distribution, coarse_geometry, new_patch_set)
    call all_ranks_accept_eb_2d( &
      new_distribution, local_ok, accepted, global_ok)
    ok = global_ok .and. accepted
  end subroutine transfer_direct_sparse_regrid_overlaps_2d

  subroutine regrid_tagged_sparse_owned_reactive_eb_patch_set_2d( &
      species, distribution, sparse_patch_set, coarse_geometry, &
      patch_set_template, criteria, refinement_ratio, geometry_builder, &
      ok, changed, local_root_transfers, local_restriction_transfers, &
      local_prolongation_transfers, local_overlap_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_distribution_2d), intent(inout) :: distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(inout) :: sparse_patch_set
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(inout) :: patch_set_template
    type(amr_eb_tagging_criteria_2d), intent(in) :: criteria
    integer, intent(in) :: refinement_ratio
    procedure(sparse_eb_geometry_builder_2d) :: geometry_builder
    logical, intent(out) :: ok, changed
    integer, intent(out), optional :: local_root_transfers
    integer, intent(out), optional :: local_restriction_transfers
    integer, intent(out), optional :: local_prolongation_transfers
    integer, intent(out), optional :: local_overlap_transfers

    type(amr_eb_regrid_plan_collection_2d) :: collection
    type(eb_geometry_2d), allocatable :: fine_geometries(:)
    real(dp), allocatable :: root_state(:, :, :)
    real(dp), allocatable :: root_temperature(:, :)
    real(dp) :: numeric_controls(3), numeric_maximum(3), numeric_minimum(3)
    logical, allocatable :: tags(:, :)
    logical :: accepted, entity_ok, global_ok, local_ok
    integer, allocatable :: metadata(:)
    integer :: child, header(4), ierr, index, integer_controls(5)
    integer :: integer_maximum(5), integer_minimum(5), patch_count
    integer :: overlap_transfers, prolongation_transfers
    integer :: restriction_transfers, root_owner, transfers

    ok = .false.
    changed = .false.
    transfers = 0
    restriction_transfers = 0
    prolongation_transfers = 0
    overlap_transfers = 0
    if (present(local_root_transfers)) local_root_transfers = 0
    if (present(local_restriction_transfers)) &
      local_restriction_transfers = 0
    if (present(local_prolongation_transfers)) &
      local_prolongation_transfers = 0
    if (present(local_overlap_transfers)) local_overlap_transfers = 0
    numeric_controls = [ &
      criteria%relative_gradient_threshold, &
      criteria%absolute_gradient_threshold, criteria%scale_floor]
    integer_controls = [ &
      refinement_ratio, criteria%buffer_cells, &
      criteria%minimum_patch_cells_x, criteria%minimum_patch_cells_y, &
      criteria%maximum_patch_gap_cells]
    local_ok = size(species) >= 1 .and. refinement_ratio >= 2 .and. &
      criteria%is_valid(coarse_geometry%nx, coarse_geometry%ny) .and. &
      sparse_patch_set%is_valid( &
        distribution, coarse_geometry, patch_set_template)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Allreduce( &
      numeric_controls, numeric_minimum, 3, MPI_DOUBLE_PRECISION, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      numeric_controls, numeric_maximum, 3, MPI_DOUBLE_PRECISION, MPI_MAX, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_controls, integer_minimum, 5, MPI_INTEGER, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_controls, integer_maximum, 5, MPI_INTEGER, MPI_MAX, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. &
        any(numeric_minimum /= numeric_maximum) .or. &
        any(integer_minimum /= integer_maximum)) return

    call gather_sparse_root_to_owner_2d( &
      distribution, sparse_patch_set, coarse_geometry, root_state, &
      root_temperature, transfers, local_ok)
    if (.not. local_ok) return
    root_owner = distribution%root_level_owner()
    entity_ok = root_owner >= 0 .and. root_owner < distribution%nranks
    if (distribution%rank == root_owner .and. entity_ok) then
      allocate(tags(coarse_geometry%nx, coarse_geometry%ny))
      call plan_reactive_eb_temperature_regrid_collection_2d( &
        root_temperature, coarse_geometry, criteria, tags, collection, &
        entity_ok)
    end if
    call all_ranks_accept_eb_2d( &
      distribution, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    if (distribution%rank == root_owner) header = [ &
      collection%coarse_nx, collection%coarse_ny, &
      collection%tagged_cell_count, collection%patch_count()]
    call MPI_Bcast( &
      header, 4, MPI_INTEGER, root_owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    patch_count = header(4)
    if (distribution%rank /= root_owner) then
      collection%coarse_nx = header(1)
      collection%coarse_ny = header(2)
      collection%tagged_cell_count = header(3)
      allocate(collection%plans(patch_count))
    end if
    allocate(metadata(12 * patch_count))
    if (distribution%rank == root_owner) then
      index = 1
      do child = 1, patch_count
        metadata(index:index + 11) = [ &
          merge(1, 0, collection%plans(child)%active), &
          collection%plans(child)%coarse_nx, &
          collection%plans(child)%coarse_ny, &
          collection%plans(child)%tagged_cell_count, &
          collection%plans(child)%tag_i_lower, &
          collection%plans(child)%tag_i_upper, &
          collection%plans(child)%tag_j_lower, &
          collection%plans(child)%tag_j_upper, &
          collection%plans(child)%coarse_i_lower, &
          collection%plans(child)%coarse_i_upper, &
          collection%plans(child)%coarse_j_lower, &
          collection%plans(child)%coarse_j_upper]
        index = index + 12
      end do
    end if
    if (patch_count > 0) then
      call MPI_Bcast( &
        metadata, size(metadata), MPI_INTEGER, root_owner, &
        distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
    end if
    if (distribution%rank /= root_owner) then
      index = 1
      do child = 1, patch_count
        collection%plans(child)%active = metadata(index) == 1
        collection%plans(child)%coarse_nx = metadata(index + 1)
        collection%plans(child)%coarse_ny = metadata(index + 2)
        collection%plans(child)%tagged_cell_count = metadata(index + 3)
        collection%plans(child)%tag_i_lower = metadata(index + 4)
        collection%plans(child)%tag_i_upper = metadata(index + 5)
        collection%plans(child)%tag_j_lower = metadata(index + 6)
        collection%plans(child)%tag_j_upper = metadata(index + 7)
        collection%plans(child)%coarse_i_lower = metadata(index + 8)
        collection%plans(child)%coarse_i_upper = metadata(index + 9)
        collection%plans(child)%coarse_j_lower = metadata(index + 10)
        collection%plans(child)%coarse_j_upper = metadata(index + 11)
        index = index + 12
      end do
    end if
    local_ok = collection%is_valid()
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    allocate(fine_geometries(patch_count))
    local_ok = .true.
    do child = 1, patch_count
      call geometry_builder( &
        coarse_geometry, collection%plans(child)%coarse_i_lower, &
        collection%plans(child)%coarse_i_upper, &
        collection%plans(child)%coarse_j_lower, &
        collection%plans(child)%coarse_j_upper, refinement_ratio, &
        fine_geometries(child), entity_ok)
      local_ok = local_ok .and. entity_ok
    end do
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call regrid_sparse_owned_reactive_eb_patch_set_2d( &
      species, distribution, sparse_patch_set, coarse_geometry, &
      patch_set_template, fine_geometries, collection, refinement_ratio, &
      local_ok, changed, restriction_transfers, prolongation_transfers, &
      overlap_transfers)
    if (.not. local_ok) then
      changed = .false.
      return
    end if

    ok = .true.
    if (present(local_root_transfers)) local_root_transfers = transfers
    if (present(local_restriction_transfers)) &
      local_restriction_transfers = restriction_transfers
    if (present(local_prolongation_transfers)) &
      local_prolongation_transfers = prolongation_transfers
    if (present(local_overlap_transfers)) &
      local_overlap_transfers = overlap_transfers
  end subroutine regrid_tagged_sparse_owned_reactive_eb_patch_set_2d

  subroutine gather_sparse_root_to_owner_2d( &
      distribution, sparse_patch_set, coarse_geometry, root_state, &
      root_temperature, local_transfers, ok)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(in) :: sparse_patch_set
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    real(dp), allocatable, intent(out) :: root_state(:, :, :)
    real(dp), allocatable, intent(out) :: root_temperature(:, :)
    integer, intent(inout) :: local_transfers
    logical, intent(out) :: ok

    type(MPI_Status) :: status
    real(dp), allocatable :: payload(:)
    logical :: accepted, global_ok, local_ok
    integer :: cell_count, ierr, j_lower, j_upper, owner, root_owner
    integer :: rows, state_count, tile

    ok = .false.
    root_owner = distribution%root_level_owner()
    local_ok = root_owner >= 0 .and. root_owner < distribution%nranks .and. &
      sparse_patch_set%nvar >= 1 .and. &
      allocated(sparse_patch_set%root_tiles) .and. &
      size(sparse_patch_set%root_tiles) == distribution%root_tile_count()
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    if (distribution%rank == root_owner) then
      allocate(root_state( &
        sparse_patch_set%nvar, coarse_geometry%nx, coarse_geometry%ny), &
        source=0.0_dp)
      allocate(root_temperature( &
        coarse_geometry%nx, coarse_geometry%ny), source=1.0_dp)
    end if
    do tile = 1, distribution%root_tile_count()
      owner = distribution%root_tiles(tile)%owner
      j_lower = distribution%root_tiles(tile)%j_lower
      j_upper = distribution%root_tiles(tile)%j_upper
      rows = j_upper - j_lower + 1
      cell_count = distribution%root_tiles(tile)%cell_count
      state_count = sparse_patch_set%nvar * cell_count
      if (owner == root_owner) then
        if (distribution%rank == root_owner) then
          root_state(:, :, j_lower:j_upper) = &
            sparse_patch_set%root_tiles(tile)%state
          root_temperature(:, j_lower:j_upper) = &
            sparse_patch_set%root_tiles(tile)%temperature
        end if
      else if (distribution%rank == owner) then
        allocate(payload(state_count + cell_count))
        payload(1:state_count) = reshape( &
          sparse_patch_set%root_tiles(tile)%state, [state_count])
        payload(state_count + 1:state_count + cell_count) = reshape( &
          sparse_patch_set%root_tiles(tile)%temperature, [cell_count])
        call MPI_Send( &
          payload, size(payload), MPI_DOUBLE_PRECISION, root_owner, &
          sparse_root_gather_tag, distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS) return
        local_transfers = local_transfers + 1
        deallocate(payload)
      else if (distribution%rank == root_owner) then
        allocate(payload(state_count + cell_count))
        call MPI_Recv( &
          payload, size(payload), MPI_DOUBLE_PRECISION, owner, &
          sparse_root_gather_tag, distribution%comm, status, ierr)
        if (ierr /= MPI_SUCCESS) return
        root_state(:, :, j_lower:j_upper) = reshape( &
          payload(1:state_count), &
          [sparse_patch_set%nvar, coarse_geometry%nx, rows])
        root_temperature(:, j_lower:j_upper) = reshape( &
          payload(state_count + 1:state_count + cell_count), &
          [coarse_geometry%nx, rows])
        deallocate(payload)
      end if
    end do
    local_ok = .true.
    if (distribution%rank == root_owner) local_ok = &
      allocated(root_state) .and. allocated(root_temperature)
    if (distribution%rank == root_owner .and. local_ok) local_ok = &
      all(ieee_is_finite(root_state)) .and. &
      all(ieee_is_finite(root_temperature))
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    ok = global_ok .and. accepted
  end subroutine gather_sparse_root_to_owner_2d

  subroutine advance_sparse_owned_reactive_eb_root_tiles_hydro_2d( &
      species, distribution, sparse_patch_set, geometry, solver, &
      reconstruction, limiter, state_redist_target_volume_fraction, &
      state_redist_max_order, dt, local_tile_states, local_tile_fluxes, &
      local_transfers, ok, local_tile_advances, local_computed_cells)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(in) :: sparse_patch_set
    type(eb_geometry_2d), intent(in) :: geometry
    character(len=*), intent(in) :: solver, reconstruction, limiter
    real(dp), intent(in) :: state_redist_target_volume_fraction, dt
    integer, intent(in) :: state_redist_max_order
    type(mpi_amr_eb_root_tile_transport_state_2d), allocatable, intent(out) :: &
      local_tile_states(:)
    type(mpi_amr_eb_root_tile_transport_flux_2d), allocatable, intent(out) :: &
      local_tile_fluxes(:)
    integer, intent(inout) :: local_transfers
    logical, intent(out) :: ok
    integer, intent(out) :: local_tile_advances, local_computed_cells

    type(eb_geometry_2d) :: band_geometry
    type(MPI_Status) :: status
    real(dp), allocatable :: band_state(:, :, :), band_temperature(:, :)
    real(dp), allocatable :: band_updated_state(:, :, :)
    real(dp), allocatable :: band_updated_temperature(:, :)
    real(dp), allocatable :: band_x_flux(:, :, :), band_y_flux(:, :, :)
    real(dp), allocatable :: payload(:)
    logical :: accepted, entity_ok, global_ok, local_ok
    integer :: band_j_lower, band_j_upper, band_rows
    integer :: computed_cells, face_j_lower, face_j_upper
    integer :: halo_cell_count, halo_state_count, halo_value_count
    integer :: ierr, local_face_j_lower, local_face_j_upper
    integer :: local_j_lower, local_j_upper, nvar, overlap_j_lower
    integer :: overlap_j_upper, source, source_owner
    integer :: target, target_owner, tile_advances

    ok = .false.
    local_tile_advances = 0
    local_computed_cells = 0
    tile_advances = 0
    computed_cells = 0
    nvar = reactive_nvar(size(species))
    local_ok = nvar >= 1 .and. sparse_patch_set%nvar == nvar .and. &
      allocated(sparse_patch_set%root_tiles) .and. &
      size(sparse_patch_set%root_tiles) == distribution%root_tile_count()
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    allocate(local_tile_states(distribution%root_tile_count()))
    allocate(local_tile_fluxes(distribution%root_tile_count()))

    do target = 1, distribution%root_tile_count()
      target_owner = distribution%root_tiles(target)%owner
      band_j_lower = max(1, distribution%root_tiles(target)%j_lower - &
        mpi_amr_eb_root_tile_hydro_halo_cells)
      band_j_upper = min(geometry%ny, &
        distribution%root_tiles(target)%j_upper + &
          mpi_amr_eb_root_tile_hydro_halo_cells)
      band_rows = band_j_upper - band_j_lower + 1
      entity_ok = .true.
      if (distribution%rank == target_owner) then
        call extract_eb_geometry_y_band_2d( &
          geometry, band_j_lower, band_j_upper, band_geometry, entity_ok)
        if (entity_ok) then
          allocate(band_state(nvar, geometry%nx, band_rows), source=0.0_dp)
          allocate(band_temperature(geometry%nx, band_rows), source=0.0_dp)
          allocate(band_updated_state, mold=band_state)
          allocate(band_updated_temperature, mold=band_temperature)
          allocate(band_x_flux(nvar, 0:geometry%nx, band_rows))
          allocate(band_y_flux(nvar, geometry%nx, 0:band_rows))
        end if
      end if
      call all_ranks_accept_eb_2d( &
        distribution, entity_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) return

      do source = 1, distribution%root_tile_count()
        overlap_j_lower = max( &
          band_j_lower, distribution%root_tiles(source)%j_lower)
        overlap_j_upper = min( &
          band_j_upper, distribution%root_tiles(source)%j_upper)
        if (overlap_j_upper < overlap_j_lower) cycle
        source_owner = distribution%root_tiles(source)%owner
        local_j_lower = overlap_j_lower - &
          distribution%root_tiles(source)%j_lower + 1
        local_j_upper = overlap_j_upper - &
          distribution%root_tiles(source)%j_lower + 1
        if (source_owner == target_owner) then
          if (distribution%rank == target_owner) then
            band_state(:, :, &
              overlap_j_lower - band_j_lower + 1: &
                overlap_j_upper - band_j_lower + 1) = &
              sparse_patch_set%root_tiles(source)%state( &
                :, :, local_j_lower:local_j_upper)
            band_temperature(:, &
              overlap_j_lower - band_j_lower + 1: &
                overlap_j_upper - band_j_lower + 1) = &
              sparse_patch_set%root_tiles(source)%temperature( &
                :, local_j_lower:local_j_upper)
          end if
          cycle
        end if

        halo_cell_count = geometry%nx * &
          (overlap_j_upper - overlap_j_lower + 1)
        halo_state_count = nvar * halo_cell_count
        halo_value_count = halo_state_count + halo_cell_count
        if (distribution%rank == source_owner) then
          allocate(payload(halo_value_count))
          payload(1:halo_state_count) = reshape( &
            sparse_patch_set%root_tiles(source)%state( &
              :, :, local_j_lower:local_j_upper), [halo_state_count])
          payload(halo_state_count + 1:halo_value_count) = reshape( &
            sparse_patch_set%root_tiles(source)%temperature( &
              :, local_j_lower:local_j_upper), [halo_cell_count])
          call MPI_Send( &
            payload, halo_value_count, MPI_DOUBLE_PRECISION, target_owner, &
            sparse_root_halo_tag, distribution%comm, ierr)
          if (ierr /= MPI_SUCCESS) return
          local_transfers = local_transfers + 1
          deallocate(payload)
        else if (distribution%rank == target_owner) then
          allocate(payload(halo_value_count))
          call MPI_Recv( &
            payload, halo_value_count, MPI_DOUBLE_PRECISION, source_owner, &
            sparse_root_halo_tag, distribution%comm, status, ierr)
          if (ierr /= MPI_SUCCESS) return
          band_state(:, :, &
            overlap_j_lower - band_j_lower + 1: &
              overlap_j_upper - band_j_lower + 1) = reshape( &
            payload(1:halo_state_count), &
            [nvar, geometry%nx, overlap_j_upper - overlap_j_lower + 1])
          band_temperature(:, &
            overlap_j_lower - band_j_lower + 1: &
              overlap_j_upper - band_j_lower + 1) = reshape( &
            payload(halo_state_count + 1:halo_value_count), &
            [geometry%nx, overlap_j_upper - overlap_j_lower + 1])
          deallocate(payload)
        end if
      end do

      entity_ok = .true.
      if (distribution%rank == target_owner) entity_ok = &
        all(ieee_is_finite(band_state)) .and. &
        all(ieee_is_finite(band_temperature))
      call all_ranks_accept_eb_2d( &
        distribution, entity_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) return
      if (distribution%rank == target_owner) then
        call advance_reactive_eb_level_2d( &
          species, band_state, band_temperature, band_geometry, solver, &
          reconstruction, limiter, state_redist_target_volume_fraction, &
          state_redist_max_order, dt, band_updated_state, &
          band_updated_temperature, band_x_flux, band_y_flux, entity_ok)
      end if
      call all_ranks_accept_eb_2d( &
        distribution, entity_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) return

      face_j_lower = distribution%root_tiles(target)%j_lower - 1
      face_j_upper = distribution%root_tiles(target)%j_upper - 1
      if (distribution%root_tiles(target)%j_upper == geometry%ny) &
        face_j_upper = geometry%ny
      local_j_lower = distribution%root_tiles(target)%j_lower - &
        band_j_lower + 1
      local_j_upper = distribution%root_tiles(target)%j_upper - &
        band_j_lower + 1
      local_face_j_lower = face_j_lower - band_j_lower + 1
      local_face_j_upper = face_j_upper - band_j_lower + 1
      if (distribution%rank == target_owner) then
        allocate(local_tile_states(target)%start_state( &
          nvar, 1:geometry%nx, &
          distribution%root_tiles(target)%j_lower: &
            distribution%root_tiles(target)%j_upper))
        allocate(local_tile_states(target)%start_temperature( &
          1:geometry%nx, distribution%root_tiles(target)%j_lower: &
            distribution%root_tiles(target)%j_upper))
        allocate(local_tile_states(target)%end_state( &
          nvar, 1:geometry%nx, &
          distribution%root_tiles(target)%j_lower: &
            distribution%root_tiles(target)%j_upper))
        allocate(local_tile_states(target)%end_temperature( &
          1:geometry%nx, distribution%root_tiles(target)%j_lower: &
            distribution%root_tiles(target)%j_upper))
        allocate(local_tile_states(target)%corrected_state( &
          nvar, 1:geometry%nx, &
          distribution%root_tiles(target)%j_lower: &
            distribution%root_tiles(target)%j_upper))
        allocate(local_tile_states(target)%corrected_temperature( &
          1:geometry%nx, distribution%root_tiles(target)%j_lower: &
            distribution%root_tiles(target)%j_upper))
        local_tile_states(target)%start_state = &
          band_state(:, :, local_j_lower:local_j_upper)
        local_tile_states(target)%start_temperature = &
          band_temperature(:, local_j_lower:local_j_upper)
        local_tile_states(target)%end_state = &
          band_updated_state(:, :, local_j_lower:local_j_upper)
        local_tile_states(target)%end_temperature = &
          band_updated_temperature(:, local_j_lower:local_j_upper)
        local_tile_states(target)%corrected_state = &
          local_tile_states(target)%end_state
        local_tile_states(target)%corrected_temperature = &
          local_tile_states(target)%end_temperature
        allocate(local_tile_fluxes(target)%x_flux( &
          nvar, 0:geometry%nx, &
          distribution%root_tiles(target)%j_lower: &
            distribution%root_tiles(target)%j_upper))
        allocate(local_tile_fluxes(target)%y_flux( &
          nvar, 1:geometry%nx, face_j_lower:face_j_upper))
        local_tile_fluxes(target)%x_flux = &
          band_x_flux(:, :, local_j_lower:local_j_upper)
        local_tile_fluxes(target)%y_flux = &
          band_y_flux(:, :, local_face_j_lower:local_face_j_upper)
      end if

      if (distribution%rank == target_owner) then
        tile_advances = tile_advances + 1
        computed_cells = computed_cells + geometry%nx * band_rows
        deallocate( &
          band_state, band_temperature, band_updated_state, &
          band_updated_temperature, band_x_flux, band_y_flux)
      end if
    end do

    local_ok = .true.
    do target = 1, distribution%root_tile_count()
      if (.not. distribution%root_tile_is_local(target)) cycle
      local_ok = local_ok .and. &
        allocated(local_tile_states(target)%start_state) .and. &
        allocated(local_tile_states(target)%start_temperature) .and. &
        allocated(local_tile_states(target)%end_state) .and. &
        allocated(local_tile_states(target)%end_temperature) .and. &
        allocated(local_tile_states(target)%corrected_state) .and. &
        allocated(local_tile_states(target)%corrected_temperature) .and. &
        allocated(local_tile_fluxes(target)%x_flux) .and. &
        allocated(local_tile_fluxes(target)%y_flux)
      if (.not. local_ok) cycle
      local_ok = all(ieee_is_finite( &
          local_tile_states(target)%start_state)) .and. &
        all(ieee_is_finite( &
          local_tile_states(target)%start_temperature)) .and. &
        all(ieee_is_finite(local_tile_states(target)%end_state)) .and. &
        all(ieee_is_finite( &
          local_tile_states(target)%end_temperature)) .and. &
        all(ieee_is_finite( &
          local_tile_states(target)%corrected_state)) .and. &
        all(ieee_is_finite( &
          local_tile_states(target)%corrected_temperature)) .and. &
        all(ieee_is_finite( &
          local_tile_fluxes(target)%x_flux)) .and. &
        all(ieee_is_finite(local_tile_fluxes(target)%y_flux))
    end do
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    ok = .true.
    local_tile_advances = tile_advances
    local_computed_cells = computed_cells
  end subroutine advance_sparse_owned_reactive_eb_root_tiles_hydro_2d

  subroutine advance_sparse_owned_reactive_eb_root_tiles_transport_euler_2d( &
      species, transport, distribution, sparse_patch_set, geometry, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      state_redist_target_volume_fraction, state_redist_max_order, &
      local_tile_states, local_tile_fluxes, minimum_theta, &
      local_transfers, ok, local_tile_advances, local_computed_cells)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(in) :: sparse_patch_set
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), intent(in) :: dt, state_redist_target_volume_fraction
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    integer, intent(in) :: state_redist_max_order
    type(mpi_amr_eb_root_tile_transport_state_2d), allocatable, intent(out) :: &
      local_tile_states(:)
    type(mpi_amr_eb_root_tile_transport_flux_2d), allocatable, intent(out) :: &
      local_tile_fluxes(:)
    real(dp), intent(out) :: minimum_theta
    integer, intent(inout) :: local_transfers
    logical, intent(out) :: ok
    integer, intent(out) :: local_tile_advances, local_computed_cells

    type(eb_geometry_2d) :: band_geometry
    type(MPI_Status) :: status
    real(dp), allocatable :: band_rhs(:, :, :), band_state(:, :, :)
    real(dp), allocatable :: band_temperature(:, :)
    real(dp), allocatable :: band_updated_state(:, :, :)
    real(dp), allocatable :: band_updated_temperature(:, :)
    real(dp), allocatable :: band_x_flux(:, :, :), band_y_flux(:, :, :)
    real(dp), allocatable :: payload(:)
    integer, allocatable :: band_source_rows(:)
    real(dp) :: band_theta, theta
    logical :: accepted, entity_ok, global_ok, local_ok, wrapped_band
    integer :: band_j_lower, band_j_upper, band_rows
    integer :: computed_cells, face_j_lower, face_j_upper
    integer :: halo_cell_count, halo_state_count, halo_value_count
    integer :: ierr, local_face_j_lower, local_face_j_upper
    integer :: local_j_lower, local_j_upper, nvar, offset, overlap_j_lower
    integer :: overlap_j_upper, rows, source, source_owner
    integer :: segment_j_lower, segment_j_upper, source_global_j_lower
    integer :: source_global_j_upper
    integer :: target, target_owner, tile_advances

    ok = .false.
    minimum_theta = 1.0_dp
    local_tile_advances = 0
    local_computed_cells = 0
    tile_advances = 0
    computed_cells = 0
    theta = 1.0_dp
    nvar = reactive_nvar(size(species))
    local_ok = nvar >= 1 .and. sparse_patch_set%nvar == nvar .and. &
      allocated(sparse_patch_set%root_tiles) .and. &
      size(sparse_patch_set%root_tiles) == distribution%root_tile_count()
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    allocate(local_tile_states(distribution%root_tile_count()))
    allocate(local_tile_fluxes(distribution%root_tile_count()))

    do target = 1, distribution%root_tile_count()
      target_owner = distribution%root_tiles(target)%owner
      rows = distribution%root_tiles(target)%j_upper - &
        distribution%root_tiles(target)%j_lower + 1
      band_j_lower = max(1, distribution%root_tiles(target)%j_lower - &
        mpi_amr_eb_root_tile_transport_halo_cells)
      band_j_upper = min(geometry%ny, &
        distribution%root_tiles(target)%j_upper + &
          mpi_amr_eb_root_tile_transport_halo_cells)
      wrapped_band = reactive_boundary_is_periodic( &
          boundaries%face(boundary_y_lower)) .and. &
        (distribution%root_tiles(target)%j_lower == 1 .or. &
         distribution%root_tiles(target)%j_upper == geometry%ny) .and. &
        rows + 2 * (mpi_amr_eb_root_tile_transport_halo_cells + 1) < &
          geometry%ny
      if (reactive_boundary_is_periodic( &
          boundaries%face(boundary_y_lower)) .and. &
          (distribution%root_tiles(target)%j_lower == 1 .or. &
           distribution%root_tiles(target)%j_upper == geometry%ny) .and. &
          .not. wrapped_band) then
        band_j_lower = 1
        band_j_upper = geometry%ny
      end if
      if (wrapped_band) then
        ! Keep one row beyond the qualified halo so the artificial gap cannot
        ! contaminate a cell in the target's dependency footprint.
        band_rows = rows + &
          2 * (mpi_amr_eb_root_tile_transport_halo_cells + 1)
      else
        band_rows = band_j_upper - band_j_lower + 1
      end if
      allocate(band_source_rows(band_rows))
      if (wrapped_band .and. &
          distribution%root_tiles(target)%j_lower == 1) then
        band_j_upper = distribution%root_tiles(target)%j_upper + &
          mpi_amr_eb_root_tile_transport_halo_cells + 1
        band_source_rows(1:band_j_upper) = &
          [(source, source = 1, band_j_upper)]
        band_source_rows(band_j_upper + 1:band_rows) = [( &
          source, source = geometry%ny - &
            mpi_amr_eb_root_tile_transport_halo_cells, geometry%ny)]
        local_j_lower = 1
        local_j_upper = rows
      else if (wrapped_band) then
        band_j_lower = distribution%root_tiles(target)%j_lower - &
          mpi_amr_eb_root_tile_transport_halo_cells - 1
        band_source_rows( &
          1:mpi_amr_eb_root_tile_transport_halo_cells + 1) = [( &
            source, source = 1, &
              mpi_amr_eb_root_tile_transport_halo_cells + 1)]
        band_source_rows( &
          mpi_amr_eb_root_tile_transport_halo_cells + 2:band_rows) = [( &
            source, source = band_j_lower, geometry%ny)]
        local_j_lower = band_rows - rows + 1
        local_j_upper = band_rows
      else
        band_source_rows = [( &
          source, source = band_j_lower, band_j_upper)]
        local_j_lower = distribution%root_tiles(target)%j_lower - &
          band_j_lower + 1
        local_j_upper = distribution%root_tiles(target)%j_upper - &
          band_j_lower + 1
      end if
      entity_ok = .true.
      if (distribution%rank == target_owner) then
        if (wrapped_band) then
          call extract_eb_geometry_y_rows_2d( &
            geometry, band_source_rows, band_geometry, entity_ok)
        else
          call extract_eb_geometry_y_band_2d( &
            geometry, band_j_lower, band_j_upper, band_geometry, entity_ok)
        end if
        if (entity_ok) then
          allocate(band_state(nvar, geometry%nx, band_rows), source=0.0_dp)
          allocate(band_temperature(geometry%nx, band_rows), source=0.0_dp)
          allocate(band_rhs, mold=band_state)
          allocate(band_updated_state, mold=band_state)
          allocate(band_updated_temperature, mold=band_temperature)
          allocate(band_x_flux(nvar, 0:geometry%nx, band_rows))
          allocate(band_y_flux(nvar, geometry%nx, 0:band_rows))
        end if
      end if
      call all_ranks_accept_eb_2d( &
        distribution, entity_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) return

      ! A cyclic edge band contains two globally contiguous row segments.
      segment_j_lower = 1
      do while (segment_j_lower <= band_rows)
        source = 0
        do offset = 1, distribution%root_tile_count()
          if (band_source_rows(segment_j_lower) < &
              distribution%root_tiles(offset)%j_lower .or. &
              band_source_rows(segment_j_lower) > &
              distribution%root_tiles(offset)%j_upper) cycle
          source = offset
          exit
        end do
        if (source == 0) return
        segment_j_upper = segment_j_lower
        do while (segment_j_upper < band_rows)
          if (band_source_rows(segment_j_upper + 1) /= &
              band_source_rows(segment_j_upper) + 1 .or. &
              band_source_rows(segment_j_upper + 1) > &
                distribution%root_tiles(source)%j_upper) exit
          segment_j_upper = segment_j_upper + 1
        end do
        source_global_j_lower = band_source_rows(segment_j_lower)
        source_global_j_upper = band_source_rows(segment_j_upper)
        source_owner = distribution%root_tiles(source)%owner
        overlap_j_lower = source_global_j_lower - &
          distribution%root_tiles(source)%j_lower + 1
        overlap_j_upper = source_global_j_upper - &
          distribution%root_tiles(source)%j_lower + 1
        if (source_owner == target_owner) then
          if (distribution%rank == target_owner) then
            band_state(:, :, segment_j_lower:segment_j_upper) = &
              sparse_patch_set%root_tiles(source)%state( &
                :, :, overlap_j_lower:overlap_j_upper)
            band_temperature(:, segment_j_lower:segment_j_upper) = &
              sparse_patch_set%root_tiles(source)%temperature( &
                :, overlap_j_lower:overlap_j_upper)
          end if
        else
          halo_cell_count = geometry%nx * &
            (segment_j_upper - segment_j_lower + 1)
          halo_state_count = nvar * halo_cell_count
          halo_value_count = halo_state_count + halo_cell_count
          if (distribution%rank == source_owner) then
            allocate(payload(halo_value_count))
            payload(1:halo_state_count) = reshape( &
              sparse_patch_set%root_tiles(source)%state( &
                :, :, overlap_j_lower:overlap_j_upper), [halo_state_count])
            payload(halo_state_count + 1:halo_value_count) = reshape( &
              sparse_patch_set%root_tiles(source)%temperature( &
                :, overlap_j_lower:overlap_j_upper), [halo_cell_count])
            call MPI_Send( &
              payload, halo_value_count, MPI_DOUBLE_PRECISION, &
              target_owner, sparse_root_halo_tag, distribution%comm, ierr)
            if (ierr /= MPI_SUCCESS) return
            local_transfers = local_transfers + 1
            deallocate(payload)
          else if (distribution%rank == target_owner) then
            allocate(payload(halo_value_count))
            call MPI_Recv( &
              payload, halo_value_count, MPI_DOUBLE_PRECISION, &
              source_owner, sparse_root_halo_tag, distribution%comm, &
              status, ierr)
            if (ierr /= MPI_SUCCESS) return
            band_state(:, :, segment_j_lower:segment_j_upper) = reshape( &
              payload(1:halo_state_count), &
              [nvar, geometry%nx, &
                segment_j_upper - segment_j_lower + 1])
            band_temperature(:, segment_j_lower:segment_j_upper) = reshape( &
              payload(halo_state_count + 1:halo_value_count), &
              [geometry%nx, segment_j_upper - segment_j_lower + 1])
            deallocate(payload)
          end if
        end if
        segment_j_lower = segment_j_upper + 1
      end do

      entity_ok = .true.
      if (distribution%rank == target_owner) entity_ok = &
        all(ieee_is_finite(band_state)) .and. &
        all(ieee_is_finite(band_temperature))
      call all_ranks_accept_eb_2d( &
        distribution, entity_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) return
      if (distribution%rank == target_owner) then
        call reactive_eb_transport_fluxes_rhs_2d( &
          species, transport, band_state, band_temperature, band_geometry, &
          dt, viscosity_enabled, thermal_conduction_enabled, &
          species_diffusion_enabled, barodiffusion_enabled, boundaries, &
          band_rhs, band_x_flux, band_y_flux, band_theta, entity_ok)
        if (entity_ok) call advance_reactive_eb_state_redistributed_2d( &
          species, band_state, band_temperature, band_geometry, band_rhs, dt, &
          band_updated_state, band_updated_temperature, entity_ok, &
          state_redist_target_volume_fraction, state_redist_max_order)
      end if
      call all_ranks_accept_eb_2d( &
        distribution, entity_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) return

      rows = distribution%root_tiles(target)%j_upper - &
        distribution%root_tiles(target)%j_lower + 1
      face_j_lower = distribution%root_tiles(target)%j_lower - 1
      face_j_upper = distribution%root_tiles(target)%j_upper - 1
      if (distribution%root_tiles(target)%j_upper == geometry%ny) &
        face_j_upper = geometry%ny
      local_face_j_lower = local_j_lower - 1
      local_face_j_upper = local_j_upper - 1
      if (distribution%root_tiles(target)%j_upper == geometry%ny) &
        local_face_j_upper = local_j_upper
      if (distribution%rank == target_owner) then
        allocate(local_tile_states(target)%start_state( &
          nvar, 1:geometry%nx, &
          distribution%root_tiles(target)%j_lower: &
            distribution%root_tiles(target)%j_upper))
        allocate(local_tile_states(target)%start_temperature( &
          1:geometry%nx, distribution%root_tiles(target)%j_lower: &
            distribution%root_tiles(target)%j_upper))
        allocate(local_tile_states(target)%corrected_state( &
          nvar, 1:geometry%nx, &
          distribution%root_tiles(target)%j_lower: &
            distribution%root_tiles(target)%j_upper))
        allocate(local_tile_states(target)%end_state( &
          nvar, 1:geometry%nx, &
          distribution%root_tiles(target)%j_lower: &
            distribution%root_tiles(target)%j_upper))
        allocate(local_tile_states(target)%end_temperature( &
          1:geometry%nx, distribution%root_tiles(target)%j_lower: &
            distribution%root_tiles(target)%j_upper))
        allocate(local_tile_states(target)%corrected_temperature( &
          1:geometry%nx, distribution%root_tiles(target)%j_lower: &
            distribution%root_tiles(target)%j_upper))
        local_tile_states(target)%start_state = &
          band_state(:, :, local_j_lower:local_j_upper)
        local_tile_states(target)%start_temperature = &
          band_temperature(:, local_j_lower:local_j_upper)
        local_tile_states(target)%end_state = &
          band_updated_state(:, :, local_j_lower:local_j_upper)
        local_tile_states(target)%end_temperature = &
          band_updated_temperature(:, local_j_lower:local_j_upper)
        local_tile_states(target)%corrected_state = &
          local_tile_states(target)%end_state
        local_tile_states(target)%corrected_temperature = &
          local_tile_states(target)%end_temperature
        allocate(local_tile_fluxes(target)%x_flux( &
          nvar, 0:geometry%nx, &
          distribution%root_tiles(target)%j_lower: &
            distribution%root_tiles(target)%j_upper))
        allocate(local_tile_fluxes(target)%y_flux( &
          nvar, 1:geometry%nx, face_j_lower:face_j_upper))
        local_tile_fluxes(target)%x_flux = &
          band_x_flux(:, :, local_j_lower:local_j_upper)
        local_tile_fluxes(target)%y_flux = &
          band_y_flux(:, :, local_face_j_lower:local_face_j_upper)
      end if

      if (distribution%rank == target_owner) then
        tile_advances = tile_advances + 1
        computed_cells = computed_cells + geometry%nx * band_rows
        theta = min(theta, band_theta)
        deallocate( &
          band_state, band_temperature, band_rhs, band_updated_state, &
          band_updated_temperature, band_x_flux, band_y_flux)
      end if
      deallocate(band_source_rows)
    end do

    local_ok = ieee_is_finite(theta) .and. theta >= 0.0_dp .and. &
      theta <= 1.0_dp
    do target = 1, distribution%root_tile_count()
      if (.not. distribution%root_tile_is_local(target)) cycle
      local_ok = local_ok .and. &
        allocated(local_tile_states(target)%start_state) .and. &
        allocated(local_tile_states(target)%start_temperature) .and. &
        allocated(local_tile_states(target)%end_state) .and. &
        allocated(local_tile_states(target)%end_temperature) .and. &
        allocated(local_tile_states(target)%corrected_state) .and. &
        allocated(local_tile_states(target)%corrected_temperature) .and. &
        allocated(local_tile_fluxes(target)%x_flux) .and. &
        allocated(local_tile_fluxes(target)%y_flux)
      if (.not. local_ok) cycle
      local_ok = all(ieee_is_finite( &
          local_tile_states(target)%start_state)) .and. &
        all(ieee_is_finite( &
          local_tile_states(target)%start_temperature)) .and. &
        all(ieee_is_finite(local_tile_states(target)%end_state)) .and. &
        all(ieee_is_finite( &
          local_tile_states(target)%end_temperature)) .and. &
        all(ieee_is_finite( &
          local_tile_states(target)%corrected_state)) .and. &
        all(ieee_is_finite( &
          local_tile_states(target)%corrected_temperature)) .and. &
        all(ieee_is_finite(local_tile_fluxes(target)%x_flux)) .and. &
        all(ieee_is_finite(local_tile_fluxes(target)%y_flux))
    end do
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    ok = .true.
    minimum_theta = theta
    local_tile_advances = tile_advances
    local_computed_cells = computed_cells
  end subroutine advance_sparse_owned_reactive_eb_root_tiles_transport_euler_2d

  subroutine reduce_sparse_root_transport_boundary_change_2d( &
      distribution, geometry, component_count, local_tile_fluxes, dt, &
      boundary_change, ok)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(eb_geometry_2d), intent(in) :: geometry
    integer, intent(in) :: component_count
    type(mpi_amr_eb_root_tile_transport_flux_2d), intent(in) :: &
      local_tile_fluxes(:)
    real(dp), intent(in) :: dt
    real(dp), allocatable, intent(out) :: boundary_change(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: local_change(:)
    logical :: accepted, global_ok, local_ok
    integer :: face_j_lower, face_j_upper, i, ierr, j
    integer :: j_lower, j_upper, target

    ok = .false.
    local_ok = component_count >= 1 .and. ieee_is_finite(dt) .and. &
      dt > 0.0_dp .and. &
      size(local_tile_fluxes) == distribution%root_tile_count()
    allocate(local_change(component_count), source=0.0_dp)
    if (local_ok) then
      do target = 1, distribution%root_tile_count()
        if (.not. distribution%root_tile_is_local(target)) cycle
        j_lower = distribution%root_tiles(target)%j_lower
        j_upper = distribution%root_tiles(target)%j_upper
        face_j_lower = j_lower - 1
        face_j_upper = j_upper - 1
        if (j_upper == geometry%ny) face_j_upper = geometry%ny
        local_ok = allocated(local_tile_fluxes(target)%x_flux) .and. &
          allocated(local_tile_fluxes(target)%y_flux)
        if (local_ok) local_ok = &
          size(local_tile_fluxes(target)%x_flux, 1) == component_count .and. &
          lbound(local_tile_fluxes(target)%x_flux, 2) == 0 .and. &
          ubound(local_tile_fluxes(target)%x_flux, 2) == geometry%nx .and. &
          lbound(local_tile_fluxes(target)%x_flux, 3) == j_lower .and. &
          ubound(local_tile_fluxes(target)%x_flux, 3) == j_upper .and. &
          size(local_tile_fluxes(target)%y_flux, 1) == component_count .and. &
          lbound(local_tile_fluxes(target)%y_flux, 2) == 1 .and. &
          ubound(local_tile_fluxes(target)%y_flux, 2) == geometry%nx .and. &
          lbound(local_tile_fluxes(target)%y_flux, 3) == face_j_lower .and. &
          ubound(local_tile_fluxes(target)%y_flux, 3) == face_j_upper
        if (local_ok) local_ok = &
          all(ieee_is_finite(local_tile_fluxes(target)%x_flux)) .and. &
          all(ieee_is_finite(local_tile_fluxes(target)%y_flux))
        if (.not. local_ok) exit
        do j = j_lower, j_upper
          local_change = local_change + dt * geometry%dy * ( &
            geometry%x_face_fraction(0, j) * &
              local_tile_fluxes(target)%x_flux(:, 0, j) - &
            geometry%x_face_fraction(geometry%nx, j) * &
              local_tile_fluxes(target)%x_flux(:, geometry%nx, j))
        end do
        if (j_lower == 1) then
          do i = 1, geometry%nx
            local_change = local_change + dt * geometry%dx * &
              geometry%y_face_fraction(i, 0) * &
                local_tile_fluxes(target)%y_flux(:, i, 0)
          end do
        end if
        if (j_upper == geometry%ny) then
          do i = 1, geometry%nx
            local_change = local_change - dt * geometry%dx * &
              geometry%y_face_fraction(i, geometry%ny) * &
                local_tile_fluxes(target)%y_flux(:, i, geometry%ny)
          end do
        end if
      end do
    end if
    local_ok = local_ok .and. all(ieee_is_finite(local_change))
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    allocate(boundary_change(component_count))
    call MPI_Allreduce( &
      local_change, boundary_change, component_count, MPI_DOUBLE_PRECISION, &
      MPI_SUM, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    local_ok = all(ieee_is_finite(boundary_change))
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    ok = .true.
  end subroutine reduce_sparse_root_transport_boundary_change_2d

  subroutine transfer_sparse_child_coarse_flux_support_2d( &
      distribution, coarse_geometry, fine_geometry, patch, component_count, &
      destination, local_tile_fluxes, x_flux, y_flux, local_transfers, ok)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    integer, intent(in) :: component_count, destination
    type(mpi_amr_eb_root_tile_transport_flux_2d), intent(in) :: &
      local_tile_fluxes(:)
    real(dp), allocatable, intent(out) :: x_flux(:, :, :)
    real(dp), allocatable, intent(out) :: y_flux(:, :, :)
    integer, intent(inout) :: local_transfers
    logical, intent(out) :: ok

    type(MPI_Status) :: status
    real(dp), allocatable :: payload(:)
    logical, allocatable :: x_covered(:, :), y_covered(:, :)
    logical :: accepted, global_ok, local_ok
    integer :: face_j_lower, face_j_upper, ierr, offset
    integer :: source_owner, target
    integer :: value_count, x_count, x_i_lower, x_i_upper
    integer :: x_j_lower, x_j_upper, x_rows
    integer :: y_count, y_i_lower, y_i_upper
    integer :: y_j_lower, y_j_upper, y_rows

    ok = .false.
    local_ok = component_count >= 1 .and. destination >= 0 .and. &
      destination < distribution%nranks .and. &
      patch%is_valid(coarse_geometry, fine_geometry) .and. &
      size(local_tile_fluxes) == distribution%root_tile_count()
    do target = 1, distribution%root_tile_count()
      if (.not. local_ok) exit
      source_owner = distribution%root_tiles(target)%owner
      if (distribution%rank /= source_owner) cycle
      face_j_lower = distribution%root_tiles(target)%j_lower - 1
      face_j_upper = distribution%root_tiles(target)%j_upper - 1
      if (distribution%root_tiles(target)%j_upper == coarse_geometry%ny) &
        face_j_upper = coarse_geometry%ny
      local_ok = allocated(local_tile_fluxes(target)%x_flux) .and. &
        allocated(local_tile_fluxes(target)%y_flux)
      if (.not. local_ok) exit
      local_ok = all(shape(local_tile_fluxes(target)%x_flux) == [ &
          component_count, coarse_geometry%nx + 1, &
          distribution%root_tiles(target)%j_upper - &
            distribution%root_tiles(target)%j_lower + 1]) .and. &
        all(shape(local_tile_fluxes(target)%y_flux) == [ &
          component_count, coarse_geometry%nx, &
          face_j_upper - face_j_lower + 1]) .and. &
        lbound(local_tile_fluxes(target)%x_flux, 2) == 0 .and. &
        lbound(local_tile_fluxes(target)%x_flux, 3) == &
          distribution%root_tiles(target)%j_lower .and. &
        lbound(local_tile_fluxes(target)%y_flux, 2) == 1 .and. &
        lbound(local_tile_fluxes(target)%y_flux, 3) == face_j_lower .and. &
        all(ieee_is_finite(local_tile_fluxes(target)%x_flux)) .and. &
        all(ieee_is_finite(local_tile_fluxes(target)%y_flux))
    end do
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    x_i_lower = patch%coarse_i_lower - 1
    x_i_upper = patch%coarse_i_upper
    y_i_lower = patch%coarse_i_lower
    y_i_upper = patch%coarse_i_upper
    if (distribution%rank == destination) then
      allocate(x_flux( &
        component_count, x_i_lower:x_i_upper, &
        patch%coarse_j_lower:patch%coarse_j_upper), source=0.0_dp)
      allocate(y_flux( &
        component_count, y_i_lower:y_i_upper, &
        patch%coarse_j_lower - 1:patch%coarse_j_upper), source=0.0_dp)
      allocate(x_covered( &
        x_i_lower:x_i_upper, &
        patch%coarse_j_lower:patch%coarse_j_upper), source=.false.)
      allocate(y_covered( &
        y_i_lower:y_i_upper, &
        patch%coarse_j_lower - 1:patch%coarse_j_upper), source=.false.)
    end if

    do target = 1, distribution%root_tile_count()
      source_owner = distribution%root_tiles(target)%owner
      x_j_lower = max( &
        patch%coarse_j_lower, distribution%root_tiles(target)%j_lower)
      x_j_upper = min( &
        patch%coarse_j_upper, distribution%root_tiles(target)%j_upper)
      face_j_lower = distribution%root_tiles(target)%j_lower - 1
      face_j_upper = distribution%root_tiles(target)%j_upper - 1
      if (distribution%root_tiles(target)%j_upper == coarse_geometry%ny) &
        face_j_upper = coarse_geometry%ny
      y_j_lower = max(patch%coarse_j_lower - 1, face_j_lower)
      y_j_upper = min(patch%coarse_j_upper, face_j_upper)
      x_rows = max(0, x_j_upper - x_j_lower + 1)
      y_rows = max(0, y_j_upper - y_j_lower + 1)
      if (x_rows == 0 .and. y_rows == 0) cycle
      x_count = component_count * (x_i_upper - x_i_lower + 1) * x_rows
      y_count = component_count * (y_i_upper - y_i_lower + 1) * y_rows
      value_count = x_count + y_count
      if (source_owner == destination) then
        if (distribution%rank == destination) then
          if (x_rows > 0) then
            x_flux(:, :, x_j_lower:x_j_upper) = &
              local_tile_fluxes(target)%x_flux( &
                :, x_i_lower:x_i_upper, x_j_lower:x_j_upper)
            x_covered(:, x_j_lower:x_j_upper) = .true.
          end if
          if (y_rows > 0) then
            y_flux(:, :, y_j_lower:y_j_upper) = &
              local_tile_fluxes(target)%y_flux( &
                :, y_i_lower:y_i_upper, y_j_lower:y_j_upper)
            y_covered(:, y_j_lower:y_j_upper) = .true.
          end if
        end if
      else if (distribution%rank == source_owner) then
        allocate(payload(value_count))
        offset = 0
        if (x_rows > 0) then
          payload(1:x_count) = reshape( &
            local_tile_fluxes(target)%x_flux( &
              :, x_i_lower:x_i_upper, x_j_lower:x_j_upper), [x_count])
          offset = x_count
        end if
        if (y_rows > 0) payload(offset + 1:value_count) = reshape( &
          local_tile_fluxes(target)%y_flux( &
            :, y_i_lower:y_i_upper, y_j_lower:y_j_upper), [y_count])
        call MPI_Send( &
          payload, value_count, MPI_DOUBLE_PRECISION, destination, &
          sparse_child_coarse_flux_support_tag, distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS) return
        local_transfers = local_transfers + 1
        deallocate(payload)
      else if (distribution%rank == destination) then
        allocate(payload(value_count))
        call MPI_Recv( &
          payload, value_count, MPI_DOUBLE_PRECISION, source_owner, &
          sparse_child_coarse_flux_support_tag, distribution%comm, status, &
          ierr)
        if (ierr /= MPI_SUCCESS) return
        offset = 0
        if (x_rows > 0) then
          x_flux(:, :, x_j_lower:x_j_upper) = reshape( &
            payload(1:x_count), &
            [component_count, x_i_upper - x_i_lower + 1, x_rows])
          x_covered(:, x_j_lower:x_j_upper) = .true.
          offset = x_count
        end if
        if (y_rows > 0) then
          y_flux(:, :, y_j_lower:y_j_upper) = reshape( &
            payload(offset + 1:value_count), &
            [component_count, y_i_upper - y_i_lower + 1, y_rows])
          y_covered(:, y_j_lower:y_j_upper) = .true.
        end if
        deallocate(payload)
      end if
    end do

    local_ok = .true.
    if (distribution%rank == destination) local_ok = &
      all(x_covered) .and. all(y_covered) .and. &
      all(ieee_is_finite(x_flux)) .and. all(ieee_is_finite(y_flux))
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    ok = global_ok .and. accepted
  end subroutine transfer_sparse_child_coarse_flux_support_2d

  subroutine transfer_sparse_child_state_support_2d( &
      distribution, coarse_geometry, fine_geometry, patch, component_count, &
      destination, local_tile_states, context, support_state, &
      support_temperature, local_transfers, ok)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    integer, intent(in) :: component_count, destination
    type(mpi_amr_eb_root_tile_transport_state_2d), intent(in) :: &
      local_tile_states(:)
    type(reactive_eb_patch_exterior_context_2d), intent(out) :: context
    real(dp), allocatable, intent(out) :: support_state(:, :, :)
    real(dp), allocatable, intent(out) :: support_temperature(:, :)
    integer, intent(inout) :: local_transfers
    logical, intent(out) :: ok

    type(MPI_Status) :: status
    real(dp), allocatable :: payload(:)
    real(dp), allocatable :: start_state(:, :, :)
    real(dp), allocatable :: start_temperature(:, :)
    real(dp), allocatable :: end_state(:, :, :)
    real(dp), allocatable :: end_temperature(:, :)
    logical, allocatable :: row_coverage(:)
    logical :: accepted, global_ok, local_ok
    integer :: cell_count, ierr, offset, overlap_j_lower, overlap_j_upper
    integer :: rows, source, state_count, support_i_lower, support_i_upper
    integer :: support_j_lower, support_j_upper, target, value_count

    ok = .false.
    local_ok = component_count >= 1 .and. destination >= 0 .and. &
      destination < distribution%nranks .and. &
      patch%is_valid(coarse_geometry, fine_geometry) .and. &
      size(local_tile_states) == distribution%root_tile_count()
    do target = 1, distribution%root_tile_count()
      if (distribution%rank /= distribution%root_tiles(target)%owner) cycle
      rows = distribution%root_tiles(target)%j_upper - &
        distribution%root_tiles(target)%j_lower + 1
      local_ok = local_ok .and. &
        allocated(local_tile_states(target)%start_state) .and. &
        allocated(local_tile_states(target)%start_temperature) .and. &
        allocated(local_tile_states(target)%end_state) .and. &
        allocated(local_tile_states(target)%end_temperature) .and. &
        allocated(local_tile_states(target)%corrected_state) .and. &
        allocated(local_tile_states(target)%corrected_temperature)
      if (.not. local_ok) cycle
      local_ok = &
        all(shape(local_tile_states(target)%start_state) == &
          [component_count, coarse_geometry%nx, rows]) .and. &
        all(shape(local_tile_states(target)%start_temperature) == &
          [coarse_geometry%nx, rows]) .and. &
        all(shape(local_tile_states(target)%end_state) == &
          shape(local_tile_states(target)%start_state)) .and. &
        all(shape(local_tile_states(target)%end_temperature) == &
          shape(local_tile_states(target)%start_temperature)) .and. &
        all(shape(local_tile_states(target)%corrected_state) == &
          shape(local_tile_states(target)%start_state)) .and. &
        all(shape(local_tile_states(target)%corrected_temperature) == &
          shape(local_tile_states(target)%start_temperature)) .and. &
        lbound(local_tile_states(target)%start_state, 2) == 1 .and. &
        lbound(local_tile_states(target)%start_state, 3) == &
          distribution%root_tiles(target)%j_lower .and. &
        lbound(local_tile_states(target)%start_temperature, 1) == 1 .and. &
        lbound(local_tile_states(target)%start_temperature, 2) == &
          distribution%root_tiles(target)%j_lower .and. &
        lbound(local_tile_states(target)%end_state, 2) == 1 .and. &
        lbound(local_tile_states(target)%end_state, 3) == &
          distribution%root_tiles(target)%j_lower .and. &
        lbound(local_tile_states(target)%end_temperature, 1) == 1 .and. &
        lbound(local_tile_states(target)%end_temperature, 2) == &
          distribution%root_tiles(target)%j_lower .and. &
        lbound(local_tile_states(target)%corrected_state, 2) == 1 .and. &
        lbound(local_tile_states(target)%corrected_state, 3) == &
          distribution%root_tiles(target)%j_lower .and. &
        lbound(local_tile_states(target)%corrected_temperature, 1) == 1 .and. &
        lbound(local_tile_states(target)%corrected_temperature, 2) == &
          distribution%root_tiles(target)%j_lower .and. &
        all(ieee_is_finite(local_tile_states(target)%start_state)) .and. &
        all(ieee_is_finite(local_tile_states(target)%start_temperature)) .and. &
        all(ieee_is_finite(local_tile_states(target)%end_state)) .and. &
        all(ieee_is_finite(local_tile_states(target)%end_temperature)) .and. &
        all(ieee_is_finite(local_tile_states(target)%corrected_state)) .and. &
        all(ieee_is_finite( &
          local_tile_states(target)%corrected_temperature))
    end do
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    support_i_lower = max(1, patch%coarse_i_lower - 2)
    support_i_upper = min(coarse_geometry%nx, patch%coarse_i_upper + 2)
    support_j_lower = max(1, patch%coarse_j_lower - 2)
    support_j_upper = min(coarse_geometry%ny, patch%coarse_j_upper + 2)
    if (distribution%rank == destination) then
      allocate(start_state( &
        component_count, support_i_lower:support_i_upper, &
        support_j_lower:support_j_upper))
      allocate(start_temperature( &
        support_i_lower:support_i_upper, support_j_lower:support_j_upper))
      allocate(end_state( &
        component_count, support_i_lower:support_i_upper, &
        support_j_lower:support_j_upper))
      allocate(end_temperature( &
        support_i_lower:support_i_upper, support_j_lower:support_j_upper))
      allocate(support_state( &
        component_count, support_i_lower:support_i_upper, &
        support_j_lower:support_j_upper))
      allocate(support_temperature( &
        support_i_lower:support_i_upper, support_j_lower:support_j_upper))
      allocate(row_coverage(support_j_lower:support_j_upper), source=.false.)
    end if

    do target = 1, distribution%root_tile_count()
      overlap_j_lower = max( &
        support_j_lower, distribution%root_tiles(target)%j_lower)
      overlap_j_upper = min( &
        support_j_upper, distribution%root_tiles(target)%j_upper)
      if (overlap_j_lower > overlap_j_upper) cycle
      source = distribution%root_tiles(target)%owner
      rows = overlap_j_upper - overlap_j_lower + 1
      cell_count = (support_i_upper - support_i_lower + 1) * rows
      state_count = component_count * cell_count
      value_count = 3 * (state_count + cell_count)
      if (source == destination) then
        if (distribution%rank == destination) then
          if (any(row_coverage(overlap_j_lower:overlap_j_upper))) return
          start_state(:, :, overlap_j_lower:overlap_j_upper) = &
            local_tile_states(target)%start_state( &
              :, support_i_lower:support_i_upper, &
              overlap_j_lower:overlap_j_upper)
          start_temperature(:, overlap_j_lower:overlap_j_upper) = &
            local_tile_states(target)%start_temperature( &
              support_i_lower:support_i_upper, &
              overlap_j_lower:overlap_j_upper)
          end_state(:, :, overlap_j_lower:overlap_j_upper) = &
            local_tile_states(target)%end_state( &
              :, support_i_lower:support_i_upper, &
              overlap_j_lower:overlap_j_upper)
          end_temperature(:, overlap_j_lower:overlap_j_upper) = &
            local_tile_states(target)%end_temperature( &
              support_i_lower:support_i_upper, &
              overlap_j_lower:overlap_j_upper)
          support_state(:, :, overlap_j_lower:overlap_j_upper) = &
            local_tile_states(target)%corrected_state( &
              :, support_i_lower:support_i_upper, &
              overlap_j_lower:overlap_j_upper)
          support_temperature(:, overlap_j_lower:overlap_j_upper) = &
            local_tile_states(target)%corrected_temperature( &
              support_i_lower:support_i_upper, &
              overlap_j_lower:overlap_j_upper)
          row_coverage(overlap_j_lower:overlap_j_upper) = .true.
        end if
      else if (distribution%rank == source) then
        allocate(payload(value_count))
        offset = 0
        payload(offset + 1:offset + state_count) = reshape( &
          local_tile_states(target)%start_state( &
            :, support_i_lower:support_i_upper, &
            overlap_j_lower:overlap_j_upper), [state_count])
        offset = offset + state_count
        payload(offset + 1:offset + cell_count) = reshape( &
          local_tile_states(target)%start_temperature( &
            support_i_lower:support_i_upper, &
            overlap_j_lower:overlap_j_upper), [cell_count])
        offset = offset + cell_count
        payload(offset + 1:offset + state_count) = reshape( &
          local_tile_states(target)%end_state( &
            :, support_i_lower:support_i_upper, &
            overlap_j_lower:overlap_j_upper), [state_count])
        offset = offset + state_count
        payload(offset + 1:offset + cell_count) = reshape( &
          local_tile_states(target)%end_temperature( &
            support_i_lower:support_i_upper, &
            overlap_j_lower:overlap_j_upper), [cell_count])
        offset = offset + cell_count
        payload(offset + 1:offset + state_count) = reshape( &
          local_tile_states(target)%corrected_state( &
            :, support_i_lower:support_i_upper, &
            overlap_j_lower:overlap_j_upper), [state_count])
        offset = offset + state_count
        payload(offset + 1:value_count) = reshape( &
          local_tile_states(target)%corrected_temperature( &
            support_i_lower:support_i_upper, &
            overlap_j_lower:overlap_j_upper), [cell_count])
        call MPI_Send( &
          payload, value_count, MPI_DOUBLE_PRECISION, destination, &
          sparse_child_state_support_tag, distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS) return
        local_transfers = local_transfers + 1
        deallocate(payload)
      else if (distribution%rank == destination) then
        if (any(row_coverage(overlap_j_lower:overlap_j_upper))) return
        allocate(payload(value_count))
        call MPI_Recv( &
          payload, value_count, MPI_DOUBLE_PRECISION, source, &
          sparse_child_state_support_tag, distribution%comm, status, ierr)
        if (ierr /= MPI_SUCCESS) return
        offset = 0
        start_state(:, :, overlap_j_lower:overlap_j_upper) = reshape( &
          payload(offset + 1:offset + state_count), &
          [component_count, support_i_upper - support_i_lower + 1, rows])
        offset = offset + state_count
        start_temperature(:, overlap_j_lower:overlap_j_upper) = reshape( &
          payload(offset + 1:offset + cell_count), &
          [support_i_upper - support_i_lower + 1, rows])
        offset = offset + cell_count
        end_state(:, :, overlap_j_lower:overlap_j_upper) = reshape( &
          payload(offset + 1:offset + state_count), &
          [component_count, support_i_upper - support_i_lower + 1, rows])
        offset = offset + state_count
        end_temperature(:, overlap_j_lower:overlap_j_upper) = reshape( &
          payload(offset + 1:offset + cell_count), &
          [support_i_upper - support_i_lower + 1, rows])
        offset = offset + cell_count
        support_state(:, :, overlap_j_lower:overlap_j_upper) = reshape( &
          payload(offset + 1:offset + state_count), &
          [component_count, support_i_upper - support_i_lower + 1, rows])
        offset = offset + state_count
        support_temperature(:, overlap_j_lower:overlap_j_upper) = reshape( &
          payload(offset + 1:value_count), &
          [support_i_upper - support_i_lower + 1, rows])
        row_coverage(overlap_j_lower:overlap_j_upper) = .true.
        deallocate(payload)
      end if
    end do

    local_ok = .true.
    if (distribution%rank == destination) then
      local_ok = all(row_coverage) .and. &
        all(ieee_is_finite(start_state)) .and. &
        all(ieee_is_finite(start_temperature)) .and. &
        all(ieee_is_finite(end_state)) .and. &
        all(ieee_is_finite(end_temperature)) .and. &
        all(ieee_is_finite(support_state)) .and. &
        all(ieee_is_finite(support_temperature))
      if (local_ok) call extract_reactive_eb_patch_exterior_context_support_2d( &
        support_i_lower, support_j_lower, start_state, start_temperature, &
        end_state, end_temperature, coarse_geometry, fine_geometry, &
        patch, component_count, context, local_ok)
      deallocate( &
        start_state, start_temperature, end_state, end_temperature, &
        row_coverage)
    end if
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    ok = global_ok .and. accepted
  end subroutine transfer_sparse_child_state_support_2d

  subroutine transfer_sparse_child_state_correction_2d( &
      distribution, coarse_geometry, fine_geometry, patch, component_count, &
      source, support_state, support_temperature, local_tile_states, &
      local_transfers, ok)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(eb_geometry_2d), intent(in) :: coarse_geometry, fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    integer, intent(in) :: component_count, source
    real(dp), allocatable, intent(in) :: support_state(:, :, :)
    real(dp), allocatable, intent(in) :: support_temperature(:, :)
    type(mpi_amr_eb_root_tile_transport_state_2d), intent(inout) :: &
      local_tile_states(:)
    integer, intent(inout) :: local_transfers
    logical, intent(out) :: ok

    type(MPI_Status) :: status
    real(dp), allocatable :: payload(:)
    logical :: accepted, global_ok, local_ok
    integer :: cell_count, destination, ierr, overlap_j_lower
    integer :: overlap_j_upper, rows, state_count, support_i_lower
    integer :: support_i_upper, support_j_lower, support_j_upper
    integer :: target, value_count

    ok = .false.
    support_i_lower = max(1, patch%coarse_i_lower - 2)
    support_i_upper = min(coarse_geometry%nx, patch%coarse_i_upper + 2)
    support_j_lower = max(1, patch%coarse_j_lower - 2)
    support_j_upper = min(coarse_geometry%ny, patch%coarse_j_upper + 2)
    local_ok = component_count >= 1 .and. source >= 0 .and. &
      source < distribution%nranks .and. &
      patch%is_valid(coarse_geometry, fine_geometry) .and. &
      size(local_tile_states) == distribution%root_tile_count()
    if (distribution%rank == source .and. local_ok) local_ok = &
      allocated(support_state) .and. allocated(support_temperature)
    if (distribution%rank == source .and. local_ok) local_ok = &
      all(shape(support_state) == &
        [component_count, support_i_upper - support_i_lower + 1, &
         support_j_upper - support_j_lower + 1]) .and. &
      all(shape(support_temperature) == &
        [support_i_upper - support_i_lower + 1, &
         support_j_upper - support_j_lower + 1]) .and. &
      lbound(support_state, 2) == support_i_lower .and. &
      lbound(support_state, 3) == support_j_lower .and. &
      lbound(support_temperature, 1) == support_i_lower .and. &
      lbound(support_temperature, 2) == support_j_lower .and. &
      all(ieee_is_finite(support_state)) .and. &
      all(ieee_is_finite(support_temperature))
    do target = 1, distribution%root_tile_count()
      if (distribution%rank /= distribution%root_tiles(target)%owner) cycle
      rows = distribution%root_tiles(target)%j_upper - &
        distribution%root_tiles(target)%j_lower + 1
      local_ok = local_ok .and. &
        allocated(local_tile_states(target)%corrected_state) .and. &
        allocated(local_tile_states(target)%corrected_temperature)
      if (.not. local_ok) cycle
      local_ok = &
        all(shape(local_tile_states(target)%corrected_state) == &
          [component_count, coarse_geometry%nx, rows]) .and. &
        all(shape(local_tile_states(target)%corrected_temperature) == &
          [coarse_geometry%nx, rows]) .and. &
        lbound(local_tile_states(target)%corrected_state, 2) == 1 .and. &
        lbound(local_tile_states(target)%corrected_state, 3) == &
          distribution%root_tiles(target)%j_lower .and. &
        lbound(local_tile_states(target)%corrected_temperature, 1) == 1 .and. &
        lbound(local_tile_states(target)%corrected_temperature, 2) == &
          distribution%root_tiles(target)%j_lower .and. &
        all(ieee_is_finite(local_tile_states(target)%corrected_state)) .and. &
        all(ieee_is_finite( &
          local_tile_states(target)%corrected_temperature))
    end do
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    do target = 1, distribution%root_tile_count()
      overlap_j_lower = max( &
        support_j_lower, distribution%root_tiles(target)%j_lower)
      overlap_j_upper = min( &
        support_j_upper, distribution%root_tiles(target)%j_upper)
      if (overlap_j_lower > overlap_j_upper) cycle
      destination = distribution%root_tiles(target)%owner
      rows = overlap_j_upper - overlap_j_lower + 1
      cell_count = (support_i_upper - support_i_lower + 1) * rows
      state_count = component_count * cell_count
      value_count = state_count + cell_count
      if (source == destination) then
        if (distribution%rank == source) then
          local_tile_states(target)%corrected_state( &
            :, support_i_lower:support_i_upper, &
            overlap_j_lower:overlap_j_upper) = &
              support_state(:, :, overlap_j_lower:overlap_j_upper)
          local_tile_states(target)%corrected_temperature( &
            support_i_lower:support_i_upper, &
            overlap_j_lower:overlap_j_upper) = &
              support_temperature(:, overlap_j_lower:overlap_j_upper)
        end if
      else if (distribution%rank == source) then
        allocate(payload(value_count))
        payload(1:state_count) = reshape( &
          support_state(:, :, overlap_j_lower:overlap_j_upper), [state_count])
        payload(state_count + 1:value_count) = reshape( &
          support_temperature(:, overlap_j_lower:overlap_j_upper), &
          [cell_count])
        call MPI_Send( &
          payload, value_count, MPI_DOUBLE_PRECISION, destination, &
          sparse_child_state_correction_tag, distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS) return
        local_transfers = local_transfers + 1
        deallocate(payload)
      else if (distribution%rank == destination) then
        allocate(payload(value_count))
        call MPI_Recv( &
          payload, value_count, MPI_DOUBLE_PRECISION, source, &
          sparse_child_state_correction_tag, distribution%comm, status, ierr)
        if (ierr /= MPI_SUCCESS) return
        local_tile_states(target)%corrected_state( &
          :, support_i_lower:support_i_upper, &
          overlap_j_lower:overlap_j_upper) = reshape( &
            payload(1:state_count), &
            [component_count, support_i_upper - support_i_lower + 1, rows])
        local_tile_states(target)%corrected_temperature( &
          support_i_lower:support_i_upper, &
          overlap_j_lower:overlap_j_upper) = reshape( &
            payload(state_count + 1:value_count), &
            [support_i_upper - support_i_lower + 1, rows])
        deallocate(payload)
      end if
    end do

    local_ok = .true.
    do target = 1, distribution%root_tile_count()
      if (distribution%rank /= distribution%root_tiles(target)%owner) cycle
      local_ok = local_ok .and. &
        all(ieee_is_finite(local_tile_states(target)%corrected_state)) .and. &
        all(ieee_is_finite( &
          local_tile_states(target)%corrected_temperature))
    end do
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    ok = global_ok .and. accepted
  end subroutine transfer_sparse_child_state_correction_2d


  pure logical function coarse_cell_is_refined_2d( &
      patch_set_template, i, j) result(refined)
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set_template
    integer, intent(in) :: i, j

    integer :: child

    refined = .false.
    do child = 1, patch_set_template%patch_count()
      refined = &
        i >= patch_set_template%children(child)%patch%coarse_i_lower .and. &
        i <= patch_set_template%children(child)%patch%coarse_i_upper .and. &
        j >= patch_set_template%children(child)%patch%coarse_j_lower .and. &
        j <= patch_set_template%children(child)%patch%coarse_j_upper
      if (refined) return
    end do
  end function coarse_cell_is_refined_2d

  subroutine composite_sparse_owned_reactive_eb_patch_set_integral_2d( &
      distribution, sparse_patch_set, coarse_geometry, patch_set_template, &
      integral, ok)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(in) :: sparse_patch_set
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set_template
    real(dp), intent(out) :: integral(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: local_integral(:)
    logical :: accepted, global_ok, local_ok
    integer :: child, global_j, i, ierr, j, local_j, tile

    integral = 0.0_dp
    ok = .false.
    local_ok = size(integral) == sparse_patch_set%nvar .and. &
      sparse_patch_set%is_valid( &
        distribution, coarse_geometry, patch_set_template)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    allocate(local_integral(size(integral)), source=0.0_dp)

    do tile = 1, distribution%root_tile_count()
      if (.not. distribution%root_tile_is_local(tile)) cycle
      do local_j = 1, size(sparse_patch_set%root_tiles(tile)%state, 3)
        global_j = distribution%root_tiles(tile)%j_lower + local_j - 1
        do i = 1, coarse_geometry%nx
          if (coarse_cell_is_refined_2d( &
              patch_set_template, i, global_j)) cycle
          local_integral = local_integral + &
            coarse_geometry%volume_fraction(i, global_j) * &
            sparse_patch_set%root_tiles(tile)%state(:, i, local_j) * &
            coarse_geometry%dx * coarse_geometry%dy
        end do
      end do
    end do
    do child = 1, distribution%child_count()
      if (.not. distribution%child_is_local(child)) cycle
      do j = 1, patch_set_template%children(child)%geometry%ny
        do i = 1, patch_set_template%children(child)%geometry%nx
          local_integral = local_integral + &
            patch_set_template%children(child)%geometry% &
              volume_fraction(i, j) * &
            sparse_patch_set%children(child)%state(:, i, j) * &
            patch_set_template%children(child)%geometry%dx * &
            patch_set_template%children(child)%geometry%dy
        end do
      end do
    end do
    call MPI_Allreduce( &
      local_integral, integral, size(integral), MPI_DOUBLE_PRECISION, &
      MPI_SUM, distribution%comm, ierr)
    local_ok = ierr == MPI_SUCCESS .and. all(ieee_is_finite(integral))
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    ok = global_ok .and. accepted
    if (.not. ok) integral = 0.0_dp
  end subroutine composite_sparse_owned_reactive_eb_patch_set_integral_2d

  subroutine close_sparse_cut_patch_set_conservation_2d( &
      species, integral_before, distribution, sparse_patch_set, &
      coarse_geometry, patch_set_template, boundary_change, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: integral_before(:)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(inout) :: sparse_patch_set
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set_template
    real(dp), intent(in) :: boundary_change(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: correction(:)
    real(dp), allocatable :: current_integral(:), primitive(:), residual(:)
    logical, allocatable :: refined(:, :), recipients(:, :)
    real(dp) :: closure_tolerance, recipient_volume, recovered_temperature
    real(dp) :: scale, sound_speed, species_residual
    logical :: accepted, entity_ok, global_ok, local_ok
    integer :: child, component, global_j, i, j, k, local_j, nvar, tile

    ok = .false.
    nvar = reactive_nvar(size(species))
    local_ok = nvar >= 1 .and. size(integral_before) == nvar .and. &
      size(boundary_change) == nvar .and. &
      all(ieee_is_finite(boundary_change)) .and. &
      sparse_patch_set%is_valid( &
        distribution, coarse_geometry, patch_set_template)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    allocate(current_integral(nvar), residual(nvar), correction(nvar))
    call composite_sparse_owned_reactive_eb_patch_set_integral_2d( &
      distribution, sparse_patch_set, coarse_geometry, patch_set_template, &
      current_integral, local_ok)
    if (.not. local_ok) return
    residual = integral_before + boundary_change - current_integral
    correction = residual
    species_residual = 0.0_dp
    do k = 1, size(species)
      component = reactive_species_component(k)
      species_residual = species_residual + residual(component)
    end do
    scale = max(1.0_dp, maxval(abs(residual)), abs(species_residual))
    closure_tolerance = 4096.0_dp * epsilon(1.0_dp) * scale
    local_ok = abs(residual(irho) - species_residual) <= closure_tolerance
    if (.not. local_ok) return
    component = reactive_species_component(size(species))
    correction(component) = correction(component) + &
      residual(irho) - species_residual

    allocate(refined(coarse_geometry%nx, coarse_geometry%ny), &
      recipients(coarse_geometry%nx, coarse_geometry%ny))
    refined = .false.
    do child = 1, patch_set_template%patch_count()
      refined( &
        patch_set_template%children(child)%patch%coarse_i_lower: &
          patch_set_template%children(child)%patch%coarse_i_upper, &
        patch_set_template%children(child)%patch%coarse_j_lower: &
          patch_set_template%children(child)%patch%coarse_j_upper) = .true.
    end do
    recipients = .false.
    local_ok = .true.
    do child = 1, patch_set_template%patch_count()
      call mark_local_coarse_fine_interface_recipients_2d( &
        coarse_geometry, patch_set_template%children(child)%geometry, &
        patch_set_template%children(child)%patch, refined, recipients, &
        entity_ok)
      local_ok = local_ok .and. entity_ok
    end do
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    recipient_volume = 0.0_dp
    do j = 1, coarse_geometry%ny
      do i = 1, coarse_geometry%nx
        if (.not. recipients(i, j)) cycle
        recipient_volume = recipient_volume + &
          coarse_geometry%volume_fraction(i, j) * &
          coarse_geometry%dx * coarse_geometry%dy
      end do
    end do
    local_ok = ieee_is_finite(recipient_volume) .and. &
      recipient_volume > tiny(1.0_dp)
    if (.not. local_ok) return
    correction = correction / recipient_volume
    if (any(.not. ieee_is_finite(correction))) return

    allocate(primitive(reactive_nprim(size(species))))
    entity_ok = .true.
    do tile = 1, distribution%root_tile_count()
      if (.not. distribution%root_tile_is_local(tile)) cycle
      do local_j = 1, size(sparse_patch_set%root_tiles(tile)%state, 3)
        global_j = distribution%root_tiles(tile)%j_lower + local_j - 1
        do i = 1, coarse_geometry%nx
          if (.not. recipients(i, global_j)) cycle
          sparse_patch_set%root_tiles(tile)%state(:, i, local_j) = &
            sparse_patch_set%root_tiles(tile)%state(:, i, local_j) + &
            correction
          call reactive_conserved_to_primitive( &
            species, sparse_patch_set%root_tiles(tile)%state(:, i, local_j), &
            sparse_patch_set%root_tiles(tile)%temperature(i, local_j), &
            primitive, recovered_temperature, sound_speed, local_ok)
          if (.not. local_ok) then
            entity_ok = .false.
            exit
          end if
          sparse_patch_set%root_tiles(tile)%temperature(i, local_j) = &
            recovered_temperature
        end do
        if (.not. entity_ok) exit
      end do
      if (.not. entity_ok) exit
    end do
    call all_ranks_accept_eb_2d( &
      distribution, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    local_ok = sparse_patch_set%is_valid( &
      distribution, coarse_geometry, patch_set_template)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    ok = global_ok .and. accepted
  end subroutine close_sparse_cut_patch_set_conservation_2d

  subroutine transfer_sparse_restriction_to_root_owners_2d( &
      distribution, child_owner, coarse_j_lower, coarse_j_upper, nvar, &
      coarse_width, restricted_state, local_transfers, ok)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    integer, intent(in) :: child_owner, coarse_j_lower, coarse_j_upper
    integer, intent(in) :: nvar, coarse_width
    real(dp), allocatable, intent(inout) :: restricted_state(:, :, :)
    integer, intent(inout) :: local_transfers
    logical, intent(out) :: ok

    type(MPI_Status) :: status
    logical, allocatable :: recipients(:)
    integer :: ierr, recipient, tile, value_count

    ok = .false.
    if (child_owner < 0 .or. child_owner >= distribution%nranks .or. &
        coarse_j_lower < 1 .or. coarse_j_upper < coarse_j_lower .or. &
        nvar < 1 .or. coarse_width < 1) return
    allocate(recipients(distribution%nranks), source=.false.)
    do tile = 1, distribution%root_tile_count()
      if (distribution%root_tiles(tile)%j_upper < coarse_j_lower .or. &
          distribution%root_tiles(tile)%j_lower > coarse_j_upper) cycle
      recipient = distribution%root_tiles(tile)%owner
      if (recipient < 0 .or. recipient >= distribution%nranks) return
      recipients(recipient + 1) = .true.
    end do
    if (.not. any(recipients)) return
    value_count = nvar * coarse_width * &
      (coarse_j_upper - coarse_j_lower + 1)

    if (distribution%rank == child_owner) then
      if (.not. allocated(restricted_state)) return
      if (size(restricted_state) /= value_count) return
      do recipient = 0, distribution%nranks - 1
        if (recipient == child_owner .or. &
            .not. recipients(recipient + 1)) cycle
        call MPI_Send( &
          restricted_state, value_count, MPI_DOUBLE_PRECISION, recipient, &
          sparse_restriction_tag, distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS) return
        local_transfers = local_transfers + 1
      end do
    else if (recipients(distribution%rank + 1)) then
      allocate(restricted_state( &
        nvar, coarse_width, coarse_j_upper - coarse_j_lower + 1))
      call MPI_Recv( &
        restricted_state, value_count, MPI_DOUBLE_PRECISION, child_owner, &
        sparse_restriction_tag, distribution%comm, status, ierr)
      if (ierr /= MPI_SUCCESS) return
    end if
    ok = .true.
  end subroutine transfer_sparse_restriction_to_root_owners_2d

  subroutine average_down_sparse_owned_reactive_eb_patch_set_2d( &
      species, distribution, sparse_patch_set, coarse_geometry, &
      patch_set_template, ok, local_restriction_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(inout) :: sparse_patch_set
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set_template
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_restriction_transfers

    type(mpi_amr_eb_sparse_patch_set_2d) :: backup, candidate
    real(dp), allocatable :: primitive(:), restricted_state(:, :, :)
    real(dp) :: fine_volume, recovered_temperature, sound_speed
    logical :: accepted, entity_ok, global_ok, local_ok
    integer :: child, coarse_i, coarse_i_lower, coarse_i_upper, coarse_j
    integer :: coarse_j_lower, coarse_j_upper, component, fine_i_lower
    integer :: fine_i_upper, fine_j_lower, fine_j_upper, ierr, local_i
    integer :: local_j, nspecies, nspecies_maximum, nspecies_minimum
    integer :: owner, ratio, tile, transfers

    ok = .false.
    transfers = 0
    if (present(local_restriction_transfers)) &
      local_restriction_transfers = 0
    local_ok = size(species) >= 1 .and. &
      sparse_patch_set%nvar == reactive_nvar(size(species)) .and. &
      sparse_patch_set%is_valid( &
        distribution, coarse_geometry, patch_set_template)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    nspecies = size(species)
    call MPI_Allreduce( &
      nspecies, nspecies_minimum, 1, MPI_INTEGER, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      nspecies, nspecies_maximum, 1, MPI_INTEGER, MPI_MAX, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. nspecies_minimum /= nspecies_maximum) return

    backup = sparse_patch_set
    candidate = sparse_patch_set
    allocate(primitive(reactive_nprim(size(species))))
    do child = 1, distribution%child_count()
      coarse_i_lower = &
        patch_set_template%children(child)%patch%coarse_i_lower
      coarse_i_upper = &
        patch_set_template%children(child)%patch%coarse_i_upper
      coarse_j_lower = &
        patch_set_template%children(child)%patch%coarse_j_lower
      coarse_j_upper = &
        patch_set_template%children(child)%patch%coarse_j_upper
      owner = distribution%child_owner(child)
      ratio = &
        patch_set_template%children(child)%patch%refinement_ratio
      if (distribution%rank == owner) then
        allocate(restricted_state( &
          candidate%nvar, coarse_i_upper - coarse_i_lower + 1, &
          coarse_j_upper - coarse_j_lower + 1), source=0.0_dp)
        do coarse_j = coarse_j_lower, coarse_j_upper
          local_j = coarse_j - coarse_j_lower + 1
          fine_j_lower = (coarse_j - coarse_j_lower) * ratio + 1
          fine_j_upper = fine_j_lower + ratio - 1
          do coarse_i = coarse_i_lower, coarse_i_upper
            local_i = coarse_i - coarse_i_lower + 1
            fine_i_lower = (coarse_i - coarse_i_lower) * ratio + 1
            fine_i_upper = fine_i_lower + ratio - 1
            fine_volume = sum( &
              patch_set_template%children(child)%geometry%volume_fraction( &
                fine_i_lower:fine_i_upper, fine_j_lower:fine_j_upper))
            if (fine_volume > tiny(1.0_dp)) then
              do component = 1, candidate%nvar
                restricted_state(component, local_i, local_j) = sum( &
                    patch_set_template%children(child)%geometry% &
                      volume_fraction( &
                        fine_i_lower:fine_i_upper, &
                        fine_j_lower:fine_j_upper) * &
                    candidate%children(child)%state( &
                      component, fine_i_lower:fine_i_upper, &
                      fine_j_lower:fine_j_upper)) / fine_volume
              end do
            else
              restricted_state(:, local_i, local_j) = &
                  candidate%children(child)%state( &
                    :, fine_i_lower, fine_j_lower)
            end if
          end do
        end do
      end if
      call transfer_sparse_restriction_to_root_owners_2d( &
        distribution, owner, coarse_j_lower, coarse_j_upper, &
        candidate%nvar, coarse_i_upper - coarse_i_lower + 1, &
        restricted_state, transfers, local_ok)
      if (.not. local_ok) then
        sparse_patch_set = backup
        return
      end if

      entity_ok = .true.
      if (allocated(restricted_state)) &
        entity_ok = all(ieee_is_finite(restricted_state))
      do tile = 1, distribution%root_tile_count()
        if (.not. distribution%root_tile_is_local(tile)) cycle
        if (distribution%root_tiles(tile)%j_upper < coarse_j_lower .or. &
            distribution%root_tiles(tile)%j_lower > coarse_j_upper) cycle
        if (.not. allocated(restricted_state)) then
          entity_ok = .false.
          cycle
        end if
        do coarse_j = max( &
            distribution%root_tiles(tile)%j_lower, &
            coarse_j_lower), &
            min( &
              distribution%root_tiles(tile)%j_upper, &
              coarse_j_upper)
          local_j = coarse_j - distribution%root_tiles(tile)%j_lower + 1
          do coarse_i = coarse_i_lower, coarse_i_upper
            local_i = coarse_i - coarse_i_lower + 1
            if (coarse_geometry%cell_type(coarse_i, coarse_j) == &
                eb_covered_cell) cycle
            if (.not. entity_ok) cycle
            if (candidate%root_tiles(tile)%temperature( &
                coarse_i, local_j) <= 0.0_dp) then
              entity_ok = .false.
              cycle
            end if
            call reactive_conserved_to_primitive( &
              species, restricted_state(:, local_i, &
                coarse_j - coarse_j_lower + 1), &
              candidate%root_tiles(tile)%temperature(coarse_i, local_j), &
              primitive, recovered_temperature, sound_speed, local_ok)
            if (.not. local_ok) then
              entity_ok = .false.
              cycle
            end if
            candidate%root_tiles(tile)%state(:, coarse_i, local_j) = &
              restricted_state( &
                :, local_i, coarse_j - coarse_j_lower + 1)
            candidate%root_tiles(tile)%temperature(coarse_i, local_j) = &
              recovered_temperature
          end do
        end do
      end do
      call all_ranks_accept_eb_2d( &
        distribution, entity_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) then
        sparse_patch_set = backup
        return
      end if
      if (allocated(restricted_state)) deallocate(restricted_state)
    end do

    local_ok = candidate%is_valid( &
      distribution, coarse_geometry, patch_set_template)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) then
      sparse_patch_set = backup
      return
    end if
    sparse_patch_set = candidate
    ok = .true.
    if (present(local_restriction_transfers)) &
      local_restriction_transfers = transfers
  end subroutine average_down_sparse_owned_reactive_eb_patch_set_2d

  subroutine advance_sparse_owned_reactive_eb_patch_set_chemistry_2d( &
      species, reactions, interval, relative_tolerance, absolute_tolerance, &
      distribution, sparse_patch_set, coarse_geometry, patch_set_template, &
      ok, local_entity_advances, local_restriction_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: interval, relative_tolerance, absolute_tolerance
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(inout) :: sparse_patch_set
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set_template
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_entity_advances
    integer, intent(out), optional :: local_restriction_transfers

    type(mpi_amr_eb_sparse_patch_set_2d) :: backup, candidate
    real(dp) :: controls(3), control_maximum(3), control_minimum(3)
    logical, allocatable :: active_mask(:, :)
    logical :: accepted, entity_ok, global_ok, local_ok
    integer :: advances, average_down_transfers, child, ierr
    integer :: integer_controls(2)
    integer :: integer_maximum(2), integer_minimum(2)
    integer :: j_lower, j_upper, tile

    ok = .false.
    advances = 0
    if (present(local_entity_advances)) local_entity_advances = 0
    if (present(local_restriction_transfers)) &
      local_restriction_transfers = 0
    controls = [interval, relative_tolerance, absolute_tolerance]
    integer_controls = [size(species), size(reactions)]
    local_ok = interval >= 0.0_dp .and. relative_tolerance > 0.0_dp .and. &
      absolute_tolerance > 0.0_dp .and. all(ieee_is_finite(controls)) .and. &
      size(species) >= 1 .and. size(reactions) >= 1 .and. &
      sparse_patch_set%is_valid( &
        distribution, coarse_geometry, patch_set_template)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Allreduce( &
      controls, control_minimum, 3, MPI_DOUBLE_PRECISION, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      controls, control_maximum, 3, MPI_DOUBLE_PRECISION, MPI_MAX, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_controls, integer_minimum, 2, MPI_INTEGER, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_controls, integer_maximum, 2, MPI_INTEGER, MPI_MAX, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. &
        any(control_minimum /= control_maximum) .or. &
        any(integer_minimum /= integer_maximum)) return

    backup = sparse_patch_set
    candidate = sparse_patch_set
    do tile = 1, distribution%root_tile_count()
      entity_ok = .true.
      if (distribution%root_tile_is_local(tile)) then
        j_lower = distribution%root_tiles(tile)%j_lower
        j_upper = distribution%root_tiles(tile)%j_upper
        allocate(active_mask(coarse_geometry%nx, j_upper - j_lower + 1))
        active_mask = coarse_geometry%cell_type(:, j_lower:j_upper) /= &
          eb_covered_cell
        call advance_reactive_chemistry_2d( &
          species, reactions, candidate%root_tiles(tile)%state, &
          candidate%root_tiles(tile)%temperature, coarse_geometry%nx, &
          j_upper - j_lower + 1, interval, relative_tolerance, &
          absolute_tolerance, entity_ok, active_mask)
        deallocate(active_mask)
        if (entity_ok) advances = advances + 1
      end if
      call all_ranks_accept_eb_2d( &
        distribution, entity_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) then
        sparse_patch_set = backup
        return
      end if
    end do
    do child = 1, distribution%child_count()
      entity_ok = .true.
      if (distribution%child_is_local(child)) then
        allocate(active_mask( &
          patch_set_template%children(child)%geometry%nx, &
          patch_set_template%children(child)%geometry%ny))
        active_mask = &
          patch_set_template%children(child)%geometry%cell_type /= &
            eb_covered_cell
        call advance_reactive_chemistry_2d( &
          species, reactions, candidate%children(child)%state, &
          candidate%children(child)%temperature, &
          patch_set_template%children(child)%geometry%nx, &
          patch_set_template%children(child)%geometry%ny, interval, &
          relative_tolerance, absolute_tolerance, entity_ok, active_mask)
        deallocate(active_mask)
        if (entity_ok) advances = advances + 1
      end if
      call all_ranks_accept_eb_2d( &
        distribution, entity_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) then
        sparse_patch_set = backup
        return
      end if
    end do

    call average_down_sparse_owned_reactive_eb_patch_set_2d( &
      species, distribution, candidate, coarse_geometry, &
      patch_set_template, local_ok, average_down_transfers)
    if (.not. local_ok) then
      sparse_patch_set = backup
      return
    end if
    sparse_patch_set = candidate
    ok = .true.
    if (present(local_entity_advances)) local_entity_advances = advances
    if (present(local_restriction_transfers)) &
      local_restriction_transfers = average_down_transfers
  end subroutine advance_sparse_owned_reactive_eb_patch_set_chemistry_2d

  subroutine compute_sparse_owned_reactive_eb_patch_set_timestep_2d( &
      species, transport, distribution, sparse_patch_set, coarse_geometry, &
      patch_set_template, hydro_cfl, transport_cfl, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, dt, ok, local_root_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(in) :: sparse_patch_set
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set_template
    real(dp), intent(in) :: hydro_cfl, transport_cfl
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_root_transfers

    type(eb_geometry_2d) :: root_tile_geometry
    real(dp) :: entity_dt, global_dt, local_dt, maximum_diffusivity
    real(dp) :: numeric_controls(2), numeric_maximum(2), numeric_minimum(2)
    logical :: accepted, entity_ok, global_ok, local_ok, transport_active
    integer :: child, ierr, integer_controls(5), integer_maximum(5)
    integer :: integer_minimum(5), ratio, tile, transfers

    dt = 0.0_dp
    ok = .false.
    transfers = 0
    if (present(local_root_transfers)) local_root_transfers = 0
    numeric_controls = [hydro_cfl, transport_cfl]
    integer_controls = [ &
      size(species), merge(1, 0, viscosity_enabled), &
      merge(1, 0, thermal_conduction_enabled), &
      merge(1, 0, species_diffusion_enabled), &
      merge(1, 0, barodiffusion_enabled)]
    local_ok = size(species) >= 1 .and. &
      all(ieee_is_finite(numeric_controls)) .and. &
      hydro_cfl > 0.0_dp .and. hydro_cfl <= 1.0_dp .and. &
      transport_cfl > 0.0_dp .and. transport_cfl <= 0.5_dp .and. &
      sparse_patch_set%nvar == reactive_nvar(size(species)) .and. &
      sparse_patch_set%is_valid( &
        distribution, coarse_geometry, patch_set_template)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Allreduce( &
      numeric_controls, numeric_minimum, 2, MPI_DOUBLE_PRECISION, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      numeric_controls, numeric_maximum, 2, MPI_DOUBLE_PRECISION, MPI_MAX, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_controls, integer_minimum, 5, MPI_INTEGER, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_controls, integer_maximum, 5, MPI_INTEGER, MPI_MAX, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. &
        any(numeric_minimum /= numeric_maximum) .or. &
        any(integer_minimum /= integer_maximum)) return
    call collective_transport_preflight_2d( &
      species, transport, distribution, coarse_geometry, patch_set_template, &
      0.0_dp, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, 2, &
      0.5_dp, local_ok)
    if (.not. local_ok) return

    transport_active = viscosity_enabled .or. thermal_conduction_enabled .or. &
      species_diffusion_enabled
    local_dt = huge(1.0_dp)
    entity_ok = .true.
    do tile = 1, distribution%root_tile_count()
      if (.not. distribution%root_tile_is_local(tile)) cycle
      call extract_eb_geometry_y_band_2d( &
        coarse_geometry, distribution%root_tiles(tile)%j_lower, &
        distribution%root_tiles(tile)%j_upper, root_tile_geometry, entity_ok)
      if (.not. entity_ok) exit
      if (count(root_tile_geometry%cell_type /= eb_covered_cell) == 0) cycle
      call compute_reactive_eb_cfl_timestep_2d( &
        species, sparse_patch_set%root_tiles(tile)%state, &
        sparse_patch_set%root_tiles(tile)%temperature, root_tile_geometry, &
        hydro_cfl, entity_dt, entity_ok)
      if (entity_ok) local_dt = min(local_dt, entity_dt)
      if (entity_ok .and. transport_active) then
        call reactive_eb_transport_timestep_2d( &
          species, transport, sparse_patch_set%root_tiles(tile)%state, &
          sparse_patch_set%root_tiles(tile)%temperature, root_tile_geometry, &
          transport_cfl, viscosity_enabled, thermal_conduction_enabled, &
          species_diffusion_enabled, entity_dt, maximum_diffusivity, entity_ok)
        if (entity_ok) local_dt = min(local_dt, entity_dt)
      end if
      if (.not. entity_ok) exit
    end do
    do child = 1, distribution%child_count()
      if (.not. distribution%child_is_local(child)) cycle
      if (count(patch_set_template%children(child)%geometry%cell_type /= &
          eb_covered_cell) == 0) cycle
      call compute_reactive_eb_cfl_timestep_2d( &
        species, sparse_patch_set%children(child)%state, &
        sparse_patch_set%children(child)%temperature, &
        patch_set_template%children(child)%geometry, hydro_cfl, entity_dt, &
        entity_ok)
      if (.not. entity_ok) exit
      ratio = patch_set_template%children(child)%patch%refinement_ratio
      local_dt = min(local_dt, real(ratio, dp) * entity_dt)
      if (transport_active) then
        call reactive_eb_transport_timestep_2d( &
          species, transport, sparse_patch_set%children(child)%state, &
          sparse_patch_set%children(child)%temperature, &
          patch_set_template%children(child)%geometry, transport_cfl, &
          viscosity_enabled, thermal_conduction_enabled, &
          species_diffusion_enabled, entity_dt, maximum_diffusivity, entity_ok)
        if (.not. entity_ok) exit
        local_dt = min(local_dt, real(ratio, dp) * entity_dt)
      end if
    end do
    local_ok = entity_ok .and. ieee_is_finite(local_dt) .and. &
      local_dt > 0.0_dp
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Allreduce( &
      local_dt, global_dt, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. .not. ieee_is_finite(global_dt) .or. &
        global_dt <= 0.0_dp .or. global_dt >= huge(1.0_dp)) return

    dt = global_dt
    ok = .true.
    if (present(local_root_transfers)) local_root_transfers = transfers
  end subroutine compute_sparse_owned_reactive_eb_patch_set_timestep_2d

  subroutine advance_sparse_owned_reactive_eb_patch_set_hydro_2d( &
      species, distribution, sparse_patch_set, coarse_geometry, &
      patch_set_template, solver, reconstruction, limiter, &
      state_redist_max_order, dt, ok, local_level_advances, &
      state_redist_target_volume_fraction, local_root_transfers, &
      local_root_hydro_cells)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(inout) :: sparse_patch_set
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set_template
    character(len=*), intent(in) :: solver, reconstruction, limiter
    integer, intent(in) :: state_redist_max_order
    real(dp), intent(in) :: dt
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_level_advances
    real(dp), intent(in), optional :: state_redist_target_volume_fraction
    integer, intent(out), optional :: local_root_transfers
    integer, intent(out), optional :: local_root_hydro_cells

    type(amr_eb_flux_register_2d) :: flux_register
    type(reactive_eb_exterior_state_2d) :: exterior
    type(reactive_eb_patch_exterior_context_2d) :: exterior_context
    type(mpi_amr_eb_sparse_patch_set_2d) :: candidate
    type(mpi_amr_eb_root_tile_transport_state_2d), allocatable :: &
      local_tile_states(:)
    type(mpi_amr_eb_root_tile_transport_flux_2d), allocatable :: &
      local_tile_fluxes(:)
    real(dp), allocatable :: coarse_support(:, :, :)
    real(dp), allocatable :: coarse_support_temperature(:, :)
    real(dp), allocatable :: coarse_work(:, :, :)
    real(dp), allocatable :: coarse_work_temperature(:, :)
    real(dp), allocatable :: coarse_x_flux_support(:, :, :)
    real(dp), allocatable :: coarse_y_flux_support(:, :, :)
    real(dp), allocatable :: fine_work(:, :, :)
    real(dp), allocatable :: fine_work_temperature(:, :)
    real(dp), allocatable :: fine_x_flux(:, :, :)
    real(dp), allocatable :: fine_y_flux(:, :, :)
    real(dp) :: alpha, fine_dt, numeric_controls(2)
    real(dp) :: numeric_maximum(2), numeric_minimum(2), selected_target
    logical :: accepted, entity_ok, global_ok, local_ok
    integer :: advances, character_index, child, ierr, integer_controls(2)
    integer :: integer_maximum(2), integer_minimum(2)
    integer :: nvar, owner, ratio, root_hydro_cells, root_tile_advances
    integer :: substep, tile, transfers
    integer :: flux_x_i_lower, flux_x_j_lower
    integer :: flux_y_i_lower, flux_y_j_lower
    integer :: support_i_lower, support_i_upper
    integer :: support_j_lower, support_j_upper
    integer :: string_codes(32, 3), string_maximum(32, 3)
    integer :: string_minimum(32, 3)

    ok = .false.
    advances = 0
    transfers = 0
    root_hydro_cells = 0
    if (present(local_level_advances)) local_level_advances = 0
    if (present(local_root_transfers)) local_root_transfers = 0
    if (present(local_root_hydro_cells)) local_root_hydro_cells = 0
    selected_target = 0.5_dp
    if (present(state_redist_target_volume_fraction)) &
      selected_target = state_redist_target_volume_fraction
    nvar = reactive_nvar(size(species))
    numeric_controls = [dt, selected_target]
    integer_controls = [state_redist_max_order, size(species)]
    string_codes = 0
    local_ok = len_trim(solver) >= 1 .and. len_trim(solver) <= 32 .and. &
      len_trim(reconstruction) >= 1 .and. &
      len_trim(reconstruction) <= 32 .and. &
      len_trim(limiter) >= 1 .and. len_trim(limiter) <= 32
    if (local_ok) then
      do character_index = 1, len_trim(solver)
        string_codes(character_index, 1) = &
          iachar(solver(character_index:character_index))
      end do
      do character_index = 1, len_trim(reconstruction)
        string_codes(character_index, 2) = &
          iachar(reconstruction(character_index:character_index))
      end do
      do character_index = 1, len_trim(limiter)
        string_codes(character_index, 3) = &
          iachar(limiter(character_index:character_index))
      end do
    end if
    local_ok = local_ok .and. size(species) >= 1 .and. &
      all(ieee_is_finite(numeric_controls)) .and. dt > 0.0_dp .and. &
      selected_target > 0.0_dp .and. selected_target <= 1.0_dp .and. &
      (state_redist_max_order == 0 .or. state_redist_max_order == 2) .and. &
      sparse_patch_set%nvar == nvar .and. &
      sparse_patch_set%is_valid( &
        distribution, coarse_geometry, patch_set_template)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Allreduce( &
      numeric_controls, numeric_minimum, 2, MPI_DOUBLE_PRECISION, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      numeric_controls, numeric_maximum, 2, MPI_DOUBLE_PRECISION, MPI_MAX, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_controls, integer_minimum, 2, MPI_INTEGER, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_controls, integer_maximum, 2, MPI_INTEGER, MPI_MAX, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. &
        any(numeric_minimum /= numeric_maximum) .or. &
        any(integer_minimum /= integer_maximum)) return
    call MPI_Allreduce( &
      string_codes, string_minimum, size(string_codes), MPI_INTEGER, &
      MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      string_codes, string_maximum, size(string_codes), MPI_INTEGER, &
      MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. &
        any(string_minimum /= string_maximum)) return

    candidate = sparse_patch_set
    call advance_sparse_owned_reactive_eb_root_tiles_hydro_2d( &
      species, distribution, candidate, coarse_geometry, trim(solver), &
      trim(reconstruction), trim(limiter), selected_target, &
      state_redist_max_order, dt, local_tile_states, local_tile_fluxes, &
      transfers, local_ok, root_tile_advances, root_hydro_cells)
    if (.not. local_ok) return
    advances = advances + root_tile_advances
    do child = 1, distribution%child_count()
      if (allocated(coarse_support)) deallocate(coarse_support)
      if (allocated(coarse_support_temperature)) &
        deallocate(coarse_support_temperature)
      if (allocated(coarse_x_flux_support)) deallocate(coarse_x_flux_support)
      if (allocated(coarse_y_flux_support)) deallocate(coarse_y_flux_support)
      if (allocated(coarse_work)) deallocate(coarse_work)
      if (allocated(coarse_work_temperature)) &
        deallocate(coarse_work_temperature)
      if (allocated(fine_work)) deallocate(fine_work)
      if (allocated(fine_work_temperature)) &
        deallocate(fine_work_temperature)
      if (allocated(fine_x_flux)) deallocate(fine_x_flux)
      if (allocated(fine_y_flux)) deallocate(fine_y_flux)
      exterior_context = reactive_eb_patch_exterior_context_2d()
      flux_register = amr_eb_flux_register_2d()
      owner = distribution%child_owner(child)
      support_i_lower = max( &
        1, patch_set_template%children(child)%patch%coarse_i_lower - 2)
      support_i_upper = min( &
        coarse_geometry%nx, &
        patch_set_template%children(child)%patch%coarse_i_upper + 2)
      support_j_lower = max( &
        1, patch_set_template%children(child)%patch%coarse_j_lower - 2)
      support_j_upper = min( &
        coarse_geometry%ny, &
        patch_set_template%children(child)%patch%coarse_j_upper + 2)
      flux_x_i_lower = &
        patch_set_template%children(child)%patch%coarse_i_lower - 1
      flux_x_j_lower = &
        patch_set_template%children(child)%patch%coarse_j_lower
      flux_y_i_lower = &
        patch_set_template%children(child)%patch%coarse_i_lower
      flux_y_j_lower = &
        patch_set_template%children(child)%patch%coarse_j_lower - 1
      call transfer_sparse_child_state_support_2d( &
        distribution, coarse_geometry, &
        patch_set_template%children(child)%geometry, &
        patch_set_template%children(child)%patch, nvar, owner, &
        local_tile_states, exterior_context, coarse_support, &
        coarse_support_temperature, transfers, local_ok)
      if (.not. local_ok) return
      call transfer_sparse_child_coarse_flux_support_2d( &
        distribution, coarse_geometry, &
        patch_set_template%children(child)%geometry, &
        patch_set_template%children(child)%patch, nvar, owner, &
        local_tile_fluxes, coarse_x_flux_support, coarse_y_flux_support, &
        transfers, local_ok)
      if (.not. local_ok) return
      entity_ok = owner >= 0 .and. owner < distribution%nranks
      if (distribution%rank == owner .and. entity_ok) then
        call initialize_amr_eb_flux_register_2d( &
          coarse_geometry, patch_set_template%children(child)%geometry, &
          patch_set_template%children(child)%patch, nvar, flux_register, &
          entity_ok)
        if (entity_ok) call accumulate_coarse_eb_fluxes_patch_support_2d( &
          flux_register, coarse_geometry, &
          patch_set_template%children(child)%geometry, &
          patch_set_template%children(child)%patch, &
          flux_x_i_lower, flux_x_j_lower, coarse_x_flux_support, &
          flux_y_i_lower, flux_y_j_lower, coarse_y_flux_support, dt, &
          entity_ok)
        allocate(coarse_work, mold=coarse_support)
        allocate(coarse_work_temperature, mold=coarse_support_temperature)
        allocate(fine_work, mold=candidate%children(child)%state)
        allocate(fine_work_temperature, &
          mold=candidate%children(child)%temperature)
        allocate(fine_x_flux(nvar, &
          0:patch_set_template%children(child)%geometry%nx, &
          patch_set_template%children(child)%geometry%ny))
        allocate(fine_y_flux(nvar, &
          patch_set_template%children(child)%geometry%nx, &
          0:patch_set_template%children(child)%geometry%ny))
        ratio = patch_set_template%children(child)%patch%refinement_ratio
        fine_dt = dt / real(ratio, dp)
        do substep = 1, ratio
          if (.not. entity_ok) exit
          if (trim(reconstruction) == "characteristic_plm") then
            alpha = (real(substep, dp) - 0.5_dp) / real(ratio, dp)
          else
            alpha = real(substep - 1, dp) / real(ratio, dp)
          end if
          call build_reactive_eb_patch_exterior_from_context_2d( &
            species, exterior_context, coarse_geometry, &
            patch_set_template%children(child)%geometry, &
            patch_set_template%children(child)%patch, alpha, exterior, &
            entity_ok, candidate%children(child)%state, &
            candidate%children(child)%temperature)
          if (.not. entity_ok) exit
          call advance_reactive_eb_level_2d( &
            species, candidate%children(child)%state, &
            candidate%children(child)%temperature, &
            patch_set_template%children(child)%geometry, trim(solver), &
            trim(reconstruction), trim(limiter), selected_target, &
            state_redist_max_order, fine_dt, fine_work, &
            fine_work_temperature, fine_x_flux, fine_y_flux, entity_ok, &
            exterior)
          if (.not. entity_ok) exit
          advances = advances + 1
          candidate%children(child)%state = fine_work
          candidate%children(child)%temperature = fine_work_temperature
          call accumulate_fine_eb_fluxes_2d( &
            flux_register, coarse_geometry, &
            patch_set_template%children(child)%geometry, &
            patch_set_template%children(child)%patch, fine_x_flux, &
            fine_y_flux, fine_dt, entity_ok)
        end do
        if (entity_ok) call reflux_reactive_eb_state_patch_support_2d( &
          species, support_i_lower, support_j_lower, coarse_support, &
          coarse_support_temperature, coarse_geometry, &
          candidate%children(child)%state, &
          candidate%children(child)%temperature, &
          patch_set_template%children(child)%geometry, &
          patch_set_template%children(child)%patch, flux_register, &
          coarse_work, coarse_work_temperature, fine_work, &
          fine_work_temperature, entity_ok)
        if (entity_ok) then
          coarse_support = coarse_work
          coarse_support_temperature = coarse_work_temperature
          candidate%children(child)%state = fine_work
          candidate%children(child)%temperature = fine_work_temperature
        end if
      end if
      call all_ranks_accept_eb_2d( &
        distribution, entity_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) return
      call transfer_sparse_child_state_correction_2d( &
        distribution, coarse_geometry, &
        patch_set_template%children(child)%geometry, &
        patch_set_template%children(child)%patch, nvar, owner, &
        coarse_support, coarse_support_temperature, local_tile_states, &
        transfers, local_ok)
      if (.not. local_ok) return
    end do
    entity_ok = .true.
    do tile = 1, distribution%root_tile_count()
      if (.not. distribution%root_tile_is_local(tile)) cycle
      entity_ok = entity_ok .and. &
        allocated(local_tile_states(tile)%corrected_state) .and. &
        allocated(local_tile_states(tile)%corrected_temperature)
      if (.not. entity_ok) cycle
      candidate%root_tiles(tile)%state = &
        local_tile_states(tile)%corrected_state
      candidate%root_tiles(tile)%temperature = &
        local_tile_states(tile)%corrected_temperature
    end do
    call all_ranks_accept_eb_2d( &
      distribution, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call average_down_sparse_owned_reactive_eb_patch_set_2d( &
      species, distribution, candidate, coarse_geometry, &
      patch_set_template, local_ok)
    if (.not. local_ok) return
    local_ok = candidate%is_valid( &
      distribution, coarse_geometry, patch_set_template)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    sparse_patch_set = candidate
    ok = .true.
    if (present(local_level_advances)) local_level_advances = advances
    if (present(local_root_transfers)) local_root_transfers = transfers
    if (present(local_root_hydro_cells)) &
      local_root_hydro_cells = root_hydro_cells
  end subroutine advance_sparse_owned_reactive_eb_patch_set_hydro_2d

  subroutine advance_sparse_owned_reactive_eb_patch_set_transport_2d( &
      species, transport, distribution, sparse_patch_set, coarse_geometry, &
      patch_set_template, interval, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, state_redist_max_order, ok, &
      local_euler_advances, minimum_theta, &
      state_redist_target_volume_fraction, local_root_transfers, &
      local_root_transport_cells)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(inout) :: sparse_patch_set
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set_template
    real(dp), intent(in) :: interval
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    integer, intent(in) :: state_redist_max_order
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_euler_advances
    real(dp), intent(out), optional :: minimum_theta
    real(dp), intent(in), optional :: state_redist_target_volume_fraction
    integer, intent(out), optional :: local_root_transfers
    integer, intent(out), optional :: local_root_transport_cells

    type(eb_geometry_2d) :: root_tile_geometry
    type(mpi_amr_eb_sparse_patch_set_2d) :: candidate, euler, stage, start
    real(dp) :: selected_target, theta_one, theta_two
    logical :: accepted, entity_ok, global_ok, local_ok
    integer :: advances_one, advances_two, child, nvar, owner, tile
    integer :: root_transport_cells_one, root_transport_cells_two
    integer :: transfers_one, transfers_two

    ok = .false.
    if (present(local_euler_advances)) local_euler_advances = 0
    if (present(minimum_theta)) minimum_theta = 1.0_dp
    if (present(local_root_transfers)) local_root_transfers = 0
    if (present(local_root_transport_cells)) local_root_transport_cells = 0
    selected_target = 0.5_dp
    if (present(state_redist_target_volume_fraction)) &
      selected_target = state_redist_target_volume_fraction
    nvar = reactive_nvar(size(species))
    local_ok = nvar >= 1 .and. sparse_patch_set%nvar == nvar .and. &
      sparse_patch_set%is_valid( &
        distribution, coarse_geometry, patch_set_template)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    call collective_transport_preflight_2d( &
      species, transport, distribution, coarse_geometry, patch_set_template, &
      interval, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      state_redist_max_order, selected_target, local_ok)
    if (.not. local_ok) return
    if (interval <= tiny(1.0_dp) .or. .not. (viscosity_enabled .or. &
        thermal_conduction_enabled .or. species_diffusion_enabled)) then
      ok = .true.
      return
    end if

    start = sparse_patch_set
    call advance_sparse_owned_reactive_eb_patch_set_transport_euler_2d( &
      species, transport, distribution, start, coarse_geometry, &
      patch_set_template, interval, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, state_redist_max_order, &
      selected_target, stage, theta_one, local_ok, advances_one, &
      transfers_one, root_transport_cells_one)
    if (.not. local_ok) return
    call advance_sparse_owned_reactive_eb_patch_set_transport_euler_2d( &
      species, transport, distribution, stage, coarse_geometry, &
      patch_set_template, interval, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, state_redist_max_order, &
      selected_target, euler, theta_two, local_ok, advances_two, &
      transfers_two, root_transport_cells_two)
    if (.not. local_ok) return

    candidate = start
    entity_ok = .true.
    do tile = 1, distribution%root_tile_count()
      if (.not. distribution%root_tile_is_local(tile)) cycle
      call extract_eb_geometry_y_band_2d( &
        coarse_geometry, distribution%root_tiles(tile)%j_lower, &
        distribution%root_tiles(tile)%j_upper, root_tile_geometry, entity_ok)
      if (.not. entity_ok) exit
      candidate%root_tiles(tile)%state = 0.5_dp * &
        (start%root_tiles(tile)%state + euler%root_tiles(tile)%state)
      call recover_transport_temperature_2d( &
        species, candidate%root_tiles(tile)%state, &
        0.5_dp * (start%root_tiles(tile)%temperature + &
          euler%root_tiles(tile)%temperature), root_tile_geometry, &
        candidate%root_tiles(tile)%temperature, entity_ok)
      if (.not. entity_ok) exit
    end do
    call all_ranks_accept_eb_2d( &
      distribution, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    do child = 1, distribution%child_count()
      owner = distribution%child_owner(child)
      entity_ok = owner >= 0 .and. owner < distribution%nranks
      if (distribution%rank == owner .and. entity_ok) then
        candidate%children(child)%state = 0.5_dp * &
          (start%children(child)%state + euler%children(child)%state)
        call recover_transport_temperature_2d( &
          species, candidate%children(child)%state, &
          0.5_dp * (start%children(child)%temperature + &
            euler%children(child)%temperature), &
          patch_set_template%children(child)%geometry, &
          candidate%children(child)%temperature, entity_ok)
      end if
      call all_ranks_accept_eb_2d( &
        distribution, entity_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) return
    end do
    call average_down_sparse_owned_reactive_eb_patch_set_2d( &
      species, distribution, candidate, coarse_geometry, &
      patch_set_template, local_ok)
    if (.not. local_ok) return
    local_ok = candidate%is_valid( &
      distribution, coarse_geometry, patch_set_template)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    sparse_patch_set = candidate
    ok = .true.
    if (present(local_euler_advances)) &
      local_euler_advances = advances_one + advances_two
    if (present(minimum_theta)) minimum_theta = min(theta_one, theta_two)
    if (present(local_root_transfers)) &
      local_root_transfers = transfers_one + transfers_two
    if (present(local_root_transport_cells)) &
      local_root_transport_cells = &
        root_transport_cells_one + root_transport_cells_two
  end subroutine advance_sparse_owned_reactive_eb_patch_set_transport_2d

  subroutine advance_sparse_owned_reactive_eb_patch_set_transport_euler_2d( &
      species, transport, distribution, sparse_patch_set, coarse_geometry, &
      patch_set_template, dt, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, state_redist_max_order, &
      target_volume_fraction, new_sparse_patch_set, minimum_theta, ok, &
      local_euler_advances, local_root_transfers, &
      local_root_transport_cells)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(in) :: sparse_patch_set
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set_template
    real(dp), intent(in) :: dt
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    integer, intent(in) :: state_redist_max_order
    real(dp), intent(in) :: target_volume_fraction
    type(mpi_amr_eb_sparse_patch_set_2d), intent(out) :: new_sparse_patch_set
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    integer, intent(out) :: local_euler_advances
    integer, intent(out) :: local_root_transfers
    integer, intent(out) :: local_root_transport_cells

    type(amr_eb_flux_register_2d) :: flux_register
    type(reactive_eb_exterior_state_2d) :: exterior
    type(mpi_amr_eb_sparse_patch_set_2d) :: candidate
    type(mpi_amr_eb_root_tile_transport_state_2d), allocatable :: &
      local_tile_states(:)
    type(mpi_amr_eb_root_tile_transport_flux_2d), allocatable :: &
      local_tile_fluxes(:)
    real(dp), allocatable :: coarse_work(:, :, :)
    real(dp), allocatable :: coarse_work_temperature(:, :)
    real(dp), allocatable :: coarse_x_flux_support(:, :, :)
    real(dp), allocatable :: coarse_y_flux_support(:, :, :)
    real(dp), allocatable :: fine_rhs(:, :, :), fine_work(:, :, :)
    real(dp), allocatable :: fine_work_temperature(:, :)
    real(dp), allocatable :: coarse_support(:, :, :)
    real(dp), allocatable :: coarse_support_temperature(:, :)
    type(reactive_eb_patch_exterior_context_2d) :: exterior_context
    real(dp), allocatable :: fine_x_flux(:, :, :), fine_y_flux(:, :, :)
    real(dp), allocatable :: boundary_change(:), integral_before(:)
    real(dp) :: alpha, coarse_theta, fine_dt, fine_theta, local_theta
    logical :: accepted, cut_interface, entity_ok, global_ok, local_ok
    integer :: advances, child, ierr, nvar, owner, ratio
    integer :: flux_x_i_lower, flux_x_j_lower
    integer :: flux_y_i_lower, flux_y_j_lower
    integer :: root_tile_advances, root_transport_cells
    integer :: substep, support_i_lower, support_i_upper, tile
    integer :: support_j_lower, support_j_upper, transfers

    new_sparse_patch_set = sparse_patch_set
    minimum_theta = 1.0_dp
    local_euler_advances = 0
    local_root_transfers = 0
    local_root_transport_cells = 0
    ok = .false.
    advances = 0
    transfers = 0
    local_theta = 1.0_dp
    nvar = reactive_nvar(size(species))
    candidate = sparse_patch_set

    allocate(integral_before(nvar))
    call composite_sparse_owned_reactive_eb_patch_set_integral_2d( &
      distribution, candidate, coarse_geometry, patch_set_template, &
      integral_before, local_ok)
    if (.not. local_ok) return

    call advance_sparse_owned_reactive_eb_root_tiles_transport_euler_2d( &
      species, transport, distribution, candidate, coarse_geometry, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, state_redist_max_order, local_tile_states, &
      local_tile_fluxes, coarse_theta, transfers, local_ok, &
      root_tile_advances, root_transport_cells)
    if (.not. local_ok) return
    advances = advances + root_tile_advances
    local_theta = min(local_theta, coarse_theta)
    cut_interface = .false.
    do child = 1, distribution%child_count()
      if (allocated(coarse_work)) deallocate(coarse_work)
      if (allocated(coarse_work_temperature)) &
        deallocate(coarse_work_temperature)
      if (allocated(coarse_support)) deallocate(coarse_support)
      if (allocated(coarse_support_temperature)) &
        deallocate(coarse_support_temperature)
      if (allocated(coarse_x_flux_support)) deallocate(coarse_x_flux_support)
      if (allocated(coarse_y_flux_support)) deallocate(coarse_y_flux_support)
      if (allocated(fine_rhs)) deallocate(fine_rhs)
      if (allocated(fine_work)) deallocate(fine_work)
      if (allocated(fine_work_temperature)) &
        deallocate(fine_work_temperature)
      if (allocated(fine_x_flux)) deallocate(fine_x_flux)
      if (allocated(fine_y_flux)) deallocate(fine_y_flux)
      exterior_context = reactive_eb_patch_exterior_context_2d()
      flux_register = amr_eb_flux_register_2d()

      cut_interface = cut_interface .or. .not. level_two_interface_is_regular( &
        patch_set_template%children(child)%geometry)
      owner = distribution%child_owner(child)
      support_i_lower = max( &
        1, patch_set_template%children(child)%patch%coarse_i_lower - 2)
      support_i_upper = min( &
        coarse_geometry%nx, &
        patch_set_template%children(child)%patch%coarse_i_upper + 2)
      support_j_lower = max( &
        1, patch_set_template%children(child)%patch%coarse_j_lower - 2)
      support_j_upper = min( &
        coarse_geometry%ny, &
        patch_set_template%children(child)%patch%coarse_j_upper + 2)
      flux_x_i_lower = &
        patch_set_template%children(child)%patch%coarse_i_lower - 1
      flux_x_j_lower = &
        patch_set_template%children(child)%patch%coarse_j_lower
      flux_y_i_lower = &
        patch_set_template%children(child)%patch%coarse_i_lower
      flux_y_j_lower = &
        patch_set_template%children(child)%patch%coarse_j_lower - 1
      call transfer_sparse_child_state_support_2d( &
        distribution, coarse_geometry, &
        patch_set_template%children(child)%geometry, &
        patch_set_template%children(child)%patch, nvar, owner, &
        local_tile_states, exterior_context, coarse_support, &
        coarse_support_temperature, transfers, local_ok)
      if (.not. local_ok) return

      call transfer_sparse_child_coarse_flux_support_2d( &
        distribution, coarse_geometry, &
        patch_set_template%children(child)%geometry, &
        patch_set_template%children(child)%patch, nvar, owner, &
        local_tile_fluxes, coarse_x_flux_support, coarse_y_flux_support, &
        transfers, local_ok)
      if (.not. local_ok) return

      entity_ok = owner >= 0 .and. owner < distribution%nranks
      if (distribution%rank == owner .and. entity_ok) then
        call initialize_amr_eb_flux_register_2d( &
          coarse_geometry, patch_set_template%children(child)%geometry, &
          patch_set_template%children(child)%patch, nvar, flux_register, &
          entity_ok)
        if (entity_ok) call accumulate_coarse_eb_fluxes_patch_support_2d( &
          flux_register, coarse_geometry, &
          patch_set_template%children(child)%geometry, &
          patch_set_template%children(child)%patch, &
          flux_x_i_lower, flux_x_j_lower, coarse_x_flux_support, &
          flux_y_i_lower, flux_y_j_lower, coarse_y_flux_support, dt, &
          entity_ok)
      end if
      call all_ranks_accept_eb_2d( &
        distribution, entity_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) return
      if (distribution%rank == owner .and. entity_ok) then
        allocate(fine_rhs, mold=candidate%children(child)%state)
        allocate(fine_work, mold=candidate%children(child)%state)
        allocate(fine_work_temperature, &
          mold=candidate%children(child)%temperature)
        allocate(fine_x_flux(nvar, &
          0:patch_set_template%children(child)%geometry%nx, &
          patch_set_template%children(child)%geometry%ny))
        allocate(fine_y_flux(nvar, &
          patch_set_template%children(child)%geometry%nx, &
          0:patch_set_template%children(child)%geometry%ny))
        ratio = patch_set_template%children(child)%patch%refinement_ratio
        fine_dt = dt / real(ratio, dp)
        do substep = 1, ratio
          if (.not. entity_ok) exit
          alpha = real(substep - 1, dp) / real(ratio, dp)
          call build_reactive_eb_patch_exterior_from_context_2d( &
            species, exterior_context, coarse_geometry, &
            patch_set_template%children(child)%geometry, &
            patch_set_template%children(child)%patch, alpha, exterior, &
            entity_ok, candidate%children(child)%state, &
            candidate%children(child)%temperature)
          if (.not. entity_ok) exit
          call reactive_eb_transport_fluxes_rhs_2d( &
            species, transport, candidate%children(child)%state, &
            candidate%children(child)%temperature, &
            patch_set_template%children(child)%geometry, fine_dt, &
            viscosity_enabled, thermal_conduction_enabled, &
            species_diffusion_enabled, barodiffusion_enabled, boundaries, &
            fine_rhs, fine_x_flux, fine_y_flux, fine_theta, entity_ok, &
            exterior)
          if (.not. entity_ok) exit
          call advance_reactive_eb_state_redistributed_2d( &
            species, candidate%children(child)%state, &
            candidate%children(child)%temperature, &
            patch_set_template%children(child)%geometry, fine_rhs, fine_dt, &
            fine_work, fine_work_temperature, entity_ok, &
            target_volume_fraction, state_redist_max_order)
          if (.not. entity_ok) exit
          advances = advances + 1
          local_theta = min(local_theta, fine_theta)
          candidate%children(child)%state = fine_work
          candidate%children(child)%temperature = fine_work_temperature
          call accumulate_fine_eb_fluxes_2d( &
            flux_register, coarse_geometry, &
            patch_set_template%children(child)%geometry, &
            patch_set_template%children(child)%patch, fine_x_flux, &
            fine_y_flux, fine_dt, entity_ok)
        end do
        if (entity_ok) then
          allocate(coarse_work( &
            nvar, support_i_lower:support_i_upper, &
            support_j_lower:support_j_upper))
          allocate(coarse_work_temperature( &
            support_i_lower:support_i_upper, &
            support_j_lower:support_j_upper))
          call reflux_reactive_eb_state_patch_support_2d( &
            species, support_i_lower, support_j_lower, coarse_support, &
            coarse_support_temperature, coarse_geometry, &
            candidate%children(child)%state, &
            candidate%children(child)%temperature, &
            patch_set_template%children(child)%geometry, &
            patch_set_template%children(child)%patch, flux_register, &
            coarse_work, coarse_work_temperature, fine_work, &
            fine_work_temperature, entity_ok)
        end if
        if (entity_ok) then
          coarse_support = coarse_work
          coarse_support_temperature = coarse_work_temperature
          candidate%children(child)%state = fine_work
          candidate%children(child)%temperature = fine_work_temperature
        end if
      end if
      call all_ranks_accept_eb_2d( &
        distribution, entity_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) return
      call transfer_sparse_child_state_correction_2d( &
        distribution, coarse_geometry, &
        patch_set_template%children(child)%geometry, &
        patch_set_template%children(child)%patch, nvar, owner, &
        coarse_support, coarse_support_temperature, local_tile_states, &
        transfers, local_ok)
      if (.not. local_ok) return
    end do
    entity_ok = .true.
    do tile = 1, distribution%root_tile_count()
      if (.not. distribution%root_tile_is_local(tile)) cycle
      entity_ok = entity_ok .and. &
        allocated(local_tile_states(tile)%corrected_state) .and. &
        allocated(local_tile_states(tile)%corrected_temperature)
      if (.not. entity_ok) cycle
      candidate%root_tiles(tile)%state = &
        local_tile_states(tile)%corrected_state
      candidate%root_tiles(tile)%temperature = &
        local_tile_states(tile)%corrected_temperature
    end do
    call all_ranks_accept_eb_2d( &
      distribution, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call average_down_sparse_owned_reactive_eb_patch_set_2d( &
      species, distribution, candidate, coarse_geometry, &
      patch_set_template, local_ok)
    if (.not. local_ok) return
    if (cut_interface) then
      call reduce_sparse_root_transport_boundary_change_2d( &
        distribution, coarse_geometry, nvar, local_tile_fluxes, dt, &
        boundary_change, local_ok)
      if (.not. local_ok) return
      call close_sparse_cut_patch_set_conservation_2d( &
        species, integral_before, distribution, candidate, coarse_geometry, &
        patch_set_template, boundary_change, local_ok)
      if (.not. local_ok) return
    end if
    call MPI_Allreduce( &
      local_theta, minimum_theta, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    local_ok = candidate%is_valid( &
      distribution, coarse_geometry, patch_set_template) .and. &
      ieee_is_finite(minimum_theta)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    new_sparse_patch_set = candidate
    local_euler_advances = advances
    local_root_transfers = transfers
    local_root_transport_cells = root_transport_cells
    ok = .true.
  end subroutine advance_sparse_owned_reactive_eb_patch_set_transport_euler_2d

  subroutine advance_sparse_owned_reactive_eb_patch_set_strang_2d( &
      species, reactions, transport, distribution, sparse_patch_set, &
      coarse_geometry, patch_set_template, solver, reconstruction, limiter, &
      state_redist_max_order, dt, rtol, atol, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, ok, local_chemistry_advances, &
      local_hydro_advances, local_transport_euler_advances, &
      minimum_transport_theta, state_redist_target_volume_fraction)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(inout) :: sparse_patch_set
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set_template
    character(len=*), intent(in) :: solver, reconstruction, limiter
    integer, intent(in) :: state_redist_max_order
    real(dp), intent(in) :: dt, rtol, atol
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_chemistry_advances
    integer, intent(out), optional :: local_hydro_advances
    integer, intent(out), optional :: local_transport_euler_advances
    real(dp), intent(out), optional :: minimum_transport_theta
    real(dp), intent(in), optional :: state_redist_target_volume_fraction

    type(mpi_amr_eb_sparse_patch_set_2d) :: candidate
    real(dp) :: selected_target, theta, theta_one, theta_two
    logical :: local_ok
    integer :: chemistry_one, chemistry_two, hydro_advances
    integer :: transport_one, transport_two

    ok = .false.
    if (present(local_chemistry_advances)) local_chemistry_advances = 0
    if (present(local_hydro_advances)) local_hydro_advances = 0
    if (present(local_transport_euler_advances)) &
      local_transport_euler_advances = 0
    if (present(minimum_transport_theta)) minimum_transport_theta = 1.0_dp
    selected_target = 0.5_dp
    if (present(state_redist_target_volume_fraction)) &
      selected_target = state_redist_target_volume_fraction
    candidate = sparse_patch_set

    call advance_sparse_owned_reactive_eb_patch_set_chemistry_2d( &
      species, reactions, 0.5_dp * dt, rtol, atol, distribution, candidate, &
      coarse_geometry, patch_set_template, local_ok, chemistry_one)
    if (.not. local_ok) return

    call advance_sparse_owned_reactive_eb_patch_set_transport_2d( &
      species, transport, distribution, candidate, coarse_geometry, &
      patch_set_template, 0.5_dp * dt, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, state_redist_max_order, local_ok, &
      transport_one, theta_one, selected_target)
    if (.not. local_ok) return
    call advance_sparse_owned_reactive_eb_patch_set_hydro_2d( &
      species, distribution, candidate, coarse_geometry, patch_set_template, &
      solver, reconstruction, limiter, state_redist_max_order, dt, local_ok, &
      hydro_advances, selected_target)
    if (.not. local_ok) return
    call advance_sparse_owned_reactive_eb_patch_set_transport_2d( &
      species, transport, distribution, candidate, coarse_geometry, &
      patch_set_template, 0.5_dp * dt, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, state_redist_max_order, local_ok, &
      transport_two, theta_two, selected_target)
    if (.not. local_ok) return

    call advance_sparse_owned_reactive_eb_patch_set_chemistry_2d( &
      species, reactions, 0.5_dp * dt, rtol, atol, distribution, candidate, &
      coarse_geometry, patch_set_template, local_ok, chemistry_two)
    if (.not. local_ok) return

    theta = min(theta_one, theta_two)
    sparse_patch_set = candidate
    ok = .true.
    if (present(local_chemistry_advances)) &
      local_chemistry_advances = chemistry_one + chemistry_two
    if (present(local_hydro_advances)) &
      local_hydro_advances = hydro_advances
    if (present(local_transport_euler_advances)) &
      local_transport_euler_advances = transport_one + transport_two
    if (present(minimum_transport_theta)) minimum_transport_theta = theta
  end subroutine advance_sparse_owned_reactive_eb_patch_set_strang_2d

  subroutine advance_sparse_owned_reactive_eb_patch_set_to_time_2d( &
      species, reactions, transport, distribution, sparse_patch_set, &
      coarse_geometry, patch_set_template, solver, reconstruction, limiter, &
      state_redist_max_order, time, final_time, steps, maximum_steps, &
      hydro_cfl, transport_cfl, rtol, atol, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, ok, minimum_dt, advanced_steps, &
      local_chemistry_advances, local_hydro_advances, &
      local_transport_euler_advances, minimum_transport_theta, &
      state_redist_target_volume_fraction, local_timestep_root_transfers, &
      regrid_evaluations, regrids, regrid_interval, regrid_criteria, &
      refinement_ratio, geometry_builder, local_regrid_root_transfers, &
      local_regrid_restriction_transfers, &
      local_regrid_prolongation_transfers, local_regrid_overlap_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(mpi_amr_eb_patch_distribution_2d), intent(inout) :: distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(inout) :: sparse_patch_set
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(inout) :: patch_set_template
    character(len=*), intent(in) :: solver, reconstruction, limiter
    integer, intent(in) :: state_redist_max_order
    real(dp), intent(inout) :: time
    real(dp), intent(in) :: final_time
    integer, intent(inout) :: steps
    integer, intent(in) :: maximum_steps
    real(dp), intent(in) :: hydro_cfl, transport_cfl, rtol, atol
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    logical, intent(out) :: ok
    real(dp), intent(out) :: minimum_dt
    integer, intent(out), optional :: advanced_steps
    integer, intent(out), optional :: local_chemistry_advances
    integer, intent(out), optional :: local_hydro_advances
    integer, intent(out), optional :: local_transport_euler_advances
    real(dp), intent(out), optional :: minimum_transport_theta
    real(dp), intent(in), optional :: state_redist_target_volume_fraction
    integer, intent(out), optional :: local_timestep_root_transfers
    integer, intent(inout), optional :: regrid_evaluations, regrids
    integer, intent(in), optional :: regrid_interval, refinement_ratio
    type(amr_eb_tagging_criteria_2d), intent(in), optional :: regrid_criteria
    procedure(sparse_eb_geometry_builder_2d), optional :: geometry_builder
    integer, intent(out), optional :: local_regrid_root_transfers
    integer, intent(out), optional :: local_regrid_restriction_transfers
    integer, intent(out), optional :: local_regrid_prolongation_transfers
    integer, intent(out), optional :: local_regrid_overlap_transfers

    type(amr_eb_tagging_criteria_2d) :: selected_regrid_criteria
    type(mpi_amr_eb_patch_distribution_2d) :: candidate_distribution
    type(mpi_amr_eb_sparse_patch_set_2d) :: candidate_sparse_patch_set
    type(reactive_eb_patch_set_2d) :: candidate_patch_set_template
    real(dp) :: dt, numeric_controls(10), numeric_maximum(10)
    real(dp) :: numeric_minimum(10), remaining, selected_target
    real(dp) :: step_theta, time_tolerance
    logical :: accepted, global_ok, local_ok, regrid_enabled
    logical :: regrid_requested, scheduled_regrid, step_changed
    integer :: chemistry_advances, hydro_advances, ierr
    integer :: integer_controls(11), integer_maximum(11)
    integer :: integer_minimum(11)
    integer :: regrid_transfers, selected_regrid_evaluations
    integer :: regrid_overlap_transfers, regrid_prolongation_transfers
    integer :: regrid_restriction_transfers
    integer :: selected_regrid_interval, selected_regrids
    integer :: selected_refinement_ratio, step_regrid_transfers
    integer :: step_regrid_overlap_transfers
    integer :: step_regrid_prolongation_transfers
    integer :: step_regrid_restriction_transfers
    integer :: step_chemistry, step_hydro, step_timestep_transfers
    integer :: step_transport, timestep_transfers, transport_advances

    ok = .false.
    minimum_dt = 0.0_dp
    chemistry_advances = 0
    hydro_advances = 0
    transport_advances = 0
    timestep_transfers = 0
    regrid_transfers = 0
    regrid_restriction_transfers = 0
    regrid_prolongation_transfers = 0
    regrid_overlap_transfers = 0
    if (present(advanced_steps)) advanced_steps = 0
    if (present(local_chemistry_advances)) local_chemistry_advances = 0
    if (present(local_hydro_advances)) local_hydro_advances = 0
    if (present(local_transport_euler_advances)) &
      local_transport_euler_advances = 0
    if (present(minimum_transport_theta)) minimum_transport_theta = 1.0_dp
    if (present(local_timestep_root_transfers)) &
      local_timestep_root_transfers = 0
    if (present(local_regrid_root_transfers)) &
      local_regrid_root_transfers = 0
    if (present(local_regrid_restriction_transfers)) &
      local_regrid_restriction_transfers = 0
    if (present(local_regrid_prolongation_transfers)) &
      local_regrid_prolongation_transfers = 0
    if (present(local_regrid_overlap_transfers)) &
      local_regrid_overlap_transfers = 0
    selected_target = 0.5_dp
    if (present(state_redist_target_volume_fraction)) &
      selected_target = state_redist_target_volume_fraction
    selected_regrid_evaluations = 0
    selected_regrids = 0
    selected_regrid_interval = 1
    selected_refinement_ratio = 2
    selected_regrid_criteria = amr_eb_tagging_criteria_2d()
    if (present(regrid_evaluations)) &
      selected_regrid_evaluations = regrid_evaluations
    if (present(regrids)) selected_regrids = regrids
    if (present(regrid_interval)) selected_regrid_interval = regrid_interval
    if (present(refinement_ratio)) &
      selected_refinement_ratio = refinement_ratio
    if (present(regrid_criteria)) &
      selected_regrid_criteria = regrid_criteria
    regrid_requested = present(regrid_evaluations) .or. present(regrids) .or. &
      present(regrid_interval) .or. present(regrid_criteria) .or. &
      present(refinement_ratio) .or. present(geometry_builder) .or. &
      present(local_regrid_root_transfers) .or. &
      present(local_regrid_restriction_transfers) .or. &
      present(local_regrid_prolongation_transfers) .or. &
      present(local_regrid_overlap_transfers)
    regrid_enabled = &
      present(regrid_evaluations) .and. present(regrids) .and. &
      present(regrid_interval) .and. present(regrid_criteria) .and. &
      present(refinement_ratio) .and. present(geometry_builder)
    numeric_controls = [ &
      time, final_time, hydro_cfl, transport_cfl, rtol, atol, selected_target, &
      selected_regrid_criteria%relative_gradient_threshold, &
      selected_regrid_criteria%absolute_gradient_threshold, &
      selected_regrid_criteria%scale_floor]
    integer_controls = [ &
      steps, maximum_steps, merge(1, 0, regrid_enabled), &
      selected_regrid_interval, selected_refinement_ratio, &
      selected_regrid_evaluations, selected_regrids, &
      selected_regrid_criteria%buffer_cells, &
      selected_regrid_criteria%minimum_patch_cells_x, &
      selected_regrid_criteria%minimum_patch_cells_y, &
      selected_regrid_criteria%maximum_patch_gap_cells]
    time_tolerance = 16.0_dp * epsilon(1.0_dp) * &
      max(tiny(1.0_dp), abs(final_time))
    local_ok = size(species) >= 1 .and. &
      all(ieee_is_finite(numeric_controls)) .and. &
      final_time >= time - time_tolerance .and. &
      hydro_cfl > 0.0_dp .and. hydro_cfl <= 1.0_dp .and. &
      transport_cfl > 0.0_dp .and. transport_cfl <= 0.5_dp .and. &
      rtol > 0.0_dp .and. atol > 0.0_dp .and. &
      selected_target > 0.0_dp .and. selected_target <= 1.0_dp .and. &
      steps >= 0 .and. maximum_steps >= steps .and. &
      (.not. regrid_requested .or. regrid_enabled) .and. &
      (.not. regrid_enabled .or. &
        (selected_regrid_interval >= 1 .and. &
        selected_refinement_ratio >= 2 .and. &
        selected_regrid_evaluations >= 0 .and. selected_regrids >= 0 .and. &
        selected_regrids <= selected_regrid_evaluations .and. &
        selected_regrid_criteria%is_valid( &
          coarse_geometry%nx, coarse_geometry%ny))) .and. &
      sparse_patch_set%nvar == reactive_nvar(size(species)) .and. &
      sparse_patch_set%is_valid( &
        distribution, coarse_geometry, patch_set_template)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Allreduce( &
      numeric_controls, numeric_minimum, 10, MPI_DOUBLE_PRECISION, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      numeric_controls, numeric_maximum, 10, MPI_DOUBLE_PRECISION, MPI_MAX, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_controls, integer_minimum, 11, MPI_INTEGER, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_controls, integer_maximum, 11, MPI_INTEGER, MPI_MAX, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. &
        any(numeric_minimum /= numeric_maximum) .or. &
        any(integer_minimum /= integer_maximum)) return

    do
      remaining = final_time - time
      if (remaining <= time_tolerance) exit
      if (steps >= maximum_steps) return
      call compute_sparse_owned_reactive_eb_patch_set_timestep_2d( &
        species, transport, distribution, sparse_patch_set, &
        coarse_geometry, patch_set_template, hydro_cfl, transport_cfl, &
        viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, boundaries, dt, &
        local_ok, step_timestep_transfers)
      if (.not. local_ok) return
      dt = min(dt, remaining)
      if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp) return
      scheduled_regrid = regrid_enabled .and. &
        modulo(steps + 1, selected_regrid_interval) == 0
      step_regrid_transfers = 0
      step_regrid_restriction_transfers = 0
      step_regrid_prolongation_transfers = 0
      step_regrid_overlap_transfers = 0
      step_changed = .false.
      if (scheduled_regrid) then
        candidate_distribution = distribution
        candidate_sparse_patch_set = sparse_patch_set
        candidate_patch_set_template = patch_set_template
        call advance_sparse_owned_reactive_eb_patch_set_strang_2d( &
          species, reactions, transport, candidate_distribution, &
          candidate_sparse_patch_set, coarse_geometry, &
          candidate_patch_set_template, solver, reconstruction, limiter, &
          state_redist_max_order, dt, rtol, atol, viscosity_enabled, &
          thermal_conduction_enabled, species_diffusion_enabled, &
          barodiffusion_enabled, boundaries, local_ok, step_chemistry, &
          step_hydro, step_transport, step_theta, selected_target)
        if (.not. local_ok) return
        call regrid_tagged_sparse_owned_reactive_eb_patch_set_2d( &
          species, candidate_distribution, candidate_sparse_patch_set, &
          coarse_geometry, candidate_patch_set_template, &
          selected_regrid_criteria, selected_refinement_ratio, &
          geometry_builder, local_ok, step_changed, step_regrid_transfers, &
          step_regrid_restriction_transfers, &
          step_regrid_prolongation_transfers, &
          step_regrid_overlap_transfers)
        if (.not. local_ok) return
        distribution = candidate_distribution
        sparse_patch_set = candidate_sparse_patch_set
        patch_set_template = candidate_patch_set_template
      else
        call advance_sparse_owned_reactive_eb_patch_set_strang_2d( &
          species, reactions, transport, distribution, sparse_patch_set, &
          coarse_geometry, patch_set_template, solver, reconstruction, &
          limiter, state_redist_max_order, dt, rtol, atol, &
          viscosity_enabled, thermal_conduction_enabled, &
          species_diffusion_enabled, barodiffusion_enabled, boundaries, &
          local_ok, step_chemistry, step_hydro, step_transport, step_theta, &
          selected_target)
        if (.not. local_ok) return
      end if

      time = time + dt
      steps = steps + 1
      chemistry_advances = chemistry_advances + step_chemistry
      hydro_advances = hydro_advances + step_hydro
      transport_advances = transport_advances + step_transport
      timestep_transfers = timestep_transfers + step_timestep_transfers
      regrid_transfers = regrid_transfers + step_regrid_transfers
      regrid_restriction_transfers = regrid_restriction_transfers + &
        step_regrid_restriction_transfers
      regrid_prolongation_transfers = regrid_prolongation_transfers + &
        step_regrid_prolongation_transfers
      regrid_overlap_transfers = regrid_overlap_transfers + &
        step_regrid_overlap_transfers
      if (scheduled_regrid) then
        selected_regrid_evaluations = selected_regrid_evaluations + 1
        if (step_changed) selected_regrids = selected_regrids + 1
      end if
      if (minimum_dt == 0.0_dp) then
        minimum_dt = dt
      else
        minimum_dt = min(minimum_dt, dt)
      end if
      if (present(advanced_steps)) advanced_steps = advanced_steps + 1
      if (present(local_chemistry_advances)) &
        local_chemistry_advances = chemistry_advances
      if (present(local_hydro_advances)) &
        local_hydro_advances = hydro_advances
      if (present(local_transport_euler_advances)) &
        local_transport_euler_advances = transport_advances
      if (present(minimum_transport_theta)) &
        minimum_transport_theta = min(minimum_transport_theta, step_theta)
      if (present(local_timestep_root_transfers)) &
        local_timestep_root_transfers = timestep_transfers
      if (present(local_regrid_root_transfers)) &
        local_regrid_root_transfers = regrid_transfers
      if (present(local_regrid_restriction_transfers)) &
        local_regrid_restriction_transfers = regrid_restriction_transfers
      if (present(local_regrid_prolongation_transfers)) &
        local_regrid_prolongation_transfers = regrid_prolongation_transfers
      if (present(local_regrid_overlap_transfers)) &
        local_regrid_overlap_transfers = regrid_overlap_transfers
      if (present(regrid_evaluations)) &
        regrid_evaluations = selected_regrid_evaluations
      if (present(regrids)) regrids = selected_regrids
    end do

    time = final_time
    ok = .true.
  end subroutine advance_sparse_owned_reactive_eb_patch_set_to_time_2d

  subroutine initialize_mpi_amr_eb_patch_distribution_2d( &
      coarse_geometry, patch_set, comm, distribution, ok, &
      subcycle_exponent)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set
    type(MPI_Comm), intent(in) :: comm
    type(mpi_amr_eb_patch_distribution_2d), intent(out) :: distribution
    logical, intent(out) :: ok
    integer, intent(in), optional :: subcycle_exponent

    integer(int64) :: cells, level_scale, work_count
    logical :: local_ok
    integer :: base_rows, child, exponent, exponent_max, exponent_min
    integer :: extra_rows, ierr, j_lower, nranks, owner, power, rank
    integer :: ratio, rows, tile, tile_count

    distribution%comm = comm
    ok = .false.
    call MPI_Comm_rank(comm, rank, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Comm_size(comm, nranks, ierr)
    if (ierr /= MPI_SUCCESS .or. nranks < 1) return
    distribution%rank = rank
    distribution%nranks = nranks
    exponent = 0
    if (present(subcycle_exponent)) exponent = subcycle_exponent
    call MPI_Allreduce( &
      exponent, exponent_min, 1, MPI_INTEGER, MPI_MIN, comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      exponent, exponent_max, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
    if (ierr /= MPI_SUCCESS .or. exponent_min /= exponent_max .or. &
        exponent < 0 .or. exponent > 2) return
    distribution%subcycle_exponent = exponent
    call replicated_reactive_eb_patch_set_matches_2d( &
      coarse_geometry, patch_set, comm, local_ok)
    if (.not. local_ok) return

    tile_count = min(coarse_geometry%ny, nranks)
    allocate(distribution%root_tiles(tile_count))
    allocate(distribution%child_owners(patch_set%patch_count()))
    allocate(distribution%child_cell_counts(patch_set%patch_count()))
    allocate(distribution%child_work_counts(patch_set%patch_count()))
    allocate(distribution%rank_cell_counts(nranks))
    allocate(distribution%rank_entity_counts(nranks))
    allocate(distribution%rank_work_counts(nranks))
    distribution%rank_cell_counts = 0
    distribution%rank_entity_counts = 0
    distribution%rank_work_counts = 0_int64

    base_rows = coarse_geometry%ny / tile_count
    extra_rows = modulo(coarse_geometry%ny, tile_count)
    j_lower = 1
    do tile = 1, tile_count
      rows = base_rows
      if (tile <= extra_rows) rows = rows + 1
      cells = int(coarse_geometry%nx, int64) * int(rows, int64)
      if (cells > int(huge(1), int64)) return
      owner = minloc(distribution%rank_work_counts, dim=1)
      distribution%root_tiles(tile)%owner = owner - 1
      distribution%root_tiles(tile)%i_lower = 1
      distribution%root_tiles(tile)%i_upper = coarse_geometry%nx
      distribution%root_tiles(tile)%j_lower = j_lower
      distribution%root_tiles(tile)%j_upper = j_lower + rows - 1
      distribution%root_tiles(tile)%cell_count = int(cells)
      distribution%root_tiles(tile)%work_count = cells
      distribution%rank_cell_counts(owner) = &
        distribution%rank_cell_counts(owner) + int(cells)
      distribution%rank_entity_counts(owner) = &
        distribution%rank_entity_counts(owner) + 1
      distribution%rank_work_counts(owner) = &
        distribution%rank_work_counts(owner) + cells
      j_lower = j_lower + rows
    end do

    do child = 1, patch_set%patch_count()
      cells = int(patch_set%children(child)%geometry%nx, int64) * &
        int(patch_set%children(child)%geometry%ny, int64)
      if (cells > int(huge(1), int64)) return
      ratio = patch_set%children(child)%patch%refinement_ratio
      level_scale = 1_int64
      do power = 1, exponent
        if (ratio < 1 .or. &
            level_scale > huge(level_scale) / int(ratio, int64)) return
        level_scale = level_scale * int(ratio, int64)
      end do
      if (cells > huge(work_count) / level_scale) return
      work_count = cells * level_scale
      owner = minloc(distribution%rank_work_counts, dim=1)
      distribution%child_owners(child) = owner - 1
      distribution%child_cell_counts(child) = int(cells)
      distribution%child_work_counts(child) = work_count
      distribution%rank_cell_counts(owner) = &
        distribution%rank_cell_counts(owner) + int(cells)
      distribution%rank_entity_counts(owner) = &
        distribution%rank_entity_counts(owner) + 1
      distribution%rank_work_counts(owner) = &
        distribution%rank_work_counts(owner) + work_count
    end do
    ok = distribution%is_valid(coarse_geometry, patch_set)
  end subroutine initialize_mpi_amr_eb_patch_distribution_2d

  pure logical function mpi_amr_eb_distribution_matches_patch_set_2d( &
      distribution, coarse_geometry, patch_set) result(matches)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set

    matches = distribution%is_valid(coarse_geometry, patch_set) .and. &
      distribution%root_tile_count() == &
        min(coarse_geometry%ny, distribution%nranks) .and. &
      distribution%child_count() == patch_set%patch_count()
  end function mpi_amr_eb_distribution_matches_patch_set_2d

  subroutine synchronize_owned_reactive_eb_patch_set_2d( &
      distribution, species_count, coarse_state, coarse_temperature, &
      coarse_geometry, patch_set, synchronized_coarse_state, &
      synchronized_coarse_temperature, synchronized_patch_set, ok)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    integer, intent(in) :: species_count
    real(dp), intent(in) :: coarse_state(:, :, :), coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set
    real(dp), intent(out) :: synchronized_coarse_state(:, :, :)
    real(dp), intent(out) :: synchronized_coarse_temperature(:, :)
    type(reactive_eb_patch_set_2d), intent(out) :: synchronized_patch_set
    logical, intent(out) :: ok

    type(reactive_eb_patch_set_2d) :: candidate_set
    real(dp), allocatable :: candidate_state(:, :, :)
    real(dp), allocatable :: candidate_temperature(:, :)
    logical :: global_ok, local_ok
    integer :: child, ierr, nvar, owner, tile

    synchronized_coarse_state = coarse_state
    synchronized_coarse_temperature = coarse_temperature
    synchronized_patch_set = patch_set
    ok = .false.
    nvar = reactive_nvar(species_count)
    local_ok = nvar >= 1 .and. &
      all(shape(coarse_state) == &
        [nvar, coarse_geometry%nx, coarse_geometry%ny]) .and. &
      all(shape(coarse_temperature) == &
        [coarse_geometry%nx, coarse_geometry%ny]) .and. &
      all(shape(synchronized_coarse_state) == shape(coarse_state)) .and. &
      all(shape(synchronized_coarse_temperature) == &
        shape(coarse_temperature)) .and. &
      patch_set%is_valid(coarse_geometry, nvar) .and. &
      distribution%is_valid(coarse_geometry, patch_set)
    call MPI_Allreduce( &
      local_ok, global_ok, 1, MPI_LOGICAL, MPI_LAND, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. .not. global_ok) return
    call replicated_reactive_eb_patch_set_matches_2d( &
      coarse_geometry, patch_set, distribution%comm, local_ok)
    if (.not. local_ok) return

    allocate(candidate_state, source=coarse_state)
    allocate(candidate_temperature, source=coarse_temperature)
    candidate_set = patch_set
    do tile = 1, distribution%root_tile_count()
      owner = distribution%root_tiles(tile)%owner
      call MPI_Bcast( &
        candidate_state(:, :, distribution%root_tiles(tile)%j_lower: &
          distribution%root_tiles(tile)%j_upper), &
        nvar * distribution%root_tiles(tile)%cell_count, &
        MPI_DOUBLE_PRECISION, owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
      call MPI_Bcast( &
        candidate_temperature(:, distribution%root_tiles(tile)%j_lower: &
          distribution%root_tiles(tile)%j_upper), &
        distribution%root_tiles(tile)%cell_count, MPI_DOUBLE_PRECISION, &
        owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
    end do
    do child = 1, distribution%child_count()
      owner = distribution%child_owners(child)
      call MPI_Bcast( &
        candidate_set%children(child)%state, &
        size(candidate_set%children(child)%state), MPI_DOUBLE_PRECISION, &
        owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
      call MPI_Bcast( &
        candidate_set%children(child)%temperature, &
        size(candidate_set%children(child)%temperature), &
        MPI_DOUBLE_PRECISION, owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
    end do
    local_ok = candidate_set%is_valid(coarse_geometry, nvar) .and. &
      all(ieee_is_finite(candidate_state)) .and. &
      all(ieee_is_finite(candidate_temperature))
    call MPI_Allreduce( &
      local_ok, global_ok, 1, MPI_LOGICAL, MPI_LAND, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. .not. global_ok) return
    synchronized_coarse_state = candidate_state
    synchronized_coarse_temperature = candidate_temperature
    synchronized_patch_set = candidate_set
    ok = .true.
  end subroutine synchronize_owned_reactive_eb_patch_set_2d

  subroutine advance_owned_reactive_eb_patch_set_chemistry_2d( &
      species, reactions, interval, rtol, atol, distribution, &
      coarse_state, coarse_temperature, coarse_geometry, patch_set, ok, &
      local_entity_advances)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: interval, rtol, atol
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    real(dp), intent(inout) :: coarse_state(:, :, :)
    real(dp), intent(inout) :: coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(inout) :: patch_set
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_entity_advances

    type(reactive_eb_patch_set_2d) :: candidate_set
    type(reactive_eb_patch_set_2d) :: synchronized_set
    real(dp), allocatable :: candidate_state(:, :, :)
    real(dp), allocatable :: candidate_temperature(:, :)
    real(dp), allocatable :: averaged_state(:, :, :)
    real(dp), allocatable :: averaged_temperature(:, :)
    logical, allocatable :: active_mask(:, :)
    real(dp) :: controls(3), control_minimum(3), control_maximum(3)
    logical :: accepted, entity_ok, global_ok, local_ok
    integer :: child, count_maximum(2), count_minimum(2), counts(2)
    integer :: ierr, j_lower, j_upper, nvar, owner, tile
    integer :: advances

    ok = .false.
    advances = 0
    if (present(local_entity_advances)) local_entity_advances = 0
    nvar = reactive_nvar(size(species))
    controls = [interval, rtol, atol]
    counts = [size(species), size(reactions)]
    local_ok = size(species) >= 1 .and. size(reactions) >= 1 .and. &
      all(ieee_is_finite(controls)) .and. interval >= 0.0_dp .and. &
      rtol > 0.0_dp .and. atol > 0.0_dp .and. &
      all(shape(coarse_state) == &
        [nvar, coarse_geometry%nx, coarse_geometry%ny]) .and. &
      all(shape(coarse_temperature) == &
        [coarse_geometry%nx, coarse_geometry%ny]) .and. &
      distribution%is_valid(coarse_geometry, patch_set)
    call MPI_Allreduce( &
      local_ok, global_ok, 1, MPI_LOGICAL, MPI_LAND, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. .not. global_ok) return
    call MPI_Allreduce( &
      controls, control_minimum, 3, MPI_DOUBLE_PRECISION, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      controls, control_maximum, 3, MPI_DOUBLE_PRECISION, MPI_MAX, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      counts, count_minimum, 2, MPI_INTEGER, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      counts, count_maximum, 2, MPI_INTEGER, MPI_MAX, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. &
        any(control_minimum /= control_maximum) .or. &
        any(count_minimum /= count_maximum)) return

    allocate(candidate_state, mold=coarse_state)
    allocate(candidate_temperature, mold=coarse_temperature)
    call synchronize_owned_reactive_eb_patch_set_2d( &
      distribution, size(species), coarse_state, coarse_temperature, &
      coarse_geometry, patch_set, candidate_state, candidate_temperature, &
      synchronized_set, local_ok)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    candidate_set = synchronized_set

    do tile = 1, distribution%root_tile_count()
      owner = distribution%root_tiles(tile)%owner
      j_lower = distribution%root_tiles(tile)%j_lower
      j_upper = distribution%root_tiles(tile)%j_upper
      entity_ok = .true.
      if (distribution%rank == owner) then
        active_mask = &
          coarse_geometry%cell_type(:, j_lower:j_upper) /= eb_covered_cell
        call advance_reactive_chemistry_2d( &
          species, reactions, candidate_state(:, :, j_lower:j_upper), &
          candidate_temperature(:, j_lower:j_upper), coarse_geometry%nx, &
          j_upper - j_lower + 1, interval, rtol, atol, entity_ok, active_mask)
        if (entity_ok) advances = advances + 1
        deallocate(active_mask)
      end if
      call all_ranks_accept_eb_2d( &
        distribution, entity_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) return
      call MPI_Bcast( &
        candidate_state(:, :, j_lower:j_upper), &
        nvar * distribution%root_tiles(tile)%cell_count, &
        MPI_DOUBLE_PRECISION, owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
      call MPI_Bcast( &
        candidate_temperature(:, j_lower:j_upper), &
        distribution%root_tiles(tile)%cell_count, MPI_DOUBLE_PRECISION, &
        owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
    end do

    do child = 1, distribution%child_count()
      owner = distribution%child_owners(child)
      entity_ok = .true.
      if (distribution%rank == owner) then
        active_mask = candidate_set%children(child)%geometry%cell_type /= &
          eb_covered_cell
        call advance_reactive_chemistry_2d( &
          species, reactions, candidate_set%children(child)%state, &
          candidate_set%children(child)%temperature, &
          candidate_set%children(child)%geometry%nx, &
          candidate_set%children(child)%geometry%ny, interval, rtol, atol, &
          entity_ok, active_mask)
        if (entity_ok) advances = advances + 1
        deallocate(active_mask)
      end if
      call all_ranks_accept_eb_2d( &
        distribution, entity_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) return
      call MPI_Bcast( &
        candidate_set%children(child)%state, &
        size(candidate_set%children(child)%state), MPI_DOUBLE_PRECISION, &
        owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
      call MPI_Bcast( &
        candidate_set%children(child)%temperature, &
        size(candidate_set%children(child)%temperature), &
        MPI_DOUBLE_PRECISION, owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
    end do

    allocate(averaged_state, mold=coarse_state)
    allocate(averaged_temperature, mold=coarse_temperature)
    call average_down_reactive_eb_patch_set_2d( &
      species, candidate_state, candidate_temperature, coarse_geometry, &
      candidate_set, averaged_state, averaged_temperature, local_ok)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    local_ok = candidate_set%is_valid(coarse_geometry, nvar) .and. &
      all(ieee_is_finite(averaged_state)) .and. &
      all(ieee_is_finite(averaged_temperature))
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    coarse_state = averaged_state
    coarse_temperature = averaged_temperature
    patch_set = candidate_set
    ok = .true.
    if (present(local_entity_advances)) local_entity_advances = advances
  end subroutine advance_owned_reactive_eb_patch_set_chemistry_2d

  subroutine advance_owned_reactive_eb_root_tiles_hydro_2d( &
      species, distribution, state, temperature, geometry, solver, &
      reconstruction, limiter, state_redist_target_volume_fraction, &
      state_redist_max_order, dt, new_state, new_temperature, x_flux, &
      y_flux, ok, local_tile_advances, local_computed_cells)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    type(eb_geometry_2d), intent(in) :: geometry
    character(len=*), intent(in) :: solver, reconstruction, limiter
    real(dp), intent(in) :: state_redist_target_volume_fraction, dt
    integer, intent(in) :: state_redist_max_order
    real(dp), intent(out) :: new_state(:, :, :), new_temperature(:, :)
    real(dp), intent(out) :: x_flux(:, 0:, :), y_flux(:, :, 0:)
    logical, intent(out) :: ok
    integer, intent(out) :: local_tile_advances, local_computed_cells

    type(eb_geometry_2d) :: band_geometry
    real(dp), allocatable :: band_state(:, :, :), band_temperature(:, :)
    real(dp), allocatable :: band_new_state(:, :, :)
    real(dp), allocatable :: band_new_temperature(:, :)
    real(dp), allocatable :: band_x_flux(:, :, :), band_y_flux(:, :, :)
    real(dp), allocatable :: state_contribution(:, :, :)
    real(dp), allocatable :: temperature_contribution(:, :)
    real(dp), allocatable :: x_flux_contribution(:, :, :)
    real(dp), allocatable :: y_flux_contribution(:, :, :)
    logical :: accepted, entity_ok, global_ok, local_ok
    integer :: band_j_lower, band_j_upper, face_j_lower, face_j_upper
    integer :: ierr, local_face_j_lower, local_face_j_upper
    integer :: local_j_lower, local_j_upper, nvar, tile

    new_state = 0.0_dp
    new_temperature = 0.0_dp
    x_flux = 0.0_dp
    y_flux = 0.0_dp
    ok = .false.
    local_tile_advances = 0
    local_computed_cells = 0
    nvar = reactive_nvar(size(species))
    local_ok = nvar >= 1 .and. &
      all(shape(state) == [nvar, geometry%nx, geometry%ny]) .and. &
      all(shape(temperature) == [geometry%nx, geometry%ny]) .and. &
      all(shape(new_state) == shape(state)) .and. &
      all(shape(new_temperature) == shape(temperature)) .and. &
      size(x_flux, 1) == nvar .and. &
      size(x_flux, 2) == geometry%nx + 1 .and. &
      size(x_flux, 3) == geometry%ny .and. &
      size(y_flux, 1) == nvar .and. size(y_flux, 2) == geometry%nx .and. &
      size(y_flux, 3) == geometry%ny + 1
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    allocate(state_contribution, mold=state)
    allocate(temperature_contribution, mold=temperature)
    allocate(x_flux_contribution(nvar, 0:geometry%nx, geometry%ny))
    allocate(y_flux_contribution(nvar, geometry%nx, 0:geometry%ny))
    state_contribution = 0.0_dp
    temperature_contribution = 0.0_dp
    x_flux_contribution = 0.0_dp
    y_flux_contribution = 0.0_dp
    entity_ok = .true.

    do tile = 1, distribution%root_tile_count()
      if (.not. distribution%root_tile_is_local(tile)) cycle
      band_j_lower = max(1, distribution%root_tiles(tile)%j_lower - &
        mpi_amr_eb_root_tile_hydro_halo_cells)
      band_j_upper = min(geometry%ny, &
        distribution%root_tiles(tile)%j_upper + &
          mpi_amr_eb_root_tile_hydro_halo_cells)
      call extract_eb_geometry_y_band_2d( &
        geometry, band_j_lower, band_j_upper, band_geometry, entity_ok)
      if (.not. entity_ok) exit
      allocate(band_state, source=state(:, :, band_j_lower:band_j_upper))
      allocate(band_temperature, &
        source=temperature(:, band_j_lower:band_j_upper))
      allocate(band_new_state, mold=band_state)
      allocate(band_new_temperature, mold=band_temperature)
      allocate(band_x_flux(nvar, 0:geometry%nx, band_geometry%ny))
      allocate(band_y_flux(nvar, geometry%nx, 0:band_geometry%ny))
      call advance_reactive_eb_level_2d( &
        species, band_state, band_temperature, band_geometry, solver, &
        reconstruction, limiter, state_redist_target_volume_fraction, &
        state_redist_max_order, dt, band_new_state, band_new_temperature, &
        band_x_flux, band_y_flux, entity_ok)
      if (.not. entity_ok) exit

      local_j_lower = distribution%root_tiles(tile)%j_lower - &
        band_j_lower + 1
      local_j_upper = distribution%root_tiles(tile)%j_upper - &
        band_j_lower + 1
      state_contribution(:, :, &
        distribution%root_tiles(tile)%j_lower: &
          distribution%root_tiles(tile)%j_upper) = &
        band_new_state(:, :, local_j_lower:local_j_upper)
      temperature_contribution(:, &
        distribution%root_tiles(tile)%j_lower: &
          distribution%root_tiles(tile)%j_upper) = &
        band_new_temperature(:, local_j_lower:local_j_upper)
      x_flux_contribution(:, :, &
        distribution%root_tiles(tile)%j_lower: &
          distribution%root_tiles(tile)%j_upper) = &
        band_x_flux(:, :, local_j_lower:local_j_upper)

      face_j_lower = distribution%root_tiles(tile)%j_lower - 1
      face_j_upper = distribution%root_tiles(tile)%j_upper - 1
      if (distribution%root_tiles(tile)%j_upper == geometry%ny) &
        face_j_upper = geometry%ny
      local_face_j_lower = face_j_lower - band_j_lower + 1
      local_face_j_upper = face_j_upper - band_j_lower + 1
      y_flux_contribution(:, :, face_j_lower:face_j_upper) = &
        band_y_flux(:, :, local_face_j_lower:local_face_j_upper)
      local_tile_advances = local_tile_advances + 1
      local_computed_cells = local_computed_cells + &
        geometry%nx * band_geometry%ny
      deallocate( &
        band_state, band_temperature, band_new_state, &
        band_new_temperature, band_x_flux, band_y_flux)
    end do
    call all_ranks_accept_eb_2d( &
      distribution, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) then
      local_tile_advances = 0
      local_computed_cells = 0
      return
    end if

    call MPI_Allreduce( &
      state_contribution, new_state, size(new_state), &
      MPI_DOUBLE_PRECISION, MPI_SUM, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      temperature_contribution, new_temperature, size(new_temperature), &
      MPI_DOUBLE_PRECISION, MPI_SUM, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      x_flux_contribution, x_flux, size(x_flux), MPI_DOUBLE_PRECISION, &
      MPI_SUM, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      y_flux_contribution, y_flux, size(y_flux), MPI_DOUBLE_PRECISION, &
      MPI_SUM, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    ok = all(ieee_is_finite(new_state)) .and. &
      all(ieee_is_finite(new_temperature)) .and. &
      all(ieee_is_finite(x_flux)) .and. all(ieee_is_finite(y_flux))
    if (.not. ok) then
      local_tile_advances = 0
      local_computed_cells = 0
    end if
  end subroutine advance_owned_reactive_eb_root_tiles_hydro_2d

  subroutine extract_eb_geometry_y_rows_2d( &
      geometry, source_rows, band, ok)
    type(eb_geometry_2d), intent(in) :: geometry
    integer, intent(in) :: source_rows(:)
    type(eb_geometry_2d), intent(out) :: band
    logical, intent(out) :: ok

    integer :: j, source_j
    real(dp) :: y_shift

    band = eb_geometry_2d()
    ok = .false.
    if (.not. geometry%is_valid() .or. size(source_rows) < 1 .or. &
        size(source_rows) > geometry%ny .or. source_rows(1) /= 1 .or. &
        source_rows(size(source_rows)) /= geometry%ny .or. &
        any(source_rows < 1) .or. any(source_rows > geometry%ny)) return
    if (size(source_rows) > 1) then
      if (any(source_rows(2:) <= source_rows(:size(source_rows) - 1))) return
    end if

    band%nx = geometry%nx
    band%ny = size(source_rows)
    band%x_lower = geometry%x_lower
    band%x_upper = geometry%x_upper
    band%y_lower = geometry%y_lower
    band%y_upper = geometry%y_lower + real(band%ny, dp) * geometry%dy
    band%dx = geometry%dx
    band%dy = geometry%dy
    allocate(band%volume_fraction(band%nx, band%ny))
    allocate(band%cell_centroid_x(band%nx, band%ny))
    allocate(band%cell_centroid_y(band%nx, band%ny))
    allocate(band%cell_type(band%nx, band%ny))
    allocate(band%x_face_fraction(0:band%nx, band%ny))
    allocate(band%y_face_fraction(band%nx, 0:band%ny))
    allocate(band%x_face_centroid_y(0:band%nx, band%ny))
    allocate(band%y_face_centroid_x(band%nx, 0:band%ny))
    allocate(band%boundary_length(band%nx, band%ny))
    allocate(band%boundary_centroid_x(band%nx, band%ny))
    allocate(band%boundary_centroid_y(band%nx, band%ny))
    allocate(band%boundary_normal_x(band%nx, band%ny))
    allocate(band%boundary_normal_y(band%nx, band%ny))
    allocate(band%boundary_normal_integral_x(band%nx, band%ny))
    allocate(band%boundary_normal_integral_y(band%nx, band%ny))

    do j = 1, band%ny
      source_j = source_rows(j)
      band%volume_fraction(:, j) = geometry%volume_fraction(:, source_j)
      band%cell_centroid_x(:, j) = geometry%cell_centroid_x(:, source_j)
      band%cell_centroid_y(:, j) = geometry%cell_centroid_y(:, source_j)
      band%cell_type(:, j) = geometry%cell_type(:, source_j)
      band%x_face_fraction(:, j) = geometry%x_face_fraction(:, source_j)
      band%x_face_centroid_y(:, j) = &
        geometry%x_face_centroid_y(:, source_j)
      band%boundary_length(:, j) = geometry%boundary_length(:, source_j)
      band%boundary_centroid_x(:, j) = &
        geometry%boundary_centroid_x(:, source_j)
      band%boundary_centroid_y(:, j) = &
        geometry%boundary_centroid_y(:, source_j)
      y_shift = real(j - source_j, dp) * geometry%dy
      where (band%boundary_length(:, j) > 0.0_dp)
        band%boundary_centroid_y(:, j) = &
          band%boundary_centroid_y(:, j) + y_shift
      end where
      band%boundary_normal_x(:, j) = &
        geometry%boundary_normal_x(:, source_j)
      band%boundary_normal_y(:, j) = &
        geometry%boundary_normal_y(:, source_j)
      band%boundary_normal_integral_x(:, j) = &
        geometry%boundary_normal_integral_x(:, source_j)
      band%boundary_normal_integral_y(:, j) = &
        geometry%boundary_normal_integral_y(:, source_j)
    end do
    band%y_face_fraction(:, 0) = geometry%y_face_fraction(:, 0)
    band%y_face_centroid_x(:, 0) = geometry%y_face_centroid_x(:, 0)
    do j = 1, band%ny
      source_j = source_rows(j)
      band%y_face_fraction(:, j) = &
        geometry%y_face_fraction(:, source_j)
      band%y_face_centroid_x(:, j) = &
        geometry%y_face_centroid_x(:, source_j)
    end do
    ok = band%is_valid()
  end subroutine extract_eb_geometry_y_rows_2d

  subroutine extract_eb_geometry_y_band_2d( &
      geometry, j_lower, j_upper, band, ok)
    type(eb_geometry_2d), intent(in) :: geometry
    integer, intent(in) :: j_lower, j_upper
    type(eb_geometry_2d), intent(out) :: band
    logical, intent(out) :: ok

    integer :: ny

    band = eb_geometry_2d()
    ok = .false.
    if (.not. geometry%is_valid() .or. j_lower < 1 .or. &
        j_upper > geometry%ny .or. j_upper < j_lower) return
    ny = j_upper - j_lower + 1
    band%nx = geometry%nx
    band%ny = ny
    band%x_lower = geometry%x_lower
    band%x_upper = geometry%x_upper
    band%y_lower = geometry%y_lower + &
      real(j_lower - 1, dp) * geometry%dy
    band%y_upper = geometry%y_lower + real(j_upper, dp) * geometry%dy
    band%dx = geometry%dx
    band%dy = geometry%dy
    allocate(band%volume_fraction(band%nx, ny))
    allocate(band%cell_centroid_x(band%nx, ny))
    allocate(band%cell_centroid_y(band%nx, ny))
    allocate(band%cell_type(band%nx, ny))
    allocate(band%x_face_fraction(0:band%nx, ny))
    allocate(band%y_face_fraction(band%nx, 0:ny))
    allocate(band%x_face_centroid_y(0:band%nx, ny))
    allocate(band%y_face_centroid_x(band%nx, 0:ny))
    allocate(band%boundary_length(band%nx, ny))
    allocate(band%boundary_centroid_x(band%nx, ny))
    allocate(band%boundary_centroid_y(band%nx, ny))
    allocate(band%boundary_normal_x(band%nx, ny))
    allocate(band%boundary_normal_y(band%nx, ny))
    allocate(band%boundary_normal_integral_x(band%nx, ny))
    allocate(band%boundary_normal_integral_y(band%nx, ny))
    band%volume_fraction = geometry%volume_fraction(:, j_lower:j_upper)
    band%cell_centroid_x = geometry%cell_centroid_x(:, j_lower:j_upper)
    band%cell_centroid_y = geometry%cell_centroid_y(:, j_lower:j_upper)
    band%cell_type = geometry%cell_type(:, j_lower:j_upper)
    band%x_face_fraction = geometry%x_face_fraction(:, j_lower:j_upper)
    band%x_face_centroid_y = &
      geometry%x_face_centroid_y(:, j_lower:j_upper)
    band%y_face_fraction = &
      geometry%y_face_fraction(:, j_lower - 1:j_upper)
    band%y_face_centroid_x = &
      geometry%y_face_centroid_x(:, j_lower - 1:j_upper)
    band%boundary_length = geometry%boundary_length(:, j_lower:j_upper)
    band%boundary_centroid_x = &
      geometry%boundary_centroid_x(:, j_lower:j_upper)
    band%boundary_centroid_y = &
      geometry%boundary_centroid_y(:, j_lower:j_upper)
    band%boundary_normal_x = &
      geometry%boundary_normal_x(:, j_lower:j_upper)
    band%boundary_normal_y = &
      geometry%boundary_normal_y(:, j_lower:j_upper)
    band%boundary_normal_integral_x = &
      geometry%boundary_normal_integral_x(:, j_lower:j_upper)
    band%boundary_normal_integral_y = &
      geometry%boundary_normal_integral_y(:, j_lower:j_upper)
    ok = band%is_valid()
  end subroutine extract_eb_geometry_y_band_2d

  subroutine advance_owned_reactive_eb_patch_set_hydro_2d( &
      species, distribution, coarse_state, coarse_temperature, &
      coarse_geometry, patch_set, solver, reconstruction, limiter, &
      state_redist_max_order, dt, ok, local_level_advances, &
      state_redist_target_volume_fraction, local_root_hydro_cells)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    real(dp), intent(inout) :: coarse_state(:, :, :)
    real(dp), intent(inout) :: coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(inout) :: patch_set
    character(len=*), intent(in) :: solver, reconstruction, limiter
    integer, intent(in) :: state_redist_max_order
    real(dp), intent(in) :: dt
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_level_advances
    real(dp), intent(in), optional :: state_redist_target_volume_fraction
    integer, intent(out), optional :: local_root_hydro_cells

    type(amr_eb_flux_register_2d) :: flux_register
    type(reactive_eb_exterior_state_2d) :: exterior
    type(reactive_eb_patch_set_2d) :: candidate_set, synchronized_set
    real(dp), allocatable :: averaged_state(:, :, :)
    real(dp), allocatable :: averaged_temperature(:, :)
    real(dp), allocatable :: coarse_corrected(:, :, :)
    real(dp), allocatable :: coarse_corrected_temperature(:, :)
    real(dp), allocatable :: coarse_work(:, :, :)
    real(dp), allocatable :: coarse_work_temperature(:, :)
    real(dp), allocatable :: coarse_x_flux(:, :, :)
    real(dp), allocatable :: coarse_y_flux(:, :, :)
    real(dp), allocatable :: fine_work(:, :, :)
    real(dp), allocatable :: fine_work_temperature(:, :)
    real(dp), allocatable :: fine_x_flux(:, :, :)
    real(dp), allocatable :: fine_y_flux(:, :, :)
    real(dp), allocatable :: root_start(:, :, :)
    real(dp), allocatable :: root_start_temperature(:, :)
    real(dp), allocatable :: root_hydro(:, :, :)
    real(dp), allocatable :: root_hydro_temperature(:, :)
    real(dp) :: alpha, fine_dt, numeric_controls(2)
    real(dp) :: numeric_maximum(2), numeric_minimum(2), selected_target
    logical :: accepted, entity_ok, global_ok, local_ok
    integer :: advances, character_index, child, ierr, integer_controls(2)
    integer :: integer_maximum(2), integer_minimum(2)
    integer :: nvar, owner, ratio, root_hydro_cells, root_owner
    integer :: root_tile_advances, substep
    integer :: string_codes(32, 3), string_maximum(32, 3)
    integer :: string_minimum(32, 3)

    ok = .false.
    advances = 0
    if (present(local_level_advances)) local_level_advances = 0
    if (present(local_root_hydro_cells)) local_root_hydro_cells = 0
    selected_target = 0.5_dp
    if (present(state_redist_target_volume_fraction)) &
      selected_target = state_redist_target_volume_fraction
    nvar = reactive_nvar(size(species))
    numeric_controls = [dt, selected_target]
    integer_controls = [state_redist_max_order, size(species)]
    string_codes = 0
    local_ok = len_trim(solver) >= 1 .and. len_trim(solver) <= 32 .and. &
      len_trim(reconstruction) >= 1 .and. &
      len_trim(reconstruction) <= 32 .and. &
      len_trim(limiter) >= 1 .and. len_trim(limiter) <= 32
    if (local_ok) then
      do character_index = 1, len_trim(solver)
        string_codes(character_index, 1) = &
          iachar(solver(character_index:character_index))
      end do
      do character_index = 1, len_trim(reconstruction)
        string_codes(character_index, 2) = &
          iachar(reconstruction(character_index:character_index))
      end do
      do character_index = 1, len_trim(limiter)
        string_codes(character_index, 3) = &
          iachar(limiter(character_index:character_index))
      end do
    end if
    local_ok = local_ok .and. size(species) >= 1 .and. &
      all(ieee_is_finite(numeric_controls)) .and. dt > 0.0_dp .and. &
      selected_target > 0.0_dp .and. selected_target <= 1.0_dp .and. &
      (state_redist_max_order == 0 .or. state_redist_max_order == 2) .and. &
      all(shape(coarse_state) == &
        [nvar, coarse_geometry%nx, coarse_geometry%ny]) .and. &
      all(shape(coarse_temperature) == &
        [coarse_geometry%nx, coarse_geometry%ny]) .and. &
      distribution%is_valid(coarse_geometry, patch_set)
    call MPI_Allreduce( &
      local_ok, global_ok, 1, MPI_LOGICAL, MPI_LAND, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. .not. global_ok) return
    call MPI_Allreduce( &
      numeric_controls, numeric_minimum, 2, MPI_DOUBLE_PRECISION, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      numeric_controls, numeric_maximum, 2, MPI_DOUBLE_PRECISION, MPI_MAX, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_controls, integer_minimum, 2, MPI_INTEGER, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_controls, integer_maximum, 2, MPI_INTEGER, MPI_MAX, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. &
        any(numeric_minimum /= numeric_maximum) .or. &
        any(integer_minimum /= integer_maximum)) return
    call MPI_Allreduce( &
      string_codes, string_minimum, size(string_codes), MPI_INTEGER, &
      MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      string_codes, string_maximum, size(string_codes), MPI_INTEGER, &
      MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. &
        any(string_minimum /= string_maximum)) return

    allocate(root_start, mold=coarse_state)
    allocate(root_start_temperature, mold=coarse_temperature)
    call synchronize_owned_reactive_eb_patch_set_2d( &
      distribution, size(species), coarse_state, coarse_temperature, &
      coarse_geometry, patch_set, root_start, root_start_temperature, &
      synchronized_set, local_ok)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    candidate_set = synchronized_set

    allocate(root_hydro, mold=coarse_state)
    allocate(root_hydro_temperature, mold=coarse_temperature)
    allocate(coarse_x_flux(nvar, 0:coarse_geometry%nx, coarse_geometry%ny))
    allocate(coarse_y_flux(nvar, coarse_geometry%nx, 0:coarse_geometry%ny))
    call advance_owned_reactive_eb_root_tiles_hydro_2d( &
      species, distribution, root_start, root_start_temperature, &
      coarse_geometry, trim(solver), trim(reconstruction), trim(limiter), &
      selected_target, state_redist_max_order, dt, root_hydro, &
      root_hydro_temperature, coarse_x_flux, coarse_y_flux, entity_ok, &
      root_tile_advances, root_hydro_cells)
    if (.not. entity_ok) return
    advances = advances + root_tile_advances
    root_owner = distribution%root_level_owner()

    allocate(coarse_corrected, source=root_hydro)
    allocate(coarse_corrected_temperature, source=root_hydro_temperature)
    allocate(coarse_work, mold=coarse_state)
    allocate(coarse_work_temperature, mold=coarse_temperature)
    do child = 1, candidate_set%patch_count()
      owner = distribution%child_owner(child)
      entity_ok = owner >= 0 .and. owner < distribution%nranks
      if (distribution%rank == owner .and. entity_ok) then
        call initialize_amr_eb_flux_register_2d( &
          coarse_geometry, candidate_set%children(child)%geometry, &
          candidate_set%children(child)%patch, nvar, flux_register, entity_ok)
        if (entity_ok) call accumulate_coarse_eb_fluxes_2d( &
          flux_register, coarse_geometry, &
          candidate_set%children(child)%geometry, &
          candidate_set%children(child)%patch, coarse_x_flux, coarse_y_flux, &
          dt, entity_ok)
        if (allocated(fine_work)) deallocate(fine_work)
        if (allocated(fine_work_temperature)) &
          deallocate(fine_work_temperature)
        if (allocated(fine_x_flux)) deallocate(fine_x_flux)
        if (allocated(fine_y_flux)) deallocate(fine_y_flux)
        allocate(fine_work, mold=candidate_set%children(child)%state)
        allocate(fine_work_temperature, &
          mold=candidate_set%children(child)%temperature)
        allocate(fine_x_flux(nvar, &
          0:candidate_set%children(child)%geometry%nx, &
          candidate_set%children(child)%geometry%ny))
        allocate(fine_y_flux(nvar, &
          candidate_set%children(child)%geometry%nx, &
          0:candidate_set%children(child)%geometry%ny))
        ratio = candidate_set%children(child)%patch%refinement_ratio
        fine_dt = dt / real(ratio, dp)
        do substep = 1, ratio
          if (.not. entity_ok) exit
          if (trim(reconstruction) == "characteristic_plm") then
            alpha = (real(substep, dp) - 0.5_dp) / real(ratio, dp)
          else
            alpha = real(substep - 1, dp) / real(ratio, dp)
          end if
          call build_reactive_eb_patch_exterior_2d( &
            species, root_start, root_start_temperature, root_hydro, &
            root_hydro_temperature, coarse_geometry, &
            candidate_set%children(child)%geometry, &
            candidate_set%children(child)%patch, alpha, exterior, entity_ok, &
            candidate_set%children(child)%state, &
            candidate_set%children(child)%temperature)
          if (.not. entity_ok) exit
          call advance_reactive_eb_level_2d( &
            species, candidate_set%children(child)%state, &
            candidate_set%children(child)%temperature, &
            candidate_set%children(child)%geometry, trim(solver), &
            trim(reconstruction), trim(limiter), selected_target, &
            state_redist_max_order, fine_dt, fine_work, &
            fine_work_temperature, fine_x_flux, fine_y_flux, entity_ok, &
            exterior)
          if (.not. entity_ok) exit
          advances = advances + 1
          candidate_set%children(child)%state = fine_work
          candidate_set%children(child)%temperature = fine_work_temperature
          call accumulate_fine_eb_fluxes_2d( &
            flux_register, coarse_geometry, &
            candidate_set%children(child)%geometry, &
            candidate_set%children(child)%patch, fine_x_flux, fine_y_flux, &
            fine_dt, entity_ok)
        end do
        if (entity_ok) call reflux_reactive_eb_state_patch_2d( &
          species, coarse_corrected, coarse_corrected_temperature, &
          coarse_geometry, candidate_set%children(child)%state, &
          candidate_set%children(child)%temperature, &
          candidate_set%children(child)%geometry, &
          candidate_set%children(child)%patch, flux_register, coarse_work, &
          coarse_work_temperature, fine_work, fine_work_temperature, entity_ok)
        if (entity_ok) then
          coarse_corrected = coarse_work
          coarse_corrected_temperature = coarse_work_temperature
          candidate_set%children(child)%state = fine_work
          candidate_set%children(child)%temperature = fine_work_temperature
        end if
      end if
      call all_ranks_accept_eb_2d( &
        distribution, entity_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) return
      call MPI_Bcast( &
        coarse_corrected, size(coarse_corrected), MPI_DOUBLE_PRECISION, &
        owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
      call MPI_Bcast( &
        coarse_corrected_temperature, size(coarse_corrected_temperature), &
        MPI_DOUBLE_PRECISION, owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
      call MPI_Bcast( &
        candidate_set%children(child)%state, &
        size(candidate_set%children(child)%state), MPI_DOUBLE_PRECISION, &
        owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
      call MPI_Bcast( &
        candidate_set%children(child)%temperature, &
        size(candidate_set%children(child)%temperature), &
        MPI_DOUBLE_PRECISION, owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
    end do

    allocate(averaged_state, mold=coarse_state)
    allocate(averaged_temperature, mold=coarse_temperature)
    entity_ok = .true.
    if (distribution%rank == root_owner) call &
      average_down_reactive_eb_patch_set_2d( &
        species, coarse_corrected, coarse_corrected_temperature, &
        coarse_geometry, candidate_set, averaged_state, &
        averaged_temperature, entity_ok)
    call all_ranks_accept_eb_2d( &
      distribution, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Bcast( &
      averaged_state, size(averaged_state), MPI_DOUBLE_PRECISION, &
      root_owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Bcast( &
      averaged_temperature, size(averaged_temperature), &
      MPI_DOUBLE_PRECISION, root_owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    local_ok = candidate_set%is_valid(coarse_geometry, nvar) .and. &
      all(ieee_is_finite(averaged_state)) .and. &
      all(ieee_is_finite(averaged_temperature))
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    coarse_state = averaged_state
    coarse_temperature = averaged_temperature
    patch_set = candidate_set
    ok = .true.
    if (present(local_level_advances)) local_level_advances = advances
    if (present(local_root_hydro_cells)) &
      local_root_hydro_cells = root_hydro_cells
  end subroutine advance_owned_reactive_eb_patch_set_hydro_2d

  subroutine advance_owned_reactive_eb_patch_set_transport_2d( &
      species, transport, distribution, coarse_state, coarse_temperature, &
      coarse_geometry, patch_set, interval, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, state_redist_max_order, ok, &
      local_euler_advances, minimum_theta, &
      state_redist_target_volume_fraction)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    real(dp), intent(inout) :: coarse_state(:, :, :)
    real(dp), intent(inout) :: coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(inout) :: patch_set
    real(dp), intent(in) :: interval
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    integer, intent(in) :: state_redist_max_order
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_euler_advances
    real(dp), intent(out), optional :: minimum_theta
    real(dp), intent(in), optional :: state_redist_target_volume_fraction

    type(reactive_eb_patch_set_2d) :: candidate_set, euler_set, stage_set
    real(dp), allocatable :: candidate_state(:, :, :)
    real(dp), allocatable :: candidate_temperature(:, :)
    real(dp), allocatable :: euler_state(:, :, :)
    real(dp), allocatable :: euler_temperature(:, :)
    real(dp), allocatable :: stage_state(:, :, :)
    real(dp), allocatable :: stage_temperature(:, :)
    real(dp), allocatable :: start_state(:, :, :)
    real(dp), allocatable :: start_temperature(:, :)
    real(dp), allocatable :: synchronized_state(:, :, :)
    real(dp), allocatable :: synchronized_temperature(:, :)
    type(reactive_eb_patch_set_2d) :: start_set
    real(dp) :: selected_target, theta_one, theta_two
    logical :: accepted, entity_ok, global_ok, local_ok
    integer :: advances_one, advances_two, child, ierr, nvar, owner
    integer :: root_owner

    ok = .false.
    if (present(local_euler_advances)) local_euler_advances = 0
    if (present(minimum_theta)) minimum_theta = 1.0_dp
    selected_target = 0.5_dp
    if (present(state_redist_target_volume_fraction)) &
      selected_target = state_redist_target_volume_fraction
    nvar = reactive_nvar(size(species))
    local_ok = all(shape(coarse_state) == &
        [nvar, coarse_geometry%nx, coarse_geometry%ny]) .and. &
      all(shape(coarse_temperature) == &
        [coarse_geometry%nx, coarse_geometry%ny])
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call collective_transport_preflight_2d( &
      species, transport, distribution, coarse_geometry, patch_set, &
      interval, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, state_redist_max_order, &
      selected_target, local_ok)
    if (.not. local_ok) return
    if (interval <= tiny(1.0_dp) .or. .not. (viscosity_enabled .or. &
        thermal_conduction_enabled .or. species_diffusion_enabled)) then
      ok = .true.
      return
    end if

    allocate(start_state, mold=coarse_state)
    allocate(start_temperature, mold=coarse_temperature)
    call synchronize_owned_reactive_eb_patch_set_2d( &
      distribution, size(species), coarse_state, coarse_temperature, &
      coarse_geometry, patch_set, start_state, start_temperature, &
      start_set, local_ok)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    allocate(stage_state, mold=coarse_state)
    allocate(stage_temperature, mold=coarse_temperature)
    call advance_owned_reactive_eb_patch_set_transport_euler_2d( &
      species, transport, distribution, start_state, start_temperature, &
      coarse_geometry, start_set, interval, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, state_redist_max_order, &
      selected_target, stage_state, stage_temperature, stage_set, theta_one, &
      local_ok, advances_one)
    if (.not. local_ok) return

    allocate(euler_state, mold=coarse_state)
    allocate(euler_temperature, mold=coarse_temperature)
    call advance_owned_reactive_eb_patch_set_transport_euler_2d( &
      species, transport, distribution, stage_state, stage_temperature, &
      coarse_geometry, stage_set, interval, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, state_redist_max_order, &
      selected_target, euler_state, euler_temperature, euler_set, theta_two, &
      local_ok, advances_two)
    if (.not. local_ok) return

    nvar = reactive_nvar(size(species))
    root_owner = distribution%root_level_owner()
    allocate(candidate_state, mold=coarse_state)
    allocate(candidate_temperature, mold=coarse_temperature)
    candidate_set = start_set
    entity_ok = root_owner >= 0 .and. root_owner < distribution%nranks
    if (distribution%rank == root_owner .and. entity_ok) then
      candidate_state = 0.5_dp * (start_state + euler_state)
      call recover_transport_temperature_2d( &
        species, candidate_state, &
        0.5_dp * (start_temperature + euler_temperature), &
        coarse_geometry, candidate_temperature, entity_ok)
    end if
    call all_ranks_accept_eb_2d( &
      distribution, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Bcast( &
      candidate_state, size(candidate_state), MPI_DOUBLE_PRECISION, &
      root_owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Bcast( &
      candidate_temperature, size(candidate_temperature), &
      MPI_DOUBLE_PRECISION, root_owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return

    do child = 1, candidate_set%patch_count()
      owner = distribution%child_owner(child)
      entity_ok = owner >= 0 .and. owner < distribution%nranks
      if (distribution%rank == owner .and. entity_ok) then
        candidate_set%children(child)%state = 0.5_dp * &
          (start_set%children(child)%state + &
           euler_set%children(child)%state)
        call recover_transport_temperature_2d( &
          species, candidate_set%children(child)%state, &
          0.5_dp * (start_set%children(child)%temperature + &
            euler_set%children(child)%temperature), &
          candidate_set%children(child)%geometry, &
          candidate_set%children(child)%temperature, entity_ok)
      end if
      call all_ranks_accept_eb_2d( &
        distribution, entity_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) return
      call MPI_Bcast( &
        candidate_set%children(child)%state, &
        size(candidate_set%children(child)%state), MPI_DOUBLE_PRECISION, &
        owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
      call MPI_Bcast( &
        candidate_set%children(child)%temperature, &
        size(candidate_set%children(child)%temperature), &
        MPI_DOUBLE_PRECISION, owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
    end do

    allocate(synchronized_state, mold=coarse_state)
    allocate(synchronized_temperature, mold=coarse_temperature)
    entity_ok = .true.
    if (distribution%rank == root_owner) call &
      average_down_reactive_eb_patch_set_2d( &
        species, candidate_state, candidate_temperature, coarse_geometry, &
        candidate_set, synchronized_state, synchronized_temperature, &
        entity_ok)
    call all_ranks_accept_eb_2d( &
      distribution, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Bcast( &
      synchronized_state, size(synchronized_state), MPI_DOUBLE_PRECISION, &
      root_owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Bcast( &
      synchronized_temperature, size(synchronized_temperature), &
      MPI_DOUBLE_PRECISION, root_owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    local_ok = candidate_set%is_valid(coarse_geometry, nvar) .and. &
      all(ieee_is_finite(synchronized_state)) .and. &
      all(ieee_is_finite(synchronized_temperature))
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    coarse_state = synchronized_state
    coarse_temperature = synchronized_temperature
    patch_set = candidate_set
    ok = .true.
    if (present(local_euler_advances)) &
      local_euler_advances = advances_one + advances_two
    if (present(minimum_theta)) minimum_theta = min(theta_one, theta_two)
  end subroutine advance_owned_reactive_eb_patch_set_transport_2d

  subroutine advance_owned_reactive_eb_patch_set_transport_euler_2d( &
      species, transport, distribution, coarse_state, coarse_temperature, &
      coarse_geometry, patch_set, dt, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, state_redist_max_order, &
      target_volume_fraction, new_coarse_state, new_coarse_temperature, &
      new_patch_set, minimum_theta, ok, local_euler_advances)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    real(dp), intent(in) :: coarse_state(:, :, :)
    real(dp), intent(in) :: coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set
    real(dp), intent(in) :: dt
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    integer, intent(in) :: state_redist_max_order
    real(dp), intent(in) :: target_volume_fraction
    real(dp), intent(out) :: new_coarse_state(:, :, :)
    real(dp), intent(out) :: new_coarse_temperature(:, :)
    type(reactive_eb_patch_set_2d), intent(out) :: new_patch_set
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    integer, intent(out) :: local_euler_advances

    type(amr_eb_flux_register_2d) :: flux_register
    type(reactive_eb_exterior_state_2d) :: exterior
    type(reactive_eb_patch_set_2d) :: candidate_set, synchronized_set
    real(dp), allocatable :: averaged_state(:, :, :)
    real(dp), allocatable :: averaged_temperature(:, :)
    real(dp), allocatable :: closed_state(:, :, :)
    real(dp), allocatable :: closed_temperature(:, :)
    real(dp), allocatable :: coarse_candidate(:, :, :)
    real(dp), allocatable :: coarse_candidate_temperature(:, :)
    real(dp), allocatable :: coarse_corrected(:, :, :)
    real(dp), allocatable :: coarse_corrected_temperature(:, :)
    real(dp), allocatable :: coarse_rhs(:, :, :)
    real(dp), allocatable :: coarse_work(:, :, :)
    real(dp), allocatable :: coarse_work_temperature(:, :)
    real(dp), allocatable :: coarse_x_flux(:, :, :)
    real(dp), allocatable :: coarse_y_flux(:, :, :)
    real(dp), allocatable :: fine_rhs(:, :, :)
    real(dp), allocatable :: fine_work(:, :, :)
    real(dp), allocatable :: fine_work_temperature(:, :)
    real(dp), allocatable :: fine_x_flux(:, :, :)
    real(dp), allocatable :: fine_y_flux(:, :, :)
    real(dp), allocatable :: integral_before(:)
    real(dp), allocatable :: root_start(:, :, :)
    real(dp), allocatable :: root_start_temperature(:, :)
    real(dp) :: alpha, coarse_theta, fine_dt, fine_theta, local_theta
    logical :: accepted, cut_interface, entity_ok, global_ok, local_ok
    integer :: advances, child, ierr, nvar, owner, ratio, root_owner
    integer :: substep

    new_coarse_state = coarse_state
    new_coarse_temperature = coarse_temperature
    new_patch_set = patch_set
    minimum_theta = 1.0_dp
    local_euler_advances = 0
    ok = .false.
    advances = 0
    local_theta = 1.0_dp
    nvar = reactive_nvar(size(species))

    allocate(root_start, mold=coarse_state)
    allocate(root_start_temperature, mold=coarse_temperature)
    call synchronize_owned_reactive_eb_patch_set_2d( &
      distribution, size(species), coarse_state, coarse_temperature, &
      coarse_geometry, patch_set, root_start, root_start_temperature, &
      synchronized_set, local_ok)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    candidate_set = synchronized_set
    root_owner = distribution%root_level_owner()

    allocate(integral_before(nvar))
    entity_ok = root_owner >= 0 .and. root_owner < distribution%nranks
    if (distribution%rank == root_owner .and. entity_ok) call &
      composite_reactive_eb_patch_set_integral_2d( &
        root_start, coarse_geometry, candidate_set, integral_before, &
        entity_ok)
    call all_ranks_accept_eb_2d( &
      distribution, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Bcast( &
      integral_before, size(integral_before), MPI_DOUBLE_PRECISION, &
      root_owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return

    allocate(coarse_candidate, mold=coarse_state)
    allocate(coarse_candidate_temperature, mold=coarse_temperature)
    allocate(coarse_rhs, mold=coarse_state)
    allocate(coarse_x_flux(nvar, 0:coarse_geometry%nx, coarse_geometry%ny))
    allocate(coarse_y_flux(nvar, coarse_geometry%nx, 0:coarse_geometry%ny))
    entity_ok = root_owner >= 0 .and. root_owner < distribution%nranks
    if (distribution%rank == root_owner .and. entity_ok) then
      call reactive_eb_transport_fluxes_rhs_2d( &
        species, transport, root_start, root_start_temperature, &
        coarse_geometry, dt, viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, boundaries, &
        coarse_rhs, coarse_x_flux, coarse_y_flux, coarse_theta, entity_ok)
      if (entity_ok) call advance_reactive_eb_state_redistributed_2d( &
        species, root_start, root_start_temperature, coarse_geometry, &
        coarse_rhs, dt, coarse_candidate, coarse_candidate_temperature, &
        entity_ok, target_volume_fraction, state_redist_max_order)
      if (entity_ok) then
        advances = advances + 1
        local_theta = min(local_theta, coarse_theta)
      end if
    end if
    call all_ranks_accept_eb_2d( &
      distribution, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Bcast( &
      coarse_candidate, size(coarse_candidate), MPI_DOUBLE_PRECISION, &
      root_owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Bcast( &
      coarse_candidate_temperature, size(coarse_candidate_temperature), &
      MPI_DOUBLE_PRECISION, root_owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Bcast( &
      coarse_x_flux, size(coarse_x_flux), MPI_DOUBLE_PRECISION, &
      root_owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Bcast( &
      coarse_y_flux, size(coarse_y_flux), MPI_DOUBLE_PRECISION, &
      root_owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return

    allocate(coarse_corrected, source=coarse_candidate)
    allocate(coarse_corrected_temperature, &
      source=coarse_candidate_temperature)
    allocate(coarse_work, mold=coarse_state)
    allocate(coarse_work_temperature, mold=coarse_temperature)
    cut_interface = .false.
    do child = 1, candidate_set%patch_count()
      cut_interface = cut_interface .or. .not. level_two_interface_is_regular( &
        candidate_set%children(child)%geometry)
      owner = distribution%child_owner(child)
      entity_ok = owner >= 0 .and. owner < distribution%nranks
      if (distribution%rank == owner .and. entity_ok) then
        call initialize_amr_eb_flux_register_2d( &
          coarse_geometry, candidate_set%children(child)%geometry, &
          candidate_set%children(child)%patch, nvar, flux_register, entity_ok)
        if (entity_ok) call accumulate_coarse_eb_fluxes_2d( &
          flux_register, coarse_geometry, &
          candidate_set%children(child)%geometry, &
          candidate_set%children(child)%patch, coarse_x_flux, coarse_y_flux, &
          dt, entity_ok)
        if (allocated(fine_rhs)) deallocate(fine_rhs)
        if (allocated(fine_work)) deallocate(fine_work)
        if (allocated(fine_work_temperature)) &
          deallocate(fine_work_temperature)
        if (allocated(fine_x_flux)) deallocate(fine_x_flux)
        if (allocated(fine_y_flux)) deallocate(fine_y_flux)
        allocate(fine_rhs, mold=candidate_set%children(child)%state)
        allocate(fine_work, mold=candidate_set%children(child)%state)
        allocate(fine_work_temperature, &
          mold=candidate_set%children(child)%temperature)
        allocate(fine_x_flux(nvar, &
          0:candidate_set%children(child)%geometry%nx, &
          candidate_set%children(child)%geometry%ny))
        allocate(fine_y_flux(nvar, &
          candidate_set%children(child)%geometry%nx, &
          0:candidate_set%children(child)%geometry%ny))
        ratio = candidate_set%children(child)%patch%refinement_ratio
        fine_dt = dt / real(ratio, dp)
        do substep = 1, ratio
          if (.not. entity_ok) exit
          alpha = real(substep - 1, dp) / real(ratio, dp)
          call build_reactive_eb_patch_exterior_2d( &
            species, root_start, root_start_temperature, coarse_candidate, &
            coarse_candidate_temperature, coarse_geometry, &
            candidate_set%children(child)%geometry, &
            candidate_set%children(child)%patch, alpha, exterior, entity_ok, &
            candidate_set%children(child)%state, &
            candidate_set%children(child)%temperature)
          if (.not. entity_ok) exit
          call reactive_eb_transport_fluxes_rhs_2d( &
            species, transport, candidate_set%children(child)%state, &
            candidate_set%children(child)%temperature, &
            candidate_set%children(child)%geometry, fine_dt, &
            viscosity_enabled, thermal_conduction_enabled, &
            species_diffusion_enabled, barodiffusion_enabled, boundaries, &
            fine_rhs, fine_x_flux, fine_y_flux, fine_theta, entity_ok, exterior)
          if (.not. entity_ok) exit
          call advance_reactive_eb_state_redistributed_2d( &
            species, candidate_set%children(child)%state, &
            candidate_set%children(child)%temperature, &
            candidate_set%children(child)%geometry, fine_rhs, fine_dt, &
            fine_work, fine_work_temperature, entity_ok, &
            target_volume_fraction, state_redist_max_order)
          if (.not. entity_ok) exit
          advances = advances + 1
          local_theta = min(local_theta, fine_theta)
          candidate_set%children(child)%state = fine_work
          candidate_set%children(child)%temperature = fine_work_temperature
          call accumulate_fine_eb_fluxes_2d( &
            flux_register, coarse_geometry, &
            candidate_set%children(child)%geometry, &
            candidate_set%children(child)%patch, fine_x_flux, fine_y_flux, &
            fine_dt, entity_ok)
        end do
        if (entity_ok) call reflux_reactive_eb_state_patch_2d( &
          species, coarse_corrected, coarse_corrected_temperature, &
          coarse_geometry, candidate_set%children(child)%state, &
          candidate_set%children(child)%temperature, &
          candidate_set%children(child)%geometry, &
          candidate_set%children(child)%patch, flux_register, coarse_work, &
          coarse_work_temperature, fine_work, fine_work_temperature, entity_ok)
        if (entity_ok) then
          coarse_corrected = coarse_work
          coarse_corrected_temperature = coarse_work_temperature
          candidate_set%children(child)%state = fine_work
          candidate_set%children(child)%temperature = fine_work_temperature
        end if
      end if
      call all_ranks_accept_eb_2d( &
        distribution, entity_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) return
      call MPI_Bcast( &
        coarse_corrected, size(coarse_corrected), MPI_DOUBLE_PRECISION, &
        owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
      call MPI_Bcast( &
        coarse_corrected_temperature, size(coarse_corrected_temperature), &
        MPI_DOUBLE_PRECISION, owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
      call MPI_Bcast( &
        candidate_set%children(child)%state, &
        size(candidate_set%children(child)%state), MPI_DOUBLE_PRECISION, &
        owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
      call MPI_Bcast( &
        candidate_set%children(child)%temperature, &
        size(candidate_set%children(child)%temperature), &
        MPI_DOUBLE_PRECISION, owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
    end do

    allocate(averaged_state, mold=coarse_state)
    allocate(averaged_temperature, mold=coarse_temperature)
    entity_ok = .true.
    if (distribution%rank == root_owner) then
      call average_down_reactive_eb_patch_set_2d( &
        species, coarse_corrected, coarse_corrected_temperature, &
        coarse_geometry, candidate_set, averaged_state, averaged_temperature, &
        entity_ok)
      if (entity_ok .and. cut_interface) then
        allocate(closed_state, mold=coarse_state)
        allocate(closed_temperature, mold=coarse_temperature)
        call close_cut_patch_set_conservation_2d( &
          species, integral_before, averaged_state, averaged_temperature, &
          coarse_geometry, candidate_set, coarse_x_flux, coarse_y_flux, dt, &
          closed_state, closed_temperature, entity_ok)
        if (entity_ok) then
          averaged_state = closed_state
          averaged_temperature = closed_temperature
        end if
      end if
    end if
    call all_ranks_accept_eb_2d( &
      distribution, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Bcast( &
      averaged_state, size(averaged_state), MPI_DOUBLE_PRECISION, &
      root_owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Bcast( &
      averaged_temperature, size(averaged_temperature), &
      MPI_DOUBLE_PRECISION, root_owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      local_theta, minimum_theta, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    local_ok = candidate_set%is_valid(coarse_geometry, nvar) .and. &
      all(ieee_is_finite(averaged_state)) .and. &
      all(ieee_is_finite(averaged_temperature)) .and. &
      ieee_is_finite(minimum_theta)
    call all_ranks_accept_eb_2d( &
      distribution, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    new_coarse_state = averaged_state
    new_coarse_temperature = averaged_temperature
    new_patch_set = candidate_set
    local_euler_advances = advances
    ok = .true.
  end subroutine advance_owned_reactive_eb_patch_set_transport_euler_2d

  subroutine collective_transport_preflight_2d( &
      species, transport, distribution, coarse_geometry, patch_set, &
      interval, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, state_redist_max_order, &
      target_volume_fraction, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set
    real(dp), intent(in) :: interval, target_volume_fraction
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    integer, intent(in) :: state_redist_max_order
    logical, intent(out) :: ok

    real(dp), allocatable :: numeric_controls(:), numeric_maximum(:)
    real(dp), allocatable :: numeric_minimum(:)
    integer, allocatable :: integer_controls(:), integer_maximum(:)
    integer, allocatable :: integer_minimum(:)
    logical :: local_ok
    integer :: character_index, face, ierr, integer_index, nprim, nsp
    integer :: numeric_index, species_index
    integer :: species_maximum, species_minimum

    ok = .false.
    nsp = size(species)
    call MPI_Allreduce( &
      nsp, species_minimum, 1, MPI_INTEGER, MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      nsp, species_maximum, 1, MPI_INTEGER, MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. species_minimum /= species_maximum .or. &
        nsp < 1) return
    nprim = reactive_nprim(nsp)
    local_ok = compatible_transport_database(species, transport) .and. &
      ieee_is_finite(interval) .and. interval >= 0.0_dp .and. &
      ieee_is_finite(target_volume_fraction) .and. &
      target_volume_fraction > 0.0_dp .and. &
      target_volume_fraction <= 1.0_dp .and. &
      (state_redist_max_order == 0 .or. state_redist_max_order == 2) .and. &
      distribution%is_valid(coarse_geometry, patch_set)
    if (local_ok) call validate_reactive_boundary_set_2d(boundaries, local_ok)
    if (local_ok) then
      do face = 1, 4
        local_ok = local_ok .and. &
          size(boundaries%face(face)%inflow_primitive) == nprim .and. &
          size(boundaries%face(face)%prescribed_species_flux) == nsp
      end do
    end if
    call MPI_Allreduce( &
      local_ok, ok, 1, MPI_LOGICAL, MPI_LAND, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. .not. ok) then
      ok = .false.
      return
    end if

    allocate(numeric_controls(2 + 5 * nsp + 4 * (5 + nprim + nsp)))
    allocate(numeric_minimum(size(numeric_controls)))
    allocate(numeric_maximum(size(numeric_controls)))
    allocate(integer_controls(6 + nsp + 24 * nsp + 4 * 3 * 24))
    allocate(integer_minimum(size(integer_controls)))
    allocate(integer_maximum(size(integer_controls)))
    numeric_controls = 0.0_dp
    integer_controls = 0
    numeric_controls(1:2) = [interval, target_volume_fraction]
    numeric_index = 2
    do species_index = 1, nsp
      numeric_controls(numeric_index + 1:numeric_index + 5) = [ &
        transport(species_index)%well_depth, &
        transport(species_index)%diameter, transport(species_index)%dipole, &
        transport(species_index)%polarizability, &
        transport(species_index)%rotational_relaxation]
      numeric_index = numeric_index + 5
    end do
    do face = 1, 4
      numeric_controls(numeric_index + 1:numeric_index + 5) = [ &
        boundaries%face(face)%wall_temperature, &
        boundaries%face(face)%wall_velocity, &
        boundaries%face(face)%inflow_temperature]
      numeric_index = numeric_index + 5
      numeric_controls(numeric_index + 1:numeric_index + nprim) = &
        boundaries%face(face)%inflow_primitive
      numeric_index = numeric_index + nprim
      numeric_controls(numeric_index + 1:numeric_index + nsp) = &
        boundaries%face(face)%prescribed_species_flux
      numeric_index = numeric_index + nsp
    end do
    integer_controls(1:6) = [ &
      state_redist_max_order, nsp, nprim, &
      merge(1, 0, viscosity_enabled), &
      merge(1, 0, thermal_conduction_enabled), &
      2 * merge(1, 0, species_diffusion_enabled) + &
        merge(1, 0, barodiffusion_enabled)]
    integer_index = 6
    do species_index = 1, nsp
      integer_index = integer_index + 1
      integer_controls(integer_index) = transport(species_index)%geometry
    end do
    do species_index = 1, nsp
      do character_index = 1, 24
        integer_index = integer_index + 1
        integer_controls(integer_index) = &
          iachar(transport(species_index)%name(character_index:character_index))
      end do
    end do
    do face = 1, 4
      do character_index = 1, 24
        integer_index = integer_index + 1
        integer_controls(integer_index) = &
          iachar(boundaries%face(face)%kind(character_index:character_index))
      end do
      do character_index = 1, 24
        integer_index = integer_index + 1
        integer_controls(integer_index) = &
          iachar(boundaries%face(face)%thermal(character_index:character_index))
      end do
      do character_index = 1, 24
        integer_index = integer_index + 1
        integer_controls(integer_index) = iachar( &
          boundaries%face(face)%wall_species(character_index:character_index))
      end do
    end do
    call MPI_Allreduce( &
      numeric_controls, numeric_minimum, size(numeric_controls), &
      MPI_DOUBLE_PRECISION, MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Allreduce( &
      numeric_controls, numeric_maximum, size(numeric_controls), &
      MPI_DOUBLE_PRECISION, MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Allreduce( &
      integer_controls, integer_minimum, size(integer_controls), &
      MPI_INTEGER, MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Allreduce( &
      integer_controls, integer_maximum, size(integer_controls), &
      MPI_INTEGER, MPI_MAX, distribution%comm, ierr)
    ok = ierr == MPI_SUCCESS .and. &
      all(numeric_minimum == numeric_maximum) .and. &
      all(integer_minimum == integer_maximum)
  end subroutine collective_transport_preflight_2d

  subroutine advance_owned_reactive_eb_patch_set_strang_2d( &
      species, reactions, transport, distribution, coarse_state, &
      coarse_temperature, coarse_geometry, patch_set, solver, &
      reconstruction, limiter, state_redist_max_order, dt, rtol, atol, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, ok, &
      local_chemistry_advances, local_hydro_advances, &
      local_transport_euler_advances, minimum_transport_theta, &
      state_redist_target_volume_fraction)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    real(dp), intent(inout) :: coarse_state(:, :, :)
    real(dp), intent(inout) :: coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(inout) :: patch_set
    character(len=*), intent(in) :: solver, reconstruction, limiter
    integer, intent(in) :: state_redist_max_order
    real(dp), intent(in) :: dt, rtol, atol
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_chemistry_advances
    integer, intent(out), optional :: local_hydro_advances
    integer, intent(out), optional :: local_transport_euler_advances
    real(dp), intent(out), optional :: minimum_transport_theta
    real(dp), intent(in), optional :: state_redist_target_volume_fraction

    type(reactive_eb_patch_set_2d) :: candidate_set
    real(dp), allocatable :: candidate_state(:, :, :)
    real(dp), allocatable :: candidate_temperature(:, :)
    real(dp) :: selected_target, theta, theta_one, theta_two
    logical :: local_ok
    integer :: chemistry_one, chemistry_two, hydro_advances
    integer :: transport_one, transport_two

    ok = .false.
    if (present(local_chemistry_advances)) local_chemistry_advances = 0
    if (present(local_hydro_advances)) local_hydro_advances = 0
    if (present(local_transport_euler_advances)) &
      local_transport_euler_advances = 0
    if (present(minimum_transport_theta)) minimum_transport_theta = 1.0_dp
    selected_target = 0.5_dp
    if (present(state_redist_target_volume_fraction)) &
      selected_target = state_redist_target_volume_fraction
    candidate_set = patch_set
    allocate(candidate_state, source=coarse_state)
    allocate(candidate_temperature, source=coarse_temperature)

    call advance_owned_reactive_eb_patch_set_chemistry_2d( &
      species, reactions, 0.5_dp * dt, rtol, atol, distribution, &
      candidate_state, candidate_temperature, coarse_geometry, &
      candidate_set, local_ok, chemistry_one)
    if (.not. local_ok) return
    call advance_owned_reactive_eb_patch_set_transport_2d( &
      species, transport, distribution, candidate_state, &
      candidate_temperature, coarse_geometry, candidate_set, 0.5_dp * dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      state_redist_max_order, local_ok, transport_one, theta_one, &
      selected_target)
    if (.not. local_ok) return
    call advance_owned_reactive_eb_patch_set_hydro_2d( &
      species, distribution, candidate_state, candidate_temperature, &
      coarse_geometry, candidate_set, solver, reconstruction, limiter, &
      state_redist_max_order, dt, local_ok, hydro_advances, selected_target)
    if (.not. local_ok) return
    call advance_owned_reactive_eb_patch_set_transport_2d( &
      species, transport, distribution, candidate_state, &
      candidate_temperature, coarse_geometry, candidate_set, 0.5_dp * dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      state_redist_max_order, local_ok, transport_two, theta_two, &
      selected_target)
    if (.not. local_ok) return
    call advance_owned_reactive_eb_patch_set_chemistry_2d( &
      species, reactions, 0.5_dp * dt, rtol, atol, distribution, &
      candidate_state, candidate_temperature, coarse_geometry, &
      candidate_set, local_ok, chemistry_two)
    if (.not. local_ok) return

    theta = min(theta_one, theta_two)
    coarse_state = candidate_state
    coarse_temperature = candidate_temperature
    patch_set = candidate_set
    ok = .true.
    if (present(local_chemistry_advances)) &
      local_chemistry_advances = chemistry_one + chemistry_two
    if (present(local_hydro_advances)) &
      local_hydro_advances = hydro_advances
    if (present(local_transport_euler_advances)) &
      local_transport_euler_advances = transport_one + transport_two
    if (present(minimum_transport_theta)) minimum_transport_theta = theta
  end subroutine advance_owned_reactive_eb_patch_set_strang_2d

  subroutine all_ranks_accept_eb_2d( &
      distribution, local_acceptance, accepted, mpi_ok)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    logical, intent(in) :: local_acceptance
    logical, intent(out) :: accepted, mpi_ok

    integer :: ierr

    call MPI_Allreduce( &
      local_acceptance, accepted, 1, MPI_LOGICAL, MPI_LAND, &
      distribution%comm, ierr)
    mpi_ok = ierr == MPI_SUCCESS
    if (.not. mpi_ok) accepted = .false.
  end subroutine all_ranks_accept_eb_2d

  subroutine replicated_reactive_eb_patch_set_matches_2d( &
      coarse_geometry, patch_set, comm, ok)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok

    integer, allocatable :: metadata(:), minimum(:), maximum(:)
    real(dp), allocatable :: geometry(:), minimum_geometry(:)
    real(dp), allocatable :: maximum_geometry(:)
    logical :: global_ok, local_ok
    integer :: child, count_max, count_min, ierr, index, patch_count

    ok = .false.
    local_ok = coarse_geometry%is_valid() .and. allocated(patch_set%children)
    if (local_ok) then
      do child = 1, patch_set%patch_count()
        local_ok = local_ok .and. &
          patch_set%children(child)%geometry%is_valid() .and. &
          patch_set%children(child)%patch%is_valid( &
            coarse_geometry, patch_set%children(child)%geometry) .and. &
          patch_set%children(child)%patch%refinement_ratio >= 2
      end do
    end if
    call MPI_Allreduce( &
      local_ok, global_ok, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
    if (ierr /= MPI_SUCCESS .or. .not. global_ok) return
    patch_count = patch_set%patch_count()
    call MPI_Allreduce( &
      patch_count, count_min, 1, MPI_INTEGER, MPI_MIN, comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      patch_count, count_max, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
    if (ierr /= MPI_SUCCESS .or. count_min /= count_max) return

    allocate(metadata(2 + 7 * patch_count))
    metadata(1:2) = [coarse_geometry%nx, coarse_geometry%ny]
    index = 3
    do child = 1, patch_count
      metadata(index:index + 6) = [ &
        patch_set%children(child)%patch%coarse_i_lower, &
        patch_set%children(child)%patch%coarse_i_upper, &
        patch_set%children(child)%patch%coarse_j_lower, &
        patch_set%children(child)%patch%coarse_j_upper, &
        patch_set%children(child)%patch%refinement_ratio, &
        patch_set%children(child)%geometry%nx, &
        patch_set%children(child)%geometry%ny]
      index = index + 7
    end do
    allocate(minimum(size(metadata)), maximum(size(metadata)))
    call MPI_Allreduce( &
      metadata, minimum, size(metadata), MPI_INTEGER, MPI_MIN, comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      metadata, maximum, size(metadata), MPI_INTEGER, MPI_MAX, comm, ierr)
    if (ierr /= MPI_SUCCESS .or. any(minimum /= maximum)) return

    allocate(geometry(7 * (patch_count + 1)))
    geometry(1:7) = [ &
      coarse_geometry%x_lower, coarse_geometry%x_upper, &
      coarse_geometry%y_lower, coarse_geometry%y_upper, &
      coarse_geometry%dx, coarse_geometry%dy, &
      sum(coarse_geometry%volume_fraction)]
    index = 8
    do child = 1, patch_count
      geometry(index:index + 6) = [ &
        patch_set%children(child)%geometry%x_lower, &
        patch_set%children(child)%geometry%x_upper, &
        patch_set%children(child)%geometry%y_lower, &
        patch_set%children(child)%geometry%y_upper, &
        patch_set%children(child)%geometry%dx, &
        patch_set%children(child)%geometry%dy, &
        sum(patch_set%children(child)%geometry%volume_fraction)]
      index = index + 7
    end do
    if (any(.not. ieee_is_finite(geometry))) return
    allocate(minimum_geometry(size(geometry)), maximum_geometry(size(geometry)))
    call MPI_Allreduce( &
      geometry, minimum_geometry, size(geometry), MPI_DOUBLE_PRECISION, &
      MPI_MIN, comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      geometry, maximum_geometry, size(geometry), MPI_DOUBLE_PRECISION, &
      MPI_MAX, comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    ok = all(minimum_geometry == maximum_geometry)
  end subroutine replicated_reactive_eb_patch_set_matches_2d

end module mpi_amr_eb_patch_2d_mod
