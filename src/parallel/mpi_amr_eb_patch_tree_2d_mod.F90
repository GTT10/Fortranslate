module mpi_amr_eb_patch_tree_2d_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mpi_f08
  use precision_mod, only: dp
  use state_indices_mod, only: irho, iet
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: &
    gas_transport_species, compatible_transport_database
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component
  use reactive_2d_mod, only: advance_reactive_chemistry_2d
  use reactive_boundary_2d_mod, only: &
    reactive_boundary_set_2d, validate_reactive_boundary_set_2d
  use reactive_eb_cfl_2d_mod, only: compute_reactive_eb_cfl_timestep_2d
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_covered_cell
  use eb_reactive_reconstruction_2d_mod, only: &
    reactive_eb_exterior_state_2d
  use eb_reactive_redistribution_2d_mod, only: &
    advance_reactive_eb_state_redistributed_2d
  use eb_reactive_transport_2d_mod, only: &
    reactive_eb_transport_timestep_2d, reactive_eb_transport_fluxes_rhs_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, average_down_reactive_eb_state_patch_2d
  use amr_eb_flux_register_2d_mod, only: &
    amr_eb_flux_register_2d, initialize_amr_eb_flux_register_2d, &
    accumulate_coarse_eb_fluxes_2d, accumulate_fine_eb_fluxes_2d, &
    reflux_reactive_eb_state_patch_2d
  use amr_eb_reactive_2d_mod, only: &
    reactive_eb_patch_exterior_context_2d, &
    extract_reactive_eb_patch_exterior_context_2d, &
    build_reactive_eb_patch_exterior_from_context_2d, &
    advance_reactive_eb_level_2d, prolong_reactive_eb_patch_pcm_2d
  use amr_eb_transport_2d_mod, only: recover_transport_temperature_2d
  use amr_eb_regrid_2d_mod, only: &
    amr_eb_tagging_criteria_2d, amr_eb_regrid_plan_collection_2d, &
    plan_reactive_eb_temperature_regrid_collection_2d
  use amr_eb_patch_tree_2d_mod, only: &
    amr_eb_patch_tree_child_plan_2d, &
    amr_eb_patch_tree_level_plan_2d, amr_eb_patch_tree_topology_2d, &
    initialize_amr_eb_patch_tree_topology_2d, &
    patch_tree_topologies_match_2d
  use amr_eb_patch_tree_reactive_2d_mod, only: &
    reactive_amr_eb_patch_tree_2d
  implicit none
  private

  integer, parameter :: sparse_tree_state_tag = 27101
  integer, parameter :: sparse_tree_temperature_tag = 27102
  integer, parameter :: sparse_tree_restriction_tag = 27103
  integer, parameter :: sparse_tree_hydro_context_tag = 27104
  integer, parameter :: sparse_tree_hydro_flux_tag = 27105
  integer, parameter :: sparse_tree_hydro_reflux_tag = 27106
  integer, parameter :: sparse_tree_hydro_correction_tag = 27107
  integer, parameter :: sparse_tree_hydro_average_tag = 27108
  integer, parameter :: sparse_tree_regrid_root_tag = 27109
  integer, parameter :: sparse_tree_regrid_prolongation_tag = 27110
  integer, parameter :: sparse_tree_regrid_overlap_tag = 27111
  integer, parameter :: sparse_tree_gather_state_tag = 27112
  integer, parameter :: sparse_tree_gather_temperature_tag = 27113
  integer, parameter :: sparse_tree_scatter_state_tag = 27114
  integer, parameter :: sparse_tree_scatter_temperature_tag = 27115
  real(dp), parameter :: sparse_regrid_geometry_tolerance = &
    5.0e3_dp * epsilon(1.0_dp)
  real(dp), parameter :: sparse_regrid_conservation_tolerance = &
    5.0e4_dp * epsilon(1.0_dp)

  type, public :: mpi_amr_eb_patch_tree_level_ownership_2d
    integer, allocatable :: owners(:)
    integer, allocatable :: cell_counts(:)
    integer(int64), allocatable :: work_counts(:)
  contains
    procedure :: patch_count => mpi_amr_eb_tree_level_patch_count
    procedure :: is_valid => mpi_amr_eb_tree_level_ownership_is_valid
  end type mpi_amr_eb_patch_tree_level_ownership_2d

  type, public :: mpi_amr_eb_patch_tree_distribution_2d
    type(MPI_Comm) :: comm = MPI_COMM_NULL
    integer :: rank = -1
    integer :: nranks = 0
    integer :: subcycle_exponent = 0
    type(mpi_amr_eb_patch_tree_level_ownership_2d), allocatable :: levels(:)
    integer, allocatable :: rank_cell_counts(:)
    integer, allocatable :: rank_patch_counts(:)
    integer(int64), allocatable :: rank_work_counts(:)
  contains
    procedure :: level_count => mpi_amr_eb_tree_distribution_level_count
    procedure :: owner_of => mpi_amr_eb_tree_distribution_owner_of
    procedure :: is_local => mpi_amr_eb_tree_distribution_is_local
    procedure :: is_valid => mpi_amr_eb_tree_distribution_is_valid
  end type mpi_amr_eb_patch_tree_distribution_2d

  type, public :: mpi_sparse_reactive_amr_eb_patch_tree_node_2d
    real(dp), allocatable :: state(:, :, :)
    real(dp), allocatable :: temperature(:, :)
  contains
    procedure :: has_data => mpi_sparse_reactive_amr_eb_node_has_data
  end type mpi_sparse_reactive_amr_eb_patch_tree_node_2d

  type, public :: mpi_sparse_reactive_amr_eb_patch_tree_level_2d
    type(mpi_sparse_reactive_amr_eb_patch_tree_node_2d), &
      allocatable :: patches(:)
  contains
    procedure :: patch_count => mpi_sparse_reactive_amr_eb_level_patch_count
  end type mpi_sparse_reactive_amr_eb_patch_tree_level_2d

  type, public :: mpi_sparse_reactive_amr_eb_patch_tree_2d
    integer :: nvar = 0
    type(amr_eb_patch_tree_topology_2d) :: topology
    type(mpi_sparse_reactive_amr_eb_patch_tree_level_2d), &
      allocatable :: levels(:)
  contains
    procedure :: level_count => mpi_sparse_reactive_amr_eb_level_count
    procedure :: is_valid => mpi_sparse_reactive_amr_eb_tree_is_valid
  end type mpi_sparse_reactive_amr_eb_patch_tree_2d

  abstract interface
    subroutine sparse_reactive_amr_eb_tree_geometry_builder_2d( &
        parent_geometry, coarse_i_lower, coarse_i_upper, coarse_j_lower, &
        coarse_j_upper, refinement_ratio, child_geometry, ok)
      import :: eb_geometry_2d
      type(eb_geometry_2d), intent(in) :: parent_geometry
      integer, intent(in) :: coarse_i_lower, coarse_i_upper
      integer, intent(in) :: coarse_j_lower, coarse_j_upper
      integer, intent(in) :: refinement_ratio
      type(eb_geometry_2d), intent(out) :: child_geometry
      logical, intent(out) :: ok
    end subroutine sparse_reactive_amr_eb_tree_geometry_builder_2d
  end interface

  public :: initialize_mpi_amr_eb_patch_tree_distribution_2d
  public :: mpi_amr_eb_patch_tree_distribution_matches_2d
  public :: synchronize_owned_reactive_amr_eb_patch_tree_2d
  public :: initialize_sparse_owned_reactive_amr_eb_patch_tree_2d
  public :: initialize_sparse_owned_reactive_amr_eb_patch_tree_root_2d
  public :: materialize_sparse_owned_reactive_amr_eb_patch_tree_2d
  public :: gather_sparse_owned_reactive_amr_eb_patch_tree_to_root_2d
  public :: scatter_root_reactive_amr_eb_patch_tree_to_sparse_2d
  public :: migrate_sparse_owned_reactive_amr_eb_patch_tree_2d
  public :: plan_tagged_sparse_owned_reactive_amr_eb_patch_tree_2d
  public :: regrid_sparse_owned_reactive_amr_eb_patch_tree_2d
  public :: regrid_tagged_sparse_owned_reactive_amr_eb_patch_tree_2d
  public :: compute_sparse_owned_reactive_amr_eb_patch_tree_timestep_2d
  public :: advance_sparse_owned_reactive_amr_eb_patch_tree_chemistry_2d
  public :: advance_sparse_owned_reactive_amr_eb_patch_tree_hydro_2d
  public :: advance_sparse_owned_reactive_amr_eb_patch_tree_transport_2d
  public :: advance_sparse_owned_reactive_amr_eb_patch_tree_full_physics_2d
  public :: advance_sparse_owned_reactive_amr_eb_patch_tree_to_time_2d
  public :: composite_sparse_amr_eb_patch_tree_integral_2d
  public :: composite_sparse_amr_eb_patch_subtree_integral_2d

contains

  pure integer function mpi_amr_eb_tree_level_patch_count(self) result(count)
    class(mpi_amr_eb_patch_tree_level_ownership_2d), intent(in) :: self

    count = 0
    if (allocated(self%owners)) count = size(self%owners)
  end function mpi_amr_eb_tree_level_patch_count

  pure logical function mpi_amr_eb_tree_level_ownership_is_valid( &
      self, nranks) result(valid)
    class(mpi_amr_eb_patch_tree_level_ownership_2d), intent(in) :: self
    integer, intent(in) :: nranks

    valid = nranks >= 1 .and. allocated(self%owners) .and. &
      allocated(self%cell_counts) .and. allocated(self%work_counts)
    if (.not. valid) return
    valid = size(self%owners) >= 1 .and. &
      size(self%cell_counts) == size(self%owners) .and. &
      size(self%work_counts) == size(self%owners) .and. &
      all(self%owners >= 0) .and. all(self%owners < nranks) .and. &
      all(self%cell_counts >= 1) .and. all(self%work_counts >= 1_int64)
  end function mpi_amr_eb_tree_level_ownership_is_valid

  pure integer function mpi_amr_eb_tree_distribution_level_count(self) &
      result(count)
    class(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: self

    count = 0
    if (allocated(self%levels)) count = size(self%levels)
  end function mpi_amr_eb_tree_distribution_level_count

  pure integer function mpi_amr_eb_tree_distribution_owner_of( &
      self, level, patch) result(owner)
    class(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: self
    integer, intent(in) :: level, patch

    owner = -1
    if (.not. allocated(self%levels)) return
    if (level < 0 .or. level >= size(self%levels)) return
    if (.not. allocated(self%levels(level + 1)%owners)) return
    if (patch < 1 .or. &
        patch > size(self%levels(level + 1)%owners)) return
    owner = self%levels(level + 1)%owners(patch)
  end function mpi_amr_eb_tree_distribution_owner_of

  pure logical function mpi_amr_eb_tree_distribution_is_local( &
      self, level, patch) result(local)
    class(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: self
    integer, intent(in) :: level, patch

    local = self%owner_of(level, patch) == self%rank
  end function mpi_amr_eb_tree_distribution_is_local

  pure logical function mpi_amr_eb_tree_distribution_is_valid(self) &
      result(valid)
    class(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: self

    integer, allocatable :: cells(:), patches(:)
    integer(int64), allocatable :: work(:)
    integer :: level, owner, patch

    valid = self%rank >= 0 .and. self%nranks >= 1 .and. &
      self%rank < self%nranks .and. allocated(self%levels) .and. &
      self%subcycle_exponent >= 0 .and. self%subcycle_exponent <= 2 .and. &
      allocated(self%rank_cell_counts) .and. &
      allocated(self%rank_patch_counts) .and. &
      allocated(self%rank_work_counts)
    if (.not. valid) return
    valid = size(self%levels) >= 1 .and. &
      size(self%rank_cell_counts) == self%nranks .and. &
      size(self%rank_patch_counts) == self%nranks .and. &
      size(self%rank_work_counts) == self%nranks
    if (.not. valid) return

    allocate(cells(self%nranks), patches(self%nranks), work(self%nranks))
    cells = 0
    patches = 0
    work = 0_int64
    do level = 1, size(self%levels)
      valid = self%levels(level)%is_valid(self%nranks)
      if (.not. valid) return
      do patch = 1, self%levels(level)%patch_count()
        owner = self%levels(level)%owners(patch) + 1
        if (cells(owner) > huge(cells(owner)) - &
            self%levels(level)%cell_counts(patch)) then
          valid = .false.
          return
        end if
        if (work(owner) > huge(work(owner)) - &
            self%levels(level)%work_counts(patch)) then
          valid = .false.
          return
        end if
        cells(owner) = cells(owner) + self%levels(level)%cell_counts(patch)
        patches(owner) = patches(owner) + 1
        work(owner) = work(owner) + self%levels(level)%work_counts(patch)
      end do
    end do
    valid = all(cells == self%rank_cell_counts) .and. &
      all(patches == self%rank_patch_counts) .and. &
      all(work == self%rank_work_counts) .and. &
      sum(self%rank_patch_counts) >= 1
  end function mpi_amr_eb_tree_distribution_is_valid

  pure logical function mpi_sparse_reactive_amr_eb_node_has_data(self) &
      result(has_data)
    class(mpi_sparse_reactive_amr_eb_patch_tree_node_2d), intent(in) :: self

    has_data = allocated(self%state) .and. allocated(self%temperature)
  end function mpi_sparse_reactive_amr_eb_node_has_data

  pure integer function mpi_sparse_reactive_amr_eb_level_patch_count(self) &
      result(count)
    class(mpi_sparse_reactive_amr_eb_patch_tree_level_2d), intent(in) :: self

    count = 0
    if (allocated(self%patches)) count = size(self%patches)
  end function mpi_sparse_reactive_amr_eb_level_patch_count

  pure integer function mpi_sparse_reactive_amr_eb_level_count(self) &
      result(count)
    class(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(in) :: self

    count = 0
    if (allocated(self%levels)) count = size(self%levels)
  end function mpi_sparse_reactive_amr_eb_level_count

  logical function mpi_sparse_reactive_amr_eb_tree_is_valid( &
      self, distribution) result(valid)
    class(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(in) :: self
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution

    type(eb_geometry_2d) :: geometry
    integer :: level, patch
    logical :: geometry_ok, local

    valid = self%nvar >= 1 .and. allocated(self%levels) .and. &
      self%topology%is_valid() .and. &
      mpi_amr_eb_patch_tree_distribution_matches_2d( &
        distribution, self%topology)
    if (.not. valid) return
    valid = size(self%levels) == self%topology%level_count()
    if (.not. valid) return
    do level = 1, size(self%levels)
      valid = allocated(self%levels(level)%patches) .and. &
        self%levels(level)%patch_count() == &
          self%topology%level_patch_count(level - 1)
      if (.not. valid) return
      do patch = 1, self%levels(level)%patch_count()
        local = distribution%is_local(level - 1, patch)
        valid = allocated(self%levels(level)%patches(patch)%state) .eqv. &
          allocated(self%levels(level)%patches(patch)%temperature)
        if (.not. valid) return
        valid = self%levels(level)%patches(patch)%has_data() .eqv. local
        if (.not. valid) return
        if (.not. local) cycle
        call topology_patch_geometry_2d( &
          self%topology, level - 1, patch, geometry, geometry_ok)
        if (.not. geometry_ok) then
          valid = .false.
          return
        end if
        valid = all(shape(self%levels(level)%patches(patch)%state) == &
            [self%nvar, geometry%nx, geometry%ny]) .and. &
          all(shape(self%levels(level)%patches(patch)%temperature) == &
            [geometry%nx, geometry%ny]) .and. &
          all(ieee_is_finite( &
            self%levels(level)%patches(patch)%state)) .and. &
          all(ieee_is_finite( &
            self%levels(level)%patches(patch)%temperature)) .and. &
          all(self%levels(level)%patches(patch)%temperature > 0.0_dp)
        if (.not. valid) return
      end do
    end do
  end function mpi_sparse_reactive_amr_eb_tree_is_valid

  subroutine initialize_mpi_amr_eb_patch_tree_distribution_2d( &
      topology, comm, distribution, ok, subcycle_exponent)
    type(amr_eb_patch_tree_topology_2d), intent(in) :: topology
    type(MPI_Comm), intent(in) :: comm
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(out) :: distribution
    logical, intent(out) :: ok
    integer, intent(in), optional :: subcycle_exponent

    type(eb_geometry_2d) :: geometry
    integer(int64) :: level_scale, patch_work
    integer :: cell_count, exponent, exponent_max, exponent_min
    integer :: ierr, level, owner, patch, power, ratio
    logical :: local_ok

    ok = .false.
    distribution%comm = comm
    call MPI_Comm_rank(comm, distribution%rank, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Comm_size(comm, distribution%nranks, ierr)
    if (ierr /= MPI_SUCCESS .or. distribution%nranks < 1) return
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

    call replicated_topology_matches_2d(topology, comm, local_ok)
    if (.not. local_ok) return
    allocate(distribution%levels(topology%level_count()))
    allocate(distribution%rank_cell_counts(distribution%nranks), source=0)
    allocate(distribution%rank_patch_counts(distribution%nranks), source=0)
    allocate(distribution%rank_work_counts(distribution%nranks), &
      source=0_int64)
    level_scale = 1_int64

    do level = 0, topology%level_count() - 1
      if (level > 0) then
        ratio = topology%relations(level)%refinement_ratio
        do power = 1, exponent
          if (level_scale > huge(level_scale) / int(ratio, int64)) return
          level_scale = level_scale * int(ratio, int64)
        end do
      end if
      allocate(distribution%levels(level + 1)%owners( &
        topology%level_patch_count(level)))
      allocate(distribution%levels(level + 1)%cell_counts( &
        topology%level_patch_count(level)))
      allocate(distribution%levels(level + 1)%work_counts( &
        topology%level_patch_count(level)))
      do patch = 1, topology%level_patch_count(level)
        call topology_patch_geometry_2d( &
          topology, level, patch, geometry, local_ok)
        if (.not. local_ok .or. &
            geometry%nx > huge(cell_count) / geometry%ny) return
        cell_count = geometry%nx * geometry%ny
        if (int(cell_count, int64) > &
            huge(patch_work) / level_scale) return
        patch_work = int(cell_count, int64) * level_scale
        owner = minloc(distribution%rank_work_counts, dim=1)
        distribution%levels(level + 1)%owners(patch) = owner - 1
        distribution%levels(level + 1)%cell_counts(patch) = cell_count
        distribution%levels(level + 1)%work_counts(patch) = patch_work
        if (distribution%rank_cell_counts(owner) > &
            huge(cell_count) - cell_count .or. &
            distribution%rank_work_counts(owner) > &
            huge(patch_work) - patch_work) return
        distribution%rank_cell_counts(owner) = &
          distribution%rank_cell_counts(owner) + cell_count
        distribution%rank_patch_counts(owner) = &
          distribution%rank_patch_counts(owner) + 1
        distribution%rank_work_counts(owner) = &
          distribution%rank_work_counts(owner) + patch_work
      end do
    end do
    ok = distribution%is_valid() .and. &
      mpi_amr_eb_patch_tree_distribution_matches_2d(distribution, topology)
  end subroutine initialize_mpi_amr_eb_patch_tree_distribution_2d

  logical function mpi_amr_eb_patch_tree_distribution_matches_2d( &
      distribution, topology) result(matches)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(amr_eb_patch_tree_topology_2d), intent(in) :: topology

    integer(int64) :: expected_work, level_scale
    integer :: cell_count, exponent, level, patch, power, ratio
    type(eb_geometry_2d) :: geometry
    logical :: geometry_ok

    matches = distribution%is_valid() .and. topology%is_valid()
    if (.not. matches) return
    matches = distribution%level_count() == topology%level_count()
    if (.not. matches) return
    exponent = distribution%subcycle_exponent
    level_scale = 1_int64
    do level = 0, topology%level_count() - 1
      if (level > 0) then
        ratio = topology%relations(level)%refinement_ratio
        do power = 1, exponent
          if (level_scale > huge(level_scale) / int(ratio, int64)) then
            matches = .false.
            return
          end if
          level_scale = level_scale * int(ratio, int64)
        end do
      end if
      matches = distribution%levels(level + 1)%patch_count() == &
        topology%level_patch_count(level)
      if (.not. matches) return
      do patch = 1, topology%level_patch_count(level)
        call topology_patch_geometry_2d( &
          topology, level, patch, geometry, geometry_ok)
        if (.not. geometry_ok .or. &
            geometry%nx > huge(cell_count) / geometry%ny) then
          matches = .false.
          return
        end if
        cell_count = geometry%nx * geometry%ny
        matches = distribution%levels(level + 1)%cell_counts(patch) == &
          cell_count
        if (.not. matches .or. int(cell_count, int64) > &
            huge(expected_work) / level_scale) then
          matches = .false.
          return
        end if
        expected_work = int(cell_count, int64) * level_scale
        matches = distribution%levels(level + 1)%work_counts(patch) == &
          expected_work
        if (.not. matches) return
      end do
    end do
  end function mpi_amr_eb_patch_tree_distribution_matches_2d

  subroutine synchronize_owned_reactive_amr_eb_patch_tree_2d( &
      distribution, solution, ok, local_entity_publications)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(reactive_amr_eb_patch_tree_2d), intent(inout) :: solution
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_entity_publications

    type(reactive_amr_eb_patch_tree_2d) :: candidate
    integer :: ierr, level, owner, patch, publications
    logical :: accepted, global_ok, local_ok

    ok = .false.
    publications = 0
    if (present(local_entity_publications)) local_entity_publications = 0
    local_ok = solution%is_valid() .and. &
      mpi_amr_eb_patch_tree_distribution_matches_2d( &
        distribution, solution%topology)
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    candidate = solution
    do level = 1, candidate%level_count()
      do patch = 1, candidate%levels(level)%patch_count()
        owner = distribution%owner_of(level - 1, patch)
        if (distribution%rank == owner) publications = publications + 1
        call MPI_Bcast( &
          candidate%levels(level)%patches(patch)%state, &
          size(candidate%levels(level)%patches(patch)%state), &
          MPI_DOUBLE_PRECISION, owner, distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS) return
        call MPI_Bcast( &
          candidate%levels(level)%patches(patch)%temperature, &
          size(candidate%levels(level)%patches(patch)%temperature), &
          MPI_DOUBLE_PRECISION, owner, distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS) return
      end do
    end do
    local_ok = candidate%is_valid()
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    solution = candidate
    ok = .true.
    if (present(local_entity_publications)) &
      local_entity_publications = publications
  end subroutine synchronize_owned_reactive_amr_eb_patch_tree_2d

  subroutine initialize_sparse_owned_reactive_amr_eb_patch_tree_2d( &
      distribution, replicated, sparse, ok, local_allocated_cells)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(reactive_amr_eb_patch_tree_2d), intent(in) :: replicated
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(out) :: sparse
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_allocated_cells

    type(mpi_sparse_reactive_amr_eb_patch_tree_2d) :: candidate
    integer :: level, patch
    logical :: accepted, global_ok, local_ok

    sparse = mpi_sparse_reactive_amr_eb_patch_tree_2d()
    ok = .false.
    if (present(local_allocated_cells)) local_allocated_cells = 0
    call replicated_distribution_matches_2d( &
      distribution, replicated%topology, local_ok)
    local_ok = local_ok .and. replicated%is_valid() .and. &
      mpi_amr_eb_patch_tree_distribution_matches_2d( &
        distribution, replicated%topology)
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    call allocate_sparse_tree_layout_2d( &
      distribution, replicated%nvar, replicated%topology, candidate, local_ok)
    if (.not. local_ok) return
    do level = 1, candidate%level_count()
      do patch = 1, candidate%levels(level)%patch_count()
        if (.not. distribution%is_local(level - 1, patch)) cycle
        candidate%levels(level)%patches(patch)%state = &
          replicated%levels(level)%patches(patch)%state
        candidate%levels(level)%patches(patch)%temperature = &
          replicated%levels(level)%patches(patch)%temperature
      end do
    end do
    local_ok = candidate%is_valid(distribution)
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    sparse = candidate
    ok = .true.
    if (present(local_allocated_cells)) &
      local_allocated_cells = &
        distribution%rank_cell_counts(distribution%rank + 1)
  end subroutine initialize_sparse_owned_reactive_amr_eb_patch_tree_2d

  subroutine initialize_sparse_owned_reactive_amr_eb_patch_tree_root_2d( &
      distribution, topology, nvar, root_state, root_temperature, sparse, &
      ok, local_allocated_cells)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(amr_eb_patch_tree_topology_2d), intent(in) :: topology
    integer, intent(in) :: nvar
    real(dp), allocatable, intent(inout) :: root_state(:, :, :)
    real(dp), allocatable, intent(inout) :: root_temperature(:, :)
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(out) :: sparse
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_allocated_cells

    type(mpi_sparse_reactive_amr_eb_patch_tree_2d) :: candidate
    integer :: ierr, nvar_maximum, nvar_minimum, root_owner
    logical :: accepted, global_ok, local, local_ok

    sparse = mpi_sparse_reactive_amr_eb_patch_tree_2d()
    ok = .false.
    if (present(local_allocated_cells)) local_allocated_cells = 0
    call replicated_distribution_matches_2d( &
      distribution, topology, local_ok)
    if (local_ok) local_ok = nvar >= 1
    if (local_ok) local_ok = topology%level_count() == 1
    if (local_ok) local_ok = topology%level_patch_count(0) == 1
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    call MPI_Allreduce( &
      nvar, nvar_minimum, 1, MPI_INTEGER, MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      nvar, nvar_maximum, 1, MPI_INTEGER, MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    if (nvar_minimum /= nvar_maximum) return

    root_owner = distribution%owner_of(0, 1)
    local = distribution%rank == root_owner
    if (local) then
      local_ok = allocated(root_state) .and. allocated(root_temperature)
      if (local_ok) then
        local_ok = all(shape(root_state) == [ &
            nvar, topology%root_geometry%nx, topology%root_geometry%ny]) &
          .and. all(shape(root_temperature) == [ &
            topology%root_geometry%nx, topology%root_geometry%ny])
      end if
      if (local_ok) then
        local_ok = all(ieee_is_finite(root_state)) .and. &
          all(ieee_is_finite(root_temperature)) .and. &
          all(root_temperature > 0.0_dp)
      end if
    else
      local_ok = .not. allocated(root_state) .and. &
        .not. allocated(root_temperature)
    end if
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    candidate%nvar = nvar
    candidate%topology = topology
    allocate(candidate%levels(1))
    allocate(candidate%levels(1)%patches(1))
    if (local) then
      call move_alloc( &
        root_state, candidate%levels(1)%patches(1)%state)
      call move_alloc( &
        root_temperature, candidate%levels(1)%patches(1)%temperature)
    end if
    local_ok = candidate%is_valid(distribution)
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) then
      if (local) then
        if (allocated(candidate%levels(1)%patches(1)%state)) &
          call move_alloc( &
            candidate%levels(1)%patches(1)%state, root_state)
        if (allocated(candidate%levels(1)%patches(1)%temperature)) &
          call move_alloc( &
            candidate%levels(1)%patches(1)%temperature, root_temperature)
      end if
      return
    end if

    sparse%nvar = candidate%nvar
    sparse%topology = candidate%topology
    call move_alloc(candidate%levels, sparse%levels)
    ok = .true.
    if (present(local_allocated_cells)) &
      local_allocated_cells = &
        distribution%rank_cell_counts(distribution%rank + 1)
  end subroutine &
    initialize_sparse_owned_reactive_amr_eb_patch_tree_root_2d

  subroutine materialize_sparse_owned_reactive_amr_eb_patch_tree_2d( &
      distribution, sparse, replicated, ok, local_entity_publications)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(in) :: sparse
    type(reactive_amr_eb_patch_tree_2d), intent(out) :: replicated
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_entity_publications

    type(reactive_amr_eb_patch_tree_2d) :: candidate
    type(eb_geometry_2d) :: geometry
    integer :: ierr, level, owner, patch, publications
    logical :: accepted, geometry_ok, global_ok, local_ok

    replicated = reactive_amr_eb_patch_tree_2d()
    ok = .false.
    publications = 0
    if (present(local_entity_publications)) local_entity_publications = 0
    call replicated_distribution_matches_2d( &
      distribution, sparse%topology, local_ok)
    local_ok = local_ok .and. sparse%is_valid(distribution)
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    candidate%nvar = sparse%nvar
    candidate%topology = sparse%topology
    allocate(candidate%levels(sparse%level_count()))
    do level = 1, candidate%level_count()
      allocate(candidate%levels(level)%patches( &
        sparse%levels(level)%patch_count()))
      do patch = 1, candidate%levels(level)%patch_count()
        call topology_patch_geometry_2d( &
          sparse%topology, level - 1, patch, geometry, geometry_ok)
        if (.not. geometry_ok) return
        allocate(candidate%levels(level)%patches(patch)%state( &
          sparse%nvar, geometry%nx, geometry%ny))
        allocate(candidate%levels(level)%patches(patch)%temperature( &
          geometry%nx, geometry%ny))
        owner = distribution%owner_of(level - 1, patch)
        if (distribution%rank == owner) then
          candidate%levels(level)%patches(patch)%state = &
            sparse%levels(level)%patches(patch)%state
          candidate%levels(level)%patches(patch)%temperature = &
            sparse%levels(level)%patches(patch)%temperature
          publications = publications + 1
        end if
        call MPI_Bcast( &
          candidate%levels(level)%patches(patch)%state, &
          size(candidate%levels(level)%patches(patch)%state), &
          MPI_DOUBLE_PRECISION, owner, distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS) return
        call MPI_Bcast( &
          candidate%levels(level)%patches(patch)%temperature, &
          size(candidate%levels(level)%patches(patch)%temperature), &
          MPI_DOUBLE_PRECISION, owner, distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS) return
      end do
    end do
    local_ok = candidate%is_valid()
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    replicated = candidate
    ok = .true.
    if (present(local_entity_publications)) &
      local_entity_publications = publications
  end subroutine materialize_sparse_owned_reactive_amr_eb_patch_tree_2d

  subroutine gather_sparse_owned_reactive_amr_eb_patch_tree_to_root_2d( &
      distribution, sparse, root, replicated, ok, local_entity_transfers)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(in) :: sparse
    integer, intent(in) :: root
    type(reactive_amr_eb_patch_tree_2d), intent(out) :: replicated
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_entity_transfers

    type(reactive_amr_eb_patch_tree_2d) :: candidate
    type(eb_geometry_2d) :: geometry
    type(MPI_Status) :: status
    logical :: accepted, geometry_ok, global_ok, local_ok
    integer :: ierr, level, owner, patch, root_maximum, root_minimum
    integer :: transfers

    replicated = reactive_amr_eb_patch_tree_2d()
    ok = .false.
    transfers = 0
    if (present(local_entity_transfers)) local_entity_transfers = 0
    call replicated_distribution_matches_2d( &
      distribution, sparse%topology, local_ok)
    local_ok = local_ok .and. sparse%is_valid(distribution) .and. &
      root >= 0 .and. root < distribution%nranks
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Allreduce( &
      root, root_minimum, 1, MPI_INTEGER, MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      root, root_maximum, 1, MPI_INTEGER, MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. root_minimum /= root_maximum) return

    if (distribution%rank == root) then
      candidate%nvar = sparse%nvar
      candidate%topology = sparse%topology
      allocate(candidate%levels(sparse%level_count()))
    end if
    do level = 1, sparse%level_count()
      if (distribution%rank == root) allocate( &
        candidate%levels(level)%patches( &
          sparse%levels(level)%patch_count()))
      do patch = 1, sparse%levels(level)%patch_count()
        call topology_patch_geometry_2d( &
          sparse%topology, level - 1, patch, geometry, geometry_ok)
        if (.not. geometry_ok) return
        if (distribution%rank == root) then
          allocate(candidate%levels(level)%patches(patch)%state( &
            sparse%nvar, geometry%nx, geometry%ny))
          allocate(candidate%levels(level)%patches(patch)%temperature( &
            geometry%nx, geometry%ny))
        end if
        owner = distribution%owner_of(level - 1, patch)
        if (owner == root) then
          if (distribution%rank == root) then
            candidate%levels(level)%patches(patch)%state = &
              sparse%levels(level)%patches(patch)%state
            candidate%levels(level)%patches(patch)%temperature = &
              sparse%levels(level)%patches(patch)%temperature
          end if
        else if (distribution%rank == owner) then
          call MPI_Send( &
            sparse%levels(level)%patches(patch)%state, &
            size(sparse%levels(level)%patches(patch)%state), &
            MPI_DOUBLE_PRECISION, root, sparse_tree_gather_state_tag, &
            distribution%comm, ierr)
          if (ierr /= MPI_SUCCESS) return
          call MPI_Send( &
            sparse%levels(level)%patches(patch)%temperature, &
            size(sparse%levels(level)%patches(patch)%temperature), &
            MPI_DOUBLE_PRECISION, root, &
            sparse_tree_gather_temperature_tag, distribution%comm, ierr)
          if (ierr /= MPI_SUCCESS) return
          transfers = transfers + 1
        else if (distribution%rank == root) then
          call MPI_Recv( &
            candidate%levels(level)%patches(patch)%state, &
            size(candidate%levels(level)%patches(patch)%state), &
            MPI_DOUBLE_PRECISION, owner, sparse_tree_gather_state_tag, &
            distribution%comm, status, ierr)
          if (ierr /= MPI_SUCCESS) return
          call MPI_Recv( &
            candidate%levels(level)%patches(patch)%temperature, &
            size(candidate%levels(level)%patches(patch)%temperature), &
            MPI_DOUBLE_PRECISION, owner, &
            sparse_tree_gather_temperature_tag, distribution%comm, &
            status, ierr)
          if (ierr /= MPI_SUCCESS) return
        end if
      end do
    end do
    local_ok = .true.
    if (distribution%rank == root) local_ok = candidate%is_valid()
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    if (distribution%rank == root) replicated = candidate
    ok = .true.
    if (present(local_entity_transfers)) local_entity_transfers = transfers
  end subroutine gather_sparse_owned_reactive_amr_eb_patch_tree_to_root_2d

  subroutine scatter_root_reactive_amr_eb_patch_tree_to_sparse_2d( &
      distribution, topology, replicated, root, sparse, ok, &
      local_entity_transfers)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(amr_eb_patch_tree_topology_2d), intent(in) :: topology
    type(reactive_amr_eb_patch_tree_2d), intent(in) :: replicated
    integer, intent(in) :: root
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(out) :: sparse
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_entity_transfers

    type(mpi_sparse_reactive_amr_eb_patch_tree_2d) :: candidate
    type(MPI_Status) :: status
    logical :: accepted, global_ok, local_ok
    integer :: ierr, level, nvar, owner, patch
    integer :: root_maximum, root_minimum, transfers

    sparse = mpi_sparse_reactive_amr_eb_patch_tree_2d()
    ok = .false.
    transfers = 0
    nvar = 0
    if (present(local_entity_transfers)) local_entity_transfers = 0
    local_ok = distribution%is_valid() .and. topology%is_valid() .and. &
      mpi_amr_eb_patch_tree_distribution_matches_2d( &
        distribution, topology) .and. root >= 0 .and. &
      root < distribution%nranks
    if (distribution%rank == root) local_ok = local_ok .and. &
      replicated%is_valid() .and. &
      patch_tree_topologies_match_2d(replicated%topology, topology)
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Allreduce( &
      root, root_minimum, 1, MPI_INTEGER, MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      root, root_maximum, 1, MPI_INTEGER, MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. root_minimum /= root_maximum) return
    if (distribution%rank == root) nvar = replicated%nvar
    call MPI_Bcast( &
      nvar, 1, MPI_INTEGER, root, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. nvar < 1) return

    call allocate_sparse_tree_layout_2d( &
      distribution, nvar, topology, candidate, local_ok)
    if (.not. local_ok) return
    do level = 1, candidate%level_count()
      do patch = 1, candidate%levels(level)%patch_count()
        owner = distribution%owner_of(level - 1, patch)
        if (distribution%rank == root) then
          if (owner == root) then
            candidate%levels(level)%patches(patch)%state = &
              replicated%levels(level)%patches(patch)%state
            candidate%levels(level)%patches(patch)%temperature = &
              replicated%levels(level)%patches(patch)%temperature
          else
            call MPI_Send( &
              replicated%levels(level)%patches(patch)%state, &
              size(replicated%levels(level)%patches(patch)%state), &
              MPI_DOUBLE_PRECISION, owner, sparse_tree_scatter_state_tag, &
              distribution%comm, ierr)
            if (ierr /= MPI_SUCCESS) return
            call MPI_Send( &
              replicated%levels(level)%patches(patch)%temperature, &
              size(replicated%levels(level)%patches(patch)%temperature), &
              MPI_DOUBLE_PRECISION, owner, &
              sparse_tree_scatter_temperature_tag, distribution%comm, ierr)
            if (ierr /= MPI_SUCCESS) return
            transfers = transfers + 1
          end if
        else if (distribution%rank == owner) then
          call MPI_Recv( &
            candidate%levels(level)%patches(patch)%state, &
            size(candidate%levels(level)%patches(patch)%state), &
            MPI_DOUBLE_PRECISION, root, sparse_tree_scatter_state_tag, &
            distribution%comm, status, ierr)
          if (ierr /= MPI_SUCCESS) return
          call MPI_Recv( &
            candidate%levels(level)%patches(patch)%temperature, &
            size(candidate%levels(level)%patches(patch)%temperature), &
            MPI_DOUBLE_PRECISION, root, &
            sparse_tree_scatter_temperature_tag, distribution%comm, &
            status, ierr)
          if (ierr /= MPI_SUCCESS) return
        end if
      end do
    end do
    local_ok = candidate%is_valid(distribution)
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    sparse = candidate
    ok = .true.
    if (present(local_entity_transfers)) local_entity_transfers = transfers
  end subroutine scatter_root_reactive_amr_eb_patch_tree_to_sparse_2d

  subroutine migrate_sparse_owned_reactive_amr_eb_patch_tree_2d( &
      old_distribution, new_distribution, sparse, ok, &
      local_entity_transfers)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: &
      old_distribution, new_distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(inout) :: sparse
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_entity_transfers

    type(mpi_sparse_reactive_amr_eb_patch_tree_2d) :: candidate
    type(MPI_Status) :: status
    integer :: comm_comparison, ierr, level, new_owner, old_owner, patch
    integer :: transfers
    logical :: accepted, global_ok, local_ok, new_matches, old_matches

    ok = .false.
    transfers = 0
    if (present(local_entity_transfers)) local_entity_transfers = 0
    call MPI_Comm_compare( &
      old_distribution%comm, new_distribution%comm, comm_comparison, ierr)
    local_ok = ierr == MPI_SUCCESS .and. &
      (comm_comparison == MPI_IDENT .or. comm_comparison == MPI_CONGRUENT) &
      .and. old_distribution%rank == new_distribution%rank .and. &
      old_distribution%nranks == new_distribution%nranks
    call replicated_distribution_matches_2d( &
      old_distribution, sparse%topology, old_matches)
    call replicated_distribution_matches_2d( &
      new_distribution, sparse%topology, new_matches)
    local_ok = local_ok .and. old_matches .and. new_matches .and. &
      sparse%is_valid(old_distribution)
    call all_ranks_accept_2d( &
      old_distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    call allocate_sparse_tree_layout_2d( &
      new_distribution, sparse%nvar, sparse%topology, candidate, local_ok)
    if (.not. local_ok) return
    do level = 1, sparse%level_count()
      do patch = 1, sparse%levels(level)%patch_count()
        old_owner = old_distribution%owner_of(level - 1, patch)
        new_owner = new_distribution%owner_of(level - 1, patch)
        if (old_owner == new_owner) then
          if (old_distribution%rank == old_owner) then
            candidate%levels(level)%patches(patch)%state = &
              sparse%levels(level)%patches(patch)%state
            candidate%levels(level)%patches(patch)%temperature = &
              sparse%levels(level)%patches(patch)%temperature
          end if
          cycle
        end if
        if (old_distribution%rank == old_owner) then
          call MPI_Send( &
            sparse%levels(level)%patches(patch)%state, &
            size(sparse%levels(level)%patches(patch)%state), &
            MPI_DOUBLE_PRECISION, new_owner, sparse_tree_state_tag, &
            old_distribution%comm, ierr)
          if (ierr /= MPI_SUCCESS) return
          call MPI_Send( &
            sparse%levels(level)%patches(patch)%temperature, &
            size(sparse%levels(level)%patches(patch)%temperature), &
            MPI_DOUBLE_PRECISION, new_owner, sparse_tree_temperature_tag, &
            old_distribution%comm, ierr)
          if (ierr /= MPI_SUCCESS) return
          transfers = transfers + 1
        else if (old_distribution%rank == new_owner) then
          call MPI_Recv( &
            candidate%levels(level)%patches(patch)%state, &
            size(candidate%levels(level)%patches(patch)%state), &
            MPI_DOUBLE_PRECISION, old_owner, sparse_tree_state_tag, &
            old_distribution%comm, status, ierr)
          if (ierr /= MPI_SUCCESS) return
          call MPI_Recv( &
            candidate%levels(level)%patches(patch)%temperature, &
            size(candidate%levels(level)%patches(patch)%temperature), &
            MPI_DOUBLE_PRECISION, old_owner, sparse_tree_temperature_tag, &
            old_distribution%comm, status, ierr)
          if (ierr /= MPI_SUCCESS) return
        end if
      end do
    end do
    local_ok = candidate%is_valid(new_distribution)
    call all_ranks_accept_2d( &
      old_distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    sparse = candidate
    ok = .true.
    if (present(local_entity_transfers)) local_entity_transfers = transfers
  end subroutine migrate_sparse_owned_reactive_amr_eb_patch_tree_2d

  subroutine regrid_tagged_sparse_owned_reactive_amr_eb_patch_tree_2d( &
      species, old_distribution, sparse, criteria, maximum_levels, &
      refinement_ratio, geometry_builder, new_distribution, ok, changed, &
      tagged_cells, transferred_cells, local_tagging_evaluations, &
      local_candidate_transfers, local_restriction_transfers, &
      local_prolongation_transfers, local_overlap_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: &
      old_distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(inout) :: sparse
    type(amr_eb_tagging_criteria_2d), intent(in) :: criteria
    integer, intent(in) :: maximum_levels, refinement_ratio
    procedure(sparse_reactive_amr_eb_tree_geometry_builder_2d) :: &
      geometry_builder
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(out) :: &
      new_distribution
    logical, intent(out) :: ok, changed
    integer, intent(out) :: tagged_cells, transferred_cells
    integer, intent(out), optional :: local_tagging_evaluations
    integer, intent(out), optional :: local_candidate_transfers
    integer, intent(out), optional :: local_restriction_transfers
    integer, intent(out), optional :: local_prolongation_transfers
    integer, intent(out), optional :: local_overlap_transfers

    type(mpi_sparse_reactive_amr_eb_patch_tree_2d) :: backup
    type(amr_eb_patch_tree_level_plan_2d), allocatable :: plans(:)
    integer :: candidate_transfers, overlap_transfers
    integer :: prolongation_transfers, regrid_restriction_transfers
    integer :: restriction_transfers
    integer :: tagging_evaluations
    logical :: local_ok

    ok = .false.
    changed = .false.
    tagged_cells = 0
    transferred_cells = 0
    tagging_evaluations = 0
    candidate_transfers = 0
    restriction_transfers = 0
    regrid_restriction_transfers = 0
    prolongation_transfers = 0
    overlap_transfers = 0
    new_distribution = old_distribution
    if (present(local_tagging_evaluations)) local_tagging_evaluations = 0
    if (present(local_candidate_transfers)) local_candidate_transfers = 0
    if (present(local_restriction_transfers)) &
      local_restriction_transfers = 0
    if (present(local_prolongation_transfers)) &
      local_prolongation_transfers = 0
    if (present(local_overlap_transfers)) local_overlap_transfers = 0
    backup = sparse

    call plan_tagged_sparse_owned_reactive_amr_eb_patch_tree_2d( &
      species, old_distribution, sparse, criteria, maximum_levels, &
      refinement_ratio, geometry_builder, plans, tagged_cells, local_ok, &
      tagging_evaluations, candidate_transfers, restriction_transfers)
    if (.not. local_ok) go to 900
    call regrid_sparse_owned_reactive_amr_eb_patch_tree_2d( &
      species, old_distribution, sparse, plans, new_distribution, local_ok, &
      changed, transferred_cells, regrid_restriction_transfers, &
      prolongation_transfers, overlap_transfers)
    if (.not. local_ok) go to 900
    restriction_transfers = &
      restriction_transfers + regrid_restriction_transfers

    ok = .true.
    if (present(local_tagging_evaluations)) &
      local_tagging_evaluations = tagging_evaluations
    if (present(local_candidate_transfers)) &
      local_candidate_transfers = candidate_transfers
    if (present(local_restriction_transfers)) &
      local_restriction_transfers = restriction_transfers
    if (present(local_prolongation_transfers)) &
      local_prolongation_transfers = prolongation_transfers
    if (present(local_overlap_transfers)) &
      local_overlap_transfers = overlap_transfers
    return

900 continue
    sparse = backup
    new_distribution = old_distribution
    changed = .false.
    tagged_cells = 0
    transferred_cells = 0
    if (present(local_tagging_evaluations)) local_tagging_evaluations = 0
    if (present(local_candidate_transfers)) local_candidate_transfers = 0
    if (present(local_restriction_transfers)) &
      local_restriction_transfers = 0
    if (present(local_prolongation_transfers)) &
      local_prolongation_transfers = 0
    if (present(local_overlap_transfers)) local_overlap_transfers = 0
  end subroutine &
    regrid_tagged_sparse_owned_reactive_amr_eb_patch_tree_2d

  subroutine plan_tagged_sparse_owned_reactive_amr_eb_patch_tree_2d( &
      species, distribution, sparse, criteria, maximum_levels, &
      refinement_ratio, geometry_builder, plans, tagged_cells, ok, &
      local_tagging_evaluations, local_candidate_transfers, &
      local_restriction_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(in) :: sparse
    type(amr_eb_tagging_criteria_2d), intent(in) :: criteria
    integer, intent(in) :: maximum_levels, refinement_ratio
    procedure(sparse_reactive_amr_eb_tree_geometry_builder_2d) :: &
      geometry_builder
    type(amr_eb_patch_tree_level_plan_2d), allocatable, intent(out) :: plans(:)
    integer, intent(out) :: tagged_cells
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_tagging_evaluations
    integer, intent(out), optional :: local_candidate_transfers
    integer, intent(out), optional :: local_restriction_transfers

    type(mpi_sparse_reactive_amr_eb_patch_tree_2d) :: source, candidate
    type(mpi_amr_eb_patch_tree_distribution_2d) :: candidate_distribution
    type(amr_eb_patch_tree_topology_2d) :: candidate_topology
    type(amr_eb_patch_tree_level_plan_2d), allocatable :: workspace(:)
    type(amr_eb_patch_tree_child_plan_2d), allocatable :: extended(:)
    type(amr_eb_regrid_plan_collection_2d) :: collection
    type(eb_geometry_2d) :: parent_geometry
    real(dp) :: numeric_maximum(3), numeric_minimum(3), numeric_values(3)
    integer, allocatable :: global_bounds(:), local_bounds(:)
    integer :: child, child_count, entry, global_header(2), ierr
    integer :: integer_maximum(7), integer_minimum(7), integer_values(7)
    integer :: local_header(2), parent, parent_count, parent_owner
    integer :: relation, relation_count, step_transfers
    integer :: candidate_transfers, restriction_transfers
    integer :: tagging_evaluations
    logical, allocatable :: tags(:, :)
    logical :: accepted, entity_ok, global_ok, local_ok, matches

    ok = .false.
    tagged_cells = 0
    tagging_evaluations = 0
    candidate_transfers = 0
    restriction_transfers = 0
    if (present(local_tagging_evaluations)) local_tagging_evaluations = 0
    if (present(local_candidate_transfers)) local_candidate_transfers = 0
    if (present(local_restriction_transfers)) &
      local_restriction_transfers = 0
    numeric_values = [ &
      criteria%relative_gradient_threshold, &
      criteria%absolute_gradient_threshold, criteria%scale_floor]
    integer_values = [ &
      size(species), maximum_levels, refinement_ratio, &
      criteria%buffer_cells, criteria%minimum_patch_cells_x, &
      criteria%minimum_patch_cells_y, criteria%maximum_patch_gap_cells]
    call replicated_distribution_matches_2d( &
      distribution, sparse%topology, matches)
    local_ok = matches .and. sparse%is_valid(distribution) .and. &
      sparse%nvar == reactive_nvar(size(species)) .and. &
      size(species) >= 1 .and. maximum_levels >= 1 .and. &
      refinement_ratio >= 2 .and. &
      criteria%is_valid( &
        sparse%topology%root_geometry%nx, sparse%topology%root_geometry%ny)
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Allreduce( &
      numeric_values, numeric_minimum, size(numeric_values), &
      MPI_DOUBLE_PRECISION, MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      numeric_values, numeric_maximum, size(numeric_values), &
      MPI_DOUBLE_PRECISION, MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_values, integer_minimum, size(integer_values), MPI_INTEGER, &
      MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_values, integer_maximum, size(integer_values), MPI_INTEGER, &
      MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. &
        any(numeric_minimum /= numeric_maximum) .or. &
        any(integer_minimum /= integer_maximum)) return

    source = sparse
    call synchronize_sparse_regrid_tree_2d( &
      species, distribution, source, restriction_transfers, local_ok)
    if (.not. local_ok) return
    candidate = source
    candidate_distribution = distribution
    allocate(workspace(maximum_levels - 1))
    relation_count = 0
    do relation = 1, maximum_levels - 1
      parent_count = candidate%topology%level_patch_count(relation - 1)
      child_count = 0
      do parent = 1, parent_count
        call topology_patch_geometry_2d( &
          candidate%topology, relation - 1, parent, parent_geometry, &
          entity_ok)
        parent_owner = candidate_distribution%owner_of(relation - 1, parent)
        entity_ok = entity_ok .and. parent_owner >= 0 .and. &
          parent_owner < distribution%nranks
        local_header = 0
        if (entity_ok .and. criteria%is_valid( &
            parent_geometry%nx, parent_geometry%ny) .and. &
            distribution%rank == parent_owner) then
          allocate(tags(parent_geometry%nx, parent_geometry%ny))
          call plan_reactive_eb_temperature_regrid_collection_2d( &
            candidate%levels(relation)%patches(parent)%temperature, &
            parent_geometry, criteria, tags, collection, entity_ok)
          deallocate(tags)
          if (entity_ok) then
            local_header = [ &
              collection%tagged_cell_count, collection%patch_count()]
            tagging_evaluations = tagging_evaluations + 1
          end if
        end if
        call all_ranks_accept_2d( &
          distribution%comm, entity_ok, accepted, global_ok)
        if (.not. global_ok .or. .not. accepted) return
        call MPI_Allreduce( &
          local_header, global_header, 2, MPI_INTEGER, MPI_SUM, &
          distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS .or. any(global_header < 0)) return
        tagged_cells = tagged_cells + global_header(1)
        if (global_header(2) == 0) cycle

        allocate(local_bounds(4 * global_header(2)), source=0)
        allocate(global_bounds(4 * global_header(2)))
        if (distribution%rank == parent_owner) then
          do child = 1, global_header(2)
            local_bounds(4 * child - 3:4 * child) = [ &
              collection%plans(child)%coarse_i_lower, &
              collection%plans(child)%coarse_i_upper, &
              collection%plans(child)%coarse_j_lower, &
              collection%plans(child)%coarse_j_upper]
          end do
        end if
        call MPI_Allreduce( &
          local_bounds, global_bounds, size(local_bounds), MPI_INTEGER, &
          MPI_SUM, distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS) return

        allocate(extended(child_count + global_header(2)))
        if (child_count > 0) &
          extended(1:child_count) = workspace(relation)%children
        do child = 1, global_header(2)
          entry = child_count + child
          extended(entry)%parent_patch = parent
          extended(entry)%coarse_i_lower = global_bounds(4 * child - 3)
          extended(entry)%coarse_i_upper = global_bounds(4 * child - 2)
          extended(entry)%coarse_j_lower = global_bounds(4 * child - 1)
          extended(entry)%coarse_j_upper = global_bounds(4 * child)
          call geometry_builder( &
            parent_geometry, extended(entry)%coarse_i_lower, &
            extended(entry)%coarse_i_upper, &
            extended(entry)%coarse_j_lower, &
            extended(entry)%coarse_j_upper, refinement_ratio, &
            extended(entry)%geometry, entity_ok)
          call all_ranks_accept_2d( &
            distribution%comm, entity_ok, accepted, global_ok)
          if (.not. global_ok .or. .not. accepted) return
        end do
        child_count = child_count + global_header(2)
        call move_alloc(extended, workspace(relation)%children)
        deallocate(local_bounds, global_bounds)
      end do
      if (child_count == 0) exit

      workspace(relation)%refinement_ratio = refinement_ratio
      relation_count = relation
      call initialize_amr_eb_patch_tree_topology_2d( &
        source%topology%root_geometry, workspace(1:relation_count), &
        candidate_topology, local_ok)
      call all_ranks_accept_2d( &
        distribution%comm, local_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) return
      call initialize_mpi_amr_eb_patch_tree_distribution_2d( &
        candidate_topology, distribution%comm, candidate_distribution, &
        local_ok, distribution%subcycle_exponent)
      if (.not. local_ok) return
      call initialize_direct_sparse_regrid_tree_2d( &
        species, distribution, candidate_distribution, source, &
        candidate_topology, candidate, step_transfers, local_ok)
      if (.not. local_ok) return
      candidate_transfers = candidate_transfers + step_transfers
    end do

    allocate(plans(relation_count))
    if (relation_count > 0) plans = workspace(1:relation_count)
    ok = .true.
    if (present(local_tagging_evaluations)) &
      local_tagging_evaluations = tagging_evaluations
    if (present(local_candidate_transfers)) &
      local_candidate_transfers = candidate_transfers
    if (present(local_restriction_transfers)) &
      local_restriction_transfers = restriction_transfers
  end subroutine &
    plan_tagged_sparse_owned_reactive_amr_eb_patch_tree_2d

  subroutine regrid_sparse_owned_reactive_amr_eb_patch_tree_2d( &
      species, old_distribution, sparse, plans, new_distribution, ok, &
      changed, transferred_cells, local_restriction_transfers, &
      local_prolongation_transfers, local_overlap_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: &
      old_distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(inout) :: sparse
    type(amr_eb_patch_tree_level_plan_2d), intent(in) :: plans(:)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(out) :: &
      new_distribution
    logical, intent(out) :: ok, changed
    integer, intent(out) :: transferred_cells
    integer, intent(out), optional :: local_restriction_transfers
    integer, intent(out), optional :: local_prolongation_transfers
    integer, intent(out), optional :: local_overlap_transfers

    type(mpi_sparse_reactive_amr_eb_patch_tree_2d) :: source, rebuilt
    type(amr_eb_patch_tree_topology_2d) :: new_topology
    real(dp), allocatable :: new_integral(:), old_integral(:)
    real(dp) :: integral_scale
    integer :: overlap_transfers, prolongation_transfers
    integer :: restriction_transfers, step_transfers
    logical :: accepted, global_ok, local_ok, matches

    ok = .false.
    changed = .false.
    transferred_cells = 0
    restriction_transfers = 0
    prolongation_transfers = 0
    overlap_transfers = 0
    new_distribution = old_distribution
    if (present(local_restriction_transfers)) &
      local_restriction_transfers = 0
    if (present(local_prolongation_transfers)) &
      local_prolongation_transfers = 0
    if (present(local_overlap_transfers)) local_overlap_transfers = 0

    call replicated_distribution_matches_2d( &
      old_distribution, sparse%topology, matches)
    local_ok = matches .and. sparse%is_valid(old_distribution) .and. &
      sparse%nvar == reactive_nvar(size(species)) .and. size(species) >= 1
    call all_ranks_accept_2d( &
      old_distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call initialize_amr_eb_patch_tree_topology_2d( &
      sparse%topology%root_geometry, plans, new_topology, local_ok)
    call all_ranks_accept_2d( &
      old_distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call replicated_topology_matches_2d( &
      new_topology, old_distribution%comm, local_ok)
    if (.not. local_ok) return

    changed = .not. patch_tree_topologies_match_2d( &
      sparse%topology, new_topology)
    if (.not. changed) then
      ok = .true.
      return
    end if

    allocate(old_integral(sparse%nvar), new_integral(sparse%nvar))
    call composite_sparse_amr_eb_patch_tree_integral_2d( &
      old_distribution, sparse, old_integral, local_ok)
    if (.not. local_ok) go to 900
    source = sparse
    call synchronize_sparse_regrid_tree_2d( &
      species, old_distribution, source, restriction_transfers, local_ok)
    if (.not. local_ok) go to 900
    call initialize_mpi_amr_eb_patch_tree_distribution_2d( &
      new_topology, old_distribution%comm, new_distribution, local_ok, &
      old_distribution%subcycle_exponent)
    if (.not. local_ok) go to 900
    call initialize_direct_sparse_regrid_tree_2d( &
      species, old_distribution, new_distribution, source, new_topology, &
      rebuilt, prolongation_transfers, local_ok)
    if (.not. local_ok) go to 900
    call transfer_direct_sparse_regrid_tree_overlaps_2d( &
      old_distribution, new_distribution, source, rebuilt, &
      transferred_cells, overlap_transfers, local_ok)
    if (.not. local_ok) go to 900
    step_transfers = 0
    call synchronize_sparse_regrid_tree_2d( &
      species, new_distribution, rebuilt, step_transfers, local_ok)
    if (.not. local_ok) go to 900
    restriction_transfers = restriction_transfers + step_transfers
    call composite_sparse_amr_eb_patch_tree_integral_2d( &
      new_distribution, rebuilt, new_integral, local_ok)
    if (.not. local_ok) go to 900
    integral_scale = max(1.0_dp, maxval(abs(old_integral)))
    local_ok = maxval(abs(new_integral - old_integral)) <= &
      sparse_regrid_conservation_tolerance * integral_scale
    call all_ranks_accept_2d( &
      old_distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) go to 900
    local_ok = rebuilt%is_valid(new_distribution)
    call all_ranks_accept_2d( &
      old_distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) go to 900

    sparse = rebuilt
    ok = .true.
    if (present(local_restriction_transfers)) &
      local_restriction_transfers = restriction_transfers
    if (present(local_prolongation_transfers)) &
      local_prolongation_transfers = prolongation_transfers
    if (present(local_overlap_transfers)) &
      local_overlap_transfers = overlap_transfers
    return

900 continue
    new_distribution = old_distribution
    changed = .false.
    transferred_cells = 0
    if (present(local_restriction_transfers)) &
      local_restriction_transfers = 0
    if (present(local_prolongation_transfers)) &
      local_prolongation_transfers = 0
    if (present(local_overlap_transfers)) local_overlap_transfers = 0
  end subroutine regrid_sparse_owned_reactive_amr_eb_patch_tree_2d

  subroutine synchronize_sparse_regrid_tree_2d( &
      species, distribution, sparse, local_transfers, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(inout) :: sparse
    integer, intent(out) :: local_transfers
    logical, intent(out) :: ok

    type(eb_geometry_2d) :: geometry
    real(dp), allocatable :: temperature_work(:, :)
    integer :: child, level, parent, patch, relation
    logical :: accepted, entity_ok, global_ok, local_ok

    ok = .false.
    local_transfers = 0
    local_ok = sparse%is_valid(distribution) .and. &
      sparse%nvar == reactive_nvar(size(species))
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    do level = 1, sparse%level_count()
      do patch = 1, sparse%levels(level)%patch_count()
        entity_ok = .true.
        if (distribution%is_local(level - 1, patch)) then
          call topology_patch_geometry_2d( &
            sparse%topology, level - 1, patch, geometry, entity_ok)
          if (entity_ok) then
            allocate(temperature_work(geometry%nx, geometry%ny))
            call recover_transport_temperature_2d( &
              species, sparse%levels(level)%patches(patch)%state, &
              sparse%levels(level)%patches(patch)%temperature, geometry, &
              temperature_work, entity_ok)
            if (entity_ok) &
              sparse%levels(level)%patches(patch)%temperature = &
                temperature_work
            deallocate(temperature_work)
          end if
        end if
        call all_ranks_accept_2d( &
          distribution%comm, entity_ok, accepted, global_ok)
        if (.not. global_ok .or. .not. accepted) return
      end do
    end do

    do relation = size(sparse%topology%relations), 1, -1
      do child = 1, sparse%topology%relations(relation)%child_patch_count()
        parent = sparse%topology%relations(relation)% &
          children(child)%parent_patch
        call average_down_sparse_patch_tree_edge_2d( &
          species, distribution, sparse, relation, parent, child, &
          local_transfers, local_ok)
        if (.not. local_ok) return
      end do
    end do
    local_ok = sparse%is_valid(distribution)
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    ok = global_ok .and. accepted
  end subroutine synchronize_sparse_regrid_tree_2d

  subroutine initialize_direct_sparse_regrid_tree_2d( &
      species, old_distribution, new_distribution, old_sparse, topology, &
      rebuilt, local_transfers, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: &
      old_distribution, new_distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(in) :: old_sparse
    type(amr_eb_patch_tree_topology_2d), intent(in) :: topology
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(out) :: rebuilt
    integer, intent(out) :: local_transfers
    logical, intent(out) :: ok

    type(eb_geometry_2d) :: root_geometry
    type(MPI_Status) :: status
    real(dp), allocatable :: payload(:), temperature_work(:, :)
    integer :: ierr, new_owner, old_owner, parent, relation, value_count
    logical :: accepted, entity_ok, global_ok, local_ok

    ok = .false.
    local_transfers = 0
    call allocate_sparse_tree_layout_2d( &
      new_distribution, old_sparse%nvar, topology, rebuilt, local_ok)
    call all_ranks_accept_2d( &
      old_distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call topology_patch_geometry_2d( &
      topology, 0, 1, root_geometry, entity_ok)
    old_owner = old_distribution%owner_of(0, 1)
    new_owner = new_distribution%owner_of(0, 1)
    entity_ok = entity_ok .and. old_owner >= 0 .and. new_owner >= 0 .and. &
      old_owner < old_distribution%nranks .and. &
      new_owner < new_distribution%nranks
    value_count = (rebuilt%nvar + 1) * root_geometry%nx * root_geometry%ny
    if (old_owner == new_owner) then
      if (old_distribution%rank == old_owner .and. entity_ok) then
        rebuilt%levels(1)%patches(1)%state = &
          old_sparse%levels(1)%patches(1)%state
        rebuilt%levels(1)%patches(1)%temperature = &
          old_sparse%levels(1)%patches(1)%temperature
      end if
    else if (old_distribution%rank == old_owner .and. entity_ok) then
      allocate(payload(value_count))
      call pack_sparse_patch_tree_node_fields_2d( &
        old_sparse%levels(1)%patches(1)%state, &
        old_sparse%levels(1)%patches(1)%temperature, payload)
      call MPI_Send( &
        payload, value_count, MPI_DOUBLE_PRECISION, new_owner, &
        sparse_tree_regrid_root_tag, old_distribution%comm, ierr)
      entity_ok = ierr == MPI_SUCCESS
      if (entity_ok) local_transfers = local_transfers + 1
    else if (old_distribution%rank == new_owner .and. entity_ok) then
      allocate(payload(value_count))
      call MPI_Recv( &
        payload, value_count, MPI_DOUBLE_PRECISION, old_owner, &
        sparse_tree_regrid_root_tag, old_distribution%comm, status, ierr)
      entity_ok = ierr == MPI_SUCCESS
      if (entity_ok) then
        rebuilt%levels(1)%patches(1)%state = reshape( &
          payload(1:rebuilt%nvar * root_geometry%nx * root_geometry%ny), &
          shape(rebuilt%levels(1)%patches(1)%state))
        rebuilt%levels(1)%patches(1)%temperature = reshape( &
          payload(rebuilt%nvar * root_geometry%nx * root_geometry%ny + 1:), &
          shape(rebuilt%levels(1)%patches(1)%temperature))
      end if
    end if
    if (old_distribution%rank == new_owner .and. entity_ok) then
      allocate(temperature_work(root_geometry%nx, root_geometry%ny))
      call recover_transport_temperature_2d( &
        species, rebuilt%levels(1)%patches(1)%state, &
        rebuilt%levels(1)%patches(1)%temperature, root_geometry, &
        temperature_work, entity_ok)
      if (entity_ok) &
        rebuilt%levels(1)%patches(1)%temperature = temperature_work
      deallocate(temperature_work)
    end if
    call all_ranks_accept_2d( &
      old_distribution%comm, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    do relation = 1, size(topology%relations)
      do parent = 1, topology%relations(relation)%parent_patch_count()
        call prolong_direct_sparse_regrid_parent_2d( &
          species, new_distribution, rebuilt, relation, parent, &
          local_transfers, local_ok)
        if (.not. local_ok) return
      end do
    end do
    local_ok = rebuilt%is_valid(new_distribution)
    call all_ranks_accept_2d( &
      old_distribution%comm, local_ok, accepted, global_ok)
    ok = global_ok .and. accepted
  end subroutine initialize_direct_sparse_regrid_tree_2d

  subroutine prolong_direct_sparse_regrid_parent_2d( &
      species, distribution, sparse, relation, parent, local_transfers, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(inout) :: sparse
    integer, intent(in) :: relation, parent
    integer, intent(inout) :: local_transfers
    logical, intent(out) :: ok

    type(eb_geometry_2d) :: child_geometry, parent_geometry
    type(MPI_Status) :: status
    real(dp), allocatable :: child_state(:, :, :), child_temperature(:, :)
    real(dp), allocatable :: payload(:)
    integer :: child, child_owner, first_child, ierr, last_child
    integer :: parent_owner, value_count
    logical :: accepted, entity_ok, global_ok

    ok = .false.
    call topology_patch_geometry_2d( &
      sparse%topology, relation - 1, parent, parent_geometry, entity_ok)
    parent_owner = distribution%owner_of(relation - 1, parent)
    entity_ok = entity_ok .and. parent_owner >= 0 .and. &
      parent_owner < distribution%nranks
    first_child = sparse%topology%relations(relation)% &
      child_offsets(parent) + 1
    last_child = sparse%topology%relations(relation)% &
      child_offsets(parent + 1)
    call all_ranks_accept_2d( &
      distribution%comm, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    do child = first_child, last_child
      child_geometry = sparse%topology%relations(relation)% &
        children(child)%geometry
      child_owner = distribution%owner_of(relation, child)
      entity_ok = child_owner >= 0 .and. child_owner < distribution%nranks
      if (distribution%rank == parent_owner .and. entity_ok) then
        allocate(child_state( &
          sparse%nvar, child_geometry%nx, child_geometry%ny))
        allocate(child_temperature(child_geometry%nx, child_geometry%ny))
        call prolong_reactive_eb_patch_pcm_2d( &
          species, sparse%levels(relation)%patches(parent)%state, &
          sparse%levels(relation)%patches(parent)%temperature, &
          parent_geometry, child_geometry, &
          sparse%topology%relations(relation)%children(child)%patch, &
          child_state, child_temperature, entity_ok)
      end if
      call all_ranks_accept_2d( &
        distribution%comm, entity_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) return

      value_count = (sparse%nvar + 1) * &
        child_geometry%nx * child_geometry%ny
      if (parent_owner == child_owner) then
        if (distribution%rank == child_owner) then
          sparse%levels(relation + 1)%patches(child)%state = child_state
          sparse%levels(relation + 1)%patches(child)%temperature = &
            child_temperature
        end if
      else if (distribution%rank == parent_owner) then
        allocate(payload(value_count))
        call pack_sparse_patch_tree_node_fields_2d( &
          child_state, child_temperature, payload)
        call MPI_Send( &
          payload, value_count, MPI_DOUBLE_PRECISION, child_owner, &
          sparse_tree_regrid_prolongation_tag, distribution%comm, ierr)
        entity_ok = ierr == MPI_SUCCESS
        if (entity_ok) local_transfers = local_transfers + 1
      else if (distribution%rank == child_owner) then
        allocate(payload(value_count))
        call MPI_Recv( &
          payload, value_count, MPI_DOUBLE_PRECISION, parent_owner, &
          sparse_tree_regrid_prolongation_tag, distribution%comm, &
          status, ierr)
        entity_ok = ierr == MPI_SUCCESS
        if (entity_ok) then
          sparse%levels(relation + 1)%patches(child)%state = reshape( &
            payload(1:sparse%nvar * child_geometry%nx * &
              child_geometry%ny), &
            shape(sparse%levels(relation + 1)%patches(child)%state))
          sparse%levels(relation + 1)%patches(child)%temperature = reshape( &
            payload(sparse%nvar * child_geometry%nx * &
              child_geometry%ny + 1:), &
            shape(sparse%levels(relation + 1)%patches(child)%temperature))
        end if
      end if
      call all_ranks_accept_2d( &
        distribution%comm, entity_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) return
      if (allocated(child_state)) deallocate(child_state)
      if (allocated(child_temperature)) deallocate(child_temperature)
      if (allocated(payload)) deallocate(payload)
    end do
    ok = .true.
  end subroutine prolong_direct_sparse_regrid_parent_2d

  subroutine transfer_direct_sparse_regrid_tree_overlaps_2d( &
      old_distribution, new_distribution, old_sparse, new_sparse, &
      transferred_cells, local_transfers, ok)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: &
      old_distribution, new_distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(in) :: old_sparse
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(inout) :: new_sparse
    integer, intent(out) :: transferred_cells, local_transfers
    logical, intent(out) :: ok

    type(eb_geometry_2d) :: new_geometry, old_geometry
    real(dp) :: new_i_offset, new_j_offset, old_i_offset, old_j_offset
    real(dp) :: overlap_x_lower, overlap_x_upper
    real(dp) :: overlap_y_lower, overlap_y_upper, scale, tolerance
    integer :: common_levels, i, j, level, new_i_lower, new_j_lower
    integer :: new_owner, new_patch, nx, ny, old_i_lower, old_j_lower
    integer :: old_owner, old_patch
    logical :: accepted, entity_ok, global_ok

    ok = .false.
    transferred_cells = 0
    local_transfers = 0
    common_levels = min(old_sparse%level_count(), new_sparse%level_count())
    do level = 2, common_levels
      do old_patch = 1, old_sparse%levels(level)%patch_count()
        call topology_patch_geometry_2d( &
          old_sparse%topology, level - 1, old_patch, old_geometry, entity_ok)
        if (.not. entity_ok) return
        old_owner = old_distribution%owner_of(level - 1, old_patch)
        do new_patch = 1, new_sparse%levels(level)%patch_count()
          call topology_patch_geometry_2d( &
            new_sparse%topology, level - 1, new_patch, new_geometry, &
            entity_ok)
          if (.not. entity_ok) return
          scale = max( &
            1.0_dp, abs(old_geometry%dx), abs(old_geometry%dy), &
            abs(new_geometry%dx), abs(new_geometry%dy))
          tolerance = sparse_regrid_geometry_tolerance * scale
          if (abs(old_geometry%dx - new_geometry%dx) > tolerance .or. &
              abs(old_geometry%dy - new_geometry%dy) > tolerance) cycle
          overlap_x_lower = max( &
            old_geometry%x_lower, new_geometry%x_lower)
          overlap_x_upper = min( &
            old_geometry%x_upper, new_geometry%x_upper)
          overlap_y_lower = max( &
            old_geometry%y_lower, new_geometry%y_lower)
          overlap_y_upper = min( &
            old_geometry%y_upper, new_geometry%y_upper)
          if (overlap_x_upper <= overlap_x_lower + tolerance .or. &
              overlap_y_upper <= overlap_y_lower + tolerance) cycle
          old_i_offset = &
            (overlap_x_lower - old_geometry%x_lower) / old_geometry%dx
          old_j_offset = &
            (overlap_y_lower - old_geometry%y_lower) / old_geometry%dy
          new_i_offset = &
            (overlap_x_lower - new_geometry%x_lower) / new_geometry%dx
          new_j_offset = &
            (overlap_y_lower - new_geometry%y_lower) / new_geometry%dy
          entity_ok = &
            abs(old_i_offset - real(nint(old_i_offset), dp)) <= &
              tolerance / old_geometry%dx .and. &
            abs(old_j_offset - real(nint(old_j_offset), dp)) <= &
              tolerance / old_geometry%dy .and. &
            abs(new_i_offset - real(nint(new_i_offset), dp)) <= &
              tolerance / new_geometry%dx .and. &
            abs(new_j_offset - real(nint(new_j_offset), dp)) <= &
              tolerance / new_geometry%dy
          if (.not. entity_ok) return
          old_i_lower = nint(old_i_offset) + 1
          old_j_lower = nint(old_j_offset) + 1
          new_i_lower = nint(new_i_offset) + 1
          new_j_lower = nint(new_j_offset) + 1
          nx = nint( &
            (overlap_x_upper - overlap_x_lower) / old_geometry%dx)
          ny = nint( &
            (overlap_y_upper - overlap_y_lower) / old_geometry%dy)
          if (nx < 1 .or. ny < 1) cycle
          entity_ok = &
            old_i_lower + nx - 1 <= old_geometry%nx .and. &
            old_j_lower + ny - 1 <= old_geometry%ny .and. &
            new_i_lower + nx - 1 <= new_geometry%nx .and. &
            new_j_lower + ny - 1 <= new_geometry%ny
          if (.not. entity_ok) return
          do j = 0, ny - 1
            do i = 0, nx - 1
              entity_ok = sparse_regrid_overlap_cell_geometry_matches_2d( &
                old_geometry, old_i_lower + i, old_j_lower + j, &
                new_geometry, new_i_lower + i, new_j_lower + j)
              if (.not. entity_ok) return
            end do
          end do
          new_owner = new_distribution%owner_of(level - 1, new_patch)
          call transfer_sparse_regrid_tree_rectangle_2d( &
            old_distribution, old_owner, new_owner, old_sparse%nvar, &
            old_sparse%levels(level)%patches(old_patch), old_i_lower, &
            old_j_lower, new_sparse%levels(level)%patches(new_patch), &
            new_i_lower, new_j_lower, nx, ny, entity_ok)
          call all_ranks_accept_2d( &
            old_distribution%comm, entity_ok, accepted, global_ok)
          if (.not. global_ok .or. .not. accepted) return
          if (nx > huge(transferred_cells) / ny) return
          if (transferred_cells > huge(transferred_cells) - nx * ny) return
          transferred_cells = transferred_cells + nx * ny
          if (old_owner /= new_owner .and. &
              old_distribution%rank == old_owner) &
            local_transfers = local_transfers + 1
        end do
      end do
    end do
    ok = .true.
  end subroutine transfer_direct_sparse_regrid_tree_overlaps_2d

  subroutine transfer_sparse_regrid_tree_rectangle_2d( &
      distribution, old_owner, new_owner, nvar, old_node, old_i_lower, &
      old_j_lower, new_node, new_i_lower, new_j_lower, nx, ny, ok)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    integer, intent(in) :: old_owner, new_owner, nvar
    type(mpi_sparse_reactive_amr_eb_patch_tree_node_2d), intent(in) :: old_node
    integer, intent(in) :: old_i_lower, old_j_lower
    type(mpi_sparse_reactive_amr_eb_patch_tree_node_2d), intent(inout) :: &
      new_node
    integer, intent(in) :: new_i_lower, new_j_lower, nx, ny
    logical, intent(out) :: ok

    type(MPI_Status) :: status
    real(dp), allocatable :: payload(:)
    integer :: ierr, state_count, value_count

    ok = old_owner >= 0 .and. old_owner < distribution%nranks .and. &
      new_owner >= 0 .and. new_owner < distribution%nranks .and. &
      nvar >= 1 .and. old_i_lower >= 1 .and. old_j_lower >= 1 .and. &
      new_i_lower >= 1 .and. new_j_lower >= 1 .and. nx >= 1 .and. ny >= 1
    if (.not. ok) return
    if (old_owner == new_owner) then
      if (distribution%rank == old_owner) then
        new_node%state(:, &
            new_i_lower:new_i_lower + nx - 1, &
            new_j_lower:new_j_lower + ny - 1) = &
          old_node%state(:, &
            old_i_lower:old_i_lower + nx - 1, &
            old_j_lower:old_j_lower + ny - 1)
        new_node%temperature( &
            new_i_lower:new_i_lower + nx - 1, &
            new_j_lower:new_j_lower + ny - 1) = &
          old_node%temperature( &
            old_i_lower:old_i_lower + nx - 1, &
            old_j_lower:old_j_lower + ny - 1)
      end if
      return
    end if
    if (distribution%rank /= old_owner .and. &
        distribution%rank /= new_owner) return

    state_count = nvar * nx * ny
    value_count = state_count + nx * ny
    allocate(payload(value_count))
    if (distribution%rank == old_owner) then
      payload(1:state_count) = reshape( &
        old_node%state(:, old_i_lower:old_i_lower + nx - 1, &
          old_j_lower:old_j_lower + ny - 1), [state_count])
      payload(state_count + 1:) = reshape( &
        old_node%temperature( &
          old_i_lower:old_i_lower + nx - 1, &
          old_j_lower:old_j_lower + ny - 1), [nx * ny])
      call MPI_Send( &
        payload, value_count, MPI_DOUBLE_PRECISION, new_owner, &
        sparse_tree_regrid_overlap_tag, distribution%comm, ierr)
      ok = ierr == MPI_SUCCESS
      return
    end if

    call MPI_Recv( &
      payload, value_count, MPI_DOUBLE_PRECISION, old_owner, &
      sparse_tree_regrid_overlap_tag, distribution%comm, status, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    new_node%state(:, &
        new_i_lower:new_i_lower + nx - 1, &
        new_j_lower:new_j_lower + ny - 1) = &
      reshape(payload(1:state_count), [nvar, nx, ny])
    new_node%temperature( &
        new_i_lower:new_i_lower + nx - 1, &
        new_j_lower:new_j_lower + ny - 1) = &
      reshape(payload(state_count + 1:), [nx, ny])
    ok = .true.
  end subroutine transfer_sparse_regrid_tree_rectangle_2d

  pure logical function sparse_regrid_overlap_cell_geometry_matches_2d( &
      first, first_i, first_j, second, second_i, second_j) result(matches)
    type(eb_geometry_2d), intent(in) :: first, second
    integer, intent(in) :: first_i, first_j, second_i, second_j

    matches = first%cell_type(first_i, first_j) == &
        second%cell_type(second_i, second_j) .and. &
      abs(first%volume_fraction(first_i, first_j) - &
        second%volume_fraction(second_i, second_j)) <= &
        sparse_regrid_geometry_tolerance .and. &
      abs(first%cell_centroid_x(first_i, first_j) - &
        second%cell_centroid_x(second_i, second_j)) <= &
        sparse_regrid_geometry_tolerance .and. &
      abs(first%cell_centroid_y(first_i, first_j) - &
        second%cell_centroid_y(second_i, second_j)) <= &
        sparse_regrid_geometry_tolerance .and. &
      abs(first%boundary_length(first_i, first_j) - &
        second%boundary_length(second_i, second_j)) <= &
        sparse_regrid_geometry_tolerance .and. &
      abs(first%boundary_centroid_x(first_i, first_j) - &
        second%boundary_centroid_x(second_i, second_j)) <= &
        sparse_regrid_geometry_tolerance .and. &
      abs(first%boundary_centroid_y(first_i, first_j) - &
        second%boundary_centroid_y(second_i, second_j)) <= &
        sparse_regrid_geometry_tolerance .and. &
      abs(first%boundary_normal_x(first_i, first_j) - &
        second%boundary_normal_x(second_i, second_j)) <= &
        sparse_regrid_geometry_tolerance .and. &
      abs(first%boundary_normal_y(first_i, first_j) - &
        second%boundary_normal_y(second_i, second_j)) <= &
        sparse_regrid_geometry_tolerance .and. &
      abs(first%boundary_normal_integral_x(first_i, first_j) - &
        second%boundary_normal_integral_x(second_i, second_j)) <= &
        sparse_regrid_geometry_tolerance .and. &
      abs(first%boundary_normal_integral_y(first_i, first_j) - &
        second%boundary_normal_integral_y(second_i, second_j)) <= &
        sparse_regrid_geometry_tolerance .and. &
      all(abs(first%x_face_fraction(first_i - 1:first_i, first_j) - &
        second%x_face_fraction(second_i - 1:second_i, second_j)) <= &
        sparse_regrid_geometry_tolerance) .and. &
      all(abs(first%x_face_centroid_y(first_i - 1:first_i, first_j) - &
        second%x_face_centroid_y(second_i - 1:second_i, second_j)) <= &
        sparse_regrid_geometry_tolerance) .and. &
      all(abs(first%y_face_fraction(first_i, first_j - 1:first_j) - &
        second%y_face_fraction(second_i, second_j - 1:second_j)) <= &
        sparse_regrid_geometry_tolerance) .and. &
      all(abs(first%y_face_centroid_x(first_i, first_j - 1:first_j) - &
        second%y_face_centroid_x(second_i, second_j - 1:second_j)) <= &
        sparse_regrid_geometry_tolerance)
  end function sparse_regrid_overlap_cell_geometry_matches_2d

  subroutine compute_sparse_owned_reactive_amr_eb_patch_tree_timestep_2d( &
      species, transport, distribution, sparse, hydro_cfl, transport_cfl, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, dt, ok, local_active_nodes)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(in) :: sparse
    real(dp), intent(in) :: hydro_cfl, transport_cfl
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_active_nodes

    type(eb_geometry_2d) :: geometry
    real(dp) :: global_dt, level_scale, local_dt, maximum_diffusivity
    real(dp) :: node_dt, numeric_maximum(2), numeric_minimum(2)
    real(dp) :: numeric_values(2), scaled_dt
    integer :: global_nodes, ierr, integer_maximum(4), integer_minimum(4)
    integer :: integer_values(4), level, local_nodes, patch
    integer :: refinement_ratio
    logical :: accepted, entity_ok, global_ok, local_ok, matches
    logical :: transport_active

    dt = 0.0_dp
    ok = .false.
    local_nodes = 0
    if (present(local_active_nodes)) local_active_nodes = 0
    transport_active = viscosity_enabled .or. &
      thermal_conduction_enabled .or. species_diffusion_enabled
    numeric_values = [hydro_cfl, transport_cfl]
    integer_values = [ &
      size(species), merge(1, 0, viscosity_enabled), &
      merge(1, 0, thermal_conduction_enabled), &
      merge(1, 0, species_diffusion_enabled)]

    call replicated_distribution_matches_2d( &
      distribution, sparse%topology, matches)
    local_ok = matches .and. sparse%is_valid(distribution) .and. &
      sparse%nvar == reactive_nvar(size(species)) .and. &
      all(ieee_is_finite(numeric_values)) .and. &
      hydro_cfl > 0.0_dp .and. hydro_cfl <= 1.0_dp .and. &
      transport_cfl > 0.0_dp .and. transport_cfl <= 0.5_dp
    if (local_ok .and. transport_active) &
      local_ok = compatible_transport_database(species, transport)
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    call MPI_Allreduce( &
      numeric_values, numeric_minimum, size(numeric_values), &
      MPI_DOUBLE_PRECISION, MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      numeric_values, numeric_maximum, size(numeric_values), &
      MPI_DOUBLE_PRECISION, MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_values, integer_minimum, size(integer_values), MPI_INTEGER, &
      MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_values, integer_maximum, size(integer_values), MPI_INTEGER, &
      MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. &
        any(numeric_minimum /= numeric_maximum) .or. &
        any(integer_minimum /= integer_maximum)) return

    local_dt = huge(1.0_dp)
    level_scale = 1.0_dp
    entity_ok = .true.
    do level = 1, sparse%level_count()
      if (level > 1) then
        refinement_ratio = &
          sparse%topology%relations(level - 1)%refinement_ratio
        if (level_scale > huge(1.0_dp) / real(refinement_ratio, dp)) then
          entity_ok = .false.
          exit
        end if
        level_scale = level_scale * real(refinement_ratio, dp)
      end if
      do patch = 1, sparse%levels(level)%patch_count()
        if (.not. distribution%is_local(level - 1, patch)) cycle
        call topology_patch_geometry_2d( &
          sparse%topology, level - 1, patch, geometry, entity_ok)
        if (.not. entity_ok) exit
        if (count(geometry%cell_type /= eb_covered_cell) == 0) cycle
        local_nodes = local_nodes + 1
        call compute_reactive_eb_cfl_timestep_2d( &
          species, sparse%levels(level)%patches(patch)%state, &
          sparse%levels(level)%patches(patch)%temperature, geometry, &
          hydro_cfl, node_dt, entity_ok)
        if (.not. entity_ok) exit
        if (node_dt > huge(1.0_dp) / level_scale) then
          scaled_dt = huge(1.0_dp)
        else
          scaled_dt = level_scale * node_dt
        end if
        local_dt = min(local_dt, scaled_dt)
        if (transport_active) then
          call reactive_eb_transport_timestep_2d( &
            species, transport, &
            sparse%levels(level)%patches(patch)%state, &
            sparse%levels(level)%patches(patch)%temperature, geometry, &
            transport_cfl, viscosity_enabled, &
            thermal_conduction_enabled, species_diffusion_enabled, &
            node_dt, maximum_diffusivity, entity_ok)
          if (.not. entity_ok) exit
          if (node_dt > huge(1.0_dp) / level_scale) then
            scaled_dt = huge(1.0_dp)
          else
            scaled_dt = level_scale * node_dt
          end if
          local_dt = min(local_dt, scaled_dt)
        end if
      end do
      if (.not. entity_ok) exit
    end do

    local_ok = entity_ok .and. ieee_is_finite(local_dt) .and. &
      local_dt > 0.0_dp
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Allreduce( &
      local_dt, global_dt, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      local_nodes, global_nodes, 1, MPI_INTEGER, MPI_SUM, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. global_nodes < 1 .or. &
        .not. ieee_is_finite(global_dt) .or. global_dt <= 0.0_dp .or. &
        global_dt >= huge(1.0_dp)) return

    dt = global_dt
    ok = .true.
    if (present(local_active_nodes)) local_active_nodes = local_nodes
  end subroutine compute_sparse_owned_reactive_amr_eb_patch_tree_timestep_2d

  subroutine advance_sparse_owned_reactive_amr_eb_patch_tree_chemistry_2d( &
      species, reactions, distribution, sparse, interval, rtol, atol, ok, &
      local_level_advances, local_restriction_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(inout) :: sparse
    real(dp), intent(in) :: interval, rtol, atol
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_level_advances(:)
    integer, intent(out), optional :: local_restriction_transfers

    type(mpi_sparse_reactive_amr_eb_patch_tree_2d) :: candidate
    type(eb_geometry_2d) :: geometry, parent_geometry
    type(MPI_Status) :: status
    real(dp), allocatable :: child_state(:, :, :)
    real(dp), allocatable :: state_work(:, :, :), temperature_work(:, :)
    real(dp) :: control_maximum(3), control_minimum(3), controls(3)
    logical, allocatable :: active_mask(:, :)
    integer, allocatable :: advances(:)
    integer :: child, child_owner, ierr, integer_maximum(2)
    integer :: integer_minimum(2), integer_values(2), level, parent
    integer :: parent_owner, patch, relation, transfers
    logical :: accepted, entity_ok, global_ok, local_ok, matches

    ok = .false.
    transfers = 0
    if (present(local_level_advances)) local_level_advances = 0
    if (present(local_restriction_transfers)) &
      local_restriction_transfers = 0
    controls = [interval, rtol, atol]
    integer_values = [size(species), size(reactions)]

    call replicated_distribution_matches_2d( &
      distribution, sparse%topology, matches)
    local_ok = matches .and. sparse%is_valid(distribution) .and. &
      sparse%nvar == reactive_nvar(size(species)) .and. &
      size(species) >= 1 .and. size(reactions) >= 1 .and. &
      all(ieee_is_finite(controls)) .and. interval >= 0.0_dp .and. &
      rtol > 0.0_dp .and. atol > 0.0_dp
    if (present(local_level_advances)) local_ok = local_ok .and. &
      size(local_level_advances) == sparse%level_count()
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    call MPI_Allreduce( &
      controls, control_minimum, size(controls), MPI_DOUBLE_PRECISION, &
      MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      controls, control_maximum, size(controls), MPI_DOUBLE_PRECISION, &
      MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_values, integer_minimum, size(integer_values), MPI_INTEGER, &
      MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_values, integer_maximum, size(integer_values), MPI_INTEGER, &
      MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. &
        any(control_minimum /= control_maximum) .or. &
        any(integer_minimum /= integer_maximum)) return

    candidate = sparse
    allocate(advances(candidate%level_count()), source=0)
    do level = 1, candidate%level_count()
      do patch = 1, candidate%levels(level)%patch_count()
        entity_ok = .true.
        if (distribution%is_local(level - 1, patch)) then
          call topology_patch_geometry_2d( &
            candidate%topology, level - 1, patch, geometry, entity_ok)
          if (entity_ok) then
            allocate(active_mask(geometry%nx, geometry%ny))
            active_mask = geometry%cell_type /= eb_covered_cell
            call advance_reactive_chemistry_2d( &
              species, reactions, &
              candidate%levels(level)%patches(patch)%state, &
              candidate%levels(level)%patches(patch)%temperature, &
              geometry%nx, geometry%ny, interval, rtol, atol, entity_ok, &
              active_mask)
            deallocate(active_mask)
          end if
          if (entity_ok) then
            allocate(temperature_work(geometry%nx, geometry%ny))
            call recover_transport_temperature_2d( &
              species, candidate%levels(level)%patches(patch)%state, &
              candidate%levels(level)%patches(patch)%temperature, geometry, &
              temperature_work, entity_ok)
            if (entity_ok) &
              candidate%levels(level)%patches(patch)%temperature = &
                temperature_work
            deallocate(temperature_work)
          end if
          if (entity_ok) advances(level) = advances(level) + 1
        end if
        call all_ranks_accept_2d( &
          distribution%comm, entity_ok, accepted, global_ok)
        if (.not. global_ok .or. .not. accepted) return
      end do
    end do

    do relation = size(candidate%topology%relations), 1, -1
      do child = 1, candidate%topology%relations(relation)% &
          child_patch_count()
        parent = candidate%topology%relations(relation)% &
          children(child)%parent_patch
        parent_owner = distribution%owner_of(relation - 1, parent)
        child_owner = distribution%owner_of(relation, child)
        entity_ok = parent_owner >= 0 .and. child_owner >= 0
        if (parent_owner /= child_owner) then
          if (distribution%rank == child_owner) then
            call MPI_Send( &
              candidate%levels(relation + 1)%patches(child)%state, &
              size(candidate%levels(relation + 1)%patches(child)%state), &
              MPI_DOUBLE_PRECISION, parent_owner, &
              sparse_tree_restriction_tag, distribution%comm, ierr)
            entity_ok = ierr == MPI_SUCCESS
            if (entity_ok) transfers = transfers + 1
          else if (distribution%rank == parent_owner) then
            geometry = candidate%topology%relations(relation)% &
              children(child)%geometry
            allocate(child_state(candidate%nvar, geometry%nx, geometry%ny))
            call MPI_Recv( &
              child_state, size(child_state), MPI_DOUBLE_PRECISION, &
              child_owner, sparse_tree_restriction_tag, distribution%comm, &
              status, ierr)
            entity_ok = ierr == MPI_SUCCESS
          end if
        end if
        if (distribution%rank == parent_owner .and. entity_ok) then
          call topology_patch_geometry_2d( &
            candidate%topology, relation - 1, parent, parent_geometry, &
            entity_ok)
          if (entity_ok) then
            allocate(state_work, mold= &
              candidate%levels(relation)%patches(parent)%state)
            allocate(temperature_work, mold= &
              candidate%levels(relation)%patches(parent)%temperature)
            if (parent_owner == child_owner) then
              call average_down_reactive_eb_state_patch_2d( &
                species, &
                candidate%levels(relation)%patches(parent)%state, &
                candidate%levels(relation)%patches(parent)%temperature, &
                parent_geometry, &
                candidate%levels(relation + 1)%patches(child)%state, &
                candidate%topology%relations(relation)% &
                  children(child)%geometry, &
                candidate%topology%relations(relation)%children(child)% &
                  patch, state_work, temperature_work, entity_ok)
            else
              call average_down_reactive_eb_state_patch_2d( &
                species, &
                candidate%levels(relation)%patches(parent)%state, &
                candidate%levels(relation)%patches(parent)%temperature, &
                parent_geometry, child_state, &
                candidate%topology%relations(relation)% &
                  children(child)%geometry, &
                candidate%topology%relations(relation)%children(child)% &
                  patch, state_work, temperature_work, entity_ok)
            end if
            if (entity_ok) then
              candidate%levels(relation)%patches(parent)%state = state_work
              candidate%levels(relation)%patches(parent)%temperature = &
                temperature_work
            end if
            deallocate(state_work, temperature_work)
          end if
        end if
        if (allocated(child_state)) deallocate(child_state)
        call all_ranks_accept_2d( &
          distribution%comm, entity_ok, accepted, global_ok)
        if (.not. global_ok .or. .not. accepted) return
      end do
    end do

    local_ok = candidate%is_valid(distribution)
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    sparse = candidate
    ok = .true.
    if (present(local_level_advances)) local_level_advances = advances
    if (present(local_restriction_transfers)) &
      local_restriction_transfers = transfers
  end subroutine &
    advance_sparse_owned_reactive_amr_eb_patch_tree_chemistry_2d

  subroutine advance_sparse_owned_reactive_amr_eb_patch_tree_hydro_2d( &
      species, distribution, sparse, solver, reconstruction, limiter, &
      state_redist_max_order, dt, ok, &
      state_redist_target_volume_fraction, failure_context, &
      local_level_advances, local_entity_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(inout) :: sparse
    character(len=*), intent(in) :: solver, reconstruction, limiter
    integer, intent(in) :: state_redist_max_order
    real(dp), intent(in) :: dt
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: state_redist_target_volume_fraction
    character(len=*), intent(out), optional :: failure_context
    integer, intent(out), optional :: local_level_advances(:)
    integer, intent(out), optional :: local_entity_transfers

    type(mpi_sparse_reactive_amr_eb_patch_tree_2d) :: candidate
    real(dp), allocatable :: x_flux(:, :, :), y_flux(:, :, :)
    real(dp) :: numeric_maximum(2), numeric_minimum(2)
    real(dp) :: numeric_values(2), selected_target
    integer, allocatable :: advances(:)
    integer :: character_index, ierr, integer_maximum(2)
    integer :: integer_minimum(2), integer_values(2), transfers
    integer :: string_codes(32, 3), string_maximum(32, 3)
    integer :: string_minimum(32, 3)
    character(len=160) :: context
    logical :: accepted, global_ok, local_ok, matches

    ok = .false.
    transfers = 0
    context = "input validation"
    if (present(failure_context)) failure_context = context
    if (present(local_level_advances)) local_level_advances = 0
    if (present(local_entity_transfers)) local_entity_transfers = 0
    selected_target = 0.5_dp
    if (present(state_redist_target_volume_fraction)) &
      selected_target = state_redist_target_volume_fraction
    numeric_values = [dt, selected_target]
    integer_values = [state_redist_max_order, size(species)]
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
    call replicated_distribution_matches_2d( &
      distribution, sparse%topology, matches)
    local_ok = local_ok .and. matches .and. sparse%is_valid(distribution) &
      .and. sparse%nvar == reactive_nvar(size(species)) .and. &
      size(species) >= 1 .and. all(ieee_is_finite(numeric_values)) .and. &
      dt > 0.0_dp .and. selected_target > 0.0_dp .and. &
      selected_target <= 1.0_dp .and. &
      (state_redist_max_order == 0 .or. state_redist_max_order == 2)
    if (present(local_level_advances)) local_ok = local_ok .and. &
      size(local_level_advances) == sparse%level_count()
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    call MPI_Allreduce( &
      numeric_values, numeric_minimum, size(numeric_values), &
      MPI_DOUBLE_PRECISION, MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      numeric_values, numeric_maximum, size(numeric_values), &
      MPI_DOUBLE_PRECISION, MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_values, integer_minimum, size(integer_values), MPI_INTEGER, &
      MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_values, integer_maximum, size(integer_values), MPI_INTEGER, &
      MPI_MAX, distribution%comm, ierr)
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

    candidate = sparse
    allocate(advances(candidate%level_count()), source=0)
    call advance_sparse_amr_eb_patch_tree_hydro_node_2d( &
      species, distribution, candidate, 1, 1, trim(solver), &
      trim(reconstruction), trim(limiter), state_redist_max_order, &
      selected_target, dt, x_flux, y_flux, advances, transfers, context, &
      local_ok)
    if (.not. local_ok) then
      if (present(failure_context)) failure_context = context
      return
    end if

    context = "final sparse validation"
    local_ok = candidate%is_valid(distribution)
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) then
      if (present(failure_context)) failure_context = context
      return
    end if
    sparse = candidate
    ok = .true.
    if (present(failure_context)) failure_context = "none"
    if (present(local_level_advances)) local_level_advances = advances
    if (present(local_entity_transfers)) local_entity_transfers = transfers
  end subroutine &
    advance_sparse_owned_reactive_amr_eb_patch_tree_hydro_2d

  subroutine advance_sparse_owned_reactive_amr_eb_patch_tree_transport_2d( &
      species, transport, distribution, sparse, interval, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, max_order, minimum_theta, ok, failure_context, &
      local_level_advances, local_entity_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(inout) :: sparse
    real(dp), intent(in) :: interval, target_volume_fraction
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    integer, intent(in) :: max_order
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    character(len=*), intent(out), optional :: failure_context
    integer, intent(out), optional :: local_level_advances(:)
    integer, intent(out), optional :: local_entity_transfers

    type(mpi_sparse_reactive_amr_eb_patch_tree_2d) :: candidate, start
    type(eb_geometry_2d) :: geometry
    real(dp), allocatable :: temperature_work(:, :)
    real(dp) :: first_theta, numeric_maximum(2), numeric_minimum(2)
    real(dp) :: numeric_values(2), second_theta
    integer, allocatable :: first_advances(:), second_advances(:)
    integer :: child, ierr, integer_maximum(6), integer_minimum(6)
    integer :: integer_values(6), level, parent, patch, relation, transfers
    character(len=160) :: context
    logical :: accepted, boundaries_match, global_ok, local_ok, matches
    logical :: transport_active

    ok = .false.
    minimum_theta = 1.0_dp
    transfers = 0
    context = "input validation"
    if (present(failure_context)) failure_context = context
    if (present(local_level_advances)) local_level_advances = 0
    if (present(local_entity_transfers)) local_entity_transfers = 0
    transport_active = viscosity_enabled .or. &
      thermal_conduction_enabled .or. species_diffusion_enabled
    numeric_values = [interval, target_volume_fraction]
    integer_values = [ &
      max_order, size(species), merge(1, 0, viscosity_enabled), &
      merge(1, 0, thermal_conduction_enabled), &
      merge(1, 0, species_diffusion_enabled), &
      merge(1, 0, barodiffusion_enabled)]

    call replicated_distribution_matches_2d( &
      distribution, sparse%topology, matches)
    call replicated_reactive_boundaries_match_2d( &
      boundaries, size(species), distribution%comm, boundaries_match)
    local_ok = matches .and. boundaries_match .and. &
      sparse%is_valid(distribution) .and. &
      sparse%nvar == reactive_nvar(size(species)) .and. &
      size(species) >= 1 .and. all(ieee_is_finite(numeric_values)) .and. &
      interval >= 0.0_dp .and. target_volume_fraction > 0.0_dp .and. &
      target_volume_fraction <= 1.0_dp .and. &
      (max_order == 0 .or. max_order == 2)
    if (transport_active) local_ok = local_ok .and. &
      compatible_transport_database(species, transport)
    if (present(local_level_advances)) local_ok = local_ok .and. &
      size(local_level_advances) == sparse%level_count()
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    call MPI_Allreduce( &
      numeric_values, numeric_minimum, size(numeric_values), &
      MPI_DOUBLE_PRECISION, MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      numeric_values, numeric_maximum, size(numeric_values), &
      MPI_DOUBLE_PRECISION, MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_values, integer_minimum, size(integer_values), MPI_INTEGER, &
      MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_values, integer_maximum, size(integer_values), MPI_INTEGER, &
      MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. &
        any(numeric_minimum /= numeric_maximum) .or. &
        any(integer_minimum /= integer_maximum)) return
    if (interval <= tiny(1.0_dp) .or. .not. transport_active) then
      ok = .true.
      if (present(failure_context)) failure_context = "none"
      return
    end if

    start = sparse
    candidate = sparse
    allocate(first_advances(candidate%level_count()), source=0)
    allocate(second_advances(candidate%level_count()), source=0)
    call advance_sparse_amr_eb_patch_tree_transport_euler_2d( &
      species, transport, distribution, candidate, interval, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, max_order, first_theta, first_advances, &
      transfers, context, local_ok)
    if (.not. local_ok) then
      if (present(failure_context)) failure_context = context
      return
    end if
    call advance_sparse_amr_eb_patch_tree_transport_euler_2d( &
      species, transport, distribution, candidate, interval, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, max_order, second_theta, second_advances, &
      transfers, context, local_ok)
    if (.not. local_ok) then
      if (present(failure_context)) failure_context = context
      return
    end if

    do level = 1, candidate%level_count()
      do patch = 1, candidate%levels(level)%patch_count()
        context = "stage blend geometry"
        call topology_patch_geometry_2d( &
          candidate%topology, level - 1, patch, geometry, local_ok)
        if (distribution%is_local(level - 1, patch) .and. local_ok) then
          candidate%levels(level)%patches(patch)%state = 0.5_dp * ( &
            start%levels(level)%patches(patch)%state + &
            candidate%levels(level)%patches(patch)%state)
          allocate(temperature_work(geometry%nx, geometry%ny))
          context = "stage blend temperature recovery"
          call recover_transport_temperature_2d( &
            species, candidate%levels(level)%patches(patch)%state, &
            0.5_dp * (start%levels(level)%patches(patch)%temperature + &
              candidate%levels(level)%patches(patch)%temperature), geometry, &
            temperature_work, local_ok)
          if (local_ok) candidate%levels(level)%patches(patch)%temperature = &
            temperature_work
          deallocate(temperature_work)
        end if
        call all_ranks_accept_2d( &
          distribution%comm, local_ok, accepted, global_ok)
        if (.not. global_ok .or. .not. accepted) then
          if (present(failure_context)) failure_context = context
          return
        end if
      end do
    end do

    context = "final hierarchy synchronization"
    do relation = size(candidate%topology%relations), 1, -1
      do child = 1, candidate%topology%relations(relation)% &
          child_patch_count()
        parent = candidate%topology%relations(relation)% &
          children(child)%parent_patch
        call average_down_sparse_patch_tree_edge_2d( &
          species, distribution, candidate, relation, parent, child, &
          transfers, local_ok)
        if (.not. local_ok) then
          if (present(failure_context)) failure_context = context
          return
        end if
      end do
    end do
    local_ok = candidate%is_valid(distribution)
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) then
      if (present(failure_context)) failure_context = context
      return
    end if

    sparse = candidate
    minimum_theta = min(first_theta, second_theta)
    ok = .true.
    if (present(failure_context)) failure_context = "none"
    if (present(local_level_advances)) &
      local_level_advances = first_advances + second_advances
    if (present(local_entity_transfers)) local_entity_transfers = transfers
  end subroutine &
    advance_sparse_owned_reactive_amr_eb_patch_tree_transport_2d

  subroutine advance_sparse_owned_reactive_amr_eb_patch_tree_full_physics_2d( &
      species, reactions, transport, distribution, sparse, solver, &
      reconstruction, limiter, max_order, dt, chemistry_enabled, rtol, atol, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, minimum_transport_theta, ok, failure_context, &
      local_chemistry_level_advances, local_transport_level_advances, &
      local_hydro_level_advances, local_chemistry_transfers, &
      local_transport_transfers, local_hydro_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(inout) :: sparse
    character(len=*), intent(in) :: solver, reconstruction, limiter
    integer, intent(in) :: max_order
    real(dp), intent(in) :: dt, rtol, atol
    logical, intent(in) :: chemistry_enabled
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    real(dp), intent(in) :: target_volume_fraction
    real(dp), intent(out) :: minimum_transport_theta
    logical, intent(out) :: ok
    character(len=*), intent(out), optional :: failure_context
    integer, intent(out), optional :: local_chemistry_level_advances(:)
    integer, intent(out), optional :: local_transport_level_advances(:)
    integer, intent(out), optional :: local_hydro_level_advances(:)
    integer, intent(out), optional :: local_chemistry_transfers
    integer, intent(out), optional :: local_transport_transfers
    integer, intent(out), optional :: local_hydro_transfers

    type(mpi_sparse_reactive_amr_eb_patch_tree_2d) :: candidate
    real(dp) :: first_transport_theta, second_transport_theta
    real(dp) :: numeric_maximum(4), numeric_minimum(4), numeric_values(4)
    integer, allocatable :: first_chemistry(:), second_chemistry(:)
    integer, allocatable :: first_transport(:), second_transport(:)
    integer, allocatable :: hydro_advances(:)
    integer :: first_chemistry_transfers, first_transport_transfers
    integer :: hydro_transfers, ierr, integer_maximum(9)
    integer :: integer_minimum(9), integer_values(9)
    integer :: second_chemistry_transfers, second_transport_transfers
    character(len=160) :: context
    logical :: accepted, global_ok, local_ok, matches, transport_active

    ok = .false.
    minimum_transport_theta = 1.0_dp
    context = "input validation"
    if (present(failure_context)) failure_context = context
    if (present(local_chemistry_level_advances)) &
      local_chemistry_level_advances = 0
    if (present(local_transport_level_advances)) &
      local_transport_level_advances = 0
    if (present(local_hydro_level_advances)) local_hydro_level_advances = 0
    if (present(local_chemistry_transfers)) local_chemistry_transfers = 0
    if (present(local_transport_transfers)) local_transport_transfers = 0
    if (present(local_hydro_transfers)) local_hydro_transfers = 0
    first_transport_theta = 1.0_dp
    second_transport_theta = 1.0_dp
    first_chemistry_transfers = 0
    second_chemistry_transfers = 0
    first_transport_transfers = 0
    second_transport_transfers = 0
    hydro_transfers = 0
    transport_active = viscosity_enabled .or. &
      thermal_conduction_enabled .or. species_diffusion_enabled
    numeric_values = [dt, rtol, atol, target_volume_fraction]
    integer_values = [ &
      size(species), size(reactions), size(transport), max_order, &
      merge(1, 0, chemistry_enabled), merge(1, 0, viscosity_enabled), &
      merge(1, 0, thermal_conduction_enabled), &
      merge(1, 0, species_diffusion_enabled), &
      merge(1, 0, barodiffusion_enabled)]

    call replicated_distribution_matches_2d( &
      distribution, sparse%topology, matches)
    local_ok = matches .and. sparse%is_valid(distribution) .and. &
      sparse%nvar == reactive_nvar(size(species)) .and. &
      size(species) >= 1 .and. all(ieee_is_finite(numeric_values)) .and. &
      dt > 0.0_dp .and. rtol > 0.0_dp .and. atol > 0.0_dp .and. &
      target_volume_fraction > 0.0_dp .and. &
      target_volume_fraction <= 1.0_dp .and. &
      (max_order == 0 .or. max_order == 2)
    if (chemistry_enabled) local_ok = local_ok .and. size(reactions) >= 1
    if (transport_active) local_ok = local_ok .and. &
      compatible_transport_database(species, transport)
    if (present(local_chemistry_level_advances)) local_ok = local_ok .and. &
      size(local_chemistry_level_advances) == sparse%level_count()
    if (present(local_transport_level_advances)) local_ok = local_ok .and. &
      size(local_transport_level_advances) == sparse%level_count()
    if (present(local_hydro_level_advances)) local_ok = local_ok .and. &
      size(local_hydro_level_advances) == sparse%level_count()
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    call MPI_Allreduce( &
      numeric_values, numeric_minimum, size(numeric_values), &
      MPI_DOUBLE_PRECISION, MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      numeric_values, numeric_maximum, size(numeric_values), &
      MPI_DOUBLE_PRECISION, MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_values, integer_minimum, size(integer_values), MPI_INTEGER, &
      MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_values, integer_maximum, size(integer_values), MPI_INTEGER, &
      MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. &
        any(numeric_minimum /= numeric_maximum) .or. &
        any(integer_minimum /= integer_maximum)) return

    candidate = sparse
    allocate(first_chemistry(candidate%level_count()), source=0)
    allocate(second_chemistry(candidate%level_count()), source=0)
    allocate(first_transport(candidate%level_count()), source=0)
    allocate(second_transport(candidate%level_count()), source=0)
    allocate(hydro_advances(candidate%level_count()), source=0)

    if (chemistry_enabled) then
      context = "first chemistry"
      call advance_sparse_owned_reactive_amr_eb_patch_tree_chemistry_2d( &
        species, reactions, distribution, candidate, 0.5_dp * dt, rtol, &
        atol, local_ok, first_chemistry, first_chemistry_transfers)
      if (.not. local_ok) then
        if (present(failure_context)) failure_context = context
        return
      end if
    end if

    context = "first transport"
    call advance_sparse_owned_reactive_amr_eb_patch_tree_transport_2d( &
      species, transport, distribution, candidate, 0.5_dp * dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, max_order, first_transport_theta, local_ok, &
      context, first_transport, first_transport_transfers)
    if (.not. local_ok) then
      if (present(failure_context)) failure_context = context
      return
    end if

    context = "hydro"
    call advance_sparse_owned_reactive_amr_eb_patch_tree_hydro_2d( &
      species, distribution, candidate, solver, reconstruction, limiter, &
      max_order, dt, local_ok, target_volume_fraction, context, &
      hydro_advances, hydro_transfers)
    if (.not. local_ok) then
      if (present(failure_context)) failure_context = context
      return
    end if

    context = "second transport"
    call advance_sparse_owned_reactive_amr_eb_patch_tree_transport_2d( &
      species, transport, distribution, candidate, 0.5_dp * dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, max_order, second_transport_theta, local_ok, &
      context, second_transport, second_transport_transfers)
    if (.not. local_ok) then
      if (present(failure_context)) failure_context = context
      return
    end if

    if (chemistry_enabled) then
      context = "second chemistry"
      call advance_sparse_owned_reactive_amr_eb_patch_tree_chemistry_2d( &
        species, reactions, distribution, candidate, 0.5_dp * dt, rtol, &
        atol, local_ok, second_chemistry, second_chemistry_transfers)
      if (.not. local_ok) then
        if (present(failure_context)) failure_context = context
        return
      end if
    end if

    context = "advance count overflow"
    local_ok = &
      all(first_chemistry <= huge(1) - second_chemistry) .and. &
      all(first_transport <= huge(1) - second_transport) .and. &
      first_chemistry_transfers <= huge(1) - second_chemistry_transfers .and. &
      first_transport_transfers <= huge(1) - second_transport_transfers
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) then
      if (present(failure_context)) failure_context = context
      return
    end if

    context = "final sparse validation"
    local_ok = candidate%is_valid(distribution)
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) then
      if (present(failure_context)) failure_context = context
      return
    end if

    sparse = candidate
    minimum_transport_theta = min( &
      first_transport_theta, second_transport_theta)
    ok = .true.
    if (present(failure_context)) failure_context = "none"
    if (present(local_chemistry_level_advances)) &
      local_chemistry_level_advances = first_chemistry + second_chemistry
    if (present(local_transport_level_advances)) &
      local_transport_level_advances = first_transport + second_transport
    if (present(local_hydro_level_advances)) &
      local_hydro_level_advances = hydro_advances
    if (present(local_chemistry_transfers)) local_chemistry_transfers = &
      first_chemistry_transfers + second_chemistry_transfers
    if (present(local_transport_transfers)) local_transport_transfers = &
      first_transport_transfers + second_transport_transfers
    if (present(local_hydro_transfers)) local_hydro_transfers = &
      hydro_transfers
  end subroutine &
    advance_sparse_owned_reactive_amr_eb_patch_tree_full_physics_2d

  subroutine advance_sparse_owned_reactive_amr_eb_patch_tree_to_time_2d( &
      species, reactions, transport, distribution, sparse, solver, &
      reconstruction, limiter, max_order, time, final_time, steps, &
      maximum_steps, hydro_cfl, transport_cfl, chemistry_enabled, rtol, atol, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, minimum_dt, minimum_transport_theta, ok, &
      failure_context, advanced_steps, local_timestep_evaluations, &
      local_chemistry_level_advances, local_transport_level_advances, &
      local_hydro_level_advances, local_chemistry_transfers, &
      local_transport_transfers, local_hydro_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(inout) :: sparse
    character(len=*), intent(in) :: solver, reconstruction, limiter
    integer, intent(in) :: max_order
    real(dp), intent(inout) :: time
    real(dp), intent(in) :: final_time
    integer, intent(inout) :: steps
    integer, intent(in) :: maximum_steps
    real(dp), intent(in) :: hydro_cfl, transport_cfl, rtol, atol
    logical, intent(in) :: chemistry_enabled
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    real(dp), intent(in) :: target_volume_fraction
    real(dp), intent(out) :: minimum_dt, minimum_transport_theta
    logical, intent(out) :: ok
    character(len=*), intent(out), optional :: failure_context
    integer, intent(out), optional :: advanced_steps
    integer, intent(out), optional :: local_timestep_evaluations
    integer, intent(out), optional :: local_chemistry_level_advances(:)
    integer, intent(out), optional :: local_transport_level_advances(:)
    integer, intent(out), optional :: local_hydro_level_advances(:)
    integer, intent(out), optional :: local_chemistry_transfers
    integer, intent(out), optional :: local_transport_transfers
    integer, intent(out), optional :: local_hydro_transfers

    type(mpi_sparse_reactive_amr_eb_patch_tree_2d) :: candidate
    real(dp) :: dt, numeric_maximum(7), numeric_minimum(7)
    real(dp) :: numeric_values(7), remaining, step_theta, time_tolerance
    integer, allocatable :: accumulated_chemistry(:)
    integer, allocatable :: accumulated_hydro(:)
    integer, allocatable :: accumulated_transport(:)
    integer, allocatable :: step_chemistry(:), step_hydro(:)
    integer, allocatable :: step_transport(:)
    integer :: accumulated_chemistry_transfers
    integer :: accumulated_hydro_transfers
    integer :: accumulated_timestep_evaluations
    integer :: accumulated_transport_transfers, completed_steps
    integer :: ierr, integer_maximum(11), integer_minimum(11)
    integer :: integer_values(11), step_chemistry_transfers
    integer :: step_hydro_transfers, step_timestep_evaluations
    integer :: step_transport_transfers
    character(len=160) :: context
    logical :: accepted, global_ok, local_ok, matches, transport_active

    ok = .false.
    minimum_dt = 0.0_dp
    minimum_transport_theta = 1.0_dp
    context = "input validation"
    if (present(failure_context)) failure_context = context
    if (present(advanced_steps)) advanced_steps = 0
    if (present(local_timestep_evaluations)) local_timestep_evaluations = 0
    if (present(local_chemistry_level_advances)) &
      local_chemistry_level_advances = 0
    if (present(local_transport_level_advances)) &
      local_transport_level_advances = 0
    if (present(local_hydro_level_advances)) local_hydro_level_advances = 0
    if (present(local_chemistry_transfers)) local_chemistry_transfers = 0
    if (present(local_transport_transfers)) local_transport_transfers = 0
    if (present(local_hydro_transfers)) local_hydro_transfers = 0
    accumulated_timestep_evaluations = 0
    accumulated_chemistry_transfers = 0
    accumulated_transport_transfers = 0
    accumulated_hydro_transfers = 0
    completed_steps = 0
    transport_active = viscosity_enabled .or. &
      thermal_conduction_enabled .or. species_diffusion_enabled
    numeric_values = [ &
      time, final_time, hydro_cfl, transport_cfl, rtol, atol, &
      target_volume_fraction]
    integer_values = [ &
      steps, maximum_steps, max_order, size(species), size(reactions), &
      size(transport), merge(1, 0, chemistry_enabled), &
      merge(1, 0, viscosity_enabled), &
      merge(1, 0, thermal_conduction_enabled), &
      merge(1, 0, species_diffusion_enabled), &
      merge(1, 0, barodiffusion_enabled)]
    time_tolerance = 16.0_dp * epsilon(1.0_dp) * &
      max(tiny(1.0_dp), abs(final_time))

    call replicated_distribution_matches_2d( &
      distribution, sparse%topology, matches)
    local_ok = matches .and. sparse%is_valid(distribution) .and. &
      sparse%nvar == reactive_nvar(size(species)) .and. &
      size(species) >= 1 .and. all(ieee_is_finite(numeric_values)) .and. &
      time >= 0.0_dp .and. final_time >= time - time_tolerance .and. &
      steps >= 0 .and. maximum_steps >= steps .and. &
      hydro_cfl > 0.0_dp .and. hydro_cfl <= 1.0_dp .and. &
      transport_cfl > 0.0_dp .and. transport_cfl <= 0.5_dp .and. &
      rtol > 0.0_dp .and. atol > 0.0_dp .and. &
      target_volume_fraction > 0.0_dp .and. &
      target_volume_fraction <= 1.0_dp .and. &
      (max_order == 0 .or. max_order == 2)
    if (chemistry_enabled) local_ok = local_ok .and. size(reactions) >= 1
    if (transport_active) local_ok = local_ok .and. &
      compatible_transport_database(species, transport)
    if (present(local_chemistry_level_advances)) local_ok = local_ok .and. &
      size(local_chemistry_level_advances) == sparse%level_count()
    if (present(local_transport_level_advances)) local_ok = local_ok .and. &
      size(local_transport_level_advances) == sparse%level_count()
    if (present(local_hydro_level_advances)) local_ok = local_ok .and. &
      size(local_hydro_level_advances) == sparse%level_count()
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    call MPI_Allreduce( &
      numeric_values, numeric_minimum, size(numeric_values), &
      MPI_DOUBLE_PRECISION, MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      numeric_values, numeric_maximum, size(numeric_values), &
      MPI_DOUBLE_PRECISION, MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_values, integer_minimum, size(integer_values), MPI_INTEGER, &
      MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_values, integer_maximum, size(integer_values), MPI_INTEGER, &
      MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. &
        any(numeric_minimum /= numeric_maximum) .or. &
        any(integer_minimum /= integer_maximum)) return

    allocate(accumulated_chemistry(sparse%level_count()), source=0)
    allocate(accumulated_transport(sparse%level_count()), source=0)
    allocate(accumulated_hydro(sparse%level_count()), source=0)
    allocate(step_chemistry(sparse%level_count()), source=0)
    allocate(step_transport(sparse%level_count()), source=0)
    allocate(step_hydro(sparse%level_count()), source=0)

    do
      remaining = final_time - time
      if (remaining <= time_tolerance) exit
      if (steps >= maximum_steps) then
        if (present(failure_context)) failure_context = "maximum steps"
        return
      end if

      context = "timestep"
      call compute_sparse_owned_reactive_amr_eb_patch_tree_timestep_2d( &
        species, transport, distribution, sparse, hydro_cfl, transport_cfl, &
        viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, dt, local_ok, &
        step_timestep_evaluations)
      if (.not. local_ok) then
        if (present(failure_context)) failure_context = context
        return
      end if
      dt = min(dt, remaining)
      if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp) then
        if (present(failure_context)) failure_context = context
        return
      end if

      candidate = sparse
      step_chemistry = 0
      step_transport = 0
      step_hydro = 0
      step_chemistry_transfers = 0
      step_transport_transfers = 0
      step_hydro_transfers = 0
      context = "full physics"
      call advance_sparse_owned_reactive_amr_eb_patch_tree_full_physics_2d( &
        species, reactions, transport, distribution, candidate, solver, &
        reconstruction, limiter, max_order, dt, chemistry_enabled, rtol, &
        atol, viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, boundaries, &
        target_volume_fraction, step_theta, local_ok, context, &
        step_chemistry, step_transport, step_hydro, &
        step_chemistry_transfers, step_transport_transfers, &
        step_hydro_transfers)
      if (.not. local_ok) then
        if (present(failure_context)) failure_context = context
        return
      end if

      context = "advance count overflow"
      local_ok = &
        all(step_chemistry <= huge(1) - accumulated_chemistry) .and. &
        all(step_transport <= huge(1) - accumulated_transport) .and. &
        all(step_hydro <= huge(1) - accumulated_hydro) .and. &
        step_timestep_evaluations <= &
          huge(1) - accumulated_timestep_evaluations .and. &
        step_chemistry_transfers <= &
          huge(1) - accumulated_chemistry_transfers .and. &
        step_transport_transfers <= &
          huge(1) - accumulated_transport_transfers .and. &
        step_hydro_transfers <= huge(1) - accumulated_hydro_transfers
      call all_ranks_accept_2d( &
        distribution%comm, local_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) then
        if (present(failure_context)) failure_context = context
        return
      end if

      sparse = candidate
      time = time + dt
      steps = steps + 1
      completed_steps = completed_steps + 1
      accumulated_timestep_evaluations = &
        accumulated_timestep_evaluations + step_timestep_evaluations
      accumulated_chemistry = accumulated_chemistry + step_chemistry
      accumulated_transport = accumulated_transport + step_transport
      accumulated_hydro = accumulated_hydro + step_hydro
      accumulated_chemistry_transfers = &
        accumulated_chemistry_transfers + step_chemistry_transfers
      accumulated_transport_transfers = &
        accumulated_transport_transfers + step_transport_transfers
      accumulated_hydro_transfers = &
        accumulated_hydro_transfers + step_hydro_transfers
      if (minimum_dt == 0.0_dp) then
        minimum_dt = dt
      else
        minimum_dt = min(minimum_dt, dt)
      end if
      minimum_transport_theta = min(minimum_transport_theta, step_theta)
      if (present(advanced_steps)) advanced_steps = completed_steps
      if (present(local_timestep_evaluations)) &
        local_timestep_evaluations = accumulated_timestep_evaluations
      if (present(local_chemistry_level_advances)) &
        local_chemistry_level_advances = accumulated_chemistry
      if (present(local_transport_level_advances)) &
        local_transport_level_advances = accumulated_transport
      if (present(local_hydro_level_advances)) &
        local_hydro_level_advances = accumulated_hydro
      if (present(local_chemistry_transfers)) local_chemistry_transfers = &
        accumulated_chemistry_transfers
      if (present(local_transport_transfers)) local_transport_transfers = &
        accumulated_transport_transfers
      if (present(local_hydro_transfers)) local_hydro_transfers = &
        accumulated_hydro_transfers
    end do

    time = final_time
    ok = .true.
    if (present(failure_context)) failure_context = "none"
  end subroutine &
    advance_sparse_owned_reactive_amr_eb_patch_tree_to_time_2d

  subroutine advance_sparse_amr_eb_patch_tree_transport_euler_2d( &
      species, transport, distribution, sparse, interval, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, max_order, minimum_theta, advances, transfers, &
      failure_context, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(inout) :: sparse
    real(dp), intent(in) :: interval, target_volume_fraction
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    integer, intent(in) :: max_order
    real(dp), intent(out) :: minimum_theta
    integer, intent(out) :: advances(:)
    integer, intent(inout) :: transfers
    character(len=*), intent(inout) :: failure_context
    logical, intent(out) :: ok

    real(dp), allocatable :: x_flux(:, :, :), y_flux(:, :, :)
    real(dp) :: local_theta
    integer :: ierr
    logical :: accepted, global_ok, local_ok

    ok = .false.
    minimum_theta = 1.0_dp
    advances = 0
    local_theta = 1.0_dp
    call advance_sparse_amr_eb_patch_tree_transport_node_2d( &
      species, transport, distribution, sparse, 1, 1, interval, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, max_order, x_flux, y_flux, local_theta, &
      advances, transfers, failure_context, local_ok)
    if (.not. local_ok) return
    call MPI_Allreduce( &
      local_theta, minimum_theta, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. .not. ieee_is_finite(minimum_theta) .or. &
        minimum_theta < 0.0_dp .or. minimum_theta > 1.0_dp) return
    local_ok = sparse%is_valid(distribution)
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    ok = global_ok .and. accepted
  end subroutine advance_sparse_amr_eb_patch_tree_transport_euler_2d

  recursive subroutine advance_sparse_amr_eb_patch_tree_transport_node_2d( &
      species, transport, distribution, sparse, level, patch, interval, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, max_order, x_flux, y_flux, minimum_theta, &
      advances, transfers, failure_context, ok, exterior)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(inout) :: sparse
    integer, intent(in) :: level, patch, max_order
    real(dp), intent(in) :: interval, target_volume_fraction
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    real(dp), allocatable, intent(out) :: x_flux(:, :, :), y_flux(:, :, :)
    real(dp), intent(inout) :: minimum_theta
    integer, intent(inout) :: advances(:), transfers
    character(len=*), intent(inout) :: failure_context
    logical, intent(out) :: ok
    type(reactive_eb_exterior_state_2d), intent(in), optional :: exterior

    type(eb_geometry_2d) :: child_geometry, geometry
    type(amr_eb_flux_register_2d), allocatable :: registers(:)
    type(reactive_eb_patch_exterior_context_2d), allocatable :: contexts(:)
    type(reactive_eb_exterior_state_2d) :: child_exterior
    real(dp), allocatable :: child_x_flux(:, :, :), child_y_flux(:, :, :)
    real(dp), allocatable :: integral_before(:), rhs(:, :, :)
    real(dp), allocatable :: parent_fine_x_flux(:, :, :)
    real(dp), allocatable :: parent_fine_y_flux(:, :, :)
    real(dp), allocatable :: state_end(:, :, :), state_start(:, :, :)
    real(dp), allocatable :: temperature_end(:, :), temperature_start(:, :)
    real(dp) :: alpha, child_interval, node_theta
    integer :: child, child_count, child_owner, first_child, global_child
    integer :: owner, ratio, substep
    logical :: accepted, entity_ok, global_ok, local_ok

    ok = .false.
    entity_ok = level >= 1 .and. level <= sparse%level_count() .and. &
      patch >= 1 .and. patch <= sparse%levels(level)%patch_count() .and. &
      size(advances) == sparse%level_count() .and. &
      ieee_is_finite(interval) .and. interval > 0.0_dp
    call all_ranks_accept_2d( &
      distribution%comm, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call topology_patch_geometry_2d( &
      sparse%topology, level - 1, patch, geometry, entity_ok)
    owner = distribution%owner_of(level - 1, patch)
    entity_ok = entity_ok .and. owner >= 0 .and. owner < distribution%nranks
    call all_ranks_accept_2d( &
      distribution%comm, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    child_count = 0
    first_child = 1
    if (level < sparse%level_count()) then
      first_child = sparse%topology%relations(level)% &
        child_offsets(patch) + 1
      child_count = sparse%topology%relations(level)% &
          child_offsets(patch + 1) - &
        sparse%topology%relations(level)%child_offsets(patch)
    end if
    if (child_count > 0) then
      allocate(integral_before(sparse%nvar))
      call reduce_sparse_amr_eb_patch_subtree_integral_2d( &
        distribution, sparse, level, patch, integral_before, local_ok)
      if (.not. local_ok) return
    end if

    entity_ok = .true.
    write(failure_context, '(a,i0,a,i0)') &
      "transport level ", level - 1, " patch ", patch
    if (distribution%rank == owner) then
      allocate(state_start, source=sparse%levels(level)%patches(patch)%state)
      allocate(temperature_start, &
        source=sparse%levels(level)%patches(patch)%temperature)
      allocate(state_end, mold=state_start)
      allocate(temperature_end, mold=temperature_start)
      allocate(rhs, mold=state_start)
      allocate(x_flux(sparse%nvar, 0:geometry%nx, geometry%ny))
      allocate(y_flux(sparse%nvar, geometry%nx, 0:geometry%ny))
      if (present(exterior)) then
        call reactive_eb_transport_fluxes_rhs_2d( &
          species, transport, state_start, temperature_start, geometry, &
          interval, viscosity_enabled, thermal_conduction_enabled, &
          species_diffusion_enabled, barodiffusion_enabled, boundaries, rhs, &
          x_flux, y_flux, node_theta, entity_ok, exterior)
      else
        call reactive_eb_transport_fluxes_rhs_2d( &
          species, transport, state_start, temperature_start, geometry, &
          interval, viscosity_enabled, thermal_conduction_enabled, &
          species_diffusion_enabled, barodiffusion_enabled, boundaries, rhs, &
          x_flux, y_flux, node_theta, entity_ok)
      end if
      if (entity_ok) then
        minimum_theta = min(minimum_theta, node_theta)
        call advance_reactive_eb_state_redistributed_2d( &
          species, state_start, temperature_start, geometry, rhs, interval, &
          state_end, temperature_end, entity_ok, target_volume_fraction, &
          max_order)
      end if
      if (entity_ok) then
        sparse%levels(level)%patches(patch)%state = state_end
        sparse%levels(level)%patches(patch)%temperature = temperature_end
        advances(level) = advances(level) + 1
      end if
    end if
    call all_ranks_accept_2d( &
      distribution%comm, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    if (child_count == 0) then
      ok = .true.
      return
    end if

    ratio = sparse%topology%relations(level)%refinement_ratio
    child_interval = interval / real(ratio, dp)
    allocate(registers(child_count), contexts(child_count))
    do child = 1, child_count
      global_child = first_child + child - 1
      child_geometry = sparse%topology%relations(level)% &
        children(global_child)%geometry
      child_owner = distribution%owner_of(level, global_child)
      write(failure_context, '(a,i0,a,i0)') &
        "transport register level ", level, " child ", global_child
      entity_ok = child_owner >= 0 .and. child_owner < distribution%nranks
      if (distribution%rank == owner .and. entity_ok) then
        call initialize_amr_eb_flux_register_2d( &
          geometry, child_geometry, sparse%topology%relations(level)% &
            children(global_child)%patch, sparse%nvar, registers(child), &
          entity_ok)
        if (entity_ok) call accumulate_coarse_eb_fluxes_2d( &
          registers(child), geometry, child_geometry, &
          sparse%topology%relations(level)%children(global_child)%patch, &
          x_flux, y_flux, interval, entity_ok)
      end if
      call all_ranks_accept_2d( &
        distribution%comm, entity_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) return
      write(failure_context, '(a,i0,a,i0)') &
        "transport context level ", level, " child ", global_child
      call transfer_sparse_patch_tree_parent_context_2d( &
        distribution, sparse%nvar, owner, child_owner, state_start, &
        temperature_start, state_end, temperature_end, geometry, &
        child_geometry, sparse%topology%relations(level)% &
          children(global_child)%patch, contexts(child), transfers, local_ok)
      if (.not. local_ok) return
    end do

    do substep = 1, ratio
      alpha = real(substep - 1, dp) / real(ratio, dp)
      do child = 1, child_count
        global_child = first_child + child - 1
        child_geometry = sparse%topology%relations(level)% &
          children(global_child)%geometry
        child_owner = distribution%owner_of(level, global_child)
        write(failure_context, '(a,i0,a,i0,a,i0)') &
          "transport exterior level ", level, " child ", global_child, &
          " substep ", substep
        child_exterior = reactive_eb_exterior_state_2d()
        entity_ok = .true.
        if (distribution%rank == child_owner) then
          call build_reactive_eb_patch_exterior_from_context_2d( &
            species, contexts(child), geometry, child_geometry, &
            sparse%topology%relations(level)%children(global_child)%patch, &
            alpha, child_exterior, entity_ok, &
            sparse%levels(level + 1)%patches(global_child)%state, &
            sparse%levels(level + 1)%patches(global_child)%temperature)
        end if
        call all_ranks_accept_2d( &
          distribution%comm, entity_ok, accepted, global_ok)
        if (.not. global_ok .or. .not. accepted) return
        call advance_sparse_amr_eb_patch_tree_transport_node_2d( &
          species, transport, distribution, sparse, level + 1, global_child, &
          child_interval, viscosity_enabled, thermal_conduction_enabled, &
          species_diffusion_enabled, barodiffusion_enabled, boundaries, &
          target_volume_fraction, max_order, child_x_flux, child_y_flux, &
          minimum_theta, advances, transfers, failure_context, local_ok, &
          child_exterior)
        if (.not. local_ok) return
        write(failure_context, '(a,i0,a,i0,a,i0)') &
          "transport fine flux level ", level, " child ", global_child, &
          " substep ", substep
        call transfer_sparse_patch_tree_flux_2d( &
          distribution, sparse%nvar, child_owner, owner, child_geometry, &
          child_x_flux, child_y_flux, parent_fine_x_flux, &
          parent_fine_y_flux, transfers, local_ok)
        if (.not. local_ok) return
        entity_ok = .true.
        if (distribution%rank == owner) call accumulate_fine_eb_fluxes_2d( &
          registers(child), geometry, child_geometry, &
          sparse%topology%relations(level)%children(global_child)%patch, &
          parent_fine_x_flux, parent_fine_y_flux, child_interval, entity_ok)
        call all_ranks_accept_2d( &
          distribution%comm, entity_ok, accepted, global_ok)
        if (.not. global_ok .or. .not. accepted) return
      end do
    end do

    do child = 1, child_count
      global_child = first_child + child - 1
      write(failure_context, '(a,i0,a,i0)') &
        "transport reflux level ", level, " child ", global_child
      call reflux_sparse_patch_tree_edge_2d( &
        species, distribution, sparse, level, patch, global_child, &
        registers(child), transfers, local_ok)
      if (.not. local_ok) return
    end do
    do child = 1, child_count
      global_child = first_child + child - 1
      write(failure_context, '(a,i0,a,i0)') &
        "transport average level ", level, " child ", global_child
      call average_down_sparse_patch_tree_edge_2d( &
        species, distribution, sparse, level, patch, global_child, &
        transfers, local_ok)
      if (.not. local_ok) return
    end do

    write(failure_context, '(a,i0,a,i0)') &
      "transport closure level ", level - 1, " patch ", patch
    call close_sparse_amr_eb_patch_subtree_conservation_2d( &
      species, distribution, sparse, level, patch, integral_before, x_flux, &
      y_flux, interval, local_ok)
    if (.not. local_ok) return
    ok = .true.
  end subroutine advance_sparse_amr_eb_patch_tree_transport_node_2d

  recursive subroutine advance_sparse_amr_eb_patch_tree_hydro_node_2d( &
      species, distribution, sparse, level, patch, solver, reconstruction, &
      limiter, state_redist_max_order, selected_target, interval, x_flux, &
      y_flux, advances, transfers, failure_context, ok, exterior)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(inout) :: sparse
    integer, intent(in) :: level, patch, state_redist_max_order
    character(len=*), intent(in) :: solver, reconstruction, limiter
    real(dp), intent(in) :: selected_target, interval
    real(dp), allocatable, intent(out) :: x_flux(:, :, :), y_flux(:, :, :)
    integer, intent(inout) :: advances(:), transfers
    character(len=*), intent(inout) :: failure_context
    logical, intent(out) :: ok
    type(reactive_eb_exterior_state_2d), intent(in), optional :: exterior

    type(eb_geometry_2d) :: child_geometry, geometry
    type(amr_eb_flux_register_2d), allocatable :: registers(:)
    type(reactive_eb_patch_exterior_context_2d), allocatable :: contexts(:)
    type(reactive_eb_exterior_state_2d) :: child_exterior
    real(dp), allocatable :: child_x_flux(:, :, :), child_y_flux(:, :, :)
    real(dp), allocatable :: integral_before(:)
    real(dp), allocatable :: parent_fine_x_flux(:, :, :)
    real(dp), allocatable :: parent_fine_y_flux(:, :, :)
    real(dp), allocatable :: state_end(:, :, :), state_start(:, :, :)
    real(dp), allocatable :: temperature_end(:, :), temperature_start(:, :)
    real(dp) :: alpha, child_interval
    integer :: child, child_count, child_owner, first_child, global_child
    integer :: owner, ratio, substep
    logical :: accepted, entity_ok, global_ok, local_ok, requires_closure

    ok = .false.
    entity_ok = level >= 1 .and. level <= sparse%level_count() .and. &
      patch >= 1 .and. patch <= sparse%levels(level)%patch_count() .and. &
      size(advances) == sparse%level_count() .and. &
      ieee_is_finite(interval) .and. interval > 0.0_dp
    call all_ranks_accept_2d( &
      distribution%comm, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call topology_patch_geometry_2d( &
      sparse%topology, level - 1, patch, geometry, entity_ok)
    owner = distribution%owner_of(level - 1, patch)
    entity_ok = entity_ok .and. owner >= 0 .and. owner < distribution%nranks
    call all_ranks_accept_2d( &
      distribution%comm, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    child_count = 0
    first_child = 1
    if (level < sparse%level_count()) then
      first_child = sparse%topology%relations(level)% &
        child_offsets(patch) + 1
      child_count = sparse%topology%relations(level)% &
          child_offsets(patch + 1) - &
        sparse%topology%relations(level)%child_offsets(patch)
    end if
    requires_closure = child_count > 0
    if (requires_closure) then
      allocate(integral_before(sparse%nvar))
      call reduce_sparse_amr_eb_patch_subtree_integral_2d( &
        distribution, sparse, level, patch, integral_before, local_ok)
      if (.not. local_ok) return
    end if

    entity_ok = .true.
    write(failure_context, '(a,i0,a,i0)') &
      "level advance level ", level - 1, " patch ", patch
    if (distribution%rank == owner) then
      allocate(state_start, source=sparse%levels(level)%patches(patch)%state)
      allocate(temperature_start, &
        source=sparse%levels(level)%patches(patch)%temperature)
      allocate(state_end, mold=state_start)
      allocate(temperature_end, mold=temperature_start)
      allocate(x_flux(sparse%nvar, 0:geometry%nx, geometry%ny))
      allocate(y_flux(sparse%nvar, geometry%nx, 0:geometry%ny))
      if (present(exterior)) then
        call advance_reactive_eb_level_2d( &
          species, state_start, temperature_start, geometry, solver, &
          reconstruction, limiter, selected_target, state_redist_max_order, &
          interval, state_end, temperature_end, x_flux, y_flux, entity_ok, &
          exterior)
      else
        call advance_reactive_eb_level_2d( &
          species, state_start, temperature_start, geometry, solver, &
          reconstruction, limiter, selected_target, state_redist_max_order, &
          interval, state_end, temperature_end, x_flux, y_flux, entity_ok)
      end if
      if (entity_ok) then
        sparse%levels(level)%patches(patch)%state = state_end
        sparse%levels(level)%patches(patch)%temperature = temperature_end
        advances(level) = advances(level) + 1
      end if
    end if
    call all_ranks_accept_2d( &
      distribution%comm, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    if (child_count == 0) then
      ok = .true.
      return
    end if

    ratio = sparse%topology%relations(level)%refinement_ratio
    child_interval = interval / real(ratio, dp)
    allocate(registers(child_count), contexts(child_count))
    do child = 1, child_count
      global_child = first_child + child - 1
      child_geometry = sparse%topology%relations(level)% &
        children(global_child)%geometry
      child_owner = distribution%owner_of(level, global_child)
      write(failure_context, '(a,i0,a,i0)') &
        "flux register level ", level, " child ", global_child
      entity_ok = child_owner >= 0 .and. child_owner < distribution%nranks
      if (distribution%rank == owner .and. entity_ok) then
        call initialize_amr_eb_flux_register_2d( &
          geometry, child_geometry, sparse%topology%relations(level)% &
            children(global_child)%patch, sparse%nvar, registers(child), &
          entity_ok)
        if (entity_ok) call accumulate_coarse_eb_fluxes_2d( &
          registers(child), geometry, child_geometry, &
          sparse%topology%relations(level)%children(global_child)%patch, &
          x_flux, y_flux, interval, entity_ok)
      end if
      call all_ranks_accept_2d( &
        distribution%comm, entity_ok, accepted, global_ok)
      if (.not. global_ok .or. .not. accepted) return
      write(failure_context, '(a,i0,a,i0)') &
        "exterior context level ", level, " child ", global_child
      call transfer_sparse_patch_tree_parent_context_2d( &
        distribution, sparse%nvar, owner, child_owner, state_start, &
        temperature_start, state_end, temperature_end, geometry, &
        child_geometry, sparse%topology%relations(level)% &
          children(global_child)%patch, contexts(child), transfers, local_ok)
      if (.not. local_ok) return
    end do

    do substep = 1, ratio
      alpha = sparse_patch_tree_substep_alpha( &
        reconstruction, substep, ratio)
      do child = 1, child_count
        global_child = first_child + child - 1
        child_geometry = sparse%topology%relations(level)% &
          children(global_child)%geometry
        child_owner = distribution%owner_of(level, global_child)
        write(failure_context, '(a,i0,a,i0,a,i0)') &
          "exterior level ", level, " child ", global_child, &
          " substep ", substep
        child_exterior = reactive_eb_exterior_state_2d()
        entity_ok = .true.
        if (distribution%rank == child_owner) then
          call build_reactive_eb_patch_exterior_from_context_2d( &
            species, contexts(child), geometry, child_geometry, &
            sparse%topology%relations(level)%children(global_child)%patch, &
            alpha, child_exterior, entity_ok, &
            sparse%levels(level + 1)%patches(global_child)%state, &
            sparse%levels(level + 1)%patches(global_child)%temperature)
        end if
        call all_ranks_accept_2d( &
          distribution%comm, entity_ok, accepted, global_ok)
        if (.not. global_ok .or. .not. accepted) return
        call advance_sparse_amr_eb_patch_tree_hydro_node_2d( &
          species, distribution, sparse, level + 1, global_child, solver, &
          reconstruction, limiter, state_redist_max_order, selected_target, &
          child_interval, child_x_flux, child_y_flux, advances, transfers, &
          failure_context, local_ok, child_exterior)
        if (.not. local_ok) return
        write(failure_context, '(a,i0,a,i0,a,i0)') &
          "fine flux level ", level, " child ", global_child, &
          " substep ", substep
        call transfer_sparse_patch_tree_flux_2d( &
          distribution, sparse%nvar, child_owner, owner, child_geometry, &
          child_x_flux, child_y_flux, parent_fine_x_flux, &
          parent_fine_y_flux, transfers, local_ok)
        if (.not. local_ok) return
        entity_ok = .true.
        if (distribution%rank == owner) call accumulate_fine_eb_fluxes_2d( &
          registers(child), geometry, child_geometry, &
          sparse%topology%relations(level)%children(global_child)%patch, &
          parent_fine_x_flux, parent_fine_y_flux, child_interval, entity_ok)
        call all_ranks_accept_2d( &
          distribution%comm, entity_ok, accepted, global_ok)
        if (.not. global_ok .or. .not. accepted) return
      end do
    end do

    do child = 1, child_count
      global_child = first_child + child - 1
      write(failure_context, '(a,i0,a,i0)') &
        "reflux level ", level, " child ", global_child
      call reflux_sparse_patch_tree_edge_2d( &
        species, distribution, sparse, level, patch, global_child, &
        registers(child), transfers, local_ok)
      if (.not. local_ok) return
    end do
    do child = 1, child_count
      global_child = first_child + child - 1
      write(failure_context, '(a,i0,a,i0)') &
        "average down level ", level, " child ", global_child
      call average_down_sparse_patch_tree_edge_2d( &
        species, distribution, sparse, level, patch, global_child, &
        transfers, local_ok)
      if (.not. local_ok) return
    end do

    write(failure_context, '(a,i0,a,i0)') &
      "cut-interface closure level ", level - 1, " patch ", patch
    call close_sparse_amr_eb_patch_subtree_conservation_2d( &
      species, distribution, sparse, level, patch, integral_before, x_flux, &
      y_flux, interval, local_ok)
    if (.not. local_ok) return
    ok = .true.
  end subroutine advance_sparse_amr_eb_patch_tree_hydro_node_2d

  subroutine transfer_sparse_patch_tree_parent_context_2d( &
      distribution, component_count, source, destination, state_start, &
      temperature_start, state_end, temperature_end, parent_geometry, &
      child_geometry, patch, context, transfers, ok)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    integer, intent(in) :: component_count, source, destination
    real(dp), allocatable, intent(in) :: state_start(:, :, :)
    real(dp), allocatable, intent(in) :: temperature_start(:, :)
    real(dp), allocatable, intent(in) :: state_end(:, :, :)
    real(dp), allocatable, intent(in) :: temperature_end(:, :)
    type(eb_geometry_2d), intent(in) :: parent_geometry, child_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    type(reactive_eb_patch_exterior_context_2d), intent(out) :: context
    integer, intent(inout) :: transfers
    logical, intent(out) :: ok

    type(reactive_eb_patch_exterior_context_2d) :: source_context
    type(MPI_Status) :: status
    real(dp), allocatable :: payload(:)
    integer :: ierr, offset, value_count
    logical :: accepted, entity_ok, global_ok

    ok = .false.
    entity_ok = component_count >= 1 .and. source >= 0 .and. &
      source < distribution%nranks .and. destination >= 0 .and. &
      destination < distribution%nranks .and. &
      patch%is_valid(parent_geometry, child_geometry)
    if (distribution%rank == source .and. entity_ok) then
      entity_ok = allocated(state_start) .and. &
        allocated(temperature_start) .and. allocated(state_end) .and. &
        allocated(temperature_end)
      if (entity_ok) entity_ok = &
        all(shape(state_start) == &
          [component_count, parent_geometry%nx, parent_geometry%ny]) .and. &
        all(shape(temperature_start) == &
          [parent_geometry%nx, parent_geometry%ny]) .and. &
        all(shape(state_end) == shape(state_start)) .and. &
        all(shape(temperature_end) == shape(temperature_start)) .and. &
        all(ieee_is_finite(state_start)) .and. &
        all(ieee_is_finite(temperature_start)) .and. &
        all(ieee_is_finite(state_end)) .and. &
        all(ieee_is_finite(temperature_end))
      if (entity_ok) call extract_reactive_eb_patch_exterior_context_2d( &
        state_start, temperature_start, state_end, temperature_end, &
        parent_geometry, child_geometry, patch, component_count, &
        source_context, entity_ok)
    end if
    call all_ranks_accept_2d( &
      distribution%comm, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    if (source == destination) then
      if (distribution%rank == destination) context = source_context
    else
      value_count = 4 * (component_count + 1) * &
        (child_geometry%nx + child_geometry%ny)
      if (distribution%rank == source) then
        allocate(payload(value_count))
        offset = 0
        call pack_sparse_patch_tree_exterior_2d( &
          source_context%start, payload, offset)
        call pack_sparse_patch_tree_exterior_2d( &
          source_context%end, payload, offset)
        call MPI_Send( &
          payload, value_count, MPI_DOUBLE_PRECISION, destination, &
          sparse_tree_hydro_context_tag, distribution%comm, ierr)
        entity_ok = ierr == MPI_SUCCESS .and. offset == value_count
        if (entity_ok) transfers = transfers + 1
      else if (distribution%rank == destination) then
        allocate(payload(value_count))
        call MPI_Recv( &
          payload, value_count, MPI_DOUBLE_PRECISION, source, &
          sparse_tree_hydro_context_tag, distribution%comm, status, ierr)
        entity_ok = ierr == MPI_SUCCESS
        if (entity_ok) then
          offset = 0
          call unpack_sparse_patch_tree_exterior_2d( &
            payload, offset, component_count, child_geometry, &
            context%start)
          call unpack_sparse_patch_tree_exterior_2d( &
            payload, offset, component_count, child_geometry, context%end)
          entity_ok = offset == value_count
        end if
      end if
    end if
    if (distribution%rank == destination .and. entity_ok) &
      entity_ok = context%is_valid(child_geometry, component_count)
    call all_ranks_accept_2d( &
      distribution%comm, entity_ok, accepted, global_ok)
    ok = global_ok .and. accepted
  end subroutine transfer_sparse_patch_tree_parent_context_2d

  subroutine pack_sparse_patch_tree_exterior_2d(exterior, payload, offset)
    type(reactive_eb_exterior_state_2d), intent(in) :: exterior
    real(dp), intent(inout) :: payload(:)
    integer, intent(inout) :: offset

    integer :: count

    count = size(exterior%x_lower_state)
    payload(offset + 1:offset + count) = reshape( &
      exterior%x_lower_state, [count])
    offset = offset + count
    count = size(exterior%x_upper_state)
    payload(offset + 1:offset + count) = reshape( &
      exterior%x_upper_state, [count])
    offset = offset + count
    count = size(exterior%y_lower_state)
    payload(offset + 1:offset + count) = reshape( &
      exterior%y_lower_state, [count])
    offset = offset + count
    count = size(exterior%y_upper_state)
    payload(offset + 1:offset + count) = reshape( &
      exterior%y_upper_state, [count])
    offset = offset + count
    count = size(exterior%x_lower_temperature)
    payload(offset + 1:offset + count) = exterior%x_lower_temperature
    offset = offset + count
    count = size(exterior%x_upper_temperature)
    payload(offset + 1:offset + count) = exterior%x_upper_temperature
    offset = offset + count
    count = size(exterior%y_lower_temperature)
    payload(offset + 1:offset + count) = exterior%y_lower_temperature
    offset = offset + count
    count = size(exterior%y_upper_temperature)
    payload(offset + 1:offset + count) = exterior%y_upper_temperature
    offset = offset + count
  end subroutine pack_sparse_patch_tree_exterior_2d

  subroutine unpack_sparse_patch_tree_exterior_2d( &
      payload, offset, component_count, geometry, exterior)
    real(dp), intent(in) :: payload(:)
    integer, intent(inout) :: offset
    integer, intent(in) :: component_count
    type(eb_geometry_2d), intent(in) :: geometry
    type(reactive_eb_exterior_state_2d), intent(out) :: exterior

    integer :: count

    allocate(exterior%x_lower_state(component_count, geometry%ny))
    count = size(exterior%x_lower_state)
    exterior%x_lower_state = reshape( &
      payload(offset + 1:offset + count), shape(exterior%x_lower_state))
    offset = offset + count
    allocate(exterior%x_upper_state(component_count, geometry%ny))
    count = size(exterior%x_upper_state)
    exterior%x_upper_state = reshape( &
      payload(offset + 1:offset + count), shape(exterior%x_upper_state))
    offset = offset + count
    allocate(exterior%y_lower_state(component_count, geometry%nx))
    count = size(exterior%y_lower_state)
    exterior%y_lower_state = reshape( &
      payload(offset + 1:offset + count), shape(exterior%y_lower_state))
    offset = offset + count
    allocate(exterior%y_upper_state(component_count, geometry%nx))
    count = size(exterior%y_upper_state)
    exterior%y_upper_state = reshape( &
      payload(offset + 1:offset + count), shape(exterior%y_upper_state))
    offset = offset + count
    allocate(exterior%x_lower_temperature(geometry%ny))
    count = size(exterior%x_lower_temperature)
    exterior%x_lower_temperature = payload(offset + 1:offset + count)
    offset = offset + count
    allocate(exterior%x_upper_temperature(geometry%ny))
    count = size(exterior%x_upper_temperature)
    exterior%x_upper_temperature = payload(offset + 1:offset + count)
    offset = offset + count
    allocate(exterior%y_lower_temperature(geometry%nx))
    count = size(exterior%y_lower_temperature)
    exterior%y_lower_temperature = payload(offset + 1:offset + count)
    offset = offset + count
    allocate(exterior%y_upper_temperature(geometry%nx))
    count = size(exterior%y_upper_temperature)
    exterior%y_upper_temperature = payload(offset + 1:offset + count)
    offset = offset + count
  end subroutine unpack_sparse_patch_tree_exterior_2d

  subroutine transfer_sparse_patch_tree_flux_2d( &
      distribution, component_count, source, destination, geometry, &
      source_x_flux, source_y_flux, destination_x_flux, destination_y_flux, &
      transfers, ok)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    integer, intent(in) :: component_count, source, destination
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), allocatable, intent(in) :: source_x_flux(:, :, :)
    real(dp), allocatable, intent(in) :: source_y_flux(:, :, :)
    real(dp), allocatable, intent(out) :: destination_x_flux(:, :, :)
    real(dp), allocatable, intent(out) :: destination_y_flux(:, :, :)
    integer, intent(inout) :: transfers
    logical, intent(out) :: ok

    type(MPI_Status) :: status
    real(dp), allocatable :: payload(:)
    integer :: ierr, value_count, x_count, y_count
    logical :: accepted, entity_ok, global_ok

    ok = .false.
    x_count = component_count * (geometry%nx + 1) * geometry%ny
    y_count = component_count * geometry%nx * (geometry%ny + 1)
    value_count = x_count + y_count
    entity_ok = component_count >= 1 .and. source >= 0 .and. &
      source < distribution%nranks .and. destination >= 0 .and. &
      destination < distribution%nranks .and. geometry%is_valid()
    if (distribution%rank == source .and. entity_ok) then
      entity_ok = allocated(source_x_flux) .and. allocated(source_y_flux)
      if (entity_ok) entity_ok = &
        all(shape(source_x_flux) == &
          [component_count, geometry%nx + 1, geometry%ny]) .and. &
        all(shape(source_y_flux) == &
          [component_count, geometry%nx, geometry%ny + 1]) .and. &
        all(ieee_is_finite(source_x_flux)) .and. &
        all(ieee_is_finite(source_y_flux))
    end if
    call all_ranks_accept_2d( &
      distribution%comm, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    if (source == destination) then
      if (distribution%rank == destination) then
        allocate(destination_x_flux, source=source_x_flux)
        allocate(destination_y_flux, source=source_y_flux)
      end if
    else if (distribution%rank == source) then
      allocate(payload(value_count))
      payload(1:x_count) = reshape(source_x_flux, [x_count])
      payload(x_count + 1:value_count) = reshape(source_y_flux, [y_count])
      call MPI_Send( &
        payload, value_count, MPI_DOUBLE_PRECISION, destination, &
        sparse_tree_hydro_flux_tag, distribution%comm, ierr)
      entity_ok = ierr == MPI_SUCCESS
      if (entity_ok) transfers = transfers + 1
    else if (distribution%rank == destination) then
      allocate(payload(value_count))
      call MPI_Recv( &
        payload, value_count, MPI_DOUBLE_PRECISION, source, &
        sparse_tree_hydro_flux_tag, distribution%comm, status, ierr)
      entity_ok = ierr == MPI_SUCCESS
      if (entity_ok) then
        allocate(destination_x_flux( &
          component_count, 0:geometry%nx, geometry%ny))
        allocate(destination_y_flux( &
          component_count, geometry%nx, 0:geometry%ny))
        destination_x_flux = reshape( &
          payload(1:x_count), shape(destination_x_flux))
        destination_y_flux = reshape( &
          payload(x_count + 1:value_count), shape(destination_y_flux))
      end if
    end if
    if (distribution%rank == destination .and. entity_ok) &
      entity_ok = allocated(destination_x_flux) .and. &
        allocated(destination_y_flux) .and. &
        all(ieee_is_finite(destination_x_flux)) .and. &
        all(ieee_is_finite(destination_y_flux))
    call all_ranks_accept_2d( &
      distribution%comm, entity_ok, accepted, global_ok)
    ok = global_ok .and. accepted
  end subroutine transfer_sparse_patch_tree_flux_2d

  subroutine pack_sparse_patch_tree_node_fields_2d( &
      state, temperature, payload)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    real(dp), intent(out) :: payload(:)

    integer :: state_count

    state_count = size(state)
    payload(1:state_count) = reshape(state, [state_count])
    payload(state_count + 1:) = reshape(temperature, [size(temperature)])
  end subroutine pack_sparse_patch_tree_node_fields_2d

  subroutine unpack_sparse_patch_tree_node_fields_2d( &
      payload, component_count, geometry, state, temperature)
    real(dp), intent(in) :: payload(:)
    integer, intent(in) :: component_count
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), allocatable, intent(out) :: state(:, :, :)
    real(dp), allocatable, intent(out) :: temperature(:, :)

    integer :: state_count

    state_count = component_count * geometry%nx * geometry%ny
    allocate(state(component_count, geometry%nx, geometry%ny))
    allocate(temperature(geometry%nx, geometry%ny))
    state = reshape(payload(1:state_count), shape(state))
    temperature = reshape(payload(state_count + 1:), shape(temperature))
  end subroutine unpack_sparse_patch_tree_node_fields_2d

  subroutine reflux_sparse_patch_tree_edge_2d( &
      species, distribution, sparse, level, parent, child, flux_register, &
      transfers, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(inout) :: sparse
    integer, intent(in) :: level, parent, child
    type(amr_eb_flux_register_2d), intent(inout) :: flux_register
    integer, intent(inout) :: transfers
    logical, intent(out) :: ok

    type(eb_geometry_2d) :: child_geometry, parent_geometry
    type(MPI_Status) :: status
    real(dp), allocatable :: child_state(:, :, :), child_temperature(:, :)
    real(dp), allocatable :: child_work(:, :, :)
    real(dp), allocatable :: child_work_temperature(:, :)
    real(dp), allocatable :: parent_work(:, :, :)
    real(dp), allocatable :: parent_work_temperature(:, :)
    real(dp), allocatable :: payload(:)
    integer :: child_owner, ierr, parent_owner, value_count
    logical :: accepted, entity_ok, global_ok

    ok = .false.
    call topology_patch_geometry_2d( &
      sparse%topology, level - 1, parent, parent_geometry, entity_ok)
    if (entity_ok) call topology_patch_geometry_2d( &
      sparse%topology, level, child, child_geometry, entity_ok)
    parent_owner = distribution%owner_of(level - 1, parent)
    child_owner = distribution%owner_of(level, child)
    entity_ok = entity_ok .and. parent_owner >= 0 .and. &
      parent_owner < distribution%nranks .and. child_owner >= 0 .and. &
      child_owner < distribution%nranks
    if (distribution%rank == parent_owner .and. entity_ok) &
      entity_ok = flux_register%is_valid( &
        parent_geometry, child_geometry, sparse%topology%relations(level)% &
          children(child)%patch)
    if (distribution%rank == child_owner .and. entity_ok) &
      entity_ok = sparse%levels(level + 1)%patches(child)%has_data()
    call all_ranks_accept_2d( &
      distribution%comm, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    value_count = (sparse%nvar + 1) * child_geometry%nx * child_geometry%ny
    if (parent_owner == child_owner) then
      if (distribution%rank == parent_owner) then
        allocate(child_state, &
          source=sparse%levels(level + 1)%patches(child)%state)
        allocate(child_temperature, &
          source=sparse%levels(level + 1)%patches(child)%temperature)
      end if
    else if (distribution%rank == child_owner) then
      allocate(payload(value_count))
      call pack_sparse_patch_tree_node_fields_2d( &
        sparse%levels(level + 1)%patches(child)%state, &
        sparse%levels(level + 1)%patches(child)%temperature, payload)
      call MPI_Send( &
        payload, value_count, MPI_DOUBLE_PRECISION, parent_owner, &
        sparse_tree_hydro_reflux_tag, distribution%comm, ierr)
      entity_ok = ierr == MPI_SUCCESS
      if (entity_ok) transfers = transfers + 1
    else if (distribution%rank == parent_owner) then
      allocate(payload(value_count))
      call MPI_Recv( &
        payload, value_count, MPI_DOUBLE_PRECISION, child_owner, &
        sparse_tree_hydro_reflux_tag, distribution%comm, status, ierr)
      entity_ok = ierr == MPI_SUCCESS
      if (entity_ok) call unpack_sparse_patch_tree_node_fields_2d( &
        payload, sparse%nvar, child_geometry, child_state, child_temperature)
    end if
    call all_ranks_accept_2d( &
      distribution%comm, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    if (distribution%rank == parent_owner) then
      allocate(parent_work, mold= &
        sparse%levels(level)%patches(parent)%state)
      allocate(parent_work_temperature, mold= &
        sparse%levels(level)%patches(parent)%temperature)
      allocate(child_work, mold=child_state)
      allocate(child_work_temperature, mold=child_temperature)
      call reflux_reactive_eb_state_patch_2d( &
        species, sparse%levels(level)%patches(parent)%state, &
        sparse%levels(level)%patches(parent)%temperature, parent_geometry, &
        child_state, child_temperature, child_geometry, &
        sparse%topology%relations(level)%children(child)%patch, &
        flux_register, parent_work, parent_work_temperature, child_work, &
        child_work_temperature, entity_ok)
      if (entity_ok) then
        sparse%levels(level)%patches(parent)%state = parent_work
        sparse%levels(level)%patches(parent)%temperature = &
          parent_work_temperature
        if (parent_owner == child_owner) then
          sparse%levels(level + 1)%patches(child)%state = child_work
          sparse%levels(level + 1)%patches(child)%temperature = &
            child_work_temperature
        end if
      end if
    end if
    call all_ranks_accept_2d( &
      distribution%comm, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    if (parent_owner /= child_owner) then
      if (distribution%rank == parent_owner) then
        if (allocated(payload)) deallocate(payload)
        allocate(payload(value_count))
        call pack_sparse_patch_tree_node_fields_2d( &
          child_work, child_work_temperature, payload)
        call MPI_Send( &
          payload, value_count, MPI_DOUBLE_PRECISION, child_owner, &
          sparse_tree_hydro_correction_tag, distribution%comm, ierr)
        entity_ok = ierr == MPI_SUCCESS
        if (entity_ok) transfers = transfers + 1
      else if (distribution%rank == child_owner) then
        if (allocated(payload)) deallocate(payload)
        allocate(payload(value_count))
        call MPI_Recv( &
          payload, value_count, MPI_DOUBLE_PRECISION, parent_owner, &
          sparse_tree_hydro_correction_tag, distribution%comm, status, ierr)
        entity_ok = ierr == MPI_SUCCESS
        if (entity_ok) then
          sparse%levels(level + 1)%patches(child)%state = reshape( &
            payload(1:sparse%nvar * child_geometry%nx * &
              child_geometry%ny), &
            shape(sparse%levels(level + 1)%patches(child)%state))
          sparse%levels(level + 1)%patches(child)%temperature = reshape( &
            payload(sparse%nvar * child_geometry%nx * &
              child_geometry%ny + 1:), &
            shape(sparse%levels(level + 1)%patches(child)%temperature))
        end if
      end if
    end if
    call all_ranks_accept_2d( &
      distribution%comm, entity_ok, accepted, global_ok)
    ok = global_ok .and. accepted
  end subroutine reflux_sparse_patch_tree_edge_2d

  subroutine average_down_sparse_patch_tree_edge_2d( &
      species, distribution, sparse, level, parent, child, transfers, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(inout) :: sparse
    integer, intent(in) :: level, parent, child
    integer, intent(inout) :: transfers
    logical, intent(out) :: ok

    type(eb_geometry_2d) :: child_geometry, parent_geometry
    type(MPI_Status) :: status
    real(dp), allocatable :: child_state(:, :, :), child_temperature(:, :)
    real(dp), allocatable :: parent_work(:, :, :)
    real(dp), allocatable :: parent_work_temperature(:, :)
    real(dp), allocatable :: payload(:)
    integer :: child_owner, ierr, parent_owner, value_count
    logical :: accepted, entity_ok, global_ok

    ok = .false.
    call topology_patch_geometry_2d( &
      sparse%topology, level - 1, parent, parent_geometry, entity_ok)
    if (entity_ok) call topology_patch_geometry_2d( &
      sparse%topology, level, child, child_geometry, entity_ok)
    parent_owner = distribution%owner_of(level - 1, parent)
    child_owner = distribution%owner_of(level, child)
    entity_ok = entity_ok .and. parent_owner >= 0 .and. &
      parent_owner < distribution%nranks .and. child_owner >= 0 .and. &
      child_owner < distribution%nranks
    if (distribution%rank == child_owner .and. entity_ok) &
      entity_ok = sparse%levels(level + 1)%patches(child)%has_data()
    call all_ranks_accept_2d( &
      distribution%comm, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    value_count = (sparse%nvar + 1) * child_geometry%nx * child_geometry%ny
    if (parent_owner == child_owner) then
      if (distribution%rank == parent_owner) then
        allocate(child_state, &
          source=sparse%levels(level + 1)%patches(child)%state)
        allocate(child_temperature, &
          source=sparse%levels(level + 1)%patches(child)%temperature)
      end if
    else if (distribution%rank == child_owner) then
      allocate(payload(value_count))
      call pack_sparse_patch_tree_node_fields_2d( &
        sparse%levels(level + 1)%patches(child)%state, &
        sparse%levels(level + 1)%patches(child)%temperature, payload)
      call MPI_Send( &
        payload, value_count, MPI_DOUBLE_PRECISION, parent_owner, &
        sparse_tree_hydro_average_tag, distribution%comm, ierr)
      entity_ok = ierr == MPI_SUCCESS
      if (entity_ok) transfers = transfers + 1
    else if (distribution%rank == parent_owner) then
      allocate(payload(value_count))
      call MPI_Recv( &
        payload, value_count, MPI_DOUBLE_PRECISION, child_owner, &
        sparse_tree_hydro_average_tag, distribution%comm, status, ierr)
      entity_ok = ierr == MPI_SUCCESS
      if (entity_ok) call unpack_sparse_patch_tree_node_fields_2d( &
        payload, sparse%nvar, child_geometry, child_state, child_temperature)
    end if
    call all_ranks_accept_2d( &
      distribution%comm, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    if (distribution%rank == parent_owner) then
      allocate(parent_work, mold= &
        sparse%levels(level)%patches(parent)%state)
      allocate(parent_work_temperature, mold= &
        sparse%levels(level)%patches(parent)%temperature)
      call average_down_reactive_eb_state_patch_2d( &
        species, sparse%levels(level)%patches(parent)%state, &
        sparse%levels(level)%patches(parent)%temperature, parent_geometry, &
        child_state, child_geometry, &
        sparse%topology%relations(level)%children(child)%patch, parent_work, &
        parent_work_temperature, entity_ok)
      if (entity_ok) then
        sparse%levels(level)%patches(parent)%state = parent_work
        sparse%levels(level)%patches(parent)%temperature = &
          parent_work_temperature
      end if
    end if
    call all_ranks_accept_2d( &
      distribution%comm, entity_ok, accepted, global_ok)
    ok = global_ok .and. accepted
  end subroutine average_down_sparse_patch_tree_edge_2d

  subroutine close_sparse_amr_eb_patch_subtree_conservation_2d( &
      species, distribution, sparse, level, patch, integral_before, x_flux, &
      y_flux, interval, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(inout) :: sparse
    integer, intent(in) :: level, patch
    real(dp), intent(in) :: integral_before(:)
    real(dp), allocatable, intent(in) :: x_flux(:, :, :), y_flux(:, :, :)
    real(dp), intent(in) :: interval
    logical, intent(out) :: ok

    type(eb_geometry_2d) :: geometry
    real(dp), allocatable :: boundary_change(:), correction(:)
    real(dp), allocatable :: current_integral(:), residual(:)
    real(dp), allocatable :: temperature_work(:, :)
    real(dp) :: closure_tolerance, recipient_volume, scale
    real(dp) :: species_residual
    integer :: component, i, ierr, j, k, owner
    logical :: accepted, entity_ok, global_ok, local_ok
    logical :: needs_correction

    ok = .false.
    call topology_patch_geometry_2d( &
      sparse%topology, level - 1, patch, geometry, entity_ok)
    owner = distribution%owner_of(level - 1, patch)
    entity_ok = entity_ok .and. owner >= 0 .and. &
      owner < distribution%nranks .and. &
      size(integral_before) == sparse%nvar .and. &
      ieee_is_finite(interval) .and. interval > 0.0_dp .and. &
      all(ieee_is_finite(integral_before))
    if (distribution%rank == owner .and. entity_ok) then
      entity_ok = allocated(x_flux) .and. allocated(y_flux)
      if (entity_ok) entity_ok = &
        all(shape(x_flux) == &
          [sparse%nvar, geometry%nx + 1, geometry%ny]) .and. &
        all(shape(y_flux) == &
          [sparse%nvar, geometry%nx, geometry%ny + 1]) .and. &
        all(ieee_is_finite(x_flux)) .and. all(ieee_is_finite(y_flux))
    end if
    call all_ranks_accept_2d( &
      distribution%comm, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    allocate(current_integral(sparse%nvar), boundary_change(sparse%nvar))
    allocate(residual(sparse%nvar), correction(sparse%nvar))
    call reduce_sparse_amr_eb_patch_subtree_integral_2d( &
      distribution, sparse, level, patch, current_integral, local_ok)
    if (.not. local_ok) return
    boundary_change = 0.0_dp
    correction = 0.0_dp
    residual = 0.0_dp
    needs_correction = .false.
    entity_ok = .true.
    if (distribution%rank == owner) then
      do j = 1, geometry%ny
        boundary_change = boundary_change + interval * geometry%dy * &
          (geometry%x_face_fraction(0, j) * x_flux(:, 0, j) - &
           geometry%x_face_fraction(geometry%nx, j) * &
             x_flux(:, geometry%nx, j))
      end do
      do i = 1, geometry%nx
        boundary_change = boundary_change + interval * geometry%dx * &
          (geometry%y_face_fraction(i, 0) * y_flux(:, i, 0) - &
           geometry%y_face_fraction(i, geometry%ny) * &
             y_flux(:, i, geometry%ny))
      end do
      residual = integral_before + boundary_change - current_integral
      correction(irho) = residual(irho)
      correction(iet) = residual(iet)
      species_residual = 0.0_dp
      do k = 1, size(species)
        component = reactive_species_component(k)
        correction(component) = residual(component)
        species_residual = species_residual + residual(component)
      end do
      scale = max(1.0_dp, maxval(abs(integral_before)), &
        maxval(abs(current_integral)), maxval(abs(boundary_change)))
      closure_tolerance = 4096.0_dp * epsilon(1.0_dp) * scale
      entity_ok = &
        abs(residual(irho) - species_residual) <= closure_tolerance
      if (entity_ok) then
        component = reactive_species_component(size(species))
        correction(component) = correction(component) + &
          residual(irho) - species_residual
        needs_correction = maxval(abs(correction)) > closure_tolerance
      end if
    end if
    call all_ranks_accept_2d( &
      distribution%comm, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Bcast( &
      needs_correction, 1, MPI_LOGICAL, owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    if (.not. needs_correction) then
      ok = .true.
      return
    end if

    entity_ok = .true.
    if (distribution%rank == owner) then
      recipient_volume = 0.0_dp
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (sparse_patch_tree_parent_cell_is_refined( &
                sparse%topology, level, patch, i, j) .or. &
              geometry%cell_type(i, j) == eb_covered_cell) cycle
          recipient_volume = recipient_volume + &
            geometry%volume_fraction(i, j) * geometry%dx * geometry%dy
        end do
      end do
      entity_ok = ieee_is_finite(recipient_volume) .and. &
        recipient_volume > tiny(1.0_dp)
      if (entity_ok) then
        correction = correction / recipient_volume
        entity_ok = all(ieee_is_finite(correction))
      end if
      if (entity_ok) then
        do j = 1, geometry%ny
          do i = 1, geometry%nx
            if (sparse_patch_tree_parent_cell_is_refined( &
                  sparse%topology, level, patch, i, j) .or. &
                geometry%cell_type(i, j) == eb_covered_cell) cycle
            sparse%levels(level)%patches(patch)%state(:, i, j) = &
              sparse%levels(level)%patches(patch)%state(:, i, j) + correction
          end do
        end do
        allocate(temperature_work(geometry%nx, geometry%ny))
        call recover_transport_temperature_2d( &
          species, sparse%levels(level)%patches(patch)%state, &
          sparse%levels(level)%patches(patch)%temperature, geometry, &
          temperature_work, entity_ok)
        if (entity_ok) sparse%levels(level)%patches(patch)%temperature = &
          temperature_work
      end if
    end if
    call all_ranks_accept_2d( &
      distribution%comm, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    call reduce_sparse_amr_eb_patch_subtree_integral_2d( &
      distribution, sparse, level, patch, current_integral, local_ok)
    if (.not. local_ok) return
    entity_ok = .true.
    if (distribution%rank == owner) then
      residual = integral_before + boundary_change - current_integral
      entity_ok = abs(residual(irho)) <= 8.0_dp * closure_tolerance .and. &
        abs(residual(iet)) <= 8.0_dp * closure_tolerance
      do k = 1, size(species)
        entity_ok = entity_ok .and. &
          abs(residual(reactive_species_component(k))) <= &
            8.0_dp * closure_tolerance
      end do
    end if
    call all_ranks_accept_2d( &
      distribution%comm, entity_ok, accepted, global_ok)
    ok = global_ok .and. accepted
  end subroutine close_sparse_amr_eb_patch_subtree_conservation_2d

  pure logical function sparse_patch_tree_parent_cell_is_refined( &
      topology, level, patch, i, j) result(refined)
    type(amr_eb_patch_tree_topology_2d), intent(in) :: topology
    integer, intent(in) :: level, patch, i, j

    type(amr_eb_patch_2d) :: child_patch
    integer :: child, first_child, last_child

    refined = .false.
    if (level >= topology%level_count()) return
    first_child = topology%relations(level)%child_offsets(patch) + 1
    last_child = topology%relations(level)%child_offsets(patch + 1)
    do child = first_child, last_child
      child_patch = topology%relations(level)%children(child)%patch
      refined = i >= child_patch%coarse_i_lower .and. &
        i <= child_patch%coarse_i_upper .and. &
        j >= child_patch%coarse_j_lower .and. &
        j <= child_patch%coarse_j_upper
      if (refined) return
    end do
  end function sparse_patch_tree_parent_cell_is_refined

  pure real(dp) function sparse_patch_tree_substep_alpha( &
      reconstruction, substep, ratio) result(alpha)
    character(len=*), intent(in) :: reconstruction
    integer, intent(in) :: substep, ratio

    if (trim(reconstruction) == "characteristic_plm") then
      alpha = (real(substep, dp) - 0.5_dp) / real(ratio, dp)
    else
      alpha = real(substep - 1, dp) / real(ratio, dp)
    end if
  end function sparse_patch_tree_substep_alpha

  subroutine composite_sparse_amr_eb_patch_tree_integral_2d( &
      distribution, sparse, integral, ok, local_nodes)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(in) :: sparse
    real(dp), intent(out) :: integral(:)
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_nodes

    call composite_sparse_amr_eb_patch_subtree_integral_2d( &
      distribution, sparse, 1, 1, integral, ok, local_nodes)
  end subroutine &
    composite_sparse_amr_eb_patch_tree_integral_2d

  subroutine &
      composite_sparse_amr_eb_patch_subtree_integral_2d( &
      distribution, sparse, level, patch, integral, ok, local_nodes)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(in) :: sparse
    integer, intent(in) :: level, patch
    real(dp), intent(out) :: integral(:)
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_nodes

    integer :: ierr, integer_maximum(3)
    integer :: integer_minimum(3), integer_values(3)
    logical :: accepted, global_ok, local_ok, matches

    integral = 0.0_dp
    ok = .false.
    if (present(local_nodes)) local_nodes = 0
    call replicated_distribution_matches_2d( &
      distribution, sparse%topology, matches)
    local_ok = matches .and. sparse%is_valid(distribution) .and. &
      size(integral) == sparse%nvar
    if (local_ok) local_ok = &
      level >= 1 .and. level <= sparse%level_count()
    if (local_ok) local_ok = &
      patch >= 1 .and. patch <= sparse%levels(level)%patch_count()
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    integer_values = [level, patch, size(integral)]
    call MPI_Allreduce( &
      integer_values, integer_minimum, size(integer_values), MPI_INTEGER, &
      MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_values, integer_maximum, size(integer_values), MPI_INTEGER, &
      MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. &
        any(integer_minimum /= integer_maximum)) return

    call reduce_sparse_amr_eb_patch_subtree_integral_2d( &
      distribution, sparse, level, patch, integral, ok, local_nodes)
  end subroutine &
    composite_sparse_amr_eb_patch_subtree_integral_2d

  subroutine reduce_sparse_amr_eb_patch_subtree_integral_2d( &
      distribution, sparse, level, patch, integral, ok, local_nodes)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(in) :: sparse
    integer, intent(in) :: level, patch
    real(dp), intent(out) :: integral(:)
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_nodes

    real(dp), allocatable :: local_integral(:)
    integer :: global_node_count, ierr, local_node_count
    logical :: accepted, global_ok, local_ok

    integral = 0.0_dp
    ok = .false.
    local_node_count = 0
    if (present(local_nodes)) local_nodes = 0
    allocate(local_integral(sparse%nvar), source=0.0_dp)
    call accumulate_sparse_subtree_integral_local_2d( &
      distribution, sparse, level, patch, local_integral, &
      local_node_count, local_ok)
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Allreduce( &
      local_integral, integral, sparse%nvar, MPI_DOUBLE_PRECISION, MPI_SUM, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      integral = 0.0_dp
      return
    end if
    call MPI_Allreduce( &
      local_node_count, global_node_count, 1, MPI_INTEGER, MPI_SUM, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. global_node_count < 1 .or. &
        any(.not. ieee_is_finite(integral))) then
      integral = 0.0_dp
      return
    end if

    ok = .true.
    if (present(local_nodes)) local_nodes = local_node_count
  end subroutine &
    reduce_sparse_amr_eb_patch_subtree_integral_2d

  recursive subroutine accumulate_sparse_subtree_integral_local_2d( &
      distribution, sparse, level, patch, integral, local_nodes, ok)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(in) :: sparse
    integer, intent(in) :: level, patch
    real(dp), intent(inout) :: integral(:)
    integer, intent(inout) :: local_nodes
    logical, intent(out) :: ok

    type(eb_geometry_2d) :: geometry
    logical, allocatable :: refined(:, :)
    integer :: child, first_child, i, j, last_child
    logical :: local_ok

    ok = .false.
    call topology_patch_geometry_2d( &
      sparse%topology, level - 1, patch, geometry, local_ok)
    if (.not. local_ok) return
    allocate(refined(geometry%nx, geometry%ny), source=.false.)
    first_child = 1
    last_child = 0
    if (level < sparse%level_count()) then
      first_child = sparse%topology%relations(level)% &
        child_offsets(patch) + 1
      last_child = sparse%topology%relations(level)% &
        child_offsets(patch + 1)
      do child = first_child, last_child
        refined( &
          sparse%topology%relations(level)%children(child)%patch% &
            coarse_i_lower: &
          sparse%topology%relations(level)%children(child)%patch% &
            coarse_i_upper, &
          sparse%topology%relations(level)%children(child)%patch% &
            coarse_j_lower: &
          sparse%topology%relations(level)%children(child)%patch% &
            coarse_j_upper) = .true.
      end do
    end if

    if (distribution%is_local(level - 1, patch)) then
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (refined(i, j)) cycle
          integral = integral + geometry%volume_fraction(i, j) * &
            sparse%levels(level)%patches(patch)%state(:, i, j) * &
            geometry%dx * geometry%dy
        end do
      end do
      local_nodes = local_nodes + 1
    end if
    do child = first_child, last_child
      call accumulate_sparse_subtree_integral_local_2d( &
        distribution, sparse, level + 1, child, integral, local_nodes, &
        local_ok)
      if (.not. local_ok) return
    end do
    ok = all(ieee_is_finite(integral))
  end subroutine accumulate_sparse_subtree_integral_local_2d

  subroutine allocate_sparse_tree_layout_2d( &
      distribution, nvar, topology, sparse, ok)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    integer, intent(in) :: nvar
    type(amr_eb_patch_tree_topology_2d), intent(in) :: topology
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(out) :: sparse
    logical, intent(out) :: ok

    type(eb_geometry_2d) :: geometry
    integer :: level, patch
    logical :: geometry_ok

    sparse = mpi_sparse_reactive_amr_eb_patch_tree_2d()
    ok = nvar >= 1 .and. &
      mpi_amr_eb_patch_tree_distribution_matches_2d(distribution, topology)
    if (.not. ok) return
    sparse%nvar = nvar
    sparse%topology = topology
    allocate(sparse%levels(topology%level_count()))
    do level = 1, sparse%level_count()
      allocate(sparse%levels(level)%patches( &
        topology%level_patch_count(level - 1)))
      do patch = 1, sparse%levels(level)%patch_count()
        if (.not. distribution%is_local(level - 1, patch)) cycle
        call topology_patch_geometry_2d( &
          topology, level - 1, patch, geometry, geometry_ok)
        if (.not. geometry_ok) then
          ok = .false.
          return
        end if
        allocate(sparse%levels(level)%patches(patch)%state( &
          nvar, geometry%nx, geometry%ny))
        allocate(sparse%levels(level)%patches(patch)%temperature( &
          geometry%nx, geometry%ny))
      end do
    end do
    ok = .true.
  end subroutine allocate_sparse_tree_layout_2d

  subroutine topology_patch_geometry_2d( &
      topology, level, patch, geometry, ok)
    type(amr_eb_patch_tree_topology_2d), intent(in) :: topology
    integer, intent(in) :: level, patch
    type(eb_geometry_2d), intent(out) :: geometry
    logical, intent(out) :: ok

    geometry = eb_geometry_2d()
    ok = topology%is_valid() .and. level >= 0 .and. &
      level < topology%level_count() .and. patch >= 1 .and. &
      patch <= topology%level_patch_count(level)
    if (.not. ok) return
    if (level == 0) then
      geometry = topology%root_geometry
    else
      geometry = topology%relations(level)%children(patch)%geometry
    end if
    ok = geometry%is_valid()
  end subroutine topology_patch_geometry_2d

  subroutine replicated_reactive_boundaries_match_2d( &
      boundaries, nspecies, comm, matches)
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    integer, intent(in) :: nspecies
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: matches

    real(dp), allocatable :: numeric_maximum(:), numeric_minimum(:)
    real(dp), allocatable :: numeric_values(:)
    integer, allocatable :: integer_maximum(:), integer_minimum(:)
    integer, allocatable :: integer_values(:)
    integer :: character_index, ierr, index, maximum_species
    integer :: minimum_species, nprimitive, side
    logical :: accepted, global_ok, local_ok

    matches = .false.
    call all_ranks_accept_2d( &
      comm, nspecies >= 1, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Allreduce( &
      nspecies, minimum_species, 1, MPI_INTEGER, MPI_MIN, comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      nspecies, maximum_species, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
    if (ierr /= MPI_SUCCESS .or. minimum_species /= maximum_species) return
    nprimitive = reactive_nprim(nspecies)

    call validate_reactive_boundary_set_2d(boundaries, local_ok)
    do side = 1, 4
      local_ok = local_ok .and. &
        allocated(boundaries%face(side)%inflow_primitive) .and. &
        allocated(boundaries%face(side)%prescribed_species_flux)
      if (.not. local_ok) cycle
      local_ok = size(boundaries%face(side)%inflow_primitive) == &
          nprimitive .and. &
        size(boundaries%face(side)%prescribed_species_flux) == nspecies .and. &
        all(ieee_is_finite(boundaries%face(side)%inflow_primitive)) .and. &
        all(ieee_is_finite( &
          boundaries%face(side)%prescribed_species_flux))
    end do
    call all_ranks_accept_2d( &
      comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return

    allocate(integer_values(24 * 3 * 4))
    allocate(integer_minimum(size(integer_values)))
    allocate(integer_maximum(size(integer_values)))
    integer_values = 0
    index = 1
    do side = 1, 4
      do character_index = 1, 24
        integer_values(index) = &
          iachar(boundaries%face(side)%kind(character_index:character_index))
        index = index + 1
      end do
      do character_index = 1, 24
        integer_values(index) = iachar( &
          boundaries%face(side)%thermal(character_index:character_index))
        index = index + 1
      end do
      do character_index = 1, 24
        integer_values(index) = iachar( &
          boundaries%face(side)% &
            wall_species(character_index:character_index))
        index = index + 1
      end do
    end do
    call MPI_Allreduce( &
      integer_values, integer_minimum, size(integer_values), MPI_INTEGER, &
      MPI_MIN, comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      integer_values, integer_maximum, size(integer_values), MPI_INTEGER, &
      MPI_MAX, comm, ierr)
    if (ierr /= MPI_SUCCESS .or. &
        any(integer_minimum /= integer_maximum)) return

    allocate(numeric_values(4 * (5 + nspecies + nprimitive)))
    allocate(numeric_minimum(size(numeric_values)))
    allocate(numeric_maximum(size(numeric_values)))
    index = 1
    do side = 1, 4
      numeric_values(index:index + 4) = [ &
        boundaries%face(side)%wall_temperature, &
        boundaries%face(side)%wall_velocity, &
        boundaries%face(side)%inflow_temperature]
      index = index + 5
      numeric_values(index:index + nspecies - 1) = &
        boundaries%face(side)%prescribed_species_flux
      index = index + nspecies
      numeric_values(index:index + nprimitive - 1) = &
        boundaries%face(side)%inflow_primitive
      index = index + nprimitive
    end do
    local_ok = index == size(numeric_values) + 1 .and. &
      all(ieee_is_finite(numeric_values))
    call all_ranks_accept_2d(comm, local_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Allreduce( &
      numeric_values, numeric_minimum, size(numeric_values), &
      MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      numeric_values, numeric_maximum, size(numeric_values), &
      MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
    matches = ierr == MPI_SUCCESS .and. &
      all(numeric_minimum == numeric_maximum)
  end subroutine replicated_reactive_boundaries_match_2d

  subroutine replicated_distribution_matches_2d( &
      distribution, topology, matches)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(amr_eb_patch_tree_topology_2d), intent(in) :: topology
    logical, intent(out) :: matches

    integer, allocatable :: integer_maximum(:), integer_minimum(:)
    integer, allocatable :: integer_values(:)
    integer(int64), allocatable :: work_maximum(:), work_minimum(:)
    integer(int64), allocatable :: work_values(:)
    integer :: comm_rank, comm_size, ierr, level, patch
    logical :: accepted, global_ok, local_ok

    call MPI_Comm_rank(distribution%comm, comm_rank, ierr)
    local_ok = ierr == MPI_SUCCESS
    if (local_ok) call MPI_Comm_size(distribution%comm, comm_size, ierr)
    local_ok = local_ok .and. ierr == MPI_SUCCESS .and. &
      distribution%rank == comm_rank .and. distribution%nranks == comm_size &
      .and. mpi_amr_eb_patch_tree_distribution_matches_2d( &
        distribution, topology)
    call all_ranks_accept_2d( &
      distribution%comm, local_ok, accepted, global_ok)
    matches = global_ok .and. accepted
    if (.not. matches) return

    allocate(integer_values(3 + 2 * distribution%nranks))
    allocate(integer_minimum(size(integer_values)))
    allocate(integer_maximum(size(integer_values)))
    integer_values(1:3) = [ &
      distribution%nranks, distribution%subcycle_exponent, &
      distribution%level_count()]
    integer_values(4:3 + distribution%nranks) = &
      distribution%rank_cell_counts
    integer_values(4 + distribution%nranks:) = &
      distribution%rank_patch_counts
    call MPI_Allreduce( &
      integer_values, integer_minimum, size(integer_values), MPI_INTEGER, &
      MPI_MIN, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      matches = .false.
      return
    end if
    call MPI_Allreduce( &
      integer_values, integer_maximum, size(integer_values), MPI_INTEGER, &
      MPI_MAX, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. any(integer_minimum /= integer_maximum)) then
      matches = .false.
      return
    end if
    allocate(work_values(distribution%nranks))
    allocate(work_minimum(distribution%nranks))
    allocate(work_maximum(distribution%nranks))
    work_values = distribution%rank_work_counts
    call MPI_Allreduce( &
      work_values, work_minimum, size(work_values), MPI_INTEGER8, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      matches = .false.
      return
    end if
    call MPI_Allreduce( &
      work_values, work_maximum, size(work_values), MPI_INTEGER8, MPI_MAX, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. any(work_minimum /= work_maximum)) then
      matches = .false.
      return
    end if

    deallocate(integer_values, integer_minimum, integer_maximum)
    deallocate(work_values, work_minimum, work_maximum)
    allocate(integer_values(2), integer_minimum(2), integer_maximum(2))
    allocate(work_values(1), work_minimum(1), work_maximum(1))
    do level = 1, distribution%level_count()
      do patch = 1, distribution%levels(level)%patch_count()
        integer_values = [ &
          distribution%levels(level)%owners(patch), &
          distribution%levels(level)%cell_counts(patch)]
        work_values(1) = distribution%levels(level)%work_counts(patch)
        call MPI_Allreduce( &
          integer_values, integer_minimum, 2, MPI_INTEGER, MPI_MIN, &
          distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS) then
          matches = .false.
          return
        end if
        call MPI_Allreduce( &
          integer_values, integer_maximum, 2, MPI_INTEGER, MPI_MAX, &
          distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS .or. &
            any(integer_minimum /= integer_maximum)) then
          matches = .false.
          return
        end if
        call MPI_Allreduce( &
          work_values, work_minimum, 1, MPI_INTEGER8, MPI_MIN, &
          distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS) then
          matches = .false.
          return
        end if
        call MPI_Allreduce( &
          work_values, work_maximum, 1, MPI_INTEGER8, MPI_MAX, &
          distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS .or. &
            work_minimum(1) /= work_maximum(1)) then
          matches = .false.
          return
        end if
      end do
    end do
    matches = .true.
  end subroutine replicated_distribution_matches_2d

  subroutine replicated_topology_matches_2d(topology, comm, matches)
    type(amr_eb_patch_tree_topology_2d), intent(in) :: topology
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: matches

    type(eb_geometry_2d) :: geometry
    integer :: integer_maximum(8), integer_minimum(8), integer_values(8)
    integer :: ierr, level, maximum_levels, minimum_levels, patch
    real(dp) :: numeric_maximum(7), numeric_minimum(7), numeric_values(7)
    logical :: accepted, global_ok, geometry_ok

    call all_ranks_accept_2d( &
      comm, topology%is_valid(), accepted, global_ok)
    matches = global_ok .and. accepted
    if (.not. matches) return
    call MPI_Allreduce( &
      topology%level_count(), minimum_levels, 1, MPI_INTEGER, MPI_MIN, &
      comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      matches = .false.
      return
    end if
    call MPI_Allreduce( &
      topology%level_count(), maximum_levels, 1, MPI_INTEGER, MPI_MAX, &
      comm, ierr)
    if (ierr /= MPI_SUCCESS .or. minimum_levels /= maximum_levels) then
      matches = .false.
      return
    end if

    do level = 0, topology%level_count() - 1
      integer_values = 0
      integer_values(1) = topology%level_patch_count(level)
      if (level > 0) &
        integer_values(2) = topology%relations(level)%refinement_ratio
      call MPI_Allreduce( &
        integer_values, integer_minimum, size(integer_values), MPI_INTEGER, &
        MPI_MIN, comm, ierr)
      if (ierr /= MPI_SUCCESS) then
        matches = .false.
        return
      end if
      call MPI_Allreduce( &
        integer_values, integer_maximum, size(integer_values), MPI_INTEGER, &
        MPI_MAX, comm, ierr)
      if (ierr /= MPI_SUCCESS .or. &
          any(integer_minimum /= integer_maximum)) then
        matches = .false.
        return
      end if
      do patch = 1, topology%level_patch_count(level)
        call topology_patch_geometry_2d( &
          topology, level, patch, geometry, geometry_ok)
        if (.not. geometry_ok) then
          matches = .false.
          return
        end if
        integer_values = 0
        integer_values(1:2) = [geometry%nx, geometry%ny]
        if (level > 0) then
          integer_values(3) = &
            topology%relations(level)%children(patch)%parent_patch
          integer_values(4:7) = [ &
            topology%relations(level)%children(patch)%patch%coarse_i_lower, &
            topology%relations(level)%children(patch)%patch%coarse_i_upper, &
            topology%relations(level)%children(patch)%patch%coarse_j_lower, &
            topology%relations(level)%children(patch)%patch%coarse_j_upper]
        end if
        call MPI_Allreduce( &
          integer_values, integer_minimum, size(integer_values), &
          MPI_INTEGER, MPI_MIN, comm, ierr)
        if (ierr /= MPI_SUCCESS) then
          matches = .false.
          return
        end if
        call MPI_Allreduce( &
          integer_values, integer_maximum, size(integer_values), &
          MPI_INTEGER, MPI_MAX, comm, ierr)
        if (ierr /= MPI_SUCCESS .or. &
            any(integer_minimum /= integer_maximum)) then
          matches = .false.
          return
        end if
        numeric_values = [ &
          geometry%x_lower, geometry%x_upper, geometry%y_lower, &
          geometry%y_upper, geometry%dx, geometry%dy, &
          sum(geometry%volume_fraction)]
        call MPI_Allreduce( &
          numeric_values, numeric_minimum, size(numeric_values), &
          MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
        if (ierr /= MPI_SUCCESS) then
          matches = .false.
          return
        end if
        call MPI_Allreduce( &
          numeric_values, numeric_maximum, size(numeric_values), &
          MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
        if (ierr /= MPI_SUCCESS .or. &
            any(numeric_minimum /= numeric_maximum)) then
          matches = .false.
          return
        end if
      end do
    end do
    matches = .true.
  end subroutine replicated_topology_matches_2d

  subroutine all_ranks_accept_2d( &
      comm, local_ok, accepted, mpi_ok)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(in) :: local_ok
    logical, intent(out) :: accepted, mpi_ok

    integer :: ierr

    call MPI_Allreduce( &
      local_ok, accepted, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
    mpi_ok = ierr == MPI_SUCCESS
    if (.not. mpi_ok) accepted = .false.
  end subroutine all_ranks_accept_2d

end module mpi_amr_eb_patch_tree_2d_mod
