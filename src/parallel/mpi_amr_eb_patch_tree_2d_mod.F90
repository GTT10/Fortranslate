module mpi_amr_eb_patch_tree_2d_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: &
    gas_transport_species, compatible_transport_database
  use reactive_1d_mod, only: reactive_nvar
  use reactive_2d_mod, only: advance_reactive_chemistry_2d
  use reactive_eb_cfl_2d_mod, only: compute_reactive_eb_cfl_timestep_2d
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_covered_cell
  use eb_reactive_transport_2d_mod, only: &
    reactive_eb_transport_timestep_2d
  use amr_eb_hierarchy_2d_mod, only: &
    average_down_reactive_eb_state_patch_2d
  use amr_eb_transport_2d_mod, only: recover_transport_temperature_2d
  use amr_eb_patch_tree_2d_mod, only: amr_eb_patch_tree_topology_2d
  use amr_eb_patch_tree_reactive_2d_mod, only: &
    reactive_amr_eb_patch_tree_2d
  implicit none
  private

  integer, parameter :: sparse_tree_state_tag = 27101
  integer, parameter :: sparse_tree_temperature_tag = 27102
  integer, parameter :: sparse_tree_restriction_tag = 27103

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

  public :: initialize_mpi_amr_eb_patch_tree_distribution_2d
  public :: mpi_amr_eb_patch_tree_distribution_matches_2d
  public :: synchronize_owned_reactive_amr_eb_patch_tree_2d
  public :: initialize_sparse_owned_reactive_amr_eb_patch_tree_2d
  public :: materialize_sparse_owned_reactive_amr_eb_patch_tree_2d
  public :: migrate_sparse_owned_reactive_amr_eb_patch_tree_2d
  public :: compute_sparse_owned_reactive_amr_eb_patch_tree_timestep_2d
  public :: advance_sparse_owned_reactive_amr_eb_patch_tree_chemistry_2d
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

    real(dp), allocatable :: local_integral(:)
    integer :: global_node_count, ierr, integer_maximum(3)
    integer :: integer_minimum(3), integer_values(3), local_node_count
    logical :: accepted, global_ok, local_ok, matches

    integral = 0.0_dp
    ok = .false.
    local_node_count = 0
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
    composite_sparse_amr_eb_patch_subtree_integral_2d

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
