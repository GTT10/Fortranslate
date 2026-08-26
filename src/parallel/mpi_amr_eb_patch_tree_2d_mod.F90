module mpi_amr_eb_patch_tree_2d_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use mpi_f08
  use precision_mod, only: dp
  use eb_geometry_2d_mod, only: eb_geometry_2d
  use amr_eb_patch_tree_2d_mod, only: amr_eb_patch_tree_topology_2d
  use amr_eb_patch_tree_reactive_2d_mod, only: &
    reactive_amr_eb_patch_tree_2d
  implicit none
  private

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

  public :: initialize_mpi_amr_eb_patch_tree_distribution_2d
  public :: mpi_amr_eb_patch_tree_distribution_matches_2d
  public :: synchronize_owned_reactive_amr_eb_patch_tree_2d

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
