module mpi_amr_patch_1d_mod
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: advance_reactive_chemistry
  use amr_reactive_1d_mod, only: recover_level_temperatures_1d
  use amr_patch_tree_1d_mod, only: &
    amr_patch_tree_hierarchy_1d, amr_patch_tree_level_fields_1d, &
    patch_tree_fields_are_valid_1d, average_down_patch_tree_1d
  use amr_patch_tree_reactive_1d_mod, only: &
    amr_patch_tree_reactive_patch_1d, &
    amr_patch_tree_reactive_solution_1d, refresh_patch_tree_ghosts
  implicit none
  private

  type, public :: mpi_amr_patch_level_ownership_1d
    integer, allocatable :: owners(:)
    integer, allocatable :: cell_counts(:)
  contains
    procedure :: patch_count => mpi_amr_level_patch_count
    procedure :: is_valid => mpi_amr_level_ownership_is_valid
  end type mpi_amr_patch_level_ownership_1d

  type, public :: mpi_amr_patch_distribution_1d
    type(MPI_Comm) :: comm = MPI_COMM_NULL
    integer :: rank = -1
    integer :: nranks = 0
    type(mpi_amr_patch_level_ownership_1d), allocatable :: levels(:)
    integer, allocatable :: rank_cell_counts(:)
    integer, allocatable :: rank_patch_counts(:)
  contains
    procedure :: level_count => mpi_amr_distribution_level_count
    procedure :: owner_of => mpi_amr_distribution_owner_of
    procedure :: is_local => mpi_amr_distribution_is_local
    procedure :: is_valid => mpi_amr_distribution_is_valid
  end type mpi_amr_patch_distribution_1d

  type, public :: mpi_amr_patch_halo_1d
    logical :: has_left = .false.
    logical :: has_right = .false.
    real(dp), allocatable :: left(:, :)
    real(dp), allocatable :: right(:, :)
  end type mpi_amr_patch_halo_1d

  type, public :: mpi_amr_level_halos_1d
    type(mpi_amr_patch_halo_1d), allocatable :: patches(:)
  end type mpi_amr_level_halos_1d

  public :: initialize_mpi_amr_patch_distribution_1d
  public :: mpi_amr_distribution_matches_hierarchy_1d
  public :: synchronize_owned_patch_tree_fields_1d
  public :: exchange_owned_adjacent_patch_halos_1d
  public :: synchronize_owned_patch_tree_reactive_1d
  public :: advance_owned_patch_tree_chemistry_1d

contains

  pure integer function mpi_amr_level_patch_count(self) result(count)
    class(mpi_amr_patch_level_ownership_1d), intent(in) :: self

    count = 0
    if (allocated(self%owners)) count = size(self%owners)
  end function mpi_amr_level_patch_count

  pure logical function mpi_amr_level_ownership_is_valid( &
      self, nranks) result(valid)
    class(mpi_amr_patch_level_ownership_1d), intent(in) :: self
    integer, intent(in) :: nranks

    valid = nranks >= 1 .and. allocated(self%owners) .and. &
      allocated(self%cell_counts)
    if (.not. valid) return
    valid = size(self%owners) >= 1 .and. &
      size(self%cell_counts) == size(self%owners) .and. &
      all(self%owners >= 0) .and. all(self%owners < nranks) .and. &
      all(self%cell_counts >= 1)
  end function mpi_amr_level_ownership_is_valid

  pure integer function mpi_amr_distribution_level_count(self) result(count)
    class(mpi_amr_patch_distribution_1d), intent(in) :: self

    count = 0
    if (allocated(self%levels)) count = size(self%levels)
  end function mpi_amr_distribution_level_count

  pure integer function mpi_amr_distribution_owner_of( &
      self, level, patch) result(owner)
    class(mpi_amr_patch_distribution_1d), intent(in) :: self
    integer, intent(in) :: level, patch

    owner = -1
    if (.not. allocated(self%levels)) return
    if (level < 0 .or. level >= size(self%levels)) return
    if (.not. allocated(self%levels(level + 1)%owners)) return
    if (patch < 1 .or. &
        patch > size(self%levels(level + 1)%owners)) return
    owner = self%levels(level + 1)%owners(patch)
  end function mpi_amr_distribution_owner_of

  pure logical function mpi_amr_distribution_is_local( &
      self, level, patch) result(local)
    class(mpi_amr_patch_distribution_1d), intent(in) :: self
    integer, intent(in) :: level, patch

    local = self%owner_of(level, patch) == self%rank
  end function mpi_amr_distribution_is_local

  pure logical function mpi_amr_distribution_is_valid(self) result(valid)
    class(mpi_amr_patch_distribution_1d), intent(in) :: self

    integer, allocatable :: cells(:), patches(:)
    integer :: level, patch, owner

    valid = self%rank >= 0 .and. self%nranks >= 1 .and. &
      self%rank < self%nranks .and. allocated(self%levels) .and. &
      allocated(self%rank_cell_counts) .and. &
      allocated(self%rank_patch_counts)
    if (.not. valid) return
    valid = size(self%levels) >= 1 .and. &
      size(self%rank_cell_counts) == self%nranks .and. &
      size(self%rank_patch_counts) == self%nranks
    if (.not. valid) return

    allocate(cells(self%nranks), patches(self%nranks))
    cells = 0
    patches = 0
    do level = 1, size(self%levels)
      valid = self%levels(level)%is_valid(self%nranks)
      if (.not. valid) return
      do patch = 1, self%levels(level)%patch_count()
        owner = self%levels(level)%owners(patch) + 1
        cells(owner) = cells(owner) + &
          self%levels(level)%cell_counts(patch)
        patches(owner) = patches(owner) + 1
      end do
    end do
    valid = all(cells == self%rank_cell_counts) .and. &
      all(patches == self%rank_patch_counts) .and. &
      sum(self%rank_patch_counts) >= 1
  end function mpi_amr_distribution_is_valid

  subroutine initialize_mpi_amr_patch_distribution_1d( &
      hierarchy, comm, distribution, ok)
    type(amr_patch_tree_hierarchy_1d), intent(in) :: hierarchy
    type(MPI_Comm), intent(in) :: comm
    type(mpi_amr_patch_distribution_1d), intent(out) :: distribution
    logical, intent(out) :: ok

    logical :: local_ok
    integer :: ierr, level, patch, owner

    ok = .false.
    distribution%comm = comm
    call MPI_Comm_rank(comm, distribution%rank, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Comm_size(comm, distribution%nranks, ierr)
    if (ierr /= MPI_SUCCESS .or. distribution%nranks < 1) return

    call replicated_hierarchy_matches_1d(hierarchy, comm, local_ok)
    if (.not. local_ok) return
    allocate(distribution%levels(hierarchy%level_count()))
    allocate(distribution%rank_cell_counts(distribution%nranks))
    allocate(distribution%rank_patch_counts(distribution%nranks))
    distribution%rank_cell_counts = 0
    distribution%rank_patch_counts = 0

    do level = 0, hierarchy%level_count() - 1
      allocate(distribution%levels(level + 1)%owners( &
        hierarchy%level_patch_count(level)))
      allocate(distribution%levels(level + 1)%cell_counts( &
        hierarchy%level_patch_count(level)))
      do patch = 1, hierarchy%level_patch_count(level)
        distribution%levels(level + 1)%cell_counts(patch) = &
          patch_cell_count_1d(hierarchy, level, patch)
        owner = minloc(distribution%rank_cell_counts, dim=1)
        distribution%levels(level + 1)%owners(patch) = owner - 1
        distribution%rank_cell_counts(owner) = &
          distribution%rank_cell_counts(owner) + &
          distribution%levels(level + 1)%cell_counts(patch)
        distribution%rank_patch_counts(owner) = &
          distribution%rank_patch_counts(owner) + 1
      end do
    end do
    ok = distribution%is_valid() .and. &
      mpi_amr_distribution_matches_hierarchy_1d(distribution, hierarchy)
  end subroutine initialize_mpi_amr_patch_distribution_1d

  pure logical function mpi_amr_distribution_matches_hierarchy_1d( &
      distribution, hierarchy) result(matches)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(amr_patch_tree_hierarchy_1d), intent(in) :: hierarchy

    integer :: level, patch

    matches = distribution%is_valid() .and. hierarchy%is_valid()
    if (.not. matches) return
    matches = distribution%level_count() == hierarchy%level_count()
    if (.not. matches) return
    do level = 0, hierarchy%level_count() - 1
      matches = distribution%levels(level + 1)%patch_count() == &
        hierarchy%level_patch_count(level)
      if (.not. matches) return
      do patch = 1, hierarchy%level_patch_count(level)
        matches = distribution%levels(level + 1)%cell_counts(patch) == &
          patch_cell_count_1d(hierarchy, level, patch)
        if (.not. matches) return
      end do
    end do
  end function mpi_amr_distribution_matches_hierarchy_1d

  subroutine synchronize_owned_patch_tree_fields_1d( &
      distribution, hierarchy, fields, ok)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(amr_patch_tree_hierarchy_1d), intent(in) :: hierarchy
    type(amr_patch_tree_level_fields_1d), intent(inout) :: fields(:)
    logical, intent(out) :: ok

    integer :: ierr, level, patch, owner

    ok = mpi_amr_distribution_matches_hierarchy_1d( &
      distribution, hierarchy) .and. &
      patch_tree_fields_are_valid_1d(fields, hierarchy)
    if (.not. ok) return
    do level = 1, size(fields)
      do patch = 1, size(fields(level)%patches)
        owner = distribution%owner_of(level - 1, patch)
        call MPI_Bcast( &
          fields(level)%patches(patch)%values, &
          size(fields(level)%patches(patch)%values), &
          MPI_DOUBLE_PRECISION, owner, distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS) then
          ok = .false.
          return
        end if
      end do
    end do
    ok = patch_tree_fields_are_valid_1d(fields, hierarchy)
  end subroutine synchronize_owned_patch_tree_fields_1d

  subroutine exchange_owned_adjacent_patch_halos_1d( &
      distribution, hierarchy, fields, ghost_width, halos, ok)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(amr_patch_tree_hierarchy_1d), intent(in) :: hierarchy
    type(amr_patch_tree_level_fields_1d), intent(in) :: fields(:)
    integer, intent(in) :: ghost_width
    type(mpi_amr_level_halos_1d), allocatable, intent(out) :: halos(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: buffer(:, :)
    logical :: adjacent
    integer :: ierr, nvar, level, patch, parent, child
    integer :: left_patch, right_patch, left_nx, right_nx, layer, owner

    ok = ghost_width >= 1 .and. &
      mpi_amr_distribution_matches_hierarchy_1d( &
        distribution, hierarchy) .and. &
      patch_tree_fields_are_valid_1d(fields, hierarchy)
    if (.not. ok) return
    nvar = size(fields(1)%patches(1)%values, 1)
    allocate(halos(hierarchy%level_count()))
    do level = 1, hierarchy%level_count()
      allocate(halos(level)%patches(hierarchy%level_patch_count(level - 1)))
      do patch = 1, size(halos(level)%patches)
        allocate(halos(level)%patches(patch)%left(nvar, ghost_width))
        allocate(halos(level)%patches(patch)%right(nvar, ghost_width))
        halos(level)%patches(patch)%left = 0.0_dp
        halos(level)%patches(patch)%right = 0.0_dp
      end do
    end do
    allocate(buffer(nvar, ghost_width))

    do level = 1, size(hierarchy%relations)
      do parent = 1, hierarchy%relations(level)%parent_patch_count()
        do child = 1, hierarchy%relations(level)% &
            child_sets(parent)%patch_count() - 1
          adjacent = hierarchy%relations(level)%child_sets(parent)% &
            patches(child)%fine%upper + 1 == hierarchy%relations(level)% &
            child_sets(parent)%patches(child + 1)%fine%lower
          if (.not. adjacent) cycle
          left_patch = hierarchy%relations(level)% &
            child_index(parent, child)
          right_patch = hierarchy%relations(level)% &
            child_index(parent, child + 1)
          left_nx = size(fields(level + 1)%patches(left_patch)%values, 2)
          right_nx = size(fields(level + 1)%patches(right_patch)%values, 2)
          if (left_nx < ghost_width .or. right_nx < ghost_width) then
            ok = .false.
            return
          end if

          buffer = 0.0_dp
          owner = distribution%owner_of(level, left_patch)
          if (distribution%rank == owner) then
            do layer = 1, ghost_width
              buffer(:, layer) = fields(level + 1)%patches(left_patch)% &
                values(:, left_nx - layer + 1)
            end do
          end if
          call MPI_Bcast(buffer, size(buffer), MPI_DOUBLE_PRECISION, owner, &
            distribution%comm, ierr)
          if (ierr /= MPI_SUCCESS) then
            ok = .false.
            return
          end if
          halos(level + 1)%patches(right_patch)%left = buffer
          halos(level + 1)%patches(right_patch)%has_left = .true.

          buffer = 0.0_dp
          owner = distribution%owner_of(level, right_patch)
          if (distribution%rank == owner) then
            do layer = 1, ghost_width
              buffer(:, layer) = fields(level + 1)%patches(right_patch)% &
                values(:, layer)
            end do
          end if
          call MPI_Bcast(buffer, size(buffer), MPI_DOUBLE_PRECISION, owner, &
            distribution%comm, ierr)
          if (ierr /= MPI_SUCCESS) then
            ok = .false.
            return
          end if
          halos(level + 1)%patches(left_patch)%right = buffer
          halos(level + 1)%patches(left_patch)%has_right = .true.
        end do
      end do
    end do
    ok = .true.
  end subroutine exchange_owned_adjacent_patch_halos_1d

  subroutine synchronize_owned_patch_tree_reactive_1d( &
      distribution, solution, ok)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(amr_patch_tree_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok

    logical :: local_ok, accepted, mpi_ok
    integer :: level, patch, owner, nvar, minimum_nvar, maximum_nvar, ierr

    local_ok = solution%is_valid() .and. &
      mpi_amr_distribution_matches_hierarchy_1d( &
        distribution, solution%hierarchy)
    call all_ranks_accept_1d( &
      distribution, local_ok, accepted, mpi_ok)
    ok = mpi_ok .and. accepted
    if (.not. ok) return
    nvar = size(solution%levels(1)%patches(1)%state, 1)
    call MPI_Allreduce( &
      nvar, minimum_nvar, 1, MPI_INTEGER, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Allreduce( &
      nvar, maximum_nvar, 1, MPI_INTEGER, MPI_MAX, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. minimum_nvar /= maximum_nvar) then
      ok = .false.
      return
    end if
    do level = 1, solution%level_count()
      do patch = 1, size(solution%levels(level)%patches)
        owner = distribution%owner_of(level - 1, patch)
        call broadcast_owned_reactive_patch_1d( &
          distribution, owner, &
          solution%levels(level)%patches(patch), local_ok)
        if (.not. local_ok) then
          ok = .false.
          return
        end if
      end do
    end do
    ok = solution%is_valid()
  end subroutine synchronize_owned_patch_tree_reactive_1d

  subroutine advance_owned_patch_tree_chemistry_1d( &
      species, reactions, config, interval, distribution, solution, ok, &
      local_patch_advances)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: interval
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(amr_patch_tree_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_patch_advances

    type(amr_patch_tree_reactive_solution_1d) :: backup
    character(len=32) :: boundary
    logical :: local_ok, patch_ok, accepted, mpi_ok
    integer :: level, patch, owner, nx, advances

    ok = .false.
    advances = 0
    if (present(local_patch_advances)) local_patch_advances = 0
    local_ok = interval >= 0.0_dp .and. size(species) >= 1 .and. &
      solution%is_valid() .and. &
      mpi_amr_distribution_matches_hierarchy_1d( &
        distribution, solution%hierarchy)
    call all_ranks_accept_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return

    call synchronize_owned_patch_tree_reactive_1d( &
      distribution, solution, local_ok)
    call all_ranks_accept_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
    backup = solution

    do level = 1, solution%level_count()
      boundary = "outflow"
      if (level == 1) boundary = config%boundary_condition
      do patch = 1, size(solution%levels(level)%patches)
        owner = distribution%owner_of(level - 1, patch)
        patch_ok = .true.
        if (distribution%rank == owner) then
          nx = size(solution%levels(level)%patches(patch)%state, 2) - 2
          call advance_reactive_chemistry( &
            species, reactions, &
            solution%levels(level)%patches(patch)%state, &
            solution%levels(level)%patches(patch)%temperature, nx, interval, &
            config%chemistry_relative_tolerance, &
            config%chemistry_absolute_tolerance, boundary, patch_ok)
          if (patch_ok) advances = advances + 1
        end if
        call all_ranks_accept_1d( &
          distribution, patch_ok, accepted, mpi_ok)
        if (.not. mpi_ok .or. .not. accepted) then
          solution = backup
          advances = 0
          return
        end if
        call broadcast_owned_reactive_patch_1d( &
          distribution, owner, &
          solution%levels(level)%patches(patch), local_ok)
        call all_ranks_accept_1d( &
          distribution, local_ok, accepted, mpi_ok)
        if (.not. mpi_ok .or. .not. accepted) then
          solution = backup
          advances = 0
          return
        end if
      end do
    end do

    call average_down_reactive_solution_1d( &
      species, config, solution, local_ok)
    call all_ranks_accept_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) then
      solution = backup
      advances = 0
      return
    end if
    local_ok = solution%is_valid()
    call all_ranks_accept_1d( &
      distribution, local_ok, accepted, mpi_ok)
    ok = mpi_ok .and. accepted
    if (.not. ok) then
      solution = backup
      advances = 0
      return
    end if
    if (present(local_patch_advances)) local_patch_advances = advances
  end subroutine advance_owned_patch_tree_chemistry_1d

  subroutine broadcast_owned_reactive_patch_1d( &
      distribution, owner, patch, ok)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    integer, intent(in) :: owner
    type(amr_patch_tree_reactive_patch_1d), intent(inout) :: patch
    logical, intent(out) :: ok

    integer :: ierr

    ok = owner >= 0 .and. owner < distribution%nranks
    if (.not. ok) return
    call MPI_Bcast(patch%state, size(patch%state), MPI_DOUBLE_PRECISION, &
      owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast( &
      patch%temperature, size(patch%temperature), MPI_DOUBLE_PRECISION, &
      owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast( &
      patch%left_ghost_state, size(patch%left_ghost_state), &
      MPI_DOUBLE_PRECISION, owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast( &
      patch%right_ghost_state, size(patch%right_ghost_state), &
      MPI_DOUBLE_PRECISION, owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast( &
      patch%left_ghost_temperature, &
      size(patch%left_ghost_temperature), MPI_DOUBLE_PRECISION, owner, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast( &
      patch%right_ghost_temperature, &
      size(patch%right_ghost_temperature), MPI_DOUBLE_PRECISION, owner, &
      distribution%comm, ierr)
    ok = ierr == MPI_SUCCESS
  end subroutine broadcast_owned_reactive_patch_1d

  subroutine average_down_reactive_solution_1d( &
      species, config, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_patch_tree_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok

    type(amr_patch_tree_level_fields_1d), allocatable :: fields(:)
    logical :: local_ok
    integer :: level, patch, nx

    ok = solution%is_valid()
    if (.not. ok) return
    allocate(fields(solution%level_count()))
    do level = 1, solution%level_count()
      allocate(fields(level)%patches(size(solution%levels(level)%patches)))
      do patch = 1, size(solution%levels(level)%patches)
        nx = size(solution%levels(level)%patches(patch)%state, 2) - 2
        fields(level)%patches(patch)%values = &
          solution%levels(level)%patches(patch)%state(:, 1:nx)
      end do
    end do
    call average_down_patch_tree_1d( &
      fields, solution%hierarchy, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    do level = 1, solution%level_count()
      do patch = 1, size(solution%levels(level)%patches)
        nx = size(solution%levels(level)%patches(patch)%state, 2) - 2
        solution%levels(level)%patches(patch)%state(:, 1:nx) = &
          fields(level)%patches(patch)%values
        call recover_level_temperatures_1d( &
          species, solution%levels(level)%patches(patch)%state, &
          solution%levels(level)%patches(patch)%temperature, nx, local_ok)
        if (.not. local_ok) then
          ok = .false.
          return
        end if
      end do
    end do
    call refresh_patch_tree_ghosts(species, config, solution, ok)
  end subroutine average_down_reactive_solution_1d

  subroutine all_ranks_accept_1d( &
      distribution, local_value, global_value, ok)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    logical, intent(in) :: local_value
    logical, intent(out) :: global_value, ok

    integer :: ierr

    call MPI_Allreduce( &
      local_value, global_value, 1, MPI_LOGICAL, MPI_LAND, &
      distribution%comm, ierr)
    ok = ierr == MPI_SUCCESS
    if (.not. ok) global_value = .false.
  end subroutine all_ranks_accept_1d

  pure integer function patch_cell_count_1d( &
      hierarchy, level, patch) result(count)
    type(amr_patch_tree_hierarchy_1d), intent(in) :: hierarchy
    integer, intent(in) :: level, patch

    integer :: parent, local_child, first_child, last_child

    count = 0
    if (level < 0 .or. level >= hierarchy%level_count()) return
    if (patch < 1 .or. patch > hierarchy%level_patch_count(level)) return
    if (level == 0) then
      count = hierarchy%base_cells
      return
    end if
    do parent = 1, hierarchy%relations(level)%parent_patch_count()
      first_child = hierarchy%relations(level)%child_offsets(parent) + 1
      last_child = hierarchy%relations(level)%child_offsets(parent + 1)
      if (patch < first_child .or. patch > last_child) cycle
      local_child = patch - hierarchy%relations(level)%child_offsets(parent)
      count = hierarchy%relations(level)%child_sets(parent)% &
        patches(local_child)%fine%cell_count()
      return
    end do
  end function patch_cell_count_1d

  subroutine replicated_hierarchy_matches_1d(hierarchy, comm, matches)
    type(amr_patch_tree_hierarchy_1d), intent(in) :: hierarchy
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: matches

    integer, allocatable :: metadata(:), minimum(:), maximum(:)
    real(dp) :: geometry(2), minimum_geometry(2), maximum_geometry(2)
    logical :: local_valid, global_valid
    integer :: ierr, metadata_size, minimum_size, maximum_size

    matches = .false.
    local_valid = hierarchy%is_valid()
    call MPI_Allreduce(local_valid, global_valid, 1, MPI_LOGICAL, MPI_LAND, &
      comm, ierr)
    if (ierr /= MPI_SUCCESS .or. .not. global_valid) return
    metadata_size = hierarchy_metadata_size_1d(hierarchy)
    call MPI_Allreduce(metadata_size, minimum_size, 1, MPI_INTEGER, MPI_MIN, &
      comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce(metadata_size, maximum_size, 1, MPI_INTEGER, MPI_MAX, &
      comm, ierr)
    if (ierr /= MPI_SUCCESS .or. minimum_size /= maximum_size) return

    allocate(metadata(metadata_size), minimum(metadata_size), &
      maximum(metadata_size))
    call pack_hierarchy_metadata_1d(hierarchy, metadata)
    call MPI_Allreduce(metadata, minimum, metadata_size, MPI_INTEGER, MPI_MIN, &
      comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce(metadata, maximum, metadata_size, MPI_INTEGER, MPI_MAX, &
      comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    geometry = [hierarchy%x_lower, hierarchy%x_upper]
    call MPI_Allreduce(geometry, minimum_geometry, 2, MPI_DOUBLE_PRECISION, &
      MPI_MIN, comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce(geometry, maximum_geometry, 2, MPI_DOUBLE_PRECISION, &
      MPI_MAX, comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    matches = all(metadata == minimum) .and. all(metadata == maximum) .and. &
      all(geometry == minimum_geometry) .and. &
      all(geometry == maximum_geometry)
  end subroutine replicated_hierarchy_matches_1d

  pure integer function hierarchy_metadata_size_1d(hierarchy) result(count)
    type(amr_patch_tree_hierarchy_1d), intent(in) :: hierarchy

    integer :: relation, parent

    count = 2
    do relation = 1, size(hierarchy%relations)
      count = count + 3
      do parent = 1, hierarchy%relations(relation)%parent_patch_count()
        count = count + 1 + 3 * hierarchy%relations(relation)% &
          child_sets(parent)%patch_count()
      end do
    end do
  end function hierarchy_metadata_size_1d

  pure subroutine pack_hierarchy_metadata_1d(hierarchy, metadata)
    type(amr_patch_tree_hierarchy_1d), intent(in) :: hierarchy
    integer, intent(out) :: metadata(:)

    integer :: relation, parent, child, offset

    metadata = 0
    offset = 2
    metadata(1) = hierarchy%base_cells
    metadata(2) = hierarchy%level_count()
    do relation = 1, size(hierarchy%relations)
      metadata(offset + 1) = hierarchy%relations(relation)%refinement_ratio
      metadata(offset + 2) = hierarchy%relations(relation)%parent_patch_count()
      metadata(offset + 3) = hierarchy%relations(relation)%child_patch_count()
      offset = offset + 3
      do parent = 1, hierarchy%relations(relation)%parent_patch_count()
        metadata(offset + 1) = hierarchy%relations(relation)% &
          child_sets(parent)%patch_count()
        offset = offset + 1
        do child = 1, hierarchy%relations(relation)% &
            child_sets(parent)%patch_count()
          metadata(offset + 1) = parent
          metadata(offset + 2) = hierarchy%relations(relation)% &
            child_sets(parent)%patches(child)%fine_coarse_lower
          metadata(offset + 3) = hierarchy%relations(relation)% &
            child_sets(parent)%patches(child)%fine_coarse_upper
          offset = offset + 3
        end do
      end do
    end do
  end subroutine pack_hierarchy_metadata_1d

end module mpi_amr_patch_1d_mod
