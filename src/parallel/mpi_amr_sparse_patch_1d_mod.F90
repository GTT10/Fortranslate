module mpi_amr_sparse_patch_1d_mod
  use mpi_f08
  use precision_mod, only: dp
  use amr_patch_tree_1d_mod, only: amr_patch_tree_hierarchy_1d
  use amr_patch_tree_reactive_1d_mod, only: &
    amr_patch_tree_reactive_patch_1d, &
    amr_patch_tree_reactive_solution_1d
  use mpi_amr_patch_1d_mod, only: &
    mpi_amr_patch_distribution_1d, &
    mpi_amr_distribution_matches_hierarchy_1d, &
    synchronize_owned_patch_tree_reactive_1d
  implicit none
  private

  type, public :: mpi_amr_sparse_reactive_level_1d
    type(amr_patch_tree_reactive_patch_1d), allocatable :: patches(:)
    logical, allocatable :: is_local(:)
  end type mpi_amr_sparse_reactive_level_1d

  type, public :: mpi_amr_sparse_reactive_solution_1d
    type(amr_patch_tree_hierarchy_1d) :: hierarchy
    type(mpi_amr_sparse_reactive_level_1d), allocatable :: levels(:)
    integer, allocatable :: level_advances(:)
    integer, allocatable :: transport_level_advances(:)
    integer :: rank = -1
    integer :: nranks = 0
    integer :: nvar = 0
    integer :: ghost_width = 0
    real(dp) :: time = 0.0_dp
    integer :: steps = 0
    integer :: regrid_evaluations = 0
    integer :: regrids = 0
    integer :: overlap_cells_transferred = 0
  contains
    procedure :: is_valid => mpi_amr_sparse_reactive_is_valid
    procedure :: local_patch_count => mpi_amr_sparse_local_patch_count
    procedure :: local_cell_count => mpi_amr_sparse_local_cell_count
    procedure :: local_value_count => mpi_amr_sparse_local_value_count
  end type mpi_amr_sparse_reactive_solution_1d

  public :: scatter_owned_patch_tree_reactive_1d
  public :: gather_owned_patch_tree_reactive_1d
  public :: migrate_owned_patch_tree_reactive_1d

contains

  logical function mpi_amr_sparse_reactive_is_valid( &
      self, distribution) result(valid)
    class(mpi_amr_sparse_reactive_solution_1d), intent(in) :: self
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution

    logical :: local
    integer :: level, patch, nx

    valid = self%rank == distribution%rank .and. &
      self%nranks == distribution%nranks .and. self%nvar >= 1 .and. &
      self%ghost_width >= 1 .and. allocated(self%levels) .and. &
      allocated(self%level_advances) .and. &
      allocated(self%transport_level_advances) .and. &
      mpi_amr_distribution_matches_hierarchy_1d( &
        distribution, self%hierarchy)
    if (.not. valid) return
    valid = size(self%levels) == self%hierarchy%level_count() .and. &
      size(self%level_advances) == size(self%levels) .and. &
      size(self%transport_level_advances) == size(self%levels) .and. &
      all(self%level_advances >= 0) .and. &
      all(self%transport_level_advances >= 0) .and. self%time >= 0.0_dp .and. &
      self%steps >= 0 .and. self%regrid_evaluations >= 0 .and. &
      self%regrids >= 0 .and. self%overlap_cells_transferred >= 0
    if (.not. valid) return
    do level = 1, size(self%levels)
      valid = allocated(self%levels(level)%patches) .and. &
        allocated(self%levels(level)%is_local)
      if (.not. valid) return
      valid = size(self%levels(level)%patches) == &
        self%hierarchy%level_patch_count(level - 1) .and. &
        size(self%levels(level)%is_local) == &
          self%hierarchy%level_patch_count(level - 1)
      if (.not. valid) return
      do patch = 1, size(self%levels(level)%patches)
        local = distribution%is_local(level - 1, patch)
        valid = self%levels(level)%is_local(patch) .eqv. local
        if (.not. valid) return
        nx = distribution%levels(level)%cell_counts(patch)
        if (local) then
          valid = sparse_patch_has_shape( &
            self%levels(level)%patches(patch), self%nvar, nx, &
            self%ghost_width)
        else
          valid = sparse_patch_is_empty(self%levels(level)%patches(patch))
        end if
        if (.not. valid) return
      end do
    end do
    valid = self%local_patch_count() == &
      distribution%rank_patch_counts(self%rank + 1) .and. &
      self%local_cell_count() == &
        distribution%rank_cell_counts(self%rank + 1)
  end function mpi_amr_sparse_reactive_is_valid

  pure integer function mpi_amr_sparse_local_patch_count(self) result(patch_count)
    class(mpi_amr_sparse_reactive_solution_1d), intent(in) :: self

    integer :: level

    patch_count = 0
    if (.not. allocated(self%levels)) return
    do level = 1, size(self%levels)
      if (allocated(self%levels(level)%is_local)) &
        patch_count = patch_count + count(self%levels(level)%is_local)
    end do
  end function mpi_amr_sparse_local_patch_count

  pure integer function mpi_amr_sparse_local_cell_count(self) result(count)
    class(mpi_amr_sparse_reactive_solution_1d), intent(in) :: self

    integer :: level, patch

    count = 0
    if (.not. allocated(self%levels)) return
    do level = 1, size(self%levels)
      if (.not. allocated(self%levels(level)%patches)) cycle
      do patch = 1, size(self%levels(level)%patches)
        if (.not. allocated(self%levels(level)%patches(patch)%state)) cycle
        count = count + size(self%levels(level)%patches(patch)%state, 2) - 2
      end do
    end do
  end function mpi_amr_sparse_local_cell_count

  pure integer function mpi_amr_sparse_local_value_count(self) result(count)
    class(mpi_amr_sparse_reactive_solution_1d), intent(in) :: self

    integer :: level, patch

    count = 0
    if (.not. allocated(self%levels)) return
    do level = 1, size(self%levels)
      if (.not. allocated(self%levels(level)%patches)) cycle
      do patch = 1, size(self%levels(level)%patches)
        if (.not. allocated(self%levels(level)%patches(patch)%state)) cycle
        count = count + sparse_patch_value_count( &
          self%levels(level)%patches(patch))
      end do
    end do
  end function mpi_amr_sparse_local_value_count

  subroutine scatter_owned_patch_tree_reactive_1d( &
      distribution, replicated, sparse, ok)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(amr_patch_tree_reactive_solution_1d), intent(inout) :: replicated
    type(mpi_amr_sparse_reactive_solution_1d), intent(out) :: sparse
    logical, intent(out) :: ok

    logical :: local_ok, accepted, mpi_ok
    integer :: level, patch

    ok = .false.
    local_ok = replicated%is_valid() .and. &
      mpi_amr_distribution_matches_hierarchy_1d( &
        distribution, replicated%hierarchy)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
    call synchronize_owned_patch_tree_reactive_1d( &
      distribution, replicated, local_ok)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return

    sparse%hierarchy = replicated%hierarchy
    sparse%rank = distribution%rank
    sparse%nranks = distribution%nranks
    sparse%nvar = size(replicated%levels(1)%patches(1)%state, 1)
    sparse%ghost_width = &
      size(replicated%levels(1)%patches(1)%left_ghost_state, 2)
    sparse%level_advances = replicated%level_advances
    sparse%transport_level_advances = replicated%transport_level_advances
    sparse%time = replicated%time
    sparse%steps = replicated%steps
    sparse%regrid_evaluations = replicated%regrid_evaluations
    sparse%regrids = replicated%regrids
    sparse%overlap_cells_transferred = replicated%overlap_cells_transferred
    allocate(sparse%levels(replicated%level_count()))
    do level = 1, size(sparse%levels)
      allocate(sparse%levels(level)%patches( &
        replicated%hierarchy%level_patch_count(level - 1)))
      allocate(sparse%levels(level)%is_local( &
        replicated%hierarchy%level_patch_count(level - 1)))
      do patch = 1, size(sparse%levels(level)%patches)
        sparse%levels(level)%is_local(patch) = &
          distribution%is_local(level - 1, patch)
        if (sparse%levels(level)%is_local(patch)) &
          sparse%levels(level)%patches(patch) = &
            replicated%levels(level)%patches(patch)
      end do
    end do
    local_ok = sparse%is_valid(distribution)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    ok = mpi_ok .and. accepted
    if (.not. ok) sparse = mpi_amr_sparse_reactive_solution_1d()
  end subroutine scatter_owned_patch_tree_reactive_1d

  subroutine gather_owned_patch_tree_reactive_1d( &
      distribution, sparse, replicated, ok)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(in) :: sparse
    type(amr_patch_tree_reactive_solution_1d), intent(inout) :: replicated
    logical, intent(out) :: ok

    logical :: local_ok, accepted, mpi_ok
    integer :: level, patch, root_owner

    ok = .false.
    local_ok = sparse%is_valid(distribution) .and. replicated%is_valid() .and. &
      mpi_amr_distribution_matches_hierarchy_1d( &
        distribution, replicated%hierarchy)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
    do level = 1, sparse%hierarchy%level_count()
      do patch = 1, sparse%hierarchy%level_patch_count(level - 1)
        if (distribution%is_local(level - 1, patch)) &
          replicated%levels(level)%patches(patch) = &
            sparse%levels(level)%patches(patch)
      end do
    end do
    root_owner = distribution%owner_of(0, 1)
    if (distribution%rank == root_owner) then
      replicated%level_advances = sparse%level_advances
      replicated%transport_level_advances = sparse%transport_level_advances
      replicated%time = sparse%time
      replicated%steps = sparse%steps
      replicated%regrid_evaluations = sparse%regrid_evaluations
      replicated%regrids = sparse%regrids
      replicated%overlap_cells_transferred = sparse%overlap_cells_transferred
    end if
    call synchronize_owned_patch_tree_reactive_1d( &
      distribution, replicated, local_ok)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    ok = mpi_ok .and. accepted
  end subroutine gather_owned_patch_tree_reactive_1d

  subroutine migrate_owned_patch_tree_reactive_1d( &
      old_distribution, new_distribution, old_sparse, new_sparse, ok)
    type(mpi_amr_patch_distribution_1d), intent(in) :: old_distribution
    type(mpi_amr_patch_distribution_1d), intent(in) :: new_distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(in) :: old_sparse
    type(mpi_amr_sparse_reactive_solution_1d), intent(out) :: new_sparse
    logical, intent(out) :: ok

    logical :: local_ok, accepted, mpi_ok
    integer :: level, patch, nx, old_owner, new_owner

    ok = .false.
    local_ok = old_sparse%is_valid(old_distribution) .and. &
      old_distribution%rank == new_distribution%rank .and. &
      old_distribution%nranks == new_distribution%nranks .and. &
      mpi_amr_distribution_matches_hierarchy_1d( &
        new_distribution, old_sparse%hierarchy)
    call all_ranks_accept_sparse_1d( &
      old_distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
    call copy_sparse_metadata(old_sparse, new_distribution, new_sparse)
    do level = 1, old_sparse%hierarchy%level_count()
      do patch = 1, old_sparse%hierarchy%level_patch_count(level - 1)
        nx = old_distribution%levels(level)%cell_counts(patch)
        old_owner = old_distribution%owner_of(level - 1, patch)
        new_owner = new_distribution%owner_of(level - 1, patch)
        call migrate_one_patch_1d( &
          old_distribution%comm, old_distribution%rank, old_owner, new_owner, &
          old_sparse%nvar, nx, old_sparse%ghost_width, &
          old_sparse%levels(level)%patches(patch), &
          new_sparse%levels(level)%patches(patch), local_ok)
        if (.not. local_ok) return
      end do
    end do
    local_ok = new_sparse%is_valid(new_distribution)
    call all_ranks_accept_sparse_1d( &
      old_distribution, local_ok, accepted, mpi_ok)
    ok = mpi_ok .and. accepted
    if (.not. ok) new_sparse = mpi_amr_sparse_reactive_solution_1d()
  end subroutine migrate_owned_patch_tree_reactive_1d

  subroutine copy_sparse_metadata(source, distribution, target)
    type(mpi_amr_sparse_reactive_solution_1d), intent(in) :: source
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(out) :: target

    integer :: level, patch_count

    target%hierarchy = source%hierarchy
    target%rank = distribution%rank
    target%nranks = distribution%nranks
    target%nvar = source%nvar
    target%ghost_width = source%ghost_width
    target%level_advances = source%level_advances
    target%transport_level_advances = source%transport_level_advances
    target%time = source%time
    target%steps = source%steps
    target%regrid_evaluations = source%regrid_evaluations
    target%regrids = source%regrids
    target%overlap_cells_transferred = source%overlap_cells_transferred
    allocate(target%levels(source%hierarchy%level_count()))
    do level = 1, size(target%levels)
      patch_count = source%hierarchy%level_patch_count(level - 1)
      allocate(target%levels(level)%patches(patch_count))
      allocate(target%levels(level)%is_local(patch_count))
      target%levels(level)%is_local = distribution%levels(level)%owners == &
        distribution%rank
    end do
  end subroutine copy_sparse_metadata

  subroutine migrate_one_patch_1d( &
      comm, rank, old_owner, new_owner, nvar, nx, ghost_width, source, &
      destination, ok)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: rank, old_owner, new_owner, nvar, nx, ghost_width
    type(amr_patch_tree_reactive_patch_1d), intent(in) :: source
    type(amr_patch_tree_reactive_patch_1d), intent(inout) :: destination
    logical, intent(out) :: ok

    type(amr_patch_tree_reactive_patch_1d) :: work
    integer :: ierr

    call allocate_sparse_patch(work, nvar, nx, ghost_width)
    if (rank == old_owner) call copy_sparse_patch_values(source, work)
    call MPI_Bcast(work%state, size(work%state), MPI_DOUBLE_PRECISION, &
      old_owner, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast(work%temperature, size(work%temperature), &
      MPI_DOUBLE_PRECISION, old_owner, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast(work%left_ghost_state, size(work%left_ghost_state), &
      MPI_DOUBLE_PRECISION, old_owner, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast(work%right_ghost_state, size(work%right_ghost_state), &
      MPI_DOUBLE_PRECISION, old_owner, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast( &
      work%left_ghost_temperature, size(work%left_ghost_temperature), &
      MPI_DOUBLE_PRECISION, old_owner, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast( &
      work%right_ghost_temperature, size(work%right_ghost_temperature), &
      MPI_DOUBLE_PRECISION, old_owner, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    if (rank == new_owner) destination = work
    ok = .true.
  end subroutine migrate_one_patch_1d

  subroutine allocate_sparse_patch(patch, nvar, nx, ghost_width)
    type(amr_patch_tree_reactive_patch_1d), intent(out) :: patch
    integer, intent(in) :: nvar, nx, ghost_width

    allocate(patch%state(nvar, 0:nx + 1))
    allocate(patch%temperature(0:nx + 1))
    allocate(patch%left_ghost_state(nvar, ghost_width))
    allocate(patch%right_ghost_state(nvar, ghost_width))
    allocate(patch%left_ghost_temperature(ghost_width))
    allocate(patch%right_ghost_temperature(ghost_width))
    patch%state = 0.0_dp
    patch%temperature = 0.0_dp
    patch%left_ghost_state = 0.0_dp
    patch%right_ghost_state = 0.0_dp
    patch%left_ghost_temperature = 0.0_dp
    patch%right_ghost_temperature = 0.0_dp
  end subroutine allocate_sparse_patch

  subroutine copy_sparse_patch_values(source, destination)
    type(amr_patch_tree_reactive_patch_1d), intent(in) :: source
    type(amr_patch_tree_reactive_patch_1d), intent(inout) :: destination

    destination%state = source%state
    destination%temperature = source%temperature
    destination%left_ghost_state = source%left_ghost_state
    destination%right_ghost_state = source%right_ghost_state
    destination%left_ghost_temperature = source%left_ghost_temperature
    destination%right_ghost_temperature = source%right_ghost_temperature
  end subroutine copy_sparse_patch_values

  pure logical function sparse_patch_has_shape( &
      patch, nvar, nx, ghost_width) result(valid)
    type(amr_patch_tree_reactive_patch_1d), intent(in) :: patch
    integer, intent(in) :: nvar, nx, ghost_width

    valid = allocated(patch%state) .and. allocated(patch%temperature) .and. &
      allocated(patch%left_ghost_state) .and. &
      allocated(patch%right_ghost_state) .and. &
      allocated(patch%left_ghost_temperature) .and. &
      allocated(patch%right_ghost_temperature)
    if (.not. valid) return
    valid = size(patch%state, 1) == nvar .and. &
      lbound(patch%state, 2) == 0 .and. &
      ubound(patch%state, 2) == nx + 1 .and. &
      lbound(patch%temperature, 1) == 0 .and. &
      ubound(patch%temperature, 1) == nx + 1 .and. &
      size(patch%left_ghost_state, 1) == nvar .and. &
      size(patch%right_ghost_state, 1) == nvar .and. &
      size(patch%left_ghost_state, 2) == ghost_width .and. &
      size(patch%right_ghost_state, 2) == ghost_width .and. &
      size(patch%left_ghost_temperature) == ghost_width .and. &
      size(patch%right_ghost_temperature) == ghost_width
  end function sparse_patch_has_shape

  pure logical function sparse_patch_is_empty(patch) result(empty)
    type(amr_patch_tree_reactive_patch_1d), intent(in) :: patch

    empty = .not. allocated(patch%state) .and. &
      .not. allocated(patch%temperature) .and. &
      .not. allocated(patch%left_ghost_state) .and. &
      .not. allocated(patch%right_ghost_state) .and. &
      .not. allocated(patch%left_ghost_temperature) .and. &
      .not. allocated(patch%right_ghost_temperature)
  end function sparse_patch_is_empty

  pure integer function sparse_patch_value_count(patch) result(count)
    type(amr_patch_tree_reactive_patch_1d), intent(in) :: patch

    count = size(patch%state) + size(patch%temperature) + &
      size(patch%left_ghost_state) + size(patch%right_ghost_state) + &
      size(patch%left_ghost_temperature) + &
      size(patch%right_ghost_temperature)
  end function sparse_patch_value_count

  subroutine all_ranks_accept_sparse_1d( &
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
  end subroutine all_ranks_accept_sparse_1d

end module mpi_amr_sparse_patch_1d_mod
