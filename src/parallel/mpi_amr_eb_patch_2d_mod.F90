module mpi_amr_eb_patch_2d_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use reactive_1d_mod, only: reactive_nvar
  use reactive_2d_mod, only: advance_reactive_chemistry_2d
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_covered_cell
  use amr_eb_regrid_2d_mod, only: &
    reactive_eb_patch_set_2d, average_down_reactive_eb_patch_set_2d
  implicit none
  private

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
    procedure :: root_tile_is_local => mpi_amr_eb_root_tile_is_local
    procedure :: child_is_local => mpi_amr_eb_child_is_local
    procedure :: is_valid => mpi_amr_eb_distribution_is_valid
  end type mpi_amr_eb_patch_distribution_2d

  public :: initialize_mpi_amr_eb_patch_distribution_2d
  public :: mpi_amr_eb_distribution_matches_patch_set_2d
  public :: synchronize_owned_reactive_eb_patch_set_2d
  public :: advance_owned_reactive_eb_patch_set_chemistry_2d

contains

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
