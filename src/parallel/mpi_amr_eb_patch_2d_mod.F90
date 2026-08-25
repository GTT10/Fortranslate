module mpi_amr_eb_patch_2d_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use reactive_2d_mod, only: advance_reactive_chemistry_2d
  use reactive_boundary_2d_mod, only: &
    reactive_boundary_set_2d, validate_reactive_boundary_set_2d
  use transport_database_mod, only: &
    gas_transport_species, compatible_transport_database
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_covered_cell
  use eb_reactive_reconstruction_2d_mod, only: &
    reactive_eb_exterior_state_2d
  use amr_eb_reactive_2d_mod, only: &
    build_reactive_eb_patch_exterior_2d, advance_reactive_eb_level_2d
  use eb_reactive_redistribution_2d_mod, only: &
    advance_reactive_eb_state_redistributed_2d
  use eb_reactive_transport_2d_mod, only: &
    reactive_eb_transport_fluxes_rhs_2d
  use amr_eb_flux_register_2d_mod, only: &
    amr_eb_flux_register_2d, initialize_amr_eb_flux_register_2d, &
    accumulate_coarse_eb_fluxes_2d, accumulate_fine_eb_fluxes_2d, &
    reflux_reactive_eb_state_patch_2d
  use amr_eb_regrid_2d_mod, only: &
    reactive_eb_patch_set_2d, average_down_reactive_eb_patch_set_2d, &
    composite_reactive_eb_patch_set_integral_2d
  use amr_eb_transport_2d_mod, only: recover_transport_temperature_2d
  use amr_eb_multilevel_reactive_2d_mod, only: &
    level_two_interface_is_regular
  use amr_eb_multipatch_transport_2d_mod, only: &
    close_cut_patch_set_conservation_2d
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

  public :: initialize_mpi_amr_eb_patch_distribution_2d
  public :: mpi_amr_eb_distribution_matches_patch_set_2d
  public :: synchronize_owned_reactive_eb_patch_set_2d
  public :: advance_owned_reactive_eb_patch_set_chemistry_2d
  public :: advance_owned_reactive_eb_patch_set_hydro_2d
  public :: advance_owned_reactive_eb_patch_set_transport_2d
  public :: advance_owned_reactive_eb_patch_set_strang_2d
  public :: scatter_owned_reactive_eb_patch_set_2d
  public :: materialize_owned_reactive_eb_patch_set_2d
  public :: average_down_sparse_owned_reactive_eb_patch_set_2d
  public :: advance_sparse_owned_reactive_eb_patch_set_chemistry_2d
  public :: advance_sparse_owned_reactive_eb_patch_set_strang_2d

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

  subroutine average_down_sparse_owned_reactive_eb_patch_set_2d( &
      species, distribution, sparse_patch_set, coarse_geometry, &
      patch_set_template, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(inout) :: sparse_patch_set
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set_template
    logical, intent(out) :: ok

    type(mpi_amr_eb_sparse_patch_set_2d) :: backup, candidate
    real(dp), allocatable :: primitive(:), restricted_state(:, :, :)
    real(dp) :: fine_volume, recovered_temperature, sound_speed
    logical :: accepted, entity_ok, global_ok, local_ok
    integer :: child, coarse_i, coarse_i_lower, coarse_i_upper, coarse_j
    integer :: coarse_j_lower, coarse_j_upper, component, fine_i_lower
    integer :: fine_i_upper, fine_j_lower, fine_j_upper, ierr, local_i
    integer :: local_j, nspecies, nspecies_maximum, nspecies_minimum
    integer :: owner, ratio, tile

    ok = .false.
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
      allocate(restricted_state( &
        candidate%nvar, coarse_i_upper - coarse_i_lower + 1, &
        coarse_j_upper - coarse_j_lower + 1))
      restricted_state = 0.0_dp
      owner = distribution%child_owner(child)
      ratio = &
        patch_set_template%children(child)%patch%refinement_ratio
      if (distribution%rank == owner) then
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
      call MPI_Bcast( &
        restricted_state, size(restricted_state), MPI_DOUBLE_PRECISION, &
        owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) then
        sparse_patch_set = backup
        return
      end if

      entity_ok = all(ieee_is_finite(restricted_state))
      do tile = 1, distribution%root_tile_count()
        if (.not. distribution%root_tile_is_local(tile)) cycle
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
      deallocate(restricted_state)
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
  end subroutine average_down_sparse_owned_reactive_eb_patch_set_2d

  subroutine advance_sparse_owned_reactive_eb_patch_set_chemistry_2d( &
      species, reactions, interval, relative_tolerance, absolute_tolerance, &
      distribution, sparse_patch_set, coarse_geometry, patch_set_template, &
      ok, local_entity_advances)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: interval, relative_tolerance, absolute_tolerance
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    type(mpi_amr_eb_sparse_patch_set_2d), intent(inout) :: sparse_patch_set
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set_template
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_entity_advances

    type(mpi_amr_eb_sparse_patch_set_2d) :: backup, candidate
    real(dp) :: controls(3), control_maximum(3), control_minimum(3)
    logical, allocatable :: active_mask(:, :)
    logical :: accepted, entity_ok, global_ok, local_ok
    integer :: advances, child, ierr, integer_controls(2)
    integer :: integer_maximum(2), integer_minimum(2)
    integer :: j_lower, j_upper, tile

    ok = .false.
    advances = 0
    if (present(local_entity_advances)) local_entity_advances = 0
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
      patch_set_template, local_ok)
    if (.not. local_ok) then
      sparse_patch_set = backup
      return
    end if
    sparse_patch_set = candidate
    ok = .true.
    if (present(local_entity_advances)) local_entity_advances = advances
  end subroutine advance_sparse_owned_reactive_eb_patch_set_chemistry_2d

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
    type(reactive_eb_patch_set_2d) :: materialized_set
    real(dp), allocatable :: fallback_state(:, :, :)
    real(dp), allocatable :: fallback_temperature(:, :)
    real(dp), allocatable :: materialized_state(:, :, :)
    real(dp), allocatable :: materialized_temperature(:, :)
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

    allocate(fallback_state( &
      candidate%nvar, coarse_geometry%nx, coarse_geometry%ny), source=0.0_dp)
    allocate(fallback_temperature( &
      coarse_geometry%nx, coarse_geometry%ny), source=1.0_dp)
    allocate(materialized_state, mold=fallback_state)
    allocate(materialized_temperature, mold=fallback_temperature)
    call materialize_owned_reactive_eb_patch_set_2d( &
      distribution, candidate, fallback_state, fallback_temperature, &
      coarse_geometry, patch_set_template, materialized_state, &
      materialized_temperature, materialized_set, local_ok)
    if (.not. local_ok) return

    call advance_owned_reactive_eb_patch_set_transport_2d( &
      species, transport, distribution, materialized_state, &
      materialized_temperature, coarse_geometry, materialized_set, &
      0.5_dp * dt, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      state_redist_max_order, local_ok, transport_one, theta_one, &
      selected_target)
    if (.not. local_ok) return
    call advance_owned_reactive_eb_patch_set_hydro_2d( &
      species, distribution, materialized_state, materialized_temperature, &
      coarse_geometry, materialized_set, solver, reconstruction, limiter, &
      state_redist_max_order, dt, local_ok, hydro_advances, selected_target)
    if (.not. local_ok) return
    call advance_owned_reactive_eb_patch_set_transport_2d( &
      species, transport, distribution, materialized_state, &
      materialized_temperature, coarse_geometry, materialized_set, &
      0.5_dp * dt, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      state_redist_max_order, local_ok, transport_two, theta_two, &
      selected_target)
    if (.not. local_ok) return

    call scatter_owned_reactive_eb_patch_set_2d( &
      distribution, size(species), materialized_state, &
      materialized_temperature, coarse_geometry, materialized_set, &
      candidate, local_ok)
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

  subroutine advance_owned_reactive_eb_patch_set_hydro_2d( &
      species, distribution, coarse_state, coarse_temperature, &
      coarse_geometry, patch_set, solver, reconstruction, limiter, &
      state_redist_max_order, dt, ok, local_level_advances, &
      state_redist_target_volume_fraction)
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
    integer :: nvar, owner, ratio, root_owner, substep
    integer :: string_codes(32, 3), string_maximum(32, 3)
    integer :: string_minimum(32, 3)

    ok = .false.
    advances = 0
    if (present(local_level_advances)) local_level_advances = 0
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
    root_owner = distribution%root_level_owner()
    entity_ok = root_owner >= 0 .and. root_owner < distribution%nranks
    if (distribution%rank == root_owner .and. entity_ok) then
      call advance_reactive_eb_level_2d( &
        species, root_start, root_start_temperature, coarse_geometry, &
        trim(solver), trim(reconstruction), trim(limiter), selected_target, &
        state_redist_max_order, dt, root_hydro, root_hydro_temperature, &
        coarse_x_flux, coarse_y_flux, entity_ok)
      if (entity_ok) advances = advances + 1
    end if
    call all_ranks_accept_eb_2d( &
      distribution, entity_ok, accepted, global_ok)
    if (.not. global_ok .or. .not. accepted) return
    call MPI_Bcast( &
      root_hydro, size(root_hydro), MPI_DOUBLE_PRECISION, root_owner, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Bcast( &
      root_hydro_temperature, size(root_hydro_temperature), &
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
    call collective_transport_preflight_2d( &
      species, transport, distribution, coarse_state, coarse_temperature, &
      coarse_geometry, patch_set, interval, viscosity_enabled, &
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
      species, transport, distribution, coarse_state, coarse_temperature, &
      coarse_geometry, patch_set, interval, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, state_redist_max_order, &
      target_volume_fraction, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(mpi_amr_eb_patch_distribution_2d), intent(in) :: distribution
    real(dp), intent(in) :: coarse_state(:, :, :)
    real(dp), intent(in) :: coarse_temperature(:, :)
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
      all(shape(coarse_state) == &
        [reactive_nvar(nsp), coarse_geometry%nx, coarse_geometry%ny]) .and. &
      all(shape(coarse_temperature) == &
        [coarse_geometry%nx, coarse_geometry%ny]) .and. &
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
