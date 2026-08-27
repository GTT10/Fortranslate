module mpi_amr_eb_patch_tree_io_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use eb_geometry_2d_mod, only: eb_geometry_2d
  use amr_eb_patch_tree_2d_mod, only: &
    amr_eb_patch_tree_level_plan_2d, amr_eb_patch_tree_topology_2d, &
    initialize_amr_eb_patch_tree_topology_2d
  use amr_eb_patch_tree_reactive_2d_mod, only: &
    reactive_amr_eb_patch_tree_2d, &
    reactive_amr_eb_patch_tree_checkpoint_fingerprint_2d, &
    write_reactive_amr_eb_patch_tree_2d_checkpoint, &
    read_reactive_amr_eb_patch_tree_2d_checkpoint, &
    write_reactive_amr_eb_patch_tree_2d_csv, &
    composite_integral_reactive_amr_eb_patch_tree_2d
  use mpi_amr_eb_patch_tree_2d_mod, only: &
    mpi_amr_eb_patch_tree_distribution_2d, &
    mpi_sparse_reactive_amr_eb_patch_tree_2d, &
    initialize_mpi_amr_eb_patch_tree_distribution_2d, &
    gather_sparse_owned_reactive_amr_eb_patch_tree_to_root_2d, &
    scatter_root_reactive_amr_eb_patch_tree_to_sparse_2d
  implicit none
  private

  integer, parameter :: maximum_checkpoint_levels = 64
  integer, parameter :: maximum_checkpoint_patches = 1000000
  integer, parameter :: maximum_checkpoint_geometry_cells = 100000000

  public :: write_sparse_owned_reactive_amr_eb_patch_tree_2d_checkpoint
  public :: read_sparse_owned_reactive_amr_eb_patch_tree_2d_checkpoint
  public :: write_sparse_owned_reactive_amr_eb_patch_tree_2d_csv

contains

  subroutine write_sparse_owned_reactive_amr_eb_patch_tree_2d_csv( &
      path, species, distribution, sparse, root, time, ok, &
      local_entity_transfers)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(in) :: sparse
    integer, intent(in) :: root
    real(dp), intent(in) :: time
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_entity_transfers

    type(reactive_amr_eb_patch_tree_2d) :: gathered
    logical :: controls_ok, gathered_ok, write_ok
    integer :: ierr, transfers

    ok = .false.
    transfers = 0
    if (present(local_entity_transfers)) local_entity_transfers = 0
    call output_write_controls_match_2d( &
      distribution%comm, root, time, controls_ok)
    if (.not. controls_ok) return
    call checkpoint_species_match_2d( &
      distribution%comm, distribution%rank, root, species, controls_ok)
    if (.not. controls_ok) return

    call gather_sparse_owned_reactive_amr_eb_patch_tree_to_root_2d( &
      distribution, sparse, root, gathered, gathered_ok, transfers)
    if (.not. gathered_ok) return
    write_ok = .true.
    if (distribution%rank == root) call &
      write_reactive_amr_eb_patch_tree_2d_csv( &
        path, species, gathered, time, write_ok)
    call MPI_Bcast( &
      write_ok, 1, MPI_LOGICAL, root, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. .not. write_ok) return

    ok = .true.
    if (present(local_entity_transfers)) local_entity_transfers = transfers
  end subroutine write_sparse_owned_reactive_amr_eb_patch_tree_2d_csv

  subroutine write_sparse_owned_reactive_amr_eb_patch_tree_2d_checkpoint( &
      path, species, distribution, sparse, root, time, steps, regrids, &
      minimum_dt, ok, local_entity_transfers, fingerprint, &
      minimum_transport_theta, initial_integrals)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(in) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(in) :: sparse
    integer, intent(in) :: root
    real(dp), intent(in) :: time, minimum_dt
    integer, intent(in) :: steps, regrids
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_entity_transfers
    type(reactive_amr_eb_patch_tree_checkpoint_fingerprint_2d), &
      intent(in), optional :: fingerprint
    real(dp), intent(in), optional :: minimum_transport_theta
    real(dp), intent(in), optional :: initial_integrals(:)

    type(reactive_amr_eb_patch_tree_2d) :: gathered
    logical :: controls_ok, gathered_ok, write_ok
    integer :: ierr, transfers
    real(dp), allocatable :: selected_initial_integrals(:)
    real(dp) :: selected_minimum_transport_theta

    ok = .false.
    transfers = 0
    selected_minimum_transport_theta = 1.0_dp
    if (present(minimum_transport_theta)) &
      selected_minimum_transport_theta = minimum_transport_theta
    if (present(local_entity_transfers)) local_entity_transfers = 0
    call checkpoint_write_controls_match_2d( &
      distribution%comm, root, time, minimum_dt, &
      selected_minimum_transport_theta, steps, regrids, sparse%nvar, &
      controls_ok, initial_integrals)
    if (.not. controls_ok) return
    call checkpoint_species_match_2d( &
      distribution%comm, distribution%rank, root, species, controls_ok)
    if (.not. controls_ok) return

    call gather_sparse_owned_reactive_amr_eb_patch_tree_to_root_2d( &
      distribution, sparse, root, gathered, gathered_ok, transfers)
    if (.not. gathered_ok) return
    write_ok = .true.
    if (distribution%rank == root) then
      allocate(selected_initial_integrals(gathered%nvar))
      if (present(initial_integrals)) then
        selected_initial_integrals = initial_integrals
      else
        call composite_integral_reactive_amr_eb_patch_tree_2d( &
          gathered, selected_initial_integrals, write_ok)
      end if
      if (write_ok .and. present(fingerprint)) then
        call write_reactive_amr_eb_patch_tree_2d_checkpoint( &
          path, species, gathered, time, steps, regrids, minimum_dt, &
          write_ok, fingerprint, selected_minimum_transport_theta, &
          selected_initial_integrals)
      else if (write_ok) then
        call write_reactive_amr_eb_patch_tree_2d_checkpoint( &
          path, species, gathered, time, steps, regrids, minimum_dt, write_ok, &
          minimum_transport_theta=selected_minimum_transport_theta, &
          initial_integrals=selected_initial_integrals)
      end if
    end if
    call MPI_Bcast( &
      write_ok, 1, MPI_LOGICAL, root, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. .not. write_ok) return

    ok = .true.
    if (present(local_entity_transfers)) local_entity_transfers = transfers
  end subroutine write_sparse_owned_reactive_amr_eb_patch_tree_2d_checkpoint

  subroutine read_sparse_owned_reactive_amr_eb_patch_tree_2d_checkpoint( &
      path, species, comm, root, maximum_levels, subcycle_exponent, &
      distribution, sparse, time, steps, regrids, minimum_dt, ok, &
      local_entity_transfers, fingerprint, minimum_transport_theta, &
      initial_integrals)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: root, maximum_levels, subcycle_exponent
    type(mpi_amr_eb_patch_tree_distribution_2d), intent(out) :: distribution
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(out) :: sparse
    real(dp), intent(out) :: time, minimum_dt
    integer, intent(out) :: steps, regrids
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_entity_transfers
    type(reactive_amr_eb_patch_tree_checkpoint_fingerprint_2d), &
      intent(in), optional :: fingerprint
    real(dp), intent(out), optional :: minimum_transport_theta
    real(dp), allocatable, intent(out), optional :: initial_integrals(:)

    type(reactive_amr_eb_patch_tree_2d) :: loaded
    type(amr_eb_patch_tree_topology_2d) :: topology
    real(dp), allocatable :: restored_initial_integrals(:)
    real(dp) :: real_metadata(3)
    integer :: ierr, integer_metadata(2), rank, transfers
    logical :: controls_ok, distributed_ok, read_ok, topology_ok

    distribution = mpi_amr_eb_patch_tree_distribution_2d()
    sparse = mpi_sparse_reactive_amr_eb_patch_tree_2d()
    time = 0.0_dp
    minimum_dt = 0.0_dp
    if (present(minimum_transport_theta)) minimum_transport_theta = 1.0_dp
    steps = 0
    regrids = 0
    ok = .false.
    transfers = 0
    real_metadata = 0.0_dp
    real_metadata(3) = 1.0_dp
    integer_metadata = 0
    if (present(local_entity_transfers)) local_entity_transfers = 0

    call checkpoint_read_controls_match_2d( &
      comm, root, maximum_levels, subcycle_exponent, rank, controls_ok)
    if (.not. controls_ok) return
    call checkpoint_species_match_2d( &
      comm, rank, root, species, controls_ok)
    if (.not. controls_ok) return

    read_ok = .true.
    if (rank == root) then
      if (present(fingerprint)) then
        call read_reactive_amr_eb_patch_tree_2d_checkpoint( &
          path, species, maximum_levels, loaded, real_metadata(1), &
          integer_metadata(1), integer_metadata(2), real_metadata(2), &
          read_ok, fingerprint, real_metadata(3), restored_initial_integrals)
      else
        call read_reactive_amr_eb_patch_tree_2d_checkpoint( &
          path, species, maximum_levels, loaded, real_metadata(1), &
          integer_metadata(1), integer_metadata(2), real_metadata(2), read_ok, &
          minimum_transport_theta=real_metadata(3), &
          initial_integrals=restored_initial_integrals)
      end if
    end if
    call MPI_Bcast(read_ok, 1, MPI_LOGICAL, root, comm, ierr)
    if (ierr /= MPI_SUCCESS .or. .not. read_ok) return
    if (rank /= root) allocate(restored_initial_integrals( &
      reactive_nvar(size(species))))
    call MPI_Bcast( &
      restored_initial_integrals, size(restored_initial_integrals), &
      MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Bcast( &
      real_metadata, size(real_metadata), MPI_DOUBLE_PRECISION, root, &
      comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Bcast( &
      integer_metadata, size(integer_metadata), MPI_INTEGER, root, comm, &
      ierr)
    if (ierr /= MPI_SUCCESS) return

    call broadcast_root_patch_tree_topology_2d( &
      comm, rank, root, loaded, topology, topology_ok)
    if (.not. topology_ok) return
    call initialize_mpi_amr_eb_patch_tree_distribution_2d( &
      topology, comm, distribution, distributed_ok, subcycle_exponent)
    if (.not. distributed_ok) then
      distribution = mpi_amr_eb_patch_tree_distribution_2d()
      return
    end if
    call scatter_root_reactive_amr_eb_patch_tree_to_sparse_2d( &
      distribution, topology, loaded, root, sparse, distributed_ok, &
      transfers)
    if (.not. distributed_ok) then
      distribution = mpi_amr_eb_patch_tree_distribution_2d()
      sparse = mpi_sparse_reactive_amr_eb_patch_tree_2d()
      return
    end if

    time = real_metadata(1)
    minimum_dt = real_metadata(2)
    if (present(minimum_transport_theta)) &
      minimum_transport_theta = real_metadata(3)
    if (present(initial_integrals)) &
      allocate(initial_integrals, source=restored_initial_integrals)
    steps = integer_metadata(1)
    regrids = integer_metadata(2)
    ok = .true.
    if (present(local_entity_transfers)) local_entity_transfers = transfers
  end subroutine read_sparse_owned_reactive_amr_eb_patch_tree_2d_checkpoint

  subroutine checkpoint_write_controls_match_2d( &
      comm, root, time, minimum_dt, minimum_transport_theta, &
      steps, regrids, nvar, ok, initial_integrals)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: root, steps, regrids, nvar
    real(dp), intent(in) :: time, minimum_dt, minimum_transport_theta
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: initial_integrals(:)

    real(dp) :: real_values(3), real_minimum(3), real_maximum(3)
    real(dp), allocatable :: integral_minimum(:), integral_maximum(:)
    integer :: ierr, integer_values(5), integer_minimum(5), nranks
    integer :: integer_maximum(5), integral_presence
    logical :: accepted, local_ok

    real_values = [time, minimum_dt, minimum_transport_theta]
    integral_presence = merge(1, 0, present(initial_integrals))
    integer_values = [root, steps, regrids, integral_presence, nvar]
    call MPI_Comm_size(comm, nranks, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    local_ok = root >= 0 .and. root < nranks .and. &
      ieee_is_finite(time) .and. ieee_is_finite(minimum_dt) .and. &
      ieee_is_finite(minimum_transport_theta) .and. &
      minimum_transport_theta >= 0.0_dp .and. &
      minimum_transport_theta <= 1.0_dp
    call MPI_Allreduce( &
      local_ok, accepted, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
    if (ierr /= MPI_SUCCESS .or. .not. accepted) then
      ok = .false.
      return
    end if
    call MPI_Allreduce( &
      real_values, real_minimum, size(real_values), MPI_DOUBLE_PRECISION, &
      MPI_MIN, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Allreduce( &
      real_values, real_maximum, size(real_values), MPI_DOUBLE_PRECISION, &
      MPI_MAX, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Allreduce( &
      integer_values, integer_minimum, size(integer_values), MPI_INTEGER, &
      MPI_MIN, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Allreduce( &
      integer_values, integer_maximum, size(integer_values), MPI_INTEGER, &
      MPI_MAX, comm, ierr)
    ok = ierr == MPI_SUCCESS .and. &
      all(real_minimum == real_maximum) .and. &
      all(integer_minimum == integer_maximum)
    if (.not. ok .or. integral_presence == 0) return
    local_ok = size(initial_integrals) == nvar
    if (local_ok) local_ok = all(ieee_is_finite(initial_integrals))
    call MPI_Allreduce( &
      local_ok, accepted, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
    if (ierr /= MPI_SUCCESS .or. .not. accepted) then
      ok = .false.
      return
    end if
    allocate(integral_minimum(nvar), integral_maximum(nvar))
    call MPI_Allreduce( &
      initial_integrals, integral_minimum, nvar, MPI_DOUBLE_PRECISION, &
      MPI_MIN, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Allreduce( &
      initial_integrals, integral_maximum, nvar, MPI_DOUBLE_PRECISION, &
      MPI_MAX, comm, ierr)
    ok = ierr == MPI_SUCCESS .and. &
      all(integral_minimum == integral_maximum)
  end subroutine checkpoint_write_controls_match_2d

  subroutine output_write_controls_match_2d(comm, root, time, ok)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: root
    real(dp), intent(in) :: time
    logical, intent(out) :: ok

    real(dp) :: maximum_time, minimum_time
    integer :: ierr, maximum_root, minimum_root, nranks
    logical :: accepted, local_ok

    call MPI_Comm_size(comm, nranks, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    local_ok = root >= 0 .and. root < nranks .and. &
      ieee_is_finite(time) .and. time >= 0.0_dp
    call MPI_Allreduce( &
      local_ok, accepted, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
    if (ierr /= MPI_SUCCESS .or. .not. accepted) then
      ok = .false.
      return
    end if
    call MPI_Allreduce( &
      time, minimum_time, 1, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Allreduce( &
      time, maximum_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Allreduce(root, minimum_root, 1, MPI_INTEGER, MPI_MIN, &
      comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Allreduce(root, maximum_root, 1, MPI_INTEGER, MPI_MAX, &
      comm, ierr)
    ok = ierr == MPI_SUCCESS .and. minimum_time == maximum_time .and. &
      minimum_root == maximum_root
  end subroutine output_write_controls_match_2d

  subroutine checkpoint_read_controls_match_2d( &
      comm, root, maximum_levels, subcycle_exponent, rank, ok)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: root, maximum_levels, subcycle_exponent
    integer, intent(out) :: rank
    logical, intent(out) :: ok

    integer :: ierr, nranks, values(3), minimum_values(3), maximum_values(3)
    logical :: accepted, local_ok

    rank = -1
    call MPI_Comm_rank(comm, rank, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Comm_size(comm, nranks, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    local_ok = root >= 0 .and. root < nranks .and. maximum_levels >= 1 &
      .and. maximum_levels <= maximum_checkpoint_levels .and. &
      subcycle_exponent >= 0 .and. subcycle_exponent <= 2
    call MPI_Allreduce( &
      local_ok, accepted, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
    if (ierr /= MPI_SUCCESS .or. .not. accepted) then
      ok = .false.
      return
    end if
    values = [root, maximum_levels, subcycle_exponent]
    call MPI_Allreduce( &
      values, minimum_values, size(values), MPI_INTEGER, MPI_MIN, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Allreduce( &
      values, maximum_values, size(values), MPI_INTEGER, MPI_MAX, comm, ierr)
    ok = ierr == MPI_SUCCESS .and. all(minimum_values == maximum_values)
  end subroutine checkpoint_read_controls_match_2d

  subroutine checkpoint_species_match_2d( &
      comm, rank, root, species, ok)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: rank, root
    type(nasa7_species), intent(in) :: species(:)
    logical, intent(out) :: ok

    character(len=128) :: root_name
    integer :: count_maximum, count_minimum, ierr, local_count, species_index
    logical :: accepted, local_ok

    local_count = size(species)
    call MPI_Allreduce( &
      local_count, count_minimum, 1, MPI_INTEGER, MPI_MIN, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Allreduce( &
      local_count, count_maximum, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
    if (ierr /= MPI_SUCCESS .or. count_minimum < 1 .or. &
        count_minimum /= count_maximum) then
      ok = .false.
      return
    end if
    local_ok = .true.
    do species_index = 1, count_minimum
      root_name = ""
      if (rank == root) root_name = trim(species(species_index)%name)
      call MPI_Bcast( &
        root_name, len(root_name), MPI_CHARACTER, root, comm, ierr)
      if (ierr /= MPI_SUCCESS) then
        ok = .false.
        return
      end if
      local_ok = local_ok .and. &
        trim(root_name) == trim(species(species_index)%name)
    end do
    call MPI_Allreduce( &
      local_ok, accepted, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
    ok = ierr == MPI_SUCCESS .and. accepted
  end subroutine checkpoint_species_match_2d

  subroutine broadcast_root_patch_tree_topology_2d( &
      comm, rank, root, root_tree, topology, ok)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: rank, root
    type(reactive_amr_eb_patch_tree_2d), intent(in) :: root_tree
    type(amr_eb_patch_tree_topology_2d), intent(out) :: topology
    logical, intent(out) :: ok

    type(amr_eb_patch_tree_level_plan_2d), allocatable :: plans(:)
    type(eb_geometry_2d) :: root_geometry, source_geometry
    integer :: child, header(2), ierr, level_count, metadata(5), relation
    logical :: accepted, geometry_ok, local_ok

    topology = amr_eb_patch_tree_topology_2d()
    level_count = 0
    if (rank == root) level_count = root_tree%level_count()
    call MPI_Bcast(level_count, 1, MPI_INTEGER, root, comm, ierr)
    if (ierr /= MPI_SUCCESS .or. level_count < 1 .or. &
        level_count > maximum_checkpoint_levels) then
      ok = .false.
      return
    end if
    source_geometry = eb_geometry_2d()
    if (rank == root) source_geometry = root_tree%topology%root_geometry
    call broadcast_root_geometry_2d( &
      comm, rank, root, source_geometry, root_geometry, geometry_ok)
    if (.not. geometry_ok) then
      ok = .false.
      return
    end if

    allocate(plans(level_count - 1))
    do relation = 1, size(plans)
      header = 0
      if (rank == root) header = [ &
        root_tree%topology%relations(relation)%refinement_ratio, &
        root_tree%topology%relations(relation)%child_patch_count()]
      call MPI_Bcast(header, size(header), MPI_INTEGER, root, comm, ierr)
      if (ierr /= MPI_SUCCESS .or. header(1) < 2 .or. header(2) < 1 .or. &
          header(2) > maximum_checkpoint_patches) then
        ok = .false.
        return
      end if
      plans(relation)%refinement_ratio = header(1)
      allocate(plans(relation)%children(header(2)))
      do child = 1, header(2)
        metadata = 0
        if (rank == root) metadata = [ &
          root_tree%topology%relations(relation)%children(child)% &
            parent_patch, &
          root_tree%topology%relations(relation)%children(child)%patch% &
            coarse_i_lower, &
          root_tree%topology%relations(relation)%children(child)%patch% &
            coarse_i_upper, &
          root_tree%topology%relations(relation)%children(child)%patch% &
            coarse_j_lower, &
          root_tree%topology%relations(relation)%children(child)%patch% &
            coarse_j_upper]
        call MPI_Bcast( &
          metadata, size(metadata), MPI_INTEGER, root, comm, ierr)
        if (ierr /= MPI_SUCCESS) then
          ok = .false.
          return
        end if
        plans(relation)%children(child)%parent_patch = metadata(1)
        plans(relation)%children(child)%coarse_i_lower = metadata(2)
        plans(relation)%children(child)%coarse_i_upper = metadata(3)
        plans(relation)%children(child)%coarse_j_lower = metadata(4)
        plans(relation)%children(child)%coarse_j_upper = metadata(5)
        source_geometry = eb_geometry_2d()
        if (rank == root) source_geometry = root_tree%topology% &
          relations(relation)%children(child)%geometry
        call broadcast_root_geometry_2d( &
          comm, rank, root, source_geometry, &
          plans(relation)%children(child)%geometry, geometry_ok)
        if (.not. geometry_ok) then
          ok = .false.
          return
        end if
      end do
    end do

    call initialize_amr_eb_patch_tree_topology_2d( &
      root_geometry, plans, topology, local_ok)
    call MPI_Allreduce( &
      local_ok, accepted, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
    ok = ierr == MPI_SUCCESS .and. accepted
    if (.not. ok) topology = amr_eb_patch_tree_topology_2d()
  end subroutine broadcast_root_patch_tree_topology_2d

  subroutine broadcast_root_geometry_2d( &
      comm, rank, root, root_geometry, geometry, ok)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: rank, root
    type(eb_geometry_2d), intent(in) :: root_geometry
    type(eb_geometry_2d), intent(out) :: geometry
    logical, intent(out) :: ok

    real(dp) :: metadata(6)
    integer :: header(2), ierr, nx, ny

    geometry = eb_geometry_2d()
    header = 0
    metadata = 0.0_dp
    if (rank == root) then
      header = [root_geometry%nx, root_geometry%ny]
      metadata = [root_geometry%x_lower, root_geometry%x_upper, &
        root_geometry%y_lower, root_geometry%y_upper, &
        root_geometry%dx, root_geometry%dy]
    end if
    call MPI_Bcast(header, size(header), MPI_INTEGER, root, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    nx = header(1)
    ny = header(2)
    if (nx < 1 .or. ny < 1) then
      ok = .false.
      return
    end if
    if (nx > maximum_checkpoint_geometry_cells / ny) then
      ok = .false.
      return
    end if
    call MPI_Bcast( &
      metadata, size(metadata), MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr /= MPI_SUCCESS .or. any(.not. ieee_is_finite(metadata))) then
      ok = .false.
      return
    end if
    if (rank == root) then
      geometry = root_geometry
    else
      call allocate_checkpoint_geometry_2d(nx, ny, metadata, geometry)
    end if

    call MPI_Bcast(geometry%volume_fraction, size( &
      geometry%volume_fraction), MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast(geometry%cell_centroid_x, size( &
      geometry%cell_centroid_x), MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast(geometry%cell_centroid_y, size( &
      geometry%cell_centroid_y), MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast(geometry%x_face_fraction, size( &
      geometry%x_face_fraction), MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast(geometry%y_face_fraction, size( &
      geometry%y_face_fraction), MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast(geometry%x_face_centroid_y, size( &
      geometry%x_face_centroid_y), MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast(geometry%y_face_centroid_x, size( &
      geometry%y_face_centroid_x), MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast(geometry%boundary_length, size( &
      geometry%boundary_length), MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast(geometry%boundary_centroid_x, size( &
      geometry%boundary_centroid_x), MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast(geometry%boundary_centroid_y, size( &
      geometry%boundary_centroid_y), MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast(geometry%boundary_normal_x, size( &
      geometry%boundary_normal_x), MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast(geometry%boundary_normal_y, size( &
      geometry%boundary_normal_y), MPI_DOUBLE_PRECISION, root, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast(geometry%boundary_normal_integral_x, size( &
      geometry%boundary_normal_integral_x), MPI_DOUBLE_PRECISION, root, &
      comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast(geometry%boundary_normal_integral_y, size( &
      geometry%boundary_normal_integral_y), MPI_DOUBLE_PRECISION, root, &
      comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast( &
      geometry%cell_type, size(geometry%cell_type), MPI_INTEGER, root, &
      comm, ierr)
    ok = ierr == MPI_SUCCESS .and. geometry%is_valid()
  end subroutine broadcast_root_geometry_2d

  subroutine allocate_checkpoint_geometry_2d( &
      nx, ny, metadata, geometry)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: metadata(6)
    type(eb_geometry_2d), intent(out) :: geometry

    geometry%nx = nx
    geometry%ny = ny
    geometry%x_lower = metadata(1)
    geometry%x_upper = metadata(2)
    geometry%y_lower = metadata(3)
    geometry%y_upper = metadata(4)
    geometry%dx = metadata(5)
    geometry%dy = metadata(6)
    allocate(geometry%volume_fraction(1:nx, 1:ny))
    allocate(geometry%cell_centroid_x(1:nx, 1:ny))
    allocate(geometry%cell_centroid_y(1:nx, 1:ny))
    allocate(geometry%x_face_fraction(0:nx, 1:ny))
    allocate(geometry%y_face_fraction(1:nx, 0:ny))
    allocate(geometry%x_face_centroid_y(0:nx, 1:ny))
    allocate(geometry%y_face_centroid_x(1:nx, 0:ny))
    allocate(geometry%boundary_length(1:nx, 1:ny))
    allocate(geometry%boundary_centroid_x(1:nx, 1:ny))
    allocate(geometry%boundary_centroid_y(1:nx, 1:ny))
    allocate(geometry%boundary_normal_x(1:nx, 1:ny))
    allocate(geometry%boundary_normal_y(1:nx, 1:ny))
    allocate(geometry%boundary_normal_integral_x(1:nx, 1:ny))
    allocate(geometry%boundary_normal_integral_y(1:nx, 1:ny))
    allocate(geometry%cell_type(1:nx, 1:ny))
  end subroutine allocate_checkpoint_geometry_2d

end module mpi_amr_eb_patch_tree_io_2d_mod
