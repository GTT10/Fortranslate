module mpi_domain_1d_mod
  use iso_fortran_env, only: real64
  use mpi_f08
  implicit none
  private

  integer, parameter :: dp = real64

  type, public :: mpi_domain_1d
    type(MPI_Comm) :: comm = MPI_COMM_NULL
    integer :: rank = -1
    integer :: nranks = 0
    integer :: global_cells = 0
    integer :: local_cells = 0
    integer :: global_first = 0
    integer :: global_last = -1
    integer :: left_rank = MPI_PROC_NULL
    integer :: right_rank = MPI_PROC_NULL
    integer, allocatable :: counts(:)
    integer, allocatable :: displacements(:)
  end type mpi_domain_1d

  public :: initialize_mpi_domain_1d
  public :: exchange_periodic_halo_1d
  public :: global_minimum_1d
  public :: global_sum_1d
  public :: gather_state_1d

contains

  subroutine initialize_mpi_domain_1d(domain, global_cells, comm, ok)
    type(mpi_domain_1d), intent(out) :: domain
    integer, intent(in) :: global_cells
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer :: ierr, r, base_cells, remainder

    ok = .false.
    if (global_cells <= 0) return
    domain%comm = comm
    domain%global_cells = global_cells
    call MPI_Comm_rank(comm, domain%rank, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Comm_size(comm, domain%nranks, ierr)
    if (ierr /= MPI_SUCCESS .or. domain%nranks <= 0) return

    allocate(domain%counts(domain%nranks), domain%displacements(domain%nranks))
    base_cells = global_cells / domain%nranks
    remainder = modulo(global_cells, domain%nranks)
    domain%displacements(1) = 0
    do r = 0, domain%nranks - 1
      domain%counts(r + 1) = base_cells
      if (r < remainder) domain%counts(r + 1) = domain%counts(r + 1) + 1
      if (r > 0) then
        domain%displacements(r + 1) = domain%displacements(r) + domain%counts(r)
      end if
    end do
    domain%local_cells = domain%counts(domain%rank + 1)
    domain%global_first = domain%displacements(domain%rank + 1) + 1
    domain%global_last = domain%global_first + domain%local_cells - 1
    if (domain%nranks == 1) then
      domain%left_rank = 0
      domain%right_rank = 0
    else
      domain%left_rank = modulo(domain%rank - 1 + domain%nranks, domain%nranks)
      domain%right_rank = modulo(domain%rank + 1, domain%nranks)
    end if
    ok = domain%local_cells > 0 .and. sum(domain%counts) == global_cells
  end subroutine initialize_mpi_domain_1d

  subroutine exchange_periodic_halo_1d(domain, state, ok)
    type(mpi_domain_1d), intent(in) :: domain
    real(dp), intent(inout), contiguous :: state(:, 0:)
    logical, intent(out) :: ok
    type(MPI_Request) :: requests(4)
    type(MPI_Status) :: statuses(4)
    integer :: ierr, nvar

    ok = .false.
    nvar = size(state, 1)
    if (ubound(state, 2) < domain%local_cells + 1) return
    if (domain%nranks == 1) then
      state(:, 0) = state(:, domain%local_cells)
      state(:, domain%local_cells + 1) = state(:, 1)
      ok = .true.
      return
    end if

    call MPI_Irecv(state(:, 0), nvar, MPI_DOUBLE_PRECISION, &
      domain%left_rank, 100, domain%comm, requests(1), ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Irecv(state(:, domain%local_cells + 1), nvar, MPI_DOUBLE_PRECISION, &
      domain%right_rank, 101, domain%comm, requests(2), ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Isend(state(:, 1), nvar, MPI_DOUBLE_PRECISION, &
      domain%left_rank, 101, domain%comm, requests(3), ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Isend(state(:, domain%local_cells), nvar, MPI_DOUBLE_PRECISION, &
      domain%right_rank, 100, domain%comm, requests(4), ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Waitall(4, requests, statuses, ierr)
    ok = ierr == MPI_SUCCESS
  end subroutine exchange_periodic_halo_1d

  subroutine global_minimum_1d(domain, local_value, global_value, ok)
    type(mpi_domain_1d), intent(in) :: domain
    real(dp), intent(in) :: local_value
    real(dp), intent(out) :: global_value
    logical, intent(out) :: ok
    integer :: ierr
    call MPI_Allreduce(local_value, global_value, 1, MPI_DOUBLE_PRECISION, &
      MPI_MIN, domain%comm, ierr)
    ok = ierr == MPI_SUCCESS
  end subroutine global_minimum_1d

  subroutine global_sum_1d(domain, local_values, global_values, ok)
    type(mpi_domain_1d), intent(in) :: domain
    real(dp), intent(in), contiguous :: local_values(:)
    real(dp), intent(out), contiguous :: global_values(:)
    logical, intent(out) :: ok
    integer :: ierr
    ok = .false.
    if (size(global_values) /= size(local_values)) return
    call MPI_Allreduce(local_values, global_values, size(local_values), &
      MPI_DOUBLE_PRECISION, MPI_SUM, domain%comm, ierr)
    ok = ierr == MPI_SUCCESS
  end subroutine global_sum_1d

  subroutine gather_state_1d(domain, local_state, global_state, root, ok)
    type(mpi_domain_1d), intent(in) :: domain
    real(dp), intent(in), contiguous :: local_state(:, :)
    real(dp), allocatable, intent(out) :: global_state(:, :)
    integer, intent(in) :: root
    logical, intent(out) :: ok
    real(dp), allocatable :: send_buffer(:), receive_buffer(:)
    integer, allocatable :: element_counts(:), element_displacements(:)
    integer :: ierr, nvar, i, v, offset

    ok = .false.
    nvar = size(local_state, 1)
    if (size(local_state, 2) /= domain%local_cells) return
    allocate(send_buffer(nvar * domain%local_cells))
    offset = 0
    do i = 1, domain%local_cells
      do v = 1, nvar
        offset = offset + 1
        send_buffer(offset) = local_state(v, i)
      end do
    end do
    allocate(element_counts(domain%nranks), element_displacements(domain%nranks))
    element_counts = nvar * domain%counts
    element_displacements = nvar * domain%displacements
    if (domain%rank == root) then
      allocate(receive_buffer(nvar * domain%global_cells))
    else
      allocate(receive_buffer(0))
    end if
    call MPI_Gatherv(send_buffer, size(send_buffer), MPI_DOUBLE_PRECISION, &
      receive_buffer, element_counts, element_displacements, MPI_DOUBLE_PRECISION, &
      root, domain%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    if (domain%rank == root) then
      allocate(global_state(nvar, domain%global_cells))
      offset = 0
      do i = 1, domain%global_cells
        do v = 1, nvar
          offset = offset + 1
          global_state(v, i) = receive_buffer(offset)
        end do
      end do
    else
      allocate(global_state(0, 0))
    end if
    ok = .true.
  end subroutine gather_state_1d

end module mpi_domain_1d_mod
