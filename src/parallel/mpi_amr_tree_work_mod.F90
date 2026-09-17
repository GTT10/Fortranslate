module mpi_amr_tree_work_mod
  use mpi_f08, only: MPI_Wtime
  use precision_mod, only: dp
  implicit none
  private
  public :: begin_tree_work,end_tree_work,tree_work_start,tree_work_stop
  ! One process-local sampling interval; inactive for ordinary library calls.
  ! Not reentrant/thread-safe. The application owns interval lifetime and retries.
  real(dp), allocatable, save :: elapsed(:)
contains
  subroutine begin_tree_work(nodes,ok)
    integer, intent(in) :: nodes
    logical, intent(out) :: ok
    integer :: status
    if (allocated(elapsed)) deallocate(elapsed)
    ok=.false.
    if (nodes<1) return
    allocate(elapsed(nodes),stat=status)
    if (status/=0) return
    elapsed=0._dp;ok=.true.
  end subroutine
  subroutine end_tree_work(cost)
    real(dp), allocatable, intent(inout) :: cost(:)
    if (allocated(cost)) deallocate(cost)
    call move_alloc(elapsed,cost)
  end subroutine
  function tree_work_start() result(stamp)
    real(dp) :: stamp
    stamp=0._dp
    if (allocated(elapsed)) stamp=MPI_Wtime()
  end function
  subroutine tree_work_stop(node,stamp)
    integer, intent(in) :: node
    real(dp), intent(in) :: stamp
    if (.not.allocated(elapsed)) return
    elapsed(node)=elapsed(node)+max(0._dp,MPI_Wtime()-stamp)
  end subroutine
end module
