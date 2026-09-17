module mpi_amr_patch_tree_restriction_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: restrict_average_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  implicit none
  private
  public :: average_down_mpi_amr_tree_3d
contains
  subroutine average_down_mpi_amr_tree_3d(tree,owners,local,nvar,comm,ok,minimum_depth)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar
    type(amr_tree_state_3d), intent(inout) :: local
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer, optional, intent(in) :: minimum_depth
    type(amr_tree_state_3d) :: candidate
    real(dp), allocatable :: restricted(:,:,:,:)
    integer :: rank,status,ios,node,parent,lo(3),hi(3),n(3),axis,count,depth,low,high,ratio
    integer(int64) :: elements
    logical :: valid,ready
    ok=.false.
    call validate_mpi_amr_tree_local_3d(tree,owners,local,nvar,comm,ready)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    depth=0
    if (present(minimum_depth)) depth=minimum_depth
    call MPI_Allreduce(depth,low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(depth,high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (low/=high.or.depth<0.or.depth>=tree%level_count()) return
    valid=.true.
    do node=1,size(tree%nodes)
      if (owners(node)==rank) then
        if (.not.all(ieee_is_finite(local%fields(node)%values))) valid=.false.
      end if
      if (node==1) cycle
      ratio=tree%nodes(node)%patch%refinement_ratio
      ! The shared restriction kernel forms ratio**3 using default integers.
      if (ratio>huge(1)/ratio/ratio) valid=.false.
      n=tree%nodes(node)%shape/ratio
      elements=int(nvar,int64)
      do axis=1,3
        if (elements>int(huge(1),int64)/int(n(axis),int64)) then
          valid=.false.;exit
        end if
        elements=elements*int(n(axis),int64)
      end do
    end do
    call agree(valid,ready)
    if (.not.ready) return
    allocate(candidate%fields(size(tree%nodes)),stat=ios)
    valid=ios==0
    if (valid) then
      do node=1,size(tree%nodes)
        if (owners(node)/=rank) cycle
        allocate(candidate%fields(node)%values,source=local%fields(node)%values,stat=ios)
        if (ios/=0) then
          valid=.false.;exit
        end if
      end do
    end if
    call agree(valid,ready)
    if (.not.ready) return
    ! Parent-before-child indexing makes reverse traversal bottom-up, including
    ! branches whose child and parent owners differ at successive levels.
    do node=size(tree%nodes),2,-1
      parent=tree%nodes(node)%parent
      if (tree%nodes(parent)%level<depth) cycle
      associate(box=>tree%nodes(node)%patch)
        lo=[box%coarse_i_lower,box%coarse_j_lower,box%coarse_k_lower]
        hi=[box%coarse_i_upper,box%coarse_j_upper,box%coarse_k_upper]
      end associate
      n=hi-lo+1
      valid=.true.
      if (rank==owners(node).or.rank==owners(parent)) then
        allocate(restricted(nvar,n(1),n(2),n(3)),stat=ios)
        valid=ios==0
      end if
      call agree(valid,ready)
      if (.not.ready) return
      if (rank==owners(node)) then
        call restrict_average_3d(candidate%fields(node)%values,tree%nodes(node)%patch,restricted,valid)
      end if
      call agree(valid,ready)
      if (.not.ready) return
      status=MPI_SUCCESS
      count=int(int(nvar,int64)*product(int(n,int64)))
      if (owners(node)/=owners(parent)) then
        if (rank==owners(node)) then
          call MPI_Send(restricted,count,MPI_DOUBLE_PRECISION,owners(parent),28403,comm,status)
        else if (rank==owners(parent)) then
          call MPI_Recv(restricted,count,MPI_DOUBLE_PRECISION,owners(node),28403,comm,MPI_STATUS_IGNORE,status)
        end if
      end if
      call agree(status==MPI_SUCCESS,ready)
      if (.not.ready) return
      if (rank==owners(parent)) then
        candidate%fields(parent)%values(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=restricted
      end if
      if (allocated(restricted)) deallocate(restricted)
    end do
    call move_alloc(candidate%fields,local%fields)
    ok=.true.
  contains
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) result=.false.
    end subroutine
  end subroutine
end module
