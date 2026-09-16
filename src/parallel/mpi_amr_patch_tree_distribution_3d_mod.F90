module mpi_amr_patch_tree_distribution_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  implicit none
  private
  public :: scatter_mpi_amr_tree_3d,gather_mpi_amr_tree_3d
contains
  subroutine scatter_mpi_amr_tree_3d(tree,owners,global,local,comm,root,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),root
    type(amr_tree_state_3d), intent(inout) :: global,local
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    call transfer(tree,owners,global,local,comm,root,.true.,ok)
  end subroutine
  subroutine gather_mpi_amr_tree_3d(tree,owners,local,global,comm,root,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),root
    type(amr_tree_state_3d), intent(inout) :: local,global
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    call transfer(tree,owners,global,local,comm,root,.false.,ok)
  end subroutine

  subroutine transfer(tree,owners,global,local,comm,root,scatter,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),root
    type(amr_tree_state_3d), intent(inout) :: global,local
    type(MPI_Comm), intent(in) :: comm
    logical, intent(in) :: scatter
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: candidate
    integer, allocatable :: records(:,:),reference(:,:)
    integer :: rank,nranks,ierr,node,nvar,local_nvar,n(3),count,ios,header(3),low(3),high(3),axis
    integer(int64) :: elements
    logical :: valid,ready
    ok=.false.
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=tree%is_valid().and.root>=0.and.root<nranks
    call agree(valid,ready)
    if (.not.ready) return
    header=[size(tree%nodes),root,merge(1,0,scatter)]
    call MPI_Allreduce(header,low,3,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(header,high,3,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (any(low/=high)) return
    valid=size(owners)==size(tree%nodes)
    if (valid) valid=all(owners>=0).and.all(owners<nranks)
    if (size(tree%nodes)>huge(1)/13) valid=.false.
    call agree(valid,ready)
    if (.not.ready) return
    allocate(records(13,size(tree%nodes)),reference(13,size(tree%nodes)),stat=ios)
    call agree(ios==0,ready)
    if (.not.ready) return
    records=0
    do node=1,size(tree%nodes)
      records(1:5,node)=[tree%nodes(node)%parent,tree%nodes(node)%level,tree%nodes(node)%shape]
      if (node>1) then
        associate(box=>tree%nodes(node)%patch)
          records(6:12,node)=[box%coarse_i_lower,box%coarse_i_upper,box%coarse_j_lower,box%coarse_j_upper, &
            box%coarse_k_lower,box%coarse_k_upper,box%refinement_ratio]
        end associate
      end if
      records(13,node)=owners(node)
    end do
    if (rank==root) reference=records
    call MPI_Bcast(reference,size(reference),MPI_INTEGER,root,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call agree(all(records==reference),ready)
    if (.not.ready) return
    valid=.true.;local_nvar=0
    if (scatter) then
      if (rank==root) then
        valid=global%is_valid(tree)
        if (valid) local_nvar=size(global%fields(1)%values,1)
      end if
    else
      valid=allocated(local%fields)
      if (valid) valid=size(local%fields)==size(tree%nodes)
      if (valid) then
        do node=1,size(tree%nodes)
          if (owners(node)/=rank) then
            if (allocated(local%fields(node)%values)) valid=.false.
            cycle
          end if
          if (.not.allocated(local%fields(node)%values)) then
            valid=.false.;cycle
          end if
          n=tree%nodes(node)%shape
          if (local_nvar==0) local_nvar=size(local%fields(node)%values,1)
          if (local_nvar<1) valid=.false.
          if (any(shape(local%fields(node)%values)/=[local_nvar,n])) valid=.false.
          if (.not.all(ieee_is_finite(local%fields(node)%values))) valid=.false.
        end do
      end if
    end if
    call agree(valid,ready)
    if (.not.ready) return
    call MPI_Allreduce(local_nvar,nvar,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=nvar>0
    if (local_nvar>0.and.local_nvar/=nvar) valid=.false.
    do node=1,size(tree%nodes)
      elements=int(nvar,int64)
      do axis=1,3
        if (elements>int(huge(1),int64)/int(tree%nodes(node)%shape(axis),int64)) then
          valid=.false.;exit
        end if
        elements=elements*int(tree%nodes(node)%shape(axis),int64)
      end do
    end do
    call agree(valid,ready)
    if (.not.ready) return
    allocate(candidate%fields(size(tree%nodes)),stat=ios)
    valid=ios==0
    if (valid) then
      do node=1,size(tree%nodes)
        if (scatter) then
          if (owners(node)/=rank) cycle
        else
          if (rank/=root) cycle
        end if
        n=tree%nodes(node)%shape
        allocate(candidate%fields(node)%values(nvar,n(1),n(2),n(3)),stat=ios)
        if (ios/=0) then
          valid=.false.;exit
        end if
      end do
    end if
    call agree(valid,ready)
    if (.not.ready) return
    ! All ranks traverse nodes in the same order. Only root and owner exchange
    ! each payload; no nonowner allocates or receives another patch's values.
    do node=1,size(tree%nodes)
      count=int(int(nvar,int64)*product(int(tree%nodes(node)%shape,int64)))
      ierr=MPI_SUCCESS
      if (scatter) then
        if (owners(node)==root) then
          if (rank==root) candidate%fields(node)%values=global%fields(node)%values
        else if (rank==root) then
          call MPI_Send(global%fields(node)%values,count,MPI_DOUBLE_PRECISION,owners(node),28401,comm,ierr)
        else if (rank==owners(node)) then
          call MPI_Recv(candidate%fields(node)%values,count,MPI_DOUBLE_PRECISION,root,28401,comm,MPI_STATUS_IGNORE,ierr)
        end if
      else
        if (owners(node)==root) then
          if (rank==root) candidate%fields(node)%values=local%fields(node)%values
        else if (rank==owners(node)) then
          call MPI_Send(local%fields(node)%values,count,MPI_DOUBLE_PRECISION,root,28402,comm,ierr)
        else if (rank==root) then
          call MPI_Recv(candidate%fields(node)%values,count,MPI_DOUBLE_PRECISION,owners(node),28402,comm,MPI_STATUS_IGNORE,ierr)
        end if
      end if
      valid=valid.and.ierr==MPI_SUCCESS
    end do
    call agree(valid,ready)
    if (.not.ready) return
    if (scatter) then
      local=candidate
    else if (rank==root) then
      global=candidate
    end if
    ok=.true.
  contains
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: status
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if (status/=MPI_SUCCESS) result=.false.
    end subroutine
  end subroutine
end module mpi_amr_patch_tree_distribution_3d_mod
