module mpi_amr_patch_tree_cells_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  implicit none
  private
  public :: fetch_mpi_amr_tree_cells_3d,validate_mpi_amr_tree_local_3d,validate_mpi_amr_tree_metadata_3d
contains
  ! Each request column is [node,i,j,k], with one-based cell indices.
  ! Request counts/order may differ on every rank, including zero requests.
  subroutine fetch_mpi_amr_tree_cells_3d(tree,owners,local,nvar,requests,values,comm,ok,maximum_depth)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar,requests(:,:)
    type(amr_tree_state_3d), intent(in) :: local
    real(dp), allocatable, intent(inout) :: values(:,:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_depth
    integer, allocatable :: send_count(:),recv_count(:),send_offset(:),recv_offset(:)
    integer, allocatable :: cursor(:),order(:),packed(:,:),received(:,:)
    real(dp), allocatable :: replies(:,:),returned(:,:),candidate(:,:)
    integer :: rank,nranks,status,ios,node,n(3),query,destination,slot,nquery,incoming,peer
    integer(int64) :: total
    logical :: valid,ready
    ok=.false.
    call validate_mpi_amr_tree_local_3d(tree,owners,local,nvar,comm,ready,maximum_depth)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,status)
    if (status/=MPI_SUCCESS) return
    call agree(size(requests,1)==4,ready)
    if (.not.ready) return
    valid=.true.
    nquery=size(requests,2)
    if (nquery>huge(1)/max(4,nvar)) valid=.false.
    do query=1,nquery
      node=requests(1,query)
      if (node<1.or.node>size(tree%nodes)) then
        valid=.false.;cycle
      end if
      if (present(maximum_depth)) then
        if (tree%nodes(node)%level>maximum_depth) valid=.false.
      end if
      if (any(requests(2:4,query)<1).or.any(requests(2:4,query)>tree%nodes(node)%shape)) valid=.false.
    end do
    call agree(valid,ready)
    if (.not.ready) return
    allocate(send_count(nranks),recv_count(nranks),send_offset(nranks),recv_offset(nranks),cursor(nranks), &
      order(nquery),packed(4,nquery),candidate(nvar,nquery),returned(nvar,nquery),stat=ios)
    call agree(ios==0,ready)
    if (.not.ready) return
    send_count=0
    do query=1,nquery
      destination=owners(requests(1,query))+1
      send_count(destination)=send_count(destination)+1
    end do
    call MPI_Alltoall(send_count,1,MPI_INTEGER,recv_count,1,MPI_INTEGER,comm,status)
    if (status/=MPI_SUCCESS) return
    total=sum(int(recv_count,int64))
    call agree(total<=int(huge(1)/max(4,nvar),int64),ready)
    if (.not.ready) return
    incoming=int(total)
    allocate(received(4,incoming),replies(nvar,incoming),stat=ios)
    call agree(ios==0,ready)
    if (.not.ready) return
    send_offset(1)=0;recv_offset(1)=0
    do peer=2,nranks
      send_offset(peer)=send_offset(peer-1)+send_count(peer-1)
      recv_offset(peer)=recv_offset(peer-1)+recv_count(peer-1)
    end do
    cursor=send_offset
    do query=1,nquery
      destination=owners(requests(1,query))+1
      cursor(destination)=cursor(destination)+1;slot=cursor(destination)
      packed(:,slot)=requests(:,query);order(slot)=query
    end do
    call MPI_Alltoallv(packed,4*send_count,4*send_offset,MPI_INTEGER, &
      received,4*recv_count,4*recv_offset,MPI_INTEGER,comm,status)
    if (status/=MPI_SUCCESS) return
    valid=.true.
    do query=1,incoming
      node=received(1,query);n=received(2:4,query)
      ! All coordinates/owners were collectively checked before exchange.
      replies(:,query)=local%fields(node)%values(:,n(1),n(2),n(3))
      if (.not.all(ieee_is_finite(replies(:,query)))) valid=.false.
    end do
    call agree(valid,ready)
    if (.not.ready) return
    call MPI_Alltoallv(replies,nvar*recv_count,nvar*recv_offset,MPI_DOUBLE_PRECISION, &
      returned,nvar*send_count,nvar*send_offset,MPI_DOUBLE_PRECISION,comm,status)
    call agree(status==MPI_SUCCESS,ready)
    if (.not.ready) return
    do query=1,nquery
      candidate(:,order(query))=returned(:,query)
    end do
    call move_alloc(candidate,values);ok=.true.
  contains
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) result=.false.
    end subroutine
  end subroutine
  ! Collective metadata and allocation/shape check; does not scan payload values.
  subroutine validate_mpi_amr_tree_local_3d(tree,owners,local,nvar,comm,ok,maximum_depth)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar
    type(amr_tree_state_3d), intent(in) :: local
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_depth
    integer :: rank,status,node,depth,low,high
    logical :: valid
    call validate_mpi_amr_tree_metadata_3d(tree,owners,nvar,comm,ok)
    if (.not.ok) return
    ok=.false.
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    depth=tree%level_count()-1
    if (present(maximum_depth)) depth=maximum_depth
    call MPI_Allreduce(depth,low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(depth,high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (low/=high.or.depth< -1.or.depth>=tree%level_count()) return
    valid=allocated(local%fields)
    if (valid) valid=size(local%fields)==size(tree%nodes)
    if (valid) then
      do node=1,size(tree%nodes)
        if (tree%nodes(node)%level>depth) cycle
        if (owners(node)/=rank) then
          if (allocated(local%fields(node)%values)) valid=.false.
        else if (.not.allocated(local%fields(node)%values)) then
          valid=.false.
        else
          if (any(shape(local%fields(node)%values)/=[nvar,tree%nodes(node)%shape])) valid=.false.
        end if
      end do
    end if
    call MPI_Allreduce(valid,ok,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if (status/=MPI_SUCCESS) ok=.false.
  end subroutine validate_mpi_amr_tree_local_3d
  ! Collective agreement on topology, owners and component count.
  subroutine validate_mpi_amr_tree_metadata_3d(tree,owners,nvar,comm,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer, allocatable :: records(:,:),reference(:,:)
    integer :: rank,nranks,status,ios,node,header(2),low(2),high(2)
    logical :: valid,ready
    ok=.false.
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,status)
    if (status/=MPI_SUCCESS) return
    valid=tree%is_valid().and.nvar>0
    call agree(valid,ready)
    if (.not.ready) return
    header=[size(tree%nodes),nvar]
    call MPI_Allreduce(header,low,2,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(header,high,2,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(low/=high)) return
    valid=size(owners)==size(tree%nodes).and.size(tree%nodes)<=huge(1)/13
    if (valid) valid=all(owners>=0).and.all(owners<nranks)
    call agree(valid,ready)
    if (.not.ready) return
    allocate(records(13,size(tree%nodes)),reference(13,size(tree%nodes)),stat=ios)
    call agree(ios==0,ready)
    if (.not.ready) return
    records=0;valid=.true.
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
    if (rank==0) reference=records
    call MPI_Bcast(reference,size(reference),MPI_INTEGER,0,comm,status)
    if (status/=MPI_SUCCESS) return
    valid=valid.and.all(records==reference)
    call agree(valid,ok)
  contains
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) result=.false.
    end subroutine
  end subroutine validate_mpi_amr_tree_metadata_3d
end module mpi_amr_patch_tree_cells_3d_mod
