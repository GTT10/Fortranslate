module mpi_amr_patch_tree_tagging_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use amr_patch_tree_tagging_3d_mod, only: hot_cell,build_tagged_tree_3d
  use amr_patch_tree_remap_3d_mod, only: gcd,covered
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  implicit none
  private
  public :: tag_mpi_amr_tree_temperature_3d
contains
  subroutine tag_mpi_amr_tree_temperature_3d(old_tree,owners,temperature,thresholds,ratios,tile_shape,buffer, &
      maximum_nodes,periodic,new_tree,comm,ok)
    type(amr_patch_tree_3d), intent(in) :: old_tree
    integer, intent(in) :: owners(:),ratios(:),tile_shape(3),buffer,maximum_nodes
    type(amr_tree_state_3d), intent(in) :: temperature
    real(dp), intent(in) :: thresholds(:)
    logical, intent(in) :: periodic(3)
    type(amr_patch_tree_3d), intent(inout) :: new_tree
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    type(amr_patch_tree_3d) :: candidate
    type(hot_cell), allocatable :: local(:),cells(:)
    integer, allocatable :: counts(:),offsets(:),controls(:),low(:),high(:)
    real(dp), allocatable :: minimum(:),maximum(:),sent_t(:),received_t(:)
    integer(int64), allocatable :: sent_boxes(:,:),received_boxes(:,:)
    integer(int64) :: common(3),total
    integer :: rank,nranks,status,ios,levels,lo,hi,node,peer,nlocal
    logical :: valid,ready
    ok=.false.
    call validate_mpi_amr_tree_local_3d(old_tree,owners,temperature,1,comm,ready)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,status)
    if (status/=MPI_SUCCESS) return
    levels=size(ratios)
    call MPI_Allreduce(levels,lo,1,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(levels,hi,1,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    call agree(comm,lo==hi.and.size(thresholds)==levels.and.levels<=huge(1)-8,ready)
    if (.not.ready) return
    call agree(comm,all(ieee_is_finite(thresholds)),ready)
    if (.not.ready) return
    allocate(controls(levels+8),low(levels+8),high(levels+8),minimum(levels),maximum(levels), &
      counts(nranks),offsets(nranks),stat=ios)
    call agree(comm,ios==0,ready)
    if (.not.ready) return
    controls=[ratios,tile_shape,buffer,maximum_nodes,merge(1,0,periodic)]
    call MPI_Allreduce(controls,low,size(controls),MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,size(controls),MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(low/=high)) return
    call MPI_Allreduce(thresholds,minimum,levels,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(thresholds,maximum,levels,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(minimum/=maximum)) return
    call collect_local(old_tree,owners,rank,temperature,thresholds,ratios,tile_shape,buffer,maximum_nodes, &
      local,common,valid)
    call agree(comm,valid,ready)
    if (.not.ready) return
    nlocal=size(local)
    call MPI_Allgather(nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,status)
    if (status/=MPI_SUCCESS) return
    total=sum(int(counts,int64))
    if (total>int(huge(1)/6,int64)) return
    offsets(1)=0
    do peer=2,nranks
      offsets(peer)=offsets(peer-1)+counts(peer-1)
    end do
    allocate(cells(int(total)),sent_boxes(6,nlocal),received_boxes(6,int(total)), &
      sent_t(nlocal),received_t(int(total)),stat=ios)
    call agree(comm,ios==0,ready)
    if (.not.ready) return
    do node=1,nlocal
      sent_boxes(:,node)=[local(node)%lower,local(node)%upper]
      sent_t(node)=local(node)%temperature
    end do
    call MPI_Allgatherv(sent_t,nlocal,MPI_DOUBLE_PRECISION,received_t,counts,offsets, &
      MPI_DOUBLE_PRECISION,comm,status)
    if (status/=MPI_SUCCESS) return
    counts=counts*6;offsets=offsets*6
    call MPI_Allgatherv(sent_boxes,6*nlocal,MPI_INTEGER8,received_boxes,counts,offsets,MPI_INTEGER8,comm,status)
    if (status/=MPI_SUCCESS) return
    do node=1,int(total)
      cells(node)%lower=received_boxes(1:3,node);cells(node)%upper=received_boxes(4:6,node)
      cells(node)%temperature=received_t(node)
    end do
    ! Tags are Boolean votes: gathered source order does not change deterministic tile order.
    call build_tagged_tree_3d(old_tree,cells,common,thresholds,ratios,tile_shape,buffer, &
      maximum_nodes,periodic,candidate,valid)
    call agree(comm,valid,ready)
    if (.not.ready) return
    call move_alloc(candidate%nodes,new_tree%nodes)
    ok=.true.
  end subroutine

  subroutine collect_local(old_tree,owners,rank,temperature,thresholds,ratios,tile_shape,buffer,maximum_nodes, &
      cells,common,ok)
    type(amr_patch_tree_3d), intent(in) :: old_tree
    integer, intent(in) :: owners(:),rank,ratios(:),tile_shape(3),buffer,maximum_nodes
    type(amr_tree_state_3d), intent(in) :: temperature
    real(dp), intent(in) :: thresholds(:)
    logical, intent(out) :: ok
    type(hot_cell), allocatable, intent(out) :: cells(:)
    integer(int64), allocatable :: old_origin(:,:),old_domain(:,:)
    integer, allocatable :: old_ratios(:)
    integer(int64), intent(out) :: common(3)
    integer(int64) :: target(3),factor(3),a,b,g,total
    integer :: axis,level,node,n(3),i,j,k,idx(3),used,pass
    logical :: valid
    ok=.false.
    if (size(thresholds)/=size(ratios)) return
    if (any(.not.ieee_is_finite(thresholds))) return
    if (any(thresholds<=0)) return
    if (any(ratios<2).or.any(tile_shape<1).or.buffer<0.or.maximum_nodes<1) return
    do level=2,size(thresholds)
      if (thresholds(level)<thresholds(level-1)) return
    end do
    do node=1,size(old_tree%nodes)
      if (owners(node)/=rank) cycle
      if (.not.all(ieee_is_finite(temperature%fields(node)%values))) return
      if (any(temperature%fields(node)%values<=0)) return
    end do
    call amr_tree_integer_layout_3d(old_tree,old_origin,old_domain,old_ratios,valid)
    if (.not.valid) return
    ! One exact integer lattice for both old cells and all requested new depths.
    ! Reserve headroom for buffered and periodically translated tile bounds.
    target=int(old_tree%nodes(1)%shape,int64)
    do level=1,size(ratios)
      if (any(target>(huge(1_int64)/4)/int(ratios(level),int64))) return
      target=target*int(ratios(level),int64)
    end do
    do axis=1,3
      a=old_domain(axis,size(old_ratios));b=target(axis);g=gcd(a,b)
      if (a/g>(huge(1_int64)/4)/b) return
      common(axis)=(a/g)*b
    end do
    ! Only exposed old cells vote. Covered parent temperatures may be stale;
    ! averaging a hot fine cell down before tagging would lose that feature.
    total=0
    do pass=1,2
      used=0
      do node=1,size(old_tree%nodes)
        if (owners(node)/=rank) cycle
        n=old_tree%nodes(node)%shape;level=old_tree%nodes(node)%level+1
        factor=common/old_domain(:,level)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              if (size(thresholds)==0) cycle
              if (temperature%fields(node)%values(1,i,j,k)<thresholds(1)) cycle
              idx=[i,j,k]
              if (covered(old_tree,node,idx)) cycle
              if (pass==1) then
                total=total+1
                if (total>int(huge(1),int64)) return
              else
                used=used+1
                cells(used)%lower=(old_origin(:,node)+int(idx-1,int64))*factor
                cells(used)%upper=cells(used)%lower+factor
                cells(used)%temperature=temperature%fields(node)%values(1,i,j,k)
              end if
            end do
          end do
        end do
      end do
      if (pass==1) then
        allocate(cells(int(total)),stat=i)
        if (i/=0) return
      end if
    end do
    ok=.true.
  end subroutine

  subroutine agree(comm,condition,result)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(in) :: condition
    logical, intent(out) :: result
    integer :: status
    call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if (status/=MPI_SUCCESS) result=.false.
  end subroutine
end module
