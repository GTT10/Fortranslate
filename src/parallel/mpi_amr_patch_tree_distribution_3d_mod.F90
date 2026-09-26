module mpi_amr_patch_tree_distribution_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d,validate_mpi_amr_tree_metadata_3d
  implicit none
  private
  public :: redistribute_mpi_amr_tree_3d,rebalance_mpi_amr_tree_3d,rebalance_mpi_amr_tree_with_field_3d
  public :: scatter_mpi_amr_tree_3d,gather_mpi_amr_tree_3d,partition_amr_tree_3d
contains
  ! Rebalance an additional signed cell field atomically with gas and temperature.
  ! Packing reuses the whole-patch transfer and never interpolates source estimates.
  subroutine rebalance_mpi_amr_tree_with_field_3d(tree,owners,nvar,state,temperature,nfield,field, &
      local_cost,minimum_gain,comm,changed,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(inout) :: owners(:)
    integer, intent(in) :: nvar,nfield
    type(amr_tree_state_3d), intent(inout) :: state,temperature,field
    real(dp), intent(in) :: local_cost(:),minimum_gain
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: changed,ok
    type(amr_tree_state_3d) :: packed,temps,gas,extra
    integer, allocatable :: next_owners(:)
    integer :: rank,status,node,n(3),ios
    logical :: valid,ready,moved
    changed=.false.;ok=.false.
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nvar,comm,ready)
    if(.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,field,nfield,comm,ready)
    if(.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,ready)
    if(.not.ready) return
    if(nfield>huge(1)-nvar) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    allocate(packed%fields(size(owners)),gas%fields(size(owners)),extra%fields(size(owners)),stat=ios)
    valid=ios==0
    call agree(valid,ready)
    if(.not.ready) return
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      if(.not.all(ieee_is_finite(field%fields(node)%values))) then
        valid=.false.;exit
      end if
      allocate(packed%fields(node)%values(nvar+nfield,n(1),n(2),n(3)),stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      end if
      packed%fields(node)%values(:nvar,:,:,:)=state%fields(node)%values
      packed%fields(node)%values(nvar+1:,:,:,:)=field%fields(node)%values
    end do
    call agree(valid,ready)
    if(.not.ready) return
    temps=temperature;next_owners=owners
    call rebalance_mpi_amr_tree_3d(tree,next_owners,nvar+nfield,packed,temps,local_cost,minimum_gain,comm,moved,ready)
    if(.not.ready) return
    if(.not.moved) then
      ok=.true.;return
    end if
    do node=1,size(owners)
      if(next_owners(node)/=rank) cycle
      allocate(gas%fields(node)%values,source=packed%fields(node)%values(:nvar,:,:,:),stat=ios)
      if(ios==0) allocate(extra%fields(node)%values,source=packed%fields(node)%values(nvar+1:,:,:,:),stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      end if
    end do
    call agree(valid,ready)
    if(.not.ready) return
    call move_alloc(gas%fields,state%fields)
    call move_alloc(temps%fields,temperature%fields)
    call move_alloc(extra%fields,field%fields)
    owners=next_owners;changed=.true.;ok=.true.
  contains
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: code
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,code)
      if(code/=MPI_SUCCESS) result=.false.
    end subroutine
  end subroutine rebalance_mpi_amr_tree_with_field_3d

  ! Owner-local costs must describe the same completed sampling interval. Only
  ! owner entries are inspected; this API neither times kernels nor schedules calls.
  subroutine rebalance_mpi_amr_tree_3d(tree,owners,nvar,state,temperature,local_cost,minimum_gain,comm,changed,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(inout) :: owners(:)
    integer, intent(in) :: nvar
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: local_cost(:),minimum_gain
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: changed,ok
    real(dp), allocatable :: contribution(:),cost(:),work(:),old_load(:),new_load(:)
    integer, allocatable :: candidate(:)
    real(dp) :: lower,upper,scale,peak,gain
    integer :: rank,nranks,ierr,ios,node,next,destination,iteration
    logical :: valid,agreed
    changed=.false.;ok=.false.
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nvar,comm,valid)
    if (.not.valid) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,valid)
    if (.not.valid) return
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=size(local_cost)==size(owners).and.ieee_is_finite(minimum_gain)
    if (valid) valid=minimum_gain>=0._dp.and.minimum_gain<1._dp
    if (valid) then
      do node=1,size(owners)
        if (owners(node)/=rank) cycle
        if (.not.ieee_is_finite(local_cost(node))) valid=.false.
        if (local_cost(node)<0._dp) valid=.false.
      end do
    end if
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(minimum_gain,lower,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(minimum_gain,upper,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (lower/=upper) return
    allocate(contribution(size(owners)),cost(size(owners)),work(size(owners)), &
      candidate(size(owners)),old_load(nranks),new_load(nranks),stat=ios)
    valid=ios==0
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    contribution=0._dp
    do node=1,size(owners)
      if (owners(node)==rank) contribution(node)=local_cost(node)
    end do
    call MPI_Allreduce(contribution,cost,size(cost),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    scale=maxval(cost)
    if (scale==0._dp) then
      ok=.true.;return
    end if
    ! Normalize before summation so finite large timings cannot overflow rank loads.
    cost=cost/scale;work=cost;old_load=0._dp;new_load=0._dp
    candidate=owners
    do node=1,size(owners)
      old_load(owners(node)+1)=old_load(owners(node)+1)+cost(node)
    end do
    do iteration=1,size(owners)
      next=maxloc(work,dim=1)
      if (work(next)<=0._dp) exit
      destination=minloc(new_load,dim=1)
      candidate(next)=destination-1
      new_load(destination)=new_load(destination)+work(next)
      work(next)=-1._dp
    end do
    if(all(candidate==owners)) then
      ok=.true.;return
    end if
    peak=maxval(old_load)
    gain=(peak-maxval(new_load))/peak
    if (gain<=0._dp.or.gain<minimum_gain) then
      ok=.true.;return
    end if
    call redistribute_mpi_amr_tree_3d(tree,owners,candidate,nvar,state,temperature,comm,valid)
    if (.not.valid) return
    owners=candidate
    changed=.true.;ok=.true.
  end subroutine

  ! Move U/T together between whole-patch owners, without interpolation or a root gather.
  ! Geometry and both owner maps must agree across the borrowed communicator.
  subroutine redistribute_mpi_amr_tree_3d(tree,old_owners,new_owners,nvar,state,temperature,comm,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: old_owners(:),new_owners(:),nvar
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: candidate,temps
    integer :: rank,status,node,axis,n(3),ios,count,temperature_count
    integer(int64) :: elements
    logical :: valid
    ok=.false.
    call validate_mpi_amr_tree_local_3d(tree,old_owners,state,nvar,comm,valid)
    if (.not.valid) return
    call validate_mpi_amr_tree_local_3d(tree,old_owners,temperature,1,comm,valid)
    if (.not.valid) return
    call validate_mpi_amr_tree_metadata_3d(tree,new_owners,nvar,comm,valid)
    if (.not.valid) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    valid=.true.
    do node=1,size(tree%nodes)
      elements=int(nvar,int64)
      do axis=1,3
        if (elements>int(huge(1),int64)/int(tree%nodes(node)%shape(axis),int64)) then
          valid=.false.;exit
        end if
        elements=elements*int(tree%nodes(node)%shape(axis),int64)
      end do
      if (old_owners(node)/=rank) cycle
      if (.not.all(ieee_is_finite(state%fields(node)%values))) valid=.false.
      if (.not.all(ieee_is_finite(temperature%fields(node)%values))) valid=.false.
      if (any(temperature%fields(node)%values<=0)) valid=.false.
    end do
    call agree_move(valid)
    if (.not.valid) return
    allocate(candidate%fields(size(tree%nodes)),temps%fields(size(tree%nodes)),stat=ios)
    valid=ios==0
    if (valid) then
      do node=1,size(tree%nodes)
        if (new_owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        allocate(candidate%fields(node)%values(nvar,n(1),n(2),n(3)), &
          temps%fields(node)%values(1,n(1),n(2),n(3)),stat=ios)
        if (ios/=0) then
          valid=.false.;exit
        end if
      end do
    end if
    call agree_move(valid)
    if (.not.valid) return
    do node=1,size(tree%nodes)
      temperature_count=int(product(int(tree%nodes(node)%shape,int64)))
      count=nvar*temperature_count
      if (old_owners(node)==new_owners(node)) then
        if (rank==old_owners(node)) then
          candidate%fields(node)%values=state%fields(node)%values
          temps%fields(node)%values=temperature%fields(node)%values
        end if
        cycle
      end if
      status=MPI_SUCCESS
      if (rank==old_owners(node)) then
        call MPI_Send(state%fields(node)%values,count,MPI_DOUBLE_PRECISION,new_owners(node),28407,comm,status)
      else if (rank==new_owners(node)) then
        call MPI_Recv(candidate%fields(node)%values,count,MPI_DOUBLE_PRECISION,old_owners(node),28407, &
          comm,MPI_STATUS_IGNORE,status)
      end if
      valid=valid.and.status==MPI_SUCCESS
      status=MPI_SUCCESS
      if (rank==old_owners(node)) then
        call MPI_Send(temperature%fields(node)%values,temperature_count,MPI_DOUBLE_PRECISION,new_owners(node), &
          28408,comm,status)
      else if (rank==new_owners(node)) then
        call MPI_Recv(temps%fields(node)%values,temperature_count,MPI_DOUBLE_PRECISION,old_owners(node), &
          28408,comm,MPI_STATUS_IGNORE,status)
      end if
      valid=valid.and.status==MPI_SUCCESS
    end do
    call agree_move(valid)
    if (.not.valid) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(temps%fields,temperature%fields)
    ok=.true.
  contains
    subroutine agree_move(condition)
      logical, intent(inout) :: condition
      logical :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      condition=.false.
      if (ierr==MPI_SUCCESS) condition=result
    end subroutine
  end subroutine

  ! Deterministic largest-work-first whole-patch assignment. Work is a cell-step
  ! proxy, not a measured runtime: exponent 0=common, 1=hydro, 2=transport subcycling.
  ! Replicated geometry/controls are a caller contract; no field payload is read.
  subroutine partition_amr_tree_3d(tree,nranks,exponent,owners,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: nranks,exponent
    integer, allocatable, intent(inout) :: owners(:)
    logical, intent(out) :: ok
    integer(int64), allocatable :: work(:),counts(:),load(:)
    integer, allocatable :: candidate(:)
    integer(int64) :: factor,cells
    integer :: node,axis,parent,rank,next,ios,iteration
    ok=.false.
    if (.not.tree%is_valid().or.nranks<1.or.exponent<0.or.exponent>2) return
    allocate(work(size(tree%nodes)),counts(size(tree%nodes)),load(nranks), &
      candidate(size(tree%nodes)),stat=ios)
    if (ios/=0) return
    counts(1)=1;load=0
    do node=1,size(tree%nodes)
      if (node>1) then
        parent=tree%nodes(node)%parent
        factor=int(tree%nodes(node)%patch%refinement_ratio,int64)**exponent
        if (counts(parent)>huge(1_int64)/factor) return
        counts(node)=counts(parent)*factor
      end if
      cells=1
      do axis=1,3
        if (cells>huge(1_int64)/int(tree%nodes(node)%shape(axis),int64)) return
        cells=cells*int(tree%nodes(node)%shape(axis),int64)
      end do
      if (cells>huge(1_int64)/counts(node)) return
      work(node)=cells*counts(node)
    end do
    ! MAXLOC/MINLOC break ties by lower node/rank index, independent of MPI size history.
    do iteration=1,size(tree%nodes)
      next=maxloc(work,dim=1);rank=minloc(load,dim=1)
      if (load(rank)>huge(1_int64)-work(next)) return
      candidate(next)=rank-1
      load(rank)=load(rank)+work(next);work(next)=-1
    end do
    call move_alloc(candidate,owners)
    ok=.true.
  end subroutine
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
