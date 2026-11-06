module mpi_amr_patch_tree_remap_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use amr_patch_tree_remap_3d_mod, only: remap_cell_covered=>covered,remap_gcd=>gcd
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d, &
    validate_mpi_amr_tree_metadata_3d,fetch_mpi_amr_tree_cells_3d
  use mpi_amr_patch_tree_restriction_3d_mod, only: average_down_mpi_amr_tree_3d
  use mpi_amr_patch_tree_thermo_3d_mod, only: recover_mpi_amr_tree_temperatures_3d
  implicit none
  private
  type :: leaf_cell
    integer(int64) :: lower(3),upper(3)
    integer :: node,index(3),slot=0
  end type
  type :: remap_plan
    type(leaf_cell), allocatable :: cells(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer(int64) :: common(3)
    integer, allocatable :: requests(:,:)
    integer :: used=0
  end type
  public :: remap_mpi_amr_tree_3d,remap_mpi_amr_tree_reactive_3d
contains
  subroutine remap_mpi_amr_tree_3d(old_tree,old_owners,old_state,new_tree,new_owners,nvar,new_state,comm,ok)
    type(amr_patch_tree_3d), intent(in) :: old_tree,new_tree
    integer, intent(in) :: old_owners(:),new_owners(:),nvar
    type(amr_tree_state_3d), intent(in) :: old_state
    type(amr_tree_state_3d), intent(inout) :: new_state
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    type(remap_plan) :: plan
    type(amr_tree_state_3d) :: candidate
    real(dp), allocatable :: values(:,:)
    integer :: rank,status,node
    logical :: valid,ready
    ok=.false.
    call validate_mpi_amr_tree_local_3d(old_tree,old_owners,old_state,nvar,comm,ready)
    if (.not.ready) return
    call validate_mpi_amr_tree_metadata_3d(new_tree,new_owners,nvar,comm,ready)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    valid=all(old_tree%nodes(1)%shape==new_tree%nodes(1)%shape)
    do node=1,size(old_tree%nodes)
      if (old_owners(node)/=rank) cycle
      if (.not.all(ieee_is_finite(old_state%fields(node)%values))) valid=.false.
    end do
    call agree(comm,valid,ready)
    if (.not.ready) return
    call build_plan(old_tree,new_tree,new_owners,rank,plan,valid)
    call agree(comm,valid,ready)
    if (.not.ready) return
    ! Only old exposed cells intersecting this rank's new exposed cells are fetched.
    ! Each source cell is requested once per destination rank, in serial source order.
    call fetch_mpi_amr_tree_cells_3d(old_tree,old_owners,old_state,nvar,plan%requests,values,comm,ready)
    if (.not.ready) return
    call evaluate_plan(new_tree,new_owners,rank,nvar,plan,values,candidate,valid)
    call agree(comm,valid,ready)
    if (.not.ready) return
    call average_down_mpi_amr_tree_3d(new_tree,new_owners,candidate,nvar,comm,ready)
    if (.not.ready) return
    call move_alloc(candidate%fields,new_state%fields)
    ok=.true.
  end subroutine

  subroutine remap_mpi_amr_tree_reactive_3d(old_tree,old_owners,old_state,old_temperature,new_tree,new_owners, &
      species,state,temperature,comm,ok,old_reaction_source,reaction_source,old_heat_release,heat_release, &
      old_internal_energy,internal_energy)
    type(amr_patch_tree_3d), intent(in) :: old_tree,new_tree
    integer, intent(in) :: old_owners(:),new_owners(:)
    type(amr_tree_state_3d), intent(in) :: old_state,old_temperature
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(in) :: old_reaction_source,old_heat_release
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source,heat_release
    type(amr_tree_state_3d), optional, intent(in) :: old_internal_energy
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy
    type(amr_tree_state_3d) :: candidate,temps,rates,heat
    type(amr_tree_state_3d), allocatable :: held
    integer :: rank,status,node,flags,low,high
    logical :: valid,ready
    ok=.false.
    flags=merge(1,0,present(old_reaction_source))+2*merge(1,0,present(reaction_source))+ &
      4*merge(1,0,present(old_heat_release))+8*merge(1,0,present(heat_release))+ &
      16*merge(1,0,present(old_internal_energy))+32*merge(1,0,present(internal_energy))
    call MPI_Allreduce(flags,low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.low/=high) return
    if(present(old_internal_energy).neqv.present(internal_energy)) return
    if(present(old_internal_energy)) then
      allocate(held)
      call remap_mpi_amr_tree_3d(old_tree,old_owners,old_internal_energy,new_tree,new_owners,1,held,comm,ready)
      if(.not.ready) return
    endif
    if(present(old_reaction_source).neqv.present(reaction_source)) return
    if(present(old_heat_release).neqv.present(heat_release)) return
    if(present(old_reaction_source)) then
      call remap_mpi_amr_tree_3d(old_tree,old_owners,old_reaction_source,new_tree,new_owners, &
        reactive_nvar(size(species))-4,rates,comm,ready)
      if(.not.ready) return
    end if
    if(present(old_heat_release)) then
      call remap_mpi_amr_tree_3d(old_tree,old_owners,old_heat_release,new_tree,new_owners,1,heat,comm,ready)
      if(.not.ready) return
    end if
    call validate_mpi_amr_tree_local_3d(old_tree,old_owners,old_temperature,1,comm,ready)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    valid=size(species)>0
    do node=1,size(old_tree%nodes)
      if (old_owners(node)/=rank) cycle
      if (.not.all(ieee_is_finite(old_temperature%fields(node)%values))) then
        valid=.false.
      else if (any(old_temperature%fields(node)%values<=0)) then
        valid=.false.
      end if
    end do
    call agree(comm,valid,ready)
    if (.not.ready) return
    call remap_mpi_amr_tree_3d(old_tree,old_owners,old_state,new_tree,new_owners,reactive_nvar(size(species)), &
      candidate,comm,ready)
    if (.not.ready) return
    call remap_mpi_amr_tree_3d(old_tree,old_owners,old_temperature,new_tree,new_owners,1,temps,comm,ready)
    if (.not.ready) return
    call recover_mpi_amr_tree_temperatures_3d(new_tree,new_owners,species,candidate,temps,comm,ready,internal_energy=held)
    if (.not.ready) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(temps%fields,temperature%fields)
    if(present(internal_energy)) call move_alloc(held%fields,internal_energy%fields)
    if(present(reaction_source)) call move_alloc(rates%fields,reaction_source%fields)
    if(present(heat_release)) call move_alloc(heat%fields,heat_release%fields)
    ok=.true.
  end subroutine

  subroutine build_plan(old_tree,new_tree,owners,rank,plan,ok)
    type(amr_patch_tree_3d), intent(in) :: old_tree,new_tree
    integer, intent(in) :: owners(:),rank
    type(remap_plan), intent(out) :: plan
    logical, intent(out) :: ok
    integer(int64), allocatable :: old_origin(:,:),old_domain(:,:)
    integer, allocatable :: old_ratios(:),new_ratios(:)
    integer(int64) :: a,b,g,total,count,factor(3),lo(3),hi(3),overlap(3)
    integer :: axis,node,n(3),level,i,j,k,idx(3),cell,ios,slot
    logical :: valid
    ok=.false.
    call amr_tree_integer_layout_3d(old_tree,old_origin,old_domain,old_ratios,valid)
    if (.not.valid) return
    call amr_tree_integer_layout_3d(new_tree,plan%origin,plan%domain,new_ratios,valid)
    if (.not.valid) return
    do axis=1,3
      a=old_domain(axis,size(old_ratios));b=plan%domain(axis,size(new_ratios))
      g=remap_gcd(a,b)
      if (a/g>huge(1_int64)/b) return
      plan%common(axis)=(a/g)*b
    end do
    total=0
    do node=1,size(old_tree%nodes)
      count=1
      do axis=1,3
        if (count>int(huge(1),int64)/old_tree%nodes(node)%shape(axis)) return
        count=count*old_tree%nodes(node)%shape(axis)
      end do
      if (count>int(huge(1),int64)-total) return
      total=total+count
    end do
    allocate(plan%cells(int(total)),stat=ios)
    if (ios/=0) return
    do node=1,size(old_tree%nodes)
      n=old_tree%nodes(node)%shape;level=old_tree%nodes(node)%level+1
      factor=plan%common/old_domain(:,level)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            idx=[i,j,k]
            if (remap_cell_covered(old_tree,node,idx)) cycle
            plan%used=plan%used+1
            associate(item=>plan%cells(plan%used))
              item%node=node;item%index=idx
              item%lower=(old_origin(:,node)+int(idx-1,int64))*factor
              item%upper=item%lower+factor
            end associate
          end do
        end do
      end do
    end do
    do node=1,size(new_tree%nodes)
      if (owners(node)/=rank) cycle
      n=new_tree%nodes(node)%shape;factor=plan%common/plan%domain(:,new_tree%nodes(node)%level+1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            idx=[i,j,k]
            if (remap_cell_covered(new_tree,node,idx)) cycle
            lo=(plan%origin(:,node)+int(idx-1,int64))*factor;hi=lo+factor
            do cell=1,plan%used
              overlap=min(hi,plan%cells(cell)%upper)-max(lo,plan%cells(cell)%lower)
              if (all(overlap>0)) plan%cells(cell)%slot=1
            end do
          end do
        end do
      end do
    end do
    slot=0
    do cell=1,plan%used
      if (plan%cells(cell)%slot==0) cycle
      slot=slot+1;plan%cells(cell)%slot=slot
    end do
    if (slot>huge(1)/4) return
    allocate(plan%requests(4,slot),stat=ios)
    if (ios/=0) return
    do cell=1,plan%used
      slot=plan%cells(cell)%slot
      if (slot==0) cycle
      plan%requests(:,slot)=[plan%cells(cell)%node,plan%cells(cell)%index]
    end do
    ok=.true.
  end subroutine

  subroutine evaluate_plan(tree,owners,rank,nvar,plan,values,candidate,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),rank,nvar
    type(remap_plan), intent(in) :: plan
    real(dp), intent(in) :: values(:,:)
    type(amr_tree_state_3d), intent(out) :: candidate
    logical, intent(out) :: ok
    real(dp), allocatable :: value(:),compensation(:),term(:),increment(:),next_value(:)
    real(dp) :: weight,coverage,coverage_comp,delta,next_coverage
    integer(int64) :: factor(3),lo(3),hi(3),overlap(3)
    integer :: node,n(3),i,j,k,cell,ios,idx(3)
    ok=.false.
    allocate(candidate%fields(size(tree%nodes)),value(nvar),compensation(nvar),term(nvar), &
      increment(nvar),next_value(nvar),stat=ios)
    if (ios/=0) return
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape;factor=plan%common/plan%domain(:,tree%nodes(node)%level+1)
      allocate(candidate%fields(node)%values(nvar,n(1),n(2),n(3)),stat=ios)
      if (ios/=0) return
      candidate%fields(node)%values=0
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            idx=[i,j,k]
            if (remap_cell_covered(tree,node,idx)) cycle
            lo=(plan%origin(:,node)+int(idx-1,int64))*factor;hi=lo+factor
            value=0;compensation=0;coverage=0;coverage_comp=0
            do cell=1,plan%used
              overlap=min(hi,plan%cells(cell)%upper)-max(lo,plan%cells(cell)%lower)
              if (any(overlap<=0)) cycle
              if (plan%cells(cell)%slot==0) return
              weight=product(real(overlap,dp)/real(factor,dp))
              term=weight*values(:,plan%cells(cell)%slot)
              increment=term-compensation;next_value=value+increment
              compensation=(next_value-value)-increment;value=next_value
              delta=weight-coverage_comp;next_coverage=coverage+delta
              coverage_comp=(next_coverage-coverage)-delta;coverage=next_coverage
            end do
            if (abs(coverage-1.0_dp)>128*epsilon(1.0_dp)) return
            candidate%fields(node)%values(:,i,j,k)=value/coverage
          end do
        end do
      end do
      if (.not.all(ieee_is_finite(candidate%fields(node)%values))) return
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
