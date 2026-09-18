module mpi_amr_patch_tree_subcycle_3d_mod
  use mpi_f08
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use mpi_amr_patch_tree_chemistry_3d_mod, only: clone_mpi_amr_tree_inputs_3d
  use mpi_amr_patch_tree_restriction_3d_mod, only: average_down_mpi_amr_tree_3d
  use mpi_amr_patch_tree_thermo_3d_mod, only: recover_mpi_amr_tree_temperatures_3d
  use mpi_amr_patch_tree_reflux_3d_mod, only: reflux_mpi_amr_tree_3d
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use reactive_3d_mod, only: compute_reactive_cfl_timestep_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use mpi_amr_patch_tree_hydro_level_3d_mod, only: advance_mpi_amr_tree_hydro_level_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  implicit none
  private
  public :: compute_mpi_amr_tree_subcycle_cfl_3d,advance_mpi_amr_tree_subcycle_hydro_3d
contains
  subroutine compute_mpi_amr_tree_subcycle_cfl_3d(tree,owners,species,state,temperature,spacing,cfl,dt,comm,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: spacing(3),cfl
    real(dp), intent(inout) :: dt
    logical, intent(out) :: ok
    integer(int64), allocatable :: origin(:,:),domain(:,:),counts(:)
    integer, allocatable :: ratios(:)
    real(dp) :: candidate,local_dt,h(3),count_real
    integer :: node,n(3),level,rank,status
    real(dp) :: low,high,global_dt
    logical :: valid
    ok=.false.
    call layout(tree,owners,species,state,temperature,spacing,origin,domain,ratios,counts,comm,valid)
    if (.not.valid) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    valid=ieee_is_finite(cfl)
    if (valid) valid=cfl>0
    call agree(comm,valid)
    if (.not.valid) return
    call MPI_Allreduce(cfl,low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(cfl,high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (low/=high) return
    candidate=huge(1._dp)
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      level=tree%nodes(node)%level+1;n=tree%nodes(node)%shape
      count_real=real(counts(level),dp);h=spacing/count_real
      call compute_reactive_cfl_timestep_3d(species,state%fields(node)%values,temperature%fields(node)%values(1,:,:,:), &
        n(1),n(2),n(3),h(1),h(2),h(3),cfl,local_dt,valid)
      if (.not.valid) exit
      if (local_dt<=huge(1._dp)/count_real) candidate=min(candidate,local_dt*count_real)
    end do
    call agree(comm,valid)
    if (.not.valid) return
    call MPI_Allreduce(candidate,global_dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    dt=global_dt;ok=.true.
  end subroutine

  subroutine advance_mpi_amr_tree_subcycle_hydro_3d(tree,owners,species,state,temperature,spacing,dt,faces, &
      reconstruction,limiter,solver,boundary_change,maximum_correction,comm,ok,inflow_state,inflow_temperature, &
      level_steps,outlet_pressure)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: spacing(3),dt
    character(len=*), intent(in) :: faces(6),reconstruction,limiter,solver
    real(dp), intent(inout) :: boundary_change(:),maximum_correction
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:)
    real(dp), optional, intent(in) :: outlet_pressure(6)
    integer(int64), optional, intent(inout) :: level_steps(:)
    type(amr_tree_state_3d) :: candidate,temps,old_endpoint,new_endpoint
    type(amr_tree_flux_3d), allocatable :: integrated(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:),counts(:),executed(:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: starts(:),finishes(:)
    real(dp) :: budget(size(boundary_change)),largest,dt_low,dt_high
    integer :: nvar,node,n(3),rank,status,ios
    logical :: valid
    ok=.false.
    call layout(tree,owners,species,state,temperature,spacing,origin,domain,ratios,counts,comm,valid)
    if (.not.valid) return
    nvar=reactive_nvar(size(species))
    valid=size(boundary_change)==nvar.and.ieee_is_finite(dt)
    if (valid) valid=dt>0
    if (present(level_steps)) valid=valid.and.size(level_steps)==size(counts)
    call agree(comm,valid)
    if (.not.valid) return
    call MPI_Allreduce(dt,dt_low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(dt,dt_high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (dt_low/=dt_high) return
    call clone_mpi_amr_tree_inputs_3d(tree,owners,size(species),state,temperature, &
      candidate,temps,comm,rank,valid)
    if (.not.valid) return
    call average_down_mpi_amr_tree_3d(tree,owners,candidate,nvar,comm,valid)
    if (.not.valid) return
    call recover_from(0,valid)
    if (.not.valid) return
    allocate(old_endpoint%fields(size(tree%nodes)),new_endpoint%fields(size(tree%nodes)),stat=ios)
    valid=ios==0;call agree(comm,valid)
    if (.not.valid) return
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      allocate(old_endpoint%fields(node)%values(nvar+1,n(1),n(2),n(3)),stat=ios)
      if (ios/=0) exit
      allocate(new_endpoint%fields(node)%values(nvar+1,n(1),n(2),n(3)),stat=ios)
      if (ios/=0) exit
    end do
    valid=ios==0;call agree(comm,valid)
    if (.not.valid) return
    allocate(starts(size(counts)),finishes(size(counts)),executed(size(counts)),stat=ios)
    valid=ios==0;call agree(comm,valid)
    if (.not.valid) return
    starts=0;finishes=1;executed=0;budget=0;largest=0
    call interval(0,0._dp,1._dp,integrated,valid)
    if (.not.valid) return
    if (any(executed/=counts)) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(temps%fields,temperature%fields)
    boundary_change=budget;maximum_correction=largest
    if (present(level_steps)) level_steps=executed
    ok=.true.
  contains
    recursive subroutine interval(depth,begin_time,end_time,effective,success)
      integer, intent(in) :: depth
      real(dp), intent(in) :: begin_time,end_time
      type(amr_tree_flux_3d), allocatable, intent(out) :: effective(:)
      logical, intent(out) :: success
      type(amr_tree_flux_3d), allocatable :: total(:),child_flux(:)
      real(dp) :: fraction0(depth),fraction1(depth),correction(nvar),max_change,sub_begin,sub_end,local_dt
      integer :: peer,ancestor,substep,ratio,allocation_status
      logical :: good
      success=.false.
      if (end_time<=begin_time) return
      local_dt=dt*(end_time-begin_time)
      if (local_dt<=0) return
      starts(depth+1)=begin_time;finishes(depth+1)=end_time
      do ancestor=1,depth
        fraction0(ancestor)=max(0._dp,min(1._dp,(begin_time-starts(ancestor))/(finishes(ancestor)-starts(ancestor))))
        fraction1(ancestor)=max(0._dp,min(1._dp,(end_time-starts(ancestor))/(finishes(ancestor)-starts(ancestor))))
      end do
      call pack_depth(depth,old_endpoint)
      call advance_mpi_amr_tree_hydro_level_3d(tree,owners,species,candidate,temps,depth,spacing,local_dt,faces, &
        reconstruction,limiter,solver,total,comm,good,old_endpoint,new_endpoint,fraction0,fraction1, &
        inflow_state,inflow_temperature,outlet_pressure)
      if (.not.good) return
      executed(depth+1)=executed(depth+1)+1
      call pack_depth(depth,new_endpoint)
      if (depth==0) then
        if (rank==owners(1)) call root_budget(total(1))
        call MPI_Bcast(budget,nvar,MPI_DOUBLE_PRECISION,owners(1),comm,status)
        if (status/=MPI_SUCCESS) return
      end if
      if (depth+1<size(counts)) then
        ratio=ratios(depth+2)
        sub_begin=begin_time
        do substep=1,ratio
          sub_end=begin_time+(end_time-begin_time)*real(substep,dp)/real(ratio,dp)
          if (substep==ratio) sub_end=end_time
          call interval(depth+1,sub_begin,sub_end,child_flux,good)
          if (.not.good) return
          allocation_status=0
          do peer=1,size(tree%nodes)
            if (owners(peer)/=rank) cycle
            if (tree%nodes(peer)%level<=depth) cycle
            if (.not.allocated(total(peer)%x)) then
              allocate(total(peer)%x,source=child_flux(peer)%x,stat=allocation_status)
              if (allocation_status/=0) exit
              allocate(total(peer)%y,source=child_flux(peer)%y,stat=allocation_status)
              if (allocation_status/=0) exit
              allocate(total(peer)%z,source=child_flux(peer)%z,stat=allocation_status)
              if (allocation_status/=0) exit
            else
              total(peer)%x=total(peer)%x+child_flux(peer)%x
              total(peer)%y=total(peer)%y+child_flux(peer)%y
              total(peer)%z=total(peer)%z+child_flux(peer)%z
            end if
          end do
          good=allocation_status==0;call agree(comm,good)
          if (.not.good) return
          sub_begin=sub_end
        end do
      end if
      correction=0;max_change=0
      call reflux_mpi_amr_tree_3d(tree,owners,candidate,nvar,total,spacing,faces,correction,max_change,comm,good, &
        minimum_depth=depth,synchronized_flux=effective)
      if (.not.good) return
      largest=max(largest,max_change)
      if (depth==0) budget=budget+correction
      call recover_from(depth,good)
      if (.not.good) return
      success=.true.
    end subroutine
    subroutine pack_depth(depth,packed)
      integer, intent(in) :: depth
      type(amr_tree_state_3d), intent(inout) :: packed
      integer :: peer
      do peer=1,size(tree%nodes)
        if (owners(peer)/=rank.or.tree%nodes(peer)%level/=depth) cycle
        packed%fields(peer)%values(:nvar,:,:,:)=candidate%fields(peer)%values
        packed%fields(peer)%values(nvar+1,:,:,:)=temps%fields(peer)%values(1,:,:,:)
      end do
    end subroutine
    subroutine recover_from(depth,success)
      integer, intent(in) :: depth
      logical, intent(out) :: success
      integer :: level
      do level=depth,tree%level_count()-1
        call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,candidate,temps,comm,success,depth=level)
        if (.not.success) return
      end do
      success=.true.
    end subroutine
    subroutine root_budget(flux)
      type(amr_tree_flux_3d), intent(in) :: flux
      integer :: i,j,k,dims(3)
      dims=tree%nodes(1)%shape
      if (faces(1)/='periodic') then
        do k=1,dims(3)
          do j=1,dims(2)
            budget=budget+(flux%x(:,1,j,k)-flux%x(:,dims(1)+1,j,k))*spacing(2)*spacing(3)
          end do
        end do
      end if
      if (faces(3)/='periodic') then
        do k=1,dims(3)
          do i=1,dims(1)
            budget=budget+(flux%y(:,i,1,k)-flux%y(:,i,dims(2)+1,k))*spacing(1)*spacing(3)
          end do
        end do
      end if
      if (faces(5)/='periodic') then
        do j=1,dims(2)
          do i=1,dims(1)
            budget=budget+(flux%z(:,i,j,1)-flux%z(:,i,j,dims(3)+1))*spacing(1)*spacing(2)
          end do
        end do
      end if
    end subroutine
  end subroutine

  subroutine layout(tree,owners,species,state,temperature,spacing,origin,domain,ratios,counts,comm,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: spacing(3)
    integer(int64), allocatable, intent(out) :: origin(:,:),domain(:,:),counts(:)
    integer, allocatable, intent(out) :: ratios(:)
    logical, intent(out) :: ok
    integer :: level,node,rank,status,ios
    real(dp) :: low(3),high(3)
    logical :: valid
    call validate_mpi_amr_tree_local_3d(tree,owners,state,reactive_nvar(size(species)),comm,ok)
    if (.not.ok) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,ok)
    if (.not.ok) return
    ok=.false.
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    valid=all(ieee_is_finite(spacing))
    if (valid) valid=all(spacing>0)
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      if (.not.all(ieee_is_finite(state%fields(node)%values))) valid=.false.
      if (.not.all(ieee_is_finite(temperature%fields(node)%values))) then
        valid=.false.
      else if (any(temperature%fields(node)%values<=0)) then
        valid=.false.
      end if
    end do
    call agree(comm,valid)
    if (.not.valid) return
    call MPI_Allreduce(spacing,low,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(spacing,high,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(low/=high)) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,ok)
    call agree(comm,ok)
    if (.not.ok) return
    ok=.false.
    allocate(counts(size(ratios)),stat=ios)
    valid=ios==0;call agree(comm,valid)
    if (.not.valid) return
    counts(1)=1
    do level=2,size(ratios)
      if (counts(level-1)>huge(1_int64)/int(ratios(level),int64)) return
      counts(level)=counts(level-1)*int(ratios(level),int64)
    end do
    ok=.true.
  end subroutine
  subroutine agree(comm,condition)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(inout) :: condition
    logical :: result
    integer :: status
    call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,status)
    condition=.false.
    if (status==MPI_SUCCESS) condition=result
  end subroutine
end module mpi_amr_patch_tree_subcycle_3d_mod
