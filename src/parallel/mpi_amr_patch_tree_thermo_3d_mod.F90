module mpi_amr_patch_tree_thermo_3d_mod
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use reactive_3d_mod, only: recover_reactive_temperatures_3d,compute_reactive_cfl_timestep_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  implicit none
  private
  public :: recover_mpi_amr_tree_temperatures_3d,compute_mpi_amr_tree_cfl_3d
contains
  subroutine recover_mpi_amr_tree_temperatures_3d(tree,owners,species,state,temperature,comm,ok,depth)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(in) :: state
    type(amr_tree_state_3d), intent(inout) :: temperature
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer, optional, intent(in) :: depth
    type(amr_tree_state_3d) :: candidate
    integer :: rank,status,node,n(3),ios,selected,low,high
    logical :: valid,ready
    ok=.false.;selected=-1
    if (present(depth)) selected=depth
    call MPI_Allreduce(selected,low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(selected,high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (low/=high) return
    call validate_inputs(tree,owners,species,state,temperature,comm,rank,ready)
    if (.not.ready) return
    if (selected< -1.or.selected>=tree%level_count()) return
    allocate(candidate%fields(size(tree%nodes)),stat=ios)
    valid=ios==0
    if (valid) then
      do node=1,size(tree%nodes)
        if (owners(node)/=rank) cycle
        if (selected>=0.and.tree%nodes(node)%level/=selected) cycle
        n=tree%nodes(node)%shape
        allocate(candidate%fields(node)%values(1,n(1),n(2),n(3)),stat=ios)
        if (ios/=0) then
          valid=.false.;exit
        end if
        call recover_reactive_temperatures_3d(species,state%fields(node)%values,temperature%fields(node)%values(1,:,:,:), &
          n(1),n(2),n(3),candidate%fields(node)%values(1,:,:,:),valid)
        if (.not.valid) exit
      end do
    end if
    call agree(comm,valid,ready)
    if (.not.ready) return
    do node=1,size(tree%nodes)
      if (.not.allocated(candidate%fields(node)%values)) cycle
      call move_alloc(candidate%fields(node)%values,temperature%fields(node)%values)
    end do
    ok=.true.
  end subroutine

  subroutine compute_mpi_amr_tree_cfl_3d(tree,owners,species,state,temperature,root_spacing,cfl,dt,comm,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: root_spacing(3),cfl
    real(dp), intent(inout) :: dt
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    real(dp), allocatable :: h(:,:)
    real(dp) :: controls(4),low(4),high(4),local_dt,node_dt,global_dt
    integer :: rank,status,node,parent,n(3),ios
    logical :: valid,ready
    ok=.false.
    call validate_inputs(tree,owners,species,state,temperature,comm,rank,ready)
    if (.not.ready) return
    controls=[root_spacing,cfl]
    valid=all(ieee_is_finite(controls))
    if (valid) valid=all(controls>0)
    call agree(comm,valid,ready)
    if (.not.ready) return
    call MPI_Allreduce(controls,low,4,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,4,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(low/=high)) return
    allocate(h(3,size(tree%nodes)),stat=ios)
    call agree(comm,ios==0,ready)
    if (.not.ready) return
    h(:,1)=root_spacing
    do node=2,size(tree%nodes)
      parent=tree%nodes(node)%parent
      h(:,node)=h(:,parent)/real(tree%nodes(node)%patch%refinement_ratio,dp)
      if (any(h(:,node)<=0)) return
    end do
    local_dt=huge(1.0_dp);valid=.true.
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      call compute_reactive_cfl_timestep_3d(species,state%fields(node)%values,temperature%fields(node)%values(1,:,:,:), &
        n(1),n(2),n(3),h(1,node),h(2,node),h(3,node),cfl,node_dt,valid)
      if (.not.valid) exit
      local_dt=min(local_dt,node_dt)
    end do
    call agree(comm,valid,ready)
    if (.not.ready) return
    call MPI_Allreduce(local_dt,global_dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    dt=global_dt;ok=.true.
  end subroutine

  subroutine validate_inputs(tree,owners,species,state,temperature,comm,rank,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    type(MPI_Comm), intent(in) :: comm
    integer, intent(out) :: rank
    logical, intent(out) :: ok
    integer :: status,node
    logical :: valid
    call validate_mpi_amr_tree_local_3d(tree,owners,state,reactive_nvar(size(species)),comm,ok)
    if (.not.ok) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,ok)
    if (.not.ok) return
    ok=.false.
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    valid=size(species)>0
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      if (.not.all(ieee_is_finite(state%fields(node)%values))) valid=.false.
      if (.not.all(ieee_is_finite(temperature%fields(node)%values))) valid=.false.
      if (any(temperature%fields(node)%values<=0)) valid=.false.
    end do
    call agree(comm,valid,ok)
  end subroutine
  subroutine agree(comm,condition,result)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(in) :: condition
    logical, intent(out) :: result
    integer :: ierr
    call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS) result=.false.
  end subroutine
end module
