module mpi_amr_patch_tree_ghost_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use amr_patch_tree_ghost_3d_mod, only: sample_amr_patch_tree_ghosts_3d,amr_tree_cell_source_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: fetch_mpi_amr_tree_cells_3d
  implicit none
  private
  public :: fill_mpi_amr_tree_ghosts_3d
  type, extends(amr_tree_cell_source_3d) :: cached_tree_source
    integer, allocatable :: requests(:,:)
    real(dp), allocatable :: values(:,:)
    integer :: pass=0,nvar=0
    integer(int64) :: count=0,expected=0
  contains
    procedure :: read=>read_cached_cell
  end type
contains
  ! Collective synchronous ghost fill. Only locally owned output fields allocate.
  subroutine fill_mpi_amr_tree_ghosts_3d(tree,owners,local,nvar,faces,ghosts,comm,ok, &
      viscous_wall,inflow_state,wall_velocity,coarse_limiter,species_count)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar
    type(amr_tree_state_3d), intent(in) :: local
    character(len=*), intent(in) :: faces(6)
    type(amr_tree_state_3d), intent(inout) :: ghosts
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    logical, optional, intent(in) :: viscous_wall
    real(dp), optional, intent(in) :: inflow_state(:,:),wall_velocity(:,:)
    character(len=*), optional, intent(in) :: coarse_limiter
    integer, optional, intent(in) :: species_count
    type(amr_tree_state_3d) :: candidate
    type(cached_tree_source) :: reader
    real(dp), allocatable :: scratch(:,:,:,:)
    integer :: rank,nranks,status,ios,node,pass
    logical :: valid,ready
    ok=.false.
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,status)
    if (status/=MPI_SUCCESS) return
    valid=tree%is_valid().and.nvar>=5
    if (valid) valid=size(owners)==size(tree%nodes)
    if (valid) valid=all(owners>=0).and.all(owners<nranks)
    call agree(valid,ready)
    if (.not.ready) return
    ! Validate boundary options on empty ranks too; planning ignores payload arithmetic.
    reader%pass=0;reader%count=0;reader%nvar=nvar
    call sample_amr_patch_tree_ghosts_3d(tree,nvar,1,faces,reader,scratch,valid, &
      viscous_wall,inflow_state,wall_velocity,coarse_limiter=coarse_limiter,species_count=species_count, &
      planning=.true.,validate_only=.true.)
    call agree(valid,ready)
    if (.not.ready) return
    if (allocated(scratch)) deallocate(scratch)
    allocate(candidate%fields(size(tree%nodes)),stat=ios)
    call agree(ios==0,ready)
    if (.not.ready) return
    ! Two metadata walks count then record the exact donor sequence, including PLM stencils.
    do pass=1,2
      reader%pass=pass;reader%count=0;valid=.true.
      do node=1,size(tree%nodes)
        if (owners(node)/=rank) cycle
        call sample_amr_patch_tree_ghosts_3d(tree,nvar,node,faces,reader,scratch,valid, &
          viscous_wall,inflow_state,wall_velocity,coarse_limiter=coarse_limiter,species_count=species_count,planning=.true.)
        if (.not.valid) exit
      end do
      call agree(valid,ready)
      if (.not.ready) return
      if (pass==1) then
        reader%expected=reader%count
        allocate(reader%requests(4,int(reader%count)),stat=ios)
        call agree(ios==0,ready)
        if (.not.ready) return
      else
        call agree(reader%count==reader%expected,ready)
        if (.not.ready) return
      end if
    end do
    call fetch_mpi_amr_tree_cells_3d(tree,owners,local,nvar,reader%requests,reader%values,comm,valid)
    if (.not.valid) return
    reader%pass=3;reader%count=0;valid=.true.
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      call sample_amr_patch_tree_ghosts_3d(tree,nvar,node,faces,reader,candidate%fields(node)%values,valid, &
        viscous_wall,inflow_state,wall_velocity,coarse_limiter=coarse_limiter,species_count=species_count)
      if (.not.valid) exit
    end do
    call agree(valid.and.reader%count==reader%expected,ready)
    if (.not.ready) return
    call move_alloc(candidate%fields,ghosts%fields)
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
  subroutine read_cached_cell(self,node,index,future,value,ok)
    class(cached_tree_source), intent(inout) :: self
    integer, intent(in) :: node,index(3)
    logical, intent(in) :: future
    real(dp), intent(out) :: value(:)
    logical, intent(out) :: ok
    integer :: slot
    ok=.false.
    if (future) return
    value=1.0_dp
    if (self%pass==0) then
      ok=.true.;return
    end if
    self%count=self%count+1
    if (self%count>int(huge(1)/max(4,self%nvar),int64)) return
    slot=int(self%count)
    if (self%pass==2) then
      if (self%count>self%expected) return
      self%requests(:,slot)=[node,index]
    else if (self%pass==3) then
      if (self%count>self%expected) return
      if (any(self%requests(:,slot)/=[node,index])) return
      value=self%values(:,slot)
    end if
    ok=.true.
  end subroutine
end module
