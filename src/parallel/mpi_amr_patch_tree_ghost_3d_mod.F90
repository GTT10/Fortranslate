module mpi_amr_patch_tree_ghost_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use rotating_frame_3d_mod, only: rotating_frame_3d
  use nasa7_thermo_mod, only: nasa7_species
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use amr_patch_tree_ghost_3d_mod, only: sample_amr_patch_tree_ghosts_3d,amr_tree_cell_source_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: fetch_mpi_amr_tree_cells_3d,validate_mpi_amr_tree_local_3d
  implicit none
  private
  public :: fill_mpi_amr_tree_ghosts_3d
  type, extends(amr_tree_cell_source_3d) :: cached_tree_source
    integer, allocatable :: requests(:,:)
    real(dp), allocatable :: values(:,:)
    logical, allocatable :: future(:)
    integer :: pass=0,nvar=0
    integer(int64) :: count=0,expected=0
  contains
    procedure :: read=>read_cached_cell
  end type
contains
  ! Collective ghost fill. Only owned selected-depth output fields allocate.
  subroutine fill_mpi_amr_tree_ghosts_3d(tree,owners,local,nvar,faces,ghosts,comm,ok, &
      viscous_wall,inflow_state,wall_velocity,coarse_limiter,species_count,coarse_future,coarse_fractions,depth,ghost_width, &
      species,inflow_mach,inflow_temperature,inflow_pressure,limit_internal_energy,stored_internal_energy, &
      frame,domain_lower,root_spacing)
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3),root_spacing(3)
    logical, optional, intent(in) :: limit_internal_energy,stored_internal_energy
    type(nasa7_species), optional, intent(in) :: species(:)
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
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
    type(amr_tree_state_3d), optional, intent(in) :: coarse_future
    real(dp), optional, intent(in) :: coarse_fractions(:)
    integer, optional, intent(in) :: depth
    integer, optional, intent(in) :: ghost_width
    type(amr_tree_state_3d) :: candidate
    type(cached_tree_source) :: reader
    real(dp), allocatable :: scratch(:,:,:,:),fractions(:),low_f(:),high_f(:),endpoint_values(:,:)
    integer, allocatable :: queries(:,:)
    integer :: controls(8),low(8),high(8),active,levels,query,slot,endpoint,width
    real(dp) :: inlet_controls(18),inlet_low(18),inlet_high(18)
    real(dp) :: frame_controls(11),frame_low(11),frame_high(11)
    logical :: use_future,include_internal,stored_internal
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
    active=-1
    if (present(depth)) active=depth
    width=2
    if(present(ghost_width)) width=ghost_width
    stored_internal=.false.
    if(present(stored_internal_energy)) stored_internal=stored_internal_energy
    include_internal=.false.
    if(present(limit_internal_energy)) include_internal=limit_internal_energy
    controls=[active,merge(1,0,present(coarse_future)),merge(1,0,present(coarse_fractions)),width, &
      merge(1,0,include_internal),merge(1,0,stored_internal),merge(1,0,present(frame)), &
      merge(1,0,present(domain_lower).and.present(root_spacing))]
    call MPI_Allreduce(controls,low,8,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,8,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(low/=high).or.controls(2)/=controls(3)) return
    valid=controls(7)==0.or.controls(8)==1
    if(present(frame)) valid=valid.and.frame%is_valid()
    if(present(domain_lower)) valid=valid.and.all(ieee_is_finite(domain_lower))
    if(present(root_spacing)) valid=valid.and.all(ieee_is_finite(root_spacing)).and.all(root_spacing>0)
    call agree(valid,ready)
    if(.not.ready) return
    frame_controls=0
    if(present(frame)) then
      frame_controls(1)=real(frame%axis,dp)
      frame_controls(2)=frame%angular_speed
      frame_controls(3:5)=frame%axis_point
      frame_controls(6:8)=domain_lower
      frame_controls(9:11)=root_spacing
    endif
    call MPI_Allreduce(frame_controls,frame_low,11,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(frame_controls,frame_high,11,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(any(frame_low/=frame_high)) return
    if(width<1.or.width>128) return
    if (active< -1.or.active>=tree%level_count()) return
    call validate_mpi_amr_tree_local_3d(tree,owners,local,nvar,comm,ready)
    if (.not.ready) return
    levels=tree%level_count()-1
    if (active>=0) levels=active
    use_future=present(coarse_future)
    if (use_future) then
      valid=size(coarse_fractions)==levels
      if (valid) valid=all(ieee_is_finite(coarse_fractions))
      if (valid) valid=all(coarse_fractions>=0).and.all(coarse_fractions<=1)
      call agree(valid,ready)
      if (.not.ready) return
      allocate(fractions(levels),low_f(levels),high_f(levels),stat=ios)
      call agree(ios==0,ready)
      if (.not.ready) return
      fractions=coarse_fractions
      call MPI_Allreduce(fractions,low_f,levels,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if (status/=MPI_SUCCESS) return
      call MPI_Allreduce(fractions,high_f,levels,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if (status/=MPI_SUCCESS) return
      if (any(low_f/=high_f)) return
      call validate_mpi_amr_tree_local_3d(tree,owners,coarse_future,nvar,comm,ready,maximum_depth=levels-1)
      if (.not.ready) return
      valid=.true.
      do node=1,size(tree%nodes)
        if (owners(node)/=rank.or.tree%nodes(node)%level>=levels) cycle
        if (.not.all(ieee_is_finite(coarse_future%fields(node)%values))) valid=.false.
      end do
      call agree(valid,ready)
      if (.not.ready) return
    end if
    ! Validate boundary options on empty ranks too; planning ignores payload arithmetic.
    reader%pass=0;reader%count=0;reader%nvar=nvar
    call sample_node(1,scratch,valid,planning=.true.,validate_only=.true.)
    call agree(valid,ready)
    if (.not.ready) return
    inlet_controls=0
    inlet_controls(:6)=-1
    if(present(inflow_mach)) inlet_controls(:6)=inflow_mach
    if(present(inflow_pressure)) inlet_controls(13:18)=inflow_pressure
    if(present(inflow_temperature)) then
      where(inlet_controls(:6)>=0) inlet_controls(7:12)=inflow_temperature
    end if
    call MPI_Allreduce(inlet_controls,inlet_low,18,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(inlet_controls,inlet_high,18,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(any(inlet_low/=inlet_high)) return
    if (allocated(scratch)) deallocate(scratch)
    allocate(candidate%fields(size(tree%nodes)),stat=ios)
    call agree(ios==0,ready)
    if (.not.ready) return
    ! Two metadata walks count then record the exact donor sequence, including PLM stencils.
    do pass=1,2
      reader%pass=pass;reader%count=0;valid=.true.
      do node=1,size(tree%nodes)
        if (owners(node)/=rank) cycle
        if (active>=0.and.tree%nodes(node)%level/=active) cycle
        call sample_node(node,scratch,valid,planning=.true.)
        if (.not.valid) exit
      end do
      call agree(valid,ready)
      if (.not.ready) return
      if (pass==1) then
        reader%expected=reader%count
        allocate(reader%requests(4,int(reader%count)),reader%future(int(reader%count)),stat=ios)
        call agree(ios==0,ready)
        if (.not.ready) return
      else
        call agree(reader%count==reader%expected,ready)
        if (.not.ready) return
      end if
    end do
    allocate(reader%values(nvar,int(reader%expected)),stat=ios)
    call agree(ios==0,ready)
    if (.not.ready) return
    do endpoint=0,merge(1,0,use_future)
      allocate(queries(4,count(reader%future.eqv.(endpoint==1))),stat=ios)
      call agree(ios==0,ready)
      if (.not.ready) return
      slot=0
      do query=1,int(reader%expected)
        if (reader%future(query).neqv.(endpoint==1)) cycle
        slot=slot+1;queries(:,slot)=reader%requests(:,query)
      end do
      if (endpoint==0) then
        call fetch_mpi_amr_tree_cells_3d(tree,owners,local,nvar,queries,endpoint_values,comm,valid)
      else
        call fetch_mpi_amr_tree_cells_3d(tree,owners,coarse_future,nvar,queries,endpoint_values,comm,valid, &
          maximum_depth=levels-1)
      end if
      if (.not.valid) return
      slot=0
      do query=1,int(reader%expected)
        if (reader%future(query).neqv.(endpoint==1)) cycle
        slot=slot+1;reader%values(:,query)=endpoint_values(:,slot)
      end do
      deallocate(queries,endpoint_values)
    end do
    reader%pass=3;reader%count=0;valid=.true.
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      if (active>=0.and.tree%nodes(node)%level/=active) cycle
      call sample_node(node,candidate%fields(node)%values,valid)
      if (.not.valid) exit
    end do
    call agree(valid.and.reader%count==reader%expected,ready)
    if (.not.ready) return
    call move_alloc(candidate%fields,ghosts%fields)
    ok=.true.
  contains
    subroutine sample_node(peer,output,success,planning,validate_only)
      integer, intent(in) :: peer
      real(dp), allocatable, intent(inout) :: output(:,:,:,:)
      logical, intent(out) :: success
      logical, optional, intent(in) :: planning,validate_only
      if (use_future) then
        call sample_amr_patch_tree_ghosts_3d(tree,nvar,peer,faces,reader,output,success, &
          viscous_wall,inflow_state,wall_velocity,coarse_fractions=fractions(:tree%nodes(peer)%level), &
          coarse_limiter=coarse_limiter,species_count=species_count,planning=planning,validate_only=validate_only, &
          ghost_width=width,species=species,inflow_mach=inflow_mach, &
          inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
          limit_internal_energy=limit_internal_energy,stored_internal_energy=stored_internal_energy, &
          frame=frame,domain_lower=domain_lower,root_spacing=root_spacing)
      else
        call sample_amr_patch_tree_ghosts_3d(tree,nvar,peer,faces,reader,output,success, &
          viscous_wall,inflow_state,wall_velocity,coarse_limiter=coarse_limiter, &
          species_count=species_count,planning=planning,validate_only=validate_only,ghost_width=width, &
          species=species,inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
          limit_internal_energy=limit_internal_energy,stored_internal_energy=stored_internal_energy, &
          frame=frame,domain_lower=domain_lower,root_spacing=root_spacing)
      end if
    end subroutine
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
      self%future(slot)=future
    else if (self%pass==3) then
      if (self%count>self%expected) return
      if (any(self%requests(:,slot)/=[node,index])) return
      if (self%future(slot).neqv.future) return
      value=self%values(:,slot)
    end if
    ok=.true.
  end subroutine
end module
