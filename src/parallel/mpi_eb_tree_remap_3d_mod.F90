module mpi_eb_tree_remap_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use rotating_frame_3d_mod, only: rotating_frame_3d
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim
  use eb_reactive_transfer_3d_mod, only: recover_reactive_eb_cell_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_two_level_3d_mod, only: matching_eb_metrics_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_metadata_3d,fetch_mpi_amr_tree_cells_3d
  use mpi_eb_tree_transfer_3d_mod, only: validate_mpi_eb_tree_geometry_3d,average_down_mpi_eb_tree_3d, &
    average_down_mpi_eb_tree_field_3d
  implicit none
  private
  public :: remap_mpi_eb_tree_3d
contains
  ! Shared root/refinement lattice and immutable geometry, matching serial PCM.
  ! Geometry and owner-local U/T commit together; caller commits topology/owners.
  subroutine remap_mpi_eb_tree_3d(species,old_tree,new_tree,old_owners,new_owners,finest,geometry,state,temperature,comm,ok, &
      reaction_source,frame,heat_release,internal_energy)
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_patch_tree_3d), intent(in) :: old_tree,new_tree
    integer, intent(in) :: old_owners(:),new_owners(:)
    type(eb_geometry_3d), intent(in) :: finest
    type(eb_geometry_3d), allocatable, intent(inout) :: geometry(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source,heat_release,internal_energy
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: snapshot,st,candidate,ct,old_source,new_source,old_heat,new_heat,old_held,new_held
    type(eb_geometry_3d), allocatable :: new_geometry(:)
    integer(int64), allocatable :: old_origin(:,:),old_domain(:,:),new_origin(:,:),new_domain(:,:)
    integer, allocatable :: old_ratios(:),new_ratios(:),requests(:,:),destination(:,:)
    real(dp), allocatable :: values(:,:),guesses(:,:),source_values(:,:),heat_values(:,:),held_values(:,:)
    real(dp) :: primitive(reactive_nprim(size(species))),recovered,speed
    integer(int64) :: cells,total,global(3),mapped(3),local(3)
    integer :: nv,rank,status,node,donor,level,depth,n(3),i,j,k,source(3),q,ios,axis,old_levels,new_levels
    integer :: supplied,minimum_supplied,maximum_supplied
    integer :: frame_axis,axis_min,axis_max
    real(dp) :: frame_values(4),frame_min(4),frame_max(4)
    logical :: valid,ready,found
    ok=.false.;nv=reactive_nvar(size(species))
    supplied=merge(1,0,present(internal_energy))
    call MPI_Allreduce(supplied,minimum_supplied,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(supplied,maximum_supplied,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.minimum_supplied/=maximum_supplied) return
    supplied=merge(1,0,present(reaction_source))
    call MPI_Allreduce(supplied,minimum_supplied,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(supplied,maximum_supplied,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.minimum_supplied/=maximum_supplied) return
    supplied=merge(1,0,present(heat_release))
    call MPI_Allreduce(supplied,minimum_supplied,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(supplied,maximum_supplied,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.minimum_supplied/=maximum_supplied) return
    supplied=merge(1,0,present(frame))
    call MPI_Allreduce(supplied,minimum_supplied,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(supplied,maximum_supplied,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.minimum_supplied/=maximum_supplied) return
    if(present(frame)) then
      call agree(frame%is_valid(),ready)
      if(.not.ready) return
      frame_axis=frame%axis
      frame_values=[frame%angular_speed,frame%axis_point]
      call MPI_Allreduce(frame_axis,axis_min,1,MPI_INTEGER,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(frame_axis,axis_max,1,MPI_INTEGER,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(frame_values,frame_min,4,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(frame_values,frame_max,4,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.axis_min/=axis_max.or.any(frame_min/=frame_max)) return
    end if
    call agree(allocated(geometry),ready)
    if (.not.ready) return
    snapshot=state;st=temperature
    if(present(internal_energy)) then
      old_held=internal_energy
      call average_down_mpi_eb_tree_3d(species,old_tree,old_owners,geometry,snapshot,st,comm,ready,frame,old_held)
    else
      call average_down_mpi_eb_tree_3d(species,old_tree,old_owners,geometry,snapshot,st,comm,ready,frame)
    endif
    if (.not.ready) return
    if(present(reaction_source)) then
      old_source=reaction_source
      call average_down_mpi_eb_tree_field_3d(old_tree,old_owners,geometry,old_source,nv-4,comm,ready)
      if(.not.ready) return
    end if
    if(present(heat_release)) then
      old_heat=heat_release
      call average_down_mpi_eb_tree_field_3d(old_tree,old_owners,geometry,old_heat,1,comm,ready)
      if(.not.ready) return
    end if
    call validate_mpi_amr_tree_metadata_3d(new_tree,new_owners,nv,comm,ready)
    if (.not.ready) return
    call amr_tree_integer_layout_3d(old_tree,old_origin,old_domain,old_ratios,valid)
    call agree(valid,ready)
    if (.not.ready) return
    call amr_tree_integer_layout_3d(new_tree,new_origin,new_domain,new_ratios,valid)
    call agree(valid,ready)
    if (.not.ready) return
    old_levels=size(old_ratios);new_levels=size(new_ratios)
    call build_geometry(valid)
    call agree(valid,ready)
    if (.not.ready) return
    call validate_mpi_eb_tree_geometry_3d(new_tree,new_owners,new_geometry,nv,comm,ready)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    total=0_int64;valid=.true.
    do node=1,size(new_owners)
      if (new_owners(node)/=rank) cycle
      cells=1_int64
      do axis=1,3
        if (cells>int(huge(1)/max(nv,4),int64)/int(new_tree%nodes(node)%shape(axis),int64)) then
          valid=.false.;exit
        end if
        cells=cells*int(new_tree%nodes(node)%shape(axis),int64)
      end do
      if (.not.valid) exit
      if (cells>int(huge(1)/max(nv,4),int64)-total) then
        valid=.false.;exit
      end if
      total=total+cells
    end do
    call agree(valid,ready)
    if (.not.ready) return
    allocate(requests(4,int(total)),destination(4,int(total)),candidate%fields(size(new_owners)), &
      ct%fields(size(new_owners)),stat=ios)
    valid=ios==0
    if (valid) then
      do node=1,size(new_owners)
        if (new_owners(node)/=rank) cycle
        n=new_tree%nodes(node)%shape
        allocate(candidate%fields(node)%values(nv,n(1),n(2),n(3)), &
          ct%fields(node)%values(1,n(1),n(2),n(3)),stat=ios)
        if (ios/=0) then
          valid=.false.;exit
        end if
      end do
    end if
    call agree(valid,ready)
    if (.not.ready) return
    if(present(reaction_source)) then
      allocate(new_source%fields(size(new_owners)),stat=ios)
      valid=ios==0
      if(valid) then
        do node=1,size(new_owners)
          if(new_owners(node)/=rank) cycle
          n=new_tree%nodes(node)%shape
          allocate(new_source%fields(node)%values(nv-4,n(1),n(2),n(3)),source=0._dp,stat=ios)
          if(ios/=0) then
            valid=.false.;exit
          end if
        end do
      end if
      call agree(valid,ready)
      if(.not.ready) return
    end if
    if(present(heat_release)) then
      allocate(new_heat%fields(size(new_owners)),stat=ios)
      valid=ios==0
      if(valid) then
        do node=1,size(new_owners)
          if(new_owners(node)/=rank) cycle
          n=new_tree%nodes(node)%shape
          allocate(new_heat%fields(node)%values(1,n(1),n(2),n(3)),source=0._dp,stat=ios)
          if(ios/=0) then
            valid=.false.;exit
          end if
        end do
      end if
      call agree(valid,ready)
      if(.not.ready) return
    end if
    if(present(internal_energy)) then
      allocate(new_held%fields(size(new_owners)),stat=ios)
      valid=ios==0
      if(valid) then
        do node=1,size(new_owners)
          if(new_owners(node)/=rank) cycle
          n=new_tree%nodes(node)%shape
          allocate(new_held%fields(node)%values(1,n(1),n(2),n(3)),source=0._dp,stat=ios)
          if(ios/=0) then
            valid=.false.;exit
          endif
        enddo
      endif
      call agree(valid,ready)
      if(.not.ready) return
    endif
    q=0
    do node=1,size(new_owners)
      if (new_owners(node)/=rank) cycle
      n=new_tree%nodes(node)%shape;level=new_tree%nodes(node)%level+1
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        global=new_origin(:,node)+int([i,j,k]-1,int64);found=.false.
        do depth=min(level,old_levels),1,-1
          mapped=global/(new_domain(:,level)/old_domain(:,depth))
          do donor=1,size(old_owners)
            if (old_tree%nodes(donor)%level+1/=depth) cycle
            local=mapped-old_origin(:,donor)+1_int64
            if (any(local<1_int64).or.any(local>int(old_tree%nodes(donor)%shape,int64))) cycle
            source=int(local)
            if (new_geometry(node)%volume_fraction(i,j,k)>0._dp) then
              if (geometry(donor)%volume_fraction(source(1),source(2),source(3))<=0._dp) valid=.false.
            end if
            q=q+1;requests(:,q)=[donor,source];destination(:,q)=[node,i,j,k]
            found=.true.;exit
          end do
          if (found) exit
        end do
        if (.not.found) valid=.false.
      end do
      end do
      end do
    end do
    call agree(valid.and.q==int(total),ready)
    if (.not.ready) return
    ! Fetch finite payloads including inactive cells, preserving exact overlap.
    call fetch_mpi_amr_tree_cells_3d(old_tree,old_owners,snapshot,nv,requests,values,comm,ready)
    if (.not.ready) return
    call fetch_mpi_amr_tree_cells_3d(old_tree,old_owners,st,1,requests,guesses,comm,ready)
    if (.not.ready) return
    if(present(reaction_source)) then
      call fetch_mpi_amr_tree_cells_3d(old_tree,old_owners,old_source,nv-4,requests,source_values,comm,ready)
      if(.not.ready) return
    end if
    if(present(heat_release)) then
      call fetch_mpi_amr_tree_cells_3d(old_tree,old_owners,old_heat,1,requests,heat_values,comm,ready)
      if(.not.ready) return
    end if
    if(present(internal_energy)) then
      call fetch_mpi_amr_tree_cells_3d(old_tree,old_owners,old_held,1,requests,held_values,comm,ready)
      if(.not.ready) return
    endif
    do q=1,int(total)
      node=destination(1,q);i=destination(2,q);j=destination(3,q);k=destination(4,q)
      if (new_geometry(node)%volume_fraction(i,j,k)>0._dp) then
        if(present(reaction_source)) new_source%fields(node)%values(:,i,j,k)=source_values(:,q)
        if(present(heat_release)) new_heat%fields(node)%values(:,i,j,k)=heat_values(:,q)
        if(present(internal_energy)) then
        call recover_reactive_eb_cell_3d(species,new_geometry(node),[i,j,k],values(:,q),guesses(1,q), &
          primitive,recovered,speed,ready,frame,held_values(1,q))
        else
        call recover_reactive_eb_cell_3d(species,new_geometry(node),[i,j,k],values(:,q),guesses(1,q), &
          primitive,recovered,speed,ready,frame)
        endif
        if (.not.ready) then
          valid=.false.;cycle
        end if
      end if
      candidate%fields(node)%values(:,i,j,k)=values(:,q)
      ct%fields(node)%values(1,i,j,k)=guesses(1,q)
      if(present(internal_energy)) then
        new_held%fields(node)%values(1,i,j,k)=held_values(1,q)
        if(new_geometry(node)%volume_fraction(i,j,k)>0) ct%fields(node)%values(1,i,j,k)=recovered
      endif
      if(present(frame).and.new_geometry(node)%volume_fraction(i,j,k)>0._dp.and. &
        new_tree%nodes(node)%level/=old_tree%nodes(requests(1,q))%level) &
        ct%fields(node)%values(1,i,j,k)=recovered
    end do
    call agree(valid,ready)
    if (.not.ready) return
    if(present(internal_energy)) then
      call average_down_mpi_eb_tree_3d(species,new_tree,new_owners,new_geometry,candidate,ct,comm,ready,frame,new_held)
    else
      call average_down_mpi_eb_tree_3d(species,new_tree,new_owners,new_geometry,candidate,ct,comm,ready,frame)
    endif
    if (.not.ready) return
    if(present(reaction_source)) then
      call average_down_mpi_eb_tree_field_3d(new_tree,new_owners,new_geometry,new_source,nv-4,comm,ready)
      if(.not.ready) return
    end if
    if(present(heat_release)) then
      call average_down_mpi_eb_tree_field_3d(new_tree,new_owners,new_geometry,new_heat,1,comm,ready)
      if(.not.ready) return
    end if
    if(present(reaction_source)) call move_alloc(new_source%fields,reaction_source%fields)
    if(present(heat_release)) call move_alloc(new_heat%fields,heat_release%fields)
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(ct%fields,temperature%fields)
    call move_alloc(new_geometry,geometry)
    if(present(internal_energy)) call move_alloc(new_held%fields,internal_energy%fields)
    ok=.true.
  contains
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) result=.false.
    end subroutine
    subroutine build_geometry(good)
      logical, intent(out) :: good
      type(eb_geometry_3d), allocatable :: levels(:)
      type(eb_geometry_3d) :: section
      integer, allocatable :: ratios(:)
      integer :: depth,p,lo(3),hi(3),level_count,allocation_status
      logical :: success
      good=.false.
      if (.not.finest%is_valid()) return
      if (any(old_domain(:,1)/=new_domain(:,1))) return
      depth=min(old_levels,new_levels)
      if (any(old_ratios(:depth)/=new_ratios(:depth))) return
      level_count=max(old_levels,new_levels)
      if (old_levels>=new_levels) then
        ratios=old_ratios
        if (any(int([finest%nx,finest%ny,finest%nz],int64)/=old_domain(:,level_count))) return
      else
        ratios=new_ratios
        if (any(int([finest%nx,finest%ny,finest%nz],int64)/=new_domain(:,level_count))) return
      end if
      allocate(levels(level_count),new_geometry(size(new_owners)),stat=allocation_status)
      if (allocation_status/=0) return
      levels(level_count)=finest
      do depth=level_count,2,-1
        call coarsen_eb_geometry_3d(levels(depth),ratios(depth),levels(depth-1),success)
        if (.not.success) return
      end do
      do p=1,size(old_owners)
        depth=old_tree%nodes(p)%level+1
        lo=int(old_origin(:,p))+1;hi=int(old_origin(:,p))+old_tree%nodes(p)%shape
        call extract_eb_geometry_3d(levels(depth),lo,hi,section,success)
        if (.not.success) return
        if (.not.matching_eb_metrics_3d(section,geometry(p))) return
      end do
      do p=1,size(new_owners)
        depth=new_tree%nodes(p)%level+1
        lo=int(new_origin(:,p))+1;hi=int(new_origin(:,p))+new_tree%nodes(p)%shape
        call extract_eb_geometry_3d(levels(depth),lo,hi,new_geometry(p),success)
        if (.not.success) return
      end do
      good=.true.
    end subroutine
  end subroutine
end module mpi_eb_tree_remap_3d_mod
