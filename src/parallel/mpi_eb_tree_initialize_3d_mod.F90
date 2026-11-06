module mpi_eb_tree_initialize_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use dual_energy_mod, only: initialize_nasa7_internal_energy_patch
  use nasa7_thermo_mod, only: nasa7_species
  use rotating_frame_3d_mod, only: rotating_frame_3d
  use reactive_1d_mod, only: reactive_nvar
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_transfer_3d_mod, only: prolong_reactive_eb_state_pcm_3d
  use mpi_eb_tree_transfer_3d_mod, only: validate_mpi_eb_tree_geometry_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_metadata_3d,fetch_mpi_amr_tree_cells_3d
  implicit none
  private
  public :: initialize_mpi_eb_tree_3d
contains
  subroutine initialize_mpi_eb_tree_3d(species,tree,owners,finest,root_state,root_temperature, &
      geometry,state,temperature,comm,ok,frame,internal_energy)
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: finest
    ! Only owners(1) supplies a root seed; other ranks may pass zero-sized arrays.
    real(dp), intent(in) :: root_state(:,:,:,:),root_temperature(:,:,:)
    type(eb_geometry_3d), allocatable, intent(inout) :: geometry(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    type(eb_geometry_3d), allocatable :: levels(:),shapes(:)
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy
    type(amr_tree_state_3d) :: candidate,ct,ch
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:),requests(:,:)
    real(dp), allocatable :: values(:,:),guesses(:,:),held_values(:,:),coarse_u(:,:,:,:),coarse_t(:,:,:)
    real(dp), allocatable :: root_held(:,:,:),coarse_held(:,:,:)
    integer :: nv,rank,status,count,node,level,d(3),lo(3),hi(3),coarse_n(3),parent,ratio,i,j,k,q,idx(3),ios
    integer :: supplied,minimum_supplied,maximum_supplied,axis_min,axis_max
    integer(int64) :: cells
    real(dp) :: frame_values(4),frame_min(4),frame_max(4)
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    supplied=merge(1,0,present(internal_energy))
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
      valid=frame%is_valid();call agree(valid)
      if(.not.valid) return
      call MPI_Allreduce(frame%axis,axis_min,1,MPI_INTEGER,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(frame%axis,axis_max,1,MPI_INTEGER,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS) return
      frame_values=[frame%angular_speed,frame%axis_point]
      call MPI_Allreduce(frame_values,frame_min,4,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(frame_values,frame_max,4,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.axis_min/=axis_max.or.any(frame_min/=frame_max)) return
    end if
    call validate_mpi_amr_tree_metadata_3d(tree,owners,nv,comm,valid)
    if (.not.valid) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid)
    if (.not.valid) return
    count=size(ratios)
    valid=finest%is_valid()
    if (valid) valid=all(int([finest%nx,finest%ny,finest%nz],int64)==domain(:,count))
    call agree(valid)
    if (.not.valid) return
    allocate(levels(count),shapes(size(owners)),candidate%fields(size(owners)),ct%fields(size(owners)))
    if(present(internal_energy)) then
      allocate(ch%fields(size(owners)),stat=ios)
      valid=ios==0;call agree(valid)
      if(.not.valid) return
    endif
    levels(count)=finest
    do level=count,2,-1
      call coarsen_eb_geometry_3d(levels(level),ratios(level),levels(level-1),valid)
      if (.not.valid) exit
    end do
    call agree(valid)
    if (.not.valid) return
    do node=1,size(owners)
      level=tree%nodes(node)%level+1;d=tree%nodes(node)%shape
      lo=int(origin(:,node))+1;hi=int(origin(:,node))+d
      call extract_eb_geometry_3d(levels(level),lo,hi,shapes(node),valid)
      if (.not.valid) exit
      if (owners(node)/=rank) cycle
      allocate(candidate%fields(node)%values(nv,d(1),d(2),d(3)),ct%fields(node)%values(1,d(1),d(2),d(3)),stat=ios)
      valid=ios==0
      if (.not.valid) exit
      candidate%fields(node)%values=0._dp;ct%fields(node)%values=0._dp
      if(present(internal_energy)) then
        allocate(ch%fields(node)%values(1,d(1),d(2),d(3)),source=0._dp,stat=ios)
        valid=ios==0
        if(.not.valid) exit
      endif
    end do
    call agree(valid)
    if (.not.valid) return
    call validate_mpi_eb_tree_geometry_3d(tree,owners,shapes,nv,comm,valid)
    if (.not.valid) return
    if (rank==owners(1)) then
      d=tree%nodes(1)%shape
      valid=all(shape(root_state)==[nv,d]).and.all(shape(root_temperature)==d)
      if (valid) then
        candidate%fields(1)%values=root_state;ct%fields(1)%values(1,:,:,:)=root_temperature
        if(present(internal_energy)) then
          allocate(root_held(d(1),d(2),d(3)),source=0._dp,stat=ios)
          valid=ios==0
          if(valid) call initialize_nasa7_internal_energy_patch(species,root_state,root_temperature, &
            shapes(1)%volume_fraction,root_held,valid)
          if(valid) call prolong_reactive_eb_state_pcm_3d(species,shapes(1),1,root_state,root_temperature, &
            candidate%fields(1)%values,ct%fields(1)%values(1,:,:,:),valid,frame,root_held,ch%fields(1)%values(1,:,:,:))
        else
          call prolong_reactive_eb_state_pcm_3d(species,shapes(1),1,root_state,root_temperature, &
            candidate%fields(1)%values,ct%fields(1)%values(1,:,:,:),valid,frame=frame)
        endif
      end if
    end if
    call agree(valid)
    if (.not.valid) return
    do node=2,size(owners)
      parent=tree%nodes(node)%parent;ratio=tree%nodes(node)%patch%refinement_ratio;d=tree%nodes(node)%shape
      associate(p=>tree%nodes(node)%patch)
      lo=[p%coarse_i_lower,p%coarse_j_lower,p%coarse_k_lower]
      hi=[p%coarse_i_upper,p%coarse_j_upper,p%coarse_k_upper]
      end associate
      coarse_n=hi-lo+1;cells=1
      do i=1,3
        if (int(coarse_n(i),int64)>int(huge(1)/max(4,nv),int64)/cells) return
        cells=cells*int(coarse_n(i),int64)
      end do
      count=0
      if (owners(node)==rank) count=int(cells)
      allocate(requests(4,count),stat=ios)
      valid=ios==0;call agree(valid)
      if (.not.valid) return
      q=0
      if (owners(node)==rank) then
        do k=lo(3),hi(3)
        do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          q=q+1;requests(:,q)=[parent,i,j,k]
        end do
        end do
        end do
      end if
      call fetch_mpi_amr_tree_cells_3d(tree,owners,candidate,nv,requests,values,comm,valid)
      if (.not.valid) return
      call fetch_mpi_amr_tree_cells_3d(tree,owners,ct,1,requests,guesses,comm,valid)
      if (.not.valid) return
      if(present(internal_energy)) then
        call fetch_mpi_amr_tree_cells_3d(tree,owners,ch,1,requests,held_values,comm,valid)
        if(.not.valid) return
      endif
      deallocate(requests)
      if (owners(node)==rank) then
        coarse_u=reshape(values,[nv,coarse_n]);coarse_t=reshape(guesses,[coarse_n])
        do k=1,d(3)
        do j=1,d(2)
        do i=1,d(1)
          idx=1+([i,j,k]-1)/ratio
          candidate%fields(node)%values(:,i,j,k)=coarse_u(:,idx(1),idx(2),idx(3))
          ct%fields(node)%values(1,i,j,k)=coarse_t(idx(1),idx(2),idx(3))
        end do
        end do
        end do
        if(present(internal_energy)) then
          coarse_held=reshape(held_values,[coarse_n])
          call prolong_reactive_eb_state_pcm_3d(species,shapes(node),ratio,coarse_u,coarse_t, &
            candidate%fields(node)%values,ct%fields(node)%values(1,:,:,:),valid,frame, &
            coarse_held,ch%fields(node)%values(1,:,:,:))
        else
          call prolong_reactive_eb_state_pcm_3d(species,shapes(node),ratio,coarse_u,coarse_t, &
            candidate%fields(node)%values,ct%fields(node)%values(1,:,:,:),valid,frame=frame)
        endif
      end if
      call agree(valid)
      if (.not.valid) return
    end do
    call move_alloc(shapes,geometry)
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(ct%fields,temperature%fields)
    if(present(internal_energy)) call move_alloc(ch%fields,internal_energy%fields)
    ok=.true.
  contains
    subroutine agree(condition)
      logical, intent(inout) :: condition
      logical :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      condition=.false.
      if (ierr==MPI_SUCCESS) condition=result
    end subroutine
  end subroutine
end module mpi_eb_tree_initialize_3d_mod
