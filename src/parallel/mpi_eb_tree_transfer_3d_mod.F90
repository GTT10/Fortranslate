module mpi_eb_tree_transfer_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_two_level_3d_mod, only: matching_eb_metrics_3d
  use eb_reactive_transfer_3d_mod, only: restrict_reactive_eb_state_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d,validate_mpi_amr_tree_metadata_3d
  implicit none
  private
  public :: validate_mpi_eb_tree_geometry_3d,average_down_mpi_eb_tree_3d,integrate_mpi_eb_tree_3d
  public :: average_down_mpi_eb_tree_field_3d
contains
  subroutine agree(condition,comm,ok)
    logical, intent(in) :: condition
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer :: status
    call MPI_Allreduce(condition,ok,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if (status/=MPI_SUCCESS) ok=.false.
  end subroutine

  subroutine bounds(tree,node,lo,hi)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: node
    integer, intent(out) :: lo(3),hi(3)
    associate(p=>tree%nodes(node)%patch)
    lo=[p%coarse_i_lower,p%coarse_j_lower,p%coarse_k_lower]
    hi=[p%coarse_i_upper,p%coarse_j_upper,p%coarse_k_upper]
    end associate
  end subroutine

  ! Geometry is replicated; U/T remain allocated only on whole-patch owners.
  ! Check every metric exactly across ranks, then the coarsened parent relation.
  subroutine validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nvar,comm,ok,padding)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer, optional, intent(in) :: padding
    integer :: width,width_min,width_max
    type(eb_geometry_3d) :: coarse,section
    real(dp), allocatable :: packed(:),reference(:)
    integer(int64) :: elements
    integer :: node,parent,lo(3),hi(3),status,ios,offset,n
    logical :: valid,ready
    ok=.false.;width=0
    if(present(padding)) width=padding
    call MPI_Allreduce(width,width_min,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(width,width_max,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.width_min/=width_max.or.width<0.or.width>huge(1)/4) return
    call validate_mpi_amr_tree_metadata_3d(tree,owners,nvar,comm,ready)
    if (.not.ready) return
    call agree(size(geometry)==size(tree%nodes),comm,ready)
    if (.not.ready) return
    valid=.true.
    do node=1,size(geometry)
      if (.not.geometry(node)%is_valid()) then
        valid=.false.;cycle
      end if
      if(any(tree%nodes(node)%shape>huge(1)-2*width)) then
        valid=.false.;cycle
      end if
      if(any([geometry(node)%nx,geometry(node)%ny,geometry(node)%nz]/=tree%nodes(node)%shape+2*width)) valid=.false.
    end do
    call agree(valid,comm,ready)
    if (.not.ready) return
    do node=1,size(geometry)
      associate(g=>geometry(node))
      ! Fixed-size optional-context record keeps collective counts equal even
      ! when allocation presence differs across ranks.
      elements=10_int64+product(int([g%nx,g%ny,g%nz],int64)+2_int64)
      elements=elements+size(g%volume_fraction,kind=int64)
      elements=elements+size(g%cell_centroid_x,kind=int64)
      elements=elements+size(g%cell_centroid_y,kind=int64)
      elements=elements+size(g%cell_centroid_z,kind=int64)
      elements=elements+size(g%x_face_fraction,kind=int64)
      elements=elements+size(g%y_face_fraction,kind=int64)
      elements=elements+size(g%z_face_fraction,kind=int64)
      elements=elements+size(g%x_face_centroid_y,kind=int64)
      elements=elements+size(g%x_face_centroid_z,kind=int64)
      elements=elements+size(g%y_face_centroid_x,kind=int64)
      elements=elements+size(g%y_face_centroid_z,kind=int64)
      elements=elements+size(g%z_face_centroid_x,kind=int64)
      elements=elements+size(g%z_face_centroid_y,kind=int64)
      elements=elements+size(g%boundary_area,kind=int64)
      elements=elements+size(g%boundary_centroid_x,kind=int64)
      elements=elements+size(g%boundary_centroid_y,kind=int64)
      elements=elements+size(g%boundary_centroid_z,kind=int64)
      elements=elements+size(g%boundary_normal_x,kind=int64)
      elements=elements+size(g%boundary_normal_y,kind=int64)
      elements=elements+size(g%boundary_normal_z,kind=int64)
      elements=elements+size(g%boundary_normal_integral_x,kind=int64)
      elements=elements+size(g%boundary_normal_integral_y,kind=int64)
      elements=elements+size(g%boundary_normal_integral_z,kind=int64)
      elements=elements+size(g%cell_type,kind=int64)
      call agree(elements<=int(huge(1),int64),comm,ready)
      if (.not.ready) return
      allocate(packed(int(elements)),reference(int(elements)),stat=ios)
      call agree(ios==0,comm,ready)
      if (.not.ready) return
      packed(:9)=[g%x_lower,g%x_upper,g%y_lower,g%y_upper,g%z_lower,g%z_upper,g%dx,g%dy,g%dz]
      packed(10)=merge(1._dp,0._dp,allocated(g%halo_volume_fraction))
      n=int(product(int([g%nx,g%ny,g%nz],int64)+2_int64))
      packed(11:10+n)=-1._dp
      if(allocated(g%halo_volume_fraction)) packed(11:10+n)=reshape(g%halo_volume_fraction,[n])
      offset=10+n
      n=size(g%volume_fraction);packed(offset+1:offset+n)=reshape(g%volume_fraction,[n]);offset=offset+n
      n=size(g%cell_centroid_x);packed(offset+1:offset+n)=reshape(g%cell_centroid_x,[n]);offset=offset+n
      n=size(g%cell_centroid_y);packed(offset+1:offset+n)=reshape(g%cell_centroid_y,[n]);offset=offset+n
      n=size(g%cell_centroid_z);packed(offset+1:offset+n)=reshape(g%cell_centroid_z,[n]);offset=offset+n
      n=size(g%x_face_fraction);packed(offset+1:offset+n)=reshape(g%x_face_fraction,[n]);offset=offset+n
      n=size(g%y_face_fraction);packed(offset+1:offset+n)=reshape(g%y_face_fraction,[n]);offset=offset+n
      n=size(g%z_face_fraction);packed(offset+1:offset+n)=reshape(g%z_face_fraction,[n]);offset=offset+n
      n=size(g%x_face_centroid_y);packed(offset+1:offset+n)=reshape(g%x_face_centroid_y,[n]);offset=offset+n
      n=size(g%x_face_centroid_z);packed(offset+1:offset+n)=reshape(g%x_face_centroid_z,[n]);offset=offset+n
      n=size(g%y_face_centroid_x);packed(offset+1:offset+n)=reshape(g%y_face_centroid_x,[n]);offset=offset+n
      n=size(g%y_face_centroid_z);packed(offset+1:offset+n)=reshape(g%y_face_centroid_z,[n]);offset=offset+n
      n=size(g%z_face_centroid_x);packed(offset+1:offset+n)=reshape(g%z_face_centroid_x,[n]);offset=offset+n
      n=size(g%z_face_centroid_y);packed(offset+1:offset+n)=reshape(g%z_face_centroid_y,[n]);offset=offset+n
      n=size(g%boundary_area);packed(offset+1:offset+n)=reshape(g%boundary_area,[n]);offset=offset+n
      n=size(g%boundary_centroid_x);packed(offset+1:offset+n)=reshape(g%boundary_centroid_x,[n]);offset=offset+n
      n=size(g%boundary_centroid_y);packed(offset+1:offset+n)=reshape(g%boundary_centroid_y,[n]);offset=offset+n
      n=size(g%boundary_centroid_z);packed(offset+1:offset+n)=reshape(g%boundary_centroid_z,[n]);offset=offset+n
      n=size(g%boundary_normal_x);packed(offset+1:offset+n)=reshape(g%boundary_normal_x,[n]);offset=offset+n
      n=size(g%boundary_normal_y);packed(offset+1:offset+n)=reshape(g%boundary_normal_y,[n]);offset=offset+n
      n=size(g%boundary_normal_z);packed(offset+1:offset+n)=reshape(g%boundary_normal_z,[n]);offset=offset+n
      n=size(g%boundary_normal_integral_x);packed(offset+1:offset+n)=reshape(g%boundary_normal_integral_x,[n]);offset=offset+n
      n=size(g%boundary_normal_integral_y);packed(offset+1:offset+n)=reshape(g%boundary_normal_integral_y,[n]);offset=offset+n
      n=size(g%boundary_normal_integral_z);packed(offset+1:offset+n)=reshape(g%boundary_normal_integral_z,[n]);offset=offset+n
      n=size(g%cell_type);packed(offset+1:offset+n)=reshape(g%cell_type,[n]);offset=offset+n
      reference=packed
      call MPI_Bcast(reference,size(reference),MPI_DOUBLE_PRECISION,0,comm,status)
      call agree(status==MPI_SUCCESS.and.all(packed==reference),comm,ready)
      if (.not.ready) return
      deallocate(packed,reference)
      end associate
      ! Padded caches are checked against physical centers by their consumer.
      if (node==1.or.width>0) cycle
      parent=tree%nodes(node)%parent
      call bounds(tree,node,lo,hi)
      call coarsen_eb_geometry_3d(geometry(node),tree%nodes(node)%patch%refinement_ratio,coarse,valid)
      if (valid) call extract_eb_geometry_3d(geometry(parent),lo,hi,section,valid)
      if (valid) valid=matching_eb_metrics_3d(coarse,section)
      call agree(valid,comm,ready)
      if (.not.ready) return
    end do
    ok=.true.
  end subroutine

  ! Restrict on the child owner; send only its coarse U/T to the parent owner.
  ! Reverse traversal handles successive remote owners without a global gather.
  subroutine average_down_mpi_eb_tree_3d(species,tree,owners,geometry,state,temperature,comm,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: candidate,temps
    real(dp), allocatable :: buffer(:,:,:,:)
    integer :: nv,node,parent,rank,status,ios,n(3),lo(3),hi(3),i,j,k,axis
    integer(int64) :: elements
    logical :: valid,ready
    ok=.false.;nv=reactive_nvar(size(species))
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready)
    if (.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nv,comm,ready)
    if (.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,ready)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    allocate(candidate%fields(size(owners)),temps%fields(size(owners)),stat=ios)
    valid=ios==0
    if (valid) then
      do node=1,size(owners)
        if (owners(node)/=rank) cycle
        allocate(candidate%fields(node)%values,source=state%fields(node)%values,stat=ios)
        if (ios/=0) then
          valid=.false.;exit
        end if
        allocate(temps%fields(node)%values,source=temperature%fields(node)%values,stat=ios)
        if (ios/=0) then
          valid=.false.;exit
        end if
      end do
    end if
    call agree(valid,comm,ready)
    if (.not.ready) return
    do node=size(owners),2,-1
      parent=tree%nodes(node)%parent
      call bounds(tree,node,lo,hi)
      n=hi-lo+1;elements=int(nv,int64)+1_int64
      valid=.true.
      do axis=1,3
        if (elements>int(huge(1),int64)/int(n(axis),int64)) then
          valid=.false.;exit
        end if
        elements=elements*int(n(axis),int64)
      end do
      call agree(valid,comm,ready)
      if (.not.ready) return
      if (rank==owners(node).or.rank==owners(parent)) then
        allocate(buffer(nv+1,n(1),n(2),n(3)),stat=ios)
        valid=ios==0
      end if
      call agree(valid,comm,ready)
      if (.not.ready) return
      if (rank==owners(node)) then
        buffer=0._dp
        call restrict_reactive_eb_state_3d(species,geometry(node),tree%nodes(node)%patch%refinement_ratio, &
          candidate%fields(node)%values,temps%fields(node)%values(1,:,:,:), &
          buffer(:nv,:,:,:),buffer(nv+1,:,:,:),valid)
      end if
      call agree(valid,comm,ready)
      if (.not.ready) return
      status=MPI_SUCCESS
      if (owners(node)/=owners(parent)) then
        if (rank==owners(node)) then
          call MPI_Send(buffer,int(elements),MPI_DOUBLE_PRECISION,owners(parent),28501,comm,status)
        else if (rank==owners(parent)) then
          call MPI_Recv(buffer,int(elements),MPI_DOUBLE_PRECISION,owners(node),28501,comm,MPI_STATUS_IGNORE,status)
        end if
      end if
      call agree(status==MPI_SUCCESS,comm,ready)
      if (.not.ready) return
      if (rank==owners(parent)) then
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          ! Solid parent payload is inactive and must remain unchanged.
          if (geometry(parent)%volume_fraction(lo(1)+i-1,lo(2)+j-1,lo(3)+k-1)<=0._dp) cycle
          candidate%fields(parent)%values(:,lo(1)+i-1,lo(2)+j-1,lo(3)+k-1)=buffer(:nv,i,j,k)
          temps%fields(parent)%values(1,lo(1)+i-1,lo(2)+j-1,lo(3)+k-1)=buffer(nv+1,i,j,k)
        end do
        end do
        end do
      end if
      if (allocated(buffer)) deallocate(buffer)
    end do
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(temps%fields,temperature%fields)
    ok=.true.
  end subroutine

  subroutine integrate_mpi_eb_tree_3d(tree,owners,geometry,state,values,comm,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: state
    real(dp), intent(inout) :: values(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    real(dp) :: candidate(size(values)),volume
    real(dp), allocatable :: local(:,:),by_node(:,:)
    logical, allocatable :: exposed(:,:,:)
    integer :: rank,status,node,child,n(3),lo(3),hi(3),i,j,k,ios
    logical :: valid,ready
    ok=.false.
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,size(values),comm,ready)
    if (.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,state,size(values),comm,ready)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    allocate(local(size(values),size(owners)),by_node(size(values),size(owners)),stat=ios)
    call agree(ios==0,comm,ready)
    if(.not.ready) return
    local=0._dp;valid=.true.
    do node=1,size(owners)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      allocate(exposed(n(1),n(2),n(3)),stat=ios)
      if (ios/=0) then
        valid=.false.;exit
      end if
      exposed=geometry(node)%volume_fraction>0._dp
      do child=node+1,size(owners)
        if (tree%nodes(child)%parent/=node) cycle
        call bounds(tree,child,lo,hi)
        exposed(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.false.
      end do
      volume=geometry(node)%dx*geometry(node)%dy*geometry(node)%dz
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if (.not.exposed(i,j,k)) cycle
        if (.not.all(ieee_is_finite(state%fields(node)%values(:,i,j,k)))) then
          valid=.false.;cycle
        end if
        local(:,node)=local(:,node)+volume*geometry(node)%volume_fraction(i,j,k)*state%fields(node)%values(:,i,j,k)
      end do
      end do
      end do
      deallocate(exposed)
    end do
    call agree(valid.and.all(ieee_is_finite(local)),comm,ready)
    if (.not.ready) return
    ! One contributor per column; then sum in tree order, independent of owners.
    call MPI_Allreduce(local,by_node,size(local),MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
    candidate=0
    if(status==MPI_SUCCESS) then
      do node=1,size(owners)
        candidate=candidate+by_node(:,node)
      end do
    end if
    call agree(status==MPI_SUCCESS.and.all(ieee_is_finite(candidate)),comm,ready)
    if (.not.ready) return
    values=candidate;ok=.true.
  end subroutine
  ! Signed arbitrary fields: fluid-volume restriction, zero for solid parents.
  subroutine average_down_mpi_eb_tree_field_3d(tree,owners,geometry,state,nv,comm,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nv
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(inout) :: state
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: candidate
    real(dp), allocatable :: buffer(:,:,:,:)
    integer :: node,parent,rank,status,ios,n(3),lo(3),hi(3),i,j,k,axis,r,first(3),last(3),c
    integer(int64) :: elements
    real(dp) :: volume
    logical :: valid,ready
    ok=.false.
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready)
    if (.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nv,comm,ready)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    allocate(candidate%fields(size(owners)),stat=ios)
    valid=ios==0
    if (valid) then
      do node=1,size(owners)
        if (owners(node)/=rank) cycle
        allocate(candidate%fields(node)%values,source=state%fields(node)%values,stat=ios)
        if (ios/=0) then
          valid=.false.;exit
        end if
        if(.not.all(ieee_is_finite(candidate%fields(node)%values))) then
          valid=.false.;exit
        end if
      end do
    end if
    call agree(valid,comm,ready)
    if (.not.ready) return
    do node=size(owners),2,-1
      parent=tree%nodes(node)%parent
      call bounds(tree,node,lo,hi)
      n=hi-lo+1;elements=int(nv,int64)
      valid=.true.
      do axis=1,3
        if (elements>int(huge(1),int64)/int(n(axis),int64)) then
          valid=.false.;exit
        end if
        elements=elements*int(n(axis),int64)
      end do
      call agree(valid,comm,ready)
      if (.not.ready) return
      if (rank==owners(node).or.rank==owners(parent)) then
        allocate(buffer(nv,n(1),n(2),n(3)),stat=ios)
        valid=ios==0
      end if
      call agree(valid,comm,ready)
      if (.not.ready) return
      if (rank==owners(node)) then
        buffer=0._dp
        r=tree%nodes(node)%patch%refinement_ratio
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          first=([i,j,k]-1)*r+1;last=first+r-1
          associate(v=>geometry(node)%volume_fraction(first(1):last(1),first(2):last(2),first(3):last(3)))
          volume=sum(v)
          if(volume>0) then
            do c=1,nv
              buffer(c,i,j,k)=sum(v*candidate%fields(node)%values(c, &
                first(1):last(1),first(2):last(2),first(3):last(3)))/volume
            end do
          end if
          end associate
        end do
        end do
        end do
        valid=all(ieee_is_finite(buffer))
      end if
      call agree(valid,comm,ready)
      if (.not.ready) return
      status=MPI_SUCCESS
      if (owners(node)/=owners(parent)) then
        if (rank==owners(node)) then
          call MPI_Send(buffer,int(elements),MPI_DOUBLE_PRECISION,owners(parent),28501,comm,status)
        else if (rank==owners(parent)) then
          call MPI_Recv(buffer,int(elements),MPI_DOUBLE_PRECISION,owners(node),28501,comm,MPI_STATUS_IGNORE,status)
        end if
      end if
      call agree(status==MPI_SUCCESS,comm,ready)
      if (.not.ready) return
      if (rank==owners(parent)) then
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          candidate%fields(parent)%values(:,lo(1)+i-1,lo(2)+j-1,lo(3)+k-1)=buffer(:nv,i,j,k)
        end do
        end do
        end do
      end if
      if (allocated(buffer)) deallocate(buffer)
    end do
    call move_alloc(candidate%fields,state%fields)
    ok=.true.
  end subroutine
end module mpi_eb_tree_transfer_3d_mod
