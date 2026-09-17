module mpi_eb_tree_halo_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use eb_reactive_transfer_3d_mod, only: face_metrics
  use mpi_eb_tree_transfer_3d_mod, only: average_down_mpi_eb_tree_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: fetch_mpi_amr_tree_cells_3d
  implicit none
  private
  public :: fill_mpi_eb_tree_halo_3d
contains
  ! Collective spatial PCM/outflow halos for all locally owned patches.
  ! Only face donors are exchanged; interior, solid faces, edges and corners are zero.
  ! Geometry is replicated. Output fields have spatial lower bounds zero.
  subroutine fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,state,temperature,ghost,ghost_temperature,comm,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    type(amr_tree_state_3d), intent(inout) :: ghost,ghost_temperature
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: snapshot,st,candidate,ct
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:),requests(:,:),destination(:,:)
    real(dp), allocatable :: values(:,:),guesses(:,:)
    real(dp) :: primitive(reactive_nprim(size(species))),recovered,speed,metrics(3)
    integer(int64) :: maximum,faces,global(3),mapped(3),cell(3),offset(3)
    integer :: nv,rank,status,ios,node,target,n(3),level,depth,axis,side,first(3),last(3),idx(3),face(3)
    integer :: i,j,k,donor,source(3),count,q
    logical :: valid,ready
    ok=.false.;nv=reactive_nvar(size(species))
    ! Average only a private owner-local snapshot. This also validates geometry,
    ! topology, owners and both input container shapes collectively.
    snapshot=state;st=temperature
    call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,snapshot,st,comm,ready)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid,ready)
    if (.not.ready) return
    maximum=0_int64;valid=.true.
    do node=1,size(owners)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      if (any(n>huge(1)-2)) then
        valid=.false.;exit
      end if
      do axis=1,3
        faces=int(n(axis),int64)*int(n(mod(axis,3)+1),int64)
        if (faces>(int(huge(1)/max(nv,4),int64)-maximum)/2_int64) then
          valid=.false.;exit
        end if
        maximum=maximum+2_int64*faces
      end do
      if (.not.valid) exit
    end do
    call agree(valid,ready)
    if (.not.ready) return
    allocate(requests(4,int(maximum)),destination(4,int(maximum)), &
      candidate%fields(size(owners)),ct%fields(size(owners)),stat=ios)
    valid=ios==0
    if (valid) then
      do node=1,size(owners)
        if (owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        allocate(candidate%fields(node)%values(nv,0:n(1)+1,0:n(2)+1,0:n(3)+1), &
          ct%fields(node)%values(1,0:n(1)+1,0:n(2)+1,0:n(3)+1),stat=ios)
        if (ios/=0) then
          valid=.false.;exit
        end if
        candidate%fields(node)%values=0._dp;ct%fields(node)%values=0._dp
      end do
    end if
    call agree(valid,ready)
    if (.not.ready) return
    count=0
    do target=1,size(owners)
      if (owners(target)/=rank) cycle
      n=tree%nodes(target)%shape;level=tree%nodes(target)%level+1
      ! Validate closed peer faces too: skipping a closed target request must
      ! not conceal an open face on its peer.
      do node=1,size(owners)
        if (node==target.or.tree%nodes(node)%level+1/=level) cycle
        offset=origin(:,node)-origin(:,target)
        if (any(offset < -int(tree%nodes(node)%shape,int64)).or.any(offset>int(n,int64))) cycle
        if (.not.peer_faces_match(geometry(target),geometry(node),int(offset))) valid=.false.
      end do
      do axis=1,3
      do side=0,1
        first=1;last=n;first(axis)=side*(n(axis)+1);last(axis)=first(axis)
        do k=first(3),last(3)
        do j=first(2),last(2)
        do i=first(1),last(1)
          idx=[i,j,k];face=idx;face(axis)=side*n(axis)
          metrics=face_metrics(geometry(target),axis,face)
          if (metrics(1)<=0._dp) cycle
          global=origin(:,target)+int(idx-1,int64);donor=0
          if (any(global<0_int64).or.any(global>=domain(:,level))) then
            donor=target;source=max(1,min(n,idx))
          else
            do depth=level,1,-1
              mapped=global/(domain(:,level)/domain(:,depth))
              do node=1,size(owners)
                if (tree%nodes(node)%level+1/=depth) cycle
                cell=mapped-origin(:,node)
                if (any(cell<0_int64).or.any(cell>=int(tree%nodes(node)%shape,int64))) cycle
                donor=node;source=int(cell)+1;exit
              end do
              if (donor/=0) exit
            end do
          end if
          if (donor==0) then
            valid=.false.;cycle
          end if
          if (geometry(donor)%volume_fraction(source(1),source(2),source(3))<=0._dp) then
            valid=.false.;cycle
          end if
          count=count+1
          requests(:,count)=[donor,source];destination(:,count)=[target,idx]
        end do
        end do
        end do
      end do
      end do
    end do
    call agree(valid,ready)
    if (.not.ready) return
    call fetch_mpi_amr_tree_cells_3d(tree,owners,snapshot,nv,requests(:,:count),values,comm,ready)
    if (.not.ready) return
    call fetch_mpi_amr_tree_cells_3d(tree,owners,st,1,requests(:,:count),guesses,comm,ready)
    if (.not.ready) return
    do q=1,count
      if (.not.ieee_is_finite(guesses(1,q))) then
        valid=.false.;cycle
      end if
      if (guesses(1,q)<=0._dp) then
        valid=.false.;cycle
      end if
      call reactive_conserved_to_primitive(species,values(:,q),guesses(1,q),primitive,recovered,speed,ready)
      if (.not.ready) then
        valid=.false.;cycle
      end if
      node=destination(1,q);i=destination(2,q);j=destination(3,q);k=destination(4,q)
      candidate%fields(node)%values(:,i,j,k)=values(:,q)
      ct%fields(node)%values(1,i,j,k)=recovered
    end do
    call agree(valid,ready)
    if (.not.ready) return
    call move_alloc(candidate%fields,ghost%fields)
    call move_alloc(ct%fields,ghost_temperature%fields)
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

  logical function peer_faces_match(target,peer,offset) result(valid)
    type(eb_geometry_3d), intent(in) :: target,peer
    integer, intent(in) :: offset(3)
    integer :: n(3),pn(3),axis,side,first(3),last(3),idx(3),source(3),face(3),other(3),i,j,k
    real(dp) :: spacing(3),origin(3),expected(3),tolerance,left(3),right(3)
    valid=.false.
    n=[target%nx,target%ny,target%nz];pn=[peer%nx,peer%ny,peer%nz]
    if (all(offset<n.and.offset+pn>0)) return
    spacing=[target%dx,target%dy,target%dz]
    if (any(abs([peer%dx,peer%dy,peer%dz]-spacing)>512._dp*epsilon(1._dp)*spacing)) return
    origin=[target%x_lower,target%y_lower,target%z_lower];expected=origin+real(offset,dp)*spacing
    tolerance=512._dp*epsilon(1._dp)*max(1._dp,maxval(abs(origin)),maxval(abs(expected)))
    if (any(abs([peer%x_lower,peer%y_lower,peer%z_lower]-expected)>tolerance)) return
    do axis=1,3
    do side=0,1
      if (side==0.and.offset(axis)+pn(axis)/=0) cycle
      if (side==1.and.offset(axis)/=n(axis)) cycle
      first=max(1,offset+1);last=min(n,offset+pn)
      first(axis)=side*(n(axis)+1);last(axis)=first(axis)
      if (any(first>last)) cycle
      do k=first(3),last(3)
      do j=first(2),last(2)
      do i=first(1),last(1)
        idx=[i,j,k];source=idx-offset
        face=idx;face(axis)=side*n(axis)
        other=source;other(axis)=source(axis)-side
        left=face_metrics(target,axis,face);right=face_metrics(peer,axis,other)
        if (abs(left(1)-right(1))>512._dp*epsilon(1._dp)) return
        if (left(1)<=0._dp) cycle
        if (any(abs(left(2:3)-right(2:3))>512._dp*epsilon(1._dp))) return
      end do
      end do
      end do
    end do
    end do
    valid=.true.
  end function
end module mpi_eb_tree_halo_3d_mod
