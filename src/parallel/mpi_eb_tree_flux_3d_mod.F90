module mpi_eb_tree_flux_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d,periodic_eb_faces_match
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d,valid_eb_integrated_flux_3d, &
    reconcile_reactive_eb_peer_flux_3d,restrict_eb_integrated_flux_3d,reflux_reactive_eb_state_3d
  use mpi_eb_tree_transfer_3d_mod, only: validate_mpi_eb_tree_geometry_3d,average_down_mpi_eb_tree_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  implicit none
  private
  public :: synchronize_mpi_eb_tree_flux_3d
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

  ! Match serial peer ordering, deepest-first replacement, and exposed-cell reflux.
  ! State, temperature and flux are owner-local; geometry is replicated.
  ! All three outputs commit only after collective success. The peers-only mode
  ! corrects provisional endpoints without restriction or descendant averaging.
  subroutine synchronize_mpi_eb_tree_flux_3d(species,tree,owners,geometry,state,temperature,flux,comm,ok, &
      minimum_level,maximum_level,peers_only,boundary_faces)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(eb_integrated_flux_3d), intent(inout) :: flux(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer, optional, intent(in) :: minimum_level,maximum_level
    logical, optional, intent(in) :: peers_only
    character(len=*), optional, intent(in) :: boundary_faces(6)
    character(len=16) :: faces(6)
    integer :: codes(6),low_codes(6),high_codes(6),seam,side
    real(dp) :: period(3),length(3)
    type(amr_tree_state_3d) :: candidate,temps
    type(eb_integrated_flux_3d), allocatable :: effective(:)
    type(eb_integrated_flux_3d) :: received,restricted
    type(eb_geometry_3d) :: section
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer(int64) :: offset(3),section_origin(3)
    logical, allocatable :: refined(:,:,:)
    integer :: nv,rank,status,node,peer,owner,receiver,axis,depth,parent,child,lo(3),hi(3),n(3),ratio,ios
    integer :: levels(3),low(3),high(3)
    logical :: valid,ready
    ok=.false.;nv=reactive_nvar(size(species))
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready)
    if (.not.ready) return
    faces='outflow'
    if(present(boundary_faces)) faces=boundary_faces
    valid=all(faces=='outflow'.or.faces=='inflow'.or.faces=='periodic'.or. &
      faces=='slip_wall'.or.faces=='no_slip_wall'.or.faces=='isothermal_wall')
    valid=valid.and.periodic_eb_faces_match(geometry(1),faces)
    call agree(valid,comm,ready)
    if(.not.ready) return
    codes=merge(1,0,faces=='inflow')+merge(2,0,faces=='periodic')+ &
      merge(3,0,faces=='slip_wall')+merge(4,0,faces=='no_slip_wall')+merge(5,0,faces=='isothermal_wall')
    call MPI_Allreduce(codes,low_codes,6,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(codes,high_codes,6,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(any(low_codes/=high_codes)) return
    associate(g=>geometry(1))
    length=[g%x_upper-g%x_lower,g%y_upper-g%y_lower,g%z_upper-g%z_lower]
    end associate
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nv,comm,ready)
    if (.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,ready)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    levels=[0,tree%level_count()-1,0]
    if (present(peers_only)) levels(3)=merge(1,0,peers_only)
    if (present(minimum_level)) levels(1)=minimum_level
    if (present(maximum_level)) levels(2)=maximum_level
    call MPI_Allreduce(levels,low,3,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(levels,high,3,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(low/=high).or.levels(1)<0.or.levels(2)<levels(1).or.levels(2)>=tree%level_count()) return
    valid=size(flux)==size(owners)
    if (valid) then
      do node=1,size(owners)
        if (owners(node)==rank) then
          if (.not.valid_eb_integrated_flux_3d(flux(node),geometry(node),nv)) valid=.false.
        else
          if (allocated(flux(node)%x).or.allocated(flux(node)%y).or.allocated(flux(node)%z)) valid=.false.
        end if
      end do
    end if
    call agree(valid,comm,ready)
    if (.not.ready) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid,comm,ready)
    if (.not.ready) return
    candidate=state;temps=temperature;effective=flux
    do node=1,size(owners)
    do peer=node+1,size(owners)
      if (tree%nodes(node)%level<levels(1).or.tree%nodes(node)%level>levels(2)) cycle
      if (tree%nodes(peer)%level/=tree%nodes(node)%level) cycle
      owner=node;receiver=peer
      do axis=1,3
        if (origin(axis,node)==origin(axis,peer)) cycle
        if (origin(axis,node)>origin(axis,peer)) then
          owner=peer;receiver=node
        end if
        exit
      end do
      offset=origin(:,receiver)-origin(:,owner)
      if (any(offset < -int(tree%nodes(receiver)%shape,int64)).or. &
        any(offset>int(tree%nodes(owner)%shape,int64))) cycle
      call send_flux(effective(owner),owners(owner),owners(receiver),tree%nodes(owner)%shape,nv,received,comm,ready)
      if (.not.ready) return
      valid=.true.
      if (rank==owners(receiver)) then
        if (levels(3)==1) then
          call reconcile_reactive_eb_peer_flux_3d(species,geometry(owner),geometry(receiver),int(offset), &
            received,effective(receiver),candidate%fields(receiver)%values,temps%fields(receiver)%values(1,:,:,:),valid)
        else
          call reconcile_reactive_eb_peer_flux_3d(species,geometry(owner),geometry(receiver),int(offset), &
            received,effective(receiver),ok=valid)
        end if
      end if
      call agree(valid,comm,ready)
      if (.not.ready) return
    end do
    end do
    do seam=1,3
      if(faces(2*seam-1)/='periodic') cycle
      period=0;period(seam)=length(seam)
      do owner=1,size(owners)
        depth=tree%nodes(owner)%level
        if(depth<levels(1).or.depth>levels(2)) cycle
        if(origin(seam,owner)+tree%nodes(owner)%shape(seam)/=domain(seam,depth+1)) cycle
        do receiver=1,size(owners)
          if(tree%nodes(receiver)%level/=depth.or.origin(seam,receiver)/=0) cycle
          offset=origin(:,receiver)-origin(:,owner)
          offset(seam)=offset(seam)+domain(seam,depth+1)
          if(.not.seam_overlap(offset,tree%nodes(owner)%shape,tree%nodes(receiver)%shape,seam)) cycle
          call send_flux(effective(owner),owners(owner),owners(receiver),tree%nodes(owner)%shape,nv,received,comm,ready)
          if(.not.ready) return
          valid=.true.
          if(rank==owners(receiver)) then
            if(levels(3)==1) then
              call reconcile_reactive_eb_peer_flux_3d(species,geometry(owner),geometry(receiver),int(offset), &
                received,effective(receiver),candidate%fields(receiver)%values, &
                temps%fields(receiver)%values(1,:,:,:),valid,periodic_length=period)
            else
              call reconcile_reactive_eb_peer_flux_3d(species,geometry(owner),geometry(receiver),int(offset), &
                received,effective(receiver),ok=valid,periodic_length=period)
            end if
          end if
          call agree(valid,comm,ready)
          if(.not.ready) return
        end do
      end do
    end do
    ! Provisional level endpoints need peer corrections before child substeps,
    ! without restriction or averaging from descendants at another time.
    if (levels(3)==1) then
      call move_alloc(candidate%fields,state%fields)
      call move_alloc(temps%fields,temperature%fields)
      flux=effective;ok=.true.
      return
    end if
    do depth=levels(2),levels(1),-1
    do node=1,size(owners)
      if (tree%nodes(node)%level/=depth) cycle
      valid=.true.
      if (rank==owners(node)) then
        n=tree%nodes(node)%shape
        allocate(refined(n(1),n(2),n(3)),stat=ios)
        valid=ios==0
        if (valid) then
          refined=.false.
          do child=node+1,size(owners)
            if (tree%nodes(child)%parent/=node) cycle
            call bounds(tree,child,lo,hi)
            refined(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.true.
          end do
          call reflux_reactive_eb_state_3d(species,geometry(node),refined,flux(node),effective(node), &
            candidate%fields(node)%values,temps%fields(node)%values(1,:,:,:),valid,all_faces=.true.)
          deallocate(refined)
        end if
      end if
      call agree(valid,comm,ready)
      if (.not.ready) return
      if (depth==levels(1)) cycle
      parent=tree%nodes(node)%parent;ratio=tree%nodes(node)%patch%refinement_ratio
      call bounds(tree,node,lo,hi)
      if (rank==owners(node)) call restrict_eb_integrated_flux_3d(geometry(node),ratio,effective(node),restricted,valid)
      call agree(valid,comm,ready)
      if (.not.ready) return
      call send_flux(restricted,owners(node),owners(parent),hi-lo+1,nv,received,comm,ready)
      if (.not.ready) return
      if (rank==owners(parent)) then
        effective(parent)%x(:,lo(1)-1:hi(1),lo(2):hi(2),lo(3):hi(3))=received%x
        effective(parent)%y(:,lo(1):hi(1),lo(2)-1:hi(2),lo(3):hi(3))=received%y
        effective(parent)%z(:,lo(1):hi(1),lo(2):hi(2),lo(3)-1:hi(3))=received%z
      end if
      call coarsen_eb_geometry_3d(geometry(node),ratio,section,valid)
      call agree(valid,comm,ready)
      if (.not.ready) return
      section_origin=origin(:,node)/int(ratio,int64)
      do peer=1,size(owners)
        if (peer==parent.or.tree%nodes(peer)%level/=depth-1) cycle
        offset=origin(:,peer)-section_origin
        if (any(offset < -int(tree%nodes(peer)%shape,int64)).or.any(offset>int(hi-lo+1,int64))) cycle
        call send_flux(restricted,owners(node),owners(peer),hi-lo+1,nv,received,comm,ready)
        if (.not.ready) return
        valid=.true.
        if (rank==owners(peer)) then
          call reconcile_reactive_eb_peer_flux_3d(species,section,geometry(peer),int(offset), &
            received,effective(peer),ok=valid)
        end if
        call agree(valid,comm,ready)
        if (.not.ready) return
      end do
      do seam=1,3
        if(faces(2*seam-1)/='periodic') cycle
        period=0;period(seam)=length(seam)
        do side=0,1
          if(side==0.and.section_origin(seam)/=0) cycle
          if(side==1.and.section_origin(seam)+hi(seam)-lo(seam)+1/=domain(seam,depth)) cycle
          do peer=1,size(owners)
            if(tree%nodes(peer)%level/=depth-1) cycle
            if(side==0.and.origin(seam,peer)+tree%nodes(peer)%shape(seam)/=domain(seam,depth)) cycle
            if(side==1.and.origin(seam,peer)/=0) cycle
            offset=origin(:,peer)-section_origin
            offset(seam)=offset(seam)+int(2*side-1,int64)*domain(seam,depth)
            if(.not.seam_overlap(offset,hi-lo+1,tree%nodes(peer)%shape,seam)) cycle
            call send_flux(restricted,owners(node),owners(peer),hi-lo+1,nv,received,comm,ready)
            if(.not.ready) return
            valid=.true.
            if(rank==owners(peer)) call reconcile_reactive_eb_peer_flux_3d(species,section,geometry(peer),int(offset), &
              received,effective(peer),ok=valid,periodic_length=period)
            call agree(valid,comm,ready)
            if(.not.ready) return
          end do
        end do
      end do
    end do
    end do
    call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,temps,comm,ready)
    if (.not.ready) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(temps%fields,temperature%fields)
    flux=effective;ok=.true.
  end subroutine

  pure logical function seam_overlap(displacement,owner_shape,receiver_shape,normal) result(overlap)
    integer(int64), intent(in) :: displacement(3)
    integer, intent(in) :: owner_shape(3),receiver_shape(3),normal
    integer :: tangent
    overlap=.true.
    do tangent=1,3
      if(tangent==normal) cycle
      if(displacement(tangent)>=owner_shape(tangent).or. &
        displacement(tangent)+receiver_shape(tangent)<=0) overlap=.false.
    end do
  end function seam_overlap

  ! A single temporary patch flux travels only from source owner to receiver.
  ! No other rank allocates its payload. Counts are checked before MPI calls.
  subroutine send_flux(source,sender,receiver,n,nv,received,comm,ok)
    type(eb_integrated_flux_3d), intent(in) :: source
    integer, intent(in) :: sender,receiver,n(3),nv
    type(eb_integrated_flux_3d), intent(inout) :: received
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer :: rank,status,counts(3),axis,d,shape4(4),ios
    integer(int64) :: elements
    logical :: valid,ready
    ok=.false.;valid=all(n>0).and.all(n<huge(1)).and.nv>0
    if (valid) then
      do axis=1,3
        shape4=[nv,n];shape4(axis+1)=shape4(axis+1)+1;elements=1_int64
        do d=1,4
          if (elements>int(huge(1),int64)/int(shape4(d),int64)) then
            valid=.false.;exit
          end if
          elements=elements*int(shape4(d),int64)
        end do
        if (.not.valid) exit
        counts(axis)=int(elements)
      end do
    end if
    call agree(valid,comm,ready)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    received=eb_integrated_flux_3d()
    if (rank==receiver) then
      allocate(received%x(nv,0:n(1),n(2),n(3)),received%y(nv,n(1),0:n(2),n(3)), &
        received%z(nv,n(1),n(2),0:n(3)),stat=ios)
      valid=ios==0
    end if
    call agree(valid,comm,ready)
    if (.not.ready) return
    if (sender==receiver) then
      if (rank==receiver) received=source
      ok=.true.;return
    end if
    do axis=1,3
      status=MPI_SUCCESS
      if (rank==sender) then
        select case(axis)
        case(1);call MPI_Send(source%x,counts(axis),MPI_DOUBLE_PRECISION,receiver,28510+axis,comm,status)
        case(2);call MPI_Send(source%y,counts(axis),MPI_DOUBLE_PRECISION,receiver,28510+axis,comm,status)
        case(3);call MPI_Send(source%z,counts(axis),MPI_DOUBLE_PRECISION,receiver,28510+axis,comm,status)
        end select
      else if (rank==receiver) then
        select case(axis)
        case(1)
          call MPI_Recv(received%x,counts(axis),MPI_DOUBLE_PRECISION,sender,28510+axis,comm,MPI_STATUS_IGNORE,status)
        case(2)
          call MPI_Recv(received%y,counts(axis),MPI_DOUBLE_PRECISION,sender,28510+axis,comm,MPI_STATUS_IGNORE,status)
        case(3)
          call MPI_Recv(received%z,counts(axis),MPI_DOUBLE_PRECISION,sender,28510+axis,comm,MPI_STATUS_IGNORE,status)
        end select
      end if
      call agree(status==MPI_SUCCESS,comm,ready)
      if (.not.ready) return
    end do
    ok=.true.
  end subroutine
end module mpi_eb_tree_flux_3d_mod
