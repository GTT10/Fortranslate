module mpi_eb_tree_halo_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use simulation_config_reactive_3d_mod, only: mach_inlet_ghost_3d
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,periodic_eb_faces_match
  use eb_reactive_transfer_3d_mod, only: face_metrics
  use mpi_eb_tree_transfer_3d_mod, only: average_down_mpi_eb_tree_3d,validate_mpi_eb_tree_geometry_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: fetch_mpi_amr_tree_cells_3d,validate_mpi_amr_tree_local_3d
  implicit none
  private
  public :: fill_mpi_eb_tree_halo_3d
contains
  ! Collective PCM halos with outflow/Mach inlet and optional ancestor time interpolation.
  ! By default only face donors are exchanged; remaining slots are zero.
  ! Optional volume fractions request complete fluid donors including corners.
  ! Without explicit target metrics, ancestor geometry uses PCM volume fractions.
  ! Geometry is replicated. Spatial bounds are 1-halo_width:N+halo_width (default width 1).
  ! Wider halos require owner-local target_volume_fraction fields on those bounds,
  ! from the target resolution's extended geometry; coarse metric substitution is forbidden.
  subroutine fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,state,temperature,ghost,ghost_temperature,comm,ok, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction, &
      inflow_mach,inflow_temperature,inflow_pressure,ghost_volume_fraction,halo_width,target_volume_fraction, &
      boundary_faces,inflow_state,restrict_snapshot)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    type(amr_tree_state_3d), intent(inout) :: ghost,ghost_temperature
    type(amr_tree_state_3d), optional, intent(inout) :: ghost_volume_fraction
    logical, optional, intent(in) :: restrict_snapshot
    integer, optional, intent(in) :: halo_width
    type(amr_tree_state_3d), optional, intent(in) :: target_volume_fraction
    type(amr_tree_state_3d) :: cv
    type(MPI_Comm), intent(in) :: comm
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    ! One fraction per ancestor depth (root first), length max tree level.
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    type(amr_tree_state_3d) :: snapshot,st,candidate,ct,interpolated,it
    real(dp), allocatable :: fraction_min(:),fraction_max(:),temporal_values(:,:),temporal_t(:,:)
    integer, allocatable :: temporal_requests(:,:),temporal_index(:)
    integer :: supplied,minimum_supplied,maximum_supplied,temporal_count
    real(dp) :: fraction
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:),requests(:,:),destination(:,:)
    real(dp), allocatable :: values(:,:),guesses(:,:)
    real(dp) :: primitive(reactive_nprim(size(species))),recovered,speed,metrics(3)
    integer(int64) :: maximum,faces,global(3),mapped(3),cell(3),offset(3)
    integer :: nv,rank,status,ios,node,target,n(3),level,depth,axis,side,first(3),last(3),idx(3),face(3)
    integer :: i,j,k,donor,source(3),count,q,width,wrap_axis
    real(dp) :: controls(30+6*reactive_nvar(size(species)))
    real(dp) :: low_controls(size(controls)),high_controls(size(controls)),inlet(reactive_nvar(size(species)))
    real(dp) :: inlet_t
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    integer :: inlet_face,momentum_sign(3)
    logical :: valid,ready,full_halo,restrict_state
    ok=.false.;nv=reactive_nvar(size(species))
    restrict_state=.true.
    if(present(restrict_snapshot)) restrict_state=restrict_snapshot
    full_halo=present(ghost_volume_fraction)
    width=1
    if(present(halo_width)) width=halo_width
    ! Agree on finite inlet controls before any rank starts donor exchange.
    controls(:6)=-1;controls(7:)=0;valid=.true.
    if (present(inflow_temperature).or.present(inflow_pressure)) valid=present(inflow_mach)
    if(present(boundary_faces)) then
      valid=.not.present(inflow_mach).and..not.present(inflow_pressure)
      valid=valid.and.all(boundary_faces=='inflow'.or.boundary_faces=='outflow'.or.boundary_faces=='periodic'.or. &
        boundary_faces=='slip_wall'.or.boundary_faces=='no_slip_wall'.or.boundary_faces=='isothermal_wall')
      controls(22)=1
      controls(23:28)=merge(1._dp,0._dp,boundary_faces=='inflow')+ &
        merge(2._dp,0._dp,boundary_faces=='periodic')+merge(3._dp,0._dp,boundary_faces=='slip_wall')+ &
        merge(4._dp,0._dp,boundary_faces=='no_slip_wall')+merge(5._dp,0._dp,boundary_faces=='isothermal_wall')
      if(any(boundary_faces=='inflow')) valid=valid.and.present(inflow_state).and.present(inflow_temperature)
    else
      valid=valid.and..not.present(inflow_state)
    end if
    if(present(inflow_state)) then
      controls(29)=1
      if(all(shape(inflow_state)==[nv,6])) then
        controls(30:29+6*nv)=reshape(inflow_state,[6*nv])
      else
        valid=.false.
      end if
    end if
    if (present(inflow_mach)) controls(:6)=inflow_mach
    if (present(inflow_temperature)) controls(7:12)=inflow_temperature
    if (present(inflow_pressure)) controls(13:18)=inflow_pressure
    controls(size(controls))=merge(1._dp,0._dp,restrict_state)
    controls(19)=merge(1._dp,0._dp,full_halo)
    controls(20)=real(width,dp);controls(21)=merge(1._dp,0._dp,present(target_volume_fraction))
    valid=valid.and.width>=1.and.width<=huge(1)/4
    if(width>1) valid=valid.and.full_halo.and.present(target_volume_fraction)
    if(present(target_volume_fraction)) valid=valid.and.full_halo
    valid=valid.and.all(ieee_is_finite(controls))
    call agree(valid,ready)
    if (.not.ready) return
    valid=all(controls(:6)<=100).and.all(controls(13:18)>=0).and.all(controls(13:18)<=1e15_dp)
    valid=valid.and..not.any(controls(13:18)>0.and.controls(:6)<0)
    valid=valid.and..not.any(controls(7:12)<=0.and.controls(:6)>=0)
    call agree(valid,ready)
    if (.not.ready) return
    call MPI_Allreduce(controls,low_controls,size(controls),MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high_controls,size(controls),MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(low_controls/=high_controls)) return
    ! SDC stage snapshots retain provisional parent states. Both paths validate
    ! geometry, topology, owners and input shapes before donor exchange.
    snapshot=state;st=temperature
    if(restrict_state) then
      call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,snapshot,st,comm,ready)
      if (.not.ready) return
    else
      call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready)
      if (.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,state,nv,comm,ready)
      if (.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,ready)
      if (.not.ready) return
    end if
    valid=.true.
    if(present(boundary_faces)) valid=periodic_eb_faces_match(geometry(1),boundary_faces)
    call agree(valid,ready)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid,ready)
    if (.not.ready) return
    supplied=sum(merge(1,0,[present(ancestor_start),present(ancestor_end),present(temperature_start), &
      present(temperature_end),present(ancestor_fraction)]))
    call MPI_Allreduce(supplied,minimum_supplied,1,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(supplied,maximum_supplied,1,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (minimum_supplied/=maximum_supplied) return
    if (supplied/=0.and.supplied/=5) return
    if (supplied==5) then
      call agree(size(ancestor_fraction)==size(ratios)-1,ready)
      if (.not.ready) return
      call agree(all(ieee_is_finite(ancestor_fraction)),ready)
      if (.not.ready) return
      call agree(all(ancestor_fraction>=0._dp).and.all(ancestor_fraction<=1._dp),ready)
      if (.not.ready) return
      allocate(fraction_min(size(ancestor_fraction)),fraction_max(size(ancestor_fraction)),stat=ios)
      call agree(ios==0,ready)
      if (.not.ready) return
      call MPI_Allreduce(ancestor_fraction,fraction_min,size(ancestor_fraction),MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if (status/=MPI_SUCCESS) return
      call MPI_Allreduce(ancestor_fraction,fraction_max,size(ancestor_fraction),MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if (status/=MPI_SUCCESS) return
      if (any(fraction_min/=fraction_max)) return
      call validate_mpi_amr_tree_local_3d(tree,owners,ancestor_start,nv,comm,ready)
      if (.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,ancestor_end,nv,comm,ready)
      if (.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,temperature_start,1,comm,ready)
      if (.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,temperature_end,1,comm,ready)
      if (.not.ready) return
      ! Endpoint fields share the supplied topology, owners and immutable geometry.
      ! Provisional ancestor updates must not be averaged from their descendants.
      interpolated=ancestor_start;it=temperature_start
      do node=1,size(owners)
        if (owners(node)/=rank) cycle
        depth=tree%nodes(node)%level+1
        if (depth>size(ancestor_fraction)) cycle
        fraction=ancestor_fraction(depth)
        if (fraction==0._dp) cycle
        if (fraction==1._dp) then
          interpolated%fields(node)%values=ancestor_end%fields(node)%values
          it%fields(node)%values=temperature_end%fields(node)%values
        else
          interpolated%fields(node)%values=(1._dp-fraction)*ancestor_start%fields(node)%values+ &
            fraction*ancestor_end%fields(node)%values
          it%fields(node)%values=(1._dp-fraction)*temperature_start%fields(node)%values+ &
            fraction*temperature_end%fields(node)%values
        end if
      end do
    end if
    if(present(target_volume_fraction)) then
      valid=allocated(target_volume_fraction%fields)
      call agree(valid,ready)
      if(.not.ready) return
      call agree(size(target_volume_fraction%fields)==size(owners),ready)
      if(.not.ready) return
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        if(any(n>huge(1)-2*width)) then
          valid=.false.;exit
        end if
        if(.not.allocated(target_volume_fraction%fields(node)%values)) then
          valid=.false.;exit
        end if
        associate(v=>target_volume_fraction%fields(node)%values)
        if(any(shape(v)/=[1,n+2*width]).or.any(lbound(v)/=[1,1-width,1-width,1-width])) then
          valid=.false.;exit
        end if
        if(.not.all(ieee_is_finite(v))) then
          valid=.false.;exit
        end if
        if(any(v<0).or.any(v>1)) then
          valid=.false.;exit
        end if
        if(any(v(1,1:n(1),1:n(2),1:n(3))/=geometry(node)%volume_fraction)) valid=.false.
        end associate
      end do
      call agree(valid,ready)
      if(.not.ready) return
    end if
    maximum=0_int64;valid=.true.
    do node=1,size(owners)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      if (any(n>huge(1)-2*width)) then
        valid=.false.;exit
      end if
      if(full_halo) then
        faces=1
        do axis=1,3
          if(faces>int(huge(1)/max(nv,4),int64)/int(n(axis)+2*width,int64)) then
            valid=.false.;exit
          end if
          faces=faces*int(n(axis)+2*width,int64)
        end do
        if(.not.valid) exit
        if(faces>int(huge(1)/max(nv,4),int64)-maximum) then
          valid=.false.;exit
        end if
        maximum=maximum+faces
        cycle
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
    if(valid.and.full_halo) then
      allocate(cv%fields(size(owners)),stat=ios)
      valid=ios==0
    end if
    if (valid) then
      do node=1,size(owners)
        if (owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        allocate(candidate%fields(node)%values(nv,1-width:n(1)+width,1-width:n(2)+width,1-width:n(3)+width), &
          ct%fields(node)%values(1,1-width:n(1)+width,1-width:n(2)+width,1-width:n(3)+width),stat=ios)
        if (ios/=0) then
          valid=.false.;exit
        end if
        candidate%fields(node)%values=0._dp;ct%fields(node)%values=0._dp
        if(full_halo) then
          allocate(cv%fields(node)%values(1,1-width:n(1)+width,1-width:n(2)+width,1-width:n(3)+width),stat=ios)
          if(ios/=0) then
            valid=.false.;exit
          end if
          cv%fields(node)%values=0
        end if
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
        if(full_halo) then
          if(axis/=1.or.side/=0) cycle
          first=1-width;last=n+width
        end if
        do k=first(3),last(3)
        do j=first(2),last(2)
        do i=first(1),last(1)
          idx=[i,j,k];face=idx;face(axis)=side*n(axis)
          if(.not.full_halo) then
            metrics=face_metrics(geometry(target),axis,face)
            if (metrics(1)<=0._dp) cycle
          end if
          global=origin(:,target)+int(idx-1,int64);donor=0
          call map_wall_donor(global,momentum_sign)
          if(full_halo) global=max(0_int64,min(domain(:,level)-1,global))
          if(full_halo.and.all(idx>=1).and.all(idx<=n)) then
            donor=target;source=idx
          else if (any(global<0_int64).or.any(global>=domain(:,level))) then
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
          if(full_halo) cv%fields(target)%values(1,i,j,k)= &
            geometry(donor)%volume_fraction(source(1),source(2),source(3))
          if(full_halo) then
            if(present(target_volume_fraction)) then
              cv%fields(target)%values(1,i,j,k)=target_volume_fraction%fields(target)%values(1,i,j,k)
            else if(allocated(geometry(target)%halo_volume_fraction)) then
              cell=global-origin(:,target)+1
              do wrap_axis=1,3
                if(controls(22+2*wrap_axis-1)>=2.or.controls(22+2*wrap_axis)>=3) &
                  cell(wrap_axis)=idx(wrap_axis)
              end do
              if(geometry(target)%halo_volume_fraction(cell(1),cell(2),cell(3))>=0) &
                cv%fields(target)%values(1,i,j,k)=geometry(target)%halo_volume_fraction(cell(1),cell(2),cell(3))
            end if
            if(cv%fields(target)%values(1,i,j,k)<=0) cycle
          end if
          if (geometry(donor)%volume_fraction(source(1),source(2),source(3))<=0._dp) then
            if(.not.full_halo.or.width>1) valid=.false.
            cycle
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
    if (supplied==5) then
      allocate(temporal_requests(4,count),temporal_index(count),stat=ios)
      call agree(ios==0,ready)
      if (.not.ready) return
      temporal_count=0
      do q=1,count
        donor=requests(1,q);target=destination(1,q)
        if (tree%nodes(donor)%level>=tree%nodes(target)%level) cycle
        temporal_count=temporal_count+1
        temporal_requests(:,temporal_count)=requests(:,q);temporal_index(temporal_count)=q
      end do
      call fetch_mpi_amr_tree_cells_3d(tree,owners,interpolated,nv, &
        temporal_requests(:,:temporal_count),temporal_values,comm,ready)
      if (.not.ready) return
      call fetch_mpi_amr_tree_cells_3d(tree,owners,it,1, &
        temporal_requests(:,:temporal_count),temporal_t,comm,ready)
      if (.not.ready) return
      do q=1,temporal_count
        values(:,temporal_index(q))=temporal_values(:,q)
        guesses(:,temporal_index(q))=temporal_t(:,q)
      end do
    end if
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
      level=tree%nodes(node)%level+1
      global=origin(:,node)+int([i,j,k]-1,int64);inlet_face=0
      do axis=1,3
        if(full_halo) then
          if(inlet_face/=0) cycle
          if(global(axis)<0.and.(controls(2*axis-1)>=0.or.controls(22+2*axis-1)==1)) inlet_face=2*axis-1
          if(global(axis)>=domain(axis,level).and.(controls(2*axis)>=0.or.controls(22+2*axis)==1)) inlet_face=2*axis
        else
          if (global(axis)<0) inlet_face=2*axis-1
          if (global(axis)>=domain(axis,level)) inlet_face=2*axis
        end if
      end do
      if (inlet_face>0) then
        if(controls(22+inlet_face)==1) then
          call reactive_conserved_to_primitive(species,inflow_state(:,inlet_face),controls(6+inlet_face), &
            primitive,recovered,speed,ready)
          if(.not.ready) then
            valid=.false.;cycle
          end if
          values(:,q)=inflow_state(:,inlet_face)
        else if (controls(inlet_face)>=0) then
          if (controls(12+inlet_face)>0) then
            call mach_inlet_ghost_3d(species,values(:,q),recovered,inlet_face,controls(6+inlet_face), &
              controls(inlet_face),inlet,inlet_t,ready,inlet_pressure=controls(12+inlet_face))
          else
            call mach_inlet_ghost_3d(species,values(:,q),recovered,inlet_face,controls(6+inlet_face), &
              controls(inlet_face),inlet,inlet_t,ready)
          end if
          if (.not.ready) then
            valid=.false.;cycle
          end if
          values(:,q)=inlet;recovered=inlet_t
        end if
      end if
      call map_wall_donor(global,momentum_sign)
      values(2:4,q)=values(2:4,q)*real(momentum_sign,dp)
      candidate%fields(node)%values(:,i,j,k)=values(:,q)
      ct%fields(node)%values(1,i,j,k)=recovered
    end do
    call agree(valid,ready)
    if (.not.ready) return
    call move_alloc(candidate%fields,ghost%fields)
    call move_alloc(ct%fields,ghost_temperature%fields)
    if(full_halo) call move_alloc(cv%fields,ghost_volume_fraction%fields)
    ok=.true.
  contains
    subroutine map_wall_donor(position,signs)
      integer(int64), intent(inout) :: position(3)
      integer, intent(out) :: signs(3)
      integer :: direction,boundary_code,boundary_face
      signs=1
      do direction=1,3
        if(controls(22+2*direction-1)==2) then
          position(direction)=modulo(position(direction),domain(direction,level))
          cycle
        end if
        do while(position(direction)<0.or.position(direction)>=domain(direction,level))
          boundary_face=2*direction
          if(position(direction)<0) boundary_face=boundary_face-1
          boundary_code=nint(controls(22+boundary_face))
          if(boundary_code<3.or.boundary_code>5) exit
          if(boundary_code>=4) then
            signs=-signs
          else
            signs(direction)=-signs(direction)
          end if
          if(position(direction)<0) then
            position(direction)=-position(direction)-1
          else
            position(direction)=2*domain(direction,level)-position(direction)-1
          end if
        end do
      end do
    end subroutine
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
