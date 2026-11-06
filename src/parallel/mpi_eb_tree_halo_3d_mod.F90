module mpi_eb_tree_halo_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use dual_energy_mod, only: initialize_nasa7_internal_energy_patch
  use nasa7_thermo_mod, only: nasa7_species
  use rotating_frame_3d_mod, only: rotating_frame_3d,convert_rotating_energy_3d
  use simulation_config_reactive_3d_mod, only: mach_inlet_ghost_3d,valid_reactive_boundaries_3d
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,periodic_eb_faces_match
  use eb_reactive_transfer_3d_mod, only: face_metrics,recover_reactive_eb_cell_3d, &
    plan_eb_linear_halo_3d,interpolate_eb_linear_halo_3d,reconstruct_eb_pressure_state_3d, &
    interpolate_eb_pressure_time_3d
  use mpi_eb_tree_transfer_3d_mod, only: average_down_mpi_eb_tree_3d,validate_mpi_eb_tree_geometry_3d, &
    average_down_mpi_eb_tree_field_3d
  use amr_patch_tree_remap_3d_mod, only: covered
  use mpi_amr_patch_tree_cells_3d_mod, only: fetch_mpi_amr_tree_cells_3d,validate_mpi_amr_tree_local_3d
  implicit none
  private
  public :: fill_mpi_eb_tree_halo_3d,fill_sdc_mpi_eb_tree_source_halo_3d
  public :: fill_sdc_mpi_eb_particle_predictor_halo_3d
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
      boundary_faces,inflow_state,restrict_snapshot,frame,internal_energy,ghost_internal_energy, &
      ancestor_internal_energy_start,ancestor_internal_energy_end,coarse_geometry,separate_momentum, &
      consistent_velocity,paired_pressure,pressure_level,connected_directions,preserve_temperature_seed)
    logical, optional, intent(in) :: preserve_temperature_seed
    integer, optional, intent(in) :: pressure_level
    type(eb_geometry_3d), optional, intent(in) :: coarse_geometry(:)
    logical, optional, intent(in) :: separate_momentum,consistent_velocity,paired_pressure,connected_directions
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy,ancestor_internal_energy_start, &
      ancestor_internal_energy_end
    type(amr_tree_state_3d), intent(inout) :: ghost,ghost_temperature
    type(amr_tree_state_3d), optional, intent(inout) :: ghost_internal_energy
    type(amr_tree_state_3d), optional, intent(inout) :: ghost_volume_fraction
    logical, optional, intent(in) :: restrict_snapshot
    integer, optional, intent(in) :: halo_width
    type(amr_tree_state_3d), optional, intent(in) :: target_volume_fraction
    type(amr_tree_state_3d) :: cv,pressure_field,temporal_pressure
    type(MPI_Comm), intent(in) :: comm
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    ! One fraction per ancestor depth (root first), length max tree level.
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    type(amr_tree_state_3d) :: snapshot,st,candidate,ct,interpolated,it,held_snapshot,held_candidate,held_interpolated
    real(dp), allocatable :: fraction_min(:),fraction_max(:),temporal_values(:,:),temporal_t(:,:), &
      held_values(:,:),temporal_held(:,:)
    integer, allocatable :: temporal_requests(:,:),temporal_index(:)
    integer :: supplied,minimum_supplied,maximum_supplied,temporal_count
    integer :: held_controls(8),held_min(8),held_max(8)
    real(dp) :: fraction
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: node_shapes(:,:)
    integer :: shape_node
    integer, allocatable :: ratios(:),requests(:,:),destination(:,:)
    real(dp), allocatable :: values(:,:),guesses(:,:),donor_held
    real(dp) :: primitive(reactive_nprim(size(species))),recovered,speed,metrics(3)
    integer(int64) :: maximum,faces,global(3),mapped(3),cell(3),offset(3)
    integer :: nv,rank,status,ios,node,target,n(3),level,depth,axis,side,first(3),last(3),idx(3),face(3)
    integer :: i,j,k,donor,source(3),count,q,width,wrap_axis
    real(dp) :: controls(30+6*reactive_nvar(size(species)))
    real(dp) :: low_controls(size(controls)),high_controls(size(controls)),inlet(reactive_nvar(size(species)))
    real(dp) :: inlet_t
    real(dp) :: position(3),energy,frame_values(4),frame_min(4),frame_max(4)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    integer :: inlet_face,momentum_sign(3)
    integer :: frame_present,frame_min_present,frame_max_present,axis_min,axis_max
    logical :: valid,ready,full_halo,restrict_state,linear
    character(len=16) :: linear_faces(6)
    integer :: stencil(4,6),signs(3,6),stencil_count,r,h
    integer, allocatable :: stencil_requests(:,:),stencil_index(:,:),stencil_signs(:,:,:)
    real(dp), allocatable :: stencil_values(:,:),stencil_t(:,:),stencil_held(:,:),stencil_offsets(:,:)
    real(dp) :: packed(reactive_nvar(size(species))+2),neighbors(reactive_nvar(size(species))+2,6),shift(3)
    real(dp), allocatable :: pressures(:,:),stencil_pressures(:,:),time_pressures(:,:)
    real(dp) :: target_pressure(1),pressure_neighbors(1,6),energy_delta
    logical :: pressure_policy
    ok=.false.;nv=reactive_nvar(size(species))
    pressure_policy=.false.
    if(present(paired_pressure)) pressure_policy=paired_pressure
    frame_present=merge(1,0,present(frame))
    call MPI_Allreduce(frame_present,frame_min_present,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(frame_present,frame_max_present,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.frame_min_present/=frame_max_present) return
    if(present(frame)) then
      call agree(frame%is_valid(),ready)
      if(.not.ready) return
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
    held_controls=[merge(1,0,present(internal_energy)),merge(1,0,present(ghost_internal_energy)), &
      merge(1,0,present(ancestor_internal_energy_start)),merge(1,0,present(ancestor_internal_energy_end)), &
      merge(1,0,present(coarse_geometry)),0,-1,0]
    if(present(preserve_temperature_seed)) held_controls(8)=merge(1,0,preserve_temperature_seed)
    if(present(pressure_level)) held_controls(7)=pressure_level
    if(present(separate_momentum)) held_controls(6)=merge(1,0,separate_momentum)
    if(present(consistent_velocity)) held_controls(6)=held_controls(6)+merge(2,0,consistent_velocity)
    if(present(paired_pressure)) held_controls(6)=held_controls(6)+merge(4,0,paired_pressure)
    if(present(connected_directions)) held_controls(6)=held_controls(6)+merge(8,0,connected_directions)
    valid=valid.and.iand(held_controls(6),3)/=2
    if(btest(held_controls(6),3)) valid=valid.and.btest(held_controls(6),0)
    if(btest(held_controls(6),2)) valid=valid.and.iand(held_controls(6),7)==7.and. &
      .not.present(frame)
    if(held_controls(6)>0) valid=valid.and.present(coarse_geometry)
    if(present(coarse_geometry)) then
      valid=valid.and.present(internal_energy).and.width>1.and.size(coarse_geometry)==size(owners)
    endif
    valid=valid.and.held_controls(1)==held_controls(2)
    valid=valid.and.held_controls(3)==held_controls(4)
    valid=valid.and.(held_controls(1)==1.or.held_controls(3)==0)
    if(held_controls(1)==1) valid=valid.and.(present(ancestor_fraction).eqv.(held_controls(3)==1))
    valid=valid.and.all(ieee_is_finite(controls))
    call agree(valid,ready)
    if (.not.ready) return
    call MPI_Allreduce(held_controls,held_min,size(held_controls),MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(held_controls,held_max,size(held_controls),MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(held_min/=held_max)) return
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
    if(present(internal_energy)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,internal_energy,1,comm,ready)
      if(.not.ready) return
      held_snapshot=internal_energy
      allocate(donor_held,stat=ios)
      call agree(ios==0,ready)
      if(.not.ready) return
    endif
    if(restrict_state) then
      if(present(internal_energy)) then
        call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,snapshot,st,comm,ready,frame,held_snapshot)
      else
        call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,snapshot,st,comm,ready,frame)
      endif
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
    allocate(node_shapes(3,size(tree%nodes)),stat=ios)
    call agree(ios==0,ready)
    if(.not.ready) return
    do shape_node=1,size(tree%nodes)
      node_shapes(:,shape_node)=tree%nodes(shape_node)%shape
    enddo
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid,ready)
    if (.not.ready) return
    if(pressure_policy) then
      pressure_field=internal_energy
      valid=.true.
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        pressure_field%fields(node)%values=0
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(geometry(node)%volume_fraction(i,j,k)<=0) cycle
          call reactive_conserved_to_primitive(species,state%fields(node)%values(:,i,j,k), &
            temperature%fields(node)%values(1,i,j,k),primitive,recovered,speed,ready, &
            internal_energy_density=internal_energy%fields(node)%values(1,i,j,k))
          if(.not.ready) then
            valid=.false.;cycle
          endif
          pressure_field%fields(node)%values(1,i,j,k)=primitive(5)
        enddo
        enddo
        enddo
      enddo
      call agree(valid,ready)
      if(.not.ready) return
      if(held_controls(7)<0.and..not.present(ancestor_fraction)) then
        call average_down_mpi_eb_tree_field_3d(tree,owners,geometry,pressure_field,1,comm,ready)
        if(.not.ready) return
      endif
    endif
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
      if(present(internal_energy)) then
        call validate_mpi_amr_tree_local_3d(tree,owners,ancestor_internal_energy_start,1,comm,ready)
        if(.not.ready) return
        call validate_mpi_amr_tree_local_3d(tree,owners,ancestor_internal_energy_end,1,comm,ready)
        if(.not.ready) return
      endif
      ! Endpoint fields share the supplied topology, owners and immutable geometry.
      ! Provisional ancestor updates must not be averaged from their descendants.
      interpolated=ancestor_start;it=temperature_start
      if(present(internal_energy)) held_interpolated=ancestor_internal_energy_start
      do node=1,size(owners)
        if (owners(node)/=rank) cycle
        depth=tree%nodes(node)%level+1
        if (depth>size(ancestor_fraction)) cycle
        fraction=ancestor_fraction(depth)
        if (fraction==0._dp) cycle
        if (fraction==1._dp) then
          interpolated%fields(node)%values=ancestor_end%fields(node)%values
          it%fields(node)%values=temperature_end%fields(node)%values
          if(present(internal_energy)) held_interpolated%fields(node)%values= &
            ancestor_internal_energy_end%fields(node)%values
        else
          interpolated%fields(node)%values=(1._dp-fraction)*ancestor_start%fields(node)%values+ &
            fraction*ancestor_end%fields(node)%values
          it%fields(node)%values=(1._dp-fraction)*temperature_start%fields(node)%values+ &
            fraction*temperature_end%fields(node)%values
          if(present(internal_energy)) held_interpolated%fields(node)%values= &
            (1._dp-fraction)*ancestor_internal_energy_start%fields(node)%values+ &
            fraction*ancestor_internal_energy_end%fields(node)%values
        end if
      end do
    end if
    if(pressure_policy.and.supplied==5) then
      temporal_pressure=pressure_field;valid=.true.
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        depth=tree%nodes(node)%level+1
        if(depth>size(ancestor_fraction)) cycle
        n=tree%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(geometry(node)%volume_fraction(i,j,k)<=0) cycle
          call interpolate_eb_pressure_time_3d(species,ancestor_start%fields(node)%values(:,i,j,k), &
            ancestor_end%fields(node)%values(:,i,j,k),ancestor_internal_energy_start%fields(node)%values(1,i,j,k), &
            ancestor_internal_energy_end%fields(node)%values(1,i,j,k),temperature_start%fields(node)%values(1,i,j,k), &
            temperature_end%fields(node)%values(1,i,j,k),ancestor_fraction(depth), &
            temporal_pressure%fields(node)%values(1,i,j,k),ready)
          valid=valid.and.ready
        enddo
        enddo
        enddo
      enddo
      call agree(valid,ready)
      if(.not.ready) return
    endif
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
      if(present(internal_energy)) then
        allocate(held_candidate%fields(size(owners)),stat=ios)
        valid=ios==0
      endif
    endif
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
        if(present(internal_energy)) then
          allocate(held_candidate%fields(node)%values(1,1-width:n(1)+width,1-width:n(2)+width, &
            1-width:n(3)+width),stat=ios)
          if(ios/=0) then
            valid=.false.;exit
          endif
          held_candidate%fields(node)%values=0
        endif
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
    if(present(coarse_geometry)) then
      call agree(count<=huge(1)/(6*(nv+2)),ready)
      if(.not.ready) return
      allocate(stencil_requests(4,6*count),stencil_index(6,count),stencil_signs(3,6,count), &
        stencil_offsets(3,count),stat=ios)
      call agree(ios==0,ready)
      if(.not.ready) return
      stencil_count=0;stencil_index=0;linear_faces='outflow'
      if(present(boundary_faces)) linear_faces=boundary_faces
      where(controls(:6)>=0) linear_faces='inflow'
      do q=1,count
        donor=requests(1,q);target=destination(1,q)
        level=tree%nodes(target)%level+1
        if(tree%nodes(donor)%level+1>=level) cycle
        global=origin(:,target)+int(destination(2:4,q)-1,int64)
        call map_wall_donor(global,momentum_sign)
        global=max(0_int64,min(domain(:,level)-1,global))
        call plan_eb_linear_halo_3d(tree%nodes%level,node_shapes,origin,domain,donor,requests(2:4,q),level,global, &
          coarse_geometry(donor),width,linear_faces,stencil,signs,shift,linear,ready,connected_directions=connected_directions)
        if(.not.ready) then
          valid=.false.;exit
        endif
        if(.not.linear) cycle
        stencil_offsets(:,q)=shift;stencil_signs(:,:,q)=signs
        do r=1,6
          stencil_count=stencil_count+1
          stencil_requests(:,stencil_count)=stencil(:,r);stencil_index(r,q)=stencil_count
        enddo
      enddo
      call agree(valid,ready)
      if(.not.ready) return
      if(pressure_policy) then
        if(supplied==5) then
          call fetch_mpi_amr_tree_cells_3d(tree,owners,temporal_pressure,1, &
            stencil_requests(:,:stencil_count),stencil_pressures,comm,ready)
        else
          call fetch_mpi_amr_tree_cells_3d(tree,owners,pressure_field,1, &
            stencil_requests(:,:stencil_count),stencil_pressures,comm,ready)
        endif
        if(.not.ready) return
      endif
      if(supplied==5) then
        call fetch_mpi_amr_tree_cells_3d(tree,owners,interpolated,nv, &
          stencil_requests(:,:stencil_count),stencil_values,comm,ready)
        if(.not.ready) return
        call fetch_mpi_amr_tree_cells_3d(tree,owners,it,1, &
          stencil_requests(:,:stencil_count),stencil_t,comm,ready)
        if(.not.ready) return
        call fetch_mpi_amr_tree_cells_3d(tree,owners,held_interpolated,1, &
          stencil_requests(:,:stencil_count),stencil_held,comm,ready)
      else
        call fetch_mpi_amr_tree_cells_3d(tree,owners,snapshot,nv, &
          stencil_requests(:,:stencil_count),stencil_values,comm,ready)
        if(.not.ready) return
        call fetch_mpi_amr_tree_cells_3d(tree,owners,st,1, &
          stencil_requests(:,:stencil_count),stencil_t,comm,ready)
        if(.not.ready) return
        call fetch_mpi_amr_tree_cells_3d(tree,owners,held_snapshot,1, &
          stencil_requests(:,:stencil_count),stencil_held,comm,ready)
      endif
      if(.not.ready) return
    endif
    if(pressure_policy) then
      call fetch_mpi_amr_tree_cells_3d(tree,owners,pressure_field,1,requests(:,:count),pressures,comm,ready)
      if(.not.ready) return
    endif
    call fetch_mpi_amr_tree_cells_3d(tree,owners,snapshot,nv,requests(:,:count),values,comm,ready)
    if (.not.ready) return
    call fetch_mpi_amr_tree_cells_3d(tree,owners,st,1,requests(:,:count),guesses,comm,ready)
    if (.not.ready) return
    if(present(internal_energy)) then
      call fetch_mpi_amr_tree_cells_3d(tree,owners,held_snapshot,1,requests(:,:count),held_values,comm,ready)
      if(.not.ready) return
    endif
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
      if(present(internal_energy)) then
        call fetch_mpi_amr_tree_cells_3d(tree,owners,held_interpolated,1, &
          temporal_requests(:,:temporal_count),temporal_held,comm,ready)
        if(.not.ready) return
      endif
      if(pressure_policy) then
        call fetch_mpi_amr_tree_cells_3d(tree,owners,temporal_pressure,1, &
          temporal_requests(:,:temporal_count),time_pressures,comm,ready)
        if(.not.ready) return
        do q=1,temporal_count
          pressures(:,temporal_index(q))=time_pressures(:,q)
        enddo
      endif
      do q=1,temporal_count
        values(:,temporal_index(q))=temporal_values(:,q)
        guesses(:,temporal_index(q))=temporal_t(:,q)
        if(present(internal_energy)) held_values(:,temporal_index(q))=temporal_held(:,q)
      end do
    end if
    do q=1,count
      if(pressure_policy) target_pressure=pressures(:,q)
      if(present(coarse_geometry)) then
        if(stencil_index(1,q)>0) then
          packed=[values(:,q),guesses(1,q),held_values(1,q)]
          do r=1,6
            h=stencil_index(r,q)
            neighbors(:,r)=[stencil_values(:,h),stencil_t(1,h),stencil_held(1,h)]
            neighbors(2:4,r)=neighbors(2:4,r)*real(stencil_signs(:,r,q),dp)
            if(pressure_policy) pressure_neighbors(:,r)=stencil_pressures(:,h)
          enddo
          call interpolate_eb_linear_halo_3d(packed,neighbors,stencil_offsets(:,q),ready,separate_momentum,consistent_velocity)
          if(.not.ready) then
            valid=.false.;cycle
          endif
          values(:,q)=packed(:nv);guesses(1,q)=packed(nv+1);held_values(1,q)=packed(nv+2)
          if(pressure_policy) then
            call interpolate_eb_linear_halo_3d(target_pressure,pressure_neighbors,stencil_offsets(:,q),ready)
            if(.not.ready) then
              valid=.false.;cycle
            endif
          endif
        endif
      endif
      if(pressure_policy) then
        donor=requests(1,q);target=destination(1,q)
        if(tree%nodes(donor)%level<tree%nodes(target)%level.or.covered(tree,donor,requests(2:4,q))) then
          call reconstruct_eb_pressure_state_3d(species,values(:,q),held_values(1,q),guesses(1,q), &
            target_pressure(1),energy_delta,ready)
          if(.not.ready) then
            valid=.false.;cycle
          endif
        endif
      endif
      if (.not.ieee_is_finite(guesses(1,q))) then
        valid=.false.;cycle
      end if
      if (guesses(1,q)<=0._dp) then
        valid=.false.;cycle
      end if
      node=destination(1,q);i=destination(2,q);j=destination(3,q);k=destination(4,q)
      if(present(internal_energy)) donor_held=held_values(1,q)
      if(present(internal_energy)) then
        call recover_reactive_eb_cell_3d(species,geometry(node),[i,j,k],values(:,q),guesses(1,q), &
          primitive,recovered,speed,ready,frame,held_values(1,q))
      else
        call recover_reactive_eb_cell_3d(species,geometry(node),[i,j,k],values(:,q),guesses(1,q), &
          primitive,recovered,speed,ready,frame)
      endif
      if (.not.ready) then
        valid=.false.;cycle
      end if
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
          if(present(internal_energy)) then
            values(:,q)=inflow_state(:,inlet_face);recovered=controls(6+inlet_face)
            call gas_internal_energy(donor_held,ready)
            if(.not.ready) then
              valid=.false.;cycle
            endif
          endif
          call reactive_conserved_to_primitive(species,inflow_state(:,inlet_face),controls(6+inlet_face), &
            primitive,recovered,speed,ready,internal_energy_density=donor_held)
          if(.not.ready) then
            valid=.false.;cycle
          end if
          values(:,q)=inflow_state(:,inlet_face)
          call boundary_energy(.true.,ready)
          if(.not.ready) then
            valid=.false.;cycle
          end if
        else if (controls(inlet_face)>=0) then
          call boundary_energy(.false.,ready)
          if(.not.ready) then
            valid=.false.;cycle
          end if
          if (controls(12+inlet_face)>0) then
            call mach_inlet_ghost_3d(species,values(:,q),recovered,inlet_face,controls(6+inlet_face), &
              controls(inlet_face),inlet,inlet_t,ready,inlet_pressure=controls(12+inlet_face), &
              internal_energy_density=donor_held)
          else
            call mach_inlet_ghost_3d(species,values(:,q),recovered,inlet_face,controls(6+inlet_face), &
              controls(inlet_face),inlet,inlet_t,ready,internal_energy_density=donor_held)
          end if
          if (.not.ready) then
            valid=.false.;cycle
          end if
          values(:,q)=inlet;recovered=inlet_t
          call boundary_energy(.true.,ready)
          if(.not.ready) then
            valid=.false.;cycle
          end if
        end if
      end if
      call map_wall_donor(global,momentum_sign)
      values(2:4,q)=values(2:4,q)*real(momentum_sign,dp)
      candidate%fields(node)%values(:,i,j,k)=values(:,q)
      ! Keep the validated interpolation seed for the native MOL clock.
      if(held_controls(8)==1.and.inlet_face==0.and..not.pressure_policy) then
        if(tree%nodes(requests(1,q))%level<tree%nodes(node)%level) recovered=guesses(1,q)
      endif
      ct%fields(node)%values(1,i,j,k)=recovered
      if(present(internal_energy)) then
        if(inlet_face>0) then
          if(controls(22+inlet_face)==1.or.controls(inlet_face)>=0) then
            call gas_internal_energy(held_values(1,q),ready)
            if(.not.ready) then
              valid=.false.;cycle
            endif
          endif
        endif
        held_candidate%fields(node)%values(1,i,j,k)=held_values(1,q)
      endif
    end do
    call agree(valid,ready)
    if (.not.ready) return
    call move_alloc(candidate%fields,ghost%fields)
    call move_alloc(ct%fields,ghost_temperature%fields)
    if(full_halo) call move_alloc(cv%fields,ghost_volume_fraction%fields)
    if(present(ghost_internal_energy)) call move_alloc(held_candidate%fields,ghost_internal_energy%fields)
    ok=.true.
  contains
    subroutine gas_internal_energy(rhoe,accepted)
      real(dp), intent(out) :: rhoe
      logical, intent(out) :: accepted
      real(dp) :: cell(nv,1,1,1),temp(1,1,1),volume(1,1,1),scalar(1,1,1)
      cell(:,1,1,1)=values(:,q);temp=recovered;volume=1;scalar=0
      call initialize_nasa7_internal_energy_patch(species,cell,temp,volume,scalar,accepted)
      rhoe=scalar(1,1,1)
    end subroutine gas_internal_energy
    subroutine boundary_energy(to_rotating,accepted)
      logical, intent(in) :: to_rotating
      logical, intent(out) :: accepted
      accepted=.true.
      if(.not.present(frame)) return
      position=[geometry(node)%x_lower,geometry(node)%y_lower,geometry(node)%z_lower]+ &
        (real([i,j,k],dp)-.5_dp)*[geometry(node)%dx,geometry(node)%dy,geometry(node)%dz]
      energy=values(5,q)
      call convert_rotating_energy_3d(frame,position,values(1,q),values(5,q),to_rotating,energy,accepted)
      if(accepted) values(5,q)=energy
    end subroutine boundary_energy
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
  ! Replace interpolated particle ghost sources by same-level particle donors
  ! plus the unmatched deposited remainder. Valid cells retain rates exactly.
  subroutine fill_sdc_mpi_eb_particle_predictor_halo_3d(tree,owners,geometry,rates,particle_rates,particle_remainder, &
      nvar,width,ghost,comm,ok,boundary_faces,target_volume_fraction, &
      ancestor_start,ancestor_end,particle_start,particle_end,ancestor_fraction)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar,width
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: rates,particle_rates,particle_remainder,target_volume_fraction
    type(amr_tree_state_3d), intent(inout) :: ghost
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), intent(in) :: boundary_faces(6)
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,particle_start,particle_end
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    type(amr_tree_state_3d) :: candidate,interpolated,peers
    integer :: flags,low,high,status,node,i,j,k,n(3),offset(3),rank
    real(dp) :: value(nvar)
    logical :: valid,agreed
    ok=.false.
    flags=sum(merge(1,0,[present(ancestor_start),present(ancestor_end),present(particle_start), &
      present(particle_end),present(ancestor_fraction)]))
    call MPI_Allreduce(flags,low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.low/=high.or.(flags/=0.and.flags/=5)) return
    call fill_sdc_mpi_eb_tree_source_halo_3d(tree,owners,geometry,rates,nvar,width,candidate,comm,valid, &
      boundary_faces,target_volume_fraction,ancestor_start,ancestor_end,ancestor_fraction)
    if(.not.valid) return
    call fill_sdc_mpi_eb_tree_source_halo_3d(tree,owners,geometry,particle_rates,nvar,width,interpolated,comm,valid, &
      boundary_faces,target_volume_fraction,particle_start,particle_end,ancestor_fraction)
    if(.not.valid) return
    call fill_sdc_mpi_eb_tree_source_halo_3d(tree,owners,geometry,particle_rates,nvar,width,peers,comm,valid, &
      boundary_faces,target_volume_fraction,same_level_only=.true.)
    if(.not.valid) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    valid=allocated(particle_remainder%fields)
    if(valid) valid=size(particle_remainder%fields)==size(owners)
    call agree();if(.not.agreed) return
    do node=1,size(owners)
      valid=allocated(particle_remainder%fields(node)%values).eqv.(owners(node)==rank)
      if(.not.valid) exit
      if(owners(node)/=rank) cycle
      valid=all(shape(particle_remainder%fields(node)%values)==shape(candidate%fields(node)%values))
      if(.not.valid) exit
      valid=all(ieee_is_finite(particle_remainder%fields(node)%values))
      if(.not.valid) exit
      n=tree%nodes(node)%shape
      offset=[lbound(particle_remainder%fields(node)%values,2),lbound(particle_remainder%fields(node)%values,3), &
        lbound(particle_remainder%fields(node)%values,4)]-1
      valid=all(particle_remainder%fields(node)%values(:,width+1+offset(1):width+n(1)+offset(1), &
        width+1+offset(2):width+n(2)+offset(2),width+1+offset(3):width+n(3)+offset(3))==0)
      if(.not.valid) exit
      do k=1,n(3)+2*width
      do j=1,n(2)+2*width
      do i=1,n(1)+2*width
        if(all([i,j,k]>width).and.all([i,j,k]<=n+width)) cycle
        block
          integer :: loc(4)
          loc=lbound(target_volume_fraction%fields(node)%values)+[0,i-1,j-1,k-1]
          if(target_volume_fraction%fields(node)%values(loc(1),loc(2),loc(3),loc(4))<=0) then
            valid=valid.and.all(particle_remainder%fields(node)%values(:,i+offset(1),j+offset(2),k+offset(3))==0)
            cycle
          end if
        end block
        value=candidate%fields(node)%values(:,i,j,k)
        call add(-interpolated%fields(node)%values(:,i,j,k))
        call add(peers%fields(node)%values(:,i,j,k))
        call add(particle_remainder%fields(node)%values(:,i+offset(1),j+offset(2),k+offset(3)))
        if(valid) candidate%fields(node)%values(:,i,j,k)=value
      end do
      end do
      end do
      if(.not.valid) exit
    end do
    call agree();if(.not.agreed) return
    call move_alloc(candidate%fields,ghost%fields);ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS) agreed=.false.
    end subroutine agree
    subroutine add(increment)
      real(dp), intent(in) :: increment(:)
      if(.not.valid) return
      valid=.not.any(value>huge(1._dp)-max(increment,0._dp)).and. &
        .not.any(value< -huge(1._dp)-min(increment,0._dp))
      if(valid) value=value+increment
    end subroutine add
  end subroutine fill_sdc_mpi_eb_particle_predictor_halo_3d

  ! Signed sources use scalar extension with mirrored wall donor positions;
  ! component values have no EOS recovery or momentum sign flip.
  ! Only fluid donors needed by owned targets are exchanged; covered values may be NaN.
  ! same_level_only leaves physical exterior and missing same-level donors zero;
  ! it cannot be combined with ancestor interpolation. Default behavior is unchanged.
  subroutine fill_sdc_mpi_eb_tree_source_halo_3d(tree,owners,geometry,rates,nvar,halo_width,ghost,comm,ok, &
      boundary_faces,target_volume_fraction,ancestor_start,ancestor_end,ancestor_fraction,same_level_only)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar,halo_width
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: rates
    type(amr_tree_state_3d), intent(inout) :: ghost
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(amr_tree_state_3d), optional, intent(in) :: target_volume_fraction
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end
    ! One fraction per ancestor depth, root first; length max tree level.
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    logical, optional, intent(in) :: same_level_only
    logical :: peers
    type(amr_tree_state_3d) :: candidate
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:),requests(:,:),destination(:,:),selected(:),queries(:,:)
    real(dp), allocatable :: values(:,:),fraction(:),fraction_low(:),fraction_high(:)
    integer(int64) :: global(3),mapped(3),local(3),maximum
    integer :: rank,status,ios,node,target,n(3),level,depth,donor,source(3),idx(3),i,j,k,axis,q
    integer :: boundary_face
    integer :: controls(5),low(5),high(5),width
    character(len=32) :: faces(6),reference(6)
    real(dp) :: volume
    logical :: valid,ready
    ok=.false.;width=halo_width;peers=.false.
    if(present(same_level_only)) peers=same_level_only
    controls=[nvar,width,merge(1,0,present(target_volume_fraction)), &
      sum(merge(1,0,[present(ancestor_start),present(ancestor_end),present(ancestor_fraction)])),merge(1,0,peers)]
    call MPI_Allreduce(controls,low,5,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,5,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    if(nvar<1.or.width<1.or.width>7) return
    if(controls(4)/=0.and.controls(4)/=3) return
    if(peers.and.controls(4)/=0) return
    faces='outflow'
    if(present(boundary_faces)) faces=boundary_faces
    reference=faces
    call MPI_Bcast(reference,len(reference)*6,MPI_CHARACTER,0,comm,status)
    call agree(status==MPI_SUCCESS.and.all(faces==reference).and.valid_reactive_boundaries_3d(faces),ready)
    if(.not.ready) return
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nvar,comm,ready)
    if(.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,rates,nvar,comm,ready)
    if(.not.ready) return
    call agree(periodic_eb_faces_match(geometry(1),faces),ready)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid,ready)
    if(.not.ready) return
    if(present(ancestor_fraction)) then
      call agree(size(ancestor_fraction)==size(ratios)-1,ready)
      if(.not.ready) return
      call agree(all(ieee_is_finite(ancestor_fraction)),ready)
      if(.not.ready) return
      call agree(all(ancestor_fraction>=0).and.all(ancestor_fraction<=1),ready)
      if(.not.ready) return
      allocate(fraction_low(size(ancestor_fraction)),fraction_high(size(ancestor_fraction)),stat=ios)
      call agree(ios==0,ready)
      if(.not.ready) return
      call MPI_Allreduce(ancestor_fraction,fraction_low,size(ancestor_fraction),MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(ancestor_fraction,fraction_high,size(ancestor_fraction),MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.any(fraction_low/=fraction_high)) return
      call validate_mpi_amr_tree_local_3d(tree,owners,ancestor_start,nvar,comm,ready)
      if(.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,ancestor_end,nvar,comm,ready)
      if(.not.ready) return
    end if
    valid=.true.
    if(present(target_volume_fraction)) then
      valid=allocated(target_volume_fraction%fields)
      if(valid) valid=size(target_volume_fraction%fields)==size(owners)
    end if
    call agree(valid,ready)
    if(.not.ready) return
    maximum=0
    do node=1,size(owners)
      n=tree%nodes(node)%shape
      if(any(n>huge(1)-2*width)) then
        valid=.false.;exit
      end if
      if(owners(node)/=rank) cycle
      maximum=maximum+product(int(n+2*width,int64))
      if(present(target_volume_fraction)) then
        if(.not.allocated(target_volume_fraction%fields(node)%values)) then
          valid=.false.;exit
        end if
        associate(v=>target_volume_fraction%fields(node)%values)
        if(any(shape(v)/=[1,n+2*width])) then
          valid=.false.;exit
        end if
        if(.not.all(ieee_is_finite(v)).or.any(v<0).or.any(v>1)) then
          valid=.false.;exit
        end if
        ! Use array positions, independent of caller lower bounds.
        end associate
      end if
    end do
    call agree(valid.and.maximum<=int(huge(1)/max(4,nvar),int64),ready)
    if(.not.ready) return
    allocate(candidate%fields(size(owners)),requests(4,int(maximum)),destination(4,int(maximum)), &
      selected(int(maximum)),queries(4,int(maximum)),fraction(int(maximum)),stat=ios)
    call agree(ios==0,ready)
    if(.not.ready) return
    q=0;fraction=-1
    targets: do target=1,size(owners)
      if(owners(target)/=rank) cycle
      n=tree%nodes(target)%shape;level=tree%nodes(target)%level+1
      allocate(candidate%fields(target)%values(nvar,n(1)+2*width,n(2)+2*width,n(3)+2*width),source=0._dp,stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      end if
      do k=1,n(3)+2*width
      do j=1,n(2)+2*width
      do i=1,n(1)+2*width
        idx=[i,j,k]-width
        global=origin(:,target)+int(idx-1,int64)
        do axis=1,3
          if(faces(2*axis-1)=='periodic') global(axis)=modulo(global(axis),domain(axis,level))
        end do
        if(peers.and.(any(global<0).or.any(global>=domain(:,level)))) cycle
        ! Match state/geometry donor locations at physical walls. Source
        ! components retain their existing scalar extension convention.
        do axis=1,3
          do while(global(axis)<0.or.global(axis)>=domain(axis,level))
            boundary_face=2*axis
            if(global(axis)<0) boundary_face=boundary_face-1
            if(faces(boundary_face)/='slip_wall'.and.faces(boundary_face)/='no_slip_wall'.and. &
              faces(boundary_face)/='isothermal_wall') exit
            if(global(axis)<0) then
              global(axis)=-global(axis)-1
            else
              global(axis)=2*domain(axis,level)-global(axis)-1
            end if
          end do
        end do
        global=max(0_int64,min(domain(:,level)-1,global))
        donor=0
        if(all(idx>=1).and.all(idx<=n)) then
          donor=target;source=idx
        else
          do depth=level,merge(level,1,peers),-1
            mapped=global/(domain(:,level)/domain(:,depth))
            do node=1,size(owners)
              if(tree%nodes(node)%level+1/=depth) cycle
              local=mapped-origin(:,node)
              if(any(local<0).or.any(local>=int(tree%nodes(node)%shape,int64))) cycle
              donor=node;source=int(local)+1;exit
            end do
            if(donor/=0) exit
          end do
        end if
        if(donor==0.and.peers) cycle
        if(donor==0) then
          valid=.false.;exit targets
        end if
        volume=geometry(donor)%volume_fraction(source(1),source(2),source(3))
        if(present(target_volume_fraction)) then
          block
            integer :: loc(4)
            loc=lbound(target_volume_fraction%fields(target)%values)+[0,i-1,j-1,k-1]
            associate(v=>target_volume_fraction%fields(target)%values)
            if(all(idx>=1).and.all(idx<=n)) then
              if(v(loc(1),loc(2),loc(3),loc(4))/=volume) then
                valid=.false.;exit targets
              end if
            end if
            if(v(loc(1),loc(2),loc(3),loc(4))<=0) cycle
            end associate
          end block
          if(volume<=0) then
            valid=.false.;exit targets
          end if
        else
          if(volume<=0) cycle
        end if
        q=q+1;requests(:,q)=[donor,source];destination(:,q)=[target,i,j,k]
        depth=tree%nodes(donor)%level+1
        if(present(ancestor_fraction).and.depth<level) fraction(q)=ancestor_fraction(depth)
      end do
      end do
      end do
    end do targets
    call agree(valid,ready)
    if(.not.ready) return
    call fetch_sources(rates,0,ready)
    if(.not.ready) return
    if(present(ancestor_fraction)) then
      call fetch_sources(ancestor_start,1,ready)
      if(.not.ready) return
      call fetch_sources(ancestor_end,2,ready)
      if(.not.ready) return
    end if
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      if(.not.all(ieee_is_finite(candidate%fields(node)%values))) valid=.false.
    end do
    call agree(valid,ready)
    if(.not.ready) return
    call move_alloc(candidate%fields,ghost%fields);ok=.true.
  contains
    subroutine fetch_sources(field,phase,success)
      type(amr_tree_state_3d), intent(in) :: field
      integer, intent(in) :: phase
      logical, intent(out) :: success
      integer :: query,used,slot,p,cell(3)
      real(dp) :: weight
      used=0
      do query=1,q
        select case(phase)
        case(0)
          if(fraction(query)>=0) cycle
        case(1)
          if(fraction(query)<0.or.fraction(query)==1) cycle
        case(2)
          if(fraction(query)<=0) cycle
        end select
        used=used+1;queries(:,used)=requests(:,query);selected(used)=query
      end do
      call fetch_mpi_amr_tree_cells_3d(tree,owners,field,nvar,queries(:,:used),values,comm,success)
      if(.not.success) return
      do slot=1,used
        query=selected(slot);p=destination(1,query);cell=destination(2:4,query)
        if(phase==0) then
          candidate%fields(p)%values(:,cell(1),cell(2),cell(3))=values(:,slot)
        else
          weight=fraction(query)
          if(phase==1) weight=1-weight
          candidate%fields(p)%values(:,cell(1),cell(2),cell(3))= &
            candidate%fields(p)%values(:,cell(1),cell(2),cell(3))+weight*values(:,slot)
        end if
      end do
    end subroutine fetch_sources
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS) result=.false.
    end subroutine
  end subroutine fill_sdc_mpi_eb_tree_source_halo_3d
end module mpi_eb_tree_halo_3d_mod
