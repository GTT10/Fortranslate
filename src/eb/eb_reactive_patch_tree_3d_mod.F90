module eb_reactive_patch_tree_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use simulation_config_reactive_3d_mod, only: mach_inlet_ghost_3d
  use reactive_3d_mod, only: advance_reactive_chemistry_3d
  use gas_transport_mod, only: gas_transport_species
  use eb_reactive_transport_3d_mod, only: advance_reactive_eb_transport_patch_euler_3d, &
    reactive_eb_transport_patch_timestep_3d,reactive_eb_halo_bounds_3d
  use reactive_eb_cfl_3d_mod, only: reactive_eb_signal_speeds_3d,compute_reactive_eb_euler_cfl_timestep_3d
  use eb_reactive_hydro_3d_mod, only: advance_reactive_eb_euler_3d,advance_reactive_eb_plm_ctu_3d, &
    compute_reactive_eb_flattening_3d,eb_plm_boundary_3d
  use eb_reactive_wall_flux_3d_mod, only: reactive_eb_slip_wall_source_3d
  use elementary_kinetics_mod, only: elementary_reaction
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_integer_layout_3d
  use amr_patch_tree_tagging_3d_mod, only: hot_cell,build_tagged_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d,extract_eb_geometry_3d,periodic_eb_faces_match
  use eb_reactive_two_level_3d_mod, only: eb_reactive_patch_3d,matching_eb_metrics_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d,valid_eb_integrated_flux_3d, &
    restrict_eb_integrated_flux_3d,reconcile_reactive_eb_peer_flux_3d,reflux_reactive_eb_state_3d
  use eb_reactive_transfer_3d_mod, only: prolong_reactive_eb_state_pcm_3d,restrict_reactive_eb_state_3d, &
    overlay_reactive_eb_peer_halo_3d
  implicit none
  private
  type, public :: eb_reactive_tree_3d
    type(amr_patch_tree_3d) :: topology
    type(eb_reactive_patch_3d), allocatable :: nodes(:)
  end type
  public :: build_eb_tree_extended_geometry_3d
  public :: tag_reactive_eb_tree_temperature_3d
  public :: remap_reactive_eb_tree_3d
  public :: reactive_eb_tree_timestep_3d
  public :: advance_reactive_eb_tree_full_3d
  public :: advance_reactive_eb_tree_transport_3d
  public :: advance_reactive_eb_tree_subcycle_transport_3d
  public :: advance_reactive_eb_tree_subcycle_hydro_3d
  public :: advance_reactive_eb_tree_hydro_3d
  public :: synchronize_reactive_eb_tree_flux_3d
  public :: fill_reactive_eb_tree_halo_3d
  public :: react_reactive_eb_tree_3d
  public :: initialize_reactive_eb_tree_3d,average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d
contains
  subroutine bounds(tree,node,lo,hi)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: node
    integer, intent(out) :: lo(3),hi(3)
    lo=[tree%nodes(node)%patch%coarse_i_lower,tree%nodes(node)%patch%coarse_j_lower, &
      tree%nodes(node)%patch%coarse_k_lower]
    hi=[tree%nodes(node)%patch%coarse_i_upper,tree%nodes(node)%patch%coarse_j_upper, &
      tree%nodes(node)%patch%coarse_k_upper]
  end subroutine bounds

  ! Build target-resolution margins from a finest domain padded by width times
  ! the total refinement ratio. Geometry remains replicated for MPI callers.
  subroutine build_eb_tree_extended_geometry_3d(topology,physical,finest,patches,ok,halo_width)
    type(amr_patch_tree_3d), intent(in) :: topology
    type(eb_geometry_3d), intent(in) :: physical(:),finest
    type(eb_geometry_3d), allocatable, intent(inout) :: patches(:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: halo_width
    type(eb_geometry_3d), allocatable :: levels(:),candidate(:)
    type(eb_geometry_3d) :: center
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer(int64) :: padding
    integer, allocatable :: ratios(:),margin(:)
    integer :: width,count,level,node,n(3),lo(3),hi(3)
    logical :: valid
    ok=.false.;width=6
    if(present(halo_width)) width=halo_width
    if(width<1.or.width>huge(1)/4) return
    call amr_tree_integer_layout_3d(topology,origin,domain,ratios,valid)
    if(.not.valid.or..not.finest%is_valid()) return
    if(size(physical)/=size(topology%nodes)) return
    count=size(ratios)
    padding=domain(1,count)/domain(1,1)
    if(padding>int(huge(1)/2,int64)/width) return
    padding=padding*width
    if(any(domain(:,count)>int(huge(1),int64)-2*padding)) return
    if(any(int([finest%nx,finest%ny,finest%nz],int64)/=domain(:,count)+2*padding)) return
    allocate(levels(count),margin(count),candidate(size(physical)))
    levels(count)=finest;margin(count)=int(padding)
    do level=count,2,-1
      if(mod(margin(level),ratios(level))/=0) return
      margin(level-1)=margin(level)/ratios(level)
      call coarsen_eb_geometry_3d(levels(level),ratios(level),levels(level-1),valid)
      if(.not.valid) return
    end do
    do node=1,size(physical)
      if(.not.physical(node)%is_valid()) return
      level=topology%nodes(node)%level+1;n=topology%nodes(node)%shape
      if(any([physical(node)%nx,physical(node)%ny,physical(node)%nz]/=n)) return
      lo=int(origin(:,node))+1+margin(level)-width;hi=lo+n+2*width-1
      call extract_eb_geometry_3d(levels(level),lo,hi,candidate(node),valid)
      if(.not.valid) return
      call extract_eb_geometry_3d(candidate(node),[width+1,width+1,width+1],n+width,center,valid)
      if(.not.valid) return
      if(.not.matching_eb_metrics_3d(center,physical(node))) return
      ! Coordinates used during clipping can differ by roundoff after padding.
      ! Keep the accepted physical patch metrics bit-for-bit, including its faces.
      associate(g=>candidate(node),p=>physical(node))
      g%volume_fraction(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%volume_fraction
      g%cell_centroid_x(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%cell_centroid_x
      g%cell_centroid_y(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%cell_centroid_y
      g%cell_centroid_z(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%cell_centroid_z
      g%boundary_area(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%boundary_area
      g%boundary_centroid_x(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%boundary_centroid_x
      g%boundary_centroid_y(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%boundary_centroid_y
      g%boundary_centroid_z(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%boundary_centroid_z
      g%boundary_normal_x(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%boundary_normal_x
      g%boundary_normal_y(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%boundary_normal_y
      g%boundary_normal_z(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%boundary_normal_z
      g%boundary_normal_integral_x(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%boundary_normal_integral_x
      g%boundary_normal_integral_y(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%boundary_normal_integral_y
      g%boundary_normal_integral_z(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%boundary_normal_integral_z
      g%cell_type(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%cell_type
      g%x_face_fraction(width:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%x_face_fraction
      g%x_face_centroid_y(width:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%x_face_centroid_y
      g%x_face_centroid_z(width:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%x_face_centroid_z
      g%y_face_fraction(width+1:width+n(1),width:width+n(2),width+1:width+n(3))=p%y_face_fraction
      g%y_face_centroid_x(width+1:width+n(1),width:width+n(2),width+1:width+n(3))=p%y_face_centroid_x
      g%y_face_centroid_z(width+1:width+n(1),width:width+n(2),width+1:width+n(3))=p%y_face_centroid_z
      g%z_face_fraction(width+1:width+n(1),width+1:width+n(2),width:width+n(3))=p%z_face_fraction
      g%z_face_centroid_x(width+1:width+n(1),width+1:width+n(2),width:width+n(3))=p%z_face_centroid_x
      g%z_face_centroid_y(width+1:width+n(1),width+1:width+n(2),width:width+n(3))=p%z_face_centroid_y
      if(allocated(g%halo_volume_fraction)) g%halo_volume_fraction(1:g%nx,1:g%ny,1:g%nz)=g%volume_fraction
      end associate
      if(.not.candidate(node)%is_valid()) return
    end do
    call move_alloc(candidate,patches);ok=.true.
  end subroutine build_eb_tree_extended_geometry_3d

  subroutine initialize_reactive_eb_tree_3d(species,topology,finest,root_state,root_temperature,tree,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: topology
    type(eb_geometry_3d), intent(in) :: finest
    real(dp), intent(in) :: root_state(:,:,:,:),root_temperature(:,:,:)
    type(eb_reactive_tree_3d), intent(inout) :: tree
    logical, intent(out) :: ok
    type(eb_reactive_tree_3d) :: candidate
    type(eb_geometry_3d), allocatable :: levels(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer :: count,level,node,parent,ratio,lo(3),hi(3),d(3),i,j,k,index(3),nv
    logical :: valid
    ok=.false.
    call amr_tree_integer_layout_3d(topology,origin,domain,ratios,valid)
    if (.not.valid.or..not.finest%is_valid()) return
    count=size(ratios);nv=reactive_nvar(size(species))
    if (any(int([finest%nx,finest%ny,finest%nz],int64)/=domain(:,count))) return
    if (any(shape(root_state)/=[nv,topology%nodes(1)%shape])) return
    if (any(shape(root_temperature)/=topology%nodes(1)%shape)) return
    allocate(levels(count));levels(count)=finest
    do level=count,2,-1
      call coarsen_eb_geometry_3d(levels(level),ratios(level),levels(level-1),valid)
      if (.not.valid) return
    end do
    candidate%topology=topology
    allocate(candidate%nodes(size(topology%nodes)))
    do node=1,size(topology%nodes)
      level=topology%nodes(node)%level+1;d=topology%nodes(node)%shape
      lo=int(origin(:,node))+1;hi=lo+d-1
      call extract_eb_geometry_3d(levels(level),lo,hi,candidate%nodes(node)%geometry,valid)
      if (.not.valid) return
      allocate(candidate%nodes(node)%state(nv,d(1),d(2),d(3)),candidate%nodes(node)%temperature(d(1),d(2),d(3)))
      if (node==1) then
        candidate%nodes(node)%state=root_state;candidate%nodes(node)%temperature=root_temperature
        call prolong_reactive_eb_state_pcm_3d(species,levels(1),1,root_state,root_temperature, &
          candidate%nodes(1)%state,candidate%nodes(1)%temperature,valid)
      else
        parent=topology%nodes(node)%parent;ratio=topology%nodes(node)%patch%refinement_ratio
        call bounds(topology,node,lo,hi)
        do k=1,d(3)
        do j=1,d(2)
        do i=1,d(1)
          index=lo+([i,j,k]-1)/ratio
          candidate%nodes(node)%state(:,i,j,k)=candidate%nodes(parent)%state(:,index(1),index(2),index(3))
          candidate%nodes(node)%temperature(i,j,k)=candidate%nodes(parent)%temperature(index(1),index(2),index(3))
        end do
        end do
        end do
        call prolong_reactive_eb_state_pcm_3d(species,candidate%nodes(node)%geometry,ratio, &
          candidate%nodes(parent)%state(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
          candidate%nodes(parent)%temperature(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
          candidate%nodes(node)%state,candidate%nodes(node)%temperature,valid)
      end if
      if (.not.valid) return
    end do
    if (.not.valid_tree(candidate,nv)) return
    tree=candidate;ok=.true.
  end subroutine initialize_reactive_eb_tree_3d

  ! Exposed fluid cells vote without averaging away fine temperature peaks.
  ! Reuse regular AMR's deterministic tile builder, then omit solid-only nodes.
  subroutine tag_reactive_eb_tree_temperature_3d(tree,finest,thresholds,ratios,tile_shape,buffer, &
      maximum_nodes,topology,ok)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: finest
    real(dp), intent(in) :: thresholds(:)
    integer, intent(in) :: ratios(:),tile_shape(3),buffer,maximum_nodes
    type(amr_patch_tree_3d), intent(inout) :: topology
    logical, intent(out) :: ok
    type(hot_cell), allocatable :: cells(:)
    type(amr_patch_tree_3d) :: candidate,pruned
    type(eb_geometry_3d), allocatable :: levels(:)
    type(eb_geometry_3d) :: section
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: old_ratios(:),shared(:),mapping(:)
    integer(int64) :: common(3),factor(3),total
    logical, allocatable :: active(:,:,:)
    integer :: nv,node,child,level,count,pass,i,j,k,n(3),lo(3),hi(3),used
    real(dp) :: temperature
    logical :: valid
    ok=.false.
    if (.not.allocated(tree%nodes)) return
    if (size(tree%nodes)<1) return
    if (.not.allocated(tree%nodes(1)%state)) return
    nv=size(tree%nodes(1)%state,1)
    if (.not.valid_tree(tree,nv).or..not.finest%is_valid()) return
    if (size(thresholds)/=size(ratios)) return
    if (.not.all(ieee_is_finite(thresholds)).or.any(thresholds<=0._dp)) return
    if (any(ratios<2).or.any(tile_shape<1).or.buffer<0.or.maximum_nodes<1) return
    do level=2,size(thresholds)
      if (thresholds(level)<thresholds(level-1)) return
    end do
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,old_ratios,valid)
    if (.not.valid) return
    level=min(size(old_ratios)-1,size(ratios))
    if (any(old_ratios(2:level+1)/=ratios(:level))) return
    count=max(size(old_ratios),size(ratios)+1)
    allocate(shared(count));shared(:size(old_ratios))=old_ratios
    shared(2:size(ratios)+1)=ratios;common=domain(:,1)
    do level=2,count
      if (any(common>(huge(1_int64)/4)/int(shared(level),int64))) return
      common=common*int(shared(level),int64)
    end do
    if (any(common/=int([finest%nx,finest%ny,finest%nz],int64))) return
    allocate(levels(count));levels(count)=finest
    do level=count,2,-1
      call coarsen_eb_geometry_3d(levels(level),shared(level),levels(level-1),valid)
      if (.not.valid) return
    end do
    do node=1,size(tree%nodes)
      level=tree%topology%nodes(node)%level+1
      lo=int(origin(:,node))+1;hi=lo+tree%topology%nodes(node)%shape-1
      call extract_eb_geometry_3d(levels(level),lo,hi,section,valid)
      if (.not.valid) return
      if (.not.matching_eb_metrics_3d(section,tree%nodes(node)%geometry)) return
    end do
    total=0
    do pass=1,2
      used=0
      do node=1,size(tree%nodes)
        n=tree%topology%nodes(node)%shape;level=tree%topology%nodes(node)%level+1
        factor=common/domain(:,level)
        allocate(active(n(1),n(2),n(3)));active=tree%nodes(node)%geometry%volume_fraction>0._dp
        do child=node+1,size(tree%nodes)
          if (tree%topology%nodes(child)%parent/=node) cycle
          call bounds(tree%topology,child,lo,hi)
          active(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.false.
        end do
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if (.not.active(i,j,k)) cycle
          temperature=tree%nodes(node)%temperature(i,j,k)
          if (.not.ieee_is_finite(temperature).or.temperature<=0._dp) return
          if (size(thresholds)==0) cycle
          if (temperature<thresholds(1)) cycle
          if (pass==1) then
            total=total+1
            if (total>int(huge(1),int64)) return
          else
            used=used+1
            cells(used)%lower=(origin(:,node)+int([i,j,k]-1,int64))*factor
            cells(used)%upper=cells(used)%lower+factor;cells(used)%temperature=temperature
          end if
        end do
        end do
        end do
        deallocate(active)
      end do
      if (pass==1) allocate(cells(int(total)))
    end do
    call build_tagged_tree_3d(tree%topology,cells,common,thresholds,ratios,tile_shape,buffer, &
      maximum_nodes,[.false.,.false.,.false.],candidate,valid)
    if (.not.valid) return
    call amr_tree_integer_layout_3d(candidate,origin,domain,old_ratios,valid)
    if (.not.valid) return
    allocate(mapping(size(candidate%nodes)));mapping=0
    pruned=candidate;used=1;mapping(1)=1
    do node=2,size(candidate%nodes)
      child=candidate%nodes(node)%parent
      if (mapping(child)==0) cycle
      level=candidate%nodes(node)%level+1
      lo=int(origin(:,node))+1;hi=lo+candidate%nodes(node)%shape-1
      call extract_eb_geometry_3d(levels(level),lo,hi,section,valid)
      if (.not.valid) return
      if (.not.any(section%volume_fraction>0._dp)) cycle
      used=used+1;mapping(node)=used
      pruned%nodes(used)=candidate%nodes(node);pruned%nodes(used)%parent=mapping(child)
    end do
    pruned%nodes=pruned%nodes(:used)
    if (.not.pruned%is_valid()) return
    topology=pruned;ok=.true.
  end subroutine tag_reactive_eb_tree_temperature_3d

  ! Conservative remap on a shared refinement lattice. finest describes the
  ! deeper of old/new trees; ratios on shared levels and root cells must match.
  ! Averaged old cells supply coarsening, exact same-level overlap or ancestor
  ! PCM supplies new refinement. Geometry is immutable wherever old nodes exist.
  subroutine remap_reactive_eb_tree_3d(species,topology,finest,tree,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: topology
    type(eb_geometry_3d), intent(in) :: finest
    type(eb_reactive_tree_3d), intent(inout) :: tree
    logical, intent(out) :: ok
    type(eb_reactive_tree_3d) :: snapshot,candidate,expected
    type(eb_geometry_3d), allocatable :: levels(:)
    integer(int64), allocatable :: old_origin(:,:),old_domain(:,:),new_origin(:,:),new_domain(:,:)
    integer, allocatable :: old_ratios(:),new_ratios(:),ratios(:)
    integer(int64) :: global(3),mapped(3),local(3)
    integer :: nv,nlevels,old_levels,new_levels,level,node,donor,depth,i,j,k,d(3),index(3)
    real(dp) :: primitive(reactive_nprim(size(species))),temperature,speed
    logical :: valid,found
    ok=.false.;nv=reactive_nvar(size(species))
    if (.not.valid_tree(tree,nv).or..not.finest%is_valid()) return
    call amr_tree_integer_layout_3d(tree%topology,old_origin,old_domain,old_ratios,valid)
    if (.not.valid) return
    call amr_tree_integer_layout_3d(topology,new_origin,new_domain,new_ratios,valid)
    if (.not.valid) return
    old_levels=size(old_ratios);new_levels=size(new_ratios);nlevels=max(old_levels,new_levels)
    if (any(old_domain(:,1)/=new_domain(:,1))) return
    level=min(old_levels,new_levels)
    if (any(old_ratios(:level)/=new_ratios(:level))) return
    if (old_levels>=new_levels) then
      ratios=old_ratios
      if (any(int([finest%nx,finest%ny,finest%nz],int64)/=old_domain(:,nlevels))) return
    else
      ratios=new_ratios
      if (any(int([finest%nx,finest%ny,finest%nz],int64)/=new_domain(:,nlevels))) return
    end if
    allocate(levels(nlevels));levels(nlevels)=finest
    do level=nlevels,2,-1
      call coarsen_eb_geometry_3d(levels(level),ratios(level),levels(level-1),valid)
      if (.not.valid) return
    end do
    snapshot=tree
    call average_down_reactive_eb_tree_3d(species,snapshot,valid)
    if (.not.valid) return
    call initialize_reactive_eb_tree_3d(species,tree%topology,levels(old_levels),snapshot%nodes(1)%state, &
      snapshot%nodes(1)%temperature,expected,valid)
    if (.not.valid) return
    do node=1,size(tree%nodes)
      if (.not.matching_eb_metrics_3d(tree%nodes(node)%geometry,expected%nodes(node)%geometry)) return
    end do
    call initialize_reactive_eb_tree_3d(species,topology,levels(new_levels),snapshot%nodes(1)%state, &
      snapshot%nodes(1)%temperature,candidate,valid)
    if (.not.valid) return
    do node=1,size(candidate%nodes)
      level=topology%nodes(node)%level+1;d=topology%nodes(node)%shape
      do k=1,d(3)
      do j=1,d(2)
      do i=1,d(1)
        global=new_origin(:,node)+int([i,j,k]-1,int64);found=.false.
        do depth=min(level,old_levels),1,-1
          mapped=global/(new_domain(:,level)/old_domain(:,depth))
          do donor=1,size(snapshot%nodes)
            if (snapshot%topology%nodes(donor)%level+1/=depth) cycle
            local=mapped-old_origin(:,donor)+1_int64
            if (any(local<1_int64).or.any(local>int(snapshot%topology%nodes(donor)%shape,int64))) cycle
            index=int(local)
            if (candidate%nodes(node)%geometry%volume_fraction(i,j,k)>0._dp) then
              if (snapshot%nodes(donor)%geometry%volume_fraction(index(1),index(2),index(3))<=0._dp) return
              call reactive_conserved_to_primitive(species,snapshot%nodes(donor)%state(:,index(1),index(2),index(3)), &
                snapshot%nodes(donor)%temperature(index(1),index(2),index(3)),primitive,temperature,speed,valid)
              if (.not.valid) return
            end if
            candidate%nodes(node)%state(:,i,j,k)=snapshot%nodes(donor)%state(:,index(1),index(2),index(3))
            candidate%nodes(node)%temperature(i,j,k)=snapshot%nodes(donor)%temperature(index(1),index(2),index(3))
            found=.true.;exit
          end do
          if (found) exit
        end do
        if (.not.found) return
      end do
      end do
      end do
    end do
    call average_down_reactive_eb_tree_3d(species,candidate,valid)
    if (.not.valid) return
    tree=candidate;ok=.true.
  end subroutine remap_reactive_eb_tree_3d

  logical function valid_tree(tree,nv) result(valid)
    type(eb_reactive_tree_3d), intent(in) :: tree
    integer, intent(in) :: nv
    type(eb_geometry_3d) :: section,restricted
    integer :: node,parent,lo(3),hi(3),d(3)
    logical :: good
    valid=.false.
    if (.not.tree%topology%is_valid().or..not.allocated(tree%nodes)) return
    if (size(tree%nodes)/=size(tree%topology%nodes)) return
    do node=1,size(tree%nodes)
      if (.not.tree%nodes(node)%geometry%is_valid()) return
      if (.not.allocated(tree%nodes(node)%state).or..not.allocated(tree%nodes(node)%temperature)) return
      d=tree%topology%nodes(node)%shape
      if (any([tree%nodes(node)%geometry%nx,tree%nodes(node)%geometry%ny,tree%nodes(node)%geometry%nz]/=d)) return
      if (any(shape(tree%nodes(node)%state)/=[nv,d]).or.any(shape(tree%nodes(node)%temperature)/=d)) return
      if (node==1) cycle
      parent=tree%topology%nodes(node)%parent
      call bounds(tree%topology,node,lo,hi)
      call coarsen_eb_geometry_3d(tree%nodes(node)%geometry,tree%topology%nodes(node)%patch%refinement_ratio,restricted,good)
      if (.not.good) return
      call extract_eb_geometry_3d(tree%nodes(parent)%geometry,lo,hi,section,good)
      if (.not.good) return
      if (.not.matching_eb_metrics_3d(section,restricted)) return
    end do
    valid=.true.
  end function valid_tree

  subroutine average_down_reactive_eb_tree_3d(species,tree,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_reactive_tree_3d), intent(inout) :: tree
    logical, intent(out) :: ok
    type(eb_reactive_tree_3d) :: candidate
    integer :: node,parent,lo(3),hi(3)
    logical :: valid
    ok=.false.
    if (.not.valid_tree(tree,reactive_nvar(size(species)))) return
    candidate=tree
    do node=size(tree%nodes),2,-1
      parent=tree%topology%nodes(node)%parent
      call bounds(tree%topology,node,lo,hi)
      call restrict_reactive_eb_state_3d(species,candidate%nodes(node)%geometry, &
        tree%topology%nodes(node)%patch%refinement_ratio,candidate%nodes(node)%state,candidate%nodes(node)%temperature, &
        candidate%nodes(parent)%state(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
        candidate%nodes(parent)%temperature(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),valid)
      if (.not.valid) return
    end do
    tree=candidate;ok=.true.
  end subroutine average_down_reactive_eb_tree_3d

  ! Bound every patch update, including hidden parent cells which are advanced
  ! before restriction. Global donor maxima cover peer and ancestor halos.
  subroutine reactive_eb_tree_timestep_3d(species,transport,tree,hydro_cfl,transport_cfl, &
      transport_enabled,viscosity,conduction,diffusion,dt,ok,subcycling, &
      inflow_mach,inflow_temperature,inflow_pressure,boundary_faces,inflow_state)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(eb_reactive_tree_3d), intent(in) :: tree
    real(dp), intent(in) :: hydro_cfl,transport_cfl
    logical, intent(in) :: transport_enabled,viscosity,conduction,diffusion
    real(dp), intent(inout) :: dt
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: subcycling
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    type(eb_reactive_tree_3d) :: snapshot
    real(dp) :: speeds(4),donors(4),diffusivity,maximum_diffusivity,interval,candidate,factor
    real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:)
    integer :: node,level,n(3)
    logical :: valid,active,recursive
    ok=.false.
    if (.not.all(ieee_is_finite([hydro_cfl,transport_cfl]))) return
    if (hydro_cfl<=0._dp.or.hydro_cfl>1._dp.or.transport_cfl<=0._dp.or.transport_cfl>.5_dp) return
    if (.not.valid_tree(tree,reactive_nvar(size(species)))) return
    recursive=.false.
    if (present(subcycling)) recursive=subcycling
    if (recursive) then
      call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
      if (.not.valid) return
    end if
    snapshot=tree
    call average_down_reactive_eb_tree_3d(species,snapshot,valid)
    if (.not.valid) return
    donors=0._dp;maximum_diffusivity=0._dp;candidate=huge(1._dp);active=.false.
    do node=1,size(tree%nodes)
      associate(g=>snapshot%nodes(node)%geometry,u=>snapshot%nodes(node)%state,t=>snapshot%nodes(node)%temperature)
      if (.not.any(g%volume_fraction>0._dp)) cycle
      active=.true.
      call reactive_eb_signal_speeds_3d(species,u,t,g,speeds,valid)
      if (.not.valid) return
      donors=max(donors,speeds)
      if (present(inflow_mach).or.present(inflow_temperature).or.present(inflow_pressure).or. &
        present(boundary_faces).or.present(inflow_state)) then
        n=tree%topology%nodes(node)%shape
        allocate(ghost(size(u,1),0:n(1)+1,0:n(2)+1,0:n(3)+1),gt(0:n(1)+1,0:n(2)+1,0:n(3)+1))
        call fill_reactive_eb_tree_halo_3d(species,snapshot,node,ghost,gt,valid, &
          inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
          boundary_faces=boundary_faces,inflow_state=inflow_state)
        if (.not.valid) return
        call reactive_eb_halo_bounds_3d(species,transport,g,ghost,gt,transport_enabled, &
          viscosity,conduction,diffusion,speeds,diffusivity,valid)
        if (.not.valid) return
        donors=max(donors,speeds);maximum_diffusivity=max(maximum_diffusivity,diffusivity)
        deallocate(ghost,gt)
      end if
      if (transport_enabled) then
        call reactive_eb_transport_patch_timestep_3d(species,transport,u,t,g,transport_cfl, &
          viscosity,conduction,diffusion,interval,diffusivity,valid)
        if (.not.valid) return
        maximum_diffusivity=max(maximum_diffusivity,diffusivity)
      end if
      end associate
    end do
    if (.not.active) return
    do node=1,size(tree%nodes)
      associate(g=>snapshot%nodes(node)%geometry,u=>snapshot%nodes(node)%state,t=>snapshot%nodes(node)%temperature)
      if (.not.any(g%volume_fraction>0._dp)) cycle
      call compute_reactive_eb_euler_cfl_timestep_3d(species,u,t,g,hydro_cfl,interval,valid,donors)
      if (.not.valid) return
      factor=1._dp
      if (recursive) then
        level=tree%topology%nodes(node)%level+1
        factor=real(domain(1,level),dp)/real(domain(1,1),dp)
      end if
      ! Compare before multiplying so an inactive huge bound cannot overflow.
      if (interval<candidate/factor) candidate=interval*factor
      if (transport_enabled) then
        call reactive_eb_transport_patch_timestep_3d(species,transport,u,t,g,transport_cfl, &
          viscosity,conduction,diffusion,interval,diffusivity,valid,maximum_diffusivity)
        if (.not.valid) return
        ! Conservatively use a full-interval diffusion bound for split halves.
        factor=factor*factor
        if (interval<candidate/factor) candidate=interval*factor
      end if
      end associate
    end do
    dt=candidate;ok=.true.
  end subroutine reactive_eb_tree_timestep_3d

  ! R/2 T/2 H T/2 R/2 composition. Each successful stage remains private
  ! until the complete step and its exchange vectors can be committed.
  subroutine advance_reactive_eb_tree_full_3d(species,reactions,transport,solver,dt,chemistry_enabled, &
      rtol,atol,transport_enabled,viscosity,conduction,diffusion,barodiffusion,tree,boundary_change, &
      wall_change,chemical_change,minimum_theta,ok,integrator,subcycling, &
      inflow_mach,inflow_temperature,inflow_pressure,difmag,extended_geometry,boundary_faces,inflow_state)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), optional, intent(in) :: extended_geometry(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt,rtol,atol
    logical, intent(in) :: chemistry_enabled,transport_enabled,viscosity,conduction,diffusion,barodiffusion
    type(eb_reactive_tree_3d), intent(inout) :: tree
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemical_change(:),minimum_theta
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: difmag
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: integrator
    logical, optional, intent(in) :: subcycling
    type(eb_reactive_tree_3d) :: candidate
    real(dp), allocatable :: boundary(:),wall(:),chemical(:),exchange(:),transport_sum(:)
    real(dp) :: theta,stage_theta
    integer :: nv,half
    logical :: valid,recursive_steps
    ok=.false.;nv=reactive_nvar(size(species))
    if (.not.ieee_is_finite(dt).or.dt<=0._dp) return
    if (size(boundary_change)/=nv.or.size(wall_change)/=nv.or.size(chemical_change)/=nv) return
    if (.not.valid_tree(tree,nv)) return
    recursive_steps=.false.
    if (present(subcycling)) recursive_steps=subcycling
    if(present(difmag)) then
      if(.not.ieee_is_finite(difmag).or.difmag<0) return
    end if
    candidate=tree
    allocate(boundary(nv),wall(nv),chemical(nv),exchange(nv),transport_sum(nv))
    chemical=0._dp;transport_sum=0._dp;theta=1._dp
    do half=1,2
      if (half==1.and.chemistry_enabled) then
        call reaction_half(valid)
        if (.not.valid) return
      end if
      if (transport_enabled) then
        if (recursive_steps) then
          call advance_reactive_eb_tree_subcycle_transport_3d(species,transport,0.5_dp*dt,viscosity,conduction, &
            diffusion,barodiffusion,candidate,exchange,stage_theta,valid, &
            inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
            boundary_faces=boundary_faces,inflow_state=inflow_state)
        else
          call advance_reactive_eb_tree_transport_3d(species,transport,0.5_dp*dt,viscosity,conduction, &
            diffusion,barodiffusion,candidate,exchange,stage_theta,valid, &
            inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
            boundary_faces=boundary_faces,inflow_state=inflow_state)
        end if
        if (.not.valid) return
        transport_sum=transport_sum+exchange;theta=min(theta,stage_theta)
      end if
      if (half==1) then
        if (recursive_steps) then
          call advance_reactive_eb_tree_subcycle_hydro_3d(species,solver,dt,candidate,boundary,wall,valid, &
            inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
            boundary_faces=boundary_faces,inflow_state=inflow_state, &
            difmag=difmag,extended_geometry=extended_geometry)
        else
          call advance_reactive_eb_tree_hydro_3d(species,solver,dt,candidate,boundary,wall,valid, &
            inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
            boundary_faces=boundary_faces,inflow_state=inflow_state, &
            difmag=difmag,extended_geometry=extended_geometry)
        end if
        if (.not.valid) return
      else if (chemistry_enabled) then
        call reaction_half(valid)
        if (.not.valid) return
      end if
    end do
    boundary=boundary+transport_sum
    if (.not.all(ieee_is_finite(boundary)).or..not.all(ieee_is_finite(chemical))) return
    tree=candidate;boundary_change=boundary;wall_change=wall
    chemical_change=chemical;minimum_theta=theta;ok=.true.
  contains
    subroutine reaction_half(good)
      logical, intent(out) :: good
      call react_reactive_eb_tree_3d(species,reactions,0.5_dp*dt,rtol,atol,candidate,exchange,good,integrator)
      if (good) chemical=chemical+exchange
    end subroutine reaction_half
  end subroutine advance_reactive_eb_tree_full_3d

  ! Common-interval molecular transport: closed outflow edges or open Mach inlets.
  subroutine advance_reactive_eb_tree_transport_3d(species,transport,dt,viscosity,conduction, &
      diffusion,barodiffusion,tree,boundary_change,minimum_theta,ok, &
      inflow_mach,inflow_temperature,inflow_pressure,boundary_faces,inflow_state)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    logical, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), intent(in) :: dt
    type(eb_reactive_tree_3d), intent(inout) :: tree
    real(dp), intent(inout) :: boundary_change(:),minimum_theta
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    logical, intent(out) :: ok
    type(eb_reactive_tree_3d) :: snapshot,candidate
    type(eb_integrated_flux_3d), allocatable :: flux(:)
    real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),boundary(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp) :: theta,stage_theta
    logical :: faces(6)
    integer :: nv,node,n(3),depth,c
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    if (.not.ieee_is_finite(dt).or.dt<=0._dp) return
    if (size(boundary_change)/=nv) return
    if (.not.valid_tree(tree,nv)) return
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if (.not.valid) return
    snapshot=tree
    call average_down_reactive_eb_tree_3d(species,snapshot,valid)
    if (.not.valid) return
    candidate=snapshot
    allocate(flux(size(tree%nodes)),boundary(nv));theta=1._dp
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      allocate(ghost(nv,0:n(1)+1,0:n(2)+1,0:n(3)+1),gt(0:n(1)+1,0:n(2)+1,0:n(3)+1))
      call fill_reactive_eb_tree_halo_3d(species,snapshot,node,ghost,gt,valid, &
        inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state)
      if (.not.valid) return
      depth=tree%topology%nodes(node)%level+1
      ! Inactive physical domain faces remain closed to diffusion. Parent patch edges
      ! inside the domain are donor interfaces, even when no peer exists.
      faces=[origin(1,node)>0,origin(1,node)+n(1)<domain(1,depth), &
        origin(2,node)>0,origin(2,node)+n(2)<domain(2,depth), &
        origin(3,node)>0,origin(3,node)+n(3)<domain(3,depth)]
      if (present(inflow_mach)) faces=faces.or.(inflow_mach>=0)
      if (present(boundary_faces)) faces=.true.
      call advance_reactive_eb_transport_patch_euler_3d(species,transport,snapshot%nodes(node)%geometry, &
        dt,viscosity,conduction,diffusion,barodiffusion,candidate%nodes(node)%state, &
        candidate%nodes(node)%temperature,flux(node),stage_theta,valid,ghost,gt,faces)
      if (.not.valid) return
      theta=min(theta,stage_theta)
      deallocate(ghost,gt)
    end do
    call synchronize_reactive_eb_tree_flux_3d(species,candidate,flux,valid,boundary_faces=boundary_faces)
    if (.not.valid) return
    n=tree%topology%nodes(1)%shape
    associate(g=>tree%nodes(1)%geometry,f=>flux(1))
    do c=1,nv
      boundary(c)=g%dy*g%dz*sum(g%x_face_fraction(0,:,:)*f%x(c,0,:,:)- &
        g%x_face_fraction(n(1),:,:)*f%x(c,n(1),:,:))+ &
        g%dx*g%dz*sum(g%y_face_fraction(:,0,:)*f%y(c,:,0,:)- &
        g%y_face_fraction(:,n(2),:)*f%y(c,:,n(2),:))+ &
        g%dx*g%dy*sum(g%z_face_fraction(:,:,0)*f%z(c,:,:,0)- &
        g%z_face_fraction(:,:,n(3))*f%z(c,:,:,n(3)))
    end do
    end associate
    if (.not.all(ieee_is_finite(boundary))) return
    tree=candidate;boundary_change=boundary;minimum_theta=theta;ok=.true.
  end subroutine advance_reactive_eb_tree_transport_3d

  ! Serial common-interval hydro with PCM donors and outflow/Mach inlet edges.
  ! Supplying extended geometry selects PLM CTU; omission retains PCM Euler.
  ! Every node reads the same averaged snapshot; all outputs commit together.
  subroutine advance_reactive_eb_tree_hydro_3d(species,solver,dt,tree,boundary_change,wall_change,ok, &
      inflow_mach,inflow_temperature,inflow_pressure,difmag,extended_geometry,boundary_faces,inflow_state)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), optional, intent(in) :: extended_geometry(:)
    type(eb_geometry_3d) :: working,center
    type(eb_plm_boundary_3d) :: physical_boundary
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: flat(:,:,:)
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    type(eb_reactive_tree_3d), intent(inout) :: tree
    real(dp), intent(inout) :: boundary_change(:),wall_change(:)
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: difmag
    real(dp) :: coefficient
    real(dp), allocatable :: vf(:,:,:)
    logical, intent(out) :: ok
    type(eb_reactive_tree_3d) :: snapshot,candidate
    type(eb_integrated_flux_3d), allocatable :: flux(:)
    real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),source(:,:,:,:),boundary(:),wall(:)
    logical, allocatable :: exposed(:,:,:)
    integer :: nv,node,child,n(3),lo(3),hi(3),c,level
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    if (.not.ieee_is_finite(dt).or.dt<=0._dp) return
    if (size(boundary_change)/=nv.or.size(wall_change)/=nv) return
    if (.not.valid_tree(tree,nv)) return
    coefficient=.1_dp
    if(present(difmag)) coefficient=difmag
    if(.not.ieee_is_finite(coefficient).or.coefficient<0) return
    if(present(extended_geometry)) then
      if(solver/='pelec'.or.size(extended_geometry)/=size(tree%nodes)) return
      call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
      if(.not.valid) return
      if(any(domain>int(huge(1)-4,int64))) return
      do node=1,size(tree%nodes)
        n=tree%topology%nodes(node)%shape
        if(any([extended_geometry(node)%nx,extended_geometry(node)%ny,extended_geometry(node)%nz]/=n+12)) return
        call extract_eb_geometry_3d(extended_geometry(node),[7,7,7],n+6,center,valid)
        if(.not.valid) return
        if(.not.matching_eb_metrics_3d(center,tree%nodes(node)%geometry)) return
      end do
    end if
    snapshot=tree
    call average_down_reactive_eb_tree_3d(species,snapshot,valid)
    if(.not.valid) return
    candidate=snapshot
    allocate(flux(size(tree%nodes)),boundary(nv),wall(nv));wall=0._dp
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      if(present(extended_geometry)) then
        allocate(ghost(nv,-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),gt(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), &
          vf(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),flat(0:n(1)+7,0:n(2)+7,0:n(3)+7))
        call fill_reactive_eb_tree_halo_3d(species,snapshot,node,ghost,gt,valid, &
          inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
          boundary_faces=boundary_faces,inflow_state=inflow_state, &
          ghost_volume_fraction=vf,halo_width=6,target_volume_fraction=extended_geometry(node)%volume_fraction)
        if(.not.valid) return
        call extract_eb_geometry_3d(extended_geometry(node),[4,4,4],n+9,working,valid)
        if(.not.valid) return
        call compute_reactive_eb_flattening_3d(species,ghost,gt,vf,n+6,flat,valid)
        if(.not.valid) return
        level=tree%topology%nodes(node)%level+1
        physical_boundary%lower=4-int(origin(:,node))
        physical_boundary%upper=3+int(domain(:,level)-origin(:,node))
        physical_boundary%faces='outflow'
        if(present(boundary_faces)) physical_boundary%faces=boundary_faces
        if(present(inflow_mach)) then
          where(inflow_mach>=0) physical_boundary%faces='inflow'
        end if
        call advance_reactive_eb_plm_ctu_3d(species,working,ghost,gt,vf,dt,flat,candidate%nodes(node)%state, &
          candidate%nodes(node)%temperature,valid,difmag=coefficient,integrated_flux=flux(node),boundary=physical_boundary)
        if(.not.valid) return
        do c=1,nv
          where(snapshot%nodes(node)%geometry%volume_fraction<=0) &
            candidate%nodes(node)%state(c,:,:,:)=snapshot%nodes(node)%state(c,:,:,:)
        end do
        where(snapshot%nodes(node)%geometry%volume_fraction<=0) &
          candidate%nodes(node)%temperature=snapshot%nodes(node)%temperature
        deallocate(flat)
      else
      allocate(ghost(nv,0:n(1)+1,0:n(2)+1,0:n(3)+1),gt(0:n(1)+1,0:n(2)+1,0:n(3)+1))
      if(coefficient>0) allocate(vf(0:n(1)+1,0:n(2)+1,0:n(3)+1))
      call fill_reactive_eb_tree_halo_3d(species,snapshot,node,ghost,gt,valid, &
        inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,ghost_volume_fraction=vf)
      if(.not.valid) return
      call advance_reactive_eb_euler_3d(species,snapshot%nodes(node)%state,snapshot%nodes(node)%temperature, &
        snapshot%nodes(node)%geometry,solver,dt,candidate%nodes(node)%state,candidate%nodes(node)%temperature, &
        valid,ghost,gt,flux(node),difmag=coefficient,ghost_volume_fraction=vf)
      end if
      if(.not.valid) return
      deallocate(ghost,gt)
      if(allocated(vf)) deallocate(vf)
      allocate(source(nv,n(1),n(2),n(3)),exposed(n(1),n(2),n(3)));exposed=.true.
      do child=node+1,size(tree%nodes)
        if (tree%topology%nodes(child)%parent/=node) cycle
        call bounds(tree%topology,child,lo,hi)
        exposed(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.false.
      end do
      call reactive_eb_slip_wall_source_3d(species,snapshot%nodes(node)%state,snapshot%nodes(node)%temperature, &
        snapshot%nodes(node)%geometry,source,valid)
      if(.not.valid) return
      associate(g=>snapshot%nodes(node)%geometry)
      do c=1,nv
        wall(c)=wall(c)+dt*g%dx*g%dy*g%dz*sum(source(c,:,:,:)*g%volume_fraction,mask=exposed)
      end do
      end associate
      deallocate(source,exposed)
    end do
    call synchronize_reactive_eb_tree_flux_3d(species,candidate,flux,valid,boundary_faces=boundary_faces)
    if(.not.valid) return
    n=tree%topology%nodes(1)%shape
    associate(g=>tree%nodes(1)%geometry,f=>flux(1))
    do c=1,nv
      boundary(c)=g%dy*g%dz*sum(g%x_face_fraction(0,:,:)*f%x(c,0,:,:)- &
        g%x_face_fraction(n(1),:,:)*f%x(c,n(1),:,:))+ &
        g%dx*g%dz*sum(g%y_face_fraction(:,0,:)*f%y(c,:,0,:)- &
        g%y_face_fraction(:,n(2),:)*f%y(c,:,n(2),:))+ &
        g%dx*g%dy*sum(g%z_face_fraction(:,:,0)*f%z(c,:,:,0)- &
        g%z_face_fraction(:,:,n(3))*f%z(c,:,:,n(3)))
    end do
    end associate
    if (.not.all(ieee_is_finite(boundary)).or..not.all(ieee_is_finite(wall))) return
    tree=candidate;boundary_change=boundary;wall_change=wall;ok=.true.
  end subroutine advance_reactive_eb_tree_hydro_3d

  ! Each level advances together, with ratio fine substeps, time-interpolated
  ! ancestor donors, per-substep peer fluxes and interval-local reflux.
  subroutine advance_reactive_eb_tree_subcycle_hydro_3d(species,solver,dt,tree,boundary_change,wall_change,ok, &
      inflow_mach,inflow_temperature,inflow_pressure,difmag,extended_geometry,boundary_faces,inflow_state)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), optional, intent(in) :: extended_geometry(:)
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    type(eb_reactive_tree_3d), intent(inout) :: tree
    real(dp), intent(inout) :: boundary_change(:),wall_change(:)
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: difmag
    logical, intent(out) :: ok
    call advance_reactive_eb_tree_subcycle_operator_3d(species,solver,dt,tree,boundary_change,wall_change,ok, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
      boundary_faces=boundary_faces,inflow_state=inflow_state, &
            difmag=difmag,extended_geometry=extended_geometry)
  end subroutine advance_reactive_eb_tree_subcycle_hydro_3d

  subroutine advance_reactive_eb_tree_subcycle_transport_3d(species,transport,dt,viscosity,conduction, &
      diffusion,barodiffusion,tree,boundary_change,minimum_theta,ok, &
      inflow_mach,inflow_temperature,inflow_pressure,boundary_faces,inflow_state)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: dt
    logical, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    type(eb_reactive_tree_3d), intent(inout) :: tree
    real(dp), intent(inout) :: boundary_change(:),minimum_theta
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    logical, intent(out) :: ok
    real(dp) :: wall(size(boundary_change))
    call advance_reactive_eb_tree_subcycle_operator_3d(species,'',dt,tree,boundary_change,wall,ok, &
      transport,viscosity,conduction,diffusion,barodiffusion,minimum_theta, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
      boundary_faces=boundary_faces,inflow_state=inflow_state)
  end subroutine advance_reactive_eb_tree_subcycle_transport_3d

  subroutine advance_reactive_eb_tree_subcycle_operator_3d(species,solver,dt,tree,boundary_change,wall_change,ok, &
      transport,viscosity,conduction,diffusion,barodiffusion,minimum_theta, &
      inflow_mach,inflow_temperature,inflow_pressure,difmag,extended_geometry,boundary_faces,inflow_state)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), optional, intent(in) :: extended_geometry(:)
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    type(eb_reactive_tree_3d), intent(inout) :: tree
    real(dp), intent(inout) :: boundary_change(:),wall_change(:)
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: difmag
    logical, intent(out) :: ok
    type(gas_transport_species), optional, intent(in) :: transport(:)
    logical, optional, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), optional, intent(inout) :: minimum_theta
    type(eb_reactive_tree_3d) :: candidate,starts,ends
    type(eb_integrated_flux_3d), allocatable :: total(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: clock_start(:),clock_dt(:),wall(:),boundary(:)
    real(dp) :: theta,coefficient
    type(eb_geometry_3d) :: center
    integer :: nv,levels,c,n(3),node
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    if (.not.ieee_is_finite(dt).or.dt<=0._dp) return
    if (size(boundary_change)/=nv.or.size(wall_change)/=nv) return
    if (.not.valid_tree(tree,nv)) return
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if (.not.valid) return
    if(present(extended_geometry)) then
      if(present(transport).or.solver/='pelec') return
      if(size(extended_geometry)/=size(tree%nodes)) return
      if(any(domain>int(huge(1)-4,int64))) return
      do node=1,size(tree%nodes)
        n=tree%topology%nodes(node)%shape
        if(any([extended_geometry(node)%nx,extended_geometry(node)%ny,extended_geometry(node)%nz]/=n+12)) return
        call extract_eb_geometry_3d(extended_geometry(node),[7,7,7],n+6,center,valid)
        if(.not.valid) return
        if(.not.matching_eb_metrics_3d(center,tree%nodes(node)%geometry)) return
      end do
    end if
    if (present(transport)) then
      if (.not.(present(viscosity).and.present(conduction).and.present(diffusion).and. &
        present(barodiffusion).and.present(minimum_theta))) return
      if (any(ratios>huge(1)/ratios)) return
    end if
    coefficient=.1_dp
    if(present(transport)) coefficient=0
    if(present(difmag)) coefficient=difmag
    if(.not.ieee_is_finite(coefficient).or.coefficient<0) return
    if(present(transport).and.coefficient>0) return
    theta=1._dp
    levels=size(ratios);candidate=tree
    call average_down_reactive_eb_tree_3d(species,candidate,valid)
    if (.not.valid) return
    starts=candidate;ends=candidate
    allocate(clock_start(levels),clock_dt(levels),wall(nv),boundary(nv));wall=0._dp
    call advance_level(0,0._dp,dt,total,valid)
    if (.not.valid) return
    n=tree%topology%nodes(1)%shape
    associate(g=>tree%nodes(1)%geometry,f=>total(1))
    do c=1,nv
      boundary(c)=g%dy*g%dz*sum(g%x_face_fraction(0,:,:)*f%x(c,0,:,:)- &
        g%x_face_fraction(n(1),:,:)*f%x(c,n(1),:,:))+ &
        g%dx*g%dz*sum(g%y_face_fraction(:,0,:)*f%y(c,:,0,:)- &
        g%y_face_fraction(:,n(2),:)*f%y(c,:,n(2),:))+ &
        g%dx*g%dy*sum(g%z_face_fraction(:,:,0)*f%z(c,:,:,0)- &
        g%z_face_fraction(:,:,n(3))*f%z(c,:,:,n(3)))
    end do
    end associate
    if (.not.all(ieee_is_finite(boundary)).or..not.all(ieee_is_finite(wall))) return
    tree=candidate;boundary_change=boundary;wall_change=wall
    if (present(minimum_theta)) minimum_theta=theta
    ok=.true.
  contains
    recursive subroutine advance_level(level,time,interval,flux,good)
      integer, intent(in) :: level
      real(dp), intent(in) :: time,interval
      type(eb_integrated_flux_3d), allocatable, intent(out) :: flux(:)
      logical, intent(out) :: good
      type(eb_integrated_flux_3d), allocatable :: fine_flux(:)
      type(eb_integrated_flux_3d) :: owned
      real(dp) :: period(3),length(3)
      integer :: seam,tangent
      logical :: overlap
      type(eb_reactive_tree_3d) :: snapshot
      type(eb_geometry_3d) :: working
      type(eb_plm_boundary_3d) :: physical_boundary
      real(dp), allocatable :: flat(:,:,:)
      real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),source(:,:,:,:),fraction(:),vf(:,:,:)
      logical, allocatable :: exposed(:,:,:)
      real(dp) :: elapsed,sub_dt,stage_theta
      integer :: node,peer,owner,receiver,axis,child,substep,d(3),lo(3),hi(3),component,substeps
      integer(int64) :: offset(3)
      logical :: valid,faces(6)
      good=.false.;snapshot=candidate
      clock_start(level+1)=time;clock_dt(level+1)=interval
      allocate(flux(size(tree%nodes)),fraction(level))
      do node=1,size(tree%nodes)
        d=tree%topology%nodes(node)%shape
        allocate(flux(node)%x(nv,0:d(1),d(2),d(3)),flux(node)%y(nv,d(1),0:d(2),d(3)), &
          flux(node)%z(nv,d(1),d(2),0:d(3)))
        flux(node)%x=0._dp;flux(node)%y=0._dp;flux(node)%z=0._dp
        if (tree%topology%nodes(node)%level/=level) cycle
        starts%nodes(node)=snapshot%nodes(node)
        fraction=max(0._dp,min(1._dp,(time-clock_start(:level))/clock_dt(:level)))
        if(present(extended_geometry)) then
          allocate(ghost(nv,-5:d(1)+6,-5:d(2)+6,-5:d(3)+6),gt(-5:d(1)+6,-5:d(2)+6,-5:d(3)+6), &
            vf(-5:d(1)+6,-5:d(2)+6,-5:d(3)+6),flat(0:d(1)+7,0:d(2)+7,0:d(3)+7))
          call fill_reactive_eb_tree_halo_3d(species,snapshot,node,ghost,gt,valid,starts,ends,fraction, &
            inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
            boundary_faces=boundary_faces,inflow_state=inflow_state, &
            ghost_volume_fraction=vf,halo_width=6,target_volume_fraction=extended_geometry(node)%volume_fraction)
          if(.not.valid) return
          call extract_eb_geometry_3d(extended_geometry(node),[4,4,4],d+9,working,valid)
          if(.not.valid) return
          call compute_reactive_eb_flattening_3d(species,ghost,gt,vf,d+6,flat,valid)
          if(.not.valid) return
          physical_boundary%lower=4-int(origin(:,node))
          physical_boundary%upper=3+int(domain(:,level+1)-origin(:,node))
          physical_boundary%faces='outflow'
          if(present(boundary_faces)) physical_boundary%faces=boundary_faces
          if(present(inflow_mach)) then
            where(inflow_mach>=0) physical_boundary%faces='inflow'
          end if
        else
        allocate(ghost(nv,0:d(1)+1,0:d(2)+1,0:d(3)+1),gt(0:d(1)+1,0:d(2)+1,0:d(3)+1))
        if(coefficient>0) allocate(vf(0:d(1)+1,0:d(2)+1,0:d(3)+1))
        call fill_reactive_eb_tree_halo_3d(species,snapshot,node,ghost,gt,valid,starts,ends,fraction, &
          inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
          boundary_faces=boundary_faces,inflow_state=inflow_state,ghost_volume_fraction=vf)
        if (.not.valid) return
        end if
        if (present(transport)) then
          faces=[origin(1,node)>0,origin(1,node)+int(d(1),int64)<domain(1,level+1), &
            origin(2,node)>0,origin(2,node)+int(d(2),int64)<domain(2,level+1), &
            origin(3,node)>0,origin(3,node)+int(d(3),int64)<domain(3,level+1)]
          if (present(inflow_mach)) faces=faces.or.(inflow_mach>=0)
          if (present(boundary_faces)) faces=.true.
          call advance_reactive_eb_transport_patch_euler_3d(species,transport,snapshot%nodes(node)%geometry, &
            interval,viscosity,conduction,diffusion,barodiffusion,candidate%nodes(node)%state, &
            candidate%nodes(node)%temperature,flux(node),stage_theta,valid,ghost,gt,faces)
          if (.not.valid) return
          theta=min(theta,stage_theta)
        else if(present(extended_geometry)) then
          call advance_reactive_eb_plm_ctu_3d(species,working,ghost,gt,vf,interval,flat,candidate%nodes(node)%state, &
            candidate%nodes(node)%temperature,valid,difmag=coefficient,integrated_flux=flux(node),boundary=physical_boundary)
          if(.not.valid) return
          do component=1,nv
            where(snapshot%nodes(node)%geometry%volume_fraction<=0) &
              candidate%nodes(node)%state(component,:,:,:)=snapshot%nodes(node)%state(component,:,:,:)
          end do
          where(snapshot%nodes(node)%geometry%volume_fraction<=0) &
            candidate%nodes(node)%temperature=snapshot%nodes(node)%temperature
          deallocate(flat)
        else
        call advance_reactive_eb_euler_3d(species,snapshot%nodes(node)%state,snapshot%nodes(node)%temperature, &
          snapshot%nodes(node)%geometry,solver,interval,candidate%nodes(node)%state,candidate%nodes(node)%temperature, &
          valid,ghost,gt,flux(node),difmag=coefficient,ghost_volume_fraction=vf)
        end if
        if (.not.valid) return
        deallocate(ghost,gt)
        if(allocated(vf)) deallocate(vf)
        if (present(transport)) cycle
        allocate(source(nv,d(1),d(2),d(3)),exposed(d(1),d(2),d(3)));exposed=.true.
        do child=node+1,size(tree%nodes)
          if (tree%topology%nodes(child)%parent/=node) cycle
          call bounds(tree%topology,child,lo,hi)
          exposed(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.false.
        end do
        call reactive_eb_slip_wall_source_3d(species,snapshot%nodes(node)%state,snapshot%nodes(node)%temperature, &
          snapshot%nodes(node)%geometry,source,valid)
        if (.not.valid) return
        associate(g=>snapshot%nodes(node)%geometry)
        do component=1,nv
          wall(component)=wall(component)+interval*g%dx*g%dy*g%dz* &
            sum(source(component,:,:,:)*g%volume_fraction,mask=exposed)
        end do
        end associate
        deallocate(source,exposed)
      end do
      do node=1,size(tree%nodes)
      do peer=node+1,size(tree%nodes)
        if (tree%topology%nodes(node)%level/=level.or.tree%topology%nodes(peer)%level/=level) cycle
        owner=node;receiver=peer
        do axis=1,3
          if (origin(axis,node)==origin(axis,peer)) cycle
          if (origin(axis,node)>origin(axis,peer)) then
            owner=peer;receiver=node
          end if
          exit
        end do
        offset=origin(:,receiver)-origin(:,owner)
        if (any(offset < -int(tree%topology%nodes(receiver)%shape,int64)).or. &
          any(offset>int(tree%topology%nodes(owner)%shape,int64))) cycle
        call reconcile_reactive_eb_peer_flux_3d(species,tree%nodes(owner)%geometry,tree%nodes(receiver)%geometry, &
          int(offset),flux(owner),flux(receiver),candidate%nodes(receiver)%state, &
          candidate%nodes(receiver)%temperature,valid)
        if (.not.valid) return
      end do
      end do
      ! Reconcile periodic provisional endpoints before descendants interpolate
      ! in time; final fine/coarse synchronization follows the recursive steps.
      if(present(boundary_faces)) then
        associate(g=>tree%nodes(1)%geometry)
        length=[g%x_upper-g%x_lower,g%y_upper-g%y_lower,g%z_upper-g%z_lower]
        end associate
        do seam=1,3
          if(boundary_faces(2*seam-1)/='periodic') cycle
          period=0;period(seam)=length(seam)
          do owner=1,size(tree%nodes)
            if(tree%topology%nodes(owner)%level/=level) cycle
            if(origin(seam,owner)+tree%topology%nodes(owner)%shape(seam)/=domain(seam,level+1)) cycle
            owned=flux(owner)
            do receiver=1,size(tree%nodes)
              if(tree%topology%nodes(receiver)%level/=level.or.origin(seam,receiver)/=0) cycle
              offset=origin(:,receiver)-origin(:,owner)
              offset(seam)=offset(seam)+domain(seam,level+1)
              overlap=.true.
              do tangent=1,3
                if(tangent==seam) cycle
                if(offset(tangent)>=tree%topology%nodes(owner)%shape(tangent).or. &
                  offset(tangent)+tree%topology%nodes(receiver)%shape(tangent)<=0) overlap=.false.
              end do
              if(.not.overlap) cycle
              call reconcile_reactive_eb_peer_flux_3d(species,tree%nodes(owner)%geometry, &
                tree%nodes(receiver)%geometry,int(offset),owned,flux(receiver),candidate%nodes(receiver)%state, &
                candidate%nodes(receiver)%temperature,valid,periodic_length=period)
              if(.not.valid) return
            end do
          end do
        end do
      end if
      do node=1,size(tree%nodes)
        if (tree%topology%nodes(node)%level==level) ends%nodes(node)=candidate%nodes(node)
      end do
      if (level+1<levels) then
        elapsed=0._dp
        substeps=ratios(level+2)
        if (present(transport)) substeps=substeps**2
        do substep=1,substeps
          sub_dt=interval/real(substeps,dp)
          if (substep==substeps) sub_dt=interval-elapsed
          if (sub_dt<=0._dp.or.time+elapsed+sub_dt<=time+elapsed) return
          call advance_level(level+1,time+elapsed,sub_dt,fine_flux,valid)
          if (.not.valid) return
          do node=1,size(tree%nodes)
            if (tree%topology%nodes(node)%level/=level+1) cycle
            flux(node)%x=flux(node)%x+fine_flux(node)%x
            flux(node)%y=flux(node)%y+fine_flux(node)%y
            flux(node)%z=flux(node)%z+fine_flux(node)%z
          end do
          elapsed=elapsed+sub_dt
        end do
      end if
      call synchronize_reactive_eb_tree_flux_3d(species,candidate,flux,valid,level,min(level+1,levels-1),boundary_faces)
      if (.not.valid) return
      good=.true.
    end subroutine advance_level
  end subroutine advance_reactive_eb_tree_subcycle_operator_3d

  ! Input states have already accepted each node's raw time-integrated flux.
  ! Replace shared transfers from finest to coarsest on exposed cells only,
  ! then average down. Publish corrected states and effective fluxes together.
  subroutine synchronize_reactive_eb_tree_flux_3d(species,tree,flux,ok,minimum_level,maximum_level,boundary_faces)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_reactive_tree_3d), intent(inout) :: tree
    type(eb_integrated_flux_3d), intent(inout) :: flux(:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: minimum_level,maximum_level
    character(len=*), optional, intent(in) :: boundary_faces(6)
    character(len=16) :: faces(6)
    real(dp) :: period(3),length(3)
    type(eb_reactive_tree_3d) :: candidate
    type(eb_integrated_flux_3d), allocatable :: effective(:)
    type(eb_integrated_flux_3d) :: restricted,owned
    type(eb_geometry_3d) :: section
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer(int64) :: offset(3),section_origin(3)
    integer :: nv,node,peer,owner,receiver,axis,depth,parent,child,lo(3),hi(3),d(3),ratio,first_level,last_level,seam,side
    logical, allocatable :: refined(:,:,:)
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    if (.not.valid_tree(tree,nv)) return
    faces='outflow'
    if(present(boundary_faces)) faces=boundary_faces
    if(any(faces/='outflow'.and.faces/='inflow'.and.faces/='periodic'.and. &
      faces/='slip_wall'.and.faces/='no_slip_wall')) return
    if(.not.periodic_eb_faces_match(tree%nodes(1)%geometry,faces)) return
    associate(g=>tree%nodes(1)%geometry)
    length=[g%x_upper-g%x_lower,g%y_upper-g%y_lower,g%z_upper-g%z_lower]
    end associate
    first_level=0;last_level=tree%topology%level_count()-1
    if (present(minimum_level)) first_level=minimum_level
    if (present(maximum_level)) last_level=maximum_level
    if (first_level<0.or.last_level<first_level.or.last_level>=tree%topology%level_count()) return
    if (size(flux)/=size(tree%nodes)) return
    do node=1,size(flux)
      if (.not.valid_eb_integrated_flux_3d(flux(node),tree%nodes(node)%geometry,nv)) return
    end do
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if (.not.valid) return
    candidate=tree;effective=flux
    ! Reconcile raw peer faces first at every depth. Descendant replacements
    ! below are mirrored into both incident parent patches and must not be
    ! overwritten later by a coarser owner decision.
    do node=1,size(flux)
    do peer=node+1,size(flux)
      if (tree%topology%nodes(node)%level<first_level.or.tree%topology%nodes(node)%level>last_level) cycle
      if (tree%topology%nodes(peer)%level/=tree%topology%nodes(node)%level) cycle
      owner=node;receiver=peer
      do axis=1,3
        if (origin(axis,node)==origin(axis,peer)) cycle
        if (origin(axis,node)>origin(axis,peer)) then
          owner=peer;receiver=node
        end if
        exit
      end do
      offset=origin(:,receiver)-origin(:,owner)
      if (any(offset < -int(tree%topology%nodes(receiver)%shape,int64)).or. &
        any(offset>int(tree%topology%nodes(owner)%shape,int64))) cycle
      call reconcile_reactive_eb_peer_flux_3d(species,tree%nodes(owner)%geometry,tree%nodes(receiver)%geometry, &
        int(offset),effective(owner),effective(receiver),ok=valid)
      if (.not.valid) return
    end do
    end do
    ! Canonical raw periodic transfer comes from the high-domain face.
    ! Copy owner flux so a patch spanning the domain can receive its own seam.
    do seam=1,3
      if(faces(2*seam-1)/='periodic') cycle
      period=0;period(seam)=length(seam)
      do owner=1,size(flux)
        depth=tree%topology%nodes(owner)%level
        if(depth<first_level.or.depth>last_level) cycle
        if(origin(seam,owner)+tree%topology%nodes(owner)%shape(seam)/=domain(seam,depth+1)) cycle
        owned=effective(owner)
        do receiver=1,size(flux)
          if(tree%topology%nodes(receiver)%level/=depth.or.origin(seam,receiver)/=0) cycle
          offset=origin(:,receiver)-origin(:,owner)
          offset(seam)=offset(seam)+domain(seam,depth+1)
          if(.not.seam_overlap(offset,tree%topology%nodes(owner)%shape, &
            tree%topology%nodes(receiver)%shape,seam)) cycle
          call reconcile_reactive_eb_peer_flux_3d(species,tree%nodes(owner)%geometry,tree%nodes(receiver)%geometry, &
            int(offset),owned,effective(receiver),ok=valid,periodic_length=period)
          if(.not.valid) return
        end do
      end do
    end do
    do depth=last_level,first_level,-1
    do node=1,size(flux)
      if (tree%topology%nodes(node)%level/=depth) cycle
      d=tree%topology%nodes(node)%shape
      allocate(refined(d(1),d(2),d(3)));refined=.false.
      do child=node+1,size(flux)
        if (tree%topology%nodes(child)%parent/=node) cycle
        call bounds(tree%topology,child,lo,hi)
        refined(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.true.
      end do
      call reflux_reactive_eb_state_3d(species,tree%nodes(node)%geometry,refined,flux(node),effective(node), &
        candidate%nodes(node)%state,candidate%nodes(node)%temperature,valid,all_faces=.true.)
      if (.not.valid) return
      deallocate(refined)
      if (depth==first_level) cycle
      parent=tree%topology%nodes(node)%parent;ratio=tree%topology%nodes(node)%patch%refinement_ratio
      call bounds(tree%topology,node,lo,hi)
      call restrict_eb_integrated_flux_3d(tree%nodes(node)%geometry,ratio,effective(node),restricted,valid)
      if (.not.valid) return
      effective(parent)%x(:,lo(1)-1:hi(1),lo(2):hi(2),lo(3):hi(3))=restricted%x
      effective(parent)%y(:,lo(1):hi(1),lo(2)-1:hi(2),lo(3):hi(3))=restricted%y
      effective(parent)%z(:,lo(1):hi(1),lo(2):hi(2),lo(3)-1:hi(3))=restricted%z
      call coarsen_eb_geometry_3d(tree%nodes(node)%geometry,ratio,section,valid)
      if (.not.valid) return
      section_origin=origin(:,node)/int(ratio,int64)
      do peer=1,size(flux)
        if (peer==parent.or.tree%topology%nodes(peer)%level/=depth-1) cycle
        offset=origin(:,peer)-section_origin
        if (any(offset < -int(tree%topology%nodes(peer)%shape,int64)).or.any(offset>int(hi-lo+1,int64))) cycle
        call reconcile_reactive_eb_peer_flux_3d(species,section,tree%nodes(peer)%geometry,int(offset), &
          restricted,effective(peer),ok=valid)
        if (.not.valid) return
      end do
      ! A restricted fine face also replaces the opposite periodic parent
      ! face, including a root that is both parent and opposite receiver.
      do seam=1,3
        if(faces(2*seam-1)/='periodic') cycle
        period=0;period(seam)=length(seam)
        do side=0,1
          if(side==0.and.section_origin(seam)/=0) cycle
          if(side==1.and.section_origin(seam)+hi(seam)-lo(seam)+1/=domain(seam,depth)) cycle
          do peer=1,size(flux)
            if(tree%topology%nodes(peer)%level/=depth-1) cycle
            if(side==0.and.origin(seam,peer)+tree%topology%nodes(peer)%shape(seam)/=domain(seam,depth)) cycle
            if(side==1.and.origin(seam,peer)/=0) cycle
            offset=origin(:,peer)-section_origin
            offset(seam)=offset(seam)+int(2*side-1,int64)*domain(seam,depth)
            if(.not.seam_overlap(offset,hi-lo+1,tree%topology%nodes(peer)%shape,seam)) cycle
            call reconcile_reactive_eb_peer_flux_3d(species,section,tree%nodes(peer)%geometry,int(offset), &
              restricted,effective(peer),ok=valid,periodic_length=period)
            if(.not.valid) return
          end do
        end do
      end do
    end do
    end do
    call average_down_reactive_eb_tree_3d(species,candidate,valid)
    if (.not.valid) return
    tree=candidate;flux=effective;ok=.true.
  contains
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
  end subroutine synchronize_reactive_eb_tree_flux_3d

  ! Six-face Mach < 0 retains outflow; active entries require inlet temperature.
  ! Pressure zero extrapolates the interior value; positive entries prescribe it.
  ! Requesting volume fractions fills interior/edges/corners too. Ancestor
  ! fractions use PCM donor values; physical corners use the first active inlet
  ! axis (x, then y, then z). No geometric fine-cell reconstruction is implied.
  ! Wider halos require explicit target-resolution volume fractions, including
  ! edges/corners. Actual output bounds correspond to 1-width:N+width; array
  ! assignment preserves caller bounds. Width one retains the legacy face route.
  subroutine fill_reactive_eb_tree_halo_3d(species,tree,target,ghost,ghost_temperature,ok, &
      ancestor_start,ancestor_end,ancestor_fraction,inflow_mach,inflow_temperature,inflow_pressure, &
      ghost_volume_fraction,halo_width,target_volume_fraction,boundary_faces,inflow_state)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_reactive_tree_3d), intent(in) :: tree
    integer, intent(in) :: target
    integer, optional, intent(in) :: halo_width
    real(dp), optional, intent(in) :: target_volume_fraction(:,:,:)
    real(dp), intent(inout) :: ghost(:,0:,0:,0:),ghost_temperature(0:,0:,0:)
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    real(dp), optional, intent(inout) :: ghost_volume_fraction(0:,0:,0:)
    real(dp), allocatable :: volumes(:,:,:)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    character(len=16) :: faces(6)
    integer :: inlet_face,momentum_sign(3)
    logical, allocatable :: filled(:,:,:)
    logical, intent(out) :: ok
    type(eb_reactive_tree_3d), optional, intent(in) :: ancestor_start,ancestor_end
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    type(eb_reactive_tree_3d) :: snapshot
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer(int64) :: global(3),mapped(3),local(3),offset(3)
    real(dp), allocatable :: candidate(:,:,:,:),temps(:,:,:),primitive(:)
    real(dp) :: area,guess,recovered,speed,fraction,value(reactive_nvar(size(species)))
    integer :: n(3),nv,level,depth,node,donor,index(3),inside(3),source(3),first(3),last(3),axis,side,i,j,k,width
    real(dp) :: mach(6),pin(6),tin(6),inlet_value(reactive_nvar(size(species))),inlet_t
    logical :: valid,physical
    ok=.false.;nv=reactive_nvar(size(species));width=1
    if(present(halo_width)) width=halo_width
    if(width<1.or.width>huge(1)/4) return
    if(width>1.and.(.not.present(ghost_volume_fraction).or..not.present(target_volume_fraction))) return
    faces='outflow'
    if(present(boundary_faces)) then
      if(any(boundary_faces/='inflow'.and.boundary_faces/='outflow'.and.boundary_faces/='periodic'.and. &
        boundary_faces/='slip_wall'.and.boundary_faces/='no_slip_wall')) return
      if(present(inflow_mach).or.present(inflow_pressure)) return
      faces=boundary_faces
      if(any(faces=='inflow')) then
        if(.not.present(inflow_state).or..not.present(inflow_temperature)) return
      end if
    else
      if(present(inflow_state)) return
    end if
    if(present(inflow_state)) then
      if(any(shape(inflow_state)/=[nv,6])) return
      if(.not.all(ieee_is_finite(inflow_state))) return
    end if
    mach=-1;pin=0;tin=0
    if (present(inflow_temperature).or.present(inflow_pressure)) then
      if (.not.present(inflow_mach).and..not.present(boundary_faces)) return
    end if
    if (present(inflow_mach)) mach=inflow_mach
    if (present(inflow_temperature)) tin=inflow_temperature
    if (present(inflow_pressure)) pin=inflow_pressure
    if (.not.all(ieee_is_finite(mach)).or.any(mach>100)) return
    if (.not.all(ieee_is_finite(pin)).or.any(pin<0).or.any(pin>1e15_dp)) return
    if (.not.all(ieee_is_finite(tin))) return
    if (any(pin>0.and.mach<0).or.any(tin<=0.and.mach>=0)) return
    if (.not.valid_tree(tree,nv)) return
    if(.not.periodic_eb_faces_match(tree%nodes(1)%geometry,faces)) return
    if (target<1.or.target>size(tree%nodes)) return
    n=tree%topology%nodes(target)%shape;level=tree%topology%nodes(target)%level+1
    if(any(n>huge(1)-2*width)) return
    if (any(shape(ghost)/=[nv,n+2*width]).or.any(shape(ghost_temperature)/=n+2*width)) return
    if(present(ghost_volume_fraction)) then
      if(any(shape(ghost_volume_fraction)/=n+2*width)) return
    end if
    if(present(target_volume_fraction)) then
      if(.not.present(ghost_volume_fraction)) return
      if(any(shape(target_volume_fraction)/=n+2*width)) return
      if(.not.all(ieee_is_finite(target_volume_fraction))) return
      if(any(target_volume_fraction<0).or.any(target_volume_fraction>1)) return
      if(any(target_volume_fraction(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))/= &
        tree%nodes(target)%geometry%volume_fraction)) return
    end if
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if (.not.valid) return
    if (present(ancestor_start).or.present(ancestor_end).or.present(ancestor_fraction)) then
      if (.not.(present(ancestor_start).and.present(ancestor_end).and.present(ancestor_fraction))) return
      if (size(ancestor_fraction)/=level-1) return
      if (.not.all(ieee_is_finite(ancestor_fraction))) return
      if (any(ancestor_fraction<0._dp).or.any(ancestor_fraction>1._dp)) return
      if (.not.matching_snapshot(ancestor_start).or..not.matching_snapshot(ancestor_end)) return
    end if
    snapshot=tree
    call average_down_reactive_eb_tree_3d(species,snapshot,valid)
    if (.not.valid) return
    allocate(candidate(nv,1-width:n(1)+width,1-width:n(2)+width,1-width:n(3)+width), &
      temps(1-width:n(1)+width,1-width:n(2)+width,1-width:n(3)+width), &
      primitive(reactive_nprim(size(species))))
    candidate=0._dp;temps=0._dp
    allocate(filled(1-width:n(1)+width,1-width:n(2)+width,1-width:n(3)+width));filled=.false.
    if(width==1) then
    do axis=1,3
    do side=0,1
      first=1;last=n;first(axis)=side*(n(axis)+1);last(axis)=first(axis)
      do k=first(3),last(3)
      do j=first(2),last(2)
      do i=first(1),last(1)
        index=[i,j,k];inside=max(1,min(n,index))
        select case(axis)
        case(1);area=tree%nodes(target)%geometry%x_face_fraction(side*n(1),j,k)
        case(2);area=tree%nodes(target)%geometry%y_face_fraction(i,side*n(2),k)
        case(3);area=tree%nodes(target)%geometry%z_face_fraction(i,j,side*n(3))
        end select
        if (area<=0._dp) cycle
        global=origin(:,target)+int(index-1,int64)
        donor=0
        if(faces(2*axis-1)=='periodic') global(axis)=modulo(global(axis),domain(axis,level))
        physical=any(global<0_int64).or.any(global>=domain(:,level))
        if (physical) then
          ! Physical boundary starts from the target's adjacent interior cell.
          donor=target;source=inside
        else
          do depth=level,1,-1
            mapped=global/(domain(:,level)/domain(:,depth))
            do node=1,size(tree%nodes)
              if (tree%topology%nodes(node)%level+1/=depth) cycle
              local=mapped-origin(:,node)
              if (any(local<0_int64).or.any(local>=int(tree%topology%nodes(node)%shape,int64))) cycle
              donor=node;source=int(local)+1;exit
            end do
            if (donor/=0) exit
          end do
        end if
        if (donor==0) return
        if (snapshot%nodes(donor)%geometry%volume_fraction(source(1),source(2),source(3))<=0._dp) return
        guess=snapshot%nodes(donor)%temperature(source(1),source(2),source(3))
        value=snapshot%nodes(donor)%state(:,source(1),source(2),source(3))
        depth=tree%topology%nodes(donor)%level+1
        if (present(ancestor_fraction)) then
          if (depth<level) then
            fraction=ancestor_fraction(depth)
            ! Endpoint states are provisional level updates. Do not average
            ! them down from descendants that may be at a different time.
            if (fraction==0._dp) then
              value=ancestor_start%nodes(donor)%state(:,source(1),source(2),source(3))
              guess=ancestor_start%nodes(donor)%temperature(source(1),source(2),source(3))
            else if (fraction==1._dp) then
              value=ancestor_end%nodes(donor)%state(:,source(1),source(2),source(3))
              guess=ancestor_end%nodes(donor)%temperature(source(1),source(2),source(3))
            else
              value=(1._dp-fraction)*ancestor_start%nodes(donor)%state(:,source(1),source(2),source(3))+ &
                fraction*ancestor_end%nodes(donor)%state(:,source(1),source(2),source(3))
              guess=(1._dp-fraction)*ancestor_start%nodes(donor)%temperature(source(1),source(2),source(3))+ &
                fraction*ancestor_end%nodes(donor)%temperature(source(1),source(2),source(3))
            end if
          end if
        end if
        if (.not.ieee_is_finite(guess)) return
        if (guess<=0._dp) return
        call reactive_conserved_to_primitive(species,value, &
          guess,primitive,recovered,speed,valid)
        if (.not.valid) return
        if (physical.and.mach(2*axis-1+side)>=0) then
          if (pin(2*axis-1+side)>0) then
            call mach_inlet_ghost_3d(species,value,recovered,2*axis-1+side,tin(2*axis-1+side), &
              mach(2*axis-1+side),inlet_value,inlet_t,valid,inlet_pressure=pin(2*axis-1+side))
          else
            call mach_inlet_ghost_3d(species,value,recovered,2*axis-1+side,tin(2*axis-1+side), &
              mach(2*axis-1+side),inlet_value,inlet_t,valid)
          end if
          if (.not.valid) return
          value=inlet_value;recovered=inlet_t
        end if
        if(physical) then
          call prescribed_inlet(2*axis-1+side,valid)
          if(.not.valid) return
          if(faces(2*axis-1+side)=='slip_wall') value(1+axis)=-value(1+axis)
          if(faces(2*axis-1+side)=='no_slip_wall') value(2:4)=-value(2:4)
        end if
        candidate(:,i,j,k)=value
        temps(i,j,k)=recovered
        filled(i,j,k)=.true.
      end do
      end do
      end do
    end do
    end do
    end if
    if(present(ghost_volume_fraction)) then
      allocate(volumes(1-width:n(1)+width,1-width:n(2)+width,1-width:n(3)+width));volumes=0
      do k=1-width,n(3)+width
      do j=1-width,n(2)+width
      do i=1-width,n(1)+width
        index=[i,j,k];global=origin(:,target)+int(index-1,int64)
        inlet_face=0
        do axis=1,3
          side=-1
          if(global(axis)<0) side=0
          if(global(axis)>=domain(axis,level)) side=1
          if(side<0) cycle
          if((mach(2*axis-1+side)>=0.or.faces(2*axis-1+side)=='inflow').and.inlet_face==0) inlet_face=2*axis-1+side
        end do
        call map_wall_donor(global,momentum_sign)
        global=max(0_int64,min(domain(:,level)-1,global))
        donor=0
        if(all(index>=1).and.all(index<=n)) then
          donor=target;source=index
        else
          do depth=level,1,-1
            mapped=global/(domain(:,level)/domain(:,depth))
            do node=1,size(tree%nodes)
              if(tree%topology%nodes(node)%level+1/=depth) cycle
              local=mapped-origin(:,node)
              if(any(local<0).or.any(local>=int(tree%topology%nodes(node)%shape,int64))) cycle
              donor=node;source=int(local)+1;exit
            end do
            if(donor/=0) exit
          end do
        end if
        if(donor==0) return
        volumes(i,j,k)=snapshot%nodes(donor)%geometry%volume_fraction(source(1),source(2),source(3))
        if(present(target_volume_fraction)) then
          volumes(i,j,k)=target_volume_fraction(i+width,j+width,k+width)
        else if(allocated(snapshot%nodes(target)%geometry%halo_volume_fraction)) then
          local=global-origin(:,target)+1
          do axis=1,3
            if(faces(2*axis-1)=='periodic'.or.faces(2*axis-1)=='slip_wall'.or. &
              faces(2*axis-1)=='no_slip_wall'.or.faces(2*axis)=='slip_wall'.or. &
              faces(2*axis)=='no_slip_wall') local(axis)=index(axis)
          end do
          if(snapshot%nodes(target)%geometry%halo_volume_fraction(local(1),local(2),local(3))>=0) &
            volumes(i,j,k)=snapshot%nodes(target)%geometry%halo_volume_fraction(local(1),local(2),local(3))
        end if
        if(volumes(i,j,k)<=0) cycle
        if(width>1.and.snapshot%nodes(donor)%geometry%volume_fraction(source(1),source(2),source(3))<=0) return
        if(filled(i,j,k)) cycle
        value=snapshot%nodes(donor)%state(:,source(1),source(2),source(3))
        guess=snapshot%nodes(donor)%temperature(source(1),source(2),source(3))
        depth=tree%topology%nodes(donor)%level+1
        if(present(ancestor_fraction)) then
          if(depth<level) then
            fraction=ancestor_fraction(depth)
            value=(1-fraction)*ancestor_start%nodes(donor)%state(:,source(1),source(2),source(3))+ &
              fraction*ancestor_end%nodes(donor)%state(:,source(1),source(2),source(3))
            guess=(1-fraction)*ancestor_start%nodes(donor)%temperature(source(1),source(2),source(3))+ &
              fraction*ancestor_end%nodes(donor)%temperature(source(1),source(2),source(3))
          end if
        end if
        call reactive_conserved_to_primitive(species,value,guess,primitive,recovered,speed,valid)
        if(.not.valid) return
        if(inlet_face>0) then
          if(faces(inlet_face)=='inflow') then
            call prescribed_inlet(inlet_face,valid)
            if(.not.valid) return
          else
          if(pin(inlet_face)>0) then
            call mach_inlet_ghost_3d(species,value,recovered,inlet_face,tin(inlet_face),mach(inlet_face), &
              inlet_value,inlet_t,valid,inlet_pressure=pin(inlet_face))
          else
            call mach_inlet_ghost_3d(species,value,recovered,inlet_face,tin(inlet_face),mach(inlet_face), &
              inlet_value,inlet_t,valid)
          end if
          if(.not.valid) return
          value=inlet_value;recovered=inlet_t
          end if
        end if
        value(2:4)=value(2:4)*real(momentum_sign,dp)
        candidate(:,i,j,k)=value;temps(i,j,k)=recovered
      end do
      end do
      end do
    end if
    ! The peer overlay also verifies open-face aperture and centroid agreement,
    ! including peers that belong to different parent patches.
    do node=1,size(tree%nodes)
      if (node==target.or.tree%topology%nodes(node)%level+1/=level) cycle
      offset=origin(:,node)-origin(:,target)
      if (any(offset < -int(tree%topology%nodes(node)%shape,int64)).or.any(offset>int(n,int64))) cycle
      call overlay_reactive_eb_peer_halo_3d(species,snapshot%nodes(target)%geometry,snapshot%nodes(node)%geometry, &
        int(offset),snapshot%nodes(node)%state,snapshot%nodes(node)%temperature, &
        candidate(:,0:n(1)+1,0:n(2)+1,0:n(3)+1),temps(0:n(1)+1,0:n(2)+1,0:n(3)+1),valid)
      if (.not.valid) return
    end do
    ghost=candidate;ghost_temperature=temps
    if(present(ghost_volume_fraction)) ghost_volume_fraction=volumes
    ok=.true.
  contains
    subroutine map_wall_donor(position,signs)
      integer(int64), intent(inout) :: position(3)
      integer, intent(out) :: signs(3)
      integer :: direction,face
      signs=1
      do direction=1,3
        if(faces(2*direction-1)=='periodic') then
          position(direction)=modulo(position(direction),domain(direction,level))
          cycle
        end if
        do while(position(direction)<0.or.position(direction)>=domain(direction,level))
          face=2*direction
          if(position(direction)<0) face=face-1
          if(faces(face)/='slip_wall'.and.faces(face)/='no_slip_wall') exit
          if(faces(face)=='no_slip_wall') then
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
    subroutine prescribed_inlet(face,accepted)
      integer, intent(in) :: face
      logical, intent(out) :: accepted
      accepted=.true.
      if(faces(face)/='inflow') return
      call reactive_conserved_to_primitive(species,inflow_state(:,face),tin(face), &
        primitive,recovered,speed,accepted)
      if(accepted) value=inflow_state(:,face)
    end subroutine
    logical function matching_snapshot(other) result(matches)
      type(eb_reactive_tree_3d), intent(in) :: other
      integer(int64), allocatable :: other_origin(:,:),other_domain(:,:)
      integer, allocatable :: other_ratios(:)
      integer :: p
      logical :: good
      matches=.false.
      if (.not.valid_tree(other,nv)) return
      if (size(other%nodes)/=size(tree%nodes)) return
      call amr_tree_integer_layout_3d(other%topology,other_origin,other_domain,other_ratios,good)
      if (.not.good) return
      if (size(other_ratios)/=size(ratios)) return
      if (any(other_ratios/=ratios).or.any(other_domain/=domain).or.any(other_origin/=origin)) return
      do p=1,size(tree%nodes)
        if (other%topology%nodes(p)%parent/=tree%topology%nodes(p)%parent.or. &
          other%topology%nodes(p)%level/=tree%topology%nodes(p)%level) return
        if (.not.matching_eb_metrics_3d(other%nodes(p)%geometry,tree%nodes(p)%geometry)) return
      end do
      matches=.true.
    end function matching_snapshot
  end subroutine fill_reactive_eb_tree_halo_3d

  subroutine react_reactive_eb_tree_3d(species,reactions,dt,rtol,atol,tree,chemical_change,ok,integrator)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: dt,rtol,atol
    type(eb_reactive_tree_3d), intent(inout) :: tree
    real(dp), intent(inout) :: chemical_change(:)
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: integrator
    type(eb_reactive_tree_3d) :: candidate
    logical, allocatable :: active(:,:,:)
    real(dp) :: before(size(chemical_change)),after(size(chemical_change))
    integer :: node,child,d(3),lo(3),hi(3)
    logical :: valid
    ok=.false.
    if (size(chemical_change)/=reactive_nvar(size(species))) return
    if (.not.all(ieee_is_finite([dt,rtol,atol]))) return
    if (min(dt,rtol,atol)<=0._dp.or.size(reactions)<1) return
    call integrate_reactive_eb_tree_3d(tree,before,valid)
    if (.not.valid) return
    candidate=tree
    do node=1,size(tree%nodes)
      d=tree%topology%nodes(node)%shape
      allocate(active(d(1),d(2),d(3)))
      active=tree%nodes(node)%geometry%volume_fraction>0._dp
      do child=node+1,size(tree%nodes)
        if (tree%topology%nodes(child)%parent/=node) cycle
        call bounds(tree%topology,child,lo,hi)
        active(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.false.
      end do
      if (any(active)) then
        call advance_reactive_chemistry_3d(species,reactions,candidate%nodes(node)%state, &
          candidate%nodes(node)%temperature,d(1),d(2),d(3),dt,rtol,atol,valid,active,integrator)
        if (.not.valid) return
      end if
      deallocate(active)
    end do
    call average_down_reactive_eb_tree_3d(species,candidate,valid)
    if (.not.valid) return
    call integrate_reactive_eb_tree_3d(candidate,after,valid)
    if (.not.valid) return
    if (.not.all(ieee_is_finite(after-before))) return
    tree=candidate;chemical_change=after-before;ok=.true.
  end subroutine react_reactive_eb_tree_3d

  subroutine integrate_reactive_eb_tree_3d(tree,values,ok)
    type(eb_reactive_tree_3d), intent(in) :: tree
    real(dp), intent(inout) :: values(:)
    logical, intent(out) :: ok
    real(dp) :: candidate(size(values)),volume
    logical, allocatable :: exposed(:,:,:)
    integer :: node,child,d(3),lo(3),hi(3),i,j,k
    ok=.false.
    if (.not.valid_tree(tree,size(values))) return
    candidate=0._dp
    do node=1,size(tree%nodes)
      d=tree%topology%nodes(node)%shape
      allocate(exposed(d(1),d(2),d(3)))
      exposed=tree%nodes(node)%geometry%volume_fraction>0._dp
      do child=node+1,size(tree%nodes)
        if (tree%topology%nodes(child)%parent/=node) cycle
        call bounds(tree%topology,child,lo,hi)
        exposed(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.false.
      end do
      volume=tree%nodes(node)%geometry%dx*tree%nodes(node)%geometry%dy*tree%nodes(node)%geometry%dz
      do k=1,d(3)
      do j=1,d(2)
      do i=1,d(1)
        if (.not.exposed(i,j,k)) cycle
        if (.not.all(ieee_is_finite(tree%nodes(node)%state(:,i,j,k)))) return
        candidate=candidate+volume*tree%nodes(node)%geometry%volume_fraction(i,j,k)*tree%nodes(node)%state(:,i,j,k)
      end do
      end do
      end do
      deallocate(exposed)
    end do
    if (.not.all(ieee_is_finite(candidate))) return
    values=candidate;ok=.true.
  end subroutine integrate_reactive_eb_tree_3d
end module eb_reactive_patch_tree_3d_mod
