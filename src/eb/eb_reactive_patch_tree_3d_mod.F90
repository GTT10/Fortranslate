module eb_reactive_patch_tree_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use reactive_3d_mod, only: advance_reactive_chemistry_3d
  use gas_transport_mod, only: gas_transport_species
  use eb_reactive_transport_3d_mod, only: advance_reactive_eb_transport_patch_euler_3d, &
    reactive_eb_transport_patch_timestep_3d
  use reactive_eb_cfl_3d_mod, only: reactive_eb_signal_speeds_3d,compute_reactive_eb_euler_cfl_timestep_3d
  use eb_reactive_hydro_3d_mod, only: advance_reactive_eb_euler_3d
  use eb_reactive_wall_flux_3d_mod, only: reactive_eb_slip_wall_source_3d
  use elementary_kinetics_mod, only: elementary_reaction
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_integer_layout_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d,extract_eb_geometry_3d
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
  public :: remap_reactive_eb_tree_3d
  public :: reactive_eb_tree_timestep_3d
  public :: advance_reactive_eb_tree_full_3d
  public :: advance_reactive_eb_tree_transport_3d
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
      transport_enabled,viscosity,conduction,diffusion,dt,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(eb_reactive_tree_3d), intent(in) :: tree
    real(dp), intent(in) :: hydro_cfl,transport_cfl
    logical, intent(in) :: transport_enabled,viscosity,conduction,diffusion
    real(dp), intent(inout) :: dt
    logical, intent(out) :: ok
    type(eb_reactive_tree_3d) :: snapshot
    real(dp) :: speeds(4),donors(4),diffusivity,maximum_diffusivity,interval,candidate
    integer :: node
    logical :: valid,active
    ok=.false.
    if (.not.all(ieee_is_finite([hydro_cfl,transport_cfl]))) return
    if (hydro_cfl<=0._dp.or.hydro_cfl>1._dp.or.transport_cfl<=0._dp.or.transport_cfl>.5_dp) return
    if (.not.valid_tree(tree,reactive_nvar(size(species)))) return
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
      candidate=min(candidate,interval)
      if (transport_enabled) then
        call reactive_eb_transport_patch_timestep_3d(species,transport,u,t,g,transport_cfl, &
          viscosity,conduction,diffusion,interval,diffusivity,valid,maximum_diffusivity)
        if (.not.valid) return
        ! Conservatively use a full-interval diffusion bound for split halves.
        candidate=min(candidate,interval)
      end if
      end associate
    end do
    dt=candidate;ok=.true.
  end subroutine reactive_eb_tree_timestep_3d

  ! R/2 T/2 H T/2 R/2 composition. Each successful stage remains private
  ! until the complete step and its exchange vectors can be committed.
  subroutine advance_reactive_eb_tree_full_3d(species,reactions,transport,solver,dt,chemistry_enabled, &
      rtol,atol,transport_enabled,viscosity,conduction,diffusion,barodiffusion,tree,boundary_change, &
      wall_change,chemical_change,minimum_theta,ok,integrator)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt,rtol,atol
    logical, intent(in) :: chemistry_enabled,transport_enabled,viscosity,conduction,diffusion,barodiffusion
    type(eb_reactive_tree_3d), intent(inout) :: tree
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemical_change(:),minimum_theta
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: integrator
    type(eb_reactive_tree_3d) :: candidate
    real(dp), allocatable :: boundary(:),wall(:),chemical(:),exchange(:),transport_sum(:)
    real(dp) :: theta,stage_theta
    integer :: nv,half
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    if (.not.ieee_is_finite(dt).or.dt<=0._dp) return
    if (size(boundary_change)/=nv.or.size(wall_change)/=nv.or.size(chemical_change)/=nv) return
    if (.not.valid_tree(tree,nv)) return
    candidate=tree
    allocate(boundary(nv),wall(nv),chemical(nv),exchange(nv),transport_sum(nv))
    chemical=0._dp;transport_sum=0._dp;theta=1._dp
    do half=1,2
      if (half==1.and.chemistry_enabled) then
        call reaction_half(valid)
        if (.not.valid) return
      end if
      if (transport_enabled) then
        call advance_reactive_eb_tree_transport_3d(species,transport,0.5_dp*dt,viscosity,conduction, &
          diffusion,barodiffusion,candidate,exchange,stage_theta,valid)
        if (.not.valid) return
        transport_sum=transport_sum+exchange;theta=min(theta,stage_theta)
      end if
      if (half==1) then
        call advance_reactive_eb_tree_hydro_3d(species,solver,dt,candidate,boundary,wall,valid)
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

  ! Common-interval molecular transport with closed physical boundaries.
  subroutine advance_reactive_eb_tree_transport_3d(species,transport,dt,viscosity,conduction, &
      diffusion,barodiffusion,tree,boundary_change,minimum_theta,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    logical, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), intent(in) :: dt
    type(eb_reactive_tree_3d), intent(inout) :: tree
    real(dp), intent(inout) :: boundary_change(:),minimum_theta
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
      call fill_reactive_eb_tree_halo_3d(species,snapshot,node,ghost,gt,valid)
      if (.not.valid) return
      depth=tree%topology%nodes(node)%level+1
      ! Physical domain faces remain closed to diffusion. Parent patch edges
      ! inside the domain are donor interfaces, even when no peer exists.
      faces=[origin(1,node)>0,origin(1,node)+n(1)<domain(1,depth), &
        origin(2,node)>0,origin(2,node)+n(2)<domain(2,depth), &
        origin(3,node)>0,origin(3,node)+n(3)<domain(3,depth)]
      call advance_reactive_eb_transport_patch_euler_3d(species,transport,snapshot%nodes(node)%geometry, &
        dt,viscosity,conduction,diffusion,barodiffusion,candidate%nodes(node)%state, &
        candidate%nodes(node)%temperature,flux(node),stage_theta,valid,ghost,gt,faces)
      if (.not.valid) return
      theta=min(theta,stage_theta)
      deallocate(ghost,gt)
    end do
    call synchronize_reactive_eb_tree_flux_3d(species,candidate,flux,valid)
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

  ! Serial common-interval Euler hydro, spatial PCM donors and outflow edges.
  ! Every node reads the same averaged snapshot; all outputs commit together.
  subroutine advance_reactive_eb_tree_hydro_3d(species,solver,dt,tree,boundary_change,wall_change,ok)
    type(nasa7_species), intent(in) :: species(:)
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    type(eb_reactive_tree_3d), intent(inout) :: tree
    real(dp), intent(inout) :: boundary_change(:),wall_change(:)
    logical, intent(out) :: ok
    type(eb_reactive_tree_3d) :: snapshot,candidate
    type(eb_integrated_flux_3d), allocatable :: flux(:)
    real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),source(:,:,:,:),boundary(:),wall(:)
    logical, allocatable :: exposed(:,:,:)
    integer :: nv,node,child,n(3),lo(3),hi(3),c
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    if (.not.ieee_is_finite(dt).or.dt<=0._dp) return
    if (size(boundary_change)/=nv.or.size(wall_change)/=nv) return
    if (.not.valid_tree(tree,nv)) return
    snapshot=tree
    call average_down_reactive_eb_tree_3d(species,snapshot,valid)
    if (.not.valid) return
    candidate=snapshot
    allocate(flux(size(tree%nodes)),boundary(nv),wall(nv));wall=0._dp
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      allocate(ghost(nv,0:n(1)+1,0:n(2)+1,0:n(3)+1),gt(0:n(1)+1,0:n(2)+1,0:n(3)+1))
      call fill_reactive_eb_tree_halo_3d(species,snapshot,node,ghost,gt,valid)
      if (.not.valid) return
      call advance_reactive_eb_euler_3d(species,snapshot%nodes(node)%state,snapshot%nodes(node)%temperature, &
        snapshot%nodes(node)%geometry,solver,dt,candidate%nodes(node)%state,candidate%nodes(node)%temperature, &
        valid,ghost,gt,flux(node))
      if (.not.valid) return
      deallocate(ghost,gt)
      allocate(source(nv,n(1),n(2),n(3)),exposed(n(1),n(2),n(3)));exposed=.true.
      do child=node+1,size(tree%nodes)
        if (tree%topology%nodes(child)%parent/=node) cycle
        call bounds(tree%topology,child,lo,hi)
        exposed(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.false.
      end do
      call reactive_eb_slip_wall_source_3d(species,snapshot%nodes(node)%state,snapshot%nodes(node)%temperature, &
        snapshot%nodes(node)%geometry,source,valid)
      if (.not.valid) return
      associate(g=>snapshot%nodes(node)%geometry)
      do c=1,nv
        wall(c)=wall(c)+dt*g%dx*g%dy*g%dz*sum(source(c,:,:,:)*g%volume_fraction,mask=exposed)
      end do
      end associate
      deallocate(source,exposed)
    end do
    call synchronize_reactive_eb_tree_flux_3d(species,candidate,flux,valid)
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
    if (.not.all(ieee_is_finite(boundary)).or..not.all(ieee_is_finite(wall))) return
    tree=candidate;boundary_change=boundary;wall_change=wall;ok=.true.
  end subroutine advance_reactive_eb_tree_hydro_3d

  ! Input states have already accepted each node's raw time-integrated flux.
  ! Replace shared transfers from finest to coarsest on exposed cells only,
  ! then average down. Publish corrected states and effective fluxes together.
  subroutine synchronize_reactive_eb_tree_flux_3d(species,tree,flux,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_reactive_tree_3d), intent(inout) :: tree
    type(eb_integrated_flux_3d), intent(inout) :: flux(:)
    logical, intent(out) :: ok
    type(eb_reactive_tree_3d) :: candidate
    type(eb_integrated_flux_3d), allocatable :: effective(:)
    type(eb_integrated_flux_3d) :: restricted
    type(eb_geometry_3d) :: section
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer(int64) :: offset(3),section_origin(3)
    integer :: nv,node,peer,owner,receiver,axis,depth,parent,child,lo(3),hi(3),d(3),ratio
    logical, allocatable :: refined(:,:,:)
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    if (.not.valid_tree(tree,nv)) return
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
    do depth=tree%topology%level_count()-1,0,-1
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
      if (node==1) cycle
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
    end do
    end do
    call average_down_reactive_eb_tree_3d(species,candidate,valid)
    if (.not.valid) return
    tree=candidate;flux=effective;ok=.true.
  end subroutine synchronize_reactive_eb_tree_flux_3d

  subroutine fill_reactive_eb_tree_halo_3d(species,tree,target,ghost,ghost_temperature,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_reactive_tree_3d), intent(in) :: tree
    integer, intent(in) :: target
    real(dp), intent(inout) :: ghost(:,0:,0:,0:),ghost_temperature(0:,0:,0:)
    logical, intent(out) :: ok
    type(eb_reactive_tree_3d) :: snapshot
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer(int64) :: global(3),mapped(3),local(3),offset(3)
    real(dp), allocatable :: candidate(:,:,:,:),temps(:,:,:),primitive(:)
    real(dp) :: area,guess,recovered,speed
    integer :: n(3),nv,level,depth,node,donor,index(3),inside(3),source(3),first(3),last(3),axis,side,i,j,k
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    if (.not.valid_tree(tree,nv)) return
    if (target<1.or.target>size(tree%nodes)) return
    n=tree%topology%nodes(target)%shape;level=tree%topology%nodes(target)%level+1
    if (any(shape(ghost)/=[nv,n+2]).or.any(shape(ghost_temperature)/=n+2)) return
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if (.not.valid) return
    snapshot=tree
    call average_down_reactive_eb_tree_3d(species,snapshot,valid)
    if (.not.valid) return
    allocate(candidate(nv,0:n(1)+1,0:n(2)+1,0:n(3)+1),temps(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
      primitive(reactive_nprim(size(species))))
    candidate=0._dp;temps=0._dp
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
        if (any(global<0_int64).or.any(global>=domain(:,level))) then
          ! Physical outflow uses the target's adjacent interior cell.
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
        if (.not.ieee_is_finite(guess)) return
        if (guess<=0._dp) return
        call reactive_conserved_to_primitive(species,snapshot%nodes(donor)%state(:,source(1),source(2),source(3)), &
          guess,primitive,recovered,speed,valid)
        if (.not.valid) return
        candidate(:,i,j,k)=snapshot%nodes(donor)%state(:,source(1),source(2),source(3))
        temps(i,j,k)=recovered
      end do
      end do
      end do
    end do
    end do
    ! The peer overlay also verifies open-face aperture and centroid agreement,
    ! including peers that belong to different parent patches.
    do node=1,size(tree%nodes)
      if (node==target.or.tree%topology%nodes(node)%level+1/=level) cycle
      offset=origin(:,node)-origin(:,target)
      if (any(offset < -int(tree%topology%nodes(node)%shape,int64)).or.any(offset>int(n,int64))) cycle
      call overlay_reactive_eb_peer_halo_3d(species,snapshot%nodes(target)%geometry,snapshot%nodes(node)%geometry, &
        int(offset),snapshot%nodes(node)%state,snapshot%nodes(node)%temperature,candidate,temps,valid)
      if (.not.valid) return
    end do
    ghost=candidate;ghost_temperature=temps;ok=.true.
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
