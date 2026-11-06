module eb_reactive_redistribution_3d_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use slope_limiter_mod, only: amrex_limited_slope,amrex_extdir_limited_slope
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d, eb_cut_cell_3d, &
    eb_regular_cell_3d,periodic_eb_faces_match,eb_cell_connectivity_3d
  use eb_reactive_reflux_3d_mod, only: assemble_eb_rereflux_3d
  implicit none
  private

  public :: state_redistribution_graph_slopes_3d,limit_state_recipient_slopes
  public :: state_redistribution_boundary_nodes_3d
  public :: post_eb_spatial_density_rate_3d
  public :: accumulate_state_redist_mol_registration_stage_3d
  public :: compute_state_redist_level_stage_3d
  public :: redistribute_sdc_eb_rate_3d
  public :: reactive_eb_flux_redistribute_3d
  public :: advance_reactive_eb_redistributed_3d
  public :: reactive_eb_weighted_state_redistribute_3d,select_state_redist_tracker_3d
  public :: advance_reactive_eb_state_redistributed_3d

  ! One-based patch coordinates in a prepared geometry/state halo. Periodic
  ! images are explicit halo cells; this box never wraps array indices.
  type, public :: state_redist_box_3d
    integer :: valid_lower(3)=0,valid_upper(3)=0
    integer(int64) :: domain_lower(3)=0_int64,domain_upper(3)=0_int64
    logical :: periodic(3)=.false.
  contains
    procedure :: is_valid => valid_state_redist_box_3d
    procedure :: contains_valid => state_redist_contains_valid_3d
    procedure :: contains_domain => state_redist_contains_domain_3d
  end type state_redist_box_3d

  integer, parameter, public :: srd_coarse_regular=0,srd_coarse_boundary=1,srd_coarse_fine_covered=2
  type, public :: state_redist_register_3d
    integer, allocatable :: coarse_flags(:,:,:)
    logical, allocatable :: fine_ghost(:,:,:)
    ! Valid box ownership is distinct from the level mask: a peer halo is
    ! neither a coarse-interpolated fine ghost nor a valid cell of this patch.
    logical, allocatable :: valid_cells(:,:,:)
    real(dp), allocatable :: coarse_delta(:,:,:,:),fine_delta(:,:,:,:)
    real(dp) :: scale=1._dp
  contains
    procedure :: is_valid => valid_state_redist_register_3d
  end type state_redist_register_3d
  public :: initialize_state_redist_register_3d,accumulate_state_redist_transfer_3d
  public :: assemble_state_redist_interface_rereflux_3d,assemble_state_redist_patch_rereflux_3d
  public :: state_redist_patch_spill_requests_3d
  public :: state_redist_patch_rereflux_budget_3d,state_redist_rereflux_budget_valid_3d
  public :: assemble_state_redist_interface_source_3d,restrict_state_redist_patch_source_3d
  public :: accumulate_state_redist_register_3d,accumulate_state_redist_coarse_register_3d
  public :: initialize_state_redist_tree_register_3d
  public :: initialize_flux_rereflux_register_3d
  public :: restrict_state_redist_fine_delta_3d,build_state_redist_coarse_flags_3d

contains
  ! EBFluxRegister redistributes the coarse/fine flux mismatch even when
  ! the spatial operator uses NoRedist. Reuse the interface ledger's owner,
  ! periodic and connectivity rules for an already integrated density change.
  ! Only coarse boundary cells enter this register; same-level face repairs
  ! and fine-covered averages are not additional coarse/fine inventory.
  subroutine initialize_flux_rereflux_register_3d(tree,node,geometry,source,reg,ok,boundary_faces)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: node
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: source(:,:,:,:)
    type(state_redist_register_3d), intent(inout) :: reg
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(state_redist_register_3d) :: candidate
    integer :: n(3),extra(3),h,i,j,k,nv
    logical :: good
    ok=.false.
    if(.not.tree%is_valid()) return
    if(node<1.or.node>size(tree%nodes)) return
    n=tree%nodes(node)%shape;nv=size(source,1)
    if(nv<1.or.any(shape(source)/=[nv,n])) return
    if(.not.geometry%is_valid().or..not.all(ieee_is_finite(source))) return
    extra=[geometry%nx,geometry%ny,geometry%nz]-n
    if(any(extra<4).or.any(modulo(extra,2)/=0).or.any(extra/=extra(1))) return
    h=extra(1)/2
    call initialize_state_redist_tree_register_3d(tree,node,nv,h,candidate,good,boundary_faces)
    if(.not.good) return
    if(allocated(candidate%coarse_flags)) then
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(candidate%coarse_flags(i+h,j+h,k+h)/=srd_coarse_boundary) cycle
        ! The redistribution mask also includes diagonal neighbors. Flux
        ! registers receive mismatch only through a shared coarse/fine face.
        if(.not.any([candidate%coarse_flags(i+h-1,j+h,k+h),candidate%coarse_flags(i+h+1,j+h,k+h), &
          candidate%coarse_flags(i+h,j+h-1,k+h),candidate%coarse_flags(i+h,j+h+1,k+h), &
          candidate%coarse_flags(i+h,j+h,k+h-1),candidate%coarse_flags(i+h,j+h,k+h+1)]== &
          srd_coarse_fine_covered)) cycle
        candidate%coarse_delta(:,i+h,j+h,k+h)=source(:,i,j,k)
      enddo
      enddo
      enddo
    endif
    reg=candidate;ok=.true.
  end subroutine initialize_flux_rereflux_register_3d

  ! Pinned EB.cpp::pc_post_eb_redistribution, before its optional species
  ! cleaner: covered rates are zero and rho rate is the ordered species sum.
  ! Species, momentum, energy and the already assembled registers are untouched.
  subroutine post_eb_spatial_density_rate_3d(geometry,rate,ok)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(inout) :: rate(:,:,:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: candidate(:,:,:,:)
    real(dp) :: density_rate,term
    integer :: i,j,k,n
    ok=.false.
    if(.not.geometry%is_valid().or.size(rate,1)<6) return
    if(any(shape(rate(1,:,:,:))/=[geometry%nx,geometry%ny,geometry%nz])) return
    candidate=rate
    do k=1,geometry%nz
    do j=1,geometry%ny
    do i=1,geometry%nx
      if(geometry%volume_fraction(i,j,k)<=0) then
        candidate(:,i,j,k)=0
        cycle
      endif
      if(.not.all(ieee_is_finite(candidate(:,i,j,k)))) return
      density_rate=0
      do n=6,size(rate,1)
        term=candidate(n,i,j,k)
        ! Reject overflow before performing the sum under trapping builds.
        if(term>0) then
          if(density_rate>huge(1._dp)-term) return
        else if(term<0) then
          if(density_rate< -huge(1._dp)-term) return
        endif
        density_rate=density_rate+term
      enddo
      candidate(1,i,j,k)=density_rate
    enddo
    enddo
    enddo
    rate=candidate;ok=.true.
  end subroutine post_eb_spatial_density_rate_3d
  logical function valid_state_redist_box_3d(box,dims) result(valid)
    class(state_redist_box_3d), intent(in) :: box
    integer, intent(in) :: dims(3)
    valid=.false.
    if(any(dims<11).or.any(box%valid_lower<6).or.any(box%valid_upper>dims-5)) return
    if(any(box%valid_upper<box%valid_lower)) return
    if(any(box%domain_lower>int(box%valid_lower,int64)).or.any(box%domain_upper<int(box%valid_upper,int64))) return
    if(any(box%domain_lower< -huge(1_int64)+5_int64).or.any(box%domain_upper>huge(1_int64)-5_int64)) return
    valid=.true.
  end function valid_state_redist_box_3d

  pure logical function state_redist_contains_valid_3d(box,index,grow) result(inside)
    class(state_redist_box_3d), intent(in) :: box
    integer, intent(in) :: index(3),grow
    inside=all(int(index,int64)>=int(box%valid_lower,int64)-int(grow,int64)).and. &
      all(int(index,int64)<=int(box%valid_upper,int64)+int(grow,int64))
  end function state_redist_contains_valid_3d

  pure logical function state_redist_contains_domain_3d(box,index,periodic_grow) result(inside)
    class(state_redist_box_3d), intent(in) :: box
    integer, intent(in) :: index(3),periodic_grow
    integer(int64) :: growth(3)
    growth=merge(int(periodic_grow,int64),0_int64,box%periodic)
    inside=all(int(index,int64)>=box%domain_lower-growth).and.all(int(index,int64)<=box%domain_upper+growth)
  end function state_redist_contains_domain_3d

  logical function valid_state_redist_register_3d(reg,nv,dims) result(valid)
    class(state_redist_register_3d), intent(in) :: reg
    integer, intent(in) :: nv,dims(3)
    valid=.false.
    if(nv<1.or.any(dims<1).or..not.ieee_is_finite(reg%scale)) return
    if(reg%scale<0) return
    if(allocated(reg%coarse_flags).neqv.allocated(reg%coarse_delta)) return
    if(allocated(reg%fine_ghost).neqv.allocated(reg%fine_delta)) return
    if(.not.allocated(reg%coarse_flags).and..not.allocated(reg%fine_ghost)) return
    if(allocated(reg%valid_cells)) then
      if(any(shape(reg%valid_cells)/=dims)) return
    endif
    if(allocated(reg%coarse_flags)) then
      if(any(shape(reg%coarse_flags)/=dims).or.any(shape(reg%coarse_delta)/=[nv,dims])) return
      if(any(reg%coarse_flags<0).or.any(reg%coarse_flags>2)) return
      if(any(.not.ieee_is_finite(reg%coarse_delta))) return
    endif
    if(allocated(reg%fine_ghost)) then
      if(any(shape(reg%fine_ghost)/=dims).or.any(shape(reg%fine_delta)/=[nv,dims])) return
      if(allocated(reg%valid_cells)) then
        if(any(reg%valid_cells.and.reg%fine_ghost)) return
      endif
      if(any(.not.ieee_is_finite(reg%fine_delta))) return
    endif
    valid=.true.
  end function valid_state_redist_register_3d

  ! Add an accepted stage/subinterval ledger exactly once. The scale is
  ! provenance of already weighted deltas; accumulation never reapplies it.
  ! Fixed MLStateRedistribute writes coarse data directly on every source
  ! evaluation. Provisional fine dm is added only when its flux is accepted.
  ! Preserve destination fine data while carrying prior coarse evaluations.
  ! Pinned Diffusion.cpp/Advance.cpp: coarse deltaR is written directly on
  ! every MOL evaluation; FineAdd receives only the old and final endpoint dm.
  ! Each source register already contains fac_for_redist=.5. Do not rescale it.
  ! Advance.cpp evaluates stages in increasing order, including zero-flux
  ! provisional stages. Keep that order rather than summing them in reverse.
  subroutine accumulate_state_redist_mol_registration_stage_3d(destination,source,stage,last_stage,nv,dims,ok)
    type(state_redist_register_3d), intent(inout) :: destination
    type(state_redist_register_3d), intent(in) :: source
    integer, intent(in) :: stage,last_stage,nv,dims(3)
    logical, intent(out) :: ok
    ok=.false.
    if(last_stage<1.or.stage<0.or.stage>last_stage) return
    if(.not.source%is_valid(nv,dims).or.source%scale/=.5_dp) return
    if(.not.standard_bounds(source)) return
    if(stage>0) then
      if(.not.destination%is_valid(nv,dims).or..not.standard_bounds(destination)) return
    endif
    if(stage==0) then
      destination=source;ok=.true.
    else if(stage==last_stage) then
      call accumulate_state_redist_register_3d(destination,source,nv,dims,ok)
    else
      call accumulate_state_redist_coarse_register_3d(destination,source,nv,dims,ok)
    endif
  contains
    logical function standard_bounds(reg) result(good)
      type(state_redist_register_3d), intent(in) :: reg
      good=.false.
      if(allocated(reg%coarse_flags)) then
        if(.not.allocated(reg%coarse_delta)) return
        if(any(lbound(reg%coarse_flags)/=1).or.any(lbound(reg%coarse_delta)/=1)) return
      endif
      if(allocated(reg%fine_ghost)) then
        if(.not.allocated(reg%fine_delta)) return
        if(any(lbound(reg%fine_ghost)/=1).or.any(lbound(reg%fine_delta)/=1)) return
      endif
      if(allocated(reg%valid_cells)) then
        if(any(lbound(reg%valid_cells)/=1)) return
      endif
      good=.true.
    end function standard_bounds
  end subroutine accumulate_state_redist_mol_registration_stage_3d

  subroutine accumulate_state_redist_coarse_register_3d(destination,source,nv,dims,ok)
    type(state_redist_register_3d), intent(inout) :: destination
    type(state_redist_register_3d), intent(in) :: source
    integer, intent(in) :: nv,dims(3)
    logical, intent(out) :: ok
    type(state_redist_register_3d) :: coarse_only
    ok=.false.
    if(.not.source%is_valid(nv,dims)) return
    coarse_only=source
    if(allocated(coarse_only%fine_delta)) coarse_only%fine_delta=0
    call accumulate_state_redist_register_3d(destination,coarse_only,nv,dims,ok)
  end subroutine accumulate_state_redist_coarse_register_3d

  subroutine accumulate_state_redist_register_3d(destination,source,nv,dims,ok)
    type(state_redist_register_3d), intent(inout) :: destination
    type(state_redist_register_3d), intent(in) :: source
    integer, intent(in) :: nv,dims(3)
    logical, intent(out) :: ok
    type(state_redist_register_3d) :: candidate
    real(dp) :: delta(nv)
    integer :: i,j,k
    logical :: good
    ok=.false.
    if(.not.destination%is_valid(nv,dims).or..not.source%is_valid(nv,dims)) return
    if(destination%scale/=source%scale) return
    if(allocated(destination%coarse_flags).neqv.allocated(source%coarse_flags)) return
    if(allocated(destination%fine_ghost).neqv.allocated(source%fine_ghost)) return
    if(allocated(destination%valid_cells).neqv.allocated(source%valid_cells)) return
    if(allocated(source%coarse_flags)) then
      if(any(destination%coarse_flags/=source%coarse_flags)) return
    endif
    if(allocated(source%fine_ghost)) then
      if(any(destination%fine_ghost.neqv.source%fine_ghost)) return
    endif
    if(allocated(source%valid_cells)) then
      if(any(destination%valid_cells.neqv.source%valid_cells)) return
    endif
    candidate=destination
    do k=1,dims(3)
    do j=1,dims(2)
    do i=1,dims(1)
      if(allocated(source%coarse_delta)) then
        delta=source%coarse_delta(:,i,j,k)
        call add_register_increment(candidate%coarse_delta(:,i,j,k),delta,good)
        if(.not.good) return
        candidate%coarse_delta(:,i,j,k)=delta
      endif
      if(allocated(source%fine_delta)) then
        delta=source%fine_delta(:,i,j,k)
        call add_register_increment(candidate%fine_delta(:,i,j,k),delta,good)
        if(.not.good) return
        candidate%fine_delta(:,i,j,k)=delta
      endif
    enddo
    enddo
    enddo
    destination=candidate;ok=.true.
  end subroutine accumulate_state_redist_register_3d

  subroutine initialize_state_redist_register_3d(nv,dims,reg,ok,coarse_flags,fine_ghost,scale,valid_cells)
    integer, intent(in) :: nv,dims(3)
    type(state_redist_register_3d), intent(inout) :: reg
    logical, intent(out) :: ok
    integer, optional, intent(in) :: coarse_flags(:,:,:)
    logical, optional, intent(in) :: fine_ghost(:,:,:),valid_cells(:,:,:)
    real(dp), optional, intent(in) :: scale
    type(state_redist_register_3d) :: candidate
    integer :: status
    ok=.false.
    if(nv<1.or.any(dims<1)) return
    if(present(scale)) candidate%scale=scale
    if(present(valid_cells)) then
      if(any(shape(valid_cells)/=dims)) return
      candidate%valid_cells=valid_cells
    endif
    if(present(coarse_flags)) then
      if(any(shape(coarse_flags)/=dims)) return
      candidate%coarse_flags=coarse_flags
      allocate(candidate%coarse_delta(nv,dims(1),dims(2),dims(3)),source=0._dp,stat=status)
      if(status/=0) return
    endif
    if(present(fine_ghost)) then
      if(any(shape(fine_ghost)/=dims)) return
      candidate%fine_ghost=fine_ghost
      allocate(candidate%fine_delta(nv,dims(1),dims(2),dims(3)),source=0._dp,stat=status)
      if(status/=0) return
    endif
    if(.not.candidate%is_valid(nv,dims)) return
    reg=candidate;ok=.true.
  end subroutine initialize_state_redist_register_3d

  ! Construct the masks for a patch at one AMR level, including its halo.
  ! Refined coverage and fine-ghost membership use the union of all level
  ! patches, so ownership changes and same-level peer halos do not create
  ! artificial coarse/fine exchanges. Global indices are zero-based int64.
  subroutine initialize_state_redist_tree_register_3d(tree,node,nv,halo_width,reg,ok,boundary_faces,scale,box)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: node,nv,halo_width
    type(state_redist_register_3d), intent(inout) :: reg
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: scale
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:),flags(:,:,:)
    logical, allocatable :: ghost(:,:,:),valid_cells(:,:,:)
    type(state_redist_box_3d), optional, intent(inout) :: box
    type(state_redist_box_3d) :: candidate_box
    type(state_redist_register_3d) :: candidate
    integer :: level,n(3),dims(3),i,j,k,ii,jj,kk,ios
    integer(int64) :: point(3),neighbor(3),h
    logical :: periodic(3),good,inside_domain,as_coarse,as_fine
    ok=.false.
    if(.not.tree%is_valid().or.nv<1.or.halo_width<0) return
    if(node<1.or.node>size(tree%nodes)) return
    n=tree%nodes(node)%shape;level=tree%nodes(node)%level
    if(halo_width>(huge(1)-maxval(n))/2) return
    dims=n+2*halo_width;h=int(halo_width,int64)
    periodic=.false.
    if(present(boundary_faces)) then
      if(any(boundary_faces/='outflow'.and.boundary_faces/='inflow'.and.boundary_faces/='periodic'.and. &
        boundary_faces/='slip_wall'.and.boundary_faces/='no_slip_wall'.and.boundary_faces/='isothermal_wall')) return
      do i=1,3
        if((boundary_faces(2*i-1)=='periodic').neqv.(boundary_faces(2*i)=='periodic')) return
        periodic(i)=boundary_faces(2*i-1)=='periodic'
      enddo
    endif
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,good)
    if(.not.good) return
    if(any(domain(:,level+1)>huge(1_int64)-h-1_int64)) return
    candidate_box%valid_lower=halo_width+1;candidate_box%valid_upper=halo_width+n
    candidate_box%domain_lower=1_int64+h-origin(:,node)
    candidate_box%domain_upper=domain(:,level+1)+h-origin(:,node)
    candidate_box%periodic=periodic
    if(present(box)) then
      if(.not.candidate_box%is_valid(dims)) return
    endif
    as_fine=level>0;as_coarse=level<tree%level_count()-1
    ! A one-level tree has no exchange but still returns a usable zero ledger.
    if(.not.as_fine) as_coarse=.true.
    allocate(valid_cells(dims(1),dims(2),dims(3)),stat=ios)
    if(ios/=0) return
    valid_cells=.false.
    valid_cells(halo_width+1:halo_width+n(1),halo_width+1:halo_width+n(2),halo_width+1:halo_width+n(3))=.true.
    if(as_fine) then
      allocate(ghost(dims(1),dims(2),dims(3)),source=.false.,stat=ios)
      if(ios/=0) return
    endif
    if(as_coarse) then
      allocate(flags(dims(1),dims(2),dims(3)),source=srd_coarse_regular,stat=ios)
      if(ios/=0) return
    endif
    do k=1,dims(3)
    do j=1,dims(2)
    do i=1,dims(1)
      point=origin(:,node)+int([i,j,k]-1,int64)-h
      call canonical_point(point,inside_domain)
      if(.not.inside_domain) cycle
      if(as_fine) ghost(i,j,k)=.not.in_level(point)
      if(.not.as_coarse) cycle
      if(is_refined(point)) then
        flags(i,j,k)=srd_coarse_fine_covered;cycle
      endif
      do kk=-1,1
      do jj=-1,1
      do ii=-1,1
        if(ii==0.and.jj==0.and.kk==0) cycle
        neighbor=point+int([ii,jj,kk],int64)
        call canonical_point(neighbor,inside_domain)
        if(.not.inside_domain) cycle
        if(is_refined(neighbor)) flags(i,j,k)=srd_coarse_boundary
      enddo
      enddo
      enddo
    enddo
    enddo
    enddo
    call initialize_state_redist_register_3d(nv,dims,candidate,good,flags,ghost,scale,valid_cells)
    if(.not.good) return
    reg=candidate
    if(present(box)) box=candidate_box
    ok=.true.
  contains
    subroutine canonical_point(index,inside)
      integer(int64), intent(inout) :: index(3)
      logical, intent(out) :: inside
      integer :: a
      inside=.true.
      do a=1,3
        if(periodic(a)) then
          index(a)=modulo(index(a),domain(a,level+1))
        else if(index(a)<0.or.index(a)>=domain(a,level+1)) then
          inside=.false.
        endif
      enddo
    end subroutine canonical_point
    logical function in_level(index) result(found)
      integer(int64), intent(in) :: index(3)
      integer :: p
      found=.false.
      do p=1,size(tree%nodes)
        if(tree%nodes(p)%level/=level) cycle
        if(all(index>=origin(:,p)).and.all(index<origin(:,p)+int(tree%nodes(p)%shape,int64))) then
          found=.true.;return
        endif
      enddo
    end function in_level
    logical function is_refined(index) result(found)
      integer(int64), intent(in) :: index(3)
      integer(int64) :: lower(3),upper(3),ratio
      integer :: p
      found=.false.
      do p=2,size(tree%nodes)
        if(tree%nodes(p)%level/=level+1) cycle
        ratio=int(tree%nodes(p)%patch%refinement_ratio,int64)
        lower=origin(:,p)/ratio;upper=lower+int(tree%nodes(p)%shape,int64)/ratio
        if(all(index>=lower).and.all(index<upper)) then
          found=.true.;return
        endif
      enddo
    end function is_refined
  end subroutine initialize_state_redist_tree_register_3d

  ! Prepared halos contain the spatial RHS evaluated at this level, including
  ! its ghost cells. This routine does not interpolate a coarse spatial RHS.
  ! Registers describe this stage only; interval accumulation belongs to the
  ! timestep caller. Rate and registers commit together after all owned patches.
  subroutine compute_state_redist_level_stage_3d(tree,geometry,old_halos,rhs_halos,nv,depth,dt,rate,registers,ok, &
      boundary_faces,target_volume_fraction,max_order,register_scale,owned,update_scales)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: old_halos,rhs_halos
    type(amr_tree_state_3d), optional, intent(in) :: update_scales
    integer, intent(in) :: nv,depth
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: rate
    type(state_redist_register_3d), allocatable, intent(inout) :: registers(:)
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: target_volume_fraction,register_scale
    integer, optional, intent(in) :: max_order
    logical, optional, intent(in) :: owned(:)
    type(amr_tree_state_3d) :: pending_rate
    type(state_redist_register_3d), allocatable :: pending_registers(:)
    type(state_redist_box_3d) :: box
    real(dp), allocatable :: patch_rate(:,:,:,:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer :: count,p,n(3),dims(3),extra(3),h,root_h,root_extra(3),ios
    real(dp) :: root_dx(3),domain_lower(3),expected_dx(3),expected_lower(3),actual_lower(3),tolerance(3)
    logical :: good
    ok=.false.
    if(.not.tree%is_valid().or.nv<1.or..not.ieee_is_finite(dt)) return
    if(dt<=0.or.depth<0.or.depth>=tree%level_count()) return
    if(present(max_order)) then
      if(max_order/=0.and.max_order/=2.and.max_order/=4) return
    endif
    if(present(target_volume_fraction)) then
      if(.not.ieee_is_finite(target_volume_fraction)) return
      if(target_volume_fraction<=0.or.target_volume_fraction>1) return
    endif
    if(present(register_scale)) then
      if(.not.ieee_is_finite(register_scale)) return
      if(register_scale<0) return
    endif
    if(present(boundary_faces)) then
      if(any(boundary_faces/='outflow'.and.boundary_faces/='inflow'.and.boundary_faces/='periodic'.and. &
        boundary_faces/='slip_wall'.and.boundary_faces/='no_slip_wall'.and.boundary_faces/='isothermal_wall')) return
      do p=1,3
        if((boundary_faces(2*p-1)=='periodic').neqv.(boundary_faces(2*p)=='periodic')) return
      enddo
    endif
    count=size(tree%nodes)
    if(size(geometry)/=count) return
    if(.not.allocated(old_halos%fields).or..not.allocated(rhs_halos%fields).or..not.allocated(rate%fields)) return
    if(size(old_halos%fields)/=count.or.size(rhs_halos%fields)/=count.or.size(rate%fields)/=count) return
    if(present(update_scales)) then
      if(.not.allocated(update_scales%fields)) return
      if(size(update_scales%fields)/=count) return
    endif
    if(present(owned)) then
      if(size(owned)/=count) return
    endif
    if(allocated(registers)) then
      if(size(registers)/=count) return
    endif
    ! Geometry is replicated metadata, including on idle MPI ranks. Check its
    ! physical origin against the integer AMR layout before using halo indices.
    if(.not.geometry(1)%is_valid()) return
    root_extra=[geometry(1)%nx,geometry(1)%ny,geometry(1)%nz]-tree%nodes(1)%shape
    if(any(root_extra<10).or.any(modulo(root_extra,2)/=0)) return
    if(any(root_extra/=root_extra(1))) return
    root_h=root_extra(1)/2
    root_dx=[geometry(1)%dx,geometry(1)%dy,geometry(1)%dz]
    domain_lower=[geometry(1)%x_lower,geometry(1)%y_lower,geometry(1)%z_lower]+real(root_h,dp)*root_dx
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,good)
    if(.not.good) return
    pending_rate=rate
    if(allocated(registers)) then
      pending_registers=registers
    else
      allocate(pending_registers(count),stat=ios)
      if(ios/=0) return
    endif
    do p=1,count
      if(tree%nodes(p)%level/=depth) cycle
      if(.not.geometry(p)%is_valid()) return
      n=tree%nodes(p)%shape;dims=[geometry(p)%nx,geometry(p)%ny,geometry(p)%nz]
      extra=dims-n
      if(any(extra<10).or.any(modulo(extra,2)/=0)) return
      if(any(extra/=extra(1))) return
      h=extra(1)/2
      expected_dx=root_dx*(real(domain(:,1),dp)/real(domain(:,depth+1),dp))
      actual_lower=[geometry(p)%x_lower,geometry(p)%y_lower,geometry(p)%z_lower]
      expected_lower=domain_lower+(real(origin(:,p),dp)-real(h,dp))*expected_dx
      tolerance=256._dp*epsilon(1._dp)*max(1._dp,abs(expected_lower),abs(domain_lower),abs(expected_dx))
      if(any(abs(actual_lower-expected_lower)>tolerance)) return
      if(any(abs([geometry(p)%dx,geometry(p)%dy,geometry(p)%dz]-expected_dx)> &
        256._dp*epsilon(1._dp)*expected_dx)) return
      if(present(owned)) then
        if(.not.owned(p)) cycle
      endif
      if(.not.allocated(old_halos%fields(p)%values).or..not.allocated(rhs_halos%fields(p)%values)) return
      if(.not.allocated(rate%fields(p)%values)) return
      if(any(shape(old_halos%fields(p)%values)/=[nv,dims])) return
      if(any(shape(rhs_halos%fields(p)%values)/=[nv,dims])) return
      if(any(shape(rate%fields(p)%values)/=[nv,n])) return
      call initialize_state_redist_tree_register_3d(tree,p,nv,h,pending_registers(p),good, &
        boundary_faces,register_scale,box)
      if(.not.good) return
      allocate(patch_rate(nv,dims(1),dims(2),dims(3)),source=0._dp,stat=ios)
      if(ios/=0) return
      if(present(update_scales)) then
        if(.not.allocated(update_scales%fields(p)%values)) return
        if(any(shape(update_scales%fields(p)%values)/=[1,dims])) return
        call redistribute_sdc_eb_rate_3d(geometry(p),old_halos%fields(p)%values,rhs_halos%fields(p)%values, &
          dt,'state_redist',patch_rate,good,target_volume_fraction,max_order,boundary_faces=boundary_faces, &
          registers=pending_registers(p),box=box,update_scale=update_scales%fields(p)%values(1,:,:,:))
      else
      call redistribute_sdc_eb_rate_3d(geometry(p),old_halos%fields(p)%values,rhs_halos%fields(p)%values, &
        dt,'state_redist',patch_rate,good,target_volume_fraction,max_order,boundary_faces=boundary_faces, &
        registers=pending_registers(p),box=box)
      endif
      if(.not.good) return
      pending_rate%fields(p)%values=patch_rate(:,h+1:h+n(1),h+1:h+n(2),h+1:h+n(3))
      deallocate(patch_rate)
    enddo
    rate=pending_rate
    call move_alloc(pending_registers,registers)
    ok=.true.
  end subroutine compute_state_redist_level_stage_3d

  ! Re-reflux accepted interface ledgers into coarse density increments.
  ! Fine-covered coarse recipients remain present: the tree flux caller
  ! injects them into fine gas before joint EOS recovery/average-down.
  subroutine assemble_state_redist_interface_rereflux_3d(tree,geometry,registers,nv,depth,increment,ok, &
      boundary_faces,spill,spill_active)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(state_redist_register_3d), intent(in) :: registers(:)
    integer, intent(in) :: nv,depth
    type(amr_tree_state_3d), intent(inout) :: increment
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(amr_tree_state_3d), optional, intent(inout) :: spill
    logical, optional, intent(out) :: spill_active
    type(amr_tree_state_3d) :: source,pending,pending_spill
    type(state_redist_box_3d) :: box
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer :: h
    logical :: periodic(3)
    real(dp), allocatable :: patch_increment(:,:,:,:),halo_spill(:,:,:,:),values(:,:)
    integer, allocatable :: requests(:,:)
    real(dp) :: budget(3,nv),patch_budget(3,nv),spill_budget(2,nv),delta(nv)
    integer :: p,n(3),q,index(3),r
    logical :: good
    ok=.false.
    if(present(spill_active)) spill_active=.false.
    call assemble_state_redist_interface_source_3d(tree,geometry,registers,nv,depth,source,good,boundary_faces)
    if(.not.good) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,good)
    if(.not.good) return
    periodic=.false.
    if(present(boundary_faces)) periodic=[boundary_faces(1)=='periodic',boundary_faces(3)=='periodic', &
      boundary_faces(5)=='periodic']
    allocate(pending%fields(size(tree%nodes)))
    if(present(spill)) then
      allocate(pending_spill%fields(size(tree%nodes)))
      do p=1,size(tree%nodes)
        n=tree%nodes(p)%shape
        allocate(pending_spill%fields(p)%values(nv,n(1),n(2),n(3)),source=0._dp)
      enddo
    endif
    budget=0._dp
    do p=1,size(tree%nodes)
      n=tree%nodes(p)%shape
      allocate(pending%fields(p)%values(nv,n(1),n(2),n(3)),source=0._dp)
      if(tree%nodes(p)%level/=depth) cycle
      h=(geometry(p)%nx-n(1))/2
      box%valid_lower=h+1;box%valid_upper=h+n
      box%domain_lower=1_int64+int(h,int64)-origin(:,p)
      box%domain_upper=domain(:,depth+1)+int(h,int64)-origin(:,p)
      box%periodic=periodic
      if(present(spill).and.depth>0) then
        call assemble_state_redist_patch_rereflux_3d(geometry(p),registers(p),n,source%fields(p)%values, &
          patch_increment,good,halo_spill,box)
      else
        call assemble_state_redist_patch_rereflux_3d(geometry(p),registers(p),n,source%fields(p)%values, &
          patch_increment,good,box=box)
      endif
      if(.not.good) return
      call state_redist_patch_rereflux_budget_3d(geometry(p),n,source%fields(p)%values, &
        patch_increment,patch_budget,good)
      if(.not.good) return
      budget=budget+patch_budget
      if(present(spill).and.depth>0) then
        call state_redist_patch_spill_requests_3d(tree,p,geometry,halo_spill,requests,values,spill_budget,good,boundary_faces)
        if(.not.good) return
        budget(2,:)=budget(2,:)+spill_budget(1,:)
        budget(3,:)=budget(3,:)+spill_budget(2,:)
        do r=1,size(requests,2)
          q=requests(1,r);index=requests(2:4,r);delta=values(:,r)
          call add_register_increment(pending_spill%fields(q)%values(:,index(1),index(2),index(3)),delta,good)
          if(.not.good) return
          pending_spill%fields(q)%values(:,index(1),index(2),index(3))=delta
        enddo
        if(present(spill_active)) spill_active=spill_active.or.size(requests,2)>0
      endif
      call move_alloc(patch_increment,pending%fields(p)%values)
    enddo
    if(.not.state_redist_rereflux_budget_valid_3d(budget)) return
    if(present(spill)) call move_alloc(pending_spill%fields,spill%fields)
    call move_alloc(pending%fields,increment%fields);ok=.true.
  end subroutine assemble_state_redist_interface_rereflux_3d

  ! Shared owner-patch operator. Connectivity follows monotone 27-neighbor
  ! paths through open EB faces; fine-covered coarse recipients are retained.
  subroutine assemble_state_redist_patch_rereflux_3d(geometry,reg,n,source,increment,ok,halo_spill,box)
    type(eb_geometry_3d), intent(in) :: geometry
    type(state_redist_box_3d), optional, intent(in) :: box
    type(state_redist_register_3d), intent(in) :: reg
    integer, intent(in) :: n(3)
    real(dp), intent(in) :: source(:,:,:,:)
    real(dp), allocatable, intent(inout) :: increment(:,:,:,:)
    logical, intent(out) :: ok
    real(dp), allocatable, optional, intent(inout) :: halo_spill(:,:,:,:)
    real(dp), allocatable :: patch_increment(:,:,:,:),pending(:,:,:,:),full_increment(:,:,:,:)
    logical, allocatable :: connected(:,:,:,:,:,:),boundary(:,:,:),owner_roots(:,:,:)
    integer :: dims(3),extra(3),h,i,j,k,nv,ii,jj,kk
    logical :: good
    ok=.false.;nv=size(source,1);dims=[geometry%nx,geometry%ny,geometry%nz];extra=dims-n
    if(nv<1.or.any(n<1).or..not.geometry%is_valid()) return
    if(any(extra<4).or.any(modulo(extra,2)/=0).or.any(extra/=extra(1))) return
    if(any(shape(source)/=[nv,dims]).or..not.reg%is_valid(nv,dims)) return
    if(.not.allocated(reg%coarse_flags).or..not.all(ieee_is_finite(source))) return
    h=extra(1)/2
    if(present(box)) then
      if(.not.box%is_valid(dims)) return
      if(any(box%valid_lower/=h+1).or.any(box%valid_upper/=h+n)) return
    endif
    allocate(connected(3,3,3,dims(1),dims(2),dims(3)),source=.false.)
    boundary=reg%coarse_flags==srd_coarse_boundary
    do k=h,h+n(3)+1
    do j=h,h+n(2)+1
    do i=h,h+n(1)+1
      if(.not.boundary(i,j,k).or.all(source(:,i,j,k)==0)) cycle
      call eb_cell_connectivity_3d(geometry,[i,j,k],connected(:,:,:,i,j,k),good)
      if(.not.good) return
      ! Physical-domain ghost cells have geometric EB connectivity, but
      ! cannot receive conserved register inventory. Retain periodic images
      ! and all peer/ancestor halo recipients; restrict only nonperiodic domain
      ! faces before the unchanged re-reflux kernel normalizes its weights.
      if(present(box)) then
        do kk=-1,1
        do jj=-1,1
        do ii=-1,1
          if(.not.box%contains_domain([i+ii,j+jj,k+kk],2)) connected(ii+2,jj+2,kk+2,i,j,k)=.false.
        enddo
        enddo
        enddo
      endif
    enddo
    enddo
    enddo
    call assemble_eb_rereflux_3d(geometry,source,boundary,connected,patch_increment,good,[h,h,h]+1,n+h)
    if(.not.good) return
    if(present(halo_spill)) then
      if(.not.allocated(reg%valid_cells)) return
      owner_roots=boundary.and.reg%valid_cells
      call assemble_eb_rereflux_3d(geometry,source,owner_roots,connected,full_increment,good)
      if(.not.good) return
      call move_alloc(full_increment,halo_spill)
    endif
    pending=patch_increment(:,h+1:h+n(1),h+1:h+n(2),h+1:h+n(3))
    call move_alloc(pending,increment);ok=.true.
  end subroutine assemble_state_redist_patch_rereflux_3d

  ! Project the one-cell halo portion of owner-root re-reflux onto an exposed
  ! cell of the next coarser level. Same-level peer and periodic recipients
  ! are already assembled by their owner patch and are excluded here.
  subroutine state_redist_patch_spill_requests_3d(tree,p,geometry,halo_spill,requests,values,budget,ok,boundary_faces)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: p
    type(eb_geometry_3d), intent(in) :: geometry(:)
    real(dp), intent(in) :: halo_spill(:,:,:,:)
    integer, allocatable, intent(inout) :: requests(:,:)
    real(dp), allocatable, intent(inout) :: values(:,:)
    real(dp), intent(out) :: budget(:,:)
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer(int64) :: cell(3),parent_cell(3)
    integer, allocatable :: ratios(:),pending_requests(:,:)
    real(dp), allocatable :: pending_values(:,:)
    real(dp) :: fine_mass,coarse_mass,scale,contribution
    integer :: n(3),dims(3),h,depth,ratio,capacity,count,i,j,k,q,axis,c,index(3),parent_h,nv
    logical :: periodic(3),inside,found,good
    ok=.false.;budget=0._dp
    if(.not.tree%is_valid().or.p<1.or.p>size(tree%nodes)) return
    if(size(geometry)/=size(tree%nodes)) return
    depth=tree%nodes(p)%level;nv=size(halo_spill,1)
    if(depth<1.or.nv<1.or.any(shape(budget)/=[2,nv])) return
    n=tree%nodes(p)%shape;dims=[geometry(p)%nx,geometry(p)%ny,geometry(p)%nz]
    if(.not.geometry(p)%is_valid().or.any(dims-n<4).or.any(modulo(dims-n,2)/=0)) return
    if(any(dims-n/=dims(1)-n(1)).or.any(shape(halo_spill)/=[nv,dims])) return
    if(any(.not.ieee_is_finite(halo_spill))) return
    h=(dims(1)-n(1))/2;ratio=tree%nodes(p)%patch%refinement_ratio
    if(ratio<1) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,good)
    if(.not.good) return
    periodic=.false.
    if(present(boundary_faces)) then
      do axis=1,3
        if((boundary_faces(2*axis-1)=='periodic').neqv.(boundary_faces(2*axis)=='periodic')) return
        periodic(axis)=boundary_faces(2*axis-1)=='periodic'
      enddo
    endif
    capacity=product(n+2)-product(n)
    allocate(pending_requests(4,capacity),pending_values(nv,capacity));count=0
    do k=h,h+n(3)+1
    do j=h,h+n(2)+1
    do i=h,h+n(1)+1
      if(all([i,j,k]>=h+1).and.all([i,j,k]<=h+n)) cycle
      if(all(halo_spill(:,i,j,k)==0._dp)) cycle
      cell=origin(:,p)+int([i,j,k]-h-1,int64);inside=.true.
      do axis=1,3
        if(periodic(axis)) then
          cell(axis)=modulo(cell(axis),domain(axis,depth+1))
        else if(cell(axis)<0.or.cell(axis)>=domain(axis,depth+1)) then
          inside=.false.
        endif
      enddo
      if(.not.inside) cycle
      found=.false.
      do q=1,size(tree%nodes)
        if(tree%nodes(q)%level/=depth) cycle
        if(any(cell<origin(:,q)).or.any(cell>=origin(:,q)+int(tree%nodes(q)%shape,int64))) cycle
        found=.true.;exit
      enddo
      if(found) cycle
      parent_cell=(cell-modulo(cell,int(ratio,int64)))/int(ratio,int64)
      q=0
      do axis=1,size(tree%nodes)
        if(tree%nodes(axis)%level/=depth-1) cycle
        if(any(parent_cell<origin(:,axis)).or. &
          any(parent_cell>=origin(:,axis)+int(tree%nodes(axis)%shape,int64))) cycle
        q=axis;exit
      enddo
      if(q==0.or..not.geometry(q)%is_valid()) return
      index=int(parent_cell-origin(:,q))+1
      dims=[geometry(q)%nx,geometry(q)%ny,geometry(q)%nz]
      if(any(modulo(dims-tree%nodes(q)%shape,2)/=0)) return
      parent_h=(dims(1)-tree%nodes(q)%shape(1))/2
      fine_mass=geometry(p)%volume_fraction(i,j,k)*geometry(p)%dx*geometry(p)%dy*geometry(p)%dz
      coarse_mass=geometry(q)%volume_fraction(index(1)+parent_h,index(2)+parent_h,index(3)+parent_h)* &
        geometry(q)%dx*geometry(q)%dy*geometry(q)%dz
      if(fine_mass<=0._dp.or.coarse_mass<=0._dp) return
      scale=fine_mass/coarse_mass
      if(.not.ieee_is_finite(scale)) return
      count=count+1;pending_requests(:,count)=[q,index]
      do c=1,nv
        contribution=halo_spill(c,i,j,k)*scale
        if(.not.ieee_is_finite(contribution)) return
        pending_values(c,count)=contribution
        budget(1,c)=budget(1,c)+fine_mass*halo_spill(c,i,j,k)
        budget(2,c)=budget(2,c)+abs(fine_mass*halo_spill(c,i,j,k))
      enddo
    enddo
    enddo
    enddo
    if(any(.not.ieee_is_finite(budget))) return
    requests=pending_requests(:,:count);values=pending_values(:,:count);ok=.true.
  end subroutine state_redist_patch_spill_requests_3d

  ! A cut-cell re-reflux may scatter into a coarse-level halo with no valid
  ! cell at this level. In that case the clipped patch increment is not a
  ! conservative interface correction and must not be accepted as a step.
  subroutine state_redist_patch_rereflux_budget_3d(geometry,n,source,increment,budget,ok)
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: n(3)
    real(dp), intent(in) :: source(:,:,:,:),increment(:,:,:,:)
    real(dp), intent(out) :: budget(:,:)
    logical, intent(out) :: ok
    real(dp) :: volume,src,delta,kappa
    integer :: nv,dims(3),h,c,i,j,k
    ok=.false.;budget=0._dp
    nv=size(source,1);dims=[geometry%nx,geometry%ny,geometry%nz]
    if(.not.geometry%is_valid().or.nv<1.or.any(n<1)) return
    if(any(dims-n<4).or.any(modulo(dims-n,2)/=0).or.any(dims-n/=dims(1)-n(1))) return
    if(any(shape(source)/=[nv,dims]).or.any(shape(increment)/=[nv,n]).or.any(shape(budget)/=[3,nv])) return
    if(any(.not.ieee_is_finite(source)).or.any(.not.ieee_is_finite(increment))) return
    h=(dims(1)-n(1))/2;volume=geometry%dx*geometry%dy*geometry%dz
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      kappa=geometry%volume_fraction(i+h,j+h,k+h)
      if(kappa<=0._dp) cycle
      do c=1,nv
        src=kappa*volume*source(c,i+h,j+h,k+h)
        delta=kappa*volume*increment(c,i,j,k)
        budget(1,c)=budget(1,c)+src
        budget(2,c)=budget(2,c)+delta
        budget(3,c)=budget(3,c)+abs(src)+abs(delta)
      enddo
    enddo
    enddo
    enddo
    ok=all(ieee_is_finite(budget))
  end subroutine state_redist_patch_rereflux_budget_3d

  pure logical function state_redist_rereflux_budget_valid_3d(budget) result(valid)
    real(dp), intent(in) :: budget(:,:)
    valid=.false.
    if(size(budget,1)/=3.or.size(budget,2)<1) return
    if(any(.not.ieee_is_finite(budget))) return
    valid=all(abs(budget(1,:)-budget(2,:))<= &
      128._dp*epsilon(1._dp)*max(tiny(1._dp),budget(3,:)))
  end function state_redist_rereflux_budget_valid_3d

  ! Assemble an interval's signed coarse interface source before EB re-reflux.
  ! Fine dm is clipped to valid+1, but the normalization uses the whole
  ! fine geometric block, as EBFluxRegister::FineAdd does. No gas is updated.
  subroutine assemble_state_redist_interface_source_3d(tree,geometry,registers,nv,depth,source_halos,ok,boundary_faces)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(state_redist_register_3d), intent(in) :: registers(:)
    integer, intent(in) :: nv,depth
    type(amr_tree_state_3d), intent(inout) :: source_halos
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(amr_tree_state_3d) :: sources,pending
    type(state_redist_register_3d) :: expected
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:),width(:),requests(:,:)
    real(dp), allocatable :: contributions(:,:)
    integer(int64) :: cell(3)
    integer :: p,q,i,j,k,c,n(3),dims(3),extra(3),index(3),level,count
    real(dp) :: delta(nv)
    logical :: good,inside,periodic(3)
    ok=.false.
    if(.not.tree%is_valid().or.nv<1.or.depth<0.or.depth+1>=tree%level_count()) return
    count=size(tree%nodes)
    if(size(geometry)/=count.or.size(registers)/=count) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,good)
    if(.not.good) return
    allocate(sources%fields(count),pending%fields(count),width(count));width=0
    periodic=.false.
    if(present(boundary_faces)) then
      do c=1,3
        if((boundary_faces(2*c-1)=='periodic').neqv.(boundary_faces(2*c)=='periodic')) return
        periodic(c)=boundary_faces(2*c-1)=='periodic'
      enddo
    endif
    do p=1,count
      level=tree%nodes(p)%level
      if(level/=depth.and.level/=depth+1) cycle
      n=tree%nodes(p)%shape;dims=[geometry(p)%nx,geometry(p)%ny,geometry(p)%nz];extra=dims-n
      if(.not.geometry(p)%is_valid().or.any(extra<2).or.any(modulo(extra,2)/=0).or.any(extra/=extra(1))) return
      width(p)=extra(1)/2
      if(.not.registers(p)%is_valid(nv,dims)) return
      call initialize_state_redist_tree_register_3d(tree,p,nv,width(p),expected,good,boundary_faces,registers(p)%scale)
      if(.not.good) return
      if(allocated(expected%coarse_flags).neqv.allocated(registers(p)%coarse_flags)) return
      if(allocated(expected%fine_ghost).neqv.allocated(registers(p)%fine_ghost)) return
      if(.not.allocated(registers(p)%valid_cells)) return
      if(any(expected%valid_cells.neqv.registers(p)%valid_cells)) return
      if(allocated(expected%coarse_flags)) then
        if(any(expected%coarse_flags/=registers(p)%coarse_flags)) return
      endif
      if(allocated(expected%fine_ghost)) then
        if(any(expected%fine_ghost.neqv.registers(p)%fine_ghost)) return
      endif
      if(level/=depth) cycle
      allocate(sources%fields(p)%values(nv,n(1),n(2),n(3)),source=0._dp)
      allocate(pending%fields(p)%values(nv,dims(1),dims(2),dims(3)),source=0._dp)
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        index=[i,j,k]+width(p)
        if(registers(p)%coarse_flags(index(1),index(2),index(3))/=srd_coarse_boundary) cycle
        sources%fields(p)%values(:,i,j,k)=registers(p)%coarse_delta(:,index(1),index(2),index(3))
      enddo
      enddo
      enddo
    enddo
    do p=1,count
      if(tree%nodes(p)%level/=depth+1) cycle
      call restrict_state_redist_patch_source_3d(tree,p,geometry(p),registers(p),nv,requests,contributions, &
        good,boundary_faces)
      if(.not.good) return
      do i=1,size(requests,2)
        q=requests(1,i);index=requests(2:4,i)
        if(registers(q)%coarse_flags(index(1)+width(q),index(2)+width(q),index(3)+width(q))/=srd_coarse_boundary) return
        delta=contributions(:,i)
        call add_register_increment(sources%fields(q)%values(:,index(1),index(2),index(3)),delta,good)
        if(.not.good) return
        sources%fields(q)%values(:,index(1),index(2),index(3))=delta
      enddo
    enddo
    ! Coarse peer/periodic roots need the complete accepted source during
    ! re-reflux. Copy valid owners into valid+1 images, never interpolate it.
    do p=1,count
      if(tree%nodes(p)%level/=depth) cycle
      n=tree%nodes(p)%shape
      do k=width(p),width(p)+n(3)+1
      do j=width(p),width(p)+n(2)+1
      do i=width(p),width(p)+n(1)+1
        cell=origin(:,p)+int([i,j,k]-width(p)-1,int64)
        call canonical(cell,depth,inside)
        if(.not.inside) cycle
        q=donor(cell,depth)
        if(q==0) cycle
        index=int(cell-origin(:,q))+1
        pending%fields(p)%values(:,i,j,k)=sources%fields(q)%values(:,index(1),index(2),index(3))
      enddo
      enddo
      enddo
    enddo
    call move_alloc(pending%fields,source_halos%fields);ok=.true.
  contains
    subroutine canonical(cell,level,inside)
      integer(int64), intent(inout) :: cell(3)
      integer, intent(in) :: level
      logical, intent(out) :: inside
      integer :: a
      inside=.true.
      do a=1,3
        if(periodic(a)) then
          cell(a)=modulo(cell(a),domain(a,level+1))
        else if(cell(a)<0.or.cell(a)>=domain(a,level+1)) then
          inside=.false.
        endif
      enddo
    end subroutine canonical
    integer function donor(cell,level) result(node)
      integer(int64), intent(in) :: cell(3)
      integer, intent(in) :: level
      integer :: q
      node=0
      do q=1,count
        if(tree%nodes(q)%level/=level) cycle
        if(any(cell<origin(:,q)).or.any(cell>=origin(:,q)+int(tree%nodes(q)%shape,int64))) cycle
        node=q;return
      enddo
    end function donor
  end subroutine assemble_state_redist_interface_source_3d

  ! One fine patch's signed FineAdd records. The complete fine block supplies
  ! kappa; only valid+1 fine ghosts supply dm. Serial and MPI use this same
  ! backend, so ownership changes neither the restriction nor its signs.
  subroutine restrict_state_redist_patch_source_3d(tree,p,geometry,reg,nv,requests,contributions,ok,boundary_faces)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: p,nv
    type(eb_geometry_3d), intent(in) :: geometry
    type(state_redist_register_3d), intent(in) :: reg
    integer, allocatable, intent(inout) :: requests(:,:)
    real(dp), allocatable, intent(inout) :: contributions(:,:)
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(state_redist_register_3d) :: expected
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:),pending(:,:)
    real(dp), allocatable :: values(:,:)
    integer(int64) :: lower(3),upper(3),first(3),last(3),cell(3),point(3),block_count,capacity,ic,jc,kc
    integer :: n(3),dims(3),extra(3),width,depth,ratio,c,ii,jj,kk,lo(3),hi(3),q,index(3),count
    real(dp) :: total(nv),delta(nv),fluid,cutoff
    logical :: good,inside,periodic(3)
    ok=.false.
    if(.not.tree%is_valid().or.nv<1) return
    if(p<1.or.p>size(tree%nodes)) return
    depth=tree%nodes(p)%level-1
    if(depth<0) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,good)
    if(.not.good) return
    n=tree%nodes(p)%shape;dims=[geometry%nx,geometry%ny,geometry%nz];extra=dims-n
    if(.not.geometry%is_valid().or.any(extra<2).or.any(modulo(extra,2)/=0).or.any(extra/=extra(1))) return
    width=extra(1)/2
    if(.not.reg%is_valid(nv,dims)) return
    call initialize_state_redist_tree_register_3d(tree,p,nv,width,expected,good,boundary_faces,reg%scale)
    if(.not.good) return
    if(.not.allocated(reg%fine_ghost).or..not.allocated(reg%valid_cells)) return
    if(any(expected%fine_ghost.neqv.reg%fine_ghost).or.any(expected%valid_cells.neqv.reg%valid_cells)) return
    periodic=.false.
    if(present(boundary_faces)) then
      do c=1,3
        if((boundary_faces(2*c-1)=='periodic').neqv.(boundary_faces(2*c)=='periodic')) return
        periodic(c)=boundary_faces(2*c-1)=='periodic'
      enddo
    endif
    ratio=ratios(depth+2);block_count=1
    do c=1,3
      if(int(ratio,int64)>huge(block_count)/block_count) return
      block_count=block_count*int(ratio,int64)
    enddo
    cutoff=1.e-14_dp*real(block_count,dp)
    first=origin(:,p)-1_int64;last=origin(:,p)+int(n,int64)
    lower=(first-modulo(first,int(ratio,int64)))/ratio
    upper=(last-modulo(last,int(ratio,int64)))/ratio
    capacity=1
    do c=1,3
      if(upper(c)-lower(c)+1>int(huge(1)/max(4,nv),int64)/capacity) return
      capacity=capacity*(upper(c)-lower(c)+1)
    enddo
    allocate(pending(4,int(capacity)),values(nv,int(capacity)));count=0
    do kc=lower(3),upper(3)
    do jc=lower(2),upper(2)
    do ic=lower(1),upper(1)
      cell=[ic,jc,kc]
      if(any(cell>huge(1_int64)/int(ratio,int64))) return
      point=cell*int(ratio,int64)-origin(:,p)+int(width+1,int64)
      if(any(point<1).or.any(point+ratio-1>int(n+2*width,int64))) return
      lo=int(point);hi=lo+ratio-1;fluid=0;total=0
      do kk=lo(3),hi(3)
      do jj=lo(2),hi(2)
      do ii=lo(1),hi(1)
        fluid=fluid+geometry%volume_fraction(ii,jj,kk)
        if(any([ii,jj,kk]<width).or.any([ii,jj,kk]>width+n+1)) cycle
        if(.not.reg%fine_ghost(ii,jj,kk)) cycle
        delta=reg%fine_delta(:,ii,jj,kk)
        call add_register_increment(total,delta,good)
        if(.not.good) return
        total=delta
      enddo
      enddo
      enddo
      if(fluid<=cutoff.or.all(total==0)) cycle
      if(fluid<1) then
        if(any(abs(total)>huge(1._dp)*fluid)) return
      endif
      delta=total/fluid
      inside=.true.
      do c=1,3
        if(periodic(c)) then
          cell(c)=modulo(cell(c),domain(c,depth+1))
        else if(cell(c)<0.or.cell(c)>=domain(c,depth+1)) then
          inside=.false.
        endif
      enddo
      if(.not.inside) cycle
      q=0
      do c=1,size(tree%nodes)
        if(tree%nodes(c)%level/=depth) cycle
        if(any(cell<origin(:,c)).or.any(cell>=origin(:,c)+int(tree%nodes(c)%shape,int64))) cycle
        q=c;exit
      enddo
      if(q==0) return
      index=int(cell-origin(:,q))+1
      count=count+1;pending(:,count)=[q,index];values(:,count)=delta
    enddo
    enddo
    enddo
    requests=pending(:,:count);contributions=values(:,:count);ok=.true.
  end subroutine restrict_state_redist_patch_source_3d

  ! A single donor-to-recipient contribution after alpha/nrs weighting.
  ! Coarse delta is a density increment; fine delta includes recipient kappa.
  ! These signs and the coarse volume ratio follow MLStateRedistribute.
  subroutine accumulate_state_redist_transfer_3d(reg,donor,recipient,contribution,donor_kappa,recipient_kappa,ok)
    type(state_redist_register_3d), intent(inout) :: reg
    integer, intent(in) :: donor(3),recipient(3)
    real(dp), intent(in) :: contribution(:),donor_kappa,recipient_kappa
    logical, intent(out) :: ok
    integer :: dims(3)
    ok=.false.
    if(allocated(reg%coarse_flags)) then
      dims=shape(reg%coarse_flags)
    else if(allocated(reg%fine_ghost)) then
      dims=shape(reg%fine_ghost)
    else
      return
    endif
    if(.not.reg%is_valid(size(contribution),dims)) return
    call record_state_redist_transfer_3d(reg,donor,recipient,contribution,donor_kappa,recipient_kappa,ok)
  end subroutine accumulate_state_redist_transfer_3d

  subroutine record_state_redist_transfer_3d(reg,donor,recipient,contribution,donor_kappa,recipient_kappa,ok)
    type(state_redist_register_3d), intent(inout) :: reg
    integer, intent(in) :: donor(3),recipient(3)
    real(dp), intent(in) :: contribution(:),donor_kappa,recipient_kappa
    logical, intent(out) :: ok
    real(dp) :: cd(size(contribution)),cr(size(contribution)),fd(size(contribution)),fr(size(contribution))
    real(dp) :: scaled(size(contribution)),ratio
    logical :: valid,donor_valid,recipient_valid
    integer :: dims(3),df,rf
    ok=.false.
    if(allocated(reg%coarse_flags)) then
      dims=shape(reg%coarse_flags)
    else if(allocated(reg%fine_ghost)) then
      dims=shape(reg%fine_ghost)
    else
      return
    endif
    if(any(donor<1).or.any(donor>dims).or.any(recipient<1).or.any(recipient>dims)) return
    if(any(.not.ieee_is_finite(contribution)).or..not.all(ieee_is_finite([donor_kappa,recipient_kappa]))) return
    if(min(donor_kappa,recipient_kappa)<=0.or.max(donor_kappa,recipient_kappa)>1) return
    if(reg%scale>1) then
      if(any(abs(contribution)>huge(1._dp)/reg%scale)) return
    endif
    scaled=reg%scale*contribution
    donor_valid=.true.;recipient_valid=.true.
    if(allocated(reg%valid_cells)) then
      donor_valid=reg%valid_cells(donor(1),donor(2),donor(3))
      recipient_valid=reg%valid_cells(recipient(1),recipient(2),recipient(3))
    endif
    cd=0;cr=0;fd=0;fr=0
    if(allocated(reg%coarse_flags)) then
      df=reg%coarse_flags(donor(1),donor(2),donor(3));rf=reg%coarse_flags(recipient(1),recipient(2),recipient(3))
      if(recipient_valid.and.rf==srd_coarse_boundary.and.df==srd_coarse_fine_covered) cr=-scaled
      if(donor_valid.and.rf==srd_coarse_fine_covered.and.df==srd_coarse_boundary) then
        if(recipient_kappa>huge(1._dp)*donor_kappa) return
        ratio=recipient_kappa/donor_kappa
        if(ratio>1) then
          if(any(abs(scaled)>huge(1._dp)/ratio)) return
        endif
        cd=scaled*ratio
      endif
      call add_register_increment(reg%coarse_delta(:,donor(1),donor(2),donor(3)),cd,valid)
      if(.not.valid) return
      call add_register_increment(reg%coarse_delta(:,recipient(1),recipient(2),recipient(3)),cr,valid)
      if(.not.valid) return
    endif
    if(allocated(reg%fine_ghost)) then
      if(reg%fine_ghost(donor(1),donor(2),donor(3)).and. &
        recipient_valid.and..not.reg%fine_ghost(recipient(1),recipient(2),recipient(3))) fd=-scaled*recipient_kappa
      if(donor_valid.and..not.reg%fine_ghost(donor(1),donor(2),donor(3)).and. &
        reg%fine_ghost(recipient(1),recipient(2),recipient(3))) fr=scaled*recipient_kappa
      call add_register_increment(reg%fine_delta(:,donor(1),donor(2),donor(3)),fd,valid)
      if(.not.valid) return
      call add_register_increment(reg%fine_delta(:,recipient(1),recipient(2),recipient(3)),fr,valid)
      if(.not.valid) return
    endif
    if(any(.not.ieee_is_finite(cd)).or.any(.not.ieee_is_finite(cr)).or. &
      any(.not.ieee_is_finite(fd)).or.any(.not.ieee_is_finite(fr))) return
    if(allocated(reg%coarse_flags)) then
      reg%coarse_delta(:,donor(1),donor(2),donor(3))=cd
      reg%coarse_delta(:,recipient(1),recipient(2),recipient(3))=cr
    endif
    if(allocated(reg%fine_ghost)) then
      reg%fine_delta(:,donor(1),donor(2),donor(3))=fd
      reg%fine_delta(:,recipient(1),recipient(2),recipient(3))=fr
    endif
    ok=.true.
  end subroutine record_state_redist_transfer_3d

  subroutine add_register_increment(old,increment,ok)
    real(dp), intent(in) :: old(:)
    real(dp), intent(inout) :: increment(:)
    logical, intent(out) :: ok
    integer :: c
    ok=.false.
    do c=1,size(old)
      if(increment(c)>0) then
        if(old(c)>huge(1._dp)-increment(c)) return
      else if(increment(c)<0) then
        if(old(c)<-huge(1._dp)-increment(c)) return
      endif
    enddo
    increment=old+increment;ok=.true.
  end subroutine add_register_increment

  subroutine build_state_redist_coarse_flags_3d(refined,flags,ok)
    logical, intent(in) :: refined(:,:,:)
    integer, allocatable, intent(inout) :: flags(:,:,:)
    logical, intent(out) :: ok
    integer, allocatable :: candidate(:,:,:)
    integer :: dims(3),i,j,k,lo(3),hi(3)
    ok=.false.;dims=shape(refined)
    if(any(dims<1)) return
    allocate(candidate(dims(1),dims(2),dims(3)));candidate=srd_coarse_regular
    do k=1,dims(3)
    do j=1,dims(2)
    do i=1,dims(1)
      if(refined(i,j,k)) then
        candidate(i,j,k)=srd_coarse_fine_covered
      else
        lo=max(1,[i,j,k]-1);hi=min(dims,[i,j,k]+1)
        if(any(refined(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))) candidate(i,j,k)=srd_coarse_boundary
      endif
    enddo
    enddo
    enddo
    call move_alloc(candidate,flags);ok=.true.
  end subroutine build_state_redist_coarse_flags_3d

  ! Match EBFluxRegister FineAdd(dm): sum fine kappa-weighted corrections,
  ! then divide by the complete coarse block's sum of fine volume fractions.
  ! Origins are zero-based global cell coordinates, including negative halos.
  ! active is the coarse/fine register mask and excludes same-level peers.
  subroutine restrict_state_redist_fine_delta_3d(fine,volume_fraction,ratio,coarse_origin,fine_origin, &
      coarse_delta,active,ok,threshold)
    type(state_redist_register_3d), intent(in) :: fine
    real(dp), intent(in) :: volume_fraction(:,:,:)
    integer, intent(in) :: ratio(3),coarse_origin(3),fine_origin(3)
    real(dp), intent(inout) :: coarse_delta(:,:,:,:)
    logical, intent(in) :: active(:,:,:)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: threshold
    real(dp), allocatable :: candidate(:,:,:,:)
    real(dp) :: cutoff,fluid,total(size(coarse_delta,1)),increment(size(coarse_delta,1))
    integer(int64) :: lower64(3),upper64(3),block_count
    integer :: dims(3),fdims(3),lo(3),hi(3),i,j,k,ii,jj,kk,nv,c
    logical :: valid
    ok=.false.;nv=size(coarse_delta,1);dims=shape(active);fdims=shape(volume_fraction)
    if(nv<1.or.any(dims<1).or.any(fdims<1).or.any(ratio<1)) return
    if(any(shape(coarse_delta)/=[nv,dims])) return
    if(.not.fine%is_valid(nv,fdims).or..not.allocated(fine%fine_delta)) return
    if(any(.not.ieee_is_finite(coarse_delta)).or.any(.not.ieee_is_finite(volume_fraction))) return
    if(any(volume_fraction<0).or.any(volume_fraction>1)) return
    cutoff=1.e-14_dp
    if(present(threshold)) cutoff=threshold
    if(.not.ieee_is_finite(cutoff).or.cutoff<0.or.cutoff>1) return
    block_count=1
    do c=1,3
      if(int(ratio(c),int64)>huge(block_count)/block_count) return
      block_count=block_count*int(ratio(c),int64)
    enddo
    cutoff=cutoff*real(block_count,dp)
    candidate=coarse_delta
    do k=1,dims(3)
    do j=1,dims(2)
    do i=1,dims(1)
      if(.not.active(i,j,k)) cycle
      lower64=(int(coarse_origin,int64)+int([i,j,k]-1,int64))*int(ratio,int64)-int(fine_origin,int64)+1
      upper64=lower64+int(ratio,int64)-1
      if(any(lower64<1).or.any(upper64>int(fdims,int64))) return
      lo=int(lower64);hi=int(upper64);fluid=0;total=0
      do kk=lo(3),hi(3)
      do jj=lo(2),hi(2)
      do ii=lo(1),hi(1)
        fluid=fluid+volume_fraction(ii,jj,kk)
        increment=fine%fine_delta(:,ii,jj,kk)
        call add_register_increment(total,increment,valid)
        if(.not.valid) return
        total=increment
      enddo
      enddo
      enddo
      ! The fixed reference uses a strict threshold comparison.
      if(fluid<=cutoff) cycle
      if(fluid<1) then
        if(any(abs(total)>huge(1._dp)*fluid)) return
      endif
      increment=total/fluid
      call add_register_increment(candidate(:,i,j,k),increment,valid)
      if(.not.valid) return
      candidate(:,i,j,k)=increment
    enddo
    enddo
    enddo
    coarse_delta=candidate;ok=.true.
  end subroutine restrict_state_redist_fine_delta_3d

  subroutine register_state_neighborhood(g,u,cells,input_weight,output_weight,delta,slope,reg,ok)
    type(eb_geometry_3d), intent(in) :: g
    real(dp), intent(in) :: u(:,:,:,:),input_weight(:),output_weight(:),delta(:,:),slope(:,:)
    integer, intent(in) :: cells(:,:)
    type(state_redist_register_3d), intent(inout) :: reg
    logical, intent(out) :: ok
    real(dp) :: contribution(size(u,1))
    integer :: d,r,di(3),ri(3)
    ok=.false.
    do r=1,size(cells,2)
      ri=cells(:,r)
      if(g%volume_fraction(ri(1),ri(2),ri(3))<=0) cycle
      do d=1,size(cells,2)
        di=cells(:,d)
        if(g%volume_fraction(di(1),di(2),di(3))<=0) cycle
        contribution=output_weight(r)*input_weight(d)*(u(:,di(1),di(2),di(3))+matmul(delta(:,r),slope))
        call record_state_redist_transfer_3d(reg,di,ri,contribution, &
          g%volume_fraction(di(1),di(2),di(3)),g%volume_fraction(ri(1),ri(2),ri(3)),ok)
        if(.not.ok) return
      enddo
    enddo
    ok=.true.
  end subroutine register_state_neighborhood

  ! Re-express redistribution as a source rate for SDC stage assembly.
  ! State redistribution acts on U + dt*R/scale, then returns
  ! scale*(redistributed-U)/dt only where the neighborhood changed.
  subroutine redistribute_sdc_eb_rate_3d(geometry,state,conservative_rhs,dt,method,rate,ok,target_volume_fraction,max_order, &
      boundary_state,boundary_geometry,boundary_faces,registers,box,update_scale)
    real(dp), optional, intent(in) :: update_scale(:,:,:)
    real(dp) :: scale_value
    type(state_redist_register_3d), optional, intent(inout) :: registers
    type(state_redist_box_3d), optional, intent(in) :: box
    type(state_redist_register_3d), allocatable :: pending_registers
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: state(:,:,:,:),conservative_rhs(:,:,:,:),dt
    character(len=*), intent(in) :: method
    real(dp), intent(inout) :: rate(:,:,:,:)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: target_volume_fraction
    integer, optional, intent(in) :: max_order
    real(dp), optional, intent(in) :: boundary_state(:,:,:,:)
    type(eb_geometry_3d), optional, intent(in) :: boundary_geometry
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), allocatable :: old(:,:,:,:),rhs(:,:,:,:),provisional(:,:,:,:),candidate(:,:,:,:)
    integer :: n(3),nc,i,j,k
    logical :: valid
    logical, allocatable :: affected(:,:,:)
    ok=.false.
    if(.not.geometry%is_valid().or..not.ieee_is_finite(dt)) return
    if(dt<=0) return
    if(present(box)) then
      if(method/='state_redist') return
    endif
    n=[geometry%nx,geometry%ny,geometry%nz];nc=size(state,1)
    if(nc<1.or.any(shape(state)/=[nc,n])) return
    if(present(box)) then
      if(.not.box%is_valid(n)) return
    endif
    if(any(shape(conservative_rhs)/=shape(state)).or.any(shape(rate)/=shape(state))) return
    if(present(update_scale)) then
      if(method/='state_redist'.or.any(shape(update_scale)/=n)) return
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(present(box)) then
          if(.not.box%contains_valid([i,j,k],3)) cycle
        endif
        if(geometry%cell_type(i,j,k)==eb_covered_cell_3d) cycle
        if(.not.ieee_is_finite(update_scale(i,j,k))) return
        if(update_scale(i,j,k)<=0) return
      enddo
      enddo
      enddo
    endif
    if(present(registers)) allocate(pending_registers,source=registers)
    old=state;rhs=conservative_rhs
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      if(present(box)) then
        ! Match ApplyMLRedistribution's grow-1 exterior reset. The second
        ! and third ghost layers retain the caller's supplied forcing.
        if(.not.box%contains_valid([i,j,k],3)) rhs(:,i,j,k)=0
        if(box%contains_valid([i,j,k],1).and..not.box%contains_domain([i,j,k],1)) rhs(:,i,j,k)=0
      endif
      if(geometry%cell_type(i,j,k)/=eb_covered_cell_3d) cycle
      old(:,i,j,k)=0;rhs(:,i,j,k)=0
    end do
    end do
    end do
    if(.not.all(ieee_is_finite(old)).or..not.all(ieee_is_finite(rhs))) return
    allocate(candidate,mold=state)
    select case(method)
    case('none')
      candidate=rhs
    case('flux_redist')
      call reactive_eb_flux_redistribute_3d(geometry,rhs,candidate,valid)
      if(.not.valid) return
    case('state_redist')
      if(present(update_scale)) then
        provisional=old
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(present(box)) then
            if(.not.box%contains_valid([i,j,k],3)) cycle
          endif
          if(geometry%cell_type(i,j,k)==eb_covered_cell_3d) cycle
          scale_value=update_scale(i,j,k)
          call scaled_sdc_state_3d(old(:,i,j,k),rhs(:,i,j,k),dt,scale_value,provisional(:,i,j,k),valid)
          if(.not.valid) return
        enddo
        enddo
        enddo
      else
        provisional=old+dt*rhs
      endif
      if(.not.all(ieee_is_finite(provisional))) return
      allocate(affected(n(1),n(2),n(3)))
      call reactive_eb_weighted_state_redistribute_3d(geometry,provisional,candidate,valid,target_volume_fraction,max_order, &
        boundary_state,boundary_geometry,boundary_faces,pending_registers,box,affected)
      if(.not.valid) return
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(affected(i,j,k)) then
          if(present(update_scale)) then
            scale_value=update_scale(i,j,k)
            call scaled_sdc_rate_3d(candidate(:,i,j,k),old(:,i,j,k),dt,scale_value,rhs(:,i,j,k),valid)
            if(.not.valid) return
            candidate(:,i,j,k)=rhs(:,i,j,k)
          else
            candidate(:,i,j,k)=(candidate(:,i,j,k)-old(:,i,j,k))/dt
          endif
        else
          candidate(:,i,j,k)=rhs(:,i,j,k)
        endif
      enddo
      enddo
      enddo
    case default
      return
    end select
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      if(geometry%cell_type(i,j,k)==eb_covered_cell_3d) candidate(:,i,j,k)=0
    end do
    end do
    end do
    if(.not.all(ieee_is_finite(candidate))) return
    rate=candidate
    if(present(registers)) registers=pending_registers
    ok=.true.
  end subroutine redistribute_sdc_eb_rate_3d


  ! Keep failure transactional even with overflow trapping enabled. These
  ! guards preserve the reference operation order, including dt*R/scale.
  subroutine scaled_sdc_state_3d(old,rhs,dt,scale,state,ok)
    real(dp), intent(in) :: old(:),rhs(:),dt,scale
    real(dp), intent(out) :: state(:)
    logical, intent(out) :: ok
    real(dp) :: increment(size(old))
    integer :: c
    ok=.false.
    if(dt>1._dp) then
      if(any(abs(rhs)>huge(1._dp)/dt)) return
    endif
    increment=dt*rhs
    if(scale<1._dp) then
      if(any(abs(increment)>huge(1._dp)*scale)) return
    endif
    increment=increment/scale
    do c=1,size(old)
      if(old(c)>0.and.increment(c)>0) then
        if(increment(c)>huge(1._dp)-old(c)) return
      else if(old(c)<0.and.increment(c)<0) then
        if(increment(c)<-huge(1._dp)-old(c)) return
      endif
    enddo
    state=old+increment
    ok=all(ieee_is_finite(state))
  end subroutine scaled_sdc_state_3d

  subroutine scaled_sdc_rate_3d(state,old,dt,scale,rate,ok)
    real(dp), intent(in) :: state(:),old(:),dt,scale
    real(dp), intent(out) :: rate(:)
    logical, intent(out) :: ok
    real(dp) :: difference(size(old))
    integer :: c
    ok=.false.
    do c=1,size(old)
      if(old(c)<0.and.state(c)>0) then
        if(state(c)>huge(1._dp)+old(c)) return
      else if(old(c)>0.and.state(c)<0) then
        if(state(c)<-huge(1._dp)+old(c)) return
      endif
    enddo
    difference=state-old
    if(scale>1._dp) then
      if(any(abs(difference)>huge(1._dp)/scale)) return
    endif
    difference=scale*difference
    if(dt<1._dp) then
      if(any(abs(difference)>huge(1._dp)*dt)) return
    endif
    rate=difference/dt
    ok=all(ieee_is_finite(rate))
  end subroutine scaled_sdc_rate_3d

  subroutine reactive_eb_flux_redistribute_3d( &
      geometry, conservative_rhs, redistributed_rhs, ok)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: conservative_rhs(:, :, :, :)
    real(dp), intent(out) :: redistributed_rhs(:, :, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate(:, :, :, :)
    real(dp), allocatable :: neighborhood_rhs(:, :, :, :)
    real(dp), allocatable :: neighbor_volume_fraction(:, :, :)
    real(dp), allocatable :: excess(:)
    real(dp) :: kappa, total_volume_fraction
    integer :: i, j, k, ncomp

    redistributed_rhs = 0.0_dp
    ok = .false.
    if (.not. geometry%is_valid()) return
    ncomp = size(conservative_rhs, 1)
    if (ncomp < 1 .or. &
        size(conservative_rhs, 2) /= geometry%nx .or. &
        size(conservative_rhs, 3) /= geometry%ny .or. &
        size(conservative_rhs, 4) /= geometry%nz .or. &
        any(shape(redistributed_rhs) /= shape(conservative_rhs)) .or. &
        any(.not. ieee_is_finite(conservative_rhs))) return

    allocate(candidate(ncomp, geometry%nx, geometry%ny, geometry%nz))
    allocate(neighborhood_rhs( &
      ncomp, geometry%nx, geometry%ny, geometry%nz))
    allocate(neighbor_volume_fraction( &
      geometry%nx, geometry%ny, geometry%nz))
    allocate(excess(ncomp))
    candidate = 0.0_dp
    neighborhood_rhs = 0.0_dp
    neighbor_volume_fraction = 0.0_dp

    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          select case (geometry%cell_type(i, j, k))
          case (eb_covered_cell_3d)
            cycle
          case (eb_regular_cell_3d)
            candidate(:, i, j, k) = conservative_rhs(:, i, j, k)
          case (eb_cut_cell_3d)
            kappa = geometry%volume_fraction(i, j, k)
            total_volume_fraction = kappa
            neighborhood_rhs(:, i, j, k) = &
              kappa * conservative_rhs(:, i, j, k)
            call accumulate_neighbors( &
              geometry, conservative_rhs, i, j, k, &
              neighborhood_rhs(:, i, j, k), total_volume_fraction, &
              neighbor_volume_fraction(i, j, k))
            if (neighbor_volume_fraction(i, j, k) <= 0.0_dp .or. &
                total_volume_fraction <= kappa) return
            neighborhood_rhs(:, i, j, k) = &
              neighborhood_rhs(:, i, j, k) / total_volume_fraction
            candidate(:, i, j, k) = &
              kappa * conservative_rhs(:, i, j, k) + &
              (1.0_dp - kappa) * neighborhood_rhs(:, i, j, k)
          case default
            return
          end select
        end do
      end do
    end do

    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) /= eb_cut_cell_3d) cycle
          kappa = geometry%volume_fraction(i, j, k)
          excess = kappa * (1.0_dp - kappa) * &
            (conservative_rhs(:, i, j, k) - &
              neighborhood_rhs(:, i, j, k))
          call distribute_to_neighbors( &
            geometry, i, j, k, &
            excess / neighbor_volume_fraction(i, j, k), candidate)
        end do
      end do
    end do
    if (any(.not. ieee_is_finite(candidate))) return

    redistributed_rhs = candidate
    ok = .true.
  end subroutine reactive_eb_flux_redistribute_3d

  subroutine accumulate_neighbors( &
      geometry, rhs, i, j, k, weighted_rhs, total_volume, &
      neighbor_volume)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: rhs(:, :, :, :)
    integer, intent(in) :: i, j, k
    real(dp), intent(inout) :: weighted_rhs(:), total_volume
    real(dp), intent(out) :: neighbor_volume

    neighbor_volume = 0.0_dp
    if (i > 1) then
      if (geometry%x_face_fraction(i - 1, j, k) > 0.0_dp) &
        call accumulate_cell(i - 1, j, k)
    end if
    if (i < geometry%nx) then
      if (geometry%x_face_fraction(i, j, k) > 0.0_dp) &
        call accumulate_cell(i + 1, j, k)
    end if
    if (j > 1) then
      if (geometry%y_face_fraction(i, j - 1, k) > 0.0_dp) &
        call accumulate_cell(i, j - 1, k)
    end if
    if (j < geometry%ny) then
      if (geometry%y_face_fraction(i, j, k) > 0.0_dp) &
        call accumulate_cell(i, j + 1, k)
    end if
    if (k > 1) then
      if (geometry%z_face_fraction(i, j, k - 1) > 0.0_dp) &
        call accumulate_cell(i, j, k - 1)
    end if
    if (k < geometry%nz) then
      if (geometry%z_face_fraction(i, j, k) > 0.0_dp) &
        call accumulate_cell(i, j, k + 1)
    end if

  contains

    subroutine accumulate_cell(neighbor_i, neighbor_j, neighbor_k)
      integer, intent(in) :: neighbor_i, neighbor_j, neighbor_k
      real(dp) :: neighbor_kappa

      if (geometry%cell_type(neighbor_i, neighbor_j, neighbor_k) == &
          eb_covered_cell_3d) return
      neighbor_kappa = &
        geometry%volume_fraction(neighbor_i, neighbor_j, neighbor_k)
      weighted_rhs = weighted_rhs + neighbor_kappa * &
        rhs(:, neighbor_i, neighbor_j, neighbor_k)
      total_volume = total_volume + neighbor_kappa
      neighbor_volume = neighbor_volume + neighbor_kappa
    end subroutine accumulate_cell

  end subroutine accumulate_neighbors

  subroutine distribute_to_neighbors( &
      geometry, i, j, k, increment, redistributed_rhs)
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: i, j, k
    real(dp), intent(in) :: increment(:)
    real(dp), intent(inout) :: redistributed_rhs(:, :, :, :)

    if (i > 1) then
      if (geometry%x_face_fraction(i - 1, j, k) > 0.0_dp) &
        call add_to_cell(i - 1, j, k)
    end if
    if (i < geometry%nx) then
      if (geometry%x_face_fraction(i, j, k) > 0.0_dp) &
        call add_to_cell(i + 1, j, k)
    end if
    if (j > 1) then
      if (geometry%y_face_fraction(i, j - 1, k) > 0.0_dp) &
        call add_to_cell(i, j - 1, k)
    end if
    if (j < geometry%ny) then
      if (geometry%y_face_fraction(i, j, k) > 0.0_dp) &
        call add_to_cell(i, j + 1, k)
    end if
    if (k > 1) then
      if (geometry%z_face_fraction(i, j, k - 1) > 0.0_dp) &
        call add_to_cell(i, j, k - 1)
    end if
    if (k < geometry%nz) then
      if (geometry%z_face_fraction(i, j, k) > 0.0_dp) &
        call add_to_cell(i, j, k + 1)
    end if

  contains

    subroutine add_to_cell(neighbor_i, neighbor_j, neighbor_k)
      integer, intent(in) :: neighbor_i, neighbor_j, neighbor_k

      if (geometry%cell_type(neighbor_i, neighbor_j, neighbor_k) == &
          eb_covered_cell_3d) return
      redistributed_rhs(:, neighbor_i, neighbor_j, neighbor_k) = &
        redistributed_rhs(:, neighbor_i, neighbor_j, neighbor_k) + &
        increment
    end subroutine add_to_cell

  end subroutine distribute_to_neighbors

  subroutine advance_reactive_eb_redistributed_3d( &
      species, state, temperature, geometry, conservative_rhs, dt, &
      new_state, new_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: conservative_rhs(:, :, :, :), dt
    real(dp), intent(out) :: new_state(:, :, :, :)
    real(dp), intent(out) :: new_temperature(:, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    real(dp), allocatable :: redistributed_rhs(:, :, :, :)
    real(dp), allocatable :: primitive(:)
    real(dp) :: recovered_temperature, sound_speed
    logical :: local_ok
    integer :: i, j, k, nvar

    new_state = 0.0_dp
    new_temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        size(state, 4) /= geometry%nz .or. &
        any(shape(temperature) /= &
          [geometry%nx, geometry%ny, geometry%nz]) .or. &
        any(shape(conservative_rhs) /= shape(state)) .or. &
        any(shape(new_state) /= shape(state)) .or. &
        any(shape(new_temperature) /= shape(temperature))) return
    new_state = state
    new_temperature = temperature
    if (.not. ieee_is_finite(dt) .or. dt < 0.0_dp .or. &
        any(.not. ieee_is_finite(state)) .or. &
        any(.not. ieee_is_finite(temperature)) .or. &
        any(.not. ieee_is_finite(conservative_rhs))) return

    allocate(redistributed_rhs(nvar, geometry%nx, geometry%ny, geometry%nz))
    call reactive_eb_flux_redistribute_3d( &
      geometry, conservative_rhs, redistributed_rhs, local_ok)
    if (.not. local_ok) return

    allocate(candidate_state, source=state)
    allocate(candidate_temperature, source=temperature)
    allocate(primitive(reactive_nprim(size(species))))
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          if (temperature(i, j, k) <= 0.0_dp) return
          candidate_state(:, i, j, k) = state(:, i, j, k) + &
            dt * redistributed_rhs(:, i, j, k)
          call reactive_conserved_to_primitive( &
            species, candidate_state(:, i, j, k), temperature(i, j, k), &
            primitive, recovered_temperature, sound_speed, local_ok)
          if (.not. local_ok) return
          candidate_temperature(i, j, k) = recovered_temperature
        end do
      end do
    end do
    if (any(.not. ieee_is_finite(candidate_state)) .or. &
        any(.not. ieee_is_finite(candidate_temperature))) return

    new_state = candidate_state
    new_temperature = candidate_temperature
    ok = .true.
  end subroutine advance_reactive_eb_redistributed_3d

  subroutine reactive_eb_weighted_state_redistribute_3d( &
      geometry, provisional_state, redistributed_state, ok, &
      target_volume_fraction,max_order, &
      boundary_state,boundary_geometry,boundary_faces,registers,box,affected)
    type(state_redist_register_3d), optional, intent(inout) :: registers
    type(state_redist_box_3d), optional, intent(in) :: box
    logical, optional, intent(inout) :: affected(:,:,:)
    type(state_redist_register_3d) :: pending_registers
    integer, allocatable :: register_cells(:,:)
    real(dp), allocatable :: register_in(:),register_out(:),register_delta(:,:),register_slope(:,:)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: provisional_state(:, :, :, :)
    real(dp), intent(out) :: redistributed_state(:, :, :, :)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: boundary_state(:,:,:,:)
    type(eb_geometry_3d), optional, intent(in) :: boundary_geometry
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), intent(in), optional :: target_volume_fraction

    integer, optional, intent(in) :: max_order
    real(dp), allocatable :: centroid(:,:,:,:),slopes(:,:,:,:,:),lower_bound(:),upper_bound(:)
    real(dp) :: position(3,2),delta(3),self_weight,other_weight
    integer :: selected_order,recipient_count,allocation_status,dims(3),image_offset(3)
    integer, allocatable :: neighbor_offset(:,:,:,:)
    logical :: periodic(3)
    character(len=16) :: patch_faces(6)

    logical, allocatable :: has_neighbor(:, :, :)
    integer, allocatable :: neighbor_i(:, :, :)
    integer, allocatable :: neighbor_j(:, :, :)
    integer, allocatable :: neighbor_k(:, :, :)
    integer, allocatable :: neighborhood_count(:, :, :)
    real(dp), allocatable :: alpha_self(:, :, :)
    real(dp), allocatable :: alpha_neighbor(:, :, :)
    real(dp), allocatable :: neighborhood_volume(:, :, :)
    real(dp), allocatable :: neighborhood_state(:, :, :, :)
    real(dp), allocatable :: candidate(:, :, :, :)
    real(dp), parameter :: weight_tolerance = &
      1024.0_dp * epsilon(1.0_dp)
    real(dp) :: target, kappa, receiving_kappa
    logical :: local_ok
    integer :: i, j, k, ni, nj, nk, ncomp

    redistributed_state = 0.0_dp
    ok = .false.
    if (.not. geometry%is_valid()) return
    ncomp = size(provisional_state, 1)
    if (ncomp < 1 .or. &
        size(provisional_state, 2) /= geometry%nx .or. &
        size(provisional_state, 3) /= geometry%ny .or. &
        size(provisional_state, 4) /= geometry%nz .or. &
        any(shape(redistributed_state) /= shape(provisional_state)) .or. &
        any(.not. ieee_is_finite(provisional_state))) return

    selected_order=0
    if(present(max_order)) selected_order=max_order
    if(selected_order/=0.and.selected_order/=2.and.selected_order/=4) return
    target = 0.5_dp
    if (present(target_volume_fraction)) target = target_volume_fraction
    if (.not. ieee_is_finite(target) .or. &
        target <= 0.0_dp .or. target > 1.0_dp) return

    periodic=.false.;dims=[geometry%nx,geometry%ny,geometry%nz]
    if(present(affected)) then
      if(any(shape(affected)/=dims)) return
    endif
    if(present(box)) then
      if(.not.box%is_valid(dims)) return
      if(present(boundary_state).or.present(boundary_geometry)) return
      patch_faces='outflow'
      do i=1,3
        if(box%periodic(i)) patch_faces(2*i-1:2*i)='periodic'
      enddo
      if(present(boundary_faces)) patch_faces=boundary_faces
      if(any(patch_faces/='outflow'.and.patch_faces/='inflow'.and.patch_faces/='periodic'.and. &
        patch_faces/='slip_wall'.and.patch_faces/='no_slip_wall'.and.patch_faces/='isothermal_wall')) return
      do i=1,3
        if((patch_faces(2*i-1)=='periodic').neqv.box%periodic(i)) return
        if((patch_faces(2*i)=='periodic').neqv.box%periodic(i)) return
      enddo
    else if(present(boundary_faces)) then
      if(.not.periodic_eb_faces_match(geometry,boundary_faces)) return
      periodic=boundary_faces(1:6:2)=='periodic'
      if(any(periodic)) then
        if(.not.present(boundary_state).or..not.present(boundary_geometry)) return
        if(any([boundary_geometry%nx,boundary_geometry%ny,boundary_geometry%nz]/=dims+4)) return
        if(any(shape(boundary_state)/=[ncomp,dims+4])) return
        if(.not.periodic_eb_faces_match(boundary_geometry,boundary_faces,halo=2)) return
      endif
    endif
    if(present(box)) then
      call redistribute_connected_neighborhoods(geometry,provisional_state,target,redistributed_state,ok,selected_order, &
        boundary_state,boundary_geometry,boundary_faces,registers,box,affected)
      return
    endif
    if(present(registers)) then
      if(.not.registers%is_valid(ncomp,dims)) return
      pending_registers=registers
      if(present(box)) then
        call prepare_state_redist_box_register_3d(box,pending_registers,ncomp,dims,local_ok)
        if(.not.local_ok) return
      endif
      allocate(register_cells(3,2),register_in(2),register_out(2),register_delta(3,2),register_slope(3,ncomp))
    endif
    allocate(neighbor_offset(3,geometry%nx,geometry%ny,geometry%nz));neighbor_offset=0
    allocate(has_neighbor(geometry%nx, geometry%ny, geometry%nz))
    allocate(neighbor_i(geometry%nx, geometry%ny, geometry%nz))
    allocate(neighbor_j(geometry%nx, geometry%ny, geometry%nz))
    allocate(neighbor_k(geometry%nx, geometry%ny, geometry%nz))
    allocate(neighborhood_count(geometry%nx, geometry%ny, geometry%nz))
    allocate(alpha_self(geometry%nx, geometry%ny, geometry%nz))
    allocate(alpha_neighbor(geometry%nx, geometry%ny, geometry%nz))
    allocate(neighborhood_volume(geometry%nx, geometry%ny, geometry%nz))
    has_neighbor = .false.
    neighbor_i = 0
    neighbor_j = 0
    neighbor_k = 0
    neighborhood_count = 1
    alpha_self = 0.0_dp
    alpha_neighbor = 0.0_dp
    neighborhood_volume = 0.0_dp

    ! Choose an open, outward neighbor that meets the target by itself.
    ! Such receivers do not need their own merge; shared receivers retain the
    ! existing conservative neighborhood weights.
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          alpha_self(i, j, k) = 1.0_dp
          kappa = geometry%volume_fraction(i, j, k)
          if (kappa >= target) cycle
          if(present(box)) then
            if(.not.box%contains_valid([i,j,k],4).or..not.box%contains_domain([i,j,k],4)) cycle
          endif
          call select_outward_state_receiver( &
            geometry, i, j, k, target, ni, nj, nk, local_ok,periodic,image_offset,box)
          if (.not. local_ok) then
            call redistribute_connected_neighborhoods( &
              geometry, provisional_state, target, redistributed_state, ok,selected_order, &
              boundary_state,boundary_geometry,boundary_faces,registers,box,affected)
            return
          end if
          receiving_kappa = geometry%volume_fraction(ni, nj, nk)
          if (kappa + receiving_kappa + weight_tolerance < target) return
          neighbor_offset(:,i,j,k)=image_offset
          has_neighbor(i, j, k) = .true.
          neighbor_i(i, j, k) = ni
          neighbor_j(i, j, k) = nj
          neighbor_k(i, j, k) = nk
          alpha_neighbor(i, j, k) = &
            (target - kappa) / receiving_kappa
          neighborhood_count(ni, nj, nk) = &
            neighborhood_count(ni, nj, nk) + 1
        end do
      end do
    end do

    ! Each receiving cell reserves an equal share for every neighborhood to
    ! which it belongs, matching the zeroth-order weighted StateRedist rule.
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (.not. has_neighbor(i, j, k)) cycle
          ni = neighbor_i(i, j, k)
          nj = neighbor_j(i, j, k)
          nk = neighbor_k(i, j, k)
          alpha_self(ni, nj, nk) = alpha_self(ni, nj, nk) - &
            alpha_neighbor(i, j, k) / &
              real(neighborhood_count(ni, nj, nk), dp)
        end do
      end do
    end do

    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          if (.not. ieee_is_finite(alpha_self(i, j, k)) .or. &
              alpha_self(i, j, k) < -weight_tolerance .or. &
              alpha_neighbor(i, j, k) < 0.0_dp) return
          neighborhood_volume(i, j, k) = alpha_self(i, j, k) * &
            geometry%volume_fraction(i, j, k)
          if (has_neighbor(i, j, k)) then
            ni = neighbor_i(i, j, k)
            nj = neighbor_j(i, j, k)
            nk = neighbor_k(i, j, k)
            neighborhood_volume(i, j, k) = &
              neighborhood_volume(i, j, k) + &
              alpha_neighbor(i, j, k) * &
                geometry%volume_fraction(ni, nj, nk) / &
                real(neighborhood_count(ni, nj, nk), dp)
          end if
          if (.not. ieee_is_finite(neighborhood_volume(i, j, k)) .or. &
              neighborhood_volume(i, j, k) <= 0.0_dp) return
        end do
      end do
    end do

    allocate(neighborhood_state( &
      ncomp, geometry%nx, geometry%ny, geometry%nz))
    allocate(candidate(ncomp, geometry%nx, geometry%ny, geometry%nz))
    neighborhood_state = 0.0_dp
    if(present(box)) neighborhood_state=provisional_state
    candidate = 0.0_dp
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          if(present(box)) then
            if(.not.box%contains_valid([i,j,k],2).or..not.box%contains_domain([i,j,k],2)) cycle
          endif
          neighborhood_state(:, i, j, k) = alpha_self(i, j, k) * &
            geometry%volume_fraction(i, j, k) * &
            provisional_state(:, i, j, k)
          if (has_neighbor(i, j, k)) then
            ni = neighbor_i(i, j, k)
            nj = neighbor_j(i, j, k)
            nk = neighbor_k(i, j, k)
            neighborhood_state(:, i, j, k) = &
              neighborhood_state(:, i, j, k) + &
              alpha_neighbor(i, j, k) * &
                geometry%volume_fraction(ni, nj, nk) * &
                provisional_state(:, ni, nj, nk) / &
                real(neighborhood_count(ni, nj, nk), dp)
          end if
          neighborhood_state(:, i, j, k) = &
            neighborhood_state(:, i, j, k) / &
              neighborhood_volume(i, j, k)
        end do
      end do
    end do

    if(selected_order>0) then
      allocate(centroid(3,geometry%nx,geometry%ny,geometry%nz), &
        slopes(3,ncomp,geometry%nx,geometry%ny,geometry%nz),lower_bound(ncomp),upper_bound(ncomp), &
        stat=allocation_status)
      if(allocation_status/=0) return
      centroid=0
      do k=1,geometry%nz
      do j=1,geometry%ny
      do i=1,geometry%nx
        if(geometry%volume_fraction(i,j,k)<=0) cycle
        if(present(box)) then
          centroid(:,i,j,k)=state_cell_centroid(geometry,[i,j,k])
          if(.not.box%contains_valid([i,j,k],3).or..not.box%contains_domain([i,j,k],5)) cycle
        endif
        self_weight=alpha_self(i,j,k)*geometry%volume_fraction(i,j,k)
        centroid(:,i,j,k)=self_weight*state_cell_centroid(geometry,[i,j,k])
        if(has_neighbor(i,j,k)) then
          ni=neighbor_i(i,j,k);nj=neighbor_j(i,j,k);nk=neighbor_k(i,j,k)
          other_weight=alpha_neighbor(i,j,k)*geometry%volume_fraction(ni,nj,nk)/ &
            real(neighborhood_count(ni,nj,nk),dp)
          centroid(:,i,j,k)=centroid(:,i,j,k)+other_weight* &
            (real(neighbor_offset(:,i,j,k),dp)+state_cell_centroid(geometry,[ni,nj,nk]))
        endif
        centroid(:,i,j,k)=centroid(:,i,j,k)/neighborhood_volume(i,j,k)
      enddo
      enddo
      enddo
      call state_redistribution_slopes_3d(geometry,neighborhood_state,centroid,selected_order,slopes,local_ok, &
        boundary_state,boundary_geometry,boundary_faces,box)
      if(.not.local_ok) return
      do ni=1,ncomp
        if(present(box)) then
          lower_bound(ni)=minval(provisional_state(ni,box%valid_lower(1)-3:box%valid_upper(1)+3, &
            box%valid_lower(2)-3:box%valid_upper(2)+3,box%valid_lower(3)-3:box%valid_upper(3)+3), &
            mask=geometry%volume_fraction(box%valid_lower(1)-3:box%valid_upper(1)+3, &
            box%valid_lower(2)-3:box%valid_upper(2)+3,box%valid_lower(3)-3:box%valid_upper(3)+3)>0)
          upper_bound(ni)=maxval(provisional_state(ni,box%valid_lower(1)-3:box%valid_upper(1)+3, &
            box%valid_lower(2)-3:box%valid_upper(2)+3,box%valid_lower(3)-3:box%valid_upper(3)+3), &
            mask=geometry%volume_fraction(box%valid_lower(1)-3:box%valid_upper(1)+3, &
            box%valid_lower(2)-3:box%valid_upper(2)+3,box%valid_lower(3)-3:box%valid_upper(3)+3)>0)
        else
          lower_bound(ni)=minval(provisional_state(ni,:,:,:),mask=geometry%volume_fraction>0)
          upper_bound(ni)=maxval(provisional_state(ni,:,:,:),mask=geometry%volume_fraction>0)
        endif
      enddo
      do k=1,geometry%nz
      do j=1,geometry%ny
      do i=1,geometry%nx
        if(geometry%volume_fraction(i,j,k)<=0) cycle
        if(present(box)) then
          if(.not.box%contains_valid([i,j,k],1).or..not.box%contains_domain([i,j,k],2)) cycle
        endif
        recipient_count=1;position(:,1)=state_cell_centroid(geometry,[i,j,k])
        if(has_neighbor(i,j,k)) then
          ni=neighbor_i(i,j,k);nj=neighbor_j(i,j,k);nk=neighbor_k(i,j,k)
          recipient_count=2
          position(:,2)=real(neighbor_offset(:,i,j,k),dp)+state_cell_centroid(geometry,[ni,nj,nk])
        endif
        call limit_state_recipient_slopes(neighborhood_state(:,i,j,k),centroid(:,i,j,k), &
          slopes(:,:,i,j,k),position(:,:recipient_count),lower_bound,upper_bound)
      enddo
      enddo
      enddo
    endif

    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          if(present(box)) then
            if(.not.box%contains_valid([i,j,k],1).or..not.box%contains_domain([i,j,k],2)) cycle
          endif
          candidate(:, i, j, k) = candidate(:, i, j, k) + &
            alpha_self(i, j, k) * &
              real(neighborhood_count(i, j, k), dp) * &
              neighborhood_state(:, i, j, k)
          if(selected_order>0) then
            delta=state_cell_centroid(geometry,[i,j,k])-centroid(:,i,j,k)
            candidate(:,i,j,k)=candidate(:,i,j,k)+alpha_self(i,j,k)*real(neighborhood_count(i,j,k),dp)* &
              matmul(delta,slopes(:,:,i,j,k))
          endif
          if (has_neighbor(i, j, k)) then
            ni = neighbor_i(i, j, k)
            nj = neighbor_j(i, j, k)
            nk = neighbor_k(i, j, k)
            candidate(:, ni, nj, nk) = candidate(:, ni, nj, nk) + &
              alpha_neighbor(i, j, k) * &
                neighborhood_state(:, i, j, k)
            if(selected_order>0) then
              delta=real(neighbor_offset(:,i,j,k),dp)+state_cell_centroid(geometry,[ni,nj,nk])-centroid(:,i,j,k)
              candidate(:,ni,nj,nk)=candidate(:,ni,nj,nk)+alpha_neighbor(i,j,k)* &
                matmul(delta,slopes(:,:,i,j,k))
            endif
          end if
        end do
      end do
    end do
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          candidate(:, i, j, k) = candidate(:, i, j, k) / &
            real(neighborhood_count(i, j, k), dp)
        end do
      end do
    end do
    if(present(box)) call preserve_state_redist_patch_halo_3d(box,provisional_state,candidate)
    if (any(.not. ieee_is_finite(candidate))) return

    if(present(registers)) then
      do k=1,geometry%nz
      do j=1,geometry%ny
      do i=1,geometry%nx
        if(geometry%volume_fraction(i,j,k)<=0) cycle
        if(present(box)) then
          if(.not.box%contains_valid([i,j,k],1).or..not.box%contains_domain([i,j,k],2)) cycle
        endif
        recipient_count=1;register_cells(:,1)=[i,j,k]
        register_in(1)=alpha_self(i,j,k)*geometry%volume_fraction(i,j,k)/neighborhood_volume(i,j,k)
        register_out(1)=alpha_self(i,j,k);register_delta=0;register_slope=0
        if(selected_order>0) then
          register_delta(:,1)=state_cell_centroid(geometry,[i,j,k])-centroid(:,i,j,k)
          register_slope=slopes(:,:,i,j,k)
        endif
        if(has_neighbor(i,j,k)) then
          recipient_count=2;ni=neighbor_i(i,j,k);nj=neighbor_j(i,j,k);nk=neighbor_k(i,j,k)
          register_cells(:,2)=[ni,nj,nk]
          register_out(2)=alpha_neighbor(i,j,k)/real(neighborhood_count(ni,nj,nk),dp)
          register_in(2)=register_out(2)*geometry%volume_fraction(ni,nj,nk)/neighborhood_volume(i,j,k)
          if(selected_order>0) register_delta(:,2)=real(neighbor_offset(:,i,j,k),dp)+ &
            state_cell_centroid(geometry,[ni,nj,nk])-centroid(:,i,j,k)
        endif
        call register_state_neighborhood(geometry,provisional_state,register_cells(:,:recipient_count), &
          register_in(:recipient_count),register_out(:recipient_count),register_delta(:,:recipient_count), &
          register_slope,pending_registers,local_ok)
        if(.not.local_ok) return
      enddo
      enddo
      enddo
      registers=pending_registers
    endif
    if(present(affected)) then
      affected=has_neighbor.or.neighborhood_count>1
      if(present(box)) then
        do k=1,dims(3)
        do j=1,dims(2)
        do i=1,dims(1)
          if(.not.box%contains_valid([i,j,k],0)) affected(i,j,k)=.false.
        enddo
        enddo
        enddo
      endif
    endif
    redistributed_state = candidate
    ok = .true.
  end subroutine reactive_eb_weighted_state_redistribute_3d

  subroutine prepare_state_redist_box_register_3d(box,reg,nv,dims,ok)
    type(state_redist_box_3d), intent(in) :: box
    type(state_redist_register_3d), intent(inout) :: reg
    integer, intent(in) :: nv,dims(3)
    logical, intent(out) :: ok
    logical, allocatable :: owned(:,:,:)
    integer :: ios
    ok=.false.
    allocate(owned(dims(1),dims(2),dims(3)),source=.false.,stat=ios)
    if(ios/=0) return
    owned(box%valid_lower(1):box%valid_upper(1),box%valid_lower(2):box%valid_upper(2), &
      box%valid_lower(3):box%valid_upper(3))=.true.
    if(allocated(reg%valid_cells)) then
      if(any(reg%valid_cells.neqv.owned)) return
    else
      call move_alloc(owned,reg%valid_cells)
    endif
    ok=reg%is_valid(nv,dims)
  end subroutine prepare_state_redist_box_register_3d

  subroutine preserve_state_redist_patch_halo_3d(box,input,result)
    type(state_redist_box_3d), intent(in) :: box
    real(dp), intent(in) :: input(:,:,:,:)
    real(dp), intent(inout) :: result(:,:,:,:)
    integer :: i,j,k
    do k=1,size(input,4)
    do j=1,size(input,3)
    do i=1,size(input,2)
      if(.not.box%contains_valid([i,j,k],0)) result(:,i,j,k)=input(:,i,j,k)
    enddo
    enddo
    enddo
  end subroutine preserve_state_redist_patch_halo_3d

  ! Qhat is updated only in grow(valid,2); its outer layer and physical
  ! exterior samples are supplied U_in. cent_hat has independent grow-3 scope.
  subroutine patch_box_state_slopes_3d(g,qhat,centroid,order,slopes,ok,box,boundary_faces)
    type(eb_geometry_3d), intent(in) :: g
    real(dp), intent(in) :: qhat(:,:,:,:),centroid(:,:,:,:)
    integer, intent(in) :: order
    real(dp), intent(out) :: slopes(:,:,:,:,:)
    logical, intent(out) :: ok
    type(state_redist_box_3d), intent(in) :: box
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp) :: values(size(qhat,1),125),points(3,125),fit_points(3,125),centers(3,125)
    real(dp) :: widths(3,125),fractions(125),root_slopes(3,size(qhat,1))
    integer :: cell(3),other(3),off(3),bounds(2,3),i,j,k,ii,jj,kk,ns,a,face_index
    integer(int64) :: relative(2,3),clip
    character(len=16) :: faces(6)
    logical :: connected(5,5,5),near(125),extdir(6),face_sample(3),valid
    slopes=0;ok=.false.;faces='outflow';widths=1;clip=int(huge(1)/4,int64)
    if(present(boundary_faces)) faces=boundary_faces
    extdir=faces=='inflow'
    do k=box%valid_lower(3)-1,box%valid_upper(3)+1
    do j=box%valid_lower(2)-1,box%valid_upper(2)+1
    do i=box%valid_lower(1)-1,box%valid_upper(1)+1
      cell=[i,j,k]
      if(g%volume_fraction(i,j,k)<=0.or..not.box%contains_domain(cell,2)) cycle
      do a=1,3
        ! Far global boundaries can exceed a local default-integer index.
        ! Clamp before subtraction so even an int64 domain near its limit
        ! cannot overflow while forming this two-cell slope stencil.
        if(box%domain_lower(a)<int(cell(a),int64)-clip) then
          relative(1,a)=-clip
        else
          relative(1,a)=box%domain_lower(a)-int(cell(a),int64)
        endif
        if(box%domain_upper(a)>int(cell(a),int64)+clip) then
          relative(2,a)=clip
        else
          relative(2,a)=box%domain_upper(a)-int(cell(a),int64)
        endif
        if(box%periodic(a)) relative(:,a)=[-clip,clip]
      enddo
      bounds=int(max(-clip,min(clip,relative)))
      call connected_state_stencil(g,cell,[1,1,1],connected)
      ns=1;values(:,1)=qhat(:,i,j,k);centers(:,1)=0
      points(:,1)=centroid(:,i,j,k);fit_points(:,1)=points(:,1)
      fractions(1)=g%volume_fraction(i,j,k);near(1)=.true.
      do kk=-2,2
      do jj=-2,2
      do ii=-2,2
        off=[ii,jj,kk]
        if(all(off==0)) cycle
        other=cell+off
        if(g%volume_fraction(other(1),other(2),other(3))<=0) cycle
        ns=ns+1;centers(:,ns)=real(off,dp)
        fractions(ns)=g%volume_fraction(other(1),other(2),other(3))
        near(ns)=all(abs(off)<=1).and.connected(ii+3,jj+3,kk+3)
        values(:,ns)=qhat(:,other(1),other(2),other(3))
        points(:,ns)=real(off,dp)+centroid(:,other(1),other(2),other(3))
        fit_points(:,ns)=points(:,ns);face_sample=.false.
        do a=1,3
          face_sample(a)=(extdir(2*a-1).and.int(cell(a),int64)==box%domain_lower(a).and.off(a)==-1).or. &
            (extdir(2*a).and.int(cell(a),int64)==box%domain_upper(a).and.off(a)==1)
        enddo
        do a=1,3
          if(.not.face_sample(a)) cycle
          fit_points(a,ns)=.5_dp*real(off(a),dp)
          if(count(face_sample)/=1) cycle
          face_index=cell(a)+min(off(a),0)
          select case(a)
          case(1)
            fit_points(2,ns)=jj+g%x_face_centroid_y(face_index,other(2),other(3))
            fit_points(3,ns)=kk+g%x_face_centroid_z(face_index,other(2),other(3))
          case(2)
            fit_points(1,ns)=ii+g%y_face_centroid_x(other(1),face_index,other(3))
            fit_points(3,ns)=kk+g%y_face_centroid_z(other(1),other(2),face_index)
          case(3)
            fit_points(1,ns)=ii+g%z_face_centroid_x(other(1),other(2),face_index)
            fit_points(2,ns)=jj+g%z_face_centroid_y(other(1),other(2),face_index)
          end select
        enddo
      enddo
      enddo
      enddo
      call state_redistribution_graph_slopes_3d(values(:,:ns),points(:,:ns),centers(:,:ns),widths(:,:ns), &
        fractions(:ns),order,root_slopes,valid,near(:ns),extdir,bounds,fit_points(:,:ns))
      if(.not.valid) return
      slopes(:,:,i,j,k)=root_slopes
    enddo
    enddo
    enddo
    ok=.true.
  end subroutine patch_box_state_slopes_3d

  function state_cell_centroid(g,index) result(center)
    type(eb_geometry_3d), intent(in) :: g
    integer, intent(in) :: index(3)
    real(dp) :: center(3)
    center=[g%cell_centroid_x(index(1),index(2),index(3)), &
      g%cell_centroid_y(index(1),index(2),index(3)),g%cell_centroid_z(index(1),index(2),index(3))]
  end function state_cell_centroid

  ! Least-squares reconstruction at neighborhood centroids, followed by the
  ! reference regular-stencil override and centroid prediction limiter.
  subroutine state_redistribution_slopes_3d(g,qhat,centroid,order,slopes,ok, &
      boundary_state,boundary_geometry,boundary_faces,box)
    type(eb_geometry_3d), intent(in) :: g
    type(state_redist_box_3d), optional, intent(in) :: box
    real(dp), intent(in) :: qhat(:,:,:,:),centroid(:,:,:,:)
    integer, intent(in) :: order
    real(dp), intent(out) :: slopes(:,:,:,:,:)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: boundary_state(:,:,:,:)
    type(eb_geometry_3d), optional, intent(in) :: boundary_geometry
    character(len=*), optional, intent(in) :: boundary_faces(6)
    logical :: connected(5,5,5),near(5,5,5),valid,regular
    real(dp) :: matrix(3,3),inverse(3,3),rhs(3),delta(3),width_low(3),width_high(3),values(-2:2)
    real(dp) :: change,low,high,small,factor,limiter(3),value
    integer :: dims(3),cell(3),other(3),i,j,k,ii,jj,kk,a,b,c,radius(3),distance,selected
    slopes=0;ok=.false.;dims=[g%nx,g%ny,g%nz]
    if(present(box)) then
      call patch_box_state_slopes_3d(g,qhat,centroid,order,slopes,ok,box,boundary_faces)
      return
    endif
    if(present(boundary_state).or.present(boundary_geometry).or.present(boundary_faces)) then
      if(.not.present(boundary_state).or..not.present(boundary_geometry).or..not.present(boundary_faces)) return
      call boundary_state_slopes(g,qhat,centroid,order,slopes,ok,boundary_state,boundary_geometry,boundary_faces)
      return
    endif
    if(any(.not.ieee_is_finite(qhat)).or.any(.not.ieee_is_finite(centroid))) return
    do k=1,g%nz
    do j=1,g%ny
    do i=1,g%nx
      if(g%volume_fraction(i,j,k)<=0) cycle
      cell=[i,j,k];radius=1
      call connected_state_stencil(g,cell,radius,connected)
      near=connected;width_low=0;width_high=0
      do kk=-1,1
      do jj=-1,1
      do ii=-1,1
        if(.not.connected(ii+3,jj+3,kk+3)) cycle
        other=cell+[ii,jj,kk]
        delta=real([ii,jj,kk],dp)+centroid(:,other(1),other(2),other(3))-centroid(:,i,j,k)
        width_low=min(width_low,delta);width_high=max(width_high,delta)
      enddo
      enddo
      enddo
      radius=merge(2,1,width_high-width_low<.5_dp.and.dims>1)
      if(any(radius==2)) then
        ! AMReX's grown least-squares stencil uses every noncovered cell
        ! in the selected rectangular range; the centroid limiter retains
        ! the original connected one-cell stencil in near.
        connected=.false.
        do kk=-radius(3),radius(3)
        do jj=-radius(2),radius(2)
        do ii=-radius(1),radius(1)
          other=cell+[ii,jj,kk]
          if(any(other<1).or.any(other>dims)) cycle
          connected(ii+3,jj+3,kk+3)=g%volume_fraction(other(1),other(2),other(3))>0
        enddo
        enddo
        enddo
      endif
      matrix=0
      do kk=-radius(3),radius(3)
      do jj=-radius(2),radius(2)
      do ii=-radius(1),radius(1)
        if(.not.connected(ii+3,jj+3,kk+3)) cycle
        other=cell+[ii,jj,kk]
        delta=real([ii,jj,kk],dp)+centroid(:,other(1),other(2),other(3))-centroid(:,i,j,k)
        do a=1,3
        do b=1,3
          matrix(a,b)=matrix(a,b)+delta(a)*delta(b)
        enddo
        enddo
      enddo
      enddo
      enddo
      call state_slope_pseudoinverse(matrix,inverse)
      do c=1,size(qhat,1)
        rhs=0
        do kk=-radius(3),radius(3)
        do jj=-radius(2),radius(2)
        do ii=-radius(1),radius(1)
          if(.not.connected(ii+3,jj+3,kk+3)) cycle
          other=cell+[ii,jj,kk]
          delta=real([ii,jj,kk],dp)+centroid(:,other(1),other(2),other(3))-centroid(:,i,j,k)
          rhs=rhs+delta*(qhat(c,other(1),other(2),other(3))-qhat(c,i,j,k))
        enddo
        enddo
        enddo
        slopes(:,c,i,j,k)=matmul(inverse,rhs)
        do a=1,3
          selected=2
          if(order==4) selected=4
          do
            distance=selected/2;regular=.true.
            if(cell(a)<=distance.or.cell(a)+distance>dims(a)) regular=.false.
            if(regular) then
              do b=-distance,distance
                other=cell;other(a)=other(a)+b
                if(g%volume_fraction(other(1),other(2),other(3))/=1) regular=.false.
              enddo
            endif
            if(regular.or.selected==2) exit
            selected=2
          enddo
          if(.not.regular) cycle
          values=0
          do b=-distance,distance
            other=cell;other(a)=other(a)+b
            values(b)=qhat(c,other(1),other(2),other(3))
          enddo
          call amrex_limited_slope(values(-2),values(-1),values(0),values(1),values(2), &
            1._dp,selected,value,valid)
          if(.not.valid) return
          slopes(a,c,i,j,k)=value
        enddo
        limiter=1
        do kk=-1,1
        do jj=-1,1
        do ii=-1,1
          if(.not.near(ii+3,jj+3,kk+3)) cycle
          other=cell+[ii,jj,kk]
          delta=real([ii,jj,kk],dp)+centroid(:,other(1),other(2),other(3))-centroid(:,i,j,k)
          change=dot_product(delta,slopes(:,c,i,j,k))
          low=min(qhat(c,i,j,k),qhat(c,other(1),other(2),other(3)))
          high=max(qhat(c,i,j,k),qhat(c,other(1),other(2),other(3)))
          small=1.e-12_dp*max(abs(low),abs(high));factor=1
          if(change>small) factor=min(1._dp,(high-qhat(c,i,j,k))/change)
          if(change< -small) factor=min(1._dp,(low-qhat(c,i,j,k))/change)
          do a=1,3
            if(abs(delta(a))>1.e-12_dp) limiter(a)=min(limiter(a),max(0._dp,factor))
          enddo
        enddo
        enddo
        enddo
        slopes(:,c,i,j,k)=limiter*slopes(:,c,i,j,k)
      enddo
    enddo
    enddo
    enddo
    ok=all(ieee_is_finite(slopes))
  end subroutine state_redistribution_slopes_3d

  ! Physical boundary ghosts retain the pre-update state, matching the
  ! upstream zero dUdt outside domain_per_grown. Neighborhood scatter still
  ! uses only physical cells; these ghosts contribute exclusively to slopes.
  subroutine boundary_state_slopes(g,qhat,centroid,order,slopes,ok,halo,hg,faces)
    type(eb_geometry_3d), intent(in) :: g,hg
    real(dp), intent(in) :: qhat(:,:,:,:),centroid(:,:,:,:),halo(:,:,:,:)
    integer, intent(in) :: order
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(out) :: slopes(:,:,:,:,:)
    logical, intent(out) :: ok
    real(dp) :: values(size(qhat,1),125),points(3,125),fit_points(3,125),centers(3,125)
    real(dp) :: widths(3,125),fractions(125),root_slopes(3,size(qhat,1))
    integer :: dims(3),cell(3),other(3),ext(3),off(3),bounds(2,3),i,j,k,ii,jj,kk,ns,a,face_index
    logical :: connected(5,5,5),near(125),extdir(6),outside,face_sample(3),valid,periodic(3)
    slopes=0;ok=.false.;dims=[g%nx,g%ny,g%nz]
    if(.not.hg%is_valid()) return
    if(any([hg%nx,hg%ny,hg%nz]/=dims+4)) return
    if(maxval(abs([hg%dx-g%dx,hg%dy-g%dy,hg%dz-g%dz, &
      hg%x_lower-(g%x_lower-2*g%dx),hg%y_lower-(g%y_lower-2*g%dy),hg%z_lower-(g%z_lower-2*g%dz), &
      hg%x_upper-(g%x_upper+2*g%dx),hg%y_upper-(g%y_upper+2*g%dy),hg%z_upper-(g%z_upper+2*g%dz)])) &
      >512*epsilon(1._dp)*max(1._dp,abs(g%x_lower),abs(g%x_upper),abs(g%y_lower),abs(g%y_upper), &
      abs(g%z_lower),abs(g%z_upper))) return
    if(any(shape(halo)/=[size(qhat,1),dims+4])) return
    if(any(.not.ieee_is_finite(halo))) return
    do a=1,6
      select case(trim(faces(a)))
      case('outflow','inflow','slip_wall','no_slip_wall','isothermal_wall','periodic')
      case default
        return
      end select
    enddo
    if(.not.periodic_eb_faces_match(g,faces)) return
    if(.not.periodic_eb_faces_match(hg,faces,halo=2)) return
    periodic=faces(1:6:2)=='periodic'
    extdir=faces=='inflow';widths=1
    do k=1,dims(3)
    do j=1,dims(2)
    do i=1,dims(1)
      if(g%volume_fraction(i,j,k)<=0) cycle
      cell=[i,j,k];bounds(1,:)=1-cell;bounds(2,:)=dims-cell
      call connected_state_stencil(hg,cell+2,[1,1,1],connected)
      ns=1;values(:,1)=qhat(:,i,j,k);centers(:,1)=0
      points(:,1)=centroid(:,i,j,k);fit_points(:,1)=points(:,1)
      fractions(1)=g%volume_fraction(i,j,k);near(1)=.true.
      do kk=-2,2
      do jj=-2,2
      do ii=-2,2
        off=[ii,jj,kk]
        if(all(off==0)) cycle
        other=cell+off
        do a=1,3
          if(periodic(a)) other(a)=1+modulo(other(a)-1,dims(a))
        enddo
        ext=other+2
        if(hg%volume_fraction(ext(1),ext(2),ext(3))<=0) cycle
        outside=any(other<1).or.any(other>dims)
        ns=ns+1;centers(:,ns)=real(off,dp)
        fractions(ns)=hg%volume_fraction(ext(1),ext(2),ext(3))
        near(ns)=all(abs(off)<=1).and.connected(ii+3,jj+3,kk+3)
        if(outside) then
          values(:,ns)=halo(:,ext(1),ext(2),ext(3))
          points(:,ns)=real(off,dp)+[hg%cell_centroid_x(ext(1),ext(2),ext(3)), &
            hg%cell_centroid_y(ext(1),ext(2),ext(3)),hg%cell_centroid_z(ext(1),ext(2),ext(3))]
        else
          values(:,ns)=qhat(:,other(1),other(2),other(3))
          points(:,ns)=real(off,dp)+centroid(:,other(1),other(2),other(3))
        endif
        fit_points(:,ns)=points(:,ns);face_sample=.false.
        do a=1,3
          face_sample(a)=(extdir(2*a-1).and.cell(a)==1.and.off(a)==-1).or. &
            (extdir(2*a).and.cell(a)==dims(a).and.off(a)==1)
        enddo
        do a=1,3
          if(.not.face_sample(a)) cycle
          fit_points(a,ns)=.5_dp*real(off(a),dp)
          if(count(face_sample)/=1) cycle
          face_index=cell(a)+2+min(off(a),0)
          select case(a)
          case(1)
            fit_points(2,ns)=jj+hg%x_face_centroid_y(face_index,ext(2),ext(3))
            fit_points(3,ns)=kk+hg%x_face_centroid_z(face_index,ext(2),ext(3))
          case(2)
            fit_points(1,ns)=ii+hg%y_face_centroid_x(ext(1),face_index,ext(3))
            fit_points(3,ns)=kk+hg%y_face_centroid_z(ext(1),face_index,ext(3))
          case(3)
            fit_points(1,ns)=ii+hg%z_face_centroid_x(ext(1),ext(2),face_index)
            fit_points(2,ns)=jj+hg%z_face_centroid_y(ext(1),ext(2),face_index)
          end select
        enddo
      enddo
      enddo
      enddo
      call state_redistribution_graph_slopes_3d(values(:,:ns),points(:,:ns),centers(:,:ns),widths(:,:ns), &
        fractions(:ns),order,root_slopes,valid,near(:ns),extdir,bounds,fit_points(:,:ns))
      if(.not.valid) return
      slopes(:,:,i,j,k)=root_slopes
    enddo
    enddo
    enddo
    ok=.true.
  end subroutine boundary_state_slopes

  ! Append physical exterior samples to a composite leaf stencil. Coordinates
  ! are relative to the root lattice center; no exterior cell is a recipient.
  subroutine state_redistribution_boundary_nodes_3d(hg,halo,cell,global_cell,domain_dims,faces, &
      values,points,fit_points,centers,fractions,near,ns,ok)
    type(eb_geometry_3d), intent(in) :: hg
    real(dp), intent(in) :: halo(:,:,:,:)
    integer, intent(in) :: cell(3),global_cell(3),domain_dims(3)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(out) :: values(:,:),points(:,:),fit_points(:,:),centers(:,:),fractions(:)
    logical, intent(out) :: near(:),ok
    integer, intent(out) :: ns
    integer :: off(3),ext(3),outside(3),i,j,k,a,face_index
    logical :: connected(5,5,5),physical,face_sample(3),extdir(6)
    ns=0;ok=.false.
    if(.not.hg%is_valid()) return
    if(any(shape(halo)/=[size(values,1),hg%nx,hg%ny,hg%nz])) return
    if(any(.not.ieee_is_finite(halo)).or.any(global_cell<1).or.any(global_cell>domain_dims)) return
    if(any(cell<3).or.any(cell>[hg%nx,hg%ny,hg%nz]-2)) return
    if(size(values,2)<125.or.any(shape(points)/=[3,125]).or.any(shape(fit_points)/=[3,125]).or. &
      any(shape(centers)/=[3,125]).or.size(fractions)<125.or.size(near)<125) return
    do a=1,6
      select case(trim(faces(a)))
      case('outflow','inflow','slip_wall','no_slip_wall','isothermal_wall','periodic')
      case default
        return
      end select
    enddo
    do a=1,3
      if((faces(2*a-1)=='periodic').neqv.(faces(2*a)=='periodic')) return
    enddo
    extdir=faces=='inflow'
    call connected_state_stencil(hg,cell,[1,1,1],connected)
    do k=-2,2
    do j=-2,2
    do i=-2,2
      off=[i,j,k];outside=global_cell+off;physical=.false.
      do a=1,3
        if(outside(a)<1.and.faces(2*a-1)/='periodic') physical=.true.
        if(outside(a)>domain_dims(a).and.faces(2*a)/='periodic') physical=.true.
      enddo
      if(.not.physical) cycle
      ext=cell+off
      if(hg%volume_fraction(ext(1),ext(2),ext(3))<=0) cycle
      ns=ns+1;values(:,ns)=halo(:,ext(1),ext(2),ext(3));centers(:,ns)=real(off,dp)
      points(:,ns)=real(off,dp)+state_cell_centroid(hg,ext)
      fractions(ns)=hg%volume_fraction(ext(1),ext(2),ext(3))
      near(ns)=all(abs(off)<=1).and.connected(i+3,j+3,k+3)
      fit_points(:,ns)=points(:,ns);face_sample=.false.
      do a=1,3
        face_sample(a)=(extdir(2*a-1).and.global_cell(a)==1.and.off(a)==-1).or. &
          (extdir(2*a).and.global_cell(a)==domain_dims(a).and.off(a)==1)
      enddo
      do a=1,3
        if(.not.face_sample(a)) cycle
        fit_points(a,ns)=.5_dp*real(off(a),dp)
        if(count(face_sample)/=1) cycle
        face_index=cell(a)+min(off(a),0)
        select case(a)
        case(1)
          fit_points(2,ns)=j+hg%x_face_centroid_y(face_index,ext(2),ext(3))
          fit_points(3,ns)=k+hg%x_face_centroid_z(face_index,ext(2),ext(3))
        case(2)
          fit_points(1,ns)=i+hg%y_face_centroid_x(ext(1),face_index,ext(3))
          fit_points(3,ns)=k+hg%y_face_centroid_z(ext(1),face_index,ext(3))
        case(3)
          fit_points(1,ns)=i+hg%z_face_centroid_x(ext(1),ext(2),face_index)
          fit_points(2,ns)=j+hg%z_face_centroid_y(ext(1),ext(2),face_index)
        end select
      enddo
    enddo
    enddo
    enddo
    ok=.true.
  end subroutine state_redistribution_boundary_nodes_3d

  subroutine connected_state_stencil(g,cell,radius,connected)
    type(eb_geometry_3d), intent(in) :: g
    integer, intent(in) :: cell(3),radius(3)
    logical, intent(out) :: connected(5,5,5)
    integer :: queue(3,125),head,tail,index(3),other(3),offset(3),dims(3),a,step
    real(dp) :: area
    dims=[g%nx,g%ny,g%nz];connected=.false.;connected(3,3,3)=.true.
    head=1;tail=1;queue(:,1)=cell
    do while(head<=tail)
      index=queue(:,head)
      do a=1,3
      do step=-1,1,2
        other=index;other(a)=other(a)+step;offset=other-cell
        if(any(abs(offset)>radius).or.any(other<1).or.any(other>dims)) cycle
        if(connected(offset(1)+3,offset(2)+3,offset(3)+3)) cycle
        if(g%volume_fraction(other(1),other(2),other(3))<=0) cycle
        select case(a)
        case(1);area=g%x_face_fraction(index(1)+min(step,0),index(2),index(3))
        case(2);area=g%y_face_fraction(index(1),index(2)+min(step,0),index(3))
        case(3);area=g%z_face_fraction(index(1),index(2),index(3)+min(step,0))
        end select
        if(area<=0) cycle
        connected(offset(1)+3,offset(2)+3,offset(3)+3)=.true.
        tail=tail+1;queue(:,tail)=other
      enddo
      enddo
      head=head+1
    enddo
  end subroutine connected_state_stencil

  ! The first column is the root. Coordinates and widths are measured in the
  ! root cell's lattice units, with periodic images already chosen by the caller.
  ! Leaf graphs let the same reconstruction cross patch and coarse/fine faces.
  ! near_connected marks the original connected stencil for extent and limiting.
  ! Grown fits use all supplied noncovered leaves inside the selected range.
  subroutine state_redistribution_graph_slopes_3d(qhat,centroid,centers,widths,fractions,order,slopes,ok,near_connected, &
      extdir_faces,extdir_bounds,fit_centroid)
    real(dp), intent(in) :: qhat(:,:),centroid(:,:),centers(:,:),widths(:,:),fractions(:)
    integer, intent(in) :: order
    real(dp), intent(out) :: slopes(:,:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: near_connected(:)
    logical, optional, intent(in) :: extdir_faces(6)
    ! Domain indices relative to the root (index zero), per axis. Supplied
    ! ghost Qhat values are boundary samples; fit_centroid can supply face
    ! positions independently of cent_hat used by the centroid limiter.
    integer, optional, intent(in) :: extdir_bounds(2,3)
    real(dp), optional, intent(in) :: fit_centroid(:,:)
    real(dp) :: matrix(3,3),inverse(3,3),rhs(3),delta(3),extent_low(3),extent_high(3),radius(3)
    real(dp) :: values(-2:2),target(3),change,low,high,small,factor,limiter(3),value,tolerance
    logical :: near(size(qhat,2)),fit(size(qhat,2)),regular,valid
    integer :: ns,nv,a,b,c,m,selected,distance,match
    slopes=0;ok=.false.;ns=size(qhat,2);nv=size(qhat,1)
    if(ns<1.or.nv<1.or.any(shape(slopes)/=[3,nv])) return
    if(any(shape(centroid)/=[3,ns]).or.any(shape(centers)/=[3,ns]).or. &
      any(shape(widths)/=[3,ns]).or.size(fractions)/=ns) return
    if(order/=2.and.order/=4) return
    if(any(.not.ieee_is_finite(qhat)).or.any(.not.ieee_is_finite(centroid)).or. &
      any(.not.ieee_is_finite(centers)).or.any(.not.ieee_is_finite(widths))) return
    if(any(widths<=0).or.any(.not.ieee_is_finite(fractions))) return
    if(present(near_connected)) then
      if(size(near_connected)/=ns) return
    endif
    if(present(extdir_faces).neqv.present(extdir_bounds)) return
    if(present(extdir_bounds)) then
      if(any(extdir_bounds(1,:)>0).or.any(extdir_bounds(2,:)<0)) return
    endif
    if(present(fit_centroid)) then
      if(any(shape(fit_centroid)/=[3,ns])) return
      if(any(.not.ieee_is_finite(fit_centroid))) return
    endif
    tolerance=512*epsilon(1._dp)
    extent_low=0;extent_high=0
    do m=1,ns
      near(m)=all(abs(centers(:,m)-centers(:,1))<=1+tolerance)
      if(present(near_connected)) near(m)=near(m).and.near_connected(m)
      if(.not.near(m)) cycle
      delta=centroid(:,m)-centroid(:,1)
      extent_low=min(extent_low,delta);extent_high=max(extent_high,delta)
    enddo
    radius=merge(2._dp,1._dp,extent_high-extent_low<.5_dp)
    matrix=0
    do m=1,ns
      fit(m)=all(abs(centers(:,m)-centers(:,1))<=radius+tolerance)
      if(.not.any(radius>1)) fit(m)=fit(m).and.near(m)
      if(.not.fit(m)) cycle
      delta=centroid(:,m)-centroid(:,1)
      if(present(fit_centroid)) delta=fit_centroid(:,m)-fit_centroid(:,1)
      do a=1,3
      do b=1,3
        matrix(a,b)=matrix(a,b)+delta(a)*delta(b)
      enddo
      enddo
    enddo
    call state_slope_pseudoinverse(matrix,inverse)
    do c=1,nv
      rhs=0
      do m=1,ns
        if(.not.fit(m)) cycle
        delta=centroid(:,m)-centroid(:,1)
        if(present(fit_centroid)) delta=fit_centroid(:,m)-fit_centroid(:,1)
        rhs=rhs+delta*(qhat(c,m)-qhat(c,1))
      enddo
      slopes(:,c)=matmul(inverse,rhs)
      do a=1,3
        selected=order
        do
          distance=selected/2;regular=fractions(1)==1.and.all(abs(widths(:,1)-1)<tolerance)
          values=0
          do b=-distance,distance
            target=centers(:,1);target(a)=target(a)+b;match=0
            do m=1,ns
              if(any(abs(centers(:,m)-target)>tolerance)) cycle
              if(fractions(m)/=1.or.any(abs(widths(:,m)-1)>tolerance)) cycle
              match=m;exit
            enddo
            if(match==0) then
              regular=.false.
            else
              values(b)=qhat(c,match)
            endif
          enddo
          if(regular.or.selected==2) exit
          selected=2
        enddo
        if(.not.regular) cycle
        if(present(extdir_faces)) then
          call amrex_extdir_limited_slope(values(-2),values(-1),values(0),values(1),values(2), &
            selected,0,extdir_bounds(1,a),extdir_bounds(2,a),extdir_faces(2*a-1),extdir_faces(2*a),value,valid)
        else
          call amrex_limited_slope(values(-2),values(-1),values(0),values(1),values(2), &
            1._dp,selected,value,valid)
        endif
        if(.not.valid) return
        slopes(a,c)=value
      enddo
      limiter=1
      do m=1,ns
        if(.not.near(m)) cycle
        delta=centroid(:,m)-centroid(:,1);change=dot_product(delta,slopes(:,c))
        low=min(qhat(c,1),qhat(c,m));high=max(qhat(c,1),qhat(c,m))
        small=1.e-12_dp*max(abs(low),abs(high));factor=1
        if(change>small) factor=min(1._dp,(high-qhat(c,1))/change)
        if(change< -small) factor=min(1._dp,(low-qhat(c,1))/change)
        do a=1,3
          if(abs(delta(a))>1.e-12_dp) limiter(a)=min(limiter(a),max(0._dp,factor))
        enddo
      enddo
      slopes(:,c)=limiter*slopes(:,c)
    enddo
    ok=all(ieee_is_finite(slopes))
  end subroutine state_redistribution_graph_slopes_3d

  ! Resolve rank-deficient extruded stencils without introducing a spurious
  ! gradient in an unobserved direction. Jacobi diagonalization is local 3x3.
  subroutine state_slope_pseudoinverse(matrix,inverse)
    real(dp), intent(in) :: matrix(3,3)
    real(dp), intent(out) :: inverse(3,3)
    real(dp) :: diagonal(3,3),vectors(3,3),rotation(3,3),angle,scale,largest
    integer :: iteration,p,q,a,b,n
    diagonal=matrix;vectors=0
    do a=1,3
      vectors(a,a)=1
    enddo
    scale=max(1._dp,maxval(abs(matrix)))
    do iteration=1,32
      p=1;q=2;largest=abs(diagonal(p,q))
      do a=1,2
      do b=a+1,3
        if(abs(diagonal(a,b))<=largest) cycle
        p=a;q=b;largest=abs(diagonal(a,b))
      enddo
      enddo
      if(largest<=64*epsilon(1._dp)*scale) exit
      angle=.5_dp*atan2(2*diagonal(p,q),diagonal(q,q)-diagonal(p,p))
      rotation=0
      do a=1,3
        rotation(a,a)=1
      enddo
      rotation(p,p)=cos(angle);rotation(q,q)=cos(angle)
      rotation(p,q)=sin(angle);rotation(q,p)=-sin(angle)
      diagonal=matmul(transpose(rotation),matmul(diagonal,rotation))
      vectors=matmul(vectors,rotation)
    enddo
    inverse=0
    do n=1,3
      if(diagonal(n,n)<=4096*epsilon(1._dp)*scale) cycle
      do a=1,3
      do b=1,3
        inverse(a,b)=inverse(a,b)+vectors(a,n)*vectors(b,n)/diagonal(n,n)
      enddo
      enddo
    enddo
  end subroutine state_slope_pseudoinverse

  ! One factor for every recipient of a neighborhood preserves its zero
  ! weighted linear moment while keeping predictions in the original bounds.
  subroutine limit_state_recipient_slopes(qhat,centroid,slopes,positions,low,high)
    real(dp), intent(in) :: qhat(:),centroid(3),positions(:,:),low(:),high(:)
    real(dp), intent(inout) :: slopes(:,:)
    real(dp) :: change,predicted,small,factor
    integer :: c,recipient
    do c=1,size(qhat)
      factor=1;small=1.e-12_dp*max(abs(low(c)),abs(high(c)))
      do recipient=1,size(positions,2)
        change=dot_product(positions(:,recipient)-centroid,slopes(:,c));predicted=qhat(c)+change
        if(predicted>high(c)+small.and.change>0) factor=min(factor,(high(c)-qhat(c))/change)
        if(predicted<low(c)-small.and.change<0) factor=min(factor,(low(c)-qhat(c))/change)
      enddo
      slopes(:,c)=max(0._dp,min(1._dp,factor))*slopes(:,c)
    enddo
  end subroutine limit_state_recipient_slopes

  ! Pinned AMReX 1546c1cc MakeITracker: normal order, symmetry expansion,
  ! and sequential physical-boundary overrides. Periodic images are explicit.
  subroutine select_state_redist_tracker_3d(g,cell,target,box,neighbors,count,ok)
    type(eb_geometry_3d), intent(in) :: g
    integer, intent(in) :: cell(3)
    real(dp), intent(in) :: target
    type(state_redist_box_3d), intent(in) :: box
    integer, intent(out) :: neighbors(3,7),count
    logical, intent(out) :: ok
    integer :: a,b,c,step(3),index(3),m
    real(dp) :: normal(3),norm,volume
    logical :: equal_norm(3,3),allowed(3,2),symmetry
    neighbors=0;count=0;ok=.false.
    if(any(cell<1).or.any(cell>[g%nx,g%ny,g%nz])) return
    if(.not.ieee_is_finite(target).or.target<=0.or.target>1) return
    volume=g%volume_fraction(cell(1),cell(2),cell(3))
    if(volume<=0.or.volume>=target) then
      ok=.true.;return
    endif
    normal=[g%x_face_fraction(cell(1),cell(2),cell(3))-g%x_face_fraction(cell(1)-1,cell(2),cell(3)), &
      g%y_face_fraction(cell(1),cell(2),cell(3))-g%y_face_fraction(cell(1),cell(2)-1,cell(3)), &
      g%z_face_fraction(cell(1),cell(2),cell(3))-g%z_face_fraction(cell(1),cell(2),cell(3)-1)]
    norm=sqrt(sum(normal**2))
    if(.not.ieee_is_finite(norm).or.norm<=0) return
    normal=normal/norm
    do a=1,3
      step(a)=merge(1,-1,normal(a)>0)
      allowed(a,1)=box%periodic(a).or.int(cell(a),int64)>box%domain_lower(a)
      allowed(a,2)=box%periodic(a).or.int(cell(a),int64)<box%domain_upper(a)
      do b=1,3
        equal_norm(a,b)=abs(normal(a)-normal(b))<1.e-8_dp.or.abs(normal(a)+normal(b))<1.e-8_dp
      enddo
    enddo
    if(abs(normal(1))>abs(normal(2)).and.abs(normal(1))>abs(normal(3))) then
      a=1
    else if(abs(normal(2))>=abs(normal(1)).and.abs(normal(2))>abs(normal(3))) then
      a=2
    else
      a=3
    endif
    if(a==1) then
      if(.not.allowed(1,merge(2,1,step(1)>0))) a=merge(2,3,abs(normal(2))>abs(normal(3)))
    endif
    if(a==2) then
      if(.not.allowed(2,merge(2,1,step(2)>0))) a=merge(1,3,abs(normal(1))>abs(normal(3)))
    endif
    if(a==3) then
      if(.not.allowed(3,merge(2,1,step(3)>0))) a=merge(1,2,abs(normal(1))>abs(normal(2)))
    endif
    neighbors(a,1)=step(a);count=1
    index=cell+neighbors(:,1)
    if(any(index<1).or.any(index>[g%nx,g%ny,g%nz])) return
    if(g%volume_fraction(index(1),index(2),index(3))<=0) return
    volume=volume+g%volume_fraction(index(1),index(2),index(3))
    symmetry=.false.
    do b=1,3
      if(b/=a) symmetry=symmetry.or.equal_norm(a,b)
    enddo
    if(volume<target.or.symmetry) then
      select case(a)
      case(1)
        if(equal_norm(1,2)) then
          b=2
        else if(equal_norm(1,3)) then
          b=3
        else
          b=merge(2,3,abs(normal(2))>abs(normal(3)))
        endif
      case(2)
        if(equal_norm(1,2)) then
          b=1
        else if(equal_norm(2,3)) then
          b=3
        else
          b=merge(1,3,abs(normal(1))>abs(normal(3)))
        endif
      case(3)
        if(equal_norm(1,3)) then
          b=1
        else if(equal_norm(2,3)) then
          b=2
        else
          b=merge(1,2,abs(normal(1))>abs(normal(2)))
        endif
      end select
      neighbors(b,2)=step(b)
      neighbors(:,3)=neighbors(:,1)+neighbors(:,2);count=3
      do m=2,3
        index=cell+neighbors(:,m)
        if(any(index<1).or.any(index>[g%nx,g%ny,g%nz])) return
        volume=volume+g%volume_fraction(index(1),index(2),index(3))
      enddo
      c=6-a-b
      symmetry=equal_norm(a,c).or.equal_norm(b,c)
      if(volume<target.or.symmetry) then
        neighbors(c,4)=step(c)
        neighbors(:,5)=neighbors(:,1)+neighbors(:,4)
        neighbors(:,6)=neighbors(:,2)+neighbors(:,4)
        ! MakeITracker orders these corners by coordinate axis, independent
        ! of which face neighbor was selected first.
        if(a>b) then
          neighbors(:,5)=neighbors(:,2)+neighbors(:,4)
          neighbors(:,6)=neighbors(:,1)+neighbors(:,4)
        endif
        neighbors(:,7)=neighbors(:,3)+neighbors(:,4);count=7
        do m=4,7
          index=cell+neighbors(:,m)
          if(any(index<1).or.any(index>[g%nx,g%ny,g%nz])) return
          volume=volume+g%volume_fraction(index(1),index(2),index(3))
        enddo
      endif
    endif
    ok=ieee_is_finite(volume).and.volume>=target
  end subroutine select_state_redist_tracker_3d

  ! Fallback for connected clusters without a target-sized outward receiver.
  ! Every active cell owns a neighborhood; large cells retain a singleton.
  ! Equal membership shares make both averaging passes conservative and convex.
  ! The growth threshold is the raw volume, not the overlap-adjusted volume;
  ! this construction alone does not establish a general small-cell CFL bound.
  subroutine redistribute_connected_neighborhoods(geometry, state, target, result, ok,order, &
      boundary_state,boundary_geometry,boundary_faces,registers,box,affected)
    type(state_redist_register_3d), optional, intent(inout) :: registers
    type(state_redist_box_3d), optional, intent(in) :: box
    logical, optional, intent(inout) :: affected(:,:,:)
    type(state_redist_register_3d) :: pending_registers
    integer, allocatable :: register_cells(:,:)
    real(dp), allocatable :: register_in(:),register_out(:),register_delta(:,:),register_slope(:,:)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: state(:,:,:,:), target
    real(dp), intent(out) :: result(:,:,:,:)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: boundary_state(:,:,:,:)
    type(eb_geometry_3d), optional, intent(in) :: boundary_geometry
    character(len=*), optional, intent(in) :: boundary_faces(6)
    integer, intent(in) :: order
    real(dp), allocatable :: qhat(:,:,:,:),centroid(:,:,:,:),slopes(:,:,:,:,:),positions(:,:),low(:),high(:)
    real(dp) :: center(3),delta(3)
    integer :: root_index(3),component,ncomp
    logical :: valid
    type :: neighborhood
      integer, allocatable :: member(:),image(:,:)
    end type neighborhood
    type(neighborhood), allocatable :: groups(:)
    integer, allocatable :: queue(:), seen(:), memberships(:),queue_image(:,:)
    logical :: periodic(3),compute_average,compute_center,use_tracker,crosses_physical,on_boundary
    real(dp), allocatable :: kappa(:), candidate(:,:,:,:), average(:),alpha_self(:),alpha_neighbor(:)
    real(dp) :: volume, weight, aperture,neighbor_volume,boundary_normal(3)
    integer :: nx, ny, nz, ncell, root, head, tail, cell, other, m, count, allocation_status
    integer :: i, j, k, axis, offset, index(3), dims(3),tracker(3,7),tracker_count

    ! Caller has checked shapes, geometry, target and finite input.
    result=0._dp
    ok=.false.
    nx=geometry%nx; ny=geometry%ny; nz=geometry%nz
    if (nx>huge(ncell)/ny) return
    ncell=nx*ny
    if (ncell>huge(ncell)/nz) return
    ncell=ncell*nz
    dims=[nx,ny,nz];periodic=.false.
    if(present(boundary_faces).and..not.present(box)) periodic=boundary_faces(1:6:2)=='periodic'
    allocate(queue_image(3,ncell))
    allocate(groups(ncell),queue(ncell),seen(ncell),memberships(ncell), &
      kappa(ncell),candidate(size(state,1),nx,ny,nz),average(size(state,1)),stat=allocation_status)
    if (allocation_status/=0) return
    if(present(registers)) then
      if(.not.registers%is_valid(size(state,1),dims)) return
      pending_registers=registers
      if(present(box)) then
        call prepare_state_redist_box_register_3d(box,pending_registers,size(state,1),dims,valid)
        if(.not.valid) return
      endif
      allocate(register_cells(3,ncell),register_in(ncell),register_out(ncell), &
        register_delta(3,ncell),register_slope(3,size(state,1)))
    endif
    kappa=reshape(geometry%volume_fraction,[ncell])
    if(present(box)) then
      do root=1,ncell
        call decode(root,i,j,k)
        if(.not.box%contains_valid([i,j,k],5).or..not.box%contains_domain([i,j,k],5)) kappa(root)=0
      enddo
    endif
    memberships=merge(1,0,kappa>0._dp)
    ! Native weighted StateRedist uses the same overlap weights for a
    ! physical-domain array and a prepared patch. Exterior cells are excluded
    ! by kappa, rather than changing the neighborhood averaging convention.
    allocate(alpha_self(ncell),alpha_neighbor(ncell))
    alpha_self=merge(1._dp,0._dp,kappa>0);alpha_neighbor=alpha_self
    seen=0
    do root=1,ncell
      if (kappa(root)<=0._dp.or.kappa(root)>=target) cycle
      if(present(box)) then
        call decode(root,i,j,k)
        if(.not.box%contains_valid([i,j,k],4).or..not.box%contains_domain([i,j,k],4)) cycle
      endif
      use_tracker=.false.
      if(present(box)) then
        call decode(root,i,j,k)
        call select_state_redist_tracker_3d(geometry,[i,j,k],target,box,tracker,tracker_count,valid)
        on_boundary=any((int([i,j,k],int64)==box%domain_lower.or. &
          int([i,j,k],int64)==box%domain_upper).and..not.box%periodic)
        if(.not.valid) then
          ! Interior native tracker failures retain their existing contract.
          ! At a physical face the normal-directed receiver can be exterior
          ! or covered, so search the connected physical volume instead.
          if(.not.on_boundary) return
          boundary_normal=[geometry%x_face_fraction(i,j,k)-geometry%x_face_fraction(i-1,j,k), &
            geometry%y_face_fraction(i,j,k)-geometry%y_face_fraction(i,j-1,k), &
            geometry%z_face_fraction(i,j,k)-geometry%z_face_fraction(i,j,k-1)]
          if(all(boundary_normal==0)) return
        else
          if(tracker_count<1) return
          crosses_physical=.false.
          do m=1,tracker_count
            index=[i,j,k]+tracker(:,m)
            if(box%contains_domain(index,5)) cycle
            if(geometry%volume_fraction(index(1),index(2),index(3))>0) crosses_physical=.true.
          end do
          tail=1;queue(1)=root;queue_image(:,1)=0
          volume=kappa(root);neighbor_volume=0
          do m=1,tracker_count
            index=[i,j,k]+tracker(:,m)
            ! A physical exterior is a slope boundary sample, never a
            ! neighborhood inventory or a redistribution recipient.
            if(crosses_physical.and..not.box%contains_domain(index,5)) cycle
            other=index(1)+nx*((index(2)-1)+ny*(index(3)-1))
            tail=tail+1;queue(tail)=other;queue_image(:,tail)=tracker(:,m)
            volume=volume+kappa(other);neighbor_volume=neighbor_volume+kappa(other)
          end do
          use_tracker=volume>=target
          if(.not.use_tracker.and..not.crosses_physical) return
        end if
      endif
      if(.not.use_tracker) then
      call decode(root,i,j,k);root_index=[i,j,k]
      head=1; tail=1; queue(1)=root; seen(root)=root;queue_image(:,1)=0
      volume=kappa(root)
      do while (volume<target.and.head<=tail)
        cell=queue(head)
        call decode(cell,i,j,k)
        ! Deterministic breadth-first order: x-, x+, y-, y+, z-, z+.
        do axis=1,3
          do offset=-1,1,2
            if (volume>=target) exit
            index=[i,j,k]; index(axis)=index(axis)+offset
            if(periodic(axis)) index(axis)=1+modulo(index(axis)-1,dims(axis))
            if (any(index<1).or.any(index>dims)) cycle
            if(present(box)) then
              ! Keep the prepared grow1 scatter/support contract. A physical
              ! boundary may use either side of each axis, but no distant
              ! root may contribute without a matching stage donor halo.
              if(any(abs(index-root_index)>1)) cycle
            end if
            other=index(1)+nx*((index(2)-1)+ny*(index(3)-1))
            if (seen(other)==root.or.kappa(other)<=0._dp) cycle
            select case(axis)
            case(1)
              aperture=geometry%x_face_fraction(i+min(offset,0),j,k)
            case(2)
              aperture=geometry%y_face_fraction(i,j+min(offset,0),k)
            case(3)
              aperture=geometry%z_face_fraction(i,j,k+min(offset,0))
            end select
            if (aperture<=0._dp) cycle
            tail=tail+1; queue(tail)=other; seen(other)=root
            queue_image(:,tail)=queue_image(:,head);queue_image(axis,tail)=queue_image(axis,tail)+offset
            volume=volume+kappa(other)
          end do
        end do
        head=head+1
      end do
      endif
      if(.not.use_tracker) then
        neighbor_volume=0
        do m=2,tail
          neighbor_volume=neighbor_volume+kappa(queue(m))
        end do
      end if
      if(neighbor_volume<=0) return
      alpha_neighbor(root)=(target-kappa(root))/neighbor_volume
      if (volume<target) return
      allocate(groups(root)%member(tail),groups(root)%image(3,tail),stat=allocation_status)
      if (allocation_status/=0) return
      groups(root)%member=queue(:tail);groups(root)%image=queue_image(:,:tail)
      do m=2,tail
        other=queue(m)
        memberships(other)=memberships(other)+1
      end do
    end do

    do root=1,ncell
      if(.not.allocated(groups(root)%member)) cycle
      do m=2,size(groups(root)%member)
        other=groups(root)%member(m)
        alpha_self(other)=alpha_self(other)-alpha_neighbor(root)/real(memberships(other),dp)
      enddo
    enddo
    candidate=0._dp
    ncomp=size(state,1)
    if(order>0) then
      allocate(qhat(ncomp,nx,ny,nz),centroid(3,nx,ny,nz),slopes(3,ncomp,nx,ny,nz), &
        positions(3,ncell),low(ncomp),high(ncomp),stat=allocation_status)
      if(allocation_status/=0) return
      qhat=0;centroid=0
      if(present(box)) then
        qhat=state
        do k=1,nz
        do j=1,ny
        do i=1,nx
          centroid(:,i,j,k)=state_cell_centroid(geometry,[i,j,k])
        enddo
        enddo
        enddo
      endif
      do component=1,ncomp
        if(present(box)) then
          low(component)=minval(state(component,box%valid_lower(1)-3:box%valid_upper(1)+3, &
            box%valid_lower(2)-3:box%valid_upper(2)+3,box%valid_lower(3)-3:box%valid_upper(3)+3), &
            mask=geometry%volume_fraction(box%valid_lower(1)-3:box%valid_upper(1)+3, &
            box%valid_lower(2)-3:box%valid_upper(2)+3,box%valid_lower(3)-3:box%valid_upper(3)+3)>0)
          high(component)=maxval(state(component,box%valid_lower(1)-3:box%valid_upper(1)+3, &
            box%valid_lower(2)-3:box%valid_upper(2)+3,box%valid_lower(3)-3:box%valid_upper(3)+3), &
            mask=geometry%volume_fraction(box%valid_lower(1)-3:box%valid_upper(1)+3, &
            box%valid_lower(2)-3:box%valid_upper(2)+3,box%valid_lower(3)-3:box%valid_upper(3)+3)>0)
        else
          low(component)=minval(state(component,:,:,:),mask=geometry%volume_fraction>0)
          high(component)=maxval(state(component,:,:,:),mask=geometry%volume_fraction>0)
        endif
      enddo
    endif
    do root=1,ncell
      if (kappa(root)<=0._dp) cycle
      count=1
      if (allocated(groups(root)%member)) count=size(groups(root)%member)
      volume=0._dp; average=0._dp;center=0
      call decode(root,i,j,k);root_index=[i,j,k]
      compute_average=.true.;compute_center=order>0
      if(present(box)) then
        if(order>0) then
          compute_average=box%contains_valid(root_index,2).and.box%contains_domain(root_index,2)
          compute_center=box%contains_valid(root_index,3).and.box%contains_domain(root_index,5)
          if(.not.compute_average.and..not.compute_center) cycle
        else
          if(.not.box%contains_valid(root_index,1).or..not.box%contains_domain(root_index,2)) cycle
        endif
      endif
      do m=1,count
        cell=root
        if (allocated(groups(root)%member)) cell=groups(root)%member(m)
        call decode(cell,i,j,k)
        weight=kappa(cell)*recipient_weight(root,cell,m)
        volume=volume+weight
        if(compute_average) average=average+weight*state(:,i,j,k)
        if(compute_center) then
          index=0
          if(allocated(groups(root)%member)) index=groups(root)%image(:,m)
          center=center+weight*(real(index,dp)+state_cell_centroid(geometry,[i,j,k]))
        endif
      end do
      if (volume<=0._dp.or..not.ieee_is_finite(volume)) return
      if(compute_average) average=average/volume
      if(order>0) then
        call decode(root,i,j,k)
        if(present(box)) then
          if(box%contains_valid([i,j,k],2).and.box%contains_domain([i,j,k],2)) qhat(:,i,j,k)=average
          if(box%contains_valid([i,j,k],3).and.box%contains_domain([i,j,k],5)) centroid(:,i,j,k)=center/volume
        else
          qhat(:,i,j,k)=average;centroid(:,i,j,k)=center/volume
        endif
        cycle
      endif
      if(present(box)) then
        if(.not.box%contains_valid(root_index,1).or..not.box%contains_domain(root_index,2)) cycle
      endif
      do m=1,count
        cell=root
        if (allocated(groups(root)%member)) cell=groups(root)%member(m)
        call decode(cell,i,j,k)
        candidate(:,i,j,k)=candidate(:,i,j,k)+average*recipient_weight(root,cell,m)
      end do
    end do
    if(order>0) then
      call state_redistribution_slopes_3d(geometry,qhat,centroid,order,slopes,valid, &
        boundary_state,boundary_geometry,boundary_faces,box)
      if(.not.valid) return
      do root=1,ncell
        if(kappa(root)<=0) cycle
        call decode(root,i,j,k);root_index=[i,j,k]
        if(present(box)) then
          if(.not.box%contains_valid(root_index,1).or..not.box%contains_domain(root_index,2)) cycle
        endif
        count=1
        if(allocated(groups(root)%member)) count=size(groups(root)%member)
        do m=1,count
          cell=root
          if(allocated(groups(root)%member)) cell=groups(root)%member(m)
          call decode(cell,i,j,k)
          index=0
          if(allocated(groups(root)%member)) index=groups(root)%image(:,m)
          positions(:,m)=real(index,dp)+state_cell_centroid(geometry,[i,j,k])
        enddo
        call limit_state_recipient_slopes(qhat(:,root_index(1),root_index(2),root_index(3)), &
          centroid(:,root_index(1),root_index(2),root_index(3)), &
          slopes(:,:,root_index(1),root_index(2),root_index(3)),positions(:,:count),low,high)
        do m=1,count
          cell=root
          if(allocated(groups(root)%member)) cell=groups(root)%member(m)
          call decode(cell,i,j,k)
          delta=positions(:,m)-centroid(:,root_index(1),root_index(2),root_index(3))
          candidate(:,i,j,k)=candidate(:,i,j,k)+(qhat(:,root_index(1),root_index(2),root_index(3))+ &
            matmul(delta,slopes(:,:,root_index(1),root_index(2),root_index(3))))*recipient_weight(root,cell,m)
        enddo
      enddo
    endif
    if(present(box)) call preserve_state_redist_patch_halo_3d(box,state,candidate)
    if (any(.not.ieee_is_finite(candidate))) return
    if(present(registers)) then
      do root=1,ncell
        if(kappa(root)<=0) cycle
        call decode(root,i,j,k)
        if(present(box)) then
          if(.not.box%contains_valid([i,j,k],1).or..not.box%contains_domain([i,j,k],2)) cycle
        endif
        count=1
        if(allocated(groups(root)%member)) count=size(groups(root)%member)
        volume=0;register_delta=0;register_slope=0
        if(order>0) then
          call decode(root,i,j,k);register_slope=slopes(:,:,i,j,k)
        endif
        do m=1,count
          cell=root
          if(allocated(groups(root)%member)) cell=groups(root)%member(m)
          call decode(cell,i,j,k);register_cells(:,m)=[i,j,k]
          register_out(m)=recipient_weight(root,cell,m)
          register_in(m)=kappa(cell)*register_out(m);volume=volume+register_in(m)
        enddo
        register_in(:count)=register_in(:count)/volume
        if(order>0) then
          call decode(root,i,j,k);root_index=[i,j,k]
          do m=1,count
            index=0
            if(allocated(groups(root)%member)) index=groups(root)%image(:,m)
            register_delta(:,m)=real(index,dp)+state_cell_centroid(geometry,register_cells(:,m))- &
              centroid(:,root_index(1),root_index(2),root_index(3))
          enddo
        endif
        call register_state_neighborhood(geometry,state,register_cells(:,:count),register_in(:count), &
          register_out(:count),register_delta(:,:count),register_slope,pending_registers,valid)
        if(.not.valid) return
      enddo
      registers=pending_registers
    endif
    if(present(affected)) then
      affected=.false.
      do root=1,ncell
        call decode(root,i,j,k)
        if(allocated(groups(root)%member).or.memberships(root)>1) affected(i,j,k)=.true.
        if(present(box)) then
          if(.not.box%contains_valid([i,j,k],0)) affected(i,j,k)=.false.
        endif
      enddo
    endif
    result=candidate
    ok=.true.
  contains
    real(dp) function recipient_weight(root,cell,member) result(weight)
      integer, intent(in) :: root,cell,member
      if(member==1) then
        weight=alpha_self(root)
      else
        weight=alpha_neighbor(root)/real(memberships(cell),dp)
      endif
    end function recipient_weight
    subroutine decode(id,i,j,k)
      integer, intent(in) :: id
      integer, intent(out) :: i,j,k
      i=mod(id-1,nx)+1
      j=mod((id-1)/nx,ny)+1
      k=(id-1)/(nx*ny)+1
    end subroutine decode
  end subroutine redistribute_connected_neighborhoods

  subroutine select_outward_state_receiver( &
      geometry,i,j,k,target,neighbor_i,neighbor_j,neighbor_k,ok,periodic,image_offset,box)
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: i,j,k
    real(dp), intent(in) :: target
    integer, intent(out) :: neighbor_i,neighbor_j,neighbor_k
    logical, intent(out) :: ok
    logical, intent(in) :: periodic(3)
    type(state_redist_box_3d), optional, intent(in) :: box
    integer, intent(out) :: image_offset(3)
    real(dp) :: normal(3),aperture,score,best
    integer :: axis,offset,index(3),dims(3)
    neighbor_i=i;neighbor_j=j;neighbor_k=k;image_offset=0;ok=.false.;best=-1._dp
    normal=[geometry%boundary_normal_x(i,j,k),geometry%boundary_normal_y(i,j,k), &
      geometry%boundary_normal_z(i,j,k)]
    dims=[geometry%nx,geometry%ny,geometry%nz]
    ! Largest outward normal component wins; equal scores prefer x, then y, z.
    ! Axis-plane geometries keep their original single receiver and weights.
    do axis=1,3
      score=abs(normal(axis))
      if (score<=512._dp*epsilon(1._dp).or.score<=best) cycle
      offset=merge(1,-1,normal(axis)>0._dp)
      index=[i,j,k];index(axis)=index(axis)+offset
      if(periodic(axis)) index(axis)=1+modulo(index(axis)-1,dims(axis))
      if (any(index<1).or.any(index>dims).or.all(index==[i,j,k])) cycle
      if(present(box)) then
        if(.not.box%contains_domain(index,5)) cycle
      endif
      if (geometry%volume_fraction(index(1),index(2),index(3))<target) cycle
      select case(axis)
      case(1)
        aperture=geometry%x_face_fraction(i+min(offset,0),j,k)
      case(2)
        aperture=geometry%y_face_fraction(i,j+min(offset,0),k)
      case(3)
        aperture=geometry%z_face_fraction(i,j,k+min(offset,0))
      end select
      if (aperture<=0._dp) cycle
      neighbor_i=index(1);neighbor_j=index(2);neighbor_k=index(3)
      image_offset=0;image_offset(axis)=offset
      best=score;ok=.true.
    end do
  end subroutine select_outward_state_receiver

  subroutine advance_reactive_eb_state_redistributed_3d( &
      species, state, temperature, geometry, conservative_rhs, dt, &
      new_state, new_temperature, ok, target_volume_fraction,max_order, &
      boundary_state,boundary_geometry,boundary_faces,registers)
    type(state_redist_register_3d), optional, intent(inout) :: registers
    type(state_redist_register_3d), allocatable :: pending_registers
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: conservative_rhs(:, :, :, :), dt
    real(dp), intent(out) :: new_state(:, :, :, :)
    real(dp), intent(out) :: new_temperature(:, :, :)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: boundary_state(:,:,:,:)
    type(eb_geometry_3d), optional, intent(in) :: boundary_geometry
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), intent(in), optional :: target_volume_fraction
    integer, optional, intent(in) :: max_order

    real(dp), allocatable :: provisional_state(:, :, :, :)
    real(dp), allocatable :: redistributed_state(:, :, :, :)
    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    real(dp), allocatable :: primitive(:)
    real(dp) :: recovered_temperature, sound_speed, target
    logical :: local_ok
    integer :: i, j, k, nvar

    new_state = 0.0_dp
    new_temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        size(state, 4) /= geometry%nz .or. &
        any(shape(temperature) /= &
          [geometry%nx, geometry%ny, geometry%nz]) .or. &
        any(shape(conservative_rhs) /= shape(state)) .or. &
        any(shape(new_state) /= shape(state)) .or. &
        any(shape(new_temperature) /= shape(temperature))) return
    new_state = state
    new_temperature = temperature
    if (.not. ieee_is_finite(dt) .or. dt < 0.0_dp .or. &
        any(.not. ieee_is_finite(state)) .or. &
        any(.not. ieee_is_finite(temperature)) .or. &
        any(.not. ieee_is_finite(conservative_rhs))) return

    allocate(provisional_state, source=state)
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          provisional_state(:, i, j, k) = state(:, i, j, k) + &
            dt * conservative_rhs(:, i, j, k)
        end do
      end do
    end do
    target = 0.5_dp
    if (present(target_volume_fraction)) target = target_volume_fraction
    if(present(registers)) allocate(pending_registers,source=registers)
    allocate(redistributed_state(nvar, geometry%nx, geometry%ny, geometry%nz))
    call reactive_eb_weighted_state_redistribute_3d( &
      geometry, provisional_state, redistributed_state, local_ok, target,max_order, &
      boundary_state,boundary_geometry,boundary_faces,pending_registers)
    if (.not. local_ok) return

    allocate(candidate_state, source=state)
    allocate(candidate_temperature, source=temperature)
    allocate(primitive(reactive_nprim(size(species))))
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          if (temperature(i, j, k) <= 0.0_dp) return
          candidate_state(:, i, j, k) = &
            redistributed_state(:, i, j, k)
          call reactive_conserved_to_primitive( &
            species, candidate_state(:, i, j, k), temperature(i, j, k), &
            primitive, recovered_temperature, sound_speed, local_ok)
          if (.not. local_ok) return
          candidate_temperature(i, j, k) = recovered_temperature
        end do
      end do
    end do
    if (any(.not. ieee_is_finite(candidate_state)) .or. &
        any(.not. ieee_is_finite(candidate_temperature))) return

    if(present(registers)) registers=pending_registers
    new_state = candidate_state
    new_temperature = candidate_temperature
    ok = .true.
  end subroutine advance_reactive_eb_state_redistributed_3d

end module eb_reactive_redistribution_3d_mod
