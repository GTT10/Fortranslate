module eb_reactive_patch_tree_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use rotating_frame_3d_mod, only: rotating_frame_3d,convert_rotating_energy_3d, &
    rotating_potential_3d,rotating_momentum_source_3d
  use precision_mod, only: dp
  use constant_volume_reactor_mod, only: forced_reactor_integrator
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use simulation_config_reactive_3d_mod, only: mach_inlet_ghost_3d,valid_reactive_boundaries_3d
  use reactive_3d_mod, only: advance_reactive_chemistry_3d,react_sdc_state_3d
  use gas_transport_mod, only: gas_transport_species
  use eb_reactive_transport_3d_mod, only: advance_reactive_eb_transport_patch_euler_3d, &
    reactive_eb_isothermal_timestep_3d,compute_sdc_eb_transport_3d, &
    reactive_eb_transport_patch_timestep_3d,reactive_eb_halo_bounds_3d
  use reactive_eb_cfl_3d_mod, only: reactive_eb_signal_speeds_3d,compute_reactive_eb_euler_cfl_timestep_3d
  use eb_reactive_hydro_3d_mod, only: advance_reactive_eb_euler_3d,advance_reactive_eb_plm_ctu_3d, &
    compute_reactive_eb_flattening_3d,eb_plm_boundary_3d,compute_sdc_eb_hydro_3d
  use eb_wall_heat_source_3d_mod, only: compute_eb_wall_transport_source_3d
  use eb_mol_hydro_3d_mod, only: compute_mol_eb_hydro_3d
  use eb_reactive_wall_flux_3d_mod, only: reactive_eb_slip_wall_source_3d
  use elementary_kinetics_mod, only: elementary_reaction
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use amr_patch_tree_tagging_3d_mod, only: hot_cell,build_tagged_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d,extract_eb_geometry_3d,periodic_eb_faces_match
  use eb_reactive_two_level_3d_mod, only: eb_reactive_patch_3d,matching_eb_metrics_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d,valid_eb_integrated_flux_3d, &
    restrict_eb_integrated_flux_3d,reconcile_reactive_eb_peer_flux_3d,reflux_reactive_eb_state_3d
  use eb_reactive_transfer_3d_mod, only: recover_reactive_eb_cell_3d, prolong_reactive_eb_state_pcm_3d, &
      restrict_reactive_eb_state_3d, &
    overlay_reactive_eb_peer_halo_3d
  use amr_patch_tree_remap_3d_mod, only: covered
  implicit none
  private
  type, public :: eb_reactive_tree_3d
    type(amr_patch_tree_3d) :: topology
    type(eb_reactive_patch_3d), allocatable :: nodes(:)
  end type
  public :: redistribute_eb_tree_state_3d,redistribute_sdc_eb_tree_state_rate_3d
  public :: redistribute_sdc_eb_tree_state_level_rate_3d
  public :: build_eb_tree_extended_geometry_3d
  public :: tag_reactive_eb_tree_temperature_3d
  public :: remap_reactive_eb_tree_3d
  public :: redistribute_sdc_eb_tree_flux_rate_3d,redistribute_sdc_eb_tree_level_rate_3d
  public :: reactive_eb_tree_timestep_3d
  public :: advance_reactive_eb_tree_full_3d,advance_reactive_eb_tree_sdc_3d
  public :: advance_reactive_eb_tree_sdc_subcycle_3d
  public :: advance_reactive_eb_tree_transport_3d,compute_sdc_eb_tree_transport_3d
  public :: advance_reactive_eb_tree_subcycle_transport_3d
  public :: advance_reactive_eb_tree_subcycle_hydro_3d
  public :: advance_reactive_eb_tree_hydro_3d,compute_sdc_eb_tree_hydro_3d
  public :: compute_eb_tree_wall_transport_3d
  public :: compute_mol_eb_tree_hydro_3d
  public :: synchronize_reactive_eb_tree_flux_3d
  public :: fill_reactive_eb_tree_halo_3d,fill_sdc_eb_tree_source_halo_3d
  public :: react_reactive_eb_tree_3d,react_sdc_eb_tree_3d
  public :: initialize_reactive_eb_tree_3d,average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d
  ! Auxiliary sources receive the old state at stage 0 and provisional new
  ! states at stages 1..iterations. Rates are per fluid volume per time.
  ! Callback-owned physical state must remain private until the driver succeeds.
  abstract interface
    subroutine eb_sdc_auxiliary_source(state,interval,stage,depth,rate,ok)
      import :: eb_reactive_tree_3d,amr_tree_state_3d,dp
      type(eb_reactive_tree_3d), intent(in) :: state
      real(dp), intent(in) :: interval
      integer, intent(in) :: stage,depth
      type(amr_tree_state_3d), intent(inout) :: rate
      logical, intent(out) :: ok
    end subroutine
    ! Time is the level interval start relative to the outer step, for all stages.
    subroutine eb_sdc_recursive_source(state,interval,stage,depth,time,rate,ok)
      import :: eb_reactive_tree_3d,amr_tree_state_3d,dp
      type(eb_reactive_tree_3d), intent(in) :: state
      real(dp), intent(in) :: interval,time
      integer, intent(in) :: stage,depth
      type(amr_tree_state_3d), intent(inout) :: rate
      logical, intent(out) :: ok
    end subroutine
  end interface
  public :: eb_sdc_auxiliary_source,eb_sdc_recursive_source
contains
  ! Composite flux redistribution. Covered parents do not donate or receive.
  ! Each cut leaf exchanges excess with face-connected fluid leaves, including
  ! peers and coarse/fine neighbors. Paired periodic domain faces wrap; other
  ! physical faces are closed for redistribution.
  subroutine redistribute_sdc_eb_tree_flux_rate_3d(tree,rhs,rate,ok,boundary_faces)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: rhs
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    character(len=32) :: faces(6)
    type(amr_tree_state_3d) :: candidate,active,delta
    integer(int64), allocatable :: origin(:,:),domain(:,:),scale(:)
    integer, allocatable :: ratios(:)
    integer(int64) :: cell_lo(3),cell_hi(3)
    integer :: nv,p,child,n(3),lo(3),hi(3),i,j,k,pass
    real(dp) :: kappa,volume,neighbors,total
    real(dp), allocatable :: average(:),excess(:)
    logical :: valid
    ok=.false.
    if(.not.allocated(rhs%fields)) return
    if(size(rhs%fields)<1) return
    if(.not.allocated(rhs%fields(1)%values)) return
    nv=size(rhs%fields(1)%values,1)
    if(nv<1) return
    if(.not.valid_tree(tree,nv)) return
    faces='outflow'
    if(present(boundary_faces)) faces=boundary_faces
    if(.not.valid_reactive_boundaries_3d(faces)) return
    if(.not.periodic_eb_faces_match(tree%nodes(1)%geometry,faces)) return
    if(size(rhs%fields)/=size(tree%nodes)) return
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if(.not.valid) return
    allocate(average(nv),excess(nv))
    candidate=rhs
    allocate(active%fields(size(tree%nodes)),delta%fields(size(tree%nodes)),scale(size(ratios)))
    scale=domain(1,size(ratios))/domain(1,:)
    do p=1,size(tree%nodes)
      n=tree%topology%nodes(p)%shape
      if(.not.allocated(rhs%fields(p)%values)) return
      if(any(shape(rhs%fields(p)%values)/=[nv,n])) return
      allocate(active%fields(p)%values(1,n(1),n(2),n(3)),source=1._dp)
      allocate(delta%fields(p)%values(nv,n(1),n(2),n(3)),source=0._dp)
      do child=p+1,size(tree%nodes)
        if(tree%topology%nodes(child)%parent/=p) cycle
        call bounds(tree%topology,child,lo,hi)
        active%fields(p)%values(1,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=0
      end do
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(tree%nodes(p)%geometry%volume_fraction(i,j,k)<=0) then
          candidate%fields(p)%values(:,i,j,k)=0
          active%fields(p)%values(1,i,j,k)=0
        else
          if(.not.all(ieee_is_finite(rhs%fields(p)%values(:,i,j,k)))) return
        end if
      end do
      end do
      end do
    end do
    do p=1,size(tree%nodes)
      n=tree%topology%nodes(p)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(active%fields(p)%values(1,i,j,k)==0) cycle
        kappa=tree%nodes(p)%geometry%volume_fraction(i,j,k)
        if(kappa>=1) cycle
        cell_lo=(origin(:,p)+[i-1,j-1,k-1])*scale(tree%topology%nodes(p)%level+1)
        cell_hi=cell_lo+scale(tree%topology%nodes(p)%level+1)
        volume=kappa*tree%nodes(p)%geometry%dx*tree%nodes(p)%geometry%dy*tree%nodes(p)%geometry%dz
        if(.not.ieee_is_finite(volume).or.volume<=0) return
        average=volume*rhs%fields(p)%values(:,i,j,k);neighbors=0
        do pass=1,2
          call visit_neighbors()
          if(pass==1) then
            total=volume+neighbors
            if(neighbors<=0.or..not.ieee_is_finite(total)) return
            average=average/total
            candidate%fields(p)%values(:,i,j,k)=kappa*rhs%fields(p)%values(:,i,j,k)+(1-kappa)*average
            excess=volume*(1-kappa)*(rhs%fields(p)%values(:,i,j,k)-average)/neighbors
          end if
        end do
      end do
      end do
      end do
    end do
    do p=1,size(tree%nodes)
      candidate%fields(p)%values=candidate%fields(p)%values+delta%fields(p)%values
      if(.not.all(ieee_is_finite(candidate%fields(p)%values))) return
    end do
    call move_alloc(candidate%fields,rate%fields)
    ok=.true.
  contains
    subroutine visit_neighbors()
      integer :: axis,side,q,direction,lower(3),upper(3),xyz(3),a,b,c,face_cell
      integer(int64) :: face,s,qlo(3),qhi(3),first(3),last(3)
      real(dp) :: neighbor_volume
      do axis=1,3
      do side=0,1
        if(.not.open_face(tree%nodes(p)%geometry,[i,j,k],axis,side)) cycle
        face=cell_lo(axis)
        if(side==1) face=cell_hi(axis)
        if(faces(2*axis-1)=='periodic') then
          if(side==0.and.face==0) face=domain(axis,size(ratios))
          if(side==1.and.face==domain(axis,size(ratios))) face=0
        end if
        do q=1,size(tree%nodes)
          s=scale(tree%topology%nodes(q)%level+1)
          ! A leaf can touch this face only on an integer cell boundary.
          if(modulo(face,s)/=0) cycle
          qlo=origin(:,q)*s;qhi=(origin(:,q)+tree%topology%nodes(q)%shape)*s
          if(side==0.and.(face<=qlo(axis).or.face>qhi(axis))) cycle
          if(side==1.and.(face<qlo(axis).or.face>=qhi(axis))) cycle
          first=0;last=0
          do direction=1,3
            if(direction==axis) cycle
            first(direction)=max(cell_lo(direction),qlo(direction))
            last(direction)=min(cell_hi(direction),qhi(direction))
          end do
          if(any([(last(direction)<=first(direction).and.direction/=axis,direction=1,3)])) cycle
          do direction=1,3
            if(direction==axis) cycle
            lower(direction)=int(first(direction)/s-origin(direction,q))+1
            upper(direction)=int((last(direction)-1)/s-origin(direction,q))+1
          end do
          face_cell=int(face/s-origin(axis,q))+side
          lower(axis)=face_cell;upper(axis)=face_cell
          do c=lower(3),upper(3)
          do b=lower(2),upper(2)
          do a=lower(1),upper(1)
            xyz=[a,b,c]
            if(active%fields(q)%values(1,a,b,c)==0) cycle
            if(.not.open_face(tree%nodes(q)%geometry,xyz,axis,1-side)) cycle
            neighbor_volume=tree%nodes(q)%geometry%volume_fraction(a,b,c)* &
              tree%nodes(q)%geometry%dx*tree%nodes(q)%geometry%dy*tree%nodes(q)%geometry%dz
            if(pass==1) then
              average=average+neighbor_volume*rhs%fields(q)%values(:,a,b,c)
              neighbors=neighbors+neighbor_volume
            else
              delta%fields(q)%values(:,a,b,c)=delta%fields(q)%values(:,a,b,c)+excess
            end if
          end do
          end do
          end do
        end do
      end do
      end do
    end subroutine visit_neighbors
    logical function open_face(g,cell,axis,side) result(opened)
      type(eb_geometry_3d), intent(in) :: g
      integer, intent(in) :: cell(3),axis,side
      integer :: f(3)
      f=cell;f(axis)=f(axis)-1+side
      select case(axis)
      case(1)
        opened=g%x_face_fraction(f(1),f(2),f(3))>0
      case(2)
        opened=g%y_face_fraction(f(1),f(2),f(3))>0
      case(3)
        opened=g%z_face_fraction(f(1),f(2),f(3))>0
      end select
    end function open_face
  end subroutine redistribute_sdc_eb_tree_flux_rate_3d

  ! Decompose the linear composite operator by source level. The selected
  ! level receives a rate; all other levels receive a time-integrated transfer
  ! for a caller-owned interval register. No gas state is updated here.
  subroutine redistribute_sdc_eb_tree_level_rate_3d(tree,rhs,depth,dt,rate,transfer,ok,boundary_faces)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: rhs
    integer, intent(in) :: depth
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: rate,transfer
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(amr_tree_state_3d) :: selected,result,exchange
    integer :: p,nv
    logical :: good
    ok=.false.
    if(.not.tree%topology%is_valid()) return
    if(depth<0.or.depth>=tree%topology%level_count()) return
    if(.not.ieee_is_finite(dt)) return
    if(dt<=0.or..not.allocated(rhs%fields)) return
    if(size(rhs%fields)/=size(tree%topology%nodes)) return
    if(.not.allocated(rhs%fields(1)%values)) return
    nv=size(rhs%fields(1)%values,1)
    selected=rhs
    do p=1,size(tree%topology%nodes)
      if(.not.allocated(rhs%fields(p)%values)) return
      if(any(shape(rhs%fields(p)%values)/=[nv,tree%topology%nodes(p)%shape])) return
      if(tree%topology%nodes(p)%level/=depth) selected%fields(p)%values=0
    end do
    call redistribute_sdc_eb_tree_flux_rate_3d(tree,selected,result,good,boundary_faces)
    if(.not.good) return
    exchange=result
    do p=1,size(tree%topology%nodes)
      if(tree%topology%nodes(p)%level==depth) then
        exchange%fields(p)%values=0
      else
        if(dt>1) then
          if(any(abs(exchange%fields(p)%values)>huge(1._dp)/dt)) return
        end if
        exchange%fields(p)%values=dt*exchange%fields(p)%values
        result%fields(p)%values=0
      end if
    end do
    call move_alloc(result%fields,rate%fields)
    call move_alloc(exchange%fields,transfer%fields)
    ok=.true.
  end subroutine redistribute_sdc_eb_tree_level_rate_3d

  ! S(P_l R) + D_l(U)/dt, where D_l contains conservative neighborhood
  ! corrections rooted on level l. Off-level increments are returned separately.
  subroutine redistribute_sdc_eb_tree_state_level_rate_3d(tree,rhs,depth,dt,rate,transfer,ok,target_volume_fraction,boundary_faces)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: rhs
    integer, intent(in) :: depth
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: rate,transfer
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: target_volume_fraction
    type(amr_tree_state_3d) :: selected,result,exchange,baseline,correction
    integer :: p,nv
    logical :: good
    ok=.false.
    if(.not.tree%topology%is_valid()) return
    if(depth<0.or.depth>=tree%topology%level_count()) return
    if(.not.ieee_is_finite(dt)) return
    if(dt<=0.or..not.allocated(rhs%fields)) return
    if(size(rhs%fields)/=size(tree%topology%nodes)) return
    if(.not.allocated(rhs%fields(1)%values)) return
    nv=size(rhs%fields(1)%values,1)
    selected=rhs
    do p=1,size(tree%topology%nodes)
      if(.not.allocated(rhs%fields(p)%values)) return
      if(any(shape(rhs%fields(p)%values)/=[nv,tree%topology%nodes(p)%shape])) return
      if(tree%topology%nodes(p)%level/=depth) selected%fields(p)%values=0
    end do
    call redistribute_eb_tree_state_3d(tree,selected,result,good,target_volume_fraction,boundary_faces)
    if(.not.good) return
    baseline=selected
    do p=1,size(tree%nodes)
      baseline%fields(p)%values=tree%nodes(p)%state
    end do
    call redistribute_eb_tree_state_3d(tree,baseline,correction,good,target_volume_fraction,boundary_faces,depth)
    if(.not.good) return
    exchange=result
    do p=1,size(tree%topology%nodes)
      if(tree%topology%nodes(p)%level==depth) then
        exchange%fields(p)%values=0
        result%fields(p)%values=result%fields(p)%values+correction%fields(p)%values/dt
      else
        if(dt>1) then
          if(any(abs(exchange%fields(p)%values)>huge(1._dp)/dt)) return
        end if
        exchange%fields(p)%values=dt*exchange%fields(p)%values+correction%fields(p)%values
        result%fields(p)%values=0
      end if
    end do
    do p=1,size(tree%nodes)
      if(.not.all(ieee_is_finite(result%fields(p)%values))) return
      if(.not.all(ieee_is_finite(exchange%fields(p)%values))) return
    end do
    call move_alloc(result%fields,rate%fields)
    call move_alloc(exchange%fields,transfer%fields)
    ok=.true.
  end subroutine redistribute_sdc_eb_tree_state_level_rate_3d

  ! Zeroth-order weighted StateRedist on composite fluid leaves. Large outward
  ! receivers use the single-grid weights; otherwise connected neighborhoods
  ! grow to the source cell's target physical volume with overlap accounting.
  ! correction_depth returns only the conservative displacement contributed by
  ! neighborhoods rooted at that depth, with zero hidden/solid entries.
  subroutine redistribute_eb_tree_state_3d(tree,input,result,ok,target_volume_fraction,boundary_faces,correction_depth)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: input
    type(amr_tree_state_3d), intent(inout) :: result
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: target_volume_fraction
    character(len=*), optional, intent(in) :: boundary_faces(6)
    integer, optional, intent(in) :: correction_depth
    type :: neighborhood
      integer, allocatable :: member(:)
    end type
    type(neighborhood), allocatable :: groups(:)
    type(amr_tree_state_3d) :: ids,candidate
    integer(int64), allocatable :: origin(:,:),domain(:,:),scale(:)
    integer, allocatable :: ratios(:),cells(:,:),receiver(:),memberships(:),queue(:),seen(:),neighbors(:),directions(:)
    real(dp), allocatable :: volume(:),target(:),beta(:),self(:),value(:,:),average(:,:),output(:,:)
    integer(int64) :: maximum,cell_lo(3),cell_hi(3)
    character(len=32) :: faces(6)
    integer :: nv,p,n(3),i,j,k,child,lo(3),hi(3),nc,root,other,m,count,head,tail,current,allocation_status
    real(dp) :: fraction,total,weight,normal(3),score,best
    logical :: valid,fallback
    ok=.false.;fraction=.5_dp;faces='outflow'
    if(present(target_volume_fraction)) fraction=target_volume_fraction
    if(present(boundary_faces)) faces=boundary_faces
    if(.not.ieee_is_finite(fraction)) return
    if(fraction<=0.or.fraction>1.or..not.valid_reactive_boundaries_3d(faces)) return
    if(.not.allocated(input%fields)) return
    if(size(input%fields)<1) return
    if(.not.allocated(input%fields(1)%values)) return
    nv=size(input%fields(1)%values,1)
    if(nv<1.or..not.valid_tree(tree,nv)) return
    if(size(input%fields)/=size(tree%nodes)) return
    if(present(correction_depth)) then
      if(correction_depth<0.or.correction_depth>=tree%topology%level_count()) return
    end if
    if(.not.periodic_eb_faces_match(tree%nodes(1)%geometry,faces)) return
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if(.not.valid) return
    maximum=0
    do p=1,size(tree%nodes)
      if(.not.allocated(input%fields(p)%values)) return
      n=tree%topology%nodes(p)%shape
      if(any(shape(input%fields(p)%values)/=[nv,n])) return
      if(int(n(1),int64)*n(2)>int(huge(1),int64)/n(3)) return
      maximum=maximum+int(n(1),int64)*n(2)*n(3)
      if(maximum>huge(1)) return
    end do
    allocate(ids%fields(size(tree%nodes)),scale(size(ratios)),cells(4,int(maximum)),stat=allocation_status)
    if(allocation_status/=0) return
    scale=domain(1,size(ratios))/domain(1,:)
    candidate=input;nc=0
    if(present(correction_depth)) then
      do p=1,size(tree%nodes)
        candidate%fields(p)%values=0
      end do
    end if
    do p=1,size(tree%nodes)
      n=tree%topology%nodes(p)%shape
      allocate(ids%fields(p)%values(1,n(1),n(2),n(3)),source=0._dp,stat=allocation_status)
      if(allocation_status/=0) return
      do child=p+1,size(tree%nodes)
        if(tree%topology%nodes(child)%parent/=p) cycle
        call bounds(tree%topology,child,lo,hi)
        ids%fields(p)%values(1,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=-1
      end do
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(tree%nodes(p)%geometry%volume_fraction(i,j,k)<=0) then
          candidate%fields(p)%values(:,i,j,k)=0
          ids%fields(p)%values(1,i,j,k)=0
          cycle
        end if
        if(.not.all(ieee_is_finite(input%fields(p)%values(:,i,j,k)))) return
        if(ids%fields(p)%values(1,i,j,k)<0) cycle
        nc=nc+1;cells(:,nc)=[p,i,j,k];ids%fields(p)%values(1,i,j,k)=real(nc,dp)
      end do
      end do
      end do
    end do
    if(nc>huge(1)/6) return
    allocate(groups(nc),receiver(nc),memberships(nc),queue(nc),seen(nc),neighbors(6*nc),directions(6*nc), &
      volume(nc),target(nc),beta(nc),self(nc),value(nv,nc),average(nv,nc),output(nv,nc),stat=allocation_status)
    if(allocation_status/=0) return
    do root=1,nc
      p=cells(1,root);i=cells(2,root);j=cells(3,root);k=cells(4,root)
      associate(g=>tree%nodes(p)%geometry)
        target(root)=fraction*g%dx*g%dy*g%dz
        volume(root)=g%volume_fraction(i,j,k)*g%dx*g%dy*g%dz
      end associate
      if(min(volume(root),target(root))<=0) return
      value(:,root)=input%fields(p)%values(:,i,j,k)
    end do
    receiver=0;memberships=1;self=1;beta=0;fallback=.false.
    do root=1,nc
      if(volume(root)>=target(root)) cycle
      current=root
      call adjacent(count)
      p=cells(1,root);i=cells(2,root);j=cells(3,root);k=cells(4,root)
      associate(g=>tree%nodes(p)%geometry)
        normal=[g%boundary_normal_x(i,j,k),g%boundary_normal_y(i,j,k),g%boundary_normal_z(i,j,k)]
      end associate
      best=-1
      do m=1,count
        other=neighbors(m)
        if(volume(other)<max(target(root),target(other))) cycle
        score=normal(abs(directions(m)))*real(sign(1,directions(m)),dp)
        if(score<=512*epsilon(1._dp).or.score<=best) cycle
        best=score;receiver(root)=other
      end do
      if(receiver(root)==0) then
        fallback=.true.;exit
      end if
      other=receiver(root)
      beta(root)=(target(root)-volume(root))/volume(other)
      memberships(other)=memberships(other)+1
    end do
    if(.not.fallback) then
      do root=1,nc
        other=receiver(root)
        if(other==0) cycle
        self(other)=self(other)-beta(root)/real(memberships(other),dp)
      end do
      do root=1,nc
        total=self(root)*volume(root);average(:,root)=total*value(:,root)
        other=receiver(root)
        if(other>0) then
          weight=beta(root)*volume(other)/real(memberships(other),dp)
          total=total+weight;average(:,root)=average(:,root)+weight*value(:,other)
        end if
        if(total<=0.or.self(root)<0) return
        average(:,root)=average(:,root)/total
      end do
      output=0
      do root=1,nc
        if(present(correction_depth)) then
          if(tree%topology%nodes(cells(1,root))%level/=correction_depth) cycle
          output(:,root)=output(:,root)+self(root)*real(memberships(root),dp)*(average(:,root)-value(:,root))
          other=receiver(root)
          if(other>0) output(:,other)=output(:,other)+beta(root)*(average(:,root)-value(:,other))
        else
          output(:,root)=output(:,root)+self(root)*real(memberships(root),dp)*average(:,root)
          other=receiver(root)
          if(other>0) output(:,other)=output(:,other)+beta(root)*average(:,root)
        end if
      end do
      do root=1,nc
        output(:,root)=output(:,root)/real(memberships(root),dp)
      end do
    else
      memberships=1;seen=0
      do root=1,nc
        if(volume(root)>=target(root)) cycle
        head=1;tail=1;queue(1)=root;seen(root)=root;total=volume(root)
        do while(total<target(root).and.head<=tail)
          current=queue(head)
          call adjacent(count)
          do m=1,count
            if(total>=target(root)) exit
            other=neighbors(m)
            if(seen(other)==root) cycle
            tail=tail+1;queue(tail)=other;seen(other)=root;total=total+volume(other)
          end do
          head=head+1
        end do
        if(total<target(root)) return
        allocate(groups(root)%member(tail),stat=allocation_status)
        if(allocation_status/=0) return
        groups(root)%member=queue(:tail)
        do m=2,tail
          other=queue(m);memberships(other)=memberships(other)+1
        end do
      end do
      output=0
      do root=1,nc
        count=1
        if(allocated(groups(root)%member)) count=size(groups(root)%member)
        total=0;average(:,root)=0
        do m=1,count
          other=root
          if(allocated(groups(root)%member)) other=groups(root)%member(m)
          weight=volume(other)/real(memberships(other),dp)
          total=total+weight;average(:,root)=average(:,root)+weight*value(:,other)
        end do
        if(total<=0) return
        average(:,root)=average(:,root)/total
        if(present(correction_depth)) then
          if(tree%topology%nodes(cells(1,root))%level/=correction_depth) cycle
        end if
        do m=1,count
          other=root
          if(allocated(groups(root)%member)) other=groups(root)%member(m)
          if(present(correction_depth)) then
            output(:,other)=output(:,other)+(average(:,root)-value(:,other))/real(memberships(other),dp)
          else
            output(:,other)=output(:,other)+average(:,root)/real(memberships(other),dp)
          end if
        end do
      end do
    end if
    if(.not.all(ieee_is_finite(output))) return
    do root=1,nc
      p=cells(1,root);i=cells(2,root);j=cells(3,root);k=cells(4,root)
      candidate%fields(p)%values(:,i,j,k)=output(:,root)
    end do
    call move_alloc(candidate%fields,result%fields);ok=.true.
  contains
    subroutine adjacent(count)
      integer, intent(out) :: count
      integer :: axis,side,q,direction,lower(3),upper(3),xyz(3),a,b,c,face_cell
      integer(int64) :: face,s,qlo(3),qhi(3),first(3),last(3)
      integer :: neighbor
      p=cells(1,current);i=cells(2,current);j=cells(3,current);k=cells(4,current)
      cell_lo=(origin(:,p)+[i-1,j-1,k-1])*scale(tree%topology%nodes(p)%level+1)
      cell_hi=cell_lo+scale(tree%topology%nodes(p)%level+1)
      count=0
      do axis=1,3
      do side=0,1
        if(.not.open_face(tree%nodes(p)%geometry,[i,j,k],axis,side)) cycle
        face=cell_lo(axis)
        if(side==1) face=cell_hi(axis)
        if(faces(2*axis-1)=='periodic') then
          if(side==0.and.face==0) face=domain(axis,size(ratios))
          if(side==1.and.face==domain(axis,size(ratios))) face=0
        end if
        do q=1,size(tree%nodes)
          s=scale(tree%topology%nodes(q)%level+1)
          ! A leaf can touch this face only on an integer cell boundary.
          if(modulo(face,s)/=0) cycle
          qlo=origin(:,q)*s;qhi=(origin(:,q)+tree%topology%nodes(q)%shape)*s
          if(side==0.and.(face<=qlo(axis).or.face>qhi(axis))) cycle
          if(side==1.and.(face<qlo(axis).or.face>=qhi(axis))) cycle
          first=0;last=0
          do direction=1,3
            if(direction==axis) cycle
            first(direction)=max(cell_lo(direction),qlo(direction))
            last(direction)=min(cell_hi(direction),qhi(direction))
          end do
          if(any([(last(direction)<=first(direction).and.direction/=axis,direction=1,3)])) cycle
          do direction=1,3
            if(direction==axis) cycle
            lower(direction)=int(first(direction)/s-origin(direction,q))+1
            upper(direction)=int((last(direction)-1)/s-origin(direction,q))+1
          end do
          face_cell=int(face/s-origin(axis,q))+side
          lower(axis)=face_cell;upper(axis)=face_cell
          do c=lower(3),upper(3)
          do b=lower(2),upper(2)
          do a=lower(1),upper(1)
            xyz=[a,b,c]
            if(ids%fields(q)%values(1,a,b,c)<=0) cycle
            if(.not.open_face(tree%nodes(q)%geometry,xyz,axis,1-side)) cycle
            neighbor=int(ids%fields(q)%values(1,a,b,c))
            count=count+1;neighbors(count)=neighbor;directions(count)=axis*(2*side-1)
          end do
          end do
          end do
        end do
      end do
      end do
    end subroutine adjacent
    logical function open_face(g,cell,axis,side) result(opened)
      type(eb_geometry_3d), intent(in) :: g
      integer, intent(in) :: cell(3),axis,side
      integer :: f(3)
      f=cell;f(axis)=f(axis)-1+side
      select case(axis)
      case(1)
        opened=g%x_face_fraction(f(1),f(2),f(3))>0
      case(2)
        opened=g%y_face_fraction(f(1),f(2),f(3))>0
      case(3)
        opened=g%z_face_fraction(f(1),f(2),f(3))>0
      end select
    end function open_face
  end subroutine redistribute_eb_tree_state_3d

  ! Match the single-grid StateRedist source convention: S(U+dt*R)-U,
  ! divided by dt. Covered cells are excluded before arithmetic.
  subroutine redistribute_sdc_eb_tree_state_rate_3d(tree,rhs,dt,rate,ok,target_volume_fraction,boundary_faces)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: rhs
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: target_volume_fraction
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(amr_tree_state_3d) :: provisional,candidate
    integer :: nv,p,n(3),i,j,k
    logical :: valid
    ok=.false.
    if(.not.ieee_is_finite(dt)) return
    if(dt<=0.or..not.allocated(rhs%fields)) return
    if(size(rhs%fields)<1) return
    if(.not.allocated(rhs%fields(1)%values)) return
    nv=size(rhs%fields(1)%values,1)
    if(nv<1.or..not.valid_tree(tree,nv)) return
    if(size(rhs%fields)/=size(tree%nodes)) return
    provisional=rhs
    do p=1,size(tree%nodes)
      n=tree%topology%nodes(p)%shape
      if(.not.allocated(rhs%fields(p)%values)) return
      if(any(shape(rhs%fields(p)%values)/=[nv,n])) return
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        provisional%fields(p)%values(:,i,j,k)=0
        if(tree%nodes(p)%geometry%volume_fraction(i,j,k)<=0) cycle
        if(.not.all(ieee_is_finite(rhs%fields(p)%values(:,i,j,k)))) return
        if(.not.all(ieee_is_finite(tree%nodes(p)%state(:,i,j,k)))) return
        provisional%fields(p)%values(:,i,j,k)=tree%nodes(p)%state(:,i,j,k)+dt*rhs%fields(p)%values(:,i,j,k)
      end do
      end do
      end do
    end do
    call redistribute_eb_tree_state_3d(tree,provisional,candidate,valid,target_volume_fraction,boundary_faces)
    if(.not.valid) return
    do p=1,size(tree%nodes)
      n=tree%topology%nodes(p)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(tree%nodes(p)%geometry%volume_fraction(i,j,k)<=0) cycle
        candidate%fields(p)%values(:,i,j,k)=(candidate%fields(p)%values(:,i,j,k)-tree%nodes(p)%state(:,i,j,k))/dt
      end do
      end do
      end do
      if(.not.all(ieee_is_finite(candidate%fields(p)%values))) return
    end do
    call move_alloc(candidate%fields,rate%fields);ok=.true.
  end subroutine redistribute_sdc_eb_tree_state_rate_3d

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

  subroutine initialize_reactive_eb_tree_3d(species,topology,finest,root_state,root_temperature,tree,ok,frame)
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
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
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
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
          candidate%nodes(1)%state,candidate%nodes(1)%temperature,valid,frame)
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
          candidate%nodes(node)%state,candidate%nodes(node)%temperature,valid,frame)
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
  subroutine remap_reactive_eb_tree_3d(species,topology,finest,tree,ok,reaction_source)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: topology
    type(eb_geometry_3d), intent(in) :: finest
    type(eb_reactive_tree_3d), intent(inout) :: tree
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source
    type(amr_tree_state_3d) :: old_source,new_source
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
    if(present(reaction_source)) then
      if(.not.reaction_source%is_valid(tree%topology)) return
      if(size(reaction_source%fields(1)%values,1)/=nv-4) return
      old_source=reaction_source
      call restrict_source(tree,old_source,valid)
      if(.not.valid) return
    end if
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
    if(present(reaction_source)) allocate(new_source%fields(size(candidate%nodes)))
    do node=1,size(candidate%nodes)
      if(present(reaction_source)) then
        d=topology%nodes(node)%shape
        allocate(new_source%fields(node)%values(nv-4,d(1),d(2),d(3)),source=0._dp)
      end if
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
            if(present(reaction_source)) then
              if(candidate%nodes(node)%geometry%volume_fraction(i,j,k)>0) &
                new_source%fields(node)%values(:,i,j,k)=old_source%fields(donor)%values(:,index(1),index(2),index(3))
            end if
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
    if(present(reaction_source)) then
      call restrict_source(candidate,new_source,valid)
      if(.not.valid) return
      reaction_source=new_source
    end if
    tree=candidate;ok=.true.
  contains
    subroutine restrict_source(layout,field,good)
      type(eb_reactive_tree_3d), intent(in) :: layout
      type(amr_tree_state_3d), intent(inout) :: field
      logical, intent(out) :: good
      integer :: child,parent,r,lo(3),hi(3),x,y,z,first(3),last(3),component
      real(dp) :: volume,value(nv-4)
      good=.false.
      do child=size(layout%nodes),2,-1
        parent=layout%topology%nodes(child)%parent
        r=layout%topology%nodes(child)%patch%refinement_ratio
        call bounds(layout%topology,child,lo,hi)
        do z=lo(3),hi(3)
        do y=lo(2),hi(2)
        do x=lo(1),hi(1)
          first=([x,y,z]-lo)*r+1;last=first+r-1
          associate(v=>layout%nodes(child)%geometry%volume_fraction(first(1):last(1),first(2):last(2),first(3):last(3)))
          volume=sum(v);value=0
          if(volume>0) then
            do component=1,nv-4
              value(component)=sum(v*field%fields(child)%values(component, &
                first(1):last(1),first(2):last(2),first(3):last(3)))/volume
            end do
          end if
          end associate
          if(.not.all(ieee_is_finite(value))) return
          field%fields(parent)%values(:,x,y,z)=value
        end do
        end do
        end do
      end do
      good=.true.
    end subroutine restrict_source
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

  subroutine average_down_reactive_eb_tree_3d(species,tree,ok,frame)
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(eb_reactive_tree_3d), intent(inout) :: tree
    logical, intent(out) :: ok
    type(eb_reactive_tree_3d) :: candidate
    integer :: node,parent,lo(3),hi(3)
    logical :: valid
    ok=.false.
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    if (.not.valid_tree(tree,reactive_nvar(size(species)))) return
    candidate=tree
    do node=size(tree%nodes),2,-1
      parent=tree%topology%nodes(node)%parent
      call bounds(tree%topology,node,lo,hi)
      call restrict_reactive_eb_state_3d(species,candidate%nodes(node)%geometry, &
        tree%topology%nodes(node)%patch%refinement_ratio,candidate%nodes(node)%state,candidate%nodes(node)%temperature, &
        candidate%nodes(parent)%state(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
        candidate%nodes(parent)%temperature(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),valid,frame)
      if (.not.valid) return
    end do
    tree=candidate;ok=.true.
  end subroutine average_down_reactive_eb_tree_3d

  ! Bound every patch update, including hidden parent cells which are advanced
  ! before restriction. Global donor maxima cover peer and ancestor halos.
  subroutine reactive_eb_tree_timestep_3d(species,transport,tree,hydro_cfl,transport_cfl, &
      transport_enabled,viscosity,conduction,diffusion,dt,ok,subcycling, &
      inflow_mach,inflow_temperature,inflow_pressure,boundary_faces,inflow_state,wall_temperature,coupled_subcycling)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(eb_reactive_tree_3d), intent(in) :: tree
    real(dp), intent(in) :: hydro_cfl,transport_cfl
    logical, intent(in) :: transport_enabled,viscosity,conduction,diffusion
    real(dp), intent(inout) :: dt
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: wall_temperature(6)
    character(len=16) :: thermal_faces(6)
    integer :: direction
    logical :: thermal_wall
    logical, intent(out) :: ok
    logical, optional, intent(in) :: subcycling,coupled_subcycling
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    type(eb_reactive_tree_3d) :: snapshot
    real(dp) :: speeds(4),donors(4),diffusivity,maximum_diffusivity,interval,candidate,factor
    real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:)
    integer :: node,level,n(3)
    logical :: valid,active,recursive,coupled
    ok=.false.
    if (.not.all(ieee_is_finite([hydro_cfl,transport_cfl]))) return
    if (hydro_cfl<=0._dp.or.hydro_cfl>1._dp.or.transport_cfl<=0._dp.or.transport_cfl>.5_dp) return
    if (.not.valid_tree(tree,reactive_nvar(size(species)))) return
    thermal_wall=.false.
    if(present(wall_temperature)) then
      if(.not.present(boundary_faces)) return
      if(.not.all(ieee_is_finite(wall_temperature)).or.any(wall_temperature<=0)) return
    end if
    if(present(boundary_faces)) thermal_wall=any(boundary_faces=='isothermal_wall')
    if(thermal_wall) then
      if(.not.present(wall_temperature).or..not.transport_enabled.or..not.conduction) return
    end if
    coupled=.false.
    if(present(coupled_subcycling)) coupled=coupled_subcycling
    recursive=.false.
    if (present(subcycling)) recursive=subcycling
    if (recursive.or.thermal_wall) then
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
        if(.not.coupled) factor=factor*factor
        if (interval<candidate/factor) candidate=interval*factor
        if(thermal_wall) then
          level=tree%topology%nodes(node)%level+1;n=tree%topology%nodes(node)%shape
          thermal_faces='outflow'
          do direction=1,3
            if(origin(direction,node)==0) thermal_faces(2*direction-1)=boundary_faces(2*direction-1)
            if(origin(direction,node)+n(direction)==domain(direction,level)) &
              thermal_faces(2*direction)=boundary_faces(2*direction)
          end do
          call reactive_eb_isothermal_timestep_3d(species,transport,u,t,g,thermal_faces, &
            wall_temperature,transport_cfl,interval,valid)
          if(.not.valid) return
          if(interval<candidate/factor) candidate=interval*factor
        end if
      end if
      end associate
    end do
    dt=candidate;ok=.true.
  end subroutine reactive_eb_tree_timestep_3d

  ! Recursive SDC with interval-owned redistribution registers. Each register
  ! closes when the coarser source/recipient level finishes its interval.
  subroutine advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,dt,iterations,chemistry, &
      transport_flags,force,rtol,atol,tree,extended_geometry,boundary_faces, &
      boundary_change,wall_change,chemical_change,source_change,minimum_theta,ok,reaction_source, &
      difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature,level_steps,redistribution, &
      target_volume_fraction,auxiliary_source,auxiliary_change,route_auxiliary, &
      evaluate_reaction,mol_mode,mol_lax_friedrichs,embedded_wall_temperature,embedded_noslip,frame)
    real(dp), optional, intent(in) :: embedded_wall_temperature
    logical, optional, intent(in) :: mol_mode,mol_lax_friedrichs,embedded_noslip
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: dt,force(3),rtol,atol
    integer, intent(in) :: iterations
    logical, intent(in) :: chemistry,transport_flags(4)
    type(eb_reactive_tree_3d), intent(inout) :: tree
    type(eb_geometry_3d), intent(in) :: extended_geometry(:)
    character(len=*), intent(in) :: boundary_faces(6)
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemical_change(:),source_change(:),minimum_theta
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source
    real(dp), optional, intent(in) :: difmag,inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: wall_temperature(6)
    integer(int64), optional, intent(inout) :: level_steps(:)
    character(len=*), optional, intent(in) :: redistribution
    real(dp), optional, intent(in) :: target_volume_fraction
    procedure(eb_sdc_recursive_source), optional :: auxiliary_source
    real(dp), optional, intent(inout) :: auxiliary_change(:)
    real(dp) :: auxiliary(size(source_change)),source_time
    logical, optional, intent(in) :: route_auxiliary
    logical :: routing
    type(amr_tree_state_3d) :: off_old,off_new
    type(amr_tree_state_3d), allocatable :: auxiliary_pending(:)
    real(dp) :: redist_target
    character(len=16) :: method
    type(amr_tree_state_3d), allocatable :: pending(:)
    type(eb_reactive_tree_3d) :: candidate,starts,ends
    type(amr_tree_state_3d), allocatable :: ir
    type(amr_tree_state_3d) :: source0,source1
    type(eb_integrated_flux_3d), allocatable :: total(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:),counts(:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: clock_start(:),clock_dt(:)
    real(dp) :: boundary(size(boundary_change)),wall(size(wall_change)),chemical(size(chemical_change))
    real(dp) :: external(size(source_change)),theta
    integer :: nv,node,n(3),levels,c
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    if(present(evaluate_reaction).and..not.chemistry) return
    routing=.false.
    if(present(route_auxiliary)) routing=route_auxiliary
    if(.not.valid_tree(tree,nv).or..not.ieee_is_finite(dt).or.dt<=0) return
    if(any([size(boundary_change),size(wall_change),size(chemical_change),size(source_change)]/=nv)) return
    if(present(reaction_source).and..not.chemistry) return
    if(routing.and..not.present(auxiliary_source)) return
    if(present(auxiliary_source).neqv.present(auxiliary_change)) return
    if(present(auxiliary_change)) then
      if(size(auxiliary_change)/=nv) return
    end if
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if(.not.valid) return
    levels=size(ratios)
    if(present(level_steps)) then
      if(size(level_steps)/=levels) return
    end if
    method='none'
    if(present(redistribution)) method=redistribution
    if(method/='none'.and.method/='flux_redist'.and.method/='state_redist') return
    redist_target=.5_dp
    if(present(target_volume_fraction)) redist_target=target_volume_fraction
    if(.not.ieee_is_finite(redist_target)) return
    if(redist_target<=0.or.redist_target>1) return
    candidate=tree
    call average_down_reactive_eb_tree_3d(species,candidate,valid,frame)
    if(.not.valid) return
    starts=candidate;ends=candidate
    allocate(clock_start(levels),clock_dt(levels),counts(levels),source0%fields(size(tree%nodes)))
    counts=0;clock_start=0;clock_dt=dt
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      allocate(source0%fields(node)%values(nv,n(1),n(2),n(3)),source=0._dp)
    end do
    source1=source0
    if(routing) then
      allocate(auxiliary_pending(0:levels-1))
      do c=0,levels-1
        auxiliary_pending(c)=source0
      end do
    end if
    if(method/='none') then
      allocate(pending(0:levels-1))
      do c=0,levels-1
        pending(c)=source0
      end do
    end if
    if(chemistry) then
      allocate(ir)
      if(present(reaction_source)) then
        ir=reaction_source
      else
        allocate(ir%fields(size(tree%nodes)))
        do node=1,size(tree%nodes)
          n=tree%topology%nodes(node)%shape
          allocate(ir%fields(node)%values(nv-4,n(1),n(2),n(3)),source=0._dp)
        end do
      end if
    end if
    auxiliary=0;source_time=0
    wall=0;chemical=0;external=0;theta=1
    call advance_level(0,0._dp,dt,total,valid)
    if(.not.valid) return
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
    if(.not.all(ieee_is_finite(boundary)).or..not.all(ieee_is_finite(wall)).or. &
      .not.all(ieee_is_finite(chemical)).or..not.all(ieee_is_finite(external))) return
    tree=candidate;boundary_change=boundary;wall_change=wall;chemical_change=chemical;source_change=external
    if(present(auxiliary_change)) auxiliary_change=auxiliary
    minimum_theta=theta
    if(present(reaction_source)) reaction_source=ir
    if(present(level_steps)) level_steps=counts
    ok=.true.
  contains
    subroutine source_adapter(state,interval,stage,depth,rate,ok)
      type(eb_reactive_tree_3d), intent(in) :: state
      real(dp), intent(in) :: interval
      integer, intent(in) :: stage,depth
      type(amr_tree_state_3d), intent(inout) :: rate
      logical, intent(out) :: ok
      type(amr_tree_state_3d) :: off
      integer :: p,x,y,z,d(3)
      logical :: valid,agreed
      call auxiliary_source(state,interval,stage,depth,source_time,rate,ok)
      if(.not.routing) return
      if(.not.ok) return
      ok=.false.
      if(.not.allocated(rate%fields)) return
      if(size(rate%fields)/=size(tree%nodes)) return
      off=source0;valid=.true.
      do p=1,size(tree%nodes)
        if(allocated(off%fields(p)%values)) off%fields(p)%values=0
      end do
      do p=1,size(tree%nodes)
        if(tree%topology%nodes(p)%level==depth) cycle
        d=tree%topology%nodes(p)%shape
        if(.not.allocated(rate%fields(p)%values)) return
        if(any(shape(rate%fields(p)%values)/=[nv,d])) return
        off%fields(p)%values=0
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(state%nodes(p)%geometry%volume_fraction(x,y,z)<=0.or.covered(tree%topology,p,[x,y,z])) cycle
          off%fields(p)%values(:,x,y,z)=rate%fields(p)%values(:,x,y,z)
        end do
        end do
        end do
        valid=valid.and.all(ieee_is_finite(off%fields(p)%values))
      end do
      if(.not.valid) return
      if(stage==0) then
        off_old=off
      else
        off_new=off
      end if
      ok=.true.

    end subroutine source_adapter
    recursive subroutine advance_level(level,time,interval,flux,good)
      integer, intent(in) :: level
      real(dp), intent(in) :: time,interval
      type(eb_integrated_flux_3d), allocatable, intent(out) :: flux(:)
      logical, intent(out) :: good
      type(eb_integrated_flux_3d), allocatable :: fine_flux(:)
      type(amr_tree_state_3d) :: op,np,correction,redistributed
      type(amr_tree_state_3d) :: off_transfer
      type(eb_reactive_tree_3d) :: off_ledger
      real(dp) :: off_inventory(nv)
      real(dp), allocatable :: step_aux(:)
      type(amr_tree_state_3d), allocatable :: transfer
      type(eb_reactive_tree_3d) :: before_sync
      real(dp) :: fstart(level),fend(level),b(nv),w(nv),ch(nv),ext(nv),limiter,elapsed,sub_dt
      integer :: p,d(3),step,target,bucket
      logical :: accepted
      good=.false.
      if(counts(level+1)==huge(1_int64)) return
      clock_start(level+1)=time;clock_dt(level+1)=interval
      fstart=max(0._dp,min(1._dp,(time-clock_start(:level))/clock_dt(:level)))
      fend=max(fstart,min(1._dp,(time+interval-clock_start(:level))/clock_dt(:level)))
      do p=1,size(tree%nodes)
        if(tree%topology%nodes(p)%level==level) starts%nodes(p)=candidate%nodes(p)
      end do
      if(method/='none') allocate(transfer)
      if(present(auxiliary_source)) allocate(step_aux(nv))
      source_time=time
      if(present(auxiliary_source)) then
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,interval,iterations,chemistry, &
        transport_flags,force,rtol,atol,candidate,extended_geometry,boundary_faces,b,w,ch,ext,limiter,accepted,ir, &
        difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature, &
        level,starts,ends,fstart,fend,source0,source1,flux,op,np,method,transfer,redist_target,source_adapter,step_aux, &
          evaluate_reaction=evaluate_reaction,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs, &
          embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip,frame=frame)
      else
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,interval,iterations,chemistry, &
        transport_flags,force,rtol,atol,candidate,extended_geometry,boundary_faces,b,w,ch,ext,limiter,accepted,ir, &
        difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature, &
        level,starts,ends,fstart,fend,source0,source1,flux,op,np,method,transfer,redist_target, &
          evaluate_reaction=evaluate_reaction,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs, &
          embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip,frame=frame)
      end if
      if(.not.accepted) return
      if(method/='none') then
        do p=1,size(tree%nodes)
          target=tree%topology%nodes(p)%level;bucket=min(level,target)
          pending(bucket)%fields(p)%values=pending(bucket)%fields(p)%values+transfer%fields(p)%values
        end do
      end if
      if(routing) then
        off_transfer=off_old
        do p=1,size(tree%nodes)
          target=tree%topology%nodes(p)%level;bucket=min(level,target)
          off_transfer%fields(p)%values=.5_dp*interval*(off_old%fields(p)%values+off_new%fields(p)%values)
          auxiliary_pending(bucket)%fields(p)%values=auxiliary_pending(bucket)%fields(p)%values+ &
            off_transfer%fields(p)%values
        end do
        off_ledger=candidate
        do p=1,size(tree%nodes)
          off_ledger%nodes(p)%state=off_transfer%fields(p)%values
        end do
        call integrate_reactive_eb_tree_3d(off_ledger,off_inventory,accepted)
        if(.not.accepted) return
        step_aux=step_aux+off_inventory
      end if
      if(present(auxiliary_source)) then
        auxiliary=auxiliary+step_aux
        if(.not.all(ieee_is_finite(auxiliary))) return
      end if
      counts(level+1)=counts(level+1)+1
      wall=wall+w;chemical=chemical+ch;external=external+ext;theta=min(theta,limiter)
      do p=1,size(tree%nodes)
        if(tree%topology%nodes(p)%level==level) then
          ends%nodes(p)=candidate%nodes(p)
          source0%fields(p)%values=op%fields(p)%values;source1%fields(p)%values=np%fields(p)%values
        else
          d=tree%topology%nodes(p)%shape
          allocate(flux(p)%x(nv,0:d(1),d(2),d(3)),flux(p)%y(nv,d(1),0:d(2),d(3)), &
            flux(p)%z(nv,d(1),d(2),0:d(3)))
          flux(p)%x=0;flux(p)%y=0;flux(p)%z=0
        end if
      end do
      if(level+1<levels) then
        elapsed=0
        do step=1,ratios(level+2)
          sub_dt=interval/real(ratios(level+2),dp)
          if(step==ratios(level+2)) sub_dt=interval-elapsed
          if(sub_dt<=0.or.time+elapsed+sub_dt<=time+elapsed) return
          call advance_level(level+1,time+elapsed,sub_dt,fine_flux,accepted)
          if(.not.accepted) return
          do p=1,size(tree%nodes)
            if(tree%topology%nodes(p)%level/=level+1) cycle
            flux(p)%x=flux(p)%x+fine_flux(p)%x
            flux(p)%y=flux(p)%y+fine_flux(p)%y
            flux(p)%z=flux(p)%z+fine_flux(p)%z
          end do
          elapsed=elapsed+sub_dt
        end do
      end if
      if(method/='none') before_sync=candidate
      call synchronize_reactive_eb_tree_flux_3d(species,candidate,flux,accepted,level, &
        min(level+1,levels-1),boundary_faces,frame)
      if(.not.accepted) return
      if(method/='none') then
        correction=source0
        do p=1,size(tree%nodes)
          correction%fields(p)%values=0
          target=tree%topology%nodes(p)%level
          if(target<level.or.target>level+1) cycle
          correction%fields(p)%values=candidate%nodes(p)%state-before_sync%nodes(p)%state
        end do
        ! This operator is linear in its argument, so an integrated reflux
        ! increment can be redistributed directly without division by dt.
        if(method=='state_redist') then
          call redistribute_eb_tree_state_3d(candidate,correction,redistributed,accepted,redist_target,boundary_faces)
        else
          call redistribute_sdc_eb_tree_flux_rate_3d(candidate,correction,redistributed,accepted,boundary_faces)
        end if
        if(.not.accepted) return
        candidate=before_sync
        do p=1,size(tree%nodes)
          target=tree%topology%nodes(p)%level
          if(target<level) then
            pending(target)%fields(p)%values=pending(target)%fields(p)%values+redistributed%fields(p)%values
          else
            candidate%nodes(p)%state=candidate%nodes(p)%state+redistributed%fields(p)%values+ &
              pending(level)%fields(p)%values
          end if
          pending(level)%fields(p)%values=0
        end do
        call recover_transfer(accepted)
        if(.not.accepted) return
        call average_down_reactive_eb_tree_3d(species,candidate,accepted,frame)
        if(.not.accepted) return
      end if
      if(routing) then
        do p=1,size(tree%nodes)
          candidate%nodes(p)%state=candidate%nodes(p)%state+auxiliary_pending(level)%fields(p)%values
          auxiliary_pending(level)%fields(p)%values=0
        end do
        call recover_transfer(accepted)
        if(.not.accepted) return
        call average_down_reactive_eb_tree_3d(species,candidate,accepted,frame)
        if(.not.accepted) return
      end if
      good=.true.
    end subroutine advance_level
    subroutine recover_transfer(good)
      logical, intent(out) :: good
      real(dp) :: primitive(reactive_nprim(size(species))),temperature,speed
      integer :: p,x,y,z,d(3)
      logical :: accepted
      good=.false.
      do p=1,size(tree%nodes)
        d=tree%topology%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(candidate%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
          if(any(candidate%nodes(p)%state(6:,x,y,z)<0._dp)) return
          call recover_reactive_eb_cell_3d(species,candidate%nodes(p)%geometry,[x,y,z],candidate%nodes(p)%state(:,x,y,z), &
            candidate%nodes(p)%temperature(x,y,z),primitive,temperature,speed,accepted,frame)
          if(.not.accepted) return
          candidate%nodes(p)%temperature(x,y,z)=temperature
        end do
        end do
        end do
      end do
      good=.true.
    end subroutine recover_transfer
  end subroutine advance_reactive_eb_tree_sdc_subcycle_3d

  subroutine advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,iterations,chemistry, &
      transport_flags,force,rtol,atol,tree,extended_geometry,boundary_faces, &
      boundary_change,wall_change,chemical_change,source_change,minimum_theta,ok,reaction_source, &
      difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature, &
      target_depth,ancestor_start,ancestor_end,fractions_start,fractions_end,source_start,source_end, &
      integrated_flux,old_predictor,new_predictor,redistribution,redistribution_transfer,target_volume_fraction, &
      auxiliary_source,auxiliary_change, &
      evaluate_reaction,mol_mode,mol_lax_friedrichs,embedded_wall_temperature,embedded_noslip,frame)
    real(dp), optional, intent(in) :: embedded_wall_temperature
    logical, optional, intent(in) :: mol_mode,mol_lax_friedrichs,embedded_noslip
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: dt,force(3),rtol,atol
    integer, intent(in) :: iterations
    logical, intent(in) :: chemistry,transport_flags(4)
    type(eb_reactive_tree_3d), intent(inout) :: tree
    type(eb_geometry_3d), intent(in) :: extended_geometry(:)
    character(len=*), intent(in) :: boundary_faces(6)
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemical_change(:),source_change(:),minimum_theta
    logical, intent(out) :: ok
    integer, optional, intent(in) :: target_depth
    type(eb_reactive_tree_3d), optional, intent(in) :: ancestor_start,ancestor_end
    real(dp), optional, intent(in) :: fractions_start(:),fractions_end(:)
    type(amr_tree_state_3d), optional, intent(in) :: source_start,source_end
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: integrated_flux(:)
    type(amr_tree_state_3d), optional, intent(inout) :: old_predictor,new_predictor
    type(amr_tree_state_3d), optional, intent(inout) :: redistribution_transfer
    type(amr_tree_state_3d) :: transfer0,transfer1,hydro_transfer,final_transfer
    real(dp), optional, intent(in) :: target_volume_fraction
    procedure(eb_sdc_auxiliary_source), optional :: auxiliary_source
    real(dp), optional, intent(inout) :: auxiliary_change(:)
    type(amr_tree_state_3d) :: a0,a1
    real(dp) :: auxiliary(size(source_change))
    real(dp) :: redist_target
    character(len=*), optional, intent(in) :: redistribution
    character(len=16) :: redistribution_method
    type(amr_tree_state_3d) :: redistributed,sync_rate
    type(eb_reactive_tree_3d) :: before_sync
    type(amr_tree_state_3d) :: op,np
    real(dp), allocatable :: midpoint(:)
    integer :: selected_depth
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source
    real(dp), optional, intent(in) :: difmag,inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: wall_temperature(6)
    type(eb_reactive_tree_3d) :: old,candidate,ledger
    type(amr_tree_state_3d) :: e0,e1,d0,d1,rates,hydro,ir,u,t
    type(eb_integrated_flux_3d), allocatable :: f0(:),f1(:),hf(:)
    real(dp) :: boundary(size(boundary_change)),wall(size(wall_change)),chemical(size(chemical_change))
    real(dp) :: external(size(source_change)),v0(size(source_change)),v1(size(source_change)),theta,stage_theta
    real(dp), allocatable :: primitive(:),wall_rate(:,:,:,:)
    real(dp) :: recovered,speed
    integer :: nv,node,n(3),it,i,j,k,c
    logical :: valid,diffusive,mol,noslip,wall_active
    type(eb_geometry_3d), allocatable :: mol_geometry(:),heat_geometry(:)
    real(dp), allocatable :: heat_stage,heat_frequency,wall_impulse(:)
    real(dp) :: heat_old,heat_new,impulse_old(3),impulse_new(3)
    real(dp) :: mol_wall_old(size(wall_change)),mol_wall_new(size(wall_change))
    ok=.false.;nv=reactive_nvar(size(species));theta=1
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    mol=.false.
    if(present(mol_mode)) mol=mol_mode
    if(present(mol_lax_friedrichs)) then
      if(mol_lax_friedrichs.and..not.mol) return
    end if
    if(present(frame)) then
      if(.not.mol) return
    end if
    mol_wall_old=0;mol_wall_new=0
    heat_old=0;heat_new=0;impulse_old=0;impulse_new=0
    noslip=.false.
    if(present(embedded_noslip)) noslip=embedded_noslip
    wall_active=noslip.or.present(embedded_wall_temperature)
    if(noslip.and..not.transport_flags(1)) return
    if(present(embedded_wall_temperature)) then
      if(.not.ieee_is_finite(embedded_wall_temperature).or.embedded_wall_temperature<=0) return
      if(.not.transport_flags(2)) return
    end if
    if(present(evaluate_reaction).and..not.chemistry) return
    if(.not.valid_tree(tree,nv)) return
    if(any([size(boundary_change),size(wall_change),size(chemical_change),size(source_change)]/=nv)) return
    if(iterations<1.or..not.all(ieee_is_finite([dt,force,rtol,atol]))) return
    if(min(dt,rtol,atol)<=0.or..not.valid_reactive_boundaries_3d(boundary_faces)) return
    if(transport_flags(4).and..not.transport_flags(3)) return
    if(present(reaction_source).and..not.chemistry) return
    if(present(auxiliary_source).neqv.present(auxiliary_change)) return
    if(present(auxiliary_change)) then
      if(size(auxiliary_change)/=nv) return
    end if
    selected_depth=-1
    if(present(target_depth)) then
      selected_depth=target_depth
      if(selected_depth<0.or.selected_depth>=tree%topology%level_count()) return
      if(.not.present(integrated_flux)) return
    end if
    redistribution_method='none'
    if(present(redistribution)) redistribution_method=redistribution
    if(redistribution_method/='none'.and.redistribution_method/='flux_redist'.and. &
      redistribution_method/='state_redist') return
    redist_target=.5_dp
    if(present(target_volume_fraction)) redist_target=target_volume_fraction
    if(.not.ieee_is_finite(redist_target)) return
    if(redist_target<=0.or.redist_target>1) return
    if(selected_depth>=0.and.redistribution_method/='none') then
      if(.not.present(redistribution_transfer)) return
    else
      if(present(redistribution_transfer)) return
    end if
    if(present(ancestor_start).or.present(ancestor_end).or.present(fractions_start).or.present(fractions_end)) then
      if(selected_depth<0) return
      if(.not.(present(ancestor_start).and.present(ancestor_end).and. &
        present(fractions_start).and.present(fractions_end))) return
      if(size(fractions_start)/=selected_depth.or.size(fractions_end)/=selected_depth) return
      if(.not.all(ieee_is_finite(fractions_start)).or..not.all(ieee_is_finite(fractions_end))) return
      if(any(fractions_start<0).or.any(fractions_end>1).or.any(fractions_end<fractions_start)) return
    end if
    if(present(source_start).neqv.present(source_end)) return
    if(present(source_start)) then
      if(.not.present(fractions_start).or..not.present(fractions_end)) return
      midpoint=.5_dp*(fractions_start+fractions_end)
    end if
    old=tree
    if(selected_depth<0) then
      call average_down_reactive_eb_tree_3d(species,old,valid,frame)
      if(.not.valid) return
    end if
    if(mol.or.wall_active) then
      if(size(extended_geometry)/=size(tree%nodes)) return
      if(mol) allocate(mol_geometry(size(tree%nodes)))
      if(wall_active) then
        allocate(heat_geometry(size(tree%nodes)),heat_stage,heat_frequency)
        if(noslip) allocate(wall_impulse(3))
      end if
      do node=1,size(tree%nodes)
        n=tree%topology%nodes(node)%shape
        if(any([extended_geometry(node)%nx,extended_geometry(node)%ny,extended_geometry(node)%nz]/=n+12)) return
        if(mol) then
          call extract_eb_geometry_3d(extended_geometry(node),[4,4,4],n+9,mol_geometry(node),valid)
          if(.not.valid) return
        end if
        if(wall_active) then
          call extract_eb_geometry_3d(extended_geometry(node),[5,5,5],n+8,heat_geometry(node),valid)
          if(.not.valid) return
        end if
      end do
    end if
    allocate(primitive(reactive_nprim(size(species))),ir%fields(size(tree%nodes)))
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      allocate(ir%fields(node)%values(nv-4,n(1),n(2),n(3)),source=0._dp)
    end do
    if(present(reaction_source)) then
      if(.not.allocated(reaction_source%fields)) return
      if(size(reaction_source%fields)/=size(tree%nodes)) return
      do node=1,size(tree%nodes)
        if(.not.allocated(reaction_source%fields(node)%values)) return
        if(any(shape(reaction_source%fields(node)%values)/=shape(ir%fields(node)%values))) return
      end do
      ir=reaction_source
    end if
    do node=1,size(tree%nodes)
      if(.not.active(node)) cycle
      do c=1,nv-4
        where(old%nodes(node)%geometry%volume_fraction<=0) ir%fields(node)%values(c,:,:,:)=0
      end do
      if(.not.all(ieee_is_finite(ir%fields(node)%values))) return
    end do
    candidate=old;diffusive=any(transport_flags(:3))
    call sources(old,e0,d0,a0,0,f0,transfer0,valid,fractions_start)
    if(.not.valid) return
    e1=e0;d1=d0;mol_wall_new=mol_wall_old;heat_new=heat_old;impulse_new=impulse_old
    if(present(auxiliary_source)) a1=a0
    call predictor_field(e0,d0,a0,op)
    allocate(u%fields(size(tree%nodes)),t%fields(size(tree%nodes)))
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      u%fields(node)%values=old%nodes(node)%state
      allocate(t%fields(node)%values(1,n(1),n(2),n(3)))
      t%fields(node)%values(1,:,:,:)=old%nodes(node)%temperature
    end do
    if(mol) then
      hydro=d0;hf=f0
      do node=1,size(tree%nodes)
        hydro%fields(node)%values=0
        if(.not.active(node)) cycle
        hf(node)%x=0;hf(node)%y=0;hf(node)%z=0
      end do
      if(present(redistribution_transfer)) hydro_transfer=hydro
    end if
    do it=1,merge(1,iterations,mol.and..not.chemistry)
      ! MOL builds an Euler predictor only for the first iteration. Later
      ! reacting iterations use the preceding reacted endpoint directly.
      if(.not.mol.or.it==1) then
      call combined_rates(.true.)
      if(.not.mol) then
      call compute_sdc_eb_tree_hydro_3d(species,old,rates,extended_geometry,dt,hydro,hf,valid, &
        difmag,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        target_depth,ancestor_start,ancestor_end,fractions_start,source_start,source_end,midpoint)
      if(.not.valid) return
      if(redistribution_method/='none') then
        if(selected_depth>=0) then
          if(redistribution_method=='state_redist') then
            call redistribute_sdc_eb_tree_state_level_rate_3d(old,hydro,selected_depth,dt,redistributed, &
              hydro_transfer,valid,redist_target,boundary_faces)
          else
            call redistribute_sdc_eb_tree_level_rate_3d(old,hydro,selected_depth,dt,redistributed, &
              hydro_transfer,valid,boundary_faces)
          end if
        else
          if(redistribution_method=='state_redist') then
            call redistribute_sdc_eb_tree_state_rate_3d(old,hydro,dt,redistributed,valid,redist_target,boundary_faces)
          else
            call redistribute_sdc_eb_tree_flux_rate_3d(old,hydro,redistributed,valid,boundary_faces)
          end if
        end if
        if(.not.valid) return
        hydro=redistributed
      end if
      end if
      do node=1,size(tree%nodes)
        if(.not.active(node)) cycle
        candidate%nodes(node)%state=old%nodes(node)%state+dt*(rates%fields(node)%values+hydro%fields(node)%values)
      end do
      call recover(valid)
      if(.not.valid) return
      end if
      call sources(candidate,e1,d1,a1,it,f1,transfer1,valid,fractions_end)
      if(.not.valid) return
      call combined_rates(.false.)
      do node=1,size(tree%nodes)
        rates%fields(node)%values=rates%fields(node)%values+hydro%fields(node)%values
      end do
      if(chemistry) then
        call react_sdc_eb_tree_3d(species,reactions,old,rates,dt,rtol,atol,u,t,ir,valid,target_depth, &
          evaluate_reaction=evaluate_reaction,frame=frame)
        if(.not.valid) return
        do node=1,size(tree%nodes)
          if(.not.active(node)) cycle
          candidate%nodes(node)%state=u%fields(node)%values
          candidate%nodes(node)%temperature=t%fields(node)%values(1,:,:,:)
        end do
      else
        do node=1,size(tree%nodes)
          if(.not.active(node)) cycle
          candidate%nodes(node)%state=old%nodes(node)%state+dt*rates%fields(node)%values
        end do
        call recover(valid)
        if(.not.valid) return
      end if
    end do
    ledger=old
    do node=1,size(tree%nodes)
      ledger%nodes(node)%state=0
      if(.not.active(node)) cycle
      ledger%nodes(node)%state=candidate%nodes(node)%state-(old%nodes(node)%state+dt*rates%fields(node)%values)
    end do
    call integrate_reactive_eb_tree_3d(ledger,chemical,valid)
    if(.not.valid) return
    call field_inventory(e0,v0,valid)
    if(.not.valid) return
    call field_inventory(e1,v1,valid)
    if(.not.valid) return
    external=.5_dp*dt*(v0+v1)
    if(present(auxiliary_source)) then
      call field_inventory(a0,v0,valid)
      if(.not.valid) return
      call field_inventory(a1,v1,valid)
      if(.not.valid) return
      auxiliary=.5_dp*dt*(v0+v1)
      if(.not.all(ieee_is_finite(auxiliary))) return
    end if
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      ledger%nodes(node)%state=0
      if(.not.active(node)) cycle
      allocate(wall_rate(nv,n(1),n(2),n(3)))
      call wall_source(old,node,wall_rate,valid)
      if(.not.valid) return
      ledger%nodes(node)%state=dt*wall_rate
      deallocate(wall_rate)
      if(diffusive.or.mol) then
        hf(node)%x=hf(node)%x+.5_dp*(f0(node)%x+f1(node)%x)
        hf(node)%y=hf(node)%y+.5_dp*(f0(node)%y+f1(node)%y)
        hf(node)%z=hf(node)%z+.5_dp*(f0(node)%z+f1(node)%z)
      end if
    end do
    call integrate_reactive_eb_tree_3d(ledger,wall,valid)
    if(.not.valid) return
    if(mol) wall=.5_dp*(mol_wall_old+mol_wall_new)
    wall(5)=wall(5)+.5_dp*(heat_old+heat_new)
    wall(2:4)=wall(2:4)+.5_dp*(impulse_old+impulse_new)
    if(.not.all(ieee_is_finite(wall))) return
    if(selected_depth<0) then
      if(redistribution_method/='none') before_sync=candidate
      call synchronize_reactive_eb_tree_flux_3d(species,candidate,hf,valid,boundary_faces=boundary_faces,frame=frame)
      if(.not.valid) return
      if(redistribution_method/='none') then
        ! Redistribute the reflux increment too: applying the same linear
        ! operator to rates and flux corrections preserves composite closure.
        sync_rate=hydro
        do node=1,size(tree%nodes)
          sync_rate%fields(node)%values=(candidate%nodes(node)%state-before_sync%nodes(node)%state)/dt
        end do
        if(redistribution_method=='state_redist') then
          call redistribute_eb_tree_state_3d(old,sync_rate,redistributed,valid,redist_target,boundary_faces)
        else
          call redistribute_sdc_eb_tree_flux_rate_3d(old,sync_rate,redistributed,valid,boundary_faces)
        end if
        if(.not.valid) return
        do node=1,size(tree%nodes)
          candidate%nodes(node)%state=before_sync%nodes(node)%state+dt*redistributed%fields(node)%values
        end do
        call recover(valid)
        if(.not.valid) return
        call average_down_reactive_eb_tree_3d(species,candidate,valid,frame)
        if(.not.valid) return
      end if
    end if
    boundary=0
    if(selected_depth<=0) then
    n=tree%topology%nodes(1)%shape
    associate(g=>tree%nodes(1)%geometry,f=>hf(1))
    do c=1,nv
      boundary(c)=g%dy*g%dz*sum(g%x_face_fraction(0,:,:)*f%x(c,0,:,:)- &
        g%x_face_fraction(n(1),:,:)*f%x(c,n(1),:,:))+ &
        g%dx*g%dz*sum(g%y_face_fraction(:,0,:)*f%y(c,:,0,:)- &
        g%y_face_fraction(:,n(2),:)*f%y(c,:,n(2),:))+ &
        g%dx*g%dy*sum(g%z_face_fraction(:,:,0)*f%z(c,:,:,0)- &
        g%z_face_fraction(:,:,n(3))*f%z(c,:,:,n(3)))
    end do
    end associate
    end if
    if(.not.all(ieee_is_finite(boundary)).or..not.all(ieee_is_finite(external))) return
    call predictor_field(e1,d1,a1,np)
    if(.not.op%is_valid(tree%topology).or..not.np%is_valid(tree%topology)) return
    if(present(redistribution_transfer)) then
      final_transfer=hydro_transfer
      do node=1,size(tree%nodes)
        final_transfer%fields(node)%values=hydro_transfer%fields(node)%values+ &
          .5_dp*(transfer0%fields(node)%values+transfer1%fields(node)%values)
      end do
      if(.not.final_transfer%is_valid(tree%topology)) return
    end if
    ! EOS recovery permits roundoff-sized negative trial densities. Published
    ! states must satisfy the strict nonnegative checkpoint contract instead.
    ! Check after reflux/redistribution and before publishing any output.
    do node=1,size(tree%nodes)
      if(.not.active(node)) cycle
      n=tree%topology%nodes(node)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(candidate%nodes(node)%geometry%volume_fraction(i,j,k)<=0) cycle
        if(any(candidate%nodes(node)%state(6:,i,j,k)<0._dp)) return
      end do
      end do
      end do
    end do
    if(present(reaction_source)) then
      do node=1,size(tree%nodes)
        do c=1,nv-4
          where(old%nodes(node)%geometry%volume_fraction<=0) &
            ir%fields(node)%values(c,:,:,:)=reaction_source%fields(node)%values(c,:,:,:)
        end do
      end do
      reaction_source=ir
    end if
    if(present(integrated_flux)) integrated_flux=hf
    if(present(old_predictor)) old_predictor=op
    if(present(new_predictor)) new_predictor=np
    if(present(auxiliary_change)) auxiliary_change=auxiliary
    if(present(redistribution_transfer)) redistribution_transfer=final_transfer
    tree=candidate;boundary_change=boundary;wall_change=wall;chemical_change=chemical
    source_change=external;minimum_theta=theta;ok=.true.
  contains
    logical function active(p)
      integer, intent(in) :: p
      active=selected_depth<0.or.tree%topology%nodes(p)%level==selected_depth
    end function active
    subroutine predictor_field(ext,diff,aux,output)
      type(amr_tree_state_3d), intent(in) :: ext,diff,aux
      type(amr_tree_state_3d), intent(out) :: output
      integer :: p
      output=ext
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        output%fields(p)%values=ext%fields(p)%values+diff%fields(p)%values
        if(present(auxiliary_source)) output%fields(p)%values=output%fields(p)%values+aux%fields(p)%values
        output%fields(p)%values(6:,:,:,:)=output%fields(p)%values(6:,:,:,:)+ir%fields(p)%values(:nv-5,:,:,:)
        output%fields(p)%values(5,:,:,:)=output%fields(p)%values(5,:,:,:)+ir%fields(p)%values(nv-4,:,:,:)
      end do
    end subroutine predictor_field
    subroutine combined_rates(with_reaction)
      logical, intent(in) :: with_reaction
      integer :: p
      rates=e0
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        rates%fields(p)%values=.5_dp*(e0%fields(p)%values+e1%fields(p)%values+ &
          d0%fields(p)%values+d1%fields(p)%values)
        if(present(auxiliary_source)) rates%fields(p)%values=rates%fields(p)%values+ &
          .5_dp*(a0%fields(p)%values+a1%fields(p)%values)
        if(.not.with_reaction) cycle
        rates%fields(p)%values(6:,:,:,:)=rates%fields(p)%values(6:,:,:,:)+ir%fields(p)%values(:nv-5,:,:,:)
        rates%fields(p)%values(5,:,:,:)=rates%fields(p)%values(5,:,:,:)+ir%fields(p)%values(nv-4,:,:,:)
      end do
    end subroutine combined_rates
    subroutine sources(state,external_rate,transport_rate,auxiliary_rate,stage,flux,exchange,good,fractions)
      type(eb_reactive_tree_3d), intent(in) :: state
      type(amr_tree_state_3d), intent(out) :: external_rate,transport_rate,auxiliary_rate,exchange
      integer, intent(in) :: stage
      type(eb_integrated_flux_3d), allocatable, intent(out) :: flux(:)
      logical, intent(out) :: good
      real(dp), optional, intent(in) :: fractions(:)
      integer :: p,d(3),x,y,z
      real(dp) :: position(3),rotation_source(3)
      type(amr_tree_state_3d) :: extra
      type(eb_integrated_flux_3d), allocatable :: extra_flux(:)
      good=.false.
      allocate(external_rate%fields(size(tree%nodes)))
      do p=1,size(tree%nodes)
        d=tree%topology%nodes(p)%shape
        allocate(external_rate%fields(p)%values(nv,d(1),d(2),d(3)),source=0._dp)
        if(.not.active(p)) cycle
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(state%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
          if(state%nodes(p)%state(1,x,y,z)<=0) return
          external_rate%fields(p)%values(2:4,x,y,z)=force
          external_rate%fields(p)%values(5,x,y,z)=dot_product(force,state%nodes(p)%state(2:4,x,y,z))/ &
            state%nodes(p)%state(1,x,y,z)
          if(present(frame)) then
            associate(g=>state%nodes(p)%geometry)
            position=[g%x_lower,g%y_lower,g%z_lower]+(real([x,y,z],dp)-.5_dp)*[g%dx,g%dy,g%dz]
            end associate
            call rotating_momentum_source_3d(frame,position,state%nodes(p)%state(1,x,y,z), &
              state%nodes(p)%state(2:4,x,y,z),rotation_source,good)
            if(.not.good) return
            good=.false.
            external_rate%fields(p)%values(2:4,x,y,z)=external_rate%fields(p)%values(2:4,x,y,z)+rotation_source
            ! Rotating total energy has no rotational source term.
          end if
        end do
        end do
        end do
      end do
      transport_rate=external_rate
      do p=1,size(tree%nodes)
        transport_rate%fields(p)%values=0
      end do
      if(present(auxiliary_source)) then
        auxiliary_rate=transport_rate
        call auxiliary_source(state,dt,stage,selected_depth,auxiliary_rate,good)
        if(.not.good) return
        good=.false.
        if(.not.allocated(auxiliary_rate%fields)) return
        if(size(auxiliary_rate%fields)/=size(tree%nodes)) return
        do p=1,size(tree%nodes)
          d=tree%topology%nodes(p)%shape
          if(.not.allocated(auxiliary_rate%fields(p)%values)) return
          if(any(shape(auxiliary_rate%fields(p)%values)/=[nv,d])) return
          if(.not.active(p)) then
            auxiliary_rate%fields(p)%values=0;cycle
          end if
          do z=1,d(3)
          do y=1,d(2)
          do x=1,d(1)
            if(state%nodes(p)%geometry%volume_fraction(x,y,z)<=0) auxiliary_rate%fields(p)%values(:,x,y,z)=0
          end do
          end do
          end do
          if(.not.all(ieee_is_finite(auxiliary_rate%fields(p)%values))) return
        end do
      end if
      if(selected_depth>=0.and.redistribution_method/='none') exchange=transport_rate
      if(diffusive) then
        call compute_sdc_eb_tree_transport_3d(species,transport,dt,transport_flags(1),transport_flags(2), &
          transport_flags(3),transport_flags(4),state,transport_rate,flux,stage_theta,good, &
          inflow_mach,inflow_temperature,inflow_pressure,boundary_faces,inflow_state,wall_temperature, &
          target_depth,ancestor_start,ancestor_end,fractions, &
          embedded_wall_temperature=embedded_wall_temperature,extended_geometry=heat_geometry, &
          embedded_wall_energy=heat_stage,embedded_frequency=heat_frequency, &
          embedded_noslip=embedded_noslip,embedded_wall_impulse=wall_impulse,frame=frame)
        if(.not.good) return
        if(wall_active) then
          if(stage==0) then
            heat_old=heat_stage
            if(noslip) impulse_old=wall_impulse
          else
            heat_new=heat_stage
            if(noslip) impulse_new=wall_impulse
          end if
        end if
      else
        stage_theta=1
      end if
      if(mol) then
        call compute_mol_eb_tree_hydro_3d(species,state,mol_geometry,dt,2,extra,extra_flux,good, &
          boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
          target_depth,ancestor_start,ancestor_end,fractions,mol_lax_friedrichs,frame)
        if(.not.good) return
        do p=1,size(tree%nodes)
          if(.not.active(p)) cycle
          transport_rate%fields(p)%values=transport_rate%fields(p)%values+extra%fields(p)%values
          if(diffusive) then
            flux(p)%x=flux(p)%x+extra_flux(p)%x
            flux(p)%y=flux(p)%y+extra_flux(p)%y
            flux(p)%z=flux(p)%z+extra_flux(p)%z
          end if
        end do
        if(.not.diffusive) call move_alloc(extra_flux,flux)
        do p=1,size(tree%nodes)
          extra%fields(p)%values=0
          if(.not.active(p)) cycle
          call wall_source(state,p,extra%fields(p)%values,good)
          if(.not.good) return
        end do
        if(stage==0) then
          call field_inventory(extra,mol_wall_old,good)
          mol_wall_old=dt*mol_wall_old
        else
          call field_inventory(extra,mol_wall_new,good)
          mol_wall_new=dt*mol_wall_new
        end if
        if(.not.good) return
      end if
      if(diffusive.or.mol) then
        if(redistribution_method/='none') then
          if(selected_depth>=0) then
            if(redistribution_method=='state_redist') then
              call redistribute_sdc_eb_tree_state_level_rate_3d(state,transport_rate,selected_depth,dt,redistributed, &
                exchange,good,redist_target,boundary_faces)
            else
              call redistribute_sdc_eb_tree_level_rate_3d(state,transport_rate,selected_depth,dt,redistributed, &
                exchange,good,boundary_faces)
            end if
          else
            if(redistribution_method=='state_redist') then
              call redistribute_sdc_eb_tree_state_rate_3d(state,transport_rate,dt,redistributed,good, &
                redist_target,boundary_faces)
            else
              call redistribute_sdc_eb_tree_flux_rate_3d(state,transport_rate,redistributed,good,boundary_faces)
            end if
          end if
          if(.not.good) return
          transport_rate=redistributed
        end if
        theta=min(theta,stage_theta)
      end if
      good=.true.
    end subroutine sources
    subroutine wall_source(state,p,output,good)
      type(eb_reactive_tree_3d), intent(in) :: state
      integer, intent(in) :: p
      real(dp), intent(out) :: output(:,:,:,:)
      logical, intent(out) :: good
      real(dp), allocatable :: potential(:,:,:)
      real(dp) :: position(3)
      integer :: x,y,z
      good=.false.;output=0
      associate(g=>state%nodes(p)%geometry)
      if(present(frame)) then
        allocate(potential(g%nx,g%ny,g%nz),source=0._dp)
        do z=1,g%nz
        do y=1,g%ny
        do x=1,g%nx
          if(g%volume_fraction(x,y,z)<=0) cycle
          position=[g%x_lower,g%y_lower,g%z_lower]+(real([x,y,z],dp)-.5_dp)*[g%dx,g%dy,g%dz]
          call rotating_potential_3d(frame,position,potential(x,y,z),good)
          if(.not.good) return
        end do
        end do
        end do
      end if
      call reactive_eb_slip_wall_source_3d(species,state%nodes(p)%state,state%nodes(p)%temperature, &
        g,output,good,potential)
      end associate
    end subroutine wall_source
    subroutine recover(good)
      logical, intent(out) :: good
      good=.false.
      do node=1,size(tree%nodes)
        if(.not.active(node)) cycle
        n=tree%topology%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(candidate%nodes(node)%geometry%volume_fraction(i,j,k)<=0) cycle
          call recover_reactive_eb_cell_3d(species,candidate%nodes(node)%geometry,[i,j,k],candidate%nodes(node)%state(:,i,j,k), &
            candidate%nodes(node)%temperature(i,j,k),primitive,recovered,speed,valid,frame)
          if(.not.valid) return
          candidate%nodes(node)%temperature(i,j,k)=recovered
        end do
        end do
        end do
      end do
      good=.true.
    end subroutine recover
    subroutine field_inventory(field,value,good)
      type(amr_tree_state_3d), intent(in) :: field
      real(dp), intent(inout) :: value(:)
      logical, intent(out) :: good
      type(eb_reactive_tree_3d) :: integral_tree
      integer :: p
      integral_tree=old
      do p=1,size(tree%nodes)
        integral_tree%nodes(p)%state=field%fields(p)%values
      end do
      call integrate_reactive_eb_tree_3d(integral_tree,value,good)
    end subroutine field_inventory
  end subroutine advance_reactive_eb_tree_sdc_3d

  subroutine advance_reactive_eb_tree_full_3d(species,reactions,transport,solver,dt,chemistry_enabled, &
      rtol,atol,transport_enabled,viscosity,conduction,diffusion,barodiffusion,tree,boundary_change, &
      wall_change,chemical_change,minimum_theta,ok,integrator,subcycling, &
      inflow_mach,inflow_temperature,inflow_pressure,difmag,extended_geometry,boundary_faces,inflow_state,wall_temperature)
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
    real(dp), optional, intent(in) :: wall_temperature(6)
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
    if(present(wall_temperature)) then
      if(.not.present(boundary_faces)) return
      if(.not.all(ieee_is_finite(wall_temperature)).or.any(wall_temperature<=0)) return
    end if
    if(present(boundary_faces)) then
      if(any(boundary_faces=='isothermal_wall')) then
        if(.not.present(wall_temperature).or..not.transport_enabled.or..not.conduction) return
      end if
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
            boundary_faces=boundary_faces,inflow_state=inflow_state,wall_temperature=wall_temperature)
        else
          call advance_reactive_eb_tree_transport_3d(species,transport,0.5_dp*dt,viscosity,conduction, &
            diffusion,barodiffusion,candidate,exchange,stage_theta,valid, &
            inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
            boundary_faces=boundary_faces,inflow_state=inflow_state,wall_temperature=wall_temperature)
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

  ! Two-cell donor margins for the quadratic embedded-wall stencil. Publish
  ! composite wall exchange only after every active patch succeeds.
  ! Transport has zero total-mass flux, so its energy flux needs no potential
  ! correction. Restore ordinary energy only in private input copies, after
  ! conserved rotating halos have been filled at their receiving positions.
  subroutine ordinary_transport_energy(frame,geometry,first,active,state,ok)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: first(3)
    logical, intent(in) :: active(:,:,:)
    real(dp), intent(inout) :: state(:,:,:,:)
    logical, intent(out) :: ok
    real(dp) :: position(3),energy
    integer :: i,j,k
    ok=.true.
    if(.not.present(frame)) return
    ok=.false.
    if(.not.frame%is_valid().or.size(state,1)<5) return
    if(any(shape(active)/=shape(state(1,:,:,:)))) return
    do k=1,size(active,3)
    do j=1,size(active,2)
    do i=1,size(active,1)
      if(.not.active(i,j,k)) cycle
      position=[geometry%x_lower,geometry%y_lower,geometry%z_lower]+ &
        (real(first+[i,j,k]-1,dp)-.5_dp)*[geometry%dx,geometry%dy,geometry%dz]
      energy=state(5,i,j,k)
      call convert_rotating_energy_3d(frame,position,state(1,i,j,k),state(5,i,j,k),.false.,energy,ok)
      if(.not.ok) return
      state(5,i,j,k)=energy
    end do
    end do
    end do
    ok=.true.
  end subroutine ordinary_transport_energy

  subroutine compute_eb_tree_wall_transport_3d(species,transport,tree,extended_geometry, &
      conduction,noslip,wall_temperature,rate,wall_exchange,frequency,ok, &
      boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,target_depth, &
      ancestor_start,ancestor_end,ancestor_fraction,frame)
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(gas_transport_species), intent(in) :: transport(:)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: extended_geometry(:)
    logical, intent(in) :: conduction,noslip
    real(dp), intent(in) :: wall_temperature
    type(amr_tree_state_3d), intent(inout) :: rate
    real(dp), intent(inout) :: wall_exchange(:),frequency
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    integer, optional, intent(in) :: target_depth
    type(eb_reactive_tree_3d), optional, intent(in) :: ancestor_start,ancestor_end
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    type(amr_tree_state_3d) :: candidate
    type(eb_reactive_tree_3d) :: ledger
    type(eb_geometry_3d) :: center
    real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),vf(:,:,:)
    real(dp) :: bound,node_frequency,inventory(reactive_nvar(size(species))),node_exchange(reactive_nvar(size(species)))
    integer :: nv,node,n(3),depth
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species));depth=-1
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    if(.not.valid_tree(tree,nv).or.size(transport)/=size(species)) return
    if(.not.ieee_is_finite(wall_temperature).or.size(wall_exchange)/=nv) return
    if(conduction.and.wall_temperature<=0) return
    if(size(extended_geometry)/=size(tree%nodes)) return
    if(present(target_depth)) then
      depth=target_depth
      if(depth<0.or.depth>=tree%topology%level_count()) return
    end if
    if(present(ancestor_start).or.present(ancestor_end).or.present(ancestor_fraction)) then
      if(depth<0) return
    end if
    allocate(candidate%fields(size(tree%nodes)));bound=0
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      if(any([extended_geometry(node)%nx,extended_geometry(node)%ny,extended_geometry(node)%nz]/=n+4)) return
      call extract_eb_geometry_3d(extended_geometry(node),[3,3,3],n+2,center,valid)
      if(.not.valid) return
      if(.not.matching_eb_metrics_3d(center,tree%nodes(node)%geometry)) return
      allocate(candidate%fields(node)%values(nv,n(1),n(2),n(3)),source=0._dp)
      if(depth>=0.and.tree%topology%nodes(node)%level/=depth) cycle
      allocate(ghost(nv,-1:n(1)+2,-1:n(2)+2,-1:n(3)+2), &
        gt(-1:n(1)+2,-1:n(2)+2,-1:n(3)+2),vf(-1:n(1)+2,-1:n(2)+2,-1:n(3)+2))
      call fill_reactive_eb_tree_halo_3d(species,tree,node,ghost,gt,valid, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_mach=inflow_mach, &
        inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,ghost_volume_fraction=vf, &
        halo_width=2,target_volume_fraction=extended_geometry(node)%volume_fraction,restrict_snapshot=.false., &
        ancestor_start=ancestor_start,ancestor_end=ancestor_end,ancestor_fraction=ancestor_fraction,frame=frame)
      if(.not.valid) return
      call ordinary_transport_energy(frame,tree%nodes(node)%geometry,[-1,-1,-1],vf>0,ghost,valid)
      if(.not.valid) return
      call compute_eb_wall_transport_source_3d(species,transport,ghost,gt,extended_geometry(node), &
        [3,3,3],n+2,conduction,noslip,wall_temperature,candidate%fields(node)%values,node_exchange,node_frequency,valid)
      if(.not.valid) return
      bound=max(bound,node_frequency)
      deallocate(ghost,gt,vf)
    end do
    ledger=tree
    do node=1,size(tree%nodes)
      ledger%nodes(node)%state=candidate%fields(node)%values
    end do
    ! Exclude covered coarse cells rather than summing per-patch ledgers.
    call integrate_reactive_eb_tree_3d(ledger,inventory,valid)
    if(.not.valid) return
    call move_alloc(candidate%fields,rate%fields)
    wall_exchange=inventory;frequency=bound;ok=.true.
  end subroutine compute_eb_tree_wall_transport_3d

  ! Common-interval molecular transport: closed outflow edges or open Mach inlets.
  ! Provisional transport rates and time-integrated face densities on every
  ! level. Parent states are retained until the coupled driver synchronizes.
  subroutine compute_sdc_eb_tree_transport_3d(species,transport,dt,viscosity,conduction, &
      diffusion,barodiffusion,tree,transport_rate,integrated_flux,minimum_theta,ok, &
      inflow_mach,inflow_temperature,inflow_pressure,boundary_faces,inflow_state,wall_temperature, &
      target_depth,ancestor_start,ancestor_end,ancestor_fraction, &
      embedded_wall_temperature,extended_geometry,embedded_wall_energy,embedded_frequency, &
      embedded_noslip,embedded_wall_impulse,frame)
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(gas_transport_species), intent(in) :: transport(:)
    logical, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), intent(in) :: dt
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(inout) :: transport_rate
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: integrated_flux(:)
    real(dp), intent(inout) :: minimum_theta
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: wall_temperature(6)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: target_depth
    type(eb_reactive_tree_3d), optional, intent(in) :: ancestor_start,ancestor_end
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    real(dp), optional, intent(in) :: embedded_wall_temperature
    type(eb_geometry_3d), optional, intent(in) :: extended_geometry(:)
    real(dp), optional, intent(inout) :: embedded_wall_energy,embedded_frequency,embedded_wall_impulse(3)
    logical, optional, intent(in) :: embedded_noslip
    logical :: noslip,wall_active
    real(dp) :: wall_values(reactive_nvar(size(species))),heat_frequency,heat_energy,wall_impulse(3),wall_temp
    type(amr_tree_state_3d) :: heat
    integer :: selected_depth
    type(amr_tree_state_3d) :: rates
    type(eb_integrated_flux_3d), allocatable :: flux(:)
    real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),ordinary(:,:,:,:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp) :: theta,stage_theta
    logical :: faces(6)
    character(len=16) :: patch_faces(6)
    integer :: nv,node,n(3),depth
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    noslip=.false.;wall_temp=0
    if(present(embedded_noslip)) noslip=embedded_noslip
    wall_active=noslip.or.present(embedded_wall_temperature)
    if(noslip.neqv.present(embedded_wall_impulse)) return
    if(noslip.and..not.viscosity) return
    if(wall_active.neqv.present(extended_geometry)) return
    if(wall_active.neqv.present(embedded_wall_energy)) return
    if(wall_active.neqv.present(embedded_frequency)) return
    if(present(embedded_wall_temperature)) then
      wall_temp=embedded_wall_temperature
      if(.not.conduction.or..not.ieee_is_finite(wall_temp).or.wall_temp<=0) return
    end if
    if(present(wall_temperature)) then
      if(.not.present(boundary_faces)) return
      if(.not.all(ieee_is_finite(wall_temperature)).or.any(wall_temperature<=0)) return
    end if
    if (.not.ieee_is_finite(dt).or.dt<=0._dp) return
    if (.not.valid_tree(tree,nv)) return
    selected_depth=-1
    if(present(target_depth)) then
      selected_depth=target_depth
      if(selected_depth<0.or.selected_depth>=tree%topology%level_count()) return
    end if
    if(present(ancestor_start).or.present(ancestor_end).or.present(ancestor_fraction)) then
      if(selected_depth<0) return
      if(.not.(present(ancestor_start).and.present(ancestor_end).and.present(ancestor_fraction))) return
    end if

    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if (.not.valid) return
    allocate(flux(size(tree%nodes)),rates%fields(size(tree%nodes)));theta=1._dp
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      allocate(rates%fields(node)%values(nv,n(1),n(2),n(3)),source=0._dp)
      if(selected_depth>=0.and.tree%topology%nodes(node)%level/=selected_depth) cycle
      allocate(ghost(nv,0:n(1)+1,0:n(2)+1,0:n(3)+1),gt(0:n(1)+1,0:n(2)+1,0:n(3)+1))
      call fill_reactive_eb_tree_halo_3d(species,tree,node,ghost,gt,valid, &
        inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,restrict_snapshot=.false., &
        ancestor_start=ancestor_start,ancestor_end=ancestor_end,ancestor_fraction=ancestor_fraction,frame=frame)
      if (.not.valid) return
      depth=tree%topology%nodes(node)%level+1
      ! Inactive physical domain faces remain closed to diffusion. Parent patch edges
      ! inside the domain are donor interfaces, even when no peer exists.
      faces=[origin(1,node)>0,origin(1,node)+n(1)<domain(1,depth), &
        origin(2,node)>0,origin(2,node)+n(2)<domain(2,depth), &
        origin(3,node)>0,origin(3,node)+n(3)<domain(3,depth)]
      if (present(inflow_mach)) faces=faces.or.(inflow_mach>=0)
      patch_faces='outflow'
      if(present(boundary_faces)) then
        where(.not.faces.and.(boundary_faces=='slip_wall'.or.boundary_faces=='no_slip_wall'.or. &
              boundary_faces=='isothermal_wall')) &
          patch_faces=boundary_faces
        faces=.true.
      end if
      ordinary=tree%nodes(node)%state
      call ordinary_transport_energy(frame,tree%nodes(node)%geometry,[1,1,1], &
        tree%nodes(node)%geometry%volume_fraction>0,ordinary,valid)
      if(.not.valid) return
      call ordinary_transport_energy(frame,tree%nodes(node)%geometry,[0,0,0],gt>0,ghost,valid)
      if(.not.valid) return
      call compute_sdc_eb_transport_3d(species,transport,ordinary,tree%nodes(node)%temperature, &
        tree%nodes(node)%geometry,dt,[viscosity,conduction,diffusion,barodiffusion],rates%fields(node)%values, &
        flux(node),stage_theta,valid,ghost,gt,faces,physical_faces=patch_faces,wall_temperature=wall_temperature)
      if (.not.valid) return
      theta=min(theta,stage_theta)
      deallocate(ghost,gt)
    end do
    if(.not.ieee_is_finite(theta)) return
    if(wall_active) then
      call compute_eb_tree_wall_transport_3d(species,transport,tree,extended_geometry, &
        present(embedded_wall_temperature),noslip,wall_temp,heat,wall_values,heat_frequency,valid, &
        boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,target_depth, &
        ancestor_start,ancestor_end,ancestor_fraction,frame)
      if(.not.valid) return
      heat_energy=dt*wall_values(5);wall_impulse=dt*wall_values(2:4)
      if(.not.all(ieee_is_finite([heat_energy,wall_impulse]))) return
      do node=1,size(tree%nodes)
        rates%fields(node)%values=rates%fields(node)%values+heat%fields(node)%values
        if(.not.all(ieee_is_finite(rates%fields(node)%values))) return
      end do
    end if
    transport_rate=rates;call move_alloc(flux,integrated_flux)
    minimum_theta=theta;ok=.true.
    if(wall_active) then
      embedded_wall_energy=heat_energy;embedded_frequency=heat_frequency
      if(noslip) embedded_wall_impulse=wall_impulse
    end if
  end subroutine compute_sdc_eb_tree_transport_3d

  subroutine advance_reactive_eb_tree_transport_3d(species,transport,dt,viscosity,conduction, &
      diffusion,barodiffusion,tree,boundary_change,minimum_theta,ok, &
      inflow_mach,inflow_temperature,inflow_pressure,boundary_faces,inflow_state,wall_temperature)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    logical, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), intent(in) :: dt
    type(eb_reactive_tree_3d), intent(inout) :: tree
    real(dp), intent(inout) :: boundary_change(:),minimum_theta
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: wall_temperature(6)
    logical, intent(out) :: ok
    type(eb_reactive_tree_3d) :: snapshot,candidate
    type(eb_integrated_flux_3d), allocatable :: flux(:)
    real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),boundary(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp) :: theta,stage_theta
    logical :: faces(6)
    character(len=16) :: patch_faces(6)
    integer :: nv,node,n(3),depth,c
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    if(present(wall_temperature)) then
      if(.not.present(boundary_faces)) return
      if(.not.all(ieee_is_finite(wall_temperature)).or.any(wall_temperature<=0)) return
    end if
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
      patch_faces='outflow'
      if(present(boundary_faces)) then
        where(.not.faces.and.(boundary_faces=='slip_wall'.or.boundary_faces=='no_slip_wall'.or. &
              boundary_faces=='isothermal_wall')) &
          patch_faces=boundary_faces
        faces=.true.
      end if
      call advance_reactive_eb_transport_patch_euler_3d(species,transport,snapshot%nodes(node)%geometry, &
        dt,viscosity,conduction,diffusion,barodiffusion,candidate%nodes(node)%state, &
        candidate%nodes(node)%temperature,flux(node),stage_theta,valid,ghost,gt,faces, &
            physical_faces=patch_faces,wall_temperature=wall_temperature)
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
  ! Provisional source-aware hydro: no restriction, redistribution or reflux.
  ! Rates are per fluid volume; fluxes are time-integrated face densities.
  subroutine compute_mol_eb_tree_hydro_3d(species,tree,extended_geometry,dt,order,hydro_rate,flux,ok, &
      boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
      target_depth,ancestor_start,ancestor_end,ancestor_fraction,lax_friedrichs,frame)
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: extended_geometry(:)
    real(dp), intent(in) :: dt
    integer, intent(in) :: order
    type(amr_tree_state_3d), intent(inout) :: hydro_rate
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: lax_friedrichs
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    integer, optional, intent(in) :: target_depth
    type(eb_reactive_tree_3d), optional, intent(in) :: ancestor_start,ancestor_end
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    type(amr_tree_state_3d) :: rates
    type(eb_integrated_flux_3d), allocatable :: candidate_flux(:)
    type(eb_geometry_3d) :: center
    real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),vf(:,:,:)
    real(dp) :: boundary(reactive_nvar(size(species))),wall(reactive_nvar(size(species)))
    integer :: nv,node,n(3),depth
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species));depth=-1
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    if(.not.valid_tree(tree,nv).or..not.ieee_is_finite(dt)) return
    if(dt<=0.or.(order/=1.and.order/=2)) return
    if(size(extended_geometry)/=size(tree%nodes)) return
    if(present(target_depth)) then
      depth=target_depth
      if(depth<0.or.depth>=tree%topology%level_count()) return
    end if
    if(present(ancestor_start).or.present(ancestor_end).or.present(ancestor_fraction)) then
      if(depth<0) return
    end if
    allocate(rates%fields(size(tree%nodes)),candidate_flux(size(tree%nodes)))
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      if(any([extended_geometry(node)%nx,extended_geometry(node)%ny,extended_geometry(node)%nz]/=n+6)) return
      call extract_eb_geometry_3d(extended_geometry(node),[4,4,4],n+3,center,valid)
      if(.not.valid) return
      if(.not.matching_eb_metrics_3d(center,tree%nodes(node)%geometry)) return
      allocate(rates%fields(node)%values(nv,n(1),n(2),n(3)),source=0._dp)
      if(depth>=0.and.tree%topology%nodes(node)%level/=depth) cycle
      allocate(ghost(nv,-2:n(1)+3,-2:n(2)+3,-2:n(3)+3), &
        gt(-2:n(1)+3,-2:n(2)+3,-2:n(3)+3),vf(-2:n(1)+3,-2:n(2)+3,-2:n(3)+3))
      call fill_reactive_eb_tree_halo_3d(species,tree,node,ghost,gt,valid, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_mach=inflow_mach, &
        inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,ghost_volume_fraction=vf, &
        halo_width=3,target_volume_fraction=extended_geometry(node)%volume_fraction,restrict_snapshot=.false., &
        ancestor_start=ancestor_start,ancestor_end=ancestor_end,ancestor_fraction=ancestor_fraction,frame=frame)
      if(.not.valid) return
      call compute_mol_eb_hydro_3d(species,extended_geometry(node),ghost,gt,dt,order, &
        rates%fields(node)%values,candidate_flux(node),boundary,wall,valid,lax_friedrichs,frame)
      if(.not.valid) return
      deallocate(ghost,gt,vf)
    end do
    call move_alloc(rates%fields,hydro_rate%fields);call move_alloc(candidate_flux,flux);ok=.true.
  end subroutine compute_mol_eb_tree_hydro_3d

  subroutine compute_sdc_eb_tree_hydro_3d(species,tree,predictor,extended_geometry,dt,hydro_rate,flux,ok, &
      difmag,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
      target_depth,ancestor_start,ancestor_end,ancestor_fraction,source_start,source_end,source_fraction)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: predictor
    type(eb_geometry_3d), intent(in) :: extended_geometry(:)
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: hydro_rate
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: target_depth
    type(eb_reactive_tree_3d), optional, intent(in) :: ancestor_start,ancestor_end
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    type(amr_tree_state_3d), optional, intent(in) :: source_start,source_end
    real(dp), optional, intent(in) :: source_fraction(:)
    integer :: selected_depth
    real(dp), optional, intent(in) :: difmag,inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(amr_tree_state_3d) :: rates
    type(eb_integrated_flux_3d), allocatable :: candidate_flux(:)
    type(eb_geometry_3d) :: working,center
    type(eb_plm_boundary_3d) :: physical_boundary
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),vf(:,:,:),flat(:,:,:),source(:,:,:,:)
    real(dp) :: coefficient,boundary(reactive_nvar(size(species))),wall(reactive_nvar(size(species)))
    character(len=32) :: faces(6)
    integer :: nv,node,n(3),level
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    if(.not.ieee_is_finite(dt).or.dt<=0) return
    if(.not.valid_tree(tree,nv)) return
    selected_depth=-1
    if(present(target_depth)) then
      selected_depth=target_depth
      if(selected_depth<0.or.selected_depth>=tree%topology%level_count()) return
    end if
    if(present(ancestor_start).or.present(ancestor_end).or.present(ancestor_fraction)) then
      if(selected_depth<0) return
      if(.not.(present(ancestor_start).and.present(ancestor_end).and.present(ancestor_fraction))) return
    end if
    if(present(source_start).or.present(source_end).or.present(source_fraction)) then
      if(selected_depth<0) return
      if(.not.(present(source_start).and.present(source_end).and.present(source_fraction))) return
    end if

    if(size(extended_geometry)/=size(tree%nodes)) return
    if(.not.allocated(predictor%fields)) return
    if(size(predictor%fields)/=size(tree%nodes)) return
    coefficient=.1_dp
    if(present(difmag)) coefficient=difmag
    if(.not.ieee_is_finite(coefficient).or.coefficient<0) return
    faces='outflow'
    if(present(boundary_faces)) faces=boundary_faces
    if(present(inflow_mach)) then
      where(inflow_mach>=0) faces='inflow'
    end if
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if(.not.valid) return
    if(any(domain>int(huge(1)-4,int64))) return
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      if(.not.allocated(predictor%fields(node)%values)) return
      if(any(shape(predictor%fields(node)%values)/=[nv,n])) return
      if(any([extended_geometry(node)%nx,extended_geometry(node)%ny,extended_geometry(node)%nz]/=n+12)) return
      call extract_eb_geometry_3d(extended_geometry(node),[7,7,7],n+6,center,valid)
      if(.not.valid) return
      if(.not.matching_eb_metrics_3d(center,tree%nodes(node)%geometry)) return
    end do
    allocate(rates%fields(size(tree%nodes)),candidate_flux(size(tree%nodes)))
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      allocate(rates%fields(node)%values(nv,n(1),n(2),n(3)),source=0._dp)
      if(selected_depth>=0.and.tree%topology%nodes(node)%level/=selected_depth) cycle
      allocate(ghost(nv,-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),gt(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), &
        vf(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),flat(0:n(1)+7,0:n(2)+7,0:n(3)+7), &
        source(nv,n(1)+6,n(2)+6,n(3)+6))
      call fill_reactive_eb_tree_halo_3d(species,tree,node,ghost,gt,valid, &
        inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,ghost_volume_fraction=vf,halo_width=6, &
        target_volume_fraction=extended_geometry(node)%volume_fraction,restrict_snapshot=.false., &
        ancestor_start=ancestor_start,ancestor_end=ancestor_end,ancestor_fraction=ancestor_fraction)
      if(.not.valid) return
      call extract_eb_geometry_3d(extended_geometry(node),[4,4,4],n+9,working,valid)
      if(.not.valid) return
      call fill_sdc_eb_tree_source_halo_3d(tree,predictor,node,3,source,valid, &
        boundary_faces=faces,target_volume_fraction=working%volume_fraction, &
        ancestor_start=source_start,ancestor_end=source_end,ancestor_fraction=source_fraction)
      if(.not.valid) return
      call compute_reactive_eb_flattening_3d(species,ghost,gt,vf,n+6,flat,valid)
      if(.not.valid) return
      level=tree%topology%nodes(node)%level+1
      physical_boundary%lower=4-int(origin(:,node))
      physical_boundary%upper=3+int(domain(:,level)-origin(:,node))
      physical_boundary%faces=faces
      call compute_sdc_eb_hydro_3d(species,working,ghost,gt,vf,dt,flat,source,rates%fields(node)%values, &
        candidate_flux(node),boundary,wall,valid,difmag=coefficient,boundary=physical_boundary)
      if(.not.valid) return
      deallocate(ghost,gt,vf,flat,source)
    end do
    hydro_rate=rates;call move_alloc(candidate_flux,flux);ok=.true.
  end subroutine compute_sdc_eb_tree_hydro_3d

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
      inflow_mach,inflow_temperature,inflow_pressure,boundary_faces,inflow_state,wall_temperature)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: dt
    logical, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    type(eb_reactive_tree_3d), intent(inout) :: tree
    real(dp), intent(inout) :: boundary_change(:),minimum_theta
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: wall_temperature(6)
    logical, intent(out) :: ok
    real(dp) :: wall(size(boundary_change))
    call advance_reactive_eb_tree_subcycle_operator_3d(species,'',dt,tree,boundary_change,wall,ok, &
      transport,viscosity,conduction,diffusion,barodiffusion,minimum_theta, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
      boundary_faces=boundary_faces,inflow_state=inflow_state,wall_temperature=wall_temperature)
  end subroutine advance_reactive_eb_tree_subcycle_transport_3d

  subroutine advance_reactive_eb_tree_subcycle_operator_3d(species,solver,dt,tree,boundary_change,wall_change,ok, &
      transport,viscosity,conduction,diffusion,barodiffusion,minimum_theta, &
      inflow_mach,inflow_temperature,inflow_pressure,difmag,extended_geometry,boundary_faces,inflow_state,wall_temperature)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), optional, intent(in) :: extended_geometry(:)
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    type(eb_reactive_tree_3d), intent(inout) :: tree
    real(dp), intent(inout) :: boundary_change(:),wall_change(:)
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: wall_temperature(6)
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
    if(present(wall_temperature)) then
      if(.not.present(boundary_faces)) return
      if(.not.all(ieee_is_finite(wall_temperature)).or.any(wall_temperature<=0)) return
    end if
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
      character(len=16) :: patch_faces(6)
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
          patch_faces='outflow'
          if(present(boundary_faces)) then
            where(.not.faces.and.(boundary_faces=='slip_wall'.or.boundary_faces=='no_slip_wall'.or. &
              boundary_faces=='isothermal_wall')) &
              patch_faces=boundary_faces
            faces=.true.
          end if
          call advance_reactive_eb_transport_patch_euler_3d(species,transport,snapshot%nodes(node)%geometry, &
            interval,viscosity,conduction,diffusion,barodiffusion,candidate%nodes(node)%state, &
            candidate%nodes(node)%temperature,flux(node),stage_theta,valid,ghost,gt,faces, &
            physical_faces=patch_faces,wall_temperature=wall_temperature)
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
  subroutine synchronize_reactive_eb_tree_flux_3d(species,tree,flux,ok,minimum_level,maximum_level,boundary_faces,frame)
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
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
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    if (.not.valid_tree(tree,nv)) return
    faces='outflow'
    if(present(boundary_faces)) faces=boundary_faces
    if(any(faces/='outflow'.and.faces/='inflow'.and.faces/='periodic'.and. &
      faces/='slip_wall'.and.faces/='no_slip_wall'.and.faces/='isothermal_wall')) return
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
        candidate%nodes(node)%state,candidate%nodes(node)%temperature,valid,all_faces=.true.,frame=frame)
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
    call average_down_reactive_eb_tree_3d(species,candidate,valid,frame)
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
  ! Signed SDC sources use PCM ancestors and scalar extension: mirrored donor
  ! positions at walls, clamped elsewhere. No EOS, momentum sign flip or restriction.
  subroutine fill_sdc_eb_tree_source_halo_3d(tree,rates,target,halo_width,ghost,ok, &
      boundary_faces,target_volume_fraction,ancestor_start,ancestor_end,ancestor_fraction)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: rates
    integer, intent(in) :: target,halo_width
    real(dp), intent(inout) :: ghost(:,:,:,:)
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: target_volume_fraction(:,:,:)
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    real(dp), allocatable :: candidate(:,:,:,:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer(int64) :: global(3),mapped(3),local(3)
    integer :: boundary_face
    integer :: n(3),index(3),source(3),nv,width,level,depth,node,donor,i,j,k,axis
    real(dp) :: fraction,volume
    character(len=32) :: faces(6)
    logical :: valid
    ok=.false.;width=halo_width;nv=size(ghost,1)
    if(nv<1.or.width<1.or.width>6) return
    if(.not.tree%topology%is_valid().or..not.allocated(tree%nodes)) return
    if(size(tree%nodes)/=size(tree%topology%nodes)) return
    if(target<1.or.target>size(tree%nodes)) return
    if(.not.allocated(tree%nodes(1)%state)) return
    if(.not.valid_tree(tree,size(tree%nodes(1)%state,1))) return
    if(.not.field_shape(rates)) return
    faces='outflow'
    if(present(boundary_faces)) faces=boundary_faces
    if(.not.valid_reactive_boundaries_3d(faces)) return
    if(.not.periodic_eb_faces_match(tree%nodes(1)%geometry,faces)) return
    n=tree%topology%nodes(target)%shape;level=tree%topology%nodes(target)%level+1
    if(any(n>huge(1)-2*width)) return
    if(any(shape(ghost)/=[nv,n+2*width])) return
    if(present(target_volume_fraction)) then
      if(any(shape(target_volume_fraction)/=n+2*width)) return
      if(.not.all(ieee_is_finite(target_volume_fraction))) return
      if(any(target_volume_fraction<0).or.any(target_volume_fraction>1)) return
      if(any(target_volume_fraction(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))/= &
        tree%nodes(target)%geometry%volume_fraction)) return
    end if
    if(present(ancestor_start).or.present(ancestor_end).or.present(ancestor_fraction)) then
      if(.not.(present(ancestor_start).and.present(ancestor_end).and.present(ancestor_fraction))) return
      if(.not.field_shape(ancestor_start).or..not.field_shape(ancestor_end)) return
      if(size(ancestor_fraction)/=level-1) return
      if(.not.all(ieee_is_finite(ancestor_fraction))) return
      if(any(ancestor_fraction<0).or.any(ancestor_fraction>1)) return
    end if
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if(.not.valid) return
    allocate(candidate(nv,n(1)+2*width,n(2)+2*width,n(3)+2*width),source=0._dp)
    do k=1,n(3)+2*width
    do j=1,n(2)+2*width
    do i=1,n(1)+2*width
      index=[i,j,k]-width
      global=origin(:,target)+int(index-1,int64)
      do axis=1,3
        if(faces(2*axis-1)=='periodic') global(axis)=modulo(global(axis),domain(axis,level))
      end do
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
      volume=tree%nodes(donor)%geometry%volume_fraction(source(1),source(2),source(3))
      if(present(target_volume_fraction)) then
        if(target_volume_fraction(i,j,k)<=0) cycle
        if(volume<=0) return
      else
        if(volume<=0) cycle
      end if
      depth=tree%topology%nodes(donor)%level+1
      if(present(ancestor_fraction).and.depth<level) then
        fraction=ancestor_fraction(depth)
        if(fraction==0) then
          candidate(:,i,j,k)=ancestor_start%fields(donor)%values(:,source(1),source(2),source(3))
        else if(fraction==1) then
          candidate(:,i,j,k)=ancestor_end%fields(donor)%values(:,source(1),source(2),source(3))
        else
          candidate(:,i,j,k)=(1-fraction)*ancestor_start%fields(donor)%values(:,source(1),source(2),source(3))+ &
            fraction*ancestor_end%fields(donor)%values(:,source(1),source(2),source(3))
        end if
      else
        candidate(:,i,j,k)=rates%fields(donor)%values(:,source(1),source(2),source(3))
      end if
      if(.not.all(ieee_is_finite(candidate(:,i,j,k)))) return
    end do
    end do
    end do
    ghost=candidate;ok=.true.
  contains
    logical function field_shape(field) result(good)
      type(amr_tree_state_3d), intent(in) :: field
      integer :: p
      good=.false.
      if(.not.allocated(field%fields)) return
      if(size(field%fields)/=size(tree%nodes)) return
      do p=1,size(field%fields)
        if(.not.allocated(field%fields(p)%values)) return
        if(any(shape(field%fields(p)%values)/=[nv,tree%topology%nodes(p)%shape])) return
      end do
      good=.true.
    end function field_shape
  end subroutine fill_sdc_eb_tree_source_halo_3d

  subroutine fill_reactive_eb_tree_halo_3d(species,tree,target,ghost,ghost_temperature,ok, &
      ancestor_start,ancestor_end,ancestor_fraction,inflow_mach,inflow_temperature,inflow_pressure, &
      ghost_volume_fraction,halo_width,target_volume_fraction,boundary_faces,inflow_state,restrict_snapshot,frame)
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
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
    logical, optional, intent(in) :: restrict_snapshot
    logical :: valid,physical,restrict_state
    ok=.false.;nv=reactive_nvar(size(species));width=1
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    if(present(halo_width)) width=halo_width
    if(width<1.or.width>huge(1)/4) return
    if(width>1.and.(.not.present(ghost_volume_fraction).or..not.present(target_volume_fraction))) return
    faces='outflow'
    if(present(boundary_faces)) then
      if(any(boundary_faces/='inflow'.and.boundary_faces/='outflow'.and.boundary_faces/='periodic'.and. &
        boundary_faces/='slip_wall'.and.boundary_faces/='no_slip_wall'.and.boundary_faces/='isothermal_wall')) return
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
    restrict_state=.true.
    if(present(restrict_snapshot)) restrict_state=restrict_snapshot
    if(restrict_state) then
      call average_down_reactive_eb_tree_3d(species,snapshot,valid,frame)
      if (.not.valid) return
    end if
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
        call recover_reactive_eb_cell_3d(species,tree%nodes(target)%geometry,[i,j,k],value, &
          guess,primitive,recovered,speed,valid,frame)
        if (.not.valid) return
        if (physical.and.mach(2*axis-1+side)>=0) then
          call boundary_energy(.false.,valid)
          if(.not.valid) return
          if (pin(2*axis-1+side)>0) then
            call mach_inlet_ghost_3d(species,value,recovered,2*axis-1+side,tin(2*axis-1+side), &
              mach(2*axis-1+side),inlet_value,inlet_t,valid,inlet_pressure=pin(2*axis-1+side))
          else
            call mach_inlet_ghost_3d(species,value,recovered,2*axis-1+side,tin(2*axis-1+side), &
              mach(2*axis-1+side),inlet_value,inlet_t,valid)
          end if
          if (.not.valid) return
          value=inlet_value;recovered=inlet_t
          call boundary_energy(.true.,valid)
          if(.not.valid) return
        end if
        if(physical) then
          call prescribed_inlet(2*axis-1+side,valid)
          if(.not.valid) return
          if(faces(2*axis-1+side)=='slip_wall') value(1+axis)=-value(1+axis)
          if(faces(2*axis-1+side)=='no_slip_wall'.or.faces(2*axis-1+side)=='isothermal_wall') value(2:4)=-value(2:4)
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
              faces(2*axis)=='no_slip_wall'.or.any(faces(2*axis-1:2*axis)=='isothermal_wall')) local(axis)=index(axis)
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
        call recover_reactive_eb_cell_3d(species,tree%nodes(target)%geometry,[i,j,k],value,guess,primitive, &
      recovered,speed,valid,frame)
        if(.not.valid) return
        if(inlet_face>0) then
          if(faces(inlet_face)=='inflow') then
            call prescribed_inlet(inlet_face,valid)
            if(.not.valid) return
          else
          call boundary_energy(.false.,valid)
          if(.not.valid) return
          if(pin(inlet_face)>0) then
            call mach_inlet_ghost_3d(species,value,recovered,inlet_face,tin(inlet_face),mach(inlet_face), &
              inlet_value,inlet_t,valid,inlet_pressure=pin(inlet_face))
          else
            call mach_inlet_ghost_3d(species,value,recovered,inlet_face,tin(inlet_face),mach(inlet_face), &
              inlet_value,inlet_t,valid)
          end if
          if(.not.valid) return
          value=inlet_value;recovered=inlet_t
          call boundary_energy(.true.,valid)
          if(.not.valid) return
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
        candidate(:,0:n(1)+1,0:n(2)+1,0:n(3)+1),temps(0:n(1)+1,0:n(2)+1,0:n(3)+1),valid,frame=frame)
      if (.not.valid) return
    end do
    ghost=candidate;ghost_temperature=temps
    if(present(ghost_volume_fraction)) ghost_volume_fraction=volumes
    ok=.true.
  contains
    ! Prescribed and Mach inlet APIs use ordinary relative total energy.
    subroutine boundary_energy(to_rotating,accepted)
      logical, intent(in) :: to_rotating
      logical, intent(out) :: accepted
      real(dp) :: position(3),energy
      accepted=.true.
      if(.not.present(frame)) return
      position=[tree%nodes(target)%geometry%x_lower,tree%nodes(target)%geometry%y_lower, &
        tree%nodes(target)%geometry%z_lower]+(real([i,j,k],dp)-.5_dp)* &
        [tree%nodes(target)%geometry%dx,tree%nodes(target)%geometry%dy,tree%nodes(target)%geometry%dz]
      energy=value(5)
      call convert_rotating_energy_3d(frame,position,value(1),value(5),to_rotating,energy,accepted)
      if(accepted) value(5)=energy
    end subroutine boundary_energy
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
          if(faces(face)/='slip_wall'.and.faces(face)/='no_slip_wall'.and.faces(face)/='isothermal_wall') exit
          if(faces(face)=='no_slip_wall'.or.faces(face)=='isothermal_wall') then
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
      if(accepted) then
        value=inflow_state(:,face)
        call boundary_energy(.true.,accepted)
      end if
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

  ! Provisional SDC chemistry evolves all fluid cells at each selected level,
  ! including refined parents. Reflux and restriction belong to the time driver.
  subroutine react_sdc_eb_tree_3d(species,reactions,old_tree,nonreacting_rates,dt,rtol,atol, &
      state,temperature,reaction_source,ok,target_depth, &
      evaluate_reaction,frame)
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), allocatable :: potential(:,:,:)
    real(dp) :: position(3)
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(eb_reactive_tree_3d), intent(in) :: old_tree
    type(amr_tree_state_3d), intent(in) :: nonreacting_rates
    real(dp), intent(in) :: dt,rtol,atol
    type(amr_tree_state_3d), intent(inout) :: state,temperature,reaction_source
    logical, intent(out) :: ok
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d) :: candidate,temps,sources
    integer :: nv,node,n(3),depth,i,j,k
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species));depth=-1
    if(present(frame)) then
      if(.not.frame%is_valid()) return
      allocate(potential(1,1,1))
    end if
    if(.not.valid_tree(old_tree,nv).or.size(reactions)<1) return
    if(.not.field_shape(nonreacting_rates,nv).or..not.field_shape(reaction_source,nv-4)) return
    if(.not.all(ieee_is_finite([dt,rtol,atol])).or.min(dt,rtol,atol)<=0) return
    if(present(target_depth)) then
      depth=target_depth
      if(depth<0.or.depth>=old_tree%topology%level_count()) return
      if(.not.field_shape(state,nv).or..not.field_shape(temperature,1)) return
      candidate=state;temps=temperature
    else
      allocate(candidate%fields(size(old_tree%nodes)),temps%fields(size(old_tree%nodes)))
      do node=1,size(old_tree%nodes)
        n=old_tree%topology%nodes(node)%shape
        candidate%fields(node)%values=old_tree%nodes(node)%state
        allocate(temps%fields(node)%values(1,n(1),n(2),n(3)))
        temps%fields(node)%values(1,:,:,:)=old_tree%nodes(node)%temperature
      end do
    end if
    sources=reaction_source
    do node=1,size(old_tree%nodes)
      if(depth>=0.and.old_tree%topology%nodes(node)%level/=depth) cycle
      n=old_tree%topology%nodes(node)%shape
      candidate%fields(node)%values=old_tree%nodes(node)%state
      temps%fields(node)%values(1,:,:,:)=old_tree%nodes(node)%temperature
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(old_tree%nodes(node)%geometry%volume_fraction(i,j,k)<=0) cycle
        if(present(frame)) then
          associate(g=>old_tree%nodes(node)%geometry)
          position=[g%x_lower,g%y_lower,g%z_lower]+(real([i,j,k],dp)-.5_dp)*[g%dx,g%dy,g%dz]
          end associate
          call rotating_potential_3d(frame,position,potential(1,1,1),valid)
          if(.not.valid) return
        end if
        call react_sdc_state_3d(species,reactions,old_tree%nodes(node)%state(:,i:i,j:j,k:k), &
          old_tree%nodes(node)%temperature(i:i,j:j,k:k),nonreacting_rates%fields(node)%values(:,i:i,j:j,k:k), &
          dt,rtol,atol,candidate%fields(node)%values(:,i:i,j:j,k:k),temps%fields(node)%values(1,i:i,j:j,k:k), &
          sources%fields(node)%values(:,i:i,j:j,k:k),valid, &
          evaluate_reaction=evaluate_reaction,rotating_potential=potential)
        if(.not.valid) return
      end do
      end do
      end do
    end do
    state=candidate;temperature=temps;reaction_source=sources;ok=.true.
  contains
    logical function field_shape(field,components) result(good)
      type(amr_tree_state_3d), intent(in) :: field
      integer, intent(in) :: components
      integer :: p
      good=.false.
      if(.not.allocated(field%fields)) return
      if(size(field%fields)/=size(old_tree%nodes)) return
      do p=1,size(field%fields)
        if(.not.allocated(field%fields(p)%values)) return
        if(any(shape(field%fields(p)%values)/=[components,old_tree%topology%nodes(p)%shape])) return
      end do
      good=.true.
    end function field_shape
  end subroutine react_sdc_eb_tree_3d

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
