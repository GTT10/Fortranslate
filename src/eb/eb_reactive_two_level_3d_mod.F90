module eb_reactive_two_level_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use eb_reactive_transport_3d_mod, only: advance_reactive_eb_transport_patch_euler_3d
  use elementary_kinetics_mod, only: elementary_reaction
  use reactive_3d_mod, only: advance_reactive_chemistry_3d
  use reactive_1d_mod, only: reactive_nvar
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_transfer_3d_mod, only: restrict_reactive_eb_state_3d,fill_reactive_eb_patch_halo_3d, &
    overlay_reactive_eb_peer_halo_3d,prolong_reactive_eb_state_pcm_3d
  use eb_reactive_hydro_3d_mod, only: advance_reactive_eb_euler_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d,restrict_eb_integrated_flux_3d,reflux_reactive_eb_state_3d, &
    reconcile_reactive_eb_peer_flux_3d
  use eb_reactive_wall_flux_3d_mod, only: reactive_eb_slip_wall_source_3d
  implicit none
  private
  type, public :: eb_reactive_patch_3d
    type(eb_geometry_3d) :: geometry
    real(dp), allocatable :: state(:,:,:,:),temperature(:,:,:)
  end type eb_reactive_patch_3d
  public :: tag_reactive_eb_temperature_3d
  public :: remap_reactive_eb_multipatch_3d
  ! Multipatch stepping accepts empty patch arrays: only the coarse grid
  ! advances, with its physical boundary/wall and chemistry ledgers retained.
  public :: advance_reactive_eb_multipatch_full_3d
  public :: advance_reactive_eb_multipatch_euler_3d,advance_reactive_eb_multipatch_transport_3d
  public :: advance_reactive_eb_two_level_subcycled_transport_3d
  public :: advance_reactive_eb_two_level_subcycled_euler_3d
  public :: advance_reactive_eb_two_level_euler_3d,advance_reactive_eb_two_level_strang_3d
  public :: advance_reactive_eb_two_level_transport_3d,advance_reactive_eb_two_level_full_3d
contains
  ! Deterministic coarse-index tiles selected by hot exposed fluid cells.
  ! Parent values hidden by fine patches and solid-cell values do not vote.
  subroutine tag_reactive_eb_temperature_3d(coarse,coarse_temperature,lower,ratio,patches, &
      threshold,tile_shape,buffer,maximum_patches,tag_lower,tag_upper,ok)
    type(eb_geometry_3d), intent(in) :: coarse
    real(dp), intent(in) :: coarse_temperature(:,:,:),threshold
    integer, intent(in) :: lower(:,:),ratio,tile_shape(3),buffer,maximum_patches
    type(eb_reactive_patch_3d), intent(in) :: patches(:)
    integer, allocatable, intent(inout) :: tag_lower(:,:),tag_upper(:,:)
    logical, intent(out) :: ok
    type(eb_geometry_3d) :: section,restricted
    logical, allocatable :: refined(:,:,:),hot(:,:,:),padded(:,:,:)
    integer, allocatable :: candidate_lower(:,:),candidate_upper(:,:)
    integer :: n(3),d(3),upper(3),index(3),lo(3),hi(3),tile(3),pad(3)
    integer :: p,i,j,k,count,pass,allocation_status
    real(dp) :: temperature
    logical :: valid
    ok=.false.
    if (.not.coarse%is_valid()) return
    if (.not.ieee_is_finite(threshold)) return
    if (threshold<=0._dp.or.ratio<1.or.any(tile_shape<1).or.buffer<0.or.maximum_patches<1) return
    n=[coarse%nx,coarse%ny,coarse%nz]
    if (any(shape(coarse_temperature)/=n).or.any(shape(lower)/=[3,size(patches)])) return
    allocate(refined(n(1),n(2),n(3)),hot(n(1),n(2),n(3)),padded(n(1),n(2),n(3)),stat=allocation_status)
    if (allocation_status/=0) return
    refined=.false.;hot=.false.;padded=.false.
    do p=1,size(patches)
      if (.not.patches(p)%geometry%is_valid().or..not.allocated(patches(p)%temperature)) return
      d=[patches(p)%geometry%nx,patches(p)%geometry%ny,patches(p)%geometry%nz]
      if (any(shape(patches(p)%temperature)/=d)) return
      if (any(mod(d,ratio)/=0).or.any(lower(:,p)<1).or.any(lower(:,p)>n)) return
      if (any(d/ratio>n-lower(:,p)+1)) return
      upper=lower(:,p)+d/ratio-1
      if (any(refined(lower(1,p):upper(1),lower(2,p):upper(2),lower(3,p):upper(3)))) return
      call coarsen_eb_geometry_3d(patches(p)%geometry,ratio,restricted,valid)
      if (.not.valid) return
      call extract_eb_geometry_3d(coarse,lower(:,p),upper,section,valid)
      if (.not.valid) return
      if (.not.matching_metrics(section,restricted)) return
      refined(lower(1,p):upper(1),lower(2,p):upper(2),lower(3,p):upper(3))=.true.
      do k=1,d(3)
      do j=1,d(2)
      do i=1,d(1)
        if (patches(p)%geometry%volume_fraction(i,j,k)<=0._dp) cycle
        temperature=patches(p)%temperature(i,j,k)
        if (.not.ieee_is_finite(temperature)) return
        if (temperature<=0._dp) return
        if (temperature<threshold) cycle
        index=lower(:,p)+([i,j,k]-1)/ratio
        hot(index(1),index(2),index(3))=.true.
      end do
      end do
      end do
    end do
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      if (refined(i,j,k).or.coarse%volume_fraction(i,j,k)<=0._dp) cycle
      temperature=coarse_temperature(i,j,k)
      if (.not.ieee_is_finite(temperature)) return
      if (temperature<=0._dp) return
      if (temperature>=threshold) hot(i,j,k)=.true.
    end do
    end do
    end do
    ! Clip before adding the halo width, including arbitrarily large buffers.
    pad=min(buffer,n)
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      if (.not.hot(i,j,k)) cycle
      index=[i,j,k];lo=index-min(pad,index-1);hi=index+min(pad,n-index)
      padded(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.true.
    end do
    end do
    end do
    tile=min(tile_shape,n)
    do pass=1,2
      count=0
      do k=1,n(3),tile(3)
      do j=1,n(2),tile(2)
      do i=1,n(1),tile(1)
        lo=[i,j,k];hi=lo+min(tile-1,n-lo)
        if (.not.any(padded(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))) cycle
        ! A buffered tile made entirely of solid cells needs no fine patch.
        if (.not.any(coarse%volume_fraction(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))>0._dp)) cycle
        if (count==maximum_patches) return
        count=count+1
        if (pass==2) then
          candidate_lower(:,count)=lo;candidate_upper(:,count)=hi
        end if
      end do
      end do
      end do
      if (pass==1) then
        allocate(candidate_lower(3,count),candidate_upper(3,count),stat=allocation_status)
        if (allocation_status/=0) return
      end if
    end do
    call move_alloc(candidate_lower,tag_lower);call move_alloc(candidate_upper,tag_upper)
    ok=.true.
  end subroutine tag_reactive_eb_temperature_3d

  ! Fixed-geometry conservative regridding. New patch geometries are supplied
  ! by the caller. Old/new arrays must be distinct; all output publishes once.
  subroutine remap_reactive_eb_multipatch_3d(species,coarse,old_lower,old_ratio,old_patches, &
      new_lower,new_ratio,new_patches,coarse_state,coarse_temperature,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: coarse
    integer, intent(in) :: old_lower(:,:),new_lower(:,:),old_ratio,new_ratio
    type(eb_reactive_patch_3d), intent(in) :: old_patches(:)
    type(eb_reactive_patch_3d), intent(inout) :: new_patches(:)
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:)
    logical, intent(out) :: ok
    type(eb_reactive_patch_3d), allocatable :: candidate(:)
    type(eb_geometry_3d) :: old_overlap,new_overlap
    real(dp), allocatable :: u(:,:,:,:),t(:,:,:),checked(:,:,:,:),checked_t(:,:,:)
    integer :: n(3),nv,p,q,i,j,k,parent(3),lo(3),hi(3),ol(3),oh(3),nl(3),nh(3),dims(3)
    integer :: old_upper(3,size(old_patches)),new_upper(3,size(new_patches))
    logical :: valid
    ok=.false.
    if (.not.coarse%is_valid()) return
    n=[coarse%nx,coarse%ny,coarse%nz];nv=reactive_nvar(size(species))
    if (any(shape(coarse_state)/=[nv,n]).or.any(shape(coarse_temperature)/=n)) return
    call validate_layout(old_patches,old_lower,old_ratio,old_upper,valid)
    if (.not.valid) return
    call validate_layout(new_patches,new_lower,new_ratio,new_upper,valid)
    if (.not.valid) return
    u=coarse_state;t=coarse_temperature
    do p=1,size(old_patches)
      if (.not.allocated(old_patches(p)%state).or..not.allocated(old_patches(p)%temperature)) return
      lo=old_lower(:,p);hi=old_upper(:,p)
      call restrict_reactive_eb_state_3d(species,old_patches(p)%geometry,old_ratio, &
        old_patches(p)%state,old_patches(p)%temperature,u(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
        t(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),valid)
      if (.not.valid) return
    end do
    ! Validate exposed parent cells as well as restricted fine cells, and
    ! recover temperatures from the conserved donor state without aliasing.
    checked=u;checked_t=t
    call prolong_reactive_eb_state_pcm_3d(species,coarse,1,u,t,checked,checked_t,valid)
    if (.not.valid) return
    u=checked;t=checked_t
    candidate=new_patches
    do p=1,size(candidate)
      dims=[candidate(p)%geometry%nx,candidate(p)%geometry%ny,candidate(p)%geometry%nz]
      if (allocated(candidate(p)%state)) deallocate(candidate(p)%state)
      if (allocated(candidate(p)%temperature)) deallocate(candidate(p)%temperature)
      allocate(candidate(p)%state(nv,dims(1),dims(2),dims(3)), &
        candidate(p)%temperature(dims(1),dims(2),dims(3)))
      ! Preserve meaningful covered-cell payload as well as fluid state.
      ! Prolongation updates fluid children only; checkpoints also retain the
      ! opaque solid-cell values and require finite, positive storage there.
      do k=1,dims(3)
      do j=1,dims(2)
      do i=1,dims(1)
        parent=new_lower(:,p)+([i,j,k]-1)/new_ratio
        candidate(p)%state(:,i,j,k)=u(:,parent(1),parent(2),parent(3))
        candidate(p)%temperature(i,j,k)=t(parent(1),parent(2),parent(3))
      end do
      end do
      end do
      lo=new_lower(:,p);hi=new_upper(:,p)
      call prolong_reactive_eb_state_pcm_3d(species,candidate(p)%geometry,new_ratio, &
        u(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),t(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
        candidate(p)%state,candidate(p)%temperature,valid)
      if (.not.valid) return
      if (old_ratio/=new_ratio) cycle
      do q=1,size(old_patches)
        lo=max(new_lower(:,p),old_lower(:,q));hi=min(new_upper(:,p),old_upper(:,q))
        if (any(lo>hi)) cycle
        nl=(lo-new_lower(:,p))*new_ratio+1;nh=(hi-new_lower(:,p)+1)*new_ratio
        ol=(lo-old_lower(:,q))*old_ratio+1;oh=(hi-old_lower(:,q)+1)*old_ratio
        call extract_eb_geometry_3d(candidate(p)%geometry,nl,nh,new_overlap,valid)
        if (.not.valid) return
        call extract_eb_geometry_3d(old_patches(q)%geometry,ol,oh,old_overlap,valid)
        if (.not.valid) return
        if (.not.matching_metrics(old_overlap,new_overlap)) return
        candidate(p)%state(:,nl(1):nh(1),nl(2):nh(2),nl(3):nh(3))= &
          old_patches(q)%state(:,ol(1):oh(1),ol(2):oh(2),ol(3):oh(3))
        candidate(p)%temperature(nl(1):nh(1),nl(2):nh(2),nl(3):nh(3))= &
          old_patches(q)%temperature(ol(1):oh(1),ol(2):oh(2),ol(3):oh(3))
      end do
    end do
    coarse_state=u;coarse_temperature=t;new_patches=candidate
    ok=.true.
  contains
    subroutine validate_layout(patches,lower,ratio,upper,valid)
      type(eb_reactive_patch_3d), intent(in) :: patches(:)
      integer, intent(in) :: lower(:,:),ratio
      integer, intent(out) :: upper(:,:)
      logical, intent(out) :: valid
      type(eb_geometry_3d) :: section,restricted
      logical :: good
      integer :: a,b,d(3)
      valid=.false.
      if (ratio<1.or.any(shape(lower)/=[3,size(patches)])) return
      do a=1,size(patches)
        if (.not.patches(a)%geometry%is_valid()) return
        d=[patches(a)%geometry%nx,patches(a)%geometry%ny,patches(a)%geometry%nz]
        if (any(mod(d,ratio)/=0).or.any(lower(:,a)<1).or.any(lower(:,a)>n)) return
        if (any(d/ratio>n-lower(:,a)+1)) return
        upper(:,a)=lower(:,a)+d/ratio-1
        do b=1,a-1
          if (all(lower(:,a)<=upper(:,b)).and.all(lower(:,b)<=upper(:,a))) return
        end do
        call coarsen_eb_geometry_3d(patches(a)%geometry,ratio,restricted,good)
        if (.not.good) return
        call extract_eb_geometry_3d(coarse,lower(:,a),upper(:,a),section,good)
        if (.not.good) return
        if (.not.matching_metrics(section,restricted)) return
      end do
      valid=.true.
    end subroutine validate_layout
  end subroutine remap_reactive_eb_multipatch_3d

  subroutine advance_reactive_eb_multipatch_full_3d(species,reactions,transport,coarse,lower,ratio,solver,dt, &
      chemistry_enabled,rtol,atol,transport_enabled,viscosity,conduction,diffusion,barodiffusion, &
      coarse_state,coarse_temperature,patches,boundary_change,wall_change,chemistry_change,minimum_theta,ok, &
      chemistry_integrator,subcycle)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(eb_geometry_3d), intent(in) :: coarse
    integer, intent(in) :: lower(:,:),ratio
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt,rtol,atol
    logical, intent(in) :: chemistry_enabled,transport_enabled,viscosity,conduction,diffusion,barodiffusion
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:)
    type(eb_reactive_patch_3d), intent(inout) :: patches(:)
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemistry_change(:),minimum_theta
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: chemistry_integrator
    logical, optional, intent(in) :: subcycle
    type(eb_reactive_patch_3d), allocatable :: children(:)
    real(dp), allocatable :: candidate(:,:,:,:),temps(:,:,:),boundary(:),wall(:),chemical(:)
    real(dp), allocatable :: exchange(:),transport_sum(:)
    logical, allocatable :: active(:,:,:),refined(:,:,:)
    integer, allocatable :: dims(:,:),upper(:,:),order(:)
    integer :: n(3),np,nvar,p,a,b,idx
    real(dp) :: theta,stage_theta
    logical :: valid
    ok=.false.;np=size(patches);nvar=reactive_nvar(size(species))
    if (size(chemistry_change)/=nvar) return
    if (.not.chemistry_enabled.and..not.transport_enabled) then
      call advance_reactive_eb_multipatch_euler_3d(species,coarse,lower,ratio,solver,dt, &
        coarse_state,coarse_temperature,patches,boundary_change,wall_change,ok,subcycle)
      if (ok) then
        chemistry_change=0._dp;minimum_theta=1._dp
      end if
      return
    end if
    if (ratio<1.or..not.ieee_is_finite(dt)) return
    if (dt<=0._dp.or..not.coarse%is_valid()) return
    if (chemistry_enabled) then
      if (.not.all(ieee_is_finite([rtol,atol]))) return
      if (min(rtol,atol)<=0._dp.or.size(reactions)<1) return
    end if
    n=[coarse%nx,coarse%ny,coarse%nz]
    if (any(shape(lower)/=[3,np])) return
    if (any(shape(coarse_state)/=[nvar,n]).or.any(shape(coarse_temperature)/=n)) return
    if (size(boundary_change)/=nvar.or.size(wall_change)/=nvar) return
    allocate(dims(3,np),upper(3,np),order(np),active(n(1),n(2),n(3)),refined(n(1),n(2),n(3)))
    refined=.false.
    do p=1,np
      if (.not.patches(p)%geometry%is_valid()) return
      if (.not.allocated(patches(p)%state).or..not.allocated(patches(p)%temperature)) return
      dims(:,p)=[patches(p)%geometry%nx,patches(p)%geometry%ny,patches(p)%geometry%nz]
      if (any(mod(dims(:,p),ratio)/=0).or.any(lower(:,p)<1).or.any(lower(:,p)>n)) return
      if (any(dims(:,p)/ratio>n-lower(:,p)+1)) return
      upper(:,p)=lower(:,p)+dims(:,p)/ratio-1
      if (any(shape(patches(p)%state)/=[nvar,dims(:,p)]).or. &
        any(shape(patches(p)%temperature)/=dims(:,p))) return
      if (any(refined(lower(1,p):upper(1,p),lower(2,p):upper(2,p),lower(3,p):upper(3,p)))) return
      refined(lower(1,p):upper(1,p),lower(2,p):upper(2,p),lower(3,p):upper(3,p))=.true.
      order(p)=p
    end do
    do a=2,np
      idx=order(a);b=a-1
      do while(b>=1)
        if (.not.box_precedes(lower(:,idx),lower(:,order(b)))) exit
        order(b+1)=order(b);b=b-1
      end do
      order(b+1)=idx
    end do
    active=coarse%volume_fraction>0._dp.and..not.refined
    candidate=coarse_state;temps=coarse_temperature;children=patches
    allocate(boundary(nvar),wall(nvar),chemical(nvar),exchange(nvar),transport_sum(nvar))
    chemical=0._dp;transport_sum=0._dp;theta=1._dp
    if (chemistry_enabled) then
      call reaction_half(valid)
      if (.not.valid) return
    end if
    if (transport_enabled) then
      call transport_half(valid)
      if (.not.valid) return
    end if
    call advance_reactive_eb_multipatch_euler_3d(species,coarse,lower,ratio,solver,dt, &
      candidate,temps,children,boundary,wall,valid,subcycle)
    if (.not.valid) return
    if (transport_enabled) then
      call transport_half(valid)
      if (.not.valid) return
    end if
    if (chemistry_enabled) then
      call reaction_half(valid)
      if (.not.valid) return
    end if
    if (transport_enabled) boundary=boundary+transport_sum
    if (.not.all(ieee_is_finite(chemical)).or..not.all(ieee_is_finite(boundary))) return
    coarse_state=candidate;coarse_temperature=temps;patches=children
    boundary_change=boundary;wall_change=wall;chemistry_change=chemical;minimum_theta=theta;ok=.true.
  contains
    subroutine transport_half(success)
      logical, intent(out) :: success
      call advance_reactive_eb_multipatch_transport_3d(species,transport,coarse,lower,ratio,.5_dp*dt, &
        viscosity,conduction,diffusion,barodiffusion,candidate,temps,children,exchange,stage_theta,success,subcycle)
      if (.not.success) return
      transport_sum=transport_sum+exchange;theta=min(theta,stage_theta)
    end subroutine transport_half

    subroutine reaction_half(success)
      logical, intent(out) :: success
      real(dp) :: before(nvar),after(nvar)
      integer :: position,child_index
      call inventory(before)
      call advance_reactive_chemistry_3d(species,reactions,candidate,temps,n(1),n(2),n(3), &
        .5_dp*dt,rtol,atol,success,active,chemistry_integrator)
      if (.not.success) return
      do position=1,np
        child_index=order(position)
        call advance_reactive_chemistry_3d(species,reactions,children(child_index)%state, &
          children(child_index)%temperature,dims(1,child_index),dims(2,child_index),dims(3,child_index), &
          .5_dp*dt,rtol,atol,success,children(child_index)%geometry%volume_fraction>0._dp,chemistry_integrator)
        if (.not.success) return
      end do
      call inventory(after)
      chemical=chemical+after-before
      do position=1,np
        child_index=order(position)
        call restrict_reactive_eb_state_3d(species,children(child_index)%geometry,ratio, &
          children(child_index)%state,children(child_index)%temperature, &
          candidate(:,lower(1,child_index):upper(1,child_index),lower(2,child_index):upper(2,child_index), &
            lower(3,child_index):upper(3,child_index)), &
          temps(lower(1,child_index):upper(1,child_index),lower(2,child_index):upper(2,child_index), &
            lower(3,child_index):upper(3,child_index)),success)
        if (.not.success) return
      end do
    end subroutine reaction_half

    subroutine inventory(values)
      real(dp), intent(out) :: values(:)
      integer :: component,position,child_index
      do component=1,nvar
        values(component)=coarse%dx*coarse%dy*coarse%dz* &
          sum(candidate(component,:,:,:)*coarse%volume_fraction,mask=active)
        do position=1,np
          child_index=order(position)
          values(component)=values(component)+children(child_index)%geometry%dx*children(child_index)%geometry%dy* &
            children(child_index)%geometry%dz*sum(children(child_index)%state(component,:,:,:)* &
            children(child_index)%geometry%volume_fraction,mask=children(child_index)%geometry%volume_fraction>0._dp)
        end do
      end do
    end subroutine inventory
  end subroutine advance_reactive_eb_multipatch_full_3d

  subroutine advance_reactive_eb_multipatch_euler_3d(species,coarse,lower,ratio,solver,dt, &
      coarse_state,coarse_temperature,patches,boundary_change,wall_change,ok,subcycle)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: coarse
    integer, intent(in) :: lower(:,:),ratio
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:),boundary_change(:),wall_change(:)
    type(eb_reactive_patch_3d), intent(inout) :: patches(:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: subcycle
    call advance_multipatch_euler(species,coarse,lower,ratio,solver,dt, &
      coarse_state,coarse_temperature,patches,boundary_change,wall_change,ok,subcycle=subcycle)
  end subroutine advance_reactive_eb_multipatch_euler_3d

  subroutine advance_reactive_eb_multipatch_transport_3d(species,transport,coarse,lower,ratio,dt, &
      viscosity,conduction,diffusion,barodiffusion,coarse_state,coarse_temperature,patches,boundary_change, &
      minimum_theta,ok,subcycle)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(eb_geometry_3d), intent(in) :: coarse
    integer, intent(in) :: lower(:,:),ratio
    real(dp), intent(in) :: dt
    logical, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:),boundary_change(:),minimum_theta
    type(eb_reactive_patch_3d), intent(inout) :: patches(:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: subcycle
    real(dp) :: wall(size(boundary_change))
    call advance_multipatch_euler(species,coarse,lower,ratio,'rusanov',dt, &
      coarse_state,coarse_temperature,patches,boundary_change,wall,ok,transport,viscosity,conduction, &
      diffusion,barodiffusion,minimum_theta,subcycle)
  end subroutine advance_reactive_eb_multipatch_transport_3d

  subroutine advance_multipatch_euler(species,coarse,lower,ratio,solver,dt, &
      coarse_state,coarse_temperature,patches,boundary_change,wall_change,ok, &
      transport,viscosity,conduction,diffusion,barodiffusion,minimum_theta,subcycle)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: coarse
    integer, intent(in) :: lower(:,:),ratio
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:),boundary_change(:),wall_change(:)
    type(eb_reactive_patch_3d), intent(inout) :: patches(:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: subcycle
    type(gas_transport_species), optional, intent(in) :: transport(:)
    logical, optional, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), optional, intent(inout) :: minimum_theta
    real(dp) :: theta,stage_theta,sub_dt,elapsed,fraction
    integer :: substeps,substep,i,j,k
    real(dp), allocatable :: donor(:,:,:,:),donor_t(:,:,:)
    logical :: faces(6)
    type(eb_reactive_patch_3d), allocatable :: candidate_patches(:),old_patches(:)
    type(eb_integrated_flux_3d), allocatable :: fine_flux(:),total_flux(:)
    type(eb_integrated_flux_3d) :: coarse_flux,effective,restricted_flux
    type(eb_geometry_3d) :: section,restricted
    integer, allocatable :: upper(:,:),dims(:,:),order(:)
    logical, allocatable :: refined(:,:,:)
    real(dp), allocatable :: initial(:,:,:,:),initial_t(:,:,:),candidate(:,:,:,:),temps(:,:,:)
    real(dp), allocatable :: ghost(:,:,:,:),ghost_t(:,:,:),source(:,:,:,:),wall(:),boundary(:)
    integer :: n(3),np,nvar,p,q,a,b,c,idx,offset(3)
    logical :: valid
    ok=.false.;np=size(patches);nvar=reactive_nvar(size(species))
    if (ratio<1.or..not.ieee_is_finite(dt)) return
    substeps=1
    if (present(subcycle)) then
      if (subcycle.and.np>0) then
        substeps=ratio
        if (present(transport)) then
          if (ratio>huge(ratio)/ratio) return
          substeps=ratio*ratio
        end if
      end if
    end if
    if (dt<=0._dp.or..not.coarse%is_valid()) return
    if (any(shape(lower)/=[3,np])) return
    n=[coarse%nx,coarse%ny,coarse%nz]
    if (any(shape(coarse_state)/=[nvar,n]).or.any(shape(coarse_temperature)/=n)) return
    if (size(boundary_change)/=nvar.or.size(wall_change)/=nvar) return
    allocate(upper(3,np),dims(3,np),order(np),refined(n(1),n(2),n(3)))
    refined=.false.
    do p=1,np
      if (.not.patches(p)%geometry%is_valid()) return
      if (.not.allocated(patches(p)%state).or..not.allocated(patches(p)%temperature)) return
      dims(:,p)=[patches(p)%geometry%nx,patches(p)%geometry%ny,patches(p)%geometry%nz]
      if (any(mod(dims(:,p),ratio)/=0).or.any(lower(:,p)<1).or.any(lower(:,p)>n)) return
      if (any(dims(:,p)/ratio>n-lower(:,p)+1)) return
      upper(:,p)=lower(:,p)+dims(:,p)/ratio-1
      if (any(shape(patches(p)%state)/=[nvar,dims(:,p)]).or. &
        any(shape(patches(p)%temperature)/=dims(:,p))) return
      if (any(refined(lower(1,p):upper(1,p),lower(2,p):upper(2,p),lower(3,p):upper(3,p)))) return
      refined(lower(1,p):upper(1,p),lower(2,p):upper(2,p),lower(3,p):upper(3,p))=.true.
      call coarsen_eb_geometry_3d(patches(p)%geometry,ratio,restricted,valid)
      if (.not.valid) return
      call extract_eb_geometry_3d(coarse,lower(:,p),upper(:,p),section,valid)
      if (.not.valid) return
      if (.not.matching_metrics(section,restricted)) return
      order(p)=p
    end do
    ! Lexicographic box origins define face ownership and summation order,
    ! independent of the caller's patch-array permutation.
    do a=2,np
      idx=order(a);b=a-1
      do while(b>=1)
        if (.not.box_precedes(lower(:,idx),lower(:,order(b)))) exit
        order(b+1)=order(b);b=b-1
      end do
      order(b+1)=idx
    end do
    initial=coarse_state;initial_t=coarse_temperature
    do a=1,np
      p=order(a)
      call restrict_reactive_eb_state_3d(species,patches(p)%geometry,ratio,patches(p)%state,patches(p)%temperature, &
        initial(:,lower(1,p):upper(1,p),lower(2,p):upper(2,p),lower(3,p):upper(3,p)), &
        initial_t(lower(1,p):upper(1,p),lower(2,p):upper(2,p),lower(3,p):upper(3,p)),valid)
      if (.not.valid) return
    end do
    candidate=initial;temps=initial_t;candidate_patches=patches
    allocate(fine_flux(np),total_flux(np),wall(nvar),boundary(nvar),source(nvar,n(1),n(2),n(3)))
    theta=1._dp;source=0._dp
    if (present(transport)) then
      call advance_reactive_eb_transport_patch_euler_3d(species,transport,coarse,dt,viscosity,conduction, &
        diffusion,barodiffusion,candidate,temps,coarse_flux,theta,valid)
      if (.not.valid) return
    else
      call advance_reactive_eb_euler_3d(species,initial,initial_t,coarse,solver,dt,candidate,temps,valid, &
        integrated_flux=coarse_flux)
      if (.not.valid) return
      call reactive_eb_slip_wall_source_3d(species,initial,initial_t,coarse,source,valid)
      if (.not.valid) return
    end if
    do c=1,nvar
      wall(c)=dt*coarse%dx*coarse%dy*coarse%dz*sum(source(c,:,:,:)*coarse%volume_fraction,mask=.not.refined)
    end do
    deallocate(source)
    donor=initial;donor_t=initial_t;elapsed=0._dp
    do substep=1,substeps
      sub_dt=dt/real(substeps,dp)
      if (substep==substeps) sub_dt=dt-elapsed
      if (sub_dt<=0._dp.or.elapsed+sub_dt<=elapsed) return
      fraction=elapsed/dt
      if (substep>1) then
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              if (coarse%volume_fraction(i,j,k)<=0._dp) cycle
              donor(:,i,j,k)=(1._dp-fraction)*initial(:,i,j,k)+fraction*candidate(:,i,j,k)
              donor_t(i,j,k)=(1._dp-fraction)*initial_t(i,j,k)+fraction*temps(i,j,k)
            end do
          end do
        end do
      end if
      ! All peers see the same accepted fine-time snapshot. Halo construction
      ! recovers EOS temperatures from interpolated conserved coarse donors.
      old_patches=candidate_patches
      do a=1,np
        p=order(a)
        allocate(ghost(nvar,0:dims(1,p)+1,0:dims(2,p)+1,0:dims(3,p)+1), &
          ghost_t(0:dims(1,p)+1,0:dims(2,p)+1,0:dims(3,p)+1))
        call fill_reactive_eb_patch_halo_3d(species,coarse,old_patches(p)%geometry,lower(:,p),ratio,donor,donor_t, &
          old_patches(p)%state,old_patches(p)%temperature,ghost,ghost_t,valid)
        if (.not.valid) return
        do b=1,np
          q=order(b)
          if (q==p) cycle
          offset=lower(:,q)-lower(:,p)
          if (any(offset < -dims(:,q)/ratio).or.any(offset>dims(:,p)/ratio)) cycle
          offset=offset*ratio
          call overlay_reactive_eb_peer_halo_3d(species,old_patches(p)%geometry,old_patches(q)%geometry,offset, &
            old_patches(q)%state,old_patches(q)%temperature,ghost,ghost_t,valid)
          if (.not.valid) return
        end do
        if (present(transport)) then
          faces=[lower(1,p)>1,upper(1,p)<n(1),lower(2,p)>1,upper(2,p)<n(2),lower(3,p)>1,upper(3,p)<n(3)]
          call advance_reactive_eb_transport_patch_euler_3d(species,transport,old_patches(p)%geometry,sub_dt,viscosity, &
            conduction,diffusion,barodiffusion,candidate_patches(p)%state,candidate_patches(p)%temperature, &
            fine_flux(p),stage_theta,valid,ghost,ghost_t,faces)
          if (.not.valid) return
          theta=min(theta,stage_theta)
        else
          call advance_reactive_eb_euler_3d(species,old_patches(p)%state,old_patches(p)%temperature,old_patches(p)%geometry, &
            solver,sub_dt,candidate_patches(p)%state,candidate_patches(p)%temperature,valid,ghost,ghost_t,fine_flux(p))
          if (.not.valid) return
        end if
        deallocate(ghost,ghost_t)
        allocate(source(nvar,dims(1,p),dims(2,p),dims(3,p)))
        source=0._dp
        if (.not.present(transport)) then
          call reactive_eb_slip_wall_source_3d(species,old_patches(p)%state,old_patches(p)%temperature, &
            old_patches(p)%geometry,source,valid)
          if (.not.valid) return
        end if
        do c=1,nvar
          wall(c)=wall(c)+sub_dt*old_patches(p)%geometry%dx*old_patches(p)%geometry%dy*old_patches(p)%geometry%dz* &
            sum(source(c,:,:,:)*old_patches(p)%geometry%volume_fraction)
        end do
        deallocate(source)
      end do
      do a=1,np
        p=order(a)
        do b=a+1,np
          q=order(b);offset=lower(:,q)-lower(:,p)
          if (any(offset < -dims(:,q)/ratio).or.any(offset>dims(:,p)/ratio)) cycle
          offset=offset*ratio
          call reconcile_reactive_eb_peer_flux_3d(species,old_patches(p)%geometry,old_patches(q)%geometry,offset, &
            fine_flux(p),fine_flux(q),candidate_patches(q)%state,candidate_patches(q)%temperature,valid)
          if (.not.valid) return
        end do
      end do
      do a=1,np
        p=order(a)
        if (substep==1) then
          total_flux(p)=fine_flux(p)
        else
          total_flux(p)%x=total_flux(p)%x+fine_flux(p)%x
          total_flux(p)%y=total_flux(p)%y+fine_flux(p)%y
          total_flux(p)%z=total_flux(p)%z+fine_flux(p)%z
        end if
      end do
      elapsed=elapsed+sub_dt
    end do
    effective=coarse_flux
    do a=1,np
      p=order(a)
      call restrict_eb_integrated_flux_3d(patches(p)%geometry,ratio,total_flux(p),restricted_flux,valid)
      if (.not.valid) return
      effective%x(:,lower(1,p)-1:upper(1,p),lower(2,p):upper(2,p),lower(3,p):upper(3,p))=restricted_flux%x
      effective%y(:,lower(1,p):upper(1,p),lower(2,p)-1:upper(2,p),lower(3,p):upper(3,p))=restricted_flux%y
      effective%z(:,lower(1,p):upper(1,p),lower(2,p):upper(2,p),lower(3,p)-1:upper(3,p))=restricted_flux%z
    end do
    if (np>0) then
      call reflux_reactive_eb_state_3d(species,coarse,refined,coarse_flux,effective,candidate,temps,valid)
      if (.not.valid) return
    end if
    do a=1,np
      p=order(a)
      call restrict_reactive_eb_state_3d(species,patches(p)%geometry,ratio, &
        candidate_patches(p)%state,candidate_patches(p)%temperature, &
        candidate(:,lower(1,p):upper(1,p),lower(2,p):upper(2,p),lower(3,p):upper(3,p)), &
        temps(lower(1,p):upper(1,p),lower(2,p):upper(2,p),lower(3,p):upper(3,p)),valid)
      if (.not.valid) return
    end do
    do c=1,nvar
      boundary(c)=coarse%dy*coarse%dz*sum( &
        coarse%x_face_fraction(0,:,:)*effective%x(c,0,:,:)-coarse%x_face_fraction(n(1),:,:)*effective%x(c,n(1),:,:)) &
        +coarse%dx*coarse%dz*sum( &
        coarse%y_face_fraction(:,0,:)*effective%y(c,:,0,:)-coarse%y_face_fraction(:,n(2),:)*effective%y(c,:,n(2),:)) &
        +coarse%dx*coarse%dy*sum( &
        coarse%z_face_fraction(:,:,0)*effective%z(c,:,:,0)-coarse%z_face_fraction(:,:,n(3))*effective%z(c,:,:,n(3)))
    end do
    if (.not.all(ieee_is_finite(boundary)).or..not.all(ieee_is_finite(wall))) return
    coarse_state=candidate;coarse_temperature=temps;patches=candidate_patches
    boundary_change=boundary;wall_change=wall
    if (present(minimum_theta)) minimum_theta=theta
    ok=.true.
  end subroutine advance_multipatch_euler

  pure logical function box_precedes(left,right) result(precedes)
    integer, intent(in) :: left(3),right(3)
    integer :: axis
    precedes=.false.
    do axis=1,3
      if (left(axis)==right(axis)) cycle
      precedes=left(axis)<right(axis);return
    end do
  end function box_precedes

  ! Symmetric reaction/hydro composition; the Euler hydro stage remains first
  ! order. Only exposed coarse cells and active fine cells own chemical solves.
  subroutine advance_reactive_eb_two_level_strang_3d(species,reactions,coarse,fine,lower,ratio,solver,dt, &
      chemistry_enabled,rtol,atol,coarse_state,coarse_temperature,fine_state,fine_temperature, &
      boundary_change,wall_change,chemistry_change,ok,chemistry_integrator,subcycle)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(eb_geometry_3d), intent(in) :: coarse,fine
    integer, intent(in) :: lower(3),ratio
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt,rtol,atol
    logical, intent(in) :: chemistry_enabled
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:)
    real(dp), intent(inout) :: fine_state(:,:,:,:),fine_temperature(:,:,:)
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemistry_change(:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: subcycle
    character(len=*), optional, intent(in) :: chemistry_integrator
    call advance_two_level_split(species,reactions,coarse,fine,lower,ratio,solver,dt,chemistry_enabled,rtol,atol, &
      coarse_state,coarse_temperature,fine_state,fine_temperature,boundary_change,wall_change,chemistry_change, &
      ok,chemistry_integrator,subcycle=subcycle)
  end subroutine advance_reactive_eb_two_level_strang_3d

  subroutine advance_reactive_eb_two_level_full_3d(species,reactions,transport,coarse,fine,lower,ratio,solver,dt, &
      chemistry_enabled,rtol,atol,transport_enabled,viscosity,conduction,diffusion,barodiffusion, &
      coarse_state,coarse_temperature,fine_state,fine_temperature,boundary_change,wall_change,chemistry_change, &
      minimum_theta,ok,chemistry_integrator,subcycle)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(eb_geometry_3d), intent(in) :: coarse,fine
    integer, intent(in) :: lower(3),ratio
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt,rtol,atol
    logical, intent(in) :: chemistry_enabled,transport_enabled,viscosity,conduction,diffusion,barodiffusion
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:),fine_state(:,:,:,:),fine_temperature(:,:,:)
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemistry_change(:),minimum_theta
    logical, intent(out) :: ok
    logical, optional, intent(in) :: subcycle
    character(len=*), optional, intent(in) :: chemistry_integrator
    if (transport_enabled) then
      call advance_two_level_split(species,reactions,coarse,fine,lower,ratio,solver,dt,chemistry_enabled,rtol,atol, &
        coarse_state,coarse_temperature,fine_state,fine_temperature,boundary_change,wall_change,chemistry_change, &
        ok,chemistry_integrator,transport,viscosity,conduction,diffusion,barodiffusion,minimum_theta,subcycle)
    else
      call advance_two_level_split(species,reactions,coarse,fine,lower,ratio,solver,dt,chemistry_enabled,rtol,atol, &
        coarse_state,coarse_temperature,fine_state,fine_temperature,boundary_change,wall_change,chemistry_change, &
        ok,chemistry_integrator,subcycle=subcycle)
      if (ok) minimum_theta=1._dp
    end if
  end subroutine advance_reactive_eb_two_level_full_3d

  subroutine advance_two_level_split(species,reactions,coarse,fine,lower,ratio,solver,dt, &
      chemistry_enabled,rtol,atol,coarse_state,coarse_temperature,fine_state,fine_temperature, &
      boundary_change,wall_change,chemistry_change,ok,chemistry_integrator, &
      transport,viscosity,conduction,diffusion,barodiffusion,minimum_theta,subcycle)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(eb_geometry_3d), intent(in) :: coarse,fine
    integer, intent(in) :: lower(3),ratio
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt,rtol,atol
    logical, intent(in) :: chemistry_enabled
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:)
    real(dp), intent(inout) :: fine_state(:,:,:,:),fine_temperature(:,:,:)
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemistry_change(:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: subcycle
    character(len=*), optional, intent(in) :: chemistry_integrator
    type(gas_transport_species), optional, intent(in) :: transport(:)
    logical, optional, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), optional, intent(inout) :: minimum_theta
    real(dp), allocatable :: transport_boundary(:),transport_sum(:)
    real(dp) :: stage_theta,theta
    real(dp), allocatable :: candidate(:,:,:,:),temps(:,:,:),child(:,:,:,:),child_t(:,:,:)
    real(dp), allocatable :: boundary(:),wall(:),chemical(:),before(:)
    logical, allocatable :: active_coarse(:,:,:),active_fine(:,:,:)
    integer :: nvar,n(3),nf(3),upper(3)
    logical :: valid
    procedure(advance_reactive_eb_two_level_euler_3d), pointer :: hydro_step
    procedure(advance_reactive_eb_two_level_transport_3d), pointer :: transport_step
    hydro_step=>advance_reactive_eb_two_level_euler_3d
    transport_step=>advance_reactive_eb_two_level_transport_3d
    if (present(subcycle)) then
      if (subcycle) then
        hydro_step=>advance_reactive_eb_two_level_subcycled_euler_3d
        transport_step=>advance_reactive_eb_two_level_subcycled_transport_3d
      end if
    end if
    ok=.false.;nvar=reactive_nvar(size(species))
    if (size(chemistry_change)/=nvar) return
    if (.not.chemistry_enabled.and..not.present(transport)) then
      call hydro_step(species,coarse,fine,lower,ratio,solver,dt, &
        coarse_state,coarse_temperature,fine_state,fine_temperature,boundary_change,wall_change,ok)
      if (ok) chemistry_change=0._dp
      return
    end if
    if (.not.ieee_is_finite(dt)) return
    if (dt<=0._dp.or.ratio<1) return
    if (chemistry_enabled) then
      if (.not.all(ieee_is_finite([rtol,atol]))) return
      if (min(rtol,atol)<=0._dp.or.size(reactions)<1) return
    end if
    if (.not.coarse%is_valid().or..not.fine%is_valid()) return
    n=[coarse%nx,coarse%ny,coarse%nz];nf=[fine%nx,fine%ny,fine%nz]
    if (any(mod(nf,ratio)/=0).or.any(lower<1).or.any(lower>n)) return
    if (any(nf/ratio>n-lower+1)) return
    upper=lower+nf/ratio-1
    if (any(shape(coarse_state)/=[nvar,n]).or.any(shape(coarse_temperature)/=n).or. &
      any(shape(fine_state)/=[nvar,nf]).or.any(shape(fine_temperature)/=nf)) return
    if (size(boundary_change)/=nvar.or.size(wall_change)/=nvar) return
    allocate(candidate,source=coarse_state);allocate(temps,source=coarse_temperature)
    allocate(child,source=fine_state);allocate(child_t,source=fine_temperature)
    allocate(active_coarse(n(1),n(2),n(3)),active_fine(nf(1),nf(2),nf(3)))
    active_coarse=coarse%volume_fraction>0._dp;active_fine=fine%volume_fraction>0._dp
    active_coarse(lower(1):upper(1),lower(2):upper(2),lower(3):upper(3))=.false.
    allocate(boundary(nvar),wall(nvar),chemical(nvar),before(nvar));chemical=0._dp
    if (chemistry_enabled) then
      call reaction_half(valid)
      if (.not.valid) return
    end if
    allocate(transport_boundary(nvar),transport_sum(nvar));transport_sum=0._dp;theta=1._dp
    if (present(transport)) then
      call transport_half(valid)
      if (.not.valid) return
    end if
    call hydro_step(species,coarse,fine,lower,ratio,solver,dt, &
      candidate,temps,child,child_t,boundary,wall,valid)
    if (.not.valid) return
    if (present(transport)) then
      call transport_half(valid)
      if (.not.valid) return
    end if
    if (chemistry_enabled) then
      call reaction_half(valid)
      if (.not.valid) return
    end if
    if (present(transport)) boundary=boundary+transport_sum
    if (.not.all(ieee_is_finite(chemical)).or..not.all(ieee_is_finite(boundary))) return
    coarse_state=candidate;coarse_temperature=temps;fine_state=child;fine_temperature=child_t
    boundary_change=boundary;wall_change=wall;chemistry_change=chemical
    if (present(minimum_theta)) minimum_theta=theta
    ok=.true.
  contains
    subroutine transport_half(success)
      logical, intent(out) :: success
      call transport_step(species,transport,coarse,fine,lower,ratio,.5_dp*dt, &
        viscosity,conduction,diffusion,barodiffusion,candidate,temps,child,child_t,transport_boundary,stage_theta,success)
      if (.not.success) return
      transport_sum=transport_sum+transport_boundary;theta=min(theta,stage_theta)
    end subroutine transport_half

    subroutine reaction_half(success)
      logical, intent(out) :: success
      integer :: c
      real(dp) :: after
      success=.false.
      do c=1,nvar
        before(c)=coarse%dx*coarse%dy*coarse%dz*sum(candidate(c,:,:,:)*coarse%volume_fraction,mask=active_coarse)+ &
          fine%dx*fine%dy*fine%dz*sum(child(c,:,:,:)*fine%volume_fraction,mask=active_fine)
      end do
      call advance_reactive_chemistry_3d(species,reactions,candidate,temps,n(1),n(2),n(3), &
        .5_dp*dt,rtol,atol,success,active_coarse,chemistry_integrator)
      if (.not.success) return
      call advance_reactive_chemistry_3d(species,reactions,child,child_t,nf(1),nf(2),nf(3), &
        .5_dp*dt,rtol,atol,success,active_fine,chemistry_integrator)
      if (.not.success) return
      do c=1,nvar
        after=coarse%dx*coarse%dy*coarse%dz*sum(candidate(c,:,:,:)*coarse%volume_fraction,mask=active_coarse)+ &
          fine%dx*fine%dy*fine%dz*sum(child(c,:,:,:)*fine%volume_fraction,mask=active_fine)
        chemical(c)=chemical(c)+after-before(c)
      end do
      call restrict_reactive_eb_state_3d(species,fine,ratio,child,child_t, &
        candidate(:,lower(1):upper(1),lower(2):upper(2),lower(3):upper(3)), &
        temps(lower(1):upper(1),lower(2):upper(2),lower(3):upper(3)),success)
    end subroutine reaction_half
  end subroutine advance_two_level_split

  ! One static fine patch, a common Euler interval and physical outflow faces.
  ! All four state arrays and both per-step inventories form one transaction.
  subroutine advance_reactive_eb_two_level_euler_3d(species,coarse,fine,lower,ratio,solver,dt, &
      coarse_state,coarse_temperature,fine_state,fine_temperature,boundary_change,wall_change,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: coarse,fine
    integer, intent(in) :: lower(3),ratio
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:)
    real(dp), intent(inout) :: fine_state(:,:,:,:),fine_temperature(:,:,:),boundary_change(:),wall_change(:)
    logical, intent(out) :: ok
    call advance_two_level_euler(species,coarse,fine,lower,ratio,solver,dt, &
      coarse_state,coarse_temperature,fine_state,fine_temperature,boundary_change,wall_change,ok)
  end subroutine advance_reactive_eb_two_level_euler_3d

  ! One coarse Euler step and ratio fine Euler steps. Coarse face donors use
  ! linear interpolation between unrefluxed coarse endpoints. Only the final
  ! synchronized state is published; this does not raise temporal order.
  subroutine advance_reactive_eb_two_level_subcycled_euler_3d(species,coarse,fine,lower,ratio,solver,dt, &
      coarse_state,coarse_temperature,fine_state,fine_temperature,boundary_change,wall_change,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: coarse,fine
    integer, intent(in) :: lower(3),ratio
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:)
    real(dp), intent(inout) :: fine_state(:,:,:,:),fine_temperature(:,:,:),boundary_change(:),wall_change(:)
    logical, intent(out) :: ok
    call advance_two_level_euler(species,coarse,fine,lower,ratio,solver,dt, &
      coarse_state,coarse_temperature,fine_state,fine_temperature,boundary_change,wall_change,ok,subcycle=.true.)
  end subroutine advance_reactive_eb_two_level_subcycled_euler_3d

  subroutine advance_reactive_eb_two_level_transport_3d(species,transport,coarse,fine,lower,ratio,dt, &
      viscosity,conduction,diffusion,barodiffusion,coarse_state,coarse_temperature,fine_state,fine_temperature, &
      boundary_change,minimum_theta,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(eb_geometry_3d), intent(in) :: coarse,fine
    integer, intent(in) :: lower(3),ratio
    real(dp), intent(in) :: dt
    logical, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:),fine_state(:,:,:,:),fine_temperature(:,:,:)
    real(dp), intent(inout) :: boundary_change(:),minimum_theta
    logical, intent(out) :: ok
    real(dp) :: wall(size(boundary_change))
    call advance_two_level_euler(species,coarse,fine,lower,ratio,'rusanov',dt, &
      coarse_state,coarse_temperature,fine_state,fine_temperature,boundary_change,wall,ok, &
      transport,viscosity,conduction,diffusion,barodiffusion,minimum_theta)
  end subroutine advance_reactive_eb_two_level_transport_3d

  ! Diffusive refinement uses ratio squared fine steps per coarse interval.
  subroutine advance_reactive_eb_two_level_subcycled_transport_3d(species,transport,coarse,fine,lower,ratio,dt, &
      viscosity,conduction,diffusion,barodiffusion,coarse_state,coarse_temperature,fine_state,fine_temperature, &
      boundary_change,minimum_theta,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(eb_geometry_3d), intent(in) :: coarse,fine
    integer, intent(in) :: lower(3),ratio
    real(dp), intent(in) :: dt
    logical, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:),fine_state(:,:,:,:),fine_temperature(:,:,:)
    real(dp), intent(inout) :: boundary_change(:),minimum_theta
    logical, intent(out) :: ok
    real(dp) :: wall(size(boundary_change))
    call advance_two_level_euler(species,coarse,fine,lower,ratio,'rusanov',dt, &
      coarse_state,coarse_temperature,fine_state,fine_temperature,boundary_change,wall,ok, &
      transport,viscosity,conduction,diffusion,barodiffusion,minimum_theta,subcycle=.true.)
  end subroutine advance_reactive_eb_two_level_subcycled_transport_3d

  subroutine advance_two_level_euler(species,coarse,fine,lower,ratio,solver,dt, &
      coarse_state,coarse_temperature,fine_state,fine_temperature,boundary_change,wall_change,ok, &
      transport,viscosity,conduction,diffusion,barodiffusion,minimum_theta,subcycle)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: coarse,fine
    integer, intent(in) :: lower(3),ratio
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:)
    real(dp), intent(inout) :: fine_state(:,:,:,:),fine_temperature(:,:,:),boundary_change(:),wall_change(:)
    logical, intent(out) :: ok
    type(gas_transport_species), optional, intent(in) :: transport(:)
    logical, optional, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), optional, intent(inout) :: minimum_theta
    logical, optional, intent(in) :: subcycle
    real(dp) :: coarse_theta,fine_theta,step_theta,sub_dt,elapsed,fraction
    integer :: substeps,substep,i,j,k
    real(dp), allocatable :: donor(:,:,:,:),donor_t(:,:,:),next_child(:,:,:,:),next_t(:,:,:)
    real(dp), allocatable :: sub_wall(:,:,:,:),fine_wall_change(:)
    type(eb_integrated_flux_3d) :: sub_flux
    logical :: faces(6)
    type(eb_geometry_3d) :: restricted,section
    type(eb_integrated_flux_3d) :: coarse_flux,fine_flux,patch_flux,effective
    real(dp), allocatable :: initial(:,:,:,:),initial_t(:,:,:),candidate(:,:,:,:),temps(:,:,:)
    real(dp), allocatable :: child(:,:,:,:),child_t(:,:,:),ghost(:,:,:,:),ghost_t(:,:,:)
    real(dp), allocatable :: coarse_wall(:,:,:,:),fine_wall(:,:,:,:),boundary(:),wall(:)
    logical, allocatable :: refined(:,:,:)
    integer :: nvar,n(3),nf(3),upper(3),c
    real(dp) :: volume
    logical :: valid
    ok=.false.
    substeps=1
    if (ratio<1) return
    if (present(subcycle)) then
      if (subcycle) then
        substeps=ratio
        if (present(transport)) then
          if (ratio>huge(ratio)/ratio) return
          substeps=ratio*ratio
        end if
      end if
    end if
    if (.not.ieee_is_finite(dt)) return
    if (dt<=0._dp.or.ratio<1) return
    if (.not.coarse%is_valid().or..not.fine%is_valid()) return
    n=[coarse%nx,coarse%ny,coarse%nz];nf=[fine%nx,fine%ny,fine%nz]
    if (any(mod(nf,ratio)/=0).or.any(lower<1).or.any(lower>n)) return
    if (any(nf/ratio>n-lower+1)) return
    upper=lower+nf/ratio-1;nvar=reactive_nvar(size(species))
    if (any(shape(coarse_state)/=[nvar,n]).or.any(shape(coarse_temperature)/=n).or. &
      any(shape(fine_state)/=[nvar,nf]).or.any(shape(fine_temperature)/=nf)) return
    if (size(boundary_change)/=nvar.or.size(wall_change)/=nvar) return
    call coarsen_eb_geometry_3d(fine,ratio,restricted,valid)
    if (.not.valid) return
    call extract_eb_geometry_3d(coarse,lower,upper,section,valid)
    if (.not.valid) return
    if (.not.matching_metrics(section,restricted)) return
    allocate(initial,source=coarse_state);allocate(initial_t,source=coarse_temperature)
    call restrict_reactive_eb_state_3d(species,fine,ratio,fine_state,fine_temperature, &
      initial(:,lower(1):upper(1),lower(2):upper(2),lower(3):upper(3)), &
      initial_t(lower(1):upper(1),lower(2):upper(2),lower(3):upper(3)),valid)
    if (.not.valid) return
    allocate(ghost(nvar,0:nf(1)+1,0:nf(2)+1,0:nf(3)+1),ghost_t(0:nf(1)+1,0:nf(2)+1,0:nf(3)+1))
    call fill_reactive_eb_patch_halo_3d(species,coarse,fine,lower,ratio,initial,initial_t, &
      fine_state,fine_temperature,ghost,ghost_t,valid)
    if (.not.valid) return
    allocate(candidate,source=initial);allocate(temps,source=initial_t)
    allocate(child,source=fine_state);allocate(child_t,source=fine_temperature)
    if (present(transport)) then
      faces=[lower(1)>1,upper(1)<n(1),lower(2)>1,upper(2)<n(2),lower(3)>1,upper(3)<n(3)]
      call advance_reactive_eb_transport_patch_euler_3d(species,transport,coarse,dt,viscosity,conduction, &
        diffusion,barodiffusion,candidate,temps,coarse_flux,coarse_theta,valid)
      if (.not.valid) return
    else
      call advance_reactive_eb_euler_3d(species,initial,initial_t,coarse,solver,dt,candidate,temps,valid, &
        integrated_flux=coarse_flux)
      if (.not.valid) return
    end if
    if (substeps==1) then
      if (present(transport)) then
        call advance_reactive_eb_transport_patch_euler_3d(species,transport,fine,dt,viscosity,conduction, &
          diffusion,barodiffusion,child,child_t,fine_flux,fine_theta,valid,ghost,ghost_t,faces)
      else
        call advance_reactive_eb_euler_3d(species,fine_state,fine_temperature,fine,solver,dt,child,child_t,valid, &
          ghost,ghost_t,fine_flux)
      end if
      if (.not.valid) return
    else
      allocate(donor,source=initial);allocate(donor_t,source=initial_t)
      allocate(next_child,source=child);allocate(next_t,source=child_t)
      allocate(sub_wall(nvar,nf(1),nf(2),nf(3)),fine_wall_change(nvar))
      fine_wall_change=0._dp;elapsed=0._dp;fine_theta=1._dp
      do substep=1,substeps
        sub_dt=dt/real(substeps,dp)
        if (substep==substeps) sub_dt=dt-elapsed
        if (sub_dt<=0._dp.or.elapsed+sub_dt<=elapsed) return
        fraction=elapsed/dt
        if (substep>1) then
          do k=1,n(3)
            do j=1,n(2)
              do i=1,n(1)
                if (coarse%volume_fraction(i,j,k)<=0._dp) cycle
                donor(:,i,j,k)=(1._dp-fraction)*initial(:,i,j,k)+fraction*candidate(:,i,j,k)
                donor_t(i,j,k)=(1._dp-fraction)*initial_t(i,j,k)+fraction*temps(i,j,k)
              end do
            end do
          end do
          ! Halo construction recovers EOS temperature from the interpolated
          ! conserved donor; donor_t is only its nonlinear solve guess.
          call fill_reactive_eb_patch_halo_3d(species,coarse,fine,lower,ratio,donor,donor_t, &
            child,child_t,ghost,ghost_t,valid)
          if (.not.valid) return
        end if
        if (present(transport)) then
          next_child=child;next_t=child_t
          call advance_reactive_eb_transport_patch_euler_3d(species,transport,fine,sub_dt,viscosity,conduction, &
            diffusion,barodiffusion,next_child,next_t,sub_flux,step_theta,valid,ghost,ghost_t,faces)
          if (.not.valid) return
          fine_theta=min(fine_theta,step_theta)
        else
          call reactive_eb_slip_wall_source_3d(species,child,child_t,fine,sub_wall,valid)
          if (.not.valid) return
          do c=1,nvar
            fine_wall_change(c)=fine_wall_change(c)+sub_dt*fine%dx*fine%dy*fine%dz* &
              sum(sub_wall(c,:,:,:)*fine%volume_fraction)
          end do
          call advance_reactive_eb_euler_3d(species,child,child_t,fine,solver,sub_dt,next_child,next_t,valid, &
            ghost,ghost_t,sub_flux)
          if (.not.valid) return
        end if
        child=next_child;child_t=next_t
        if (substep==1) then
          fine_flux=sub_flux
        else
          fine_flux%x=fine_flux%x+sub_flux%x
          fine_flux%y=fine_flux%y+sub_flux%y
          fine_flux%z=fine_flux%z+sub_flux%z
        end if
        elapsed=elapsed+sub_dt
      end do
    end if
    call restrict_eb_integrated_flux_3d(fine,ratio,fine_flux,patch_flux,valid)
    if (.not.valid) return
    effective=coarse_flux
    effective%x(:,lower(1)-1:upper(1),lower(2):upper(2),lower(3):upper(3))=patch_flux%x
    effective%y(:,lower(1):upper(1),lower(2)-1:upper(2),lower(3):upper(3))=patch_flux%y
    effective%z(:,lower(1):upper(1),lower(2):upper(2),lower(3)-1:upper(3))=patch_flux%z
    allocate(refined(n(1),n(2),n(3)));refined=.false.
    refined(lower(1):upper(1),lower(2):upper(2),lower(3):upper(3))=.true.
    call reflux_reactive_eb_state_3d(species,coarse,refined,coarse_flux,effective,candidate,temps,valid)
    if (.not.valid) return
    call restrict_reactive_eb_state_3d(species,fine,ratio,child,child_t, &
      candidate(:,lower(1):upper(1),lower(2):upper(2),lower(3):upper(3)), &
      temps(lower(1):upper(1),lower(2):upper(2),lower(3):upper(3)),valid)
    if (.not.valid) return
    allocate(coarse_wall(nvar,n(1),n(2),n(3)),fine_wall(nvar,nf(1),nf(2),nf(3)),boundary(nvar),wall(nvar))
    coarse_wall=0._dp;fine_wall=0._dp
    if (.not.present(transport)) then
    call reactive_eb_slip_wall_source_3d(species,initial,initial_t,coarse,coarse_wall,valid)
    if (.not.valid) return
    call reactive_eb_slip_wall_source_3d(species,fine_state,fine_temperature,fine,fine_wall,valid)
    if (.not.valid) return
    end if
    volume=coarse%dx*coarse%dy*coarse%dz
    do c=1,nvar
      boundary(c)=coarse%dy*coarse%dz*sum( &
        coarse%x_face_fraction(0,:,:)*effective%x(c,0,:,:)-coarse%x_face_fraction(n(1),:,:)*effective%x(c,n(1),:,:)) &
        +coarse%dx*coarse%dz*sum( &
        coarse%y_face_fraction(:,0,:)*effective%y(c,:,0,:)-coarse%y_face_fraction(:,n(2),:)*effective%y(c,:,n(2),:)) &
        +coarse%dx*coarse%dy*sum( &
        coarse%z_face_fraction(:,:,0)*effective%z(c,:,:,0)-coarse%z_face_fraction(:,:,n(3))*effective%z(c,:,:,n(3)))
      if (substeps>1) then
        wall(c)=dt*volume*sum(coarse_wall(c,:,:,:)*coarse%volume_fraction,mask=.not.refined)+fine_wall_change(c)
      else
      wall(c)=dt*(volume*sum(coarse_wall(c,:,:,:)*coarse%volume_fraction,mask=.not.refined) &
        +fine%dx*fine%dy*fine%dz*sum(fine_wall(c,:,:,:)*fine%volume_fraction))
      end if
    end do
    if (.not.all(ieee_is_finite(boundary)).or..not.all(ieee_is_finite(wall))) return
    coarse_state=candidate;coarse_temperature=temps;fine_state=child;fine_temperature=child_t
    boundary_change=boundary;wall_change=wall
    if (present(minimum_theta)) minimum_theta=min(coarse_theta,fine_theta)
    ok=.true.
  end subroutine advance_two_level_euler

  logical function matching_metrics(a,b) result(valid)
    type(eb_geometry_3d), intent(in) :: a,b
    real(dp), parameter :: tol=2048._dp*epsilon(1._dp)
    real(dp) :: length_scale,area_scale
    valid=.false.
    if (any([a%nx,a%ny,a%nz]/=[b%nx,b%ny,b%nz])) return
    length_scale=max(1._dp,abs(a%x_lower),abs(a%x_upper),abs(a%y_lower),abs(a%y_upper), &
      abs(a%z_lower),abs(a%z_upper))
    area_scale=max(a%dx*a%dy,a%dx*a%dz,a%dy*a%dz)
    if (abs(a%x_lower-b%x_lower)>tol*length_scale.or.abs(a%x_upper-b%x_upper)>tol*length_scale) return
    if (abs(a%dx-b%dx)>tol*a%dx) return
    if (abs(a%y_lower-b%y_lower)>tol*length_scale.or.abs(a%y_upper-b%y_upper)>tol*length_scale) return
    if (abs(a%dy-b%dy)>tol*a%dy) return
    if (abs(a%z_lower-b%z_lower)>tol*length_scale.or.abs(a%z_upper-b%z_upper)>tol*length_scale) return
    if (abs(a%dz-b%dz)>tol*a%dz) return
    if (any(abs(a%volume_fraction-b%volume_fraction)>tol*1._dp)) return
    if (any(abs(a%cell_centroid_x-b%cell_centroid_x)>tol*1._dp)) return
    if (any(abs(a%cell_centroid_y-b%cell_centroid_y)>tol*1._dp)) return
    if (any(abs(a%cell_centroid_z-b%cell_centroid_z)>tol*1._dp)) return
    if (any(abs(a%x_face_fraction-b%x_face_fraction)>tol*1._dp)) return
    if (any(abs(a%y_face_fraction-b%y_face_fraction)>tol*1._dp)) return
    if (any(abs(a%z_face_fraction-b%z_face_fraction)>tol*1._dp)) return
    if (any(abs(a%x_face_centroid_y-b%x_face_centroid_y)>tol*1._dp)) return
    if (any(abs(a%x_face_centroid_z-b%x_face_centroid_z)>tol*1._dp)) return
    if (any(abs(a%y_face_centroid_x-b%y_face_centroid_x)>tol*1._dp)) return
    if (any(abs(a%y_face_centroid_z-b%y_face_centroid_z)>tol*1._dp)) return
    if (any(abs(a%z_face_centroid_x-b%z_face_centroid_x)>tol*1._dp)) return
    if (any(abs(a%z_face_centroid_y-b%z_face_centroid_y)>tol*1._dp)) return
    if (any(abs(a%boundary_area-b%boundary_area)>tol*area_scale)) return
    if (any(abs(a%boundary_centroid_x-b%boundary_centroid_x)>tol*length_scale)) return
    if (any(abs(a%boundary_centroid_y-b%boundary_centroid_y)>tol*length_scale)) return
    if (any(abs(a%boundary_centroid_z-b%boundary_centroid_z)>tol*length_scale)) return
    if (any(abs(a%boundary_normal_x-b%boundary_normal_x)>tol*1._dp)) return
    if (any(abs(a%boundary_normal_y-b%boundary_normal_y)>tol*1._dp)) return
    if (any(abs(a%boundary_normal_z-b%boundary_normal_z)>tol*1._dp)) return
    if (any(abs(a%boundary_normal_integral_x-b%boundary_normal_integral_x)>tol*area_scale)) return
    if (any(abs(a%boundary_normal_integral_y-b%boundary_normal_integral_y)>tol*area_scale)) return
    if (any(abs(a%boundary_normal_integral_z-b%boundary_normal_integral_z)>tol*area_scale)) return
    if (any(a%cell_type/=b%cell_type)) return
    valid=.true.
  end function matching_metrics
end module eb_reactive_two_level_3d_mod
