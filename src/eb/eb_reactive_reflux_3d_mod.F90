module eb_reactive_reflux_3d_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use constants_mod, only: density_floor
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use eb_reactive_transfer_3d_mod, only: recover_reactive_eb_cell_3d
  use nasa7_thermo_mod, only: nasa7_species,nasa7_eos_bounds
  use mixture_thermo_mod, only: mixture_temperature_bounds,mixture_mass_properties,valid_mixture_composition
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d,eb_cut_cell_3d,eb_regular_cell_3d
  implicit none
  private
  type, public :: eb_integrated_flux_3d
    ! Time-integrated flux per OPEN face area. Normal indices are 0:n;
    ! component and tangential indices start at 1. Apertures are applied once.
    real(dp), allocatable :: x(:,:,:,:),y(:,:,:,:),z(:,:,:,:)
  end type eb_integrated_flux_3d
  public :: assemble_eb_rereflux_3d,add_eb_rereflux_to_fine_3d,apply_reactive_eb_rereflux_pair_3d
  public :: restrict_eb_integrated_flux_3d,reflux_reactive_eb_state_3d
  public :: reactive_eb_exchange_interval_3d
  public :: compute_reactive_eb_held_projection,compute_accepted_eb_held_clock_prefixes_3d
  public :: supply_reactive_eb_interface_donors_3d
  public :: supply_reactive_eb_endpoint_donor_3d,supply_reactive_eb_held_eos_donor_3d
  public :: supply_reactive_eb_gas_eos_donor_3d,compute_reactive_eb_reactor_bias_3d,admit_reactive_eb_reactor_energy_3d
  public :: reconcile_reactive_eb_interface_budget_3d
  public :: compute_reactive_eb_reflux_increment_3d
  public :: valid_eb_integrated_flux_3d
  public :: reconcile_reactive_eb_peer_flux_3d
  ! Admission of a private donor proposal through the caller's actual
  ! delivery expression. The callback must not mutate public state or invoke
  ! a reactor provider. MPI callers make the result uniform on all ranks.
  abstract interface
    subroutine reactive_eb_funding_admission(node,cell,pre,require_post,ok)
      import dp
      integer, intent(in) :: node,cell(3)
      real(dp), intent(in) :: pre(:)
      logical, intent(in) :: require_post
      logical, intent(out) :: ok
    end subroutine reactive_eb_funding_admission
  end interface
  public :: reactive_eb_funding_admission
contains
  ! Difference between accepted local UEINT and its transported mean,
  ! including an initial scalar reset. This is accounting, not a gas-energy
  ! source or a modification of the pinned reactor's unpacked state.
  subroutine compute_reactive_eb_held_projection(initial,reset_initial,transported,local,projection,ok)
    real(dp), intent(in) :: initial,reset_initial,transported,local
    real(dp), intent(out) :: projection
    logical, intent(out) :: ok
    real(dp) :: difference(1),reset(1)
    projection=0;ok=.false.
    if(.not.all(ieee_is_finite([initial,reset_initial,transported,local]))) return
    difference=local
    call add_rereflux_increment(difference,[-transported],ok)
    if(.not.ok) return
    reset=reset_initial
    call add_rereflux_increment(reset,[-initial],ok)
    if(.not.ok) return
    call add_rereflux_increment(difference,reset,ok)
    if(ok) projection=difference(1)
  end subroutine compute_reactive_eb_held_projection

  ! EBFluxRegister's coarse re-reflux, using the caller's actual 27-point
  ! EB connectivity. Geometry/source may include a halo; only the selected
  ! valid box receives increments. Exterior donors need complete metrics.
  subroutine assemble_eb_rereflux_3d(geometry,source,boundary,connected,increment,ok,valid_lower,valid_upper)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: source(:,:,:,:)
    logical, intent(in) :: boundary(:,:,:),connected(:,:,:,:,:,:)
    real(dp), allocatable, intent(inout) :: increment(:,:,:,:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: valid_lower(3),valid_upper(3)
    real(dp), allocatable :: candidate(:,:,:,:)
    real(dp) :: dm(size(source,1)),drho(size(source,1)),total,kappa,factor
    integer :: dims(3),lo(3),hi(3),i,j,k,ii,jj,kk,other(3)
    logical :: valid
    ok=.false.
    if(.not.geometry%is_valid().or.size(source,1)<1) return
    dims=[geometry%nx,geometry%ny,geometry%nz];lo=1;hi=dims
    if(present(valid_lower).neqv.present(valid_upper)) return
    if(present(valid_lower)) then
      lo=valid_lower;hi=valid_upper
    endif
    if(any(lo<1).or.any(hi>dims).or.any(hi<lo)) return
    if(any(shape(source)/=[size(source,1),dims]).or.any(shape(boundary)/=dims)) return
    if(any(shape(connected)/=[3,3,3,dims]).or.any(.not.ieee_is_finite(source))) return
    allocate(candidate,mold=source);candidate=0
    do k=max(1,lo(3)-1),min(dims(3),hi(3)+1)
    do j=max(1,lo(2)-1),min(dims(2),hi(2)+1)
    do i=max(1,lo(1)-1),min(dims(1),hi(1)+1)
      if(.not.boundary(i,j,k).or.all(source(:,i,j,k)==0)) cycle
      select case(geometry%cell_type(i,j,k))
      case(eb_regular_cell_3d)
        if(all([i,j,k]>=lo).and.all([i,j,k]<=hi)) then
          call add_rereflux_increment(candidate(:,i,j,k),source(:,i,j,k),valid)
          if(.not.valid) return
        endif
      case(eb_cut_cell_3d)
        kappa=geometry%volume_fraction(i,j,k);dm=source(:,i,j,k)*kappa;total=0
        do kk=-1,1
        do jj=-1,1
        do ii=-1,1
          if(ii==0.and.jj==0.and.kk==0) cycle
          if(.not.connected(ii+2,jj+2,kk+2,i,j,k)) cycle
          other=[i+ii,j+jj,k+kk]
          if(any(other<1).or.any(other>dims)) return
          if(geometry%volume_fraction(other(1),other(2),other(3))<=0) return
          total=total+geometry%volume_fraction(other(1),other(2),other(3))
        enddo
        enddo
        enddo
        if(total<=0) return
        if(1-kappa>huge(1._dp)*min(1._dp,total)) return
        factor=(1-kappa)/total
        if(factor>1) then
          if(any(abs(dm)>huge(1._dp)/factor)) return
        endif
        drho=dm*factor
        if(all([i,j,k]>=lo).and.all([i,j,k]<=hi)) then
          call add_rereflux_increment(candidate(:,i,j,k),dm,valid)
          if(.not.valid) return
        endif
        do kk=-1,1
        do jj=-1,1
        do ii=-1,1
          if(ii==0.and.jj==0.and.kk==0) cycle
          if(.not.connected(ii+2,jj+2,kk+2,i,j,k)) cycle
          other=[i+ii,j+jj,k+kk]
          if(any(other<lo).or.any(other>hi)) cycle
          call add_rereflux_increment(candidate(:,other(1),other(2),other(3)),drho,valid)
          if(.not.valid) return
        enddo
        enddo
        enddo
      end select
    enddo
    enddo
    enddo
    call move_alloc(candidate,increment);ok=.true.
  end subroutine assemble_eb_rereflux_3d

  ! Coarse corrections in fine-covered cells are injected as density
  ! increments into every fine child selected by the inside mask.
  subroutine add_eb_rereflux_to_fine_3d(coarse,inside,ratio,coarse_origin,fine_origin,fine,ok)
    real(dp), intent(in) :: coarse(:,:,:,:)
    logical, intent(in) :: inside(:,:,:)
    integer, intent(in) :: ratio(3),coarse_origin(3),fine_origin(3)
    real(dp), intent(inout) :: fine(:,:,:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: candidate(:,:,:,:)
    integer(int64) :: global(3),index64(3)
    integer :: n(3),i,j,k,index(3)
    logical :: valid
    ok=.false.;n=shape(inside)
    if(size(fine,1)<1.or.any(ratio<1).or.any(n<1).or.any(shape(fine(1,:,:,:))<1)) return
    if(any(shape(coarse)/=[size(fine,1),n])) return
    if(any(.not.ieee_is_finite(coarse)).or.any(.not.ieee_is_finite(fine))) return
    candidate=fine
    do k=1,size(fine,4)
    do j=1,size(fine,3)
    do i=1,size(fine,2)
      global=int(fine_origin,int64)+int([i,j,k]-1,int64)
      ! Integer coarsening rounds toward minus infinity, including halo indices.
      index64=(global-modulo(global,int(ratio,int64)))/int(ratio,int64)-int(coarse_origin,int64)+1
      if(any(index64<1).or.any(index64>int(n,int64))) return
      index=int(index64)
      if(.not.inside(index(1),index(2),index(3))) cycle
      call add_rereflux_increment(candidate(:,i,j,k),coarse(:,index(1),index(2),index(3)),valid)
      if(.not.valid) return
    enddo
    enddo
    enddo
    fine=candidate;ok=.true.
  end subroutine add_eb_rereflux_to_fine_3d

  ! Publish exposed coarse and fine-covered corrections together after EOS
  ! acceptance. Covered coarse cells are left for the caller's average-down.
  subroutine apply_reactive_eb_rereflux_pair_3d(species,coarse_geometry,fine_geometry,increment,inside,ratio, &
      coarse_origin,fine_origin,coarse_state,coarse_temperature,fine_state,fine_temperature,ok,frame)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: coarse_geometry,fine_geometry
    real(dp), intent(in) :: increment(:,:,:,:)
    logical, intent(in) :: inside(:,:,:)
    integer, intent(in) :: ratio(3),coarse_origin(3),fine_origin(3)
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:),fine_state(:,:,:,:),fine_temperature(:,:,:)
    logical, intent(out) :: ok
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), allocatable :: cs(:,:,:,:),ct(:,:,:),fs(:,:,:,:),ft(:,:,:),fd(:,:,:,:),primitive(:)
    real(dp) :: expected(3),actual(3),spacing(3),speed
    integer :: cn(3),fn(3),nv,i,j,k
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    if(.not.coarse_geometry%is_valid().or..not.fine_geometry%is_valid().or.any(ratio<1)) return
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    endif
    cn=[coarse_geometry%nx,coarse_geometry%ny,coarse_geometry%nz]
    fn=[fine_geometry%nx,fine_geometry%ny,fine_geometry%nz]
    if(any(shape(coarse_state)/=[nv,cn]).or.any(shape(increment)/=[nv,cn]).or.any(shape(inside)/=cn)) return
    if(any(shape(fine_state)/=[nv,fn]).or.any(shape(coarse_temperature)/=cn).or.any(shape(fine_temperature)/=fn)) return
    if(any(.not.ieee_is_finite(increment))) return
    spacing=[coarse_geometry%dx,coarse_geometry%dy,coarse_geometry%dz]
    expected=spacing/real(ratio,dp);actual=[fine_geometry%dx,fine_geometry%dy,fine_geometry%dz]
    if(any(abs(expected-actual)>512*epsilon(1._dp)*spacing)) return
    expected=[coarse_geometry%x_lower,coarse_geometry%y_lower,coarse_geometry%z_lower]+ &
      (real(fine_origin,dp)/real(ratio,dp)-real(coarse_origin,dp))*spacing
    actual=[fine_geometry%x_lower,fine_geometry%y_lower,fine_geometry%z_lower]
    if(any(abs(expected-actual)>512*epsilon(1._dp)*max(1._dp,abs(expected)))) return
    cs=coarse_state;ct=coarse_temperature;fs=fine_state;ft=fine_temperature
    allocate(fd,mold=fine_state);fd=0
    call add_eb_rereflux_to_fine_3d(increment,inside,ratio,coarse_origin,fine_origin,fd,valid)
    if(.not.valid) return
    allocate(primitive(reactive_nprim(size(species))))
    do k=1,cn(3)
    do j=1,cn(2)
    do i=1,cn(1)
      if(inside(i,j,k).or.coarse_geometry%volume_fraction(i,j,k)<=0.or.all(increment(:,i,j,k)==0)) cycle
      call add_rereflux_increment(cs(:,i,j,k),increment(:,i,j,k),valid)
      if(.not.valid.or.any(cs(6:,i,j,k)<0)) return
      call recover_reactive_eb_cell_3d(species,coarse_geometry,[i,j,k],cs(:,i,j,k),coarse_temperature(i,j,k), &
        primitive,ct(i,j,k),speed,valid,frame)
      if(.not.valid) return
    enddo
    enddo
    enddo
    do k=1,fn(3)
    do j=1,fn(2)
    do i=1,fn(1)
      if(fine_geometry%volume_fraction(i,j,k)<=0.or.all(fd(:,i,j,k)==0)) cycle
      call add_rereflux_increment(fs(:,i,j,k),fd(:,i,j,k),valid)
      if(.not.valid.or.any(fs(6:,i,j,k)<0)) return
      call recover_reactive_eb_cell_3d(species,fine_geometry,[i,j,k],fs(:,i,j,k),fine_temperature(i,j,k), &
        primitive,ft(i,j,k),speed,valid,frame)
      if(.not.valid) return
    enddo
    enddo
    enddo
    coarse_state=cs;coarse_temperature=ct;fine_state=fs;fine_temperature=ft;ok=.true.
  end subroutine apply_reactive_eb_rereflux_pair_3d

  ! Rebuild one complete accepted clock forest in its actual entry/exit order.
  ! A parent reset occurs after all descendants, even though its record is
  ! before them. Only accepted scalar impulses belong here; a replacement
  ! forest uses the immutable prefix before its first entry, not trial totals.
  ! All four arrays publish together, after shape/order/finite/overflow gates.
  subroutine compute_accepted_eb_held_clock_prefixes_3d(depth,start_reset,local_projection,local_reset,end_reset, &
      initial_projection,initial_reset,projection_before,reset_before,projection_after,reset_after,ok)
    integer, intent(in) :: depth(:)
    real(dp), intent(in) :: start_reset(:),local_projection(:),local_reset(:),end_reset(:)
    real(dp), intent(in) :: initial_projection,initial_reset
    real(dp), allocatable, intent(inout) :: projection_before(:),reset_before(:),projection_after(:),reset_after(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: pb(:),rb(:),pa(:),ra(:)
    integer, allocatable :: stack(:)
    real(dp) :: current(2)
    integer :: n,r,top,node,ios
    logical :: valid
    ok=.false.;n=size(depth)
    if(n<1) return
    if(size(start_reset)/=n.or.size(local_projection)/=n.or.size(local_reset)/=n.or.size(end_reset)/=n) return
    if(any(depth<depth(1)).or.depth(1)<0) return
    do r=2,n
      if(depth(r)>depth(r-1)) then
        if(depth(r)-depth(r-1)/=1) return
      endif
    enddo
    if(.not.all(ieee_is_finite([initial_projection,initial_reset]))) return
    if(.not.all(ieee_is_finite(start_reset)).or..not.all(ieee_is_finite(local_projection)).or. &
      .not.all(ieee_is_finite(local_reset)).or..not.all(ieee_is_finite(end_reset))) return
    allocate(pb(n),rb(n),pa(n),ra(n),stack(n),stat=ios)
    if(ios/=0) return
    current=[initial_projection,initial_reset];top=0
    do r=1,n
      ! Close the previous sibling and its descendants before the next entry.
      do while(top>0)
        node=stack(top)
        if(depth(node)<depth(r)) exit
        call add_rereflux_increment(current,[0._dp,end_reset(node)],valid)
        if(.not.valid) return
        pa(node)=current(1);ra(node)=current(2);top=top-1
      enddo
      call add_rereflux_increment(current,[0._dp,start_reset(r)],valid)
      if(.not.valid) return
      pb(r)=current(1);rb(r)=current(2)
      call add_rereflux_increment(current,[local_projection(r),local_reset(r)],valid)
      if(.not.valid) return
      top=top+1;stack(top)=r
    enddo
    ! Accepted parent closures are deepest first at the end of the forest.
    do while(top>0)
      node=stack(top)
      call add_rereflux_increment(current,[0._dp,end_reset(node)],valid)
      if(.not.valid) return
      pa(node)=current(1);ra(node)=current(2);top=top-1
    enddo
    call move_alloc(pb,projection_before);call move_alloc(rb,reset_before)
    call move_alloc(pa,projection_after);call move_alloc(ra,reset_after)
    ok=.true.
  end subroutine compute_accepted_eb_held_clock_prefixes_3d

  subroutine add_rereflux_increment(value,increment,ok)
    real(dp), intent(inout) :: value(:)
    real(dp), intent(in) :: increment(:)
    logical, intent(out) :: ok
    integer :: c
    ok=.false.
    if(size(value)/=size(increment)) return
    if(any(.not.ieee_is_finite(value)).or.any(.not.ieee_is_finite(increment))) return
    do c=1,size(value)
      if(increment(c)>0) then
        if(value(c)>huge(1._dp)-increment(c)) return
      else if(increment(c)<0) then
        if(value(c)<-huge(1._dp)-increment(c)) return
      endif
    enddo
    value=value+increment;ok=.true.
  end subroutine add_rereflux_increment

  ! Find a represented admissible interval for baseline + factor*direction.
  ! Baselines may be rejected predictors with negative density/species. This
  ! permits a donor exchange to repair several physical endpoints together.
  ! Density is strictly positive; species may reach exact zero. The caller
  ! still owns EOS recovery and conservative face/state publication.
  subroutine reactive_eb_exchange_interval_3d(baseline,direction,lower,upper,ok, &
      passive_components)
    real(dp), intent(in) :: baseline(:,:),direction(:,:)
    real(dp), intent(inout) :: lower,upper
    logical, intent(out) :: ok
    ! Trailing independent scalars share every conservative packet and face
    ! blend, but do not impose species nonnegativity on its gas allowance.
    integer, optional, intent(in) :: passive_components
    integer :: active_components,passives
    real(dp) :: lo,hi,b,g,crossing,middle
    integer :: nv,q,c,attempt
    ok=.false.;nv=size(baseline,1)
    passives=0
    if(present(passive_components)) passives=passive_components
    if(passives<0.or.passives>nv-6) return
    active_components=nv-passives
    if(nv<6.or.any(shape(direction)/=shape(baseline))) return
    if(any(.not.ieee_is_finite(baseline)).or.any(.not.ieee_is_finite(direction))) return
    lo=0;hi=1
    do q=1,size(baseline,2)
    do c=1,active_components
      if(c>=2.and.c<=5) cycle
      b=baseline(c,q);g=direction(c,q)
      if(g==0) then
        if(b<0.or.(c==1.and.b==0)) return
      else if(g>0) then
        if(b<0.or.(c==1.and.b==0)) then
          if(-b>g) return
          crossing=-b/g
          if(c==1) crossing=nearest(crossing,1._dp)
          lo=max(lo,crossing)
        endif
        if(b>=0.and.g>huge(1._dp)-b) hi=min(hi,(huge(1._dp)-b)/g)
      else
        if(b<0.or.(c==1.and.b==0)) return
        if(b<=-g) then
          crossing=b/(-g)
          if(c==1.and.crossing>0) crossing=nearest(crossing,-1._dp)
          hi=min(hi,crossing)
        endif
      endif
    enddo
    enddo
    if(lo>hi) return
    ! Analytic ratios can round onto a forbidden represented endpoint. Move
    ! inward without clipping a state or relaxing the accepted-state contract.
    do attempt=1,16
      if(represented(lo)) exit
      if(lo==hi) return
      lo=nearest(lo,1._dp)
      if(lo>hi) return
    enddo
    do attempt=1,16
      if(represented(hi)) exit
      if(lo==hi) return
      hi=nearest(hi,-1._dp)
      if(lo>hi) return
    enddo
    if(.not.represented(lo).or..not.represented(hi)) then
      middle=lo+.5_dp*(hi-lo)
      if(.not.represented(middle)) return
      if(.not.represented(lo)) lo=middle
      if(.not.represented(hi)) hi=middle
    endif
    lower=lo;upper=hi;ok=.true.
  contains
    logical function represented(factor) result(good)
      real(dp), intent(in) :: factor
      real(dp) :: amount,value
      integer :: column,component
      good=.false.
      do column=1,size(baseline,2)
      do component=1,nv
        if(component>=2.and.component<=5) cycle
        amount=factor*direction(component,column)
        if(amount>0) then
          if(baseline(component,column)>huge(1._dp)-amount) return
        else if(amount<0) then
          if(baseline(component,column)<-huge(1._dp)-amount) return
        endif
        value=baseline(component,column)+amount
        if(.not.ieee_is_finite(value)) return
        if(component<=active_components) then
          if(value<0.or.(component==1.and.value==0)) return
        endif
      enddo
      enddo
      good=.true.
    end function represented
  end subroutine reactive_eb_exchange_interval_3d

  ! Move one conservative packet from an exposed same-level donor toward a
  ! reflux endpoint. The caller publishes the matching interior face packet
  ! and tests EOS before committing either state. A partial packet can be
  ! followed by another donor; failure leaves all outputs untouched.
  subroutine supply_reactive_eb_endpoint_donor_3d(target,correction,target_volume,donor,donor_volume, &
      funded_target,funded_donor,packet,complete,ok, &
      passive_components,allow_inventory_packet)
    real(dp), intent(in) :: target(:),correction(:),target_volume,donor(:),donor_volume
    real(dp), intent(inout) :: funded_target(:),funded_donor(:),packet(:)
    logical, intent(out) :: complete,ok
    ! Trailing independent scalars share every conservative packet and face
    ! blend, but do not impose species nonnegativity on its gas allowance.
    integer, optional, intent(in) :: passive_components
    logical, optional, intent(in) :: allow_inventory_packet
    integer :: active_components,passives
    real(dp) :: total_volume,scale,amount,factor,required,before,after,pair_scale,margin
    real(dp) :: difference(size(target))
    real(dp) :: transfer(size(target)),candidate_target(size(target)),candidate_donor(size(target))
    integer :: nv,c
    logical :: progress,finished
    ok=.false.;complete=.false.;nv=size(target)
    passives=0
    if(present(passive_components)) passives=passive_components
    if(passives<0.or.passives>nv-6) return
    active_components=nv-passives
    if(nv<6.or.size(correction)/=nv.or.size(donor)/=nv.or.size(funded_target)/=nv.or. &
      size(funded_donor)/=nv.or.size(packet)/=nv) return
    if(any(.not.ieee_is_finite(target)).or.any(.not.ieee_is_finite(correction)).or. &
      any(.not.ieee_is_finite(donor)).or..not.ieee_is_finite(target_volume).or. &
      .not.ieee_is_finite(donor_volume)) return
    if(target_volume<=0.or.donor_volume<=0.or.target_volume>huge(1._dp)-donor_volume) return
    if(target(1)<=0.or.donor(1)<=0.or.any(target(6:active_components)<0).or.any(donor(6:active_components)<0)) return
    total_volume=target_volume+donor_volume
    scale=donor_volume/total_volume
    pair_scale=target_volume*scale
    do c=1,nv
      if(target(c)<0) then
        if(donor(c)>huge(1._dp)+target(c)) return
      else if(target(c)>0) then
        if(donor(c)<-huge(1._dp)+target(c)) return
      endif
      difference(c)=donor(c)-target(c)
      if(abs(difference(c))>1) then
        if(pair_scale>huge(1._dp)/abs(difference(c))) return
      endif
    enddo
    transfer=pair_scale*difference
    if(any(.not.ieee_is_finite(transfer))) return
    factor=0;progress=.false.
    do c=1,active_components
      if(c>=2.and.c<=5) cycle
      if(correction(c)>0) then
        if(target(c)>huge(1._dp)-correction(c)) return
      else if(correction(c)<0) then
        if(target(c)<-huge(1._dp)-correction(c)) return
      endif
      before=target(c)+correction(c)
      if(.not.ieee_is_finite(before)) return
      if(before>=0.and.(c/=1.or.before>0)) cycle
      amount=transfer(c)/target_volume
      if(amount<=0) cycle
      progress=.true.
      if(c==1.and.before==0) then
        margin=32*epsilon(1._dp)*target(1)
        if(margin>=amount) then
          required=1
        else
          required=margin/amount
        endif
      else if(-before>=amount) then
        required=1
      else
        required=-before/amount
      endif
      if(.not.ieee_is_finite(required)) return
      factor=max(factor,required)
    enddo
    if(.not.progress) then
      call supply_inventory_packet()
      return
    endif
    factor=min(1._dp,factor)
    factor=min(1._dp,factor*(1+32*epsilon(1._dp)))
    if(factor<=0) return
    transfer=factor*transfer
    candidate_target=target+transfer/target_volume
    candidate_donor=donor-transfer/donor_volume
    if(any(.not.ieee_is_finite(candidate_target)).or.any(.not.ieee_is_finite(candidate_donor))) return
    if(candidate_target(1)<=0.or.candidate_donor(1)<=0.or. &
      any(candidate_target(6:active_components)<0).or.any(candidate_donor(6:active_components)<0)) return
    finished=.true.;progress=.false.
    do c=1,active_components
      if(c>=2.and.c<=5) cycle
      before=target(c)+correction(c)
      if(correction(c)>0) then
        if(candidate_target(c)>huge(1._dp)-correction(c)) return
      else if(correction(c)<0) then
        if(candidate_target(c)<-huge(1._dp)-correction(c)) return
      endif
      after=candidate_target(c)+correction(c)
      if(.not.ieee_is_finite(after)) return
      if(before<0.or.(c==1.and.before==0)) then
        if(after<before) return
        if(after>before) progress=.true.
      else if(after<0.or.(c==1.and.after==0)) then
        return
      endif
      if(after<0.or.(c==1.and.after==0)) finished=.false.
    enddo
    if(.not.progress) then
      call supply_inventory_packet()
      return
    endif
    funded_target=candidate_target;funded_donor=candidate_donor;packet=transfer
    complete=finished;ok=.true.
  contains
    subroutine supply_inventory_packet()
      real(dp) :: ratio,rate,needed,fraction,loss(nv),gain(nv),before_value,after_value,density_margin
      integer :: component
      logical :: deficit
      ! Registered debt can exhaust a uniform endpoint even though mean
      ! contrast is zero. Spend one complete physical donor packet only if
      ! it closes every gas deficit; partial contrast pooling stays unchanged.
      if(.not.present(allow_inventory_packet)) return
      if(.not.allow_inventory_packet) return
      if(target_volume<1) then
        if(donor_volume>huge(1._dp)*target_volume) return
      endif
      ratio=donor_volume/target_volume
      if(.not.ieee_is_finite(ratio).or.ratio<=0) return
      fraction=0;deficit=.false.
      density_margin=density_floor+32*epsilon(1._dp)*max(target(1),abs(correction(1)))
      do component=1,active_components
        if(component>=2.and.component<=5) cycle
        before_value=target(component)+correction(component)
        if(before_value>=0.and.(component/=1.or.before_value>0)) cycle
        deficit=.true.
        if(donor(component)<=0) return
        if(ratio>1) then
          if(donor(component)>huge(1._dp)/ratio) return
        endif
        rate=donor(component)*ratio
        if(.not.ieee_is_finite(rate).or.rate<=0) return
        needed=-before_value
        if(component==1) then
          if(needed>huge(1._dp)-density_margin) return
          needed=needed+density_margin
        endif
        if(needed>=rate) return
        fraction=max(fraction,needed/rate)
      enddo
      if(.not.deficit.or.fraction<=0) return
      fraction=fraction*(1+32*epsilon(1._dp))
      if(.not.ieee_is_finite(fraction).or.fraction>=1-32*epsilon(1._dp)) return
      loss=fraction*donor
      do component=1,nv
        if(ratio>1) then
          if(abs(loss(component))>huge(1._dp)/ratio) return
        endif
        gain(component)=loss(component)*ratio
        if(abs(target(component))>huge(1._dp)-abs(gain(component))) return
        if(donor_volume>1) then
          if(abs(loss(component))>huge(1._dp)/donor_volume) return
        endif
      enddo
      candidate_target=target+gain;candidate_donor=donor-loss
      if(any(.not.ieee_is_finite(candidate_target)).or.any(.not.ieee_is_finite(candidate_donor))) return
      if(candidate_donor(1)<=0.or.any(candidate_donor(6:active_components)<0)) return
      do component=1,active_components
        if(component>=2.and.component<=5) cycle
        if(correction(component)>0) then
          if(candidate_target(component)>huge(1._dp)-correction(component)) return
        else if(correction(component)<0) then
          if(candidate_target(component)<-huge(1._dp)-correction(component)) return
        endif
        after_value=candidate_target(component)+correction(component)
        if(.not.ieee_is_finite(after_value).or.after_value<0) return
        if(component==1.and.after_value<=density_floor) return
      enddo
      transfer=loss*donor_volume
      if(any(.not.ieee_is_finite(transfer))) return
      funded_target=candidate_target;funded_donor=candidate_donor;packet=transfer
      complete=.true.;ok=.true.
    end subroutine supply_inventory_packet
  end subroutine supply_reactive_eb_endpoint_donor_3d

  ! A complete physical packet can repair actual-UEINT EOS debt even when
  ! density/species are already nonnegative. NASA bounds are signed caloric
  ! constraints, not an H>=0 rule. Four pre/post endpoints share one fraction.
  subroutine supply_reactive_eb_held_eos_donor_3d(species,target_geometry,target_cell,target_guess, &
      donor_geometry,donor_cell,donor_guess,target,correction,target_volume,donor,donor_correction,donor_volume, &
      funded_target,funded_donor,packet,complete,ok,frame,allow_raw_target)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: target_geometry,donor_geometry
    integer, intent(in) :: target_cell(3),donor_cell(3)
    real(dp), intent(in) :: target_guess,donor_guess,target(:),correction(:),target_volume
    real(dp), intent(in) :: donor(:),donor_correction(:),donor_volume
    real(dp), intent(inout) :: funded_target(:),funded_donor(:),packet(:)
    logical, intent(out) :: complete,ok
    type(rotating_frame_3d), optional, intent(in) :: frame
    logical, optional, intent(in) :: allow_raw_target
    real(dp) :: baseline(size(target),4),direction(size(target),4),bounds(size(species),2)
    real(dp) :: trial(size(target),4),transfer(size(target)),rho_bounds(size(target),4)
    real(dp) :: pv(reactive_nprim(size(species))),lo,hi,fraction,ratio
    real(dp) :: tg,sound,b,g,crossing,needed,after_margin
    integer :: nv,active,q,side,attempt
    integer :: sides
    logical :: valid,partial,post_valid,progress,raw_target
    ok=.false.;complete=.false.;active=reactive_nvar(size(species));nv=active+1
    if(size(target)/=nv.or.size(donor)/=nv.or.size(correction)/=nv.or.size(donor_correction)/=nv) return
    if(size(funded_target)/=nv.or.size(funded_donor)/=nv.or.size(packet)/=nv) return
    if(.not.all(ieee_is_finite([target_volume,donor_volume,target_guess,donor_guess]))) return
    if(target_volume<=0.or.donor_volume<=0) return
    if(any(.not.ieee_is_finite(target)).or.any(.not.ieee_is_finite(donor)).or. &
      any(.not.ieee_is_finite(correction)).or.any(.not.ieee_is_finite(donor_correction))) return
    raw_target=.false.
    if(present(allow_raw_target)) raw_target=allow_raw_target.and.all(correction==0)
    if(raw_target) then
      if(target(1)<=density_floor) return
      if(target(1)<1) then
        if(any(abs(target(6:active))>huge(1._dp)*target(1))) return
      endif
      if(.not.valid_mixture_composition(species,target(6:active)/target(1))) return
    endif
    baseline(:,1)=target;baseline(:,2)=target;baseline(:,3)=donor;baseline(:,4)=donor
    call add_rereflux_increment(baseline(:,2),correction,valid)
    if(.not.valid) return
    call add_rereflux_increment(baseline(:,4),-donor_correction,valid)
    if(.not.valid) return
    call recover_reactive_eb_cell_3d(species,target_geometry,target_cell,baseline(:active,2), &
      target_guess,pv,tg,sound,valid,frame,baseline(nv,2))
    ! EOS composition tolerances do not settle signed material debt. A
    ! negative trace still needs a conservative packet before strict closure.
    if(valid.and.baseline(1,2)>density_floor.and.all(baseline(6:active,2)>=0)) return
    do q=1,4
      if(q==2) cycle
      if(baseline(1,q)<=density_floor.or.any(baseline(6:active,q)<0)) return
      if(q==1.and.raw_target) cycle
      if(q==1) then
        call recover_reactive_eb_cell_3d(species,target_geometry,target_cell,baseline(:active,q), &
          target_guess,pv,tg,sound,valid,frame,baseline(nv,q))
      else
        call recover_reactive_eb_cell_3d(species,donor_geometry,donor_cell,baseline(:active,q), &
          donor_guess,pv,tg,sound,valid,frame,baseline(nv,q))
      endif
      if(.not.valid) return
    enddo
    if(target_volume<1) then
      if(donor_volume>huge(1._dp)*target_volume) return
    endif
    ratio=donor_volume/target_volume
    if(.not.ieee_is_finite(ratio).or.ratio<=0) return
    if(ratio>1) then
      if(any(abs(donor)>huge(1._dp)/ratio)) return
    endif
    direction(:,1)=ratio*donor;direction(:,2)=direction(:,1)
    direction(:,3)=-donor;direction(:,4)=-donor
    rho_bounds=baseline
    do q=1,4
      if(rho_bounds(1,q)<-huge(1._dp)+density_floor) return
      rho_bounds(1,q)=rho_bounds(1,q)-density_floor
    enddo
    lo=0;hi=1
    call reactive_eb_exchange_interval_3d(rho_bounds,direction,lo,hi,valid,passive_components=1)
    if(.not.valid) return
    call reactive_eb_eos_search_bounds_3d(species,bounds,sides,valid)
    if(.not.valid) return
    needed=lo
    do q=1,4
    do side=1,sides
      call caloric_margin(baseline(:,q),bounds(:,side),b,valid)
      if(.not.valid) return
      call caloric_margin(direction(:,q),bounds(:,side),g,valid)
      if(.not.valid) return
      if(side==2) then
        b=-b;g=-g
      endif
      if((q==2.or.(q==1.and.raw_target)).and.b<0) then
        if(g<=0) return
        if(-b>g) then
          needed=2
        else
          needed=max(needed,-b/g)
        endif
        cycle
      endif
      if(g==0) then
        if(b<0) return
      else if(g>0.and.b<0) then
        if(-b>g) return
        crossing=-b/g;lo=max(lo,crossing)
      else if(g<0) then
        if(b<0) return
        if(b<=-g) hi=min(hi,b/(-g))
      endif
      if(lo>hi) return
    enddo
    enddo
    ! Check actual EOS inside the analytically admissible interval. An
    ! inward represented point avoids an energy-boundary cancellation; the
    ! midpoint and upper interior also handle a rounded narrow lower bound.
    partial=needed>hi
    if(.not.partial) lo=max(lo,needed)
    do attempt=1,3
      if(partial) then
        fraction=hi-(hi-lo)*32*epsilon(1._dp)
      else
      select case(attempt)
      case(1)
        fraction=lo+(hi-lo)*32*epsilon(1._dp)
      case(2)
        fraction=lo+.5_dp*(hi-lo)
      case(3)
        fraction=hi-(hi-lo)*32*epsilon(1._dp)
      end select
      endif
      if(fraction<=0.or.fraction>=1) cycle
      trial=baseline
      call add_rereflux_increment(trial(:,1),fraction*direction(:,1),valid)
      if(.not.valid) return
      trial(:,2)=trial(:,1)
      call add_rereflux_increment(trial(:,2),correction,valid)
      if(.not.valid) return
      ! Scale the remaining whole material vector together. Subtracting a
      ! nearly full packet componentwise destroys rho=sum(rhoY) at tiny stock.
      trial(:,3)=(1-fraction)*donor;trial(:,4)=trial(:,3)
      call add_rereflux_increment(trial(:,4),-donor_correction,valid)
      if(.not.valid) return
      do q=1,4
        if(partial.and.(q==2.or.(q==1.and.raw_target))) cycle
        if(q<=2) then
          call recover_reactive_eb_cell_3d(species,target_geometry,target_cell,trial(:active,q), &
            target_guess,pv,tg,sound,valid,frame,trial(nv,q))
        else
          call recover_reactive_eb_cell_3d(species,donor_geometry,donor_cell,trial(:active,q), &
            donor_guess,pv,tg,sound,valid,frame,trial(nv,q))
        endif
        if(.not.valid) exit
      enddo
      if(.not.valid) then
        if(partial) return
        cycle
      endif
      post_valid=.true.
      if(partial) then
        ! As with complementary gas donors, a healthy pre-register target
        ! may pool several complete physical packets before debt closes.
        ! Do not worsen any satisfied caloric bound or a remaining deficit.
        progress=.false.
        do side=1,sides
          call caloric_margin(baseline(:,2),bounds(:,side),b,valid)
          if(.not.valid) return
          call caloric_margin(trial(:,2),bounds(:,side),after_margin,valid)
          if(.not.valid) return
          if(side==2) then
            b=-b;after_margin=-after_margin
          endif
          if(b<0) then
            if(after_margin<=b) return
            progress=.true.
          else if(after_margin<0) then
            return
          endif
        enddo
        if(.not.progress) return
        call recover_reactive_eb_cell_3d(species,target_geometry,target_cell,trial(:active,2), &
          target_guess,pv,tg,sound,post_valid,frame,trial(nv,2))
      endif
      transfer=fraction*donor
      if(donor_volume>1) then
        if(any(abs(transfer)>huge(1._dp)/donor_volume)) return
      endif
      transfer=transfer*donor_volume
      if(any(.not.ieee_is_finite(transfer))) return
      funded_target=trial(:,1);funded_donor=trial(:,3);packet=transfer
      complete=post_valid;ok=.true.;return
    enddo
  contains
    subroutine caloric_margin(value,energies,result,good)
      real(dp), intent(in) :: value(nv),energies(size(species))
      real(dp), intent(out) :: result
      logical, intent(out) :: good
      real(dp) :: term
      integer :: component
      good=.false.;result=value(nv)
      do component=1,size(species)
        if(abs(energies(component))>1) then
          if(abs(value(5+component))>huge(1._dp)/abs(energies(component))) return
        endif
        term=-value(5+component)*energies(component)
        if(term>0) then
          if(result>huge(1._dp)-term) return
        else if(term<0) then
          if(result<-huge(1._dp)-term) return
        endif
        result=result+term
      enddo
      good=ieee_is_finite(result)
    end subroutine caloric_margin
  end subroutine supply_reactive_eb_held_eos_donor_3d

  ! Gas-only EOS uses total energy minus kinetic/rotating energy. Along one
  ! complete material packet each NASA bound is quadratic after multiplying
  ! by positive density. Keep disconnected admissible intervals: mixing two
  ! healthy streams can exceed the temperature ceiling in between them.
  ! Strict EOS keeps both fitted caloric bounds. Pinned extrapolation has
  ! no fitted-temperature ceiling; its low branch is e(Tlo)+Cv(Tlo)*(T-Tlo).
  ! Therefore positive T requires e>e(Tlo)-Cv(Tlo)*Tlo, an affine species
  ! bound even for mixed material packets. These are necessary search bounds;
  ! every proposal still passes the actual gas/UEINT/biased-reactor EOS,
  ! including positive Cp/Cv and finite properties, before publication.
  subroutine reactive_eb_eos_search_bounds_3d(species,bounds,sides,ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(out) :: bounds(size(species),2)
    integer, intent(out) :: sides
    logical, intent(out) :: ok
    real(dp) :: y(size(species)),lo,hi,mw,rg,cp,cv,gamma,h,e,entropy,edge(1),delta(1)
    integer :: sp
    logical :: active,valid
    ok=.false.;bounds=0;sides=2
    if(size(species)<1) return
    active=species(1)%pinned_eos_extrapolation
    if(any(species%pinned_eos_extrapolation.neqv.active)) return
    if(active) then
      call nasa7_eos_bounds(species,lo,hi,valid)
      sides=1
    else
      call mixture_temperature_bounds(species,lo,hi,valid)
    endif
    if(.not.valid) return
    do sp=1,size(species)
      y=0;y(sp)=1
      call mixture_mass_properties(species,y,lo,mw,rg,cp,cv,gamma,h,e,entropy,valid)
      if(.not.valid) return
      if(active) then
        if(lo>1) then
          if(abs(cv)>huge(1._dp)/lo) return
        endif
        edge=e;delta=-cv*lo
        call add_rereflux_increment(edge,delta,valid)
        if(.not.valid) return
        bounds(sp,1)=edge(1)
      else
        bounds(sp,1)=e
        call mixture_mass_properties(species,y,hi,mw,rg,cp,cv,gamma,h,e,entropy,valid)
        if(.not.valid) return
        bounds(sp,2)=e
      endif
    enddo
    ok=.true.
  end subroutine reactive_eb_eos_search_bounds_3d

  ! React.cpp packs old UEINT separately from the gas projection used for
  ! its external energy. This read-only offset converts a post-gas EOS query
  ! into the actual forced reactor-energy query, including rotating potential.
  subroutine compute_reactive_eb_reactor_bias_3d(geometry,cell,gas,held,bias,ok,frame)
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: cell(3)
    real(dp), intent(in) :: gas(:),held
    real(dp), intent(inout) :: bias
    logical, intent(out) :: ok
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp) :: scale,velocity,kinetic,shape_factor,projected(1),delta(1),phi,position(3),candidate(1)
    logical :: valid
    ok=.false.
    if(size(gas)<6.or..not.ieee_is_finite(held)) return
    if(any(.not.ieee_is_finite(gas)).or.gas(1)<=density_floor) return
    scale=maxval(abs(gas(2:4)));kinetic=0
    if(scale>0) then
      if(gas(1)<1) then
        if(scale>huge(1._dp)*gas(1)) return
      endif
      velocity=scale/gas(1)
      if(velocity>1) then
        if(scale>huge(1._dp)/velocity) return
      endif
      kinetic=scale*velocity;shape_factor=.5_dp*sum((gas(2:4)/scale)**2)
      if(shape_factor>1) then
        if(kinetic>huge(1._dp)/shape_factor) return
      endif
      kinetic=kinetic*shape_factor
    endif
    projected=gas(5);delta=-kinetic
    call add_rereflux_increment(projected,delta,valid)
    if(.not.valid) return
    if(present(frame)) then
      position=[geometry%x_lower,geometry%y_lower,geometry%z_lower]+ &
        (real(cell,dp)-.5_dp)*[geometry%dx,geometry%dy,geometry%dz]
      call rotating_potential_3d(frame,position,phi,valid)
      if(.not.valid) return
      if(phi>1) then
        if(gas(1)>huge(1._dp)/phi) return
      endif
      delta=gas(1)*phi
      call add_rereflux_increment(projected,delta,valid)
      if(.not.valid) return
    endif
    candidate=held
    call add_rereflux_increment(candidate,-projected,valid)
    if(.not.valid) return
    bias=candidate(1);ok=.true.
  end subroutine compute_reactive_eb_reactor_bias_3d

  subroutine admit_reactive_eb_reactor_energy_3d(species,geometry,cell,gas,bias,guess,ok,frame)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: cell(3)
    real(dp), intent(in) :: gas(:),bias,guess
    logical, intent(out) :: ok
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp) :: shifted(size(gas)),delta(size(gas)),qp(reactive_nprim(size(species))),tg,sound
    ok=.false.
    if(size(gas)/=reactive_nvar(size(species))) return
    call recover_reactive_eb_cell_3d(species,geometry,cell,gas,guess,qp,tg,sound,ok,frame)
    if(.not.ok) return
    shifted=gas;delta=0;delta(5)=bias
    call add_rereflux_increment(shifted,delta,ok)
    if(.not.ok) return
    call recover_reactive_eb_cell_3d(species,geometry,cell,shifted,guess,qp,tg,sound,ok,frame)
  end subroutine admit_reactive_eb_reactor_energy_3d

  subroutine supply_reactive_eb_gas_eos_donor_3d(species,target_geometry,target_cell,target_guess, &
      donor_geometry,donor_cell,donor_guess,target,correction,target_volume,donor,donor_correction,donor_volume, &
      funded_target,funded_donor,packet,complete,ok,frame,allow_raw_target,reactor_energy_bias, &
      native_projection,inventory_nodes,check_native_target)
    logical, optional, intent(in) :: check_native_target
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: target_geometry,donor_geometry
    integer, intent(in) :: target_cell(3),donor_cell(3)
    real(dp), intent(in) :: target_guess,donor_guess,target(:),correction(:),target_volume
    real(dp), intent(in) :: donor(:),donor_correction(:),donor_volume
    real(dp), intent(inout) :: funded_target(:),funded_donor(:),packet(:)
    logical, intent(out) :: complete,ok
    type(rotating_frame_3d), optional, intent(in) :: frame
    logical, optional, intent(in) :: allow_raw_target
    ! With a reactor bias, the final component is independent UEINT. The
    ! read-only density-energy offsets are old UEINT minus old projected gas
    ! internal energy at the target and donor. They constrain admission only;
    ! every transported component remains the actual complete material packet.
    real(dp), optional, intent(in) :: reactor_energy_bias(2)
    procedure(reactive_eb_funding_admission), optional :: native_projection
    integer, optional, intent(in) :: inventory_nodes(2)
    integer, parameter :: capacity=64
    real(dp) :: baseline(size(target),4),direction(size(target),4),trial(size(target),4)
    real(dp) :: normalized(size(target),4),slope(size(target),4),rho_bounds(size(target),4)
    real(dp) :: bounds(size(species),2),scaled_bounds(size(species)),transfer(size(target))
    real(dp) :: pv(reactive_nprim(size(species))),phi(4),position(3)
    real(dp) :: tg,sound,lo,hi,ratio,rho_scale,energy_scale,momentum_scale
    real(dp) :: physical_lo(capacity),physical_hi(capacity),full_lo(capacity),full_hi(capacity)
    real(dp) :: coefficients(3,2,4,3),bias(4),normalized_bias(4),full_lower,full_upper,b0,b1,r0,r1,m0(3),m1(3),fraction,before,after
    integer :: nv,gas_nv,passives,families,family,sp,q,side,physical_count,full_count,interval,attempt,retained
    real(dp) :: refused_fraction,accepted_fraction,probe_fraction
    integer :: refinement
    integer :: sides
    logical :: valid,post_valid,progress,raw_target,joint,native_target_required
    ok=.false.;complete=.false.;gas_nv=reactive_nvar(size(species))
    native_target_required=.false.
    if(present(check_native_target)) native_target_required=check_native_target
    if(native_target_required.and..not.present(native_projection)) return
    if(present(native_projection).and..not.present(inventory_nodes)) return
    if(present(inventory_nodes)) then
      if(any(inventory_nodes<1)) return
    endif
    joint=present(reactor_energy_bias);passives=0;families=1;bias=0
    if(joint) then
      passives=1;families=3
      if(any(.not.ieee_is_finite(reactor_energy_bias))) return
      bias=[reactor_energy_bias(1),reactor_energy_bias(1),reactor_energy_bias(2),reactor_energy_bias(2)]
    endif
    nv=gas_nv+passives
    if(gas_nv<6) return
    if(size(target)/=nv.or.size(donor)/=nv.or.size(correction)/=nv.or.size(donor_correction)/=nv) return
    if(size(funded_target)/=nv.or.size(funded_donor)/=nv.or.size(packet)/=nv) return
    if(.not.all(ieee_is_finite([target_volume,donor_volume,target_guess,donor_guess]))) return
    if(target_volume<=0.or.donor_volume<=0) return
    if(any(.not.ieee_is_finite(target)).or.any(.not.ieee_is_finite(donor)).or. &
      any(.not.ieee_is_finite(correction)).or.any(.not.ieee_is_finite(donor_correction))) return
    raw_target=.false.
    if(present(allow_raw_target)) raw_target=allow_raw_target.and.all(correction==0)
    if(raw_target) then
      if(target(1)<=density_floor) return
      if(target(1)<1) then
        if(any(abs(target(6:gas_nv))>huge(1._dp)*target(1))) return
      endif
      if(.not.valid_mixture_composition(species,target(6:gas_nv)/target(1))) return
    endif
    baseline(:,1)=target;baseline(:,2)=target;baseline(:,3)=donor;baseline(:,4)=donor
    call add_rereflux_increment(baseline(:,2),correction,valid)
    if(.not.valid) return
    call add_rereflux_increment(baseline(:,4),-donor_correction,valid)
    if(.not.valid) return
    call recover_endpoint(baseline(:,2),2,valid)
    ! The downstream gas/species gate is strict even when caloric recovery
    ! tolerates a tiny negative mass fraction. Continue the donor search.
    if(valid.and.baseline(1,2)>density_floor.and.all(baseline(6:gas_nv,2)>=0)) then
      if(.not.native_target_required) return
      ! A statically valid endpoint can still fail its actual forced reactor.
      ! Query the pre-debt target: the caller adds the original source once.
      call native_projection(inventory_nodes(1),target_cell,baseline(:,1),.true.,valid)
      if(valid) return
    endif
    do q=1,4
      if(q==2) cycle
      if(baseline(1,q)<=density_floor.or.any(baseline(6:gas_nv,q)<0)) return
      if(q==1.and.raw_target) cycle
      call recover_endpoint(baseline(:,q),q,valid)
      if(.not.valid) return
    enddo
    if(target_volume<1) then
      if(donor_volume>huge(1._dp)*target_volume) return
    endif
    ratio=donor_volume/target_volume
    if(.not.ieee_is_finite(ratio).or.ratio<=0) return
    if(ratio>1) then
      if(any(abs(donor)>huge(1._dp)/ratio)) return
    endif
    direction(:,1)=ratio*donor;direction(:,2)=direction(:,1)
    direction(:,3)=-donor;direction(:,4)=-donor
    rho_bounds=baseline
    do q=1,4
      if(rho_bounds(1,q)<-huge(1._dp)+density_floor) return
      rho_bounds(1,q)=rho_bounds(1,q)-density_floor
    enddo
    lo=0;hi=1
    if(joint) then
      ! A single healthy donor may supply only part of a material deficit.
      ! Keep the pending receiver out of the physical interval; final closure
      ! still requires its represented density, species and all three EOSs.
      call reactive_eb_exchange_interval_3d(rho_bounds(:,[1,3,4]),direction(:,[1,3,4]),lo,hi,valid,1)
    else
      call reactive_eb_exchange_interval_3d(rho_bounds,direction,lo,hi,valid)
    endif
    if(.not.valid) return
    call reactive_eb_eos_search_bounds_3d(species,bounds,sides,valid)
    if(.not.valid) return
    phi=0
    if(present(frame)) then
      position=[target_geometry%x_lower,target_geometry%y_lower,target_geometry%z_lower]+ &
        (real(target_cell,dp)-.5_dp)*[target_geometry%dx,target_geometry%dy,target_geometry%dz]
      call rotating_potential_3d(frame,position,phi(1),valid)
      if(.not.valid) return
      position=[donor_geometry%x_lower,donor_geometry%y_lower,donor_geometry%z_lower]+ &
        (real(donor_cell,dp)-.5_dp)*[donor_geometry%dx,donor_geometry%dy,donor_geometry%dz]
      call rotating_potential_3d(frame,position,phi(3),valid)
      if(.not.valid) return
      phi(2)=phi(1);phi(4)=phi(3)
    endif
    ! Separate density and energy scales keep the quadratic products bounded
    ! without ever forming rho_scale**2*energy_scale or squared dimensional M.
    rho_scale=max(maxval(abs(baseline(1,:))),maxval(abs(direction(1,:))), &
      maxval(abs(baseline(6:gas_nv,:))),maxval(abs(direction(6:gas_nv,:))))
    if(rho_scale<=0.or..not.ieee_is_finite(rho_scale)) return
    if(rho_scale<1) then
      if(any(abs(baseline(2:5,:))>huge(1._dp)*rho_scale).or. &
        any(abs(direction(2:5,:))>huge(1._dp)*rho_scale)) return
    endif
    if(joint.and.rho_scale<1) then
      if(any(abs(baseline(nv,:))>huge(1._dp)*rho_scale).or. &
        any(abs(direction(nv,:))>huge(1._dp)*rho_scale).or.any(abs(bias)>huge(1._dp)*rho_scale)) return
    endif
    normalized=baseline/rho_scale;slope=direction/rho_scale
    normalized_bias=bias/rho_scale
    momentum_scale=max(maxval(abs(normalized(2:4,:))),maxval(abs(slope(2:4,:))))
    if(momentum_scale>sqrt(huge(1._dp)/16)) return
    energy_scale=max(1._dp,maxval(abs(bounds)),maxval(abs(phi)), &
      maxval(abs(normalized(5,:))),maxval(abs(slope(5,:))),momentum_scale**2)
    if(joint) energy_scale=max(energy_scale,maxval(abs(normalized(nv,:))), &
      maxval(abs(slope(nv,:))),maxval(abs(normalized_bias)))
    normalized(2:4,:)=normalized(2:4,:)/sqrt(energy_scale)
    slope(2:4,:)=slope(2:4,:)/sqrt(energy_scale)
    normalized(5,:)=normalized(5,:)/energy_scale;slope(5,:)=slope(5,:)/energy_scale
    phi=phi/energy_scale;normalized_bias=normalized_bias/energy_scale
    if(joint) then
      normalized(nv,:)=normalized(nv,:)/energy_scale;slope(nv,:)=slope(nv,:)/energy_scale
    endif
    do q=1,4
      r0=normalized(1,q);r1=slope(1,q);m0=normalized(2:4,q);m1=slope(2:4,q)
      do family=1,families
      do side=1,sides
        scaled_bounds=bounds(:,side)/energy_scale
        b0=normalized(5,q)+phi(q)*r0-sum(normalized(6:gas_nv,q)*scaled_bounds)
        b1=slope(5,q)+phi(q)*r1-sum(slope(6:gas_nv,q)*scaled_bounds)
        if(family==2) b0=b0+normalized_bias(q)
        if(family==3) then
          coefficients(:,side,q,family)=[normalized(nv,q)-sum(normalized(6:gas_nv,q)*scaled_bounds), &
            slope(nv,q)-sum(slope(6:gas_nv,q)*scaled_bounds),0._dp]
        else
          coefficients(1,side,q,family)=2*r0*b0-sum(m0*m0)
          coefficients(2,side,q,family)=2*(r0*b1+r1*b0)-2*sum(m0*m1)
          coefficients(3,side,q,family)=2*r1*b1-sum(m1*m1)
        endif
        if(side==2) coefficients(:,side,q,family)=-coefficients(:,side,q,family)
      enddo
      enddo
    enddo
    physical_lo=0;physical_hi=0;full_lo=0;full_hi=0
    physical_count=1;physical_lo(1)=lo;physical_hi(1)=hi
    do q=1,4
      if(q==2.or.(q==1.and.raw_target)) cycle
      do family=1,families
      do side=1,sides
        call clip_quadratic(coefficients(:,side,q,family),physical_lo,physical_hi,physical_count,valid)
        if(.not.valid.or.physical_count==0) return
      enddo
      enddo
    enddo
    full_count=physical_count;full_lo=physical_lo;full_hi=physical_hi
    if(joint) then
      full_lower=lo;full_upper=hi
      call reactive_eb_exchange_interval_3d(rho_bounds(:,2:2),direction(:,2:2),full_lower,full_upper,valid,1)
      if(valid) then
        ! The material interval kernel returns its own complete interval.
        ! Intersect each physical piece rather than replacing donor density
        ! bounds: multiplied caloric polynomials alone also admit rho <= 0.
        retained=0
        do interval=1,physical_count
          if(max(physical_lo(interval),full_lower)>min(physical_hi(interval),full_upper)) cycle
          retained=retained+1
          full_lo(retained)=max(physical_lo(interval),full_lower)
          full_hi(retained)=min(physical_hi(interval),full_upper)
        enddo
        full_count=retained
      else
        full_count=0
      endif
    endif
    do family=1,families
    do side=1,sides
      call clip_quadratic(coefficients(:,side,2,family),full_lo,full_hi,full_count,valid)
      if(.not.valid) return
    enddo
    enddo
    do interval=1,full_count
    refused_fraction=full_lo(interval)
    do attempt=1,3
      select case(attempt)
      case(1)
        fraction=full_lo(interval)+(full_hi(interval)-full_lo(interval))*32*epsilon(1._dp)
      case(2)
        fraction=full_lo(interval)+.5_dp*(full_hi(interval)-full_lo(interval))
      case(3)
        fraction=full_hi(interval)-(full_hi(interval)-full_lo(interval))*32*epsilon(1._dp)
      end select
      call propose(fraction,.true.,valid,post_valid)
      if(valid) then
        if(present(native_projection).and.attempt>1) then
          ! The EOS interval edge may fail after the caller converts the
          ! impulse to a rate and back. Retain an admitted interior packet
          ! while reducing excess borrowing. No monotonic EOS assumption is
          ! needed: every retained proposal passes both private queries and
          ! the actual delivery expression again before publication.
          accepted_fraction=fraction
          do refinement=1,48
            probe_fraction=refused_fraction+.5_dp*(accepted_fraction-refused_fraction)
            if(probe_fraction<=refused_fraction.or.probe_fraction>=accepted_fraction) exit
            call propose(probe_fraction,.true.,valid,post_valid)
            if(valid) then
              accepted_fraction=probe_fraction
            else
              refused_fraction=probe_fraction
            endif
          enddo
          fraction=accepted_fraction
          call propose(fraction,.true.,valid,post_valid)
          if(.not.valid) return
        endif
        call publish(fraction,post_valid)
        return
      endif
      refused_fraction=fraction
    enddo
    enddo
    if(full_count>0) return
    ! A partial physical loan can pool additional providers. Its pending
    ! material or thermal deficits improve; joint reactor funding may still
    ! leave negative pending density. Final closure owns complete post-debt
    ! EOS admission and whole-transaction rollback.
    if(.not.joint) then
      if(baseline(1,2)<=density_floor.or.any(baseline(6:gas_nv,2)<0)) return
    endif
    do interval=physical_count,1,-1
      fraction=physical_hi(interval)-(physical_hi(interval)-physical_lo(interval))*32*epsilon(1._dp)
      call propose(fraction,.false.,valid,post_valid)
      if(.not.valid) cycle
      progress=.false.
      if(joint.and.(baseline(1,2)<=density_floor.or.any(baseline(6:gas_nv,2)<0))) then
        do sp=1,gas_nv
          if(sp>=2.and.sp<=5) cycle
          if(baseline(sp,2)<0.or.(sp==1.and.baseline(sp,2)<=density_floor)) then
            if(trial(sp,2)<=baseline(sp,2)) then
              valid=.false.;exit
            endif
            progress=.true.
          else if(trial(sp,2)<0) then
            valid=.false.;exit
          endif
        enddo
      else
        do family=1,families
        do side=1,sides
          call margin(baseline(:,2),side,family,before,valid)
          if(.not.valid) exit
          call margin(trial(:,2),side,family,after,valid)
          if(.not.valid) exit
          if(before<0) then
            if(after<=before) then
              valid=.false.;exit
            endif
            progress=.true.
          else if(after<0) then
            valid=.false.;exit
          endif
        enddo
        if(.not.valid) exit
        enddo
      endif
      if(valid.and.progress) then
        call publish(fraction,post_valid)
        return
      endif
    enddo
  contains
    subroutine recover_endpoint(value,column,good)
      real(dp), intent(in) :: value(nv)
      integer, intent(in) :: column
      logical, intent(out) :: good
      real(dp) :: shifted(gas_nv),shift(gas_nv)
      if(column<=2) then
        call recover_reactive_eb_cell_3d(species,target_geometry,target_cell,value(:gas_nv), &
          target_guess,pv,tg,sound,good,frame)
      else
        call recover_reactive_eb_cell_3d(species,donor_geometry,donor_cell,value(:gas_nv), &
          donor_guess,pv,tg,sound,good,frame)
      endif
      if(.not.good.or..not.joint) return
      if(column<=2) then
        call recover_reactive_eb_cell_3d(species,target_geometry,target_cell,value(:gas_nv), &
          target_guess,pv,tg,sound,good,frame,internal_energy_density=value(nv))
      else
        call recover_reactive_eb_cell_3d(species,donor_geometry,donor_cell,value(:gas_nv), &
          donor_guess,pv,tg,sound,good,frame,internal_energy_density=value(nv))
      endif
      if(.not.good) return
      shifted=value(:gas_nv);shift=0;shift(5)=bias(column)
      call add_rereflux_increment(shifted,shift,good)
      if(.not.good) return
      if(column<=2) then
        call recover_reactive_eb_cell_3d(species,target_geometry,target_cell,shifted, &
          target_guess,pv,tg,sound,good,frame)
      else
        call recover_reactive_eb_cell_3d(species,donor_geometry,donor_cell,shifted, &
          donor_guess,pv,tg,sound,good,frame)
      endif
    end subroutine recover_endpoint
    subroutine propose(factor,require_post,good,post_good)
      real(dp), intent(in) :: factor
      logical, intent(in) :: require_post
      logical, intent(out) :: good,post_good
      integer :: column
      good=.false.;post_good=.false.
      if(factor<=0.or.factor>=1) return
      trial=baseline
      call add_rereflux_increment(trial(:,1),factor*direction(:,1),good)
      if(.not.good) return
      trial(:,2)=trial(:,1)
      call add_rereflux_increment(trial(:,2),correction,good)
      if(.not.good) return
      trial(:,3)=(1-factor)*donor;trial(:,4)=trial(:,3)
      call add_rereflux_increment(trial(:,4),-donor_correction,good)
      if(.not.good) return
      do column=1,4
        if(joint.and.column==2.and..not.require_post) cycle
        if(trial(1,column)<=density_floor.or.any(trial(6:gas_nv,column)<0)) then
          good=.false.;return
        endif
        if(.not.require_post.and.(column==2.or.(column==1.and.raw_target))) cycle
        call recover_endpoint(trial(:,column),column,good)
        if(.not.good) return
      enddo
      call recover_endpoint(trial(:,2),2,post_good)
      if(present(native_projection)) then
        call native_projection(inventory_nodes(1),target_cell,trial(:,1),require_post,good)
        if(.not.good) return
        call native_projection(inventory_nodes(2),donor_cell,trial(:,4),.true.,good)
        if(.not.good) return
      endif
      good=.true.
    end subroutine propose
    subroutine publish(factor,post_good)
      real(dp), intent(in) :: factor
      logical, intent(in) :: post_good
      transfer=factor*donor
      if(donor_volume>1) then
        if(any(abs(transfer)>huge(1._dp)/donor_volume)) return
      endif
      transfer=transfer*donor_volume
      if(any(.not.ieee_is_finite(transfer))) return
      funded_target=trial(:,1);funded_donor=trial(:,3);packet=transfer
      complete=post_good;ok=.true.
    end subroutine publish
    subroutine margin(value,side,family,result,good)
      real(dp), intent(in) :: value(nv)
      integer, intent(in) :: side,family
      real(dp), intent(out) :: result
      logical, intent(out) :: good
      real(dp) :: scaled(nv),kinetic
      good=.false.;result=0
      if(value(1)<=density_floor) return
      if(rho_scale<1) then
        if(any(abs(value)>huge(1._dp)*rho_scale)) return
      endif
      scaled=value/rho_scale;scaled(2:4)=scaled(2:4)/sqrt(energy_scale)
      scaled(5)=scaled(5)/energy_scale
      if(scaled(1)<=0) return
      kinetic=.5_dp*sum(scaled(2:4)**2)
      if(scaled(1)<1) then
        if(kinetic>huge(1._dp)*scaled(1)) return
      endif
      result=scaled(5)+phi(2)*scaled(1)-sum(scaled(6:gas_nv)*(bounds(:,side)/energy_scale))-kinetic/scaled(1)
      if(family==2) result=result+normalized_bias(2)
      if(family==3) result=scaled(nv)/energy_scale-sum(scaled(6:gas_nv)*(bounds(:,side)/energy_scale))
      if(side==2) result=-result
      good=ieee_is_finite(result)
    end subroutine margin
    subroutine clip_quadratic(coef,starts,ends,count,good)
      real(dp), intent(in) :: coef(3)
      real(dp), intent(inout) :: starts(capacity),ends(capacity)
      integer, intent(inout) :: count
      logical, intent(out) :: good
      real(dp) :: a,b,c,scale,disc,root,other,first,last,qroot
      real(dp) :: allowed_lo(2),allowed_hi(2),next_lo(capacity),next_hi(capacity)
      integer :: allowed,piece,part,next_count
      good=.false.
      if(any(.not.ieee_is_finite(coef))) return
      scale=maxval(abs(coef));allowed=0
      if(scale==0) then
        good=.true.;return
      endif
      a=coef(3)/scale;b=coef(2)/scale;c=coef(1)/scale
      if(a==0) then
        if(b==0) then
          if(c>=0) then
            good=.true.;return
          endif
        else
          root=bounded_ratio(-c,b);allowed=1
          if(b>0) then
            allowed_lo(1)=max(0._dp,root);allowed_hi(1)=1
          else
            allowed_lo(1)=0;allowed_hi(1)=min(1._dp,root)
          endif
        endif
      else
        disc=b*b-4*a*c
        if(disc==0.and.a>0) then
          good=.true.;return
        endif
        if(disc<0) then
          if(a>0) then
            good=.true.;return
          endif
        else
          qroot=-.5_dp*(b+sign(sqrt(disc),b))
          if(qroot==0) then
            root=0;other=0
          else
            root=bounded_ratio(qroot,a);other=bounded_ratio(c,qroot)
          endif
          first=min(root,other);last=max(root,other)
          if(a>0) then
            allowed=2;allowed_lo=[0._dp,max(0._dp,last)];allowed_hi=[min(1._dp,first),1._dp]
          else
            allowed=1;allowed_lo(1)=max(0._dp,first);allowed_hi(1)=min(1._dp,last)
          endif
        endif
      endif
      next_count=0
      do piece=1,count
      do part=1,allowed
        first=max(starts(piece),allowed_lo(part));last=min(ends(piece),allowed_hi(part))
        if(first>last) cycle
        if(next_count==capacity) return
        next_count=next_count+1;next_lo(next_count)=first;next_hi(next_count)=last
      enddo
      enddo
      count=next_count
      if(count>0) then
        starts(:count)=next_lo(:count);ends(:count)=next_hi(:count)
      endif
      good=.true.
    end subroutine clip_quadratic
    real(dp) function bounded_ratio(numerator,denominator) result(value)
      real(dp), intent(in) :: numerator,denominator
      if(abs(numerator)>2*abs(denominator)) then
        value=2*sign(1._dp,numerator)*sign(1._dp,denominator)
      else
        value=numerator/denominator
      endif
    end function bounded_ratio
  end subroutine supply_reactive_eb_gas_eos_donor_3d

  ! The caller selects a deterministic owner for a shared same-level face.
  ! Replace the receiver's integrated face flux with the owner's physical
  ! transfer and correct the already advanced receiver state by the difference.
  ! Owner data is read-only; receiver flux, conserved state and T commit together.
  subroutine reconcile_reactive_eb_peer_flux_3d(species,owner,receiver,receiver_offset,owner_flux, &
      receiver_flux,state,temperature,ok,periodic_length,frame,component_count)
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(eb_geometry_3d), intent(in) :: owner,receiver
    integer, intent(in) :: receiver_offset(3)
    ! Optional domain length on exactly one axis; zero on tangential axes.
    real(dp), optional, intent(in) :: periodic_length(3)
    integer, optional, intent(in) :: component_count
    type(eb_integrated_flux_3d), intent(in) :: owner_flux
    type(eb_integrated_flux_3d), intent(inout) :: receiver_flux
    real(dp), optional, intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    logical, intent(out) :: ok
    type(eb_integrated_flux_3d) :: candidate_flux
    real(dp), allocatable :: candidate(:,:,:,:),temps(:,:,:),replacement(:)
    real(dp) :: spacing(3),rspacing(3),origin(3),expected(3),tol,oa,ra,owner_area,receiver_area,volume,sound_speed
    real(dp) :: primitive(reactive_nprim(size(species)))
    integer :: n(3),rn(3),nvar,axis,side,first(3),last(3),index(3),cell(3),oface(3),rface(3),i,j,k
    logical :: valid
    ok=.false.
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    if (.not.owner%is_valid().or..not.receiver%is_valid()) return
    n=[owner%nx,owner%ny,owner%nz];rn=[receiver%nx,receiver%ny,receiver%nz]
    if (any(receiver_offset < -rn).or.any(receiver_offset>n)) return
    if (all(receiver_offset<n.and.receiver_offset+rn>0)) return
    spacing=[owner%dx,owner%dy,owner%dz];rspacing=[receiver%dx,receiver%dy,receiver%dz]
    tol=512._dp*epsilon(1._dp)
    if (any(abs(rspacing-spacing)>tol*spacing)) return
    origin=[owner%x_lower,owner%y_lower,owner%z_lower]
    expected=origin+real(receiver_offset,dp)*spacing
    if(present(periodic_length)) then
      if(.not.all(ieee_is_finite(periodic_length))) return
      if(any(periodic_length<0).or.count(periodic_length>0)/=1) return
      axis=maxloc(periodic_length,dim=1)
      if(receiver_offset(axis)/=n(axis).and.receiver_offset(axis)/=-rn(axis)) return
      if(periodic_length(axis)<max(n(axis)*spacing(axis),rn(axis)*rspacing(axis))) return
      expected(axis)=expected(axis)-sign(periodic_length(axis),real(receiver_offset(axis),dp))
    end if
    if (any(abs([receiver%x_lower,receiver%y_lower,receiver%z_lower]-expected)> &
      tol*max(1._dp,maxval(abs(origin)),maxval(abs(expected))))) return
    nvar=reactive_nvar(size(species))
    if(present(component_count)) then
      nvar=component_count
      if(nvar<1) return
      if(nvar/=reactive_nvar(size(species)).and.present(state)) return
    endif
    if (present(state).neqv.present(temperature)) return
    if (present(state)) then
      if (any(shape(state)/=[nvar,rn]).or.any(shape(temperature)/=rn)) return
    end if
    if (.not.valid_eb_integrated_flux_3d(owner_flux,owner,nvar).or. &
      .not.valid_eb_integrated_flux_3d(receiver_flux,receiver,nvar)) return
    allocate(replacement(nvar))
    candidate_flux=receiver_flux
    if (present(state)) then
      candidate=state;temps=temperature
    end if
    do axis=1,3
      do side=0,1
        if (side==0.and.receiver_offset(axis)+rn(axis)/=0) cycle
        if (side==1.and.receiver_offset(axis)/=n(axis)) cycle
        first=max(1,receiver_offset+1);last=min(n,receiver_offset+rn)
        first(axis)=side*(n(axis)+1);last(axis)=first(axis)
        if (any(first>last)) cycle
        do k=first(3),last(3)
          do j=first(2),last(2)
            do i=first(1),last(1)
              index=[i,j,k];cell=index-receiver_offset
              oface=index;oface(axis)=side*n(axis)
              rface=cell;rface(axis)=cell(axis)-side
              oa=aperture(owner,axis,oface);ra=aperture(receiver,axis,rface)
              if (abs(oa-ra)>tol) return
              if ((oa>0._dp).neqv.(ra>0._dp)) return
              if (oa<=0._dp) cycle
              if(present(periodic_length)) then
                if(any(abs(face_centroid(owner,axis,oface)-face_centroid(receiver,axis,rface))>tol)) return
              end if
              if (receiver%volume_fraction(cell(1),cell(2),cell(3))<=0._dp) return
              owner_area=oa*product(spacing)/spacing(axis)
              receiver_area=ra*product(rspacing)/rspacing(axis)
              volume=receiver%volume_fraction(cell(1),cell(2),cell(3))*product(rspacing)
              replacement=(owner_area/receiver_area)*face(owner_flux,axis,oface)
              if (present(state)) then
                candidate(:,cell(1),cell(2),cell(3))=state(:,cell(1),cell(2),cell(3))+ &
                  real(2*side-1,dp)*(receiver_area/volume)*(replacement-face(receiver_flux,axis,rface))
                if (.not.ieee_is_finite(temperature(cell(1),cell(2),cell(3)))) return
                if (temperature(cell(1),cell(2),cell(3))<=0._dp) return
                ! Accepted states must satisfy the same nonnegative species contract as restart.
                if(any(candidate(6:,cell(1),cell(2),cell(3))<0)) return
                call recover_reactive_eb_cell_3d(species,receiver,cell,candidate(:,cell(1),cell(2),cell(3)), &
                  temperature(cell(1),cell(2),cell(3)),primitive,temps(cell(1),cell(2),cell(3)),sound_speed,valid,frame)
                if (.not.valid) return
              end if
              select case(axis)
              case(1);candidate_flux%x(:,rface(1),rface(2),rface(3))=replacement
              case(2);candidate_flux%y(:,rface(1),rface(2),rface(3))=replacement
              case(3);candidate_flux%z(:,rface(1),rface(2),rface(3))=replacement
              end select
            end do
          end do
        end do
      end do
    end do
    receiver_flux=candidate_flux
    if (present(state)) then
      state=candidate;temperature=temps
    end if
    ok=.true.
  end subroutine reconcile_reactive_eb_peer_flux_3d

  pure function face_centroid(g,axis,index) result(value)
    type(eb_geometry_3d), intent(in) :: g
    integer, intent(in) :: axis,index(3)
    real(dp) :: value(2)
    integer :: i,j,k
    i=index(1);j=index(2);k=index(3)
    select case(axis)
    case(1);value=[g%x_face_centroid_y(i,j,k),g%x_face_centroid_z(i,j,k)]
    case(2);value=[g%y_face_centroid_x(i,j,k),g%y_face_centroid_z(i,j,k)]
    case(3);value=[g%z_face_centroid_x(i,j,k),g%z_face_centroid_y(i,j,k)]
    end select
  end function face_centroid

  subroutine restrict_eb_integrated_flux_3d(fine,ratio,fine_flux,coarse_flux,ok)
    type(eb_geometry_3d), intent(in) :: fine
    integer, intent(in) :: ratio
    type(eb_integrated_flux_3d), intent(in) :: fine_flux
    type(eb_integrated_flux_3d), intent(inout) :: coarse_flux
    logical, intent(out) :: ok
    type(eb_geometry_3d) :: coarse
    type(eb_integrated_flux_3d) :: candidate
    integer :: axis,i,j,k,ii,jj,kk,lo(3),hi(3),first(3),dims(3),index(3),nvar
    real(dp) :: area_fraction,weight
    real(dp), allocatable :: average(:)
    logical :: valid
    ok=.false.
    call coarsen_eb_geometry_3d(fine,ratio,coarse,valid)
    if (.not.valid) return
    if (.not.allocated(fine_flux%x)) return
    nvar=size(fine_flux%x,1)
    if (.not.valid_eb_integrated_flux_3d(fine_flux,fine,nvar)) return
    dims=[coarse%nx,coarse%ny,coarse%nz]
    allocate(candidate%x(nvar,0:dims(1),dims(2),dims(3)), &
      candidate%y(nvar,dims(1),0:dims(2),dims(3)),candidate%z(nvar,dims(1),dims(2),0:dims(3)),average(nvar))
    do axis=1,3
      first=1;first(axis)=0
      do k=first(3),dims(3)
        do j=first(2),dims(2)
          do i=first(1),dims(1)
            index=[i,j,k];lo=(index-1)*ratio+1;hi=index*ratio;lo(axis)=hi(axis)
            average=0._dp;area_fraction=0._dp
            do kk=lo(3),hi(3)
              do jj=lo(2),hi(2)
                do ii=lo(1),hi(1)
                  weight=aperture(fine,axis,[ii,jj,kk])
                  area_fraction=area_fraction+weight
                  average=average+weight*face(fine_flux,axis,[ii,jj,kk])
                end do
              end do
            end do
            if (area_fraction>0._dp) average=average/area_fraction
            select case(axis)
            case(1)
              candidate%x(:,i,j,k)=average
            case(2)
              candidate%y(:,i,j,k)=average
            case(3)
              candidate%z(:,i,j,k)=average
            end select
          end do
        end do
      end do
    end do
    if (.not.valid_eb_integrated_flux_3d(candidate,coarse,nvar)) return
    coarse_flux=candidate;ok=.true.
  end subroutine restrict_eb_integrated_flux_3d

  ! Reconcile one coarse interface and its fine subfaces using both endpoint
  ! inventories. Transfers are time-integrated OPEN-AREA inventories oriented
  ! into the exposed coarse cell; the fine cells see the opposite transfer.
  ! weights partition the original coarse transfer over open fine subfaces.
  ! One theta applies to every conserved component and every subface, so the
  ! corrected coarse/fine states and face transfers remain conservative.
  ! An empty admissible interval rejects without publishing any output. That
  ! case needs a wider donor network; a two-cell adjustment cannot repair it.
  ! EOS acceptance and interactions with other faces remain caller duties.
  subroutine reconcile_reactive_eb_interface_budget_3d(coarse_state,coarse_volume,coarse_transfer, &
      fine_state,fine_volume,fine_transfer,weights,corrected_transfer,coarse_increment,fine_increment,theta,ok, &
      minimum_theta,maximum_theta, &
      passive_components)
    real(dp), intent(in) :: coarse_state(:),coarse_volume,coarse_transfer(:)
    real(dp), intent(in) :: fine_state(:,:),fine_volume(:),fine_transfer(:,:),weights(:)
    real(dp), allocatable, intent(inout) :: corrected_transfer(:,:),coarse_increment(:),fine_increment(:,:)
    real(dp), intent(inout) :: theta
    logical, intent(out) :: ok
    ! Trailing independent scalars share every conservative packet and face
    ! blend, but do not impose species nonnegativity on its gas allowance.
    integer, optional, intent(in) :: passive_components
    integer :: active_components,passives
    ! Descendant inventories may impose a stronger floor than a parent average.
    real(dp), optional, intent(in) :: minimum_theta,maximum_theta
    real(dp), allocatable :: baseline(:,:),difference(:,:),net(:),faces(:,:),dc(:),df(:,:),endpoint(:)
    real(dp) :: lower,upper,budget,ratio,factor,floor,ceiling
    integer :: nv,nf,c,q,ios
    logical :: admissible
    ok=.false.;nv=size(coarse_state);nf=size(weights)
    passives=0
    if(present(passive_components)) passives=passive_components
    if(passives<0.or.passives>nv-6) return
    active_components=nv-passives
    floor=0;ceiling=1
    if(present(minimum_theta)) then
      if(.not.ieee_is_finite(minimum_theta)) return
      if(minimum_theta<0.or.minimum_theta>1) return
      floor=minimum_theta
    endif
    if(present(maximum_theta)) then
      if(.not.ieee_is_finite(maximum_theta)) return
      if(maximum_theta<0.or.maximum_theta>1) return
      ceiling=maximum_theta
    end if
    if(ceiling<floor) return
    if(nv<6.or.nf<1) return
    if(size(coarse_transfer)/=nv.or.any(shape(fine_state)/=[nv,nf]).or. &
      any(shape(fine_transfer)/=[nv,nf]).or.size(fine_volume)/=nf) return
    if(any(.not.ieee_is_finite(coarse_state)).or.any(.not.ieee_is_finite(coarse_transfer)).or. &
      any(.not.ieee_is_finite(fine_state)).or.any(.not.ieee_is_finite(fine_transfer))) return
    if(.not.ieee_is_finite(coarse_volume).or.coarse_volume<=0) return
    if(any(.not.ieee_is_finite(fine_volume)).or.any(fine_volume<=0)) return
    if(any(.not.ieee_is_finite(weights)).or.any(weights<0)) return
    if(abs(sum(weights)-1._dp)>32*epsilon(1._dp)*real(nf,dp)) return
    if(coarse_state(1)<=0.or.any(coarse_state(6:active_components)<0)) return
    if(any(fine_state(1,:)<=0).or.any(fine_state(6:active_components,:)<0)) return
    do q=1,nf
      if(weights(q)==0.and.any(fine_transfer(:,q)/=0)) return
    enddo
    allocate(baseline(nv,nf),difference(nv,nf),net(nv),faces(nv,nf),dc(nv),df(nv,nf),endpoint(nv),stat=ios)
    if(ios/=0) return
    do q=1,nf
      baseline(:,q)=weights(q)*coarse_transfer
    enddo
    difference=fine_transfer-baseline;net=sum(difference,dim=2)
    if(any(.not.ieee_is_finite(difference)).or.any(.not.ieee_is_finite(net))) return
    factor=ceiling
    call represented_budget(factor,admissible)
    if(.not.admissible) then
      lower=floor;upper=ceiling
      do c=1,active_components
        if(c>=2.and.c<=5) cycle
        if(net(c)<0) then
          budget=coarse_volume*coarse_state(c)
          if(.not.ieee_is_finite(budget)) return
          if(budget<=-net(c)) upper=min(upper,budget/(-net(c)))
        endif
        do q=1,nf
          if(difference(c,q)>=0) cycle
          budget=fine_volume(q)*fine_state(c,q)
          if(.not.ieee_is_finite(budget)) return
          if(budget<=-difference(c,q)) then
            ratio=budget/(-difference(c,q))
            lower=max(lower,1-ratio)
          endif
        enddo
      enddo
      if(lower>upper) return
      factor=upper
      if(upper<1.and.lower<upper) factor=upper-min(16*epsilon(1._dp)*upper,.25_dp*(upper-lower))
      call represented_budget(factor,admissible)
      if(.not.admissible.and.lower<upper) then
        ! Stay strictly inside both analytic donor bounds when the face
        ! representation loses the margin at an interval endpoint.
        factor=lower+.5_dp*(upper-lower)
        call represented_budget(factor,admissible)
      endif
      if(.not.admissible) return
    endif
    call move_alloc(faces,corrected_transfer);call move_alloc(dc,coarse_increment)
    call move_alloc(df,fine_increment);theta=factor;ok=.true.
  contains
    subroutine represented_budget(a,good)
      real(dp), intent(in) :: a
      logical, intent(out) :: good
      if(a==1) then
        faces=fine_transfer
      else if(a==0) then
        faces=baseline
      else
        faces=baseline+a*difference
      endif
      ! Compute actual published differences, rather than assuming the
      ! analytic blend survived rounding in the integrated face arrays.
      dc=(sum(faces,dim=2)-coarse_transfer)/coarse_volume
      do q=1,nf
        df(:,q)=(fine_transfer(:,q)-faces(:,q))/fine_volume(q)
      enddo
      good=a>=floor.and.all(ieee_is_finite(faces)).and.all(ieee_is_finite(dc)).and.all(ieee_is_finite(df))
      endpoint=coarse_state+dc
      good=good.and.all(ieee_is_finite(endpoint)).and.endpoint(1)>0.and.all(endpoint(6:active_components)>=0)
      do q=1,nf
        endpoint=fine_state(:,q)+df(:,q)
        good=good.and.all(ieee_is_finite(endpoint)).and.endpoint(1)>0.and.all(endpoint(6:active_components)>=0)
      enddo
    end subroutine represented_budget
  end subroutine reconcile_reactive_eb_interface_budget_3d

  ! Supply an interface deficit by a conservative exchange with one connected
  ! interior donor per fine subface. The smallest convex mixing fraction that
  ! funds the coarse-admissible interface transfer is selected for each pair.
  ! All conserved components move together, with the returned inventory used
  ! as an internal face transfer. Inputs and all outputs roll back on failure.
  ! A donor that cannot fund its pair requires a wider neighborhood.
  ! target_subface funds one endpoint for sequential searches. The caller must
  ! finish funding and reconcile the complete interface before accepting it.
  ! allow_partial accepts an improving convex exchange while an endpoint
  ! deficit remains. The full interface still requires strict reconciliation.
  subroutine supply_reactive_eb_interface_donors_3d(coarse_state,coarse_volume,coarse_transfer, &
      fine_state,fine_volume,fine_transfer,weights,donor_state,donor_volume,funded_fine,funded_donor,exchange,ok, &
      target_subface,allow_partial,force_exchange, &
      passive_components)
    real(dp), intent(in) :: coarse_state(:),coarse_volume,coarse_transfer(:),fine_state(:,:),fine_volume(:)
    real(dp), intent(in) :: fine_transfer(:,:),weights(:),donor_state(:,:),donor_volume(:)
    real(dp), allocatable, intent(inout) :: funded_fine(:,:),funded_donor(:,:),exchange(:,:)
    logical, intent(out) :: ok
    ! Trailing independent scalars share every conservative packet and face
    ! blend, but do not impose species nonnegativity on its gas allowance.
    integer, optional, intent(in) :: passive_components
    integer :: active_components,passives
    integer, optional, intent(in) :: target_subface
    logical, optional, intent(in) :: allow_partial,force_exchange
    logical :: partial,improved,unresolved
    real(dp), allocatable :: baseline(:,:),difference(:,:),net(:),target(:,:),candidate(:,:),donors(:,:),transfer(:,:)
    real(dp), allocatable :: desired(:),gain(:),pair(:),before(:)
    real(dp) :: factor,budget,lambda,required,volume_sum,cap
    integer :: nv,nf,c,q,ios
    ok=.false.;nv=size(coarse_state);nf=size(weights)
    passives=0
    if(present(passive_components)) passives=passive_components
    if(passives<0.or.passives>nv-6) return
    active_components=nv-passives
    partial=.false.
    if(present(allow_partial)) partial=allow_partial
    if(nv<6.or.nf<1) return
    if(present(target_subface)) then
      if(target_subface<1.or.target_subface>nf) return
    endif
    if(size(coarse_transfer)/=nv.or.any(shape(fine_state)/=[nv,nf]).or.any(shape(fine_transfer)/=[nv,nf]).or. &
      any(shape(donor_state)/=[nv,nf]).or.size(fine_volume)/=nf.or.size(donor_volume)/=nf) return
    if(.not.ieee_is_finite(coarse_volume).or.coarse_volume<=0.or.any(.not.ieee_is_finite(fine_volume)).or. &
      any(fine_volume<=0).or.any(.not.ieee_is_finite(donor_volume)).or.any(donor_volume<=0)) return
    if(any(.not.ieee_is_finite(coarse_state)).or.any(.not.ieee_is_finite(coarse_transfer)).or. &
      any(.not.ieee_is_finite(fine_state)).or.any(.not.ieee_is_finite(fine_transfer)).or. &
      any(.not.ieee_is_finite(donor_state)).or.any(.not.ieee_is_finite(weights))) return
    if(any(weights<0).or.abs(sum(weights)-1._dp)>32*epsilon(1._dp)*real(nf,dp)) return
    if(coarse_state(1)<=0.or.any(coarse_state(6:active_components)<0).or.any(fine_state(1,:)<=0).or. &
      any(fine_state(6:active_components,:)<0).or.any(donor_state(1,:)<=0).or.any(donor_state(6:active_components,:)<0)) return
    allocate(baseline(nv,nf),difference(nv,nf),net(nv),target(nv,nf),candidate(nv,nf),donors(nv,nf), &
      transfer(nv,nf),desired(nv),gain(nv),pair(nv),before(nv),stat=ios)
    if(ios/=0) return
    do q=1,nf
      baseline(:,q)=weights(q)*coarse_transfer
    enddo
    difference=fine_transfer-baseline;net=sum(difference,dim=2)
    if(any(.not.ieee_is_finite(difference)).or.any(.not.ieee_is_finite(net))) return
    factor=1
    do c=1,active_components
      if(c>=2.and.c<=5) cycle
      if(net(c)>=0) cycle
      budget=coarse_volume*coarse_state(c)
      if(.not.ieee_is_finite(budget)) return
      if(budget<=-net(c)) factor=min(factor,budget/(-net(c)))
    enddo
    if(factor<1.and.factor>0) factor=factor*(1-32*epsilon(1._dp))
    target=baseline+factor*difference;candidate=fine_state;donors=donor_state;transfer=0
    do q=1,nf
      if(present(target_subface)) then
        if(q/=target_subface) cycle
      endif
      desired=fine_state(:,q)+(fine_transfer(:,q)-target(:,q))/fine_volume(q)
      before=desired
      volume_sum=fine_volume(q)+donor_volume(q)
      if(.not.ieee_is_finite(volume_sum)) return
      ! Inventory of the full equalization exchange, evaluated without the
      ! product of two cell volumes to avoid unnecessary under/overflow.
      pair=(fine_volume(q)/volume_sum)*donor_volume(q)*(donor_state(:,q)-fine_state(:,q))
      gain=pair/fine_volume(q);lambda=0;unresolved=.false.
      ! A covered endpoint can have an admissible average while a physical
      ! descendant needs supply. The tree caller checks that full exchange
      ! direction against its physical leaves and clips it if necessary.
      if(present(force_exchange)) then
        if(force_exchange) lambda=1
      endif
      do c=1,active_components
        if(c>=2.and.c<=5) cycle
        if(desired(c)>=0.and.(c/=1.or.desired(c)>0)) cycle
        if(gain(c)<=0) then
          if(.not.partial.or.gain(c)<0) return
          unresolved=.true.;cycle
        endif
        required=-desired(c)/gain(c)
        if(c==1.and.desired(c)==0) required=32*epsilon(1._dp)*fine_state(c,q)/gain(c)
        lambda=max(lambda,required)
      enddo
      if(.not.ieee_is_finite(lambda)) return
      if(lambda>1) then
        if(.not.partial) return
        lambda=1
      endif
      if(lambda>0) lambda=min(1._dp,lambda*(1+32*epsilon(1._dp)))
      if(partial) then
        cap=1
        do c=1,active_components
          if(c>=2.and.c<=5) cycle
          if(before(c)<0.or.gain(c)>=0) cycle
          if(before(c)>-gain(c)) cycle
          required=before(c)/(-gain(c))
          if(c==1.and.required>0) required=nearest(required,-1._dp)
          cap=min(cap,max(0._dp,required))
        enddo
        if(unresolved) lambda=1
        lambda=min(lambda,cap)
      endif
      transfer(:,q)=lambda*pair
      candidate(:,q)=fine_state(:,q)+transfer(:,q)/fine_volume(q)
      donors(:,q)=donor_state(:,q)-transfer(:,q)/donor_volume(q)
      desired=candidate(:,q)+(fine_transfer(:,q)-target(:,q))/fine_volume(q)
      if(any(.not.ieee_is_finite(candidate(:,q))).or.any(.not.ieee_is_finite(donors(:,q))).or. &
        any(.not.ieee_is_finite(desired)).or.candidate(1,q)<=0.or.any(candidate(6:active_components,q)<0).or. &
        donors(1,q)<=0.or.any(donors(6:active_components,q)<0)) return
      if(desired(1)<=0.or.any(desired(6:active_components)<0)) then
        if(.not.partial) return
        improved=.false.
        do c=1,active_components
          if(c>=2.and.c<=5) cycle
          if(before(c)<0.or.(c==1.and.before(c)==0)) then
            if(desired(c)<before(c)) return
            if(desired(c)>before(c)) improved=.true.
          else if(desired(c)<0.or.(c==1.and.desired(c)==0)) then
            return
          endif
        enddo
        if(.not.improved) return
      endif
    enddo
    call move_alloc(candidate,funded_fine);call move_alloc(donors,funded_donor);call move_alloc(transfer,exchange)
    ok=.true.
  end subroutine supply_reactive_eb_interface_donors_3d

  ! Apply fine-minus-coarse integrated flux only at interior interfaces between
  ! refined and unrefined coarse cells. Fine states are updated/averaged separately.
  ! With all_faces, apply the full flux difference on exposed cells, including
  ! outer faces. The caller then owns the synchronized outer-boundary budget.
  ! Assemble signed conservative reflux increments without applying them to
  ! a thermodynamic state. Composite donor budgets and redistribution must
  ! operate on these increments before final state acceptance. Neither flux
  ! nor input state is modified; output storage commits only on success.
  subroutine compute_reactive_eb_reflux_increment_3d(geometry,refined,coarse_flux,fine_flux,increment,ok,all_faces,touched, &
      redistribution_delta)
    real(dp), optional, intent(in) :: redistribution_delta(:,:,:,:)
    type(eb_geometry_3d), intent(in) :: geometry
    logical, intent(in) :: refined(:,:,:)
    type(eb_integrated_flux_3d), intent(in) :: coarse_flux,fine_flux
    real(dp), allocatable, intent(inout) :: increment(:,:,:,:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: all_faces
    logical, allocatable, optional, intent(inout) :: touched(:,:,:)
    real(dp), allocatable :: candidate(:,:,:,:),correction(:)
    logical, allocatable :: mask(:,:,:)
    real(dp) :: spacing(3),volume,scale
    integer :: i,j,k,axis,step,idx(3),neighbor(3),face_idx(3),dims(3),nvar
    logical :: complete,changed
    ok=.false.;complete=.false.
    if(present(all_faces)) complete=all_faces
    if(.not.geometry%is_valid().or..not.allocated(coarse_flux%x)) return
    dims=[geometry%nx,geometry%ny,geometry%nz];nvar=size(coarse_flux%x,1)
    if(nvar<1.or.any(shape(refined)/=dims)) return
    if(.not.valid_eb_integrated_flux_3d(coarse_flux,geometry,nvar).or. &
      .not.valid_eb_integrated_flux_3d(fine_flux,geometry,nvar)) return
    if(present(redistribution_delta)) then
      if(any(shape(redistribution_delta)/=[nvar,dims])) return
      if(any(.not.ieee_is_finite(redistribution_delta))) return
      do k=1,dims(3)
      do j=1,dims(2)
      do i=1,dims(1)
        ! Fine-covered corrections must be exported to fine states separately.
        ! Reject them here instead of dropping a register contribution.
        if(.not.refined(i,j,k).and.geometry%volume_fraction(i,j,k)>0) cycle
        if(any(redistribution_delta(:,i,j,k)/=0)) return
      enddo
      enddo
      enddo
    endif
    allocate(candidate(nvar,dims(1),dims(2),dims(3)),mask(dims(1),dims(2),dims(3)),correction(nvar))
    candidate=0;mask=.false.
    spacing=[geometry%dx,geometry%dy,geometry%dz]
    do k=1,dims(3)
      do j=1,dims(2)
        do i=1,dims(1)
          if (refined(i,j,k).or.geometry%volume_fraction(i,j,k)<=0._dp) cycle
          idx=[i,j,k];correction=0._dp;changed=.false.
          volume=geometry%volume_fraction(i,j,k)*product(spacing)
          if(.not.ieee_is_finite(volume).or.volume<=0) return
          do axis=1,3
            do step=-1,1,2
              neighbor=idx;neighbor(axis)=neighbor(axis)+step
              if (.not.complete) then
                if (any(neighbor<1).or.any(neighbor>dims)) cycle
                if (.not.refined(neighbor(1),neighbor(2),neighbor(3))) cycle
              end if
              face_idx=idx;face_idx(axis)=face_idx(axis)+min(step,0)
              if (complete) then
                if (all(face(fine_flux,axis,face_idx)==face(coarse_flux,axis,face_idx))) cycle
              end if
              scale=aperture(geometry,axis,face_idx)
              if (scale<=0._dp) cycle
              scale=-real(step,dp)*scale*(product(spacing)/spacing(axis))/volume
              correction=correction+scale*(face(fine_flux,axis,face_idx)-face(coarse_flux,axis,face_idx))
              changed=.true.
            end do
          end do
          if(present(redistribution_delta)) then
            do axis=1,nvar
              if(redistribution_delta(axis,i,j,k)>0) then
                if(correction(axis)>huge(1._dp)-redistribution_delta(axis,i,j,k)) return
              else if(redistribution_delta(axis,i,j,k)<0) then
                if(correction(axis)<-huge(1._dp)-redistribution_delta(axis,i,j,k)) return
              endif
            enddo
            correction=correction+redistribution_delta(:,i,j,k)
            changed=changed.or.any(redistribution_delta(:,i,j,k)/=0)
          endif
          candidate(:,i,j,k)=correction
          mask(i,j,k)=changed
        end do
      end do
    end do
    if(any(.not.ieee_is_finite(candidate))) return
    call move_alloc(candidate,increment)
    if(present(touched)) call move_alloc(mask,touched)
    ok=.true.
  end subroutine compute_reactive_eb_reflux_increment_3d

  subroutine reflux_reactive_eb_state_3d(species,geometry,refined,coarse_flux,fine_flux,state,temperature,ok,all_faces,frame, &
      redistribution_delta)
    real(dp), optional, intent(in) :: redistribution_delta(:,:,:,:)
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(eb_geometry_3d), intent(in) :: geometry
    logical, intent(in) :: refined(:,:,:)
    type(eb_integrated_flux_3d), intent(in) :: coarse_flux,fine_flux
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: all_faces
    real(dp), allocatable :: candidate(:,:,:,:),temps(:,:,:),primitive(:),increment(:,:,:,:)
    logical, allocatable :: touched(:,:,:)
    real(dp) :: sound_speed
    integer :: i,j,k,dims(3),nvar
    logical :: valid,complete
    ok=.false.;complete=.false.
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    if (present(all_faces)) complete=all_faces
    if (.not.geometry%is_valid()) return
    dims=[geometry%nx,geometry%ny,geometry%nz];nvar=reactive_nvar(size(species))
    if (any(shape(state)/=[nvar,dims]).or.any(shape(temperature)/=dims).or.any(shape(refined)/=dims)) return
    if (.not.valid_eb_integrated_flux_3d(coarse_flux,geometry,nvar).or. &
      .not.valid_eb_integrated_flux_3d(fine_flux,geometry,nvar)) return
    allocate(candidate,source=state);allocate(temps,source=temperature)
    allocate(primitive(reactive_nprim(size(species))))
    call compute_reactive_eb_reflux_increment_3d(geometry,refined,coarse_flux,fine_flux,increment,valid,complete, &
      touched,redistribution_delta)
    if(.not.valid) return
    do k=1,dims(3)
      do j=1,dims(2)
        do i=1,dims(1)
          if(.not.touched(i,j,k)) cycle
          if (.not.ieee_is_finite(temperature(i,j,k))) return
          if (temperature(i,j,k)<=0._dp) return
          candidate(:,i,j,k)=state(:,i,j,k)+increment(:,i,j,k)
          if(any(candidate(6:,i,j,k)<0)) return
          call recover_reactive_eb_cell_3d(species,geometry,[i,j,k],candidate(:,i,j,k),temperature(i,j,k), &
            primitive,temps(i,j,k),sound_speed,valid,frame)
          if (.not.valid) return
        end do
      end do
    end do
    state=candidate;temperature=temps;ok=.true.
  end subroutine reflux_reactive_eb_state_3d

  logical function valid_eb_integrated_flux_3d(flux,geometry,nvar) result(valid)
    type(eb_integrated_flux_3d), intent(in) :: flux
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: nvar
    valid=.false.
    if (nvar<1.or..not.allocated(flux%x).or..not.allocated(flux%y).or..not.allocated(flux%z)) return
    if (any(shape(flux%x)/=[nvar,geometry%nx+1,geometry%ny,geometry%nz]).or. &
      any(shape(flux%y)/=[nvar,geometry%nx,geometry%ny+1,geometry%nz]).or. &
      any(shape(flux%z)/=[nvar,geometry%nx,geometry%ny,geometry%nz+1])) return
    if (any(lbound(flux%x)/=[1,0,1,1]).or.any(lbound(flux%y)/=[1,1,0,1]).or. &
      any(lbound(flux%z)/=[1,1,1,0])) return
    valid=all(ieee_is_finite(flux%x)).and.all(ieee_is_finite(flux%y)).and.all(ieee_is_finite(flux%z))
  end function valid_eb_integrated_flux_3d

  function face(flux,axis,index) result(value)
    type(eb_integrated_flux_3d), intent(in) :: flux
    integer, intent(in) :: axis,index(3)
    real(dp) :: value(size(flux%x,1))
    select case(axis)
    case(1)
      value=flux%x(:,index(1),index(2),index(3))
    case(2)
      value=flux%y(:,index(1),index(2),index(3))
    case(3)
      value=flux%z(:,index(1),index(2),index(3))
    end select
  end function face

  real(dp) function aperture(geometry,axis,index) result(value)
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: axis,index(3)
    select case(axis)
    case(1)
      value=geometry%x_face_fraction(index(1),index(2),index(3))
    case(2)
      value=geometry%y_face_fraction(index(1),index(2),index(3))
    case(3)
      value=geometry%z_face_fraction(index(1),index(2),index(3))
    end select
  end function aperture
end module eb_reactive_reflux_3d_mod
