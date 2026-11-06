module eb_reactive_transfer_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive, &
    reactive_primitive_to_conserved
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d,eb_cell_connectivity_3d
  implicit none
  private
  public :: face_metrics,recover_reactive_eb_cell_3d
  public :: restrict_reactive_eb_state_3d,prolong_reactive_eb_state_pcm_3d
  public :: fill_reactive_eb_patch_halo_3d,overlay_reactive_eb_peer_halo_3d
  public :: plan_eb_linear_halo_3d,interpolate_eb_linear_halo_3d
  public :: reconstruct_eb_pressure_state_3d,interpolate_eb_pressure_time_3d
contains
  ! Nonrotating auxiliary pressure representation for halo/covered states.
  ! This is not a conservative update of a valid leaf cell. The caller owns
  ! that distinction and any energy ledger. Preserve rho, momentum, species
  ! and the total-minus-held energy gap; publish only after EOS admission.
  ! All inout arguments, including the energy correction, survive refusal.
  subroutine reconstruct_eb_pressure_state_3d(species,state,held,temperature,pressure,energy_delta,ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:),held,temperature,energy_delta
    real(dp), intent(in) :: pressure
    logical, intent(out) :: ok
    real(dp) :: q(reactive_nprim(size(species))),rebuilt(size(state)),candidate(size(state))
    real(dp) :: t,cs,kinetic,numerator,delta,new_held,readback_t
    real(dp), parameter :: bound=huge(1._dp)/16
    logical :: good
    ok=.false.
    if(size(state)/=reactive_nvar(size(species))) return
    if(.not.all(ieee_is_finite(state))) return
    if(.not.all(ieee_is_finite([held,temperature,pressure]))) return
    if(pressure<=0.or.temperature<=0.or.state(1)<=0) return
    ! Bound intermediates before arithmetic, including under FPE trapping.
    if(any(abs(state)>bound).or.abs(held)>bound) return
    if(any(abs(state(2:4))>sqrt(bound))) return
    numerator=.5_dp*sum(state(2:4)**2)
    if(state(1)<1) then
      if(numerator>bound*state(1)) return
    endif
    kinetic=numerator/state(1)
    call reactive_conserved_to_primitive(species,state,temperature,q,t,cs,good,internal_energy_density=held)
    if(.not.good) return
    q(5)=pressure
    call reactive_primitive_to_conserved(species,q,rebuilt,t,cs,good)
    if(.not.good) return
    if(abs(rebuilt(5))>bound) return
    delta=rebuilt(5)-kinetic-held
    candidate=state;candidate(5)=state(5)+delta;new_held=held+delta
    call reactive_conserved_to_primitive(species,candidate,t,q,readback_t,cs,good, &
      internal_energy_density=new_held)
    if(.not.good) return
    state(5)=candidate(5);held=new_held;temperature=t;energy_delta=delta;ok=.true.
  end subroutine reconstruct_eb_pressure_state_3d

  ! Interpolate pressure on the parent clock, not EOS of mixed conserved
  ! states. Endpoint gas and held energies must describe the same clock.
  subroutine interpolate_eb_pressure_time_3d(species,left,right,left_held,right_held, &
      left_t,right_t,fraction,pressure,ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: left(:),right(:),left_held,right_held,left_t,right_t,fraction
    real(dp), intent(inout) :: pressure
    logical, intent(out) :: ok
    real(dp) :: q(reactive_nprim(size(species))),t,cs,p0,p1
    logical :: good
    ok=.false.
    if(.not.ieee_is_finite(fraction)) return
    if(fraction<0.or.fraction>1) return
    call reactive_conserved_to_primitive(species,left,left_t,q,t,cs,good,internal_energy_density=left_held)
    if(.not.good) return
    p0=q(5)
    call reactive_conserved_to_primitive(species,right,right_t,q,t,cs,good,internal_energy_density=right_held)
    if(.not.good) return
    p1=q(5)
    pressure=(1-fraction)*p0+fraction*p1;ok=.true.
  end subroutine interpolate_eb_pressure_time_3d

  ! Plan the six state reads for the pinned EB mf_lincc interpolater. The
  ! geometry includes a halo, so patch edges do not artificially flatten a
  ! connected stencil. State requests remain at the donor's own level.
  ! A missing same-level neighbor or inlet stencil retains PCM; recursive
  ! spatial interpolation and ext_dir slopes are not provided by this path.
  subroutine plan_eb_linear_halo_3d(levels,shapes,origin,domain,donor,source,fine_level,fine_index, &
      geometry,width,faces,requests,signs,offset,linear,ok,connected_directions)
    integer, intent(in) :: levels(:),shapes(:,:)
    integer(int64), intent(in) :: origin(:,:),domain(:,:),fine_index(3)
    integer, intent(in) :: donor,source(3),fine_level,width
    type(eb_geometry_3d), intent(in) :: geometry
    character(len=*), intent(in) :: faces(6)
    integer, intent(out) :: requests(4,6),signs(3,6)
    real(dp), intent(out) :: offset(3)
    logical, intent(out) :: linear,ok
    logical, optional, intent(in) :: connected_directions
    logical :: connected(3,3,3),good,directional,enabled(3)
    integer(int64) :: query(3),local(3),factor(3)
    integer :: depth,a,s,q,b,face,p,minus_index(3),plus_index(3),ii,jj,kk
    ok=.false.;linear=.false.;requests=0;signs=1;offset=0
    if(donor<1.or.donor>size(levels)) return
    if(any(shape(shapes)/=[3,size(levels)]).or.any(shape(origin)/=[3,size(levels)])) return
    if(size(domain,1)/=3.or.fine_level<1.or.fine_level>size(domain,2)) return
    if(any(shapes<=0).or.any(domain<=0)) return
    if(any(source<1).or.any(source>shapes(:,donor))) return
    depth=levels(donor)+1
    if(depth<1.or.depth>size(domain,2)) return
    if(depth>=fine_level) then
      ok=.true.;return
    endif
    if(width<1) return
    if(any([geometry%nx,geometry%ny,geometry%nz]/=shapes(:,donor)+2*width)) return
    call eb_cell_connectivity_3d(geometry,source+width,connected,good)
    if(.not.good) return
    ok=.true.
    directional=.false.
    if(present(connected_directions)) directional=connected_directions
    if(.not.directional.and..not.all(connected)) return
    enabled=.true.
    if(directional) then
      if(.not.connected(2,2,2)) return
      do a=1,3
        minus_index=2;plus_index=2;minus_index(a)=1;plus_index(a)=3
        enabled(a)=connected(minus_index(1),minus_index(2),minus_index(3)).and. &
          connected(plus_index(1),plus_index(2),plus_index(3))
      enddo
      if(.not.any(enabled)) return
      ! Admit only a connected stencil in the retained subspace. Axial reads
      ! alone must not bypass a disconnected diagonal around a corner.
      do kk=1,3
      do jj=1,3
      do ii=1,3
        if(.not.enabled(1).and.ii/=2) cycle
        if(.not.enabled(2).and.jj/=2) cycle
        if(.not.enabled(3).and.kk/=2) cycle
        if(.not.connected(ii,jj,kk)) return
      enddo
      enddo
      enddo
    endif
    factor=domain(:,fine_level)/domain(:,depth)
    offset=(real(modulo(fine_index,factor),dp)+.5_dp)/real(factor,dp)-.5_dp
    do a=1,3
    if(directional) then
      ! Opt-in: a blocked axis contributes zero slope without sampling solid
      ! state. Connected tangential axes retain their existing donor reads.
      if(.not.enabled(a)) then
        requests(:,2*a-1)=[donor,source];requests(:,2*a)=[donor,source]
        cycle
      endif
    endif
    do s=1,2
      q=2*a-2+s
      query=origin(:,donor)+int(source-1,int64);query(a)=query(a)+2*s-3
      do b=1,3
        if(faces(2*b-1)=='periodic') then
          query(b)=modulo(query(b),domain(b,depth));cycle
        endif
        if(query(b)>=0.and.query(b)<domain(b,depth)) cycle
        face=2*b
        if(query(b)<0) face=face-1
        if(faces(face)=='inflow') return
        if(faces(face)=='slip_wall') signs(b,q)=-1
        if(faces(face)=='no_slip_wall'.or.faces(face)=='isothermal_wall') signs(:,q)=-1
        ! Only one adjacent coarse cell is requested. Reflection and
        ! foextrap therefore select the same boundary-adjacent cell.
        query(b)=max(0_int64,min(domain(b,depth)-1,query(b)))
      enddo
      do p=1,size(levels)
        if(levels(p)+1/=depth) cycle
        local=query-origin(:,p)
        if(any(local<0).or.any(local>=int(shapes(:,p),int64))) cycle
        requests(:,q)=[p,int(local)+1];exit
      enddo
      if(requests(1,q)==0) return
    enddo
    enddo
    linear=.true.
  end subroutine plan_eb_linear_halo_3d

  ! Default: one limiter per direction, shared by gas, temperature and stored rho e.
  ! The opt-in candidate limits momentum components separately, while retaining
  ! one scalar factor for density, species, energies and temperature. This avoids
  ! roundoff-sized momentum extrema flattening resolved scalar gradients.
  ! It requires the conserved layout rho,momentum(3),energy,...; callers still
  ! own EOS admission. consistent_velocity is a further opt-in;
  ! it couples advected momentum to the limited density slope and limits
  ! momentum relative to donor velocity. Applications select it explicitly.
  ! Preserve the caller's center on refusal, including nonfinite neighbors.
  subroutine interpolate_eb_linear_halo_3d(value,neighbors,offset,ok,separate_momentum,consistent_velocity)
    real(dp), intent(inout) :: value(:)
    real(dp), intent(in) :: neighbors(:,:),offset(3)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: separate_momentum,consistent_velocity
    real(dp) :: velocity,rm,rp
    real(dp) :: slopes(size(value),3),candidate(size(value)),dc,df,db,limited,scale
    integer :: a,c
    logical :: separate,consistent
    ok=.false.
    separate=.false.
    if(present(separate_momentum)) separate=separate_momentum
    consistent=.false.
    if(present(consistent_velocity)) consistent=consistent_velocity
    if(consistent.and..not.separate) return
    if(separate.and.size(value)<5) return
    if(any(shape(neighbors)/=[size(value),6])) return
    if(.not.all(ieee_is_finite(value)).or..not.all(ieee_is_finite(neighbors))) return
    if(.not.all(ieee_is_finite(offset)).or.any(abs(offset)>.5_dp)) return
    if(consistent) then
      if(value(1)<=0.or.any(neighbors(1,:)<=0)) return
    endif
    do a=1,3
      scale=1
      do c=1,size(value)
        dc=.5_dp*(neighbors(c,2*a)-neighbors(c,2*a-1))
        df=2*(neighbors(c,2*a)-value(c));db=2*(value(c)-neighbors(c,2*a-1))
        limited=0
        if(df*db>=0) limited=min(abs(df),abs(db))
        limited=sign(min(limited,abs(dc)),dc)
        if(separate.and.c>=2.and.c<=4) then
          slopes(c,a)=limited
        else
          if(dc/=0) scale=min(scale,limited/dc)
          slopes(c,a)=dc
        endif
      enddo
      if(separate) then
        slopes(1,a)=scale*slopes(1,a)
        slopes(5:,a)=scale*slopes(5:,a)
        if(consistent) then
          ! Limit momentum relative to the donor velocity, then add the
          ! density-coupled advective part. Constant velocity survives a
          ! scalar limiter clipping density; symmetric children conserve m.
          do c=2,4
            velocity=value(c)/value(1)
            rm=neighbors(c,2*a-1)-velocity*neighbors(1,2*a-1)
            rp=neighbors(c,2*a)-velocity*neighbors(1,2*a)
            dc=.5_dp*(rp-rm);df=2*rp;db=-2*rm
            if(.not.all(ieee_is_finite([velocity,rm,rp,dc,df,db]))) return
            limited=0
            if((df>=0.and.db>=0).or.(df<=0.and.db<=0)) limited=min(abs(df),abs(db))
            limited=sign(min(limited,abs(dc)),dc)
            slopes(c,a)=velocity*slopes(1,a)+limited
          enddo
        endif
      else
        slopes(:,a)=scale*slopes(:,a)
      endif
    enddo
    candidate=value+matmul(slopes,offset)
    if(.not.all(ieee_is_finite(candidate))) return
    value=candidate;ok=.true.
  end subroutine interpolate_eb_linear_halo_3d

  ! Overlay one same-resolution peer on previously filled face halos.
  ! peer_offset is the peer origin minus the target origin in cell indices.
  ! Only their shared open face is replaced; edges/corners and other donors
  ! remain untouched unless volume fractions request all peer-covered slots.
  ! The caller owns collection-wide overlap validation.
  subroutine overlay_reactive_eb_peer_halo_3d(species,target,peer,peer_offset,peer_state,peer_temperature, &
      ghost_state,ghost_temperature,ok,ghost_volume_fraction,frame,peer_internal_energy,ghost_internal_energy)
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(eb_geometry_3d), intent(in) :: target,peer
    integer, intent(in) :: peer_offset(3)
    real(dp), intent(in) :: peer_state(:,:,:,:),peer_temperature(:,:,:)
    real(dp), intent(inout) :: ghost_state(:,0:,0:,0:),ghost_temperature(0:,0:,0:)
    real(dp), optional, intent(inout) :: ghost_volume_fraction(0:,0:,0:)
    real(dp), optional, intent(in) :: peer_internal_energy(:,:,:)
    real(dp), optional, intent(inout) :: ghost_internal_energy(0:,0:,0:)
    real(dp), allocatable :: held(:,:,:),cell_held
    real(dp), allocatable :: volumes(:,:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: candidate(:,:,:,:),temps(:,:,:)
    real(dp) :: primitive(reactive_nprim(size(species))),recovered,speed,guess
    real(dp) :: spacing(3),origin(3),expected(3),tolerance,left_metrics(3),right_metrics(3)
    integer :: n(3),pn(3),nvar,axis,side,first(3),last(3),idx(3),source(3),face(3),peer_face(3),i,j,k
    logical :: valid
    ok=.false.
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    if (.not.target%is_valid().or..not.peer%is_valid()) return
    n=[target%nx,target%ny,target%nz];pn=[peer%nx,peer%ny,peer%nz]
    ! Only touching boxes are meaningful to this pair operation. Bound offsets
    ! before doing integer index arithmetic, including for rejected inputs.
    if (any(peer_offset < -pn).or.any(peer_offset>n)) return
    if (all(peer_offset<n.and.peer_offset+pn>0)) return
    spacing=[target%dx,target%dy,target%dz]
    if (any(abs([peer%dx,peer%dy,peer%dz]-spacing)>512._dp*epsilon(1._dp)*spacing)) return
    origin=[target%x_lower,target%y_lower,target%z_lower]
    expected=origin+real(peer_offset,dp)*spacing
    tolerance=512._dp*epsilon(1._dp)*max(1._dp,maxval(abs(origin)),maxval(abs(expected)))
    if (any(abs([peer%x_lower,peer%y_lower,peer%z_lower]-expected)>tolerance)) return
    nvar=reactive_nvar(size(species))
    if (.not.fields_match(peer_state,peer_temperature,peer,nvar)) return
    if (any(shape(ghost_state)/=[nvar,n+2]).or.any(shape(ghost_temperature)/=n+2)) return
    if(present(peer_internal_energy).neqv.present(ghost_internal_energy)) return
    if(present(peer_internal_energy)) then
      if(any(shape(peer_internal_energy)/=pn).or.any(shape(ghost_internal_energy)/=n+2)) return
      allocate(held(0:n(1)+1,0:n(2)+1,0:n(3)+1),cell_held)
      held=ghost_internal_energy
    endif
    if(present(ghost_volume_fraction)) then
      if(any(shape(ghost_volume_fraction)/=n+2)) return
      allocate(volumes(0:n(1)+1,0:n(2)+1,0:n(3)+1))
      volumes=ghost_volume_fraction
    end if
    allocate(candidate(nvar,0:n(1)+1,0:n(2)+1,0:n(3)+1),temps(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    candidate=ghost_state;temps=ghost_temperature
    do axis=1,3
      do side=0,1
        if (side==0.and.peer_offset(axis)+pn(axis)/=0) cycle
        if (side==1.and.peer_offset(axis)/=n(axis)) cycle
        first=max(1,peer_offset+1);last=min(n,peer_offset+pn)
        first(axis)=side*(n(axis)+1);last(axis)=first(axis)
        if (any(first>last)) cycle
        do k=first(3),last(3)
          do j=first(2),last(2)
            do i=first(1),last(1)
              idx=[i,j,k];source=idx-peer_offset
              face=idx;face(axis)=side*n(axis)
              peer_face=source;peer_face(axis)=source(axis)-side
              left_metrics=face_metrics(target,axis,face)
              right_metrics=face_metrics(peer,axis,peer_face)
              if (abs(left_metrics(1)-right_metrics(1))>512._dp*epsilon(1._dp)) return
              if (left_metrics(1)<=0._dp) cycle
              if (any(abs(left_metrics(2:3)-right_metrics(2:3))>512._dp*epsilon(1._dp))) return
              if (peer%volume_fraction(source(1),source(2),source(3))<=0._dp) return
              guess=peer_temperature(source(1),source(2),source(3))
              if (.not.ieee_is_finite(guess)) return
              if (guess<=0._dp) return
              if(present(peer_internal_energy)) cell_held=peer_internal_energy(source(1),source(2),source(3))
              call recover_reactive_eb_cell_3d(species,target,[i,j,k],peer_state(:,source(1),source(2),source(3)),guess, &
                primitive,recovered,speed,valid,frame,cell_held)
              if (.not.valid) return
              candidate(:,i,j,k)=peer_state(:,source(1),source(2),source(3));temps(i,j,k)=recovered
              if(present(peer_internal_energy)) held(i,j,k)=cell_held
            end do
          end do
        end do
      end do
    end do
    if(present(ghost_volume_fraction)) then
      first=max(0,peer_offset+1);last=min(n+1,peer_offset+pn)
      do k=first(3),last(3)
        do j=first(2),last(2)
          do i=first(1),last(1)
            source=[i,j,k]-peer_offset
            volumes(i,j,k)=peer%volume_fraction(source(1),source(2),source(3))
            candidate(:,i,j,k)=0;temps(i,j,k)=0
            if(present(peer_internal_energy)) held(i,j,k)=0
            if(volumes(i,j,k)<=0) cycle
            guess=peer_temperature(source(1),source(2),source(3))
            if(.not.ieee_is_finite(guess).or.guess<=0) return
            if(present(peer_internal_energy)) cell_held=peer_internal_energy(source(1),source(2),source(3))
            call recover_reactive_eb_cell_3d(species,target,[i,j,k],peer_state(:,source(1),source(2),source(3)),guess, &
              primitive,recovered,speed,valid,frame,cell_held)
            if(.not.valid) return
            candidate(:,i,j,k)=peer_state(:,source(1),source(2),source(3));temps(i,j,k)=recovered
            if(present(peer_internal_energy)) held(i,j,k)=cell_held
          end do
        end do
      end do
      ghost_volume_fraction=volumes
    end if
    if(present(ghost_internal_energy)) ghost_internal_energy=held
    ghost_state=candidate;ghost_temperature=temps;ok=.true.
  end subroutine overlay_reactive_eb_peer_halo_3d

  pure function face_metrics(geometry,axis,index) result(values)
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: axis,index(3)
    real(dp) :: values(3)
    integer :: i,j,k
    i=index(1);j=index(2);k=index(3)
    select case(axis)
    case(1)
      values=[geometry%x_face_fraction(i,j,k),geometry%x_face_centroid_y(i,j,k),geometry%x_face_centroid_z(i,j,k)]
    case(2)
      values=[geometry%y_face_fraction(i,j,k),geometry%y_face_centroid_x(i,j,k),geometry%y_face_centroid_z(i,j,k)]
    case(3)
      values=[geometry%z_face_fraction(i,j,k),geometry%z_face_centroid_x(i,j,k),geometry%z_face_centroid_y(i,j,k)]
    end select
  end function face_metrics

  ! PCM parent donors at internal patch faces, fine-state extrapolation at
  ! physical domain faces. Only open face slots are populated; all other halo
  ! storage is zero. Optional volume fractions request complete halos, with PCM
  ! parent geometry at coarse/fine boundaries. This is not full metric parity.
  subroutine fill_reactive_eb_patch_halo_3d(species,coarse,fine,lower,ratio,coarse_state,coarse_temperature, &
      fine_state,fine_temperature,ghost_state,ghost_temperature,ok,ghost_volume_fraction,frame)
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(eb_geometry_3d), intent(in) :: coarse,fine
    integer, intent(in) :: lower(3),ratio
    real(dp), intent(in) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:)
    real(dp), intent(in) :: fine_state(:,:,:,:),fine_temperature(:,:,:)
    real(dp), intent(inout) :: ghost_state(:,0:,0:,0:),ghost_temperature(0:,0:,0:)
    real(dp), optional, intent(inout) :: ghost_volume_fraction(0:,0:,0:)
    real(dp), allocatable :: volumes(:,:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: candidate(:,:,:,:),temps(:,:,:),primitive(:),donor(:)
    integer :: dims(3),parent_dims(3),upper(3),axis,side,i,j,k,start(3),finish(3),index(3),parent(3),inside(3),nvar
    real(dp) :: spacing(3),origin(3),fine_origin(3),fine_end(3),expected_end(3),tolerance,area,guess,recovered,speed
    logical :: valid
    ok=.false.
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    if (.not.coarse%is_valid().or..not.fine%is_valid().or.ratio<1) return
    dims=[fine%nx,fine%ny,fine%nz];parent_dims=[coarse%nx,coarse%ny,coarse%nz]
    if (any(mod(dims,ratio)/=0).or.any(lower<1).or.any(lower>parent_dims)) return
    if (any(dims/ratio>parent_dims-lower+1)) return
    upper=lower+dims/ratio-1
    spacing=[coarse%dx,coarse%dy,coarse%dz]
    origin=[coarse%x_lower,coarse%y_lower,coarse%z_lower]+real(lower-1,dp)*spacing
    expected_end=[coarse%x_lower,coarse%y_lower,coarse%z_lower]+real(upper,dp)*spacing
    fine_origin=[fine%x_lower,fine%y_lower,fine%z_lower];fine_end=[fine%x_upper,fine%y_upper,fine%z_upper]
    tolerance=512._dp*epsilon(1._dp)*max(1._dp,maxval(abs(origin)),maxval(abs(expected_end)))
    if (any(abs(fine_origin-origin)>tolerance).or.any(abs(fine_end-expected_end)>tolerance)) return
    if (any(abs([fine%dx,fine%dy,fine%dz]*ratio-spacing)> &
      512._dp*epsilon(1._dp)*spacing)) return
    nvar=reactive_nvar(size(species))
    if (.not.fields_match(coarse_state,coarse_temperature,coarse,nvar)) return
    if (.not.fields_match(fine_state,fine_temperature,fine,nvar)) return
    if (any(shape(ghost_state)/=[nvar,dims+2]).or.any(shape(ghost_temperature)/=dims+2)) return
    if(present(ghost_volume_fraction)) then
      if(any(shape(ghost_volume_fraction)/=dims+2)) return
      allocate(volumes(0:dims(1)+1,0:dims(2)+1,0:dims(3)+1))
      volumes=0
    end if
    allocate(candidate(nvar,0:dims(1)+1,0:dims(2)+1,0:dims(3)+1), &
      temps(0:dims(1)+1,0:dims(2)+1,0:dims(3)+1),primitive(reactive_nprim(size(species))),donor(nvar))
    candidate=0._dp;temps=0._dp
    do axis=1,3
      do side=0,1
        start=1;finish=dims;start(axis)=side*(dims(axis)+1);finish(axis)=start(axis)
        do k=start(3),finish(3)
          do j=start(2),finish(2)
            do i=start(1),finish(1)
              index=[i,j,k];inside=max(1,min(dims,index))
              select case(axis)
              case(1)
                area=fine%x_face_fraction(side*dims(1),j,k)
              case(2)
                area=fine%y_face_fraction(i,side*dims(2),k)
              case(3)
                area=fine%z_face_fraction(i,j,side*dims(3))
              end select
              if (area<=0._dp) cycle
              parent=lower+(inside-1)/ratio
              parent(axis)=lower(axis)-1
              if (side==1) parent(axis)=upper(axis)+1
              if (parent(axis)<1.or.parent(axis)>parent_dims(axis)) then
                donor=fine_state(:,inside(1),inside(2),inside(3))
                guess=fine_temperature(inside(1),inside(2),inside(3))
              else
                if (coarse%volume_fraction(parent(1),parent(2),parent(3))<=0._dp) return
                donor=coarse_state(:,parent(1),parent(2),parent(3))
                guess=coarse_temperature(parent(1),parent(2),parent(3))
              end if
              if (.not.ieee_is_finite(guess)) return
              if (guess<=0._dp) return
              call recover_reactive_eb_cell_3d(species,fine,[i,j,k],donor,guess,primitive,recovered,speed,valid,frame)
              if (.not.valid) return
              candidate(:,i,j,k)=donor;temps(i,j,k)=recovered
            end do
          end do
        end do
      end do
    end do
    if(present(ghost_volume_fraction)) then
      do k=0,dims(3)+1
        do j=0,dims(2)+1
          do i=0,dims(1)+1
            index=[i,j,k]
            ! Clamp physical-domain ghosts before selecting a fine or parent donor.
            where(lower==1) index=max(1,index)
            where(upper==parent_dims) index=min(dims,index)
            if(all(index>=1.and.index<=dims)) then
              volumes(i,j,k)=fine%volume_fraction(index(1),index(2),index(3))
              donor=fine_state(:,index(1),index(2),index(3))
              guess=fine_temperature(index(1),index(2),index(3))
            else
              parent=lower+(index+ratio-1)/ratio-1
              volumes(i,j,k)=coarse%volume_fraction(parent(1),parent(2),parent(3))
              donor=coarse_state(:,parent(1),parent(2),parent(3))
              guess=coarse_temperature(parent(1),parent(2),parent(3))
            end if
            if(allocated(fine%halo_volume_fraction)) then
              if(fine%halo_volume_fraction(index(1),index(2),index(3))>=0) &
                volumes(i,j,k)=fine%halo_volume_fraction(index(1),index(2),index(3))
            end if
            if(volumes(i,j,k)<=0) then
              candidate(:,i,j,k)=0;temps(i,j,k)=0
              cycle
            end if
            ! Open faces were already checked above. Retain their existing values.
            if(temps(i,j,k)>0) cycle
            if(.not.ieee_is_finite(guess).or.guess<=0) return
            call recover_reactive_eb_cell_3d(species,fine,[i,j,k],donor,guess,primitive,recovered,speed,valid,frame)
            if(.not.valid) return
            candidate(:,i,j,k)=donor;temps(i,j,k)=recovered
          end do
        end do
      end do
      ghost_volume_fraction=volumes
    end if
    ghost_state=candidate;ghost_temperature=temps;ok=.true.
  end subroutine fill_reactive_eb_patch_halo_3d

  ! Full-domain restriction onto the geometry implied by fine and ratio.
  ! Covered destination storage is retained. Publish U and T together only
  ! after all active donor and destination states pass the existing EOS checks.
  subroutine restrict_reactive_eb_state_3d(species,fine,ratio,fine_state,fine_temperature, &
      coarse_state,coarse_temperature,ok,frame,fine_internal_energy,coarse_internal_energy)
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(eb_geometry_3d), intent(in) :: fine
    integer, intent(in) :: ratio
    real(dp), intent(in) :: fine_state(:,:,:,:),fine_temperature(:,:,:)
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:)
    real(dp), optional, intent(in) :: fine_internal_energy(:,:,:)
    real(dp), optional, intent(inout) :: coarse_internal_energy(:,:,:)
    logical, intent(out) :: ok
    type(eb_geometry_3d) :: coarse
    real(dp), allocatable :: candidate(:,:,:,:),temps(:,:,:),primitive(:)
    real(dp), allocatable :: held(:,:,:)
    real(dp) :: weight,volume,recovered,sound_speed,guess
    integer :: i,j,k,ii,jj,kk,nvar,allocation_status,lo(3),hi(3)
    logical :: valid
    ok=.false.
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    call coarsen_eb_geometry_3d(fine,ratio,coarse,valid)
    if (.not.valid) return
    nvar=reactive_nvar(size(species))
    if (.not.fields_match(fine_state,fine_temperature,fine,nvar)) return
    if (.not.fields_match(coarse_state,coarse_temperature,coarse,nvar)) return
    allocate(candidate(nvar,coarse%nx,coarse%ny,coarse%nz), &
      temps(coarse%nx,coarse%ny,coarse%nz),primitive(reactive_nprim(size(species))),stat=allocation_status)
    if (allocation_status/=0) return
    if(present(fine_internal_energy).neqv.present(coarse_internal_energy)) return
    if(present(fine_internal_energy)) then
      if(any(shape(fine_internal_energy)/=shape(fine_temperature))) return
      if(any(shape(coarse_internal_energy)/=shape(coarse_temperature))) return
      allocate(held,source=coarse_internal_energy,stat=allocation_status)
      if(allocation_status/=0) return
    endif
    candidate=coarse_state; temps=coarse_temperature
    do k=1,coarse%nz
      do j=1,coarse%ny
        do i=1,coarse%nx
          if (coarse%volume_fraction(i,j,k)<=0._dp) cycle
          lo=([i,j,k]-1)*ratio+1; hi=[i,j,k]*ratio
          volume=sum(fine%volume_fraction(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
          if (volume<=0._dp) return
          candidate(:,i,j,k)=0._dp; guess=0._dp
          if(present(fine_internal_energy)) held(i,j,k)=0
          do kk=lo(3),hi(3)
            do jj=lo(2),hi(2)
              do ii=lo(1),hi(1)
                if (fine%volume_fraction(ii,jj,kk)<=0._dp) cycle
                if (.not.ieee_is_finite(fine_temperature(ii,jj,kk))) return
                if (fine_temperature(ii,jj,kk)<=0._dp) return
                if(present(fine_internal_energy)) then
                call recover_reactive_eb_cell_3d(species,fine,[ii,jj,kk],fine_state(:,ii,jj,kk),fine_temperature(ii,jj,kk), &
                  primitive,recovered,sound_speed,valid,frame, &
                  fine_internal_energy(ii,jj,kk))
                else
                call recover_reactive_eb_cell_3d(species,fine,[ii,jj,kk],fine_state(:,ii,jj,kk),fine_temperature(ii,jj,kk), &
                  primitive,recovered,sound_speed,valid,frame)
                endif
                if (.not.valid) return
                weight=fine%volume_fraction(ii,jj,kk)/volume
                candidate(:,i,j,k)=candidate(:,i,j,k)+weight*fine_state(:,ii,jj,kk)
                guess=guess+weight*recovered
                if(present(fine_internal_energy)) held(i,j,k)=held(i,j,k)+weight*fine_internal_energy(ii,jj,kk)
              end do
            end do
          end do
          if(present(fine_internal_energy)) then
          call recover_reactive_eb_cell_3d(species,coarse,[i,j,k],candidate(:,i,j,k),guess, &
            primitive,temps(i,j,k),sound_speed,valid,frame,held(i,j,k))
          else
          call recover_reactive_eb_cell_3d(species,coarse,[i,j,k],candidate(:,i,j,k),guess, &
            primitive,temps(i,j,k),sound_speed,valid,frame)
          endif
          if (.not.valid) return
        end do
      end do
    end do
    coarse_state=candidate; coarse_temperature=temps
    if(present(coarse_internal_energy)) coarse_internal_energy=held
    ok=.true.
  end subroutine restrict_reactive_eb_state_3d

  ! Piecewise constant prolongation preserves the kappa-weighted inventory
  ! because the coarse geometry is restricted from the supplied fine geometry.
  subroutine prolong_reactive_eb_state_pcm_3d(species,fine,ratio,coarse_state,coarse_temperature, &
      fine_state,fine_temperature,ok,frame,coarse_internal_energy,fine_internal_energy)
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(eb_geometry_3d), intent(in) :: fine
    integer, intent(in) :: ratio
    real(dp), intent(in) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:)
    real(dp), intent(inout) :: fine_state(:,:,:,:),fine_temperature(:,:,:)
    real(dp), optional, intent(in) :: coarse_internal_energy(:,:,:)
    real(dp), optional, intent(inout) :: fine_internal_energy(:,:,:)
    logical, intent(out) :: ok
    type(eb_geometry_3d) :: coarse
    real(dp), allocatable :: candidate(:,:,:,:),temps(:,:,:),primitive(:)
    real(dp), allocatable :: held(:,:,:)
    real(dp) :: recovered,sound_speed
    integer :: i,j,k,ii,jj,kk,nvar,allocation_status,lo(3),hi(3)
    logical :: valid
    ok=.false.
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    call coarsen_eb_geometry_3d(fine,ratio,coarse,valid)
    if (.not.valid) return
    nvar=reactive_nvar(size(species))
    if (.not.fields_match(fine_state,fine_temperature,fine,nvar)) return
    if (.not.fields_match(coarse_state,coarse_temperature,coarse,nvar)) return
    allocate(candidate(nvar,fine%nx,fine%ny,fine%nz),temps(fine%nx,fine%ny,fine%nz), &
      primitive(reactive_nprim(size(species))),stat=allocation_status)
    if (allocation_status/=0) return
    if(present(coarse_internal_energy).neqv.present(fine_internal_energy)) return
    if(present(coarse_internal_energy)) then
      if(any(shape(coarse_internal_energy)/=shape(coarse_temperature))) return
      if(any(shape(fine_internal_energy)/=shape(fine_temperature))) return
      allocate(held,source=fine_internal_energy,stat=allocation_status)
      if(allocation_status/=0) return
    endif
    candidate=fine_state; temps=fine_temperature
    do k=1,coarse%nz
      do j=1,coarse%ny
        do i=1,coarse%nx
          if (coarse%volume_fraction(i,j,k)<=0._dp) cycle
          if (.not.ieee_is_finite(coarse_temperature(i,j,k))) return
          if (coarse_temperature(i,j,k)<=0._dp) return
          if(present(coarse_internal_energy)) then
          call recover_reactive_eb_cell_3d(species,coarse,[i,j,k],coarse_state(:,i,j,k),coarse_temperature(i,j,k), &
            primitive,recovered,sound_speed,valid,frame,coarse_internal_energy(i,j,k))
          else
          call recover_reactive_eb_cell_3d(species,coarse,[i,j,k],coarse_state(:,i,j,k),coarse_temperature(i,j,k), &
            primitive,recovered,sound_speed,valid,frame)
          endif
          if (.not.valid) return
          lo=([i,j,k]-1)*ratio+1; hi=[i,j,k]*ratio
          do kk=lo(3),hi(3)
            do jj=lo(2),hi(2)
              do ii=lo(1),hi(1)
                if (fine%volume_fraction(ii,jj,kk)<=0._dp) cycle
                candidate(:,ii,jj,kk)=coarse_state(:,i,j,k)
                temps(ii,jj,kk)=recovered
                if(present(coarse_internal_energy)) held(ii,jj,kk)=coarse_internal_energy(i,j,k)
                if(present(frame)) then
                  if(present(coarse_internal_energy)) then
                  call recover_reactive_eb_cell_3d(species,fine,[ii,jj,kk],candidate(:,ii,jj,kk),recovered, &
                    primitive,temps(ii,jj,kk),sound_speed,valid,frame,held(ii,jj,kk))
                  else
                  call recover_reactive_eb_cell_3d(species,fine,[ii,jj,kk],candidate(:,ii,jj,kk),recovered, &
                    primitive,temps(ii,jj,kk),sound_speed,valid,frame)
                  endif
                  if(.not.valid) return
                end if
              end do
            end do
          end do
        end do
      end do
    end do
    fine_state=candidate; fine_temperature=temps
    if(present(fine_internal_energy)) fine_internal_energy=held
    ok=.true.
  end subroutine prolong_reactive_eb_state_pcm_3d

  ! Conserved transfer keeps rotating energy unchanged. Recover thermodynamics
  ! at the receiving cell center, including halo indices outside the patch.
  subroutine recover_reactive_eb_cell_3d(species,geometry,index,state,guess,primitive,temperature,speed,ok,frame, &
      internal_energy_density)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: index(3)
    real(dp), intent(in) :: state(:),guess
    real(dp), intent(out) :: primitive(:),temperature,speed
    logical, intent(out) :: ok
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: internal_energy_density
    real(dp), allocatable :: potential
    real(dp) :: position(3)
    primitive=0;temperature=0;speed=0;ok=.false.
    if(present(frame)) then
      allocate(potential)
      position=[geometry%x_lower,geometry%y_lower,geometry%z_lower]+ &
        (real(index,dp)-.5_dp)*[geometry%dx,geometry%dy,geometry%dz]
      call rotating_potential_3d(frame,position,potential,ok)
      if(.not.ok) return
    end if
    ! Preserve the independent UEINT all the way into EOS recovery. Building
    ! UEDEN from it and subtracting kinetic/potential energy loses the scalar
    ! again precisely where the dual-energy formulation is required.
    call reactive_conserved_to_primitive(species,state,guess,primitive,temperature,speed,ok, &
      rotating_potential=potential,internal_energy_density=internal_energy_density)
  end subroutine recover_reactive_eb_cell_3d

  pure logical function fields_match(state,temperature,geometry,nvar) result(valid)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:)
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: nvar
    valid=nvar>0.and.all(shape(state)==[nvar,geometry%nx,geometry%ny,geometry%nz]).and. &
      all(shape(temperature)==[geometry%nx,geometry%ny,geometry%nz])
  end function fields_match
end module eb_reactive_transfer_3d_mod
