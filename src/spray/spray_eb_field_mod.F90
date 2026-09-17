module spray_eb_field_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use spray_eb_intersection_mod, only: intersect_spray_level_set
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar,reactive_conserved_to_primitive
  implicit none
  private
  public :: eb_spray_stencil,sample_eb_spray_field,deposit_eb_spray_gas
contains
  ! Volume-weighted CIC, restricted to face-connected fluid cells in the local
  ! eight-cell support. The caller binds phi to the geometry built from it.
  subroutine eb_spray_stencil(g,phi,position,cells,weights,volumes,ok)
    type(eb_geometry_3d), intent(in) :: g
    real(dp), intent(in) :: phi(0:,0:,0:),position(3)
    integer, intent(inout) :: cells(3,8)
    real(dp), intent(inout) :: weights(8),volumes(8)
    logical, intent(out) :: ok
    integer :: n(3),base(3),owner(3),idx(3),c(3,8),i,j,k,q,r,pass,axis,d(3),face(3),dummy(3)
    real(dp) :: lo(3),hi(3),h(3),u(3),w(8),v(8),fraction,point(3),normal(3),aperture
    logical :: hit,valid,connected(8)
    ok=.false.
    if (.not.g%is_valid()) return
    n=[g%nx,g%ny,g%nz]
    if (any(shape(phi)/=n+1)) return
    lo=[g%x_lower,g%y_lower,g%z_lower];hi=[g%x_upper,g%y_upper,g%z_upper];h=[g%dx,g%dy,g%dz]
    hit=.false.;fraction=0;point=0;normal=0;dummy=0
    call intersect_spray_level_set(phi,lo,hi,position,position,hit,fraction,point,normal,dummy,valid)
    if (.not.valid) return
    u=(position-lo)/h+.5_dp;base=floor(u);owner=min(n,floor((position-lo)/h)+1)
    c=0;w=0;v=0;connected=.false.;q=0
    do k=0,1
    do j=0,1
    do i=0,1
      q=q+1;idx=base+[i,j,k]
      if (any(idx<1).or.any(idx>n)) cycle
      if (g%volume_fraction(idx(1),idx(2),idx(3))<=0) cycle
      c(:,q)=idx
      v(q)=g%volume_fraction(idx(1),idx(2),idx(3))*product(h)
      w(q)=product(max(0._dp,1-abs(u-real(idx,dp))))*g%volume_fraction(idx(1),idx(2),idx(3))
      connected(q)=all(idx==owner)
    end do
    end do
    end do
    if (.not.any(connected)) return
    do pass=1,7
      do q=1,8
        if (connected(q).or.w(q)<=0) cycle
        do r=1,8
          if (.not.connected(r)) cycle
          d=abs(c(:,q)-c(:,r))
          if (sum(d)/=1) cycle
          axis=maxloc(d,dim=1);face=min(c(:,q),c(:,r))
          select case(axis)
          case(1)
            aperture=g%x_face_fraction(face(1),face(2),face(3))
          case(2)
            aperture=g%y_face_fraction(face(1),face(2),face(3))
          case(3)
            aperture=g%z_face_fraction(face(1),face(2),face(3))
          end select
          if (aperture>0) connected(q)=.true.
        end do
      end do
    end do
    where(.not.connected) w=0
    if (sum(w)<=0) return
    w=w/sum(w)
    cells=c;weights=w;volumes=v;ok=.true.
  end subroutine

  subroutine sample_eb_spray_field(g,phi,field,positions,samples,ok)
    type(eb_geometry_3d), intent(in) :: g
    real(dp), intent(in) :: phi(0:,0:,0:),field(:,:,:,:),positions(:,:)
    real(dp), allocatable, intent(inout) :: samples(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: trial(:,:)
    real(dp) :: w(8),v(8)
    integer :: c(3,8),p,q,ios
    logical :: valid
    ok=.false.
    if (.not.g%is_valid().or.size(positions,1)/=3) return
    if (any(shape(phi)/=[g%nx+1,g%ny+1,g%nz+1])) return
    if (.not.all(ieee_is_finite(phi))) return
    if (any(shape(field)/=[size(field,1),g%nx,g%ny,g%nz]).or.size(field,1)<1) return
    if (.not.all(ieee_is_finite(field))) return
    if (maxval(abs(field))>1e100_dp) return
    allocate(trial(size(field,1),size(positions,2)),stat=ios)
    if (ios/=0) return
    trial=0;c=0;w=0;v=0
    do p=1,size(positions,2)
      call eb_spray_stencil(g,phi,positions(:,p),c,w,v,valid)
      if (.not.valid) return
      do q=1,8
        if (w(q)<=0) cycle
        trial(:,p)=trial(:,p)+w(q)*field(:,c(1,q),c(2,q),c(3,q))
      end do
    end do
    call move_alloc(trial,samples);ok=.true.
  end subroutine

  ! Contributions are integrated mass/momentum/energy/species, while state holds
  ! densities per fluid volume. Covered cells and their temperature are untouched.
  subroutine deposit_eb_spray_gas(species,g,phi,positions,contributions,state,temperature,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: g
    real(dp), intent(in) :: phi(0:,0:,0:),positions(:,:),contributions(:,:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: trial(:,:,:,:),temp(:,:,:),primitive(:),term(:)
    real(dp) :: w(8),v(8),sound
    integer :: c(3,8),nv,p,q,i,j,k,ios
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    if (.not.g%is_valid().or.size(positions,1)/=3) return
    if (any(shape(phi)/=[g%nx+1,g%ny+1,g%nz+1])) return
    if (.not.all(ieee_is_finite(phi))) return
    if (any(shape(state)/=[nv,g%nx,g%ny,g%nz])) return
    if (any(shape(temperature)/=[g%nx,g%ny,g%nz])) return
    if (any(shape(contributions)/=[nv,size(positions,2)])) return
    if (.not.all(ieee_is_finite(state)).or..not.all(ieee_is_finite(contributions))) return
    if (maxval(abs(state))>1e100_dp.or.maxval(abs(contributions))>1e100_dp) return
    allocate(trial(nv,g%nx,g%ny,g%nz),temp(g%nx,g%ny,g%nz),primitive(nv),term(nv),stat=ios)
    if (ios/=0) return
    trial=state;temp=temperature;c=0;w=0;v=0
    ! Reject invalid incoming active-cell states as well as invalid final states.
    call recover(valid)
    if (.not.valid) return
    do p=1,size(positions,2)
      call eb_spray_stencil(g,phi,positions(:,p),c,w,v,valid)
      if (.not.valid) return
      do q=1,8
        if (w(q)<=0) cycle
        term=w(q)*contributions(:,p)
        if (any(abs(term)>1e100_dp*min(1._dp,v(q)))) return
        term=term/v(q)
        i=c(1,q);j=c(2,q);k=c(3,q)
        trial(:,i,j,k)=trial(:,i,j,k)+term
        if (any(abs(trial(:,i,j,k))>1e100_dp)) return
      end do
    end do
    call recover(valid)
    if (.not.valid) return
    state=trial;temperature=temp;ok=.true.
  contains
    subroutine recover(accepted)
      logical, intent(out) :: accepted
      integer :: ii,jj,kk
      real(dp) :: t
      logical :: cell_ok
      accepted=.false.
      do kk=1,g%nz
      do jj=1,g%ny
      do ii=1,g%nx
        if (g%volume_fraction(ii,jj,kk)<=0) cycle
        call reactive_conserved_to_primitive(species,trial(:,ii,jj,kk),temperature(ii,jj,kk),primitive,t,sound,cell_ok)
        if (.not.cell_ok) return
        temp(ii,jj,kk)=t
      end do
      end do
      end do
      accepted=.true.
    end subroutine
  end subroutine
end module
