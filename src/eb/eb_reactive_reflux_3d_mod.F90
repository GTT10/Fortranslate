module eb_reactive_reflux_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d
  implicit none
  private
  type, public :: eb_integrated_flux_3d
    ! Time-integrated flux per OPEN face area. Normal indices are 0:n;
    ! component and tangential indices start at 1. Apertures are applied once.
    real(dp), allocatable :: x(:,:,:,:),y(:,:,:,:),z(:,:,:,:)
  end type eb_integrated_flux_3d
  public :: restrict_eb_integrated_flux_3d,reflux_reactive_eb_state_3d
  public :: valid_eb_integrated_flux_3d
  public :: reconcile_reactive_eb_peer_flux_3d
contains
  ! The caller selects a deterministic owner for a shared same-level face.
  ! Replace the receiver's integrated face flux with the owner's physical
  ! transfer and correct the already advanced receiver state by the difference.
  ! Owner data is read-only; receiver flux, conserved state and T commit together.
  subroutine reconcile_reactive_eb_peer_flux_3d(species,owner,receiver,receiver_offset,owner_flux, &
      receiver_flux,state,temperature,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: owner,receiver
    integer, intent(in) :: receiver_offset(3)
    type(eb_integrated_flux_3d), intent(in) :: owner_flux
    type(eb_integrated_flux_3d), intent(inout) :: receiver_flux
    real(dp), optional, intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    logical, intent(out) :: ok
    type(eb_integrated_flux_3d) :: candidate_flux
    real(dp), allocatable :: candidate(:,:,:,:),temps(:,:,:)
    real(dp) :: spacing(3),rspacing(3),origin(3),expected(3),tol,oa,ra,owner_area,receiver_area,volume,sound_speed
    real(dp) :: primitive(reactive_nprim(size(species))),replacement(reactive_nvar(size(species)))
    integer :: n(3),rn(3),nvar,axis,side,first(3),last(3),index(3),cell(3),oface(3),rface(3),i,j,k
    logical :: valid
    ok=.false.
    if (.not.owner%is_valid().or..not.receiver%is_valid()) return
    n=[owner%nx,owner%ny,owner%nz];rn=[receiver%nx,receiver%ny,receiver%nz]
    if (any(receiver_offset < -rn).or.any(receiver_offset>n)) return
    if (all(receiver_offset<n.and.receiver_offset+rn>0)) return
    spacing=[owner%dx,owner%dy,owner%dz];rspacing=[receiver%dx,receiver%dy,receiver%dz]
    tol=512._dp*epsilon(1._dp)
    if (any(abs(rspacing-spacing)>tol*spacing)) return
    origin=[owner%x_lower,owner%y_lower,owner%z_lower]
    expected=origin+real(receiver_offset,dp)*spacing
    if (any(abs([receiver%x_lower,receiver%y_lower,receiver%z_lower]-expected)> &
      tol*max(1._dp,maxval(abs(origin)),maxval(abs(expected))))) return
    nvar=reactive_nvar(size(species))
    if (present(state).neqv.present(temperature)) return
    if (present(state)) then
      if (any(shape(state)/=[nvar,rn]).or.any(shape(temperature)/=rn)) return
    end if
    if (.not.valid_eb_integrated_flux_3d(owner_flux,owner,nvar).or. &
      .not.valid_eb_integrated_flux_3d(receiver_flux,receiver,nvar)) return
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
                call reactive_conserved_to_primitive(species,candidate(:,cell(1),cell(2),cell(3)), &
                  temperature(cell(1),cell(2),cell(3)),primitive,temps(cell(1),cell(2),cell(3)),sound_speed,valid)
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

  ! Apply fine-minus-coarse integrated flux only at interior interfaces between
  ! refined and unrefined coarse cells. Fine states are updated/averaged separately.
  ! With all_faces, apply the full flux difference on exposed cells, including
  ! outer faces. The caller then owns the synchronized outer-boundary budget.
  subroutine reflux_reactive_eb_state_3d(species,geometry,refined,coarse_flux,fine_flux,state,temperature,ok,all_faces)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: geometry
    logical, intent(in) :: refined(:,:,:)
    type(eb_integrated_flux_3d), intent(in) :: coarse_flux,fine_flux
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: all_faces
    real(dp), allocatable :: candidate(:,:,:,:),temps(:,:,:),primitive(:),correction(:)
    real(dp) :: spacing(3),volume,scale,sound_speed
    integer :: i,j,k,axis,step,idx(3),neighbor(3),face_idx(3),dims(3),nvar
    logical :: valid,touched,complete
    ok=.false.;complete=.false.
    if (present(all_faces)) complete=all_faces
    if (.not.geometry%is_valid()) return
    dims=[geometry%nx,geometry%ny,geometry%nz];nvar=reactive_nvar(size(species))
    if (any(shape(state)/=[nvar,dims]).or.any(shape(temperature)/=dims).or.any(shape(refined)/=dims)) return
    if (.not.valid_eb_integrated_flux_3d(coarse_flux,geometry,nvar).or. &
      .not.valid_eb_integrated_flux_3d(fine_flux,geometry,nvar)) return
    allocate(candidate,source=state);allocate(temps,source=temperature)
    allocate(primitive(reactive_nprim(size(species))),correction(nvar))
    spacing=[geometry%dx,geometry%dy,geometry%dz]
    do k=1,dims(3)
      do j=1,dims(2)
        do i=1,dims(1)
          if (refined(i,j,k).or.geometry%volume_fraction(i,j,k)<=0._dp) cycle
          idx=[i,j,k];correction=0._dp;touched=.false.
          volume=geometry%volume_fraction(i,j,k)*product(spacing)
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
              touched=.true.
            end do
          end do
          if (.not.touched) cycle
          if (.not.ieee_is_finite(temperature(i,j,k))) return
          if (temperature(i,j,k)<=0._dp) return
          candidate(:,i,j,k)=state(:,i,j,k)+correction
          call reactive_conserved_to_primitive(species,candidate(:,i,j,k),temperature(i,j,k), &
            primitive,temps(i,j,k),sound_speed,valid)
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
