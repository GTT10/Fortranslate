module eb_reactive_transfer_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d
  implicit none
  private
  public :: face_metrics
  public :: restrict_reactive_eb_state_3d,prolong_reactive_eb_state_pcm_3d
  public :: fill_reactive_eb_patch_halo_3d,overlay_reactive_eb_peer_halo_3d
contains
  ! Overlay one same-resolution peer on previously filled face halos.
  ! peer_offset is the peer origin minus the target origin in cell indices.
  ! Only their shared open face is replaced; edges/corners and other donors
  ! remain untouched. The caller owns collection-wide overlap validation.
  subroutine overlay_reactive_eb_peer_halo_3d(species,target,peer,peer_offset,peer_state,peer_temperature, &
      ghost_state,ghost_temperature,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: target,peer
    integer, intent(in) :: peer_offset(3)
    real(dp), intent(in) :: peer_state(:,:,:,:),peer_temperature(:,:,:)
    real(dp), intent(inout) :: ghost_state(:,0:,0:,0:),ghost_temperature(0:,0:,0:)
    logical, intent(out) :: ok
    real(dp), allocatable :: candidate(:,:,:,:),temps(:,:,:)
    real(dp) :: primitive(reactive_nprim(size(species))),recovered,speed,guess
    real(dp) :: spacing(3),origin(3),expected(3),tolerance,left_metrics(3),right_metrics(3)
    integer :: n(3),pn(3),nvar,axis,side,first(3),last(3),idx(3),source(3),face(3),peer_face(3),i,j,k
    logical :: valid
    ok=.false.
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
              call reactive_conserved_to_primitive(species,peer_state(:,source(1),source(2),source(3)),guess, &
                primitive,recovered,speed,valid)
              if (.not.valid) return
              candidate(:,i,j,k)=peer_state(:,source(1),source(2),source(3));temps(i,j,k)=recovered
            end do
          end do
        end do
      end do
    end do
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
  ! storage is zero. This checks Cartesian alignment, not full metric parity.
  subroutine fill_reactive_eb_patch_halo_3d(species,coarse,fine,lower,ratio,coarse_state,coarse_temperature, &
      fine_state,fine_temperature,ghost_state,ghost_temperature,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: coarse,fine
    integer, intent(in) :: lower(3),ratio
    real(dp), intent(in) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:)
    real(dp), intent(in) :: fine_state(:,:,:,:),fine_temperature(:,:,:)
    real(dp), intent(inout) :: ghost_state(:,0:,0:,0:),ghost_temperature(0:,0:,0:)
    logical, intent(out) :: ok
    real(dp), allocatable :: candidate(:,:,:,:),temps(:,:,:),primitive(:),donor(:)
    integer :: dims(3),parent_dims(3),upper(3),axis,side,i,j,k,start(3),finish(3),index(3),parent(3),inside(3),nvar
    real(dp) :: spacing(3),origin(3),fine_origin(3),fine_end(3),expected_end(3),tolerance,area,guess,recovered,speed
    logical :: valid
    ok=.false.
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
              call reactive_conserved_to_primitive(species,donor,guess,primitive,recovered,speed,valid)
              if (.not.valid) return
              candidate(:,i,j,k)=donor;temps(i,j,k)=recovered
            end do
          end do
        end do
      end do
    end do
    ghost_state=candidate;ghost_temperature=temps;ok=.true.
  end subroutine fill_reactive_eb_patch_halo_3d

  ! Full-domain restriction onto the geometry implied by fine and ratio.
  ! Covered destination storage is retained. Publish U and T together only
  ! after all active donor and destination states pass the existing EOS checks.
  subroutine restrict_reactive_eb_state_3d(species,fine,ratio,fine_state,fine_temperature, &
      coarse_state,coarse_temperature,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: fine
    integer, intent(in) :: ratio
    real(dp), intent(in) :: fine_state(:,:,:,:),fine_temperature(:,:,:)
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:)
    logical, intent(out) :: ok
    type(eb_geometry_3d) :: coarse
    real(dp), allocatable :: candidate(:,:,:,:),temps(:,:,:),primitive(:)
    real(dp) :: weight,volume,recovered,sound_speed,guess
    integer :: i,j,k,ii,jj,kk,nvar,allocation_status,lo(3),hi(3)
    logical :: valid
    ok=.false.
    call coarsen_eb_geometry_3d(fine,ratio,coarse,valid)
    if (.not.valid) return
    nvar=reactive_nvar(size(species))
    if (.not.fields_match(fine_state,fine_temperature,fine,nvar)) return
    if (.not.fields_match(coarse_state,coarse_temperature,coarse,nvar)) return
    allocate(candidate(nvar,coarse%nx,coarse%ny,coarse%nz), &
      temps(coarse%nx,coarse%ny,coarse%nz),primitive(reactive_nprim(size(species))),stat=allocation_status)
    if (allocation_status/=0) return
    candidate=coarse_state; temps=coarse_temperature
    do k=1,coarse%nz
      do j=1,coarse%ny
        do i=1,coarse%nx
          if (coarse%volume_fraction(i,j,k)<=0._dp) cycle
          lo=([i,j,k]-1)*ratio+1; hi=[i,j,k]*ratio
          volume=sum(fine%volume_fraction(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
          if (volume<=0._dp) return
          candidate(:,i,j,k)=0._dp; guess=0._dp
          do kk=lo(3),hi(3)
            do jj=lo(2),hi(2)
              do ii=lo(1),hi(1)
                if (fine%volume_fraction(ii,jj,kk)<=0._dp) cycle
                if (.not.ieee_is_finite(fine_temperature(ii,jj,kk))) return
                if (fine_temperature(ii,jj,kk)<=0._dp) return
                call reactive_conserved_to_primitive(species,fine_state(:,ii,jj,kk),fine_temperature(ii,jj,kk), &
                  primitive,recovered,sound_speed,valid)
                if (.not.valid) return
                weight=fine%volume_fraction(ii,jj,kk)/volume
                candidate(:,i,j,k)=candidate(:,i,j,k)+weight*fine_state(:,ii,jj,kk)
                guess=guess+weight*recovered
              end do
            end do
          end do
          call reactive_conserved_to_primitive(species,candidate(:,i,j,k),guess, &
            primitive,temps(i,j,k),sound_speed,valid)
          if (.not.valid) return
        end do
      end do
    end do
    coarse_state=candidate; coarse_temperature=temps
    ok=.true.
  end subroutine restrict_reactive_eb_state_3d

  ! Piecewise constant prolongation preserves the kappa-weighted inventory
  ! because the coarse geometry is restricted from the supplied fine geometry.
  subroutine prolong_reactive_eb_state_pcm_3d(species,fine,ratio,coarse_state,coarse_temperature, &
      fine_state,fine_temperature,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: fine
    integer, intent(in) :: ratio
    real(dp), intent(in) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:)
    real(dp), intent(inout) :: fine_state(:,:,:,:),fine_temperature(:,:,:)
    logical, intent(out) :: ok
    type(eb_geometry_3d) :: coarse
    real(dp), allocatable :: candidate(:,:,:,:),temps(:,:,:),primitive(:)
    real(dp) :: recovered,sound_speed
    integer :: i,j,k,ii,jj,kk,nvar,allocation_status,lo(3),hi(3)
    logical :: valid
    ok=.false.
    call coarsen_eb_geometry_3d(fine,ratio,coarse,valid)
    if (.not.valid) return
    nvar=reactive_nvar(size(species))
    if (.not.fields_match(fine_state,fine_temperature,fine,nvar)) return
    if (.not.fields_match(coarse_state,coarse_temperature,coarse,nvar)) return
    allocate(candidate(nvar,fine%nx,fine%ny,fine%nz),temps(fine%nx,fine%ny,fine%nz), &
      primitive(reactive_nprim(size(species))),stat=allocation_status)
    if (allocation_status/=0) return
    candidate=fine_state; temps=fine_temperature
    do k=1,coarse%nz
      do j=1,coarse%ny
        do i=1,coarse%nx
          if (coarse%volume_fraction(i,j,k)<=0._dp) cycle
          if (.not.ieee_is_finite(coarse_temperature(i,j,k))) return
          if (coarse_temperature(i,j,k)<=0._dp) return
          call reactive_conserved_to_primitive(species,coarse_state(:,i,j,k),coarse_temperature(i,j,k), &
            primitive,recovered,sound_speed,valid)
          if (.not.valid) return
          lo=([i,j,k]-1)*ratio+1; hi=[i,j,k]*ratio
          do kk=lo(3),hi(3)
            do jj=lo(2),hi(2)
              do ii=lo(1),hi(1)
                if (fine%volume_fraction(ii,jj,kk)<=0._dp) cycle
                candidate(:,ii,jj,kk)=coarse_state(:,i,j,k)
                temps(ii,jj,kk)=recovered
              end do
            end do
          end do
        end do
      end do
    end do
    fine_state=candidate; fine_temperature=temps
    ok=.true.
  end subroutine prolong_reactive_eb_state_pcm_3d

  pure logical function fields_match(state,temperature,geometry,nvar) result(valid)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:)
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: nvar
    valid=nvar>0.and.all(shape(state)==[nvar,geometry%nx,geometry%ny,geometry%nz]).and. &
      all(shape(temperature)==[geometry%nx,geometry%ny,geometry%nz])
  end function fields_match
end module eb_reactive_transfer_3d_mod
