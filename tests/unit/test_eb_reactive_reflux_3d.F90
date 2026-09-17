program test_eb_reactive_reflux_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_is_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d, &
    extract_eb_geometry_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d,restrict_eb_integrated_flux_3d, &
    reflux_reactive_eb_state_3d,reconcile_reactive_eb_peer_flux_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(eb_geometry_3d) :: fine,coarse
  type(eb_integrated_flux_3d) :: ff,cf,restricted,saved_flux
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: state(:,:,:,:),temperature(:,:,:),saved(:,:,:,:),saved_t(:,:,:),reference(:),delta(:)
  real(dp) :: fine_change,coarse_change,volume,scale
  logical :: refined(3,3,3),ok
  integer :: ratio,n,nvar,i,j,k,c
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo load')
  nvar=reactive_nvar(size(species))
  allocate(state(nvar,3,3,3),temperature(3,3,3),saved(nvar,3,3,3),saved_t(3,3,3),reference(nvar),delta(nvar))
  allocate(cf%x(nvar,0:3,3,3),cf%y(nvar,3,0:3,3),cf%z(nvar,3,3,0:3))
  config%nx=3;config%ny=3;config%nz=3
  refined=.false.;refined(2,2,2)=.true.
  do ratio=2,3
    n=3*ratio
    call build_axis_plane_eb_geometry_3d(n,n,n,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'z',.37_dp,fine,ok)
    call require(ok,'fine cut geometry')
    call coarsen_eb_geometry_3d(fine,ratio,coarse,ok)
    call require(ok,'coarse cut geometry')
    call initialize_reactive_eb_density_sheet_3d(species,config,coarse,state,temperature,ok)
    call require(ok,'initial coarse state')
    saved=state;saved_t=temperature;reference=state(:,2,2,2)
    allocate(ff%x(nvar,0:n,n,n),ff%y(nvar,n,0:n,n),ff%z(nvar,n,n,0:n))
    do k=1,n
      do j=1,n
        do i=0,n
          ff%x(:,i,j,k)=1.e-6_dp*(i+2*j+3*k)*reference
        end do
      end do
    end do
    do k=1,n
      do j=0,n
        do i=1,n
          ff%y(:,i,j,k)=2.e-6_dp*(i+2*j+3*k)*reference
        end do
      end do
    end do
    do k=0,n
      do j=1,n
        do i=1,n
          ff%z(:,i,j,k)=3.e-6_dp*(i+2*j+3*k)*reference
        end do
      end do
    end do
    call restrict_eb_integrated_flux_3d(fine,ratio,ff,restricted,ok)
    call require(ok,'aperture weighted flux restriction')
    call require(all(restricted%z(:,:, :,0)==0._dp),'closed face flux zero')
    cf%x=0._dp;cf%y=0._dp;cf%z=0._dp
    call reflux_reactive_eb_state_3d(species,coarse,refined,cf,restricted,state,temperature,ok)
    call require(ok,'cut cell reflux')
    call require(maxval(abs(state-saved))>1.e-10_dp,'reflux is nontrivial')
    call require(all(state(:,2,2,2)==saved(:,2,2,2)).and.temperature(2,2,2)==saved_t(2,2,2), &
      'refined coarse storage unchanged')
    call require(all(state(:,1,1,3)==saved(:,1,1,3)).and.temperature(1,1,3)==saved_t(1,1,3), &
      'noninterface storage unchanged')
    call require(all(state(:,:,:,1)==saved(:,:,:,1)),'covered storage unchanged')
    ! Independent fine-volume divergence over the refined cube. Its integrated
    ! inventory change must cancel the unrefined coarse-cell corrections.
    delta=0._dp
    do k=ratio+1,2*ratio
      do j=ratio+1,2*ratio
        do i=ratio+1,2*ratio
          delta=delta-fine%dy*fine%dz*(fine%x_face_fraction(i,j,k)*ff%x(:,i,j,k)- &
            fine%x_face_fraction(i-1,j,k)*ff%x(:,i-1,j,k)) &
            -fine%dx*fine%dz*(fine%y_face_fraction(i,j,k)*ff%y(:,i,j,k)- &
            fine%y_face_fraction(i,j-1,k)*ff%y(:,i,j-1,k)) &
            -fine%dx*fine%dy*(fine%z_face_fraction(i,j,k)*ff%z(:,i,j,k)- &
            fine%z_face_fraction(i,j,k-1)*ff%z(:,i,j,k-1))
        end do
      end do
    end do
    volume=coarse%dx*coarse%dy*coarse%dz
    do c=1,nvar
      coarse_change=sum((state(c,:,:,:)-saved(c,:,:,:))*coarse%volume_fraction)*volume
      fine_change=delta(c)
      scale=max(1._dp,sum(abs(saved(c,:,:,:))*coarse%volume_fraction)*volume)
      call require(abs(coarse_change+fine_change)<=2.e-12_dp*scale,'composite flux conservation')
    end do
    state=saved;temperature=saved_t
    call reflux_reactive_eb_state_3d(species,coarse,refined,restricted,restricted,state,temperature,ok)
    call require(ok.and.all(state==saved),'matching flux no state correction')
    state=saved;temperature=saved_t
    saved_flux=restricted
    restricted%x=restricted%x*1.e12_dp
    restricted%y=restricted%y*1.e12_dp
    restricted%z=restricted%z*1.e12_dp
    call reflux_reactive_eb_state_3d(species,coarse,refined,cf,restricted,state,temperature,ok)
    call require(.not.ok.and.all(state==saved).and.all(temperature==saved_t),'nonphysical reflux rollback')
    restricted=saved_flux
    call restrict_eb_integrated_flux_3d(fine,0,ff,restricted,ok)
    call require(.not.ok.and.all(restricted%x==saved_flux%x).and.all(restricted%y==saved_flux%y).and. &
      all(restricted%z==saved_flux%z),'invalid restriction rollback')
    call reflux_reactive_eb_state_3d(species,coarse,refined(:2,:,:),cf,restricted,state,temperature,ok)
    call require(.not.ok.and.all(state==saved).and.all(temperature==saved_t),'bad mask rollback')
    temperature(1,2,2)=ieee_value(0._dp,ieee_quiet_nan)
    call reflux_reactive_eb_state_3d(species,coarse,refined,cf,restricted,state,temperature,ok)
    call require(.not.ok.and.all(state==saved).and.ieee_is_nan(temperature(1,2,2)), &
      'nonfinite interface temperature rejected without floating point trap')
    temperature=saved_t
    ff%x(1,0,1,1)=ieee_value(0._dp,ieee_quiet_nan)
    call restrict_eb_integrated_flux_3d(fine,ratio,ff,restricted,ok)
    call require(.not.ok.and.all(restricted%x==saved_flux%x),'nonfinite flux rollback')
    deallocate(ff%x,ff%y,ff%z)
  end do
  call check_peer_flux()
  write(*,'(a)') 'test_eb_reactive_reflux_3d: PASS'
contains
  subroutine check_peer_flux()
    type(eb_geometry_3d) :: domain,owner,receiver
    type(eb_integrated_flux_3d) :: owned,received,kept_flux
    real(dp), allocatable :: u(:,:,:,:),t(:,:,:),kept(:,:,:,:),kept_t(:,:,:)
    real(dp) :: vector(nvar),expected(nvar),actual(nvar),area,h(3),rvolume
    integer :: axis,side,ol(3),oh(3),rl(3),rh(3),offset(3),rn(3),first(3),last(3),component
    logical :: valid
    call build_axis_plane_eb_geometry_3d(6,6,6,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'z',.37_dp,domain,valid)
    call require(valid,'peer flux domain geometry')
    allocate(owned%x(nvar,0:3,3,3),owned%y(nvar,3,0:3,3),owned%z(nvar,3,3,0:3))
    do axis=1,3
      do side=0,1
        ol=1;oh=3;rl=2;rh=3
        if (side==1) then
          rl(axis)=4;rh(axis)=6
        else
          ol(axis)=4;oh(axis)=6;rl(axis)=1;rh(axis)=3
        end if
        offset=rl-ol;rn=rh-rl+1
        call extract_eb_geometry_3d(domain,ol,oh,owner,valid)
        call require(valid,'owner geometry')
        call extract_eb_geometry_3d(domain,rl,rh,receiver,valid)
        call require(valid,'receiver partial face geometry')
        allocate(u(nvar,rn(1),rn(2),rn(3)),t(rn(1),rn(2),rn(3)))
        allocate(received%x(nvar,0:rn(1),rn(2),rn(3)),received%y(nvar,rn(1),0:rn(2),rn(3)), &
          received%z(nvar,rn(1),rn(2),0:rn(3)))
        config%nx=rn(1);config%ny=rn(2);config%nz=rn(3)
        call initialize_reactive_eb_density_sheet_3d(species,config,receiver,u,t,valid)
        call require(valid,'receiver physical state')
        vector=1.e-7_dp*u(:,rn(1),rn(2),rn(3))
        do component=1,nvar
          owned%x(component,:,:,:)=vector(component)
          owned%y(component,:,:,:)=vector(component)
          owned%z(component,:,:,:)=vector(component)
        end do
        received%x=0._dp;received%y=0._dp;received%z=0._dp
        kept=u;kept_t=t
        call reconcile_reactive_eb_peer_flux_3d(species,owner,receiver,offset,owned,received,u,t,valid)
        call require(valid,'peer flux synchronization')
        first=offset+1;last=offset+rn;first(axis)=side*3;last(axis)=first(axis)
        select case(axis)
        case(1);area=sum(owner%x_face_fraction(first(1):last(1),first(2):last(2),first(3):last(3)))
        case(2);area=sum(owner%y_face_fraction(first(1):last(1),first(2):last(2),first(3):last(3)))
        case(3);area=sum(owner%z_face_fraction(first(1):last(1),first(2):last(2),first(3):last(3)))
        end select
        h=[owner%dx,owner%dy,owner%dz];area=area*product(h)/h(axis)
        expected=real(2*side-1,dp)*area*vector
        rvolume=receiver%dx*receiver%dy*receiver%dz
        do component=1,nvar
          actual(component)=rvolume*sum((u(component,:,:,:)-kept(component,:,:,:))*receiver%volume_fraction)
        end do
        call require(maxval(abs(expected))>1.e-10_dp,'peer correction is nonzero')
        call require(all(abs(actual-expected)<3.e-12_dp*max(1._dp,abs(expected))),'receiver inventory equals face difference')
        call require(all(t>0._dp),'corrected peer temperatures positive')
        kept=u;kept_t=t;kept_flux=received
        ! Matching fluxes are idempotent (up to EOS temperature recovery).
        call reconcile_reactive_eb_peer_flux_3d(species,owner,receiver,offset,owned,received,u,t,valid)
        call require(valid.and.all(u==kept),'synchronized conserved state is idempotent')
        kept_t=t
        do component=1,nvar
          owned%x(component,:,:,:)=-real(2*side-1,dp)*1.e13_dp*vector(component)
          owned%y(component,:,:,:)=-real(2*side-1,dp)*1.e13_dp*vector(component)
          owned%z(component,:,:,:)=-real(2*side-1,dp)*1.e13_dp*vector(component)
        end do
        call reconcile_reactive_eb_peer_flux_3d(species,owner,receiver,offset,owned,received,u,t,valid)
        call require(.not.valid.and.all(u==kept).and.all(t==kept_t),'nonphysical correction rolls back state')
        call require(all(received%x==kept_flux%x).and.all(received%y==kept_flux%y).and. &
          all(received%z==kept_flux%z),'failed correction rolls back fluxes')
        call reconcile_reactive_eb_peer_flux_3d(species,owner,receiver,[0,0,0],owned,received,u,t,valid)
        call require(.not.valid.and.all(u==kept).and.all(t==kept_t),'overlapping peer flux boxes rejected')
        deallocate(u,t,received%x,received%y,received%z)
      end do
    end do
  end subroutine check_peer_flux

  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not.condition) error stop message
  end subroutine require
end program test_eb_reactive_reflux_3d
