module eb_reactive_hydro_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use simulation_config_reactive_3d_mod, only: map_reactive_boundary_donor_3d,valid_reactive_boundaries_3d
  use reactive_3d_mod, only: adjust_reactive_artificial_viscosity_3d
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive, &
    reactive_riemann_flux_x,reconstruct_pelec_plm_cell,compute_pelec_plm_slopes,reactive_pelec_primitive_flux_x, &
    reactive_ppm_flattening_coefficient,reactive_pelec_source_to_primitive
  use reactive_2d_mod, only: reactive_riemann_flux_y
  use reactive_directional_flux_3d_mod, only: reactive_riemann_flux_z
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d,eb_regular_cell_3d,extract_eb_geometry_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d,valid_eb_integrated_flux_3d
  use eb_reactive_wall_flux_3d_mod, only: &
    reactive_eb_flux_divergence_3d,reactive_eb_slip_wall_source_3d,interpolate_reactive_eb_centroid_fluxes_3d
  use eb_reactive_redistribution_3d_mod, only: &
    advance_reactive_eb_redistributed_3d, &
    advance_reactive_eb_state_redistributed_3d
  implicit none
  private

  type, public :: eb_plm_boundary_3d
    ! Physical cell bounds in working-geometry coordinates, not cropped indices.
    integer :: lower(3)=1,upper(3)=1
    character(len=24) :: faces(6)='interior'
  end type eb_plm_boundary_3d

  type, public :: eb_plm_face_prediction_3d
    ! Normal face index starts at zero, tangential indices at one.
    real(dp), allocatable :: left(:,:,:,:),right(:,:,:,:),flux(:,:,:,:)
    real(dp), allocatable :: energy_left(:,:,:),energy_right(:,:,:),pressure(:,:,:),velocity(:,:,:)
    ! Independent intensive UEINT Riemann flux, captured only on request.
    real(dp), allocatable :: internal_energy_flux(:,:,:)
  end type eb_plm_face_prediction_3d
  public :: compute_reactive_eb_plm_normal_fluxes_3d,apply_reactive_eb_transverse_3d
  public :: compute_reactive_eb_plm_ctu_fluxes_3d,advance_reactive_eb_plm_ctu_3d
  public :: compute_eb_godunov_pressure_work_3d
  public :: compute_sdc_eb_hydro_3d,compute_sdc_eb_boundary_hydro_3d
  public :: advance_reactive_eb_boundary_ctu_3d
  public :: advance_reactive_eb_outflow_ctu_3d
  public :: compute_reactive_eb_flattening_3d
  public :: limit_reactive_eb_flux_correction_3d,correct_reactive_eb_spatial_hydro_flux_3d
  public :: predict_reactive_eb_plm_cell_3d
  public :: reactive_eb_outflow_riemann_fluxes_3d,reactive_eb_riemann_fluxes_3d
  public :: advance_reactive_eb_euler_3d
  public :: advance_reactive_eb_redistributed_euler_3d
  public :: advance_reactive_eb_state_redistributed_euler_3d

  interface reactive_eb_riemann_fluxes_3d
    module procedure reactive_eb_outflow_riemann_fluxes_3d
  end interface

contains

  ! Correct the whole spatial CTU region before its valid crop and ML
  ! redistribution. The caller supplies actual grown state/source/geometry;
  ! the retained hydro residual includes the EB wall pressure contribution.
  subroutine correct_reactive_eb_spatial_hydro_flux_3d(geometry,state,source,dt,hydro_rate,flux,baseline,ok, &
      internal_energy_rate,internal_energy_flux,internal_energy_baseline)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: state(:,:,:,:),source(:,:,:,:),dt
    real(dp), intent(inout) :: hydro_rate(:,:,:,:)
    type(eb_integrated_flux_3d), intent(inout) :: flux
    type(eb_integrated_flux_3d), intent(in) :: baseline
    logical, intent(out) :: ok
    real(dp), optional, intent(inout) :: internal_energy_rate(:,:,:)
    type(eb_integrated_flux_3d), optional, intent(inout) :: internal_energy_flux
    type(eb_integrated_flux_3d), optional, intent(in) :: internal_energy_baseline
    type(eb_integrated_flux_3d) :: held_original,held_low,held_limited,held_delta
    real(dp), allocatable :: held_candidate(:,:,:)
    real(dp) :: held_increment
    integer :: held_options
    logical :: held
    type(eb_integrated_flux_3d) :: original,low,limited,delta
    real(dp), allocatable :: candidate(:,:,:,:)
    real(dp) :: areas(3),volume,increment(size(state,1)),endpoint(size(state,1))
    integer :: n(3),nv,i,j,k,c
    logical :: good
    ok=.false.;n=[geometry%nx,geometry%ny,geometry%nz];nv=size(state,1)
    held_options=count([present(internal_energy_rate),present(internal_energy_flux),present(internal_energy_baseline)])
    if(held_options/=0.and.held_options/=3) return
    held=held_options==3
    if(.not.geometry%is_valid().or.nv<6.or..not.ieee_is_finite(dt).or.dt<=0) return
    if(any(shape(state)/=[nv,n]).or.any(shape(source)/=[nv,n]).or.any(shape(hydro_rate)/=[nv,n])) return
    if(held) then
      if(any(shape(internal_energy_rate)/=n)) return
      if(.not.all(ieee_is_finite(internal_energy_rate))) return
      if(.not.valid_eb_integrated_flux_3d(internal_energy_flux,geometry,1).or. &
        .not.valid_eb_integrated_flux_3d(internal_energy_baseline,geometry,1)) return
      if(dt<1) then
        if(maxval(abs(internal_energy_flux%x))>huge(1._dp)*dt.or. &
          maxval(abs(internal_energy_flux%y))>huge(1._dp)*dt.or. &
          maxval(abs(internal_energy_flux%z))>huge(1._dp)*dt) return
        if(maxval(abs(internal_energy_baseline%x))>huge(1._dp)*dt.or. &
          maxval(abs(internal_energy_baseline%y))>huge(1._dp)*dt.or. &
          maxval(abs(internal_energy_baseline%z))>huge(1._dp)*dt) return
      endif
    endif
    if(.not.allocated(flux%x).or..not.allocated(flux%y).or..not.allocated(flux%z)) return
    if(.not.allocated(baseline%x).or..not.allocated(baseline%y).or..not.allocated(baseline%z)) return
    if(any(shape(flux%x)/=[nv,n(1)+1,n(2),n(3)]).or. &
      any(shape(flux%y)/=[nv,n(1),n(2)+1,n(3)]).or.any(shape(flux%z)/=[nv,n(1),n(2),n(3)+1])) return
    if(any(shape(baseline%x)/=shape(flux%x)).or.any(shape(baseline%y)/=shape(flux%y)).or. &
      any(shape(baseline%z)/=shape(flux%z))) return
    if(.not.all(ieee_is_finite(state)).or..not.all(ieee_is_finite(source)).or. &
      .not.all(ieee_is_finite(hydro_rate))) return
    if(.not.all(ieee_is_finite(flux%x)).or..not.all(ieee_is_finite(flux%y)).or. &
      .not.all(ieee_is_finite(flux%z))) return
    if(.not.all(ieee_is_finite(baseline%x)).or..not.all(ieee_is_finite(baseline%y)).or. &
      .not.all(ieee_is_finite(baseline%z))) return
    if(dt<1) then
      if(maxval(abs(flux%x))>huge(1._dp)*dt.or.maxval(abs(flux%y))>huge(1._dp)*dt.or. &
        maxval(abs(flux%z))>huge(1._dp)*dt) return
      if(maxval(abs(baseline%x))>huge(1._dp)*dt.or.maxval(abs(baseline%y))>huge(1._dp)*dt.or. &
        maxval(abs(baseline%z))>huge(1._dp)*dt) return
    endif
    original=flux;low=baseline
    original%x=original%x/dt;original%y=original%y/dt;original%z=original%z/dt
    low%x=low%x/dt;low%y=low%y/dt;low%z=low%z/dt;limited=original
    if(held) then
      held_original=internal_energy_flux;held_low=internal_energy_baseline
      held_original%x=held_original%x/dt;held_original%y=held_original%y/dt;held_original%z=held_original%z/dt
      held_low%x=held_low%x/dt;held_low%y=held_low%y/dt;held_low%z=held_low%z/dt;held_limited=held_original
      call limit_reactive_eb_flux_correction_3d(geometry,state,low,limited,dt,good,source=source, &
        internal_energy_baseline=held_low,internal_energy_flux=held_limited)
    else
      call limit_reactive_eb_flux_correction_3d(geometry,state,low,limited,dt,good,source=source)
    endif
    if(.not.good) return
    where(limited%x==original%x)
      limited%x=flux%x
    elsewhere
      limited%x=dt*limited%x
    endwhere
    where(limited%y==original%y)
      limited%y=flux%y
    elsewhere
      limited%y=dt*limited%y
    endwhere
    where(limited%z==original%z)
      limited%z=flux%z
    elsewhere
      limited%z=dt*limited%z
    endwhere
    delta=limited;delta%x=delta%x-flux%x;delta%y=delta%y-flux%y;delta%z=delta%z-flux%z
    if(held) then
      where(held_limited%x==held_original%x)
        held_limited%x=internal_energy_flux%x
      elsewhere
        held_limited%x=dt*held_limited%x
      endwhere
      where(held_limited%y==held_original%y)
        held_limited%y=internal_energy_flux%y
      elsewhere
        held_limited%y=dt*held_limited%y
      endwhere
      where(held_limited%z==held_original%z)
        held_limited%z=internal_energy_flux%z
      elsewhere
        held_limited%z=dt*held_limited%z
      endwhere
      held_delta=held_limited
      held_delta%x=held_delta%x-internal_energy_flux%x
      held_delta%y=held_delta%y-internal_energy_flux%y
      held_delta%z=held_delta%z-internal_energy_flux%z
      held_candidate=internal_energy_rate
    endif
    candidate=hydro_rate;areas=[geometry%dy*geometry%dz,geometry%dx*geometry%dz,geometry%dx*geometry%dy]
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      if(geometry%volume_fraction(i,j,k)<=0) cycle
      volume=dt*geometry%dx*geometry%dy*geometry%dz*geometry%volume_fraction(i,j,k)
      if(.not.ieee_is_finite(volume).or.volume<=0) return
      increment=( &
        areas(1)*(geometry%x_face_fraction(i-1,j,k)*delta%x(:,i-1,j,k)- &
          geometry%x_face_fraction(i,j,k)*delta%x(:,i,j,k))+ &
        areas(2)*(geometry%y_face_fraction(i,j-1,k)*delta%y(:,i,j-1,k)- &
          geometry%y_face_fraction(i,j,k)*delta%y(:,i,j,k))+ &
        areas(3)*(geometry%z_face_fraction(i,j,k-1)*delta%z(:,i,j,k-1)- &
          geometry%z_face_fraction(i,j,k)*delta%z(:,i,j,k)))/volume
      if(any(increment/=0)) candidate(:,i,j,k)=hydro_rate(:,i,j,k)+increment
      if(held) then
        held_increment=( &
          areas(1)*(geometry%x_face_fraction(i-1,j,k)*held_delta%x(1,i-1,j,k)- &
            geometry%x_face_fraction(i,j,k)*held_delta%x(1,i,j,k))+ &
          areas(2)*(geometry%y_face_fraction(i,j-1,k)*held_delta%y(1,i,j-1,k)- &
            geometry%y_face_fraction(i,j,k)*held_delta%y(1,i,j,k))+ &
          areas(3)*(geometry%z_face_fraction(i,j,k-1)*held_delta%z(1,i,j,k-1)- &
            geometry%z_face_fraction(i,j,k)*held_delta%z(1,i,j,k)))/volume
        ! Only the advective increment changes. Preserve the captured pdivu
        ! residual, including exact identity when this face family is unchanged.
        if(held_increment/=0) held_candidate(i,j,k)=internal_energy_rate(i,j,k)+held_increment
      endif
      do c=1,nv
        if(c>=2.and.c<=5) cycle
        if(increment(c)==0) cycle
        candidate(c,i,j,k)=( &
          areas(1)*(geometry%x_face_fraction(i-1,j,k)*limited%x(c,i-1,j,k)- &
            geometry%x_face_fraction(i,j,k)*limited%x(c,i,j,k))+ &
          areas(2)*(geometry%y_face_fraction(i,j-1,k)*limited%y(c,i,j-1,k)- &
            geometry%y_face_fraction(i,j,k)*limited%y(c,i,j,k))+ &
          areas(3)*(geometry%z_face_fraction(i,j,k-1)*limited%z(c,i,j,k-1)- &
            geometry%z_face_fraction(i,j,k)*limited%z(c,i,j,k)))/volume
      enddo
      endpoint=state(:,i,j,k)+dt*(candidate(:,i,j,k)+source(:,i,j,k))
      if(any(.not.ieee_is_finite(endpoint)).or.endpoint(1)<=0.or.any(endpoint(6:)<0)) return
    enddo
    enddo
    enddo
    if(held) then
      if(.not.all(ieee_is_finite(held_candidate)).or..not.all(ieee_is_finite(held_limited%x)).or. &
        .not.all(ieee_is_finite(held_limited%y)).or..not.all(ieee_is_finite(held_limited%z))) return
    endif
    hydro_rate=candidate;flux=limited
    if(held) then
      internal_energy_rate=held_candidate;internal_energy_flux=held_limited
    endif
    ok=.true.
  end subroutine correct_reactive_eb_spatial_hydro_flux_3d

  ! Baseline and corrected fluxes are intensive and have identical physical
  ! face bounds. One coefficient per face scales every conserved component.
  ! The baseline forward state, including the supplied source, must already
  ! have nonnegative density and species. This routine does not repair it.
  ! peer_alpha supplies the opposite boundary donors in one halo layer,
  ! indexed 0:N+1 on each axis. Obtain both patches' allowances from the
  ! original baseline/corrected pair, exchange them, and reapply to that
  ! original pair. Shared fluxes and apertures must agree before this call.
  ! Without peer_alpha boundary faces only constrain the local donor.
  ! alpha_output component 1 is density; components 2: are species in order.
  subroutine limit_reactive_eb_flux_correction_3d(geometry,state,baseline,flux,dt,ok,source,alpha_output,peer_alpha, &
      internal_energy_baseline,internal_energy_flux)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: state(:,:,:,:),dt
    type(eb_integrated_flux_3d), intent(in) :: baseline
    type(eb_integrated_flux_3d), intent(inout) :: flux
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: source(:,:,:,:)
    real(dp), allocatable, optional, intent(inout) :: alpha_output(:,:,:,:)
    real(dp), optional, intent(in) :: peer_alpha(:,0:,0:,0:)
    type(eb_integrated_flux_3d), optional, intent(in) :: internal_energy_baseline
    type(eb_integrated_flux_3d), optional, intent(inout) :: internal_energy_flux
    type(eb_integrated_flux_3d) :: candidate,held_candidate
    real(dp), allocatable :: alpha(:,:,:,:),available(:,:,:,:)
    real(dp) :: volume,cell_volume,areas(3),budget,outgoing,delta,theta,face_area
    integer :: n(3),nv,nc,i,j,k,c,component,axis,face(3),left(3),right(3),lo(3),hi(3)
    ok=.false.;n=[geometry%nx,geometry%ny,geometry%nz];nv=size(state,1);nc=nv-4
    if(nv<6.or..not.geometry%is_valid().or.any(shape(state)/=[nv,n])) return
    if(.not.ieee_is_finite(dt).or.dt<=0) return
    if(present(source)) then
      if(any(shape(source)/=[nv,n])) return
    endif
    if(present(peer_alpha)) then
      if(any(shape(peer_alpha)/=[nc,n+2])) return
      if(any(.not.ieee_is_finite(peer_alpha))) return
      if(any(peer_alpha<0).or.any(peer_alpha>1)) return
    endif
    if(.not.valid_flux(baseline).or..not.valid_flux(flux)) return
    volume=geometry%dx*geometry%dy*geometry%dz
    areas=[geometry%dy*geometry%dz,geometry%dx*geometry%dz,geometry%dx*geometry%dy]
    if(.not.ieee_is_finite(volume).or.volume<=0.or.any(.not.ieee_is_finite(areas))) return
    if(present(internal_energy_baseline).neqv.present(internal_energy_flux)) return
    if(present(internal_energy_flux)) then
      if(.not.valid_eb_integrated_flux_3d(internal_energy_baseline,geometry,1).or. &
        .not.valid_eb_integrated_flux_3d(internal_energy_flux,geometry,1)) return
      held_candidate=internal_energy_flux
    endif
    candidate=flux
    allocate(alpha(nc,n(1),n(2),n(3)),available(nc,n(1),n(2),n(3)))
    alpha=1;available=0
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      if(geometry%volume_fraction(i,j,k)<=0) cycle
      if(any(.not.ieee_is_finite(state(:,i,j,k)))) return
      if(present(source)) then
        if(any(.not.ieee_is_finite(source(:,i,j,k)))) return
      endif
      cell_volume=volume*geometry%volume_fraction(i,j,k)
      do c=1,nc
        component=c+4
        if(c==1) component=1
        budget=cell_volume*state(component,i,j,k)
        if(present(source)) budget=budget+dt*cell_volume*source(component,i,j,k)
        outgoing=0
        do axis=1,3
          face=[i,j,k];face(axis)=face(axis)-1
          face_area=areas(axis)*aperture(axis,face)
          budget=budget+dt*face_area*value(baseline,axis,component,face)
          delta=value(flux,axis,component,face)-value(baseline,axis,component,face)
          outgoing=outgoing+dt*face_area*max(0._dp,-delta)
          face(axis)=face(axis)+1
          face_area=areas(axis)*aperture(axis,face)
          budget=budget-dt*face_area*value(baseline,axis,component,face)
          delta=value(flux,axis,component,face)-value(baseline,axis,component,face)
          outgoing=outgoing+dt*face_area*max(0._dp,delta)
        enddo
        if(.not.ieee_is_finite(budget).or..not.ieee_is_finite(outgoing)) return
        if(budget<0) return
        available(c,i,j,k)=budget
        ! A smaller outgoing face can disappear when added to a saturated
        ! budget. Apply the existing margin at equality too, so cancelling
        ! the baseline inflow cannot leave that represented face unbalanced.
        if(outgoing>=budget.and.outgoing>0) alpha(c,i,j,k)=(budget/outgoing)*(1-16*epsilon(1._dp))
      enddo
    enddo
    enddo
    enddo
    do axis=1,3
      lo=1;lo(axis)=0;hi=n
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        face=[i,j,k]
        if(aperture(axis,face)<=0) cycle
        left=face;right=face;right(axis)=right(axis)+1
        theta=1
        do c=1,nc
          component=c+4
          if(c==1) component=1
          delta=value(flux,axis,component,face)-value(baseline,axis,component,face)
          if(delta>0.and.left(axis)>0) theta=min(theta,alpha(c,left(1),left(2),left(3)))
          if(delta<0.and.right(axis)<=n(axis)) theta=min(theta,alpha(c,right(1),right(2),right(3)))
          if(present(peer_alpha)) then
            if(delta>0.and.left(axis)==0) theta=min(theta,peer_alpha(c,left(1),left(2),left(3)))
            if(delta<0.and.right(axis)==n(axis)+1) theta=min(theta,peer_alpha(c,right(1),right(2),right(3)))
          endif
        enddo
        if(theta>=1) cycle
        select case(axis)
        case(1);candidate%x(:,i,j,k)=baseline%x(:,i,j,k)+theta*(flux%x(:,i,j,k)-baseline%x(:,i,j,k))
        case(2);candidate%y(:,i,j,k)=baseline%y(:,i,j,k)+theta*(flux%y(:,i,j,k)-baseline%y(:,i,j,k))
        case(3);candidate%z(:,i,j,k)=baseline%z(:,i,j,k)+theta*(flux%z(:,i,j,k)-baseline%z(:,i,j,k))
        end select
        ! The extra held face uses exactly the gas face allowance. No
        ! internal-energy clipping or second acceptance coefficient is used.
        if(present(internal_energy_flux)) then
          select case(axis)
          case(1);held_candidate%x(:,i,j,k)=internal_energy_baseline%x(:,i,j,k)+ &
            theta*(internal_energy_flux%x(:,i,j,k)-internal_energy_baseline%x(:,i,j,k))
          case(2);held_candidate%y(:,i,j,k)=internal_energy_baseline%y(:,i,j,k)+ &
            theta*(internal_energy_flux%y(:,i,j,k)-internal_energy_baseline%y(:,i,j,k))
          case(3);held_candidate%z(:,i,j,k)=internal_energy_baseline%z(:,i,j,k)+ &
            theta*(internal_energy_flux%z(:,i,j,k)-internal_energy_baseline%z(:,i,j,k))
          end select
        endif
      enddo
      enddo
      enddo
    enddo
    if(.not.valid_flux(candidate)) return
    if(present(internal_energy_flux)) then
      if(.not.valid_eb_integrated_flux_3d(held_candidate,geometry,1)) return
    endif
    ! Check represented budgets as well as the analytic allowance; a tiny
    ! negative roundoff is a rejected transaction, never a species clip.
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      if(geometry%volume_fraction(i,j,k)<=0) cycle
      do c=1,nc
        component=c+4
        if(c==1) component=1
        budget=available(c,i,j,k)
        do axis=1,3
          face=[i,j,k];face(axis)=face(axis)-1
          budget=budget+dt*areas(axis)*aperture(axis,face)* &
            (value(candidate,axis,component,face)-value(baseline,axis,component,face))
          face(axis)=face(axis)+1
          budget=budget-dt*areas(axis)*aperture(axis,face)* &
            (value(candidate,axis,component,face)-value(baseline,axis,component,face))
        enddo
        if(.not.ieee_is_finite(budget).or.budget<0) return
      enddo
    enddo
    enddo
    enddo
    call move_alloc(candidate%x,flux%x);call move_alloc(candidate%y,flux%y);call move_alloc(candidate%z,flux%z)
    if(present(alpha_output)) call move_alloc(alpha,alpha_output)
    if(present(internal_energy_flux)) then
      call move_alloc(held_candidate%x,internal_energy_flux%x)
      call move_alloc(held_candidate%y,internal_energy_flux%y)
      call move_alloc(held_candidate%z,internal_energy_flux%z)
    endif
    ok=.true.
  contains
    logical function valid_flux(f) result(valid)
      type(eb_integrated_flux_3d), intent(in) :: f
      valid=.false.
      if(.not.allocated(f%x).or..not.allocated(f%y).or..not.allocated(f%z)) return
      if(any(shape(f%x)/=[nv,n(1)+1,n(2),n(3)]).or. &
        any(shape(f%y)/=[nv,n(1),n(2)+1,n(3)]).or.any(shape(f%z)/=[nv,n(1),n(2),n(3)+1])) return
      if(any(lbound(f%x)/=[1,0,1,1]).or.any(lbound(f%y)/=[1,1,0,1]).or. &
        any(lbound(f%z)/=[1,1,1,0])) return
      valid=all(ieee_is_finite(f%x)).and.all(ieee_is_finite(f%y)).and.all(ieee_is_finite(f%z))
    end function valid_flux
    real(dp) function value(f,d,index_component,index) result(v)
      type(eb_integrated_flux_3d), intent(in) :: f
      integer, intent(in) :: d,index_component,index(3)
      select case(d)
      case(1);v=f%x(index_component,index(1),index(2),index(3))
      case(2);v=f%y(index_component,index(1),index(2),index(3))
      case(3);v=f%z(index_component,index(1),index(2),index(3))
      end select
    end function value
    real(dp) function aperture(d,index) result(v)
      integer, intent(in) :: d,index(3)
      select case(d)
      case(1);v=geometry%x_face_fraction(index(1),index(2),index(3))
      case(2);v=geometry%y_face_fraction(index(1),index(2),index(3))
      case(3);v=geometry%z_face_fraction(index(1),index(2),index(3))
      end select
    end function aperture
  end subroutine limit_reactive_eb_flux_correction_3d

  logical function valid_plm_boundary(boundary) result(valid)
    type(eb_plm_boundary_3d), optional, intent(in) :: boundary
    integer :: face
    valid=.true.
    if(.not.present(boundary)) return
    valid=.false.
    if(any(boundary%lower>boundary%upper)) return
    do face=1,6
      select case(trim(boundary%faces(face)))
      case('interior','periodic','inflow','outflow','slip_wall','no_slip_wall','isothermal_wall','symmetry')
      case default
        return
      end select
    end do
    valid=.true.
  end function valid_plm_boundary

  ! pc_cmpflx corrects only the local Riemann arguments. Reconstructed states,
  ! species fractions and independently traced internal energy remain unchanged.
  subroutine correct_plm_boundary(left,right,axis,face,boundary)
    real(dp), intent(inout) :: left(:),right(:)
    integer, intent(in) :: axis,face
    type(eb_plm_boundary_3d), optional, intent(in) :: boundary
    integer :: side
    logical :: wall
    if(.not.present(boundary)) return
    if(face==boundary%lower(axis)-1) then
      side=2*axis-1
    else if(face==boundary%upper(axis)) then
      side=2*axis
    else
      return
    end if
    select case(trim(boundary%faces(side)))
    case('slip_wall','no_slip_wall','isothermal_wall','symmetry');wall=.true.
    case('outflow');wall=.false.
    case default;return
    end select
    if(mod(side,2)==1) then
      left(1:5)=right(1:5)
      if(wall) left(axis+1)=-right(axis+1)
    else
      right(1:5)=left(1:5)
      if(wall) right(axis+1)=-left(axis+1)
    end if
  end subroutine correct_plm_boundary

  ! Cell flattening for the normal predictor. Inputs contain three halos;
  ! outputs include its one donor layer. Missing pencils and cut-cell pencils
  ! retain unity, following flatten_eb. Covered storage is never inspected.
  subroutine compute_reactive_eb_flattening_3d(species,ghost,temperature,volume_fraction,n,flattening,ok)
    type(nasa7_species), intent(in) :: species(:)
    integer, intent(in) :: n(3)
    real(dp), intent(in) :: ghost(:,-2:,-2:,-2:),temperature(-2:,-2:,-2:),volume_fraction(-2:,-2:,-2:)
    real(dp), intent(inout) :: flattening(0:,0:,0:)
    logical, intent(out) :: ok
    real(dp), allocatable :: candidate(:,:,:)
    real(dp) :: pressure(-3:3),velocity(-3:3),vf(-3:3),q(reactive_nvar(size(species))),recovered,speed
    integer :: i,j,k,axis,d,idx(3),point(3)
    logical :: valid
    ok=.false.
    if(any(n<1)) return
    if(any(shape(ghost)/=[reactive_nvar(size(species)),n+6]).or.any(shape(temperature)/=n+6).or. &
      any(shape(volume_fraction)/=n+6).or.any(shape(flattening)/=n+2)) return
    if(.not.all(ieee_is_finite(volume_fraction))) return
    if(any(volume_fraction<0).or.any(volume_fraction>1)) return
    allocate(candidate(0:n(1)+1,0:n(2)+1,0:n(3)+1));candidate=1
    do k=0,n(3)+1
    do j=0,n(2)+1
    do i=0,n(1)+1
      idx=[i,j,k]
      ! Normal prediction uses face donors, never exterior edge/corner donors.
      if(count(idx<1.or.idx>n)>1.or.volume_fraction(i,j,k)<=0) cycle
      do axis=1,3
        if(idx(axis)-3< -2.or.idx(axis)+3>n(axis)+3) cycle
        do d=-3,3
          point=idx;point(axis)=point(axis)+d
          vf(d)=volume_fraction(point(1),point(2),point(3))
        end do
        ! Cut cells disable the reference shock sensor in this direction.
        ! Also avoid undefined thermodynamics in a fully covered pencil.
        if(any(vf<1)) cycle
        do d=-3,3
          point=idx;point(axis)=point(axis)+d
          call reactive_conserved_to_primitive(species,ghost(:,point(1),point(2),point(3)), &
            temperature(point(1),point(2),point(3)),q,recovered,speed,valid)
          if(.not.valid) return
          pressure(d)=q(5);velocity(d)=q(axis+1)
        end do
        candidate(i,j,k)=min(candidate(i,j,k),reactive_ppm_flattening_coefficient(pressure,velocity))
      end do
    end do
    end do
    end do
    flattening=candidate;ok=.true.
  end subroutine compute_reactive_eb_flattening_3d

  ! Application-facing outflow step. Extended geometry contains six physical
  ! halo cells; use domain-face extension of the source geometry, not cell-metric
  ! replication. State donors are nearest physical cells in each direction.
  subroutine advance_reactive_eb_outflow_ctu_3d(species,state,temperature,extended,dt,new_state,new_temperature,ok, &
      redistribution,difmag,integrated_flux,boundary_change,wall_change,target_volume_fraction,order,shock_flattening)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),dt
    type(eb_geometry_3d), intent(in) :: extended
    real(dp), intent(inout) :: new_state(:,:,:,:),new_temperature(:,:,:)
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: redistribution
    real(dp), optional, intent(in) :: difmag,target_volume_fraction
    integer, optional, intent(in) :: order
    logical, optional, intent(in) :: shock_flattening
    type(eb_integrated_flux_3d), optional, intent(inout) :: integrated_flux
    real(dp), optional, intent(inout) :: boundary_change(:),wall_change(:)
    call advance_reactive_eb_boundary_ctu_3d(species,state,temperature,extended,dt,new_state,new_temperature,ok, &
      [character(len=7)::'outflow','outflow','outflow','outflow','outflow','outflow'], &
      redistribution,difmag,integrated_flux,boundary_change,wall_change,target_volume_fraction,order,shock_flattening)
  end subroutine advance_reactive_eb_outflow_ctu_3d

  subroutine compute_sdc_eb_boundary_hydro_3d(species,state,temperature,extended,dt,conserved_source, &
      hydro_rate,integrated_flux,boundary_change,wall_change,ok,boundary_faces,difmag,order,shock_flattening, &
      inflow_state,inflow_temperature,baseline_flux,internal_energy_state,internal_energy_rate, &
      internal_energy_integrated_flux,internal_energy_baseline_flux,internal_energy_pressure_work,inflow_internal_energy)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: extended
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),dt,conserved_source(:,:,:,:)
    real(dp), intent(inout) :: hydro_rate(:,:,:,:),boundary_change(:),wall_change(:)
    type(eb_integrated_flux_3d), intent(inout) :: integrated_flux
    logical, intent(out) :: ok
    character(len=*), intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: difmag,inflow_state(:,:),inflow_temperature(:)
    integer, optional, intent(in) :: order
    logical, optional, intent(in) :: shock_flattening
    type(eb_integrated_flux_3d), optional, intent(inout) :: baseline_flux,internal_energy_integrated_flux, &
      internal_energy_baseline_flux
    real(dp), optional, intent(in) :: internal_energy_state(:,:,:),inflow_internal_energy(:)
    real(dp), optional, intent(inout) :: internal_energy_rate(:,:,:),internal_energy_pressure_work(:,:,:)
    real(dp), allocatable :: unused(:,:,:,:),unused_t(:,:,:)
    allocate(unused,mold=state);allocate(unused_t,mold=temperature)
    call advance_reactive_eb_boundary_ctu_3d(species,state,temperature,extended,dt,unused,unused_t,ok, &
      boundary_faces,difmag=difmag,integrated_flux=integrated_flux,boundary_change=boundary_change, &
      wall_change=wall_change,order=order,shock_flattening=shock_flattening,inflow_state=inflow_state, &
      inflow_temperature=inflow_temperature,conserved_source=conserved_source,hydro_rate=hydro_rate, &
      baseline_flux=baseline_flux,internal_energy_state=internal_energy_state,internal_energy_rate=internal_energy_rate, &
      internal_energy_integrated_flux=internal_energy_integrated_flux, &
      internal_energy_baseline_flux=internal_energy_baseline_flux,internal_energy_pressure_work=internal_energy_pressure_work, &
      inflow_internal_energy=inflow_internal_energy)
  end subroutine compute_sdc_eb_boundary_hydro_3d

  subroutine advance_reactive_eb_boundary_ctu_3d(species,state,temperature,extended,dt, &
      new_state,new_temperature,ok,boundary_faces, &
      redistribution,difmag,integrated_flux,boundary_change,wall_change,target_volume_fraction,order,shock_flattening, &
      inflow_state,inflow_temperature,conserved_source,hydro_rate, &
      state_redist_max_order,baseline_flux,internal_energy_state,internal_energy_rate,internal_energy_integrated_flux, &
      internal_energy_baseline_flux,internal_energy_pressure_work,inflow_internal_energy)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),dt
    type(eb_geometry_3d), intent(in) :: extended
    real(dp), intent(inout) :: new_state(:,:,:,:),new_temperature(:,:,:)
    integer, optional, intent(in) :: state_redist_max_order
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: redistribution
    real(dp), optional, intent(in) :: difmag,target_volume_fraction
    integer, optional, intent(in) :: order
    logical, optional, intent(in) :: shock_flattening
    type(eb_integrated_flux_3d), optional, intent(inout) :: integrated_flux
    real(dp), optional, intent(inout) :: boundary_change(:),wall_change(:)
    character(len=*), intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:)
    real(dp), optional, intent(in) :: conserved_source(:,:,:,:)
    real(dp), optional, intent(inout) :: hydro_rate(:,:,:,:)
    type(eb_integrated_flux_3d), optional, intent(inout) :: baseline_flux,internal_energy_integrated_flux, &
      internal_energy_baseline_flux
    real(dp), optional, intent(in) :: internal_energy_state(:,:,:),inflow_internal_energy(:)
    real(dp), optional, intent(inout) :: internal_energy_rate(:,:,:),internal_energy_pressure_work(:,:,:)
    real(dp), allocatable :: held_ghost(:,:,:)
    integer :: held_options
    logical :: held
    real(dp), allocatable :: source_ghost(:,:,:,:)
    character(len=32) :: source_faces(6)
    real(dp) :: shift(3)
    logical :: reflect(3)
    type(eb_geometry_3d) :: working
    type(eb_plm_boundary_3d) :: boundary
    real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),vf(:,:,:),flat(:,:,:),candidate(:,:,:,:),ct(:,:,:)
    integer :: n(3),nv,i,j,k,donor(3),inlet,axis
    logical :: valid,use_flattening
    ok=.false.;nv=reactive_nvar(size(species));n=shape(temperature)
    held_options=count([present(internal_energy_state),present(internal_energy_rate),present(internal_energy_integrated_flux)])
    if(held_options/=0.and.held_options/=3) return
    held=held_options==3
    if(held) then
      if(.not.present(hydro_rate)) return
      if(any(shape(internal_energy_state)/=n).or.any(shape(internal_energy_rate)/=n)) return
      if(present(internal_energy_pressure_work)) then
        if(any(shape(internal_energy_pressure_work)/=n)) return
      endif
      if(present(internal_energy_baseline_flux).and..not.present(baseline_flux)) return
      if(present(inflow_internal_energy)) then
        if(size(inflow_internal_energy)/=6) return
        if(.not.all(ieee_is_finite(inflow_internal_energy))) return
      endif
    else
      if(present(internal_energy_pressure_work).or.present(internal_energy_baseline_flux).or. &
        present(inflow_internal_energy)) return
    endif
    if(.not.valid_reactive_boundaries_3d(boundary_faces)) return
    if(any(boundary_faces=='inflow')) then
      if(.not.present(inflow_state).or..not.present(inflow_temperature)) return
      if(any(shape(inflow_state)/=[nv,6]).or.size(inflow_temperature)/=6) return
      if(.not.all(ieee_is_finite(inflow_state)).or..not.all(ieee_is_finite(inflow_temperature))) return
    end if
    if(any(n<1).or.any(shape(state)/=[nv,n])) return
    if(present(conserved_source)) then
      if(any(shape(conserved_source)/=[nv,n])) return
    end if
    if(present(hydro_rate)) then
      if(.not.present(conserved_source).or.any(shape(hydro_rate)/=[nv,n])) return
    end if
    if(any(shape(new_state)/=[nv,n]).or.any(shape(new_temperature)/=n)) return
    if(any([extended%nx,extended%ny,extended%nz]/=n+12)) return
    call extract_eb_geometry_3d(extended,[4,4,4],n+9,working,valid)
    if(.not.valid) return
    allocate(ghost(nv,-2:n(1)+9,-2:n(2)+9,-2:n(3)+9),gt(-2:n(1)+9,-2:n(2)+9,-2:n(3)+9), &
      vf(-2:n(1)+9,-2:n(2)+9,-2:n(3)+9),flat(0:n(1)+7,0:n(2)+7,0:n(3)+7))
    vf=extended%volume_fraction;ghost=0;gt=0;flat=1
    if(held) allocate(held_ghost(-2:n(1)+9,-2:n(2)+9,-2:n(3)+9),source=0._dp)
    do k=-2,n(3)+9
    do j=-2,n(2)+9
    do i=-2,n(1)+9
      if(vf(i,j,k)<=0) cycle
      call map_reactive_boundary_donor_3d([i,j,k]-3,n,boundary_faces,.false.,donor,reflect,shift,inlet,valid)
      if(.not.valid) return
      if(inlet>0) then
        ghost(:,i,j,k)=inflow_state(:,inlet);gt(i,j,k)=inflow_temperature(inlet)
        if(held) then
          if(present(inflow_internal_energy)) then
            held_ghost(i,j,k)=inflow_internal_energy(inlet)
          else
            ! A prescribed gas inlet initializes its UEINT consistently.
            ! Explicit independent inlet UEINT may instead be supplied.
            if(inflow_state(1,inlet)<=0) return
            held_ghost(i,j,k)=inflow_state(5,inlet)-.5_dp*sum(inflow_state(2:4,inlet)**2)/inflow_state(1,inlet)
          endif
        endif
      else
        ! Open target geometry cannot obtain fluid state from a covered donor.
        if(extended%volume_fraction(donor(1)+6,donor(2)+6,donor(3)+6)<=0) return
        ghost(:,i,j,k)=state(:,donor(1),donor(2),donor(3))
        gt(i,j,k)=temperature(donor(1),donor(2),donor(3))
        if(held) held_ghost(i,j,k)=internal_energy_state(donor(1),donor(2),donor(3))
      end if
      do axis=1,3
        if(reflect(axis)) ghost(axis+1,i,j,k)=-ghost(axis+1,i,j,k)
      end do
    end do
    end do
    end do
    if(present(conserved_source)) then
      allocate(source_ghost(nv,n(1)+6,n(2)+6,n(3)+6));source_ghost=0
      source_faces=boundary_faces
      where(source_faces/='periodic') source_faces='outflow'
      do k=1,n(3)+6
      do j=1,n(2)+6
      do i=1,n(1)+6
        if(working%volume_fraction(i,j,k)<=0) cycle
        call map_reactive_boundary_donor_3d([i,j,k]-3,n,source_faces,.false.,donor,reflect,shift,inlet,valid)
        if(.not.valid) return
        if(extended%volume_fraction(donor(1)+6,donor(2)+6,donor(3)+6)<=0) return
        source_ghost(:,i,j,k)=conserved_source(:,donor(1),donor(2),donor(3))
      end do
      end do
      end do
    end if
    use_flattening=.true.
    if(present(shock_flattening)) use_flattening=shock_flattening
    if(use_flattening) then
      call compute_reactive_eb_flattening_3d(species,ghost,gt,vf,n+6,flat,valid)
      if(.not.valid) return
    end if
    boundary%lower=4;boundary%upper=n+3;boundary%faces=boundary_faces
    candidate=state;ct=temperature
    call advance_reactive_eb_plm_ctu_3d(species,working,ghost,gt,vf,dt,flat,candidate,ct,ok, &
      redistribution,difmag,integrated_flux,boundary_change,wall_change,target_volume_fraction,order,boundary, &
      source_ghost,hydro_rate, &
      state_redist_max_order=state_redist_max_order,baseline_flux=baseline_flux, &
      internal_energy_state=held_ghost,internal_energy_rate=internal_energy_rate, &
      internal_energy_integrated_flux=internal_energy_integrated_flux, &
      internal_energy_baseline_flux=internal_energy_baseline_flux,internal_energy_pressure_work=internal_energy_pressure_work)
    if(.not.ok.or.present(hydro_rate)) return
    ! Existing application and split operators preserve covered-cell storage.
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      if(vf(i+3,j+3,k+3)>0) cycle
      candidate(:,i,j,k)=state(:,i,j,k);ct(i,j,k)=temperature(i,j,k)
    end do
    end do
    end do
    new_state=candidate;new_temperature=ct
  end subroutine advance_reactive_eb_boundary_ctu_3d

  ! Complete hydro update for the central cells 4:N-3 of a prepared working
  ! geometry. Physical boundary filling/flattening remain caller responsibilities.
  ! Uses the same extended geometry/ghost layout as the CTU advance.
  ! Returns conservative rates before redistribution and time-integrated flux densities.
  ! Pinned Hydro.H pc_eb_pdivu, with one geometry/face halo around output.
  ! Regular cells use face-pressure averages; cut cells use cell pressure and
  ! owner-cell connectivity for the centroid velocity stencil on BOTH faces.
  ! This is a non-face source in raw UEINT, before ML redistribution.
  subroutine compute_eb_godunov_pressure_work_3d(geometry,prediction,cell_pressure,pressure_work,ok)
    type(eb_geometry_3d), intent(in) :: geometry
    type(eb_plm_face_prediction_3d), intent(in) :: prediction(3)
    real(dp), intent(in) :: cell_pressure(:,:,:)
    real(dp), intent(inout) :: pressure_work(:,:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: candidate(:,:,:)
    real(dp) :: dxinv(3),components(3),value,vm,vp,pm,pp,apm,app,weighted,term,factor
    integer :: n(3),axis,i,j,k,owner(3),minus(3),plus(3),expected(3),lower(3)
    logical :: good
    ok=.false.;n=[geometry%nx,geometry%ny,geometry%nz]
    if(.not.geometry%is_valid().or.any(n<3)) return
    if(any(shape(cell_pressure)/=n).or.any(shape(pressure_work)/=n-2)) return
    if(min(geometry%dx,geometry%dy,geometry%dz)<1._dp/huge(1._dp)) return
    dxinv=1._dp/[geometry%dx,geometry%dy,geometry%dz]
    do axis=1,3
      if(.not.allocated(prediction(axis)%pressure).or..not.allocated(prediction(axis)%velocity)) return
      expected=n;expected(axis)=expected(axis)+1;lower=1;lower(axis)=0
      if(any(shape(prediction(axis)%pressure)/=expected).or.any(shape(prediction(axis)%velocity)/=expected)) return
      if(any(lbound(prediction(axis)%pressure)/=lower).or.any(lbound(prediction(axis)%velocity)/=lower)) return
    enddo
    allocate(candidate(n(1)-2,n(2)-2,n(3)-2),source=0._dp)
    do k=2,n(3)-1
    do j=2,n(2)-1
    do i=2,n(1)-1
      if(geometry%volume_fraction(i,j,k)<=0) cycle
      owner=[i,j,k];components=0
      if(geometry%cell_type(i,j,k)==eb_regular_cell_3d) then
        do axis=1,3
          minus=owner;minus(axis)=minus(axis)-1;plus=owner
          vm=prediction(axis)%velocity(minus(1),minus(2),minus(3))
          vp=prediction(axis)%velocity(plus(1),plus(2),plus(3))
          pm=prediction(axis)%pressure(minus(1),minus(2),minus(3))
          pp=prediction(axis)%pressure(plus(1),plus(2),plus(3))
          if(.not.all(ieee_is_finite([vm,vp,pm,pp]))) return
          call add_value(vp,-vm,value,good);if(.not.good) return
          call add_value(pm,pp,weighted,good);if(.not.good) return
          call multiply_value(weighted,value,term,good);if(.not.good) return
          call multiply_value(term,dxinv(axis),components(axis),good);if(.not.good) return
        enddo
        call add_value(components(1),components(2),value,good);if(.not.good) return
        call add_value(value,components(3),weighted,good);if(.not.good) return
        value=.5_dp*weighted
      else
        if(.not.ieee_is_finite(cell_pressure(i,j,k)).or.cell_pressure(i,j,k)<=0) return
        do axis=1,3
          minus=owner;minus(axis)=minus(axis)-1;plus=owner
          call face_velocity(axis,minus,vm,apm,good);if(.not.good) return
          call face_velocity(axis,plus,vp,app,good);if(.not.good) return
          call add_value(app*vp,-apm*vm,value,good);if(.not.good) return
          call multiply_value(dxinv(axis),value,components(axis),good);if(.not.good) return
        enddo
        call add_value(components(1),components(2),value,good);if(.not.good) return
        call add_value(value,components(3),weighted,good);if(.not.good) return
        if(geometry%volume_fraction(i,j,k)<1._dp/huge(1._dp)) return
        call multiply_value(1._dp/geometry%volume_fraction(i,j,k),cell_pressure(i,j,k),factor,good)
        if(.not.good) return
        call multiply_value(factor,weighted,value,good);if(.not.good) return
      endif
      if(.not.ieee_is_finite(value)) return
      candidate(i-1,j-1,k-1)=value
    enddo
    enddo
    enddo
    pressure_work=candidate;ok=.true.
  contains
    subroutine face_velocity(direction,face,value,aperture,good)
      integer, intent(in) :: direction,face(3)
      real(dp), intent(out) :: value,aperture
      logical, intent(out) :: good
      real(dp) :: offset(2),fraction(2),weight,next_value,accumulated
      integer :: b,c,shift(2),cell(3),donor(3),s,t
      good=.false.;value=0;fraction=0
      b=merge(2,1,direction==1);c=merge(2,3,direction==3)
      select case(direction)
      case(1)
        aperture=geometry%x_face_fraction(face(1),face(2),face(3))
        offset=[geometry%x_face_centroid_y(face(1),face(2),face(3)), &
          geometry%x_face_centroid_z(face(1),face(2),face(3))]
      case(2)
        aperture=geometry%y_face_fraction(face(1),face(2),face(3))
        offset=[geometry%y_face_centroid_x(face(1),face(2),face(3)), &
          geometry%y_face_centroid_z(face(1),face(2),face(3))]
      case(3)
        aperture=geometry%z_face_fraction(face(1),face(2),face(3))
        offset=[geometry%z_face_centroid_x(face(1),face(2),face(3)), &
          geometry%z_face_centroid_y(face(1),face(2),face(3))]
      end select
      if(aperture<=0) then
        good=.true.;return
      endif
      if(aperture<1) then
        ! The pinned 3-D kernel uses the positive-side cell flag for each
        ! face, including the upper face. Only axial tangential connections
        ! enter this stencil; their bits equal that cell's open apertures.
        cell=face;cell(direction)=cell(direction)+1
        shift=int(sign(1._dp,offset))
        if(tangent_connected(cell,b,shift(1))) fraction(1)=abs(offset(1))
        if(tangent_connected(cell,c,shift(2))) fraction(2)=abs(offset(2))
      else
        shift=0
      endif
      do t=0,1
      do s=0,1
        weight=merge(fraction(1),1-fraction(1),s==1)*merge(fraction(2),1-fraction(2),t==1)
        if(weight==0) cycle
        donor=face;donor(b)=donor(b)+s*shift(1);donor(c)=donor(c)+t*shift(2)
        next_value=prediction(direction)%velocity(donor(1),donor(2),donor(3))
        if(.not.ieee_is_finite(next_value)) return
        call add_value(value,weight*next_value,accumulated,good);if(.not.good) return
        value=accumulated
      enddo
      enddo
      good=.true.
    end subroutine face_velocity
    logical function tangent_connected(cell,direction,shift) result(connected)
      integer, intent(in) :: cell(3),direction,shift
      integer :: face(3)
      connected=.false.
      if(geometry%volume_fraction(cell(1),cell(2),cell(3))<=0) return
      face=cell
      if(shift<0) face(direction)=face(direction)-1
      select case(direction)
      case(1);connected=geometry%x_face_fraction(face(1),face(2),face(3))>0
      case(2);connected=geometry%y_face_fraction(face(1),face(2),face(3))>0
      case(3);connected=geometry%z_face_fraction(face(1),face(2),face(3))>0
      end select
    end function tangent_connected
    subroutine add_value(a,b,c,good)
      real(dp), intent(in) :: a,b
      real(dp), intent(out) :: c
      logical, intent(out) :: good
      good=.false.
      if(b>0) then
        if(a>huge(1._dp)-b) return
      else if(b<0) then
        if(a< -huge(1._dp)-b) return
      endif
      c=a+b;good=ieee_is_finite(c)
    end subroutine add_value
    subroutine multiply_value(a,b,c,good)
      real(dp), intent(in) :: a,b
      real(dp), intent(out) :: c
      logical, intent(out) :: good
      good=.false.
      if(abs(a)>1) then
        if(abs(b)>huge(1._dp)/abs(a)) return
      else if(abs(b)>1) then
        if(abs(a)>huge(1._dp)/abs(b)) return
      endif
      c=a*b;good=ieee_is_finite(c)
    end subroutine multiply_value
  end subroutine compute_eb_godunov_pressure_work_3d

  subroutine compute_sdc_eb_hydro_3d(species,geometry,ghost,temperature,volume_fraction,dt,flattening, &
      conserved_source,hydro_rate,integrated_flux,boundary_change,wall_change,ok,difmag,order,boundary,baseline_flux, &
      internal_energy_state,internal_energy_rate,internal_energy_integrated_flux, &
      internal_energy_baseline_flux,internal_energy_pressure_work,spatial_flux_weight)
    real(dp), optional, intent(in) :: spatial_flux_weight
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: ghost(:,-2:,-2:,-2:),temperature(-2:,-2:,-2:),volume_fraction(-2:,-2:,-2:)
    real(dp), intent(in) :: dt,flattening(0:,0:,0:),conserved_source(:,:,:,:)
    real(dp), intent(inout) :: hydro_rate(:,:,:,:),boundary_change(:),wall_change(:)
    type(eb_integrated_flux_3d), intent(inout) :: integrated_flux
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: difmag
    integer, optional, intent(in) :: order
    type(eb_plm_boundary_3d), optional, intent(in) :: boundary
    type(eb_integrated_flux_3d), optional, intent(inout) :: baseline_flux
    real(dp), optional, intent(in) :: internal_energy_state(-2:,-2:,-2:)
    real(dp), optional, intent(inout) :: internal_energy_rate(:,:,:),internal_energy_pressure_work(:,:,:)
    type(eb_integrated_flux_3d), optional, intent(inout) :: internal_energy_integrated_flux,internal_energy_baseline_flux
    real(dp), allocatable :: unused_state(:,:,:,:),unused_temperature(:,:,:)
    integer :: n(3),nv
    ok=.false.;n=[geometry%nx,geometry%ny,geometry%nz];nv=reactive_nvar(size(species))
    if(any(n<7)) return
    allocate(unused_state(nv,n(1)-6,n(2)-6,n(3)-6),unused_temperature(n(1)-6,n(2)-6,n(3)-6))
    call advance_reactive_eb_plm_ctu_3d(species,geometry,ghost,temperature,volume_fraction,dt,flattening, &
      unused_state,unused_temperature,ok,difmag=difmag,integrated_flux=integrated_flux, &
      boundary_change=boundary_change,wall_change=wall_change,order=order,boundary=boundary, &
      conserved_source=conserved_source,hydro_rate=hydro_rate,baseline_flux=baseline_flux, &
      internal_energy_state=internal_energy_state,internal_energy_rate=internal_energy_rate, &
      internal_energy_integrated_flux=internal_energy_integrated_flux, &
      internal_energy_baseline_flux=internal_energy_baseline_flux,internal_energy_pressure_work=internal_energy_pressure_work, &
      spatial_flux_weight=spatial_flux_weight)
  end subroutine compute_sdc_eb_hydro_3d

  subroutine advance_reactive_eb_plm_ctu_3d(species,geometry,ghost,temperature,volume_fraction,dt,flattening, &
      new_state,new_temperature,ok,redistribution,difmag,integrated_flux,boundary_change,wall_change, &
      target_volume_fraction,order,boundary,conserved_source,hydro_rate,baseline_flux, &
      state_redist_max_order,internal_energy_state,internal_energy_rate,internal_energy_integrated_flux, &
      internal_energy_baseline_flux,internal_energy_pressure_work,spatial_flux_weight)
    real(dp), optional, intent(in) :: spatial_flux_weight
    type(nasa7_species), intent(in) :: species(:)
    type(eb_plm_boundary_3d), optional, intent(in) :: boundary
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: ghost(:,-2:,-2:,-2:),temperature(-2:,-2:,-2:),volume_fraction(-2:,-2:,-2:)
    real(dp), intent(in) :: dt,flattening(0:,0:,0:)
    real(dp), intent(inout) :: new_state(:,:,:,:),new_temperature(:,:,:)
    integer, optional, intent(in) :: state_redist_max_order
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: redistribution
    real(dp), optional, intent(in) :: difmag,target_volume_fraction
    integer, optional, intent(in) :: order
    type(eb_integrated_flux_3d), optional, intent(inout) :: integrated_flux
    real(dp), optional, intent(inout) :: boundary_change(:),wall_change(:)
    real(dp), optional, intent(in) :: conserved_source(:,:,:,:)
    real(dp), optional, intent(inout) :: hydro_rate(:,:,:,:)
    type(eb_integrated_flux_3d), optional, intent(inout) :: baseline_flux
    real(dp), optional, intent(in) :: internal_energy_state(-2:,-2:,-2:)
    real(dp), optional, intent(inout) :: internal_energy_rate(:,:,:),internal_energy_pressure_work(:,:,:)
    type(eb_integrated_flux_3d), optional, intent(inout) :: internal_energy_integrated_flux,internal_energy_baseline_flux
    type(eb_integrated_flux_3d) :: held_raw,held_flux,held_low
    real(dp), allocatable :: held_rhs(:,:,:),pwork(:,:,:),cell_pressure(:,:,:)
    real(dp) :: held_volume
    integer :: held_options
    logical :: held
    type(eb_plm_face_prediction_3d) :: prediction(3),low_prediction(3)
    type(eb_geometry_3d) :: margin,physical,redist_halo
    type(eb_integrated_flux_3d) :: flux,low_flux,blended_flux,blended_held_flux
    real(dp) :: face_weight
    real(dp), allocatable :: state(:,:,:,:),temp(:,:,:),candidate(:,:,:,:),ct(:,:,:),rhs(:,:,:,:)
    real(dp) :: coefficient,target,total,recovered,speed,q(reactive_nvar(size(species)))
    real(dp) :: inventory_boundary(reactive_nvar(size(species))),wall(reactive_nvar(size(species)))
    integer :: n(3),m(3),nv,i,j,k,a,lo(4),hi(4),wall_lo(3),wall_hi(3)
    character(len=32) :: redist_faces(6)
    character(len=16) :: method
    logical :: valid,no_slip(6)
    ok=.false.;nv=reactive_nvar(size(species));n=[geometry%nx,geometry%ny,geometry%nz]
    if(any(n<7).or..not.geometry%is_valid()) return
    face_weight=1
    if(present(spatial_flux_weight)) then
      if(.not.present(hydro_rate)) return
      face_weight=spatial_flux_weight
      if(.not.ieee_is_finite(face_weight)) return
      if(face_weight<0.or.face_weight>1) return
    endif
    held_options=count([present(internal_energy_state),present(internal_energy_rate),present(internal_energy_integrated_flux)])
    if(held_options/=0.and.held_options/=3) return
    held=held_options==3
    if(held) then
      ! The tree integrator owns the stored scalar and energy reset. This
      ! prepared patch operator returns its raw spatial family atomically.
      if(.not.present(hydro_rate)) return
      if(any(shape(internal_energy_state)/=n+6).or.any(shape(internal_energy_rate)/=n-6)) return
      if(present(internal_energy_pressure_work)) then
        if(any(shape(internal_energy_pressure_work)/=n-6)) return
      endif
      if(present(internal_energy_baseline_flux).and..not.present(baseline_flux)) return
    else
      if(present(internal_energy_pressure_work).or.present(internal_energy_baseline_flux)) return
    endif
    if(any(shape(new_state)/=[nv,n-6]).or.any(shape(new_temperature)/=n-6)) return
    if(.not.ieee_is_finite(dt).or.dt<=0) return
    coefficient=.1_dp;target=.5_dp;method='none'
    if(present(difmag)) coefficient=difmag
    if(present(target_volume_fraction)) target=target_volume_fraction
    if(present(redistribution)) method=redistribution
    if(.not.all(ieee_is_finite([coefficient,target])).or.coefficient<0.or.target<=0.or.target>1) return
    if(method/='none'.and.method/='flux_redist'.and.method/='state_redist') return
    if(present(hydro_rate)) then
      if(method/='none'.or..not.present(conserved_source)) return
      if(any(shape(hydro_rate)/=[nv,n-6])) return
    end if
    if(present(baseline_flux).and..not.present(hydro_rate)) return
    if(present(boundary_change)) then
      if(size(boundary_change)/=nv) return
    end if
    if(present(wall_change)) then
      if(size(wall_change)/=nv) return
    end if
    call compute_reactive_eb_plm_ctu_fluxes_3d(species,geometry,ghost,temperature,volume_fraction, &
      dt,flattening,prediction,valid,order,boundary=boundary,conserved_source=conserved_source, &
      capture_internal_energy_flux=held)
    if(.not.valid) return
    call extract_eb_geometry_3d(geometry,[3,3,3],n-2,margin,valid)
    if(.not.valid) return
    m=n-4
    no_slip=.false.;wall_lo=0;wall_hi=m
    if(present(boundary)) then
      no_slip=boundary%faces=='no_slip_wall'.or.boundary%faces=='isothermal_wall'
      wall_lo=boundary%lower-3;wall_hi=boundary%upper-2
    end if
    if(held) then
    call adjust_reactive_artificial_viscosity_3d(ghost(:,2:n(1)-1,2:n(2)-1,2:n(3)-1), &
      m,[geometry%dx,geometry%dy,geometry%dz],coefficient,no_slip, &
      prediction(1)%flux,prediction(2)%flux,prediction(3)%flux,valid, &
      volume_fraction(2:n(1)-1,2:n(2)-1,2:n(3)-1), &
      margin%x_face_fraction>0,margin%y_face_fraction>0,margin%z_face_fraction>0,wall_lo,wall_hi, &
      internal_energy_state(2:n(1)-1,2:n(2)-1,2:n(3)-1), &
      prediction(1)%internal_energy_flux,prediction(2)%internal_energy_flux,prediction(3)%internal_energy_flux)
    else
    call adjust_reactive_artificial_viscosity_3d(ghost(:,2:n(1)-1,2:n(2)-1,2:n(3)-1), &
      m,[geometry%dx,geometry%dy,geometry%dz],coefficient,no_slip, &
      prediction(1)%flux,prediction(2)%flux,prediction(3)%flux,valid, &
      volume_fraction(2:n(1)-1,2:n(2)-1,2:n(3)-1), &
      margin%x_face_fraction>0,margin%y_face_fraction>0,margin%z_face_fraction>0,wall_lo,wall_hi)
    endif
    if(.not.valid) return
    ! The viscosity API deliberately treats zero as identity. Godunov still
    ! normalizes species flux when viscosity is disabled.
    if(coefficient==0) then
      do a=1,3
        lo=lbound(prediction(a)%flux);hi=ubound(prediction(a)%flux)
        do k=lo(4),hi(4)
        do j=lo(3),hi(3)
        do i=lo(2),hi(2)
          total=sum(prediction(a)%flux(6:,i,j,k))
          if(total/=0) prediction(a)%flux(6:,i,j,k)= &
            prediction(a)%flux(6:,i,j,k)*(prediction(a)%flux(1,i,j,k)/total)
        end do
        end do
        end do
      end do
    end if
    m=n-6
    allocate(flux%x(nv,0:m(1),m(2),m(3)),flux%y(nv,m(1),0:m(2),m(3)),flux%z(nv,m(1),m(2),0:m(3)))
    call interpolate_reactive_eb_centroid_fluxes_3d(margin,prediction(1)%flux,prediction(2)%flux,prediction(3)%flux, &
      flux%x,flux%y,flux%z,valid)
    if(.not.valid) return
    if(held) then
      allocate(held_raw%x(1,0:margin%nx,margin%ny,margin%nz), &
        held_raw%y(1,margin%nx,0:margin%ny,margin%nz),held_raw%z(1,margin%nx,margin%ny,0:margin%nz))
      held_raw%x(1,:,:,:)=prediction(1)%internal_energy_flux
      held_raw%y(1,:,:,:)=prediction(2)%internal_energy_flux
      held_raw%z(1,:,:,:)=prediction(3)%internal_energy_flux
      allocate(held_flux%x(1,0:m(1),m(2),m(3)),held_flux%y(1,m(1),0:m(2),m(3)), &
        held_flux%z(1,m(1),m(2),0:m(3)))
      call interpolate_reactive_eb_centroid_fluxes_3d(margin,held_raw%x,held_raw%y,held_raw%z, &
        held_flux%x,held_flux%y,held_flux%z,valid)
      if(.not.valid) return
    endif
    if(present(baseline_flux).or.face_weight<1) then
      ! Piecewise-constant normal fluxes are the candidate low-order baseline.
      ! Omit transverse forcing, viscosity and centroid corrections; those
      ! belong to the full correction. Its forward density/species budget
      ! still must be checked with the actual interval and complete source.
      call compute_reactive_eb_plm_normal_fluxes_3d(species,geometry,ghost,temperature,volume_fraction, &
        dt,flattening,low_prediction,valid,order=1,boundary=boundary,capture_internal_energy_flux=held)
      if(.not.valid) return
      allocate(low_flux%x(nv,0:m(1),m(2),m(3)),low_flux%y(nv,m(1),0:m(2),m(3)), &
        low_flux%z(nv,m(1),m(2),0:m(3)))
      low_flux%x=low_prediction(1)%flux(:,3:n(1)-3,4:n(2)-3,4:n(3)-3)
      low_flux%y=low_prediction(2)%flux(:,4:n(1)-3,3:n(2)-3,4:n(3)-3)
      low_flux%z=low_prediction(3)%flux(:,4:n(1)-3,4:n(2)-3,3:n(3)-3)
      low_flux%x=dt*low_flux%x;low_flux%y=dt*low_flux%y;low_flux%z=dt*low_flux%z
      if(any(.not.ieee_is_finite(low_flux%x)).or.any(.not.ieee_is_finite(low_flux%y)).or. &
        any(.not.ieee_is_finite(low_flux%z))) return
      if(present(internal_energy_baseline_flux).or.(held.and.face_weight<1)) then
        allocate(held_low%x(1,0:m(1),m(2),m(3)),held_low%y(1,m(1),0:m(2),m(3)),held_low%z(1,m(1),m(2),0:m(3)))
        held_low%x(1,:,:,:)=dt*low_prediction(1)%internal_energy_flux(3:n(1)-3,4:n(2)-3,4:n(3)-3)
        held_low%y(1,:,:,:)=dt*low_prediction(2)%internal_energy_flux(4:n(1)-3,3:n(2)-3,4:n(3)-3)
        held_low%z(1,:,:,:)=dt*low_prediction(3)%internal_energy_flux(4:n(1)-3,4:n(2)-3,3:n(3)-3)
        if(.not.all(ieee_is_finite(held_low%x)).or..not.all(ieee_is_finite(held_low%y)).or. &
          .not.all(ieee_is_finite(held_low%z))) return
      endif
    endif
    if(face_weight<1) then
      ! A single native interval coefficient acts on both face families.
      ! Work with integrated faces so zero returns the exact PCM baseline;
      ! their intensity determines advection, while pdivu stays native.
      blended_flux=low_flux
      if(face_weight>0) then
        blended_flux%x=(1-face_weight)*low_flux%x+face_weight*(dt*flux%x)
        blended_flux%y=(1-face_weight)*low_flux%y+face_weight*(dt*flux%y)
        blended_flux%z=(1-face_weight)*low_flux%z+face_weight*(dt*flux%z)
      endif
      flux%x=blended_flux%x/dt;flux%y=blended_flux%y/dt;flux%z=blended_flux%z/dt
      if(.not.all(ieee_is_finite(flux%x)).or..not.all(ieee_is_finite(flux%y)).or. &
        .not.all(ieee_is_finite(flux%z))) return
      if(held) then
        blended_held_flux=held_low
        if(face_weight>0) then
          blended_held_flux%x=(1-face_weight)*held_low%x+face_weight*(dt*held_flux%x)
          blended_held_flux%y=(1-face_weight)*held_low%y+face_weight*(dt*held_flux%y)
          blended_held_flux%z=(1-face_weight)*held_low%z+face_weight*(dt*held_flux%z)
        endif
        held_flux%x=blended_held_flux%x/dt;held_flux%y=blended_held_flux%y/dt;held_flux%z=blended_held_flux%z/dt
        if(.not.all(ieee_is_finite(held_flux%x)).or..not.all(ieee_is_finite(held_flux%y)).or. &
          .not.all(ieee_is_finite(held_flux%z))) return
      endif
    endif
    call extract_eb_geometry_3d(geometry,[4,4,4],n-3,physical,valid)
    if(.not.valid) return
    state=ghost(:,4:n(1)-3,4:n(2)-3,4:n(3)-3);temp=temperature(4:n(1)-3,4:n(2)-3,4:n(3)-3)
    do k=1,m(3)
    do j=1,m(2)
    do i=1,m(1)
      if(physical%volume_fraction(i,j,k)>0) cycle
      state(:,i,j,k)=0;temp(i,j,k)=0
    end do
    end do
    end do
    allocate(rhs(nv,m(1),m(2),m(3)))
    call reactive_eb_flux_divergence_3d(species,state,temp,physical,flux%x,flux%y,flux%z,rhs,valid)
    if(.not.valid) return
    if(held) then
      allocate(cell_pressure(margin%nx,margin%ny,margin%nz),source=0._dp)
      allocate(pwork(m(1),m(2),m(3)),held_rhs(m(1),m(2),m(3)),source=0._dp)
      do k=2,margin%nz-1
      do j=2,margin%ny-1
      do i=2,margin%nx-1
        if(margin%volume_fraction(i,j,k)<=0) cycle
        if(margin%cell_type(i,j,k)==eb_regular_cell_3d) cycle
        call reactive_conserved_to_primitive(species,ghost(:,i+2,j+2,k+2),temperature(i+2,j+2,k+2), &
          q,recovered,speed,valid)
        if(.not.valid) return
        cell_pressure(i,j,k)=q(5)
      enddo
      enddo
      enddo
      call compute_eb_godunov_pressure_work_3d(margin,prediction,cell_pressure,pwork,valid)
      if(.not.valid) return
      held_volume=physical%dx*physical%dy*physical%dz
      do k=1,m(3)
      do j=1,m(2)
      do i=1,m(1)
        if(physical%volume_fraction(i,j,k)<=0) cycle
        held_rhs(i,j,k)=( &
          physical%dy*physical%dz*(physical%x_face_fraction(i-1,j,k)*held_flux%x(1,i-1,j,k)- &
            physical%x_face_fraction(i,j,k)*held_flux%x(1,i,j,k))+ &
          physical%dx*physical%dz*(physical%y_face_fraction(i,j-1,k)*held_flux%y(1,i,j-1,k)- &
            physical%y_face_fraction(i,j,k)*held_flux%y(1,i,j,k))+ &
          physical%dx*physical%dy*(physical%z_face_fraction(i,j,k-1)*held_flux%z(1,i,j,k-1)- &
            physical%z_face_fraction(i,j,k)*held_flux%z(1,i,j,k)))/ &
          (held_volume*physical%volume_fraction(i,j,k))-pwork(i,j,k)
      enddo
      enddo
      enddo
      if(face_weight<1) then
        held_flux=blended_held_flux
      else
      held_flux%x=dt*held_flux%x;held_flux%y=dt*held_flux%y;held_flux%z=dt*held_flux%z
      endif
      if(.not.all(ieee_is_finite(held_rhs)).or..not.all(ieee_is_finite(held_flux%x)).or. &
        .not.all(ieee_is_finite(held_flux%y)).or..not.all(ieee_is_finite(held_flux%z))) return
    endif
    call hydro_inventory_changes(species,state,temp,physical,dt,flux,inventory_boundary,wall,valid)
    if(.not.valid) return
    if(present(hydro_rate)) then
      if(.not.all(ieee_is_finite(rhs))) return
      if(face_weight<1) then
        flux=blended_flux
      else
      flux%x=dt*flux%x;flux%y=dt*flux%y;flux%z=dt*flux%z
      endif
      if(.not.all(ieee_is_finite(flux%x)).or..not.all(ieee_is_finite(flux%y)).or. &
        .not.all(ieee_is_finite(flux%z))) return
      hydro_rate=rhs
      if(present(integrated_flux)) integrated_flux=flux
      if(present(boundary_change)) boundary_change=inventory_boundary
      if(present(wall_change)) wall_change=wall
      if(present(baseline_flux)) baseline_flux=low_flux
      if(held) then
        internal_energy_rate=held_rhs;internal_energy_integrated_flux=held_flux
        if(present(internal_energy_baseline_flux)) internal_energy_baseline_flux=held_low
        if(present(internal_energy_pressure_work)) internal_energy_pressure_work=pwork
      endif
      ok=.true.;return
    end if
    allocate(candidate,mold=state);allocate(ct,mold=temp)
    select case(method)
    case('flux_redist')
      call advance_reactive_eb_redistributed_3d(species,state,temp,physical,rhs,dt,candidate,ct,valid)
    case('state_redist')
      ! A prepared patch may have peer/ancestor halos rather than physical
      ! ghosts. Only a complete physical domain can use zero exterior dUdt.
      valid=.false.
      if(present(boundary)) valid=all(boundary%lower==4).and.all(boundary%upper==n-3)
      if(valid) then
        redist_faces=boundary%faces
        call extract_eb_geometry_3d(geometry,[2,2,2],n-1,redist_halo,valid)
        if(.not.valid) return
        call advance_reactive_eb_state_redistributed_3d(species,state,temp,physical,rhs,dt,candidate,ct,valid,target, &
          max_order=state_redist_max_order,boundary_state=ghost(:,2:n(1)-1,2:n(2)-1,2:n(3)-1), &
          boundary_geometry=redist_halo,boundary_faces=redist_faces)
      else
        call advance_reactive_eb_state_redistributed_3d(species,state,temp,physical,rhs,dt,candidate,ct,valid,target, &
          max_order=state_redist_max_order)
      endif
    case('none')
      candidate=state;ct=temp
      do k=1,m(3)
      do j=1,m(2)
      do i=1,m(1)
        if(physical%volume_fraction(i,j,k)<=0) cycle
        candidate(:,i,j,k)=state(:,i,j,k)+dt*rhs(:,i,j,k)
        call reactive_conserved_to_primitive(species,candidate(:,i,j,k),temp(i,j,k),q,recovered,speed,valid)
        if(.not.valid) return
        ct(i,j,k)=recovered
      end do
      end do
      end do
    end select
    if(.not.valid) return
    if(.not.all(ieee_is_finite(candidate)).or..not.all(ieee_is_finite(ct))) return
    flux%x=dt*flux%x;flux%y=dt*flux%y;flux%z=dt*flux%z
    if(.not.all(ieee_is_finite(flux%x)).or..not.all(ieee_is_finite(flux%y)).or. &
      .not.all(ieee_is_finite(flux%z))) return
    new_state=candidate;new_temperature=ct
    if(present(integrated_flux)) integrated_flux=flux
    if(present(boundary_change)) boundary_change=inventory_boundary
    if(present(wall_change)) wall_change=wall
    ok=.true.
  end subroutine advance_reactive_eb_plm_ctu_3d

  ! Twelve Riemann solves: three normal, six single-transverse, three final.
  ! Geometry includes a two-cell margin around the returned region. The normal
  ! predictor additionally needs three state/metric halos outside this geometry.
  ! Output indices describe cells 3:N-2, shifted down by two. Fluxes are intensive
  ! face-center values, before viscosity or centroid interpolation.
  subroutine compute_reactive_eb_plm_ctu_fluxes_3d(species,geometry,ghost,temperature,volume_fraction, &
      dt,flattening,prediction,ok,order,source,source_energy,boundary,conserved_source,normal_prediction, &
      capture_internal_energy_flux)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: ghost(:,-2:,-2:,-2:),temperature(-2:,-2:,-2:),volume_fraction(-2:,-2:,-2:)
    real(dp), intent(in) :: dt,flattening(0:,0:,0:)
    type(eb_plm_face_prediction_3d), intent(inout) :: prediction(3)
    ! Uncorrected normal predictors use the same cropped face indices as prediction.
    ! This aligns the two stages for conservative face correction budgets.
    type(eb_plm_face_prediction_3d), optional, intent(inout) :: normal_prediction(3)
    type(eb_plm_boundary_3d), optional, intent(in) :: boundary
    logical, optional, intent(in) :: capture_internal_energy_flux
    logical :: capture_held
    real(dp) :: held_flux
    logical, intent(out) :: ok
    integer, optional, intent(in) :: order
    real(dp), optional, intent(in) :: source(:,:,:,:),source_energy(:,:,:),conserved_source(:,:,:,:)
    real(dp), allocatable :: forcing(:,:,:,:),forcing_energy(:,:,:)
    type(eb_plm_face_prediction_3d) :: normal(3),corner(3,3),final(3),candidate(3),normal_candidate(3)
    real(dp), allocatable :: gamma(:,:,:),sound(:,:,:)
    real(dp) :: q(reactive_nvar(size(species))),recovered,speed,spacing(3)
    integer :: n(3),nv,i,j,k,axis,t,lo(3),hi(3),other(2),m
    logical :: valid
    capture_held=.false.
    if(present(capture_internal_energy_flux)) capture_held=capture_internal_energy_flux
    ok=.false.;n=[geometry%nx,geometry%ny,geometry%nz];nv=reactive_nvar(size(species))
    if(any(n<5)) return
    if(.not.valid_plm_boundary(boundary)) return
    if(present(source).neqv.present(source_energy)) return
    if(present(conserved_source)) then
      if(present(source).or.present(source_energy)) return
      if(any(shape(conserved_source)/=[nv,n])) return
    end if
    if(present(source)) then
      if(any(shape(source)/=[nv,n]).or.any(shape(source_energy)/=n)) return
    end if
    call compute_reactive_eb_plm_normal_fluxes_3d(species,geometry,ghost,temperature,volume_fraction, &
      dt,flattening,normal,valid,order,boundary,capture_internal_energy_flux)
    if(.not.valid) return
    spacing=[geometry%dx,geometry%dy,geometry%dz]
    allocate(gamma(n(1),n(2),n(3)),sound(n(1),n(2),n(3)));gamma=0;sound=0
    if(present(source)) then
      forcing=source;forcing_energy=source_energy
    else if(present(conserved_source)) then
      allocate(forcing(nv,n(1),n(2),n(3)),forcing_energy(n(1),n(2),n(3)))
      forcing=0;forcing_energy=0
    end if
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      if(geometry%volume_fraction(i,j,k)<=0) cycle
      call reactive_conserved_to_primitive(species,ghost(:,i,j,k),temperature(i,j,k),q,recovered,speed,valid)
      if(.not.valid) return
      gamma(i,j,k)=speed**2*q(1)/q(5);sound(i,j,k)=speed
      if(present(conserved_source)) then
        call reactive_pelec_source_to_primitive(species,q,conserved_source(:,i,j,k), &
          forcing(:,i,j,k),forcing_energy(i,j,k),valid)
        if(.not.valid) return
      end if
    end do
    end do
    end do
    do axis=1,3
      do t=1,3
        if(t==axis) cycle
        corner(axis,t)=normal(axis)
        lo=2;lo(axis)=1;hi=n-1
        do k=lo(3),hi(3)
        do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          call correct_face(axis,[t],[i,j,k],corner(axis,t),valid)
          if(.not.valid) return
        end do
        end do
        end do
      end do
    end do
    do axis=1,3
      m=0
      do t=1,3
        if(t==axis) cycle
        m=m+1;other(m)=t
      end do
      final(axis)=normal(axis)
      lo=3;lo(axis)=2;hi=n-2
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        call correct_face(axis,other,[i,j,k],final(axis),valid)
        if(.not.valid) return
      end do
      end do
      end do
      call crop_prediction(final(axis),lo,hi,candidate(axis))
      if(present(normal_prediction)) call crop_prediction(normal(axis),lo,hi,normal_candidate(axis))
    end do
    if(present(normal_prediction)) normal_prediction=normal_candidate
    prediction=candidate;ok=.true.
  contains
    ! Shift the normal and final stages together, including auxiliary Riemann
    ! states used by viscosity and centroid interpolation. Publish only on success.
    subroutine crop_prediction(f,lower,upper,c)
      type(eb_plm_face_prediction_3d), intent(in) :: f
      integer, intent(in) :: lower(3),upper(3)
      type(eb_plm_face_prediction_3d), intent(out) :: c
      integer :: idx(3)
      idx=lower-2
      allocate(c%left(nv,idx(1):upper(1)-2,idx(2):upper(2)-2,idx(3):upper(3)-2), &
        c%right(nv,idx(1):upper(1)-2,idx(2):upper(2)-2,idx(3):upper(3)-2), &
        c%flux(nv,idx(1):upper(1)-2,idx(2):upper(2)-2,idx(3):upper(3)-2), &
        c%energy_left(idx(1):upper(1)-2,idx(2):upper(2)-2,idx(3):upper(3)-2), &
        c%energy_right(idx(1):upper(1)-2,idx(2):upper(2)-2,idx(3):upper(3)-2), &
        c%pressure(idx(1):upper(1)-2,idx(2):upper(2)-2,idx(3):upper(3)-2), &
        c%velocity(idx(1):upper(1)-2,idx(2):upper(2)-2,idx(3):upper(3)-2))
      c%left=f%left(:,lower(1):upper(1),lower(2):upper(2),lower(3):upper(3))
      c%right=f%right(:,lower(1):upper(1),lower(2):upper(2),lower(3):upper(3))
      c%flux=f%flux(:,lower(1):upper(1),lower(2):upper(2),lower(3):upper(3))
      c%energy_left=f%energy_left(lower(1):upper(1),lower(2):upper(2),lower(3):upper(3))
      c%energy_right=f%energy_right(lower(1):upper(1),lower(2):upper(2),lower(3):upper(3))
      c%pressure=f%pressure(lower(1):upper(1),lower(2):upper(2),lower(3):upper(3))
      c%velocity=f%velocity(lower(1):upper(1),lower(2):upper(2),lower(3):upper(3))
      if(allocated(f%internal_energy_flux)) then
        allocate(c%internal_energy_flux(idx(1):upper(1)-2,idx(2):upper(2)-2,idx(3):upper(3)-2))
        c%internal_energy_flux=f%internal_energy_flux(lower(1):upper(1),lower(2):upper(2),lower(3):upper(3))
      endif
    end subroutine crop_prediction

    real(dp) function area(direction,index) result(value)
      integer, intent(in) :: direction,index(3)
      select case(direction)
      case(1);value=geometry%x_face_fraction(index(1),index(2),index(3))
      case(2);value=geometry%y_face_fraction(index(1),index(2),index(3))
      case(3);value=geometry%z_face_fraction(index(1),index(2),index(3))
      end select
    end function area
    subroutine correct_face(direction,transverse,index,result,success)
      integer, intent(in) :: direction,transverse(:),index(3)
      type(eb_plm_face_prediction_3d), intent(inout) :: result
      logical, intent(out) :: success
      real(dp) :: states(nv,2),energies(2),ql(nv),qr(nv),f(nv),el,er
      real(dp) :: tf(nv,2,size(transverse)),tp(2,size(transverse)),tv(2,size(transverse))
      real(dp) :: ta(2,size(transverse)),coef(size(transverse)),na(2),ap,swap,pface,uface
      integer :: side,d,face,cell(3),point(3),nt,normal_side,second
      success=.false.;nt=size(transverse);ap=area(direction,index)
      if(ap<=0) then
        success=.true.;return
      end if
      do side=0,1
        cell=index;cell(direction)=cell(direction)+side
        states=0;energies=0;na=0;normal_side=2-side;na(normal_side)=ap
        if(side==0) then
          states(:,2)=normal(direction)%left(:,index(1),index(2),index(3))
          energies(2)=normal(direction)%energy_left(index(1),index(2),index(3))
        else
          states(:,1)=normal(direction)%right(:,index(1),index(2),index(3))
          energies(1)=normal(direction)%energy_right(index(1),index(2),index(3))
        end if
        do d=1,nt
          coef(d)=dt/(3*spacing(transverse(d)))
          if(nt==2) coef(d)=.5_dp*dt/spacing(transverse(d))
          do face=1,2
            point=cell;point(transverse(d))=point(transverse(d))+face-2
            ta(face,d)=area(transverse(d),point)
            if(nt==1) then
              tf(:,face,d)=normal(transverse(d))%flux(:,point(1),point(2),point(3))
              tp(face,d)=normal(transverse(d))%pressure(point(1),point(2),point(3))
              tv(face,d)=normal(transverse(d))%velocity(point(1),point(2),point(3))
            else
              second=transverse(3-d)
              tf(:,face,d)=corner(transverse(d),second)%flux(:,point(1),point(2),point(3))
              tp(face,d)=corner(transverse(d),second)%pressure(point(1),point(2),point(3))
              tv(face,d)=corner(transverse(d),second)%velocity(point(1),point(2),point(3))
            end if
          end do
        end do
        if(nt==2.and.allocated(forcing)) then
          call apply_reactive_eb_transverse_3d(states,energies,gamma(cell(1),cell(2),cell(3)),na,ta,tf,tp,tv, &
            coef,success,forcing(:,cell(1),cell(2),cell(3)),forcing_energy(cell(1),cell(2),cell(3)),.5_dp*dt)
        else
          call apply_reactive_eb_transverse_3d(states,energies,gamma(cell(1),cell(2),cell(3)),na,ta,tf,tp,tv,coef,success)
        end if
        if(.not.success) return
        if(side==0) then
          ql=states(:,2);el=energies(2)
        else
          qr=states(:,1);er=energies(1)
        end if
      end do
      result%left(:,index(1),index(2),index(3))=ql;result%right(:,index(1),index(2),index(3))=qr
      result%energy_left(index(1),index(2),index(3))=el;result%energy_right(index(1),index(2),index(3))=er
      call correct_plm_boundary(ql,qr,direction,index(direction),boundary)
      swap=ql(2);ql(2)=ql(direction+1);ql(direction+1)=swap
      swap=qr(2);qr(2)=qr(direction+1);qr(direction+1)=swap
      cell=index;cell(direction)=cell(direction)+1
      if(capture_held) then
      call reactive_pelec_primitive_flux_x(species,ql,qr, &
        .5_dp*(sound(index(1),index(2),index(3))+sound(cell(1),cell(2),cell(3))), &
        f,success,interface_pressure=pface,interface_velocity=uface,passive_upwind=.true.,internal_energy_flux=held_flux)
      else
      call reactive_pelec_primitive_flux_x(species,ql,qr, &
        .5_dp*(sound(index(1),index(2),index(3))+sound(cell(1),cell(2),cell(3))), &
        f,success,interface_pressure=pface,interface_velocity=uface,passive_upwind=.true.)
      endif
      if(.not.success) return
      swap=f(2);f(2)=f(direction+1);f(direction+1)=swap
      result%flux(:,index(1),index(2),index(3))=f
      result%pressure(index(1),index(2),index(3))=pface;result%velocity(index(1),index(2),index(3))=uface
      if(capture_held) result%internal_energy_flux(index(1),index(2),index(3))=held_flux
    end subroutine correct_face
  end subroutine compute_reactive_eb_plm_ctu_fluxes_3d

  ! Pinned pc_transdo/pc_transdd for a cell's two normal face states.
  ! Fluxes are intensive and in global momentum order. For two transverse
  ! directions, all four faces must be open, as in PeleC's joint EB gate.
  subroutine apply_reactive_eb_transverse_3d(states,internal_energy,gamma,normal_area, &
      transverse_area,flux,face_pressure,face_velocity,coefficient,ok,source,source_energy,source_dt)
    real(dp), intent(inout) :: states(:,:),internal_energy(:)
    real(dp), intent(in) :: gamma,normal_area(2),transverse_area(:,:),flux(:,:,:)
    real(dp), intent(in) :: face_pressure(:,:),face_velocity(:,:),coefficient(:)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: source(:),source_energy,source_dt
    real(dp) :: candidate(size(states,1),2),energy(2),delta(size(states,1)),pressure_work(2)
    real(dp) :: forcing(size(states,1)),forcing_energy,hdt,rho,new_rho,momentum(3),total,pnew
    integer :: nv,nt,d,side,supplied
    logical :: use_flux
    ok=.false.;nv=size(states,1);nt=size(coefficient)
    if(nv<6.or.size(states,2)/=2.or.size(internal_energy)/=2.or.nt<1.or.nt>2) return
    if(any(shape(flux)/=[nv,2,nt]).or.any(shape(transverse_area)/=[2,nt])) return
    if(any(shape(face_pressure)/=[2,nt]).or.any(shape(face_velocity)/=[2,nt])) return
    if(.not.all(ieee_is_finite([gamma,normal_area,coefficient]))) return
    if(.not.all(ieee_is_finite(transverse_area))) return
    if(gamma<=0.or.any(normal_area<0).or.any(normal_area>1).or.any(coefficient<0)) return
    if(any(transverse_area<0).or.any(transverse_area>1)) return
    forcing=0;forcing_energy=0;hdt=0
    supplied=count([present(source),present(source_energy),present(source_dt)])
    if(supplied/=0.and.supplied/=3) return
    if(supplied==3) then
      if(nt/=2.or.size(source)/=nv) return
      if(.not.all(ieee_is_finite(source)).or..not.all(ieee_is_finite([source_energy,source_dt]))) return
      if(source_dt<0) return
      forcing=source;forcing_energy=source_energy;hdt=source_dt
    end if
    use_flux=all(transverse_area>0)
    delta=0;pressure_work=0
    if(use_flux) then
      if(.not.all(ieee_is_finite(flux)).or..not.all(ieee_is_finite(face_pressure)).or. &
        .not.all(ieee_is_finite(face_velocity))) return
      do d=1,nt
        delta=delta+coefficient(d)*(flux(:,2,d)-flux(:,1,d))
        pressure_work(d)=coefficient(d)*(face_pressure(2,d)*face_velocity(2,d)- &
          face_pressure(1,d)*face_velocity(1,d)+.5_dp*(face_pressure(2,d)+face_pressure(1,d))* &
          (face_velocity(2,d)-face_velocity(1,d))*(gamma-1))
      end do
    end if
    candidate=states;energy=internal_energy
    do side=1,2
      if(normal_area(side)<=0) cycle
      if(.not.all(ieee_is_finite(states(:,side))).or..not.ieee_is_finite(internal_energy(side))) return
      rho=states(1,side);new_rho=rho-delta(1)
      if(rho<=0.or.new_rho<=0) return
      momentum=rho*states(2:4,side)-delta(2:4)
      total=internal_energy(side)+.5_dp*rho*sum(states(2:4,side)**2)-delta(5)
      candidate(1,side)=new_rho+hdt*forcing(1)
      candidate(2:4,side)=momentum/new_rho+hdt*forcing(2:4)
      energy(side)=total-.5_dp*sum(momentum**2)/new_rho+hdt*forcing_energy
      ! Preserve the reference's ordered subtraction for the two-direction stage.
      pnew=states(5,side)-pressure_work(1)-pressure_work(2)+hdt*forcing(5)
      candidate(5,side)=max(tiny(1._dp),pnew)
      candidate(6:,side)=(rho*states(6:,side)-delta(6:))/new_rho+hdt*forcing(6:)
      if(.not.all(ieee_is_finite(candidate(:,side))).or..not.ieee_is_finite(energy(side))) return
      ! Reference transverse passives may be signed before the Riemann solve.
      if(candidate(1,side)<=0.or.sum(candidate(6:,side))<=0) return
    end do
    states=candidate;internal_energy=energy;ok=.true.
  end subroutine apply_reactive_eb_transverse_3d

  ! First Godunov stage only: intensive normal fluxes before transverse,
  ! viscosity, centroid interpolation or aperture/area integration.
  subroutine compute_reactive_eb_plm_normal_fluxes_3d(species,geometry,ghost,temperature,volume_fraction, &
      dt,flattening,prediction,ok,order,boundary,capture_internal_energy_flux)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: ghost(:,-2:,-2:,-2:),temperature(-2:,-2:,-2:),volume_fraction(-2:,-2:,-2:)
    real(dp), intent(in) :: dt,flattening(0:,0:,0:)
    type(eb_plm_face_prediction_3d), intent(inout) :: prediction(3)
    type(eb_plm_boundary_3d), optional, intent(in) :: boundary
    logical, optional, intent(in) :: capture_internal_energy_flux
    logical :: capture_held
    real(dp) :: held_flux
    logical, intent(out) :: ok
    integer, optional, intent(in) :: order
    type(eb_plm_face_prediction_3d) :: candidate(3)
    real(dp), allocatable :: primitive(:,:,:,:),sound(:,:,:),energy(:,:,:)
    real(dp) :: stencil(reactive_nvar(size(species)),-2:2),epencil(-2:2)
    real(dp) :: ql(reactive_nvar(size(species))),qr(reactive_nvar(size(species))),unused(reactive_nvar(size(species)))
    real(dp) :: f(reactive_nvar(size(species))),el,er,ignored,aperture,spacing(3),recovered,speed,swap,pface,uface
    integer :: n(3),nv,axis,i,j,k,idx(3),lo(3),hi(3),cell(3),point(3),offset,side,selected
    logical :: covered(-2:2),valid
    capture_held=.false.
    if(present(capture_internal_energy_flux)) capture_held=capture_internal_energy_flux
    ok=.false.;nv=reactive_nvar(size(species));selected=4
    if(present(order)) selected=order
    if(.not.valid_plm_boundary(boundary)) return
    if(.not.geometry%is_valid().or.(selected/=1.and.selected/=2.and.selected/=4)) return
    n=[geometry%nx,geometry%ny,geometry%nz];spacing=[geometry%dx,geometry%dy,geometry%dz]
    if(any(shape(ghost)/=[nv,n+6]).or.any(shape(temperature)/=n+6).or.any(shape(volume_fraction)/=n+6)) return
    if(any(shape(flattening)/=n+2)) return
    if(.not.ieee_is_finite(dt).or.dt<0) return
    if(.not.all(ieee_is_finite(flattening)).or..not.all(ieee_is_finite(volume_fraction))) return
    if(any(flattening<0).or.any(flattening>1).or.any(volume_fraction<0).or.any(volume_fraction>1)) return
    if(any(volume_fraction(1:n(1),1:n(2),1:n(3))/=geometry%volume_fraction)) return
    allocate(primitive(nv,-2:n(1)+3,-2:n(2)+3,-2:n(3)+3), &
      sound(-2:n(1)+3,-2:n(2)+3,-2:n(3)+3),energy(-2:n(1)+3,-2:n(2)+3,-2:n(3)+3))
    primitive=0;sound=0;energy=0
    do k=-2,n(3)+3
    do j=-2,n(2)+3
    do i=-2,n(1)+3
      idx=[i,j,k]
      ! Only face-normal pencils enter this stage; unused edge/corner data may be poisoned.
      if(count(idx<1.or.idx>n)>1.or.volume_fraction(i,j,k)==0) cycle
      call reactive_conserved_to_primitive(species,ghost(:,i,j,k),temperature(i,j,k), &
        primitive(:,i,j,k),recovered,speed,valid)
      if(.not.valid) return
      sound(i,j,k)=speed
      energy(i,j,k)=ghost(5,i,j,k)-.5_dp*sum(ghost(2:4,i,j,k)*primitive(2:4,i,j,k))
    end do
    end do
    end do
    do axis=1,3
      lo=1;lo(axis)=0;hi=n
      associate(result=>candidate(axis))
      allocate(result%left(nv,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
        result%right(nv,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
        result%flux(nv,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
        result%energy_left(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
        result%energy_right(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
        result%pressure(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
        result%velocity(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      result%left=0;result%right=0;result%flux=0;result%energy_left=0;result%energy_right=0
      result%pressure=0;result%velocity=0
      if(capture_held) then
        allocate(result%internal_energy_flux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),source=0._dp)
      endif
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        idx=[i,j,k]
        select case(axis)
        case(1);aperture=geometry%x_face_fraction(i,j,k)
        case(2);aperture=geometry%y_face_fraction(i,j,k)
        case(3);aperture=geometry%z_face_fraction(i,j,k)
        end select
        if(aperture<=0) cycle
        do side=0,1
          cell=idx;cell(axis)=cell(axis)+side
          if(volume_fraction(cell(1),cell(2),cell(3))<=0) return
          do offset=-2,2
            point=cell;point(axis)=point(axis)+offset
            covered(offset)=volume_fraction(point(1),point(2),point(3))<=0
            stencil(:,offset)=primitive(:,point(1),point(2),point(3))
            epencil(offset)=energy(point(1),point(2),point(3))
          end do
          if(side==0) then
            call predict_reactive_eb_plm_cell_3d(stencil,epencil,sound(cell(1),cell(2),cell(3)),axis,spacing(axis), &
              dt,flattening(cell(1),cell(2),cell(3)),covered,[0._dp,aperture],unused,ql,ignored,el,valid,selected)
          else
            call predict_reactive_eb_plm_cell_3d(stencil,epencil,sound(cell(1),cell(2),cell(3)),axis,spacing(axis), &
              dt,flattening(cell(1),cell(2),cell(3)),covered,[aperture,0._dp],qr,unused,er,ignored,valid,selected)
          end if
          if(.not.valid) return
        end do
        result%left(:,i,j,k)=ql;result%right(:,i,j,k)=qr
        result%energy_left(i,j,k)=el;result%energy_right(i,j,k)=er
        call correct_plm_boundary(ql,qr,axis,idx(axis),boundary)
        swap=ql(2);ql(2)=ql(axis+1);ql(axis+1)=swap
        swap=qr(2);qr(2)=qr(axis+1);qr(axis+1)=swap
        cell=idx;cell(axis)=cell(axis)+1
        if(capture_held) then
        call reactive_pelec_primitive_flux_x(species,ql,qr,.5_dp*(sound(i,j,k)+sound(cell(1),cell(2),cell(3))), &
          f,valid,interface_pressure=pface,interface_velocity=uface,passive_upwind=.true.,internal_energy_flux=held_flux)
        else
        call reactive_pelec_primitive_flux_x(species,ql,qr,.5_dp*(sound(i,j,k)+sound(cell(1),cell(2),cell(3))), &
          f,valid,interface_pressure=pface,interface_velocity=uface,passive_upwind=.true.)
        endif
        if(.not.valid) return
        swap=f(2);f(2)=f(axis+1);f(axis+1)=swap
        result%flux(:,i,j,k)=f;result%pressure(i,j,k)=pface;result%velocity(i,j,k)=uface
        if(capture_held) result%internal_energy_flux(i,j,k)=held_flux
      end do
      end do
      end do
      end associate
    end do
    prediction=candidate;ok=.true.
  end subroutine compute_reactive_eb_plm_normal_fluxes_3d

  ! Pinned pc_plm_d_eb normal prediction. Primitive layout is rho,u,v,w,p,Y;
  ! rhoe is kept separately for subsequent EB transverse pressure/energy work.
  subroutine predict_reactive_eb_plm_cell_3d(stencil,rhoe,sound,axis,dx,dt,flat,covered,area, &
      minus,plus,energy_minus,energy_plus,ok,order)
    real(dp), intent(in) :: stencil(:,-2:),rhoe(-2:),sound,dx,dt,flat,area(2)
    integer, intent(in) :: axis
    logical, intent(in) :: covered(-2:)
    real(dp), intent(out) :: minus(:),plus(:),energy_minus,energy_plus
    logical, intent(out) :: ok
    integer, optional, intent(in) :: order
    real(dp) :: rotated(size(stencil,1),-2:2),qm(size(stencil,1)),qp(size(stencil,1))
    real(dp) :: energy_stencil(1,-2:2),energy_slope(1),slope(size(stencil,1))
    real(dp) :: w(3),rho,c2,enth,am,ap,a0e,contact,acoustic,dtdx,em,ep
    logical :: opened(2),valid
    integer :: selected,nv
    minus=0;plus=0;energy_minus=0;energy_plus=0;ok=.false.
    nv=size(stencil,1);selected=4
    if(present(order)) selected=order
    if(nv<6.or.size(stencil,2)/=5.or.size(rhoe)/=5.or.size(covered)/=5) return
    if(size(minus)/=nv.or.size(plus)/=nv.or.axis<1.or.axis>3) return
    if(.not.all(ieee_is_finite([sound,dx,dt,flat,area]))) return
    if(sound<=0.or.dx<=0.or.dt<0.or.flat<0.or.flat>1.or.any(area<0).or.any(area>1)) return
    if(selected/=1.and.selected/=2.and.selected/=4) return
    opened=area>0
    if(.not.any(opened)) then
      ok=.true.;return
    end if
    if(covered(0).or..not.ieee_is_finite(rhoe(0))) return
    rotated=stencil
    rotated(2,:)=stencil(axis+1,:);rotated(axis+1,:)=stencil(2,:)
    dtdx=dt/dx
    ! pc_plm_d_eb_passive preserves independently limited face weights.
    call reconstruct_pelec_plm_cell(rotated,sound,dtdx,flat,qm,qp,valid,selected,covered,opened, &
      normalize_species=.false.)
    if(.not.valid) return
    call compute_pelec_plm_slopes(rotated,flat,selected,slope,valid,covered)
    if(.not.valid) return
    energy_stencil(1,:)=rhoe
    call compute_pelec_plm_slopes(energy_stencil,flat,selected,energy_slope,valid,covered)
    if(.not.valid) return
    rho=rotated(1,0);c2=sound**2
    enth=((rhoe(0)+rotated(5,0))/rho)/c2
    w=[rotated(2,0)-sound,rotated(2,0),rotated(2,0)+sound]
    am=.5_dp*(slope(5)/(rho*sound)-slope(2))*rho/sound
    ap=.5_dp*(slope(5)/(rho*sound)+slope(2))*rho/sound
    a0e=energy_slope(1)-slope(5)*enth
    em=0;ep=0
    if(opened(1)) then
      acoustic=.25_dp*dtdx*(w(1)-w(3))*(1-sign(1._dp,w(3)))*ap
      contact=.25_dp*dtdx*(w(1)-w(2))*(1-sign(1._dp,w(2)))*a0e
      em=rhoe(0)-.5_dp*(1+dtdx*min(w(1),0._dp))*energy_slope(1)+acoustic*enth*c2+contact
    end if
    if(opened(2)) then
      acoustic=.25_dp*dtdx*(w(3)-w(1))*(1+sign(1._dp,w(1)))*am
      contact=.25_dp*dtdx*(w(3)-w(2))*(1+sign(1._dp,w(2)))*a0e
      ep=rhoe(0)+.5_dp*(1-dtdx*max(w(3),0._dp))*energy_slope(1)+acoustic*enth*c2+contact
    end if
    if(.not.all(ieee_is_finite([em,ep]))) return
    minus=qm;plus=qp
    minus(2)=qm(axis+1);minus(axis+1)=qm(2)
    plus(2)=qp(axis+1);plus(axis+1)=qp(2)
    energy_minus=em;energy_plus=ep;ok=.true.
  end subroutine predict_reactive_eb_plm_cell_3d

  subroutine reactive_eb_outflow_riemann_fluxes_3d( &
      species, state, temperature, geometry, solver, x_flux, y_flux, &
      z_flux, ok, ghost_state, ghost_temperature, difmag, ghost_volume_fraction, no_slip_faces)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    character(len=*), intent(in) :: solver
    real(dp), intent(out) :: x_flux(:, 0:, :, :)
    real(dp), intent(out) :: y_flux(:, :, 0:, :)
    real(dp), intent(out) :: z_flux(:, :, :, 0:)
    logical, intent(out) :: ok
    ! Optional face halo, indexed 0:n+1 in each spatial direction. Only open
    ! boundary-face donors are read for explicit zero viscosity. Positive difmag
    ! also needs edge/corner state donors and a matching volume-fraction halo.
    ! With no halo, retain the historical extrapolated-outflow boundary states.
    real(dp), optional, intent(in) :: ghost_state(:,0:,0:,0:),ghost_temperature(0:,0:,0:)
    real(dp), optional, intent(in) :: difmag,ghost_volume_fraction(0:,0:,0:)
    logical, optional, intent(in) :: no_slip_faces(6)
    real(dp), allocatable :: halo(:,:,:,:),vf(:,:,:)
    real(dp) :: coefficient
    logical :: walls(6)
    integer :: n(3),q(3)
    real(dp), allocatable :: left_state(:),right_state(:)
    real(dp) :: left_temperature,right_temperature

    real(dp), allocatable :: candidate_x(:, :, :, :)
    real(dp), allocatable :: candidate_y(:, :, :, :)
    real(dp), allocatable :: candidate_z(:, :, :, :)
    logical :: face_ok
    integer :: i, j, k, lower_i, upper_i, lower_j, upper_j
    integer :: lower_k, upper_k, nvar

    x_flux = 0.0_dp
    y_flux = 0.0_dp
    z_flux = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (.not. valid_riemann_solver(solver)) return
    if (.not. valid_state_shapes(state, temperature, geometry, nvar) .or. &
        size(x_flux, 1) /= nvar .or. &
        size(x_flux, 2) /= geometry%nx + 1 .or. &
        size(x_flux, 3) /= geometry%ny .or. &
        size(x_flux, 4) /= geometry%nz .or. &
        size(y_flux, 1) /= nvar .or. &
        size(y_flux, 2) /= geometry%nx .or. &
        size(y_flux, 3) /= geometry%ny + 1 .or. &
        size(y_flux, 4) /= geometry%nz .or. &
        size(z_flux, 1) /= nvar .or. &
        size(z_flux, 2) /= geometry%nx .or. &
        size(z_flux, 3) /= geometry%ny .or. &
        size(z_flux, 4) /= geometry%nz + 1) return

    if (present(ghost_state).neqv.present(ghost_temperature)) return
    if (present(ghost_state)) then
      if (any(shape(ghost_state)/=[nvar,geometry%nx+2,geometry%ny+2,geometry%nz+2])) return
      if (any(shape(ghost_temperature)/=[geometry%nx+2,geometry%ny+2,geometry%nz+2])) return
    end if
    coefficient=.1_dp
    if(present(difmag)) coefficient=difmag
    if(.not.ieee_is_finite(coefficient).or.coefficient<0) return
    n=[geometry%nx,geometry%ny,geometry%nz]
    if(coefficient>0) then
      ! Positive viscosity needs a complete corner halo, not just face donors.
      if(present(ghost_state).neqv.present(ghost_volume_fraction)) return
      if(present(ghost_volume_fraction)) then
        if(any(shape(ghost_volume_fraction)/=n+2)) return
        if(any(ghost_volume_fraction(1:n(1),1:n(2),1:n(3))/=geometry%volume_fraction)) return
      end if
    end if
    allocate(left_state(nvar),right_state(nvar))

    allocate(candidate_x(nvar, 0:geometry%nx, geometry%ny, geometry%nz))
    allocate(candidate_y(nvar, geometry%nx, 0:geometry%ny, geometry%nz))
    allocate(candidate_z(nvar, geometry%nx, geometry%ny, 0:geometry%nz))
    candidate_x = 0.0_dp
    candidate_y = 0.0_dp
    candidate_z = 0.0_dp

    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 0, geometry%nx
          if (geometry%x_face_fraction(i, j, k) <= 0.0_dp) cycle
          lower_i = max(1, i)
          upper_i = min(geometry%nx, i + 1)
          if (geometry%cell_type(lower_i, j, k) == eb_covered_cell_3d &
              .or. geometry%cell_type(upper_i, j, k) == &
                eb_covered_cell_3d) return
          left_state=state(:,lower_i, j, k);right_state=state(:,upper_i, j, k)
          left_temperature=temperature(lower_i, j, k);right_temperature=temperature(upper_i, j, k)
          if (present(ghost_state)) then
            if (i==0) then
              left_state=ghost_state(:,0,j,k);left_temperature=ghost_temperature(0,j,k)
            end if
            if (i==geometry%nx) then
              right_state=ghost_state(:,geometry%nx+1,j,k);right_temperature=ghost_temperature(geometry%nx+1,j,k)
            end if
          end if
          call reactive_riemann_flux_x( &
            species, left_state, right_state, &
            left_temperature, right_temperature, &
            solver, candidate_x(:, i, j, k), face_ok)
          if (.not. face_ok) return
        end do
      end do
    end do

    do k = 1, geometry%nz
      do j = 0, geometry%ny
        lower_j = max(1, j)
        upper_j = min(geometry%ny, j + 1)
        do i = 1, geometry%nx
          if (geometry%y_face_fraction(i, j, k) <= 0.0_dp) cycle
          if (geometry%cell_type(i, lower_j, k) == eb_covered_cell_3d &
              .or. geometry%cell_type(i, upper_j, k) == &
                eb_covered_cell_3d) return
          left_state=state(:,i, lower_j, k);right_state=state(:,i, upper_j, k)
          left_temperature=temperature(i, lower_j, k);right_temperature=temperature(i, upper_j, k)
          if (present(ghost_state)) then
            if (j==0) then
              left_state=ghost_state(:,i,0,k);left_temperature=ghost_temperature(i,0,k)
            end if
            if (j==geometry%ny) then
              right_state=ghost_state(:,i,geometry%ny+1,k);right_temperature=ghost_temperature(i,geometry%ny+1,k)
            end if
          end if
          call reactive_riemann_flux_y( &
            species, left_state, right_state, &
            left_temperature, right_temperature, &
            solver, candidate_y(:, i, j, k), face_ok)
          if (.not. face_ok) return
        end do
      end do
    end do

    do k = 0, geometry%nz
      lower_k = max(1, k)
      upper_k = min(geometry%nz, k + 1)
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%z_face_fraction(i, j, k) <= 0.0_dp) cycle
          if (geometry%cell_type(i, j, lower_k) == eb_covered_cell_3d &
              .or. geometry%cell_type(i, j, upper_k) == &
                eb_covered_cell_3d) return
          left_state=state(:,i, j, lower_k);right_state=state(:,i, j, upper_k)
          left_temperature=temperature(i, j, lower_k);right_temperature=temperature(i, j, upper_k)
          if (present(ghost_state)) then
            if (k==0) then
              left_state=ghost_state(:,i,j,0);left_temperature=ghost_temperature(i,j,0)
            end if
            if (k==geometry%nz) then
              right_state=ghost_state(:,i,j,geometry%nz+1);right_temperature=ghost_temperature(i,j,geometry%nz+1)
            end if
          end if
          call reactive_riemann_flux_z( &
            species, left_state, right_state, &
            left_temperature, right_temperature, &
            solver, candidate_z(:, i, j, k), face_ok)
          if (.not. face_ok) return
        end do
      end do
    end do
    if (any(.not. ieee_is_finite(candidate_x)) .or. &
        any(.not. ieee_is_finite(candidate_y)) .or. &
        any(.not. ieee_is_finite(candidate_z))) return

    if(coefficient>0) then
      allocate(halo(nvar,0:n(1)+1,0:n(2)+1,0:n(3)+1),vf(0:n(1)+1,0:n(2)+1,0:n(3)+1))
      if(present(ghost_state)) then
        halo=ghost_state;vf=ghost_volume_fraction
      else
        do k=0,n(3)+1
        do j=0,n(2)+1
        do i=0,n(1)+1
          q=max(1,min(n,[i,j,k]))
          halo(:,i,j,k)=state(:,q(1),q(2),q(3))
          vf(i,j,k)=geometry%volume_fraction(q(1),q(2),q(3))
        end do
        end do
        end do
      end if
      halo(:,1:n(1),1:n(2),1:n(3))=state
      walls=.false.
      if(present(no_slip_faces)) walls=no_slip_faces
      call adjust_reactive_artificial_viscosity_3d(halo,n,[geometry%dx,geometry%dy,geometry%dz], &
        coefficient,walls,candidate_x,candidate_y,candidate_z,face_ok,vf, &
        geometry%x_face_fraction>0,geometry%y_face_fraction>0,geometry%z_face_fraction>0)
      if(.not.face_ok) return
    end if

    x_flux = candidate_x
    y_flux = candidate_y
    z_flux = candidate_z
    ok = .true.
  end subroutine reactive_eb_outflow_riemann_fluxes_3d

  subroutine build_reactive_eb_conservative_rhs_3d( &
      species, state, temperature, geometry, solver, rhs, ok, ghost_state, ghost_temperature, flux, &
      difmag, ghost_volume_fraction, no_slip_faces)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    character(len=*), intent(in) :: solver
    real(dp), intent(out) :: rhs(:, :, :, :)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: ghost_state(:,0:,0:,0:),ghost_temperature(0:,0:,0:)
    real(dp), optional, intent(in) :: difmag,ghost_volume_fraction(0:,0:,0:)
    logical, optional, intent(in) :: no_slip_faces(6)
    type(eb_integrated_flux_3d), optional, intent(inout) :: flux

    real(dp), allocatable :: x_flux(:, :, :, :), y_flux(:, :, :, :)
    real(dp), allocatable :: z_flux(:, :, :, :)
    logical :: local_ok
    integer :: nvar

    rhs = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (.not. valid_state_shapes(state, temperature, geometry, nvar) .or. &
        size(rhs, 1) /= nvar .or. &
        size(rhs, 2) /= geometry%nx .or. &
        size(rhs, 3) /= geometry%ny .or. &
        size(rhs, 4) /= geometry%nz) return

    allocate(x_flux(nvar, 0:geometry%nx, geometry%ny, geometry%nz))
    allocate(y_flux(nvar, geometry%nx, 0:geometry%ny, geometry%nz))
    allocate(z_flux(nvar, geometry%nx, geometry%ny, 0:geometry%nz))
    call reactive_eb_outflow_riemann_fluxes_3d( &
      species, state, temperature, geometry, solver, x_flux, y_flux, &
      z_flux, local_ok, ghost_state, ghost_temperature,difmag,ghost_volume_fraction,no_slip_faces)
    if (.not. local_ok) return
    call reactive_eb_flux_divergence_3d( &
      species, state, temperature, geometry, x_flux, y_flux, z_flux, &
      rhs, local_ok)
    if (.not. local_ok) then
      rhs = 0.0_dp
      return
    end if
    if (present(flux)) then
      flux%x=x_flux;flux%y=y_flux;flux%z=z_flux
    end if
    ok = .true.
  end subroutine build_reactive_eb_conservative_rhs_3d

  subroutine advance_reactive_eb_euler_3d( &
      species, state, temperature, geometry, solver, dt, new_state, &
      new_temperature, ok, ghost_state, ghost_temperature, integrated_flux, &
      difmag, ghost_volume_fraction, no_slip_faces)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    real(dp), intent(out) :: new_state(:, :, :, :)
    real(dp), intent(out) :: new_temperature(:, :, :)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: ghost_state(:,0:,0:,0:),ghost_temperature(0:,0:,0:)
    real(dp), optional, intent(in) :: difmag,ghost_volume_fraction(0:,0:,0:)
    logical, optional, intent(in) :: no_slip_faces(6)
    type(eb_integrated_flux_3d), optional, intent(inout) :: integrated_flux
    type(eb_integrated_flux_3d) :: flux

    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    real(dp), allocatable :: rhs(:, :, :, :)
    real(dp), allocatable :: primitive(:)
    real(dp) :: recovered_temperature, sound_speed
    logical :: local_ok
    integer :: i, j, k, nvar

    new_state = 0.0_dp
    new_temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (.not. valid_state_shapes(state, temperature, geometry, nvar) .or. &
        any(shape(new_state) /= shape(state)) .or. &
        any(shape(new_temperature) /= shape(temperature))) return
    new_state = state
    new_temperature = temperature
    if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp) return

    allocate(rhs(nvar, geometry%nx, geometry%ny, geometry%nz))
    if (present(integrated_flux)) then
      call build_reactive_eb_conservative_rhs_3d( &
        species,state,temperature,geometry,solver,rhs,local_ok,ghost_state,ghost_temperature,flux, &
        difmag,ghost_volume_fraction,no_slip_faces)
    else
      call build_reactive_eb_conservative_rhs_3d( &
        species,state,temperature,geometry,solver,rhs,local_ok,ghost_state,ghost_temperature, &
        difmag=difmag,ghost_volume_fraction=ghost_volume_fraction,no_slip_faces=no_slip_faces)
    end if
    if (.not. local_ok) return

    allocate(candidate_state, source=state)
    allocate(candidate_temperature, source=temperature)
    allocate(primitive(reactive_nprim(size(species))))
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          candidate_state(:, i, j, k) = &
            state(:, i, j, k) + dt * rhs(:, i, j, k)
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

    if (present(integrated_flux)) then
      flux%x=dt*flux%x;flux%y=dt*flux%y;flux%z=dt*flux%z
      if (.not.all(ieee_is_finite(flux%x)).or..not.all(ieee_is_finite(flux%y)).or. &
        .not.all(ieee_is_finite(flux%z))) return
      integrated_flux=flux
    end if
    new_state = candidate_state
    new_temperature = candidate_temperature
    ok = .true.
  end subroutine advance_reactive_eb_euler_3d

  subroutine advance_reactive_eb_redistributed_euler_3d( &
      species, state, temperature, geometry, solver, dt, new_state, &
      new_temperature, ok, boundary_change, wall_change, difmag)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    real(dp), intent(out) :: new_state(:, :, :, :)
    real(dp), intent(out) :: new_temperature(:, :, :)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: difmag
    real(dp), optional, intent(inout) :: boundary_change(:),wall_change(:)
    type(eb_integrated_flux_3d) :: flux
    real(dp) :: boundary(reactive_nvar(size(species))),wall(reactive_nvar(size(species)))

    real(dp), allocatable :: rhs(:, :, :, :)
    logical :: local_ok
    integer :: nvar

    new_state = 0.0_dp
    new_temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (.not. valid_state_shapes(state, temperature, geometry, nvar) .or. &
        any(shape(new_state) /= shape(state)) .or. &
        any(shape(new_temperature) /= shape(temperature))) return
    new_state = state
    new_temperature = temperature
    if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp) return

    if (present(boundary_change)) then
      if (size(boundary_change)/=nvar) return
    end if
    if (present(wall_change)) then
      if (size(wall_change)/=nvar) return
    end if
    allocate(rhs(nvar, geometry%nx, geometry%ny, geometry%nz))
    if (present(boundary_change).or.present(wall_change)) then
      call build_reactive_eb_conservative_rhs_3d( &
        species,state,temperature,geometry,solver,rhs,local_ok,flux=flux,difmag=difmag)
      if (.not.local_ok) return
      call hydro_inventory_changes(species,state,temperature,geometry,dt,flux,boundary,wall,local_ok)
      if (.not.local_ok) return
    else
      call build_reactive_eb_conservative_rhs_3d( &
        species,state,temperature,geometry,solver,rhs,local_ok,difmag=difmag)
      if (.not.local_ok) return
    end if

    call advance_reactive_eb_redistributed_3d( &
      species, state, temperature, geometry, rhs, dt, new_state, &
      new_temperature, local_ok)
    if (.not. local_ok) return
    if (present(boundary_change)) boundary_change=boundary
    if (present(wall_change)) wall_change=wall
    ok = .true.
  end subroutine advance_reactive_eb_redistributed_euler_3d

  subroutine advance_reactive_eb_state_redistributed_euler_3d( &
      species, state, temperature, geometry, solver, dt, new_state, &
      new_temperature, ok, target_volume_fraction, boundary_change, wall_change, difmag, &
      state_redist_max_order)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    real(dp), intent(out) :: new_state(:, :, :, :)
    real(dp), intent(out) :: new_temperature(:, :, :)
    integer, optional, intent(in) :: state_redist_max_order
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: difmag
    real(dp), optional, intent(inout) :: boundary_change(:),wall_change(:)
    type(eb_integrated_flux_3d) :: flux
    real(dp) :: boundary(reactive_nvar(size(species))),wall(reactive_nvar(size(species)))
    real(dp), intent(in), optional :: target_volume_fraction

    real(dp), allocatable :: rhs(:, :, :, :)
    real(dp) :: target
    logical :: local_ok
    integer :: nvar

    new_state = 0.0_dp
    new_temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (.not. valid_state_shapes(state, temperature, geometry, nvar) .or. &
        any(shape(new_state) /= shape(state)) .or. &
        any(shape(new_temperature) /= shape(temperature))) return
    new_state = state
    new_temperature = temperature
    if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp) return

    if (present(boundary_change)) then
      if (size(boundary_change)/=nvar) return
    end if
    if (present(wall_change)) then
      if (size(wall_change)/=nvar) return
    end if
    allocate(rhs(nvar, geometry%nx, geometry%ny, geometry%nz))
    if (present(boundary_change).or.present(wall_change)) then
      call build_reactive_eb_conservative_rhs_3d( &
        species,state,temperature,geometry,solver,rhs,local_ok,flux=flux,difmag=difmag)
      if (.not.local_ok) return
      call hydro_inventory_changes(species,state,temperature,geometry,dt,flux,boundary,wall,local_ok)
      if (.not.local_ok) return
    else
      call build_reactive_eb_conservative_rhs_3d( &
        species,state,temperature,geometry,solver,rhs,local_ok,difmag=difmag)
      if (.not.local_ok) return
    end if

    target = 0.5_dp
    if (present(target_volume_fraction)) target = target_volume_fraction
    call advance_reactive_eb_state_redistributed_3d( &
      species, state, temperature, geometry, rhs, dt, new_state, &
      new_temperature, local_ok, target, &
      max_order=state_redist_max_order)
    if (.not. local_ok) return
    if (present(boundary_change)) boundary_change=boundary
    if (present(wall_change)) wall_change=wall
    ok = .true.
  end subroutine advance_reactive_eb_state_redistributed_euler_3d

  subroutine hydro_inventory_changes(species,state,temperature,geometry,dt,flux,boundary,wall,ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),dt
    type(eb_geometry_3d), intent(in) :: geometry
    type(eb_integrated_flux_3d), intent(in) :: flux
    real(dp), intent(out) :: boundary(:),wall(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: source(:,:,:,:)
    integer :: c,n(3)
    n=[geometry%nx,geometry%ny,geometry%nz]
    allocate(source(size(state,1),n(1),n(2),n(3)))
    call reactive_eb_slip_wall_source_3d(species,state,temperature,geometry,source,ok)
    if (.not.ok) return
    ! These are the actual Riemann faces used in the conservative RHS.
    ! Redistribution transfers inventory internally, leaving this exchange unchanged.
    do c=1,size(boundary)
      boundary(c)=dt*(geometry%dy*geometry%dz*sum( &
        geometry%x_face_fraction(0,:,:)*flux%x(c,0,:,:)-geometry%x_face_fraction(n(1),:,:)*flux%x(c,n(1),:,:)) &
        +geometry%dx*geometry%dz*sum( &
        geometry%y_face_fraction(:,0,:)*flux%y(c,:,0,:)-geometry%y_face_fraction(:,n(2),:)*flux%y(c,:,n(2),:)) &
        +geometry%dx*geometry%dy*sum( &
        geometry%z_face_fraction(:,:,0)*flux%z(c,:,:,0)-geometry%z_face_fraction(:,:,n(3))*flux%z(c,:,:,n(3))))
      wall(c)=dt*geometry%dx*geometry%dy*geometry%dz*sum(source(c,:,:,:)*geometry%volume_fraction)
    end do
    ok=all(ieee_is_finite(boundary)).and.all(ieee_is_finite(wall))
  end subroutine hydro_inventory_changes

  pure logical function valid_state_shapes( &
      state, temperature, geometry, nvar) result(valid)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: nvar

    valid = size(state, 1) == nvar .and. &
      size(state, 2) == geometry%nx .and. &
      size(state, 3) == geometry%ny .and. &
      size(state, 4) == geometry%nz .and. &
      all(shape(temperature) == &
        [geometry%nx, geometry%ny, geometry%nz])
  end function valid_state_shapes

  pure logical function valid_riemann_solver(solver) result(valid)
    character(len=*), intent(in) :: solver

    select case (trim(solver))
    case ("rusanov", "hllc", "pelec")
      valid = .true.
    case default
      valid = .false.
    end select
  end function valid_riemann_solver

end module eb_reactive_hydro_3d_mod
