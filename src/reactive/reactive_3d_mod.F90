module reactive_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d, fill_reactive_boundary_ghosts_3d, &
    fill_pressure_outlet_ghosts_3d,fill_mach_inlet_boundary_ghosts_3d
  use precision_mod, only: dp
  use constants_mod, only: density_floor, pressure_floor
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved, reactive_conserved_to_primitive, reconstructed_pelec_state, &
    reactive_riemann_flux_x, reactive_pelec_primitive_flux_x, characteristic_limited_slope, &
    reactive_ppm_flattening_coefficient, reconstruct_pelec_plm_cell
  use reactive_2d_mod, only: &
    reactive_riemann_flux_y, advance_reactive_chemistry_2d, reconstruct_reactive_characteristic_ppm_cell_2d, &
    apply_reactive_transverse_correction_2d
  use reactive_directional_flux_3d_mod, only: reactive_riemann_flux_z
  use reactive_transport_3d_mod, only: advance_reactive_transport_3d
  implicit none
  private

  public :: adjust_reactive_artificial_viscosity_3d
  public :: reactive_ppm_transverse_pressure_3d
  public :: compute_reactive_ppm_halo_ctu_fluxes_3d
  public :: advance_reactive_euler_ppm_ctu_3d
  public :: compute_reactive_ppm_ctu_fluxes_3d
  public :: compute_reactive_ppm_normal_fluxes_3d
  public :: advance_reactive_boundary_hydro_3d
  public :: compute_reactive_cfl_timestep_3d
  public :: compute_reactive_euler_rhs_3d
  public :: advance_reactive_euler_ssprk2_3d
  public :: compute_reactive_face_fluxes_3d
  public :: compute_reactive_plm_face_fluxes_3d
  public :: compute_reactive_plm_slab_face_fluxes_3d
  public :: advance_reactive_euler_ssprk2_with_fluxes_3d
  public :: advance_reactive_euler_ssprk2_plm_3d
  public :: advance_reactive_euler_ssprk2_plm_with_fluxes_3d
  public :: advance_reactive_chemistry_3d
  public :: advance_reactive_strang_3d
  public :: advance_reactive_full_3d
  public :: recover_reactive_temperatures_3d
  public :: reactive_integrals_3d
  public :: reactive_extrema_3d

contains

  ! Regular-cell pc_divu/pc_artif_visc/pc_norm_spec_flx. One ghost cell per side.
  ! Optional EB volume fractions use eb_divu's eight-cell regularity cutoff.
  ! Positive coefficients zero closed faces; zero retains the identity policy.
  ! Fluxes remain intensive; the finite-volume caller supplies face areas.
  subroutine adjust_reactive_artificial_viscosity_3d(state,n,spacing,difmag, &
      no_slip_faces,xflux,yflux,zflux,ok,volume_fraction,open_x,open_y,open_z,boundary_lower,boundary_upper)
    integer, intent(in) :: n(3)
    ! Physical face indices may lie inside an extended working region.
    integer, optional, intent(in) :: boundary_lower(3),boundary_upper(3)
    real(dp), intent(in) :: state(:,0:,0:,0:),spacing(3),difmag
    logical, intent(in) :: no_slip_faces(6)
    real(dp), intent(inout) :: xflux(:,0:,:,:),yflux(:,:,0:,:),zflux(:,:,:,0:)
    logical, intent(out) :: ok
    ! EB inputs are an all-or-none group; covered states are never read.
    real(dp), optional, intent(in) :: volume_fraction(0:,0:,0:)
    logical, optional, intent(in) :: open_x(0:,:,:),open_y(:,0:,:),open_z(:,:,0:)
    real(dp), allocatable :: velocity(:,:,:,:),divergence(:,:,:),fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:)
    real(dp) :: derivatives(3),factor,total,candidate(size(state,1)),jump(size(state,1))
    integer :: nv,i,j,k,a,b,c,s,t,point(3),low(3),high(3),limit(3),normal,wall_lo(3),wall_hi(3)
    logical :: wall,eb,opened
    ok=.false.;nv=size(state,1)
    if (nv<6.or.any(n<1).or.any(shape(state)/=[nv,n+2])) return
    if (any(shape(xflux)/=[nv,n(1)+1,n(2),n(3)]).or. &
        any(shape(yflux)/=[nv,n(1),n(2)+1,n(3)]).or.any(shape(zflux)/=[nv,n(1),n(2),n(3)+1])) return
    if (.not.all(ieee_is_finite([spacing,difmag])).or.any(spacing<=0).or.difmag<0) return
    if(present(boundary_lower).neqv.present(boundary_upper)) return
    wall_lo=0;wall_hi=n
    if(present(boundary_lower)) then
      wall_lo=boundary_lower;wall_hi=boundary_upper
      if(any(wall_lo>=wall_hi)) return
    end if
    eb=present(volume_fraction)
    if(eb.neqv.present(open_x)) return
    if(eb.neqv.present(open_y)) return
    if(eb.neqv.present(open_z)) return
    if(eb) then
      if(any(shape(volume_fraction)/=n+2)) return
      if(any(shape(open_x)/=[n(1)+1,n(2),n(3)]).or. &
        any(shape(open_y)/=[n(1),n(2)+1,n(3)]).or.any(shape(open_z)/=[n(1),n(2),n(3)+1])) return
      if(.not.all(ieee_is_finite(volume_fraction))) return
      if(any(volume_fraction<0).or.any(volume_fraction>1)) return
      do k=0,n(3)+1
      do j=0,n(2)+1
      do i=0,n(1)+1
        if(volume_fraction(i,j,k)==0) cycle
        if(.not.all(ieee_is_finite(state(:,i,j,k))).or.state(1,i,j,k)<=density_floor) return
      end do
      end do
      end do
    else
      if (.not.all(ieee_is_finite(state)).or.any(state(1,:,:,:)<=density_floor)) return
    end if
    if (.not.all(ieee_is_finite(xflux)).or..not.all(ieee_is_finite(yflux)).or. &
        .not.all(ieee_is_finite(zflux))) return
    if (difmag==0) then
      ok=.true.;return
    end if
    allocate(velocity(3,0:n(1)+1,0:n(2)+1,0:n(3)+1),divergence(0:n(1),0:n(2),0:n(3)))
    velocity=0
    do k=0,n(3)+1
    do j=0,n(2)+1
    do i=0,n(1)+1
      if(eb) then
        if(volume_fraction(i,j,k)==0) cycle
      end if
      velocity(:,i,j,k)=state(2:4,i,j,k)/state(1,i,j,k)
    end do
    end do
    end do
    do k=0,n(3)
    do j=0,n(2)
    do i=0,n(1)
      divergence(i,j,k)=0
      if(eb) then
        if(sum(volume_fraction(i:i+1,j:j+1,k:k+1))<8._dp*(1._dp-1e-12_dp)) cycle
      end if
      do a=1,3
        b=merge(2,1,a==1);c=merge(2,3,a==3)
        derivatives(a)=0
        do t=0,1
        do s=0,1
          low=[i,j,k];low(b)=low(b)+s;low(c)=low(c)+t
          high=low;high(a)=high(a)+1
          derivatives(a)=derivatives(a)+velocity(a,high(1),high(2),high(3))-velocity(a,low(1),low(2),low(3))
        end do
        end do
        derivatives(a)=.25_dp*derivatives(a)/spacing(a)
      end do
      divergence(i,j,k)=sum(derivatives)
    end do
    end do
    end do
    if (.not.all(ieee_is_finite(divergence))) return
    allocate(fx(nv,0:n(1),n(2),n(3)),fy(nv,n(1),0:n(2),n(3)),fz(nv,n(1),n(2),0:n(3)))
    fx=xflux;fy=yflux;fz=zflux
    do a=1,3
      b=merge(2,1,a==1);c=merge(2,3,a==3)
      limit=n;limit(a)=n(a)+1
      do k=1,limit(3)
      do j=1,limit(2)
      do i=1,limit(1)
        low=[i,j,k];low(a)=low(a)-1;high=low;high(a)=high(a)+1
        if(eb) then
          select case(a)
          case(1);opened=open_x(low(1),low(2),low(3))
          case(2);opened=open_y(low(1),low(2),low(3))
          case(3);opened=open_z(low(1),low(2),low(3))
          end select
          if(.not.opened) then
            select case(a)
            case(1);fx(:,low(1),low(2),low(3))=0
            case(2);fy(:,low(1),low(2),low(3))=0
            case(3);fz(:,low(1),low(2),low(3))=0
            end select
            cycle
          end if
          if(volume_fraction(low(1),low(2),low(3))==0.or.volume_fraction(high(1),high(2),high(3))==0) return
        end if
        point=low;point(b)=point(b)-1;point(c)=point(c)-1
        total=0
        do t=0,1
        do s=0,1
          high=point;high(b)=high(b)+s;high(c)=high(c)+t
          total=total+divergence(high(1),high(2),high(3))
        end do
        end do
        factor=spacing(a)*difmag*min(0._dp,.25_dp*total)
        high=low;high(a)=high(a)+1
        jump=factor*(state(:,high(1),high(2),high(3))-state(:,low(1),low(2),low(3)))
        normal=low(a)
        wall=(normal==wall_lo(a).and.no_slip_faces(2*a-1)).or.(normal==wall_hi(a).and.no_slip_faces(2*a))
        if (wall) then
          jump(b+1)=0;jump(c+1)=0
        end if
        select case(a)
        case(1);candidate=fx(:,low(1),low(2),low(3))+jump
        case(2);candidate=fy(:,low(1),low(2),low(3))+jump
        case(3);candidate=fz(:,low(1),low(2),low(3))+jump
        end select
        total=sum(candidate(6:))
        if (total/=0) candidate(6:)=candidate(6:)*(candidate(1)/total)
        if (.not.all(ieee_is_finite(candidate))) return
        select case(a)
        case(1);fx(:,low(1),low(2),low(3))=candidate
        case(2);fy(:,low(1),low(2),low(3))=candidate
        case(3);fz(:,low(1),low(2),low(3))=candidate
        end select
      end do
      end do
      end do
    end do
    xflux=fx;yflux=fy;zflux=fz;ok=.true.
  end subroutine adjust_reactive_artificial_viscosity_3d


  ! Regular-cell pc_transdo/pc_transdd pressure work (no forcing sources).
  ! Conserved flux differences and primitive pressure are evolved separately.
  subroutine reactive_ppm_transverse_pressure_3d(normal,normal_pressure,delta,gamma, &
      lower_pressure,upper_pressure,lower_velocity,upper_velocity,scale,primitive,internal_energy,ok)
    real(dp), intent(in) :: normal(:),normal_pressure,delta(:),gamma
    real(dp), intent(in) :: lower_pressure(:),upper_pressure(:),lower_velocity(:),upper_velocity(:),scale(:)
    real(dp), intent(out) :: primitive(:),internal_energy
    logical, intent(out) :: ok
    real(dp) :: candidate(size(normal)),state(size(normal)),pressure,work,pav,du,floor_pressure
    integer :: axis,nd
    primitive=0;internal_energy=0;ok=.false.;nd=size(scale)
    if (nd<1.or.nd>2.or.size(normal)<6.or.size(delta)/=size(normal).or.size(primitive)/=size(normal)) return
    if (size(lower_pressure)/=nd.or.size(upper_pressure)/=nd.or. &
        size(lower_velocity)/=nd.or.size(upper_velocity)/=nd) return
    if (.not.all(ieee_is_finite(normal)).or..not.all(ieee_is_finite(delta))) return
    if (.not.all(ieee_is_finite([normal_pressure,gamma,lower_pressure,upper_pressure, &
        lower_velocity,upper_velocity,scale]))) return
    if (normal(1)<=density_floor.or.normal_pressure<=0.or.gamma<=0.or.any(scale<0).or. &
        any(lower_pressure<=0).or.any(upper_pressure<=0)) return
    state=normal-delta
    if (state(1)<=density_floor) return
    ! pc_transdo_passive/pc_transdd_passive retain signed predictor weights.
    ! They are consumed by the primitive Riemann path; accepted cell species
    ! remain subject to the ordinary conserved-state/EOS validation.
    candidate(1)=state(1);candidate(2:4)=state(2:4)/state(1);candidate(6:)=state(6:)/state(1)
    pressure=normal_pressure
    do axis=1,nd
      work=upper_pressure(axis)*upper_velocity(axis)-lower_pressure(axis)*lower_velocity(axis)
      pav=.5_dp*(upper_pressure(axis)+lower_pressure(axis))
      du=upper_velocity(axis)-lower_velocity(axis)
      pressure=pressure-scale(axis)*(work+pav*du*(gamma-1))
    end do
    floor_pressure=tiny(1._dp)
    if (nd==1) floor_pressure=epsilon(1._dp)*1e-100_dp
    candidate(5)=max(pressure,floor_pressure)
    work=state(5)-.5_dp*sum(state(2:4)**2)/state(1)
    if (.not.all(ieee_is_finite(candidate)).or..not.ieee_is_finite(work)) return
    primitive=candidate;internal_energy=work;ok=.true.
  end subroutine reactive_ppm_transverse_pressure_3d

  ! Periodic time-centered normal PPM fluxes. These are predictors for a CTU
  ! update; transverse/corner corrections must be applied by the caller.
  subroutine compute_reactive_ppm_normal_fluxes_3d(species,state,temperature,spacing,dt, &
      steepening,flattening,solver,xflux,yflux,zflux,ok,normal_minus,normal_plus,normal_tm,normal_tp)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),spacing(3),dt
    logical, intent(in) :: steepening,flattening
    character(len=*), intent(in) :: solver
    real(dp), intent(out) :: xflux(:,0:,:,:),yflux(:,:,0:,:),zflux(:,:,:,0:)
    logical, intent(out) :: ok
    real(dp), allocatable, optional, intent(out) :: normal_minus(:,:,:,:,:),normal_plus(:,:,:,:,:)
    real(dp), allocatable, optional, intent(out) :: normal_tm(:,:,:,:),normal_tp(:,:,:,:)
    real(dp), allocatable :: nm(:,:,:,:,:),np(:,:,:,:,:),ntm(:,:,:,:),ntp(:,:,:,:)
    logical :: save_states
    real(dp), allocatable :: q(:,:,:,:),sound(:,:,:),minus(:,:,:,:),plus(:,:,:,:),tm(:,:,:),tp(:,:,:)
    real(dp), allocatable :: fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:)
    real(dp) :: stencil(reactive_nprim(size(species)),-3:3)
    real(dp) :: qm(reactive_nprim(size(species))),qp(reactive_nprim(size(species)))
    real(dp) :: recovered,speed,swap,flux(reactive_nvar(size(species)))
    integer :: n(3),nv,axis,i,j,k,offset,index(3),donor(3),left(3),face(3)
    logical :: valid
    xflux=0;yflux=0;zflux=0;ok=.false.;nv=reactive_nvar(size(species))
    save_states=present(normal_minus).or.present(normal_plus).or.present(normal_tm).or.present(normal_tp)
    if (save_states) then
      if (.not.(present(normal_minus).and.present(normal_plus).and.present(normal_tm).and.present(normal_tp))) return
    end if
    n=shape(temperature)
    if (any(n<2).or.any(shape(state)/=[nv,n])) return
    if (any(shape(xflux)/=[nv,n(1)+1,n(2),n(3)])) return
    if (any(shape(yflux)/=[nv,n(1),n(2)+1,n(3)])) return
    if (any(shape(zflux)/=[nv,n(1),n(2),n(3)+1])) return
    if (.not.all(ieee_is_finite(spacing)).or..not.ieee_is_finite(dt)) return
    if (any(spacing<=0).or.dt<0) return
    allocate(fx(nv,0:n(1),n(2),n(3)),fy(nv,n(1),0:n(2),n(3)),fz(nv,n(1),n(2),0:n(3)))
    allocate(q(nv,n(1),n(2),n(3)),sound(n(1),n(2),n(3)))
    allocate(minus(nv,n(1),n(2),n(3)),plus(nv,n(1),n(2),n(3)),tm(n(1),n(2),n(3)),tp(n(1),n(2),n(3)))
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      call reactive_conserved_to_primitive(species,state(:,i,j,k),temperature(i,j,k), &
        q(:,i,j,k),recovered,sound(i,j,k),valid)
      if (.not.valid) return
    end do
    end do
    end do
    if (save_states) allocate(nm(nv,n(1),n(2),n(3),3),np(nv,n(1),n(2),n(3),3), &
      ntm(n(1),n(2),n(3),3),ntp(n(1),n(2),n(3),3))
    do axis=1,3
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        index=[i,j,k]
        do offset=-3,3
          donor=index;donor(axis)=modulo(index(axis)-1+offset,n(axis))+1
          stencil(:,offset)=q(:,donor(1),donor(2),donor(3))
          swap=stencil(2,offset);stencil(2,offset)=stencil(axis+1,offset);stencil(axis+1,offset)=swap
        end do
        if (dt>spacing(axis)/(abs(stencil(2,0))+sound(i,j,k))) return
        call reconstruct_reactive_characteristic_ppm_cell_2d(species,stencil,sound(i,j,k),dt/spacing(axis), &
          steepening,flattening,qm,qp,valid)
        if (.not.valid) return
        swap=qm(2);qm(2)=qm(axis+1);qm(axis+1)=swap
        swap=qp(2);qp(2)=qp(axis+1);qp(axis+1)=swap
        call reactive_primitive_to_conserved(species,qm,minus(:,i,j,k),tm(i,j,k),speed,valid)
        if (.not.valid) return
        call reactive_primitive_to_conserved(species,qp,plus(:,i,j,k),tp(i,j,k),speed,valid)
        if (.not.valid) return
      end do
      end do
      end do
      if (save_states) then
        nm(:,:,:,:,axis)=minus;np(:,:,:,:,axis)=plus;ntm(:,:,:,axis)=tm;ntp(:,:,:,axis)=tp
      end if
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        index=[i,j,k];left=index;left(axis)=modulo(index(axis)-2,n(axis))+1
        select case(axis)
        case(1)
          call reactive_riemann_flux_x(species,plus(:,left(1),left(2),left(3)),minus(:,i,j,k), &
            tp(left(1),left(2),left(3)),tm(i,j,k),solver,flux,valid)
        case(2)
          call reactive_riemann_flux_y(species,plus(:,left(1),left(2),left(3)),minus(:,i,j,k), &
            tp(left(1),left(2),left(3)),tm(i,j,k),solver,flux,valid)
        case(3)
          call reactive_riemann_flux_z(species,plus(:,left(1),left(2),left(3)),minus(:,i,j,k), &
            tp(left(1),left(2),left(3)),tm(i,j,k),solver,flux,valid)
        end select
        if (.not.valid) return
        face=index;face(axis)=face(axis)-1
        select case(axis)
        case(1);fx(:,face(1),j,k)=flux
        case(2);fy(:,i,face(2),k)=flux
        case(3);fz(:,i,j,face(3))=flux
        end select
      end do
      end do
      end do
    end do
    fx(:,n(1),:,:)=fx(:,0,:,:)
    fy(:,:,n(2),:)=fy(:,:,0,:)
    fz(:,:,:,n(3))=fz(:,:,:,0)
    xflux=fx;yflux=fy;zflux=fz
    if (save_states) then
      call move_alloc(nm,normal_minus);call move_alloc(np,normal_plus)
      call move_alloc(ntm,normal_tm);call move_alloc(ntp,normal_tp)
    end if
    ok=.true.
  end subroutine compute_reactive_ppm_normal_fluxes_3d

  ! Twelve-solve CTU ordering: 3 normal, 6 single-transverse (dt/3),
  ! then 3 double-transverse (dt/2) solves.
  ! The pelec solver uses separate primitive pressure work; other solvers use
  ! conservative/EOS correction. An explicit override supports kernel audits.
  subroutine compute_reactive_ppm_ctu_fluxes_3d(species,state,temperature,spacing,dt, &
      steepening,flattening,solver,xflux,yflux,zflux,minimum_theta,ok,pressure_tracing,difmag)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),spacing(3),dt
    logical, intent(in) :: steepening,flattening
    character(len=*), intent(in) :: solver
    real(dp), intent(out) :: xflux(:,0:,:,:),yflux(:,:,0:,:),zflux(:,:,:,0:),minimum_theta
    logical, intent(out) :: ok
    logical, optional, intent(in) :: pressure_tracing
    real(dp), optional, intent(in) :: difmag
    real(dp) :: viscosity
    logical :: trace_pressure
    real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:)
    character(len=32) :: faces(6)
    real(dp), allocatable :: fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:),nm(:,:,:,:,:),np(:,:,:,:,:)
    real(dp), allocatable :: ntm(:,:,:,:),ntp(:,:,:,:),normal(:,:,:,:,:),corner(:,:,:,:,:,:),final(:,:,:,:,:)
    real(dp), allocatable :: minus(:,:,:,:),plus(:,:,:,:),tm(:,:,:),tp(:,:,:),delta(:,:,:,:)
    real(dp) :: theta,theta_m,theta_p,zero(reactive_nvar(size(species)))
    integer :: n(3),nv,a,b,c,i,j,k,index(3),lo(3),hi(3)
    logical :: valid
    xflux=0;yflux=0;zflux=0;minimum_theta=0;ok=.false.;n=shape(temperature);nv=reactive_nvar(size(species))
    if (any(n<2)) return
    if (any(shape(xflux)/=[nv,n(1)+1,n(2),n(3)])) return
    if (any(shape(yflux)/=[nv,n(1),n(2)+1,n(3)])) return
    if (any(shape(zflux)/=[nv,n(1),n(2),n(3)+1])) return
    trace_pressure=solver=="pelec"
    if (present(pressure_tracing)) trace_pressure=pressure_tracing
    viscosity=.1_dp
    if (present(difmag)) viscosity=difmag
    if (.not.ieee_is_finite(viscosity).or.viscosity<0) return
    if (trace_pressure.or.viscosity>0) then
      if (any(shape(state)/=[nv,n])) return
      faces='periodic'
      allocate(ghost(nv,n(1)+8,n(2)+8,n(3)+8),gt(n(1)+8,n(2)+8,n(3)+8))
      call fill_reactive_boundary_ghosts_3d(state,temperature,faces,ghost,gt,valid,ghost_width=4)
      if (.not.valid) return
      call compute_reactive_ppm_halo_ctu_fluxes_3d(species,ghost,gt,n,spacing,dt,steepening,flattening,solver, &
        xflux,yflux,zflux,minimum_theta,ok,pressure_tracing=trace_pressure,difmag=viscosity)
      return
    end if
    allocate(fx(nv,0:n(1),n(2),n(3)),fy(nv,n(1),0:n(2),n(3)),fz(nv,n(1),n(2),0:n(3)))
    call compute_reactive_ppm_normal_fluxes_3d(species,state,temperature,spacing,dt,steepening,flattening, &
      solver,fx,fy,fz,valid,nm,np,ntm,ntp)
    if (.not.valid) return
    allocate(normal(nv,n(1),n(2),n(3),3),corner(nv,n(1),n(2),n(3),3,3),final(nv,n(1),n(2),n(3),3))
    allocate(minus(nv,n(1),n(2),n(3)),plus(nv,n(1),n(2),n(3)),delta(nv,n(1),n(2),n(3)))
    allocate(tm(n(1),n(2),n(3)),tp(n(1),n(2),n(3)))
    normal(:,:,:,:,1)=fx(:,1:,:,:);normal(:,:,:,:,2)=fy(:,:,1:,:);normal(:,:,:,:,3)=fz(:,:,:,1:)
    theta=1;zero=0;corner=0
    do a=1,3
      do b=1,3
        if (a==b) cycle
        delta=0
        call add_difference(normal(:,:,:,:,b),b,(dt/spacing(b))/3)
        call correct_and_solve(a,corner(:,:,:,:,a,b),valid)
        if (.not.valid) return
      end do
    end do
    do a=1,3
      b=modulo(a,3)+1;c=modulo(a+1,3)+1
      delta=0
      call add_difference(corner(:,:,:,:,b,c),b,.5_dp*dt/spacing(b))
      call add_difference(corner(:,:,:,:,c,b),c,.5_dp*dt/spacing(c))
      call correct_and_solve(a,final(:,:,:,:,a),valid)
      if (.not.valid) return
    end do
    xflux(:,1:,:,:)=final(:,:,:,:,1);xflux(:,0,:,:)=final(:,n(1),:,:,1)
    yflux(:,:,1:,:)=final(:,:,:,:,2);yflux(:,:,0,:)=final(:,:,n(2),:,2)
    zflux(:,:,:,1:)=final(:,:,:,:,3);zflux(:,:,:,0)=final(:,:,:,n(3),3)
    minimum_theta=theta;ok=.true.
  contains
    subroutine add_difference(flux,direction,scale)
      real(dp), intent(in) :: flux(:,:,:,:),scale
      integer, intent(in) :: direction
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        index=[i,j,k];lo=index;lo(direction)=modulo(index(direction)-2,n(direction))+1
        delta(:,i,j,k)=delta(:,i,j,k)+scale*(flux(:,i,j,k)-flux(:,lo(1),lo(2),lo(3)))
      end do
      end do
      end do
    end subroutine
    subroutine correct_and_solve(direction,flux,good)
      integer, intent(in) :: direction
      real(dp), intent(out) :: flux(:,:,:,:)
      logical, intent(out) :: good
      good=.false.
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        call apply_reactive_transverse_correction_2d(species,nm(:,i,j,k,direction),ntm(i,j,k,direction), &
          delta(:,i,j,k),zero,1._dp,minus(:,i,j,k),tm(i,j,k),theta_m,valid)
        if (.not.valid) return
        call apply_reactive_transverse_correction_2d(species,np(:,i,j,k,direction),ntp(i,j,k,direction), &
          delta(:,i,j,k),zero,1._dp,plus(:,i,j,k),tp(i,j,k),theta_p,valid)
        if (.not.valid) return
        theta=min(theta,theta_m,theta_p)
      end do
      end do
      end do
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        index=[i,j,k];hi=index;hi(direction)=modulo(index(direction),n(direction))+1
        select case(direction)
        case(1)
          call reactive_riemann_flux_x(species,plus(:,i,j,k),minus(:,hi(1),hi(2),hi(3)), &
            tp(i,j,k),tm(hi(1),hi(2),hi(3)),solver,flux(:,i,j,k),valid)
        case(2)
          call reactive_riemann_flux_y(species,plus(:,i,j,k),minus(:,hi(1),hi(2),hi(3)), &
            tp(i,j,k),tm(hi(1),hi(2),hi(3)),solver,flux(:,i,j,k),valid)
        case(3)
          call reactive_riemann_flux_z(species,plus(:,i,j,k),minus(:,hi(1),hi(2),hi(3)), &
            tp(i,j,k),tm(hi(1),hi(2),hi(3)),solver,flux(:,i,j,k),valid)
        end select
        if (.not.valid) return
      end do
      end do
      end do
      good=.true.
    end subroutine
  end subroutine compute_reactive_ppm_ctu_fluxes_3d

  ! Four supplied ghost layers; physical cells are 1:n. No periodic wrapping
  ! or ghost filling occurs here. Corner solves are bounded to their consumers.
  subroutine compute_reactive_ppm_halo_ctu_fluxes_3d(species,ghost,gt,n,spacing,dt, &
      steepening,flattening,solver,xflux,yflux,zflux,minimum_theta,ok,pressure_tracing,inflow_faces,difmag,no_slip_faces, &
      reflecting_faces,outflow_faces,outlet_pressure)
    type(nasa7_species), intent(in) :: species(:)
    integer, intent(in) :: n(3)
    real(dp), intent(in) :: ghost(:,-3:,-3:,-3:),gt(-3:,-3:,-3:),spacing(3),dt
    logical, intent(in) :: steepening,flattening
    character(len=*), intent(in) :: solver
    real(dp), intent(out) :: xflux(:,0:,:,:),yflux(:,:,0:,:),zflux(:,:,:,0:),minimum_theta
    logical, intent(out) :: ok
    logical, optional, intent(in) :: pressure_tracing,inflow_faces(6),no_slip_faces(6)
    logical, optional, intent(in) :: reflecting_faces(6),outflow_faces(6)
    real(dp), optional, intent(in) :: difmag,outlet_pressure(6)
    real(dp) :: viscosity
    logical :: walls(6),reflecting(6),outflow(6)
    real(dp), allocatable :: cell_flat(:,:,:)
    real(dp) :: pressure_pencil(-3:3),velocity_pencil(-3:3)
    logical :: trace_pressure
    real(dp), allocatable :: pnminus(:,:,:,:),pnplus(:,:,:,:),normal_face(:,:,:,:,:),corner_face(:,:,:,:,:,:)
    real(dp), allocatable :: face_result(:,:,:,:)
    real(dp) :: pl(2),ph(2),vl(2),vh(2),weights(2),gamma,internal_energy
    real(dp), allocatable :: q(:,:,:,:),sound(:,:,:),nm(:,:,:,:,:),np(:,:,:,:,:)
    real(dp), allocatable :: ntm(:,:,:,:),ntp(:,:,:,:),normal(:,:,:,:,:),corner(:,:,:,:,:,:)
    real(dp), allocatable :: minus(:,:,:,:),plus(:,:,:,:),tm(:,:,:),tp(:,:,:),final(:,:,:,:,:)
    real(dp) :: stencil(reactive_nprim(size(species)),-3:3),qm(reactive_nprim(size(species)))
    real(dp) :: qp(reactive_nprim(size(species))),delta(reactive_nvar(size(species)))
    real(dp) :: zero(reactive_nvar(size(species))),theta,theta_m,theta_p,recovered,speed,swap
    integer :: nv,i,j,k,a,b,c,offset,idx(3),donor(3),lo(3),hi(3),lower(3),upper(3)
    logical :: valid
    ok=.false.;xflux=0;yflux=0;zflux=0;minimum_theta=0
    trace_pressure=solver=="pelec"
    if (present(pressure_tracing)) trace_pressure=pressure_tracing
    if (trace_pressure.and.solver/="pelec") return
    viscosity=.1_dp
    if (present(difmag)) viscosity=difmag
    if (.not.ieee_is_finite(viscosity).or.viscosity<0) return
    reflecting=.false.;outflow=.false.
    if(present(reflecting_faces)) reflecting=reflecting_faces
    if(present(outflow_faces)) outflow=outflow_faces
    if(any(reflecting.and.outflow)) return
    if(present(outlet_pressure)) then
      if(.not.all(ieee_is_finite(outlet_pressure))) return
      ! Prescribed-pressure outlets retain their EOS-adjusted exterior donor.
      outflow=outflow.and.outlet_pressure<=0
    end if
    walls=.false.
    if (present(no_slip_faces)) walls=no_slip_faces
    nv=reactive_nvar(size(species))
    ! Supplied halos support a single owned cell on an MPI slab.
    if (any(n<1)) return
    if (any(shape(ghost)/=[nv,n+8]).or.any(shape(gt)/=n+8)) return
    if (any(shape(xflux)/=[nv,n(1)+1,n(2),n(3)])) return
    if (any(shape(yflux)/=[nv,n(1),n(2)+1,n(3)])) return
    if (any(shape(zflux)/=[nv,n(1),n(2),n(3)+1])) return
    if (.not.all(ieee_is_finite(spacing)).or..not.ieee_is_finite(dt)) return
    if (any(spacing<=0).or.dt<0) return
    allocate(q(nv,-3:n(1)+4,-3:n(2)+4,-3:n(3)+4),sound(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    ! Only the union of normal reconstruction pencils needs EOS recovery.
    q=0
    do k=-3,n(3)+4
    do j=-3,n(2)+4
    do i=-3,n(1)+4
      idx=[i,j,k]
      if (count(idx<0.or.idx>n+1)>1) cycle
      call reactive_conserved_to_primitive(species,ghost(:,i,j,k),gt(i,j,k),q(:,i,j,k),recovered,speed,valid)
      if (.not.valid) return
      if (all(idx>=0.and.idx<=n+1)) sound(i,j,k)=speed
    end do
    end do
    end do
    allocate(nm(nv,0:n(1)+1,0:n(2)+1,0:n(3)+1,3),np(nv,0:n(1)+1,0:n(2)+1,0:n(3)+1,3))
    allocate(ntm(0:n(1)+1,0:n(2)+1,0:n(3)+1,3),ntp(0:n(1)+1,0:n(2)+1,0:n(3)+1,3))
    allocate(normal(nv,0:n(1)+1,0:n(2)+1,0:n(3)+1,3), &
      corner(nv,0:n(1)+1,0:n(2)+1,0:n(3)+1,3,3),final(nv,0:n(1)+1,0:n(2)+1,0:n(3)+1,3))
    allocate(minus(nv,0:n(1)+1,0:n(2)+1,0:n(3)+1),plus(nv,0:n(1)+1,0:n(2)+1,0:n(3)+1))
    allocate(tm(0:n(1)+1,0:n(2)+1,0:n(3)+1),tp(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    if (trace_pressure) then
      allocate(pnminus(0:n(1)+1,0:n(2)+1,0:n(3)+1,3),pnplus(0:n(1)+1,0:n(2)+1,0:n(3)+1,3))
      allocate(normal_face(2,0:n(1)+1,0:n(2)+1,0:n(3)+1,3), &
        corner_face(2,0:n(1)+1,0:n(2)+1,0:n(3)+1,3,3),face_result(2,0:n(1)+1,0:n(2)+1,0:n(3)+1))
      normal_face=0;corner_face=0
      allocate(cell_flat(0:n(1)+1,0:n(2)+1,0:n(3)+1))
      cell_flat=1
      ! PeleC trace_ppm uses the minimum flattening over all spatial axes.
      if (flattening) then
        do k=0,n(3)+1
        do j=0,n(2)+1
        do i=0,n(1)+1
          do b=1,3
            do offset=-3,3
              donor=[i,j,k];donor(b)=donor(b)+offset
              pressure_pencil(offset)=q(5,donor(1),donor(2),donor(3))
              velocity_pencil(offset)=q(b+1,donor(1),donor(2),donor(3))
            end do
            cell_flat(i,j,k)=min(cell_flat(i,j,k), &
              reactive_ppm_flattening_coefficient(pressure_pencil,velocity_pencil))
          end do
        end do
        end do
        end do
      end if
    end if
    normal=0;corner=0;final=0;zero=0;theta=1
    do a=1,3
      do k=0,n(3)+1
      do j=0,n(2)+1
      do i=0,n(1)+1
        idx=[i,j,k]
        do offset=-3,3
          donor=idx;donor(a)=donor(a)+offset
          stencil(:,offset)=q(:,donor(1),donor(2),donor(3))
          swap=stencil(2,offset);stencil(2,offset)=stencil(a+1,offset);stencil(a+1,offset)=swap
        end do
        if (dt>spacing(a)/(abs(stencil(2,0))+sound(i,j,k))) return
        if (trace_pressure) then
          call reconstruct_reactive_characteristic_ppm_cell_2d(species,stencil,sound(i,j,k),dt/spacing(a), &
            steepening,flattening,qm,qp,valid,flattening_override=cell_flat(i,j,k),pelec_projection=.true.)
        else
          call reconstruct_reactive_characteristic_ppm_cell_2d(species,stencil,sound(i,j,k),dt/spacing(a), &
            steepening,flattening,qm,qp,valid)
        end if
        if (.not.valid) return
        swap=qm(2);qm(2)=qm(a+1);qm(a+1)=swap
        swap=qp(2);qp(2)=qp(a+1);qp(a+1)=swap
        if (trace_pressure) then
          pnminus(i,j,k,a)=qm(5);pnplus(i,j,k,a)=qp(5)
        end if
        if(trace_pressure) then
          call reconstructed_pelec_state(species,qm,nm(:,i,j,k,a),ntm(i,j,k,a),speed,valid)
          if (.not.valid) return
          call reconstructed_pelec_state(species,qp,np(:,i,j,k,a),ntp(i,j,k,a),speed,valid)
        else
          call reactive_primitive_to_conserved(species,qm,nm(:,i,j,k,a),ntm(i,j,k,a),speed,valid)
          if (.not.valid) return
          call reactive_primitive_to_conserved(species,qp,np(:,i,j,k,a),ntp(i,j,k,a),speed,valid)
        end if
        if (.not.valid) return
      end do
      end do
      end do
      minus=nm(:,:,:,:,a);plus=np(:,:,:,:,a);tm=ntm(:,:,:,a);tp=ntp(:,:,:,a)
      if (trace_pressure) then
        call normal_primitives(a,valid)
        if (.not.valid) return
      end if
      lower=0;upper=n+1;upper(a)=n(a)
      call solve_faces(a,normal(:,:,:,:,a),valid)
      if (.not.valid) return
      if (trace_pressure) normal_face(:,:,:,:,a)=face_result
    end do
    do a=1,3
      do b=1,3
        if (a==b) cycle
        lower=0;upper=n+1;lower(b)=1;upper(b)=n(b)
        call correct_states(a,b,0,valid)
        if (.not.valid) return
        upper(a)=n(a)
        call solve_faces(a,corner(:,:,:,:,a,b),valid)
        if (.not.valid) return
        if (trace_pressure) corner_face(:,:,:,:,a,b)=face_result
      end do
    end do
    do a=1,3
      b=modulo(a,3)+1;c=modulo(a+1,3)+1
      if (trace_pressure) then
        ! pc_transdd subtracts transverse work in ascending spatial order:
        ! X <- Y,Z; Y <- X,Z; Z <- X,Y. The pressure subtraction is sequential.
        b=merge(2,1,a==1);c=merge(2,3,a==3)
      end if
      lower=1;upper=n;lower(a)=0;upper(a)=n(a)+1
      call correct_states(a,b,c,valid)
      if (.not.valid) return
      upper(a)=n(a)
      call solve_faces(a,final(:,:,:,:,a),valid)
      if (.not.valid) return
    end do
    if (trace_pressure.and.present(inflow_faces)) then
      call replace_inflow_fluxes(valid)
      if (.not.valid) return
    end if
    if (viscosity>0) then
      call adjust_reactive_artificial_viscosity_3d(ghost(:,0:n(1)+1,0:n(2)+1,0:n(3)+1), &
        n,spacing,viscosity,walls,final(:,0:n(1),1:n(2),1:n(3),1), &
        final(:,1:n(1),0:n(2),1:n(3),2),final(:,1:n(1),1:n(2),0:n(3),3),valid)
      if (.not.valid) return
    end if
    if(trace_pressure.and.viscosity==0) then
      do a=1,3
        lower=1;upper=n;lower(a)=0
        do k=lower(3),upper(3)
        do j=lower(2),upper(2)
        do i=lower(1),upper(1)
          ! pc_norm_spec_flx runs even when artificial viscosity is disabled.
          speed=sum(final(6:,i,j,k,a))
          if(speed/=0) final(6:,i,j,k,a)=final(6:,i,j,k,a)*(final(1,i,j,k,a)/speed)
        end do
        end do
        end do
      end do
      if(.not.all(ieee_is_finite(final))) return
    end if
    xflux=final(:,0:n(1),1:n(2),1:n(3),1)
    yflux=final(:,1:n(1),0:n(2),1:n(3),2)
    zflux=final(:,1:n(1),1:n(2),0:n(3),3)
    minimum_theta=theta;ok=.true.
  contains
    subroutine correct_boundary_states(direction,point,left,right)
      integer, intent(in) :: direction,point(3)
      real(dp), intent(inout) :: left(:),right(:)
      integer :: face
      face=0
      if(point(direction)==0) face=2*direction-1
      if(point(direction)==n(direction)) face=2*direction
      if(face==0) return
      if(.not.(reflecting(face).or.outflow(face))) return
      ! pc_cmpflx corrects only hydro variables; passive face weights persist.
      if(mod(face,2)==1) then
        left(1:5)=right(1:5)
        if(reflecting(face)) left(direction+1)=-right(direction+1)
      else
        right(1:5)=left(1:5)
        if(reflecting(face)) right(direction+1)=-left(direction+1)
      end if
    end subroutine correct_boundary_states
    subroutine replace_inflow_fluxes(good)
      logical, intent(out) :: good
      integer :: face,side,direction,cell(3),point(3),ii,jj,kk,shift
      real(dp) :: pencil(nv,-2:2),left(nv),right(nv),unused(nv),rotated(nv)
      good=.false.
      do face=1,6
        if (.not.inflow_faces(face)) cycle
        direction=(face+1)/2
        lower=1;upper=n
        lower(direction)=merge(0,n(direction),modulo(face,2)==1)
        upper(direction)=lower(direction)
        do kk=lower(3),upper(3)
        do jj=lower(2),upper(2)
        do ii=lower(1),upper(1)
          do side=0,1
            cell=[ii,jj,kk];cell(direction)=cell(direction)+side
            do shift=-2,2
              point=cell;point(direction)=point(direction)+shift
              pencil(:,shift)=q(:,point(1),point(2),point(3))
              swap=pencil(2,shift);pencil(2,shift)=pencil(direction+1,shift);pencil(direction+1,shift)=swap
            end do
            if (side==0) then
              call reconstruct_pelec_plm_cell(pencil,sound(cell(1),cell(2),cell(3)),dt/spacing(direction), &
                cell_flat(cell(1),cell(2),cell(3)),unused,left,good)
            else
              call reconstruct_pelec_plm_cell(pencil,sound(cell(1),cell(2),cell(3)),dt/spacing(direction), &
                cell_flat(cell(1),cell(2),cell(3)),right,unused,good)
            end if
            if (.not.good) return
          end do
          point=[ii,jj,kk];point(direction)=point(direction)+1
          call reactive_pelec_primitive_flux_x(species,left,right, &
            .5_dp*(sound(ii,jj,kk)+sound(point(1),point(2),point(3))),rotated,good,passive_upwind=.true.)
          if (.not.good) return
          swap=rotated(2);rotated(2)=rotated(direction+1);rotated(direction+1)=swap
          final(:,ii,jj,kk,direction)=rotated
        end do
        end do
        end do
      end do
      good=.true.
    end subroutine replace_inflow_fluxes
    subroutine normal_primitives(direction,good)
      integer, intent(in) :: direction
      logical, intent(out) :: good
      do k=0,n(3)+1
      do j=0,n(2)+1
      do i=0,n(1)+1
        minus(1,i,j,k)=nm(1,i,j,k,direction)
        minus(2:4,i,j,k)=nm(2:4,i,j,k,direction)/nm(1,i,j,k,direction)
        minus(5,i,j,k)=pnminus(i,j,k,direction)
        minus(6:,i,j,k)=nm(6:,i,j,k,direction)/nm(1,i,j,k,direction)
        plus(1,i,j,k)=np(1,i,j,k,direction)
        plus(2:4,i,j,k)=np(2:4,i,j,k,direction)/np(1,i,j,k,direction)
        plus(5,i,j,k)=pnplus(i,j,k,direction)
        plus(6:,i,j,k)=np(6:,i,j,k,direction)/np(1,i,j,k,direction)
      end do
      end do
      end do
      good=all(ieee_is_finite(minus)).and.all(ieee_is_finite(plus))
    end subroutine normal_primitives
    subroutine correct_states(direction,transverse,second,good)
      integer, intent(in) :: direction,transverse,second
      logical, intent(out) :: good
      good=.false.
      do k=lower(3),upper(3)
      do j=lower(2),upper(2)
      do i=lower(1),upper(1)
        idx=[i,j,k];lo=idx;lo(transverse)=lo(transverse)-1
        if (second==0) then
          delta=(dt/spacing(transverse))/3*(normal(:,i,j,k,transverse)- &
            normal(:,lo(1),lo(2),lo(3),transverse))
        else
          delta=(.5_dp*dt/spacing(transverse))*(corner(:,i,j,k,transverse,second)- &
            corner(:,lo(1),lo(2),lo(3),transverse,second))
          lo=idx;lo(second)=lo(second)-1
          delta=delta+(.5_dp*dt/spacing(second))*(corner(:,i,j,k,second,transverse)- &
            corner(:,lo(1),lo(2),lo(3),second,transverse))
        end if
        if (trace_pressure) then
          gamma=sound(i,j,k)**2*q(1,i,j,k)/q(5,i,j,k)
          lo=idx;lo(transverse)=lo(transverse)-1
          if (second==0) then
            pl(1)=normal_face(1,lo(1),lo(2),lo(3),transverse)
            ph(1)=normal_face(1,i,j,k,transverse)
            vl(1)=normal_face(2,lo(1),lo(2),lo(3),transverse)
            vh(1)=normal_face(2,i,j,k,transverse)
            weights(1)=(dt/spacing(transverse))/3
          else
            pl(1)=corner_face(1,lo(1),lo(2),lo(3),transverse,second)
            ph(1)=corner_face(1,i,j,k,transverse,second)
            vl(1)=corner_face(2,lo(1),lo(2),lo(3),transverse,second)
            vh(1)=corner_face(2,i,j,k,transverse,second)
            lo=idx;lo(second)=lo(second)-1
            pl(2)=corner_face(1,lo(1),lo(2),lo(3),second,transverse)
            ph(2)=corner_face(1,i,j,k,second,transverse)
            vl(2)=corner_face(2,lo(1),lo(2),lo(3),second,transverse)
            vh(2)=corner_face(2,i,j,k,second,transverse)
            weights=[.5_dp*dt/spacing(transverse),.5_dp*dt/spacing(second)]
          end if
          offset=merge(1,2,second==0)
          call reactive_ppm_transverse_pressure_3d(nm(:,i,j,k,direction),pnminus(i,j,k,direction),delta,gamma, &
            pl(:offset),ph(:offset),vl(:offset),vh(:offset),weights(:offset),minus(:,i,j,k),internal_energy,good)
          if (.not.good) return
          call reactive_ppm_transverse_pressure_3d(np(:,i,j,k,direction),pnplus(i,j,k,direction),delta,gamma, &
            pl(:offset),ph(:offset),vl(:offset),vh(:offset),weights(:offset),plus(:,i,j,k),internal_energy,good)
          if (.not.good) return
          cycle
        end if
        call apply_reactive_transverse_correction_2d(species,nm(:,i,j,k,direction),ntm(i,j,k,direction), &
          delta,zero,1._dp,minus(:,i,j,k),tm(i,j,k),theta_m,good)
        if (.not.good) return
        call apply_reactive_transverse_correction_2d(species,np(:,i,j,k,direction),ntp(i,j,k,direction), &
          delta,zero,1._dp,plus(:,i,j,k),tp(i,j,k),theta_p,good)
        if (.not.good) return
        theta=min(theta,theta_m,theta_p)
      end do
      end do
      end do
      good=.true.
    end subroutine
    subroutine solve_faces(direction,flux,good)
      integer, intent(in) :: direction
      real(dp), intent(out) :: flux(:,0:,0:,0:)
      logical, intent(out) :: good
      good=.false.;flux=0
      if (trace_pressure) face_result=0
      do k=lower(3),upper(3)
      do j=lower(2),upper(2)
      do i=lower(1),upper(1)
        hi=[i,j,k];hi(direction)=hi(direction)+1
        if (trace_pressure) then
          qm=plus(:,i,j,k);qp=minus(:,hi(1),hi(2),hi(3))
          call correct_boundary_states(direction,[i,j,k],qm,qp)
          swap=qm(2);qm(2)=qm(direction+1);qm(direction+1)=swap
          swap=qp(2);qp(2)=qp(direction+1);qp(direction+1)=swap
          call reactive_pelec_primitive_flux_x(species,qm,qp,.5_dp*(sound(i,j,k)+sound(hi(1),hi(2),hi(3))), &
            flux(:,i,j,k),good,interface_velocity=face_result(2,i,j,k), &
            interface_pressure=face_result(1,i,j,k),passive_upwind=.true.)
          if (.not.good) return
          swap=flux(2,i,j,k);flux(2,i,j,k)=flux(direction+1,i,j,k);flux(direction+1,i,j,k)=swap
          cycle
        end if
        select case(direction)
        case(1)
          call reactive_riemann_flux_x(species,plus(:,i,j,k),minus(:,hi(1),hi(2),hi(3)), &
            tp(i,j,k),tm(hi(1),hi(2),hi(3)),solver,flux(:,i,j,k),good)
        case(2)
          call reactive_riemann_flux_y(species,plus(:,i,j,k),minus(:,hi(1),hi(2),hi(3)), &
            tp(i,j,k),tm(hi(1),hi(2),hi(3)),solver,flux(:,i,j,k),good)
        case(3)
          call reactive_riemann_flux_z(species,plus(:,i,j,k),minus(:,hi(1),hi(2),hi(3)), &
            tp(i,j,k),tm(hi(1),hi(2),hi(3)),solver,flux(:,i,j,k),good)
        end select
        if (.not.good) return
      end do
      end do
      end do
      good=.true.
    end subroutine
  end subroutine compute_reactive_ppm_halo_ctu_fluxes_3d

  ! One conservative CTU step. A failed final EOS recovery rolls back the
  ! complete field and limiter diagnostic; the caller can retry a smaller dt.
  subroutine advance_reactive_euler_ppm_ctu_3d(species,state,temperature,spacing,dt, &
      steepening,flattening,solver,minimum_theta,ok,difmag)
    real(dp), optional, intent(in) :: difmag
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:),minimum_theta
    real(dp), intent(in) :: spacing(3),dt
    logical, intent(in) :: steepening,flattening
    character(len=*), intent(in) :: solver
    logical, intent(out) :: ok
    real(dp), allocatable :: fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:),candidate(:,:,:,:),ct(:,:,:)
    real(dp) :: theta
    integer :: n(3),nv,i,j,k
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species));n=shape(temperature)
    if (any(n<2).or.any(shape(state)/=[nv,n])) return
    if (.not.all(ieee_is_finite(spacing)).or..not.ieee_is_finite(dt)) return
    if (any(spacing<=0).or.dt<0) return
    allocate(ct(n(1),n(2),n(3)))
    if (dt==0) then
      call recover_reactive_temperatures_3d(species,state,temperature,n(1),n(2),n(3),ct,valid)
      if (.not.valid) return
      minimum_theta=1;ok=.true.;return
    end if
    allocate(fx(nv,0:n(1),n(2),n(3)),fy(nv,n(1),0:n(2),n(3)),fz(nv,n(1),n(2),0:n(3)))
    call compute_reactive_ppm_ctu_fluxes_3d(species,state,temperature,spacing,dt,steepening,flattening, &
      solver,fx,fy,fz,theta,valid,difmag=difmag)
    if (.not.valid) return
    candidate=state
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      candidate(:,i,j,k)=state(:,i,j,k)-(dt/spacing(1))*(fx(:,i,j,k)-fx(:,i-1,j,k))- &
        (dt/spacing(2))*(fy(:,i,j,k)-fy(:,i,j-1,k))-(dt/spacing(3))*(fz(:,i,j,k)-fz(:,i,j,k-1))
    end do
    end do
    end do
    call recover_reactive_temperatures_3d(species,candidate,temperature,n(1),n(2),n(3),ct,valid)
    if (.not.valid) return
    state=candidate;temperature=ct;minimum_theta=theta;ok=.true.
  end subroutine advance_reactive_euler_ppm_ctu_3d

  subroutine compute_reactive_cfl_timestep_3d( &
      species, state, temperature, nx, ny, nz, dx, dy, dz, cfl, dt, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz, cfl
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: local_temperature, sound_speed, rate, maximum_rate
    logical :: cell_ok
    integer :: i, j, k, nvar

    dt = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return
    if (.not. all(ieee_is_finite([dx, dy, dz, cfl]))) return
    if (nx < 2 .or. ny < 2 .or. nz < 2 .or. &
        dx <= 0.0_dp .or. dy <= 0.0_dp .or. dz <= 0.0_dp .or. &
        cfl <= 0.0_dp .or. cfl > 1.0_dp) return

    allocate(primitive(reactive_nprim(size(species))))
    maximum_rate = 0.0_dp
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), primitive, &
            local_temperature, sound_speed, cell_ok)
          if (.not. cell_ok) return
          rate = (abs(primitive(2)) + sound_speed) / dx + &
            (abs(primitive(3)) + sound_speed) / dy + &
            (abs(primitive(4)) + sound_speed) / dz
          maximum_rate = max(maximum_rate, rate)
        end do
      end do
    end do
    if (.not. ieee_is_finite(maximum_rate)) return
    if (maximum_rate <= 0.0_dp) return
    dt = cfl / maximum_rate
    ok = ieee_is_finite(dt)
    if (ok) ok = dt > 0.0_dp
  end subroutine compute_reactive_cfl_timestep_3d

  subroutine compute_reactive_euler_rhs_3d( &
      species, state, temperature, nx, ny, nz, dx, dy, dz, &
      riemann_solver, rhs, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz
    character(len=*), intent(in) :: riemann_solver
    real(dp), intent(out) :: rhs(:, :, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: flux_x(:, :, :, :)
    real(dp), allocatable :: flux_y(:, :, :, :)
    real(dp), allocatable :: flux_z(:, :, :, :)
    logical :: face_ok
    integer :: i, j, k, next_i, next_j, next_k
    integer :: previous_i, previous_j, previous_k, nvar

    rhs = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return
    if (size(rhs, 1) /= nvar .or. size(rhs, 2) /= nx .or. &
        size(rhs, 3) /= ny .or. size(rhs, 4) /= nz) return
    if (.not. all(ieee_is_finite([dx, dy, dz]))) return
    if (nx < 2 .or. ny < 2 .or. nz < 2 .or. &
        dx <= 0.0_dp .or. dy <= 0.0_dp .or. dz <= 0.0_dp) return

    allocate(flux_x(nvar, nx, ny, nz))
    allocate(flux_y(nvar, nx, ny, nz))
    allocate(flux_z(nvar, nx, ny, nz))

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          next_i = modulo(i, nx) + 1
          call reactive_riemann_flux_x( &
            species, state(:, i, j, k), state(:, next_i, j, k), &
            temperature(i, j, k), temperature(next_i, j, k), &
            riemann_solver, flux_x(:, i, j, k), face_ok)
          if (.not. face_ok) return
        end do
      end do
    end do

    do k = 1, nz
      do j = 1, ny
        next_j = modulo(j, ny) + 1
        do i = 1, nx
          call reactive_riemann_flux_y( &
            species, state(:, i, j, k), state(:, i, next_j, k), &
            temperature(i, j, k), temperature(i, next_j, k), &
            riemann_solver, flux_y(:, i, j, k), face_ok)
          if (.not. face_ok) return
        end do
      end do
    end do

    do k = 1, nz
      next_k = modulo(k, nz) + 1
      do j = 1, ny
        do i = 1, nx
          call reactive_riemann_flux_z( &
            species, state(:, i, j, k), state(:, i, j, next_k), &
            temperature(i, j, k), temperature(i, j, next_k), &
            riemann_solver, flux_z(:, i, j, k), face_ok)
          if (.not. face_ok) return
        end do
      end do
    end do

    do k = 1, nz
      previous_k = modulo(k - 2, nz) + 1
      do j = 1, ny
        previous_j = modulo(j - 2, ny) + 1
        do i = 1, nx
          previous_i = modulo(i - 2, nx) + 1
          rhs(:, i, j, k) = &
            -(flux_x(:, i, j, k) - &
              flux_x(:, previous_i, j, k)) / dx &
            -(flux_y(:, i, j, k) - &
              flux_y(:, i, previous_j, k)) / dy &
            -(flux_z(:, i, j, k) - &
              flux_z(:, i, j, previous_k)) / dz
        end do
      end do
    end do
    ok = all(ieee_is_finite(rhs))
  end subroutine compute_reactive_euler_rhs_3d

  subroutine advance_reactive_euler_ssprk2_3d( &
      species, state, temperature, nx, ny, nz, dx, dy, dz, dt, &
      riemann_solver, ok, difmag)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz, dt
    character(len=*), intent(in) :: riemann_solver
    real(dp), optional, intent(in) :: difmag
    logical, intent(out) :: ok

    real(dp), allocatable :: old_state(:, :, :, :)
    real(dp), allocatable :: stage_state(:, :, :, :)
    real(dp), allocatable :: updated_state(:, :, :, :)
    real(dp), allocatable :: rhs(:, :, :, :)
    real(dp), allocatable :: old_temperature(:, :, :)
    real(dp), allocatable :: stage_temperature(:, :, :)
    real(dp), allocatable :: updated_temperature(:, :, :)
    logical :: local_ok
    real(dp) :: coefficient, boundary_change(size(state,1))
    character(len=8) :: faces(6)
    integer :: nvar

    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return
    if (.not. ieee_is_finite(dt)) return
    if (nx < 2 .or. ny < 2 .or. nz < 2 .or. dt <= 0.0_dp) return

    coefficient=.1_dp
    if (present(difmag)) coefficient=difmag
    if (.not.ieee_is_finite(coefficient).or.coefficient<0) return
    if (coefficient>0) then
      faces='periodic'
      call advance_reactive_boundary_hydro_3d(species,state,temperature,[dx,dy,dz],dt, &
        faces,'pcm','mc',riemann_solver,boundary_change,ok, &
        difmag=coefficient)
      return
    end if

    allocate(old_state(nvar, nx, ny, nz))
    allocate(stage_state(nvar, nx, ny, nz))
    allocate(updated_state(nvar, nx, ny, nz))
    allocate(rhs(nvar, nx, ny, nz))
    allocate(old_temperature(nx, ny, nz))
    allocate(stage_temperature(nx, ny, nz))
    allocate(updated_temperature(nx, ny, nz))
    old_state = state
    call recover_reactive_temperatures_3d( &
      species, old_state, temperature, nx, ny, nz, old_temperature, local_ok)
    if (.not. local_ok) return

    call compute_reactive_euler_rhs_3d( &
      species, old_state, old_temperature, nx, ny, nz, dx, dy, dz, &
      riemann_solver, rhs, local_ok)
    if (.not. local_ok) return
    stage_state = old_state + dt * rhs
    call recover_reactive_temperatures_3d( &
      species, stage_state, old_temperature, nx, ny, nz, &
      stage_temperature, local_ok)
    if (.not. local_ok) return

    call compute_reactive_euler_rhs_3d( &
      species, stage_state, stage_temperature, nx, ny, nz, dx, dy, dz, &
      riemann_solver, rhs, local_ok)
    if (.not. local_ok) return
    updated_state = 0.5_dp * old_state + &
      0.5_dp * (stage_state + dt * rhs)
    call recover_reactive_temperatures_3d( &
      species, updated_state, stage_temperature, nx, ny, nz, &
      updated_temperature, local_ok)
    if (.not. local_ok) return

    state = updated_state
    temperature = updated_temperature
    ok = .true.
  end subroutine advance_reactive_euler_ssprk2_3d

  subroutine compute_reactive_face_fluxes_3d( &
      species, state, temperature, nx, ny, nz, riemann_solver, &
      flux_x, flux_y, flux_z, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    character(len=*), intent(in) :: riemann_solver
    real(dp), intent(out) :: flux_x(:, :, :, :)
    real(dp), intent(out) :: flux_y(:, :, :, :)
    real(dp), intent(out) :: flux_z(:, :, :, :)
    logical, intent(out) :: ok

    logical :: face_ok
    integer :: i, j, k, next_i, next_j, next_k, nvar

    flux_x = 0.0_dp
    flux_y = 0.0_dp
    flux_z = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return
    if (nx < 2 .or. ny < 2 .or. nz < 2 .or. &
        any(shape(flux_x) /= shape(state)) .or. &
        any(shape(flux_y) /= shape(state)) .or. &
        any(shape(flux_z) /= shape(state))) return

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          next_i = modulo(i, nx) + 1
          call reactive_riemann_flux_x( &
            species, state(:, i, j, k), state(:, next_i, j, k), &
            temperature(i, j, k), temperature(next_i, j, k), &
            riemann_solver, flux_x(:, i, j, k), face_ok)
          if (.not. face_ok) return
          next_j = modulo(j, ny) + 1
          call reactive_riemann_flux_y( &
            species, state(:, i, j, k), state(:, i, next_j, k), &
            temperature(i, j, k), temperature(i, next_j, k), &
            riemann_solver, flux_y(:, i, j, k), face_ok)
          if (.not. face_ok) return
          next_k = modulo(k, nz) + 1
          call reactive_riemann_flux_z( &
            species, state(:, i, j, k), state(:, i, j, next_k), &
            temperature(i, j, k), temperature(i, j, next_k), &
            riemann_solver, flux_z(:, i, j, k), face_ok)
          if (.not. face_ok) return
        end do
      end do
    end do
    ok = all(ieee_is_finite(flux_x)) .and. &
      all(ieee_is_finite(flux_y)) .and. all(ieee_is_finite(flux_z))
  end subroutine compute_reactive_face_fluxes_3d

  subroutine compute_reactive_plm_face_fluxes_3d( &
      species, state, temperature, nx, ny, nz, limiter, riemann_solver, &
      flux_x, flux_y, flux_z, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    character(len=*), intent(in) :: limiter, riemann_solver
    real(dp), intent(out) :: flux_x(:, :, :, :)
    real(dp), intent(out) :: flux_y(:, :, :, :)
    real(dp), intent(out) :: flux_z(:, :, :, :)
    logical, intent(out) :: ok

    call compute_reactive_plm_slab_face_fluxes_3d( &
      species, state, temperature, nx, ny, nz, 1, nx, limiter, &
      riemann_solver, flux_x, flux_y, flux_z, ok)
  end subroutine compute_reactive_plm_face_fluxes_3d

  subroutine compute_reactive_plm_slab_face_fluxes_3d( &
      species, state, temperature, nx, ny, nz, first_i, last_i, &
      limiter, riemann_solver, flux_x, flux_y, flux_z, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz, first_i, last_i
    character(len=*), intent(in) :: limiter, riemann_solver
    real(dp), intent(out) :: flux_x(:, :, :, :)
    real(dp), intent(out) :: flux_y(:, :, :, :)
    real(dp), intent(out) :: flux_z(:, :, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:, :, :, :), sound_speed(:, :, :)
    real(dp), allocatable :: x_minus(:, :, :, :), x_plus(:, :, :, :)
    real(dp), allocatable :: y_minus(:, :, :, :), y_plus(:, :, :, :)
    real(dp), allocatable :: z_minus(:, :, :, :), z_plus(:, :, :, :)
    real(dp), allocatable :: candidate_x(:, :, :, :)
    real(dp), allocatable :: candidate_y(:, :, :, :)
    real(dp), allocatable :: candidate_z(:, :, :, :)
    real(dp), allocatable :: center(:), dl(:), dr(:), slope(:)
    real(dp), allocatable :: rotated_center(:), rotated_dl(:), rotated_dr(:)
    real(dp), allocatable :: rotated_slope(:), rotated_minus(:)
    real(dp), allocatable :: rotated_plus(:)
    real(dp), allocatable :: left_state(:), right_state(:)
    logical, allocatable :: primitive_needed(:), reconstruction_needed(:)
    real(dp) :: recovered_temperature, theta
    real(dp) :: left_temperature, right_temperature
    logical :: local_ok
    integer :: i, j, k, im, ip, jm, jp, km, kp
    integer :: nvar, nprimitive

    flux_x = 0.0_dp
    flux_y = 0.0_dp
    flux_z = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    nprimitive = reactive_nprim(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return
    if (nx < 3 .or. ny < 3 .or. nz < 3 .or. &
        any(shape(flux_x) /= shape(state)) .or. &
        any(shape(flux_y) /= shape(state)) .or. &
        any(shape(flux_z) /= shape(state))) return
    if (first_i < 1 .or. last_i > nx .or. last_i < first_i) return
    if (trim(limiter) /= "minmod" .and. trim(limiter) /= "mc") return

    allocate(primitive(nprimitive, nx, ny, nz), sound_speed(nx, ny, nz))
    allocate(x_minus(nprimitive, nx, ny, nz))
    allocate(x_plus(nprimitive, nx, ny, nz))
    allocate(y_minus(nprimitive, nx, ny, nz))
    allocate(y_plus(nprimitive, nx, ny, nz))
    allocate(z_minus(nprimitive, nx, ny, nz))
    allocate(z_plus(nprimitive, nx, ny, nz))
    allocate(candidate_x(nvar, nx, ny, nz))
    allocate(candidate_y(nvar, nx, ny, nz))
    allocate(candidate_z(nvar, nx, ny, nz))
    allocate(center(nprimitive), dl(nprimitive), dr(nprimitive))
    allocate(slope(nprimitive), rotated_center(nprimitive))
    allocate(rotated_dl(nprimitive), rotated_dr(nprimitive))
    allocate(rotated_slope(nprimitive), rotated_minus(nprimitive))
    allocate(rotated_plus(nprimitive))
    allocate(left_state(nvar), right_state(nvar))
    allocate(primitive_needed(nx), reconstruction_needed(nx))

    primitive_needed = .false.
    reconstruction_needed = .false.
    do i = first_i, last_i
      reconstruction_needed(i) = .true.
      reconstruction_needed(periodic_index_3d(i + 1, nx)) = .true.
    end do
    do i = 1, nx
      if (.not. reconstruction_needed(i)) cycle
      primitive_needed(i) = .true.
      primitive_needed(periodic_index_3d(i - 1, nx)) = .true.
      primitive_needed(periodic_index_3d(i + 1, nx)) = .true.
    end do

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (.not. primitive_needed(i)) cycle
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), &
            primitive(:, i, j, k), recovered_temperature, &
            sound_speed(i, j, k), local_ok)
          if (.not. local_ok) return
        end do
      end do
    end do

    do k = 1, nz
      km = periodic_index_3d(k - 1, nz)
      kp = periodic_index_3d(k + 1, nz)
      do j = 1, ny
        jm = periodic_index_3d(j - 1, ny)
        jp = periodic_index_3d(j + 1, ny)
        do i = 1, nx
          if (.not. reconstruction_needed(i)) cycle
          im = periodic_index_3d(i - 1, nx)
          ip = periodic_index_3d(i + 1, nx)
          center = primitive(:, i, j, k)
          dl = center - primitive(:, im, j, k)
          dr = primitive(:, ip, j, k) - center
          call characteristic_limited_slope( &
            center, dl, dr, sound_speed(i, j, k), limiter, slope, local_ok)
          if (.not. local_ok) return
          theta = primitive_slope_scale_3d(center, slope, size(species))
          slope = theta * slope
          x_minus(:, i, j, k) = center - 0.5_dp * slope
          x_plus(:, i, j, k) = center + 0.5_dp * slope
          call sanitize_primitive_3d( &
            x_minus(:, i, j, k), center, size(species))
          call sanitize_primitive_3d( &
            x_plus(:, i, j, k), center, size(species))

          call rotate_primitive_y_to_x_3d(center, rotated_center)
          call rotate_primitive_y_to_x_3d( &
            center - primitive(:, i, jm, k), rotated_dl)
          call rotate_primitive_y_to_x_3d( &
            primitive(:, i, jp, k) - center, rotated_dr)
          call characteristic_limited_slope( &
            rotated_center, rotated_dl, rotated_dr, sound_speed(i, j, k), &
            limiter, rotated_slope, local_ok)
          if (.not. local_ok) return
          theta = primitive_slope_scale_3d( &
            rotated_center, rotated_slope, size(species))
          rotated_slope = theta * rotated_slope
          rotated_minus = rotated_center - 0.5_dp * rotated_slope
          rotated_plus = rotated_center + 0.5_dp * rotated_slope
          call rotate_primitive_x_to_y_3d( &
            rotated_minus, y_minus(:, i, j, k))
          call rotate_primitive_x_to_y_3d( &
            rotated_plus, y_plus(:, i, j, k))
          call sanitize_primitive_3d( &
            y_minus(:, i, j, k), center, size(species))
          call sanitize_primitive_3d( &
            y_plus(:, i, j, k), center, size(species))

          call rotate_primitive_z_to_x_3d(center, rotated_center)
          call rotate_primitive_z_to_x_3d( &
            center - primitive(:, i, j, km), rotated_dl)
          call rotate_primitive_z_to_x_3d( &
            primitive(:, i, j, kp) - center, rotated_dr)
          call characteristic_limited_slope( &
            rotated_center, rotated_dl, rotated_dr, sound_speed(i, j, k), &
            limiter, rotated_slope, local_ok)
          if (.not. local_ok) return
          theta = primitive_slope_scale_3d( &
            rotated_center, rotated_slope, size(species))
          rotated_slope = theta * rotated_slope
          rotated_minus = rotated_center - 0.5_dp * rotated_slope
          rotated_plus = rotated_center + 0.5_dp * rotated_slope
          call rotate_primitive_x_to_z_3d( &
            rotated_minus, z_minus(:, i, j, k))
          call rotate_primitive_x_to_z_3d( &
            rotated_plus, z_plus(:, i, j, k))
          call sanitize_primitive_3d( &
            z_minus(:, i, j, k), center, size(species))
          call sanitize_primitive_3d( &
            z_plus(:, i, j, k), center, size(species))
        end do
      end do
    end do

    candidate_x = 0.0_dp
    candidate_y = 0.0_dp
    candidate_z = 0.0_dp
    do k = 1, nz
      kp = periodic_index_3d(k + 1, nz)
      do j = 1, ny
        jp = periodic_index_3d(j + 1, ny)
        do i = first_i, last_i
          ip = periodic_index_3d(i + 1, nx)
          call primitive_face_to_state_3d( &
            species, x_plus(:, i, j, k), primitive(:, i, j, k), &
            left_state, left_temperature, local_ok)
          if (.not. local_ok) return
          call primitive_face_to_state_3d( &
            species, x_minus(:, ip, j, k), primitive(:, ip, j, k), &
            right_state, right_temperature, local_ok)
          if (.not. local_ok) return
          call reactive_riemann_flux_x( &
            species, left_state, right_state, left_temperature, &
            right_temperature, riemann_solver, candidate_x(:, i, j, k), &
            local_ok)
          if (.not. local_ok) return

          call primitive_face_to_state_3d( &
            species, y_plus(:, i, j, k), primitive(:, i, j, k), &
            left_state, left_temperature, local_ok)
          if (.not. local_ok) return
          call primitive_face_to_state_3d( &
            species, y_minus(:, i, jp, k), primitive(:, i, jp, k), &
            right_state, right_temperature, local_ok)
          if (.not. local_ok) return
          call reactive_riemann_flux_y( &
            species, left_state, right_state, left_temperature, &
            right_temperature, riemann_solver, candidate_y(:, i, j, k), &
            local_ok)
          if (.not. local_ok) return

          call primitive_face_to_state_3d( &
            species, z_plus(:, i, j, k), primitive(:, i, j, k), &
            left_state, left_temperature, local_ok)
          if (.not. local_ok) return
          call primitive_face_to_state_3d( &
            species, z_minus(:, i, j, kp), primitive(:, i, j, kp), &
            right_state, right_temperature, local_ok)
          if (.not. local_ok) return
          call reactive_riemann_flux_z( &
            species, left_state, right_state, left_temperature, &
            right_temperature, riemann_solver, candidate_z(:, i, j, k), &
            local_ok)
          if (.not. local_ok) return
        end do
      end do
    end do
    if (.not. all(ieee_is_finite(candidate_x)) .or. &
        .not. all(ieee_is_finite(candidate_y)) .or. &
        .not. all(ieee_is_finite(candidate_z))) return
    flux_x = candidate_x
    flux_y = candidate_y
    flux_z = candidate_z
    ok = .true.
  end subroutine compute_reactive_plm_slab_face_fluxes_3d

  pure integer function periodic_index_3d(index, extent) result(wrapped)
    integer, intent(in) :: index, extent

    wrapped = 1 + modulo(index - 1, extent)
  end function periodic_index_3d

  pure subroutine rotate_primitive_y_to_x_3d(input, output)
    real(dp), intent(in) :: input(:)
    real(dp), intent(out) :: output(:)

    output = input
    if (size(input) /= size(output) .or. size(input) < 5) return
    output(2) = input(3)
    output(3) = input(2)
  end subroutine rotate_primitive_y_to_x_3d

  pure subroutine rotate_primitive_x_to_y_3d(input, output)
    real(dp), intent(in) :: input(:)
    real(dp), intent(out) :: output(:)

    call rotate_primitive_y_to_x_3d(input, output)
  end subroutine rotate_primitive_x_to_y_3d

  pure subroutine rotate_primitive_z_to_x_3d(input, output)
    real(dp), intent(in) :: input(:)
    real(dp), intent(out) :: output(:)

    output = input
    if (size(input) /= size(output) .or. size(input) < 5) return
    output(2) = input(4)
    output(4) = input(2)
  end subroutine rotate_primitive_z_to_x_3d

  pure subroutine rotate_primitive_x_to_z_3d(input, output)
    real(dp), intent(in) :: input(:)
    real(dp), intent(out) :: output(:)

    call rotate_primitive_z_to_x_3d(input, output)
  end subroutine rotate_primitive_x_to_z_3d

  pure real(dp) function primitive_lower_scale_3d( &
      center, slope, lower) result(theta)
    real(dp), intent(in) :: center, slope, lower
    real(dp) :: magnitude

    magnitude = abs(slope)
    if (magnitude <= tiny(1.0_dp) .or. center - magnitude > lower) then
      theta = 1.0_dp
    else
      theta = max(0.0_dp, min(1.0_dp, (center - lower) / magnitude))
    end if
  end function primitive_lower_scale_3d

  pure real(dp) function primitive_upper_scale_3d( &
      center, slope, upper) result(theta)
    real(dp), intent(in) :: center, slope, upper
    real(dp) :: magnitude

    magnitude = abs(slope)
    if (magnitude <= tiny(1.0_dp) .or. center + magnitude < upper) then
      theta = 1.0_dp
    else
      theta = max(0.0_dp, min(1.0_dp, (upper - center) / magnitude))
    end if
  end function primitive_upper_scale_3d

  pure real(dp) function primitive_slope_scale_3d( &
      center, slope, nspecies) result(theta)
    real(dp), intent(in) :: center(:), slope(:)
    integer, intent(in) :: nspecies
    integer :: species_index, component

    theta = 1.0_dp
    if (size(center) /= size(slope)) then
      theta = 0.0_dp
      return
    end if
    theta = min(theta, &
      primitive_lower_scale_3d(center(1), slope(1), density_floor))
    theta = min(theta, &
      primitive_lower_scale_3d(center(5), slope(5), pressure_floor))
    do species_index = 1, nspecies
      component = reactive_mass_fraction_component(species_index)
      theta = min(theta, &
        primitive_lower_scale_3d(center(component), slope(component), 0.0_dp))
      theta = min(theta, &
        primitive_upper_scale_3d(center(component), slope(component), 1.0_dp))
    end do
  end function primitive_slope_scale_3d

  pure subroutine sanitize_primitive_3d(q, fallback, nspecies)
    real(dp), intent(inout) :: q(:)
    real(dp), intent(in) :: fallback(:)
    integer, intent(in) :: nspecies

    real(dp) :: total
    integer :: species_index, component

    if (size(q) /= size(fallback) .or. size(q) < 5) return
    if (.not. all(ieee_is_finite(q))) then
      q = fallback
      return
    end if
    if (q(1) <= density_floor .or. q(5) <= pressure_floor) then
      q = fallback
      return
    end if
    total = 0.0_dp
    do species_index = 1, nspecies
      component = reactive_mass_fraction_component(species_index)
      if (q(component) < -1.0e-12_dp) then
        q = fallback
        return
      end if
      q(component) = max(0.0_dp, q(component))
      total = total + q(component)
    end do
    if (total <= tiny(1.0_dp)) then
      q = fallback
      return
    end if
    do species_index = 1, nspecies
      component = reactive_mass_fraction_component(species_index)
      q(component) = q(component) / total
    end do
  end subroutine sanitize_primitive_3d

  subroutine primitive_face_to_state_3d( &
      species, face_primitive, fallback_primitive, state, temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: face_primitive(:), fallback_primitive(:)
    real(dp), intent(out) :: state(:), temperature
    logical, intent(out) :: ok

    real(dp) :: work(size(face_primitive))
    real(dp) :: sound_speed

    work = face_primitive
    call sanitize_primitive_3d(work, fallback_primitive, size(species))
    call reactive_primitive_to_conserved( &
      species, work, state, temperature, sound_speed, ok)
    if (ok) return
    call reactive_primitive_to_conserved( &
      species, fallback_primitive, state, temperature, sound_speed, ok)
  end subroutine primitive_face_to_state_3d

  subroutine reactive_flux_divergence_3d( &
      flux_x, flux_y, flux_z, nx, ny, nz, dx, dy, dz, rhs, ok)
    real(dp), intent(in) :: flux_x(:, :, :, :)
    real(dp), intent(in) :: flux_y(:, :, :, :)
    real(dp), intent(in) :: flux_z(:, :, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz
    real(dp), intent(out) :: rhs(:, :, :, :)
    logical, intent(out) :: ok

    integer :: i, j, k, previous_i, previous_j, previous_k

    rhs = 0.0_dp
    ok = .false.
    if (.not. all(ieee_is_finite([dx, dy, dz]))) return
    ok = nx >= 2 .and. ny >= 2 .and. nz >= 2 .and. &
      dx > 0.0_dp .and. dy > 0.0_dp .and. dz > 0.0_dp
    if (.not. ok) return
    ok = all(shape(flux_x) == shape(rhs)) .and. &
      all(shape(flux_y) == shape(rhs)) .and. &
      all(shape(flux_z) == shape(rhs)) .and. &
      size(rhs, 2) == nx .and. size(rhs, 3) == ny .and. &
      size(rhs, 4) == nz
    if (.not. ok) return
    do k = 1, nz
      previous_k = modulo(k - 2, nz) + 1
      do j = 1, ny
        previous_j = modulo(j - 2, ny) + 1
        do i = 1, nx
          previous_i = modulo(i - 2, nx) + 1
          rhs(:, i, j, k) = &
            -(flux_x(:, i, j, k) - &
              flux_x(:, previous_i, j, k)) / dx &
            -(flux_y(:, i, j, k) - &
              flux_y(:, i, previous_j, k)) / dy &
            -(flux_z(:, i, j, k) - &
              flux_z(:, i, j, previous_k)) / dz
        end do
      end do
    end do
    ok = all(ieee_is_finite(rhs))
  end subroutine reactive_flux_divergence_3d

  subroutine advance_reactive_euler_ssprk2_with_fluxes_3d( &
      species, state, temperature, nx, ny, nz, dx, dy, dz, dt, &
      riemann_solver, flux_x, flux_y, flux_z, ok, difmag)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz, dt
    character(len=*), intent(in) :: riemann_solver
    real(dp), intent(out) :: flux_x(:, :, :, :)
    real(dp), intent(out) :: flux_y(:, :, :, :)
    real(dp), intent(out) :: flux_z(:, :, :, :)
    real(dp), optional, intent(in) :: difmag
    logical, intent(out) :: ok

    real(dp), allocatable :: old_state(:, :, :, :)
    real(dp), allocatable :: stage_state(:, :, :, :)
    real(dp), allocatable :: updated_state(:, :, :, :)
    real(dp), allocatable :: rhs(:, :, :, :)
    real(dp), allocatable :: first_flux_x(:, :, :, :)
    real(dp), allocatable :: first_flux_y(:, :, :, :)
    real(dp), allocatable :: first_flux_z(:, :, :, :)
    real(dp), allocatable :: second_flux_x(:, :, :, :)
    real(dp), allocatable :: second_flux_y(:, :, :, :)
    real(dp), allocatable :: second_flux_z(:, :, :, :)
    real(dp), allocatable :: old_temperature(:, :, :)
    real(dp), allocatable :: stage_temperature(:, :, :)
    real(dp), allocatable :: updated_temperature(:, :, :)
    logical :: local_ok
    real(dp) :: coefficient, boundary_change(size(state,1))
    character(len=8) :: faces(6)
    integer :: nvar

    flux_x = 0.0_dp
    flux_y = 0.0_dp
    flux_z = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return
    if (.not. ieee_is_finite(dt)) return
    if (nx < 2 .or. ny < 2 .or. nz < 2 .or. dt <= 0.0_dp) return
    if (any(shape(flux_x) /= shape(state)) .or. &
        any(shape(flux_y) /= shape(state)) .or. &
        any(shape(flux_z) /= shape(state))) return

    coefficient=.1_dp
    if (present(difmag)) coefficient=difmag
    if (.not.ieee_is_finite(coefficient).or.coefficient<0) return
    if (coefficient>0) then
      faces='periodic'
      call advance_reactive_boundary_hydro_3d(species,state,temperature,[dx,dy,dz],dt, &
        faces,'pcm','mc',riemann_solver,boundary_change,ok, flux_x=flux_x,flux_y=flux_y,flux_z=flux_z, &
        difmag=coefficient)
      return
    end if

    allocate(old_state(nvar, nx, ny, nz))
    allocate(stage_state(nvar, nx, ny, nz))
    allocate(updated_state(nvar, nx, ny, nz))
    allocate(rhs(nvar, nx, ny, nz))
    allocate(first_flux_x(nvar, nx, ny, nz))
    allocate(first_flux_y(nvar, nx, ny, nz))
    allocate(first_flux_z(nvar, nx, ny, nz))
    allocate(second_flux_x(nvar, nx, ny, nz))
    allocate(second_flux_y(nvar, nx, ny, nz))
    allocate(second_flux_z(nvar, nx, ny, nz))
    allocate(old_temperature(nx, ny, nz))
    allocate(stage_temperature(nx, ny, nz))
    allocate(updated_temperature(nx, ny, nz))
    old_state = state
    call recover_reactive_temperatures_3d( &
      species, old_state, temperature, nx, ny, nz, old_temperature, local_ok)
    if (.not. local_ok) return

    call compute_reactive_face_fluxes_3d( &
      species, old_state, old_temperature, nx, ny, nz, riemann_solver, &
      first_flux_x, first_flux_y, first_flux_z, local_ok)
    if (.not. local_ok) return
    call reactive_flux_divergence_3d( &
      first_flux_x, first_flux_y, first_flux_z, nx, ny, nz, &
      dx, dy, dz, rhs, local_ok)
    if (.not. local_ok) return
    stage_state = old_state + dt * rhs
    call recover_reactive_temperatures_3d( &
      species, stage_state, old_temperature, nx, ny, nz, &
      stage_temperature, local_ok)
    if (.not. local_ok) return

    call compute_reactive_face_fluxes_3d( &
      species, stage_state, stage_temperature, nx, ny, nz, riemann_solver, &
      second_flux_x, second_flux_y, second_flux_z, local_ok)
    if (.not. local_ok) return
    call reactive_flux_divergence_3d( &
      second_flux_x, second_flux_y, second_flux_z, nx, ny, nz, &
      dx, dy, dz, rhs, local_ok)
    if (.not. local_ok) return
    updated_state = 0.5_dp * old_state + &
      0.5_dp * (stage_state + dt * rhs)
    call recover_reactive_temperatures_3d( &
      species, updated_state, stage_temperature, nx, ny, nz, &
      updated_temperature, local_ok)
    if (.not. local_ok) return

    flux_x = 0.5_dp * (first_flux_x + second_flux_x)
    flux_y = 0.5_dp * (first_flux_y + second_flux_y)
    flux_z = 0.5_dp * (first_flux_z + second_flux_z)
    state = updated_state
    temperature = updated_temperature
    ok = .true.
  end subroutine advance_reactive_euler_ssprk2_with_fluxes_3d

  subroutine advance_reactive_euler_ssprk2_plm_3d( &
      species, state, temperature, nx, ny, nz, dx, dy, dz, dt, &
      limiter, riemann_solver, ok, difmag)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz, dt
    character(len=*), intent(in) :: limiter, riemann_solver
    real(dp), optional, intent(in) :: difmag
    logical, intent(out) :: ok

    real(dp), allocatable :: flux_x(:, :, :, :)
    real(dp), allocatable :: flux_y(:, :, :, :)
    real(dp), allocatable :: flux_z(:, :, :, :)
    integer :: nvar

    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return
    if (nx < 3 .or. ny < 3 .or. nz < 3) return
    allocate(flux_x(nvar, nx, ny, nz))
    allocate(flux_y(nvar, nx, ny, nz))
    allocate(flux_z(nvar, nx, ny, nz))
    call advance_reactive_euler_ssprk2_plm_with_fluxes_3d( &
      species, state, temperature, nx, ny, nz, dx, dy, dz, dt, &
      limiter, riemann_solver, flux_x, flux_y, flux_z, ok, difmag=difmag)
  end subroutine advance_reactive_euler_ssprk2_plm_3d

  subroutine advance_reactive_euler_ssprk2_plm_with_fluxes_3d( &
      species, state, temperature, nx, ny, nz, dx, dy, dz, dt, &
      limiter, riemann_solver, flux_x, flux_y, flux_z, ok, difmag)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz, dt
    character(len=*), intent(in) :: limiter, riemann_solver
    real(dp), intent(out) :: flux_x(:, :, :, :)
    real(dp), intent(out) :: flux_y(:, :, :, :)
    real(dp), intent(out) :: flux_z(:, :, :, :)
    real(dp), optional, intent(in) :: difmag
    logical, intent(out) :: ok

    real(dp), allocatable :: old_state(:, :, :, :)
    real(dp), allocatable :: stage_state(:, :, :, :)
    real(dp), allocatable :: updated_state(:, :, :, :)
    real(dp), allocatable :: rhs(:, :, :, :)
    real(dp), allocatable :: first_flux_x(:, :, :, :)
    real(dp), allocatable :: first_flux_y(:, :, :, :)
    real(dp), allocatable :: first_flux_z(:, :, :, :)
    real(dp), allocatable :: second_flux_x(:, :, :, :)
    real(dp), allocatable :: second_flux_y(:, :, :, :)
    real(dp), allocatable :: second_flux_z(:, :, :, :)
    real(dp), allocatable :: old_temperature(:, :, :)
    real(dp), allocatable :: stage_temperature(:, :, :)
    real(dp), allocatable :: updated_temperature(:, :, :)
    logical :: local_ok
    real(dp) :: coefficient, boundary_change(size(state,1))
    character(len=8) :: faces(6)
    integer :: nvar

    flux_x = 0.0_dp
    flux_y = 0.0_dp
    flux_z = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return
    if (.not. ieee_is_finite(dt)) return
    if (nx < 3 .or. ny < 3 .or. nz < 3 .or. dt <= 0.0_dp) return
    if (any(shape(flux_x) /= shape(state)) .or. &
        any(shape(flux_y) /= shape(state)) .or. &
        any(shape(flux_z) /= shape(state))) return

    coefficient=.1_dp
    if (present(difmag)) coefficient=difmag
    if (.not.ieee_is_finite(coefficient).or.coefficient<0) return
    if (coefficient>0) then
      faces='periodic'
      call advance_reactive_boundary_hydro_3d(species,state,temperature,[dx,dy,dz],dt, &
        faces,'characteristic_plm',limiter,riemann_solver,boundary_change,ok, flux_x=flux_x,flux_y=flux_y,flux_z=flux_z, &
        difmag=coefficient)
      return
    end if

    allocate(old_state, source=state)
    allocate(stage_state(nvar, nx, ny, nz))
    allocate(updated_state(nvar, nx, ny, nz))
    allocate(rhs(nvar, nx, ny, nz))
    allocate(first_flux_x(nvar, nx, ny, nz))
    allocate(first_flux_y(nvar, nx, ny, nz))
    allocate(first_flux_z(nvar, nx, ny, nz))
    allocate(second_flux_x(nvar, nx, ny, nz))
    allocate(second_flux_y(nvar, nx, ny, nz))
    allocate(second_flux_z(nvar, nx, ny, nz))
    allocate(old_temperature(nx, ny, nz))
    allocate(stage_temperature(nx, ny, nz))
    allocate(updated_temperature(nx, ny, nz))
    call recover_reactive_temperatures_3d( &
      species, old_state, temperature, nx, ny, nz, old_temperature, local_ok)
    if (.not. local_ok) return

    call compute_reactive_plm_face_fluxes_3d( &
      species, old_state, old_temperature, nx, ny, nz, limiter, &
      riemann_solver, first_flux_x, first_flux_y, first_flux_z, local_ok)
    if (.not. local_ok) return
    call reactive_flux_divergence_3d( &
      first_flux_x, first_flux_y, first_flux_z, nx, ny, nz, &
      dx, dy, dz, rhs, local_ok)
    if (.not. local_ok) return
    stage_state = old_state + dt * rhs
    call recover_reactive_temperatures_3d( &
      species, stage_state, old_temperature, nx, ny, nz, &
      stage_temperature, local_ok)
    if (.not. local_ok) return

    call compute_reactive_plm_face_fluxes_3d( &
      species, stage_state, stage_temperature, nx, ny, nz, limiter, &
      riemann_solver, second_flux_x, second_flux_y, second_flux_z, local_ok)
    if (.not. local_ok) return
    call reactive_flux_divergence_3d( &
      second_flux_x, second_flux_y, second_flux_z, nx, ny, nz, &
      dx, dy, dz, rhs, local_ok)
    if (.not. local_ok) return
    updated_state = 0.5_dp * old_state + &
      0.5_dp * (stage_state + dt * rhs)
    call recover_reactive_temperatures_3d( &
      species, updated_state, stage_temperature, nx, ny, nz, &
      updated_temperature, local_ok)
    if (.not. local_ok) return

    flux_x = 0.5_dp * (first_flux_x + second_flux_x)
    flux_y = 0.5_dp * (first_flux_y + second_flux_y)
    flux_z = 0.5_dp * (first_flux_z + second_flux_z)
    state = updated_state
    temperature = updated_temperature
    ok = .true.
  end subroutine advance_reactive_euler_ssprk2_plm_with_fluxes_3d

  subroutine advance_reactive_chemistry_3d( &
      species, reactions, state, temperature, nx, ny, nz, interval, &
      rtol, atol, ok, active_mask, chemistry_integrator)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: interval, rtol, atol
    logical, intent(out) :: ok
    logical, intent(in), optional :: active_mask(:, :, :)
    character(len=*), intent(in), optional :: chemistry_integrator

    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    logical :: plane_ok
    integer :: plane, nvar

    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return
    if (present(active_mask)) then
      if (any(shape(active_mask) /= [nx, ny, nz])) return
    end if
    if (.not. all(ieee_is_finite([interval, rtol, atol]))) return
    if (nx < 1 .or. ny < 1 .or. nz < 1 .or. interval < 0.0_dp .or. &
        rtol <= 0.0_dp .or. atol <= 0.0_dp) return
    if (interval <= tiny(1.0_dp)) then
      ok = .true.
      return
    end if

    allocate(candidate_state, source=state)
    allocate(candidate_temperature, source=temperature)
    do plane = 1, nz
      if (present(active_mask)) then
        call advance_reactive_chemistry_2d( &
          species, reactions, candidate_state(:, :, :, plane), &
          candidate_temperature(:, :, plane), nx, ny, interval, rtol, atol, &
          plane_ok, active_mask(:, :, plane), &
          chemistry_integrator=chemistry_integrator)
      else
        call advance_reactive_chemistry_2d( &
          species, reactions, candidate_state(:, :, :, plane), &
          candidate_temperature(:, :, plane), nx, ny, interval, rtol, atol, &
          plane_ok, chemistry_integrator=chemistry_integrator)
      end if
      if (.not. plane_ok) return
    end do
    state = candidate_state
    temperature = candidate_temperature
    ok = .true.
  end subroutine advance_reactive_chemistry_3d

  subroutine advance_reactive_strang_3d( &
      species, reactions, state, temperature, nx, ny, nz, dx, dy, dz, dt, &
      riemann_solver, chemistry_enabled, rtol, atol, ok, chemistry_integrator, difmag)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz, dt, rtol, atol
    character(len=*), intent(in) :: riemann_solver
    logical, intent(in) :: chemistry_enabled
    real(dp), optional, intent(in) :: difmag
    logical, intent(out) :: ok
    character(len=*), intent(in), optional :: chemistry_integrator

    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    logical :: local_ok
    integer :: nvar

    ok = .false.
    if (present(difmag)) then
      if (.not.ieee_is_finite(difmag).or.difmag<0) return
    end if
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return
    if (.not. ieee_is_finite(dt)) return
    if (dt <= 0.0_dp) return
    if (chemistry_enabled) then
      if (.not. all(ieee_is_finite([rtol, atol]))) return
      if (rtol <= 0.0_dp .or. atol <= 0.0_dp) return
    end if

    allocate(candidate_state, source=state)
    allocate(candidate_temperature, source=temperature)
    if (chemistry_enabled) then
      call advance_reactive_chemistry_3d( &
        species, reactions, candidate_state, candidate_temperature, &
        nx, ny, nz, 0.5_dp * dt, rtol, atol, local_ok, &
        chemistry_integrator=chemistry_integrator)
      if (.not. local_ok) return
    end if
    call advance_reactive_euler_ssprk2_3d( &
      species, candidate_state, candidate_temperature, nx, ny, nz, &
      dx, dy, dz, dt, riemann_solver, local_ok, difmag=difmag)
    if (.not. local_ok) return
    if (chemistry_enabled) then
      call advance_reactive_chemistry_3d( &
        species, reactions, candidate_state, candidate_temperature, &
        nx, ny, nz, 0.5_dp * dt, rtol, atol, local_ok, &
        chemistry_integrator=chemistry_integrator)
      if (.not. local_ok) return
    end if
    state = candidate_state
    temperature = candidate_temperature
    ok = .true.
  end subroutine advance_reactive_strang_3d

  subroutine advance_reactive_full_3d( &
      species, reactions, transport, state, temperature, nx, ny, nz, &
      dx, dy, dz, dt, riemann_solver, chemistry_enabled, rtol, atol, &
      transport_enabled, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, &
      minimum_transport_theta, ok, reconstruction, limiter, chemistry_integrator, boundary_faces, &
      boundary_change, inflow_state, inflow_temperature, wall_temperature,wall_velocity,outlet_pressure, &
      inflow_mach,inflow_pressure,difmag)
    real(dp), optional, intent(in) :: difmag
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz, dt, rtol, atol
    character(len=*), intent(in) :: riemann_solver
    logical, intent(in) :: chemistry_enabled, transport_enabled
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: minimum_transport_theta
    logical, intent(out) :: ok
    character(len=*), intent(in), optional :: reconstruction, limiter
    character(len=*), intent(in), optional :: chemistry_integrator

    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    real(dp), intent(in), optional :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6)
    character(len=*), intent(in), optional :: boundary_faces(6)
    real(dp), intent(out), optional :: boundary_change(:)
    character(len=32) :: faces(6)
    real(dp) :: change(size(state,1)),hydro_change(size(state,1)),transport_change(size(state,1))
    real(dp) :: stage_transport_theta,transverse_theta,explicit_viscosity
    logical :: local_ok
    integer :: nvar
    character(len=32) :: selected_reconstruction, selected_limiter

    ok = .false.
    minimum_transport_theta = 1.0_dp
    faces = "periodic"
    if (present(boundary_faces)) faces = boundary_faces
    if(present(inflow_pressure).and..not.present(inflow_mach)) return
    if(present(inflow_mach)) then
      if(.not.all(ieee_is_finite(inflow_mach))) return
      if(any(inflow_mach>100).or.any(inflow_mach>=0.and.faces/='inflow')) return
      if(present(inflow_pressure)) then
        if(.not.all(ieee_is_finite(inflow_pressure))) return
        if(any(inflow_pressure<0).or.any(inflow_pressure>1e15_dp)) return
        if(any(inflow_pressure>0.and.inflow_mach<0)) return
      end if
    end if
    if(present(outlet_pressure)) then
      if(.not.all(ieee_is_finite(outlet_pressure))) return
      if(any(outlet_pressure<0).or.any(outlet_pressure>1e15_dp)) return
      if(any(outlet_pressure>0.and.faces/='outflow')) return
    end if
    if (present(boundary_change)) then
      boundary_change = 0.0_dp
      if (size(boundary_change) /= size(state,1)) return
    end if
    change = 0.0_dp
    hydro_change = 0.0_dp
    if (.not. valid_reactive_boundaries_3d(faces)) return
    if (any(faces=='isothermal_wall')) then
      if (.not.transport_enabled.or..not.thermal_conduction_enabled) return
    end if
    selected_reconstruction = "pcm"
    selected_limiter = "mc"
    if (present(reconstruction)) selected_reconstruction = &
      trim(reconstruction)
    if (present(limiter)) selected_limiter = trim(limiter)
    explicit_viscosity=.1_dp
    if (present(difmag)) then
      if (.not.ieee_is_finite(difmag).or.difmag<0) return
      explicit_viscosity=difmag
    end if
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return
    if (.not. ieee_is_finite(dt)) return
    if (dt <= 0.0_dp) return
    if (chemistry_enabled) then
      if (.not. all(ieee_is_finite([rtol, atol]))) return
      if (rtol <= 0.0_dp .or. atol <= 0.0_dp) return
    end if
    if (transport_enabled .and. size(transport) /= size(species)) return
    if (trim(selected_reconstruction) /= "pcm" .and. &
        trim(selected_reconstruction) /= "characteristic_plm" .and. &
        trim(selected_reconstruction) /= "characteristic_ppm") return
    if (trim(selected_limiter) /= "minmod" .and. &
        trim(selected_limiter) /= "mc") return

    allocate(candidate_state, source=state)
    allocate(candidate_temperature, source=temperature)
    if (chemistry_enabled) then
      call advance_reactive_chemistry_3d( &
        species, reactions, candidate_state, candidate_temperature, &
        nx, ny, nz, 0.5_dp * dt, rtol, atol, local_ok, &
        chemistry_integrator=chemistry_integrator)
      if (.not. local_ok) return
    end if
    if (transport_enabled) then
      call advance_reactive_transport_3d( &
        species, transport, candidate_state, candidate_temperature, &
        nx, ny, nz, dx, dy, dz, 0.5_dp * dt, viscosity_enabled, &
        thermal_conduction_enabled, species_diffusion_enabled, &
        barodiffusion_enabled, stage_transport_theta, local_ok, faces, transport_change, &
        inflow_state, inflow_temperature, wall_temperature,wall_velocity,outlet_pressure, &
        inflow_mach,inflow_pressure)
      if (.not. local_ok) return
      change = change + transport_change
      minimum_transport_theta = min( &
        minimum_transport_theta, stage_transport_theta)
    end if
    if (any(faces /= "periodic").or. &
        (explicit_viscosity>0.and.selected_reconstruction/="characteristic_ppm")) then
      call advance_reactive_boundary_hydro_3d(species, candidate_state, candidate_temperature, &
        [dx,dy,dz], dt, faces, selected_reconstruction, selected_limiter, riemann_solver, hydro_change, local_ok, &
        inflow_state, inflow_temperature,outlet_pressure=outlet_pressure, &
        inflow_mach=inflow_mach,inflow_pressure=inflow_pressure,difmag=difmag)
    else
      select case (trim(selected_reconstruction))
      case ("pcm")
        call advance_reactive_euler_ssprk2_3d( &
          species, candidate_state, candidate_temperature, nx, ny, nz, &
          dx, dy, dz, dt, riemann_solver, local_ok, difmag=explicit_viscosity)
      case ("characteristic_ppm")
        call advance_reactive_euler_ppm_ctu_3d(species,candidate_state,candidate_temperature,[dx,dy,dz],dt, &
          .false.,.true.,riemann_solver,transverse_theta,local_ok,difmag=difmag)
      case ("characteristic_plm")
        call advance_reactive_euler_ssprk2_plm_3d( &
          species, candidate_state, candidate_temperature, nx, ny, nz, &
          dx, dy, dz, dt, selected_limiter, riemann_solver, local_ok, difmag=explicit_viscosity)
      end select
    end if
    if (.not. local_ok) return
    change = change + hydro_change
    if (transport_enabled) then
      call advance_reactive_transport_3d( &
        species, transport, candidate_state, candidate_temperature, &
        nx, ny, nz, dx, dy, dz, 0.5_dp * dt, viscosity_enabled, &
        thermal_conduction_enabled, species_diffusion_enabled, &
        barodiffusion_enabled, stage_transport_theta, local_ok, faces, transport_change, &
        inflow_state, inflow_temperature, wall_temperature,wall_velocity,outlet_pressure, &
        inflow_mach,inflow_pressure)
      if (.not. local_ok) return
      change = change + transport_change
      minimum_transport_theta = min( &
        minimum_transport_theta, stage_transport_theta)
    end if
    if (chemistry_enabled) then
      call advance_reactive_chemistry_3d( &
        species, reactions, candidate_state, candidate_temperature, &
        nx, ny, nz, 0.5_dp * dt, rtol, atol, local_ok, &
        chemistry_integrator=chemistry_integrator)
      if (.not. local_ok) return
    end if
    if (present(boundary_change)) boundary_change = change
    state = candidate_state
    temperature = candidate_temperature
    ok = .true.
  end subroutine advance_reactive_full_3d


  ! Physical face ordering: x-low, x-high, y-low, y-high, z-low, z-high.
  ! PPM uses four supplied layers and one CTU update. For PCM/PLM,
  ! ghosts are reconstructed at each RK stage. Two layers contain the entire
  ! PLM stencil at the physical face, so periodic wrapping outside the padded
  ! array cannot enter a physical-cell flux. Only physical cells are advanced.
  subroutine advance_reactive_boundary_hydro_3d(species,state,temperature,spacing,dt, &
      faces,reconstruction,limiter,riemann_solver,boundary_change,ok,inflow_state,inflow_temperature, &
      flux_x,flux_y,flux_z,outlet_pressure,inflow_mach,inflow_pressure,difmag)
    real(dp), optional, intent(in) :: difmag
    real(dp), intent(in), optional :: inflow_state(:,:),inflow_temperature(:),outlet_pressure(6)
    real(dp), intent(in), optional :: inflow_mach(6),inflow_pressure(6)
    ! Optional time-centered fluxes: upper faces, or all faces including index zero.
    real(dp), optional, intent(out) :: flux_x(:,:,:,:),flux_y(:,:,:,:),flux_z(:,:,:,:)
    real(dp), allocatable :: mean_x(:,:,:,:),mean_y(:,:,:,:),mean_z(:,:,:,:)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    real(dp), intent(in) :: spacing(3),dt
    character(len=*), intent(in) :: faces(6),reconstruction,limiter,riemann_solver
    real(dp), intent(out) :: boundary_change(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: candidate(:,:,:,:),temp(:,:,:),recovered(:,:,:),rhs(:,:,:,:)
    real(dp) :: rate0(size(state,1)),rate1(size(state,1)),explicit_viscosity
    integer :: n(3)
    logical :: local_ok,full_faces
    ok=.false.; boundary_change=0.0_dp
    if (present(flux_x)) flux_x=0
    if (present(flux_y)) flux_y=0
    if (present(flux_z)) flux_z=0
    explicit_viscosity=.1_dp
    if (present(difmag)) then
      if (.not.ieee_is_finite(difmag).or.difmag<0) return
      explicit_viscosity=difmag
    end if
    if ((present(flux_x).neqv.present(flux_y)).or.(present(flux_x).neqv.present(flux_z))) return
    if(present(outlet_pressure)) then
      if(.not.all(ieee_is_finite(outlet_pressure))) return
      if(any(outlet_pressure<0).or.any(outlet_pressure>1e15_dp)) return
      if(any(outlet_pressure>0.and.faces/='outflow')) return
    end if
    if(present(inflow_pressure).and..not.present(inflow_mach)) return
    full_faces=.false.
    n=shape(temperature)
    if (present(flux_x)) then
      full_faces=all(shape(flux_x)==[size(state,1),n(1)+1,n(2),n(3)]).and. &
        all(shape(flux_y)==[size(state,1),n(1),n(2)+1,n(3)]).and. &
        all(shape(flux_z)==[size(state,1),n(1),n(2),n(3)+1])
      if (.not.full_faces) then
        if (any(shape(flux_x)/=shape(state)).or.any(shape(flux_y)/=shape(state)).or.any(shape(flux_z)/=shape(state))) return
      end if
      allocate(mean_x,mold=flux_x);allocate(mean_y,mold=flux_y);allocate(mean_z,mold=flux_z)
      mean_x=0;mean_y=0;mean_z=0
    end if
    n=shape(temperature)
    if (size(boundary_change) /= size(state,1)) return
    if (.not. valid_reactive_3d_shapes(state,temperature,reactive_nvar(size(species)),n(1),n(2),n(3))) return
    if (any(n<2) .or. .not. valid_reactive_boundaries_3d(faces)) return
    if (.not. all(ieee_is_finite([spacing,dt]))) return
    if (minval(spacing)<=0.0_dp .or. dt<=0.0_dp) return
    if (reconstruction/="pcm" .and. reconstruction/="characteristic_plm" .and. &
        reconstruction/="characteristic_ppm") return
    if (limiter/="mc" .and. limiter/="minmod") return
    if (reconstruction=="characteristic_ppm") then
      call boundary_ppm(local_ok)
      ok=local_ok
      return
    end if
    allocate(candidate,source=state); allocate(temp,mold=temperature)
    allocate(recovered,mold=temperature); allocate(rhs,mold=state)
    call recover_reactive_temperatures_3d(species,state,temperature,n(1),n(2),n(3),temp,local_ok)
    if (.not.local_ok) return
    call boundary_rhs(state,temp,rhs,rate0,local_ok)
    if (.not.local_ok) return
    candidate=state+dt*rhs
    call recover_reactive_temperatures_3d(species,candidate,temp,n(1),n(2),n(3),recovered,local_ok)
    if (.not.local_ok) return
    temp=recovered
    call boundary_rhs(candidate,temp,rhs,rate1,local_ok)
    if (.not.local_ok) return
    candidate=0.5_dp*state+0.5_dp*(candidate+dt*rhs)
    call recover_reactive_temperatures_3d(species,candidate,temp,n(1),n(2),n(3),recovered,local_ok)
    if (.not.local_ok) return
    boundary_change=0.5_dp*dt*(rate0+rate1)
    if (present(flux_x)) then
      flux_x=mean_x;flux_y=mean_y;flux_z=mean_z
    end if
    state=candidate; temperature=recovered; ok=.true.
  contains
    subroutine boundary_ppm(success)
      logical, intent(out) :: success
      real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:)
      real(dp) :: theta,change(size(state,1))
      integer :: i,j,k,nv
      success=.false.;nv=size(state,1)
      allocate(ghost(nv,n(1)+8,n(2)+8,n(3)+8),gt(n(1)+8,n(2)+8,n(3)+8))
      allocate(fx(nv,0:n(1),n(2),n(3)),fy(nv,n(1),0:n(2),n(3)),fz(nv,n(1),n(2),0:n(3)))
      if (present(inflow_mach)) then
        call fill_mach_inlet_boundary_ghosts_3d(species,state,temperature,faces,inflow_mach,ghost,gt,success, &
          inflow_state=inflow_state,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,ghost_width=4)
      else
        call fill_reactive_boundary_ghosts_3d(state,temperature,faces,ghost,gt,success, &
          inflow_state=inflow_state,inflow_temperature=inflow_temperature,ghost_width=4)
      end if
      if (.not.success) return
      call fill_pressure_outlet_ghosts_3d(species,n,faces,outlet_pressure,ghost,gt,success,ghost_width=4)
      if (.not.success) return
      call compute_reactive_ppm_halo_ctu_fluxes_3d(species,ghost,gt,n,spacing,dt,.false.,.true., &
        riemann_solver,fx,fy,fz,theta,success,inflow_faces=faces=='inflow', &
        no_slip_faces=faces=='no_slip_wall'.or.faces=='isothermal_wall',difmag=difmag, &
        reflecting_faces=faces=='slip_wall'.or.faces=='no_slip_wall'.or.faces=='isothermal_wall', &
        outflow_faces=faces=='outflow',outlet_pressure=outlet_pressure)
      if (.not.success) return
      success=.false.
      candidate=state
      allocate(recovered,mold=temperature)
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        candidate(:,i,j,k)=state(:,i,j,k)-(dt/spacing(1))*(fx(:,i,j,k)-fx(:,i-1,j,k))- &
          (dt/spacing(2))*(fy(:,i,j,k)-fy(:,i,j-1,k))-(dt/spacing(3))*(fz(:,i,j,k)-fz(:,i,j,k-1))
      end do
      end do
      end do
      call recover_reactive_temperatures_3d(species,candidate,temperature,n(1),n(2),n(3),recovered,success)
      if (.not.success) return
      change=0
      do k=1,n(3)
      do j=1,n(2)
        change=change+(fx(:,0,j,k)-fx(:,n(1),j,k))*spacing(2)*spacing(3)
      end do
      end do
      do k=1,n(3)
      do i=1,n(1)
        change=change+(fy(:,i,0,k)-fy(:,i,n(2),k))*spacing(1)*spacing(3)
      end do
      end do
      do j=1,n(2)
      do i=1,n(1)
        change=change+(fz(:,i,j,0)-fz(:,i,j,n(3)))*spacing(1)*spacing(2)
      end do
      end do
      change=dt*change
      success=all(ieee_is_finite(change))
      if (.not.success) return
      if (present(flux_x)) then
        if (full_faces) then
          flux_x=fx;flux_y=fy;flux_z=fz
        else
          flux_x=fx(:,1:,:,:);flux_y=fy(:,:,1:,:);flux_z=fz(:,:,:,1:)
        end if
      end if
      state=candidate;temperature=recovered;boundary_change=change
    end subroutine boundary_ppm
    subroutine boundary_rhs(u,t,du,rate,success)
      real(dp), intent(in) :: u(:,:,:,:),t(:,:,:)
      real(dp), intent(out) :: du(:,:,:,:),rate(:)
      logical, intent(out) :: success
      real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:)
      integer :: dims(3),i,j,k,nv
      success=.false.; du=0.0_dp; rate=0.0_dp
      dims=n+4; nv=size(u,1)
      allocate(ghost(nv,dims(1),dims(2),dims(3)),gt(dims(1),dims(2),dims(3)))
      allocate(fx,mold=ghost); allocate(fy,mold=ghost); allocate(fz,mold=ghost)
      if(present(inflow_mach)) then
        call fill_mach_inlet_boundary_ghosts_3d(species,u,t,faces,inflow_mach,ghost,gt,success, &
          inflow_state=inflow_state,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure)
      else
      call fill_reactive_boundary_ghosts_3d(u,t,faces,ghost,gt,success, &
        inflow_state=inflow_state,inflow_temperature=inflow_temperature)
      end if
      if (.not.success) return
      call fill_pressure_outlet_ghosts_3d(species,n,faces,outlet_pressure,ghost,gt,success)
      if(.not.success) return
      if (reconstruction=='pcm') then
        call compute_reactive_face_fluxes_3d(species,ghost,gt,dims(1),dims(2),dims(3), &
          riemann_solver,fx,fy,fz,success)
      else
        call compute_reactive_plm_face_fluxes_3d(species,ghost,gt,dims(1),dims(2),dims(3), &
          limiter,riemann_solver,fx,fy,fz,success)
      end if
      if (.not.success) return
      if(explicit_viscosity>0) then
        call adjust_reactive_artificial_viscosity_3d(ghost(:,2:n(1)+3,2:n(2)+3,2:n(3)+3), &
          n,spacing,explicit_viscosity,faces=='no_slip_wall'.or.faces=='isothermal_wall', &
          fx(:,2:n(1)+2,3:n(2)+2,3:n(3)+2),fy(:,3:n(1)+2,2:n(2)+2,3:n(3)+2), &
          fz(:,3:n(1)+2,3:n(2)+2,2:n(3)+2),success)
        if(.not.success) return
      end if
      do k=3,n(3)+2
        do j=3,n(2)+2
          do i=3,n(1)+2
            du(:,i-2,j-2,k-2)=-(fx(:,i,j,k)-fx(:,i-1,j,k))/spacing(1) &
              -(fy(:,i,j,k)-fy(:,i,j-1,k))/spacing(2)-(fz(:,i,j,k)-fz(:,i,j,k-1))/spacing(3)
          end do
        end do
      end do
      if (present(flux_x)) then
        if (full_faces) then
          mean_x=mean_x+.5_dp*fx(:,2:n(1)+2,3:n(2)+2,3:n(3)+2)
          mean_y=mean_y+.5_dp*fy(:,3:n(1)+2,2:n(2)+2,3:n(3)+2)
          mean_z=mean_z+.5_dp*fz(:,3:n(1)+2,3:n(2)+2,2:n(3)+2)
        else
        mean_x=mean_x+.5_dp*fx(:,3:n(1)+2,3:n(2)+2,3:n(3)+2)
        mean_y=mean_y+.5_dp*fy(:,3:n(1)+2,3:n(2)+2,3:n(3)+2)
        mean_z=mean_z+.5_dp*fz(:,3:n(1)+2,3:n(2)+2,3:n(3)+2)
        end if
      end if
      ! Surface integral, independent of volume summation used by diagnostics.
      do k=3,n(3)+2
        do j=3,n(2)+2
          rate=rate+(fx(:,2,j,k)-fx(:,n(1)+2,j,k))*spacing(2)*spacing(3)
        end do
      end do
      do k=3,n(3)+2
        do i=3,n(1)+2
          rate=rate+(fy(:,i,2,k)-fy(:,i,n(2)+2,k))*spacing(1)*spacing(3)
        end do
      end do
      do j=3,n(2)+2
        do i=3,n(1)+2
          rate=rate+(fz(:,i,j,2)-fz(:,i,j,n(3)+2))*spacing(1)*spacing(2)
        end do
      end do
      success=all(ieee_is_finite(du)).and.all(ieee_is_finite(rate))
    end subroutine boundary_rhs
  end subroutine advance_reactive_boundary_hydro_3d

  subroutine recover_reactive_temperatures_3d( &
      species, state, temperature_guess, nx, ny, nz, temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature_guess(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(out) :: temperature(:, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: sound_speed
    logical :: cell_ok
    integer :: i, j, k, nvar

    temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature_guess, nvar, nx, ny, nz)) return
    if (size(temperature, 1) /= nx .or. size(temperature, 2) /= ny .or. &
        size(temperature, 3) /= nz) return

    allocate(primitive(reactive_nprim(size(species))))
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature_guess(i, j, k), &
            primitive, temperature(i, j, k), sound_speed, cell_ok)
          if (.not. cell_ok) then
            temperature = 0.0_dp
            return
          end if
        end do
      end do
    end do
    ok = all(ieee_is_finite(temperature))
    if (.not. ok) then
      temperature = 0.0_dp
      return
    end if
    ok = minval(temperature) > 0.0_dp
    if (.not. ok) temperature = 0.0_dp
  end subroutine recover_reactive_temperatures_3d

  subroutine reactive_integrals_3d( &
      state, nx, ny, nz, dx, dy, dz, integrals, ok)
    real(dp), intent(in) :: state(:, :, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz
    real(dp), intent(out) :: integrals(:)
    logical, intent(out) :: ok

    integer :: component

    integrals = 0.0_dp
    ok = .false.
    if (.not. all(ieee_is_finite([dx, dy, dz]))) return
    if (.not. all(ieee_is_finite(state))) return
    ok = size(state, 1) == size(integrals) .and. &
      size(state, 2) == nx .and. size(state, 3) == ny .and. &
      size(state, 4) == nz .and. nx > 0 .and. ny > 0 .and. nz > 0 .and. &
      dx > 0.0_dp .and. dy > 0.0_dp .and. dz > 0.0_dp
    if (.not. ok) return
    do component = 1, size(integrals)
      integrals(component) = sum(state(component, :, :, :)) * dx * dy * dz
    end do
    ok = all(ieee_is_finite(integrals))
    if (.not. ok) integrals = 0.0_dp
  end subroutine reactive_integrals_3d

  subroutine reactive_extrema_3d( &
      species, state, temperature, nx, ny, nz, &
      minimum_density, maximum_density, minimum_pressure, maximum_pressure, &
      minimum_temperature, maximum_temperature, maximum_speed, &
      maximum_closure_error, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(out) :: minimum_density, maximum_density
    real(dp), intent(out) :: minimum_pressure, maximum_pressure
    real(dp), intent(out) :: minimum_temperature, maximum_temperature
    real(dp), intent(out) :: maximum_speed, maximum_closure_error
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: local_temperature, sound_speed, closure, speed
    logical :: cell_ok
    integer :: i, j, k, species_index, nvar

    minimum_density = huge(1.0_dp)
    maximum_density = -huge(1.0_dp)
    minimum_pressure = huge(1.0_dp)
    maximum_pressure = -huge(1.0_dp)
    minimum_temperature = huge(1.0_dp)
    maximum_temperature = -huge(1.0_dp)
    maximum_speed = 0.0_dp
    maximum_closure_error = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. valid_reactive_3d_shapes( &
          state, temperature, nvar, nx, ny, nz)) return

    allocate(primitive(reactive_nprim(size(species))))
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), primitive, &
            local_temperature, sound_speed, cell_ok)
          if (.not. cell_ok) return
          minimum_density = min(minimum_density, primitive(1))
          maximum_density = max(maximum_density, primitive(1))
          minimum_pressure = min(minimum_pressure, primitive(5))
          maximum_pressure = max(maximum_pressure, primitive(5))
          minimum_temperature = min(minimum_temperature, local_temperature)
          maximum_temperature = max(maximum_temperature, local_temperature)
          speed = sqrt(sum(primitive(2:4)**2))
          maximum_speed = max(maximum_speed, speed)
          closure = 0.0_dp
          do species_index = 1, size(species)
            closure = closure + primitive( &
              reactive_mass_fraction_component(species_index))
          end do
          maximum_closure_error = max( &
            maximum_closure_error, abs(closure - 1.0_dp))
        end do
      end do
    end do
    ok = all(ieee_is_finite([ &
      minimum_density, maximum_density, minimum_pressure, maximum_pressure, &
      minimum_temperature, maximum_temperature, maximum_speed, &
      maximum_closure_error]))
  end subroutine reactive_extrema_3d

  pure logical function valid_reactive_3d_shapes( &
      state, temperature, nvar, nx, ny, nz) result(valid)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nvar, nx, ny, nz

    valid = nvar > 0 .and. size(state, 1) == nvar .and. &
      size(state, 2) == nx .and. size(state, 3) == ny .and. &
      size(state, 4) == nz .and. size(temperature, 1) == nx .and. &
      size(temperature, 2) == ny .and. size(temperature, 3) == nz
  end function valid_reactive_3d_shapes

end module reactive_3d_mod
