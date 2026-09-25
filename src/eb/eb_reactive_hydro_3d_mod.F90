module eb_reactive_hydro_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use reactive_3d_mod, only: adjust_reactive_artificial_viscosity_3d
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive, &
    reactive_riemann_flux_x,reconstruct_pelec_plm_cell,compute_pelec_plm_slopes,reactive_pelec_primitive_flux_x, &
    reactive_ppm_flattening_coefficient
  use reactive_2d_mod, only: reactive_riemann_flux_y
  use reactive_directional_flux_3d_mod, only: reactive_riemann_flux_z
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d,extract_eb_geometry_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
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
  end type eb_plm_face_prediction_3d
  public :: compute_reactive_eb_plm_normal_fluxes_3d,apply_reactive_eb_transverse_3d
  public :: compute_reactive_eb_plm_ctu_fluxes_3d,advance_reactive_eb_plm_ctu_3d
  public :: advance_reactive_eb_outflow_ctu_3d
  public :: compute_reactive_eb_flattening_3d
  public :: predict_reactive_eb_plm_cell_3d
  public :: reactive_eb_outflow_riemann_fluxes_3d,reactive_eb_riemann_fluxes_3d
  public :: advance_reactive_eb_euler_3d
  public :: advance_reactive_eb_redistributed_euler_3d
  public :: advance_reactive_eb_state_redistributed_euler_3d

  interface reactive_eb_riemann_fluxes_3d
    module procedure reactive_eb_outflow_riemann_fluxes_3d
  end interface

contains

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
    type(eb_geometry_3d) :: working
    type(eb_plm_boundary_3d) :: boundary
    real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),vf(:,:,:),flat(:,:,:),candidate(:,:,:,:),ct(:,:,:)
    integer :: n(3),nv,i,j,k,donor(3)
    logical :: valid,use_flattening
    ok=.false.;nv=reactive_nvar(size(species));n=shape(temperature)
    if(any(n<1).or.any(shape(state)/=[nv,n])) return
    if(any(shape(new_state)/=[nv,n]).or.any(shape(new_temperature)/=n)) return
    if(any([extended%nx,extended%ny,extended%nz]/=n+12)) return
    call extract_eb_geometry_3d(extended,[4,4,4],n+9,working,valid)
    if(.not.valid) return
    allocate(ghost(nv,-2:n(1)+9,-2:n(2)+9,-2:n(3)+9),gt(-2:n(1)+9,-2:n(2)+9,-2:n(3)+9), &
      vf(-2:n(1)+9,-2:n(2)+9,-2:n(3)+9),flat(0:n(1)+7,0:n(2)+7,0:n(3)+7))
    vf=extended%volume_fraction;ghost=0;gt=0;flat=1
    do k=-2,n(3)+9
    do j=-2,n(2)+9
    do i=-2,n(1)+9
      if(vf(i,j,k)<=0) cycle
      donor=max(1,min(n,[i,j,k]-3))
      ! Do not invent a fluid state if the supplied geometry needs a covered donor.
      if(extended%volume_fraction(donor(1)+6,donor(2)+6,donor(3)+6)<=0) return
      ghost(:,i,j,k)=state(:,donor(1),donor(2),donor(3))
      gt(i,j,k)=temperature(donor(1),donor(2),donor(3))
    end do
    end do
    end do
    use_flattening=.true.
    if(present(shock_flattening)) use_flattening=shock_flattening
    if(use_flattening) then
      call compute_reactive_eb_flattening_3d(species,ghost,gt,vf,n+6,flat,valid)
      if(.not.valid) return
    end if
    boundary%lower=4;boundary%upper=n+3;boundary%faces='outflow'
    candidate=state;ct=temperature
    call advance_reactive_eb_plm_ctu_3d(species,working,ghost,gt,vf,dt,flat,candidate,ct,ok, &
      redistribution,difmag,integrated_flux,boundary_change,wall_change,target_volume_fraction,order,boundary)
    if(.not.ok) return
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
  end subroutine advance_reactive_eb_outflow_ctu_3d

  ! Complete hydro update for the central cells 4:N-3 of a prepared working
  ! geometry. Physical boundary filling/flattening remain caller responsibilities.
  subroutine advance_reactive_eb_plm_ctu_3d(species,geometry,ghost,temperature,volume_fraction,dt,flattening, &
      new_state,new_temperature,ok,redistribution,difmag,integrated_flux,boundary_change,wall_change, &
      target_volume_fraction,order,boundary)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_plm_boundary_3d), optional, intent(in) :: boundary
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: ghost(:,-2:,-2:,-2:),temperature(-2:,-2:,-2:),volume_fraction(-2:,-2:,-2:)
    real(dp), intent(in) :: dt,flattening(0:,0:,0:)
    real(dp), intent(inout) :: new_state(:,:,:,:),new_temperature(:,:,:)
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: redistribution
    real(dp), optional, intent(in) :: difmag,target_volume_fraction
    integer, optional, intent(in) :: order
    type(eb_integrated_flux_3d), optional, intent(inout) :: integrated_flux
    real(dp), optional, intent(inout) :: boundary_change(:),wall_change(:)
    type(eb_plm_face_prediction_3d) :: prediction(3)
    type(eb_geometry_3d) :: margin,physical
    type(eb_integrated_flux_3d) :: flux
    real(dp), allocatable :: state(:,:,:,:),temp(:,:,:),candidate(:,:,:,:),ct(:,:,:),rhs(:,:,:,:)
    real(dp) :: coefficient,target,total,recovered,speed,q(reactive_nvar(size(species)))
    real(dp) :: inventory_boundary(reactive_nvar(size(species))),wall(reactive_nvar(size(species)))
    integer :: n(3),m(3),nv,i,j,k,a,lo(4),hi(4),wall_lo(3),wall_hi(3)
    character(len=16) :: method
    logical :: valid,no_slip(6)
    ok=.false.;nv=reactive_nvar(size(species));n=[geometry%nx,geometry%ny,geometry%nz]
    if(any(n<7).or..not.geometry%is_valid()) return
    if(any(shape(new_state)/=[nv,n-6]).or.any(shape(new_temperature)/=n-6)) return
    if(.not.ieee_is_finite(dt).or.dt<=0) return
    coefficient=.1_dp;target=.5_dp;method='none'
    if(present(difmag)) coefficient=difmag
    if(present(target_volume_fraction)) target=target_volume_fraction
    if(present(redistribution)) method=redistribution
    if(.not.all(ieee_is_finite([coefficient,target])).or.coefficient<0.or.target<=0.or.target>1) return
    if(method/='none'.and.method/='flux_redist'.and.method/='state_redist') return
    if(present(boundary_change)) then
      if(size(boundary_change)/=nv) return
    end if
    if(present(wall_change)) then
      if(size(wall_change)/=nv) return
    end if
    call compute_reactive_eb_plm_ctu_fluxes_3d(species,geometry,ghost,temperature,volume_fraction, &
      dt,flattening,prediction,valid,order,boundary=boundary)
    if(.not.valid) return
    call extract_eb_geometry_3d(geometry,[3,3,3],n-2,margin,valid)
    if(.not.valid) return
    m=n-4
    no_slip=.false.;wall_lo=0;wall_hi=m
    if(present(boundary)) then
      no_slip=boundary%faces=='no_slip_wall'.or.boundary%faces=='isothermal_wall'
      wall_lo=boundary%lower-3;wall_hi=boundary%upper-2
    end if
    call adjust_reactive_artificial_viscosity_3d(ghost(:,2:n(1)-1,2:n(2)-1,2:n(3)-1), &
      m,[geometry%dx,geometry%dy,geometry%dz],coefficient,no_slip, &
      prediction(1)%flux,prediction(2)%flux,prediction(3)%flux,valid, &
      volume_fraction(2:n(1)-1,2:n(2)-1,2:n(3)-1), &
      margin%x_face_fraction>0,margin%y_face_fraction>0,margin%z_face_fraction>0,wall_lo,wall_hi)
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
    call hydro_inventory_changes(species,state,temp,physical,dt,flux,inventory_boundary,wall,valid)
    if(.not.valid) return
    allocate(candidate,mold=state);allocate(ct,mold=temp)
    select case(method)
    case('flux_redist')
      call advance_reactive_eb_redistributed_3d(species,state,temp,physical,rhs,dt,candidate,ct,valid)
    case('state_redist')
      call advance_reactive_eb_state_redistributed_3d(species,state,temp,physical,rhs,dt,candidate,ct,valid,target)
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
      dt,flattening,prediction,ok,order,source,source_energy,boundary)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: ghost(:,-2:,-2:,-2:),temperature(-2:,-2:,-2:),volume_fraction(-2:,-2:,-2:)
    real(dp), intent(in) :: dt,flattening(0:,0:,0:)
    type(eb_plm_face_prediction_3d), intent(inout) :: prediction(3)
    type(eb_plm_boundary_3d), optional, intent(in) :: boundary
    logical, intent(out) :: ok
    integer, optional, intent(in) :: order
    real(dp), optional, intent(in) :: source(:,:,:,:),source_energy(:,:,:)
    type(eb_plm_face_prediction_3d) :: normal(3),corner(3,3),final(3),candidate(3)
    real(dp), allocatable :: gamma(:,:,:),sound(:,:,:)
    real(dp) :: q(reactive_nvar(size(species))),recovered,speed,spacing(3)
    integer :: n(3),nv,i,j,k,axis,t,lo(3),hi(3),idx(3),other(2),m
    logical :: valid
    ok=.false.;n=[geometry%nx,geometry%ny,geometry%nz];nv=reactive_nvar(size(species))
    if(any(n<5)) return
    if(.not.valid_plm_boundary(boundary)) return
    if(present(source).neqv.present(source_energy)) return
    if(present(source)) then
      if(any(shape(source)/=[nv,n]).or.any(shape(source_energy)/=n)) return
    end if
    call compute_reactive_eb_plm_normal_fluxes_3d(species,geometry,ghost,temperature,volume_fraction, &
      dt,flattening,normal,valid,order,boundary)
    if(.not.valid) return
    spacing=[geometry%dx,geometry%dy,geometry%dz]
    allocate(gamma(n(1),n(2),n(3)),sound(n(1),n(2),n(3)));gamma=0;sound=0
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      if(geometry%volume_fraction(i,j,k)<=0) cycle
      call reactive_conserved_to_primitive(species,ghost(:,i,j,k),temperature(i,j,k),q,recovered,speed,valid)
      if(.not.valid) return
      gamma(i,j,k)=speed**2*q(1)/q(5);sound(i,j,k)=speed
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
      associate(c=>candidate(axis),f=>final(axis))
      idx=lo-2
      allocate(c%left(nv,idx(1):hi(1)-2,idx(2):hi(2)-2,idx(3):hi(3)-2), &
        c%right(nv,idx(1):hi(1)-2,idx(2):hi(2)-2,idx(3):hi(3)-2), &
        c%flux(nv,idx(1):hi(1)-2,idx(2):hi(2)-2,idx(3):hi(3)-2), &
        c%energy_left(idx(1):hi(1)-2,idx(2):hi(2)-2,idx(3):hi(3)-2), &
        c%energy_right(idx(1):hi(1)-2,idx(2):hi(2)-2,idx(3):hi(3)-2), &
        c%pressure(idx(1):hi(1)-2,idx(2):hi(2)-2,idx(3):hi(3)-2), &
        c%velocity(idx(1):hi(1)-2,idx(2):hi(2)-2,idx(3):hi(3)-2))
      c%left=f%left(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      c%right=f%right(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      c%flux=f%flux(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      c%energy_left=f%energy_left(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      c%energy_right=f%energy_right(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      c%pressure=f%pressure(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      c%velocity=f%velocity(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      end associate
    end do
    prediction=candidate;ok=.true.
  contains
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
        if(nt==2.and.present(source)) then
          call apply_reactive_eb_transverse_3d(states,energies,gamma(cell(1),cell(2),cell(3)),na,ta,tf,tp,tv, &
            coef,success,source(:,cell(1),cell(2),cell(3)),source_energy(cell(1),cell(2),cell(3)),.5_dp*dt)
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
      call reactive_pelec_primitive_flux_x(species,ql,qr, &
        .5_dp*(sound(index(1),index(2),index(3))+sound(cell(1),cell(2),cell(3))), &
        f,success,interface_pressure=pface,interface_velocity=uface)
      if(.not.success) return
      swap=f(2);f(2)=f(direction+1);f(direction+1)=swap
      result%flux(:,index(1),index(2),index(3))=f
      result%pressure(index(1),index(2),index(3))=pface;result%velocity(index(1),index(2),index(3))=uface
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
      dt,flattening,prediction,ok,order,boundary)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: ghost(:,-2:,-2:,-2:),temperature(-2:,-2:,-2:),volume_fraction(-2:,-2:,-2:)
    real(dp), intent(in) :: dt,flattening(0:,0:,0:)
    type(eb_plm_face_prediction_3d), intent(inout) :: prediction(3)
    type(eb_plm_boundary_3d), optional, intent(in) :: boundary
    logical, intent(out) :: ok
    integer, optional, intent(in) :: order
    type(eb_plm_face_prediction_3d) :: candidate(3)
    real(dp), allocatable :: primitive(:,:,:,:),sound(:,:,:),energy(:,:,:)
    real(dp) :: stencil(reactive_nvar(size(species)),-2:2),epencil(-2:2)
    real(dp) :: ql(reactive_nvar(size(species))),qr(reactive_nvar(size(species))),unused(reactive_nvar(size(species)))
    real(dp) :: f(reactive_nvar(size(species))),el,er,ignored,aperture,spacing(3),recovered,speed,swap,pface,uface
    integer :: n(3),nv,axis,i,j,k,idx(3),lo(3),hi(3),cell(3),point(3),offset,side,selected
    logical :: covered(-2:2),valid
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
        call reactive_pelec_primitive_flux_x(species,ql,qr,.5_dp*(sound(i,j,k)+sound(cell(1),cell(2),cell(3))), &
          f,valid,interface_pressure=pface,interface_velocity=uface)
        if(.not.valid) return
        swap=f(2);f(2)=f(axis+1);f(axis+1)=swap
        result%flux(:,i,j,k)=f;result%pressure(i,j,k)=pface;result%velocity(i,j,k)=uface
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
    call reconstruct_pelec_plm_cell(rotated,sound,dtdx,flat,qm,qp,valid,selected,covered,opened)
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
      new_temperature, ok, target_volume_fraction, boundary_change, wall_change, difmag)
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
      new_temperature, local_ok, target)
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
