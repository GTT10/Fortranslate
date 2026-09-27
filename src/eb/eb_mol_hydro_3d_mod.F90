module eb_mol_hydro_3d_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use mol_riemann_3d_mod, only: mol_riemann_flux_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_reactive_wall_flux_3d_mod, only: interpolate_reactive_eb_centroid_fluxes_3d, &
    reactive_eb_flux_divergence_3d,reactive_eb_slip_wall_source_3d
  implicit none
  private
  public :: compute_mol_eb_hydro_3d
contains
  ! Prepared geometry/state contain three same-resolution halo cells on each
  ! side. Return central spatial rates and dt-integrated intensive face fluxes.
  ! Caller supplies physical boundary ghosts and performs redistribution.
  subroutine compute_mol_eb_hydro_3d(species,geometry,state,temperature,dt,order, &
      hydro_rate,integrated_flux,boundary_change,wall_change,ok,lax_friedrichs,frame)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),dt
    integer, intent(in) :: order
    real(dp), intent(inout) :: hydro_rate(:,:,:,:),boundary_change(:),wall_change(:)
    type(eb_integrated_flux_3d), intent(inout) :: integrated_flux
    logical, intent(out) :: ok
    logical, optional, intent(in) :: lax_friedrichs
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), allocatable :: potential(:,:,:),wall_potential(:,:,:),potential_cell,potential_face(:)
    real(dp) :: position(3)
    type(eb_geometry_3d) :: margin,physical
    type(eb_integrated_flux_3d) :: raw,flux
    real(dp), allocatable :: primitive(:,:,:,:),sound(:,:,:),rates(:,:,:,:),walls(:,:,:,:)
    real(dp) :: stencil(reactive_nprim(size(species)),-1:2),cs(2),f(reactive_nvar(size(species)))
    real(dp) :: boundary(reactive_nvar(size(species))),wall(reactive_nvar(size(species))),t,volume
    integer :: n(3),m(3),out(3),nv,np,axis,i,j,k,c,cell(3),peer(3),idx(3),lo(3),hi(3),s
    logical :: valid,connected(2,2)
    ok=.false.;n=[geometry%nx,geometry%ny,geometry%nz];out=n-6;m=n-4
    nv=reactive_nvar(size(species));np=reactive_nprim(size(species))
    if(any(out<1).or..not.geometry%is_valid()) return
    if(order/=1.and.order/=2) return
    if(.not.ieee_is_finite(dt)) return
    if(dt<=0) return
    if(any(shape(state)/=[nv,n]).or.any(shape(temperature)/=n)) return
    if(any(shape(hydro_rate)/=[nv,out]).or.size(boundary_change)/=nv.or.size(wall_change)/=nv) return
    allocate(primitive(np,n(1),n(2),n(3)),sound(n(1),n(2),n(3)))
    primitive=0;sound=0
    if(present(frame)) then
      if(.not.frame%is_valid()) return
      allocate(potential(n(1),n(2),n(3)),wall_potential(out(1),out(2),out(3)),potential_cell,potential_face(2))
      potential=0
    end if
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      if(geometry%volume_fraction(i,j,k)<=0) cycle
      if(present(frame)) then
        position=[geometry%x_lower+(real(i,dp)-.5_dp)*geometry%dx, &
          geometry%y_lower+(real(j,dp)-.5_dp)*geometry%dy,geometry%z_lower+(real(k,dp)-.5_dp)*geometry%dz]
        call rotating_potential_3d(frame,position,potential_cell,valid)
        if(.not.valid) return
        potential(i,j,k)=potential_cell
      end if
      call reactive_conserved_to_primitive(species,state(:,i,j,k),temperature(i,j,k), &
        primitive(:,i,j,k),t,sound(i,j,k),valid,potential_cell)
      if(.not.valid) return
    end do
    end do
    end do
    if(allocated(wall_potential)) wall_potential=potential(4:n(1)-3,4:n(2)-3,4:n(3)-3)
    call extract_eb_geometry_3d(geometry,[3,3,3],n-2,margin,valid)
    if(.not.valid) return
    call extract_eb_geometry_3d(geometry,[4,4,4],n-3,physical,valid)
    if(.not.valid) return
    allocate(raw%x(nv,0:m(1),m(2),m(3)),raw%y(nv,m(1),0:m(2),m(3)),raw%z(nv,m(1),m(2),0:m(3)))
    raw%x=0;raw%y=0;raw%z=0
    do axis=1,3
      lo=1;lo(axis)=0;hi=m
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        idx=[i,j,k];cell=idx+2
        if(aperture(axis,cell)<=0) cycle
        peer=cell;peer(axis)=peer(axis)+1
        if(geometry%volume_fraction(cell(1),cell(2),cell(3))<=0.or. &
          geometry%volume_fraction(peer(1),peer(2),peer(3))<=0) return
        cs=[sound(cell(1),cell(2),cell(3)),sound(peer(1),peer(2),peer(3))]
        if(allocated(potential_face)) potential_face= &
          [potential(cell(1),cell(2),cell(3)),potential(peer(1),peer(2),peer(3))]
        connected=.true.
        peer=cell;peer(axis)=peer(axis)-1
        connected(1,1)=aperture(axis,peer)>0
        peer=cell;peer(axis)=peer(axis)+1
        connected(2,2)=aperture(axis,peer)>0
        do s=-1,2
          peer=cell;peer(axis)=peer(axis)+s
          stencil(:,s)=primitive(:,peer(1),peer(2),peer(3))
        end do
        call mol_riemann_flux_3d(species,stencil,cs,axis,connected,order,f,valid,lax_friedrichs,potential_face)
        if(.not.valid) return
        select case(axis)
        case(1);raw%x(:,i,j,k)=f
        case(2);raw%y(:,i,j,k)=f
        case(3);raw%z(:,i,j,k)=f
        end select
      end do
      end do
      end do
    end do
    allocate(flux%x(nv,0:out(1),out(2),out(3)),flux%y(nv,out(1),0:out(2),out(3)), &
      flux%z(nv,out(1),out(2),0:out(3)),rates(nv,out(1),out(2),out(3)),walls(nv,out(1),out(2),out(3)))
    call interpolate_reactive_eb_centroid_fluxes_3d(margin,raw%x,raw%y,raw%z,flux%x,flux%y,flux%z,valid,mol_stencil=.true.)
    if(.not.valid) return
    call reactive_eb_flux_divergence_3d(species,state(:,4:n(1)-3,4:n(2)-3,4:n(3)-3), &
      temperature(4:n(1)-3,4:n(2)-3,4:n(3)-3),physical,flux%x,flux%y,flux%z,rates,valid,wall_potential)
    if(.not.valid) return
    call reactive_eb_slip_wall_source_3d(species,state(:,4:n(1)-3,4:n(2)-3,4:n(3)-3), &
      temperature(4:n(1)-3,4:n(2)-3,4:n(3)-3),physical,walls,valid,wall_potential)
    if(.not.valid) return
    volume=physical%dx*physical%dy*physical%dz
    do c=1,nv
      boundary(c)=dt*(physical%dy*physical%dz*sum( &
        physical%x_face_fraction(0,:,:)*flux%x(c,0,:,:)-physical%x_face_fraction(out(1),:,:)*flux%x(c,out(1),:,:))+ &
        physical%dx*physical%dz*sum(physical%y_face_fraction(:,0,:)*flux%y(c,:,0,:)- &
        physical%y_face_fraction(:,out(2),:)*flux%y(c,:,out(2),:))+ &
        physical%dx*physical%dy*sum(physical%z_face_fraction(:,:,0)*flux%z(c,:,:,0)- &
        physical%z_face_fraction(:,:,out(3))*flux%z(c,:,:,out(3))))
      wall(c)=dt*volume*sum(physical%volume_fraction*walls(c,:,:,:))
    end do
    flux%x=dt*flux%x;flux%y=dt*flux%y;flux%z=dt*flux%z
    if(.not.all(ieee_is_finite(rates)).or..not.all(ieee_is_finite(boundary)).or..not.all(ieee_is_finite(wall))) return
    if(.not.all(ieee_is_finite(flux%x)).or..not.all(ieee_is_finite(flux%y)).or. &
      .not.all(ieee_is_finite(flux%z))) return
    hydro_rate=rates;boundary_change=boundary;wall_change=wall
    call move_alloc(flux%x,integrated_flux%x);call move_alloc(flux%y,integrated_flux%y)
    call move_alloc(flux%z,integrated_flux%z);ok=.true.
  contains
    real(dp) function aperture(direction,index) result(a)
      integer, intent(in) :: direction,index(3)
      select case(direction)
      case(1);a=geometry%x_face_fraction(index(1),index(2),index(3))
      case(2);a=geometry%y_face_fraction(index(1),index(2),index(3))
      case(3);a=geometry%z_face_fraction(index(1),index(2),index(3))
      end select
    end function
  end subroutine
end module eb_mol_hydro_3d_mod
