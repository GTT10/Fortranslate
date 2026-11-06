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
  public :: compute_mol_eb_hydro_3d,limit_mol_centroid_species_fluxes_3d
  ! Diagnostic snapshots are instantaneous (not dt-integrated). pressure and
  ! raw_flux include one transverse margin: face index zero has input index 2.
  ! centroid_flux/rate/wall_rate use the physical box (input cells 4:n-3),
  ! including spatial-rate halos when requested. Closed-face pressures are zero;
  ! consult geometry apertures before interpreting them. Caller owns stage IDs.
  type, public :: eb_mol_diagnostics_3d
    type(eb_integrated_flux_3d) :: pressure,kinematics,raw_flux,centroid_flux
    real(dp), allocatable :: cell_pressure(:,:,:),rate(:,:,:,:),wall_rate(:,:,:,:)
    real(dp), allocatable :: cell_primitive(:,:,:,:),cell_sound(:,:,:)
  end type
contains
  ! Prepared geometry/state contain three same-resolution halo cells on each
  ! side. Return central spatial rates and dt-integrated intensive face fluxes.
  ! Caller supplies physical boundary ghosts and performs redistribution.
  subroutine compute_mol_eb_hydro_3d(species,geometry,state,temperature,dt,order, &
      hydro_rate,integrated_flux,boundary_change,wall_change,ok,lax_friedrichs,frame, &
      raw_face_flux,limiter_alpha,include_boundary_correction,spatial_rate_halo, &
      internal_energy_rate,internal_energy_integrated_flux,internal_energy_raw_face_flux, &
      internal_energy_spatial_rate_halo,internal_energy_density,diagnostics,primitive_from_total_energy, &
      preserve_species_ratios)
    ! Select one energy basis for both face primitives and EB wall pressure.
    ! Held-energy transport remains independent; default preserves held recovery.
    logical, optional, intent(in) :: primitive_from_total_energy
    logical, optional, intent(in) :: preserve_species_ratios
    logical :: use_total_energy
    type(eb_mol_diagnostics_3d), optional, intent(inout) :: diagnostics
    type(eb_mol_diagnostics_3d) :: diagnostic_candidate
    real(dp) :: face_pressure(2),face_kinematics(6)
    real(dp), optional, intent(in) :: internal_energy_density(:,:,:)
    real(dp), optional, intent(inout) :: internal_energy_rate(:,:,:)
    type(eb_integrated_flux_3d), optional, intent(inout) :: internal_energy_integrated_flux
    type(eb_integrated_flux_3d), optional, intent(inout) :: internal_energy_raw_face_flux
    real(dp), allocatable, optional, intent(inout) :: internal_energy_spatial_rate_halo(:,:,:)
    real(dp), allocatable, optional, intent(inout) :: spatial_rate_halo(:,:,:,:)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),dt
    integer, intent(in) :: order
    real(dp), intent(inout) :: hydro_rate(:,:,:,:),boundary_change(:),wall_change(:)
    type(eb_integrated_flux_3d), intent(inout) :: integrated_flux
    logical, intent(out) :: ok
    logical, optional, intent(in) :: lax_friedrichs
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(eb_integrated_flux_3d), optional, intent(inout) :: raw_face_flux
    real(dp), allocatable, optional, intent(inout) :: limiter_alpha(:,:,:,:)
    logical, optional, intent(in) :: include_boundary_correction
    real(dp), allocatable :: potential(:,:,:),wall_potential(:,:,:),potential_cell,potential_face(:),held_cell,wall_held(:,:,:)
    real(dp) :: position(3)
    type(eb_geometry_3d) :: margin,physical,valid_geometry
    type(eb_integrated_flux_3d) :: raw,flux
    type(eb_integrated_flux_3d) :: raw_candidate,valid_flux,valid_raw,auxiliary_flux,auxiliary_raw,valid_auxiliary_raw
    real(dp), allocatable :: auxiliary_rate(:,:,:),auxiliary_halo_candidate(:,:,:)
    real(dp) :: auxiliary_face,area(3)
    real(dp), allocatable :: alpha_candidate(:,:,:,:),halo_candidate(:,:,:,:)
    real(dp), allocatable :: primitive(:,:,:,:),sound(:,:,:),rates(:,:,:,:),walls(:,:,:,:)
    real(dp) :: stencil(reactive_nprim(size(species)),-1:2),cs(2),f(reactive_nvar(size(species)))
    real(dp) :: boundary(reactive_nvar(size(species))),wall(reactive_nvar(size(species))),t,volume
    integer :: n(3),m(3),out(3),valid_n(3),g,nv,nf,np,axis,i,j,k,c,cell(3),peer(3),idx(3),lo(3),hi(3),s
    logical :: valid,connected(2,2),use_lax,include_boundary
    use_total_energy=.false.
    if(present(primitive_from_total_energy)) use_total_energy=primitive_from_total_energy
    ok=.false.;n=[geometry%nx,geometry%ny,geometry%nz];out=n-6;m=n-4
    nv=reactive_nvar(size(species));np=reactive_nprim(size(species));nf=nv
    if(present(internal_energy_rate).neqv.present(internal_energy_integrated_flux)) return
    if(present(internal_energy_raw_face_flux).and. &
      (.not.present(internal_energy_rate).or..not.present(raw_face_flux))) return
    if(present(internal_energy_spatial_rate_halo).and..not.present(internal_energy_rate)) return
    if(present(internal_energy_rate)) nf=nv+1
    ! With a halo output, six input ghost layers provide the three layers
    ! of spatial RHS used by StateRedist. The same operator computes every
    ! valid and ghost rate; face fluxes/diagnostics remain on the valid box.
    g=0
    if(present(spatial_rate_halo).or.present(internal_energy_spatial_rate_halo)) g=3
    valid_n=out-2*g
    if(any(valid_n<1).or..not.geometry%is_valid()) return
    if(order/=1.and.order/=2) return
    use_lax=.false.
    if(present(lax_friedrichs)) use_lax=lax_friedrichs
    include_boundary=.false.
    if(present(include_boundary_correction)) include_boundary=include_boundary_correction
    if((present(raw_face_flux).or.present(limiter_alpha).or.include_boundary).and..not.use_lax) return
    if(present(raw_face_flux).neqv.present(limiter_alpha)) return
    if(include_boundary.and..not.present(limiter_alpha)) return
    if(.not.ieee_is_finite(dt)) return
    if(dt<=0) return
    if(any(shape(state)/=[nv,n]).or.any(shape(temperature)/=n)) return
    if(any(shape(hydro_rate)/=[nv,valid_n]).or.size(boundary_change)/=nv.or.size(wall_change)/=nv) return
    if(present(internal_energy_rate)) then
      if(any(shape(internal_energy_rate)/=valid_n)) return
    endif
    if(present(internal_energy_density)) then
      if(any(shape(internal_energy_density)/=n)) return
      if(.not.use_total_energy) then
        allocate(held_cell)
        wall_held=internal_energy_density(4:n(1)-3,4:n(2)-3,4:n(3)-3)
      endif
    endif
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
      if(allocated(held_cell)) held_cell=internal_energy_density(i,j,k)
      call reactive_conserved_to_primitive(species,state(:,i,j,k),temperature(i,j,k), &
        primitive(:,i,j,k),t,sound(i,j,k),valid,potential_cell,held_cell,preserve_species_ratios)
      if(.not.valid) return
    end do
    end do
    end do
    if(allocated(wall_potential)) wall_potential=potential(4:n(1)-3,4:n(2)-3,4:n(3)-3)
    call extract_eb_geometry_3d(geometry,[3,3,3],n-2,margin,valid)
    if(.not.valid) return
    call extract_eb_geometry_3d(geometry,[4,4,4],n-3,physical,valid)
    if(.not.valid) return
    call extract_eb_geometry_3d(geometry,[4+g,4+g,4+g],n-3-g,valid_geometry,valid)
    if(.not.valid) return
    allocate(raw%x(nf,0:m(1),m(2),m(3)),raw%y(nf,m(1),0:m(2),m(3)),raw%z(nf,m(1),m(2),0:m(3)))
    raw%x=0;raw%y=0;raw%z=0
    if(present(diagnostics)) then
      allocate(diagnostic_candidate%pressure%x(2,0:m(1),m(2),m(3)), &
        diagnostic_candidate%pressure%y(2,m(1),0:m(2),m(3)), &
        diagnostic_candidate%pressure%z(2,m(1),m(2),0:m(3)))
      diagnostic_candidate%pressure%x=0;diagnostic_candidate%pressure%y=0;diagnostic_candidate%pressure%z=0
      allocate(diagnostic_candidate%kinematics%x(6,0:m(1),m(2),m(3)), &
        diagnostic_candidate%kinematics%y(6,m(1),0:m(2),m(3)), &
        diagnostic_candidate%kinematics%z(6,m(1),m(2),0:m(3)))
      diagnostic_candidate%kinematics%x=0;diagnostic_candidate%kinematics%y=0;diagnostic_candidate%kinematics%z=0
      diagnostic_candidate%cell_pressure=primitive(5,:,:,:)
      diagnostic_candidate%cell_primitive=primitive
      diagnostic_candidate%cell_sound=sound
    endif
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
        call mol_riemann_flux_3d(species,stencil,cs,axis,connected,order,f,valid, &
          lax_friedrichs,potential_face,auxiliary_face,face_pressure,face_kinematics)
        if(.not.valid) return
        if(present(diagnostics)) then
          select case(axis)
          case(1);diagnostic_candidate%pressure%x(:,i,j,k)=face_pressure
            diagnostic_candidate%kinematics%x(:,i,j,k)=face_kinematics
          case(2);diagnostic_candidate%pressure%y(:,i,j,k)=face_pressure
            diagnostic_candidate%kinematics%y(:,i,j,k)=face_kinematics
          case(3);diagnostic_candidate%pressure%z(:,i,j,k)=face_pressure
            diagnostic_candidate%kinematics%z(:,i,j,k)=face_kinematics
          end select
        endif
        select case(axis)
        case(1);raw%x(1:nv,i,j,k)=f
          if(nf>nv) raw%x(nf,i,j,k)=auxiliary_face
        case(2);raw%y(1:nv,i,j,k)=f
          if(nf>nv) raw%y(nf,i,j,k)=auxiliary_face
        case(3);raw%z(1:nv,i,j,k)=f
          if(nf>nv) raw%z(nf,i,j,k)=auxiliary_face
        end select
      end do
      end do
      end do
    end do
    allocate(flux%x(nf,0:out(1),out(2),out(3)),flux%y(nf,out(1),0:out(2),out(3)), &
      flux%z(nf,out(1),out(2),0:out(3)),rates(nv,out(1),out(2),out(3)),walls(nv,out(1),out(2),out(3)))
    call interpolate_reactive_eb_centroid_fluxes_3d(margin,raw%x,raw%y,raw%z,flux%x,flux%y,flux%z,valid,mol_stencil=.true.)
    if(.not.valid) return
    if(use_lax) then
      if(present(limiter_alpha)) then
        call limit_mol_centroid_species_fluxes_3d(physical,state(:,4:n(1)-3,4:n(2)-3,4:n(3)-3), &
          raw,flux,dt,valid,include_boundary,alpha_candidate)
      else
        call limit_mol_centroid_species_fluxes_3d(physical,state(:,4:n(1)-3,4:n(2)-3,4:n(3)-3), &
          raw,flux,dt,valid)
      end if
      if(.not.valid) return
    end if
    call reactive_eb_flux_divergence_3d(species,state(:,4:n(1)-3,4:n(2)-3,4:n(3)-3), &
      temperature(4:n(1)-3,4:n(2)-3,4:n(3)-3),physical,flux%x(1:nv,:,:,:),flux%y(1:nv,:,:,:),flux%z(1:nv,:,:,:), &
      rates,valid,wall_potential,wall_held)
    if(.not.valid) return
    call reactive_eb_slip_wall_source_3d(species,state(:,4:n(1)-3,4:n(2)-3,4:n(3)-3), &
      temperature(4:n(1)-3,4:n(2)-3,4:n(3)-3),physical,walls,valid,wall_potential,wall_held)
    if(.not.valid) return
    if(nf>nv) then
      allocate(auxiliary_rate(out(1),out(2),out(3)),source=0._dp)
      volume=physical%dx*physical%dy*physical%dz
      area=[physical%dy*physical%dz,physical%dx*physical%dz,physical%dx*physical%dy]
      ! Pinned MOL transports UEINT without the CTU pressure-work source.
      ! A stationary slip wall has zero internal-energy advective flux.
      do k=1,out(3)
      do j=1,out(2)
      do i=1,out(1)
        if(physical%volume_fraction(i,j,k)<=0) cycle
        auxiliary_rate(i,j,k)=( &
          area(1)*(physical%x_face_fraction(i-1,j,k)*flux%x(nf,i-1,j,k)- &
                   physical%x_face_fraction(i,j,k)*flux%x(nf,i,j,k))+ &
          area(2)*(physical%y_face_fraction(i,j-1,k)*flux%y(nf,i,j-1,k)- &
                   physical%y_face_fraction(i,j,k)*flux%y(nf,i,j,k))+ &
          area(3)*(physical%z_face_fraction(i,j,k-1)*flux%z(nf,i,j,k-1)- &
                   physical%z_face_fraction(i,j,k)*flux%z(nf,i,j,k)))/(volume*physical%volume_fraction(i,j,k))
      enddo
      enddo
      enddo
      if(.not.all(ieee_is_finite(auxiliary_rate))) return
    endif
    volume=valid_geometry%dx*valid_geometry%dy*valid_geometry%dz
    do c=1,nv
      boundary(c)=dt*(valid_geometry%dy*valid_geometry%dz*sum( &
        valid_geometry%x_face_fraction(0,:,:)*flux%x(c,g,g+1:g+valid_n(2),g+1:g+valid_n(3))- &
        valid_geometry%x_face_fraction(valid_n(1),:,:)*flux%x(c,g+valid_n(1),g+1:g+valid_n(2),g+1:g+valid_n(3)))+ &
        valid_geometry%dx*valid_geometry%dz*sum( &
        valid_geometry%y_face_fraction(:,0,:)*flux%y(c,g+1:g+valid_n(1),g,g+1:g+valid_n(3))- &
        valid_geometry%y_face_fraction(:,valid_n(2),:)*flux%y(c,g+1:g+valid_n(1),g+valid_n(2),g+1:g+valid_n(3)))+ &
        valid_geometry%dx*valid_geometry%dy*sum( &
        valid_geometry%z_face_fraction(:,:,0)*flux%z(c,g+1:g+valid_n(1),g+1:g+valid_n(2),g)- &
        valid_geometry%z_face_fraction(:,:,valid_n(3))*flux%z(c,g+1:g+valid_n(1),g+1:g+valid_n(2),g+valid_n(3))))
      wall(c)=dt*volume*sum(valid_geometry%volume_fraction*walls(c,g+1:g+valid_n(1),g+1:g+valid_n(2),g+1:g+valid_n(3)))
    end do
    if(present(diagnostics)) then
      diagnostic_candidate%raw_flux=raw
      diagnostic_candidate%centroid_flux=flux
      diagnostic_candidate%rate=rates
      diagnostic_candidate%wall_rate=walls
    endif
    flux%x=dt*flux%x;flux%y=dt*flux%y;flux%z=dt*flux%z
    if(.not.all(ieee_is_finite(rates)).or..not.all(ieee_is_finite(boundary)).or..not.all(ieee_is_finite(wall))) return
    if(.not.all(ieee_is_finite(flux%x)).or..not.all(ieee_is_finite(flux%y)).or. &
      .not.all(ieee_is_finite(flux%z))) return
    if(present(raw_face_flux)) then
      allocate(raw_candidate%x(nv,0:out(1),out(2),out(3)), &
        raw_candidate%y(nv,out(1),0:out(2),out(3)),raw_candidate%z(nv,out(1),out(2),0:out(3)))
      raw_candidate%x=dt*raw%x(1:nv,1:out(1)+1,2:out(2)+1,2:out(3)+1)
      raw_candidate%y=dt*raw%y(1:nv,2:out(1)+1,1:out(2)+1,2:out(3)+1)
      raw_candidate%z=dt*raw%z(1:nv,2:out(1)+1,2:out(2)+1,1:out(3)+1)
      if(.not.all(ieee_is_finite(raw_candidate%x)).or..not.all(ieee_is_finite(raw_candidate%y)).or. &
        .not.all(ieee_is_finite(raw_candidate%z))) return
      if(present(internal_energy_raw_face_flux)) then
        allocate(auxiliary_raw%x(1,0:out(1),out(2),out(3)), &
          auxiliary_raw%y(1,out(1),0:out(2),out(3)),auxiliary_raw%z(1,out(1),out(2),0:out(3)))
        auxiliary_raw%x(1,:,:,:)=dt*raw%x(nf,1:out(1)+1,2:out(2)+1,2:out(3)+1)
        auxiliary_raw%y(1,:,:,:)=dt*raw%y(nf,2:out(1)+1,1:out(2)+1,2:out(3)+1)
        auxiliary_raw%z(1,:,:,:)=dt*raw%z(nf,2:out(1)+1,2:out(2)+1,1:out(3)+1)
        if(.not.all(ieee_is_finite(auxiliary_raw%x)).or..not.all(ieee_is_finite(auxiliary_raw%y)).or. &
          .not.all(ieee_is_finite(auxiliary_raw%z))) return
      endif
    end if
    call crop_valid_flux(flux,valid_flux,nv)
    if(nf>nv) call crop_valid_flux(flux,auxiliary_flux,1,nf)
    if(present(raw_face_flux)) call crop_valid_flux(raw_candidate,valid_raw)
    if(present(internal_energy_raw_face_flux)) call crop_valid_flux(auxiliary_raw,valid_auxiliary_raw)
    if(present(spatial_rate_halo)) then
      allocate(halo_candidate(nv,n(1),n(2),n(3)),source=0._dp)
      halo_candidate(:,4:n(1)-3,4:n(2)-3,4:n(3)-3)=rates
    endif
    if(present(internal_energy_spatial_rate_halo)) then
      allocate(auxiliary_halo_candidate(n(1),n(2),n(3)),source=0._dp)
      auxiliary_halo_candidate(4:n(1)-3,4:n(2)-3,4:n(3)-3)=auxiliary_rate
    endif
    if(nf>nv) then
      internal_energy_rate=auxiliary_rate(g+1:g+valid_n(1),g+1:g+valid_n(2),g+1:g+valid_n(3))
      internal_energy_integrated_flux=auxiliary_flux
    endif
    hydro_rate=rates(:,g+1:g+valid_n(1),g+1:g+valid_n(2),g+1:g+valid_n(3))
    boundary_change=boundary;wall_change=wall
    call move_alloc(valid_flux%x,integrated_flux%x);call move_alloc(valid_flux%y,integrated_flux%y)
    call move_alloc(valid_flux%z,integrated_flux%z)
    if(present(raw_face_flux)) then
      raw_face_flux=valid_raw
      limiter_alpha=alpha_candidate(:,g+1:g+valid_n(1),g+1:g+valid_n(2),g+1:g+valid_n(3))
      if(present(internal_energy_raw_face_flux)) internal_energy_raw_face_flux=valid_auxiliary_raw
    endif
    if(present(spatial_rate_halo)) call move_alloc(halo_candidate,spatial_rate_halo)
    if(present(internal_energy_spatial_rate_halo)) &
      call move_alloc(auxiliary_halo_candidate,internal_energy_spatial_rate_halo)
    if(present(diagnostics)) diagnostics=diagnostic_candidate
    ok=.true.
  contains
    subroutine crop_valid_flux(source,destination,components,first_component)
      type(eb_integrated_flux_3d), intent(in) :: source
      type(eb_integrated_flux_3d), intent(out) :: destination
      integer, optional, intent(in) :: components,first_component
      integer :: count,start
      count=size(source%x,1);start=1
      if(present(components)) count=components
      if(present(first_component)) start=first_component
      allocate(destination%x(count,0:valid_n(1),valid_n(2),valid_n(3)), &
        destination%y(count,valid_n(1),0:valid_n(2),valid_n(3)), &
        destination%z(count,valid_n(1),valid_n(2),0:valid_n(3)))
      destination%x=source%x(start:start+count-1,g:g+valid_n(1),g+1:g+valid_n(2),g+1:g+valid_n(3))
      destination%y=source%y(start:start+count-1,g+1:g+valid_n(1),g:g+valid_n(2),g+1:g+valid_n(3))
      destination%z=source%z(start:start+count-1,g+1:g+valid_n(1),g+1:g+valid_n(2),g:g+valid_n(3))
    end subroutine crop_valid_flux
    real(dp) function aperture(direction,index) result(a)
      integer, intent(in) :: direction,index(3)
      select case(direction)
      case(1);a=geometry%x_face_fraction(index(1),index(2),index(3))
      case(2);a=geometry%y_face_fraction(index(1),index(2),index(3))
      case(3);a=geometry%z_face_fraction(index(1),index(2),index(3))
      end select
    end function
  end subroutine

  ! Bound the centroid correction on interior faces by the species mass left
  ! after the uncorrected LF flux. One coefficient scales every conserved
  ! component of a face correction, preserving its shared conservative flux.
  ! A caller may include boundary corrections in the cell allowance and
  ! receive the allowances for a later cross-patch face exchange.
  subroutine limit_mol_centroid_species_fluxes_3d(geometry,state,raw,flux,dt,ok,include_boundary,alpha_output)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: state(:,:,:,:),dt
    type(eb_integrated_flux_3d), intent(in) :: raw
    type(eb_integrated_flux_3d), intent(inout) :: flux
    logical, intent(out) :: ok
    logical, optional, intent(in) :: include_boundary
    real(dp), allocatable, optional, intent(inout) :: alpha_output(:,:,:,:)
    type(eb_integrated_flux_3d) :: candidate
    real(dp), allocatable :: alpha(:,:,:,:)
    real(dp) :: volume,area(3),available,outgoing,low_flux,high_flux,low_area,high_area
    real(dp) :: delta,theta
    integer :: n(3),nv,nf,ns,s,i,j,k,axis,face(3),left(3),right(3)
    logical :: boundary_correction
    ok=.false.;n=[geometry%nx,geometry%ny,geometry%nz];nv=size(state,1);ns=nv-5
    boundary_correction=.false.
    if(present(include_boundary)) boundary_correction=include_boundary
    if(ns<1.or.any(n<1).or.any(shape(state)/=[nv,n]).or..not.geometry%is_valid()) return
    if(.not.ieee_is_finite(dt).or.dt<=0) return
    if(.not.allocated(raw%x).or..not.allocated(raw%y).or..not.allocated(raw%z)) return
    if(.not.allocated(flux%x).or..not.allocated(flux%y).or..not.allocated(flux%z)) return
    nf=size(raw%x,1)
    if(nf<nv) return
    if(any(shape(raw%x)/=[nf,n(1)+3,n(2)+2,n(3)+2]).or. &
      any(shape(raw%y)/=[nf,n(1)+2,n(2)+3,n(3)+2]).or. &
      any(shape(raw%z)/=[nf,n(1)+2,n(2)+2,n(3)+3])) return
    if(any(shape(flux%x)/=[nf,n(1)+1,n(2),n(3)]).or. &
      any(shape(flux%y)/=[nf,n(1),n(2)+1,n(3)]).or. &
      any(shape(flux%z)/=[nf,n(1),n(2),n(3)+1])) return
    if(any(.not.ieee_is_finite(raw%x)).or.any(.not.ieee_is_finite(raw%y)).or. &
      any(.not.ieee_is_finite(raw%z))) return
    if(any(.not.ieee_is_finite(flux%x)).or.any(.not.ieee_is_finite(flux%y)).or. &
      any(.not.ieee_is_finite(flux%z))) return
    candidate=flux
    volume=geometry%dx*geometry%dy*geometry%dz
    area=[geometry%dy*geometry%dz,geometry%dx*geometry%dz,geometry%dx*geometry%dy]
    if(.not.ieee_is_finite(volume).or.volume<=0.or.any(.not.ieee_is_finite(area))) return
    allocate(alpha(ns,n(1),n(2),n(3)));alpha=1
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      if(geometry%volume_fraction(i,j,k)<=0) cycle
      do s=1,ns
        available=geometry%volume_fraction(i,j,k)*state(s+5,i,j,k)*volume
        outgoing=0
        do axis=1,3
          face=[i,j,k];face(axis)=face(axis)-1
          low_area=area(axis)*aperture(axis,face)
          low_flux=base_flux(axis,s+5,face)
          if(face(axis)>0.or.boundary_correction) &
            outgoing=outgoing+dt*low_area*max(0._dp,-correction(axis,s+5,face))
          face(axis)=face(axis)+1
          high_area=area(axis)*aperture(axis,face)
          high_flux=base_flux(axis,s+5,face)
          if(face(axis)<n(axis).or.boundary_correction) &
            outgoing=outgoing+dt*high_area*max(0._dp,correction(axis,s+5,face))
          available=available+dt*(low_area*low_flux-high_area*high_flux)
        end do
        if(.not.ieee_is_finite(available).or..not.ieee_is_finite(outgoing)) return
        if(available<0) cycle
        if(outgoing>available) alpha(s,i,j,k)=max(0._dp,(available/outgoing)*(1-16*epsilon(1._dp)))
      end do
    end do
    end do
    end do
    do axis=1,3
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        face=[i,j,k]
        if(face(axis)>=n(axis)) cycle
        left=face;right=face;right(axis)=right(axis)+1
        theta=1
        do s=1,ns
          delta=correction(axis,s+5,face)
          if(delta>0.and.geometry%volume_fraction(left(1),left(2),left(3))>0) &
            theta=min(theta,alpha(s,left(1),left(2),left(3)))
          if(delta<0.and.geometry%volume_fraction(right(1),right(2),right(3))>0) &
            theta=min(theta,alpha(s,right(1),right(2),right(3)))
        end do
        if(theta>=1) cycle
        select case(axis)
        case(1)
          candidate%x(:,i,j,k)=raw%x(:,i+1,j+1,k+1)+theta*(flux%x(:,i,j,k)-raw%x(:,i+1,j+1,k+1))
        case(2)
          candidate%y(:,i,j,k)=raw%y(:,i+1,j+1,k+1)+theta*(flux%y(:,i,j,k)-raw%y(:,i+1,j+1,k+1))
        case(3)
          candidate%z(:,i,j,k)=raw%z(:,i+1,j+1,k+1)+theta*(flux%z(:,i,j,k)-raw%z(:,i+1,j+1,k+1))
        end select
      end do
      end do
      end do
    end do
    if(any(.not.ieee_is_finite(candidate%x)).or.any(.not.ieee_is_finite(candidate%y)).or. &
      any(.not.ieee_is_finite(candidate%z))) return
    call move_alloc(candidate%x,flux%x);call move_alloc(candidate%y,flux%y)
    call move_alloc(candidate%z,flux%z)
    if(present(alpha_output)) alpha_output=alpha
    ok=.true.
  contains
    real(dp) function aperture(direction,index) result(value)
      integer, intent(in) :: direction,index(3)
      select case(direction)
      case(1);value=geometry%x_face_fraction(index(1),index(2),index(3))
      case(2);value=geometry%y_face_fraction(index(1),index(2),index(3))
      case(3);value=geometry%z_face_fraction(index(1),index(2),index(3))
      end select
    end function
    real(dp) function base_flux(direction,component,index) result(value)
      integer, intent(in) :: direction,component,index(3)
      if(.not.boundary_correction.and.(index(direction)==0.or.index(direction)==n(direction))) then
        select case(direction)
        case(1);value=flux%x(component,index(1),index(2),index(3))
        case(2);value=flux%y(component,index(1),index(2),index(3))
        case(3);value=flux%z(component,index(1),index(2),index(3))
        end select
      else
        select case(direction)
        case(1);value=raw%x(component,index(1)+1,index(2)+1,index(3)+1)
        case(2);value=raw%y(component,index(1)+1,index(2)+1,index(3)+1)
        case(3);value=raw%z(component,index(1)+1,index(2)+1,index(3)+1)
        end select
      end if
    end function
    real(dp) function correction(direction,component,index) result(value)
      integer, intent(in) :: direction,component,index(3)
      value=0
      if(.not.boundary_correction.and.(index(direction)==0.or.index(direction)==n(direction))) return
      select case(direction)
      case(1);value=flux%x(component,index(1),index(2),index(3))- &
        raw%x(component,index(1)+1,index(2)+1,index(3)+1)
      case(2);value=flux%y(component,index(1),index(2),index(3))- &
        raw%y(component,index(1)+1,index(2)+1,index(3)+1)
      case(3);value=flux%z(component,index(1),index(2),index(3))- &
        raw%z(component,index(1)+1,index(2)+1,index(3)+1)
      end select
    end function
  end subroutine limit_mol_centroid_species_fluxes_3d
end module eb_mol_hydro_3d_mod
