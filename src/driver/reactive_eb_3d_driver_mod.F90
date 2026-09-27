module reactive_eb_3d_driver_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use constant_volume_reactor_mod, only: forced_reactor_integrator
  use state_indices_mod, only: imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_mass_fraction_component, reactive_conserved_to_primitive, &
    reactive_primitive_to_conserved
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d, eb_cut_cell_3d,extract_eb_geometry_3d,periodic_eb_faces_match
  use eb_reactive_hydro_3d_mod, only: &
    advance_reactive_eb_redistributed_euler_3d,advance_reactive_eb_boundary_ctu_3d, &
    advance_reactive_eb_state_redistributed_euler_3d
  use reactive_3d_mod, only: advance_reactive_chemistry_3d,react_sdc_state_3d
  use eb_reactive_transport_3d_mod, only: advance_reactive_eb_transport_3d,eb_transport_halo_provider_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d,valid_eb_integrated_flux_3d
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use eb_reactive_hydro_3d_mod, only: compute_sdc_eb_boundary_hydro_3d
  use eb_reactive_transport_3d_mod, only: compute_sdc_eb_transport_3d
  use eb_reactive_redistribution_3d_mod, only: redistribute_sdc_eb_rate_3d
  use simulation_config_reactive_3d_mod, only: fill_reactive_boundary_ghosts_3d
  implicit none
  private

  public :: initialize_reactive_eb_density_sheet_3d
  public :: advance_reactive_eb_strang_3d
  public :: advance_reactive_eb_full_3d
  public :: advance_reactive_eb_sdc_3d
  public :: react_sdc_eb_state_3d
  public :: reactive_eb_integrals_3d
  public :: reactive_eb_element_amounts_3d
  public :: reactive_eb_composition_integrals_3d
  public :: reactive_eb_element_integrals_3d
  public :: reactive_eb_element_species_supported_3d
  public :: reactive_eb_extrema_3d
  public :: write_reactive_eb_3d_csv

contains
  subroutine advance_reactive_eb_sdc_3d(species,reactions,transport,state,temperature,geometry,extended,dt, &
      iterations,chemistry,transport_flags,force,rtol,atol,redistribution,boundary_faces, &
      boundary_change,wall_change,source_change,minimum_theta,ok,reaction_source,integrated_flux, &
      difmag,target_volume_fraction,inflow_state,inflow_temperature,wall_temperature,chemistry_change, &
      evaluate_reaction)
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(eb_geometry_3d), intent(in) :: geometry,extended
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    real(dp), intent(in) :: dt,force(3),rtol,atol
    integer, intent(in) :: iterations
    logical, intent(in) :: chemistry,transport_flags(4)
    character(len=*), intent(in) :: redistribution,boundary_faces(6)
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),source_change(:),minimum_theta
    logical, intent(out) :: ok
    real(dp), optional, intent(inout) :: chemistry_change(:)
    real(dp), optional, intent(inout) :: reaction_source(:,:,:,:)
    type(eb_integrated_flux_3d), optional, intent(inout) :: integrated_flux
    real(dp), optional, intent(in) :: difmag,target_volume_fraction,inflow_state(:,:),inflow_temperature(:),wall_temperature(6)
    real(dp), allocatable :: old(:,:,:,:),ot(:,:,:),u(:,:,:,:),t(:,:,:),ir(:,:,:,:),hydro(:,:,:,:),rates(:,:,:,:)
    real(dp), allocatable :: e0(:,:,:,:),e1(:,:,:,:),d0(:,:,:,:),d1(:,:,:,:),raw(:,:,:,:)
    type(eb_integrated_flux_3d) :: hf,f0,f1
    real(dp) :: boundary(size(boundary_change)),wall(size(wall_change)),external(size(source_change))
    real(dp) :: v0(size(source_change)),v1(size(source_change)),theta,stheta
    integer :: nv,n(3),i,j,k,it
    logical :: valid,diffusive
    ok=.false.;nv=reactive_nvar(size(species));n=[geometry%nx,geometry%ny,geometry%nz]
    if(present(evaluate_reaction).and..not.chemistry) return
    if(.not.geometry%is_valid().or..not.valid_ctu_geometry(geometry,extended,'pelec')) return
    if(any(shape(state)/=[nv,n]).or.any(shape(temperature)/=n)) return
    if(size(boundary_change)/=nv.or.size(wall_change)/=nv.or.size(source_change)/=nv) return
    if(iterations<1.or..not.all(ieee_is_finite([dt,force,rtol,atol]))) return
    if(min(dt,rtol,atol)<=0.or..not.valid_reactive_boundaries_3d(boundary_faces)) return
    if(transport_flags(4).and..not.transport_flags(3)) return
    if(present(chemistry_change)) then
      if(size(chemistry_change)/=nv) return
    end if
    if(present(reaction_source)) then
      if(.not.chemistry.or.any(shape(reaction_source)/=[nv-4,n])) return
    end if
    old=state;ot=temperature
    allocate(ir(nv-4,n(1),n(2),n(3)),hydro(nv,n(1),n(2),n(3)),raw(nv,n(1),n(2),n(3)));ir=0
    if(present(reaction_source)) ir=reaction_source
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      if(geometry%cell_type(i,j,k)/=eb_covered_cell_3d) cycle
      old(:,i,j,k)=0;ot(i,j,k)=0;ir(:,i,j,k)=0
    end do
    end do
    end do
    if(.not.all(ieee_is_finite(ir))) return
    call recover(old,ot,valid)
    if(.not.valid) return
    diffusive=any(transport_flags(:3));theta=1
    call sources(old,ot,e0,d0,f0,valid)
    if(.not.valid) return
    e1=e0;d1=d0;u=old;t=ot
    do it=1,iterations
      rates=.5_dp*e1+.5_dp*e0+.5_dp*d1+.5_dp*d0
      rates(6:,:,:,:)=rates(6:,:,:,:)+ir(:nv-5,:,:,:)
      rates(5,:,:,:)=rates(5,:,:,:)+ir(nv-4,:,:,:)
      call compute_sdc_eb_boundary_hydro_3d(species,old,ot,extended,dt,rates,raw,hf,boundary,wall,valid, &
        boundary_faces,difmag=difmag,inflow_state=inflow_state,inflow_temperature=inflow_temperature)
      if(.not.valid) return
      call redistribute_sdc_eb_rate_3d(geometry,old,raw,dt,redistribution,hydro,valid,target_volume_fraction)
      if(.not.valid) return
      u=old+dt*(rates+hydro)
      call recover(u,t,valid)
      if(.not.valid) return
      call sources(u,t,e1,d1,f1,valid)
      if(.not.valid) return
      rates=.5_dp*e1+.5_dp*e0+.5_dp*d1+.5_dp*d0+hydro
      if(chemistry) then
        call react_sdc_eb_state_3d(species,reactions,old,ot,geometry,rates,dt,rtol,atol,u,t,ir,valid, &
          evaluate_reaction=evaluate_reaction)
      else
        u=old+dt*rates
        call recover(u,t,valid)
      end if
      if(.not.valid) return
    end do
    if(diffusive) then
      hf%x=hf%x+.5_dp*(f0%x+f1%x);hf%y=hf%y+.5_dp*(f0%y+f1%y);hf%z=hf%z+.5_dp*(f0%z+f1%z)
      call transport_inventory(f0,v0)
      call transport_inventory(f1,v1)
      boundary=boundary+.5_dp*(v0+v1)
    end if
    call reactive_eb_integrals_3d(e0,geometry,v0,valid)
    if(.not.valid) return
    call reactive_eb_integrals_3d(e1,geometry,v1,valid)
    if(.not.valid) return
    external=.5_dp*dt*(v0+v1)
    if(.not.all(ieee_is_finite(boundary)).or..not.all(ieee_is_finite(wall)).or. &
      .not.all(ieee_is_finite(external))) return
    if(.not.all(ieee_is_finite(hf%x)).or..not.all(ieee_is_finite(hf%y)).or..not.all(ieee_is_finite(hf%z))) return
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      if(geometry%cell_type(i,j,k)/=eb_covered_cell_3d) cycle
      u(:,i,j,k)=state(:,i,j,k);t(i,j,k)=temperature(i,j,k)
      if(present(reaction_source)) ir(:,i,j,k)=reaction_source(:,i,j,k)
    end do
    end do
    end do
    if(present(chemistry_change)) then
      raw=0
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(geometry%cell_type(i,j,k)==eb_covered_cell_3d) cycle
        if(chemistry) raw(:,i,j,k)=u(:,i,j,k)-(old(:,i,j,k)+dt*rates(:,i,j,k))
      end do
      end do
      end do
      call reactive_eb_integrals_3d(raw,geometry,v0,valid)
      if(.not.valid) return
    end if
    state=u;temperature=t;boundary_change=boundary;wall_change=wall;source_change=external;minimum_theta=theta
    if(present(reaction_source)) reaction_source=ir
    if(present(chemistry_change)) chemistry_change=v0
    if(present(integrated_flux)) integrated_flux=hf
    ok=.true.
  contains
    subroutine transport_inventory(flux,result)
      type(eb_integrated_flux_3d), intent(in) :: flux
      real(dp), intent(out) :: result(:)
      integer :: c
      do c=1,nv
        result(c)=geometry%dy*geometry%dz*sum(geometry%x_face_fraction(0,:,:)*flux%x(c,0,:,:)- &
          geometry%x_face_fraction(n(1),:,:)*flux%x(c,n(1),:,:))+ &
          geometry%dx*geometry%dz*sum(geometry%y_face_fraction(:,0,:)*flux%y(c,:,0,:)- &
          geometry%y_face_fraction(:,n(2),:)*flux%y(c,:,n(2),:))+ &
          geometry%dx*geometry%dy*sum(geometry%z_face_fraction(:,:,0)*flux%z(c,:,:,0)- &
          geometry%z_face_fraction(:,:,n(3))*flux%z(c,:,:,n(3)))
      end do
    end subroutine transport_inventory
    subroutine recover(q,temp,good)
      real(dp), intent(in) :: q(:,:,:,:)
      real(dp), intent(inout) :: temp(:,:,:)
      logical, intent(out) :: good
      real(dp) :: primitive(nv),recovered,speed
      integer :: ii,jj,kk
      good=.false.
      do kk=1,n(3)
      do jj=1,n(2)
      do ii=1,n(1)
        if(geometry%cell_type(ii,jj,kk)==eb_covered_cell_3d) cycle
        call reactive_conserved_to_primitive(species,q(:,ii,jj,kk),temp(ii,jj,kk),primitive,recovered,speed,good)
        if(.not.good) return
        temp(ii,jj,kk)=recovered
      end do
      end do
      end do
      good=.true.
    end subroutine recover
    subroutine sources(q,temp,e,d,flux,good)
      real(dp), intent(in) :: q(:,:,:,:),temp(:,:,:)
      real(dp), allocatable, intent(out) :: e(:,:,:,:),d(:,:,:,:)
      type(eb_integrated_flux_3d), intent(out) :: flux
      logical, intent(out) :: good
      real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),unredistributed(:,:,:,:)
      integer :: ii,jj,kk
      allocate(e,mold=q);allocate(d,mold=q);e=0;d=0;good=.false.
      do kk=1,n(3)
      do jj=1,n(2)
      do ii=1,n(1)
        if(geometry%cell_type(ii,jj,kk)==eb_covered_cell_3d) cycle
        if(q(1,ii,jj,kk)<=0) return
        e(2:4,ii,jj,kk)=force
        e(5,ii,jj,kk)=dot_product(q(2:4,ii,jj,kk)/q(1,ii,jj,kk),force)
      end do
      end do
      end do
      if(.not.all(ieee_is_finite(e))) return
      if(diffusive) then
        allocate(ghost(nv,0:n(1)+1,0:n(2)+1,0:n(3)+1),gt(0:n(1)+1,0:n(2)+1,0:n(3)+1))
        allocate(unredistributed,mold=q)
        call fill_reactive_boundary_ghosts_3d(q,temp,boundary_faces,ghost,gt,good,viscous_wall=.true., &
          inflow_state=inflow_state,inflow_temperature=inflow_temperature,ghost_width=1)
        if(.not.good) return
        call compute_sdc_eb_transport_3d(species,transport,q,temp,geometry,dt,transport_flags, &
          unredistributed,flux,stheta,good,ghost,gt,[.true.,.true.,.true.,.true.,.true.,.true.], &
          boundary_faces,wall_temperature)
        if(.not.good) return
        theta=min(theta,stheta)
        call redistribute_sdc_eb_rate_3d(geometry,q,unredistributed,dt,redistribution,d,good,target_volume_fraction)
        if(.not.good) return
      end if
      good=.true.
    end subroutine sources
  end subroutine advance_reactive_eb_sdc_3d


  ! Rates and states are per fluid volume; geometry supplies the active mask.
  subroutine react_sdc_eb_state_3d(species,reactions,old_state,old_temperature,geometry,nonreacting_rates, &
      dt,rtol,atol,state,temperature,reaction_source,ok, &
      evaluate_reaction)
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: old_state(:,:,:,:),old_temperature(:,:,:),nonreacting_rates(:,:,:,:),dt,rtol,atol
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:),reaction_source(:,:,:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: candidate(:,:,:,:),ct(:,:,:),ir(:,:,:,:)
    integer :: nv,n(3),i,j,k
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    if(.not.geometry%is_valid().or.size(species)<1.or.size(reactions)<1) return
    n=[geometry%nx,geometry%ny,geometry%nz]
    if(any(shape(old_state)/=[nv,n]).or.any(shape(nonreacting_rates)/=[nv,n])) return
    if(any(shape(state)/=[nv,n]).or.any(shape(old_temperature)/=n).or.any(shape(temperature)/=n)) return
    if(any(shape(reaction_source)/=[nv-4,n])) return
    if(.not.all(ieee_is_finite([dt,rtol,atol])).or.min(dt,rtol,atol)<=0) return
    candidate=old_state;ct=old_temperature;ir=reaction_source
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      if(geometry%cell_type(i,j,k)==eb_covered_cell_3d) cycle
      call react_sdc_state_3d(species,reactions,old_state(:,i:i,j:j,k:k),old_temperature(i:i,j:j,k:k), &
        nonreacting_rates(:,i:i,j:j,k:k),dt,rtol,atol,candidate(:,i:i,j:j,k:k),ct(i:i,j:j,k:k), &
        ir(:,i:i,j:j,k:k),valid, &
          evaluate_reaction=evaluate_reaction)
      if(.not.valid) return
    end do
    end do
    end do
    state=candidate;temperature=ct;reaction_source=ir;ok=.true.
  end subroutine react_sdc_eb_state_3d

  logical function valid_ctu_geometry(physical,extended,solver) result(valid)
    type(eb_geometry_3d), intent(in) :: physical
    type(eb_geometry_3d), optional, intent(in) :: extended
    character(len=*), intent(in) :: solver
    type(eb_geometry_3d) :: cropped
    integer :: n(3)
    real(dp), parameter :: tolerance=256*epsilon(1._dp)
    real(dp) :: bounds(6),other(6)
    valid=.true.
    if(.not.present(extended)) return
    valid=.false.
    if(solver/='pelec') return
    n=[physical%nx,physical%ny,physical%nz]
    if(any([extended%nx,extended%ny,extended%nz]/=n+12)) return
    call extract_eb_geometry_3d(extended,[7,7,7],n+6,cropped,valid)
    if(.not.valid) return
    valid=.false.
    bounds=[physical%x_lower,physical%x_upper,physical%y_lower,physical%y_upper,physical%z_lower,physical%z_upper]
    other=[cropped%x_lower,cropped%x_upper,cropped%y_lower,cropped%y_upper,cropped%z_lower,cropped%z_upper]
    if(any(abs(bounds-other)>tolerance*max(1._dp,abs(bounds)))) return
    if(any(cropped%cell_type/=physical%cell_type)) return
    if(any(abs(cropped%volume_fraction-physical%volume_fraction)> &
      tolerance*max(1._dp,abs(physical%volume_fraction)))) return
    if(any(abs(cropped%cell_centroid_x-physical%cell_centroid_x)> &
      tolerance*max(1._dp,abs(physical%cell_centroid_x)))) return
    if(any(abs(cropped%cell_centroid_y-physical%cell_centroid_y)> &
      tolerance*max(1._dp,abs(physical%cell_centroid_y)))) return
    if(any(abs(cropped%cell_centroid_z-physical%cell_centroid_z)> &
      tolerance*max(1._dp,abs(physical%cell_centroid_z)))) return
    if(any(abs(cropped%x_face_fraction-physical%x_face_fraction)> &
      tolerance*max(1._dp,abs(physical%x_face_fraction)))) return
    if(any(abs(cropped%y_face_fraction-physical%y_face_fraction)> &
      tolerance*max(1._dp,abs(physical%y_face_fraction)))) return
    if(any(abs(cropped%z_face_fraction-physical%z_face_fraction)> &
      tolerance*max(1._dp,abs(physical%z_face_fraction)))) return
    if(any(abs(cropped%x_face_centroid_y-physical%x_face_centroid_y)> &
      tolerance*max(1._dp,abs(physical%x_face_centroid_y)))) return
    if(any(abs(cropped%x_face_centroid_z-physical%x_face_centroid_z)> &
      tolerance*max(1._dp,abs(physical%x_face_centroid_z)))) return
    if(any(abs(cropped%y_face_centroid_x-physical%y_face_centroid_x)> &
      tolerance*max(1._dp,abs(physical%y_face_centroid_x)))) return
    if(any(abs(cropped%y_face_centroid_z-physical%y_face_centroid_z)> &
      tolerance*max(1._dp,abs(physical%y_face_centroid_z)))) return
    if(any(abs(cropped%z_face_centroid_x-physical%z_face_centroid_x)> &
      tolerance*max(1._dp,abs(physical%z_face_centroid_x)))) return
    if(any(abs(cropped%z_face_centroid_y-physical%z_face_centroid_y)> &
      tolerance*max(1._dp,abs(physical%z_face_centroid_y)))) return
    if(any(abs(cropped%boundary_area-physical%boundary_area)> &
      tolerance*max(1._dp,abs(physical%boundary_area)))) return
    if(any(abs(cropped%boundary_centroid_x-physical%boundary_centroid_x)> &
      tolerance*max(1._dp,abs(physical%boundary_centroid_x)))) return
    if(any(abs(cropped%boundary_centroid_y-physical%boundary_centroid_y)> &
      tolerance*max(1._dp,abs(physical%boundary_centroid_y)))) return
    if(any(abs(cropped%boundary_centroid_z-physical%boundary_centroid_z)> &
      tolerance*max(1._dp,abs(physical%boundary_centroid_z)))) return
    if(any(abs(cropped%boundary_normal_x-physical%boundary_normal_x)> &
      tolerance*max(1._dp,abs(physical%boundary_normal_x)))) return
    if(any(abs(cropped%boundary_normal_y-physical%boundary_normal_y)> &
      tolerance*max(1._dp,abs(physical%boundary_normal_y)))) return
    if(any(abs(cropped%boundary_normal_z-physical%boundary_normal_z)> &
      tolerance*max(1._dp,abs(physical%boundary_normal_z)))) return
    if(any(abs(cropped%boundary_normal_integral_x-physical%boundary_normal_integral_x)> &
      tolerance*max(1._dp,abs(physical%boundary_normal_integral_x)))) return
    if(any(abs(cropped%boundary_normal_integral_y-physical%boundary_normal_integral_y)> &
      tolerance*max(1._dp,abs(physical%boundary_normal_integral_y)))) return
    if(any(abs(cropped%boundary_normal_integral_z-physical%boundary_normal_integral_z)> &
      tolerance*max(1._dp,abs(physical%boundary_normal_integral_z)))) return
    valid=.true.
  end function valid_ctu_geometry

  subroutine initialize_reactive_eb_density_sheet_3d( &
      species, config, geometry, state, temperature, ok, &
      base_mole_fractions)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_eb_3d_config), intent(in) :: config
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(out) :: state(:, :, :, :), temperature(:, :, :)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: base_mole_fractions(:)

    real(dp), parameter :: default_mole_fractions(7) = [ &
      0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
      1.0e-5_dp, 0.0_dp, 0.55643_dp]
    real(dp), allocatable :: primitive(:), mass_fractions(:)
    real(dp), allocatable :: mole_fractions(:)
    real(dp), allocatable :: regular_state(:), cut_state(:)
    real(dp) :: regular_temperature, cut_temperature, sound_speed
    logical :: local_ok
    integer :: i, j, k, nvar, species_index

    state = 0.0_dp
    temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        size(state, 4) /= geometry%nz .or. &
        any(shape(temperature) /= &
          [geometry%nx, geometry%ny, geometry%nz])) return

    allocate(primitive(reactive_nprim(size(species))))
    allocate(mass_fractions(size(species)))
    allocate(mole_fractions(size(species)))
    allocate(regular_state(nvar), cut_state(nvar))
    if (present(base_mole_fractions)) then
      if (size(base_mole_fractions) /= size(species)) return
      mole_fractions = base_mole_fractions
    else
      if (size(species) /= size(default_mole_fractions)) return
      mole_fractions = default_mole_fractions
    end if
    call mass_fractions_from_mole_fractions( &
      species, mole_fractions, mass_fractions, local_ok)
    if (.not. local_ok) return
    primitive(1:5) = [ &
      config%regular_density, config%initial_velocity_x, &
      config%initial_velocity_y, config%initial_velocity_z, &
      config%initial_pressure]
    do species_index = 1, size(species)
      primitive(reactive_mass_fraction_component(species_index)) = &
        mass_fractions(species_index)
    end do
    call reactive_primitive_to_conserved( &
      species, primitive, regular_state, regular_temperature, sound_speed, &
      local_ok)
    if (.not. local_ok) return
    primitive(1) = config%cut_density
    call reactive_primitive_to_conserved( &
      species, primitive, cut_state, cut_temperature, sound_speed, local_ok)
    if (.not. local_ok) return

    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_cut_cell_3d) then
            state(:, i, j, k) = cut_state
            temperature(i, j, k) = cut_temperature
          else
            state(:, i, j, k) = regular_state
            temperature(i, j, k) = regular_temperature
          end if
        end do
      end do
    end do
    ok = all(ieee_is_finite(state)) .and. &
      all(ieee_is_finite(temperature))
  end subroutine initialize_reactive_eb_density_sheet_3d

  subroutine advance_reactive_eb_strang_3d( &
      species, reactions, state, temperature, geometry, solver, redistribution, &
      dt, chemistry_enabled, rtol, atol, new_state, new_temperature, ok, &
      target_volume_fraction, chemistry_integrator, boundary_change, wall_change, chemistry_change, difmag,extended_geometry, &
      boundary_faces,inflow_state,inflow_temperature)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    type(eb_geometry_3d), optional, intent(in) :: extended_geometry
    character(len=*), intent(in) :: solver, redistribution
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:)
    character(len=32) :: hydro_faces(6)
    real(dp), intent(in) :: dt, rtol, atol
    logical, intent(in) :: chemistry_enabled
    real(dp), intent(out) :: new_state(:, :, :, :), new_temperature(:, :, :)
    real(dp), optional, intent(in) :: difmag
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: target_volume_fraction
    character(len=*), intent(in), optional :: chemistry_integrator
    real(dp), optional, intent(inout) :: boundary_change(:),wall_change(:),chemistry_change(:)
    real(dp), allocatable :: exchange_boundary(:),exchange_wall(:),exchange_chemistry(:),before(:),after(:)

    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    real(dp), allocatable :: hydro_state(:, :, :, :)
    real(dp), allocatable :: hydro_temperature(:, :, :)
    logical, allocatable :: active_mask(:, :, :)
    real(dp) :: target
    logical :: local_ok
    integer :: nvar

    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (size(state, 1) /= nvar .or. size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        size(state, 4) /= geometry%nz .or. &
        any(shape(temperature) /= [geometry%nx, geometry%ny, geometry%nz]) .or. &
        any(shape(new_state) /= shape(state)) .or. &
        any(shape(new_temperature) /= shape(temperature))) return
    new_state = state
    new_temperature = temperature
    if(.not.valid_ctu_geometry(geometry,extended_geometry,solver)) return
    hydro_faces='outflow'
    if(present(boundary_faces)) then
      if(.not.present(extended_geometry)) return
      if(.not.valid_reactive_boundaries_3d(boundary_faces)) return
      hydro_faces=boundary_faces
      if(.not.periodic_eb_faces_match(extended_geometry,hydro_faces,6)) return
    else
      if(present(inflow_state).or.present(inflow_temperature)) return
    end if
    if(present(difmag)) then
      if(.not.ieee_is_finite(difmag).or.difmag<0) return
    end if
    if (present(boundary_change)) then
      if (size(boundary_change)/=nvar) return
    end if
    if (present(wall_change)) then
      if (size(wall_change)/=nvar) return
    end if
    if (present(chemistry_change)) then
      if (size(chemistry_change)/=nvar) return
    end if
    if (present(boundary_change).or.present(wall_change).or.present(chemistry_change)) then
      allocate(exchange_boundary(nvar),exchange_wall(nvar),exchange_chemistry(nvar),before(nvar),after(nvar))
      exchange_chemistry=0._dp
    end if
    if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp .or. &
        .not. ieee_is_finite(rtol) .or. .not. ieee_is_finite(atol) .or. &
        rtol <= 0.0_dp .or. atol <= 0.0_dp) return
    if (chemistry_enabled .and. size(reactions) < 1) return
    target = 0.5_dp
    if (present(target_volume_fraction)) target = target_volume_fraction
    if (.not. ieee_is_finite(target) .or. target <= 0.0_dp .or. &
        target > 1.0_dp) return
    select case (trim(redistribution))
    case ("flux_redist", "state_redist")
    case default
      return
    end select

    allocate(candidate_state, source=state)
    allocate(candidate_temperature, source=temperature)
    allocate(hydro_state, source=state)
    allocate(hydro_temperature, source=temperature)
    allocate(active_mask(geometry%nx, geometry%ny, geometry%nz))
    active_mask = geometry%cell_type /= eb_covered_cell_3d

    if (chemistry_enabled) then
      if (allocated(exchange_chemistry)) then
        call reactive_eb_integrals_3d(candidate_state,geometry,before,local_ok)
        if (.not.local_ok) return
      end if
      call advance_reactive_chemistry_3d( &
        species, reactions, candidate_state, candidate_temperature, &
        geometry%nx, geometry%ny, geometry%nz, 0.5_dp * dt, rtol, atol, &
        local_ok, active_mask, chemistry_integrator)
      if (.not. local_ok) return
      call covered_cells_match_3d( &
        state, temperature, candidate_state, candidate_temperature, geometry, &
        local_ok)
      if (.not. local_ok) return
      if (allocated(exchange_chemistry)) then
        call reactive_eb_integrals_3d(candidate_state,geometry,after,local_ok)
        if (.not.local_ok) return
        exchange_chemistry=exchange_chemistry+after-before
      end if
    end if

    if(present(extended_geometry)) then
      call advance_reactive_eb_boundary_ctu_3d(species,candidate_state,candidate_temperature,extended_geometry,dt, &
        hydro_state,hydro_temperature,local_ok,hydro_faces,redistribution,difmag,boundary_change=exchange_boundary, &
        wall_change=exchange_wall,target_volume_fraction=target,inflow_state=inflow_state,inflow_temperature=inflow_temperature)
    else
      select case (trim(redistribution))
      case ("flux_redist")
        call advance_reactive_eb_redistributed_euler_3d( &
          species, candidate_state, candidate_temperature, geometry, solver, dt, &
          hydro_state, hydro_temperature, local_ok,exchange_boundary,exchange_wall,difmag=difmag)
      case ("state_redist")
        call advance_reactive_eb_state_redistributed_euler_3d( &
          species, candidate_state, candidate_temperature, geometry, solver, dt, &
          hydro_state, hydro_temperature, local_ok, target,exchange_boundary,exchange_wall,difmag=difmag)
      end select
    end if
    if (.not. local_ok) return
    call covered_cells_match_3d( &
      state, temperature, hydro_state, hydro_temperature, geometry, local_ok)
    if (.not. local_ok) return
    candidate_state = hydro_state
    candidate_temperature = hydro_temperature

    if (chemistry_enabled) then
      if (allocated(exchange_chemistry)) then
        call reactive_eb_integrals_3d(candidate_state,geometry,before,local_ok)
        if (.not.local_ok) return
      end if
      call advance_reactive_chemistry_3d( &
        species, reactions, candidate_state, candidate_temperature, &
        geometry%nx, geometry%ny, geometry%nz, 0.5_dp * dt, rtol, atol, &
        local_ok, active_mask, chemistry_integrator)
      if (.not. local_ok) return
      call covered_cells_match_3d( &
        state, temperature, candidate_state, candidate_temperature, geometry, &
        local_ok)
      if (.not. local_ok) return
      if (allocated(exchange_chemistry)) then
        call reactive_eb_integrals_3d(candidate_state,geometry,after,local_ok)
        if (.not.local_ok) return
        exchange_chemistry=exchange_chemistry+after-before
      end if
    end if
    if (any(.not. ieee_is_finite(candidate_state)) .or. &
        any(.not. ieee_is_finite(candidate_temperature))) return
    if (allocated(exchange_chemistry)) then
      if (.not.all(ieee_is_finite(exchange_chemistry))) return
    end if
    if (present(boundary_change)) boundary_change=exchange_boundary
    if (present(wall_change)) wall_change=exchange_wall
    if (present(chemistry_change)) chemistry_change=exchange_chemistry
    new_state = candidate_state
    new_temperature = candidate_temperature
    ok = .true.
  end subroutine advance_reactive_eb_strang_3d

  subroutine advance_reactive_eb_full_3d( &
      species, reactions, transport, state, temperature, geometry, solver, &
      redistribution, dt, chemistry_enabled, rtol, atol, &
      transport_enabled, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, new_state, &
      new_temperature, minimum_transport_theta, ok, target_volume_fraction, &
      chemistry_integrator, boundary_change, wall_change, chemistry_change, difmag,extended_geometry, &
      boundary_faces,inflow_state,inflow_temperature,prepare_transport_halo,wall_temperature)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    type(eb_geometry_3d), optional, intent(in) :: extended_geometry
    character(len=*), intent(in) :: solver, redistribution
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_temperature(6)
    character(len=32) :: hydro_faces(6)
    procedure(eb_transport_halo_provider_3d), optional :: prepare_transport_halo
    type(eb_integrated_flux_3d), allocatable :: transport_flux
    real(dp), allocatable :: transport_change(:),half_change(:)
    character(len=32), allocatable :: transport_faces(:)
    real(dp), allocatable :: transport_inflow(:,:),transport_inflow_temperature(:)
    real(dp), intent(in) :: dt, rtol, atol
    logical, intent(in) :: chemistry_enabled, transport_enabled
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: new_state(:, :, :, :), new_temperature(:, :, :)
    real(dp), intent(out) :: minimum_transport_theta
    real(dp), optional, intent(in) :: difmag
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: target_volume_fraction
    character(len=*), intent(in), optional :: chemistry_integrator
    real(dp), optional, intent(inout) :: boundary_change(:),wall_change(:),chemistry_change(:)
    real(dp), allocatable :: exchange_boundary(:),exchange_wall(:),exchange_chemistry(:),before(:),after(:)

    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    real(dp), allocatable :: hydro_state(:, :, :, :)
    real(dp), allocatable :: hydro_temperature(:, :, :)
    logical, allocatable :: active_mask(:, :, :)
    real(dp) :: target, theta_one, theta_two
    logical :: local_ok
    integer :: nvar

    minimum_transport_theta = 1.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        size(state, 4) /= geometry%nz .or. &
        any(shape(temperature) /= &
          [geometry%nx, geometry%ny, geometry%nz]) .or. &
        any(shape(new_state) /= shape(state)) .or. &
        any(shape(new_temperature) /= shape(temperature))) return
    new_state = state
    new_temperature = temperature
    if(.not.valid_ctu_geometry(geometry,extended_geometry,solver)) return
    hydro_faces='outflow'
    if(present(boundary_faces)) then
      if(.not.present(extended_geometry)) return
      if(.not.valid_reactive_boundaries_3d(boundary_faces)) return
      hydro_faces=boundary_faces
      if(.not.periodic_eb_faces_match(extended_geometry,hydro_faces,6)) return
    else
      if(present(inflow_state).or.present(inflow_temperature)) return
    end if
    ! A custom stage provider must identify the matching hydro boundaries.
    ! Without one, the common boundary filler supplies inflow/outflow donors.
    if(present(prepare_transport_halo).and..not.present(boundary_faces)) return
    if(present(difmag)) then
      if(.not.ieee_is_finite(difmag).or.difmag<0) return
    end if
    if (present(boundary_change)) then
      if (size(boundary_change)/=nvar) return
    end if
    if (present(wall_change)) then
      if (size(wall_change)/=nvar) return
    end if
    if (present(chemistry_change)) then
      if (size(chemistry_change)/=nvar) return
    end if
    if (present(boundary_change).or.present(wall_change).or.present(chemistry_change)) then
      allocate(exchange_boundary(nvar),exchange_wall(nvar),exchange_chemistry(nvar),before(nvar),after(nvar))
      exchange_chemistry=0._dp
    end if
    if (.not. transport_enabled) then
      if(present(wall_temperature)) return
      call advance_reactive_eb_strang_3d( &
        species, reactions, state, temperature, geometry, solver, &
        redistribution, dt, chemistry_enabled, rtol, atol, new_state, &
        new_temperature, ok, target_volume_fraction, chemistry_integrator, &
        boundary_change,wall_change,chemistry_change,difmag=difmag,extended_geometry=extended_geometry, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_temperature=inflow_temperature)
      return
    end if

    if(present(boundary_faces)) then
      if(present(prepare_transport_halo)) then
        if(any(hydro_faces/='outflow'.and.hydro_faces/='inflow')) return
      end if
      if(any(hydro_faces=='isothermal_wall')) then
        if(.not.thermal_conduction_enabled.or..not.present(wall_temperature)) return
      end if
    else
      if(present(wall_temperature)) return
    end if
    if (size(transport) /= size(species)) return
    if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp .or. &
        .not. ieee_is_finite(rtol) .or. .not. ieee_is_finite(atol) .or. &
        rtol <= 0.0_dp .or. atol <= 0.0_dp) return
    if (trim(redistribution) /= "state_redist") return
    if (.not. (viscosity_enabled .or. thermal_conduction_enabled .or. &
        species_diffusion_enabled)) return
    if (barodiffusion_enabled .and. .not. species_diffusion_enabled) return
    if (chemistry_enabled .and. size(reactions) < 1) return
    target = 0.5_dp
    if (present(target_volume_fraction)) target = target_volume_fraction
    if (.not. ieee_is_finite(target) .or. target <= 0.0_dp .or. &
        target > 1.0_dp) return

    if(present(boundary_faces).and..not.present(prepare_transport_halo)) then
      allocate(transport_faces,source=hydro_faces)
      if(present(inflow_state)) allocate(transport_inflow,source=inflow_state)
      if(present(inflow_temperature)) allocate(transport_inflow_temperature,source=inflow_temperature)
    end if
    if(present(boundary_faces).and.allocated(exchange_boundary)) then
      allocate(transport_flux,transport_change(nvar),half_change(nvar))
      transport_change=0
    end if
    allocate(candidate_state, source=state)
    allocate(candidate_temperature, source=temperature)
    allocate(hydro_state, source=state)
    allocate(hydro_temperature, source=temperature)
    allocate(active_mask(geometry%nx, geometry%ny, geometry%nz))
    active_mask = geometry%cell_type /= eb_covered_cell_3d

    if (chemistry_enabled) then
      if (allocated(exchange_chemistry)) then
        call reactive_eb_integrals_3d(candidate_state,geometry,before,local_ok)
        if (.not.local_ok) return
      end if
      call advance_reactive_chemistry_3d( &
        species, reactions, candidate_state, candidate_temperature, &
        geometry%nx, geometry%ny, geometry%nz, 0.5_dp * dt, rtol, atol, &
        local_ok, active_mask, chemistry_integrator)
      if (.not. local_ok) return
      call covered_cells_match_3d( &
        state, temperature, candidate_state, candidate_temperature, &
        geometry, local_ok)
      if (.not. local_ok) return
      if (allocated(exchange_chemistry)) then
        call reactive_eb_integrals_3d(candidate_state,geometry,after,local_ok)
        if (.not.local_ok) return
        exchange_chemistry=exchange_chemistry+after-before
      end if
    end if

    call advance_reactive_eb_transport_3d( &
      species, transport, candidate_state, candidate_temperature, &
      geometry, 0.5_dp * dt, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, target, theta_one, local_ok,prepare_transport_halo,transport_flux, &
      transport_faces,transport_inflow,transport_inflow_temperature,wall_temperature)
    if (.not. local_ok) return
    if(allocated(transport_flux)) then
      call transport_boundary_inventory(geometry,transport_flux,half_change,local_ok)
      if(.not.local_ok) return
      transport_change=transport_change+half_change
    end if
    call covered_cells_match_3d( &
      state, temperature, candidate_state, candidate_temperature, &
      geometry, local_ok)
    if (.not. local_ok) return

    if(present(extended_geometry)) then
      call advance_reactive_eb_boundary_ctu_3d(species,candidate_state,candidate_temperature,extended_geometry,dt, &
        hydro_state,hydro_temperature,local_ok,hydro_faces,redistribution,difmag,boundary_change=exchange_boundary, &
        wall_change=exchange_wall,target_volume_fraction=target,inflow_state=inflow_state,inflow_temperature=inflow_temperature)
    else
      call advance_reactive_eb_state_redistributed_euler_3d( &
        species, candidate_state, candidate_temperature, geometry, solver, &
        dt, hydro_state, hydro_temperature, local_ok, target,exchange_boundary,exchange_wall,difmag=difmag)
    end if
    if (.not. local_ok) return
    call covered_cells_match_3d( &
      state, temperature, hydro_state, hydro_temperature, geometry, &
      local_ok)
    if (.not. local_ok) return
    candidate_state = hydro_state
    candidate_temperature = hydro_temperature

    call advance_reactive_eb_transport_3d( &
      species, transport, candidate_state, candidate_temperature, &
      geometry, 0.5_dp * dt, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, target, theta_two, local_ok,prepare_transport_halo,transport_flux, &
      transport_faces,transport_inflow,transport_inflow_temperature,wall_temperature)
    if (.not. local_ok) return
    if(allocated(transport_flux)) then
      call transport_boundary_inventory(geometry,transport_flux,half_change,local_ok)
      if(.not.local_ok) return
      transport_change=transport_change+half_change
    end if
    call covered_cells_match_3d( &
      state, temperature, candidate_state, candidate_temperature, &
      geometry, local_ok)
    if (.not. local_ok) return

    if (chemistry_enabled) then
      if (allocated(exchange_chemistry)) then
        call reactive_eb_integrals_3d(candidate_state,geometry,before,local_ok)
        if (.not.local_ok) return
      end if
      call advance_reactive_chemistry_3d( &
        species, reactions, candidate_state, candidate_temperature, &
        geometry%nx, geometry%ny, geometry%nz, 0.5_dp * dt, rtol, atol, &
        local_ok, active_mask, chemistry_integrator)
      if (.not. local_ok) return
      call covered_cells_match_3d( &
        state, temperature, candidate_state, candidate_temperature, &
        geometry, local_ok)
      if (.not. local_ok) return
      if (allocated(exchange_chemistry)) then
        call reactive_eb_integrals_3d(candidate_state,geometry,after,local_ok)
        if (.not.local_ok) return
        exchange_chemistry=exchange_chemistry+after-before
      end if
    end if
    if (any(.not. ieee_is_finite(candidate_state)) .or. &
        any(.not. ieee_is_finite(candidate_temperature))) return

    if (allocated(exchange_chemistry)) then
      if (.not.all(ieee_is_finite(exchange_chemistry))) return
    end if
    if(allocated(transport_change)) then
      exchange_boundary=exchange_boundary+transport_change
      if(.not.all(ieee_is_finite(exchange_boundary))) return
    end if
    if (present(boundary_change)) boundary_change=exchange_boundary
    if (present(wall_change)) wall_change=exchange_wall
    if (present(chemistry_change)) chemistry_change=exchange_chemistry
    new_state = candidate_state
    new_temperature = candidate_temperature
    minimum_transport_theta = min(theta_one, theta_two)
    ok = .true.
  end subroutine advance_reactive_eb_full_3d

  ! Flux is already integrated in time and per open face area.
  subroutine transport_boundary_inventory(geometry,flux,change,ok)
    type(eb_geometry_3d), intent(in) :: geometry
    type(eb_integrated_flux_3d), intent(in) :: flux
    real(dp), intent(out) :: change(:)
    logical, intent(out) :: ok
    integer :: c,n(3)
    ok=.false.
    if(.not.valid_eb_integrated_flux_3d(flux,geometry,size(change))) return
    n=[geometry%nx,geometry%ny,geometry%nz]
    do c=1,size(change)
      change(c)=geometry%dy*geometry%dz*sum( &
        geometry%x_face_fraction(0,:,:)*flux%x(c,0,:,:)-geometry%x_face_fraction(n(1),:,:)*flux%x(c,n(1),:,:)) &
        +geometry%dx*geometry%dz*sum( &
        geometry%y_face_fraction(:,0,:)*flux%y(c,:,0,:)-geometry%y_face_fraction(:,n(2),:)*flux%y(c,:,n(2),:)) &
        +geometry%dx*geometry%dy*sum( &
        geometry%z_face_fraction(:,:,0)*flux%z(c,:,:,0)-geometry%z_face_fraction(:,:,n(3))*flux%z(c,:,:,n(3)))
    end do
    ok=all(ieee_is_finite(change))
  end subroutine transport_boundary_inventory

  subroutine covered_cells_match_3d( &
      reference_state, reference_temperature, candidate_state, &
      candidate_temperature, geometry, ok)
    real(dp), intent(in) :: reference_state(:, :, :, :)
    real(dp), intent(in) :: reference_temperature(:, :, :)
    real(dp), intent(in) :: candidate_state(:, :, :, :)
    real(dp), intent(in) :: candidate_temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    logical, intent(out) :: ok

    integer :: i, j, k

    ok = .false.
    if (any(shape(candidate_state) /= shape(reference_state)) .or. &
        any(shape(candidate_temperature) /= shape(reference_temperature))) return
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) /= eb_covered_cell_3d) cycle
          if (.not. all(same_real_bits(candidate_state(:, i, j, k), &
              reference_state(:, i, j, k))) .or. &
              .not. same_real_bits(candidate_temperature(i, j, k), &
                reference_temperature(i, j, k))) return
        end do
      end do
    end do
    ok = .true.
  end subroutine covered_cells_match_3d

  pure elemental logical function same_real_bits(left, right)
    real(dp), intent(in) :: left, right

    same_real_bits = transfer(left, 0_int64) == transfer(right, 0_int64)
  end function same_real_bits

  subroutine reactive_eb_integrals_3d( &
      field, geometry, integrals, ok, l1_integrals)
    real(dp), intent(in) :: field(:, :, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(out) :: integrals(:)
    logical, intent(out) :: ok
    real(dp), intent(out), optional :: l1_integrals(:)

    real(dp) :: cell_volume
    integer :: component

    integrals = 0.0_dp
    if (present(l1_integrals)) l1_integrals = 0.0_dp
    ok = .false.
    if (.not. geometry%is_valid() .or. &
        size(field, 2) /= geometry%nx .or. &
        size(field, 3) /= geometry%ny .or. &
        size(field, 4) /= geometry%nz .or. &
        size(integrals) /= size(field, 1) .or. &
        any(.not. ieee_is_finite(field))) return
    if (present(l1_integrals)) then
      if (size(l1_integrals) /= size(field, 1)) return
    end if
    cell_volume = geometry%dx * geometry%dy * geometry%dz
    do component = 1, size(field, 1)
      integrals(component) = sum(geometry%volume_fraction * &
        field(component, :, :, :)) * cell_volume
      if (present(l1_integrals)) then
        l1_integrals(component) = sum(geometry%volume_fraction * &
          abs(field(component, :, :, :))) * cell_volume
      end if
    end do
    ok = all(ieee_is_finite(integrals))
    if (present(l1_integrals)) then
      ok = ok .and. all(ieee_is_finite(l1_integrals))
    end if
  end subroutine reactive_eb_integrals_3d

  subroutine reactive_eb_composition_integrals_3d(species,state,geometry,element_atoms,integrals,ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),element_atoms(:,:)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(inout) :: integrals(:)
    logical, intent(out) :: ok
    real(dp) :: candidate(size(integrals)),amount,volume
    integer :: i,j,k,s,nvar
    ok=.false.
    nvar=reactive_nvar(size(species))
    if (size(integrals)<1.or.size(species)<1.or..not.geometry%is_valid()) return
    if (any(shape(element_atoms)/=[size(integrals),size(species)])) return
    if (any(shape(state)/=[nvar,geometry%nx,geometry%ny,geometry%nz])) return
    if (.not.all(ieee_is_finite(element_atoms)).or.any(element_atoms<0._dp)) return
    if (any(maxval(element_atoms,dim=1)<=0._dp)) return
    if (.not.all(ieee_is_finite(species%molecular_weight))) return
    if (any(species%molecular_weight<=0._dp)) return
    candidate=0._dp;volume=geometry%dx*geometry%dy*geometry%dz
    do k=1,geometry%nz
      do j=1,geometry%ny
        do i=1,geometry%nx
          if (geometry%cell_type(i,j,k)==eb_covered_cell_3d) cycle
          if (.not.all(ieee_is_finite(state(:,i,j,k)))) return
          do s=1,size(species)
            amount=state(reactive_species_component(s),i,j,k)
            if (amount<0._dp) return
            amount=amount/species(s)%molecular_weight*geometry%volume_fraction(i,j,k)*volume
            candidate=candidate+element_atoms(:,s)*amount
          end do
        end do
      end do
    end do
    if (.not.all(ieee_is_finite(candidate))) return
    integrals=candidate;ok=.true.
  end subroutine reactive_eb_composition_integrals_3d

  subroutine reactive_eb_element_integrals_3d( &
      species, state, geometry, element_integrals, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(out) :: element_integrals(:)
    logical, intent(out) :: ok

    real(dp) :: atoms(3), cell_volume
    logical :: local_ok
    integer :: i, j, k, species_index, nvar

    element_integrals = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. size(element_integrals) /= 3 .or. &
        .not. geometry%is_valid() .or. size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. size(state, 3) /= geometry%ny .or. &
        size(state, 4) /= geometry%nz) return
    cell_volume = geometry%dx * geometry%dy * geometry%dz
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          if (any(.not. ieee_is_finite(state(:, i, j, k)))) return
          do species_index = 1, size(species)
            call element_atom_counts( &
              species(species_index)%name, atoms, local_ok)
            if (.not. local_ok .or. &
                .not. ieee_is_finite(species(species_index)%molecular_weight) .or. &
                species(species_index)%molecular_weight <= 0.0_dp) return
            element_integrals = element_integrals + &
              geometry%volume_fraction(i, j, k) * &
              state(reactive_species_component(species_index), i, j, k) * &
              atoms / species(species_index)%molecular_weight * cell_volume
          end do
        end do
      end do
    end do
    ok = all(ieee_is_finite(element_integrals))
  end subroutine reactive_eb_element_integrals_3d

  pure logical function reactive_eb_element_species_supported_3d( &
      species) result(supported)
    type(nasa7_species), intent(in) :: species(:)

    real(dp) :: atoms(3)
    logical :: ok
    integer :: species_index

    supported = size(species) >= 1
    do species_index = 1, size(species)
      call element_atom_counts(species(species_index)%name, atoms, ok)
      if (.not. ok) then
        supported = .false.
        return
      end if
    end do
  end function reactive_eb_element_species_supported_3d

  subroutine reactive_eb_element_amounts_3d(species,amounts,elements,ok,element_atoms)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: amounts(:)
    real(dp), intent(inout) :: elements(:)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: element_atoms(:,:)
    real(dp) :: candidate(size(elements)),atoms(3)
    integer :: k
    ok=.false.
    if (size(species)<1.or.size(elements)<1) return
    if (size(amounts)/=reactive_nvar(size(species))) return
    if (.not.all(ieee_is_finite(amounts))) return
    if (.not.all(ieee_is_finite(species%molecular_weight))) return
    if (any(species%molecular_weight<=0._dp)) return
    if (present(element_atoms)) then
      if (any(shape(element_atoms)/=[size(elements),size(species)])) return
      if (.not.all(ieee_is_finite(element_atoms))) return
      if (any(element_atoms<0._dp)) return
      if (any(sum(element_atoms,dim=1)<=0._dp)) return
      candidate=matmul(element_atoms,amounts(6:)/species%molecular_weight)
    else
      if (size(elements)/=3) return
      candidate=0._dp
      do k=1,size(species)
        call element_atom_counts(species(k)%name,atoms,ok)
        if (.not.ok) return
        candidate=candidate+atoms*amounts(k+5)/species(k)%molecular_weight
      end do
    end if
    ok=.false.
    if (.not.all(ieee_is_finite(candidate))) return
    elements=candidate
    ok=.true.
  end subroutine reactive_eb_element_amounts_3d

  pure subroutine element_atom_counts(name, atoms, ok)
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: atoms(3)
    logical, intent(out) :: ok

    atoms = 0.0_dp
    select case (trim(name))
    case ("H2")
      atoms = [2.0_dp, 0.0_dp, 0.0_dp]
    case ("H")
      atoms = [1.0_dp, 0.0_dp, 0.0_dp]
    case ("O")
      atoms = [0.0_dp, 1.0_dp, 0.0_dp]
    case ("O2")
      atoms = [0.0_dp, 2.0_dp, 0.0_dp]
    case ("OH")
      atoms = [1.0_dp, 1.0_dp, 0.0_dp]
    case ("H2O")
      atoms = [2.0_dp, 1.0_dp, 0.0_dp]
    case ("HO2")
      atoms = [1.0_dp, 2.0_dp, 0.0_dp]
    case ("H2O2")
      atoms = [2.0_dp, 2.0_dp, 0.0_dp]
    case ("N2")
      atoms = [0.0_dp, 0.0_dp, 2.0_dp]
    case ("AR")
      atoms = 0.0_dp
    case default
      ok = .false.
      return
    end select
    ok = .true.
  end subroutine element_atom_counts

  subroutine reactive_eb_extrema_3d( &
      species, state, temperature, geometry, minimum_density, &
      maximum_density, minimum_pressure, maximum_pressure, &
      minimum_temperature, maximum_temperature, maximum_speed, &
      maximum_closure_error, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(out) :: minimum_density, maximum_density
    real(dp), intent(out) :: minimum_pressure, maximum_pressure
    real(dp), intent(out) :: minimum_temperature, maximum_temperature
    real(dp), intent(out) :: maximum_speed, maximum_closure_error
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: local_temperature, sound_speed, speed, closure_error
    logical :: local_ok
    integer :: i, j, k, species_index, active_cells, nvar

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
    if (nvar <= 0 .or. .not. geometry%is_valid() .or. &
        size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        size(state, 4) /= geometry%nz .or. &
        any(shape(temperature) /= &
          [geometry%nx, geometry%ny, geometry%nz])) return

    allocate(primitive(reactive_nprim(size(species))))
    active_cells = 0
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          active_cells = active_cells + 1
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), &
            primitive, local_temperature, sound_speed, local_ok)
          if (.not. local_ok) return
          minimum_density = min(minimum_density, primitive(1))
          maximum_density = max(maximum_density, primitive(1))
          minimum_pressure = min(minimum_pressure, primitive(5))
          maximum_pressure = max(maximum_pressure, primitive(5))
          minimum_temperature = min(minimum_temperature, local_temperature)
          maximum_temperature = max(maximum_temperature, local_temperature)
          speed = sqrt(sum(primitive(2:4)**2))
          maximum_speed = max(maximum_speed, speed)
          closure_error = 0.0_dp
          do species_index = 1, size(species)
            closure_error = closure_error + &
              primitive(reactive_mass_fraction_component(species_index))
          end do
          maximum_closure_error = max( &
            maximum_closure_error, abs(closure_error - 1.0_dp))
        end do
      end do
    end do
    ok = active_cells > 0 .and. all(ieee_is_finite([ &
      minimum_density, maximum_density, minimum_pressure, &
      maximum_pressure, minimum_temperature, maximum_temperature, &
      maximum_speed, maximum_closure_error]))
  end subroutine reactive_eb_extrema_3d

  subroutine write_reactive_eb_3d_csv( &
      path, species, geometry, state, temperature, time, ok, message)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :), time
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    real(dp), allocatable :: primitive(:)
    real(dp) :: x, y, z, centroid_x, centroid_y, centroid_z
    real(dp) :: local_temperature, sound_speed
    logical :: cell_ok
    integer :: unit, io_status, i, j, k, species_index, nvar

    ok = .false.
    message = ""
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid() .or. &
        .not. ieee_is_finite(time) .or. &
        size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        size(state, 4) /= geometry%nz .or. &
        any(shape(temperature) /= &
          [geometry%nx, geometry%ny, geometry%nz])) then
      message = "Reactive EB 3D CSV input contract is invalid"
      return
    end if
    allocate(primitive(reactive_nprim(size(species))))
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), &
            primitive, local_temperature, sound_speed, cell_ok)
          if (.not. cell_ok) then
            write(message, '(a,i0,a,i0,a,i0,a)') &
              "Invalid output state at (", i, ",", j, ",", k, ")"
            return
          end if
        end do
      end do
    end do

    open(newunit=unit, file=trim(path), status="replace", action="write", &
      iostat=io_status)
    if (io_status /= 0) then
      write(message, '(a,1x,a)') "Could not create output file:", trim(path)
      return
    end if
    write(unit, '(a)', advance='no') &
      "time,i,j,k,x,y,z,cell_type,volume_fraction,fluid_centroid_x," // &
      "fluid_centroid_y,fluid_centroid_z,rho,u,v,w,pressure," // &
      "temperature,rhoE,rhou,rhov,rhow"
    do species_index = 1, size(species)
      write(unit, '(a)', advance='no') &
        ",Y_" // trim(species(species_index)%name)
    end do
    do species_index = 1, size(species)
      write(unit, '(a)', advance='no') &
        ",rhoY_" // trim(species(species_index)%name)
    end do
    write(unit, '(a)') ""

    do k = 1, geometry%nz
      z = geometry%z_lower + (real(k, dp) - 0.5_dp) * geometry%dz
      do j = 1, geometry%ny
        y = geometry%y_lower + (real(j, dp) - 0.5_dp) * geometry%dy
        do i = 1, geometry%nx
          x = geometry%x_lower + (real(i, dp) - 0.5_dp) * geometry%dx
          centroid_x = x + geometry%cell_centroid_x(i, j, k) * geometry%dx
          centroid_y = y + geometry%cell_centroid_y(i, j, k) * geometry%dy
          centroid_z = z + geometry%cell_centroid_z(i, j, k) * geometry%dz
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), &
            primitive, local_temperature, sound_speed, cell_ok)
          if (.not. cell_ok) then
            close(unit)
            message = "Reactive EB 3D state changed during CSV output"
            return
          end if
          write(unit, '(*(es25.16e3,:,","))') &
            time, real(i, dp), real(j, dp), real(k, dp), x, y, z, &
            real(geometry%cell_type(i, j, k), dp), &
            geometry%volume_fraction(i, j, k), &
            centroid_x, centroid_y, centroid_z, primitive(1:5), &
            local_temperature, state(iet, i, j, k), &
            state(imx, i, j, k), state(imy, i, j, k), &
            state(imz, i, j, k), &
            (primitive(reactive_mass_fraction_component(species_index)), &
              species_index = 1, size(species)), &
            (state(reactive_species_component(species_index), i, j, k), &
              species_index = 1, size(species))
        end do
      end do
    end do
    close(unit)
    ok = .true.
  end subroutine write_reactive_eb_3d_csv

end module reactive_eb_3d_driver_mod
