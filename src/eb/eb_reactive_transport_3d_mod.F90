module eb_reactive_transport_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use simulation_config_reactive_3d_mod, only: fill_reactive_boundary_ghosts_3d,valid_reactive_boundaries_3d
  use state_indices_mod, only: iet
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_conserved_to_primitive
  use reactive_transport_3d_mod, only: &
    reactive_transport_face_flux_3d, reactive_transport_timestep_3d,wall_heat_transfer_3d
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d,periodic_eb_faces_match
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_reactive_redistribution_3d_mod, only: &
    advance_reactive_eb_state_redistributed_3d
  implicit none
  private

  ! A stage donor provider fills one halo layer from the current stage state.
  ! Physical labels, prescribed values and patch exchange belong to the caller.
  abstract interface
    subroutine eb_transport_halo_provider_3d(state,temperature,geometry,ghost,ghost_temperature,faces,ok)
      import dp,eb_geometry_3d
      real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:)
      type(eb_geometry_3d), intent(in) :: geometry
      real(dp), intent(out) :: ghost(:,0:,0:,0:),ghost_temperature(0:,0:,0:)
      logical, intent(out) :: faces(6),ok
    end subroutine eb_transport_halo_provider_3d
  end interface
  public :: eb_transport_halo_provider_3d

  real(dp), parameter :: species_safety = 0.90_dp

  public :: reactive_eb_halo_bounds_3d,reactive_eb_isothermal_timestep_3d
  public :: reactive_eb_transport_timestep_3d,reactive_eb_transport_patch_timestep_3d
  public :: reactive_eb_transport_fluxes_rhs_3d
  public :: reactive_eb_transport_rhs_3d
  public :: reactive_eb_transport_euler_update_3d
  public :: advance_reactive_eb_transport_3d,advance_reactive_eb_transport_patch_euler_3d

contains
  ! Rate from the same half-cell heat-transfer law used by wall face fluxes.
  ! Multiplying by aperture/kappa accounts for the actual fluid heat capacity.
  subroutine reactive_eb_isothermal_timestep_3d(species,transport,state,temperature,geometry,faces,wall_t,cfl,dt,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),wall_t(6),cfl
    type(eb_geometry_3d), intent(in) :: geometry
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok
    real(dp) :: spacing(3),heat,rate,area,weighted,cell_rate,maximum_rate,candidate,vf
    integer :: n(3),idx(3),i,j,k,face,axis,side
    logical :: valid
    dt=0;ok=.false.
    if(.not.valid_transport_shapes(species,transport,state,temperature,geometry)) return
    if(.not.valid_reactive_boundaries_3d(faces).or..not.ieee_is_finite(cfl)) return
    if(cfl<=0.or.cfl>.5_dp) return
    n=[geometry%nx,geometry%ny,geometry%nz];spacing=[geometry%dx,geometry%dy,geometry%dz]
    maximum_rate=0
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      if(geometry%cell_type(i,j,k)==eb_covered_cell_3d) cycle
      idx=[i,j,k];cell_rate=0;vf=geometry%volume_fraction(i,j,k)
      if(vf<=0) return
      do face=1,6
        if(faces(face)/='isothermal_wall') cycle
        axis=(face+1)/2;side=modulo(face+1,2)
        if(side==0.and.idx(axis)/=1) cycle
        if(side==1.and.idx(axis)/=n(axis)) cycle
        select case(axis)
        case(1);area=geometry%x_face_fraction(side*n(1),j,k)
        case(2);area=geometry%y_face_fraction(i,side*n(2),k)
        case(3);area=geometry%z_face_fraction(i,j,side*n(3))
        end select
        if(area<=0) cycle
        call wall_heat_transfer_3d(species,transport,state(:,i,j,k),temperature(i,j,k), &
          wall_t(face),spacing(axis),heat,rate,valid)
        if(.not.valid) return
        weighted=rate*area
        if(weighted>huge(1._dp)*vf) return
        weighted=weighted/vf
        if(weighted>huge(1._dp)-cell_rate) return
        cell_rate=cell_rate+weighted
      end do
      maximum_rate=max(maximum_rate,cell_rate)
    end do
    end do
    end do
    candidate=huge(1._dp)
    if(maximum_rate>0) then
      if(maximum_rate<1) then
        if(cfl>huge(1._dp)*maximum_rate) return
      end if
      candidate=cfl/maximum_rate
    end if
    if(.not.ieee_is_finite(candidate).or.candidate<=0) return
    dt=candidate;ok=.true.
  end subroutine reactive_eb_isothermal_timestep_3d

  ! Bounds from open face donors only. Edge/corner and covered-face padding may
  ! be zero or nonfinite and must never enter the physical timestep estimate.
  subroutine reactive_eb_halo_bounds_3d(species,transport,geometry,ghost,temperature, &
      transport_enabled,viscosity,conduction,diffusion,speeds,maximum_diffusivity,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: ghost(:,0:,0:,0:),temperature(0:,0:,0:)
    logical, intent(in) :: transport_enabled,viscosity,conduction,diffusion
    real(dp), intent(out) :: speeds(4),maximum_diffusivity
    logical, intent(out) :: ok
    real(dp) :: q(reactive_nprim(size(species))),recovered,sound,area,dt,dmax,candidate(4),maximum
    real(dp) :: u(reactive_nvar(size(species)),1,1,1),t(1,1,1)
    integer :: n(3),axis,side,first(3),last(3),i,j,k
    logical :: valid
    speeds=0;maximum_diffusivity=0;ok=.false.;candidate=0;maximum=0
    if (.not.geometry%is_valid()) return
    n=[geometry%nx,geometry%ny,geometry%nz]
    if (any(shape(ghost)/=[size(u,1),n+2]).or.any(shape(temperature)/=n+2)) return
    do axis=1,3
    do side=0,1
      first=1;last=n;first(axis)=side*(n(axis)+1);last(axis)=first(axis)
      do k=first(3),last(3)
      do j=first(2),last(2)
      do i=first(1),last(1)
        select case(axis)
        case(1);area=geometry%x_face_fraction(side*n(1),j,k)
        case(2);area=geometry%y_face_fraction(i,side*n(2),k)
        case(3);area=geometry%z_face_fraction(i,j,side*n(3))
        end select
        if (area<=0) cycle
        call reactive_conserved_to_primitive(species,ghost(:,i,j,k),temperature(i,j,k),q,recovered,sound,valid)
        if (.not.valid) return
        candidate(:3)=max(candidate(:3),abs(q(2:4))+sound)
        candidate(4)=max(candidate(4),sum(abs(q(2:4)))+sound)
        if (transport_enabled) then
          u(:,1,1,1)=ghost(:,i,j,k);t(1,1,1)=recovered
          call reactive_transport_timestep_3d(species,transport,u,t,1,1,1,1._dp,1._dp,1._dp,.5_dp, &
            viscosity,conduction,diffusion,dt,dmax,valid)
          if (.not.valid) return
          maximum=max(maximum,dmax)
        end if
      end do
      end do
      end do
    end do
    end do
    speeds=candidate;maximum_diffusivity=maximum;ok=.true.
  end subroutine reactive_eb_halo_bounds_3d
  ! Conservative Euler update without redistribution, for flux-synchronized
  ! patches. The caller supplies a small-cell-limited diffusion timestep.
  subroutine advance_reactive_eb_transport_patch_euler_3d(species,transport,geometry,dt,viscosity,conduction, &
      diffusion,barodiffusion,state,temperature,integrated_flux,minimum_theta,ok,ghost_state,ghost_temperature,faces,physical_faces)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: dt
    logical, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:),minimum_theta
    type(eb_integrated_flux_3d), intent(inout) :: integrated_flux
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: ghost_state(:,0:,0:,0:),ghost_temperature(0:,0:,0:)
    logical, optional, intent(in) :: faces(6)
    character(len=*), optional, intent(in) :: physical_faces(6)
    type(eb_integrated_flux_3d) :: flux
    real(dp), allocatable :: rhs(:,:,:,:),candidate(:,:,:,:),temps(:,:,:),primitive(:)
    real(dp) :: theta,speed
    integer :: i,j,k,nvar
    logical :: valid
    ok=.false.
    if (.not.ieee_is_finite(dt)) return
    if (dt<=0._dp) return
    if (.not.valid_transport_shapes(species,transport,state,temperature,geometry)) return
    nvar=size(state,1)
    allocate(rhs(nvar,geometry%nx,geometry%ny,geometry%nz),primitive(reactive_nprim(size(species))))
    allocate(flux%x(nvar,0:geometry%nx,geometry%ny,geometry%nz), &
      flux%y(nvar,geometry%nx,0:geometry%ny,geometry%nz),flux%z(nvar,geometry%nx,geometry%ny,0:geometry%nz))
    call reactive_eb_transport_fluxes_rhs_3d(species,transport,state,temperature,geometry,dt, &
      viscosity,conduction,diffusion,barodiffusion,rhs,flux%x,flux%y,flux%z,theta,valid, &
      ghost_state,ghost_temperature,faces,physical_faces)
    if (.not.valid) return
    allocate(candidate,source=state);allocate(temps,source=temperature)
    do k=1,geometry%nz
      do j=1,geometry%ny
        do i=1,geometry%nx
          if (geometry%cell_type(i,j,k)==eb_covered_cell_3d) cycle
          candidate(:,i,j,k)=state(:,i,j,k)+dt*rhs(:,i,j,k)
          call reactive_conserved_to_primitive(species,candidate(:,i,j,k),temperature(i,j,k), &
            primitive,temps(i,j,k),speed,valid)
          if (.not.valid) return
        end do
      end do
    end do
    flux%x=dt*flux%x;flux%y=dt*flux%y;flux%z=dt*flux%z
    if (.not.all(ieee_is_finite(flux%x)).or..not.all(ieee_is_finite(flux%y)).or. &
      .not.all(ieee_is_finite(flux%z))) return
    state=candidate;temperature=temps;integrated_flux=flux;minimum_theta=theta;ok=.true.
  end subroutine advance_reactive_eb_transport_patch_euler_3d


  pure logical function valid_transport_shapes( &
      species, transport, state, temperature, geometry) result(valid)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry

    integer :: nvar

    nvar = reactive_nvar(size(species))
    valid = nvar > 0 .and. size(transport) == size(species) .and. &
      geometry%is_valid() .and. size(state, 1) == nvar .and. &
      size(state, 2) == geometry%nx .and. &
      size(state, 3) == geometry%ny .and. &
      size(state, 4) == geometry%nz .and. &
      all(shape(temperature) == &
        [geometry%nx, geometry%ny, geometry%nz])
  end function valid_transport_shapes

  subroutine recover_active_primitives( &
      species, state, temperature, geometry, primitive, &
      checked_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(out) :: primitive(:, :, :, :)
    real(dp), intent(out) :: checked_temperature(:, :, :)
    logical, intent(out) :: ok

    real(dp) :: sound_speed
    logical :: local_ok
    integer :: i, j, k, nprim

    primitive = 0.0_dp
    checked_temperature = 0.0_dp
    ok = .false.
    nprim = reactive_nprim(size(species))
    if (size(primitive, 1) /= nprim .or. &
        size(primitive, 2) /= geometry%nx .or. &
        size(primitive, 3) /= geometry%ny .or. &
        size(primitive, 4) /= geometry%nz .or. &
        any(shape(checked_temperature) /= &
          [geometry%nx, geometry%ny, geometry%nz])) return
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), &
            primitive(:, i, j, k), checked_temperature(i, j, k), &
            sound_speed, local_ok)
          if (.not. local_ok) return
        end do
      end do
    end do
    ok = .true.
  end subroutine recover_active_primitives

  pure real(dp) function active_velocity_derivative(primitive,geometry,i,j,k,velocity_component,direction, &
      periodic) result(derivative)
    real(dp), intent(in) :: primitive(:,:,:,:)
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: i,j,k,velocity_component,direction
    logical, optional, intent(in) :: periodic(3)
    integer :: lower(3),upper(3),idx(3),n(3),component
    real(dp) :: spacing(3),lower_area,upper_area
    logical :: has_lower,has_upper,wrap
    derivative=0
    if(direction<1.or.direction>3) return
    idx=[i,j,k];lower=idx;upper=idx;n=[geometry%nx,geometry%ny,geometry%nz]
    spacing=[geometry%dx,geometry%dy,geometry%dz];component=velocity_component+1
    lower(direction)=idx(direction)-1;upper(direction)=idx(direction)+1
    select case(direction)
    case(1);lower_area=geometry%x_face_fraction(i-1,j,k);upper_area=geometry%x_face_fraction(i,j,k)
    case(2);lower_area=geometry%y_face_fraction(i,j-1,k);upper_area=geometry%y_face_fraction(i,j,k)
    case(3);lower_area=geometry%z_face_fraction(i,j,k-1);upper_area=geometry%z_face_fraction(i,j,k)
    end select
    wrap=.false.
    if(present(periodic)) wrap=periodic(direction)
    has_lower=(idx(direction)>1.or.wrap).and.lower_area>0
    has_upper=(idx(direction)<n(direction).or.wrap).and.upper_area>0
    if(wrap) then
      lower(direction)=modulo(lower(direction)-1,n(direction))+1
      upper(direction)=modulo(upper(direction)-1,n(direction))+1
    end if
    if(has_lower) has_lower=geometry%cell_type(lower(1),lower(2),lower(3))/=eb_covered_cell_3d
    if(has_upper) has_upper=geometry%cell_type(upper(1),upper(2),upper(3))/=eb_covered_cell_3d
    if(has_lower.and.has_upper) then
      derivative=(primitive(component,upper(1),upper(2),upper(3))- &
        primitive(component,lower(1),lower(2),lower(3)))/(2*spacing(direction))
    else if(has_upper) then
      derivative=(primitive(component,upper(1),upper(2),upper(3))-primitive(component,i,j,k))/spacing(direction)
    else if(has_lower) then
      derivative=(primitive(component,i,j,k)-primitive(component,lower(1),lower(2),lower(3)))/spacing(direction)
    end if
  end function active_velocity_derivative

  pure subroutine active_face_velocity_gradient( &
      primitive, geometry, left_i, left_j, left_k, right_i, right_j, &
      right_k, normal, gradient,periodic)
    real(dp), intent(in) :: primitive(:, :, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: left_i, left_j, left_k
    integer, intent(in) :: right_i, right_j, right_k, normal
    real(dp), intent(out) :: gradient(3, 3)
    logical, optional, intent(in) :: periodic(3)

    real(dp) :: spacing(3)
    integer :: velocity_component, direction

    spacing = [geometry%dx, geometry%dy, geometry%dz]
    gradient = 0.0_dp
    do velocity_component = 1, 3
      gradient(velocity_component, normal) = ( &
        primitive(velocity_component + 1, right_i, right_j, right_k) - &
        primitive(velocity_component + 1, left_i, left_j, left_k)) / &
        spacing(normal)
      do direction = 1, 3
        if (direction == normal) cycle
        gradient(velocity_component, direction) = 0.5_dp * ( &
          active_velocity_derivative( &
            primitive, geometry, left_i, left_j, left_k, &
            velocity_component, direction,periodic) + &
          active_velocity_derivative( &
            primitive, geometry, right_i, right_j, right_k, &
            velocity_component, direction,periodic))
      end do
    end do
  end subroutine active_face_velocity_gradient

  subroutine reactive_eb_transport_timestep_3d( &
      species, transport, state, temperature, geometry, transport_cfl, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, dt, maximum_diffusivity, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: transport_cfl
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled
    real(dp), intent(out) :: dt, maximum_diffusivity
    logical, intent(out) :: ok

    real(dp), allocatable :: active_state(:, :, :, :)
    real(dp), allocatable :: active_temperature(:, :, :)
    logical :: found_active
    integer :: i, j, k, reference_i, reference_j, reference_k

    dt = 0.0_dp
    maximum_diffusivity = 0.0_dp
    ok = .false.
    if (.not. valid_transport_shapes( &
          species, transport, state, temperature, geometry)) return
    found_active = .false.
    reference_i = 0
    reference_j = 0
    reference_k = 0
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          reference_i = i
          reference_j = j
          reference_k = k
          found_active = .true.
          exit
        end do
        if (found_active) exit
      end do
      if (found_active) exit
    end do
    if (.not. found_active) return
    allocate(active_state, source=state)
    allocate(active_temperature, source=temperature)
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) /= eb_covered_cell_3d) cycle
          active_state(:, i, j, k) = &
            state(:, reference_i, reference_j, reference_k)
          active_temperature(i, j, k) = &
            temperature(reference_i, reference_j, reference_k)
        end do
      end do
    end do
    call reactive_transport_timestep_3d( &
      species, transport, active_state, active_temperature, &
      geometry%nx, geometry%ny, geometry%nz, geometry%dx, geometry%dy, &
      geometry%dz, transport_cfl, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, dt, &
      maximum_diffusivity, ok)
  end subroutine reactive_eb_transport_timestep_3d

  ! Small-cell estimate for a plain Euler patch. Unlike the redistributed
  ! estimate, open-face conductances are divided by fluid volume. Physical
  ! outer faces are included conservatively even when their transport flux is
  ! disabled. Embedded walls have zero molecular flux and add no conductance.
  subroutine reactive_eb_transport_patch_timestep_3d(species,transport,state,temperature,geometry,transport_cfl, &
      viscosity_enabled,thermal_conduction_enabled,species_diffusion_enabled,dt,maximum_diffusivity,ok, &
      donor_diffusivity)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),transport_cfl
    type(eb_geometry_3d), intent(in) :: geometry
    logical, intent(in) :: viscosity_enabled,thermal_conduction_enabled,species_diffusion_enabled
    real(dp), intent(out) :: dt,maximum_diffusivity
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: donor_diffusivity
    real(dp) :: base_dt,diffusivity,metric,maximum_metric,rate,candidate
    integer :: i,j,k
    logical :: valid
    dt=0._dp;maximum_diffusivity=0._dp;ok=.false.
    if (present(donor_diffusivity)) then
      if (.not.ieee_is_finite(donor_diffusivity)) return
      if (donor_diffusivity<0._dp) return
    end if
    call reactive_eb_transport_timestep_3d(species,transport,state,temperature,geometry,transport_cfl, &
      viscosity_enabled,thermal_conduction_enabled,species_diffusion_enabled,base_dt,diffusivity,valid)
    if (.not.valid) return
    if (.not.(viscosity_enabled.or.thermal_conduction_enabled.or.species_diffusion_enabled)) then
      dt=huge(1._dp);ok=.true.;return
    end if
    if (present(donor_diffusivity)) diffusivity=max(diffusivity,donor_diffusivity)
    if (diffusivity==0._dp) then
      dt=huge(1._dp);ok=.true.;return
    end if
    maximum_metric=0._dp
    do k=1,geometry%nz
      do j=1,geometry%ny
        do i=1,geometry%nx
          if (geometry%cell_type(i,j,k)==eb_covered_cell_3d) cycle
          metric=0.5_dp*(geometry%x_face_fraction(i-1,j,k)+geometry%x_face_fraction(i,j,k))/geometry%dx**2 &
            +0.5_dp*(geometry%y_face_fraction(i,j-1,k)+geometry%y_face_fraction(i,j,k))/geometry%dy**2 &
            +0.5_dp*(geometry%z_face_fraction(i,j,k-1)+geometry%z_face_fraction(i,j,k))/geometry%dz**2
          if (metric>huge(1._dp)*geometry%volume_fraction(i,j,k)) return
          metric=metric/geometry%volume_fraction(i,j,k)
          if (.not.ieee_is_finite(metric).or.metric<=0._dp) return
          maximum_metric=max(maximum_metric,metric)
        end do
      end do
    end do
    if (maximum_metric<=0._dp) return
    if (maximum_metric>1._dp) then
      if (diffusivity>huge(1._dp)/maximum_metric) return
    end if
    rate=diffusivity*maximum_metric
    if (rate<=0._dp) return
    if (rate<1._dp) then
      if (transport_cfl>huge(1._dp)*rate) return
    end if
    candidate=transport_cfl/rate
    if (.not.ieee_is_finite(candidate).or.candidate<=0._dp) return
    dt=candidate;maximum_diffusivity=diffusivity;ok=.true.
  end subroutine reactive_eb_transport_patch_timestep_3d

  pure subroutine limit_species_face( &
      flux, species_energy, nspecies, theta)
    real(dp), intent(inout) :: flux(:)
    real(dp), intent(in) :: species_energy, theta
    integer, intent(in) :: nspecies

    integer :: species_index, component

    do species_index = 1, nspecies
      component = reactive_species_component(species_index)
      flux(component) = theta * flux(component)
    end do
    flux(iet) = flux(iet) + (theta - 1.0_dp) * species_energy
  end subroutine limit_species_face

  subroutine limit_reactive_eb_transport_fluxes_3d( &
      species, state, geometry, dt, x_flux, y_flux, z_flux, &
      species_energy_x, species_energy_y, species_energy_z, &
      minimum_theta, ok,periodic)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), dt
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(inout) :: x_flux(:, 0:, :, :)
    real(dp), intent(inout) :: y_flux(:, :, 0:, :)
    real(dp), intent(inout) :: z_flux(:, :, :, 0:)
    real(dp), intent(in) :: species_energy_x(0:, :, :)
    real(dp), intent(in) :: species_energy_y(:, 0:, :)
    real(dp), intent(in) :: species_energy_z(:, :, 0:)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    logical, intent(in) :: periodic(3)

    real(dp), allocatable :: theta_cell(:, :, :)
    real(dp) :: outgoing, available, candidate, theta_face
    real(dp) :: x_area, y_area, z_area, cell_volume
    integer :: i, j, k, face_i, face_j, face_k
    integer :: species_index, component

    minimum_theta = 1.0_dp
    ok = .false.
    if (.not. ieee_is_finite(dt) .or. dt < 0.0_dp) return
    if (dt <= tiny(1.0_dp)) then
      ok = .true.
      return
    end if
    allocate(theta_cell(geometry%nx, geometry%ny, geometry%nz))
    theta_cell = 1.0_dp
    x_area = geometry%dy * geometry%dz
    y_area = geometry%dx * geometry%dz
    z_area = geometry%dx * geometry%dy
    cell_volume = geometry%dx * geometry%dy * geometry%dz
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          do species_index = 1, size(species)
            component = reactive_species_component(species_index)
            outgoing = x_area * ( &
              max(geometry%x_face_fraction(i, j, k) * &
                x_flux(component, i, j, k), 0.0_dp) + &
              max(-geometry%x_face_fraction(i - 1, j, k) * &
                x_flux(component, i - 1, j, k), 0.0_dp)) + &
              y_area * ( &
              max(geometry%y_face_fraction(i, j, k) * &
                y_flux(component, i, j, k), 0.0_dp) + &
              max(-geometry%y_face_fraction(i, j - 1, k) * &
                y_flux(component, i, j - 1, k), 0.0_dp)) + &
              z_area * ( &
              max(geometry%z_face_fraction(i, j, k) * &
                z_flux(component, i, j, k), 0.0_dp) + &
              max(-geometry%z_face_fraction(i, j, k - 1) * &
                z_flux(component, i, j, k - 1), 0.0_dp))
            available = max(0.0_dp, state(component, i, j, k)) * &
              geometry%volume_fraction(i, j, k) * cell_volume
            if (outgoing > 0.0_dp) then
              candidate = species_safety * available / (dt * outgoing)
              theta_cell(i, j, k) = min(theta_cell(i, j, k), &
                max(0.0_dp, min(1.0_dp, candidate)))
            end if
          end do
        end do
      end do
    end do

    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do face_i = 0, geometry%nx
          theta_face = 1.0_dp
          if (face_i >= 1) then
            if (geometry%cell_type(face_i, j, k) /= eb_covered_cell_3d) &
              theta_face = min(theta_face, theta_cell(face_i, j, k))
          end if
          if (face_i + 1 <= geometry%nx) then
            if (geometry%cell_type(face_i + 1, j, k) /= &
                eb_covered_cell_3d) theta_face = min( &
                  theta_face, theta_cell(face_i + 1, j, k))
          end if
          if(periodic(1).and.(face_i==0.or.face_i==geometry%nx)) &
            theta_face=min(theta_face,theta_cell(1,j,k),theta_cell(geometry%nx,j,k))
          call limit_species_face( &
            x_flux(:, face_i, j, k), species_energy_x(face_i, j, k), &
            size(species), theta_face)
          minimum_theta = min(minimum_theta, theta_face)
        end do
      end do
    end do
    do k = 1, geometry%nz
      do face_j = 0, geometry%ny
        do i = 1, geometry%nx
          theta_face = 1.0_dp
          if (face_j >= 1) then
            if (geometry%cell_type(i, face_j, k) /= eb_covered_cell_3d) &
              theta_face = min(theta_face, theta_cell(i, face_j, k))
          end if
          if (face_j + 1 <= geometry%ny) then
            if (geometry%cell_type(i, face_j + 1, k) /= &
                eb_covered_cell_3d) theta_face = min( &
                  theta_face, theta_cell(i, face_j + 1, k))
          end if
          if(periodic(2).and.(face_j==0.or.face_j==geometry%ny)) &
            theta_face=min(theta_face,theta_cell(i,1,k),theta_cell(i,geometry%ny,k))
          call limit_species_face( &
            y_flux(:, i, face_j, k), species_energy_y(i, face_j, k), &
            size(species), theta_face)
          minimum_theta = min(minimum_theta, theta_face)
        end do
      end do
    end do
    do face_k = 0, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          theta_face = 1.0_dp
          if (face_k >= 1) then
            if (geometry%cell_type(i, j, face_k) /= eb_covered_cell_3d) &
              theta_face = min(theta_face, theta_cell(i, j, face_k))
          end if
          if (face_k + 1 <= geometry%nz) then
            if (geometry%cell_type(i, j, face_k + 1) /= &
                eb_covered_cell_3d) theta_face = min( &
                  theta_face, theta_cell(i, j, face_k + 1))
          end if
          if(periodic(3).and.(face_k==0.or.face_k==geometry%nz)) &
            theta_face=min(theta_face,theta_cell(i,j,1),theta_cell(i,j,geometry%nz))
          call limit_species_face( &
            z_flux(:, i, j, face_k), species_energy_z(i, j, face_k), &
            size(species), theta_face)
          minimum_theta = min(minimum_theta, theta_face)
        end do
      end do
    end do
    ok = all(ieee_is_finite(x_flux)) .and. &
      all(ieee_is_finite(y_flux)) .and. &
      all(ieee_is_finite(z_flux)) .and. &
      ieee_is_finite(minimum_theta) .and. &
      minimum_theta >= 0.0_dp .and. minimum_theta <= 1.0_dp
  end subroutine limit_reactive_eb_transport_fluxes_3d

  subroutine reactive_eb_transport_fluxes_rhs_3d( &
      species, transport, state, temperature, geometry, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, rhs, &
      x_flux, y_flux, z_flux, minimum_theta, ok, ghost_state, ghost_temperature, boundary_faces,physical_faces,wall_temperature)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: dt
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: rhs(:, :, :, :)
    real(dp), intent(out) :: x_flux(:, 0:, :, :)
    real(dp), intent(out) :: y_flux(:, :, 0:, :)
    real(dp), intent(out) :: z_flux(:, :, :, 0:)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: ghost_state(:,0:,0:,0:),ghost_temperature(0:,0:,0:)
    logical, optional, intent(in) :: boundary_faces(6)
    character(len=*), optional, intent(in) :: physical_faces(6)
    real(dp), optional, intent(in) :: wall_temperature(6)
    logical :: enabled_faces(6),periodic(3)
    real(dp), allocatable :: donor_primitive(:)
    real(dp) :: donor_temperature,sound_speed,spacing(3),face_area,species_energy
    integer :: axis,side,start(3),finish(3),inside(3),outside(3),v,d
    real(dp), allocatable :: left(:),right(:),boundary_flux(:)
    real(dp) :: left_t,right_t,wall_heat,wall_rate,other_area,area_ratio
    character(len=32) :: label
    integer :: face,opposite(3)

    real(dp), allocatable :: primitive(:, :, :, :)
    real(dp), allocatable :: checked_temperature(:, :, :)
    real(dp), allocatable :: species_energy_x(:, :, :)
    real(dp), allocatable :: species_energy_y(:, :, :)
    real(dp), allocatable :: species_energy_z(:, :, :)
    real(dp) :: gradient(3, 3), fluid_volume
    real(dp) :: x_area, y_area, z_area
    logical :: local_ok
    integer :: i, j, k, nvar, nprim

    rhs = 0.0_dp
    x_flux = 0.0_dp
    y_flux = 0.0_dp
    z_flux = 0.0_dp
    minimum_theta = 1.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    nprim = reactive_nprim(size(species))
    if (.not. valid_transport_shapes( &
          species, transport, state, temperature, geometry) .or. &
        any(shape(rhs) /= shape(state)) .or. &
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
        size(z_flux, 4) /= geometry%nz + 1 .or. &
        .not. ieee_is_finite(dt) .or. dt < 0.0_dp .or. &
        (barodiffusion_enabled .and. &
          .not. species_diffusion_enabled)) return
    if (present(ghost_state).neqv.present(ghost_temperature)) return
    enabled_faces=present(ghost_state)
    if (present(boundary_faces)) enabled_faces=boundary_faces
    if (any(enabled_faces).and..not.present(ghost_state)) return
    if (present(ghost_state)) then
      if (any(shape(ghost_state)/=[nvar,geometry%nx+2,geometry%ny+2,geometry%nz+2])) return
      if (any(shape(ghost_temperature)/=[geometry%nx+2,geometry%ny+2,geometry%nz+2])) return
    end if
    periodic=.false.
    if(present(physical_faces)) then
      if(.not.valid_reactive_boundaries_3d(physical_faces)) return
      if(.not.periodic_eb_faces_match(geometry,physical_faces)) return
      do axis=1,3
        periodic(axis)=physical_faces(2*axis-1)=='periodic'
        if(periodic(axis).and..not.all(enabled_faces(2*axis-1:2*axis))) return
      end do
      if(any(physical_faces=='isothermal_wall')) then
        if(.not.thermal_conduction_enabled.or..not.present(wall_temperature)) return
      end if
    else
      if(present(wall_temperature)) return
    end if
    if (.not. (viscosity_enabled .or. thermal_conduction_enabled .or. &
        species_diffusion_enabled)) then
      ok = .true.
      return
    end if

    allocate(primitive(nprim, geometry%nx, geometry%ny, geometry%nz))
    allocate(checked_temperature(geometry%nx, geometry%ny, geometry%nz))
    allocate(species_energy_x(0:geometry%nx, geometry%ny, geometry%nz))
    allocate(species_energy_y(geometry%nx, 0:geometry%ny, geometry%nz))
    allocate(species_energy_z(geometry%nx, geometry%ny, 0:geometry%nz))
    species_energy_x = 0.0_dp
    species_energy_y = 0.0_dp
    species_energy_z = 0.0_dp
    call recover_active_primitives( &
      species, state, temperature, geometry, primitive, &
      checked_temperature, local_ok)
    if (.not. local_ok) return

    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx - 1
          if (geometry%x_face_fraction(i, j, k) <= 0.0_dp) cycle
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d .or. &
              geometry%cell_type(i + 1, j, k) == &
                eb_covered_cell_3d) return
          call active_face_velocity_gradient( &
            primitive, geometry, i, j, k, i + 1, j, k, 1, gradient,periodic)
          call reactive_transport_face_flux_3d( &
            species, transport, primitive(:, i, j, k), &
            primitive(:, i + 1, j, k), checked_temperature(i, j, k), &
            checked_temperature(i + 1, j, k), geometry%dx, 1, gradient, &
            viscosity_enabled, thermal_conduction_enabled, &
            species_diffusion_enabled, barodiffusion_enabled, &
            x_flux(:, i, j, k), species_energy_x(i, j, k), local_ok)
          if (.not. local_ok) return
        end do
      end do
    end do
    do k = 1, geometry%nz
      do j = 1, geometry%ny - 1
        do i = 1, geometry%nx
          if (geometry%y_face_fraction(i, j, k) <= 0.0_dp) cycle
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d .or. &
              geometry%cell_type(i, j + 1, k) == &
                eb_covered_cell_3d) return
          call active_face_velocity_gradient( &
            primitive, geometry, i, j, k, i, j + 1, k, 2, gradient,periodic)
          call reactive_transport_face_flux_3d( &
            species, transport, primitive(:, i, j, k), &
            primitive(:, i, j + 1, k), checked_temperature(i, j, k), &
            checked_temperature(i, j + 1, k), geometry%dy, 2, gradient, &
            viscosity_enabled, thermal_conduction_enabled, &
            species_diffusion_enabled, barodiffusion_enabled, &
            y_flux(:, i, j, k), species_energy_y(i, j, k), local_ok)
          if (.not. local_ok) return
        end do
      end do
    end do
    do k = 1, geometry%nz - 1
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%z_face_fraction(i, j, k) <= 0.0_dp) cycle
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d .or. &
              geometry%cell_type(i, j, k + 1) == &
                eb_covered_cell_3d) return
          call active_face_velocity_gradient( &
            primitive, geometry, i, j, k, i, j, k + 1, 3, gradient,periodic)
          call reactive_transport_face_flux_3d( &
            species, transport, primitive(:, i, j, k), &
            primitive(:, i, j, k + 1), checked_temperature(i, j, k), &
            checked_temperature(i, j, k + 1), geometry%dz, 3, gradient, &
            viscosity_enabled, thermal_conduction_enabled, &
            species_diffusion_enabled, barodiffusion_enabled, &
            z_flux(:, i, j, k), species_energy_z(i, j, k), local_ok)
          if (.not. local_ok) return
        end do
      end do
    end do
    ! Parent-provided face donors. Normal gradients use the ghost-center
    ! spacing; tangential velocity gradients use the adjacent interior stencil.
    ! Disabled physical faces retain the existing zero transport flux.
    if (any(enabled_faces)) then
      allocate(donor_primitive(nprim),left(nprim),right(nprim),boundary_flux(nvar))
      spacing=[geometry%dx,geometry%dy,geometry%dz]
      do axis=1,3
        do side=0,1
          if (.not.enabled_faces(2*axis-1+side)) cycle
          if(periodic(axis).and.side==1) cycle
          start=1;finish=[geometry%nx,geometry%ny,geometry%nz]
          start(axis)=side*finish(axis);finish(axis)=start(axis)
          do k=start(3),finish(3)
            do j=start(2),finish(2)
              do i=start(1),finish(1)
                inside=max(1,[i,j,k]);outside=inside;outside(axis)=inside(axis)+2*side-1
                select case(axis)
                case(1)
                  face_area=geometry%x_face_fraction(i,j,k)
                case(2)
                  face_area=geometry%y_face_fraction(i,j,k)
                case(3)
                  face_area=geometry%z_face_fraction(i,j,k)
                end select
                if (face_area<=0._dp) cycle
                if (geometry%cell_type(inside(1),inside(2),inside(3))==eb_covered_cell_3d) return
                if(periodic(axis)) then
                  opposite=inside
                  select case(axis)
                  case(1);opposite(axis)=geometry%nx;other_area=geometry%x_face_fraction(geometry%nx,j,k)
                  case(2);opposite(axis)=geometry%ny;other_area=geometry%y_face_fraction(i,geometry%ny,k)
                  case(3);opposite(axis)=geometry%nz;other_area=geometry%z_face_fraction(i,j,geometry%nz)
                  end select
                  if(other_area<=0) return
                  if(geometry%cell_type(opposite(1),opposite(2),opposite(3))==eb_covered_cell_3d) return
                  donor_primitive=primitive(:,opposite(1),opposite(2),opposite(3))
                  donor_temperature=checked_temperature(opposite(1),opposite(2),opposite(3))
                  area_ratio=face_area/other_area
                else
                  call reactive_conserved_to_primitive(species,ghost_state(:,outside(1),outside(2),outside(3)), &
                    ghost_temperature(outside(1),outside(2),outside(3)),donor_primitive,donor_temperature, &
                    sound_speed,local_ok)
                  if (.not.local_ok) return
                end if
                left=primitive(:,inside(1),inside(2),inside(3));right=left
                left_t=checked_temperature(inside(1),inside(2),inside(3));right_t=left_t
                if (side==0) then
                  left=donor_primitive;left_t=donor_temperature
                else
                  right=donor_primitive;right_t=donor_temperature
                end if
                gradient=0._dp
                do v=1,3
                  gradient(v,axis)=(right(v+1)-left(v+1))/spacing(axis)
                  do d=1,3
                    if (d==axis) cycle
                    gradient(v,d)=active_velocity_derivative(primitive,geometry,inside(1),inside(2),inside(3),v,d,periodic)
                  end do
                end do
                if(periodic(axis)) then
                  call active_face_velocity_gradient(primitive,geometry,opposite(1),opposite(2),opposite(3), &
                    inside(1),inside(2),inside(3),axis,gradient,periodic)
                end if
                label=''
                face=2*axis-1+side
                if(present(physical_faces)) label=physical_faces(face)
                select case(label)
                case('slip_wall')
                  do d=1,3
                    if(d/=axis) gradient(axis,d)=0
                  end do
                case('no_slip_wall','isothermal_wall')
                  do d=1,3
                    if(d/=axis) gradient(:,d)=0
                  end do
                end select
                call reactive_transport_face_flux_3d(species,transport,left,right,left_t,right_t,spacing(axis),axis, &
                  gradient,viscosity_enabled,thermal_conduction_enabled,species_diffusion_enabled, &
                  barodiffusion_enabled,boundary_flux,species_energy,local_ok)
                if (.not.local_ok) return
                select case(label)
                case('slip_wall','no_slip_wall','isothermal_wall')
                  boundary_flux(1)=0;boundary_flux(6:)=0;species_energy=0
                  boundary_flux(iet)=0
                  if(label=='slip_wall') then
                    do v=1,3
                      if(v/=axis) boundary_flux(v+1)=0
                    end do
                  else if(label=='isothermal_wall') then
                    call wall_heat_transfer_3d(species,transport,state(:,inside(1),inside(2),inside(3)), &
                      checked_temperature(inside(1),inside(2),inside(3)),wall_temperature(face), &
                      spacing(axis),wall_heat,wall_rate,local_ok)
                    if(.not.local_ok) return
                    boundary_flux(iet)=real(1-2*side,dp)*wall_heat
                  end if
                end select
                select case(axis)
                case(1)
                  x_flux(:,i,j,k)=boundary_flux;species_energy_x(i,j,k)=species_energy
                  if(periodic(axis)) then
                    x_flux(:,geometry%nx,j,k)=area_ratio*boundary_flux
                    species_energy_x(geometry%nx,j,k)=area_ratio*species_energy
                  end if
                case(2)
                  y_flux(:,i,j,k)=boundary_flux;species_energy_y(i,j,k)=species_energy
                  if(periodic(axis)) then
                    y_flux(:,i,geometry%ny,k)=area_ratio*boundary_flux
                    species_energy_y(i,geometry%ny,k)=area_ratio*species_energy
                  end if
                case(3)
                  z_flux(:,i,j,k)=boundary_flux;species_energy_z(i,j,k)=species_energy
                  if(periodic(axis)) then
                    z_flux(:,i,j,geometry%nz)=area_ratio*boundary_flux
                    species_energy_z(i,j,geometry%nz)=area_ratio*species_energy
                  end if
                end select
              end do
            end do
          end do
        end do
      end do
    end if
    call limit_reactive_eb_transport_fluxes_3d( &
      species, state, geometry, dt, x_flux, y_flux, z_flux, &
      species_energy_x, species_energy_y, species_energy_z, &
      minimum_theta, local_ok,periodic)
    if (.not. local_ok) return

    x_area = geometry%dy * geometry%dz
    y_area = geometry%dx * geometry%dz
    z_area = geometry%dx * geometry%dy
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          fluid_volume = geometry%volume_fraction(i, j, k) * &
            geometry%dx * geometry%dy * geometry%dz
          if (fluid_volume <= 0.0_dp) return
          rhs(:, i, j, k) = -( &
            x_area * ( &
              geometry%x_face_fraction(i, j, k) * &
                x_flux(:, i, j, k) - &
              geometry%x_face_fraction(i - 1, j, k) * &
                x_flux(:, i - 1, j, k)) + &
            y_area * ( &
              geometry%y_face_fraction(i, j, k) * &
                y_flux(:, i, j, k) - &
              geometry%y_face_fraction(i, j - 1, k) * &
                y_flux(:, i, j - 1, k)) + &
            z_area * ( &
              geometry%z_face_fraction(i, j, k) * &
                z_flux(:, i, j, k) - &
              geometry%z_face_fraction(i, j, k - 1) * &
                z_flux(:, i, j, k - 1))) / fluid_volume
        end do
      end do
    end do
    ok = all(ieee_is_finite(rhs))
  end subroutine reactive_eb_transport_fluxes_rhs_3d

  ! Optional caller-prepared donors use one halo layer (indices 0:n+1).
  ! Returned flux is dt times the limited intensive face flux; apply apertures
  ! and face areas once when forming a boundary inventory.
  subroutine reactive_eb_transport_rhs_3d( &
      species, transport, state, temperature, geometry, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, rhs, &
      minimum_theta, ok, ghost_state, ghost_temperature, boundary_faces, integrated_flux,physical_faces,wall_temperature)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: dt
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: rhs(:, :, :, :)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: ghost_state(:,0:,0:,0:),ghost_temperature(0:,0:,0:)
    logical, optional, intent(in) :: boundary_faces(6)
    character(len=*), optional, intent(in) :: physical_faces(6)
    real(dp), optional, intent(in) :: wall_temperature(6)
    type(eb_integrated_flux_3d), optional, intent(inout) :: integrated_flux

    real(dp), allocatable :: x_flux(:, :, :, :)
    real(dp), allocatable :: y_flux(:, :, :, :)
    real(dp), allocatable :: z_flux(:, :, :, :)

    allocate(x_flux(size(state, 1), 0:geometry%nx, &
      geometry%ny, geometry%nz))
    allocate(y_flux(size(state, 1), geometry%nx, &
      0:geometry%ny, geometry%nz))
    allocate(z_flux(size(state, 1), geometry%nx, &
      geometry%ny, 0:geometry%nz))
    call reactive_eb_transport_fluxes_rhs_3d( &
      species, transport, state, temperature, geometry, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, rhs, &
      x_flux, y_flux, z_flux, minimum_theta, ok, ghost_state, ghost_temperature, boundary_faces, &
      physical_faces,wall_temperature)
    if (.not.ok) return
    if (present(integrated_flux)) then
      x_flux=dt*x_flux;y_flux=dt*y_flux;z_flux=dt*z_flux
      ok=all(ieee_is_finite(x_flux)).and.all(ieee_is_finite(y_flux)).and.all(ieee_is_finite(z_flux))
      if (.not.ok) return
      call move_alloc(x_flux,integrated_flux%x)
      call move_alloc(y_flux,integrated_flux%y)
      call move_alloc(z_flux,integrated_flux%z)
    end if
  end subroutine reactive_eb_transport_rhs_3d

  ! The caller must rebuild donors for each time stage. Flux publication waits
  ! for successful redistribution and thermodynamic recovery; failure restores
  ! the input state/temperature and leaves the caller's flux object untouched.
  subroutine reactive_eb_transport_euler_update_3d( &
      species, transport, state, temperature, geometry, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, &
      target_volume_fraction, new_state, new_temperature, &
      minimum_theta, ok, ghost_state, ghost_temperature, boundary_faces, integrated_flux,physical_faces,wall_temperature)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: dt, target_volume_fraction
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: new_state(:, :, :, :)
    real(dp), intent(out) :: new_temperature(:, :, :)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: ghost_state(:,0:,0:,0:),ghost_temperature(0:,0:,0:)
    logical, optional, intent(in) :: boundary_faces(6)
    character(len=*), optional, intent(in) :: physical_faces(6)
    real(dp), optional, intent(in) :: wall_temperature(6)
    type(eb_integrated_flux_3d), optional, intent(inout) :: integrated_flux

    real(dp), allocatable :: rhs(:, :, :, :),candidate(:,:,:,:),candidate_temperature(:,:,:)
    type(eb_integrated_flux_3d) :: flux
    logical :: local_ok

    minimum_theta = 1.0_dp
    ok = .false.
    if (.not. valid_transport_shapes( &
          species, transport, state, temperature, geometry) .or. &
        any(shape(new_state) /= shape(state)) .or. &
        any(shape(new_temperature) /= shape(temperature)) .or. &
        .not. ieee_is_finite(dt) .or. dt < 0.0_dp .or. &
        (barodiffusion_enabled .and. &
          .not. species_diffusion_enabled)) return
    new_state = state
    new_temperature = temperature
    if (dt <= tiny(1.0_dp) .or. .not. (viscosity_enabled .or. &
        thermal_conduction_enabled .or. species_diffusion_enabled)) then
      if (present(integrated_flux)) then
        allocate(flux%x(size(state,1),0:geometry%nx,geometry%ny,geometry%nz), &
          flux%y(size(state,1),geometry%nx,0:geometry%ny,geometry%nz), &
          flux%z(size(state,1),geometry%nx,geometry%ny,0:geometry%nz))
        flux%x=0;flux%y=0;flux%z=0;integrated_flux=flux
      end if
      ok = .true.
      return
    end if
    allocate(rhs, mold=state)
    call reactive_eb_transport_rhs_3d( &
      species, transport, state, temperature, geometry, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, rhs, &
      minimum_theta, local_ok, ghost_state, ghost_temperature, boundary_faces, flux,physical_faces,wall_temperature)
    if (.not. local_ok) return
    allocate(candidate,mold=state);allocate(candidate_temperature,mold=temperature)
    call advance_reactive_eb_state_redistributed_3d( &
      species, state, temperature, geometry, rhs, dt, candidate, &
      candidate_temperature, local_ok, target_volume_fraction)
    if (.not. local_ok) return
    new_state=candidate;new_temperature=candidate_temperature
    if (present(integrated_flux)) integrated_flux=flux
    ok = .true.
  end subroutine reactive_eb_transport_euler_update_3d

  subroutine advance_reactive_eb_transport_3d( &
      species, transport, state, temperature, geometry, interval, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, &
      target_volume_fraction, minimum_theta, ok, prepare_halo, integrated_flux, &
      boundary_faces,inflow_state,inflow_temperature,wall_temperature)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: interval, target_volume_fraction
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    procedure(eb_transport_halo_provider_3d), optional :: prepare_halo
    type(eb_integrated_flux_3d), optional, intent(inout) :: integrated_flux
    type(eb_integrated_flux_3d), allocatable :: flux_one,flux_two
    real(dp), allocatable :: ghost(:,:,:,:),ghost_temperature(:,:,:)
    logical :: faces(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_temperature(6)

    real(dp), allocatable :: initial_state(:, :, :, :)
    real(dp), allocatable :: initial_temperature(:, :, :)
    real(dp), allocatable :: stage_state(:, :, :, :)
    real(dp), allocatable :: stage_temperature(:, :, :)
    real(dp), allocatable :: euler_state(:, :, :, :)
    real(dp), allocatable :: euler_temperature(:, :, :)
    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    real(dp), allocatable :: primitive(:)
    real(dp) :: theta_one, theta_two, recovered_temperature, sound_speed
    logical :: local_ok
    integer :: i, j, k

    minimum_theta = 1.0_dp
    ok = .false.
    if (.not. valid_transport_shapes( &
          species, transport, state, temperature, geometry) .or. &
        .not. ieee_is_finite(interval) .or. interval < 0.0_dp .or. &
        .not. ieee_is_finite(target_volume_fraction) .or. &
        target_volume_fraction <= 0.0_dp .or. &
        target_volume_fraction > 1.0_dp .or. &
        (barodiffusion_enabled .and. &
          .not. species_diffusion_enabled)) return
    if(present(boundary_faces)) then
      if(present(prepare_halo)) return
      if(.not.valid_reactive_boundaries_3d(boundary_faces)) return
      if(any(boundary_faces=='isothermal_wall')) then
        if(.not.thermal_conduction_enabled.or..not.present(wall_temperature)) return
      end if
    else
      if(present(inflow_state).or.present(inflow_temperature).or.present(wall_temperature)) return
    end if
    if (present(integrated_flux)) then
      allocate(flux_one,flux_two)
      allocate(flux_one%x(size(state,1),0:geometry%nx,geometry%ny,geometry%nz), &
        flux_one%y(size(state,1),geometry%nx,0:geometry%ny,geometry%nz), &
        flux_one%z(size(state,1),geometry%nx,geometry%ny,0:geometry%nz))
      flux_one%x=0;flux_one%y=0;flux_one%z=0
    end if
    if (interval <= tiny(1.0_dp) .or. .not. (viscosity_enabled .or. &
        thermal_conduction_enabled .or. species_diffusion_enabled)) then
      if (present(integrated_flux)) integrated_flux=flux_one
      ok = .true.
      return
    end if
    allocate(initial_state, source=state)
    allocate(initial_temperature, source=temperature)
    allocate(stage_state, mold=state)
    allocate(stage_temperature, mold=temperature)
    allocate(euler_state, mold=state)
    allocate(euler_temperature, mold=temperature)
    faces=.false.
    if (present(prepare_halo).or.present(boundary_faces)) then
      allocate(ghost(size(state,1),0:geometry%nx+1,0:geometry%ny+1,0:geometry%nz+1), &
        ghost_temperature(0:geometry%nx+1,0:geometry%ny+1,0:geometry%nz+1))
      if(present(prepare_halo)) then
        call prepare_halo(initial_state,initial_temperature,geometry,ghost,ghost_temperature,faces,local_ok)
      else
        faces=.true.
        call fill_reactive_boundary_ghosts_3d(initial_state,initial_temperature,boundary_faces, &
          ghost,ghost_temperature,local_ok,viscous_wall=.true., &
          inflow_state=inflow_state,inflow_temperature=inflow_temperature,ghost_width=1)
      end if
      if (.not.local_ok) return
    end if
    call reactive_eb_transport_euler_update_3d( &
      species, transport, initial_state, initial_temperature, geometry, &
      interval, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, &
      target_volume_fraction, stage_state, stage_temperature, theta_one, &
      local_ok,ghost,ghost_temperature,faces,flux_one,boundary_faces,wall_temperature)
    if (.not. local_ok) return
    if (present(prepare_halo).or.present(boundary_faces)) then
      if(present(prepare_halo)) then
        call prepare_halo(stage_state,stage_temperature,geometry,ghost,ghost_temperature,faces,local_ok)
      else
        faces=.true.
        call fill_reactive_boundary_ghosts_3d(stage_state,stage_temperature,boundary_faces, &
          ghost,ghost_temperature,local_ok,viscous_wall=.true., &
          inflow_state=inflow_state,inflow_temperature=inflow_temperature,ghost_width=1)
      end if
      if (.not.local_ok) return
    end if
    call reactive_eb_transport_euler_update_3d( &
      species, transport, stage_state, stage_temperature, geometry, &
      interval, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, &
      target_volume_fraction, euler_state, euler_temperature, theta_two, &
      local_ok,ghost,ghost_temperature,faces,flux_two,boundary_faces,wall_temperature)
    if (.not. local_ok) return

    allocate(candidate_state, source=0.5_dp * (initial_state + euler_state))
    allocate(candidate_temperature, source=initial_temperature)
    allocate(primitive(reactive_nprim(size(species))))
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) then
            candidate_state(:, i, j, k) = initial_state(:, i, j, k)
            candidate_temperature(i, j, k) = initial_temperature(i, j, k)
            cycle
          end if
          call reactive_conserved_to_primitive( &
            species, candidate_state(:, i, j, k), &
            0.5_dp * (initial_temperature(i, j, k) + &
              euler_temperature(i, j, k)), primitive, &
            recovered_temperature, sound_speed, local_ok)
          if (.not. local_ok) return
          candidate_temperature(i, j, k) = recovered_temperature
        end do
      end do
    end do
    if (any(.not. ieee_is_finite(candidate_state)) .or. &
        any(.not. ieee_is_finite(candidate_temperature))) return
    if (present(integrated_flux)) then
      flux_one%x=.5_dp*flux_one%x+.5_dp*flux_two%x
      flux_one%y=.5_dp*flux_one%y+.5_dp*flux_two%y
      flux_one%z=.5_dp*flux_one%z+.5_dp*flux_two%z
      if (.not.all(ieee_is_finite(flux_one%x)).or..not.all(ieee_is_finite(flux_one%y)).or. &
        .not.all(ieee_is_finite(flux_one%z))) return
      integrated_flux=flux_one
    end if
    state = candidate_state
    temperature = candidate_temperature
    minimum_theta = min(theta_one, theta_two)
    ok = .true.
  end subroutine advance_reactive_eb_transport_3d

end module eb_reactive_transport_3d_mod
