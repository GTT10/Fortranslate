module eb_reactive_hydro_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use reactive_3d_mod, only: adjust_reactive_artificial_viscosity_3d
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive, &
    reactive_riemann_flux_x,reconstruct_pelec_plm_cell,compute_pelec_plm_slopes
  use reactive_2d_mod, only: reactive_riemann_flux_y
  use reactive_directional_flux_3d_mod, only: reactive_riemann_flux_z
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_reactive_wall_flux_3d_mod, only: &
    reactive_eb_flux_divergence_3d,reactive_eb_slip_wall_source_3d
  use eb_reactive_redistribution_3d_mod, only: &
    advance_reactive_eb_redistributed_3d, &
    advance_reactive_eb_state_redistributed_3d
  implicit none
  private

  public :: predict_reactive_eb_plm_cell_3d
  public :: reactive_eb_outflow_riemann_fluxes_3d,reactive_eb_riemann_fluxes_3d
  public :: advance_reactive_eb_euler_3d
  public :: advance_reactive_eb_redistributed_euler_3d
  public :: advance_reactive_eb_state_redistributed_euler_3d

  interface reactive_eb_riemann_fluxes_3d
    module procedure reactive_eb_outflow_riemann_fluxes_3d
  end interface

contains

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
