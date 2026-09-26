program test_reactive_eb_hydro_3d
  use precision_mod, only: dp
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config,build_configured_eb_geometry_3d
  use state_indices_mod, only: irho
  use reactive_3d_mod, only: adjust_reactive_artificial_viscosity_3d
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_species_component, reactive_primitive_to_conserved,reactive_pelec_primitive_flux_x
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d, build_axis_plane_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_wall_flux_3d_mod, only: &
    reactive_eb_flux_divergence_3d
  use eb_reactive_hydro_3d_mod, only: &
    reactive_eb_outflow_riemann_fluxes_3d,predict_reactive_eb_plm_cell_3d, &
    compute_reactive_eb_plm_normal_fluxes_3d,compute_reactive_eb_plm_ctu_fluxes_3d,advance_reactive_eb_plm_ctu_3d, &
    eb_plm_face_prediction_3d,eb_plm_boundary_3d,apply_reactive_eb_transverse_3d,compute_reactive_eb_flattening_3d, &
    advance_reactive_eb_euler_3d,advance_reactive_eb_outflow_ctu_3d,advance_reactive_eb_boundary_ctu_3d, &
    advance_reactive_eb_redistributed_euler_3d, &
    advance_reactive_eb_state_redistributed_euler_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  implicit none

  integer, parameter :: nx = 10
  integer, parameter :: ny = 8
  integer, parameter :: nz = 6
  real(dp), parameter :: pressure = 135000.0_dp
  type(nasa7_species), allocatable :: species(:)
  type(eb_geometry_3d) :: geometry
  real(dp), allocatable :: primitive(:), state_cell(:), perturbed_cell(:)
  real(dp), allocatable :: mass_fractions(:)
  real(dp), allocatable :: state(:, :, :, :), new_state(:, :, :, :)
  real(dp), allocatable :: temperature(:, :, :), new_temperature(:, :, :)
  real(dp), allocatable :: x_flux(:, :, :, :), y_flux(:, :, :, :)
  real(dp), allocatable :: z_flux(:, :, :, :)
  real(dp), allocatable :: before_integrals(:), after_integrals(:)
  real(dp) :: mole_fractions(7), base_temperature, perturbed_temperature
  real(dp) :: sound_speed, perturbed_sound_speed
  real(dp) :: maximum_uniform_state_error
  real(dp) :: maximum_uniform_temperature_error
  real(dp) :: conservation_error, species_closure_error
  real(dp) :: small_cell_full_step, redistributed_cut_density
  real(dp) :: state_redistributed_cut_density
  logical :: ok
  integer :: i, j, k, nvar, first_species, last_species

  maximum_uniform_state_error = 0.0_dp
  maximum_uniform_temperature_error = 0.0_dp

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamic database load")
  nvar = reactive_nvar(size(species))
  first_species = reactive_species_component(1)
  last_species = reactive_species_component(size(species))
  allocate(primitive(reactive_nprim(size(species))))
  allocate(state_cell(nvar), perturbed_cell(nvar), mass_fractions(size(species)))
  allocate(state(nvar, nx, ny, nz), new_state(nvar, nx, ny, nz))
  allocate(temperature(nx, ny, nz), new_temperature(nx, ny, nz))
  allocate(x_flux(nvar, 0:nx, ny, nz))
  allocate(y_flux(nvar, nx, 0:ny, nz))
  allocate(z_flux(nvar, nx, ny, 0:nz))
  allocate(before_integrals(nvar), after_integrals(nvar))

  mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.55643_dp]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  call require(ok, "composition conversion")

  call verify_uniform_tangent_update("x", "hllc", &
    [0.0_dp, 18.0_dp, -11.0_dp])
  call verify_uniform_tangent_update("y", "rusanov", &
    [18.0_dp, 0.0_dp, -11.0_dp])
  call verify_uniform_tangent_update("z", "pelec", &
    [18.0_dp, -11.0_dp, 0.0_dp])
  call verify_small_cell_hydro_relief()
  call verify_viscosity()
  call verify_normal_prediction()
  call verify_normal_fluxes()
  call verify_transverse_prediction()
  call verify_ctu_fluxes()
  call verify_ctu_update()
  call verify_flattening()
  call verify_physical_fluxes()
  call verify_interior_wall_viscosity()
  call verify_outflow_ctu()
  call verify_boundary_ctu()

  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", 1.1_dp, geometry, ok)
  call require(ok, "covered solver-validation geometry")
  call reactive_eb_outflow_riemann_fluxes_3d( &
    species, state, temperature, geometry, "unknown", &
    x_flux, y_flux, z_flux, ok)
  call require(.not. ok .and. maxval(abs(x_flux)) == 0.0_dp .and. &
    maxval(abs(y_flux)) == 0.0_dp .and. maxval(abs(z_flux)) == 0.0_dp, &
    "covered geometry invalid solver rejection")

  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", 0.37_dp, geometry, ok)
  call require(ok, "perturbation geometry")
  call build_cell_state(0.31_dp, [0.0_dp, 0.0_dp, 0.0_dp], &
    state_cell, base_temperature, sound_speed)
  call fill_state(state_cell, base_temperature)
  call build_cell_state(0.34_dp, [0.0_dp, 0.0_dp, 0.0_dp], &
    perturbed_cell, perturbed_temperature, perturbed_sound_speed)
  state(:, 7, 4, 3) = perturbed_cell
  temperature(7, 4, 3) = perturbed_temperature

  call reactive_eb_outflow_riemann_fluxes_3d( &
    species, state, temperature, geometry, "rusanov", &
    x_flux, y_flux, z_flux, ok)
  call require(ok, "perturbation face fluxes")
  call require(maxval(abs(x_flux(:, 3, :, :))) == 0.0_dp, &
    "closed embedded-normal face flux")
  call require(maxval(abs(x_flux(:, 4, :, :))) > 0.0_dp, &
    "open embedded-normal face flux")

  call fluid_integrals(state, geometry, before_integrals)
  call advance_reactive_eb_euler_3d( &
    species, state, temperature, geometry, "rusanov", 1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(ok, "perturbation Euler update")
  call fluid_integrals(new_state, geometry, after_integrals)
  conservation_error = maxval(abs(after_integrals - before_integrals) / &
    max(1.0_dp, abs(before_integrals)))
  call require(conservation_error <= 3.0e-13_dp, &
    "fluid-volume conservation")
  species_closure_error = maxval(abs( &
      sum(new_state(first_species:last_species, :, :, :), dim=1) - &
        new_state(irho, :, :, :)), &
      mask=geometry%cell_type /= eb_covered_cell_3d)
  call require(species_closure_error <= 3.0e-13_dp, &
    "species-density closure")
  call require(minval(new_temperature, &
    mask=geometry%cell_type /= eb_covered_cell_3d) > 0.0_dp, &
    "positive recovered temperature")
  call require(all(new_state(:, 1:3, :, :) == state(:, 1:3, :, :)) .and. &
    all(new_temperature(1:3, :, :) == temperature(1:3, :, :)), &
    "covered cells unchanged")

  call advance_reactive_eb_redistributed_euler_3d( &
    species, state, temperature, geometry, "rusanov", 1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(ok, "perturbation redistributed Euler update")
  call fluid_integrals(new_state, geometry, after_integrals)
  conservation_error = max(conservation_error, &
    maxval(abs(after_integrals - before_integrals) / &
      max(1.0_dp, abs(before_integrals))))
  call require(conservation_error <= 3.0e-13_dp, &
    "redistributed fluid-volume conservation")
  species_closure_error = max(species_closure_error, maxval(abs( &
      sum(new_state(first_species:last_species, :, :, :), dim=1) - &
        new_state(irho, :, :, :)), &
      mask=geometry%cell_type /= eb_covered_cell_3d))
  call require(species_closure_error <= 3.0e-13_dp, &
    "redistributed species-density closure")
  call require(minval(new_temperature, &
    mask=geometry%cell_type /= eb_covered_cell_3d) > 0.0_dp, &
    "positive redistributed temperature")
  call require(all(new_state(:, 1:3, :, :) == state(:, 1:3, :, :)) .and. &
    all(new_temperature(1:3, :, :) == temperature(1:3, :, :)), &
    "redistributed covered cells unchanged")

  call advance_reactive_eb_state_redistributed_euler_3d( &
    species, state, temperature, geometry, "rusanov", 1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(ok, "perturbation StateRedist Euler update")
  call fluid_integrals(new_state, geometry, after_integrals)
  conservation_error = max(conservation_error, &
    maxval(abs(after_integrals - before_integrals) / &
      max(1.0_dp, abs(before_integrals))))
  call require(conservation_error <= 3.0e-13_dp, &
    "StateRedist fluid-volume conservation")
  species_closure_error = max(species_closure_error, maxval(abs( &
      sum(new_state(first_species:last_species, :, :, :), dim=1) - &
        new_state(irho, :, :, :)), &
      mask=geometry%cell_type /= eb_covered_cell_3d))
  call require(species_closure_error <= 3.0e-13_dp, &
    "StateRedist species-density closure")
  call require(minval(new_temperature, &
    mask=geometry%cell_type /= eb_covered_cell_3d) > 0.0_dp, &
    "positive StateRedist temperature")
  call require(all(new_state(:, 1:3, :, :) == state(:, 1:3, :, :)) .and. &
    all(new_temperature(1:3, :, :) == temperature(1:3, :, :)), &
    "StateRedist covered cells unchanged")

  call advance_reactive_eb_euler_3d( &
    species, state, temperature, geometry, "unknown", 1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(.not. ok .and. all(new_state == state) .and. &
    all(new_temperature == temperature), "invalid solver rollback")
  call advance_reactive_eb_euler_3d( &
    species, state, temperature, geometry, "rusanov", -1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(.not. ok .and. all(new_state == state) .and. &
    all(new_temperature == temperature), "invalid timestep rollback")
  call advance_reactive_eb_redistributed_euler_3d( &
    species, state, temperature, geometry, "unknown", 1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(.not. ok .and. all(new_state == state) .and. &
    all(new_temperature == temperature), &
    "redistributed invalid solver rollback")
  call advance_reactive_eb_redistributed_euler_3d( &
    species, state, temperature, geometry, "rusanov", -1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(.not. ok .and. all(new_state == state) .and. &
    all(new_temperature == temperature), &
    "redistributed invalid timestep rollback")
  call advance_reactive_eb_state_redistributed_euler_3d( &
    species, state, temperature, geometry, "unknown", 1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(.not. ok .and. all(new_state == state) .and. &
    all(new_temperature == temperature), &
    "StateRedist invalid solver rollback")
  call advance_reactive_eb_state_redistributed_euler_3d( &
    species, state, temperature, geometry, "rusanov", -1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(.not. ok .and. all(new_state == state) .and. &
    all(new_temperature == temperature), &
    "StateRedist invalid timestep rollback")
  call advance_reactive_eb_state_redistributed_euler_3d( &
    species, state, temperature, geometry, "rusanov", 1.0e-8_dp, &
    new_state, new_temperature, ok, 0.0_dp)
  call require(.not. ok .and. all(new_state == state) .and. &
    all(new_temperature == temperature), &
    "StateRedist invalid target rollback")

  state(irho, 4, 2, 3) = -1.0_dp
  call advance_reactive_eb_euler_3d( &
    species, state, temperature, geometry, "rusanov", 1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(.not. ok .and. all(new_state == state) .and. &
    all(new_temperature == temperature), "invalid cut-state rollback")
  call advance_reactive_eb_redistributed_euler_3d( &
    species, state, temperature, geometry, "rusanov", 1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(.not. ok .and. all(new_state == state) .and. &
    all(new_temperature == temperature), &
    "redistributed invalid cut-state rollback")
  call advance_reactive_eb_state_redistributed_euler_3d( &
    species, state, temperature, geometry, "rusanov", 1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(.not. ok .and. all(new_state == state) .and. &
    all(new_temperature == temperature), &
    "StateRedist invalid cut-state rollback")

  write(*, '(a,es24.16)') &
    "maximum uniform state error: ", maximum_uniform_state_error
  write(*, '(a,es24.16)') &
    "maximum uniform temperature error: ", &
    maximum_uniform_temperature_error
  write(*, '(a,es24.16)') &
    "maximum conservation error: ", conservation_error
  write(*, '(a,es24.16)') &
    "maximum species closure error: ", species_closure_error
  write(*, '(a,es24.16)') &
    "small-cell full-grid timestep: ", small_cell_full_step
  write(*, '(a,es24.16)') &
    "FluxRedist cut-cell density: ", redistributed_cut_density
  write(*, '(a,es24.16)') &
    "StateRedist cut-cell density: ", state_redistributed_cut_density
  write(*, '(a)') "test_reactive_eb_hydro_3d: PASS"

contains
  subroutine verify_outflow_ctu()
    type(reactive_eb_3d_config) :: cfg
    type(eb_geometry_3d) :: wide,physical
    real(dp) :: u(nvar,4,4,4),t(4,4,4),updated(nvar,4,4,4),ut(4,4,4),boundary(nvar),wall(nvar),cs,vel(3),change,scale
    integer :: axis,mode,ii,jj,kk,c,tangent
    character(len=24) :: faces(6)
    character(len=1), parameter :: axes(3)=['x','y','z']
    character(len=12), parameter :: modes(3)=[character(len=12)::'none','flux_redist','state_redist']
    logical :: valid
    cfg=reactive_eb_3d_config();cfg%nx=4;cfg%ny=4;cfg%nz=4;cfg%plane_position=.395_dp
    do axis=1,3
      cfg%plane_axis=axes(axis)
      call build_configured_eb_geometry_3d(cfg,wide,valid,halo=6)
      call require(valid,'outflow extended geometry')
      call extract_eb_geometry_3d(wide,[7,7,7],[10,10,10],physical,valid)
      call require(valid,'outflow physical geometry')
      u=0;t=0;vel=[1._dp,2._dp,3._dp];vel(axis)=0
      do kk=1,4
      do jj=1,4
      do ii=1,4
        if(physical%volume_fraction(ii,jj,kk)<=0) cycle
        call build_cell_state(.31_dp+.002_dp*(ii+jj+kk),vel,u(:,ii,jj,kk),t(ii,jj,kk),cs)
      end do
      end do
      end do
      do mode=1,3
        call advance_reactive_eb_outflow_ctu_3d(species,u,t,wide,1e-7_dp,updated,ut,valid, &
          redistribution=trim(modes(mode)),boundary_change=boundary,wall_change=wall)
        call require(valid.and.any(updated/=u),'outflow prepared CTU update')
        do c=1,nvar
          change=sum((updated(c,:,:,:)-u(c,:,:,:))*physical%volume_fraction)/4._dp**3
          scale=max(1._dp,sum(abs(u(c,:,:,:))*physical%volume_fraction)/4._dp**3)
          call require(abs(change-boundary(c)-wall(c))<2e-12_dp*scale,'outflow CTU inventory closure')
        end do
      end do
      faces='outflow';tangent=modulo(axis,3)+1
      faces(2*tangent-1:2*tangent)='slip_wall'
      call advance_reactive_eb_boundary_ctu_3d(species,u,t,wide,1e-7_dp,updated,ut,valid,faces, &
        boundary_change=boundary,wall_change=wall)
      call require(valid,'cut geometry with tangential physical walls')
      do c=1,nvar
        change=sum((updated(c,:,:,:)-u(c,:,:,:))*physical%volume_fraction)/64
        scale=max(1._dp,sum(abs(u(c,:,:,:))*physical%volume_fraction)/64)
        call require(abs(change-boundary(c)-wall(c))<2e-12_dp*scale,'cut physical-wall CTU conservation')
      end do
      updated=-17;ut=-18;boundary=-19;wall=-20
      call advance_reactive_eb_outflow_ctu_3d(species,u,t,wide,-1._dp,updated,ut,valid, &
        boundary_change=boundary,wall_change=wall)
      call require(.not.valid.and.all(updated==-17).and.all(ut==-18).and.all(boundary==-19).and.all(wall==-20), &
        'outflow CTU transactional failure')
    end do
  end subroutine verify_outflow_ctu

  subroutine verify_boundary_ctu()
    type(reactive_eb_3d_config) :: cfg
    type(eb_geometry_3d) :: wide
    type(eb_integrated_flux_3d) :: flux
    real(dp) :: u(nvar,4,4,4),t(4,4,4),updated(nvar,4,4,4),ut(4,4,4),boundary(nvar),wall(nvar)
    real(dp) :: inlet(nvar,6),tin(6),cs,change,scale
    character(len=24) :: faces(6)
    character(len=24), parameter :: kinds(5)=[character(len=24):: &
      'periodic','slip_wall','no_slip_wall','isothermal_wall','inflow']
    integer :: axis,mode,i,j,k,c
    logical :: valid
    cfg=reactive_eb_3d_config();cfg%nx=4;cfg%ny=4;cfg%nz=4;cfg%plane_position=-1
    call build_configured_eb_geometry_3d(cfg,wide,valid,halo=6)
    call require(valid.and.all(wide%volume_fraction==1),'physical boundary regular geometry')
    do k=1,4
    do j=1,4
    do i=1,4
      call build_cell_state(.31_dp+.002_dp*(i+j+k),[1._dp,2._dp,3._dp],u(:,i,j,k),t(i,j,k),cs)
    end do
    end do
    end do
    do c=1,6
      inlet(:,c)=u(:,2,2,2);tin(c)=t(2,2,2)
    end do
    do axis=1,3
    do mode=1,5
      faces='outflow';faces(2*axis-1:2*axis)=kinds(mode)
      call advance_reactive_eb_boundary_ctu_3d(species,u,t,wide,1e-7_dp,updated,ut,valid,faces, &
        integrated_flux=flux,boundary_change=boundary,wall_change=wall,inflow_state=inlet,inflow_temperature=tin)
      call require(valid,'physical boundary CTU step')
      do c=1,nvar
        change=sum(updated(c,:,:,:)-u(c,:,:,:))/64
        scale=max(1._dp,sum(abs(u(c,:,:,:)))/64)
        call require(abs(change-boundary(c)-wall(c))<2e-12_dp*scale,'physical boundary inventory closure')
      end do
      if(mode>=2.and.mode<=4) then
        select case(axis)
        case(1)
          call require(all(flux%x(1,0,:,:)==0).and.all(flux%x(1,4,:,:)==0),'x wall mass closure')
        case(2)
          call require(all(flux%y(1,:,0,:)==0).and.all(flux%y(1,:,4,:)==0),'y wall mass closure')
        case(3)
          call require(all(flux%z(1,:,:,0)==0).and.all(flux%z(1,:,:,4)==0),'z wall mass closure')
        end select
      end if
    end do
    end do
    updated=-17;ut=-18;boundary=-19;wall=-20
    faces='outflow';faces(1)='inflow'
    call advance_reactive_eb_boundary_ctu_3d(species,u,t,wide,1e-7_dp,updated,ut,valid,faces, &
      boundary_change=boundary,wall_change=wall)
    call require(.not.valid.and.all(updated==-17).and.all(ut==-18).and.all(boundary==-19).and.all(wall==-20), &
      'missing inflow transactional rejection')
    faces(1)='periodic'
    call advance_reactive_eb_boundary_ctu_3d(species,u,t,wide,1e-7_dp,updated,ut,valid,faces)
    call require(.not.valid.and.all(updated==-17).and.all(ut==-18),'unpaired periodic rejection')
  end subroutine verify_boundary_ctu

  subroutine verify_interior_wall_viscosity()
    real(dp) :: u(nvar,0:5,0:5,0:5),fx(nvar,0:4,4,4),fy(nvar,4,0:4,4),fz(nvar,4,4,0:4)
    integer :: ii,jj,kk,face
    logical :: valid
    u=0;u(1,:,:,:)=1;u(5,:,:,:)=1000;u(6,:,:,:)=1
    do kk=0,5
    do jj=0,5
    do ii=0,5
      u(2:4,ii,jj,kk)=-real(ii+jj+kk,dp)*[1._dp,2._dp,3._dp]
    end do
    end do
    end do
    fx=0;fy=0;fz=0
    call adjust_reactive_artificial_viscosity_3d(u,[4,4,4],[1._dp,1._dp,1._dp],.1_dp, &
      [.true.,.true.,.true.,.true.,.true.,.true.],fx,fy,fz,valid,boundary_lower=[1,1,1],boundary_upper=[3,3,3])
    call require(valid,'interior physical wall viscosity')
    do face=1,3,2
      call require(all(fx(3:4,face,:,:)==0).and.all(fy(2,:,face,:)==0).and.all(fy(4,:,face,:)==0).and. &
        all(fz(2:3,:,:,face)==0),'six interior no-slip faces exclude tangential viscosity')
    end do
    call require(all(fx(3:4,2,:,:)/=0).and.all(fy(2,:,2,:)/=0).and.all(fz(2:3,:,:,2)/=0), &
      'ordinary interior faces retain tangential viscosity')
    call require(all(fx(3:4,0,:,:)/=0).and.all(fx(3:4,4,:,:)/=0), &
      'working-region edges are not physical walls')
  end subroutine verify_interior_wall_viscosity

  subroutine verify_physical_fluxes()
    type(eb_geometry_3d) :: working
    type(eb_plm_boundary_3d) :: boundary
    type(eb_plm_face_prediction_3d) :: result(3),saved(3)
    real(dp) :: ghost(nvar,-2:11,-2:11,-2:11),gt(-2:11,-2:11,-2:11),vf(-2:11,-2:11,-2:11)
    real(dp) :: sound(-2:11,-2:11,-2:11),flat(0:9,0:9,0:9),q(nvar),ql(nvar),qr(nvar),expected(nvar)
    real(dp) :: swap,pface,uface,cav,scale,updated(nvar,2,2,2),updated_t(2,2,2),ledger(nvar),wall(nvar),change
    integer :: ii,jj,kk,mode,stage,axis,side,idx(3),cell(3),output(3)
    character(len=24), parameter :: modes(5)=[character(len=24):: &
      'outflow','slip_wall','no_slip_wall','isothermal_wall','symmetry']
    logical :: valid
    call build_axis_plane_eb_geometry_3d(8,8,8,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
      'x',-.1_dp,working,valid)
    call require(valid,'physical boundary working geometry')
    vf=1;flat=1;boundary%lower=4;boundary%upper=5
    do kk=-2,11
    do jj=-2,11
    do ii=-2,11
      q(1)=.31_dp+.001_dp*(ii+jj+kk)
      q(2:4)=[7._dp,-3._dp,5._dp]+.05_dp*[ii,jj,kk];q(5)=pressure
      q(6:)=mass_fractions
      q(6)=q(6)+.0001_dp*(ii+jj+kk);q(nvar)=q(nvar)-.0001_dp*(ii+jj+kk)
      call reactive_primitive_to_conserved(species,q,ghost(:,ii,jj,kk),gt(ii,jj,kk),sound(ii,jj,kk),valid)
      call require(valid,'physical flux input')
    end do
    end do
    end do
    do mode=1,size(modes)
      boundary%faces=modes(mode)
      do stage=1,2
        if(stage==1) then
          call compute_reactive_eb_plm_normal_fluxes_3d(species,working,ghost,gt,vf,1e-7_dp,flat, &
            result,valid,boundary=boundary)
        else
          call compute_reactive_eb_plm_ctu_fluxes_3d(species,working,ghost,gt,vf,1e-7_dp,flat, &
            result,valid,boundary=boundary)
        end if
        call require(valid,'physical normal and corrected CTU fluxes')
        do axis=1,3
        do side=1,3
          idx=4
          if(side==1) idx(axis)=3
          if(side==2) idx(axis)=5
          output=idx-2*(stage-1)
          ql=result(axis)%left(:,output(1),output(2),output(3))
          qr=result(axis)%right(:,output(1),output(2),output(3))
          ! Independently apply the pinned local-argument rule. Species remain
          ! donor-specific; interior faces (side 3) are deliberately untouched.
          if(side==1) then
            ql(1:5)=qr(1:5)
            if(mode>1) ql(axis+1)=-qr(axis+1)
          else if(side==2) then
            qr(1:5)=ql(1:5)
            if(mode>1) qr(axis+1)=-ql(axis+1)
          end if
          swap=ql(2);ql(2)=ql(axis+1);ql(axis+1)=swap
          swap=qr(2);qr(2)=qr(axis+1);qr(axis+1)=swap
          cell=idx;cell(axis)=cell(axis)+1
          cav=.5_dp*(sound(idx(1),idx(2),idx(3))+sound(cell(1),cell(2),cell(3)))
          call reactive_pelec_primitive_flux_x(species,ql,qr,cav,expected,valid, &
            interface_pressure=pface,interface_velocity=uface)
          call require(valid,'physical reference Riemann arguments')
          swap=expected(2);expected(2)=expected(axis+1);expected(axis+1)=swap
          scale=max(1._dp,maxval(abs(expected)))
          call require(maxval(abs(result(axis)%flux(:,output(1),output(2),output(3))-expected))<1e-13_dp*scale, &
            'six-face Riemann correction preserves species and interior fluxes')
          call require(abs(result(axis)%pressure(output(1),output(2),output(3))-pface)<1e-10_dp.and. &
            abs(result(axis)%velocity(output(1),output(2),output(3))-uface)<1e-12_dp, &
            'physical Godunov pressure and velocity')
        end do
        end do
      end do
      call advance_reactive_eb_plm_ctu_3d(species,working,ghost,gt,vf,1e-7_dp,flat,updated,updated_t,valid, &
        boundary=boundary,boundary_change=ledger,wall_change=wall)
      call require(valid,'physical boundary complete CTU update')
      do ii=1,nvar
        change=sum(updated(ii,:,:,:)-ghost(ii,4:5,4:5,4:5))/8._dp**3
        scale=max(1._dp,sum(abs(ghost(ii,4:5,4:5,4:5)))/8._dp**3)
        call require(abs(change-ledger(ii)-wall(ii))<2e-12_dp*scale,'physical update inventory balance')
      end do
    end do
    saved=result;boundary%faces(1)='invalid'
    call compute_reactive_eb_plm_ctu_fluxes_3d(species,working,ghost,gt,vf,1e-7_dp,flat,result,valid,boundary=boundary)
    call require(.not.valid,'invalid physical boundary rejection')
    do axis=1,3
      call require(all(result(axis)%flux==saved(axis)%flux),'physical boundary failure preserves output')
    end do
  end subroutine verify_physical_fluxes

  subroutine verify_flattening()
    real(dp) :: ghost(nvar,-2:11,-2:11,-2:11),gt(-2:11,-2:11,-2:11),vf(-2:11,-2:11,-2:11)
    real(dp) :: flat(0:9,0:9,0:9),q(nvar),cs,nan_value,vel(3),p
    integer :: axis,ii,jj,kk,idx(3),offset,point(3)
    logical :: valid
    nan_value=ieee_value(0._dp,ieee_quiet_nan)
    do axis=1,3
      vf=1
      do kk=-2,11
      do jj=-2,11
      do ii=-2,11
        idx=[ii,jj,kk];p=100000._dp
        if(idx(axis)>=5) p=200000._dp
        vel=0;vel(axis)=real(5-idx(axis),dp)
        q(1)=.31_dp;q(2:4)=vel;q(5)=p;q(6:)=mass_fractions
        call reactive_primitive_to_conserved(species,q,ghost(:,ii,jj,kk),gt(ii,jj,kk),cs,valid)
        call require(valid,'flattening input state')
      end do
      end do
      end do
      call compute_reactive_eb_flattening_3d(species,ghost,gt,vf,[8,8,8],flat,valid)
      call require(valid.and.flat(4,4,4)==0,'compressive shock flattening in each axis')
      call require(flat(1,1,1)==1.and.flat(0,0,0)==1,'smooth and unused donors retain unity')
      do offset=-3,3
        point=4;point(axis)=point(axis)+offset
        vf(point(1),point(2),point(3))=.5_dp
        call compute_reactive_eb_flattening_3d(species,ghost,gt,vf,[8,8,8],flat,valid)
        call require(valid.and.flat(4,4,4)==1,'each cut neighbor disables directional shock sensor')
        vf(point(1),point(2),point(3))=1
      end do
      ! Invalid active thermodynamics must preserve the previous output.
      ghost(:,4,4,4)=nan_value;flat=-7
      call compute_reactive_eb_flattening_3d(species,ghost,gt,vf,[8,8,8],flat,valid)
      call require(.not.valid.and.all(flat==-7),'flattening transactional invalid state')
      vf=0;ghost=nan_value;gt=nan_value
      call compute_reactive_eb_flattening_3d(species,ghost,gt,vf,[8,8,8],flat,valid)
      call require(valid.and.all(flat==1),'covered flattening ignores poisoned storage')
    end do
  end subroutine verify_flattening

  subroutine verify_ctu_update()
    type(eb_geometry_3d) :: working,physical
    type(eb_integrated_flux_3d) :: integral,saved_flux
    real(dp) :: ghost(nvar,-2:11,-2:11,-2:11),gt(-2:11,-2:11,-2:11),vf(-2:11,-2:11,-2:11),flat(0:9,0:9,0:9)
    real(dp) :: initial(nvar,2,2,2),q(nvar,2,2,2),t(2,2,2),saved(nvar,2,2,2),saved_t(2,2,2)
    real(dp) :: boundary(nvar),wall(nvar),change,scale,vel(3),cs,position,nan_value
    integer :: cut,mode,ii,jj,kk,c,idx(3),inside(3)
    character(len=1), parameter :: axes(3)=['x','y','z']
    character(len=12), parameter :: methods(3)=[character(len=12)::'none','flux_redist','state_redist']
    logical :: valid
    flat=1;nan_value=ieee_value(0._dp,ieee_quiet_nan)
    do cut=0,3
      position=.46_dp
      if(cut==0) position=-.1_dp
      call build_axis_plane_eb_geometry_3d(8,8,8,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
        axes(max(1,cut)),position,working,valid)
      call require(valid,'CTU update working geometry')
      call extract_eb_geometry_3d(working,[4,4,4],[5,5,5],physical,valid)
      call require(valid,'CTU physical geometry')
      ghost=nan_value;gt=nan_value
      do kk=-2,11
      do jj=-2,11
      do ii=-2,11
        idx=[ii,jj,kk];inside=max(1,min(8,idx))
        vf(ii,jj,kk)=working%volume_fraction(inside(1),inside(2),inside(3))
        if(vf(ii,jj,kk)==0) cycle
        vel=[4._dp,2._dp,-1._dp]-.1_dp*inside
        call build_cell_state(.31_dp+.001_dp*sum(inside),vel,ghost(:,ii,jj,kk),gt(ii,jj,kk),cs)
      end do
      end do
      end do
      call compute_reactive_eb_flattening_3d(species,ghost,gt,vf,[8,8,8],flat,valid)
      call require(valid,'CTU update flattening preparation')
      initial=ghost(:,4:5,4:5,4:5)
      do kk=1,2
      do jj=1,2
      do ii=1,2
        if(physical%volume_fraction(ii,jj,kk)==0) initial(:,ii,jj,kk)=0
      end do
      end do
      end do
      do mode=1,3
        q=-17;t=-18;boundary=19;wall=20
        call advance_reactive_eb_plm_ctu_3d(species,working,ghost,gt,vf,1e-8_dp,flat,q,t,valid, &
          redistribution=trim(methods(mode)),integrated_flux=integral,boundary_change=boundary,wall_change=wall)
        call require(valid.and.any(q/=initial),'complete CTU hydro update and redistribution')
        do c=1,nvar
          change=sum((q(c,:,:,:)-initial(c,:,:,:))*physical%volume_fraction)*physical%dx*physical%dy*physical%dz
          scale=max(1._dp,sum(abs(initial(c,:,:,:))*physical%volume_fraction)*physical%dx*physical%dy*physical%dz)
          call require(abs(change-boundary(c)-wall(c))<2e-12_dp*scale,'CTU boundary/wall inventory closure')
        end do
        call require(maxval(abs(sum(q(6:,:,:,:),dim=1)-q(1,:,:,:)))<2e-13_dp,'CTU species closure')
        saved=q;saved_t=t;saved_flux=integral
        if(cut==0.and.mode==1) then
          call advance_reactive_eb_plm_ctu_3d(species,working,ghost,gt,vf,1e-8_dp,flat,q,t,valid,difmag=.1_dp)
          call require(valid.and.all(q==saved).and.all(t==saved_t),'CTU omitted viscosity reference policy')
          call advance_reactive_eb_plm_ctu_3d(species,working,ghost,gt,vf,1e-8_dp,flat,q,t,valid,difmag=0._dp)
          call require(valid.and.any(q/=saved),'CTU zero viscosity disables compression correction')
          q=saved;t=saved_t
        end if
        boundary=19;wall=20
        call advance_reactive_eb_plm_ctu_3d(species,working,ghost,gt,vf,1._dp,flat,q,t,valid, &
          redistribution=trim(methods(mode)),integrated_flux=integral,boundary_change=boundary,wall_change=wall)
        call require(.not.valid.and.all(q==saved).and.all(t==saved_t).and.all(boundary==19).and.all(wall==20), &
          'CTU failed update preserves state and ledgers')
        call require(all(integral%x==saved_flux%x).and.all(integral%y==saved_flux%y).and. &
          all(integral%z==saved_flux%z),'CTU failed update preserves integrated flux')
      end do
    end do
  end subroutine verify_ctu_update
  subroutine verify_ctu_fluxes()
    type(eb_geometry_3d) :: shape
    type(eb_plm_face_prediction_3d) :: result(3),saved(3)
    real(dp) :: ghost(nvar,-2:11,-2:11,-2:11),gt(-2:11,-2:11,-2:11),vf(-2:11,-2:11,-2:11)
    real(dp) :: flat(0:9,0:9,0:9),source(nvar,8,8,8),esource(8,8,8),conserved(nvar,8,8,8)
    real(dp) :: vel(3),rho,cs,tc,position,expected(nvar),u(nvar),dt,face_rho,coord(3),nan_value
    integer :: cut,axis,ii,jj,kk,inside(3),idx(3),lo(3),hi(3)
    logical :: valid
    character(len=1), parameter :: axes(3)=['x','y','z']
    vel=[7._dp,-3._dp,5._dp];dt=1e-6_dp;flat=1;nan_value=ieee_value(0._dp,ieee_quiet_nan)
    do cut=0,3
      position=.37_dp
      if(cut==0) position=-.1_dp
      call build_axis_plane_eb_geometry_3d(8,8,8,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
        axes(max(1,cut)),position,shape,valid)
      call require(valid,'CTU working geometry')
      ghost=nan_value;gt=nan_value
      do kk=-2,11
      do jj=-2,11
      do ii=-2,11
        idx=[ii,jj,kk];inside=max(1,min(8,idx))
        vf(ii,jj,kk)=shape%volume_fraction(inside(1),inside(2),inside(3))
        if(vf(ii,jj,kk)==0.or.count(idx<1.or.idx>8)>1) cycle
        rho=.31_dp
        if(cut==0) rho=rho+.002_dp*sum(inside*[1,2,3])
        call build_cell_state(rho,vel,ghost(:,ii,jj,kk),gt(ii,jj,kk),cs)
      end do
      end do
      end do
      call compute_reactive_eb_plm_ctu_fluxes_3d(species,shape,ghost,gt,vf,dt,flat,result,valid)
      call require(valid,'EB twelve-solve CTU flux pipeline')
      do axis=1,3
        lo=1;lo(axis)=0;hi=4
        call require(all(lbound(result(axis)%flux)==[1,lo]).and.all(ubound(result(axis)%flux)==[nvar,hi]), &
          'CTU cropped output bounds')
        if(cut==0) then
          coord=4;coord(axis)=4.5_dp
          face_rho=.31_dp+.002_dp*sum(coord*[1._dp,2._dp,3._dp])-.5_dp*dt*.002_dp*8*sum(vel*[1._dp,2._dp,3._dp])
          call require(abs(result(axis)%flux(1,2,2,2)-vel(axis)*face_rho)<2e-12_dp, &
            'unsplit linear density advection through transverse corner fluxes')
        else
          call build_cell_state(.31_dp,vel,u,tc,cs)
          expected=u*vel(axis);expected(axis+1)=expected(axis+1)+pressure;expected(5)=expected(5)+pressure*vel(axis)
          do kk=lo(3),hi(3)
          do jj=lo(2),hi(2)
          do ii=lo(1),hi(1)
            if(result(axis)%pressure(ii,jj,kk)==0) then
              call require(all(result(axis)%flux(:,ii,jj,kk)==0),'closed final CTU flux')
            else
              call require(maxval(abs(result(axis)%flux(:,ii,jj,kk)-expected)/max(1._dp,abs(expected)))<3e-13_dp, &
                'cut-cell uniform CTU intensive Euler flux')
            end if
          end do
          end do
          end do
        end if
      end do
      saved=result;source=0;esource=0
      call compute_reactive_eb_plm_ctu_fluxes_3d(species,shape,ghost,gt,vf,dt,flat,result,valid, &
        source=source,source_energy=esource)
      call require(valid,'CTU zero primitive source')
      do axis=1,3
        call require(all(result(axis)%flux==saved(axis)%flux),'zero-source CTU exact parity')
      end do
      source=0;source(2,:,:,:)=8;esource=0;conserved=nan_value
      do kk=1,8
      do jj=1,8
      do ii=1,8
        if(vf(ii,jj,kk)==0) cycle
        conserved(:,ii,jj,kk)=0
        conserved(2,ii,jj,kk)=8*ghost(1,ii,jj,kk)
        conserved(5,ii,jj,kk)=vel(1)*conserved(2,ii,jj,kk)
      end do
      end do
      end do
      call compute_reactive_eb_plm_ctu_fluxes_3d(species,shape,ghost,gt,vf,dt,flat,result,valid, &
        source=source,source_energy=esource)
      call require(valid,'analytic acceleration primitive source')
      call require(any(result(1)%flux/=saved(1)%flux),'source changes CTU flux')
      saved=result
      call compute_reactive_eb_plm_ctu_fluxes_3d(species,shape,ghost,gt,vf,dt,flat,result,valid,conserved_source=conserved)
      call require(valid,'conserved source with poisoned covered cells')
      do axis=1,3
        call require(maxval(abs(result(axis)%flux-saved(axis)%flux)/max(1._dp,abs(saved(axis)%flux)))<2e-12_dp, &
          'conserved acceleration matches analytic primitive source flux')
      end do
      saved=result
      conserved(2,8,8,8)=nan_value
      call compute_reactive_eb_plm_ctu_fluxes_3d(species,shape,ghost,gt,vf,dt,flat,result,valid,conserved_source=conserved)
      call require(.not.valid,'invalid active conserved source rejected')
      call compute_reactive_eb_plm_ctu_fluxes_3d(species,shape,ghost,gt,vf,1._dp,flat,result,valid)
      call require(.not.valid,'CTU failed stage rejects')
      do axis=1,3
        call require(all(result(axis)%flux==saved(axis)%flux).and.all(result(axis)%left==saved(axis)%left).and. &
          all(result(axis)%right==saved(axis)%right).and.all(result(axis)%pressure==saved(axis)%pressure), &
          'CTU stage failure rollback')
      end do
    end do
  end subroutine verify_ctu_fluxes
  subroutine verify_transverse_prediction()
    real(dp) :: initial(nvar,2),q(nvar,2),e(2),original_energy(2),flx(nvar,2,2),p(2,2),v(2,2)
    real(dp) :: ta(2,2),na(2),coef(2),delta(nvar),wanted(nvar),source(nvar),base(nvar,2),base_e(2),nan_value
    integer :: nt,mask,normal_mask,d,side
    logical :: valid,active
    initial=0;initial(1,:)=2;initial(2,:)=1;initial(3,:)=2;initial(4,:)=3;initial(5,:)=100
    initial(6,:)=.25_dp;initial(nvar,:)=.75_dp;original_energy=50
    delta=0;delta(1:5)=[.2_dp,.1_dp,.2_dp,.3_dp,.4_dp];delta(6)=.05_dp;delta(nvar)=.15_dp
    coef=[.1_dp,.2_dp];nan_value=ieee_value(0._dp,ieee_quiet_nan)
    do nt=1,2
      do mask=0,2**(2*nt)-1
        do d=1,nt
          ta(:,d)=[merge(.3_dp,0._dp,btest(mask,2*d-2)),merge(.7_dp,0._dp,btest(mask,2*d-1))]
        end do
        active=all(ta(:,:nt)>0)
        flx=0;flx(:,2,1)=10*delta;flx(:,2,2)=10*delta
        p(:,1)=[100._dp,110._dp];p(:,2)=[80._dp,90._dp]
        v(:,1)=[2._dp,3._dp];v(:,2)=[1._dp,2._dp]
        if(.not.active) then
          flx=nan_value;p=nan_value;v=nan_value
        end if
        do normal_mask=0,3
          na=[merge(.2_dp,0._dp,btest(normal_mask,0)),merge(.8_dp,0._dp,btest(normal_mask,1))]
          q=initial;e=original_energy
          call apply_reactive_eb_transverse_3d(q,e,1.4_dp,na,ta(:,:nt),flx(:,:,:nt),p(:,:nt),v(:,:nt), &
            coef(:nt),valid)
          call require(valid,'single/double EB transverse correction')
          do side=1,2
            if(na(side)==0.or..not.active) then
              call require(all(q(:,side)==initial(:,side)).and.e(side)==50,'closed normal or transverse gate')
            else
              wanted=delta*merge(1._dp,3._dp,nt==1)
              call require(abs(q(1,side)-(2-wanted(1)))<1e-14_dp,'transverse density conservation')
              call require(maxval(abs(q(1,side)*q(2:4,side)-([2._dp,4._dp,6._dp]-wanted(2:4))))<1e-14_dp, &
                'transverse momentum conservation')
              call require(abs(e(side)+.5_dp*q(1,side)*sum(q(2:4,side)**2)-(64-wanted(5)))<2e-14_dp, &
                'transverse total energy conservation')
              call require(abs(q(5,side)-merge(82.8_dp,56._dp,nt==1))<2e-14_dp,'ordered transverse pressure work')
              call require(maxval(abs(q(1,side)*q(6:,side)-(2*initial(6:,side)-wanted(6:))))<1e-14_dp, &
                'transverse species conservation')
            end if
          end do
        end do
      end do
    end do
    q=initial;e=original_energy
    call apply_reactive_eb_transverse_3d(q,e,1.4_dp,[1._dp,1._dp],ta,flx,p,v,coef,valid)
    call require(valid,'source-free transverse reference')
    base=q;base_e=e
    source=0;source(1:5)=[.2_dp,.4_dp,.6_dp,.8_dp,2._dp];source(6)=.01_dp;source(nvar)=-.01_dp
    q=initial;e=original_energy
    call apply_reactive_eb_transverse_3d(q,e,1.4_dp,[1._dp,1._dp],ta,flx,p,v,coef,valid,source,3._dp,.5_dp)
    call require(valid.and.maxval(abs(q(:,1)-base(:,1)-.5_dp*source))<1e-14_dp.and. &
      all(e==base_e+1.5_dp),'double transverse primitive source integration')
    q=initial;e=original_energy;flx(1,2,1)=1e6_dp
    call apply_reactive_eb_transverse_3d(q,e,1.4_dp,[1._dp,1._dp],ta,flx,p,v,coef,valid)
    call require(.not.valid.and.all(q==initial).and.all(e==original_energy),'nonpositive transverse density rollback')
  end subroutine verify_transverse_prediction
  subroutine verify_normal_fluxes()
    type(eb_geometry_3d) :: shape
    type(eb_plm_face_prediction_3d) :: result(3),saved(3)
    real(dp) :: ghost(nvar,-2:7,-2:7,-2:7),gt(-2:7,-2:7,-2:7),vf(-2:7,-2:7,-2:7),flat(0:5,0:5,0:5)
    real(dp) :: fx(nvar,0:4,4,4),fy(nvar,4,0:4,4),fz(nvar,4,4,0:4),expected(nvar),cell_state(nvar)
    real(dp) :: tc,cs,vel(3),area,nan_value
    integer :: cut,mode,axis,ii,jj,kk,idx(3),inside(3),lo(3),hi(3)
    logical :: valid
    character(len=1), parameter :: axes(3)=['x','y','z']
    nan_value=ieee_value(0._dp,ieee_quiet_nan);flat=1;vel=[7._dp,-3._dp,5._dp]
    do cut=1,3
      call build_axis_plane_eb_geometry_3d(4,4,4,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
        axes(cut),.37_dp,shape,valid)
      call require(valid,'normal flux geometry')
      do mode=1,2
        ghost=nan_value;gt=nan_value
        do kk=-2,7
        do jj=-2,7
        do ii=-2,7
          idx=[ii,jj,kk];inside=max(1,min(4,idx))
          vf(ii,jj,kk)=shape%volume_fraction(inside(1),inside(2),inside(3))
          if(vf(ii,jj,kk)==0.or.count(idx<1.or.idx>4)>1) cycle
          call build_cell_state(.31_dp+.002_dp*(mode-1)*sum(inside),vel,ghost(:,ii,jj,kk),gt(ii,jj,kk),cs)
        end do
        end do
        end do
        call compute_reactive_eb_plm_normal_fluxes_3d(species,shape,ghost,gt,vf,1e-7_dp,flat,result,valid)
        call require(valid,'normal flux pipeline ignores covered and unused corner states')
        if(mode==1) then
          call build_cell_state(.31_dp,vel,cell_state,tc,cs)
          do axis=1,3
            expected=cell_state*vel(axis);expected(axis+1)=expected(axis+1)+pressure
            expected(5)=expected(5)+pressure*vel(axis)
            lo=1;lo(axis)=0;hi=4
            do kk=lo(3),hi(3)
            do jj=lo(2),hi(2)
            do ii=lo(1),hi(1)
              select case(axis)
              case(1);area=shape%x_face_fraction(ii,jj,kk)
              case(2);area=shape%y_face_fraction(ii,jj,kk)
              case(3);area=shape%z_face_fraction(ii,jj,kk)
              end select
              if(area==0) then
                call require(all(result(axis)%flux(:,ii,jj,kk)==0).and.result(axis)%pressure(ii,jj,kk)==0.and. &
                  result(axis)%velocity(ii,jj,kk)==0,'closed normal flux zero')
              else
                call require(maxval(abs(result(axis)%flux(:,ii,jj,kk)-expected)/max(1._dp,abs(expected)))<3e-13_dp, &
                  'uniform intensive flux including partially open faces')
                call require(abs(result(axis)%pressure(ii,jj,kk)/pressure-1)<3e-13_dp.and. &
                  abs(result(axis)%velocity(ii,jj,kk)-vel(axis))<3e-13_dp,'Godunov pressure and velocity')
              end if
            end do
            end do
            end do
          end do
        else
          saved=result
          call compute_reactive_eb_plm_normal_fluxes_3d(species,shape,ghost,gt,vf,1e-7_dp,flat,result,valid,order=1)
          call require(valid.and.any(result(1)%flux/=saved(1)%flux),'PLM changes nonuniform normal flux')
          call reactive_eb_outflow_riemann_fluxes_3d(species,ghost(:,1:4,1:4,1:4),gt(1:4,1:4,1:4),shape,'pelec', &
            fx,fy,fz,valid,ghost(:,0:5,0:5,0:5),gt(0:5,0:5,0:5),difmag=0._dp)
          call require(valid.and.maxval(abs(result(1)%flux-fx)/max(1._dp,abs(fx)))<3e-13_dp.and. &
            maxval(abs(result(2)%flux-fy)/max(1._dp,abs(fy)))<3e-13_dp.and. &
            maxval(abs(result(3)%flux-fz)/max(1._dp,abs(fz)))<3e-13_dp,'first-order PCM flux limit')
        end if
        saved=result;ghost(1,5,4,4)=-1
        call compute_reactive_eb_plm_normal_fluxes_3d(species,shape,ghost,gt,vf,1e-7_dp,flat,result,valid)
        call require(.not.valid,'bad open donor rejected')
        do axis=1,3
          call require(all(result(axis)%flux==saved(axis)%flux).and.all(result(axis)%left==saved(axis)%left).and. &
            all(result(axis)%right==saved(axis)%right).and.all(result(axis)%energy_left==saved(axis)%energy_left).and. &
            all(result(axis)%energy_right==saved(axis)%energy_right).and. &
            all(result(axis)%pressure==saved(axis)%pressure).and.all(result(axis)%velocity==saved(axis)%velocity), &
            'normal prediction and flux transaction rollback')
        end do
      end do
    end do
  end subroutine verify_normal_fluxes
  subroutine verify_normal_prediction()
    real(dp) :: q(nvar,-2:2),e(-2:2),qm(nvar),qp(nvar),em,ep,area(2),nan_value
    logical :: covered(-2:2),valid
    integer :: axis,mask,j
    nan_value=ieee_value(0._dp,ieee_quiet_nan)
    do axis=1,3
      do mask=0,3
        do j=-2,2
          q(:,j)=0;q(1,j)=1+.02_dp*j;q(5,j)=135000
          q(2:4,j)=[3._dp,5._dp,7._dp]+j*[.1_dp,.2_dp,.3_dp]
          q(nvar,j)=1;e(j)=300000+200*j
        end do
        covered=.false.;area=[merge(.7_dp,0._dp,btest(mask,0)),merge(.4_dp,0._dp,btest(mask,1))]
        call predict_reactive_eb_plm_cell_3d(q,e,400._dp,axis,.01_dp,0._dp,1._dp,covered,area, &
          qm,qp,em,ep,valid)
        call require(valid,'EB PLM normal predictor')
        if(area(1)>0) then
          call require(abs(qm(1)-.99_dp)<1e-14_dp.and.abs(em-299900)<1e-10_dp,'lower density/internal energy trace')
          call require(maxval(abs(qm(2:4)-([3._dp,5._dp,7._dp]-.5_dp*[.1_dp,.2_dp,.3_dp])))<1e-14_dp, &
            'lower rotated velocities')
        else
          call require(all(qm==0).and.em==0,'closed lower face zero')
        end if
        if(area(2)>0) then
          call require(abs(qp(1)-1.01_dp)<1e-14_dp.and.abs(ep-300100)<1e-10_dp,'upper density/internal energy trace')
          call require(maxval(abs(qp(2:4)-([3._dp,5._dp,7._dp]+.5_dp*[.1_dp,.2_dp,.3_dp])))<1e-14_dp, &
            'upper rotated velocities')
        else
          call require(all(qp==0).and.ep==0,'closed upper face zero')
        end if
      end do
      covered(-2)=.true.;q(:,-2)=nan_value;e(-2)=nan_value
      call predict_reactive_eb_plm_cell_3d(q,e,400._dp,axis,.01_dp,0._dp,1._dp,covered,[1._dp,1._dp], &
        qm,qp,em,ep,valid)
      call require(valid.and.abs(em-(300000._dp-100._dp*4/3))<1e-10_dp.and. &
        abs(ep-(300000._dp+100._dp*4/3))<1e-10_dp,'EB energy slope order fallback')
      call predict_reactive_eb_plm_cell_3d(q,e,400._dp,axis,.01_dp,1e-6_dp,0._dp,covered,[1._dp,1._dp], &
        qm,qp,em,ep,valid)
      call require(valid.and.all(qm==q(:,0)).and.all(qp==q(:,0)).and.em==e(0).and.ep==e(0), &
        'flattened finite-time EB predictor')
      call predict_reactive_eb_plm_cell_3d(q,e,400._dp,axis,.01_dp,1._dp,1._dp,covered,[1._dp,1._dp], &
        qm,qp,em,ep,valid)
      call require(.not.valid.and.all(qm==0).and.all(qp==0).and.em==0.and.ep==0,'normal CFL zero outputs')
      covered=.false.
      do j=-2,2
        q(:,j)=0;q(1,j)=1+.02_dp*j;q(5,j)=135000;q(nvar,j)=1
        q(axis+1,j)=-20;e(j)=300000+200*j
      end do
      call predict_reactive_eb_plm_cell_3d(q,e,400._dp,axis,.01_dp,1e-6_dp,1._dp,covered,[1._dp,1._dp], &
        qm,qp,em,ep,valid)
      call require(valid.and.abs(em-299900.2_dp)<1e-10_dp.and.abs(ep-300096.2_dp)<1e-10_dp, &
        'finite-time internal energy contact-wave integration')
      ! Independent limited species slopes need not sum to zero. Pinned
      ! pc_plm_d_eb_passive retains the resulting nonunit face weights.
      q=0;q(1,:)=1;q(5,:)=135000;e=300000
      q(6:8,-2)=[.1_dp,.3_dp,.6_dp];q(6:8,-1)=q(6:8,-2)
      q(6:8,0)=[.2_dp,.2_dp,.6_dp]
      q(6:8,1)=[.4_dp,.1_dp,.5_dp];q(6:8,2)=q(6:8,1)
      do mask=0,3
        area=[merge(.7_dp,0._dp,btest(mask,0)),merge(.4_dp,0._dp,btest(mask,1))]
        call predict_reactive_eb_plm_cell_3d(q,e,400._dp,axis,.01_dp,0._dp,1._dp,covered,area, &
          qm,qp,em,ep,valid,order=2)
        call require(valid,'independent species EB prediction')
        if(area(1)>0) then
          call require(maxval(abs(qm(6:8)-[.1_dp,4._dp/15,.6_dp]))<1e-14_dp, &
            'lower EB species retain independent slopes')
        else
          call require(all(qm==0).and.em==0,'independent species lower closed face')
        end if
        if(area(2)>0) then
          call require(maxval(abs(qp(6:8)-[.3_dp,2._dp/15,.6_dp]))<1e-14_dp, &
            'upper EB species retain independent slopes')
        else
          call require(all(qp==0).and.ep==0,'independent species upper closed face')
        end if
      end do
    end do
    q=nan_value;e=nan_value;covered=.true.
    call predict_reactive_eb_plm_cell_3d(q,e,400._dp,1,.01_dp,1e-6_dp,1._dp,covered,[0._dp,0._dp], &
      qm,qp,em,ep,valid)
    call require(valid.and.all(qm==0).and.all(qp==0).and.em==0.and.ep==0,'fully closed predictor ignores covered state')
  end subroutine verify_normal_prediction

  subroutine verify_viscosity()
    real(dp), allocatable :: base(:,:,:,:),rhs(:,:,:,:),saved(:,:,:,:),saved_t(:,:,:)
    real(dp), parameter :: step=1e-8_dp
    call build_axis_plane_eb_geometry_3d(nx,ny,nz,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
      'x',.37_dp,geometry,ok)
    call require(ok,'viscosity geometry')
    state=0;temperature=0
    do k=1,nz
    do j=1,ny
    do i=1,nx
      if(geometry%volume_fraction(i,j,k)==0) cycle
      call build_cell_state(.31_dp+.001_dp*i,[-real(i,dp),-2._dp*j,-3._dp*k], &
        state_cell,base_temperature,sound_speed)
      state(:,i,j,k)=state_cell;temperature(i,j,k)=base_temperature
    end do
    end do
    end do
    saved=state;saved_t=temperature
    call reactive_eb_outflow_riemann_fluxes_3d(species,state,temperature,geometry,'pelec', &
      x_flux,y_flux,z_flux,ok,difmag=0._dp)
    call require(ok,'zero EB coefficient')
    base=x_flux
    call reactive_eb_outflow_riemann_fluxes_3d(species,state,temperature,geometry,'pelec', &
      x_flux,y_flux,z_flux,ok,difmag=.037_dp)
    call require(ok.and.any(x_flux/=base),'EB viscosity changes open flux')
    call require(all(x_flux(:,0:3,:,:)==0),'EB viscosity preserves closed flux')
    allocate(rhs,mold=state)
    call reactive_eb_flux_divergence_3d(species,state,temperature,geometry,x_flux,y_flux,z_flux,rhs,ok)
    call require(ok,'corrected EB divergence')
    call advance_reactive_eb_euler_3d(species,state,temperature,geometry,'pelec',step,new_state,new_temperature,ok, &
      difmag=.037_dp)
    call require(ok,'EB viscosity Euler step')
    call require(maxval(abs(new_state-(state+step*rhs))/max(1._dp,abs(new_state)))<2e-14_dp, &
      'Euler uses corrected intensive flux with geometry divergence')
    call advance_reactive_eb_redistributed_euler_3d(species,state,temperature,geometry,'pelec', &
      step,new_state,new_temperature,ok,difmag=.037_dp)
    call require(ok,'EB viscosity flux redistribution')
    call advance_reactive_eb_state_redistributed_euler_3d(species,state,temperature,geometry,'pelec', &
      step,new_state,new_temperature,ok,difmag=.037_dp)
    call require(ok,'EB viscosity state redistribution')
    call advance_reactive_eb_euler_3d(species,state,temperature,geometry,'pelec',step,new_state,new_temperature,ok, &
      difmag=-1._dp)
    call require(.not.ok.and.all(new_state==saved).and.all(new_temperature==saved_t),'negative EB coefficient rollback')
  end subroutine verify_viscosity

  subroutine verify_small_cell_hydro_relief()
    real(dp), allocatable :: raw_rhs(:, :, :, :)
    real(dp) :: stable_dt, local_conservation_error
    logical :: local_ok

    call build_axis_plane_eb_geometry_3d( &
      nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      "x", 0.395_dp, geometry, local_ok)
    call require(local_ok, "small-cell hydro geometry")
    call assert_close(geometry%volume_fraction(4, 4, 3), 0.05_dp, &
      3.0e-13_dp, "hydro small-cell fraction")
    call build_cell_state(0.31_dp, [0.0_dp, 0.0_dp, 0.0_dp], &
      state_cell, base_temperature, sound_speed)
    call fill_state(state_cell, base_temperature)
    call build_cell_state(0.93_dp, [0.0_dp, 0.0_dp, 0.0_dp], &
      perturbed_cell, perturbed_temperature, perturbed_sound_speed)
    state(:, 4, 4, 3) = perturbed_cell
    temperature(4, 4, 3) = perturbed_temperature
    stable_dt = 0.5_dp / (max(sound_speed, perturbed_sound_speed) * &
      (1.0_dp / geometry%dx + 1.0_dp / geometry%dy + &
        1.0_dp / geometry%dz))
    small_cell_full_step = stable_dt

    call reactive_eb_outflow_riemann_fluxes_3d( &
      species, state, temperature, geometry, "rusanov", &
      x_flux, y_flux, z_flux, local_ok)
    call require(local_ok, "small-cell raw face fluxes")
    allocate(raw_rhs(nvar, nx, ny, nz))
    call reactive_eb_flux_divergence_3d( &
      species, state, temperature, geometry, x_flux, y_flux, z_flux, &
      raw_rhs, local_ok)
    call require(local_ok, "small-cell raw flux divergence")
    call require(state(irho, 4, 4, 3) + &
      stable_dt * raw_rhs(irho, 4, 4, 3) < 0.0_dp, &
      "raw full-cell timestep produces negative cut density")

    call advance_reactive_eb_euler_3d( &
      species, state, temperature, geometry, "rusanov", stable_dt, &
      new_state, new_temperature, local_ok)
    call require(.not. local_ok .and. all(new_state == state) .and. &
      all(new_temperature == temperature), &
      "raw small-cell hydro rejects full-cell timestep")

    call fluid_integrals(state, geometry, before_integrals)
    call advance_reactive_eb_redistributed_euler_3d( &
      species, state, temperature, geometry, "rusanov", stable_dt, &
      new_state, new_temperature, local_ok)
    call require(local_ok, &
      "redistribution accepts full-cell-scale hydro timestep")
    call require(new_state(irho, 4, 4, 3) > 0.0_dp, &
      "redistributed small-cell density remains positive")
    redistributed_cut_density = new_state(irho, 4, 4, 3)
    call fluid_integrals(new_state, geometry, after_integrals)
    local_conservation_error = maxval( &
      abs(after_integrals - before_integrals) / &
        max(1.0_dp, abs(before_integrals)))
    call require(local_conservation_error <= 3.0e-13_dp, &
      "small-cell full-step conservation")

    call advance_reactive_eb_state_redistributed_euler_3d( &
      species, state, temperature, geometry, "rusanov", stable_dt, &
      new_state, new_temperature, local_ok)
    call require(local_ok, &
      "StateRedist accepts full-cell-scale hydro timestep")
    call require(new_state(irho, 4, 4, 3) > 0.0_dp, &
      "StateRedist small-cell density remains positive")
    state_redistributed_cut_density = new_state(irho, 4, 4, 3)
    call fluid_integrals(new_state, geometry, after_integrals)
    local_conservation_error = maxval( &
      abs(after_integrals - before_integrals) / &
        max(1.0_dp, abs(before_integrals)))
    call require(local_conservation_error <= 3.0e-13_dp, &
      "StateRedist small-cell full-step conservation")
  end subroutine verify_small_cell_hydro_relief

  subroutine verify_uniform_tangent_update(axis, solver, velocity)
    character(len=*), intent(in) :: axis, solver
    real(dp), intent(in) :: velocity(3)

    real(dp) :: local_temperature, local_sound_speed
    real(dp) :: local_scale, local_error
    logical :: local_ok

    call build_axis_plane_eb_geometry_3d( &
      nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      axis, 0.37_dp, geometry, local_ok)
    call require(local_ok, trim(axis)//"-plane uniform geometry")
    call build_cell_state( &
      0.31_dp, velocity, state_cell, local_temperature, local_sound_speed)
    call fill_state(state_cell, local_temperature)
    call advance_reactive_eb_euler_3d( &
      species, state, temperature, geometry, solver, 1.0e-7_dp, &
      new_state, new_temperature, local_ok)
    call require(local_ok, trim(axis)//"-plane uniform Euler update")
    local_scale = max(1.0_dp, maxval(abs(state)))
    local_error = maxval(abs(new_state - state)) / local_scale
    maximum_uniform_state_error = max( &
      maximum_uniform_state_error, local_error)
    call require(local_error <= 2.0e-14_dp, &
      trim(axis)//"-plane uniform state invariance")
    local_error = maxval(abs(new_temperature - temperature)) / &
      max(1.0_dp, local_temperature)
    maximum_uniform_temperature_error = max( &
      maximum_uniform_temperature_error, local_error)
    call require(local_error <= 2.0e-11_dp, &
      trim(axis)//"-plane uniform temperature invariance")
    call advance_reactive_eb_redistributed_euler_3d( &
      species, state, temperature, geometry, solver, 1.0e-7_dp, &
      new_state, new_temperature, local_ok)
    call require(local_ok, &
      trim(axis)//"-plane uniform redistributed Euler update")
    local_error = maxval(abs(new_state - state)) / local_scale
    maximum_uniform_state_error = max( &
      maximum_uniform_state_error, local_error)
    call require(local_error <= 2.0e-14_dp, &
      trim(axis)//"-plane redistributed state invariance")
    local_error = maxval(abs(new_temperature - temperature)) / &
      max(1.0_dp, local_temperature)
    maximum_uniform_temperature_error = max( &
      maximum_uniform_temperature_error, local_error)
    call require(local_error <= 2.0e-11_dp, &
      trim(axis)//"-plane redistributed temperature invariance")
    call advance_reactive_eb_state_redistributed_euler_3d( &
      species, state, temperature, geometry, solver, 1.0e-7_dp, &
      new_state, new_temperature, local_ok)
    call require(local_ok, &
      trim(axis)//"-plane uniform StateRedist Euler update")
    local_error = maxval(abs(new_state - state)) / local_scale
    maximum_uniform_state_error = max( &
      maximum_uniform_state_error, local_error)
    call require(local_error <= 2.0e-14_dp, &
      trim(axis)//"-plane StateRedist state invariance")
    local_error = maxval(abs(new_temperature - temperature)) / &
      max(1.0_dp, local_temperature)
    maximum_uniform_temperature_error = max( &
      maximum_uniform_temperature_error, local_error)
    call require(local_error <= 2.0e-11_dp, &
      trim(axis)//"-plane StateRedist temperature invariance")
  end subroutine verify_uniform_tangent_update

  subroutine build_cell_state( &
      density, velocity, conserved, recovered_temperature, recovered_sound)
    real(dp), intent(in) :: density, velocity(3)
    real(dp), intent(out) :: conserved(:)
    real(dp), intent(out) :: recovered_temperature, recovered_sound

    logical :: local_ok
    integer :: species_index

    primitive(1:5) = [density, velocity, pressure]
    do species_index = 1, size(species)
      primitive(reactive_mass_fraction_component(species_index)) = &
        mass_fractions(species_index)
    end do
    call reactive_primitive_to_conserved( &
      species, primitive, conserved, recovered_temperature, &
      recovered_sound, local_ok)
    call require(local_ok, "cell-state construction")
  end subroutine build_cell_state

  subroutine fill_state(conserved, cell_temperature)
    real(dp), intent(in) :: conserved(:), cell_temperature

    integer :: local_i, local_j, local_k

    do local_k = 1, nz
      do local_j = 1, ny
        do local_i = 1, nx
          state(:, local_i, local_j, local_k) = conserved
          temperature(local_i, local_j, local_k) = cell_temperature
        end do
      end do
    end do
  end subroutine fill_state

  subroutine fluid_integrals(field, local_geometry, integrals)
    real(dp), intent(in) :: field(:, :, :, :)
    type(eb_geometry_3d), intent(in) :: local_geometry
    real(dp), intent(out) :: integrals(:)

    real(dp) :: cell_volume
    integer :: component

    cell_volume = local_geometry%dx * local_geometry%dy * local_geometry%dz
    do component = 1, size(field, 1)
      integrals(component) = sum( &
        field(component, :, :, :) * local_geometry%volume_fraction) * &
        cell_volume
    end do
  end subroutine fluid_integrals

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message

    if (abs(actual - expected) > tolerance) error stop message
  end subroutine assert_close

end program test_reactive_eb_hydro_3d
