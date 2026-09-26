program test_reactive_eb_transport_3d
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use eb_geometry_2d_mod, only: eb_geometry_2d,build_eb_geometry_2d
  use state_indices_mod, only: irho, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_species_component, reactive_conserved_to_primitive, &
    reactive_primitive_to_conserved
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d, eb_cut_cell_3d, &
    build_axis_plane_eb_geometry_3d,extrude_eb_geometry_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: &
    initialize_reactive_eb_density_sheet_3d, reactive_eb_integrals_3d, &
    reactive_eb_element_integrals_3d, advance_reactive_eb_strang_3d, &
    advance_reactive_eb_full_3d
  use eb_reactive_transport_3d_mod, only: &
    reactive_eb_transport_timestep_3d, reactive_eb_transport_rhs_3d, reactive_eb_transport_patch_timestep_3d, &
    reactive_eb_transport_fluxes_rhs_3d, reactive_eb_transport_euler_update_3d, &
    advance_reactive_eb_transport_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_reactive_transport_3d_mod, only: compute_sdc_eb_transport_3d
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  logical :: ok

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "transport thermodynamics load")
  call load_h2o2_elementary_mechanism(reactions, ok)
  call require(ok, "transport chemistry load")
  call load_h2o2_elementary_transport(transport, ok)
  call require(ok, "transport database load")
  call check_regular_patch()
  call check_axis(species, transport, "x")
  call check_axis(species, transport, "y")
  call check_axis(species, transport, "z")
  call check_full_driver(species, reactions, transport)
  call check_curved_driver(species,reactions,transport,.false.)
  call check_curved_driver(species,reactions,transport,.true.)
  write(*, '(a)') "test_reactive_eb_transport_3d: PASS"

contains

  subroutine check_regular_patch()
    type(eb_geometry_3d) :: geometry
    type(reactive_eb_3d_config) :: config
    real(dp), allocatable :: state(:,:,:,:),temperature(:,:,:)
    real(dp) :: reference_dt,reference_d,patch_dt,patch_d
    logical :: valid
    config%nx=4;config%ny=4;config%nz=4
    call build_axis_plane_eb_geometry_3d(4,4,4,0._dp,1._dp,0._dp,2._dp,0._dp,3._dp,'x',-1._dp,geometry,valid)
    call require(valid,'regular transport patch geometry')
    allocate(state(reactive_nvar(size(species)),4,4,4),temperature(4,4,4))
    call initialize_reactive_eb_density_sheet_3d(species,config,geometry,state,temperature,valid)
    call require(valid,'regular transport patch state')
    call reactive_eb_transport_timestep_3d(species,transport,state,temperature,geometry,.35_dp, &
      .true.,.true.,.true.,reference_dt,reference_d,valid)
    call require(valid,'regular reference transport bound')
    call reactive_eb_transport_patch_timestep_3d(species,transport,state,temperature,geometry,.35_dp, &
      .true.,.true.,.true.,patch_dt,patch_d,valid)
    call require(valid.and.patch_d==reference_d,'regular patch diffusivity matches existing calculation')
    call require(abs(patch_dt-reference_dt)<2.e-14_dp*reference_dt,'anisotropic regular patch interval')
    state(irho,1,1,1)=-1._dp
    call reactive_eb_transport_patch_timestep_3d(species,transport,state,temperature,geometry,.35_dp, &
      .true.,.true.,.true.,patch_dt,patch_d,valid)
    call require(.not.valid.and.patch_dt==0._dp.and.patch_d==0._dp,'invalid active patch state resets outputs')
  end subroutine check_regular_patch

  subroutine check_curved_driver(species,reactions,transport,cluster)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(eb_geometry_2d) :: section
    type(eb_geometry_3d) :: geometry
    type(reactive_eb_3d_config) :: config
    logical, intent(in) :: cluster
    real(dp) :: phi(0:10,0:8),theta,initial_elements(3),final_elements(3),error
    real(dp) :: primitive(reactive_nprim(size(species))),recovered_temperature,sound_speed
    real(dp), allocatable :: state(:,:,:,:),temperature(:,:,:),result(:,:,:,:),temps(:,:,:)
    real(dp), allocatable :: hydro(:,:,:,:),diffused(:,:,:,:),before(:),after(:)
    integer :: i,j,nvar,mode
    logical :: valid
    do j=0,8
      do i=0,10
        phi(i,j)=(real(i,dp)/10-.5_dp)**2+(real(j,dp)/8-.5_dp)**2-.29_dp**2
      end do
    end do
    call build_eb_geometry_2d(phi,0._dp,1._dp,0._dp,1._dp,section,valid)
    call require(valid,'curved full-driver section')
    call extrude_eb_geometry_3d(section,3,3,0._dp,1._dp,geometry,valid)
    call require(valid,'curved full-driver geometry')
    if (cluster) then
      call build_axis_plane_eb_geometry_3d(10,8,3,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
        'x',.99_dp,geometry,valid)
      call require(valid,'cluster full-driver geometry')
      call require(maxval(geometry%volume_fraction)<.5_dp,'full-driver cluster has no large receiver')
    end if
    config%nx=10;config%ny=8;config%nz=3
    nvar=reactive_nvar(size(species))
    allocate(state(nvar,10,8,3),temperature(10,8,3),result(nvar,10,8,3),temps(10,8,3))
    allocate(hydro(nvar,10,8,3),diffused(nvar,10,8,3),before(nvar),after(nvar))
    call initialize_reactive_eb_density_sheet_3d(species,config,geometry,state,temperature,valid)
    call require(valid,'curved full-driver initial state')
    if (cluster) then
      call reactive_conserved_to_primitive(species,state(:,10,4,2),temperature(10,4,2), &
        primitive,recovered_temperature,sound_speed,valid)
      call require(valid,'cluster primitive recovery')
      primitive(irho)=1.1_dp*primitive(irho)
      call reactive_primitive_to_conserved(species,primitive,state(:,10,4,2),temperature(10,4,2),sound_speed,valid)
      call require(valid,'cluster density perturbation at fixed pressure')
    end if
    call reactive_eb_integrals_3d(state,geometry,before,valid)
    call require(valid,'curved initial inventories')
    call reactive_eb_element_integrals_3d(species,state,geometry,initial_elements,valid)
    call require(valid,'curved initial elements')
    do mode=1,3
      call advance_reactive_eb_full_3d(species,reactions,transport,state,temperature,geometry, &
        'rusanov','state_redist',1.e-8_dp,mode==3,2.e-7_dp,1.e-12_dp, &
        mode>=2,.true.,.true.,.true.,.true.,result,temps,theta,valid,.5_dp)
      call require(valid,'curved hydro/transport/reaction full driver')
      call require(all(temps>0._dp).and.all(result(irho,:,:,:)>0._dp),'curved positive state')
      call require(all(same_real_bits(result,state).or.spread(geometry%volume_fraction>0,1,nvar)), &
        'curved covered storage unchanged')
      call reactive_eb_integrals_3d(result,geometry,after,valid)
      call require(valid,'curved final inventories')
      error=maxval(abs(after([irho,iet])-before([irho,iet]))/max(1.e-20_dp,abs(before([irho,iet]))))
      call require(error<=2.e-9_dp,'curved mass and energy conservation')
      call reactive_eb_element_integrals_3d(species,result,geometry,final_elements,valid)
      call require(valid,'curved final elements')
      call require(maxval(abs(final_elements-initial_elements)/max(1.e-30_dp,abs(initial_elements)))<=2.e-9_dp, &
        'curved element conservation')
      if (mode==1) hydro=result
      if (mode==2) then
        diffused=result
        call require(maxval(abs(diffused-hydro))>1.e-12_dp,'curved transport response')
      end if
      if (mode==3) call require(maxval(abs(result-diffused))>1.e-12_dp,'curved chemistry response')
    end do
  end subroutine check_curved_driver

  subroutine check_axis(species, transport, axis)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    character(len=*), intent(in) :: axis

    type(eb_integrated_flux_3d) :: stage_flux,saved_flux
    type(gas_transport_species) :: bad_transport(size(transport))
    real(dp), allocatable :: stage_rate(:,:,:,:)
    real(dp) :: stage_theta,rate_sum(5+size(species)),rate_scale(5+size(species))
    type(reactive_eb_3d_config) :: config
    type(eb_geometry_3d) :: geometry
    real(dp), allocatable :: state(:, :, :, :), saved_state(:, :, :, :)
    real(dp), allocatable :: temperature(:, :, :), saved_temperature(:, :, :)
    real(dp), allocatable :: initial_integrals(:), final_integrals(:)
    real(dp), allocatable :: rhs(:, :, :, :), primitive(:)
    real(dp), allocatable :: x_flux(:, :, :, :), y_flux(:, :, :, :)
    real(dp), allocatable :: z_flux(:, :, :, :)
    real(dp), allocatable :: bad_state(:, :, :, :), bad_temperature(:, :, :)
    real(dp) :: initial_elements(3), final_elements(3)
    real(dp) :: timestep, modified_timestep, maximum_diffusivity
    real(dp) :: modified_diffusivity, interval, minimum_theta
    real(dp) :: conserved_error, element_error, response
    real(dp) :: recovered_temperature, sound_speed
    logical :: local_ok
    integer :: i, j, k, nvar

    config%nx = 4
    config%ny = 4
    config%nz = 4
    config%x_upper = 1.0e-4_dp
    config%y_upper = 1.0e-4_dp
    config%z_upper = 1.0e-4_dp
    config%plane_axis = axis
    config%plane_position = 4.875e-5_dp
    call build_axis_plane_eb_geometry_3d( &
      config%nx, config%ny, config%nz, config%x_lower, config%x_upper, &
      config%y_lower, config%y_upper, config%z_lower, config%z_upper, &
      axis, config%plane_position, geometry, local_ok)
    call require(local_ok, "transport geometry")
    call require(abs(minval(geometry%volume_fraction, &
      mask=geometry%volume_fraction > 0.0_dp) - 0.05_dp) < 1.0e-12_dp, &
      "transport small-cell geometry")

    nvar = reactive_nvar(size(species))
    allocate(state(nvar, 4, 4, 4), saved_state(nvar, 4, 4, 4))
    allocate(temperature(4, 4, 4), saved_temperature(4, 4, 4))
    allocate(initial_integrals(nvar), final_integrals(nvar))
    allocate(rhs(nvar, 4, 4, 4), primitive(reactive_nprim(size(species))))
    allocate(x_flux(nvar, 0:geometry%nx, geometry%ny, geometry%nz))
    allocate(y_flux(nvar, geometry%nx, 0:geometry%ny, geometry%nz))
    allocate(z_flux(nvar, geometry%nx, geometry%ny, 0:geometry%nz))
    call initialize_reactive_eb_density_sheet_3d( &
      species, config, geometry, state, temperature, local_ok)
    call require(local_ok, "transport initial state")
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), primitive, &
            recovered_temperature, sound_speed, local_ok)
          call require(local_ok, "transport primitive recovery")
          primitive(2:4) = [ &
            2.0_dp + real(i, dp), 3.0_dp + real(j, dp), &
            4.0_dp + real(k, dp)]
          if (geometry%cell_type(i, j, k) == eb_cut_cell_3d) then
            primitive(reactive_mass_fraction_component(1)) = &
              primitive(reactive_mass_fraction_component(1)) + 0.02_dp
            primitive(reactive_mass_fraction_component(size(species))) = &
              primitive(reactive_mass_fraction_component(size(species))) - &
              0.02_dp
          end if
          call reactive_primitive_to_conserved( &
            species, primitive, state(:, i, j, k), temperature(i, j, k), &
            sound_speed, local_ok)
          call require(local_ok, "transport state reconstruction")
        end do
      end do
    end do
    saved_state = state
    saved_temperature = temperature

    call reactive_eb_transport_timestep_3d( &
      species, transport, state, temperature, geometry, 0.35_dp, &
      .true., .true., .true., timestep, maximum_diffusivity, local_ok)
    call require(local_ok .and. timestep > 0.0_dp .and. &
      maximum_diffusivity > 0.0_dp, "transport timestep")
    call reactive_eb_transport_patch_timestep_3d(species,transport,state,temperature,geometry,.35_dp, &
      .true.,.true.,.true.,modified_timestep,modified_diffusivity,local_ok)
    call require(local_ok.and.modified_diffusivity==maximum_diffusivity,'plain patch transport estimate')
    ! kappa=.05: one normal opening, four tangential openings with area .05;
    ! no embedded-wall conductance. Half the open-face sum gives 12/h**2.
    call require(abs(modified_timestep-timestep/4._dp)<1.e-12_dp*timestep,'thin slab diffusion analytical rate')
    call reactive_eb_transport_patch_timestep_3d(species,transport,state,temperature,geometry,.35_dp, &
      .true.,.true.,.true.,modified_timestep,modified_diffusivity,local_ok,2._dp*maximum_diffusivity)
    call require(local_ok.and.modified_diffusivity==2._dp*maximum_diffusivity,'external diffusion bound')
    call require(abs(modified_timestep-timestep/8._dp)<1.e-12_dp*timestep,'external diffusion bound halves interval')
    call reactive_eb_transport_patch_timestep_3d(species,transport,state,temperature,geometry,.35_dp, &
      .false.,.false.,.false.,modified_timestep,modified_diffusivity,local_ok,maximum_diffusivity)
    call require(local_ok.and.modified_timestep==huge(1._dp).and.modified_diffusivity==0._dp, &
      'disabled transport imposes no bound')
    call reactive_eb_transport_patch_timestep_3d(species,transport,state,temperature,geometry,.35_dp, &
      .true.,.true.,.true.,modified_timestep,modified_diffusivity,local_ok,-1._dp)
    call require(.not.local_ok.and.modified_timestep==0._dp.and.modified_diffusivity==0._dp, &
      'invalid external diffusion bound fails with reset outputs')
    call overwrite_covered_storage(state, temperature, geometry)
    call reactive_eb_transport_timestep_3d( &
      species, transport, state, temperature, geometry, 0.35_dp, &
      .true., .true., .true., modified_timestep, modified_diffusivity, &
      local_ok)
    call require(local_ok .and. modified_timestep == timestep .and. &
      modified_diffusivity == maximum_diffusivity, &
      "covered storage excluded from transport timestep")
    call reactive_eb_transport_patch_timestep_3d(species,transport,state,temperature,geometry,.35_dp, &
      .true.,.true.,.true.,modified_timestep,modified_diffusivity,local_ok)
    call require(local_ok.and.modified_diffusivity==maximum_diffusivity,'covered patch storage ignored')
    call require(abs(modified_timestep-timestep/4._dp)<1.e-12_dp*timestep,'covered storage retains patch bound')
    state = saved_state
    temperature = saved_temperature
    interval = 0.2_dp * timestep
    allocate(bad_state(nvar - 1, 4, 4, 4), bad_temperature(4, 4, 3))
    bad_state = -71.0_dp
    bad_temperature = -72.0_dp
    call reactive_eb_transport_euler_update_3d( &
      species, transport, state, temperature, geometry, interval, .true., &
      .true., .true., .true., 0.5_dp, bad_state, bad_temperature, &
      minimum_theta, local_ok)
    call require(.not. local_ok, "invalid EB transport output shapes")

    call reactive_eb_transport_fluxes_rhs_3d( &
      species, transport, state, temperature, geometry, interval, &
      .true., .true., .true., .true., rhs, x_flux, y_flux, z_flux, &
      minimum_theta, local_ok)
    call require(local_ok, "transport face fluxes")
    call require_closed_transport_faces(geometry, x_flux, y_flux, z_flux)
    allocate(stage_rate,source=rhs);saved_state=state;saved_temperature=temperature
    call compute_sdc_eb_transport_3d(species,transport,state,temperature,geometry,interval, &
      [.true.,.true.,.true.,.true.],stage_rate,stage_flux,stage_theta,local_ok)
    call require(local_ok.and.stage_theta==minimum_theta,'SDC transport stage limiter')
    call require(all(stage_rate==rhs),'SDC transport rates exact')
    call require(all(stage_flux%x==interval*x_flux).and.all(stage_flux%y==interval*y_flux).and. &
      all(stage_flux%z==interval*z_flux),'SDC time integrated flux density exact')
    call require(all(same_real_bits(state,saved_state)).and.all(same_real_bits(temperature,saved_temperature)), &
      'SDC transport input unchanged')
    call reactive_eb_integrals_3d(stage_rate,geometry,rate_sum,local_ok)
    call require(local_ok,'SDC volume weighted rate integral')
    call reactive_eb_integrals_3d(abs(stage_rate),geometry,rate_scale,local_ok)
    call require(local_ok.and.all(abs(rate_sum)<=5e-12_dp*max(tiny(1._dp),rate_scale)), &
      'closed EB source conserves volume weighted inventory')
    saved_flux=stage_flux;bad_transport=transport;bad_transport(1)%diameter=-1
    call compute_sdc_eb_transport_3d(species,bad_transport,state,temperature,geometry,interval, &
      [.true.,.true.,.true.,.true.],stage_rate,stage_flux,stage_theta,local_ok)
    call require(.not.local_ok.and.stage_theta==minimum_theta,'bad SDC transport coefficient rejected')
    call require(all(stage_rate==rhs).and.all(stage_flux%x==saved_flux%x).and.all(stage_flux%y==saved_flux%y).and. &
      all(stage_flux%z==saved_flux%z),'SDC transport outputs rollback')

    call reactive_eb_integrals_3d( &
      state, geometry, initial_integrals, local_ok)
    call require(local_ok, "transport initial integrals")
    call reactive_eb_element_integrals_3d( &
      species, state, geometry, initial_elements, local_ok)
    call require(local_ok, "transport initial elements")
    call advance_reactive_eb_transport_3d( &
      species, transport, state, temperature, geometry, interval, &
      .true., .true., .true., .true., 0.5_dp, minimum_theta, local_ok)
    call require(local_ok, "transport SSPRK2 advance")
    call reactive_eb_integrals_3d( &
      state, geometry, final_integrals, local_ok)
    call require(local_ok, "transport final integrals")
    call reactive_eb_element_integrals_3d( &
      species, state, geometry, final_elements, local_ok)
    call require(local_ok, "transport final elements")
    conserved_error = maxval(abs(final_integrals - initial_integrals) / &
      max(1.0e-20_dp, abs(initial_integrals)))
    element_error = maxval(abs(final_elements - initial_elements) / &
      max(1.0e-30_dp, abs(initial_elements)))
    response = max(maxval(abs(state - saved_state)), &
      maxval(abs(temperature - saved_temperature)))
    call require(conserved_error <= 2.0e-11_dp, &
      "transport component conservation")
    call require(element_error <= 2.0e-11_dp, &
      "transport element conservation")
    call require(response > 1.0e-10_dp, "transport resolved response")
    call require(minimum_theta >= 0.0_dp .and. minimum_theta <= 1.0_dp, &
      "transport limiter range")
    call require_covered_identity( &
      state, temperature, saved_state, saved_temperature, geometry)

    state = saved_state
    temperature = saved_temperature
    call advance_reactive_eb_transport_3d( &
      species, transport, state, temperature, geometry, interval, &
      .false., .false., .false., .false., 0.5_dp, minimum_theta, local_ok)
    call require(local_ok .and. all(same_real_bits(state, saved_state)) .and. &
      all(same_real_bits(temperature, saved_temperature)), &
      "disabled transport bitwise identity")
    call advance_reactive_eb_transport_3d( &
      species, transport, state, temperature, geometry, -interval, &
      .true., .true., .true., .true., 0.5_dp, minimum_theta, local_ok)
    call require(.not. local_ok .and. &
      all(same_real_bits(state, saved_state)) .and. &
      all(same_real_bits(temperature, saved_temperature)), &
      "negative interval rollback")
    call advance_reactive_eb_transport_3d( &
      species, transport, state, temperature, geometry, interval, &
      .true., .true., .false., .true., 0.5_dp, minimum_theta, local_ok)
    call require(.not. local_ok .and. &
      all(same_real_bits(state, saved_state)) .and. &
      all(same_real_bits(temperature, saved_temperature)), &
      "invalid barodiffusion rollback")
    call advance_reactive_eb_transport_3d( &
      species, transport(1:size(transport) - 1), state, temperature, &
      geometry, interval, .true., .true., .true., .true., 0.5_dp, &
      minimum_theta, local_ok)
    call require(.not. local_ok .and. &
      all(same_real_bits(state, saved_state)) .and. &
      all(same_real_bits(temperature, saved_temperature)), &
      "invalid transport table rollback")

    call reactive_eb_transport_rhs_3d( &
      species, transport, state, temperature, geometry, &
      1000.0_dp * timestep, .false., .false., .true., .true., rhs, &
      minimum_theta, local_ok)
    call require(local_ok .and. minimum_theta < 1.0_dp .and. &
      minimum_theta >= 0.0_dp, "species outflow limiter activation")
    write(*, '(a,1x,a,5(a,es12.4))') &
      "transport axis", axis, ", dt=", timestep, &
      ", diffusivity=", maximum_diffusivity, ", theta=", minimum_theta, &
      ", conserved=", conserved_error, ", elements=", element_error
  end subroutine check_axis

  subroutine check_full_driver(species, reactions, transport)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)

    type(reactive_eb_3d_config) :: config
    type(eb_geometry_3d) :: geometry
    real(dp), allocatable :: state(:, :, :, :), saved_state(:, :, :, :)
    real(dp), allocatable :: temperature(:, :, :), saved_temperature(:, :, :)
    real(dp), allocatable :: legacy_state(:, :, :, :)
    real(dp), allocatable :: legacy_temperature(:, :, :)
    real(dp), allocatable :: disabled_state(:, :, :, :)
    real(dp), allocatable :: disabled_temperature(:, :, :)
    real(dp), allocatable :: transport_state(:, :, :, :)
    real(dp), allocatable :: transport_temperature(:, :, :)
    real(dp), allocatable :: coupled_state(:, :, :, :)
    real(dp), allocatable :: coupled_temperature(:, :, :)
    real(dp), allocatable :: bad_state(:, :, :, :), bad_temperature(:, :, :)
    real(dp), allocatable :: initial_integrals(:), transport_integrals(:)
    real(dp), allocatable :: coupled_integrals(:)
    real(dp) :: initial_elements(3), transport_elements(3), coupled_elements(3)
    real(dp) :: dt, minimum_theta, transport_error, conserved_error
    real(dp) :: element_error, transport_change, chemistry_change
    logical :: local_ok
    integer :: nvar

    config%nx = 4
    config%ny = 4
    config%nz = 4
    config%x_upper = 1.0_dp
    config%y_upper = 1.0_dp
    config%z_upper = 1.0_dp
    config%plane_axis = "x"
    config%plane_position = 0.4875_dp
    config%riemann_solver = "rusanov"
    config%redistribution = "state_redist"
    call build_axis_plane_eb_geometry_3d( &
      config%nx, config%ny, config%nz, config%x_lower, config%x_upper, &
      config%y_lower, config%y_upper, config%z_lower, config%z_upper, &
      config%plane_axis, config%plane_position, geometry, local_ok)
    call require(local_ok, "full driver geometry")
    nvar = reactive_nvar(size(species))
    allocate(state(nvar, 4, 4, 4), saved_state(nvar, 4, 4, 4))
    allocate(temperature(4, 4, 4), saved_temperature(4, 4, 4))
    allocate(legacy_state(nvar, 4, 4, 4), disabled_state(nvar, 4, 4, 4))
    allocate(legacy_temperature(4, 4, 4), disabled_temperature(4, 4, 4))
    allocate(transport_state(nvar, 4, 4, 4), coupled_state(nvar, 4, 4, 4))
    allocate(transport_temperature(4, 4, 4), coupled_temperature(4, 4, 4))
    allocate(bad_state(nvar - 1, 4, 4, 4), bad_temperature(4, 4, 3))
    allocate(initial_integrals(nvar), transport_integrals(nvar))
    allocate(coupled_integrals(nvar))
    call initialize_reactive_eb_density_sheet_3d( &
      species, config, geometry, state, temperature, local_ok)
    call require(local_ok, "full driver initial state")
    dt = 1.0e-7_dp

    bad_state = -81.0_dp
    bad_temperature = -82.0_dp
    call advance_reactive_eb_full_3d( &
      species, reactions, transport, state, temperature, geometry, &
      config%riemann_solver, config%redistribution, dt, .false., &
      2.0e-7_dp, 1.0e-12_dp, .false., .false., .false., .false., .false., &
      bad_state, bad_temperature, minimum_theta, local_ok, &
      config%state_redist_target_volume_fraction)
    call require(.not. local_ok .and. minimum_theta == 1.0_dp, &
      "invalid full driver output shapes")

    call advance_reactive_eb_strang_3d( &
      species, reactions, state, temperature, geometry, config%riemann_solver, &
      config%redistribution, dt, .false., 2.0e-7_dp, 1.0e-12_dp, &
      legacy_state, legacy_temperature, local_ok, &
      config%state_redist_target_volume_fraction)
    call require(local_ok, "legacy disabled split")
    call advance_reactive_eb_full_3d( &
      species, reactions, transport, state, temperature, geometry, &
      config%riemann_solver, config%redistribution, dt, .false., &
      2.0e-7_dp, 1.0e-12_dp, .false., .false., .false., .false., .false., &
      disabled_state, disabled_temperature, minimum_theta, local_ok, &
      config%state_redist_target_volume_fraction)
    call require(local_ok .and. minimum_theta == 1.0_dp, &
      "disabled full driver")
    call require(all(same_real_bits(legacy_state, disabled_state)) .and. &
      all(same_real_bits(legacy_temperature, disabled_temperature)), &
      "transport-disabled full driver parity")

    call apply_full_driver_profile( &
      species, state, temperature, geometry, config%initial_pressure)
    saved_state = state
    saved_temperature = temperature
    call reactive_eb_integrals_3d( &
      state, geometry, initial_integrals, local_ok)
    call require(local_ok, "full driver initial integrals")
    call reactive_eb_element_integrals_3d( &
      species, state, geometry, initial_elements, local_ok)
    call require(local_ok, "full driver initial elements")
    call advance_reactive_eb_full_3d( &
      species, reactions, transport, state, temperature, geometry, &
      config%riemann_solver, config%redistribution, dt, .false., &
      2.0e-7_dp, 1.0e-12_dp, .false., .false., .false., .false., .false., &
      disabled_state, disabled_temperature, minimum_theta, local_ok, &
      config%state_redist_target_volume_fraction)
    call require(local_ok .and. minimum_theta == 1.0_dp, &
      "profile disabled full driver")

    call advance_reactive_eb_full_3d( &
      species, reactions, transport, state, temperature, geometry, &
      config%riemann_solver, config%redistribution, dt, .false., &
      2.0e-7_dp, 1.0e-12_dp, .true., .true., .true., .true., .true., &
      transport_state, transport_temperature, minimum_theta, local_ok, &
      config%state_redist_target_volume_fraction)
    call require(local_ok, "transport-only full driver")
    call reactive_eb_integrals_3d( &
      transport_state, geometry, transport_integrals, local_ok)
    call require(local_ok, "transport-only integrals")
    call reactive_eb_element_integrals_3d( &
      species, transport_state, geometry, transport_elements, local_ok)
    call require(local_ok, "transport-only elements")
    transport_error = maxval(abs( &
      transport_integrals([irho, imy, imz, iet]) - &
      initial_integrals([irho, imy, imz, iet])) / &
      max(1.0e-20_dp, abs(initial_integrals([irho, imy, imz, iet]))))
    transport_error = max(transport_error, maxval(abs( &
      transport_integrals(reactive_species_component(1): &
        reactive_species_component(size(species))) - &
      initial_integrals(reactive_species_component(1): &
        reactive_species_component(size(species)))) / &
      max(1.0e-20_dp, abs(initial_integrals(reactive_species_component(1): &
        reactive_species_component(size(species)))))))
    call require(transport_error <= 2.0e-9_dp, &
      "transport-only full conservation")
    call require(maxval(abs(transport_elements - initial_elements) / &
      max(1.0e-30_dp, abs(initial_elements))) <= 2.0e-9_dp, &
      "transport-only full elements")
    transport_change = maxval(abs( &
      transport_state(reactive_species_component(1): &
        reactive_species_component(size(species)), :, :, :) - &
      disabled_state(reactive_species_component(1): &
        reactive_species_component(size(species)), :, :, :)))
    call require(transport_change > 1.0e-12_dp, &
      "transport-only species response")
    call require_covered_identity( &
      transport_state, transport_temperature, state, temperature, geometry)

    call advance_reactive_eb_full_3d( &
      species, reactions, transport, state, temperature, geometry, &
      config%riemann_solver, config%redistribution, dt, .true., &
      2.0e-7_dp, 1.0e-12_dp, .true., .true., .true., .true., .true., &
      coupled_state, coupled_temperature, minimum_theta, local_ok, &
      config%state_redist_target_volume_fraction)
    call require(local_ok, "transport-chemistry full driver")
    call reactive_eb_integrals_3d( &
      coupled_state, geometry, coupled_integrals, local_ok)
    call require(local_ok, "coupled full integrals")
    call reactive_eb_element_integrals_3d( &
      species, coupled_state, geometry, coupled_elements, local_ok)
    call require(local_ok, "coupled full elements")
    conserved_error = maxval(abs( &
      coupled_integrals([irho, imy, imz, iet]) - &
      initial_integrals([irho, imy, imz, iet])) / &
      max(1.0e-20_dp, abs(initial_integrals([irho, imy, imz, iet]))))
    element_error = maxval(abs(coupled_elements - initial_elements) / &
      max(1.0e-30_dp, abs(initial_elements)))
    chemistry_change = maxval(abs( &
      coupled_state(reactive_species_component(1): &
        reactive_species_component(size(species)), :, :, :) - &
      transport_state(reactive_species_component(1): &
        reactive_species_component(size(species)), :, :, :)))
    call require(conserved_error <= 5.0e-8_dp, &
      "coupled tangential conservation")
    call require(element_error <= 5.0e-8_dp, "coupled element conservation")
    call require(chemistry_change > 1.0e-12_dp, &
      "coupled chemistry species response")
    call require(minimum_theta >= 0.0_dp .and. minimum_theta <= 1.0_dp, &
      "coupled limiter range")
    call require_covered_identity( &
      coupled_state, coupled_temperature, state, temperature, geometry)

    disabled_state = -101.0_dp
    disabled_temperature = -102.0_dp
    call advance_reactive_eb_full_3d( &
      species, reactions, transport, state, temperature, geometry, "invalid", &
      config%redistribution, dt, .true., 2.0e-7_dp, 1.0e-12_dp, .true., &
      .true., .true., .true., .true., disabled_state, disabled_temperature, &
      minimum_theta, local_ok, config%state_redist_target_volume_fraction)
    call require(.not. local_ok .and. minimum_theta == 1.0_dp .and. &
      all(same_real_bits( &
      disabled_state, state)) .and. all(same_real_bits( &
      disabled_temperature, temperature)), "invalid solver rollback")
    disabled_state = -103.0_dp
    disabled_temperature = -104.0_dp
    call advance_reactive_eb_full_3d( &
      species, reactions, transport, state, temperature, geometry, &
      config%riemann_solver, "flux_redist", dt, .true., 2.0e-7_dp, &
      1.0e-12_dp, .true., .true., .true., .true., .true., disabled_state, &
      disabled_temperature, minimum_theta, local_ok, &
      config%state_redist_target_volume_fraction)
    call require(.not. local_ok .and. minimum_theta == 1.0_dp .and. &
      all(same_real_bits( &
      disabled_state, state)) .and. all(same_real_bits( &
      disabled_temperature, temperature)), "invalid redistribution rollback")
    disabled_state = -105.0_dp
    disabled_temperature = -106.0_dp
    call advance_reactive_eb_full_3d( &
      species, reactions, transport(1:size(transport) - 1), state, &
      temperature, geometry, config%riemann_solver, config%redistribution, dt, &
      .true., 2.0e-7_dp, 1.0e-12_dp, .true., .true., .true., .true., .true., &
      disabled_state, disabled_temperature, minimum_theta, local_ok, &
      config%state_redist_target_volume_fraction)
    call require(.not. local_ok .and. minimum_theta == 1.0_dp .and. &
      all(same_real_bits( &
      disabled_state, state)) .and. all(same_real_bits( &
      disabled_temperature, temperature)), "invalid transport rollback")
    write(*, '(a,es12.4,4(a,es12.4))') "full driver dt=", dt, &
      ", transport_error=", transport_error, ", conserved=", conserved_error, &
      ", elements=", element_error, ", chemistry_change=", chemistry_change
  end subroutine check_full_driver

  subroutine apply_full_driver_profile( &
      species, state, temperature, geometry, pressure)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: pressure

    real(dp), allocatable :: primitive(:)
    real(dp) :: recovered_temperature, sound_speed, delta
    logical :: local_ok
    integer :: i, j, k, nprim

    nprim = reactive_nprim(size(species))
    allocate(primitive(nprim))
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), primitive, &
            recovered_temperature, sound_speed, local_ok)
          call require(local_ok, "full driver primitive recovery")
          primitive(2:4) = 0.0_dp
          primitive(5) = pressure
          delta = 0.02_dp * real(i - 1, dp) / real(geometry%nx, dp)
          primitive(reactive_mass_fraction_component(1)) = &
            primitive(reactive_mass_fraction_component(1)) + delta
          primitive(reactive_mass_fraction_component(size(species))) = &
            primitive(reactive_mass_fraction_component(size(species))) - delta
          call reactive_primitive_to_conserved( &
            species, primitive, state(:, i, j, k), temperature(i, j, k), &
            sound_speed, local_ok)
          call require(local_ok, "full driver profile reconstruction")
        end do
      end do
    end do
  end subroutine apply_full_driver_profile

  subroutine require_closed_transport_faces(geometry, x_flux, y_flux, z_flux)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: x_flux(:, 0:, :, :)
    real(dp), intent(in) :: y_flux(:, :, 0:, :)
    real(dp), intent(in) :: z_flux(:, :, :, 0:)

    integer :: i, j, k

    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 0, geometry%nx
          if (i == 0 .or. i == geometry%nx .or. &
              geometry%x_face_fraction(i, j, k) <= 0.0_dp) then
            call require(all(same_real_bits( &
              x_flux(:, i, j, k), 0.0_dp)), "zero x transport face")
          end if
        end do
      end do
    end do
    do k = 1, geometry%nz
      do i = 1, geometry%nx
        do j = 0, geometry%ny
          if (j == 0 .or. j == geometry%ny .or. &
              geometry%y_face_fraction(i, j, k) <= 0.0_dp) then
            call require(all(same_real_bits( &
              y_flux(:, i, j, k), 0.0_dp)), "zero y transport face")
          end if
        end do
      end do
    end do
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        do k = 0, geometry%nz
          if (k == 0 .or. k == geometry%nz .or. &
              geometry%z_face_fraction(i, j, k) <= 0.0_dp) then
            call require(all(same_real_bits( &
              z_flux(:, i, j, k), 0.0_dp)), "zero z transport face")
          end if
        end do
      end do
    end do
  end subroutine require_closed_transport_faces

  subroutine overwrite_covered_storage(state, temperature, geometry)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry

    integer :: i, j, k

    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) /= eb_covered_cell_3d) cycle
          state(:, i, j, k) = -17.0_dp
          temperature(i, j, k) = -23.0_dp
        end do
      end do
    end do
  end subroutine overwrite_covered_storage

  subroutine require_covered_identity( &
      state, temperature, saved_state, saved_temperature, geometry)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(in) :: saved_state(:, :, :, :)
    real(dp), intent(in) :: saved_temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry

    integer :: i, j, k

    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) /= eb_covered_cell_3d) cycle
          call require(all(same_real_bits( &
            state(:, i, j, k), saved_state(:, i, j, k))) .and. &
            same_real_bits(temperature(i, j, k), &
              saved_temperature(i, j, k)), "covered transport identity")
        end do
      end do
    end do
  end subroutine require_covered_identity

  pure elemental logical function same_real_bits(left, right)
    real(dp), intent(in) :: left, right

    same_real_bits = transfer(left, 0_int64) == transfer(right, 0_int64)
  end function same_real_bits

  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) error stop label
  end subroutine require

end program test_reactive_eb_transport_3d
