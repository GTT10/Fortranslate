program test_reactive_eb_chemistry_3d
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use reactive_1d_mod, only: reactive_nvar, reactive_species_component
  use reactive_3d_mod, only: advance_reactive_chemistry_3d,react_sdc_state_3d
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d, eb_regular_cell_3d, &
    build_axis_plane_eb_geometry_3d
  use eb_reactive_hydro_3d_mod, only: &
    advance_reactive_eb_state_redistributed_euler_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config,build_configured_eb_geometry_3d
  use reactive_eb_3d_driver_mod, only: &
    initialize_reactive_eb_density_sheet_3d, &
    advance_reactive_eb_strang_3d, reactive_eb_integrals_3d, &
    reactive_eb_element_integrals_3d,react_sdc_eb_state_3d
  use reactive_eb_3d_driver_mod, only: advance_reactive_eb_sdc_3d
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  logical :: ok

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "elementary thermodynamics load")
  call load_h2o2_elementary_mechanism(reactions, ok)
  call require(ok, "elementary mechanism load")

  call check_mask_axis(species, reactions, "x")
  call check_mask_axis(species, reactions, "y")
  call check_mask_axis(species, reactions, "z")
  call check_eb_split(species, reactions)
  call check_sdc()
  call check_sdc_loop()
  write(*, '(a)') "test_reactive_eb_chemistry_3d: PASS"

contains
  subroutine check_sdc_loop()
    type(reactive_eb_3d_config) :: config
    type(eb_geometry_3d) :: geometry,extended
    type(gas_transport_species), allocatable :: transport(:)
    type(eb_integrated_flux_3d) :: flux,saved_flux
    real(dp) :: u(12,2,2,2),t(2,2,2),ir(8,2,2,2),initial(12),final(12),scale(12)
    real(dp) :: boundary(12),wall(12),external(12),theta,saved(12,2,2,2),st(2,2,2),si(8,2,2,2)
    integer :: mode,i,j,k
    character(len=12), parameter :: methods(3)=[character(len=12)::'none','flux_redist','state_redist']
    character(len=7) :: faces(6)
    config%nx=2;config%ny=2;config%nz=2;config%plane_position=.37_dp
    call build_configured_eb_geometry_3d(config,geometry,ok)
    call require(ok,'coupled SDC geometry')
    call build_configured_eb_geometry_3d(config,extended,ok,halo=6)
    call require(ok,'coupled SDC extended geometry')
    call load_h2o2_elementary_transport(transport,ok)
    call require(ok,'coupled SDC transport')
    faces='outflow'
    do mode=1,3
      call initialize_reactive_eb_density_sheet_3d(species,config,geometry,u,t,ok)
      call require(ok,'coupled SDC initial state')
      call reactive_eb_integrals_3d(u,geometry,initial,ok)
      call require(ok,'coupled SDC initial inventory')
      call reactive_eb_integrals_3d(abs(u),geometry,scale,ok)
      ir=0;boundary=-1;wall=-2;external=-3;theta=-4
      call advance_reactive_eb_sdc_3d(species,reactions,transport,u,t,geometry,extended,1e-9_dp,2,.true., &
        [.true.,.true.,.true.,.true.],[20000._dp,-30000._dp,40000._dp],1e-6_dp,1e-12_dp,trim(methods(mode)), &
        faces,boundary,wall,external,theta,ok,ir,flux)
      call require(ok,'coupled reacting transporting forced EB SDC')
      call reactive_eb_integrals_3d(u,geometry,final,ok)
      call require(ok.and.all(abs(final(:5)-initial(:5)-boundary(:5)-wall(:5)-external(:5)) &
        <=2e-10_dp*max(1._dp,scale(:5))),'coupled EB SDC Euler inventory closure')
      do k=1,2
      do j=1,2
      do i=1,2
        if(geometry%cell_type(i,j,k)==eb_covered_cell_3d) cycle
        call require(abs(sum(u(6:,i,j,k))-u(1,i,j,k))<2e-12_dp,'coupled EB SDC species closure')
      end do
      end do
      end do
      saved=u;st=t;si=ir;saved_flux=flux
      boundary=-1;wall=-2;external=-3;theta=-4
      call advance_reactive_eb_sdc_3d(species,reactions,transport,u,t,geometry,extended,1e-9_dp,2,.true., &
        [.true.,.true.,.true.,.true.],[20000._dp,-30000._dp,40000._dp],1e-6_dp,1e-12_dp,trim(methods(mode)), &
        faces,boundary,wall,external,theta,ok,ir,flux,difmag=-1._dp)
      call require(.not.ok.and.all(u==saved).and.all(t==st).and.all(ir==si),'coupled EB SDC state rollback')
      call require(all(boundary==-1).and.all(wall==-2).and.all(external==-3).and.theta==-4,'coupled SDC ledger rollback')
      call require(all(flux%x==saved_flux%x).and.all(flux%y==saved_flux%y).and.all(flux%z==saved_flux%z), &
        'coupled SDC integrated flux rollback')
    end do
  end subroutine check_sdc_loop


  subroutine check_sdc()
    type(reactive_eb_3d_config) :: config
    type(eb_geometry_3d) :: geometry
    real(dp), allocatable :: old(:,:,:,:),t0(:,:,:),rates(:,:,:,:),u(:,:,:,:),t(:,:,:),ir(:,:,:,:)
    real(dp), allocatable :: reference(:,:,:,:),rt(:,:,:),ri(:,:,:,:),saved(:,:,:,:),st(:,:,:),si(:,:,:,:)
    integer :: nv,i
    config%nx=4;config%ny=1;config%nz=1
    call build_axis_plane_eb_geometry_3d(4,1,1,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',.37_dp,geometry,ok)
    call require(ok,'SDC EB geometry')
    nv=reactive_nvar(size(species))
    allocate(old(nv,4,1,1),t0(4,1,1),rates(nv,4,1,1),ir(nv-4,4,1,1))
    call initialize_reactive_eb_density_sheet_3d(species,config,geometry,old,t0,ok)
    call require(ok,'SDC EB initial state')
    reference=old;rt=t0
    do i=1,4
      if(geometry%cell_type(i,1,1)/=eb_covered_cell_3d) cycle
      reference(:,i,1,1)=old(:,2,1,1);rt(i,1,1)=t0(2,1,1)
    end do
    rates=0;rates(2,:,:,:)=20000;rates(5,:,:,:)=30000
    u=reference;t=rt;ir=-7;ri=ir
    call react_sdc_state_3d(species,reactions,reference,rt,rates,1e-9_dp,1e-6_dp,1e-12_dp,u,t,ri,ok)
    call require(ok,'regular forced reaction reference')
    reference=u;rt=t
    do i=1,4
      if(geometry%cell_type(i,1,1)/=eb_covered_cell_3d) cycle
      old(:,i,1,1)=ieee_value(0._dp,ieee_quiet_nan)
      t0(i,1,1)=-99;rates(:,i,1,1)=ieee_value(0._dp,ieee_quiet_nan)
    end do
    u=old;t=t0
    call react_sdc_eb_state_3d(species,reactions,old,t0,geometry,rates,1e-9_dp,1e-6_dp,1e-12_dp,u,t,ir,ok)
    call require(ok,'EB forced reaction accepts covered sentinel')
    do i=1,4
      if(geometry%cell_type(i,1,1)==eb_covered_cell_3d) then
        call require(all(same_real_bits(u(:,i,1,1),old(:,i,1,1))).and.t(i,1,1)==t0(i,1,1), &
          'covered state and temperature preserved')
        call require(all(ir(:,i,1,1)==-7),'covered reaction estimate preserved')
      else
        call require(all(u(:,i,1,1)==reference(:,i,1,1)).and.t(i,1,1)==rt(i,1,1),'active forced state exact')
        call require(all(ir(:,i,1,1)==ri(:,i,1,1)),'active forced reaction estimate exact')
      end if
    end do
    saved=u;st=t;si=ir;rates(1,4,1,1)=-1e30_dp
    call react_sdc_eb_state_3d(species,reactions,old,t0,geometry,rates,1e-9_dp,1e-6_dp,1e-12_dp,u,t,ir,ok)
    call require(.not.ok,'late active cell rejects EB reaction')
    call require(all(same_real_bits(u,saved)).and.all(same_real_bits(t,st)).and.all(same_real_bits(ir,si)), &
      'failed EB reaction retains all outputs')
  end subroutine check_sdc

  subroutine check_mask_axis(species, reactions, axis)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    character(len=*), intent(in) :: axis

    type(reactive_eb_3d_config) :: config
    type(eb_geometry_3d) :: geometry
    real(dp), allocatable :: state(:, :, :, :), saved_state(:, :, :, :)
    real(dp), allocatable :: temperature(:, :, :), saved_temperature(:, :, :)
    real(dp), allocatable :: before_integrals(:), after_integrals(:)
    real(dp) :: before_elements(3), after_elements(3)
    real(dp) :: conserved_error, element_error, reaction_change
    logical, allocatable :: active_mask(:, :, :)
    logical :: local_ok, found_disabled
    integer :: i, j, k, nvar, disabled_i, disabled_j, disabled_k

    config%nx = 4
    config%ny = 4
    config%nz = 4
    config%plane_axis = axis
    config%plane_position = 0.37_dp
    call build_axis_plane_eb_geometry_3d( &
      config%nx, config%ny, config%nz, config%x_lower, config%x_upper, &
      config%y_lower, config%y_upper, config%z_lower, config%z_upper, &
      axis, config%plane_position, geometry, local_ok)
    call require(local_ok, "mask geometry")

    nvar = reactive_nvar(size(species))
    allocate(state(nvar, config%nx, config%ny, config%nz))
    allocate(saved_state(nvar, config%nx, config%ny, config%nz))
    allocate(temperature(config%nx, config%ny, config%nz))
    allocate(saved_temperature(config%nx, config%ny, config%nz))
    allocate(active_mask(config%nx, config%ny, config%nz))
    allocate(before_integrals(nvar), after_integrals(nvar))
    call initialize_reactive_eb_density_sheet_3d( &
      species, config, geometry, state, temperature, local_ok)
    call require(local_ok, "mask state")
    saved_state = state
    saved_temperature = temperature
    active_mask = geometry%cell_type /= eb_covered_cell_3d
    found_disabled = .false.
    disabled_i = 0
    disabled_j = 0
    disabled_k = 0
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (.not. found_disabled .and. &
              geometry%cell_type(i, j, k) == eb_regular_cell_3d) then
            active_mask(i, j, k) = .false.
            disabled_i = i
            disabled_j = j
            disabled_k = k
            found_disabled = .true.
          end if
        end do
      end do
    end do
    call require(found_disabled, "disabled active cell")
    call reactive_eb_integrals_3d( &
      state, geometry, before_integrals, local_ok)
    call require(local_ok, "mask conserved integrals")
    call reactive_eb_element_integrals_3d( &
      species, state, geometry, before_elements, local_ok)
    call require(local_ok, "mask element integrals")

    call advance_reactive_chemistry_3d( &
      species, reactions, state, temperature, config%nx, config%ny, config%nz, &
      5.0e-7_dp, 2.0e-7_dp, 1.0e-12_dp, local_ok, active_mask)
    call require(local_ok, "active-mask chemistry")
    call reactive_eb_integrals_3d( &
      state, geometry, after_integrals, local_ok)
    call require(local_ok, "mask final conserved integrals")
    call reactive_eb_element_integrals_3d( &
      species, state, geometry, after_elements, local_ok)
    call require(local_ok, "mask final element integrals")

    conserved_error = maxval(abs(after_integrals([irho, imx, imy, imz, iet]) - &
      before_integrals([irho, imx, imy, imz, iet])) / &
      max(1.0_dp, abs(before_integrals([irho, imx, imy, imz, iet]))))
    element_error = maxval(abs(after_elements - before_elements) / &
      max(1.0e-30_dp, abs(before_elements)))
    reaction_change = maxval(abs( &
      state(reactive_species_component(1): &
        reactive_species_component(size(species)), :, :, :) - &
      saved_state(reactive_species_component(1): &
        reactive_species_component(size(species)), :, :, :)))
    call require(conserved_error <= 3.0e-12_dp, "mask conserved quantities")
    call require(element_error <= 3.0e-10_dp, "mask element conservation")
    call require(reaction_change > 1.0e-12_dp, "mask chemistry change")
    call require(all(same_real_bits( &
      state(:, disabled_i, disabled_j, disabled_k), &
      saved_state(:, disabled_i, disabled_j, disabled_k))) .and. &
      same_real_bits(temperature(disabled_i, disabled_j, disabled_k), &
        saved_temperature(disabled_i, disabled_j, disabled_k)), &
      "disabled cell identity")
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) /= eb_covered_cell_3d) cycle
          call require(all(same_real_bits( &
            state(:, i, j, k), saved_state(:, i, j, k))) .and. &
            same_real_bits(temperature(i, j, k), &
              saved_temperature(i, j, k)), &
            "covered cell identity")
        end do
      end do
    end do

    state = saved_state
    temperature = saved_temperature
    active_mask = .true.
    call advance_reactive_chemistry_3d( &
      species, reactions, state, temperature, config%nx, config%ny, config%nz, &
      5.0e-7_dp, 2.0e-7_dp, 1.0e-12_dp, local_ok, active_mask(1:3, :, :))
    call require(.not. local_ok .and. &
      all(same_real_bits(state, saved_state)) .and. &
      all(same_real_bits(temperature, saved_temperature)), &
      "invalid mask rollback")

    state = saved_state
    temperature = saved_temperature
    active_mask = geometry%cell_type /= eb_covered_cell_3d
    call require(active_mask(config%nx, config%ny, config%nz), &
      "late-plane failure cell is active")
    state(irho, config%nx, config%ny, config%nz) = &
      -abs(state(irho, config%nx, config%ny, config%nz))
    saved_state = state
    saved_temperature = temperature
    call advance_reactive_chemistry_3d( &
      species, reactions, state, temperature, config%nx, config%ny, config%nz, &
      5.0e-7_dp, 2.0e-7_dp, 1.0e-12_dp, local_ok, active_mask)
    call require(.not. local_ok .and. &
      all(same_real_bits(state, saved_state)) .and. &
      all(same_real_bits(temperature, saved_temperature)), &
      "late-plane chemistry rollback")
    write(*, '(a,1x,a,3(a,es12.4))') &
      "mask axis", axis, ", conserved=", conserved_error, &
      ", elements=", element_error, ", change=", reaction_change
  end subroutine check_mask_axis

  subroutine check_eb_split(species, reactions)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)

    type(reactive_eb_3d_config) :: config
    type(eb_geometry_3d) :: geometry
    real(dp), allocatable :: state(:, :, :, :), inert_state(:, :, :, :)
    real(dp), allocatable :: reactive_state(:, :, :, :), saved_state(:, :, :, :)
    real(dp), allocatable :: temperature(:, :, :), inert_temperature(:, :, :)
    real(dp), allocatable :: reactive_temperature(:, :, :), saved_temperature(:, :, :)
    real(dp), allocatable :: initial_integrals(:), final_integrals(:)
    real(dp), allocatable :: direct_state(:, :, :, :), direct_temperature(:, :, :)
    real(dp) :: initial_elements(3), final_elements(3)
    real(dp) :: conserved_error, element_error, chemistry_change
    logical :: local_ok
    integer :: nvar, i, j, k

    config%nx = 4
    config%ny = 4
    config%nz = 4
    config%plane_axis = "x"
    config%plane_position = 0.4875_dp
    config%redistribution = "state_redist"
    call build_axis_plane_eb_geometry_3d( &
      config%nx, config%ny, config%nz, config%x_lower, config%x_upper, &
      config%y_lower, config%y_upper, config%z_lower, config%z_upper, &
      config%plane_axis, config%plane_position, geometry, local_ok)
    call require(local_ok, "split geometry")
    call require(minval(geometry%volume_fraction, &
      mask=geometry%volume_fraction > 0.0_dp) < &
      config%state_redist_target_volume_fraction, &
      "split activates StateRedist merge")
    nvar = reactive_nvar(size(species))
    allocate(state(nvar, 4, 4, 4), inert_state(nvar, 4, 4, 4))
    allocate(reactive_state(nvar, 4, 4, 4), saved_state(nvar, 4, 4, 4))
    allocate(temperature(4, 4, 4), inert_temperature(4, 4, 4))
    allocate(reactive_temperature(4, 4, 4), saved_temperature(4, 4, 4))
    allocate(direct_state(nvar, 4, 4, 4), direct_temperature(4, 4, 4))
    allocate(initial_integrals(nvar), final_integrals(nvar))
    call initialize_reactive_eb_density_sheet_3d( &
      species, config, geometry, state, temperature, local_ok)
    call require(local_ok, "split state")
    call advance_reactive_eb_state_redistributed_euler_3d( &
      species, state, temperature, geometry, "rusanov", 1.0e-7_dp, &
      direct_state, direct_temperature, local_ok, &
      config%state_redist_target_volume_fraction)
    call require(local_ok, "direct inert hydro")
    call advance_reactive_eb_strang_3d( &
      species, reactions, state, temperature, geometry, "rusanov", &
      config%redistribution, 1.0e-7_dp, .false., 2.0e-7_dp, 1.0e-12_dp, &
      inert_state, inert_temperature, local_ok, &
      config%state_redist_target_volume_fraction)
    call require(local_ok, "inert EB split")
    call require(all(same_real_bits(inert_state, direct_state)) .and. &
      all(same_real_bits(inert_temperature, direct_temperature)), &
      "inert hydro parity")

    call reactive_eb_integrals_3d( &
      state, geometry, initial_integrals, local_ok)
    call require(local_ok, "split initial integrals")
    call reactive_eb_element_integrals_3d( &
      species, state, geometry, initial_elements, local_ok)
    call require(local_ok, "split initial elements")
    call advance_reactive_eb_strang_3d( &
      species, reactions, state, temperature, geometry, "rusanov", &
      config%redistribution, 1.0e-7_dp, .true., 2.0e-7_dp, 1.0e-12_dp, &
      reactive_state, reactive_temperature, local_ok, &
      config%state_redist_target_volume_fraction)
    call require(local_ok, "reactive EB split")
    call reactive_eb_integrals_3d( &
      reactive_state, geometry, final_integrals, local_ok)
    call require(local_ok, "split final integrals")
    call reactive_eb_element_integrals_3d( &
      species, reactive_state, geometry, final_elements, local_ok)
    call require(local_ok, "split final elements")
    conserved_error = maxval(abs(final_integrals([irho, imy, imz, iet]) - &
      initial_integrals([irho, imy, imz, iet])) / &
      max(1.0_dp, abs(initial_integrals([irho, imy, imz, iet]))))
    element_error = maxval(abs(final_elements - initial_elements) / &
      max(1.0e-30_dp, abs(initial_elements)))
    chemistry_change = maxval(abs(reactive_state(reactive_species_component(1): &
      reactive_species_component(size(species)), :, :, :) - &
      inert_state(reactive_species_component(1): &
        reactive_species_component(size(species)), :, :, :)))
    call require(conserved_error <= 5.0e-10_dp, "split conserved quantities")
    call require(element_error <= 5.0e-9_dp, "split element conservation")
    call require(chemistry_change > 1.0e-12_dp, "split chemistry effect")
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) /= eb_covered_cell_3d) cycle
          call require(all(same_real_bits( &
            reactive_state(:, i, j, k), state(:, i, j, k))) .and. &
            same_real_bits(reactive_temperature(i, j, k), &
              temperature(i, j, k)), &
            "split covered identity")
        end do
      end do
    end do

    saved_state = state
    saved_temperature = temperature
    call advance_reactive_eb_strang_3d( &
      species, reactions, state, temperature, geometry, "invalid", &
      config%redistribution, 1.0e-7_dp, .true., 2.0e-7_dp, 1.0e-12_dp, &
      inert_state, inert_temperature, local_ok, &
      config%state_redist_target_volume_fraction)
    call require(.not. local_ok .and. &
      all(same_real_bits(inert_state, saved_state)) .and. &
      all(same_real_bits(inert_temperature, saved_temperature)), &
      "split invalid solver rollback")
    call advance_reactive_eb_strang_3d( &
      species, reactions, state, temperature, geometry, "rusanov", &
      config%redistribution, 1.0e-7_dp, .true., -1.0_dp, 1.0e-12_dp, &
      inert_state, inert_temperature, local_ok, &
      config%state_redist_target_volume_fraction)
    call require(.not. local_ok .and. &
      all(same_real_bits(inert_state, saved_state)) .and. &
      all(same_real_bits(inert_temperature, saved_temperature)), &
      "split invalid tolerance rollback")
    call advance_reactive_eb_strang_3d( &
      species, reactions(1:0), state, temperature, geometry, "rusanov", &
      config%redistribution, 1.0e-7_dp, .true., 2.0e-7_dp, 1.0e-12_dp, &
      inert_state, inert_temperature, local_ok, &
      config%state_redist_target_volume_fraction)
    call require(.not. local_ok .and. &
      all(same_real_bits(inert_state, saved_state)) .and. &
      all(same_real_bits(inert_temperature, saved_temperature)), &
      "split empty mechanism rollback")
    write(*, '(3(a,es12.4))') &
      "split conserved=", conserved_error, ", elements=", element_error, &
      ", chemistry_change=", chemistry_change
  end subroutine check_eb_split

  pure elemental logical function same_real_bits(left, right)
    real(dp), intent(in) :: left, right

    same_real_bits = transfer(left, 0_int64) == transfer(right, 0_int64)
  end function same_real_bits

  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) error stop label
  end subroutine require

end program test_reactive_eb_chemistry_3d
