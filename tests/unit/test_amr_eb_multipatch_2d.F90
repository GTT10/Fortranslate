program test_amr_eb_multipatch_2d
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use precision_mod, only: dp
  use state_indices_mod, only: irho, iet
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_mass_properties
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_mass_fraction_component, reactive_primitive_to_conserved
  use reactive_boundary_2d_mod, only: &
    reactive_boundary_set_2d, build_reactive_boundary_set_2d
  use simulation_config_reactive_2d_mod, only: reactive_2d_config
  use eb_geometry_2d_mod, only: &
    eb_geometry_2d, eb_covered_cell, build_eb_geometry_2d
  use eb_reactive_transport_2d_mod, only: &
    reactive_eb_transport_timestep_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, build_amr_eb_patch_2d
  use amr_eb_regrid_2d_mod, only: &
    amr_eb_tagging_criteria_2d, amr_eb_regrid_plan_collection_2d, &
    reactive_eb_patch_set_2d, build_amr_eb_regrid_plan_collection_2d, &
    initialize_reactive_eb_patch_set_2d, &
    average_down_reactive_eb_patch_set_2d, &
    composite_reactive_eb_patch_set_integral_2d, &
    regrid_reactive_eb_patch_set_2d, &
    advance_reactive_eb_patch_set_hydro_2d
  use reactive_eb_amr_2d_driver_mod, only: &
    advance_reactive_eb_patch_set_strang_2d
  use amr_eb_multipatch_transport_2d_mod, only: &
    advance_reactive_eb_patch_set_transport_2d
  implicit none

  integer, parameter :: coarse_nx = 14, coarse_ny = 14, ratio = 2
  type(eb_geometry_2d) :: coarse_geometry
  type(eb_geometry_2d), allocatable :: old_geometries(:)
  type(eb_geometry_2d), allocatable :: new_geometries(:)
  type(amr_eb_patch_2d) :: geometry_patch
  type(amr_eb_tagging_criteria_2d) :: criteria
  type(amr_eb_regrid_plan_collection_2d) :: old_collection
  type(amr_eb_regrid_plan_collection_2d) :: new_collection
  type(reactive_eb_patch_set_2d) :: old_set, new_set, removed_set
  type(reactive_eb_patch_set_2d) :: hydro_set, failed_set
  type(reactive_eb_patch_set_2d) :: chemistry_set, chemistry_failed_set
  type(reactive_eb_patch_set_2d) :: transport_set, transport_new_set
  type(reactive_eb_patch_set_2d) :: transport_failed_set
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(reactive_2d_config) :: transport_config
  type(reactive_boundary_set_2d) :: boundaries
  real(dp) :: coarse_level_set(0:coarse_nx, 0:coarse_ny)
  real(dp), allocatable :: primitive(:), mass_fractions(:), state_cell(:)
  real(dp), allocatable :: coarse_state(:, :, :), coarse_temperature(:, :)
  real(dp), allocatable :: averaged_state(:, :, :)
  real(dp), allocatable :: averaged_temperature(:, :)
  real(dp), allocatable :: new_coarse_state(:, :, :)
  real(dp), allocatable :: new_coarse_temperature(:, :)
  real(dp), allocatable :: chemistry_coarse_state(:, :, :)
  real(dp), allocatable :: chemistry_coarse_temperature(:, :)
  real(dp), allocatable :: transport_coarse_state(:, :, :)
  real(dp), allocatable :: transport_coarse_temperature(:, :)
  real(dp), allocatable :: transport_new_coarse_state(:, :, :)
  real(dp), allocatable :: transport_new_coarse_temperature(:, :)
  real(dp), allocatable :: integral_before(:), integral_after(:)
  real(dp) :: mole_fractions(7), x, y, temperature_cell, sound_speed
  real(dp) :: dt, factor, integral_scale, reaction_change, state_scale
  real(dp) :: coarse_transport_dt, fine_transport_dt, transport_interval
  real(dp) :: maximum_diffusivity, initial_span, final_span, minimum_theta
  character(len=64) :: hydro_failure_context
  logical :: old_tags(coarse_nx, coarse_ny)
  logical :: new_tags(coarse_nx, coarse_ny), ok
  integer :: child, component, i, j, k, nvar

  do j = 0, coarse_ny
    y = real(j, dp) / real(coarse_ny, dp)
    do i = 0, coarse_nx
      x = real(i, dp) / real(coarse_nx, dp)
      coarse_level_set(i, j) = x + y - 0.78_dp
    end do
  end do
  call build_eb_geometry_2d( &
    coarse_level_set, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    coarse_geometry, ok)
  call require(ok, "multipatch coarse EB geometry")

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "multipatch thermodynamic database")
  call load_h2o2_elementary_transport(transport, ok)
  call require(ok, "multipatch transport database")
  call load_h2o2_elementary_mechanism(reactions, ok)
  call require(ok, "multipatch elementary kinetics")
  nvar = reactive_nvar(size(species))
  allocate(primitive(reactive_nprim(size(species))))
  allocate(mass_fractions(size(species)), state_cell(nvar))
  mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.55643_dp]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  call require(ok, "multipatch composition conversion")
  primitive(1:5) = [0.31_dp, 0.0_dp, 0.0_dp, 0.0_dp, 135000.0_dp]
  do i = 1, size(species)
    primitive(reactive_mass_fraction_component(i)) = mass_fractions(i)
  end do
  call reactive_primitive_to_conserved( &
    species, primitive, state_cell, temperature_cell, sound_speed, ok)
  call require(ok, "multipatch reference state")

  allocate(coarse_state(nvar, coarse_nx, coarse_ny))
  allocate(coarse_temperature(coarse_nx, coarse_ny))
  allocate(averaged_state(nvar, coarse_nx, coarse_ny))
  allocate(averaged_temperature(coarse_nx, coarse_ny))
  allocate(new_coarse_state(nvar, coarse_nx, coarse_ny))
  allocate(new_coarse_temperature(coarse_nx, coarse_ny))
  allocate(chemistry_coarse_state(nvar, coarse_nx, coarse_ny))
  allocate(chemistry_coarse_temperature(coarse_nx, coarse_ny))
  allocate(integral_before(nvar), integral_after(nvar))
  coarse_state = spread(spread(state_cell, 2, coarse_nx), 3, coarse_ny)
  coarse_temperature = temperature_cell

  criteria%buffer_cells = 0
  criteria%minimum_patch_cells_x = 5
  criteria%minimum_patch_cells_y = 5
  criteria%maximum_patch_gap_cells = 0
  old_tags = .false.
  old_tags(1:5, 6:10) = .true.
  old_tags(9:13, 9:13) = .true.
  call build_amr_eb_regrid_plan_collection_2d( &
    old_tags, criteria, old_collection, ok)
  call require(ok .and. old_collection%patch_count() == 2, &
    "old two-patch plan")
  allocate(old_geometries(old_collection%patch_count()))
  do child = 1, old_collection%patch_count()
    call build_patch_geometry( &
      coarse_geometry, old_collection%plans(child)%coarse_i_lower, &
      old_collection%plans(child)%coarse_i_upper, &
      old_collection%plans(child)%coarse_j_lower, &
      old_collection%plans(child)%coarse_j_upper, ratio, &
      old_geometries(child), geometry_patch, ok)
    call require(ok, "old child EB geometry")
  end do
  call initialize_reactive_eb_patch_set_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, &
    old_geometries, old_collection, ratio, old_set, ok)
  call require(ok .and. old_set%is_valid(coarse_geometry, nvar) .and. &
    old_set%patch_count() == 2, "initialize reactive EB patch set")
  old_set%children(1)%state = 1.01_dp * old_set%children(1)%state
  old_set%children(2)%state = 0.99_dp * old_set%children(2)%state
  call require(old_set%is_valid(coarse_geometry, nvar), &
    "nonmatching multipatch hydro state")

  allocate(transport_coarse_state, source=coarse_state)
  allocate(transport_coarse_temperature, source=coarse_temperature)
  allocate(transport_new_coarse_state, mold=coarse_state)
  allocate(transport_new_coarse_temperature, mold=coarse_temperature)
  transport_set = old_set
  call impose_double_hotspot( &
    species, mass_fractions, coarse_geometry, transport_coarse_state, &
    transport_coarse_temperature, ok)
  call require(ok, "multipatch coarse transport hotspot")
  do child = 1, transport_set%patch_count()
    call impose_double_hotspot( &
      species, mass_fractions, transport_set%children(child)%geometry, &
      transport_set%children(child)%state, &
      transport_set%children(child)%temperature, ok)
    call require(ok, "multipatch fine transport hotspot")
  end do
  transport_config%initial_temperature = temperature_cell
  transport_config%initial_pressure = 135000.0_dp
  transport_config%initial_velocity_x = 0.0_dp
  transport_config%initial_velocity_y = 0.0_dp
  transport_config%boundary_x_lower = "outflow"
  transport_config%boundary_x_upper = "outflow"
  transport_config%boundary_y_lower = "outflow"
  transport_config%boundary_y_upper = "outflow"
  transport_config%transport_enabled = .true.
  transport_config%viscosity_enabled = .false.
  transport_config%thermal_conduction_enabled = .true.
  transport_config%species_diffusion_enabled = .false.
  transport_config%barodiffusion_enabled = .false.
  call build_reactive_boundary_set_2d( &
    species, transport_config, boundaries, ok)
  call require(ok, "multipatch transport boundaries")
  call composite_reactive_eb_patch_set_integral_2d( &
    transport_coarse_state, coarse_geometry, transport_set, &
    integral_before, ok)
  call require(ok, "initial multipatch transport integral")
  initial_span = patch_set_temperature_span( &
    transport_coarse_temperature, coarse_geometry, transport_set)
  call reactive_eb_transport_timestep_2d( &
    species, transport, transport_coarse_state, &
    transport_coarse_temperature, coarse_geometry, 0.30_dp, .false., &
    .true., .false., coarse_transport_dt, maximum_diffusivity, ok)
  call require(ok .and. coarse_transport_dt > 0.0_dp, &
    "multipatch coarse transport timestep")
  transport_interval = 0.25_dp * coarse_transport_dt
  do child = 1, transport_set%patch_count()
    call reactive_eb_transport_timestep_2d( &
      species, transport, transport_set%children(child)%state, &
      transport_set%children(child)%temperature, &
      transport_set%children(child)%geometry, 0.30_dp, .false., .true., &
      .false., fine_transport_dt, maximum_diffusivity, ok)
    call require(ok .and. fine_transport_dt > 0.0_dp, &
      "multipatch fine transport timestep")
    transport_interval = min(transport_interval, 0.25_dp * &
      real(transport_set%children(child)%patch%refinement_ratio, dp) * &
      fine_transport_dt)
  end do
  call advance_reactive_eb_patch_set_transport_2d( &
    species, transport, transport_coarse_state, &
    transport_coarse_temperature, coarse_geometry, transport_set, &
    transport_interval, .false., .true., .false., .false., boundaries, &
    0.5_dp, 2, transport_new_coarse_state, &
    transport_new_coarse_temperature, transport_new_set, minimum_theta, &
    ok, failure_context=hydro_failure_context)
  if (.not. ok) write(*, '(a)') trim(hydro_failure_context)
  call require(ok .and. minimum_theta > 0.999999999_dp .and. &
    transport_new_set%is_valid(coarse_geometry, nvar), &
    "subcycled multipatch EB transport")
  call composite_reactive_eb_patch_set_integral_2d( &
    transport_new_coarse_state, coarse_geometry, transport_new_set, &
    integral_after, ok)
  integral_scale = max(1.0_dp, maxval(abs(integral_before)))
  call require(ok .and. maxval(abs(integral_after - integral_before)) <= &
    5.0e-10_dp * integral_scale, "multipatch transport conservation")
  final_span = patch_set_temperature_span( &
    transport_new_coarse_temperature, coarse_geometry, transport_new_set)
  call require(final_span < initial_span, &
    "multipatch conduction reduces hierarchy temperature span")
  call assert_patch_set_covered_unchanged( &
    transport_coarse_state, transport_set, transport_new_coarse_state, &
    transport_new_set, coarse_geometry, "multipatch covered state unchanged")

  call advance_reactive_eb_patch_set_transport_2d( &
    species, transport, transport_coarse_state, &
    transport_coarse_temperature, coarse_geometry, transport_set, &
    -transport_interval, .false., .true., .false., .false., boundaries, &
    0.5_dp, 2, transport_new_coarse_state, &
    transport_new_coarse_temperature, transport_failed_set, minimum_theta, &
    ok)
  call require(.not. ok .and. &
    all(transport_new_coarse_state == transport_coarse_state) .and. &
    all(transport_new_coarse_temperature == transport_coarse_temperature), &
    "invalid multipatch transport coarse rollback")
  do child = 1, transport_set%patch_count()
    call require(all(transport_failed_set%children(child)%state == &
      transport_set%children(child)%state) .and. &
      all(transport_failed_set%children(child)%temperature == &
        transport_set%children(child)%temperature), &
      "invalid multipatch transport fine rollback")
  end do

  call composite_reactive_eb_patch_set_integral_2d( &
    coarse_state, coarse_geometry, old_set, integral_before, ok)
  call require(ok, "initial multipatch hydro integral")
  dt = 0.02_dp * min(coarse_geometry%dx, coarse_geometry%dy) / sound_speed
  call advance_reactive_eb_patch_set_hydro_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, old_set, &
    "hllc", "pcm", "mc", 2, dt, new_coarse_state, &
    new_coarse_temperature, hydro_set, ok, &
    failure_context=hydro_failure_context)
  if (.not. ok) write(*, '(a)') trim(hydro_failure_context)
  call require(ok .and. hydro_set%is_valid(coarse_geometry, nvar), &
    "subcycled multipatch EB hydro")
  call composite_reactive_eb_patch_set_integral_2d( &
    new_coarse_state, coarse_geometry, hydro_set, integral_after, ok)
  integral_scale = max(1.0_dp, maxval(abs(integral_before)))
  if (.not. ok .or. &
      abs(integral_after(irho) - integral_before(irho)) > &
        5.0e-11_dp * integral_scale .or. &
      abs(integral_after(iet) - integral_before(iet)) > &
        5.0e-11_dp * integral_scale) then
    write(*, '(a,4(es25.16e3,1x))') "rho conservation: ", &
      integral_before(irho), integral_after(irho), &
      integral_after(irho) - integral_before(irho), &
      5.0e-11_dp * integral_scale
    write(*, '(a,4(es25.16e3,1x))') "energy conservation: ", &
      integral_before(iet), integral_after(iet), &
      integral_after(iet) - integral_before(iet), &
      5.0e-11_dp * integral_scale
  end if
  call require(ok .and. abs(integral_after(irho) - integral_before(irho)) <= &
    5.0e-11_dp * integral_scale .and. &
    abs(integral_after(iet) - integral_before(iet)) <= &
    5.0e-11_dp * integral_scale, &
    "multipatch hydro mass and energy conservation")
  do k = 1, size(species)
    component = reactive_species_component(k)
    call require(abs(integral_after(component) - &
      integral_before(component)) <= 5.0e-11_dp * integral_scale, &
      "multipatch hydro species conservation")
  end do
  state_scale = max(1.0_dp, maxval(abs(state_cell)))
  call average_down_reactive_eb_patch_set_2d( &
    species, new_coarse_state, new_coarse_temperature, coarse_geometry, &
    hydro_set, averaged_state, averaged_temperature, ok)
  call require(ok .and. maxval(abs(averaged_state - &
    new_coarse_state)) <= 4.0e-13_dp * state_scale, &
    "multipatch hydro synchronized hierarchy")

  call advance_reactive_eb_patch_set_strang_2d( &
    species, reactions, coarse_state, coarse_temperature, coarse_geometry, &
    old_set, "hllc", "pcm", "mc", 2, dt, .true., 1.0e-8_dp, &
    1.0e-14_dp, chemistry_coarse_state, chemistry_coarse_temperature, &
    chemistry_set, ok, failure_context=hydro_failure_context)
  if (.not. ok) write(*, '(a)') trim(hydro_failure_context)
  call require(ok .and. chemistry_set%is_valid(coarse_geometry, nvar), &
    "multipatch chemistry-hydro transaction")
  call composite_reactive_eb_patch_set_integral_2d( &
    chemistry_coarse_state, coarse_geometry, chemistry_set, &
    integral_after, ok)
  call require(ok .and. abs(integral_after(irho) - integral_before(irho)) <= &
    5.0e-10_dp * integral_scale .and. &
    abs(integral_after(iet) - integral_before(iet)) <= &
    5.0e-10_dp * integral_scale, &
    "multipatch chemistry mass and energy conservation")
  reaction_change = 0.0_dp
  do k = 1, size(species)
    component = reactive_species_component(k)
    reaction_change = max(reaction_change, maxval(abs( &
      chemistry_coarse_state(component, :, :) - &
      new_coarse_state(component, :, :))))
    do child = 1, chemistry_set%patch_count()
      reaction_change = max(reaction_change, maxval(abs( &
        chemistry_set%children(child)%state(component, :, :) - &
        hydro_set%children(child)%state(component, :, :))))
    end do
  end do
  call require(reaction_change > 1.0e-14_dp * state_scale, &
    "multipatch chemistry changes species")
  call average_down_reactive_eb_patch_set_2d( &
    species, chemistry_coarse_state, chemistry_coarse_temperature, &
    coarse_geometry, chemistry_set, averaged_state, averaged_temperature, ok)
  call require(ok .and. maxval(abs(averaged_state - &
    chemistry_coarse_state)) <= 4.0e-13_dp * state_scale, &
    "multipatch chemistry synchronized hierarchy")

  call advance_reactive_eb_patch_set_strang_2d( &
    species, reactions, coarse_state, coarse_temperature, coarse_geometry, &
    old_set, "invalid", "pcm", "mc", 2, dt, .true., 1.0e-8_dp, &
    1.0e-14_dp, new_coarse_state, new_coarse_temperature, &
    chemistry_failed_set, ok)
  call require(.not. ok .and. &
    maxval(abs(new_coarse_state - coarse_state)) == 0.0_dp .and. &
    maxval(abs(new_coarse_temperature - coarse_temperature)) == 0.0_dp .and. &
    chemistry_failed_set%patch_count() == old_set%patch_count(), &
    "multipatch chemistry coarse rollback")
  do child = 1, old_set%patch_count()
    call require(maxval(abs( &
      chemistry_failed_set%children(child)%state - &
      old_set%children(child)%state)) == 0.0_dp .and. &
      maxval(abs(chemistry_failed_set%children(child)%temperature - &
        old_set%children(child)%temperature)) == 0.0_dp, &
      "multipatch chemistry fine rollback")
  end do

  call advance_reactive_eb_patch_set_hydro_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, old_set, &
    "invalid", "pcm", "mc", 2, dt, new_coarse_state, &
    new_coarse_temperature, failed_set, ok)
  call require(.not. ok .and. &
    maxval(abs(new_coarse_state - coarse_state)) == 0.0_dp .and. &
    maxval(abs(new_coarse_temperature - coarse_temperature)) == 0.0_dp .and. &
    failed_set%patch_count() == old_set%patch_count(), &
    "multipatch hydro coarse rollback")
  do child = 1, old_set%patch_count()
    call require(maxval(abs(failed_set%children(child)%state - &
      old_set%children(child)%state)) == 0.0_dp .and. &
      maxval(abs(failed_set%children(child)%temperature - &
        old_set%children(child)%temperature)) == 0.0_dp, &
      "multipatch hydro fine rollback")
  end do

  do child = 1, old_set%patch_count()
    do j = 1, old_set%children(child)%geometry%ny
      do i = 1, old_set%children(child)%geometry%nx
        factor = 1.0_dp + 1.0e-3_dp * &
          real(10 * child + i + 2 * j, dp)
        old_set%children(child)%state(:, i, j) = factor * &
          old_set%children(child)%state(:, i, j)
      end do
    end do
  end do
  call require(old_set%is_valid(coarse_geometry, nvar), &
    "perturbed reactive EB patch set")
  call composite_reactive_eb_patch_set_integral_2d( &
    coarse_state, coarse_geometry, old_set, integral_before, ok)
  call require(ok, "old multipatch composite integral")
  call average_down_reactive_eb_patch_set_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, old_set, &
    averaged_state, averaged_temperature, ok)
  call require(ok, "multipatch reactive average down")
  integral_after = 0.0_dp
  do i = 1, nvar
    integral_after(i) = sum(coarse_geometry%volume_fraction * &
      averaged_state(i, :, :)) * coarse_geometry%dx * coarse_geometry%dy
  end do
  integral_scale = max(1.0_dp, maxval(abs(integral_before)))
  call require(maxval(abs(integral_after - integral_before)) <= &
    8.0e-12_dp * integral_scale, "multipatch average-down conservation")

  new_tags = .false.
  new_tags(2:6, 6:10) = .true.
  new_tags(9:13, 9:13) = .true.
  call build_amr_eb_regrid_plan_collection_2d( &
    new_tags, criteria, new_collection, ok)
  call require(ok .and. new_collection%patch_count() == 2, &
    "moved two-patch plan")
  allocate(new_geometries(new_collection%patch_count()))
  do child = 1, new_collection%patch_count()
    call build_patch_geometry( &
      coarse_geometry, new_collection%plans(child)%coarse_i_lower, &
      new_collection%plans(child)%coarse_i_upper, &
      new_collection%plans(child)%coarse_j_lower, &
      new_collection%plans(child)%coarse_j_upper, ratio, &
      new_geometries(child), geometry_patch, ok)
    call require(ok, "new child EB geometry")
  end do
  call regrid_reactive_eb_patch_set_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, old_set, &
    new_geometries, new_collection, ratio, new_coarse_state, &
    new_coarse_temperature, new_set, ok)
  call require(ok .and. new_set%is_valid(coarse_geometry, nvar) .and. &
    new_set%patch_count() == 2, "transactional multipatch regrid")
  call composite_reactive_eb_patch_set_integral_2d( &
    new_coarse_state, coarse_geometry, new_set, integral_after, ok)
  call require(ok .and. maxval(abs(integral_after - integral_before)) <= &
    8.0e-12_dp * integral_scale, "multipatch regrid conservation")
  call require(maxval(abs(new_set%children(1)%state(:, 1:8, :) - &
    old_set%children(1)%state(:, 3:10, :))) == 0.0_dp, &
    "moved patch exact fine overlap")
  call require(maxval(abs(new_set%children(2)%state - &
    old_set%children(2)%state)) == 0.0_dp, &
    "unchanged patch exact state retention")
  do j = 1, new_set%children(1)%geometry%ny
    call require(maxval(abs(new_set%children(1)%state(:, 9:10, j) - &
      spread(new_coarse_state(:, 6, 6 + (j - 1) / ratio), 2, 2))) <= &
      5.0e-14_dp * state_scale, &
      "new patch cells use synchronized coarse PCM")
  end do

  old_set%children(1)%state(:, 1, 1) = &
    ieee_value(0.0_dp, ieee_quiet_nan)
  call regrid_reactive_eb_patch_set_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, old_set, &
    new_geometries, new_collection, ratio, new_coarse_state, &
    new_coarse_temperature, removed_set, ok)
  call require(.not. ok .and. .not. allocated(removed_set%children) .and. &
    maxval(abs(new_coarse_state - coarse_state)) == 0.0_dp .and. &
    maxval(abs(new_coarse_temperature - coarse_temperature)) == 0.0_dp, &
    "invalid multipatch regrid rollback")
  old_set%children(1)%state(:, 1, 1) = &
    1.01_dp * (1.0_dp + 1.0e-3_dp * 13.0_dp) * state_cell

  new_tags = .false.
  call build_amr_eb_regrid_plan_collection_2d( &
    new_tags, criteria, new_collection, ok)
  call require(ok .and. new_collection%patch_count() == 0, &
    "empty multipatch removal plan")
  deallocate(new_geometries)
  allocate(new_geometries(0))
  call regrid_reactive_eb_patch_set_2d( &
    species, coarse_state, coarse_temperature, coarse_geometry, old_set, &
    new_geometries, new_collection, ratio, new_coarse_state, &
    new_coarse_temperature, removed_set, ok)
  call require(ok .and. removed_set%is_valid(coarse_geometry, nvar) .and. &
    removed_set%patch_count() == 0, "remove complete EB patch set")
  integral_after = 0.0_dp
  do i = 1, nvar
    integral_after(i) = sum(coarse_geometry%volume_fraction * &
      new_coarse_state(i, :, :)) * coarse_geometry%dx * coarse_geometry%dy
  end do
  call require(maxval(abs(integral_after - integral_before)) <= &
    8.0e-12_dp * integral_scale, "multipatch removal conservation")

  write(*, '(a)') "test_amr_eb_multipatch_2d: PASS"

contains

  subroutine impose_double_hotspot( &
      thermo_species, composition, geometry, state, temperature, valid)
    type(nasa7_species), intent(in) :: thermo_species(:)
    real(dp), intent(in) :: composition(:)
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), intent(inout) :: state(:, :, :), temperature(:, :)
    logical, intent(out) :: valid

    real(dp) :: molecular_weight, gas_constant, cp, cv, gamma
    real(dp) :: enthalpy, internal_energy, entropy, local_temperature
    real(dp) :: local_x, local_y
    logical :: properties_ok
    integer :: local_i, local_j

    valid = .false.
    do local_j = 1, geometry%ny
      local_y = geometry%y_lower + &
        (real(local_j, dp) - 0.5_dp) * geometry%dy
      do local_i = 1, geometry%nx
        local_x = geometry%x_lower + &
          (real(local_i, dp) - 0.5_dp) * geometry%dx
        local_temperature = temperature_cell + &
          180.0_dp * exp(-((local_x - 0.26_dp)**2 + &
            (local_y - 0.58_dp)**2) / 0.008_dp) + &
          140.0_dp * exp(-((local_x - 0.78_dp)**2 + &
            (local_y - 0.77_dp)**2) / 0.006_dp)
        call mixture_mass_properties( &
          thermo_species, composition, local_temperature, molecular_weight, &
          gas_constant, cp, cv, gamma, enthalpy, internal_energy, entropy, &
          properties_ok)
        if (.not. properties_ok) return
        state(iet, local_i, local_j) = &
          state(irho, local_i, local_j) * internal_energy
        temperature(local_i, local_j) = local_temperature
      end do
    end do
    valid = .true.
  end subroutine impose_double_hotspot

  real(dp) function patch_set_temperature_span( &
      root_temperature, root_geometry, patch_set) result(span)
    real(dp), intent(in) :: root_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: root_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set

    real(dp) :: minimum_temperature, maximum_temperature
    integer :: local_child, local_i, local_j

    minimum_temperature = huge(1.0_dp)
    maximum_temperature = -huge(1.0_dp)
    do local_j = 1, root_geometry%ny
      do local_i = 1, root_geometry%nx
        if (root_geometry%cell_type(local_i, local_j) == eb_covered_cell .or. &
            coarse_cell_is_refined(patch_set, local_i, local_j)) cycle
        minimum_temperature = min( &
          minimum_temperature, root_temperature(local_i, local_j))
        maximum_temperature = max( &
          maximum_temperature, root_temperature(local_i, local_j))
      end do
    end do
    do local_child = 1, patch_set%patch_count()
      do local_j = 1, patch_set%children(local_child)%geometry%ny
        do local_i = 1, patch_set%children(local_child)%geometry%nx
          if (patch_set%children(local_child)%geometry%cell_type( &
              local_i, local_j) == eb_covered_cell) cycle
          minimum_temperature = min(minimum_temperature, &
            patch_set%children(local_child)%temperature(local_i, local_j))
          maximum_temperature = max(maximum_temperature, &
            patch_set%children(local_child)%temperature(local_i, local_j))
        end do
      end do
    end do
    span = maximum_temperature - minimum_temperature
  end function patch_set_temperature_span

  pure logical function coarse_cell_is_refined(patch_set, i, j) &
      result(refined)
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set
    integer, intent(in) :: i, j
    integer :: local_child

    refined = .false.
    do local_child = 1, patch_set%patch_count()
      refined = &
        i >= patch_set%children(local_child)%patch%coarse_i_lower .and. &
        i <= patch_set%children(local_child)%patch%coarse_i_upper .and. &
        j >= patch_set%children(local_child)%patch%coarse_j_lower .and. &
        j <= patch_set%children(local_child)%patch%coarse_j_upper
      if (refined) return
    end do
  end function coarse_cell_is_refined

  subroutine assert_patch_set_covered_unchanged( &
      original_root, original_set, candidate_root, candidate_set, &
      root_geometry, message)
    real(dp), intent(in) :: original_root(:, :, :), candidate_root(:, :, :)
    type(reactive_eb_patch_set_2d), intent(in) :: original_set, candidate_set
    type(eb_geometry_2d), intent(in) :: root_geometry
    character(len=*), intent(in) :: message
    integer :: local_child, local_i, local_j

    do local_j = 1, root_geometry%ny
      do local_i = 1, root_geometry%nx
        if (root_geometry%cell_type(local_i, local_j) /= &
            eb_covered_cell) cycle
        call require(all(candidate_root(:, local_i, local_j) == &
          original_root(:, local_i, local_j)), message)
      end do
    end do
    do local_child = 1, original_set%patch_count()
      do local_j = 1, original_set%children(local_child)%geometry%ny
        do local_i = 1, original_set%children(local_child)%geometry%nx
          if (original_set%children(local_child)%geometry%cell_type( &
              local_i, local_j) /= eb_covered_cell) cycle
          call require(all(candidate_set%children(local_child)%state( &
            :, local_i, local_j) == &
            original_set%children(local_child)%state(:, local_i, local_j)), &
            message)
        end do
      end do
    end do
  end subroutine assert_patch_set_covered_unchanged

  subroutine build_patch_geometry( &
      root_geometry, i_lower, i_upper, j_lower, j_upper, refinement_ratio, &
      fine_geometry, patch, valid)
    type(eb_geometry_2d), intent(in) :: root_geometry
    integer, intent(in) :: i_lower, i_upper, j_lower, j_upper
    integer, intent(in) :: refinement_ratio
    type(eb_geometry_2d), intent(out) :: fine_geometry
    type(amr_eb_patch_2d), intent(out) :: patch
    logical, intent(out) :: valid

    real(dp), allocatable :: level_set(:, :)
    real(dp) :: x_lower, x_upper, y_lower, y_upper, local_x, local_y
    integer :: fine_nx, fine_ny, local_i, local_j

    fine_nx = (i_upper - i_lower + 1) * refinement_ratio
    fine_ny = (j_upper - j_lower + 1) * refinement_ratio
    x_lower = root_geometry%x_lower + real(i_lower - 1, dp) * &
      root_geometry%dx
    x_upper = root_geometry%x_lower + real(i_upper, dp) * root_geometry%dx
    y_lower = root_geometry%y_lower + real(j_lower - 1, dp) * &
      root_geometry%dy
    y_upper = root_geometry%y_lower + real(j_upper, dp) * root_geometry%dy
    allocate(level_set(0:fine_nx, 0:fine_ny))
    do local_j = 0, fine_ny
      local_y = y_lower + real(local_j, dp) * &
        (y_upper - y_lower) / real(fine_ny, dp)
      do local_i = 0, fine_nx
        local_x = x_lower + real(local_i, dp) * &
          (x_upper - x_lower) / real(fine_nx, dp)
        level_set(local_i, local_j) = local_x + local_y - 0.78_dp
      end do
    end do
    call build_eb_geometry_2d( &
      level_set, x_lower, x_upper, y_lower, y_upper, fine_geometry, valid)
    if (.not. valid) return
    call build_amr_eb_patch_2d( &
      root_geometry, fine_geometry, i_lower, i_upper, j_lower, j_upper, &
      refinement_ratio, patch, valid)
  end subroutine build_patch_geometry

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_amr_eb_multipatch_2d
