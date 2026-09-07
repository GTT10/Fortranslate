program test_reactive_eb_chemistry_3d
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use reactive_1d_mod, only: reactive_nvar, reactive_species_component
  use reactive_3d_mod, only: advance_reactive_chemistry_3d
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d, eb_regular_cell_3d, &
    build_axis_plane_eb_geometry_3d
  use eb_reactive_hydro_3d_mod, only: &
    advance_reactive_eb_state_redistributed_euler_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: &
    initialize_reactive_eb_density_sheet_3d, &
    advance_reactive_eb_strang_3d, reactive_eb_integrals_3d, &
    reactive_eb_element_integrals_3d
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
  write(*, '(a)') "test_reactive_eb_chemistry_3d: PASS"

contains

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
