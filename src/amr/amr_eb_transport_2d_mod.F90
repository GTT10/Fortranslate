module amr_eb_transport_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use transport_database_mod, only: gas_transport_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use reactive_boundary_2d_mod, only: reactive_boundary_set_2d
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_covered_cell
  use eb_reactive_reconstruction_2d_mod, only: &
    reactive_eb_exterior_state_2d
  use eb_reactive_redistribution_2d_mod, only: &
    advance_reactive_eb_state_redistributed_2d
  use eb_reactive_transport_2d_mod, only: &
    reactive_eb_transport_fluxes_rhs_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, average_down_reactive_eb_state_patch_2d
  use amr_eb_flux_register_2d_mod, only: &
    amr_eb_flux_register_2d, initialize_amr_eb_flux_register_2d, &
    accumulate_coarse_eb_fluxes_2d, accumulate_fine_eb_fluxes_2d, &
    reflux_reactive_eb_state_patch_2d
  use amr_eb_reactive_2d_mod, only: build_reactive_eb_patch_exterior_2d
  implicit none
  private

  public :: advance_two_level_reactive_eb_transport_euler_2d
  public :: advance_two_level_reactive_eb_transport_2d

contains

  subroutine advance_two_level_reactive_eb_transport_euler_2d( &
      species, transport, coarse_state, coarse_temperature, &
      coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
      dt, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, max_order, new_coarse_state, &
      new_coarse_temperature, new_fine_state, new_fine_temperature, &
      minimum_theta, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: coarse_state(:, :, :), coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    real(dp), intent(in) :: fine_state(:, :, :), fine_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    real(dp), intent(in) :: dt
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    real(dp), intent(in) :: target_volume_fraction
    integer, intent(in) :: max_order
    real(dp), intent(out) :: new_coarse_state(:, :, :)
    real(dp), intent(out) :: new_coarse_temperature(:, :)
    real(dp), intent(out) :: new_fine_state(:, :, :)
    real(dp), intent(out) :: new_fine_temperature(:, :)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok

    type(amr_eb_flux_register_2d) :: flux_register
    type(reactive_eb_exterior_state_2d) :: exterior
    real(dp), allocatable :: coarse_candidate(:, :, :)
    real(dp), allocatable :: coarse_candidate_temperature(:, :)
    real(dp), allocatable :: coarse_rhs(:, :, :)
    real(dp), allocatable :: coarse_x_flux(:, :, :)
    real(dp), allocatable :: coarse_y_flux(:, :, :)
    real(dp), allocatable :: coarse_work(:, :, :)
    real(dp), allocatable :: coarse_work_temperature(:, :)
    real(dp), allocatable :: fine_candidate(:, :, :)
    real(dp), allocatable :: fine_candidate_temperature(:, :)
    real(dp), allocatable :: fine_rhs(:, :, :)
    real(dp), allocatable :: fine_x_flux(:, :, :)
    real(dp), allocatable :: fine_y_flux(:, :, :)
    real(dp), allocatable :: fine_work(:, :, :)
    real(dp), allocatable :: fine_work_temperature(:, :)
    real(dp) :: alpha, coarse_theta, fine_theta, fine_dt
    logical :: local_ok
    integer :: nvar, ratio, substep

    new_coarse_state = coarse_state
    new_coarse_temperature = coarse_temperature
    new_fine_state = fine_state
    new_fine_temperature = fine_temperature
    minimum_theta = 1.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar < 1 .or. size(transport) /= size(species) .or. &
        .not. ieee_is_finite(dt) .or. dt <= 0.0_dp .or. &
        .not. patch%is_valid(coarse_geometry, fine_geometry) .or. &
        any(shape(new_coarse_state) /= shape(coarse_state)) .or. &
        any(shape(new_coarse_temperature) /= shape(coarse_temperature)) .or. &
        any(shape(new_fine_state) /= shape(fine_state)) .or. &
        any(shape(new_fine_temperature) /= shape(fine_temperature))) return

    allocate(coarse_candidate, mold=coarse_state)
    allocate(coarse_candidate_temperature, mold=coarse_temperature)
    allocate(coarse_rhs, mold=coarse_state)
    allocate(coarse_x_flux(nvar, 0:coarse_geometry%nx, coarse_geometry%ny))
    allocate(coarse_y_flux(nvar, coarse_geometry%nx, 0:coarse_geometry%ny))
    call reactive_eb_transport_fluxes_rhs_2d( &
      species, transport, coarse_state, coarse_temperature, coarse_geometry, &
      dt, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      coarse_rhs, coarse_x_flux, coarse_y_flux, coarse_theta, local_ok)
    if (.not. local_ok) return
    call advance_reactive_eb_state_redistributed_2d( &
      species, coarse_state, coarse_temperature, coarse_geometry, &
      coarse_rhs, dt, coarse_candidate, coarse_candidate_temperature, &
      local_ok, target_volume_fraction, max_order)
    if (.not. local_ok) return

    call initialize_amr_eb_flux_register_2d( &
      coarse_geometry, fine_geometry, patch, nvar, flux_register, local_ok)
    if (.not. local_ok) return
    call accumulate_coarse_eb_fluxes_2d( &
      flux_register, coarse_geometry, fine_geometry, patch, &
      coarse_x_flux, coarse_y_flux, dt, local_ok)
    if (.not. local_ok) return

    allocate(fine_candidate, source=fine_state)
    allocate(fine_candidate_temperature, source=fine_temperature)
    allocate(fine_rhs, mold=fine_state)
    allocate(fine_x_flux(nvar, 0:fine_geometry%nx, fine_geometry%ny))
    allocate(fine_y_flux(nvar, fine_geometry%nx, 0:fine_geometry%ny))
    allocate(fine_work, mold=fine_state)
    allocate(fine_work_temperature, mold=fine_temperature)
    ratio = patch%refinement_ratio
    fine_dt = dt / real(ratio, dp)
    do substep = 1, ratio
      alpha = real(substep - 1, dp) / real(ratio, dp)
      call build_reactive_eb_patch_exterior_2d( &
        species, coarse_state, coarse_temperature, coarse_candidate, &
        coarse_candidate_temperature, coarse_geometry, fine_geometry, &
        patch, alpha, exterior, local_ok, fine_candidate, &
        fine_candidate_temperature)
      if (.not. local_ok) return
      call reactive_eb_transport_fluxes_rhs_2d( &
        species, transport, fine_candidate, fine_candidate_temperature, &
        fine_geometry, fine_dt, viscosity_enabled, &
        thermal_conduction_enabled, species_diffusion_enabled, &
        barodiffusion_enabled, boundaries, fine_rhs, fine_x_flux, &
        fine_y_flux, fine_theta, local_ok, exterior)
      if (.not. local_ok) return
      minimum_theta = min(minimum_theta, fine_theta)
      call advance_reactive_eb_state_redistributed_2d( &
        species, fine_candidate, fine_candidate_temperature, fine_geometry, &
        fine_rhs, fine_dt, fine_work, fine_work_temperature, local_ok, &
        target_volume_fraction, max_order)
      if (.not. local_ok) return
      fine_candidate = fine_work
      fine_candidate_temperature = fine_work_temperature
      call accumulate_fine_eb_fluxes_2d( &
        flux_register, coarse_geometry, fine_geometry, patch, &
        fine_x_flux, fine_y_flux, fine_dt, local_ok)
      if (.not. local_ok) return
    end do
    minimum_theta = min(minimum_theta, coarse_theta)

    allocate(coarse_work, mold=coarse_state)
    allocate(coarse_work_temperature, mold=coarse_temperature)
    call reflux_reactive_eb_state_patch_2d( &
      species, coarse_candidate, coarse_candidate_temperature, &
      coarse_geometry, fine_candidate, fine_candidate_temperature, &
      fine_geometry, patch, flux_register, coarse_work, &
      coarse_work_temperature, fine_work, fine_work_temperature, local_ok)
    if (.not. local_ok) return
    call average_down_reactive_eb_state_patch_2d( &
      species, coarse_work, coarse_work_temperature, coarse_geometry, &
      fine_work, fine_geometry, patch, coarse_candidate, &
      coarse_candidate_temperature, local_ok)
    if (.not. local_ok) return

    new_coarse_state = coarse_candidate
    new_coarse_temperature = coarse_candidate_temperature
    new_fine_state = fine_work
    new_fine_temperature = fine_work_temperature
    ok = .true.
  end subroutine advance_two_level_reactive_eb_transport_euler_2d

  subroutine recover_transport_temperature_2d( &
      species, state, temperature_guess, geometry, temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :), temperature_guess(:, :)
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), intent(out) :: temperature(:, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: sound_speed
    logical :: local_ok
    integer :: i, j

    temperature = temperature_guess
    ok = .false.
    if (any(shape(state) /= &
        [reactive_nvar(size(species)), geometry%nx, geometry%ny]) .or. &
        any(shape(temperature_guess) /= [geometry%nx, geometry%ny]) .or. &
        any(shape(temperature) /= shape(temperature_guess))) return
    allocate(primitive(reactive_nprim(size(species))))
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) == eb_covered_cell) cycle
        call reactive_conserved_to_primitive( &
          species, state(:, i, j), temperature_guess(i, j), primitive, &
          temperature(i, j), sound_speed, local_ok)
        if (.not. local_ok) return
      end do
    end do
    ok = .true.
  end subroutine recover_transport_temperature_2d

  subroutine advance_two_level_reactive_eb_transport_2d( &
      species, transport, coarse_state, coarse_temperature, &
      coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
      interval, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, max_order, new_coarse_state, &
      new_coarse_temperature, new_fine_state, new_fine_temperature, &
      minimum_theta, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: coarse_state(:, :, :), coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    real(dp), intent(in) :: fine_state(:, :, :), fine_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    real(dp), intent(in) :: interval
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    real(dp), intent(in) :: target_volume_fraction
    integer, intent(in) :: max_order
    real(dp), intent(out) :: new_coarse_state(:, :, :)
    real(dp), intent(out) :: new_coarse_temperature(:, :)
    real(dp), intent(out) :: new_fine_state(:, :, :)
    real(dp), intent(out) :: new_fine_temperature(:, :)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok

    real(dp), allocatable :: stage_coarse_state(:, :, :)
    real(dp), allocatable :: stage_coarse_temperature(:, :)
    real(dp), allocatable :: stage_fine_state(:, :, :)
    real(dp), allocatable :: stage_fine_temperature(:, :)
    real(dp), allocatable :: euler_coarse_state(:, :, :)
    real(dp), allocatable :: euler_coarse_temperature(:, :)
    real(dp), allocatable :: euler_fine_state(:, :, :)
    real(dp), allocatable :: euler_fine_temperature(:, :)
    real(dp), allocatable :: candidate_coarse_state(:, :, :)
    real(dp), allocatable :: candidate_coarse_temperature(:, :)
    real(dp), allocatable :: candidate_fine_state(:, :, :)
    real(dp), allocatable :: candidate_fine_temperature(:, :)
    real(dp) :: theta_one, theta_two
    logical :: local_ok

    new_coarse_state = coarse_state
    new_coarse_temperature = coarse_temperature
    new_fine_state = fine_state
    new_fine_temperature = fine_temperature
    minimum_theta = 1.0_dp
    ok = .false.
    if (.not. ieee_is_finite(interval) .or. interval < 0.0_dp) return
    if (interval <= tiny(1.0_dp) .or. .not. (viscosity_enabled .or. &
        thermal_conduction_enabled .or. species_diffusion_enabled)) then
      ok = .true.
      return
    end if

    allocate(stage_coarse_state, mold=coarse_state)
    allocate(stage_coarse_temperature, mold=coarse_temperature)
    allocate(stage_fine_state, mold=fine_state)
    allocate(stage_fine_temperature, mold=fine_temperature)
    call advance_two_level_reactive_eb_transport_euler_2d( &
      species, transport, coarse_state, coarse_temperature, &
      coarse_geometry, fine_state, fine_temperature, fine_geometry, patch, &
      interval, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, max_order, stage_coarse_state, &
      stage_coarse_temperature, stage_fine_state, stage_fine_temperature, &
      theta_one, local_ok)
    if (.not. local_ok) return

    allocate(euler_coarse_state, mold=coarse_state)
    allocate(euler_coarse_temperature, mold=coarse_temperature)
    allocate(euler_fine_state, mold=fine_state)
    allocate(euler_fine_temperature, mold=fine_temperature)
    call advance_two_level_reactive_eb_transport_euler_2d( &
      species, transport, stage_coarse_state, stage_coarse_temperature, &
      coarse_geometry, stage_fine_state, stage_fine_temperature, &
      fine_geometry, patch, interval, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, target_volume_fraction, max_order, &
      euler_coarse_state, euler_coarse_temperature, euler_fine_state, &
      euler_fine_temperature, theta_two, local_ok)
    if (.not. local_ok) return

    allocate(candidate_coarse_state, source= &
      0.5_dp * (coarse_state + euler_coarse_state))
    allocate(candidate_fine_state, source= &
      0.5_dp * (fine_state + euler_fine_state))
    allocate(candidate_coarse_temperature, mold=coarse_temperature)
    allocate(candidate_fine_temperature, mold=fine_temperature)
    call recover_transport_temperature_2d( &
      species, candidate_coarse_state, &
      0.5_dp * (coarse_temperature + euler_coarse_temperature), &
      coarse_geometry, candidate_coarse_temperature, local_ok)
    if (.not. local_ok) return
    call recover_transport_temperature_2d( &
      species, candidate_fine_state, &
      0.5_dp * (fine_temperature + euler_fine_temperature), fine_geometry, &
      candidate_fine_temperature, local_ok)
    if (.not. local_ok) return
    call average_down_reactive_eb_state_patch_2d( &
      species, candidate_coarse_state, candidate_coarse_temperature, &
      coarse_geometry, candidate_fine_state, fine_geometry, patch, &
      new_coarse_state, new_coarse_temperature, local_ok)
    if (.not. local_ok) return
    new_fine_state = candidate_fine_state
    new_fine_temperature = candidate_fine_temperature
    minimum_theta = min(theta_one, theta_two)
    ok = .true.
  end subroutine advance_two_level_reactive_eb_transport_2d

end module amr_eb_transport_2d_mod
