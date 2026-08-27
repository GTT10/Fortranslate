module amr_eb_multilevel_transport_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use transport_database_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar
  use reactive_boundary_2d_mod, only: reactive_boundary_set_2d
  use eb_geometry_2d_mod, only: eb_geometry_2d
  use eb_reactive_reconstruction_2d_mod, only: &
    reactive_eb_exterior_state_2d
  use eb_reactive_redistribution_2d_mod, only: &
    advance_reactive_eb_state_redistributed_2d
  use eb_reactive_transport_2d_mod, only: &
    reactive_eb_transport_fluxes_rhs_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, average_down_reactive_eb_state_patch_2d, &
    composite_eb_integral_2d
  use amr_eb_multilevel_2d_mod, only: &
    average_down_three_level_reactive_eb_state_2d, &
    composite_three_level_eb_integral_2d
  use amr_eb_flux_register_2d_mod, only: &
    amr_eb_flux_register_2d, initialize_amr_eb_flux_register_2d, &
    accumulate_coarse_eb_fluxes_2d, accumulate_fine_eb_fluxes_2d, &
    reflux_reactive_eb_state_patch_2d
  use amr_eb_reactive_2d_mod, only: build_reactive_eb_patch_exterior_2d
  use amr_eb_transport_2d_mod, only: recover_transport_temperature_2d
  use amr_eb_multilevel_reactive_2d_mod, only: &
    level_two_interface_is_regular, close_cut_interface_conservation_2d
  implicit none
  private

  public :: advance_three_level_reactive_eb_transport_euler_2d
  public :: advance_three_level_reactive_eb_transport_2d

contains

  subroutine advance_three_level_reactive_eb_transport_euler_2d( &
      species, transport, root_state, root_temperature, root_geometry, &
      level_one_state, level_one_temperature, level_one_geometry, root_patch, &
      level_two_state, level_two_temperature, level_two_geometry, &
      level_one_patch, dt, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, max_order, new_root_state, new_root_temperature, &
      new_level_one_state, new_level_one_temperature, new_level_two_state, &
      new_level_two_temperature, minimum_theta, ok, failure_context)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: root_state(:, :, :), root_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: root_geometry
    real(dp), intent(in) :: level_one_state(:, :, :)
    real(dp), intent(in) :: level_one_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: level_one_geometry
    type(amr_eb_patch_2d), intent(in) :: root_patch
    real(dp), intent(in) :: level_two_state(:, :, :)
    real(dp), intent(in) :: level_two_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: level_two_geometry
    type(amr_eb_patch_2d), intent(in) :: level_one_patch
    real(dp), intent(in) :: dt, target_volume_fraction
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    integer, intent(in) :: max_order
    real(dp), intent(out) :: new_root_state(:, :, :)
    real(dp), intent(out) :: new_root_temperature(:, :)
    real(dp), intent(out) :: new_level_one_state(:, :, :)
    real(dp), intent(out) :: new_level_one_temperature(:, :)
    real(dp), intent(out) :: new_level_two_state(:, :, :)
    real(dp), intent(out) :: new_level_two_temperature(:, :)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    character(len=*), intent(out), optional :: failure_context

    type(amr_eb_flux_register_2d) :: root_register, level_one_register
    type(reactive_eb_exterior_state_2d) :: level_one_exterior
    type(reactive_eb_exterior_state_2d) :: level_two_exterior
    real(dp), allocatable :: root_candidate(:, :, :), root_temperature_work(:, :)
    real(dp), allocatable :: root_rhs(:, :, :), root_x_flux(:, :, :)
    real(dp), allocatable :: root_y_flux(:, :, :), root_refluxed(:, :, :)
    real(dp), allocatable :: root_refluxed_temperature(:, :)
    real(dp), allocatable :: root_closed(:, :, :)
    real(dp), allocatable :: root_closed_temperature(:, :)
    real(dp), allocatable :: level_one_candidate(:, :, :)
    real(dp), allocatable :: level_one_candidate_temperature(:, :)
    real(dp), allocatable :: level_one_start(:, :, :)
    real(dp), allocatable :: level_one_start_temperature(:, :)
    real(dp), allocatable :: level_one_rhs(:, :, :)
    real(dp), allocatable :: level_one_uncorrected(:, :, :)
    real(dp), allocatable :: level_one_uncorrected_temperature(:, :)
    real(dp), allocatable :: level_one_refluxed(:, :, :)
    real(dp), allocatable :: level_one_refluxed_temperature(:, :)
    real(dp), allocatable :: level_one_closed(:, :, :)
    real(dp), allocatable :: level_one_closed_temperature(:, :)
    real(dp), allocatable :: level_one_x_flux(:, :, :)
    real(dp), allocatable :: level_one_y_flux(:, :, :)
    real(dp), allocatable :: level_two_candidate(:, :, :)
    real(dp), allocatable :: level_two_candidate_temperature(:, :)
    real(dp), allocatable :: level_two_rhs(:, :, :)
    real(dp), allocatable :: level_two_work(:, :, :)
    real(dp), allocatable :: level_two_work_temperature(:, :)
    real(dp), allocatable :: level_two_refluxed(:, :, :)
    real(dp), allocatable :: level_two_refluxed_temperature(:, :)
    real(dp), allocatable :: level_two_x_flux(:, :, :)
    real(dp), allocatable :: level_two_y_flux(:, :, :)
    real(dp), allocatable :: level_one_integral_before(:)
    real(dp), allocatable :: root_integral_before(:)
    real(dp) :: level_one_dt, level_two_dt, alpha
    real(dp) :: root_theta, level_one_theta, level_two_theta
    logical :: local_ok
    integer :: nvar, level_one_ratio, level_two_ratio
    integer :: level_one_substep, level_two_substep

    new_root_state = root_state
    new_root_temperature = root_temperature
    new_level_one_state = level_one_state
    new_level_one_temperature = level_one_temperature
    new_level_two_state = level_two_state
    new_level_two_temperature = level_two_temperature
    minimum_theta = 1.0_dp
    ok = .false.
    if (present(failure_context)) failure_context = "input validation"
    nvar = reactive_nvar(size(species))
    if (nvar < 1 .or. size(transport) /= size(species) .or. &
        .not. ieee_is_finite(dt) .or. dt <= 0.0_dp .or. &
        .not. root_patch%is_valid(root_geometry, level_one_geometry) .or. &
        .not. level_one_patch%is_valid( &
          level_one_geometry, level_two_geometry) .or. &
        level_one_patch%coarse_i_lower < 3 .or. &
        level_one_patch%coarse_i_upper > level_one_geometry%nx - 2 .or. &
        level_one_patch%coarse_j_lower < 3 .or. &
        level_one_patch%coarse_j_upper > level_one_geometry%ny - 2 .or. &
        any(shape(new_root_state) /= shape(root_state)) .or. &
        any(shape(new_root_temperature) /= shape(root_temperature)) .or. &
        any(shape(new_level_one_state) /= shape(level_one_state)) .or. &
        any(shape(new_level_one_temperature) /= &
          shape(level_one_temperature)) .or. &
        any(shape(new_level_two_state) /= shape(level_two_state)) .or. &
        any(shape(new_level_two_temperature) /= &
          shape(level_two_temperature))) return

    allocate(root_integral_before(nvar))
    if (present(failure_context)) failure_context = "initial composite integral"
    call composite_three_level_eb_integral_2d( &
      root_state, root_geometry, level_one_state, level_one_geometry, &
      root_patch, level_two_state, level_two_geometry, level_one_patch, &
      root_integral_before, local_ok)
    if (.not. local_ok) return

    allocate(root_candidate, mold=root_state)
    allocate(root_temperature_work, mold=root_temperature)
    allocate(root_rhs, mold=root_state)
    allocate(root_x_flux(nvar, 0:root_geometry%nx, root_geometry%ny))
    allocate(root_y_flux(nvar, root_geometry%nx, 0:root_geometry%ny))
    if (present(failure_context)) failure_context = "root transport flux"
    call reactive_eb_transport_fluxes_rhs_2d( &
      species, transport, root_state, root_temperature, root_geometry, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, root_rhs, &
      root_x_flux, root_y_flux, root_theta, local_ok)
    if (.not. local_ok) return
    if (present(failure_context)) failure_context = "root redistribution"
    call advance_reactive_eb_state_redistributed_2d( &
      species, root_state, root_temperature, root_geometry, root_rhs, dt, &
      root_candidate, root_temperature_work, local_ok, &
      target_volume_fraction, max_order)
    if (.not. local_ok) return
    if (present(failure_context)) failure_context = "root register initialization"
    call initialize_amr_eb_flux_register_2d( &
      root_geometry, level_one_geometry, root_patch, nvar, root_register, &
      local_ok)
    if (.not. local_ok) return
    if (present(failure_context)) &
      failure_context = "root coarse flux accumulation"
    call accumulate_coarse_eb_fluxes_2d( &
      root_register, root_geometry, level_one_geometry, root_patch, &
      root_x_flux, root_y_flux, dt, local_ok)
    if (.not. local_ok) return

    allocate(level_one_candidate, source=level_one_state)
    allocate(level_one_candidate_temperature, source=level_one_temperature)
    allocate(level_one_start, mold=level_one_state)
    allocate(level_one_start_temperature, mold=level_one_temperature)
    allocate(level_one_rhs, mold=level_one_state)
    allocate(level_one_uncorrected, mold=level_one_state)
    allocate(level_one_uncorrected_temperature, mold=level_one_temperature)
    allocate(level_one_refluxed, mold=level_one_state)
    allocate(level_one_refluxed_temperature, mold=level_one_temperature)
    allocate(level_one_closed, mold=level_one_state)
    allocate(level_one_closed_temperature, mold=level_one_temperature)
    allocate(level_one_x_flux( &
      nvar, 0:level_one_geometry%nx, level_one_geometry%ny))
    allocate(level_one_y_flux( &
      nvar, level_one_geometry%nx, 0:level_one_geometry%ny))
    allocate(level_two_candidate, source=level_two_state)
    allocate(level_two_candidate_temperature, source=level_two_temperature)
    allocate(level_two_rhs, mold=level_two_state)
    allocate(level_two_work, mold=level_two_state)
    allocate(level_two_work_temperature, mold=level_two_temperature)
    allocate(level_two_refluxed, mold=level_two_state)
    allocate(level_two_refluxed_temperature, mold=level_two_temperature)
    allocate(level_two_x_flux( &
      nvar, 0:level_two_geometry%nx, level_two_geometry%ny))
    allocate(level_two_y_flux( &
      nvar, level_two_geometry%nx, 0:level_two_geometry%ny))
    allocate(level_one_integral_before(nvar))

    level_one_ratio = root_patch%refinement_ratio
    level_two_ratio = level_one_patch%refinement_ratio
    level_one_dt = dt / real(level_one_ratio, dp)
    level_two_dt = level_one_dt / real(level_two_ratio, dp)
    do level_one_substep = 1, level_one_ratio
      if (present(failure_context)) &
        failure_context = "level-one composite integral"
      call composite_eb_integral_2d( &
        level_one_candidate, level_one_geometry, level_two_candidate, &
        level_two_geometry, level_one_patch, level_one_integral_before, &
        local_ok)
      if (.not. local_ok) return
      level_one_start = level_one_candidate
      level_one_start_temperature = level_one_candidate_temperature
      alpha = real(level_one_substep - 1, dp) / real(level_one_ratio, dp)
      if (present(failure_context)) failure_context = "level-one exterior"
      call build_reactive_eb_patch_exterior_2d( &
        species, root_state, root_temperature, root_candidate, &
        root_temperature_work, root_geometry, level_one_geometry, root_patch, &
        alpha, level_one_exterior, local_ok, level_one_candidate, &
        level_one_candidate_temperature)
      if (.not. local_ok) return
      if (present(failure_context)) failure_context = "level-one transport flux"
      call reactive_eb_transport_fluxes_rhs_2d( &
        species, transport, level_one_candidate, &
        level_one_candidate_temperature, level_one_geometry, level_one_dt, &
        viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, boundaries, &
        level_one_rhs, level_one_x_flux, level_one_y_flux, level_one_theta, &
        local_ok, level_one_exterior)
      if (.not. local_ok) return
      minimum_theta = min(minimum_theta, level_one_theta)
      if (present(failure_context)) failure_context = "level-one redistribution"
      call advance_reactive_eb_state_redistributed_2d( &
        species, level_one_candidate, level_one_candidate_temperature, &
        level_one_geometry, level_one_rhs, level_one_dt, &
        level_one_uncorrected, level_one_uncorrected_temperature, local_ok, &
        target_volume_fraction, max_order)
      if (.not. local_ok) return
      if (present(failure_context)) &
        failure_context = "root fine flux accumulation"
      call accumulate_fine_eb_fluxes_2d( &
        root_register, root_geometry, level_one_geometry, root_patch, &
        level_one_x_flux, level_one_y_flux, level_one_dt, local_ok)
      if (.not. local_ok) return

      if (present(failure_context)) &
        failure_context = "level-two register initialization"
      call initialize_amr_eb_flux_register_2d( &
        level_one_geometry, level_two_geometry, level_one_patch, nvar, &
        level_one_register, local_ok)
      if (.not. local_ok) return
      if (present(failure_context)) &
        failure_context = "level-two coarse flux accumulation"
      call accumulate_coarse_eb_fluxes_2d( &
        level_one_register, level_one_geometry, level_two_geometry, &
        level_one_patch, level_one_x_flux, level_one_y_flux, level_one_dt, &
        local_ok)
      if (.not. local_ok) return

      do level_two_substep = 1, level_two_ratio
        alpha = real(level_two_substep - 1, dp) / real(level_two_ratio, dp)
        if (present(failure_context)) failure_context = "level-two exterior"
        call build_reactive_eb_patch_exterior_2d( &
          species, level_one_start, level_one_start_temperature, &
          level_one_uncorrected, level_one_uncorrected_temperature, &
          level_one_geometry, level_two_geometry, level_one_patch, alpha, &
          level_two_exterior, local_ok, level_two_candidate, &
          level_two_candidate_temperature)
        if (.not. local_ok) return
        if (present(failure_context)) &
          failure_context = "level-two transport flux"
        call reactive_eb_transport_fluxes_rhs_2d( &
          species, transport, level_two_candidate, &
          level_two_candidate_temperature, level_two_geometry, level_two_dt, &
          viscosity_enabled, thermal_conduction_enabled, &
          species_diffusion_enabled, barodiffusion_enabled, boundaries, &
          level_two_rhs, level_two_x_flux, level_two_y_flux, level_two_theta, &
          local_ok, level_two_exterior)
        if (.not. local_ok) return
        minimum_theta = min(minimum_theta, level_two_theta)
        if (present(failure_context)) &
          failure_context = "level-two redistribution"
        call advance_reactive_eb_state_redistributed_2d( &
          species, level_two_candidate, level_two_candidate_temperature, &
          level_two_geometry, level_two_rhs, level_two_dt, level_two_work, &
          level_two_work_temperature, local_ok, target_volume_fraction, &
          max_order)
        if (.not. local_ok) return
        level_two_candidate = level_two_work
        level_two_candidate_temperature = level_two_work_temperature
        if (present(failure_context)) &
          failure_context = "level-two fine flux accumulation"
        call accumulate_fine_eb_fluxes_2d( &
          level_one_register, level_one_geometry, level_two_geometry, &
          level_one_patch, level_two_x_flux, level_two_y_flux, level_two_dt, &
          local_ok)
        if (.not. local_ok) return
      end do

      if (present(failure_context)) failure_context = "level-two reflux"
      call reflux_reactive_eb_state_patch_2d( &
        species, level_one_uncorrected, &
        level_one_uncorrected_temperature, level_one_geometry, &
        level_two_candidate, level_two_candidate_temperature, &
        level_two_geometry, level_one_patch, level_one_register, &
        level_one_refluxed, level_one_refluxed_temperature, &
        level_two_refluxed, level_two_refluxed_temperature, local_ok)
      if (.not. local_ok) return
      if (present(failure_context)) failure_context = "level-two average-down"
      call average_down_reactive_eb_state_patch_2d( &
        species, level_one_refluxed, level_one_refluxed_temperature, &
        level_one_geometry, level_two_refluxed, level_two_geometry, &
        level_one_patch, level_one_candidate, &
        level_one_candidate_temperature, local_ok)
      if (.not. local_ok) return
      level_two_candidate = level_two_refluxed
      level_two_candidate_temperature = level_two_refluxed_temperature
      if (.not. level_two_interface_is_regular(level_two_geometry)) then
        if (present(failure_context)) &
          failure_context = "level-two conservation closure"
        call close_cut_interface_conservation_2d( &
          species, level_one_integral_before, level_one_candidate, &
          level_one_candidate_temperature, level_one_geometry, &
          level_two_candidate, level_two_geometry, level_one_patch, &
          level_one_x_flux, level_one_y_flux, level_one_dt, &
          level_one_closed, level_one_closed_temperature, local_ok)
        if (.not. local_ok) return
        level_one_candidate = level_one_closed
        level_one_candidate_temperature = level_one_closed_temperature
      end if
    end do

    minimum_theta = min(minimum_theta, root_theta)
    allocate(root_refluxed, mold=root_state)
    allocate(root_refluxed_temperature, mold=root_temperature)
    if (present(failure_context)) failure_context = "root reflux"
    call reflux_reactive_eb_state_patch_2d( &
      species, root_candidate, root_temperature_work, root_geometry, &
      level_one_candidate, level_one_candidate_temperature, &
      level_one_geometry, root_patch, root_register, root_refluxed, &
      root_refluxed_temperature, level_one_refluxed, &
      level_one_refluxed_temperature, local_ok)
    if (.not. local_ok) return
    if (present(failure_context)) failure_context = "root average-down"
    call average_down_three_level_reactive_eb_state_2d( &
      species, root_refluxed, root_refluxed_temperature, root_geometry, &
      level_one_refluxed, level_one_refluxed_temperature, &
      level_one_geometry, root_patch, level_two_candidate, &
      level_two_candidate_temperature, level_two_geometry, level_one_patch, &
      root_candidate, root_temperature_work, level_one_candidate, &
      level_one_candidate_temperature, local_ok)
    if (.not. local_ok) return
    if (.not. level_two_interface_is_regular(level_one_geometry)) then
      allocate(root_closed, mold=root_state)
      allocate(root_closed_temperature, mold=root_temperature)
      if (present(failure_context)) &
        failure_context = "root conservation closure"
      call close_cut_interface_conservation_2d( &
        species, root_integral_before, root_candidate, &
        root_temperature_work, root_geometry, level_one_candidate, &
        level_one_geometry, root_patch, root_x_flux, root_y_flux, dt, &
        root_closed, root_closed_temperature, local_ok)
      if (.not. local_ok) return
      root_candidate = root_closed
      root_temperature_work = root_closed_temperature
    end if

    new_root_state = root_candidate
    new_root_temperature = root_temperature_work
    new_level_one_state = level_one_candidate
    new_level_one_temperature = level_one_candidate_temperature
    new_level_two_state = level_two_candidate
    new_level_two_temperature = level_two_candidate_temperature
    ok = .true.
    if (present(failure_context)) failure_context = "none"
  end subroutine advance_three_level_reactive_eb_transport_euler_2d

  subroutine advance_three_level_reactive_eb_transport_2d( &
      species, transport, root_state, root_temperature, root_geometry, &
      level_one_state, level_one_temperature, level_one_geometry, root_patch, &
      level_two_state, level_two_temperature, level_two_geometry, &
      level_one_patch, interval, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, target_volume_fraction, max_order, &
      new_root_state, new_root_temperature, new_level_one_state, &
      new_level_one_temperature, new_level_two_state, &
      new_level_two_temperature, minimum_theta, ok, failure_context)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: root_state(:, :, :), root_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: root_geometry
    real(dp), intent(in) :: level_one_state(:, :, :)
    real(dp), intent(in) :: level_one_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: level_one_geometry
    type(amr_eb_patch_2d), intent(in) :: root_patch
    real(dp), intent(in) :: level_two_state(:, :, :)
    real(dp), intent(in) :: level_two_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: level_two_geometry
    type(amr_eb_patch_2d), intent(in) :: level_one_patch
    real(dp), intent(in) :: interval, target_volume_fraction
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    integer, intent(in) :: max_order
    real(dp), intent(out) :: new_root_state(:, :, :)
    real(dp), intent(out) :: new_root_temperature(:, :)
    real(dp), intent(out) :: new_level_one_state(:, :, :)
    real(dp), intent(out) :: new_level_one_temperature(:, :)
    real(dp), intent(out) :: new_level_two_state(:, :, :)
    real(dp), intent(out) :: new_level_two_temperature(:, :)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    character(len=*), intent(out), optional :: failure_context

    real(dp), allocatable :: stage_root(:, :, :), stage_root_temperature(:, :)
    real(dp), allocatable :: stage_one(:, :, :), stage_one_temperature(:, :)
    real(dp), allocatable :: stage_two(:, :, :), stage_two_temperature(:, :)
    real(dp), allocatable :: euler_root(:, :, :), euler_root_temperature(:, :)
    real(dp), allocatable :: euler_one(:, :, :), euler_one_temperature(:, :)
    real(dp), allocatable :: euler_two(:, :, :), euler_two_temperature(:, :)
    real(dp), allocatable :: candidate_root(:, :, :)
    real(dp), allocatable :: candidate_root_temperature(:, :)
    real(dp), allocatable :: candidate_one(:, :, :)
    real(dp), allocatable :: candidate_one_temperature(:, :)
    real(dp), allocatable :: candidate_two(:, :, :)
    real(dp), allocatable :: candidate_two_temperature(:, :)
    real(dp) :: theta_one, theta_two
    logical :: local_ok

    new_root_state = root_state
    new_root_temperature = root_temperature
    new_level_one_state = level_one_state
    new_level_one_temperature = level_one_temperature
    new_level_two_state = level_two_state
    new_level_two_temperature = level_two_temperature
    minimum_theta = 1.0_dp
    ok = .false.
    if (present(failure_context)) failure_context = "input validation"
    if (.not. ieee_is_finite(interval) .or. interval < 0.0_dp) return
    if (interval <= tiny(1.0_dp) .or. .not. (viscosity_enabled .or. &
        thermal_conduction_enabled .or. species_diffusion_enabled)) then
      ok = .true.
      return
    end if

    allocate(stage_root, mold=root_state)
    allocate(stage_root_temperature, mold=root_temperature)
    allocate(stage_one, mold=level_one_state)
    allocate(stage_one_temperature, mold=level_one_temperature)
    allocate(stage_two, mold=level_two_state)
    allocate(stage_two_temperature, mold=level_two_temperature)
    if (present(failure_context)) failure_context = "first Euler stage"
    call advance_three_level_reactive_eb_transport_euler_2d( &
      species, transport, root_state, root_temperature, root_geometry, &
      level_one_state, level_one_temperature, level_one_geometry, root_patch, &
      level_two_state, level_two_temperature, level_two_geometry, &
      level_one_patch, interval, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, target_volume_fraction, max_order, &
      stage_root, stage_root_temperature, stage_one, stage_one_temperature, &
      stage_two, stage_two_temperature, theta_one, local_ok, failure_context)
    if (.not. local_ok) return

    allocate(euler_root, mold=root_state)
    allocate(euler_root_temperature, mold=root_temperature)
    allocate(euler_one, mold=level_one_state)
    allocate(euler_one_temperature, mold=level_one_temperature)
    allocate(euler_two, mold=level_two_state)
    allocate(euler_two_temperature, mold=level_two_temperature)
    if (present(failure_context)) failure_context = "second Euler stage"
    call advance_three_level_reactive_eb_transport_euler_2d( &
      species, transport, stage_root, stage_root_temperature, root_geometry, &
      stage_one, stage_one_temperature, level_one_geometry, root_patch, &
      stage_two, stage_two_temperature, level_two_geometry, level_one_patch, &
      interval, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, max_order, euler_root, euler_root_temperature, &
      euler_one, euler_one_temperature, euler_two, euler_two_temperature, &
      theta_two, local_ok, failure_context)
    if (.not. local_ok) return

    allocate(candidate_root, source=0.5_dp * (root_state + euler_root))
    allocate(candidate_one, source=0.5_dp * (level_one_state + euler_one))
    allocate(candidate_two, source=0.5_dp * (level_two_state + euler_two))
    allocate(candidate_root_temperature, mold=root_temperature)
    allocate(candidate_one_temperature, mold=level_one_temperature)
    allocate(candidate_two_temperature, mold=level_two_temperature)
    if (present(failure_context)) failure_context = "root blend recovery"
    call recover_transport_temperature_2d( &
      species, candidate_root, &
      0.5_dp * (root_temperature + euler_root_temperature), root_geometry, &
      candidate_root_temperature, local_ok)
    if (.not. local_ok) return
    if (present(failure_context)) failure_context = "level-one blend recovery"
    call recover_transport_temperature_2d( &
      species, candidate_one, &
      0.5_dp * (level_one_temperature + euler_one_temperature), &
      level_one_geometry, candidate_one_temperature, local_ok)
    if (.not. local_ok) return
    if (present(failure_context)) failure_context = "level-two blend recovery"
    call recover_transport_temperature_2d( &
      species, candidate_two, &
      0.5_dp * (level_two_temperature + euler_two_temperature), &
      level_two_geometry, candidate_two_temperature, local_ok)
    if (.not. local_ok) return
    if (present(failure_context)) failure_context = "final average-down"
    call average_down_three_level_reactive_eb_state_2d( &
      species, candidate_root, candidate_root_temperature, root_geometry, &
      candidate_one, candidate_one_temperature, level_one_geometry, &
      root_patch, candidate_two, candidate_two_temperature, &
      level_two_geometry, level_one_patch, new_root_state, &
      new_root_temperature, new_level_one_state, new_level_one_temperature, &
      local_ok)
    if (.not. local_ok) return
    new_level_two_state = candidate_two
    new_level_two_temperature = candidate_two_temperature
    minimum_theta = min(theta_one, theta_two)
    ok = .true.
    if (present(failure_context)) failure_context = "none"
  end subroutine advance_three_level_reactive_eb_transport_2d

end module amr_eb_multilevel_transport_2d_mod
