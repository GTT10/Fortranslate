module reactive_eb_amr_2d_driver_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_2d_mod, only: initialize_reactive_2d
  use eb_geometry_2d_mod, only: eb_geometry_2d
  use simulation_config_reactive_eb_amr_2d_mod, only: &
    reactive_eb_amr_2d_config
  use reactive_eb_2d_driver_mod, only: &
    build_configured_eb_geometry_2d, &
    build_configured_eb_geometry_region_2d, &
    compute_reactive_eb_cfl_timestep_2d, reactive_eb_integrals_2d
  use eb_reactive_hydro_2d_mod, only: advance_reactive_eb_hydro_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, build_amr_eb_patch_2d, composite_eb_integral_2d
  use amr_eb_reactive_2d_mod, only: &
    prolong_reactive_eb_patch_pcm_2d, &
    advance_two_level_reactive_eb_hydro_2d
  use amr_eb_regrid_2d_mod, only: &
    amr_eb_tagging_criteria_2d, amr_eb_regrid_plan_2d, &
    plan_reactive_eb_temperature_regrid_2d, &
    collapse_two_level_reactive_eb_patch_2d, &
    regrid_two_level_reactive_eb_patch_2d
  implicit none
  private

  public :: compute_reactive_eb_amr_cfl_timestep_2d
  public :: regrid_reactive_eb_amr_hierarchy_2d
  public :: simulate_reactive_eb_amr_2d

contains

  pure logical function supported_reactive_eb_amr_config(config) &
      result(supported)
    type(reactive_eb_amr_2d_config), intent(in) :: config

    supported = config%eb%flow%nx >= 4 .and. &
      config%eb%flow%ny >= 4 .and. &
      config%eb%flow%maximum_steps >= 1 .and. &
      config%eb%flow%x_upper > config%eb%flow%x_lower .and. &
      config%eb%flow%y_upper > config%eb%flow%y_lower .and. &
      ieee_is_finite(config%eb%flow%final_time) .and. &
      config%eb%flow%final_time > 0.0_dp .and. &
      ieee_is_finite(config%eb%flow%cfl) .and. &
      config%eb%flow%cfl > 0.0_dp .and. config%eb%flow%cfl <= 0.8_dp .and. &
      .not. config%eb%flow%chemistry_enabled .and. &
      .not. config%eb%flow%transport_enabled .and. &
      .not. config%eb%flow%use_transverse_correction .and. &
      (trim(config%eb%flow%reconstruction) == "pcm" .or. &
       trim(config%eb%flow%reconstruction) == "characteristic_plm") .and. &
      trim(config%eb%flow%boundary_x_lower) == "outflow" .and. &
      trim(config%eb%flow%boundary_x_upper) == "outflow" .and. &
      trim(config%eb%flow%boundary_y_lower) == "outflow" .and. &
      trim(config%eb%flow%boundary_y_upper) == "outflow" .and. &
      ieee_is_finite(config%eb%state_redist_target_volume_fraction) .and. &
      config%eb%state_redist_target_volume_fraction > 0.0_dp .and. &
      config%eb%state_redist_target_volume_fraction <= 1.0_dp .and. &
      (config%eb%state_redist_max_order == 0 .or. &
       config%eb%state_redist_max_order == 2) .and. &
      config%coarse_i_lower > 1 .and. &
      config%coarse_i_upper < config%eb%flow%nx .and. &
      config%coarse_j_lower > 1 .and. &
      config%coarse_j_upper < config%eb%flow%ny .and. &
      config%coarse_i_upper >= config%coarse_i_lower .and. &
      config%coarse_j_upper >= config%coarse_j_lower .and. &
      config%refinement_ratio >= 2 .and. config%regrid_interval >= 1 .and. &
      (.not. config%remove_fine_patch_when_untagged .or. &
       config%dynamic_regridding) .and. &
      ieee_is_finite(config%regrid_relative_temperature_gradient) .and. &
      config%regrid_relative_temperature_gradient >= 0.0_dp .and. &
      ieee_is_finite(config%regrid_absolute_temperature_gradient) .and. &
      config%regrid_absolute_temperature_gradient >= 0.0_dp .and. &
      ieee_is_finite(config%regrid_temperature_scale_floor) .and. &
      config%regrid_temperature_scale_floor > 0.0_dp .and. &
      config%regrid_buffer_cells >= 0 .and. &
      config%regrid_minimum_patch_cells_x >= 1 .and. &
      config%regrid_minimum_patch_cells_x <= config%eb%flow%nx - 2 .and. &
      config%regrid_minimum_patch_cells_y >= 1 .and. &
      config%regrid_minimum_patch_cells_y <= config%eb%flow%ny - 2
  end function supported_reactive_eb_amr_config

  subroutine compute_reactive_eb_amr_cfl_timestep_2d( &
      species, coarse_state, coarse_temperature, coarse_geometry, &
      fine_state, fine_temperature, fine_geometry, refinement_ratio, &
      cfl, dt, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: coarse_state(:, :, :), coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    real(dp), intent(in) :: fine_state(:, :, :), fine_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: fine_geometry
    integer, intent(in) :: refinement_ratio
    real(dp), intent(in) :: cfl
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok

    real(dp) :: coarse_dt, fine_dt
    logical :: local_ok

    dt = 0.0_dp
    ok = .false.
    if (refinement_ratio < 2) return
    call compute_reactive_eb_cfl_timestep_2d( &
      species, coarse_state, coarse_temperature, coarse_geometry, &
      cfl, coarse_dt, local_ok)
    if (.not. local_ok) return
    call compute_reactive_eb_cfl_timestep_2d( &
      species, fine_state, fine_temperature, fine_geometry, &
      cfl, fine_dt, local_ok)
    if (.not. local_ok) return
    dt = min(coarse_dt, real(refinement_ratio, dp) * fine_dt)
    ok = ieee_is_finite(dt) .and. dt > 0.0_dp
  end subroutine compute_reactive_eb_amr_cfl_timestep_2d

  subroutine compute_reactive_eb_amr_integrals_2d( &
      coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
      fine_active, integrals, ok)
    real(dp), intent(in) :: coarse_state(:, :, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    real(dp), allocatable, intent(in) :: fine_state(:, :, :)
    type(eb_geometry_2d), intent(in) :: fine_geometry
    type(amr_eb_patch_2d), intent(in) :: patch
    logical, intent(in) :: fine_active
    real(dp), intent(out) :: integrals(:)
    logical, intent(out) :: ok

    if (fine_active) then
      ok = allocated(fine_state)
      if (.not. ok) return
      call composite_eb_integral_2d( &
        coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
        integrals, ok)
      return
    end if
    if (allocated(fine_state)) then
      integrals = 0.0_dp
      ok = .false.
      return
    end if
    call reactive_eb_integrals_2d( &
      coarse_state, coarse_geometry, integrals, ok)
  end subroutine compute_reactive_eb_amr_integrals_2d

  subroutine build_reactive_eb_amr_patch_geometry_2d( &
      config, coarse_geometry, coarse_i_lower, coarse_i_upper, &
      coarse_j_lower, coarse_j_upper, fine_geometry, patch, ok)
    type(reactive_eb_amr_2d_config), intent(in) :: config
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    integer, intent(in) :: coarse_i_lower, coarse_i_upper
    integer, intent(in) :: coarse_j_lower, coarse_j_upper
    type(eb_geometry_2d), intent(out) :: fine_geometry
    type(amr_eb_patch_2d), intent(out) :: patch
    logical, intent(out) :: ok

    real(dp) :: fine_x_lower, fine_x_upper, fine_y_lower, fine_y_upper
    integer :: fine_nx, fine_ny

    ok = .false.
    if (coarse_i_lower <= 1 .or. &
        coarse_i_upper >= coarse_geometry%nx .or. &
        coarse_j_lower <= 1 .or. &
        coarse_j_upper >= coarse_geometry%ny .or. &
        coarse_i_upper < coarse_i_lower .or. &
        coarse_j_upper < coarse_j_lower) return
    fine_x_lower = coarse_geometry%x_lower + &
      real(coarse_i_lower - 1, dp) * coarse_geometry%dx
    fine_x_upper = coarse_geometry%x_lower + &
      real(coarse_i_upper, dp) * coarse_geometry%dx
    fine_y_lower = coarse_geometry%y_lower + &
      real(coarse_j_lower - 1, dp) * coarse_geometry%dy
    fine_y_upper = coarse_geometry%y_lower + &
      real(coarse_j_upper, dp) * coarse_geometry%dy
    fine_nx = (coarse_i_upper - coarse_i_lower + 1) * &
      config%refinement_ratio
    fine_ny = (coarse_j_upper - coarse_j_lower + 1) * &
      config%refinement_ratio
    call build_configured_eb_geometry_region_2d( &
      config%eb, fine_nx, fine_ny, fine_x_lower, fine_x_upper, &
      fine_y_lower, fine_y_upper, fine_geometry, ok)
    if (.not. ok) return
    call build_amr_eb_patch_2d( &
      coarse_geometry, fine_geometry, coarse_i_lower, coarse_i_upper, &
      coarse_j_lower, coarse_j_upper, config%refinement_ratio, patch, ok)
  end subroutine build_reactive_eb_amr_patch_geometry_2d

  subroutine regrid_reactive_eb_amr_hierarchy_2d( &
      species, config, coarse_state, coarse_temperature, coarse_geometry, &
      fine_state, fine_temperature, fine_geometry, patch, fine_active, &
      changed, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_eb_amr_2d_config), intent(in) :: config
    real(dp), allocatable, intent(inout) :: coarse_state(:, :, :)
    real(dp), allocatable, intent(inout) :: coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    real(dp), allocatable, intent(inout) :: fine_state(:, :, :)
    real(dp), allocatable, intent(inout) :: fine_temperature(:, :)
    type(eb_geometry_2d), intent(inout) :: fine_geometry
    type(amr_eb_patch_2d), intent(inout) :: patch
    logical, intent(inout) :: fine_active
    logical, intent(out) :: changed, ok

    type(amr_eb_tagging_criteria_2d) :: criteria
    type(amr_eb_regrid_plan_2d) :: plan
    type(eb_geometry_2d) :: candidate_fine_geometry
    type(amr_eb_patch_2d) :: candidate_patch
    real(dp), allocatable :: candidate_coarse_state(:, :, :)
    real(dp), allocatable :: candidate_coarse_temperature(:, :)
    real(dp), allocatable :: candidate_fine_state(:, :, :)
    real(dp), allocatable :: candidate_fine_temperature(:, :)
    logical, allocatable :: tags(:, :)
    logical :: local_ok
    integer :: nvar

    changed = .false.
    ok = .false.
    if (fine_active) then
      if (.not. allocated(fine_state) .or. &
          .not. allocated(fine_temperature) .or. &
          .not. patch%is_valid(coarse_geometry, fine_geometry)) return
    else
      if (allocated(fine_state) .or. allocated(fine_temperature) .or. &
          fine_geometry%is_valid() .or. patch%refinement_ratio /= 0) return
    end if
    criteria%relative_gradient_threshold = &
      config%regrid_relative_temperature_gradient
    criteria%absolute_gradient_threshold = &
      config%regrid_absolute_temperature_gradient
    criteria%scale_floor = config%regrid_temperature_scale_floor
    criteria%buffer_cells = config%regrid_buffer_cells
    criteria%minimum_patch_cells_x = config%regrid_minimum_patch_cells_x
    criteria%minimum_patch_cells_y = config%regrid_minimum_patch_cells_y
    allocate(tags(coarse_geometry%nx, coarse_geometry%ny))
    call plan_reactive_eb_temperature_regrid_2d( &
      coarse_temperature, coarse_geometry, criteria, tags, plan, local_ok)
    if (.not. local_ok) return
    if (.not. plan%active) then
      if (fine_active .and. config%remove_fine_patch_when_untagged) then
        allocate(candidate_coarse_state, mold=coarse_state)
        allocate(candidate_coarse_temperature, mold=coarse_temperature)
        call collapse_two_level_reactive_eb_patch_2d( &
          species, coarse_state, coarse_temperature, coarse_geometry, &
          fine_state, fine_geometry, patch, candidate_coarse_state, &
          candidate_coarse_temperature, local_ok)
        if (.not. local_ok) return
        call move_alloc(candidate_coarse_state, coarse_state)
        call move_alloc(candidate_coarse_temperature, coarse_temperature)
        deallocate(fine_state, fine_temperature)
        fine_geometry = eb_geometry_2d()
        patch = amr_eb_patch_2d()
        fine_active = .false.
        changed = .true.
      end if
      ok = .true.
      return
    end if
    if (fine_active .and. &
        plan%coarse_i_lower == patch%coarse_i_lower .and. &
        plan%coarse_i_upper == patch%coarse_i_upper .and. &
        plan%coarse_j_lower == patch%coarse_j_lower .and. &
        plan%coarse_j_upper == patch%coarse_j_upper) then
      ok = .true.
      return
    end if

    call build_reactive_eb_amr_patch_geometry_2d( &
      config, coarse_geometry, plan%coarse_i_lower, &
      plan%coarse_i_upper, plan%coarse_j_lower, plan%coarse_j_upper, &
      candidate_fine_geometry, candidate_patch, local_ok)
    if (.not. local_ok) return
    nvar = size(coarse_state, 1)
    allocate(candidate_coarse_state, mold=coarse_state)
    allocate(candidate_coarse_temperature, mold=coarse_temperature)
    allocate(candidate_fine_state( &
      nvar, candidate_fine_geometry%nx, candidate_fine_geometry%ny))
    allocate(candidate_fine_temperature( &
      candidate_fine_geometry%nx, candidate_fine_geometry%ny))
    if (fine_active) then
      call regrid_two_level_reactive_eb_patch_2d( &
        species, coarse_state, coarse_temperature, coarse_geometry, &
        fine_state, fine_temperature, fine_geometry, patch, &
        candidate_fine_geometry, candidate_patch, candidate_coarse_state, &
        candidate_coarse_temperature, candidate_fine_state, &
        candidate_fine_temperature, local_ok)
    else
      candidate_coarse_state = coarse_state
      candidate_coarse_temperature = coarse_temperature
      call prolong_reactive_eb_patch_pcm_2d( &
        species, coarse_state, coarse_temperature, coarse_geometry, &
        candidate_fine_geometry, candidate_patch, candidate_fine_state, &
        candidate_fine_temperature, local_ok)
    end if
    if (.not. local_ok) return
    call move_alloc(candidate_coarse_state, coarse_state)
    call move_alloc(candidate_coarse_temperature, coarse_temperature)
    call move_alloc(candidate_fine_state, fine_state)
    call move_alloc(candidate_fine_temperature, fine_temperature)
    fine_geometry = candidate_fine_geometry
    patch = candidate_patch
    fine_active = .true.
    changed = .true.
    ok = .true.
  end subroutine regrid_reactive_eb_amr_hierarchy_2d

  subroutine simulate_reactive_eb_amr_2d( &
      species, config, coarse_state, coarse_temperature, coarse_geometry, &
      fine_state, fine_temperature, fine_geometry, patch, fine_active, time, &
      steps, regrids, initial_integrals, final_integrals, minimum_dt, &
      base_density, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_eb_amr_2d_config), intent(in) :: config
    real(dp), allocatable, intent(out) :: coarse_state(:, :, :)
    real(dp), allocatable, intent(out) :: coarse_temperature(:, :)
    type(eb_geometry_2d), intent(out) :: coarse_geometry
    real(dp), allocatable, intent(out) :: fine_state(:, :, :)
    real(dp), allocatable, intent(out) :: fine_temperature(:, :)
    type(eb_geometry_2d), intent(out) :: fine_geometry
    type(amr_eb_patch_2d), intent(out) :: patch
    logical, intent(out) :: fine_active
    real(dp), intent(out) :: time, minimum_dt, base_density
    integer, intent(out) :: steps, regrids
    real(dp), allocatable, intent(out) :: initial_integrals(:)
    real(dp), allocatable, intent(out) :: final_integrals(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: coarse_candidate(:, :, :)
    real(dp), allocatable :: coarse_candidate_temperature(:, :)
    real(dp), allocatable :: fine_candidate(:, :, :)
    real(dp), allocatable :: fine_candidate_temperature(:, :)
    real(dp) :: coarse_dx, coarse_dy, dt, remaining, time_tolerance
    logical :: changed, local_ok
    integer :: fine_nx, fine_ny, nvar

    ok = .false.
    time = 0.0_dp
    steps = 0
    regrids = 0
    fine_active = .false.
    minimum_dt = 0.0_dp
    base_density = 0.0_dp
    if (.not. supported_reactive_eb_amr_config(config)) return
    call build_configured_eb_geometry_2d( &
      config%eb, coarse_geometry, local_ok)
    if (.not. local_ok) return
    call build_reactive_eb_amr_patch_geometry_2d( &
      config, coarse_geometry, config%coarse_i_lower, &
      config%coarse_i_upper, config%coarse_j_lower, &
      config%coarse_j_upper, fine_geometry, patch, local_ok)
    if (.not. local_ok) return
    fine_nx = fine_geometry%nx
    fine_ny = fine_geometry%ny

    call initialize_reactive_2d( &
      species, config%eb%flow, coarse_state, coarse_temperature, &
      coarse_dx, coarse_dy, base_density, local_ok)
    if (.not. local_ok) return
    if (abs(coarse_dx - coarse_geometry%dx) > &
        8.0_dp * epsilon(1.0_dp) * coarse_geometry%dx .or. &
        abs(coarse_dy - coarse_geometry%dy) > &
        8.0_dp * epsilon(1.0_dp) * coarse_geometry%dy) return
    nvar = size(coarse_state, 1)
    allocate(fine_state(nvar, fine_nx, fine_ny))
    allocate(fine_temperature(fine_nx, fine_ny))
    call prolong_reactive_eb_patch_pcm_2d( &
      species, coarse_state, coarse_temperature, coarse_geometry, &
      fine_geometry, patch, fine_state, fine_temperature, local_ok)
    if (.not. local_ok) return
    fine_active = .true.
    if (config%dynamic_regridding .and. config%regrid_at_initialization) then
      call regrid_reactive_eb_amr_hierarchy_2d( &
        species, config, coarse_state, coarse_temperature, coarse_geometry, &
        fine_state, fine_temperature, fine_geometry, patch, fine_active, &
        changed, local_ok)
      if (.not. local_ok) return
      if (changed) regrids = regrids + 1
    end if

    allocate(initial_integrals(nvar), final_integrals(nvar))
    call compute_reactive_eb_amr_integrals_2d( &
      coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
      fine_active, initial_integrals, local_ok)
    if (.not. local_ok) return
    minimum_dt = huge(1.0_dp)
    time_tolerance = 16.0_dp * epsilon(1.0_dp) * &
      max(tiny(1.0_dp), abs(config%eb%flow%final_time))

    do
      remaining = config%eb%flow%final_time - time
      if (remaining <= time_tolerance) exit
      if (steps >= config%eb%flow%maximum_steps) return
      if (allocated(coarse_candidate)) deallocate(coarse_candidate)
      if (allocated(coarse_candidate_temperature)) &
        deallocate(coarse_candidate_temperature)
      if (allocated(fine_candidate)) deallocate(fine_candidate)
      if (allocated(fine_candidate_temperature)) &
        deallocate(fine_candidate_temperature)
      allocate(coarse_candidate, mold=coarse_state)
      allocate(coarse_candidate_temperature, mold=coarse_temperature)
      if (fine_active) then
        allocate(fine_candidate, mold=fine_state)
        allocate(fine_candidate_temperature, mold=fine_temperature)
        call compute_reactive_eb_amr_cfl_timestep_2d( &
          species, coarse_state, coarse_temperature, coarse_geometry, &
          fine_state, fine_temperature, fine_geometry, &
          config%refinement_ratio, config%eb%flow%cfl, dt, local_ok)
      else
        call compute_reactive_eb_cfl_timestep_2d( &
          species, coarse_state, coarse_temperature, coarse_geometry, &
          config%eb%flow%cfl, dt, local_ok)
      end if
      if (.not. local_ok) return
      dt = min(dt, remaining)
      if (fine_active) then
        call advance_two_level_reactive_eb_hydro_2d( &
          species, coarse_state, coarse_temperature, coarse_geometry, &
          fine_state, fine_temperature, fine_geometry, patch, &
          config%eb%flow%riemann_solver, config%eb%flow%reconstruction, &
          config%eb%flow%limiter, config%eb%state_redist_max_order, dt, &
          coarse_candidate, coarse_candidate_temperature, fine_candidate, &
          fine_candidate_temperature, local_ok, &
          config%eb%state_redist_target_volume_fraction)
      else
        call advance_reactive_eb_hydro_2d( &
          species, coarse_state, coarse_temperature, coarse_geometry, &
          config%eb%flow%riemann_solver, dt, coarse_candidate, &
          coarse_candidate_temperature, local_ok, &
          config%eb%state_redist_target_volume_fraction, &
          config%eb%flow%reconstruction, config%eb%flow%limiter, &
          config%eb%state_redist_max_order)
      end if
      if (.not. local_ok) return
      coarse_state = coarse_candidate
      coarse_temperature = coarse_candidate_temperature
      if (fine_active) then
        fine_state = fine_candidate
        fine_temperature = fine_candidate_temperature
      end if
      time = time + dt
      minimum_dt = min(minimum_dt, dt)
      steps = steps + 1
      if (config%dynamic_regridding .and. &
          modulo(steps, config%regrid_interval) == 0) then
        call regrid_reactive_eb_amr_hierarchy_2d( &
          species, config, coarse_state, coarse_temperature, &
          coarse_geometry, fine_state, fine_temperature, fine_geometry, &
          patch, fine_active, changed, local_ok)
        if (.not. local_ok) return
        if (changed) regrids = regrids + 1
      end if
    end do
    time = config%eb%flow%final_time
    call compute_reactive_eb_amr_integrals_2d( &
      coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
      fine_active, final_integrals, local_ok)
    if (.not. local_ok) return
    ok = steps > 0 .and. ieee_is_finite(minimum_dt) .and. minimum_dt > 0.0_dp
  end subroutine simulate_reactive_eb_amr_2d

end module reactive_eb_amr_2d_driver_mod
