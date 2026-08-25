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
    compute_reactive_eb_cfl_timestep_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, build_amr_eb_patch_2d, composite_eb_integral_2d
  use amr_eb_reactive_2d_mod, only: &
    prolong_reactive_eb_patch_pcm_2d, &
    advance_two_level_reactive_eb_hydro_2d
  implicit none
  private

  public :: compute_reactive_eb_amr_cfl_timestep_2d
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
      config%refinement_ratio >= 2
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

  subroutine simulate_reactive_eb_amr_2d( &
      species, config, coarse_state, coarse_temperature, coarse_geometry, &
      fine_state, fine_temperature, fine_geometry, patch, time, steps, &
      initial_integrals, final_integrals, minimum_dt, base_density, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_eb_amr_2d_config), intent(in) :: config
    real(dp), allocatable, intent(out) :: coarse_state(:, :, :)
    real(dp), allocatable, intent(out) :: coarse_temperature(:, :)
    type(eb_geometry_2d), intent(out) :: coarse_geometry
    real(dp), allocatable, intent(out) :: fine_state(:, :, :)
    real(dp), allocatable, intent(out) :: fine_temperature(:, :)
    type(eb_geometry_2d), intent(out) :: fine_geometry
    type(amr_eb_patch_2d), intent(out) :: patch
    real(dp), intent(out) :: time, minimum_dt, base_density
    integer, intent(out) :: steps
    real(dp), allocatable, intent(out) :: initial_integrals(:)
    real(dp), allocatable, intent(out) :: final_integrals(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: coarse_candidate(:, :, :)
    real(dp), allocatable :: coarse_candidate_temperature(:, :)
    real(dp), allocatable :: fine_candidate(:, :, :)
    real(dp), allocatable :: fine_candidate_temperature(:, :)
    real(dp) :: coarse_dx, coarse_dy, fine_x_lower, fine_x_upper
    real(dp) :: fine_y_lower, fine_y_upper, dt, remaining, time_tolerance
    logical :: local_ok
    integer :: fine_nx, fine_ny, nvar

    ok = .false.
    time = 0.0_dp
    steps = 0
    minimum_dt = 0.0_dp
    base_density = 0.0_dp
    if (.not. supported_reactive_eb_amr_config(config)) return
    call build_configured_eb_geometry_2d( &
      config%eb, coarse_geometry, local_ok)
    if (.not. local_ok) return
    fine_x_lower = coarse_geometry%x_lower + &
      real(config%coarse_i_lower - 1, dp) * coarse_geometry%dx
    fine_x_upper = coarse_geometry%x_lower + &
      real(config%coarse_i_upper, dp) * coarse_geometry%dx
    fine_y_lower = coarse_geometry%y_lower + &
      real(config%coarse_j_lower - 1, dp) * coarse_geometry%dy
    fine_y_upper = coarse_geometry%y_lower + &
      real(config%coarse_j_upper, dp) * coarse_geometry%dy
    fine_nx = (config%coarse_i_upper - config%coarse_i_lower + 1) * &
      config%refinement_ratio
    fine_ny = (config%coarse_j_upper - config%coarse_j_lower + 1) * &
      config%refinement_ratio
    call build_configured_eb_geometry_region_2d( &
      config%eb, fine_nx, fine_ny, fine_x_lower, fine_x_upper, &
      fine_y_lower, fine_y_upper, fine_geometry, local_ok)
    if (.not. local_ok) return
    call build_amr_eb_patch_2d( &
      coarse_geometry, fine_geometry, config%coarse_i_lower, &
      config%coarse_i_upper, config%coarse_j_lower, &
      config%coarse_j_upper, config%refinement_ratio, patch, local_ok)
    if (.not. local_ok) return

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

    allocate(initial_integrals(nvar), final_integrals(nvar))
    call composite_eb_integral_2d( &
      coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
      initial_integrals, local_ok)
    if (.not. local_ok) return
    allocate(coarse_candidate, mold=coarse_state)
    allocate(coarse_candidate_temperature, mold=coarse_temperature)
    allocate(fine_candidate, mold=fine_state)
    allocate(fine_candidate_temperature, mold=fine_temperature)
    minimum_dt = huge(1.0_dp)
    time_tolerance = 16.0_dp * epsilon(1.0_dp) * &
      max(tiny(1.0_dp), abs(config%eb%flow%final_time))

    do
      remaining = config%eb%flow%final_time - time
      if (remaining <= time_tolerance) exit
      if (steps >= config%eb%flow%maximum_steps) return
      call compute_reactive_eb_amr_cfl_timestep_2d( &
        species, coarse_state, coarse_temperature, coarse_geometry, &
        fine_state, fine_temperature, fine_geometry, &
        config%refinement_ratio, config%eb%flow%cfl, dt, local_ok)
      if (.not. local_ok) return
      dt = min(dt, remaining)
      call advance_two_level_reactive_eb_hydro_2d( &
        species, coarse_state, coarse_temperature, coarse_geometry, &
        fine_state, fine_temperature, fine_geometry, patch, &
        config%eb%flow%riemann_solver, config%eb%flow%reconstruction, &
        config%eb%flow%limiter, config%eb%state_redist_max_order, dt, &
        coarse_candidate, coarse_candidate_temperature, fine_candidate, &
        fine_candidate_temperature, local_ok, &
        config%eb%state_redist_target_volume_fraction)
      if (.not. local_ok) return
      coarse_state = coarse_candidate
      coarse_temperature = coarse_candidate_temperature
      fine_state = fine_candidate
      fine_temperature = fine_candidate_temperature
      time = time + dt
      minimum_dt = min(minimum_dt, dt)
      steps = steps + 1
    end do
    time = config%eb%flow%final_time
    call composite_eb_integral_2d( &
      coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
      final_integrals, local_ok)
    if (.not. local_ok) return
    ok = steps > 0 .and. ieee_is_finite(minimum_dt) .and. minimum_dt > 0.0_dp
  end subroutine simulate_reactive_eb_amr_2d

end module reactive_eb_amr_2d_driver_mod
