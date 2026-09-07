module amr_reactive_transport_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar, reactive_nprim
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use reactive_transport_3d_mod, only: &
    reactive_transport_fluxes_3d, reactive_transport_ghosted_fluxes_3d, &
    reactive_transport_interface_theta_3d, reactive_transport_timestep_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d, average_down_3d
  use amr_reactive_3d_mod, only: &
    advance_amr_reactive_chemistry_3d, advance_amr_reactive_hydro_3d, &
    advance_amr_reactive_strang_3d, accumulate_fine_interface_fluxes_3d, &
    reflux_coarse_3d, interpolate_coarse_fine_ghost_plm_3d
  implicit none
  private

  public :: compute_amr_reactive_transport_timestep_3d
  public :: advance_amr_reactive_transport_3d
  public :: advance_amr_reactive_full_3d

contains

  subroutine compute_amr_reactive_transport_timestep_3d( &
      species, transport, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, dx, dy, dz, transport_cfl, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, dt, maximum_diffusivity, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_state(:, :, :, :)
    real(dp), intent(in) :: coarse_temperature(:, :, :)
    real(dp), intent(in) :: fine_state(:, :, :, :)
    real(dp), intent(in) :: fine_temperature(:, :, :)
    real(dp), intent(in) :: dx, dy, dz, transport_cfl
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled
    real(dp), intent(out) :: dt, maximum_diffusivity
    logical, intent(out) :: ok

    real(dp) :: coarse_dt, fine_dt, coarse_diffusivity, fine_diffusivity
    real(dp) :: ratio
    logical :: local_ok

    dt = 0.0_dp
    maximum_diffusivity = 0.0_dp
    ok = .false.
    if (.not. valid_amr_transport_shapes_3d( &
          species, transport, patch, coarse_state, coarse_temperature, &
          fine_state, fine_temperature)) return
    if (.not. all(ieee_is_finite([dx, dy, dz, transport_cfl]))) return
    if (dx <= 0.0_dp .or. dy <= 0.0_dp .or. dz <= 0.0_dp) return

    call reactive_transport_timestep_3d( &
      species, transport, coarse_state, coarse_temperature, &
      patch%coarse_nx, patch%coarse_ny, patch%coarse_nz, dx, dy, dz, &
      transport_cfl, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, coarse_dt, coarse_diffusivity, local_ok)
    if (.not. local_ok) return
    ratio = real(patch%refinement_ratio, dp)
    call reactive_transport_timestep_3d( &
      species, transport, fine_state, fine_temperature, &
      patch%fine_nx(), patch%fine_ny(), patch%fine_nz(), &
      dx / ratio, dy / ratio, dz / ratio, transport_cfl, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, fine_dt, fine_diffusivity, local_ok)
    if (.not. local_ok) return

    dt = min(coarse_dt, ratio * ratio * fine_dt)
    maximum_diffusivity = max(coarse_diffusivity, fine_diffusivity)
    if (.not. ieee_is_finite(dt)) return
    if (.not. ieee_is_finite(maximum_diffusivity)) return
    ok = dt > 0.0_dp
  end subroutine compute_amr_reactive_transport_timestep_3d

  subroutine advance_amr_reactive_transport_3d( &
      species, transport, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, dx, dy, dz, interval, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, limiter, &
      minimum_theta, maximum_reflux_correction, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(inout) :: coarse_state(:, :, :, :)
    real(dp), intent(inout) :: coarse_temperature(:, :, :)
    real(dp), intent(inout) :: fine_state(:, :, :, :)
    real(dp), intent(inout) :: fine_temperature(:, :, :)
    real(dp), intent(in) :: dx, dy, dz, interval
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    character(len=*), intent(in) :: limiter
    real(dp), intent(out) :: minimum_theta, maximum_reflux_correction
    logical, intent(out) :: ok

    real(dp), allocatable :: stage_coarse(:, :, :, :)
    real(dp), allocatable :: stage_coarse_temperature(:, :, :)
    real(dp), allocatable :: stage_fine(:, :, :, :)
    real(dp), allocatable :: stage_fine_temperature(:, :, :)
    real(dp), allocatable :: euler_coarse(:, :, :, :)
    real(dp), allocatable :: euler_coarse_temperature(:, :, :)
    real(dp), allocatable :: euler_fine(:, :, :, :)
    real(dp), allocatable :: euler_fine_temperature(:, :, :)
    real(dp), allocatable :: candidate_coarse(:, :, :, :)
    real(dp), allocatable :: candidate_coarse_temperature(:, :, :)
    real(dp), allocatable :: candidate_fine(:, :, :, :)
    real(dp), allocatable :: candidate_fine_temperature(:, :, :)
    real(dp), allocatable :: recovered_temperature(:, :, :)
    real(dp) :: theta_one, theta_two, reflux_one, reflux_two
    logical :: local_ok

    minimum_theta = 1.0_dp
    maximum_reflux_correction = 0.0_dp
    ok = .false.
    if (.not. valid_amr_transport_shapes_3d( &
          species, transport, patch, coarse_state, coarse_temperature, &
          fine_state, fine_temperature)) return
    if (.not. all(ieee_is_finite([dx, dy, dz, interval]))) return
    if (interval < 0.0_dp .or. &
        dx <= 0.0_dp .or. dy <= 0.0_dp .or. dz <= 0.0_dp .or. &
        (barodiffusion_enabled .and. .not. species_diffusion_enabled) .or. &
        (trim(limiter) /= "minmod" .and. trim(limiter) /= "mc")) return
    if (interval <= tiny(1.0_dp) .or. .not. (viscosity_enabled .or. &
        thermal_conduction_enabled .or. species_diffusion_enabled)) then
      ok = .true.
      return
    end if

    allocate(stage_coarse, mold=coarse_state)
    allocate(stage_coarse_temperature, mold=coarse_temperature)
    allocate(stage_fine, mold=fine_state)
    allocate(stage_fine_temperature, mold=fine_temperature)
    call advance_amr_transport_euler_3d( &
      species, transport, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, dx, dy, dz, interval, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, limiter, &
      stage_coarse, stage_coarse_temperature, stage_fine, &
      stage_fine_temperature, theta_one, reflux_one, local_ok)
    if (.not. local_ok) return

    allocate(euler_coarse, mold=coarse_state)
    allocate(euler_coarse_temperature, mold=coarse_temperature)
    allocate(euler_fine, mold=fine_state)
    allocate(euler_fine_temperature, mold=fine_temperature)
    call advance_amr_transport_euler_3d( &
      species, transport, patch, stage_coarse, stage_coarse_temperature, &
      stage_fine, stage_fine_temperature, dx, dy, dz, interval, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, limiter, &
      euler_coarse, euler_coarse_temperature, euler_fine, &
      euler_fine_temperature, theta_two, reflux_two, local_ok)
    if (.not. local_ok) return

    allocate(candidate_coarse, source=0.5_dp * (coarse_state + euler_coarse))
    allocate(candidate_fine, source=0.5_dp * (fine_state + euler_fine))
    allocate(candidate_coarse_temperature, mold=coarse_temperature)
    allocate(candidate_fine_temperature, mold=fine_temperature)
    call recover_reactive_temperatures_3d( &
      species, candidate_fine, &
      0.5_dp * (fine_temperature + euler_fine_temperature), &
      patch%fine_nx(), patch%fine_ny(), patch%fine_nz(), &
      candidate_fine_temperature, local_ok)
    if (.not. local_ok) return
    call average_down_3d(candidate_coarse, candidate_fine, patch, local_ok)
    if (.not. local_ok) return
    allocate(recovered_temperature, mold=coarse_temperature)
    call recover_reactive_temperatures_3d( &
      species, candidate_coarse, &
      0.5_dp * (coarse_temperature + euler_coarse_temperature), &
      patch%coarse_nx, patch%coarse_ny, patch%coarse_nz, &
      recovered_temperature, local_ok)
    if (.not. local_ok) return
    candidate_coarse_temperature = recovered_temperature
    if (.not. all(ieee_is_finite(candidate_coarse)) .or. &
        .not. all(ieee_is_finite(candidate_coarse_temperature)) .or. &
        .not. all(ieee_is_finite(candidate_fine)) .or. &
        .not. all(ieee_is_finite(candidate_fine_temperature))) return

    coarse_state = candidate_coarse
    coarse_temperature = candidate_coarse_temperature
    fine_state = candidate_fine
    fine_temperature = candidate_fine_temperature
    minimum_theta = min(theta_one, theta_two)
    maximum_reflux_correction = max(reflux_one, reflux_two)
    ok = .true.
  end subroutine advance_amr_reactive_transport_3d

  subroutine advance_amr_reactive_full_3d( &
      species, reactions, transport, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, dx, dy, dz, dt, riemann_solver, &
      chemistry_enabled, rtol, atol, transport_enabled, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, minimum_transport_theta, &
      maximum_reflux_correction, ok, reconstruction, limiter, &
      chemistry_integrator)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(inout) :: coarse_state(:, :, :, :)
    real(dp), intent(inout) :: coarse_temperature(:, :, :)
    real(dp), intent(inout) :: fine_state(:, :, :, :)
    real(dp), intent(inout) :: fine_temperature(:, :, :)
    real(dp), intent(in) :: dx, dy, dz, dt, rtol, atol
    character(len=*), intent(in) :: riemann_solver
    logical, intent(in) :: chemistry_enabled, transport_enabled
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: minimum_transport_theta
    real(dp), intent(out) :: maximum_reflux_correction
    logical, intent(out) :: ok
    character(len=*), intent(in), optional :: reconstruction, limiter
    character(len=*), intent(in), optional :: chemistry_integrator

    real(dp), allocatable :: coarse_candidate(:, :, :, :)
    real(dp), allocatable :: coarse_candidate_temperature(:, :, :)
    real(dp), allocatable :: fine_candidate(:, :, :, :)
    real(dp), allocatable :: fine_candidate_temperature(:, :, :)
    real(dp) :: stage_theta, stage_reflux, hydro_reflux
    real(dp) :: candidate_theta, candidate_reflux
    logical :: local_ok
    character(len=32) :: selected_reconstruction, selected_limiter

    minimum_transport_theta = 1.0_dp
    maximum_reflux_correction = 0.0_dp
    ok = .false.
    candidate_theta = 1.0_dp
    candidate_reflux = 0.0_dp
    selected_reconstruction = "pcm"
    selected_limiter = "mc"
    if (present(reconstruction)) selected_reconstruction = trim(reconstruction)
    if (present(limiter)) selected_limiter = trim(limiter)
    if (.not. transport_enabled) then
      call advance_amr_reactive_strang_3d( &
        species, reactions, patch, coarse_state, coarse_temperature, &
        fine_state, fine_temperature, dx, dy, dz, dt, riemann_solver, &
        chemistry_enabled, rtol, atol, candidate_reflux, ok, &
        selected_reconstruction, selected_limiter, chemistry_integrator)
      if (ok) maximum_reflux_correction = candidate_reflux
      return
    end if
    if (.not. valid_amr_transport_shapes_3d( &
          species, transport, patch, coarse_state, coarse_temperature, &
          fine_state, fine_temperature)) return
    if (.not. ieee_is_finite(dt)) return
    if (dt <= 0.0_dp) return
    if (chemistry_enabled) then
      if (.not. all(ieee_is_finite([rtol, atol]))) return
      if (rtol <= 0.0_dp .or. atol <= 0.0_dp) return
    end if
    if ((barodiffusion_enabled .and. .not. species_diffusion_enabled) .or. &
        (trim(selected_reconstruction) /= "pcm" .and. &
          trim(selected_reconstruction) /= "characteristic_plm") .or. &
        (trim(selected_limiter) /= "minmod" .and. &
          trim(selected_limiter) /= "mc")) return

    allocate(coarse_candidate, source=coarse_state)
    allocate(coarse_candidate_temperature, source=coarse_temperature)
    allocate(fine_candidate, source=fine_state)
    allocate(fine_candidate_temperature, source=fine_temperature)
    if (chemistry_enabled) then
      call advance_amr_reactive_chemistry_3d( &
        species, reactions, patch, coarse_candidate, &
        coarse_candidate_temperature, fine_candidate, &
        fine_candidate_temperature, 0.5_dp * dt, rtol, atol, local_ok, &
        chemistry_integrator)
      if (.not. local_ok) return
    end if
    call advance_amr_reactive_transport_3d( &
      species, transport, patch, coarse_candidate, &
      coarse_candidate_temperature, fine_candidate, fine_candidate_temperature, &
      dx, dy, dz, 0.5_dp * dt, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, selected_limiter, stage_theta, stage_reflux, &
      local_ok)
    if (.not. local_ok) return
    candidate_theta = min(candidate_theta, stage_theta)
    candidate_reflux = max(candidate_reflux, stage_reflux)

    call advance_amr_reactive_hydro_3d( &
      species, patch, coarse_candidate, coarse_candidate_temperature, &
      fine_candidate, fine_candidate_temperature, dx, dy, dz, dt, &
      riemann_solver, hydro_reflux, local_ok, selected_reconstruction, &
      selected_limiter)
    if (.not. local_ok) return
    candidate_reflux = max(candidate_reflux, hydro_reflux)

    call advance_amr_reactive_transport_3d( &
      species, transport, patch, coarse_candidate, &
      coarse_candidate_temperature, fine_candidate, fine_candidate_temperature, &
      dx, dy, dz, 0.5_dp * dt, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, selected_limiter, stage_theta, stage_reflux, &
      local_ok)
    if (.not. local_ok) return
    candidate_theta = min(candidate_theta, stage_theta)
    candidate_reflux = max(candidate_reflux, stage_reflux)
    if (chemistry_enabled) then
      call advance_amr_reactive_chemistry_3d( &
        species, reactions, patch, coarse_candidate, &
        coarse_candidate_temperature, fine_candidate, &
        fine_candidate_temperature, 0.5_dp * dt, rtol, atol, local_ok, &
        chemistry_integrator)
      if (.not. local_ok) return
    end if

    coarse_state = coarse_candidate
    coarse_temperature = coarse_candidate_temperature
    fine_state = fine_candidate
    fine_temperature = fine_candidate_temperature
    minimum_transport_theta = candidate_theta
    maximum_reflux_correction = candidate_reflux
    ok = .true.
  end subroutine advance_amr_reactive_full_3d

  subroutine advance_amr_transport_euler_3d( &
      species, transport, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, dx, dy, dz, dt, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, limiter, new_coarse_state, &
      new_coarse_temperature, new_fine_state, new_fine_temperature, &
      minimum_theta, maximum_reflux_correction, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_state(:, :, :, :)
    real(dp), intent(in) :: coarse_temperature(:, :, :)
    real(dp), intent(in) :: fine_state(:, :, :, :)
    real(dp), intent(in) :: fine_temperature(:, :, :)
    real(dp), intent(in) :: dx, dy, dz, dt
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    character(len=*), intent(in) :: limiter
    real(dp), intent(out) :: new_coarse_state(:, :, :, :)
    real(dp), intent(out) :: new_coarse_temperature(:, :, :)
    real(dp), intent(out) :: new_fine_state(:, :, :, :)
    real(dp), intent(out) :: new_fine_temperature(:, :, :)
    real(dp), intent(out) :: minimum_theta, maximum_reflux_correction
    logical, intent(out) :: ok

    real(dp), allocatable :: coarse_flux_x(:, :, :, :)
    real(dp), allocatable :: coarse_flux_y(:, :, :, :)
    real(dp), allocatable :: coarse_flux_z(:, :, :, :)
    real(dp), allocatable :: fine_flux_x(:, :, :, :)
    real(dp), allocatable :: fine_flux_y(:, :, :, :)
    real(dp), allocatable :: fine_flux_z(:, :, :, :)
    real(dp), allocatable :: ghost_state(:, :, :, :)
    real(dp), allocatable :: ghost_temperature(:, :, :)
    real(dp), allocatable :: recovered_temperature(:, :, :)
    real(dp), allocatable :: fine_x_lower(:, :, :), fine_x_upper(:, :, :)
    real(dp), allocatable :: fine_y_lower(:, :, :), fine_y_upper(:, :, :)
    real(dp), allocatable :: fine_z_lower(:, :, :), fine_z_upper(:, :, :)
    real(dp), allocatable :: boundary_x_lower(:, :, :)
    real(dp), allocatable :: boundary_x_upper(:, :, :)
    real(dp), allocatable :: boundary_y_lower(:, :, :)
    real(dp), allocatable :: boundary_y_upper(:, :, :)
    real(dp), allocatable :: boundary_z_lower(:, :, :)
    real(dp), allocatable :: boundary_z_upper(:, :, :)
    real(dp) :: fine_dt, ratio_real, coarse_theta, fine_theta, alpha
    real(dp) :: interface_theta
    logical :: local_ok
    integer :: i, j, k, previous_i, previous_j, previous_k
    integer :: ratio, subcycles, substep, nvar
    integer :: covered_nx, covered_ny, covered_nz

    new_coarse_state = coarse_state
    new_coarse_temperature = coarse_temperature
    new_fine_state = fine_state
    new_fine_temperature = fine_temperature
    minimum_theta = 1.0_dp
    maximum_reflux_correction = 0.0_dp
    ok = .false.
    if (.not. valid_amr_transport_shapes_3d( &
          species, transport, patch, coarse_state, coarse_temperature, &
          fine_state, fine_temperature)) return
    if (.not. all(ieee_is_finite([dx, dy, dz, dt]))) return
    if (dt <= 0.0_dp .or. dx <= 0.0_dp .or. dy <= 0.0_dp .or. &
        dz <= 0.0_dp) return

    nvar = reactive_nvar(size(species))
    allocate(coarse_flux_x, mold=coarse_state)
    allocate(coarse_flux_y, mold=coarse_state)
    allocate(coarse_flux_z, mold=coarse_state)
    call reactive_transport_fluxes_3d( &
      species, transport, coarse_state, coarse_temperature, &
      patch%coarse_nx, patch%coarse_ny, patch%coarse_nz, dx, dy, dz, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, &
      coarse_flux_x, coarse_flux_y, coarse_flux_z, coarse_theta, local_ok)
    if (.not. local_ok) return
    do k = 1, patch%coarse_nz
      previous_k = 1 + modulo(k - 2, patch%coarse_nz)
      do j = 1, patch%coarse_ny
        previous_j = 1 + modulo(j - 2, patch%coarse_ny)
        do i = 1, patch%coarse_nx
          previous_i = 1 + modulo(i - 2, patch%coarse_nx)
          new_coarse_state(:, i, j, k) = coarse_state(:, i, j, k) - &
            dt / dx * (coarse_flux_x(:, i, j, k) - &
              coarse_flux_x(:, previous_i, j, k)) - &
            dt / dy * (coarse_flux_y(:, i, j, k) - &
              coarse_flux_y(:, i, previous_j, k)) - &
            dt / dz * (coarse_flux_z(:, i, j, k) - &
              coarse_flux_z(:, i, j, previous_k))
        end do
      end do
    end do
    allocate(recovered_temperature, mold=coarse_temperature)
    call recover_reactive_temperatures_3d( &
      species, new_coarse_state, coarse_temperature, patch%coarse_nx, &
      patch%coarse_ny, patch%coarse_nz, recovered_temperature, local_ok)
    if (.not. local_ok) return
    new_coarse_temperature = recovered_temperature
    deallocate(recovered_temperature)

    ratio = patch%refinement_ratio
    ratio_real = real(ratio, dp)
    subcycles = ratio * ratio
    fine_dt = dt / real(subcycles, dp)
    covered_nx = patch%coarse_i_upper - patch%coarse_i_lower + 1
    covered_ny = patch%coarse_j_upper - patch%coarse_j_lower + 1
    covered_nz = patch%coarse_k_upper - patch%coarse_k_lower + 1
    allocate(fine_flux_x(nvar, 0:patch%fine_nx(), &
      patch%fine_ny(), patch%fine_nz()))
    allocate(fine_flux_y(nvar, patch%fine_nx(), &
      0:patch%fine_ny(), patch%fine_nz()))
    allocate(fine_flux_z(nvar, patch%fine_nx(), &
      patch%fine_ny(), 0:patch%fine_nz()))
    allocate(ghost_state(nvar, 0:patch%fine_nx() + 1, &
      0:patch%fine_ny() + 1, 0:patch%fine_nz() + 1))
    allocate(ghost_temperature(0:patch%fine_nx() + 1, &
      0:patch%fine_ny() + 1, 0:patch%fine_nz() + 1))
    allocate(fine_x_lower(nvar, covered_ny, covered_nz))
    allocate(fine_x_upper(nvar, covered_ny, covered_nz))
    allocate(fine_y_lower(nvar, covered_nx, covered_nz))
    allocate(fine_y_upper(nvar, covered_nx, covered_nz))
    allocate(fine_z_lower(nvar, covered_nx, covered_ny))
    allocate(fine_z_upper(nvar, covered_nx, covered_ny))
    allocate(boundary_x_lower(nvar, patch%fine_ny(), patch%fine_nz()))
    allocate(boundary_x_upper(nvar, patch%fine_ny(), patch%fine_nz()))
    allocate(boundary_y_lower(nvar, patch%fine_nx(), patch%fine_nz()))
    allocate(boundary_y_upper(nvar, patch%fine_nx(), patch%fine_nz()))
    allocate(boundary_z_lower(nvar, patch%fine_nx(), patch%fine_ny()))
    allocate(boundary_z_upper(nvar, patch%fine_nx(), patch%fine_ny()))
    fine_x_lower = 0.0_dp
    fine_x_upper = 0.0_dp
    fine_y_lower = 0.0_dp
    fine_y_upper = 0.0_dp
    fine_z_lower = 0.0_dp
    fine_z_upper = 0.0_dp
    boundary_x_lower = 0.0_dp
    boundary_x_upper = 0.0_dp
    boundary_y_lower = 0.0_dp
    boundary_y_upper = 0.0_dp
    boundary_z_lower = 0.0_dp
    boundary_z_upper = 0.0_dp
    minimum_theta = coarse_theta

    do substep = 1, subcycles
      alpha = real(substep - 1, dp) / real(subcycles, dp)
      call build_fine_transport_ghosts_3d( &
        species, patch, coarse_state, coarse_temperature, &
        new_coarse_state, new_coarse_temperature, new_fine_state, &
        new_fine_temperature, alpha, limiter, ghost_state, &
        ghost_temperature, local_ok)
      if (.not. local_ok) return
      call reactive_transport_ghosted_fluxes_3d( &
        species, transport, ghost_state, ghost_temperature, &
        patch%fine_nx(), patch%fine_ny(), patch%fine_nz(), &
        dx / ratio_real, dy / ratio_real, dz / ratio_real, fine_dt, &
        viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, &
        fine_flux_x, fine_flux_y, fine_flux_z, fine_theta, local_ok)
      if (.not. local_ok) return
      minimum_theta = min(minimum_theta, fine_theta)
      do k = 1, patch%fine_nz()
        do j = 1, patch%fine_ny()
          do i = 1, patch%fine_nx()
            new_fine_state(:, i, j, k) = new_fine_state(:, i, j, k) - &
              fine_dt / (dx / ratio_real) * &
                (fine_flux_x(:, i, j, k) - &
                 fine_flux_x(:, i - 1, j, k)) - &
              fine_dt / (dy / ratio_real) * &
                (fine_flux_y(:, i, j, k) - &
                 fine_flux_y(:, i, j - 1, k)) - &
              fine_dt / (dz / ratio_real) * &
                (fine_flux_z(:, i, j, k) - &
                 fine_flux_z(:, i, j, k - 1))
          end do
        end do
      end do
      allocate(recovered_temperature, mold=new_fine_temperature)
      call recover_reactive_temperatures_3d( &
        species, new_fine_state, new_fine_temperature, patch%fine_nx(), &
        patch%fine_ny(), patch%fine_nz(), recovered_temperature, local_ok)
      if (.not. local_ok) return
      new_fine_temperature = recovered_temperature
      deallocate(recovered_temperature)
      call accumulate_fine_interface_fluxes_3d( &
        patch, fine_flux_x, fine_flux_y, fine_flux_z, &
        fine_x_lower, fine_x_upper, fine_y_lower, fine_y_upper, &
        fine_z_lower, fine_z_upper, local_ok)
      if (.not. local_ok) return
      boundary_x_lower = boundary_x_lower + fine_flux_x(:, 0, :, :)
      boundary_x_upper = boundary_x_upper + &
        fine_flux_x(:, patch%fine_nx(), :, :)
      boundary_y_lower = boundary_y_lower + fine_flux_y(:, :, 0, :)
      boundary_y_upper = boundary_y_upper + &
        fine_flux_y(:, :, patch%fine_ny(), :)
      boundary_z_lower = boundary_z_lower + fine_flux_z(:, :, :, 0)
      boundary_z_upper = boundary_z_upper + &
        fine_flux_z(:, :, :, patch%fine_nz())
    end do
    fine_x_lower = fine_x_lower / real(subcycles, dp)
    fine_x_upper = fine_x_upper / real(subcycles, dp)
    fine_y_lower = fine_y_lower / real(subcycles, dp)
    fine_y_upper = fine_y_upper / real(subcycles, dp)
    fine_z_lower = fine_z_lower / real(subcycles, dp)
    fine_z_upper = fine_z_upper / real(subcycles, dp)
    boundary_x_lower = boundary_x_lower / real(subcycles, dp)
    boundary_x_upper = boundary_x_upper / real(subcycles, dp)
    boundary_y_lower = boundary_y_lower / real(subcycles, dp)
    boundary_y_upper = boundary_y_upper / real(subcycles, dp)
    boundary_z_lower = boundary_z_lower / real(subcycles, dp)
    boundary_z_upper = boundary_z_upper / real(subcycles, dp)

    call limit_fine_transport_interfaces_3d( &
      size(species), patch, new_coarse_state, coarse_flux_x, coarse_flux_y, &
      coarse_flux_z, fine_x_lower, fine_x_upper, fine_y_lower, fine_y_upper, &
      fine_z_lower, fine_z_upper, boundary_x_lower, boundary_x_upper, &
      boundary_y_lower, boundary_y_upper, boundary_z_lower, &
      boundary_z_upper, new_fine_state, dx, dy, dz, dt, &
      interface_theta, local_ok)
    if (.not. local_ok) return
    minimum_theta = min(minimum_theta, interface_theta)
    allocate(recovered_temperature, mold=new_fine_temperature)
    call recover_reactive_temperatures_3d( &
      species, new_fine_state, new_fine_temperature, patch%fine_nx(), &
      patch%fine_ny(), patch%fine_nz(), recovered_temperature, local_ok)
    if (.not. local_ok) return
    new_fine_temperature = recovered_temperature
    deallocate(recovered_temperature)

    call reflux_coarse_3d( &
      patch, new_coarse_state, coarse_flux_x, coarse_flux_y, coarse_flux_z, &
      fine_x_lower, fine_x_upper, fine_y_lower, fine_y_upper, &
      fine_z_lower, fine_z_upper, dx, dy, dz, dt, &
      maximum_reflux_correction, local_ok)
    if (.not. local_ok) return
    call average_down_3d(new_coarse_state, new_fine_state, patch, local_ok)
    if (.not. local_ok) return
    allocate(recovered_temperature, mold=new_coarse_temperature)
    call recover_reactive_temperatures_3d( &
      species, new_coarse_state, new_coarse_temperature, patch%coarse_nx, &
      patch%coarse_ny, patch%coarse_nz, recovered_temperature, local_ok)
    if (.not. local_ok) return
    new_coarse_temperature = recovered_temperature
    ok = all(ieee_is_finite(new_coarse_state)) .and. &
      all(ieee_is_finite(new_coarse_temperature)) .and. &
      all(ieee_is_finite(new_fine_state)) .and. &
      all(ieee_is_finite(new_fine_temperature))
  end subroutine advance_amr_transport_euler_3d

  subroutine limit_fine_transport_interfaces_3d( &
      nspecies, patch, coarse_state, coarse_flux_x, coarse_flux_y, &
      coarse_flux_z, fine_x_lower, fine_x_upper, fine_y_lower, &
      fine_y_upper, fine_z_lower, fine_z_upper, boundary_x_lower, &
      boundary_x_upper, boundary_y_lower, boundary_y_upper, &
      boundary_z_lower, boundary_z_upper, fine_state, dx, dy, dz, dt, &
      minimum_theta, ok)
    integer, intent(in) :: nspecies
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_state(:, :, :, :)
    real(dp), intent(in) :: coarse_flux_x(:, :, :, :)
    real(dp), intent(in) :: coarse_flux_y(:, :, :, :)
    real(dp), intent(in) :: coarse_flux_z(:, :, :, :)
    real(dp), intent(inout) :: fine_x_lower(:, :, :)
    real(dp), intent(inout) :: fine_x_upper(:, :, :)
    real(dp), intent(inout) :: fine_y_lower(:, :, :)
    real(dp), intent(inout) :: fine_y_upper(:, :, :)
    real(dp), intent(inout) :: fine_z_lower(:, :, :)
    real(dp), intent(inout) :: fine_z_upper(:, :, :)
    real(dp), intent(in) :: boundary_x_lower(:, :, :)
    real(dp), intent(in) :: boundary_x_upper(:, :, :)
    real(dp), intent(in) :: boundary_y_lower(:, :, :)
    real(dp), intent(in) :: boundary_y_upper(:, :, :)
    real(dp), intent(in) :: boundary_z_lower(:, :, :)
    real(dp), intent(in) :: boundary_z_upper(:, :, :)
    real(dp), intent(inout) :: fine_state(:, :, :, :)
    real(dp), intent(in) :: dx, dy, dz, dt
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok

    real(dp), allocatable :: base_state(:)
    real(dp) :: theta, fine_dx, fine_dy, fine_dz
    integer :: i, j, k, fine_i, fine_j, fine_k, ratio
    logical :: local_ok

    minimum_theta = 1.0_dp
    ok = .false.
    if (nspecies < 1) return
    if (.not. patch%is_strictly_interior()) return
    if (.not. all(ieee_is_finite([dx, dy, dz, dt]))) return
    if (dx <= 0.0_dp .or. dy <= 0.0_dp .or. dz <= 0.0_dp .or. &
        dt <= 0.0_dp) return
    if (.not. all(ieee_is_finite(coarse_state))) return
    if (.not. all(ieee_is_finite(fine_state))) return

    ratio = patch%refinement_ratio
    fine_dx = dx / real(ratio, dp)
    fine_dy = dy / real(ratio, dp)
    fine_dz = dz / real(ratio, dp)
    allocate(base_state(size(coarse_state, 1)))

    do k = 1, size(fine_x_lower, 3)
      do j = 1, size(fine_x_lower, 2)
        base_state = coarse_state(:, patch%coarse_i_lower - 1, &
          patch%coarse_j_lower + j - 1, patch%coarse_k_lower + k - 1) + &
          dt / dx * coarse_flux_x(:, patch%coarse_i_lower - 1, &
            patch%coarse_j_lower + j - 1, patch%coarse_k_lower + k - 1)
        call reactive_transport_interface_theta_3d( &
          nspecies, base_state, fine_x_lower(:, j, k), -dt / dx, &
          theta, local_ok)
        if (.not. local_ok) return
        minimum_theta = min(minimum_theta, theta)
        if (theta < 1.0_dp) then
          fine_x_lower(:, j, k) = theta * fine_x_lower(:, j, k)
          do fine_k = (k - 1) * ratio + 1, k * ratio
            do fine_j = (j - 1) * ratio + 1, j * ratio
              fine_state(:, 1, fine_j, fine_k) = &
                fine_state(:, 1, fine_j, fine_k) + &
                dt / fine_dx * (theta - 1.0_dp) * &
                  boundary_x_lower(:, fine_j, fine_k)
            end do
          end do
        end if

        base_state = coarse_state(:, patch%coarse_i_upper + 1, &
          patch%coarse_j_lower + j - 1, patch%coarse_k_lower + k - 1) - &
          dt / dx * coarse_flux_x(:, patch%coarse_i_upper, &
            patch%coarse_j_lower + j - 1, patch%coarse_k_lower + k - 1)
        call reactive_transport_interface_theta_3d( &
          nspecies, base_state, fine_x_upper(:, j, k), dt / dx, &
          theta, local_ok)
        if (.not. local_ok) return
        minimum_theta = min(minimum_theta, theta)
        if (theta < 1.0_dp) then
          fine_x_upper(:, j, k) = theta * fine_x_upper(:, j, k)
          do fine_k = (k - 1) * ratio + 1, k * ratio
            do fine_j = (j - 1) * ratio + 1, j * ratio
              fine_state(:, patch%fine_nx(), fine_j, fine_k) = &
                fine_state(:, patch%fine_nx(), fine_j, fine_k) - &
                dt / fine_dx * (theta - 1.0_dp) * &
                  boundary_x_upper(:, fine_j, fine_k)
            end do
          end do
        end if
      end do
    end do

    do k = 1, size(fine_y_lower, 3)
      do i = 1, size(fine_y_lower, 2)
        base_state = coarse_state(:, patch%coarse_i_lower + i - 1, &
          patch%coarse_j_lower - 1, patch%coarse_k_lower + k - 1) + &
          dt / dy * coarse_flux_y(:, patch%coarse_i_lower + i - 1, &
            patch%coarse_j_lower - 1, patch%coarse_k_lower + k - 1)
        call reactive_transport_interface_theta_3d( &
          nspecies, base_state, fine_y_lower(:, i, k), -dt / dy, &
          theta, local_ok)
        if (.not. local_ok) return
        minimum_theta = min(minimum_theta, theta)
        if (theta < 1.0_dp) then
          fine_y_lower(:, i, k) = theta * fine_y_lower(:, i, k)
          do fine_k = (k - 1) * ratio + 1, k * ratio
            do fine_i = (i - 1) * ratio + 1, i * ratio
              fine_state(:, fine_i, 1, fine_k) = &
                fine_state(:, fine_i, 1, fine_k) + &
                dt / fine_dy * (theta - 1.0_dp) * &
                  boundary_y_lower(:, fine_i, fine_k)
            end do
          end do
        end if

        base_state = coarse_state(:, patch%coarse_i_lower + i - 1, &
          patch%coarse_j_upper + 1, patch%coarse_k_lower + k - 1) - &
          dt / dy * coarse_flux_y(:, patch%coarse_i_lower + i - 1, &
            patch%coarse_j_upper, patch%coarse_k_lower + k - 1)
        call reactive_transport_interface_theta_3d( &
          nspecies, base_state, fine_y_upper(:, i, k), dt / dy, &
          theta, local_ok)
        if (.not. local_ok) return
        minimum_theta = min(minimum_theta, theta)
        if (theta < 1.0_dp) then
          fine_y_upper(:, i, k) = theta * fine_y_upper(:, i, k)
          do fine_k = (k - 1) * ratio + 1, k * ratio
            do fine_i = (i - 1) * ratio + 1, i * ratio
              fine_state(:, fine_i, patch%fine_ny(), fine_k) = &
                fine_state(:, fine_i, patch%fine_ny(), fine_k) - &
                dt / fine_dy * (theta - 1.0_dp) * &
                  boundary_y_upper(:, fine_i, fine_k)
            end do
          end do
        end if
      end do
    end do

    do j = 1, size(fine_z_lower, 3)
      do i = 1, size(fine_z_lower, 2)
        base_state = coarse_state(:, patch%coarse_i_lower + i - 1, &
          patch%coarse_j_lower + j - 1, patch%coarse_k_lower - 1) + &
          dt / dz * coarse_flux_z(:, patch%coarse_i_lower + i - 1, &
            patch%coarse_j_lower + j - 1, patch%coarse_k_lower - 1)
        call reactive_transport_interface_theta_3d( &
          nspecies, base_state, fine_z_lower(:, i, j), -dt / dz, &
          theta, local_ok)
        if (.not. local_ok) return
        minimum_theta = min(minimum_theta, theta)
        if (theta < 1.0_dp) then
          fine_z_lower(:, i, j) = theta * fine_z_lower(:, i, j)
          do fine_j = (j - 1) * ratio + 1, j * ratio
            do fine_i = (i - 1) * ratio + 1, i * ratio
              fine_state(:, fine_i, fine_j, 1) = &
                fine_state(:, fine_i, fine_j, 1) + &
                dt / fine_dz * (theta - 1.0_dp) * &
                  boundary_z_lower(:, fine_i, fine_j)
            end do
          end do
        end if

        base_state = coarse_state(:, patch%coarse_i_lower + i - 1, &
          patch%coarse_j_lower + j - 1, patch%coarse_k_upper + 1) - &
          dt / dz * coarse_flux_z(:, patch%coarse_i_lower + i - 1, &
            patch%coarse_j_lower + j - 1, patch%coarse_k_upper)
        call reactive_transport_interface_theta_3d( &
          nspecies, base_state, fine_z_upper(:, i, j), dt / dz, &
          theta, local_ok)
        if (.not. local_ok) return
        minimum_theta = min(minimum_theta, theta)
        if (theta < 1.0_dp) then
          fine_z_upper(:, i, j) = theta * fine_z_upper(:, i, j)
          do fine_j = (j - 1) * ratio + 1, j * ratio
            do fine_i = (i - 1) * ratio + 1, i * ratio
              fine_state(:, fine_i, fine_j, patch%fine_nz()) = &
                fine_state(:, fine_i, fine_j, patch%fine_nz()) - &
                dt / fine_dz * (theta - 1.0_dp) * &
                  boundary_z_upper(:, fine_i, fine_j)
            end do
          end do
        end if
      end do
    end do

    if (.not. all(ieee_is_finite(fine_x_lower))) return
    if (.not. all(ieee_is_finite(fine_x_upper))) return
    if (.not. all(ieee_is_finite(fine_y_lower))) return
    if (.not. all(ieee_is_finite(fine_y_upper))) return
    if (.not. all(ieee_is_finite(fine_z_lower))) return
    if (.not. all(ieee_is_finite(fine_z_upper))) return
    if (.not. all(ieee_is_finite(fine_state))) return
    if (.not. ieee_is_finite(minimum_theta)) return
    ok = minimum_theta >= 0.0_dp .and. minimum_theta <= 1.0_dp
  end subroutine limit_fine_transport_interfaces_3d

  subroutine build_fine_transport_ghosts_3d( &
      species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, fine_state, fine_temperature, &
      alpha, limiter, ghost_state, ghost_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_start(:, :, :, :)
    real(dp), intent(in) :: coarse_start_temperature(:, :, :)
    real(dp), intent(in) :: coarse_end(:, :, :, :)
    real(dp), intent(in) :: coarse_end_temperature(:, :, :)
    real(dp), intent(in) :: fine_state(:, :, :, :)
    real(dp), intent(in) :: fine_temperature(:, :, :)
    real(dp), intent(in) :: alpha
    character(len=*), intent(in) :: limiter
    real(dp), intent(out) :: ghost_state(:, 0:, 0:, 0:)
    real(dp), intent(out) :: ghost_temperature(0:, 0:, 0:)
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:), interpolated_state(:)
    real(dp) :: interpolated_temperature
    logical :: local_ok
    integer :: i, j, k, coarse_i, coarse_j, coarse_k, ratio
    integer :: nx, ny, nz

    ghost_state = 0.0_dp
    ghost_temperature = 0.0_dp
    ok = .false.
    if (size(species) < 1) return
    if (.not. patch%is_strictly_interior()) return
    if (.not. ieee_is_finite(alpha)) return
    if (alpha < 0.0_dp .or. alpha > 1.0_dp) return
    nx = patch%fine_nx()
    ny = patch%fine_ny()
    nz = patch%fine_nz()
    ratio = patch%refinement_ratio
    if (size(ghost_state, 1) /= reactive_nvar(size(species)) .or. &
        size(ghost_state, 2) /= nx + 2 .or. &
        size(ghost_state, 3) /= ny + 2 .or. &
        size(ghost_state, 4) /= nz + 2 .or. &
        size(ghost_temperature, 1) /= nx + 2 .or. &
        size(ghost_temperature, 2) /= ny + 2 .or. &
        size(ghost_temperature, 3) /= nz + 2) return
    if (size(fine_state, 1) /= reactive_nvar(size(species)) .or. &
        size(fine_state, 2) /= nx .or. size(fine_state, 3) /= ny .or. &
        size(fine_state, 4) /= nz) return
    if (any(shape(fine_temperature) /= [nx, ny, nz])) return
    ghost_state(:, 1:nx, 1:ny, 1:nz) = fine_state
    ghost_temperature(1:nx, 1:ny, 1:nz) = fine_temperature
    allocate(primitive(reactive_nprim(size(species))))
    allocate(interpolated_state(reactive_nvar(size(species))))
    do k = 0, nz + 1
      coarse_k = patch%coarse_k_lower + floor_divide(k - 1, ratio)
      do j = 0, ny + 1
        coarse_j = patch%coarse_j_lower + floor_divide(j - 1, ratio)
        do i = 0, nx + 1
          if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. &
              k >= 1 .and. k <= nz) cycle
          coarse_i = patch%coarse_i_lower + floor_divide(i - 1, ratio)
          call interpolate_coarse_fine_ghost_plm_3d( &
            species, coarse_start, coarse_start_temperature, &
            coarse_end, coarse_end_temperature, coarse_i, coarse_j, &
            coarse_k, i, j, k, ratio, alpha, limiter, primitive, &
            interpolated_state, interpolated_temperature, local_ok)
          if (.not. local_ok) return
          ghost_state(:, i, j, k) = interpolated_state
          ghost_temperature(i, j, k) = interpolated_temperature
        end do
      end do
    end do
    if (.not. all(ieee_is_finite(ghost_state))) return
    if (.not. all(ieee_is_finite(ghost_temperature))) return
    ok = minval(ghost_temperature) > 0.0_dp
  end subroutine build_fine_transport_ghosts_3d

  pure integer function floor_divide(value, divisor) result(quotient)
    integer, intent(in) :: value, divisor

    if (value >= 0) then
      quotient = value / divisor
    else
      quotient = -((-value + divisor - 1) / divisor)
    end if
  end function floor_divide

  pure logical function valid_amr_transport_shapes_3d( &
      species, transport, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature) result(valid)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_state(:, :, :, :)
    real(dp), intent(in) :: coarse_temperature(:, :, :)
    real(dp), intent(in) :: fine_state(:, :, :, :)
    real(dp), intent(in) :: fine_temperature(:, :, :)

    integer :: nvar

    valid = .false.
    if (size(species) < 1) return
    if (size(transport) /= size(species)) return
    if (.not. patch%is_strictly_interior()) return
    nvar = reactive_nvar(size(species))
    if (any(shape(coarse_state) /= &
        [nvar, patch%coarse_nx, patch%coarse_ny, patch%coarse_nz])) return
    if (any(shape(coarse_temperature) /= &
        [patch%coarse_nx, patch%coarse_ny, patch%coarse_nz])) return
    if (any(shape(fine_state) /= &
        [nvar, patch%fine_nx(), patch%fine_ny(), patch%fine_nz()])) return
    if (any(shape(fine_temperature) /= &
        [patch%fine_nx(), patch%fine_ny(), patch%fine_nz()])) return
    if (.not. all(ieee_is_finite(coarse_state))) return
    if (.not. all(ieee_is_finite(coarse_temperature))) return
    if (.not. all(ieee_is_finite(fine_state))) return
    if (.not. all(ieee_is_finite(fine_temperature))) return
    if (minval(coarse_temperature) <= 0.0_dp) return
    if (minval(fine_temperature) <= 0.0_dp) return
    valid = .true.
  end function valid_amr_transport_shapes_3d

end module amr_reactive_transport_3d_mod
