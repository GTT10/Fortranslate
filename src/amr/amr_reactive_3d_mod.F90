module amr_reactive_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use constants_mod, only: density_floor, pressure_floor
  use elementary_kinetics_mod, only: &
    elementary_reaction, valid_elementary_reaction
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_mass_fraction_component, &
    reactive_conserved_to_primitive, reactive_primitive_to_conserved, &
    reactive_riemann_flux_x
  use slope_limiter_mod, only: limited_slope
  use reactive_2d_mod, only: reactive_riemann_flux_y
  use reactive_directional_flux_3d_mod, only: reactive_riemann_flux_z
  use reactive_3d_mod, only: &
    compute_reactive_cfl_timestep_3d, &
    advance_reactive_euler_ssprk2_with_fluxes_3d, &
    advance_reactive_euler_ssprk2_plm_with_fluxes_3d, &
    compute_reactive_plm_slab_face_fluxes_3d, &
    advance_reactive_chemistry_3d, &
    recover_reactive_temperatures_3d
  use amr_hierarchy_3d_mod, only: &
    amr_patch_3d, average_down_3d, composite_integrals_amr_3d
  implicit none
  private

  public :: compute_amr_reactive_cfl_timestep_3d
  public :: composite_element_integrals_amr_3d
  public :: element_integrals_from_reactive_integrals_3d
  public :: advance_amr_reactive_chemistry_3d
  public :: advance_amr_reactive_strang_3d
  public :: advance_amr_reactive_hydro_3d
  public :: accumulate_fine_interface_fluxes_3d
  public :: reflux_coarse_3d
  public :: compute_fine_patch_plm_slab_face_fluxes_3d
  public :: interpolate_coarse_fine_ghost_plm_3d

contains

  subroutine composite_element_integrals_amr_3d( &
      species, patch, coarse_state, fine_state, dx, dy, dz, elements, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_state(:, :, :, :)
    real(dp), intent(in) :: fine_state(:, :, :, :)
    real(dp), intent(in) :: dx, dy, dz
    real(dp), intent(out) :: elements(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: integrals(:)
    elements = 0.0_dp
    ok = size(elements) == 3 .and. size(species) > 0 .and. &
      size(coarse_state, 1) == reactive_nvar(size(species)) .and. &
      size(fine_state, 1) == reactive_nvar(size(species))
    if (.not. ok) return
    allocate(integrals(reactive_nvar(size(species))))
    call composite_integrals_amr_3d( &
      coarse_state, fine_state, patch, dx, dy, dz, integrals, ok)
    if (.not. ok) return
    call element_integrals_from_reactive_integrals_3d( &
      species, integrals, elements, ok)
  end subroutine composite_element_integrals_amr_3d

  pure subroutine element_integrals_from_reactive_integrals_3d( &
      species, integrals, elements, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: integrals(:)
    real(dp), intent(out) :: elements(:)
    logical, intent(out) :: ok

    real(dp) :: atom_counts(3), amount
    integer :: species_index

    elements = 0.0_dp
    ok = size(elements) == 3 .and. size(species) > 0 .and. &
      size(integrals) == reactive_nvar(size(species)) .and. &
      all(ieee_is_finite(integrals))
    if (.not. ok) return
    do species_index = 1, size(species)
      call fixed_h2o2_atom_counts( &
        species(species_index)%name, atom_counts, ok)
      if (.not. ok) then
        elements = 0.0_dp
        return
      end if
      ok = ieee_is_finite(species(species_index)%molecular_weight)
      if (ok) ok = species(species_index)%molecular_weight > 0.0_dp
      if (.not. ok) then
        elements = 0.0_dp
        return
      end if
      amount = integrals(reactive_species_component(species_index)) / &
        species(species_index)%molecular_weight
      elements = elements + amount * atom_counts
    end do
    ok = all(ieee_is_finite(elements))
    if (.not. ok) elements = 0.0_dp
  end subroutine element_integrals_from_reactive_integrals_3d

  pure subroutine fixed_h2o2_atom_counts(name, atom_counts, ok)
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: atom_counts(3)
    logical, intent(out) :: ok

    atom_counts = 0.0_dp
    ok = .true.
    select case (trim(name))
    case ("H2")
      atom_counts = [2.0_dp, 0.0_dp, 0.0_dp]
    case ("H")
      atom_counts = [1.0_dp, 0.0_dp, 0.0_dp]
    case ("O")
      atom_counts = [0.0_dp, 1.0_dp, 0.0_dp]
    case ("O2")
      atom_counts = [0.0_dp, 2.0_dp, 0.0_dp]
    case ("OH")
      atom_counts = [1.0_dp, 1.0_dp, 0.0_dp]
    case ("H2O")
      atom_counts = [2.0_dp, 1.0_dp, 0.0_dp]
    case ("HO2")
      atom_counts = [1.0_dp, 2.0_dp, 0.0_dp]
    case ("H2O2")
      atom_counts = [2.0_dp, 2.0_dp, 0.0_dp]
    case ("N2")
      atom_counts = [0.0_dp, 0.0_dp, 2.0_dp]
    case ("AR")
      continue
    case default
      ok = .false.
    end select
  end subroutine fixed_h2o2_atom_counts

  subroutine compute_amr_reactive_cfl_timestep_3d( &
      species, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, dx, dy, dz, cfl, dt, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_state(:, :, :, :)
    real(dp), intent(in) :: coarse_temperature(:, :, :)
    real(dp), intent(in) :: fine_state(:, :, :, :)
    real(dp), intent(in) :: fine_temperature(:, :, :)
    real(dp), intent(in) :: dx, dy, dz, cfl
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok

    real(dp) :: coarse_dt, fine_dt
    logical :: local_ok
    integer :: ratio

    dt = 0.0_dp
    ok = patch%is_strictly_interior() .and. valid_amr_reactive_shapes( &
      species, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature)
    if (.not. ok) return
    ratio = patch%refinement_ratio
    call compute_reactive_cfl_timestep_3d( &
      species, coarse_state, coarse_temperature, patch%coarse_nx, &
      patch%coarse_ny, patch%coarse_nz, dx, dy, dz, cfl, &
      coarse_dt, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call compute_reactive_cfl_timestep_3d( &
      species, fine_state, fine_temperature, patch%fine_nx(), &
      patch%fine_ny(), patch%fine_nz(), dx / real(ratio, dp), &
      dy / real(ratio, dp), dz / real(ratio, dp), cfl, fine_dt, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    dt = min(coarse_dt, real(ratio, dp) * fine_dt)
    ok = ieee_is_finite(dt)
    if (ok) ok = dt > 0.0_dp
  end subroutine compute_amr_reactive_cfl_timestep_3d

  subroutine advance_amr_reactive_chemistry_3d( &
      species, reactions, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, interval, rtol, atol, ok, &
      chemistry_integrator)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(inout) :: coarse_state(:, :, :, :)
    real(dp), intent(inout) :: coarse_temperature(:, :, :)
    real(dp), intent(inout) :: fine_state(:, :, :, :)
    real(dp), intent(inout) :: fine_temperature(:, :, :)
    real(dp), intent(in) :: interval, rtol, atol
    logical, intent(out) :: ok
    character(len=*), intent(in), optional :: chemistry_integrator

    real(dp), allocatable :: coarse_candidate(:, :, :, :)
    real(dp), allocatable :: coarse_candidate_temperature(:, :, :)
    real(dp), allocatable :: fine_candidate(:, :, :, :)
    real(dp), allocatable :: fine_candidate_temperature(:, :, :)
    real(dp), allocatable :: synchronized_temperature(:, :, :)
    logical :: local_ok

    ok = .false.
    if (.not. valid_amr_reactive_chemistry_inputs( &
          species, reactions, patch, coarse_state, coarse_temperature, &
          fine_state, fine_temperature, interval, rtol, atol)) return
    if (present(chemistry_integrator)) then
      if (.not. valid_amr_reactive_chemistry_integrator( &
            chemistry_integrator)) return
    end if

    ! Keep both levels private until chemistry and the synchronization have
    ! succeeded.  This makes the source phase transactional as a hierarchy.
    allocate(coarse_candidate, source=coarse_state)
    allocate(coarse_candidate_temperature, source=coarse_temperature)
    allocate(fine_candidate, source=fine_state)
    allocate(fine_candidate_temperature, source=fine_temperature)
    if (present(chemistry_integrator)) then
      call advance_reactive_chemistry_3d( &
        species, reactions, coarse_candidate, coarse_candidate_temperature, &
        patch%coarse_nx, patch%coarse_ny, patch%coarse_nz, interval, rtol, &
        atol, local_ok, chemistry_integrator=chemistry_integrator)
    else
      call advance_reactive_chemistry_3d( &
        species, reactions, coarse_candidate, coarse_candidate_temperature, &
        patch%coarse_nx, patch%coarse_ny, patch%coarse_nz, interval, rtol, &
        atol, local_ok)
    end if
    if (.not. local_ok) return

    if (present(chemistry_integrator)) then
      call advance_reactive_chemistry_3d( &
        species, reactions, fine_candidate, fine_candidate_temperature, &
        patch%fine_nx(), patch%fine_ny(), patch%fine_nz(), interval, rtol, &
        atol, local_ok, chemistry_integrator=chemistry_integrator)
    else
      call advance_reactive_chemistry_3d( &
        species, reactions, fine_candidate, fine_candidate_temperature, &
        patch%fine_nx(), patch%fine_ny(), patch%fine_nz(), interval, rtol, &
        atol, local_ok)
    end if
    if (.not. local_ok) return

    call average_down_3d( &
      coarse_candidate, fine_candidate, patch, local_ok)
    if (.not. local_ok) return
    allocate(synchronized_temperature, mold=coarse_temperature)
    call recover_reactive_temperatures_3d( &
      species, coarse_candidate, coarse_candidate_temperature, &
      patch%coarse_nx, patch%coarse_ny, patch%coarse_nz, &
      synchronized_temperature, local_ok)
    if (.not. local_ok) return
    if (.not. all(ieee_is_finite(coarse_candidate)) .or. &
        .not. all(ieee_is_finite(synchronized_temperature)) .or. &
        .not. all(ieee_is_finite(fine_candidate)) .or. &
        .not. all(ieee_is_finite(fine_candidate_temperature))) return

    coarse_state = coarse_candidate
    coarse_temperature = synchronized_temperature
    fine_state = fine_candidate
    fine_temperature = fine_candidate_temperature
    ok = .true.
  end subroutine advance_amr_reactive_chemistry_3d

  subroutine advance_amr_reactive_strang_3d( &
      species, reactions, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, dx, dy, dz, dt, riemann_solver, &
      chemistry_enabled, rtol, atol, maximum_reflux_correction, ok, &
      reconstruction, limiter, chemistry_integrator)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(inout) :: coarse_state(:, :, :, :)
    real(dp), intent(inout) :: coarse_temperature(:, :, :)
    real(dp), intent(inout) :: fine_state(:, :, :, :)
    real(dp), intent(inout) :: fine_temperature(:, :, :)
    real(dp), intent(in) :: dx, dy, dz, dt
    character(len=*), intent(in) :: riemann_solver
    logical, intent(in) :: chemistry_enabled
    real(dp), intent(in) :: rtol, atol
    real(dp), intent(out) :: maximum_reflux_correction
    logical, intent(out) :: ok
    character(len=*), intent(in), optional :: reconstruction, limiter
    character(len=*), intent(in), optional :: chemistry_integrator

    real(dp), allocatable :: coarse_candidate(:, :, :, :)
    real(dp), allocatable :: coarse_candidate_temperature(:, :, :)
    real(dp), allocatable :: fine_candidate(:, :, :, :)
    real(dp), allocatable :: fine_candidate_temperature(:, :, :)
    real(dp) :: candidate_reflux
    logical :: local_ok

    maximum_reflux_correction = 0.0_dp
    ok = .false.
    if (.not. patch%is_strictly_interior()) return
    if (.not. valid_amr_reactive_shapes( &
          species, patch, coarse_state, coarse_temperature, &
          fine_state, fine_temperature)) return
    if (.not. all(ieee_is_finite([dx, dy, dz, dt]))) return
    if (dx <= 0.0_dp .or. dy <= 0.0_dp .or. dz <= 0.0_dp .or. &
        dt <= 0.0_dp) return

    ! Dispatch the old operation directly when chemistry is disabled.  In
    ! addition to avoiding unnecessary copies, this preserves its exact bit
    ! path and its existing optional reconstruction/limiter behavior.
    if (.not. chemistry_enabled) then
      call dispatch_amr_reactive_hydro_3d( &
        species, patch, coarse_state, coarse_temperature, fine_state, &
        fine_temperature, dx, dy, dz, dt, riemann_solver, &
        maximum_reflux_correction, ok, reconstruction, limiter)
      return
    end if

    if (.not. valid_amr_reactive_chemistry_inputs( &
          species, reactions, patch, coarse_state, coarse_temperature, &
          fine_state, fine_temperature, 0.5_dp * dt, rtol, atol)) return
    if (present(chemistry_integrator)) then
      if (.not. valid_amr_reactive_chemistry_integrator( &
            chemistry_integrator)) return
    end if

    ! The full split is also transactional: neither public hierarchy level is
    ! assigned until both source phases and the AMR hydro phase succeed.
    allocate(coarse_candidate, source=coarse_state)
    allocate(coarse_candidate_temperature, source=coarse_temperature)
    allocate(fine_candidate, source=fine_state)
    allocate(fine_candidate_temperature, source=fine_temperature)
    if (present(chemistry_integrator)) then
      call advance_amr_reactive_chemistry_3d( &
        species, reactions, patch, coarse_candidate, &
        coarse_candidate_temperature, fine_candidate, &
        fine_candidate_temperature, 0.5_dp * dt, rtol, atol, local_ok, &
        chemistry_integrator=chemistry_integrator)
    else
      call advance_amr_reactive_chemistry_3d( &
        species, reactions, patch, coarse_candidate, &
        coarse_candidate_temperature, fine_candidate, &
        fine_candidate_temperature, 0.5_dp * dt, rtol, atol, local_ok)
    end if
    if (.not. local_ok) return

    candidate_reflux = 0.0_dp
    call dispatch_amr_reactive_hydro_3d( &
      species, patch, coarse_candidate, coarse_candidate_temperature, &
      fine_candidate, fine_candidate_temperature, dx, dy, dz, dt, &
      riemann_solver, candidate_reflux, local_ok, reconstruction, limiter)
    if (.not. local_ok) return

    if (present(chemistry_integrator)) then
      call advance_amr_reactive_chemistry_3d( &
        species, reactions, patch, coarse_candidate, &
        coarse_candidate_temperature, fine_candidate, &
        fine_candidate_temperature, 0.5_dp * dt, rtol, atol, local_ok, &
        chemistry_integrator=chemistry_integrator)
    else
      call advance_amr_reactive_chemistry_3d( &
        species, reactions, patch, coarse_candidate, &
        coarse_candidate_temperature, fine_candidate, &
        fine_candidate_temperature, 0.5_dp * dt, rtol, atol, local_ok)
    end if
    if (.not. local_ok) return
    if (.not. all(ieee_is_finite(coarse_candidate)) .or. &
        .not. all(ieee_is_finite(coarse_candidate_temperature)) .or. &
        .not. all(ieee_is_finite(fine_candidate)) .or. &
        .not. all(ieee_is_finite(fine_candidate_temperature))) return

    coarse_state = coarse_candidate
    coarse_temperature = coarse_candidate_temperature
    fine_state = fine_candidate
    fine_temperature = fine_candidate_temperature
    maximum_reflux_correction = candidate_reflux
    ok = .true.
  end subroutine advance_amr_reactive_strang_3d

  subroutine advance_amr_reactive_hydro_3d( &
      species, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, dx, dy, dz, dt, riemann_solver, &
      maximum_reflux_correction, ok, reconstruction, limiter)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(inout) :: coarse_state(:, :, :, :)
    real(dp), intent(inout) :: coarse_temperature(:, :, :)
    real(dp), intent(inout) :: fine_state(:, :, :, :)
    real(dp), intent(inout) :: fine_temperature(:, :, :)
    real(dp), intent(in) :: dx, dy, dz, dt
    character(len=*), intent(in) :: riemann_solver
    real(dp), intent(out) :: maximum_reflux_correction
    logical, intent(out) :: ok
    character(len=*), intent(in), optional :: reconstruction, limiter

    real(dp), allocatable :: coarse_start(:, :, :, :)
    real(dp), allocatable :: coarse_candidate(:, :, :, :)
    real(dp), allocatable :: fine_candidate(:, :, :, :)
    real(dp), allocatable :: coarse_start_temperature(:, :, :)
    real(dp), allocatable :: coarse_candidate_temperature(:, :, :)
    real(dp), allocatable :: synchronized_coarse_temperature(:, :, :)
    real(dp), allocatable :: fine_candidate_temperature(:, :, :)
    real(dp), allocatable :: coarse_flux_x(:, :, :, :)
    real(dp), allocatable :: coarse_flux_y(:, :, :, :)
    real(dp), allocatable :: coarse_flux_z(:, :, :, :)
    real(dp), allocatable :: fine_flux_x(:, :, :, :)
    real(dp), allocatable :: fine_flux_y(:, :, :, :)
    real(dp), allocatable :: fine_flux_z(:, :, :, :)
    real(dp), allocatable :: fine_x_lower(:, :, :), fine_x_upper(:, :, :)
    real(dp), allocatable :: fine_y_lower(:, :, :), fine_y_upper(:, :, :)
    real(dp), allocatable :: fine_z_lower(:, :, :), fine_z_upper(:, :, :)
    real(dp) :: fine_dt, alpha_start, alpha_end
    logical :: local_ok
    integer :: nvar, ratio, substep
    integer :: covered_nx, covered_ny, covered_nz
    character(len=32) :: selected_reconstruction, selected_limiter

    maximum_reflux_correction = 0.0_dp
    ok = .false.
    selected_reconstruction = "pcm"
    selected_limiter = "mc"
    if (present(reconstruction)) then
      selected_reconstruction = trim(reconstruction)
    end if
    if (present(limiter)) selected_limiter = trim(limiter)
    if (.not. valid_amr_reactive_shapes( &
          species, patch, coarse_state, coarse_temperature, &
          fine_state, fine_temperature)) return
    if (.not. patch%is_strictly_interior()) return
    if (.not. all(ieee_is_finite([dx, dy, dz, dt]))) return
    if (dx <= 0.0_dp .or. dy <= 0.0_dp .or. dz <= 0.0_dp .or. &
        dt <= 0.0_dp) return
    if (trim(selected_reconstruction) /= "pcm" .and. &
        trim(selected_reconstruction) /= "characteristic_plm") return
    if (trim(selected_limiter) /= "minmod" .and. &
        trim(selected_limiter) /= "mc") return

    nvar = reactive_nvar(size(species))
    ratio = patch%refinement_ratio
    covered_nx = patch%coarse_i_upper - patch%coarse_i_lower + 1
    covered_ny = patch%coarse_j_upper - patch%coarse_j_lower + 1
    covered_nz = patch%coarse_k_upper - patch%coarse_k_lower + 1
    coarse_start = coarse_state
    coarse_candidate = coarse_state
    fine_candidate = fine_state
    allocate(coarse_start_temperature, mold=coarse_temperature)
    coarse_candidate_temperature = coarse_temperature
    fine_candidate_temperature = fine_temperature
    call recover_reactive_temperatures_3d( &
      species, coarse_start, coarse_temperature, patch%coarse_nx, &
      patch%coarse_ny, patch%coarse_nz, coarse_start_temperature, local_ok)
    if (.not. local_ok) return

    allocate(coarse_flux_x, mold=coarse_state)
    allocate(coarse_flux_y, mold=coarse_state)
    allocate(coarse_flux_z, mold=coarse_state)
    select case (trim(selected_reconstruction))
    case ("pcm")
      call advance_reactive_euler_ssprk2_with_fluxes_3d( &
        species, coarse_candidate, coarse_candidate_temperature, &
        patch%coarse_nx, patch%coarse_ny, patch%coarse_nz, &
        dx, dy, dz, dt, riemann_solver, coarse_flux_x, coarse_flux_y, &
        coarse_flux_z, local_ok)
    case ("characteristic_plm")
      call advance_reactive_euler_ssprk2_plm_with_fluxes_3d( &
        species, coarse_candidate, coarse_candidate_temperature, &
        patch%coarse_nx, patch%coarse_ny, patch%coarse_nz, &
        dx, dy, dz, dt, selected_limiter, riemann_solver, &
        coarse_flux_x, coarse_flux_y, coarse_flux_z, local_ok)
    end select
    if (.not. local_ok) return

    allocate(fine_flux_x(nvar, 0:patch%fine_nx(), &
      patch%fine_ny(), patch%fine_nz()))
    allocate(fine_flux_y(nvar, patch%fine_nx(), &
      0:patch%fine_ny(), patch%fine_nz()))
    allocate(fine_flux_z(nvar, patch%fine_nx(), &
      patch%fine_ny(), 0:patch%fine_nz()))
    allocate(fine_x_lower(nvar, covered_ny, covered_nz), &
      fine_x_upper(nvar, covered_ny, covered_nz))
    allocate(fine_y_lower(nvar, covered_nx, covered_nz), &
      fine_y_upper(nvar, covered_nx, covered_nz))
    allocate(fine_z_lower(nvar, covered_nx, covered_ny), &
      fine_z_upper(nvar, covered_nx, covered_ny))
    fine_x_lower = 0.0_dp
    fine_x_upper = 0.0_dp
    fine_y_lower = 0.0_dp
    fine_y_upper = 0.0_dp
    fine_z_lower = 0.0_dp
    fine_z_upper = 0.0_dp

    fine_dt = dt / real(ratio, dp)
    do substep = 1, ratio
      alpha_start = real(substep - 1, dp) / real(ratio, dp)
      alpha_end = real(substep, dp) / real(ratio, dp)
      call advance_fine_patch_ssprk2_3d( &
        species, patch, coarse_start, coarse_start_temperature, &
        coarse_candidate, coarse_candidate_temperature, &
        fine_candidate, fine_candidate_temperature, &
        dx / real(ratio, dp), dy / real(ratio, dp), &
        dz / real(ratio, dp), fine_dt, alpha_start, alpha_end, &
        riemann_solver, selected_reconstruction, selected_limiter, &
        fine_flux_x, fine_flux_y, fine_flux_z, local_ok)
      if (.not. local_ok) return
      call accumulate_fine_interface_fluxes_3d( &
        patch, fine_flux_x, fine_flux_y, fine_flux_z, &
        fine_x_lower, fine_x_upper, fine_y_lower, fine_y_upper, &
        fine_z_lower, fine_z_upper, local_ok)
      if (.not. local_ok) return
    end do
    fine_x_lower = fine_x_lower / real(ratio, dp)
    fine_x_upper = fine_x_upper / real(ratio, dp)
    fine_y_lower = fine_y_lower / real(ratio, dp)
    fine_y_upper = fine_y_upper / real(ratio, dp)
    fine_z_lower = fine_z_lower / real(ratio, dp)
    fine_z_upper = fine_z_upper / real(ratio, dp)

    call reflux_coarse_3d( &
      patch, coarse_candidate, coarse_flux_x, coarse_flux_y, coarse_flux_z, &
      fine_x_lower, fine_x_upper, fine_y_lower, fine_y_upper, &
      fine_z_lower, fine_z_upper, dx, dy, dz, dt, &
      maximum_reflux_correction, local_ok)
    if (.not. local_ok) return
    call average_down_3d(coarse_candidate, fine_candidate, patch, local_ok)
    if (.not. local_ok) return
    allocate(synchronized_coarse_temperature, mold=coarse_temperature)
    call recover_reactive_temperatures_3d( &
      species, coarse_candidate, coarse_candidate_temperature, &
      patch%coarse_nx, patch%coarse_ny, patch%coarse_nz, &
      synchronized_coarse_temperature, local_ok)
    if (.not. local_ok) return
    coarse_candidate_temperature = synchronized_coarse_temperature

    coarse_state = coarse_candidate
    coarse_temperature = coarse_candidate_temperature
    fine_state = fine_candidate
    fine_temperature = fine_candidate_temperature
    ok = .true.
  end subroutine advance_amr_reactive_hydro_3d

  subroutine advance_fine_patch_ssprk2_3d( &
      species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, state, temperature, &
      dx, dy, dz, dt, alpha_start, alpha_end, riemann_solver, &
      reconstruction, limiter, face_flux_x, face_flux_y, face_flux_z, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_start(:, :, :, :)
    real(dp), intent(in) :: coarse_start_temperature(:, :, :)
    real(dp), intent(in) :: coarse_end(:, :, :, :)
    real(dp), intent(in) :: coarse_end_temperature(:, :, :)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(in) :: dx, dy, dz, dt, alpha_start, alpha_end
    character(len=*), intent(in) :: riemann_solver
    character(len=*), intent(in) :: reconstruction, limiter
    real(dp), intent(out) :: face_flux_x(:, 0:, :, :)
    real(dp), intent(out) :: face_flux_y(:, :, 0:, :)
    real(dp), intent(out) :: face_flux_z(:, :, :, 0:)
    logical, intent(out) :: ok

    real(dp), allocatable :: old_state(:, :, :, :)
    real(dp), allocatable :: stage_state(:, :, :, :)
    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: rhs(:, :, :, :)
    real(dp), allocatable :: first_flux_x(:, :, :, :)
    real(dp), allocatable :: first_flux_y(:, :, :, :)
    real(dp), allocatable :: first_flux_z(:, :, :, :)
    real(dp), allocatable :: second_flux_x(:, :, :, :)
    real(dp), allocatable :: second_flux_y(:, :, :, :)
    real(dp), allocatable :: second_flux_z(:, :, :, :)
    real(dp), allocatable :: old_temperature(:, :, :)
    real(dp), allocatable :: stage_temperature(:, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    logical :: local_ok
    integer :: nvar, nx, ny, nz

    face_flux_x = 0.0_dp
    face_flux_y = 0.0_dp
    face_flux_z = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    nx = patch%fine_nx()
    ny = patch%fine_ny()
    nz = patch%fine_nz()
    if (.not. valid_fine_advance_shapes( &
          species, patch, coarse_start, coarse_start_temperature, &
          coarse_end, coarse_end_temperature, state, temperature, &
          face_flux_x, face_flux_y, face_flux_z)) return
    if (.not. all(ieee_is_finite( &
          [dx, dy, dz, dt, alpha_start, alpha_end]))) return
    if (dx <= 0.0_dp .or. dy <= 0.0_dp .or. dz <= 0.0_dp .or. &
        dt <= 0.0_dp .or. alpha_start < 0.0_dp .or. alpha_end > 1.0_dp .or. &
        alpha_end <= alpha_start) return

    old_state = state
    allocate(old_temperature(nx, ny, nz))
    call recover_reactive_temperatures_3d( &
      species, old_state, temperature, nx, ny, nz, old_temperature, local_ok)
    if (.not. local_ok) return
    allocate(stage_state(nvar, nx, ny, nz))
    allocate(candidate_state(nvar, nx, ny, nz))
    allocate(rhs(nvar, nx, ny, nz))
    allocate(stage_temperature(nx, ny, nz))
    allocate(candidate_temperature(nx, ny, nz))
    allocate(first_flux_x(nvar, 0:nx, ny, nz))
    allocate(first_flux_y(nvar, nx, 0:ny, nz))
    allocate(first_flux_z(nvar, nx, ny, 0:nz))
    allocate(second_flux_x(nvar, 0:nx, ny, nz))
    allocate(second_flux_y(nvar, nx, 0:ny, nz))
    allocate(second_flux_z(nvar, nx, ny, 0:nz))

    call compute_fine_patch_face_fluxes_3d( &
      species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, old_state, old_temperature, &
      alpha_start, riemann_solver, reconstruction, limiter, &
      first_flux_x, first_flux_y, &
      first_flux_z, local_ok)
    if (.not. local_ok) return
    call fine_flux_divergence_3d( &
      first_flux_x, first_flux_y, first_flux_z, dx, dy, dz, rhs, local_ok)
    if (.not. local_ok) return
    stage_state = old_state + dt * rhs
    call recover_reactive_temperatures_3d( &
      species, stage_state, old_temperature, nx, ny, nz, &
      stage_temperature, local_ok)
    if (.not. local_ok) return

    call compute_fine_patch_face_fluxes_3d( &
      species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, stage_state, stage_temperature, &
      alpha_end, riemann_solver, reconstruction, limiter, &
      second_flux_x, second_flux_y, &
      second_flux_z, local_ok)
    if (.not. local_ok) return
    call fine_flux_divergence_3d( &
      second_flux_x, second_flux_y, second_flux_z, dx, dy, dz, rhs, local_ok)
    if (.not. local_ok) return
    candidate_state = 0.5_dp * old_state + &
      0.5_dp * (stage_state + dt * rhs)
    call recover_reactive_temperatures_3d( &
      species, candidate_state, stage_temperature, nx, ny, nz, &
      candidate_temperature, local_ok)
    if (.not. local_ok) return

    face_flux_x = 0.5_dp * (first_flux_x + second_flux_x)
    face_flux_y = 0.5_dp * (first_flux_y + second_flux_y)
    face_flux_z = 0.5_dp * (first_flux_z + second_flux_z)
    state = candidate_state
    temperature = candidate_temperature
    ok = .true.
  end subroutine advance_fine_patch_ssprk2_3d

  subroutine compute_fine_patch_face_fluxes_3d( &
      species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, state, temperature, alpha, &
      riemann_solver, reconstruction, limiter, &
      face_flux_x, face_flux_y, face_flux_z, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_start(:, :, :, :)
    real(dp), intent(in) :: coarse_start_temperature(:, :, :)
    real(dp), intent(in) :: coarse_end(:, :, :, :)
    real(dp), intent(in) :: coarse_end_temperature(:, :, :)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(in) :: alpha
    character(len=*), intent(in) :: riemann_solver
    character(len=*), intent(in) :: reconstruction, limiter
    real(dp), intent(out) :: face_flux_x(:, 0:, :, :)
    real(dp), intent(out) :: face_flux_y(:, :, 0:, :)
    real(dp), intent(out) :: face_flux_z(:, :, :, 0:)
    logical, intent(out) :: ok

    real(dp), allocatable :: ghost_state(:), primitive(:)
    real(dp) :: ghost_temperature
    logical :: face_ok
    integer :: i, j, k, coarse_i, coarse_j, coarse_k
    integer :: nx, ny, nz, ratio

    face_flux_x = 0.0_dp
    face_flux_y = 0.0_dp
    face_flux_z = 0.0_dp
    ok = .false.
    nx = patch%fine_nx()
    ny = patch%fine_ny()
    nz = patch%fine_nz()
    ratio = patch%refinement_ratio
    if (.not. valid_fine_advance_shapes( &
          species, patch, coarse_start, coarse_start_temperature, &
          coarse_end, coarse_end_temperature, state, temperature, &
          face_flux_x, face_flux_y, face_flux_z)) return
    if (.not. ieee_is_finite(alpha)) return
    if (alpha < 0.0_dp .or. alpha > 1.0_dp) return
    if (trim(reconstruction) /= "pcm" .and. &
        trim(reconstruction) /= "characteristic_plm") return
    if (trim(limiter) /= "minmod" .and. trim(limiter) /= "mc") return
    if (trim(reconstruction) == "characteristic_plm") then
      call compute_fine_patch_plm_face_fluxes_3d( &
        species, patch, coarse_start, coarse_start_temperature, &
        coarse_end, coarse_end_temperature, state, temperature, alpha, &
        limiter, riemann_solver, face_flux_x, face_flux_y, face_flux_z, ok)
      return
    end if
    allocate(ghost_state(reactive_nvar(size(species))))
    allocate(primitive(reactive_nprim(size(species))))

    do k = 1, nz
      coarse_k = patch%coarse_k_lower + (k - 1) / ratio
      do j = 1, ny
        coarse_j = patch%coarse_j_lower + (j - 1) / ratio
        call interpolate_coarse_cell( &
          species, coarse_start, coarse_start_temperature, &
          coarse_end, coarse_end_temperature, &
          patch%coarse_i_lower - 1, coarse_j, coarse_k, alpha, &
          primitive, ghost_state, ghost_temperature, face_ok)
        if (.not. face_ok) return
        call reactive_riemann_flux_x( &
          species, ghost_state, state(:, 1, j, k), ghost_temperature, &
          temperature(1, j, k), riemann_solver, &
          face_flux_x(:, 0, j, k), face_ok)
        if (.not. face_ok) return
        call interpolate_coarse_cell( &
          species, coarse_start, coarse_start_temperature, &
          coarse_end, coarse_end_temperature, &
          patch%coarse_i_upper + 1, coarse_j, coarse_k, alpha, &
          primitive, ghost_state, ghost_temperature, face_ok)
        if (.not. face_ok) return
        call reactive_riemann_flux_x( &
          species, state(:, nx, j, k), ghost_state, &
          temperature(nx, j, k), ghost_temperature, riemann_solver, &
          face_flux_x(:, nx, j, k), face_ok)
        if (.not. face_ok) return
        do i = 1, nx - 1
          call reactive_riemann_flux_x( &
            species, state(:, i, j, k), state(:, i + 1, j, k), &
            temperature(i, j, k), temperature(i + 1, j, k), &
            riemann_solver, face_flux_x(:, i, j, k), face_ok)
          if (.not. face_ok) return
        end do
      end do
    end do

    do k = 1, nz
      coarse_k = patch%coarse_k_lower + (k - 1) / ratio
      do i = 1, nx
        coarse_i = patch%coarse_i_lower + (i - 1) / ratio
        call interpolate_coarse_cell( &
          species, coarse_start, coarse_start_temperature, &
          coarse_end, coarse_end_temperature, coarse_i, &
          patch%coarse_j_lower - 1, coarse_k, alpha, &
          primitive, ghost_state, ghost_temperature, face_ok)
        if (.not. face_ok) return
        call reactive_riemann_flux_y( &
          species, ghost_state, state(:, i, 1, k), ghost_temperature, &
          temperature(i, 1, k), riemann_solver, &
          face_flux_y(:, i, 0, k), face_ok)
        if (.not. face_ok) return
        call interpolate_coarse_cell( &
          species, coarse_start, coarse_start_temperature, &
          coarse_end, coarse_end_temperature, coarse_i, &
          patch%coarse_j_upper + 1, coarse_k, alpha, &
          primitive, ghost_state, ghost_temperature, face_ok)
        if (.not. face_ok) return
        call reactive_riemann_flux_y( &
          species, state(:, i, ny, k), ghost_state, &
          temperature(i, ny, k), ghost_temperature, riemann_solver, &
          face_flux_y(:, i, ny, k), face_ok)
        if (.not. face_ok) return
        do j = 1, ny - 1
          call reactive_riemann_flux_y( &
            species, state(:, i, j, k), state(:, i, j + 1, k), &
            temperature(i, j, k), temperature(i, j + 1, k), &
            riemann_solver, face_flux_y(:, i, j, k), face_ok)
          if (.not. face_ok) return
        end do
      end do
    end do

    do j = 1, ny
      coarse_j = patch%coarse_j_lower + (j - 1) / ratio
      do i = 1, nx
        coarse_i = patch%coarse_i_lower + (i - 1) / ratio
        call interpolate_coarse_cell( &
          species, coarse_start, coarse_start_temperature, &
          coarse_end, coarse_end_temperature, coarse_i, coarse_j, &
          patch%coarse_k_lower - 1, alpha, primitive, ghost_state, &
          ghost_temperature, face_ok)
        if (.not. face_ok) return
        call reactive_riemann_flux_z( &
          species, ghost_state, state(:, i, j, 1), ghost_temperature, &
          temperature(i, j, 1), riemann_solver, &
          face_flux_z(:, i, j, 0), face_ok)
        if (.not. face_ok) return
        call interpolate_coarse_cell( &
          species, coarse_start, coarse_start_temperature, &
          coarse_end, coarse_end_temperature, coarse_i, coarse_j, &
          patch%coarse_k_upper + 1, alpha, primitive, ghost_state, &
          ghost_temperature, face_ok)
        if (.not. face_ok) return
        call reactive_riemann_flux_z( &
          species, state(:, i, j, nz), ghost_state, &
          temperature(i, j, nz), ghost_temperature, riemann_solver, &
          face_flux_z(:, i, j, nz), face_ok)
        if (.not. face_ok) return
        do k = 1, nz - 1
          call reactive_riemann_flux_z( &
            species, state(:, i, j, k), state(:, i, j, k + 1), &
            temperature(i, j, k), temperature(i, j, k + 1), &
            riemann_solver, face_flux_z(:, i, j, k), face_ok)
          if (.not. face_ok) return
        end do
      end do
    end do
    ok = all(ieee_is_finite(face_flux_x)) .and. &
      all(ieee_is_finite(face_flux_y)) .and. &
      all(ieee_is_finite(face_flux_z))
  end subroutine compute_fine_patch_face_fluxes_3d

  subroutine compute_fine_patch_plm_face_fluxes_3d( &
      species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, state, temperature, alpha, &
      limiter, riemann_solver, face_flux_x, face_flux_y, face_flux_z, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_start(:, :, :, :)
    real(dp), intent(in) :: coarse_start_temperature(:, :, :)
    real(dp), intent(in) :: coarse_end(:, :, :, :)
    real(dp), intent(in) :: coarse_end_temperature(:, :, :)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(in) :: alpha
    character(len=*), intent(in) :: limiter, riemann_solver
    real(dp), intent(out) :: face_flux_x(:, 0:, :, :)
    real(dp), intent(out) :: face_flux_y(:, :, 0:, :)
    real(dp), intent(out) :: face_flux_z(:, :, :, 0:)
    logical, intent(out) :: ok

    call compute_fine_patch_plm_slab_face_fluxes_3d( &
      species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, state, temperature, alpha, &
      1, patch%fine_nx(), limiter, riemann_solver, face_flux_x, &
      face_flux_y, face_flux_z, ok)
  end subroutine compute_fine_patch_plm_face_fluxes_3d

  subroutine compute_fine_patch_plm_slab_face_fluxes_3d( &
      species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, state, temperature, alpha, &
      first_i, last_i, limiter, riemann_solver, face_flux_x, face_flux_y, &
      face_flux_z, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_start(:, :, :, :)
    real(dp), intent(in) :: coarse_start_temperature(:, :, :)
    real(dp), intent(in) :: coarse_end(:, :, :, :)
    real(dp), intent(in) :: coarse_end_temperature(:, :, :)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(in) :: alpha
    integer, intent(in) :: first_i, last_i
    character(len=*), intent(in) :: limiter, riemann_solver
    real(dp), intent(out) :: face_flux_x(:, 0:, :, :)
    real(dp), intent(out) :: face_flux_y(:, :, 0:, :)
    real(dp), intent(out) :: face_flux_z(:, :, :, 0:)
    logical, intent(out) :: ok

    real(dp), allocatable :: extended_state(:, :, :, :)
    real(dp), allocatable :: extended_temperature(:, :, :)
    real(dp), allocatable :: extended_flux_x(:, :, :, :)
    real(dp), allocatable :: extended_flux_y(:, :, :, :)
    real(dp), allocatable :: extended_flux_z(:, :, :, :)
    real(dp), allocatable :: primitive(:)
    real(dp) :: ghost_temperature
    logical :: local_ok
    integer :: i, j, k, local_i, local_j, local_k
    integer :: coarse_i, coarse_j, coarse_k
    integer :: nx, ny, nz, ex, ey, ez, nvar, ratio
    integer :: extended_first, extended_last, fill_first, fill_last
    integer :: first_x_face

    face_flux_x = 0.0_dp
    face_flux_y = 0.0_dp
    face_flux_z = 0.0_dp
    ok = .false.
    nx = patch%fine_nx()
    ny = patch%fine_ny()
    nz = patch%fine_nz()
    ex = nx + 4
    ey = ny + 4
    ez = nz + 4
    nvar = reactive_nvar(size(species))
    ratio = patch%refinement_ratio
    if (.not. valid_fine_advance_shapes( &
          species, patch, coarse_start, coarse_start_temperature, &
          coarse_end, coarse_end_temperature, state, temperature, &
          face_flux_x, face_flux_y, face_flux_z)) return
    if (.not. ieee_is_finite(alpha)) return
    if (alpha < 0.0_dp .or. alpha > 1.0_dp) return
    if (first_i < 1 .or. last_i > nx .or. last_i < first_i) return
    if (trim(limiter) /= "minmod" .and. trim(limiter) /= "mc") return

    allocate(extended_state(nvar, ex, ey, ez))
    allocate(extended_temperature(ex, ey, ez))
    allocate(extended_flux_x(nvar, ex, ey, ez))
    allocate(extended_flux_y(nvar, ex, ey, ez))
    allocate(extended_flux_z(nvar, ex, ey, ez))
    allocate(primitive(reactive_nprim(size(species))))
    extended_state = 0.0_dp
    extended_temperature = 0.0_dp
    extended_first = first_i + 2
    if (first_i == 1) extended_first = 2
    extended_last = last_i + 2
    fill_first = max(1, extended_first - 1)
    fill_last = min(ex, extended_last + 2)
    do k = 1, ez
      local_k = k - 2
      coarse_k = patch%coarse_k_lower + &
        floor(real(local_k - 1, dp) / real(ratio, dp))
      do j = 1, ey
        local_j = j - 2
        coarse_j = patch%coarse_j_lower + &
          floor(real(local_j - 1, dp) / real(ratio, dp))
        do i = fill_first, fill_last
          local_i = i - 2
          if (i >= 3 .and. i <= nx + 2 .and. &
              j >= 3 .and. j <= ny + 2 .and. &
              k >= 3 .and. k <= nz + 2) then
            extended_state(:, i, j, k) = &
              state(:, local_i, local_j, local_k)
            extended_temperature(i, j, k) = &
              temperature(local_i, local_j, local_k)
            cycle
          end if
          coarse_i = patch%coarse_i_lower + &
            floor(real(local_i - 1, dp) / real(ratio, dp))
          call interpolate_coarse_fine_ghost_plm_3d( &
            species, coarse_start, coarse_start_temperature, &
            coarse_end, coarse_end_temperature, coarse_i, coarse_j, &
            coarse_k, local_i, local_j, local_k, ratio, alpha, limiter, &
            primitive, extended_state(:, i, j, k), ghost_temperature, &
            local_ok)
          if (.not. local_ok) return
          extended_temperature(i, j, k) = ghost_temperature
        end do
      end do
    end do
    call compute_reactive_plm_slab_face_fluxes_3d( &
      species, extended_state, extended_temperature, ex, ey, ez, &
      extended_first, extended_last, limiter, riemann_solver, &
      extended_flux_x, extended_flux_y, extended_flux_z, local_ok)
    if (.not. local_ok) return
    first_x_face = first_i
    if (first_i == 1) first_x_face = 0
    do k = 1, nz
      do j = 1, ny
        do i = first_x_face, last_i
          face_flux_x(:, i, j, k) = &
            extended_flux_x(:, i + 2, j + 2, k + 2)
        end do
      end do
    end do
    do k = 1, nz
      do j = 0, ny
        do i = first_i, last_i
          face_flux_y(:, i, j, k) = &
            extended_flux_y(:, i + 2, j + 2, k + 2)
        end do
      end do
    end do
    do k = 0, nz
      do j = 1, ny
        do i = first_i, last_i
          face_flux_z(:, i, j, k) = &
            extended_flux_z(:, i + 2, j + 2, k + 2)
        end do
      end do
    end do
    ok = all(ieee_is_finite(face_flux_x)) .and. &
      all(ieee_is_finite(face_flux_y)) .and. &
      all(ieee_is_finite(face_flux_z))
  end subroutine compute_fine_patch_plm_slab_face_fluxes_3d

  subroutine interpolate_coarse_fine_ghost_plm_3d( &
      species, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, coarse_i, coarse_j, coarse_k, &
      local_i, local_j, local_k, ratio, alpha, limiter, primitive, &
      state, temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: coarse_start(:, :, :, :)
    real(dp), intent(in) :: coarse_start_temperature(:, :, :)
    real(dp), intent(in) :: coarse_end(:, :, :, :)
    real(dp), intent(in) :: coarse_end_temperature(:, :, :)
    integer, intent(in) :: coarse_i, coarse_j, coarse_k
    integer, intent(in) :: local_i, local_j, local_k, ratio
    real(dp), intent(in) :: alpha
    character(len=*), intent(in) :: limiter
    real(dp), intent(out) :: primitive(:), state(:), temperature
    logical, intent(out) :: ok

    real(dp) :: center(size(primitive)), minus(size(primitive))
    real(dp) :: plus(size(primitive)), slope_x(size(primitive))
    real(dp) :: slope_y(size(primitive)), slope_z(size(primitive))
    real(dp) :: delta(size(primitive))
    real(dp) :: offset_x, offset_y, offset_z, theta, candidate_value
    real(dp) :: sound_speed, mass_fraction_sum
    logical :: local_ok
    integer :: component, species_index, im, ip, jm, jp, km, kp

    ok = .false.
    if (.not. valid_ghost_plm_inputs_3d( &
          species, coarse_start, coarse_start_temperature, coarse_end, &
          coarse_end_temperature, coarse_i, coarse_j, coarse_k, local_i, &
          local_j, local_k, ratio, alpha, limiter, primitive, state)) return
    im = 1 + modulo(coarse_i - 2, size(coarse_start, 2))
    ip = 1 + modulo(coarse_i, size(coarse_start, 2))
    jm = 1 + modulo(coarse_j - 2, size(coarse_start, 3))
    jp = 1 + modulo(coarse_j, size(coarse_start, 3))
    km = 1 + modulo(coarse_k - 2, size(coarse_start, 4))
    kp = 1 + modulo(coarse_k, size(coarse_start, 4))
    call time_interpolated_coarse_primitive_3d( &
      species, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, coarse_i, coarse_j, coarse_k, &
      alpha, center, local_ok)
    if (.not. local_ok) return
    call time_interpolated_coarse_primitive_3d( &
      species, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, im, coarse_j, coarse_k, &
      alpha, minus, local_ok)
    if (.not. local_ok) return
    call time_interpolated_coarse_primitive_3d( &
      species, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, ip, coarse_j, coarse_k, &
      alpha, plus, local_ok)
    if (.not. local_ok) return
    do component = 1, size(primitive)
      call limited_slope( &
        center(component) - minus(component), &
        plus(component) - center(component), limiter, &
        slope_x(component), local_ok)
      if (.not. local_ok) return
    end do
    call time_interpolated_coarse_primitive_3d( &
      species, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, coarse_i, jm, coarse_k, &
      alpha, minus, local_ok)
    if (.not. local_ok) return
    call time_interpolated_coarse_primitive_3d( &
      species, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, coarse_i, jp, coarse_k, &
      alpha, plus, local_ok)
    if (.not. local_ok) return
    do component = 1, size(primitive)
      call limited_slope( &
        center(component) - minus(component), &
        plus(component) - center(component), limiter, &
        slope_y(component), local_ok)
      if (.not. local_ok) return
    end do
    call time_interpolated_coarse_primitive_3d( &
      species, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, coarse_i, coarse_j, km, &
      alpha, minus, local_ok)
    if (.not. local_ok) return
    call time_interpolated_coarse_primitive_3d( &
      species, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, coarse_i, coarse_j, kp, &
      alpha, plus, local_ok)
    if (.not. local_ok) return
    do component = 1, size(primitive)
      call limited_slope( &
        center(component) - minus(component), &
        plus(component) - center(component), limiter, &
        slope_z(component), local_ok)
      if (.not. local_ok) return
    end do

    offset_x = (real(modulo(local_i - 1, ratio), dp) + 0.5_dp) / &
      real(ratio, dp) - 0.5_dp
    offset_y = (real(modulo(local_j - 1, ratio), dp) + 0.5_dp) / &
      real(ratio, dp) - 0.5_dp
    offset_z = (real(modulo(local_k - 1, ratio), dp) + 0.5_dp) / &
      real(ratio, dp) - 0.5_dp
    delta = offset_x * slope_x + offset_y * slope_y + offset_z * slope_z
    theta = positive_increment_scale_3d( &
      center(1), delta(1), density_floor)
    theta = min(theta, positive_increment_scale_3d( &
      center(5), delta(5), pressure_floor))
    do species_index = 1, size(species)
      component = reactive_mass_fraction_component(species_index)
      theta = min(theta, positive_increment_scale_3d( &
        center(component), delta(component), 0.0_dp))
    end do
    primitive = center + max(0.0_dp, min(1.0_dp, theta)) * delta
    mass_fraction_sum = 0.0_dp
    do species_index = 1, size(species)
      component = reactive_mass_fraction_component(species_index)
      primitive(component) = max(0.0_dp, primitive(component))
      mass_fraction_sum = mass_fraction_sum + primitive(component)
    end do
    if (mass_fraction_sum <= tiny(1.0_dp)) return
    do species_index = 1, size(species)
      component = reactive_mass_fraction_component(species_index)
      primitive(component) = primitive(component) / mass_fraction_sum
    end do
    call reactive_primitive_to_conserved( &
      species, primitive, state, temperature, sound_speed, local_ok)
    if (.not. local_ok) then
      primitive = center
      call reactive_primitive_to_conserved( &
        species, primitive, state, temperature, sound_speed, local_ok)
    end if
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    candidate_value = temperature + sound_speed
    ok = ieee_is_finite(candidate_value)
  end subroutine interpolate_coarse_fine_ghost_plm_3d

  subroutine time_interpolated_coarse_primitive_3d( &
      species, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, i, j, k, alpha, primitive, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: coarse_start(:, :, :, :)
    real(dp), intent(in) :: coarse_start_temperature(:, :, :)
    real(dp), intent(in) :: coarse_end(:, :, :, :)
    real(dp), intent(in) :: coarse_end_temperature(:, :, :)
    integer, intent(in) :: i, j, k
    real(dp), intent(in) :: alpha
    real(dp), intent(out) :: primitive(:)
    logical, intent(out) :: ok

    real(dp) :: interpolated_state(size(coarse_start, 1))
    real(dp) :: temperature_guess, temperature, sound_speed

    interpolated_state = &
      (1.0_dp - alpha) * coarse_start(:, i, j, k) + &
      alpha * coarse_end(:, i, j, k)
    temperature_guess = &
      (1.0_dp - alpha) * coarse_start_temperature(i, j, k) + &
      alpha * coarse_end_temperature(i, j, k)
    call reactive_conserved_to_primitive( &
      species, interpolated_state, temperature_guess, primitive, &
      temperature, sound_speed, ok)
  end subroutine time_interpolated_coarse_primitive_3d

  pure real(dp) function positive_increment_scale_3d( &
      center, delta, lower_bound) result(theta)
    real(dp), intent(in) :: center, delta, lower_bound

    theta = 1.0_dp
    if (delta < 0.0_dp) then
      theta = min(1.0_dp, max(0.0_dp, &
        (center - lower_bound) / max(-delta, tiny(1.0_dp))))
    end if
  end function positive_increment_scale_3d

  subroutine interpolate_coarse_cell( &
      species, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, i, j, k, alpha, primitive, &
      state, temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: coarse_start(:, :, :, :)
    real(dp), intent(in) :: coarse_start_temperature(:, :, :)
    real(dp), intent(in) :: coarse_end(:, :, :, :)
    real(dp), intent(in) :: coarse_end_temperature(:, :, :)
    integer, intent(in) :: i, j, k
    real(dp), intent(in) :: alpha
    real(dp), intent(out) :: primitive(:), state(:), temperature
    logical, intent(out) :: ok

    real(dp) :: sound_speed, temperature_guess

    state = (1.0_dp - alpha) * coarse_start(:, i, j, k) + &
      alpha * coarse_end(:, i, j, k)
    temperature_guess = &
      (1.0_dp - alpha) * coarse_start_temperature(i, j, k) + &
      alpha * coarse_end_temperature(i, j, k)
    call reactive_conserved_to_primitive( &
      species, state, temperature_guess, primitive, temperature, &
      sound_speed, ok)
  end subroutine interpolate_coarse_cell

  subroutine fine_flux_divergence_3d( &
      face_flux_x, face_flux_y, face_flux_z, dx, dy, dz, rhs, ok)
    real(dp), intent(in) :: face_flux_x(:, 0:, :, :)
    real(dp), intent(in) :: face_flux_y(:, :, 0:, :)
    real(dp), intent(in) :: face_flux_z(:, :, :, 0:)
    real(dp), intent(in) :: dx, dy, dz
    real(dp), intent(out) :: rhs(:, :, :, :)
    logical, intent(out) :: ok

    integer :: i, j, k, nx, ny, nz

    rhs = 0.0_dp
    nx = size(rhs, 2)
    ny = size(rhs, 3)
    nz = size(rhs, 4)
    ok = all(ieee_is_finite([dx, dy, dz]))
    if (.not. ok) return
    ok = dx > 0.0_dp .and. dy > 0.0_dp .and. dz > 0.0_dp .and. &
      size(face_flux_x, 1) == size(rhs, 1) .and. &
      size(face_flux_x, 2) == nx + 1 .and. &
      size(face_flux_x, 3) == ny .and. size(face_flux_x, 4) == nz .and. &
      size(face_flux_y, 1) == size(rhs, 1) .and. &
      size(face_flux_y, 2) == nx .and. &
      size(face_flux_y, 3) == ny + 1 .and. size(face_flux_y, 4) == nz .and. &
      size(face_flux_z, 1) == size(rhs, 1) .and. &
      size(face_flux_z, 2) == nx .and. size(face_flux_z, 3) == ny .and. &
      size(face_flux_z, 4) == nz + 1
    if (.not. ok) return
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          rhs(:, i, j, k) = &
            -(face_flux_x(:, i, j, k) - &
              face_flux_x(:, i - 1, j, k)) / dx &
            -(face_flux_y(:, i, j, k) - &
              face_flux_y(:, i, j - 1, k)) / dy &
            -(face_flux_z(:, i, j, k) - &
              face_flux_z(:, i, j, k - 1)) / dz
        end do
      end do
    end do
    ok = all(ieee_is_finite(rhs))
  end subroutine fine_flux_divergence_3d

  subroutine accumulate_fine_interface_fluxes_3d( &
      patch, face_flux_x, face_flux_y, face_flux_z, &
      x_lower, x_upper, y_lower, y_upper, z_lower, z_upper, ok)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: face_flux_x(:, 0:, :, :)
    real(dp), intent(in) :: face_flux_y(:, :, 0:, :)
    real(dp), intent(in) :: face_flux_z(:, :, :, 0:)
    real(dp), intent(inout) :: x_lower(:, :, :), x_upper(:, :, :)
    real(dp), intent(inout) :: y_lower(:, :, :), y_upper(:, :, :)
    real(dp), intent(inout) :: z_lower(:, :, :), z_upper(:, :, :)
    logical, intent(out) :: ok

    real(dp) :: inverse_face_children
    real(dp), allocatable :: candidate_x_lower(:, :, :), candidate_x_upper(:, :, :)
    real(dp), allocatable :: candidate_y_lower(:, :, :), candidate_y_upper(:, :, :)
    real(dp), allocatable :: candidate_z_lower(:, :, :), candidate_z_upper(:, :, :)
    integer :: i, j, k, fine_i, fine_j, fine_k
    integer :: nx, ny, nz, ratio

    ok = .false.
    ratio = patch%refinement_ratio
    nx = patch%fine_nx()
    ny = patch%fine_ny()
    nz = patch%fine_nz()
    if (.not. valid_interface_flux_shapes( &
          patch, face_flux_x, face_flux_y, face_flux_z, &
          x_lower, x_upper, y_lower, y_upper, z_lower, z_upper)) return
    allocate(candidate_x_lower, source=x_lower)
    allocate(candidate_x_upper, source=x_upper)
    allocate(candidate_y_lower, source=y_lower)
    allocate(candidate_y_upper, source=y_upper)
    allocate(candidate_z_lower, source=z_lower)
    allocate(candidate_z_upper, source=z_upper)
    inverse_face_children = 1.0_dp / real(ratio**2, dp)
    do k = 1, size(x_lower, 3)
      do j = 1, size(x_lower, 2)
        do fine_k = (k - 1) * ratio + 1, k * ratio
          do fine_j = (j - 1) * ratio + 1, j * ratio
            candidate_x_lower(:, j, k) = candidate_x_lower(:, j, k) + &
              inverse_face_children * face_flux_x(:, 0, fine_j, fine_k)
            candidate_x_upper(:, j, k) = candidate_x_upper(:, j, k) + &
              inverse_face_children * face_flux_x(:, nx, fine_j, fine_k)
          end do
        end do
      end do
    end do
    do k = 1, size(y_lower, 3)
      do i = 1, size(y_lower, 2)
        do fine_k = (k - 1) * ratio + 1, k * ratio
          do fine_i = (i - 1) * ratio + 1, i * ratio
            candidate_y_lower(:, i, k) = candidate_y_lower(:, i, k) + &
              inverse_face_children * face_flux_y(:, fine_i, 0, fine_k)
            candidate_y_upper(:, i, k) = candidate_y_upper(:, i, k) + &
              inverse_face_children * face_flux_y(:, fine_i, ny, fine_k)
          end do
        end do
      end do
    end do
    do j = 1, size(z_lower, 3)
      do i = 1, size(z_lower, 2)
        do fine_j = (j - 1) * ratio + 1, j * ratio
          do fine_i = (i - 1) * ratio + 1, i * ratio
            candidate_z_lower(:, i, j) = candidate_z_lower(:, i, j) + &
              inverse_face_children * face_flux_z(:, fine_i, fine_j, 0)
            candidate_z_upper(:, i, j) = candidate_z_upper(:, i, j) + &
              inverse_face_children * face_flux_z(:, fine_i, fine_j, nz)
          end do
        end do
      end do
    end do
    ok = all(ieee_is_finite(candidate_x_lower)) .and. &
      all(ieee_is_finite(candidate_x_upper)) .and. &
      all(ieee_is_finite(candidate_y_lower)) .and. &
      all(ieee_is_finite(candidate_y_upper)) .and. &
      all(ieee_is_finite(candidate_z_lower)) .and. &
      all(ieee_is_finite(candidate_z_upper))
    if (.not. ok) return
    x_lower = candidate_x_lower
    x_upper = candidate_x_upper
    y_lower = candidate_y_lower
    y_upper = candidate_y_upper
    z_lower = candidate_z_lower
    z_upper = candidate_z_upper
  end subroutine accumulate_fine_interface_fluxes_3d

  subroutine reflux_coarse_3d( &
      patch, coarse_state, coarse_flux_x, coarse_flux_y, coarse_flux_z, &
      fine_x_lower, fine_x_upper, fine_y_lower, fine_y_upper, &
      fine_z_lower, fine_z_upper, dx, dy, dz, dt, &
      maximum_correction, ok)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(inout) :: coarse_state(:, :, :, :)
    real(dp), intent(in) :: coarse_flux_x(:, :, :, :)
    real(dp), intent(in) :: coarse_flux_y(:, :, :, :)
    real(dp), intent(in) :: coarse_flux_z(:, :, :, :)
    real(dp), intent(in) :: fine_x_lower(:, :, :), fine_x_upper(:, :, :)
    real(dp), intent(in) :: fine_y_lower(:, :, :), fine_y_upper(:, :, :)
    real(dp), intent(in) :: fine_z_lower(:, :, :), fine_z_upper(:, :, :)
    real(dp), intent(in) :: dx, dy, dz, dt
    real(dp), intent(out) :: maximum_correction
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate_state(:, :, :, :), correction(:)
    real(dp) :: candidate_maximum_correction
    integer :: i, j, k, local_i, local_j, local_k

    maximum_correction = 0.0_dp
    ok = .false.
    if (.not. valid_reflux_inputs_3d( &
          patch, coarse_state, coarse_flux_x, coarse_flux_y, coarse_flux_z, &
          fine_x_lower, fine_x_upper, fine_y_lower, fine_y_upper, &
          fine_z_lower, fine_z_upper, dx, dy, dz, dt)) return
    allocate(candidate_state, source=coarse_state)
    allocate(correction(size(coarse_state, 1)))
    candidate_maximum_correction = 0.0_dp

    do k = patch%coarse_k_lower, patch%coarse_k_upper
      local_k = k - patch%coarse_k_lower + 1
      do j = patch%coarse_j_lower, patch%coarse_j_upper
        local_j = j - patch%coarse_j_lower + 1
        correction = -dt / dx * (fine_x_lower(:, local_j, local_k) - &
          coarse_flux_x(:, patch%coarse_i_lower - 1, j, k))
        candidate_state(:, patch%coarse_i_lower - 1, j, k) = &
          candidate_state(:, patch%coarse_i_lower - 1, j, k) + correction
        candidate_maximum_correction = max(candidate_maximum_correction, &
          maxval(abs(correction)))
        correction = dt / dx * (fine_x_upper(:, local_j, local_k) - &
          coarse_flux_x(:, patch%coarse_i_upper, j, k))
        candidate_state(:, patch%coarse_i_upper + 1, j, k) = &
          candidate_state(:, patch%coarse_i_upper + 1, j, k) + correction
        candidate_maximum_correction = max(candidate_maximum_correction, &
          maxval(abs(correction)))
      end do
    end do
    do k = patch%coarse_k_lower, patch%coarse_k_upper
      local_k = k - patch%coarse_k_lower + 1
      do i = patch%coarse_i_lower, patch%coarse_i_upper
        local_i = i - patch%coarse_i_lower + 1
        correction = -dt / dy * (fine_y_lower(:, local_i, local_k) - &
          coarse_flux_y(:, i, patch%coarse_j_lower - 1, k))
        candidate_state(:, i, patch%coarse_j_lower - 1, k) = &
          candidate_state(:, i, patch%coarse_j_lower - 1, k) + correction
        candidate_maximum_correction = max(candidate_maximum_correction, &
          maxval(abs(correction)))
        correction = dt / dy * (fine_y_upper(:, local_i, local_k) - &
          coarse_flux_y(:, i, patch%coarse_j_upper, k))
        candidate_state(:, i, patch%coarse_j_upper + 1, k) = &
          candidate_state(:, i, patch%coarse_j_upper + 1, k) + correction
        candidate_maximum_correction = max(candidate_maximum_correction, &
          maxval(abs(correction)))
      end do
    end do
    do j = patch%coarse_j_lower, patch%coarse_j_upper
      local_j = j - patch%coarse_j_lower + 1
      do i = patch%coarse_i_lower, patch%coarse_i_upper
        local_i = i - patch%coarse_i_lower + 1
        correction = -dt / dz * (fine_z_lower(:, local_i, local_j) - &
          coarse_flux_z(:, i, j, patch%coarse_k_lower - 1))
        candidate_state(:, i, j, patch%coarse_k_lower - 1) = &
          candidate_state(:, i, j, patch%coarse_k_lower - 1) + correction
        candidate_maximum_correction = max(candidate_maximum_correction, &
          maxval(abs(correction)))
        correction = dt / dz * (fine_z_upper(:, local_i, local_j) - &
          coarse_flux_z(:, i, j, patch%coarse_k_upper))
        candidate_state(:, i, j, patch%coarse_k_upper + 1) = &
          candidate_state(:, i, j, patch%coarse_k_upper + 1) + correction
        candidate_maximum_correction = max(candidate_maximum_correction, &
          maxval(abs(correction)))
      end do
    end do
    ok = ieee_is_finite(candidate_maximum_correction) .and. &
      all(ieee_is_finite(candidate_state))
    if (.not. ok) return
    coarse_state = candidate_state
    maximum_correction = candidate_maximum_correction
  end subroutine reflux_coarse_3d

  subroutine dispatch_amr_reactive_hydro_3d( &
      species, patch, coarse_state, coarse_temperature, fine_state, &
      fine_temperature, dx, dy, dz, dt, riemann_solver, &
      maximum_reflux_correction, ok, reconstruction, limiter)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(inout) :: coarse_state(:, :, :, :)
    real(dp), intent(inout) :: coarse_temperature(:, :, :)
    real(dp), intent(inout) :: fine_state(:, :, :, :)
    real(dp), intent(inout) :: fine_temperature(:, :, :)
    real(dp), intent(in) :: dx, dy, dz, dt
    character(len=*), intent(in) :: riemann_solver
    real(dp), intent(out) :: maximum_reflux_correction
    logical, intent(out) :: ok
    character(len=*), intent(in), optional :: reconstruction, limiter

    if (present(reconstruction)) then
      if (present(limiter)) then
        call advance_amr_reactive_hydro_3d( &
          species, patch, coarse_state, coarse_temperature, fine_state, &
          fine_temperature, dx, dy, dz, dt, riemann_solver, &
          maximum_reflux_correction, ok, reconstruction=reconstruction, &
          limiter=limiter)
      else
        call advance_amr_reactive_hydro_3d( &
          species, patch, coarse_state, coarse_temperature, fine_state, &
          fine_temperature, dx, dy, dz, dt, riemann_solver, &
          maximum_reflux_correction, ok, reconstruction=reconstruction)
      end if
    else if (present(limiter)) then
      call advance_amr_reactive_hydro_3d( &
        species, patch, coarse_state, coarse_temperature, fine_state, &
        fine_temperature, dx, dy, dz, dt, riemann_solver, &
        maximum_reflux_correction, ok, limiter=limiter)
    else
      call advance_amr_reactive_hydro_3d( &
        species, patch, coarse_state, coarse_temperature, fine_state, &
        fine_temperature, dx, dy, dz, dt, riemann_solver, &
        maximum_reflux_correction, ok)
    end if
  end subroutine dispatch_amr_reactive_hydro_3d

  pure logical function valid_ghost_plm_inputs_3d( &
      species, coarse_start, coarse_start_temperature, coarse_end, &
      coarse_end_temperature, coarse_i, coarse_j, coarse_k, local_i, &
      local_j, local_k, ratio, alpha, limiter, primitive, state) result(valid)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: coarse_start(:, :, :, :)
    real(dp), intent(in) :: coarse_start_temperature(:, :, :)
    real(dp), intent(in) :: coarse_end(:, :, :, :)
    real(dp), intent(in) :: coarse_end_temperature(:, :, :)
    integer, intent(in) :: coarse_i, coarse_j, coarse_k
    integer, intent(in) :: local_i, local_j, local_k, ratio
    real(dp), intent(in) :: alpha
    character(len=*), intent(in) :: limiter
    real(dp), intent(in) :: primitive(:), state(:)

    integer :: nvar, nprim, nx, ny, nz

    valid = .false.
    if (size(species) < 1 .or. ratio < 2) return
    nvar = reactive_nvar(size(species))
    nprim = reactive_nprim(size(species))
    nx = size(coarse_start, 2)
    ny = size(coarse_start, 3)
    nz = size(coarse_start, 4)
    if (nx < 2 .or. ny < 2 .or. nz < 2) return
    if (size(coarse_start, 1) /= nvar .or. &
        size(coarse_end, 1) /= nvar .or. &
        size(coarse_end, 2) /= nx .or. size(coarse_end, 3) /= ny .or. &
        size(coarse_end, 4) /= nz .or. &
        any(shape(coarse_start_temperature) /= [nx, ny, nz]) .or. &
        any(shape(coarse_end_temperature) /= [nx, ny, nz]) .or. &
        size(primitive) /= nprim .or. size(state) /= nvar) return
    if (coarse_i < 1 .or. coarse_i > nx .or. &
        coarse_j < 1 .or. coarse_j > ny .or. &
        coarse_k < 1 .or. coarse_k > nz) return
    if (local_i < -1 .or. local_i > ratio * nx + 2 .or. &
        local_j < -1 .or. local_j > ratio * ny + 2 .or. &
        local_k < -1 .or. local_k > ratio * nz + 2) return
    if (.not. ieee_is_finite(alpha)) return
    if (alpha < 0.0_dp .or. alpha > 1.0_dp) return
    if (trim(limiter) /= "minmod" .and. trim(limiter) /= "mc") return
    valid = .true.
  end function valid_ghost_plm_inputs_3d

  pure logical function valid_reflux_inputs_3d( &
      patch, coarse_state, coarse_flux_x, coarse_flux_y, coarse_flux_z, &
      fine_x_lower, fine_x_upper, fine_y_lower, fine_y_upper, &
      fine_z_lower, fine_z_upper, dx, dy, dz, dt) result(valid)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_state(:, :, :, :)
    real(dp), intent(in) :: coarse_flux_x(:, :, :, :)
    real(dp), intent(in) :: coarse_flux_y(:, :, :, :)
    real(dp), intent(in) :: coarse_flux_z(:, :, :, :)
    real(dp), intent(in) :: fine_x_lower(:, :, :), fine_x_upper(:, :, :)
    real(dp), intent(in) :: fine_y_lower(:, :, :), fine_y_upper(:, :, :)
    real(dp), intent(in) :: fine_z_lower(:, :, :), fine_z_upper(:, :, :)
    real(dp), intent(in) :: dx, dy, dz, dt

    integer :: nvar, covered_nx, covered_ny, covered_nz

    valid = .false.
    if (.not. patch%is_strictly_interior()) return
    nvar = size(coarse_state, 1)
    covered_nx = patch%coarse_i_upper - patch%coarse_i_lower + 1
    covered_ny = patch%coarse_j_upper - patch%coarse_j_lower + 1
    covered_nz = patch%coarse_k_upper - patch%coarse_k_lower + 1
    if (nvar < 1 .or. size(coarse_state, 2) /= patch%coarse_nx .or. &
        size(coarse_state, 3) /= patch%coarse_ny .or. &
        size(coarse_state, 4) /= patch%coarse_nz .or. &
        any(shape(coarse_flux_x) /= shape(coarse_state)) .or. &
        any(shape(coarse_flux_y) /= shape(coarse_state)) .or. &
        any(shape(coarse_flux_z) /= shape(coarse_state)) .or. &
        any(shape(fine_x_lower) /= [nvar, covered_ny, covered_nz]) .or. &
        any(shape(fine_x_upper) /= shape(fine_x_lower)) .or. &
        any(shape(fine_y_lower) /= [nvar, covered_nx, covered_nz]) .or. &
        any(shape(fine_y_upper) /= shape(fine_y_lower)) .or. &
        any(shape(fine_z_lower) /= [nvar, covered_nx, covered_ny]) .or. &
        any(shape(fine_z_upper) /= shape(fine_z_lower))) return
    if (.not. all(ieee_is_finite([dx, dy, dz, dt]))) return
    if (dx <= 0.0_dp .or. dy <= 0.0_dp .or. dz <= 0.0_dp .or. &
        dt <= 0.0_dp) return
    if (.not. all(ieee_is_finite(coarse_state)) .or. &
        .not. all(ieee_is_finite(coarse_flux_x)) .or. &
        .not. all(ieee_is_finite(coarse_flux_y)) .or. &
        .not. all(ieee_is_finite(coarse_flux_z)) .or. &
        .not. all(ieee_is_finite(fine_x_lower)) .or. &
        .not. all(ieee_is_finite(fine_x_upper)) .or. &
        .not. all(ieee_is_finite(fine_y_lower)) .or. &
        .not. all(ieee_is_finite(fine_y_upper)) .or. &
        .not. all(ieee_is_finite(fine_z_lower)) .or. &
        .not. all(ieee_is_finite(fine_z_upper))) return
    valid = .true.
  end function valid_reflux_inputs_3d

  logical function valid_amr_reactive_chemistry_inputs( &
      species, reactions, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, interval, rtol, atol) result(valid)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_state(:, :, :, :)
    real(dp), intent(in) :: coarse_temperature(:, :, :)
    real(dp), intent(in) :: fine_state(:, :, :, :)
    real(dp), intent(in) :: fine_temperature(:, :, :)
    real(dp), intent(in) :: interval, rtol, atol

    integer :: reaction_index

    valid = .false.
    if (.not. patch%is_strictly_interior()) return
    if (.not. valid_amr_reactive_shapes( &
          species, patch, coarse_state, coarse_temperature, &
          fine_state, fine_temperature)) return
    if (size(reactions) < 1) return
    if (.not. all(ieee_is_finite([interval, rtol, atol]))) return
    if (interval < 0.0_dp .or. rtol <= 0.0_dp .or. atol <= 0.0_dp) return
    do reaction_index = 1, size(reactions)
      if (.not. valid_elementary_reaction( &
            reactions(reaction_index), size(species))) return
    end do
    valid = .true.
  end function valid_amr_reactive_chemistry_inputs

  logical function valid_amr_reactive_chemistry_integrator( &
      chemistry_integrator) result(valid)
    character(len=*), intent(in) :: chemistry_integrator

    valid = trim(chemistry_integrator) == "explicit" .or. &
      trim(chemistry_integrator) == "implicit"
  end function valid_amr_reactive_chemistry_integrator

  pure logical function valid_amr_reactive_shapes( &
      species, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature) result(valid)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_state(:, :, :, :)
    real(dp), intent(in) :: coarse_temperature(:, :, :)
    real(dp), intent(in) :: fine_state(:, :, :, :)
    real(dp), intent(in) :: fine_temperature(:, :, :)

    integer :: nvar

    nvar = reactive_nvar(size(species))
    valid = patch%is_valid() .and. size(coarse_state, 1) == nvar .and. &
      size(coarse_state, 2) == patch%coarse_nx .and. &
      size(coarse_state, 3) == patch%coarse_ny .and. &
      size(coarse_state, 4) == patch%coarse_nz .and. &
      size(coarse_temperature, 1) == patch%coarse_nx .and. &
      size(coarse_temperature, 2) == patch%coarse_ny .and. &
      size(coarse_temperature, 3) == patch%coarse_nz .and. &
      size(fine_state, 1) == nvar .and. &
      size(fine_state, 2) == patch%fine_nx() .and. &
      size(fine_state, 3) == patch%fine_ny() .and. &
      size(fine_state, 4) == patch%fine_nz() .and. &
      size(fine_temperature, 1) == patch%fine_nx() .and. &
      size(fine_temperature, 2) == patch%fine_ny() .and. &
      size(fine_temperature, 3) == patch%fine_nz() .and. &
      all(ieee_is_finite(coarse_state)) .and. &
      all(ieee_is_finite(coarse_temperature)) .and. &
      all(ieee_is_finite(fine_state)) .and. &
      all(ieee_is_finite(fine_temperature))
  end function valid_amr_reactive_shapes

  pure logical function valid_fine_advance_shapes( &
      species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, state, temperature, &
      face_flux_x, face_flux_y, face_flux_z) result(valid)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_start(:, :, :, :)
    real(dp), intent(in) :: coarse_start_temperature(:, :, :)
    real(dp), intent(in) :: coarse_end(:, :, :, :)
    real(dp), intent(in) :: coarse_end_temperature(:, :, :)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(in) :: face_flux_x(:, 0:, :, :)
    real(dp), intent(in) :: face_flux_y(:, :, 0:, :)
    real(dp), intent(in) :: face_flux_z(:, :, :, 0:)

    valid = valid_amr_reactive_shapes( &
      species, patch, coarse_start, coarse_start_temperature, &
      state, temperature) .and. &
      all(shape(coarse_end) == shape(coarse_start)) .and. &
      all(shape(coarse_end_temperature) == &
        shape(coarse_start_temperature)) .and. &
      all(ieee_is_finite(coarse_end)) .and. &
      all(ieee_is_finite(coarse_end_temperature)) .and. &
      size(face_flux_x, 1) == size(state, 1) .and. &
      size(face_flux_x, 2) == size(state, 2) + 1 .and. &
      size(face_flux_x, 3) == size(state, 3) .and. &
      size(face_flux_x, 4) == size(state, 4) .and. &
      size(face_flux_y, 1) == size(state, 1) .and. &
      size(face_flux_y, 2) == size(state, 2) .and. &
      size(face_flux_y, 3) == size(state, 3) + 1 .and. &
      size(face_flux_y, 4) == size(state, 4) .and. &
      size(face_flux_z, 1) == size(state, 1) .and. &
      size(face_flux_z, 2) == size(state, 2) .and. &
      size(face_flux_z, 3) == size(state, 3) .and. &
      size(face_flux_z, 4) == size(state, 4) + 1
  end function valid_fine_advance_shapes

  pure logical function valid_interface_flux_shapes( &
      patch, face_flux_x, face_flux_y, face_flux_z, &
      x_lower, x_upper, y_lower, y_upper, z_lower, z_upper) result(valid)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: face_flux_x(:, 0:, :, :)
    real(dp), intent(in) :: face_flux_y(:, :, 0:, :)
    real(dp), intent(in) :: face_flux_z(:, :, :, 0:)
    real(dp), intent(in) :: x_lower(:, :, :), x_upper(:, :, :)
    real(dp), intent(in) :: y_lower(:, :, :), y_upper(:, :, :)
    real(dp), intent(in) :: z_lower(:, :, :), z_upper(:, :, :)

    integer :: covered_nx, covered_ny, covered_nz, nvar

    covered_nx = patch%coarse_i_upper - patch%coarse_i_lower + 1
    covered_ny = patch%coarse_j_upper - patch%coarse_j_lower + 1
    covered_nz = patch%coarse_k_upper - patch%coarse_k_lower + 1
    nvar = size(face_flux_x, 1)
    valid = patch%is_valid() .and. nvar > 0 .and. &
      size(face_flux_x, 2) == patch%fine_nx() + 1 .and. &
      size(face_flux_x, 3) == patch%fine_ny() .and. &
      size(face_flux_x, 4) == patch%fine_nz() .and. &
      size(face_flux_y, 1) == nvar .and. &
      size(face_flux_y, 2) == patch%fine_nx() .and. &
      size(face_flux_y, 3) == patch%fine_ny() + 1 .and. &
      size(face_flux_y, 4) == patch%fine_nz() .and. &
      size(face_flux_z, 1) == nvar .and. &
      size(face_flux_z, 2) == patch%fine_nx() .and. &
      size(face_flux_z, 3) == patch%fine_ny() .and. &
      size(face_flux_z, 4) == patch%fine_nz() + 1 .and. &
      all(shape(x_lower) == [nvar, covered_ny, covered_nz]) .and. &
      all(shape(x_upper) == shape(x_lower)) .and. &
      all(shape(y_lower) == [nvar, covered_nx, covered_nz]) .and. &
      all(shape(y_upper) == shape(y_lower)) .and. &
      all(shape(z_lower) == [nvar, covered_nx, covered_ny]) .and. &
      all(shape(z_upper) == shape(z_lower)) .and. &
      all(ieee_is_finite(face_flux_x)) .and. &
      all(ieee_is_finite(face_flux_y)) .and. &
      all(ieee_is_finite(face_flux_z)) .and. &
      all(ieee_is_finite(x_lower)) .and. &
      all(ieee_is_finite(x_upper)) .and. &
      all(ieee_is_finite(y_lower)) .and. &
      all(ieee_is_finite(y_upper)) .and. &
      all(ieee_is_finite(z_lower)) .and. &
      all(ieee_is_finite(z_upper))
  end function valid_interface_flux_shapes

end module amr_reactive_3d_mod
