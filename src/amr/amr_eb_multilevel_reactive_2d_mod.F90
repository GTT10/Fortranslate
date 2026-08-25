module amr_eb_multilevel_reactive_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use state_indices_mod, only: irho, iet
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_conserved_to_primitive
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_covered_cell
  use eb_reactive_reconstruction_2d_mod, only: &
    reactive_eb_exterior_state_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, average_down_reactive_eb_state_patch_2d, &
    composite_eb_integral_2d
  use amr_eb_multilevel_2d_mod, only: &
    average_down_three_level_reactive_eb_state_2d
  use amr_eb_flux_register_2d_mod, only: &
    amr_eb_flux_register_2d, initialize_amr_eb_flux_register_2d, &
    accumulate_coarse_eb_fluxes_2d, accumulate_fine_eb_fluxes_2d, &
    reflux_reactive_eb_state_patch_2d
  use amr_eb_reactive_2d_mod, only: &
    build_reactive_eb_patch_exterior_2d, advance_reactive_eb_level_2d
  implicit none
  private

  public :: advance_three_level_reactive_eb_hydro_2d
  public :: level_two_interface_is_regular
  public :: close_cut_interface_conservation_2d

contains

  subroutine advance_three_level_reactive_eb_hydro_2d( &
      species, root_state, root_temperature, root_geometry, &
      level_one_state, level_one_temperature, level_one_geometry, &
      root_patch, level_two_state, level_two_temperature, &
      level_two_geometry, level_one_patch, solver, reconstruction, limiter, &
      state_redist_max_order, dt, new_root_state, new_root_temperature, &
      new_level_one_state, new_level_one_temperature, new_level_two_state, &
      new_level_two_temperature, ok, state_redist_target_volume_fraction)
    type(nasa7_species), intent(in) :: species(:)
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
    character(len=*), intent(in) :: solver, reconstruction, limiter
    integer, intent(in) :: state_redist_max_order
    real(dp), intent(in) :: dt
    real(dp), intent(out) :: new_root_state(:, :, :)
    real(dp), intent(out) :: new_root_temperature(:, :)
    real(dp), intent(out) :: new_level_one_state(:, :, :)
    real(dp), intent(out) :: new_level_one_temperature(:, :)
    real(dp), intent(out) :: new_level_two_state(:, :, :)
    real(dp), intent(out) :: new_level_two_temperature(:, :)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: state_redist_target_volume_fraction

    type(amr_eb_flux_register_2d) :: root_register, level_one_register
    type(reactive_eb_exterior_state_2d) :: level_one_exterior
    type(reactive_eb_exterior_state_2d) :: level_two_exterior
    real(dp), allocatable :: root_candidate(:, :, :)
    real(dp), allocatable :: root_candidate_temperature(:, :)
    real(dp), allocatable :: root_refluxed(:, :, :)
    real(dp), allocatable :: root_refluxed_temperature(:, :)
    real(dp), allocatable :: level_one_candidate(:, :, :)
    real(dp), allocatable :: level_one_candidate_temperature(:, :)
    real(dp), allocatable :: level_one_start(:, :, :)
    real(dp), allocatable :: level_one_start_temperature(:, :)
    real(dp), allocatable :: level_one_uncorrected(:, :, :)
    real(dp), allocatable :: level_one_uncorrected_temperature(:, :)
    real(dp), allocatable :: level_one_refluxed(:, :, :)
    real(dp), allocatable :: level_one_refluxed_temperature(:, :)
    real(dp), allocatable :: level_one_closed(:, :, :)
    real(dp), allocatable :: level_one_closed_temperature(:, :)
    real(dp), allocatable :: level_two_candidate(:, :, :)
    real(dp), allocatable :: level_two_candidate_temperature(:, :)
    real(dp), allocatable :: level_two_work(:, :, :)
    real(dp), allocatable :: level_two_work_temperature(:, :)
    real(dp), allocatable :: level_two_refluxed(:, :, :)
    real(dp), allocatable :: level_two_refluxed_temperature(:, :)
    real(dp), allocatable :: root_x_flux(:, :, :), root_y_flux(:, :, :)
    real(dp), allocatable :: level_one_x_flux(:, :, :)
    real(dp), allocatable :: level_one_y_flux(:, :, :)
    real(dp), allocatable :: level_two_x_flux(:, :, :)
    real(dp), allocatable :: level_two_y_flux(:, :, :)
    real(dp), allocatable :: level_one_integral_before(:)
    real(dp) :: level_one_dt, level_two_dt, alpha, selected_target
    logical :: local_ok
    integer :: nvar, level_one_ratio, level_two_ratio
    integer :: level_one_substep, level_two_substep

    new_root_state = root_state
    new_root_temperature = root_temperature
    new_level_one_state = level_one_state
    new_level_one_temperature = level_one_temperature
    new_level_two_state = level_two_state
    new_level_two_temperature = level_two_temperature
    ok = .false.
    nvar = reactive_nvar(size(species))
    selected_target = 0.5_dp
    if (present(state_redist_target_volume_fraction)) &
      selected_target = state_redist_target_volume_fraction
    if (nvar < 1 .or. .not. ieee_is_finite(dt) .or. dt <= 0.0_dp .or. &
        .not. ieee_is_finite(selected_target) .or. &
        selected_target <= 0.0_dp .or. selected_target > 1.0_dp .or. &
        .not. root_patch%is_valid(root_geometry, level_one_geometry) .or. &
        .not. level_one_patch%is_valid( &
          level_one_geometry, level_two_geometry) .or. &
        .not. level_two_patch_is_separated( &
          level_one_patch, level_one_geometry) .or. &
        any(shape(root_state) /= &
          [nvar, root_geometry%nx, root_geometry%ny]) .or. &
        any(shape(root_temperature) /= &
          [root_geometry%nx, root_geometry%ny]) .or. &
        any(shape(level_one_state) /= &
          [nvar, level_one_geometry%nx, level_one_geometry%ny]) .or. &
        any(shape(level_one_temperature) /= &
          [level_one_geometry%nx, level_one_geometry%ny]) .or. &
        any(shape(level_two_state) /= &
          [nvar, level_two_geometry%nx, level_two_geometry%ny]) .or. &
        any(shape(level_two_temperature) /= &
          [level_two_geometry%nx, level_two_geometry%ny]) .or. &
        any(shape(new_root_state) /= shape(root_state)) .or. &
        any(shape(new_root_temperature) /= shape(root_temperature)) .or. &
        any(shape(new_level_one_state) /= shape(level_one_state)) .or. &
        any(shape(new_level_one_temperature) /= &
          shape(level_one_temperature)) .or. &
        any(shape(new_level_two_state) /= shape(level_two_state)) .or. &
        any(shape(new_level_two_temperature) /= &
          shape(level_two_temperature)) .or. &
        any(.not. ieee_is_finite(root_state)) .or. &
        any(.not. ieee_is_finite(root_temperature)) .or. &
        any(.not. ieee_is_finite(level_one_state)) .or. &
        any(.not. ieee_is_finite(level_one_temperature)) .or. &
        any(.not. ieee_is_finite(level_two_state)) .or. &
        any(.not. ieee_is_finite(level_two_temperature)) .or. &
        any(root_temperature <= 0.0_dp) .or. &
        any(level_one_temperature <= 0.0_dp) .or. &
        any(level_two_temperature <= 0.0_dp)) return

    allocate(root_candidate, mold=root_state)
    allocate(root_candidate_temperature, mold=root_temperature)
    allocate(root_x_flux(nvar, 0:root_geometry%nx, root_geometry%ny))
    allocate(root_y_flux(nvar, root_geometry%nx, 0:root_geometry%ny))
    call advance_reactive_eb_level_2d( &
      species, root_state, root_temperature, root_geometry, solver, &
      reconstruction, limiter, selected_target, state_redist_max_order, dt, &
      root_candidate, root_candidate_temperature, root_x_flux, root_y_flux, &
      local_ok)
    if (.not. local_ok) return
    call initialize_amr_eb_flux_register_2d( &
      root_geometry, level_one_geometry, root_patch, nvar, root_register, &
      local_ok)
    if (.not. local_ok) return
    call accumulate_coarse_eb_fluxes_2d( &
      root_register, root_geometry, level_one_geometry, root_patch, &
      root_x_flux, root_y_flux, dt, local_ok)
    if (.not. local_ok) return

    allocate(level_one_candidate, source=level_one_state)
    allocate(level_one_candidate_temperature, source=level_one_temperature)
    allocate(level_one_start, mold=level_one_state)
    allocate(level_one_start_temperature, mold=level_one_temperature)
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
      call composite_eb_integral_2d( &
        level_one_candidate, level_one_geometry, level_two_candidate, &
        level_two_geometry, level_one_patch, level_one_integral_before, &
        local_ok)
      if (.not. local_ok) return
      level_one_start = level_one_candidate
      level_one_start_temperature = level_one_candidate_temperature
      alpha = substep_time_alpha( &
        reconstruction, level_one_substep, level_one_ratio)
      call build_reactive_eb_patch_exterior_2d( &
        species, root_state, root_temperature, root_candidate, &
        root_candidate_temperature, root_geometry, level_one_geometry, &
        root_patch, alpha, level_one_exterior, local_ok, level_one_candidate, &
        level_one_candidate_temperature)
      if (.not. local_ok) return
      call advance_reactive_eb_level_2d( &
        species, level_one_candidate, level_one_candidate_temperature, &
        level_one_geometry, solver, reconstruction, limiter, selected_target, &
        state_redist_max_order, level_one_dt, level_one_uncorrected, &
        level_one_uncorrected_temperature, level_one_x_flux, &
        level_one_y_flux, local_ok, level_one_exterior)
      if (.not. local_ok) return
      call accumulate_fine_eb_fluxes_2d( &
        root_register, root_geometry, level_one_geometry, root_patch, &
        level_one_x_flux, level_one_y_flux, level_one_dt, local_ok)
      if (.not. local_ok) return

      call initialize_amr_eb_flux_register_2d( &
        level_one_geometry, level_two_geometry, level_one_patch, nvar, &
        level_one_register, local_ok)
      if (.not. local_ok) return
      call accumulate_coarse_eb_fluxes_2d( &
        level_one_register, level_one_geometry, level_two_geometry, &
        level_one_patch, level_one_x_flux, level_one_y_flux, level_one_dt, &
        local_ok)
      if (.not. local_ok) return

      do level_two_substep = 1, level_two_ratio
        alpha = substep_time_alpha( &
          reconstruction, level_two_substep, level_two_ratio)
        call build_reactive_eb_patch_exterior_2d( &
          species, level_one_start, level_one_start_temperature, &
          level_one_uncorrected, level_one_uncorrected_temperature, &
          level_one_geometry, level_two_geometry, level_one_patch, alpha, &
          level_two_exterior, local_ok, level_two_candidate, &
          level_two_candidate_temperature)
        if (.not. local_ok) return
        call advance_reactive_eb_level_2d( &
          species, level_two_candidate, level_two_candidate_temperature, &
          level_two_geometry, solver, reconstruction, limiter, &
          selected_target, state_redist_max_order, level_two_dt, &
          level_two_work, level_two_work_temperature, level_two_x_flux, &
          level_two_y_flux, local_ok, level_two_exterior)
        if (.not. local_ok) return
        level_two_candidate = level_two_work
        level_two_candidate_temperature = level_two_work_temperature
        call accumulate_fine_eb_fluxes_2d( &
          level_one_register, level_one_geometry, level_two_geometry, &
          level_one_patch, level_two_x_flux, level_two_y_flux, level_two_dt, &
          local_ok)
        if (.not. local_ok) return
      end do

      call reflux_reactive_eb_state_patch_2d( &
        species, level_one_uncorrected, &
        level_one_uncorrected_temperature, level_one_geometry, &
        level_two_candidate, level_two_candidate_temperature, &
        level_two_geometry, level_one_patch, level_one_register, &
        level_one_refluxed, level_one_refluxed_temperature, &
        level_two_refluxed, level_two_refluxed_temperature, local_ok)
      if (.not. local_ok) return
      call average_down_reactive_eb_state_patch_2d( &
        species, level_one_refluxed, level_one_refluxed_temperature, &
        level_one_geometry, level_two_refluxed, level_two_geometry, &
        level_one_patch, level_one_candidate, &
        level_one_candidate_temperature, local_ok)
      if (.not. local_ok) return
      level_two_candidate = level_two_refluxed
      level_two_candidate_temperature = level_two_refluxed_temperature
      if (.not. level_two_interface_is_regular(level_two_geometry)) then
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

    allocate(root_refluxed, mold=root_state)
    allocate(root_refluxed_temperature, mold=root_temperature)
    call reflux_reactive_eb_state_patch_2d( &
      species, root_candidate, root_candidate_temperature, root_geometry, &
      level_one_candidate, level_one_candidate_temperature, &
      level_one_geometry, root_patch, root_register, root_refluxed, &
      root_refluxed_temperature, level_one_refluxed, &
      level_one_refluxed_temperature, local_ok)
    if (.not. local_ok) return
    call average_down_three_level_reactive_eb_state_2d( &
      species, root_refluxed, root_refluxed_temperature, root_geometry, &
      level_one_refluxed, level_one_refluxed_temperature, &
      level_one_geometry, root_patch, level_two_candidate, &
      level_two_candidate_temperature, level_two_geometry, level_one_patch, &
      root_candidate, root_candidate_temperature, level_one_candidate, &
      level_one_candidate_temperature, local_ok)
    if (.not. local_ok) return

    new_root_state = root_candidate
    new_root_temperature = root_candidate_temperature
    new_level_one_state = level_one_candidate
    new_level_one_temperature = level_one_candidate_temperature
    new_level_two_state = level_two_candidate
    new_level_two_temperature = level_two_candidate_temperature
    ok = .true.
  end subroutine advance_three_level_reactive_eb_hydro_2d

  pure real(dp) function substep_time_alpha( &
      reconstruction, substep, ratio) result(alpha)
    character(len=*), intent(in) :: reconstruction
    integer, intent(in) :: substep, ratio

    if (trim(reconstruction) == "characteristic_plm") then
      alpha = (real(substep, dp) - 0.5_dp) / real(ratio, dp)
    else
      alpha = real(substep - 1, dp) / real(ratio, dp)
    end if
  end function substep_time_alpha

  pure logical function level_two_patch_is_separated( &
      patch, level_one_geometry) result(separated)
    type(amr_eb_patch_2d), intent(in) :: patch
    type(eb_geometry_2d), intent(in) :: level_one_geometry

    separated = patch%coarse_i_lower >= 3 .and. &
      patch%coarse_i_upper <= level_one_geometry%nx - 2 .and. &
      patch%coarse_j_lower >= 3 .and. &
      patch%coarse_j_upper <= level_one_geometry%ny - 2
  end function level_two_patch_is_separated

  pure logical function level_two_interface_is_regular( &
      geometry) result(regular)
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), parameter :: tolerance = 64.0_dp * epsilon(1.0_dp)

    regular = .false.
    if (geometry%nx < 1 .or. geometry%ny < 1 .or. &
        .not. allocated(geometry%x_face_fraction) .or. &
        .not. allocated(geometry%y_face_fraction)) return
    regular = &
      all(abs(geometry%x_face_fraction(0, :) - 1.0_dp) <= tolerance) .and. &
      all(abs(geometry%x_face_fraction(geometry%nx, :) - 1.0_dp) <= &
        tolerance) .and. &
      all(abs(geometry%y_face_fraction(:, 0) - 1.0_dp) <= tolerance) .and. &
      all(abs(geometry%y_face_fraction(:, geometry%ny) - 1.0_dp) <= &
        tolerance)
  end function level_two_interface_is_regular

  subroutine close_cut_interface_conservation_2d( &
      species, integral_before, parent_state, parent_temperature, &
      parent_geometry, child_state, child_geometry, patch, x_flux, y_flux, &
      dt, closed_state, closed_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: integral_before(:)
    real(dp), intent(in) :: parent_state(:, :, :)
    real(dp), intent(in) :: parent_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: parent_geometry, child_geometry
    real(dp), intent(in) :: child_state(:, :, :)
    type(amr_eb_patch_2d), intent(in) :: patch
    real(dp), intent(in) :: x_flux(:, 0:, :), y_flux(:, :, 0:), dt
    real(dp), intent(out) :: closed_state(:, :, :)
    real(dp), intent(out) :: closed_temperature(:, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: current_integral(:), boundary_change(:)
    real(dp), allocatable :: residual(:), correction(:), primitive(:)
    real(dp) :: recipient_volume, recovered_temperature, sound_speed
    real(dp) :: scale, closure_tolerance, species_residual
    logical :: local_ok
    integer :: i, j, k, nvar, component

    closed_state = parent_state
    closed_temperature = parent_temperature
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar < 1 .or. size(integral_before) /= nvar .or. &
        any(shape(parent_state) /= &
          [nvar, parent_geometry%nx, parent_geometry%ny]) .or. &
        any(shape(parent_temperature) /= &
          [parent_geometry%nx, parent_geometry%ny]) .or. &
        any(shape(child_state) /= &
          [nvar, child_geometry%nx, child_geometry%ny]) .or. &
        any(shape(x_flux) /= &
          [nvar, parent_geometry%nx + 1, parent_geometry%ny]) .or. &
        any(shape(y_flux) /= &
          [nvar, parent_geometry%nx, parent_geometry%ny + 1]) .or. &
        any(shape(closed_state) /= shape(parent_state)) .or. &
        any(shape(closed_temperature) /= shape(parent_temperature)) .or. &
        .not. ieee_is_finite(dt) .or. dt <= 0.0_dp .or. &
        any(.not. ieee_is_finite(integral_before)) .or. &
        any(.not. ieee_is_finite(x_flux)) .or. &
        any(.not. ieee_is_finite(y_flux))) return

    allocate(current_integral(nvar), boundary_change(nvar))
    allocate(residual(nvar), correction(nvar))
    call composite_eb_integral_2d( &
      parent_state, parent_geometry, child_state, child_geometry, patch, &
      current_integral, local_ok)
    if (.not. local_ok) return
    boundary_change = 0.0_dp
    do j = 1, parent_geometry%ny
      boundary_change = boundary_change + dt * parent_geometry%dy * &
        (parent_geometry%x_face_fraction(0, j) * x_flux(:, 0, j) - &
         parent_geometry%x_face_fraction(parent_geometry%nx, j) * &
           x_flux(:, parent_geometry%nx, j))
    end do
    do i = 1, parent_geometry%nx
      boundary_change = boundary_change + dt * parent_geometry%dx * &
        (parent_geometry%y_face_fraction(i, 0) * y_flux(:, i, 0) - &
         parent_geometry%y_face_fraction(i, parent_geometry%ny) * &
           y_flux(:, i, parent_geometry%ny))
    end do
    residual = integral_before + boundary_change - current_integral
    correction = 0.0_dp
    correction(irho) = residual(irho)
    correction(iet) = residual(iet)
    species_residual = 0.0_dp
    do k = 1, size(species)
      component = reactive_species_component(k)
      correction(component) = residual(component)
      species_residual = species_residual + residual(component)
    end do
    scale = max(1.0_dp, abs(residual(irho)), abs(species_residual))
    closure_tolerance = 4096.0_dp * epsilon(1.0_dp) * scale
    if (abs(residual(irho) - species_residual) > closure_tolerance) return
    component = reactive_species_component(size(species))
    correction(component) = correction(component) + &
      residual(irho) - species_residual

    recipient_volume = 0.0_dp
    do j = 1, parent_geometry%ny
      do i = 1, parent_geometry%nx
        if (cell_is_inside_patch(patch, i, j) .or. &
            parent_geometry%cell_type(i, j) == eb_covered_cell) cycle
        recipient_volume = recipient_volume + &
          parent_geometry%volume_fraction(i, j) * &
          parent_geometry%dx * parent_geometry%dy
      end do
    end do
    if (.not. ieee_is_finite(recipient_volume) .or. &
        recipient_volume <= tiny(1.0_dp)) return
    correction = correction / recipient_volume
    if (any(.not. ieee_is_finite(correction))) return

    allocate(primitive(reactive_nprim(size(species))))
    do j = 1, parent_geometry%ny
      do i = 1, parent_geometry%nx
        if (cell_is_inside_patch(patch, i, j) .or. &
            parent_geometry%cell_type(i, j) == eb_covered_cell) cycle
        closed_state(:, i, j) = closed_state(:, i, j) + correction
        call reactive_conserved_to_primitive( &
          species, closed_state(:, i, j), closed_temperature(i, j), &
          primitive, recovered_temperature, sound_speed, local_ok)
        if (.not. local_ok) then
          closed_state = parent_state
          closed_temperature = parent_temperature
          return
        end if
        closed_temperature(i, j) = recovered_temperature
      end do
    end do
    ok = .true.
  end subroutine close_cut_interface_conservation_2d

  pure logical function cell_is_inside_patch(patch, i, j) result(inside)
    type(amr_eb_patch_2d), intent(in) :: patch
    integer, intent(in) :: i, j

    inside = i >= patch%coarse_i_lower .and. &
      i <= patch%coarse_i_upper .and. &
      j >= patch%coarse_j_lower .and. j <= patch%coarse_j_upper
  end function cell_is_inside_patch

end module amr_eb_multilevel_reactive_2d_mod
