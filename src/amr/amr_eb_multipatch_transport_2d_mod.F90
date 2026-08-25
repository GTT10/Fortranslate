module amr_eb_multipatch_transport_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use state_indices_mod, only: irho
  use nasa7_thermo_mod, only: nasa7_species
  use transport_database_mod, only: gas_transport_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_conserved_to_primitive
  use reactive_boundary_2d_mod, only: reactive_boundary_set_2d
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_covered_cell
  use eb_reactive_reconstruction_2d_mod, only: &
    reactive_eb_exterior_state_2d
  use eb_reactive_redistribution_2d_mod, only: &
    advance_reactive_eb_state_redistributed_2d
  use eb_reactive_transport_2d_mod, only: &
    reactive_eb_transport_fluxes_rhs_2d
  use amr_eb_flux_register_2d_mod, only: &
    amr_eb_flux_register_2d, initialize_amr_eb_flux_register_2d, &
    accumulate_coarse_eb_fluxes_2d, accumulate_fine_eb_fluxes_2d, &
    reflux_reactive_eb_state_patch_2d
  use amr_eb_reactive_2d_mod, only: build_reactive_eb_patch_exterior_2d
  use amr_eb_regrid_2d_mod, only: &
    reactive_eb_patch_set_2d, average_down_reactive_eb_patch_set_2d, &
    composite_reactive_eb_patch_set_integral_2d
  use amr_eb_transport_2d_mod, only: recover_transport_temperature_2d
  use amr_eb_multilevel_reactive_2d_mod, only: &
    level_two_interface_is_regular
  implicit none
  private

  public :: advance_reactive_eb_patch_set_transport_euler_2d
  public :: advance_reactive_eb_patch_set_transport_2d
  public :: close_cut_patch_set_conservation_2d

contains

  subroutine advance_reactive_eb_patch_set_transport_euler_2d( &
      species, transport, coarse_state, coarse_temperature, coarse_geometry, &
      patch_set, dt, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, max_order, new_coarse_state, &
      new_coarse_temperature, new_patch_set, minimum_theta, ok, &
      failure_context)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: coarse_state(:, :, :), coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set
    real(dp), intent(in) :: dt, target_volume_fraction
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    integer, intent(in) :: max_order
    real(dp), intent(out) :: new_coarse_state(:, :, :)
    real(dp), intent(out) :: new_coarse_temperature(:, :)
    type(reactive_eb_patch_set_2d), intent(out) :: new_patch_set
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    character(len=*), intent(out), optional :: failure_context

    type(amr_eb_flux_register_2d) :: flux_register
    type(reactive_eb_exterior_state_2d) :: exterior
    type(reactive_eb_patch_set_2d) :: candidate_set
    real(dp), allocatable :: coarse_candidate(:, :, :)
    real(dp), allocatable :: coarse_candidate_temperature(:, :)
    real(dp), allocatable :: coarse_corrected(:, :, :)
    real(dp), allocatable :: coarse_corrected_temperature(:, :)
    real(dp), allocatable :: coarse_work(:, :, :)
    real(dp), allocatable :: coarse_work_temperature(:, :)
    real(dp), allocatable :: coarse_closed(:, :, :)
    real(dp), allocatable :: coarse_closed_temperature(:, :)
    real(dp), allocatable :: coarse_rhs(:, :, :)
    real(dp), allocatable :: coarse_x_flux(:, :, :)
    real(dp), allocatable :: coarse_y_flux(:, :, :)
    real(dp), allocatable :: fine_rhs(:, :, :)
    real(dp), allocatable :: fine_work(:, :, :)
    real(dp), allocatable :: fine_work_temperature(:, :)
    real(dp), allocatable :: fine_x_flux(:, :, :)
    real(dp), allocatable :: fine_y_flux(:, :, :)
    real(dp), allocatable :: integral_before(:)
    real(dp) :: coarse_theta, fine_theta, fine_dt, alpha
    logical :: local_ok, cut_interface
    integer :: child, nvar, ratio, substep

    new_coarse_state = coarse_state
    new_coarse_temperature = coarse_temperature
    new_patch_set = patch_set
    minimum_theta = 1.0_dp
    ok = .false.
    if (present(failure_context)) failure_context = "input validation"
    nvar = reactive_nvar(size(species))
    if (nvar < 1 .or. size(transport) /= size(species) .or. &
        .not. ieee_is_finite(dt) .or. dt <= 0.0_dp .or. &
        .not. ieee_is_finite(target_volume_fraction) .or. &
        target_volume_fraction <= 0.0_dp .or. &
        target_volume_fraction > 1.0_dp .or. &
        any(shape(coarse_state) /= &
          [nvar, coarse_geometry%nx, coarse_geometry%ny]) .or. &
        any(shape(coarse_temperature) /= &
          [coarse_geometry%nx, coarse_geometry%ny]) .or. &
        any(shape(new_coarse_state) /= shape(coarse_state)) .or. &
        any(shape(new_coarse_temperature) /= shape(coarse_temperature)) .or. &
        .not. patch_set%is_valid(coarse_geometry, nvar)) return

    allocate(integral_before(nvar))
    call composite_reactive_eb_patch_set_integral_2d( &
      coarse_state, coarse_geometry, patch_set, integral_before, local_ok)
    if (.not. local_ok) return
    allocate(coarse_candidate, mold=coarse_state)
    allocate(coarse_candidate_temperature, mold=coarse_temperature)
    allocate(coarse_rhs, mold=coarse_state)
    allocate(coarse_x_flux(nvar, 0:coarse_geometry%nx, coarse_geometry%ny))
    allocate(coarse_y_flux(nvar, coarse_geometry%nx, 0:coarse_geometry%ny))
    if (present(failure_context)) failure_context = "coarse transport"
    call reactive_eb_transport_fluxes_rhs_2d( &
      species, transport, coarse_state, coarse_temperature, coarse_geometry, &
      dt, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      coarse_rhs, coarse_x_flux, coarse_y_flux, coarse_theta, local_ok)
    if (.not. local_ok) return
    call advance_reactive_eb_state_redistributed_2d( &
      species, coarse_state, coarse_temperature, coarse_geometry, coarse_rhs, &
      dt, coarse_candidate, coarse_candidate_temperature, local_ok, &
      target_volume_fraction, max_order)
    if (.not. local_ok) return
    minimum_theta = min(minimum_theta, coarse_theta)

    allocate(coarse_corrected, source=coarse_candidate)
    allocate(coarse_corrected_temperature, &
      source=coarse_candidate_temperature)
    allocate(coarse_work, mold=coarse_state)
    allocate(coarse_work_temperature, mold=coarse_temperature)
    candidate_set = patch_set
    cut_interface = .false.
    do child = 1, candidate_set%patch_count()
      cut_interface = cut_interface .or. .not. level_two_interface_is_regular( &
        candidate_set%children(child)%geometry)
      if (present(failure_context)) &
        failure_context = "child flux-register initialization"
      call initialize_amr_eb_flux_register_2d( &
        coarse_geometry, candidate_set%children(child)%geometry, &
        candidate_set%children(child)%patch, nvar, flux_register, local_ok)
      if (.not. local_ok) return
      call accumulate_coarse_eb_fluxes_2d( &
        flux_register, coarse_geometry, &
        candidate_set%children(child)%geometry, &
        candidate_set%children(child)%patch, coarse_x_flux, coarse_y_flux, &
        dt, local_ok)
      if (.not. local_ok) return

      if (allocated(fine_rhs)) deallocate(fine_rhs)
      if (allocated(fine_work)) deallocate(fine_work)
      if (allocated(fine_work_temperature)) deallocate(fine_work_temperature)
      if (allocated(fine_x_flux)) deallocate(fine_x_flux)
      if (allocated(fine_y_flux)) deallocate(fine_y_flux)
      allocate(fine_rhs, mold=candidate_set%children(child)%state)
      allocate(fine_work, mold=candidate_set%children(child)%state)
      allocate(fine_work_temperature, &
        mold=candidate_set%children(child)%temperature)
      allocate(fine_x_flux( &
        nvar, 0:candidate_set%children(child)%geometry%nx, &
        candidate_set%children(child)%geometry%ny))
      allocate(fine_y_flux( &
        nvar, candidate_set%children(child)%geometry%nx, &
        0:candidate_set%children(child)%geometry%ny))
      ratio = candidate_set%children(child)%patch%refinement_ratio
      fine_dt = dt / real(ratio, dp)
      do substep = 1, ratio
        alpha = real(substep - 1, dp) / real(ratio, dp)
        if (present(failure_context)) failure_context = "child exterior fill"
        call build_reactive_eb_patch_exterior_2d( &
          species, coarse_state, coarse_temperature, coarse_candidate, &
          coarse_candidate_temperature, coarse_geometry, &
          candidate_set%children(child)%geometry, &
          candidate_set%children(child)%patch, alpha, exterior, local_ok, &
          candidate_set%children(child)%state, &
          candidate_set%children(child)%temperature)
        if (.not. local_ok) return
        if (present(failure_context)) write(failure_context, '(a,i0,a,i0)') &
          "fine transport child ", child, " substep ", substep
        call reactive_eb_transport_fluxes_rhs_2d( &
          species, transport, candidate_set%children(child)%state, &
          candidate_set%children(child)%temperature, &
          candidate_set%children(child)%geometry, fine_dt, &
          viscosity_enabled, thermal_conduction_enabled, &
          species_diffusion_enabled, barodiffusion_enabled, boundaries, &
          fine_rhs, fine_x_flux, fine_y_flux, fine_theta, local_ok, exterior)
        if (.not. local_ok) return
        minimum_theta = min(minimum_theta, fine_theta)
        call advance_reactive_eb_state_redistributed_2d( &
          species, candidate_set%children(child)%state, &
          candidate_set%children(child)%temperature, &
          candidate_set%children(child)%geometry, fine_rhs, fine_dt, &
          fine_work, fine_work_temperature, local_ok, target_volume_fraction, &
          max_order)
        if (.not. local_ok) return
        candidate_set%children(child)%state = fine_work
        candidate_set%children(child)%temperature = fine_work_temperature
        call accumulate_fine_eb_fluxes_2d( &
          flux_register, coarse_geometry, &
          candidate_set%children(child)%geometry, &
          candidate_set%children(child)%patch, fine_x_flux, fine_y_flux, &
          fine_dt, local_ok)
        if (.not. local_ok) return
      end do

      if (present(failure_context)) failure_context = "child reflux"
      call reflux_reactive_eb_state_patch_2d( &
        species, coarse_corrected, coarse_corrected_temperature, &
        coarse_geometry, candidate_set%children(child)%state, &
        candidate_set%children(child)%temperature, &
        candidate_set%children(child)%geometry, &
        candidate_set%children(child)%patch, flux_register, coarse_work, &
        coarse_work_temperature, fine_work, fine_work_temperature, local_ok)
      if (.not. local_ok) return
      coarse_corrected = coarse_work
      coarse_corrected_temperature = coarse_work_temperature
      candidate_set%children(child)%state = fine_work
      candidate_set%children(child)%temperature = fine_work_temperature
    end do

    if (present(failure_context)) failure_context = "patch-set average down"
    call average_down_reactive_eb_patch_set_2d( &
      species, coarse_corrected, coarse_corrected_temperature, &
      coarse_geometry, candidate_set, coarse_work, &
      coarse_work_temperature, local_ok)
    if (.not. local_ok) return
    if (cut_interface) then
      allocate(coarse_closed, mold=coarse_state)
      allocate(coarse_closed_temperature, mold=coarse_temperature)
      if (present(failure_context)) &
        failure_context = "patch-set cut-interface conservation"
      call close_cut_patch_set_conservation_2d( &
        species, integral_before, coarse_work, coarse_work_temperature, &
        coarse_geometry, candidate_set, coarse_x_flux, coarse_y_flux, dt, &
        coarse_closed, coarse_closed_temperature, local_ok)
      if (.not. local_ok) return
      coarse_work = coarse_closed
      coarse_work_temperature = coarse_closed_temperature
    end if
    if (.not. candidate_set%is_valid(coarse_geometry, nvar)) return
    new_coarse_state = coarse_work
    new_coarse_temperature = coarse_work_temperature
    new_patch_set = candidate_set
    ok = .true.
    if (present(failure_context)) failure_context = "none"
  end subroutine advance_reactive_eb_patch_set_transport_euler_2d

  subroutine advance_reactive_eb_patch_set_transport_2d( &
      species, transport, coarse_state, coarse_temperature, coarse_geometry, &
      patch_set, interval, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, max_order, new_coarse_state, &
      new_coarse_temperature, new_patch_set, minimum_theta, ok, &
      failure_context)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: coarse_state(:, :, :), coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set
    real(dp), intent(in) :: interval, target_volume_fraction
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    integer, intent(in) :: max_order
    real(dp), intent(out) :: new_coarse_state(:, :, :)
    real(dp), intent(out) :: new_coarse_temperature(:, :)
    type(reactive_eb_patch_set_2d), intent(out) :: new_patch_set
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    character(len=*), intent(out), optional :: failure_context

    type(reactive_eb_patch_set_2d) :: stage_set, euler_set, candidate_set
    real(dp), allocatable :: stage_coarse(:, :, :)
    real(dp), allocatable :: stage_coarse_temperature(:, :)
    real(dp), allocatable :: euler_coarse(:, :, :)
    real(dp), allocatable :: euler_coarse_temperature(:, :)
    real(dp), allocatable :: candidate_coarse(:, :, :)
    real(dp), allocatable :: candidate_coarse_temperature(:, :)
    real(dp), allocatable :: synchronized_coarse(:, :, :)
    real(dp), allocatable :: synchronized_coarse_temperature(:, :)
    real(dp) :: theta_one, theta_two
    logical :: local_ok
    integer :: child, nvar

    new_coarse_state = coarse_state
    new_coarse_temperature = coarse_temperature
    new_patch_set = patch_set
    minimum_theta = 1.0_dp
    ok = .false.
    if (present(failure_context)) failure_context = "input validation"
    if (.not. ieee_is_finite(interval) .or. interval < 0.0_dp) return
    if (interval <= tiny(1.0_dp) .or. .not. (viscosity_enabled .or. &
        thermal_conduction_enabled .or. species_diffusion_enabled)) then
      ok = .true.
      if (present(failure_context)) failure_context = "none"
      return
    end if

    allocate(stage_coarse, mold=coarse_state)
    allocate(stage_coarse_temperature, mold=coarse_temperature)
    call advance_reactive_eb_patch_set_transport_euler_2d( &
      species, transport, coarse_state, coarse_temperature, coarse_geometry, &
      patch_set, interval, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, max_order, stage_coarse, &
      stage_coarse_temperature, stage_set, theta_one, local_ok, &
      failure_context)
    if (.not. local_ok) return

    allocate(euler_coarse, mold=coarse_state)
    allocate(euler_coarse_temperature, mold=coarse_temperature)
    call advance_reactive_eb_patch_set_transport_euler_2d( &
      species, transport, stage_coarse, stage_coarse_temperature, &
      coarse_geometry, stage_set, interval, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, target_volume_fraction, max_order, &
      euler_coarse, euler_coarse_temperature, euler_set, theta_two, local_ok, &
      failure_context)
    if (.not. local_ok) return

    nvar = reactive_nvar(size(species))
    if (.not. euler_set%is_valid(coarse_geometry, nvar) .or. &
        euler_set%patch_count() /= patch_set%patch_count()) return
    allocate(candidate_coarse, source=0.5_dp * &
      (coarse_state + euler_coarse))
    allocate(candidate_coarse_temperature, mold=coarse_temperature)
    candidate_set = patch_set
    do child = 1, candidate_set%patch_count()
      candidate_set%children(child)%state = 0.5_dp * &
        (patch_set%children(child)%state + &
         euler_set%children(child)%state)
      call recover_transport_temperature_2d( &
        species, candidate_set%children(child)%state, &
        0.5_dp * (patch_set%children(child)%temperature + &
          euler_set%children(child)%temperature), &
        candidate_set%children(child)%geometry, &
        candidate_set%children(child)%temperature, local_ok)
      if (.not. local_ok) return
    end do
    call recover_transport_temperature_2d( &
      species, candidate_coarse, &
      0.5_dp * (coarse_temperature + euler_coarse_temperature), &
      coarse_geometry, candidate_coarse_temperature, local_ok)
    if (.not. local_ok) return
    allocate(synchronized_coarse, mold=coarse_state)
    allocate(synchronized_coarse_temperature, mold=coarse_temperature)
    call average_down_reactive_eb_patch_set_2d( &
      species, candidate_coarse, candidate_coarse_temperature, &
      coarse_geometry, candidate_set, synchronized_coarse, &
      synchronized_coarse_temperature, local_ok)
    if (.not. local_ok) return
    if (.not. candidate_set%is_valid(coarse_geometry, nvar)) return
    new_coarse_state = synchronized_coarse
    new_coarse_temperature = synchronized_coarse_temperature
    new_patch_set = candidate_set
    minimum_theta = min(theta_one, theta_two)
    ok = .true.
    if (present(failure_context)) failure_context = "none"
  end subroutine advance_reactive_eb_patch_set_transport_2d

  subroutine close_cut_patch_set_conservation_2d( &
      species, integral_before, coarse_state, coarse_temperature, &
      coarse_geometry, patch_set, x_flux, y_flux, dt, closed_state, &
      closed_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: integral_before(:)
    real(dp), intent(in) :: coarse_state(:, :, :), coarse_temperature(:, :)
    type(eb_geometry_2d), intent(in) :: coarse_geometry
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set
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

    closed_state = coarse_state
    closed_temperature = coarse_temperature
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar < 1 .or. size(integral_before) /= nvar .or. &
        any(shape(x_flux) /= &
          [nvar, coarse_geometry%nx + 1, coarse_geometry%ny]) .or. &
        any(shape(y_flux) /= &
          [nvar, coarse_geometry%nx, coarse_geometry%ny + 1]) .or. &
        .not. ieee_is_finite(dt) .or. dt <= 0.0_dp) return
    allocate(current_integral(nvar), boundary_change(nvar))
    allocate(residual(nvar), correction(nvar))
    call composite_reactive_eb_patch_set_integral_2d( &
      coarse_state, coarse_geometry, patch_set, current_integral, local_ok)
    if (.not. local_ok) return
    boundary_change = 0.0_dp
    do j = 1, coarse_geometry%ny
      boundary_change = boundary_change + dt * coarse_geometry%dy * &
        (coarse_geometry%x_face_fraction(0, j) * x_flux(:, 0, j) - &
         coarse_geometry%x_face_fraction(coarse_geometry%nx, j) * &
           x_flux(:, coarse_geometry%nx, j))
    end do
    do i = 1, coarse_geometry%nx
      boundary_change = boundary_change + dt * coarse_geometry%dx * &
        (coarse_geometry%y_face_fraction(i, 0) * y_flux(:, i, 0) - &
         coarse_geometry%y_face_fraction(i, coarse_geometry%ny) * &
           y_flux(:, i, coarse_geometry%ny))
    end do
    residual = integral_before + boundary_change - current_integral
    correction = residual
    species_residual = 0.0_dp
    do k = 1, size(species)
      component = reactive_species_component(k)
      species_residual = species_residual + residual(component)
    end do
    scale = max(1.0_dp, maxval(abs(residual)), abs(species_residual))
    closure_tolerance = 4096.0_dp * epsilon(1.0_dp) * scale
    if (abs(residual(irho) - species_residual) > closure_tolerance) return
    component = reactive_species_component(size(species))
    correction(component) = correction(component) + &
      residual(irho) - species_residual

    recipient_volume = 0.0_dp
    do j = 1, coarse_geometry%ny
      do i = 1, coarse_geometry%nx
        if (cell_is_inside_any_patch(patch_set, i, j) .or. &
            coarse_geometry%cell_type(i, j) == eb_covered_cell) cycle
        recipient_volume = recipient_volume + &
          coarse_geometry%volume_fraction(i, j) * &
          coarse_geometry%dx * coarse_geometry%dy
      end do
    end do
    if (.not. ieee_is_finite(recipient_volume) .or. &
        recipient_volume <= tiny(1.0_dp)) return
    correction = correction / recipient_volume
    if (any(.not. ieee_is_finite(correction))) return

    allocate(primitive(reactive_nprim(size(species))))
    do j = 1, coarse_geometry%ny
      do i = 1, coarse_geometry%nx
        if (cell_is_inside_any_patch(patch_set, i, j) .or. &
            coarse_geometry%cell_type(i, j) == eb_covered_cell) cycle
        closed_state(:, i, j) = closed_state(:, i, j) + correction
        call reactive_conserved_to_primitive( &
          species, closed_state(:, i, j), closed_temperature(i, j), &
          primitive, recovered_temperature, sound_speed, local_ok)
        if (.not. local_ok) then
          closed_state = coarse_state
          closed_temperature = coarse_temperature
          return
        end if
        closed_temperature(i, j) = recovered_temperature
      end do
    end do
    ok = .true.
  end subroutine close_cut_patch_set_conservation_2d

  pure logical function cell_is_inside_any_patch(patch_set, i, j) &
      result(inside)
    type(reactive_eb_patch_set_2d), intent(in) :: patch_set
    integer, intent(in) :: i, j
    integer :: child

    inside = .false.
    do child = 1, patch_set%patch_count()
      inside = &
        i >= patch_set%children(child)%patch%coarse_i_lower .and. &
        i <= patch_set%children(child)%patch%coarse_i_upper .and. &
        j >= patch_set%children(child)%patch%coarse_j_lower .and. &
        j <= patch_set%children(child)%patch%coarse_j_upper
      if (inside) return
    end do
  end function cell_is_inside_any_patch

end module amr_eb_multipatch_transport_2d_mod
