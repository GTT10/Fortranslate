module eb_reactive_transport_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use transport_database_mod, only: gas_transport_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_conserved_to_primitive
  use reactive_boundary_2d_mod, only: reactive_boundary_set_2d
  use reactive_transport_2d_mod, only: &
    reactive_transport_fluxes_2d_faces, reactive_transport_timestep_2d
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_covered_cell
  use eb_reactive_reconstruction_2d_mod, only: &
    interpolate_reactive_eb_face_centroid_fluxes_2d
  use eb_reactive_redistribution_2d_mod, only: &
    advance_reactive_eb_state_redistributed_2d
  implicit none
  private

  real(dp), parameter :: eb_species_safety = 0.90_dp

  public :: reactive_eb_transport_timestep_2d
  public :: reactive_eb_transport_rhs_2d
  public :: reactive_eb_transport_euler_update_2d
  public :: advance_reactive_eb_transport_2d

contains

  subroutine reactive_eb_transport_timestep_2d( &
      species, transport, state, temperature, geometry, transport_cfl, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, dt, maximum_diffusivity, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), intent(in) :: transport_cfl
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled
    real(dp), intent(out) :: dt, maximum_diffusivity
    logical, intent(out) :: ok

    dt = 0.0_dp
    maximum_diffusivity = 0.0_dp
    ok = .false.
    if (.not. geometry%is_valid()) return
    if (size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        any(shape(temperature) /= [geometry%nx, geometry%ny])) return
    call reactive_transport_timestep_2d( &
      species, transport, state, temperature, geometry%nx, geometry%ny, &
      geometry%dx, geometry%dy, transport_cfl, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, dt, &
      maximum_diffusivity, ok)
  end subroutine reactive_eb_transport_timestep_2d

  subroutine limit_reactive_eb_transport_fluxes_2d( &
      species, state, geometry, dt, flux_x, flux_y, minimum_theta, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :)
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), intent(in) :: dt
    real(dp), intent(inout) :: flux_x(:, 0:, :), flux_y(:, :, 0:)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok

    real(dp), allocatable :: theta_cell(:, :)
    real(dp) :: outgoing, available, candidate, theta_face
    integer :: i, j, face_i, face_j, species_index, component

    minimum_theta = 1.0_dp
    ok = .false.
    if (.not. geometry%is_valid() .or. .not. ieee_is_finite(dt) .or. &
        dt < 0.0_dp .or. size(state, 1) /= reactive_nvar(size(species)) .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        size(flux_x, 1) /= size(state, 1) .or. &
        size(flux_x, 2) /= geometry%nx + 1 .or. &
        size(flux_x, 3) /= geometry%ny .or. &
        size(flux_y, 1) /= size(state, 1) .or. &
        size(flux_y, 2) /= geometry%nx .or. &
        size(flux_y, 3) /= geometry%ny + 1) return
    if (dt <= tiny(1.0_dp)) then
      ok = .true.
      return
    end if

    allocate(theta_cell(geometry%nx, geometry%ny))
    theta_cell = 1.0_dp
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) == eb_covered_cell) cycle
        do species_index = 1, size(species)
          component = reactive_species_component(species_index)
          outgoing = geometry%dy * ( &
            max(geometry%x_face_fraction(i, j) * &
              flux_x(component, i, j), 0.0_dp) + &
            max(-geometry%x_face_fraction(i - 1, j) * &
              flux_x(component, i - 1, j), 0.0_dp)) + &
            geometry%dx * ( &
            max(geometry%y_face_fraction(i, j) * &
              flux_y(component, i, j), 0.0_dp) + &
            max(-geometry%y_face_fraction(i, j - 1) * &
              flux_y(component, i, j - 1), 0.0_dp))
          available = max(0.0_dp, state(component, i, j)) * &
            geometry%volume_fraction(i, j) * geometry%dx * geometry%dy
          if (outgoing > 0.0_dp) then
            candidate = eb_species_safety * available / (dt * outgoing)
            theta_cell(i, j) = min(theta_cell(i, j), &
              max(0.0_dp, min(1.0_dp, candidate)))
          end if
        end do
      end do
    end do

    do j = 1, geometry%ny
      do face_i = 0, geometry%nx
        theta_face = 1.0_dp
        if (face_i >= 1) then
          if (geometry%cell_type(face_i, j) /= eb_covered_cell) &
            theta_face = min(theta_face, theta_cell(face_i, j))
        end if
        if (face_i + 1 <= geometry%nx) then
          if (geometry%cell_type(face_i + 1, j) /= eb_covered_cell) &
            theta_face = min(theta_face, theta_cell(face_i + 1, j))
        end if
        flux_x(:, face_i, j) = theta_face * flux_x(:, face_i, j)
        minimum_theta = min(minimum_theta, theta_face)
      end do
    end do
    do face_j = 0, geometry%ny
      do i = 1, geometry%nx
        theta_face = 1.0_dp
        if (face_j >= 1) then
          if (geometry%cell_type(i, face_j) /= eb_covered_cell) &
            theta_face = min(theta_face, theta_cell(i, face_j))
        end if
        if (face_j + 1 <= geometry%ny) then
          if (geometry%cell_type(i, face_j + 1) /= eb_covered_cell) &
            theta_face = min(theta_face, theta_cell(i, face_j + 1))
        end if
        flux_y(:, i, face_j) = theta_face * flux_y(:, i, face_j)
        minimum_theta = min(minimum_theta, theta_face)
      end do
    end do
    ok = all(ieee_is_finite(flux_x)) .and. &
      all(ieee_is_finite(flux_y)) .and. &
      ieee_is_finite(minimum_theta) .and. minimum_theta >= 0.0_dp .and. &
      minimum_theta <= 1.0_dp
  end subroutine limit_reactive_eb_transport_fluxes_2d

  subroutine reactive_eb_transport_rhs_2d( &
      species, transport, state, temperature, geometry, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, rhs, &
      minimum_theta, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), intent(in) :: dt
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    real(dp), intent(out) :: rhs(:, :, :)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok

    real(dp), allocatable :: center_x(:, :, :), center_y(:, :, :)
    real(dp), allocatable :: centroid_x(:, :, :), centroid_y(:, :, :)
    real(dp) :: regular_theta, fluid_volume
    logical :: local_ok
    integer :: i, j, nvar

    rhs = 0.0_dp
    minimum_theta = 1.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid() .or. &
        any(shape(rhs) /= shape(state))) return
    allocate(center_x(nvar, 0:geometry%nx, geometry%ny))
    allocate(center_y(nvar, geometry%nx, 0:geometry%ny))
    allocate(centroid_x(nvar, 0:geometry%nx, geometry%ny))
    allocate(centroid_y(nvar, geometry%nx, 0:geometry%ny))
    call reactive_transport_fluxes_2d_faces( &
      species, transport, state, temperature, geometry%nx, geometry%ny, &
      geometry%dx, geometry%dy, dt, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, center_x, center_y, regular_theta, local_ok, &
      boundaries)
    if (.not. local_ok) return
    call interpolate_reactive_eb_face_centroid_fluxes_2d( &
      geometry, center_x, center_y, centroid_x, centroid_y, local_ok)
    if (.not. local_ok) return
    call limit_reactive_eb_transport_fluxes_2d( &
      species, state, geometry, dt, centroid_x, centroid_y, minimum_theta, &
      local_ok)
    if (.not. local_ok) return
    minimum_theta = min(minimum_theta, regular_theta)

    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) == eb_covered_cell) cycle
        fluid_volume = geometry%volume_fraction(i, j) * &
          geometry%dx * geometry%dy
        if (fluid_volume <= 0.0_dp) return
        rhs(:, i, j) = -( &
          geometry%dy * ( &
            geometry%x_face_fraction(i, j) * centroid_x(:, i, j) - &
            geometry%x_face_fraction(i - 1, j) * &
              centroid_x(:, i - 1, j)) + &
          geometry%dx * ( &
            geometry%y_face_fraction(i, j) * centroid_y(:, i, j) - &
            geometry%y_face_fraction(i, j - 1) * &
              centroid_y(:, i, j - 1))) / fluid_volume
      end do
    end do
    ok = all(ieee_is_finite(rhs))
  end subroutine reactive_eb_transport_rhs_2d

  subroutine reactive_eb_transport_euler_update_2d( &
      species, transport, state, temperature, geometry, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, max_order, new_state, new_temperature, &
      minimum_theta, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), intent(in) :: dt, target_volume_fraction
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    integer, intent(in) :: max_order
    real(dp), intent(out) :: new_state(:, :, :), new_temperature(:, :)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok

    real(dp), allocatable :: rhs(:, :, :)
    logical :: local_ok

    new_state = state
    new_temperature = temperature
    minimum_theta = 1.0_dp
    ok = .false.
    if (any(shape(new_state) /= shape(state)) .or. &
        any(shape(new_temperature) /= shape(temperature))) return
    allocate(rhs, mold=state)
    call reactive_eb_transport_rhs_2d( &
      species, transport, state, temperature, geometry, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, rhs, &
      minimum_theta, local_ok)
    if (.not. local_ok) return
    call advance_reactive_eb_state_redistributed_2d( &
      species, state, temperature, geometry, rhs, dt, new_state, &
      new_temperature, local_ok, target_volume_fraction, max_order)
    if (.not. local_ok) return
    ok = .true.
  end subroutine reactive_eb_transport_euler_update_2d

  subroutine advance_reactive_eb_transport_2d( &
      species, transport, state, temperature, geometry, interval, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, max_order, minimum_theta, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(inout) :: state(:, :, :), temperature(:, :)
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), intent(in) :: interval, target_volume_fraction
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    integer, intent(in) :: max_order
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok

    real(dp), allocatable :: initial_state(:, :, :)
    real(dp), allocatable :: initial_temperature(:, :)
    real(dp), allocatable :: stage_state(:, :, :), stage_temperature(:, :)
    real(dp), allocatable :: euler_state(:, :, :), euler_temperature(:, :)
    real(dp), allocatable :: candidate_state(:, :, :)
    real(dp), allocatable :: candidate_temperature(:, :), primitive(:)
    real(dp) :: theta_one, theta_two, recovered_temperature, sound_speed
    logical :: local_ok
    integer :: i, j

    minimum_theta = 1.0_dp
    ok = .false.
    if (.not. ieee_is_finite(interval) .or. interval < 0.0_dp) return
    if (interval <= tiny(1.0_dp) .or. .not. (viscosity_enabled .or. &
        thermal_conduction_enabled .or. species_diffusion_enabled)) then
      ok = .true.
      return
    end if
    allocate(initial_state, source=state)
    allocate(initial_temperature, source=temperature)
    allocate(stage_state, mold=state)
    allocate(stage_temperature, mold=temperature)
    allocate(euler_state, mold=state)
    allocate(euler_temperature, mold=temperature)
    call reactive_eb_transport_euler_update_2d( &
      species, transport, initial_state, initial_temperature, geometry, &
      interval, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, max_order, stage_state, stage_temperature, &
      theta_one, local_ok)
    if (.not. local_ok) return
    call reactive_eb_transport_euler_update_2d( &
      species, transport, stage_state, stage_temperature, geometry, interval, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, max_order, euler_state, euler_temperature, &
      theta_two, local_ok)
    if (.not. local_ok) return

    allocate(candidate_state, source=0.5_dp * (initial_state + euler_state))
    allocate(candidate_temperature, source=initial_temperature)
    allocate(primitive(reactive_nprim(size(species))))
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) == eb_covered_cell) then
          candidate_state(:, i, j) = initial_state(:, i, j)
          cycle
        end if
        call reactive_conserved_to_primitive( &
          species, candidate_state(:, i, j), &
          0.5_dp * (initial_temperature(i, j) + &
            euler_temperature(i, j)), primitive, recovered_temperature, &
          sound_speed, local_ok)
        if (.not. local_ok) return
        candidate_temperature(i, j) = recovered_temperature
      end do
    end do
    state = candidate_state
    temperature = candidate_temperature
    minimum_theta = min(theta_one, theta_two)
    ok = .true.
  end subroutine advance_reactive_eb_transport_2d

end module eb_reactive_transport_2d_mod
