module eb_reactive_transport_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use transport_database_mod, only: gas_transport_species
  use mixture_transport_mod, only: mixture_transport_coefficients
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_mass_fraction_component, reactive_conserved_to_primitive
  use reactive_boundary_2d_mod, only: &
    reactive_boundary_face_2d, reactive_boundary_set_2d
  use reactive_transport_2d_mod, only: &
    reactive_transport_exterior_2d, reactive_transport_fluxes_2d_faces, &
    reactive_transport_timestep_2d
  use eb_geometry_2d_mod, only: &
    eb_geometry_2d, eb_covered_cell, eb_cut_cell
  use eb_reactive_reconstruction_2d_mod, only: &
    reactive_eb_exterior_state_2d, &
    interpolate_reactive_eb_face_centroid_fluxes_2d
  use eb_reactive_redistribution_2d_mod, only: &
    advance_reactive_eb_state_redistributed_2d
  implicit none
  private

  real(dp), parameter :: eb_species_safety = 0.90_dp

  public :: reactive_eb_transport_timestep_2d
  public :: reactive_eb_wall_transport_flux_2d
  public :: reactive_eb_transport_fluxes_rhs_2d
  public :: reactive_eb_transport_rhs_2d
  public :: reactive_eb_transport_euler_update_2d
  public :: advance_reactive_eb_transport_2d

contains

  subroutine reactive_eb_wall_transport_flux_2d( &
      species, transport, conserved, temperature_guess, fluid_normal, &
      normal_distance, wall, viscosity_enabled, &
      thermal_conduction_enabled, flux, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: conserved(:), temperature_guess
    real(dp), intent(in) :: fluid_normal(2), normal_distance
    type(reactive_boundary_face_2d), intent(in) :: wall
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    real(dp), intent(out) :: flux(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:), mass_fractions(:), diffusion(:)
    real(dp) :: recovered_temperature, sound_speed, normal_magnitude
    real(dp) :: viscosity, conductivity, normal_velocity_difference
    real(dp) :: velocity_difference(3), momentum_flux(3)
    logical :: local_ok
    integer :: k, nspecies

    flux = 0.0_dp
    ok = .false.
    nspecies = size(species)
    if (nspecies < 1 .or. size(transport) /= nspecies .or. &
        size(conserved) /= reactive_nvar(nspecies) .or. &
        size(flux) /= size(conserved) .or. &
        .not. ieee_is_finite(temperature_guess) .or. &
        .not. all(ieee_is_finite(fluid_normal)) .or. &
        .not. ieee_is_finite(normal_distance) .or. &
        normal_distance <= tiny(1.0_dp) .or. &
        .not. all(ieee_is_finite(wall%wall_velocity)) .or. &
        .not. ieee_is_finite(wall%wall_temperature) .or. &
        wall%wall_temperature <= 0.0_dp) return
    if (trim(wall%kind) /= "slip_wall" .and. &
        trim(wall%kind) /= "no_slip_wall") return
    if (trim(wall%thermal) /= "adiabatic" .and. &
        trim(wall%thermal) /= "isothermal") return
    if (trim(wall%wall_species) /= "impermeable") return
    normal_magnitude = sqrt(sum(fluid_normal**2))
    if (.not. ieee_is_finite(normal_magnitude) .or. &
        abs(normal_magnitude - 1.0_dp) > &
          128.0_dp * epsilon(1.0_dp)) return
    if ((.not. viscosity_enabled .or. trim(wall%kind) == "slip_wall") .and. &
        (.not. thermal_conduction_enabled .or. &
         trim(wall%thermal) == "adiabatic")) then
      ok = .true.
      return
    end if

    allocate(primitive(reactive_nprim(nspecies)))
    allocate(mass_fractions(nspecies), diffusion(nspecies))
    call reactive_conserved_to_primitive( &
      species, conserved, temperature_guess, primitive, &
      recovered_temperature, sound_speed, local_ok)
    if (.not. local_ok) return
    do k = 1, nspecies
      mass_fractions(k) = &
        primitive(reactive_mass_fraction_component(k))
    end do
    call mixture_transport_coefficients( &
      species, transport, mass_fractions, recovered_temperature, &
      primitive(5), viscosity, conductivity, diffusion, local_ok)
    if (.not. local_ok) return

    if (viscosity_enabled .and. trim(wall%kind) == "no_slip_wall") then
      velocity_difference = primitive(2:4) - wall%wall_velocity
      normal_velocity_difference = &
        velocity_difference(1) * fluid_normal(1) + &
        velocity_difference(2) * fluid_normal(2)
      momentum_flux(1) = -viscosity * ( &
        velocity_difference(1) + &
        normal_velocity_difference * fluid_normal(1) / 3.0_dp) / &
        normal_distance
      momentum_flux(2) = -viscosity * ( &
        velocity_difference(2) + &
        normal_velocity_difference * fluid_normal(2) / 3.0_dp) / &
        normal_distance
      momentum_flux(3) = &
        -viscosity * velocity_difference(3) / normal_distance
      flux(imx:imz) = momentum_flux
      flux(iet) = dot_product(momentum_flux, wall%wall_velocity)
    end if
    if (thermal_conduction_enabled .and. &
        trim(wall%thermal) == "isothermal") then
      flux(iet) = flux(iet) - conductivity * &
        (recovered_temperature - wall%wall_temperature) / normal_distance
    end if
    flux(irho) = 0.0_dp
    do k = 1, nspecies
      flux(reactive_species_component(k)) = 0.0_dp
    end do
    ok = all(ieee_is_finite(flux))
  end subroutine reactive_eb_wall_transport_flux_2d

  subroutine recover_transport_exterior_cell( &
      species, exterior_state, exterior_temperature, fallback_state, &
      fallback_temperature, primitive, temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: exterior_state(:), exterior_temperature
    real(dp), intent(in) :: fallback_state(:), fallback_temperature
    real(dp), intent(out) :: primitive(:), temperature
    logical, intent(out) :: ok

    real(dp) :: sound_speed

    call reactive_conserved_to_primitive( &
      species, exterior_state, exterior_temperature, primitive, &
      temperature, sound_speed, ok)
    if (ok) return
    call reactive_conserved_to_primitive( &
      species, fallback_state, fallback_temperature, primitive, &
      temperature, sound_speed, ok)
  end subroutine recover_transport_exterior_cell

  subroutine build_reactive_transport_exterior_2d( &
      species, state, temperature, geometry, eb_exterior, exterior, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    type(eb_geometry_2d), intent(in) :: geometry
    type(reactive_eb_exterior_state_2d), intent(in) :: eb_exterior
    type(reactive_transport_exterior_2d), intent(out) :: exterior
    logical, intent(out) :: ok

    type(reactive_transport_exterior_2d) :: candidate
    logical :: local_ok
    integer :: i, j, nprim, nvar

    ok = .false.
    nvar = reactive_nvar(size(species))
    nprim = reactive_nprim(size(species))
    if (nvar < 1 .or. .not. geometry%is_valid() .or. &
        any(shape(state) /= [nvar, geometry%nx, geometry%ny]) .or. &
        any(shape(temperature) /= [geometry%nx, geometry%ny]) .or. &
        .not. eb_exterior%is_valid(geometry, nvar)) return
    allocate(candidate%x_lower_primitive(nprim, geometry%ny))
    allocate(candidate%x_upper_primitive(nprim, geometry%ny))
    allocate(candidate%y_lower_primitive(nprim, geometry%nx))
    allocate(candidate%y_upper_primitive(nprim, geometry%nx))
    allocate(candidate%x_lower_temperature(geometry%ny))
    allocate(candidate%x_upper_temperature(geometry%ny))
    allocate(candidate%y_lower_temperature(geometry%nx))
    allocate(candidate%y_upper_temperature(geometry%nx))
    do j = 1, geometry%ny
      call recover_transport_exterior_cell( &
        species, eb_exterior%x_lower_state(:, j), &
        eb_exterior%x_lower_temperature(j), state(:, 1, j), &
        temperature(1, j), candidate%x_lower_primitive(:, j), &
        candidate%x_lower_temperature(j), local_ok)
      if (.not. local_ok) return
      call recover_transport_exterior_cell( &
        species, eb_exterior%x_upper_state(:, j), &
        eb_exterior%x_upper_temperature(j), state(:, geometry%nx, j), &
        temperature(geometry%nx, j), candidate%x_upper_primitive(:, j), &
        candidate%x_upper_temperature(j), local_ok)
      if (.not. local_ok) return
    end do
    do i = 1, geometry%nx
      call recover_transport_exterior_cell( &
        species, eb_exterior%y_lower_state(:, i), &
        eb_exterior%y_lower_temperature(i), state(:, i, 1), &
        temperature(i, 1), candidate%y_lower_primitive(:, i), &
        candidate%y_lower_temperature(i), local_ok)
      if (.not. local_ok) return
      call recover_transport_exterior_cell( &
        species, eb_exterior%y_upper_state(:, i), &
        eb_exterior%y_upper_temperature(i), state(:, i, geometry%ny), &
        temperature(i, geometry%ny), candidate%y_upper_primitive(:, i), &
        candidate%y_upper_temperature(i), local_ok)
      if (.not. local_ok) return
    end do
    if (.not. candidate%is_valid(nprim, geometry%nx, geometry%ny)) return
    exterior = candidate
    ok = .true.
  end subroutine build_reactive_transport_exterior_2d

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

  subroutine reactive_eb_transport_fluxes_rhs_2d( &
      species, transport, state, temperature, geometry, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, rhs, &
      x_flux, y_flux, minimum_theta, ok, exterior)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), intent(in) :: dt
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    real(dp), intent(out) :: rhs(:, :, :)
    real(dp), intent(out) :: x_flux(:, 0:, :), y_flux(:, :, 0:)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    type(reactive_eb_exterior_state_2d), intent(in), optional :: exterior

    real(dp), allocatable :: center_x(:, :, :), center_y(:, :, :)
    real(dp), allocatable :: wall_flux(:)
    type(reactive_transport_exterior_2d) :: transport_exterior
    real(dp) :: regular_theta, fluid_volume, normal_distance
    real(dp) :: cell_centroid_x, cell_centroid_y
    real(dp) :: fluid_normal(2)
    logical :: local_ok
    integer :: i, j, nvar

    rhs = 0.0_dp
    x_flux = 0.0_dp
    y_flux = 0.0_dp
    minimum_theta = 1.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid() .or. &
        any(shape(rhs) /= shape(state)) .or. &
        size(x_flux, 1) /= nvar .or. &
        size(x_flux, 2) /= geometry%nx + 1 .or. &
        size(x_flux, 3) /= geometry%ny .or. &
        size(y_flux, 1) /= nvar .or. &
        size(y_flux, 2) /= geometry%nx .or. &
        size(y_flux, 3) /= geometry%ny + 1) return
    allocate(center_x(nvar, 0:geometry%nx, geometry%ny))
    allocate(center_y(nvar, geometry%nx, 0:geometry%ny))
    allocate(wall_flux(nvar))
    if (present(exterior)) then
      call build_reactive_transport_exterior_2d( &
        species, state, temperature, geometry, exterior, &
        transport_exterior, local_ok)
      if (.not. local_ok) return
      call reactive_transport_fluxes_2d_faces( &
        species, transport, state, temperature, geometry%nx, geometry%ny, &
        geometry%dx, geometry%dy, dt, viscosity_enabled, &
        thermal_conduction_enabled, species_diffusion_enabled, &
        barodiffusion_enabled, center_x, center_y, regular_theta, local_ok, &
        boundaries, transport_exterior)
    else
      call reactive_transport_fluxes_2d_faces( &
        species, transport, state, temperature, geometry%nx, geometry%ny, &
        geometry%dx, geometry%dy, dt, viscosity_enabled, &
        thermal_conduction_enabled, species_diffusion_enabled, &
        barodiffusion_enabled, center_x, center_y, regular_theta, local_ok, &
        boundaries)
    end if
    if (.not. local_ok) return
    call interpolate_reactive_eb_face_centroid_fluxes_2d( &
      geometry, center_x, center_y, x_flux, y_flux, local_ok)
    if (.not. local_ok) return
    call limit_reactive_eb_transport_fluxes_2d( &
      species, state, geometry, dt, x_flux, y_flux, minimum_theta, &
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
            geometry%x_face_fraction(i, j) * x_flux(:, i, j) - &
            geometry%x_face_fraction(i - 1, j) * &
              x_flux(:, i - 1, j)) + &
          geometry%dx * ( &
            geometry%y_face_fraction(i, j) * y_flux(:, i, j) - &
            geometry%y_face_fraction(i, j - 1) * &
              y_flux(:, i, j - 1))) / fluid_volume
        if (geometry%cell_type(i, j) == eb_cut_cell .and. &
            (viscosity_enabled .or. thermal_conduction_enabled)) then
          cell_centroid_x = geometry%x_lower + &
            (real(i, dp) - 0.5_dp + geometry%cell_centroid_x(i, j)) * &
            geometry%dx
          cell_centroid_y = geometry%y_lower + &
            (real(j, dp) - 0.5_dp + geometry%cell_centroid_y(i, j)) * &
            geometry%dy
          fluid_normal = [geometry%boundary_normal_x(i, j), &
            geometry%boundary_normal_y(i, j)]
          normal_distance = &
            (cell_centroid_x - geometry%boundary_centroid_x(i, j)) * &
              fluid_normal(1) + &
            (cell_centroid_y - geometry%boundary_centroid_y(i, j)) * &
              fluid_normal(2)
          call reactive_eb_wall_transport_flux_2d( &
            species, transport, state(:, i, j), temperature(i, j), &
            fluid_normal, normal_distance, boundaries%embedded_wall, &
            viscosity_enabled, thermal_conduction_enabled, wall_flux, &
            local_ok)
          if (.not. local_ok) then
            rhs = 0.0_dp
            return
          end if
          rhs(:, i, j) = rhs(:, i, j) + &
            geometry%boundary_length(i, j) * wall_flux / fluid_volume
        end if
      end do
    end do
    ok = all(ieee_is_finite(rhs))
  end subroutine reactive_eb_transport_fluxes_rhs_2d

  subroutine reactive_eb_transport_rhs_2d( &
      species, transport, state, temperature, geometry, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, rhs, &
      minimum_theta, ok, exterior)
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
    type(reactive_eb_exterior_state_2d), intent(in), optional :: exterior

    real(dp), allocatable :: x_flux(:, :, :), y_flux(:, :, :)

    allocate(x_flux(size(state, 1), 0:geometry%nx, geometry%ny))
    allocate(y_flux(size(state, 1), geometry%nx, 0:geometry%ny))
    call reactive_eb_transport_fluxes_rhs_2d( &
      species, transport, state, temperature, geometry, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, rhs, &
      x_flux, y_flux, minimum_theta, ok, exterior)
  end subroutine reactive_eb_transport_rhs_2d

  subroutine reactive_eb_transport_euler_update_2d( &
      species, transport, state, temperature, geometry, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, max_order, new_state, new_temperature, &
      minimum_theta, ok, exterior)
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
    type(reactive_eb_exterior_state_2d), intent(in), optional :: exterior

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
      minimum_theta, local_ok, exterior)
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
