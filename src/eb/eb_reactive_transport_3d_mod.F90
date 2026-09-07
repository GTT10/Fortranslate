module eb_reactive_transport_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use state_indices_mod, only: iet
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_conserved_to_primitive
  use reactive_transport_3d_mod, only: &
    reactive_transport_face_flux_3d, reactive_transport_timestep_3d
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d
  use eb_reactive_redistribution_3d_mod, only: &
    advance_reactive_eb_state_redistributed_3d
  implicit none
  private

  real(dp), parameter :: species_safety = 0.90_dp

  public :: reactive_eb_transport_timestep_3d
  public :: reactive_eb_transport_fluxes_rhs_3d
  public :: reactive_eb_transport_rhs_3d
  public :: reactive_eb_transport_euler_update_3d
  public :: advance_reactive_eb_transport_3d

contains

  pure logical function valid_transport_shapes( &
      species, transport, state, temperature, geometry) result(valid)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry

    integer :: nvar

    nvar = reactive_nvar(size(species))
    valid = nvar > 0 .and. size(transport) == size(species) .and. &
      geometry%is_valid() .and. size(state, 1) == nvar .and. &
      size(state, 2) == geometry%nx .and. &
      size(state, 3) == geometry%ny .and. &
      size(state, 4) == geometry%nz .and. &
      all(shape(temperature) == &
        [geometry%nx, geometry%ny, geometry%nz])
  end function valid_transport_shapes

  subroutine recover_active_primitives( &
      species, state, temperature, geometry, primitive, &
      checked_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(out) :: primitive(:, :, :, :)
    real(dp), intent(out) :: checked_temperature(:, :, :)
    logical, intent(out) :: ok

    real(dp) :: sound_speed
    logical :: local_ok
    integer :: i, j, k, nprim

    primitive = 0.0_dp
    checked_temperature = 0.0_dp
    ok = .false.
    nprim = reactive_nprim(size(species))
    if (size(primitive, 1) /= nprim .or. &
        size(primitive, 2) /= geometry%nx .or. &
        size(primitive, 3) /= geometry%ny .or. &
        size(primitive, 4) /= geometry%nz .or. &
        any(shape(checked_temperature) /= &
          [geometry%nx, geometry%ny, geometry%nz])) return
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), &
            primitive(:, i, j, k), checked_temperature(i, j, k), &
            sound_speed, local_ok)
          if (.not. local_ok) return
        end do
      end do
    end do
    ok = .true.
  end subroutine recover_active_primitives

  pure real(dp) function active_velocity_derivative( &
      primitive, geometry, i, j, k, velocity_component, direction) &
      result(derivative)
    real(dp), intent(in) :: primitive(:, :, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: i, j, k, velocity_component, direction

    integer :: lower_i, lower_j, lower_k, upper_i, upper_j, upper_k
    integer :: primitive_component
    real(dp) :: spacing
    logical :: has_lower, has_upper

    derivative = 0.0_dp
    primitive_component = velocity_component + 1
    lower_i = i
    lower_j = j
    lower_k = k
    upper_i = i
    upper_j = j
    upper_k = k
    select case (direction)
    case (1)
      lower_i = i - 1
      upper_i = i + 1
      spacing = geometry%dx
      has_lower = i > 1 .and. &
        geometry%x_face_fraction(i - 1, j, k) > 0.0_dp
      has_upper = i < geometry%nx .and. &
        geometry%x_face_fraction(i, j, k) > 0.0_dp
    case (2)
      lower_j = j - 1
      upper_j = j + 1
      spacing = geometry%dy
      has_lower = j > 1 .and. &
        geometry%y_face_fraction(i, j - 1, k) > 0.0_dp
      has_upper = j < geometry%ny .and. &
        geometry%y_face_fraction(i, j, k) > 0.0_dp
    case (3)
      lower_k = k - 1
      upper_k = k + 1
      spacing = geometry%dz
      has_lower = k > 1 .and. &
        geometry%z_face_fraction(i, j, k - 1) > 0.0_dp
      has_upper = k < geometry%nz .and. &
        geometry%z_face_fraction(i, j, k) > 0.0_dp
    case default
      return
    end select
    if (has_lower) has_lower = geometry%cell_type( &
      lower_i, lower_j, lower_k) /= eb_covered_cell_3d
    if (has_upper) has_upper = geometry%cell_type( &
      upper_i, upper_j, upper_k) /= eb_covered_cell_3d
    if (has_lower .and. has_upper) then
      derivative = ( &
        primitive(primitive_component, upper_i, upper_j, upper_k) - &
        primitive(primitive_component, lower_i, lower_j, lower_k)) / &
        (2.0_dp * spacing)
    else if (has_upper) then
      derivative = ( &
        primitive(primitive_component, upper_i, upper_j, upper_k) - &
        primitive(primitive_component, i, j, k)) / spacing
    else if (has_lower) then
      derivative = ( &
        primitive(primitive_component, i, j, k) - &
        primitive(primitive_component, lower_i, lower_j, lower_k)) / spacing
    end if
  end function active_velocity_derivative

  pure subroutine active_face_velocity_gradient( &
      primitive, geometry, left_i, left_j, left_k, right_i, right_j, &
      right_k, normal, gradient)
    real(dp), intent(in) :: primitive(:, :, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: left_i, left_j, left_k
    integer, intent(in) :: right_i, right_j, right_k, normal
    real(dp), intent(out) :: gradient(3, 3)

    real(dp) :: spacing(3)
    integer :: velocity_component, direction

    spacing = [geometry%dx, geometry%dy, geometry%dz]
    gradient = 0.0_dp
    do velocity_component = 1, 3
      gradient(velocity_component, normal) = ( &
        primitive(velocity_component + 1, right_i, right_j, right_k) - &
        primitive(velocity_component + 1, left_i, left_j, left_k)) / &
        spacing(normal)
      do direction = 1, 3
        if (direction == normal) cycle
        gradient(velocity_component, direction) = 0.5_dp * ( &
          active_velocity_derivative( &
            primitive, geometry, left_i, left_j, left_k, &
            velocity_component, direction) + &
          active_velocity_derivative( &
            primitive, geometry, right_i, right_j, right_k, &
            velocity_component, direction))
      end do
    end do
  end subroutine active_face_velocity_gradient

  subroutine reactive_eb_transport_timestep_3d( &
      species, transport, state, temperature, geometry, transport_cfl, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, dt, maximum_diffusivity, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: transport_cfl
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled
    real(dp), intent(out) :: dt, maximum_diffusivity
    logical, intent(out) :: ok

    real(dp), allocatable :: active_state(:, :, :, :)
    real(dp), allocatable :: active_temperature(:, :, :)
    logical :: found_active
    integer :: i, j, k, reference_i, reference_j, reference_k

    dt = 0.0_dp
    maximum_diffusivity = 0.0_dp
    ok = .false.
    if (.not. valid_transport_shapes( &
          species, transport, state, temperature, geometry)) return
    found_active = .false.
    reference_i = 0
    reference_j = 0
    reference_k = 0
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          reference_i = i
          reference_j = j
          reference_k = k
          found_active = .true.
          exit
        end do
        if (found_active) exit
      end do
      if (found_active) exit
    end do
    if (.not. found_active) return
    allocate(active_state, source=state)
    allocate(active_temperature, source=temperature)
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) /= eb_covered_cell_3d) cycle
          active_state(:, i, j, k) = &
            state(:, reference_i, reference_j, reference_k)
          active_temperature(i, j, k) = &
            temperature(reference_i, reference_j, reference_k)
        end do
      end do
    end do
    call reactive_transport_timestep_3d( &
      species, transport, active_state, active_temperature, &
      geometry%nx, geometry%ny, geometry%nz, geometry%dx, geometry%dy, &
      geometry%dz, transport_cfl, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, dt, &
      maximum_diffusivity, ok)
  end subroutine reactive_eb_transport_timestep_3d

  pure subroutine limit_species_face( &
      flux, species_energy, nspecies, theta)
    real(dp), intent(inout) :: flux(:)
    real(dp), intent(in) :: species_energy, theta
    integer, intent(in) :: nspecies

    integer :: species_index, component

    do species_index = 1, nspecies
      component = reactive_species_component(species_index)
      flux(component) = theta * flux(component)
    end do
    flux(iet) = flux(iet) + (theta - 1.0_dp) * species_energy
  end subroutine limit_species_face

  subroutine limit_reactive_eb_transport_fluxes_3d( &
      species, state, geometry, dt, x_flux, y_flux, z_flux, &
      species_energy_x, species_energy_y, species_energy_z, &
      minimum_theta, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), dt
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(inout) :: x_flux(:, 0:, :, :)
    real(dp), intent(inout) :: y_flux(:, :, 0:, :)
    real(dp), intent(inout) :: z_flux(:, :, :, 0:)
    real(dp), intent(in) :: species_energy_x(0:, :, :)
    real(dp), intent(in) :: species_energy_y(:, 0:, :)
    real(dp), intent(in) :: species_energy_z(:, :, 0:)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok

    real(dp), allocatable :: theta_cell(:, :, :)
    real(dp) :: outgoing, available, candidate, theta_face
    real(dp) :: x_area, y_area, z_area, cell_volume
    integer :: i, j, k, face_i, face_j, face_k
    integer :: species_index, component

    minimum_theta = 1.0_dp
    ok = .false.
    if (.not. ieee_is_finite(dt) .or. dt < 0.0_dp) return
    if (dt <= tiny(1.0_dp)) then
      ok = .true.
      return
    end if
    allocate(theta_cell(geometry%nx, geometry%ny, geometry%nz))
    theta_cell = 1.0_dp
    x_area = geometry%dy * geometry%dz
    y_area = geometry%dx * geometry%dz
    z_area = geometry%dx * geometry%dy
    cell_volume = geometry%dx * geometry%dy * geometry%dz
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          do species_index = 1, size(species)
            component = reactive_species_component(species_index)
            outgoing = x_area * ( &
              max(geometry%x_face_fraction(i, j, k) * &
                x_flux(component, i, j, k), 0.0_dp) + &
              max(-geometry%x_face_fraction(i - 1, j, k) * &
                x_flux(component, i - 1, j, k), 0.0_dp)) + &
              y_area * ( &
              max(geometry%y_face_fraction(i, j, k) * &
                y_flux(component, i, j, k), 0.0_dp) + &
              max(-geometry%y_face_fraction(i, j - 1, k) * &
                y_flux(component, i, j - 1, k), 0.0_dp)) + &
              z_area * ( &
              max(geometry%z_face_fraction(i, j, k) * &
                z_flux(component, i, j, k), 0.0_dp) + &
              max(-geometry%z_face_fraction(i, j, k - 1) * &
                z_flux(component, i, j, k - 1), 0.0_dp))
            available = max(0.0_dp, state(component, i, j, k)) * &
              geometry%volume_fraction(i, j, k) * cell_volume
            if (outgoing > 0.0_dp) then
              candidate = species_safety * available / (dt * outgoing)
              theta_cell(i, j, k) = min(theta_cell(i, j, k), &
                max(0.0_dp, min(1.0_dp, candidate)))
            end if
          end do
        end do
      end do
    end do

    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do face_i = 0, geometry%nx
          theta_face = 1.0_dp
          if (face_i >= 1) then
            if (geometry%cell_type(face_i, j, k) /= eb_covered_cell_3d) &
              theta_face = min(theta_face, theta_cell(face_i, j, k))
          end if
          if (face_i + 1 <= geometry%nx) then
            if (geometry%cell_type(face_i + 1, j, k) /= &
                eb_covered_cell_3d) theta_face = min( &
                  theta_face, theta_cell(face_i + 1, j, k))
          end if
          call limit_species_face( &
            x_flux(:, face_i, j, k), species_energy_x(face_i, j, k), &
            size(species), theta_face)
          minimum_theta = min(minimum_theta, theta_face)
        end do
      end do
    end do
    do k = 1, geometry%nz
      do face_j = 0, geometry%ny
        do i = 1, geometry%nx
          theta_face = 1.0_dp
          if (face_j >= 1) then
            if (geometry%cell_type(i, face_j, k) /= eb_covered_cell_3d) &
              theta_face = min(theta_face, theta_cell(i, face_j, k))
          end if
          if (face_j + 1 <= geometry%ny) then
            if (geometry%cell_type(i, face_j + 1, k) /= &
                eb_covered_cell_3d) theta_face = min( &
                  theta_face, theta_cell(i, face_j + 1, k))
          end if
          call limit_species_face( &
            y_flux(:, i, face_j, k), species_energy_y(i, face_j, k), &
            size(species), theta_face)
          minimum_theta = min(minimum_theta, theta_face)
        end do
      end do
    end do
    do face_k = 0, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          theta_face = 1.0_dp
          if (face_k >= 1) then
            if (geometry%cell_type(i, j, face_k) /= eb_covered_cell_3d) &
              theta_face = min(theta_face, theta_cell(i, j, face_k))
          end if
          if (face_k + 1 <= geometry%nz) then
            if (geometry%cell_type(i, j, face_k + 1) /= &
                eb_covered_cell_3d) theta_face = min( &
                  theta_face, theta_cell(i, j, face_k + 1))
          end if
          call limit_species_face( &
            z_flux(:, i, j, face_k), species_energy_z(i, j, face_k), &
            size(species), theta_face)
          minimum_theta = min(minimum_theta, theta_face)
        end do
      end do
    end do
    ok = all(ieee_is_finite(x_flux)) .and. &
      all(ieee_is_finite(y_flux)) .and. &
      all(ieee_is_finite(z_flux)) .and. &
      ieee_is_finite(minimum_theta) .and. &
      minimum_theta >= 0.0_dp .and. minimum_theta <= 1.0_dp
  end subroutine limit_reactive_eb_transport_fluxes_3d

  subroutine reactive_eb_transport_fluxes_rhs_3d( &
      species, transport, state, temperature, geometry, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, rhs, &
      x_flux, y_flux, z_flux, minimum_theta, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: dt
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: rhs(:, :, :, :)
    real(dp), intent(out) :: x_flux(:, 0:, :, :)
    real(dp), intent(out) :: y_flux(:, :, 0:, :)
    real(dp), intent(out) :: z_flux(:, :, :, 0:)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:, :, :, :)
    real(dp), allocatable :: checked_temperature(:, :, :)
    real(dp), allocatable :: species_energy_x(:, :, :)
    real(dp), allocatable :: species_energy_y(:, :, :)
    real(dp), allocatable :: species_energy_z(:, :, :)
    real(dp) :: gradient(3, 3), fluid_volume
    real(dp) :: x_area, y_area, z_area
    logical :: local_ok
    integer :: i, j, k, nvar, nprim

    rhs = 0.0_dp
    x_flux = 0.0_dp
    y_flux = 0.0_dp
    z_flux = 0.0_dp
    minimum_theta = 1.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    nprim = reactive_nprim(size(species))
    if (.not. valid_transport_shapes( &
          species, transport, state, temperature, geometry) .or. &
        any(shape(rhs) /= shape(state)) .or. &
        size(x_flux, 1) /= nvar .or. &
        size(x_flux, 2) /= geometry%nx + 1 .or. &
        size(x_flux, 3) /= geometry%ny .or. &
        size(x_flux, 4) /= geometry%nz .or. &
        size(y_flux, 1) /= nvar .or. &
        size(y_flux, 2) /= geometry%nx .or. &
        size(y_flux, 3) /= geometry%ny + 1 .or. &
        size(y_flux, 4) /= geometry%nz .or. &
        size(z_flux, 1) /= nvar .or. &
        size(z_flux, 2) /= geometry%nx .or. &
        size(z_flux, 3) /= geometry%ny .or. &
        size(z_flux, 4) /= geometry%nz + 1 .or. &
        .not. ieee_is_finite(dt) .or. dt < 0.0_dp .or. &
        (barodiffusion_enabled .and. &
          .not. species_diffusion_enabled)) return
    if (.not. (viscosity_enabled .or. thermal_conduction_enabled .or. &
        species_diffusion_enabled)) then
      ok = .true.
      return
    end if

    allocate(primitive(nprim, geometry%nx, geometry%ny, geometry%nz))
    allocate(checked_temperature(geometry%nx, geometry%ny, geometry%nz))
    allocate(species_energy_x(0:geometry%nx, geometry%ny, geometry%nz))
    allocate(species_energy_y(geometry%nx, 0:geometry%ny, geometry%nz))
    allocate(species_energy_z(geometry%nx, geometry%ny, 0:geometry%nz))
    species_energy_x = 0.0_dp
    species_energy_y = 0.0_dp
    species_energy_z = 0.0_dp
    call recover_active_primitives( &
      species, state, temperature, geometry, primitive, &
      checked_temperature, local_ok)
    if (.not. local_ok) return

    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx - 1
          if (geometry%x_face_fraction(i, j, k) <= 0.0_dp) cycle
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d .or. &
              geometry%cell_type(i + 1, j, k) == &
                eb_covered_cell_3d) return
          call active_face_velocity_gradient( &
            primitive, geometry, i, j, k, i + 1, j, k, 1, gradient)
          call reactive_transport_face_flux_3d( &
            species, transport, primitive(:, i, j, k), &
            primitive(:, i + 1, j, k), checked_temperature(i, j, k), &
            checked_temperature(i + 1, j, k), geometry%dx, 1, gradient, &
            viscosity_enabled, thermal_conduction_enabled, &
            species_diffusion_enabled, barodiffusion_enabled, &
            x_flux(:, i, j, k), species_energy_x(i, j, k), local_ok)
          if (.not. local_ok) return
        end do
      end do
    end do
    do k = 1, geometry%nz
      do j = 1, geometry%ny - 1
        do i = 1, geometry%nx
          if (geometry%y_face_fraction(i, j, k) <= 0.0_dp) cycle
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d .or. &
              geometry%cell_type(i, j + 1, k) == &
                eb_covered_cell_3d) return
          call active_face_velocity_gradient( &
            primitive, geometry, i, j, k, i, j + 1, k, 2, gradient)
          call reactive_transport_face_flux_3d( &
            species, transport, primitive(:, i, j, k), &
            primitive(:, i, j + 1, k), checked_temperature(i, j, k), &
            checked_temperature(i, j + 1, k), geometry%dy, 2, gradient, &
            viscosity_enabled, thermal_conduction_enabled, &
            species_diffusion_enabled, barodiffusion_enabled, &
            y_flux(:, i, j, k), species_energy_y(i, j, k), local_ok)
          if (.not. local_ok) return
        end do
      end do
    end do
    do k = 1, geometry%nz - 1
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%z_face_fraction(i, j, k) <= 0.0_dp) cycle
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d .or. &
              geometry%cell_type(i, j, k + 1) == &
                eb_covered_cell_3d) return
          call active_face_velocity_gradient( &
            primitive, geometry, i, j, k, i, j, k + 1, 3, gradient)
          call reactive_transport_face_flux_3d( &
            species, transport, primitive(:, i, j, k), &
            primitive(:, i, j, k + 1), checked_temperature(i, j, k), &
            checked_temperature(i, j, k + 1), geometry%dz, 3, gradient, &
            viscosity_enabled, thermal_conduction_enabled, &
            species_diffusion_enabled, barodiffusion_enabled, &
            z_flux(:, i, j, k), species_energy_z(i, j, k), local_ok)
          if (.not. local_ok) return
        end do
      end do
    end do
    call limit_reactive_eb_transport_fluxes_3d( &
      species, state, geometry, dt, x_flux, y_flux, z_flux, &
      species_energy_x, species_energy_y, species_energy_z, &
      minimum_theta, local_ok)
    if (.not. local_ok) return

    x_area = geometry%dy * geometry%dz
    y_area = geometry%dx * geometry%dz
    z_area = geometry%dx * geometry%dy
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          fluid_volume = geometry%volume_fraction(i, j, k) * &
            geometry%dx * geometry%dy * geometry%dz
          if (fluid_volume <= 0.0_dp) return
          rhs(:, i, j, k) = -( &
            x_area * ( &
              geometry%x_face_fraction(i, j, k) * &
                x_flux(:, i, j, k) - &
              geometry%x_face_fraction(i - 1, j, k) * &
                x_flux(:, i - 1, j, k)) + &
            y_area * ( &
              geometry%y_face_fraction(i, j, k) * &
                y_flux(:, i, j, k) - &
              geometry%y_face_fraction(i, j - 1, k) * &
                y_flux(:, i, j - 1, k)) + &
            z_area * ( &
              geometry%z_face_fraction(i, j, k) * &
                z_flux(:, i, j, k) - &
              geometry%z_face_fraction(i, j, k - 1) * &
                z_flux(:, i, j, k - 1))) / fluid_volume
        end do
      end do
    end do
    ok = all(ieee_is_finite(rhs))
  end subroutine reactive_eb_transport_fluxes_rhs_3d

  subroutine reactive_eb_transport_rhs_3d( &
      species, transport, state, temperature, geometry, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, rhs, &
      minimum_theta, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: dt
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: rhs(:, :, :, :)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok

    real(dp), allocatable :: x_flux(:, :, :, :)
    real(dp), allocatable :: y_flux(:, :, :, :)
    real(dp), allocatable :: z_flux(:, :, :, :)

    allocate(x_flux(size(state, 1), 0:geometry%nx, &
      geometry%ny, geometry%nz))
    allocate(y_flux(size(state, 1), geometry%nx, &
      0:geometry%ny, geometry%nz))
    allocate(z_flux(size(state, 1), geometry%nx, &
      geometry%ny, 0:geometry%nz))
    call reactive_eb_transport_fluxes_rhs_3d( &
      species, transport, state, temperature, geometry, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, rhs, &
      x_flux, y_flux, z_flux, minimum_theta, ok)
  end subroutine reactive_eb_transport_rhs_3d

  subroutine reactive_eb_transport_euler_update_3d( &
      species, transport, state, temperature, geometry, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, &
      target_volume_fraction, new_state, new_temperature, &
      minimum_theta, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: dt, target_volume_fraction
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: new_state(:, :, :, :)
    real(dp), intent(out) :: new_temperature(:, :, :)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok

    real(dp), allocatable :: rhs(:, :, :, :)
    logical :: local_ok

    minimum_theta = 1.0_dp
    ok = .false.
    if (.not. valid_transport_shapes( &
          species, transport, state, temperature, geometry) .or. &
        any(shape(new_state) /= shape(state)) .or. &
        any(shape(new_temperature) /= shape(temperature)) .or. &
        .not. ieee_is_finite(dt) .or. dt < 0.0_dp .or. &
        (barodiffusion_enabled .and. &
          .not. species_diffusion_enabled)) return
    new_state = state
    new_temperature = temperature
    if (dt <= tiny(1.0_dp) .or. .not. (viscosity_enabled .or. &
        thermal_conduction_enabled .or. species_diffusion_enabled)) then
      ok = .true.
      return
    end if
    allocate(rhs, mold=state)
    call reactive_eb_transport_rhs_3d( &
      species, transport, state, temperature, geometry, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, rhs, &
      minimum_theta, local_ok)
    if (.not. local_ok) return
    call advance_reactive_eb_state_redistributed_3d( &
      species, state, temperature, geometry, rhs, dt, new_state, &
      new_temperature, local_ok, target_volume_fraction)
    if (.not. local_ok) return
    ok = .true.
  end subroutine reactive_eb_transport_euler_update_3d

  subroutine advance_reactive_eb_transport_3d( &
      species, transport, state, temperature, geometry, interval, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, &
      target_volume_fraction, minimum_theta, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: interval, target_volume_fraction
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok

    real(dp), allocatable :: initial_state(:, :, :, :)
    real(dp), allocatable :: initial_temperature(:, :, :)
    real(dp), allocatable :: stage_state(:, :, :, :)
    real(dp), allocatable :: stage_temperature(:, :, :)
    real(dp), allocatable :: euler_state(:, :, :, :)
    real(dp), allocatable :: euler_temperature(:, :, :)
    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    real(dp), allocatable :: primitive(:)
    real(dp) :: theta_one, theta_two, recovered_temperature, sound_speed
    logical :: local_ok
    integer :: i, j, k

    minimum_theta = 1.0_dp
    ok = .false.
    if (.not. valid_transport_shapes( &
          species, transport, state, temperature, geometry) .or. &
        .not. ieee_is_finite(interval) .or. interval < 0.0_dp .or. &
        .not. ieee_is_finite(target_volume_fraction) .or. &
        target_volume_fraction <= 0.0_dp .or. &
        target_volume_fraction > 1.0_dp .or. &
        (barodiffusion_enabled .and. &
          .not. species_diffusion_enabled)) return
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
    call reactive_eb_transport_euler_update_3d( &
      species, transport, initial_state, initial_temperature, geometry, &
      interval, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, &
      target_volume_fraction, stage_state, stage_temperature, theta_one, &
      local_ok)
    if (.not. local_ok) return
    call reactive_eb_transport_euler_update_3d( &
      species, transport, stage_state, stage_temperature, geometry, &
      interval, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, &
      target_volume_fraction, euler_state, euler_temperature, theta_two, &
      local_ok)
    if (.not. local_ok) return

    allocate(candidate_state, source=0.5_dp * (initial_state + euler_state))
    allocate(candidate_temperature, source=initial_temperature)
    allocate(primitive(reactive_nprim(size(species))))
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) then
            candidate_state(:, i, j, k) = initial_state(:, i, j, k)
            candidate_temperature(i, j, k) = initial_temperature(i, j, k)
            cycle
          end if
          call reactive_conserved_to_primitive( &
            species, candidate_state(:, i, j, k), &
            0.5_dp * (initial_temperature(i, j, k) + &
              euler_temperature(i, j, k)), primitive, &
            recovered_temperature, sound_speed, local_ok)
          if (.not. local_ok) return
          candidate_temperature(i, j, k) = recovered_temperature
        end do
      end do
    end do
    if (any(.not. ieee_is_finite(candidate_state)) .or. &
        any(.not. ieee_is_finite(candidate_temperature))) return
    state = candidate_state
    temperature = candidate_temperature
    minimum_theta = min(theta_one, theta_two)
    ok = .true.
  end subroutine advance_reactive_eb_transport_3d

end module eb_reactive_transport_3d_mod
