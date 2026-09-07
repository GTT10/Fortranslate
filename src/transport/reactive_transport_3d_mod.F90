module reactive_transport_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use mixture_thermo_mod, only: mixture_mass_properties
  use gas_transport_mod, only: gas_transport_species
  use mixture_transport_mod, only: mixture_transport_coefficients
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_mass_fraction_component, reactive_conserved_to_primitive
  use reactive_transport_2d_mod, only: &
    face_transport_data, species_face_flux
  implicit none
  private

  real(dp), parameter :: species_safety = 0.90_dp

  public :: reactive_transport_fluxes_3d
  public :: reactive_transport_ghosted_fluxes_3d
  public :: reactive_transport_face_flux_3d
  public :: reactive_transport_interface_theta_3d
  public :: reactive_transport_timestep_3d
  public :: reactive_transport_euler_update_3d
  public :: advance_reactive_transport_3d

contains

  pure logical function finite_product_fits(left, right) result(fits)
    real(dp), intent(in) :: left, right

    fits = .false.
    if (.not. ieee_is_finite(left)) return
    if (.not. ieee_is_finite(right)) return
    if (abs(left) <= 1.0_dp .or. abs(right) <= 1.0_dp) then
      fits = .true.
      return
    end if
    fits = abs(left) <= huge(1.0_dp) / abs(right)
  end function finite_product_fits

  pure logical function finite_quotient_fits(numerator, denominator) &
      result(fits)
    real(dp), intent(in) :: numerator, denominator

    fits = .false.
    if (.not. ieee_is_finite(numerator)) return
    if (.not. ieee_is_finite(denominator)) return
    if (abs(denominator) <= 0.0_dp) return
    if (abs(denominator) >= 1.0_dp) then
      fits = .true.
      return
    end if
    fits = abs(numerator) / huge(1.0_dp) <= abs(denominator)
  end function finite_quotient_fits

  pure logical function finite_sum_fits(left, right) result(fits)
    real(dp), intent(in) :: left, right

    fits = .false.
    if (.not. ieee_is_finite(left)) return
    if (.not. ieee_is_finite(right)) return
    if (right > 0.0_dp) then
      if (left > huge(1.0_dp) - right) return
    else if (right < 0.0_dp) then
      if (left < -huge(1.0_dp) - right) return
    end if
    fits = .true.
  end function finite_sum_fits

  pure subroutine checked_sum(left, right, value, ok)
    real(dp), intent(in) :: left, right
    real(dp), intent(out) :: value
    logical, intent(out) :: ok

    value = 0.0_dp
    ok = .false.
    if (.not. finite_sum_fits(left, right)) return
    value = left + right
    if (.not. ieee_is_finite(value)) then
      value = 0.0_dp
      return
    end if
    ok = .true.
  end subroutine checked_sum

  pure subroutine checked_difference(left, right, value, ok)
    real(dp), intent(in) :: left, right
    real(dp), intent(out) :: value
    logical, intent(out) :: ok

    value = 0.0_dp
    ok = .false.
    if (.not. ieee_is_finite(left)) return
    if (.not. ieee_is_finite(right)) return
    if (.not. finite_sum_fits(left, -right)) return
    value = left - right
    if (.not. ieee_is_finite(value)) then
      value = 0.0_dp
      return
    end if
    ok = .true.
  end subroutine checked_difference

  pure subroutine checked_product(left, right, value, ok)
    real(dp), intent(in) :: left, right
    real(dp), intent(out) :: value
    logical, intent(out) :: ok

    value = 0.0_dp
    ok = .false.
    if (.not. finite_product_fits(left, right)) return
    value = left * right
    if (.not. ieee_is_finite(value)) then
      value = 0.0_dp
      return
    end if
    ok = .true.
  end subroutine checked_product

  pure subroutine checked_quotient(numerator, denominator, value, ok)
    real(dp), intent(in) :: numerator, denominator
    real(dp), intent(out) :: value
    logical, intent(out) :: ok

    value = 0.0_dp
    ok = .false.
    if (.not. finite_quotient_fits(numerator, denominator)) return
    value = numerator / denominator
    if (.not. ieee_is_finite(value)) then
      value = 0.0_dp
      return
    end if
    ok = .true.
  end subroutine checked_quotient

  pure subroutine checked_midpoint(left, right, value, ok)
    real(dp), intent(in) :: left, right
    real(dp), intent(out) :: value
    logical, intent(out) :: ok

    real(dp) :: left_half, right_half

    value = 0.0_dp
    ok = .false.
    if (finite_sum_fits(left, right)) then
      value = 0.5_dp * (left + right)
      if (.not. ieee_is_finite(value)) then
        value = 0.0_dp
        return
      end if
      ok = .true.
      return
    end if
    if (.not. finite_product_fits(0.5_dp, left)) return
    if (.not. finite_product_fits(0.5_dp, right)) return
    left_half = 0.5_dp * left
    right_half = 0.5_dp * right
    call checked_sum(left_half, right_half, value, ok)
  end subroutine checked_midpoint

  pure subroutine checked_dot_product(left, right, value, ok)
    real(dp), intent(in) :: left(:), right(:)
    real(dp), intent(out) :: value
    logical, intent(out) :: ok

    real(dp) :: absolute_bound, product
    integer :: i

    value = 0.0_dp
    ok = .false.
    if (size(left) /= size(right)) return
    absolute_bound = 0.0_dp
    do i = 1, size(left)
      if (.not. finite_product_fits(left(i), right(i))) return
      product = left(i) * right(i)
      if (.not. finite_sum_fits(absolute_bound, abs(product))) return
      absolute_bound = absolute_bound + abs(product)
    end do
    value = dot_product(left, right)
    if (.not. ieee_is_finite(value)) then
      value = 0.0_dp
      return
    end if
    ok = .true.
  end subroutine checked_dot_product

  pure subroutine checked_species_outgoing( &
      positive_x, negative_x, positive_y, negative_y, positive_z, negative_z, &
      dx, dy, dz, outgoing, ok)
    real(dp), intent(in) :: positive_x, negative_x, positive_y, negative_y
    real(dp), intent(in) :: positive_z, negative_z, dx, dy, dz
    real(dp), intent(out) :: outgoing
    logical, intent(out) :: ok

    real(dp) :: term, partial
    logical :: local_ok

    outgoing = 0.0_dp
    ok = .false.
    if (.not. all(ieee_is_finite([positive_x, negative_x, positive_y, &
        negative_y, positive_z, negative_z, dx, dy, dz]))) return
    if (dx <= 0.0_dp .or. dy <= 0.0_dp .or. dz <= 0.0_dp) return

    call checked_quotient(max(positive_x, 0.0_dp), dx, term, local_ok)
    if (.not. local_ok) return
    outgoing = term
    call checked_quotient(max(-negative_x, 0.0_dp), dx, term, local_ok)
    if (.not. local_ok) return
    call checked_sum(outgoing, term, partial, local_ok)
    if (.not. local_ok) return
    outgoing = partial
    call checked_quotient(max(positive_y, 0.0_dp), dy, term, local_ok)
    if (.not. local_ok) return
    call checked_sum(outgoing, term, partial, local_ok)
    if (.not. local_ok) return
    outgoing = partial
    call checked_quotient(max(-negative_y, 0.0_dp), dy, term, local_ok)
    if (.not. local_ok) return
    call checked_sum(outgoing, term, partial, local_ok)
    if (.not. local_ok) return
    outgoing = partial
    call checked_quotient(max(positive_z, 0.0_dp), dz, term, local_ok)
    if (.not. local_ok) return
    call checked_sum(outgoing, term, partial, local_ok)
    if (.not. local_ok) return
    outgoing = partial
    call checked_quotient(max(-negative_z, 0.0_dp), dz, term, local_ok)
    if (.not. local_ok) return
    call checked_sum(outgoing, term, partial, local_ok)
    if (.not. local_ok) return
    outgoing = partial
    ok = .true.
  end subroutine checked_species_outgoing

  subroutine reset_transport_flux_outputs_3d( &
      flux_x, flux_y, flux_z, minimum_theta)
    real(dp), intent(out) :: flux_x(:, :, :, :)
    real(dp), intent(out) :: flux_y(:, :, :, :)
    real(dp), intent(out) :: flux_z(:, :, :, :)
    real(dp), intent(out) :: minimum_theta

    flux_x = 0.0_dp
    flux_y = 0.0_dp
    flux_z = 0.0_dp
    minimum_theta = 1.0_dp
  end subroutine reset_transport_flux_outputs_3d

  subroutine reset_transport_face_outputs_3d(flux, species_energy)
    real(dp), intent(out) :: flux(:)
    real(dp), intent(out) :: species_energy

    flux = 0.0_dp
    species_energy = 0.0_dp
  end subroutine reset_transport_face_outputs_3d

  subroutine reset_transport_timestep_outputs_3d(dt, maximum_diffusivity)
    real(dp), intent(out) :: dt, maximum_diffusivity

    dt = 0.0_dp
    maximum_diffusivity = 0.0_dp
  end subroutine reset_transport_timestep_outputs_3d

  subroutine reset_transport_euler_outputs_3d( &
      output_state, output_temperature, minimum_theta)
    real(dp), intent(out) :: output_state(:, :, :, :)
    real(dp), intent(out) :: output_temperature(:, :, :)
    real(dp), intent(out) :: minimum_theta

    output_state = 0.0_dp
    output_temperature = 0.0_dp
    minimum_theta = 1.0_dp
  end subroutine reset_transport_euler_outputs_3d

  subroutine reset_transport_ghosted_outputs_3d( &
      flux_x, flux_y, flux_z, minimum_theta, theta_cell_output)
    real(dp), intent(out) :: flux_x(:, 0:, :, :)
    real(dp), intent(out) :: flux_y(:, :, 0:, :)
    real(dp), intent(out) :: flux_z(:, :, :, 0:)
    real(dp), intent(out) :: minimum_theta
    real(dp), intent(out), optional :: theta_cell_output(:, :, :)

    flux_x = 0.0_dp
    flux_y = 0.0_dp
    flux_z = 0.0_dp
    minimum_theta = 1.0_dp
    if (present(theta_cell_output)) theta_cell_output = 1.0_dp
  end subroutine reset_transport_ghosted_outputs_3d

  pure subroutine reactive_transport_interface_theta_3d( &
      nspecies, base_state, face_flux, signed_coefficient, theta, ok)
    integer, intent(in) :: nspecies
    real(dp), intent(in) :: base_state(:), face_flux(:)
    real(dp), intent(in) :: signed_coefficient
    real(dp), intent(out) :: theta
    logical, intent(out) :: ok

    real(dp) :: loss, loss_numerator, candidate_theta
    integer :: component, species_index

    theta = 0.0_dp
    ok = .false.
    if (nspecies < 1 .or. nspecies > 32 .or. &
        size(base_state) /= reactive_nvar(nspecies) .or. &
        size(face_flux) /= reactive_nvar(nspecies) .or. &
        .not. ieee_is_finite(signed_coefficient) .or. &
        .not. all(ieee_is_finite(base_state)) .or. &
        .not. all(ieee_is_finite(face_flux))) return

    candidate_theta = 1.0_dp
    do species_index = 1, nspecies
      component = reactive_species_component(species_index)
      if (base_state(component) < 0.0_dp) return
      if (.not. finite_product_fits( &
          signed_coefficient, face_flux(component))) return
      loss = -signed_coefficient * face_flux(component)
      if (.not. ieee_is_finite(loss)) return
      loss = max(loss, 0.0_dp)
      if (loss > 0.0_dp) then
        loss_numerator = species_safety * base_state(component)
        if (.not. finite_quotient_fits(loss_numerator, loss)) return
        loss = loss_numerator / loss
        if (.not. ieee_is_finite(loss)) return
        candidate_theta = min( &
          candidate_theta, max(0.0_dp, min(1.0_dp, loss)))
      end if
    end do
    if (.not. ieee_is_finite(candidate_theta)) return
    if (candidate_theta < 0.0_dp .or. candidate_theta > 1.0_dp) return
    theta = candidate_theta
    ok = .true.
  end subroutine reactive_transport_interface_theta_3d

  pure logical function valid_transport_shapes_3d( &
      state, temperature, nvar, nx, ny, nz) result(valid)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nvar, nx, ny, nz

    valid = nvar >= 1 .and. nx >= 1 .and. ny >= 1 .and. nz >= 1 .and. &
      size(state, 1) == nvar .and. size(state, 2) == nx .and. &
      size(state, 3) == ny .and. size(state, 4) == nz .and. &
      size(temperature, 1) == nx .and. size(temperature, 2) == ny .and. &
      size(temperature, 3) == nz
  end function valid_transport_shapes_3d

  pure integer function periodic_index(index, extent) result(wrapped)
    integer, intent(in) :: index, extent
    wrapped = 1 + modulo(index - 1, extent)
  end function periodic_index

  subroutine recover_transport_primitives_3d( &
      species, state, temperature, nx, ny, nz, primitive, &
      checked_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(out) :: primitive(:, :, :, :)
    real(dp), intent(out) :: checked_temperature(:, :, :)
    logical, intent(out) :: ok

    real(dp) :: sound_speed
    logical :: local_ok
    integer :: i, j, k, nvar, nprim

    primitive = 0.0_dp
    checked_temperature = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    nprim = reactive_nprim(size(species))
    if (.not. valid_transport_shapes_3d( &
          state, temperature, nvar, nx, ny, nz)) return
    if (size(primitive, 1) /= nprim .or. &
        size(primitive, 2) /= nx .or. size(primitive, 3) /= ny .or. &
        size(primitive, 4) /= nz .or. &
        size(checked_temperature, 1) /= nx .or. &
        size(checked_temperature, 2) /= ny .or. &
        size(checked_temperature, 3) /= nz) return
    if (.not. all(ieee_is_finite(state)) .or. &
        .not. all(ieee_is_finite(temperature))) return

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), &
            primitive(:, i, j, k), checked_temperature(i, j, k), &
            sound_speed, local_ok)
          if (.not. local_ok) return
        end do
      end do
    end do
    if (.not. all(ieee_is_finite(primitive)) .or. &
        .not. all(ieee_is_finite(checked_temperature))) return
    ok = .true.
  end subroutine recover_transport_primitives_3d

  pure subroutine centered_velocity_derivative( &
      primitive, nx, ny, nz, i, j, k, velocity_component, direction, &
      dx, dy, dz, derivative, ok)
    real(dp), intent(in) :: primitive(:, :, :, :)
    integer, intent(in) :: nx, ny, nz, i, j, k
    integer, intent(in) :: velocity_component, direction
    real(dp), intent(in) :: dx, dy, dz
    real(dp), intent(out) :: derivative
    logical, intent(out) :: ok

    real(dp) :: numerator, denominator
    integer :: lower, upper, primitive_component

    derivative = 0.0_dp
    ok = .false.
    if (velocity_component < 1 .or. velocity_component > 3) return
    if (direction < 1 .or. direction > 3) return
    primitive_component = velocity_component + 1
    select case (direction)
    case (1)
      lower = periodic_index(i - 1, nx)
      upper = periodic_index(i + 1, nx)
      call checked_difference( &
        primitive(primitive_component, upper, j, k), &
        primitive(primitive_component, lower, j, k), numerator, ok)
      if (.not. ok) return
      call checked_product(2.0_dp, dx, denominator, ok)
      if (.not. ok) return
      call checked_quotient(numerator, denominator, derivative, ok)
      if (.not. ok) return
      derivative = (primitive(primitive_component, upper, j, k) - &
        primitive(primitive_component, lower, j, k)) / (2.0_dp * dx)
    case (2)
      lower = periodic_index(j - 1, ny)
      upper = periodic_index(j + 1, ny)
      call checked_difference( &
        primitive(primitive_component, i, upper, k), &
        primitive(primitive_component, i, lower, k), numerator, ok)
      if (.not. ok) return
      call checked_product(2.0_dp, dy, denominator, ok)
      if (.not. ok) return
      call checked_quotient(numerator, denominator, derivative, ok)
      if (.not. ok) return
      derivative = (primitive(primitive_component, i, upper, k) - &
        primitive(primitive_component, i, lower, k)) / (2.0_dp * dy)
    case (3)
      lower = periodic_index(k - 1, nz)
      upper = periodic_index(k + 1, nz)
      call checked_difference( &
        primitive(primitive_component, i, j, upper), &
        primitive(primitive_component, i, j, lower), numerator, ok)
      if (.not. ok) return
      call checked_product(2.0_dp, dz, denominator, ok)
      if (.not. ok) return
      call checked_quotient(numerator, denominator, derivative, ok)
      if (.not. ok) return
      derivative = (primitive(primitive_component, i, j, upper) - &
        primitive(primitive_component, i, j, lower)) / (2.0_dp * dz)
    end select
  end subroutine centered_velocity_derivative

  pure subroutine face_velocity_gradient( &
      primitive, nx, ny, nz, left_i, left_j, left_k, &
      right_i, right_j, right_k, normal, dx, dy, dz, gradient, ok)
    real(dp), intent(in) :: primitive(:, :, :, :)
    integer, intent(in) :: nx, ny, nz
    integer, intent(in) :: left_i, left_j, left_k
    integer, intent(in) :: right_i, right_j, right_k, normal
    real(dp), intent(in) :: dx, dy, dz
    real(dp), intent(out) :: gradient(3, 3)
    logical, intent(out) :: ok

    real(dp) :: normal_difference, left_derivative, right_derivative
    real(dp) :: spacing(3)
    logical :: local_ok
    integer :: velocity_component, direction

    spacing = [dx, dy, dz]
    gradient = 0.0_dp
    ok = .false.
    if (normal < 1 .or. normal > 3) return
    do velocity_component = 1, 3
      call checked_difference( &
        primitive(velocity_component + 1, right_i, right_j, right_k), &
        primitive(velocity_component + 1, left_i, left_j, left_k), &
        normal_difference, local_ok)
      if (.not. local_ok) return
      call checked_quotient( &
        normal_difference, spacing(normal), &
        gradient(velocity_component, normal), local_ok)
      if (.not. local_ok) return
      gradient(velocity_component, normal) = &
        (primitive(velocity_component + 1, right_i, right_j, right_k) - &
         primitive(velocity_component + 1, left_i, left_j, left_k)) / &
        spacing(normal)
      do direction = 1, 3
        if (direction == normal) cycle
        call centered_velocity_derivative( &
          primitive, nx, ny, nz, left_i, left_j, left_k, &
          velocity_component, direction, dx, dy, dz, left_derivative, local_ok)
        if (.not. local_ok) return
        call centered_velocity_derivative( &
          primitive, nx, ny, nz, right_i, right_j, right_k, &
          velocity_component, direction, dx, dy, dz, right_derivative, local_ok)
        if (.not. local_ok) return
        call checked_midpoint( &
          left_derivative, right_derivative, &
          gradient(velocity_component, direction), local_ok)
        if (.not. local_ok) return
        gradient(velocity_component, direction) = &
          0.5_dp * (left_derivative + right_derivative)
      end do
    end do
    ok = .true.
  end subroutine face_velocity_gradient

  subroutine assemble_periodic_transport_face( &
      species, transport, left_primitive, right_primitive, &
      left_temperature, right_temperature, spacing, normal, gradient, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, flux, &
      species_energy, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: left_primitive(:), right_primitive(:)
    real(dp), intent(in) :: left_temperature, right_temperature, spacing
    integer, intent(in) :: normal
    real(dp), intent(in) :: gradient(3, 3)
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: flux(:), species_energy
    logical, intent(out) :: ok

    real(dp), allocatable :: diffusion(:), yleft(:), yright(:), yface(:)
    real(dp), allocatable :: xleft(:), xright(:), hface(:), species_flux(:)
    real(dp) :: viscosity, conductivity, density_face, pressure_face
    real(dp) :: divergence, tau_normal(3), velocity_face(3)
    real(dp) :: viscous_power
    real(dp) :: divergence_partial, strain_sum, correction_factor
    real(dp) :: normal_correction, temperature_difference, conductive_term
    logical :: local_ok
    integer :: component, species_index, nspecies

    flux = 0.0_dp
    species_energy = 0.0_dp
    ok = .false.
    nspecies = size(species)
    if (.not. all(ieee_is_finite(left_primitive)) .or. &
        .not. all(ieee_is_finite(right_primitive)) .or. &
        .not. ieee_is_finite(left_temperature) .or. &
        .not. ieee_is_finite(right_temperature) .or. &
        .not. ieee_is_finite(spacing) .or. &
        .not. all(ieee_is_finite(gradient))) return
    if (size(transport) /= nspecies .or. &
        size(left_primitive) /= reactive_nprim(nspecies) .or. &
        size(right_primitive) /= reactive_nprim(nspecies) .or. &
        size(flux) /= reactive_nvar(nspecies) .or. normal < 1 .or. &
        normal > 3 .or. spacing <= 0.0_dp .or. &
        (barodiffusion_enabled .and. &
          .not. species_diffusion_enabled)) return
    allocate(diffusion(nspecies), yleft(nspecies), yright(nspecies))
    allocate(yface(nspecies), xleft(nspecies), xright(nspecies))
    allocate(hface(nspecies), species_flux(nspecies))
    call face_transport_data( &
      species, transport, left_primitive, right_primitive, &
      left_temperature, right_temperature, viscosity, conductivity, &
      diffusion, yleft, yright, yface, xleft, xright, hface, local_ok)
    if (.not. local_ok) then
      call reset_transport_face_outputs_3d(flux, species_energy)
      return
    end if
    if (.not. ieee_is_finite(viscosity) .or. &
        .not. ieee_is_finite(conductivity) .or. &
        .not. all(ieee_is_finite(diffusion)) .or. &
        .not. all(ieee_is_finite(yleft)) .or. &
        .not. all(ieee_is_finite(yright)) .or. &
        .not. all(ieee_is_finite(yface)) .or. &
        .not. all(ieee_is_finite(xleft)) .or. &
        .not. all(ieee_is_finite(xright)) .or. &
        .not. all(ieee_is_finite(hface))) then
      call reset_transport_face_outputs_3d(flux, species_energy)
      return
    end if

    if (.not. finite_sum_fits( &
        left_primitive(1), right_primitive(1))) then
      call reset_transport_face_outputs_3d(flux, species_energy)
      return
    end if
    if (.not. finite_sum_fits( &
        left_primitive(5), right_primitive(5))) then
      call reset_transport_face_outputs_3d(flux, species_energy)
      return
    end if
    do component = 2, 4
      if (.not. finite_sum_fits( &
          left_primitive(component), right_primitive(component))) then
        call reset_transport_face_outputs_3d(flux, species_energy)
        return
      end if
    end do
    density_face = 0.5_dp * (left_primitive(1) + right_primitive(1))
    pressure_face = 0.5_dp * (left_primitive(5) + right_primitive(5))
    velocity_face = 0.5_dp * &
      (left_primitive(2:4) + right_primitive(2:4))
    if (.not. all(ieee_is_finite([density_face, pressure_face])) .or. &
        .not. all(ieee_is_finite(velocity_face))) then
      call reset_transport_face_outputs_3d(flux, species_energy)
      return
    end if
    if (viscosity_enabled) then
      if (.not. finite_sum_fits(gradient(1, 1), gradient(2, 2))) then
        call reset_transport_face_outputs_3d(flux, species_energy)
        return
      end if
      divergence_partial = gradient(1, 1) + gradient(2, 2)
      if (.not. finite_sum_fits(divergence_partial, gradient(3, 3))) then
        call reset_transport_face_outputs_3d(flux, species_energy)
        return
      end if
      divergence = gradient(1, 1) + gradient(2, 2) + gradient(3, 3)
      if (.not. ieee_is_finite(divergence)) then
        call reset_transport_face_outputs_3d(flux, species_energy)
        return
      end if
      do component = 1, 3
        if (.not. finite_sum_fits( &
            gradient(normal, component), gradient(component, normal))) then
          call reset_transport_face_outputs_3d(flux, species_energy)
          return
        end if
        strain_sum = gradient(normal, component) + gradient(component, normal)
        if (.not. finite_product_fits(viscosity, strain_sum)) then
          call reset_transport_face_outputs_3d(flux, species_energy)
          return
        end if
        tau_normal(component) = viscosity * &
          (gradient(normal, component) + gradient(component, normal))
        if (component == normal) then
          correction_factor = (2.0_dp / 3.0_dp) * viscosity
          if (.not. finite_product_fits( &
              correction_factor, divergence)) then
            call reset_transport_face_outputs_3d(flux, species_energy)
            return
          end if
          normal_correction = correction_factor * divergence
          if (.not. finite_sum_fits( &
              tau_normal(component), -normal_correction)) then
            call reset_transport_face_outputs_3d(flux, species_energy)
            return
          end if
          tau_normal(component) = tau_normal(component) - &
            (2.0_dp / 3.0_dp) * viscosity * divergence
        end if
      end do
      if (.not. all(ieee_is_finite(tau_normal))) then
        call reset_transport_face_outputs_3d(flux, species_energy)
        return
      end if
      flux(imx:imz) = -tau_normal
      call checked_dot_product(tau_normal, velocity_face, viscous_power, local_ok)
      if (.not. local_ok) then
        call reset_transport_face_outputs_3d(flux, species_energy)
        return
      end if
      flux(iet) = -viscous_power
    end if
    if (thermal_conduction_enabled) then
      if (.not. finite_sum_fits( &
          right_temperature, -left_temperature)) then
        call reset_transport_face_outputs_3d(flux, species_energy)
        return
      end if
      temperature_difference = right_temperature - left_temperature
      if (.not. finite_product_fits( &
          conductivity, temperature_difference)) then
        call reset_transport_face_outputs_3d(flux, species_energy)
        return
      end if
      conductive_term = conductivity * temperature_difference
      if (.not. finite_quotient_fits(conductive_term, spacing)) then
        call reset_transport_face_outputs_3d(flux, species_energy)
        return
      end if
      conductive_term = conductive_term / spacing
      if (.not. finite_sum_fits(flux(iet), -conductive_term)) then
        call reset_transport_face_outputs_3d(flux, species_energy)
        return
      end if
      flux(iet) = flux(iet) - conductivity * &
        (right_temperature - left_temperature) / spacing
    end if
    if (species_diffusion_enabled) then
      call species_face_flux( &
        species, diffusion, yleft, yright, yface, xleft, xright, hface, &
        density_face, left_primitive(5), right_primitive(5), pressure_face, &
        spacing, barodiffusion_enabled, species_flux, species_energy, local_ok)
      if (.not. local_ok) then
        call reset_transport_face_outputs_3d(flux, species_energy)
        return
      end if
      if (.not. all(ieee_is_finite(species_flux)) .or. &
          .not. ieee_is_finite(species_energy)) then
        call reset_transport_face_outputs_3d(flux, species_energy)
        return
      end if
      do species_index = 1, nspecies
        flux(reactive_species_component(species_index)) = &
          species_flux(species_index)
      end do
      if (.not. finite_sum_fits(flux(iet), species_energy)) then
        call reset_transport_face_outputs_3d(flux, species_energy)
        return
      end if
      flux(iet) = flux(iet) + species_energy
    end if
    flux(irho) = 0.0_dp
    if (.not. all(ieee_is_finite(flux)) .or. &
        .not. ieee_is_finite(species_energy)) then
      call reset_transport_face_outputs_3d(flux, species_energy)
      return
    end if
    ok = .true.
  end subroutine assemble_periodic_transport_face

  subroutine reactive_transport_face_flux_3d( &
      species, transport, left_primitive, right_primitive, &
      left_temperature, right_temperature, spacing, normal, gradient, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, flux, &
      species_energy, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: left_primitive(:), right_primitive(:)
    real(dp), intent(in) :: left_temperature, right_temperature, spacing
    integer, intent(in) :: normal
    real(dp), intent(in) :: gradient(3, 3)
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: flux(:), species_energy
    logical, intent(out) :: ok

    call assemble_periodic_transport_face( &
      species, transport, left_primitive, right_primitive, &
      left_temperature, right_temperature, spacing, normal, gradient, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, flux, &
      species_energy, ok)
  end subroutine reactive_transport_face_flux_3d

  subroutine reactive_transport_fluxes_3d( &
      species, transport, state, temperature, nx, ny, nz, dx, dy, dz, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, &
      flux_x, flux_y, flux_z, minimum_theta, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz, dt
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: flux_x(:, :, :, :)
    real(dp), intent(out) :: flux_y(:, :, :, :)
    real(dp), intent(out) :: flux_z(:, :, :, :)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:, :, :, :)
    real(dp), allocatable :: checked_temperature(:, :, :)
    real(dp), allocatable :: species_energy_x(:, :, :)
    real(dp), allocatable :: species_energy_y(:, :, :)
    real(dp), allocatable :: species_energy_z(:, :, :)
    real(dp), allocatable :: theta_cell(:, :, :)
    real(dp) :: gradient(3, 3), outgoing, mass, candidate, theta_face
    real(dp) :: theta_numerator, theta_denominator
    logical :: local_ok
    integer :: i, j, k, next_i, next_j, next_k
    integer :: previous_i, previous_j, previous_k
    integer :: component, species_index, nvar, nprim

    flux_x = 0.0_dp
    flux_y = 0.0_dp
    flux_z = 0.0_dp
    minimum_theta = 1.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    nprim = reactive_nprim(size(species))
    if (.not. all(ieee_is_finite([dx, dy, dz, dt]))) return
    if (barodiffusion_enabled .and. .not. species_diffusion_enabled) return
    if (size(transport) /= size(species) .or. nx < 2 .or. ny < 2 .or. &
        nz < 2) return
    if (dx <= 0.0_dp .or. dy <= 0.0_dp .or. dz <= 0.0_dp .or. dt < 0.0_dp) return
    if (.not. valid_transport_shapes_3d( &
          state, temperature, nvar, nx, ny, nz)) return
    if (any(shape(flux_x) /= shape(state)) .or. &
        any(shape(flux_y) /= shape(state)) .or. &
        any(shape(flux_z) /= shape(state))) return
    if (.not. all(ieee_is_finite(state)) .or. &
        .not. all(ieee_is_finite(temperature))) return
    if (.not. (viscosity_enabled .or. thermal_conduction_enabled .or. &
        species_diffusion_enabled)) then
      ok = .true.
      return
    end if

    allocate(primitive(nprim, nx, ny, nz))
    allocate(checked_temperature(nx, ny, nz))
    allocate(species_energy_x(nx, ny, nz))
    allocate(species_energy_y(nx, ny, nz))
    allocate(species_energy_z(nx, ny, nz))
    allocate(theta_cell(nx, ny, nz))
    species_energy_x = 0.0_dp
    species_energy_y = 0.0_dp
    species_energy_z = 0.0_dp
    call recover_transport_primitives_3d( &
      species, state, temperature, nx, ny, nz, primitive, &
      checked_temperature, local_ok)
    if (.not. local_ok) then
      call reset_transport_flux_outputs_3d( &
        flux_x, flux_y, flux_z, minimum_theta)
      return
    end if
    if (.not. all(ieee_is_finite(primitive)) .or. &
        .not. all(ieee_is_finite(checked_temperature))) then
      call reset_transport_flux_outputs_3d( &
        flux_x, flux_y, flux_z, minimum_theta)
      return
    end if

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          next_i = periodic_index(i + 1, nx)
          call face_velocity_gradient( &
            primitive, nx, ny, nz, i, j, k, next_i, j, k, 1, &
            dx, dy, dz, gradient, local_ok)
          if (.not. local_ok) then
            call reset_transport_flux_outputs_3d( &
              flux_x, flux_y, flux_z, minimum_theta)
            return
          end if
          call assemble_periodic_transport_face( &
            species, transport, primitive(:, i, j, k), &
            primitive(:, next_i, j, k), checked_temperature(i, j, k), &
            checked_temperature(next_i, j, k), dx, 1, gradient, &
            viscosity_enabled, thermal_conduction_enabled, &
            species_diffusion_enabled, barodiffusion_enabled, &
            flux_x(:, i, j, k), species_energy_x(i, j, k), local_ok)
          if (.not. local_ok) then
            call reset_transport_flux_outputs_3d( &
              flux_x, flux_y, flux_z, minimum_theta)
            return
          end if

          next_j = periodic_index(j + 1, ny)
          call face_velocity_gradient( &
            primitive, nx, ny, nz, i, j, k, i, next_j, k, 2, &
            dx, dy, dz, gradient, local_ok)
          if (.not. local_ok) then
            call reset_transport_flux_outputs_3d( &
              flux_x, flux_y, flux_z, minimum_theta)
            return
          end if
          call assemble_periodic_transport_face( &
            species, transport, primitive(:, i, j, k), &
            primitive(:, i, next_j, k), checked_temperature(i, j, k), &
            checked_temperature(i, next_j, k), dy, 2, gradient, &
            viscosity_enabled, thermal_conduction_enabled, &
            species_diffusion_enabled, barodiffusion_enabled, &
            flux_y(:, i, j, k), species_energy_y(i, j, k), local_ok)
          if (.not. local_ok) then
            call reset_transport_flux_outputs_3d( &
              flux_x, flux_y, flux_z, minimum_theta)
            return
          end if

          next_k = periodic_index(k + 1, nz)
          call face_velocity_gradient( &
            primitive, nx, ny, nz, i, j, k, i, j, next_k, 3, &
            dx, dy, dz, gradient, local_ok)
          if (.not. local_ok) then
            call reset_transport_flux_outputs_3d( &
              flux_x, flux_y, flux_z, minimum_theta)
            return
          end if
          call assemble_periodic_transport_face( &
            species, transport, primitive(:, i, j, k), &
            primitive(:, i, j, next_k), checked_temperature(i, j, k), &
            checked_temperature(i, j, next_k), dz, 3, gradient, &
            viscosity_enabled, thermal_conduction_enabled, &
            species_diffusion_enabled, barodiffusion_enabled, &
            flux_z(:, i, j, k), species_energy_z(i, j, k), local_ok)
          if (.not. local_ok) then
            call reset_transport_flux_outputs_3d( &
              flux_x, flux_y, flux_z, minimum_theta)
            return
          end if
        end do
      end do
    end do

    if (.not. all(ieee_is_finite(flux_x)) .or. &
        .not. all(ieee_is_finite(flux_y)) .or. &
        .not. all(ieee_is_finite(flux_z)) .or. &
        .not. all(ieee_is_finite(species_energy_x)) .or. &
        .not. all(ieee_is_finite(species_energy_y)) .or. &
        .not. all(ieee_is_finite(species_energy_z))) then
      call reset_transport_flux_outputs_3d( &
        flux_x, flux_y, flux_z, minimum_theta)
      return
    end if
    if (species_diffusion_enabled .and. dt > 0.0_dp) then
      theta_cell = 1.0_dp
      do k = 1, nz
        previous_k = periodic_index(k - 1, nz)
        do j = 1, ny
          previous_j = periodic_index(j - 1, ny)
          do i = 1, nx
            previous_i = periodic_index(i - 1, nx)
            do species_index = 1, size(species)
              component = reactive_species_component(species_index)
              call checked_species_outgoing( &
                flux_x(component, i, j, k), &
                flux_x(component, previous_i, j, k), &
                flux_y(component, i, j, k), &
                flux_y(component, i, previous_j, k), &
                flux_z(component, i, j, k), &
                flux_z(component, i, j, previous_k), dx, dy, dz, &
                outgoing, local_ok)
              if (.not. local_ok) then
                call reset_transport_flux_outputs_3d( &
                  flux_x, flux_y, flux_z, minimum_theta)
                return
              end if
              mass = max(0.0_dp, state(component, i, j, k))
              if (.not. ieee_is_finite(mass)) then
                call reset_transport_flux_outputs_3d( &
                  flux_x, flux_y, flux_z, minimum_theta)
                return
              end if
              if (outgoing > 0.0_dp) then
                call checked_product( &
                  species_safety, mass, theta_numerator, local_ok)
                if (.not. local_ok) then
                  call reset_transport_flux_outputs_3d( &
                    flux_x, flux_y, flux_z, minimum_theta)
                  return
                end if
                call checked_product(dt, outgoing, theta_denominator, local_ok)
                if (.not. local_ok) then
                  call reset_transport_flux_outputs_3d( &
                    flux_x, flux_y, flux_z, minimum_theta)
                  return
                end if
                call checked_quotient( &
                  theta_numerator, theta_denominator, candidate, local_ok)
                if (.not. local_ok) then
                  call reset_transport_flux_outputs_3d( &
                    flux_x, flux_y, flux_z, minimum_theta)
                  return
                end if
                theta_cell(i, j, k) = min(theta_cell(i, j, k), &
                  max(0.0_dp, min(1.0_dp, candidate)))
                if (.not. ieee_is_finite(theta_cell(i, j, k))) then
                  call reset_transport_flux_outputs_3d( &
                    flux_x, flux_y, flux_z, minimum_theta)
                  return
                end if
              end if
            end do
          end do
      end do
      end do
      if (.not. all(ieee_is_finite(theta_cell))) then
        call reset_transport_flux_outputs_3d( &
          flux_x, flux_y, flux_z, minimum_theta)
        return
      end if

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            next_i = periodic_index(i + 1, nx)
            theta_face = min(theta_cell(i, j, k), &
              theta_cell(next_i, j, k))
            if (.not. ieee_is_finite(theta_face)) then
              call reset_transport_flux_outputs_3d( &
                flux_x, flux_y, flux_z, minimum_theta)
              return
            end if
            minimum_theta = min(minimum_theta, theta_face)
            call limit_species_face( &
              flux_x(:, i, j, k), species_energy_x(i, j, k), &
              size(species), theta_face, local_ok)
            if (.not. local_ok) then
              call reset_transport_flux_outputs_3d( &
                flux_x, flux_y, flux_z, minimum_theta)
              return
            end if

            next_j = periodic_index(j + 1, ny)
            theta_face = min(theta_cell(i, j, k), &
              theta_cell(i, next_j, k))
            if (.not. ieee_is_finite(theta_face)) then
              call reset_transport_flux_outputs_3d( &
                flux_x, flux_y, flux_z, minimum_theta)
              return
            end if
            minimum_theta = min(minimum_theta, theta_face)
            call limit_species_face( &
              flux_y(:, i, j, k), species_energy_y(i, j, k), &
              size(species), theta_face, local_ok)
            if (.not. local_ok) then
              call reset_transport_flux_outputs_3d( &
                flux_x, flux_y, flux_z, minimum_theta)
              return
            end if

            next_k = periodic_index(k + 1, nz)
            theta_face = min(theta_cell(i, j, k), &
              theta_cell(i, j, next_k))
            if (.not. ieee_is_finite(theta_face)) then
              call reset_transport_flux_outputs_3d( &
                flux_x, flux_y, flux_z, minimum_theta)
              return
            end if
            minimum_theta = min(minimum_theta, theta_face)
            call limit_species_face( &
              flux_z(:, i, j, k), species_energy_z(i, j, k), &
              size(species), theta_face, local_ok)
            if (.not. local_ok) then
              call reset_transport_flux_outputs_3d( &
                flux_x, flux_y, flux_z, minimum_theta)
              return
            end if
          end do
        end do
      end do
    end if

    if (.not. all(ieee_is_finite(flux_x)) .or. &
        .not. all(ieee_is_finite(flux_y)) .or. &
        .not. all(ieee_is_finite(flux_z)) .or. &
        .not. ieee_is_finite(minimum_theta)) then
      call reset_transport_flux_outputs_3d( &
        flux_x, flux_y, flux_z, minimum_theta)
      return
    end if
    if (minimum_theta < 0.0_dp .or. minimum_theta > 1.0_dp) then
      call reset_transport_flux_outputs_3d( &
        flux_x, flux_y, flux_z, minimum_theta)
      return
    end if
    ok = .true.
  end subroutine reactive_transport_fluxes_3d

  pure subroutine ghosted_centered_velocity_derivative_3d( &
      primitive, i, j, k, velocity_component, direction, dx, dy, dz, &
      derivative, ok)
    real(dp), intent(in) :: primitive(:, 0:, 0:, 0:)
    integer, intent(in) :: i, j, k, velocity_component, direction
    real(dp), intent(in) :: dx, dy, dz
    real(dp), intent(out) :: derivative
    logical, intent(out) :: ok

    real(dp) :: numerator, denominator
    integer :: primitive_component

    derivative = 0.0_dp
    ok = .false.
    if (velocity_component < 1 .or. velocity_component > 3) return
    if (direction < 1 .or. direction > 3) return
    primitive_component = velocity_component + 1
    select case (direction)
    case (1)
      call checked_difference( &
        primitive(primitive_component, i + 1, j, k), &
        primitive(primitive_component, i - 1, j, k), numerator, ok)
      if (.not. ok) return
      call checked_product(2.0_dp, dx, denominator, ok)
      if (.not. ok) return
      call checked_quotient(numerator, denominator, derivative, ok)
      if (.not. ok) return
      derivative = (primitive(primitive_component, i + 1, j, k) - &
        primitive(primitive_component, i - 1, j, k)) / (2.0_dp * dx)
    case (2)
      call checked_difference( &
        primitive(primitive_component, i, j + 1, k), &
        primitive(primitive_component, i, j - 1, k), numerator, ok)
      if (.not. ok) return
      call checked_product(2.0_dp, dy, denominator, ok)
      if (.not. ok) return
      call checked_quotient(numerator, denominator, derivative, ok)
      if (.not. ok) return
      derivative = (primitive(primitive_component, i, j + 1, k) - &
        primitive(primitive_component, i, j - 1, k)) / (2.0_dp * dy)
    case (3)
      call checked_difference( &
        primitive(primitive_component, i, j, k + 1), &
        primitive(primitive_component, i, j, k - 1), numerator, ok)
      if (.not. ok) return
      call checked_product(2.0_dp, dz, denominator, ok)
      if (.not. ok) return
      call checked_quotient(numerator, denominator, derivative, ok)
      if (.not. ok) return
      derivative = (primitive(primitive_component, i, j, k + 1) - &
        primitive(primitive_component, i, j, k - 1)) / (2.0_dp * dz)
    case default
      return
    end select
  end subroutine ghosted_centered_velocity_derivative_3d

  pure subroutine ghosted_face_velocity_gradient_3d( &
      primitive, left_i, left_j, left_k, right_i, right_j, right_k, &
      normal, dx, dy, dz, gradient, ok)
    real(dp), intent(in) :: primitive(:, 0:, 0:, 0:)
    integer, intent(in) :: left_i, left_j, left_k
    integer, intent(in) :: right_i, right_j, right_k, normal
    real(dp), intent(in) :: dx, dy, dz
    real(dp), intent(out) :: gradient(3, 3)
    logical, intent(out) :: ok

    real(dp) :: normal_difference, left_derivative, right_derivative
    real(dp) :: spacing(3)
    logical :: local_ok
    integer :: velocity_component, direction

    spacing = [dx, dy, dz]
    gradient = 0.0_dp
    ok = .false.
    if (normal < 1 .or. normal > 3) return
    do velocity_component = 1, 3
      call checked_difference( &
        primitive(velocity_component + 1, right_i, right_j, right_k), &
        primitive(velocity_component + 1, left_i, left_j, left_k), &
        normal_difference, local_ok)
      if (.not. local_ok) return
      call checked_quotient( &
        normal_difference, spacing(normal), &
        gradient(velocity_component, normal), local_ok)
      if (.not. local_ok) return
      gradient(velocity_component, normal) = &
        (primitive(velocity_component + 1, right_i, right_j, right_k) - &
         primitive(velocity_component + 1, left_i, left_j, left_k)) / &
        spacing(normal)
      do direction = 1, 3
        if (direction == normal) cycle
        call ghosted_centered_velocity_derivative_3d( &
          primitive, left_i, left_j, left_k, velocity_component, &
          direction, dx, dy, dz, left_derivative, local_ok)
        if (.not. local_ok) return
        call ghosted_centered_velocity_derivative_3d( &
          primitive, right_i, right_j, right_k, velocity_component, &
          direction, dx, dy, dz, right_derivative, local_ok)
        if (.not. local_ok) return
        call checked_midpoint( &
          left_derivative, right_derivative, &
          gradient(velocity_component, direction), local_ok)
        if (.not. local_ok) return
        gradient(velocity_component, direction) = &
          0.5_dp * (left_derivative + right_derivative)
      end do
    end do
    ok = .true.
  end subroutine ghosted_face_velocity_gradient_3d

  subroutine reactive_transport_ghosted_fluxes_3d( &
      species, transport, state, temperature, nx, ny, nz, dx, dy, dz, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, &
      flux_x, flux_y, flux_z, minimum_theta, ok, theta_cell_output, &
      exterior_theta_x_lower, exterior_theta_x_upper, &
      periodic_theta_y, periodic_theta_z)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, 0:, 0:, 0:), temperature(0:, 0:, 0:)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz, dt
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: flux_x(:, 0:, :, :)
    real(dp), intent(out) :: flux_y(:, :, 0:, :)
    real(dp), intent(out) :: flux_z(:, :, :, 0:)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    real(dp), intent(out), optional :: theta_cell_output(:, :, :)
    real(dp), intent(in), optional :: exterior_theta_x_lower(:, :)
    real(dp), intent(in), optional :: exterior_theta_x_upper(:, :)
    logical, intent(in), optional :: periodic_theta_y, periodic_theta_z

    real(dp), allocatable :: primitive(:, :, :, :)
    real(dp), allocatable :: checked_temperature(:, :, :)
    real(dp), allocatable :: species_energy_x(:, :, :)
    real(dp), allocatable :: species_energy_y(:, :, :)
    real(dp), allocatable :: species_energy_z(:, :, :)
    real(dp), allocatable :: theta_cell(:, :, :)
    real(dp) :: gradient(3, 3), outgoing, mass, candidate, theta_face
    real(dp) :: sound_speed
    real(dp) :: theta_numerator, theta_denominator
    logical :: local_ok
    logical :: use_periodic_theta_y, use_periodic_theta_z
    integer :: i, j, k, component, species_index, nvar, nprim

    flux_x = 0.0_dp
    flux_y = 0.0_dp
    flux_z = 0.0_dp
    minimum_theta = 1.0_dp
    ok = .false.
    if (present(theta_cell_output)) theta_cell_output = 1.0_dp
    use_periodic_theta_y = .false.
    use_periodic_theta_z = .false.
    if (present(periodic_theta_y)) use_periodic_theta_y = periodic_theta_y
    if (present(periodic_theta_z)) use_periodic_theta_z = periodic_theta_z
    nvar = reactive_nvar(size(species))
    nprim = reactive_nprim(size(species))
    if (.not. all(ieee_is_finite([dx, dy, dz, dt]))) return
    if (size(transport) /= size(species) .or. &
        nx < 1 .or. ny < 1 .or. nz < 1) return
    if (dx <= 0.0_dp .or. dy <= 0.0_dp .or. dz <= 0.0_dp .or. dt < 0.0_dp) return
    if (barodiffusion_enabled .and. .not. species_diffusion_enabled) return
    if (size(state, 1) /= nvar .or. size(state, 2) /= nx + 2 .or. &
        size(state, 3) /= ny + 2 .or. size(state, 4) /= nz + 2 .or. &
        size(temperature, 1) /= nx + 2 .or. &
        size(temperature, 2) /= ny + 2 .or. &
        size(temperature, 3) /= nz + 2) return
    if (size(flux_x, 1) /= nvar .or. size(flux_x, 2) /= nx + 1 .or. &
        size(flux_x, 3) /= ny .or. size(flux_x, 4) /= nz .or. &
        size(flux_y, 1) /= nvar .or. size(flux_y, 2) /= nx .or. &
        size(flux_y, 3) /= ny + 1 .or. size(flux_y, 4) /= nz .or. &
        size(flux_z, 1) /= nvar .or. size(flux_z, 2) /= nx .or. &
        size(flux_z, 3) /= ny .or. size(flux_z, 4) /= nz + 1) return
    if (.not. all(ieee_is_finite(state)) .or. &
        .not. all(ieee_is_finite(temperature))) return
    if (present(theta_cell_output)) then
      if (any(shape(theta_cell_output) /= [nx, ny, nz])) return
    end if
    if (present(exterior_theta_x_lower)) then
      if (any(shape(exterior_theta_x_lower) /= [ny, nz])) return
      if (.not. all(ieee_is_finite(exterior_theta_x_lower))) return
      if (minval(exterior_theta_x_lower) < 0.0_dp .or. &
          maxval(exterior_theta_x_lower) > 1.0_dp) return
    end if
    if (present(exterior_theta_x_upper)) then
      if (any(shape(exterior_theta_x_upper) /= [ny, nz])) return
      if (.not. all(ieee_is_finite(exterior_theta_x_upper))) return
      if (minval(exterior_theta_x_upper) < 0.0_dp .or. &
          maxval(exterior_theta_x_upper) > 1.0_dp) return
    end if
    if (.not. (viscosity_enabled .or. thermal_conduction_enabled .or. &
        species_diffusion_enabled)) then
      ok = .true.
      return
    end if

    allocate(primitive(nprim, 0:nx + 1, 0:ny + 1, 0:nz + 1))
    allocate(checked_temperature(0:nx + 1, 0:ny + 1, 0:nz + 1))
    do k = 0, nz + 1
      do j = 0, ny + 1
        do i = 0, nx + 1
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), &
            primitive(:, i, j, k), checked_temperature(i, j, k), &
            sound_speed, local_ok)
          if (.not. local_ok) then
            call reset_transport_ghosted_outputs_3d( &
              flux_x, flux_y, flux_z, minimum_theta, theta_cell_output)
            return
          end if
        end do
      end do
    end do

    allocate(species_energy_x(0:nx, ny, nz))
    allocate(species_energy_y(nx, 0:ny, nz))
    allocate(species_energy_z(nx, ny, 0:nz))
    species_energy_x = 0.0_dp
    species_energy_y = 0.0_dp
    species_energy_z = 0.0_dp
    do k = 1, nz
      do j = 1, ny
        do i = 0, nx
          call ghosted_face_velocity_gradient_3d( &
            primitive, i, j, k, i + 1, j, k, 1, dx, dy, dz, gradient, local_ok)
          if (.not. local_ok) then
            call reset_transport_ghosted_outputs_3d( &
              flux_x, flux_y, flux_z, minimum_theta, theta_cell_output)
            return
          end if
          call assemble_periodic_transport_face( &
            species, transport, primitive(:, i, j, k), &
            primitive(:, i + 1, j, k), checked_temperature(i, j, k), &
            checked_temperature(i + 1, j, k), dx, 1, gradient, &
            viscosity_enabled, thermal_conduction_enabled, &
            species_diffusion_enabled, barodiffusion_enabled, &
            flux_x(:, i, j, k), species_energy_x(i, j, k), local_ok)
          if (.not. local_ok) then
            call reset_transport_ghosted_outputs_3d( &
              flux_x, flux_y, flux_z, minimum_theta, theta_cell_output)
            return
          end if
        end do
      end do
    end do
    do k = 1, nz
      do j = 0, ny
        do i = 1, nx
          call ghosted_face_velocity_gradient_3d( &
            primitive, i, j, k, i, j + 1, k, 2, dx, dy, dz, gradient, local_ok)
          if (.not. local_ok) then
            call reset_transport_ghosted_outputs_3d( &
              flux_x, flux_y, flux_z, minimum_theta, theta_cell_output)
            return
          end if
          call assemble_periodic_transport_face( &
            species, transport, primitive(:, i, j, k), &
            primitive(:, i, j + 1, k), checked_temperature(i, j, k), &
            checked_temperature(i, j + 1, k), dy, 2, gradient, &
            viscosity_enabled, thermal_conduction_enabled, &
            species_diffusion_enabled, barodiffusion_enabled, &
            flux_y(:, i, j, k), species_energy_y(i, j, k), local_ok)
          if (.not. local_ok) then
            call reset_transport_ghosted_outputs_3d( &
              flux_x, flux_y, flux_z, minimum_theta, theta_cell_output)
            return
          end if
        end do
      end do
    end do
    do k = 0, nz
      do j = 1, ny
        do i = 1, nx
          call ghosted_face_velocity_gradient_3d( &
            primitive, i, j, k, i, j, k + 1, 3, dx, dy, dz, gradient, local_ok)
          if (.not. local_ok) then
            call reset_transport_ghosted_outputs_3d( &
              flux_x, flux_y, flux_z, minimum_theta, theta_cell_output)
            return
          end if
          call assemble_periodic_transport_face( &
            species, transport, primitive(:, i, j, k), &
            primitive(:, i, j, k + 1), checked_temperature(i, j, k), &
            checked_temperature(i, j, k + 1), dz, 3, gradient, &
            viscosity_enabled, thermal_conduction_enabled, &
            species_diffusion_enabled, barodiffusion_enabled, &
            flux_z(:, i, j, k), species_energy_z(i, j, k), local_ok)
          if (.not. local_ok) then
            call reset_transport_ghosted_outputs_3d( &
              flux_x, flux_y, flux_z, minimum_theta, theta_cell_output)
            return
          end if
        end do
      end do
    end do

    if (.not. all(ieee_is_finite(flux_x)) .or. &
        .not. all(ieee_is_finite(flux_y)) .or. &
        .not. all(ieee_is_finite(flux_z)) .or. &
        .not. all(ieee_is_finite(species_energy_x)) .or. &
        .not. all(ieee_is_finite(species_energy_y)) .or. &
        .not. all(ieee_is_finite(species_energy_z))) then
      call reset_transport_ghosted_outputs_3d( &
        flux_x, flux_y, flux_z, minimum_theta, theta_cell_output)
      return
    end if
    if (species_diffusion_enabled .and. dt > 0.0_dp) then
      allocate(theta_cell(nx, ny, nz))
      theta_cell = 1.0_dp
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            do species_index = 1, size(species)
              component = reactive_species_component(species_index)
              call checked_species_outgoing( &
                flux_x(component, i, j, k), &
                flux_x(component, i - 1, j, k), &
                flux_y(component, i, j, k), &
                flux_y(component, i, j - 1, k), &
                flux_z(component, i, j, k), &
                flux_z(component, i, j, k - 1), dx, dy, dz, &
                outgoing, local_ok)
              if (.not. local_ok) then
                call reset_transport_ghosted_outputs_3d( &
                  flux_x, flux_y, flux_z, minimum_theta, theta_cell_output)
                return
              end if
              mass = max(0.0_dp, state(component, i, j, k))
              if (.not. ieee_is_finite(mass)) then
                call reset_transport_ghosted_outputs_3d( &
                  flux_x, flux_y, flux_z, minimum_theta, theta_cell_output)
                return
              end if
              if (outgoing > 0.0_dp) then
                call checked_product( &
                  species_safety, mass, theta_numerator, local_ok)
                if (.not. local_ok) then
                  call reset_transport_ghosted_outputs_3d( &
                    flux_x, flux_y, flux_z, minimum_theta, theta_cell_output)
                  return
                end if
                call checked_product(dt, outgoing, theta_denominator, local_ok)
                if (.not. local_ok) then
                  call reset_transport_ghosted_outputs_3d( &
                    flux_x, flux_y, flux_z, minimum_theta, theta_cell_output)
                  return
                end if
                call checked_quotient( &
                  theta_numerator, theta_denominator, candidate, local_ok)
                if (.not. local_ok) then
                  call reset_transport_ghosted_outputs_3d( &
                    flux_x, flux_y, flux_z, minimum_theta, theta_cell_output)
                  return
                end if
                theta_cell(i, j, k) = min(theta_cell(i, j, k), &
                  max(0.0_dp, min(1.0_dp, candidate)))
                if (.not. ieee_is_finite(theta_cell(i, j, k))) then
                  call reset_transport_ghosted_outputs_3d( &
                    flux_x, flux_y, flux_z, minimum_theta, theta_cell_output)
                  return
                end if
              end if
            end do
          end do
      end do
      end do
      if (.not. all(ieee_is_finite(theta_cell))) then
        call reset_transport_ghosted_outputs_3d( &
          flux_x, flux_y, flux_z, minimum_theta, theta_cell_output)
        return
      end if
      if (present(theta_cell_output)) theta_cell_output = theta_cell
      do k = 1, nz
        do j = 1, ny
          do i = 0, nx
            if (i == 0) then
              theta_face = theta_cell(1, j, k)
              if (present(exterior_theta_x_lower)) &
                theta_face = min( &
                  theta_face, exterior_theta_x_lower(j, k))
            else if (i == nx) then
              theta_face = theta_cell(nx, j, k)
              if (present(exterior_theta_x_upper)) &
                theta_face = min( &
                  theta_face, exterior_theta_x_upper(j, k))
            else
              theta_face = min(theta_cell(i, j, k), &
                theta_cell(i + 1, j, k))
            end if
            minimum_theta = min(minimum_theta, theta_face)
            call limit_species_face( &
              flux_x(:, i, j, k), species_energy_x(i, j, k), &
              size(species), theta_face, local_ok)
            if (.not. local_ok) then
              call reset_transport_ghosted_outputs_3d( &
                flux_x, flux_y, flux_z, minimum_theta, theta_cell_output)
              return
            end if
          end do
        end do
      end do
      do k = 1, nz
        do j = 0, ny
          do i = 1, nx
            if (j == 0) then
              theta_face = theta_cell(i, 1, k)
              if (use_periodic_theta_y) &
                theta_face = min(theta_face, theta_cell(i, ny, k))
            else if (j == ny) then
              theta_face = theta_cell(i, ny, k)
              if (use_periodic_theta_y) &
                theta_face = min(theta_face, theta_cell(i, 1, k))
            else
              theta_face = min(theta_cell(i, j, k), &
                theta_cell(i, j + 1, k))
            end if
            minimum_theta = min(minimum_theta, theta_face)
            call limit_species_face( &
              flux_y(:, i, j, k), species_energy_y(i, j, k), &
              size(species), theta_face, local_ok)
            if (.not. local_ok) then
              call reset_transport_ghosted_outputs_3d( &
                flux_x, flux_y, flux_z, minimum_theta, theta_cell_output)
              return
            end if
          end do
        end do
      end do
      do k = 0, nz
        do j = 1, ny
          do i = 1, nx
            if (k == 0) then
              theta_face = theta_cell(i, j, 1)
              if (use_periodic_theta_z) &
                theta_face = min(theta_face, theta_cell(i, j, nz))
            else if (k == nz) then
              theta_face = theta_cell(i, j, nz)
              if (use_periodic_theta_z) &
                theta_face = min(theta_face, theta_cell(i, j, 1))
            else
              theta_face = min(theta_cell(i, j, k), &
                theta_cell(i, j, k + 1))
            end if
            minimum_theta = min(minimum_theta, theta_face)
            call limit_species_face( &
              flux_z(:, i, j, k), species_energy_z(i, j, k), &
              size(species), theta_face, local_ok)
            if (.not. local_ok) then
              call reset_transport_ghosted_outputs_3d( &
                flux_x, flux_y, flux_z, minimum_theta, theta_cell_output)
              return
            end if
          end do
        end do
      end do
    end if

    if (.not. all(ieee_is_finite(flux_x)) .or. &
        .not. all(ieee_is_finite(flux_y)) .or. &
        .not. all(ieee_is_finite(flux_z)) .or. &
        .not. ieee_is_finite(minimum_theta)) then
      call reset_transport_ghosted_outputs_3d( &
        flux_x, flux_y, flux_z, minimum_theta, theta_cell_output)
      return
    end if
    if (minimum_theta < 0.0_dp .or. minimum_theta > 1.0_dp) then
      call reset_transport_ghosted_outputs_3d( &
        flux_x, flux_y, flux_z, minimum_theta, theta_cell_output)
      return
    end if
    ok = .true.
  end subroutine reactive_transport_ghosted_fluxes_3d

  pure subroutine limit_species_face( &
      flux, species_energy, nspecies, theta, ok)
    real(dp), intent(inout) :: flux(:)
    real(dp), intent(in) :: species_energy, theta
    integer, intent(in) :: nspecies
    logical, intent(out) :: ok

    real(dp) :: correction, energy_correction
    real(dp) :: candidate_flux(size(flux))
    logical :: local_ok
    integer :: species_index, component

    ok = .false.
    candidate_flux = 0.0_dp
    if (nspecies < 1 .or. size(flux) < 1) return
    if (iet < 1 .or. iet > size(flux)) return
    if (.not. ieee_is_finite(theta) .or. &
        .not. ieee_is_finite(species_energy) .or. &
        .not. all(ieee_is_finite(flux))) return
    if (theta < 0.0_dp .or. theta > 1.0_dp) return
    candidate_flux = flux
    do species_index = 1, nspecies
      component = reactive_species_component(species_index)
      if (component < 1 .or. component > size(flux)) return
      call checked_product( &
        theta, flux(component), candidate_flux(component), local_ok)
      if (.not. local_ok) return
    end do
    call checked_difference(theta, 1.0_dp, correction, local_ok)
    if (.not. local_ok) return
    call checked_product( &
      correction, species_energy, energy_correction, local_ok)
    if (.not. local_ok) return
    call checked_sum( &
      flux(iet), energy_correction, candidate_flux(iet), local_ok)
    if (.not. local_ok) return
    if (.not. all(ieee_is_finite(candidate_flux))) return
    flux = candidate_flux
    ok = .true.
  end subroutine limit_species_face

  subroutine reactive_transport_timestep_3d( &
      species, transport, state, temperature, nx, ny, nz, dx, dy, dz, &
      transport_cfl, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, dt, maximum_diffusivity, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz, transport_cfl
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled
    real(dp), intent(out) :: dt, maximum_diffusivity
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:), mass_fractions(:), diffusion(:)
    real(dp) :: checked_temperature, sound_speed, viscosity, conductivity
    real(dp) :: molecular_weight, gas_constant, cp, cv, gamma
    real(dp) :: enthalpy, internal_energy, entropy, candidate, denominator
    real(dp) :: numerator, thermal_denominator
    real(dp) :: inverse_dx, inverse_dy, inverse_dz
    real(dp) :: inverse_dx_squared, inverse_dy_squared, inverse_dz_squared
    real(dp) :: metric_partial, metric_sum
    logical :: local_ok
    integer :: i, j, k, species_index, nspecies, nvar

    dt = 0.0_dp
    maximum_diffusivity = 0.0_dp
    ok = .false.
    nspecies = size(species)
    nvar = reactive_nvar(nspecies)
    if (.not. all(ieee_is_finite([dx, dy, dz, transport_cfl]))) return
    if (size(transport) /= nspecies .or. nx < 1 .or. ny < 1 .or. nz < 1) return
    if (dx <= 0.0_dp .or. dy <= 0.0_dp .or. dz <= 0.0_dp .or. &
        transport_cfl <= 0.0_dp .or. transport_cfl > 0.5_dp) return
    if (.not. valid_transport_shapes_3d( &
          state, temperature, nvar, nx, ny, nz)) return
    if (.not. all(ieee_is_finite(state)) .or. &
        .not. all(ieee_is_finite(temperature))) return
    if (.not. (viscosity_enabled .or. thermal_conduction_enabled .or. &
        species_diffusion_enabled)) then
      dt = huge(1.0_dp)
      ok = .true.
      return
    end if
    allocate(primitive(reactive_nprim(nspecies)))
    allocate(mass_fractions(nspecies), diffusion(nspecies))
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), primitive, &
            checked_temperature, sound_speed, local_ok)
          if (.not. local_ok) then
            call reset_transport_timestep_outputs_3d( &
              dt, maximum_diffusivity)
            return
          end if
          if (.not. all(ieee_is_finite(primitive)) .or. &
              .not. all(ieee_is_finite([checked_temperature, sound_speed]))) then
            call reset_transport_timestep_outputs_3d( &
              dt, maximum_diffusivity)
            return
          end if
          do species_index = 1, nspecies
            mass_fractions(species_index) = &
              primitive(reactive_mass_fraction_component(species_index))
          end do
          call mixture_transport_coefficients( &
            species, transport, mass_fractions, checked_temperature, &
            primitive(5), viscosity, conductivity, diffusion, local_ok)
          if (.not. local_ok) then
            call reset_transport_timestep_outputs_3d( &
              dt, maximum_diffusivity)
            return
          end if
          if (.not. ieee_is_finite(viscosity) .or. &
              .not. ieee_is_finite(conductivity) .or. &
              .not. all(ieee_is_finite(diffusion))) then
            call reset_transport_timestep_outputs_3d( &
              dt, maximum_diffusivity)
            return
          end if
          call mixture_mass_properties( &
            species, mass_fractions, checked_temperature, molecular_weight, &
            gas_constant, cp, cv, gamma, enthalpy, internal_energy, entropy, &
            local_ok)
          if (.not. local_ok) then
            call reset_transport_timestep_outputs_3d( &
              dt, maximum_diffusivity)
            return
          end if
          if (.not. all(ieee_is_finite([molecular_weight, gas_constant, cp, &
              cv, gamma, enthalpy, internal_energy, entropy]))) then
            call reset_transport_timestep_outputs_3d( &
              dt, maximum_diffusivity)
            return
          end if
          if (primitive(1) <= 0.0_dp .or. cv <= 0.0_dp) then
            call reset_transport_timestep_outputs_3d( &
              dt, maximum_diffusivity)
            return
          end if
          if (viscosity_enabled) then
            if (.not. finite_product_fits( &
                4.0_dp / 3.0_dp, viscosity)) then
              call reset_transport_timestep_outputs_3d( &
                dt, maximum_diffusivity)
              return
            end if
            numerator = (4.0_dp / 3.0_dp) * viscosity
            if (.not. finite_quotient_fits( &
                numerator, primitive(1))) then
              call reset_transport_timestep_outputs_3d( &
                dt, maximum_diffusivity)
              return
            end if
            candidate = (4.0_dp / 3.0_dp) * viscosity / primitive(1)
            if (.not. ieee_is_finite(candidate)) then
              call reset_transport_timestep_outputs_3d( &
                dt, maximum_diffusivity)
              return
            end if
            maximum_diffusivity = max(maximum_diffusivity, candidate)
          end if
          if (thermal_conduction_enabled) then
            if (.not. finite_product_fits(primitive(1), cv)) then
              call reset_transport_timestep_outputs_3d( &
                dt, maximum_diffusivity)
              return
            end if
            thermal_denominator = primitive(1) * cv
            if (.not. finite_quotient_fits( &
                conductivity, thermal_denominator)) then
              call reset_transport_timestep_outputs_3d( &
                dt, maximum_diffusivity)
              return
            end if
            candidate = conductivity / (primitive(1) * cv)
            if (.not. ieee_is_finite(candidate)) then
              call reset_transport_timestep_outputs_3d( &
                dt, maximum_diffusivity)
              return
            end if
            maximum_diffusivity = max(maximum_diffusivity, candidate)
          end if
          if (species_diffusion_enabled) then
            candidate = maxval(diffusion)
            if (.not. ieee_is_finite(candidate)) then
              call reset_transport_timestep_outputs_3d( &
                dt, maximum_diffusivity)
              return
            end if
            maximum_diffusivity = max(maximum_diffusivity, candidate)
          end if
        end do
      end do
    end do
    if (.not. ieee_is_finite(maximum_diffusivity)) then
      call reset_transport_timestep_outputs_3d(dt, maximum_diffusivity)
      return
    end if
    if (maximum_diffusivity <= 0.0_dp) then
      dt = huge(1.0_dp)
    else
      call checked_quotient(1.0_dp, dx, inverse_dx, local_ok)
      if (.not. local_ok) then
        call reset_transport_timestep_outputs_3d(dt, maximum_diffusivity)
        return
      end if
      call checked_quotient(1.0_dp, dy, inverse_dy, local_ok)
      if (.not. local_ok) then
        call reset_transport_timestep_outputs_3d(dt, maximum_diffusivity)
        return
      end if
      call checked_quotient(1.0_dp, dz, inverse_dz, local_ok)
      if (.not. local_ok) then
        call reset_transport_timestep_outputs_3d(dt, maximum_diffusivity)
        return
      end if
      if (.not. finite_product_fits(dx, dx) .or. &
          .not. finite_product_fits(dy, dy) .or. &
          .not. finite_product_fits(dz, dz)) then
        call reset_transport_timestep_outputs_3d(dt, maximum_diffusivity)
        return
      end if
      call checked_product( &
        inverse_dx, inverse_dx, inverse_dx_squared, local_ok)
      if (.not. local_ok) then
        call reset_transport_timestep_outputs_3d(dt, maximum_diffusivity)
        return
      end if
      call checked_product( &
        inverse_dy, inverse_dy, inverse_dy_squared, local_ok)
      if (.not. local_ok) then
        call reset_transport_timestep_outputs_3d(dt, maximum_diffusivity)
        return
      end if
      call checked_product( &
        inverse_dz, inverse_dz, inverse_dz_squared, local_ok)
      if (.not. local_ok) then
        call reset_transport_timestep_outputs_3d(dt, maximum_diffusivity)
        return
      end if
      call checked_sum( &
        inverse_dx_squared, inverse_dy_squared, metric_partial, local_ok)
      if (.not. local_ok) then
        call reset_transport_timestep_outputs_3d(dt, maximum_diffusivity)
        return
      end if
      call checked_sum( &
        metric_partial, inverse_dz_squared, metric_sum, local_ok)
      if (.not. local_ok) then
        call reset_transport_timestep_outputs_3d(dt, maximum_diffusivity)
        return
      end if
      call checked_product( &
        maximum_diffusivity, metric_sum, denominator, local_ok)
      if (.not. local_ok) then
        call reset_transport_timestep_outputs_3d(dt, maximum_diffusivity)
        return
      end if
      if (.not. finite_quotient_fits(transport_cfl, denominator)) then
        call reset_transport_timestep_outputs_3d(dt, maximum_diffusivity)
        return
      end if
      dt = transport_cfl / denominator
    end if
    if (.not. ieee_is_finite(dt) .or. &
        .not. ieee_is_finite(maximum_diffusivity)) then
      call reset_transport_timestep_outputs_3d(dt, maximum_diffusivity)
      return
    end if
    if (dt <= 0.0_dp) then
      call reset_transport_timestep_outputs_3d(dt, maximum_diffusivity)
      return
    end if
    ok = .true.
  end subroutine reactive_transport_timestep_3d

  subroutine reactive_transport_euler_update_3d( &
      species, transport, input_state, input_temperature, nx, ny, nz, &
      dx, dy, dz, dt, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, output_state, &
      output_temperature, minimum_theta, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: input_state(:, :, :, :)
    real(dp), intent(in) :: input_temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz, dt
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: output_state(:, :, :, :)
    real(dp), intent(out) :: output_temperature(:, :, :)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok

    real(dp), allocatable :: flux_x(:, :, :, :)
    real(dp), allocatable :: flux_y(:, :, :, :)
    real(dp), allocatable :: flux_z(:, :, :, :), primitive(:), candidate_cell(:)
    real(dp) :: checked_temperature, sound_speed
    real(dp) :: coefficient_x, coefficient_y, coefficient_z
    real(dp) :: difference, term, updated
    logical :: local_ok
    integer :: i, j, k, previous_i, previous_j, previous_k, component, nvar

    output_state = 0.0_dp
    output_temperature = 0.0_dp
    minimum_theta = 1.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. all(ieee_is_finite([dx, dy, dz, dt]))) return
    if (barodiffusion_enabled .and. .not. species_diffusion_enabled) return
    if (.not. valid_transport_shapes_3d( &
          input_state, input_temperature, nvar, nx, ny, nz)) return
    if (any(shape(output_state) /= shape(input_state)) .or. &
        any(shape(output_temperature) /= shape(input_temperature))) return
    if (.not. all(ieee_is_finite(input_state)) .or. &
        .not. all(ieee_is_finite(input_temperature))) return
    allocate(flux_x(nvar, nx, ny, nz), flux_y(nvar, nx, ny, nz))
    allocate(flux_z(nvar, nx, ny, nz))
    allocate(primitive(reactive_nprim(size(species))))
    allocate(candidate_cell(nvar))
    call reactive_transport_fluxes_3d( &
      species, transport, input_state, input_temperature, nx, ny, nz, &
      dx, dy, dz, dt, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, &
      flux_x, flux_y, flux_z, minimum_theta, local_ok)
    if (.not. local_ok) then
      call reset_transport_euler_outputs_3d( &
        output_state, output_temperature, minimum_theta)
      return
    end if

    call checked_quotient(dt, dx, coefficient_x, local_ok)
    if (.not. local_ok) then
      call reset_transport_euler_outputs_3d( &
        output_state, output_temperature, minimum_theta)
      return
    end if
    call checked_quotient(dt, dy, coefficient_y, local_ok)
    if (.not. local_ok) then
      call reset_transport_euler_outputs_3d( &
        output_state, output_temperature, minimum_theta)
      return
    end if
    call checked_quotient(dt, dz, coefficient_z, local_ok)
    if (.not. local_ok) then
      call reset_transport_euler_outputs_3d( &
        output_state, output_temperature, minimum_theta)
      return
    end if

    do k = 1, nz
      previous_k = periodic_index(k - 1, nz)
      do j = 1, ny
        previous_j = periodic_index(j - 1, ny)
        do i = 1, nx
          previous_i = periodic_index(i - 1, nx)
          candidate_cell = input_state(:, i, j, k)
          do component = 1, nvar
            call checked_difference( &
              flux_x(component, i, j, k), &
              flux_x(component, previous_i, j, k), difference, local_ok)
            if (.not. local_ok) then
              call reset_transport_euler_outputs_3d( &
                output_state, output_temperature, minimum_theta)
              return
            end if
            call checked_product(coefficient_x, difference, term, local_ok)
            if (.not. local_ok) then
              call reset_transport_euler_outputs_3d( &
                output_state, output_temperature, minimum_theta)
              return
            end if
            call checked_difference(candidate_cell(component), term, updated, local_ok)
            if (.not. local_ok) then
              call reset_transport_euler_outputs_3d( &
                output_state, output_temperature, minimum_theta)
              return
            end if
            candidate_cell(component) = updated

            call checked_difference( &
              flux_y(component, i, j, k), &
              flux_y(component, i, previous_j, k), difference, local_ok)
            if (.not. local_ok) then
              call reset_transport_euler_outputs_3d( &
                output_state, output_temperature, minimum_theta)
              return
            end if
            call checked_product(coefficient_y, difference, term, local_ok)
            if (.not. local_ok) then
              call reset_transport_euler_outputs_3d( &
                output_state, output_temperature, minimum_theta)
              return
            end if
            call checked_difference(candidate_cell(component), term, updated, local_ok)
            if (.not. local_ok) then
              call reset_transport_euler_outputs_3d( &
                output_state, output_temperature, minimum_theta)
              return
            end if
            candidate_cell(component) = updated

            call checked_difference( &
              flux_z(component, i, j, k), &
              flux_z(component, i, j, previous_k), difference, local_ok)
            if (.not. local_ok) then
              call reset_transport_euler_outputs_3d( &
                output_state, output_temperature, minimum_theta)
              return
            end if
            call checked_product(coefficient_z, difference, term, local_ok)
            if (.not. local_ok) then
              call reset_transport_euler_outputs_3d( &
                output_state, output_temperature, minimum_theta)
              return
            end if
            call checked_difference(candidate_cell(component), term, updated, local_ok)
            if (.not. local_ok) then
              call reset_transport_euler_outputs_3d( &
                output_state, output_temperature, minimum_theta)
              return
            end if
            candidate_cell(component) = updated
          end do
          candidate_cell = input_state(:, i, j, k) - &
            dt / dx * (flux_x(:, i, j, k) - &
              flux_x(:, previous_i, j, k)) - &
            dt / dy * (flux_y(:, i, j, k) - &
              flux_y(:, i, previous_j, k)) - &
            dt / dz * (flux_z(:, i, j, k) - &
              flux_z(:, i, j, previous_k))
          if (.not. all(ieee_is_finite(candidate_cell))) then
            call reset_transport_euler_outputs_3d( &
              output_state, output_temperature, minimum_theta)
            return
          end if
          call reactive_conserved_to_primitive( &
            species, candidate_cell, input_temperature(i, j, k), &
            primitive, checked_temperature, sound_speed, local_ok)
          if (.not. local_ok) then
            call reset_transport_euler_outputs_3d( &
              output_state, output_temperature, minimum_theta)
            return
          end if
          if (.not. all(ieee_is_finite(primitive)) .or. &
              .not. all(ieee_is_finite([checked_temperature, sound_speed]))) then
            call reset_transport_euler_outputs_3d( &
              output_state, output_temperature, minimum_theta)
            return
          end if
          output_state(:, i, j, k) = candidate_cell
          output_temperature(i, j, k) = checked_temperature
        end do
      end do
    end do
    if (.not. all(ieee_is_finite(output_state)) .or. &
        .not. all(ieee_is_finite(output_temperature)) .or. &
        .not. ieee_is_finite(minimum_theta)) then
      call reset_transport_euler_outputs_3d( &
        output_state, output_temperature, minimum_theta)
      return
    end if
    if (minimum_theta < 0.0_dp .or. minimum_theta > 1.0_dp) then
      call reset_transport_euler_outputs_3d( &
        output_state, output_temperature, minimum_theta)
      return
    end if
    ok = .true.
  end subroutine reactive_transport_euler_update_3d

  subroutine advance_reactive_transport_3d( &
      species, transport, state, temperature, nx, ny, nz, dx, dy, dz, &
      interval, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, minimum_theta, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: dx, dy, dz, interval
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
    real(dp) :: theta1, theta2, checked_temperature, sound_speed
    real(dp) :: temperature_guess, midpoint_value
    logical :: local_ok
    integer :: i, j, k, component, nvar

    minimum_theta = 1.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. all(ieee_is_finite([dx, dy, dz, interval]))) return
    if (dx <= 0.0_dp .or. dy <= 0.0_dp .or. dz <= 0.0_dp) return
    if (barodiffusion_enabled .and. .not. species_diffusion_enabled) return
    if (.not. valid_transport_shapes_3d( &
          state, temperature, nvar, nx, ny, nz)) return
    if (.not. all(ieee_is_finite(state)) .or. &
        .not. all(ieee_is_finite(temperature))) return
    if (interval < 0.0_dp) return
    if (interval <= tiny(1.0_dp) .or. .not. (viscosity_enabled .or. &
        thermal_conduction_enabled .or. species_diffusion_enabled)) then
      ok = .true.
      return
    end if
    allocate(initial_state, source=state)
    allocate(initial_temperature, source=temperature)
    allocate(stage_state(nvar, nx, ny, nz))
    allocate(stage_temperature(nx, ny, nz))
    allocate(euler_state(nvar, nx, ny, nz))
    allocate(euler_temperature(nx, ny, nz))
    allocate(candidate_state(nvar, nx, ny, nz))
    allocate(candidate_temperature(nx, ny, nz))
    allocate(primitive(reactive_nprim(size(species))))

    call reactive_transport_euler_update_3d( &
      species, transport, initial_state, initial_temperature, nx, ny, nz, &
      dx, dy, dz, interval, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, stage_state, &
      stage_temperature, theta1, local_ok)
    if (.not. local_ok) return
    call reactive_transport_euler_update_3d( &
      species, transport, stage_state, stage_temperature, nx, ny, nz, &
      dx, dy, dz, interval, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, euler_state, &
      euler_temperature, theta2, local_ok)
    if (.not. local_ok) return
    if (.not. all(ieee_is_finite(stage_state)) .or. &
        .not. all(ieee_is_finite(stage_temperature)) .or. &
        .not. all(ieee_is_finite(euler_state)) .or. &
        .not. all(ieee_is_finite(euler_temperature))) return
    if (.not. ieee_is_finite(theta1) .or. &
        .not. ieee_is_finite(theta2)) return

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          do component = 1, nvar
            call checked_midpoint( &
              initial_state(component, i, j, k), &
              euler_state(component, i, j, k), &
              midpoint_value, local_ok)
            if (.not. local_ok) return
          end do
        end do
      end do
    end do
    candidate_state = 0.5_dp * (initial_state + euler_state)
    if (.not. all(ieee_is_finite(candidate_state))) return
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          call checked_midpoint( &
            initial_temperature(i, j, k), euler_temperature(i, j, k), &
            midpoint_value, local_ok)
          if (.not. local_ok) return
          temperature_guess = 0.5_dp * &
            (initial_temperature(i, j, k) + euler_temperature(i, j, k))
          call reactive_conserved_to_primitive( &
            species, candidate_state(:, i, j, k), temperature_guess, &
            primitive, checked_temperature, sound_speed, local_ok)
          if (.not. local_ok) return
          if (.not. all(ieee_is_finite(primitive)) .or. &
              .not. all(ieee_is_finite([checked_temperature, sound_speed]))) return
          candidate_temperature(i, j, k) = checked_temperature
        end do
      end do
    end do
    if (.not. all(ieee_is_finite(candidate_temperature))) return
    minimum_theta = min(theta1, theta2)
    if (.not. ieee_is_finite(minimum_theta)) then
      minimum_theta = 1.0_dp
      return
    end if
    if (minimum_theta < 0.0_dp .or. minimum_theta > 1.0_dp) then
      minimum_theta = 1.0_dp
      return
    end if
    state = candidate_state
    temperature = candidate_temperature
    ok = .true.
  end subroutine advance_reactive_transport_3d

end module reactive_transport_3d_mod
