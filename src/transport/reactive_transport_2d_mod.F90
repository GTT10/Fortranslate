module reactive_transport_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species, nasa7_mass_properties
  use mixture_thermo_mod, only: &
    mole_fractions_from_mass_fractions, mixture_mass_properties
  use gas_transport_mod, only: gas_transport_species
  use mixture_transport_mod, only: mixture_transport_coefficients
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_mass_fraction_component, reactive_conserved_to_primitive
  use reactive_boundary_2d_mod, only: &
    reactive_boundary_face_2d, reactive_boundary_set_2d, &
    initialize_periodic_boundary_set_2d, validate_reactive_boundary_set_2d, &
    sample_reactive_primitive_2d, reactive_boundary_is_periodic, &
    reactive_boundary_is_wall, reactive_boundary_is_inflow, &
    reactive_boundary_has_prescribed_species_flux
  implicit none
  private

  real(dp), parameter :: species_safety = 0.90_dp

  type, public :: reactive_transport_exterior_2d
    real(dp), allocatable :: x_lower_primitive(:, :)
    real(dp), allocatable :: x_upper_primitive(:, :)
    real(dp), allocatable :: y_lower_primitive(:, :)
    real(dp), allocatable :: y_upper_primitive(:, :)
    real(dp), allocatable :: x_lower_temperature(:)
    real(dp), allocatable :: x_upper_temperature(:)
    real(dp), allocatable :: y_lower_temperature(:)
    real(dp), allocatable :: y_upper_temperature(:)
  contains
    procedure :: is_valid => reactive_transport_exterior_is_valid
  end type reactive_transport_exterior_2d

  public :: reactive_transport_timestep_2d
  public :: reactive_transport_fluxes_2d_faces
  public :: reactive_transport_fluxes_2d
  public :: reactive_transport_euler_update_2d
  public :: advance_reactive_transport_2d
  public :: face_transport_data
  public :: species_face_flux

contains

  pure logical function reactive_transport_exterior_is_valid( &
      self, nprimitive, nx, ny) result(valid)
    class(reactive_transport_exterior_2d), intent(in) :: self
    integer, intent(in) :: nprimitive, nx, ny

    valid = nprimitive >= 1 .and. nx >= 1 .and. ny >= 1 .and. &
      allocated(self%x_lower_primitive) .and. &
      allocated(self%x_upper_primitive) .and. &
      allocated(self%y_lower_primitive) .and. &
      allocated(self%y_upper_primitive) .and. &
      allocated(self%x_lower_temperature) .and. &
      allocated(self%x_upper_temperature) .and. &
      allocated(self%y_lower_temperature) .and. &
      allocated(self%y_upper_temperature)
    if (.not. valid) return
    valid = all(shape(self%x_lower_primitive) == [nprimitive, ny]) .and. &
      all(shape(self%x_upper_primitive) == [nprimitive, ny]) .and. &
      all(shape(self%y_lower_primitive) == [nprimitive, nx]) .and. &
      all(shape(self%y_upper_primitive) == [nprimitive, nx]) .and. &
      size(self%x_lower_temperature) == ny .and. &
      size(self%x_upper_temperature) == ny .and. &
      size(self%y_lower_temperature) == nx .and. &
      size(self%y_upper_temperature) == nx
    if (.not. valid) return
    if (.not. all(ieee_is_finite(self%x_lower_primitive))) then
      valid = .false.
      return
    end if
    if (.not. all(ieee_is_finite(self%x_upper_primitive))) then
      valid = .false.
      return
    end if
    if (.not. all(ieee_is_finite(self%y_lower_primitive))) then
      valid = .false.
      return
    end if
    if (.not. all(ieee_is_finite(self%y_upper_primitive))) then
      valid = .false.
      return
    end if
    if (.not. all(ieee_is_finite(self%x_lower_temperature))) then
      valid = .false.
      return
    end if
    if (.not. all(ieee_is_finite(self%x_upper_temperature))) then
      valid = .false.
      return
    end if
    if (.not. all(ieee_is_finite(self%y_lower_temperature))) then
      valid = .false.
      return
    end if
    if (.not. all(ieee_is_finite(self%y_upper_temperature))) then
      valid = .false.
      return
    end if
    if (any(self%x_lower_temperature <= 0.0_dp)) then
      valid = .false.
      return
    end if
    if (any(self%x_upper_temperature <= 0.0_dp)) then
      valid = .false.
      return
    end if
    if (any(self%y_lower_temperature <= 0.0_dp)) then
      valid = .false.
      return
    end if
    if (any(self%y_upper_temperature <= 0.0_dp)) then
      valid = .false.
      return
    end if
    valid = .true.
  end function reactive_transport_exterior_is_valid

  subroutine sample_transport_primitive_2d( &
      primitive, temperature, nx, ny, i, j, boundaries, sampled, &
      sampled_temperature, ok, exterior)
    real(dp), intent(in) :: primitive(:, :, :), temperature(:, :)
    integer, intent(in) :: nx, ny, i, j
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    real(dp), intent(out) :: sampled(:), sampled_temperature
    logical, intent(out) :: ok
    type(reactive_transport_exterior_2d), intent(in), optional :: exterior

    integer :: mapped_i, mapped_j

    if (.not. present(exterior) .or. &
        (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny)) then
      call sample_reactive_primitive_2d( &
        primitive, temperature, nx, ny, i, j, boundaries, sampled, &
        sampled_temperature, ok)
      return
    end if
    sampled = 0.0_dp
    sampled_temperature = 0.0_dp
    ok = .false.
    if (.not. exterior%is_valid(size(primitive, 1), nx, ny)) return
    mapped_i = min(nx, max(1, i))
    mapped_j = min(ny, max(1, j))
    if (i < 1) then
      sampled = exterior%x_lower_primitive(:, mapped_j)
      sampled_temperature = exterior%x_lower_temperature(mapped_j)
    else if (i > nx) then
      sampled = exterior%x_upper_primitive(:, mapped_j)
      sampled_temperature = exterior%x_upper_temperature(mapped_j)
    else if (j < 1) then
      sampled = exterior%y_lower_primitive(:, mapped_i)
      sampled_temperature = exterior%y_lower_temperature(mapped_i)
    else if (j > ny) then
      sampled = exterior%y_upper_primitive(:, mapped_i)
      sampled_temperature = exterior%y_upper_temperature(mapped_i)
    else
      return
    end if
    if (.not. all(ieee_is_finite(sampled))) return
    if (.not. ieee_is_finite(sampled_temperature)) return
    if (sampled_temperature <= 0.0_dp) return
    ok = .true.
  end subroutine sample_transport_primitive_2d

  pure integer function periodic_index(index, extent) result(wrapped)
    integer, intent(in) :: index, extent
    wrapped = 1 + modulo(index - 1, extent)
  end function periodic_index

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
    fits = abs(numerator) <= huge(1.0_dp) * abs(denominator)
  end function finite_quotient_fits

  pure subroutine add_finite(left, right, result, ok)
    real(dp), intent(in) :: left, right
    real(dp), intent(out) :: result
    logical, intent(out) :: ok

    result = 0.0_dp
    ok = .false.
    if (.not. finite_sum_fits(left, right)) return
    result = left + right
    if (.not. ieee_is_finite(result)) then
      result = 0.0_dp
      return
    end if
    ok = .true.
  end subroutine add_finite

  pure subroutine subtract_finite(left, right, result, ok)
    real(dp), intent(in) :: left, right
    real(dp), intent(out) :: result
    logical, intent(out) :: ok

    call add_finite(left, -right, result, ok)
  end subroutine subtract_finite

  pure subroutine multiply_finite(left, right, result, ok)
    real(dp), intent(in) :: left, right
    real(dp), intent(out) :: result
    logical, intent(out) :: ok

    result = 0.0_dp
    ok = .false.
    if (.not. finite_product_fits(left, right)) return
    result = left * right
    if (.not. ieee_is_finite(result)) then
      result = 0.0_dp
      return
    end if
    ok = .true.
  end subroutine multiply_finite

  pure subroutine divide_finite(numerator, denominator, result, ok)
    real(dp), intent(in) :: numerator, denominator
    real(dp), intent(out) :: result
    logical, intent(out) :: ok

    result = 0.0_dp
    ok = .false.
    if (.not. finite_quotient_fits(numerator, denominator)) return
    result = numerator / denominator
    if (.not. ieee_is_finite(result)) then
      result = 0.0_dp
      return
    end if
    ok = .true.
  end subroutine divide_finite

  pure subroutine midpoint_finite(left, right, result, ok)
    real(dp), intent(in) :: left, right
    real(dp), intent(out) :: result
    logical, intent(out) :: ok

    result = 0.0_dp
    ok = .false.
    if (.not. ieee_is_finite(left)) return
    if (.not. ieee_is_finite(right)) return
    if (finite_sum_fits(left, right)) then
      result = 0.5_dp * (left + right)
    else
      result = 0.5_dp * left + 0.5_dp * right
    end if
    if (.not. ieee_is_finite(result)) then
      result = 0.0_dp
      return
    end if
    ok = .true.
  end subroutine midpoint_finite

  pure subroutine finite_sum(values, total, ok)
    real(dp), intent(in) :: values(:)
    real(dp), intent(out) :: total
    logical, intent(out) :: ok

    real(dp) :: candidate
    logical :: local_ok
    integer :: i

    total = 0.0_dp
    ok = .false.
    do i = 1, size(values)
      call add_finite(total, values(i), candidate, local_ok)
      if (.not. local_ok) then
        total = 0.0_dp
        return
      end if
      total = candidate
    end do
    ok = .true.
  end subroutine finite_sum

  pure logical function finite_absolute_sum_fits(values) result(fits)
    real(dp), intent(in) :: values(:)

    real(dp) :: bound, next_bound
    integer :: i

    fits = .false.
    bound = 0.0_dp
    do i = 1, size(values)
      if (.not. ieee_is_finite(values(i))) return
      if (.not. finite_sum_fits(bound, abs(values(i)))) return
      next_bound = bound + abs(values(i))
      if (.not. ieee_is_finite(next_bound)) return
      bound = next_bound
    end do
    fits = .true.
  end function finite_absolute_sum_fits

  pure subroutine checked_difference_quotient_2d( &
      left, right, spacing, value, ok)
    real(dp), intent(in) :: left, right, spacing
    real(dp), intent(out) :: value
    logical, intent(out) :: ok

    real(dp) :: difference

    value = 0.0_dp
    ok = .false.
    call subtract_finite(left, right, difference, ok)
    if (.not. ok) return
    call divide_finite(difference, spacing, value, ok)
  end subroutine checked_difference_quotient_2d

  pure subroutine checked_half_sum_2d(left, right, value, ok)
    real(dp), intent(in) :: left, right
    real(dp), intent(out) :: value
    logical, intent(out) :: ok

    real(dp) :: sum

    value = 0.0_dp
    ok = .false.
    call add_finite(left, right, sum, ok)
    if (.not. ok) return
    call multiply_finite(0.5_dp, sum, value, ok)
  end subroutine checked_half_sum_2d

  pure subroutine checked_centered_difference_2d( &
      left, right, spacing, value, ok)
    real(dp), intent(in) :: left, right, spacing
    real(dp), intent(out) :: value
    logical, intent(out) :: ok

    real(dp) :: difference, half_difference, twice_spacing

    value = 0.0_dp
    ok = .false.
    call subtract_finite(left, right, difference, ok)
    if (.not. ok) return
    call multiply_finite(0.5_dp, difference, half_difference, ok)
    if (.not. ok) return
    call multiply_finite(2.0_dp, spacing, twice_spacing, ok)
    if (.not. ok) return
    call divide_finite(half_difference, twice_spacing, value, ok)
  end subroutine checked_centered_difference_2d

  pure subroutine validate_face_stress_arithmetic_2d( &
      viscosity, normal_gradient, cross_gradient_one, cross_gradient_two, &
      transverse_gradient, divergence, uleft, uright, vleft, vright, &
      wleft, wright, ok)
    real(dp), intent(in) :: viscosity, normal_gradient
    real(dp), intent(in) :: cross_gradient_one, cross_gradient_two
    real(dp), intent(in) :: transverse_gradient, divergence
    real(dp), intent(in) :: uleft, uright, vleft, vright, wleft, wright
    logical, intent(out) :: ok

    real(dp) :: twice_normal, divergence_term, normal_term
    real(dp) :: cross_term, tau_normal, tau_cross, tau_transverse
    real(dp) :: velocity_sum, uface, vface, wface
    real(dp) :: energy_normal, energy_cross, energy_transverse
    real(dp) :: energy_sum, energy_total
    logical :: local_ok

    ok = .false.
    call multiply_finite(2.0_dp, normal_gradient, twice_normal, local_ok)
    if (.not. local_ok) return
    call multiply_finite(2.0_dp / 3.0_dp, divergence, divergence_term, local_ok)
    if (.not. local_ok) return
    call subtract_finite(twice_normal, divergence_term, normal_term, local_ok)
    if (.not. local_ok) return
    call multiply_finite(viscosity, normal_term, tau_normal, local_ok)
    if (.not. local_ok) return
    call add_finite(cross_gradient_one, cross_gradient_two, cross_term, local_ok)
    if (.not. local_ok) return
    call multiply_finite(viscosity, cross_term, tau_cross, local_ok)
    if (.not. local_ok) return
    call multiply_finite(viscosity, transverse_gradient, tau_transverse, local_ok)
    if (.not. local_ok) return

    call add_finite(uleft, uright, velocity_sum, local_ok)
    if (.not. local_ok) return
    call multiply_finite(0.5_dp, velocity_sum, uface, local_ok)
    if (.not. local_ok) return
    call add_finite(vleft, vright, velocity_sum, local_ok)
    if (.not. local_ok) return
    call multiply_finite(0.5_dp, velocity_sum, vface, local_ok)
    if (.not. local_ok) return
    call add_finite(wleft, wright, velocity_sum, local_ok)
    if (.not. local_ok) return
    call multiply_finite(0.5_dp, velocity_sum, wface, local_ok)
    if (.not. local_ok) return

    call multiply_finite(tau_normal, uface, energy_normal, local_ok)
    if (.not. local_ok) return
    call multiply_finite(tau_cross, vface, energy_cross, local_ok)
    if (.not. local_ok) return
    call multiply_finite(tau_transverse, wface, energy_transverse, local_ok)
    if (.not. local_ok) return
    call add_finite(energy_normal, energy_cross, energy_sum, local_ok)
    if (.not. local_ok) return
    call add_finite(energy_sum, energy_transverse, energy_total, local_ok)
    ok = local_ok
  end subroutine validate_face_stress_arithmetic_2d

  pure subroutine validate_euler_component_arithmetic_2d( &
      input_value, x_flux_right, x_flux_left, y_flux_right, y_flux_left, &
      dx, dy, dt, ok)
    real(dp), intent(in) :: input_value, x_flux_right, x_flux_left
    real(dp), intent(in) :: y_flux_right, y_flux_left, dx, dy, dt
    logical, intent(out) :: ok

    real(dp) :: x_difference, y_difference, x_rate, y_rate
    real(dp) :: x_update, y_update, candidate, candidate_next
    logical :: local_ok

    ok = .false.
    call subtract_finite(x_flux_right, x_flux_left, x_difference, local_ok)
    if (.not. local_ok) return
    call subtract_finite(y_flux_right, y_flux_left, y_difference, local_ok)
    if (.not. local_ok) return
    call divide_finite(dt, dx, x_rate, local_ok)
    if (.not. local_ok) return
    call divide_finite(dt, dy, y_rate, local_ok)
    if (.not. local_ok) return
    call multiply_finite(x_rate, x_difference, x_update, local_ok)
    if (.not. local_ok) return
    call multiply_finite(y_rate, y_difference, y_update, local_ok)
    if (.not. local_ok) return
    call subtract_finite(input_value, x_update, candidate, local_ok)
    if (.not. local_ok) return
    call subtract_finite(candidate, y_update, candidate_next, local_ok)
    ok = local_ok
  end subroutine validate_euler_component_arithmetic_2d

  pure subroutine validate_timestep_arithmetic_2d( &
      dx, dy, transport_cfl, maximum_diffusivity, ok)
    real(dp), intent(in) :: dx, dy, transport_cfl, maximum_diffusivity
    logical, intent(out) :: ok

    real(dp) :: dx_squared, dy_squared, inverse_dx_squared
    real(dp) :: inverse_dy_squared, metric_sum, denominator, candidate_dt
    logical :: local_ok

    ok = .false.
    call multiply_finite(dx, dx, dx_squared, local_ok)
    if (.not. local_ok) return
    call multiply_finite(dy, dy, dy_squared, local_ok)
    if (.not. local_ok) return
    call divide_finite(1.0_dp, dx_squared, inverse_dx_squared, local_ok)
    if (.not. local_ok) return
    call divide_finite(1.0_dp, dy_squared, inverse_dy_squared, local_ok)
    if (.not. local_ok) return
    call add_finite( &
      inverse_dx_squared, inverse_dy_squared, metric_sum, local_ok)
    if (.not. local_ok) return
    call multiply_finite( &
      maximum_diffusivity, metric_sum, denominator, local_ok)
    if (.not. local_ok) return
    call divide_finite(transport_cfl, denominator, candidate_dt, local_ok)
    ok = local_ok
  end subroutine validate_timestep_arithmetic_2d

  pure subroutine validate_limiter_candidate_2d( &
      flux_x_right, flux_x_left, flux_y_right, flux_y_left, mass, dx, dy, dt, &
      ok)
    real(dp), intent(in) :: flux_x_right, flux_x_left, flux_y_right, flux_y_left
    real(dp), intent(in) :: mass, dx, dy, dt
    logical, intent(out) :: ok

    real(dp) :: outgoing_x_right, outgoing_x_left
    real(dp) :: outgoing_y_right, outgoing_y_left, outgoing
    real(dp) :: denominator, numerator, candidate
    real(dp) :: partial_outgoing, next_partial_outgoing
    logical :: local_ok

    ok = .false.
    call divide_finite(max(flux_x_right, 0.0_dp), dx, &
      outgoing_x_right, local_ok)
    if (.not. local_ok) return
    call divide_finite(max(-flux_x_left, 0.0_dp), dx, &
      outgoing_x_left, local_ok)
    if (.not. local_ok) return
    call divide_finite(max(flux_y_right, 0.0_dp), dy, &
      outgoing_y_right, local_ok)
    if (.not. local_ok) return
    call divide_finite(max(-flux_y_left, 0.0_dp), dy, &
      outgoing_y_left, local_ok)
    if (.not. local_ok) return
    call add_finite(outgoing_x_right, outgoing_x_left, &
      partial_outgoing, local_ok)
    if (.not. local_ok) return
    call add_finite(partial_outgoing, outgoing_y_right, &
      next_partial_outgoing, local_ok)
    if (.not. local_ok) return
    partial_outgoing = next_partial_outgoing
    call add_finite(partial_outgoing, outgoing_y_left, outgoing, local_ok)
    if (.not. local_ok) return
    if (outgoing <= 0.0_dp) then
      ok = .true.
      return
    end if
    call multiply_finite(dt, outgoing, denominator, local_ok)
    if (.not. local_ok) return
    call multiply_finite(species_safety, mass, numerator, local_ok)
    if (.not. local_ok) return
    call divide_finite(numerator, denominator, candidate, local_ok)
    ok = local_ok
  end subroutine validate_limiter_candidate_2d

  subroutine recover_primitives_2d( &
      species, state, temperature, nx, ny, primitive, checked_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    integer, intent(in) :: nx, ny
    real(dp), intent(out) :: primitive(:, :, :)
    real(dp), intent(out) :: checked_temperature(:, :)
    logical, intent(out) :: ok

    real(dp) :: sound_speed
    logical :: local_ok
    integer :: i, j

    ok = .false.
    if (size(state, 1) /= reactive_nvar(size(species)) .or. &
        size(state, 2) < nx .or. size(state, 3) < ny .or. &
        size(temperature, 1) < nx .or. size(temperature, 2) < ny .or. &
        size(primitive, 1) /= reactive_nprim(size(species)) .or. &
        size(primitive, 2) < nx .or. size(primitive, 3) < ny .or. &
        size(checked_temperature, 1) < nx .or. &
        size(checked_temperature, 2) < ny) return
    do j = 1, ny
      do i = 1, nx
        call reactive_conserved_to_primitive( &
          species, state(:, i, j), temperature(i, j), primitive(:, i, j), &
          checked_temperature(i, j), sound_speed, local_ok)
        if (.not. local_ok) return
      end do
    end do
    ok = .true.
  end subroutine recover_primitives_2d

  subroutine face_transport_data( &
      species, transport, left_primitive, right_primitive, left_temperature, &
      right_temperature, viscosity, conductivity, diffusion, yleft, yright, &
      yface, xleft, xright, hface, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: left_primitive(:), right_primitive(:)
    real(dp), intent(in) :: left_temperature, right_temperature
    real(dp), intent(out) :: viscosity, conductivity, diffusion(:)
    real(dp), intent(out) :: yleft(:), yright(:), yface(:)
    real(dp), intent(out) :: xleft(:), xright(:), hface(:)
    logical, intent(out) :: ok

    real(dp) :: tface, pface, cp, cv, hleft, hright, eint, entropy
    real(dp) :: candidate_viscosity, candidate_conductivity
    real(dp) :: midpoint, composition_total
    real(dp), allocatable :: candidate_diffusion(:), candidate_yleft(:)
    real(dp), allocatable :: candidate_yright(:), candidate_yface(:)
    real(dp), allocatable :: candidate_xleft(:), candidate_xright(:)
    real(dp), allocatable :: candidate_hface(:)
    logical :: local_ok
    integer :: k, nspecies

    viscosity = 0.0_dp
    conductivity = 0.0_dp
    diffusion = 0.0_dp
    yleft = 0.0_dp
    yright = 0.0_dp
    yface = 0.0_dp
    xleft = 0.0_dp
    xright = 0.0_dp
    hface = 0.0_dp
    ok = .false.
    nspecies = size(species)
    if (nspecies < 1) return
    if (size(transport) /= nspecies .or. &
        size(left_primitive) /= reactive_nprim(nspecies) .or. &
        size(right_primitive) /= reactive_nprim(nspecies) .or. &
        size(diffusion) /= nspecies .or. size(yleft) /= nspecies .or. &
        size(yright) /= nspecies .or. size(yface) /= nspecies .or. &
        size(xleft) /= nspecies .or. size(xright) /= nspecies .or. &
        size(hface) /= nspecies) return

    if (.not. all(ieee_is_finite(left_primitive))) return
    if (.not. all(ieee_is_finite(right_primitive))) return
    if (.not. ieee_is_finite(left_temperature)) return
    if (.not. ieee_is_finite(right_temperature)) return
    if (left_temperature <= 0.0_dp) return
    if (right_temperature <= 0.0_dp) return
    if (left_primitive(1) <= 0.0_dp) return
    if (right_primitive(1) <= 0.0_dp) return
    if (left_primitive(5) <= 0.0_dp) return
    if (right_primitive(5) <= 0.0_dp) return

    allocate(candidate_diffusion(nspecies), candidate_yleft(nspecies))
    allocate(candidate_yright(nspecies), candidate_yface(nspecies))
    allocate(candidate_xleft(nspecies), candidate_xright(nspecies))
    allocate(candidate_hface(nspecies))
    candidate_diffusion = 0.0_dp
    candidate_yleft = 0.0_dp
    candidate_yright = 0.0_dp
    candidate_yface = 0.0_dp
    candidate_xleft = 0.0_dp
    candidate_xright = 0.0_dp
    candidate_hface = 0.0_dp
    candidate_viscosity = 0.0_dp
    candidate_conductivity = 0.0_dp

    do k = 1, nspecies
      candidate_yleft(k) = left_primitive(reactive_mass_fraction_component(k))
      candidate_yright(k) = right_primitive(reactive_mass_fraction_component(k))
    end do
    if (.not. all(ieee_is_finite(candidate_yleft))) return
    if (.not. all(ieee_is_finite(candidate_yright))) return
    do k = 1, nspecies
      call midpoint_finite( &
        candidate_yleft(k), candidate_yright(k), midpoint, local_ok)
      if (.not. local_ok) return
      candidate_yface(k) = max(0.0_dp, midpoint)
    end do
    call finite_sum(candidate_yface, composition_total, local_ok)
    if (.not. local_ok) return
    if (composition_total <= tiny(1.0_dp)) return
    do k = 1, nspecies
      call divide_finite( &
        candidate_yface(k), composition_total, midpoint, local_ok)
      if (.not. local_ok) return
      candidate_yface(k) = midpoint
    end do
    call midpoint_finite(left_temperature, right_temperature, tface, local_ok)
    if (.not. local_ok) return
    call midpoint_finite(left_primitive(5), right_primitive(5), pface, local_ok)
    if (.not. local_ok) return
    call mixture_transport_coefficients( &
      species, transport, candidate_yface, tface, pface, candidate_viscosity, &
      candidate_conductivity, candidate_diffusion, local_ok)
    if (.not. local_ok) return
    call mole_fractions_from_mass_fractions( &
      species, candidate_yleft, candidate_xleft, local_ok)
    if (.not. local_ok) return
    call mole_fractions_from_mass_fractions( &
      species, candidate_yright, candidate_xright, local_ok)
    if (.not. local_ok) return
    if (.not. all(ieee_is_finite(candidate_xleft))) return
    if (.not. all(ieee_is_finite(candidate_xright))) return
    do k = 1, nspecies
      call nasa7_mass_properties( &
        species(k), left_temperature, cp, cv, hleft, eint, entropy, local_ok)
      if (.not. local_ok) return
      call nasa7_mass_properties( &
        species(k), right_temperature, cp, cv, hright, eint, entropy, local_ok)
      if (.not. local_ok) return
      if (.not. ieee_is_finite(hleft)) return
      if (.not. ieee_is_finite(hright)) return
      call midpoint_finite(hleft, hright, candidate_hface(k), local_ok)
      if (.not. local_ok) return
    end do
    if (.not. ieee_is_finite(candidate_viscosity)) return
    if (candidate_viscosity <= 0.0_dp) return
    if (.not. ieee_is_finite(candidate_conductivity)) return
    if (candidate_conductivity <= 0.0_dp) return
    if (.not. all(ieee_is_finite(candidate_diffusion))) return
    if (any(candidate_diffusion < 0.0_dp)) return
    if (.not. all(ieee_is_finite(candidate_hface))) return
    viscosity = candidate_viscosity
    conductivity = candidate_conductivity
    diffusion = candidate_diffusion
    yleft = candidate_yleft
    yright = candidate_yright
    yface = candidate_yface
    xleft = candidate_xleft
    xright = candidate_xright
    hface = candidate_hface
    ok = .true.
  end subroutine face_transport_data

  subroutine species_face_flux( &
      species, diffusion, yleft, yright, yface, xleft, xright, hface, &
      density_face, pressure_left, pressure_right, pressure_face, spacing, &
      barodiffusion_enabled, species_flux, enthalpy_flux, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: diffusion(:), yleft(:), yright(:), yface(:)
    real(dp), intent(in) :: xleft(:), xright(:), hface(:)
    real(dp), intent(in) :: density_face, pressure_left, pressure_right
    real(dp), intent(in) :: pressure_face, spacing
    logical, intent(in) :: barodiffusion_enabled
    real(dp), intent(out) :: species_flux(:), enthalpy_flux
    logical, intent(out) :: ok

    real(dp), allocatable :: raw_flux(:), candidate_flux(:)
    real(dp) :: dlnp, correction, enthalpy_candidate, enthalpy_next
    real(dp) :: pressure_difference, denominator, midpoint
    real(dp) :: mole_difference, gradient, composition_difference
    real(dp) :: barodiffusion_term, driving_term, density_diffusion
    real(dp) :: product, closure, scale, tolerance, tolerance_factor
    logical :: local_ok
    integer :: k, nspecies

    species_flux = 0.0_dp
    enthalpy_flux = 0.0_dp
    ok = .false.
    nspecies = size(species)
    if (nspecies < 1) return
    if (size(diffusion) /= nspecies .or. size(yleft) /= nspecies .or. &
        size(yright) /= nspecies .or. size(yface) /= nspecies .or. &
        size(xleft) /= nspecies .or. size(xright) /= nspecies .or. &
        size(hface) /= nspecies .or. size(species_flux) /= nspecies) return
    if (.not. ieee_is_finite(density_face)) return
    if (.not. ieee_is_finite(pressure_left)) return
    if (.not. ieee_is_finite(pressure_right)) return
    if (.not. ieee_is_finite(pressure_face)) return
    if (.not. ieee_is_finite(spacing)) return
    if (density_face <= 0.0_dp) return
    if (pressure_left <= 0.0_dp) return
    if (pressure_right <= 0.0_dp) return
    if (pressure_face <= 0.0_dp) return
    if (spacing <= 0.0_dp) return
    if (.not. all(ieee_is_finite(diffusion))) return
    if (any(diffusion < 0.0_dp)) return
    if (.not. all(ieee_is_finite(yleft))) return
    if (.not. all(ieee_is_finite(yright))) return
    if (.not. all(ieee_is_finite(yface))) return
    if (.not. all(ieee_is_finite(xleft))) return
    if (.not. all(ieee_is_finite(xright))) return
    if (.not. all(ieee_is_finite(hface))) return
    allocate(raw_flux(nspecies))
    allocate(candidate_flux(nspecies))
    raw_flux = 0.0_dp
    candidate_flux = 0.0_dp
    dlnp = 0.0_dp
    if (barodiffusion_enabled) then
      call subtract_finite( &
        pressure_right, pressure_left, pressure_difference, local_ok)
      if (.not. local_ok) return
      call multiply_finite(spacing, pressure_face, denominator, local_ok)
      if (.not. local_ok) return
      call divide_finite(pressure_difference, denominator, dlnp, local_ok)
      if (.not. local_ok) return
    end if
    do k = 1, nspecies
      call subtract_finite( &
        xright(k), xleft(k), mole_difference, local_ok)
      if (.not. local_ok) return
      call divide_finite(mole_difference, spacing, gradient, local_ok)
      if (.not. local_ok) return
      if (.not. finite_sum_fits(xleft(k), xright(k))) return
      call midpoint_finite(xleft(k), xright(k), midpoint, local_ok)
      if (.not. local_ok) return
      call subtract_finite( &
        midpoint, yface(k), composition_difference, local_ok)
      if (.not. local_ok) return
      call multiply_finite( &
        composition_difference, dlnp, barodiffusion_term, local_ok)
      if (.not. local_ok) return
      call add_finite(gradient, barodiffusion_term, driving_term, local_ok)
      if (.not. local_ok) return
      call multiply_finite( &
        density_face, diffusion(k), density_diffusion, local_ok)
      if (.not. local_ok) return
      call multiply_finite( &
        density_diffusion, driving_term, raw_flux(k), local_ok)
      if (.not. local_ok) return
      raw_flux(k) = -raw_flux(k)
    end do
    if (.not. finite_absolute_sum_fits(raw_flux)) return
    call finite_sum(raw_flux, correction, local_ok)
    if (.not. local_ok) return
    do k = 1, nspecies
      call multiply_finite(yface(k), correction, product, local_ok)
      if (.not. local_ok) return
      call subtract_finite(raw_flux(k), product, candidate_flux(k), local_ok)
      if (.not. local_ok) return
    end do
    if (nspecies > 1) then
      call finite_sum( &
        candidate_flux(1:nspecies - 1), correction, local_ok)
      if (.not. local_ok) return
      candidate_flux(nspecies) = -correction
    else
      candidate_flux(1) = 0.0_dp
    end if
    if (.not. all(ieee_is_finite(candidate_flux))) return
    enthalpy_candidate = 0.0_dp
    do k = 1, nspecies
      call multiply_finite( &
        hface(k), candidate_flux(k), product, local_ok)
      if (.not. local_ok) return
      call add_finite(enthalpy_candidate, product, enthalpy_next, local_ok)
      if (.not. local_ok) return
      enthalpy_candidate = enthalpy_next
    end do

    ! The checks above are a fail-closed preflight.  Keep the established
    ! normal-input operation order so Debug and Release checkpoint contracts
    ! retain their independently frozen byte streams.
    dlnp = 0.0_dp
    if (barodiffusion_enabled) then
      dlnp = (pressure_right - pressure_left) / (spacing * pressure_face)
    end if
    do k = 1, nspecies
      raw_flux(k) = -density_face * diffusion(k) * &
        ((xright(k) - xleft(k)) / spacing + &
          (0.5_dp * (xleft(k) + xright(k)) - yface(k)) * dlnp)
    end do
    if (.not. all(ieee_is_finite(raw_flux))) return
    if (.not. finite_absolute_sum_fits(raw_flux)) return
    correction = sum(raw_flux)
    candidate_flux = raw_flux - yface * correction
    if (.not. all(ieee_is_finite(candidate_flux))) return
    if (nspecies > 1) then
      if (.not. finite_absolute_sum_fits( &
          candidate_flux(1:nspecies - 1))) return
      candidate_flux(nspecies) = -sum(candidate_flux(1:nspecies - 1))
    else
      candidate_flux(1) = 0.0_dp
    end if
    if (.not. all(ieee_is_finite(candidate_flux))) return
    do k = 1, nspecies
      if (.not. finite_product_fits(hface(k), candidate_flux(k))) return
      raw_flux(k) = hface(k) * candidate_flux(k)
    end do
    if (.not. finite_absolute_sum_fits(raw_flux)) return
    enthalpy_candidate = sum(hface * candidate_flux)
    if (.not. ieee_is_finite(enthalpy_candidate)) return
    call finite_sum(candidate_flux, closure, local_ok)
    if (.not. local_ok) return
    scale = max(1.0_dp, maxval(abs(candidate_flux)))
    tolerance_factor = 2.0e3_dp * epsilon(1.0_dp)
    if (.not. ieee_is_finite(tolerance_factor)) return
    if (tolerance_factor <= 1.0_dp) then
      tolerance = tolerance_factor * scale
      if (.not. ieee_is_finite(tolerance)) return
      if (abs(closure) > tolerance) return
    else if (scale <= huge(1.0_dp) / tolerance_factor) then
      tolerance = tolerance_factor * scale
      if (.not. ieee_is_finite(tolerance)) return
      if (abs(closure) > tolerance) return
    end if
    species_flux = candidate_flux
    enthalpy_flux = enthalpy_candidate
    ok = .true.
  end subroutine species_face_flux

  subroutine prescribed_wall_species_flux( &
      face, coordinate_sign, hface, species_flux, enthalpy_flux, ok)
    type(reactive_boundary_face_2d), intent(in) :: face
    real(dp), intent(in) :: coordinate_sign, hface(:)
    real(dp), intent(out) :: species_flux(:), enthalpy_flux
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate_flux(:)
    real(dp) :: scale, prescribed_sum, closure, tolerance, tolerance_factor
    real(dp) :: enthalpy_candidate, enthalpy_next, product
    logical :: local_ok
    integer :: k, nspecies

    species_flux = 0.0_dp
    enthalpy_flux = 0.0_dp
    ok = .false.
    nspecies = size(species_flux)
    if (nspecies < 1) return
    if (.not. reactive_boundary_has_prescribed_species_flux(face)) return
    if (.not. allocated(face%prescribed_species_flux)) return
    if (size(face%prescribed_species_flux) /= size(species_flux) .or. &
        size(hface) /= size(species_flux)) return
    if (.not. ieee_is_finite(coordinate_sign)) return
    if (abs(coordinate_sign) < 1.0_dp .or. &
        abs(coordinate_sign) > 1.0_dp) return
    if (.not. all(ieee_is_finite(face%prescribed_species_flux))) return
    if (.not. all(ieee_is_finite(hface))) return
    call finite_sum(face%prescribed_species_flux, prescribed_sum, local_ok)
    if (.not. local_ok) return
    scale = max(1.0_dp, maxval(abs(face%prescribed_species_flux)))
    tolerance_factor = 2.0e3_dp * epsilon(1.0_dp)
    if (.not. ieee_is_finite(tolerance_factor)) return
    if (tolerance_factor <= 1.0_dp) then
      tolerance = tolerance_factor * scale
      if (.not. ieee_is_finite(tolerance)) return
      if (abs(prescribed_sum) > tolerance) return
    else if (scale <= huge(1.0_dp) / tolerance_factor) then
      tolerance = tolerance_factor * scale
      if (.not. ieee_is_finite(tolerance)) return
      if (abs(prescribed_sum) > tolerance) return
    end if
    allocate(candidate_flux(nspecies))
    candidate_flux = 0.0_dp
    do k = 1, nspecies
      call multiply_finite( &
        coordinate_sign, face%prescribed_species_flux(k), &
        candidate_flux(k), local_ok)
      if (.not. local_ok) return
    end do
    if (nspecies > 1) then
      call finite_sum(candidate_flux(1:nspecies - 1), closure, local_ok)
      if (.not. local_ok) return
      candidate_flux(nspecies) = -closure
    else
      candidate_flux(1) = 0.0_dp
    end if
    enthalpy_candidate = 0.0_dp
    do k = 1, nspecies
      call multiply_finite(hface(k), candidate_flux(k), product, local_ok)
      if (.not. local_ok) return
      call add_finite(enthalpy_candidate, product, enthalpy_next, local_ok)
      if (.not. local_ok) return
      enthalpy_candidate = enthalpy_next
    end do
    if (.not. all(ieee_is_finite(candidate_flux))) return
    if (.not. ieee_is_finite(enthalpy_candidate)) return
    species_flux = candidate_flux
    enthalpy_flux = enthalpy_candidate
    ok = .true.
  end subroutine prescribed_wall_species_flux

  subroutine transport_face_flux_x( &
      species, transport, primitive, checked_temperature, nx, ny, face_i, j, &
      dx, dy, boundaries, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, flux, &
      species_energy, ok, exterior)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: primitive(:, :, :), checked_temperature(:, :)
    integer, intent(in) :: nx, ny, face_i, j
    real(dp), intent(in) :: dx, dy
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: flux(:), species_energy
    logical, intent(out) :: ok
    type(reactive_transport_exterior_2d), intent(in), optional :: exterior

    real(dp), allocatable :: qleft(:), qright(:), qtmp(:), qtmp2(:)
    real(dp), allocatable :: yleft(:), yright(:), yface(:)
    real(dp), allocatable :: xleft(:), xright(:), diffusion(:), hface(:)
    real(dp), allocatable :: species_flux(:)
    real(dp) :: tleft, tright, ttmp, ttmp2, spacing, coordinate_sign
    real(dp) :: viscosity, conductivity, density_face, pressure_face
    real(dp) :: dudx, dvdx, dwdx, dudy, dvdy, divu
    real(dp) :: tau_xx, tau_xy, tau_xz, dtdx, uface, vface, wface
    real(dp) :: guard_value, guard_value_two, guard_sum, centered_value
    logical :: local_ok, wall_face, slip_face
    integer :: nspecies, nprim, k, side

    flux = 0.0_dp
    species_energy = 0.0_dp
    ok = .false.
    nspecies = size(species)
    nprim = reactive_nprim(nspecies)
    if (face_i < 0 .or. face_i > nx .or. j < 1 .or. j > ny .or. &
        size(flux) /= reactive_nvar(nspecies)) return
    allocate(qleft(nprim), qright(nprim), qtmp(nprim), qtmp2(nprim))
    allocate(yleft(nspecies), yright(nspecies), yface(nspecies))
    allocate(xleft(nspecies), xright(nspecies), diffusion(nspecies))
    allocate(hface(nspecies), species_flux(nspecies))

    call sample_transport_primitive_2d( &
      primitive, checked_temperature, nx, ny, face_i, j, boundaries, &
      qleft, tleft, local_ok, exterior)
    if (.not. local_ok) return
    call sample_transport_primitive_2d( &
      primitive, checked_temperature, nx, ny, face_i + 1, j, boundaries, &
      qright, tright, local_ok, exterior)
    if (.not. local_ok) return

    spacing = dx
    wall_face = .false.
    slip_face = .false.
    side = 0
    if (face_i == 0 .and. .not. &
        reactive_boundary_is_periodic(boundaries%face(1))) side = 1
    if (face_i == nx .and. .not. &
        reactive_boundary_is_periodic(boundaries%face(2))) side = 2
    if (side > 0) then
      wall_face = reactive_boundary_is_wall(boundaries%face(side))
      slip_face = trim(boundaries%face(side)%kind) == 'slip_wall'
      if (reactive_boundary_is_inflow(boundaries%face(side))) spacing = 0.5_dp * dx
    end if

    call face_transport_data( &
      species, transport, qleft, qright, tleft, tright, viscosity, &
      conductivity, diffusion, yleft, yright, yface, xleft, xright, hface, &
      local_ok)
    if (.not. local_ok) return
    call checked_half_sum_2d(qleft(1), qright(1), guard_value, local_ok)
    if (.not. local_ok) return
    density_face = 0.5_dp * (qleft(1) + qright(1))
    call checked_half_sum_2d(qleft(5), qright(5), guard_value, local_ok)
    if (.not. local_ok) return
    pressure_face = 0.5_dp * (qleft(5) + qright(5))
    call checked_difference_quotient_2d( &
      qright(2), qleft(2), spacing, guard_value, local_ok)
    if (.not. local_ok) return
    dudx = (qright(2) - qleft(2)) / spacing
    call checked_difference_quotient_2d( &
      qright(3), qleft(3), spacing, guard_value, local_ok)
    if (.not. local_ok) return
    dvdx = (qright(3) - qleft(3)) / spacing
    call checked_difference_quotient_2d( &
      qright(4), qleft(4), spacing, guard_value, local_ok)
    if (.not. local_ok) return
    dwdx = (qright(4) - qleft(4)) / spacing

    call sample_transport_primitive_2d( &
      primitive, checked_temperature, nx, ny, face_i, j + 1, boundaries, &
      qtmp, ttmp, local_ok, exterior)
    if (.not. local_ok) return
    call sample_transport_primitive_2d( &
      primitive, checked_temperature, nx, ny, face_i, j - 1, boundaries, &
      qtmp2, ttmp2, local_ok, exterior)
    if (.not. local_ok) return
    call checked_centered_difference_2d( &
      qtmp(2), qtmp2(2), dy, centered_value, local_ok)
    if (.not. local_ok) return
    dudy = 0.5_dp * (qtmp(2) - qtmp2(2)) / (2.0_dp * dy)
    call checked_centered_difference_2d( &
      qtmp(3), qtmp2(3), dy, centered_value, local_ok)
    if (.not. local_ok) return
    dvdy = 0.5_dp * (qtmp(3) - qtmp2(3)) / (2.0_dp * dy)
    call sample_transport_primitive_2d( &
      primitive, checked_temperature, nx, ny, face_i + 1, j + 1, &
      boundaries, qtmp, ttmp, local_ok, exterior)
    if (.not. local_ok) return
    call sample_transport_primitive_2d( &
      primitive, checked_temperature, nx, ny, face_i + 1, j - 1, &
      boundaries, qtmp2, ttmp2, local_ok, exterior)
    if (.not. local_ok) return
    call checked_centered_difference_2d( &
      qtmp(2), qtmp2(2), dy, centered_value, local_ok)
    if (.not. local_ok) return
    call add_finite(dudy, centered_value, guard_sum, local_ok)
    if (.not. local_ok) return
    dudy = dudy + 0.5_dp * (qtmp(2) - qtmp2(2)) / (2.0_dp * dy)
    call checked_centered_difference_2d( &
      qtmp(3), qtmp2(3), dy, centered_value, local_ok)
    if (.not. local_ok) return
    call add_finite(dvdy, centered_value, guard_sum, local_ok)
    if (.not. local_ok) return
    dvdy = dvdy + 0.5_dp * (qtmp(3) - qtmp2(3)) / (2.0_dp * dy)
    call add_finite(dudx, dvdy, guard_value, local_ok)
    if (.not. local_ok) return
    divu = dudx + dvdy

    if (viscosity_enabled) then
      call validate_face_stress_arithmetic_2d( &
        viscosity, dudx, dudy, dvdx, dwdx, divu, qleft(2), qright(2), &
        qleft(3), qright(3), qleft(4), qright(4), local_ok)
      if (.not. local_ok) return
      tau_xx = viscosity * (2.0_dp * dudx - (2.0_dp / 3.0_dp) * divu)
      tau_xy = viscosity * (dudy + dvdx)
      tau_xz = viscosity * dwdx
      if (slip_face) then
        tau_xy = 0.0_dp
        tau_xz = 0.0_dp
      end if
      flux(imx) = -tau_xx
      flux(imy) = -tau_xy
      flux(imz) = -tau_xz
      uface = 0.5_dp * (qleft(2) + qright(2))
      vface = 0.5_dp * (qleft(3) + qright(3))
      wface = 0.5_dp * (qleft(4) + qright(4))
      flux(iet) = -(tau_xx * uface + tau_xy * vface + tau_xz * wface)
    end if
    if (thermal_conduction_enabled) then
      call checked_difference_quotient_2d( &
        tright, tleft, spacing, guard_value, local_ok)
      if (.not. local_ok) return
      call multiply_finite(conductivity, guard_value, guard_value_two, local_ok)
      if (.not. local_ok) return
      call add_finite(flux(iet), -guard_value_two, guard_sum, local_ok)
      if (.not. local_ok) return
      dtdx = (tright - tleft) / spacing
      flux(iet) = flux(iet) - conductivity * dtdx
    end if
    if (species_diffusion_enabled) then
      if (wall_face) then
        species_flux = 0.0_dp
        species_energy = 0.0_dp
        local_ok = .true.
        if (reactive_boundary_has_prescribed_species_flux( &
            boundaries%face(side))) then
          coordinate_sign = 1.0_dp
          if (side == 2) coordinate_sign = -1.0_dp
          call prescribed_wall_species_flux( &
            boundaries%face(side), coordinate_sign, hface, species_flux, &
            species_energy, local_ok)
        end if
      else
        call species_face_flux( &
          species, diffusion, yleft, yright, yface, xleft, xright, hface, &
          density_face, qleft(5), qright(5), pressure_face, spacing, &
          barodiffusion_enabled, species_flux, species_energy, local_ok)
      end if
      if (.not. local_ok) return
      call add_finite(flux(iet), species_energy, guard_sum, local_ok)
      if (.not. local_ok) return
      do k = 1, nspecies
        flux(reactive_species_component(k)) = species_flux(k)
      end do
      flux(iet) = flux(iet) + species_energy
    end if
    flux(irho) = 0.0_dp
    ok = all(ieee_is_finite(flux)) .and. ieee_is_finite(species_energy)
  end subroutine transport_face_flux_x

  subroutine transport_face_flux_y( &
      species, transport, primitive, checked_temperature, nx, ny, i, face_j, &
      dx, dy, boundaries, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, flux, &
      species_energy, ok, exterior)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: primitive(:, :, :), checked_temperature(:, :)
    integer, intent(in) :: nx, ny, i, face_j
    real(dp), intent(in) :: dx, dy
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: flux(:), species_energy
    logical, intent(out) :: ok
    type(reactive_transport_exterior_2d), intent(in), optional :: exterior

    real(dp), allocatable :: qlower(:), qupper(:), qtmp(:), qtmp2(:)
    real(dp), allocatable :: ylower(:), yupper(:), yface(:)
    real(dp), allocatable :: xlower(:), xupper(:), diffusion(:), hface(:)
    real(dp), allocatable :: species_flux(:)
    real(dp) :: tlower, tupper, ttmp, ttmp2, spacing, coordinate_sign
    real(dp) :: viscosity, conductivity, density_face, pressure_face
    real(dp) :: dudy, dvdy, dwdy, dudx, dvdx, divu
    real(dp) :: tau_yx, tau_yy, tau_yz, dtdy, uface, vface, wface
    real(dp) :: guard_value, guard_value_two, guard_sum, centered_value
    logical :: local_ok, wall_face, slip_face
    integer :: nspecies, nprim, k, side

    flux = 0.0_dp
    species_energy = 0.0_dp
    ok = .false.
    nspecies = size(species)
    nprim = reactive_nprim(nspecies)
    if (i < 1 .or. i > nx .or. face_j < 0 .or. face_j > ny .or. &
        size(flux) /= reactive_nvar(nspecies)) return
    allocate(qlower(nprim), qupper(nprim), qtmp(nprim), qtmp2(nprim))
    allocate(ylower(nspecies), yupper(nspecies), yface(nspecies))
    allocate(xlower(nspecies), xupper(nspecies), diffusion(nspecies))
    allocate(hface(nspecies), species_flux(nspecies))

    call sample_transport_primitive_2d( &
      primitive, checked_temperature, nx, ny, i, face_j, boundaries, &
      qlower, tlower, local_ok, exterior)
    if (.not. local_ok) return
    call sample_transport_primitive_2d( &
      primitive, checked_temperature, nx, ny, i, face_j + 1, boundaries, &
      qupper, tupper, local_ok, exterior)
    if (.not. local_ok) return

    spacing = dy
    wall_face = .false.
    slip_face = .false.
    side = 0
    if (face_j == 0 .and. .not. &
        reactive_boundary_is_periodic(boundaries%face(3))) side = 3
    if (face_j == ny .and. .not. &
        reactive_boundary_is_periodic(boundaries%face(4))) side = 4
    if (side > 0) then
      wall_face = reactive_boundary_is_wall(boundaries%face(side))
      slip_face = trim(boundaries%face(side)%kind) == 'slip_wall'
      if (reactive_boundary_is_inflow(boundaries%face(side))) spacing = 0.5_dp * dy
    end if

    call face_transport_data( &
      species, transport, qlower, qupper, tlower, tupper, viscosity, &
      conductivity, diffusion, ylower, yupper, yface, xlower, xupper, hface, &
      local_ok)
    if (.not. local_ok) return
    call checked_half_sum_2d(qlower(1), qupper(1), guard_value, local_ok)
    if (.not. local_ok) return
    density_face = 0.5_dp * (qlower(1) + qupper(1))
    call checked_half_sum_2d(qlower(5), qupper(5), guard_value, local_ok)
    if (.not. local_ok) return
    pressure_face = 0.5_dp * (qlower(5) + qupper(5))
    call checked_difference_quotient_2d( &
      qupper(2), qlower(2), spacing, guard_value, local_ok)
    if (.not. local_ok) return
    dudy = (qupper(2) - qlower(2)) / spacing
    call checked_difference_quotient_2d( &
      qupper(3), qlower(3), spacing, guard_value, local_ok)
    if (.not. local_ok) return
    dvdy = (qupper(3) - qlower(3)) / spacing
    call checked_difference_quotient_2d( &
      qupper(4), qlower(4), spacing, guard_value, local_ok)
    if (.not. local_ok) return
    dwdy = (qupper(4) - qlower(4)) / spacing

    call sample_transport_primitive_2d( &
      primitive, checked_temperature, nx, ny, i + 1, face_j, boundaries, &
      qtmp, ttmp, local_ok, exterior)
    if (.not. local_ok) return
    call sample_transport_primitive_2d( &
      primitive, checked_temperature, nx, ny, i - 1, face_j, boundaries, &
      qtmp2, ttmp2, local_ok, exterior)
    if (.not. local_ok) return
    call checked_centered_difference_2d( &
      qtmp(2), qtmp2(2), dx, centered_value, local_ok)
    if (.not. local_ok) return
    dudx = 0.5_dp * (qtmp(2) - qtmp2(2)) / (2.0_dp * dx)
    call checked_centered_difference_2d( &
      qtmp(3), qtmp2(3), dx, centered_value, local_ok)
    if (.not. local_ok) return
    dvdx = 0.5_dp * (qtmp(3) - qtmp2(3)) / (2.0_dp * dx)
    call sample_transport_primitive_2d( &
      primitive, checked_temperature, nx, ny, i + 1, face_j + 1, &
      boundaries, qtmp, ttmp, local_ok, exterior)
    if (.not. local_ok) return
    call sample_transport_primitive_2d( &
      primitive, checked_temperature, nx, ny, i - 1, face_j + 1, &
      boundaries, qtmp2, ttmp2, local_ok, exterior)
    if (.not. local_ok) return
    call checked_centered_difference_2d( &
      qtmp(2), qtmp2(2), dx, centered_value, local_ok)
    if (.not. local_ok) return
    call add_finite(dudx, centered_value, guard_sum, local_ok)
    if (.not. local_ok) return
    dudx = dudx + 0.5_dp * (qtmp(2) - qtmp2(2)) / (2.0_dp * dx)
    call checked_centered_difference_2d( &
      qtmp(3), qtmp2(3), dx, centered_value, local_ok)
    if (.not. local_ok) return
    call add_finite(dvdx, centered_value, guard_sum, local_ok)
    if (.not. local_ok) return
    dvdx = dvdx + 0.5_dp * (qtmp(3) - qtmp2(3)) / (2.0_dp * dx)
    call add_finite(dudx, dvdy, guard_value, local_ok)
    if (.not. local_ok) return
    divu = dudx + dvdy

    if (viscosity_enabled) then
      call validate_face_stress_arithmetic_2d( &
        viscosity, dvdy, dudy, dvdx, dwdy, divu, qlower(2), qupper(2), &
        qlower(3), qupper(3), qlower(4), qupper(4), local_ok)
      if (.not. local_ok) return
      tau_yx = viscosity * (dudy + dvdx)
      tau_yy = viscosity * (2.0_dp * dvdy - (2.0_dp / 3.0_dp) * divu)
      tau_yz = viscosity * dwdy
      if (slip_face) then
        tau_yx = 0.0_dp
        tau_yz = 0.0_dp
      end if
      flux(imx) = -tau_yx
      flux(imy) = -tau_yy
      flux(imz) = -tau_yz
      uface = 0.5_dp * (qlower(2) + qupper(2))
      vface = 0.5_dp * (qlower(3) + qupper(3))
      wface = 0.5_dp * (qlower(4) + qupper(4))
      flux(iet) = -(tau_yx * uface + tau_yy * vface + tau_yz * wface)
    end if
    if (thermal_conduction_enabled) then
      call checked_difference_quotient_2d( &
        tupper, tlower, spacing, guard_value, local_ok)
      if (.not. local_ok) return
      call multiply_finite(conductivity, guard_value, guard_value_two, local_ok)
      if (.not. local_ok) return
      call add_finite(flux(iet), -guard_value_two, guard_sum, local_ok)
      if (.not. local_ok) return
      dtdy = (tupper - tlower) / spacing
      flux(iet) = flux(iet) - conductivity * dtdy
    end if
    if (species_diffusion_enabled) then
      if (wall_face) then
        species_flux = 0.0_dp
        species_energy = 0.0_dp
        local_ok = .true.
        if (reactive_boundary_has_prescribed_species_flux( &
            boundaries%face(side))) then
          coordinate_sign = 1.0_dp
          if (side == 4) coordinate_sign = -1.0_dp
          call prescribed_wall_species_flux( &
            boundaries%face(side), coordinate_sign, hface, species_flux, &
            species_energy, local_ok)
        end if
      else
        call species_face_flux( &
          species, diffusion, ylower, yupper, yface, xlower, xupper, hface, &
          density_face, qlower(5), qupper(5), pressure_face, spacing, &
          barodiffusion_enabled, species_flux, species_energy, local_ok)
      end if
      if (.not. local_ok) return
      call add_finite(flux(iet), species_energy, guard_sum, local_ok)
      if (.not. local_ok) return
      do k = 1, nspecies
        flux(reactive_species_component(k)) = species_flux(k)
      end do
      flux(iet) = flux(iet) + species_energy
    end if
    flux(irho) = 0.0_dp
    ok = all(ieee_is_finite(flux)) .and. ieee_is_finite(species_energy)
  end subroutine transport_face_flux_y

  subroutine reactive_transport_fluxes_2d_faces( &
      species, transport, state, temperature, nx, ny, dx, dy, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, flux_x, flux_y, &
      minimum_theta, ok, boundaries, exterior)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: dx, dy, dt
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: flux_x(:, 0:, :), flux_y(:, :, 0:)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    type(reactive_boundary_set_2d), intent(in), optional :: boundaries
    type(reactive_transport_exterior_2d), intent(in), optional :: exterior

    type(reactive_boundary_set_2d) :: active_boundaries
    real(dp), allocatable :: primitive(:, :, :), checked_temperature(:, :)
    real(dp), allocatable :: species_energy_x(:, :), species_energy_y(:, :)
    real(dp), allocatable :: theta_cell(:, :)
    real(dp), allocatable :: candidate_flux_x(:, :, :), candidate_flux_y(:, :, :)
    real(dp) :: outgoing, mass, candidate, theta_face
    real(dp) :: candidate_minimum_theta
    logical :: local_ok, periodic_x, periodic_y
    integer :: i, j, face_i, face_j, k, component
    integer :: nspecies, nvar, nprim, left_i, right_i, lower_j, upper_j

    flux_x = 0.0_dp
    flux_y = 0.0_dp
    minimum_theta = 1.0_dp
    ok = .false.
    nspecies = size(species)
    if (nspecies < 1 .or. nspecies > 32) return
    nvar = reactive_nvar(nspecies)
    nprim = reactive_nprim(nspecies)
    if (.not. all(ieee_is_finite([dx, dy, dt]))) return
    if (size(transport) /= nspecies .or. nx < 2 .or. ny < 2) return
    if (dx <= 0.0_dp .or. dy <= 0.0_dp .or. dt < 0.0_dp) return
    if (size(state, 1) /= nvar .or. size(state, 2) < nx .or. &
        size(state, 3) < ny .or. size(temperature, 1) < nx .or. &
        size(temperature, 2) < ny) return
    if (size(flux_x, 1) /= nvar .or. ubound(flux_x, 2) < nx .or. &
        size(flux_x, 3) < ny .or. size(flux_y, 1) /= nvar .or. &
        size(flux_y, 2) < nx .or. ubound(flux_y, 3) < ny) return
    if (.not. all(ieee_is_finite(state(:, 1:nx, 1:ny)))) return
    if (.not. all(ieee_is_finite(temperature(1:nx, 1:ny)))) return
    if (present(boundaries)) then
      active_boundaries = boundaries
    else
      call initialize_periodic_boundary_set_2d(nprim, active_boundaries)
    end if
    call validate_reactive_boundary_set_2d(active_boundaries, local_ok)
    if (.not. local_ok) return
    if (present(exterior)) then
      if (.not. exterior%is_valid(nprim, nx, ny)) return
    end if
    if (.not. species_diffusion_enabled .and. &
        (reactive_boundary_has_prescribed_species_flux( &
           active_boundaries%face(1)) .or. &
         reactive_boundary_has_prescribed_species_flux( &
           active_boundaries%face(2)) .or. &
         reactive_boundary_has_prescribed_species_flux( &
           active_boundaries%face(3)) .or. &
         reactive_boundary_has_prescribed_species_flux( &
           active_boundaries%face(4)))) return
    periodic_x = reactive_boundary_is_periodic(active_boundaries%face(1))
    periodic_y = reactive_boundary_is_periodic(active_boundaries%face(3))
    if (.not. (viscosity_enabled .or. thermal_conduction_enabled .or. &
        species_diffusion_enabled)) then
      ok = .true.
      return
    end if

    allocate(candidate_flux_x(nvar, 0:nx, ny), candidate_flux_y(nvar, nx, 0:ny))
    candidate_flux_x = 0.0_dp
    candidate_flux_y = 0.0_dp
    candidate_minimum_theta = 1.0_dp
    allocate(primitive(nprim, nx, ny), checked_temperature(nx, ny))
    allocate(species_energy_x(0:nx, ny), species_energy_y(nx, 0:ny))
    allocate(theta_cell(nx, ny))
    species_energy_x = 0.0_dp
    species_energy_y = 0.0_dp
    call recover_primitives_2d( &
      species, state, temperature, nx, ny, primitive, checked_temperature, &
      local_ok)
    if (.not. local_ok) return

    do j = 1, ny
      do face_i = 0, nx
        call transport_face_flux_x( &
          species, transport, primitive, checked_temperature, nx, ny, face_i, &
          j, dx, dy, active_boundaries, viscosity_enabled, &
          thermal_conduction_enabled, species_diffusion_enabled, &
          barodiffusion_enabled, candidate_flux_x(:, face_i, j), &
          species_energy_x(face_i, j), local_ok, exterior)
        if (.not. local_ok) return
      end do
    end do
    do face_j = 0, ny
      do i = 1, nx
        call transport_face_flux_y( &
          species, transport, primitive, checked_temperature, nx, ny, i, &
          face_j, dx, dy, active_boundaries, viscosity_enabled, &
          thermal_conduction_enabled, species_diffusion_enabled, &
          barodiffusion_enabled, candidate_flux_y(:, i, face_j), &
          species_energy_y(i, face_j), local_ok, exterior)
        if (.not. local_ok) return
      end do
    end do

    if (species_diffusion_enabled .and. dt > 0.0_dp) then
      theta_cell = 1.0_dp
      do j = 1, ny
        do i = 1, nx
          do k = 1, nspecies
            component = reactive_species_component(k)
            call validate_limiter_candidate_2d( &
              candidate_flux_x(component, i, j), &
              candidate_flux_x(component, i - 1, j), &
              candidate_flux_y(component, i, j), &
              candidate_flux_y(component, i, j - 1), &
              max(0.0_dp, state(component, i, j)), dx, dy, dt, local_ok)
            if (.not. local_ok) return
            outgoing = max(candidate_flux_x(component, i, j), 0.0_dp) / dx + &
              max(-candidate_flux_x(component, i - 1, j), 0.0_dp) / dx + &
              max(candidate_flux_y(component, i, j), 0.0_dp) / dy + &
              max(-candidate_flux_y(component, i, j - 1), 0.0_dp) / dy
            if (.not. ieee_is_finite(outgoing)) return
            mass = max(0.0_dp, state(component, i, j))
            if (outgoing > 0.0_dp) then
              candidate = species_safety * mass / (dt * outgoing)
              if (.not. ieee_is_finite(candidate)) return
              theta_cell(i, j) = min(theta_cell(i, j), &
                max(0.0_dp, min(1.0_dp, candidate)))
            end if
          end do
        end do
      end do

      do j = 1, ny
        do face_i = 0, nx
          theta_face = 1.0_dp
          left_i = face_i
          right_i = face_i + 1
          if (face_i == 0) then
            if (periodic_x) left_i = nx
          else if (face_i == nx) then
            if (periodic_x) right_i = 1
          end if
          if (left_i >= 1 .and. left_i <= nx) &
            theta_face = min(theta_face, theta_cell(left_i, j))
          if (right_i >= 1 .and. right_i <= nx) &
            theta_face = min(theta_face, theta_cell(right_i, j))
          candidate_minimum_theta = min(candidate_minimum_theta, theta_face)
          do k = 1, nspecies
            component = reactive_species_component(k)
            if (.not. finite_product_fits(theta_face, &
                candidate_flux_x(component, face_i, j))) return
            candidate_flux_x(component, face_i, j) = &
              theta_face * candidate_flux_x(component, face_i, j)
          end do
          if (.not. finite_product_fits(theta_face - 1.0_dp, &
              species_energy_x(face_i, j))) return
          if (.not. finite_sum_fits(candidate_flux_x(iet, face_i, j), &
              (theta_face - 1.0_dp) * species_energy_x(face_i, j))) return
          candidate_flux_x(iet, face_i, j) = candidate_flux_x(iet, face_i, j) + &
            (theta_face - 1.0_dp) * species_energy_x(face_i, j)
        end do
      end do
      do face_j = 0, ny
        do i = 1, nx
          theta_face = 1.0_dp
          lower_j = face_j
          upper_j = face_j + 1
          if (face_j == 0) then
            if (periodic_y) lower_j = ny
          else if (face_j == ny) then
            if (periodic_y) upper_j = 1
          end if
          if (lower_j >= 1 .and. lower_j <= ny) &
            theta_face = min(theta_face, theta_cell(i, lower_j))
          if (upper_j >= 1 .and. upper_j <= ny) &
            theta_face = min(theta_face, theta_cell(i, upper_j))
          candidate_minimum_theta = min(candidate_minimum_theta, theta_face)
          do k = 1, nspecies
            component = reactive_species_component(k)
            if (.not. finite_product_fits(theta_face, &
                candidate_flux_y(component, i, face_j))) return
            candidate_flux_y(component, i, face_j) = &
              theta_face * candidate_flux_y(component, i, face_j)
          end do
          if (.not. finite_product_fits(theta_face - 1.0_dp, &
              species_energy_y(i, face_j))) return
          if (.not. finite_sum_fits(candidate_flux_y(iet, i, face_j), &
              (theta_face - 1.0_dp) * species_energy_y(i, face_j))) return
          candidate_flux_y(iet, i, face_j) = candidate_flux_y(iet, i, face_j) + &
            (theta_face - 1.0_dp) * species_energy_y(i, face_j)
        end do
      end do
    end if

    if (.not. all(ieee_is_finite(candidate_flux_x))) return
    if (.not. all(ieee_is_finite(candidate_flux_y))) return
    if (.not. ieee_is_finite(candidate_minimum_theta)) return
    if (candidate_minimum_theta < 0.0_dp) return
    if (candidate_minimum_theta > 1.0_dp) return
    flux_x(:, 0:nx, 1:ny) = candidate_flux_x
    flux_y(:, 1:nx, 0:ny) = candidate_flux_y
    minimum_theta = candidate_minimum_theta
    ok = .true.
  end subroutine reactive_transport_fluxes_2d_faces

  subroutine reactive_transport_fluxes_2d( &
      species, transport, state, temperature, nx, ny, dx, dy, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, flux_x, flux_y, &
      minimum_theta, ok, boundaries)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: dx, dy, dt
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: flux_x(:, :, :), flux_y(:, :, :)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    type(reactive_boundary_set_2d), intent(in), optional :: boundaries

    real(dp), allocatable :: face_x(:, :, :), face_y(:, :, :)
    logical :: local_ok
    integer :: nvar, nspecies, i, j

    flux_x = 0.0_dp
    flux_y = 0.0_dp
    minimum_theta = 1.0_dp
    ok = .false.
    nspecies = size(species)
    if (nspecies < 1 .or. nspecies > 32) return
    nvar = reactive_nvar(nspecies)
    if (.not. all(ieee_is_finite([dx, dy, dt]))) return
    if (size(transport) /= nspecies .or. nx < 2 .or. ny < 2) return
    if (dx <= 0.0_dp .or. dy <= 0.0_dp .or. dt < 0.0_dp) return
    if (size(state, 1) /= nvar .or. size(state, 2) < nx .or. &
        size(state, 3) < ny .or. size(temperature, 1) < nx .or. &
        size(temperature, 2) < ny) return
    if (size(flux_x, 1) /= nvar .or. size(flux_x, 2) < nx .or. &
        size(flux_x, 3) < ny .or. size(flux_y, 1) /= nvar .or. &
        size(flux_y, 2) < nx .or. size(flux_y, 3) < ny) return
    if (.not. all(ieee_is_finite(state(:, 1:nx, 1:ny)))) return
    if (.not. all(ieee_is_finite(temperature(1:nx, 1:ny)))) return
    allocate(face_x(nvar, 0:nx, ny), face_y(nvar, nx, 0:ny))
    if (present(boundaries)) then
      call reactive_transport_fluxes_2d_faces( &
        species, transport, state, temperature, nx, ny, dx, dy, dt, &
        viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, face_x, face_y, &
        minimum_theta, local_ok, boundaries)
    else
      call reactive_transport_fluxes_2d_faces( &
        species, transport, state, temperature, nx, ny, dx, dy, dt, &
        viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, face_x, face_y, &
        minimum_theta, local_ok)
    end if
    if (.not. local_ok) return
    do j = 1, ny
      do i = 1, nx
        flux_x(:, i, j) = face_x(:, i, j)
        flux_y(:, i, j) = face_y(:, i, j)
      end do
    end do
    ok = .true.
  end subroutine reactive_transport_fluxes_2d

  subroutine reactive_transport_timestep_2d( &
      species, transport, state, temperature, nx, ny, dx, dy, transport_cfl, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, dt, maximum_diffusivity, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: dx, dy, transport_cfl
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled
    real(dp), intent(out) :: dt, maximum_diffusivity
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:), y(:), diffusion(:)
    real(dp) :: checked_temperature, sound_speed, viscosity, conductivity
    real(dp) :: molecular_weight, gas_constant, cp, cv, gamma
    real(dp) :: enthalpy, internal_energy, entropy, candidate, denominator
    real(dp) :: guard_value, guard_value_two, guard_denominator
    real(dp) :: candidate_maximum_diffusivity, candidate_dt
    logical :: local_ok
    integer :: i, j, k, nspecies

    dt = 0.0_dp
    maximum_diffusivity = 0.0_dp
    ok = .false.
    nspecies = size(species)
    if (nspecies < 1 .or. nspecies > 32) return
    if (.not. all(ieee_is_finite([dx, dy, transport_cfl]))) return
    if (size(transport) /= nspecies .or. nx < 1 .or. ny < 1) return
    if (dx <= 0.0_dp .or. dy <= 0.0_dp .or. transport_cfl <= 0.0_dp .or. &
        transport_cfl > 0.5_dp) return
    if (size(state, 1) /= reactive_nvar(nspecies) .or. &
        size(state, 2) < nx .or. size(state, 3) < ny .or. &
        size(temperature, 1) < nx .or. size(temperature, 2) < ny) return
    if (.not. all(ieee_is_finite(state(:, 1:nx, 1:ny)))) return
    if (.not. all(ieee_is_finite(temperature(1:nx, 1:ny)))) return
    if (.not. (viscosity_enabled .or. thermal_conduction_enabled .or. &
        species_diffusion_enabled)) then
      dt = huge(1.0_dp)
      ok = .true.
      return
    end if
    candidate_maximum_diffusivity = 0.0_dp
    allocate(primitive(reactive_nprim(nspecies)), y(nspecies))
    allocate(diffusion(nspecies))
    do j = 1, ny
      do i = 1, nx
        call reactive_conserved_to_primitive( &
          species, state(:, i, j), temperature(i, j), primitive, &
          checked_temperature, sound_speed, local_ok)
        if (.not. local_ok) return
        do k = 1, nspecies
          y(k) = primitive(reactive_mass_fraction_component(k))
        end do
        call mixture_transport_coefficients( &
          species, transport, y, checked_temperature, primitive(5), viscosity, &
          conductivity, diffusion, local_ok)
        if (.not. local_ok) return
        call mixture_mass_properties( &
          species, y, checked_temperature, molecular_weight, gas_constant, &
          cp, cv, gamma, enthalpy, internal_energy, entropy, local_ok)
        if (.not. local_ok) return
        if (.not. ieee_is_finite(primitive(1))) return
        if (.not. ieee_is_finite(cv)) return
        if (primitive(1) <= 0.0_dp) return
        if (cv <= 0.0_dp) return
        if (viscosity_enabled) then
          call multiply_finite(4.0_dp / 3.0_dp, viscosity, &
            guard_value, local_ok)
          if (.not. local_ok) return
          call divide_finite(guard_value, primitive(1), guard_value_two, local_ok)
          if (.not. local_ok) return
          candidate = (4.0_dp / 3.0_dp) * viscosity / primitive(1)
          if (.not. ieee_is_finite(candidate)) return
          candidate_maximum_diffusivity = &
            max(candidate_maximum_diffusivity, candidate)
        end if
        if (thermal_conduction_enabled) then
          call multiply_finite(primitive(1), cv, guard_denominator, local_ok)
          if (.not. local_ok) return
          call divide_finite(conductivity, guard_denominator, &
            guard_value, local_ok)
          if (.not. local_ok) return
          candidate = conductivity / (primitive(1) * cv)
          if (.not. ieee_is_finite(candidate)) return
          candidate_maximum_diffusivity = &
            max(candidate_maximum_diffusivity, candidate)
        end if
        if (species_diffusion_enabled) then
          if (.not. all(ieee_is_finite(diffusion))) return
          if (any(diffusion < 0.0_dp)) return
          candidate_maximum_diffusivity = max( &
            candidate_maximum_diffusivity, maxval(diffusion))
        end if
      end do
    end do
    if (.not. ieee_is_finite(candidate_maximum_diffusivity)) return
    if (candidate_maximum_diffusivity <= 0.0_dp) then
      candidate_dt = huge(1.0_dp)
    else
      call validate_timestep_arithmetic_2d( &
        dx, dy, transport_cfl, candidate_maximum_diffusivity, local_ok)
      if (.not. local_ok) return
      denominator = candidate_maximum_diffusivity * &
        (1.0_dp / dx**2 + 1.0_dp / dy**2)
      candidate_dt = transport_cfl / denominator
    end if
    if (.not. ieee_is_finite(candidate_dt)) return
    if (candidate_dt <= 0.0_dp) return
    dt = candidate_dt
    maximum_diffusivity = candidate_maximum_diffusivity
    ok = .true.
  end subroutine reactive_transport_timestep_2d

  subroutine reactive_transport_euler_update_2d( &
      species, transport, input_state, input_temperature, nx, ny, dx, dy, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, output_state, &
      output_temperature, minimum_theta, ok, boundaries)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: input_state(:, :, :), input_temperature(:, :)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: dx, dy, dt
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: output_state(:, :, :), output_temperature(:, :)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    type(reactive_boundary_set_2d), intent(in), optional :: boundaries

    real(dp), allocatable :: flux_x(:, :, :), flux_y(:, :, :), primitive(:)
    real(dp), allocatable :: candidate_output_state(:, :, :)
    real(dp), allocatable :: candidate_output_temperature(:, :)
    real(dp) :: checked_temperature, sound_speed, candidate_minimum_theta
    logical :: local_ok
    integer :: i, j, k, nvar, nspecies

    minimum_theta = 1.0_dp
    ok = .false.
    nspecies = size(species)
    if (nspecies < 1 .or. nspecies > 32) return
    nvar = reactive_nvar(nspecies)
    if (size(output_state, 1) /= nvar .or. size(output_state, 2) < nx .or. &
        size(output_state, 3) < ny .or. size(output_temperature, 1) < nx .or. &
        size(output_temperature, 2) < ny) return
    output_state = 0.0_dp
    output_temperature = 0.0_dp
    if (.not. all(ieee_is_finite([dx, dy, dt]))) return
    if (size(transport) /= nspecies .or. nx < 2 .or. ny < 2) return
    if (dx <= 0.0_dp .or. dy <= 0.0_dp .or. dt < 0.0_dp) return
    if (size(input_state, 1) /= nvar .or. size(input_state, 2) < nx .or. &
        size(input_state, 3) < ny .or. size(input_temperature, 1) < nx .or. &
        size(input_temperature, 2) < ny) return
    if (.not. all(ieee_is_finite(input_state(:, 1:nx, 1:ny)))) return
    if (.not. all(ieee_is_finite(input_temperature(1:nx, 1:ny)))) return
    allocate(flux_x(nvar, 0:nx, ny), flux_y(nvar, nx, 0:ny))
    allocate(primitive(reactive_nprim(nspecies)))
    allocate(candidate_output_state(nvar, nx, ny), &
      candidate_output_temperature(nx, ny))
    if (present(boundaries)) then
      call reactive_transport_fluxes_2d_faces( &
        species, transport, input_state, input_temperature, nx, ny, dx, dy, &
        dt, viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, flux_x, flux_y, &
        candidate_minimum_theta, local_ok, boundaries)
    else
      call reactive_transport_fluxes_2d_faces( &
        species, transport, input_state, input_temperature, nx, ny, dx, dy, &
        dt, viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, flux_x, flux_y, &
        candidate_minimum_theta, local_ok)
    end if
    if (.not. local_ok) return
    candidate_output_state = input_state(:, 1:nx, 1:ny)
    candidate_output_temperature = input_temperature(1:nx, 1:ny)
    do j = 1, ny
      do i = 1, nx
        do k = 1, nvar
          call validate_euler_component_arithmetic_2d( &
            input_state(k, i, j), flux_x(k, i, j), flux_x(k, i - 1, j), &
            flux_y(k, i, j), flux_y(k, i, j - 1), dx, dy, dt, local_ok)
          if (.not. local_ok) return
        end do
        candidate_output_state(:, i, j) = input_state(:, i, j) - &
          dt / dx * (flux_x(:, i, j) - flux_x(:, i - 1, j)) - &
          dt / dy * (flux_y(:, i, j) - flux_y(:, i, j - 1))
        call reactive_conserved_to_primitive( &
          species, candidate_output_state(:, i, j), input_temperature(i, j), primitive, &
          checked_temperature, sound_speed, local_ok)
        if (.not. local_ok) return
        candidate_output_temperature(i, j) = checked_temperature
      end do
    end do
    output_state(:, 1:nx, 1:ny) = candidate_output_state
    output_temperature(1:nx, 1:ny) = candidate_output_temperature
    minimum_theta = candidate_minimum_theta
    ok = .true.
  end subroutine reactive_transport_euler_update_2d

  subroutine advance_reactive_transport_2d( &
      species, transport, state, temperature, nx, ny, dx, dy, interval, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, minimum_theta, ok, &
      boundaries)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(inout) :: state(:, :, :), temperature(:, :)
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: dx, dy, interval
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    type(reactive_boundary_set_2d), intent(in), optional :: boundaries

    real(dp), allocatable :: initial_state(:, :, :), initial_temperature(:, :)
    real(dp), allocatable :: stage1_state(:, :, :), stage1_temperature(:, :)
    real(dp), allocatable :: euler2_state(:, :, :), euler2_temperature(:, :)
    real(dp), allocatable :: candidate_state(:, :, :), candidate_temperature(:, :)
    real(dp), allocatable :: primitive(:)
    real(dp) :: theta1, theta2, checked_temperature, sound_speed
    real(dp) :: temperature_guess, guard_value
    logical :: local_ok
    integer :: i, j, k, nvar, nspecies

    ok = .false.
    minimum_theta = 1.0_dp
    nspecies = size(species)
    if (nspecies < 1 .or. nspecies > 32) return
    nvar = reactive_nvar(nspecies)
    if (.not. all(ieee_is_finite([dx, dy, interval]))) return
    if (size(transport) /= nspecies .or. nx < 2 .or. ny < 2) return
    if (dx <= 0.0_dp .or. dy <= 0.0_dp .or. interval < 0.0_dp) return
    if (size(state, 1) /= nvar .or. size(state, 2) < nx .or. &
        size(state, 3) < ny .or. size(temperature, 1) < nx .or. &
        size(temperature, 2) < ny) return
    if (.not. all(ieee_is_finite(state(:, 1:nx, 1:ny)))) return
    if (.not. all(ieee_is_finite(temperature(1:nx, 1:ny)))) return
    if (interval <= tiny(1.0_dp) .or. .not. (viscosity_enabled .or. &
        thermal_conduction_enabled .or. species_diffusion_enabled)) then
      ok = .true.
      return
    end if
    allocate(initial_state(nvar, nx, ny), initial_temperature(nx, ny))
    allocate(stage1_state(nvar, nx, ny), stage1_temperature(nx, ny))
    allocate(euler2_state(nvar, nx, ny), euler2_temperature(nx, ny))
    allocate(candidate_state(nvar, nx, ny), candidate_temperature(nx, ny))
    allocate(primitive(reactive_nprim(nspecies)))
    initial_state = state(:, 1:nx, 1:ny)
    initial_temperature = temperature(1:nx, 1:ny)

    if (present(boundaries)) then
      call reactive_transport_euler_update_2d( &
        species, transport, initial_state, initial_temperature, nx, ny, dx, &
        dy, interval, viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, stage1_state, &
        stage1_temperature, theta1, local_ok, boundaries)
    else
      call reactive_transport_euler_update_2d( &
        species, transport, initial_state, initial_temperature, nx, ny, dx, &
        dy, interval, viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, stage1_state, &
        stage1_temperature, theta1, local_ok)
    end if
    if (.not. local_ok) return
    if (present(boundaries)) then
      call reactive_transport_euler_update_2d( &
        species, transport, stage1_state, stage1_temperature, nx, ny, dx, dy, &
        interval, viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, euler2_state, &
        euler2_temperature, theta2, local_ok, boundaries)
    else
      call reactive_transport_euler_update_2d( &
        species, transport, stage1_state, stage1_temperature, nx, ny, dx, dy, &
        interval, viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, euler2_state, &
        euler2_temperature, theta2, local_ok)
    end if
    if (.not. local_ok) return

    candidate_state = initial_state
    candidate_temperature = initial_temperature
    do j = 1, ny
      do i = 1, nx
        do k = 1, nvar
          call checked_half_sum_2d( &
            initial_state(k, i, j), euler2_state(k, i, j), guard_value, &
            local_ok)
          if (.not. local_ok) return
        end do
        call checked_half_sum_2d( &
          initial_temperature(i, j), euler2_temperature(i, j), &
          guard_value, local_ok)
        if (.not. local_ok) return
        candidate_state(:, i, j) = 0.5_dp * &
          (initial_state(:, i, j) + euler2_state(:, i, j))
        temperature_guess = 0.5_dp * &
          (initial_temperature(i, j) + euler2_temperature(i, j))
        call reactive_conserved_to_primitive( &
          species, candidate_state(:, i, j), temperature_guess, primitive, &
          checked_temperature, sound_speed, local_ok)
        if (.not. local_ok) return
        candidate_temperature(i, j) = checked_temperature
      end do
    end do
    state(:, 1:nx, 1:ny) = candidate_state
    temperature(1:nx, 1:ny) = candidate_temperature
    minimum_theta = min(theta1, theta2)
    ok = .true.
  end subroutine advance_reactive_transport_2d


end module reactive_transport_2d_mod
