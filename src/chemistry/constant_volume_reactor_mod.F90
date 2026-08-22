module constant_volume_reactor_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species, nasa7_mass_properties
  use mixture_thermo_mod, only: &
    valid_mixture_composition, mixture_mass_properties, &
    temperature_from_internal_energy
  use elementary_kinetics_mod, only: &
    elementary_reaction, elementary_mass_fraction_rhs, &
    elementary_mass_fraction_jacobian
  implicit none
  private

  integer, parameter, public :: reactor_max_step_attempts = 24

  public :: reactor_specific_internal_energy
  public :: reactor_rhs
  public :: advance_constant_volume_adaptive
  public :: reactor_reduced_jacobian
  public :: backward_euler_trial
  public :: advance_constant_volume_implicit_adaptive

contains

  real(dp) function reactor_specific_internal_energy( &
      species, mass_fractions, temperature, ok) result(internal_energy)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: mass_fractions(:), temperature
    logical, intent(out) :: ok

    real(dp) :: molecular_weight, gas_constant, cp, cv, gamma
    real(dp) :: enthalpy, entropy

    internal_energy = 0.0_dp
    call mixture_mass_properties( &
      species, mass_fractions, temperature, molecular_weight, gas_constant, &
      cp, cv, gamma, enthalpy, internal_energy, entropy, ok)
  end function reactor_specific_internal_energy

  subroutine reactor_rhs( &
      species, reactions, density, target_internal_energy, mass_fractions, &
      temperature_guess, derivative, temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: density, target_internal_energy
    real(dp), intent(in) :: mass_fractions(:), temperature_guess
    real(dp), intent(out) :: derivative(:), temperature
    logical, intent(out) :: ok

    derivative = 0.0_dp
    temperature = 0.0_dp
    ok = density > 0.0_dp .and. &
      size(derivative) == size(species) .and. &
      valid_mixture_composition(species, mass_fractions)
    if (.not. ok) return

    call temperature_from_internal_energy( &
      species, mass_fractions, target_internal_energy, temperature_guess, &
      temperature, ok)
    if (.not. ok) return
    call elementary_mass_fraction_rhs( &
      species, reactions, temperature, density, mass_fractions, derivative, ok)
  end subroutine reactor_rhs

  subroutine advance_constant_volume_adaptive( &
      species, reactions, density, target_internal_energy, &
      requested_time_step, relative_tolerance, absolute_tolerance, &
      mass_fractions, temperature, accepted_time_step, next_time_step, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: density, target_internal_energy
    real(dp), intent(in) :: requested_time_step
    real(dp), intent(in) :: relative_tolerance, absolute_tolerance
    real(dp), intent(inout) :: mass_fractions(:), temperature
    real(dp), intent(out) :: accepted_time_step, next_time_step
    logical, intent(out) :: ok

    real(dp), allocatable :: full_state(:), half_state(:), intermediate(:)
    real(dp) :: full_temperature, half_temperature, intermediate_temperature
    real(dp) :: trial_time_step, error_norm, factor
    logical :: full_ok, half_ok
    integer :: attempt

    accepted_time_step = 0.0_dp
    next_time_step = 0.0_dp
    ok = requested_time_step > 0.0_dp .and. density > 0.0_dp .and. &
      relative_tolerance > 0.0_dp .and. absolute_tolerance > 0.0_dp .and. &
      size(mass_fractions) == size(species) .and. &
      valid_mixture_composition(species, mass_fractions)
    if (.not. ok) return

    allocate(full_state(size(species)), half_state(size(species)))
    allocate(intermediate(size(species)))
    trial_time_step = requested_time_step

    do attempt = 1, reactor_max_step_attempts
      call rk4_trial( &
        mass_fractions, temperature, trial_time_step, &
        full_state, full_temperature, full_ok)
      if (full_ok) then
        call rk4_trial( &
          mass_fractions, temperature, 0.5_dp * trial_time_step, &
          intermediate, intermediate_temperature, half_ok)
      else
        half_ok = .false.
      end if
      if (half_ok) then
        call rk4_trial( &
          intermediate, intermediate_temperature, &
          0.5_dp * trial_time_step, half_state, half_temperature, half_ok)
      end if

      if (.not. full_ok .or. .not. half_ok) then
        trial_time_step = 0.25_dp * trial_time_step
        cycle
      end if

      error_norm = normalized_error( &
        full_state, full_temperature, half_state, half_temperature, &
        relative_tolerance, absolute_tolerance)
      if (.not. ieee_is_finite(error_norm)) then
        trial_time_step = 0.25_dp * trial_time_step
        cycle
      end if

      if (error_norm <= 1.0_dp) then
        call enforce_roundoff_composition(half_state, half_ok)
        if (.not. half_ok) then
          trial_time_step = 0.25_dp * trial_time_step
          cycle
        end if
        mass_fractions = half_state
        call temperature_from_internal_energy( &
          species, mass_fractions, target_internal_energy, half_temperature, &
          temperature, half_ok)
        if (.not. half_ok) then
          trial_time_step = 0.25_dp * trial_time_step
          cycle
        end if

        if (error_norm <= 1.0e-12_dp) then
          factor = 2.0_dp
        else
          factor = min(2.0_dp, max(0.25_dp, &
            0.9_dp * error_norm**(-0.2_dp)))
        end if
        accepted_time_step = trial_time_step
        next_time_step = trial_time_step * factor
        ok = .true.
        return
      end if

      factor = max(0.1_dp, 0.9_dp * error_norm**(-0.25_dp))
      trial_time_step = trial_time_step * factor
    end do

    ok = .false.

  contains

    subroutine rk4_trial( &
        initial_state, initial_temperature, time_step, &
        updated_state, updated_temperature, trial_ok)
      real(dp), intent(in) :: initial_state(:), initial_temperature, time_step
      real(dp), intent(out) :: updated_state(:), updated_temperature
      logical, intent(out) :: trial_ok

      real(dp), allocatable :: k1(:), k2(:), k3(:), k4(:), stage(:)
      real(dp) :: temperature_1, temperature_2, temperature_3, temperature_4

      trial_ok = time_step > 0.0_dp
      updated_state = 0.0_dp
      updated_temperature = 0.0_dp
      if (.not. trial_ok) return
      allocate(k1(size(species)), k2(size(species)))
      allocate(k3(size(species)), k4(size(species)), stage(size(species)))

      call reactor_rhs( &
        species, reactions, density, target_internal_energy, initial_state, &
        initial_temperature, k1, temperature_1, trial_ok)
      if (.not. trial_ok) return

      stage = initial_state + 0.5_dp * time_step * k1
      call reactor_rhs( &
        species, reactions, density, target_internal_energy, stage, &
        temperature_1, k2, temperature_2, trial_ok)
      if (.not. trial_ok) return

      stage = initial_state + 0.5_dp * time_step * k2
      call reactor_rhs( &
        species, reactions, density, target_internal_energy, stage, &
        temperature_2, k3, temperature_3, trial_ok)
      if (.not. trial_ok) return

      stage = initial_state + time_step * k3
      call reactor_rhs( &
        species, reactions, density, target_internal_energy, stage, &
        temperature_3, k4, temperature_4, trial_ok)
      if (.not. trial_ok) return

      updated_state = initial_state + time_step * &
        (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4) / 6.0_dp
      call enforce_roundoff_composition(updated_state, trial_ok)
      if (.not. trial_ok) return
      call temperature_from_internal_energy( &
        species, updated_state, target_internal_energy, temperature_4, &
        updated_temperature, trial_ok)
    end subroutine rk4_trial

  end subroutine advance_constant_volume_adaptive

  real(dp) function normalized_error( &
      first_state, first_temperature, second_state, second_temperature, &
      relative_tolerance, absolute_tolerance) result(error_norm)
    real(dp), intent(in) :: first_state(:), first_temperature
    real(dp), intent(in) :: second_state(:), second_temperature
    real(dp), intent(in) :: relative_tolerance, absolute_tolerance

    real(dp) :: scale
    integer :: i

    error_norm = 0.0_dp
    do i = 1, size(first_state)
      scale = absolute_tolerance + relative_tolerance * &
        max(abs(first_state(i)), abs(second_state(i)))
      error_norm = max(error_norm, &
        abs(second_state(i) - first_state(i)) / scale)
    end do
    scale = max(1.0e-6_dp, &
      relative_tolerance * max(abs(first_temperature), &
        abs(second_temperature)))
    error_norm = max(error_norm, &
      abs(second_temperature - first_temperature) / scale)
  end function normalized_error

  subroutine enforce_roundoff_composition(mass_fractions, ok)
    real(dp), intent(inout) :: mass_fractions(:)
    logical, intent(out) :: ok

    real(dp), parameter :: tolerance = 2.0e-10_dp
    real(dp) :: total

    ok = .false.
    if (any(.not. ieee_is_finite(mass_fractions))) return
    if (any(mass_fractions < -tolerance)) return
    mass_fractions = max(0.0_dp, mass_fractions)
    total = sum(mass_fractions)
    if (total <= 0.0_dp .or. abs(total - 1.0_dp) > tolerance) return
    mass_fractions = mass_fractions / total
    ok = .true.
  end subroutine enforce_roundoff_composition



  subroutine reactor_reduced_jacobian( &
      species, reactions, density, target_internal_energy, mass_fractions, &
      temperature_guess, jacobian, temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: density, target_internal_energy
    real(dp), intent(in) :: mass_fractions(:), temperature_guess
    real(dp), intent(out) :: jacobian(:, :), temperature
    logical, intent(out) :: ok

    real(dp), allocatable :: fixed_temperature_jacobian(:, :)
    real(dp), allocatable :: rhs_lower(:), rhs_upper(:)
    real(dp), allocatable :: d_rhs_d_temperature(:)
    real(dp), allocatable :: species_internal_energy(:)
    real(dp) :: molecular_weight, gas_constant, cp, cv, gamma
    real(dp) :: enthalpy, mixture_internal_energy, entropy
    real(dp) :: species_cp, species_cv, species_h, species_u, species_s
    real(dp) :: lower_temperature, upper_temperature, temperature_delta
    real(dp) :: common_minimum, common_maximum, d_temperature_d_variable
    logical :: property_ok
    integer :: i, j, nspecies, reduced_size

    jacobian = 0.0_dp
    temperature = 0.0_dp
    nspecies = size(species)
    reduced_size = nspecies - 1
    ok = nspecies >= 2 .and. size(mass_fractions) == nspecies .and. &
      size(jacobian, 1) == reduced_size .and. &
      size(jacobian, 2) == reduced_size .and. density > 0.0_dp
    if (.not. ok) return
    ok = valid_mixture_composition(species, mass_fractions)
    if (.not. ok) return

    call temperature_from_internal_energy( &
      species, mass_fractions, target_internal_energy, temperature_guess, &
      temperature, ok)
    if (.not. ok) return

    allocate(fixed_temperature_jacobian(nspecies, nspecies))
    allocate(rhs_lower(nspecies), rhs_upper(nspecies))
    allocate(d_rhs_d_temperature(nspecies), species_internal_energy(nspecies))
    call elementary_mass_fraction_jacobian( &
      species, reactions, temperature, density, mass_fractions, &
      fixed_temperature_jacobian, ok)
    if (.not. ok) return

    common_minimum = maxval(species%temperature_min)
    common_maximum = minval(species%temperature_max)
    temperature_delta = max(1.0e-4_dp, 1.0e-6_dp * temperature)
    lower_temperature = max(common_minimum, temperature - temperature_delta)
    upper_temperature = min(common_maximum, temperature + temperature_delta)
    if (upper_temperature <= lower_temperature) then
      ok = .false.
      return
    end if
    call elementary_mass_fraction_rhs( &
      species, reactions, lower_temperature, density, mass_fractions, &
      rhs_lower, ok)
    if (.not. ok) return
    call elementary_mass_fraction_rhs( &
      species, reactions, upper_temperature, density, mass_fractions, &
      rhs_upper, ok)
    if (.not. ok) return
    d_rhs_d_temperature = (rhs_upper - rhs_lower) / &
      (upper_temperature - lower_temperature)

    call mixture_mass_properties( &
      species, mass_fractions, temperature, molecular_weight, gas_constant, &
      cp, cv, gamma, enthalpy, mixture_internal_energy, entropy, ok)
    if (.not. ok .or. cv <= 0.0_dp) then
      ok = .false.
      return
    end if
    do i = 1, nspecies
      call nasa7_mass_properties( &
        species(i), temperature, species_cp, species_cv, species_h, &
        species_u, species_s, property_ok)
      if (.not. property_ok) then
        ok = .false.
        return
      end if
      species_internal_energy(i) = species_u
    end do

    do j = 1, reduced_size
      d_temperature_d_variable = -(species_internal_energy(j) - &
        species_internal_energy(nspecies)) / cv
      do i = 1, reduced_size
        jacobian(i, j) = fixed_temperature_jacobian(i, j) - &
          fixed_temperature_jacobian(i, nspecies) + &
          d_rhs_d_temperature(i) * d_temperature_d_variable
      end do
    end do
    ok = all(ieee_is_finite(jacobian))
  end subroutine reactor_reduced_jacobian

  subroutine backward_euler_trial( &
      species, reactions, density, target_internal_energy, initial_state, &
      initial_temperature, time_step, relative_tolerance, absolute_tolerance, &
      updated_state, updated_temperature, newton_iterations, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: density, target_internal_energy
    real(dp), intent(in) :: initial_state(:), initial_temperature, time_step
    real(dp), intent(in) :: relative_tolerance, absolute_tolerance
    real(dp), intent(out) :: updated_state(:), updated_temperature
    integer, intent(out) :: newton_iterations
    logical, intent(out) :: ok

    integer, parameter :: maximum_newton_iterations = 20
    integer, parameter :: maximum_line_search_iterations = 24
    real(dp), allocatable :: reduced_initial(:), reduced_state(:)
    real(dp), allocatable :: trial_reduced(:), full_state(:), trial_full(:)
    real(dp), allocatable :: rhs(:), trial_rhs(:), residual(:), trial_residual(:)
    real(dp), allocatable :: jacobian(:, :), matrix(:, :), correction(:)
    real(dp) :: temperature, evaluated_temperature, trial_temperature
    real(dp) :: trial_evaluated_temperature, residual_norm
    real(dp) :: trial_norm, line_factor
    logical :: solve_ok, trial_ok
    integer :: i, iteration, line_iteration, nspecies, reduced_size

    updated_state = 0.0_dp
    updated_temperature = 0.0_dp
    newton_iterations = 0
    nspecies = size(species)
    reduced_size = nspecies - 1
    ok = nspecies >= 2 .and. size(initial_state) == nspecies .and. &
      size(updated_state) == nspecies .and. density > 0.0_dp .and. &
      time_step > 0.0_dp .and. relative_tolerance > 0.0_dp .and. &
      absolute_tolerance > 0.0_dp
    if (.not. ok) return
    ok = valid_mixture_composition(species, initial_state)
    if (.not. ok) return

    allocate(reduced_initial(reduced_size), reduced_state(reduced_size))
    allocate(trial_reduced(reduced_size), full_state(nspecies))
    allocate(trial_full(nspecies), rhs(nspecies), trial_rhs(nspecies))
    allocate(residual(reduced_size), trial_residual(reduced_size))
    allocate(jacobian(reduced_size, reduced_size))
    allocate(matrix(reduced_size, reduced_size), correction(reduced_size))
    reduced_initial = initial_state(1:reduced_size)
    reduced_state = reduced_initial
    temperature = initial_temperature

    do iteration = 1, maximum_newton_iterations
      call reconstruct_full_state(reduced_state, full_state, ok)
      if (.not. ok) return
      call reactor_rhs( &
        species, reactions, density, target_internal_energy, full_state, &
        temperature, rhs, evaluated_temperature, ok)
      if (.not. ok) return
      temperature = evaluated_temperature
      residual = reduced_state - reduced_initial - &
        time_step * rhs(1:reduced_size)
      residual_norm = implicit_residual_norm( &
        residual, reduced_state, reduced_initial, &
        relative_tolerance, absolute_tolerance)
      if (.not. ieee_is_finite(residual_norm)) then
        ok = .false.
        return
      end if
      if (residual_norm <= 5.0e-2_dp) then
        updated_state = full_state
        call enforce_roundoff_composition(updated_state, ok)
        if (.not. ok) return
        call temperature_from_internal_energy( &
          species, updated_state, target_internal_energy, temperature, &
          updated_temperature, ok)
        newton_iterations = iteration
        return
      end if

      call reactor_reduced_jacobian( &
        species, reactions, density, target_internal_energy, full_state, &
        temperature, jacobian, evaluated_temperature, ok)
      if (.not. ok) return
      temperature = evaluated_temperature
      matrix = -time_step * jacobian
      do i = 1, reduced_size
        matrix(i, i) = matrix(i, i) + 1.0_dp
      end do
      call solve_dense_linear_system( &
        matrix, -residual, correction, solve_ok)
      if (.not. solve_ok) then
        ok = .false.
        return
      end if

      line_factor = 1.0_dp
      trial_ok = .false.
      do line_iteration = 1, maximum_line_search_iterations
        trial_reduced = reduced_state + line_factor * correction
        call reconstruct_full_state(trial_reduced, trial_full, trial_ok)
        if (trial_ok) then
          trial_temperature = temperature
          call reactor_rhs( &
            species, reactions, density, target_internal_energy, trial_full, &
            trial_temperature, trial_rhs, trial_evaluated_temperature, &
            trial_ok)
          if (trial_ok) trial_temperature = trial_evaluated_temperature
        end if
        if (trial_ok) then
          trial_residual = trial_reduced - reduced_initial - &
            time_step * trial_rhs(1:reduced_size)
          trial_norm = implicit_residual_norm( &
            trial_residual, trial_reduced, reduced_initial, &
            relative_tolerance, absolute_tolerance)
          if (ieee_is_finite(trial_norm) .and. trial_norm < residual_norm) then
            reduced_state = trial_reduced
            temperature = trial_temperature
            exit
          end if
        end if
        line_factor = 0.5_dp * line_factor
      end do
      if (line_iteration > maximum_line_search_iterations) then
        ok = .false.
        return
      end if
    end do
    ok = .false.
  end subroutine backward_euler_trial

  subroutine advance_constant_volume_implicit_adaptive( &
      species, reactions, density, target_internal_energy, &
      requested_time_step, relative_tolerance, absolute_tolerance, &
      mass_fractions, temperature, accepted_time_step, next_time_step, &
      newton_iterations, rejected_attempts, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: density, target_internal_energy
    real(dp), intent(in) :: requested_time_step
    real(dp), intent(in) :: relative_tolerance, absolute_tolerance
    real(dp), intent(inout) :: mass_fractions(:), temperature
    real(dp), intent(out) :: accepted_time_step, next_time_step
    integer, intent(out) :: newton_iterations, rejected_attempts
    logical, intent(out) :: ok

    real(dp), allocatable :: full_state(:), half_state(:), intermediate(:)
    real(dp), allocatable :: extrapolated_state(:)
    real(dp) :: full_temperature, half_temperature, intermediate_temperature
    real(dp) :: extrapolated_temperature
    real(dp) :: trial_time_step, error_norm, factor
    logical :: full_ok, half_ok, extrapolated_ok
    integer :: attempt, full_iterations, half_iterations_1, half_iterations_2

    accepted_time_step = 0.0_dp
    next_time_step = 0.0_dp
    newton_iterations = 0
    rejected_attempts = 0
    ok = requested_time_step > 0.0_dp .and. density > 0.0_dp .and. &
      relative_tolerance > 0.0_dp .and. absolute_tolerance > 0.0_dp .and. &
      size(mass_fractions) == size(species)
    if (.not. ok) return
    ok = valid_mixture_composition(species, mass_fractions)
    if (.not. ok) return

    allocate(full_state(size(species)), half_state(size(species)))
    allocate(intermediate(size(species)), extrapolated_state(size(species)))
    trial_time_step = requested_time_step

    do attempt = 1, reactor_max_step_attempts
      call backward_euler_trial( &
        species, reactions, density, target_internal_energy, mass_fractions, &
        temperature, trial_time_step, relative_tolerance, absolute_tolerance, &
        full_state, full_temperature, full_iterations, full_ok)
      if (full_ok) then
        call backward_euler_trial( &
          species, reactions, density, target_internal_energy, mass_fractions, &
          temperature, 0.5_dp * trial_time_step, relative_tolerance, &
          absolute_tolerance, intermediate, intermediate_temperature, &
          half_iterations_1, half_ok)
      else
        half_ok = .false.
        half_iterations_1 = 0
      end if
      if (half_ok) then
        call backward_euler_trial( &
          species, reactions, density, target_internal_energy, intermediate, &
          intermediate_temperature, 0.5_dp * trial_time_step, &
          relative_tolerance, absolute_tolerance, half_state, &
          half_temperature, half_iterations_2, half_ok)
      else
        half_iterations_2 = 0
      end if

      if (.not. full_ok .or. .not. half_ok) then
        rejected_attempts = rejected_attempts + 1
        trial_time_step = 0.25_dp * trial_time_step
        cycle
      end if
      error_norm = normalized_error( &
        full_state, full_temperature, half_state, half_temperature, &
        relative_tolerance, absolute_tolerance)
      if (.not. ieee_is_finite(error_norm)) then
        rejected_attempts = rejected_attempts + 1
        trial_time_step = 0.25_dp * trial_time_step
        cycle
      end if

      if (error_norm <= 1.0_dp) then
        ! Backward Euler is first order.  Two half steps have half the leading
        ! local truncation error of one full step, so Richardson extrapolation
        ! cancels that term and provides the accepted second-order state.
        ! Near a positivity boundary the extrapolated composition can leave the
        ! simplex; in that case retain the already verified two-half-step state.
        extrapolated_state = 2.0_dp * half_state - full_state
        call enforce_roundoff_composition(extrapolated_state, extrapolated_ok)
        if (extrapolated_ok) then
          call temperature_from_internal_energy( &
            species, extrapolated_state, target_internal_energy, &
            half_temperature, extrapolated_temperature, extrapolated_ok)
        end if

        if (extrapolated_ok) then
          mass_fractions = extrapolated_state
          temperature = extrapolated_temperature
        else
          call enforce_roundoff_composition(half_state, half_ok)
          if (.not. half_ok) then
            rejected_attempts = rejected_attempts + 1
            trial_time_step = 0.25_dp * trial_time_step
            cycle
          end if
          mass_fractions = half_state
          call temperature_from_internal_energy( &
            species, mass_fractions, target_internal_energy, half_temperature, &
            temperature, half_ok)
          if (.not. half_ok) then
            rejected_attempts = rejected_attempts + 1
            trial_time_step = 0.25_dp * trial_time_step
            cycle
          end if
        end if
        if (error_norm <= 1.0e-12_dp) then
          factor = 2.0_dp
        else
          factor = min(2.0_dp, max(0.2_dp, &
            0.9_dp * error_norm**(-0.5_dp)))
        end if
        accepted_time_step = trial_time_step
        next_time_step = trial_time_step * factor
        newton_iterations = full_iterations + half_iterations_1 + &
          half_iterations_2
        ok = .true.
        return
      end if

      rejected_attempts = rejected_attempts + 1
      factor = max(0.1_dp, 0.9_dp * error_norm**(-0.5_dp))
      trial_time_step = trial_time_step * factor
    end do
    ok = .false.
  end subroutine advance_constant_volume_implicit_adaptive

  subroutine reconstruct_full_state(reduced_state, full_state, ok)
    real(dp), intent(in) :: reduced_state(:)
    real(dp), intent(out) :: full_state(:)
    logical, intent(out) :: ok
    real(dp), parameter :: tolerance = 5.0e-13_dp

    full_state = 0.0_dp
    ok = size(full_state) == size(reduced_state) + 1 .and. &
      all(ieee_is_finite(reduced_state))
    if (.not. ok) return
    full_state(1:size(reduced_state)) = reduced_state
    full_state(size(full_state)) = 1.0_dp - sum(reduced_state)
    if (any(full_state < -tolerance)) then
      ok = .false.
      return
    end if
    full_state = max(0.0_dp, full_state)
    if (abs(sum(full_state) - 1.0_dp) > 5.0_dp * tolerance) then
      ok = .false.
      return
    end if
    full_state = full_state / sum(full_state)
    ok = .true.
  end subroutine reconstruct_full_state

  real(dp) function implicit_residual_norm( &
      residual, state, reference_state, relative_tolerance, &
      absolute_tolerance) result(norm)
    real(dp), intent(in) :: residual(:), state(:), reference_state(:)
    real(dp), intent(in) :: relative_tolerance, absolute_tolerance
    real(dp) :: scale
    integer :: i

    norm = 0.0_dp
    do i = 1, size(residual)
      scale = absolute_tolerance + relative_tolerance * &
        max(abs(state(i)), abs(reference_state(i)))
      norm = max(norm, abs(residual(i)) / scale)
    end do
  end function implicit_residual_norm

  subroutine solve_dense_linear_system(matrix, right_hand_side, solution, ok)
    real(dp), intent(in) :: matrix(:, :), right_hand_side(:)
    real(dp), intent(out) :: solution(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: work(:, :), rhs(:), temporary_row(:)
    real(dp) :: pivot_value, factor, threshold, temporary_value
    integer :: i, j, k, pivot, n

    solution = 0.0_dp
    n = size(right_hand_side)
    ok = n >= 1 .and. size(matrix, 1) == n .and. &
      size(matrix, 2) == n .and. size(solution) == n
    if (.not. ok) return
    allocate(work(n, n), rhs(n), temporary_row(n))
    work = matrix
    rhs = right_hand_side
    threshold = 100.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, maxval(abs(work)))

    do k = 1, n - 1
      pivot = k - 1 + maxloc(abs(work(k:n, k)), dim=1)
      pivot_value = abs(work(pivot, k))
      if (pivot_value <= threshold) then
        ok = .false.
        return
      end if
      if (pivot /= k) then
        temporary_row = work(k, :)
        work(k, :) = work(pivot, :)
        work(pivot, :) = temporary_row
        temporary_value = rhs(k)
        rhs(k) = rhs(pivot)
        rhs(pivot) = temporary_value
      end if
      do i = k + 1, n
        factor = work(i, k) / work(k, k)
        work(i, k) = 0.0_dp
        do j = k + 1, n
          work(i, j) = work(i, j) - factor * work(k, j)
        end do
        rhs(i) = rhs(i) - factor * rhs(k)
      end do
    end do
    if (abs(work(n, n)) <= threshold) then
      ok = .false.
      return
    end if

    solution(n) = rhs(n) / work(n, n)
    do i = n - 1, 1, -1
      solution(i) = (rhs(i) - dot_product( &
        work(i, i + 1:n), solution(i + 1:n))) / work(i, i)
    end do
    ok = all(ieee_is_finite(solution))
  end subroutine solve_dense_linear_system

end module constant_volume_reactor_mod
