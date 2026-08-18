module implicit_reactor_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use pressure_dependent_kinetics_mod, only: pressure_dependent_reaction
  use full_reactor_rhs_mod, only: evaluate_full_reactor_rhs
  use h2o2_full_jacobian_mod, only: &
    h2o2_full_active_species, evaluate_h2o2_full_jacobian
  implicit none
  private

  public :: backward_euler_newton_step
  public :: advance_implicit_adaptive

contains

  subroutine backward_euler_newton_step( &
      species, reactions, density, target_internal_energy, time_step, &
      relative_tolerance, absolute_tolerance, maximum_iterations, &
      old_mass_fractions, old_temperature, new_mass_fractions, &
      new_temperature, iterations, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(pressure_dependent_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: density, target_internal_energy, time_step
    real(dp), intent(in) :: relative_tolerance, absolute_tolerance
    integer, intent(in) :: maximum_iterations
    real(dp), intent(in) :: old_mass_fractions(:), old_temperature
    real(dp), intent(out) :: new_mass_fractions(:), new_temperature
    integer, intent(out) :: iterations
    logical, intent(out) :: ok
    real(dp), allocatable :: trial(:), candidate(:), derivative(:)
    real(dp), allocatable :: residual(:), matrix(:, :), delta(:), jacobian(:, :)
    real(dp), allocatable :: candidate_derivative(:), candidate_residual(:)
    real(dp) :: residual_norm, candidate_norm, line_scale, temperature
    real(dp) :: candidate_temperature, reactive_total
    logical :: local_ok
    integer :: row, line_iteration

    new_mass_fractions = old_mass_fractions
    new_temperature = old_temperature
    iterations = 0
    ok = size(species) == size(old_mass_fractions)
    ok = ok .and. size(new_mass_fractions) == size(species)
    ok = ok .and. size(species) == h2o2_full_active_species + 1
    ok = ok .and. density > 0.0_dp .and. time_step > 0.0_dp
    ok = ok .and. relative_tolerance > 0.0_dp
    ok = ok .and. absolute_tolerance > 0.0_dp
    ok = ok .and. maximum_iterations > 0
    if (.not. ok) return

    allocate(trial(size(species)), candidate(size(species)))
    allocate(derivative(size(species)), candidate_derivative(size(species)))
    allocate(residual(h2o2_full_active_species))
    allocate(candidate_residual(h2o2_full_active_species))
    allocate(matrix(h2o2_full_active_species, h2o2_full_active_species))
    allocate(jacobian(h2o2_full_active_species, h2o2_full_active_species))
    allocate(delta(h2o2_full_active_species))

    trial = max(old_mass_fractions, 0.0_dp)
    trial = trial / sum(trial)
    reactive_total = 1.0_dp - old_mass_fractions(size(species))
    temperature = old_temperature

    do iterations = 1, maximum_iterations
      call evaluate_h2o2_full_jacobian( &
        species, reactions, density, target_internal_energy, trial, &
        temperature, derivative, jacobian, temperature, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      residual = trial(1:h2o2_full_active_species) - &
        old_mass_fractions(1:h2o2_full_active_species) - &
        time_step * derivative(1:h2o2_full_active_species)
      residual_norm = scaled_vector_norm( &
        residual, trial(1:h2o2_full_active_species), &
        old_mass_fractions(1:h2o2_full_active_species), &
        relative_tolerance, absolute_tolerance)
      if (residual_norm <= 0.2_dp) then
        new_mass_fractions = trial
        new_temperature = temperature
        ok = all(new_mass_fractions >= -1.0e-13_dp)
        return
      end if

      matrix = -time_step * jacobian
      do row = 1, h2o2_full_active_species
        matrix(row, row) = matrix(row, row) + 1.0_dp
      end do
      delta = -residual
      call solve_dense_linear_system(matrix, delta, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if

      line_scale = 1.0_dp
      do line_iteration = 1, 18
        candidate = trial
        candidate(1:h2o2_full_active_species) = &
          trial(1:h2o2_full_active_species) + line_scale * delta
        candidate(size(species)) = &
          1.0_dp - sum(candidate(1:h2o2_full_active_species))
        if (minval(candidate) >= -1.0e-13_dp .and. &
            abs(sum(candidate(1:h2o2_full_active_species)) - &
            reactive_total) <= 5.0e-7_dp) then
          candidate = max(candidate, 0.0_dp)
          candidate = candidate / sum(candidate)
          call evaluate_full_reactor_rhs( &
            species, reactions, density, target_internal_energy, candidate, &
            temperature, candidate_derivative, candidate_temperature, &
            ok=local_ok)
          if (local_ok) then
            candidate_residual = &
              candidate(1:h2o2_full_active_species) - &
              old_mass_fractions(1:h2o2_full_active_species) - &
              time_step * &
              candidate_derivative(1:h2o2_full_active_species)
            candidate_norm = scaled_vector_norm( &
              candidate_residual, &
              candidate(1:h2o2_full_active_species), &
              old_mass_fractions(1:h2o2_full_active_species), &
              relative_tolerance, absolute_tolerance)
            if (candidate_norm < residual_norm .or. line_scale <= 1.0e-4_dp) then
              trial = candidate
              temperature = candidate_temperature
              exit
            end if
          end if
        end if
        line_scale = 0.5_dp * line_scale
      end do
      if (line_iteration > 18) then
        ok = .false.
        return
      end if
    end do

    ok = .false.
  end subroutine backward_euler_newton_step

  subroutine advance_implicit_adaptive( &
      species, reactions, density, target_internal_energy, requested_step, &
      minimum_step, maximum_step, relative_tolerance, absolute_tolerance, &
      maximum_newton_iterations, mass_fractions, temperature, accepted_step, &
      suggested_step, newton_iterations, rejected_steps, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(pressure_dependent_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: density, target_internal_energy, requested_step
    real(dp), intent(in) :: minimum_step, maximum_step
    real(dp), intent(in) :: relative_tolerance, absolute_tolerance
    integer, intent(in) :: maximum_newton_iterations
    real(dp), intent(inout) :: mass_fractions(:), temperature
    real(dp), intent(out) :: accepted_step, suggested_step
    integer, intent(out) :: newton_iterations, rejected_steps
    logical, intent(out) :: ok
    real(dp), allocatable :: full_state(:), half_state(:), second_half_state(:)
    real(dp) :: full_temperature, half_temperature, second_half_temperature
    real(dp) :: trial_step, error, factor, temperature_scale
    integer :: full_iterations, half_iterations, second_iterations, attempt
    logical :: full_ok, half_ok, second_ok

    accepted_step = 0.0_dp
    suggested_step = requested_step
    newton_iterations = 0
    rejected_steps = 0
    ok = requested_step > 0.0_dp .and. minimum_step > 0.0_dp
    ok = ok .and. maximum_step >= minimum_step
    ok = ok .and. size(species) == size(mass_fractions)
    if (.not. ok) return

    allocate(full_state(size(species)), half_state(size(species)))
    allocate(second_half_state(size(species)))
    trial_step = min(maximum_step, max(minimum_step, requested_step))

    do attempt = 1, 30
      call backward_euler_newton_step( &
        species, reactions, density, target_internal_energy, trial_step, &
        relative_tolerance, absolute_tolerance, maximum_newton_iterations, &
        mass_fractions, temperature, full_state, full_temperature, &
        full_iterations, full_ok)
      call backward_euler_newton_step( &
        species, reactions, density, target_internal_energy, &
        0.5_dp * trial_step, relative_tolerance, absolute_tolerance, &
        maximum_newton_iterations, mass_fractions, temperature, half_state, &
        half_temperature, half_iterations, half_ok)
      second_ok = .false.
      second_iterations = 0
      if (half_ok) then
        call backward_euler_newton_step( &
          species, reactions, density, target_internal_energy, &
          0.5_dp * trial_step, relative_tolerance, absolute_tolerance, &
          maximum_newton_iterations, half_state, half_temperature, &
          second_half_state, second_half_temperature, second_iterations, &
          second_ok)
      end if

      if (full_ok .and. half_ok .and. second_ok) then
        error = scaled_state_difference( &
          second_half_state, full_state, relative_tolerance, &
          absolute_tolerance)
        temperature_scale = absolute_tolerance + relative_tolerance * &
          max(abs(second_half_temperature), abs(full_temperature))
        error = max(error, abs(second_half_temperature - full_temperature) / &
          max(temperature_scale, 1.0e-12_dp))
      else
        error = huge(1.0_dp)
      end if

      if (error <= 1.0_dp .and. second_ok) then
        mass_fractions = max(second_half_state, 0.0_dp)
        mass_fractions = mass_fractions / sum(mass_fractions)
        temperature = second_half_temperature
        accepted_step = trial_step
        newton_iterations = full_iterations + half_iterations + &
          second_iterations
        factor = 0.9_dp / sqrt(max(error, 1.0e-8_dp))
        factor = min(2.0_dp, max(0.25_dp, factor))
        suggested_step = min(maximum_step, max(minimum_step, &
          trial_step * factor))
        ok = all(ieee_is_finite(mass_fractions)) .and. &
          ieee_is_finite(temperature)
        return
      end if

      rejected_steps = rejected_steps + 1
      if (trial_step <= minimum_step * (1.0_dp + 10.0_dp * epsilon(1.0_dp))) then
        ok = .false.
        return
      end if
      if (ieee_is_finite(error) .and. error > 0.0_dp) then
        factor = 0.8_dp / sqrt(error)
        factor = min(0.5_dp, max(0.1_dp, factor))
      else
        factor = 0.25_dp
      end if
      trial_step = max(minimum_step, trial_step * factor)
    end do
    ok = .false.
  end subroutine advance_implicit_adaptive

  real(dp) function scaled_vector_norm( &
      vector, current, reference, relative_tolerance, absolute_tolerance) &
      result(norm)
    real(dp), intent(in) :: vector(:), current(:), reference(:)
    real(dp), intent(in) :: relative_tolerance, absolute_tolerance
    real(dp) :: scale
    integer :: index

    norm = 0.0_dp
    do index = 1, size(vector)
      scale = absolute_tolerance + relative_tolerance * &
        max(abs(current(index)), abs(reference(index)))
      norm = max(norm, abs(vector(index)) / max(scale, 1.0e-30_dp))
    end do
  end function scaled_vector_norm

  real(dp) function scaled_state_difference( &
      left, right, relative_tolerance, absolute_tolerance) result(error)
    real(dp), intent(in) :: left(:), right(:)
    real(dp), intent(in) :: relative_tolerance, absolute_tolerance
    real(dp) :: scale
    integer :: index

    error = 0.0_dp
    do index = 1, size(left)
      scale = absolute_tolerance + relative_tolerance * &
        max(abs(left(index)), abs(right(index)))
      error = max(error, abs(left(index) - right(index)) / &
        max(scale, 1.0e-30_dp))
    end do
  end function scaled_state_difference

  subroutine solve_dense_linear_system(matrix, right_hand_side, ok)
    real(dp), intent(inout) :: matrix(:, :)
    real(dp), intent(inout) :: right_hand_side(:)
    logical, intent(out) :: ok
    real(dp) :: pivot_value, multiplier, temporary
    real(dp), allocatable :: temporary_row(:)
    integer :: dimension, pivot, row, column, pivot_row

    dimension = size(right_hand_side)
    ok = size(matrix, 1) == dimension .and. size(matrix, 2) == dimension
    if (.not. ok) return
    allocate(temporary_row(dimension))

    do pivot = 1, dimension - 1
      pivot_row = pivot - 1 + maxloc(abs(matrix(pivot:dimension, pivot)), dim=1)
      pivot_value = matrix(pivot_row, pivot)
      if (.not. ieee_is_finite(pivot_value) .or. &
          abs(pivot_value) <= 100.0_dp * tiny(1.0_dp)) then
        ok = .false.
        return
      end if
      if (pivot_row /= pivot) then
        temporary_row = matrix(pivot, :)
        matrix(pivot, :) = matrix(pivot_row, :)
        matrix(pivot_row, :) = temporary_row
        temporary = right_hand_side(pivot)
        right_hand_side(pivot) = right_hand_side(pivot_row)
        right_hand_side(pivot_row) = temporary
      end if
      do row = pivot + 1, dimension
        multiplier = matrix(row, pivot) / matrix(pivot, pivot)
        matrix(row, pivot) = 0.0_dp
        do column = pivot + 1, dimension
          matrix(row, column) = matrix(row, column) - &
            multiplier * matrix(pivot, column)
        end do
        right_hand_side(row) = right_hand_side(row) - &
          multiplier * right_hand_side(pivot)
      end do
    end do
    if (abs(matrix(dimension, dimension)) <= 100.0_dp * tiny(1.0_dp)) then
      ok = .false.
      return
    end if

    right_hand_side(dimension) = right_hand_side(dimension) / &
      matrix(dimension, dimension)
    do row = dimension - 1, 1, -1
      temporary = right_hand_side(row)
      do column = row + 1, dimension
        temporary = temporary - matrix(row, column) * right_hand_side(column)
      end do
      if (abs(matrix(row, row)) <= 100.0_dp * tiny(1.0_dp)) then
        ok = .false.
        return
      end if
      right_hand_side(row) = temporary / matrix(row, row)
    end do
    ok = all(ieee_is_finite(right_hand_side))
  end subroutine solve_dense_linear_system

end module implicit_reactor_mod
