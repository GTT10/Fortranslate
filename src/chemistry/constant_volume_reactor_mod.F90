module constant_volume_reactor_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use mixture_thermo_mod, only: &
    valid_mixture_composition, mixture_mass_properties, &
    temperature_from_internal_energy
  use elementary_kinetics_mod, only: &
    elementary_reaction, elementary_mass_fraction_rhs
  implicit none
  private

  integer, parameter, public :: reactor_max_step_attempts = 24

  public :: reactor_specific_internal_energy
  public :: reactor_rhs
  public :: advance_constant_volume_adaptive

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

end module constant_volume_reactor_mod
