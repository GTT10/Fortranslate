program test_sundials_constant_volume_multi_context
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density, &
    mixture_mass_properties
  use elementary_kinetics_mod, only: elementary_reaction
  use sundials_constant_volume_reactor_mod, only: &
    cvode_reactor_context, cvode_reactor_context_capacity, &
    cvode_reactor_statistics, &
    initialize_constant_volume_cvode, advance_constant_volume_cvode, &
    get_constant_volume_cvode_statistics, finalize_constant_volume_cvode
  use fixture_mechanism_mod, only: &
    load_fixture_thermo_data, load_fixture_mechanism
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  logical :: ok

  call load_fixture_thermo_data(species, ok)
  if (.not. ok) error stop "Could not load fixture thermodynamics"
  call load_fixture_mechanism(reactions, ok)
  if (.not. ok) error stop "Could not load fixture reactions"

  call check_interleaved_contexts(species, reactions)
  call check_failure_isolation(species, reactions)
  call check_stale_handle_isolation(species, reactions)
  call check_context_capacity(species, reactions)

  deallocate(reactions, species)
  write(*, '(a)') "test_sundials_constant_volume_multi_context: PASS"

contains

  subroutine check_interleaved_contexts(species, reactions)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)

    real(dp), parameter :: mole_a(2) = [0.8_dp, 0.2_dp]
    real(dp), parameter :: mole_b(2) = [0.1_dp, 0.9_dp]
    real(dp), parameter :: targets_a(2) = [2.5e-8_dp, 5.0e-8_dp]
    real(dp), parameter :: targets_b(3) = &
      [1.25e-8_dp, 2.5e-8_dp, 5.0e-8_dp]
    type(cvode_reactor_context) :: context_a, context_b
    type(cvode_reactor_statistics) :: baseline_stats_a, baseline_stats_b
    type(cvode_reactor_statistics) :: interleaved_stats_a, interleaved_stats_b
    real(dp) :: baseline_y_a(2), baseline_y_b(2)
    real(dp) :: y_a(2), y_b(2), baseline_t_a, baseline_t_b
    real(dp) :: temperature_a, temperature_b, density_a, density_b
    real(dp) :: energy_a, energy_b
    character(len=1024) :: message
    logical :: ok

    call run_standalone( &
      species, reactions, mole_a, 1000.0_dp, 101325.0_dp, targets_a, &
      baseline_y_a, baseline_t_a, baseline_stats_a)
    call run_standalone( &
      species, reactions, mole_b, 1200.0_dp, 202650.0_dp, targets_b, &
      baseline_y_b, baseline_t_b, baseline_stats_b)

    call prepare_state( &
      species, mole_a, 1000.0_dp, 101325.0_dp, y_a, temperature_a, &
      density_a, energy_a)
    call prepare_state( &
      species, mole_b, 1200.0_dp, 202650.0_dp, y_b, temperature_b, &
      density_b, energy_b)
    if (maxloc(y_a, dim=1) == maxloc(y_b, dim=1)) then
      error stop "Multi-context fixture did not exercise distinct closures"
    end if

    call initialize_fixture_context( &
      context_a, species, reactions, density_a, energy_a, y_a, &
      temperature_a, 100000, ok, message)
    if (.not. ok) error stop trim(message)
    call initialize_fixture_context( &
      context_b, species, reactions, density_b, energy_b, y_b, &
      temperature_b, 100000, ok, message)
    if (.not. ok) error stop trim(message)

    call advance_constant_volume_cvode( &
      context_a, targets_a(1), y_a, temperature_a, ok, message)
    if (.not. ok) error stop trim(message)
    call advance_constant_volume_cvode( &
      context_b, targets_b(1), y_b, temperature_b, ok, message)
    if (.not. ok) error stop trim(message)
    call advance_constant_volume_cvode( &
      context_b, targets_b(2), y_b, temperature_b, ok, message)
    if (.not. ok) error stop trim(message)
    call advance_constant_volume_cvode( &
      context_a, targets_a(2), y_a, temperature_a, ok, message)
    if (.not. ok) error stop trim(message)
    call get_constant_volume_cvode_statistics( &
      context_a, interleaved_stats_a, ok, message)
    if (.not. ok) error stop trim(message)
    if (any(y_a /= baseline_y_a) .or. temperature_a /= baseline_t_a .or. &
        .not. same_statistics(interleaved_stats_a, baseline_stats_a)) then
      error stop "Interleaved context A differed from standalone execution"
    end if

    call finalize_constant_volume_cvode(context_a, ok, message)
    if (.not. ok) error stop trim(message)
    call finalize_constant_volume_cvode(context_a, ok, message)
    if (.not. ok) error stop "Context A double finalization failed"

    call advance_constant_volume_cvode( &
      context_b, targets_b(3), y_b, temperature_b, ok, message)
    if (.not. ok) error stop trim(message)
    call get_constant_volume_cvode_statistics( &
      context_b, interleaved_stats_b, ok, message)
    if (.not. ok) error stop trim(message)
    if (any(y_b /= baseline_y_b) .or. temperature_b /= baseline_t_b .or. &
        .not. same_statistics(interleaved_stats_b, baseline_stats_b)) then
      error stop "Interleaved context B differed from standalone execution"
    end if
    call finalize_constant_volume_cvode(context_b, ok, message)
    if (.not. ok) error stop trim(message)
    call finalize_constant_volume_cvode(context_b, ok, message)
    if (.not. ok) error stop "Context B double finalization failed"
  end subroutine check_interleaved_contexts

  subroutine check_failure_isolation(species, reactions)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)

    type(cvode_reactor_context) :: failing_context, surviving_context
    type(cvode_reactor_statistics) :: statistics
    real(dp) :: y_failing(2), y_surviving(2), saved_y(2)
    real(dp) :: temperature_failing, temperature_surviving
    real(dp) :: saved_temperature, density, energy
    character(len=1024) :: message
    logical :: ok

    call prepare_state( &
      species, [0.8_dp, 0.2_dp], 1000.0_dp, 101325.0_dp, y_failing, &
      temperature_failing, density, energy)
    call initialize_fixture_context( &
      failing_context, species, reactions, density, energy, y_failing, &
      temperature_failing, 1, ok, message)
    if (.not. ok) error stop trim(message)

    call prepare_state( &
      species, [0.1_dp, 0.9_dp], 1200.0_dp, 202650.0_dp, y_surviving, &
      temperature_surviving, density, energy)
    call initialize_fixture_context( &
      surviving_context, species, reactions, density, energy, y_surviving, &
      temperature_surviving, 100000, ok, message)
    if (.not. ok) error stop trim(message)

    saved_y = y_failing
    saved_temperature = temperature_failing
    call advance_constant_volume_cvode( &
      failing_context, 1.0e-7_dp, y_failing, temperature_failing, ok, message)
    if (ok .or. index(message, "FCVode failed") == 0) then
      error stop "Context-local CVODE failure was not reported"
    end if
    if (any(y_failing /= saved_y) .or. &
        temperature_failing /= saved_temperature) then
      error stop "Failed context modified its caller state"
    end if
    call advance_constant_volume_cvode( &
      failing_context, 1.25e-7_dp, y_failing, temperature_failing, ok, message)
    if (ok .or. index(message, "failed and must be finalized") == 0) then
      error stop "Failed context was not quarantined"
    end if

    call advance_constant_volume_cvode( &
      surviving_context, 2.5e-8_dp, y_surviving, temperature_surviving, &
      ok, message)
    if (.not. ok) error stop "A failed context contaminated its peer"
    call get_constant_volume_cvode_statistics( &
      surviving_context, statistics, ok, message)
    if (.not. ok .or. statistics%internal_steps <= 0) then
      error stop "Surviving context did not retain independent statistics"
    end if

    call finalize_constant_volume_cvode(surviving_context, ok, message)
    if (.not. ok) error stop trim(message)
    call finalize_constant_volume_cvode(failing_context, ok, message)
    if (.not. ok) error stop trim(message)
    call finalize_constant_volume_cvode(failing_context, ok, message)
    if (.not. ok) error stop "Failed context double finalization failed"
  end subroutine check_failure_isolation

  subroutine check_stale_handle_isolation(species, reactions)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)

    type(cvode_reactor_context) :: owner, stale_copy, replacement
    real(dp) :: mass_fractions(2), temperature, density, energy
    character(len=1024) :: message
    logical :: ok

    call prepare_state( &
      species, [0.8_dp, 0.2_dp], 1000.0_dp, 101325.0_dp, &
      mass_fractions, temperature, density, energy)
    call initialize_fixture_context( &
      owner, species, reactions, density, energy, mass_fractions, &
      temperature, 100000, ok, message)
    if (.not. ok) error stop trim(message)
    stale_copy = owner
    call finalize_constant_volume_cvode(owner, ok, message)
    if (.not. ok) error stop trim(message)

    call initialize_fixture_context( &
      replacement, species, reactions, density, energy, mass_fractions, &
      temperature, 100000, ok, message)
    if (.not. ok) error stop trim(message)
    call finalize_constant_volume_cvode(stale_copy, ok, message)
    if (ok .or. index(message, "stale") == 0) then
      error stop "Stale CVODE handle was not rejected"
    end if
    call advance_constant_volume_cvode( &
      replacement, 2.5e-8_dp, mass_fractions, temperature, ok, message)
    if (.not. ok) error stop "Stale finalization released a replacement context"
    call finalize_constant_volume_cvode(replacement, ok, message)
    if (.not. ok) error stop trim(message)
    call finalize_constant_volume_cvode(stale_copy, ok, message)
    if (.not. ok) error stop "Cleared stale handle was not idempotent"
  end subroutine check_stale_handle_isolation

  subroutine check_context_capacity(species, reactions)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)

    type(cvode_reactor_context) :: &
      contexts(cvode_reactor_context_capacity + 1)
    real(dp) :: mass_fractions(2), temperature, density, energy
    character(len=1024) :: message
    logical :: ok
    integer :: context_index

    call prepare_state( &
      species, [0.8_dp, 0.2_dp], 1000.0_dp, 101325.0_dp, &
      mass_fractions, temperature, density, energy)
    do context_index = 1, cvode_reactor_context_capacity
      call initialize_fixture_context( &
        contexts(context_index), species, reactions, density, energy, &
        mass_fractions, temperature, 100000, ok, message)
      if (.not. ok) error stop "CVODE context capacity was underfilled"
    end do
    call initialize_fixture_context( &
      contexts(cvode_reactor_context_capacity + 1), species, reactions, &
      density, energy, mass_fractions, temperature, 100000, ok, message)
    if (ok) error stop "CVODE context capacity exhaustion was not rejected"
    if (index(message, "capacity is exhausted") == 0) error stop trim(message)
    call advance_constant_volume_cvode( &
      contexts(cvode_reactor_context_capacity), 2.5e-8_dp, mass_fractions, &
      temperature, ok, message)
    if (.not. ok) error stop "Capacity rejection contaminated a live context"
    do context_index = cvode_reactor_context_capacity, 1, -1
      call finalize_constant_volume_cvode( &
        contexts(context_index), ok, message)
      if (.not. ok) error stop "Reverse-order CVODE finalization failed"
    end do
    call finalize_constant_volume_cvode( &
      contexts(cvode_reactor_context_capacity + 1), ok, message)
    if (.not. ok) error stop "Rejected capacity handle was not empty"
  end subroutine check_context_capacity

  subroutine run_standalone( &
      species, reactions, mole_fractions, initial_temperature, pressure, &
      target_times, final_mass_fractions, final_temperature, statistics)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: mole_fractions(:), initial_temperature, pressure
    real(dp), intent(in) :: target_times(:)
    real(dp), intent(out) :: final_mass_fractions(:), final_temperature
    type(cvode_reactor_statistics), intent(out) :: statistics

    type(cvode_reactor_context) :: context
    real(dp) :: density, energy
    character(len=1024) :: message
    logical :: ok
    integer :: output_index

    call prepare_state( &
      species, mole_fractions, initial_temperature, pressure, &
      final_mass_fractions, final_temperature, density, energy)
    call initialize_fixture_context( &
      context, species, reactions, density, energy, final_mass_fractions, &
      final_temperature, 100000, ok, message)
    if (.not. ok) error stop trim(message)
    do output_index = 1, size(target_times)
      call advance_constant_volume_cvode( &
        context, target_times(output_index), final_mass_fractions, &
        final_temperature, ok, message)
      if (.not. ok) error stop trim(message)
    end do
    call get_constant_volume_cvode_statistics( &
      context, statistics, ok, message)
    if (.not. ok) error stop trim(message)
    call finalize_constant_volume_cvode(context, ok, message)
    if (.not. ok) error stop trim(message)
  end subroutine run_standalone

  subroutine initialize_fixture_context( &
      context, species, reactions, density, energy, mass_fractions, &
      temperature, maximum_steps, ok, message)
    type(cvode_reactor_context), intent(inout) :: context
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: density, energy, mass_fractions(:), temperature
    integer, intent(in) :: maximum_steps
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    call initialize_constant_volume_cvode( &
      context, species, reactions, density, energy, 0.0_dp, 1.0e-12_dp, &
      1.0e-16_dp, 2.5e-8_dp, 1.0e-7_dp, 1.0e-13_dp, maximum_steps, &
      mass_fractions, temperature, ok, message)
  end subroutine initialize_fixture_context

  subroutine prepare_state( &
      species, mole_fractions, initial_temperature, pressure, &
      mass_fractions, temperature, density, energy)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: mole_fractions(:), initial_temperature, pressure
    real(dp), intent(out) :: mass_fractions(:), temperature, density, energy

    real(dp) :: molecular_weight, gas_constant, cp, cv, gamma
    real(dp) :: enthalpy, entropy
    logical :: ok

    call mass_fractions_from_mole_fractions( &
      species, mole_fractions, mass_fractions, ok)
    if (.not. ok) error stop "Could not create fixture mass fractions"
    temperature = initial_temperature
    density = mixture_density( &
      species, mass_fractions, pressure, temperature, ok)
    if (.not. ok) error stop "Could not create fixture density"
    call mixture_mass_properties( &
      species, mass_fractions, temperature, molecular_weight, gas_constant, &
      cp, cv, gamma, enthalpy, energy, entropy, ok)
    if (.not. ok) error stop "Could not create fixture energy"
  end subroutine prepare_state

  logical function same_statistics(left, right) result(same)
    type(cvode_reactor_statistics), intent(in) :: left, right

    same = left%internal_steps == right%internal_steps .and. &
      left%rhs_evaluations == right%rhs_evaluations .and. &
      left%jacobian_evaluations == right%jacobian_evaluations .and. &
      left%nonlinear_iterations == right%nonlinear_iterations .and. &
      left%nonlinear_convergence_failures == &
        right%nonlinear_convergence_failures .and. &
      left%error_test_failures == right%error_test_failures
  end function same_statistics

end program test_sundials_constant_volume_multi_context
