program test_implicit_reactor
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_full_thermo
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_full_mechanism_mod, only: &
    h2o2_full_nspecies, load_h2o2_full_mechanism
  use constant_volume_reactor_mod, only: &
    reactor_specific_internal_energy, reactor_rhs, reactor_reduced_jacobian, &
    backward_euler_trial, advance_constant_volume_implicit_adaptive
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp) :: mole_fractions(h2o2_full_nspecies)
  real(dp) :: mass_fractions(h2o2_full_nspecies)
  real(dp) :: plus_state(h2o2_full_nspecies)
  real(dp) :: minus_state(h2o2_full_nspecies)
  real(dp) :: rhs_plus(h2o2_full_nspecies)
  real(dp) :: rhs_minus(h2o2_full_nspecies)
  real(dp) :: updated_state(h2o2_full_nspecies)
  real(dp) :: jacobian(h2o2_full_nspecies - 1, h2o2_full_nspecies - 1)
  real(dp) :: finite_difference(h2o2_full_nspecies - 1)
  real(dp) :: density, target_energy, recovered_temperature
  real(dp) :: plus_temperature, minus_temperature, updated_temperature
  real(dp) :: perturbation, scale, maximum_scaled_error
  real(dp) :: accepted_dt, next_dt, current_energy
  logical :: ok
  integer :: i, j, newton_iterations, rejected_attempts

  call load_h2o2_full_thermo(species, ok)
  if (.not. ok) error stop "Failed to load full H2/O2 thermodynamics"
  call load_h2o2_full_mechanism(reactions, ok)
  if (.not. ok) error stop "Failed to load full H2/O2 mechanism"

  mole_fractions = [ &
    0.280_dp, 1.0e-3_dp, 1.0e-3_dp, 0.140_dp, 1.0e-3_dp, &
    1.0e-3_dp, 1.0e-4_dp, 1.0e-4_dp, 0.020_dp, 0.5558_dp ]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  if (.not. ok) error stop "Failed to initialize Jacobian test composition"
  density = mixture_density( &
    species, mass_fractions, 101325.0_dp, 1200.0_dp, ok)
  if (.not. ok) error stop "Failed to initialize Jacobian test density"
  target_energy = reactor_specific_internal_energy( &
    species, mass_fractions, 1200.0_dp, ok)
  if (.not. ok) error stop "Failed to initialize Jacobian test energy"

  call reactor_reduced_jacobian( &
    species, reactions, density, target_energy, mass_fractions, 1200.0_dp, &
    jacobian, recovered_temperature, ok)
  if (.not. ok) error stop "Reduced reactor Jacobian evaluation failed"
  call assert_close(recovered_temperature, 1200.0_dp, 2.0e-11_dp, &
    "recovered Jacobian temperature")

  maximum_scaled_error = 0.0_dp
  do j = 1, h2o2_full_nspecies - 1
    perturbation = min(1.0e-7_dp, 1.0e-4_dp * mass_fractions(j))
    perturbation = max(perturbation, 1.0e-11_dp)
    if (perturbation >= 0.25_dp * mass_fractions(j)) &
      perturbation = 0.1_dp * mass_fractions(j)
    plus_state = mass_fractions
    minus_state = mass_fractions
    plus_state(j) = plus_state(j) + perturbation
    plus_state(h2o2_full_nspecies) = &
      plus_state(h2o2_full_nspecies) - perturbation
    minus_state(j) = minus_state(j) - perturbation
    minus_state(h2o2_full_nspecies) = &
      minus_state(h2o2_full_nspecies) + perturbation
    call reactor_rhs( &
      species, reactions, density, target_energy, plus_state, 1200.0_dp, &
      rhs_plus, plus_temperature, ok)
    if (.not. ok) error stop "Positive reduced finite difference failed"
    call reactor_rhs( &
      species, reactions, density, target_energy, minus_state, 1200.0_dp, &
      rhs_minus, minus_temperature, ok)
    if (.not. ok) error stop "Negative reduced finite difference failed"
    finite_difference = (rhs_plus(1:h2o2_full_nspecies - 1) - &
      rhs_minus(1:h2o2_full_nspecies - 1)) / (2.0_dp * perturbation)
    do i = 1, h2o2_full_nspecies - 1
      scale = 5.0e-3_dp + 2.0e-3_dp * &
        max(abs(jacobian(i, j)), abs(finite_difference(i)))
      maximum_scaled_error = max(maximum_scaled_error, &
        abs(jacobian(i, j) - finite_difference(i)) / scale)
    end do
  end do
  if (maximum_scaled_error > 1.0_dp) then
    write(*, '(a,es24.16)') &
      "maximum reduced Jacobian scaled error: ", maximum_scaled_error
    error stop "Reduced reactor Jacobian disagrees with finite difference"
  end if

  call backward_euler_trial( &
    species, reactions, density, target_energy, mass_fractions, 1200.0_dp, &
    1.0e-6_dp, 1.0e-6_dp, 1.0e-12_dp, updated_state, &
    updated_temperature, newton_iterations, ok)
  if (.not. ok) error stop "Backward Euler Newton trial failed"
  if (newton_iterations < 1 .or. any(updated_state < -1.0e-13_dp)) then
    error stop "Backward Euler returned an invalid state"
  end if
  call assert_close(sum(updated_state), 1.0_dp, 5.0e-13_dp, &
    "Backward Euler composition closure")

  mass_fractions = updated_state
  recovered_temperature = updated_temperature
  call advance_constant_volume_implicit_adaptive( &
    species, reactions, density, target_energy, 5.0e-6_dp, &
    1.0e-5_dp, 1.0e-12_dp, mass_fractions, recovered_temperature, &
    accepted_dt, next_dt, newton_iterations, rejected_attempts, ok)
  if (.not. ok) error stop "Adaptive implicit reactor step failed"
  if (accepted_dt <= 0.0_dp .or. next_dt <= 0.0_dp .or. &
      newton_iterations < 1 .or. any(mass_fractions < -1.0e-13_dp)) then
    error stop "Adaptive implicit reactor returned invalid diagnostics"
  end if
  current_energy = reactor_specific_internal_energy( &
    species, mass_fractions, recovered_temperature, ok)
  if (.not. ok) error stop "Implicit energy reconstruction failed"
  call assert_close(current_energy, target_energy, 2.0e-10_dp, &
    "implicit internal-energy conservation")

  write(*, '(a,es24.16)') &
    "maximum reduced Jacobian scaled error: ", maximum_scaled_error
  write(*, '(a,i0)') "implicit Newton iterations: ", newton_iterations
  write(*, '(a)') "test_implicit_reactor: PASS"

contains

  subroutine assert_close(actual, expected, relative_tolerance, label)
    real(dp), intent(in) :: actual, expected, relative_tolerance
    character(len=*), intent(in) :: label
    real(dp) :: relative_error

    relative_error = abs(actual - expected) / max(1.0_dp, abs(expected))
    if (relative_error > relative_tolerance) then
      write(*, '(a,2(1x,es24.16),1x,es12.4)') &
        trim(label), actual, expected, relative_error
      error stop "Implicit reactor reference mismatch"
    end if
  end subroutine assert_close

end program test_implicit_reactor
