program test_full_h2o2_jacobian
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density, &
    mixture_mass_properties
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_full_thermo_mod, only: full_nspecies, load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: &
    load_h2o2_full_mechanism, h2o2_full_mass_fraction_jacobian
  use constant_volume_reactor_mod, only: &
    reactor_rhs, reactor_reduced_jacobian
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp) :: x(full_nspecies), y(full_nspecies), y_plus(full_nspecies)
  real(dp) :: y_minus(full_nspecies), rhs_plus(full_nspecies)
  real(dp) :: rhs_minus(full_nspecies)
  real(dp) :: jacobian(full_nspecies, full_nspecies), density
  real(dp) :: reduced_jacobian(full_nspecies - 1, full_nspecies - 1)
  real(dp) :: direction(full_nspecies - 1), projected(full_nspecies - 1)
  real(dp) :: finite_difference(full_nspecies - 1)
  real(dp) :: directions(full_nspecies - 1, full_nspecies)
  real(dp) :: molecular_weight, gas_constant, cp, cv, gamma
  real(dp) :: enthalpy, target_energy, entropy
  real(dp) :: temperature, evaluated_temperature
  real(dp) :: temperature_plus, temperature_minus
  real(dp) :: perturbation, directional_error, directional_scale
  real(dp) :: last_species_direction, maximum_directional_error
  real(dp) :: maximum_temperature_span
  logical :: ok
  integer :: direction_index, species_index

  call load_h2o2_full_thermo(species, ok)
  if (.not. ok) error stop "Failed to load full thermo"
  call load_h2o2_full_mechanism(reactions, ok)
  if (.not. ok) error stop "Failed to load full mechanism"
  x = [2.0_dp, 1.0e-5_dp, 2.0e-6_dp, 1.0_dp, 3.0e-6_dp, &
    1.0e-4_dp, 1.0e-5_dp, 1.0e-6_dp, 0.1_dp, 3.0_dp]
  x = x / sum(x)
  call mass_fractions_from_mole_fractions(species, x, y, ok)
  if (.not. ok) error stop "Composition conversion failed"
  density = mixture_density(species, y, 202650.0_dp, 1100.0_dp, ok)
  if (.not. ok) error stop "Density evaluation failed"
  call h2o2_full_mass_fraction_jacobian( &
    species, reactions, 1100.0_dp, density, y, jacobian, ok)
  if (.not. ok) error stop "Full mass-fraction Jacobian failed"
  if (maxval(abs(jacobian)) <= 1.0_dp) then
    error stop "Full Jacobian did not capture stiff scales"
  end if
  if (any(jacobian /= jacobian)) error stop "Full Jacobian contains NaN"

  temperature = 1100.0_dp
  call mixture_mass_properties( &
    species, y, temperature, molecular_weight, gas_constant, cp, cv, gamma, &
    enthalpy, target_energy, entropy, ok)
  if (.not. ok) error stop "Initial energy evaluation failed"
  call reactor_reduced_jacobian( &
    species, reactions, density, target_energy, y, temperature, &
    reduced_jacobian, evaluated_temperature, ok)
  if (.not. ok) error stop "Reduced reactor Jacobian failed"
  temperature = evaluated_temperature

  directions = 0.0_dp
  do direction_index = 1, full_nspecies - 1
    directions(direction_index, direction_index) = 1.0_dp
  end do
  directions(1, full_nspecies) = 0.75_dp
  directions(4, full_nspecies) = -0.25_dp
  directions(6, full_nspecies) = 0.10_dp
  maximum_directional_error = 0.0_dp
  maximum_temperature_span = 0.0_dp

  do direction_index = 1, full_nspecies
    direction = directions(:, direction_index)
    last_species_direction = -sum(direction)
    perturbation = 1.0e-7_dp
    do species_index = 1, full_nspecies - 1
      if (abs(direction(species_index)) > 0.0_dp) then
        perturbation = min( &
          perturbation, 0.05_dp * y(species_index) / &
          abs(direction(species_index)))
      end if
    end do
    if (abs(last_species_direction) > 0.0_dp) then
      perturbation = min( &
        perturbation, 0.05_dp * y(full_nspecies) / &
        abs(last_species_direction))
    end if
    if (perturbation <= 10.0_dp * epsilon(1.0_dp)) then
      error stop "Directional perturbation is too small"
    end if

    y_plus(1:full_nspecies - 1) = y(1:full_nspecies - 1) + &
      perturbation * direction
    y_minus(1:full_nspecies - 1) = y(1:full_nspecies - 1) - &
      perturbation * direction
    y_plus(full_nspecies) = 1.0_dp - sum(y_plus(1:full_nspecies - 1))
    y_minus(full_nspecies) = 1.0_dp - sum(y_minus(1:full_nspecies - 1))
    call reactor_rhs( &
      species, reactions, density, target_energy, y_plus, temperature, &
      rhs_plus, temperature_plus, ok)
    if (.not. ok) error stop "Positive directional state failed"
    call reactor_rhs( &
      species, reactions, density, target_energy, y_minus, temperature, &
      rhs_minus, temperature_minus, ok)
    if (.not. ok) error stop "Negative directional state failed"

    finite_difference = ( &
      rhs_plus(1:full_nspecies - 1) - rhs_minus(1:full_nspecies - 1)) / &
      (2.0_dp * perturbation)
    projected = matmul(reduced_jacobian, direction)
    directional_scale = max( &
      1.0_dp, maxval(abs(finite_difference)), maxval(abs(projected)))
    directional_error = maxval(abs(projected - finite_difference)) / &
      directional_scale
    maximum_directional_error = max( &
      maximum_directional_error, directional_error)
    maximum_temperature_span = max( &
      maximum_temperature_span, abs(temperature_plus - temperature_minus))
  end do

  ! The CVODE adapter eliminates the initially largest species rather than
  ! assuming that the final mechanism species is dependent.  Exercise that
  ! generalized reduced-state map independently of SUNDIALS.
  call reactor_reduced_jacobian( &
    species, reactions, density, target_energy, y, temperature, &
    reduced_jacobian, evaluated_temperature, ok, 1)
  if (.not. ok) error stop "Nonfinal-dependent reactor Jacobian failed"
  direction = 0.0_dp
  direction(3) = 1.0_dp
  perturbation = min(1.0e-7_dp, 0.05_dp * y(1), 0.05_dp * y(4))
  y_plus = y
  y_minus = y
  y_plus(4) = y_plus(4) + perturbation
  y_plus(1) = y_plus(1) - perturbation
  y_minus(4) = y_minus(4) - perturbation
  y_minus(1) = y_minus(1) + perturbation
  call reactor_rhs( &
    species, reactions, density, target_energy, y_plus, temperature, &
    rhs_plus, temperature_plus, ok)
  if (.not. ok) error stop "Positive nonfinal-dependent state failed"
  call reactor_rhs( &
    species, reactions, density, target_energy, y_minus, temperature, &
    rhs_minus, temperature_minus, ok)
  if (.not. ok) error stop "Negative nonfinal-dependent state failed"
  finite_difference = (rhs_plus(2:full_nspecies) - &
    rhs_minus(2:full_nspecies)) / (2.0_dp * perturbation)
  projected = matmul(reduced_jacobian, direction)
  directional_scale = max( &
    1.0_dp, maxval(abs(finite_difference)), maxval(abs(projected)))
  directional_error = maxval(abs(projected - finite_difference)) / &
    directional_scale
  maximum_directional_error = max( &
    maximum_directional_error, directional_error)
  maximum_temperature_span = max( &
    maximum_temperature_span, abs(temperature_plus - temperature_minus))
  if (maximum_directional_error > 2.0e-4_dp) then
    error stop "Reduced Jacobian failed directional finite difference"
  end if
  if (maximum_temperature_span <= 1.0e-10_dp) then
    error stop "Directional test did not exercise energy-constrained temperature"
  end if

  write(*, '(a,es24.16,a,es24.16)') &
    "test_full_h2o2_jacobian: PASS, max|J|=", maxval(abs(jacobian)), &
    ", directional error=", maximum_directional_error
end program test_full_h2o2_jacobian
