program test_pressure_dependent_kinetics
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_full_thermo
  use elementary_kinetics_mod, only: &
    elementary_reaction, reaction_kind_elementary, reaction_kind_three_body, &
    reaction_kind_falloff, effective_third_body_concentration, &
    troe_falloff_factor, elementary_production_rates_from_concentrations, &
    elementary_production_jacobian
  use h2o2_full_mechanism_mod, only: &
    h2o2_full_nspecies, h2o2_full_nreactions, load_h2o2_full_mechanism
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp) :: concentrations(h2o2_full_nspecies)
  real(dp) :: plus_concentrations(h2o2_full_nspecies)
  real(dp) :: minus_concentrations(h2o2_full_nspecies)
  real(dp) :: rates_plus(h2o2_full_nspecies)
  real(dp) :: rates_minus(h2o2_full_nspecies)
  real(dp) :: finite_difference(h2o2_full_nspecies)
  real(dp) :: jacobian(h2o2_full_nspecies, h2o2_full_nspecies)
  real(dp) :: third_body, expected_third_body, falloff, perturbation
  real(dp) :: scale, maximum_scaled_error
  logical :: ok
  integer :: i, j, elementary_count, third_body_count, falloff_count

  call load_h2o2_full_thermo(species, ok)
  if (.not. ok) error stop "Failed to load full H2/O2 thermodynamics"
  call load_h2o2_full_mechanism(reactions, ok)
  if (.not. ok) error stop "Failed to load full H2/O2 mechanism"
  if (size(species) /= h2o2_full_nspecies) error stop "Species count mismatch"
  if (size(reactions) /= h2o2_full_nreactions) error stop "Reaction count mismatch"

  elementary_count = 0
  third_body_count = 0
  falloff_count = 0
  do i = 1, size(reactions)
    select case (reactions(i)%kind)
    case (reaction_kind_elementary)
      elementary_count = elementary_count + 1
    case (reaction_kind_three_body)
      third_body_count = third_body_count + 1
    case (reaction_kind_falloff)
      falloff_count = falloff_count + 1
    case default
      error stop "Unknown generated reaction kind"
    end select
  end do
  if (elementary_count /= 23 .or. third_body_count /= 5 .or. &
      falloff_count /= 1) then
    error stop "Generated pressure-dependent reaction counts are incorrect"
  end if

  concentrations = 1.0e-3_dp * [ &
    1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, &
    6.0_dp, 7.0_dp, 8.0_dp, 9.0_dp, 10.0_dp ]
  third_body = effective_third_body_concentration( &
    reactions(1), concentrations, ok)
  if (.not. ok) error stop "Third-body concentration evaluation failed"
  expected_third_body = dot_product( &
    reactions(1)%third_body_efficiencies, concentrations)
  call assert_close(third_body, expected_third_body, 2.0e-14_dp, &
    "third-body concentration")

  falloff = troe_falloff_factor( &
    reactions(22)%troe, 1000.0_dp, 0.1_dp, ok)
  if (.not. ok) error stop "Troe factor evaluation failed"
  call assert_close(falloff, 0.5986808444383703_dp, 2.0e-13_dp, &
    "Troe factor")

  concentrations = [ &
    1.0e-2_dp, 1.0e-5_dp, 2.0e-5_dp, 5.0e-3_dp, 1.5e-5_dp, &
    1.0e-3_dp, 2.0e-6_dp, 5.0e-7_dp, 3.0e-3_dp, 2.0e-2_dp ]
  call elementary_production_jacobian( &
    species, reactions, 1200.0_dp, concentrations, jacobian, ok)
  if (.not. ok) error stop "Analytic production Jacobian evaluation failed"

  maximum_scaled_error = 0.0_dp
  do j = 1, h2o2_full_nspecies
    perturbation = max(1.0e-11_dp, 2.0e-6_dp * concentrations(j))
    if (perturbation >= 0.25_dp * concentrations(j)) then
      perturbation = 0.1_dp * concentrations(j)
    end if
    plus_concentrations = concentrations
    minus_concentrations = concentrations
    plus_concentrations(j) = plus_concentrations(j) + perturbation
    minus_concentrations(j) = minus_concentrations(j) - perturbation
    call elementary_production_rates_from_concentrations( &
      species, reactions, 1200.0_dp, plus_concentrations, rates_plus, ok)
    if (.not. ok) error stop "Positive finite-difference rate failed"
    call elementary_production_rates_from_concentrations( &
      species, reactions, 1200.0_dp, minus_concentrations, rates_minus, ok)
    if (.not. ok) error stop "Negative finite-difference rate failed"
    finite_difference = (rates_plus - rates_minus) / (2.0_dp * perturbation)
    do i = 1, h2o2_full_nspecies
      scale = 1.0e-5_dp + 5.0e-5_dp * &
        max(abs(jacobian(i, j)), abs(finite_difference(i)))
      maximum_scaled_error = max(maximum_scaled_error, &
        abs(jacobian(i, j) - finite_difference(i)) / scale)
    end do
  end do
  if (maximum_scaled_error > 1.0_dp) then
    write(*, '(a,es24.16)') &
      "maximum production Jacobian scaled error: ", maximum_scaled_error
    error stop "Analytic production Jacobian disagrees with finite difference"
  end if

  write(*, '(a,es24.16)') &
    "maximum production Jacobian scaled error: ", maximum_scaled_error
  write(*, '(a)') "test_pressure_dependent_kinetics: PASS"

contains

  subroutine assert_close(actual, expected, relative_tolerance, label)
    real(dp), intent(in) :: actual, expected, relative_tolerance
    character(len=*), intent(in) :: label
    real(dp) :: relative_error

    relative_error = abs(actual - expected) / max(1.0_dp, abs(expected))
    if (relative_error > relative_tolerance) then
      write(*, '(a,2(1x,es24.16),1x,es12.4)') &
        trim(label), actual, expected, relative_error
      error stop "Pressure-dependent kinetics reference mismatch"
    end if
  end subroutine assert_close

end program test_pressure_dependent_kinetics
