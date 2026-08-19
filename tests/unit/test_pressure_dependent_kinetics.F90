program test_pressure_dependent_kinetics
  use precision_mod, only: dp
  use elementary_kinetics_mod, only: &
    elementary_reaction, effective_third_body_concentration, &
    troe_falloff_factor, reaction_effective_forward_rate
  use h2o2_full_mechanism_mod, only: &
    h2o2_full_nspecies, load_h2o2_full_mechanism
  implicit none

  type(elementary_reaction), allocatable :: reactions(:)
  real(dp) :: concentrations(h2o2_full_nspecies), derivative(h2o2_full_nspecies)
  real(dp) :: collider, factor, rate_low, rate_high
  logical :: ok

  call load_h2o2_full_mechanism(reactions, ok)
  if (.not. ok) error stop "Failed to load full mechanism"

  concentrations = 1.0e-3_dp
  collider = effective_third_body_concentration(reactions(1), concentrations, ok)
  if (.not. ok) error stop "Third-body concentration failed"
  call assert_close(collider, sum(reactions(1)%third_body_efficiencies) * &
    1.0e-3_dp, 2.0e-13_dp, "third-body collider")

  factor = troe_falloff_factor(reactions(22)%troe, 1000.0_dp, 1.0_dp, ok)
  if (.not. ok .or. factor <= 0.0_dp .or. factor > 1.0_dp) then
    error stop "Troe factor is outside (0,1]"
  end if

  concentrations = 1.0e-12_dp
  call reaction_effective_forward_rate( &
    reactions(22), 1000.0_dp, concentrations, rate_low, derivative, ok)
  if (.not. ok .or. rate_low <= 0.0_dp) error stop "Low-pressure rate failed"
  concentrations = 1.0e3_dp
  call reaction_effective_forward_rate( &
    reactions(22), 1000.0_dp, concentrations, rate_high, derivative, ok)
  if (.not. ok .or. rate_high <= rate_low) then
    error stop "Falloff rate did not increase with collider concentration"
  end if
  if (any(derivative < 0.0_dp)) then
    error stop "Falloff concentration derivative became negative"
  end if

  write(*, '(a)') "test_pressure_dependent_kinetics: PASS"

contains

  subroutine assert_close(actual, expected, relative_tolerance, label)
    real(dp), intent(in) :: actual, expected, relative_tolerance
    character(len=*), intent(in) :: label
    real(dp) :: error
    error = abs(actual - expected) / max(1.0_dp, abs(expected))
    if (error > relative_tolerance) then
      write(*, '(a,3(1x,es24.16))') trim(label), actual, expected, error
      error stop "Pressure-dependent kinetics mismatch"
    end if
  end subroutine assert_close
end program test_pressure_dependent_kinetics
