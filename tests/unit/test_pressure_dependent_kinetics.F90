program test_pressure_dependent_kinetics
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use elementary_kinetics_mod, only: &
    elementary_reaction, effective_third_body_concentration, &
    troe_falloff_factor, reaction_effective_forward_rate
  use h2o2_full_mechanism_mod, only: &
    h2o2_full_nspecies, load_h2o2_full_mechanism
  implicit none

  type(elementary_reaction), allocatable :: reactions(:)
  type(elementary_reaction) :: extreme
  real(dp) :: concentrations(h2o2_full_nspecies), derivative(h2o2_full_nspecies)
  real(dp) :: collider, factor, rate_low, rate_high, expected
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

  ! Lindemann high-pressure limit: a finite derivative must not overflow
  ! merely because the intermediate (1+Pr)**2 is unrepresentable.
  extreme=reactions(22)
  extreme%troe%enabled=.false.
  extreme%low_pressure_rate%pre_exponential=1.e180_dp
  extreme%low_pressure_rate%temperature_exponent=0._dp
  extreme%low_pressure_rate%activation_energy=0._dp
  extreme%high_pressure_rate%pre_exponential=1.e10_dp
  extreme%high_pressure_rate%temperature_exponent=0._dp
  extreme%high_pressure_rate%activation_energy=0._dp
  extreme%third_body_efficiencies=1._dp
  concentrations=1._dp/real(h2o2_full_nspecies,dp)
  collider=sum(concentrations)
  call reaction_effective_forward_rate(extreme,1000._dp,concentrations,rate_high,derivative,ok)
  if(.not.ok) error stop 'Representable high-pressure derivative refused'
  if(.not.ieee_is_finite(rate_high).or.any(.not.ieee_is_finite(derivative))) &
    error stop 'Nonfinite high-pressure coefficient or derivative'
  expected=1.e-160_dp/collider**2
  if(any(abs(derivative/expected-1._dp)>2.e-12_dp)) &
    error stop 'High-pressure analytic Lindemann derivative mismatch'
  call assert_close(rate_high,1.e10_dp,2.e-13_dp,'high-pressure limiting coefficient')

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
