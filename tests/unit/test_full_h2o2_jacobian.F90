program test_full_h2o2_jacobian
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species,nasa7_dimensionless_properties,universal_gas_constant
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density, &
    mixture_mass_properties
  use elementary_kinetics_mod, only: elementary_reaction,reaction_equilibrium_constant,reaction_progress_rate, &
    reaction_effective_forward_rate,elementary_production_rates_from_concentrations,elementary_production_jacobian
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
  call check_inverse_equilibrium_range()
contains
  subroutine check_inverse_equilibrium_range()
    type(nasa7_species), allocatable :: table(:)
    type(elementary_reaction), allocatable :: mechanism(:)
    type(elementary_reaction) :: counter(1)
    real(dp) :: y(10),c(10),plus(10),minus(10),rates(10),fp(10),fm(10),jac(10,10),fd(10),expected(10)
    real(dp) :: temperatures(3)=[100._dp,80._dp,20._dp],qf,qr,q,kc,kf,unused(10),g(10)
    real(dp) :: cp,h,s,log_kc,log_coefficient,wanted,residual,mass_residual
    real(dp), parameter :: step=2e-7_dp
    integer :: thermal,signed,sp,cases
    logical :: good,kc_good
    call load_h2o2_full_thermo(table,good);call range_require(good,'pinned thermo')
    call load_h2o2_full_mechanism(mechanism,good);call range_require(good,'pinned mechanism')
    table%pinned_eos_extrapolation=.true.;cases=0
    y=1e-9_dp;y(1)=.028_dp;y(4)=.224_dp;y(10)=1-sum(y(:9))
    do thermal=1,size(temperatures)
      do signed=0,1
        c=y/table%molecular_weight
        if(signed==1) c(2)=-1e-14_dp
        call elementary_production_rates_from_concentrations(table,mechanism,temperatures(thermal),c,rates,good, &
          signed_concentrations=signed==1)
        call range_require(good,'finite cold physical/signed progress despite unrepresentable Kc')
        call elementary_production_jacobian(table,mechanism,temperatures(thermal),c,jac,good, &
          signed_concentrations=signed==1)
        call range_require(good,'finite cold physical/signed progress Jacobian')
        plus=c*(1+step);minus=c*(1-step)
        call elementary_production_rates_from_concentrations(table,mechanism,temperatures(thermal),plus,fp,good, &
          signed_concentrations=signed==1)
        call range_require(good,'positive concentration perturbation')
        call elementary_production_rates_from_concentrations(table,mechanism,temperatures(thermal),minus,fm,good, &
          signed_concentrations=signed==1)
        call range_require(good,'negative concentration perturbation')
        fd=(fp-fm)/(2*step);expected=matmul(jac,c)
        residual=maxval(abs(fd-expected))/max(maxval(abs(fd)),maxval(abs(expected)),1e-30_dp)
        mass_residual=abs(sum(table%molecular_weight*rates))/max(sum(abs(table%molecular_weight*rates)),1e-30_dp)
        print *, 'Inverse-equilibrium cold coordinate:',temperatures(thermal),signed,residual,mass_residual
        call range_require(residual<=2e-6_dp,'independent directional RHS/Jacobian comparison')
        call range_require(mass_residual<=3e-12_dp,'physical/signed reaction mass balance')
        cases=cases+1
      enddo
    enddo
    c=y/table%molecular_weight
    do sp=1,10
      call nasa7_dimensionless_properties(table(sp),80._dp,cp,h,s,good)
      call range_require(good,'independent Gibbs query');g(sp)=h-s
    enddo
    call range_require(mechanism(1)%reactant_stoich(3)==2.and.mechanism(1)%product_stoich(4)==1, &
      'counterfactual uses declared2O+M=O2+M ordering')
    log_kc=2*g(3)-g(4)-log(101325._dp/(universal_gas_constant*80._dp))
    kc=reaction_equilibrium_constant(table,mechanism(1),80._dp,kc_good)
    call range_require(.not.kc_good.and.kc==0,'standalone overflowing Kc still refuses')
    call reaction_effective_forward_rate(mechanism(1),80._dp,c,kf,unused,good)
    call range_require(good.and.kf>0,'independent forward coefficient')
    call reaction_progress_rate(table,mechanism(1),80._dp,c,qf,qr,q,good)
    call range_require(good,'finite inverse-equilibrium progress')
    wanted=exp(log(kf)+log(c(4))-log_kc)
    call range_require(wanted>0.and.abs(qr/wanted-1)<=2e-11_dp,'nonzero small reverse contribution retained')
    c(4)=0
    call reaction_progress_rate(table,mechanism(1),80._dp,c,qf,qr,q,good)
    call range_require(good.and.qr==0,'absent reverse product remains exact zero')
    c(4)=-1e-14_dp
    call reaction_progress_rate(table,mechanism(1),80._dp,c,qf,qr,q,good,signed_concentrations=.true.)
    call range_require(good.and.qr==0,'true small reverse underflow remains zero')
    c(4)=-.5_dp*y(4)/table(4)%molecular_weight
    call reaction_progress_rate(table,mechanism(1),80._dp,c,qf,qr,q,good,signed_concentrations=.true.)
    call range_require(good.and.qr<0,'signed reverse product sign retained')
    call reaction_progress_rate(table,mechanism(1),80._dp,c,qf,qr,q,good)
    call range_require(.not.good.and.qf==0.and.qr==0.and.q==0,'ordinary negative concentrations still refuse')
    ! Reverse the declared stoichiometry solely as an arithmetic countercase.
    ! Kc and its reciprocal cannot be formed here, and c(O)**2 underflows,
    ! while the actual reverse progress and its derivative are representable.
    counter(1)=mechanism(1)
    counter(1)%reactant_stoich=mechanism(1)%product_stoich
    counter(1)%product_stoich=mechanism(1)%reactant_stoich
    c=y/table%molecular_weight
    call reaction_effective_forward_rate(counter(1),80._dp,c,kf,unused,good)
    call range_require(good.and.kf>0,'countercase forward coefficient')
    log_coefficient=log(kf)+log_kc
    c(3)=exp(-.55_dp*log_coefficient)
    call reaction_effective_forward_rate(counter(1),80._dp,c,kf,unused,good)
    call range_require(good.and.kf>0,'countercase changed collider')
    log_coefficient=log(kf)+log_kc
    wanted=exp(log_coefficient+2*log(c(3)))
    call range_require(c(3)*c(3)==0.and.wanted>0,'representable progress with underflowing activity product')
    kc=reaction_equilibrium_constant(table,counter(1),80._dp,kc_good)
    call range_require(.not.kc_good.and.kc==0,'standalone microscopic Kc still refuses')
    call reaction_progress_rate(table,counter(1),80._dp,c,qf,qr,q,good)
    call range_require(good.and.abs(qr/wanted-1)<=2e-11_dp,'log action preserves finite reverse progress')
    call elementary_production_jacobian(table,counter,80._dp,c,jac,good)
    call range_require(good,'log activity derivative remains finite')
    c(3)=exp(.5_dp*(log(huge(1._dp))-2-log_coefficient))
    call reaction_progress_rate(table,counter(1),80._dp,c,qf,qr,q,good)
    call range_require(good,'large countercase progress is finite')
    call elementary_production_jacobian(table,counter,80._dp,c,jac,good)
    call range_require(.not.good,'true derivative overflow refuses without an FPE or clipping')
    print *, 'PASS: inverse-equilibrium finite progress/Jacobian:',cases,'cold physical/signed cases and range controls'
  end subroutine check_inverse_equilibrium_range
  subroutine range_require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL:',label
      error stop 1
    endif
  end subroutine range_require
end program test_full_h2o2_jacobian
