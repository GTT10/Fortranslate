module elementary_kinetics_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: &
    nasa7_species, nasa7_dimensionless_properties, universal_gas_constant
  use mixture_thermo_mod, only: valid_mixture_composition
  implicit none
  private

  real(dp), parameter, public :: kinetics_standard_pressure = 101325.0_dp
  real(dp), parameter :: logarithm_limit = 700.0_dp
  real(dp), parameter :: minimum_reduced_pressure = 1.0e-300_dp

  integer, parameter, public :: reaction_kind_elementary = 1
  integer, parameter, public :: reaction_kind_three_body = 2
  integer, parameter, public :: reaction_kind_falloff = 3
  integer, parameter, public :: reaction_kind_sri = 4
  integer, parameter, public :: reaction_kind_activated = 5
  integer, parameter, public :: reaction_kind_sri_activated = 6
  integer, parameter, public :: reaction_kind_tsang = 7
  integer, parameter, public :: reaction_kind_tsang_activated = 8
  integer, parameter, public :: reaction_kind_plog = 9
  integer, parameter, public :: plog_capacity = 64
  integer, parameter, public :: reaction_kind_chebyshev = 10
  integer, parameter, public :: chebyshev_capacity = 8

  type, public :: arrhenius_rate
    real(dp) :: pre_exponential = 0.0_dp
    real(dp) :: temperature_exponent = 0.0_dp
    real(dp) :: activation_energy = 0.0_dp ! J / kmol
  end type arrhenius_rate

  type, public :: troe_parameters
    logical :: enabled = .false.
    real(dp) :: alpha = 0.0_dp
    real(dp) :: temperature_3 = 0.0_dp
    real(dp) :: temperature_1 = 0.0_dp
    real(dp) :: temperature_2 = 0.0_dp
  end type troe_parameters

  type, public :: elementary_reaction
    character(len=128) :: equation = ""
    integer :: kind = reaction_kind_elementary
    real(dp), allocatable :: reactant_stoich(:)
    real(dp), allocatable :: product_stoich(:)
    type(arrhenius_rate) :: forward_rate
    type(arrhenius_rate) :: low_pressure_rate
    type(arrhenius_rate) :: high_pressure_rate
    real(dp), allocatable :: third_body_efficiencies(:)
    type(troe_parameters) :: troe
    logical :: reversible = .true.
    ! Explicit kinetic reference-concentration convention [J/(kmol K)].
    ! This does not change Arrhenius activation units or the thermodynamic EOS.
    real(dp) :: equilibrium_reference_gas_constant = universal_gas_constant
    real(dp) :: sri(5) = [1.0_dp,0.0_dp,1.0_dp,1.0_dp,0.0_dp]
    real(dp) :: tsang(2) = [1.0_dp,0.0_dp]
    integer :: plog_count=0
    real(dp) :: plog(4,plog_capacity)=0.0_dp ! pressure [Pa], A, b, Ea [J/kmol]
    integer :: chebyshev_shape(2)=0
    real(dp) :: chebyshev_bounds(4)=0.0_dp ! Tmin, Tmax [K], Pmin, Pmax [Pa]
    real(dp) :: chebyshev(chebyshev_capacity,chebyshev_capacity)=0.0_dp
  end type elementary_reaction

  public :: valid_arrhenius_rate
  public :: valid_troe_parameters
  public :: valid_elementary_reaction
  public :: arrhenius_rate_constant
  public :: reaction_equilibrium_constant
  public :: effective_third_body_concentration
  public :: troe_falloff_factor
  public :: reaction_effective_forward_rate
  public :: reaction_progress_rate
  public :: molar_concentrations_from_mass_fractions
  public :: elementary_production_rates
  public :: elementary_production_rates_from_concentrations
  public :: elementary_production_jacobian
  public :: elementary_mass_fraction_rhs
  public :: elementary_mass_fraction_jacobian

contains

  logical function valid_arrhenius_rate(rate) result(valid)
    type(arrhenius_rate), intent(in) :: rate

    valid = rate%pre_exponential >= 0.0_dp .and. &
      all(ieee_is_finite([rate%pre_exponential, &
        rate%temperature_exponent, rate%activation_energy]))
  end function valid_arrhenius_rate

  logical function valid_troe_parameters(parameters) result(valid)
    type(troe_parameters), intent(in) :: parameters

    if (.not. parameters%enabled) then
      valid = .true.
      return
    end if
    valid = parameters%temperature_3 > 0.0_dp .and. &
      parameters%temperature_1 > 0.0_dp .and. &
      all(ieee_is_finite([parameters%alpha, parameters%temperature_3, &
        parameters%temperature_1, parameters%temperature_2]))
  end function valid_troe_parameters

  logical function valid_elementary_reaction(reaction, nspecies) result(valid)
    type(elementary_reaction), intent(in) :: reaction
    integer, intent(in) :: nspecies

    valid = .false.
    if (nspecies < 1) return
    if (.not.ieee_is_finite(reaction%equilibrium_reference_gas_constant)) return
    if (reaction%equilibrium_reference_gas_constant <= 0._dp) return
    if (.not. allocated(reaction%reactant_stoich)) return
    if (.not. allocated(reaction%product_stoich)) return
    if (size(reaction%reactant_stoich) /= nspecies) return
    if (size(reaction%product_stoich) /= nspecies) return
    if (any(.not. ieee_is_finite(reaction%reactant_stoich))) return
    if (any(.not. ieee_is_finite(reaction%product_stoich))) return
    if (any(reaction%reactant_stoich < 0.0_dp)) return
    if (any(reaction%product_stoich < 0.0_dp)) return
    if (sum(reaction%reactant_stoich) <= 0.0_dp) return
    if (sum(reaction%product_stoich) <= 0.0_dp) return
    if (maxval(abs(reaction%product_stoich - &
        reaction%reactant_stoich)) <= 0.0_dp) return

    select case (reaction%kind)
    case (reaction_kind_chebyshev)
      if (.not.valid_chebyshev_table(reaction)) return
    case (reaction_kind_plog)
      if (.not.valid_plog_table(reaction)) return
    case (reaction_kind_elementary)
      if (.not. valid_arrhenius_rate(reaction%forward_rate)) return
    case (reaction_kind_three_body)
      if (.not. valid_arrhenius_rate(reaction%forward_rate)) return
      if (.not. valid_efficiencies(reaction, nspecies)) return
    case (reaction_kind_falloff,reaction_kind_sri,reaction_kind_activated,reaction_kind_sri_activated, &
        reaction_kind_tsang,reaction_kind_tsang_activated)
      if (.not. valid_arrhenius_rate(reaction%low_pressure_rate)) return
      if (.not. valid_arrhenius_rate(reaction%high_pressure_rate)) return
      if (.not. valid_efficiencies(reaction, nspecies)) return
      if (.not. valid_troe_parameters(reaction%troe)) return
      if (reaction%kind==reaction_kind_tsang.or.reaction%kind==reaction_kind_tsang_activated) then
        if (reaction%troe%enabled.or..not.all(ieee_is_finite(reaction%tsang))) return
        if (all(reaction%tsang<=0)) return
      end if
      if (reaction%kind==reaction_kind_sri.or.reaction%kind==reaction_kind_sri_activated) then
        if (reaction%troe%enabled.or..not.all(ieee_is_finite(reaction%sri))) return
        if (reaction%sri(1)<0.or.reaction%sri(3)<0.or.reaction%sri(4)<=0) return
        if (reaction%sri(1)==0.and.reaction%sri(3)==0) return
      end if
    case default
      return
    end select
    valid = .true.
  end function valid_elementary_reaction

  logical function valid_efficiencies(reaction, nspecies) result(valid)
    type(elementary_reaction), intent(in) :: reaction
    integer, intent(in) :: nspecies

    valid = allocated(reaction%third_body_efficiencies)
    if (.not. valid) return
    valid = size(reaction%third_body_efficiencies) == nspecies
    if (.not. valid) return
    valid = all(ieee_is_finite(reaction%third_body_efficiencies)) .and. &
      all(reaction%third_body_efficiencies >= 0.0_dp)
  end function valid_efficiencies

  logical function valid_chebyshev_table(reaction) result(ok)
    type(elementary_reaction), intent(in) :: reaction
    integer :: nt,np
    ok=.false.
    if (any(reaction%chebyshev_shape<1).or.any(reaction%chebyshev_shape>chebyshev_capacity)) return
    nt=reaction%chebyshev_shape(1);np=reaction%chebyshev_shape(2)
    if (.not.all(ieee_is_finite(reaction%chebyshev_bounds))) return
    if (any(reaction%chebyshev_bounds<=0)) return
    if (reaction%chebyshev_bounds(2)<=reaction%chebyshev_bounds(1)) return
    if (reaction%chebyshev_bounds(4)<=reaction%chebyshev_bounds(3)) return
    if (.not.all(ieee_is_finite(reaction%chebyshev(1:nt,1:np)))) return
    ok=.true.
  end function

  subroutine chebyshev_forward_rate(reaction,temperature,total_concentration,rate,derivative,ok)
    type(elementary_reaction), intent(in) :: reaction
    real(dp), intent(in) :: temperature,total_concentration
    real(dp), intent(out) :: rate,derivative
    logical, intent(out) :: ok
    real(dp) :: tx(chebyshev_capacity),px(chebyshev_capacity),dx(chebyshev_capacity)
    real(dp) :: x,y,logp,span,lograte,slope,bounds(4)
    integer :: i,j,nt,np
    ok=.false.;rate=0;derivative=0
    if (.not.valid_chebyshev_table(reaction)) return
    if (.not.ieee_is_finite(temperature).or..not.ieee_is_finite(total_concentration)) return
    if (temperature<=0.or.total_concentration<=0) return
    bounds=reaction%chebyshev_bounds
    nt=reaction%chebyshev_shape(1);np=reaction%chebyshev_shape(2)
    ! Log pressure avoids forming R*T*C, which can overflow outside a fit range.
    logp=log10(universal_gas_constant)+log10(temperature)+log10(total_concentration)
    span=log10(bounds(4))-log10(bounds(3))
    if (span<=0.or.1/bounds(2)==1/bounds(1)) return
    x=(2/temperature-1/bounds(1)-1/bounds(2))/(1/bounds(2)-1/bounds(1))
    y=(2*logp-log10(bounds(3))-log10(bounds(4)))/span
    tx=0;px=0;dx=0;tx(1)=1;px(1)=1
    tx(2)=x;px(2)=y;dx(2)=1
    do i=3,max(nt,np)
      if (i<=nt) tx(i)=2*x*tx(i-1)-tx(i-2)
      if (i<=np) then
        px(i)=2*y*px(i-1)-px(i-2)
        dx(i)=2*px(i-1)+2*y*dx(i-1)-dx(i-2)
      end if
    end do
    lograte=0;slope=0
    do j=1,np
      do i=1,nt
        lograte=lograte+reaction%chebyshev(i,j)*tx(i)*px(j)
        slope=slope+reaction%chebyshev(i,j)*tx(i)*dx(j)
      end do
    end do
    lograte=lograte*log(10.0_dp)
    if (.not.ieee_is_finite(lograte).or..not.ieee_is_finite(slope)) return
    if (abs(lograte)>logarithm_limit) return
    rate=exp(lograte)
    derivative=2*rate*slope/(total_concentration*span)
    ok=ieee_is_finite(rate).and.ieee_is_finite(derivative)
  end subroutine

  logical function valid_plog_table(reaction) result(ok)
    type(elementary_reaction), intent(in) :: reaction
    integer :: i,n
    ok=.false.;n=reaction%plog_count
    if (n<1.or.n>plog_capacity) return
    if (.not.all(ieee_is_finite(reaction%plog(:,1:n)))) return
    if (any(reaction%plog(1,1:n)<=0)) return
    do i=2,n
      if (reaction%plog(1,i)<reaction%plog(1,i-1)) return
    end do
    ok=.true.
  end function

  subroutine plog_forward_rate(reaction,temperature,total_concentration,rate,derivative,ok)
    type(elementary_reaction), intent(in) :: reaction
    real(dp), intent(in) :: temperature,total_concentration
    real(dp), intent(out) :: rate,derivative
    logical, intent(out) :: ok
    real(dp) :: pressures(plog_capacity),logs(plog_capacity),terms(plog_capacity),weights(plog_capacity)
    real(dp) :: pivot,total,pressure,slope,log_rate,a
    integer :: i,j,first,n,groups,g
    rate=0;derivative=0;ok=.false.
    if (.not.valid_plog_table(reaction)) return
    if (.not.all(ieee_is_finite([temperature,total_concentration]))) return
    if (temperature<=0.or.total_concentration<0) return
    n=reaction%plog_count;groups=0;i=1
    do while(i<=n)
      first=i;groups=groups+1;pressures(groups)=reaction%plog(1,i)
      do
        a=reaction%plog(2,i);weights(i)=0;terms(i)=-huge(1.0_dp)/4
        if (a/=0) then
          weights(i)=sign(1.0_dp,a)
          terms(i)=log(abs(a))+reaction%plog(3,i)*log(temperature)- &
            reaction%plog(4,i)/(universal_gas_constant*temperature)
        end if
        i=i+1
        if (i>n) exit
        if (reaction%plog(1,i)/=pressures(groups)) exit
      end do
      j=i-1;pivot=maxval(terms(first:j))
      total=sum(weights(first:j)*exp(terms(first:j)-pivot))
      if (total<=0.or..not.ieee_is_finite(total)) return
      logs(groups)=pivot+log(total)
    end do
    pressure=universal_gas_constant*temperature*total_concentration
    if (.not.ieee_is_finite(pressure)) return
    slope=0
    if (pressure<=pressures(1)) then
      log_rate=logs(1)
    else if (pressure>=pressures(groups)) then
      log_rate=logs(groups)
    else
      g=1
      do while(g<groups-1)
        if (pressure<pressures(g+1)) exit
        g=g+1
      end do
      slope=(logs(g+1)-logs(g))/(log(pressures(g+1))-log(pressures(g)))
      log_rate=logs(g)+slope*(log(pressure)-log(pressures(g)))
    end if
    if (.not.ieee_is_finite(log_rate).or.abs(log_rate)>logarithm_limit) return
    rate=exp(log_rate)
    if (slope/=0.and.total_concentration>0) derivative=rate*slope/total_concentration
    ok=ieee_is_finite(rate).and.ieee_is_finite(derivative)
  end subroutine

  real(dp) function arrhenius_rate_constant(rate, temperature, ok) &
      result(rate_constant)
    type(arrhenius_rate), intent(in) :: rate
    real(dp), intent(in) :: temperature
    logical, intent(out) :: ok

    real(dp) :: logarithm

    rate_constant = 0.0_dp
    ok = valid_arrhenius_rate(rate) .and. temperature > 0.0_dp .and. &
      ieee_is_finite(temperature)
    if (.not. ok) return
    if (rate%pre_exponential <= tiny(1.0_dp)) return

    logarithm = log(rate%pre_exponential) + &
      rate%temperature_exponent * log(temperature) - &
      rate%activation_energy / (universal_gas_constant * temperature)
    if (.not. ieee_is_finite(logarithm) .or. &
        logarithm > logarithm_limit) then
      ok = .false.
      return
    end if
    if (logarithm < -logarithm_limit) then
      rate_constant = 0.0_dp
    else
      rate_constant = exp(logarithm)
    end if
    ok = ieee_is_finite(rate_constant) .and. rate_constant >= 0.0_dp
  end function arrhenius_rate_constant

  real(dp) function reaction_equilibrium_constant( &
      species,reaction,temperature,ok) result(equilibrium_constant)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reaction
    real(dp), intent(in) :: temperature
    logical, intent(out) :: ok
    real(dp) :: log_kc
    equilibrium_constant=0
    log_kc=reaction_log_equilibrium_constant(species,reaction,temperature,ok)
    if(.not.ok) return
    if(abs(log_kc)>logarithm_limit) then
      ok=.false.;return
    endif
    equilibrium_constant=exp(log_kc)
    ok=ieee_is_finite(equilibrium_constant).and.equilibrium_constant>0
  end function reaction_equilibrium_constant

  ! A finite inverse-equilibrium contribution does not require a
  ! representable standalone Kc. Keep the ordinary division path exactly,
  ! and combine logarithms only when that intermediate is out of range.
  ! Signed Newton products retain their sign; true output overflow refuses.
  subroutine scaled_inverse_equilibrium_action(rate,concentrations,stoich,log_kc,result,ok,gradient_species)
    real(dp), intent(in) :: rate,concentrations(:),stoich(:),log_kc
    real(dp), intent(out) :: result
    logical, intent(out) :: ok
    integer, optional, intent(in) :: gradient_species
    real(dp) :: logarithm,powers(size(stoich)),factor,term,result_sign
    integer :: i,order
    result=0;ok=.false.
    if(size(concentrations)/=size(stoich)) return
    if(.not.all(ieee_is_finite([rate,log_kc])).or.any(.not.ieee_is_finite(concentrations)).or. &
      any(.not.ieee_is_finite(stoich)).or.any(stoich<0)) return
    powers=stoich;factor=1
    if(present(gradient_species)) then
      if(gradient_species<1.or.gradient_species>size(stoich)) return
      factor=stoich(gradient_species)
      if(factor==0) then
        ok=.true.;return
      endif
      powers(gradient_species)=powers(gradient_species)-1
    endif
    ! Validate signs and fractional zero-coordinate derivatives before an
    ! absent factor can return zero. Signed integer mass-action polynomials
    ! keep their existing opt-in behavior; fractional negative powers refuse.
    result_sign=sign(1._dp,rate)
    do i=1,size(stoich)
      if(concentrations(i)==0.and.powers(i)<0) return
      if(concentrations(i)>=0.or.powers(i)==0) cycle
      if(abs(powers(i))>real(huge(order),dp)) return
      order=nint(powers(i))
      if(real(order,dp)/=powers(i)) return
      if(mod(abs(order),2)==1) result_sign=-result_sign
    enddo
    if(rate==0.or.any(concentrations==0.and.powers>0)) then
      ok=.true.;return
    endif
    ! Combine the coefficient and activities before exponentiation. Neither
    ! an overflowing reciprocal Kc nor an underflowing concentration product
    ! is required to exist separately when the actual contribution is finite.
    logarithm=log(abs(rate))+log(factor)-log_kc
    do i=1,size(stoich)
      if(powers(i)==0) cycle
      term=log(abs(concentrations(i)))
      if(abs(term)>1) then
        if(abs(powers(i))>huge(1._dp)/abs(term)) return
      endif
      term=term*powers(i)
      if(term>0) then
        if(logarithm>huge(1._dp)-term) return
      else if(term<0) then
        if(logarithm< -huge(1._dp)-term) return
      endif
      logarithm=logarithm+term
    enddo
    if(.not.ieee_is_finite(logarithm).or.logarithm>log(huge(1._dp))) return
    result=result_sign*exp(logarithm)
    ok=ieee_is_finite(result)
  end subroutine scaled_inverse_equilibrium_action

  subroutine subtract_inverse_equilibrium_term(value,term,ok)
    real(dp), intent(inout) :: value
    real(dp), intent(in) :: term
    logical, intent(out) :: ok
    ok=.false.
    if(.not.all(ieee_is_finite([value,term]))) return
    if(term<0) then
      if(value>huge(1._dp)+term) return
    else if(term>0) then
      if(value< -huge(1._dp)+term) return
    endif
    value=value-term;ok=ieee_is_finite(value)
  end subroutine subtract_inverse_equilibrium_term

  real(dp) function reaction_log_equilibrium_constant( &
      species, reaction, temperature, ok) result(log_kc)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reaction
    real(dp), intent(in) :: temperature
    logical, intent(out) :: ok

    real(dp) :: cp_over_r, h_over_rt, s_over_r
    real(dp) :: delta_g_over_rt, delta_nu
    real(dp) :: reference_concentration, net_stoich
    logical :: species_ok
    integer :: i

    log_kc = 0.0_dp
    ok = valid_elementary_reaction(reaction, size(species)) .and. &
      temperature > 0.0_dp
    if (.not. ok) return

    delta_g_over_rt = 0.0_dp
    delta_nu = 0.0_dp
    do i = 1, size(species)
      call nasa7_dimensionless_properties( &
        species(i), temperature, cp_over_r, h_over_rt, s_over_r, species_ok)
      if (.not. species_ok) then
        ok = .false.
        return
      end if
      net_stoich = reaction%product_stoich(i) - &
        reaction%reactant_stoich(i)
      delta_g_over_rt = delta_g_over_rt + &
        net_stoich * (h_over_rt - s_over_r)
      delta_nu = delta_nu + net_stoich
    end do

    ! Reject unrepresentable reference concentrations before arithmetic traps.
    ok = .false.
    if (temperature > 1._dp) then
      if (reaction%equilibrium_reference_gas_constant > huge(1._dp)/temperature) return
    endif
    if (reaction%equilibrium_reference_gas_constant < &
        (kinetics_standard_pressure/huge(1._dp))/temperature) return
    reference_concentration = kinetics_standard_pressure / &
      (reaction%equilibrium_reference_gas_constant * temperature)
    log_kc = -delta_g_over_rt + delta_nu * log(reference_concentration)
    ok=ieee_is_finite(log_kc)
  end function reaction_log_equilibrium_constant

  real(dp) function effective_third_body_concentration( &
      reaction, concentrations, ok,signed_concentrations) result(effective_concentration)
    type(elementary_reaction), intent(in) :: reaction
    real(dp), intent(in) :: concentrations(:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: signed_concentrations
    logical :: signed_trial

    signed_trial=.false.
    if(present(signed_concentrations)) signed_trial=signed_concentrations
    effective_concentration = 0.0_dp
    ok = allocated(reaction%third_body_efficiencies)
    if (.not. ok) return
    ok = size(reaction%third_body_efficiencies) == size(concentrations) .and. &
      all(ieee_is_finite(concentrations)) .and. &
      (signed_trial.or.all(concentrations >= 0.0_dp))
    if (.not. ok) return
    effective_concentration = sum( &
      reaction%third_body_efficiencies * concentrations)
    ok = ieee_is_finite(effective_concentration) .and. &
      (signed_trial.or.effective_concentration >= 0.0_dp)
  end function effective_third_body_concentration

  ! Cantera 3.2 SRI form: d*T**e*(a*exp(-b/T)+exp(-T/c))**X.
  real(dp) function sri_falloff_factor(parameters,temperature,reduced_pressure,ok,derivative) result(factor)
    real(dp), intent(in) :: parameters(5),temperature,reduced_pressure
    logical, intent(out) :: ok
    real(dp), intent(out) :: derivative
    real(dp) :: terms(2),center,lp,x,lf,pivot
    factor=0;derivative=0;ok=.false.
    if (.not.all(ieee_is_finite([parameters,temperature,reduced_pressure]))) return
    if (temperature<=0.or.reduced_pressure<0) return
    if (parameters(1)<0.or.parameters(3)<0.or.parameters(4)<=0) return
    if (parameters(1)==0.and.parameters(3)==0) return
    ! Log-sum-exp avoids underflow of the temperature-dependent center.
    terms=-huge(1.0_dp)/4
    if (parameters(1)>0) terms(1)=log(parameters(1))-parameters(2)/temperature
    if (parameters(3)>0) terms(2)=-temperature/parameters(3)
    pivot=maxval(terms);center=pivot+log(sum(exp(terms-pivot)))
    lp=log10(max(reduced_pressure,minimum_reduced_pressure));x=1/(1+lp*lp)
    lf=log(parameters(4))+parameters(5)*log(temperature)+x*center
    if (.not.ieee_is_finite(lf).or.abs(lf)>logarithm_limit) return
    factor=exp(lf)
    derivative=-2*lp*x*x*center/log(10.0_dp)
    ok=ieee_is_finite(factor).and.ieee_is_finite(derivative)
  end function sri_falloff_factor

  real(dp) function troe_falloff_factor( &
      parameters, temperature, reduced_pressure, ok, &
      logarithmic_derivative,center_override) result(falloff_factor)
    type(troe_parameters), intent(in) :: parameters
    real(dp), intent(in), optional :: center_override
    real(dp), intent(in) :: temperature, reduced_pressure
    logical, intent(out) :: ok
    real(dp), intent(out), optional :: logarithmic_derivative

    real(dp) :: f_center, log_f_center, log_pressure
    real(dp) :: c_parameter, n_parameter, numerator, denominator
    real(dp) :: broadening_argument, log_factor
    real(dp) :: derivative_argument, derivative_log_factor

    falloff_factor = 1.0_dp
    if (present(logarithmic_derivative)) logarithmic_derivative = 0.0_dp
    ok = temperature > 0.0_dp .and. reduced_pressure >= 0.0_dp .and. &
      ieee_is_finite(temperature) .and. ieee_is_finite(reduced_pressure) .and. &
      valid_troe_parameters(parameters)
    if (.not. ok) return
    if (.not.parameters%enabled.and..not.present(center_override)) return

    if (present(center_override)) then
      f_center=center_override
    else
    f_center = (1.0_dp - parameters%alpha) * &
      exp(-temperature / parameters%temperature_3) + &
      parameters%alpha * exp(-temperature / parameters%temperature_1)
    if (parameters%temperature_2 /= 0.0_dp) then
      if (-parameters%temperature_2/temperature > logarithm_limit) then
        ok=.false.
        return
      end if
      f_center = f_center + exp(-parameters%temperature_2 / temperature)
    end if
    end if
    if (f_center <= 0.0_dp .or. .not. ieee_is_finite(f_center)) then
      ok = .false.
      return
    end if

    log_f_center = log10(f_center)
    log_pressure = log10(max(reduced_pressure, minimum_reduced_pressure))
    c_parameter = -0.4_dp - 0.67_dp * log_f_center
    n_parameter = 0.75_dp - 1.27_dp * log_f_center
    numerator = log_pressure + c_parameter
    denominator = n_parameter - 0.14_dp * numerator
    if (abs(denominator) <= tiny(1.0_dp)) then
      ok = .false.
      return
    end if
    broadening_argument = numerator / denominator
    log_factor = log_f_center / &
      (1.0_dp + broadening_argument * broadening_argument)
    falloff_factor = 10.0_dp**log_factor

    if (present(logarithmic_derivative)) then
      derivative_argument = n_parameter / (denominator * denominator)
      derivative_log_factor = -2.0_dp * log_f_center * &
        broadening_argument * derivative_argument / &
        (1.0_dp + broadening_argument * broadening_argument)**2
      logarithmic_derivative = derivative_log_factor
    end if
    ok = ieee_is_finite(falloff_factor) .and. falloff_factor > 0.0_dp
  end function troe_falloff_factor

  subroutine reaction_effective_forward_rate( &
      reaction, temperature, concentrations, rate_constant, &
      concentration_derivative, ok, &
      signed_concentrations)
    type(elementary_reaction), intent(in) :: reaction
    real(dp), intent(in) :: temperature, concentrations(:)
    real(dp), intent(out) :: rate_constant
    real(dp), intent(out), optional :: concentration_derivative(:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: signed_concentrations
    logical :: signed_trial

    real(dp) :: base_rate, low_rate, high_rate
    real(dp) :: third_body, reduced_pressure, falloff_factor
    real(dp) :: log_factor_derivative, derivative_wrt_pressure
    real(dp) :: derivative_wrt_third_body
    integer :: i

    signed_trial=.false.
    if(present(signed_concentrations)) signed_trial=signed_concentrations
    rate_constant = 0.0_dp
    if (present(concentration_derivative)) concentration_derivative = 0.0_dp
    ok = size(concentrations) > 0 .and. &
      all(ieee_is_finite(concentrations)) .and. (signed_trial.or.all(concentrations >= 0.0_dp))
    if (present(concentration_derivative)) then
      ok = ok .and. size(concentration_derivative) == size(concentrations)
    end if
    if (.not. ok) return

    select case (reaction%kind)
    case (reaction_kind_chebyshev)
      call chebyshev_forward_rate(reaction,temperature,sum(concentrations),rate_constant,derivative_wrt_pressure,ok)
      if (present(concentration_derivative)) concentration_derivative=derivative_wrt_pressure
    case (reaction_kind_plog)
      call plog_forward_rate(reaction,temperature,sum(concentrations),rate_constant,derivative_wrt_pressure,ok)
      if (present(concentration_derivative)) concentration_derivative=derivative_wrt_pressure
    case (reaction_kind_elementary)
      rate_constant = arrhenius_rate_constant( &
        reaction%forward_rate, temperature, ok)
    case (reaction_kind_three_body)
      if (.not. valid_efficiencies(reaction, size(concentrations))) then
        ok = .false.
        return
      end if
      base_rate = arrhenius_rate_constant( &
        reaction%forward_rate, temperature, ok)
      if (.not. ok) return
      third_body = effective_third_body_concentration( &
        reaction, concentrations, ok,signed_concentrations=signed_trial)
      if (.not. ok) return
      rate_constant = base_rate * third_body
      if (present(concentration_derivative)) then
        concentration_derivative = base_rate * &
          reaction%third_body_efficiencies
      end if
    case (reaction_kind_falloff,reaction_kind_sri,reaction_kind_activated,reaction_kind_sri_activated, &
        reaction_kind_tsang,reaction_kind_tsang_activated)
      if (.not. valid_efficiencies(reaction, size(concentrations))) then
        ok = .false.
        return
      end if
      low_rate = arrhenius_rate_constant( &
        reaction%low_pressure_rate, temperature, ok)
      if (.not. ok) return
      high_rate = arrhenius_rate_constant( &
        reaction%high_pressure_rate, temperature, ok)
      if (.not. ok .or. high_rate <= 0.0_dp) then
        ok = .false.
        return
      end if
      third_body = effective_third_body_concentration( &
        reaction, concentrations, ok,signed_concentrations=signed_trial)
      if (.not. ok) return
      reduced_pressure = low_rate * third_body / high_rate
      if (reaction%kind==reaction_kind_sri.or.reaction%kind==reaction_kind_sri_activated) then
        falloff_factor=sri_falloff_factor(reaction%sri,temperature,reduced_pressure,ok,log_factor_derivative)
      else if (reaction%kind==reaction_kind_tsang.or.reaction%kind==reaction_kind_tsang_activated) then
        falloff_factor=troe_falloff_factor(reaction%troe,temperature,reduced_pressure,ok, &
          log_factor_derivative,reaction%tsang(1)+reaction%tsang(2)*temperature)
      else
      falloff_factor = troe_falloff_factor( &
        reaction%troe, temperature, reduced_pressure, ok, &
        log_factor_derivative)
      end if
      if (.not. ok) return
      if (reaction%kind==reaction_kind_activated.or.reaction%kind==reaction_kind_sri_activated.or. &
          reaction%kind==reaction_kind_tsang_activated) then
        rate_constant=low_rate/(1+reduced_pressure)*falloff_factor
      else
        rate_constant = high_rate * reduced_pressure / &
          (1.0_dp + reduced_pressure) * falloff_factor
      end if
      if (present(concentration_derivative)) then
        if (reaction%kind==reaction_kind_activated.or.reaction%kind==reaction_kind_sri_activated.or. &
          reaction%kind==reaction_kind_tsang_activated) then
          derivative_wrt_third_body=-rate_constant/(1+reduced_pressure)*(low_rate/high_rate)
          if (reduced_pressure>minimum_reduced_pressure) &
            derivative_wrt_third_body=derivative_wrt_third_body+rate_constant*log_factor_derivative/third_body
        else if (third_body <= tiny(1.0_dp)) then
          derivative_wrt_third_body = low_rate * falloff_factor
        else if (reduced_pressure > 0.5_dp*sqrt(huge(1.0_dp))) then
          ! d(k)/dM = k/M * [1/(1+Pr) + d(log F)/d(log Pr)].
          ! Avoid overflowing (1+Pr)**2 when the derivative is finite.
          ! Keep the ordinary-range arithmetic unchanged.
          derivative_wrt_third_body = rate_constant / third_body * &
            (1.0_dp / (1.0_dp + reduced_pressure) + log_factor_derivative)
        else
          derivative_wrt_pressure = high_rate * falloff_factor * &
            (1.0_dp / (1.0_dp + reduced_pressure)**2 + &
             log_factor_derivative / (1.0_dp + reduced_pressure))
          derivative_wrt_third_body = derivative_wrt_pressure * &
            low_rate / high_rate
        end if
        do i = 1, size(concentrations)
          concentration_derivative(i) = derivative_wrt_third_body * &
            reaction%third_body_efficiencies(i)
        end do
      end if
    case default
      ok = .false.
    end select
    ok = ok .and. ieee_is_finite(rate_constant) .and. &
      (rate_constant >= 0.0_dp.or.(signed_trial.and.reaction%kind==reaction_kind_three_body))
    if (present(concentration_derivative)) then
      ok = ok .and. all(ieee_is_finite(concentration_derivative))
    end if
  end subroutine reaction_effective_forward_rate

  subroutine reaction_progress_rate( &
      species, reaction, temperature, concentrations, forward_progress, &
      reverse_progress, net_progress, ok, &
      signed_concentrations)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reaction
    real(dp), intent(in) :: temperature, concentrations(:)
    real(dp), intent(out) :: forward_progress, reverse_progress, net_progress
    logical, intent(out) :: ok
    logical, optional, intent(in) :: signed_concentrations
    logical :: signed_trial

    real(dp) :: forward_constant, reverse_constant, equilibrium_constant,log_kc
    real(dp) :: reactant_product, product_product
    real(dp), allocatable :: unused_derivative(:)

    signed_trial=.false.
    if(present(signed_concentrations)) signed_trial=signed_concentrations
    forward_progress = 0.0_dp
    reverse_progress = 0.0_dp
    net_progress = 0.0_dp
    ok = valid_elementary_reaction(reaction, size(species)) .and. &
      size(concentrations) == size(species) .and. &
      all(ieee_is_finite(concentrations)) .and. &
      (signed_trial.or.all(concentrations >= 0.0_dp))
    if (.not. ok) return

    allocate(unused_derivative(size(species)))
    call reaction_effective_forward_rate( &
      reaction, temperature, concentrations, forward_constant, &
      unused_derivative, ok,signed_concentrations=signed_trial)
    if (.not. ok) return
    call concentration_product( &
      concentrations, reaction%reactant_stoich, reactant_product, ok,signed_concentrations=signed_trial)
    if (.not. ok) return
    forward_progress = forward_constant * reactant_product

    if (reaction%reversible) then
      log_kc = reaction_log_equilibrium_constant(species,reaction,temperature,ok)
      if (.not. ok) return
      if(abs(log_kc)<=logarithm_limit) then
        call concentration_product( &
          concentrations, reaction%product_stoich, product_product, ok,signed_concentrations=signed_trial)
        if (.not. ok) return
        equilibrium_constant=exp(log_kc)
        reverse_constant=forward_constant/equilibrium_constant
        reverse_progress=reverse_constant*product_product
      else
        call scaled_inverse_equilibrium_action(forward_constant,concentrations,reaction%product_stoich, &
          log_kc,reverse_progress,ok)
        if(.not.ok) return
      endif
    end if

    net_progress = forward_progress - reverse_progress
    ok = all(ieee_is_finite( &
      [forward_progress, reverse_progress, net_progress]))
  end subroutine reaction_progress_rate

  subroutine concentration_product(concentrations, stoich, product, ok, &
      signed_concentrations)
    real(dp), intent(in) :: concentrations(:), stoich(:)
    real(dp), intent(out) :: product
    logical, intent(out) :: ok
    logical, optional, intent(in) :: signed_concentrations
    logical :: signed_trial
    integer :: i,order

    signed_trial=.false.
    if(present(signed_concentrations)) signed_trial=signed_concentrations
    product = 1.0_dp
    ok = size(concentrations) == size(stoich) .and. &
      all(ieee_is_finite(concentrations)) .and. &
      all(ieee_is_finite(stoich)) .and. &
      (signed_trial.or.all(concentrations >= 0.0_dp)) .and. all(stoich >= 0.0_dp)
    if (.not. ok) return
    do i = 1, size(concentrations)
      if(stoich(i)<=0) cycle
      if(concentrations(i)<0) then
        ! A real exponent on a negative activity can trap. Integer reaction
        ! orders retain their signed polynomial; noninteger orders refuse.
        ok=stoich(i)<=real(huge(order),dp)
        if(.not.ok) return
        order=nint(stoich(i))
        ok=real(order,dp)==stoich(i)
        if(.not.ok) return
        product=product*concentrations(i)**order
      else
        product=product*concentrations(i)**stoich(i)
      endif
    end do
    ok = ieee_is_finite(product) .and. (signed_trial.or.product >= 0.0_dp)
  end subroutine concentration_product

  subroutine concentration_product_gradient( &
      concentrations, stoich, product, gradient, ok, &
      signed_concentrations)
    real(dp), intent(in) :: concentrations(:), stoich(:)
    real(dp), intent(out) :: product, gradient(:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: signed_concentrations
    logical :: signed_trial
    real(dp) :: factor,powers(size(stoich))
    integer :: i, j

    signed_trial=.false.
    if(present(signed_concentrations)) signed_trial=signed_concentrations
    gradient = 0.0_dp
    call concentration_product(concentrations, stoich, product, ok,signed_concentrations=signed_trial)
    if (.not. ok .or. size(gradient) /= size(concentrations)) then
      ok = .false.
      return
    end if

    do j = 1, size(concentrations)
      if (stoich(j) <= 0.0_dp) cycle
      if (concentrations(j) > 0.0_dp.or.(signed_trial.and.concentrations(j)<0)) then
        gradient(j) = product * stoich(j) / concentrations(j)
      else
        if (stoich(j) < 1.0_dp) then
          ok = .false.
          return
        else if (abs(stoich(j) - 1.0_dp) <= epsilon(1.0_dp)) then
          if(signed_trial) then
            powers=stoich;powers(j)=0
            call concentration_product(concentrations,powers,factor,ok,signed_concentrations=.true.)
            if(.not.ok) return
          else
          factor = 1.0_dp
          do i = 1, size(concentrations)
            if (i == j .or. stoich(i) <= 0.0_dp) cycle
            factor = factor * concentrations(i)**stoich(i)
          end do
          endif
          gradient(j) = factor
        else
          gradient(j) = 0.0_dp
        end if
      end if
    end do
    ok = all(ieee_is_finite(gradient))
  end subroutine concentration_product_gradient

  subroutine molar_concentrations_from_mass_fractions( &
      species, density, mass_fractions, concentrations, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: density, mass_fractions(:)
    real(dp), intent(out) :: concentrations(:)
    logical, intent(out) :: ok

    integer :: i

    concentrations = 0.0_dp
    ok = density > 0.0_dp .and. size(concentrations) == size(species)
    if (.not. ok) return
    ok = valid_mixture_composition(species, mass_fractions)
    if (.not. ok) return
    do i = 1, size(species)
      concentrations(i) = density * max(0.0_dp, mass_fractions(i)) / &
        species(i)%molecular_weight
    end do
    ok = all(ieee_is_finite(concentrations)) .and. &
      all(concentrations >= 0.0_dp)
  end subroutine molar_concentrations_from_mass_fractions

  subroutine elementary_production_rates( &
      species, reactions, temperature, density, mass_fractions, &
      molar_production_rates, ok, forward_progress, reverse_progress)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: temperature, density, mass_fractions(:)
    real(dp), intent(out) :: molar_production_rates(:)
    logical, intent(out) :: ok
    real(dp), intent(out), optional :: forward_progress(:), reverse_progress(:)

    real(dp), allocatable :: concentrations(:)

    molar_production_rates = 0.0_dp
    ok = size(molar_production_rates) == size(species) .and. &
      size(reactions) >= 1
    if (.not. ok) return
    allocate(concentrations(size(species)))
    call molar_concentrations_from_mass_fractions( &
      species, density, mass_fractions, concentrations, ok)
    if (.not. ok) return
    call elementary_production_rates_from_concentrations( &
      species, reactions, temperature, concentrations, &
      molar_production_rates, ok, forward_progress, reverse_progress)
  end subroutine elementary_production_rates

  subroutine elementary_production_rates_from_concentrations( &
      species, reactions, temperature, concentrations, &
      molar_production_rates, ok, forward_progress, reverse_progress, &
      signed_concentrations)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: temperature, concentrations(:)
    real(dp), intent(out) :: molar_production_rates(:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: signed_concentrations
    logical :: signed_trial
    real(dp), intent(out), optional :: forward_progress(:), reverse_progress(:)

    real(dp) :: q_forward, q_reverse, q_net
    logical :: reaction_ok
    integer :: i, reaction_index

    signed_trial=.false.
    if(present(signed_concentrations)) signed_trial=signed_concentrations
    molar_production_rates = 0.0_dp
    ok = size(molar_production_rates) == size(species) .and. &
      size(concentrations) == size(species) .and. size(reactions) >= 1 .and. &
      all(ieee_is_finite(concentrations)) .and. &
      (signed_trial.or.all(concentrations >= 0.0_dp))
    if (.not. ok) return
    if (present(forward_progress)) then
      if (size(forward_progress) /= size(reactions)) then
        ok = .false.
        return
      end if
      forward_progress = 0.0_dp
    end if
    if (present(reverse_progress)) then
      if (size(reverse_progress) /= size(reactions)) then
        ok = .false.
        return
      end if
      reverse_progress = 0.0_dp
    end if

    do reaction_index = 1, size(reactions)
      call reaction_progress_rate( &
        species, reactions(reaction_index), temperature, concentrations, &
        q_forward, q_reverse, q_net, reaction_ok,signed_concentrations=signed_trial)
      if (.not. reaction_ok) then
        ok = .false.
        return
      end if
      if (present(forward_progress)) &
        forward_progress(reaction_index) = q_forward
      if (present(reverse_progress)) &
        reverse_progress(reaction_index) = q_reverse
      do i = 1, size(species)
        molar_production_rates(i) = molar_production_rates(i) + &
          (reactions(reaction_index)%product_stoich(i) - &
           reactions(reaction_index)%reactant_stoich(i)) * q_net
      end do
    end do
    ok = all(ieee_is_finite(molar_production_rates))
  end subroutine elementary_production_rates_from_concentrations

  subroutine elementary_production_jacobian( &
      species, reactions, temperature, concentrations, jacobian, ok, &
      signed_concentrations)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: temperature, concentrations(:)
    real(dp), intent(out) :: jacobian(:, :)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: signed_concentrations
    logical :: signed_trial

    real(dp), allocatable :: rate_derivative(:)
    real(dp), allocatable :: reactant_gradient(:), product_gradient(:)
    real(dp), allocatable :: progress_derivative(:)
    real(dp) :: forward_constant, reverse_constant, equilibrium_constant,log_kc
    real(dp) :: reactant_product, product_product, net_stoich,reverse_term
    logical :: reaction_ok
    integer :: i, j, reaction_index, nspecies

    signed_trial=.false.
    if(present(signed_concentrations)) signed_trial=signed_concentrations
    jacobian = 0.0_dp
    nspecies = size(species)
    ok = size(concentrations) == nspecies .and. &
      size(jacobian, 1) == nspecies .and. size(jacobian, 2) == nspecies .and. &
      all(ieee_is_finite(concentrations)) .and. (signed_trial.or.all(concentrations >= 0.0_dp))
    if (.not. ok) return
    allocate(rate_derivative(nspecies), reactant_gradient(nspecies))
    allocate(product_gradient(nspecies), progress_derivative(nspecies))

    do reaction_index = 1, size(reactions)
      if (.not. valid_elementary_reaction( &
          reactions(reaction_index), nspecies)) then
        ok = .false.
        return
      end if
      call reaction_effective_forward_rate( &
        reactions(reaction_index), temperature, concentrations, &
        forward_constant, rate_derivative, reaction_ok,signed_concentrations=signed_trial)
      if (.not. reaction_ok) then
        ok = .false.
        return
      end if
      call concentration_product_gradient( &
        concentrations, reactions(reaction_index)%reactant_stoich, &
        reactant_product, reactant_gradient, reaction_ok,signed_concentrations=signed_trial)
      if (.not. reaction_ok) then
        ok = .false.
        return
      end if
      progress_derivative = rate_derivative * reactant_product + &
        forward_constant * reactant_gradient

      if (reactions(reaction_index)%reversible) then
        log_kc = reaction_log_equilibrium_constant(species,reactions(reaction_index),temperature,reaction_ok)
        if (.not. reaction_ok) then
          ok = .false.
          return
        end if
        if(abs(log_kc)<=logarithm_limit) then
          call concentration_product_gradient( &
            concentrations, reactions(reaction_index)%product_stoich, &
            product_product, product_gradient, reaction_ok,signed_concentrations=signed_trial)
          if (.not. reaction_ok) then
            ok = .false.
            return
          end if
          equilibrium_constant=exp(log_kc)
          reverse_constant=forward_constant/equilibrium_constant
          progress_derivative = progress_derivative - &
            rate_derivative * product_product / equilibrium_constant - &
            reverse_constant * product_gradient
        else
          do j=1,nspecies
            call scaled_inverse_equilibrium_action(rate_derivative(j),concentrations, &
              reactions(reaction_index)%product_stoich,log_kc,reverse_term,reaction_ok)
            if(reaction_ok) call subtract_inverse_equilibrium_term(progress_derivative(j),reverse_term,reaction_ok)
            if(.not.reaction_ok) then
              ok=.false.;return
            endif
            call scaled_inverse_equilibrium_action(forward_constant,concentrations,reactions(reaction_index)%product_stoich, &
              log_kc,reverse_term,reaction_ok,gradient_species=j)
            if(reaction_ok) call subtract_inverse_equilibrium_term(progress_derivative(j),reverse_term,reaction_ok)
            if(.not.reaction_ok) then
              ok=.false.;return
            endif
          enddo
        endif
      end if

      do i = 1, nspecies
        net_stoich = reactions(reaction_index)%product_stoich(i) - &
          reactions(reaction_index)%reactant_stoich(i)
        if (abs(net_stoich) <= tiny(1.0_dp)) cycle
        do j = 1, nspecies
          jacobian(i, j) = jacobian(i, j) + &
            net_stoich * progress_derivative(j)
        end do
      end do
    end do
    ok = all(ieee_is_finite(jacobian))
  end subroutine elementary_production_jacobian

  subroutine elementary_mass_fraction_rhs( &
      species, reactions, temperature, density, mass_fractions, derivative, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: temperature, density, mass_fractions(:)
    real(dp), intent(out) :: derivative(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: molar_production_rates(:)
    integer :: i

    derivative = 0.0_dp
    ok = size(derivative) == size(species) .and. density > 0.0_dp
    if (.not. ok) return
    allocate(molar_production_rates(size(species)))
    call elementary_production_rates( &
      species, reactions, temperature, density, mass_fractions, &
      molar_production_rates, ok)
    if (.not. ok) return
    do i = 1, size(species)
      derivative(i) = species(i)%molecular_weight * &
        molar_production_rates(i) / density
    end do
    ok = all(ieee_is_finite(derivative))
  end subroutine elementary_mass_fraction_rhs

  subroutine elementary_mass_fraction_jacobian( &
      species, reactions, temperature, density, mass_fractions, jacobian, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: temperature, density, mass_fractions(:)
    real(dp), intent(out) :: jacobian(:, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: concentrations(:), molar_jacobian(:, :)
    integer :: i, j, nspecies

    jacobian = 0.0_dp
    nspecies = size(species)
    ok = density > 0.0_dp .and. size(jacobian, 1) == nspecies .and. &
      size(jacobian, 2) == nspecies
    if (.not. ok) return
    allocate(concentrations(nspecies), molar_jacobian(nspecies, nspecies))
    call molar_concentrations_from_mass_fractions( &
      species, density, mass_fractions, concentrations, ok)
    if (.not. ok) return
    call elementary_production_jacobian( &
      species, reactions, temperature, concentrations, molar_jacobian, ok)
    if (.not. ok) return
    do i = 1, nspecies
      do j = 1, nspecies
        jacobian(i, j) = species(i)%molecular_weight / &
          species(j)%molecular_weight * molar_jacobian(i, j)
      end do
    end do
    ok = all(ieee_is_finite(jacobian))
  end subroutine elementary_mass_fraction_jacobian

end module elementary_kinetics_mod
