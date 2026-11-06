module constant_volume_reactor_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species, nasa7_mass_properties,nasa7_eos_bounds
  use mixture_thermo_mod, only: &
    valid_mixture_composition, mixture_mass_properties, &
    temperature_from_internal_energy
  use elementary_kinetics_mod, only: &
    elementary_reaction, elementary_mass_fraction_rhs, valid_elementary_reaction, &
    elementary_mass_fraction_jacobian,elementary_production_rates_from_concentrations,elementary_production_jacobian
  implicit none
  private

  integer, parameter, public :: reactor_max_step_attempts = 24

  public :: reactor_specific_internal_energy
  public :: reaction_source_heat_release
  public :: reactor_rhs
  public :: forced_reactor_reduced_jacobian
  public :: forced_reactor_rhs
  public :: advance_forced_reactor_implicit_adaptive
  public :: advance_forced_reactor_interval
  public :: advance_constant_volume_adaptive
  public :: reactor_reduced_jacobian
  public :: backward_euler_trial
  public :: advance_constant_volume_implicit_adaptive

  public :: forced_reactor_integrator
  abstract interface
    subroutine forced_reactor_integrator(species,reactions,initial_energy_density,external_species, &
        external_energy,duration,rtol,atol,species_density,temperature,completed_steps,ok,maximum_steps)
      import :: dp,nasa7_species,elementary_reaction
      type(nasa7_species), intent(in) :: species(:)
      type(elementary_reaction), intent(in) :: reactions(:)
      real(dp), intent(in) :: initial_energy_density,external_species(:),external_energy,duration,rtol,atol
      real(dp), intent(inout) :: species_density(:),temperature
      integer, intent(out) :: completed_steps
      logical, intent(out) :: ok
      integer, optional, intent(in) :: maximum_steps
    end subroutine
  end interface
contains

  ! PeleC React.cpp stores -sum(h_k I_R_k), using accepted-time species
  ! mass enthalpies and chemical rhoY rates (external forcing excluded).
  ! SI output is W/m^3. Endothermic rates are valid negative outputs.
  subroutine reaction_source_heat_release(species,temperature,species_rates,heat_release,ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: temperature,species_rates(:)
    real(dp), intent(inout) :: heat_release
    logical, intent(out) :: ok
    real(dp) :: cp,cv,enthalpy,energy,entropy,term,candidate
    logical :: valid
    integer :: k
    ok=.false.;candidate=0
    if(size(species)<1.or.size(species_rates)/=size(species)) return
    if(.not.ieee_is_finite(temperature).or.any(.not.ieee_is_finite(species_rates))) return
    do k=1,size(species)
      call nasa7_mass_properties(species(k),temperature,cp,cv,enthalpy,energy,entropy,valid)
      if(.not.valid) return
      if(abs(enthalpy)>1._dp) then
        if(abs(species_rates(k))>huge(1._dp)/abs(enthalpy)) return
      endif
      term=enthalpy*species_rates(k)
      if(term>0) then
        if(candidate< -huge(1._dp)+term) return
      elseif(term<0) then
        if(candidate>huge(1._dp)+term) return
      endif
      candidate=candidate-term
    enddo
    if(.not.ieee_is_finite(candidate)) return
    heat_release=candidate;ok=.true.
  end subroutine reaction_source_heat_release

  ! Reactor state for SDC: rhoY evolves with chemistry plus fixed external
  ! species rates; rhoe(t)=rhoe(0)+t*external_energy. Density follows sum(rhoY).
  subroutine forced_reactor_rhs(species,reactions,initial_energy_density,elapsed,species_density, &
      external_species,external_energy,temperature_guess,derivative,temperature,ok, &
      permit_signed_trials)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: initial_energy_density,elapsed,species_density(:),external_species(:),external_energy
    real(dp), intent(in) :: temperature_guess
    real(dp), intent(out) :: derivative(:),temperature
    logical, intent(out) :: ok
    logical, optional, intent(in) :: permit_signed_trials
    real(dp) :: density,energy,candidate_temperature,y(size(species)),rates(size(species)),candidate(size(species))
    logical :: valid,signed_trial
    signed_trial=.false.
    if(present(permit_signed_trials)) signed_trial=permit_signed_trials
    derivative=0;temperature=0;ok=.false.
    if(size(species)<1.or.size(species_density)/=size(species).or. &
      size(external_species)/=size(species).or.size(derivative)/=size(species)) return
    if(.not.all(ieee_is_finite([initial_energy_density,elapsed,external_energy,temperature_guess]))) return
    if(.not.all(ieee_is_finite(species_density)).or..not.all(ieee_is_finite(external_species))) return
    if(elapsed<0) return
    if(.not.signed_trial.and.any(species_density<0)) return
    density=sum(species_density)
    if(.not.ieee_is_finite(density).or.density<=0) return
    y=species_density/density
    energy=(initial_energy_density+elapsed*external_energy)/density
    if(.not.ieee_is_finite(energy)) return
    if(size(reactions)==0) then
      ! An empty mechanism is the inert limit of the forced reactor.
      rates=0
      call recover_trial_temperature(species,y,energy,temperature_guess,candidate_temperature,valid,signed_trial)
    else
      call reactor_rhs(species,reactions,density,energy,y,temperature_guess,rates,candidate_temperature,valid, &
        permit_signed_trials=signed_trial)
    end if
    if(.not.valid) return
    candidate=density*rates+external_species
    if(.not.all(ieee_is_finite(candidate))) return
    derivative=candidate;temperature=candidate_temperature;ok=.true.
  end subroutine forced_reactor_rhs

  ! Accept one adaptive substep. The caller retains rhoe(0), advances elapsed
  ! by accepted_time_step, and repeats until its requested interval is covered.
  ! absolute_tolerance is in species-density units, not mass-fraction units.
  ! At fixed time, the forced total density and internal-energy density are
  ! analytic constants. For independent rhoY coordinates, the factors rho in
  ! rho*dY/dt and 1/rho in dY/d(rhoY) cancel. The reduced fixed-energy Jacobian
  ! therefore applies directly, including the dependent-species closure.
  subroutine forced_reactor_reduced_jacobian(species,reactions,initial_energy_density,elapsed,species_density, &
      external_species,external_energy,temperature_guess,jacobian,temperature,ok,dependent_species,permit_signed_trials)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: initial_energy_density,elapsed,species_density(:),external_species(:),external_energy
    real(dp), intent(in) :: temperature_guess
    real(dp), intent(out) :: jacobian(:,:),temperature
    logical, intent(out) :: ok
    integer, optional, intent(in) :: dependent_species
    logical, optional, intent(in) :: permit_signed_trials
    integer :: dependent
    real(dp) :: density,energy,rates(size(species)),recovered
    logical :: good,signed_trial
    signed_trial=.false.
    if(present(permit_signed_trials)) signed_trial=permit_signed_trials
    jacobian=0;temperature=0;ok=.false.
    if(size(species)<2) return
    if(any(shape(jacobian)/=[size(species)-1,size(species)-1])) return
    if(present(dependent_species)) then
      if(dependent_species<1.or.dependent_species>size(species)) return
    endif
    call forced_reactor_rhs(species,reactions,initial_energy_density,elapsed,species_density, &
      external_species,external_energy,temperature_guess,rates,recovered,good,permit_signed_trials=signed_trial)
    if(.not.good) return
    if(size(reactions)==0) then
      ! Constant external forcing has an exactly zero state Jacobian.
      temperature=recovered;ok=.true.;return
    endif
    density=sum(species_density)
    energy=(initial_energy_density+elapsed*external_energy)/density
    if(signed_trial.and.any(species_density<0)) then
      dependent=size(species)
      if(present(dependent_species)) dependent=dependent_species
      call signed_trial_jacobian(species,reactions,density,energy,species_density/density, &
        recovered,jacobian,temperature,ok,dependent)
      return
    endif
    call reactor_reduced_jacobian(species,reactions,density,energy,species_density/density, &
      recovered,jacobian,temperature,ok,dependent_species)
  end subroutine forced_reactor_reduced_jacobian

  subroutine advance_forced_reactor_implicit_adaptive(species,reactions,initial_energy_density,elapsed, &
      external_species,external_energy,requested_time_step,relative_tolerance,absolute_tolerance, &
      species_density,temperature,accepted_time_step,next_time_step,ok, &
      permit_signed_trials)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: initial_energy_density,elapsed,external_species(:),external_energy
    real(dp), intent(in) :: requested_time_step,relative_tolerance,absolute_tolerance
    real(dp), intent(inout) :: species_density(:),temperature
    real(dp), intent(out) :: accepted_time_step,next_time_step
    logical, intent(out) :: ok
    logical, optional, intent(in) :: permit_signed_trials
    real(dp) :: full(size(species)),half(size(species)),middle(size(species)),candidate(size(species))
    real(dp) :: rates(size(species)),tf,th,tm,tc,h,err,factor
    logical :: good,extrapolated_ok,signed_trial
    integer :: attempt

    signed_trial=.false.
    if(present(permit_signed_trials)) signed_trial=permit_signed_trials
    accepted_time_step=0;next_time_step=0;ok=.false.
    if(.not.all(ieee_is_finite([requested_time_step,relative_tolerance,absolute_tolerance]))) return
    if(min(requested_time_step,relative_tolerance,absolute_tolerance)<=0) return
    call forced_reactor_rhs(species,reactions,initial_energy_density,elapsed,species_density, &
      external_species,external_energy,temperature,rates,tc,good,permit_signed_trials=signed_trial)
    if(.not.good) return
    h=requested_time_step
    do attempt=1,reactor_max_step_attempts
      if(elapsed+h<=elapsed) return
      call trial(species_density,temperature,elapsed,h,full,tf,good)
      if(good) call trial(species_density,temperature,elapsed,.5_dp*h,middle,tm,good)
      if(good) call trial(middle,tm,elapsed+.5_dp*h,.5_dp*h,half,th,good)
      if(.not.good) then
        h=.25_dp*h
        cycle
      end if
      err=normalized_error(full,tf,half,th,relative_tolerance,absolute_tolerance)
      if(.not.ieee_is_finite(err)) then
        h=.25_dp*h
        cycle
      end if
      if(err>1) then
        h=h*max(.1_dp,.9_dp*err**(-.5_dp))
        cycle
      end if
      candidate=2*half-full
      ! Richardson extrapolation can reverse the sign of a newly produced
      ! high-order trace even when both backward-Euler endpoints are physical.
      ! Keep its physical admission check; signed Newton coordinates remain
      ! available in the half-step result when affine forcing requires them.
      call forced_reactor_rhs(species,reactions,initial_energy_density,elapsed+h,candidate, &
        external_species,external_energy,th,rates,tc,extrapolated_ok)
      if(.not.extrapolated_ok) then
        candidate=half
        tc=th
      end if
      factor=2
      if(err>1e-12_dp) factor=min(2._dp,max(.2_dp,.9_dp*err**(-.5_dp)))
      species_density=candidate;temperature=tc
      accepted_time_step=h;next_time_step=h*factor;ok=.true.
      return
    end do
  contains
    subroutine trial(start,start_temperature,time,step,finish,finish_temperature,valid)
      real(dp), intent(in) :: start(:),start_temperature,time,step
      real(dp), intent(out) :: finish(:),finish_temperature
      logical, intent(out) :: valid
      real(dp) :: density,energy,reference(size(species)),y(size(species)),end_y(size(species))
      real(dp) :: unused(size(species))
      integer :: iterations
      valid=.false.;finish=0;finish_temperature=0
      density=sum(start)+step*sum(external_species)
      if(.not.ieee_is_finite(density).or.density<=0) return
      energy=(initial_energy_density+(time+step)*external_energy)/density
      if(.not.ieee_is_finite(energy)) return
      reference=(start+step*external_species)/density
      if(size(reactions)==0) then
        finish=start+step*external_species
        call forced_reactor_rhs(species,reactions,initial_energy_density,time+step,finish, &
          external_species,external_energy,start_temperature,unused,finish_temperature,valid,permit_signed_trials=signed_trial)
        return
      end if
      y=start/sum(start)
      ! Carry admissible affine forcing into the Newton starting guess. An
      ! unchanged guess can satisfy atol while erasing a newly advected tracer.
      if(all(reference>=0)) y=reference
      call backward_euler_trial(species,reactions,density,energy,y,start_temperature,step, &
        relative_tolerance,absolute_tolerance/density,end_y,finish_temperature,iterations,valid,reference, &
        preserve_signed_trials=signed_trial)
      if(valid) finish=density*end_y
    end subroutine trial
  end subroutine advance_forced_reactor_implicit_adaptive

  ! Transactional full interval: publish only after every substep succeeds.
  subroutine advance_forced_reactor_interval(species,reactions,initial_energy_density,external_species, &
      external_energy,duration,relative_tolerance,absolute_tolerance,species_density,temperature, &
      completed_steps,ok,initial_time_step,maximum_steps)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: initial_energy_density,external_species(:),external_energy,duration
    real(dp), intent(in) :: relative_tolerance,absolute_tolerance
    real(dp), intent(inout) :: species_density(:),temperature
    integer, intent(out) :: completed_steps
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: initial_time_step
    integer, intent(in), optional :: maximum_steps
    real(dp) :: candidate(size(species)),rates(size(species)),candidate_temperature,recovered_temperature
    real(dp) :: elapsed,step,accepted,next_step,remaining
    integer :: limit
    logical :: valid,inactive_forcing

    ok=.false.;completed_steps=0;limit=100000
    if(present(maximum_steps)) limit=maximum_steps
    if(limit<1) return
    if(.not.all(ieee_is_finite([duration,relative_tolerance,absolute_tolerance]))) return
    if(duration<0.or.min(relative_tolerance,absolute_tolerance)<=0) return
    step=duration
    if(present(initial_time_step)) then
      if(.not.ieee_is_finite(initial_time_step)) return
      if(initial_time_step<=0) return
      step=min(step,initial_time_step)
    end if
    inactive_forcing=.false.
    if(.not.present(initial_time_step)) &
      inactive_forcing=chemically_inactive_forcing(reactions,species_density,external_species,size(species))
    if(inactive_forcing) then
      call forced_reactor_rhs(species,reactions(:0),initial_energy_density,0._dp,species_density, &
        external_species,external_energy,temperature,rates,recovered_temperature,valid)
    else
    call forced_reactor_rhs(species,reactions,initial_energy_density,0._dp,species_density, &
      external_species,external_energy,temperature,rates,recovered_temperature,valid)
    endif
    if(.not.valid) return
    if(duration==0) then
      ok=.true.
      return
    end if
    if((size(reactions)==0.or.inactive_forcing).and..not.present(initial_time_step)) then
      ! With no chemical reactions and constant forcing, rhoY is linear in
      ! time. Solve the full interval directly so adaptive half-step rounding
      ! cannot be reported by SDC as a fictitious chemical source. Explicit
      ! substep requests retain the existing step-count/limit contract.
      candidate=species_density+duration*external_species
      call forced_reactor_rhs(species,reactions(:0),initial_energy_density,duration,candidate, &
        external_species,external_energy,recovered_temperature,rates,candidate_temperature,valid)
      if(.not.valid) return
      species_density=candidate;temperature=candidate_temperature;completed_steps=1;ok=.true.
      return
    endif
    ! A negative conserved molar pool cannot reach any nonnegative endpoint.
    ! Do not consume the adaptive step limit approaching its signed-trial bound
    ! before returning the same transactional refusal. This is a sufficient
    ! certificate only; affine-negative species that chemistry can replenish
    ! retain the ordinary signed Newton integration.
    if(negative_conserved_molar_pool(species,reactions,species_density,external_species,duration)) return
    candidate=species_density;candidate_temperature=recovered_temperature;elapsed=0
    do while(elapsed<duration)
      if(completed_steps>=limit) return
      remaining=duration-elapsed
      step=min(step,remaining)
      call advance_forced_reactor_implicit_adaptive(species,reactions,initial_energy_density,elapsed, &
        external_species,external_energy,step,relative_tolerance,absolute_tolerance,candidate, &
        candidate_temperature,accepted,next_step,valid,permit_signed_trials=.true.)
      if(.not.valid) return
      if(accepted<=0.or.accepted>remaining.or.elapsed+accepted<=elapsed) return
      completed_steps=completed_steps+1
      if(accepted==remaining) then
        elapsed=duration
      else
        elapsed=elapsed+accepted
      end if
      step=next_step
    end do
    ! Signed internal trial coordinates never bypass final physical admission.
    call forced_reactor_rhs(species,reactions,initial_energy_density,duration,candidate, &
      external_species,external_energy,candidate_temperature,rates,recovered_temperature,valid)
    if(.not.valid) return
    species_density=candidate;temperature=recovered_temperature;ok=.true.
  end subroutine advance_forced_reactor_interval

  logical function chemically_inactive_forcing(reactions,densities,forcing,ns) result(inactive)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: densities(:),forcing(:)
    integer, intent(in) :: ns
    integer :: r,i
    ! Sufficient structural certificate, never a small-rate approximation:
    ! every species in either reaction direction stays identically absent.
    ! Third-body-only species may remain present, but cannot create reactants.
    inactive=.false.
    if(ns<1.or.size(densities)/=ns.or.size(forcing)/=ns) return
    if(.not.all(ieee_is_finite(densities)).or..not.all(ieee_is_finite(forcing))) return
    do r=1,size(reactions)
      if(.not.valid_elementary_reaction(reactions(r),ns)) return
      do i=1,ns
        if(reactions(r)%reactant_stoich(i)>0.or.reactions(r)%product_stoich(i)>0) then
          if(densities(i)/=0.or.forcing(i)/=0) return
        endif
      enddo
    enddo
    inactive=.true.
  end function chemically_inactive_forcing

  logical function negative_conserved_molar_pool(species,reactions,densities,forcing,duration) result(impossible)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: densities(:),forcing(:),duration
    integer :: labels(size(species)),ns,r,i,pool,leader,previous
    real(dp) :: net(size(species)),increment(size(species)),endpoint(size(species)),bound
    real(dp) :: old_moles,forced_moles,total,roundoff_scale,term,compensation,shift,updated,absolute_term
    logical :: conserved(size(species))
    impossible=.false.;ns=size(species)
    ! The caller has validated the physical input. An overflow or unsupported
    ! stoichiometry leaves admission to the existing integration path.
    if(size(reactions)==0.or.ns==0) return
    do i=1,ns
      if(abs(forcing(i))>huge(1._dp)/max(1._dp,duration)) return
      increment(i)=duration*forcing(i)
      if(densities(i)>huge(1._dp)-max(increment(i),0._dp)) return
      if(densities(i)< -huge(1._dp)-min(increment(i),0._dp)) return
      endpoint(i)=densities(i)+increment(i)
    enddo
    if(all(endpoint>=0)) return
    labels=[(i,i=1,ns)]
    bound=1._dp/(epsilon(1._dp)*real(ns,dp))
    do r=1,size(reactions)
      net=reactions(r)%product_stoich-reactions(r)%reactant_stoich
      ! Integral net coefficients with a bounded sum admit an exact binary64
      ! zero test. Fractional or very large coefficients provide no certificate.
      if(any(abs(net)>bound)) return
      if(any(net/=anint(net))) return
      if(all(net==0)) cycle
      leader=minval(labels,mask=net/=0)
      do i=1,ns
        if(net(i)==0) cycle
        previous=labels(i)
        where(labels==previous) labels=leader
      enddo
    enddo
    conserved=.true.
    do r=1,size(reactions)
      net=reactions(r)%product_stoich-reactions(r)%reactant_stoich
      if(all(net==0)) cycle
      pool=minval(labels,mask=net/=0)
      if(sum(net)/=0) conserved(pool)=.false.
    enddo
    do pool=1,ns
      if(.not.conserved(pool).or..not.any(labels==pool)) cycle
      total=0;roundoff_scale=0;compensation=0
      do i=1,ns
        if(labels(i)/=pool) cycle
        if(species(i)%molecular_weight<1) then
          if(abs(densities(i))>huge(1._dp)*species(i)%molecular_weight) return
          if(abs(increment(i))>huge(1._dp)*species(i)%molecular_weight) return
        endif
        old_moles=densities(i)/species(i)%molecular_weight
        forced_moles=increment(i)/species(i)%molecular_weight
        if(max(abs(old_moles),abs(forced_moles))>huge(1._dp)/(4*real(ns,dp))) return
        if(abs(old_moles)>huge(1._dp)-abs(forced_moles)) return
        absolute_term=abs(old_moles)+abs(forced_moles)
        if(roundoff_scale>huge(1._dp)-absolute_term) return
        roundoff_scale=roundoff_scale+absolute_term
        term=old_moles+forced_moles
        if(total>huge(1._dp)-max(term,0._dp)) return
        if(total< -huge(1._dp)-min(term,0._dp)) return
        shift=term-compensation;updated=total+shift
        if(.not.ieee_is_finite(updated)) return
        compensation=(updated-total)-shift;total=updated
      enddo
      ! Cancellation near a zero inventory is inconclusive, so the existing
      ! reactor/final physical gate decides it. No source or state is clipped.
      if(total< -64*epsilon(1._dp)*roundoff_scale) then
        impossible=.true.;return
      endif
    enddo
  end function negative_conserved_molar_pool

  ! The same existing Newton composition bound applies to signed trials.
  ! Preserve their coordinates and affine balance instead of projecting them.
  subroutine admit_trial_composition(y,ok,signed_trial)
    real(dp), intent(inout) :: y(:)
    logical, intent(out) :: ok
    logical, intent(in) :: signed_trial
    if(.not.signed_trial) then
      call enforce_roundoff_composition(y,ok)
      return
    endif
    ok=.false.
    if(.not.all(ieee_is_finite(y))) return
    if(any(y< -5e-13_dp).or.any(y>1+5e-13_dp)) return
    ok=abs(sum(y)-1)<=2.5e-12_dp
  end subroutine admit_trial_composition

  subroutine recover_trial_temperature(species,y,energy,guess,temperature,ok,signed_trial)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: y(:),energy,guess
    real(dp), intent(out) :: temperature
    logical, intent(out) :: ok
    logical, intent(in) :: signed_trial
    real(dp) :: lower,upper,elo,ehi,current,cv,residual,candidate,t,tolerance,cv_lower,cv_upper
    integer :: iteration
    temperature=0;ok=.false.
    if(.not.signed_trial.or..not.any(y<0)) then
      call temperature_from_internal_energy(species,y,energy,guess,temperature,ok)
      return
    endif
    if(.not.valid_mixture_composition(species,y)) return
    if(.not.all(ieee_is_finite([energy,guess]))) return
    call nasa7_eos_bounds(species,lower,upper,ok)
    if(.not.ok) return
    call properties(lower,elo,cv,ok);if(.not.ok) return
    cv_lower=cv
    call properties(upper,ehi,cv,ok);if(.not.ok) return
    cv_upper=cv
    if(species(1)%pinned_eos_extrapolation) then
      call pinned_signed_inverse(candidate,ok)
      if(ok) temperature=candidate
      return
    endif
    tolerance=1e-11_dp*max(1._dp,abs(energy))
    if(energy<elo-tolerance.or.energy>ehi+tolerance) then
      ok=.false.;return
    endif
    t=min(upper,max(lower,guess));ok=.false.
    do iteration=1,100
      call properties(t,current,cv,ok);if(.not.ok) return
      residual=current-energy
      if(abs(residual)<=tolerance) then
        temperature=t;ok=.true.;return
      endif
      if(residual>0) then
        upper=t;ehi=current
      else
        lower=t;elo=current
      endif
      if(nearest(lower,1._dp)>=upper.and. &
        .not.any(lower<species%temperature_mid.and.upper>=species%temperature_mid)) then
        temperature=lower
        if(abs(ehi-energy)<abs(elo-energy)) temperature=upper
        ok=.true.;return
      endif
      candidate=t-residual/cv
      if(.not.ieee_is_finite(candidate)) then
        ok=.false.;return
      endif
      if(candidate<=lower.or.candidate>=upper) candidate=.5_dp*lower+.5_dp*upper
      t=candidate
    enddo
    ok=.false.
  contains
    subroutine properties(t,e,cv,good)
      real(dp), intent(in) :: t
      real(dp), intent(out) :: e,cv
      logical, intent(out) :: good
      real(dp) :: cp_s,cv_s,h_s,e_s,s_s,term
      integer :: k
      e=0;cv=0;good=.false.
      do k=1,size(species)
        call nasa7_mass_properties(species(k),t,cp_s,cv_s,h_s,e_s,s_s,good)
        if(.not.good) return
        if(abs(e_s)>huge(1._dp)/4.or.abs(cv_s)>huge(1._dp)/4) then
          good=.false.;return
        endif
        term=y(k)*e_s
        if(abs(e)>huge(1._dp)/2-abs(term)) then
          good=.false.;return
        endif
        e=e+term;term=y(k)*cv_s
        if(abs(cv)>huge(1._dp)/2-abs(term)) then
          good=.false.;return
        endif
        cv=cv+term
      enddo
      good=all(ieee_is_finite([e,cv])).and.cv>0
    end subroutine properties
    subroutine pinned_signed_inverse(result,good)
      real(dp), intent(out) :: result
      logical, intent(out) :: good
      real(dp) :: value,heat_capacity,delta,delta_energy,denominator,ratio,working,anchor,temperature_tolerance
      integer :: iteration,maximum_iterations
      result=0;good=.false.
      maximum_iterations=200;temperature_tolerance=1e-6_dp
      if(species(1)%pinned_eos_convergence) then
        maximum_iterations=5000;temperature_tolerance=1e-12_dp
      endif
      ! Use the same boundary-Cv linear extension and native200-update
      ! inverse as the physical Fuego EOS. Signed trial coordinates retain
      ! their own caloric sums rather than clipping a temporary population.
      if(energy<elo.or.energy>ehi) then
        anchor=lower;value=elo;heat_capacity=cv_lower
        if(energy>ehi) then
          anchor=upper;value=ehi;heat_capacity=cv_upper
        endif
        call reactor_energy_difference(energy,value,delta_energy,good)
        if(.not.good) return
        if(heat_capacity<1) then
          if(abs(delta_energy)>huge(1._dp)*heat_capacity) then
            good=.false.;return
          endif
        endif
        delta=delta_energy/heat_capacity
        call reactor_energy_difference(anchor,-delta,working,good)
        if(.not.good) return
        if(working<=0) then
          good=.false.;return
        endif
        call properties(working,value,heat_capacity,good)
        if(good) result=working
        return
      endif
      working=guess
      if(working<lower.or.working>upper) then
        call reactor_energy_difference(ehi,elo,denominator,good)
        if(.not.good) return
        if(denominator<=0) then
          good=.false.;return
        endif
        if(denominator<1) then
          if(upper-lower>huge(1._dp)*denominator) then
            good=.false.;return
          endif
        endif
        ratio=(upper-lower)/denominator
        call reactor_energy_difference(energy,elo,delta_energy,good)
        if(.not.good) return
        if(ratio>1) then
          if(abs(delta_energy)>huge(1._dp)/ratio) then
            good=.false.;return
          endif
        endif
        delta=ratio*delta_energy
        call reactor_energy_difference(lower,-delta,working,good)
        if(.not.good) return
      endif
      do iteration=1,maximum_iterations
        call properties(working,value,heat_capacity,good)
        if(.not.good) return
        call reactor_energy_difference(energy,value,delta_energy,good)
        if(.not.good) return
        if(heat_capacity<1) then
          if(abs(delta_energy)>huge(1._dp)*heat_capacity) then
            good=.false.;return
          endif
        endif
        delta=delta_energy/heat_capacity
        if(delta>100._dp) then
          delta=100._dp
        else if(delta< -100._dp) then
          delta=-100._dp
        else if(abs(delta)<temperature_tolerance) then
          result=working;good=.true.;return
        endif
        call reactor_energy_difference(working,-delta,anchor,good)
        if(.not.good) return
        working=anchor
        if(working<=0) then
          good=.false.;return
        endif
      enddo
      good=.false.
    end subroutine pinned_signed_inverse
  end subroutine recover_trial_temperature

  subroutine signed_trial_rhs(species,reactions,density,energy,y,guess,rhs,temperature,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: density,energy,y(:),guess
    real(dp), intent(out) :: rhs(:),temperature
    logical, intent(out) :: ok
    real(dp) :: concentrations(size(species)),production(size(species)),value
    integer :: k
    rhs=0;temperature=0;ok=.false.
    if(size(rhs)/=size(species).or..not.ieee_is_finite(density)) return
    if(density<=0.or..not.valid_mixture_composition(species,y)) return
    call recover_trial_temperature(species,y,energy,guess,temperature,ok,.true.)
    if(.not.ok) return
    do k=1,size(species)
      if(density>huge(1._dp)/max(1._dp,abs(y(k)))) then
        ok=.false.;return
      endif
      value=density*y(k)
      if(species(k)%molecular_weight<1) then
        if(abs(value)>huge(1._dp)*species(k)%molecular_weight) then
          ok=.false.;return
        endif
      endif
      concentrations(k)=value/species(k)%molecular_weight
    enddo
    call elementary_production_rates_from_concentrations(species,reactions,temperature,concentrations, &
      production,ok,signed_concentrations=.true.)
    if(.not.ok) return
    do k=1,size(species)
      if(abs(production(k))>huge(1._dp)/max(1._dp,species(k)%molecular_weight)) then
        ok=.false.;return
      endif
      value=production(k)*species(k)%molecular_weight
      if(density<1) then
        if(abs(value)>huge(1._dp)*density) then
          ok=.false.;return
        endif
      endif
      rhs(k)=value/density
    enddo
    ok=all(ieee_is_finite(rhs))
  end subroutine signed_trial_rhs

  subroutine signed_trial_jacobian(species,reactions,density,energy,y,guess,jacobian,temperature,ok,dependent)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: density,energy,y(:),guess
    real(dp), intent(out) :: jacobian(:,:),temperature
    logical, intent(out) :: ok
    integer, intent(in) :: dependent
    real(dp) :: concentrations(size(species)),molar_jacobian(size(species),size(species))
    real(dp) :: low_rates(size(species)),high_rates(size(species)),thermal_derivative(size(species)),ei(size(species))
    real(dp) :: lower,upper,h,cv,cp_s,cv_s,h_s,e_s,s_s,dt_dy,mw_i,reference_temperature
    integer :: i,j,row,column,coordinates(size(species)-1)
    jacobian=0;temperature=0;ok=.false.
    coordinates=[(reduced_species_index(i,dependent),i=1,size(species)-1)]
    call recover_trial_temperature(species,y,energy,guess,temperature,ok,.true.);if(.not.ok) return
    concentrations=density*y/species%molecular_weight
    call elementary_production_jacobian(species,reactions,temperature,concentrations,molar_jacobian,ok, &
      signed_concentrations=.true.)
    if(.not.ok) return
    h=max(1e-4_dp,1e-6_dp*temperature)
    call reactor_temperature_difference_interval(species,temperature,h,lower,upper,reference_temperature,ok)
    if(.not.ok) return
    if(upper<=lower) then
      ok=.false.;return
    endif
    call elementary_production_rates_from_concentrations(species,reactions,lower,concentrations,low_rates,ok, &
      signed_concentrations=.true.)
    if(.not.ok) return
    call elementary_production_rates_from_concentrations(species,reactions,upper,concentrations,high_rates,ok, &
      signed_concentrations=.true.)
    if(.not.ok) return
    thermal_derivative=(high_rates-low_rates)/(upper-lower)*species%molecular_weight/density
    cv=0
    do i=1,size(species)
      call nasa7_mass_properties(species(i),reference_temperature,cp_s,cv_s,h_s,e_s,s_s,ok)
      if(.not.ok) return
      cv=cv+y(i)*cv_s
      call reactor_linearized_species_energy(e_s,cv_s,temperature-reference_temperature,ei(i),ok)
      if(.not.ok) return
    enddo
    if(.not.ieee_is_finite(cv).or.cv<=0) then
      ok=.false.;return
    endif
    ! Exact signed concentration derivatives avoid a finite species perturbation
    ! overwhelming a trace product. Temperature differentiation keeps the same
    ! step as the preceding physical Jacobian, including formation energies.
    do j=1,size(coordinates)
      column=coordinates(j);dt_dy=-(ei(column)-ei(dependent))/cv
      do i=1,size(coordinates)
        row=coordinates(i);mw_i=species(row)%molecular_weight
        jacobian(i,j)=mw_i/species(column)%molecular_weight*molar_jacobian(row,column)- &
          mw_i/species(dependent)%molecular_weight*molar_jacobian(row,dependent)+thermal_derivative(row)*dt_dy
      enddo
    enddo
    ok=all(ieee_is_finite(jacobian))
  end subroutine signed_trial_jacobian

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
      temperature_guess, derivative, temperature, ok, &
      permit_signed_trials)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: density, target_internal_energy
    real(dp), intent(in) :: mass_fractions(:), temperature_guess
    real(dp), intent(out) :: derivative(:), temperature
    logical, intent(out) :: ok
    logical, optional, intent(in) :: permit_signed_trials

    derivative=0;temperature=0;ok=.false.
    if(size(mass_fractions)/=size(species).or.size(derivative)/=size(species)) return
    if(.not.all(ieee_is_finite(mass_fractions))) return
    if(present(permit_signed_trials)) then
      if(permit_signed_trials.and.any(mass_fractions<0)) then
        call signed_trial_rhs(species,reactions,density,target_internal_energy,mass_fractions, &
          temperature_guess,derivative,temperature,ok)
        return
      endif
    endif
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



  subroutine reactor_reduced_jacobian( &
      species, reactions, density, target_internal_energy, mass_fractions, &
      temperature_guess, jacobian, temperature, ok, dependent_species)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: density, target_internal_energy
    real(dp), intent(in) :: mass_fractions(:), temperature_guess
    real(dp), intent(out) :: jacobian(:, :), temperature
    logical, intent(out) :: ok
    integer, intent(in), optional :: dependent_species

    real(dp), allocatable :: fixed_temperature_jacobian(:, :)
    real(dp), allocatable :: rhs_lower(:), rhs_upper(:)
    real(dp), allocatable :: d_rhs_d_temperature(:)
    real(dp), allocatable :: species_internal_energy(:)
    real(dp) :: molecular_weight, gas_constant, cp, cv, gamma
    real(dp) :: enthalpy, mixture_internal_energy, entropy
    real(dp) :: species_cp, species_cv, species_h, species_u, species_s
    real(dp) :: lower_temperature, upper_temperature, temperature_delta
    real(dp) :: reference_temperature,d_temperature_d_variable
    logical :: property_ok
    integer :: dependent, i, j, nspecies, reduced_size
    integer :: row_species, column_species

    jacobian = 0.0_dp
    temperature = 0.0_dp
    nspecies = size(species)
    reduced_size = nspecies - 1
    dependent = nspecies
    if (present(dependent_species)) dependent = dependent_species
    ok = nspecies >= 2 .and. size(mass_fractions) == nspecies .and. &
      size(jacobian, 1) == reduced_size .and. &
      size(jacobian, 2) == reduced_size .and. density > 0.0_dp .and. &
      dependent >= 1 .and. dependent <= nspecies
    if (.not. ok) return
    ok = valid_mixture_composition(species, mass_fractions)
    if (.not. ok) return

    call temperature_from_internal_energy( &
      species, mass_fractions, target_internal_energy, temperature_guess, &
      temperature, ok)
    if (.not. ok) return

    allocate(fixed_temperature_jacobian(nspecies, nspecies))
    allocate(rhs_lower(nspecies), rhs_upper(nspecies))
    allocate(d_rhs_d_temperature(nspecies), species_internal_energy(nspecies))
    call elementary_mass_fraction_jacobian( &
      species, reactions, temperature, density, mass_fractions, &
      fixed_temperature_jacobian, ok)
    if (.not. ok) return

    temperature_delta = max(1.0e-4_dp, 1.0e-6_dp * temperature)
    call reactor_temperature_difference_interval(species,temperature,temperature_delta, &
      lower_temperature,upper_temperature,reference_temperature,ok)
    if(.not.ok) return
    if (upper_temperature <= lower_temperature) then
      ok = .false.
      return
    end if
    call elementary_mass_fraction_rhs( &
      species, reactions, lower_temperature, density, mass_fractions, &
      rhs_lower, ok)
    if (.not. ok) return
    call elementary_mass_fraction_rhs( &
      species, reactions, upper_temperature, density, mass_fractions, &
      rhs_upper, ok)
    if (.not. ok) return
    d_rhs_d_temperature = (rhs_upper - rhs_lower) / &
      (upper_temperature - lower_temperature)

    call mixture_mass_properties( &
      species, mass_fractions, temperature, molecular_weight, gas_constant, &
      cp, cv, gamma, enthalpy, mixture_internal_energy, entropy, ok)
    if (.not. ok .or. cv <= 0.0_dp) then
      ok = .false.
      return
    end if
    if(reference_temperature/=temperature) then
      call mixture_mass_properties(species,mass_fractions,reference_temperature,molecular_weight,gas_constant, &
        cp,cv,gamma,enthalpy,mixture_internal_energy,entropy,ok)
      if(.not.ok) return
    endif
    do i = 1, nspecies
      call nasa7_mass_properties( &
        species(i), reference_temperature, species_cp, species_cv, species_h, &
        species_u, species_s, property_ok)
      if (.not. property_ok) then
        ok = .false.
        return
      end if
      call reactor_linearized_species_energy(species_u,species_cv,temperature-reference_temperature, &
        species_internal_energy(i),property_ok)
      if(.not.property_ok) then
        ok=.false.;return
      endif
    end do

    do j = 1, reduced_size
      column_species = reduced_species_index(j, dependent)
      d_temperature_d_variable = -( &
        species_internal_energy(column_species) - &
        species_internal_energy(dependent)) / cv
      do i = 1, reduced_size
        row_species = reduced_species_index(i, dependent)
        jacobian(i, j) = &
          fixed_temperature_jacobian(row_species, column_species) - &
          fixed_temperature_jacobian(row_species, dependent) + &
          d_rhs_d_temperature(row_species) * d_temperature_d_variable
      end do
    end do
    ok = all(ieee_is_finite(jacobian))
  end subroutine reactor_reduced_jacobian

  pure integer function reduced_species_index( &
      reduced_index, dependent_species) result(species_index)
    integer, intent(in) :: reduced_index, dependent_species

    species_index = reduced_index
    if (reduced_index >= dependent_species) species_index = reduced_index + 1
  end function reduced_species_index

  subroutine backward_euler_trial( &
      species, reactions, density, target_internal_energy, initial_state, &
      initial_temperature, time_step, relative_tolerance, absolute_tolerance, &
      updated_state, updated_temperature, newton_iterations, ok, residual_reference, &
      preserve_signed_trials)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: density, target_internal_energy
    real(dp), intent(in) :: initial_state(:), initial_temperature, time_step
    real(dp), intent(in) :: relative_tolerance, absolute_tolerance
    real(dp), intent(out) :: updated_state(:), updated_temperature
    integer, intent(out) :: newton_iterations
    logical, intent(out) :: ok

    ! Optional externally forced reference need not lie in the simplex;
    ! initial_state remains a physical Newton starting guess.
    real(dp), intent(in), optional :: residual_reference(:)
    logical, optional, intent(in) :: preserve_signed_trials

    integer, parameter :: maximum_newton_iterations = 20
    integer, parameter :: maximum_line_search_iterations = 24
    real(dp), allocatable :: reduced_initial(:), reduced_state(:)
    real(dp), allocatable :: trial_reduced(:), full_state(:), trial_full(:)
    real(dp), allocatable :: rhs(:), trial_rhs(:), residual(:), trial_residual(:)
    real(dp), allocatable :: jacobian(:, :), matrix(:, :), correction(:)
    real(dp) :: temperature, evaluated_temperature, trial_temperature
    real(dp) :: trial_evaluated_temperature, residual_norm
    real(dp) :: trial_norm, line_factor
    logical :: solve_ok, trial_ok, forced_balance,signed_trial
    integer :: i, iteration, line_iteration, nspecies, reduced_size,dependent
    integer, allocatable :: coordinates(:)

    signed_trial=.false.
    if(present(preserve_signed_trials)) signed_trial=preserve_signed_trials
    signed_trial=signed_trial.and.present(residual_reference)
    updated_state = 0.0_dp
    updated_temperature = 0.0_dp
    newton_iterations = 0
    nspecies = size(species)
    reduced_size = nspecies - 1
    ok = nspecies >= 2 .and. size(initial_state) == nspecies .and. &
      size(updated_state) == nspecies .and. density > 0.0_dp .and. &
      time_step > 0.0_dp .and. relative_tolerance > 0.0_dp .and. &
      absolute_tolerance > 0.0_dp
    if (.not. ok) return
    ok = valid_mixture_composition(species, initial_state)
    if (.not. ok) return

    allocate(reduced_initial(reduced_size), reduced_state(reduced_size))
    allocate(trial_reduced(reduced_size), full_state(nspecies))
    allocate(trial_full(nspecies), rhs(nspecies), trial_rhs(nspecies))
    allocate(residual(reduced_size), trial_residual(reduced_size))
    allocate(jacobian(reduced_size, reduced_size))
    allocate(matrix(reduced_size, reduced_size), correction(reduced_size))
    ! The largest forced component carries the mass closure. Eliminating a
    ! trace product as 1-sum(other species) would erase it through cancellation.
    ! Ordinary unforced trials retain their existing coordinate convention.
    dependent=nspecies
    if(present(residual_reference)) dependent=maxloc(initial_state,dim=1)
    allocate(coordinates(reduced_size))
    coordinates=[(reduced_species_index(i,dependent),i=1,reduced_size)]
    reduced_initial = initial_state(coordinates)
    if(present(residual_reference)) then
      ok=size(residual_reference)==nspecies
      if(.not.ok) return
      ok=all(ieee_is_finite(residual_reference))
      if(.not.ok) return
      reduced_initial=residual_reference(coordinates)
    end if
    reduced_state = initial_state(coordinates)
    temperature = initial_temperature

    do iteration = 1, maximum_newton_iterations
      call reconstruct_full_state(reduced_state, full_state, ok,dependent,signed_trial)
      if (.not. ok) return
      if(present(residual_reference)) reduced_state=full_state(coordinates)
      call reactor_rhs( &
        species, reactions, density, target_internal_energy, full_state, &
        temperature, rhs, evaluated_temperature, ok,permit_signed_trials=signed_trial)
      if (.not. ok) return
      temperature = evaluated_temperature
      residual = reduced_state - reduced_initial - &
        time_step * rhs(coordinates)
      residual_norm = implicit_residual_norm( &
        residual, reduced_state, reduced_initial, &
        relative_tolerance, absolute_tolerance)
      if (.not. ieee_is_finite(residual_norm)) then
        ok = .false.
        return
      end if
      forced_balance=.true.
      if(present(residual_reference)) then
        ! SDC exports the chemical increment separately from affine transport.
        ! An atol-sized bulk residual can hide an entire newly produced trace
        ! species. Resolve its balance relative to its represented inventory
        ! and source, while allowing ordinary floating-point cancellation.
        forced_balance=all(abs(residual)<=max(relative_tolerance,64*epsilon(1._dp))* &
          max(abs(reduced_state),abs(reduced_initial),abs(time_step*rhs(coordinates)),tiny(1._dp)))
      endif
      ! Forced updates feed an SDC reaction source and face predictor. Resolve
      ! a nonzero initial residual once, including newly produced trace species.
      if (residual_norm <= 5.0e-2_dp.and.forced_balance.and. &
          .not.(present(residual_reference).and.iteration==1.and.residual_norm>0)) then
        updated_state = full_state
        call admit_trial_composition(updated_state,ok,signed_trial)
        if (.not. ok) return
        call recover_trial_temperature(species,updated_state,target_internal_energy,temperature, &
          updated_temperature,ok,signed_trial)
        newton_iterations = iteration
        return
      end if

      if(signed_trial.and.any(full_state<0)) then
        call signed_trial_jacobian(species,reactions,density,target_internal_energy,full_state, &
          temperature,jacobian,evaluated_temperature,ok,dependent)
      else
      call reactor_reduced_jacobian( &
        species, reactions, density, target_internal_energy, full_state, &
        temperature, jacobian, evaluated_temperature, ok,dependent)
      endif
      if (.not. ok) return
      temperature = evaluated_temperature
      matrix = -time_step * jacobian
      do i = 1, reduced_size
        matrix(i, i) = matrix(i, i) + 1.0_dp
      end do
      call solve_dense_linear_system( &
        matrix, -residual, correction, solve_ok)
      if (.not. solve_ok) then
        ok = .false.
        return
      end if
      if(present(residual_reference).and.residual_norm<=5.0e-2_dp) then
        ! In a stiff equilibrium, rate cancellation can exceed the trace's
        ! inventory-relative residual. The solved Newton correction measures
        ! the state error after that stiffness is accounted for. An absent
        ! product still needs a representable update before this can converge.
        if(all(abs(correction)<=max(relative_tolerance,64*epsilon(1._dp))* &
          max(abs(reduced_state),abs(reduced_initial),tiny(1._dp)))) then
          updated_state=full_state
          call admit_trial_composition(updated_state,ok,signed_trial)
          if(.not.ok) return
          call recover_trial_temperature(species,updated_state,target_internal_energy, &
            temperature,updated_temperature,ok,signed_trial)
          newton_iterations=iteration
          return
        endif
      endif
      ! The tolerance is already satisfied and no component can represent
      ! this correction. Recheck convergence without a futile line search.
      if(residual_norm<=5.0e-2_dp.and.all(reduced_state+correction==reduced_state)) cycle

      line_factor = 1.0_dp
      trial_ok = .false.
      do line_iteration = 1, maximum_line_search_iterations
        trial_reduced = reduced_state + line_factor * correction
        call reconstruct_full_state(trial_reduced, trial_full, trial_ok,dependent,signed_trial)
        if(trial_ok.and.present(residual_reference)) trial_reduced=trial_full(coordinates)
        if (trial_ok) then
          trial_temperature = temperature
          call reactor_rhs( &
            species, reactions, density, target_internal_energy, trial_full, &
            trial_temperature, trial_rhs, trial_evaluated_temperature, &
            trial_ok,permit_signed_trials=signed_trial)
          if (trial_ok) trial_temperature = trial_evaluated_temperature
        end if
        if (trial_ok) then
          trial_residual = trial_reduced - reduced_initial - &
            time_step * trial_rhs(coordinates)
          trial_norm = implicit_residual_norm( &
            trial_residual, trial_reduced, reduced_initial, &
            relative_tolerance, absolute_tolerance)
          ! A forced tracer update may change represented trace components
          ! while its maximum residual stays on the same rounding plateau.
          ! Both endpoints must satisfy the existing convergence criterion.
          if (ieee_is_finite(trial_norm).and.(trial_norm<residual_norm.or. &
              (present(residual_reference).and.residual_norm<=5.0e-2_dp.and.trial_norm<=5.0e-2_dp))) then
            reduced_state = trial_reduced
            temperature = trial_temperature
            exit
          end if
        end if
        line_factor = 0.5_dp * line_factor
      end do
      if (line_iteration > maximum_line_search_iterations) then
        ok = .false.
        return
      end if
    end do
    ok = .false.
  end subroutine backward_euler_trial

  subroutine advance_constant_volume_implicit_adaptive( &
      species, reactions, density, target_internal_energy, &
      requested_time_step, relative_tolerance, absolute_tolerance, &
      mass_fractions, temperature, accepted_time_step, next_time_step, &
      newton_iterations, rejected_attempts, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: density, target_internal_energy
    real(dp), intent(in) :: requested_time_step
    real(dp), intent(in) :: relative_tolerance, absolute_tolerance
    real(dp), intent(inout) :: mass_fractions(:), temperature
    real(dp), intent(out) :: accepted_time_step, next_time_step
    integer, intent(out) :: newton_iterations, rejected_attempts
    logical, intent(out) :: ok

    real(dp), allocatable :: full_state(:), half_state(:), intermediate(:)
    real(dp), allocatable :: extrapolated_state(:)
    real(dp) :: full_temperature, half_temperature, intermediate_temperature
    real(dp) :: extrapolated_temperature
    real(dp) :: trial_time_step, error_norm, factor
    logical :: full_ok, half_ok, extrapolated_ok
    integer :: attempt, full_iterations, half_iterations_1, half_iterations_2

    accepted_time_step = 0.0_dp
    next_time_step = 0.0_dp
    newton_iterations = 0
    rejected_attempts = 0
    ok = requested_time_step > 0.0_dp .and. density > 0.0_dp .and. &
      relative_tolerance > 0.0_dp .and. absolute_tolerance > 0.0_dp .and. &
      size(mass_fractions) == size(species)
    if (.not. ok) return
    ok = valid_mixture_composition(species, mass_fractions)
    if (.not. ok) return

    allocate(full_state(size(species)), half_state(size(species)))
    allocate(intermediate(size(species)), extrapolated_state(size(species)))
    trial_time_step = requested_time_step

    do attempt = 1, reactor_max_step_attempts
      call backward_euler_trial( &
        species, reactions, density, target_internal_energy, mass_fractions, &
        temperature, trial_time_step, relative_tolerance, absolute_tolerance, &
        full_state, full_temperature, full_iterations, full_ok)
      if (full_ok) then
        call backward_euler_trial( &
          species, reactions, density, target_internal_energy, mass_fractions, &
          temperature, 0.5_dp * trial_time_step, relative_tolerance, &
          absolute_tolerance, intermediate, intermediate_temperature, &
          half_iterations_1, half_ok)
      else
        half_ok = .false.
        half_iterations_1 = 0
      end if
      if (half_ok) then
        call backward_euler_trial( &
          species, reactions, density, target_internal_energy, intermediate, &
          intermediate_temperature, 0.5_dp * trial_time_step, &
          relative_tolerance, absolute_tolerance, half_state, &
          half_temperature, half_iterations_2, half_ok)
      else
        half_iterations_2 = 0
      end if

      if (.not. full_ok .or. .not. half_ok) then
        rejected_attempts = rejected_attempts + 1
        trial_time_step = 0.25_dp * trial_time_step
        cycle
      end if
      error_norm = normalized_error( &
        full_state, full_temperature, half_state, half_temperature, &
        relative_tolerance, absolute_tolerance)
      if (.not. ieee_is_finite(error_norm)) then
        rejected_attempts = rejected_attempts + 1
        trial_time_step = 0.25_dp * trial_time_step
        cycle
      end if

      if (error_norm <= 1.0_dp) then
        ! Backward Euler is first order.  Two half steps have half the leading
        ! local truncation error of one full step, so Richardson extrapolation
        ! cancels that term and provides the accepted second-order state.
        ! Near a positivity boundary the extrapolated composition can leave the
        ! simplex; in that case retain the already verified two-half-step state.
        extrapolated_state = 2.0_dp * half_state - full_state
        call enforce_roundoff_composition(extrapolated_state, extrapolated_ok)
        if (extrapolated_ok) then
          call temperature_from_internal_energy( &
            species, extrapolated_state, target_internal_energy, &
            half_temperature, extrapolated_temperature, extrapolated_ok)
        end if

        if (extrapolated_ok) then
          mass_fractions = extrapolated_state
          temperature = extrapolated_temperature
        else
          call enforce_roundoff_composition(half_state, half_ok)
          if (.not. half_ok) then
            rejected_attempts = rejected_attempts + 1
            trial_time_step = 0.25_dp * trial_time_step
            cycle
          end if
          mass_fractions = half_state
          call temperature_from_internal_energy( &
            species, mass_fractions, target_internal_energy, half_temperature, &
            temperature, half_ok)
          if (.not. half_ok) then
            rejected_attempts = rejected_attempts + 1
            trial_time_step = 0.25_dp * trial_time_step
            cycle
          end if
        end if
        if (error_norm <= 1.0e-12_dp) then
          factor = 2.0_dp
        else
          factor = min(2.0_dp, max(0.2_dp, &
            0.9_dp * error_norm**(-0.5_dp)))
        end if
        accepted_time_step = trial_time_step
        next_time_step = trial_time_step * factor
        newton_iterations = full_iterations + half_iterations_1 + &
          half_iterations_2
        ok = .true.
        return
      end if

      rejected_attempts = rejected_attempts + 1
      factor = max(0.1_dp, 0.9_dp * error_norm**(-0.5_dp))
      trial_time_step = trial_time_step * factor
    end do
    ok = .false.
  end subroutine advance_constant_volume_implicit_adaptive

  subroutine reconstruct_full_state(reduced_state, full_state, ok,dependent_species,preserve_signed)
    real(dp), intent(in) :: reduced_state(:)
    real(dp), intent(out) :: full_state(:)
    logical, intent(out) :: ok
    real(dp), parameter :: tolerance = 5.0e-13_dp
    integer, optional, intent(in) :: dependent_species
    logical, optional, intent(in) :: preserve_signed
    integer :: dependent,i

    full_state = 0.0_dp
    ok = size(full_state) == size(reduced_state) + 1 .and. &
      all(ieee_is_finite(reduced_state))
    if (.not. ok) return
    dependent=size(full_state)
    if(present(dependent_species)) dependent=dependent_species
    if(dependent<1.or.dependent>size(full_state)) then
      ok=.false.;return
    endif
    do i=1,size(reduced_state)
      full_state(reduced_species_index(i,dependent))=reduced_state(i)
    enddo
    full_state(dependent) = 1.0_dp - sum(reduced_state)
    if (any(full_state < -tolerance)) then
      ok = .false.
      return
    end if
    if(present(preserve_signed)) then
      if(preserve_signed) then
        ok=abs(sum(full_state)-1)<=5*tolerance
        return
      endif
    endif
    full_state = max(0.0_dp, full_state)
    if (abs(sum(full_state) - 1.0_dp) > 5.0_dp * tolerance) then
      ok = .false.
      return
    end if
    full_state = full_state / sum(full_state)
    ok = .true.
  end subroutine reconstruct_full_state

  real(dp) function implicit_residual_norm( &
      residual, state, reference_state, relative_tolerance, &
      absolute_tolerance) result(norm)
    real(dp), intent(in) :: residual(:), state(:), reference_state(:)
    real(dp), intent(in) :: relative_tolerance, absolute_tolerance
    real(dp) :: scale
    integer :: i

    norm = 0.0_dp
    do i = 1, size(residual)
      scale = absolute_tolerance + relative_tolerance * &
        max(abs(state(i)), abs(reference_state(i)))
      norm = max(norm, abs(residual(i)) / scale)
    end do
  end function implicit_residual_norm

  subroutine solve_dense_linear_system(matrix,right_hand_side,solution,ok)
    real(dp), intent(in) :: matrix(:,:),right_hand_side(:)
    real(dp), intent(out) :: solution(:)
    logical, intent(out) :: ok
    ! Keep ordinary accepted Newton corrections byte-for-byte unchanged.
    ! A global pivot threshold can reject an independently scaled equation;
    ! retry the same system in equilibrated coordinates, with the same pivot
    ! tolerance and an independent componentwise backward-error check.
    call solve_unscaled_dense_linear_system(matrix,right_hand_side,solution,ok)
    if(.not.ok) call solve_equilibrated_dense_linear_system(matrix,right_hand_side,solution,ok)
  end subroutine solve_dense_linear_system

  subroutine solve_unscaled_dense_linear_system(matrix, right_hand_side, solution, ok)
    real(dp), intent(in) :: matrix(:, :), right_hand_side(:)
    real(dp), intent(out) :: solution(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: work(:, :), rhs(:), temporary_row(:)
    real(dp) :: pivot_value, factor, threshold, temporary_value
    integer :: i, j, k, pivot, n

    solution = 0.0_dp
    n = size(right_hand_side)
    ok = n >= 1 .and. size(matrix, 1) == n .and. &
      size(matrix, 2) == n .and. size(solution) == n
    if (.not. ok) return
    allocate(work(n, n), rhs(n), temporary_row(n))
    work = matrix
    rhs = right_hand_side
    threshold = 100.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, maxval(abs(work)))

    do k = 1, n - 1
      pivot = k - 1 + maxloc(abs(work(k:n, k)), dim=1)
      pivot_value = abs(work(pivot, k))
      if (pivot_value <= threshold) then
        ok = .false.
        return
      end if
      if (pivot /= k) then
        temporary_row = work(k, :)
        work(k, :) = work(pivot, :)
        work(pivot, :) = temporary_row
        temporary_value = rhs(k)
        rhs(k) = rhs(pivot)
        rhs(pivot) = temporary_value
      end if
      do i = k + 1, n
        factor = work(i, k) / work(k, k)
        work(i, k) = 0.0_dp
        do j = k + 1, n
          work(i, j) = work(i, j) - factor * work(k, j)
        end do
        rhs(i) = rhs(i) - factor * rhs(k)
      end do
    end do
    if (abs(work(n, n)) <= threshold) then
      ok = .false.
      return
    end if

    solution(n) = rhs(n) / work(n, n)
    do i = n - 1, 1, -1
      solution(i) = (rhs(i) - dot_product( &
        work(i, i + 1:n), solution(i + 1:n))) / work(i, i)
    end do
    ok = all(ieee_is_finite(solution))
  end subroutine solve_unscaled_dense_linear_system

  subroutine solve_equilibrated_dense_linear_system(matrix,right_hand_side,solution,ok)
    real(dp), intent(in) :: matrix(:,:),right_hand_side(:)
    real(dp), intent(out) :: solution(:)
    logical, intent(out) :: ok
    real(dp) :: scaled(size(right_hand_side),size(right_hand_side)),rhs(size(right_hand_side))
    real(dp) :: column_scale(size(right_hand_side)),answer(size(right_hand_side)),row_scale
    real(dp) :: terms(size(right_hand_side)),magnitude,error,normalizer,normalized_rhs
    integer :: i,j,n
    logical :: valid
    solution=0;ok=.false.;n=size(right_hand_side)
    if(n<1.or.size(matrix,1)/=n.or.size(matrix,2)/=n.or.size(solution)/=n) return
    if(any(.not.ieee_is_finite(matrix)).or.any(.not.ieee_is_finite(right_hand_side))) return
    do i=1,n
      row_scale=maxval(abs(matrix(i,:)))
      if(row_scale==0) return
      if(row_scale<1) then
        if(abs(right_hand_side(i))>huge(1._dp)*row_scale) return
      endif
      scaled(i,:)=matrix(i,:)/row_scale
      rhs(i)=right_hand_side(i)/row_scale
    enddo
    do j=1,n
      column_scale(j)=maxval(abs(scaled(:,j)))
      if(column_scale(j)==0) return
      scaled(:,j)=scaled(:,j)/column_scale(j)
    enddo
    call solve_unscaled_dense_linear_system(scaled,rhs,answer,valid)
    if(.not.valid) return
    ! Check before undoing column scaling, using the original equilibrated
    ! equations rather than the factorized matrix. Normalize each dot product
    ! so neither its absolute sum nor its residual needs to overflow.
    do i=1,n
      do j=1,n
        if(abs(scaled(i,j))>1) then
          if(abs(answer(j))>huge(1._dp)/abs(scaled(i,j))) return
        endif
      enddo
      terms=scaled(i,:)*answer
      normalizer=max(abs(rhs(i)),maxval(abs(terms)),tiny(1._dp))
      terms=terms/normalizer;normalized_rhs=rhs(i)/normalizer
      error=abs(sum(terms)-normalized_rhs)
      magnitude=sum(abs(terms))+abs(normalized_rhs)
      if(error>100*epsilon(1._dp)*magnitude) return
    enddo
    do j=1,n
      if(column_scale(j)<1) then
        if(abs(answer(j))>huge(1._dp)*column_scale(j)) return
      endif
      answer(j)=answer(j)/column_scale(j)
    enddo
    if(any(.not.ieee_is_finite(answer))) return
    solution=answer;ok=.true.
  end subroutine solve_equilibrated_dense_linear_system


  ! The inverse's reference bounds are linear-extension anchors when the
  ! pinned EOS is enabled, not a temperature clamp on a physical rate query.
  ! Its out-of-bound fixed-energy derivative uses anchor Cv and affine species
  ! energies e_i(anchor)+(T-anchor)*Cv_i(anchor), matching that same inverse.
  subroutine reactor_temperature_difference_interval(species,temperature,step,lower,upper,reference,ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: temperature,step
    real(dp), intent(out) :: lower,upper,reference
    logical, intent(out) :: ok
    real(dp) :: anchor_lower,anchor_upper
    lower=0;upper=0;reference=temperature;ok=.false.
    call nasa7_eos_bounds(species,anchor_lower,anchor_upper,ok)
    if(.not.ok) return
    ok=.false.
    if(.not.all(ieee_is_finite([temperature,step])).or.min(temperature,step)<=0) return
    if(temperature>huge(1._dp)-step) return
    if(species(1)%pinned_eos_extrapolation) then
      lower=max(tiny(1._dp),temperature-step);upper=temperature+step
      reference=min(anchor_upper,max(anchor_lower,temperature))
    else
      lower=max(anchor_lower,temperature-step);upper=min(anchor_upper,temperature+step)
    endif
    ok=upper>lower
  end subroutine reactor_temperature_difference_interval

  subroutine reactor_linearized_species_energy(energy,cv,delta_temperature,result,ok)
    real(dp), intent(in) :: energy,cv,delta_temperature
    real(dp), intent(out) :: result
    logical, intent(out) :: ok
    real(dp) :: increment
    result=energy;ok=.true.
    if(delta_temperature==0) return
    ok=.false.
    if(.not.all(ieee_is_finite([energy,cv,delta_temperature]))) return
    if(abs(delta_temperature)>1) then
      if(abs(cv)>huge(1._dp)/abs(delta_temperature)) return
    endif
    increment=cv*delta_temperature
    call reactor_energy_difference(energy,-increment,result,ok)
  end subroutine reactor_linearized_species_energy

  subroutine reactor_energy_difference(left,right,result,ok)
    real(dp), intent(in) :: left,right
    real(dp), intent(out) :: result
    logical, intent(out) :: ok
    result=0;ok=.false.
    if(.not.all(ieee_is_finite([left,right]))) return
    if(right<0) then
      if(left>huge(1._dp)+right) return
    else if(right>0) then
      if(left< -huge(1._dp)+right) return
    endif
    result=left-right;ok=ieee_is_finite(result)
  end subroutine reactor_energy_difference

end module constant_volume_reactor_mod
