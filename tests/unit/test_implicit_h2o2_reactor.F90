program test_implicit_h2o2_reactor
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species,universal_gas_constant
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density, mixture_mass_properties
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_full_thermo_mod, only: full_nspecies, load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use constant_volume_reactor_mod, only: &
    advance_constant_volume_implicit_adaptive,forced_reactor_rhs,reactor_rhs,forced_reactor_reduced_jacobian, &
    advance_forced_reactor_implicit_adaptive,advance_forced_reactor_interval,reaction_source_heat_release, &
    backward_euler_trial
  use reactive_3d_mod, only: react_sdc_state_3d
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp) :: x(full_nspecies), y(full_nspecies), density, temperature
  real(dp) :: target_energy, final_energy, molecular_weight, gas_constant
  real(dp) :: cp, cv, gamma, enthalpy, entropy, accepted, next_step
  integer :: newton_iterations, rejected,steps,trace_index
  logical :: ok
  type(elementary_reaction), allocatable :: inert(:)
  real(dp) :: rhoy(full_nspecies),source(full_nspecies),rates(full_nspecies),baseline(full_nspecies),source_t
  real(dp) :: forced_t,forced_accepted,forced_next,saved(full_nspecies),trace_energy

  call check_forced_jacobian()
  call check_inactive_mechanism()
  call load_h2o2_full_thermo(species, ok)
  if (.not. ok) error stop "Failed to load full thermo"
  call load_h2o2_full_mechanism(reactions, ok)
  if (.not. ok) error stop "Failed to load full mechanism"
  x = [2.0_dp, 1.0e-6_dp, 1.0e-12_dp, 1.0_dp, 1.0e-12_dp, &
    0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 3.0_dp]
  x = x / sum(x)
  call mass_fractions_from_mole_fractions(species, x, y, ok)
  if (.not. ok) error stop "Composition conversion failed"
  temperature = 1000.0_dp
  density = mixture_density(species, y, 101325.0_dp, temperature, ok)
  if (.not. ok) error stop "Density evaluation failed"
  call mixture_mass_properties(species, y, temperature, molecular_weight, &
    gas_constant, cp, cv, gamma, enthalpy, target_energy, entropy, ok)
  if (.not. ok) error stop "Initial energy evaluation failed"
  allocate(inert(0));rhoy=density*y;source=0
  call reactor_rhs(species,reactions,density,target_energy,y,temperature,baseline,source_t,ok)
  if(.not.ok) error stop 'unforced RHS baseline'
  call forced_reactor_rhs(species,reactions,density*target_energy,0._dp,rhoy,source,0._dp, &
    temperature,rates,source_t,ok)
  if(.not.ok.or.maxval(abs(rates-density*baseline)/max(1._dp,abs(density*baseline)))>1e-12_dp) &
    error stop 'forced RHS zero-source chemistry limit'
  source=.2_dp*y;rhoy=(density+.05_dp)*y
  call forced_reactor_rhs(species,inert,density*target_energy,.25_dp,rhoy,source,.2_dp*target_energy, &
    temperature,rates,source_t,ok)
  if(.not.ok.or.any(rates/=source).or.abs(source_t-temperature)>1e-8_dp) &
    error stop 'forced RHS variable density and energy'
  call reactor_rhs(species,reactions,density+.05_dp,target_energy,y,temperature,baseline,source_t,ok)
  if(.not.ok) error stop 'forced chemistry comparison baseline'
  call forced_reactor_rhs(species,reactions,density*target_energy,.25_dp,rhoy,source,.2_dp*target_energy, &
    temperature,rates,source_t,ok)
  if(.not.ok.or.maxval(abs(rates-((density+.05_dp)*baseline+source))/ &
    max(1._dp,abs((density+.05_dp)*baseline+source)))>1e-12_dp) error stop 'chemistry and external species rates'
  rhoy(1)=-1
  call forced_reactor_rhs(species,reactions,density*target_energy,.25_dp,rhoy,source,.2_dp*target_energy, &
    temperature,rates,source_t,ok)
  if(ok.or.any(rates/=0).or.source_t/=0) error stop 'forced RHS invalid species zero outputs'
  rhoy=density*y;source=.2_dp*y;forced_t=temperature
  call advance_forced_reactor_implicit_adaptive(species,inert,density*target_energy,0._dp,source, &
    .2_dp*target_energy,.25_dp,1e-6_dp,1e-12_dp,rhoy,forced_t,forced_accepted,forced_next,ok)
  if(.not.ok.or.forced_accepted/=.25_dp.or.maxval(abs(rhoy-(density+.05_dp)*y))>1e-14_dp.or. &
    abs(forced_t-temperature)>1e-8_dp) error stop 'forced inert integration'
  saved=rhoy
  call advance_forced_reactor_implicit_adaptive(species,inert,density*target_energy,0._dp,source, &
    .2_dp*target_energy,-.25_dp,1e-6_dp,1e-12_dp,rhoy,forced_t,forced_accepted,forced_next,ok)
  if(ok.or.any(rhoy/=saved).or.forced_accepted/=0.or.forced_next/=0) error stop 'forced step rollback'
  block
    type(elementary_reaction) :: inactive(1)
    real(dp) :: tracer(full_nspecies),forcing(full_nspecies),tracer_t,tracer_energy
    ! A nonempty mechanism exercises Newton even when chemistry is inactive.
    inactive(1)=reactions(1)
    inactive(1)%forward_rate%pre_exponential=0._dp
    tracer=0;tracer(full_nspecies)=1;tracer_t=300
    call mixture_mass_properties(species,tracer,tracer_t,molecular_weight,gas_constant, &
      cp,cv,gamma,enthalpy,tracer_energy,entropy,ok)
    if(.not.ok) error stop 'trace forcing initial energy'
    forcing=0;forcing(1)=1e-10_dp
    call advance_forced_reactor_implicit_adaptive(species,inactive,tracer_energy,0._dp,forcing, &
      0._dp,1e-6_dp,1e-6_dp,1e-12_dp,tracer,tracer_t,forced_accepted,forced_next,ok)
    if(.not.ok.or.forced_accepted/=1e-6_dp) error stop 'trace forcing adaptive step'
    if(abs(tracer(1)-forcing(1)*forced_accepted)>1e-28_dp.or.any(tracer<0)) &
      error stop 'sub-atol affine tracer must survive Newton convergence'
  end block
  block
    type(elementary_reaction) :: slow(1)
    real(dp) :: tracer(full_nspecies),forcing(full_nspecies),initial_rates(full_nspecies)
    real(dp) :: tracer_t,tracer_energy,rhs_t,expected
    ! Balanced O recombination creates O2 below atol, with no affine forcing.
    slow(1)=reactions(1)
    slow(1)%forward_rate%pre_exponential=1e-4_dp
    slow(1)%reversible=.false.
    tracer=0;tracer(3)=.1_dp;tracer(full_nspecies)=.9_dp;tracer_t=300
    call mixture_mass_properties(species,tracer,tracer_t,molecular_weight,gas_constant, &
      cp,cv,gamma,enthalpy,tracer_energy,entropy,ok)
    if(.not.ok) error stop 'trace chemistry initial energy'
    forcing=0
    call forced_reactor_rhs(species,slow,tracer_energy,0._dp,tracer,forcing,0._dp, &
      tracer_t,initial_rates,rhs_t,ok)
    if(.not.ok.or.initial_rates(4)<=0) error stop 'trace chemistry production baseline'
    call advance_forced_reactor_implicit_adaptive(species,slow,tracer_energy,0._dp,forcing, &
      0._dp,1e-6_dp,1e-6_dp,1e-12_dp,tracer,tracer_t,forced_accepted,forced_next,ok)
    if(.not.ok.or.forced_accepted/=1e-6_dp) error stop 'trace chemistry adaptive step'
    expected=forced_accepted*initial_rates(4)
    if(expected>=1e-12_dp.or.abs(tracer(4)-expected)>1e-8_dp*expected.or.any(tracer<0)) &
      error stop 'sub-atol chemical product must survive Newton convergence'
    ! A correction smaller than every populated component's ulp is resolved.
    slow(1)%forward_rate%pre_exponential=1e-100_dp
    tracer=0;tracer(3)=.1_dp;tracer(4)=.1_dp;tracer(full_nspecies)=.8_dp;tracer_t=300
    saved=tracer
    call mixture_mass_properties(species,tracer,tracer_t,molecular_weight,gas_constant, &
      cp,cv,gamma,enthalpy,tracer_energy,entropy,ok)
    if(.not.ok) error stop 'unrepresentable correction initial energy'
    call advance_forced_reactor_implicit_adaptive(species,slow,tracer_energy,0._dp,forcing, &
      0._dp,1e-6_dp,1e-6_dp,1e-12_dp,tracer,tracer_t,forced_accepted,forced_next,ok)
    if(.not.ok.or.forced_accepted/=1e-6_dp.or.any(tracer/=saved)) &
      error stop 'unrepresentable Newton correction must not reject a converged step'
  end block
  block
    real(dp) :: cell_species(full_nspecies),cell_source(full_nspecies),cell_temperature,cell_energy
    real(dp) :: rho_before,rho_after,expected_water,element_actual(3),element_expected(3)
    real(dp) :: atoms(3,full_nspecies)
    ! Extracted selected MPI EB recursive-SDC cell: N2 removal plus water injection.
    cell_species=0;cell_species(full_nspecies)=.93000000000000016_dp
    cell_source=0;cell_source(6)=7.4339030865166767e-5_dp
    cell_source(full_nspecies)=-6.0195696059402293_dp
    cell_temperature=489.09335141132806_dp;rho_before=sum(cell_species)
    cell_energy=50859.780854659184_dp
    call advance_forced_reactor_interval(species,reactions,cell_energy,cell_source, &
      80690783.882164389_dp,1e-8_dp,2e-7_dp,9.3e-13_dp,cell_species,cell_temperature,steps,ok)
    if(.not.ok) error stop 'N2-water forced cell must complete its interval'
    rho_after=rho_before+1e-8_dp*sum(cell_source)
    expected_water=1e-8_dp*cell_source(6)
    if(any(cell_species<0).or.abs(sum(cell_species)-rho_after)>3e-15_dp*rho_after) &
      error stop 'N2-water forced cell positivity and density balance'
    if(abs(cell_species(6)-expected_water)>1e-10_dp*expected_water) &
      error stop 'N2-water injected sub-atol component must survive'
    atoms(1,:)=[2._dp,1._dp,0._dp,0._dp,1._dp,2._dp,1._dp,2._dp,0._dp,0._dp]
    atoms(2,:)=[0._dp,0._dp,1._dp,2._dp,1._dp,1._dp,2._dp,2._dp,0._dp,0._dp]
    atoms(3,:)=0;atoms(3,full_nspecies)=2
    element_actual=matmul(atoms,cell_species/species%molecular_weight)
    element_expected=atoms(:,6)*expected_water/species(6)%molecular_weight+ &
      atoms(:,full_nspecies)*(rho_before+1e-8_dp*cell_source(full_nspecies))/species(full_nspecies)%molecular_weight
    if(any(abs(element_actual-element_expected)>1e-10_dp*element_expected)) &
      error stop 'N2-water forced cell H-O-N inventory'
    call mixture_mass_properties(species,cell_species/sum(cell_species),cell_temperature,molecular_weight, &
      gas_constant,cp,cv,gamma,enthalpy,final_energy,entropy,ok)
    if(.not.ok.or.abs(sum(cell_species)*final_energy-(cell_energy+1e-8_dp*80690783.882164389_dp))> &
      5e-10_dp*abs(cell_energy)) error stop 'N2-water forced cell energy balance'
  end block
  block
    real(dp) :: initial(full_nspecies),reference(full_nspecies),finish(full_nspecies),rhs(full_nspecies)
    real(dp) :: dt,rho,energy,finished_t,rhs_t,forcing(full_nspecies),scale
    integer :: k,iterations
    ! Water in nitrogen produces several coupled intermediates far below atol.
    ! A single Newton correction can resolve the first product while erasing a
    ! later one. Check the actual backward-Euler chemical balance of each trace.
    dt=1e-8_dp;initial=0;initial(full_nspecies)=.93000000000000016_dp
    forcing=0;forcing(6)=7.4339030865166767e-5_dp
    forcing(full_nspecies)=-6.0195696059402293_dp
    rho=sum(initial)+dt*sum(forcing)
    reference=(initial+dt*forcing)/rho
    energy=(50859.780854659184_dp+dt*80690783.882164389_dp)/rho
    call backward_euler_trial(species,reactions,rho,energy,reference,489.09335141132806_dp,dt, &
      2e-7_dp,9.3e-13_dp/rho,finish,finished_t,iterations,ok,reference)
    if(.not.ok.or.any(finish<0)) error stop 'coupled sub-atol trace trial must complete'
    call reactor_rhs(species,reactions,rho,energy,finish,finished_t,rhs,rhs_t,ok)
    if(.not.ok.or.finish(8)<=0) error stop 'coupled H2O2 trace must be represented'
    do k=1,full_nspecies-1
      if(k==6) cycle
      scale=max(abs(finish(k)),abs(reference(k)),abs(dt*rhs(k)),tiny(1._dp))
      if(abs(finish(k)-reference(k)-dt*rhs(k))>3e-7_dp*scale) &
        error stop 'sub-atol coupled chemical balance'
    enddo
  end block
  block
    type(nasa7_species) :: caloric(2),bad_species(2)
    real(dp) :: qdot,expected,temps(3),cp_coefficient,production(2)
    integer :: k
    caloric=species(1:2)
    caloric(1)%molecular_weight=2;caloric(2)%molecular_weight=4
    caloric%temperature_min=200;caloric%temperature_mid=1000;caloric%temperature_max=5000
    do k=1,2
      caloric(k)%low_coefficients=0;caloric(k)%high_coefficients=0
      caloric(k)%low_coefficients(1)=3.5_dp;caloric(k)%high_coefficients(1)=4.5_dp
      caloric(k)%low_coefficients(6)=100._dp-600._dp*(k-1)
      caloric(k)%high_coefficients(6)=caloric(k)%low_coefficients(6)
    enddo
    temps=[300._dp,nearest(1000._dp,-1._dp),1000._dp]
    production=[-1._dp,1._dp]
    do k=1,size(temps)
      cp_coefficient=3.5_dp
      if(temps(k)>=1000._dp) cp_coefficient=4.5_dp
      expected=universal_gas_constant*((cp_coefficient*temps(k)+100._dp)/2- &
        (cp_coefficient*temps(k)-500._dp)/4)
      qdot=17
      call reaction_source_heat_release(caloric,temps(k),production,qdot,ok)
      if(.not.ok.or.abs(qdot-expected)>8*epsilon(1._dp)*abs(expected)) &
        error stop 'reaction heat release SI enthalpy sum and midpoint branch'
      call reaction_source_heat_release(caloric,temps(k),-production,qdot,ok)
      if(.not.ok.or.abs(qdot+expected)>8*epsilon(1._dp)*abs(expected)) &
        error stop 'endothermic reaction heat release keeps negative sign'
    enddo
    call reaction_source_heat_release(caloric,300._dp,[0._dp,0._dp],qdot,ok)
    if(.not.ok.or.qdot/=0) error stop 'zero chemical source gives zero heat release'
    qdot=17
    call reaction_source_heat_release(caloric,300._dp,production(:1),qdot,ok)
    if(ok.or.qdot/=17) error stop 'heat release source-shape rollback'
    call reaction_source_heat_release(caloric,ieee_value(0._dp,ieee_quiet_nan),production,qdot,ok)
    if(ok.or.qdot/=17) error stop 'heat release temperature rollback'
    bad_species=caloric;bad_species(2)%molecular_weight=0
    call reaction_source_heat_release(bad_species,300._dp,production,qdot,ok)
    if(ok.or.qdot/=17) error stop 'heat release late species failure rollback'
    call reaction_source_heat_release(caloric,300._dp,[huge(1._dp),0._dp],qdot,ok)
    if(ok.or.qdot/=17) error stop 'heat release product overflow rejects without trap'
  end block
  rhoy=density*y;source=0;forced_t=temperature
  call advance_forced_reactor_implicit_adaptive(species,reactions,density*target_energy,0._dp,source, &
    0._dp,1e-7_dp,1e-6_dp,density*1e-12_dp,rhoy,forced_t,forced_accepted,forced_next,ok)
  if(.not.ok) error stop 'forced zero-source integration'
  call advance_constant_volume_implicit_adaptive( &
    species, reactions, density, target_energy, 1.0e-7_dp, 1.0e-6_dp, &
    1.0e-12_dp, y, temperature, accepted, next_step, newton_iterations, &
    rejected, ok)
  if (.not. ok) error stop "Implicit adaptive step failed"
  if(abs(forced_accepted-accepted)>1e-15_dp.or.maxval(abs(rhoy/density-y))>1e-10_dp.or. &
    abs(forced_t-temperature)>1e-6_dp) error stop 'forced zero-source integrator limit'
  if (accepted <= 0.0_dp .or. next_step <= 0.0_dp) then
    error stop "Implicit adaptive step returned invalid time steps"
  end if
  if (minval(y) < -1.0e-13_dp .or. abs(sum(y) - 1.0_dp) > 1.0e-12_dp) then
    error stop "Implicit step violated composition bounds"
  end if
  call mixture_mass_properties(species, y, temperature, molecular_weight, &
    gas_constant, cp, cv, gamma, enthalpy, final_energy, entropy, ok)
  if (.not. ok) error stop "Final energy evaluation failed"
  if (abs(final_energy - target_energy) / max(1.0_dp, abs(target_energy)) &
      > 5.0e-10_dp) then
    error stop "Implicit step violated constant-volume energy"
  end if
  rhoy=density*y;source=.2_dp*y;forced_t=temperature
  call advance_forced_reactor_implicit_adaptive(species,reactions,density*target_energy,0._dp,source, &
    .2_dp*target_energy,1e-7_dp,1e-6_dp,1e-12_dp,rhoy,forced_t,forced_accepted,forced_next,ok)
  if(.not.ok.or.forced_accepted<=0) error stop 'forced reacting integration'
  if(abs(sum(rhoy)-(density+.2_dp*forced_accepted))>1e-13_dp) error stop 'forced reacting mass inventory'
  call mixture_mass_properties(species,rhoy/sum(rhoy),forced_t,molecular_weight,gas_constant, &
    cp,cv,gamma,enthalpy,final_energy,entropy,ok)
  if(.not.ok.or.abs(sum(rhoy)*final_energy-(density+.2_dp*forced_accepted)*target_energy)/ &
    max(1._dp,abs(density*target_energy))>5e-10_dp) error stop 'forced reacting energy inventory'
  rhoy=density*y;source=.2_dp*y;forced_t=temperature;saved=rhoy
  call advance_forced_reactor_interval(species,inert,density*target_energy,source,.2_dp*target_energy, &
    .25_dp,1e-6_dp,1e-12_dp,rhoy,forced_t,steps,ok,initial_time_step=.03125_dp,maximum_steps=1)
  if(ok.or.steps/=1.or.any(rhoy/=saved).or.forced_t/=temperature) &
    error stop 'forced interval rollback after accepted substep'
  call advance_forced_reactor_interval(species,inert,density*target_energy,source,.2_dp*target_energy, &
    .25_dp,1e-6_dp,1e-12_dp,rhoy,forced_t,steps,ok,initial_time_step=.03125_dp)
  if(.not.ok.or.steps/=4.or.maxval(abs(rhoy-(density+.05_dp)*y))>1e-14_dp.or. &
    abs(forced_t-temperature)>1e-8_dp) error stop 'forced interval inert elapsed-time accounting'
  ! No-reaction full intervals must not invent an SDC chemical remainder
  ! at a newly introduced trace species. Verify the exact linear increment
  ! and rejection/rollback when the same zero species is instead depleted.
  rhoy=density*y;trace_index=minloc(y,dim=1);rhoy(trace_index)=0
  call mixture_mass_properties(species,rhoy/sum(rhoy),temperature,molecular_weight,gas_constant, &
    cp,cv,gamma,enthalpy,trace_energy,entropy,ok)
  if(.not.ok) error stop 'trace forced interval properties'
  source=0;source(trace_index)=1e-46_dp;forced_t=temperature;saved=rhoy
  call advance_forced_reactor_interval(species,inert,sum(saved)*trace_energy,source,0._dp, &
    1e-8_dp,1e-6_dp,1e-12_dp,rhoy,forced_t,steps,ok)
  if(.not.ok.or.steps/=1.or.any(rhoy/=saved+1e-8_dp*source)) &
    error stop 'trace forced interval exact inert solution'
  block
    real(dp) :: old_cell(5+full_nspecies,1,1,1),new_cell(5+full_nspecies,1,1,1)
    real(dp) :: cell_rate(5+full_nspecies,1,1,1),old_t(1,1,1),new_t(1,1,1),ir(full_nspecies+1,1,1,1)
    old_cell=0;old_cell(1,1,1,1)=sum(saved);old_cell(5,1,1,1)=sum(saved)*trace_energy
    old_cell(6:,1,1,1)=saved;old_t=temperature;new_cell=old_cell;new_t=old_t;ir=53
    cell_rate=0;cell_rate(1,1,1,1)=sum(source);cell_rate(6:,1,1,1)=source
    call react_sdc_state_3d(species,inert,old_cell,old_t,cell_rate,1e-8_dp,1e-6_dp,1e-12_dp, &
      new_cell,new_t,ir,ok)
    if(.not.ok.or.any(ir/=0).or.any(new_cell(6:,1,1,1)/=saved+1e-8_dp*source)) &
      error stop 'inert SDC has exactly zero chemical species and energy sources'
  end block
  rhoy=saved;source(trace_index)=-1e-46_dp;forced_t=temperature
  call advance_forced_reactor_interval(species,inert,sum(saved)*trace_energy,source,0._dp, &
    1e-8_dp,1e-6_dp,1e-12_dp,rhoy,forced_t,steps,ok)
  if(ok.or.steps/=0.or.any(rhoy/=saved).or.forced_t/=temperature) &
    error stop 'trace forced interval depletion rollback'
  rhoy=density*y;source=.2_dp*y;forced_t=temperature
  call advance_forced_reactor_interval(species,reactions,density*target_energy,source,.2_dp*target_energy, &
    1e-6_dp,1e-6_dp,1e-12_dp,rhoy,forced_t,steps,ok,initial_time_step=1e-8_dp)
  if(.not.ok.or.steps<2) error stop 'forced reacting full interval'
  if(abs(sum(rhoy)-(density+.2_dp*1e-6_dp))>1e-13_dp) error stop 'forced interval mass inventory'
  call mixture_mass_properties(species,rhoy/sum(rhoy),forced_t,molecular_weight,gas_constant, &
    cp,cv,gamma,enthalpy,final_energy,entropy,ok)
  if(.not.ok.or.abs(sum(rhoy)*final_energy-(density+.2_dp*1e-6_dp)*target_energy)/ &
    max(1._dp,abs(density*target_energy))>5e-10_dp) error stop 'forced interval energy inventory'
  write(*, '(a)') "test_implicit_h2o2_reactor: PASS"
contains
  subroutine check_inactive_mechanism()
    type(nasa7_species), allocatable :: table(:)
    type(elementary_reaction), allocatable :: mechanism(:),empty(:)
    real(dp) :: state(full_nspecies),control(full_nspecies),forcing(full_nspecies),t,tc
    integer :: n,idx,count,control_count
    logical :: good
    call load_h2o2_full_thermo(table,good)
    call require(good,'inactive mechanism thermo')
    call load_h2o2_full_mechanism(mechanism,good)
    call require(good,'inactive mechanism rates')
    allocate(empty(0));table%pinned_eos_extrapolation=.true.
    state=0;forcing=0;n=0
    do idx=1,size(table)
      if(table(idx)%name=='N2') n=idx
    enddo
    call require(n>0,'inactive bath species exists')
    state(n)=1;forcing(n)=-67777554.259668991_dp;control=state
    t=700;tc=t
    call advance_forced_reactor_interval(table,empty,218595.25696862090_dp,forcing,-31829032696484.051_dp, &
      1.e-8_dp,1.e-8_dp,1.e-12_dp,control,tc,control_count,good)
    call require(good.and.control_count==1,'affine inert endpoint admitted by existing EOS policy')
    call advance_forced_reactor_interval(table,mechanism,218595.25696862090_dp,forcing,-31829032696484.051_dp, &
      1.e-8_dp,1.e-8_dp,1.e-12_dp,state,t,count,good)
    call require(good.and.count==1,'structurally inactive full mechanism takes exact interval')
    call require(all(state==control).and.t==tc,'full inactive mechanism equals independent empty-mechanism control')
    state=0;state(n)=1;t=700;forcing(n)=-2.e8_dp
    call advance_forced_reactor_interval(table,mechanism,218595.25696862090_dp,forcing,0._dp, &
      1.e-8_dp,1.e-8_dp,1.e-12_dp,state,t,count,good)
    call require(.not.good.and.count==0.and.state(n)==1.and.sum(state)==1.and.t==700, &
      'inactive negative endpoint refuses transactionally')
  end subroutine check_inactive_mechanism
  subroutine check_forced_jacobian()
    type(nasa7_species), allocatable :: thermo(:)
    type(elementary_reaction), allocatable :: mechanism(:)
    type(elementary_reaction) :: empty_mechanism(0)
    real(dp), allocatable :: state(:),source(:),trial(:),plus(:),minus(:),jac(:,:),control(:,:)
    integer, allocatable :: independent(:)
    real(dp) :: mw,r,cp,cv,gamma,h,e,s,t,delta,error
    integer :: n,dep,column,idx
    logical :: good
    call load_h2o2_full_thermo(thermo,good)
    call require(good,'forced Jacobian thermo')
    call load_h2o2_full_mechanism(mechanism,good)
    call require(good,'forced Jacobian mechanism')
    n=size(thermo)
    allocate(state(n),source(n),trial(n),plus(n),minus(n),jac(n-1,n-1),control(n-1,n-1))
    state=real([(idx,idx=1,n)],dp);state=2.3_dp*state/sum(state)
    source=.01_dp*state;source(1)=-.002_dp
    call mixture_mass_properties(thermo,state/sum(state),1200._dp,mw,r,cp,cv,gamma,h,e,s,good)
    call require(good,'forced Jacobian initial energy')
    e=e*sum(state)
    call forced_reactor_reduced_jacobian(thermo,empty_mechanism,e,7.e-5_dp,state,source,1.e6_dp, &
      1200._dp,jac,t,good,1)
    call require(good.and.all(jac==0),'inert forced Jacobian exactly zero')
    call forced_reactor_reduced_jacobian(thermo,empty_mechanism,e,7.e-5_dp,state,source,1.e6_dp, &
      1200._dp,jac,t,good,n+1)
    call require(.not.good.and.all(jac==0),'invalid forced closure rejected before inert Jacobian')
    do dep=1,n,n-1
      independent=pack([(idx,idx=1,n)],[(idx/=dep,idx=1,n)])
      call forced_reactor_reduced_jacobian(thermo,mechanism,e,7.e-5_dp,state,source,1.e6_dp, &
        1200._dp,jac,t,good,dep)
      call require(good,'forced Jacobian analytic density closure')
      do column=1,n-1
        idx=independent(column);delta=1.e-6_dp*min(state(idx),state(dep))
        trial=state;trial(idx)=trial(idx)+delta;trial(dep)=trial(dep)-delta
        call forced_reactor_rhs(thermo,mechanism,e,7.e-5_dp,trial,source,1.e6_dp,t,plus,s,good)
        call require(good,'forced Jacobian upper difference')
        trial=state;trial(idx)=trial(idx)-delta;trial(dep)=trial(dep)+delta
        call forced_reactor_rhs(thermo,mechanism,e,7.e-5_dp,trial,source,1.e6_dp,t,minus,s,good)
        call require(good,'forced Jacobian lower difference')
        control(:,column)=(plus(independent)-minus(independent))/(2*delta)
      enddo
      error=maxval(abs(jac-control))/max(1._dp,maxval(abs(control)))
      write(*,*) 'forced Jacobian dependent species/scaled error',dep,error
      call require(error<2.e-5_dp,'forced Jacobian independent finite-difference control')
    enddo
  end subroutine check_forced_jacobian
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if(.not.condition) error stop message
  end subroutine require
end program test_implicit_h2o2_reactor
