program test_equilibrium_reference_convention
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species,universal_gas_constant
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use elementary_kinetics_mod, only: elementary_reaction,valid_elementary_reaction, &
    reaction_equilibrium_constant,elementary_production_jacobian,elementary_production_rates_from_concentrations, &
    reaction_effective_forward_rate
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism,h2o2_h2_index,h2o2_h_index
  implicit none
  type(nasa7_species),allocatable :: species(:)
  type(elementary_reaction),allocatable :: reactions(:)
  type(elementary_reaction) :: probe
  real(dp),allocatable :: c(:),plus(:),minus(:),jac(:,:),cp(:),cm(:)
  real(dp) :: kc,baseline,expected,delta,bad(4),forward_before,forward_after
  integer :: mode,i,n
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok)
  if(.not.ok) error stop 'thermo'
  call load_h2o2_elementary_mechanism(reactions,ok)
  if(.not.ok) error stop 'mechanism'
  n=size(species)
  allocate(c(n));c=.001_dp
  do i=1,size(reactions)
    if(reactions(i)%equilibrium_reference_gas_constant/=universal_gas_constant) error stop 'changed default'
  enddo
  do mode=-1,1
    probe=reactions(1)
    if(mode/=0) then
      probe%reactant_stoich=0;probe%product_stoich=0
      if(mode==1) then
        probe%reactant_stoich(h2o2_h2_index)=1;probe%product_stoich(h2o2_h_index)=2
      else
        probe%reactant_stoich(h2o2_h_index)=2;probe%product_stoich(h2o2_h2_index)=1
      endif
    endif
    baseline=reaction_equilibrium_constant(species,probe,1200._dp,ok)
    if(.not.ok.or.baseline<=0) error stop 'baseline Kc'
    call reaction_effective_forward_rate(probe,1200._dp,c,forward_before,ok=ok)
    if(.not.ok) error stop 'baseline forward rate'
    probe%equilibrium_reference_gas_constant=8314.46_dp
    call reaction_effective_forward_rate(probe,1200._dp,c,forward_after,ok=ok)
    if(.not.ok.or.forward_before/=forward_after) error stop 'forward rate changed'
    kc=reaction_equilibrium_constant(species,probe,1200._dp,ok)
    expected=(universal_gas_constant/8314.46_dp)**mode
    if(.not.ok.or.abs(kc/baseline-expected)>3.e-14_dp) error stop 'reference concentration power'
    if(mode==0.and.kc/=baseline) error stop 'zero net stoichiometry changed'
  enddo
  bad=[0._dp,-1._dp,ieee_value(0._dp,ieee_quiet_nan),ieee_value(0._dp,ieee_positive_inf)]
  do i=1,size(bad)
    probe%equilibrium_reference_gas_constant=bad(i)
    if(valid_elementary_reaction(probe,n)) error stop 'invalid convention admitted'
    kc=reaction_equilibrium_constant(species,probe,1200._dp,ok)
    if(ok.or.kc/=0) error stop 'invalid convention evaluated'
  enddo
  probe%equilibrium_reference_gas_constant=huge(1._dp)
  kc=reaction_equilibrium_constant(species,probe,1200._dp,ok)
  if(ok.or.kc/=0) error stop 'overflowing reference concentration'
  probe%equilibrium_reference_gas_constant=tiny(1._dp)
  kc=reaction_equilibrium_constant(species,probe,1200._dp,ok)
  if(ok.or.kc/=0) error stop 'unrepresentable reference concentration'
  ! The same convention must enter both reverse rates and their density Jacobian.
  do i=1,size(reactions)
    reactions(i)%equilibrium_reference_gas_constant=8314.46_dp
  enddo
  allocate(plus(n),minus(n),jac(n,n),cp(n),cm(n))
  call elementary_production_jacobian(species,reactions,1200._dp,c,jac,ok)
  if(.not.ok) error stop 'Jacobian'
  do i=1,n
    delta=1.e-7_dp*c(i);cp=c;cm=c;cp(i)=cp(i)+delta;cm(i)=cm(i)-delta
    call elementary_production_rates_from_concentrations(species,reactions,1200._dp,cp,plus,ok)
    if(.not.ok) error stop 'plus rates'
    call elementary_production_rates_from_concentrations(species,reactions,1200._dp,cm,minus,ok)
    if(.not.ok) error stop 'minus rates'
    if(maxval(abs((plus-minus)/(2*delta)-jac(:,i)))> &
       2.e-6_dp*max(1._dp,maxval(abs(jac(:,i))))) error stop 'convention Jacobian mismatch'
  enddo
  print *, 'test_equilibrium_reference_convention: PASS'
end program
