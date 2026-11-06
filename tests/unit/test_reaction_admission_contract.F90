program test_reaction_admission_contract
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use elementary_kinetics_mod, only: elementary_reaction
  use constant_volume_reactor_mod, only: check_forced_reactor_admission, &
    reactor_specific_internal_energy,advance_forced_reactor_interval, &
    reaction_admission_not_attempted,reaction_admission_accepted,reaction_admission_rejected, &
    reaction_admission_unknown,reaction_admission_provider_refused
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(elementary_reaction) :: inert(0)
  real(dp) :: original(10),forcing(10),saved(10),saved_forcing(10),candidate(10)
  real(dp) :: temperature,energy,external_energy,source_energy,accepted_temperature,dt,old_gas,new_gas
  real(dp) :: gas(15),gas_rate(15),gas_endpoint(15)
  integer :: status,steps,labels(5),i,j
  logical :: ok
  call load_h2o2_full_thermo(species,ok);call require(ok,'full thermo')
  call load_h2o2_full_mechanism(reactions,ok);call require(ok,'full mechanism')
  call require(size(species)==10,'packet species count')
  call require(all(species%name==[character(len=16)::'H2','H','O','O2','OH','H2O','HO2','H2O2','AR','N2']), &
    'packet species order')
  species%pinned_eos_extrapolation=.true.
  labels=[reaction_admission_not_attempted,reaction_admission_accepted,reaction_admission_rejected, &
    reaction_admission_unknown,reaction_admission_provider_refused]
  do i=1,5
    do j=i+1,5
      call require(labels(i)/=labels(j),'distinct status meanings')
    enddo
  enddo

  ! Whole-interval inert acceptance is a copied-state preview, not publication.
  original=0;original(10)=1;forcing=0;temperature=1000;dt=1.e-8_dp
  energy=reactor_specific_internal_energy(species,original,temperature,ok)
  call require(ok,'inert initial caloric energy')
  saved=original;saved_forcing=forcing
  call check_forced_reactor_admission(species,inert,energy,forcing,0._dp,dt,1.e-8_dp,1.e-12_dp, &
    original,temperature,status,maximum_steps=1)
  call require(status==reaction_admission_accepted,'exact inert interval certifies within one step')
  call require(all(original==saved).and.all(forcing==saved_forcing).and.temperature==1000, &
    'accepted preview leaves caller packet untouched')

  call check_forced_reactor_admission(species,inert,energy,forcing,0._dp,dt,1.e-8_dp,1.e-12_dp, &
    original,temperature,status,maximum_steps=0)
  call require(status==reaction_admission_unknown,'zero explicit work budget never adopts the default unbounded work')
  call require(all(original==saved).and.all(forcing==saved_forcing).and.temperature==1000, &
    'invalid work budget preserves packet')

  ! Actual phase2 physical fine donor (node2,4,1,2), not the imposed sink.
  ! Preserve the different old held and old gas energies. These 17-digit
  ! constants come from complete actual caller assembly, not a tuned fixture.
  ! cell-rank0.dat SHA edee78405c6f1d15a8606c9f1b40950b0db98c0fd7129aaf03900a6818180b4f
  ! assembly-rank0.dat SHA 12e7f58394d71aab827cb39d739644bb9eb82629060bdd9e4641dfc6ce023e9f
  original=1.00000000000000006e-9_dp
  original(1)=2.80000000000000006e-2_dp;original(4)=2.24000000000000005e-1_dp
  original(10)=7.47999992999999974e-1_dp
  forcing=-7.49411718181922509e-2_dp
  forcing(1)=-2.09835281090938300e6_dp;forcing(4)=-1.67868224872750640e7_dp
  forcing(10)=-5.60559959954196066e7_dp
  energy=1.09727798565058387e6_dp;temperature=1.40000000000000023e3_dp
  gas(:5)=[1._dp,1._dp,2._dp,3._dp,1.59841135239832290e6_dp];gas(6:)=original
  gas_rate(:5)=[-7.49411718181922585e7_dp,-7.49411718181918561e7_dp, &
    -1.49882343636384517e8_dp,-2.24823515454576761e8_dp,-1.19788009731249406e14_dp]
  gas_rate(6:)=forcing;gas_endpoint=gas+dt*gas_rate
  old_gas=gas(5)-.5_dp*sum(gas(2:4)**2)/gas(1)
  new_gas=gas_endpoint(5)-.5_dp*sum(gas_endpoint(2:4)**2)/gas_endpoint(1)
  external_energy=(new_gas-old_gas)/dt
  call require(old_gas/=energy,'independent held and gas initial energies retained')
  call require(external_energy==-1.19787485143046688e14_dp,'actual gas-projected packet source identity')
  call require(all(original+dt*forcing>0),'actual affine species endpoint is positive')
  saved=original;saved_forcing=forcing;source_energy=external_energy
  call check_forced_reactor_admission(species,reactions,energy,forcing,external_energy,dt,1.e-8_dp,1.e-12_dp, &
    original,temperature,status,maximum_steps=512)
  call require(status==reaction_admission_unknown,'original full donor work limit is unknown, not infeasible')
  call require(all(original==saved).and.all(forcing==saved_forcing).and.external_energy==source_energy.and. &
    temperature==1.40000000000000023e3_dp.and.energy==1.09727798565058387e6_dp, &
    'unknown preview preserves exact actual packet')

  ! A negative affine trace is replenished by chemistry; preview cannot use
  ! affine species sign or an exhausted small work budget as impossibility.
  original=0;original(10)=1;temperature=1000.0000000000001_dp
  energy=469600.82041824982_dp;external_energy=-34759358186.641472_dp
  forcing=[-8.8294917300898273e-27_dp,-1.4683862885727810e-20_dp,-2.8938113755454018e-18_dp, &
    -2.9736855926687663e-20_dp,-3.5592099773271806e-7_dp,899.85617208188751_dp, &
    -2.9820151528746120e-12_dp,2750.1162925888539_dp,0._dp,-7621.7142497948162_dp]
  call require(any(original+dt*forcing<0),'replenishment negative-affine control active')
  saved=original;saved_forcing=forcing
  call check_forced_reactor_admission(species,reactions,energy,forcing,external_energy,dt,2.e-7_dp,1.e-12_dp, &
    original,temperature,status,maximum_steps=1)
  call require(status==reaction_admission_unknown,'insufficient replenishment work is unknown')
  call check_forced_reactor_admission(species,reactions,energy,forcing,external_energy,dt,2.e-7_dp,1.e-12_dp, &
    original,temperature,status,maximum_steps=512)
  call require(status==reaction_admission_accepted,'signed replenishment remains certifiable')
  call require(all(original==saved).and.all(forcing==saved_forcing).and.temperature==1000.0000000000001_dp, &
    'both replenishment previews are immutable')
  candidate=original;accepted_temperature=temperature
  call advance_forced_reactor_interval(species,reactions,energy,forcing,external_energy,dt,2.e-7_dp,1.e-12_dp, &
    candidate,accepted_temperature,steps,ok,maximum_steps=512)
  call require(ok.and.all(candidate>=0),'certified replenishment has real physical endpoint')

  ! Even a legacy refusal backed by a structural conserved-pool certificate
  ! is conservatively UNKNOWN in this built-in helper, which exposes no
  ! certificate. Never infer a stronger reason solely from boolean false.
  original=0;original(10)=1;forcing=0;forcing(9)=-1.e-10_dp;temperature=1000
  saved=original;saved_forcing=forcing
  call check_forced_reactor_admission(species,reactions,energy,forcing,0._dp,dt,2.e-7_dp,1.e-12_dp, &
    original,temperature,status,maximum_steps=512)
  call require(status==reaction_admission_unknown,'plain builtin refusal provides no certified status')
  call require(all(original==saved).and.all(forcing==saved_forcing).and.temperature==1000, &
    'refused preview leaves source and state untouched')
  print *, 'PASS: immutable admission status, actual donor energy packet, bounded unknown and signed replenishment'
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    endif
  end subroutine require
end program test_reaction_admission_contract
