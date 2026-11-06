program test_forced_signed_trials
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use nasa7_thermo_mod, only: nasa7_species,nasa7_mass_properties,configure_nasa7_eos_policy,universal_gas_constant
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use elementary_kinetics_mod, only: elementary_reaction,arrhenius_rate
  use constant_volume_reactor_mod, only: advance_forced_reactor_interval,forced_reactor_rhs,reactor_rhs, &
    reactor_specific_internal_energy,reactor_reduced_jacobian,forced_reactor_reduced_jacobian,backward_euler_trial
#ifdef PELEF_TEST_CVODE_SIGNED
  use sundials_forced_reactor_mod, only: advance_forced_reactor_cvode
#endif
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp) :: y(10),old(10),source(10),rate(10),t,restored_t,recovered,element_before(4),element_after(4),element_source(4)
  real(dp) :: atoms(4,10),mw(10),mass_error,residual,original_source(10)
  real(dp), parameter :: energy=469600.82041824982_dp,external_energy=-34759358186.641472_dp,dt=1e-8_dp
  integer :: steps,i
  logical :: ok
  call load_h2o2_full_thermo(species,ok);call require(ok,'thermo')
  call load_h2o2_full_mechanism(reactions,ok);call require(ok,'mechanism')
  call require(size(species)==10,'failed-input mechanism species identity')
  call require(all(species%name==[character(len=16)::'H2','H','O','O2','OH','H2O','HO2','H2O2','AR','N2']), &
    'failed-input species order')
  mw=species%molecular_weight
  atoms=0
  atoms(1,:)=[2._dp,1._dp,0._dp,0._dp,1._dp,2._dp,1._dp,2._dp,0._dp,0._dp]
  atoms(2,:)=[0._dp,0._dp,1._dp,2._dp,1._dp,1._dp,2._dp,2._dp,0._dp,0._dp]
  atoms(3,9)=1;atoms(4,10)=2
  source=[-8.8294917300898273e-27_dp,-1.4683862885727810e-20_dp,-2.8938113755454018e-18_dp, &
    -2.9736855926687663e-20_dp,-3.5592099773271806e-7_dp,899.85617208188751_dp, &
    -2.9820151528746120e-12_dp,2750.1162925888539_dp,0._dp,-7621.7142497948162_dp]
  original_source=source;old=0;old(10)=1;y=old;t=1000.0000000000001_dp
  element_before=matmul(atoms,old/mw);element_source=dt*matmul(atoms,source/mw)
  call advance_forced_reactor_interval(species,reactions,energy,source,external_energy,dt, &
    2e-7_dp,1e-12_dp,y,t,steps,ok)
  call require(ok.and.steps>0,'actual absent-species negative-forcing input accepts its physical endpoint')
  call require(all(y>=0).and.all(source==original_source),'accepted density is physical and external forcing is untouched')
  mass_error=abs(sum(y)-sum(old)-dt*sum(source))
  call require(mass_error<=5e-12_dp*sum(old),'forced mass balance')
  element_after=matmul(atoms,y/mw)
  residual=maxval(abs(element_after-element_before-element_source)/ &
    max(abs(element_before),abs(element_after),abs(element_source),1e-30_dp))
  call require(residual<=5e-12_dp,'H/O/Ar/N balance includes exact affine source')
  call forced_reactor_rhs(species,reactions,energy,dt,y,source,external_energy,t,rate,recovered,ok)
  call require(ok.and.abs(t-recovered)<=1e-11_dp*max(1._dp,abs(t)),'physical endpoint recovers ordinary EOS')
  print *, 'Signed trial endpoint:',steps,t,mass_error,residual
  do i=1,10
    print *, trim(species(i)%name),y(i)
  enddo
#ifdef PELEF_TEST_CVODE_SIGNED
  y=old;t=1000.0000000000001_dp
  call advance_forced_reactor_cvode(species,reactions,energy,source,external_energy,dt, &
    2e-7_dp,1e-12_dp,y,t,steps,ok)
  call require(ok.and.steps>0.and.all(y>=0),'CVODE signed trials retain physical endpoint')
  mass_error=abs(sum(y)-sum(old)-dt*sum(source))
  call require(mass_error<=5e-12_dp*sum(old),'CVODE forced mass balance')
  element_after=matmul(atoms,y/mw)
  residual=maxval(abs(element_after-element_before-element_source)/ &
    max(abs(element_before),abs(element_after),abs(element_source),1e-30_dp))
  call require(residual<=5e-12_dp,'CVODE signed evaluation retains H/O/Ar/N source balance')
  print *, 'CVODE signed trial endpoint:',steps,t,mass_error,residual
#endif
  ! Ordinary point API admission remains strict; a signed temporary flag must
  ! never make an invalid initial population acceptable implicitly.
  y=old;y(1)=-1e-15_dp;y(10)=y(10)+1e-15_dp
  call forced_reactor_rhs(species,reactions,energy,0._dp,y,source,external_energy,t,rate,recovered,ok)
  call require(.not.ok,'ordinary RHS rejects a negative initial species')
  y=old;y(1)=ieee_value(0._dp,ieee_quiet_nan);restored_t=t
  call reactor_rhs(species,reactions,1._dp,energy,y,t,rate,recovered,ok,permit_signed_trials=.true.)
  call require(.not.ok.and.all(rate==0),'signed point evaluator refuses NaN before a sign comparison')
  call advance_forced_reactor_interval(species,reactions,energy,source,external_energy,dt, &
    2e-7_dp,1e-12_dp,y,t,steps,ok)
  call require(.not.ok.and.t==restored_t,'nonfinite initial state refuses before publishing temperature')
  call check_negative_molar_pool()
  print *, 'PASS: signed forced Newton trials retain physical endpoint and elemental source balance'
  call check_runtime_eos_jacobians()
  call check_equilibrated_newton()
  call check_converged_signed_inverse()
contains
  subroutine check_converged_signed_inverse()
    type(nasa7_species) :: table(2)
    type(elementary_reaction) :: empty(0)
    real(dp) :: density(2),zero(2),rates(2),recovered,guess
    logical :: good
    table%name='constant';table%molecular_weight=universal_gas_constant
    table%temperature_min=300;table%temperature_mid=1000;table%temperature_max=3500
    table(1)%low_coefficients=[2._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp]
    table(1)%high_coefficients=table(1)%low_coefficients;table(2)=table(1)
    density=[-1e-14_dp,1._dp+1e-14_dp];zero=0;guess=2000._dp+5e-7_dp
    call configure_nasa7_eos_policy(table,.true.,good)
    call require(good,'default signed caloric policy')
    call forced_reactor_rhs(table,empty,2000._dp,0._dp,density,zero,0._dp,guess,rates,recovered,good, &
      permit_signed_trials=.true.)
    call require(good.and.recovered==guess,'default signed inverse retains sub-tolerance warm start')
    call configure_nasa7_eos_policy(table,.true.,good,converged_inverse=.true.)
    call require(good,'converged signed caloric policy')
    call forced_reactor_rhs(table,empty,2000._dp,0._dp,density,zero,0._dp,guess,rates,recovered,good, &
      permit_signed_trials=.true.)
    call require(good.and.abs(recovered-2000)<1e-12_dp.and.all(rates==0),'signed inverse analytic convergence')
  end subroutine check_converged_signed_inverse
  subroutine check_equilibrated_newton()
    type(nasa7_species), allocatable :: table(:)
    type(nasa7_species) :: identical_thermo(3)
    type(elementary_reaction) :: conversion(1)
    real(dp) :: initial(3),endpoint(3),energy,temperature,expected,residual
    integer :: iterations
    logical :: valid
    call load_h2o2_full_thermo(table,valid);call require(valid,'equilibrated Newton thermo')
    ! Synthetic isomers share thermodynamics and molecular weight. This
    ! isolates linear-system scaling from chemistry or EOS modelling error.
    identical_thermo=table(1)
    identical_thermo(1)%name='FAST_A';identical_thermo(2)%name='FAST_B';identical_thermo(3)%name='BULK'
    conversion(1)%reactant_stoich=[1._dp,0._dp,0._dp]
    conversion(1)%product_stoich=[0._dp,1._dp,0._dp]
    conversion(1)%forward_rate=arrhenius_rate(1e25_dp,0._dp,0._dp)
    conversion(1)%reversible=.false.
    initial=[1e-6_dp,0._dp,1._dp-1e-6_dp]
    energy=reactor_specific_internal_energy(identical_thermo,initial,1000._dp,valid)
    call require(valid,'synthetic isomer internal energy')
    call backward_euler_trial(identical_thermo,conversion,1._dp,energy,initial,1000._dp,1e-5_dp, &
      1e-8_dp,1e-12_dp,endpoint,temperature,iterations,valid,residual_reference=initial)
    call require(valid.and.iterations>0,'Newton resolves stiff independently scaled trace equations')
    expected=initial(1)/(1+1e20_dp)
    residual=abs(endpoint(1)-expected)/expected
    call require(residual<=1e-10_dp,'stiff trace agrees with independent analytic backward-Euler solution')
    call require(all(endpoint>=0).and.abs(sum(endpoint)-1)<=1e-14_dp,'stiff trace physical mass balance')
    call require(abs(endpoint(2)-(initial(1)-expected))<=1e-16_dp,'stiff trace converted material retained')
    call require(abs(endpoint(3)-initial(3))<=1e-14_dp,'inert bulk is retained')
    call require(abs(temperature-1000)<=1e-7_dp,'isomer conversion conserves energy')
    print *, 'PASS: equilibrated Newton analytic stiff trace:',iterations,residual
  end subroutine check_equilibrated_newton

  subroutine check_runtime_eos_jacobians()
    type(nasa7_species), allocatable :: table(:)
    type(elementary_reaction), allocatable :: mechanism(:)
    real(dp) :: y(10),plus(10),minus(10),jac(9,9),fp(10),fm(10),fd(9),zero(10)
    real(dp) :: wanted,anchor,energy,cv,cp_s,cv_s,h_s,e_s,s_s,recovered,tp,tm,residual
    real(dp), parameter :: step=2e-7_dp
    integer :: policy,thermal,signed,sp,cases
    logical :: ok
    call load_h2o2_full_thermo(table,ok);call require(ok,'pinned full thermo')
    call load_h2o2_full_mechanism(mechanism,ok);call require(ok,'pinned full mechanism')
    zero=0;cases=0
    do policy=1,2
      call configure_nasa7_eos_policy(table,.true.,ok)
      if(policy==2) call configure_nasa7_eos_policy(table,.true.,ok,600._dp,2400._dp)
      call require(ok,'explicit EOS policy')
      do thermal=1,2
        wanted=200;anchor=300
        if(thermal==2) then
          wanted=4000;anchor=3500
        endif
        if(policy==2) anchor=merge(600._dp,2400._dp,thermal==1)
        do signed=0,1
          y=1e-9_dp;y(1)=.028_dp;y(4)=.224_dp;y(10)=1-sum(y(:9))
          if(signed==1) then
            y(2)=-1e-15_dp;y(10)=1-sum(y(:9))
          endif
          energy=0;cv=0
          do sp=1,10
            call nasa7_mass_properties(table(sp),anchor,cp_s,cv_s,h_s,e_s,s_s,ok)
            call require(ok,'independent anchor thermodynamics')
            energy=energy+y(sp)*e_s;cv=cv+y(sp)*cv_s
          enddo
          energy=energy+(wanted-anchor)*cv
          if(signed==0) then
            call reactor_reduced_jacobian(table,mechanism,1._dp,energy,y,wanted,jac,recovered,ok)
          else
            call forced_reactor_reduced_jacobian(table,mechanism,energy,0._dp,y,zero,0._dp,wanted,jac,recovered,ok, &
              permit_signed_trials=.true.)
          endif
          call require(ok,'physical/signed fixed-energy Jacobian outside the fitted bounds')
          call require(abs(recovered-wanted)<=2e-10_dp*max(1._dp,wanted),'anchor inverse temperature')
          plus=y;minus=y;plus(1)=plus(1)+step;plus(10)=plus(10)-step
          minus(1)=minus(1)-step;minus(10)=minus(10)+step
          if(signed==0) then
            call reactor_rhs(table,mechanism,1._dp,energy,plus,wanted,fp,tp,ok)
            call require(ok,'positive finite-difference coordinate')
            call reactor_rhs(table,mechanism,1._dp,energy,minus,wanted,fm,tm,ok)
          else
            call forced_reactor_rhs(table,mechanism,energy,0._dp,plus,zero,0._dp,wanted,fp,tp,ok, &
              permit_signed_trials=.true.)
            call require(ok,'positive signed finite-difference coordinate')
            call forced_reactor_rhs(table,mechanism,energy,0._dp,minus,zero,0._dp,wanted,fm,tm,ok, &
              permit_signed_trials=.true.)
          endif
          call require(ok,'negative finite-difference coordinate')
          fd=(fp(:9)-fm(:9))/(2*step)
          residual=maxval(abs(jac(:,1)-fd))/max(maxval(abs(jac(:,1))),maxval(abs(fd)),1e-30_dp)
          print *, 'EOS thermal coordinate:',policy,thermal,signed,recovered,residual
          call require(residual<=2e-6_dp,'Jacobian agrees with independent fixed-energy RHS finite difference')
          cases=cases+1
        enddo
      enddo
    enddo
    call configure_nasa7_eos_policy(table,.false.,ok)
    call require(ok,'strict policy restored')
    y=1e-9_dp;y(1)=.028_dp;y(4)=.224_dp;y(10)=1-sum(y(:9))
    call reactor_reduced_jacobian(table,mechanism,1._dp,energy,y,4000._dp,jac,recovered,ok)
    call require(.not.ok.and.all(jac==0),'strict policy still refuses out-of-fit energy')
    print *, 'PASS: runtime EOS physical/signed Jacobian:',cases,'cold/hot default/explicit anchor cases; strict refusal'
  end subroutine check_runtime_eos_jacobians
  subroutine check_negative_molar_pool()
    type(nasa7_species), allocatable :: sp(:)
    type(elementary_reaction), allocatable :: rx(:),association(:)
    real(dp) :: initial(7),rho_y(7),affine(7),saved_affine(7),temperature,e0,saved_t
    real(dp) :: molar_before,molar_after,molar_source,initial_10(10),affine_10(10)
    real(dp), parameter :: captured_e=177815.7378503666550_dp,captured_source_e=1085066919.401288033_dp
    integer :: count,i
    logical :: good
    call load_h2o2_elementary_thermo(sp,good);call require(good,'pool thermodynamics')
    call load_h2o2_elementary_mechanism(rx,good);call require(good,'pool mechanism')
    do i=1,size(rx)
      call require(sum(rx(i)%product_stoich-rx(i)%reactant_stoich)==0,'captured pool conserves mole count')
    enddo
    initial=0;initial(7)=.3_dp;affine=0
    affine(6)=-.005982936912456951473_dp;affine(7)=672.2400835299458777_dp
    saved_affine=affine;rho_y=initial;temperature=1139.107472721268550_dp;saved_t=temperature
    call advance_forced_reactor_interval(sp,rx,captured_e,affine,captured_source_e,dt,1e-6_dp,3e-13_dp, &
      rho_y,temperature,count,good)
    call require(.not.good.and.count==0,'actual failed LF/Euler input rejects before adaptive substeps')
    call require(all(rho_y==initial).and.temperature==saved_t.and.all(affine==saved_affine), &
      'negative conserved-pool refusal preserves caller state and forcing')
    ! A dry pool deficit has no positive reservoir against which to hide a
    ! tracer. It remains impossible even far below the absolute solver tolerance.
    affine=0;affine(6)=-1e-20_dp;rho_y=initial;temperature=saved_t
    call advance_forced_reactor_interval(sp,rx,captured_e,affine,0._dp,dt,1e-6_dp,3e-13_dp, &
      rho_y,temperature,count,good)
    call require(.not.good.and.count==0.and.all(rho_y==initial),'dry trace pool deficit is not clipped or stepped')
    ! Negative affine water can be replenished by physical reacting H2/OH.
    ! The preflight must not reject a single affine-negative species by itself.
    initial=0;initial(1)=.05_dp;initial(5)=.05_dp;initial(7)=.2_dp
    e0=sum(initial)*reactor_specific_internal_energy(sp,initial/sum(initial),1000._dp,good)
    call require(good,'replenishing-pool initial energy');rho_y=initial;temperature=1000
    affine=0;affine(6)=-1e-4_dp
    molar_before=sum(initial(:6)/sp(:6)%molecular_weight)
    molar_source=dt*sum(affine(:6)/sp(:6)%molecular_weight)
    call advance_forced_reactor_interval(sp,rx,e0,affine,0._dp,dt,1e-6_dp,3e-13_dp, &
      rho_y,temperature,count,good)
    call require(good.and.count>0.and.rho_y(6)>0.and.all(rho_y>=0),'chemistry replenishes affine-negative water')
    molar_after=sum(rho_y(:6)/sp(:6)%molecular_weight)
    print *, 'Replenished pool residual:',count,molar_before,molar_after,molar_source, &
      molar_after-molar_before-molar_source
    ! Retain the original test's 5e-12 conserved-inventory gate. This
    ! replenishing case takes many adaptive substeps, rather than one sum.
    call require(abs(molar_after-molar_before-molar_source)<=5e-12_dp*molar_before, &
      'replenished endpoint preserves conserved pool and forcing')
    ! Real association changes mole count; a negative affine molar inventory
    ! cannot be used as a certificate for that component.
    allocate(association(1));association(1)=rx(1)
    association(1)%equation='2 H => H2';association(1)%reversible=.false.
    association(1)%reactant_stoich=0;association(1)%product_stoich=0
    association(1)%reactant_stoich(2)=2;association(1)%product_stoich(1)=1
    association(1)%forward_rate%pre_exponential=100
    association(1)%forward_rate%temperature_exponent=0;association(1)%forward_rate%activation_energy=0
    initial=0;initial(2)=.01_dp;initial(7)=.29_dp;affine=0;affine(1)=-1e-5_dp
    e0=sum(initial)*reactor_specific_internal_energy(sp,initial/sum(initial),1000._dp,good)
    call require(good,'association initial energy');rho_y=initial;temperature=1000
    call advance_forced_reactor_interval(sp,association,e0,affine,0._dp,dt,1e-6_dp,3e-13_dp, &
      rho_y,temperature,count,good)
    call require(good.and.count>0.and.rho_y(1)>0.and.all(rho_y>=0),'nonconserved molar component integrates normally')
    ! A pool cancellation within its arithmetic error bound is inconclusive.
    initial=0;initial(1)=.1_dp;initial(7)=.2_dp;affine=0
    affine(1)=-(initial(1)/dt)*(1+epsilon(1._dp))
    e0=sum(initial)*reactor_specific_internal_energy(sp,initial/sum(initial),1000._dp,good)
    call require(good,'roundoff-cancellation initial energy');rho_y=initial;temperature=1000
    call advance_forced_reactor_interval(sp,rx,e0,affine,0._dp,dt,1e-6_dp,3e-13_dp, &
      rho_y,temperature,count,good,maximum_steps=1)
    call require(count>0,'roundoff pool cancellation retains ordinary integration/final admission')
    ! A nonconserved active component does not erase a separate invariant for
    ! inert Ar, even when it participates in third-body rate coefficients.
    initial_10=old;affine_10=0;affine_10(9)=-1e-10_dp;y=initial_10;t=1000;restored_t=t
    call advance_forced_reactor_interval(species,reactions,energy,affine_10,0._dp,dt,2e-7_dp,1e-12_dp,y,t,count,good)
    call require(.not.good.and.count==0.and.all(y==initial_10).and.t==restored_t,'separate inert pool refuses atomically')
    print *, 'PASS: conserved molar-pool certificate, replenishment, association, cancellation and inert rollback'
  end subroutine check_negative_molar_pool

  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL:',label
      error stop 1
    endif
  end subroutine require
end program test_forced_signed_trials
