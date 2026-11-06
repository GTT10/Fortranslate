program test_dual_energy
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_is_nan
  use precision_mod, only: dp
  use dual_energy_mod, only: dual_energy_policy,reset_dual_energy_cell,reset_dual_energy_cells, &
    reset_nasa7_dual_energy_cell,reset_nasa7_dual_energy_patch,initialize_nasa7_internal_energy_patch
  use nasa7_thermo_mod, only: nasa7_species,universal_gas_constant
  implicit none
  type(nasa7_species) :: species(1)
  type(dual_energy_policy) :: policy
  real(dp) :: state(6),aux,delta,retained(6),cells(6,2),energies(2),temperatures(2),changes(2)
  real(dp) :: saved_cells(6,2),saved_energies(2),saved_changes(2)
  real(dp) :: patch(6,3,1,1),held(3,1,1),patch_temperature(3,1,1),fluid_volume(3,1,1),patch_change
  real(dp) :: formation=0
  logical :: ok,fail_eos=.false.
  call init(200._dp,20._dp)
  call reset_dual_energy_cell(state,aux,3._dp,policy,delta,ok)
  call check(ok.and.aux==100.and.delta==0,'default total-minus-kinetic reset without EOS')
  call init(200._dp,20._dp)
  call reset_dual_energy_cell(state,aux,3._dp,policy,delta,ok,potential=25._dp)
  call check(ok.and.aux==150.and.delta==0,'rotating potential included in internal energy')
  policy%allow_negative_energy=.false.
  call init(100._dp,20._dp)
  call reset_dual_energy_cell(state,aux,3._dp,policy,delta,ok)
  call check(ok.and.state(5)==120.and.aux==20.and.delta==20,'retained auxiliary repairs unresolved kinetic subtraction')
  policy%update_total_from_internal=.false.
  call init(100._dp,20._dp)
  call reset_dual_energy_cell(state,aux,3._dp,policy,delta,ok)
  call check(ok.and.state(5)==100.and.aux==20.and.delta==0,'disabled optional total repair')
  call init(0._dp,20._dp)
  call reset_dual_energy_cell(state,aux,3._dp,policy,delta,ok)
  call check(ok.and.state(5)==120.and.delta==120,'negative total branch repairs regardless of optional flag')
  policy%update_total_from_internal=.true.;policy%eta2=.1_dp
  call init(102._dp,20._dp)
  call reset_dual_energy_cell(state,aux,3._dp,policy,delta,ok)
  call check(ok.and.state(5)==120.and.aux==20,'eta threshold retains auxiliary')
  policy%eta2=.01_dp
  call init(102._dp,20._dp)
  call reset_dual_energy_cell(state,aux,3._dp,policy,delta,ok)
  call check(ok.and.state(5)==102.and.aux==2,'eta threshold accepts total-derived internal energy')
  call init(100._dp,-1._dp)
  call reset_dual_energy_cell(state,aux,3._dp,policy,delta,ok,evaluate_energy=eos)
  call check(ok.and.state(5)==101.and.aux>0.and.delta==1,'tiny-temperature replacement uses auxiliary difference')
  policy%allow_small_energy=.false.;policy%update_total_from_internal=.false.;formation=-100
  call init(-210._dp,-220._dp)
  call reset_dual_energy_cell(state,aux,3._dp,policy,delta,ok,evaluate_energy=eos)
  call check(ok.and.state(5)==-40.and.aux==-140.and.delta==170,'negative formation energy and thermal fallback')
  policy%eta2=.1_dp;formation=0
  call init(100._dp,20._dp)
  call reset_dual_energy_cell(state,aux,3._dp,policy,delta,ok,evaluate_energy=eos)
  call check(ok.and.state(5)==100.and.aux==20,'small-energy branch preserves reliable retained auxiliary')
  call init(100._dp,-1._dp);retained=state;fail_eos=.true.
  call reset_dual_energy_cell(state,aux,3._dp,policy,delta,ok,evaluate_energy=eos)
  call check(.not.ok.and.all(state==retained).and.aux==-1.and.delta==-77,'failed EOS preserves cell and correction')
  fail_eos=.false.;policy=dual_energy_policy()
  call init(200._dp,20._dp);cells(:,1)=state;cells(:,2)=state
  energies=20;temperatures=3;changes=-77
  cells(1,2)=ieee_value(0._dp,ieee_quiet_nan)
  saved_cells=cells;saved_cells(1,2)=0;saved_energies=energies;saved_changes=changes
  call reset_dual_energy_cells(cells,energies,temperatures,policy,changes,ok)
  call check(.not.ok.and.all(cells(:,1)==saved_cells(:,1)).and.all(energies==saved_energies).and. &
    all(changes==saved_changes),'late cell failure rolls back entire batch')
  cells(1,2)=2
  call reset_dual_energy_cells(cells,energies,temperatures,policy,changes,ok)
  call check(ok.and.all(energies==100).and.all(changes==0),'default batch reset')
  species(1)%name='X';species(1)%molecular_weight=10
  species(1)%temperature_min=200;species(1)%temperature_mid=1000;species(1)%temperature_max=5000
  species(1)%low_coefficients=0
  species(1)%low_coefficients(1)=1+100/universal_gas_constant
  species(1)%low_coefficients(6)=-1000/universal_gas_constant
  species(1)%high_coefficients=species(1)%low_coefficients
  policy%allow_small_energy=.false.;policy%update_total_from_internal=.false.
  call init(-210._dp,-220._dp)
  call reset_nasa7_dual_energy_cell(species,state,aux,3._dp,policy,delta,ok)
  call check(ok.and.abs(state(5)+40)<1e-12_dp.and.abs(aux+140)<1e-12_dp, &
    'NASA7 reset evaluates formation energy and explicit transient temperature')
  retained=state;species(1)%molecular_weight=0;delta=-77
  call reset_nasa7_dual_energy_cell(species,state,aux,3._dp,policy,delta,ok)
  call check(.not.ok.and.all(state==retained).and.delta==-77,'invalid NASA7 rollback')
  species(1)%molecular_weight=10;policy=dual_energy_policy();policy%allow_negative_energy=.false.
  call init(100._dp,20._dp)
  patch(:,1,1,1)=ieee_value(0._dp,ieee_quiet_nan)
  patch(:,2,1,1)=state;patch(:,3,1,1)=state
  patch(1,3,1,1)=ieee_value(0._dp,ieee_quiet_nan)
  held(:,1,1)=[ieee_value(0._dp,ieee_quiet_nan),20._dp,20._dp]
  patch_temperature(:,1,1)=[ieee_value(0._dp,ieee_quiet_nan),3._dp,3._dp]
  fluid_volume(:,1,1)=[0._dp,.25_dp,1._dp];patch_change=7
  call reset_nasa7_dual_energy_patch(species,patch,held,patch_temperature,fluid_volume,policy,patch_change,ok)
  call check(.not.ok.and.patch(5,2,1,1)==100.and.held(2,1,1)==20.and.patch_change==7, &
    'late fluid failure rolls back EB patch')
  patch(1,3,1,1)=2
  call reset_nasa7_dual_energy_patch(species,patch,held,patch_temperature,fluid_volume,policy,patch_change,ok)
  call check(ok.and.patch(5,2,1,1)==120.and.patch(5,3,1,1)==120.and.patch_change==32.and. &
    ieee_is_nan(patch(5,1,1,1)).and.ieee_is_nan(held(1,1,1)), &
    'EB patch reset preserves covered sentinels and weights total-energy repair')
  call check_native_initial_energy()
  print *, 'PASS: dual-energy policy branches, frame potential and transactional failure'
contains
  subroutine check_native_initial_energy()
    real(dp) :: gas(6,3,1,1),temps(3,1,1),mask(3,1,1),stored(3,1,1)
    gas=0;gas(:,1,1,1)=ieee_value(0._dp,ieee_quiet_nan)
    gas(:,2,1,1)=[2._dp,2.e10_dp,0._dp,0._dp,1.e20_dp,2._dp]
    gas(:,3,1,1)=gas(:,2,1,1);gas(5,3,1,1)=-1.e20_dp
    temps(:,1,1)=[ieee_value(0._dp,ieee_quiet_nan),300._dp,5001._dp]
    mask(:,1,1)=[0._dp,.25_dp,1._dp];stored(:,1,1)=[ieee_value(0._dp,ieee_quiet_nan),17._dp,19._dp]
    call initialize_nasa7_internal_energy_patch(species,gas,temps,mask,stored,ok)
    call check(.not.ok.and.stored(2,1,1)==17.and.stored(3,1,1)==19.and.ieee_is_nan(stored(1,1,1)), &
      'native initialization late EOS failure preserves all held fields')
    temps(3,1,1)=300
    call initialize_nasa7_internal_energy_patch(species,gas,temps,mask,stored,ok)
    call check(ok.and.all(abs(stored(2:3,1,1)-5800._dp)<1.e-10_dp).and.ieee_is_nan(stored(1,1,1)), &
      'native initialization uses analytic rho*e(T,Y) despite unresolved or modified total energy')
    call check(gas(5,2,1,1)-.5_dp*gas(2,2,1,1)**2/gas(1,2,1,1)==0, &
      'initialization fixture loses internal energy in kinetic subtraction')
    gas(6,3,1,1)=1;stored(2:3,1,1)=[23._dp,29._dp]
    call initialize_nasa7_internal_energy_patch(species,gas,temps,mask,stored,ok)
    call check(.not.ok.and.all(stored(2:3,1,1)==[23._dp,29._dp]), &
      'native initialization rejects inconsistent species mass atomically')
    gas(:,3,1,1)=0;gas(1,3,1,1)=huge(1._dp)/4;gas(6,3,1,1)=gas(1,3,1,1)
    call initialize_nasa7_internal_energy_patch(species,gas,temps,mask,stored,ok)
    call check(.not.ok.and.all(stored(2:3,1,1)==[23._dp,29._dp]), &
      'native initialization refuses density-energy overflow before traps or partial publication')
    gas(1,3,1,1)=tiny(1._dp);gas(6,3,1,1)=huge(1._dp)
    call initialize_nasa7_internal_energy_patch(species,gas,temps,mask,stored,ok)
    call check(.not.ok.and.all(stored(2:3,1,1)==[23._dp,29._dp]), &
      'native initialization refuses composition-ratio overflow before traps or partial publication')
  end subroutine check_native_initial_energy
  subroutine init(total,internal)
    real(dp), intent(in) :: total,internal
    state=[2._dp,20._dp,0._dp,0._dp,total,2._dp];aux=internal;delta=-77
  end subroutine
  subroutine eos(t,y,e,valid)
    real(dp), intent(in) :: t,y(:)
    real(dp), intent(out) :: e
    logical, intent(out) :: valid
    e=10*t+formation;valid=.not.fail_eos.and.size(y)==1.and.all(y==1)
  end subroutine
  subroutine check(valid,message)
    logical, intent(in) :: valid
    character(len=*), intent(in) :: message
    if(valid) return
    print *, 'FAIL: ',message
    error stop 1
  end subroutine
end program test_dual_energy
