program test_regular_gas_checkpoint
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved,reactive_conserved_to_primitive
  use reactive_csv_io_3d_mod, only: write_reactive_3d_csv
  use reactive_3d_checkpoint_mod, only: write_regular_gas_checkpoint,read_regular_gas_checkpoint
  implicit none
  type(nasa7_species), allocatable :: species(:)
  real(dp) :: state(12,1,1,1),temp(1,1,1),ir(8,1,1,1),heat(1,1,1),q(12),sound,t
  real(dp) :: seed(12,1,1,1),seed_t(1,1,1),initial(12),wall(12),source(12),metrics(2),time,previous_timestep
  integer :: steps,unit,other,status
  logical :: ok
  character(len=1024) :: message,line
  call load_h2o2_elementary_thermo(species,ok);call require(ok,'thermo')
  q=0;q(1)=.4_dp;q(5)=101325;q(12)=1
  call reactive_primitive_to_conserved(species,q,state(:,1,1,1),t,sound,ok);call require(ok,'seed')
  temp=t;seed=state;seed_t=temp;ir=.125_dp;heat=37;initial=1;wall=2;source=3;metrics=[4._dp,.5_dp];time=1;steps=2
  call write_regular_gas_checkpoint('regular-checkpoint-test.chk',identity,species,state,temp,time,steps, &
    initial,wall,source,metrics,ok,message,ir,heat,previous_timestep=.25_dp)
  call require(ok,'atomic checkpoint write')
  state(6,1,1,1)=-1
  call write_regular_gas_checkpoint('regular-checkpoint-test.chk',identity,species,state,temp,time,steps, &
    initial,wall,source,metrics,ok,message,ir,heat)
  call require(.not.ok,'failed writer retains previous file')
  state=0;temp=0;ir=0;heat=0;initial=0;wall=0;source=0;metrics=0;time=0;steps=0;previous_timestep=17
  call read_regular_gas_checkpoint('regular-checkpoint-test.chk',identity,species,state,temp,time,steps, &
    initial,wall,source,metrics,ok,message,ir,heat,previous_timestep=previous_timestep)
  call require(ok.and.all(state==seed).and.all(temp==seed_t),'state exact restore after rejected overwrite')
  call require(previous_timestep==.25_dp,'accepted timestep exact restore')
  call require(all(ir==.125_dp).and.all(heat==37),'IR heat exact restore')
  call require(time==1.and.steps==2.and.all(initial==1).and.all(wall==2).and.all(source==3).and. &
    all(metrics==[4._dp,.5_dp]),'progress inventories and diagnostics restore')
  open(newunit=unit,file='regular-checkpoint-test.chk',status='old',action='read')
  open(newunit=other,file='regular-checkpoint-truncated.chk',status='replace',action='write')
  do
    read(unit,'(a)',iostat=status) line
    if(status/=0) exit
    if(trim(line)=='END_PELEF_REGULAR_GAS') exit
    write(other,'(a)') trim(line)
  enddo
  close(unit);close(other)
  state=17;temp=19;ir=23;heat=29;time=31;steps=7;initial=37;wall=41;source=43;metrics=47;previous_timestep=53
  call read_regular_gas_checkpoint('regular-checkpoint-truncated.chk',identity,species,state,temp,time,steps, &
    initial,wall,source,metrics,ok,message,ir,heat,previous_timestep=previous_timestep)
  call require(.not.ok.and.all(state==17).and.all(temp==19).and.all(ir==23).and.all(heat==29), &
    'late reader failure preserves all state')
  call require(time==31.and.steps==7.and.all(initial==37).and.all(wall==41).and.all(source==43).and.all(metrics==47), &
    'late reader failure preserves progress inventories and diagnostics')
  call require(previous_timestep==53,'late reader failure preserves accepted interval')
  state=seed;temp=seed_t;time=1;steps=2;metrics=[4._dp,.5_dp]
  call write_regular_gas_checkpoint('regular-checkpoint-test.chk',identity,species,state,temp,time,steps, &
    initial,wall,source,metrics,ok,message,ir,heat,previous_timestep=-1._dp)
  call require(.not.ok,'invalid accepted interval refuses writer')
  call write_regular_gas_checkpoint('regular-checkpoint-legacy.chk',identity,species,state,temp,time,steps, &
    initial,wall,source,metrics,ok,message,ir,heat)
  call require(ok,'legacy kernel snapshot write')
  call read_regular_gas_checkpoint('regular-checkpoint-legacy.chk',identity,species,state,temp,time,steps, &
    initial,wall,source,metrics,ok,message,ir,heat,previous_timestep=previous_timestep)
  call require(.not.ok.and.previous_timestep==53,'required history cannot use legacy snapshot')
  call read_regular_gas_checkpoint('regular-checkpoint-legacy.chk',identity,species,state,temp,time,steps, &
    initial,wall,source,metrics,ok,message,ir,heat)
  call require(ok,'explicit legacy kernel read')
  call check_auxiliary_checkpoint()
  call check_eos_identity()
  print *, 'PASS: regular gas checkpoint exact restore, rejected overwrite and transactional late read failure'
contains
  subroutine check_eos_identity()
    state=seed;temp=seed_t;time=1;steps=2;metrics=[4._dp,.5_dp]
    species%pinned_eos_extrapolation=.true.
    call read_regular_gas_checkpoint('regular-checkpoint-legacy.chk',identity,species,state,temp,time,steps, &
      initial,wall,source,metrics,ok,message,ir,heat)
    call require(.not.ok.and.all(state==seed).and.all(temp==seed_t),'strict to reference EOS rollback')
    call write_regular_gas_checkpoint('regular-checkpoint-eos.chk',identity,species,state,temp,time,steps, &
      initial,wall,source,metrics,ok,message,ir,heat)
    call require(ok,'reference EOS checkpoint write')
    call read_regular_gas_checkpoint('regular-checkpoint-eos.chk',identity,species,state,temp,time,steps, &
      initial,wall,source,metrics,ok,message,ir,heat)
    call require(ok,'reference EOS same-mode read')
    species%pinned_eos_convergence=.true.
    call read_regular_gas_checkpoint('regular-checkpoint-eos.chk',identity,species,state,temp,time,steps, &
      initial,wall,source,metrics,ok,message,ir,heat)
    call require(.not.ok.and.all(state==seed).and.all(temp==seed_t),'changed caloric precision rollback')
    call write_regular_gas_checkpoint('regular-checkpoint-converged.chk',identity,species,state,temp,time,steps, &
      initial,wall,source,metrics,ok,message,ir,heat)
    call require(ok,'converged EOS checkpoint write')
    call read_regular_gas_checkpoint('regular-checkpoint-converged.chk',identity,species,state,temp,time,steps, &
      initial,wall,source,metrics,ok,message,ir,heat)
    call require(ok.and.all(state==seed),'converged EOS same-mode read')
    species%pinned_eos_convergence=.false.
    call read_regular_gas_checkpoint('regular-checkpoint-converged.chk',identity,species,state,temp,time,steps, &
      initial,wall,source,metrics,ok,message,ir,heat)
    call require(.not.ok.and.all(state==seed),'converged to default precision rollback')
    species%pinned_eos_temperature_min=300;species%pinned_eos_temperature_max=4000
    call read_regular_gas_checkpoint('regular-checkpoint-eos.chk',identity,species,state,temp,time,steps, &
      initial,wall,source,metrics,ok,message,ir,heat)
    call require(.not.ok.and.all(state==seed).and.all(temp==seed_t),'changed EOS bounds rollback')
    species%pinned_eos_temperature_min=0;species%pinned_eos_temperature_max=0
    species%pinned_eos_extrapolation=.false.
    call read_regular_gas_checkpoint('regular-checkpoint-eos.chk',identity,species,state,temp,time,steps, &
      initial,wall,source,metrics,ok,message,ir,heat)
    call require(.not.ok.and.all(state==seed).and.all(temp==seed_t),'reference to strict EOS rollback')
  end subroutine check_eos_identity
  subroutine check_auxiliary_checkpoint()
    real(dp) :: u(12,1,1,1),temperature(1,1,1),auxiliary(1,1,1),kept(12,1,1,1), &
      kept_auxiliary(1,1,1),totals(12),boundaries(12),forces(12),limits(2),stamp,interval,correction
    real(dp) :: primitive(12),recovered(12),cs,tt,csv_row(28)
    integer :: count,infile,outfile,ios
    logical :: good
    primitive=q;primitive(2)=1.e8_dp
    call reactive_primitive_to_conserved(species,q,u(:,1,1,1),tt,cs,good)
    call require(good,'auxiliary internal-energy seed')
    auxiliary=u(5,1,1,1);temperature=tt
    call reactive_primitive_to_conserved(species,primitive,u(:,1,1,1),tt,cs,good)
    call require(good,'high kinetic-energy seed')
    ! The retained scalar must avoid roundoff in E - KE at high Mach number.
    call reactive_conserved_to_primitive(species,u(:,1,1,1),temperature(1,1,1),recovered,tt,cs,good, &
      internal_energy_density=auxiliary(1,1,1))
    call require(good.and.abs(tt-temperature(1,1,1))<1e-10_dp,'direct auxiliary EOS at high kinetic energy')
    ! A differing total energy is legitimate when the reset policy retains rho e.
    u=seed;u(5,1,1,1)=u(5,1,1,1)+333._dp
    totals=seed(:,1,1,1);boundaries=0;forces=0;limits=[0._dp,1._dp]
    correction=333;stamp=1;count=2;interval=.25_dp
    kept=u;kept_auxiliary=auxiliary
    call write_regular_gas_checkpoint('regular-checkpoint-auxiliary.chk',identity,species,u,temperature,stamp,count, &
      totals,boundaries,forces,limits,good,message,previous_timestep=interval, &
      internal_energy_density=auxiliary,energy_correction=correction)
    call require(good,'V3 auxiliary and correction write')
    u=0;temperature=0;auxiliary=0;correction=0
    call read_regular_gas_checkpoint('regular-checkpoint-auxiliary.chk',identity,species,u,temperature,stamp,count, &
      totals,boundaries,forces,limits,good,message,spacing=[1._dp,1._dp,1._dp],previous_timestep=interval, &
      internal_energy_density=auxiliary,energy_correction=correction)
    call require(good.and.all(u==kept).and.all(auxiliary==kept_auxiliary).and.correction==333, &
      'V3 exact scalar restore and energy-correction conservation')
    call require(all(temperature==seed_t),'V3 temperature uses auxiliary instead of differing total')
    call write_reactive_3d_csv('regular-checkpoint-auxiliary.csv',species,[0._dp],[0._dp],[0._dp], &
      u,temperature,1,1,1,stamp,good,message,internal_energy_density=auxiliary)
    call require(good,'auxiliary CSV output')
    open(newunit=infile,file='regular-checkpoint-auxiliary.csv',status='old',action='read')
    read(infile,'(a)') line
    read(infile,'(a)') line
    close(infile)
    read(line,*,iostat=ios) csv_row
    call require(ios==0.and.abs(csv_row(9)-q(5))<1e-7_dp.and.abs(csv_row(10)-seed_t(1,1,1))<1e-10_dp, &
      'CSV pressure and temperature use auxiliary thermodynamic state')
    call read_regular_gas_checkpoint('regular-checkpoint-auxiliary.chk',identity,species,u,temperature,stamp,count, &
      totals,boundaries,forces,limits,good,message,previous_timestep=interval)
    call require(.not.good,'V3 cannot silently discard auxiliary state')
    call read_regular_gas_checkpoint('regular-checkpoint-test.chk',identity,species,u,temperature,stamp,count, &
      totals,boundaries,forces,limits,good,message,previous_timestep=interval, &
      internal_energy_density=auxiliary,energy_correction=correction)
    call require(.not.good,'V2 cannot silently manufacture auxiliary state')
    call write_regular_gas_checkpoint('regular-checkpoint-auxiliary-nohistory.chk',identity,species,u,temperature,stamp,count, &
      totals,boundaries,forces,limits,good,message,internal_energy_density=auxiliary,energy_correction=correction)
    call require(good,'V3 optional timestep write')
    call read_regular_gas_checkpoint('regular-checkpoint-auxiliary-nohistory.chk',identity,species,u,temperature,stamp,count, &
      totals,boundaries,forces,limits,good,message,internal_energy_density=auxiliary,energy_correction=correction)
    call require(good,'V3 optional timestep read')
    call read_regular_gas_checkpoint('regular-checkpoint-auxiliary-nohistory.chk',identity,species,u,temperature,stamp,count, &
      totals,boundaries,forces,limits,good,message,previous_timestep=interval, &
      internal_energy_density=auxiliary,energy_correction=correction)
    call require(.not.good,'V3 required history cannot use omitted history')
    open(newunit=infile,file='regular-checkpoint-auxiliary.chk',status='old',action='read')
    open(newunit=outfile,file='regular-checkpoint-auxiliary-truncated.chk',status='replace',action='write')
    do
      read(infile,'(a)',iostat=ios) line
      if(ios/=0.or.trim(line)=='END_PELEF_REGULAR_GAS') exit
      write(outfile,'(a)') trim(line)
    enddo
    close(infile);close(outfile)
    u=17;temperature=19;auxiliary=23;correction=29;interval=31
    call read_regular_gas_checkpoint('regular-checkpoint-auxiliary-truncated.chk',identity,species,u,temperature,stamp,count, &
      totals,boundaries,forces,limits,good,message,previous_timestep=interval, &
      internal_energy_density=auxiliary,energy_correction=correction)
    call require(.not.good.and.all(u==17).and.all(temperature==19).and.all(auxiliary==23).and.correction==29.and. &
      interval==31,'V3 late failure keeps scalar correction and other fields')
  end subroutine check_auxiliary_checkpoint
  subroutine identity(unit,status)
    integer, intent(in) :: unit
    integer, intent(out) :: status
    write(unit,'(a)',iostat=status) 'TEST_REGULAR_CHECKPOINT_V1'
  end subroutine
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if(.not.condition) then
      print *, 'FAIL: ',message
      error stop 1
    endif
  end subroutine
end program
