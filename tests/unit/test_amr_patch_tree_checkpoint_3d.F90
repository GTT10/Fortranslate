program test_amr_patch_tree_checkpoint_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use simulation_config_reactive_3d_mod, only: reactive_3d_config,reactive_3d_mole_fractions
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: initialize_amr_patch_tree_3d,prolong_amr_patch_tree_3d, &
    average_down_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
  use amr_patch_tree_checkpoint_3d_mod, only: amr_tree_checkpoint_3d,build_amr_tree_checkpoint_context_3d, &
    write_amr_tree_checkpoint_3d,read_amr_tree_checkpoint_3d,read_amr_tree_checkpoint_header_3d
  implicit none
  type(amr_tree_checkpoint_3d) :: source,loaded,saved
  type(reactive_3d_config) :: config,changed
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_3d) :: patches(2)
  character(len=:), allocatable :: context,other
  real(dp) :: x(7),u(12,4,4,4),t(1,4,4,4),q(12),sound
  integer :: c,unit,input,status,mode,node
  character(len=65536) :: line
  logical :: ok,exists
  character(len=16) :: test_mode
  call get_command_argument(1,test_mode)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_mechanism(reactions,ok)
  call require(ok,'reactions')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  config%nx=4;config%ny=4;config%nz=4;config%thermo_model='elementary'
  call reactive_3d_mole_fractions(config,7,x,ok)
  call require(ok,'composition')
  call build_context(config,context)
  changed=config;changed%difmag=0
  call build_context(changed,other)
  call require(context/=other,'automatic and zero non-PPM coefficient identity differ')
  changed%difmag=.1_dp
  call build_context(changed,other)
  call require(context==other,'automatic and explicit reference coefficient identity')
  changed%difmag=.037_dp
  call build_context(changed,other)
  call require(context/=other,'positive non-PPM coefficient identity')

  changed=config;changed%final_time=2e-6_dp;changed%maximum_steps=2000;changed%output_file='elsewhere.csv'
  call build_context(changed,other)
  call require(context==other,'continuation controls excluded')
  changed=config;changed%cfl=.1_dp
  call build_context(changed,other)
  call require(context/=other,'flow controls captured')
  q=0;q(1)=.3_dp;q(5)=101325;q(6)=.03_dp;q(9)=.22_dp;q(12)=.75_dp
  call reactive_primitive_to_conserved(species,q,u(:,1,1,1),t(1,1,1,1),sound,ok)
  call require(ok,'primitive')
  do c=1,12
    u(c,:,:,:)=u(c,1,1,1)
  end do
  t=t(1,1,1,1)
  call initialize_amr_patch_3d(4,4,4,1,2,1,2,1,2,2,patches(1),ok)
  call initialize_amr_patch_3d(4,4,4,2,3,1,2,1,2,3,patches(2),ok)
  call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,source%tree,ok)
  call require(ok,'three-level topology')
  call prolong_amr_patch_tree_3d(source%tree,u,source%state,ok)
  call require(ok,'state')
  call prolong_amr_patch_tree_3d(source%tree,t,source%temperature,ok)
  call require(ok,'temperature')
  call average_down_amr_patch_tree_3d(source%tree,source%state,ok)
  call require(ok,'synchronization')
  allocate(source%initial(12),source%boundary(12))
  call composite_integrals_amr_patch_tree_3d(source%tree,source%state,[.001_dp,.001_dp,.001_dp],source%initial,ok)
  call require(ok,'inventory')
  source%boundary=.001_dp*source%initial
  source%time=1e-9_dp;source%steps=3;source%regrids=2
  source%maximum_reflux=.125_dp;source%minimum_theta=.75_dp;source%maximum_remap_error=1e-16_dp
  if(test_mode=='publication') then
    open(newunit=unit,file='publication.chk',status='replace');write(unit,'(a)') 'preserve';close(unit)
    open(newunit=unit,file='publication.chk.tmp',status='new');write(unit,'(a)') 'temporary';close(unit)
    call write_amr_tree_checkpoint_3d('publication.chk',species,context,source,ok)
    call require(.not.ok,'temporary collision rejected')
    open(newunit=unit,file='publication.chk',status='old');read(unit,'(a)') line;close(unit)
    call require(trim(line)=='preserve','collision preserves final')
    inquire(file='publication.chk',size=status)
    call require(status==9,'collision preserves final size')
    open(newunit=unit,file='publication.chk.tmp',status='old');read(unit,'(a)') line
    call require(trim(line)=='temporary','collision preserves other temporary')
    inquire(unit=unit,size=status)
    call require(status==10,'collision preserves other temporary size')
    close(unit,status='delete')
    call execute_command_line('mkdir publication-directory',exitstat=status)
    call require(status==0,'create rename-rejection fixture')
    call write_amr_tree_checkpoint_3d('publication-directory',species,context,source,ok)
    call require(.not.ok,'rename failure rejected')
    inquire(file='publication-directory.tmp',exist=exists)
    call require(.not.exists,'rename failure removes own temporary')
    call execute_command_line('rmdir publication-directory',exitstat=status)
    call require(status==0,'rename failure preserves directory')
    print *, 'PASS: serial publication collision and rename rollback'
    stop
  endif
  call write_amr_tree_checkpoint_3d('tree.chk',species,context,source,ok)
  call require(ok,'write snapshot')
  loaded=source;loaded%state%fields(1)%values=-1;loaded%time=0
  call read_amr_tree_checkpoint_3d('tree.chk',species,context,[4,4,4],loaded,ok)
  call require(ok,'read snapshot')
  call same(loaded,source)
  source%previous_timestep=2e-10_dp;source%post_regrid=.true.
  call write_amr_tree_checkpoint_3d('tree_history.chk',species,context,source,ok,indexed=.true.)
  call require(ok,'indexed history write')
  call read_amr_tree_checkpoint_3d('tree_history.chk',species,context,[4,4,4],loaded,ok)
  call require(ok.and.loaded%previous_timestep==2e-10_dp.and.loaded%post_regrid,'history exact read')
  loaded%previous_timestep=53
  call read_amr_tree_checkpoint_3d('tree_history.chk',species,other,[4,4,4],loaded,ok)
  call require(.not.ok.and.loaded%previous_timestep==53.and.loaded%post_regrid,'history rollback on context failure')
  source%previous_timestep=-1
  call write_amr_tree_checkpoint_3d('tree_history_invalid.chk',species,context,source,ok)
  call require(.not.ok,'invalid history rejects writer')
  deallocate(source%previous_timestep);source%post_regrid=.false.;loaded=source
  saved=loaded
  call read_amr_tree_checkpoint_3d('tree.chk',species,other,[4,4,4],loaded,ok)
  call require(.not.ok,'context mismatch rejection')
  call same(loaded,saved)
  transport(1)%diameter=transport(1)%diameter*1.01_dp
  call build_context(config,other)
  call require(other/=context,'transport provenance captured')
  call read_amr_tree_checkpoint_3d('tree.chk',species,other,[4,4,4],loaded,ok)
  call require(.not.ok,'transport mismatch rejection')
  call same(loaded,saved)
  source%temperature%fields(3)%values=1
  call write_amr_tree_checkpoint_3d('tree.chk',species,context,source,ok)
  call require(.not.ok,'invalid write rejected before replacement')
  call read_amr_tree_checkpoint_3d('tree.chk',species,context,[4,4,4],loaded,ok)
  call require(ok,'existing checkpoint survived invalid write')
  call same(loaded,saved)
  open(newunit=input,file='tree.chk',status='old')
  open(newunit=unit,file='truncated.chk',status='replace')
  do
    read(input,'(a)',iostat=status) line
    call require(status==0,'read complete source checkpoint')
    if (trim(line)=='END_AMR_TREE_CHECKPOINT') exit
    write(unit,'(a)') trim(line)
  end do
  close(input);close(unit)
  call read_amr_tree_checkpoint_3d('truncated.chk',species,context,[4,4,4],loaded,ok)
  call require(.not.ok,'truncated read rejection')
  call same(loaded,saved)
  source=saved
  source%source=-.002_dp*source%initial
  call prolong_amr_patch_tree_3d(source%tree,-u(:8,:,:,:),source%reaction_source,ok)
  call require(ok,'signed reaction estimates')
  do mode=0,3
    if(mode==2) then
      allocate(source%heat_release%fields(size(source%tree%nodes)))
      do node=1,size(source%tree%nodes)
        source%heat_release%fields(node)%values=source%temperature%fields(node)%values
        source%heat_release%fields(node)%values=real((-1)**node*node,dp)*1.e9_dp
      end do
    end if
    call write_amr_tree_checkpoint_3d('sdc.chk',species,context,source,ok,indexed=mod(mode,2)==1)
    call require(ok,'SDC write')
    open(newunit=input,file='sdc.chk',status='old')
    call read_amr_tree_checkpoint_header_3d(input,species,context,[4,4,4],loaded,ok)
    close(input)
    call require(.not.ok,'legacy header callers reject SDC payload')
    call read_amr_tree_checkpoint_3d('sdc.chk',species,context,[4,4,4],loaded,ok)
    call require(ok,'SDC read')
    call same(loaded,source)
    call read_amr_tree_checkpoint_3d('sdc.chk',species,other,[4,4,4],loaded,ok)
    call require(.not.ok,'SDC context rejection')
    call same(loaded,source)
  end do
  deallocate(source%reaction_source%fields)
  call write_amr_tree_checkpoint_3d('sdc.chk',species,context,source,ok)
  call require(ok,'nonreacting SDC write')
  call read_amr_tree_checkpoint_3d('sdc.chk',species,context,[4,4,4],loaded,ok)
  call require(ok,'nonreacting SDC read')
  call same(loaded,source)
  call check_internal_energy()
  call check_mol_covered_reset()
contains
  subroutine check_mol_covered_reset()
    type(amr_tree_checkpoint_3d) :: payload,actual,prior
    real(dp), allocatable :: recovered(:,:,:)
    integer :: peer,n(3),i,variant
    payload=source
    payload%internal_energy=payload%temperature;payload%energy_correction=0
    peer=size(payload%tree%nodes)
    do i=1,payload%tree%nodes(peer)%shape(1)
      payload%state%fields(peer)%values(2,i,:,:)=real(i,dp)*payload%state%fields(peer)%values(1,i,:,:)
    enddo
    call average_down_amr_patch_tree_3d(payload%tree,payload%state,ok)
    call require(ok,'MOL covered conserved restriction')
    do peer=1,size(payload%tree%nodes)
      n=payload%tree%nodes(peer)%shape
      payload%internal_energy%fields(peer)%values(1,:,:,:)=payload%state%fields(peer)%values(5,:,:,:) - &
        .5_dp*sum(payload%state%fields(peer)%values(2:4,:,:,:)**2,dim=1)/payload%state%fields(peer)%values(1,:,:,:)
      allocate(recovered(n(1),n(2),n(3)))
      call recover_reactive_temperatures_3d(species,payload%state%fields(peer)%values, &
        payload%temperature%fields(peer)%values(1,:,:,:),n(1),n(2),n(3),recovered,ok, &
        internal_energy_density=payload%internal_energy%fields(peer)%values(1,:,:,:))
      call require(ok,'MOL reset temperature');payload%temperature%fields(peer)%values(1,:,:,:)=recovered
      deallocate(recovered)
    enddo
    call write_amr_tree_checkpoint_3d('mol_reset.chk',species,context,payload,ok)
    call require(.not.ok,'covered reset requires explicit MOL mode')
    do variant=0,1
      call write_amr_tree_checkpoint_3d('mol_reset.chk',species,context,payload,ok, &
        indexed=variant==1,mol_energy_reset=.true.)
      call require(ok,'MOL reset checkpoint write')
      actual=source;prior=actual
      call read_amr_tree_checkpoint_3d('mol_reset.chk',species,context,[4,4,4],actual,ok,allow_internal_energy=.true.)
      call require(.not.ok,'MOL reset strict reader rejects');call same(actual,prior)
      call read_amr_tree_checkpoint_3d('mol_reset.chk',species,context,[4,4,4],actual,ok, &
        allow_internal_energy=.true.,mol_energy_reset=.true.)
      call require(ok,'MOL reset checkpoint read');call same(actual,payload)
      do peer=1,size(payload%tree%nodes)
        call require(all(actual%internal_energy%fields(peer)%values==payload%internal_energy%fields(peer)%values), &
          'MOL reset retained values exact')
      enddo
    enddo
    prior=payload
    payload%internal_energy%fields(1)%values(1,1,1,1)=payload%internal_energy%fields(1)%values(1,1,1,1)+1
    call write_amr_tree_checkpoint_3d('mol_reset.chk',species,context,payload,ok,mol_energy_reset=.true.)
    call require(.not.ok,'MOL mode rejects arbitrary covered held change')
    payload=prior;payload%state%fields(1)%values(5,1,1,1)=payload%state%fields(1)%values(5,1,1,1)+1
    call write_amr_tree_checkpoint_3d('mol_reset.chk',species,context,payload,ok,mol_energy_reset=.true.)
    call require(.not.ok,'MOL mode preserves conserved restriction')
    call read_amr_tree_checkpoint_3d('mol_reset.chk',species,context,[4,4,4],actual,ok, &
      allow_internal_energy=.true.,mol_energy_reset=.true.)
    call require(ok,'MOL valid file survives rejected writes');call same(actual,prior)
  end subroutine check_mol_covered_reset

  subroutine check_internal_energy()
    type(amr_tree_checkpoint_3d) :: payload,actual,prior
    real(dp), allocatable :: recovered(:,:,:)
    integer :: variant,peer,n(3),from,to,stat
    character(len=65536) :: record
    do variant=0,7
      payload=source
      if(variant>=4) payload%pressure_work_change=-.25_dp
      if(variant<2) then
        deallocate(payload%source,payload%heat_release%fields)
      endif
      payload%previous_timestep=1e-10_dp
      payload%internal_energy=payload%temperature;payload%energy_correction=-7._dp
      do peer=1,size(payload%tree%nodes)
        n=payload%tree%nodes(peer)%shape
        payload%internal_energy%fields(peer)%values(1,:,:,:)=.8_dp*payload%state%fields(peer)%values(5,:,:,:)
        allocate(recovered(n(1),n(2),n(3)))
        call recover_reactive_temperatures_3d(species,payload%state%fields(peer)%values, &
          payload%temperature%fields(peer)%values(1,:,:,:),n(1),n(2),n(3),recovered,ok, &
          internal_energy_density=payload%internal_energy%fields(peer)%values(1,:,:,:))
        call require(ok,'retained EOS fixture')
        payload%temperature%fields(peer)%values(1,:,:,:)=recovered
        deallocate(recovered)
      enddo
      call write_amr_tree_checkpoint_3d('retained.chk',species,context,payload,ok,indexed=mod(variant,2)==1)
      call require(ok,'retained write with history/source/heat/index variants')
      actual=source
      call read_amr_tree_checkpoint_3d('retained.chk',species,context,[4,4,4],actual,ok)
      call require(.not.ok,'unsupported reader rejects retained payload')
      call read_amr_tree_checkpoint_3d('retained.chk',species,context,[4,4,4],actual,ok,allow_internal_energy=.true.)
      call require(ok,'retained checkpoint read')
      call same(actual,payload)
      call require(actual%energy_correction==-7._dp.and.actual%previous_timestep==1e-10_dp,'retained ledgers exact')
      call require(allocated(actual%pressure_work_change).eqv.(variant>=4),'pressure-work history presence exact')
      if(variant>=4) call require(actual%pressure_work_change==-.25_dp,'signed pressure-work history exact')
      do peer=1,size(payload%tree%nodes)
        call require(all(actual%internal_energy%fields(peer)%values==payload%internal_energy%fields(peer)%values), &
          'retained scalar roundtrip exact')
      enddo
      prior=actual
      open(newunit=from,file='retained.chk',status='old')
      open(newunit=to,file='retained_bad.chk',status='replace')
      read(from,'(a)') record;write(to,'(a)') trim(record)
      read(from,'(a)') record
      if(variant>=4) then
        write(to,'(a)') '-7.0 NaN'
      else
        write(to,'(a)') 'NaN'
      endif
      do
        read(from,'(a)',iostat=stat) record
        if(stat<0) exit
        call require(stat==0,'retained corrupt fixture copy')
        write(to,'(a)') trim(record)
      enddo
      close(from);close(to)
      call read_amr_tree_checkpoint_3d('retained_bad.chk',species,context,[4,4,4],actual,ok,allow_internal_energy=.true.)
      call require(.not.ok,'nonfinite retained correction rejected')
      call same(actual,prior)
      call require(actual%energy_correction==prior%energy_correction,'retained correction rollback')
      if(variant>=4) call require(actual%pressure_work_change==prior%pressure_work_change,'pressure-work read rollback')
      payload%temperature=source%temperature
      call write_amr_tree_checkpoint_3d('retained.chk',species,context,payload,ok)
      call require(.not.ok,'retained temperature inconsistency rejects writer')
      call read_amr_tree_checkpoint_3d('retained.chk',species,context,[4,4,4],actual,ok,allow_internal_energy=.true.)
      call require(ok,'retained valid checkpoint survives failed writer')
    enddo
  end subroutine check_internal_energy

  subroutine build_context(flow,text)
    type(reactive_3d_config), intent(in) :: flow
    character(len=:), allocatable, intent(inout) :: text
    call build_amr_tree_checkpoint_context_3d(flow,species,reactions,transport,x,'fixed-elementary','explicit', &
      'schedule-v1',text,ok)
    call require(ok,'canonical context')
  end subroutine
  subroutine same(a,b)
    type(amr_tree_checkpoint_3d), intent(in) :: a,b
    integer :: peer
    call require(a%time==b%time.and.a%steps==b%steps.and.a%regrids==b%regrids,'clock and counters')
    call require(a%maximum_reflux==b%maximum_reflux.and.a%minimum_theta==b%minimum_theta,'diagnostics')
    call require(a%maximum_remap_error==b%maximum_remap_error,'remap diagnostic')
    call require(all(a%initial==b%initial).and.all(a%boundary==b%boundary),'budgets')
    call require(allocated(a%source).eqv.allocated(b%source),'source presence')
    if (allocated(a%source)) call require(all(a%source==b%source),'source ledger exact')
    call require(allocated(a%reaction_source%fields).eqv.allocated(b%reaction_source%fields),'reaction presence')
    call require(allocated(a%heat_release%fields).eqv.allocated(b%heat_release%fields),'heat presence')
    call require(size(a%tree%nodes)==size(b%tree%nodes),'node count')
    do peer=1,size(a%tree%nodes)
      call require(a%tree%nodes(peer)%parent==b%tree%nodes(peer)%parent,'parent')
      call require(all(a%tree%nodes(peer)%shape==b%tree%nodes(peer)%shape),'shape')
      if(allocated(a%heat_release%fields)) &
        call require(all(a%heat_release%fields(peer)%values==b%heat_release%fields(peer)%values),'heat exact')
      if (allocated(a%reaction_source%fields)) &
        call require(all(a%reaction_source%fields(peer)%values==b%reaction_source%fields(peer)%values),'reaction exact')
      call require(all(a%state%fields(peer)%values==b%state%fields(peer)%values),'state exact')
      call require(all(a%temperature%fields(peer)%values==b%temperature%fields(peer)%values),'temperature exact')
    end do
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
