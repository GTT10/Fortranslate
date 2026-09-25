program test_amr_patch_tree_checkpoint_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use simulation_config_reactive_3d_mod, only: reactive_3d_config,reactive_3d_mole_fractions
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: initialize_amr_patch_tree_3d,prolong_amr_patch_tree_3d, &
    average_down_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
  use amr_patch_tree_checkpoint_3d_mod, only: amr_tree_checkpoint_3d,build_amr_tree_checkpoint_context_3d, &
    write_amr_tree_checkpoint_3d,read_amr_tree_checkpoint_3d
  implicit none
  type(amr_tree_checkpoint_3d) :: source,loaded,saved
  type(reactive_3d_config) :: config,changed
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_3d) :: patches(2)
  character(len=:), allocatable :: context,other
  real(dp) :: x(7),u(12,4,4,4),t(1,4,4,4),q(12),sound
  integer :: c,unit,input,status
  character(len=65536) :: line
  logical :: ok
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
  call require(context==other,'automatic and zero non-PPM coefficient identity')
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
  call write_amr_tree_checkpoint_3d('tree.chk',species,context,source,ok)
  call require(ok,'write snapshot')
  loaded=source;loaded%state%fields(1)%values=-1;loaded%time=0
  call read_amr_tree_checkpoint_3d('tree.chk',species,context,[4,4,4],loaded,ok)
  call require(ok,'read snapshot')
  call same(loaded,source)
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
contains
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
    call require(size(a%tree%nodes)==size(b%tree%nodes),'node count')
    do peer=1,size(a%tree%nodes)
      call require(a%tree%nodes(peer)%parent==b%tree%nodes(peer)%parent,'parent')
      call require(all(a%tree%nodes(peer)%shape==b%tree%nodes(peer)%shape),'shape')
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
