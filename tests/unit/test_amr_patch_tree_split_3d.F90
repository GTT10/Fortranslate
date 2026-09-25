program test_amr_patch_tree_split_3d
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use amr_patch_tree_split_3d_mod, only: advance_amr_patch_tree_split_3d,advance_amr_patch_tree_sdc_3d
  use reactive_3d_mod, only: advance_reactive_sdc_3d
  use amr_patch_tree_hydro_3d_mod, only: advance_amr_patch_tree_hydro_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d,average_down_amr_patch_tree_3d
  use amr_patch_tree_chemistry_3d_mod, only: advance_amr_patch_tree_chemistry_3d,advance_amr_patch_tree_strang_3d
  use amr_reactive_3d_mod, only: element_integrals_from_reactive_integrals_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction), allocatable :: reactions(:),bad(:)
  type(amr_patch_tree_3d) :: tree
  type(amr_tree_state_3d) :: state,temperature,saved,saved_t,synchronized,oracle
  type(amr_patch_3d) :: patches(3)
  real(dp) :: root(12,2,2,2),rt(1,2,2,2),q(12),sound
  real(dp) :: initial(12),final(12),elements0(3),elements1(3),boundary_elements(3),budget(12),largest,theta, &
      reference_budget(12),reference_max
  real(dp), parameter :: spacing(3)=[.001_dp,.001_dp,.001_dp],dt=1e-8_dp
  real(dp) :: outlet_pressure(6)
  integer :: node,mode
  character(len=32) :: faces(6),integrator
  character(len=16) :: argument
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  call load_h2o2_elementary_mechanism(reactions,ok)
  call require(ok,'reactions')
  q=0;q(1)=.3_dp;q(5)=202650;q(2:4)=[5.0_dp,2.0_dp,-1.0_dp]
  q(6)=.028_dp;q(7)=1e-5_dp;q(9)=.22_dp;q(12)=1-q(6)-q(7)-q(9)
  call reactive_primitive_to_conserved(species,q,root(:,1,1,1),rt(1,1,1,1),sound,ok)
  call require(ok,'primitive')
  do node=1,12
    root(node,:,:,:)=root(node,1,1,1)
  end do
  rt=rt(1,1,1,1)
  call initialize_amr_patch_3d(2,2,2,1,1,1,2,1,2,2,patches(1),ok)
  call initialize_amr_patch_3d(2,2,2,2,2,1,2,1,2,2,patches(2),ok)
  call initialize_amr_patch_3d(2,4,4,1,1,1,1,1,1,2,patches(3),ok)
  call initialize_amr_patch_tree_3d([2,2,2],[1,1,2],patches,tree,ok)
  call require(ok,'three-level branched tree')
  call get_command_argument(1,argument)
  if(trim(argument)=='sdc') then
    call check_sdc()
    print *, 'test_amr_patch_tree_sdc_3d: PASS'
    stop
  end if
  do mode=1,4
    integrator='explicit'
    if (mod(mode,2)==0) integrator='implicit'
    call prolong_amr_patch_tree_3d(tree,root,state,ok)
    call require(ok,'state')
    call prolong_amr_patch_tree_3d(tree,rt,temperature,ok)
    call require(ok,'temperature')
    call average_down_amr_patch_tree_3d(tree,state,ok)
    call require(ok,'initial source restriction')
    saved=state;saved_t=temperature
    ! Make one exposed fine cell nonuniform before the coupled flow/source step.
    state=saved;temperature=saved_t
    q(5)=210000
    call reactive_primitive_to_conserved(species,q,state%fields(4)%values(:,1,1,1), &
      temperature%fields(4)%values(1,1,1,1),sound,ok)
    call require(ok,'nonuniform leaf')
    call average_down_amr_patch_tree_3d(tree,state,ok)
    call require(ok,'initial restriction')
    call inventory(initial,elements0)
    faces='periodic'
    if (mode>=2) faces='outflow'
    outlet_pressure=0
    if (mode>=3) outlet_pressure([2,6])=203650
    if (mode==1) then
      oracle=state;synchronized=temperature
      call advance_amr_patch_tree_strang_3d(tree,species,reactions,oracle,synchronized,spacing,dt,faces, &
        'characteristic_plm','mc','hllc',2e-7_dp,1e-12_dp,reference_budget,reference_max,ok,integrator)
      call require(ok,'no-transport comparison')
    end if
    call advance_amr_patch_tree_split_3d(tree,species,reactions,transport,state,temperature,spacing,dt,faces, &
      'characteristic_plm','mc','hllc',.true.,.true.,.true.,.true.,.true.,2e-7_dp,1e-12_dp,budget,largest,theta,ok, &
      integrator,subcycling=mod(mode,2)==0,outlet_pressure=outlet_pressure)
    call require(ok,'Strang step')
    if (mode==1) call require(maxval(abs(state%fields(4)%values-oracle%fields(4)%values))>1e-10_dp, &
      'transport changes the coupled solution')
    call inventory(final,elements1)
    call element_integrals_from_reactive_integrals_3d(species,budget,boundary_elements,ok)
    call require(ok,'element boundary ledger')
    call require(maxval(abs(final(1:5)-initial(1:5)-budget(1:5))/ &
      max(1e-20_dp,abs(initial(1:5))))<2e-10_dp,'Euler inventory plus boundary')
    call require(maxval(abs(elements1-elements0-boundary_elements)/max(1e-30_dp,abs(elements0))) &
      <2e-10_dp,'element inventory plus boundary')
    synchronized=state
    call average_down_amr_patch_tree_3d(tree,synchronized,ok)
    call require(ok,'final restriction')
    do node=1,size(tree%nodes)
      call require(all(synchronized%fields(node)%values==state%fields(node)%values),'synchronized state')
      call require(maxval(abs(sum(state%fields(node)%values(6:,:,:,:),dim=1)- &
        state%fields(node)%values(1,:,:,:)))<1e-12_dp,'species closure')
    end do
    saved=state;saved_t=temperature;budget=-17;largest=-19;theta=-23
    ! A hydro failure occurs after the first successful half reaction step.
    call advance_amr_patch_tree_split_3d(tree,species,reactions,transport,state,temperature,spacing,dt,faces, &
      'invalid','mc','hllc',.true.,.true.,.true.,.true.,.true.,2e-7_dp,1e-12_dp,budget,largest,theta,ok, &
      integrator,subcycling=mod(mode,2)==0,outlet_pressure=outlet_pressure)
    call require(.not.ok,'post-source hydro rejection')
    call unchanged()
    bad=reactions;bad(1)%reactant_stoich(1)=-1
    call advance_amr_patch_tree_chemistry_3d(tree,species,bad,state,temperature,dt,2e-7_dp,1e-12_dp,ok,integrator)
    call require(.not.ok,'invalid reactions rejected')
    call unchanged()
    ! With transport disabled the integrated entry point must reduce exactly
    ! to the pre-existing reaction/hydro composition.
    if (mode<=2) then
      oracle=state;synchronized=temperature
      call advance_amr_patch_tree_strang_3d(tree,species,reactions,oracle,synchronized,spacing,dt,faces, &
        'pcm','mc','hllc',2e-7_dp,1e-12_dp,reference_budget,reference_max,ok,integrator)
      call require(ok,'reaction/hydro reference')
      call advance_amr_patch_tree_split_3d(tree,species,reactions,transport,state,temperature,spacing,dt,faces, &
        'pcm','mc','hllc',.true.,.false.,.false.,.false.,.false.,2e-7_dp,1e-12_dp,budget,largest,theta,ok,integrator)
      call require(ok,'disabled transport')
      call same_as_reference()
      call require(theta==1,'inactive transport limiter')
    end if
    oracle=state;synchronized=temperature
    call advance_amr_patch_tree_hydro_3d(tree,species,oracle,synchronized,spacing,dt,faces, &
      'pcm','mc','hllc',reference_budget,reference_max,ok,outlet_pressure=outlet_pressure)
    call require(ok,'hydro reference')
    call advance_amr_patch_tree_split_3d(tree,species,bad,transport(:0),state,temperature,spacing,dt,faces, &
      'pcm','mc','hllc',.false.,.false.,.false.,.false.,.false.,2e-7_dp,1e-12_dp,budget,largest,theta,ok, &
      outlet_pressure=outlet_pressure)
    call require(ok,'disabled chemistry and transport')
    call same_as_reference()
  end do
contains
  subroutine check_sdc()
    type(amr_tree_state_3d) :: u,t,previous,previous_t,sync
    real(dp) :: reference(12,2,2,2),ref_t(2,2,2),sources(12),bound(12),maximum,minimum
    real(dp) :: before(12),after(12),e0(3),e1(3),ext(3),source_ref(12),boundary_ref(12)
    real(dp), parameter :: h=1e-9_dp,force(3)=[20000._dp,-30000._dp,40000._dp]
    integer :: peer,trial
    faces='periodic'
    call prolong_amr_patch_tree_3d(tree,root,u,ok)
    call require(ok,'SDC tree state')
    call prolong_amr_patch_tree_3d(tree,rt,t,ok)
    call require(ok,'SDC tree temperature')
    reference=root;ref_t=rt(1,:,:,:)
    call advance_reactive_sdc_3d(species,reference,ref_t,spacing,h,0._dp,2,2,ok=ok, &
      external_force=force,transport=transport,reactions=reactions,chemistry_relative_tolerance=2e-7_dp, &
      chemistry_absolute_tolerance=1e-12_dp,source_change=source_ref,boundary_change=boundary_ref)
    call require(ok,'SDC single-grid reacting reference')
    do trial=1,2
      if(trial==2) u%fields(4)%values(5,1,1,1)=1.001_dp*u%fields(4)%values(5,1,1,1)
      call average_down_amr_patch_tree_3d(tree,u,ok)
      call require(ok,'SDC input restriction')
      call composite_integrals_amr_patch_tree_3d(tree,u,spacing,before,ok)
      call require(ok,'SDC initial composite inventory')
      call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,u,t,spacing,h,faces,2,.true., &
        [.true.,.true.,.true.,.true.],force,2e-7_dp,1e-12_dp,bound,sources,maximum,minimum,ok)
      call require(ok,'SDC coupled hierarchy step')
      if(trial==1) then
        do peer=1,size(tree%nodes)
          call require(maxval(abs(u%fields(peer)%values(:,1,1,1)-reference(:,1,1,1))/ &
            max(1._dp,abs(reference(:,1,1,1))))<1e-11_dp,'SDC hierarchy uniform reacting reduction')
          call require(maxval(abs(t%fields(peer)%values-ref_t(1,1,1)))<1e-6_dp,'SDC hierarchy uniform temperature')
        end do
      end if
      call composite_integrals_amr_patch_tree_3d(tree,u,spacing,after,ok)
      call require(ok,'SDC final composite inventory')
      call require(maxval(abs(after(:5)-before(:5)-bound(:5)-sources(:5))/max(1e-20_dp,abs(before(:5)))) &
        <2e-10_dp,'SDC Euler inventory with forcing and reflux')
      call element_integrals_from_reactive_integrals_3d(species,before,e0,ok)
      call require(ok,'SDC initial elements')
      call element_integrals_from_reactive_integrals_3d(species,after,e1,ok)
      call require(ok,'SDC final elements')
      call element_integrals_from_reactive_integrals_3d(species,bound+sources,ext,ok)
      call require(ok.and.maxval(abs(e1-e0-ext)/max(1e-30_dp,abs(e0)))<2e-10_dp,'SDC element inventory')
      sync=u
      call average_down_amr_patch_tree_3d(tree,sync,ok)
      call require(ok,'SDC final restriction check')
      do peer=1,size(tree%nodes)
        call require(all(sync%fields(peer)%values==u%fields(peer)%values),'SDC synchronized parents')
      end do
    end do
    previous=u;previous_t=t;bound=-3;sources=-4;maximum=-5;minimum=-6
    call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,u,t,spacing,h,faces,2,.true., &
      [.true.,.true.,.true.,.true.],force,2e-7_dp,1e-12_dp,bound,sources,maximum,minimum,ok,difmag=-1._dp)
    call require(.not.ok.and.all(bound==-3).and.all(sources==-4).and.maximum==-5.and.minimum==-6, &
      'SDC failed hydro stage preserves inventories')
    do peer=1,size(tree%nodes)
      call require(all(u%fields(peer)%values==previous%fields(peer)%values).and. &
        all(t%fields(peer)%values==previous_t%fields(peer)%values),'SDC whole-step rollback')
    end do
  end subroutine check_sdc
  subroutine same_as_reference()
    integer :: peer
    do peer=1,size(tree%nodes)
      call require(all(state%fields(peer)%values==oracle%fields(peer)%values),'disabled operator state parity')
      call require(all(temperature%fields(peer)%values==synchronized%fields(peer)%values),'disabled operator T parity')
    end do
    call require(all(budget==reference_budget).and.largest==reference_max,'disabled operator diagnostics')
  end subroutine
  subroutine inventory(total,elements)
    real(dp), intent(out) :: total(12),elements(3)
    call composite_integrals_amr_patch_tree_3d(tree,state,spacing,total,ok)
    call require(ok,'composite inventory')
    call element_integrals_from_reactive_integrals_3d(species,total,elements,ok)
    call require(ok,'element inventory')
  end subroutine
  subroutine unchanged()
    integer :: peer
    do peer=1,size(tree%nodes)
      call require(all(state%fields(peer)%values==saved%fields(peer)%values),'state rollback')
      call require(all(temperature%fields(peer)%values==saved_t%fields(peer)%values),'temperature rollback')
    end do
    call require(all(budget==-17).and.largest==-19.and.theta==-23,'diagnostic rollback')
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
