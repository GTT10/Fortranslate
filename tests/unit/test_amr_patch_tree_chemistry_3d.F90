program test_amr_patch_tree_chemistry_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: advance_reactive_chemistry_3d,react_sdc_state_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d,average_down_amr_patch_tree_3d
  use amr_patch_tree_chemistry_3d_mod, only: advance_amr_patch_tree_chemistry_3d,advance_amr_patch_tree_strang_3d
  use amr_patch_tree_chemistry_3d_mod, only: react_sdc_amr_patch_tree_3d
  use amr_reactive_3d_mod, only: element_integrals_from_reactive_integrals_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:),bad(:)
  type(amr_patch_tree_3d) :: tree
  type(amr_tree_state_3d) :: state,temperature,saved,saved_t,synchronized,oracle
  type(amr_patch_3d) :: patches(3)
  real(dp) :: root(12,2,2,2),rt(1,2,2,2),reference(12,2,2,2),ref_t(2,2,2),q(12),sound
  real(dp) :: initial(12),final(12),elements0(3),elements1(3),boundary_elements(3),budget(12),largest
  real(dp), parameter :: spacing(3)=[.001_dp,.001_dp,.001_dp],dt=1e-8_dp
  integer :: node,mode
  character(len=32) :: faces(6),integrator
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
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
  call check_sdc_stage()
  do mode=1,2
    integrator='explicit'
    if (mode==2) integrator='implicit'
    call prolong_amr_patch_tree_3d(tree,root,state,ok)
    call require(ok,'state')
    call prolong_amr_patch_tree_3d(tree,rt,temperature,ok)
    call require(ok,'temperature')
    call average_down_amr_patch_tree_3d(tree,state,ok)
    call require(ok,'initial source restriction')
    saved=state;saved_t=temperature
    reference=root;ref_t=rt(1,:,:,:)
    call advance_reactive_chemistry_3d(species,reactions,reference,ref_t,2,2,2,dt,2e-7_dp,1e-12_dp, &
      ok,chemistry_integrator=integrator)
    call require(ok,'single grid chemistry oracle')
    call prolong_amr_patch_tree_3d(tree,reference,oracle,ok)
    call require(ok,'oracle prolongation')
    call average_down_amr_patch_tree_3d(tree,oracle,ok)
    call require(ok,'oracle restriction')
    call advance_amr_patch_tree_chemistry_3d(tree,species,reactions,state,temperature,dt,2e-7_dp,1e-12_dp,ok,integrator)
    call require(ok,'tree chemistry')
    call require(maxval(abs(reference(6:,:,:,:)-root(6:,:,:,:)))>1e-12_dp,'active chemistry')
    do node=1,size(tree%nodes)
      call require(maxval(abs(state%fields(node)%values(1:5,:,:,:)-saved%fields(node)%values(1:5,:,:,:))) &
        <1e-10_dp,'Euler source invariants')
      call require(maxval(abs(state%fields(node)%values-oracle%fields(node)%values))<1e-10_dp,'source oracle')
      call require(maxval(abs(temperature%fields(node)%values-ref_t(1,1,1)))<1e-7_dp,'temperature oracle')
    end do
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
    if (mode==2) faces='outflow'
    call advance_amr_patch_tree_strang_3d(tree,species,reactions,state,temperature,spacing,dt,faces, &
      'pcm','mc','hllc',2e-7_dp,1e-12_dp,budget,largest,ok,integrator)
    call require(ok,'Strang step')
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
    saved=state;saved_t=temperature;budget=-17;largest=-19
    ! A hydro failure occurs after the first successful half reaction step.
    call advance_amr_patch_tree_strang_3d(tree,species,reactions,state,temperature,spacing,dt,faces, &
      'invalid','mc','hllc',2e-7_dp,1e-12_dp,budget,largest,ok,integrator)
    call require(.not.ok,'post-source hydro rejection')
    call unchanged()
    bad=reactions;bad(1)%reactant_stoich(1)=-1
    call advance_amr_patch_tree_chemistry_3d(tree,species,bad,state,temperature,dt,2e-7_dp,1e-12_dp,ok,integrator)
    call require(.not.ok,'invalid reactions rejected')
    call unchanged()
  end do
contains
  subroutine check_sdc_stage()
    type(amr_tree_state_3d) :: old,ot,rates,result,temps,ir,bad_old,saved_ir,prior,prior_t,prior_ir
    real(dp) :: expected(12,2,2,2),expected_t(2,2,2),expected_ir(8,2,2,2)
    integer :: peer
    call prolong_amr_patch_tree_3d(tree,root,old,ok)
    call require(ok,'SDC old tree')
    call prolong_amr_patch_tree_3d(tree,rt,ot,ok)
    call require(ok,'SDC temperature tree')
    call prolong_amr_patch_tree_3d(tree,.1_dp*root,rates,ok)
    call require(ok,'SDC external rates tree')
    expected=root;expected_t=rt(1,:,:,:);expected_ir=0
    call react_sdc_state_3d(species,reactions,root,rt(1,:,:,:),.1_dp*root,dt,2e-7_dp,1e-12_dp, &
      expected,expected_t,expected_ir,ok)
    call require(ok,'SDC uniform-grid reaction stage')
    call react_sdc_amr_patch_tree_3d(tree,species,reactions,old,ot,rates,dt,2e-7_dp,1e-12_dp,result,temps,ir,ok)
    call require(ok,'SDC three-level reaction stage')
    do peer=1,size(tree%nodes)
      call require(maxval(abs(result%fields(peer)%values(:,1,1,1)-expected(:,1,1,1)))<1e-10_dp, &
        'SDC hierarchy matches independent uniform-grid stage')
      call require(maxval(abs(temps%fields(peer)%values-expected_t(1,1,1)))<1e-7_dp,'SDC hierarchy temperatures')
      call require(maxval(abs(ir%fields(peer)%values(:,1,1,1)-expected_ir(:,1,1,1)))<1e-8_dp, &
        'SDC hierarchy reaction source estimates')
    end do
    prior=result;prior_t=temps;prior_ir=ir
    do peer=1,size(tree%nodes)
      result%fields(peer)%values(2,:,:,:)=result%fields(peer)%values(2,:,:,:)+1
      temps%fields(peer)%values=temps%fields(peer)%values+1
      ir%fields(peer)%values=ir%fields(peer)%values+1
    end do
    saved=result;saved_t=temps;saved_ir=ir
    call react_sdc_amr_patch_tree_3d(tree,species,reactions,old,ot,rates,dt,2e-7_dp,1e-12_dp, &
      result,temps,ir,ok,target_depth=tree%level_count()-1)
    call require(ok,'level-selective reaction stage')
    do peer=1,size(tree%nodes)
      if(tree%nodes(peer)%level==tree%level_count()-1) then
        call require(all(result%fields(peer)%values==prior%fields(peer)%values),'selected reaction state exact')
        call require(all(temps%fields(peer)%values==prior_t%fields(peer)%values),'selected reaction temperature exact')
        call require(all(ir%fields(peer)%values==prior_ir%fields(peer)%values),'selected reaction source exact')
      else
        call require(all(result%fields(peer)%values==saved%fields(peer)%values),'inactive reaction state retained')
        call require(all(temps%fields(peer)%values==saved_t%fields(peer)%values),'inactive reaction temperature retained')
        call require(all(ir%fields(peer)%values==saved_ir%fields(peer)%values),'inactive reaction source retained')
      end if
    end do
    bad_old=old;bad_old%fields(size(tree%nodes))%values(1,1,1,1)=-1
    saved=result;saved_t=temps;saved_ir=ir
    call react_sdc_amr_patch_tree_3d(tree,species,reactions,bad_old,ot,rates,dt,2e-7_dp,1e-12_dp,result,temps,ir,ok)
    call require(.not.ok,'SDC invalid late-patch state rejected')
    do peer=1,size(tree%nodes)
      call require(all(result%fields(peer)%values==saved%fields(peer)%values).and. &
        all(temps%fields(peer)%values==saved_t%fields(peer)%values).and. &
        all(ir%fields(peer)%values==saved_ir%fields(peer)%values),'SDC whole-tree rollback')
    end do
  end subroutine check_sdc_stage
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
    call require(all(budget==-17).and.largest==-19,'diagnostic rollback')
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
