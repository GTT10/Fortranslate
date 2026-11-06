program test_amr_patch_tree_split_3d
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_is_nan
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use amr_patch_tree_split_3d_mod, only: advance_amr_patch_tree_split_3d,advance_amr_patch_tree_sdc_3d, &
    advance_amr_patch_tree_sdc_subcycle_3d
  use reactive_3d_mod, only: advance_reactive_sdc_3d
  use amr_patch_tree_hydro_3d_mod, only: advance_amr_patch_tree_hydro_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use precision_mod, only: dp
  use dual_energy_mod, only: dual_energy_policy
  use reactive_mol_3d_mod, only: advance_reactive_mol_3d
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
  logical :: ok,fail_auxiliary_stage=.false.
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
  if(trim(argument)=='auxiliary') then
    call check_auxiliary_mol_time()
    print *, 'test_amr_stored_auxiliary_mol_time: PASS'
    stop
  endif
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
  subroutine check_auxiliary_mol_time()
    type(amr_patch_tree_3d) :: single
    type(amr_tree_state_3d) :: u,t,held,old_u,old_t,old_held,source0,source1
    type(amr_tree_flux_3d), allocatable :: flux(:),aflux(:),old_flux(:),old_aflux(:)
    type(dual_energy_policy) :: policy
    real(dp) :: seed(12,2,2,2),seed_t(1,2,2,2),prim(12),c,aux_seed(1,2,2,2)
    real(dp) :: regular(12,2,2,2),regular_t(2,2,2),regular_aux(2,2,2),ledger,regular_ledger
    real(dp) :: bound(12),external(12),aux_change(12),maximum,minimum,before(12),after(12)
    real(dp) :: expected_change,scalar_before(1),scalar_after(1)
    integer(int64) :: steps(3)
    integer :: peer,mode,cycle
    logical :: cycling,transport_flags(4)
    procedure(inert_reactor), pointer :: reactor
    character(len=32) :: periodic(6)
    periodic='periodic';prim=0;prim(1)=.3_dp;prim(5)=202650;prim(12)=1
    call reactive_primitive_to_conserved(species,prim,seed(:,1,1,1),seed_t(1,1,1,1),c,ok)
    call require(ok,'auxiliary uniform N2 initial state')
    do peer=1,12
      seed(peer,:,:,:)=seed(peer,1,1,1)
    enddo
    seed_t=seed_t(1,1,1,1);aux_seed(1,:,:,:)=seed(5,:,:,:)-100
    call initialize_amr_patch_tree_3d([2,2,2],[integer::],patches(:0),single,ok)
    call require(ok,'auxiliary root-only topology')
    call prolong_amr_patch_tree_3d(single,seed,u,ok)
    call require(ok,'auxiliary root state')
    call prolong_amr_patch_tree_3d(single,seed_t,t,ok)
    call require(ok,'auxiliary root temperature')
    call prolong_amr_patch_tree_3d(single,aux_seed,held,ok)
    call require(ok,'auxiliary root retained field')
    ! Ordinary MOL must reject invalid old thermodynamics without committing
    ! state, temperature, held energy or scalar diagnostics.
    old_u=u;old_t=t;old_held=held
    policy%eta2=2;policy%update_total_from_internal=.true.;policy%allow_negative_energy=.false.
    do mode=1,4
      u=old_u;t=old_t;held=old_held
      select case(mode)
      case(1)
        u%fields(1)%values(1,1,1,1)=-1
      case(2)
        u%fields(1)%values(5,1,1,1)=ieee_value(0._dp,ieee_quiet_nan)
      case(3)
        held%fields(1)%values(1,1,1,1)=ieee_value(0._dp,ieee_quiet_nan)
      case(4)
        t%fields(1)%values(1,1,1,1)=ieee_value(0._dp,ieee_quiet_nan)
      end select
      ledger=-7;bound=-3;external=-4;maximum=-5;minimum=-6
      call advance_amr_patch_tree_sdc_subcycle_3d(single,species,reactions,transport,u,t,spacing,dt,periodic, &
        1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp, &
        bound,external,maximum,minimum,ok,mol_mode=.true.,internal_energy=held, &
        energy_correction=ledger,energy_policy=policy)
      call require(.not.ok,'ordinary MOL rejects invalid old state')
      call require(ledger==-7.and.all(bound==-3).and.all(external==-4).and.maximum==-5.and.minimum==-6, &
        'ordinary MOL admission diagnostic rollback')
      select case(mode)
      case(1)
        call require(u%fields(1)%values(1,1,1,1)==-1,'invalid density preserved on refusal')
        u%fields(1)%values(1,1,1,1)=old_u%fields(1)%values(1,1,1,1)
      case(2)
        call require(ieee_is_nan(u%fields(1)%values(5,1,1,1)),'invalid total energy preserved on refusal')
        u%fields(1)%values(5,1,1,1)=old_u%fields(1)%values(5,1,1,1)
      case(3)
        call require(ieee_is_nan(held%fields(1)%values(1,1,1,1)),'invalid held energy preserved on refusal')
        held=old_held
      case(4)
        call require(ieee_is_nan(t%fields(1)%values(1,1,1,1)),'invalid temperature preserved on refusal')
        t=old_t
      end select
      call require(all(u%fields(1)%values==old_u%fields(1)%values).and. &
        all(t%fields(1)%values==old_t%fields(1)%values).and. &
        all(held%fields(1)%values==old_held%fields(1)%values),'ordinary MOL admission field rollback')
    enddo
    u=old_u;t=old_t;held=old_held
    policy%eta2=1;policy%update_total_from_internal=.false.;policy%allow_negative_energy=.false.
    regular=seed;regular_t=seed_t(1,:,:,:);regular_aux=aux_seed(1,:,:,:);regular_ledger=0;ledger=0
    call advance_reactive_mol_3d(species,regular,regular_t,spacing,dt,0._dp,1,.false.,2e-7_dp,1e-12_dp,ok, &
      internal_energy_density=regular_aux,energy_correction=regular_ledger,energy_policy=policy)
    call require(ok,'auxiliary standalone regular MOL oracle')
    call advance_amr_patch_tree_sdc_3d(single,species,reactions,transport,u,t,spacing,dt,periodic,1,.false., &
      [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp,bound,external,maximum,minimum,ok, &
      mol_mode=.true.,internal_energy=held,energy_correction=ledger,energy_policy=policy, &
      integrated_flux=flux,internal_energy_integrated_flux=aflux)
    call require(ok,'auxiliary root MOL time integration')
    call require(all(u%fields(1)%values==regular),'auxiliary root state matches regular MOL')
    call require(all(t%fields(1)%values(1,:,:,:)==regular_t),'auxiliary root temperature matches regular MOL')
    call require(all(held%fields(1)%values(1,:,:,:)==regular_aux).and.ledger==regular_ledger, &
      'auxiliary root retained scalar and correction match regular MOL')
    call require(all(held%fields(1)%values==aux_seed),'independent scalar survives MOL stage reset policy')
    ! A pressure perturbation crosses the finest coarse/fine interface. The
    ! independently stored scalar must conserve through actual nonzero fluxes.
    call prolong_amr_patch_tree_3d(tree,seed,u,ok);call require(ok,'nonuniform auxiliary state')
    call prolong_amr_patch_tree_3d(tree,seed_t,t,ok);call require(ok,'nonuniform auxiliary temperature')
    call prolong_amr_patch_tree_3d(tree,aux_seed,held,ok);call require(ok,'nonuniform auxiliary scalar')
    prim(5)=1.02_dp*prim(5)
    peer=size(tree%nodes)
    call reactive_primitive_to_conserved(species,prim,u%fields(peer)%values(:,1,1,1), &
      t%fields(peer)%values(1,1,1,1),c,ok)
    call require(ok,'finest pressure perturbation')
    held%fields(peer)%values(1,1,1,1)=u%fields(peer)%values(5,1,1,1)-100
    old_held=held
    call composite_integrals_amr_patch_tree_3d(tree,held,spacing,scalar_before,ok)
    call require(ok,'nonuniform scalar initial integral')
    ledger=0
    call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,u,t,spacing,dt,periodic, &
      1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp, &
      bound,external,maximum,minimum,ok,mol_mode=.true.,subcycling=.true., &
      internal_energy=held,energy_correction=ledger,energy_policy=policy,internal_energy_integrated_flux=aflux)
    call require(ok,'nonuniform scalar recursive transport and reflux')
    call composite_integrals_amr_patch_tree_3d(tree,held,spacing,scalar_after,ok)
    call require(ok.and.abs(scalar_after(1)-scalar_before(1))<1e-13_dp,'nonuniform scalar composite conservation')
    call require(maxval(abs(held%fields(peer)%values-old_held%fields(peer)%values))>1e-8_dp, &
      'nonuniform scalar flux changes finest state')
    call require(ledger==0,'nonuniform scalar retained policy makes no total energy repair')
    policy%update_total_from_internal=.true.;aux_seed(1,:,:,:)=seed(5,:,:,:)
    expected_change=10*8*product(spacing)
    do cycle=0,1
      cycling=cycle==1
      do mode=1,3
        call prolong_amr_patch_tree_3d(tree,seed,u,ok);call require(ok,'auxiliary nested state')
        call prolong_amr_patch_tree_3d(tree,seed_t,t,ok);call require(ok,'auxiliary nested temperature')
        call prolong_amr_patch_tree_3d(tree,aux_seed,held,ok);call require(ok,'auxiliary nested held field')
        call composite_integrals_amr_patch_tree_3d(tree,u,spacing,before,ok)
        call require(ok,'auxiliary nested initial inventory')
        call composite_integrals_amr_patch_tree_3d(tree,held,spacing,scalar_before,ok)
        call require(ok,'auxiliary scalar initial inventory')
        ledger=0;steps=-1;nullify(reactor)
        transport_flags=cycling.and.mode==1
        if(mode>1) reactor=>inert_reactor
        call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,u,t,spacing,dt,periodic, &
          merge(2,1,mode==3),mode>1,transport_flags,[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp, &
          bound,external,maximum,minimum,ok,mol_mode=.true.,subcycling=cycling,level_steps=steps, &
          internal_energy=held,energy_correction=ledger,energy_policy=policy, &
          auxiliary_source=energy_source_recursive,auxiliary_change=aux_change,evaluate_reaction=reactor, &
          integrated_flux=flux,internal_energy_integrated_flux=aflux)
        call require(ok,'auxiliary nested MOL common/recursive forced integration')
        call require(all(steps==merge([1_int64,2_int64,4_int64],[1_int64,1_int64,1_int64],cycling)), &
          'auxiliary nested MOL level schedule')
        call composite_integrals_amr_patch_tree_3d(tree,u,spacing,after,ok)
        call require(ok,'auxiliary final conserved inventory')
        call composite_integrals_amr_patch_tree_3d(tree,held,spacing,scalar_after,ok)
        call require(ok,'auxiliary final scalar inventory')
        call require(abs(after(5)-before(5)-bound(5)-external(5)-aux_change(5)-ledger)<1e-13_dp, &
          'auxiliary composite total energy plus correction ledger')
        if(mode==1) then
          call require(abs(ledger-expected_change)<1e-13_dp,'leaf-volume weighted predictor/final corrections')
          call require(abs(scalar_after(1)-scalar_before(1))<1e-13_dp,'auxiliary scalar conservative reflux')
        else
          call require(abs(scalar_after(1)-after(5))<1e-13_dp,'reactor endpoint scalar and total energy agree')
        endif
        old_u=u;old_t=t;old_held=held;old_flux=flux;old_aflux=aflux
        fail_auxiliary_stage=.true.;bound=-3;external=-4;maximum=-5;minimum=-6;ledger=-7;aux_change=-8;steps=-9
        call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,u,t,spacing,dt,periodic, &
          1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp, &
          bound,external,maximum,minimum,ok,mol_mode=.true.,subcycling=cycling,level_steps=steps, &
          internal_energy=held,energy_correction=ledger,energy_policy=policy, &
          auxiliary_source=energy_source_recursive,auxiliary_change=aux_change, &
          integrated_flux=flux,internal_energy_integrated_flux=aflux)
        fail_auxiliary_stage=.false.
        call require(.not.ok.and.ledger==-7.and.all(bound==-3).and.all(external==-4).and. &
          maximum==-5.and.minimum==-6.and.all(aux_change==-8).and.all(steps==-9),'auxiliary late stage ledger rollback')
        do peer=1,size(tree%nodes)
          call require(all(u%fields(peer)%values==old_u%fields(peer)%values).and. &
            all(t%fields(peer)%values==old_t%fields(peer)%values).and. &
            all(held%fields(peer)%values==old_held%fields(peer)%values),'auxiliary late stage field rollback')
          call require(all(flux(peer)%x==old_flux(peer)%x).and.all(aflux(peer)%x==old_aflux(peer)%x), &
            'auxiliary late stage flux rollback')
        enddo
      enddo
    enddo
  end subroutine check_auxiliary_mol_time

  subroutine energy_source_recursive(gas,temp,interval,stage,depth,time,rate,success)
    type(amr_tree_state_3d), intent(in) :: gas,temp
    type(amr_tree_state_3d), intent(inout) :: rate
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    logical, intent(out) :: success
    integer :: peer
    success=.false.
    if(fail_auxiliary_stage.and.stage>0) return
    rate=gas
    do peer=1,size(tree%nodes)
      rate%fields(peer)%values=0
      if(tree%nodes(peer)%level/=depth) cycle
      rate%fields(peer)%values(5,:,:,:)=-10/dt
    enddo
    success=.true.
  end subroutine energy_source_recursive

  subroutine inert_reactor(sp,rx,rhoe,external_y,external_e,duration,rtol,atol,densities,temp,count,good,maximum_steps)
    type(nasa7_species), intent(in) :: sp(:)
    type(elementary_reaction), intent(in) :: rx(:)
    real(dp), intent(in) :: rhoe,external_y(:),external_e,duration,rtol,atol
    real(dp), intent(inout) :: densities(:),temp
    integer, intent(out) :: count
    logical, intent(out) :: good
    integer, optional, intent(in) :: maximum_steps
    densities=densities+duration*external_y
    ! The shared SDC adapter recovers accepted temperature from the forced
    ! energy endpoint after this inert callback returns unchanged chemistry.
    count=1;good=.true.
  end subroutine inert_reactor

  subroutine check_sdc()
    type(amr_tree_state_3d) :: u,t,previous,previous_t,sync,coarse0,coarse1,op,np,ir
    type(amr_tree_flux_3d), allocatable :: flux(:)
    real(dp) :: reference(12,2,2,2),ref_t(2,2,2),sources(12),bound(12),maximum,minimum
    real(dp) :: before(12),after(12),e0(3),e1(3),ext(3),source_ref(12),boundary_ref(12)
    real(dp), parameter :: h=1e-9_dp,force(3)=[20000._dp,-30000._dp,40000._dp]
    integer :: peer,trial,dims(3)
    integer(int64) :: steps(3)
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
    previous=u;previous_t=t
    allocate(coarse0%fields(size(tree%nodes)),ir%fields(size(tree%nodes)))
    do peer=1,size(tree%nodes)
      dims=tree%nodes(peer)%shape
      allocate(coarse0%fields(peer)%values(13,dims(1),dims(2),dims(3)))
      coarse0%fields(peer)%values(:12,:,:,:)=u%fields(peer)%values
      coarse0%fields(peer)%values(13,:,:,:)=t%fields(peer)%values(1,:,:,:)
      allocate(ir%fields(peer)%values(8,dims(1),dims(2),dims(3)));ir%fields(peer)%values=0
    end do
    call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,u,t,spacing,h,faces,2,.true., &
      [.true.,.true.,.true.,.true.],force,2e-7_dp,1e-12_dp,bound,sources,maximum,minimum,ok, &
      reaction_source=ir,target_depth=0,integrated_flux=flux,old_predictor=op,new_predictor=np)
    call require(ok,'root level coupled SDC iteration')
    call require(maxval(abs(u%fields(1)%values(:,1,1,1)-reference(:,1,1,1))/ &
      max(1._dp,abs(reference(:,1,1,1))))<1e-11_dp,'root level SDC uniform reference')
    do peer=2,size(tree%nodes)
      call require(all(u%fields(peer)%values==previous%fields(peer)%values),'root step retains finer state')
      call require(all(t%fields(peer)%values==previous_t%fields(peer)%values),'root step retains finer temperature')
      call require(.not.allocated(flux(peer)%x),'root step exports no fine flux')
    end do
    coarse1=coarse0
    coarse1%fields(1)%values(:12,:,:,:)=u%fields(1)%values
    coarse1%fields(1)%values(13,:,:,:)=t%fields(1)%values(1,:,:,:)
    previous=u;previous_t=t
    call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,u,t,spacing,.5_dp*h,faces,2,.true., &
      [.true.,.true.,.true.,.true.],force,2e-7_dp,1e-12_dp,bound,sources,maximum,minimum,ok, &
      reaction_source=ir,target_depth=1,coarse_start=coarse0,coarse_end=coarse1, &
      fractions_start=[0._dp],fractions_end=[.5_dp],source_start=op,source_end=np,integrated_flux=flux)
    call require(ok,'fine coupled SDC interval with parent temporal state and sources')
    call require(all(bound==0).and.maximum==0,'fine interval defers physical budget and reflux')
    do peer=1,size(tree%nodes)
      if(tree%nodes(peer)%level==1) cycle
      call require(all(u%fields(peer)%values==previous%fields(peer)%values),'fine interval retains other levels')
      call require(all(t%fields(peer)%values==previous_t%fields(peer)%values),'fine interval retains other temperatures')
    end do
    call prolong_amr_patch_tree_3d(tree,root,u,ok)
    call prolong_amr_patch_tree_3d(tree,rt,t,ok)
    do trial=1,3
      if(trial==2) u%fields(4)%values(5,1,1,1)=1.001_dp*u%fields(4)%values(5,1,1,1)
      call average_down_amr_patch_tree_3d(tree,u,ok)
      call require(ok,'SDC input restriction')
      call composite_integrals_amr_patch_tree_3d(tree,u,spacing,before,ok)
      call require(ok,'SDC initial composite inventory')
      if(trial==3) then
        call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,u,t,spacing,h,faces,2,.true., &
          [.true.,.true.,.true.,.true.],force,2e-7_dp,1e-12_dp,bound,sources,maximum,minimum,ok,level_steps=steps)
        call require(ok.and.all(steps==[1_int64,2_int64,4_int64]),'recursive SDC level counts')
      else
      call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,u,t,spacing,h,faces,2,.true., &
        [.true.,.true.,.true.,.true.],force,2e-7_dp,1e-12_dp,bound,sources,maximum,minimum,ok)
      end if
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
    steps=-7
    call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,u,t,spacing,h,faces,2,.true., &
      [.true.,.true.,.true.,.true.],force,2e-7_dp,1e-12_dp,bound,sources,maximum,minimum,ok, &
      difmag=-1._dp,level_steps=steps)
    call require(.not.ok.and.all(bound==-3).and.all(sources==-4).and.maximum==-5.and.minimum==-6.and.all(steps==-7), &
      'recursive SDC failure preserves diagnostics and counts')
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
