program test_eb_tree_checkpoint_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use dual_energy_mod, only: dual_energy_policy
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,advance_reactive_eb_tree_hydro_3d
  use eb_reactive_tree_checkpoint_3d_mod, only: eb_tree_checkpoint_3d,reactive_eb_tree_checkpoint_record_3d
  use eb_reactive_tree_checkpoint_3d_mod, only: checkpoint_held_restriction_matches_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,middle,coarse
  type(eb_reactive_tree_3d) :: tree,saved,reference
  type(reactive_eb_3d_config) :: config
  type(eb_tree_checkpoint_3d) :: record,restored
  integer :: stream
  character(len=512) :: message
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:)
  real(dp), allocatable :: before(:),boundary(:),primitive(:),wall(:),ref_b(:),ref_w(:)
  real(dp) :: speed
  character(len=1), parameter :: axes(3)=['x','y','z']
  integer :: nv,node,axis,count,i,j,k
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  call load_h2o2_elementary_mechanism(reactions,ok)
  call require(ok,'reactions')
  nv=reactive_nvar(size(species))
  allocate(before(nv),boundary(nv),primitive(reactive_nprim(size(species))))
  allocate(wall(nv),ref_b(nv),ref_w(nv))
  do axis=3,3
  do count=5,5
  call initialize_amr_patch_3d(4,4,4,1,2,1,2,2,3,2,boxes(1),ok)
  call require(ok,'child left')
  call initialize_amr_patch_3d(4,4,4,3,4,1,2,2,3,2,boxes(2),ok)
  call require(ok,'child right')
  call initialize_amr_patch_3d(4,4,4,3,4,3,4,2,3,2,boxes(3),ok)
  call require(ok,'grandchild left')
  call initialize_amr_patch_3d(4,4,4,1,2,3,4,2,3,2,boxes(4),ok)
  call require(ok,'grandchild right')
  if (count==5) then
    call initialize_amr_patch_tree_3d([4,4,4],[1,1,2,3],boxes,topology,ok)
  else
    call initialize_amr_patch_tree_3d([4,4,4],[1,1,2],boxes(:3),topology,ok)
  end if
  call require(ok,'topology')
  call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,axes(axis),.31_dp,finest,ok)
  call require(ok,'finest geometry')
  call coarsen_eb_geometry_3d(finest,2,middle,ok)
  call require(ok,'middle geometry')
  call coarsen_eb_geometry_3d(middle,2,coarse,ok)
  call require(ok,'coarse geometry')
  allocate(u(nv,4,4,4),t(4,4,4))
  config%nx=4;config%ny=4;config%nz=4
  config%redistribution='none';config%plane_axis='z';config%plane_position=.31_dp
  call initialize_reactive_eb_density_sheet_3d(species,config,coarse,u,t,ok)
  call require(ok,'root state')
  call initialize_reactive_eb_tree_3d(species,topology,finest,u,t,tree,ok)
  call require(ok,'tree initialization')
  do node=2,count
    tree%nodes(node)%state=real(node,dp)*tree%nodes(node)%state
    tree%nodes(node)%state(5,:,:,:)=1.1_dp*tree%nodes(node)%state(5,:,:,:)
    do k=1,4
    do j=1,4
    do i=1,4
      if (tree%nodes(node)%geometry%volume_fraction(i,j,k)<=0._dp) cycle
      call reactive_conserved_to_primitive(species,tree%nodes(node)%state(:,i,j,k),t(i,j,k),primitive, &
        tree%nodes(node)%temperature(i,j,k),speed,ok)
      call require(ok,'temperature perturbation EOS')
    end do
    end do
    end do
  end do
  call average_down_reactive_eb_tree_3d(species,tree,ok)
  call require(ok,'initial average down')
  call integrate_reactive_eb_tree_3d(tree,before,ok)
  call require(ok,'initial inventory')
  call advance_reactive_eb_tree_hydro_3d(species,'hllc',1.e-8_dp,tree,boundary,wall,ok)
  call require(ok,'first hydro step')
  record%time=1.e-8_dp;record%steps=1;record%minimum_dt=1.e-8_dp;record%minimum_theta=1._dp
  record%initial=before;record%scale=abs(before);record%boundary=boundary;record%wall=wall
  record%chemical=0._dp*before
  saved=tree
  open(newunit=stream,status='scratch',form='formatted',action='readwrite')
  call reactive_eb_tree_checkpoint_record_3d(stream,.false.,species,reactions,transport,config, &
    tree,record,ok,message)
  if (.not.ok) print *,trim(message)
  call require(ok,'write tree record')
  call unchanged_tree()
  rewind(stream)
  tree%nodes(1)%geometry%halo_volume_fraction(0,0,0)=.5_dp
  restored%time=17._dp
  call reactive_eb_tree_checkpoint_record_3d(stream,.true.,species,reactions,transport,config, &
    tree,restored,ok,message)
  call require(.not.ok,'changed geometry halo rejects checkpoint')
  call unchanged_tree()
  call require(tree%nodes(1)%geometry%halo_volume_fraction(0,0,0)==.5_dp.and.restored%time==17._dp, &
    'geometry halo mismatch rollback')
  tree=saved
  rewind(stream)
  tree%nodes(1)%state=0.9_dp*tree%nodes(1)%state
  restored%time=17._dp
  call reactive_eb_tree_checkpoint_record_3d(stream,.true.,species,reactions,transport,config, &
    tree,restored,ok,message)
  if (.not.ok) print *,trim(message)
  call require(ok,'read tree record')
  do node=1,count
    call require(allocated(tree%nodes(node)%geometry%halo_volume_fraction),'restored geometry halo present')
    call require(all(tree%nodes(node)%geometry%halo_volume_fraction== &
      saved%nodes(node)%geometry%halo_volume_fraction),'exact geometry halo roundtrip')
  end do
  call unchanged_tree()
  call require(restored%time==record%time.and.restored%steps==record%steps.and. &
    all(restored%boundary==record%boundary).and.all(restored%wall==record%wall).and. &
    all(restored%chemical==record%chemical).and.all(restored%initial==record%initial).and. &
    all(restored%scale==record%scale).and.restored%minimum_dt==record%minimum_dt.and. &
    restored%minimum_theta==record%minimum_theta,'exact metadata roundtrip')
  reference=saved
  call advance_reactive_eb_tree_hydro_3d(species,'hllc',1.e-8_dp,reference,ref_b,ref_w,ok)
  call require(ok,'uninterrupted second step')
  call advance_reactive_eb_tree_hydro_3d(species,'hllc',1.e-8_dp,tree,boundary,wall,ok)
  call require(ok,'restarted second step')
  do node=1,count
    call require(all(tree%nodes(node)%state==reference%nodes(node)%state).and. &
      all(tree%nodes(node)%temperature==reference%nodes(node)%temperature),'exact restart continuation')
  end do
  call require(all(boundary==ref_b).and.all(wall==ref_w),'exact continued exchange')
  ! Overwrite only the final envelope marker; every patch parses first.
  backspace(stream)
  write(stream,'(a)') 'BROKEN_END'
  endfile(stream)
  rewind(stream)
  saved=tree;restored%time=17._dp
  call reactive_eb_tree_checkpoint_record_3d(stream,.true.,species,reactions,transport,config, &
    tree,restored,ok,message)
  call require(.not.ok,'late malformed envelope rejected')
  call unchanged_tree()
  call require(restored%time==17._dp,'late metadata rollback')
  close(stream)
  ! A wrong expected topology fails before reading any node payload.
  open(newunit=stream,status='scratch',form='formatted',action='readwrite')
  write(stream,'(a)') 'PELEF_EB_TREE_CHECKPOINT_V1'
  write(stream,*) count+1
  rewind(stream)
  call reactive_eb_tree_checkpoint_record_3d(stream,.true.,species,reactions,transport,config, &
    tree,restored,ok,message)
  call require(.not.ok,'node count mismatch rejected')
  call unchanged_tree()
  call require(restored%time==17._dp,'topology metadata rollback')
  close(stream)
  call check_direct_child_restriction()
  call check_held_restriction_policy()
  call check_sdc_record()
  call check_eos_record()
  deallocate(u,t)
  end do
  end do
  print *, 'test_eb_tree_checkpoint_3d: PASS'
contains
  subroutine check_direct_child_restriction()
    type(eb_reactive_tree_3d) :: direct,cascade,original
    real(dp) :: expected,weight
    integer :: ii,jj,kk
    logical :: valid
    original=tree
    ! Node 2 is itself covered by node 4. A saved reset of its donor energy
    ! must reach the root without first being overwritten from grandchildren.
    original%nodes(2)%state(5,3:4,3:4,2:3)=original%nodes(2)%state(5,3:4,3:4,2:3)+100._dp
    direct=original;cascade=original
    expected=0;weight=0
    do kk=1,2
    do jj=3,4
    do ii=3,4
      weight=weight+original%nodes(2)%geometry%volume_fraction(ii,jj,kk)
      expected=expected+original%nodes(2)%geometry%volume_fraction(ii,jj,kk)* &
        original%nodes(2)%state(5,ii,jj,kk)
    enddo
    enddo
    enddo
    call require(weight>0,'direct-child test fluid weight')
    expected=expected/weight
    call average_down_reactive_eb_tree_3d(species,direct,valid,direct_children=.true.)
    call require(valid,'direct-child multilevel restriction')
    call require(abs(direct%nodes(1)%state(5,2,2,2)-expected)<=2.e-12_dp*max(1._dp,abs(expected)), &
      'root uses saved immediate child energy')
    call average_down_reactive_eb_tree_3d(species,cascade,valid)
    call require(valid,'ordinary cascading multilevel restriction')
    call require(abs(cascade%nodes(1)%state(5,2,2,2)-direct%nodes(1)%state(5,2,2,2))>1._dp, &
      'fixture distinguishes saved children from cascading restriction')
  end subroutine check_direct_child_restriction

  subroutine check_held_restriction_policy()
    type(dual_energy_policy) :: policy
    real(dp) :: actual(nv,4,4,4),restricted(nv,4,4,4),held(1,4,4,4),original_held(1,4,4,4)
    integer :: ii,jj,kk
    restricted=tree%nodes(1)%state
    original_held=0
    do kk=1,4
    do jj=1,4
    do ii=1,4
      if(restricted(1,ii,jj,kk)<=0) cycle
      original_held(1,ii,jj,kk)=restricted(5,ii,jj,kk)- &
        sum(restricted(2:4,ii,jj,kk)**2)/(2*restricted(1,ii,jj,kk))
    enddo
    enddo
    enddo
    ! A covered fluid donor has positive unresolved kinetic variance after restriction.
    call require(tree%nodes(1)%geometry%volume_fraction(1,1,2)>0,'policy test covered fluid')
    ! The cold density-sheet fixture has negative formation energy. Select a
    ! positive held energy explicitly to exercise the eta2 reset branch.
    original_held(1,1,1,2)=1.e6_dp
    restricted(5,1,1,2)=original_held(1,1,1,2)+ &
      sum(restricted(2:4,1,1,2)**2)/(2*restricted(1,1,1,2))+100._dp
    actual=restricted;held=original_held
    policy%eta2=2;policy%allow_negative_energy=.false.
    call require(checkpoint_held_restriction_matches_3d(species,tree%topology,1,tree%nodes(1)%geometry, &
      actual,held,restricted,original_held,tree%nodes(1)%temperature,policy),'raw restriction accepted')
    actual(5,1,1,2)=actual(5,1,1,2)-100._dp
    call require(checkpoint_held_restriction_matches_3d(species,tree%topology,1,tree%nodes(1)%geometry, &
      actual,held,restricted,original_held,tree%nodes(1)%temperature,policy),'covered total energy reset accepted')
    policy%update_total_from_internal=.false.
    call require(.not.checkpoint_held_restriction_matches_3d(species,tree%topology,1,tree%nodes(1)%geometry, &
      actual,held,restricted,original_held,tree%nodes(1)%temperature,policy),'disabled total energy reset rejects changed energy')
    policy%update_total_from_internal=.true.
    actual(3,1,1,2)=actual(3,1,1,2)+1._dp
    call require(.not.checkpoint_held_restriction_matches_3d(species,tree%topology,1,tree%nodes(1)%geometry, &
      actual,held,restricted,original_held,tree%nodes(1)%temperature,policy),'covered momentum corruption rejected')
    actual(3,1,1,2)=restricted(3,1,1,2)
    held(1,1,1,2)=held(1,1,1,2)+1._dp
    call require(.not.checkpoint_held_restriction_matches_3d(species,tree%topology,1,tree%nodes(1)%geometry, &
      actual,held,restricted,original_held,tree%nodes(1)%temperature,policy),'covered held energy corruption rejected')
    held=original_held
    actual(5,1,1,4)=actual(5,1,1,4)+100._dp
    call require(.not.checkpoint_held_restriction_matches_3d(species,tree%topology,1,tree%nodes(1)%geometry, &
      actual,held,restricted,original_held,tree%nodes(1)%temperature,policy),'uncovered energy corruption rejected')
  end subroutine check_held_restriction_policy

  subroutine check_eos_record()
    integer :: file,mode
    type(eb_tree_checkpoint_3d) :: source_record
    real(dp) :: current(nv)
    saved=tree
    call integrate_reactive_eb_tree_3d(tree,current,ok)
    call require(ok,'EOS restart inventory')
    source_record=record;source_record%initial=current;source_record%scale=abs(current)
    source_record%boundary=0*current;source_record%wall=0*current;source_record%chemical=0*current
    do mode=0,1
      species%pinned_eos_extrapolation=mode==1;config%pinned_eos_extrapolation=mode==1
      open(newunit=file,status='scratch',form='formatted',action='readwrite')
      call reactive_eb_tree_checkpoint_record_3d(file,.false.,species,reactions,transport,config, &
        tree,source_record,ok,message)
      call require(ok,'EOS-mode tree record write')
      rewind(file)
      call reactive_eb_tree_checkpoint_record_3d(file,.true.,species,reactions,transport,config, &
        tree,restored,ok,message)
      call require(ok,'EOS-mode tree record same-mode read')
      call unchanged_tree()
      if(mode==1) then
        species%pinned_eos_temperature_min=300;species%pinned_eos_temperature_max=4000
        rewind(file);restored%time=17
        call reactive_eb_tree_checkpoint_record_3d(file,.true.,species,reactions,transport,config, &
          tree,restored,ok,message)
        call require(.not.ok.and.restored%time==17,'changed tree EOS bounds rollback')
        call unchanged_tree()
        species%pinned_eos_temperature_min=0;species%pinned_eos_temperature_max=0
      endif
      species%pinned_eos_extrapolation=mode==0;config%pinned_eos_extrapolation=mode==0
      rewind(file);restored%time=17
      call reactive_eb_tree_checkpoint_record_3d(file,.true.,species,reactions,transport,config, &
        tree,restored,ok,message)
      call require(.not.ok.and.restored%time==17,'cross-mode tree EOS transactional refusal')
      call unchanged_tree();close(file)
    enddo
    species%pinned_eos_extrapolation=.false.;config%pinned_eos_extrapolation=.false.
  end subroutine check_eos_record
  subroutine check_sdc_record()
    type(eb_tree_checkpoint_3d) :: source_record,target_record
    type(reactive_eb_3d_config) :: cfg
    real(dp) :: current(nv),bad_ir(nv-4)
    integer :: unit,p,mode
    cfg=config;cfg%time_integrator='sdc';cfg%external_forcing=[1._dp,2._dp,3._dp]
    call integrate_reactive_eb_tree_3d(tree,current,ok)
    call require(ok,'SDC checkpoint current inventory')
    source_record%initial=current;source_record%scale=abs(current)
    source_record%boundary=0*current;source_record%wall=0*current;source_record%chemical=0*current
    source_record%forcing=1.e-8_dp+0*current
    source_record%initial=current-source_record%forcing
    source_record%time=1.e-8_dp;source_record%minimum_dt=1.e-8_dp;source_record%steps=1
    do mode=1,4
      cfg%state_redist_max_order=2*(1+mod(mode-1,2))
      cfg%chemistry_enabled=mod(mode,2)==1
      if(cfg%chemistry_enabled) then
        allocate(source_record%reaction_source%fields(count))
        do p=1,count
          allocate(source_record%reaction_source%fields(p)%values(nv-4,4,4,4),source=-real(p,dp))
        end do
      else
        deallocate(source_record%reaction_source%fields)
      end if
      if(allocated(source_record%heat_release%fields)) deallocate(source_record%heat_release%fields)
      if(mode>=3) then
        allocate(source_record%heat_release%fields(count))
        do p=1,count
          allocate(source_record%heat_release%fields(p)%values(1,4,4,4),source=(-1._dp)**p*real(p,dp)*1.e9_dp)
        enddo
      endif
      open(newunit=unit,status='scratch',form='formatted',action='readwrite')
      call reactive_eb_tree_checkpoint_record_3d(unit,.false.,species,reactions,transport,cfg, &
        tree,source_record,ok,message)
      call require(ok,'write SDC tree: '//trim(message))
      rewind(unit)
      cfg%state_redist_max_order=0;target_record%time=17
      saved=tree
      call reactive_eb_tree_checkpoint_record_3d(unit,.true.,species,reactions,transport,cfg, &
        tree,target_record,ok,message)
      call require(.not.ok.and.target_record%time==17,'changed high-order checkpoint rejects transactionally')
      call unchanged_tree()
      cfg%state_redist_max_order=2*(1+mod(mode-1,2))
      rewind(unit)
      call reactive_eb_tree_checkpoint_record_3d(unit,.true.,species,reactions,transport,cfg, &
        tree,target_record,ok,message)
      call require(ok,'read SDC tree: '//trim(message))
      call require(all(target_record%forcing==source_record%forcing),'SDC forcing exact round trip')
      call require(allocated(target_record%heat_release%fields).eqv.(mode>=3),'heat presence round trip')
      if(mode>=3) then
        do p=1,count
          call require(all(target_record%heat_release%fields(p)%values== &
            source_record%heat_release%fields(p)%values),'signed stage heat exact round trip')
        enddo
      endif
      if(cfg%chemistry_enabled) then
        do p=1,count
          call require(all(target_record%reaction_source%fields(p)%values== &
            source_record%reaction_source%fields(p)%values),'SDC reaction exact round trip')
        end do
        backspace(unit);backspace(unit)
        bad_ir=0;bad_ir(1)=ieee_value(0._dp,ieee_quiet_nan)
        write(unit,'(*(es27.18e3,1x))') bad_ir
        write(unit,'(a)') 'END_EB_TREE_CHECKPOINT'
        endfile(unit);rewind(unit)
        target_record%forcing=17;target_record%reaction_source%fields(1)%values=19
        if(mode>=3) target_record%heat_release%fields(1)%values=23
        saved=tree
        call reactive_eb_tree_checkpoint_record_3d(unit,.true.,species,reactions,transport,cfg, &
          tree,target_record,ok,message)
        call require(.not.ok,'late SDC diagnostic NaN rejects')
        if(mode>=3) call require(all(target_record%heat_release%fields(1)%values==23),'late heat metadata rollback')
        call unchanged_tree()
        call require(all(target_record%forcing==17).and.all(target_record%reaction_source%fields(1)%values==19), &
          'SDC source metadata rollback')
      else
        call require(.not.allocated(target_record%reaction_source%fields),'nonreacting SDC has no source array')
      end if
      rewind(unit);cfg%external_forcing(1)=4
      call reactive_eb_tree_checkpoint_record_3d(unit,.true.,species,reactions,transport,cfg, &
        tree,target_record,ok,message)
      call require(.not.ok,'SDC forcing configuration mismatch')
      cfg%external_forcing(1)=1
      close(unit)
    end do
  end subroutine check_sdc_record

  subroutine unchanged_tree()
    integer :: n
    do n=1,size(tree%nodes)
      call require(all(tree%nodes(n)%state==saved%nodes(n)%state).and. &
        all(tree%nodes(n)%temperature==saved%nodes(n)%temperature),'all state and temperature preservation')
    end do
  end subroutine unchanged_tree
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine require
end program test_eb_tree_checkpoint_3d
