program test_eb_patch_tree_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,react_reactive_eb_tree_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d,reactive_eb_integrals_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp), allocatable :: chemical(:)
  type(amr_patch_3d) :: boxes(3)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,middle,coarse
  type(eb_reactive_tree_3d) :: tree,saved
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:),initial(:),before(:),after(:),grandchild(:)
  integer :: nv,axis,node,ratio
  logical :: ok
  character(len=1), parameter :: axes(3)=['x','y','z']
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_mechanism(reactions,ok)
  call require(ok,'reactions')
  nv=reactive_nvar(size(species))
  allocate(u(nv,4,4,4),t(4,4,4),initial(nv),before(nv),after(nv),grandchild(nv),chemical(nv))
  do ratio=2,3
  call initialize_amr_patch_3d(4,4,4,1,2,1,2,2,3,ratio,boxes(1),ok)
  call require(ok,'first child')
  call initialize_amr_patch_3d(4,4,4,3,4,3,4,2,3,ratio,boxes(2),ok)
  call require(ok,'second child')
  call initialize_amr_patch_3d(2*ratio,2*ratio,2*ratio,ratio,ratio+1,ratio,ratio+1,ratio,ratio+1,2,boxes(3),ok)
  call require(ok,'grandchild')
  call initialize_amr_patch_tree_3d([4,4,4],[1,1,2],boxes,topology,ok)
  call require(ok,'three-level topology')
  do axis=1,3
    call build_axis_plane_eb_geometry_3d(8*ratio,8*ratio,8*ratio,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
      axes(axis),.31_dp,finest,ok)
    call require(ok,'finest geometry')
    call coarsen_eb_geometry_3d(finest,2,middle,ok)
    call require(ok,'middle geometry')
    call coarsen_eb_geometry_3d(middle,ratio,coarse,ok)
    call require(ok,'coarse geometry')
    config%nx=4;config%ny=4;config%nz=4
    call initialize_reactive_eb_density_sheet_3d(species,config,coarse,u,t,ok)
    call require(ok,'root state')
    call reactive_eb_integrals_3d(u,coarse,initial,ok)
    call require(ok,'root inventory')
    call initialize_reactive_eb_tree_3d(species,topology,finest,u,t,tree,ok)
    call require(ok,'EB tree initialization')
    call integrate_reactive_eb_tree_3d(tree,before,ok)
    call require(ok.and.all(abs(before-initial)<3.e-14_dp*max(1._dp,abs(initial))),'PCM composite conservation')
    call reactive_eb_integrals_3d(tree%nodes(4)%state,tree%nodes(4)%geometry,grandchild,ok)
    call require(ok,'grandchild inventory')
    tree%nodes(4)%state=1.25_dp*tree%nodes(4)%state
    ! Hidden parent values must neither contribute nor block bottom-up repair.
    tree%nodes(1)%state(:,1:2,1:2,2:3)=-99._dp;tree%nodes(1)%temperature(1:2,1:2,2:3)=-99._dp
    tree%nodes(2)%state(:,ratio:ratio+1,ratio:ratio+1,ratio:ratio+1)=-99._dp
    tree%nodes(2)%temperature(ratio:ratio+1,ratio:ratio+1,ratio:ratio+1)=-99._dp
    call integrate_reactive_eb_tree_3d(tree,before,ok)
    call require(ok,'exposed inventory despite stale parents')
    call require(all(abs(before-initial-.25_dp*grandchild)<3.e-14_dp*max(1._dp,abs(initial))), &
      'grandchild contributes exactly once')
    call average_down_reactive_eb_tree_3d(species,tree,ok)
    call require(ok,'bottom-up restriction')
    call integrate_reactive_eb_tree_3d(tree,after,ok)
    call require(ok.and.all(after==before),'restriction preserves composite inventory exactly')
    call reactive_eb_integrals_3d(tree%nodes(1)%state,tree%nodes(1)%geometry,after,ok)
    call require(ok.and.all(abs(after-before)<3.e-14_dp*max(1._dp,abs(before))),'root receives deep inventories')
    ! Hidden parents are deliberately invalid when chemistry starts. Only
    ! exposed cells may react; subsequent restriction repairs hidden values.
    tree%nodes(1)%state(:,1:2,1:2,2:3)=-99._dp;tree%nodes(1)%temperature(1:2,1:2,2:3)=-99._dp
    tree%nodes(2)%state(:,ratio:ratio+1,ratio:ratio+1,ratio:ratio+1)=-99._dp
    tree%nodes(2)%temperature(ratio:ratio+1,ratio:ratio+1,ratio:ratio+1)=-99._dp
    call react_reactive_eb_tree_3d(species,reactions,1.e-9_dp,2.e-7_dp,1.e-12_dp,tree,chemical,ok,'explicit')
    call require(ok,'exposed-cell tree chemistry')
    call require(any(abs(chemical(6:))>1.e-15_dp),'nonzero tree species reaction')
    call integrate_reactive_eb_tree_3d(tree,after,ok)
    call require(ok.and.all(after-before==chemical),'chemical inventory accounting')
    call require(all(abs(chemical(1:5))<3.e-14_dp*max(1._dp,abs(before(1:5)))),'reaction invariants')
    tree%nodes(4)%state(1,4,4,4)=-1._dp
    saved=tree
    call average_down_reactive_eb_tree_3d(species,tree,ok)
    call require(.not.ok,'invalid grandchild rejection')
    chemical=17._dp
    call react_reactive_eb_tree_3d(species,reactions,1.e-9_dp,2.e-7_dp,1.e-12_dp,tree,chemical,ok,'explicit')
    call require(.not.ok.and.all(chemical==17._dp),'late-node reaction rejection and ledger rollback')
    do node=1,4
      call require(all(tree%nodes(node)%state==saved%nodes(node)%state).and. &
        all(tree%nodes(node)%temperature==saved%nodes(node)%temperature),'whole-tree rollback')
    end do
    tree=saved
    tree%nodes(4)%geometry%x_lower=tree%nodes(4)%geometry%x_lower+.01_dp
    after=17._dp
    call integrate_reactive_eb_tree_3d(tree,after,ok)
    call require(.not.ok.and.all(after==17._dp),'geometry mismatch integral rollback')
    tree=saved;u(1,4,4,4)=-1._dp
    call initialize_reactive_eb_tree_3d(species,topology,finest,u,t,tree,ok)
    call require(.not.ok,'invalid root rejection')
    call require(all(tree%nodes(1)%state==saved%nodes(1)%state),'initialization rollback')
  end do
  end do
  print *, 'test_eb_patch_tree_3d: PASS'
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine require
end program test_eb_patch_tree_3d
