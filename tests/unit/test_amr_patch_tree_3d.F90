program test_amr_patch_tree_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,average_down_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
  implicit none
  type(amr_patch_tree_3d) :: tree,saved_tree,root_only
  type(amr_tree_state_3d) :: state,saved_state,single
  type(amr_patch_3d) :: patches(5),bad(5),empty(0)
  integer :: parents(5),bad_parents(5),i,j,k,node
  real(dp) :: root(2,8,6,4),spacing(3),reference(2),integrals(2),after(2),leaf_volume
  logical :: ok
  parents=[1,1,2,3,4];spacing=[.5_dp,.25_dp,.125_dp]
  ! Two adjacent root children, a child on each branch, then a fourth level.
  call initialize_amr_patch_3d(8,6,4,1,2,2,3,1,2,2,patches(1),ok)
  call require(ok,'first child')
  call initialize_amr_patch_3d(8,6,4,3,4,2,3,1,2,2,patches(2),ok)
  call require(ok,'adjacent sibling')
  call initialize_amr_patch_3d(4,4,4,1,2,1,2,1,2,3,patches(3),ok)
  call require(ok,'second-level child')
  call initialize_amr_patch_3d(4,4,4,2,3,2,3,2,3,3,patches(4),ok)
  call require(ok,'second branch')
  call initialize_amr_patch_3d(6,6,6,1,2,2,3,3,4,2,patches(5),ok)
  call require(ok,'third-level child')
  call initialize_amr_patch_tree_3d([8,6,4],parents,patches,tree,ok)
  call require(ok.and.tree%is_valid(),'tree topology')
  call require(tree%level_count()==4.and.size(tree%nodes)==6,'four levels and six patches')
  do k=1,4
    do j=1,6
      do i=1,8
        root(:,i,j,k)=[real(i+10*j+100*k,dp),real(2*i+j+3*k,dp)]
      end do
    end do
  end do
  reference=[sum(root(1,:,:,:)),sum(root(2,:,:,:))]*product(spacing)
  call prolong_amr_patch_tree_3d(tree,root,state,ok)
  call require(ok.and.state%is_valid(tree),'PCM initialization throughout tree')
  integrals=-1
  call composite_integrals_amr_patch_tree_3d(tree,state,spacing,integrals,ok)
  if (maxval(abs(integrals-reference))>=1.e-11_dp) print *, 'PCM inventory difference: ',integrals-reference
  call require(ok.and.maxval(abs(integrals-reference))<1.e-11_dp,'PCM composite conservation without double counting')
  call require(all(state%fields(6)%values(1,:,:,:)==root(1,1,2,1)),'deep child samples correct ancestor')
  leaf_volume=product(spacing)/8.0_dp/27.0_dp/8.0_dp
  state%fields(6)%values(1,:,:,:)=state%fields(6)%values(1,:,:,:)+3.0_dp
  state%fields(6)%values(2,:,:,:)=state%fields(6)%values(2,:,:,:)-2.0_dp
  reference=reference+[3.0_dp,-2.0_dp]*real(product(tree%nodes(6)%shape),dp)*leaf_volume
  call composite_integrals_amr_patch_tree_3d(tree,state,spacing,integrals,ok)
  call require(ok.and.maxval(abs(integrals-reference))<1.e-11_dp,'deep leaf increment has its physical volume')
  call average_down_amr_patch_tree_3d(tree,state,ok)
  call require(ok,'bottom-up restriction')
  call composite_integrals_amr_patch_tree_3d(tree,state,spacing,after,ok)
  call require(ok.and.maxval(abs(after-integrals))<1.e-11_dp,'restriction preserves leaf composite inventory')
  call require(maxval(abs([sum(state%fields(1)%values(1,:,:,:)), &
    sum(state%fields(1)%values(2,:,:,:))]*product(spacing)-reference))<1.e-11_dp,'restriction reaches root')
  saved_tree=tree;saved_state=state
  bad=patches;bad(2)=bad(1)
  call initialize_amr_patch_tree_3d([8,6,4],parents,bad,tree,ok)
  call rejected_topology('overlapping siblings')
  bad_parents=parents;bad_parents(3)=4
  call initialize_amr_patch_tree_3d([8,6,4],bad_parents,patches,tree,ok)
  call rejected_topology('cyclic parent')
  bad=patches;bad(4)%refinement_ratio=2
  call initialize_amr_patch_tree_3d([8,6,4],parents,bad,tree,ok)
  call rejected_topology('inconsistent same-level refinement')
  bad=patches;bad(5)%coarse_nx=7
  call initialize_amr_patch_tree_3d([8,6,4],parents,bad,tree,ok)
  call rejected_topology('wrong parent shape')
  root(1,1,1,1)=ieee_value(0.0_dp,ieee_quiet_nan)
  call prolong_amr_patch_tree_3d(tree,root,state,ok)
  call require(.not.ok,'invalid root rejected')
  do node=1,size(state%fields)
    call require(all(state%fields(node)%values==saved_state%fields(node)%values),'failed prolongation rollback')
  end do
  state%fields(6)%values(1,1,1,1)=ieee_value(0.0_dp,ieee_quiet_nan)
  call average_down_amr_patch_tree_3d(tree,state,ok)
  call require(.not.ok.and.all(state%fields(1)%values==saved_state%fields(1)%values),'invalid descendant leaves root intact')
  after=-17
  call composite_integrals_amr_patch_tree_3d(tree,state,spacing,after,ok)
  call require(.not.ok.and.all(after==-17),'invalid field preserves output inventory')
  call initialize_amr_patch_tree_3d([8,6,4],[integer::],empty,root_only,ok)
  call require(ok.and.root_only%level_count()==1,'unrefined tree')
  call prolong_amr_patch_tree_3d(root_only,saved_state%fields(1)%values,single,ok)
  call require(ok,'unrefined initialization')
  call composite_integrals_amr_patch_tree_3d(root_only,single,spacing,after,ok)
  call require(ok.and.maxval(abs(after-reference))<1.e-11_dp,'unrefined inventory')
  print '(a)','test_amr_patch_tree_3d: PASS'
contains
  subroutine rejected_topology(label)
    character(len=*), intent(in) :: label
    call require(.not.ok,'reject '//label)
    call require(tree%is_valid().and.all(tree%nodes%parent==saved_tree%nodes%parent),'topology rollback '//label)
    do node=1,size(tree%nodes)
      call require(all(tree%nodes(node)%shape==saved_tree%nodes(node)%shape),'shape rollback '//label)
    end do
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (condition) return
    print '(a)','FAIL: '//label
    error stop 1
  end subroutine
end program
