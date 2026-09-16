program test_amr_patch_tree_tagging_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
  use amr_patch_tree_tagging_3d_mod, only: tag_amr_patch_tree_temperature_3d
  use amr_patch_tree_remap_3d_mod, only: remap_amr_patch_tree_3d
  implicit none
  type(amr_patch_tree_3d) :: root,tree,old,saved
  type(amr_tree_state_3d) :: temperatures,source,destination
  type(amr_patch_3d) :: empty(0),patches(1)
  real(dp) :: uniform(1,8,4,4),small(1,4,4,4),before(1),after(1),nan
  logical :: ok
  integer :: node
  call initialize_amr_patch_tree_3d([8,4,4],[integer::],empty,root,ok)
  call require(ok,'root geometry')
  uniform=1000
  call prolong_amr_patch_tree_3d(root,uniform,temperatures,ok)
  call require(ok,'root temperature')
  temperatures%fields(1)%values(1,1,1,1)=2000
  temperatures%fields(1)%values(1,8,4,4)=2200
  call tag_amr_patch_tree_temperature_3d(root,temperatures,[1500._dp,2100._dp],[2,3],[2,2,2],0,32, &
    [.false.,.false.,.false.],tree,ok)
  call require(ok,'disjoint multilevel tagging')
  call require(size(tree%nodes)==4.and.tree%level_count()==3,'two islands and one nested hotter island')
  call require(tree%nodes(2)%parent==1.and.tree%nodes(3)%parent==1.and.tree%nodes(4)%parent==3,'parent order')
  call require(tree%nodes(2)%patch%coarse_i_lower==1.and.tree%nodes(3)%patch%coarse_i_lower==7,'island bounds')
  call require(tree%nodes(4)%patch%coarse_i_lower==3.and.all(tree%nodes(4)%shape==6),'ratio-3 nested bounds')
  source=temperatures
  call composite_integrals_amr_patch_tree_3d(root,source,[1._dp,1._dp,1._dp],before,ok)
  call require(ok,'original inventory')
  call remap_amr_patch_tree_3d(root,source,tree,destination,ok)
  call require(ok,'generated geometry remap')
  call composite_integrals_amr_patch_tree_3d(tree,destination,[1._dp,1._dp,1._dp],after,ok)
  call require(ok.and.maxval(abs(after-before))<1e-9_dp,'generated geometry conservation')

  temperatures%fields(1)%values=1000
  temperatures%fields(1)%values(1,1,1,1)=2000
  call tag_amr_patch_tree_temperature_3d(root,temperatures,[1500._dp],[2],[2,2,2],1,32, &
    [.false.,.false.,.false.],tree,ok)
  call require(ok.and.size(tree%nodes)==2,'physical boundary clips buffer')
  call tag_amr_patch_tree_temperature_3d(root,temperatures,[1500._dp],[2],[2,2,2],1,32, &
    [.true.,.true.,.true.],tree,ok)
  call require(ok.and.size(tree%nodes)==9,'periodic corner buffer selects eight tiles')
  call require(tree%nodes(9)%patch%coarse_i_lower==7,'periodic opposite corner')

  ! A fine-only hot cell survives a ratio change even though its parent is cold.
  call initialize_amr_patch_3d(4,4,4,1,1,1,1,1,1,2,patches(1),ok)
  call require(ok,'old patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1],patches,old,ok)
  call require(ok,'old tree')
  small=1000
  call prolong_amr_patch_tree_3d(old,small,temperatures,ok)
  call require(ok,'old temperature')
  temperatures%fields(2)%values(1,2,2,2)=2000
  call tag_amr_patch_tree_temperature_3d(old,temperatures,[1500._dp],[3],[1,1,1],0,32, &
    [.false.,.false.,.false.],tree,ok)
  call require(ok.and.size(tree%nodes)==2,'fine-only tag')
  call require(all(tree%nodes(2)%shape==3),'changed refinement ratio')
  temperatures%fields(1)%values(1,1,1,1)=3000
  temperatures%fields(2)%values=1000
  call tag_amr_patch_tree_temperature_3d(old,temperatures,[1500._dp],[3],[1,1,1],0,32, &
    [.false.,.false.,.false.],tree,ok)
  call require(ok.and.size(tree%nodes)==1,'covered hot parent ignored; full coarsening')

  ! Transactional failures preserve the destination tree.
  saved=old;tree=saved
  temperatures%fields(2)%values=2000
  call tag_amr_patch_tree_temperature_3d(old,temperatures,[1500._dp],[3],[1,1,1],0,1, &
    [.false.,.false.,.false.],tree,ok)
  call require(.not.ok,'node budget rejection')
  call unchanged()
  call tag_amr_patch_tree_temperature_3d(old,temperatures,[1500._dp,1400._dp],[2,2],[1,1,1],0,32, &
    [.false.,.false.,.false.],tree,ok)
  call require(.not.ok,'decreasing thresholds rejection')
  call unchanged()
  call tag_amr_patch_tree_temperature_3d(old,temperatures,[1500._dp,1500._dp], &
    [huge(1),huge(1)],[1,1,1],0,32,[.false.,.false.,.false.],tree,ok)
  call require(.not.ok,'integer lattice overflow rejection')
  call unchanged()
  nan=ieee_value(1._dp,ieee_quiet_nan)
  call tag_amr_patch_tree_temperature_3d(old,temperatures,[nan],[2],[1,1,1],0,32, &
    [.false.,.false.,.false.],tree,ok)
  call require(.not.ok,'nonfinite threshold rejection')
  call unchanged()
  call tag_amr_patch_tree_temperature_3d(old,temperatures,[real(dp)::],[integer::],[1,1,1],0,1, &
    [.false.,.false.,.false.],tree,ok)
  call require(ok.and.size(tree%nodes)==1,'zero refinement levels')
  print *, 'PASS: tree temperature tagging, multilevel tiles, periodic buffers, fine-only tags, coarsening and rollback'
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: '//label
      error stop 1
    end if
  end subroutine
  subroutine unchanged()
    call require(size(tree%nodes)==size(saved%nodes),'rollback count')
    do node=1,size(saved%nodes)
      call require(tree%nodes(node)%parent==saved%nodes(node)%parent,'rollback parent')
      call require(all(tree%nodes(node)%shape==saved%nodes(node)%shape),'rollback shape')
      call require(tree%nodes(node)%patch%refinement_ratio==saved%nodes(node)%patch%refinement_ratio,'rollback ratio')
    end do
  end subroutine
end program
