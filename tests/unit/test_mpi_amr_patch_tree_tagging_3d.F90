program test_mpi_amr_patch_tree_tagging_3d
  use mpi_f08
  use mpi_amr_patch_tree_tagging_3d_mod, only: tag_mpi_amr_tree_temperature_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
  use amr_patch_tree_tagging_3d_mod, only: serial_tag=>tag_amr_patch_tree_temperature_3d
  use amr_patch_tree_remap_3d_mod, only: remap_amr_patch_tree_3d
  implicit none
  type(amr_patch_tree_3d) :: root,tree,old,saved
  type(amr_tree_state_3d) :: temperatures,source,destination
  type(amr_patch_3d) :: empty(0),patches(1)
  real(dp) :: uniform(1,8,4,4),small(1,4,4,4),before(1),after(1),nan
  logical :: ok
  integer :: node,rank,nranks,ierr
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks>=2,"rank count")
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
  if (rank==0) print *, 'PASS: MPI tagging matches serial, periodic/fine-only tags, empty ranks and collective rollback'
  call MPI_Finalize(ierr)
contains
  subroutine tag_amr_patch_tree_temperature_3d(old_tree,temperature,thresholds,ratios,tile_shape,buffer, &
      maximum_nodes,periodic,new_tree,ok)
    type(amr_patch_tree_3d), intent(in) :: old_tree
    type(amr_tree_state_3d), intent(in) :: temperature
    real(dp), intent(in) :: thresholds(:)
    integer, intent(in) :: ratios(:),tile_shape(3),buffer,maximum_nodes
    logical, intent(in) :: periodic(3)
    type(amr_patch_tree_3d), intent(inout) :: new_tree
    logical, intent(out) :: ok
    type(amr_patch_tree_3d) :: oracle,candidate,original,retained
    type(amr_tree_state_3d) :: local
    integer :: owners(size(old_tree%nodes)),peer,mode
    real(dp) :: changed(size(thresholds))
    logical :: valid,expected
    original=new_tree;oracle=new_tree
    call serial_tag(old_tree,temperature,thresholds,ratios,tile_shape,buffer,maximum_nodes,periodic,oracle,expected)
    do mode=1,2
      owners=[(mod(peer-1,nranks),peer=1,size(old_tree%nodes))]
      if (mode==2) owners=nranks-1
      if (allocated(local%fields)) deallocate(local%fields)
      allocate(local%fields(size(owners)))
      do peer=1,size(owners)
        if (owners(peer)==rank) local%fields(peer)=temperature%fields(peer)
      end do
      candidate=original
      call tag_mpi_amr_tree_temperature_3d(old_tree,owners,local,thresholds,ratios,tile_shape,buffer, &
        maximum_nodes,periodic,candidate,MPI_COMM_WORLD,valid)
      call require(valid.eqv.expected,'MPI/serial success agreement')
      call same_tree(candidate,oracle)
      if (valid.and.size(thresholds)>0) then
        retained=candidate;changed=thresholds
        if (rank==0) changed(1)=changed(1)+1
        call tag_mpi_amr_tree_temperature_3d(old_tree,owners,local,changed,ratios,tile_shape,buffer, &
          maximum_nodes,periodic,candidate,MPI_COMM_WORLD,valid)
        call require(.not.valid,'different thresholds rejected collectively')
        call same_tree(candidate,retained)
        if (rank==owners(1)) local%fields(1)%values(1,1,1,1)=ieee_value(0.0_dp,ieee_quiet_nan)
        call tag_mpi_amr_tree_temperature_3d(old_tree,owners,local,thresholds,ratios,tile_shape,buffer, &
          maximum_nodes,periodic,candidate,MPI_COMM_WORLD,valid)
        call require(.not.valid,'one-owner NaN rejected collectively')
        call same_tree(candidate,retained)
      end if
    end do
    new_tree=oracle;ok=expected
  end subroutine
  subroutine same_tree(actual,expected)
    type(amr_patch_tree_3d), intent(in) :: actual,expected
    integer :: peer
    call require(allocated(actual%nodes).eqv.allocated(expected%nodes),'tree allocation')
    if (.not.allocated(expected%nodes)) return
    call require(size(actual%nodes)==size(expected%nodes),'tree node count')
    do peer=1,size(expected%nodes)
      call require(actual%nodes(peer)%parent==expected%nodes(peer)%parent,'tree parent')
      call require(actual%nodes(peer)%level==expected%nodes(peer)%level,'tree level')
      call require(all(actual%nodes(peer)%shape==expected%nodes(peer)%shape),'tree shape')
      if (peer==1) cycle
      associate(a=>actual%nodes(peer)%patch,b=>expected%nodes(peer)%patch)
        call require(all([a%coarse_i_lower,a%coarse_i_upper,a%coarse_j_lower,a%coarse_j_upper, &
          a%coarse_k_lower,a%coarse_k_upper,a%refinement_ratio]== &
          [b%coarse_i_lower,b%coarse_i_upper,b%coarse_j_lower,b%coarse_j_upper, &
          b%coarse_k_lower,b%coarse_k_upper,b%refinement_ratio]),'tree patch')
      end associate
    end do
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: status
    if (.not.condition) then
      print *, 'FAIL: '//label
      call MPI_Abort(MPI_COMM_WORLD,1,status)
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
