program test_mpi_amr_patch_tree_restriction_3d
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_is_nan
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d,average_down_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    average_down_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
  use mpi_amr_patch_tree_restriction_3d_mod, only: average_down_mpi_amr_tree_3d
  implicit none
  type(amr_patch_tree_3d) :: tree,changed
  type(amr_patch_3d) :: patches(3)
  type(amr_tree_state_3d) :: global,local,oracle,saved,broken
  real(dp) :: before(2),after(2)
  integer :: rank,nranks,ierr,node,parent,n(3),i,j,k,mode,depth,owners(4),other(4),bad_depth
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks>=2,'rank count')
  call initialize_amr_patch_3d(4,4,4,1,2,1,2,1,2,2,patches(1),ok)
  call require(ok,'child')
  call initialize_amr_patch_3d(4,4,4,3,4,1,2,1,2,2,patches(2),ok)
  call require(ok,'peer')
  call initialize_amr_patch_3d(4,4,4,2,2,2,2,2,2,3,patches(3),ok)
  call require(ok,'grandchild')
  call initialize_amr_patch_tree_3d([4,4,4],[1,1,2],patches,tree,ok)
  call require(ok,'tree')
  allocate(global%fields(4))
  do node=1,4
    n=tree%nodes(node)%shape
    allocate(global%fields(node)%values(2,n(1),n(2),n(3)))
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          global%fields(node)%values(:,i,j,k)=[2.0_dp,3.0_dp]*node+.01_dp*i+.02_dp*j+.03_dp*k
        end do
      end do
    end do
  end do
  call composite_integrals_amr_patch_tree_3d(tree,global,[.1_dp,.2_dp,.3_dp],before,ok)
  call require(ok,'initial composite inventory')
  do mode=1,3
    owners=[0,1,1,0]
    if (mode==2) owners=0
    if (mode==3) owners=mod([0,1,2,0],nranks)
    do depth=0,2
      oracle=global
      if (depth==0) then
        call average_down_amr_patch_tree_3d(tree,oracle,ok)
        call require(ok,'full serial restriction')
      else
        do node=4,2,-1
          parent=tree%nodes(node)%parent
          if (tree%nodes(parent)%level<depth) cycle
          call average_down_3d(oracle%fields(parent)%values,oracle%fields(node)%values,tree%nodes(node)%patch,ok)
          call require(ok,'interval serial restriction')
        end do
      end if
      call composite_integrals_amr_patch_tree_3d(tree,oracle,[.1_dp,.2_dp,.3_dp],after,ok)
      call require(ok.and.all(after==before),'restriction preserves exposed inventory')
      allocate(local%fields(4))
      do node=1,4
        if (owners(node)==rank) local%fields(node)=global%fields(node)
      end do
      if (depth==0) then
        call average_down_mpi_amr_tree_3d(tree,owners,local,2,MPI_COMM_WORLD,ok)
      else
        call average_down_mpi_amr_tree_3d(tree,owners,local,2,MPI_COMM_WORLD,ok,minimum_depth=depth)
      end if
      call require(ok,'collective restriction')
      do node=1,4
        call require(allocated(local%fields(node)%values).eqv.(owners(node)==rank),'only owned fields')
        if (owners(node)/=rank) cycle
        call require(all(local%fields(node)%values==oracle%fields(node)%values),'exact serial parity')
      end do
      saved=local
      bad_depth=0
      if (rank==nranks-1) bad_depth=1
      call average_down_mpi_amr_tree_3d(tree,owners,local,2,MPI_COMM_WORLD,ok,minimum_depth=bad_depth)
      call require(.not.ok,'different valid depths rejected')
      call unchanged()
      call average_down_mpi_amr_tree_3d(tree,owners,local,2,MPI_COMM_WORLD,ok,minimum_depth=3)
      call require(.not.ok,'out of range depth rejected')
      call unchanged()
      other=owners
      if (rank==0) other(1)=1
      call average_down_mpi_amr_tree_3d(tree,other,local,2,MPI_COMM_WORLD,ok)
      call require(.not.ok,'different ownership rejected')
      call unchanged()
      changed=tree
      if (rank==0) then
        changed%nodes(4)%patch%coarse_i_lower=1
        changed%nodes(4)%patch%coarse_i_upper=1
      end if
      call average_down_mpi_amr_tree_3d(changed,owners,local,2,MPI_COMM_WORLD,ok)
      call require(.not.ok,'different valid topology rejected')
      call unchanged()
      broken=local
      if (rank==owners(4)) deallocate(broken%fields(4)%values)
      call average_down_mpi_amr_tree_3d(tree,owners,broken,2,MPI_COMM_WORLD,ok)
      call require(.not.ok,'missing child field rejected')
      if (rank==owners(4)) call require(.not.allocated(broken%fields(4)%values),'missing allocation unchanged')
      if (rank==owners(4)) local%fields(4)%values(1,1,1,1)=ieee_value(0.0_dp,ieee_quiet_nan)
      call average_down_mpi_amr_tree_3d(tree,owners,local,2,MPI_COMM_WORLD,ok)
      call require(.not.ok,'nonfinite child rejected collectively')
      if (rank==owners(4)) then
        call require(ieee_is_nan(local%fields(4)%values(1,1,1,1)),'invalid source preserved')
        local%fields(4)%values(1,1,1,1)=saved%fields(4)%values(1,1,1,1)
      end if
      call unchanged()
      deallocate(local%fields)
    end do
  end do
  if (rank==0) print *, 'PASS: distributed bottom-up restriction, depth limits, serial parity and rollback'
  call MPI_Finalize(ierr)
contains
  subroutine unchanged()
    integer :: peer
    do peer=1,4
      call require(allocated(local%fields(peer)%values).eqv.allocated(saved%fields(peer)%values),'rollback allocation')
      if (.not.allocated(saved%fields(peer)%values)) cycle
      call require(all(local%fields(peer)%values==saved%fields(peer)%values),'rollback values')
    end do
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: status
    if (condition) return
    print *, 'FAIL rank ',rank,': '//label
    call MPI_Abort(MPI_COMM_WORLD,1,status)
  end subroutine
end program
