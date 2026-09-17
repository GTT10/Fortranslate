program test_mpi_amr_patch_tree_distribution_3d
  use mpi_f08
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use mpi_amr_patch_tree_distribution_3d_mod, only: &
    scatter_mpi_amr_tree_3d,gather_mpi_amr_tree_3d,partition_amr_tree_3d,redistribute_mpi_amr_tree_3d, &
    rebalance_mpi_amr_tree_3d
  implicit none
  type(amr_patch_tree_3d) :: tree,different_tree
  type(amr_patch_3d) :: patches(2)
  type(amr_tree_state_3d) :: global,local,result,saved,broken
  integer :: ierr,rank,nranks,root,node,n(3),owners(3),other_owners(3),mode,i,j,k
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks>=2,'two or more ranks')
  root=nranks-1
  call initialize_amr_patch_3d(4,4,4,1,2,1,2,1,2,2,patches(1),ok)
  call require(ok,'child geometry')
  call initialize_amr_patch_3d(4,4,4,1,1,1,1,1,1,3,patches(2),ok)
  call require(ok,'grandchild geometry')
  call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,tree,ok)
  call require(ok,'tree')
  block
    integer, allocatable :: balanced(:),previous(:)
    type(amr_patch_tree_3d) :: oversized
    call partition_amr_tree_3d(tree,2,0,balanced,ok)
    call require(ok,'common partition')
    call require(all(balanced==[0,1,0]),'stable cell-weight ties')
    call partition_amr_tree_3d(tree,2,2,balanced,ok)
    call require(ok,'recursive transport partition')
    call require(all(balanced==[1,1,0]),'largest subcycled workload first')
    call partition_amr_tree_3d(tree,4,1,balanced,ok)
    call require(ok,'empty-rank partition')
    call require(all(balanced==[2,1,0]),'hydro weights and rank ties')
    previous=balanced
    call partition_amr_tree_3d(tree,0,0,balanced,ok)
    call require(.not.ok,'invalid rank count rejected')
    call require(all(balanced==previous),'partition rollback')
    call initialize_amr_patch_tree_3d([huge(1),huge(1),huge(1)],[integer::],patches(:0),oversized,ok)
    call require(ok,'metadata-only oversized tree')
    call partition_amr_tree_3d(oversized,2,0,balanced,ok)
    call require(.not.ok,'cell-work overflow rejected')
    call require(all(balanced==previous),'overflow rollback')
  end block
  if (rank==root) then
    allocate(global%fields(3))
    do node=1,3
      n=tree%nodes(node)%shape
      allocate(global%fields(node)%values(2,n(1),n(2),n(3)))
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            global%fields(node)%values(:,i,j,k)=values(node,i,j,k)
          end do
        end do
      end do
    end do
  end if
  do mode=1,2
    owners=[0,1,0]
    if (mode==2) owners=0
    call scatter_mpi_amr_tree_3d(tree,owners,global,local,MPI_COMM_WORLD,root,ok)
    call require(ok,'scatter')
    do node=1,3
      call require(allocated(local%fields(node)%values).eqv.(owners(node)==rank),'only owner stores payload')
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            call require(all(local%fields(node)%values(:,i,j,k)==values(node,i,j,k)),'received component/cell order')
          end do
        end do
      end do
      local%fields(node)%values=local%fields(node)%values+10*rank
    end do
    call gather_mpi_amr_tree_3d(tree,owners,local,result,MPI_COMM_WORLD,root,ok)
    call require(ok,'gather')
    if (rank==root) then
      do node=1,3
        call require(all(result%fields(node)%values==global%fields(node)%values+10*owners(node)),'exact gathered values')
      end do
    else
      call require(.not.allocated(result%fields),'nonroot gather output untouched')
    end if
    saved=local;different_tree=tree
    if (rank==0) then
      different_tree%nodes(3)%patch%coarse_i_lower=2
      different_tree%nodes(3)%patch%coarse_i_upper=2
    end if
    call scatter_mpi_amr_tree_3d(different_tree,owners,global,local,MPI_COMM_WORLD,root,ok)
    call require(.not.ok,'different valid topologies rejected collectively')
    call unchanged()
    other_owners=owners
    if (rank==0) other_owners(1)=1
    call scatter_mpi_amr_tree_3d(tree,other_owners,global,local,MPI_COMM_WORLD,root,ok)
    call require(.not.ok,'inconsistent owners rejected collectively')
    call unchanged()
    call scatter_mpi_amr_tree_3d(tree,owners,global,local,MPI_COMM_WORLD,rank,ok)
    call require(.not.ok,'inconsistent root rejected collectively')
    call unchanged()
    broken=local
    if (rank==0) deallocate(broken%fields(1)%values)
    call gather_mpi_amr_tree_3d(tree,owners,broken,result,MPI_COMM_WORLD,root,ok)
    call require(.not.ok,'missing owner data rejected collectively')
    if (rank==root) then
      do node=1,3
        call require(all(result%fields(node)%values==global%fields(node)%values+10*owners(node)),'failed gather preserves output')
      end do
    end if
    block
      type(amr_tree_state_3d) :: moving,moving_t,backup,backup_t
      integer :: destinations(3),invalid(3),peer
      real(dp) :: costs(3),threshold
      logical :: changed
      moving=local
      allocate(moving_t%fields(3))
      do peer=1,3
        if (owners(peer)/=rank) cycle
        n=tree%nodes(peer)%shape
        allocate(moving_t%fields(peer)%values(1,n(1),n(2),n(3)))
        moving_t%fields(peer)%values=800._dp+peer
      end do
      destinations=mod(owners+1,nranks)
      call redistribute_mpi_amr_tree_3d(tree,owners,destinations,2,moving,moving_t,MPI_COMM_WORLD,ok)
      call require(ok,'direct U/T owner migration')
      do peer=1,3
        call require(allocated(moving%fields(peer)%values).eqv.(destinations(peer)==rank),'new owner state only')
        call require(allocated(moving_t%fields(peer)%values).eqv.(destinations(peer)==rank),'new owner T only')
        if (destinations(peer)/=rank) cycle
        n=tree%nodes(peer)%shape
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              call require(all(moving%fields(peer)%values(:,i,j,k)==values(peer,i,j,k)+10*owners(peer)), &
                'exact migrated cell/component values')
            end do
          end do
        end do
        call require(all(moving_t%fields(peer)%values==800._dp+peer),'exact migrated temperatures')
      end do
      backup=moving;backup_t=moving_t;invalid=owners
      if (rank==0) invalid(1)=mod(owners(1)+1,nranks)
      call redistribute_mpi_amr_tree_3d(tree,destinations,invalid,2,moving,moving_t,MPI_COMM_WORLD,ok)
      call require(.not.ok,'inconsistent new owner map rejected')
      do peer=1,3
        if (destinations(peer)/=rank) cycle
        call require(all(moving%fields(peer)%values==backup%fields(peer)%values),'failed migration U rollback')
        call require(all(moving_t%fields(peer)%values==backup_t%fields(peer)%values),'failed migration T rollback')
      end do
      if (rank==destinations(1)) moving_t%fields(1)%values(1,1,1,1)=-1
      call redistribute_mpi_amr_tree_3d(tree,destinations,owners,2,moving,moving_t,MPI_COMM_WORLD,ok)
      call require(.not.ok,'bad temperature rejects entire migration')
      if (rank==destinations(1)) moving_t%fields(1)%values(1,1,1,1)=801._dp
      call redistribute_mpi_amr_tree_3d(tree,destinations,owners,2,moving,moving_t,MPI_COMM_WORLD,ok)
      call require(ok,'reverse migration')
      do peer=1,3
        if (owners(peer)/=rank) cycle
        call require(all(moving%fields(peer)%values==local%fields(peer)%values),'roundtrip exact state')
      end do
      call redistribute_mpi_amr_tree_3d(tree,owners,owners,2,moving,moving_t,MPI_COMM_WORLD,ok)
      call require(ok,'identity owner migration')
      destinations=0
      call redistribute_mpi_amr_tree_3d(tree,owners,destinations,2,moving,moving_t,MPI_COMM_WORLD,ok)
      call require(ok,'concentrate before cost balancing')
      costs=-1._dp ! Nonowner entries are deliberately invalid and must be ignored.
      if (rank==0) costs=[9._dp,6._dp,3._dp]
      call rebalance_mpi_amr_tree_3d(tree,destinations,2,moving,moving_t,costs,0.75_dp,MPI_COMM_WORLD,changed,ok)
      call require(ok.and..not.changed.and.all(destinations==0),'gain threshold retains owners')
      if (rank==0) costs(1)=-1._dp
      call rebalance_mpi_amr_tree_3d(tree,destinations,2,moving,moving_t,costs,0.25_dp,MPI_COMM_WORLD,changed,ok)
      call require(.not.ok.and..not.changed.and.all(destinations==0),'negative owner cost rejected')
      if (rank==0) costs(1)=9._dp
      threshold=0.25_dp
      if (rank==0) threshold=0.5_dp
      call rebalance_mpi_amr_tree_3d(tree,destinations,2,moving,moving_t,costs,threshold,MPI_COMM_WORLD,changed,ok)
      call require(.not.ok.and..not.changed,'rank-dependent gain rejected')
      if (rank==0) moving_t%fields(1)%values(1,1,1,1)=-1._dp
      call rebalance_mpi_amr_tree_3d(tree,destinations,2,moving,moving_t,costs,0.25_dp,MPI_COMM_WORLD,changed,ok)
      call require(.not.ok.and..not.changed.and.all(destinations==0),'migration failure preserves owners')
      if (rank==0) then
        call require(moving_t%fields(1)%values(1,1,1,1)==-1._dp,'failed balance preserves temperature')
        moving_t%fields(1)%values(1,1,1,1)=801._dp
      end if
      call rebalance_mpi_amr_tree_3d(tree,destinations,2,moving,moving_t,costs,0.25_dp,MPI_COMM_WORLD,changed,ok)
      call require(ok.and.changed,'cost balancing migrates fields')
      invalid=[0,1,1]
      if (nranks>=3) invalid=[0,1,2]
      call require(all(destinations==invalid),'largest-cost-first owners')
      do peer=1,3
        call require(allocated(moving%fields(peer)%values).eqv.(destinations(peer)==rank),'balanced owner only')
        if (destinations(peer)/=rank) cycle
        n=tree%nodes(peer)%shape
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              call require(all(moving%fields(peer)%values(:,i,j,k)==values(peer,i,j,k)+10*owners(peer)), &
                'cost balance exact conserved values')
            end do
          end do
        end do
        call require(all(moving_t%fields(peer)%values==800._dp+peer),'cost balance exact temperature')
      end do
      costs=[9._dp,6._dp,3._dp]
      call rebalance_mpi_amr_tree_3d(tree,destinations,2,moving,moving_t,costs,0._dp,MPI_COMM_WORLD,changed,ok)
      call require(ok.and..not.changed.and.all(destinations==invalid),'no improvement retains owners')
      costs=0._dp
      call rebalance_mpi_amr_tree_3d(tree,destinations,2,moving,moving_t,costs,0._dp,MPI_COMM_WORLD,changed,ok)
      call require(ok.and..not.changed.and.all(destinations==invalid),'zero measured work retains owners')
    end block
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
  end do
  if (rank==root) print *, 'PASS: distributed tree payloads, exact gather, empty ranks and collective rollback'
  call MPI_Finalize(ierr)
contains
  function values(peer,i,j,k) result(value)
    integer, intent(in) :: peer,i,j,k
    real(dp) :: value(2)
    value=[real(peer,dp)+.25_dp,real(peer,dp)-.5_dp]+.01_dp*i+.001_dp*j+.0001_dp*k
  end function
  subroutine unchanged()
    integer :: peer
    do peer=1,3
      call require(allocated(local%fields(peer)%values).eqv.allocated(saved%fields(peer)%values),'allocation rollback')
      if (.not.allocated(saved%fields(peer)%values)) cycle
      call require(all(local%fields(peer)%values==saved%fields(peer)%values),'payload rollback')
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
