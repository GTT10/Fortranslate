program test_mpi_amr_patch_tree_cells_3d
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: fetch_mpi_amr_tree_cells_3d,add_mpi_amr_tree_cells_3d
  implicit none
  type(amr_patch_tree_3d) :: tree,changed
  type(amr_patch_3d) :: patches(2)
  type(amr_tree_state_3d) :: local,checkpoint
  integer :: rank,nranks,ierr,node,n(3),i,j,k,q,mode,owners(3),other(3),components
  integer, allocatable :: requests(:,:),invalid(:,:)
  real(dp), allocatable :: result(:,:),saved(:,:),increments(:,:)
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks>=2,'rank count')
  call initialize_amr_patch_3d(4,4,4,1,2,1,2,1,2,2,patches(1),ok)
  call require(ok,'child')
  call initialize_amr_patch_3d(4,4,4,1,1,1,1,1,1,3,patches(2),ok)
  call require(ok,'grandchild')
  call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,tree,ok)
  call require(ok,'tree')
  do mode=1,2
    owners=[0,1,0]
    if (mode==2) owners=0
    allocate(local%fields(3))
    do node=1,3
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      allocate(local%fields(node)%values(2,n(1),n(2),n(3)))
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            local%fields(node)%values(:,i,j,k)=expected(node,i,j,k)
          end do
        end do
      end do
    end do
    allocate(requests(4,5+rank))
    do q=1,size(requests,2)
      node=1+mod(q+rank,3)
      requests(:,q)=[node,1+mod(q,tree%nodes(node)%shape(1)),1,2]
    end do
    requests(:,1)=[2,2,3,4];requests(:,2)=requests(:,1)
    if (rank==nranks-1) then
      deallocate(requests)
      allocate(requests(4,0))
    end if
    call fetch_mpi_amr_tree_cells_3d(tree,owners,local,2,requests,result,MPI_COMM_WORLD,ok)
    call require(ok,'fetch mixed owners and empty request rank')
    call require(all(shape(result)==[2,size(requests,2)]),'result shape')
    do q=1,size(requests,2)
      call require(all(result(:,q)==expected(requests(1,q),requests(2,q),requests(3,q),requests(4,q))), &
        'exact original request order')
    end do
    saved=result
    invalid=requests
    if (rank==0) invalid(2,1)=0
    call fetch_mpi_amr_tree_cells_3d(tree,owners,local,2,invalid,result,MPI_COMM_WORLD,ok)
    call rejected('invalid coordinate')
    if (rank==0) invalid(1,1)=4
    call fetch_mpi_amr_tree_cells_3d(tree,owners,local,2,invalid,result,MPI_COMM_WORLD,ok)
    call rejected('invalid node')
    other=owners
    if (rank==0) other(1)=1
    call fetch_mpi_amr_tree_cells_3d(tree,other,local,2,requests,result,MPI_COMM_WORLD,ok)
    call rejected('inconsistent ownership')
    components=2
    if (rank==0) components=3
    call fetch_mpi_amr_tree_cells_3d(tree,owners,local,components,requests,result,MPI_COMM_WORLD,ok)
    call rejected('inconsistent component count')
    changed=tree
    if (rank==0) then
      changed%nodes(3)%patch%coarse_i_lower=2
      changed%nodes(3)%patch%coarse_i_upper=2
    end if
    call fetch_mpi_amr_tree_cells_3d(changed,owners,local,2,requests,result,MPI_COMM_WORLD,ok)
    call rejected('inconsistent valid topology')
    if (rank==owners(2)) local%fields(2)%values(1,2,3,4)=ieee_value(0.0_dp,ieee_quiet_nan)
    call fetch_mpi_amr_tree_cells_3d(tree,owners,local,2,requests,result,MPI_COMM_WORLD,ok)
    call rejected('nonfinite requested payload')
    if (rank==owners(2)) local%fields(2)%values(:,2,3,4)=expected(2,2,3,4)
    do node=1,3
      call require(allocated(local%fields(node)%values).eqv.(owners(node)==rank),'owned allocation only')
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            call require(all(local%fields(node)%values(:,i,j,k)==expected(node,i,j,k)),'source unchanged')
          end do
        end do
      end do
    end do
    deallocate(requests)
    allocate(requests(4,0))
    call fetch_mpi_amr_tree_cells_3d(tree,owners,local,2,requests,result,MPI_COMM_WORLD,ok)
    call require(ok,'all ranks empty requests')
    call require(all(shape(result)==[2,0]),'empty output')
    checkpoint=local
    allocate(increments(2,0))
    call add_mpi_amr_tree_cells_3d(tree,owners,local,2,requests,increments,MPI_COMM_WORLD,ok)
    call require(ok,'empty deposits')
    call unchanged('empty deposits unchanged')
    deallocate(requests,increments)
    ! Rank zero sends two signed contributions to one cell; rank one also
    ! contributes. Rank two, when present, participates without any records.
    allocate(requests(4,merge(2,0,rank==0.or.rank==1)))
    allocate(increments(2,size(requests,2)))
    do q=1,size(requests,2)
      requests(:,q)=[2,2,3,4]
    end do
    if (rank==0) then
      increments(:,1)=[1._dp,-2._dp];increments(:,2)=[-.25_dp,.5_dp]
    else if (rank==1) then
      increments(:,1)=[2._dp,-1._dp];increments(:,2)=[-.75_dp,.25_dp]
    end if
    call add_mpi_amr_tree_cells_3d(tree,owners,local,2,requests,increments,MPI_COMM_WORLD,ok)
    call require(ok,'duplicate remote signed deposits')
    do node=1,3
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if (node==2.and.i==2.and.j==3.and.k==4) then
          call require(maxval(abs(local%fields(node)%values(:,i,j,k)- &
            (expected(node,i,j,k)+[2._dp,-2.25_dp])))<1.e-14_dp,'summed owner value')
        else
          call require(all(local%fields(node)%values(:,i,j,k)==expected(node,i,j,k)),'other cells unchanged')
        end if
      end do
      end do
      end do
    end do
    checkpoint=local
    invalid=requests
    if (rank==0) invalid(3,2)=0
    call add_mpi_amr_tree_cells_3d(tree,owners,local,2,invalid,increments,MPI_COMM_WORLD,ok)
    call require(.not.ok,'invalid deposit coordinate rejected')
    call unchanged('invalid deposit rollback')
    if (rank==1) increments(2,2)=ieee_value(0._dp,ieee_quiet_nan)
    call add_mpi_amr_tree_cells_3d(tree,owners,local,2,requests,increments,MPI_COMM_WORLD,ok)
    call require(.not.ok,'nonfinite deposit rejected')
    call unchanged('nonfinite deposit rollback')
    increments=0
    ! The first record succeeds on a different cell. A late second component
    ! overflows on another owner, requiring collective rollback of both.
    if (rank==0) then
      requests(:,1)=[1,1,1,1];increments(:,1)=1
      requests(:,2)=[2,2,3,4];increments(2,2)=huge(1._dp)
    end if
    if (rank==owners(2)) local%fields(2)%values(2,2,3,4)=huge(1._dp)
    checkpoint=local
    call add_mpi_amr_tree_cells_3d(tree,owners,local,2,requests,increments,MPI_COMM_WORLD,ok)
    call require(.not.ok,'late positive overflow rejected')
    call unchanged('positive overflow collective rollback')
    if (rank==0) increments(2,2)=-huge(1._dp)
    if (rank==owners(2)) local%fields(2)%values(2,2,3,4)=-huge(1._dp)
    checkpoint=local
    call add_mpi_amr_tree_cells_3d(tree,owners,local,2,requests,increments,MPI_COMM_WORLD,ok)
    call require(.not.ok,'late negative overflow rejected')
    call unchanged('negative overflow collective rollback')
    deallocate(increments)
    deallocate(requests,local%fields)
  end do
  if (rank==0) print *, 'PASS: sparse owner fetch/deposit, duplicates, empty ranks and collective rollback'
  call MPI_Finalize(ierr)
contains
  subroutine unchanged(label)
    character(len=*), intent(in) :: label
    integer :: p
    do p=1,size(tree%nodes)
      call require(allocated(local%fields(p)%values).eqv.allocated(checkpoint%fields(p)%values),label)
      if (.not.allocated(local%fields(p)%values)) cycle
      call require(all(local%fields(p)%values==checkpoint%fields(p)%values),label)
    end do
  end subroutine
  function expected(node,i,j,k) result(value)
    integer, intent(in) :: node,i,j,k
    real(dp) :: value(2)
    value=[real(node,dp)+.25_dp,real(node,dp)-.5_dp]+.01_dp*i+.001_dp*j+.0001_dp*k
  end function
  subroutine rejected(label)
    character(len=*), intent(in) :: label
    call require(.not.ok,label)
    call require(all(shape(result)==shape(saved)),'rollback shape')
    call require(all(result==saved),'rollback values')
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
