program test_mpi_amr_patch_tree_flux_peers_3d
  use mpi_f08
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,initialize_amr_patch_tree_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d,face_marks
  use mpi_amr_patch_tree_flux_peers_3d_mod, only: synchronize_mpi_amr_tree_flux_peers_3d
  implicit none
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patches(2)
  type(amr_tree_flux_3d), allocatable :: flux(:),saved(:)
  type(face_marks), allocatable :: marks(:),saved_marks(:)
  integer :: rank,nranks,ierr,node,n(3),owners(3),mode,depth,j,k
  character(len=24) :: faces(6)
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks>=2,'rank count')
  call initialize_amr_patch_3d(4,4,4,1,2,1,3,1,4,2,patches(1),ok)
  call require(ok,'left child')
  call initialize_amr_patch_3d(4,4,4,3,4,2,4,1,4,2,patches(2),ok)
  call require(ok,'offset right child')
  call initialize_amr_patch_tree_3d([4,4,4],[1,1],patches,tree,ok)
  call require(ok,'tree')
  do mode=1,4
    owners=[0,0,1]
    if (mode==4) owners=0
    faces='periodic'
    if (mode==2) faces='outflow'
    allocate(flux(3),marks(3))
    do node=1,3
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      allocate(flux(node)%x(2,n(1)+1,n(2),n(3)),flux(node)%y(2,n(1),n(2)+1,n(3)), &
        flux(node)%z(2,n(1),n(2),n(3)+1))
      allocate(marks(node)%x(n(1)+1,n(2),n(3)),marks(node)%y(n(1),n(2)+1,n(3)),marks(node)%z(n(1),n(2),n(3)+1))
      flux(node)%x=1;flux(node)%y=2;flux(node)%z=3
      marks(node)%x=.false.;marks(node)%y=.false.;marks(node)%z=.false.
    end do
    ! Partial transverse overlap: left j=3:6 corresponds to right j=1:4.
    ! Test both priority directions, including same-owner communication bypass.
    if (mode/=3.and.rank==owners(2)) then
      flux(2)%x(:,5,3:6,:)=7
      marks(2)%x(5,3:6,:)=.true.
    else if (mode==3.and.rank==owners(3)) then
      flux(3)%x(:,1,1:4,:)=7
      marks(3)%x(1,1:4,:)=.true.
    end if
    saved=flux;saved_marks=marks
    call synchronize_mpi_amr_tree_flux_peers_3d(tree,owners,2,1,faces,flux,marks,MPI_COMM_WORLD,ok)
    call require(ok,'reconcile shared faces')
    do node=1,3
      if (owners(node)/=rank) then
        call require(.not.allocated(flux(node)%x),'nonowner flux absent')
        cycle
      end if
      if (node==2) then
        saved(node)%x(:,5,3:6,:)=7;saved_marks(node)%x(5,3:6,:)=.true.
      else if (node==3) then
        saved(node)%x(:,1,1:4,:)=7;saved_marks(node)%x(1,1:4,:)=.true.
      end if
    end do
    call unchanged()
    ! A mismatching equal-priority interface must roll back all candidates.
    if (rank==owners(3)) flux(3)%x(:,1,1,1)=8
    saved=flux;saved_marks=marks
    call synchronize_mpi_amr_tree_flux_peers_3d(tree,owners,2,1,faces,flux,marks,MPI_COMM_WORLD,ok)
    call require(.not.ok,'inconsistent refined faces rejected')
    call unchanged()
    if (rank==owners(3)) flux(3)%x(:,1,1,1)=7
    ! Root opposite periodic faces are reconciled on the same node/owner.
    if (rank==owners(1)) then
      flux(1)%y(:,:,1,:)=11;marks(1)%y(:,1,:)=.true.
      flux(1)%z(:,:,:,5)=13;marks(1)%z(:,:,5)=.true.
    end if
    saved=flux;saved_marks=marks
    call synchronize_mpi_amr_tree_flux_peers_3d(tree,owners,2,0,faces,flux,marks,MPI_COMM_WORLD,ok)
    call require(ok,'root periodic self or physical faces')
    if (rank==owners(1).and.mode/=2) then
      saved(1)%y(:,:,5,:)=11;saved_marks(1)%y(:,5,:)=.true.
      saved(1)%z(:,:,:,1)=13;saved_marks(1)%z(:,:,1)=.true.
    end if
    call unchanged()
    depth=1
    if (rank==nranks-1) depth=0
    call synchronize_mpi_amr_tree_flux_peers_3d(tree,owners,2,depth,faces,flux,marks,MPI_COMM_WORLD,ok)
    call require(.not.ok,'inconsistent active depth')
    call unchanged()
    if (rank==nranks-1) faces='outflow'
    if (mode==2.and.rank==nranks-1) faces='periodic'
    call synchronize_mpi_amr_tree_flux_peers_3d(tree,owners,2,1,faces,flux,marks,MPI_COMM_WORLD,ok)
    call require(.not.ok,'inconsistent periodic geometry')
    call unchanged()
    deallocate(flux,marks)
  end do
  if (rank==0) print *, 'PASS: owner peer flux exchange, partial overlaps, refined priority, periodic self and rollback'
  call MPI_Finalize(ierr)
contains
  subroutine unchanged()
    integer :: peer
    do peer=1,3
      if (owners(peer)/=rank) cycle
      call require(all(flux(peer)%x==saved(peer)%x).and.all(flux(peer)%y==saved(peer)%y).and. &
        all(flux(peer)%z==saved(peer)%z),'exact expected flux including untouched entries')
      call require(all(marks(peer)%x.eqv.saved_marks(peer)%x).and.all(marks(peer)%y.eqv.saved_marks(peer)%y).and. &
        all(marks(peer)%z.eqv.saved_marks(peer)%z),'exact refined masks')
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
