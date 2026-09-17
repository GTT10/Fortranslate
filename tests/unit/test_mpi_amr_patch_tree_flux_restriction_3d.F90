program test_mpi_amr_patch_tree_flux_restriction_3d
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d,amr_tree_integer_layout_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d,face_marks,reflux_amr_patch_tree_3d
  use mpi_amr_patch_tree_reflux_3d_mod, only: reflux_mpi_amr_tree_3d
  use mpi_amr_patch_tree_flux_restriction_3d_mod, only: restrict_mpi_amr_tree_flux_level_3d
  use mpi_amr_patch_tree_flux_peers_3d_mod, only: synchronize_mpi_amr_tree_flux_peers_3d
  implicit none
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patches(3)
  type(amr_tree_state_3d) :: state
  type(amr_tree_flux_3d), allocatable :: global(:),local(:),oracle(:),saved(:)
  type(face_marks), allocatable :: marks(:),saved_marks(:)
  integer(int64), allocatable :: origin(:,:),domain(:,:)
  integer, allocatable :: ratios(:)
  integer :: rank,nranks,ierr,node,n(3),i,j,k,mode,depth,owners(4),bad_depth
  real(dp) :: budget(2),maximum
  character(len=24) :: faces(6)
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks>=2,'rank count')
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patches(1),ok)
  call require(ok,'left child')
  call initialize_amr_patch_3d(4,4,4,3,4,1,4,1,4,2,patches(2),ok)
  call require(ok,'right child')
  call initialize_amr_patch_3d(4,8,8,2,2,2,2,2,2,3,patches(3),ok)
  call require(ok,'grandchild')
  call initialize_amr_patch_tree_3d([4,4,4],[1,1,2],patches,tree,ok)
  call require(ok,'tree')
  call amr_tree_integer_layout_3d(tree,origin,domain,ratios,ok)
  call require(ok,'layout')
  do mode=1,3
    owners=[0,1,0,0]
    if (mode==2) owners=mod([0,1,2,0],nranks)
    if (mode==3) owners=0
    faces='periodic'
    if (mode==2) faces='outflow'
    allocate(global(4),local(4),marks(4),state%fields(4))
    do node=1,4
      n=tree%nodes(node)%shape
      allocate(global(node)%x(2,n(1)+1,n(2),n(3)),global(node)%y(2,n(1),n(2)+1,n(3)), &
        global(node)%z(2,n(1),n(2),n(3)+1),state%fields(node)%values(2,n(1),n(2),n(3)))
      state%fields(node)%values=100
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)+1
            global(node)%x(:,i,j,k)=value(node,[i,j,k],1)
          end do
        end do
      end do
      do k=1,n(3)
        do j=1,n(2)+1
          do i=1,n(1)
            global(node)%y(:,i,j,k)=value(node,[i,j,k],2)
          end do
        end do
      end do
      do k=1,n(3)+1
        do j=1,n(2)
          do i=1,n(1)
            global(node)%z(:,i,j,k)=value(node,[i,j,k],3)
          end do
        end do
      end do
      if (owners(node)/=rank) cycle
      local(node)=global(node)
      allocate(marks(node)%x(n(1)+1,n(2),n(3)),marks(node)%y(n(1),n(2)+1,n(3)),marks(node)%z(n(1),n(2),n(3)+1))
      marks(node)%x=.false.;marks(node)%y=.false.;marks(node)%z=.false.
    end do
    call reflux_amr_patch_tree_3d(tree,state,global,[1.0_dp,1.0_dp,1.0_dp],faces,budget,maximum,ok, &
      synchronized_flux=oracle)
    call require(ok,'serial full reflux flux oracle')
    do depth=2,1,-1
      call synchronize_mpi_amr_tree_flux_peers_3d(tree,owners,2,depth,faces,local,marks,MPI_COMM_WORLD,ok)
      call require(ok,'child peer synchronization')
      call restrict_mpi_amr_tree_flux_level_3d(tree,owners,2,depth,local,marks,MPI_COMM_WORLD,ok)
      call require(ok,'owner to parent flux averaging')
    end do
    call synchronize_mpi_amr_tree_flux_peers_3d(tree,owners,2,0,faces,local,marks,MPI_COMM_WORLD,ok)
    call require(ok,'root synchronization')
    do node=1,4
      if (owners(node)/=rank) then
        call require(.not.allocated(local(node)%x),'only owner stores flux')
        cycle
      end if
      call require(all(local(node)%x==oracle(node)%x).and.all(local(node)%y==oracle(node)%y).and. &
        all(local(node)%z==oracle(node)%z),'exact serial synchronized flux parity')
    end do
    if (rank==owners(1)) call require(all(marks(1)%x).and.all(marks(1)%y).and.all(marks(1)%z),'root coverage marks')
    if (rank==owners(4)) call require(.not.any(marks(4)%x).and..not.any(marks(4)%y).and. &
      .not.any(marks(4)%z),'leaf marks unchanged')
    call verify_full_reflux()
    saved=local;saved_marks=marks
    call restrict_mpi_amr_tree_flux_level_3d(tree,owners,2,1,local,marks,MPI_COMM_WORLD,ok)
    call require(ok,'matching previously restricted faces accepted')
    call unchanged()
    if (rank==owners(1)) local(1)%z(:,4,4,4)=99
    saved=local
    call restrict_mpi_amr_tree_flux_level_3d(tree,owners,2,1,local,marks,MPI_COMM_WORLD,ok)
    call require(.not.ok,'late parent marked-face mismatch rejected')
    call unchanged()
    bad_depth=1
    if (rank==nranks-1) bad_depth=2
    call restrict_mpi_amr_tree_flux_level_3d(tree,owners,2,bad_depth,local,marks,MPI_COMM_WORLD,ok)
    call require(.not.ok,'different valid child depths rejected')
    call unchanged()
    call restrict_mpi_amr_tree_flux_level_3d(tree,owners,2,0,local,marks,MPI_COMM_WORLD,ok)
    call require(.not.ok,'root has no parent')
    call unchanged()
    deallocate(global,local,marks,state%fields)
  end do
  if (rank==0) print *, 'PASS: owner flux restriction, recursive peer composition, serial parity and rollback'
  call MPI_Finalize(ierr)
contains
  subroutine verify_full_reflux()
    type(amr_tree_state_3d) :: serial_state,local_state,saved_state
    type(amr_tree_flux_3d), allocatable :: raw(:),expected(:),effective(:),saved_effective(:)
    real(dp) :: serial_budget(2),mpi_budget(2),serial_max,mpi_max
    integer :: floor_depth,peer
    allocate(raw(4))
    do peer=1,4
      if (owners(peer)==rank) raw(peer)=global(peer)
    end do
    do floor_depth=0,2
      serial_state=state
      allocate(local_state%fields(4))
      do peer=1,4
        serial_state%fields(peer)%values=100
        if (owners(peer)==rank) local_state%fields(peer)=serial_state%fields(peer)
      end do
      call reflux_amr_patch_tree_3d(tree,serial_state,global,[1.0_dp,1.0_dp,1.0_dp],faces,serial_budget,serial_max,ok, &
        minimum_depth=floor_depth,synchronized_flux=expected)
      call require(ok,'serial interval reflux oracle')
      call reflux_mpi_amr_tree_3d(tree,owners,local_state,2,raw,[1.0_dp,1.0_dp,1.0_dp],faces,mpi_budget,mpi_max, &
        MPI_COMM_WORLD,ok,minimum_depth=floor_depth,synchronized_flux=effective)
      call require(ok,'complete distributed reflux')
      call require(all(mpi_budget==serial_budget).and.mpi_max==serial_max,'exact global boundary and maximum correction')
      if (floor_depth==0) call require(serial_max>0,'fixture exercises nonzero conserved-state correction')
      if (floor_depth==0.and.mode==2) call require(any(abs(serial_budget)>0),'nonzero physical boundary correction')
      do peer=1,4
        if (owners(peer)/=rank) then
          call require(.not.allocated(local_state%fields(peer)%values),'only owned corrected states')
          call require(.not.allocated(effective(peer)%x),'only owned synchronized fluxes')
          cycle
        end if
        call require(all(local_state%fields(peer)%values==serial_state%fields(peer)%values),'exact corrected state parity')
        if (tree%nodes(peer)%level<floor_depth) then
          call require(.not.allocated(effective(peer)%x),'inactive effective flux absent')
        else
          call require(all(effective(peer)%x==expected(peer)%x).and.all(effective(peer)%y==expected(peer)%y).and. &
            all(effective(peer)%z==expected(peer)%z),'complete reflux effective flux parity')
        end if
      end do
      if (floor_depth==0) then
        saved_state=local_state;saved_effective=effective
        if (rank==owners(3)) raw(3)%x(:,1,1,1)=raw(3)%x(:,1,1,1)+1
        call reflux_mpi_amr_tree_3d(tree,owners,local_state,2,raw,[1.0_dp,1.0_dp,1.0_dp],faces,mpi_budget,mpi_max, &
          MPI_COMM_WORLD,ok,synchronized_flux=effective)
        call require(.not.ok,'later peer failure rolls back full reflux')
        call require(all(mpi_budget==serial_budget).and.mpi_max==serial_max,'diagnostic rollback')
        do peer=1,4
          if (owners(peer)/=rank) cycle
          call require(all(local_state%fields(peer)%values==saved_state%fields(peer)%values),'complete state rollback')
          call require(all(effective(peer)%x==saved_effective(peer)%x).and. &
            all(effective(peer)%y==saved_effective(peer)%y).and.all(effective(peer)%z==saved_effective(peer)%z), &
            'effective flux rollback')
        end do
        if (rank==owners(3)) raw(3)=global(3)
      end if
      deallocate(local_state%fields)
    end do
  end subroutine
  function value(peer,index,axis) result(v)
    integer, intent(in) :: peer,index(3),axis
    real(dp) :: v(2)
    v=[1.0_dp,2.0_dp]*real(axis+tree%nodes(peer)%level,dp)
    if (mode==2) v=v+.01_dp*real(tree%nodes(peer)%level+1,dp)*sum(real(origin(:,peer)+int(index-1,int64),dp)/ &
      real(domain(:,tree%nodes(peer)%level+1),dp))
  end function
  subroutine unchanged()
    integer :: peer
    do peer=1,4
      if (owners(peer)/=rank) cycle
      call require(all(local(peer)%x==saved(peer)%x).and.all(local(peer)%y==saved(peer)%y).and. &
        all(local(peer)%z==saved(peer)%z),'flux rollback')
      call require(all(marks(peer)%x.eqv.saved_marks(peer)%x).and.all(marks(peer)%y.eqv.saved_marks(peer)%y).and. &
        all(marks(peer)%z.eqv.saved_marks(peer)%z),'mark rollback')
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
