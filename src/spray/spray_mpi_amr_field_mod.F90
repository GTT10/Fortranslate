module spray_mpi_amr_field_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d,fetch_mpi_amr_tree_cells_3d
  use mpi_amr_patch_tree_restriction_3d_mod, only: average_down_mpi_amr_tree_3d
  use spray_amr_tree_partition_mod, only: locate_amr_tree_spray_cells
  use spray_amr_tree_coupling_mod, only: amr_tree_spray_stencil
  implicit none
  private
  type :: particle_stencil
    integer, allocatable :: cells(:,:)
    real(dp), allocatable :: weights(:),volumes(:)
  end type
  public :: sample_mpi_amr_tree_spray_field,deposit_mpi_amr_tree_spray_field
contains
  ! Each rank supplies its own particle positions, possibly empty. Only requested
  ! exposed-cell values are exchanged; the owner-local field is never gathered.
  subroutine sample_mpi_amr_tree_spray_field(comm,tree,owners,state,nvar,positions,lower,spacing,periodic, &
      samples,ok,maximum_requests)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar
    type(amr_tree_state_3d), intent(in) :: state
    real(dp), intent(in) :: positions(:,:),lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    real(dp), allocatable, intent(inout) :: samples(:,:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_requests
    type(particle_stencil), allocatable :: stencils(:)
    integer, allocatable :: requests(:,:)
    real(dp), allocatable :: values(:,:),candidate(:,:),term(:)
    integer :: ierr,p,q,c,total,ios
    logical :: valid,agreed
    ok=.false.
    call prepare_requests(comm,tree,owners,state,nvar,positions,lower,spacing,periodic, &
      stencils,requests,valid,maximum_requests)
    if (.not.valid) return
    allocate(candidate(nvar,size(positions,2)),term(nvar),stat=ios)
    valid=ios==0;call agree()
    if (.not.agreed) return
    call fetch_mpi_amr_tree_cells_3d(tree,owners,state,nvar,requests,values,comm,valid)
    if (.not.valid) return
    total=0;candidate=0._dp
    do p=1,size(stencils)
      do q=1,size(stencils(p)%weights)
        total=total+1;term=stencils(p)%weights(q)*values(:,total)
        do c=1,nvar
          if (term(c)>0._dp) then
            if (candidate(c,p)>huge(1._dp)-term(c)) valid=.false.
          else if (term(c)<0._dp) then
            if (candidate(c,p)< -huge(1._dp)-term(c)) valid=.false.
          end if
        end do
        if (.not.valid) exit
        candidate(:,p)=candidate(:,p)+term
      end do
      if (.not.valid) exit
    end do
    call agree()
    if (.not.agreed) return
    call move_alloc(candidate,samples);ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine
  ! Contributions are integrated quantities. Destination bins preserve source
  ! particle/stencil order; receives are applied in increasing source-rank order.
  subroutine deposit_mpi_amr_tree_spray_field(comm,tree,owners,state,nvar,positions,contributions, &
      lower,spacing,periodic,ok,maximum_requests,maximum_received)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar
    type(amr_tree_state_3d), intent(inout) :: state
    real(dp), intent(in) :: positions(:,:),contributions(:,:),lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_requests,maximum_received
    type(particle_stencil), allocatable :: stencils(:)
    type(amr_tree_state_3d) :: candidate
    integer, allocatable :: requests(:,:),send_cells(:,:),recv_cells(:,:)
    integer, allocatable :: counts(:),received(:),offsets(:),displacements(:),cursor(:)
    real(dp), allocatable :: send_values(:,:),recv_values(:,:),term(:)
    integer :: rank,nranks,ierr,ios,capacity,total,incoming,p,q,c,r,slot,node,idx(3)
    logical :: valid,agreed,has_particles,any_particles
    ok=.false.;capacity=1000000
    if (present(maximum_received)) capacity=maximum_received
    call prepare_requests(comm,tree,owners,state,nvar,positions,lower,spacing,periodic, &
      stencils,requests,valid,maximum_requests)
    if (.not.valid) return
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=all(shape(contributions)==[nvar,size(positions,2)]).and.capacity>=0.and.capacity<=1000000
    valid=valid.and.all(ieee_is_finite(contributions))
    do node=1,size(owners)
      if (owners(node)/=rank) cycle
      if (.not.all(ieee_is_finite(state%fields(node)%values))) valid=.false.
    end do
    call agree()
    if (.not.agreed) return
    has_particles=size(positions,2)>0
    call MPI_Allreduce(has_particles,any_particles,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (.not.any_particles) then
      ok=.true.;return
    end if
    total=size(requests,2)
    allocate(counts(nranks),received(nranks),offsets(nranks),displacements(nranks),cursor(nranks), &
      send_cells(4,total),send_values(nvar,total),term(nvar),stat=ios)
    valid=ios==0;call agree()
    if (.not.agreed) return
    counts=0
    do q=1,total
      r=owners(requests(1,q))+1;counts(r)=counts(r)+1
    end do
    call MPI_Alltoall(counts,1,MPI_INTEGER,received,1,MPI_INTEGER,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=sum(int(received,int64))<=int(min(capacity,huge(1)/max(4,nvar)),int64)
    call agree()
    if (.not.agreed) return
    incoming=sum(received);offsets(1)=0;displacements(1)=0
    do r=2,nranks
      offsets(r)=offsets(r-1)+counts(r-1)
      displacements(r)=displacements(r-1)+received(r-1)
    end do
    cursor=offsets
    do p=1,size(stencils)
      do q=1,size(stencils(p)%weights)
        node=stencils(p)%cells(1,q);r=owners(node)+1
        slot=cursor(r)+1;cursor(r)=slot
        send_cells(:,slot)=stencils(p)%cells(:,q)
        term=stencils(p)%weights(q)*contributions(:,p)
        if (any(abs(term)>huge(1._dp)*min(1._dp,stencils(p)%volumes(q)))) then
          valid=.false.;exit
        end if
        send_values(:,slot)=term/stencils(p)%volumes(q)
      end do
      if (.not.valid) exit
    end do
    call agree()
    if (.not.agreed) return
    allocate(recv_cells(4,incoming),recv_values(nvar,incoming),stat=ios)
    valid=ios==0;call agree()
    if (.not.agreed) return
    call MPI_Alltoallv(send_cells,4*counts,4*offsets,MPI_INTEGER, &
      recv_cells,4*received,4*displacements,MPI_INTEGER,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Alltoallv(send_values,nvar*counts,nvar*offsets,MPI_DOUBLE_PRECISION, &
      recv_values,nvar*received,nvar*displacements,MPI_DOUBLE_PRECISION,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    candidate=state
    do q=1,incoming
      node=recv_cells(1,q);idx=recv_cells(2:4,q);term=recv_values(:,q)
      do c=1,nvar
        associate(value=>candidate%fields(node)%values(c,idx(1),idx(2),idx(3)))
        if (term(c)>0._dp) then
          if (value>huge(1._dp)-term(c)) valid=.false.
        else if (term(c)<0._dp) then
          if (value< -huge(1._dp)-term(c)) valid=.false.
        end if
        end associate
      end do
      if (.not.valid) exit
      candidate%fields(node)%values(:,idx(1),idx(2),idx(3))= &
        candidate%fields(node)%values(:,idx(1),idx(2),idx(3))+term
    end do
    call agree()
    if (.not.agreed) return
    call average_down_mpi_amr_tree_3d(tree,owners,candidate,nvar,comm,valid)
    if (.not.valid) return
    call move_alloc(candidate%fields,state%fields);ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine

  subroutine prepare_requests(comm,tree,owners,state,nvar,positions,lower,spacing,periodic, &
      stencils,requests,ok,maximum_requests)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar
    type(amr_tree_state_3d), intent(in) :: state
    real(dp), intent(in) :: positions(:,:),lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    type(particle_stencil), allocatable, intent(out) :: stencils(:)
    integer, allocatable, intent(out) :: requests(:,:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_requests
    integer, allocatable :: lookup(:,:)
    real(dp) :: bounds(6),reference_bounds(6)
    integer :: flags(3),reference_flags(3),ierr,limit,p,total,n,ios
    logical :: valid,agreed
    ok=.false.;limit=1000000
    if (present(maximum_requests)) limit=maximum_requests
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nvar,comm,valid)
    if (.not.valid) return
    valid=size(positions,1)==3.and.limit>=0.and.limit<=1000000
    valid=valid.and.all(ieee_is_finite([lower,spacing]))
    if (valid) valid=size(positions,2)<=limit
    call agree()
    if (.not.agreed) return
    bounds=[lower,spacing];reference_bounds=bounds
    call MPI_Bcast(reference_bounds,6,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    flags=merge(1,0,periodic);reference_flags=flags
    call MPI_Bcast(reference_flags,3,MPI_INTEGER,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=all(bounds==reference_bounds).and.all(flags==reference_flags)
    call agree()
    if (.not.agreed) return
    call locate_amr_tree_spray_cells(tree,positions,lower,spacing,periodic,lookup,valid)
    call agree()
    if (.not.agreed) return
    allocate(stencils(size(positions,2)),stat=ios)
    valid=ios==0;call agree()
    if (.not.agreed) return
    total=0
    do p=1,size(stencils)
      if (total==limit) then
        valid=.false.;exit
      end if
      call amr_tree_spray_stencil(tree,positions(:,p),lower,spacing,periodic, &
        stencils(p)%cells,stencils(p)%weights,stencils(p)%volumes,valid,maximum_cells=limit-total)
      if (.not.valid) exit
      total=total+size(stencils(p)%weights)
    end do
    valid=valid.and.total<=huge(1)/max(4,nvar)
    call agree()
    if (.not.agreed) return
    allocate(requests(4,total),stat=ios)
    valid=ios==0;call agree()
    if (.not.agreed) return
    total=0
    do p=1,size(stencils)
      n=size(stencils(p)%weights)
      requests(:,total+1:total+n)=stencils(p)%cells;total=total+n
    end do
    ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine
end module
