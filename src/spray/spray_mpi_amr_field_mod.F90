module spray_mpi_amr_field_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d,fetch_mpi_amr_tree_cells_3d
  use spray_amr_tree_partition_mod, only: locate_amr_tree_spray_cells
  use spray_amr_tree_coupling_mod, only: amr_tree_spray_stencil
  implicit none
  private
  type :: particle_stencil
    integer, allocatable :: cells(:,:)
    real(dp), allocatable :: weights(:),volumes(:)
  end type
  public :: sample_mpi_amr_tree_spray_field
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
    integer, allocatable :: lookup(:,:),requests(:,:)
    real(dp), allocatable :: values(:,:),candidate(:,:),term(:)
    real(dp) :: bounds(6),reference_bounds(6)
    integer :: flags(3),reference_flags(3),ierr,limit,p,q,c,total,n,ios
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
    allocate(requests(4,total),candidate(nvar,size(positions,2)),term(nvar),stat=ios)
    valid=ios==0;call agree()
    if (.not.agreed) return
    total=0
    do p=1,size(stencils)
      n=size(stencils(p)%weights)
      requests(:,total+1:total+n)=stencils(p)%cells;total=total+n
    end do
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
end module
