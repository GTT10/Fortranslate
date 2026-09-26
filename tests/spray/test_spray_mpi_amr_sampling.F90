program test_spray_mpi_amr_sampling
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use spray_amr_tree_coupling_mod, only: sample_amr_tree_spray_field
  use spray_mpi_amr_field_mod, only: sample_mpi_amr_tree_spray_field
  implicit none
  type(amr_patch_3d) :: patches(2)
  type(amr_patch_tree_3d) :: tree
  type(amr_tree_state_3d) :: full,local,snapshot
  real(dp), allocatable :: positions(:,:),samples(:,:),reference(:,:),saved(:,:),good_positions(:,:)
  real(dp) :: points(3,4),zero(3)=0._dp,spacing(3)=1._dp,bounds(3)
  integer :: rank,nranks,ierr,node,n(3),i,j,k,np,owners(3),wrong(3),level
  logical :: ok,periodic(3)=[.true.,.true.,.true.]
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks==2.or.nranks==4,'rank count')
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patches(1),ok)
  call require(ok,'fine patch')
  call initialize_amr_patch_3d(4,8,8,1,2,1,2,1,2,2,patches(2),ok)
  call require(ok,'nested patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,tree,ok)
  call require(ok,'tree')
  allocate(full%fields(3))
  do node=1,3
    n=tree%nodes(node)%shape
    allocate(full%fields(node)%values(2,n(1),n(2),n(3)))
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      full%fields(node)%values(:,i,j,k)=[real(10*node+i+3*j+7*k,dp),real(node*i*j+k,dp)]
    end do
    end do
    end do
  end do
  owners=[0,1,nranks-1];call distribute()
  points(:,1)=[1.9_dp,1.1_dp,1.1_dp]
  points(:,2)=[2.1_dp,1.1_dp,1.1_dp]
  points(:,3)=[.9_dp,.9_dp,.9_dp]
  points(:,4)=[3.9_dp,1.1_dp,1.1_dp]
  np=2
  if (rank==0) np=0
  if (rank==1) np=4
  positions=points(:,:np)
  snapshot=full;positions=.5_dp
  do node=1,3
    full%fields(node)%values=real(10*node,dp)
  end do
  call distribute()
  do level=0,1
    call sample_amr_tree_spray_field(tree,full,positions,zero,spacing,periodic,reference,ok,level)
    call require(ok.and.all(abs(reference-real(10*(level+1),dp))<1e-13_dp),'analytic serial level view')
    call sample_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,local,2,positions,zero,spacing,periodic, &
      samples,ok,maximum_level=level)
    call require(ok.and.all(samples==reference),'exact MPI level view')
  end do
  saved=samples
  call sample_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,local,2,positions,zero,spacing,periodic, &
    samples,ok,maximum_level=merge(0,1,rank==0))
  call require(.not.ok.and.all(samples==saved),'level disagreement rollback')
  call sample_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,local,2,positions,zero,spacing,periodic, &
    samples,ok,maximum_level=3)
  call require(.not.ok.and.all(samples==saved),'invalid level rollback')
  if(rank==owners(3)) local%fields(3)%values=ieee_value(0._dp,ieee_quiet_nan)
  call sample_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,local,2,positions,zero,spacing,periodic, &
    samples,ok,maximum_level=0)
  call require(ok.and.all(abs(samples-10._dp)<1e-13_dp),'unrequested deeper payload ignored')
  full=snapshot;positions=points(:,:np);call distribute()
  call sample_amr_tree_spray_field(tree,full,positions,zero,spacing,periodic,reference,ok)
  call require(ok,'serial reference')
  call sample()
  call require(ok.and.all(samples==reference),'exact distributed composite sample')
  saved=samples
  ! Covered parent payload is not part of any composite interpolation request.
  if (rank==owners(1)) local%fields(1)%values(:,1,1,1)=ieee_value(0._dp,ieee_quiet_nan)
  call sample()
  call require(ok.and.all(samples==reference),'unrequested covered payload ignored')
  call distribute()
  if (rank==owners(2)) local%fields(2)%values=ieee_value(0._dp,ieee_quiet_nan)
  call sample()
  call require(.not.ok.and.all(samples==saved),'requested NaN collective rollback')
  call distribute()
  bounds=zero
  if (rank==0) bounds(1)=1._dp
  call sample_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,local,2,positions,bounds,spacing,periodic,samples,ok)
  call require(.not.ok.and.all(samples==saved),'physical bounds disagreement rollback')
  wrong=owners
  if (rank==0) wrong(1)=1
  call sample_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,wrong,local,2,positions,zero,spacing,periodic,samples,ok)
  call require(.not.ok.and.all(samples==saved),'owner disagreement rollback')
  if (rank==owners(2)) deallocate(local%fields(2)%values)
  call sample()
  call require(.not.ok.and.all(samples==saved),'missing owned field rollback')
  call distribute()
  good_positions=positions
  if (rank==nranks-1) positions(1,1)=ieee_value(0._dp,ieee_quiet_nan)
  call sample()
  call require(.not.ok.and.all(samples==saved),'bad local position collective rollback')
  positions=good_positions
  call sample_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,local,2,positions,zero,spacing,periodic, &
    samples,ok,maximum_requests=merge(0,1000000,rank==1))
  call require(.not.ok.and.all(samples==saved),'local request capacity rollback')
  owners=modulo(owners+1,nranks);call distribute()
  call sample()
  call require(ok.and.all(samples==reference),'ownership changes preserve samples')
  deallocate(positions);allocate(positions(3,0))
  call sample_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,local,2,positions,zero,spacing,periodic, &
    samples,ok,maximum_requests=0)
  call require(ok.and.size(samples,1)==2.and.size(samples,2)==0,'globally empty query')
  call sample_amr_tree_spray_field(tree,full,positions,zero,spacing,periodic,reference,ok,0)
  call require(ok.and.size(reference,2)==0,'empty serial bounded query')
  call sample_amr_tree_spray_field(tree,full,positions,zero,spacing,periodic,reference,ok,3)
  call require(.not.ok.and.size(reference,2)==0,'empty serial invalid level rollback')
  if (rank==0) print *, 'PASS: owner-local AMR spray field sampling'
  call MPI_Finalize(ierr)
contains
  subroutine sample()
    call sample_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,local,2,positions,zero,spacing,periodic,samples,ok)
  end subroutine
  subroutine distribute()
    integer :: q
    if (allocated(local%fields)) deallocate(local%fields)
    allocate(local%fields(3))
    do q=1,3
      if (owners(q)==rank) local%fields(q)%values=full%fields(q)%values
    end do
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL rank ',rank,': ',label
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
  end subroutine
end program
