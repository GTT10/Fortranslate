program test_les_mpi_filter
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use les_filter_3d_mod, only: les_filter_weights,filter_les_field_3d
  use les_mpi_filter_mod, only: filter_mpi_les_field_3d,sample_mpi_les_x_planes_3d
  use spray_partition_mod, only: spray_cartesian_block
  implicit none
  integer, parameter :: cells(3)=[5,3,2]
  real(dp) :: full(2,5,3,2),reference(2,5,3,2)
  real(dp), allocatable :: local(:,:,:,:),filtered(:,:,:,:),saved(:,:,:,:),weights(:)
  integer, allocatable :: coordinates(:)
  real(dp), allocatable :: planes(:,:,:,:),saved_planes(:,:,:,:)
  integer :: rank,nranks,ierr,i,j,k,first(3),extent(3),lo,hi,m,ratio
  character(len=8) :: model
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call spray_cartesian_block(cells,[nranks,1,1],rank,first,extent,ok)
  call require(ok,'partition')
  lo=first(1)+1;hi=lo+extent(1)-1
  do k=1,2
    do j=1,3
      do i=1,5
        full(:,i,j,k)=[sin(real(i+2*j+3*k,dp)),real(i*j*k,dp)]
      end do
    end do
  end do
  local=full(:,lo:hi,:,:);allocate(filtered,mold=local)
  ! Thirteen-layer periodic/reflected requests cross multiple uneven slabs.
  allocate(coordinates(26))
  do m=1,2
    do i=1,26
      j=lo-14+i
      if(m==1) then
        coordinates(i)=modulo(j-1,cells(1))+1
      else
        k=modulo(j-1,2*cells(1))
        coordinates(i)=min(k+1,2*cells(1)-k)
      end if
    end do
    call sample_mpi_les_x_planes_3d(MPI_COMM_WORLD,local,cells,coordinates,planes,ok)
    call require(ok,'wide x donor exchange')
    do i=1,size(coordinates)
      call require(all(planes(:,i,:,:)==full(:,coordinates(i),:,:)),'exact mapped donor planes')
    end do
  end do
  saved_planes=planes
  if(rank==nranks-1) coordinates(1)=0
  call sample_mpi_les_x_planes_3d(MPI_COMM_WORLD,local,cells,coordinates,planes,ok)
  call require(.not.ok.and.all(planes==saved_planes),'remote invalid coordinate preserves planes')
  deallocate(coordinates)
  allocate(coordinates(rank))
  coordinates=cells(1)
  call sample_mpi_les_x_planes_3d(MPI_COMM_WORLD,local,cells,coordinates,planes,ok)
  call require(ok.and.size(planes,2)==rank,'unequal and zero request counts')
  do i=1,rank
    call require(all(planes(:,i,:,:)==full(:,cells(1),:,:)),'repeated last-owner donors')
  end do
  do m=1,4
    model='box';ratio=2
    if (m==2) ratio=6
    if (m==3) then
      ratio=8;model='gaussian'
    end if
    if (m==4) model='none'
    call les_filter_weights(model,ratio,weights,ok)
    call require(ok,'weights')
    call filter_les_field_3d(full,weights,reference,ok)
    call require(ok,'serial filter')
    call filter_mpi_les_field_3d(MPI_COMM_WORLD,local,cells,weights,filtered,ok)
    call require(ok.and.all(filtered==reference(:,lo:hi,:,:)),'exact serial filter parity')
  end do
  saved=filtered
  call les_filter_weights('box',2,weights,ok)
  if (rank==nranks-1) weights=[.2_dp,.6_dp,.2_dp]
  call filter_mpi_les_field_3d(MPI_COMM_WORLD,local,cells,weights,filtered,ok)
  call require(.not.ok.and.all(filtered==saved),'different valid weights collective rollback')
  call les_filter_weights('box',2,weights,ok)
  if (rank==nranks-1) call les_filter_weights('box',6,weights,ok)
  call filter_mpi_les_field_3d(MPI_COMM_WORLD,local,cells,weights,filtered,ok)
  call require(.not.ok.and.all(filtered==saved),'different widths collective rollback')
  call les_filter_weights('box',6,weights,ok)
  if (rank==nranks-1) local(1,1,1,1)=ieee_value(1._dp,ieee_quiet_nan)
  call filter_mpi_les_field_3d(MPI_COMM_WORLD,local,cells,weights,filtered,ok)
  call require(.not.ok.and.all(filtered==saved),'nonfinite remote field collective rollback')
  if (rank==0) print *, 'PASS: distributed LES filtering'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    integer :: abort_error
    if (condition) return
    print *,rank,'FAIL: ',message
    call MPI_Abort(MPI_COMM_WORLD,1,abort_error)
  end subroutine
end program
