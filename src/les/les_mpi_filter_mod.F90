module les_mpi_filter_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use spray_partition_mod, only: spray_cartesian_block
  use les_filter_3d_mod, only: filter_les_field_3d
  implicit none
  private
  public :: filter_mpi_les_field_3d
contains
  ! Periodic x slabs. One-plane propagation handles widths exceeding local or
  ! global x extent without gathering the global field on any rank.
  subroutine filter_mpi_les_field_3d(comm,field,cells,weights,filtered,ok)
    type(MPI_Comm), intent(in) :: comm
    real(dp), intent(in) :: field(:,:,:,:),weights(:)
    integer, intent(in) :: cells(3)
    real(dp), intent(inout) :: filtered(:,:,:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: left_stream(:,:,:,:),right_stream(:,:,:,:),xhalo(:,:,:,:),halo(:,:,:,:),candidate(:,:,:,:)
    real(dp), allocatable :: send_left(:,:,:),send_right(:,:,:),recv_left(:,:,:),recv_right(:,:,:),wmin(:),wmax(:)
    integer :: rank,nranks,ierr,nc,nx,ny,nz,radius,s,j,k,jp,kp,left,right,count,ios
    integer :: config(5),low(5),high(5),first(3),extent(3)
    logical :: valid,global_valid
    ok=.false.
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    nc=size(field,1);nx=size(field,2);ny=size(field,3);nz=size(field,4)
    config=[cells,nc,size(weights)]
    call MPI_Allreduce(config,low,5,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(config,high,5,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=all(low==high).and.nc>=1.and.nc<=1024.and.size(weights)>=1.and.size(weights)<=65
    valid=valid.and.all(cells>=1).and.mod(size(weights),2)==1.and.all(shape(filtered)==shape(field))
    call agree()
    if (.not.global_valid) return
    call spray_cartesian_block(cells,[nranks,1,1],rank,first,extent,valid)
    valid=valid.and.all(extent==[nx,ny,nz])
    valid=valid.and.all(ieee_is_finite(field)).and.all(ieee_is_finite(weights))
    call agree()
    if (.not.global_valid) return
    valid=maxval(abs(field))<=1e100_dp.and.maxval(abs(weights))<=1e6_dp
    call agree()
    if (.not.global_valid) return
    valid=sum(abs(weights))<=1e6_dp.and.abs(sum(weights)-1)<=1e-12_dp
    valid=valid.and.int(nc,int64)*ny*nz<=huge(count)
    call agree()
    if (.not.global_valid) return
    allocate(wmin(size(weights)),wmax(size(weights)),stat=ios)
    valid=ios==0;call agree()
    if (.not.global_valid) return
    call MPI_Allreduce(weights,wmin,size(weights),MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(weights,wmax,size(weights),MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (any(wmin/=wmax)) return
    radius=size(weights)/2;count=nc*ny*nz
    allocate(left_stream(nc,nx,ny,nz),right_stream(nc,nx,ny,nz),xhalo(nc,nx+2*radius,ny,nz), &
      halo(nc,nx+2*radius,ny+2*radius,nz+2*radius),candidate(nc,nx,ny,nz), &
      send_left(nc,ny,nz),send_right(nc,ny,nz),recv_left(nc,ny,nz),recv_right(nc,ny,nz),stat=ios)
    valid=ios==0;call agree()
    if (.not.global_valid) return
    left=modulo(rank-1,nranks);right=modulo(rank+1,nranks)
    left_stream=field;right_stream=field;xhalo(:,radius+1:radius+nx,:,:)=field
    do s=1,radius
      send_left=right_stream(:,1,:,:);send_right=left_stream(:,nx,:,:)
      call MPI_Sendrecv(send_left,count,MPI_DOUBLE_PRECISION,left,471,recv_right,count,MPI_DOUBLE_PRECISION, &
        right,471,comm,MPI_STATUS_IGNORE,ierr)
      if (ierr/=MPI_SUCCESS) return
      call MPI_Sendrecv(send_right,count,MPI_DOUBLE_PRECISION,right,472,recv_left,count,MPI_DOUBLE_PRECISION, &
        left,472,comm,MPI_STATUS_IGNORE,ierr)
      if (ierr/=MPI_SUCCESS) return
      xhalo(:,radius+1-s,:,:)=recv_left;xhalo(:,radius+nx+s,:,:)=recv_right
      left_stream(:,2:nx,:,:)=left_stream(:,1:nx-1,:,:);left_stream(:,1,:,:)=recv_left
      right_stream(:,1:nx-1,:,:)=right_stream(:,2:nx,:,:);right_stream(:,nx,:,:)=recv_right
    end do
    do k=1,nz+2*radius
      kp=modulo(k-radius-1,nz)+1
      do j=1,ny+2*radius
        jp=modulo(j-radius-1,ny)+1
        halo(:,:,j,k)=xhalo(:,:,jp,kp)
      end do
    end do
    call filter_les_field_3d(halo,weights,candidate,valid,[radius+1,radius+1,radius+1])
    call agree()
    if (.not.global_valid) return
    filtered=candidate;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,global_valid,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) global_valid=.false.
    end subroutine
  end subroutine
end module
