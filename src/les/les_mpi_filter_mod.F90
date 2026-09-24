module les_mpi_filter_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use spray_partition_mod, only: spray_cartesian_block
  use les_filter_3d_mod, only: filter_les_field_3d
  implicit none
  private
  public :: filter_mpi_les_field_3d,sample_mpi_les_x_planes_3d
contains
  ! Fetch requested global (one-based) x planes from their slab owners.
  ! Coordinates may repeat and request counts may differ, including zero.
  ! Physical mapping/reflection belongs to the caller; no global field gather.
  subroutine sample_mpi_les_x_planes_3d(comm,field,cells,coordinates,planes,ok)
    type(MPI_Comm), intent(in) :: comm
    real(dp), intent(in) :: field(:,:,:,:)
    integer, intent(in) :: cells(3),coordinates(:)
    real(dp), allocatable, intent(inout) :: planes(:,:,:,:)
    logical, intent(out) :: ok
    integer :: rank,nranks,ierr,ios,nc,nrequest,first(3),extent(3),config(4),low(4),high(4)
    integer :: peer,i,j,owner,slot,received,plane_size,local_x
    integer, allocatable :: send_count(:),recv_count(:),send_offset(:),recv_offset(:),cursor(:)
    integer, allocatable :: requested(:),incoming(:),permutation(:),starts(:),ends(:)
    real(dp), allocatable :: answers(:),returned(:),candidate(:,:,:,:)
    logical :: valid,agreed
    ok=.false.;nc=size(field,1);nrequest=size(coordinates)
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if(ierr/=MPI_SUCCESS) return
    config=[cells,nc]
    call MPI_Allreduce(config,low,4,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(config,high,4,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    valid=all(low==high).and.all(cells>=1).and.nc>=1.and.nc<=1024.and.nrequest<=256
    valid=valid.and.all(coordinates>=1).and.all(coordinates<=cells(1))
    call agree()
    if(.not.agreed) return
    call spray_cartesian_block(cells,[nranks,1,1],rank,first,extent,valid)
    valid=valid.and.all(extent>=1).and.all(shape(field)==[nc,extent])
    valid=valid.and.all(ieee_is_finite(field))
    ! Bound products before forming default-integer MPI counts/displacements.
    valid=valid.and.int(nc,int64)*cells(2)<=huge(1)/int(cells(3),int64)
    call agree()
    if(.not.agreed) return
    plane_size=nc*cells(2)*cells(3)
    valid=int(nrequest,int64)*plane_size<=huge(1)
    call agree()
    if(.not.agreed) return
    allocate(send_count(nranks),recv_count(nranks),send_offset(nranks),recv_offset(nranks), &
      cursor(nranks),starts(nranks),ends(nranks),requested(nrequest),permutation(nrequest),stat=ios)
    valid=ios==0;call agree()
    if(.not.agreed) return
    send_count=0
    do peer=1,nranks
      call spray_cartesian_block(cells,[nranks,1,1],peer-1,first,extent,valid)
      if(.not.valid) exit
      starts(peer)=first(1)+1;ends(peer)=first(1)+extent(1)
    end do
    call agree()
    if(.not.agreed) return
    do i=1,nrequest
      do peer=1,nranks
        if(coordinates(i)<starts(peer).or.coordinates(i)>ends(peer)) cycle
        send_count(peer)=send_count(peer)+1
        exit
      end do
    end do
    call MPI_Alltoall(send_count,1,MPI_INTEGER,recv_count,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    valid=sum(int(recv_count,int64))<=huge(1)/int(plane_size,int64)
    call agree()
    if(.not.agreed) return
    send_offset(1)=0;recv_offset(1)=0
    do peer=2,nranks
      send_offset(peer)=send_offset(peer-1)+send_count(peer-1)
      recv_offset(peer)=recv_offset(peer-1)+recv_count(peer-1)
    end do
    received=sum(recv_count);cursor=send_offset
    do i=1,nrequest
      do owner=1,nranks
        if(coordinates(i)<starts(owner).or.coordinates(i)>ends(owner)) cycle
        cursor(owner)=cursor(owner)+1;slot=cursor(owner)
        requested(slot)=coordinates(i);permutation(slot)=i
        exit
      end do
    end do
    allocate(incoming(received),answers(received*plane_size),returned(nrequest*plane_size), &
      candidate(nc,nrequest,cells(2),cells(3)),stat=ios)
    valid=ios==0;call agree()
    if(.not.agreed) return
    call MPI_Alltoallv(requested,send_count,send_offset,MPI_INTEGER, &
      incoming,recv_count,recv_offset,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    valid=all(incoming>=starts(rank+1)).and.all(incoming<=ends(rank+1))
    call agree()
    if(.not.agreed) return
    do i=1,received
      local_x=incoming(i)-starts(rank+1)+1
      answers((i-1)*plane_size+1:i*plane_size)=reshape(field(:,local_x,:,:),[plane_size])
    end do
    send_count=send_count*plane_size;recv_count=recv_count*plane_size
    send_offset=send_offset*plane_size;recv_offset=recv_offset*plane_size
    call MPI_Alltoallv(answers,recv_count,recv_offset,MPI_DOUBLE_PRECISION, &
      returned,send_count,send_offset,MPI_DOUBLE_PRECISION,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    do i=1,nrequest
      j=permutation(i)
      candidate(:,j,:,:)=reshape(returned((i-1)*plane_size+1:i*plane_size),[nc,cells(2),cells(3)])
    end do
    call move_alloc(candidate,planes);ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine sample_mpi_les_x_planes_3d

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
