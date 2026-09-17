module spray_mpi_exchange_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use spray_parcel_mod, only: spray_parcel,valid_spray_parcel
  use spray_mixture_mod, only: valid_component_inventory
  implicit none
  private
  public :: exchange_spray_parcel_bins
contains
  ! Exchange contiguous destination bins, then sort received IDs. Caller checks
  ! geometric ownership on the private result before publishing application state.
  subroutine exchange_spray_parcel_bins(comm,parcels,counts,offsets,maximum_parcels,ok,component_mass)
    type(MPI_Comm), intent(in) :: comm
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    integer, intent(in) :: counts(:),offsets(:),maximum_parcels
    logical, intent(out) :: ok
    real(dp), allocatable, optional, intent(inout) :: component_mass(:,:)
    type(spray_parcel), allocatable :: incoming(:)
    real(dp), allocatable :: masses(:,:),send_real(:,:),recv_real(:,:)
    integer(int64), allocatable :: send_int(:,:),recv_int(:,:)
    integer, allocatable :: received(:),displacements(:)
    integer :: nranks,ierr,nc,nreal,nrecv,p,r,low,high,ios
    integer(int64) :: cursor
    logical :: valid,global_valid
    ok=.false.
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=allocated(parcels).and.maximum_parcels>=0.and.maximum_parcels<=1000000
    if (allocated(parcels)) valid=valid.and.size(parcels)<=1000000
    valid=valid.and.size(counts)==nranks.and.size(offsets)==nranks
    nc=0
    if (present(component_mass)) then
      valid=valid.and.allocated(component_mass)
      if (allocated(component_mass)) then
        nc=size(component_mass,1)
        if (allocated(parcels)) valid=valid.and.size(component_mass,2)==size(parcels)
      end if
      valid=valid.and.nc>=1.and.nc<=16
    end if
    call agree()
    if (.not.global_valid) return
    call MPI_Allreduce(nc,low,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(nc,high,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=low==high;cursor=0
    do r=1,nranks
      if (counts(r)<0.or.int(offsets(r),int64)/=cursor) valid=.false.
      cursor=cursor+int(counts(r),int64)
    end do
    valid=valid.and.cursor==size(parcels)
    do p=1,size(parcels)
      if (.not.valid_spray_parcel(parcels(p))) valid=.false.
      if (present(component_mass)) then
        if (.not.valid_component_inventory(parcels(p),component_mass(:,p))) valid=.false.
      end if
    end do
    call agree()
    if (.not.global_valid) return
    allocate(received(nranks),displacements(nranks),stat=ios)
    valid=ios==0;call agree()
    if (.not.global_valid) return
    call MPI_Alltoall(counts,1,MPI_INTEGER,received,1,MPI_INTEGER,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=sum(int(received,int64))<=int(maximum_parcels,int64)
    call agree()
    if (.not.global_valid) return
    nrecv=sum(received);displacements(1)=0
    do r=2,nranks
      displacements(r)=displacements(r-1)+received(r-1)
    end do
    nreal=14+nc
    allocate(send_int(3,size(parcels)),send_real(nreal,size(parcels)),recv_int(3,nrecv),recv_real(nreal,nrecv),stat=ios)
    valid=ios==0;call agree()
    if (.not.global_valid) return
    do p=1,size(parcels)
      send_int(:,p)=[parcels(p)%id,parcels(p)%parent_id,int(parcels(p)%material_id,int64)]
      send_real(1:14,p)=[parcels(p)%position,parcels(p)%velocity,parcels(p)%mass,parcels(p)%temperature, &
        parcels(p)%multiplicity,parcels(p)%deformation,parcels(p)%deformation_rate,parcels(p)%rt_time, &
        parcels(p)%khrt_initial_count,parcels(p)%khrt_shed_volume]
      if (nc>0) send_real(15:,p)=component_mass(:,p)
    end do
    call MPI_Alltoallv(send_int,3*counts,3*offsets,MPI_INTEGER8,recv_int,3*received,3*displacements,MPI_INTEGER8,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Alltoallv(send_real,nreal*counts,nreal*offsets,MPI_DOUBLE_PRECISION, &
      recv_real,nreal*received,nreal*displacements,MPI_DOUBLE_PRECISION,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (nrecv>1) call sort_records(1,nrecv)
    allocate(incoming(nrecv),stat=ios)
    valid=ios==0;call agree()
    if (.not.global_valid) return
    if (allocated(masses)) deallocate(masses)
    if (nc>0) allocate(masses(nc,nrecv),stat=ios)
    valid=ios==0;call agree()
    if (.not.global_valid) return
    valid=.true.
    do p=1,nrecv
      incoming(p)=spray_parcel(recv_int(1,p),recv_real(1:3,p),recv_real(4:6,p),recv_real(7,p), &
        recv_real(8,p),recv_real(9,p),recv_real(10,p),recv_real(11,p),recv_int(2,p),recv_real(12,p), &
        int(recv_int(3,p)),recv_real(13,p),recv_real(14,p))
      valid=valid.and.valid_spray_parcel(incoming(p))
      if (p>1) valid=valid.and.incoming(p)%id>incoming(p-1)%id
      if (nc>0) then
        masses(:,p)=recv_real(15:,p)
        valid=valid.and.valid_component_inventory(incoming(p),masses(:,p))
      end if
    end do
    call agree()
    if (.not.global_valid) return
    call move_alloc(incoming,parcels)
    if (present(component_mass)) call move_alloc(masses,component_mass)
    ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,global_valid,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) global_valid=.false.
    end subroutine
    recursive subroutine sort_records(first,last)
      integer, intent(in) :: first,last
      integer :: left,right
      integer(int64) :: pivot,tmp_int(3)
      real(dp) :: tmp_real(nreal)
      left=first;right=last;pivot=recv_int(1,first+(last-first)/2)
      do
        do while(recv_int(1,left)<pivot)
          left=left+1
        end do
        do while(recv_int(1,right)>pivot)
          right=right-1
        end do
        if (left<=right) then
          tmp_int=recv_int(:,left);recv_int(:,left)=recv_int(:,right);recv_int(:,right)=tmp_int
          tmp_real=recv_real(:,left);recv_real(:,left)=recv_real(:,right);recv_real(:,right)=tmp_real
          left=left+1;right=right-1
        end if
        if (left>right) exit
      end do
      if (first<right) call sort_records(first,right)
      if (left<last) call sort_records(left,last)
    end subroutine
  end subroutine
end module
