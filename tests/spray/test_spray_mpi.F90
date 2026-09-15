program test_spray_mpi
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use spray_parcel_mod
  use spray_mpi_mod
  implicit none
  type(spray_parcel), allocatable :: parcels(:),saved(:)
  real(dp), allocatable :: masses(:,:),saved_mass(:,:)
  integer :: rank,nranks,ierr,p,target
  integer(int64) :: last_id,first_id,saved_last
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  allocate(parcels(2),masses(2,2))
  do p=1,2
    parcels(p)=spray_parcel(int(2*rank+p,int64),[.25_dp,.25_dp,.25_dp], &
      [1.0_dp,2.0_dp,3.0_dp],1e-9_dp,350.0_dp,100.0_dp)
    parcels(p)%khrt_initial_count=100;parcels(p)%khrt_shed_volume=real(parcels(p)%id,dp)*1e-15_dp
    masses(:,p)=parcels(p)%mass*[.25_dp,.75_dp]
  end do
  saved=parcels;saved_mass=masses
  call migrate_spray_parcels(MPI_COMM_WORLD,parcels,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp], &
    [nranks,1,1],[nranks,1,1],[.true.,.false.,.false.],2,ok,masses)
  call require(.not.ok.and.all(parcels%id==saved%id).and.all(masses==saved_mass),'capacity collective rollback')
  call migrate_spray_parcels(MPI_COMM_WORLD,parcels,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp], &
    [nranks,1,1],[nranks,1,1],[.true.,.false.,.false.],100,ok,masses)
  call require(ok,'gather migration')
  call require(size(parcels)==merge(2*nranks,0,rank==0),'empty destinations')
  do p=1,size(parcels)
    call require(parcels(p)%id==p,'received ID ordering')
    call require(parcels(p)%khrt_shed_volume==real(p,dp)*1e-15_dp,'KH payload')
    target=modulo(p-1,nranks);parcels(p)%position(1)=target+.25_dp
  end do
  call migrate_spray_parcels(MPI_COMM_WORLD,parcels,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp], &
    [nranks,1,1],[nranks,1,1],[.true.,.false.,.false.],100,ok,masses)
  call require(ok.and.size(parcels)==2,'redistribute from empty sources')
  do p=1,2
    call require(parcels(p)%id==rank+1+(p-1)*nranks,'destination ID ordering')
    call require(all(masses(:,p)==parcels(p)%mass*[.25_dp,.75_dp]),'component payload')
    call require(all(parcels(p)%velocity==[1.0_dp,2.0_dp,3.0_dp]),'velocity payload')
  end do
  saved=parcels;saved_mass=masses
  if (rank==0) masses(1,1)=-1
  call migrate_spray_parcels(MPI_COMM_WORLD,parcels,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp], &
    [nranks,1,1],[nranks,1,1],[.true.,.false.,.false.],100,ok,masses)
  call require(.not.ok.and.all(parcels%id==saved%id),'one-rank invalid input rejects collectively')
  if (rank/=0) call require(all(masses==saved_mass),'other ranks preserve mass')
  masses=saved_mass;last_id=int(2*nranks,int64)
  call split_mpi_spray_parcels(MPI_COMM_WORLD,parcels,60.0_dp,merge(2,100,rank==0),last_id,ok,masses)
  call require(.not.ok.and.size(parcels)==2.and.last_id==2*nranks,'split capacity collective rollback')
  call require(all(parcels%id==saved%id).and.all(masses==saved_mass),'split failure preserves payload')
  call split_mpi_spray_parcels(MPI_COMM_WORLD,parcels,60.0_dp,100,last_id,ok,masses)
  call require(ok.and.size(parcels)==4.and.last_id==4*nranks,'collective split ID allocation')
  call require(all(parcels(3:4)%id==[2*nranks+2*rank+1,2*nranks+2*rank+2]),'split globally disjoint IDs')
  call require(all(parcels(3:4)%parent_id==saved%id),'split lineage')
  call require(sum(parcels%multiplicity)==sum(saved%multiplicity),'split multiplicity conserved')
  call require(sum(parcels%khrt_initial_count)==sum(saved%khrt_initial_count),'split KH count conserved')
  call require(abs(sum(parcels%khrt_shed_volume)-sum(saved%khrt_shed_volume))<1e-28_dp,'split KH volume conserved')
  do p=1,4
    call require(all(masses(:,p)==parcels(p)%mass*[.25_dp,.75_dp]),'split component inheritance')
    call require(parcels(p)%multiplicity<=60,'split count cap')
  end do
  last_id=100;first_id=-7
  call reserve_spray_ids(MPI_COMM_WORLD,rank,last_id,first_id,ok)
  call require(ok.and.last_id==100+nranks*(nranks-1)/2,'collective ID high-water mark')
  if (rank==0) then
    call require(first_id==0,'zero request has no ID range')
  else
    call require(first_id==101+rank*(rank-1)/2,'disjoint ID prefix range')
  end if
  saved_last=last_id
  call reserve_spray_ids(MPI_COMM_WORLD,1,last_id,first_id,ok)
  call require(ok.and.first_id==saved_last+rank+1.and.last_id==saved_last+nranks,'successive ID reservations')
  saved_last=last_id;first_id=-9
  call reserve_spray_ids(MPI_COMM_WORLD,merge(-1,1,rank==0),last_id,first_id,ok)
  call require(.not.ok.and.last_id==saved_last.and.first_id==-9,'invalid ID request collective rollback')
  last_id=huge(last_id)-1;first_id=-11
  call reserve_spray_ids(MPI_COMM_WORLD,1,last_id,first_id,ok)
  call require(.not.ok.and.last_id==huge(last_id)-1.and.first_id==-11,'ID overflow rollback')
  last_id=int(rank,int64);first_id=-13
  call reserve_spray_ids(MPI_COMM_WORLD,1,last_id,first_id,ok)
  call require(.not.ok.and.last_id==rank.and.first_id==-13,'inconsistent ID history rejected')
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (condition) return
    print *,rank,': ',message
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  end subroutine
end program
