program test_spray_mpi
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use spray_parcel_mod
  use spray_mpi_mod
  use spray_partition_mod, only: spray_cartesian_block
  use spray_coupling_3d_mod, only: particle_boundary_stencil
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(spray_liquid) :: liquids(2)
  real(dp) :: fractions(2),ledger(15),direction(3)
  type(spray_parcel), allocatable :: parcels(:),saved(:)
  real(dp), allocatable :: masses(:,:),saved_mass(:,:)
  integer :: rank,nranks,ierr,p,target,first(3),dims(3),cells(3),i,j,k,n,indices(3,8)
  real(dp), allocatable :: field(:,:,:,:),positions(:,:),samples(:,:),saved_samples(:,:)
  real(dp) :: weights(8),expected_value
  character(len=16) :: faces(6)
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
  call load_h2o2_full_thermo(species,ok)
  call require(ok,'injection thermo')
  call initialize_spray_liquid(species,'H2O',1000.0_dp,4184.0_dp,373.15_dp,101325.0_dp,2.257e6_dp,liquids(1),ok)
  call require(ok,'injection water')
  call initialize_spray_liquid(species,'H2O2',800.0_dp,2000.0_dp,350.0_dp,101325.0_dp,1e6_dp,liquids(2),ok)
  call require(ok,'injection peroxide')
  deallocate(parcels,masses);allocate(parcels(0),masses(2,0));last_id=0
  fractions=[.2_dp+.1_dp*rank,.8_dp-.1_dp*rank]
  call inject_mpi_cone_parcels(MPI_COMM_WORLD,liquids(1),1e-10_dp*rank,rank, &
    [rank+.25_dp,.25_dp,.25_dp],[0.0_dp,0.0_dp,1.0_dp],10.0_dp,0.0_dp,2e-5_dp,330.0_dp, &
    100,last_id,parcels,ledger,ok,liquids,fractions,masses)
  call require(ok.and.size(parcels)==rank.and.last_id==nranks*(nranks-1)/2,'collective mixture injection')
  call require(abs(ledger(1)-1e-10_dp*rank)<1e-24_dp,'injected total mass')
  call require(abs(ledger(5+liquids(1)%vapor_index)-ledger(1)*fractions(1))<1e-24_dp,'injected component ledger')
  do p=1,size(parcels)
    call require(all(masses(:,p)==parcels(p)%mass*fractions),'injected composition')
    call require(parcels(p)%id==rank*(rank-1)/2+p,'injected unique IDs')
  end do
  saved=parcels;saved_mass=masses;saved_last=last_id;direction=[0.0_dp,0.0_dp,1.0_dp]
  if (rank==nranks-1) direction=0
  call inject_mpi_cone_parcels(MPI_COMM_WORLD,liquids(1),1e-10_dp,1, &
    [rank+.25_dp,.25_dp,.25_dp],direction,10.0_dp,0.0_dp,2e-5_dp,330.0_dp, &
    100,last_id,parcels,ledger,ok,liquids,fractions,masses)
  call require(.not.ok.and.last_id==saved_last.and.all(ledger==0),'injection failure restores global IDs and ledger')
  call require(all(parcels%id==saved%id).and.all(masses==saved_mass),'injection failure preserves local inventories')
  cells=[2*nranks+1,3,2];faces='slip_wall';faces(1:2)='periodic'
  call spray_cartesian_block(cells,[nranks,1,1],rank,first,dims,ok)
  call require(ok,'distributed field block')
  allocate(field(2,dims(1),dims(2),dims(3)),positions(3,3))
  do k=1,dims(3)
    do j=1,dims(2)
      do i=1,dims(1)
        field(1,i,j,k)=real(i+first(1)+10*j+100*k,dp);field(2,i,j,k)=-field(1,i,j,k)
      end do
    end do
  end do
  positions(:,1)=[2.9_dp,1.3_dp,.7_dp]
  positions(:,2)=[real(cells(1),dp)-.25_dp,2.5_dp,1.0_dp]
  positions(:,3)=[-.25_dp,.1_dp,.1_dp]
  call sample_mpi_spray_field(MPI_COMM_WORLD,field,positions,[0.0_dp,0.0_dp,0.0_dp], &
    [1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],faces,samples,ok)
  call require(ok,'distributed field sampling')
  do p=1,3
    call particle_boundary_stencil(positions(:,p),[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp], &
      cells,faces,indices,weights,ok)
    call require(ok,'serial reference stencil');expected_value=0
    do n=1,8
      expected_value=expected_value+weights(n)*real(indices(1,n)+10*indices(2,n)+100*indices(3,n),dp)
    end do
    call require(samples(1,p)==expected_value.and.samples(2,p)==-expected_value,'distributed versus serial interpolation')
  end do
  if (rank==0) positions=positions(:,:0)
  call sample_mpi_spray_field(MPI_COMM_WORLD,field,positions,[0.0_dp,0.0_dp,0.0_dp], &
    [1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],faces,samples,ok)
  call require(ok.and.size(samples,2)==merge(0,3,rank==0),'empty sampling source')
  saved_samples=samples
  if (rank==nranks-1) positions(2,1)=-1
  call sample_mpi_spray_field(MPI_COMM_WORLD,field,positions,[0.0_dp,0.0_dp,0.0_dp], &
    [1.0_dp,1.0_dp,1.0_dp],cells,[nranks,1,1],faces,samples,ok)
  call require(.not.ok.and.all(samples==saved_samples),'one-rank invalid sample preserves every result')
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
