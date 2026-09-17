program test_spray_mpi_eb_impact
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use spray_parcel_mod
  use spray_wall_impact_mod, only: spray_wall_film,wall_film_totals
  use spray_eb_impact_mod, only: spray_eb_wall_options
  use spray_mpi_eb_impact_mod, only: impact_mpi_level_set
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(spray_liquid) :: liquids(1)
  type(spray_options) :: options
  type(spray_eb_wall_options) :: walls
  type(spray_parcel), allocatable :: drops(:),saved(:)
  type(spray_wall_film), allocatable :: films(:),saved_films(:)
  real(dp), allocatable :: masses(:,:),film_mass(:,:),fractions(:),previous(:,:),pressure(:)
  logical, allocatable :: hit(:)
  real(dp) :: phi(0:2,0:2,0:2),normal(3),mw(1),wall(5),before(5),after(5),d
  integer(int64) :: last_id
  integer :: rank,nranks,ierr,i,j,k,p
  logical :: ok,rejected
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks==2.or.nranks==4,'ranks')
  call load_h2o2_full_thermo(species,ok)
  call require(ok,'thermo')
  call initialize_spray_liquid(species,'H2O',1000._dp,4184._dp,373.15_dp,101325._dp,2.257e6_dp,liquids(1),ok)
  call require(ok,'liquid')
  mw=species(liquids(1)%vapor_index)%molecular_weight
  do k=0,2
  do j=0,2
  do i=0,2
    phi(i,j,k)=.5_dp*(i+j)-1
  end do
  end do
  end do
  normal=[1._dp,1._dp,0._dp]/sqrt(2._dp)
  allocate(drops(1),masses(1,1),previous(3,1),pressure(1),hit(1),fractions(1))
  drops(1)=spray_parcel(int(rank+1,int64),[.4_dp,.4_dp,.5_dp],-200*normal,1e-12_dp,350._dp,10._dp)
  masses=1e-12_dp;previous(:,1)=[.6_dp,.6_dp,.5_dp];pressure=101325
  d=parcel_diameter(drops(1),liquids(1))
  allocate(films(merge(1,0,rank==nranks-1)),film_mass(1,merge(1,0,rank==nranks-1)))
  if(rank==nranks-1) then
    films(1)=spray_wall_film([.5_dp,.5_dp,.5_dp],normal,1e-10_dp,350._dp,.001_dp,.03_dp*d,1)
    film_mass=1e-10_dp
  end if
  saved=drops;saved_films=films;last_id=nranks
  options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp;options%maximum_parcels=1
  wall=42;fractions=43;hit=.false.
  call impact()
  call require(.not.ok.and..not.rejected,'capacity rejection')
  call unchanged()
  options%maximum_parcels=8
  if(rank==0) phi(0,0,0)=-2
  call impact()
  call require(.not.ok,'geometry rank disagreement');call unchanged()
  phi(0,0,0)=-1
  if(rank==nranks-1) films(1)%height=.2_dp*d
  call impact()
  call require(.not.ok.and.rejected,'remote thick-film model rejection')
  if(rank==nranks-1) films(1)=saved_films(1)
  call unchanged()
  call inventory(before)
  call impact()
  call require(ok.and..not.rejected.and.all(hit),'distributed oblique contact')
  call require(abs(fractions(1)-.5_dp)<1e-14_dp,'contact fraction')
  call require(size(drops)==5.and.drops(1)%mass==0,'four wet-splash children')
  call require(last_id==5*nranks,'global watermark')
  do p=2,5
    call require(drops(p)%id==nranks+4*rank+p-1,'rank-prefix child IDs')
    call require(dot_product(drops(p)%velocity,normal)>=0,'outgoing children')
  end do
  call inventory(after)
  call require(all(abs(after+wall-before)<=3e-12_dp*max(abs(before),1e-20_dp)),'global inventory')
  deallocate(drops,masses,films,film_mass,previous,pressure)
  allocate(drops(merge(1,0,rank==nranks-1)),masses(1,merge(1,0,rank==nranks-1)))
  allocate(previous(3,size(drops)),pressure(size(drops)),films(0),film_mass(1,0))
  if(rank==nranks-1) then
    drops=saved;masses=1e-12_dp;previous(:,1)=[.6_dp,.6_dp,.5_dp];pressure=101325
  end if
  last_id=nranks
  call impact()
  call require(ok.and.last_id==nranks+4,'empty ranks participate in ID reservation')
  call require(size(drops)==merge(5,0,rank==nranks-1),'empty rank payload')
  if(rank==0) print *, 'PASS: distributed EB contact, remote films, IDs and rollback'
  call MPI_Finalize(ierr)
contains
  subroutine impact()
    call impact_mpi_level_set(MPI_COMM_WORLD,liquids,options,walls,phi,mw,pressure,previous, &
      [0._dp,0._dp,0._dp],[1._dp,1._dp,1._dp],drops,masses,films,film_mass,last_id,hit,fractions,wall,ok,rejected)
  end subroutine
  subroutine inventory(value)
    real(dp), intent(out) :: value(5)
    integer :: q
    value=0
    do q=1,size(drops)
      value=value+parcel_totals(drops(q),liquids(1))
      call require(masses(1,q)==drops(q)%mass,'parcel component inventory')
    end do
    do q=1,size(films)
      value=value+wall_film_totals(films(q),liquids(1))
      call require(film_mass(1,q)==films(q)%mass,'film component inventory')
    end do
    call MPI_Allreduce(MPI_IN_PLACE,value,5,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  end subroutine
  subroutine unchanged()
    call require(size(drops)==1.and.all(masses==1e-12_dp).and.last_id==nranks,'population rollback')
    call require(all(drops(1)%position==saved(1)%position).and.drops(1)%mass==saved(1)%mass,'state rollback')
    call require(size(films)==size(saved_films).and.all(wall==42).and.all(fractions==43),'film and ledger rollback')
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL rank ',rank,': ',label
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
  end subroutine
end program
