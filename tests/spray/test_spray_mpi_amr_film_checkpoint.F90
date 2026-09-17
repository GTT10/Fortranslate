program test_spray_mpi_amr_film_checkpoint
  use mpi_f08
  use mpi_amr_patch_tree_thermo_3d_mod, only: recover_mpi_amr_tree_temperatures_3d
  use spray_mpi_amr_checkpoint_mod, only: mpi_amr_spray_checkpoint
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options
  use spray_amr_tree_partition_mod, only: locate_amr_tree_spray_cells
  use spray_parcel_mod
  use spray_mixture_mod, only: liquid_mixture_properties
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_full_transport
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: patches(2)
  type(amr_patch_tree_3d) :: tree
  type(amr_tree_state_3d) :: gas,temp,saved,saved_t,expected,expected_t,local,local_t
  type(gas_transport_species), allocatable :: transport(:)
  type(spray_liquid) :: liquids(2)
  type(spray_options) :: options
  character(len=32) :: faces(6)
  real(dp), allocatable :: boundary(:)
  type(spray_wall_film), allocatable :: drops(:),saved_drops(:)
  type(spray_wall_options) :: walls
  real(dp), allocatable :: masses(:,:),saved_mass(:,:)
  real(dp), allocatable :: root(:,:,:,:),primitive(:),cell(:),y(:),before(:),after(:),scale(:)
  real(dp), allocatable :: positions(:,:)
  integer, allocatable :: locations(:,:)
  real(dp) :: spacing(3)=.001_dp,zero(3)=0._dp,t,sound,r,interval
  integer :: nv,node,n(3),i,j,k,p,water,nitrogen,steps
  integer :: rank,nranks,ierr,q,unit,event,capacity
  integer, allocatable :: owners(:)
  character(len=1024) :: path,base
  integer(int64) :: tick
  logical :: regridded
  real(dp) :: time
  type(spray_parcel), allocatable :: parcels(:)
  real(dp), allocatable :: parcel_mass(:,:),injected(:)
  logical :: ok
  integer(int64) :: last_id
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks==2.or.nranks==4,'ranks')
  call load_h2o2_full_thermo(species,ok)
  call require(ok,'thermodynamics')
  nv=reactive_nvar(size(species))
  allocate(root(nv,4,4,4),primitive(nv),cell(nv),y(size(species)),before(nv),after(nv),scale(nv))
  water=0;nitrogen=0
  do p=1,size(species)
    if (species(p)%name=='H2O') water=p
    if (species(p)%name=='N2') nitrogen=p
  end do
  call require(water>0.and.nitrogen>0,'species identity')
  y=0._dp;y(nitrogen)=1._dp
  r=mixture_specific_gas_constant(species,y,ok)
  call require(ok,'gas constant')
  primitive=0._dp;primitive(1)=1._dp;primitive(5)=r*1000._dp;primitive(6:)=y
  call reactive_primitive_to_conserved(species,primitive,cell,t,sound,ok)
  call require(ok,'initial gas')
  do k=1,4
  do j=1,4
  do i=1,4
    root(:,i,j,k)=cell
  end do
  end do
  end do
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patches(1),ok)
  call require(ok,'fine patch')
  call initialize_amr_patch_3d(4,8,8,1,2,1,2,1,2,2,patches(2),ok)
  call require(ok,'nested patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,tree,ok)
  call require(ok,'tree')
  call prolong_amr_patch_tree_3d(tree,root,gas,ok)
  call require(ok,'gas fields')
  allocate(temp%fields(3))
  do node=1,3
    n=tree%nodes(node)%shape
    allocate(temp%fields(node)%values(1,n(1),n(2),n(3)))
    temp%fields(node)%values=1000._dp
  end do
  call load_h2o2_full_transport(transport,ok)
  call require(ok,'transport')
  call initialize_spray_liquid(species,'H2O',1000._dp,4184._dp,373.15_dp,101325._dp,2.257e6_dp,liquids(1),ok)
  call require(ok,'water liquid')
  call initialize_spray_liquid(species,'H2O2',800._dp,2000._dp,350._dp,101325._dp,1e6_dp,liquids(2),ok)
  call require(ok,'peroxide liquid')
  owners=[0,1,nranks-1]
  allocate(local%fields(3),local_t%fields(3))
  do node=1,3
    if (owners(node)/=rank) cycle
    local%fields(node)%values=gas%fields(node)%values
    local_t%fields(node)%values=temp%fields(node)%values
  end do
  allocate(drops(merge(2,0,rank==nranks-1)),masses(2,merge(2,0,rank==nranks-1)))
  do p=1,size(drops)
    drops(p)=spray_wall_film([0._dp,.0005_dp,.0005_dp],[1._dp,0._dp,0._dp], &
      1e-9_dp,330._dp,.001_dp,1e-6_dp,p)
    if (p==2) then
      drops(p)%position(1)=.004_dp;drops(p)%normal(1)=-1
    end if
    masses(:,p)=drops(p)%mass*[.25_dp,.75_dp]
  end do
  call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,local,local_t,MPI_COMM_WORLD,ok)
  call require(ok,'source EOS recovery')
  do node=1,size(tree%nodes)
    if (rank==owners(node)) temp%fields(node)%values=local_t%fields(node)%values
    call MPI_Bcast(temp%fields(node)%values,size(temp%fields(node)%values),MPI_DOUBLE_PRECISION, &
      owners(node),MPI_COMM_WORLD,ierr)
  end do
  call get_command_argument(1,base)
  call require(len_trim(base)>0,'checkpoint artifact prefix')
  path=''
  if (rank==0) then
    call system_clock(tick)
    write(path,'(a,"_",i0)') trim(base),tick
  end if
  call MPI_Bcast(path,len(path),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  allocate(parcels(0),parcel_mass(2,0),injected(nv),boundary(nv))
  injected=0;boundary=0;boundary(5)=.125_dp
  time=1e-8_dp;steps=2;event=3;regridded=.false.;last_id=17;capacity=100
  saved=local;saved_t=local_t;saved_drops=drops;saved_mass=masses
  call checkpoint(MPI_COMM_WORLD,.true.,.true.)
  call require(ok,'distributed film checkpoint write')
  call checkpoint(MPI_COMM_WORLD,.true.,.true.)
  call require(.not.ok,'immutable overwrite rejection');call unchanged()
  capacity=0;time=99;boundary=123
  call checkpoint(MPI_COMM_WORLD,.false.,.true.)
  call require(.not.ok.and.time==99.and.all(boundary==123),'late film capacity read preserves metadata')
  call unchanged()
  capacity=100
  call checkpoint(MPI_COMM_WORLD,.false.,rank/=nranks-1)
  call require(.not.ok.and.time==99,'rank-disagreed film presence rejects');call unchanged()
  call checkpoint(MPI_COMM_WORLD,.false.,.false.)
  call require(.not.ok.and.time==99,'omitted film mode rejects');call unchanged()
  ! Each rank independently reads the shared checkpoint as a one-rank run.
  call checkpoint(MPI_COMM_SELF,.false.,.true.)
  call require(ok.and.all(owners==0).and.size(drops)==2,'changed-rank film restoration')
  call require(time==1e-8_dp.and.steps==2.and.event==3.and.last_id==17,'restart metadata')
  call require(boundary(5)==.125_dp.and.count(boundary/=0)==1,'cumulative boundary restoration')
  do node=1,size(tree%nodes)
    call require(all(local%fields(node)%values==gas%fields(node)%values),'exact restored gas')
    call require(all(local_t%fields(node)%values==temp%fields(node)%values),'exact restored temperature')
  end do
  do p=1,2
    call require(drops(p)%mass==1e-9_dp.and.drops(p)%temperature==330.and. &
      drops(p)%diameter==.001_dp.and.drops(p)%height==1e-6_dp.and.drops(p)%material_id==p,'film physical payload')
    call require(drops(p)%position(1)==merge(0._dp,.004_dp,p==1),'exact lower/upper wall positions')
    call require(drops(p)%normal(1)==merge(1._dp,-1._dp,p==1),'film normals')
    call require(all(masses(:,p)==1e-9_dp*[.25_dp,.75_dp]),'film component payload')
  end do
  call checkpoint(MPI_COMM_WORLD,.false.,.true.)
  call require(ok,'restore distributed owners and film payload')
  saved=local;saved_t=local_t;saved_drops=drops;saved_mass=masses
  ! Trailing data must not be silently accepted after the film records/footer.
  if (rank==0) then
    open(newunit=unit,file=trim(path)//'.particles',status='old',position='append')
    write(unit,'(a)') 'TRAILING_CORRUPTION'
    close(unit)
  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  time=77;boundary=456
  call checkpoint(MPI_COMM_WORLD,.false.,.true.)
  call require(.not.ok.and.time==77.and.all(boundary==456),'corrupt film archive rejects atomically')
  call unchanged()
  if (rank==0) print *, 'PASS: MPI AMR film checkpoint; artifacts ',trim(path)
  call MPI_Finalize(ierr)
contains
  subroutine checkpoint(comm,writing,with_films)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(in) :: writing,with_films
    if (with_films) then
      call mpi_amr_spray_checkpoint(comm,writing,trim(path),'film-test',species,tree,owners,local,local_t, &
        parcels,parcel_mass,time,steps,event,regridded,last_id,injected,zero,spacing,100,ok, &
        [.false.,.true.,.true.],boundary,drops,masses,capacity)
    else
      call mpi_amr_spray_checkpoint(comm,writing,trim(path),'film-test',species,tree,owners,local,local_t, &
        parcels,parcel_mass,time,steps,event,regridded,last_id,injected,zero,spacing,100,ok, &
        [.false.,.true.,.true.],boundary)
    end if
  end subroutine
  subroutine unchanged()
    integer :: a
    call require(size(drops)==size(saved_drops),'film count rollback')
    call require(all(masses==saved_mass),'component rollback')
    do a=1,size(drops)
      call require(all(drops(a)%position==saved_drops(a)%position).and. &
        drops(a)%diameter==saved_drops(a)%diameter.and.drops(a)%height==saved_drops(a)%height.and. &
        drops(a)%mass==saved_drops(a)%mass.and. &
        drops(a)%temperature==saved_drops(a)%temperature,'particle rollback')
    end do
    do a=1,3
      if (owners(a)/=rank) cycle
      call require(all(local%fields(a)%values==saved%fields(a)%values),'gas rollback')
      call require(all(local_t%fields(a)%values==saved_t%fields(a)%values),'temperature rollback')
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
