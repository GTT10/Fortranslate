program test_spray_mpi_eb_film_reactive
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d,coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use spray_mpi_eb_reactive_mod, only: advance_mpi_eb_spray_reactive
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use spray_parcel_mod
  use spray_mixture_mod, only: liquid_mixture_properties
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_full_transport
  use, intrinsic :: iso_fortran_env, only: int64
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options
  implicit none
  type(spray_wall_film), allocatable :: films(:),saved_films(:)
  real(dp), allocatable :: film_mass(:,:),saved_film_mass(:,:)
  type(spray_wall_options) :: walls
  integer :: film_steps
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp), allocatable :: boundary(:),chemical(:)
  real(dp) :: theta
  character(len=8) :: solver
  type(spray_liquid) :: liquids(2)
  type(spray_options) :: options
  type(spray_parcel), allocatable :: drops(:),saved_drops(:)
  real(dp), allocatable :: masses(:,:),saved_mass(:,:),before(:),after(:),ledger(:)
  integer :: steps,p
  real(dp) :: rho,cp,e,amount
  type(amr_patch_3d) :: patch(1)
  type(amr_patch_tree_3d) :: tree
  type(eb_geometry_3d) :: finest,g(2)
  type(amr_tree_state_3d) :: field,temp,saved,saved_t
  type(nasa7_species), allocatable :: species(:)
  real(dp), allocatable :: positions(:,:),feedback(:,:),cell(:),primitive(:),y(:),total(:),parent(:),fine(:)
  real(dp) :: phi(0:8,0:8,0:8),zero(3)=0,spacing(3)=.001_dp,t,sound,r,volume
  integer :: rank,nranks,ierr,owners(2),dims(3),node,i,j,k,nv,nitrogen
  logical :: ok,include_films=.true.
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks==2.or.nranks==4,'rank count')
  owners=[0,nranks-1]
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patch(1),ok)
  call require(ok,'patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1],patch,tree,ok)
  call require(ok,'topology')
  do j=0,8
    phi(:,j,:)=real(j,dp)*.0005_dp-.00125_dp
  end do
  call build_level_set_eb_geometry_3d(phi,zero,[.004_dp,.004_dp,.004_dp],finest,ok)
  call require(ok,'fine geometry')
  call coarsen_eb_geometry_3d(finest,2,g(1),ok)
  call require(ok,'root geometry')
  call extract_eb_geometry_3d(finest,[1,1,1],[4,8,8],g(2),ok)
  call require(ok,'child geometry')
  call load_h2o2_full_thermo(species,ok)
  call require(ok,'thermo');nv=reactive_nvar(size(species))
  allocate(cell(nv),primitive(nv),y(size(species)),total(nv),parent(nv),fine(nv))
  nitrogen=0
  do i=1,size(species)
    if(species(i)%name=='N2') nitrogen=i
  end do
  call require(nitrogen>0,'species')
  y=0;y(nitrogen)=1;r=mixture_specific_gas_constant(species,y,ok)
  call require(ok,'gas constant')
  primitive=0;primitive(1)=1;primitive(2:4)=[1._dp,2._dp,3._dp];primitive(5)=r*1000;primitive(6:)=y
  call reactive_primitive_to_conserved(species,primitive,cell,t,sound,ok)
  call require(ok,'gas state')
  allocate(field%fields(2),temp%fields(2),positions(3,merge(1,0,rank==0)),feedback(nv,merge(1,0,rank==0)))
  do node=1,2
    if(owners(node)/=rank) cycle
    dims=tree%nodes(node)%shape
    allocate(field%fields(node)%values(nv,dims(1),dims(2),dims(3)),temp%fields(node)%values(1,dims(1),dims(2),dims(3)))
    do k=1,dims(3)
    do j=1,dims(2)
    do i=1,dims(1)
      field%fields(node)%values(:,i,j,k)=cell;temp%fields(node)%values(1,i,j,k)=t
      if(g(node)%volume_fraction(i,j,k)==0) then
        field%fields(node)%values(:,i,j,k)=0;temp%fields(node)%values(1,i,j,k)=42
      end if
    end do
    end do
    end do
  end do
  call load_h2o2_full_transport(transport,ok)
  call require(ok,'transport')
  call load_h2o2_full_mechanism(reactions,ok)
  call require(ok,'mechanism')
  allocate(boundary(nv),chemical(nv))
  call initialize_spray_liquid(species,'H2O',1000._dp,4184._dp,373.15_dp,101325._dp,2.257e6_dp,liquids(1),ok)
  call require(ok,'water')
  call initialize_spray_liquid(species,'H2O2',800._dp,2000._dp,350._dp,101325._dp,1e6_dp,liquids(2),ok)
  call require(ok,'peroxide')
  allocate(drops(merge(1,0,rank==nranks-1)),masses(2,merge(1,0,rank==nranks-1)),before(nv),after(nv),ledger(nv))
  if(rank==nranks-1) then
    drops(1)=spray_parcel(1_int64,[.00199_dp,.001250001_dp,.0015_dp], &
      [20._dp,-20._dp,5._dp],1e-12_dp,350._dp,100._dp)
    masses(:,1)=drops(1)%mass*[.25_dp,.75_dp]
  end if
  allocate(films(merge(2,0,rank==0)),film_mass(2,merge(2,0,rank==0)))
  do p=1,size(films)
    films(p)=spray_wall_film([.001_dp,.00125_dp,.0015_dp],[0._dp,1._dp,0._dp], &
      1e-10_dp,330._dp,.0001_dp,1e-6_dp,p)
    if(p==2) films(p)%position(1)=.003_dp
    film_mass(:,p)=films(p)%mass*[.25_dp,.75_dp]
  end do
  walls%enabled=.true.;walls%evolve=.true.;walls%temperature=400;walls%liquid_conductivity=.6_dp
  walls%minimum_height=1e-10_dp
  saved_films=films;saved_film_mass=film_mass
  saved=field;saved_t=temp;saved_drops=drops;saved_mass=masses
  ledger=42;boundary=43;chemical=44;theta=45;solver='invalid'
  include_films=rank/=0
  call advance()
  call require(.not.ok.and.steps==0.and.film_steps==0,'rank-disagreed film presence rejects')
  include_films=.true.
  call advance()
  call require(.not.ok.and.steps==0.and.all(ledger==42),'whole interval failure and ledger rollback')
  call require(all(masses==saved_mass).and.size(drops)==size(saved_drops),'particle population rollback')
  if(rank==nranks-1) then
    call require(all(drops(1)%position==saved_drops(1)%position),'particle position rollback')
    call require(all(drops(1)%velocity==saved_drops(1)%velocity),'particle velocity rollback')
    call require(drops(1)%temperature==saved_drops(1)%temperature,'particle temperature rollback')
    call require(drops(1)%mass==saved_drops(1)%mass,'particle mass rollback')
  end if
  do node=1,2
    if(owners(node)/=rank) cycle
    call require(all(field%fields(node)%values==saved%fields(node)%values),'gas rollback')
    call require(all(temp%fields(node)%values==saved_t%fields(node)%values),'temperature rollback')
  end do
  call require(all(boundary==43).and.all(chemical==44).and.theta==45,'gas diagnostics rollback')
  call require(film_steps==0.and.size(films)==size(saved_films),'film population rollback')
  call require(all(film_mass==saved_film_mass),'film component rollback')
  do p=1,size(films)
    call require(films(p)%mass==saved_films(p)%mass.and.films(p)%temperature==saved_films(p)%temperature, &
      'film physical rollback')
  end do
  solver='hllc'
  call inventory(before)
  call advance()
  call require(ok.and.steps>1,'distributed EB interval')
  call inventory(after)
  call require(all(abs(after-before-ledger-boundary-chemical)<=3e-12_dp*max(abs(before),abs(after),1e-30_dp)), &
    'composite gas parcel wall conservation')
  call require(any(chemical(6:)/=0),'chemical species exchange')
  call require(size(drops)==merge(1,0,rank==nranks-1),'final fine owner')
  if(rank==nranks-1) then
    call require(drops(1)%position(1)<.002_dp.and.drops(1)%position(2)>.00125_dp,'fine cell and fluid position')
    call require(drops(1)%temperature>350.and.drops(1)%mass<1e-12_dp,'heating and evaporation')
    call require(drops(1)%velocity(2)>0.and.drops(1)%velocity(1)<20,'reflection and drag')
  end if
  call require(film_steps>=2,'both film source halves')
  call require(size(films)==merge(1,0,rank==0.or.rank==nranks-1),'film finest ownership')
  do p=1,size(films)
    call require(films(p)%mass<1e-10_dp.and.films(p)%temperature>330,'film heat and evaporation')
  end do
  if(rank==0) print *, 'PASS: MPI EB reactive particle and film composition with rollback'
  call MPI_Finalize(ierr)
contains
  subroutine advance()
    if(include_films) then
    call advance_mpi_eb_spray_reactive(MPI_COMM_WORLD,tree,owners,g,phi,species,reactions,transport, &
      liquids,options,field,temp,drops,masses,zero,spacing,1e-8_dp,solver, &
      .true.,.true.,.true.,.true.,.true.,.true.,2e-7_dp,1e-12_dp,steps,theta,boundary,ledger,chemical,ok,'implicit', &
      films=films,film_mass=film_mass,film_settings=walls,film_steps=film_steps)
    else
    call advance_mpi_eb_spray_reactive(MPI_COMM_WORLD,tree,owners,g,phi,species,reactions,transport, &
      liquids,options,field,temp,drops,masses,zero,spacing,1e-8_dp,solver, &
      .true.,.true.,.true.,.true.,.true.,.true.,2e-7_dp,1e-12_dp,steps,theta,boundary,ledger,chemical,ok,'implicit', &
      film_steps=film_steps)
    end if
  end subroutine
  subroutine inventory(result)
    real(dp), intent(out) :: result(:)
    result=0
    do node=1,2
      if(owners(node)/=rank) cycle
      dims=tree%nodes(node)%shape
      do k=1,dims(3)
      do j=1,dims(2)
      do i=1,dims(1)
        if(node==1.and.i<=2) cycle
        volume=g(node)%volume_fraction(i,j,k)*g(node)%dx*g(node)%dy*g(node)%dz
        result=result+volume*field%fields(node)%values(:,i,j,k)
      end do
      end do
      end do
    end do
    do p=1,size(drops)
      call liquid_mixture_properties(liquids,masses(:,p)/drops(p)%mass,drops(p)%temperature,rho,cp,e,ok)
      call require(ok,'liquid properties')
      amount=drops(p)%mass*drops(p)%multiplicity
      result(1)=result(1)+amount;result(2:4)=result(2:4)+amount*drops(p)%velocity
      result(5)=result(5)+amount*(e+.5_dp*sum(drops(p)%velocity**2))
      do i=1,2
        result(5+liquids(i)%vapor_index)=result(5+liquids(i)%vapor_index)+masses(i,p)*drops(p)%multiplicity
      end do
    end do
    do p=1,size(films)
      call liquid_mixture_properties(liquids,film_mass(:,p)/films(p)%mass,films(p)%temperature,rho,cp,e,ok)
      call require(ok,'film properties')
      result(1)=result(1)+films(p)%mass;result(5)=result(5)+films(p)%mass*e
      do i=1,2
        result(5+liquids(i)%vapor_index)=result(5+liquids(i)%vapor_index)+film_mass(i,p)
      end do
    end do
    call MPI_Allreduce(MPI_IN_PLACE,result,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL: ',rank,': ',label
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
  end subroutine
end program
