program test_spray_mpi_amr_reactive
  use mpi_f08
  use spray_mpi_amr_reactive_mod, only: advance_mpi_amr_spray_reactive
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use spray_parcel_mod
  use spray_mixture_mod, only: liquid_mixture_properties
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_full_transport
  use, intrinsic :: iso_fortran_env, only: int64
  use spray_mpi_amr_step_mod, only: advance_mpi_amr_mixture_coupling
  use spray_amr_tree_partition_mod, only: locate_amr_tree_spray_cells
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
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp), allocatable :: boundary(:)
  character(len=32) :: faces(6)
  real(dp) :: correction,theta
  real(dp) :: sgs(4),accel(3)=0
  real(dp), allocatable :: source(:)
  character(len=16) :: argument
  integer :: mode
  character(len=8) :: solver
  logical :: recursive,chemistry
  type(spray_parcel), allocatable :: drops(:),saved_drops(:)
  real(dp), allocatable :: masses(:,:),saved_mass(:,:)
  real(dp), allocatable :: root(:,:,:,:),primitive(:),cell(:),y(:),before(:),after(:),scale(:)
  real(dp), allocatable :: positions(:,:)
  integer, allocatable :: locations(:,:)
  real(dp) :: spacing(3)=.001_dp,zero(3)=0._dp,t,sound,r,interval
  integer :: nv,node,n(3),i,j,k,p,water,nitrogen,steps
  integer :: rank,nranks,ierr,owners(3),q
  logical :: ok
  integer(int64) :: last_id,created
  call get_command_argument(1,argument)
  sgs=[0._dp,.7_dp,.7_dp,0._dp]
  if (argument=='les') sgs=[.17_dp,.7_dp,.7_dp,.1_dp]
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks==2.or.nranks==4,'ranks')
  call load_h2o2_full_thermo(species,ok)
  call require(ok,'thermodynamics')
  nv=reactive_nvar(size(species))
  allocate(root(nv,4,4,4),primitive(nv),cell(nv),y(size(species)),before(nv),after(nv),scale(nv))
  allocate(boundary(nv));boundary=123;faces='periodic'
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
  ! All parcels begin on the last rank, close to coarse/fine and periodic interfaces.
  allocate(drops(merge(3,0,rank==nranks-1)),masses(2,merge(3,0,rank==nranks-1)))
  do p=1,size(drops)
    drops(p)=spray_parcel(int(p,int64),spacing*[1.99_dp,.8_dp,.8_dp], &
      [20._dp,10._dp,5._dp],1e-12_dp,350._dp,100._dp)
    if (p==2) drops(p)%position=spacing*[.99_dp,.9_dp,.9_dp]
    if (p==3) drops(p)%position=spacing*[3.99_dp,1.1_dp,1.1_dp]
    masses(:,p)=drops(p)%mass*[.25_dp,.75_dp]
  end do
  call load_h2o2_full_mechanism(reactions,ok)
  call require(ok,'reaction table')
  saved=local;saved_t=local_t;saved_drops=drops;saved_mass=masses
  last_id=3;interval=1.e-8_dp;correction=123;theta=456
  solver='invalid';recursive=.false.;chemistry=.false.
  if (argument=='les') then
    solver='hllc'
    if (rank==nranks-1) sgs(4)=.2_dp
    call coupled()
    call require(.not.ok.and.steps==0.and.correction==123.and.theta==456,'LES disagreement preserves diagnostics')
    call unchanged()
    sgs(4)=.1_dp
    solver='invalid'
  end if
  call coupled()
  call require(.not.ok.and.steps==0.and.correction==123.and.theta==456,'gas-stage failure preserves diagnostics')
  call require(all(boundary==123),'gas failure preserves boundary ledger')
  call unchanged()
  solver='hllc'
  do mode=1,2
    local=saved;local_t=saved_t;drops=saved_drops;masses=saved_mass
    recursive=mode==2
    call totals(before)
    call coupled()
    call require(ok.and.steps>=2.and.theta>=0.and.theta<=1,'coupled common/recursive gas step')
    call totals(after)
    scale=max(abs(before),abs(after),1.e-30_dp)
    call require(all(abs(after-before)<=3.e-12_dp*scale),'coupled composite conservation')
  end do
  chemistry=.true.
  call coupled()
  call require(ok,'reacting coupled continuation')
  call totals(after)
  call require(all(abs(after(:5)-before(:5))<=3.e-12_dp*scale(:5)),'reacting mass momentum energy conservation')
  ! Open gas/particle boundary and specular wall, with gas transport enabled.
  chemistry=.false.;faces(1)='slip_wall';faces(2)='outflow'
  do mode=1,2
    local=saved;local_t=saved_t;drops=saved_drops;masses=saved_mass
    recursive=mode==2
    do p=1,size(drops)
      drops(p)%velocity=[0._dp,10._dp,5._dp]
      if (p==1) then
        drops(p)%position(1)=3.999_dp*spacing(1);drops(p)%velocity(1)=1000
      else if (p==3) then
        drops(p)%position(1)=.001_dp*spacing(1);drops(p)%velocity(1)=-1000
      end if
    end do
    call totals(before)
    call coupled()
    call require(ok,'physical boundary common/recursive coupled step')
    call totals(after)
    scale=max(abs(before),abs(after),abs(boundary),1.e-30_dp)
    call require(all(abs(after-before-boundary)<=3.e-12_dp*scale),'combined gas and parcel boundary budget')
    q=size(drops)
    call MPI_Allreduce(MPI_IN_PLACE,q,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(q==2,'coupled outflow removes parcel')
  end do
  if (argument=='split') then
    local=saved;local_t=saved_t;drops=saved_drops;masses=saved_mass
    faces='periodic';recursive=.false.;chemistry=.false.;last_id=3
    options%breakup_model='etab';options%surface_tension=1e3_dp;options%liquid_viscosity=.001_dp
    options%maximum_multiplicity=50;options%maximum_parcels=3
    call coupled()
    call require(.not.ok.and.created==0,'split capacity rejects coupled candidate')
    call unchanged()
    options%maximum_parcels=100
    if (rank==nranks-1) options%maximum_multiplicity=60
    call coupled()
    call require(.not.ok.and.created==0,'split controls agree collectively')
    call unchanged()
    options%maximum_multiplicity=50
    call totals(before)
    call coupled()
    call require(ok.and.created==3.and.last_id==6,'split reserves three global IDs')
    call totals(after)
    scale=max(abs(before),abs(after),1.e-30_dp)
    call require(all(abs(after-before)<=3.e-12_dp*scale),'split composite conservation')
    call require(all(drops%multiplicity<=50),'split multiplicity bound')
    q=size(drops)
    call MPI_Allreduce(MPI_IN_PLACE,q,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(q==6,'split doubles three parcels')
  end if
  if (argument=='les') then
    local=saved;local_t=saved_t;drops=saved_drops;masses=saved_mass
    faces='periodic';recursive=rank==nranks-1;chemistry=.false.;last_id=3
    call coupled()
    call require(.not.ok,'rank-disagreeing subcycling rejects before LES dispatch')
    call unchanged()
  end if
  if (argument=='force') then
    allocate(source(nv));source=123
    local=saved;local_t=saved_t;drops=saved_drops;masses=saved_mass
    faces='periodic';recursive=.false.;chemistry=.false.;last_id=3
    accel=[100._dp,-200._dp,300._dp]
    if (rank==nranks-1) accel(1)=101
    call coupled()
    call require(.not.ok.and.all(source==123),'force disagreement preserves source ledger')
    call unchanged()
    accel(1)=100
    solver='invalid'
    call coupled()
    call require(.not.ok.and.all(source==123),'later failure preserves force ledger')
    call unchanged()
    solver='hllc'
    call totals(before)
    call coupled()
    call require(ok,'uniform acceleration coupled step')
    call totals(after)
    scale=max(abs(before),abs(after),abs(source),1.e-30_dp)
    call require(all(abs(after-before-source)<=3.e-12_dp*scale),'external impulse and work balance')
    call require(all(abs(source(2:4)-before(1)*interval*accel)<= &
      3.e-12_dp*abs(before(1)*interval*accel)),'analytic total impulse')
  end if
  if (rank==0) print *, 'PASS: MPI AMR reacting spray periodic and physical boundaries'
  call MPI_Finalize(ierr)
contains
  subroutine coupled()
    call advance_mpi_amr_spray_reactive(MPI_COMM_WORLD,tree,owners,species,reactions,transport,liquids,options, &
      local,local_t,drops,masses,zero,spacing,interval,'pcm','minmod',solver,chemistry, &
      .true.,.true.,.true.,.false.,1.e-6_dp,1.e-12_dp,steps,correction,theta,ok,last_id, &
      chemistry_integrator='implicit',subcycling=recursive,boundary_faces=faces,boundary_change=boundary, &
      les_settings=sgs,les_model='vreman',acceleration=accel,source_change=source,split_created=created)
  end subroutine
  subroutine unchanged()
    integer :: a
    call require(size(drops)==size(saved_drops).and.last_id==3,'particle count and IDs rollback')
    call require(all(masses==saved_mass),'component rollback')
    do a=1,size(drops)
      call require(all(drops(a)%position==saved_drops(a)%position).and. &
        all(drops(a)%velocity==saved_drops(a)%velocity).and.drops(a)%mass==saved_drops(a)%mass.and. &
        drops(a)%temperature==saved_drops(a)%temperature,'particle rollback')
    end do
    do a=1,3
      if (owners(a)/=rank) cycle
      call require(all(local%fields(a)%values==saved%fields(a)%values),'gas rollback')
      call require(all(local_t%fields(a)%values==saved_t%fields(a)%values),'temperature rollback')
    end do
  end subroutine
  subroutine totals(total)
    real(dp), intent(out) :: total(:)
    type(amr_tree_state_3d) :: assembled
    real(dp) :: liquid_total(nv),rho,cp,e,fractions(2),amount
    integer :: a,b,c
    assembled=gas
    do a=1,3
      assembled%fields(a)%values=0
      if (owners(a)==rank) assembled%fields(a)%values=local%fields(a)%values
      call MPI_Allreduce(MPI_IN_PLACE,assembled%fields(a)%values,size(assembled%fields(a)%values), &
        MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    end do
    call composite_integrals_amr_patch_tree_3d(tree,assembled,spacing,total,ok)
    call require(ok,'gas composite integral')
    liquid_total=0
    do a=1,size(drops)
      fractions=masses(:,a)/drops(a)%mass
      call liquid_mixture_properties(liquids,fractions,drops(a)%temperature,rho,cp,e,ok)
      call require(ok,'liquid energy')
      amount=drops(a)%mass*drops(a)%multiplicity
      liquid_total(1)=liquid_total(1)+amount
      liquid_total(2:4)=liquid_total(2:4)+amount*drops(a)%velocity
      liquid_total(5)=liquid_total(5)+amount*(e+.5_dp*sum(drops(a)%velocity**2))
      do b=1,2
        c=5+liquids(b)%vapor_index
        liquid_total(c)=liquid_total(c)+masses(b,a)*drops(a)%multiplicity
      end do
    end do
    call MPI_Allreduce(MPI_IN_PLACE,liquid_total,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    total=total+liquid_total
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
