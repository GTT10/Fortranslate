program test_spray_mpi_amr_film_step
  use mpi_f08
  use spray_mpi_amr_film_step_mod, only: advance_mpi_amr_film_coupling
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
  type(amr_tree_state_3d) :: exported
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
  integer :: rank,nranks,ierr,owners(3),q
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
  faces='periodic';faces(1:2)='slip_wall'
  walls%enabled(1:2)=.true.;walls%evolve=.true.;walls%temperature=400
  walls%liquid_conductivity=.6_dp;walls%minimum_height=1e-10_dp
  allocate(boundary(nv));boundary=123
  saved=local;saved_t=local_t;saved_drops=drops;saved_mass=masses
  exported=local
  do node=1,3
    if(allocated(exported%fields(node)%values)) exported%fields(node)%values=123
  end do
  interval=1e-3_dp;options%maximum_substeps=1
  call advance()
  call require(.not.ok.and.steps==0.and.all(boundary==123),'interval exhaustion preserves ledger')
  call unchanged()
  interval=1e-8_dp;options%maximum_substeps=100000
  call totals(before)
  call advance()
  call require(ok.and.steps>0.and.boundary(5)>0,'AMR film heating and wall energy')
  do node=1,3
    call require(allocated(exported%fields(node)%values).eqv.(owners(node)==rank),'source ownership')
    if(owners(node)/=rank) cycle
    call require(all(abs(local%fields(node)%values-saved%fields(node)%values-exported%fields(node)%values)<= &
      3e-13_dp*max(1._dp,abs(local%fields(node)%values))),'integrated film source reconstruction')
  end do
  call totals(after)
  scale=max(1e-12_dp,abs(before))
  call require(all(abs(after-before-boundary)<=2e-12_dp*scale),'gas film wall inventory')
  t=sum(drops%mass)
  call MPI_Allreduce(MPI_IN_PLACE,t,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call require(t<2e-9_dp,'evaporation transfers film mass to gas')
  do p=1,size(drops)
    call require(drops(p)%temperature>330,'film heats')
    call require(drops(p)%position(1)==0.or.drops(p)%position(1)==.004_dp,'exact wall positions retained')
  end do
  ! Empty populations still validate gas and return a zero exchange.
  deallocate(drops,masses);allocate(drops(0),masses(2,0))
  call advance()
  call require(ok.and.all(boundary==0),'globally empty films')
  do node=1,3
    if(owners(node)/=rank) cycle
    call require(all(exported%fields(node)%values==0),'empty film source replaces previous field')
  end do
  if (rank==0) print *, 'PASS: MPI AMR film evaporation and gas/wall energy coupling'
  call MPI_Finalize(ierr)
contains
  subroutine advance()
    call advance_mpi_amr_film_coupling(MPI_COMM_WORLD,species,transport,liquids,options,walls,local,local_t, &
      tree,owners,drops,masses,zero,spacing,interval,faces,steps,boundary,ok,gas_change=exported)
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
      call require(all(exported%fields(a)%values==123),'source rollback')
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
      amount=drops(a)%mass
      liquid_total(1)=liquid_total(1)+amount
      liquid_total(5)=liquid_total(5)+amount*e
      do b=1,2
        c=5+liquids(b)%vapor_index
        liquid_total(c)=liquid_total(c)+masses(b,a)
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
