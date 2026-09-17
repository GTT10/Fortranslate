program test_spray_mpi_amr_boundaries
  use mpi_f08
  use spray_mpi_amr_step_mod, only: advance_mpi_amr_mixture_coupling
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
  type(spray_parcel), allocatable :: drops(:),saved_drops(:)
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
  ! All parcels begin on the last rank, including a coarse/fine crossing and wrap.
  allocate(drops(merge(3,0,rank==nranks-1)),masses(2,merge(3,0,rank==nranks-1)))
  do p=1,size(drops)
    drops(p)=spray_parcel(int(p,int64),spacing*[1.99_dp,.8_dp,.8_dp], &
      [20._dp,10._dp,5._dp],1e-12_dp,350._dp,100._dp)
    if (p==2) drops(p)%position=spacing*[.99_dp,.9_dp,.9_dp]
    if (p==3) drops(p)%position=spacing*[3.99_dp,1.1_dp,1.1_dp]
    masses(:,p)=drops(p)%mass*[.25_dp,.75_dp]
  end do
  faces='periodic';faces(1)='slip_wall';faces(2)='outflow'
  allocate(boundary(nv));boundary=123
  do p=1,size(drops)
    drops(p)%velocity=0
    select case(p)
    case(1)
      drops(p)%position=spacing*[3.999_dp,.8_dp,.8_dp];drops(p)%velocity(1)=1000
    case(2)
      drops(p)%position=spacing*[.5_dp,.8_dp,.8_dp]
    case(3)
      drops(p)%position=spacing*[.001_dp,.8_dp,.8_dp];drops(p)%velocity(1)=-1000
    end select
  end do
  saved=local;saved_t=local_t;saved_drops=drops;saved_mass=masses
  last_id=3;options%drag=.false.;options%heat=.false.;options%evaporation=.false.
  interval=2.e-7_dp;options%maximum_substeps=1
  call advance()
  call require(.not.ok.and.steps==0.and.all(boundary==123),'boundary interval exhaustion keeps ledger')
  call unchanged()
  options%maximum_substeps=100000;interval=2.e-8_dp
  if (rank==nranks-1) faces(1)='outflow'
  call advance()
  call require(.not.ok.and.all(boundary==123),'rank-disagreed face type rejects');call unchanged()
  faces(1)='slip_wall';options%maximum_parcels=1
  call advance()
  call require(.not.ok.and.all(boundary==123),'late migration failure keeps boundary ledger');call unchanged()
  options%maximum_parcels=100
  call totals(before)
  call advance()
  call require(ok.and.steps>0,'mixed reflection/outflow step')
  call totals(after)
  scale=max(abs(before),abs(boundary),1.e-30_dp)
  call require(all(abs(after-before-boundary)<=3.e-12_dp*scale),'gas liquid boundary mass momentum energy species budget')
  call require(abs(boundary(1)+1.e-10_dp)<=3.e-24_dp,'exact outgoing represented mass')
  call require(abs(boundary(2)-1.e-7_dp)<=3.e-21_dp,'outflow and wall impulse')
  q=size(drops)
  call MPI_Allreduce(MPI_IN_PLACE,q,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  call require(q==2,'outflow parcel removed')
  do p=1,size(drops)
    call require(drops(p)%id/=1,'outflow ID absent')
    if (drops(p)%id==3) call require(drops(p)%velocity(1)==1000.and.drops(p)%position(1)>0,'specular reflection')
  end do
  if (rank==0) print *, 'PASS: MPI AMR particle physical boundaries and ledger'
  call MPI_Finalize(ierr)
contains
  subroutine advance()
    call advance_mpi_amr_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,options,local,local_t, &
      tree,owners,drops,masses,zero,spacing,interval,steps,ok,last_id,faces,boundary)
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
