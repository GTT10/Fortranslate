program test_spray_mpi_amr_balance
  use mpi_f08
  use spray_mpi_amr_balance_mod, only: rebalance_mpi_amr_spray
  use spray_mpi_amr_remap_mod, only: remap_mpi_amr_spray
  use amr_patch_tree_remap_3d_mod, only: remap_amr_patch_tree_reactive_3d
  use spray_amr_tree_partition_mod, only: locate_amr_tree_spray_cells
  use spray_wall_impact_mod, only: spray_wall_film
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
  type(amr_patch_tree_3d) :: tree,target,original_tree,refined_tree,root_tree
  type(amr_tree_state_3d) :: gas,temp,saved,saved_t,expected,expected_t,local,local_t
  type(gas_transport_species), allocatable :: transport(:)
  type(spray_liquid) :: liquids(2)
  type(spray_options) :: options
  type(spray_wall_film), allocatable :: films(:),saved_films(:)
  real(dp), allocatable :: film_mass(:,:),saved_film_mass(:,:)
  type(spray_parcel), allocatable :: drops(:),saved_drops(:)
  real(dp), allocatable :: masses(:,:),saved_mass(:,:)
  real(dp), allocatable :: root(:,:,:,:),primitive(:),cell(:),y(:),before(:),after(:),scale(:)
  real(dp), allocatable :: positions(:,:)
  integer, allocatable :: locations(:,:)
  real(dp) :: spacing(3)=.001_dp,zero(3)=0._dp,t,sound,r,interval
  integer :: nv,node,n(3),i,j,k,p,water,nitrogen,steps
  integer :: rank,nranks,ierr,q,mode,capacity
  integer, allocatable :: owners(:),target_owners(:),saved_owners(:)
  type(spray_parcel) :: all_drops(3)
  real(dp) :: all_mass(2,3)
  logical :: ok,changed
  real(dp) :: costs(3)
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
  original_tree=tree
  call initialize_amr_patch_tree_3d([4,4,4],[integer::],patches(:0),root_tree,ok)
  call require(ok,'root-only tree')
  call initialize_amr_patch_3d(4,4,4,3,4,1,4,1,4,2,patches(1),ok)
  call require(ok,'moved child')
  call initialize_amr_patch_3d(4,8,8,3,4,1,2,1,2,2,patches(2),ok)
  call require(ok,'moved nested child')
  call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,refined_tree,ok)
  call require(ok,'moved tree')
  ! Nonuniform positive source gas makes coarsening a real volume average.
  do node=1,3
    gas%fields(node)%values=gas%fields(node)%values*(1._dp+.1_dp*node)
  end do
  call composite_integrals_amr_patch_tree_3d(tree,gas,spacing,before,ok)
  call require(ok,'initial gas inventory')
  owners=[0,1,nranks-1]
  allocate(local%fields(3),local_t%fields(3))
  do node=1,3
    if (owners(node)/=rank) cycle
    local%fields(node)%values=gas%fields(node)%values
    local_t%fields(node)%values=temp%fields(node)%values
  end do
  do p=1,3
    all_drops(p)=spray_parcel(int(p,int64),spacing*[real(p,dp)+.2_dp,.8_dp,.8_dp], &
      [20._dp,10._dp,5._dp],p*1e-12_dp,350._dp,100._dp)
    all_drops(p)%parent_id=p-1;all_drops(p)%material_id=p
    all_drops(p)%deformation=.01_dp*p;all_drops(p)%deformation_rate=.02_dp*p
    all_drops(p)%rt_time=.03_dp*p;all_drops(p)%khrt_initial_count=100
    all_drops(p)%khrt_shed_volume=1.e-20_dp*p
    all_mass(:,p)=all_drops(p)%mass*[.25_dp,.75_dp]
    call require(valid_spray_parcel(all_drops(p)),'source particle valid')
  end do
  allocate(drops(0),masses(2,0))
  if (rank==nranks-1) then
    drops=all_drops;masses=all_mass
  end if
  saved=local;saved_t=local_t;saved_drops=drops;saved_mass=masses;saved_owners=owners
  target=root_tree;target_owners=[0];capacity=1
  call remap()
  call require(.not.ok,'particle capacity rejects coarsening');call unchanged()
  capacity=100
  if (rank==owners(3)) local%fields(3)%values(1,1,1,1)=-1
  saved=local
  call remap()
  call require(.not.ok,'invalid source gas rejects before coarsening');call unchanged()
  if (rank==owners(3)) local%fields(3)%values=gas%fields(3)%values
  saved=local
  ! Valid target with a different root resolution: particle migration can succeed,
  ! but gas remapping must fail and preserve the original topology and payloads.
  call initialize_amr_patch_tree_3d([3,4,4],[integer::],patches(:0),target,ok)
  call require(ok,'incompatible root')
  call remap()
  call require(.not.ok,'late gas failure after private particle migration');call unchanged()
  do mode=1,3
    if (mode==1) target=root_tree
    if (mode==2) target=refined_tree
    if (mode==3) target=original_tree
    target_owners=[(modulo(node,nranks),node=1,size(target%nodes))]
    call remap_amr_patch_tree_reactive_3d(tree,gas,temp,target,species,expected,expected_t,ok)
    call require(ok,'serial gas reference')
    call remap()
    call require(ok,'coupled topology change')
    call require(size(tree%nodes)==size(target%nodes).and.all(owners==target_owners),'published topology and owners')
    do node=1,size(tree%nodes)
      if (owners(node)==rank) then
        call require(all(local%fields(node)%values==expected%fields(node)%values),'serial conserved remap')
        call require(all(local_t%fields(node)%values==expected_t%fields(node)%values),'serial temperature recovery')
      else
        call require(.not.allocated(local%fields(node)%values),'nonowner storage absent')
      end if
    end do
    if (allocated(positions)) deallocate(positions)
    allocate(positions(3,size(drops)))
    do p=1,size(drops)
      q=int(drops(p)%id);positions(:,p)=drops(p)%position
      call require(all(drops(p)%position==all_drops(q)%position).and. &
        all(drops(p)%velocity==all_drops(q)%velocity),'particle phase-space preserved')
      call require(drops(p)%mass==all_drops(q)%mass.and.drops(p)%temperature==all_drops(q)%temperature.and. &
        drops(p)%multiplicity==all_drops(q)%multiplicity.and.all(masses(:,p)==all_mass(:,q)),'liquid inventory preserved')
      call require(drops(p)%parent_id==all_drops(q)%parent_id.and.drops(p)%material_id==all_drops(q)%material_id.and. &
        drops(p)%deformation==all_drops(q)%deformation.and.drops(p)%deformation_rate==all_drops(q)%deformation_rate.and. &
        drops(p)%rt_time==all_drops(q)%rt_time.and.drops(p)%khrt_initial_count==all_drops(q)%khrt_initial_count.and. &
        drops(p)%khrt_shed_volume==all_drops(q)%khrt_shed_volume,'lineage and breakup state preserved')
    end do
    call locate_amr_tree_spray_cells(tree,positions,zero,spacing,[.true.,.true.,.true.],locations,ok)
    call require(ok,'new particle locations')
    do p=1,size(drops)
      call require(owners(locations(1,p))==rank,'finest new owner')
    end do
    q=size(drops)
    call MPI_Allreduce(MPI_IN_PLACE,q,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(q==3,'global parcel count')
    call composite_integrals_amr_patch_tree_3d(tree,expected,spacing,after,ok)
    call require(ok.and.all(abs(after-before)<=3.e-13_dp*max(abs(before),1.e-30_dp)),'gas conservation across regrids')
    gas=expected;temp=expected_t
  end do
  target=tree;target_owners=0
  call remap()
  call require(ok,'concentrate gas and particles on rank zero')
  saved=local;saved_t=local_t;saved_drops=drops;saved_mass=masses;saved_owners=owners
  allocate(films(merge(1,0,rank==0)),film_mass(2,merge(1,0,rank==0)))
  if (rank==0) then
    films(1)=spray_wall_film([0._dp,.00025_dp,.00025_dp],[1._dp,0._dp,0._dp], &
      1e-9_dp,330._dp,.001_dp,1e-6_dp,2)
    film_mass(:,1)=films(1)%mass*[.25_dp,.75_dp]
  end if
  saved_films=films;saved_film_mass=film_mass
  costs=[100._dp,80._dp,60._dp]
  call rebalance_mpi_amr_spray(MPI_COMM_WORLD,tree,owners,nv,local,local_t,drops,masses, &
    zero,spacing,[.true.,.true.,.true.],costs,0._dp,0,changed,ok,films,film_mass,100)
  call require(.not.ok.and..not.changed,'late particle capacity failure rolls back balancing');call unchanged()
  call rebalance_mpi_amr_spray(MPI_COMM_WORLD,tree,owners,nv,local,local_t,drops,masses, &
    zero,spacing,[.true.,.true.,.true.],costs,0._dp,100,changed,ok,films,film_mass,0)
  call require(.not.ok.and..not.changed,'late film capacity failure rolls back balancing');call unchanged()
  call rebalance_mpi_amr_spray(MPI_COMM_WORLD,tree,owners,nv,local,local_t,drops,masses, &
    zero,spacing,[.true.,.true.,.true.],costs,0._dp,100,changed,ok,films,film_mass,100)
  call require(ok.and.changed.and.any(owners/=saved_owners),'forced cost imbalance redistributes')
  do node=1,3
    if (owners(node)==rank) then
      call require(all(local%fields(node)%values==gas%fields(node)%values),'gas redistribution exact')
      call require(all(local_t%fields(node)%values==temp%fields(node)%values),'temperature redistribution exact')
    else
      call require(.not.allocated(local%fields(node)%values),'nonowner gas absent')
    end if
  end do
  deallocate(positions);allocate(positions(3,size(drops)))
  do p=1,size(drops)
    q=int(drops(p)%id);positions(:,p)=drops(p)%position
    call require(all(masses(:,p)==all_mass(:,q)),'balanced particle composition')
  end do
  call locate_amr_tree_spray_cells(tree,positions,zero,spacing,[.true.,.true.,.true.],locations,ok)
  call require(ok,'balanced locations')
  do p=1,size(drops)
    call require(owners(locations(1,p))==rank,'balanced finest owner')
  end do
  call require(size(films)==merge(1,0,rank==owners(3)),'film follows nested wall owner')
  do p=1,size(films)
    call require(all(films(p)%position==[0._dp,.00025_dp,.00025_dp]).and. &
      all(films(p)%normal==[1._dp,0._dp,0._dp]),'exact film wall geometry')
    call require(films(p)%mass==1e-9_dp.and.films(p)%temperature==330.and. &
      films(p)%diameter==.001_dp.and.films(p)%height==1e-6_dp.and.films(p)%material_id==2,'film payload')
    call require(all(film_mass(:,p)==1e-9_dp*[.25_dp,.75_dp]),'film components')
  end do
  saved_films=films;saved_film_mass=film_mass
  saved=local;saved_t=local_t;saved_drops=drops;saved_mass=masses;saved_owners=owners;costs=0
  call rebalance_mpi_amr_spray(MPI_COMM_WORLD,tree,owners,nv,local,local_t,drops,masses, &
    zero,spacing,[.true.,.true.,.true.],costs,0._dp,100,changed,ok,films,film_mass,100)
  call require(ok.and..not.changed,'zero measured work leaves ownership unchanged');call unchanged()
  if (rank==0) print *, 'PASS: transactional MPI AMR gas and particle balancing'
  call MPI_Finalize(ierr)
contains
  subroutine remap()
    call remap_mpi_amr_spray(MPI_COMM_WORLD,species,tree,owners,local,local_t,drops,masses,target,target_owners, &
      zero,spacing,[.true.,.true.,.true.],capacity,ok)
  end subroutine
  subroutine unchanged()
    integer :: a
    if (allocated(films)) then
      call require(size(films)==size(saved_films).and.all(film_mass==saved_film_mass),'film inventory rollback')
      do a=1,size(films)
        call require(all(films(a)%position==saved_films(a)%position).and. &
          films(a)%mass==saved_films(a)%mass.and.films(a)%temperature==saved_films(a)%temperature,'film rollback')
      end do
    end if
    call require(size(tree%nodes)==size(original_tree%nodes).and.all(owners==saved_owners),'topology rollback')
    call require(size(drops)==size(saved_drops).and.all(masses==saved_mass),'particle count and components rollback')
    do a=1,size(drops)
      call require(drops(a)%id==saved_drops(a)%id.and.all(drops(a)%position==saved_drops(a)%position),'particle rollback')
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
