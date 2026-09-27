program test_spray_mpi_eb_tree
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,initialize_amr_patch_tree_3d
  use spray_parcel_mod, only: spray_parcel
  use spray_mpi_amr_tree_mod, only: migrate_eb_tree_spray_parcels,create_mpi_amr_virtual_particles, &
    create_mpi_amr_ghost_particles,initialize_mpi_amr_sdc_populations,collect_mpi_amr_level_particles,spray_level_population
  implicit none
  type(amr_patch_3d) :: patch(1)
  type(amr_patch_tree_3d) :: tree
  type(spray_parcel), allocatable :: drops(:),saved(:),virtuals(:),ghosts(:)
  real(dp), allocatable :: masses(:,:),saved_mass(:,:),virtual_mass(:,:),ghost_mass(:,:)
  real(dp) :: phi(0:8,0:8,0:8),zero(3)=0,spacing(3)=1
  integer :: rank,nranks,ierr,i,owners(2),p
  integer, allocatable :: ghost_nodes(:)
  type(spray_level_population), allocatable :: physical(:),coarse_copies(:),incoming(:)
  logical :: ok,wrap(3)
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks==2.or.nranks==4,'rank count')
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patch(1),ok)
  call require(ok,'patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1],patch,tree,ok)
  call require(ok,'tree')
  owners=[0,nranks-1]
  do i=0,8
    phi(i,:,:)=real(i,dp)*.5_dp-.75_dp
  end do
  allocate(drops(1),masses(2,1))
  drops(1)=spray_parcel(int(rank+1,int64),[1._dp,1._dp,1._dp],[1._dp,2._dp,3._dp],1e-9_dp,350._dp,100._dp)
  masses(:,1)=drops(1)%mass*[.25_dp,.75_dp]
  saved=drops;saved_mass=masses
  drops(1)%position(1)=2.125_dp
  call create_mpi_amr_ghost_particles(MPI_COMM_WORLD,tree,owners,drops,masses,zero,spacing, &
    [.false.,.false.,.false.],1,1,100,ghosts,ghost_mass,ghost_nodes,ok,phi)
  call require(ok.and.size(ghosts)==merge(nranks,0,rank==nranks-1).and.all(ghost_nodes==2),'EB ghost halo')
  if(rank==0) drops(1)%position(1)=.5_dp
  call create_mpi_amr_ghost_particles(MPI_COMM_WORLD,tree,owners,drops,masses,zero,spacing, &
    [.false.,.false.,.false.],1,1,100,ghosts,ghost_mass,ghost_nodes,ok,phi)
  call require(.not.ok.and.size(ghosts)==merge(nranks,0,rank==nranks-1),'solid ghost rollback')
  do p=1,size(ghosts)
    call require(ghosts(p)%position(1)==2.125_dp.and.ghosts(p)%id==p,'retained ghost payload')
  end do
  drops=saved
  call create_mpi_amr_virtual_particles(MPI_COMM_WORLD,tree,owners,drops,masses,zero,spacing, &
    [.false.,.false.,.false.],0,100,virtuals,virtual_mass,ok,phi)
  call require(ok.and.size(virtuals)==merge(nranks,0,rank==0),'EB virtual coarse ownership')
  call unchanged()
  if(rank==0) drops(1)%position(1)=.5_dp
  call create_mpi_amr_virtual_particles(MPI_COMM_WORLD,tree,owners,drops,masses,zero,spacing, &
    [.false.,.false.,.false.],0,100,virtuals,virtual_mass,ok,phi)
  call require(.not.ok.and.size(virtuals)==merge(nranks,0,rank==0),'solid virtual rejection rollback')
  do p=1,size(virtuals)
    call require(virtuals(p)%id==p.and.virtuals(p)%position(1)==1,'virtual output retained after rejection')
  end do
  drops=saved
  if(rank==0) phi(0,0,0)=phi(0,0,0)+.01_dp
  call migrate(100)
  call require(.not.ok,'geometry disagreement')
  call unchanged()
  phi(0,0,0)=-.75_dp
  if(rank==0) drops(1)%position(1)=.5_dp
  saved=drops
  call migrate(100)
  call require(.not.ok,'solid point after private migration')
  call unchanged()
  drops=saved;drops(1)%position(1)=1;saved=drops
  call migrate(0)
  call require(.not.ok,'capacity rollback');call unchanged()
  call migrate_eb_tree_spray_parcels(MPI_COMM_WORLD,tree,owners,phi,drops,zero,spacing,100,ok,masses,0)
  call require(ok.and.size(drops)==merge(nranks,0,rank==0),'covered coarse EB owner')
  do p=1,size(drops)
    call require(drops(p)%id==p.and.all(masses(:,p)==1e-9_dp*[.25_dp,.75_dp]),'coarse EB payload')
  end do
  call migrate(100)
  call require(ok.and.size(drops)==merge(nranks,0,rank==nranks-1),'finest fluid owner')
  do p=1,size(drops)
    call require(drops(p)%id==p.and.all(drops(p)%velocity==[1._dp,2._dp,3._dp]),'ID and velocity payload')
    call require(all(masses(:,p)==1e-9_dp*[.25_dp,.75_dp]),'component payload')
    drops(p)%position(1)=3
  end do
  call migrate(100)
  call require(ok.and.size(drops)==merge(nranks,0,rank==0),'coarse fluid owner')
  owners=[nranks-1,0]
  call migrate(100)
  call require(ok.and.size(drops)==merge(nranks,0,rank==nranks-1),'changed gas ownership')
  ! Periodic tangential axes wrap before nodal fluid checks and level routing.
  wrap=[.false.,.true.,.true.]
  do p=1,size(drops)
    drops(p)%position=[1._dp,-3._dp,9._dp]
  end do
  call migrate_eb_tree_spray_parcels(MPI_COMM_WORLD,tree,owners,phi,drops,zero,spacing,100,ok,masses,periodic=wrap)
  call require(ok.and.size(drops)==merge(nranks,0,rank==owners(2)),'periodic EB migration selects fine owner')
  do p=1,size(drops)
    call require(all(drops(p)%position==1._dp).and.drops(p)%id==p,'periodic EB wrap and ID order')
    call require(all(masses(:,p)==1e-9_dp*[.25_dp,.75_dp]),'periodic EB component inventory')
  end do
  call initialize_mpi_amr_sdc_populations(MPI_COMM_WORLD,tree,owners,drops,masses,zero,spacing,wrap,1,100, &
    physical,coarse_copies,incoming,ok,phi)
  call require(ok.and.size(physical(2)%parcels)==size(drops),'periodic EB physical hierarchy')
  call require(size(coarse_copies(1)%parcels)==merge(nranks,0,rank==owners(1)), &
    'periodic EB coarse virtual ownership')
  call collect_mpi_amr_level_particles(MPI_COMM_WORLD,tree,owners,physical,zero,spacing,wrap,100,drops,masses,ok,phi)
  call require(ok.and.size(drops)==merge(nranks,0,rank==owners(2)),'periodic EB physical collection')
  do p=1,size(drops)
    drops(p)%position=[2.125_dp,-3._dp,9._dp]
  end do
  call create_mpi_amr_ghost_particles(MPI_COMM_WORLD,tree,owners,drops,masses,zero,spacing,wrap,1,1,100, &
    ghosts,ghost_mass,ghost_nodes,ok,phi)
  call require(ok.and.size(ghosts)==merge(nranks,0,rank==owners(2)),'periodic EB grown ghost routing')
  do p=1,size(ghosts)
    call require(all(ghosts(p)%position==[2.125_dp,1._dp,1._dp]).and.ghosts(p)%id==p,'wrapped ghost payload')
  end do
  ! Geometry must be periodic on every requested axis, even for empty ranks.
  call migrate_eb_tree_spray_parcels(MPI_COMM_WORLD,tree,owners,phi,drops,zero,spacing,100,ok,masses, &
    periodic=[.true.,.false.,.false.])
  call require(.not.ok,'nonperiodic nodal seam rejects')
  do p=1,size(drops)
    call require(all(drops(p)%position==[2.125_dp,-3._dp,9._dp]),'seam rejection retains unwrapped coordinates')
  end do
  if(rank==0) wrap(3)=.false.
  call migrate_eb_tree_spray_parcels(MPI_COMM_WORLD,tree,owners,phi,drops,zero,spacing,100,ok,masses,periodic=wrap)
  call require(.not.ok,'rank-disagreed periodic controls reject')
  ! A coarse particle near the right seam supplies a left fine ghost box.
  do i=0,8
    phi(:,i,:)=real(i,dp)*.5_dp-.75_dp
  end do
  wrap=[.true.,.false.,.true.]
  do p=1,size(drops)
    drops(p)%position=[3.875_dp,1._dp,1._dp]
  end do
  call create_mpi_amr_ghost_particles(MPI_COMM_WORLD,tree,owners,drops,masses,zero,spacing,wrap,1,1,100, &
    ghosts,ghost_mass,ghost_nodes,ok,phi)
  call require(ok.and.size(ghosts)==merge(nranks,0,rank==owners(2)).and.all(ghost_nodes==2), &
    'EB ghost ownership reaches fine box across periodic seam')
  do p=1,size(ghosts)
    call require(all(ghosts(p)%position==[3.875_dp,1._dp,1._dp]),'periodic ghost keeps canonical position')
  end do
  deallocate(drops,masses);allocate(drops(0),masses(2,0))
  call migrate(100)
  call require(ok,'empty population')
  if(rank==0) print *, 'PASS: transactional MPI EB tree particle ownership'
  call MPI_Finalize(ierr)
contains
  subroutine migrate(capacity)
    integer, intent(in) :: capacity
    call migrate_eb_tree_spray_parcels(MPI_COMM_WORLD,tree,owners,phi,drops,zero,spacing,capacity,ok,masses)
  end subroutine
  subroutine unchanged()
    call require(size(drops)==size(saved),'rollback count')
    call require(all(drops%id==saved%id).and.all(masses==saved_mass),'rollback identity and components')
    call require(all(drops(1)%position==saved(1)%position),'rollback position')
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
