program test_spray_mpi_amr_tree
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,initialize_amr_patch_tree_3d
  use spray_parcel_mod, only: spray_parcel
  use spray_mpi_amr_tree_mod, only: migrate_amr_tree_spray_parcels,create_mpi_amr_virtual_particles, &
    create_mpi_amr_ghost_particles
  implicit none
  type(amr_patch_3d) :: patches(3)
  type(amr_patch_tree_3d) :: tree,coarse
  type(spray_parcel), allocatable :: parcels(:),saved(:),virtuals(:),saved_virtuals(:),ghosts(:)
  real(dp), allocatable :: masses(:,:),saved_mass(:,:),virtual_mass(:,:),saved_virtual_mass(:,:),ghost_mass(:,:)
  type(spray_parcel) :: expected
  integer :: rank,nranks,ierr,p,id,owners(4),bad_owners(4),ceiling
  real(dp) :: zero(3)=0._dp,spacing(3)=1._dp,bad_lower(3)
  integer, allocatable :: ghost_nodes(:)
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks==2.or.nranks==4,'rank fixture')
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patches(1),ok)
  call require(ok,'left patch')
  call initialize_amr_patch_3d(4,4,4,3,4,1,4,1,4,2,patches(2),ok)
  call require(ok,'right patch')
  call initialize_amr_patch_3d(4,8,8,1,2,1,4,1,4,2,patches(3),ok)
  call require(ok,'nested patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1,1,2],patches,tree,ok)
  call require(ok,'tree')
  owners=[0,0,1,nranks-1]
  allocate(parcels(2),masses(2,2))
  do p=1,2
    parcels(p)=make_parcel(2*rank+p)
    masses(:,p)=parcels(p)%mass*[.25_dp,.75_dp]
  end do
  saved=parcels;saved_mass=masses
  do p=1,size(parcels)
    parcels(p)%position=[merge(1.125_dp,3.875_dp,mod(parcels(p)%id,2_int64)==1),.25_dp,.25_dp]
  end do
  call create_mpi_amr_ghost_particles(MPI_COMM_WORLD,tree,owners,parcels,masses,zero,spacing, &
    [.true.,.false.,.false.],2,1,100,ghosts,ghost_mass,ghost_nodes,ok)
  call require(ok.and.size(ghosts)==merge(2*nranks,0,rank==nranks-1),'fine and periodic ghost halo')
  call require(all(ghost_nodes==4),'ghost target node')
  do p=1,size(ghosts)
    expected=make_parcel(int(ghosts(p)%id))
    call require(all(ghost_mass(:,p)==expected%mass*[.25_dp,.75_dp]),'ghost components')
    call require(ghosts(p)%position(1)==merge(1.125_dp,3.875_dp,mod(ghosts(p)%id,2_int64)==1), &
      'ghost position unchanged across periodic image')
  end do
  call create_mpi_amr_ghost_particles(MPI_COMM_WORLD,tree,owners,parcels,masses,zero,spacing, &
    [.true.,.false.,.false.],2,1,0,ghosts,ghost_mass,ghost_nodes,ok)
  call require(.not.ok.and.size(ghosts)==merge(2*nranks,0,rank==nranks-1).and.all(ghost_nodes==4), &
    'ghost capacity rollback')
  call create_mpi_amr_ghost_particles(MPI_COMM_WORLD,tree,owners,parcels,masses,zero,spacing, &
    [.true.,.false.,.false.],2,merge(0,1,rank==0),100,ghosts,ghost_mass,ghost_nodes,ok)
  call require(.not.ok.and.size(ghosts)==merge(2*nranks,0,rank==nranks-1),'ghost width agreement')
  call create_mpi_amr_ghost_particles(MPI_COMM_WORLD,tree,owners,parcels,masses,zero,spacing, &
    [.true.,.false.,.false.],2,0,100,ghosts,ghost_mass,ghost_nodes,ok)
  call require(ok.and.size(ghosts)==0.and.size(ghost_nodes)==0,'zero halo excludes coarse particles')
  call create_mpi_amr_ghost_particles(MPI_COMM_WORLD,tree,owners,parcels,masses,zero,spacing, &
    [.true.,.false.,.false.],1,8,100,ghosts,ghost_mass,ghost_nodes,ok)
  call require(ok.and.size(ghosts)==0,'ghosts exclude particles already at target level')
  parcels=saved
  call create_mpi_amr_virtual_particles(MPI_COMM_WORLD,tree,owners,parcels,masses,zero,spacing, &
    [.true.,.false.,.false.],1,100,virtuals,virtual_mass,ok)
  call require(ok.and.size(virtuals)==merge(nranks,0,rank==0),'virtuals select only deeper levels')
  do p=1,size(virtuals)
    call require(virtuals(p)%id==2*p-1,'virtual source identity')
    expected=make_parcel(int(virtuals(p)%id))
    call require(all(virtuals(p)%position==expected%position).and. &
      all(virtual_mass(:,p)==expected%mass*[.25_dp,.75_dp]),'virtual unaggregated payload')
  end do
  saved_virtuals=virtuals;saved_virtual_mass=virtual_mass
  call create_mpi_amr_virtual_particles(MPI_COMM_WORLD,tree,owners,parcels,masses,zero,spacing, &
    [.true.,.false.,.false.],0,0,virtuals,virtual_mass,ok)
  call require(.not.ok.and.all(virtuals%id==saved_virtuals%id).and. &
    all(virtual_mass==saved_virtual_mass),'virtual capacity rollback')
  call create_mpi_amr_virtual_particles(MPI_COMM_WORLD,tree,owners,parcels,masses,zero,spacing, &
    [.true.,.false.,.false.],merge(0,1,rank==0),100,virtuals,virtual_mass,ok)
  call require(.not.ok.and.all(virtuals%id==saved_virtuals%id).and. &
    all(virtual_mass==saved_virtual_mass),'virtual level disagreement rollback')
  call create_mpi_amr_virtual_particles(MPI_COMM_WORLD,tree,owners,parcels,masses,zero,spacing, &
    [.true.,.false.,.false.],0,100,virtuals,virtual_mass,ok)
  call require(ok.and.size(virtuals)==merge(2*nranks,0,rank==0),'all descendants represented at root')
  call create_mpi_amr_virtual_particles(MPI_COMM_WORLD,tree,owners,parcels,masses,zero,spacing, &
    [.true.,.false.,.false.],2,100,virtuals,virtual_mass,ok)
  call require(ok.and.size(virtuals)==0.and.size(virtual_mass,2)==0,'finest virtual population empty')
  call require(all(parcels%id==saved%id).and.all(masses==saved_mass),'physical population unchanged')
  ceiling=merge(0,1,rank==0)
  call migrate_amr_tree_spray_parcels(MPI_COMM_WORLD,tree,owners,parcels,zero,spacing, &
    [.true.,.false.,.false.],100,ok,masses,ceiling)
  call require(.not.ok.and.all(parcels%id==saved%id).and.all(masses==saved_mass),'level disagreement rollback')
  call migrate_amr_tree_spray_parcels(MPI_COMM_WORLD,tree,owners,parcels,zero,spacing, &
    [.true.,.false.,.false.],100,ok,masses,3)
  call require(.not.ok.and.all(parcels%id==saved%id).and.all(masses==saved_mass),'invalid level rollback')
  do ceiling=0,1
    call migrate_amr_tree_spray_parcels(MPI_COMM_WORLD,tree,owners,parcels,zero,spacing, &
      [.true.,.false.,.false.],100,ok,masses,ceiling)
    call require(ok,'bounded level migration')
    if(ceiling==0) then
      call require(size(parcels)==merge(2*nranks,0,rank==0),'covered root owns all copies')
    else
      call require(size(parcels)==merge(nranks,0,rank<=1),'level-one owner distribution')
    end if
    do p=1,size(parcels)
      expected=make_parcel(int(parcels(p)%id))
      call require(all(parcels(p)%position==expected%position),'level view retains position')
      call require(all(masses(:,p)==expected%mass*[.25_dp,.75_dp]),'level view retains components')
      if(ceiling==1) call require(mod(parcels(p)%id,2_int64)==int(1-rank,int64),'level-one ID ownership')
    end do
  end do
  parcels=saved;masses=saved_mass
  call migrate(nranks-1)
  call require(.not.ok.and.all(parcels%id==saved%id).and.all(masses==saved_mass),'capacity rollback')
  bad_owners=owners
  if (rank==0) bad_owners(1)=1
  call migrate_amr_tree_spray_parcels(MPI_COMM_WORLD,tree,bad_owners,parcels,zero,spacing, &
    [.true.,.false.,.false.],100,ok,masses)
  call require(.not.ok.and.all(parcels%id==saved%id).and.all(masses==saved_mass),'owner-map disagreement')
  bad_lower=zero
  if (rank==0) bad_lower(1)=1._dp
  call migrate_amr_tree_spray_parcels(MPI_COMM_WORLD,tree,owners,parcels,bad_lower,spacing, &
    [.true.,.false.,.false.],100,ok,masses)
  call require(.not.ok.and.all(parcels%id==saved%id),'physical domain disagreement')
  call migrate(100)
  call require(ok,'finest owner migration')
  if (nranks==2) then
    call require(size(parcels)==merge(4,0,rank==1),'two-rank ownership')
  else
    call require(size(parcels)==merge(4,0,rank==1.or.rank==3),'four-rank empty destinations')
  end if
  do p=1,size(parcels)
    id=int(parcels(p)%id);expected=make_parcel(id)
    if (p>1) call require(parcels(p)%id>parcels(p-1)%id,'receive ID order')
    call require(all(masses(:,p)==expected%mass*[.25_dp,.75_dp]),'component payload')
    call require(all(parcels(p)%position==expected%position).and. &
      all(parcels(p)%velocity==expected%velocity),'kinematic payload')
    call require(parcels(p)%parent_id==expected%parent_id.and. &
      parcels(p)%material_id==expected%material_id,'identity payload')
    call require(all([parcels(p)%mass,parcels(p)%temperature,parcels(p)%multiplicity, &
      parcels(p)%deformation,parcels(p)%deformation_rate,parcels(p)%rt_time, &
      parcels(p)%khrt_initial_count,parcels(p)%khrt_shed_volume]== &
      [expected%mass,expected%temperature,expected%multiplicity,expected%deformation, &
      expected%deformation_rate,expected%rt_time,expected%khrt_initial_count,expected%khrt_shed_volume]), &
      'complete physical and breakup payload')
    parcels(p)%position=[1.25_dp,.25_dp,.25_dp]
  end do
  call migrate(100)
  call require(ok.and.size(parcels)==merge(2*nranks,0,rank==0),'particle crossing')
  owners=modulo(owners+1,nranks)
  call migrate(100)
  call require(ok.and.size(parcels)==merge(2*nranks,0,rank==1),'gas ownership reassignment')
  call initialize_amr_patch_tree_3d([4,4,4],[integer::],patches(:0),coarse,ok)
  call require(ok,'coarsened layout')
  call migrate_amr_tree_spray_parcels(MPI_COMM_WORLD,coarse,[0],parcels,zero,spacing, &
    [.true.,.false.,.false.],100,ok,masses)
  call require(ok.and.size(parcels)==merge(2*nranks,0,rank==0),'regrid migration')
  call migrate_amr_tree_spray_parcels(MPI_COMM_WORLD,coarse,[nranks-1],parcels,zero,spacing, &
    [.true.,.false.,.false.],100,ok)
  call require(ok.and.size(parcels)==merge(2*nranks,0,rank==nranks-1),'single-component payload')
  deallocate(parcels,masses);allocate(parcels(2),masses(2,2))
  do p=1,2
    parcels(p)=make_parcel(p);masses(:,p)=parcels(p)%mass*[.25_dp,.75_dp]
  end do
  saved=parcels;saved_mass=masses
  call migrate_amr_tree_spray_parcels(MPI_COMM_WORLD,coarse,[0],parcels,zero,spacing, &
    [.true.,.false.,.false.],100,ok,masses)
  call require(.not.ok.and.all(parcels%id==saved%id).and.all(masses==saved_mass),'duplicate receive rollback')
  deallocate(parcels,masses);allocate(parcels(0),masses(2,0))
  call migrate(100)
  call require(ok.and.size(parcels)==0.and.size(masses,2)==0,'globally empty population')
  if (rank==0) print *, 'PASS: AMR spray MPI migration'
  call MPI_Finalize(ierr)
contains
  function make_parcel(number) result(parcel)
    integer, intent(in) :: number
    type(spray_parcel) :: parcel
    parcel=spray_parcel(int(number,int64),[merge(4.25_dp,3.25_dp,mod(number,2)==1),.25_dp,.25_dp], &
      [1._dp,2._dp,3._dp],number*1.e-9_dp,350._dp,100._dp)
    parcel%parent_id=number-1;parcel%material_id=1+mod(number,2)
    parcel%deformation=.1_dp*number;parcel%deformation_rate=.2_dp*number
    parcel%rt_time=.01_dp*number;parcel%khrt_initial_count=100._dp;parcel%khrt_shed_volume=number*1.e-15_dp
  end function
  subroutine migrate(capacity)
    integer, intent(in) :: capacity
    call migrate_amr_tree_spray_parcels(MPI_COMM_WORLD,tree,owners,parcels,zero,spacing, &
      [.true.,.false.,.false.],capacity,ok,masses)
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
