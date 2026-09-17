program test_spray_mpi_amr_tree
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,initialize_amr_patch_tree_3d
  use spray_parcel_mod, only: spray_parcel
  use spray_mpi_amr_tree_mod, only: migrate_amr_tree_spray_parcels
  implicit none
  type(amr_patch_3d) :: patches(3)
  type(amr_patch_tree_3d) :: tree,coarse
  type(spray_parcel), allocatable :: parcels(:),saved(:)
  real(dp), allocatable :: masses(:,:),saved_mass(:,:)
  type(spray_parcel) :: expected
  integer :: rank,nranks,ierr,p,id,owners(4),bad_owners(4)
  real(dp) :: zero(3)=0._dp,spacing(3)=1._dp,bad_lower(3)
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
