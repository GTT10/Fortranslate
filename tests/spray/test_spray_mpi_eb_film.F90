program test_spray_mpi_eb_film
  use mpi_f08
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,initialize_amr_patch_tree_3d
  use spray_wall_impact_mod, only: spray_wall_film
  use spray_mpi_amr_film_mod, only: migrate_mpi_amr_wall_films,initialize_mpi_eb_film_levels,spray_level_film_population
  implicit none
  type(amr_patch_3d) :: patches(3)
  type(amr_patch_tree_3d) :: tree,coarse
  type(spray_wall_film), allocatable :: films(:),saved(:)
  type(spray_wall_film) :: expected
  real(dp), allocatable :: masses(:,:),saved_mass(:,:)
  real(dp) :: zero(3)=0._dp,spacing(3)=1._dp,total(2),initial(2)
  real(dp) :: phi(0:16,0:16,0:16)
  type(spray_level_film_population), allocatable :: physical(:),kept(:)
  integer, allocatable :: nodes(:)
  real(dp), allocatable :: probes(:,:)
  real(dp) :: probe_offset
  integer :: level,c
  integer :: i
  integer :: rank,nranks,ierr,p,id,owners(4),bad_owners(4),expected_count
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
  do i=0,16
    phi(i,:,:)=abs(.25_dp*i-2)-1
  end do
  allocate(films(3),masses(2,3))
  do p=1,3
    films(p)=make_film(3*rank+p)
    masses(:,p)=films(p)%mass*[.25_dp,.75_dp]
  end do
  saved=films;saved_mass=masses
  initial=sum(masses,dim=2)
  call MPI_Allreduce(MPI_IN_PLACE,initial,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call migrate(0)
  call require(.not.ok,'receive capacity rejects');call unchanged()
  bad_owners=owners
  if (rank==0) bad_owners(1)=1
  call migrate_mpi_amr_wall_films(MPI_COMM_WORLD,tree,bad_owners,films,masses,zero,spacing,100,ok,phi)
  call require(.not.ok,'rank-disagreed owner map rejects');call unchanged()
  if (rank==0) films(1)%position(1)=-.01_dp
  call migrate(100)
  call require(.not.ok,'exterior film rejects collectively')
  if (rank==0) films(1)=saved(1)
  call unchanged()
  if(rank==0) films(1)%normal=-films(1)%normal
  call migrate(100)
  call require(.not.ok,'reversed EB normal rejects collectively')
  if(rank==0) films(1)=saved(1)
  call unchanged()
  if(rank==0) phi(0,0,0)=2
  call migrate(100)
  call require(.not.ok,'rank-disagreed nodal geometry rejects')
  if(rank==0) phi(0,0,0)=1
  call unchanged()
  if(rank==0) then
    call migrate_mpi_amr_wall_films(MPI_COMM_WORLD,tree,owners,films,masses,zero,spacing,100,ok)
  else
    call migrate(100)
  end if
  call require(.not.ok,'rank-disagreed EB presence rejects');call unchanged()
  call migrate(100)
  call require(ok,'finest fluid-side owner migration')
  expected_count=0
  if (rank==0) expected_count=expected_count+nranks
  if (rank==1) expected_count=expected_count+nranks
  if (rank==nranks-1) expected_count=expected_count+nranks
  call require(size(films)==expected_count,'closed lower and upper wall ownership')
  do p=1,size(films)
    id=nint(films(p)%temperature)-300
    select case(mod(id-1,3))
    case(0);call require(rank==nranks-1,'nested lower wall')
    case(1);call require(rank==1,'right upper wall')
    case(2);call require(rank==0,'left upper transverse wall')
    end select
    if (p>1) call require(films(p)%temperature>films(p-1)%temperature,'stable source rank/local ordering')
  end do
  call payload()
  allocate(nodes(1),probes(3,1));nodes=17;probes=19
  call migrate_mpi_amr_wall_films(MPI_COMM_WORLD,tree,owners,films,masses,zero,spacing,0,ok,phi, &
    nodes=nodes,probes=probes)
  call require(.not.ok.and.all(nodes==17).and.all(probes==19),'migration metadata rolls back on capacity failure')
  if(rank==0) then
    call migrate_mpi_amr_wall_films(MPI_COMM_WORLD,tree,owners,films,masses,zero,spacing,100,ok,phi)
  else
    call migrate_mpi_amr_wall_films(MPI_COMM_WORLD,tree,owners,films,masses,zero,spacing,100,ok,phi, &
      nodes=nodes,probes=probes)
  end if
  call require(.not.ok.and.all(nodes==17).and.all(probes==19),'metadata presence agreed across ranks')
  call initialize_mpi_eb_film_levels(MPI_COMM_WORLD,tree,owners,films,masses,zero,spacing,100,phi,physical,ok)
  call require(ok.and.size(physical)==3,'physical film hierarchy initialized')
  probe_offset=max(1e-8_dp*.25_dp,128*epsilon(1._dp)*4)
  total=0;c=0
  do level=1,3
    c=c+size(physical(level)%films)
    call require(size(physical(level)%component_mass,2)==size(physical(level)%films).and. &
      size(physical(level)%probes,2)==size(physical(level)%films).and. &
      size(physical(level)%nodes)==size(physical(level)%films),'film level payload dimensions')
    do p=1,size(physical(level)%films)
      id=nint(physical(level)%films(p)%temperature)-300;expected=make_film(id)
      call require(tree%nodes(physical(level)%nodes(p))%level==level-1.and. &
        owners(physical(level)%nodes(p))==rank,'physical film level and owner agree')
      call require(all(physical(level)%films(p)%position==expected%position).and. &
        all(physical(level)%films(p)%normal==expected%normal),'level initialization preserves wall geometry')
      call require(all(physical(level)%probes(:,p)==expected%position+probe_offset*expected%normal), &
        'film hierarchy exports exact fluid-side ownership probe')
      call require(physical(level)%films(p)%material_id==expected%material_id.and. &
        all(physical(level)%component_mass(:,p)==expected%mass*[.25_dp,.75_dp]),'film hierarchy material and components')
    end do
    total=total+sum(physical(level)%component_mass,dim=2)
  end do
  call require(c==size(films),'every local film assigned to exactly one level')
  call require(size(physical(1)%films)==0,'covered root has allocated empty physical film population')
  call MPI_Allreduce(MPI_IN_PLACE,total,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call require(all(abs(total-initial)<=2e-12_dp*initial),'level film inventory')
  kept=physical
  call initialize_mpi_eb_film_levels(MPI_COMM_WORLD,tree,owners,films,masses,zero,spacing,0,phi,physical,ok)
  call require(.not.ok,'film hierarchy capacity failure')
  do level=1,3
    call require(all(physical(level)%component_mass==kept(level)%component_mass).and. &
      all(physical(level)%nodes==kept(level)%nodes).and.all(physical(level)%probes==kept(level)%probes), &
      'film hierarchy failure preserves complete level metadata')
    do p=1,size(physical(level)%films)
      call require(all(physical(level)%films(p)%position==kept(level)%films(p)%position).and. &
        physical(level)%films(p)%temperature==kept(level)%films(p)%temperature,'film hierarchy payload rollback')
    end do
  end do
  call payload()
  owners=modulo(owners+1,nranks)
  call migrate(100)
  call require(ok,'changed gas owners');call payload()
  call initialize_amr_patch_tree_3d([4,4,4],[integer::],patches(:0),coarse,ok)
  call require(ok,'coarse tree')
  call migrate_mpi_amr_wall_films(MPI_COMM_WORLD,coarse,[0],films,masses,zero,spacing,100,ok,phi)
  call require(ok.and.size(films)==merge(3*nranks,0,rank==0),'coarsened film migration')
  call payload()
  call migrate(100)
  call require(ok,'refinement reassignment');call payload()
  total=sum(masses,dim=2)
  call MPI_Allreduce(MPI_IN_PLACE,total,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call require(all(abs(total-initial)<=2e-12_dp*initial),'component inventory')
  deallocate(films,masses);allocate(films(0),masses(2,0))
  call migrate(0)
  call require(ok.and.size(films)==0.and.size(masses,2)==0,'globally empty films')
  call initialize_mpi_eb_film_levels(MPI_COMM_WORLD,tree,owners,films,masses,zero,spacing,0,phi,physical,ok)
  call require(ok.and.size(physical)==3,'globally empty film hierarchy')
  do level=1,3
    call require(size(physical(level)%films)==0.and.size(physical(level)%nodes)==0.and. &
      size(physical(level)%component_mass,2)==0.and.size(physical(level)%probes,2)==0,'empty film level arrays')
  end do
  if (rank==0) print *, 'PASS: MPI EB fluid-side stationary film migration'
  call MPI_Finalize(ierr)
contains
  function make_film(number) result(film)
    integer, intent(in) :: number
    type(spray_wall_film) :: film
    film=spray_wall_film([1._dp,.25_dp,.25_dp],[-1._dp,0._dp,0._dp], &
      number*1e-9_dp,300._dp+number,.001_dp*number,1e-6_dp*number,number)
    if (mod(number-1,3)==1) then
      film%position(1)=3;film%normal(1)=1
    else if (mod(number-1,3)==2) then
      film%position=[1._dp,3._dp,.25_dp];film%normal=[-1._dp,0._dp,0._dp]
    end if
  end function
  subroutine unchanged()
    call require(size(films)==size(saved).and.all(masses==saved_mass),'component/count rollback')
    call payload()
  end subroutine
  subroutine payload()
    integer :: q
    do q=1,size(films)
      expected=make_film(nint(films(q)%temperature)-300)
      call require(all(films(q)%position==expected%position).and.all(films(q)%normal==expected%normal), &
        'exact wall position and normal')
      call require(all([films(q)%mass,films(q)%temperature,films(q)%diameter,films(q)%height]== &
        [expected%mass,expected%temperature,expected%diameter,expected%height]),'film physical payload')
      call require(films(q)%material_id==expected%material_id,'film material')
      call require(all(masses(:,q)==expected%mass*[.25_dp,.75_dp]),'exact component payload')
    end do
  end subroutine
  subroutine migrate(capacity)
    integer, intent(in) :: capacity
    call migrate_mpi_amr_wall_films(MPI_COMM_WORLD,tree,owners,films,masses,zero,spacing,capacity,ok,phi)
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
