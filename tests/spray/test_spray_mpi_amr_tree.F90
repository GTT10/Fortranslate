program test_spray_mpi_amr_tree
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,initialize_amr_patch_tree_3d
  use spray_parcel_mod, only: spray_parcel
  use spray_mpi_amr_tree_mod, only: migrate_amr_tree_spray_parcels,create_mpi_amr_virtual_particles, &
    create_mpi_amr_ghost_particles,redistribute_mpi_amr_ghost_particles,create_mpi_amr_level_particles, &
    spray_level_population,collect_mpi_amr_level_particles,initialize_mpi_amr_sdc_populations
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
  call check_physical_levels()
  call check_sdc_populations()
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
  call check_ghost_redistribution()
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
  subroutine check_sdc_populations()
    type(spray_level_population), allocatable :: active(:),coarse_copies(:),fine_copies(:)
    type(spray_level_population), allocatable :: old_active(:),old_coarse(:),old_fine(:)
    integer :: level,a,seen(2*nranks)
    call initialize_mpi_amr_sdc_populations(MPI_COMM_WORLD,tree,owners,saved,saved_mass,zero,spacing, &
      [.true.,.false.,.false.],8,100,active,coarse_copies,fine_copies,ok)
    call require(ok,'initialize synchronized SDC populations')
    call require(size(active)==3.and.size(coarse_copies)==3.and.size(fine_copies)==3,'population level count')
    seen=0
    do level=1,3
      do a=1,size(active(level)%parcels)
        seen(int(active(level)%parcels(a)%id))=seen(int(active(level)%parcels(a)%id))+1
      end do
      call check_population(active(level),level-1)
      call check_population(coarse_copies(level),level-1)
      call check_population(fine_copies(level),level-1)
    end do
    call MPI_Allreduce(MPI_IN_PLACE,seen,size(seen),MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(all(seen==1),'initialized physical population has each ID once')
    call require(size(active(1)%parcels)==0.and.size(fine_copies(1)%parcels)==0,'root physical and ghosts empty')
    call require(size(coarse_copies(1)%parcels)==merge(2*nranks,0,rank==0),'all physical descendants copied to root')
    call require(size(coarse_copies(2)%parcels)==merge(nranks,0,rank==0),'fine descendants copied to parent')
    call require(size(coarse_copies(3)%parcels)==0,'no virtuals at finest level')
    call require(size(fine_copies(3)%parcels)==merge(nranks,0,rank==nranks-1),'immediate coarse ghost population')
    old_active=active;old_coarse=coarse_copies;old_fine=fine_copies
    ! Physical level capacity fits; root virtual copies exceed it afterwards.
    call initialize_mpi_amr_sdc_populations(MPI_COMM_WORLD,tree,owners,saved,saved_mass,zero,spacing, &
      [.true.,.false.,.false.],8,nranks,active,coarse_copies,fine_copies,ok)
    call require(.not.ok,'late virtual capacity rejects initialization')
    do level=1,3
      call compare_population(active(level),old_active(level))
      call compare_population(coarse_copies(level),old_coarse(level))
      call compare_population(fine_copies(level),old_fine(level))
    end do
    call require(all(parcels%id==saved%id).and.all(masses==saved_mass),'SDC initialization retains physical input')
  end subroutine check_sdc_populations
  subroutine check_population(population,level)
    type(spray_level_population), intent(in) :: population
    integer, intent(in) :: level
    integer :: a
    do a=1,size(population%parcels)
      expected=make_parcel(int(population%parcels(a)%id))
      call require(tree%nodes(population%nodes(a))%level==level.and.owners(population%nodes(a))==rank, &
        'initialized population level and owner')
      call require(all(population%parcels(a)%position==expected%position).and. &
        all(population%parcels(a)%velocity==expected%velocity).and. &
        population%parcels(a)%temperature==expected%temperature,'population snapshot before motion')
      call require(all(population%component_mass(:,a)==expected%mass*[.25_dp,.75_dp]),'population components')
    end do
  end subroutine check_population
  subroutine compare_population(population,old)
    type(spray_level_population), intent(in) :: population,old
    call require(size(population%parcels)==size(old%parcels),'population count rollback')
    call require(all(population%parcels%id==old%parcels%id).and.all(population%nodes==old%nodes).and. &
      all(population%component_mass==old%component_mass),'population identity and components rollback')
  end subroutine compare_population
  subroutine check_physical_levels()
    type(spray_parcel), allocatable :: selected(:),before(:)
    type(spray_level_population) :: populations(3)
    real(dp), allocatable :: selected_mass(:,:),before_mass(:,:)
    integer, allocatable :: nodes(:),before_nodes(:)
    integer :: seen(2*nranks),level,a
    seen=0
    do level=0,2
      call create_mpi_amr_level_particles(MPI_COMM_WORLD,tree,owners,saved,saved_mass,zero,spacing, &
        [.true.,.false.,.false.],level,100,selected,selected_mass,nodes,ok)
      call require(ok,'physical level selection')
      populations(level+1)%parcels=selected;populations(level+1)%component_mass=selected_mass
      populations(level+1)%nodes=nodes
      if(level==0) call require(size(selected)==0,'covered root contains no physical population')
      do a=1,size(selected)
        seen(int(selected(a)%id))=seen(int(selected(a)%id))+1
        call require(tree%nodes(nodes(a))%level==level.and.owners(nodes(a))==rank,'exact physical level and owner')
        expected=make_parcel(int(selected(a)%id))
        call require(all(selected(a)%position==expected%position).and. &
          all(selected(a)%velocity==expected%velocity).and.selected(a)%temperature==expected%temperature, &
          'level selection preserves physical state')
        call require(all(selected_mass(:,a)==expected%mass*[.25_dp,.75_dp]),'physical level component inventory')
      end do
    end do
    call MPI_Allreduce(MPI_IN_PLACE,seen,size(seen),MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(all(seen==1),'every physical ID occurs in exactly one level population')
    before=selected;before_mass=selected_mass;before_nodes=nodes
    call create_mpi_amr_level_particles(MPI_COMM_WORLD,tree,owners,saved,saved_mass,zero,spacing, &
      [.true.,.false.,.false.],1,0,selected,selected_mass,nodes,ok)
    call require(.not.ok.and.size(selected)==size(before),'physical level capacity failure')
    call require(all(selected%id==before%id).and.all(selected_mass==before_mass).and.all(nodes==before_nodes), &
      'physical level outputs roll back together')
    call require(all(parcels%id==saved%id).and.all(masses==saved_mass),'physical source population unchanged')
    ! Move live populations across levels before synchronized collection.
    do level=1,3
      do a=1,size(populations(level)%parcels)
        populations(level)%parcels(a)%position(1)=merge(3.5_dp,.25_dp,level==3)
      end do
    end do
    call collect_mpi_amr_level_particles(MPI_COMM_WORLD,tree,owners,populations,zero,spacing, &
      [.true.,.false.,.false.],0,selected,selected_mass,ok)
    call require(.not.ok.and.all(selected%id==before%id).and.all(selected_mass==before_mass),'collection capacity rollback')
    call collect_mpi_amr_level_particles(MPI_COMM_WORLD,tree,owners,populations,zero,spacing, &
      [.true.,.false.,.false.],100,selected,selected_mass,ok)
    call require(ok,'collect synchronized physical level populations')
    seen=0
    do a=1,size(selected)
      seen(int(selected(a)%id))=seen(int(selected(a)%id))+1
      expected=make_parcel(int(selected(a)%id))
      call require(all(selected_mass(:,a)==expected%mass*[.25_dp,.75_dp]),'collection component conservation')
      if(mod(selected(a)%id,2_int64)==1) then
        call require(rank==owners(3).and.selected(a)%position(1)==3.5_dp,'collected fine-to-coarse crossing')
      else
        call require(rank==owners(4).and.selected(a)%position(1)==.25_dp,'collected coarse-to-fine crossing')
      end if
    end do
    call MPI_Allreduce(MPI_IN_PLACE,seen,size(seen),MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(all(seen==1),'collection retains every physical ID once')
    before=selected;before_mass=selected_mass
    if(rank==0) then
      deallocate(populations(1)%parcels,populations(1)%component_mass)
      allocate(populations(1)%parcels(1),populations(1)%component_mass(2,1))
      populations(1)%parcels(1)=make_parcel(1)
      populations(1)%component_mass(:,1)=populations(1)%parcels(1)%mass*[.25_dp,.75_dp]
    end if
    call collect_mpi_amr_level_particles(MPI_COMM_WORLD,tree,owners,populations,zero,spacing, &
      [.true.,.false.,.false.],100,selected,selected_mass,ok)
    call require(.not.ok.and.all(selected%id==before%id).and.all(selected_mass==before_mass), &
      'duplicate physical ID at different destinations rejects atomically')

  end subroutine check_physical_levels
  subroutine check_ghost_redistribution()
    type(spray_parcel), allocatable :: moved(:),before_move(:)
    real(dp), allocatable :: moved_mass(:,:),before_mass(:,:)
    integer, allocatable :: nodes(:)
    integer :: map(4),a
    moved=ghosts;moved_mass=ghost_mass;nodes=ghost_nodes;map=owners;map(4)=0
    do a=1,size(moved)
      select case(mod(moved(a)%id,4_int64))
      case(1)
        moved(a)%mass=0;moved_mass(:,a)=0;moved(a)%position=4._dp
      case(2)
        moved(a)%position(1)=2._dp
      case(3)
        moved(a)%position(1)=.875_dp
      end select
    end do
    before_move=moved;before_mass=moved_mass
    call redistribute_mpi_amr_ghost_particles(MPI_COMM_WORLD,tree,map,moved,moved_mass,nodes,zero,spacing, &
      [.true.,.false.,.false.],2,1,0,ok)
    call require(.not.ok.and.size(moved)==size(before_move),'ghost redistribution capacity rollback')
    call require(all(moved_mass==before_mass).and.all(nodes==4),'ghost redistribution preserves components/nodes on failure')
    do a=1,size(moved)
      call require(all(moved(a)%position==before_move(a)%position),'ghost failure preserves coordinates')
    end do
    call redistribute_mpi_amr_ghost_particles(MPI_COMM_WORLD,tree,map,moved,moved_mass,nodes,zero,spacing, &
      [.true.,.false.,.false.],2,1,100,ok)
    call require(ok.and.size(moved)==merge(nranks,0,rank==0),'existing ghosts reassign to changed fine owner')
    call require(all(nodes==4),'redistributed fine node indices')
    do a=1,size(moved)
      call require(mod(moved(a)%id,4_int64)==0.or.mod(moved(a)%id,4_int64)==3,'dead and escaped ghosts removed')
      expected=make_parcel(int(moved(a)%id))
      call require(all(moved_mass(:,a)==expected%mass*[.25_dp,.75_dp]),'surviving ghost components retained')
      if(a>1) call require(moved(a)%id>moved(a-1)%id,'redistributed ghost IDs sorted')
    end do
    call redistribute_mpi_amr_ghost_particles(MPI_COMM_WORLD,tree,map,moved,moved_mass,nodes,zero,spacing, &
      [.true.,.false.,.false.],2,0,100,ok)
    call require(ok.and.size(moved)==merge(nranks/2,0,rank==0),'zero-grow retains ghosts now inside fine valid cells')
    do a=1,size(moved)
      call require(moved(a)%position(1)==.875_dp,'valid-cell ghost coordinates retained')
    end do
  end subroutine check_ghost_redistribution
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
