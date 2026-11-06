module spray_mpi_amr_tree_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use spray_eb_tree_stencil_mod, only: compatible_eb_nodal_lattice
  use spray_eb_intersection_mod, only: intersect_spray_level_set
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_integer_layout_3d
  use spray_parcel_mod, only: spray_parcel
  use spray_mixture_mod, only: valid_component_inventory
  use spray_amr_tree_partition_mod, only: partition_amr_tree_spray_parcels,locate_amr_tree_spray_cells, &
    locate_amr_tree_ghost_cells
  use spray_mpi_exchange_mod, only: exchange_spray_parcel_bins
  implicit none
  private
  public :: migrate_amr_tree_spray_parcels,migrate_eb_tree_spray_parcels
  public :: create_mpi_amr_virtual_particles,create_mpi_amr_ghost_particles,redistribute_mpi_amr_ghost_particles
  public :: create_mpi_amr_level_particles,collect_mpi_amr_level_particles
  public :: initialize_mpi_amr_sdc_populations,redistribute_mpi_amr_level_particles
  public :: refresh_mpi_amr_level_copies,redistribute_mpi_amr_subtree_particles
  type, public :: spray_level_population
    type(spray_parcel), allocatable :: parcels(:)
    real(dp), allocatable :: component_mass(:,:)
    integer, allocatable :: nodes(:)
  end type spray_level_population
contains
  ! Redistribute only synchronized physical levels >= minimum_level. A live
  ! parcel outside those valid boxes stays in the grown minimum-level boxes;
  ! lower logical levels are neither inspected nor changed.
  subroutine redistribute_mpi_amr_subtree_particles(comm,tree,owners,physical,lower,spacing,periodic, &
      minimum_level,grow,maximum_parcels,ok,nodal_phi)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),minimum_level,grow,maximum_parcels
    type(spray_level_population), intent(inout) :: physical(:)
    real(dp), intent(in) :: lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    type(spray_level_population), allocatable :: selected(:),candidate(:)
    type(spray_parcel), allocatable :: merged(:)
    real(dp), allocatable :: mass(:,:),positions(:,:)
    integer, allocatable :: locations(:,:),destination(:)
    integer :: levels,level,nc,p,q,status,flags(4),low(4),high(4)
    logical :: valid,agreed
    ok=.false.;valid=tree%is_valid()
    call agree();if(.not.agreed) return
    levels=tree%level_count()
    valid=size(physical)==levels.and.minimum_level>=0.and.minimum_level<levels.and. &
      grow>=0.and.grow<=1000000.and.maximum_parcels>=0
    call agree();if(.not.agreed) return
    flags=[minimum_level,grow,maximum_parcels,merge(1,0,present(nodal_phi))]
    call MPI_Allreduce(flags,low,4,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high,4,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    valid=allocated(physical(minimum_level+1)%component_mass)
    call agree();if(.not.agreed) return
    nc=size(physical(minimum_level+1)%component_mass,1)
    allocate(selected(levels));candidate=physical
    do level=1,levels
      if(level>minimum_level) then
        selected(level)=physical(level)
      else
        allocate(selected(level)%parcels(0),selected(level)%component_mass(nc,0))
      end if
    end do
    ! Global ID validation and temporary geometric migration are private.
    call collect_mpi_amr_level_particles(comm,tree,owners,selected,lower,spacing,periodic,1000000, &
      merged,mass,valid,nodal_phi)
    if(.not.valid) return
    allocate(positions(3,size(merged)),destination(size(merged)))
    do p=1,size(merged)
      positions(:,p)=merged(p)%position
    end do
    call locate_amr_tree_spray_cells(tree,positions,lower,spacing,periodic,locations,valid)
    call agree();if(.not.agreed) return
    do p=1,size(merged)
      destination(p)=max(minimum_level,tree%nodes(locations(1,p))%level)
    end do
    do level=minimum_level,levels-1
      candidate(level+1)%parcels=pack(merged,destination==level)
      if(allocated(candidate(level+1)%component_mass)) deallocate(candidate(level+1)%component_mass)
      allocate(candidate(level+1)%component_mass(nc,count(destination==level)));q=0
      do p=1,size(merged)
        if(destination(p)/=level) cycle
        q=q+1;candidate(level+1)%component_mass(:,q)=mass(:,p)
      end do
      call redistribute_mpi_amr_level_particles(comm,tree,owners,candidate(level+1)%parcels, &
        candidate(level+1)%component_mass,candidate(level+1)%nodes,lower,spacing,periodic,level, &
        merge(grow,0,level==minimum_level),maximum_parcels,valid,nodal_phi)
      if(.not.valid) return
    end do
    physical=candidate;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS) agreed=.false.
    end subroutine agree
  end subroutine redistribute_mpi_amr_subtree_particles

  ! Refresh copies before this level moves. Physical arrays carry logical
  ! level ownership until synchronization, even across valid-box boundaries.
  ! Descendants must be at this interval's start; next-level ghosts are made
  ! from the current level only and retained through its child substeps.
  subroutine refresh_mpi_amr_level_copies(comm,tree,owners,physical,lower,spacing,periodic,target_level,grow, &
      maximum_parcels,virtuals,next_ghosts,ok,nodal_phi)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),target_level,grow,maximum_parcels
    type(spray_level_population), intent(in) :: physical(:)
    real(dp), intent(in) :: lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    type(spray_level_population), intent(inout) :: virtuals,next_ghosts
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    type(spray_level_population), allocatable :: descendants(:)
    type(spray_level_population) :: coarse_copy,fine_copy
    integer :: levels,level,nc,status,ios,flags(5),low(5),high(5)
    logical :: valid,agreed
    ok=.false.;valid=tree%is_valid()
    call agree();if(.not.agreed) return
    levels=tree%level_count()
    valid=size(physical)==levels.and.target_level>=0.and.target_level<levels.and.grow>=0.and.maximum_parcels>=0
    call agree();if(.not.agreed) return
    flags=[levels,target_level,grow,maximum_parcels,merge(1,0,present(nodal_phi))]
    call MPI_Allreduce(flags,low,5,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high,5,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    nc=0
    do level=1,levels
      valid=allocated(physical(level)%parcels).and.allocated(physical(level)%component_mass)
      if(.not.valid) exit
      if(level==1) nc=size(physical(level)%component_mass,1)
      valid=nc>=1.and.nc<=16
      valid=valid.and.all(shape(physical(level)%component_mass)==[nc,size(physical(level)%parcels)])
      if(.not.valid) exit
    end do
    call agree();if(.not.agreed) return
    allocate(descendants(levels),stat=ios)
    valid=ios==0;call agree();if(.not.agreed) return
    do level=1,levels
      if(level>target_level+1) then
        descendants(level)=physical(level)
      else
        allocate(descendants(level)%parcels(0),descendants(level)%component_mass(nc,0),stat=ios)
        valid=ios==0;call agree();if(.not.agreed) return
      end if
    end do
    ! Collection validates unique live IDs; selection was by logical level,
    ! not by current geometric containment. Its temporary finest owners do
    ! not change which parcels are represented at the requested parent level.
    call collect_mpi_amr_level_particles(comm,tree,owners,descendants,lower,spacing,periodic,1000000, &
      coarse_copy%parcels,coarse_copy%component_mass,valid,nodal_phi)
    if(.not.valid) return
    call redistribute_mpi_amr_level_particles(comm,tree,owners,coarse_copy%parcels,coarse_copy%component_mass, &
      coarse_copy%nodes,lower,spacing,periodic,target_level,grow,maximum_parcels,valid,nodal_phi)
    if(.not.valid) return
    if(target_level+1<levels) then
      call route_mpi_amr_ghost_particles(comm,tree,owners,physical(target_level+1)%parcels, &
        physical(target_level+1)%component_mass,lower,spacing,periodic,target_level+1,grow,maximum_parcels, &
        fine_copy%parcels,fine_copy%component_mass,fine_copy%nodes,valid,.false.,nodal_phi)
      if(.not.valid) return
    else
      allocate(fine_copy%parcels(0),fine_copy%component_mass(nc,0),fine_copy%nodes(0),stat=ios)
      valid=ios==0;call agree();if(.not.agreed) return
    end if
    virtuals=coarse_copy;next_ghosts=fine_copy;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS) agreed=.false.
    end subroutine agree
  end subroutine refresh_mpi_amr_level_copies

  ! Initialize an outer step from a synchronized physical population. Every
  ! copy is made before any level moves. Outputs use level+1 indexing; only
  ! physical populations may be collected into the final physical inventory.
  subroutine initialize_mpi_amr_sdc_populations(comm,tree,owners,parcels,component_mass,lower,spacing,periodic, &
      grow,maximum_parcels,physical,virtuals,ghosts,ok,nodal_phi)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),grow,maximum_parcels
    type(spray_parcel), allocatable, intent(in) :: parcels(:)
    real(dp), allocatable, intent(in) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    type(spray_level_population), allocatable, intent(inout) :: physical(:),virtuals(:),ghosts(:)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    type(spray_level_population), allocatable :: active(:),coarse_copies(:),fine_copies(:)
    real(dp), allocatable :: positions(:,:)
    integer, allocatable :: locations(:,:)
    integer :: levels,level,p,nc,status,ios,flags(4),low(4),high(4)
    logical :: valid,agreed
    ok=.false.;valid=tree%is_valid().and.grow>=0.and.maximum_parcels>=0
    call agree();if(.not.agreed) return
    levels=tree%level_count()
    flags=[levels,grow,maximum_parcels,merge(1,0,present(nodal_phi))]
    call MPI_Allreduce(flags,low,4,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high,4,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    allocate(active(levels),coarse_copies(levels),fine_copies(levels),stat=ios)
    valid=ios==0;call agree();if(.not.agreed) return
    ! Select every physical level first; none of these calls changes the source.
    do level=0,levels-1
      call create_mpi_amr_level_particles(comm,tree,owners,parcels,component_mass,lower,spacing,periodic, &
        level,maximum_parcels,active(level+1)%parcels,active(level+1)%component_mass, &
        active(level+1)%nodes,valid,nodal_phi)
      if(.not.valid) return
    end do
    nc=size(component_mass,1)
    do level=0,levels-1
      call create_mpi_amr_virtual_particles(comm,tree,owners,parcels,component_mass,lower,spacing,periodic, &
        level,maximum_parcels,coarse_copies(level+1)%parcels,coarse_copies(level+1)%component_mass,valid,nodal_phi)
      if(.not.valid) return
      allocate(positions(3,size(coarse_copies(level+1)%parcels)),stat=ios)
      valid=ios==0;call agree();if(.not.agreed) return
      do p=1,size(positions,2)
        positions(:,p)=coarse_copies(level+1)%parcels(p)%position
      end do
      call locate_amr_tree_spray_cells(tree,positions,lower,spacing,periodic,locations,valid,maximum_level=level)
      call agree();if(.not.agreed) return
      coarse_copies(level+1)%nodes=locations(1,:)
      deallocate(positions)
      if(level==0) then
        allocate(fine_copies(1)%parcels(0),fine_copies(1)%component_mass(nc,0),fine_copies(1)%nodes(0),stat=ios)
        valid=ios==0;call agree();if(.not.agreed) return
      else
        call create_mpi_amr_ghost_particles(comm,tree,owners,parcels,component_mass,lower,spacing,periodic, &
          level,grow,maximum_parcels,fine_copies(level+1)%parcels,fine_copies(level+1)%component_mass, &
          fine_copies(level+1)%nodes,valid,nodal_phi)
        if(.not.valid) return
      end if
    end do
    call move_alloc(active,physical);call move_alloc(coarse_copies,virtuals);call move_alloc(fine_copies,ghosts)
    ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS) agreed=.false.
    end subroutine agree
  end subroutine initialize_mpi_amr_sdc_populations

  ! Merge physical level snapshots after their clocks synchronize. Ghost and
  ! virtual copies must not be supplied. ID-bin exchange detects duplicates
  ! globally before the private result migrates to current geometric owners.
  subroutine collect_mpi_amr_level_particles(comm,tree,owners,levels,lower,spacing,periodic,maximum_parcels, &
      parcels,component_mass,ok,nodal_phi)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),maximum_parcels
    type(spray_level_population), intent(in) :: levels(:)
    real(dp), intent(in) :: lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    type(spray_parcel), allocatable :: candidate(:)
    real(dp), allocatable :: mass(:,:)
    integer, allocatable :: counts(:),offsets(:),cursor(:)
    integer :: status,nranks,nc,level,p,r,q,flag,low,high
    integer(int64) :: total
    logical :: valid,agreed
    ok=.false.;nc=0;total=0
    call MPI_Comm_size(comm,nranks,status)
    if(status/=MPI_SUCCESS) return
    valid=tree%is_valid()
    if(valid) valid=size(levels)==tree%level_count()
    call agree();if(.not.agreed) return
    do level=1,size(levels)
      valid=allocated(levels(level)%parcels).and.allocated(levels(level)%component_mass)
      if(.not.valid) exit
      if(level==1) nc=size(levels(level)%component_mass,1)
      valid=nc>=1.and.nc<=16
      valid=valid.and.all(shape(levels(level)%component_mass)==[nc,size(levels(level)%parcels)])
      if(.not.valid) exit
      do p=1,size(levels(level)%parcels)
        valid=valid_component_inventory(levels(level)%parcels(p),levels(level)%component_mass(:,p))
        if(.not.valid) exit
      end do
      if(.not.valid) exit
      total=total+count(levels(level)%parcels%mass>0,kind=int64)
    end do
    valid=valid.and.total<=1000000
    call agree();if(.not.agreed) return
    flag=merge(1,0,present(nodal_phi))
    call MPI_Allreduce(flag,low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flag,high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.low/=high) return
    allocate(counts(nranks),offsets(nranks),cursor(nranks));counts=0
    do level=1,size(levels)
      do p=1,size(levels(level)%parcels)
        if(levels(level)%parcels(p)%mass<=0) cycle
        r=int(modulo(levels(level)%parcels(p)%id,int(nranks,int64)))+1
        counts(r)=counts(r)+1
      end do
    end do
    offsets(1)=0
    do r=2,nranks
      offsets(r)=offsets(r-1)+counts(r-1)
    end do
    cursor=offsets;allocate(candidate(int(total)),mass(nc,int(total)))
    do level=1,size(levels)
      do p=1,size(levels(level)%parcels)
        if(levels(level)%parcels(p)%mass<=0) cycle
        r=int(modulo(levels(level)%parcels(p)%id,int(nranks,int64)))+1
        cursor(r)=cursor(r)+1;q=cursor(r)
        candidate(q)=levels(level)%parcels(p);mass(:,q)=levels(level)%component_mass(:,p)
      end do
    end do
    call exchange_spray_parcel_bins(comm,candidate,counts,offsets,1000000,valid,mass)
    if(.not.valid) return
    if(present(nodal_phi)) then
      call migrate_eb_tree_spray_parcels(comm,tree,owners,nodal_phi,candidate,lower,spacing, &
        maximum_parcels,valid,mass,periodic=periodic)
    else
      call migrate_amr_tree_spray_parcels(comm,tree,owners,candidate,lower,spacing,periodic,maximum_parcels,valid,mass)
    end if
    if(.not.valid) return
    call move_alloc(candidate,parcels);call move_alloc(mass,component_mass);ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS) agreed=.false.
    end subroutine agree
  end subroutine collect_mpi_amr_level_particles

  ! Copies immediate-coarse physical parcels into the grown next-level boxes.
  ! Each source has one deterministic destination, even when grown boxes overlap.
  subroutine create_mpi_amr_ghost_particles(comm,tree,owners,parcels,component_mass,lower,spacing,periodic, &
      target_level,grow,maximum_parcels,ghosts,ghost_mass,ghost_nodes,ok,nodal_phi)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),target_level,grow,maximum_parcels
    type(spray_parcel), allocatable, intent(in) :: parcels(:)
    real(dp), allocatable, intent(in) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    type(spray_parcel), allocatable, intent(inout) :: ghosts(:)
    real(dp), allocatable, intent(inout) :: ghost_mass(:,:)
    integer, allocatable, intent(inout) :: ghost_nodes(:)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    call route_mpi_amr_ghost_particles(comm,tree,owners,parcels,component_mass,lower,spacing,periodic, &
      target_level,grow,maximum_parcels,ghosts,ghost_mass,ghost_nodes,ok,.true.,nodal_phi)
  end subroutine create_mpi_amr_ghost_particles

  ! Reassign an existing ghost population, retaining members now inside fine
  ! valid boxes and removing dead copies or copies outside all grown boxes.
  subroutine redistribute_mpi_amr_ghost_particles(comm,tree,owners,ghosts,ghost_mass,ghost_nodes,lower,spacing, &
      periodic,target_level,grow,maximum_parcels,ok,nodal_phi,tab_tangent_history)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),target_level,grow,maximum_parcels
    type(spray_parcel), allocatable, intent(inout) :: ghosts(:)
    real(dp), allocatable, intent(inout) :: ghost_mass(:,:)
    integer, allocatable, intent(inout) :: ghost_nodes(:)
    real(dp), intent(in) :: lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    real(dp), allocatable, optional, intent(inout) :: tab_tangent_history(:)
    type(spray_parcel), allocatable :: candidate(:)
    real(dp), allocatable :: mass(:,:)
    integer, allocatable :: nodes(:)
    call route_mpi_amr_ghost_particles(comm,tree,owners,ghosts,ghost_mass,lower,spacing,periodic, &
      target_level,grow,maximum_parcels,candidate,mass,nodes,ok,.false.,nodal_phi,tab_tangent_history=tab_tangent_history)
    if(.not.ok) return
    call move_alloc(candidate,ghosts);call move_alloc(mass,ghost_mass);call move_alloc(nodes,ghost_nodes)
  end subroutine redistribute_mpi_amr_ghost_particles

  ! Keep live physical parcels on their current level until synchronization.
  ! Crossing a valid-box edge is allowed inside its grow region; escaping all
  ! grown boxes rejects collectively rather than deleting physical inventory.
  subroutine redistribute_mpi_amr_level_particles(comm,tree,owners,parcels,component_mass,nodes,lower,spacing, &
      periodic,target_level,grow,maximum_parcels,ok,nodal_phi,tab_tangent_history)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),target_level,grow,maximum_parcels
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    integer, allocatable, intent(inout) :: nodes(:)
    real(dp), intent(in) :: lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    real(dp), allocatable, optional, intent(inout) :: tab_tangent_history(:)
    type(spray_parcel), allocatable :: candidate(:)
    real(dp), allocatable :: mass(:,:)
    integer, allocatable :: locations(:)
    call route_mpi_amr_ghost_particles(comm,tree,owners,parcels,component_mass,lower,spacing,periodic, &
      target_level,grow,maximum_parcels,candidate,mass,locations,ok,.false.,nodal_phi,retain_all=.true., &
      tab_tangent_history=tab_tangent_history)
    if(.not.ok) return
    call move_alloc(candidate,parcels);call move_alloc(mass,component_mass);call move_alloc(locations,nodes)
  end subroutine redistribute_mpi_amr_level_particles

  subroutine route_mpi_amr_ghost_particles(comm,tree,owners,parcels,component_mass,lower,spacing,periodic, &
      target_level,grow,maximum_parcels,ghosts,ghost_mass,ghost_nodes,ok,coarse_only,nodal_phi,retain_all,tab_tangent_history)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),target_level,grow,maximum_parcels
    type(spray_parcel), allocatable, intent(in) :: parcels(:)
    real(dp), allocatable, intent(in) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    type(spray_parcel), allocatable, intent(inout) :: ghosts(:)
    real(dp), allocatable, intent(inout) :: ghost_mass(:,:)
    integer, allocatable, intent(inout) :: ghost_nodes(:)
    logical, intent(in) :: coarse_only
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    real(dp), allocatable, optional, intent(inout) :: tab_tangent_history(:)
    real(dp), allocatable :: history(:),copies_history(:)
    type(spray_parcel), allocatable :: physical(:),copies(:)
    real(dp), allocatable :: mass(:,:),masses(:,:),positions(:,:)
    integer, allocatable :: found(:,:),active(:,:),counts(:),offsets(:),next(:),nodes(:)
    integer :: p,r,q,status,rank,nranks,flags(6),low(6),high(6)
    logical, optional, intent(in) :: retain_all
    logical :: keep_live
    logical :: valid,agreed
    ok=.false.;keep_live=.false.
    if(present(retain_all)) keep_live=retain_all
    valid=allocated(parcels).and.allocated(component_mass).and.tree%is_valid()
    if(valid) valid=target_level>=merge(0,1,keep_live.and..not.coarse_only).and. &
      target_level<tree%level_count().and.grow>=0.and.grow<=1000000
    call agree()
    if(.not.agreed) return
    flags=[target_level,grow,merge(1,0,present(nodal_phi)),merge(1,0,coarse_only),merge(1,0,keep_live), &
      merge(1,0,present(tab_tangent_history))]
    call MPI_Allreduce(flags,low,6,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high,6,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    ! Private physical migration validates replicated topology, bounds, payload
    ! and EB geometry before any selection or output publication.
    physical=parcels;mass=component_mass
    if(.not.coarse_only) then
      valid=size(component_mass,2)==size(parcels)
      if(valid) then
        do p=1,size(parcels)
          valid=valid_component_inventory(parcels(p),component_mass(:,p))
          if(.not.valid) exit
        end do
      end if
      call agree();if(.not.agreed) return
      physical=pack(parcels,parcels%mass>0)
      deallocate(mass);allocate(mass(size(component_mass,1),size(physical)));q=0
      do p=1,size(parcels)
        if(parcels(p)%mass<=0) cycle
        q=q+1;mass(:,q)=component_mass(:,p)
      end do
    end if
    if(present(tab_tangent_history)) then
      valid=allocated(tab_tangent_history)
      call agree();if(.not.agreed) return
      valid=size(tab_tangent_history)==size(parcels)
      if(valid) valid=all(ieee_is_finite(tab_tangent_history)).and.all(tab_tangent_history>=0)
      call agree();if(.not.agreed) return
      history=tab_tangent_history
      if(.not.coarse_only) history=pack(tab_tangent_history,parcels%mass>0)
    if(present(nodal_phi)) then
      call migrate_eb_tree_spray_parcels(comm,tree,owners,nodal_phi,physical,lower,spacing,1000000,valid,mass, &
        periodic=periodic,tab_tangent_history=history)
    else
      call migrate_amr_tree_spray_parcels(comm,tree,owners,physical,lower,spacing,periodic,1000000,valid,mass, &
        tab_tangent_history=history)
    end if
    else
    if(present(nodal_phi)) then
      call migrate_eb_tree_spray_parcels(comm,tree,owners,nodal_phi,physical,lower,spacing,1000000,valid,mass,periodic=periodic)
    else
      call migrate_amr_tree_spray_parcels(comm,tree,owners,physical,lower,spacing,periodic,1000000,valid,mass)
    end if
    endif
    if(.not.valid) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,status)
    if(status/=MPI_SUCCESS) return
    allocate(positions(3,size(physical)),counts(nranks),offsets(nranks),next(nranks));counts=0
    do p=1,size(physical)
      positions(:,p)=physical(p)%position
    end do
    call locate_amr_tree_spray_cells(tree,positions,lower,spacing,periodic,active,valid)
    call agree()
    if(.not.agreed) return
    call locate_level()
    if(valid.and.keep_live) valid=all(found(1,:)>0)
    call agree()
    if(.not.agreed) return
    do p=1,size(physical)
      if(coarse_only) then
        if(tree%nodes(active(1,p))%level/=target_level-1) found(1,p)=0
      end if
      if(found(1,p)==0) cycle
      r=owners(found(1,p))+1;counts(r)=counts(r)+1
    end do
    offsets(1)=0
    do r=2,nranks
      offsets(r)=offsets(r-1)+counts(r-1)
    end do
    next=offsets+1;allocate(copies(sum(counts)),masses(size(mass,1),sum(counts)))
    if(present(tab_tangent_history)) allocate(copies_history(sum(counts)))
    do p=1,size(physical)
      if(found(1,p)==0) cycle
      r=owners(found(1,p))+1;q=next(r);next(r)=q+1
      if(present(tab_tangent_history)) copies_history(q)=history(p)
      copies(q)=physical(p);masses(:,q)=mass(:,p)
    end do
    if(present(tab_tangent_history)) then
      call exchange_spray_parcel_bins(comm,copies,counts,offsets,maximum_parcels,valid,masses,copies_history)
    else
      call exchange_spray_parcel_bins(comm,copies,counts,offsets,maximum_parcels,valid,masses)
    endif
    if(.not.valid) return
    deallocate(positions);allocate(positions(3,size(copies)))
    do p=1,size(copies)
      positions(:,p)=copies(p)%position
    end do
    call locate_level()
    if(valid) then
      do p=1,size(copies)
        if(found(1,p)==0) then
          valid=.false.;exit
        end if
        if(owners(found(1,p))/=rank) valid=.false.
      end do
    end if
    call agree()
    if(.not.agreed) return
    if(present(tab_tangent_history)) call move_alloc(copies_history,tab_tangent_history)
    nodes=found(1,:)
    call move_alloc(copies,ghosts);call move_alloc(masses,ghost_mass);call move_alloc(nodes,ghost_nodes)
    ok=.true.
  contains
    subroutine locate_level()
      if(target_level==0) then
        call locate_amr_tree_spray_cells(tree,positions,lower,spacing,periodic,found,valid,maximum_level=0)
      else
        call locate_amr_tree_ghost_cells(tree,positions,lower,spacing,periodic,target_level,grow,found,valid)
      end if
    end subroutine locate_level
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine route_mpi_amr_ghost_particles

  ! Build nonaggregated coarse representations of all finer active parcels.
  ! The separate output retains source IDs; never merge it into physical parcels.
  ! Caller publishes/clears this level-local population with its stage lifecycle.
  subroutine create_mpi_amr_virtual_particles(comm,tree,owners,parcels,component_mass,lower,spacing,periodic, &
      target_level,maximum_parcels,virtuals,virtual_mass,ok,nodal_phi)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),target_level,maximum_parcels
    type(spray_parcel), allocatable, intent(in) :: parcels(:)
    real(dp), allocatable, intent(in) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    type(spray_parcel), allocatable, intent(inout) :: virtuals(:)
    real(dp), allocatable, intent(inout) :: virtual_mass(:,:)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    call select_mpi_amr_level_copies(comm,tree,owners,parcels,component_mass,lower,spacing,periodic, &
      target_level,maximum_parcels,virtuals,virtual_mass,ok,.true.,nodal_phi)
  end subroutine create_mpi_amr_virtual_particles

  ! Select physical parcels whose finest containing node is exactly this level.
  ! Source data remain read-only; the caller must not merge these copies with it.
  subroutine create_mpi_amr_level_particles(comm,tree,owners,parcels,component_mass,lower,spacing,periodic, &
      target_level,maximum_parcels,level_parcels,level_mass,level_nodes,ok,nodal_phi)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),target_level,maximum_parcels
    type(spray_parcel), allocatable, intent(in) :: parcels(:)
    real(dp), allocatable, intent(in) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    type(spray_parcel), allocatable, intent(inout) :: level_parcels(:)
    real(dp), allocatable, intent(inout) :: level_mass(:,:)
    integer, allocatable, intent(inout) :: level_nodes(:)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    type(spray_parcel), allocatable :: candidate(:)
    real(dp), allocatable :: mass(:,:),positions(:,:)
    integer, allocatable :: found(:,:),nodes(:)
    integer :: p,status
    logical :: valid,agreed
    ok=.false.
    call select_mpi_amr_level_copies(comm,tree,owners,parcels,component_mass,lower,spacing,periodic, &
      target_level,maximum_parcels,candidate,mass,valid,.false.,nodal_phi)
    if(.not.valid) return
    allocate(positions(3,size(candidate)))
    do p=1,size(candidate)
      positions(:,p)=candidate(p)%position
    end do
    call locate_amr_tree_spray_cells(tree,positions,lower,spacing,periodic,found,valid)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    nodes=found(1,:)
    call move_alloc(candidate,level_parcels);call move_alloc(mass,level_mass);call move_alloc(nodes,level_nodes)
    ok=.true.
  end subroutine create_mpi_amr_level_particles

  subroutine select_mpi_amr_level_copies(comm,tree,owners,parcels,component_mass,lower,spacing,periodic, &
      target_level,maximum_parcels,virtuals,virtual_mass,ok,finer_copies,nodal_phi)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),target_level,maximum_parcels
    type(spray_parcel), allocatable, intent(in) :: parcels(:)
    real(dp), allocatable, intent(in) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    type(spray_parcel), allocatable, intent(inout) :: virtuals(:)
    real(dp), allocatable, intent(inout) :: virtual_mass(:,:)
    logical, intent(in) :: finer_copies
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    type(spray_parcel), allocatable :: copies(:)
    real(dp), allocatable :: masses(:,:),positions(:,:)
    integer, allocatable :: locations(:,:)
    logical, allocatable :: selected(:)
    integer :: p,q,nc,status,flags(3),low(3),high(3)
    logical :: valid,agreed
    ok=.false.
    valid=tree%is_valid().and.allocated(parcels).and.allocated(component_mass)
    if(valid) valid=target_level>=0.and.target_level<tree%level_count()
    if(valid) valid=size(component_mass,2)==size(parcels)
    if(valid) valid=size(component_mass,1)>=1.and.size(component_mass,1)<=16
    call agree()
    if(.not.agreed) return
    flags=[target_level,merge(1,0,present(nodal_phi)),merge(1,0,finer_copies)]
    call MPI_Allreduce(flags,low,3,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high,3,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    nc=size(component_mass,1)
    allocate(positions(3,size(parcels)),selected(size(parcels)))
    do p=1,size(parcels)
      valid=valid.and.valid_component_inventory(parcels(p),component_mass(:,p))
      if(p>1) valid=valid.and.parcels(p)%id>parcels(p-1)%id
      positions(:,p)=parcels(p)%position
    end do
    call agree()
    if(.not.agreed) return
    call locate_amr_tree_spray_cells(tree,positions,lower,spacing,periodic,locations,valid)
    call agree()
    if(.not.agreed) return
    do p=1,size(parcels)
      if(finer_copies) then
        selected(p)=tree%nodes(locations(1,p))%level>target_level
      else
        selected(p)=tree%nodes(locations(1,p))%level==target_level
      end if
    end do
    copies=pack(parcels,selected);allocate(masses(nc,count(selected)));q=0
    do p=1,size(parcels)
      if(.not.selected(p)) cycle
      q=q+1;masses(:,q)=component_mass(:,p)
    end do
    if(present(nodal_phi)) then
      call migrate_eb_tree_spray_parcels(comm,tree,owners,nodal_phi,copies,lower,spacing, &
        maximum_parcels,valid,masses,target_level,periodic)
    else
      call migrate_amr_tree_spray_parcels(comm,tree,owners,copies,lower,spacing,periodic, &
        maximum_parcels,valid,masses,target_level)
    end if
    if(.not.valid) return
    call move_alloc(copies,virtuals);call move_alloc(masses,virtual_mass)
    ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine select_mpi_amr_level_copies
  ! The replicated nodal field describes the root domain at the finest tree
  ! resolution or a common integer refinement of it, as used by the gas EB hierarchy builder.
  ! Optional periodic axes wrap before EB checking and require matching nodal seams.
  subroutine migrate_eb_tree_spray_parcels(comm,tree,owners,phi,parcels,lower,spacing, &
      maximum_parcels,ok,component_mass,maximum_level,periodic,tab_tangent_history)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, optional, intent(in) :: maximum_level
    integer, intent(in) :: owners(:),maximum_parcels
    real(dp), intent(in) :: phi(0:,0:,0:),lower(3),spacing(3)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, optional, intent(inout) :: component_mass(:,:)
    logical, optional, intent(in) :: periodic(3)
    logical :: wrap(3)
    logical, intent(out) :: ok
    real(dp), allocatable, optional, intent(inout) :: tab_tangent_history(:)
    real(dp), allocatable :: history(:)
    type(spray_parcel), allocatable :: candidate(:)
    real(dp), allocatable :: masses(:,:),reference(:,:,:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer :: dims(3),low(3),high(3),ios,ierr,p,axis,cell(3),history_flag,history_min,history_max
    real(dp) :: upper(3),fraction,position(3),normal(3)
    logical :: valid,agreed,hit
    ok=.false.;dims=shape(phi)-1;wrap=.false.
    if(present(periodic)) wrap=periodic
    valid=allocated(parcels).and.all(dims>=1).and.all(dims<=1000000)
    if(present(component_mass)) valid=valid.and.allocated(component_mass)
    call agree()
    if(.not.agreed) return
    valid=product(int(dims,int64))<=1000000_int64
    valid=valid.and.all(ieee_is_finite(phi))
    if(valid) valid=maxval(abs(phi))<=1e100_dp
    call agree()
    if(.not.agreed) return
    valid=.true.
    if(wrap(1)) valid=valid.and.all(phi(0,:,:)==phi(dims(1),:,:))
    if(wrap(2)) valid=valid.and.all(phi(:,0,:)==phi(:,dims(2),:))
    if(wrap(3)) valid=valid.and.all(phi(:,:,0)==phi(:,:,dims(3)))
    call agree();if(.not.agreed) return
    history_flag=merge(1,0,present(tab_tangent_history))
    call MPI_Allreduce(history_flag,history_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(history_flag,history_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.history_min/=history_max) return
    if(present(tab_tangent_history)) then
      valid=allocated(tab_tangent_history)
      call agree();if(.not.agreed) return
    candidate=parcels
    history=tab_tangent_history
    if(present(component_mass)) then
      masses=component_mass
      call migrate_amr_tree_spray_parcels(comm,tree,owners,candidate,lower,spacing, &
        wrap,maximum_parcels,valid,masses,maximum_level,history)
    else
      call migrate_amr_tree_spray_parcels(comm,tree,owners,candidate,lower,spacing, &
        wrap,maximum_parcels,valid,maximum_level=maximum_level,tab_tangent_history=history)
    end if
    else
    candidate=parcels
    if(present(component_mass)) then
      masses=component_mass
      call migrate_amr_tree_spray_parcels(comm,tree,owners,candidate,lower,spacing, &
        wrap,maximum_parcels,valid,masses,maximum_level)
    else
      call migrate_amr_tree_spray_parcels(comm,tree,owners,candidate,lower,spacing, &
        wrap,maximum_parcels,valid,maximum_level=maximum_level)
    end if
    endif
    if(.not.valid) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    if(valid) valid=compatible_eb_nodal_lattice(dims,domain(:,size(domain,2)))
    call agree()
    if(.not.agreed) return
    call MPI_Allreduce(dims,low,3,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(dims,high,3,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(low/=high)) return
    allocate(reference(dims(1)+1,dims(2)+1,dims(3)+1),stat=ios)
    valid=ios==0;call agree()
    if(.not.agreed) return
    reference=phi
    call MPI_Bcast(reference,size(reference),MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    valid=all(reference==phi)
    upper=lower+spacing*real(tree%nodes(1)%shape,dp)
    valid=valid.and.all((upper-lower)/real(dims,dp)>=1e-12_dp).and.maxval(abs([lower,upper]))<=1e9_dp
    call agree()
    if(.not.agreed) return
    hit=.false.;fraction=0;position=0;normal=0;cell=0
    do p=1,size(candidate)
      ! General AMR migration wraps the ownership lookup but preserves payload
      ! coordinates. EB intersection requires the canonical physical domain.
      do axis=1,3
        if(.not.wrap(axis)) cycle
        candidate(p)%position(axis)=lower(axis)+modulo(candidate(p)%position(axis)-lower(axis),upper(axis)-lower(axis))
      end do
      call intersect_spray_level_set(phi,lower,upper,candidate(p)%position,candidate(p)%position, &
        hit,fraction,position,normal,cell,valid)
      if(.not.valid) exit
    end do
    call agree()
    if(.not.agreed) return
    if(present(tab_tangent_history)) call move_alloc(history,tab_tangent_history)
    call move_alloc(candidate,parcels)
    if(present(component_mass)) call move_alloc(masses,component_mass)
    ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine
  subroutine migrate_amr_tree_spray_parcels(comm,tree,owners,parcels,lower,spacing,periodic, &
      maximum_parcels,ok,component_mass,maximum_level,tab_tangent_history)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, optional, intent(in) :: maximum_level
    integer, intent(in) :: owners(:),maximum_parcels
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), intent(in) :: lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    logical, intent(out) :: ok
    real(dp), allocatable, optional, intent(inout) :: component_mass(:,:)
    real(dp), allocatable, optional, intent(inout) :: tab_tangent_history(:)
    real(dp), allocatable :: history(:)
    type(spray_parcel), allocatable :: candidate(:)
    real(dp), allocatable :: masses(:,:),positions(:,:)
    integer, allocatable :: counts(:),offsets(:),locations(:,:)
    real(dp) :: bounds(6),reference_bounds(6)
    integer :: rank,nranks,ierr,nodes,low,high,node,p,ios,descriptor(13),reference(13),flags(6),reference_flags(6),ceiling
    logical :: valid,agreed
    ok=.false.
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=tree%is_valid().and.allocated(parcels)
    valid=valid.and.all(ieee_is_finite([lower,spacing]))
    if (valid) valid=size(owners)==size(tree%nodes).and.size(parcels)<=1000000
    if (valid) valid=all(owners>=0).and.all(owners<nranks)
    if (present(component_mass)) valid=valid.and.allocated(component_mass)
    if(present(tab_tangent_history)) then
      valid=valid.and.allocated(tab_tangent_history)
      if(allocated(tab_tangent_history).and.allocated(parcels)) then
        valid=valid.and.size(tab_tangent_history)==size(parcels)
        valid=valid.and.all(ieee_is_finite(tab_tangent_history)).and.all(tab_tangent_history>=0)
      endif
    endif
    call agree()
    if (.not.agreed) return
    ceiling=tree%level_count()-1
    if(present(maximum_level)) ceiling=maximum_level
    valid=ceiling>=0.and.ceiling<tree%level_count()
    call agree()
    if(.not.agreed) return
    nodes=size(tree%nodes)
    call MPI_Allreduce(nodes,low,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(nodes,high,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.low/=high) return
    bounds=[lower,spacing];reference_bounds=bounds
    call MPI_Bcast(reference_bounds,6,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    flags=[merge(1,0,periodic),merge(1,0,present(component_mass)),ceiling, &
      merge(1,0,present(tab_tangent_history))]
    reference_flags=flags
    call MPI_Bcast(reference_flags,6,MPI_INTEGER,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=all(bounds==reference_bounds).and.all(flags==reference_flags)
    call agree()
    if (.not.agreed) return
    do node=1,nodes
      associate(n=>tree%nodes(node),patch=>tree%nodes(node)%patch)
      descriptor=[owners(node),n%parent,n%level,n%shape,patch%refinement_ratio, &
        patch%coarse_i_lower,patch%coarse_i_upper,patch%coarse_j_lower,patch%coarse_j_upper, &
        patch%coarse_k_lower,patch%coarse_k_upper]
      end associate
      reference=descriptor
      call MPI_Bcast(reference,13,MPI_INTEGER,0,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      valid=valid.and.all(reference==descriptor)
    end do
    call agree()
    if (.not.agreed) return
    if(present(tab_tangent_history)) then
    candidate=parcels
    history=tab_tangent_history
    if (present(component_mass)) then
      masses=component_mass
      call partition_amr_tree_spray_parcels(tree,owners,nranks,candidate,lower,spacing,periodic, &
        counts,offsets,locations,valid,masses,maximum_level,history)
    else
      call partition_amr_tree_spray_parcels(tree,owners,nranks,candidate,lower,spacing,periodic, &
        counts,offsets,locations,valid,maximum_level=maximum_level,tab_tangent_history=history)
    end if
    call agree()
    if (.not.agreed) return
    if (present(component_mass)) then
      call exchange_spray_parcel_bins(comm,candidate,counts,offsets,maximum_parcels,valid,masses,history)
    else
      call exchange_spray_parcel_bins(comm,candidate,counts,offsets,maximum_parcels,valid,tab_tangent_history=history)
    end if
    else
    candidate=parcels
    if (present(component_mass)) then
      masses=component_mass
      call partition_amr_tree_spray_parcels(tree,owners,nranks,candidate,lower,spacing,periodic, &
        counts,offsets,locations,valid,masses,maximum_level)
    else
      call partition_amr_tree_spray_parcels(tree,owners,nranks,candidate,lower,spacing,periodic, &
        counts,offsets,locations,valid,maximum_level=maximum_level)
    end if
    call agree()
    if (.not.agreed) return
    if (present(component_mass)) then
      call exchange_spray_parcel_bins(comm,candidate,counts,offsets,maximum_parcels,valid,masses)
    else
      call exchange_spray_parcel_bins(comm,candidate,counts,offsets,maximum_parcels,valid)
    end if
    endif
    if (.not.valid) return
    allocate(positions(3,size(candidate)),stat=ios)
    valid=ios==0;call agree()
    if (.not.agreed) return
    do p=1,size(candidate)
      positions(:,p)=candidate(p)%position
    end do
    call locate_amr_tree_spray_cells(tree,positions,lower,spacing,periodic,locations,valid,maximum_level)
    if (valid) then
      do p=1,size(candidate)
        if (owners(locations(1,p))/=rank) valid=.false.
      end do
    end if
    call agree()
    if (.not.agreed) return
    if(present(tab_tangent_history)) call move_alloc(history,tab_tangent_history)
    call move_alloc(candidate,parcels)
    if (present(component_mass)) call move_alloc(masses,component_mass)
    ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine
end module
