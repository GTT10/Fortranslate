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
  use spray_amr_tree_partition_mod, only: partition_amr_tree_spray_parcels,locate_amr_tree_spray_cells
  use spray_mpi_exchange_mod, only: exchange_spray_parcel_bins
  implicit none
  private
  public :: migrate_amr_tree_spray_parcels,migrate_eb_tree_spray_parcels
  public :: create_mpi_amr_virtual_particles
contains
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
    type(spray_parcel), allocatable :: copies(:)
    real(dp), allocatable :: masses(:,:),positions(:,:)
    integer, allocatable :: locations(:,:)
    logical, allocatable :: selected(:)
    integer :: p,q,nc,status,flags(2),low(2),high(2)
    logical :: valid,agreed
    ok=.false.
    valid=tree%is_valid().and.allocated(parcels).and.allocated(component_mass)
    if(valid) valid=target_level>=0.and.target_level<tree%level_count()
    if(valid) valid=size(component_mass,2)==size(parcels)
    if(valid) valid=size(component_mass,1)>=1.and.size(component_mass,1)<=16
    if(present(nodal_phi)) valid=valid.and..not.any(periodic)
    call agree()
    if(.not.agreed) return
    flags=[target_level,merge(1,0,present(nodal_phi))]
    call MPI_Allreduce(flags,low,2,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high,2,MPI_INTEGER,MPI_MAX,comm,status)
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
      selected(p)=tree%nodes(locations(1,p))%level>target_level
    end do
    copies=pack(parcels,selected);allocate(masses(nc,count(selected)));q=0
    do p=1,size(parcels)
      if(.not.selected(p)) cycle
      q=q+1;masses(:,q)=component_mass(:,p)
    end do
    if(present(nodal_phi)) then
      call migrate_eb_tree_spray_parcels(comm,tree,owners,nodal_phi,copies,lower,spacing, &
        maximum_parcels,valid,masses,target_level)
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
  end subroutine create_mpi_amr_virtual_particles
  ! The replicated nodal field describes the root domain at the finest tree
  ! resolution or a common integer refinement of it, as used by the gas EB hierarchy builder. No periodic wrapping.
  subroutine migrate_eb_tree_spray_parcels(comm,tree,owners,phi,parcels,lower,spacing, &
      maximum_parcels,ok,component_mass,maximum_level)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, optional, intent(in) :: maximum_level
    integer, intent(in) :: owners(:),maximum_parcels
    real(dp), intent(in) :: phi(0:,0:,0:),lower(3),spacing(3)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, optional, intent(inout) :: component_mass(:,:)
    logical, intent(out) :: ok
    type(spray_parcel), allocatable :: candidate(:)
    real(dp), allocatable :: masses(:,:),reference(:,:,:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer :: dims(3),low(3),high(3),ios,ierr,p,cell(3)
    real(dp) :: upper(3),fraction,position(3),normal(3)
    logical :: valid,agreed,hit
    ok=.false.;dims=shape(phi)-1
    valid=allocated(parcels).and.all(dims>=1).and.all(dims<=1000000)
    if(present(component_mass)) valid=valid.and.allocated(component_mass)
    call agree()
    if(.not.agreed) return
    valid=product(int(dims,int64))<=1000000_int64
    valid=valid.and.all(ieee_is_finite(phi))
    if(valid) valid=maxval(abs(phi))<=1e100_dp
    call agree()
    if(.not.agreed) return
    candidate=parcels
    if(present(component_mass)) then
      masses=component_mass
      call migrate_amr_tree_spray_parcels(comm,tree,owners,candidate,lower,spacing, &
        [.false.,.false.,.false.],maximum_parcels,valid,masses,maximum_level)
    else
      call migrate_amr_tree_spray_parcels(comm,tree,owners,candidate,lower,spacing, &
        [.false.,.false.,.false.],maximum_parcels,valid,maximum_level=maximum_level)
    end if
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
      call intersect_spray_level_set(phi,lower,upper,candidate(p)%position,candidate(p)%position, &
        hit,fraction,position,normal,cell,valid)
      if(.not.valid) exit
    end do
    call agree()
    if(.not.agreed) return
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
      maximum_parcels,ok,component_mass,maximum_level)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, optional, intent(in) :: maximum_level
    integer, intent(in) :: owners(:),maximum_parcels
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), intent(in) :: lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    logical, intent(out) :: ok
    real(dp), allocatable, optional, intent(inout) :: component_mass(:,:)
    type(spray_parcel), allocatable :: candidate(:)
    real(dp), allocatable :: masses(:,:),positions(:,:)
    integer, allocatable :: counts(:),offsets(:),locations(:,:)
    real(dp) :: bounds(6),reference_bounds(6)
    integer :: rank,nranks,ierr,nodes,low,high,node,p,ios,descriptor(13),reference(13),flags(5),reference_flags(5),ceiling
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
    flags=[merge(1,0,periodic),merge(1,0,present(component_mass)),ceiling];reference_flags=flags
    call MPI_Bcast(reference_flags,5,MPI_INTEGER,0,comm,ierr)
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
