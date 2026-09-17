module spray_mpi_amr_tree_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d
  use spray_parcel_mod, only: spray_parcel
  use spray_amr_tree_partition_mod, only: partition_amr_tree_spray_parcels,locate_amr_tree_spray_cells
  use spray_mpi_exchange_mod, only: exchange_spray_parcel_bins
  implicit none
  private
  public :: migrate_amr_tree_spray_parcels
contains
  subroutine migrate_amr_tree_spray_parcels(comm,tree,owners,parcels,lower,spacing,periodic, &
      maximum_parcels,ok,component_mass)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
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
    integer :: rank,nranks,ierr,nodes,low,high,node,p,ios,descriptor(13),reference(13),flags(4),reference_flags(4)
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
    nodes=size(tree%nodes)
    call MPI_Allreduce(nodes,low,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(nodes,high,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.low/=high) return
    bounds=[lower,spacing];reference_bounds=bounds
    call MPI_Bcast(reference_bounds,6,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    flags=[merge(1,0,periodic),merge(1,0,present(component_mass))];reference_flags=flags
    call MPI_Bcast(reference_flags,4,MPI_INTEGER,0,comm,ierr)
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
        counts,offsets,locations,valid,masses)
    else
      call partition_amr_tree_spray_parcels(tree,owners,nranks,candidate,lower,spacing,periodic, &
        counts,offsets,locations,valid)
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
    call locate_amr_tree_spray_cells(tree,positions,lower,spacing,periodic,locations,valid)
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
