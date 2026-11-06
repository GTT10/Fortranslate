module spray_amr_tree_partition_mod
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_integer_layout_3d
  use spray_parcel_mod, only: spray_parcel,valid_spray_parcel
  use spray_mixture_mod, only: valid_component_inventory
  implicit none
  private
  public :: locate_amr_tree_spray_cells,partition_amr_tree_spray_parcels
  public :: locate_amr_tree_ghost_cells
contains
  ! Half-open cells/patches; finest covering node wins. Periodic lookup wraps
  ! coordinates without modifying parcel positions. Outputs publish together.
  ! maximum_level (zero based) selects the finest covering node at or below
  ! that level, including covered parents, for coarse particle representations.
  subroutine locate_amr_tree_spray_cells(tree,positions,lower,spacing,periodic,locations,ok,maximum_level)
    type(amr_patch_tree_3d), intent(in) :: tree
    real(dp), intent(in) :: positions(:,:),lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    integer, optional, intent(in) :: maximum_level
    integer, allocatable, intent(inout) :: locations(:,:)
    logical, intent(out) :: ok
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:),candidate(:,:)
    integer(int64) :: cell(3)
    real(dp) :: length(3),fraction(3),coordinate(3)
    integer :: p,node,level,best,axis,ios,ceiling
    logical :: valid
    ok=.false.
    if (size(positions,1)/=3) return
    if (.not.all(ieee_is_finite(positions)).or..not.all(ieee_is_finite([lower,spacing]))) return
    if (any(abs(lower)>1.e9_dp).or.any(abs(positions)>1.e9_dp)) return
    if (any(spacing<1.e-12_dp).or.any(spacing>1.e9_dp)) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    if (.not.valid) return
    ceiling=tree%level_count()-1
    if(present(maximum_level)) ceiling=maximum_level
    if(ceiling<0.or.ceiling>=tree%level_count()) return
    ! Keep integer cell boundaries representable in the floating lookup.
    if (any(domain>2_int64**52)) return
    length=spacing*real(tree%nodes(1)%shape,dp)
    allocate(candidate(4,size(positions,2)),stat=ios)
    if (ios/=0) return
    do p=1,size(positions,2)
      coordinate=positions(:,p)-lower
      do axis=1,3
        if (periodic(axis)) then
          coordinate(axis)=modulo(coordinate(axis),length(axis))
        else
          if (coordinate(axis)<0._dp.or.coordinate(axis)>=length(axis)) return
        end if
      end do
      fraction=coordinate/length;best=-1
      do node=1,size(tree%nodes)
        if(tree%nodes(node)%level>ceiling) cycle
        level=tree%nodes(node)%level+1
        if (level<=best) cycle
        cell=min(domain(:,level)-1,floor(fraction*real(domain(:,level),dp),kind=int64))
        if (any(cell<origin(:,node)).or.any(cell-origin(:,node)>=int(tree%nodes(node)%shape,int64))) cycle
        best=level;candidate(:,p)=[node,int(cell-origin(:,node))+1]
      end do
      if (best<1) return
    end do
    call move_alloc(candidate,locations);ok=.true.
  end subroutine

  ! Fine-level boxes grown in fine cells. Node zero means no matching box.
  ! Valid cells take precedence over grown boxes; ties use tree order.
  ! Local cell indices may lie outside the valid patch. Coordinates are unchanged.
  subroutine locate_amr_tree_ghost_cells(tree,positions,lower,spacing,periodic,target_level,grow,locations,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    real(dp), intent(in) :: positions(:,:),lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    integer, intent(in) :: target_level,grow
    integer, allocatable, intent(inout) :: locations(:,:)
    logical, intent(out) :: ok
    integer, allocatable :: candidate(:,:),physical(:,:),ratios(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer(int64) :: cell(3),shifted(3),lo(3),hi(3),width
    real(dp) :: fraction(3),length(3)
    integer :: p,node,axis,pass
    logical :: valid
    ok=.false.
    if(target_level<1.or.grow<0.or.grow>1000000) return
    call locate_amr_tree_spray_cells(tree,positions,lower,spacing,periodic,physical,valid)
    if(.not.valid) return
    if(target_level>=tree%level_count()) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    if(.not.valid) return
    length=spacing*real(tree%nodes(1)%shape,dp)
    allocate(candidate(4,size(positions,2)));candidate=0
    do p=1,size(positions,2)
      fraction=(positions(:,p)-lower)/length
      do axis=1,3
        if(periodic(axis)) fraction(axis)=modulo(fraction(axis),1._dp)
      end do
      cell=min(domain(:,target_level+1)-1,floor(fraction*real(domain(:,target_level+1),dp),kind=int64))
      do pass=0,1
        width=int(grow,int64)*pass
        do node=1,size(tree%nodes)
          if(tree%nodes(node)%level/=target_level) cycle
          lo=origin(:,node)-width;hi=origin(:,node)+int(tree%nodes(node)%shape,int64)-1+width
          shifted=cell
          do axis=1,3
            if(.not.periodic(axis)) cycle
            if(shifted(axis)<lo(axis)) shifted(axis)=shifted(axis)+domain(axis,target_level+1)
            if(shifted(axis)>hi(axis)) shifted(axis)=shifted(axis)-domain(axis,target_level+1)
          end do
          if(any(shifted<lo).or.any(shifted>hi)) cycle
          candidate(:,p)=[node,int(shifted-origin(:,node))+1];exit
        end do
        if(candidate(1,p)>0) exit
      end do
    end do
    call move_alloc(candidate,locations);ok=.true.
  end subroutine locate_amr_tree_ghost_cells

  ! Input IDs are strictly increasing, as after MPI receive sorting. Stable
  ! rank bins preserve every parcel field and optional component association.
  subroutine partition_amr_tree_spray_parcels(tree,owners,nranks,parcels,lower,spacing,periodic, &
      counts,offsets,locations,ok,component_mass,maximum_level,tab_tangent_history)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, optional, intent(in) :: maximum_level
    integer, intent(in) :: owners(:),nranks
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), intent(in) :: lower(3),spacing(3)
    logical, intent(in) :: periodic(3)
    integer, allocatable, intent(inout) :: counts(:),offsets(:),locations(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable, optional, intent(inout) :: component_mass(:,:)
    real(dp), allocatable, optional, intent(inout) :: tab_tangent_history(:)
    real(dp), allocatable :: history(:)
    type(spray_parcel), allocatable :: ordered(:)
    real(dp), allocatable :: positions(:,:),masses(:,:)
    integer, allocatable :: found(:,:),bins(:),starts(:),next(:),ordered_locations(:,:)
    integer :: p,r,destination,ios
    logical :: valid
    ok=.false.
    if (.not.allocated(parcels).or..not.tree%is_valid()) return
    if (nranks<1.or.nranks>1000000.or.size(owners)/=size(tree%nodes)) return
    if (any(owners<0).or.any(owners>=nranks)) return
    if (present(component_mass)) then
      if (.not.allocated(component_mass)) return
      if (size(component_mass,2)/=size(parcels)) return
      if (size(component_mass,1)<1.or.size(component_mass,1)>16) return
    end if
    allocate(positions(3,size(parcels)),stat=ios)
    if (ios/=0) return
    if(present(tab_tangent_history)) then
      if(.not.allocated(tab_tangent_history)) return
      if(size(tab_tangent_history)/=size(parcels)) return
      if(.not.all(ieee_is_finite(tab_tangent_history)).or.any(tab_tangent_history<0)) return
      allocate(history(size(parcels)))
    endif
    do p=1,size(parcels)
      if (.not.valid_spray_parcel(parcels(p))) return
      if (p>1) then
        if (parcels(p)%id<=parcels(p-1)%id) return
      end if
      if (present(component_mass)) then
        if (.not.valid_component_inventory(parcels(p),component_mass(:,p))) return
      end if
      positions(:,p)=parcels(p)%position
    end do
    call locate_amr_tree_spray_cells(tree,positions,lower,spacing,periodic,found,valid,maximum_level)
    if (.not.valid) return
    allocate(bins(nranks),starts(nranks),next(nranks),ordered(size(parcels)), &
      ordered_locations(4,size(parcels)),stat=ios)
    if (ios/=0) return
    bins=0
    do p=1,size(parcels)
      r=owners(found(1,p))+1;bins(r)=bins(r)+1
    end do
    starts(1)=0
    do r=2,nranks
      starts(r)=starts(r-1)+bins(r-1)
    end do
    next=starts
    if (present(component_mass)) then
      allocate(masses(size(component_mass,1),size(parcels)),stat=ios)
      if (ios/=0) return
    end if
    do p=1,size(parcels)
      r=owners(found(1,p))+1;destination=next(r)+1;next(r)=destination
      ordered(destination)=parcels(p);ordered_locations(:,destination)=found(:,p)
      if (present(component_mass)) masses(:,destination)=component_mass(:,p)
      if(present(tab_tangent_history)) history(destination)=tab_tangent_history(p)
    end do
    call move_alloc(ordered,parcels);call move_alloc(bins,counts);call move_alloc(starts,offsets)
    call move_alloc(ordered_locations,locations)
    if (present(component_mass)) call move_alloc(masses,component_mass)
    if(present(tab_tangent_history)) call move_alloc(history,tab_tangent_history)
    ok=.true.
  end subroutine
end module
