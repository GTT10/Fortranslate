module amr_patch_tree_tagging_3d_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  implicit none
  private
  type, public :: hot_cell
    integer(int64) :: lower(3),upper(3)
    real(dp) :: temperature
  end type
  public :: tag_amr_patch_tree_temperature_3d,build_tagged_tree_3d
contains
  subroutine tag_amr_patch_tree_temperature_3d(old_tree,temperature,thresholds,ratios,tile_shape,buffer, &
      maximum_nodes,periodic,new_tree,ok)
    type(amr_patch_tree_3d), intent(in) :: old_tree
    type(amr_tree_state_3d), intent(in) :: temperature
    real(dp), intent(in) :: thresholds(:)
    integer, intent(in) :: ratios(:),tile_shape(3),buffer,maximum_nodes
    logical, intent(in) :: periodic(3)
    type(amr_patch_tree_3d), intent(inout) :: new_tree
    logical, intent(out) :: ok
    type(hot_cell), allocatable :: cells(:)
    integer(int64), allocatable :: old_origin(:,:),old_domain(:,:)
    integer, allocatable :: old_ratios(:)
    integer(int64) :: common(3),target(3),factor(3),a,b,g,total
    integer :: axis,level,node,n(3),i,j,k,idx(3),used,pass
    logical :: valid
    ok=.false.
    if (.not.temperature%is_valid(old_tree)) return
    if (size(temperature%fields(1)%values,1)/=1) return
    if (size(thresholds)/=size(ratios)) return
    if (any(.not.ieee_is_finite(thresholds))) return
    if (any(thresholds<=0)) return
    if (any(ratios<2).or.any(tile_shape<1).or.buffer<0.or.maximum_nodes<1) return
    do level=2,size(thresholds)
      if (thresholds(level)<thresholds(level-1)) return
    end do
    do node=1,size(old_tree%nodes)
      if (any(temperature%fields(node)%values<=0)) return
    end do
    call amr_tree_integer_layout_3d(old_tree,old_origin,old_domain,old_ratios,valid)
    if (.not.valid) return
    ! One exact integer lattice for both old cells and all requested new depths.
    ! Reserve headroom for buffered and periodically translated tile bounds.
    target=int(old_tree%nodes(1)%shape,int64)
    do level=1,size(ratios)
      if (any(target>(huge(1_int64)/4)/int(ratios(level),int64))) return
      target=target*int(ratios(level),int64)
    end do
    do axis=1,3
      a=old_domain(axis,size(old_ratios));b=target(axis);g=gcd(a,b)
      if (a/g>(huge(1_int64)/4)/b) return
      common(axis)=(a/g)*b
    end do
    ! Only exposed old cells vote. Covered parent temperatures may be stale;
    ! averaging a hot fine cell down before tagging would lose that feature.
    total=0
    do pass=1,2
      used=0
      do node=1,size(old_tree%nodes)
        n=old_tree%nodes(node)%shape;level=old_tree%nodes(node)%level+1
        factor=common/old_domain(:,level)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              if (size(thresholds)==0) cycle
              if (temperature%fields(node)%values(1,i,j,k)<thresholds(1)) cycle
              idx=[i,j,k]
              if (covered(old_tree,node,idx)) cycle
              if (pass==1) then
                total=total+1
                if (total>int(huge(1),int64)) return
              else
                used=used+1
                cells(used)%lower=(old_origin(:,node)+int(idx-1,int64))*factor
                cells(used)%upper=cells(used)%lower+factor
                cells(used)%temperature=temperature%fields(node)%values(1,i,j,k)
              end if
            end do
          end do
        end do
      end do
      if (pass==1) then
        allocate(cells(int(total)),stat=i)
        if (i/=0) return
      end if
    end do
    call build_tagged_tree_3d(old_tree,cells,common,thresholds,ratios,tile_shape,buffer, &
      maximum_nodes,periodic,new_tree,ok)
  end subroutine

  ! Internal shared builder: callers validate controls and supply hot exposed-cell boxes
  ! on the checked common lattice. No field payload is needed here.
  subroutine build_tagged_tree_3d(old_tree,cells,common,thresholds,ratios,tile_shape,buffer, &
      maximum_nodes,periodic,new_tree,ok)
    type(amr_patch_tree_3d), intent(in) :: old_tree
    type(hot_cell), intent(in) :: cells(:)
    integer(int64), intent(in) :: common(3)
    real(dp), intent(in) :: thresholds(:)
    integer, intent(in) :: ratios(:),tile_shape(3),buffer,maximum_nodes
    logical, intent(in) :: periodic(3)
    type(amr_patch_tree_3d), intent(inout) :: new_tree
    logical, intent(out) :: ok
    type(amr_patch_tree_3d) :: candidate
    integer(int64), allocatable :: origin(:,:),width(:,:)
    integer(int64) :: pad(3),lo(3),hi(3)
    integer :: used,parent,level,n(3),ratio,axis,i,j,k,lower(3),upper(3),child,ios
    logical :: tagged,valid
    ok=.false.
    allocate(candidate%nodes(maximum_nodes),origin(3,maximum_nodes),width(3,maximum_nodes),stat=ios)
    if (ios/=0) return
    candidate%nodes(1)=old_tree%nodes(1)
    origin(:,1)=0;width(:,1)=common/int(candidate%nodes(1)%shape,int64)
    used=1;parent=1
    ! Breadth-first ordering: deterministic k/j/i tiles, independent of old
    ! patch order or decomposition. Tiles partition each parent without overlap.
    do while(parent<=used)
      level=candidate%nodes(parent)%level+1
      if (level>size(ratios)) then
        parent=parent+1;cycle
      end if
      n=candidate%nodes(parent)%shape;ratio=ratios(level)
      if (any(n>huge(1)-tile_shape)) return
      ! A buffer wider than the domain is equivalent to the whole domain.
      do axis=1,3
        if (int(buffer,int64)>common(axis)/width(axis,parent)) then
          pad(axis)=common(axis)
        else
          pad(axis)=int(buffer,int64)*width(axis,parent)
        end if
      end do
      do k=1,n(3),tile_shape(3)
        do j=1,n(2),tile_shape(2)
          do i=1,n(1),tile_shape(1)
            lower=[i,j,k];upper=min(n,lower+tile_shape-1)
            lo=origin(:,parent)+int(lower-1,int64)*width(:,parent)-pad
            hi=origin(:,parent)+int(upper,int64)*width(:,parent)+pad
            tagged=.false.
            do child=1,size(cells)
              if (cells(child)%temperature<thresholds(level)) cycle
              if (intersects(lo,hi,cells(child)%lower,cells(child)%upper,common,periodic)) then
                tagged=.true.;exit
              end if
            end do
            if (.not.tagged) cycle
            if (used==maximum_nodes) return
            if (any(int(upper-lower+1,int64)*int(ratio,int64)>int(huge(1),int64))) return
            child=used+1
            call initialize_amr_patch_3d(n(1),n(2),n(3),lower(1),upper(1),lower(2),upper(2), &
              lower(3),upper(3),ratio,candidate%nodes(child)%patch,valid)
            if (.not.valid) return
            candidate%nodes(child)%parent=parent;candidate%nodes(child)%level=level
            candidate%nodes(child)%shape=(upper-lower+1)*ratio
            origin(:,child)=origin(:,parent)+int(lower-1,int64)*width(:,parent)
            width(:,child)=width(:,parent)/int(ratio,int64)
            used=child
          end do
        end do
      end do
      parent=parent+1
    end do
    candidate%nodes=candidate%nodes(:used)
    if (.not.candidate%is_valid()) return
    new_tree=candidate;ok=.true.
  end subroutine

  pure logical function intersects(lo,hi,source_lo,source_hi,domain,periodic) result(hit)
    integer(int64), intent(in) :: lo(3),hi(3),source_lo(3),source_hi(3),domain(3)
    logical, intent(in) :: periodic(3)
    integer :: axis
    logical :: overlap
    hit=.false.
    do axis=1,3
      overlap=lo(axis)<source_hi(axis).and.hi(axis)>source_lo(axis)
      if (periodic(axis)) then
        overlap=overlap.or.(lo(axis)<source_hi(axis)-domain(axis).and.hi(axis)>source_lo(axis)-domain(axis))
        overlap=overlap.or.(lo(axis)<source_hi(axis)+domain(axis).and.hi(axis)>source_lo(axis)+domain(axis))
      end if
      if (.not.overlap) return
    end do
    hit=.true.
  end function

  pure integer(int64) function gcd(a,b) result(g)
    integer(int64), intent(in) :: a,b
    integer(int64) :: x,y,r
    x=a;y=b
    do while(y/=0)
      r=modulo(x,y);x=y;y=r
    end do
    g=x
  end function

  pure logical function covered(tree,node,idx)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: node,idx(3)
    integer :: child,lo(3),hi(3)
    covered=.false.
    do child=node+1,size(tree%nodes)
      if (tree%nodes(child)%parent/=node) cycle
      associate(box=>tree%nodes(child)%patch)
        lo=[box%coarse_i_lower,box%coarse_j_lower,box%coarse_k_lower]
        hi=[box%coarse_i_upper,box%coarse_j_upper,box%coarse_k_upper]
      end associate
      if (all(idx>=lo).and.all(idx<=hi)) then
        covered=.true.;return
      end if
    end do
  end function
end module amr_patch_tree_tagging_3d_mod
