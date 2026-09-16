module amr_patch_tree_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,prolong_pcm_3d,average_down_3d
  implicit none
  private

  type, public :: amr_tree_node_3d
    integer :: parent=0,level=0
    integer :: shape(3)=0
    type(amr_patch_3d) :: patch
  end type
  type, public :: amr_patch_tree_3d
    type(amr_tree_node_3d), allocatable :: nodes(:)
  contains
    procedure :: is_valid => tree_valid
    procedure :: level_count => tree_levels
  end type
  type, public :: amr_tree_field_3d
    real(dp), allocatable :: values(:,:,:,:)
  end type
  type, public :: amr_tree_state_3d
    type(amr_tree_field_3d), allocatable :: fields(:)
  contains
    procedure :: is_valid => state_valid
  end type
  public :: amr_tree_integer_layout_3d
  public :: initialize_amr_patch_tree_3d,prolong_amr_patch_tree_3d
  public :: average_down_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
contains
  subroutine amr_tree_integer_layout_3d(tree,origin,domain,ratios,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer(int64), allocatable, intent(out) :: origin(:,:),domain(:,:)
    integer, allocatable, intent(out) :: ratios(:)
    logical, intent(out) :: ok
    integer :: node,parent,level
    ok=.false.
    if (.not.tree%is_valid()) return
    allocate(origin(3,size(tree%nodes)),domain(3,tree%level_count()),ratios(tree%level_count()))
    origin=0;ratios=1;domain(:,1)=int(tree%nodes(1)%shape,int64)
    do node=2,size(tree%nodes)
      ratios(tree%nodes(node)%level+1)=tree%nodes(node)%patch%refinement_ratio
    end do
    do level=2,size(ratios)
      if (any(domain(:,level-1)>(huge(1_int64)-2)/int(ratios(level),int64))) return
      domain(:,level)=domain(:,level-1)*int(ratios(level),int64)
    end do
    do node=2,size(tree%nodes)
      parent=tree%nodes(node)%parent
      origin(:,node)=(origin(:,parent)+int(lower(tree%nodes(node)%patch)-1,int64))* &
        int(tree%nodes(node)%patch%refinement_ratio,int64)
    end do
    ok=.true.
  end subroutine amr_tree_integer_layout_3d

  pure function lower(patch) result(bounds)
    type(amr_patch_3d), intent(in) :: patch
    integer :: bounds(3)
    bounds=[patch%coarse_i_lower,patch%coarse_j_lower,patch%coarse_k_lower]
  end function
  pure function upper(patch) result(bounds)
    type(amr_patch_3d), intent(in) :: patch
    integer :: bounds(3)
    bounds=[patch%coarse_i_upper,patch%coarse_j_upper,patch%coarse_k_upper]
  end function
  pure logical function tree_valid(self) result(valid)
    class(amr_patch_tree_3d), intent(in) :: self
    integer :: node,parent,other
    integer(int64) :: extent(3)
    valid=.false.
    if (.not.allocated(self%nodes)) return
    if (size(self%nodes)<1) return
    if (self%nodes(1)%parent/=0.or.self%nodes(1)%level/=0.or.any(self%nodes(1)%shape<2)) return
    do node=2,size(self%nodes)
      parent=self%nodes(node)%parent
      ! Parent-before-child ordering permits a single reverse traversal for restriction.
      if (parent<1.or.parent>=node) return
      if (self%nodes(node)%level/=self%nodes(parent)%level+1) return
      if (.not.self%nodes(node)%patch%is_valid()) return
      if (any([self%nodes(node)%patch%coarse_nx,self%nodes(node)%patch%coarse_ny, &
          self%nodes(node)%patch%coarse_nz]/=self%nodes(parent)%shape)) return
      extent=int(upper(self%nodes(node)%patch),int64)-int(lower(self%nodes(node)%patch),int64)+1
      extent=extent*int(self%nodes(node)%patch%refinement_ratio,int64)
      if (any(extent>int(huge(1),int64))) return
      if (any(int(extent)/=self%nodes(node)%shape)) return
      do other=2,node-1
        if (self%nodes(other)%level==self%nodes(node)%level) then
          if (self%nodes(other)%patch%refinement_ratio/=self%nodes(node)%patch%refinement_ratio) return
        end if
        if (self%nodes(other)%parent/=parent) cycle
        ! Face/edge/corner adjacency is allowed; sharing a parent cell is not.
        if (all(lower(self%nodes(node)%patch)<=upper(self%nodes(other)%patch)).and. &
            all(lower(self%nodes(other)%patch)<=upper(self%nodes(node)%patch))) return
      end do
    end do
    valid=.true.
  end function
  pure integer function tree_levels(self) result(count)
    class(amr_patch_tree_3d), intent(in) :: self
    count=0
    if (self%is_valid()) count=maxval(self%nodes%level)+1
  end function
  subroutine initialize_amr_patch_tree_3d(root_shape,parents,patches,tree,ok)
    integer, intent(in) :: root_shape(3),parents(:)
    type(amr_patch_3d), intent(in) :: patches(:)
    type(amr_patch_tree_3d), intent(inout) :: tree
    logical, intent(out) :: ok
    type(amr_patch_tree_3d) :: candidate
    integer :: node,parent
    integer(int64) :: extent(3)
    ok=.false.
    if (size(parents)/=size(patches).or.any(root_shape<2)) return
    allocate(candidate%nodes(size(parents)+1))
    candidate%nodes(1)%shape=root_shape
    do node=2,size(candidate%nodes)
      parent=parents(node-1)
      if (parent<1.or.parent>=node) return
      if (.not.patches(node-1)%is_valid()) return
      extent=int(upper(patches(node-1)),int64)-int(lower(patches(node-1)),int64)+1
      extent=extent*int(patches(node-1)%refinement_ratio,int64)
      if (any(extent>int(huge(1),int64))) return
      candidate%nodes(node)%parent=parent
      candidate%nodes(node)%level=candidate%nodes(parent)%level+1
      candidate%nodes(node)%shape=int(extent)
      candidate%nodes(node)%patch=patches(node-1)
    end do
    if (.not.candidate%is_valid()) return
    tree=candidate;ok=.true.
  end subroutine
  pure logical function state_valid(self,tree) result(valid)
    class(amr_tree_state_3d), intent(in) :: self
    type(amr_patch_tree_3d), intent(in) :: tree
    integer :: node,nvar
    valid=.false.
    if (.not.tree%is_valid().or..not.allocated(self%fields)) return
    if (size(self%fields)/=size(tree%nodes)) return
    if (.not.allocated(self%fields(1)%values)) return
    nvar=size(self%fields(1)%values,1)
    if (nvar<1) return
    do node=1,size(tree%nodes)
      if (.not.allocated(self%fields(node)%values)) return
      if (any(shape(self%fields(node)%values)/=[nvar,tree%nodes(node)%shape])) return
      if (.not.all(ieee_is_finite(self%fields(node)%values))) return
    end do
    valid=.true.
  end function
  subroutine prolong_amr_patch_tree_3d(tree,root,state,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    real(dp), intent(in) :: root(:,:,:,:)
    type(amr_tree_state_3d), intent(inout) :: state
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: candidate
    integer :: node,parent,nvar
    ok=.false.
    if (.not.tree%is_valid()) return
    nvar=size(root,1)
    if (nvar<1.or.any(shape(root)/=[nvar,tree%nodes(1)%shape])) return
    if (.not.all(ieee_is_finite(root))) return
    allocate(candidate%fields(size(tree%nodes)))
    candidate%fields(1)%values=root
    do node=2,size(tree%nodes)
      parent=tree%nodes(node)%parent
      allocate(candidate%fields(node)%values(nvar,tree%nodes(node)%shape(1), &
        tree%nodes(node)%shape(2),tree%nodes(node)%shape(3)))
      call prolong_pcm_3d(candidate%fields(parent)%values,tree%nodes(node)%patch,candidate%fields(node)%values,ok)
      if (.not.ok) return
    end do
    state=candidate;ok=.true.
  end subroutine
  subroutine average_down_amr_patch_tree_3d(tree,state,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(inout) :: state
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: candidate
    integer :: node,parent
    ok=.false.
    if (.not.state%is_valid(tree)) return
    candidate=state
    do node=size(tree%nodes),2,-1
      parent=tree%nodes(node)%parent
      call average_down_3d(candidate%fields(parent)%values,candidate%fields(node)%values,tree%nodes(node)%patch,ok)
      if (.not.ok) return
    end do
    state=candidate;ok=.true.
  end subroutine
  subroutine composite_integrals_amr_patch_tree_3d(tree,state,root_spacing,integrals,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: state
    real(dp), intent(in) :: root_spacing(3)
    real(dp), intent(inout) :: integrals(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: candidate(:),volume(:),correction(:),term(:),updated(:)
    logical, allocatable :: covered(:,:,:)
    integer :: node,child,parent,i,j,k,lo(3),hi(3)
    real(dp) :: ratio
    ok=.false.
    if (.not.state%is_valid(tree)) return
    if (size(integrals)/=size(state%fields(1)%values,1)) return
    if (.not.all(ieee_is_finite(root_spacing))) return
    if (any(root_spacing<=0)) return
    allocate(candidate(size(integrals)),volume(size(tree%nodes)),correction(size(integrals)), &
      term(size(integrals)),updated(size(integrals)))
    candidate=0;correction=0
    volume(1)=product(root_spacing)
    if (.not.ieee_is_finite(volume(1)).or.volume(1)<=0) return
    do node=1,size(tree%nodes)
      if (node>1) then
        parent=tree%nodes(node)%parent;ratio=real(tree%nodes(node)%patch%refinement_ratio,dp)
        volume(node)=volume(parent)/ratio/ratio/ratio
        if (volume(node)<=0) return
      end if
      allocate(covered(tree%nodes(node)%shape(1),tree%nodes(node)%shape(2),tree%nodes(node)%shape(3)))
      covered=.false.
      do child=node+1,size(tree%nodes)
        if (tree%nodes(child)%parent/=node) cycle
        lo=lower(tree%nodes(child)%patch);hi=upper(tree%nodes(child)%patch)
        covered(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.true.
      end do
      do k=1,size(covered,3)
        do j=1,size(covered,2)
          do i=1,size(covered,1)
            if (covered(i,j,k)) cycle
            ! Compensate summation across levels with different cell volumes.
            term=volume(node)*state%fields(node)%values(:,i,j,k)-correction
            updated=candidate+term
            correction=(updated-candidate)-term
            candidate=updated
          end do
        end do
      end do
      deallocate(covered)
    end do
    if (.not.all(ieee_is_finite(candidate))) return
    integrals=candidate;ok=.true.
  end subroutine
end module amr_patch_tree_3d_mod
