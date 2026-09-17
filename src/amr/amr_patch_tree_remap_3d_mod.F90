module amr_patch_tree_remap_3d_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d, &
    average_down_amr_patch_tree_3d
  implicit none
  private
  type :: leaf_cell
    integer(int64) :: lower(3),upper(3)
    integer :: node,index(3)
  end type
  public :: remap_amr_patch_tree_3d,remap_amr_patch_tree_reactive_3d
  ! Shared geometry predicates for distributed conservative remapping.
  public :: covered,gcd
contains
  subroutine remap_amr_patch_tree_3d(old_tree,old_state,new_tree,new_state,ok)
    type(amr_patch_tree_3d), intent(in) :: old_tree,new_tree
    type(amr_tree_state_3d), intent(in) :: old_state
    type(amr_tree_state_3d), intent(inout) :: new_state
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: candidate
    type(leaf_cell), allocatable :: cells(:)
    integer(int64), allocatable :: old_origin(:,:),old_domain(:,:),new_origin(:,:),new_domain(:,:)
    integer, allocatable :: old_ratios(:),new_ratios(:)
    integer(int64) :: common(3),a,b,g,total,factor(3),lo(3),hi(3),overlap(3),count
    real(dp), allocatable :: value(:),compensation(:),term(:),increment(:),next_value(:)
    real(dp) :: weight,coverage,coverage_comp,delta,next_coverage
    integer :: node,n(3),i,j,k,axis,level,nvar,cell,used,idx(3),source
    logical :: valid
    ok=.false.
    if (.not.old_state%is_valid(old_tree).or..not.new_tree%is_valid()) return
    if (any(old_tree%nodes(1)%shape/=new_tree%nodes(1)%shape)) return
    call amr_tree_integer_layout_3d(old_tree,old_origin,old_domain,old_ratios,valid)
    if (.not.valid) return
    call amr_tree_integer_layout_3d(new_tree,new_origin,new_domain,new_ratios,valid)
    if (.not.valid) return
    ! Integer least-common-multiple coordinates represent every old/new face
    ! exactly, even when the refinement ratios change. Reject overflow.
    do axis=1,3
      a=old_domain(axis,size(old_ratios));b=new_domain(axis,size(new_ratios))
      g=gcd(a,b)
      if (a/g>huge(1_int64)/b) return
      common(axis)=(a/g)*b
    end do
    nvar=size(old_state%fields(1)%values,1);total=0
    do node=1,size(old_tree%nodes)
      count=size(old_state%fields(node)%values,kind=int64)/int(nvar,int64)
      if (count>int(huge(1),int64)-total) return
      total=total+count
    end do
    allocate(cells(int(total)));used=0
    do node=1,size(old_tree%nodes)
      n=old_tree%nodes(node)%shape;level=old_tree%nodes(node)%level+1
      factor=common/old_domain(:,level)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            idx=[i,j,k]
            if (covered(old_tree,node,idx)) cycle
            used=used+1;cells(used)%node=node;cells(used)%index=idx
            cells(used)%lower=(old_origin(:,node)+int(idx-1,int64))*factor
            cells(used)%upper=cells(used)%lower+factor
          end do
        end do
      end do
    end do
    allocate(candidate%fields(size(new_tree%nodes)))
    allocate(value(nvar),compensation(nvar),term(nvar),increment(nvar),next_value(nvar))
    do node=1,size(new_tree%nodes)
      n=new_tree%nodes(node)%shape;level=new_tree%nodes(node)%level+1
      factor=common/new_domain(:,level)
      allocate(candidate%fields(node)%values(nvar,n(1),n(2),n(3)));candidate%fields(node)%values=0
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            idx=[i,j,k]
            if (covered(new_tree,node,idx)) cycle
            lo=(new_origin(:,node)+int(idx-1,int64))*factor;hi=lo+factor
            value=0;compensation=0;coverage=0;coverage_comp=0
            do cell=1,used
              overlap=min(hi,cells(cell)%upper)-max(lo,cells(cell)%lower)
              if (any(overlap<=0)) cycle
              weight=product(real(overlap,dp)/real(factor,dp))
              source=cells(cell)%node;idx=cells(cell)%index
              term=weight*old_state%fields(source)%values(:,idx(1),idx(2),idx(3))
              increment=term-compensation;next_value=value+increment
              compensation=(next_value-value)-increment;value=next_value
              delta=weight-coverage_comp;next_coverage=coverage+delta
              coverage_comp=(next_coverage-coverage)-delta;coverage=next_coverage
            end do
            if (abs(coverage-1.0_dp)>128*epsilon(1.0_dp)) return
            candidate%fields(node)%values(:,i,j,k)=value/coverage
          end do
        end do
      end do
    end do
    call average_down_amr_patch_tree_3d(new_tree,candidate,valid)
    if (.not.valid) return
    if (.not.candidate%is_valid(new_tree)) return
    new_state=candidate;ok=.true.
  end subroutine

  subroutine remap_amr_patch_tree_reactive_3d(old_tree,old_state,old_temperature,new_tree,species,state,temperature,ok)
    type(amr_patch_tree_3d), intent(in) :: old_tree,new_tree
    type(amr_tree_state_3d), intent(in) :: old_state,old_temperature
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: candidate,temps
    real(dp), allocatable :: recovered(:,:,:)
    integer :: node,n(3)
    logical :: valid
    ok=.false.
    if (.not.old_state%is_valid(old_tree).or..not.old_temperature%is_valid(old_tree)) return
    if (size(old_state%fields(1)%values,1)/=reactive_nvar(size(species))) return
    if (size(old_temperature%fields(1)%values,1)/=1) return
    do node=1,size(old_tree%nodes)
      if (any(old_temperature%fields(node)%values<=0)) return
    end do
    call remap_amr_patch_tree_3d(old_tree,old_state,new_tree,candidate,valid)
    if (.not.valid) return
    ! Temperature is only a starting estimate. Recover thermodynamic state
    ! from the volume-averaged conserved mixture before publishing either.
    call remap_amr_patch_tree_3d(old_tree,old_temperature,new_tree,temps,valid)
    if (.not.valid) return
    do node=1,size(new_tree%nodes)
      n=new_tree%nodes(node)%shape;allocate(recovered(n(1),n(2),n(3)))
      call recover_reactive_temperatures_3d(species,candidate%fields(node)%values,temps%fields(node)%values(1,:,:,:), &
        n(1),n(2),n(3),recovered,valid)
      if (.not.valid) return
      temps%fields(node)%values(1,:,:,:)=recovered
      deallocate(recovered)
    end do
    state=candidate;temperature=temps;ok=.true.
  end subroutine

  pure integer(int64) function gcd(a,b) result(g)
    integer(int64), intent(in) :: a,b
    integer(int64) :: x,y,remainder
    x=a;y=b
    do while(y/=0)
      remainder=modulo(x,y);x=y;y=remainder
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
      lo=[tree%nodes(child)%patch%coarse_i_lower,tree%nodes(child)%patch%coarse_j_lower, &
        tree%nodes(child)%patch%coarse_k_lower]
      hi=[tree%nodes(child)%patch%coarse_i_upper,tree%nodes(child)%patch%coarse_j_upper, &
        tree%nodes(child)%patch%coarse_k_upper]
      if (all(idx>=lo).and.all(idx<=hi)) then
        covered=.true.;return
      end if
    end do
  end function
end module amr_patch_tree_remap_3d_mod
