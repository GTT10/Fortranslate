module amr_patch_tree_reflux_3d_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: average_down_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d, &
    amr_tree_integer_layout_3d,average_down_amr_patch_tree_3d
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d
  implicit none
  private
  type, public :: amr_tree_flux_3d
    ! Time-integrated flux per unit face area; normal indices run from 1 to n+1.
    real(dp), allocatable :: x(:,:,:,:),y(:,:,:,:),z(:,:,:,:)
  end type
  type :: face_marks
    logical, allocatable :: x(:,:,:),y(:,:,:),z(:,:,:)
  end type
  public :: reflux_amr_patch_tree_3d
contains
  subroutine reflux_amr_patch_tree_3d(tree,state,flux,root_spacing,faces,boundary_correction,maximum_correction,ok, &
      minimum_depth,synchronized_flux)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(inout) :: state
    type(amr_tree_flux_3d), intent(in) :: flux(:)
    real(dp), intent(in) :: root_spacing(3)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(inout) :: boundary_correction(:),maximum_correction
    logical, intent(out) :: ok
    integer, optional, intent(in) :: minimum_depth
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: synchronized_flux(:)
    type(amr_tree_state_3d) :: candidate
    type(amr_tree_flux_3d), allocatable :: work(:)
    type(face_marks), allocatable :: refined(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: spacing(:,:),budget(:),correction(:),delta_low(:),delta_high(:)
    logical, allocatable :: covered(:,:,:)
    real(dp) :: largest,area
    integer :: node,child,parent,level,nvar,n(3),lo(3),hi(3),i,j,k,axis,idx(3),upper_idx(3),floor_depth
    logical :: valid
    ok=.false.
    if (.not.state%is_valid(tree)) return
    floor_depth=0
    if (present(minimum_depth)) floor_depth=minimum_depth
    if (floor_depth<0.or.floor_depth>=tree%level_count()) return
    if (size(flux)/=size(tree%nodes).or..not.valid_reactive_boundaries_3d(faces)) return
    nvar=size(state%fields(1)%values,1)
    if (size(boundary_correction)/=nvar) return
    if (.not.all(ieee_is_finite(root_spacing))) return
    if (any(root_spacing<=0)) return
    do node=1,size(flux)
      if (tree%nodes(node)%level<floor_depth) cycle
      n=tree%nodes(node)%shape
      if (any(n==huge(1))) return
      if (.not.allocated(flux(node)%x).or..not.allocated(flux(node)%y).or..not.allocated(flux(node)%z)) return
      if (any(shape(flux(node)%x)/=[nvar,n(1)+1,n(2),n(3)]).or. &
          any(shape(flux(node)%y)/=[nvar,n(1),n(2)+1,n(3)]).or. &
          any(shape(flux(node)%z)/=[nvar,n(1),n(2),n(3)+1])) return
      if (.not.all(ieee_is_finite(flux(node)%x)).or..not.all(ieee_is_finite(flux(node)%y)).or. &
          .not.all(ieee_is_finite(flux(node)%z))) return
    end do
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    if (.not.valid) return
    allocate(spacing(3,size(ratios)),budget(nvar),correction(nvar),delta_low(nvar),delta_high(nvar))
    spacing(:,1)=root_spacing
    do level=2,size(ratios)
      spacing(:,level)=spacing(:,level-1)/real(ratios(level),dp)
      if (any(spacing(:,level)<=0)) return
    end do
    candidate=state;work=flux;budget=0;largest=0
    allocate(refined(size(flux)))
    do node=1,size(flux)
      if (tree%nodes(node)%level<floor_depth) cycle
      n=tree%nodes(node)%shape
      allocate(refined(node)%x(n(1)+1,n(2),n(3)),refined(node)%y(n(1),n(2)+1,n(3)), &
        refined(node)%z(n(1),n(2),n(3)+1))
      refined(node)%x=.false.;refined(node)%y=.false.;refined(node)%z=.false.
    end do
    do level=size(ratios),floor_depth+1,-1
      call synchronize_peers(level,valid)
      if (.not.valid) return
      do node=1,size(flux)
        if (tree%nodes(node)%level+1/=level) cycle
        n=tree%nodes(node)%shape
        allocate(covered(n(1),n(2),n(3)));covered=.false.
        do child=node+1,size(flux)
          if (tree%nodes(child)%parent/=node) cycle
          lo=[tree%nodes(child)%patch%coarse_i_lower,tree%nodes(child)%patch%coarse_j_lower, &
            tree%nodes(child)%patch%coarse_k_lower]
          hi=[tree%nodes(child)%patch%coarse_i_upper,tree%nodes(child)%patch%coarse_j_upper, &
            tree%nodes(child)%patch%coarse_k_upper]
          covered(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.true.
        end do
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              if (covered(i,j,k)) cycle
              idx=[i,j,k];correction=0
              do axis=1,3
                upper_idx=idx;upper_idx(axis)=upper_idx(axis)+1
                delta_low=get_face(work(node),axis,idx)-get_face(flux(node),axis,idx)
                delta_high=get_face(work(node),axis,upper_idx)-get_face(flux(node),axis,upper_idx)
                correction=correction+(delta_low-delta_high)/spacing(axis,level)
              end do
              candidate%fields(node)%values(:,i,j,k)=candidate%fields(node)%values(:,i,j,k)+correction
              largest=max(largest,maxval(abs(correction)))
            end do
          end do
        end do
        deallocate(covered)
      end do
      ! Updated child fluxes include deeper corrections before being area-averaged.
      if (level>floor_depth+1) then
        do node=2,size(flux)
          if (tree%nodes(node)%level+1/=level) cycle
          parent=tree%nodes(node)%parent
          call restrict_faces(node,parent,valid)
          if (.not.valid) return
        end do
      end if
    end do
    n=tree%nodes(1)%shape
    do axis=1,3
      if (floor_depth>0) exit
      if (faces(2*axis-1)=='periodic') cycle
      area=product(root_spacing)/root_spacing(axis)
      hi=n;hi(axis)=1
      do k=1,hi(3)
        do j=1,hi(2)
          do i=1,hi(1)
            idx=[i,j,k];upper_idx=idx;upper_idx(axis)=n(axis)+1
            budget=budget+area*(get_face(work(1),axis,idx)-get_face(flux(1),axis,idx)- &
              get_face(work(1),axis,upper_idx)+get_face(flux(1),axis,upper_idx))
          end do
        end do
      end do
    end do
    if (.not.all(ieee_is_finite(budget)).or..not.ieee_is_finite(largest)) return
    if (floor_depth==0) then
      call average_down_amr_patch_tree_3d(tree,candidate,valid)
      if (.not.valid) return
    else
      ! Do not overwrite coarse endpoints that belong to a longer interval.
      do node=size(tree%nodes),2,-1
        parent=tree%nodes(node)%parent
        if (tree%nodes(parent)%level<floor_depth) cycle
        call average_down_3d(candidate%fields(parent)%values,candidate%fields(node)%values,tree%nodes(node)%patch,valid)
        if (.not.valid) return
      end do
    end if
    if (.not.candidate%is_valid(tree)) return
    if (present(synchronized_flux)) then
      do node=1,size(work)
        if (tree%nodes(node)%level>=floor_depth) cycle
        if (allocated(work(node)%x)) deallocate(work(node)%x)
        if (allocated(work(node)%y)) deallocate(work(node)%y)
        if (allocated(work(node)%z)) deallocate(work(node)%z)
      end do
      call move_alloc(work,synchronized_flux)
    end if
    state=candidate;boundary_correction=budget;maximum_correction=largest;ok=.true.
  contains
    subroutine restrict_faces(child_node,parent_node,success)
      integer, intent(in) :: child_node,parent_node
      logical, intent(out) :: success
      integer :: a,b,c,r,low(3),high(3),last(3),p(3),f(3),i1,j1,k1,s,t
      real(dp) :: average(nvar)
      success=.false.;r=tree%nodes(child_node)%patch%refinement_ratio
      low=[tree%nodes(child_node)%patch%coarse_i_lower,tree%nodes(child_node)%patch%coarse_j_lower, &
        tree%nodes(child_node)%patch%coarse_k_lower]
      high=[tree%nodes(child_node)%patch%coarse_i_upper,tree%nodes(child_node)%patch%coarse_j_upper, &
        tree%nodes(child_node)%patch%coarse_k_upper]
      do a=1,3
        b=modulo(a,3)+1;c=modulo(a+1,3)+1;last=high;last(a)=last(a)+1
        do k1=low(3),last(3)
          do j1=low(2),last(2)
            do i1=low(1),last(1)
              p=[i1,j1,k1];average=0
              do t=0,r-1
                do s=0,r-1
                  f=(p-low)*r+1;f(b)=f(b)+s;f(c)=f(c)+t
                  average=average+get_face(work(child_node),a,f)/real(r,dp)/real(r,dp)
                end do
              end do
              if (marked(refined(parent_node),a,p)) then
                if (.not.matches(average,get_face(work(parent_node),a,p))) return
              else
                call set_face(work(parent_node),a,p,average)
                call set_mark(refined(parent_node),a,p)
              end if
            end do
          end do
        end do
      end do
      success=.true.
    end subroutine
    subroutine synchronize_peers(depth,success)
      integer, intent(in) :: depth
      logical, intent(out) :: success
      integer :: left,right,a,side,other_side,p(3),q(3),last(3),i1,j1,k1
      integer(int64) :: key(3),normal
      real(dp) :: chosen(nvar),other_value(nvar)
      logical :: left_refined,right_refined
      success=.false.
      do left=1,size(flux)
        if (tree%nodes(left)%level+1/=depth) cycle
        do a=1,3
          last=tree%nodes(left)%shape;last(a)=1
          do side=1,2
            do k1=1,last(3)
              do j1=1,last(2)
                do i1=1,last(1)
                  p=[i1,j1,k1];p(a)=1+(side-1)*tree%nodes(left)%shape(a)
                  key=origin(:,left)+int(p-1,int64)
                  if (faces(2*a-1)=='periodic') key(a)=modulo(key(a),domain(a,depth))
                  do right=1,size(flux)
                    if (tree%nodes(right)%level+1/=depth) cycle
                    q=1
                    do other_side=1,3
                      if (other_side==a) cycle
                      if (key(other_side)<origin(other_side,right).or. &
                          key(other_side)>=origin(other_side,right)+tree%nodes(right)%shape(other_side)) exit
                      q(other_side)=int(key(other_side)-origin(other_side,right))+1
                    end do
                    if (other_side<=3) cycle
                    do other_side=1,2
                      q(a)=1+(other_side-1)*tree%nodes(right)%shape(a)
                      if (right==left.and.other_side==side) cycle
                      normal=origin(a,right)+int(q(a)-1,int64)
                      if (faces(2*a-1)=='periodic') normal=modulo(normal,domain(a,depth))
                      if (normal/=key(a)) cycle
                      chosen=get_face(work(left),a,p);other_value=get_face(work(right),a,q)
                      left_refined=marked(refined(left),a,p);right_refined=marked(refined(right),a,q)
                      if (left_refined.eqv.right_refined) then
                        if (.not.matches(chosen,other_value)) return
                      else if (right_refined) then
                        chosen=other_value
                      end if
                      call set_face(work(left),a,p,chosen);call set_face(work(right),a,q,chosen)
                      if (left_refined.or.right_refined) then
                        call set_mark(refined(left),a,p);call set_mark(refined(right),a,q)
                      end if
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
      success=.true.
    end subroutine
  end subroutine reflux_amr_patch_tree_3d
  function get_face(field,axis,index) result(value)
    type(amr_tree_flux_3d), intent(in) :: field
    integer, intent(in) :: axis,index(3)
    real(dp) :: value(size(field%x,1))
    select case(axis)
    case(1)
      value=field%x(:,index(1),index(2),index(3))
    case(2)
      value=field%y(:,index(1),index(2),index(3))
    case(3)
      value=field%z(:,index(1),index(2),index(3))
    end select
  end function
  subroutine set_face(field,axis,index,value)
    type(amr_tree_flux_3d), intent(inout) :: field
    integer, intent(in) :: axis,index(3)
    real(dp), intent(in) :: value(:)
    select case(axis)
    case(1)
      field%x(:,index(1),index(2),index(3))=value
    case(2)
      field%y(:,index(1),index(2),index(3))=value
    case(3)
      field%z(:,index(1),index(2),index(3))=value
    end select
  end subroutine
  logical function marked(field,axis,index) result(value)
    type(face_marks), intent(in) :: field
    integer, intent(in) :: axis,index(3)
    select case(axis)
    case(1)
      value=field%x(index(1),index(2),index(3))
    case(2)
      value=field%y(index(1),index(2),index(3))
    case(3)
      value=field%z(index(1),index(2),index(3))
    end select
  end function
  subroutine set_mark(field,axis,index)
    type(face_marks), intent(inout) :: field
    integer, intent(in) :: axis,index(3)
    select case(axis)
    case(1)
      field%x(index(1),index(2),index(3))=.true.
    case(2)
      field%y(index(1),index(2),index(3))=.true.
    case(3)
      field%z(index(1),index(2),index(3))=.true.
    end select
  end subroutine
  logical function matches(a,b) result(ok)
    real(dp), intent(in) :: a(:),b(:)
    ok=all(ieee_is_finite(a)).and.all(ieee_is_finite(b))
    if (ok) ok=maxval(abs(a-b))<=64*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a)),maxval(abs(b)))
  end function
end module amr_patch_tree_reflux_3d_mod
