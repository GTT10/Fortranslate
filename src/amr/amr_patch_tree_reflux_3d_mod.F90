module amr_patch_tree_reflux_3d_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: average_down_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d, &
    amr_tree_integer_layout_3d,average_down_amr_patch_tree_3d
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d
  use amr_patch_tree_remap_3d_mod, only: cell_covered => covered
  use amr_tree_reflux_graph_3d_mod, only: build_amr_reflux_graph_3d
  use conservative_reflux_budget_mod, only: spread_reflux_return
  implicit none
  private
  type, public :: amr_tree_flux_3d
    ! Time-integrated flux per unit face area; normal indices run from 1 to n+1.
    real(dp), allocatable :: x(:,:,:,:),y(:,:,:,:),z(:,:,:,:)
  end type
  type, public :: face_marks
    logical, allocatable :: x(:,:,:),y(:,:,:),z(:,:,:)
  end type
  public :: reflux_amr_patch_tree_3d
  public :: get_face,set_face,marked,set_mark,matches
contains
  subroutine reflux_amr_patch_tree_3d(tree,state,flux,root_spacing,faces,boundary_correction,maximum_correction,ok, &
      minimum_depth,synchronized_flux,first_species,last_species, &
      auxiliary_state,auxiliary_flux,auxiliary_synchronized_flux)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(inout) :: state
    type(amr_tree_flux_3d), intent(in) :: flux(:)
    real(dp), intent(in) :: root_spacing(3)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(inout) :: boundary_correction(:),maximum_correction
    logical, intent(out) :: ok
    integer, optional, intent(in) :: minimum_depth
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: synchronized_flux(:)
    integer, optional, intent(in) :: first_species,last_species
    type(amr_tree_state_3d), optional, intent(inout) :: auxiliary_state
    type(amr_tree_flux_3d), optional, intent(in) :: auxiliary_flux(:)
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: auxiliary_synchronized_flux(:)
    ! Keep the maximum correction scoped to the caller's gas components.
    ok=.false.
    if(.not.state%is_valid(tree)) return
    call reflux_amr_patch_tree_core(tree,state,flux,root_spacing,faces,boundary_correction,maximum_correction,ok, &
      minimum_depth,synchronized_flux,first_species,last_species, &
      auxiliary_state,auxiliary_flux,auxiliary_synchronized_flux,size(state%fields(1)%values,1))
  end subroutine reflux_amr_patch_tree_3d

  recursive subroutine reflux_amr_patch_tree_core(tree,state,flux,root_spacing,faces,boundary_correction,maximum_correction,ok, &
      minimum_depth,synchronized_flux,first_species,last_species, &
      auxiliary_state,auxiliary_flux,auxiliary_synchronized_flux,diagnostic_components)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(inout) :: state
    type(amr_tree_flux_3d), intent(in) :: flux(:)
    real(dp), intent(in) :: root_spacing(3)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(inout) :: boundary_correction(:),maximum_correction
    logical, intent(out) :: ok
    integer, optional, intent(in) :: minimum_depth
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: synchronized_flux(:)
    integer, optional, intent(in) :: first_species,last_species
    type(amr_tree_state_3d), optional, intent(inout) :: auxiliary_state
    type(amr_tree_flux_3d), optional, intent(in) :: auxiliary_flux(:)
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: auxiliary_synchronized_flux(:)
    integer, intent(in) :: diagnostic_components
    type(amr_tree_state_3d) :: candidate,coupled,coupled_auxiliary
    type(amr_tree_flux_3d), allocatable :: coupled_flux(:),coupled_effective(:),gas_effective(:),scalar_effective(:)
    real(dp), allocatable :: coupled_budget(:)
    integer :: bounds(2),coupled_flags,ios
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
    if(present(first_species).neqv.present(last_species)) return
    bounds=0
    if(present(first_species).and.present(last_species)) bounds=[first_species,last_species]
    if(any(bounds/=0)) then
      if(bounds(1)<2.or.bounds(2)<bounds(1).or.bounds(2)>nvar) return
    endif
    coupled_flags=0
    if(present(auxiliary_state)) coupled_flags=coupled_flags+1
    if(present(auxiliary_flux)) coupled_flags=coupled_flags+2
    if(present(auxiliary_synchronized_flux)) coupled_flags=coupled_flags+4
    if(coupled_flags/=0.and.coupled_flags/=3.and.coupled_flags/=7) return
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
    if(present(auxiliary_state)) then
      if(.not.auxiliary_state%is_valid(tree)) return
      if(size(auxiliary_state%fields(1)%values,1)/=1) return
      valid=nvar<huge(1).and.size(auxiliary_flux)==size(flux)
      if(.not.valid) return
      allocate(coupled%fields(size(flux)),coupled_auxiliary%fields(size(flux)),coupled_flux(size(flux)), &
        gas_effective(size(flux)),scalar_effective(size(flux)),coupled_budget(nvar+1),stat=ios)
      if(ios/=0) return
      valid=.true.
      do node=1,size(flux)
        n=tree%nodes(node)%shape
        allocate(coupled%fields(node)%values(nvar+1,n(1),n(2),n(3)),stat=ios)
        if(ios/=0) then
          valid=.false.;exit
        endif
        coupled%fields(node)%values(:nvar,:,:,:)=state%fields(node)%values
        coupled%fields(node)%values(nvar+1,:,:,:)=auxiliary_state%fields(node)%values(1,:,:,:)
        if(tree%nodes(node)%level<floor_depth) cycle
        if(.not.allocated(flux(node)%x).or..not.allocated(flux(node)%y).or..not.allocated(flux(node)%z)) then
          valid=.false.;exit
        endif
        if(.not.allocated(auxiliary_flux(node)%x).or..not.allocated(auxiliary_flux(node)%y).or. &
            .not.allocated(auxiliary_flux(node)%z)) then
          valid=.false.;exit
        endif
        if(any(shape(flux(node)%x)/=[nvar,n(1)+1,n(2),n(3)]).or. &
          any(shape(flux(node)%y)/=[nvar,n(1),n(2)+1,n(3)]).or. &
          any(shape(flux(node)%z)/=[nvar,n(1),n(2),n(3)+1]).or. &
          any(shape(auxiliary_flux(node)%x)/=[1,n(1)+1,n(2),n(3)]).or. &
          any(shape(auxiliary_flux(node)%y)/=[1,n(1),n(2)+1,n(3)]).or. &
          any(shape(auxiliary_flux(node)%z)/=[1,n(1),n(2),n(3)+1])) then
          valid=.false.;exit
        endif
        allocate(coupled_flux(node)%x(nvar+1,n(1)+1,n(2),n(3)), &
          coupled_flux(node)%y(nvar+1,n(1),n(2)+1,n(3)),coupled_flux(node)%z(nvar+1,n(1),n(2),n(3)+1),stat=ios)
        if(ios/=0) then
          valid=.false.;exit
        endif
        coupled_flux(node)%x(:nvar,:,:,:)=flux(node)%x;coupled_flux(node)%x(nvar+1,:,:,:)=auxiliary_flux(node)%x(1,:,:,:)
        coupled_flux(node)%y(:nvar,:,:,:)=flux(node)%y;coupled_flux(node)%y(nvar+1,:,:,:)=auxiliary_flux(node)%y(1,:,:,:)
        coupled_flux(node)%z(:nvar,:,:,:)=flux(node)%z;coupled_flux(node)%z(nvar+1,:,:,:)=auxiliary_flux(node)%z(1,:,:,:)
      enddo
      if(.not.valid) return
      coupled_budget=0;largest=0
      call reflux_amr_patch_tree_core(tree,coupled,coupled_flux,root_spacing,faces,coupled_budget, &
        largest,valid,minimum_depth=floor_depth,synchronized_flux=coupled_effective, &
        first_species=first_species,last_species=last_species,diagnostic_components=diagnostic_components)
      if(.not.valid) return
      ! Allocate all output copies before publishing either conserved field.
      candidate=state;coupled_auxiliary=auxiliary_state
      do node=1,size(flux)
        candidate%fields(node)%values=coupled%fields(node)%values(:nvar,:,:,:)
        coupled_auxiliary%fields(node)%values(1,:,:,:)=coupled%fields(node)%values(nvar+1,:,:,:)
        if(tree%nodes(node)%level<floor_depth) cycle
        gas_effective(node)%x=coupled_effective(node)%x(:nvar,:,:,:)
        gas_effective(node)%y=coupled_effective(node)%y(:nvar,:,:,:)
        gas_effective(node)%z=coupled_effective(node)%z(:nvar,:,:,:)
        scalar_effective(node)%x=coupled_effective(node)%x(nvar+1:nvar+1,:,:,:)
        scalar_effective(node)%y=coupled_effective(node)%y(nvar+1:nvar+1,:,:,:)
        scalar_effective(node)%z=coupled_effective(node)%z(nvar+1:nvar+1,:,:,:)
      enddo
      call move_alloc(candidate%fields,state%fields)
      call move_alloc(coupled_auxiliary%fields,auxiliary_state%fields)
      if(present(synchronized_flux)) call move_alloc(gas_effective,synchronized_flux)
      if(present(auxiliary_synchronized_flux)) call move_alloc(scalar_effective,auxiliary_synchronized_flux)
      boundary_correction=coupled_budget(:nvar);maximum_correction=largest;ok=.true.;return
    endif
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
              largest=max(largest,maxval(abs(correction(:diagnostic_components))))
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
    if(bounds(1)>0) then
      call balance_serial_reflux_species(tree,state,candidate,nvar,work,root_spacing,faces, &
        bounds(1),bounds(2),floor_depth,valid)
      if(.not.valid) return
      ! Match the MPI diagnostic without subtracting large energy values.
      do node=1,size(flux)
        if(tree%nodes(node)%level<floor_depth) cycle
        n=tree%nodes(node)%shape
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              largest=max(largest,maxval(abs(candidate%fields(node)%values(bounds(1):bounds(2),i,j,k)- &
                state%fields(node)%values(bounds(1):bounds(2),i,j,k))))
            enddo
          enddo
        enddo
      enddo
    endif
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
  end subroutine reflux_amr_patch_tree_core
  ! A reflux correction is limited as a complete conserved vector. The part
  ! returned to finer cells is funded by their connected composite inventory,
  ! with identical transfers recorded on every copy of every affected face.
  ! This serial traversal mirrors the MPI donor/face ordering and reuses the
  ! same topology graph and complete-vector funding kernel.
  subroutine balance_serial_reflux_species(tree,original,state,nvar,flux,root_spacing, &
      faces,first_species,last_species,minimum_depth,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: nvar,first_species,last_species,minimum_depth
    type(amr_tree_state_3d), intent(in) :: original
    type(amr_tree_state_3d), intent(inout) :: state
    type(amr_tree_flux_3d), intent(inout) :: flux(:)
    real(dp), intent(in) :: root_spacing(3)
    character(len=*), intent(in) :: faces(6)
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: candidate
    type(amr_tree_flux_3d), allocatable :: work(:)
    integer, allocatable :: targets(:,:),requests(:,:),edges(:,:),edge_side(:),ratios(:)
    integer(int64), allocatable :: lower(:,:),upper(:,:),origin(:,:),domain(:,:)
    real(dp), allocatable :: values(:,:),accepted(:,:),integrals(:,:),volumes(:),seed(:)
    real(dp) :: before(nvar),current(nvar),coarse(nvar),delta(nvar),returned(nvar),theta,volume,spacing(3),lattice(3)
    integer :: node,n(3),i,j,k,depth,nt,target(4),t,c,axis,side,attempt,e,u,v,peer,ios
    integer(int64) :: lo(3),hi(3),plane
    logical :: good,negative,funded,success
    ok=.false.
    good=first_species>=2.and.last_species>=first_species.and.last_species<=nvar
    if(.not.good) return
    negative=.false.
    do node=1,size(tree%nodes)
      if(tree%nodes(node)%level<minimum_depth) cycle
      n=tree%nodes(node)%shape
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            if(cell_covered(tree,node,[i,j,k])) cycle
            if(bad(state%fields(node)%values(:,i,j,k))) negative=.true.
          enddo
        enddo
      enddo
    enddo
    if(.not.negative) then
      ok=.true.;return
    endif
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,good)
    if(.not.good) return
    candidate=state;work=flux
    do depth=tree%level_count()-1,minimum_depth,-1
      spacing=root_spacing*real(domain(:,1),dp)/real(domain(:,depth+1),dp)
      volume=product(spacing)
      if(volume<=0.or..not.ieee_is_finite(volume)) return
      do node=1,size(tree%nodes)
        if(tree%nodes(node)%level/=depth) cycle
        n=tree%nodes(node)%shape;nt=0
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              if(cell_covered(tree,node,[i,j,k])) cycle
              if(bad(candidate%fields(node)%values(:,i,j,k))) nt=nt+1
            enddo
          enddo
        enddo
        if(nt>huge(1)/4) return
        if(nt==0) cycle
        allocate(targets(4,nt),stat=ios)
        if(ios/=0) return
        t=0
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              if(cell_covered(tree,node,[i,j,k])) cycle
              if(.not.bad(candidate%fields(node)%values(:,i,j,k))) cycle
              t=t+1;targets(:,t)=[node,i,j,k]
            enddo
          enddo
        enddo
        do t=1,nt
          target=targets(:,t);i=target(2);j=target(3);k=target(4)
          before=original%fields(node)%values(:,i,j,k);current=candidate%fields(node)%values(:,i,j,k)
          if(bad(before).or..not.all(ieee_is_finite(current))) return
          delta=current-before;theta=1
          do c=1,nvar
            if(c/=1.and.(c<first_species.or.c>last_species)) cycle
            if(delta(c)>=0) cycle
            theta=min(theta,before(c)/(-delta(c)))
          enddo
          if(theta<0.or.theta>1) return
          ! Recompute, rather than clip a component. Rounding at an exhausted
          ! budget is handled by reducing the common vector coefficient.
          do attempt=1,16
            coarse=before+theta*delta
            if(.not.bad(coarse)) exit
            if(theta==0) return
            theta=max(0._dp,nearest(theta,-1._dp))
          enddo
          if(bad(coarse)) return
          returned=(current-coarse)*volume;funded=.false.
          do axis=1,3
            do side=-1,1,2
              call build_amr_reflux_graph_3d(tree,root_spacing,faces,target,axis,side,requests,volumes, &
                edges,edge_side,seed,lower,upper,lattice,good)
              if(.not.good) cycle
              if(size(requests,2)>huge(1)/nvar.or.size(edges,2)>huge(1)/nvar) return
              if(allocated(values)) deallocate(values)
              allocate(values(nvar,size(requests,2)),stat=ios)
              if(ios/=0) return
              do c=1,size(requests,2)
                peer=requests(1,c)
                values(:,c)=candidate%fields(peer)%values(:,requests(2,c),requests(3,c),requests(4,c))
              enddo
              allocate(accepted(nvar,size(requests,2)),integrals(nvar,size(edges,2)),stat=ios)
              if(ios/=0) return
              success=.false.
              call spread_reflux_return(values,volumes,edges,seed,returned, &
                first_species,last_species,accepted,integrals,success)
              if(.not.success) then
                deallocate(accepted,integrals);cycle
              endif
              do c=1,size(requests,2)
                peer=requests(1,c)
                candidate%fields(peer)%values(:,requests(2,c),requests(3,c),requests(4,c))=accepted(:,c)
              enddo
              candidate%fields(node)%values(:,i,j,k)=coarse
              ! External supply at each seed plus graph's internal transfers.
              seed=seed/maxval(seed);seed=seed/sum(seed)
              c=maxloc(seed,dim=1);seed(c)=0;seed(c)=1-sum(seed)
              do c=1,size(requests,2)
                if(seed(c)==0) cycle
                lo=lower(:,c);hi=upper(:,c);plane=lo(axis)
                if(side<0) plane=hi(axis)
                call transfer_face(lo,hi,axis,plane,real(side,dp)*seed(c)*returned,good)
                if(.not.good) return
              enddo
              do e=1,size(edges,2)
                u=edges(1,e);v=edges(2,e);peer=abs(edge_side(e))
                lo=max(lower(:,u),lower(:,v));hi=min(upper(:,u),upper(:,v))
                plane=lower(peer,u)
                if(edge_side(e)>0) plane=upper(peer,u)
                call transfer_face(lo,hi,peer,plane,sign(1._dp,real(edge_side(e),dp))*integrals(:,e),good)
                if(.not.good) return
              enddo
              deallocate(accepted,integrals);funded=.true.;exit
            enddo
            if(funded) exit
          enddo
          if(.not.funded) return
        enddo
        deallocate(targets)
      enddo
    enddo
    state=candidate;flux=work;ok=.true.
  contains
    logical function bad(value)
      real(dp), intent(in) :: value(:)
      bad=value(1)<=0.or.any(value(first_species:last_species)<0)
    end function bad
    subroutine transfer_face(face_lo,face_hi,normal,position,amount,ready)
      integer(int64), intent(in) :: face_lo(3),face_hi(3),position
      integer, intent(in) :: normal
      real(dp), intent(in) :: amount(:)
      logical, intent(out) :: ready
      integer(int64) :: base(3),last(3),w,pl,low_cell(3),high_cell(3),overlap(3),face_base(3)
      integer :: p,lev,b,d,ii,jj,kk,ix(3),wrap
      real(dp) :: region_area,cell_area,weight,next(nvar)
      logical :: valid
      b=modulo(normal,3)+1;d=modulo(normal+1,3)+1
      region_area=real(face_hi(b)-face_lo(b),dp)*lattice(b)*real(face_hi(d)-face_lo(d),dp)*lattice(d)
      ready=region_area>0.and.ieee_is_finite(region_area)
      if(.not.ready) return
      valid=.true.
      do wrap=0,1
        pl=position
        if(wrap==1) then
          if(faces(2*normal-1)/='periodic') cycle
          if(position==0) then
            pl=domain(normal,size(ratios))
          else if(position==domain(normal,size(ratios))) then
            pl=0
          else
            cycle
          endif
        endif
        do p=1,size(tree%nodes)
          if(tree%nodes(p)%level<minimum_depth) cycle
          lev=tree%nodes(p)%level+1;w=domain(1,size(ratios))/domain(1,lev)
          base=origin(:,p)*w;last=base+int(tree%nodes(p)%shape,int64)*w
          if(pl<base(normal).or.pl>last(normal)) cycle
          if(modulo(pl-base(normal),w)/=0) cycle
          if(face_hi(b)<=base(b).or.face_lo(b)>=last(b).or.face_hi(d)<=base(d).or.face_lo(d)>=last(d)) cycle
          low_cell=(max(face_lo,base)-base)/w+1;high_cell=(min(face_hi,last)-1-base)/w+1
          low_cell(normal)=(pl-base(normal))/w+1;high_cell(normal)=low_cell(normal)
          cell_area=real(w,dp)*lattice(b)*real(w,dp)*lattice(d)
          if(cell_area<=0.or..not.ieee_is_finite(cell_area)) then
            valid=.false.;cycle
          endif
          do kk=int(low_cell(3)),int(high_cell(3))
            do jj=int(low_cell(2)),int(high_cell(2))
              do ii=int(low_cell(1)),int(high_cell(1))
                ix=[ii,jj,kk];face_base=base+int(ix-1,int64)*w
                overlap=min(face_hi,face_base+w)-max(face_lo,face_base)
                weight=(real(overlap(b),dp)*lattice(b)*real(overlap(d),dp)*lattice(d)/region_area)/cell_area
                next=get_face(work(p),normal,ix)+weight*amount
                if(.not.all(ieee_is_finite(next))) valid=.false.
                call set_face(work(p),normal,ix,next)
              enddo
            enddo
          enddo
        enddo
      enddo
      ready=valid
    end subroutine transfer_face
  end subroutine balance_serial_reflux_species

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
