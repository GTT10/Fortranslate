program test_amr_patch_tree_reflux_3d
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d,amr_tree_integer_layout_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d,reflux_amr_patch_tree_3d
  implicit none
  type(amr_patch_tree_3d) :: tree
  type(amr_tree_state_3d) :: state,saved,uncorrected,staged
  type(amr_patch_3d) :: patches(5)
  type(amr_tree_flux_3d), allocatable :: flux(:),merged(:),scoped(:),effective(:),saved_effective(:)
  integer(int64), allocatable :: origin(:,:),domain(:,:)
  integer, allocatable :: ratios(:)
  real(dp) :: root(2,4,4,4),h(3,4),spacing(3),initial(2),final(2),expected(2),root_budget(2),budget(2),largest
  real(dp) :: q(3),v(2),factor,area,old,scoped_budget(2),scoped_max,reference_budget(2)
  integer :: floor_depth,node,n(3),level,i,j,k,mode,axis,idx(3),hi(3),child,lo_child(3),hi_child(3)
  logical :: ok,covered
  character(len=32) :: faces(6)
  call initialize_amr_patch_3d(4,4,4,1,2,2,3,1,4,2,patches(1),ok)
  call initialize_amr_patch_3d(4,4,4,3,4,2,3,1,4,2,patches(2),ok)
  call initialize_amr_patch_3d(4,4,8,3,4,1,2,1,2,2,patches(3),ok)
  call initialize_amr_patch_3d(4,4,8,1,2,1,2,1,2,2,patches(4),ok)
  call initialize_amr_patch_3d(4,4,4,3,4,1,2,1,2,2,patches(5),ok)
  call initialize_amr_patch_tree_3d([4,4,4],[1,1,2,3,4],patches,tree,ok)
  call require(ok,'four-level branched topology')
  call amr_tree_integer_layout_3d(tree,origin,domain,ratios,ok)
  call require(ok,'global geometry')
  spacing=[.25_dp,.5_dp,.75_dp];h(:,1)=spacing
  do level=2,4
    h(:,level)=h(:,level-1)/ratios(level)
  end do
  root(1,:,:,:)=100;root(2,:,:,:)=200
  initial=[sum(root(1,:,:,:)),sum(root(2,:,:,:))]*product(spacing)
  allocate(flux(6))
  do mode=1,2
    faces='outflow'
    if (mode==2) faces='periodic'
    call prolong_amr_patch_tree_3d(tree,root,state,ok)
    call require(ok,'initialize conserved state')
    do node=1,6
      n=tree%nodes(node)%shape;level=tree%nodes(node)%level+1
      if (allocated(flux(node)%x)) deallocate(flux(node)%x,flux(node)%y,flux(node)%z)
      allocate(flux(node)%x(2,n(1)+1,n(2),n(3)),flux(node)%y(2,n(1),n(2)+1,n(3)), &
        flux(node)%z(2,n(1),n(2),n(3)+1))
      do axis=1,3
        hi=n;hi(axis)=hi(axis)+1
        do k=1,hi(3)
          do j=1,hi(2)
            do i=1,hi(1)
              idx=[i,j,k];q=(real(origin(:,node),dp)+real(idx,dp)-.5_dp)*h(:,level)
              q(axis)=q(axis)-.5_dp*h(axis,level)
              factor=.001_dp*level*(1+sum(q**2))
              if (mode==2) factor=.001_dp*level*(1+q(axis)*(4*spacing(axis)-q(axis)))
              v=[factor,-.7_dp*factor]
              select case(axis)
              case(1)
                flux(node)%x(:,i,j,k)=v
              case(2)
                flux(node)%y(:,i,j,k)=v
              case(3)
                flux(node)%z(:,i,j,k)=v
              end select
            end do
          end do
        end do
      end do
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            state%fields(node)%values(:,i,j,k)=state%fields(node)%values(:,i,j,k)- &
              (flux(node)%x(:,i+1,j,k)-flux(node)%x(:,i,j,k))/h(1,level)- &
              (flux(node)%y(:,i,j+1,k)-flux(node)%y(:,i,j,k))/h(2,level)- &
              (flux(node)%z(:,i,j,k+1)-flux(node)%z(:,i,j,k))/h(3,level)
          end do
        end do
      end do
    end do
    ! Independently integrate only exposed leaf-cell faces on the physical domain.
    expected=0;root_budget=0
    if (mode==1) then
      do node=1,6
        n=tree%nodes(node)%shape;level=tree%nodes(node)%level+1
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              idx=[i,j,k];covered=.false.
              do child=node+1,6
                if (tree%nodes(child)%parent/=node) cycle
                lo_child=[tree%nodes(child)%patch%coarse_i_lower,tree%nodes(child)%patch%coarse_j_lower, &
                  tree%nodes(child)%patch%coarse_k_lower]
                hi_child=[tree%nodes(child)%patch%coarse_i_upper,tree%nodes(child)%patch%coarse_j_upper, &
                  tree%nodes(child)%patch%coarse_k_upper]
                if (all(idx>=lo_child).and.all(idx<=hi_child)) covered=.true.
              end do
              v=0
              do axis=1,3
                area=product(h(:,level))/h(axis,level)
                if (origin(axis,node)+idx(axis)-1==0) v=v+area*cell_face(node,axis,idx,.false.)
                if (origin(axis,node)+idx(axis)==domain(axis,level)) v=v-area*cell_face(node,axis,idx,.true.)
              end do
              if (node==1) root_budget=root_budget+v
              if (.not.covered) expected=expected+v
            end do
          end do
        end do
      end do
    end if
    uncorrected=state
    budget=-17;largest=-19
    call reflux_amr_patch_tree_3d(tree,state,flux,spacing,faces,budget,largest,ok)
    call require(ok,'tree flux correction')
    call require(largest>0,'nontrivial corrections exercised')
    call composite_integrals_amr_patch_tree_3d(tree,state,spacing,final,ok)
    call require(ok,'corrected composite inventory')
    if (maxval(abs(final-initial-expected))>2.e-11_dp) print *,final-initial-expected
    call require(maxval(abs(final-initial-expected))<2.e-11_dp,'composite balance matches exposed finest boundary flux')
    call require(maxval(abs(budget-(expected-root_budget)))<2.e-13_dp,'root boundary ledger replacement')
    call require(maxval(abs([sum(state%fields(1)%values(1,:,:,:)),sum(state%fields(1)%values(2,:,:,:))]* &
      product(spacing)-final))<2.e-11_dp,'all corrections and restriction reach the root')
    reference_budget=budget;staged=uncorrected;merged=flux
    do floor_depth=2,0,-1
      scoped=merged
      do node=1,size(tree%nodes)
        if (tree%nodes(node)%level>=floor_depth) cycle
        deallocate(scoped(node)%x,scoped(node)%y,scoped(node)%z)
      end do
      scoped_budget=-17;scoped_max=-19
      call reflux_amr_patch_tree_3d(tree,staged,scoped,spacing,faces,scoped_budget,scoped_max,ok, &
        minimum_depth=floor_depth,synchronized_flux=effective)
      call require(ok,'partial interval synchronization')
      do node=1,size(tree%nodes)
        if (tree%nodes(node)%level<floor_depth) then
          call require(all(staged%fields(node)%values==uncorrected%fields(node)%values),'coarse interval untouched')
          call require(.not.allocated(effective(node)%x),'inactive effective flux absent')
        else
          merged(node)=effective(node)
        end if
      end do
      if (floor_depth>0) call require(all(scoped_budget==0),'root budget deferred')
    end do
    call require(maxval(abs(scoped_budget-reference_budget))<2e-13_dp,'composed interval boundary budget')
    do node=1,size(tree%nodes)
      call require(maxval(abs(staged%fields(node)%values-state%fields(node)%values))<2e-13_dp,'composed interval state')
    end do
    call reflux_amr_patch_tree_3d(tree,staged,merged,spacing,faces,scoped_budget,scoped_max,ok)
    call require(ok,'reuse synchronized flux')
    call require(scoped_max<2e-13_dp.and.maxval(abs(scoped_budget))<2e-13_dp,'no duplicate correction')

  end do
  saved=state;saved_effective=effective;budget=17;largest=19
  call reflux_amr_patch_tree_3d(tree,state,flux,spacing,faces,budget,largest,ok, &
    minimum_depth=-1,synchronized_flux=effective)
  call rollback("invalid minimum depth")
  old=flux(5)%x(1,1,4,4);flux(5)%x(1,1,4,4)=old+1
  call reflux_amr_patch_tree_3d(tree,state,flux,spacing,faces,budget,largest,ok,synchronized_flux=effective)
  call rollback('inconsistent shared peer flux')
  flux(5)%x(1,1,4,4)=old
  flux(1)%z(1,1,1,1)=ieee_value(0.0_dp,ieee_quiet_nan)
  call reflux_amr_patch_tree_3d(tree,state,flux,spacing,faces,budget,largest,ok,synchronized_flux=effective)
  call rollback('nonfinite flux')
  print '(a)','test_amr_patch_tree_reflux_3d: PASS'
contains
  function cell_face(which,a,index,upper) result(value)
    integer, intent(in) :: which,a,index(3)
    logical, intent(in) :: upper
    real(dp) :: value(2)
    integer :: p(3)
    p=index
    if (upper) p(a)=p(a)+1
    select case(a)
    case(1)
      value=flux(which)%x(:,p(1),p(2),p(3))
    case(2)
      value=flux(which)%y(:,p(1),p(2),p(3))
    case(3)
      value=flux(which)%z(:,p(1),p(2),p(3))
    end select
  end function
  subroutine rollback(label)
    character(len=*), intent(in) :: label
    call require(.not.ok.and.all(budget==17).and.largest==19,'reject and preserve diagnostics: '//label)
    do node=1,6
      call require(all(state%fields(node)%values==saved%fields(node)%values),'preserve all levels: '//label)
      call require(all(effective(node)%x==saved_effective(node)%x),'preserve synchronized x flux')
      call require(all(effective(node)%y==saved_effective(node)%y),'preserve synchronized y flux')
      call require(all(effective(node)%z==saved_effective(node)%z),'preserve synchronized z flux')
    end do
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (condition) return
    print '(a)','FAIL: '//label
    error stop 1
  end subroutine
end program
