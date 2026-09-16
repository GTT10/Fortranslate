program test_amr_patch_tree_ghost_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d,prolong_pcm_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d, &
    initialize_amr_patch_tree_3d,prolong_amr_patch_tree_3d
  use amr_patch_tree_ghost_3d_mod, only: fill_amr_patch_tree_ghosts_3d
  use simulation_config_reactive_3d_mod, only: fill_reactive_boundary_ghosts_3d
  implicit none
  type(amr_patch_tree_3d) :: tree,deep
  type(amr_tree_state_3d) :: state,deep_state
  type(amr_patch_3d) :: patches(4),whole,chain(63)
  real(dp) :: root(6,4,4,4),inlet(6,6),inlet_t(6),velocity(3,6)
  real(dp), allocatable :: dense(:,:,:,:),next(:,:,:,:),oracle(:,:,:,:),ghost(:,:,:,:),saved(:,:,:,:),t(:,:,:),gt(:,:,:)
  integer :: origins(3,5),node,i,j,k,level,n,lo(3),hi(3),mode
  character(len=32) :: faces(6)
  logical :: ok,viscous
  call initialize_amr_patch_3d(4,4,4,1,2,2,3,1,4,2,patches(1),ok)
  call require(ok,'first root child')
  call initialize_amr_patch_3d(4,4,4,3,4,2,3,1,4,2,patches(2),ok)
  call require(ok,'second root child')
  call initialize_amr_patch_3d(4,4,8,3,4,1,2,2,3,2,patches(3),ok)
  call require(ok,'left grandchild')
  call initialize_amr_patch_3d(4,4,8,1,2,1,2,2,3,2,patches(4),ok)
  call require(ok,'right grandchild')
  call initialize_amr_patch_tree_3d([4,4,4],[1,1,2,3],patches,tree,ok)
  call require(ok,'tree')
  do k=1,4
    do j=1,4
      do i=1,4
        root(:,i,j,k)=[2.0_dp,real(i,dp),real(j,dp),real(k,dp),1000.0_dp,2.0_dp]
      end do
    end do
  end do
  call prolong_amr_patch_tree_3d(tree,root,state,ok)
  call require(ok,'state')
  origins(:,1)=[0,0,0];origins(:,2)=[0,2,0];origins(:,3)=[4,2,0]
  origins(:,4)=[4,4,2];origins(:,5)=[8,4,2]
  ! Make peer values distinct from their coarse ancestors.
  do node=2,5
    do k=1,tree%nodes(node)%shape(3)
      do j=1,tree%nodes(node)%shape(2)
        do i=1,tree%nodes(node)%shape(1)
          state%fields(node)%values(:,i,j,k)=[2.0_dp,real(100*node+i,dp), &
            real(10*node+j,dp),real(node+k,dp),1.e6_dp,2.0_dp]
        end do
      end do
    end do
  end do
  inlet=spread([3.0_dp,4.0_dp,5.0_dp,6.0_dp,2.e6_dp,3.0_dp],2,6);inlet_t=300
  do mode=1,3
    faces='periodic';velocity=0;viscous=mode==3
    if (mode>1) then
      faces=[character(len=32)::'slip_wall','outflow','inflow','no_slip_wall','isothermal_wall','inflow']
      velocity(1,4)=.2_dp;velocity(2,5)=.3_dp
    end if
    dense=root
    do level=0,2
      n=4*2**level
      ! Independent oracle: prolong the complete dense grid, then overlay peer boxes.
      if (level>0) then
        call initialize_amr_patch_3d(n/2,n/2,n/2,1,n/2,1,n/2,1,n/2,2,whole,ok)
        call require(ok,'dense oracle box')
        allocate(next(6,n,n,n))
        call prolong_pcm_3d(dense,whole,next,ok)
        call require(ok,'dense oracle prolongation')
        call move_alloc(next,dense)
        do node=2,5
          if (tree%nodes(node)%level/=level) cycle
          lo=origins(:,node)+1;hi=origins(:,node)+tree%nodes(node)%shape
          dense(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=state%fields(node)%values
        end do
      end if
      allocate(t(n,n,n),oracle(6,n+4,n+4,n+4),gt(n+4,n+4,n+4));t=300
      call fill_reactive_boundary_ghosts_3d(dense,t,faces,oracle,gt,ok,viscous,inlet,inlet_t,velocity)
      call require(ok,'dense physical/periodic oracle')
      do node=1,5
        if (tree%nodes(node)%level/=level) cycle
        call fill_amr_patch_tree_ghosts_3d(tree,state,node,faces,ghost,ok,viscous,inlet,velocity)
        call require(ok,'tree ghost fill')
        lo=origins(:,node)+1;hi=origins(:,node)+tree%nodes(node)%shape+4
        call require(all(ghost==oracle(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))), &
          'all halo/interior values equal dense oracle including cross-parent peers and ancestor fallback')
      end do
      deallocate(t,oracle,gt)
    end do
    deallocate(dense)
  end do
  saved=ghost
  call fill_amr_patch_tree_ghosts_3d(tree,state,5,faces,ghost,ok)
  call require(.not.ok.and.all(ghost==saved),'missing reservoir preserves output')
  faces='periodic';faces(1)='outflow'
  call fill_amr_patch_tree_ghosts_3d(tree,state,5,faces,ghost,ok)
  call require(.not.ok.and.all(ghost==saved),'unpaired periodic face rejected')
  faces='periodic'
  call fill_amr_patch_tree_ghosts_3d(tree,state,6,faces,ghost,ok)
  call require(.not.ok.and.all(ghost==saved),'invalid node rejected')
  state%fields(1)%values(1,1,1,1)=ieee_value(0.0_dp,ieee_quiet_nan)
  call fill_amr_patch_tree_ghosts_3d(tree,state,5,faces,ghost,ok)
  call require(.not.ok.and.all(ghost==saved),'invalid ancestor preserves output')
  do i=1,63
    call initialize_amr_patch_3d(2,2,2,1,1,1,1,1,1,2,chain(i),ok)
  end do
  call initialize_amr_patch_tree_3d([2,2,2],[(i,i=1,63)],chain,deep,ok)
  call require(ok,'deep compact topology')
  call prolong_amr_patch_tree_3d(deep,root(:,1:2,1:2,1:2),deep_state,ok)
  call require(ok,'deep compact fields')
  call fill_amr_patch_tree_ghosts_3d(deep,deep_state,64,faces,ghost,ok)
  call require(.not.ok.and.all(ghost==saved),'global integer-coordinate overflow rejects before allocation')
  print '(a)','test_amr_patch_tree_ghost_3d: PASS'
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (condition) return
    print '(a)','FAIL: '//label
    error stop 1
  end subroutine
end program
