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
        if(mode==1) then
          call fill_amr_patch_tree_ghosts_3d(tree,state,node,faces,ghost,ok,ghost_width=7)
          call require(ok.and.all(shape(ghost)==[6,tree%nodes(node)%shape+14]),'wide periodic ghost shape')
          do k=1,size(ghost,4)
            do j=1,size(ghost,3)
              do i=1,size(ghost,2)
                lo=modulo(origins(:,node)+[i,j,k]-8,n)+1
                call require(all(ghost(:,i,j,k)==dense(:,lo(1),lo(2),lo(3))), &
                  'wide periodic donors match dense oracle across peer and ancestor boundaries')
              end do
            end do
          end do
          saved=ghost
          call fill_amr_patch_tree_ghosts_3d(tree,state,node,faces,ghost,ok,ghost_width=0)
          call require(.not.ok.and.all(ghost==saved),'invalid ghost width preserves output')
        end if
      end do
      deallocate(t,oracle,gt)
    end do
    deallocate(dense)
  end do
  call wide_wall_checks()
  call temporal_checks()
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
  subroutine wide_wall_checks()
    type(amr_patch_tree_3d) :: box
    type(amr_tree_state_3d) :: field
    real(dp), allocatable :: halo(:,:,:,:),dense_halo(:,:,:,:),dense_t(:,:,:)
    real(dp) :: temperature(4,4,4)
    real(dp) :: walls(3,6),expected(6),momentum(3),offset
    character(len=32) :: bc(6)
    integer :: x,donor,slab,r,ix,iy,iz,kind
    logical :: success
    call initialize_amr_patch_tree_3d([4,4,4],[integer::],patches(:0),box,success)
    call require(success,'wide wall root geometry')
    call prolong_amr_patch_tree_3d(box,root,field,success)
    call require(success,'wide wall root state')
    temperature=800
    allocate(dense_halo(6,30,30,30),dense_t(30,30,30))
    do kind=1,3
      bc='periodic';bc(1:2)='no_slip_wall';walls=0
      walls(2,1)=1;walls(2,2)=3
      if(kind==2) then
        bc(1:2)='slip_wall';walls=0
      else if(kind==3) then
        bc(2)='outflow';walls(:,2)=0
      end if
      call fill_amr_patch_tree_ghosts_3d(box,field,1,bc,halo,success, &
        viscous_wall=.true.,wall_velocity=walls,ghost_width=13)
      call require(success,'wide wall halo crosses domain repeatedly')
      call fill_reactive_boundary_ghosts_3d(root,temperature,bc,dense_halo,dense_t,success, &
        viscous_wall=.true.,wall_velocity=walls,ghost_width=13)
      call require(success,'wide regular wall ghosts')
      call require(all(dense_halo==halo).and.all(dense_t==800),'wide regular and tree donors agree exactly')
      do iz=1,30
        do iy=1,30
          do ix=1,30
            x=ix-14;slab=floor(real(x,dp)/4);r=modulo(x,4)
            donor=r+1
            if(modulo(slab,2)==1) donor=4-r
            if(kind==3) donor=max(1,min(4,x+1))
            if(kind==3.and.x<0) donor=min(4,-x)
            expected=root(:,donor,modulo(iy-14,4)+1,modulo(iz-14,4)+1)
            momentum=expected(2:4)
            if(kind<=2) then
              if(modulo(slab,2)==1) then
                expected(2)=-expected(2)
                if(kind==1) expected(3:4)=-expected(3:4)
              end if
              if(kind==1) then
                ! Analytic extension between walls v_y=1 and v_y=3:
                ! every pair of reflections translates v_y by four.
                offset=4._dp*floor(real(slab,dp)/2)
                if(modulo(slab,2)==1) offset=offset+6
                expected(3)=expected(3)+expected(1)*offset
              end if
            else if(x<0) then
              expected(2:4)=-expected(2:4)
              expected(3)=expected(3)+2*expected(1)
            end if
            expected(5)=expected(5)+(sum(expected(2:4)**2)-sum(momentum**2))/(2*expected(1))
            call require(all(halo(:,ix,iy,iz)==expected),'analytic repeated wall donor, momentum and energy')
          end do
        end do
      end do
    end do
    bc='periodic';walls=0
    call fill_reactive_boundary_ghosts_3d(root,temperature,bc,dense_halo,dense_t,success,ghost_width=0)
    call require(.not.success,'zero regular ghost width rejected')
    call fill_reactive_boundary_ghosts_3d(root,temperature,bc,dense_halo,dense_t,success,ghost_width=129)
    call require(.not.success,'excessive regular ghost width rejected')
  end subroutine wide_wall_checks
  subroutine temporal_checks()
    type(amr_tree_state_3d) :: future,incomplete
    real(dp), allocatable :: full(:,:,:,:),refined(:,:,:,:),expected(:,:,:,:),actual(:,:,:,:),snapshot(:,:,:,:)
    real(dp), allocatable :: temperatures(:,:,:),expected_t(:,:,:)
    real(dp) :: alpha(2),bad_alpha(2),wall(3,6)
    character(len=32) :: boundary_faces(6)
    integer :: endpoint,physical,depth,peer,side,target_depth,lower_box(3),upper_box(3)
    logical :: success
    allocate(future%fields(5))
    do peer=1,3
      future%fields(peer)%values=state%fields(peer)%values
      future%fields(peer)%values(1,:,:,:)=future%fields(peer)%values(1,:,:,:)+.5_dp
      future%fields(peer)%values(2:4,:,:,:)=future%fields(peer)%values(2:4,:,:,:)+16.0_dp
      future%fields(peer)%values(5,:,:,:)=future%fields(peer)%values(5,:,:,:)+256.0_dp
      future%fields(peer)%values(6,:,:,:)=future%fields(peer)%values(6,:,:,:)+.5_dp
    end do
    ! Future peer/finer fields 4 and 5 deliberately remain unallocated.
    do endpoint=1,3
      select case(endpoint)
      case(1)
        alpha=0
      case(2)
        alpha=1
      case(3)
        alpha=[.25_dp,.75_dp]
      end select
      do physical=1,2
        boundary_faces='periodic';wall=0
        if (physical==2) then
          boundary_faces=[character(len=32)::'slip_wall','outflow','inflow','no_slip_wall','isothermal_wall','inflow']
          wall=velocity
        end if
        do target_depth=1,2
        full=(1-alpha(1))*root+alpha(1)*future%fields(1)%values
        do depth=1,target_depth
          side=4*2**depth
          call initialize_amr_patch_3d(side/2,side/2,side/2,1,side/2,1,side/2,1,side/2,2,whole,success)
          call require(success,'temporal dense geometry')
          allocate(refined(6,side,side,side))
          call prolong_pcm_3d(full,whole,refined,success)
          call require(success,'temporal dense prolongation')
          call move_alloc(refined,full)
          do peer=2,5
            if (tree%nodes(peer)%level/=depth) cycle
            lower_box=origins(:,peer)+1;upper_box=origins(:,peer)+tree%nodes(peer)%shape
            if (depth<target_depth) then
              full(:,lower_box(1):upper_box(1),lower_box(2):upper_box(2),lower_box(3):upper_box(3))= &
                (1-alpha(2))*state%fields(peer)%values+alpha(2)*future%fields(peer)%values
            else
              full(:,lower_box(1):upper_box(1),lower_box(2):upper_box(2),lower_box(3):upper_box(3))= &
                state%fields(peer)%values
            end if
          end do
        end do
        allocate(temperatures(side,side,side),expected(6,side+4,side+4,side+4), &
          expected_t(side+4,side+4,side+4));temperatures=300
        call fill_reactive_boundary_ghosts_3d(full,temperatures,boundary_faces,expected,expected_t, &
          success,.true.,inlet,inlet_t,wall)
        call require(success,'temporal dense boundary oracle')
        do peer=2,5
          if (tree%nodes(peer)%level/=target_depth) cycle
          call fill_amr_patch_tree_ghosts_3d(tree,state,peer,boundary_faces,actual,success, &
            .true.,inlet,wall,future,alpha(:target_depth))
          call require(success,'time-interpolated tree fill')
          lower_box=origins(:,peer)+1;upper_box=origins(:,peer)+tree%nodes(peer)%shape+4
          call require(all(actual==expected(:,lower_box(1):upper_box(1),lower_box(2):upper_box(2), &
            lower_box(3):upper_box(3))),'temporal endpoints and per-donor-level fractions match dense oracle')
          call require(all(actual(:,3:tree%nodes(peer)%shape(1)+2,3:tree%nodes(peer)%shape(2)+2, &
            3:tree%nodes(peer)%shape(3)+2)==state%fields(peer)%values), &
            'temporal interpolation leaves current fine interior exact')
        end do
        deallocate(full,temperatures,expected,expected_t)
        end do
      end do
    end do
    snapshot=actual;boundary_faces='periodic'
    call fill_amr_patch_tree_ghosts_3d(tree,state,5,boundary_faces,actual,success,coarse_future=future)
    call require(.not.success.and.all(actual==snapshot),'missing time fractions preserve output')
    call fill_amr_patch_tree_ghosts_3d(tree,state,5,boundary_faces,actual,success,coarse_fractions=alpha)
    call require(.not.success.and.all(actual==snapshot),'missing future state preserves output')
    call fill_amr_patch_tree_ghosts_3d(tree,state,5,boundary_faces,actual,success, &
      coarse_future=future,coarse_fractions=alpha(:1))
    call require(.not.success.and.all(actual==snapshot),'wrong fraction count rejected')
    bad_alpha=[-.1_dp,.5_dp]
    call fill_amr_patch_tree_ghosts_3d(tree,state,5,boundary_faces,actual,success, &
      coarse_future=future,coarse_fractions=bad_alpha)
    call require(.not.success.and.all(actual==snapshot),'time extrapolation rejected')
    bad_alpha=[.5_dp,ieee_value(0.0_dp,ieee_quiet_nan)]
    call fill_amr_patch_tree_ghosts_3d(tree,state,5,boundary_faces,actual,success, &
      coarse_future=future,coarse_fractions=bad_alpha)
    call require(.not.success.and.all(actual==snapshot),'nonfinite fraction rejected')
    incomplete=future;deallocate(incomplete%fields(3)%values)
    call fill_amr_patch_tree_ghosts_3d(tree,state,5,boundary_faces,actual,success, &
      coarse_future=incomplete,coarse_fractions=alpha)
    call require(.not.success.and.all(actual==snapshot),'missing coarse peer future rejected')
  end subroutine temporal_checks
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (condition) return
    print '(a)','FAIL: '//label
    error stop 1
  end subroutine
end program
