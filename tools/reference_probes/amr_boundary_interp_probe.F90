! Standalone production-halo probe; output is input to the pinned AMReX probe.
program amr_boundary_interp_probe
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d, &
    initialize_amr_patch_tree_3d,prolong_amr_patch_tree_3d
  use amr_patch_tree_ghost_3d_mod, only: fill_amr_patch_tree_ghosts_3d
  implicit none
  type(amr_patch_tree_3d) :: tree
  type(amr_tree_state_3d) :: state,future
  type(amr_patch_3d) :: patch(2)
  real(dp) :: root(8,4,4,4),future_root(8,4,4,4),inlet(8,6),base(7),gradient(7),curve(7),r,blend
  real(dp), allocatable :: ghost(:,:,:,:)
  integer :: axis,side,layout,profile,id,x,y,z,index(3),lo(3),hi(3),shape_parent(3)
  integer :: target_node,factor,origin(3),fine_index(3),local_index(3),parent_lo(3),count,time_case
  character(len=32) :: faces(6)
  logical :: ok,include_internal
  character(len=32) :: argument
  base=[2._dp,3._dp,4._dp,5._dp,1000._dp,.6_dp,1.4_dp]
  gradient=[1._dp,.2_dp,.3_dp,.4_dp,2._dp,.3_dp,.7_dp]
  call get_command_argument(1,argument)
  include_internal=trim(argument)/='no_aux'
  id=0
  write(*,*) 48
  do axis=1,3
    do side=1,2
      do layout=1,2
        do profile=1,2
          do time_case=1,2
          id=id+1;curve=0
          if(profile==2) curve=[.03_dp,.01_dp,-.002_dp,.004_dp,.4_dp,.02_dp,.01_dp]
          do z=1,4
            do y=1,4
              do x=1,4
                index=[x,y,z];r=real(index(axis),dp)-.5_dp
                root(:7,x,y,z)=base+gradient*r+curve*r*r
                root(8,x,y,z)=300+.2_dp*r
              enddo
            enddo
          enddo
          lo=2;hi=3;lo(axis)=1;hi(axis)=2
          if(side==2) then
            lo(axis)=3;hi(axis)=4
          endif
          if(layout==1) then
            count=1;target_node=2;factor=2;origin=(lo-1)*factor
            call initialize_amr_patch_3d(4,4,4,lo(1),hi(1),lo(2),hi(2),lo(3),hi(3),2,patch(1),ok)
            if(.not.ok) error stop 'Invalid wide stencil geometry'
          else
            count=2;target_node=3;factor=6
            hi(axis)=lo(axis)
            if(side==2) then
              lo(axis)=4;hi(axis)=4
            endif
            parent_lo=lo;shape_parent=3*(hi-lo+1)
            call initialize_amr_patch_3d(4,4,4,lo(1),hi(1),lo(2),hi(2),lo(3),hi(3),3,patch(1),ok)
            if(.not.ok) error stop 'Invalid narrow parent geometry'
            lo=1;hi=1
            if(side==2) then
              lo(axis)=3;hi(axis)=3
            endif
            origin=(parent_lo-1)*6+(lo-1)*2
            call initialize_amr_patch_3d(shape_parent(1),shape_parent(2),shape_parent(3), &
              lo(1),hi(1),lo(2),hi(2),lo(3),hi(3),2,patch(2),ok)
            if(.not.ok) error stop 'Invalid narrow child geometry'
          endif
          call initialize_amr_patch_tree_3d([4,4,4],[(x,x=1,count)],patch(:count),tree,ok)
          if(.not.ok) error stop 'Invalid tree'
          call prolong_amr_patch_tree_3d(tree,root,state,ok)
          if(.not.ok) error stop 'Invalid state'
          future=state;blend=0
          if(time_case==2) blend=.35_dp
          future_root=root
          do z=1,4
            do y=1,4
              do x=1,4
                index=[x,y,z];r=real(index(axis),dp)-.5_dp
                future_root(2,x,y,z)=root(2,x,y,z)+5*r*(4-r)
              enddo
            enddo
          enddo
          call prolong_amr_patch_tree_3d(tree,future_root,future,ok)
          if(.not.ok) error stop 'Invalid future endpoint state'
          faces='outflow' ;faces(2*axis-1:2*axis)='inflow'
          inlet(:7,:)=spread(base,2,6);inlet(:7,2*axis)=base+4*gradient+16*curve
          inlet(8,:)=300;inlet(8,2*axis)=300.8_dp
          call fill_amr_patch_tree_ghosts_3d(tree,state,target_node,faces,ghost,ok, &
            coarse_limiter='linear',species_count=2,inflow_state=inlet,ghost_width=1, &
            coarse_future=future,coarse_fractions=spread(blend,1,count),limit_internal_energy=include_internal)
          if(.not.ok) error stop 'Halo rejected'
          local_index=2;local_index(modulo(axis,3)+1)=1
          if(side==2) local_index(axis)=tree%nodes(target_node)%shape(axis)+1
          fine_index=origin+local_index-2
          write(*,*) id,axis,side,layout,profile,factor,blend
          write(*,*) origin,tree%nodes(target_node)%shape,fine_index
          do z=1,4
            do y=1,4
              do x=1,4
                write(*,'(8(es25.17e3,1x))') root(:,x,y,z)
              enddo
            enddo
          enddo
          do z=1,4
            do y=1,4
              do x=1,4
                write(*,'(8(es25.17e3,1x))') future%fields(1)%values(:,x,y,z)
              enddo
            enddo
          enddo
          write(*,'(8(es25.17e3,1x))') inlet(:,2*axis-1)
          write(*,'(8(es25.17e3,1x))') inlet(:,2*axis)
          write(*,'(8(es25.17e3,1x))') ghost(:,local_index(1),local_index(2),local_index(3))
          enddo
        enddo
      enddo
    enddo
  enddo
end program amr_boundary_interp_probe
