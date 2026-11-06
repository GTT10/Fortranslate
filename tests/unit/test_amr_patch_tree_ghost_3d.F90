program test_amr_patch_tree_ghost_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use precision_mod, only: dp
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d,prolong_pcm_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d, &
    initialize_amr_patch_tree_3d,prolong_amr_patch_tree_3d
  use amr_patch_tree_ghost_3d_mod, only: fill_amr_patch_tree_ghosts_3d,apply_amr_tree_pressure_ghosts_3d
  use simulation_config_reactive_3d_mod, only: fill_reactive_boundary_ghosts_3d, &
    fill_mach_inlet_boundary_ghosts_3d,mach_inlet_ghost_3d,characteristic_pressure_ghost_3d
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
  call mach_inlet_checks()
  call rotating_inlet_checks()
  call rotating_pressure_checks()
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
  call linear_inlet_checks()
  call linear_short_inlet_checks()
  print '(a)','test_amr_patch_tree_ghost_3d: PASS'
contains
  subroutine rotating_pressure_checks()
    type(nasa7_species), allocatable :: species(:)
    type(amr_tree_state_3d) :: packed
    type(rotating_frame_3d) :: frame
    real(dp), allocatable :: buffer(:,:,:,:),kept(:,:,:,:)
    real(dp) :: lower(3),spacing(3),position(3),potential,q(12),u(12),expected(12),donor(12)
    real(dp) :: temperature,sound,pressure(6),expected_t
    character(len=32) :: saved_faces(6)
    integer :: nv,p,x,y,z,layer
    saved_faces=faces
    call load_h2o2_elementary_thermo(species,ok)
    call require(ok,'rotating pressure thermo')
    nv=5+size(species)
    q=[.3_dp,10._dp,3._dp,-2._dp,1.e5_dp,.01_dp,.001_dp,.001_dp,.2_dp,.001_dp,.05_dp,.737_dp]
    call reactive_primitive_to_conserved(species,q,u,temperature,sound,ok)
    call require(ok,'rotating pressure initial state')
    frame%axis=3;frame%angular_speed=400._dp
    frame%axis_point=[.0004_dp,-.0007_dp,0._dp]
    lower=[-.002_dp,.003_dp,.001_dp];spacing=[.001_dp,.0015_dp,.002_dp]
    allocate(packed%fields(5))
    do p=1,5
      allocate(packed%fields(p)%values(nv+2,tree%nodes(p)%shape(1),tree%nodes(p)%shape(2),tree%nodes(p)%shape(3)))
      do z=1,tree%nodes(p)%shape(3)
      do y=1,tree%nodes(p)%shape(2)
      do x=1,tree%nodes(p)%shape(1)
        position=lower+(real(origins(:,p)+[x,y,z],dp)-.5_dp)*spacing/real(2**tree%nodes(p)%level,dp)
        call rotating_potential_3d(frame,position,potential,ok)
        call require(ok,'rotating pressure cell potential')
        packed%fields(p)%values(:nv,x,y,z)=u
        packed%fields(p)%values(5,x,y,z)=u(5)-u(1)*potential
        packed%fields(p)%values(nv+1,x,y,z)=temperature
        packed%fields(p)%values(nv+2,x,y,z)=u(5)-.5_dp*sum(u(2:4)**2)/u(1)
      enddo
      enddo
      enddo
    enddo
    faces='outflow';pressure=0;pressure(2)=1.01e5_dp
    do p=1,3,2
      call fill_amr_patch_tree_ghosts_3d(tree,packed,p,faces,buffer,ok,ghost_width=3, &
        frame=frame,domain_lower=lower,root_spacing=spacing)
      call require(ok,'rotating pressure donor halo')
      kept=buffer
      call apply_amr_tree_pressure_ghosts_3d(tree,species,p,faces,pressure,buffer,ok,ghost_width=3, &
        stored_internal_energy=.true.,frame=frame,domain_lower=lower,root_spacing=spacing)
      call require(ok,'rotating pressure application')
      x=4*2**tree%nodes(p)%level-origins(1,p)+3
      y=4;z=4
      do layer=1,size(buffer,2)-x
        position=lower+(real(origins(:,p)+[x,y,z]-3,dp)-.5_dp)*spacing/real(2**tree%nodes(p)%level,dp)
        call rotating_potential_3d(frame,position,potential,ok)
        donor=kept(:nv,x,y,z);donor(5)=donor(5)+donor(1)*potential
        call characteristic_pressure_ghost_3d(species,donor,kept(nv+1,x,y,z),2,pressure(2),real(2*layer,dp), &
          expected,expected_t,ok)
        call require(ok,'independent ordinary pressure EOS')
        call require(abs(buffer(nv+2,x+layer,y,z)-(expected(5)-.5_dp*sum(expected(2:4)**2)/expected(1)))<1e-9_dp, &
          'pressure boundary retained internal energy excludes potential')
        position=lower+(real(origins(:,p)+[x+layer,y,z]-3,dp)-.5_dp)*spacing/real(2**tree%nodes(p)%level,dp)
        call rotating_potential_3d(frame,position,potential,ok)
        expected(5)=expected(5)-expected(1)*potential
        call require(maxval(abs(buffer(:nv,x+layer,y,z)-expected))<1e-9_dp, &
          'pressure ghost uses its own potential at root and fine levels')
        call require(buffer(nv+1,x+layer,y,z)==expected_t,'pressure ghost temperature oracle')
      enddo
      call require(all(buffer(:,:x,:,:)==kept(:,:x,:,:)),'untouched pressure cells retain exact bits')
      kept=buffer
      call apply_amr_tree_pressure_ghosts_3d(tree,species,p,faces,pressure,buffer,ok,ghost_width=3, &
        stored_internal_energy=.true.,frame=frame,domain_lower=lower)
      call require(.not.ok.and.all(buffer==kept),'missing rotating pressure spacing rollback')
    enddo
    faces='outflow';faces(1)='inflow'
    call fill_amr_patch_tree_ghosts_3d(tree,packed,1,faces,buffer,ok,ghost_width=3, &
      frame=frame,domain_lower=lower,root_spacing=spacing,species=species, &
      inflow_mach=[.2_dp,-1._dp,-1._dp,-1._dp,-1._dp,-1._dp], &
      inflow_temperature=[temperature,temperature,temperature,temperature,temperature,temperature], &
      limit_internal_energy=.true.,stored_internal_energy=.true.,coarse_limiter='linear',species_count=size(species))
    call require(ok,'rotating Mach inlet with retained internal energy')
    position=lower+[-.5_dp,.5_dp,.5_dp]*spacing
    call rotating_potential_3d(frame,position,potential,ok)
    call require(abs(buffer(nv+2,3,4,4)-(buffer(5,3,4,4)+buffer(1,3,4,4)*potential- &
      .5_dp*sum(buffer(2:4,3,4,4)**2)/buffer(1,3,4,4)))<1e-9_dp, &
      'rotating Mach inlet rho e excludes centrifugal potential')
    faces=saved_faces
  end subroutine rotating_pressure_checks

  subroutine rotating_inlet_checks()
    type(rotating_frame_3d) :: frame,bad_frame
    real(dp) :: lower(3),spacing(3),position(3),potential
    frame%axis=3;frame%angular_speed=400._dp
    frame%axis_point=[.0004_dp,-.0007_dp,0._dp]
    lower=[-.002_dp,.003_dp,.001_dp];spacing=[.001_dp,.0015_dp,.002_dp]
    faces=['inflow  ','outflow ','periodic','periodic','periodic','periodic']
    call fill_amr_patch_tree_ghosts_3d(tree,state,1,faces,ghost,ok,inflow_state=inlet, &
      frame=frame,domain_lower=lower,root_spacing=spacing)
    call require(ok,'rotating reservoir inlet ghost')
    position=lower+[-.5_dp,.5_dp,.5_dp]*spacing
    call rotating_potential_3d(frame,position,potential,ok)
    call require(ok,'rotating reservoir position')
    call require(all(ghost(1:4,2,3,3)==inlet(1:4,1)).and.ghost(6,2,3,3)==inlet(6,1), &
      'rotating inlet preserves density, momentum and species')
    call require(abs(ghost(5,2,3,3)-(inlet(5,1)-inlet(1,1)*potential))<1e-9_dp, &
      'rotating reservoir energy uses receiving ghost position')
    saved=ghost;bad_frame=frame;bad_frame%axis=0
    call fill_amr_patch_tree_ghosts_3d(tree,state,1,faces,ghost,ok,inflow_state=inlet, &
      frame=bad_frame,domain_lower=lower,root_spacing=spacing)
    call require(.not.ok.and.all(ghost==saved),'invalid rotating inlet rolls back ghost output')
  end subroutine rotating_inlet_checks

  subroutine linear_short_inlet_checks()
    type(amr_patch_tree_3d) :: geometry
    type(amr_tree_state_3d) :: gas
    type(amr_patch_3d) :: patch(2)
    real(dp) :: packed(7,4,4,4),reservoir(7,6),base(7),gradient(7),profile(4),expected(7),coordinate
    real(dp), allocatable :: halo(:,:,:,:)
    character(len=32) :: bc(6)
    integer :: axis,side,lower(3),upper(3),parent_shape(3),index(3),x,y,z
    logical :: success
    base=[2._dp,3._dp,4._dp,5._dp,1000._dp,.6_dp,1.4_dp]
    gradient=[1._dp,.2_dp,.3_dp,.4_dp,2._dp,.3_dp,.7_dp]
    do axis=1,3
      do side=1,2
        profile=[.5_dp,1.5_dp,3._dp,3.5_dp]
        if(side==2) profile=[.5_dp,1._dp,2.5_dp,3.5_dp]
        do z=1,4
          do y=1,4
            do x=1,4
              index=[x,y,z];packed(:,x,y,z)=base+gradient*profile(index(axis))
            enddo
          enddo
        enddo
        lower=2;upper=3;lower(axis)=1;upper(axis)=1
        if(side==2) then
          lower(axis)=4;upper(axis)=4
        endif
        call initialize_amr_patch_3d(4,4,4,lower(1),upper(1),lower(2),upper(2),lower(3),upper(3),3,patch(1),success)
        call require(success,'narrow parent geometry')
        parent_shape=(upper-lower+1)*3;lower=1;upper=1
        if(side==2) then
          lower(axis)=3;upper(axis)=3
        endif
        call initialize_amr_patch_3d(parent_shape(1),parent_shape(2),parent_shape(3), &
          lower(1),upper(1),lower(2),upper(2),lower(3),upper(3),2,patch(2),success)
        call require(success,'narrow grandchild geometry')
        call initialize_amr_patch_tree_3d([4,4,4],[1,2],patch,geometry,success)
        call require(success,'narrow donor buffer tree')
        call prolong_amr_patch_tree_3d(geometry,packed,gas,success)
        call require(success,'narrow donor buffer state')
        bc='outflow';bc(2*axis-1:2*axis)='inflow'
        reservoir=spread(base,2,6);reservoir(:,2*axis)=base+4*gradient
        call fill_amr_patch_tree_ghosts_3d(geometry,gas,3,bc,halo,success, &
          coarse_limiter='linear',species_count=2,inflow_state=reservoir,ghost_width=1)
        call require(success,'short prescribed-inlet donor stencil')
        index=2;index(modulo(axis,3)+1)=1
        coordinate=1._dp/12
        if(side==2) then
          index(axis)=3;coordinate=4._dp-1._dp/12
        endif
        expected=base+gradient*coordinate
        call require(maxval(abs(halo(:,index(1),index(2),index(3))-expected))< &
          64*epsilon(1._dp)*maxval(abs(expected)),'short buffer uses three-point rather than global four-point slope')
      enddo
    enddo
  end subroutine linear_short_inlet_checks
  subroutine linear_inlet_checks()
    type(amr_patch_tree_3d) :: geometry
    type(amr_tree_state_3d) :: gas
    type(amr_patch_3d) :: patch(1)
    real(dp) :: packed(7,4,4,4),reservoir(7,6),base(7),gradient(7),expected(7),coordinate
    real(dp), allocatable :: halo(:,:,:,:)
    character(len=32) :: bc(6)
    integer :: axis,side,cells,variant,shape_root(3),lower(3),upper(3),index(3),x,y,z
    logical :: success
    base=[2._dp,3._dp,4._dp,5._dp,1000._dp,.6_dp,1.4_dp]
    gradient=[1._dp,.2_dp,.3_dp,.4_dp,2._dp,.3_dp,.7_dp]
    do axis=1,3
      do variant=1,2
        cells=4
        if(variant==2) cells=2
        shape_root=4;shape_root(axis)=cells
        do z=1,shape_root(3)
          do y=1,shape_root(2)
            do x=1,shape_root(1)
              index=[x,y,z];coordinate=real(index(axis),dp)-.5_dp
              packed(:,x,y,z)=base+gradient*coordinate
            enddo
          enddo
        enddo
        do side=1,2
          lower=2;upper=3
          lower(axis)=1;upper(axis)=min(2,cells)
          if(side==2) then
            lower(axis)=max(1,cells-1);upper(axis)=cells
          endif
          call initialize_amr_patch_3d(shape_root(1),shape_root(2),shape_root(3), &
            lower(1),upper(1),lower(2),upper(2),lower(3),upper(3),2,patch(1),success)
          call require(success,'linear inflow child')
          call initialize_amr_patch_tree_3d(shape_root,[1],patch,geometry,success)
          call require(success,'linear inflow geometry')
          call prolong_amr_patch_tree_3d(geometry,packed(:,:shape_root(1),:shape_root(2),:shape_root(3)),gas,success)
          call require(success,'linear inflow state')
          bc='outflow';bc(2*axis-1:2*axis)='inflow'
          reservoir=spread(base,2,6)
          reservoir(:,2*axis)=base+gradient*real(cells,dp)
          call fill_amr_patch_tree_ghosts_3d(geometry,gas,2,bc,halo,success, &
            coarse_limiter='linear',species_count=2,inflow_state=reservoir,ghost_width=1)
          call require(success,'one-sided linear inflow interpolation')
          index=3;index(modulo(axis,3)+1)=1
          index(axis)=2
          coordinate=.25_dp
          if(side==2) then
            index(axis)=geometry%nodes(2)%shape(axis)+1
            coordinate=real(cells,dp)-.25_dp
          endif
          expected=base+gradient*coordinate
          call require(maxval(abs(halo(:,index(1),index(2),index(3))-expected))< &
            64*epsilon(1._dp)*maxval(abs(expected)),'affine coarse donor next to prescribed inlet')
          call require(abs(sum(halo(6:7,index(1),index(2),index(3)))-halo(1,index(1),index(2),index(3)))< &
            16*epsilon(1._dp)*abs(halo(1,index(1),index(2),index(3))),'linear species closure at inlet')
        enddo
      enddo
    enddo
  end subroutine linear_inlet_checks
  subroutine mach_inlet_checks()
    type(nasa7_species), allocatable :: thermo(:)
    type(amr_patch_tree_3d) :: geometry
    type(amr_tree_state_3d) :: gas,future
    type(amr_patch_3d) :: empty(0),patch(1)
    real(dp) :: packed(13,4,4,4),q(12),sound,ma(6),pin(6),tin(6),walls(3,6),donor(13),expected(13)
    real(dp) :: regular(12,18,18,18),regular_t(18,18,18)
    real(dp), allocatable :: halo(:,:,:,:),previous(:,:,:,:)
    character(len=32) :: bc(6)
    integer :: i,j,k,face,wall_axis
    logical :: success
    call load_h2o2_elementary_thermo(thermo,success)
    call require(success,'Mach tree thermo')
    do k=1,4
      do j=1,4
        do i=1,4
          q=0;q(1)=.3_dp+.01_dp*i;q(5)=101325+100*j;q(6)=.03_dp;q(9)=.23_dp;q(12)=.74_dp
          q(2:4)=[real(i,dp),real(j,dp),real(k,dp)]
          call reactive_primitive_to_conserved(thermo,q,packed(:12,i,j,k),packed(13,i,j,k),sound,success)
          call require(success,'Mach tree initial gas')
        end do
      end do
    end do
    call initialize_amr_patch_tree_3d([4,4,4],[integer::],empty,geometry,success)
    call require(success,'Mach tree root')
    call prolong_amr_patch_tree_3d(geometry,packed,gas,success)
    call require(success,'Mach tree packed gas')
    tin=700
    do face=1,6
      bc='outflow';bc(face)='inflow';ma=-1;ma(face)=.2_dp;pin=0;walls=0
      if(modulo(face,2)==0) pin(face)=120000
      wall_axis=modulo((face+1)/2,3)+1
      bc(2*wall_axis-1:2*wall_axis)='no_slip_wall'
      walls((face+1)/2,2*wall_axis-1)=2
      call fill_amr_patch_tree_ghosts_3d(geometry,gas,1,bc,halo,success,viscous_wall=.true.,wall_velocity=walls, &
        ghost_width=7,species=thermo,inflow_mach=ma,inflow_temperature=tin,inflow_pressure=pin)
      call require(success,'Mach tree wide root donors')
      call fill_mach_inlet_boundary_ghosts_3d(thermo,packed(:12,:,:,:),packed(13,:,:,:),bc,ma,regular,regular_t, &
        success,viscous_wall=.true.,wall_velocity=walls,ghost_width=7,inflow_temperature=tin,inflow_pressure=pin)
      call require(success,'Mach regular root oracle')
      call require(all(halo(:12,:,:,:)==regular).and.all(halo(13,:,:,:)==regular_t),'exact Mach root including walls')
    end do
    donor=packed(:,1,1,1)
    do k=1,4
      do j=1,4
        do i=1,4
          packed(:,i,j,k)=donor
        end do
      end do
    end do
    call initialize_amr_patch_3d(4,4,4,2,3,2,3,2,3,2,patch(1),success)
    call require(success,'Mach interior refinement')
    call initialize_amr_patch_tree_3d([4,4,4],[1],patch,geometry,success)
    call require(success,'Mach refined geometry')
    call prolong_amr_patch_tree_3d(geometry,packed,gas,success)
    call require(success,'Mach refined gas')
    future=gas;future%fields(1)%values(:12,:,:,:)=1.2_dp*gas%fields(1)%values(:12,:,:,:)
    bc='outflow';bc(1)='inflow';ma=-1;ma(1)=.2_dp
    call fill_amr_patch_tree_ghosts_3d(geometry,gas,2,bc,halo,success,coarse_future=future, &
      coarse_fractions=[.5_dp],coarse_limiter='mc',species_count=7,ghost_width=7, &
      species=thermo,inflow_mach=ma,inflow_temperature=tin)
    call require(success,'Mach temporally interpolated coarse donors')
    donor(:12)=.5_dp*donor(:12)+.5_dp*(1.2_dp*donor(:12))
    expected=donor
    call mach_inlet_ghost_3d(thermo,donor(:12),donor(13),1,700._dp,.2_dp,expected(:12),expected(13),success)
    call require(success,'Mach coarse endpoint oracle')
    call require(all(halo(:,1,9,9)==expected),'Mach inlet uses blended coarse pressure and composition')
    previous=halo;tin(1)=-1
    call fill_amr_patch_tree_ghosts_3d(geometry,gas,2,bc,halo,success,species=thermo,inflow_mach=ma, &
      inflow_temperature=tin,ghost_width=7)
    call require(.not.success.and.all(halo==previous),'Mach tree failure preserves output')
  end subroutine mach_inlet_checks

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
