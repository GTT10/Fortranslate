program test_eb_tree_remap_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,remap_reactive_eb_tree_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(rotating_frame_3d), allocatable :: frame
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology,moved,root_only
  type(eb_geometry_3d) :: finest,middle,coarse
  type(eb_reactive_tree_3d) :: tree,saved,changed
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:)
  real(dp), allocatable :: before(:),after(:)
  real(dp) :: residual,position(3),potential,worst
  real(dp), allocatable :: primitive(:)
  character(len=1), parameter :: axes(3)=['x','y','z']
  integer :: nv,node,axis,count,i,j,k,rotation
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  nv=reactive_nvar(size(species))
  allocate(before(nv),after(nv),primitive(reactive_nprim(size(species))))
  worst=0
  do rotation=0,1
  if(rotation==1) then
    allocate(frame)
    frame%angular_speed=100;frame%axis_point=[10._dp,10._dp,10._dp]
  end if
  do axis=1,3
  if(allocated(frame)) frame%axis=axis
  do count=4,5
  call initialize_amr_patch_3d(4,4,4,1,2,1,2,2,3,2,boxes(1),ok)
  call require(ok,'child left')
  call initialize_amr_patch_3d(4,4,4,3,4,1,2,2,3,2,boxes(2),ok)
  call require(ok,'child right')
  call initialize_amr_patch_3d(4,4,4,3,4,3,4,2,3,2,boxes(3),ok)
  call require(ok,'grandchild left')
  call initialize_amr_patch_3d(4,4,4,1,2,3,4,2,3,2,boxes(4),ok)
  call require(ok,'grandchild right')
  if (count==5) then
    call initialize_amr_patch_tree_3d([4,4,4],[1,1,2,3],boxes,topology,ok)
  else
    call initialize_amr_patch_tree_3d([4,4,4],[1,1,2],boxes(:3),topology,ok)
  end if
  call require(ok,'topology')
  call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,axes(axis),.31_dp,finest,ok)
  call require(ok,'finest geometry')
  call coarsen_eb_geometry_3d(finest,2,middle,ok)
  call require(ok,'middle geometry')
  call coarsen_eb_geometry_3d(middle,2,coarse,ok)
  call require(ok,'coarse geometry')
  allocate(u(nv,4,4,4),t(4,4,4))
  config%nx=4;config%ny=4;config%nz=4
  call initialize_reactive_eb_density_sheet_3d(species,config,coarse,u,t,ok)
  call require(ok,'root state')
  if(allocated(frame)) then
    do k=1,4
    do j=1,4
    do i=1,4
      if(coarse%volume_fraction(i,j,k)<=0) cycle
      position=(real([i,j,k],dp)-.5_dp)/4
      call rotating_potential_3d(frame,position,potential,ok);call require(ok,'root potential')
      u(5,i,j,k)=u(5,i,j,k)-u(1,i,j,k)*potential
    end do
    end do
    end do
    call require(minval(u(5,:,:,:),mask=coarse%volume_fraction>0)<0,'negative rotating energy fixture')
  end if
  call initialize_reactive_eb_tree_3d(species,topology,finest,u,t,tree,ok,frame)
  call require(ok,'tree initialization')
  do node=2,count
    tree%nodes(node)%state=real(node,dp)*tree%nodes(node)%state
  end do
  do i=1,4
    tree%nodes(4)%state(:,i,:,:)=tree%nodes(4)%state(:,i,:,:)*(1._dp+.01_dp*real(i,dp))
  end do
  call average_down_reactive_eb_tree_3d(species,tree,ok,frame)
  call require(ok,'initial average down')
  call integrate_reactive_eb_tree_3d(tree,before,ok)
  call require(ok,'initial inventory')
  saved=tree
  call initialize_amr_patch_3d(4,4,4,2,3,3,4,2,3,2,boxes(3),ok)
  call require(ok,'shifted grandchild')
  call initialize_amr_patch_tree_3d([4,4,4],[1,1,2],boxes(:3),moved,ok)
  call require(ok,'new topology')
  call remap_reactive_eb_tree_3d(species,moved,finest,tree,ok,frame=frame)
  call require(ok,'move grandchild and remove other refinement')
  call check_inventory()
  call require(all(tree%nodes(4)%state(:,3:4,:,:)==saved%nodes(4)%state(:,1:2,:,:)), &
    'exact same-level overlap states')
  call require(all(tree%nodes(4)%temperature(3:4,:,:)==saved%nodes(4)%temperature(1:2,:,:)), &
    'exact same-level overlap temperatures')
  ! Root-only target coarsens all old descendants; the input finest geometry
  ! still describes the deeper of the source and target trees.
  call initialize_amr_patch_tree_3d([4,4,4],[integer::],boxes(:0),root_only,ok)
  call require(ok,'root-only topology')
  call remap_reactive_eb_tree_3d(species,root_only,finest,tree,ok,frame=frame)
  call require(ok.and.size(tree%nodes)==1,'complete coarsening')
  call check_inventory()
  call remap_reactive_eb_tree_3d(species,topology,finest,tree,ok,frame=frame)
  call require(ok,'refine root-only state to three levels')
  call check_inventory()
  call check_source_remap()
  ! Source EOS and immutable geometry failures preserve the original tree.
  tree=saved
  tree%nodes(count)%state(1,4,4,4)=-1._dp
  changed=tree
  call remap_reactive_eb_tree_3d(species,moved,finest,tree,ok,frame=frame)
  call require(.not.ok,'invalid source rejected')
  call unchanged()
  tree=saved;changed=tree
  call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
    axes(axis),.32_dp,finest,ok)
  call require(ok,'changed geometry fixture')
  call remap_reactive_eb_tree_3d(species,moved,finest,tree,ok,frame=frame)
  call require(.not.ok,'geometry change rejected')
  call unchanged()
  deallocate(u,t)
  end do
  end do
  end do
  print *, 'test_eb_tree_remap_3d: PASS; worst conservation residual:',worst
contains
  subroutine check_source_remap()
    type(amr_tree_state_3d) :: source,previous,heat,previous_heat
    type(eb_reactive_tree_3d) :: candidate,old
    real(dp) :: first(nv),last(nv),first_heat(nv),last_heat(nv)
    integer :: p,c
    candidate=saved
    allocate(source%fields(count),heat%fields(count))
    do p=1,count
      allocate(source%fields(p)%values(nv-4,4,4,4))
      allocate(heat%fields(p)%values(1,4,4,4),source=real((-1)**p*p,dp))
      do c=1,nv-4
        source%fields(p)%values(c,:,:,:)=real((-1)**c*c*p,dp)
      end do
    end do
    call source_inventory(candidate,source,first)
    call source_inventory(candidate,heat,first_heat)
    call remap_reactive_eb_tree_3d(species,moved,finest,candidate,ok,source,frame,heat)
    call require(ok,'signed source moved refinement')
    call source_inventory(candidate,source,last)
    call source_inventory(candidate,heat,last_heat)
    call require(all(abs(last_heat-first_heat)<=3.e-14_dp*max(1._dp,abs(first_heat))), &
      'signed stage heat remap conservation')
    call require(all(abs(last-first)<=3.e-14_dp*max(1._dp,abs(first))),'signed source movement conservation')
    call remap_reactive_eb_tree_3d(species,root_only,finest,candidate,ok,source,frame,heat)
    call require(ok,'signed source complete coarsening')
    call remap_reactive_eb_tree_3d(species,topology,finest,candidate,ok,source,frame,heat)
    call require(ok,'signed source refinement')
    call source_inventory(candidate,source,last)
    call source_inventory(candidate,heat,last_heat)
    call require(all(abs(last_heat-first_heat)<=3.e-14_dp*max(1._dp,abs(first_heat))), &
      'signed stage heat remap conservation')
    call require(all(abs(last-first)<=3.e-14_dp*max(1._dp,abs(first))),'signed source coarsen/refine conservation')
    candidate%nodes(count)%state(1,4,4,4)=-1;old=candidate;previous=source;previous_heat=heat
    call remap_reactive_eb_tree_3d(species,moved,finest,candidate,ok,source,frame,heat)
    call require(.not.ok,'source remap bad gas rejects')
    do p=1,count
      call require(all(candidate%nodes(p)%state==old%nodes(p)%state).and. &
        all(source%fields(p)%values==previous%fields(p)%values),'gas and signed source rollback')
      call require(all(heat%fields(p)%values==previous_heat%fields(p)%values),'stage heat rollback')
    end do
  end subroutine check_source_remap
  subroutine source_inventory(layout,source,inventory)
    type(eb_reactive_tree_3d), intent(in) :: layout
    type(amr_tree_state_3d), intent(in) :: source
    real(dp), intent(inout) :: inventory(:)
    type(eb_reactive_tree_3d) :: temp
    integer :: p
    temp=layout
    do p=1,size(layout%nodes)
      temp%nodes(p)%state=0
      temp%nodes(p)%state(:size(source%fields(p)%values,1),:,:,:)=source%fields(p)%values
    end do
    call integrate_reactive_eb_tree_3d(temp,inventory,ok)
    call require(ok,'signed source composite inventory')
  end subroutine source_inventory

  subroutine check_inventory()
    integer :: p,x,y,z
    real(dp) :: recovered,speed
    call integrate_reactive_eb_tree_3d(tree,after,ok)
    call require(ok,'remapped inventory')
    residual=maxval(abs(after-before)/max(1._dp,abs(before)))
    print *, 'axis, source nodes, residual: ',axes(axis),count,residual
    call require(residual<3.e-14_dp,'remap composite conservation')
    worst=max(worst,residual)
    if(allocated(frame)) then
      do p=1,size(tree%nodes)
      do z=1,4
      do y=1,4
      do x=1,4
        if(tree%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
        associate(g=>tree%nodes(p)%geometry)
        position=[g%x_lower,g%y_lower,g%z_lower]+(real([x,y,z],dp)-.5_dp)*[g%dx,g%dy,g%dz]
        end associate
        call rotating_potential_3d(frame,position,potential,ok);call require(ok,'receiving potential')
        call reactive_conserved_to_primitive(species,tree%nodes(p)%state(:,x,y,z), &
          tree%nodes(p)%temperature(x,y,z),primitive,recovered,speed,ok,potential)
        call require(ok,'receiving EOS')
        call require(abs(recovered/tree%nodes(p)%temperature(x,y,z)-1)<1e-12_dp,'receiving temperature')
      end do
      end do
      end do
      end do
    end if
  end subroutine check_inventory
  subroutine unchanged()
    integer :: n
    call require(size(tree%nodes)==size(changed%nodes),'node count rollback')
    do n=1,size(tree%nodes)
      call require(tree%topology%nodes(n)%parent==changed%topology%nodes(n)%parent,'parent rollback')
      call require(all(tree%nodes(n)%state==changed%nodes(n)%state).and. &
        all(tree%nodes(n)%temperature==changed%nodes(n)%temperature),'all state and temperature rollback')
    end do
  end subroutine unchanged
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine require
end program test_eb_tree_remap_3d
