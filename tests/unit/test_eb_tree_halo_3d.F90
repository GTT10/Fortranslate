program test_eb_tree_halo_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    fill_reactive_eb_tree_halo_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,middle,coarse
  type(eb_reactive_tree_3d) :: tree,saved,start,finish
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:),ghost(:,:,:,:),gt(:,:,:)
  real(dp), allocatable :: inlet(:,:),wide(:,:,:,:),wide_t(:,:,:),vf(:,:,:),target_vf(:,:,:)
  real(dp) :: inlet_t(6)
  character(len=16) :: faces(6)
  integer :: i,j,k
  real(dp) :: expected_fraction
  integer :: nv,node,endpoint
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  nv=reactive_nvar(size(species))
  call initialize_amr_patch_3d(4,4,4,1,2,1,2,2,3,2,boxes(1),ok)
  call require(ok,'child left')
  call initialize_amr_patch_3d(4,4,4,3,4,1,2,2,3,2,boxes(2),ok)
  call require(ok,'child right')
  call initialize_amr_patch_3d(4,4,4,3,4,3,4,2,3,2,boxes(3),ok)
  call require(ok,'grandchild left')
  call initialize_amr_patch_3d(4,4,4,1,2,3,4,2,3,2,boxes(4),ok)
  call require(ok,'grandchild right')
  call initialize_amr_patch_tree_3d([4,4,4],[1,1,2,3],boxes,topology,ok)
  call require(ok,'topology')
  call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'z',.31_dp,finest,ok)
  call require(ok,'finest geometry')
  call coarsen_eb_geometry_3d(finest,2,middle,ok)
  call require(ok,'middle geometry')
  call coarsen_eb_geometry_3d(middle,2,coarse,ok)
  call require(ok,'coarse geometry')
  allocate(u(nv,4,4,4),t(4,4,4),ghost(nv,0:5,0:5,0:5),gt(0:5,0:5,0:5))
  config%nx=4;config%ny=4;config%nz=4
  call initialize_reactive_eb_density_sheet_3d(species,config,coarse,u,t,ok)
  call require(ok,'root state')
  call initialize_reactive_eb_tree_3d(species,topology,finest,u,t,tree,ok)
  call require(ok,'tree initialization')
  do node=2,5
    tree%nodes(node)%state=real(node,dp)*tree%nodes(node)%state
  end do
  saved=tree
  call fill_reactive_eb_tree_halo_3d(species,tree,4,ghost,gt,ok)
  call require(ok,'nested halo')
  call require(all(ghost(:,5,2,2)==tree%nodes(5)%state(:,1,2,2)),'same-level peer under different parent')
  call require(all(ghost(:,2,5,2)==tree%nodes(1)%state(:,2,3,2)),'grandparent donor beyond parent edge')
  call require(all(ghost(:,2,2,0)==tree%nodes(2)%state(:,3,3,1)),'immediate parent donor')
  call require(gt(5,2,2)>0._dp.and.gt(2,5,2)>0._dp.and.gt(2,2,0)>0._dp,'donor EOS temperatures')
  call require(all(ghost(:,0,0,0)==0._dp),'unused corner remains empty')
  call fill_reactive_eb_tree_halo_3d(species,tree,2,ghost,gt,ok)
  call require(ok,'boundary-touching parent halo')
  call require(all(ghost(:,0,2,2)==tree%nodes(2)%state(:,1,2,2)),'physical outflow uses adjacent target cell')
  call fill_reactive_eb_tree_halo_3d(species,tree,1,ghost,gt,ok)
  call require(ok,'root halo')
  print *, 'root donor relative roundoff: ',maxval(abs(ghost(:,0,1,3)-2._dp*u(:,1,1,3))/ &
    max(1._dp,abs(2._dp*u(:,1,1,3))))
  ! Eight-child restriction sums weighted terms; its analytic uniform value
  ! is compared with an eight-term roundoff bound, unlike exact donor copies.
  call require(all(abs(ghost(:,0,1,3)-2._dp*u(:,1,1,3))<= &
    8._dp*epsilon(1._dp)*max(1._dp,abs(2._dp*u(:,1,1,3)))),'root uses restricted child inventory')
  call require(all(ghost(:,2,2,0)==0._dp).and.gt(2,2,0)==0._dp,'solid face remains empty')
  do node=1,5
    call require(all(tree%nodes(node)%state==saved%nodes(node)%state).and. &
      all(tree%nodes(node)%temperature==saved%nodes(node)%temperature),'halo leaves source snapshot unchanged')
  end do
  start=tree;finish=tree
  do node=1,5
    finish%nodes(node)%state=2._dp*start%nodes(node)%state
  end do
  do endpoint=0,2
    expected_fraction=.5_dp*real(endpoint,dp)
    call fill_reactive_eb_tree_halo_3d(species,tree,4,ghost,gt,ok,start,finish, &
      [expected_fraction,expected_fraction])
    call require(ok,'time-interpolated halo')
    call require(all(abs(ghost(:,2,5,2)-(1._dp+expected_fraction)*start%nodes(1)%state(:,2,3,2))<= &
      4._dp*epsilon(1._dp)*max(1._dp,abs(ghost(:,2,5,2)))),'root time donor')
    call require(all(abs(ghost(:,2,2,0)-(1._dp+expected_fraction)*start%nodes(2)%state(:,3,3,1))<= &
      4._dp*epsilon(1._dp)*max(1._dp,abs(ghost(:,2,2,0)))),'parent time donor')
    call require(all(ghost(:,5,2,2)==tree%nodes(5)%state(:,1,2,2)),'peer stays at fine time')
  end do
  call fill_reactive_eb_tree_halo_3d(species,tree,4,ghost,gt,ok,start,finish,[0._dp,1._dp])
  call require(ok,'independent ancestor intervals')
  call require(all(ghost(:,2,5,2)==start%nodes(1)%state(:,2,3,2)),'root left endpoint exact')
  call require(all(ghost(:,2,2,0)==finish%nodes(2)%state(:,3,3,1)),'parent right endpoint exact')
  ghost=17._dp;gt=19._dp
  call fill_reactive_eb_tree_halo_3d(species,tree,4,ghost,gt,ok,start,finish,[0._dp,1.1_dp])
  call require(.not.ok.and.all(ghost==17._dp).and.all(gt==19._dp),'invalid fraction rollback')
  call fill_reactive_eb_tree_halo_3d(species,tree,4,ghost,gt,ok,ancestor_start=start)
  call require(.not.ok.and.all(ghost==17._dp).and.all(gt==19._dp),'incomplete time context rollback')
  finish%nodes(2)%state(1,3,3,1)=-1._dp
  call fill_reactive_eb_tree_halo_3d(species,tree,4,ghost,gt,ok,start,finish,[0._dp,1._dp])
  call require(.not.ok.and.all(ghost==17._dp).and.all(gt==19._dp),'invalid time donor rollback')
  tree%nodes(5)%state(1,1,2,2)=-1._dp
  ghost=17._dp;gt=19._dp
  call fill_reactive_eb_tree_halo_3d(species,tree,4,ghost,gt,ok)
  call require(.not.ok.and.all(ghost==17._dp).and.all(gt==19._dp),'invalid peer rollback')
  tree=saved
  call fill_reactive_eb_tree_halo_3d(species,tree,0,ghost,gt,ok)
  call require(.not.ok.and.all(ghost==17._dp).and.all(gt==19._dp),'invalid target rollback')
  allocate(inlet(nv,6))
  do i=1,6
    inlet(:,i)=real(i+1,dp)*u(:,1,1,3)
  end do
  inlet_t=t(1,1,3);faces='outflow';faces(1)='inflow';faces(3)='inflow'
  call fill_reactive_eb_tree_halo_3d(species,tree,2,ghost,gt,ok, &
    boundary_faces=faces,inflow_state=inlet,inflow_temperature=inlet_t)
  call require(ok.and.all(ghost(:,0,2,2)==inlet(:,1)),'prescribed inlet at boundary child')
  call require(all(ghost(:,2,0,2)==inlet(:,3)),'second prescribed face')
  call fill_reactive_eb_tree_halo_3d(species,tree,4,ghost,gt,ok, &
    boundary_faces=faces,inflow_state=inlet,inflow_temperature=inlet_t)
  call require(ok.and.all(ghost(:,5,2,2)==tree%nodes(5)%state(:,1,2,2)), &
    'prescribed outer faces leave interior peer donor unchanged')
  allocate(wide(nv,0:15,0:15,0:15),wide_t(0:15,0:15,0:15),vf(0:15,0:15,0:15),target_vf(16,16,16))
  do k=0,15
  do j=0,15
  do i=0,15
    target_vf(i+1,j+1,k+1)=coarse%volume_fraction(max(1,min(4,i-5)), &
      max(1,min(4,j-5)),max(1,min(4,k-5)))
  end do
  end do
  end do
  wide=17;wide_t=19;vf=23
  call fill_reactive_eb_tree_halo_3d(species,tree,1,wide,wide_t,ok, &
    ghost_volume_fraction=vf,halo_width=6,target_volume_fraction=target_vf, &
    boundary_faces=faces,inflow_state=inlet,inflow_temperature=inlet_t)
  call require(ok.and.all(wide(:,0,6,9)==inlet(:,1)),'six-layer prescribed inlet')
  call require(all(wide(:,0,0,9)==inlet(:,1)),'corner uses first inflow axis')
  call require(all(wide(:,6,0,9)==inlet(:,3)),'wide second face inlet')
  call require(all(vf==target_vf),'wide inlet volume fractions')
  inlet(1,1)=-1;wide=17;wide_t=19;vf=23
  call fill_reactive_eb_tree_halo_3d(species,tree,1,wide,wide_t,ok, &
    ghost_volume_fraction=vf,halo_width=6,target_volume_fraction=target_vf, &
    boundary_faces=faces,inflow_state=inlet,inflow_temperature=inlet_t)
  call require(.not.ok.and.all(wide==17).and.all(wide_t==19).and.all(vf==23),'invalid wide inlet rollback')
  ghost=17;gt=19
  call fill_reactive_eb_tree_halo_3d(species,tree,2,ghost,gt,ok,boundary_faces=faces)
  call require(.not.ok.and.all(ghost==17).and.all(gt==19),'missing inlet data rollback')
  call fill_reactive_eb_tree_halo_3d(species,tree,2,ghost,gt,ok, &
    inflow_mach=[.2_dp,-1._dp,-1._dp,-1._dp,-1._dp,-1._dp],inflow_temperature=inlet_t)
  call require(ok.and.ghost(2,0,2,2)>0,'legacy Mach inlet remains available')
  ghost=17;gt=19
  call fill_reactive_eb_tree_halo_3d(species,tree,2,ghost,gt,ok, &
    boundary_faces=faces,inflow_state=inlet,inflow_temperature=inlet_t,inflow_mach=[(0._dp,i=1,6)])
  call require(.not.ok.and.all(ghost==17).and.all(gt==19),'conflicting inlet modes reject')
  print *, 'test_eb_tree_halo_3d: PASS'
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine require
end program test_eb_tree_halo_3d
