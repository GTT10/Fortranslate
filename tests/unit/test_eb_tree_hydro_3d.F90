program test_eb_tree_hydro_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use simulation_config_reactive_3d_mod, only: mach_inlet_ghost_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar,reactive_conserved_to_primitive,reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,advance_reactive_eb_tree_hydro_3d, &
    fill_reactive_eb_tree_halo_3d,build_eb_tree_extended_geometry_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,middle,coarse,extended_finest
  type(eb_geometry_3d), allocatable :: extended(:),physical(:)
  type(eb_reactive_tree_3d) :: tree,saved,ctu_start
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:)
  real(dp), allocatable :: before(:),after(:),boundary(:),wall(:),ctu_before(:)
  real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),expected(:)
  real(dp) :: mach(6),tin(6),pin(6),et
  integer :: face,idx(3),inside(3),dim
  real(dp) :: coefficient,sound
  character(len=16) :: periodic_faces(6)
  character(len=32) :: mode
  real(dp) :: residual
  character(len=1), parameter :: axes(3)=['x','y','z']
  integer :: nv,node,axis,count,i,j,k
  logical :: ok
  call get_command_argument(1,mode)
  coefficient=0
  if(mode=='viscosity') coefficient=.037_dp
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  nv=reactive_nvar(size(species))
  allocate(before(nv),after(nv),boundary(nv),wall(nv),expected(nv),ctu_before(nv))
  allocate(ghost(nv,0:5,0:5,0:5),gt(0:5,0:5,0:5))
  do axis=1,3
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
  call initialize_reactive_eb_tree_3d(species,topology,finest,u,t,tree,ok)
  call require(ok,'tree initialization')
  do node=2,count
    tree%nodes(node)%state=real(node,dp)*tree%nodes(node)%state
  end do
  call average_down_reactive_eb_tree_3d(species,tree,ok)
  call require(ok,'initial average down')
  call integrate_reactive_eb_tree_3d(tree,before,ok)
  call require(ok,'initial inventory')
  saved=tree
  boundary=17._dp;wall=19._dp
  call advance_reactive_eb_tree_hydro_3d(species,'hllc',1.e-8_dp,tree,boundary,wall,ok,difmag=coefficient)
  call require(ok,'tree hydro step')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'final inventory')
  residual=maxval(abs(after-before-boundary-wall)/max(1._dp,abs(before)))
  print *, 'axis, nodes, inventory residual: ',axes(axis),count,residual
  call require(residual<3.e-14_dp,'hydro conservation including boundary and exposed wall')
  call require(any(tree%nodes(count)%state/=saved%nodes(count)%state),'nontrivial deep hydro update')
  tree=saved
  allocate(physical(count))
  do node=1,count
    physical(node)=tree%nodes(node)%geometry
  end do
  call build_axis_plane_eb_geometry_3d(64,64,64,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
    axes(axis),.31_dp,extended_finest,ok)
  call require(ok,'padded finest CTU geometry')
  call build_eb_tree_extended_geometry_3d(topology,physical,extended_finest,extended,ok)
  call require(ok,'coarsened CTU patch geometry cache')
  do node=1,count
    call require(all(extended(node)%volume_fraction(7:10,7:10,7:10)==physical(node)%volume_fraction), &
      'CTU cache retains exact physical volume fractions')
  end do
  if(axis==1.and.count==4) then
    boundary=17;wall=19
    call advance_reactive_eb_tree_hydro_3d(species,'pelec',1.e-8_dp,tree,boundary,wall,ok, &
      difmag=coefficient,extended_geometry=extended)
    call require(.not.ok.and.all(boundary==17).and.all(wall==19),'cold Riemann star rejects below mechanism range')
    call unchanged_tree()
  end if
  ! Retain discontinuous densities with a uniform 1000 K thermodynamic fixture.
  ! The original cold case above remains an explicit rejection test.
  do node=1,count
  do k=1,4
  do j=1,4
  do i=1,4
    if(tree%nodes(node)%geometry%volume_fraction(i,j,k)<=0) cycle
    call reactive_conserved_to_primitive(species,tree%nodes(node)%state(:,i,j,k), &
      tree%nodes(node)%temperature(i,j,k),expected,et,sound,ok)
    call require(ok,'CTU fixture primitive state')
    expected(5)=expected(5)*(1000._dp/et)
    call reactive_primitive_to_conserved(species,expected,tree%nodes(node)%state(:,i,j,k), &
      tree%nodes(node)%temperature(i,j,k),sound,ok)
    call require(ok,'CTU fixture within mechanism range')
  end do
  end do
  end do
  end do
  call average_down_reactive_eb_tree_3d(species,tree,ok)
  call require(ok,'CTU initial restriction')
  ctu_start=tree
  call integrate_reactive_eb_tree_3d(tree,ctu_before,ok)
  call require(ok,'CTU initial inventory')
  call advance_reactive_eb_tree_hydro_3d(species,'pelec',1.e-8_dp,tree,boundary,wall,ok, &
    difmag=coefficient,extended_geometry=extended)
  call require(ok,'CTU common tree hydro with reflux')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'CTU composite inventory')
  residual=maxval(abs(after-ctu_before-boundary-wall)/max(1._dp,abs(ctu_before)))
  print *, 'CTU axis, nodes, inventory residual: ',axes(axis),count,residual
  call require(residual<3.e-14_dp,'CTU conservative peer and coarse/fine flux synchronization')
  call require(any(tree%nodes(count)%state/=ctu_start%nodes(count)%state),'nontrivial deep CTU update')
  periodic_faces='periodic';periodic_faces(2*axis-1:2*axis)='outflow'
  tree=ctu_start
  call advance_reactive_eb_tree_hydro_3d(species,'pelec',1.e-8_dp,tree,boundary,wall,ok, &
    difmag=coefficient,extended_geometry=extended,boundary_faces=periodic_faces)
  call require(ok,'periodic CTU common tree step')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'periodic CTU inventory')
  residual=maxval(abs(after-ctu_before-boundary-wall)/max(1._dp,abs(ctu_before)))
  call require(residual<3.e-14_dp,'periodic CTU composite conservation')
  tree=saved
  call advance_reactive_eb_tree_hydro_3d(species,'hllc',1.e-8_dp,tree,boundary,wall,ok, &
    difmag=coefficient,boundary_faces=periodic_faces)
  call require(ok,'periodic PCM common tree step')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'periodic PCM inventory')
  residual=maxval(abs(after-before-boundary-wall)/max(1._dp,abs(before)))
  call require(residual<3.e-14_dp,'periodic PCM composite conservation')
  periodic_faces='no_slip_wall';periodic_faces(2*axis-1:2*axis)='outflow'
  periodic_faces(2*modulo(axis,3)+1:2*modulo(axis,3)+2)='slip_wall'
  tree=ctu_start
  call advance_reactive_eb_tree_hydro_3d(species,'pelec',1.e-8_dp,tree,boundary,wall,ok, &
    difmag=coefficient,extended_geometry=extended,boundary_faces=periodic_faces)
  call require(ok,'mixed wall CTU common tree step')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'wall CTU inventory')
  residual=maxval(abs(after-ctu_before-boundary-wall)/max(1._dp,abs(ctu_before)))
  call require(residual<3.e-14_dp,'wall CTU composite conservation')
  tree=saved
  call advance_reactive_eb_tree_hydro_3d(species,'hllc',1.e-8_dp,tree,boundary,wall,ok, &
    difmag=coefficient,boundary_faces=periodic_faces)
  call require(ok,'mixed wall PCM common tree step')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'wall PCM inventory')
  residual=maxval(abs(after-before-boundary-wall)/max(1._dp,abs(before)))
  call require(residual<3.e-14_dp,'wall PCM composite conservation')
  tree=saved;boundary=17;wall=19
  call advance_reactive_eb_tree_hydro_3d(species,'hllc',1.e-8_dp,tree,boundary,wall,ok,extended_geometry=extended)
  call require(.not.ok.and.all(boundary==17).and.all(wall==19),'CTU incompatible solver rollback')
  call unchanged_tree()
  deallocate(physical)
  tree=saved;mach=.2_dp;tin=400;pin=0
  if (count==5) pin=120000
  ghost=-17;gt=-19
  call fill_reactive_eb_tree_halo_3d(species,tree,1,ghost,gt,ok, &
    inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
  call require(ok,'Mach root halo')
  do face=1,6
    dim=(face+1)/2
    ! The plane covers the low face in its normal direction.
    if (dim==axis.and.modulo(face,2)==1) cycle
    idx=3;idx(dim)=merge(0,5,modulo(face,2)==1);inside=max(1,min(4,idx))
    if (pin(face)>0) then
      call mach_inlet_ghost_3d(species,tree%nodes(1)%state(:,inside(1),inside(2),inside(3)), &
        tree%nodes(1)%temperature(inside(1),inside(2),inside(3)),face,tin(face),mach(face),expected,et,ok, &
        inlet_pressure=pin(face))
    else
      call mach_inlet_ghost_3d(species,tree%nodes(1)%state(:,inside(1),inside(2),inside(3)), &
        tree%nodes(1)%temperature(inside(1),inside(2),inside(3)),face,tin(face),mach(face),expected,et,ok)
    end if
    call require(ok.and.all(ghost(:,idx(1),idx(2),idx(3))==expected).and.gt(idx(1),idx(2),idx(3))==et, &
      'physical Mach face orientation and pressure')
  end do
  call advance_reactive_eb_tree_hydro_3d(species,'hllc',1.e-8_dp,tree,boundary,wall,ok,mach,tin,pin,difmag=coefficient)
  call require(ok,'Mach EB tree hydro')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'Mach final inventory')
  residual=maxval(abs(after-before-boundary-wall)/max(1._dp,abs(before)))
  call require(residual<3.e-14_dp,'Mach EB conservation with root boundary replacement')
  tree=saved;boundary=17;wall=19;tin(1)=ieee_value(0._dp,ieee_quiet_nan)
  call advance_reactive_eb_tree_hydro_3d(species,'hllc',1.e-8_dp,tree,boundary,wall,ok,mach,tin,pin,difmag=coefficient)
  call require(.not.ok,'invalid Mach temperature rejected even on covered face')
  call unchanged_tree()
  call require(all(boundary==17).and.all(wall==19),'Mach diagnostics rollback')
  ghost=-17;gt=-19
  call fill_reactive_eb_tree_halo_3d(species,tree,1,ghost,gt,ok,inflow_mach=mach)
  call require(.not.ok.and.all(ghost==-17).and.all(gt==-19),'missing inlet temperature halo rollback')
  tree=saved;boundary=17._dp;wall=19._dp
  call advance_reactive_eb_tree_hydro_3d(species,'invalid',1.e-8_dp,tree,boundary,wall,ok,difmag=coefficient)
  call require(.not.ok,'invalid solver rejected')
  call unchanged_tree()
  call require(all(boundary==17._dp).and.all(wall==19._dp),'diagnostic rollback')
  tree%nodes(count)%state(1,4,4,4)=-1._dp;saved=tree
  call advance_reactive_eb_tree_hydro_3d(species,'hllc',1.e-8_dp,tree,boundary,wall,ok,difmag=coefficient)
  call require(.not.ok,'invalid deep state rejected')
  call unchanged_tree()
  call require(all(boundary==17._dp).and.all(wall==19._dp),'deep failure diagnostic rollback')
  deallocate(u,t)
  end do
  end do
  print *, 'test_eb_tree_hydro_3d: PASS'
contains
  subroutine unchanged_tree()
    integer :: n
    do n=1,size(tree%nodes)
      call require(all(tree%nodes(n)%state==saved%nodes(n)%state).and. &
        all(tree%nodes(n)%temperature==saved%nodes(n)%temperature),'all state and temperature rollback')
    end do
  end subroutine unchanged_tree
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine require
end program test_eb_tree_hydro_3d
