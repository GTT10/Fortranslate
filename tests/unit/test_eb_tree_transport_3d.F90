program test_eb_tree_transport_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,advance_reactive_eb_tree_transport_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,middle,coarse
  type(eb_reactive_tree_3d) :: tree,saved
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:)
  real(dp), allocatable :: before(:),after(:),boundary(:),primitive(:)
  real(dp) :: residual,theta,speed
  character(len=1), parameter :: axes(3)=['x','y','z']
  integer :: nv,node,axis,count,i,j,k
  character(len=16) :: periodic_faces(6)
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  nv=reactive_nvar(size(species))
  allocate(before(nv),after(nv),boundary(nv),primitive(reactive_nprim(size(species))))
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
    tree%nodes(node)%state(5,:,:,:)=1.1_dp*tree%nodes(node)%state(5,:,:,:)
    do k=1,4
    do j=1,4
    do i=1,4
      if (tree%nodes(node)%geometry%volume_fraction(i,j,k)<=0._dp) cycle
      call reactive_conserved_to_primitive(species,tree%nodes(node)%state(:,i,j,k),t(i,j,k),primitive, &
        tree%nodes(node)%temperature(i,j,k),speed,ok)
      call require(ok,'temperature perturbation EOS')
    end do
    end do
    end do
  end do
  call average_down_reactive_eb_tree_3d(species,tree,ok)
  call require(ok,'initial average down')
  call integrate_reactive_eb_tree_3d(tree,before,ok)
  call require(ok,'initial inventory')
  saved=tree
  boundary=17._dp;theta=19._dp
  call advance_reactive_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
    tree,boundary,theta,ok)
  call require(ok,'tree transport step')
  call require(all(boundary==0._dp),'closed physical diffusion faces')
  call require(theta>=0._dp.and.theta<=1._dp,'limiter range')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'final inventory')
  residual=maxval(abs(after-before-boundary)/max(1._dp,abs(before)))
  print *, 'axis, nodes, inventory residual: ',axes(axis),count,residual
  call require(residual<3.e-14_dp,'transport conservation')
  call require(any(tree%nodes(count)%state/=saved%nodes(count)%state),'nontrivial deep transport update')
  periodic_faces='no_slip_wall'
  tree=saved
  call advance_reactive_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
    tree,boundary,theta,ok,boundary_faces=periodic_faces)
  call require(ok,'wall transport update')
  call require(boundary(1)==0.and.boundary(5)==0.and.all(boundary(6:)==0),'impermeable adiabatic wall flux')
  call require(theta>=0.and.theta<=1,'periodic limiter range')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'wall inventory')
  residual=maxval(abs(after-before-boundary)/max(1._dp,abs(before)))
  print *, 'wall axis, nodes, residual: ',axes(axis),count,residual
  call require(residual<3.e-14_dp,'wall composite transport conservation')
  periodic_faces='periodic';periodic_faces(2*axis-1:2*axis)='outflow'
  tree=saved
  call advance_reactive_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
    tree,boundary,theta,ok,boundary_faces=periodic_faces)
  call require(ok,'periodic transport update')
  call require(theta>=0.and.theta<=1,'periodic limiter range')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'periodic inventory')
  residual=maxval(abs(after-before-boundary)/max(1._dp,abs(before)))
  print *, 'periodic axis, nodes, residual: ',axes(axis),count,residual
  call require(residual<3.e-14_dp,'periodic composite transport conservation')
  tree=saved;boundary=17;theta=19
  periodic_faces(2*mod(axis,3)+1)='outflow'
  call advance_reactive_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
    tree,boundary,theta,ok,boundary_faces=periodic_faces)
  call require(.not.ok.and.all(boundary==17).and.theta==19,'unpaired periodic transport rollback')
  call unchanged_tree()
  tree=saved;boundary=17._dp;theta=19._dp
  call advance_reactive_eb_tree_transport_3d(species,transport,-1.e-8_dp,.true.,.true.,.true.,.true., &
    tree,boundary,theta,ok)
  call require(.not.ok,'invalid interval rejected')
  call unchanged_tree()
  call require(all(boundary==17._dp).and.theta==19._dp,'diagnostic rollback')
  tree%nodes(count)%state(1,4,4,4)=-1._dp;saved=tree
  call advance_reactive_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
    tree,boundary,theta,ok)
  call require(.not.ok,'invalid deep state rejected')
  call unchanged_tree()
  call require(all(boundary==17._dp).and.theta==19._dp,'deep failure diagnostic rollback')
  deallocate(u,t)
  end do
  end do
  print *, 'test_eb_tree_transport_3d: PASS'
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
end program test_eb_tree_transport_3d
