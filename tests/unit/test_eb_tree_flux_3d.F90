program test_eb_tree_flux_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,synchronize_reactive_eb_tree_flux_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,middle,coarse
  type(eb_reactive_tree_3d) :: tree,saved
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:)
  type(eb_integrated_flux_3d), allocatable :: flux(:),raw(:),bad(:)
  real(dp), allocatable :: before(:),after(:),boundary(:),basis(:)
  real(dp) :: volume,residual
  character(len=1), parameter :: axes(3)=['x','y','z']
  integer :: nv,node,axis,count,i,j,k
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  nv=reactive_nvar(size(species))
  allocate(before(nv),after(nv),boundary(nv),basis(nv))
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
  allocate(u(nv,4,4,4),t(4,4,4),flux(count))
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
  ! Manufactured time-integrated transfers, deliberately unequal between nodes.
  ! A single physical conserved vector keeps species closure and positive EOS.
  basis=u(:,4,4,4)
  do node=1,count
    allocate(flux(node)%x(nv,0:4,4,4),flux(node)%y(nv,4,0:4,4),flux(node)%z(nv,4,4,0:4))
    do k=1,4
    do j=1,4
    do i=0,4
      flux(node)%x(:,i,j,k)=basis*1.e-6_dp*real(node*(1+i+2*j+3*k),dp)
      flux(node)%y(:,j,i,k)=basis*2.e-6_dp*real(node*(2+2*i+j+k),dp)
      flux(node)%z(:,j,k,i)=basis*3.e-6_dp*real(node*(3+i+j+2*k),dp)
    end do
    end do
    end do
    associate(g=>tree%nodes(node)%geometry,s=>tree%nodes(node)%state,f=>flux(node))
    do k=1,4
    do j=1,4
    do i=1,4
      if (g%volume_fraction(i,j,k)<=0._dp) cycle
      volume=g%volume_fraction(i,j,k)*g%dx*g%dy*g%dz
      s(:,i,j,k)=s(:,i,j,k)+( &
        (g%x_face_fraction(i-1,j,k)*f%x(:,i-1,j,k)-g%x_face_fraction(i,j,k)*f%x(:,i,j,k))*g%dy*g%dz+ &
        (g%y_face_fraction(i,j-1,k)*f%y(:,i,j-1,k)-g%y_face_fraction(i,j,k)*f%y(:,i,j,k))*g%dx*g%dz+ &
        (g%z_face_fraction(i,j,k-1)*f%z(:,i,j,k-1)-g%z_face_fraction(i,j,k)*f%z(:,i,j,k))*g%dx*g%dy)/volume
    end do
    end do
    end do
    end associate
  end do
  saved=tree;raw=flux
  call synchronize_reactive_eb_tree_flux_3d(species,tree,flux,ok)
  call require(ok,'flux synchronization')
  call require(any(flux(1)%x/=raw(1)%x),'nontrivial coarse replacement')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'final inventory')
  boundary=0._dp
  associate(g=>tree%nodes(1)%geometry,f=>flux(1))
  do j=1,4
  do i=1,4
    boundary=boundary+(g%x_face_fraction(0,i,j)*f%x(:,0,i,j)-g%x_face_fraction(4,i,j)*f%x(:,4,i,j))*g%dy*g%dz &
      +(g%y_face_fraction(i,0,j)*f%y(:,i,0,j)-g%y_face_fraction(i,4,j)*f%y(:,i,4,j))*g%dx*g%dz &
      +(g%z_face_fraction(i,j,0)*f%z(:,i,j,0)-g%z_face_fraction(i,j,4)*f%z(:,i,j,4))*g%dx*g%dy
  end do
  end do
  end associate
  residual=maxval(abs(after-before-boundary)/max(1._dp,abs(before)))
  print *, 'axis, nodes, inventory residual: ',axes(axis),count,residual
  call require(residual<3.e-14_dp,'composite conservation including physical boundaries')
  ! A malformed deep flux must preserve both inputs, including earlier nodes.
  tree=saved;bad=raw
  bad(count)%x(1,0,1,1)=ieee_value(1._dp,ieee_quiet_nan)
  flux=bad
  call synchronize_reactive_eb_tree_flux_3d(species,tree,flux,ok)
  call require(.not.ok,'nonfinite flux rejected')
  call unchanged_tree()
  do node=1,count-1
    call require(all(flux(node)%x==bad(node)%x),'earlier flux rollback')
  end do
  ! Finite but enormous density transfer fails EOS after private changes.
  bad=raw;bad(count)%x(1,4,:,:)=1.e20_dp;flux=bad
  call synchronize_reactive_eb_tree_flux_3d(species,tree,flux,ok)
  call require(.not.ok,'nonphysical correction rejected')
  call unchanged_tree()
  do node=1,count
    call require(all(flux(node)%x==bad(node)%x).and.all(flux(node)%y==bad(node)%y).and. &
      all(flux(node)%z==bad(node)%z),'all flux rollback')
  end do
  deallocate(u,t,flux)
  end do
  end do
  print *, 'test_eb_tree_flux_3d: PASS'
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
end program test_eb_tree_flux_3d
