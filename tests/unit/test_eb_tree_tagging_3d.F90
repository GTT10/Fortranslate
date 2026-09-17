program test_eb_tree_tagging_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    tag_reactive_eb_tree_temperature_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology,selected,kept
  type(eb_geometry_3d) :: finest,middle,coarse
  type(eb_reactive_tree_3d) :: tree,generated
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:)
  integer :: nv,node,child,lo(3),hi(3)
  real(dp) :: nan
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
  allocate(u(nv,4,4,4),t(4,4,4))
  config%nx=4;config%ny=4;config%nz=4
  call initialize_reactive_eb_density_sheet_3d(species,config,coarse,u,t,ok)
  call require(ok,'root state')
  call initialize_reactive_eb_tree_3d(species,topology,finest,u,t,tree,ok)
  call require(ok,'tree initialization')
  do node=2,5
    tree%nodes(node)%state=real(node,dp)*tree%nodes(node)%state
  end do
  nan=ieee_value(1._dp,ieee_quiet_nan)
  do node=1,5
    tree%nodes(node)%temperature=300._dp
    where(tree%nodes(node)%geometry%volume_fraction<=0._dp) tree%nodes(node)%temperature=nan
    do child=node+1,5
      if (topology%nodes(child)%parent/=node) cycle
      associate(p=>topology%nodes(child)%patch)
        lo=[p%coarse_i_lower,p%coarse_j_lower,p%coarse_k_lower]
        hi=[p%coarse_i_upper,p%coarse_j_upper,p%coarse_k_upper]
      end associate
      tree%nodes(node)%temperature(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=nan
    end do
  end do
  tree%nodes(4)%temperature(4,4,4)=1500._dp
  call tag_reactive_eb_tree_temperature_3d(tree,finest,[1000._dp,1400._dp],[2,2],[2,2,2],0,64,selected,ok)
  call require(ok,'fine-only hot cell and inactive NaNs')
  call require(size(selected%nodes)==3.and.selected%level_count()==3,'hot descendant selects both levels')
  call require(selected%nodes(2)%patch%coarse_i_lower==1.and. &
    selected%nodes(2)%patch%coarse_k_lower==3,'root tile identifies fine hot position')
  call require(selected%nodes(3)%patch%coarse_i_lower==3.and. &
    selected%nodes(3)%patch%coarse_k_lower==1,'second-level hot tile')
  kept=selected
  call tag_reactive_eb_tree_temperature_3d(tree,finest,[1000._dp,1600._dp],[2,2],[2,2,2],0,64,selected,ok)
  call require(ok.and.size(selected%nodes)==2,'per-level temperature thresholds')
  call tag_reactive_eb_tree_temperature_3d(tree,finest,[2000._dp,3000._dp],[2,2],[2,2,2],0,64,selected,ok)
  call require(ok.and.size(selected%nodes)==1,'cold tree coarsens to root')
  call tag_reactive_eb_tree_temperature_3d(tree,finest,[1000._dp,1400._dp],[2,2],[2,2,2],huge(1),256,selected,ok)
  call require(ok,'large buffer clipped safely')
  call initialize_reactive_eb_tree_3d(species,selected,finest,u,t,generated,ok)
  call require(ok,'tagged geometry initializes')
  do node=1,size(generated%nodes)
    call require(any(generated%nodes(node)%geometry%volume_fraction>0._dp),'solid-only tiles pruned')
  end do
  call require(size(selected%nodes)<73,'buffered solid tiles omitted')
  selected=kept
  call tag_reactive_eb_tree_temperature_3d(tree,finest,[1000._dp,1400._dp],[2,2],[2,2,2],0,1,selected,ok)
  call require(.not.ok,'capacity failure')
  call unchanged()
  tree%nodes(4)%temperature(4,4,4)=nan
  call tag_reactive_eb_tree_temperature_3d(tree,finest,[1000._dp,1400._dp],[2,2],[2,2,2],0,64,selected,ok)
  call require(.not.ok,'invalid exposed temperature')
  call unchanged()
  print *, 'test_eb_tree_tagging_3d: PASS'
contains
  subroutine unchanged()
    integer :: n
    call require(size(selected%nodes)==size(kept%nodes),'node count rollback')
    do n=1,size(selected%nodes)
      call require(selected%nodes(n)%parent==kept%nodes(n)%parent.and. &
        selected%nodes(n)%level==kept%nodes(n)%level.and. &
        all(selected%nodes(n)%shape==kept%nodes(n)%shape),'topology rollback')
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
end program test_eb_tree_tagging_3d
