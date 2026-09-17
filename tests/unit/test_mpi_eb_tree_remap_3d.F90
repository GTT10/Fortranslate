program test_mpi_eb_tree_remap_3d
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,remap_reactive_eb_tree_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use mpi_eb_tree_remap_3d_mod, only: remap_mpi_eb_tree_3d
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d
  use eb_reactive_two_level_3d_mod, only: matching_eb_metrics_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology,moved,root_only,current_tree
  type(eb_geometry_3d) :: finest,middle,coarse
  type(eb_reactive_tree_3d) :: tree,saved
  type(reactive_eb_3d_config) :: config
  type(amr_tree_state_3d) :: local,lt,original,original_t
  type(eb_geometry_3d), allocatable :: geometry(:),original_geometry(:)
  integer, allocatable :: owners(:),next_owners(:)
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:)
  real(dp), allocatable :: before(:),after(:)
  real(dp) :: residual
  character(len=1), parameter :: axes(3)=['x','y','z']
  integer :: nv,node,axis,count,i,rank,nranks,status,stage
  logical :: ok
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  nv=reactive_nvar(size(species))
  allocate(before(nv),after(nv))
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
  do i=1,4
    tree%nodes(4)%state(:,i,:,:)=tree%nodes(4)%state(:,i,:,:)*(1._dp+.01_dp*real(i,dp))
  end do
  call average_down_reactive_eb_tree_3d(species,tree,ok)
  call require(ok,'initial average down')
  call integrate_reactive_eb_tree_3d(tree,before,ok)
  call require(ok,'initial inventory')
  saved=tree;current_tree=topology;stage=0
  allocate(geometry(count),owners(count),local%fields(count),lt%fields(count))
  owners=[(mod(node-1,nranks),node=1,count)]
  do node=1,count
    geometry(node)=tree%nodes(node)%geometry
    if (owners(node)/=rank) cycle
    local%fields(node)%values=tree%nodes(node)%state
    allocate(lt%fields(node)%values(1,4,4,4))
    lt%fields(node)%values(1,:,:,:)=tree%nodes(node)%temperature
  end do
  call initialize_amr_patch_3d(4,4,4,2,3,3,4,2,3,2,boxes(3),ok)
  call require(ok,'shifted grandchild')
  call initialize_amr_patch_tree_3d([4,4,4],[1,1,2],boxes(:3),moved,ok)
  call require(ok,'new topology')
  call remap_reactive_eb_tree_3d(species,moved,finest,tree,ok)
  call require(ok,'move grandchild and remove other refinement')
  call check_inventory()
  call check_mpi_remap(tree%topology)
  call require(all(tree%nodes(4)%state(:,3:4,:,:)==saved%nodes(4)%state(:,1:2,:,:)), &
    'exact same-level overlap states')
  call require(all(tree%nodes(4)%temperature(3:4,:,:)==saved%nodes(4)%temperature(1:2,:,:)), &
    'exact same-level overlap temperatures')
  ! Root-only target coarsens all old descendants; the input finest geometry
  ! still describes the deeper of the source and target trees.
  call initialize_amr_patch_tree_3d([4,4,4],[integer::],boxes(:0),root_only,ok)
  call require(ok,'root-only topology')
  call remap_reactive_eb_tree_3d(species,root_only,finest,tree,ok)
  call require(ok.and.size(tree%nodes)==1,'complete coarsening')
  call check_inventory()
  call check_mpi_remap(tree%topology)
  call remap_reactive_eb_tree_3d(species,topology,finest,tree,ok)
  call require(ok,'refine root-only state to three levels')
  call check_inventory()
  call check_mpi_remap(tree%topology)
  original=local;original_t=lt;original_geometry=geometry
  if (owners(count)==rank) local%fields(count)%values(1,4,4,4)=-1._dp
  original=local
  next_owners=[(mod(node,nranks),node=1,size(moved%nodes))]
  call remap_mpi_eb_tree_3d(species,current_tree,moved,owners,next_owners,finest,geometry,local,lt,MPI_COMM_WORLD,ok)
  call require(.not.ok,'MPI bad source rejected')
  call unchanged_mpi()
  if (owners(count)==rank) local%fields(count)%values(1,4,4,4)=tree%nodes(count)%state(1,4,4,4)
  original=local
  call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
    axes(axis),.32_dp,finest,ok)
  call require(ok,'changed geometry fixture')
  call remap_mpi_eb_tree_3d(species,current_tree,moved,owners,next_owners,finest,geometry,local,lt,MPI_COMM_WORLD,ok)
  call require(.not.ok,'MPI geometry change rejected')
  call unchanged_mpi()
  deallocate(geometry,owners,local%fields,lt%fields)
  deallocate(u,t)
  end do
  end do
  if (rank==0) print *, 'test_mpi_eb_tree_remap_3d: PASS'
  call MPI_Finalize(status)
contains
  subroutine check_inventory()
    call integrate_reactive_eb_tree_3d(tree,after,ok)
    call require(ok,'remapped inventory')
    residual=maxval(abs(after-before)/max(1._dp,abs(before)))
    print *, 'axis, source nodes, residual: ',axes(axis),count,residual
    call require(residual<3.e-14_dp,'remap composite conservation')
  end subroutine check_inventory
  subroutine check_mpi_remap(target)
    type(amr_patch_tree_3d), intent(in) :: target
    real(dp) :: inventory(nv),error
    integer :: p
    stage=stage+1
    next_owners=[(mod(p+stage,nranks),p=1,size(target%nodes))]
    call remap_mpi_eb_tree_3d(species,current_tree,target,owners,next_owners,finest,geometry,local,lt,MPI_COMM_WORLD,ok)
    call require(ok,'MPI remap')
    current_tree=target;owners=next_owners
    do p=1,size(target%nodes)
      call require(allocated(local%fields(p)%values).eqv.(owners(p)==rank),'new owner-local allocation')
      call require(matching_eb_metrics_3d(geometry(p),tree%nodes(p)%geometry),'serial target geometry')
      if (owners(p)/=rank) cycle
      call require(all(local%fields(p)%values==tree%nodes(p)%state),'exact serial remapped U')
      call require(all(lt%fields(p)%values(1,:,:,:)==tree%nodes(p)%temperature),'exact serial remapped T')
    end do
    call integrate_mpi_eb_tree_3d(target,owners,geometry,local,inventory,MPI_COMM_WORLD,ok)
    call require(ok,'MPI inventory')
    error=maxval(abs(inventory-before)/max(1._dp,abs(before)))
    call require(error<3.e-14_dp,'MPI remap conservation')
    if (rank==0) print *, 'MPI axis, nodes, stage, residual: ',axes(axis),count,stage,error
  end subroutine
  subroutine unchanged_mpi()
    integer :: p
    do p=1,size(owners)
      call require(matching_eb_metrics_3d(geometry(p),original_geometry(p)),'geometry rollback')
      if (owners(p)/=rank) cycle
      call require(all(local%fields(p)%values==original%fields(p)%values),'U rollback')
      call require(all(lt%fields(p)%values==original_t%fields(p)%values),'T rollback')
    end do
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: ierr
    if (.not.condition) then
      print *, 'FAIL: ',label
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
  end subroutine require
end program test_mpi_eb_tree_remap_3d
