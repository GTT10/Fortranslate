program test_mpi_eb_tree_transfer_3d
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use mpi_amr_patch_tree_distribution_3d_mod, only: scatter_mpi_amr_tree_3d,gather_mpi_amr_tree_3d
  use mpi_eb_tree_transfer_3d_mod, only: average_down_mpi_eb_tree_3d,integrate_mpi_eb_tree_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,middle,coarse,geometry(5),original(5)
  type(eb_reactive_tree_3d) :: tree,serial
  type(reactive_eb_3d_config) :: config
  type(amr_tree_state_3d) :: global,gt,local,lt,result,rt,saved,st
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:),expected(:),actual(:)
  real(dp) :: error
  integer :: nv,node,count,rank,nranks,status,owners(5),root,mode
  logical :: ok
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  root=nranks-1;count=5
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  nv=reactive_nvar(size(species))
  allocate(expected(nv),actual(nv))
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
  serial=tree
  call average_down_reactive_eb_tree_3d(species,serial,ok)
  call require(ok,'serial average')
  call integrate_reactive_eb_tree_3d(serial,expected,ok)
  call require(ok,'serial inventory')
  do node=1,5
    geometry(node)=tree%nodes(node)%geometry
  end do
  original=geometry
  if (rank==root) then
    allocate(global%fields(5),gt%fields(5))
    do node=1,5
      global%fields(node)%values=tree%nodes(node)%state
      allocate(gt%fields(node)%values(1,4,4,4))
      gt%fields(node)%values(1,:,:,:)=tree%nodes(node)%temperature
    end do
  end if
  do mode=1,2
    owners=[(mod(node-1,nranks),node=1,5)]
    if (mode==2) owners=0
    call scatter_mpi_amr_tree_3d(topology,owners,global,local,MPI_COMM_WORLD,root,ok)
    call require(ok,'scatter U')
    call scatter_mpi_amr_tree_3d(topology,owners,gt,lt,MPI_COMM_WORLD,root,ok)
    call require(ok,'scatter T')
    call average_down_mpi_eb_tree_3d(species,topology,owners,geometry,local,lt,MPI_COMM_WORLD,ok)
    call require(ok,'MPI average')
    do node=1,5
      call require(allocated(local%fields(node)%values).eqv.(rank==owners(node)),'owner-local U')
      call require(allocated(lt%fields(node)%values).eqv.(rank==owners(node)),'owner-local T')
    end do
    call gather_mpi_amr_tree_3d(topology,owners,local,result,MPI_COMM_WORLD,root,ok)
    call require(ok,'gather U')
    call gather_mpi_amr_tree_3d(topology,owners,lt,rt,MPI_COMM_WORLD,root,ok)
    call require(ok,'gather T')
    if (rank==root) then
      do node=1,5
        call require(all(result%fields(node)%values==serial%nodes(node)%state),'exact serial U including solids')
        call require(all(rt%fields(node)%values(1,:,:,:)==serial%nodes(node)%temperature),'exact serial T')
      end do
    end if
    call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,actual,MPI_COMM_WORLD,ok)
    call require(ok,'MPI inventory')
    error=maxval(abs(actual-expected)/max(1._dp,abs(expected)))
    call require(error<3.e-14_dp,'composite inventory')
    if (rank==0) print *, 'ranks, mode, inventory error: ',nranks,mode,error
    saved=local;st=lt
    ! Reject a late child after another child has already restricted privately.
    if (rank==owners(4)) local%fields(4)%values(1,4,4,4)=-1._dp
    saved=local
    call average_down_mpi_eb_tree_3d(species,topology,owners,geometry,local,lt,MPI_COMM_WORLD,ok)
    call require(.not.ok,'bad fine EOS rejected collectively')
    call unchanged()
    local=saved
    if (rank==owners(4)) local%fields(4)%values(1,4,4,4)=serial%nodes(4)%state(1,4,4,4)
    saved=local
    if (rank==0) then
      do node=1,5
        geometry(node)%x_lower=geometry(node)%x_lower+.01_dp
        geometry(node)%x_upper=geometry(node)%x_upper+.01_dp
        call require(geometry(node)%is_valid(),'locally valid translated geometry')
      end do
    end if
    call average_down_mpi_eb_tree_3d(species,topology,owners,geometry,local,lt,MPI_COMM_WORLD,ok)
    call require(.not.ok,'geometry mismatch rejected collectively')
    call unchanged()
    geometry=original
    if (rank==owners(5)) local%fields(5)%values(1,4,4,4)=ieee_value(0._dp,ieee_quiet_nan)
    actual=17._dp
    call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,actual,MPI_COMM_WORLD,ok)
    call require(.not.ok.and.all(actual==17._dp),'invalid exposed payload preserves inventory')
  end do
  if (rank==0) print *, 'PASS: MPI EB restriction and composite inventory'
  call MPI_Finalize(status)
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: ierr
    if (condition) return
    print *, 'FAIL rank ',rank,': ',label
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  end subroutine
  subroutine unchanged()
    integer :: n
    do n=1,5
      if (owners(n)/=rank) cycle
      call require(all(local%fields(n)%values==saved%fields(n)%values),'state rollback')
      call require(all(lt%fields(n)%values==st%fields(n)%values),'temperature rollback')
    end do
  end subroutine
end program
