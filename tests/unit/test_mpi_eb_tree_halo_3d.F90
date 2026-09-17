program test_mpi_eb_tree_halo_3d
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    fill_reactive_eb_tree_halo_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use mpi_amr_patch_tree_distribution_3d_mod, only: scatter_mpi_amr_tree_3d
  use mpi_eb_tree_halo_3d_mod, only: fill_mpi_eb_tree_halo_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,middle,coarse,geometry(5),original(5)
  type(eb_reactive_tree_3d) :: tree
  type(reactive_eb_3d_config) :: config
  type(amr_tree_state_3d) :: global,gt,local,lt,ghost,gh_t,saved,st,saved_ghost,saved_gt
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:)
  real(dp), allocatable :: reference(:,:,:,:),reference_t(:,:,:)
  integer :: nv,node,count,rank,nranks,status,owners(5),root,mode
  logical :: ok
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  root=nranks-1;count=5
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
  allocate(reference(nv,0:5,0:5,0:5),reference_t(0:5,0:5,0:5))
  do mode=1,2
    owners=[(mod(node-1,nranks),node=1,5)]
    if (mode==2) owners=0
    call scatter_mpi_amr_tree_3d(topology,owners,global,local,MPI_COMM_WORLD,root,ok)
    call require(ok,'scatter U')
    call scatter_mpi_amr_tree_3d(topology,owners,gt,lt,MPI_COMM_WORLD,root,ok)
    call require(ok,'scatter T')
    saved=local;st=lt
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,ghost,gh_t,MPI_COMM_WORLD,ok)
    call require(ok,'MPI halos')
    call unchanged_inputs()
    do node=1,5
      call require(allocated(ghost%fields(node)%values).eqv.(rank==owners(node)),'owner-local halo U')
      call require(allocated(gh_t%fields(node)%values).eqv.(rank==owners(node)),'owner-local halo T')
      if (rank/=owners(node)) cycle
      call require(all(lbound(ghost%fields(node)%values)==[1,0,0,0]),'halo lower bounds')
      call require(all(lbound(gh_t%fields(node)%values)==[1,0,0,0]),'temperature halo lower bounds')
      call fill_reactive_eb_tree_halo_3d(species,tree,node,reference,reference_t,ok)
      call require(ok,'serial halos')
      call require(all(ghost%fields(node)%values==reference),'exact serial face U and inactive zeros')
      call require(all(gh_t%fields(node)%values(1,:,:,:)==reference_t),'exact serial face T')
    end do
    ! Root unrefined boundary survives averaging, then fails donor EOS after exchange.
    if (rank==owners(1)) local%fields(1)%values(1,4,4,4)=-1._dp
    saved=local;saved_ghost=ghost;saved_gt=gh_t
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,ghost,gh_t,MPI_COMM_WORLD,ok)
    call require(.not.ok,'bad donor EOS rejected collectively')
    call unchanged_inputs()
    call unchanged_halos()
    if (rank==owners(1)) local%fields(1)%values(1,4,4,4)=tree%nodes(1)%state(1,4,4,4)
    saved=local
    if (rank==0) then
      do node=1,5
        geometry(node)%x_lower=geometry(node)%x_lower+.01_dp
        geometry(node)%x_upper=geometry(node)%x_upper+.01_dp
      end do
    end if
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,ghost,gh_t,MPI_COMM_WORLD,ok)
    call require(.not.ok,'geometry mismatch rejected collectively')
    call unchanged_inputs()
    call unchanged_halos()
    geometry=original
  end do
  if (rank==0) print *, 'PASS: MPI EB face halos, exact serial parity and rollback'
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
  subroutine unchanged_inputs()
    integer :: n
    do n=1,5
      if (owners(n)/=rank) cycle
      call require(all(local%fields(n)%values==saved%fields(n)%values),'input U unchanged')
      call require(all(lt%fields(n)%values==st%fields(n)%values),'input T unchanged')
    end do
  end subroutine
  subroutine unchanged_halos()
    integer :: n
    do n=1,5
      if (owners(n)/=rank) cycle
      call require(all(ghost%fields(n)%values==saved_ghost%fields(n)%values),'halo U rollback')
      call require(all(gh_t%fields(n)%values==saved_gt%fields(n)%values),'halo T rollback')
    end do
  end subroutine
end program
