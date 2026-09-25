program test_mpi_eb_tree_hydro_3d
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar,reactive_conserved_to_primitive,reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,advance_reactive_eb_tree_hydro_3d, &
    build_eb_tree_extended_geometry_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use mpi_eb_tree_hydro_3d_mod, only: advance_mpi_eb_tree_hydro_3d
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,middle,coarse,extended_finest
  type(eb_geometry_3d), allocatable :: extended(:),physical(:)
  type(eb_reactive_tree_3d) :: tree,saved
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:)
  real(dp), allocatable :: before(:),after(:),boundary(:),wall(:),primitive(:)
  real(dp) :: residual,coefficient,et,sound
  character(len=32) :: test_mode,solver
  character(len=1), parameter :: axes(3)=['x','y','z']
  integer :: nv,node,axis,count,rank,nranks,status,i,j,k
  logical :: ok
  call get_command_argument(1,test_mode)
  solver='hllc'
  if(test_mode=='ctu') solver='pelec'
  coefficient=0
  if(test_mode=='viscosity'.or.test_mode=='ctu') coefficient=.037_dp
  if(test_mode=='default') coefficient=.1_dp
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  nv=reactive_nvar(size(species))
  allocate(before(nv),after(nv),boundary(nv),wall(nv),primitive(nv))
  do axis=1,3
  do count=4,5
  if(test_mode=='ctu'.and.(axis/=3.or.count/=5)) cycle
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
  if(test_mode=='ctu') then
    allocate(physical(count))
    do node=1,count
      physical(node)=tree%nodes(node)%geometry
      do k=1,4
      do j=1,4
      do i=1,4
        if(physical(node)%volume_fraction(i,j,k)<=0) cycle
        call reactive_conserved_to_primitive(species,tree%nodes(node)%state(:,i,j,k), &
          tree%nodes(node)%temperature(i,j,k),primitive,et,sound,ok)
        call require(ok,'CTU primitive fixture')
        primitive(5)=primitive(5)*1000._dp/et
        call reactive_primitive_to_conserved(species,primitive,tree%nodes(node)%state(:,i,j,k), &
          tree%nodes(node)%temperature(i,j,k),sound,ok)
        call require(ok,'CTU fixture EOS')
      end do
      end do
      end do
    end do
    call build_axis_plane_eb_geometry_3d(64,64,64,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
      axes(axis),.31_dp,extended_finest,ok)
    call require(ok,'extended finest')
    call build_eb_tree_extended_geometry_3d(topology,physical,extended_finest,extended,ok)
    call require(ok,'extended cache')
  end if
  call average_down_reactive_eb_tree_3d(species,tree,ok)
  call require(ok,'initial average down')
  call integrate_reactive_eb_tree_3d(tree,before,ok)
  call require(ok,'initial inventory')
  saved=tree
  boundary=17._dp;wall=19._dp
  call advance_reactive_eb_tree_hydro_3d(species,solver,1.e-8_dp,tree,boundary,wall,ok, &
    difmag=coefficient,extended_geometry=extended)
  call require(ok,'tree hydro step')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'final inventory')
  residual=maxval(abs(after-before-boundary-wall)/max(1._dp,abs(before)))
  print *, 'axis, nodes, inventory residual: ',axes(axis),count,residual
  call require(residual<3.e-14_dp,'hydro conservation including boundary and exposed wall')
  call require(any(tree%nodes(count)%state/=saved%nodes(count)%state),'nontrivial deep hydro update')
  call check_mpi_hydro()
  deallocate(u,t)
  end do
  end do
  if (rank==0) print *, 'test_mpi_eb_tree_hydro_3d: PASS'
  call MPI_Finalize(status)
contains
  subroutine check_mpi_hydro()
    type(amr_tree_state_3d) :: local,lt,initial,it
    type(eb_geometry_3d) :: geometry(count)
    type(eb_geometry_3d), allocatable :: trial_extended(:)
    real(dp), allocatable :: mpi_coefficient
    real(dp) :: b(nv),w(nv),inventory(nv),error,dt,local_coefficient
    integer :: owners(count),p,mode,failure
    logical :: good
    do p=1,count
      geometry(p)=saved%nodes(p)%geometry
    end do
    do mode=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if (mode==2) owners=nranks-1
      allocate(local%fields(count),lt%fields(count))
      do p=1,count
        if (owners(p)/=rank) cycle
        local%fields(p)%values=saved%nodes(p)%state
        allocate(lt%fields(p)%values(1,4,4,4))
        lt%fields(p)%values(1,:,:,:)=saved%nodes(p)%temperature
      end do
      initial=local;it=lt;b=17._dp;w=19._dp;local_coefficient=coefficient
      if(test_mode/='default') mpi_coefficient=coefficient
      call advance_mpi_eb_tree_hydro_3d(species,topology,owners,geometry,solver,1.e-8_dp, &
        local,lt,b,w,MPI_COMM_WORLD,good,difmag=mpi_coefficient,extended_geometry=extended)
      call require(good,'MPI hydro')
      do p=1,count
        call require(allocated(local%fields(p)%values).eqv.(owners(p)==rank),'owner-local U')
        if (owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==tree%nodes(p)%state),'exact serial hydro U')
        call require(all(lt%fields(p)%values(1,:,:,:)==tree%nodes(p)%temperature),'exact serial hydro T')
      end do
      call require(all(b==boundary),'exact boundary ledger')
      call require(maxval(abs(w-wall)/max(1._dp,abs(wall)))<3.e-14_dp,'wall ledger')
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,inventory,MPI_COMM_WORLD,good)
      call require(good,'MPI inventory')
      error=maxval(abs(inventory-before-b-w)/max(1._dp,abs(before)))
      call require(error<3.e-14_dp,'MPI hydro conservation')
      if (rank==0) print *, 'MPI axis, nodes, mode, residual: ',axes(axis),count,mode,error
      do failure=1,8
        if(failure>=7.and.test_mode/='ctu') cycle
        if(allocated(extended)) trial_extended=extended
        if(failure==7.and.rank==0) deallocate(trial_extended)
        if(failure==8.and.rank==0) trial_extended(1)%x_face_centroid_y(0,1,1)=.1_dp
        if(failure==4.and.coefficient==0) cycle
        do p=1,count
          geometry(p)=saved%nodes(p)%geometry
        end do
        if(failure==5.and.rank==0) geometry(1)%halo_volume_fraction(0,0,0)=.5_dp
        if(failure==6.and.rank==0) deallocate(geometry(1)%halo_volume_fraction)
        local_coefficient=coefficient
        if(failure==4.and.rank==0) local_coefficient=.05_dp
        local=initial;lt=it;b=17._dp;w=19._dp;dt=1.e-8_dp
        if (failure==1) then
          call advance_mpi_eb_tree_hydro_3d(species,topology,owners,geometry,'invalid',dt, &
            local,lt,b,w,MPI_COMM_WORLD,good,difmag=local_coefficient,extended_geometry=trial_extended)
        else
          if (failure==2.and.rank==0) dt=2.e-8_dp
          if (failure==3.and.rank==owners(count)) local%fields(count)%values(1,4,4,4)=-1._dp
          call advance_mpi_eb_tree_hydro_3d(species,topology,owners,geometry,solver,dt, &
            local,lt,b,w,MPI_COMM_WORLD,good,difmag=local_coefficient,extended_geometry=trial_extended)
        end if
        call require(.not.good,'collective hydro rejection')
        if (failure==3.and.rank==owners(count)) then
          call require(local%fields(count)%values(1,4,4,4)==-1._dp,'bad density preserved')
          local%fields(count)%values(1,4,4,4)=initial%fields(count)%values(1,4,4,4)
        end if
        do p=1,count
          if (owners(p)/=rank) cycle
          call require(all(local%fields(p)%values==initial%fields(p)%values),'U rollback')
          call require(all(lt%fields(p)%values==it%fields(p)%values),'T rollback')
        end do
        call require(all(b==17._dp).and.all(w==19._dp),'ledger rollback')
      end do
      do p=1,count
        geometry(p)=saved%nodes(p)%geometry
      end do
      deallocate(local%fields,lt%fields)
    end do
  end subroutine check_mpi_hydro
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: ierr
    if (.not.condition) then
      print *, 'FAIL: ',label
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
  end subroutine require
end program test_mpi_eb_tree_hydro_3d
