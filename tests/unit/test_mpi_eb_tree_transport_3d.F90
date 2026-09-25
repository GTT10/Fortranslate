program test_mpi_eb_tree_transport_3d
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,advance_reactive_eb_tree_transport_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use mpi_eb_tree_timestep_3d_mod, only: mpi_eb_tree_timestep_3d
  use eb_reactive_patch_tree_3d_mod, only: reactive_eb_tree_timestep_3d
  use mpi_eb_tree_transport_3d_mod, only: advance_mpi_eb_tree_transport_3d
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d
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
  integer :: nv,node,axis,count,i,j,k,rank,nranks,status
  character(len=16) :: test_mode
  character(len=16), allocatable :: labels(:)
  real(dp), allocatable :: wall_t(:)
  logical :: ok
  call get_command_argument(1,test_mode)
  if(test_mode=='isothermal') then
    allocate(labels(6),wall_t(6));labels='isothermal_wall';wall_t=2000
  end if
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  nv=reactive_nvar(size(species))
  allocate(before(nv),after(nv),boundary(nv),primitive(reactive_nprim(size(species))))
  do axis=1,3
  do count=4,5
  if(allocated(labels).and.(axis/=3.or.count/=5)) cycle
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
    tree,boundary,theta,ok,boundary_faces=labels,wall_temperature=wall_t)
  call require(ok,'tree transport step')
  if(allocated(labels)) then
    call require(boundary(1)==0.and.all(boundary(6:)==0).and.boundary(5)>0,'heated impermeable walls')
  else
    call require(all(boundary==0._dp),'closed physical diffusion faces')
  end if
  call require(theta>=0._dp.and.theta<=1._dp,'limiter range')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'final inventory')
  residual=maxval(abs(after-before-boundary)/max(1._dp,abs(before)))
  print *, 'axis, nodes, inventory residual: ',axes(axis),count,residual
  call require(residual<3.e-14_dp,'transport conservation')
  call require(any(tree%nodes(count)%state/=saved%nodes(count)%state),'nontrivial deep transport update')
  call check_mpi_transport()
  deallocate(u,t)
  end do
  end do
  if (rank==0) print *, 'test_mpi_eb_tree_transport_3d: PASS'
  call MPI_Finalize(status)
contains
  subroutine check_mpi_transport()
    type(amr_tree_state_3d) :: local,lt,initial,it
    type(eb_geometry_3d) :: geometry(count)
    real(dp) :: mach(6),tin(6),pin(6),old_dt
    real(dp) :: b(nv),limiter,inventory(nv),error,dt,expected_dt,actual_dt,hcfl,tcfl
    integer :: owners(count),p,mode,failure,estimate_mode,recursive_mode
    logical :: good,visc
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
      initial=local;it=lt;b=17._dp;limiter=19._dp
      if(allocated(wall_t)) then
        do recursive_mode=0,1
          call reactive_eb_tree_timestep_3d(species,transport,saved,.5_dp,1.e-12_dp,.true., &
            .false.,.true.,.false.,old_dt,good,subcycling=recursive_mode==1)
          call require(good,'thermal baseline timestep')
          call reactive_eb_tree_timestep_3d(species,transport,saved,.5_dp,1.e-12_dp,.true., &
            .false.,.true.,.false.,expected_dt,good,subcycling=recursive_mode==1, &
            boundary_faces=labels,wall_temperature=wall_t)
          call require(good.and.expected_dt<=old_dt,'serial thermal bound')
          if(recursive_mode==1) call require(expected_dt<old_dt,'thermal bound is active')
          call mpi_eb_tree_timestep_3d(species,transport,topology,owners,geometry,local,lt,.5_dp,1.e-12_dp,.true., &
            .false.,.true.,.false.,actual_dt,MPI_COMM_WORLD,good,subcycling=recursive_mode==1, &
            boundary_faces=labels,wall_temperature=wall_t)
          call require(good.and.actual_dt==expected_dt,'exact MPI thermal timestep')
          if(rank==0) print *, 'thermal dt baseline, serial, MPI: ',old_dt,expected_dt,actual_dt
          actual_dt=17
          if(rank==0) wall_t(1)=2001
          call mpi_eb_tree_timestep_3d(species,transport,topology,owners,geometry,local,lt,.5_dp,1.e-12_dp,.true., &
            .false.,.true.,.false.,actual_dt,MPI_COMM_WORLD,good,subcycling=recursive_mode==1, &
            boundary_faces=labels,wall_temperature=wall_t)
          call require(.not.good.and.actual_dt==17,'rank-dependent wall temperature preserves timestep')
          wall_t=2000
        end do
      end if

      do recursive_mode=0,1
      do estimate_mode=1,3
        hcfl=.5_dp;tcfl=.35_dp
        if (estimate_mode==3) tcfl=1.e-10_dp
        call reactive_eb_tree_timestep_3d(species,transport,saved,hcfl,tcfl,estimate_mode/=1, &
          .true.,.true.,.true.,expected_dt,good,subcycling=recursive_mode==1)
        call require(good,'serial estimate')
        call mpi_eb_tree_timestep_3d(species,transport,topology,owners,geometry,local,lt,hcfl,tcfl,estimate_mode/=1, &
          .true.,.true.,.true.,actual_dt,MPI_COMM_WORLD,good,subcycling=recursive_mode==1)
        call require(good.and.actual_dt==expected_dt,'exact MPI hydro/diffusion bound')
        if (axis==3.and.count==5) then
          old_dt=expected_dt;mach=10;tin=2000;pin=101325
          call reactive_eb_tree_timestep_3d(species,transport,saved,hcfl,tcfl,estimate_mode/=1, &
            .true.,.true.,.true.,expected_dt,good,subcycling=recursive_mode==1, &
            inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
          call require(good.and.expected_dt<old_dt,'hot fast inlet tightens hydro and diffusion bound')
          call mpi_eb_tree_timestep_3d(species,transport,topology,owners,geometry,local,lt,hcfl,tcfl,estimate_mode/=1, &
            .true.,.true.,.true.,actual_dt,MPI_COMM_WORLD,good,subcycling=recursive_mode==1, &
            inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
          call require(good.and.actual_dt==expected_dt,'exact MPI inlet timestep')
          actual_dt=17
          if (rank==0) pin(1)=101326
          call mpi_eb_tree_timestep_3d(species,transport,topology,owners,geometry,local,lt,hcfl,tcfl,estimate_mode/=1, &
            .true.,.true.,.true.,actual_dt,MPI_COMM_WORLD,good,subcycling=recursive_mode==1, &
            inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
          call require(.not.good.and.actual_dt==17,'rank-dependent inlet preserves timestep')
        end if
      end do
      end do
      actual_dt=17._dp
      if (rank==0) hcfl=.25_dp
      call mpi_eb_tree_timestep_3d(species,transport,topology,owners,geometry,local,lt,hcfl,tcfl,.true., &
        .true.,.true.,.true.,actual_dt,MPI_COMM_WORLD,good)
      call require(.not.good.and.actual_dt==17._dp,'inconsistent CFL preserves estimate')
      call mpi_eb_tree_timestep_3d(species,transport,topology,owners,geometry,local,lt,.5_dp,tcfl,.true., &
        .true.,.true.,.true.,actual_dt,MPI_COMM_WORLD,good,subcycling=rank==0)
      call require(.not.good.and.actual_dt==17._dp,'inconsistent subcycling preserves estimate')
      do p=1,count
        if (owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==initial%fields(p)%values),'estimate preserves U')
        call require(all(lt%fields(p)%values==it%fields(p)%values),'estimate preserves T')
      end do
      call advance_mpi_eb_tree_transport_3d(species,transport,topology,owners,geometry,1.e-8_dp, &
        .true.,.true.,.true.,.true.,local,lt,b,limiter,MPI_COMM_WORLD,good, &
          boundary_faces=labels,wall_temperature=wall_t)
      call require(good,'MPI transport')
      do p=1,count
        call require(allocated(local%fields(p)%values).eqv.(owners(p)==rank),'owner-local U')
        if (owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==tree%nodes(p)%state),'exact serial transport U')
        call require(all(lt%fields(p)%values(1,:,:,:)==tree%nodes(p)%temperature),'exact serial transport T')
      end do
      call require(all(b==boundary),'exact boundary ledger')
      call require(limiter==theta,'exact serial minimum limiter')
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,inventory,MPI_COMM_WORLD,good)
      call require(good,'MPI inventory')
      error=maxval(abs(inventory-before-b)/max(1._dp,abs(before)))
      call require(error<3.e-14_dp,'MPI transport conservation')
      if (rank==0) print *, 'MPI axis, nodes, mode, residual: ',axes(axis),count,mode,error
      do failure=1,5
        if(failure==5.and..not.allocated(wall_t)) cycle
        local=initial;lt=it;b=17._dp;limiter=19._dp;dt=1.e-8_dp
        visc=.true.
        if (failure==1) dt=-1.e-8_dp
        if (failure==2.and.rank==0) dt=2.e-8_dp
        if (failure==3.and.rank==owners(count)) local%fields(count)%values(1,4,4,4)=-1._dp
        if (failure==4.and.rank==0) visc=.false.
        if(failure==5.and.rank==0) wall_t(1)=2001
        call advance_mpi_eb_tree_transport_3d(species,transport,topology,owners,geometry,dt, &
          visc,.true.,.true.,.true.,local,lt,b,limiter,MPI_COMM_WORLD,good, &
          boundary_faces=labels,wall_temperature=wall_t)
        if(allocated(wall_t)) wall_t=2000
        call require(.not.good,'collective transport rejection')
        if (failure==3.and.rank==owners(count)) then
          call require(local%fields(count)%values(1,4,4,4)==-1._dp,'bad density preserved')
          local%fields(count)%values(1,4,4,4)=initial%fields(count)%values(1,4,4,4)
        end if
        do p=1,count
          if (owners(p)/=rank) cycle
          call require(all(local%fields(p)%values==initial%fields(p)%values),'U rollback')
          call require(all(lt%fields(p)%values==it%fields(p)%values),'T rollback')
        end do
        call require(all(b==17._dp).and.limiter==19._dp,'ledger rollback')
      end do
      deallocate(local%fields,lt%fields)
    end do
  end subroutine check_mpi_transport
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: ierr
    if (.not.condition) then
      print *, 'FAIL: ',label
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
  end subroutine require
end program test_mpi_eb_tree_transport_3d
