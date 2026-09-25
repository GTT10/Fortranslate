program test_mpi_eb_tree_subcycle_full_3d
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive, &
    reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,advance_reactive_eb_tree_full_3d, &
    reactive_eb_tree_timestep_3d,build_eb_tree_extended_geometry_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use mpi_eb_tree_timestep_3d_mod, only: mpi_eb_tree_timestep_3d
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_full_3d
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,middle,coarse,extended_finest
  type(eb_geometry_3d), allocatable :: extended(:),physical(:)
  type(eb_reactive_tree_3d) :: tree,saved
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:)
  real(dp), allocatable :: before(:),after(:),boundary(:),primitive(:),wall(:),chemical(:)
  real(dp) :: residual,theta,speed,et
  character(len=1), parameter :: axes(3)=['x','y','z']
  integer :: nv,node,axis,count,i,j,k,mode,integrator_index,rank,nranks,status
  logical :: ok,chemistry_on,transport_on,recursive_steps
  real(dp), allocatable :: mach(:),pin(:),inlet(:,:),wall_temperature(:)
  real(dp) :: tin(6)
  character(len=16), allocatable :: labels(:)
  character(len=32) :: test_mode,viscosity_mode,solver
  real(dp) :: coefficient
  character(len=8), parameter :: integrators(2)=[character(len=8)::'explicit','implicit']
  call get_command_argument(2,viscosity_mode)
  coefficient=0
  if(viscosity_mode=='viscosity') coefficient=.037_dp
  call get_command_argument(1,test_mode)
  recursive_steps=index(test_mode,'common')==0
  solver='hllc'
  if(index(test_mode,'ctu')==1) solver='pelec'
  allocate(mach(6),pin(6));mach=-1;tin=400;pin=0
  if(index(test_mode,'inflow')>0) then
    deallocate(mach,pin)
    allocate(labels(6));labels='outflow';labels(1)='inflow';labels(4)='inflow'
  end if
  if(index(test_mode,'periodic')>0) then
    deallocate(mach,pin)
    allocate(labels(6));labels='periodic';labels(5:6)='outflow'
  end if
  if(index(test_mode,'wall')>0) then
    deallocate(mach,pin)
    allocate(labels(6));labels='outflow';labels(1:2)='slip_wall';labels(3:4)='no_slip_wall'
  end if
  if(index(test_mode,'isothermal')>0) then
    deallocate(mach,pin)
    allocate(labels(6),wall_temperature(6));labels='isothermal_wall';labels(5:6)='outflow'
    wall_temperature=2000
  end if
  if (index(test_mode,'mach')==1) then
    mach(1)=.2_dp;mach(4)=.1_dp;pin(4)=120000
  end if
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  call load_h2o2_elementary_mechanism(reactions,ok)
  call require(ok,'reactions')
  nv=reactive_nvar(size(species))
  allocate(before(nv),after(nv),boundary(nv),primitive(reactive_nprim(size(species))))
  allocate(wall(nv),chemical(nv))
  do axis=3,3
  do count=5,5
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
  if(solver=='pelec') then
    allocate(physical(count))
    do node=1,count
      physical(node)=tree%nodes(node)%geometry
      do k=1,4
      do j=1,4
      do i=1,4
        if(physical(node)%volume_fraction(i,j,k)<=0) cycle
        call reactive_conserved_to_primitive(species,tree%nodes(node)%state(:,i,j,k), &
          tree%nodes(node)%temperature(i,j,k),primitive,et,speed,ok)
        call require(ok,'CTU fixture primitive')
        primitive(5)=primitive(5)*1000._dp/et
        call reactive_primitive_to_conserved(species,primitive,tree%nodes(node)%state(:,i,j,k), &
          tree%nodes(node)%temperature(i,j,k),speed,ok)
        call require(ok,'CTU fixture EOS')
      end do
      end do
      end do
    end do
    call build_axis_plane_eb_geometry_3d(64,64,64,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
      axes(axis),.31_dp,extended_finest,ok)
    call require(ok,'extended finest')
    call build_eb_tree_extended_geometry_3d(topology,physical,extended_finest,extended,ok)
    call require(ok,'extended tree geometry')
  end if
  call average_down_reactive_eb_tree_3d(species,tree,ok)
  call require(ok,'initial average down')
  call integrate_reactive_eb_tree_3d(tree,before,ok)
  call require(ok,'initial inventory')
  saved=tree
  if(index(test_mode,'inflow')>0) then
    allocate(inlet(nv,6))
    call reactive_conserved_to_primitive(species,tree%nodes(1)%state(:,1,1,3), &
      tree%nodes(1)%temperature(1,1,3),primitive,et,speed,ok)
    call require(ok,'prescribed inlet primitive')
    primitive(2:4)=[100._dp,-50._dp,0._dp];primitive(5)=.8_dp*primitive(5)
    do i=1,6
      call reactive_primitive_to_conserved(species,primitive,inlet(:,i),tin(i),speed,ok)
      call require(ok,'prescribed inlet EOS')
    end do
  end if
  do integrator_index=1,2
  do mode=merge(3,0,index(test_mode,'mach')==1),3
    if (integrator_index==2.and.mode/=3) cycle
    if(solver=='pelec'.and.mode/=3) cycle
    chemistry_on=btest(mode,0);transport_on=btest(mode,1)
    tree=saved
    call advance_reactive_eb_tree_full_3d(species,reactions,transport,solver,1.e-8_dp,chemistry_on, &
      2.e-7_dp,1.e-12_dp,transport_on,.true.,.true.,.true.,.true.,tree,boundary,wall,chemical,theta,ok, &
      integrators(integrator_index),subcycling=recursive_steps, &
      inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin,difmag=coefficient,extended_geometry=extended, &
      boundary_faces=labels,inflow_state=inlet,wall_temperature=wall_temperature)
    call require(ok,'combined step')
    call integrate_reactive_eb_tree_3d(tree,after,ok)
    call require(ok,'final inventory')
    residual=maxval(abs(after-before-boundary-wall-chemical)/max(1._dp,abs(before)))
    print *, 'integrator, mode, residual: ',integrators(integrator_index),mode,residual
    call require(residual<3.e-14_dp,'combined conservation')
    call require(theta>=0._dp.and.theta<=1._dp,'limiter range')
    if (.not.transport_on) call require(theta==1._dp,'disabled transport limiter')
    if (.not.chemistry_on) call require(all(chemical==0._dp),'disabled chemistry ledger')
    if (chemistry_on) then
      call require(any(chemical(6:)/=0._dp),'nonzero chemical species exchange')
      call require(all(abs(chemical(:5))<=3.e-14_dp*max(1._dp,abs(before(:5)))),'chemical invariants')
    end if
    call check_mpi_full()
  end do
  end do
  deallocate(u,t)
  end do
  end do
  if (rank==0) print *, 'test_mpi_eb_tree_subcycle_full_3d: PASS'
  call MPI_Finalize(status)
contains
  subroutine check_mpi_full()
    type(amr_tree_state_3d) :: local,lt,initial,it
    type(eb_geometry_3d) :: geometry(count)
    real(dp) :: b(nv),w(nv),chem(nv),limiter,inventory(nv),error,local_coefficient,serial_dt,mpi_dt,saved_tin
    integer :: owners(count),p
    logical :: good,enabled
    owners=[(mod(p-1,nranks),p=1,count)]
    if (integrator_index==2) owners=nranks-1
    allocate(local%fields(count),lt%fields(count))
    do p=1,count
      geometry(p)=saved%nodes(p)%geometry
      if (owners(p)/=rank) cycle
      local%fields(p)%values=saved%nodes(p)%state
      allocate(lt%fields(p)%values(1,4,4,4))
      lt%fields(p)%values(1,:,:,:)=saved%nodes(p)%temperature
    end do
    initial=local;it=lt
    if(allocated(labels)) then
      call reactive_eb_tree_timestep_3d(species,transport,saved,.1_dp,.1_dp, &
        .true.,.true.,.true.,.true.,serial_dt,good,subcycling=recursive_steps, &
        boundary_faces=labels,inflow_state=inlet,inflow_temperature=tin,wall_temperature=wall_temperature)
      call require(good,'serial prescribed inlet timestep')
      call mpi_eb_tree_timestep_3d(species,transport,topology,owners,geometry,local,lt,.1_dp,.1_dp, &
        .true.,.true.,.true.,.true.,mpi_dt,MPI_COMM_WORLD,good,subcycling=recursive_steps, &
        boundary_faces=labels,inflow_state=inlet,inflow_temperature=tin,wall_temperature=wall_temperature)
      call require(good.and.abs(mpi_dt-serial_dt)<=1e-14_dp*serial_dt,'MPI inlet timestep parity')
    end if
    call advance_mpi_eb_tree_full_3d(species,reactions,transport,topology,owners,geometry,solver,1.e-8_dp, &
      chemistry_on,2.e-7_dp,1.e-12_dp,transport_on,.true.,.true.,.true.,.true., &
      local,lt,b,w,chem,limiter,MPI_COMM_WORLD,good,integrators(integrator_index),subcycling=recursive_steps, &
      inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin,difmag=coefficient,extended_geometry=extended, &
      boundary_faces=labels,inflow_state=inlet,wall_temperature=wall_temperature)
    call require(good,'MPI full step')
    do p=1,count
      call require(allocated(local%fields(p)%values).eqv.(owners(p)==rank),'owner-local full U')
      if (owners(p)/=rank) cycle
      call require(all(local%fields(p)%values==tree%nodes(p)%state),'exact serial full U')
      call require(all(lt%fields(p)%values(1,:,:,:)==tree%nodes(p)%temperature),'exact serial full T')
    end do
    call require(all(b==boundary).and.limiter==theta,'exact boundary and limiter')
    call require(maxval(abs(w-wall)/max(1._dp,abs(before)))<3.e-14_dp,'wall parity')
    call require(maxval(abs(chem-chemical)/max(1._dp,abs(before)))<3.e-14_dp,'chemical ledger parity')
    call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,inventory,MPI_COMM_WORLD,good)
    call require(good,'MPI inventory')
    error=maxval(abs(inventory-before-b-w-chem)/max(1._dp,abs(before)))
    call require(error<3.e-14_dp,'MPI full conservation')
    if (rank==0) print *, 'MPI integrator, mode, residual: ',integrators(integrator_index),mode,error
    if (integrator_index==1.and.mode==3) then
      local=initial;lt=it;b=17._dp;w=19._dp;chem=23._dp;limiter=29._dp
      call advance_mpi_eb_tree_full_3d(species,reactions,transport,topology,owners,geometry,'invalid',1.e-8_dp, &
        .true.,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true., &
        local,lt,b,w,chem,limiter,MPI_COMM_WORLD,good,'explicit',subcycling=recursive_steps, &
        inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin,difmag=coefficient,extended_geometry=extended, &
      boundary_faces=labels,inflow_state=inlet,wall_temperature=wall_temperature)
      call require(.not.good,'late hydro rejection after reaction and transport')
      do p=1,count
        if (owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==initial%fields(p)%values),'full U rollback')
        call require(all(lt%fields(p)%values==it%fields(p)%values),'full T rollback')
      end do
      call require(all(b==17._dp).and.all(w==19._dp).and.all(chem==23._dp).and.limiter==29._dp,'full ledger rollback')
      if (index(test_mode,'mach')==1) then
        if (rank==0) pin(4)=120001
        call advance_mpi_eb_tree_full_3d(species,reactions,transport,topology,owners,geometry,solver,1.e-8_dp, &
          .false.,2.e-7_dp,1.e-12_dp,.false.,.true.,.true.,.true.,.true., &
          local,lt,b,w,chem,limiter,MPI_COMM_WORLD,good,'explicit',subcycling=recursive_steps, &
          inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin,difmag=coefficient,extended_geometry=extended, &
      boundary_faces=labels,inflow_state=inlet,wall_temperature=wall_temperature)
        call require(.not.good,'rank-dependent inlet pressure rejected')
        pin(4)=120000
        if (rank==0) tin(1)=ieee_value(0._dp,ieee_quiet_nan)
        call advance_mpi_eb_tree_full_3d(species,reactions,transport,topology,owners,geometry,solver,1.e-8_dp, &
          .false.,2.e-7_dp,1.e-12_dp,.false.,.true.,.true.,.true.,.true., &
          local,lt,b,w,chem,limiter,MPI_COMM_WORLD,good,'explicit',subcycling=recursive_steps, &
          inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin,difmag=coefficient,extended_geometry=extended, &
      boundary_faces=labels,inflow_state=inlet,wall_temperature=wall_temperature)
        call require(.not.good,'one-rank invalid inlet collectively rejected')
        tin=400
        do p=1,count
          if (owners(p)/=rank) cycle
          call require(all(local%fields(p)%values==initial%fields(p)%values),'inlet U rollback')
          call require(all(lt%fields(p)%values==it%fields(p)%values),'inlet T rollback')
        end do
        call require(all(b==17).and.all(w==19).and.all(chem==23).and.limiter==29,'inlet ledger rollback')
      end if
      if(allocated(labels)) then
        saved_tin=tin(1)
        if(rank==0) then
          if(allocated(wall_temperature)) then
            wall_temperature(1)=2001
          else if(index(test_mode,'periodic')>0.or.index(test_mode,'wall')>0) then
            labels(1)='outflow'
          else
            tin(1)=ieee_value(0._dp,ieee_quiet_nan)
          end if
        end if
        call advance_mpi_eb_tree_full_3d(species,reactions,transport,topology,owners,geometry,solver,1.e-8_dp, &
          .true.,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true., &
          local,lt,b,w,chem,limiter,MPI_COMM_WORLD,good,'explicit',subcycling=recursive_steps, &
          inflow_temperature=tin,difmag=coefficient,extended_geometry=extended, &
          boundary_faces=labels,inflow_state=inlet,wall_temperature=wall_temperature)
        call require(.not.good,'invalid boundary after private reaction')
        tin(1)=saved_tin
        if(allocated(wall_temperature)) wall_temperature=2000
        if(index(test_mode,'periodic')>0) labels(1)='periodic'
        if(index(test_mode,'wall')>0) labels(1)='slip_wall'
        do p=1,count
          if(owners(p)/=rank) cycle
          call require(all(local%fields(p)%values==initial%fields(p)%values).and. &
            all(lt%fields(p)%values==it%fields(p)%values),'boundary full rollback')
        end do
        call require(all(b==17).and.all(w==19).and.all(chem==23).and.limiter==29,'boundary ledger rollback')
      end if
      enabled=rank==0
      call advance_mpi_eb_tree_full_3d(species,reactions,transport,topology,owners,geometry,solver,1.e-8_dp, &
        enabled,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true., &
        local,lt,b,w,chem,limiter,MPI_COMM_WORLD,good,'explicit',subcycling=recursive_steps, &
        inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin,difmag=coefficient,extended_geometry=extended, &
      boundary_faces=labels,inflow_state=inlet,wall_temperature=wall_temperature)
      call require(.not.good,'different chemistry flags rejected')
      call advance_mpi_eb_tree_full_3d(species,reactions,transport,topology,owners,geometry,solver,1.e-8_dp, &
        .true.,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true., &
        local,lt,b,w,chem,limiter,MPI_COMM_WORLD,good,'explicit',subcycling=enabled)
      call require(.not.good,'different subcycling flags rejected')
      local_coefficient=coefficient
      if(rank==0) local_coefficient=.05_dp
      call advance_mpi_eb_tree_full_3d(species,reactions,transport,topology,owners,geometry,solver,1.e-8_dp, &
        .true.,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true., &
        local,lt,b,w,chem,limiter,MPI_COMM_WORLD,good,'explicit',subcycling=recursive_steps,difmag=local_coefficient)
      call require(.not.good,'different viscosity coefficients rejected')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==initial%fields(p)%values),'coefficient U rollback')
        call require(all(lt%fields(p)%values==it%fields(p)%values),'coefficient T rollback')
      end do
      call require(all(b==17._dp).and.all(w==19._dp).and.all(chem==23._dp).and.limiter==29._dp,'mode ledger rollback')
    end if
  end subroutine check_mpi_full
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: ierr
    if (.not.condition) then
      print *, 'FAIL: ',label
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
  end subroutine require
end program test_mpi_eb_tree_subcycle_full_3d
