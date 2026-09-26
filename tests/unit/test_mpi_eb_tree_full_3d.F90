program test_mpi_eb_tree_full_3d
  use iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use mpi_f08
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
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,advance_reactive_eb_tree_full_3d,react_sdc_eb_tree_3d
  use eb_reactive_patch_tree_3d_mod, only: advance_reactive_eb_tree_sdc_3d,build_eb_tree_extended_geometry_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_full_3d,react_sdc_mpi_eb_tree_3d,advance_mpi_eb_tree_sdc_3d
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_sdc_subcycle_3d,mpi_eb_sdc_recursive_source
  use eb_reactive_patch_tree_3d_mod, only: advance_reactive_eb_tree_sdc_subcycle_3d,eb_sdc_recursive_source
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,middle,coarse
  type(eb_reactive_tree_3d) :: tree,saved
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:)
  real(dp), allocatable :: before(:),after(:),boundary(:),primitive(:),wall(:),chemical(:)
  real(dp) :: residual,theta,speed
  character(len=1), parameter :: axes(3)=['x','y','z']
  integer :: nv,node,axis,count,i,j,k,mode,integrator_index,rank,nranks,status
  integer :: serial_source_calls=0,mpi_source_calls=0,source_failure=0
  character(len=16) :: test_mode
  logical :: ok,chemistry_on,transport_on
  character(len=8), parameter :: integrators(2)=[character(len=8)::'explicit','implicit']
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  call load_h2o2_elementary_mechanism(reactions,ok)
  call require(ok,'reactions')
  call get_command_argument(1,test_mode)
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
  call average_down_reactive_eb_tree_3d(species,tree,ok)
  call require(ok,'initial average down')
  call integrate_reactive_eb_tree_3d(tree,before,ok)
  call require(ok,'initial inventory')
  saved=tree
  call check_sdc_loop()
  if(test_mode=="subcycle".or.test_mode=="redist_level".or.test_mode=="redist_recursive".or. &
     test_mode=="state_redist".or.test_mode=="state_recursive".or.test_mode=="aux_source".or. &
     test_mode=="aux_recursive".or.test_mode=="aux_route") then
    deallocate(u,t)
    cycle
  end if
  call check_sdc_reaction()
  do integrator_index=1,2
  do mode=0,3
    if (integrator_index==2.and.mode/=3) cycle
    chemistry_on=btest(mode,0);transport_on=btest(mode,1)
    tree=saved
    call advance_reactive_eb_tree_full_3d(species,reactions,transport,'hllc',1.e-8_dp,chemistry_on, &
      2.e-7_dp,1.e-12_dp,transport_on,.true.,.true.,.true.,.true.,tree,boundary,wall,chemical,theta,ok, &
      integrators(integrator_index))
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
  if (rank==0) print *, 'test_mpi_eb_tree_full_3d: PASS'
  call MPI_Finalize(status)
contains
  subroutine check_sdc_loop()
    type(eb_reactive_tree_3d) :: initial,reference
    type(eb_geometry_3d) :: geometry(count),extended_finest
    type(eb_geometry_3d), allocatable :: extended(:)
    type(amr_tree_state_3d) :: local,lt,ir,serial_ir,old_ir,old,ot
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),sb(nv),sw(nv),sch(nv),sext(nv),limiter,slimiter
    real(dp) :: initial_inventory(nv),final_inventory(nv),error,force(3),viscosity_coefficient,recovered,sound
    integer :: owners(count),p,x,y,z,layout,failure
    character(len=16) :: faces(6),method
    logical :: good
    method='none'
    if(test_mode=='state_redist') method='state_redist'
    initial=saved;faces='outflow';force=[2.e4_dp,-3.e4_dp,4.e4_dp]
    do p=1,count
      geometry(p)=initial%nodes(p)%geometry
      do z=1,4
      do y=1,4
      do x=1,4
        if(geometry(p)%volume_fraction(x,y,z)<=0) cycle
        call reactive_conserved_to_primitive(species,initial%nodes(p)%state(:,x,y,z), &
          initial%nodes(p)%temperature(x,y,z),primitive,recovered,sound,good)
        call require(good,'SDC fixture primitive')
        primitive(5)=primitive(5)*1000._dp/recovered
        call reactive_primitive_to_conserved(species,primitive,initial%nodes(p)%state(:,x,y,z), &
          initial%nodes(p)%temperature(x,y,z),sound,good)
        call require(good,'SDC fixture EOS')
      end do
      end do
      end do
    end do
    call average_down_reactive_eb_tree_3d(species,initial,good)
    call require(good,'SDC fixture restriction')
    call build_axis_plane_eb_geometry_3d(64,64,64,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
      'z',.31_dp,extended_finest,good)
    call require(good,'SDC extended finest')
    call build_eb_tree_extended_geometry_3d(topology,geometry,extended_finest,extended,good)
    call require(good,'SDC extended geometry')
    if(test_mode=="aux_source") then
      call check_auxiliary_sources(initial,geometry,extended,faces,force)
      return
    end if
    if(test_mode=="redist_level") then
      call check_redist_level(initial,geometry,extended,faces,force)
      return
    end if
    if(test_mode=="subcycle".or.test_mode=="redist_recursive".or.test_mode=="state_recursive".or. &
     test_mode=="aux_recursive".or.test_mode=="aux_route") then
      call check_sdc_recursive(initial,geometry,extended,faces,force)
      return
    end if
    allocate(serial_ir%fields(count))
    do p=1,count
      allocate(serial_ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
    end do
    reference=initial
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,.true., &
      [.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,reference,extended,faces, &
      sb,sw,sch,sext,slimiter,good,serial_ir,redistribution=method)
    call require(good,'serial coupled SDC')
    call integrate_reactive_eb_tree_3d(initial,initial_inventory,good)
    call require(good,'SDC initial inventory')
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if(layout==2) owners=0
      allocate(local%fields(count),lt%fields(count),ir%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        local%fields(p)%values=initial%nodes(p)%state
        allocate(lt%fields(p)%values(1,4,4,4),ir%fields(p)%values(nv-4,4,4,4))
        lt%fields(p)%values(1,:,:,:)=initial%nodes(p)%temperature
        ir%fields(p)%values=0
      end do
      old=local;ot=lt;old_ir=ir
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
        1.e-10_dp,2,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,faces, &
        b,w,ch,ext,limiter,MPI_COMM_WORLD,good,ir,redistribution=method)
      call require(good,'MPI coupled SDC')
      do p=1,count
        call require(allocated(local%fields(p)%values).eqv.(owners(p)==rank),'coupled SDC owner-local state')
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==reference%nodes(p)%state),'exact serial coupled SDC state')
        call require(all(lt%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature),'exact serial coupled SDC T')
        call require(all(ir%fields(p)%values==serial_ir%fields(p)%values),'exact serial coupled SDC IR')
      end do
      call require(limiter==slimiter,'exact coupled SDC limiter')
      error=max(maxval(abs(b-sb)/max(1._dp,abs(sb))),maxval(abs(w-sw)/max(1._dp,abs(sw))), &
        maxval(abs(ch-sch)/max(1._dp,abs(sch))),maxval(abs(ext-sext)/max(1._dp,abs(sext))))
      call require(error<2.e-12_dp,'serial coupled SDC inventories')
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,final_inventory,MPI_COMM_WORLD,good)
      call require(good,'SDC final inventory')
      error=maxval(abs(final_inventory-initial_inventory-b-w-ch-ext)/max(1._dp,abs(initial_inventory)))
      if(rank==0) print *, 'MPI SDC layout, composite residual: ',layout,error
      call require(error<3.e-14_dp,'SDC composite closure')
      if(test_mode/='state_redist') call check_sdc_level(initial,geometry,extended,owners,old,ot)
      do failure=1,2
        local=old;lt=ot;ir=old_ir;b=17;w=19;ch=23;ext=29;limiter=31
        viscosity_coefficient=.1_dp
        if(failure==1) viscosity_coefficient=-1
        if(failure==2.and.rank==nranks-1) force(1)=force(1)+1
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
          1.e-10_dp,2,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,faces, &
          b,w,ch,ext,limiter,MPI_COMM_WORLD,good,ir,difmag=viscosity_coefficient,redistribution=method)
        force=[2.e4_dp,-3.e4_dp,4.e4_dp]
        call require(.not.good,'collective coupled SDC rejection')
        do p=1,count
          if(owners(p)/=rank) cycle
          call require(all(local%fields(p)%values==old%fields(p)%values),'coupled SDC state rollback')
          call require(all(lt%fields(p)%values==ot%fields(p)%values),'coupled SDC T rollback')
          call require(all(ir%fields(p)%values==old_ir%fields(p)%values),'coupled SDC IR rollback')
        end do
        call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.limiter==31, &
          'coupled SDC inventory rollback')
      end do
      deallocate(local%fields,lt%fields,ir%fields)
    end do
  end subroutine check_sdc_loop
  subroutine check_auxiliary_sources(initial,geometry,extended,faces,force)
    type(eb_reactive_tree_3d), intent(in) :: initial
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: force(3)
    type(eb_reactive_tree_3d) :: reference,profile
    type(amr_tree_state_3d) :: local,lt,ir,serial_ir,old_ir,old,ot
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),sb(nv),sw(nv),sch(nv),sext(nv),limiter,slimiter
    real(dp) :: aux(nv),saux(nv),expected_aux(nv),initial_inventory(nv),final_inventory(nv),error
    integer :: owners(count),p,x,y,z,layout,failure
    logical :: good
    serial_source_calls=0;source_failure=0
    allocate(serial_ir%fields(count))
    do p=1,count
      allocate(serial_ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
    end do
    reference=initial
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,.true., &
      [.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,reference,extended,faces, &
      sb,sw,sch,sext,slimiter,good,serial_ir,auxiliary_source=serial_auxiliary,auxiliary_change=saux)
    call require(good,'serial coupled SDC')
    call require(serial_source_calls==3,'one old and two new serial source calls')
    profile=initial
    do p=1,count
      do z=1,4
      do y=1,4
      do x=1,4
        profile%nodes(p)%state(:,x,y,z)=source_profile(p,x,y,z,0)
        if(profile%nodes(p)%geometry%volume_fraction(x,y,z)<=0) profile%nodes(p)%state(:,x,y,z)=0
      end do
      end do
      end do
    end do
    call integrate_reactive_eb_tree_3d(profile,expected_aux,good)
    call require(good,'source profile inventory')
    expected_aux=2.e-10_dp*expected_aux
    call require(maxval(abs(saux-expected_aux)/max(1.e-30_dp,abs(expected_aux)))<3.e-14_dp, &
      'accepted old/new source inventory excludes intermediate iterations')
    call integrate_reactive_eb_tree_3d(initial,initial_inventory,good)
    call require(good,'SDC initial inventory')
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if(layout==2) owners=0
      allocate(local%fields(count),lt%fields(count),ir%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        local%fields(p)%values=initial%nodes(p)%state
        allocate(lt%fields(p)%values(1,4,4,4),ir%fields(p)%values(nv-4,4,4,4))
        lt%fields(p)%values(1,:,:,:)=initial%nodes(p)%temperature
        ir%fields(p)%values=0
      end do
      old=local;ot=lt;old_ir=ir;mpi_source_calls=0;source_failure=0
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
        1.e-10_dp,2,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,faces, &
        b,w,ch,ext,limiter,MPI_COMM_WORLD,good,ir,auxiliary_source=mpi_auxiliary,auxiliary_change=aux)
      call require(good,'MPI coupled SDC')
      call require(mpi_source_calls==3,'one old and two new MPI source calls')
      call require(maxval(abs(aux-saux)/max(1.e-30_dp,abs(saux)))<3.e-14_dp,'separate serial MPI source ledger')
      do p=1,count
        call require(allocated(local%fields(p)%values).eqv.(owners(p)==rank),'coupled SDC owner-local state')
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==reference%nodes(p)%state),'exact serial coupled SDC state')
        call require(all(lt%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature),'exact serial coupled SDC T')
        call require(all(ir%fields(p)%values==serial_ir%fields(p)%values),'exact serial coupled SDC IR')
      end do
      call require(limiter==slimiter,'exact coupled SDC limiter')
      error=max(maxval(abs(b-sb)/max(1._dp,abs(sb))),maxval(abs(w-sw)/max(1._dp,abs(sw))), &
        maxval(abs(ch-sch)/max(1._dp,abs(sch))),maxval(abs(ext-sext)/max(1._dp,abs(sext))))
      call require(error<2.e-12_dp,'serial coupled SDC inventories')
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,final_inventory,MPI_COMM_WORLD,good)
      call require(good,'SDC final inventory')
      error=maxval(abs(final_inventory-initial_inventory-b-w-ch-ext-aux)/max(1._dp,abs(initial_inventory)))
      if(rank==0) print *, 'MPI SDC layout, composite residual: ',layout,error
      call require(error<3.e-14_dp,'SDC composite closure')
      do failure=1,3
        local=old;lt=ot;ir=old_ir;b=17;w=19;ch=23;ext=29;limiter=31
        source_failure=failure;mpi_source_calls=0;aux=37
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
          1.e-10_dp,2,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,faces, &
          b,w,ch,ext,limiter,MPI_COMM_WORLD,good,ir,auxiliary_source=mpi_auxiliary,auxiliary_change=aux)
        call require(.not.good,'collective late auxiliary source rejection')
        call require(mpi_source_calls==3,'failure follows two successful source stages')
        call require(all(aux==37),'auxiliary ledger rollback')
        do p=1,count
          if(owners(p)/=rank) cycle
          call require(all(local%fields(p)%values==old%fields(p)%values),'coupled SDC state rollback')
          call require(all(lt%fields(p)%values==ot%fields(p)%values),'coupled SDC T rollback')
          call require(all(ir%fields(p)%values==old_ir%fields(p)%values),'coupled SDC IR rollback')
        end do
        call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.limiter==31, &
          'coupled SDC inventory rollback')
      end do
      deallocate(local%fields,lt%fields,ir%fields)
    end do
  end subroutine check_auxiliary_sources
  subroutine serial_auxiliary(state,interval,stage,depth,rate,good)
    type(eb_reactive_tree_3d), intent(in) :: state
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: good
    integer :: p,x,y,z,d(3)
    call require(interval==1.e-10_dp.and.depth==-1,'serial source interval/depth')
    call require(stage==serial_source_calls,'serial source stage sequence')
    serial_source_calls=serial_source_calls+1
    do p=1,size(state%nodes)
      d=state%topology%nodes(p)%shape
      do z=1,d(3)
      do y=1,d(2)
      do x=1,d(1)
        rate%fields(p)%values(:,x,y,z)=source_profile(p,x,y,z,stage)
      end do
      end do
      end do
    end do
    good=.true.
  end subroutine
  subroutine mpi_auxiliary(state,temperature,interval,stage,depth,rate,comm,good)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    integer :: p,x,y,z,d(3),local_rank,ierr
    call MPI_Comm_rank(comm,local_rank,ierr)
    call require(ierr==MPI_SUCCESS.and.local_rank==rank,'source communicator')
    call require(interval==1.e-10_dp.and.depth==-1,'MPI source interval/depth')
    call require(stage==mpi_source_calls,'MPI source stage sequence')
    mpi_source_calls=mpi_source_calls+1;good=.true.
    do p=1,size(state%fields)
      if(.not.allocated(state%fields(p)%values)) cycle
      call require(allocated(temperature%fields(p)%values),'source temperature ownership')
      d=shape(temperature%fields(p)%values(1,:,:,:))
      do z=1,d(3)
      do y=1,d(2)
      do x=1,d(1)
        rate%fields(p)%values(:,x,y,z)=source_profile(p,x,y,z,stage)
      end do
      end do
      end do
    end do
    if(stage/=2) return
    if(source_failure==1.and.rank==nranks-1) good=.false.
    if(allocated(state%fields(1)%values)) then
      if(source_failure==2) rate%fields(1)%values(:,4,4,4)=ieee_value(0._dp,ieee_quiet_nan)
      if(source_failure==3) deallocate(rate%fields(1)%values)
    end if
  end subroutine
  function source_profile(p,x,y,z,stage) result(value)
    integer, intent(in) :: p,x,y,z,stage
    real(dp) :: value(nv),amount
    amount=1000._dp*real(stage+1,dp)*(p+.01_dp*x+.02_dp*y+.03_dp*z)
    value=0;value(1)=amount;value(2:4)=amount*[2._dp,-3._dp,4._dp]
    value(5)=amount*2.e6_dp;value(6)=amount
    if(tree%nodes(p)%geometry%volume_fraction(x,y,z)<=0) value=ieee_value(0._dp,ieee_quiet_nan)
  end function
  subroutine check_redist_level(initial,geometry,extended,faces,force)
    type(eb_reactive_tree_3d), intent(in) :: initial
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: force(3)
    type(eb_reactive_tree_3d) :: reference,ledger
    type(amr_tree_state_3d) :: local,lt,ir,serial_ir,transfer,serial_transfer,old,ot,old_ir,old_transfer
    type(eb_integrated_flux_3d), allocatable :: flux(:),serial_flux(:)
    integer :: owners(count),p,layout,x,y,z,c
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),sb(nv),sw(nv),sch(nv),sext(nv),limiter,slimiter
    real(dp) :: flux_inventory(nv),change(nv),error
    logical :: good
    reference=initial
    allocate(serial_ir%fields(count))
    do p=1,count
      allocate(serial_ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
    end do
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,.true., &
      [.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,reference,extended,faces, &
      sb,sw,sch,sext,slimiter,good,serial_ir,target_depth=1,integrated_flux=serial_flux, &
      redistribution='flux_redist',redistribution_transfer=serial_transfer)
    call require(good,'serial redistributed selected SDC')
    ledger=initial
    do p=1,count
      ledger%nodes(p)%state=0
      if(topology%nodes(p)%level/=1) cycle
      associate(g=>geometry(p),f=>serial_flux(p))
      do z=1,4
      do y=1,4
      do x=1,4
        if(g%volume_fraction(x,y,z)<=0) cycle
        do c=1,nv
          ledger%nodes(p)%state(c,x,y,z)=( &
            g%dy*g%dz*(g%x_face_fraction(x-1,y,z)*f%x(c,x-1,y,z)-g%x_face_fraction(x,y,z)*f%x(c,x,y,z))+ &
            g%dx*g%dz*(g%y_face_fraction(x,y-1,z)*f%y(c,x,y-1,z)-g%y_face_fraction(x,y,z)*f%y(c,x,y,z))+ &
            g%dx*g%dy*(g%z_face_fraction(x,y,z-1)*f%z(c,x,y,z-1)-g%z_face_fraction(x,y,z)*f%z(c,x,y,z)))/ &
            (g%volume_fraction(x,y,z)*g%dx*g%dy*g%dz)
        end do
      end do
      end do
      end do
      end associate
    end do
    call integrate_reactive_eb_tree_3d(ledger,flux_inventory,good)
    call require(good,'selected exported flux divergence')
    do p=1,count
      ledger%nodes(p)%state=reference%nodes(p)%state-initial%nodes(p)%state+serial_transfer%fields(p)%values
      if(topology%nodes(p)%level==1) then
        call require(all(serial_transfer%fields(p)%values==0),'selected level has no deferred transfer')
      else
        call require(all(reference%nodes(p)%state==initial%nodes(p)%state),'other gas waits for transfer')
      end if
    end do
    call integrate_reactive_eb_tree_3d(ledger,change,good)
    call require(good,'selected state plus transfer inventory')
    error=maxval(abs(change-flux_inventory-sw-sch-sext)/max(1._dp,abs(flux_inventory)))
    print *, 'Selected redistribution flux/ledger residual: ',error
    call require(error<2.e-9_dp,'selected state plus transfer closes exported fluxes')
    call require(any(serial_transfer%fields(1)%values/=0),'coarse transfer is active')
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if(layout==2) owners=0
      allocate(local%fields(count),lt%fields(count),ir%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        local%fields(p)%values=initial%nodes(p)%state
        allocate(lt%fields(p)%values(1,4,4,4),ir%fields(p)%values(nv-4,4,4,4))
        lt%fields(p)%values(1,:,:,:)=initial%nodes(p)%temperature
        ir%fields(p)%values=0
      end do
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
        1.e-10_dp,2,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,faces, &
        b,w,ch,ext,limiter,MPI_COMM_WORLD,good,ir,target_depth=1,integrated_flux=flux, &
        redistribution='flux_redist',redistribution_transfer=transfer)
      call require(good,'MPI redistributed selected SDC')
      do p=1,count
        call require(allocated(transfer%fields(p)%values).eqv.(owners(p)==rank),'owner-local stage transfer')
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==reference%nodes(p)%state),'exact selected redistributed gas')
        call require(all(lt%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature),'exact selected redistributed T')
        call require(all(ir%fields(p)%values==serial_ir%fields(p)%values),'exact selected redistributed IR')
        call require(all(transfer%fields(p)%values==serial_transfer%fields(p)%values),'exact stage transfer')
      end do
      error=max(maxval(abs(b-sb)/max(1._dp,abs(sb))),maxval(abs(w-sw)/max(1._dp,abs(sw))), &
        maxval(abs(ch-sch)/max(1._dp,abs(sch))),maxval(abs(ext-sext)/max(1._dp,abs(sext))))
      call require(error<2.e-12_dp.and.limiter==slimiter,'selected redistributed diagnostics')
      old=local;ot=lt;old_ir=ir;old_transfer=transfer;b=17;w=19;ch=23;ext=29;limiter=31
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,local,lt,extended, &
        1.e-10_dp,2,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,faces, &
        b,w,ch,ext,limiter,MPI_COMM_WORLD,good,ir,difmag=-1._dp,target_depth=1,integrated_flux=flux, &
        redistribution='flux_redist',redistribution_transfer=transfer)
      call require(.not.good,'selected redistribution failed stage')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==old%fields(p)%values),'selected redistributed gas rollback')
        call require(all(lt%fields(p)%values==ot%fields(p)%values),'selected redistributed T rollback')
        call require(all(ir%fields(p)%values==old_ir%fields(p)%values),'selected redistributed IR rollback')
        call require(all(transfer%fields(p)%values==old_transfer%fields(p)%values),'stage transfer rollback')
      end do
      call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.limiter==31, &
        'selected redistributed diagnostics rollback')
      deallocate(local%fields,lt%fields,ir%fields)
    end do
  end subroutine check_redist_level
  subroutine recursive_source_schedule(interval,stage,depth,time,calls)
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    integer, intent(inout) :: calls
    integer, parameter :: depths(7)=[0,1,2,2,1,2,2],ticks(7)=[0,0,0,1,2,2,3]
    integer :: n
    n=calls/3+1
    call require(n<=7,'recursive source interval count')
    call require(stage==mod(calls,3).and.depth==depths(n),'recursive source stage/depth order')
    call require(abs(time-2.5e-11_dp*ticks(n))<1e-25_dp,'recursive source start time')
    call require(abs(interval-1e-10_dp/2**depth)<1e-25_dp,'recursive source interval')
    calls=calls+1
  end subroutine
  subroutine serial_recursive_aux(state,interval,stage,depth,time,rate,good)
    type(eb_reactive_tree_3d), intent(in) :: state
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: good
    integer :: p,x,y,z
    call recursive_source_schedule(interval,stage,depth,time,serial_source_calls)
    do p=1,size(state%nodes)
      do z=1,4
      do y=1,4
      do x=1,4
        rate%fields(p)%values(:,x,y,z)=(1+time/1e-10_dp)*source_profile(p,x,y,z,stage)
      end do
      end do
      end do
    end do
    good=.true.
  end subroutine
  subroutine mpi_recursive_aux(state,temperature,interval,stage,depth,time,rate,comm,good)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    integer :: p,x,y,z
    call recursive_source_schedule(interval,stage,depth,time,mpi_source_calls)
    do p=1,size(state%fields)
      if(.not.allocated(state%fields(p)%values)) cycle
      do z=1,4
      do y=1,4
      do x=1,4
        rate%fields(p)%values(:,x,y,z)=(1+time/1e-10_dp)*source_profile(p,x,y,z,stage)
      end do
      end do
      end do
    end do
    good=.not.(source_failure==1.and.mpi_source_calls==21.and.rank==nranks-1)
  end subroutine
  subroutine check_sdc_recursive(initial,geometry,extended,faces,force)
    type(eb_reactive_tree_3d), intent(in) :: initial
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: force(3)
    procedure(eb_sdc_recursive_source), pointer :: serial_provider=>null()
    procedure(mpi_eb_sdc_recursive_source), pointer :: mpi_provider=>null()
    real(dp), allocatable :: aux(:),saux(:)
    real(dp) :: coefficient,expected_aux(nv),source_factor
    integer :: x,y,z
    type(eb_reactive_tree_3d) :: reference,profile
    type(amr_tree_state_3d) :: local,lt,ir,serial_ir,old,ot,old_ir
    integer(int64) :: counts(3),serial_counts(3)
    integer :: owners(count),p,layout
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),sb(nv),sw(nv),sch(nv),sext(nv),limiter,slimiter,error
    real(dp) :: initial_inventory(nv),final_inventory(nv)
    character(len=16) :: method
    logical :: good
    method="none"
    if(test_mode=="redist_recursive") method="flux_redist"
    if(test_mode=="state_recursive") method="state_redist"
    nullify(serial_provider,mpi_provider)
    if(test_mode=='aux_recursive'.or.test_mode=='aux_route') then
      serial_provider=>serial_recursive_aux;mpi_provider=>mpi_recursive_aux
      allocate(aux(nv),saux(nv))
    end if
    serial_source_calls=0;source_failure=0
    reference=initial
    allocate(serial_ir%fields(count))
    do p=1,count
      allocate(serial_ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
    end do
    call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1.e-10_dp,2,.true., &
      [.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,reference,extended,faces, &
      sb,sw,sch,sext,slimiter,good,serial_ir,level_steps=serial_counts,redistribution=method, &
      auxiliary_source=serial_provider,auxiliary_change=saux,route_auxiliary=test_mode=='aux_route')
    call require(good,'serial recursive SDC')
    if(allocated(saux)) then
      call require(serial_source_calls==21,'serial recursive source calls')
      profile=initial
      do p=1,count
        source_factor=3._dp-1._dp/2**topology%nodes(p)%level
        if(test_mode=='aux_route') source_factor=7.25_dp
        do z=1,4
        do y=1,4
        do x=1,4
          profile%nodes(p)%state(:,x,y,z)=1e-10_dp*source_factor* &
            source_profile(p,x,y,z,0)
          if(profile%nodes(p)%geometry%volume_fraction(x,y,z)<=0) profile%nodes(p)%state(:,x,y,z)=0
        end do
        end do
        end do
      end do
      call integrate_reactive_eb_tree_3d(profile,expected_aux,good)
      call require(good,'recursive expected inventory')
      call require(maxval(abs(saux-expected_aux)/max(1e-30_dp,abs(expected_aux)))<3e-14_dp, &
        'recursive sources integrate accepted stages and leaf time intervals')
    end if
    call require(all(serial_counts==[1_int64,2_int64,4_int64]),'serial level counts')
    call integrate_reactive_eb_tree_3d(initial,initial_inventory,good)
    call require(good,'recursive initial inventory')
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if(layout==2) owners=0
      allocate(local%fields(count),lt%fields(count),ir%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        local%fields(p)%values=initial%nodes(p)%state
        allocate(lt%fields(p)%values(1,4,4,4),ir%fields(p)%values(nv-4,4,4,4))
        lt%fields(p)%values(1,:,:,:)=initial%nodes(p)%temperature
        ir%fields(p)%values=0
      end do
      mpi_source_calls=0;source_failure=0
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,local,lt, &
        1.e-10_dp,2,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
        b,w,ch,ext,limiter,MPI_COMM_WORLD,good,ir,level_steps=counts,redistribution=method, &
        auxiliary_source=mpi_provider,auxiliary_change=aux,route_auxiliary=test_mode=='aux_route')
      call require(good,'MPI recursive SDC')
      if(allocated(aux)) then
        call require(mpi_source_calls==21,'MPI recursive source calls')
        call require(maxval(abs(aux-saux)/max(1e-30_dp,abs(saux)))<3e-14_dp,'recursive auxiliary serial MPI')
      end if
      call require(all(counts==serial_counts),'MPI recursive level counts')
      do p=1,count
        call require(allocated(local%fields(p)%values).eqv.(owners(p)==rank),'recursive owner-local fields')
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==reference%nodes(p)%state),'exact recursive state')
        call require(all(lt%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature),'exact recursive T')
        call require(all(ir%fields(p)%values==serial_ir%fields(p)%values),'exact recursive IR')
      end do
      call require(limiter==slimiter,'recursive limiter')
      error=max(maxval(abs(b-sb)/max(1._dp,abs(sb))),maxval(abs(w-sw)/max(1._dp,abs(sw))), &
        maxval(abs(ch-sch)/max(1._dp,abs(sch))),maxval(abs(ext-sext)/max(1._dp,abs(sext))))
      call require(error<2.e-12_dp,'recursive serial ledgers')
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,final_inventory,MPI_COMM_WORLD,good)
      call require(good,'recursive final inventory')
      if(allocated(aux)) final_inventory=final_inventory-aux
      error=maxval(abs(final_inventory-initial_inventory-b-w-ch-ext)/max(1._dp,abs(initial_inventory)))
      if(rank==0) print *, 'Recursive SDC layout, residual: ',layout,error
      call require(error<3.e-14_dp,'recursive closure')
      coefficient=-1;mpi_source_calls=0
      if(allocated(aux)) then
        coefficient=.1_dp;source_failure=1;aux=41
      end if
      old=local;ot=lt;old_ir=ir;b=17;w=19;ch=23;ext=29;limiter=31;counts=37
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,local,lt, &
        1.e-10_dp,2,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,extended,faces, &
        b,w,ch,ext,limiter,MPI_COMM_WORLD,good,ir,difmag=coefficient,level_steps=counts,redistribution=method, &
        auxiliary_source=mpi_provider,auxiliary_change=aux,route_auxiliary=test_mode=='aux_route')
      call require(.not.good,'recursive rejected stage')
      if(allocated(aux)) then
        call require(mpi_source_calls==21.and.all(aux==41),'last descendant source failure and ledger rollback')
      end if
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==old%fields(p)%values),'recursive state rollback')
        call require(all(lt%fields(p)%values==ot%fields(p)%values),'recursive T rollback')
        call require(all(ir%fields(p)%values==old_ir%fields(p)%values),'recursive IR rollback')
      end do
      call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.limiter==31.and. &
        all(counts==37),'recursive diagnostics rollback')
      deallocate(local%fields,lt%fields,ir%fields)
    end do
  end subroutine check_sdc_recursive
  subroutine check_sdc_level(initial,geometry,extended,owners,old,ot)
    type(eb_reactive_tree_3d), intent(in) :: initial
    type(eb_geometry_3d), intent(in) :: geometry(:),extended(:)
    integer, intent(in) :: owners(:)
    type(amr_tree_state_3d), intent(in) :: old,ot
    type(eb_reactive_tree_3d) :: reference,ending
    type(amr_tree_state_3d) :: state,temp,ir,serial_ir,finish,source0,source1,local0,local1,op,np,sop,snp
    type(eb_integrated_flux_3d), allocatable :: flux(:),serial_flux(:)
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),sb(nv),sw(nv),sch(nv),sext(nv),theta,stheta,error
    real(dp) :: force(3)
    integer :: p,attempt
    character(len=16) :: faces(6)
    logical :: good
    reference=initial;ending=initial;state=old;temp=ot;finish=old
    faces='outflow';force=[2.e4_dp,-3.e4_dp,4.e4_dp]
    allocate(serial_ir%fields(count),ir%fields(count),source0%fields(count),source1%fields(count), &
      local0%fields(count),local1%fields(count))
    do p=1,count
      ending%nodes(p)%state=1.01_dp*initial%nodes(p)%state
      allocate(serial_ir%fields(p)%values(nv-4,4,4,4),source=0._dp)
      if(topology%nodes(p)%level/=2) serial_ir%fields(p)%values=7
      allocate(source0%fields(p)%values(nv,4,4,4),source=0._dp)
      allocate(source1%fields(p)%values(nv,4,4,4),source=0._dp)
      source0%fields(p)%values(2,:,:,:)=1.e4_dp;source1%fields(p)%values(2,:,:,:)=2.e4_dp
      if(owners(p)/=rank) cycle
      finish%fields(p)%values=ending%nodes(p)%state
      ir%fields(p)%values=serial_ir%fields(p)%values
      local0%fields(p)%values=source0%fields(p)%values
      local1%fields(p)%values=source1%fields(p)%values
    end do
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,.true., &
      [.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,reference,extended,faces, &
      sb,sw,sch,sext,stheta,good,serial_ir,target_depth=2,ancestor_start=initial,ancestor_end=ending, &
      fractions_start=[.25_dp,.25_dp],fractions_end=[.75_dp,.75_dp],source_start=source0,source_end=source1, &
      integrated_flux=serial_flux,old_predictor=sop,new_predictor=snp)
    call require(good,'serial selected coupled SDC')
    do attempt=1,2
      if(attempt==1) then
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,state,temp,extended, &
          1.e-10_dp,2,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,faces, &
          b,w,ch,ext,theta,MPI_COMM_WORLD,good,ir,target_depth=2,ancestor_start=old,ancestor_end=finish, &
          temperature_start=ot,temperature_end=ot,fractions_start=[.25_dp,.25_dp],fractions_end=[.75_dp,.75_dp], &
          source_start=local0,source_end=local1,integrated_flux=flux,old_predictor=op,new_predictor=np)
        call require(good,'MPI selected coupled SDC')
      else
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,state,temp,extended, &
          1.e-10_dp,2,.true.,[.true.,.true.,.true.,.true.],force,2.e-7_dp,1.e-12_dp,faces, &
          b,w,ch,ext,theta,MPI_COMM_WORLD,good,ir,difmag=-1._dp,target_depth=2, &
          integrated_flux=flux,old_predictor=op,new_predictor=np)
        call require(.not.good,'selected loop rejects after old transport')
      end if
      call require(theta==stheta,'selected loop limiter retained')
      error=max(maxval(abs(b-sb)/max(1._dp,abs(sb))),maxval(abs(w-sw)/max(1._dp,abs(sw))), &
        maxval(abs(ch-sch)/max(1._dp,abs(sch))),maxval(abs(ext-sext)/max(1._dp,abs(sext))))
      call require(error<2.e-12_dp,'selected loop inventories match or roll back')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(state%fields(p)%values==reference%nodes(p)%state),'selected loop exact state')
        call require(all(temp%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature),'selected loop exact T')
        call require(all(ir%fields(p)%values==serial_ir%fields(p)%values),'selected loop exact IR and inactive retention')
        call require(all(op%fields(p)%values==sop%fields(p)%values),'selected loop exact old predictor')
        call require(all(np%fields(p)%values==snp%fields(p)%values),'selected loop exact new predictor')
        if(topology%nodes(p)%level/=2) then
          call require(.not.allocated(flux(p)%x),'selected loop inactive flux absent')
          cycle
        end if
        call require(all(flux(p)%x==serial_flux(p)%x).and.all(flux(p)%y==serial_flux(p)%y).and. &
          all(flux(p)%z==serial_flux(p)%z),'selected loop exact flux')
      end do
    end do
  end subroutine check_sdc_level
  subroutine check_sdc_reaction()
    type(amr_tree_state_3d) :: rates,ref_state,ref_temp,ref_ir,local,lt,lr,ir,old,ot,kept,kt,ki
    type(eb_geometry_3d) :: geometry(count)
    real(dp) :: tolerance
    integer :: p,layout,owners(count)
    logical :: good
    allocate(rates%fields(count),ref_ir%fields(count))
    do p=1,count
      geometry(p)=saved%nodes(p)%geometry
      allocate(rates%fields(p)%values(nv,4,4,4),source=0._dp)
      allocate(ref_ir%fields(p)%values(nv-4,4,4,4),source=7._dp)
      rates%fields(p)%values(2,:,:,:)=2.e4_dp
      rates%fields(p)%values(5,:,:,:)=2.e4_dp*saved%nodes(p)%state(2,:,:,:)/saved%nodes(p)%state(1,:,:,:)
    end do
    call react_sdc_eb_tree_3d(species,reactions,saved,rates,1.e-9_dp,2.e-7_dp,1.e-12_dp, &
      ref_state,ref_temp,ref_ir,good)
    call require(good,'serial SDC reaction reference')
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if(layout==2) owners=0
      allocate(old%fields(count),ot%fields(count),lr%fields(count),ir%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        old%fields(p)%values=saved%nodes(p)%state
        allocate(ot%fields(p)%values(1,4,4,4))
        ot%fields(p)%values(1,:,:,:)=saved%nodes(p)%temperature
        lr%fields(p)%values=rates%fields(p)%values
        allocate(ir%fields(p)%values(nv-4,4,4,4),source=7._dp)
      end do
      call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
        1.e-9_dp,2.e-7_dp,1.e-12_dp,local,lt,ir,MPI_COMM_WORLD,good)
      call require(good,'MPI SDC reaction')
      do p=1,count
        call require(allocated(local%fields(p)%values).eqv.(owners(p)==rank),'SDC owner-local allocation')
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==ref_state%fields(p)%values).and. &
          all(lt%fields(p)%values==ref_temp%fields(p)%values).and.all(ir%fields(p)%values==ref_ir%fields(p)%values), &
          'MPI SDC exact serial state temperature and reaction')
        local%fields(p)%values=-17;lt%fields(p)%values=-19;ir%fields(p)%values=7
      end do
      call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
        1.e-9_dp,2.e-7_dp,1.e-12_dp,local,lt,ir,MPI_COMM_WORLD,good,target_depth=2)
      call require(good,'MPI selected level SDC reaction')
      do p=1,count
        if(owners(p)/=rank) cycle
        if(topology%nodes(p)%level==2) then
          call require(all(local%fields(p)%values==ref_state%fields(p)%values).and. &
            all(ir%fields(p)%values==ref_ir%fields(p)%values),'MPI selected SDC exact reaction')
        else
          call require(all(local%fields(p)%values==-17).and.all(lt%fields(p)%values==-19).and. &
            all(ir%fields(p)%values==7),'MPI SDC inactive level retained')
        end if
      end do
      kept=local;kt=lt;ki=ir
      if(owners(count)==rank) old%fields(count)%values(1,4,4,4)=-1
      call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
        1.e-9_dp,2.e-7_dp,1.e-12_dp,local,lt,ir,MPI_COMM_WORLD,good)
      call require(.not.good,'MPI failing owner rejects globally')
      if(owners(count)==rank) old%fields(count)%values=saved%nodes(count)%state
      tolerance=2.e-7_dp
      if(rank==nranks-1) tolerance=-1
      call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
        1.e-9_dp,tolerance,1.e-12_dp,local,lt,ir,MPI_COMM_WORLD,good)
      call require(.not.good,'MPI invalid control including idle rank rejects')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==kept%fields(p)%values).and. &
          all(lt%fields(p)%values==kt%fields(p)%values).and.all(ir%fields(p)%values==ki%fields(p)%values), &
          'MPI SDC collective rollback')
      end do
      deallocate(old%fields,ot%fields,lr%fields,ir%fields)
    end do
  end subroutine check_sdc_reaction

  subroutine check_mpi_full()
    type(amr_tree_state_3d) :: local,lt,initial,it
    type(eb_geometry_3d) :: geometry(count)
    real(dp) :: b(nv),w(nv),chem(nv),limiter,inventory(nv),error
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
    call advance_mpi_eb_tree_full_3d(species,reactions,transport,topology,owners,geometry,'hllc',1.e-8_dp, &
      chemistry_on,2.e-7_dp,1.e-12_dp,transport_on,.true.,.true.,.true.,.true., &
      local,lt,b,w,chem,limiter,MPI_COMM_WORLD,good,integrators(integrator_index))
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
        local,lt,b,w,chem,limiter,MPI_COMM_WORLD,good,'explicit')
      call require(.not.good,'late hydro rejection after reaction and transport')
      do p=1,count
        if (owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==initial%fields(p)%values),'full U rollback')
        call require(all(lt%fields(p)%values==it%fields(p)%values),'full T rollback')
      end do
      call require(all(b==17._dp).and.all(w==19._dp).and.all(chem==23._dp).and.limiter==29._dp,'full ledger rollback')
      enabled=rank==0
      call advance_mpi_eb_tree_full_3d(species,reactions,transport,topology,owners,geometry,'hllc',1.e-8_dp, &
        enabled,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true., &
        local,lt,b,w,chem,limiter,MPI_COMM_WORLD,good,'explicit')
      call require(.not.good,'different chemistry flags rejected')
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
end program test_mpi_eb_tree_full_3d
