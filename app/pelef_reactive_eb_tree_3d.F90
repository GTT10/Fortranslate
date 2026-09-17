program pelef_reactive_eb_tree_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64,iostat_end,iostat_eor
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use reactive_1d_mod, only: reactive_nvar
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,initialize_amr_patch_tree_3d,amr_tree_integer_layout_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config,build_configured_eb_geometry_3d
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    advance_reactive_eb_tree_full_3d,integrate_reactive_eb_tree_3d,reactive_eb_tree_timestep_3d, &
    remap_reactive_eb_tree_3d
  use eb_reactive_tree_checkpoint_3d_mod, only: eb_tree_checkpoint_3d,reactive_eb_tree_checkpoint_record_3d
  implicit none
  integer :: root_shape(3)=[4,4,4],child_count=0,parent(32)=1,ratio(32)=2,lower(3,32)=1,upper(3,32)=2
  integer :: regrid_count=0,regrid_stage=0,regrid_child_count(16)=0
  integer :: regrid_parent(32,16)=1,regrid_ratio(32,16)=2,regrid_lower(3,32,16)=1,regrid_upper(3,32,16)=2
  real(dp) :: regrid_times(16)=0._dp,next_time
  integer :: maximum_steps=100,maximum_step_retries=12,checkpoint_interval_steps=0
  real(dp) :: domain_lower(3)=0._dp,domain_upper(3)=1._dp,plane_position=.31_dp
  logical :: automatic_timestep=.false.
  real(dp) :: hydro_cfl=.5_dp,transport_cfl=.35_dp
  real(dp) :: time_step=1.e-9_dp,final_time=2.e-9_dp,rtol=2.e-7_dp,atol=1.e-12_dp
  logical :: chemistry_enabled=.true.,transport_enabled=.true.,stop_after_checkpoint=.false.
  character(len=16) :: geometry_type='plane',plane_axis='z',extrusion_axis='z'
  character(len=16) :: riemann_solver='rusanov',chemistry_integrator='explicit'
  character(len=512) :: level_set_file='',section_file='',output_file='eb_tree.csv',checkpoint_file='',restart_file=''
  namelist /eb_tree/ root_shape,child_count,parent,ratio,lower,upper,domain_lower,domain_upper, &
    geometry_type,plane_axis,plane_position,extrusion_axis,level_set_file,section_file,time_step,final_time, &
    maximum_steps,maximum_step_retries,rtol,atol,chemistry_enabled,transport_enabled,riemann_solver, &
    automatic_timestep,hydro_cfl,transport_cfl,chemistry_integrator,output_file,checkpoint_file,restart_file, &
    checkpoint_interval_steps,stop_after_checkpoint,regrid_count,regrid_times,regrid_child_count, &
    regrid_parent,regrid_ratio,regrid_lower,regrid_upper
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_tree_3d) :: topology
  type(amr_patch_tree_3d), allocatable :: layouts(:)
  type(eb_geometry_3d), allocatable :: geometry_levels(:)
  type(eb_reactive_tree_3d) :: tree
  type(eb_tree_checkpoint_3d) :: record
  type(eb_geometry_3d) :: finest,coarse
  type(reactive_eb_3d_config) :: config,root_config
  integer(int64), allocatable :: origin(:,:),domain(:,:)
  integer, allocatable :: ratios(:),shared_ratios(:)
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:),current(:),boundary(:),wall(:),chemical(:)
  real(dp) :: dt,theta,error,estimated_dt
  integer :: unit,ios,status,node,n(3),nv,retry,depth,event,common
  logical :: ok,stopped=.false.
  character(len=1024) :: input_path,message,budget_path
  if (command_argument_count()/=1) error stop 'Usage: pelef_reactive_eb_tree_3d input.nml'
  call get_command_argument(1,input_path)
  open(newunit=unit,file=trim(input_path),status='old',action='read',iostat=ios)
  if (ios/=0) error stop 'Cannot open input'
  read(unit,nml=eb_tree,iostat=ios)
  close(unit,iostat=status)
  if (ios/=0.or.status/=0) error stop 'Invalid eb_tree input'
  if (child_count<0.or.child_count>32.or.any(root_shape<2)) error stop 'Invalid tree dimensions'
  if (maximum_steps<1.or.maximum_step_retries<0.or.maximum_step_retries>60) error stop 'Invalid step limits'
  if (.not.all(ieee_is_finite([domain_lower,domain_upper,plane_position,time_step,final_time,rtol,atol,hydro_cfl,transport_cfl]))) &
    error stop 'Nonfinite controls'
  if (any(domain_upper<=domain_lower).or.min(time_step,final_time,rtol,atol)<=0._dp) error stop 'Invalid controls'
  if (hydro_cfl<=0._dp.or.hydro_cfl>1._dp.or.transport_cfl<=0._dp.or.transport_cfl>.5_dp) &
    error stop 'Invalid CFL controls'
  if (chemistry_integrator/='explicit'.and.chemistry_integrator/='implicit') error stop 'Invalid chemistry integrator'
  if (riemann_solver/='rusanov'.and.riemann_solver/='hllc'.and.riemann_solver/='pelec') error stop 'Invalid solver'
  if (len_trim(output_file)==0.or.checkpoint_interval_steps<0) error stop 'Invalid output policy'
  if (checkpoint_interval_steps>0.and.len_trim(checkpoint_file)==0) error stop 'Missing checkpoint path'
  if (stop_after_checkpoint.and.checkpoint_interval_steps==0) error stop 'Missing checkpoint interval'
  budget_path=trim(output_file)//'.budget.csv'
  if (output_file==input_path.or.budget_path==input_path.or.output_file==restart_file.or.budget_path==restart_file.or. &
    output_file==level_set_file.or.budget_path==level_set_file.or.output_file==section_file.or.budget_path==section_file) &
    error stop 'Output aliases input'
  if (len_trim(checkpoint_file)>0) then
    if (checkpoint_file==input_path.or.checkpoint_file==output_file.or.checkpoint_file==budget_path.or. &
      checkpoint_file==restart_file.or.checkpoint_file==level_set_file.or.checkpoint_file==section_file) &
      error stop 'Checkpoint aliases another file'
  end if
  if (regrid_count<0.or.regrid_count>16) error stop 'Invalid regrid count'
  allocate(layouts(0:regrid_count))
  call make_layout(child_count,parent,ratio,lower,upper,layouts(0))
  do event=1,regrid_count
    if (.not.ieee_is_finite(regrid_times(event))) error stop 'Nonfinite regrid time'
    if (regrid_times(event)<=0._dp) error stop 'Invalid regrid time'
    if (event>1) then
      if (regrid_times(event)<=regrid_times(event-1)) error stop 'Regrid times must increase'
    end if
    call make_layout(regrid_child_count(event),regrid_parent(:,event),regrid_ratio(:,event), &
      regrid_lower(:,:,event),regrid_upper(:,:,event),layouts(event))
  end do
  do event=0,regrid_count
    call amr_tree_integer_layout_3d(layouts(event),origin,domain,ratios,ok)
    if (.not.ok) error stop 'Invalid integer layout'
    if (.not.allocated(shared_ratios)) shared_ratios=ratios
    common=min(size(shared_ratios),size(ratios))
    if (any(shared_ratios(:common)/=ratios(:common))) error stop 'Schedule changes refinement lattice'
    if (size(ratios)>size(shared_ratios)) shared_ratios=ratios
  end do
  n=root_shape
  do depth=2,size(shared_ratios)
    if (any(n>huge(1)/shared_ratios(depth))) error stop 'Finest lattice too large'
    n=n*shared_ratios(depth)
  end do
  topology=layouts(0)
  config%nx=n(1);config%ny=n(2);config%nz=n(3)
  config%x_lower=domain_lower(1);config%y_lower=domain_lower(2);config%z_lower=domain_lower(3)
  config%x_upper=domain_upper(1);config%y_upper=domain_upper(2);config%z_upper=domain_upper(3)
  config%geometry_type=geometry_type;config%plane_axis=plane_axis;config%plane_position=plane_position
  config%extrusion_axis=extrusion_axis;config%section_file=section_file;config%level_set_file=level_set_file
  config%cfl=hydro_cfl;config%transport_cfl=transport_cfl
  config%redistribution='none';config%riemann_solver=riemann_solver
  config%chemistry_enabled=chemistry_enabled;config%transport_enabled=transport_enabled
  config%chemistry_relative_tolerance=rtol;config%chemistry_absolute_tolerance=atol
  config%viscosity_enabled=.true.;config%thermal_conduction_enabled=.true.
  config%species_diffusion_enabled=.true.;config%barodiffusion_enabled=.true.
  config%final_time=final_time;config%maximum_steps=maximum_steps
  call build_configured_eb_geometry_3d(config,finest,ok)
  if (.not.ok) error stop 'Cannot build finest geometry'
  allocate(geometry_levels(size(shared_ratios)));geometry_levels(size(shared_ratios))=finest
  do depth=size(shared_ratios),2,-1
    call coarsen_eb_geometry_3d(geometry_levels(depth),shared_ratios(depth),geometry_levels(depth-1),ok)
    if (.not.ok) error stop 'Cannot coarsen geometry'
  end do
  coarse=geometry_levels(1)
  call load_h2o2_elementary_thermo(species,ok)
  if (.not.ok) error stop 'Thermo load failed'
  call load_h2o2_elementary_mechanism(reactions,ok)
  if (.not.ok) error stop 'Reaction load failed'
  call load_h2o2_elementary_transport(transport,ok)
  if (.not.ok) error stop 'Transport load failed'
  nv=reactive_nvar(size(species));n=root_shape
  allocate(u(nv,n(1),n(2),n(3)),t(n(1),n(2),n(3)),current(nv),boundary(nv),wall(nv),chemical(nv))
  root_config=config;root_config%nx=n(1);root_config%ny=n(2);root_config%nz=n(3)
  call initialize_reactive_eb_density_sheet_3d(species,root_config,coarse,u,t,ok)
  if (.not.ok) error stop 'Initial state failed'
  call initialize_reactive_eb_tree_3d(species,topology,geometry_levels(topology%level_count()),u,t,tree,ok)
  if (.not.ok) error stop 'Tree initialization failed'
  call integrate_reactive_eb_tree_3d(tree,current,ok)
  if (.not.ok) error stop 'Initial inventory failed'
  record%initial=current;record%scale=abs(current);record%boundary=0._dp*current
  record%wall=0._dp*current;record%chemical=0._dp*current
  if (len_trim(restart_file)>0) call restore()
  do while(record%time<final_time.and.record%steps<maximum_steps)
    dt=min(time_step,final_time-record%time)
    if (automatic_timestep) then
      call reactive_eb_tree_timestep_3d(species,transport,tree,hydro_cfl,transport_cfl,transport_enabled, &
        .true.,.true.,.true.,estimated_dt,ok)
      if (.not.ok) error stop 'Tree timestep estimate failed'
      dt=min(dt,estimated_dt)
    end if
    if (regrid_stage<regrid_count) dt=min(dt,regrid_times(regrid_stage+1)-record%time)
    do retry=0,maximum_step_retries
      if (record%time+dt<=record%time) error stop 'Time increment underflow'
      call advance_reactive_eb_tree_full_3d(species,reactions,transport,riemann_solver,dt,chemistry_enabled, &
        rtol,atol,transport_enabled,.true.,.true.,.true.,.true.,tree,boundary,wall,chemical,theta,ok,chemistry_integrator)
      if (ok) exit
      if (retry==maximum_step_retries) error stop 'Tree step retry limit exhausted'
      dt=.5_dp*dt
    end do
    record%boundary=record%boundary+boundary;record%wall=record%wall+wall;record%chemical=record%chemical+chemical
    next_time=record%time+dt
    if (regrid_stage<regrid_count) then
      if (dt==regrid_times(regrid_stage+1)-record%time) next_time=regrid_times(regrid_stage+1)
    end if
    record%time=next_time;record%steps=record%steps+1
    if (record%steps==1) record%minimum_dt=dt
    record%minimum_dt=min(record%minimum_dt,dt);record%minimum_theta=min(record%minimum_theta,theta)
    if (regrid_stage<regrid_count) then
      if (record%time>=regrid_times(regrid_stage+1)) then
        call activate(regrid_stage+1,.true.)
        write(*,'(a,i0,a,i0)') 'Regrid stage: ',regrid_stage,' nodes: ',size(tree%nodes)
      end if
    end if
    call balance()
    if (checkpoint_interval_steps>0) then
      if (mod(record%steps,checkpoint_interval_steps)==0) then
        call save()
        if (stop_after_checkpoint) then
          stopped=.true.
          exit
        end if
      end if
    end if
  end do
  if (record%time<final_time.and..not.stopped) error stop 'Maximum steps reached before final time'
  call balance()
  call output()
  write(*,'(a,i0,a,es24.16)') 'Completed steps: ',record%steps,' time: ',record%time
  write(*,'(a,es24.16)') 'Composite inventory error: ',error
  write(*,'(a,es24.16)') 'Minimum accepted interval: ',record%minimum_dt
contains
  subroutine make_layout(number,parents,refinement,lo,hi,result)
    integer, intent(in) :: number,parents(:),refinement(:),lo(:,:),hi(:,:)
    type(amr_patch_tree_3d), intent(out) :: result
    type(amr_patch_3d), allocatable :: patches(:)
    integer :: node,parent_node,dims(3)
    logical :: valid
    if (number<0.or.number>32) error stop 'Invalid scheduled child count'
    allocate(patches(number))
    do node=1,number
      parent_node=parents(node)
      if (parent_node<1.or.parent_node>node) error stop 'Parents must precede children'
      dims=root_shape
      if (parent_node>1) dims=(hi(:,parent_node-1)-lo(:,parent_node-1)+1)*refinement(parent_node-1)
      call initialize_amr_patch_3d(dims(1),dims(2),dims(3),lo(1,node),hi(1,node),lo(2,node),hi(2,node), &
        lo(3,node),hi(3,node),refinement(node),patches(node),valid)
      if (.not.valid) error stop 'Invalid scheduled box'
    end do
    call initialize_amr_patch_tree_3d(root_shape,parents(:number),patches,result,valid)
    if (.not.valid) error stop 'Invalid scheduled topology'
  end subroutine
  subroutine activate(stage,remapping)
    integer, intent(in) :: stage
    logical, intent(in) :: remapping
    integer :: level
    logical :: valid
    if (remapping) then
      level=max(topology%level_count(),layouts(stage)%level_count())
      call remap_reactive_eb_tree_3d(species,layouts(stage),geometry_levels(level),tree,valid)
    else
      level=layouts(stage)%level_count()
      call initialize_reactive_eb_tree_3d(species,layouts(stage),geometry_levels(level),u,t,tree,valid)
    end if
    if (.not.valid) error stop 'Cannot activate scheduled tree'
    topology=layouts(stage);regrid_stage=stage
  end subroutine
  subroutine schedule_record(stream,reading)
    integer, intent(in) :: stream
    logical, intent(in) :: reading
    integer :: event,node,stored_count,stored_stage,node_count,descriptor(12),stored(12)
    real(dp) :: event_time
    if (reading) then
      read(stream,*,iostat=ios) stored_count,stored_stage
      if (ios/=0) error stop 'Missing restart schedule'
      if (stored_count/=regrid_count.or.stored_stage<0.or.stored_stage>regrid_count) &
        error stop 'Invalid restart schedule count or stage'
    else
      write(stream,*,iostat=ios) regrid_count,regrid_stage
      if (ios/=0) error stop 'Cannot write schedule count'
    end if
    do event=0,regrid_count
      if (event>0) then
        if (reading) then
          read(stream,*,iostat=ios) event_time
          if (ios/=0) error stop 'Missing event time'
          if (event_time/=regrid_times(event)) error stop 'Restart event time mismatch'
        else
          write(stream,'(es27.18e3)',iostat=ios) regrid_times(event)
          if (ios/=0) error stop 'Cannot write event time'
        end if
      end if
      node_count=size(layouts(event)%nodes)
      if (reading) then
        read(stream,*,iostat=ios) stored_count
        if (ios/=0) error stop 'Missing event topology'
        if (stored_count/=node_count) error stop 'Restart event node count mismatch'
      else
        write(stream,*,iostat=ios) node_count
        if (ios/=0) error stop 'Cannot write node count'
      end if
      do node=1,node_count
        associate(n=>layouts(event)%nodes(node),p=>layouts(event)%nodes(node)%patch)
          descriptor=[n%parent,n%level,n%shape,p%refinement_ratio,p%coarse_i_lower,p%coarse_i_upper, &
            p%coarse_j_lower,p%coarse_j_upper,p%coarse_k_lower,p%coarse_k_upper]
        end associate
        if (reading) then
          read(stream,*,iostat=ios) stored
          if (ios/=0) error stop 'Missing event node'
          if (any(stored/=descriptor)) error stop 'Restart event topology mismatch'
        else
          write(stream,*,iostat=ios) descriptor
          if (ios/=0) error stop 'Cannot write event node'
        end if
      end do
    end do
    if (reading) call activate(stored_stage,.false.)
  end subroutine
  subroutine balance()
    call integrate_reactive_eb_tree_3d(tree,current,ok)
    if (.not.ok) error stop 'Inventory failed'
    error=maxval(abs(current-record%initial-record%boundary-record%wall-record%chemical)/ &
      max(1._dp,record%scale,abs(current),abs(record%boundary)+abs(record%wall)+abs(record%chemical)))
    if (.not.ieee_is_finite(error).or.error>5.e-10_dp) error stop 'Composite budget failed'
  end subroutine
  subroutine envelope(stream,reading)
    integer, intent(in) :: stream
    logical, intent(in) :: reading
    character(len=64) :: magic,integrator
    logical :: automatic
    real(dp) :: interval,stored_hydro_cfl,stored_transport_cfl
    if (reading) then
      read(stream,'(a)',iostat=ios) magic
      if (ios/=0) error stop 'Missing tree application header'
      if (magic/='PELEF_EB_TREE_APPLICATION_V1'.and.magic/='PELEF_EB_TREE_APPLICATION_V2'.and. &
        magic/='PELEF_EB_TREE_APPLICATION_V3') &
        error stop 'Wrong tree application header'
      if (magic=='PELEF_EB_TREE_APPLICATION_V1'.and.automatic_timestep) error stop 'Restart policy mismatch'
      read(stream,*,iostat=ios) interval,integrator
      if (ios/=0) error stop 'Missing restart policy'
      if (interval/=time_step.or.integrator/=chemistry_integrator) error stop 'Restart policy mismatch'
      if (magic/='PELEF_EB_TREE_APPLICATION_V1') then
        read(stream,*,iostat=ios) automatic,stored_hydro_cfl,stored_transport_cfl
        if (ios/=0) error stop 'Missing automatic interval policy'
        if (automatic.neqv.automatic_timestep) error stop 'Restart interval mode mismatch'
        if (stored_hydro_cfl/=hydro_cfl.or.stored_transport_cfl/=transport_cfl) error stop 'Restart CFL mismatch'
      end if
    else
      magic='PELEF_EB_TREE_APPLICATION_V2'
      if (regrid_count>0) magic='PELEF_EB_TREE_APPLICATION_V3'
      write(stream,'(a)',iostat=ios) trim(magic)
      if (ios/=0) error stop 'Cannot write application header'
      write(stream,'(es27.18e3,1x,a)',iostat=ios) time_step,trim(chemistry_integrator)
      if (ios/=0) error stop 'Cannot write policy'
      write(stream,'(l1,2(1x,es27.18e3))',iostat=ios) automatic_timestep,hydro_cfl,transport_cfl
      if (ios/=0) error stop 'Cannot write automatic interval policy'
    end if
    if (magic=='PELEF_EB_TREE_APPLICATION_V3') then
      call schedule_record(stream,reading)
    else if (regrid_count/=0) then
      error stop 'Restart schedule missing'
    end if
    call reactive_eb_tree_checkpoint_record_3d(stream,reading,species,reactions,transport,config,tree,record,ok,message)
    if (.not.ok) then
      print *,trim(message)
      error stop 'Invalid tree checkpoint record'
    end if
    if (reading) then
      if (regrid_stage/=count(regrid_times(:regrid_count)<=record%time)) error stop 'Restart stage/time mismatch'
    end if
  end subroutine
  subroutine restore()
    character(len=1) :: extra
    open(newunit=unit,file=trim(restart_file),status='old',action='read',iostat=ios)
    if (ios/=0) error stop 'Cannot open restart'
    call envelope(unit,.true.)
    read(unit,'(a)',iostat=ios) extra
    if (ios/=iostat_end) error stop 'Trailing restart data'
    close(unit,iostat=ios)
    if (ios/=0) error stop 'Restart close failed'
    if (record%time>final_time.or.record%steps>maximum_steps) error stop 'Restart beyond requested limits'
    call balance()
  end subroutine
  subroutine save()
    integer :: scratch,destination,amount
    character(len=4096) :: chunk
    open(newunit=scratch,status='scratch',form='formatted',action='readwrite',iostat=ios)
    if (ios/=0) error stop 'Cannot stage checkpoint'
    call envelope(scratch,.false.)
    rewind(scratch,iostat=ios)
    if (ios/=0) error stop 'Checkpoint rewind failed'
    open(newunit=destination,file=trim(checkpoint_file),status='replace',action='write',iostat=ios)
    if (ios/=0) error stop 'Cannot open checkpoint destination'
    do
      read(scratch,'(a)',advance='no',size=amount,iostat=ios) chunk
      if (ios==iostat_end) exit
      if (ios/=0.and.ios/=iostat_eor) error stop 'Checkpoint staging read failed'
      write(destination,'(a)',advance='no',iostat=status) chunk(:amount)
      if (status/=0) error stop 'Checkpoint publication failed'
      if (ios==iostat_eor) then
        write(destination,'(a)',iostat=status) ''
        if (status/=0) error stop 'Checkpoint newline failed'
      end if
    end do
    close(destination,iostat=ios)
    if (ios/=0) error stop 'Checkpoint close failed'
    close(scratch,iostat=ios)
    if (ios/=0) error stop 'Checkpoint staging close failed'
  end subroutine
  subroutine output()
    integer :: i,j,k,c,child,lo(3),hi(3)
    logical, allocatable :: exposed(:,:,:)
    open(newunit=unit,file=trim(output_file),status='replace',action='write',iostat=ios)
    if (ios/=0) error stop 'Cannot open CSV'
    write(unit,'(a)',advance='no') 'node,level,i,j,k,fluid_volume,temperature,rho,rhou,rhov,rhow,rhoE'
    do c=1,size(species)
      write(unit,'(a)',advance='no') ',rhoY_'//trim(species(c)%name)
    end do
    write(unit,*)
    do node=1,size(tree%nodes)
      n=topology%nodes(node)%shape
      allocate(exposed(n(1),n(2),n(3)));exposed=tree%nodes(node)%geometry%volume_fraction>0._dp
      do child=node+1,size(tree%nodes)
        if (topology%nodes(child)%parent/=node) cycle
        associate(box=>topology%nodes(child)%patch)
          lo=[box%coarse_i_lower,box%coarse_j_lower,box%coarse_k_lower]
          hi=[box%coarse_i_upper,box%coarse_j_upper,box%coarse_k_upper]
        end associate
        exposed(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.false.
      end do
      associate(g=>tree%nodes(node)%geometry)
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if (.not.exposed(i,j,k)) cycle
        write(unit,'(i0,4(",",i0),*(:,",",es26.17e3))') node,topology%nodes(node)%level,i,j,k, &
          g%volume_fraction(i,j,k)*g%dx*g%dy*g%dz,tree%nodes(node)%temperature(i,j,k),tree%nodes(node)%state(:,i,j,k)
      end do
      end do
      end do
      end associate
      deallocate(exposed)
    end do
    close(unit,iostat=ios)
    if (ios/=0) error stop 'CSV close failed'
    open(newunit=unit,file=trim(budget_path),status='replace',action='write',iostat=ios)
    if (ios/=0) error stop 'Cannot open budget CSV'
    write(unit,'(a)') 'component,initial,final,boundary,wall,chemistry,residual'
    do c=1,nv
      write(unit,'(i0,6(",",es26.17e3))') c,record%initial(c),current(c),record%boundary(c), &
        record%wall(c),record%chemical(c),current(c)-record%initial(c)-record%boundary(c)-record%wall(c)-record%chemical(c)
    end do
    close(unit,iostat=ios)
    if (ios/=0) error stop 'Budget close failed'
  end subroutine
end program pelef_reactive_eb_tree_3d
