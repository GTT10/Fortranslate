module mpi_amr_patch_tree_application_3d_mod
  use amr_patch_tree_checkpoint_3d_mod, only: amr_tree_checkpoint_3d,build_amr_tree_checkpoint_context_3d
  use mpi_amr_patch_tree_remap_3d_mod, only: remap_mpi_amr_tree_reactive_3d,remap_mpi_amr_tree_3d
  use amr_patch_tree_remap_3d_mod, only: cell_covered=>covered
  use mpi_amr_patch_tree_tagging_3d_mod, only: tag_mpi_amr_tree_temperature_3d
  use mpi_f08
  use, intrinsic :: iso_c_binding, only: c_char,c_size_t,c_ptr,c_associated,c_null_char
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use simulation_config_selected_reactor_mod, only: selected_reactor_paths_alias_lexically
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use mixture_thermo_mod, only: mixture_density,mass_fractions_from_mole_fractions
  use selected_composition_mod, only: resolve_selected_composition
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved, &
    resolve_reactive_chemistry_integrator
  use reactive_3d_mod, only: compute_reactive_cfl_timestep_3d
  use mpi_amr_patch_tree_ghost_3d_mod, only: fill_mpi_amr_tree_ghosts_3d
  use reactive_transport_3d_mod, only: reactive_transport_timestep_3d
  use simulation_config_reactive_3d_mod, only: reactive_3d_config,resolved_reactive_boundaries_3d
  use reactive_entropy_wave_3d_problem_mod, only: initialize_reactive_problem_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    amr_tree_integer_layout_3d
  use mpi_amr_patch_tree_thermo_3d_mod, only: compute_mpi_amr_tree_cfl_3d,recover_mpi_amr_tree_temperatures_3d
  use mpi_amr_patch_tree_transport_3d_mod, only: compute_mpi_amr_tree_transport_timestep_3d
  use mpi_amr_patch_tree_subcycle_3d_mod, only: compute_mpi_amr_tree_subcycle_cfl_3d
  use mpi_amr_patch_tree_transport_subcycle_3d_mod, only: compute_mpi_amr_tree_subcycle_transport_timestep_3d
  use mpi_amr_patch_tree_split_3d_mod, only: advance_mpi_amr_tree_split_3d,advance_mpi_amr_tree_sdc_3d, &
    advance_mpi_amr_tree_sdc_subcycle_3d
  use mpi_amr_patch_tree_restriction_3d_mod, only: average_down_mpi_amr_tree_3d
  use mpi_amr_patch_tree_distribution_3d_mod, only: &
    partition_amr_tree_3d,rebalance_mpi_amr_tree_3d
  use mpi_amr_patch_tree_checkpoint_3d_mod, only: &
    write_mpi_amr_tree_checkpoint_3d,read_mpi_amr_tree_checkpoint_3d
  use mpi_amr_patch_tree_csv_3d_mod, only: write_mpi_amr_tree_csv_3d
  use mpi_amr_tree_work_mod, only: begin_tree_work,end_tree_work
  use simulation_config_reactive_3d_mod, only: resolved_reactive_difmag_3d
  implicit none
  private
  public :: run_mpi_amr_patch_tree_application_3d
  interface
    function system_getcwd(buffer,length) bind(C,name='getcwd') result(pointer)
      import c_char,c_size_t,c_ptr
      character(kind=c_char), intent(out) :: buffer(*)
      integer(c_size_t), value :: length
      type(c_ptr) :: pointer
    end function
  end interface
contains
  subroutine run_mpi_amr_patch_tree_application_3d(topology_path,input_path,config,species,reactions,transport, &
      mole_fractions,comm,chemistry_integrator,mechanism_identity,schedule_path)
    character(len=*), intent(in) :: topology_path,input_path
    character(len=*), optional, intent(in) :: chemistry_integrator,mechanism_identity,schedule_path
    type(reactive_3d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: mole_fractions(:)
    type(MPI_Comm), intent(in) :: comm
    type(reactive_3d_config) :: local
    type(amr_patch_tree_3d) :: tree,tagged_tree
    type(amr_patch_tree_3d), allocatable :: events(:)
    real(dp), allocatable :: event_times(:),source_budget(:),source_change(:)
    real(dp) :: max_remap_error
    integer :: event,nevents,slash,regrids,last_checkpoint_step
    character(len=1024) :: event_path
    type(amr_tree_state_3d) :: state,temperature
    type(amr_tree_state_3d), allocatable :: reaction_source
    integer, allocatable :: ratios(:),owners(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    real(dp), allocatable :: x(:),y(:),z(:),root_h(:),h(:,:),mass_fractions(:),inlet_x(:)
    real(dp), allocatable :: initial(:),final(:),budget(:),change(:),inlet(:,:),q(:),inlet_cube(:,:,:,:)
    real(dp) :: time,dt,transport_dt,diffusivity,base_density,correction,theta,max_correction,min_theta,stop_time
    real(dp) :: inlet_t(6),inlet_temperature(2,2,2),inlet_dt,face_dt,sound,min_h(3),mach_pressure(6)
    integer :: count,unit,status,node,n(3),i,j,k,step,attempt,face,c,rank,nranks
    character(len=32) :: faces(6),partitioning
    character(len=1024) :: message
    character(len=1024) :: checkpoint_file,restart_file
    character(len=:), allocatable :: context,identity,integrator,schedule_identity,expected_tree
    integer, parameter :: maximum_tag_levels=16
    integer :: regrid_interval_steps,refinement_levels,refinement_ratios(maximum_tag_levels)
    integer :: tile_shape(3),tag_buffer,maximum_nodes,last_tag_step
    real(dp) :: temperature_thresholds(maximum_tag_levels),tag_spacing(3)
    logical :: periodic(3)
    integer :: rebalance_interval_steps
    real(dp) :: rebalance_minimum_gain
    real(dp), allocatable :: measured_cost(:)
    logical :: sample_work,redistributed
    integer :: checkpoint_interval_steps
    logical :: stop_after_checkpoint,use_implicit,subcycling
    logical :: ok,viscosity,conduction,diffusion,barodiffusion
    logical :: checkpoint_indexed,sdc
    namelist /amr_tree_checkpoint/ checkpoint_file,restart_file,checkpoint_interval_steps,stop_after_checkpoint,checkpoint_indexed
    namelist /amr_tree_adaptation/ regrid_interval_steps,refinement_levels,refinement_ratios, &
      temperature_thresholds,tile_shape,tag_buffer,maximum_nodes
    namelist /amr_tree_distribution/ partitioning,rebalance_interval_steps,rebalance_minimum_gain
    namelist /amr_tree_timestepping/ subcycling
    sdc=trim(config%time_integrator)=='sdc'
    partitioning='round_robin'
    rebalance_interval_steps=0;rebalance_minimum_gain=0.1_dp
    subcycling=.false.
    regrid_interval_steps=0;refinement_levels=1;refinement_ratios=2
    temperature_thresholds=1200;tile_shape=4;tag_buffer=1;maximum_nodes=4096;last_tag_step=0
    checkpoint_indexed=.false.
    checkpoint_file='';restart_file='';checkpoint_interval_steps=0;stop_after_checkpoint=.false.
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) call fatal('Cannot read MPI rank')
    call MPI_Comm_size(comm,nranks,status)
    if (status/=MPI_SUCCESS) call fatal('Cannot read MPI size')
    call check_paths(input_path)
    call check_paths(topology_path)
    ! Read tree-specific controls before the shared regular-grid configuration.
    open(newunit=unit,file=input_path,status='old',action='read',iostat=status)
    if (status/=0) call fatal('Cannot open MPI tree input')
    do
      read(unit,'(a)',iostat=status) message
      if (status<0) exit
      if (status/=0) call fatal('Cannot scan MPI tree input')
      message=adjustl(message)
      do i=1,len_trim(message)
        c=iachar(message(i:i))
        if (c>=iachar('A').and.c<=iachar('Z')) message(i:i)=achar(c+32)
      end do
      if (index(message,'&amr_tree_distribution')==1.or.index(message,'$amr_tree_distribution')==1) then
        backspace(unit,iostat=status)
        if (status/=0) call fatal('Cannot locate tree distribution configuration')
        read(unit,nml=amr_tree_distribution,iostat=status)
        if (status/=0) call fatal('Invalid tree distribution configuration')
        cycle
      end if
      if (index(message,'&amr_tree_timestepping')==1.or.index(message,'$amr_tree_timestepping')==1) then
        backspace(unit,iostat=status)
        if (status/=0) call fatal('Cannot locate tree timestepping configuration')
        read(unit,nml=amr_tree_timestepping,iostat=status)
        if (status/=0) call fatal('Invalid tree timestepping configuration')
        cycle
      end if
      if (index(message,'&amr_tree_adaptation')==1.or.index(message,'$amr_tree_adaptation')==1) then
        backspace(unit,iostat=status)
        if (status/=0) call fatal('Cannot locate tree adaptation configuration')
        read(unit,nml=amr_tree_adaptation,iostat=status)
        if (status/=0) call fatal('Invalid tree adaptation configuration')
        cycle
      end if
      if (index(message,'&amr_tree_checkpoint')==1.or.index(message,'$amr_tree_checkpoint')==1) then
        backspace(unit,iostat=status)
        if (status/=0) call fatal('Cannot locate MPI tree checkpoint group')
        read(unit,nml=amr_tree_checkpoint,iostat=status)
        if (status/=0) call fatal('Invalid MPI tree checkpoint configuration')
      else if (index(message,'&amr_tree_')==1.or.index(message,'$amr_tree_')==1) then
        call fatal('Unsupported MPI tree configuration group')
      end if
    end do
    close(unit)
    if (partitioning/='round_robin'.and.partitioning/='weighted') call fatal('Invalid tree partitioning')
    if(.not.sdc.and.any(config%external_forcing/=0)) call fatal('MPI tree forcing requires SDC')
    if (rebalance_interval_steps<0) call fatal('Invalid rebalance interval')
    if (.not.ieee_is_finite(rebalance_minimum_gain)) call fatal('Invalid rebalance gain')
    if (rebalance_minimum_gain<0.or.rebalance_minimum_gain>=1) call fatal('Invalid rebalance gain')
    if (regrid_interval_steps<0) call fatal('Negative automatic tree regrid interval')
    if (regrid_interval_steps>0) then
      if (present(schedule_path)) call fatal('Automatic tree adaptation cannot use a prescribed schedule')
      if (refinement_levels<0.or.refinement_levels>maximum_tag_levels) call fatal('Invalid tree refinement depth')
      if (any(refinement_ratios(:refinement_levels)<2)) call fatal('Invalid tree refinement ratios')
      if (.not.all(ieee_is_finite(temperature_thresholds(:refinement_levels)))) call fatal('Nonfinite tree thresholds')
      if (any(temperature_thresholds(:refinement_levels)<=0)) call fatal('Invalid tree temperature thresholds')
      do i=2,refinement_levels
        if (temperature_thresholds(i)<temperature_thresholds(i-1)) call fatal('Tree thresholds must not decrease')
      end do
      if (any(tile_shape<1).or.tag_buffer<0.or.maximum_nodes<1) call fatal('Invalid tree tagging bounds')
    end if
    if (checkpoint_interval_steps<0) call fatal('Negative MPI tree checkpoint interval')
    if (stop_after_checkpoint.and.checkpoint_interval_steps==0) call fatal('Checkpoint stop requires positive interval')
    if (checkpoint_interval_steps>0.and.len_trim(checkpoint_file)==0) call fatal('Missing MPI tree checkpoint path')
    call check_checkpoint_paths(input_path)
    call check_checkpoint_paths(topology_path)
    call check_checkpoint_paths(config%output_file)
    call check_checkpoint_paths(trim(config%output_file)//'.budget.csv')
    if (present(schedule_path)) then
      call check_paths(schedule_path)
      call check_checkpoint_paths(schedule_path)
    end if
    call read_topology(topology_path,tree)
    schedule_identity=tree_key(tree)
    if (regrid_interval_steps>0) then
      write(message,'(*(i0,1x))') regrid_interval_steps,refinement_levels,tile_shape,tag_buffer,maximum_nodes
      schedule_identity=schedule_identity//'automatic-v1 '//trim(message)//new_line('a')
      do i=1,refinement_levels
        write(message,'(i0,1x,es27.18e3)') refinement_ratios(i),temperature_thresholds(i)
        schedule_identity=schedule_identity//trim(message)//new_line('a')
      end do
    end if
    nevents=0
    if (present(schedule_path)) then
      open(newunit=unit,file=schedule_path,status='old',action='read',iostat=status)
      if (status/=0) call fatal('Cannot open tree regrid schedule')
      read(unit,*,iostat=status) nevents
      if (status/=0) call fatal('Invalid tree regrid event count')
      if (nevents<0) call fatal('Invalid tree regrid event count')
    end if
    allocate(events(nevents),event_times(nevents),stat=status)
    if (status/=0) call fatal('Cannot allocate tree regrid schedule')
    do event=1,nevents
      read(unit,*,iostat=status) event_times(event),event_path
      if (status/=0) call fatal('Invalid tree regrid event')
      if (.not.ieee_is_finite(event_times(event))) call fatal('Nonfinite tree regrid time')
      if (event_times(event)<=0.or.event_times(event)>config%final_time) call fatal('Tree regrid time outside run')
      if (event>1) then
        if (event_times(event)<=event_times(event-1)) call fatal('Tree regrid times must increase')
      end if
      if (len_trim(event_path)==0) call fatal('Empty regrid topology path')
      if (event_path(1:1)/='/') then
        slash=scan(schedule_path,'/',back=.true.)
        event_path=schedule_path(:slash)//trim(event_path)
      end if
      call read_topology(trim(event_path),events(event))
      write(message,'(es27.18e3)') event_times(event)
      schedule_identity=schedule_identity//trim(message)//new_line('a')//tree_key(events(event))
    end do
    if (present(schedule_path)) close(unit)
    if (subcycling) schedule_identity=schedule_identity//'subcycling-v1'//new_line('a')

    if (checkpoint_interval_steps>0.or.len_trim(restart_file)>0) then
      identity='fixed:'//trim(config%thermo_model)
      if (present(mechanism_identity)) identity=mechanism_identity
      call resolve_reactive_chemistry_integrator(size(species),chemistry_integrator,use_implicit,ok)
      if (.not.ok) call fatal('Invalid MPI tree chemistry integrator')
      integrator='explicit'
      if (use_implicit) integrator='implicit'
      call build_amr_tree_checkpoint_context_3d(config,species,reactions,transport,mole_fractions, &
        identity,integrator,schedule_identity,context,ok)
      if (.not.ok) call fatal('Cannot build MPI tree checkpoint context')
    end if
    count=size(tree%nodes)-1
    call assign_owners(tree,owners)
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,ok)
    if (.not.ok) call fatal('Tree coordinate overflow')
    root_h=[(config%x_upper-config%x_lower)/config%nx,(config%y_upper-config%y_lower)/config%ny, &
      (config%z_upper-config%z_lower)/config%nz]
    allocate(h(3,count+1),state%fields(count+1),temperature%fields(count+1))
    allocate(mass_fractions(size(species)),inlet_x(size(species)),q(reactive_nvar(size(species))))
    allocate(initial(size(q)),final(size(q)),budget(size(q)),change(size(q)),inlet(size(q),6),inlet_cube(size(q),2,2,2))
    allocate(source_budget(size(q)),source_change(size(q)));source_budget=0;source_change=0
    if(sdc.and.config%chemistry_enabled) then
      allocate(reaction_source)
      allocate(reaction_source%fields(count+1))
      do node=1,count+1
        if(owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        allocate(reaction_source%fields(node)%values(size(q)-4,n(1),n(2),n(3)))
        reaction_source%fields(node)%values=0
      end do
    end if
    do node=1,count+1
      h(:,node)=root_h*real(domain(:,1),dp)/real(domain(:,tree%nodes(node)%level+1),dp)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      call coordinates(node)
      allocate(state%fields(node)%values(size(q),n(1),n(2),n(3)),temperature%fields(node)%values(1,n(1),n(2),n(3)))
      local=config;local%nx=n(1);local%ny=n(2);local%nz=n(3)
      call initialize_reactive_problem_3d(species,local,x,y,z,state%fields(node)%values, &
        temperature%fields(node)%values(1,:,:,:),base_density,mass_fractions,ok,mole_fractions)
      if (.not.ok) call fatal('Cannot initialize tree node')
    end do
    call average_down_mpi_amr_tree_3d(tree,owners,state,size(q),comm,ok)
    if (.not.ok) call fatal('Cannot synchronize initial MPI tree')
    call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,state,temperature,comm,ok)
    if (.not.ok) call fatal('Invalid initial MPI tree temperature')
    ! Empty owners still need the default reservoir composition.
    call mass_fractions_from_mole_fractions(species,mole_fractions,mass_fractions,ok)
    if (.not.ok) call fatal('Invalid initial mass fractions')
    faces=resolved_reactive_boundaries_3d(config)
    periodic=[faces(1)=="periodic",faces(3)=="periodic",faces(5)=="periodic"]
    viscosity=config%transport_enabled.and.config%viscosity_enabled
    conduction=config%transport_enabled.and.config%thermal_conduction_enabled
    diffusion=config%transport_enabled.and.config%species_diffusion_enabled
    barodiffusion=config%transport_enabled.and.config%barodiffusion_enabled
    min_h=minval(h,dim=2)
    if (regrid_interval_steps>0) then
      tag_spacing=root_h
      do i=1,refinement_levels
        tag_spacing=tag_spacing/real(refinement_ratios(i),dp)
      end do
      if (any(tag_spacing<=tiny(1.0_dp))) call fatal('Automatic tree spacing underflow')
      min_h=min(min_h,tag_spacing)
      ! Validate the requested integer lattice before any evolution/output.
      call tag_mpi_amr_tree_temperature_3d(tree,owners,temperature,temperature_thresholds(:refinement_levels), &
        refinement_ratios(:refinement_levels),tile_shape,tag_buffer,maximum_nodes,periodic,tagged_tree,comm,ok)
      if (.not.ok) call fatal('Invalid initial automatic tree proposal')
    end if
    do event=1,nevents
      call include_event_spacing(events(event))
    end do
    inlet=0;inlet_t=0;inlet_dt=huge(1.0_dp);mach_pressure=0
    where(config%inflow_mach>=0) mach_pressure=config%inflow_pressure
    if (any(faces=='inflow').and.config%inflow_composition_count/=0) then
      call resolve_selected_composition('Inflow',config%inflow_composition_count,config%inflow_composition_species, &
        config%inflow_composition_mole_fractions,species,inlet_x,ok,message)
      if (.not.ok) call fatal('Invalid inflow composition')
      call mass_fractions_from_mole_fractions(species,inlet_x,mass_fractions,ok)
      if (.not.ok) call fatal('Invalid inflow mass fractions')
    end if
    do face=1,6
      if (faces(face)/='inflow') cycle
      if(config%inflow_mach(face)>=0) then
        inlet_t(face)=config%inflow_temperature(face)
        cycle
      end if
      q=0;q(1)=mixture_density(species,mass_fractions,config%inflow_pressure(face),config%inflow_temperature(face),ok)
      if (.not.ok) call fatal('Invalid inflow density')
      q(2:4)=config%inflow_velocity(:,face);q(5)=config%inflow_pressure(face);q(6:)=mass_fractions
      call reactive_primitive_to_conserved(species,q,inlet(:,face),inlet_t(face),sound,ok)
      if (.not.ok) call fatal('Invalid inflow state')
      if (subcycling) min_h=root_h
      inlet_dt=min(inlet_dt,config%cfl/sum((abs(q(2:4))+sound)/min_h))
      if (config%transport_enabled) then
        do k=1,2
          do j=1,2
            do i=1,2
              inlet_cube(:,i,j,k)=inlet(:,face)
            end do
          end do
        end do
        inlet_temperature=inlet_t(face)
        call reactive_transport_timestep_3d(species,transport,inlet_cube,inlet_temperature,2,2,2, &
          min_h(1),min_h(2),min_h(3),config%transport_cfl,viscosity,conduction,diffusion,face_dt,diffusivity,ok)
        if (.not.ok) call fatal('Invalid inflow transport bound')
        inlet_dt=min(inlet_dt,face_dt)
      end if
    end do
    call distributed_inventory(tree,owners,state,initial,.false.)
    time=0;step=0;budget=0;max_correction=0;min_theta=1
    event=1;regrids=0;max_remap_error=0;last_checkpoint_step=0
    if (len_trim(restart_file)>0) call restore_checkpoint()
    if (rank==0) then
      write(*,'(a)') 'PeleF '//pelef_version//' experimental MPI regular 3D patch tree'
      write(*,'(a,l1)') 'Subcycling: ',subcycling
      write(*,'(a,a)') 'Partitioning: ',trim(partitioning)
      write(*,'(a,i0,a,i0,a,i0)') 'Ranks: ',nranks,' Nodes: ',count+1,' Levels: ',tree%level_count()
    end if
    do
      if (event<=nevents) then
        if (time==event_times(event)) then
          call apply_event(events(event))
          regrids=regrids+1;event=event+1
        end if
      end if
      if (regrid_interval_steps>0) then
        if (step>last_tag_step) then
          if (mod(step,regrid_interval_steps)==0) then
            call tag_mpi_amr_tree_temperature_3d(tree,owners,temperature,temperature_thresholds(:refinement_levels), &
              refinement_ratios(:refinement_levels),tile_shape,tag_buffer,maximum_nodes,periodic,tagged_tree,comm,ok)
            if (.not.ok) call fatal('Automatic tree tagging failed')
            if (tree_key(tree)/=tree_key(tagged_tree)) call apply_event(tagged_tree)
            regrids=regrids+1;last_tag_step=step
            if (rank==0) write(*,'(a,i0,a,i0)') 'Automatic regrid check at step ',step,' Nodes: ',size(tree%nodes)
          end if
        end if
      end if
      if (checkpoint_interval_steps>0) then
        if (step>last_checkpoint_step.and.mod(step,checkpoint_interval_steps)==0) then
          call save_checkpoint()
          last_checkpoint_step=step
          if (stop_after_checkpoint) then
            if (rank==0) write(*,'(a,i0,a,es24.16)') 'Stopped after checkpoint at step ',step,' time ',time
            return
          end if
        end if
      end if
      if (time>=config%final_time) exit
      stop_time=config%final_time
      if (event<=nevents) stop_time=event_times(event)
      if (step>=config%maximum_steps) call fatal('Tree maximum_steps reached')
      if (subcycling) then
        call compute_mpi_amr_tree_subcycle_cfl_3d(tree,owners,species,state,temperature,root_h,config%cfl,dt,comm,ok)
      else
        call compute_mpi_amr_tree_cfl_3d(tree,owners,species,state,temperature,root_h,config%cfl,dt,comm,ok)
      end if
      if (.not.ok) call fatal('Invalid MPI tree hydro timestep')
      dt=min(dt,inlet_dt,stop_time-time)
      if (config%transport_enabled) then
        if (subcycling) then
          call compute_mpi_amr_tree_subcycle_transport_timestep_3d(tree,owners,species,transport,state,temperature,root_h, &
            config%transport_cfl,viscosity,conduction,diffusion,transport_dt,diffusivity,comm,ok,faces,config%wall_temperature)
        else
          call compute_mpi_amr_tree_transport_timestep_3d(tree,owners,species,transport,state,temperature,root_h, &
            config%transport_cfl,viscosity,conduction,diffusion,transport_dt,diffusivity,comm,ok,faces,config%wall_temperature)
        end if
        if (.not.ok) call fatal('Invalid MPI tree transport timestep')
        dt=min(dt,transport_dt)
      end if
      if(any(config%inflow_mach>=0)) call limit_mach_step()
      sample_work=.false.
      if (rebalance_interval_steps>0) sample_work=mod(step+1,rebalance_interval_steps)==0
      do attempt=1,20
        if (time+dt<=time) call fatal('Tree timestep cannot advance time')
        if (sample_work) then
          call begin_tree_work(size(tree%nodes),ok)
          if (.not.ok) call fatal('Cannot allocate tree work counters')
        end if
        if(sdc) then
          if(subcycling) then
          call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,state,temperature,root_h,dt,faces, &
            config%sdc_iterations,config%chemistry_enabled,[viscosity,conduction,diffusion,barodiffusion], &
            config%external_forcing,config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance, &
            change,source_change,correction,theta,comm,ok,reaction_source=reaction_source, &
            inflow_state=inlet,inflow_temperature=inlet_t,wall_temperature=config%wall_temperature, &
            wall_velocity=config%wall_velocity,outlet_pressure=config%outlet_pressure, &
            inflow_mach=config%inflow_mach,inflow_pressure=mach_pressure,difmag=resolved_reactive_difmag_3d(config))
          else
          call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,state,temperature,root_h,dt,faces, &
            config%sdc_iterations,config%chemistry_enabled,[viscosity,conduction,diffusion,barodiffusion], &
            config%external_forcing,config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance, &
            change,source_change,correction,theta,comm,ok,reaction_source=reaction_source, &
            inflow_state=inlet,inflow_temperature=inlet_t,wall_temperature=config%wall_temperature, &
            wall_velocity=config%wall_velocity,outlet_pressure=config%outlet_pressure, &
            inflow_mach=config%inflow_mach,inflow_pressure=mach_pressure,difmag=resolved_reactive_difmag_3d(config))
          end if
        else
        call advance_mpi_amr_tree_split_3d(tree,owners,species,reactions,transport,state,temperature,root_h,dt,faces, &
          config%reconstruction,config%limiter,config%riemann_solver,config%chemistry_enabled, &
          viscosity,conduction,diffusion,barodiffusion,config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance, &
          change,correction,theta,comm,ok,chemistry_integrator=chemistry_integrator,inflow_state=inlet,inflow_temperature=inlet_t, &
          wall_temperature=config%wall_temperature,wall_velocity=config%wall_velocity,subcycling=subcycling, &
          outlet_pressure=config%outlet_pressure,inflow_mach=config%inflow_mach,inflow_pressure=mach_pressure, &
          difmag=resolved_reactive_difmag_3d(config))
        end if
        if (sample_work) call end_tree_work(measured_cost)
        if (ok) exit
        dt=.5_dp*dt
      end do
      if (.not.ok) call fatal('Tree step failed after 20 attempts')
      if (dt==stop_time-time) then
        time=stop_time
      else
        time=time+dt
      end if
      step=step+1;budget=budget+change;source_budget=source_budget+source_change
      if (sample_work) then
        call rebalance_fields()
        if (.not.ok) call fatal('Cannot rebalance MPI tree')
        if (rank==0) write(*,'(a,i0,a,l1)') 'Rebalance step: ',step,' migrated: ',redistributed
      end if
      max_correction=max(max_correction,correction);min_theta=min(min_theta,theta)
    end do
    call distributed_inventory(tree,owners,state,final,.false.)
    call write_mpi_amr_tree_csv_3d(trim(config%output_file),tree,owners,species,state,temperature, &
      [config%x_lower,config%y_lower,config%z_lower],root_h,time,comm,ok)
    if (.not.ok) call fatal('Cannot write distributed MPI tree CSV')
    if (rank==0) then
    open(newunit=unit,file=trim(config%output_file)//'.budget.csv',status='replace',action='write',iostat=status)
    if (status/=0) call fatal('Cannot open tree budget CSV')
    if(sdc) then
      write(unit,'(a)',iostat=status) 'component,initial,final,boundary_change,source_change'
    else
      write(unit,'(a)',iostat=status) 'component,initial,final,boundary_change'
    end if
    if (status/=0) call fatal('Cannot write tree budget header')
    do c=1,size(q)
      if(sdc) then
        write(unit,'(i0,4(",",es25.16e3))',iostat=status) c,initial(c),final(c),budget(c),source_budget(c)
      else
        write(unit,'(i0,3(",",es25.16e3))',iostat=status) c,initial(c),final(c),budget(c)
      end if
      if (status/=0) call fatal('Cannot write tree budget row')
    end do
    close(unit,iostat=status)
    if (status/=0) call fatal('Cannot close tree budget CSV')
    write(*,'(a,i0,a,es12.4)') 'Regrids: ',regrids,' Maximum remap error: ',max_remap_error
    write(*,'(a,i0,a,es24.16)') 'Steps: ',step,' Final time: ',time
    write(*,'(a,es12.4,a,es12.4)') 'Maximum reflux: ',max_correction,' Minimum transport theta: ',min_theta
    end if
    call MPI_Barrier(comm,status)
    if (status/=MPI_SUCCESS) call fatal('Cannot complete MPI tree output')
  contains
    subroutine limit_mach_step()
      type(amr_tree_state_3d) :: packed,halos
      real(dp) :: bound,dmax,factor,local_bound,global_bound
      integer :: peer,dims(3),nv,width
      character(len=:), allocatable :: donor_limiter
      nv=size(q)
      if(config%reconstruction/='pcm') donor_limiter=config%limiter
      width=2
      if(config%reconstruction=='characteristic_ppm') width=4
      allocate(packed%fields(size(tree%nodes)))
      do peer=1,size(tree%nodes)
        if(owners(peer)/=rank) cycle
        dims=tree%nodes(peer)%shape
        allocate(packed%fields(peer)%values(nv+1,dims(1),dims(2),dims(3)))
        packed%fields(peer)%values(:nv,:,:,:)=state%fields(peer)%values
        packed%fields(peer)%values(nv+1,:,:,:)=temperature%fields(peer)%values(1,:,:,:)
      end do
      call fill_mpi_amr_tree_ghosts_3d(tree,owners,packed,nv+1,faces,halos,comm,ok,viscous_wall=.true., &
        inflow_state=inlet_with_temperature(),wall_velocity=config%wall_velocity,coarse_limiter=donor_limiter, &
        species_count=size(species),ghost_width=width,species=species,inflow_mach=config%inflow_mach, &
        inflow_temperature=inlet_t,inflow_pressure=mach_pressure)
      if(.not.ok) call fatal('Invalid MPI tree Mach inlet donors')
      local_bound=dt
      do peer=1,size(tree%nodes)
        if(owners(peer)/=rank) cycle
        associate(ghost=>halos%fields(peer)%values)
        dims=shape(ghost(1,:,:,:));factor=1
        if(subcycling) factor=root_h(1)/h(1,peer)
        call compute_reactive_cfl_timestep_3d(species,ghost(:nv,:,:,:),ghost(nv+1,:,:,:), &
          dims(1),dims(2),dims(3),h(1,peer),h(2,peer),h(3,peer),config%cfl,bound,ok)
        if(.not.ok) call fatal('Invalid MPI tree Mach hydro timestep')
        local_bound=min(local_bound,factor*bound)
        if(config%transport_enabled) then
          call reactive_transport_timestep_3d(species,transport,ghost(:nv,:,:,:),ghost(nv+1,:,:,:), &
            dims(1),dims(2),dims(3),h(1,peer),h(2,peer),h(3,peer),config%transport_cfl, &
            viscosity,conduction,diffusion,bound,dmax,ok)
          if(.not.ok) call fatal('Invalid MPI tree Mach transport timestep')
          local_bound=min(local_bound,factor**2*bound)
        end if
        end associate
      end do
      call MPI_Allreduce(local_bound,global_bound,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) call fatal('Cannot reduce MPI tree Mach timestep')
      dt=global_bound
    end subroutine
    function inlet_with_temperature() result(packed_inlet)
      real(dp) :: packed_inlet(size(q)+1,6)
      packed_inlet(:size(q),:)=inlet;packed_inlet(size(q)+1,:)=inlet_t
    end function
    subroutine fatal(reason)
      character(len=*), intent(in) :: reason
      integer :: ierr
      write(*,'(a,i0,2a)') 'MPI tree rank ',rank,': ',reason
      call MPI_Abort(comm,2,ierr)
      error stop 2
    end subroutine
    subroutine check_checkpoint_paths(path)
      character(len=*), intent(in) :: path
      if (checkpoint_interval_steps>0) then
        if (selected_reactor_paths_alias_lexically(absolute_path(checkpoint_file),absolute_path(path))) &
          call fatal('MPI tree checkpoint aliases input/output')
      end if
      if (len_trim(restart_file)>0) then
        if (selected_reactor_paths_alias_lexically(absolute_path(restart_file),absolute_path(path))) &
          call fatal('MPI tree restart aliases input/output')
      end if
    end subroutine
    function tree_key(geometry) result(key)
      type(amr_patch_tree_3d), intent(in) :: geometry
      character(len=:), allocatable :: key
      character(len=512) :: line
      integer :: peer
      write(line,'(*(i0,1x))') size(geometry%nodes),geometry%nodes(1)%shape
      key=trim(line)//new_line('a')
      do peer=2,size(geometry%nodes)
        associate(item=>geometry%nodes(peer),box=>geometry%nodes(peer)%patch)
          write(line,'(*(i0,1x))') item%parent,box%coarse_i_lower,box%coarse_i_upper, &
            box%coarse_j_lower,box%coarse_j_upper,box%coarse_k_lower,box%coarse_k_upper,box%refinement_ratio
        end associate
        key=key//trim(line)//new_line('a')
      end do
    end function
    subroutine rebalance_fields()
      type(amr_tree_state_3d) :: packed,temps,moved_state,moved_source
      integer, allocatable :: assigned(:)
      integer :: peer,nv,dims(3)
      if(.not.allocated(reaction_source)) then
        call rebalance_mpi_amr_tree_3d(tree,owners,size(q),state,temperature, &
          measured_cost,rebalance_minimum_gain,comm,redistributed,ok)
        return
      end if
      nv=size(q);assigned=owners;temps=temperature
      allocate(packed%fields(size(tree%nodes)))
      do peer=1,size(tree%nodes)
        if(owners(peer)/=rank) cycle
        dims=tree%nodes(peer)%shape
        allocate(packed%fields(peer)%values(2*nv-4,dims(1),dims(2),dims(3)))
        packed%fields(peer)%values(:nv,:,:,:)=state%fields(peer)%values
        packed%fields(peer)%values(nv+1:,:,:,:)=reaction_source%fields(peer)%values
      end do
      call rebalance_mpi_amr_tree_3d(tree,assigned,2*nv-4,packed,temps, &
        measured_cost,rebalance_minimum_gain,comm,redistributed,ok)
      if(.not.ok.or..not.redistributed) return
      allocate(moved_state%fields(size(tree%nodes)),moved_source%fields(size(tree%nodes)))
      do peer=1,size(tree%nodes)
        if(assigned(peer)/=rank) cycle
        moved_state%fields(peer)%values=packed%fields(peer)%values(:nv,:,:,:)
        moved_source%fields(peer)%values=packed%fields(peer)%values(nv+1:,:,:,:)
      end do
      owners=assigned
      call move_alloc(moved_state%fields,state%fields)
      call move_alloc(moved_source%fields,reaction_source%fields)
      call move_alloc(temps%fields,temperature%fields)
    end subroutine rebalance_fields
    subroutine save_checkpoint()
      type(amr_tree_checkpoint_3d) :: snapshot
      snapshot%tree=tree;snapshot%initial=initial;snapshot%boundary=budget
      if(sdc) snapshot%source=source_budget
      if(allocated(reaction_source)) snapshot%reaction_source=reaction_source
      snapshot%time=time;snapshot%steps=step;snapshot%regrids=regrids
      snapshot%maximum_remap_error=max_remap_error
      snapshot%maximum_reflux=max_correction;snapshot%minimum_theta=min_theta
      call write_mpi_amr_tree_checkpoint_3d(trim(checkpoint_file),species,context,snapshot,owners,state,temperature,comm,ok, &
        indexed=checkpoint_indexed)
      if (.not.ok) call fatal('Cannot write MPI tree checkpoint')
    end subroutine
    subroutine restore_checkpoint()
      type(amr_tree_checkpoint_3d) :: snapshot
      integer, allocatable :: restored_owners(:)
      real(dp) :: original_scale(size(q)),restart_scale(size(q)),inventory(size(q))
      integer :: nconserved,applied,exponent
      ! The initialized owner fields supply the cancellation-safe ledger scale.
      call distributed_inventory(tree,owners,state,original_scale,.true.)
      exponent=0
      if (subcycling) then
        exponent=1
        if (config%transport_enabled) exponent=2
      end if
      call read_mpi_amr_tree_checkpoint_3d(trim(restart_file),species,context, &
        [config%nx,config%ny,config%nz],partitioning,exponent,snapshot,restored_owners,comm,ok)
      if (.not.ok) call fatal('Invalid MPI tree restart checkpoint/context')
      if(sdc.neqv.allocated(snapshot%source)) call fatal('MPI tree source ledger presence mismatch')
      if(allocated(snapshot%reaction_source%fields).neqv.(sdc.and.config%chemistry_enabled)) &
        call fatal('MPI tree reaction estimate presence mismatch')
      if(sdc) source_budget=snapshot%source
      if (snapshot%time>config%final_time.or.snapshot%steps>config%maximum_steps) &
        call fatal('MPI tree restart exceeds run limits')
      applied=count_events(snapshot%time)
      if (regrid_interval_steps>0) applied=snapshot%steps/regrid_interval_steps
      if (snapshot%regrids/=applied) call fatal('MPI tree restart adaptation/schedule position mismatch')
      if (regrid_interval_steps>0.and.applied>0) then
        call validate_adaptive_geometry(snapshot%tree)
      else
        expected_tree=tree_key(tree)
        if (applied>0) expected_tree=tree_key(events(applied))
        if (tree_key(snapshot%tree)/=expected_tree) call fatal('MPI tree restart geometry mismatch')
      end if
      if (any(abs(snapshot%initial-initial)>64*epsilon(1.0_dp)*max(tiny(1.0_dp),abs(initial)))) &
        call fatal('Tree restart initial inventory mismatch')
      ! Use absolute cell inventories so cancelling momentum has a meaningful scale.
      call distributed_inventory(snapshot%tree,restored_owners,snapshot%state,restart_scale,.true.)
      call distributed_inventory(snapshot%tree,restored_owners,snapshot%state,inventory,.false.)
      ! Momentum is a vector: a nearly zero transverse component must use
      ! the full momentum magnitude as its roundoff scale.
      original_scale(2:4)=maxval(original_scale(2:4))
      restart_scale(2:4)=maxval(restart_scale(2:4))
      nconserved=size(q)
      if (config%chemistry_enabled) nconserved=5
      if (any(abs(inventory(:nconserved)-snapshot%initial(:nconserved)- &
          snapshot%boundary(:nconserved)-source_budget(:nconserved))> &
          2e-10_dp*max(tiny(1.0_dp),original_scale(:nconserved),restart_scale(:nconserved), &
          abs(snapshot%boundary(:nconserved)),abs(source_budget(:nconserved))))) call fatal('Tree restart boundary ledger mismatch')
      initial=snapshot%initial;budget=snapshot%boundary;step=snapshot%steps
      time=snapshot%time;max_correction=snapshot%maximum_reflux
      min_theta=snapshot%minimum_theta;max_remap_error=snapshot%maximum_remap_error
      regrids=count_events(time);event=regrids+1;last_checkpoint_step=step;last_tag_step=step
      if (regrid_interval_steps>0) regrids=step/regrid_interval_steps
      tree=snapshot%tree
      call reset_layout()
      if (any(owners/=restored_owners)) call fatal('Restart owner assignment mismatch')
      ! All context, EOS, average and ledger checks completed before publication.
      call move_alloc(snapshot%state%fields,state%fields)
      call move_alloc(snapshot%temperature%fields,temperature%fields)
      if(allocated(reaction_source)) call move_alloc(snapshot%reaction_source%fields,reaction_source%fields)
    end subroutine
    subroutine check_paths(path)
      character(len=*), intent(in) :: path
      if (selected_reactor_paths_alias_lexically(absolute_path(path),absolute_path(config%output_file)).or. &
          selected_reactor_paths_alias_lexically(absolute_path(path),absolute_path(trim(config%output_file)//'.budget.csv'))) &
        call fatal('MPI tree output aliases an input')
    end subroutine
    function absolute_path(path) result(value)
      character(len=*), intent(in) :: path
      character(len=:), allocatable :: value
      character(kind=c_char) :: buffer(4096)
      type(c_ptr) :: pointer
      integer :: index
      if (len_trim(path)==0) call fatal('Empty MPI tree path')
      if (path(1:1)=='/') then
        value=trim(path)
      else
        pointer=system_getcwd(buffer,int(size(buffer),c_size_t))
        if (.not.c_associated(pointer)) call fatal('Cannot resolve MPI tree working directory')
        value=''
        do index=1,size(buffer)
          if (buffer(index)==c_null_char) exit
          value=value//buffer(index)
        end do
        value=value//'/'//trim(path)
      end if
      if (len(value)>1024) call fatal('MPI tree path exceeds lexical comparison capacity')
    end function
    subroutine read_topology(path,result)
      character(len=*), intent(in) :: path
      type(amr_patch_tree_3d), intent(out) :: result
      type(amr_patch_3d), allocatable :: patches(:)
      integer, allocatable :: parents(:),shapes(:,:)
      integer :: input,ios,nchild,peer,parent,dims(3),lower(3),upper(3),ratio
      logical :: valid
      call check_paths(path)
      call check_checkpoint_paths(path)
      open(newunit=input,file=path,status='old',action='read',iostat=ios)
      if (ios/=0) call fatal('Cannot open tree topology')
      read(input,*,iostat=ios) nchild
      if (ios/=0) call fatal('Invalid tree child count')
      if (nchild<0.or.nchild==huge(nchild)) call fatal('Invalid tree child count')
      allocate(patches(nchild),parents(nchild),shapes(3,nchild+1),stat=ios)
      if (ios/=0) call fatal('Cannot allocate tree topology')
      shapes(:,1)=[config%nx,config%ny,config%nz]
      do peer=1,nchild
        read(input,*,iostat=ios) parent,lower(1),upper(1),lower(2),upper(2),lower(3),upper(3),ratio
        if (ios/=0) call fatal('Invalid tree patch record')
        if (parent<1.or.parent>peer) call fatal('Tree parents must precede children')
        parents(peer)=parent;dims=shapes(:,parent)
        call initialize_amr_patch_3d(dims(1),dims(2),dims(3),lower(1),upper(1),lower(2),upper(2), &
          lower(3),upper(3),ratio,patches(peer),valid)
        if (.not.valid) call fatal('Invalid tree patch geometry')
        if (any(int(upper-lower+1,int64)*int(ratio,int64)>int(huge(1),int64))) call fatal('Tree patch extent overflow')
        shapes(:,peer+1)=[patches(peer)%fine_nx(),patches(peer)%fine_ny(),patches(peer)%fine_nz()]
      end do
      close(input)
      call initialize_amr_patch_tree_3d(shapes(:,1),parents,patches,result,valid)
      if (.not.valid) call fatal('Invalid tree topology')
    end subroutine
    subroutine validate_adaptive_geometry(geometry)
      type(amr_patch_tree_3d), intent(in) :: geometry
      integer :: peer,depth,parent,lower(3),upper(3),dims(3)
      if (size(geometry%nodes)>maximum_nodes.or.geometry%level_count()>refinement_levels+1) &
        call fatal('Tree restart exceeds adaptation bounds')
      do peer=2,size(geometry%nodes)
        depth=geometry%nodes(peer)%level;parent=geometry%nodes(peer)%parent
        dims=geometry%nodes(parent)%shape
        associate(box=>geometry%nodes(peer)%patch)
          lower=[box%coarse_i_lower,box%coarse_j_lower,box%coarse_k_lower]
          upper=[box%coarse_i_upper,box%coarse_j_upper,box%coarse_k_upper]
          if (box%refinement_ratio/=refinement_ratios(depth)) call fatal('Tree restart adaptation ratio mismatch')
        end associate
        if (any(mod(lower-1,tile_shape)/=0)) call fatal('Tree restart tile alignment mismatch')
        if (any(int(upper,int64)/=min(int(dims,int64),int(lower,int64)+int(tile_shape,int64)-1))) &
          call fatal('Tree restart tile extent mismatch')
      end do
    end subroutine
    integer function count_events(at_time) result(applied)
      real(dp), intent(in) :: at_time
      integer :: peer
      applied=0
      do peer=1,nevents
        if (event_times(peer)<=at_time) applied=applied+1
      end do
    end function
    subroutine include_event_spacing(next_tree)
      type(amr_patch_tree_3d), intent(in) :: next_tree
      integer(int64), allocatable :: o(:,:),d(:,:)
      integer, allocatable :: r(:)
      logical :: valid
      call amr_tree_integer_layout_3d(next_tree,o,d,r,valid)
      if (.not.valid) call fatal('Scheduled tree coordinate overflow')
      min_h=min(min_h,root_h*real(d(:,1),dp)/real(d(:,size(r)),dp))
    end subroutine
    subroutine assign_owners(geometry,assigned)
      type(amr_patch_tree_3d), intent(in) :: geometry
      integer, allocatable, intent(inout) :: assigned(:)
      integer :: peer,exponent
      logical :: valid
      if (partitioning=='round_robin') then
        assigned=[(mod(peer-1,nranks),peer=1,size(geometry%nodes))]
        return
      end if
      exponent=0
      if (subcycling) then
        exponent=1
        if (config%transport_enabled) exponent=2
      end if
      call partition_amr_tree_3d(geometry,nranks,exponent,assigned,valid)
      if (.not.valid) call fatal('Cannot partition MPI tree')
    end subroutine
    subroutine reset_layout()
      integer :: peer
      call amr_tree_integer_layout_3d(tree,origin,domain,ratios,ok)
      if (.not.ok) call fatal('Regridded tree coordinate overflow')
      count=size(tree%nodes)-1
      call assign_owners(tree,owners)
      if (allocated(h)) deallocate(h)
      allocate(h(3,count+1))
      do peer=1,count+1
        h(:,peer)=root_h*real(domain(:,1),dp)/real(domain(:,tree%nodes(peer)%level+1),dp)
      end do
    end subroutine
    subroutine distributed_inventory(geometry,assigned,fields,result,absolute)
      type(amr_patch_tree_3d), intent(in) :: geometry
      integer, intent(in) :: assigned(:)
      type(amr_tree_state_3d), intent(in) :: fields
      real(dp), intent(out) :: result(:)
      logical, intent(in) :: absolute
      real(dp) :: accum(size(result),2),term(size(result)),updated(size(result)),value(size(result))
      real(dp) :: volume(size(geometry%nodes)),ratio
      integer :: peer,parent,ii,jj,kk,dims(3),ierr
      accum=0;volume(1)=product(root_h)
      if (.not.ieee_is_finite(volume(1)).or.volume(1)<=0) call fatal('Invalid MPI inventory volume')
      do peer=1,size(geometry%nodes)
        if (peer>1) then
          parent=geometry%nodes(peer)%parent;ratio=real(geometry%nodes(peer)%patch%refinement_ratio,dp)
          volume(peer)=volume(parent)/ratio/ratio/ratio
          if (volume(peer)<=0) call fatal('MPI inventory volume underflow')
        end if
        if (assigned(peer)==rank) then
          dims=geometry%nodes(peer)%shape
          do kk=1,dims(3)
            do jj=1,dims(2)
              do ii=1,dims(1)
                if (cell_covered(geometry,peer,[ii,jj,kk])) cycle
                value=fields%fields(peer)%values(:,ii,jj,kk)
                if (absolute) value=abs(value)
                term=volume(peer)*value-accum(:,2)
                updated=accum(:,1)+term
                accum(:,2)=(updated-accum(:,1))-term;accum(:,1)=updated
              end do
            end do
          end do
        end if
        ! Pass only the compensated accumulator, preserving the serial cell order.
        call MPI_Bcast(accum,size(accum),MPI_DOUBLE_PRECISION,assigned(peer),comm,ierr)
        if (ierr/=MPI_SUCCESS) call fatal('Cannot transfer MPI inventory accumulator')
      end do
      if (.not.all(ieee_is_finite(accum))) call fatal('Invalid MPI inventory')
      result=accum(:,1)
    end subroutine
    subroutine apply_event(next_tree)
      type(amr_patch_tree_3d), intent(in) :: next_tree
      type(amr_tree_state_3d) :: candidate,temps,remapped_reaction
      integer, allocatable :: assigned(:)
      real(dp) :: before(size(q)),after(size(q)),scale(size(q)),error
      call assign_owners(next_tree,assigned)
      call remap_mpi_amr_tree_reactive_3d(tree,owners,state,temperature,next_tree,assigned, &
        species,candidate,temps,comm,ok)
      if (.not.ok) call fatal('Scheduled MPI tree remap failed')
      if(allocated(reaction_source)) then
        call remap_mpi_amr_tree_3d(tree,owners,reaction_source,next_tree,assigned,size(q)-4,remapped_reaction,comm,ok)
        if(.not.ok) call fatal('MPI SDC reaction-source remap failed')
      end if
      call distributed_inventory(tree,owners,state,before,.false.)
      call distributed_inventory(next_tree,assigned,candidate,after,.false.)
      call distributed_inventory(tree,owners,state,scale,.true.)
      error=maxval(abs(after-before)/max(tiny(1.0_dp),scale))
      if (error>5e-12_dp) call fatal('Scheduled MPI tree remap lost conservation')
      tree=next_tree
      call move_alloc(candidate%fields,state%fields)
      call move_alloc(temps%fields,temperature%fields)
      if(allocated(reaction_source)) call move_alloc(remapped_reaction%fields,reaction_source%fields)
      call reset_layout()
      max_remap_error=max(max_remap_error,error)
      if (rank==0) write(*,'(a,es24.16,a,i0,a,i0)') &
        'Regrid time: ',time,' Nodes: ',count+1,' Levels: ',tree%level_count()
    end subroutine
    subroutine coordinates(peer)
      integer, intent(in) :: peer
      integer :: a,dims(3)
      dims=tree%nodes(peer)%shape
      if (allocated(x)) deallocate(x,y,z)
      allocate(x(dims(1)),y(dims(2)),z(dims(3)))
      x=[(config%x_lower+(real(origin(1,peer),dp)+real(a,dp)-.5_dp)*h(1,peer),a=1,dims(1))]
      y=[(config%y_lower+(real(origin(2,peer),dp)+real(a,dp)-.5_dp)*h(2,peer),a=1,dims(2))]
      z=[(config%z_lower+(real(origin(3,peer),dp)+real(a,dp)-.5_dp)*h(3,peer),a=1,dims(3))]
    end subroutine
  end subroutine
end module mpi_amr_patch_tree_application_3d_mod
