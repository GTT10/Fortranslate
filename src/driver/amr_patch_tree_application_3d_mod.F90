module amr_patch_tree_application_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use amr_patch_tree_tagging_3d_mod, only: tag_amr_patch_tree_temperature_3d
  use amr_patch_tree_remap_3d_mod, only: remap_amr_patch_tree_reactive_3d
  use, intrinsic :: iso_fortran_env, only: int64
  use amr_patch_tree_checkpoint_3d_mod, only: amr_tree_checkpoint_3d, build_amr_tree_checkpoint_context_3d, &
    read_amr_tree_checkpoint_3d,write_amr_tree_checkpoint_3d
  use simulation_config_selected_reactor_mod, only: selected_reactor_paths_alias_lexically
  use precision_mod, only: dp
  use amr_patch_tree_subcycle_3d_mod, only: compute_amr_patch_tree_subcycle_cfl_3d
  use amr_patch_tree_transport_subcycle_3d_mod, only: compute_amr_patch_tree_subcycle_transport_timestep_3d
  use constants_mod, only: pelef_version
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use mixture_thermo_mod, only: mixture_density,mass_fractions_from_mole_fractions
  use selected_composition_mod, only: resolve_selected_composition
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved, &
    resolve_reactive_chemistry_integrator
  use reactive_3d_mod, only: recover_reactive_temperatures_3d,compute_reactive_cfl_timestep_3d
  use amr_patch_tree_ghost_3d_mod, only: fill_amr_patch_tree_ghosts_3d
  use reactive_transport_3d_mod, only: reactive_transport_timestep_3d
  use simulation_config_reactive_3d_mod, only: reactive_3d_config,resolved_reactive_boundaries_3d
  use reactive_entropy_wave_3d_problem_mod, only: initialize_reactive_problem_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    amr_tree_integer_layout_3d,average_down_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
  use amr_patch_tree_hydro_3d_mod, only: compute_amr_patch_tree_cfl_3d
  use amr_patch_tree_transport_3d_mod, only: compute_amr_patch_tree_transport_timestep_3d
  use amr_patch_tree_split_3d_mod, only: advance_amr_patch_tree_split_3d,advance_amr_patch_tree_sdc_3d
  use reactive_csv_io_3d_mod, only: reactive_3d_csv_header,format_reactive_3d_csv_row
  use simulation_config_reactive_3d_mod, only: resolved_reactive_difmag_3d
  implicit none
  private
  public :: run_amr_patch_tree_application_3d
contains
  subroutine run_amr_patch_tree_application_3d(topology_path,config,species,reactions,transport,mole_fractions, &
      schedule_path,chemistry_integrator,input_path,mechanism_identity)
    character(len=*), intent(in) :: topology_path
    character(len=*), optional, intent(in) :: schedule_path,chemistry_integrator,input_path,mechanism_identity
    type(reactive_3d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: mole_fractions(:)
    type(reactive_3d_config) :: local
    type(amr_patch_tree_3d) :: tree,tagged_tree
    type(amr_tree_checkpoint_3d) :: snapshot
    type(amr_patch_tree_3d), allocatable :: events(:)
    type(amr_tree_state_3d) :: state,temperature
    type(amr_tree_state_3d), allocatable :: reaction_source
    integer, allocatable :: ratios(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    real(dp), allocatable :: x(:),y(:),z(:),root_h(:),h(:,:),recovered(:,:,:),mass_fractions(:),inlet_x(:)
    real(dp), allocatable :: initial(:),final(:),budget(:),change(:),inlet(:,:),q(:),inlet_cube(:,:,:,:)
    real(dp), allocatable :: event_times(:),source_budget(:),source_change(:)
    real(dp) :: stop_time,max_remap_error
    real(dp) :: time,dt,transport_dt,diffusivity,base_density,correction,theta,max_correction,min_theta
    real(dp) :: inlet_t(6),inlet_temperature(2,2,2),inlet_dt,face_dt,sound,min_h(3),mach_pressure(6)
    integer :: count,unit,status,node,n(3),lo(3),hi(3),i,j,k,step,attempt,face,c,child,nevents,event,slash,regrids
    character(len=32) :: faces(6)
    character(len=1024) :: message,event_path
    character(len=:), allocatable :: row
    logical :: ok,covered,viscosity,conduction,diffusion,barodiffusion,use_implicit
    character(len=1024) :: checkpoint_file,restart_file
    character(len=:), allocatable :: context,schedule_identity,identity,integrator,expected_tree
    integer :: checkpoint_interval_steps,last_checkpoint_step
    integer, parameter :: maximum_tag_levels=16
    integer :: regrid_interval_steps,refinement_levels,refinement_ratios(maximum_tag_levels)
    integer :: tile_shape(3),tag_buffer,maximum_nodes,last_tag_step
    real(dp) :: temperature_thresholds(maximum_tag_levels),tag_spacing(3)
    logical :: periodic(3),stop_after_checkpoint,subcycling
    logical :: checkpoint_indexed,sdc
    namelist /amr_tree_checkpoint/ checkpoint_file,restart_file,checkpoint_interval_steps,stop_after_checkpoint,checkpoint_indexed
    namelist /amr_tree_adaptation/ regrid_interval_steps,refinement_levels,refinement_ratios, &
      temperature_thresholds,tile_shape,tag_buffer,maximum_nodes
    namelist /amr_tree_timestepping/ subcycling
    sdc=trim(config%time_integrator)=='sdc'
    subcycling=.false.
    regrid_interval_steps=0;refinement_levels=1;refinement_ratios=2
    temperature_thresholds=1200;tile_shape=4;tag_buffer=1;maximum_nodes=4096;last_tag_step=0
    checkpoint_indexed=.false.
    checkpoint_file='';restart_file='';checkpoint_interval_steps=0;stop_after_checkpoint=.false.
    last_checkpoint_step=0
    if (present(input_path)) then
      open(newunit=unit,file=input_path,status='old',action='read',iostat=status)
      if (status/=0) error stop 'Cannot open tree checkpoint configuration'
      do
        read(unit,'(a)',iostat=status) message
        if (status<0) exit
        if (status/=0) error stop 'Cannot read tree checkpoint configuration'
        message=adjustl(message)
        do i=1,len_trim(message)
          c=iachar(message(i:i))
          if (c>=iachar('A').and.c<=iachar('Z')) message(i:i)=achar(c+32)
        end do
        if (index(message,'&amr_tree_timestepping')==1.or.index(message,'$amr_tree_timestepping')==1) then
          backspace(unit,iostat=status)
          if (status/=0) error stop 'Cannot locate tree timestepping configuration'
          read(unit,nml=amr_tree_timestepping,iostat=status)
          if (status/=0) error stop 'Invalid tree timestepping configuration'
          cycle
        end if
        if (index(message,'&amr_tree_adaptation')==1.or.index(message,'$amr_tree_adaptation')==1) then
          backspace(unit,iostat=status)
          if (status/=0) error stop 'Cannot locate tree adaptation configuration'
          read(unit,nml=amr_tree_adaptation,iostat=status)
          if (status/=0) error stop 'Invalid tree adaptation configuration'
          cycle
        end if
        if (index(message,'&amr_tree_checkpoint')/=1.and.index(message,'$amr_tree_checkpoint')/=1) cycle
        backspace(unit,iostat=status)
        if (status/=0) error stop 'Cannot locate tree checkpoint configuration'
        read(unit,nml=amr_tree_checkpoint,iostat=status)
        if (status/=0) error stop 'Invalid tree checkpoint configuration'
      end do
      close(unit)
    end if
    if (sdc) then
      if (subcycling) error stop 'SDC tree subcycling is not implemented'
      if (present(schedule_path).or.regrid_interval_steps>0) error stop 'SDC tree regridding is not implemented'
    else
      if (any(config%external_forcing/=0)) error stop 'Tree external forcing requires SDC'
    end if
    if (regrid_interval_steps<0) error stop 'Negative automatic tree regrid interval'
    if (regrid_interval_steps>0) then
      if (present(schedule_path)) error stop 'Automatic tree adaptation cannot use a prescribed schedule'
      if (refinement_levels<0.or.refinement_levels>maximum_tag_levels) error stop 'Invalid tree refinement depth'
      if (any(refinement_ratios(:refinement_levels)<2)) error stop 'Invalid tree refinement ratios'
      if (.not.all(ieee_is_finite(temperature_thresholds(:refinement_levels)))) error stop 'Nonfinite tree thresholds'
      if (any(temperature_thresholds(:refinement_levels)<=0)) error stop 'Invalid tree temperature thresholds'
      do i=2,refinement_levels
        if (temperature_thresholds(i)<temperature_thresholds(i-1)) error stop 'Tree thresholds must not decrease'
      end do
      if (any(tile_shape<1).or.tag_buffer<0.or.maximum_nodes<1) error stop 'Invalid tree tagging bounds'
    end if
    if (checkpoint_interval_steps<0) error stop 'Negative tree checkpoint interval'
    if (stop_after_checkpoint.and.checkpoint_interval_steps==0) error stop 'Checkpoint stop requires positive interval'
    if (checkpoint_interval_steps>0.and.len_trim(checkpoint_file)==0) error stop 'Missing tree checkpoint path'
    if (present(input_path)) call check_paths(input_path)
    if (present(schedule_path)) call check_paths(schedule_path)
    call check_paths(config%output_file)
    call check_paths(trim(config%output_file)//'.budget.csv')
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
    count=size(tree%nodes)-1
    nevents=0
    if (present(schedule_path)) then
      open(newunit=unit,file=schedule_path,status='old',action='read',iostat=status)
      if (status/=0) error stop 'Cannot open tree regrid schedule'
      read(unit,*,iostat=status) nevents
      if (status/=0) error stop 'Invalid tree regrid event count'
      if (nevents<0) error stop 'Invalid tree regrid event count'
    end if
    allocate(events(nevents),event_times(nevents),stat=status)
    if (status/=0) error stop 'Cannot allocate tree regrid schedule'
    do event=1,nevents
      read(unit,*,iostat=status) event_times(event),event_path
      if (status/=0) error stop 'Invalid tree regrid event'
      if (.not.ieee_is_finite(event_times(event))) error stop 'Nonfinite tree regrid time'
      if (event_times(event)<=0.or.event_times(event)>config%final_time) error stop 'Tree regrid time outside run'
      if (event>1) then
        if (event_times(event)<=event_times(event-1)) error stop 'Tree regrid times must increase'
      end if
      if (len_trim(event_path)==0) error stop 'Empty regrid topology path'
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
      if (.not.ok) error stop 'Invalid tree chemistry integrator'
      integrator='explicit'
      if (use_implicit) integrator='implicit'
      call build_amr_tree_checkpoint_context_3d(config,species,reactions,transport,mole_fractions, &
        identity,integrator,schedule_identity,context,ok)
      if (.not.ok) error stop 'Cannot build tree checkpoint context'
    end if
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,ok)
    if (.not.ok) error stop 'Tree coordinate overflow'
    root_h=[(config%x_upper-config%x_lower)/config%nx,(config%y_upper-config%y_lower)/config%ny, &
      (config%z_upper-config%z_lower)/config%nz]
    allocate(h(3,count+1),state%fields(count+1),temperature%fields(count+1))
    allocate(mass_fractions(size(species)),inlet_x(size(species)),q(reactive_nvar(size(species))))
    allocate(initial(size(q)),final(size(q)),budget(size(q)),change(size(q)),inlet(size(q),6),inlet_cube(size(q),2,2,2))
    allocate(source_budget(size(q)),source_change(size(q)))
    source_budget=0;source_change=0
    if (sdc.and.config%chemistry_enabled) then
      allocate(reaction_source)
      allocate(reaction_source%fields(count+1))
      do node=1,count+1
        n=tree%nodes(node)%shape
        allocate(reaction_source%fields(node)%values(size(q)-4,n(1),n(2),n(3)))
        reaction_source%fields(node)%values=0
      end do
    end if
    do node=1,count+1
      h(:,node)=root_h*real(domain(:,1),dp)/real(domain(:,tree%nodes(node)%level+1),dp)
      n=tree%nodes(node)%shape
      call coordinates(node)
      allocate(state%fields(node)%values(size(q),n(1),n(2),n(3)),temperature%fields(node)%values(1,n(1),n(2),n(3)))
      local=config;local%nx=n(1);local%ny=n(2);local%nz=n(3)
      call initialize_reactive_problem_3d(species,local,x,y,z,state%fields(node)%values, &
        temperature%fields(node)%values(1,:,:,:),base_density,mass_fractions,ok,mole_fractions)
      if (.not.ok) error stop 'Cannot initialize tree node'
    end do
    call average_down_amr_patch_tree_3d(tree,state,ok)
    if (.not.ok) error stop 'Cannot synchronize initial tree'
    do node=1,count+1
      n=tree%nodes(node)%shape;allocate(recovered(n(1),n(2),n(3)))
      call recover_reactive_temperatures_3d(species,state%fields(node)%values,temperature%fields(node)%values(1,:,:,:), &
        n(1),n(2),n(3),recovered,ok)
      if (.not.ok) error stop 'Invalid initial tree temperature'
      temperature%fields(node)%values(1,:,:,:)=recovered;deallocate(recovered)
    end do
    faces=resolved_reactive_boundaries_3d(config)
    viscosity=config%transport_enabled.and.config%viscosity_enabled
    conduction=config%transport_enabled.and.config%thermal_conduction_enabled
    diffusion=config%transport_enabled.and.config%species_diffusion_enabled
    barodiffusion=config%transport_enabled.and.config%barodiffusion_enabled
    periodic=[faces(1)=='periodic',faces(3)=='periodic',faces(5)=='periodic']
    min_h=minval(h,dim=2)
    if (regrid_interval_steps>0) then
      tag_spacing=root_h
      do i=1,refinement_levels
        tag_spacing=tag_spacing/real(refinement_ratios(i),dp)
      end do
      if (any(tag_spacing<=tiny(1.0_dp))) error stop 'Automatic tree spacing underflow'
      min_h=min(min_h,tag_spacing)
      ! Validate the requested integer lattice before any evolution/output.
      call tag_amr_patch_tree_temperature_3d(tree,temperature,temperature_thresholds(:refinement_levels), &
        refinement_ratios(:refinement_levels),tile_shape,tag_buffer,maximum_nodes,periodic,tagged_tree,ok)
      if (.not.ok) error stop 'Invalid initial automatic tree proposal'
    end if
    ! Reservoir bounds remain conservative for every scheduled refinement.
    do event=1,nevents
      call include_event_spacing(events(event))
    end do
    inlet=0;inlet_t=0;inlet_dt=huge(1.0_dp);mach_pressure=0
    where(config%inflow_mach>=0) mach_pressure=config%inflow_pressure
    if (any(faces=='inflow').and.config%inflow_composition_count/=0) then
      call resolve_selected_composition('Inflow',config%inflow_composition_count,config%inflow_composition_species, &
        config%inflow_composition_mole_fractions,species,inlet_x,ok,message)
      if (.not.ok) error stop 'Invalid inflow composition'
      call mass_fractions_from_mole_fractions(species,inlet_x,mass_fractions,ok)
      if (.not.ok) error stop 'Invalid inflow mass fractions'
    end if
    do face=1,6
      if (faces(face)/='inflow') cycle
      if(config%inflow_mach(face)>=0) then
        inlet_t(face)=config%inflow_temperature(face)
        cycle
      end if
      q=0;q(1)=mixture_density(species,mass_fractions,config%inflow_pressure(face),config%inflow_temperature(face),ok)
      if (.not.ok) error stop 'Invalid inflow density'
      q(2:4)=config%inflow_velocity(:,face);q(5)=config%inflow_pressure(face);q(6:)=mass_fractions
      call reactive_primitive_to_conserved(species,q,inlet(:,face),inlet_t(face),sound,ok)
      if (.not.ok) error stop 'Invalid inflow state'
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
        if (.not.ok) error stop 'Invalid inflow transport bound'
        inlet_dt=min(inlet_dt,face_dt)
      end if
    end do
    call composite_integrals_amr_patch_tree_3d(tree,state,root_h,initial,ok)
    if (.not.ok) error stop 'Invalid initial composite inventory'
    time=0;step=0;budget=0;max_correction=0;min_theta=1;event=1;regrids=0;max_remap_error=0
    if (len_trim(restart_file)>0) call restore_checkpoint()
    write(*,'(a)') 'PeleF '//pelef_version//' experimental serial regular 3D patch tree'
    write(*,'(a,l1)') 'Subcycling: ',subcycling
    write(*,'(a,i0,a,i0)') 'Nodes: ',count+1,' Levels: ',tree%level_count()
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
            call tag_amr_patch_tree_temperature_3d(tree,temperature,temperature_thresholds(:refinement_levels), &
              refinement_ratios(:refinement_levels),tile_shape,tag_buffer,maximum_nodes,periodic,tagged_tree,ok)
            if (.not.ok) error stop 'Automatic tree tagging failed'
            if (tree_key(tree)/=tree_key(tagged_tree)) call apply_event(tagged_tree)
            regrids=regrids+1;last_tag_step=step
            write(*,'(a,i0,a,i0)') 'Automatic regrid check at step ',step,' Nodes: ',size(tree%nodes)
          end if
        end if
      end if
      if (checkpoint_interval_steps>0) then
        if (step>last_checkpoint_step) then
          if (mod(step,checkpoint_interval_steps)==0) then
            call save_checkpoint()
            last_checkpoint_step=step
            if (stop_after_checkpoint) then
              write(*,'(a,i0,a,es24.16)') 'Stopped after checkpoint at step ',step,' time ',time
              return
            end if
          end if
        end if
      end if
      if (time>=config%final_time) exit
      stop_time=config%final_time
      if (event<=nevents) stop_time=event_times(event)
      if (step>=config%maximum_steps) error stop 'Tree maximum_steps reached'
      if (subcycling) then
        call compute_amr_patch_tree_subcycle_cfl_3d(tree,species,state,temperature,root_h,config%cfl,dt,ok)
      else
        call compute_amr_patch_tree_cfl_3d(tree,species,state,temperature,root_h,config%cfl,dt,ok)
      end if
      if (.not.ok) error stop 'Invalid tree hydro timestep'
      dt=min(dt,inlet_dt,stop_time-time)
      if (config%transport_enabled) then
        if (subcycling) then
          call compute_amr_patch_tree_subcycle_transport_timestep_3d(tree,species,transport,state,temperature,root_h, &
            config%transport_cfl,viscosity,conduction,diffusion,transport_dt,diffusivity,ok,faces,config%wall_temperature)
        else
          call compute_amr_patch_tree_transport_timestep_3d(tree,species,transport,state,temperature,root_h, &
            config%transport_cfl,viscosity,conduction,diffusion,transport_dt,diffusivity,ok,faces,config%wall_temperature)
        end if
        if (.not.ok) error stop 'Invalid tree transport timestep'
        dt=min(dt,transport_dt)
      end if
      if(any(config%inflow_mach>=0)) call limit_mach_step()
      do attempt=1,20
        if (time+dt<=time) error stop 'Tree timestep cannot advance time'
        if (sdc) then
          call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,state,temperature,root_h,dt,faces, &
            config%sdc_iterations,config%chemistry_enabled,[viscosity,conduction,diffusion,barodiffusion], &
            config%external_forcing,config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance, &
            change,source_change,correction,theta,ok,reaction_source=reaction_source, &
            inflow_state=inlet,inflow_temperature=inlet_t,wall_temperature=config%wall_temperature, &
            wall_velocity=config%wall_velocity,outlet_pressure=config%outlet_pressure, &
            inflow_mach=config%inflow_mach,inflow_pressure=mach_pressure,difmag=resolved_reactive_difmag_3d(config))
        else
          call advance_amr_patch_tree_split_3d(tree,species,reactions,transport,state,temperature,root_h,dt,faces, &
            config%reconstruction,config%limiter,config%riemann_solver,config%chemistry_enabled, &
            viscosity,conduction,diffusion,barodiffusion, &
            config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance, &
            change,correction,theta,ok,chemistry_integrator=chemistry_integrator,inflow_state=inlet,inflow_temperature=inlet_t, &
            wall_temperature=config%wall_temperature,wall_velocity=config%wall_velocity,subcycling=subcycling, &
            outlet_pressure=config%outlet_pressure,inflow_mach=config%inflow_mach,inflow_pressure=mach_pressure, &
            difmag=resolved_reactive_difmag_3d(config))
        end if
        if (ok) exit
        dt=.5_dp*dt
      end do
      if (.not.ok) error stop 'Tree step failed after 20 attempts'
      if (dt==stop_time-time) then
        time=stop_time
      else
        time=time+dt
      end if
      step=step+1;budget=budget+change
      source_budget=source_budget+source_change
      max_correction=max(max_correction,correction);min_theta=min(min_theta,theta)
    end do
    call composite_integrals_amr_patch_tree_3d(tree,state,root_h,final,ok)
    if (.not.ok) error stop 'Invalid final composite inventory'
    open(newunit=unit,file=trim(config%output_file),status='replace',action='write',iostat=status)
    if (status/=0) error stop 'Cannot open tree CSV'
    write(unit,'(a)',iostat=status) 'node,level,dx,dy,dz,'//reactive_3d_csv_header(species)
    if (status/=0) error stop 'Cannot write tree CSV header'
    allocate(character(len=26*(14+2*size(species))-1)::row)
    do node=1,count+1
      n=tree%nodes(node)%shape;call coordinates(node)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            covered=.false.
            do child=node+1,count+1
              if (tree%nodes(child)%parent/=node) cycle
              lo=[tree%nodes(child)%patch%coarse_i_lower,tree%nodes(child)%patch%coarse_j_lower, &
                tree%nodes(child)%patch%coarse_k_lower]
              hi=[tree%nodes(child)%patch%coarse_i_upper,tree%nodes(child)%patch%coarse_j_upper, &
                tree%nodes(child)%patch%coarse_k_upper]
              if (all([i,j,k]>=lo).and.all([i,j,k]<=hi)) covered=.true.
            end do
            if (covered) cycle
            call format_reactive_3d_csv_row(species,[x(i),y(j),z(k)],time,state%fields(node)%values(:,i,j,k), &
              temperature%fields(node)%values(1,i,j,k),row,ok)
            if (.not.ok) error stop 'Invalid tree output state'
            write(unit,'(i0,",",i0,3(",",es25.16e3),",",a)',iostat=status) node,tree%nodes(node)%level,h(:,node),row
            if (status/=0) error stop 'Cannot write tree CSV row'
          end do
        end do
      end do
    end do
    close(unit,iostat=status)
    if (status/=0) error stop 'Cannot close tree CSV'
    open(newunit=unit,file=trim(config%output_file)//'.budget.csv',status='replace',action='write',iostat=status)
    if (status/=0) error stop 'Cannot open tree budget CSV'
    if (sdc) then
      write(unit,'(a)',iostat=status) 'component,initial,final,boundary_change,source_change'
    else
      write(unit,'(a)',iostat=status) 'component,initial,final,boundary_change'
    end if
    if (status/=0) error stop 'Cannot write tree budget header'
    do c=1,size(q)
      if (sdc) then
        write(unit,'(i0,4(",",es25.16e3))',iostat=status) c,initial(c),final(c),budget(c),source_budget(c)
      else
        write(unit,'(i0,3(",",es25.16e3))',iostat=status) c,initial(c),final(c),budget(c)
      end if
      if (status/=0) error stop 'Cannot write tree budget row'
    end do
    close(unit,iostat=status)
    if (status/=0) error stop 'Cannot close tree budget CSV'
    write(*,'(a,i0,a,es12.4)') 'Regrids: ',regrids,' Maximum remap error: ',max_remap_error
    write(*,'(a,i0,a,es24.16)') 'Steps: ',step,' Final time: ',time
    write(*,'(a,es12.4,a,es12.4)') 'Maximum reflux: ',max_correction,' Minimum transport theta: ',min_theta
  contains
    subroutine limit_mach_step()
      type(amr_tree_state_3d) :: packed
      real(dp), allocatable :: ghost(:,:,:,:)
      real(dp) :: bound,dmax,factor
      integer :: peer,dims(3),nv,width
      character(len=:), allocatable :: donor_limiter
      nv=size(q)
      if(config%reconstruction/='pcm') donor_limiter=config%limiter
      width=2
      if(config%reconstruction=='characteristic_ppm') width=4
      allocate(packed%fields(size(tree%nodes)))
      do peer=1,size(tree%nodes)
        dims=tree%nodes(peer)%shape
        allocate(packed%fields(peer)%values(nv+1,dims(1),dims(2),dims(3)))
        packed%fields(peer)%values(:nv,:,:,:)=state%fields(peer)%values
        packed%fields(peer)%values(nv+1,:,:,:)=temperature%fields(peer)%values(1,:,:,:)
      end do
      do peer=1,size(tree%nodes)
        call fill_amr_patch_tree_ghosts_3d(tree,packed,peer,faces,ghost,ok,viscous_wall=.true., &
          inflow_state=inlet_with_temperature(),wall_velocity=config%wall_velocity,coarse_limiter=donor_limiter, &
          species_count=size(species),ghost_width=width,species=species,inflow_mach=config%inflow_mach, &
          inflow_temperature=inlet_t,inflow_pressure=mach_pressure)
        if(.not.ok) error stop 'Invalid tree Mach inlet donors'
        dims=shape(ghost(1,:,:,:));factor=1
        if(subcycling) factor=root_h(1)/h(1,peer)
        call compute_reactive_cfl_timestep_3d(species,ghost(:nv,:,:,:),ghost(nv+1,:,:,:), &
          dims(1),dims(2),dims(3),h(1,peer),h(2,peer),h(3,peer),config%cfl,bound,ok)
        if(.not.ok) error stop 'Invalid tree Mach hydro timestep'
        dt=min(dt,factor*bound)
        if(config%transport_enabled) then
          call reactive_transport_timestep_3d(species,transport,ghost(:nv,:,:,:),ghost(nv+1,:,:,:), &
            dims(1),dims(2),dims(3),h(1,peer),h(2,peer),h(3,peer),config%transport_cfl, &
            viscosity,conduction,diffusion,bound,dmax,ok)
          if(.not.ok) error stop 'Invalid tree Mach transport timestep'
          dt=min(dt,factor**2*bound)
        end if
      end do
    end subroutine
    function inlet_with_temperature() result(packed_inlet)
      real(dp) :: packed_inlet(size(q)+1,6)
      packed_inlet(:size(q),:)=inlet;packed_inlet(size(q)+1,:)=inlet_t
    end function
    subroutine check_paths(path)
      character(len=*), intent(in) :: path
      if (checkpoint_interval_steps>0) then
        if (selected_reactor_paths_alias_lexically(checkpoint_file,path)) error stop 'Tree checkpoint aliases input/output'
      end if
      if (len_trim(restart_file)>0) then
        if (selected_reactor_paths_alias_lexically(restart_file,path)) error stop 'Tree restart aliases input/output'
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
    subroutine save_checkpoint()
      snapshot%tree=tree;snapshot%state=state;snapshot%temperature=temperature
      snapshot%initial=initial;snapshot%boundary=budget
      if (sdc) snapshot%source=source_budget
      if (allocated(reaction_source)) snapshot%reaction_source=reaction_source
      snapshot%time=time;snapshot%steps=step;snapshot%regrids=regrids
      snapshot%maximum_reflux=max_correction;snapshot%minimum_theta=min_theta
      snapshot%maximum_remap_error=max_remap_error
      call write_amr_tree_checkpoint_3d(trim(checkpoint_file),species,context,snapshot,ok,indexed=checkpoint_indexed)
      if (.not.ok) error stop 'Cannot write tree checkpoint'
    end subroutine
    subroutine restore_checkpoint()
      integer :: applied,peer,nconserved
      type(amr_tree_state_3d) :: absolute_state
      real(dp) :: original_scale(size(q)),restart_scale(size(q)),inventory(size(q))
      call read_amr_tree_checkpoint_3d(trim(restart_file),species,context, &
        [config%nx,config%ny,config%nz],snapshot,ok)
      if (.not.ok) error stop 'Invalid tree restart checkpoint/context'
      if (sdc.neqv.allocated(snapshot%source)) error stop 'Tree restart source ledger mismatch'
      if (allocated(snapshot%reaction_source%fields).neqv.(sdc.and.config%chemistry_enabled)) &
        error stop 'Tree restart reaction source mismatch'
      if (sdc) source_budget=snapshot%source
      if (allocated(reaction_source)) reaction_source=snapshot%reaction_source
      if (snapshot%time>config%final_time.or.snapshot%steps>config%maximum_steps) &
        error stop 'Tree restart exceeds run limits'
      applied=count_events(snapshot%time)
      if (regrid_interval_steps>0) applied=snapshot%steps/regrid_interval_steps
      if (snapshot%regrids/=applied) error stop 'Tree restart schedule position mismatch'
      if (regrid_interval_steps>0.and.applied>0) then
        call validate_adaptive_geometry(snapshot%tree)
      else
        expected_tree=tree_key(tree)
        if (applied>0) expected_tree=tree_key(events(applied))
        if (tree_key(snapshot%tree)/=expected_tree) error stop 'Tree restart geometry mismatch'
      end if
      if (any(abs(snapshot%initial-initial)>64*epsilon(1.0_dp)*max(tiny(1.0_dp),abs(initial)))) &
        error stop 'Tree restart initial inventory mismatch'
      ! Use absolute cell inventories so cancelling momentum has a meaningful scale.
      absolute_state=state
      do peer=1,size(tree%nodes)
        absolute_state%fields(peer)%values=abs(state%fields(peer)%values)
      end do
      call composite_integrals_amr_patch_tree_3d(tree,absolute_state,root_h,original_scale,ok)
      if (.not.ok) error stop 'Invalid initial restart inventory scale'
      absolute_state=snapshot%state
      do peer=1,size(snapshot%tree%nodes)
        absolute_state%fields(peer)%values=abs(snapshot%state%fields(peer)%values)
      end do
      call composite_integrals_amr_patch_tree_3d(snapshot%tree,absolute_state,root_h,restart_scale,ok)
      if (.not.ok) error stop 'Invalid restart inventory scale'
      call composite_integrals_amr_patch_tree_3d(snapshot%tree,snapshot%state,root_h,inventory,ok)
      if (.not.ok) error stop 'Invalid restart inventory'
      ! Momentum is a vector: a nearly zero transverse component must use
      ! the full momentum magnitude as its roundoff scale.
      original_scale(2:4)=maxval(original_scale(2:4))
      restart_scale(2:4)=maxval(restart_scale(2:4))
      nconserved=size(q)
      if (config%chemistry_enabled) nconserved=5
      if (any(abs(inventory(:nconserved)-snapshot%initial(:nconserved)- &
          snapshot%boundary(:nconserved)-source_budget(:nconserved))> &
          2e-10_dp*max(tiny(1.0_dp),original_scale(:nconserved),restart_scale(:nconserved), &
          abs(snapshot%boundary(:nconserved)),abs(source_budget(:nconserved))))) error stop 'Tree restart boundary ledger mismatch'
      tree=snapshot%tree;state=snapshot%state;temperature=snapshot%temperature
      initial=snapshot%initial;budget=snapshot%boundary
      time=snapshot%time;step=snapshot%steps;regrids=snapshot%regrids
      max_correction=snapshot%maximum_reflux;min_theta=snapshot%minimum_theta
      max_remap_error=snapshot%maximum_remap_error
      event=applied+1;last_checkpoint_step=step;last_tag_step=step;count=size(tree%nodes)-1
      call amr_tree_integer_layout_3d(tree,origin,domain,ratios,ok)
      if (.not.ok) error stop 'Tree restart coordinate overflow'
      deallocate(h);allocate(h(3,count+1))
      do peer=1,count+1
        h(:,peer)=root_h*real(domain(:,1),dp)/real(domain(:,tree%nodes(peer)%level+1),dp)
      end do
    end subroutine
    subroutine validate_adaptive_geometry(geometry)
      type(amr_patch_tree_3d), intent(in) :: geometry
      integer :: peer,depth,parent,lower(3),upper(3),dims(3)
      if (size(geometry%nodes)>maximum_nodes.or.geometry%level_count()>refinement_levels+1) &
        error stop 'Tree restart exceeds adaptation bounds'
      do peer=2,size(geometry%nodes)
        depth=geometry%nodes(peer)%level;parent=geometry%nodes(peer)%parent
        dims=geometry%nodes(parent)%shape
        associate(box=>geometry%nodes(peer)%patch)
          lower=[box%coarse_i_lower,box%coarse_j_lower,box%coarse_k_lower]
          upper=[box%coarse_i_upper,box%coarse_j_upper,box%coarse_k_upper]
          if (box%refinement_ratio/=refinement_ratios(depth)) error stop 'Tree restart adaptation ratio mismatch'
        end associate
        if (any(mod(lower-1,tile_shape)/=0)) error stop 'Tree restart tile alignment mismatch'
        if (any(int(upper,int64)/=min(int(dims,int64),int(lower,int64)+int(tile_shape,int64)-1))) &
          error stop 'Tree restart tile extent mismatch'
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
    subroutine read_topology(path,result)
      character(len=*), intent(in) :: path
      type(amr_patch_tree_3d), intent(out) :: result
      type(amr_patch_3d), allocatable :: patches(:)
      integer, allocatable :: parents(:),shapes(:,:)
      integer :: input,ios,nchild,peer,parent,dims(3),lower(3),upper(3),ratio
      logical :: valid
      call check_paths(path)
      open(newunit=input,file=path,status='old',action='read',iostat=ios)
      if (ios/=0) error stop 'Cannot open tree topology'
      read(input,*,iostat=ios) nchild
      if (ios/=0) error stop 'Invalid tree child count'
      if (nchild<0.or.nchild==huge(nchild)) error stop 'Invalid tree child count'
      allocate(patches(nchild),parents(nchild),shapes(3,nchild+1),stat=ios)
      if (ios/=0) error stop 'Cannot allocate tree topology'
      shapes(:,1)=[config%nx,config%ny,config%nz]
      do peer=1,nchild
        read(input,*,iostat=ios) parent,lower(1),upper(1),lower(2),upper(2),lower(3),upper(3),ratio
        if (ios/=0) error stop 'Invalid tree patch record'
        if (parent<1.or.parent>peer) error stop 'Tree parents must precede children'
        parents(peer)=parent;dims=shapes(:,parent)
        call initialize_amr_patch_3d(dims(1),dims(2),dims(3),lower(1),upper(1),lower(2),upper(2), &
          lower(3),upper(3),ratio,patches(peer),valid)
        if (.not.valid) error stop 'Invalid tree patch geometry'
        if (any(int(upper-lower+1,int64)*int(ratio,int64)>int(huge(1),int64))) error stop 'Tree patch extent overflow'
        shapes(:,peer+1)=[patches(peer)%fine_nx(),patches(peer)%fine_ny(),patches(peer)%fine_nz()]
      end do
      close(input)
      call initialize_amr_patch_tree_3d(shapes(:,1),parents,patches,result,valid)
      if (.not.valid) error stop 'Invalid tree topology'
    end subroutine
    subroutine include_event_spacing(next_tree)
      type(amr_patch_tree_3d), intent(in) :: next_tree
      integer(int64), allocatable :: o(:,:),d(:,:)
      integer, allocatable :: r(:)
      logical :: valid
      call amr_tree_integer_layout_3d(next_tree,o,d,r,valid)
      if (.not.valid) error stop 'Scheduled tree coordinate overflow'
      min_h=min(min_h,root_h*real(d(:,1),dp)/real(d(:,size(r)),dp))
    end subroutine
    subroutine apply_event(next_tree)
      type(amr_patch_tree_3d), intent(in) :: next_tree
      type(amr_tree_state_3d) :: candidate,temps,absolute_state
      integer(int64), allocatable :: o(:,:),d(:,:)
      integer, allocatable :: r(:)
      real(dp), allocatable :: spacing(:,:)
      real(dp) :: before(size(q)),after(size(q)),scale(size(q)),error
      integer :: peer
      logical :: valid
      call remap_amr_patch_tree_reactive_3d(tree,state,temperature,next_tree,species,candidate,temps,valid)
      if (.not.valid) error stop 'Scheduled tree remap failed'
      call composite_integrals_amr_patch_tree_3d(tree,state,root_h,before,valid)
      if (.not.valid) error stop 'Invalid pre-remap inventory'
      call composite_integrals_amr_patch_tree_3d(next_tree,candidate,root_h,after,valid)
      if (.not.valid) error stop 'Invalid post-remap inventory'
      absolute_state=state
      do peer=1,size(tree%nodes)
        absolute_state%fields(peer)%values=abs(state%fields(peer)%values)
      end do
      call composite_integrals_amr_patch_tree_3d(tree,absolute_state,root_h,scale,valid)
      if (.not.valid) error stop 'Invalid remap inventory scale'
      error=maxval(abs(after-before)/max(tiny(1.0_dp),scale))
      if (error>5e-12_dp) error stop 'Scheduled tree remap lost conservation'
      call amr_tree_integer_layout_3d(next_tree,o,d,r,valid)
      if (.not.valid) error stop 'Scheduled tree coordinate overflow'
      allocate(spacing(3,size(next_tree%nodes)))
      do peer=1,size(next_tree%nodes)
        spacing(:,peer)=root_h*real(d(:,1),dp)/real(d(:,next_tree%nodes(peer)%level+1),dp)
      end do
      tree=next_tree;state=candidate;temperature=temps
      origin=o;domain=d;ratios=r;h=spacing;count=size(tree%nodes)-1
      max_remap_error=max(max_remap_error,error)
      write(*,'(a,es24.16,a,i0,a,i0)') 'Regrid time: ',time,' Nodes: ',count+1,' Levels: ',tree%level_count()
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
end module amr_patch_tree_application_3d_mod
