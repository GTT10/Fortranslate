module amr_patch_tree_application_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use dual_energy_mod, only: dual_energy_policy
  use amr_patch_tree_tagging_3d_mod, only: tag_amr_patch_tree_temperature_3d
  use amr_patch_tree_remap_3d_mod, only: remap_amr_patch_tree_reactive_3d
  use, intrinsic :: iso_fortran_env, only: int64
  use amr_patch_tree_checkpoint_3d_mod, only: amr_tree_checkpoint_3d, build_amr_tree_checkpoint_context_3d, &
    read_amr_tree_checkpoint_3d,write_amr_tree_checkpoint_3d
  use simulation_config_selected_reactor_mod, only: selected_reactor_paths_alias_lexically
  use precision_mod, only: dp
  use rotating_frame_3d_mod, only: rotating_frame_3d,convert_rotating_energy_3d
  use reactive_mol_3d_mod, only: convert_regular_energy_3d
  use constant_volume_reactor_mod, only: forced_reactor_integrator
  use amr_patch_tree_chemistry_3d_mod, only: initialize_amr_tree_reaction_source_3d
  use amr_patch_tree_subcycle_3d_mod, only: compute_amr_patch_tree_subcycle_cfl_3d
  use amr_patch_tree_transport_subcycle_3d_mod, only: compute_amr_patch_tree_subcycle_transport_timestep_3d
  use constants_mod, only: pelef_version
  use nasa7_thermo_mod, only: nasa7_species,configure_nasa7_eos_policy
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
  use amr_patch_tree_hydro_3d_mod, only: compute_amr_patch_tree_cfl_3d,convert_amr_patch_tree_rotating_energy_3d
  use amr_patch_tree_transport_3d_mod, only: compute_amr_patch_tree_transport_timestep_3d
  use amr_patch_tree_split_3d_mod, only: advance_amr_patch_tree_split_3d,advance_amr_patch_tree_sdc_3d, &
    advance_amr_patch_tree_sdc_subcycle_3d
  use reactive_csv_io_3d_mod, only: reactive_3d_csv_header,format_reactive_3d_csv_row
  use simulation_config_reactive_3d_mod, only: resolved_reactive_difmag_3d
  use amr_patch_tree_split_3d_mod, only: amr_sdc_diffusion_evaluator
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use reactive_3d_application_mod, only: read_constant_transport,read_reactive_rotation_controls
  implicit none
  private
  public :: run_amr_patch_tree_application_3d
  abstract interface
    subroutine tree_application_diffusion(tree,u,t,h,faces,inlet,inlet_t,depth,fractions, &
        coarse_start,coarse_end,rates,flux,diffusivity,ok,frame,domain_lower,internal_energy)
      import amr_patch_tree_3d,amr_tree_state_3d,amr_tree_flux_3d,dp,rotating_frame_3d
      type(amr_tree_state_3d), optional, intent(in) :: internal_energy
      type(rotating_frame_3d), optional, intent(in) :: frame
      real(dp), optional, intent(in) :: domain_lower(3)
      type(amr_patch_tree_3d), intent(in) :: tree
      integer, intent(in) :: depth
      type(amr_tree_state_3d), intent(in) :: u,t
      real(dp), intent(in) :: h(3),inlet(:,:),inlet_t(:)
      character(len=*), intent(in) :: faces(6)
      real(dp), optional, intent(in) :: fractions(:)
      type(amr_tree_state_3d), optional, intent(in) :: coarse_start,coarse_end
      type(amr_tree_state_3d), intent(inout) :: rates
      type(amr_tree_flux_3d), allocatable, intent(inout) :: flux(:)
      real(dp), intent(inout) :: diffusivity
      logical, intent(out) :: ok
    end subroutine
  end interface
contains
  subroutine run_amr_patch_tree_application_3d(topology_path,config,species,reactions,transport,mole_fractions, &
      schedule_path,chemistry_integrator,input_path,mechanism_identity,evaluate_diffusion,diffusion_identity,evaluate_reaction)
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    procedure(tree_application_diffusion), optional :: evaluate_diffusion
    character(len=*), optional, intent(in) :: diffusion_identity
    type(amr_tree_state_3d) :: diffusion_rates
    type(amr_tree_flux_3d), allocatable :: diffusion_flux(:)
    character(len=*), intent(in) :: topology_path
    character(len=*), optional, intent(in) :: schedule_path,chemistry_integrator,input_path,mechanism_identity
    type(reactive_3d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    type(nasa7_species), allocatable :: eos_species(:)
    logical :: eos_policy_ok
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: mole_fractions(:)
    type(reactive_3d_config) :: local
    type(amr_patch_tree_3d) :: tree,tagged_tree
    type(amr_tree_checkpoint_3d) :: snapshot
    type(amr_patch_tree_3d), allocatable :: events(:)
    type(amr_tree_state_3d) :: state,temperature,ordinary_state
    type(rotating_frame_3d), allocatable :: frame
    real(dp), allocatable :: rotation_lower(:)
    type(amr_tree_state_3d), allocatable :: reaction_source,heat_release,internal_energy
    type(dual_energy_policy), allocatable :: energy_policy
    real(dp), allocatable :: energy_correction,pressure_work_change,step_pressure_work
    real(dp), allocatable :: held_cell
    real(dp) :: repair
    integer, allocatable :: ratios(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    real(dp), allocatable :: x(:),y(:),z(:),root_h(:),h(:,:),recovered(:,:,:),mass_fractions(:),inlet_x(:)
    real(dp), allocatable :: initial(:),final(:),budget(:),change(:),inlet(:,:),q(:),inlet_cube(:,:,:,:)
    real(dp), allocatable :: event_times(:),source_budget(:),source_change(:)
    real(dp) :: stop_time,max_remap_error
    real(dp) :: time,dt,transport_dt,diffusivity,base_density,correction,theta,max_correction,min_theta,ordinary_energy
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
    logical :: isothermal_slip_faces(6)
    logical :: periodic(3),stop_after_checkpoint,subcycling
    logical :: checkpoint_indexed,sdc,heat_known,mol,initial_reactions_pending
    real(dp), allocatable :: constant_coefficients(:),hit_transport_parameters(:)
    character(len=:), allocatable :: constant_identity
    namelist /amr_tree_checkpoint/ checkpoint_file,restart_file,checkpoint_interval_steps,stop_after_checkpoint,checkpoint_indexed
    namelist /amr_tree_adaptation/ regrid_interval_steps,refinement_levels,refinement_ratios, &
      temperature_thresholds,tile_shape,tag_buffer,maximum_nodes
    real(dp) :: maximum_dt,initial_dt,fixed_dt,initial_timestep_factor,maximum_timestep_growth,previous_timestep
    logical :: post_regrid
    namelist /amr_tree_timestepping/ subcycling,maximum_dt,initial_dt,fixed_dt, &
      initial_timestep_factor,maximum_timestep_growth,isothermal_slip_faces
    eos_species=species
    call configure_nasa7_eos_policy(eos_species,config%pinned_eos_extrapolation,eos_policy_ok, &
      config%pinned_eos_temperature_min,config%pinned_eos_temperature_max)
    if(.not.eos_policy_ok) error stop 'Invalid application EOS policy'
    mol=trim(config%time_integrator)=='mol'
    sdc=trim(config%time_integrator)=='sdc'.or.mol
    if(present(input_path)) call read_reactive_rotation_controls(input_path,config%time_integrator,frame)
    if(allocated(frame)) then
      allocate(rotation_lower(3))
      rotation_lower=[config%x_lower,config%y_lower,config%z_lower]
    endif
    heat_known=.false.
    if((config%sdc_reactor=='cvode').neqv.present(evaluate_reaction)) error stop 'SDC reactor selection does not match backend'
    if(present(evaluate_diffusion).neqv.present(diffusion_identity)) error stop 'Diffusion callback requires identity'
    if(present(evaluate_diffusion).and..not.sdc) error stop 'Additional tree diffusion requires SDC'
    constant_identity=''
    if(present(input_path)) &
      call read_constant_transport(input_path,constant_coefficients,hit_transport_parameters)
    if(allocated(hit_transport_parameters)) error stop 'HIT transport scales require initialized regular HIT input'
    if(allocated(constant_coefficients)) then
      if(.not.sdc.or..not.config%transport_enabled) error stop 'Tree constant transport requires enabled SDC transport'
      constant_identity='constant-transport:'
      do event=1,size(constant_coefficients)
        write(message,'(es27.18e3)') constant_coefficients(event)
        constant_identity=constant_identity//trim(message)//new_line('a')
      enddo
    endif
    isothermal_slip_faces=.false.
    subcycling=.false.
    maximum_dt=huge(1._dp);initial_dt=0;fixed_dt=0;initial_timestep_factor=1;maximum_timestep_growth=1.1_dp
    previous_timestep=0;post_regrid=.false.
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
    if(any(isothermal_slip_faces)) then
      faces=resolved_reactive_boundaries_3d(config)
      if(.not.sdc.or..not.config%transport_enabled.or..not.config%thermal_conduction_enabled) &
        error stop 'Thermal slip requires SDC/MOL conduction'
      if(any(isothermal_slip_faces.and.faces/='slip_wall')) error stop 'Thermal slip requires slip_wall faces'
      if(any(isothermal_slip_faces.and.(.not.ieee_is_finite(config%wall_temperature).or.config%wall_temperature<=0))) &
        error stop 'Invalid thermal slip temperature'
    endif
    if (.not.sdc.and.any(config%external_forcing/=0)) error stop 'Tree external forcing requires SDC'
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
    if(.not.all(ieee_is_finite([maximum_dt,initial_dt,fixed_dt,initial_timestep_factor,maximum_timestep_growth]))) &
      error stop 'Nonfinite tree timestep controls'
    if(maximum_dt<=0.or.fixed_dt>maximum_dt.or.initial_timestep_factor<=0.or.initial_timestep_factor>1.or. &
      maximum_timestep_growth<1) error stop 'Invalid tree timestep controls'
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
    if(mol) then
      write(message,'(a,l1)') 'regular-tree-mol-v4 ',config%mol_lax_friedrichs
      schedule_identity=schedule_identity//trim(message)//new_line('a')
    endif
    write(message,'(a,5(es27.18e3,1x))') 'tree-timestep-v1 ',maximum_dt,initial_dt,fixed_dt, &
      initial_timestep_factor,maximum_timestep_growth
    schedule_identity=schedule_identity//trim(message)//new_line('a')
    if(mol) schedule_identity=schedule_identity//'AMR_MOL_AUXILIARY_ENERGY_V1'//new_line('a')
    if(sdc.and..not.mol) schedule_identity=schedule_identity//'AMR_GODUNOV_AUXILIARY_ENERGY_V1'//new_line('a')
    if(allocated(frame)) then
      write(message,'(i0,4(1x,es27.18e3))') frame%axis,frame%angular_speed,frame%axis_point
      schedule_identity=schedule_identity//'tree-rotation-v1 '//trim(message)//new_line('a')
    endif
    if (checkpoint_interval_steps>0.or.len_trim(restart_file)>0) then
      identity='fixed:'//trim(config%thermo_model)
      if (present(mechanism_identity)) identity=mechanism_identity
      call resolve_reactive_chemistry_integrator(size(eos_species),chemistry_integrator,use_implicit,ok)
      if (.not.ok) error stop 'Invalid tree chemistry integrator'
      integrator='explicit'
      if (use_implicit) integrator='implicit'
      call build_amr_tree_checkpoint_context_3d(config,eos_species,reactions,transport,mole_fractions, &
        identity,integrator,schedule_identity,context,ok)
      if (.not.ok) error stop 'Cannot build tree checkpoint context'
      if(present(diffusion_identity)) context=context//'diffusion:'//diffusion_identity//new_line('a')
      context=context//constant_identity
      if(any(isothermal_slip_faces)) then
        write(message,'(6l1)') isothermal_slip_faces
        context=context//'isothermal-slip:'//trim(message)//new_line('a')
      endif
    end if
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,ok)
    if (.not.ok) error stop 'Tree coordinate overflow'
    root_h=[(config%x_upper-config%x_lower)/config%nx,(config%y_upper-config%y_lower)/config%ny, &
      (config%z_upper-config%z_lower)/config%nz]
    allocate(h(3,count+1),state%fields(count+1),temperature%fields(count+1))
    allocate(mass_fractions(size(eos_species)),inlet_x(size(eos_species)),q(reactive_nvar(size(eos_species))))
    allocate(initial(size(q)),final(size(q)),budget(size(q)),change(size(q)),inlet(size(q),6),inlet_cube(size(q),2,2,2))
    allocate(source_budget(size(q)),source_change(size(q)))
    source_budget=0;source_change=0
    if(sdc) then
      allocate(heat_release)
      allocate(heat_release%fields(count+1))
      do node=1,count+1
        n=tree%nodes(node)%shape
        allocate(heat_release%fields(node)%values(1,n(1),n(2),n(3)),source=0._dp)
      end do
      heat_known=.true.
    end if
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
      call initialize_reactive_problem_3d(eos_species,local,x,y,z,state%fields(node)%values, &
        temperature%fields(node)%values(1,:,:,:),base_density,mass_fractions,ok,mole_fractions)
      if (.not.ok) error stop 'Cannot initialize tree node'
    end do
    call average_down_amr_patch_tree_3d(tree,state,ok)
    if (.not.ok) error stop 'Cannot synchronize initial tree'
    do node=1,count+1
      n=tree%nodes(node)%shape;allocate(recovered(n(1),n(2),n(3)))
      call recover_reactive_temperatures_3d(eos_species,state%fields(node)%values,temperature%fields(node)%values(1,:,:,:), &
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
        config%inflow_composition_mole_fractions,eos_species,inlet_x,ok,message)
      if (.not.ok) error stop 'Invalid inflow composition'
      call mass_fractions_from_mole_fractions(eos_species,inlet_x,mass_fractions,ok)
      if (.not.ok) error stop 'Invalid inflow mass fractions'
    end if
    do face=1,6
      if (faces(face)/='inflow') cycle
      if(config%inflow_mach(face)>=0) then
        inlet_t(face)=config%inflow_temperature(face)
        cycle
      end if
      q=0;q(1)=mixture_density(eos_species,mass_fractions,config%inflow_pressure(face),config%inflow_temperature(face),ok)
      if (.not.ok) error stop 'Invalid inflow density'
      q(2:4)=config%inflow_velocity(:,face);q(5)=config%inflow_pressure(face);q(6:)=mass_fractions
      call reactive_primitive_to_conserved(eos_species,q,inlet(:,face),inlet_t(face),sound,ok)
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
        call reactive_transport_timestep_3d(eos_species,transport,inlet_cube,inlet_temperature,2,2,2, &
          min_h(1),min_h(2),min_h(3),config%transport_cfl,viscosity,conduction,diffusion,face_dt,diffusivity,ok, &
          constant_coefficients=constant_coefficients)
        if (.not.ok) error stop 'Invalid inflow transport bound'
        inlet_dt=min(inlet_dt,face_dt)
      end if
    end do
    if(sdc) then
      allocate(internal_energy,energy_policy,energy_correction,pressure_work_change,step_pressure_work)
      internal_energy=temperature;energy_correction=0;pressure_work_change=0;step_pressure_work=0
      energy_policy%allow_small_energy=config%allow_small_energy
      energy_policy%allow_negative_energy=config%allow_negative_energy
      energy_policy%update_total_from_internal=config%dual_energy_update_E_from_e
      energy_policy%eta2=config%dual_energy_eta2
      do node=1,size(tree%nodes)
        internal_energy%fields(node)%values(1,:,:,:)=state%fields(node)%values(5,:,:,:)- &
          .5_dp*sum(state%fields(node)%values(2:4,:,:,:)**2,dim=1)/state%fields(node)%values(1,:,:,:)
      enddo
    endif
    if(allocated(frame)) then
      call convert_amr_patch_tree_rotating_energy_3d(tree,state,rotation_lower,root_h,frame,.true.,ok)
      if(.not.ok) error stop 'Invalid initial rotating tree energy'
    endif
    call composite_integrals_amr_patch_tree_3d(tree,state,root_h,initial,ok)
    if (.not.ok) error stop 'Invalid initial composite inventory'
    time=0;step=0;budget=0;max_correction=0;min_theta=1;event=1;regrids=0;max_remap_error=0
    if (len_trim(restart_file)>0) call restore_checkpoint()
    initial_reactions_pending=mol.and.config%chemistry_enabled.and.len_trim(restart_file)==0
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
            post_regrid=.true.
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
      if(step==0.and.initial_dt>0) then
        dt=initial_dt
      else if(fixed_dt>0) then
        dt=fixed_dt
      else
        ordinary_state=state
        if(allocated(frame)) then
          call convert_amr_patch_tree_rotating_energy_3d(tree,ordinary_state,rotation_lower,root_h,frame,.false.,ok)
          if(.not.ok) error stop 'Invalid rotating tree timestep state'
        endif
        if (subcycling) then
          call compute_amr_patch_tree_subcycle_cfl_3d(tree,eos_species,ordinary_state,temperature,root_h,config%cfl,dt,ok, &
            internal_energy=internal_energy)
        else
          call compute_amr_patch_tree_cfl_3d(tree,eos_species,ordinary_state,temperature,root_h,config%cfl,dt,ok, &
            internal_energy=internal_energy)
        end if
        if (.not.ok) error stop 'Invalid tree hydro timestep'
        dt=min(dt,inlet_dt,maximum_dt)
        if (config%transport_enabled) then
          if (subcycling) then
            call compute_amr_patch_tree_subcycle_transport_timestep_3d(tree,eos_species,transport, &
      ordinary_state,temperature,root_h, &
              config%transport_cfl,viscosity,conduction,diffusion,transport_dt,diffusivity,ok,faces,config%wall_temperature, &
            constant_coefficients=constant_coefficients,internal_energy=internal_energy, &
            isothermal_slip_faces=isothermal_slip_faces)
          else
            call compute_amr_patch_tree_transport_timestep_3d(tree,eos_species,transport,ordinary_state,temperature,root_h, &
              config%transport_cfl,viscosity,conduction,diffusion,transport_dt,diffusivity,ok,faces,config%wall_temperature, &
            constant_coefficients=constant_coefficients,internal_energy=internal_energy, &
            isothermal_slip_faces=isothermal_slip_faces)
          end if
          if (.not.ok) error stop 'Invalid tree transport timestep'
          dt=min(dt,transport_dt)
        end if
        if(present(evaluate_diffusion)) then
          call evaluate_diffusion(tree,state,temperature,root_h,faces,inlet,inlet_t,-1, &
            rates=diffusion_rates,flux=diffusion_flux,diffusivity=diffusivity,ok=ok, &
            frame=frame,domain_lower=rotation_lower,internal_energy=internal_energy)
          if(.not.ok) error stop 'Invalid additional tree diffusion timestep'
          if(.not.ieee_is_finite(diffusivity).or.diffusivity<0) error stop 'Invalid additional diffusivity'
          if(diffusivity>0) dt=min(dt,.1_dp/(diffusivity*sum(1._dp/minval(h,dim=2)**2)))
        end if
        if(any(config%inflow_mach>=0)) call limit_mach_step()
      endif
      if(step==0) then
        if(initial_dt<=0) dt=initial_timestep_factor*dt
      else if(fixed_dt<=0) then
        if(post_regrid) then
          dt=min(dt,previous_timestep)
        else if(previous_timestep<dt/maximum_timestep_growth) then
          dt=maximum_timestep_growth*previous_timestep
        endif
      endif
      if(config%final_time-time-dt<.001_dp*dt) dt=config%final_time-time
      dt=min(dt,stop_time-time)
      do attempt=1,20
        if (time+dt<=time) error stop 'Tree timestep cannot advance time'
        if (sdc) then
          if(present(evaluate_diffusion)) then
            call advance_sdc(evaluate_stage_diffusion)
          else
            call advance_sdc()
          end if
        else
          call advance_amr_patch_tree_split_3d(tree,eos_species,reactions,transport,state,temperature,root_h,dt,faces, &
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
      if(initial_reactions_pending) write(*,'(a)') 'Initialized MOL reaction history for first accepted step'
      initial_reactions_pending=.false.
      previous_timestep=dt;post_regrid=.false.
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
    if(heat_known) then
    write(unit,'(a)',iostat=status) 'node,level,dx,dy,dz,'//reactive_3d_csv_header(eos_species)//',heat_release_W_m3'
    else
    write(unit,'(a)',iostat=status) 'node,level,dx,dy,dz,'//reactive_3d_csv_header(eos_species)
    end if
    if (status/=0) error stop 'Cannot write tree CSV header'
    allocate(character(len=26*(14+2*size(eos_species))-1+26*merge(1,0,heat_known))::row)
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
            if(sdc) held_cell=internal_energy%fields(node)%values(1,i,j,k)
            q=state%fields(node)%values(:,i,j,k)
            if(allocated(frame)) then
              call convert_rotating_energy_3d(frame,[x(i),y(j),z(k)],q(1),q(5),.false.,ordinary_energy,ok)
              if(.not.ok) error stop 'Invalid rotating tree output state'
              q(5)=ordinary_energy
            endif
            call format_reactive_3d_csv_row(eos_species,[x(i),y(j),z(k)],time,q, &
              temperature%fields(node)%values(1,i,j,k),row(:len(row)-26*merge(1,0,heat_known)),ok,internal_energy_density=held_cell)
            if (.not.ok) error stop 'Invalid tree output state'
            if(heat_known) then
              write(row(len(row)-25:),'(",",es25.16e3)',iostat=status) heat_release%fields(node)%values(1,i,j,k)
              if(status/=0) error stop 'Cannot format stage heat'
            end if
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
      if(sdc) then
        write(unit,'(a)',iostat=status) 'component,initial,final,boundary_change,source_change,energy_correction'
      else
      write(unit,'(a)',iostat=status) 'component,initial,final,boundary_change,source_change'
      endif
    else
      write(unit,'(a)',iostat=status) 'component,initial,final,boundary_change'
    end if
    if (status/=0) error stop 'Cannot write tree budget header'
    do c=1,size(q)
      if (sdc) then
        if(sdc) then
          repair=0
          if(c==5) repair=energy_correction
          write(unit,'(i0,5(",",es25.16e3))',iostat=status) c,initial(c),final(c),budget(c),source_budget(c),repair
        else
        write(unit,'(i0,4(",",es25.16e3))',iostat=status) c,initial(c),final(c),budget(c),source_budget(c)
        endif
      else
        write(unit,'(i0,3(",",es25.16e3))',iostat=status) c,initial(c),final(c),budget(c)
      end if
      if (status/=0) error stop 'Cannot write tree budget row'
    end do
    close(unit,iostat=status)
    if (status/=0) error stop 'Cannot close tree budget CSV'
    if(sdc) write(*,'(a)') 'SDC reactor: '//trim(config%sdc_reactor)
    if(allocated(pressure_work_change)) write(*,'(a,es24.16)') 'Internal energy pressure work: ',pressure_work_change
    write(*,'(a,i0,a,es12.4)') 'Regrids: ',regrids,' Maximum remap error: ',max_remap_error
    write(*,'(a,i0,a,es24.16)') 'Steps: ',step,' Final time: ',time
    write(*,'(a,es12.4,a,es12.4)') 'Maximum reflux: ',max_correction,' Minimum transport theta: ',min_theta
  contains
    subroutine advance_sdc(callback)
      procedure(amr_sdc_diffusion_evaluator), optional :: callback
      if(initial_reactions_pending) then
        call initialize_amr_tree_reaction_source_3d(tree,eos_species,reactions,state,temperature,dt,subcycling, &
          config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance,reaction_source,ok, &
          evaluate_reaction=evaluate_reaction,internal_energy=internal_energy,frame=frame, &
          domain_lower=rotation_lower,root_spacing=root_h)
        if(.not.ok) return
      endif
          if (subcycling) then
          call advance_amr_patch_tree_sdc_subcycle_3d(tree,eos_species,reactions,transport,state,temperature,root_h,dt,faces, &
            config%sdc_iterations,config%chemistry_enabled,[viscosity,conduction,diffusion,barodiffusion], &
            config%external_forcing,config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance, &
            change,source_change,correction,theta,ok,reaction_source=reaction_source, &
            inflow_state=inlet,inflow_temperature=inlet_t,wall_temperature=config%wall_temperature, &
            wall_velocity=config%wall_velocity,outlet_pressure=config%outlet_pressure, &
            inflow_mach=config%inflow_mach,inflow_pressure=mach_pressure,difmag=resolved_reactive_difmag_3d(config), &
            evaluate_diffusion=callback,evaluate_reaction=evaluate_reaction, &
          constant_coefficients=constant_coefficients,heat_release=heat_release, &
          mol_mode=mol,mol_lax_friedrichs=config%mol_lax_friedrichs, &
          internal_energy=internal_energy,energy_correction=energy_correction,energy_policy=energy_policy, &
          frame=frame,domain_lower=rotation_lower,internal_energy_pressure_work_change=step_pressure_work, &
            isothermal_slip_faces=isothermal_slip_faces)
          else
          call advance_amr_patch_tree_sdc_3d(tree,eos_species,reactions,transport,state,temperature,root_h,dt,faces, &
            config%sdc_iterations,config%chemistry_enabled,[viscosity,conduction,diffusion,barodiffusion], &
            config%external_forcing,config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance, &
            change,source_change,correction,theta,ok,reaction_source=reaction_source, &
            inflow_state=inlet,inflow_temperature=inlet_t,wall_temperature=config%wall_temperature, &
            wall_velocity=config%wall_velocity,outlet_pressure=config%outlet_pressure, &
            inflow_mach=config%inflow_mach,inflow_pressure=mach_pressure,difmag=resolved_reactive_difmag_3d(config), &
            evaluate_diffusion=callback,evaluate_reaction=evaluate_reaction, &
          constant_coefficients=constant_coefficients,heat_release=heat_release, &
          mol_mode=mol,mol_lax_friedrichs=config%mol_lax_friedrichs, &
          internal_energy=internal_energy,energy_correction=energy_correction,energy_policy=energy_policy, &
          frame=frame,domain_lower=rotation_lower,internal_energy_pressure_work_change=step_pressure_work, &
            isothermal_slip_faces=isothermal_slip_faces)
          end if
      if(ok) then
        heat_known=.true.
        if(allocated(pressure_work_change)) then
          if(step_pressure_work>0) then
            if(pressure_work_change>huge(1._dp)-step_pressure_work) error stop 'Pressure-work history overflow'
          else if(step_pressure_work<0) then
            if(pressure_work_change< -huge(1._dp)-step_pressure_work) error stop 'Pressure-work history overflow'
          endif
          pressure_work_change=pressure_work_change+step_pressure_work
        endif
      endif
    end subroutine advance_sdc
    subroutine evaluate_stage_diffusion(u,t,depth,fractions,coarse_start,coarse_end,rates,flux,good, &
      internal_energy)
      type(amr_tree_state_3d), optional, intent(in) :: internal_energy
      type(amr_tree_state_3d), intent(in) :: u,t
      type(amr_tree_state_3d), optional, intent(in) :: coarse_start,coarse_end
      integer, intent(in) :: depth
      real(dp), optional, intent(in) :: fractions(:)
      type(amr_tree_state_3d), intent(inout) :: rates
      type(amr_tree_flux_3d), allocatable, intent(inout) :: flux(:)
      logical, intent(out) :: good
      real(dp) :: maximum_diffusion
      call evaluate_diffusion(tree,u,t,root_h,faces,inlet,inlet_t,depth,fractions,coarse_start,coarse_end, &
        rates,flux,maximum_diffusion,good,frame=frame,domain_lower=rotation_lower, &
        internal_energy=internal_energy)
    end subroutine
    subroutine limit_mach_step()
      type(amr_tree_state_3d) :: packed
      real(dp), allocatable :: ghost(:,:,:,:)
      real(dp) :: bound,dmax,factor
      integer :: peer,dims(3),nv,width,np
      character(len=:), allocatable :: donor_limiter
      real(dp), allocatable :: sampled_h(:,:,:)
      nv=size(q)
      np=nv+1+merge(1,0,allocated(internal_energy))
      if(config%reconstruction/='pcm') donor_limiter=config%limiter
      if(allocated(internal_energy)) donor_limiter='linear'
      width=2
      if(config%reconstruction=='characteristic_ppm') width=4
      allocate(packed%fields(size(tree%nodes)))
      do peer=1,size(tree%nodes)
        dims=tree%nodes(peer)%shape
        allocate(packed%fields(peer)%values(np,dims(1),dims(2),dims(3)))
        packed%fields(peer)%values(:nv,:,:,:)=state%fields(peer)%values
        packed%fields(peer)%values(nv+1,:,:,:)=temperature%fields(peer)%values(1,:,:,:)
        if(allocated(internal_energy)) packed%fields(peer)%values(np:np,:,:,:)=internal_energy%fields(peer)%values
      end do
      do peer=1,size(tree%nodes)
        call fill_amr_patch_tree_ghosts_3d(tree,packed,peer,faces,ghost,ok,viscous_wall=.true., &
          inflow_state=inlet_with_temperature(),wall_velocity=config%wall_velocity,coarse_limiter=donor_limiter, &
          species_count=size(eos_species),ghost_width=width,species=eos_species,inflow_mach=config%inflow_mach, &
          inflow_temperature=inlet_t,inflow_pressure=mach_pressure,frame=frame,domain_lower=rotation_lower, &
          root_spacing=root_h,stored_internal_energy=allocated(internal_energy), &
          limit_internal_energy=allocated(internal_energy))
        if(.not.ok) error stop 'Invalid tree Mach inlet donors'
        if(allocated(frame)) then
          call convert_regular_energy_3d(frame,rotation_lower+(real(origin(:,peer),dp)-real(width,dp))*h(:,peer), &
            h(:,peer),ghost,0,.false.,ok)
          if(.not.ok) error stop 'Invalid rotating tree Mach timestep state'
        endif
        if(allocated(internal_energy)) sampled_h=ghost(np,:,:,:)
        dims=shape(ghost(1,:,:,:));factor=1
        if(subcycling) factor=root_h(1)/h(1,peer)
        call compute_reactive_cfl_timestep_3d(eos_species,ghost(:nv,:,:,:),ghost(nv+1,:,:,:), &
          dims(1),dims(2),dims(3),h(1,peer),h(2,peer),h(3,peer),config%cfl,bound,ok,internal_energy_density=sampled_h)
        if(.not.ok) error stop 'Invalid tree Mach hydro timestep'
        dt=min(dt,factor*bound)
        if(config%transport_enabled) then
          call reactive_transport_timestep_3d(eos_species,transport,ghost(:nv,:,:,:),ghost(nv+1,:,:,:), &
            dims(1),dims(2),dims(3),h(1,peer),h(2,peer),h(3,peer),config%transport_cfl, &
            viscosity,conduction,diffusion,bound,dmax,ok, &
          constant_coefficients=constant_coefficients,internal_energy_density=sampled_h)
          if(.not.ok) error stop 'Invalid tree Mach transport timestep'
          dt=min(dt,factor**2*bound)
        end if
      end do
    end subroutine
    function inlet_with_temperature() result(packed_inlet)
      real(dp), allocatable :: packed_inlet(:,:)
      integer :: face,np
      np=size(q)+1+merge(1,0,allocated(internal_energy))
      allocate(packed_inlet(np,6))
      packed_inlet(:size(q),:)=inlet;packed_inlet(size(q)+1,:)=inlet_t
      if(allocated(internal_energy)) then
        packed_inlet(np,:)=0
        do face=1,6
          if(inlet(1,face)<=0) cycle
          packed_inlet(np,face)=inlet(5,face)-.5_dp*sum(inlet(2:4,face)**2)/inlet(1,face)
        enddo
      endif
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
      if(sdc) then
        snapshot%internal_energy=internal_energy;snapshot%energy_correction=energy_correction
        snapshot%pressure_work_change=pressure_work_change
      endif
      snapshot%tree=tree;snapshot%state=state;snapshot%temperature=temperature
      snapshot%initial=initial;snapshot%boundary=budget
      if (sdc) snapshot%source=source_budget
      if (allocated(reaction_source)) snapshot%reaction_source=reaction_source
      if(heat_known) snapshot%heat_release=heat_release
      snapshot%time=time;snapshot%steps=step;snapshot%regrids=regrids
      snapshot%previous_timestep=previous_timestep;snapshot%post_regrid=post_regrid
      snapshot%maximum_reflux=max_correction;snapshot%minimum_theta=min_theta
      snapshot%maximum_remap_error=max_remap_error
      call write_amr_tree_checkpoint_3d(trim(checkpoint_file),eos_species,context,snapshot,ok,indexed=checkpoint_indexed, &
        mol_energy_reset=sdc.and..not.allocated(frame))
      if (.not.ok) error stop 'Cannot write tree checkpoint'
    end subroutine
    subroutine restore_checkpoint()
      integer :: applied,peer,nconserved
      type(amr_tree_state_3d) :: absolute_state
      real(dp) :: original_scale(size(q)),restart_scale(size(q)),inventory(size(q))
      call read_amr_tree_checkpoint_3d(trim(restart_file),eos_species,context, &
        [config%nx,config%ny,config%nz],snapshot,ok,allow_internal_energy=sdc, &
        mol_energy_reset=sdc.and..not.allocated(frame))
      if (.not.ok) error stop 'Invalid tree restart checkpoint/context'
      if(sdc.neqv.allocated(snapshot%internal_energy%fields)) error stop 'Tree retained energy presence mismatch'
      if (sdc.neqv.allocated(snapshot%source)) error stop 'Tree restart source ledger mismatch'
      if (allocated(snapshot%reaction_source%fields).neqv.(sdc.and.config%chemistry_enabled)) &
        error stop 'Tree restart reaction source mismatch'
      if (sdc) source_budget=snapshot%source
      if (allocated(reaction_source)) reaction_source=snapshot%reaction_source
      if(.not.allocated(snapshot%previous_timestep)) error stop 'Tree restart missing timestep history'
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
      if(sdc) inventory(5)=inventory(5)-snapshot%energy_correction
      nconserved=size(q)
      if (config%chemistry_enabled) nconserved=5
      if (any(abs(inventory(:nconserved)-snapshot%initial(:nconserved)- &
          snapshot%boundary(:nconserved)-source_budget(:nconserved))> &
          2e-10_dp*max(tiny(1.0_dp),original_scale(:nconserved),restart_scale(:nconserved), &
          abs(snapshot%boundary(:nconserved)),abs(source_budget(:nconserved))))) error stop 'Tree restart boundary ledger mismatch'
      tree=snapshot%tree;state=snapshot%state;temperature=snapshot%temperature
      if(sdc) then
        heat_known=allocated(snapshot%heat_release%fields)
        if(heat_known) then
          heat_release=snapshot%heat_release
        else
          heat_release=snapshot%temperature
          do peer=1,size(snapshot%tree%nodes)
            heat_release%fields(peer)%values=0
          end do
        end if
      end if
      if(sdc) then
        internal_energy=snapshot%internal_energy;energy_correction=snapshot%energy_correction
        if(allocated(snapshot%pressure_work_change)) then
          pressure_work_change=snapshot%pressure_work_change
        else if(mol) then
          ! Previous MOL snapshots have the same exactly zero pdivu source.
          pressure_work_change=0
        else
          error stop 'Godunov restart missing pressure-work history'
        endif
      endif
      initial=snapshot%initial;budget=snapshot%boundary
      time=snapshot%time;step=snapshot%steps;regrids=snapshot%regrids
      previous_timestep=snapshot%previous_timestep;post_regrid=snapshot%post_regrid
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
      type(amr_tree_state_3d), allocatable :: remapped_reaction,remapped_heat,remapped_internal
      real(dp) :: held_before(1),held_after(1),held_scale(1)
      integer(int64), allocatable :: o(:,:),d(:,:)
      integer, allocatable :: r(:)
      real(dp), allocatable :: spacing(:,:)
      real(dp) :: before(size(q)),after(size(q)),scale(size(q)),error
      integer :: peer
      logical :: valid
      if (allocated(reaction_source)) allocate(remapped_reaction)
      if(allocated(heat_release)) allocate(remapped_heat)
      if(sdc) allocate(remapped_internal)
      call remap_amr_patch_tree_reactive_3d(tree,state,temperature,next_tree,eos_species,candidate,temps,valid, &
        old_reaction_source=reaction_source,reaction_source=remapped_reaction, &
        old_heat_release=heat_release,heat_release=remapped_heat, &
        old_internal_energy=internal_energy,internal_energy=remapped_internal)
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
      if(sdc) then
        call composite_integrals_amr_patch_tree_3d(tree,internal_energy,root_h,held_before,valid)
        if(.not.valid) error stop 'Invalid retained pre-remap inventory'
        call composite_integrals_amr_patch_tree_3d(next_tree,remapped_internal,root_h,held_after,valid)
        if(.not.valid) error stop 'Invalid retained post-remap inventory'
        absolute_state=internal_energy
        do peer=1,size(tree%nodes)
          absolute_state%fields(peer)%values=abs(internal_energy%fields(peer)%values)
        enddo
        call composite_integrals_amr_patch_tree_3d(tree,absolute_state,root_h,held_scale,valid)
        if(.not.valid) error stop 'Invalid retained remap scale'
        if(any(abs(held_after-held_before)>5e-12_dp*max(tiny(1._dp),held_scale))) &
          error stop 'Scheduled remap lost retained energy'
        internal_energy=remapped_internal
      endif
      tree=next_tree;state=candidate;temperature=temps
      if (allocated(reaction_source)) reaction_source=remapped_reaction
      if(allocated(heat_release)) heat_release=remapped_heat
      origin=o;domain=d;ratios=r;h=spacing;count=size(tree%nodes)-1
      max_remap_error=max(max_remap_error,error)
      write(*,'(a,es24.16,a,i0,a,i0)') 'Regrid time: ',time,' Nodes: ',count+1,' Levels: ',tree%level_count()
      post_regrid=.true.
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
