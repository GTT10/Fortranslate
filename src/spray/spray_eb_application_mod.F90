module spray_eb_application_mod
  use precision_mod, only: dp
  use, intrinsic :: iso_c_binding, only: c_char,c_size_t,c_ptr,c_associated,c_null_char
  use simulation_config_selected_reactor_mod, only: selected_reactor_paths_alias_lexically
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config,build_configured_eb_geometry_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d,write_reactive_eb_3d_csv, &
    reactive_eb_integrals_3d
  use reactive_eb_cfl_3d_mod, only: compute_reactive_eb_cfl_timestep_3d
  use eb_reactive_transport_3d_mod, only: reactive_eb_transport_timestep_3d
  use spray_parcel_mod
  use spray_mixture_mod, only: liquid_mixture_properties
  use spray_eb_reactive_mod, only: advance_eb_spray_reactive_allocated
  use spray_eb_impact_mod, only: spray_eb_wall_options
  use spray_wall_impact_mod, only: spray_wall_film
  use spray_eb_film_step_mod, only: advance_eb_film_coupling
  use spray_eb_step_mod, only: advance_eb_mixture_coupling_allocated
  use spray_eb_field_mod, only: eb_spray_stencil
  use spray_checkpoint_mod, only: write_spray_checkpoint,read_spray_checkpoint,spray_checkpoint_context_capacity
  use amr_reactive_3d_checkpoint_mod, only: write_selected_checkpoint_context
  implicit none
  private
  public :: run_eb_spray_application
  interface
    function system_getcwd(buffer,length) bind(C,name='getcwd') result(pointer)
      import c_char,c_size_t,c_ptr
      character(kind=c_char), intent(out) :: buffer(*)
      integer(c_size_t), value :: length
      type(c_ptr) :: pointer
    end function
  end interface
contains
  subroutine run_eb_spray_application(input_path,spray_path,config,species,reactions,transport,mole_fractions,integrator)
    character(len=*), intent(in) :: input_path,spray_path
    type(reactive_eb_3d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), optional, intent(in) :: mole_fractions(:)
    character(len=*), optional, intent(in) :: integrator
    integer :: nc=1,number=1,injection_events=1
    real(dp) :: injection_times(1024)=0,injection_scales(1024)=1
    character(len=24) :: vapor_names(16)='H2O'
    real(dp) :: densities(16)=1000,heat_capacities(16)=4184,reference_temperatures(16)=373.15_dp
    real(dp) :: reference_pressures(16)=101325,latent_heats(16)=2.257e6_dp,fractions(16)=0
    real(dp) :: total_mass=1e-10_dp,center(3)=0,axis(3)=[1._dp,0._dp,0._dp],speed=1,half_angle=0
    real(dp) :: diameter=1e-5_dp,liquid_temperature=350,maximum_dt=1e-8_dp
    type(spray_options) :: options
    type(spray_eb_wall_options) :: walls
    logical :: wall_impacts=.false.
    type(spray_wall_film), allocatable :: films(:)
    real(dp), allocatable :: film_mass(:,:)
    namelist /eb_spray/ nc,number,vapor_names,densities,heat_capacities,reference_temperatures, &
      reference_pressures,latent_heats,fractions,total_mass,center,axis,speed,half_angle,diameter, &
      liquid_temperature,maximum_dt,options,injection_events,injection_times,injection_scales,wall_impacts,walls
    type(eb_geometry_3d) :: g
    type(reactive_eb_3d_config) :: canonical
    type(spray_liquid), allocatable :: liquids(:)
    type(spray_liquid) :: effective
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: phi(:,:,:),state(:,:,:,:),temperature(:,:,:),masses(:,:)
    real(dp), allocatable :: boundary(:),wall(:),chemical(:),sb(:),sw(:),sc(:),injected(:),initial(:)
    real(dp) :: time,dt,hydro_dt,transport_dt,diffusivity,theta,rho,cp,e
    integer :: unit,status,close_status,nv,n,p,step,attempt,particle_steps,next_event
    integer(int64) :: last_id,step_split_created
    integer(int64), allocatable :: split_created_total
    logical :: ok,model_rejected
    integer :: last_checkpoint_step
    character(len=:), allocatable :: context,checkpoint_path
    character(len=1024) :: message
    character(len=24) :: component_label
    character(len=:), allocatable :: header
    ! Configuration and namelist objects are local to one invocation.
    nc=1;number=1;fractions=0;fractions(1)=1
    vapor_names='H2O';densities=1000;heat_capacities=4184;reference_temperatures=373.15_dp
    reference_pressures=101325;latent_heats=2.257e6_dp;total_mass=1e-10_dp;center=0
    axis=[1._dp,0._dp,0._dp];speed=1;half_angle=0;diameter=1e-5_dp;liquid_temperature=350
    maximum_dt=1e-8_dp;options=spray_options();walls=spray_eb_wall_options();wall_impacts=.false.
    injection_events=1;injection_times=0;injection_scales=1;next_event=1;last_id=0
    call check_output(config%output_file)
    call check_output(trim(config%output_file)//'.particles.csv')
    call check_output(trim(config%output_file)//'.films.csv')
    if (config%geometry_type/='level_set') error stop 'EB spray requires nodal level-set geometry'
    if(len_trim(config%restart_file)>0) call check_checkpoint_path(config%restart_file)
    if(config%checkpoint_interval_steps>0) call check_checkpoint_path(config%checkpoint_file)
    open(newunit=unit,file=spray_path,status='old',action='read',iostat=status)
    if(status/=0) error stop 'Cannot open EB spray input'
    read(unit,nml=eb_spray,iostat=status)
    close(unit,iostat=close_status)
    if(status/=0.or.close_status/=0) error stop 'Invalid EB spray namelist'
    if(nc<1.or.nc>16.or.number<1.or.number>125000) error stop 'Invalid EB spray population or components'
    if(.not.ieee_is_finite(maximum_dt)) error stop 'Invalid EB spray maximum_dt'
    if(maximum_dt<=0.or.maximum_dt>1e6_dp) error stop 'Invalid EB spray maximum_dt'
    if(.not.all(ieee_is_finite(fractions)).or.any(fractions<0)) error stop 'Invalid liquid fractions'
    if(abs(sum(fractions(:nc))-1)>1e-12_dp.or.any(fractions(nc+1:)/=0)) error stop 'Liquid fractions must sum to one'
    if(.not.valid_spray_options(options).or.options%khrt_children) error stop 'Unsupported EB spray options'
    if(options%maximum_multiplicity>0) then
      allocate(split_created_total);split_created_total=0
    end if
    if(injection_events<1.or.injection_events>1024) error stop 'Invalid injection_events'
    if(.not.all(ieee_is_finite(injection_times)).or..not.all(ieee_is_finite(injection_scales))) &
      error stop 'Nonfinite injection schedule'
    if(any(injection_times<0).or.any(injection_times>1e6_dp).or.any(injection_scales<0).or. &
       any(injection_scales>1e6_dp)) error stop 'Invalid injection times or scales'
    do n=2,injection_events
      if(injection_times(n)<injection_times(n-1)) error stop 'Injection times must be nondecreasing'
    end do
    call build_configured_eb_geometry_3d(config,g,ok,phi)
    if(.not.ok) error stop 'Cannot build EB spray geometry'
    allocate(liquids(nc))
    do n=1,nc
      call initialize_spray_liquid(species,vapor_names(n),densities(n),heat_capacities(n),reference_temperatures(n), &
        reference_pressures(n),latent_heats(n),liquids(n),ok)
      if(.not.ok) error stop 'Invalid EB spray liquid table'
    end do
    call liquid_mixture_properties(liquids,fractions(:nc),liquid_temperature,rho,cp,e,ok)
    if(.not.ok) error stop 'Invalid injected liquid mixture'
    effective=liquids(1);effective%density=rho;effective%heat_capacity=cp
    effective%reference_temperature=liquid_temperature;effective%reference_energy=e
    allocate(drops(0),masses(nc,0),films(0),film_mass(nc,0))
    nv=reactive_nvar(size(species))
    allocate(state(nv,g%nx,g%ny,g%nz),temperature(g%nx,g%ny,g%nz))
    allocate(boundary(nv),wall(nv),chemical(nv),sb(nv),sw(nv),sc(nv),injected(nv),initial(nv))
    call initialize_reactive_eb_density_sheet_3d(species,config,g,state,temperature,ok,mole_fractions)
    if(.not.ok) error stop 'Invalid EB initial gas'
    if (walls%evolve.and..not.wall_impacts) error stop 'Film evolution requires wall_impacts'
    if (wall_impacts) then
      if (options%surface_tension<1e-8_dp.or.liquid_viscosity_at_temperature(options,liquid_temperature)<=0) &
        error stop 'EB wall impacts require positive surface tension and viscosity'
      call advance_eb_mixture_coupling_allocated(species,transport,liquids,options,g,phi,state,temperature, &
        drops,masses,0._dp,particle_steps,sw(:5),ok,walls,films,film_mass,last_id)
      if (.not.ok) error stop 'Invalid EB impact configuration'
    end if
    call reactive_eb_integrals_3d(state,g,initial,ok)
    if(.not.ok) error stop 'Invalid initial EB gas inventory'
    time=0;step=0;boundary=0;wall=0;chemical=0;injected=0;last_checkpoint_step=0
    if(config%checkpoint_interval_steps>0.or.len_trim(config%restart_file)>0) call build_context()
    if(len_trim(config%restart_file)>0) then
      call read_spray_checkpoint(config%restart_file,context,species,state,temperature,drops,time,step,next_event, &
        initial,injected,ok,message,boundary_totals=boundary,source_totals=wall,chemistry_totals=chemical, &
        breakup_state=.true.,split_state=.true.,khrt_state=.true.,material_state=.true.,shedding_state=.true., &
        component_mass=masses,last_id=last_id,film_state=.true.,wall_films=films,film_component_mass=film_mass, &
        retired_parent_state=wall_impacts,split_created_total=split_created_total)
      if(.not.ok) then
        write(*,'(a)') trim(message)
        error stop 'EB spray restart rejected'
      end if
      call validate_restart()
      last_checkpoint_step=step
    end if
    do
      do while(next_event<=injection_events)
        if(injection_times(next_event)>time) exit
        call inject_event()
        next_event=next_event+1
      end do
      if(config%checkpoint_interval_steps>0.and.step>last_checkpoint_step) then
        if(modulo(step,config%checkpoint_interval_steps)==0) then
          write(component_label,'(i0)') step
          checkpoint_path=trim(config%checkpoint_file)//'.step'//trim(component_label)
          call check_checkpoint_path(checkpoint_path)
          call write_spray_checkpoint(checkpoint_path,context,species,state,temperature,drops,time,step,next_event, &
            initial,injected,ok,message,boundary_totals=boundary,source_totals=wall,chemistry_totals=chemical, &
            breakup_state=.true.,split_state=.true.,khrt_state=.true.,material_state=.true.,shedding_state=.true., &
            component_mass=masses,last_id=last_id,film_state=.true.,wall_films=films,film_component_mass=film_mass, &
        retired_parent_state=wall_impacts,split_created_total=split_created_total)
          if(.not.ok) then
            write(*,'(a)') trim(message)
            error stop 'EB spray checkpoint rejected'
          end if
          last_checkpoint_step=step
          if(config%stop_after_checkpoint) exit
        end if
      end if
      if(time>=config%final_time) exit
      if(step>=config%maximum_steps) error stop 'EB spray maximum_steps exhausted'
      call compute_reactive_eb_cfl_timestep_3d(species,state,temperature,g,config%cfl,hydro_dt,ok)
      if(.not.ok) error stop 'Invalid EB gas CFL bound'
      dt=min(hydro_dt,maximum_dt,config%final_time-time)
      if(next_event<=injection_events) dt=min(dt,injection_times(next_event)-time)
      if(config%transport_enabled) then
        call reactive_eb_transport_timestep_3d(species,transport,state,temperature,g,config%transport_cfl, &
          config%viscosity_enabled,config%thermal_conduction_enabled,config%species_diffusion_enabled, &
          transport_dt,diffusivity,ok)
        if(.not.ok) error stop 'Invalid EB transport bound'
        dt=min(dt,transport_dt)
      end if
      do attempt=1,30
        ok=.false.
        if(dt<=0.or.time+dt<=time) exit
        theta=1;sb=0;sw=0;sc=0
        if (wall_impacts) then
        call advance_eb_spray_reactive_allocated(species,reactions,transport,liquids,options,g,phi,state,temperature, &
          drops,masses,dt,config%riemann_solver,config%redistribution,config%chemistry_enabled,config%transport_enabled, &
          config%viscosity_enabled,config%thermal_conduction_enabled,config%species_diffusion_enabled, &
          config%barodiffusion_enabled,config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance, &
          particle_steps,theta,sb,sw,sc,ok,config%state_redist_target_volume_fraction,integrator, &
          walls,films,film_mass,last_id,model_rejected,step_split_created)
          if(model_rejected) error stop 'EB impact correlation outside admissible mass/velocity range'
        else
        call advance_eb_spray_reactive_allocated(species,reactions,transport,liquids,options,g,phi,state,temperature, &
          drops,masses,dt,config%riemann_solver,config%redistribution,config%chemistry_enabled,config%transport_enabled, &
          config%viscosity_enabled,config%thermal_conduction_enabled,config%species_diffusion_enabled, &
          config%barodiffusion_enabled,config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance, &
          particle_steps,theta,sb,sw,sc,ok,config%state_redist_target_volume_fraction,integrator, &
          last_id=last_id,split_created=step_split_created)
        end if
        if(ok) exit
        dt=.5_dp*dt
      end do
      if(.not.ok) error stop 'EB spray coupled step rejected'
      if(allocated(split_created_total)) split_created_total=split_created_total+step_split_created
      time=time+dt;step=step+1;boundary=boundary+sb;wall=wall+sw;chemical=chemical+sc
    end do
    call write_reactive_eb_3d_csv(config%output_file,species,g,state,temperature,time,ok,message)
    if(.not.ok) error stop 'Cannot write EB spray gas CSV'
    open(newunit=unit,file=trim(config%output_file)//'.particles.csv',status='replace',action='write',iostat=status)
    if(status/=0) error stop 'Cannot write EB spray particle CSV'
    header='id,x,y,z,u,v,w,mass,temperature,multiplicity'
    do n=1,nc
      write(component_label,'(i0)') n
      header=header//',component_mass_'//trim(component_label)
    end do
    write(unit,'(a)',iostat=status) header
    do p=1,size(drops)
      if(status/=0) exit
      write(unit,'(i0,*(",",es25.16,:))',iostat=status) drops(p)%id,drops(p)%position,drops(p)%velocity, &
        drops(p)%mass,drops(p)%temperature,drops(p)%multiplicity,masses(:,p)
    end do
    close(unit,iostat=close_status)
    if(status/=0.or.close_status/=0) error stop 'EB particle CSV write failed'
    open(newunit=unit,file=trim(config%output_file)//'.films.csv',status='replace',action='write',iostat=status)
    if(status/=0) error stop 'Cannot write EB film CSV'
    header='material_id,x,y,z,nx,ny,nz,mass,temperature,diameter,height'
    do n=1,nc
      write(component_label,'(i0)') n
      header=header//',component_mass_'//trim(component_label)
    end do
    write(unit,'(a)',iostat=status) header
    do p=1,size(films)
      if(status/=0) exit
      write(unit,'(i0,*(",",es25.16,:))',iostat=status) films(p)%material_id,films(p)%position,films(p)%normal, &
        films(p)%mass,films(p)%temperature,films(p)%diameter,films(p)%height,film_mass(:,p)
    end do
    close(unit,iostat=close_status)
    if(status/=0.or.close_status/=0) error stop 'EB film CSV write failed'
    write(*,'(a,i0,a,es25.16)') 'EB spray steps: ',step,' time: ',time
    write(*,'(a,i0)') 'Injection events processed: ',next_event-1
    write(*,'(a,*(es25.16,1x))') 'Injected inventory: ',injected
    write(*,'(a,*(es25.16,1x))') 'Boundary inventory: ',boundary
    write(*,'(a,*(es25.16,1x))') 'Wall inventory: ',wall
    write(*,'(a,*(es25.16,1x))') 'Chemistry inventory: ',chemical
  contains
    subroutine build_context()
      integer :: scratch,ios,j,k,length
      character(len=spray_checkpoint_context_capacity) :: line
      character(len=32) :: method
      real(dp) :: composition(size(species))
      canonical=config
      canonical%final_time=0;canonical%maximum_steps=0;canonical%output_file=''
      canonical%checkpoint_interval_steps=0;canonical%checkpoint_file='';canonical%restart_file=''
      canonical%stop_after_checkpoint=.false.;canonical%level_set_file=''
      method='default';if(present(integrator)) method=integrator
      composition=0
      if(present(mole_fractions)) composition=mole_fractions
      open(newunit=scratch,status='scratch',form='formatted',delim='quote',iostat=ios)
      if(ios/=0) error stop 'Cannot prepare EB checkpoint context'
      write(scratch,'(a)',iostat=ios) 'EB_SPRAY_CONTEXT_V2'
      if(ios==0) write(scratch,*,iostat=ios) canonical
      if(ios==0) write(scratch,nml=eb_spray,iostat=ios)
      do j=1,size(species)
        if(ios==0) write(scratch,*,iostat=ios) species(j)
        if(ios==0) write(scratch,*,iostat=ios) transport(j)
      end do
      if(ios==0) call write_selected_checkpoint_context(scratch,reactions,'EB_SPRAY_V1',method,composition, &
        config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance,ios)
      do k=0,g%nz
      do j=0,g%ny
        if(ios==0) write(scratch,'(*(es26.17e3,1x))',iostat=ios) phi(:,j,k)
      end do
      end do
      if(ios/=0) error stop 'Cannot serialize EB checkpoint context'
      rewind(scratch,iostat=ios)
      if(ios/=0) error stop 'Cannot rewind EB checkpoint context'
      context=''
      do
        read(scratch,'(a)',iostat=ios) line
        if(ios<0) exit
        if(ios/=0.or.len_trim(line)==len(line)) error stop 'EB checkpoint context line exceeds capacity'
        length=len(context)+len_trim(line)+1
        if(length>spray_checkpoint_context_capacity) error stop 'EB checkpoint context exceeds capacity'
        context=context//trim(line)//'|'
      end do
      close(scratch)
    end subroutine
    subroutine validate_restart()
      integer :: q,expected_count,cells(3,8)
      real(dp) :: weights(8),volumes(8)
      logical :: valid
      if(time>config%final_time.or.step>config%maximum_steps) error stop 'Restart exceeds continuation limits'
      if(next_event<1.or.next_event>injection_events+1) error stop 'Invalid restart event index'
      if(next_event<=injection_events) then
        if(injection_times(next_event)<=time) error stop 'Unprocessed restart injection event'
      end if
      if(next_event>1) then
        if(injection_times(next_event-1)>time) error stop 'Restart injection ahead of time'
      end if
      expected_count=count(injection_scales(:next_event-1)>0)*number
      step_split_created=0
      if(allocated(split_created_total)) step_split_created=split_created_total
      if(last_id<int(expected_count,int64)) error stop 'Invalid restart injection watermark'
      if(step_split_created<0.or.step_split_created>last_id-int(expected_count,int64)) &
        error stop 'Invalid restart split count'
      if (wall_impacts) then
        if(last_id<int(expected_count,int64)) error stop 'Invalid restart impact watermark'
        if(modulo(last_id-int(expected_count,int64)-step_split_created,4_int64)/=0) error stop 'Invalid restart splash IDs'
      else
        if(size(films)/=0) error stop 'Unexpected restart films'
        if(int(size(drops),int64)/=last_id.or.last_id-int(expected_count,int64)/=step_split_created) &
          error stop 'Invalid restart injection IDs'
      end if
      if(size(drops)>min(125000,options%maximum_parcels)) error stop 'Restart exceeds population capacity'
      if(any(shape(masses)/=[nc,size(drops)])) error stop 'Restart component shape mismatch'
      cells=0;weights=0;volumes=0
      do q=1,size(drops)
        if (wall_impacts) then
          if(drops(q)%id>last_id) error stop 'Restart ID exceeds watermark'
          if(q>1) then
            if(drops(q)%id<=drops(q-1)%id) error stop 'Restart IDs not ordered'
          end if
        else
          if(drops(q)%id/=int(q,int64)) error stop 'Restart parcel ID sequence mismatch'
        end if
        call eb_spray_stencil(g,phi,drops(q)%position,cells,weights,volumes,valid)
        if(.not.valid) error stop 'Restart parcel outside fluid'
      end do
      if (wall_impacts) then
        call advance_eb_film_coupling(species,transport,liquids,options,walls,g,phi,state,temperature, &
          films,film_mass,0._dp,particle_steps,sw(:5),valid)
        if(.not.valid) error stop 'Invalid restart EB films'
      end if
    end subroutine
    subroutine check_checkpoint_path(path)
      character(len=*), intent(in) :: path
      call check_output(path)
      if(selected_reactor_paths_alias_lexically(absolute_path(path),absolute_path(config%output_file)).or. &
         selected_reactor_paths_alias_lexically(absolute_path(path), &
           absolute_path(trim(config%output_file)//'.particles.csv')).or. &
         selected_reactor_paths_alias_lexically(absolute_path(path), &
           absolute_path(trim(config%output_file)//'.films.csv'))) error stop 'Checkpoint aliases CSV output'
    end subroutine
    subroutine inject_event()
      type(spray_parcel), allocatable :: added(:),candidate(:)
      real(dp), allocatable :: candidate_mass(:,:)
      real(dp) :: change(nv),totals(5),weights(8),volumes(8)
      integer :: cells(3,8),q,j,old_count,ios
      logical :: valid
      if(injection_scales(next_event)==0) return
      old_count=size(drops)
      if(old_count>min(125000,options%maximum_parcels)-number) error stop 'EB injection capacity exceeded'
      call inject_cone_parcels(effective,total_mass*injection_scales(next_event),number,last_id+1_int64, &
        center,axis,speed,half_angle,diameter,liquid_temperature,added,valid)
      if(.not.valid) error stop 'Invalid scheduled EB injection'
      change=0;cells=0;weights=0;volumes=0
      allocate(candidate(old_count+number),candidate_mass(nc,old_count+number),stat=ios)
      if(ios/=0) error stop 'Cannot allocate scheduled injection'
      candidate(:old_count)=drops;candidate_mass(:,:old_count)=masses
      do q=1,number
        call eb_spray_stencil(g,phi,added(q)%position,cells,weights,volumes,valid)
        if(.not.valid) error stop 'EB injection point is not in fluid'
        candidate(old_count+q)=added(q)
        candidate_mass(:,old_count+q)=added(q)%mass*fractions(:nc)
        totals=parcel_totals(added(q),effective);change(1:5)=change(1:5)+totals
        do j=1,nc
          change(5+liquids(j)%vapor_index)=change(5+liquids(j)%vapor_index)+ &
            candidate_mass(j,old_count+q)*added(q)%multiplicity
        end do
      end do
      if(.not.all(ieee_is_finite(change))) error stop 'Invalid injection inventory'
      call move_alloc(candidate,drops);call move_alloc(candidate_mass,masses)
      last_id=last_id+int(number,int64);injected=injected+change
    end subroutine
    function absolute_path(path) result(value)
      character(len=*), intent(in) :: path
      character(len=:), allocatable :: value
      character(kind=c_char) :: buffer(4096)
      type(c_ptr) :: pointer
      integer :: index
      if (len_trim(path)==0) error stop 'Empty EB path'
      if (path(1:1)=='/') then
        value=trim(path)
      else
        pointer=system_getcwd(buffer,int(size(buffer),c_size_t))
        if (.not.c_associated(pointer)) error stop 'Cannot resolve EB working directory'
        value=''
        do index=1,size(buffer)
          if (buffer(index)==c_null_char) exit
          value=value//buffer(index)
        end do
        value=value//'/'//trim(path)
      end if
      if (len(value)>1024) error stop 'EB path too long'
    end function
    subroutine check_output(path)
      character(len=*), intent(in) :: path
      if(selected_reactor_paths_alias_lexically(absolute_path(path),absolute_path(input_path)).or. &
         selected_reactor_paths_alias_lexically(absolute_path(path),absolute_path(spray_path)).or. &
         selected_reactor_paths_alias_lexically(absolute_path(path),absolute_path(config%level_set_file))) &
        error stop 'EB output aliases input'
    end subroutine
  end subroutine
end module
