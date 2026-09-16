module spray_mpi_application_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: compute_reactive_cfl_timestep_3d
  use reactive_transport_3d_mod, only: reactive_transport_timestep_3d
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d,valid_wall_velocity_3d
  use spray_wall_impact_mod, only: spray_wall_film
  use spray_parcel_mod
  use spray_partition_mod, only: spray_cartesian_block
  use spray_coupling_3d_mod, only: mixture_system_totals
  use spray_mpi_checkpoint_mod, only: mpi_spray_checkpoint
  use spray_checkpoint_mod, only: spray_checkpoint_context_capacity
  use spray_injection_profile_mod, only: integrate_spray_injection_profile
  use spray_mixture_mod, only: mixture_evolution_id
  use spray_mpi_mod
  implicit none
  private
  public :: run_mpi_spray_application
contains
  subroutine run_mpi_spray_application(input_path,species,reactions,transport,mechanism_hash)
  character(len=*), intent(in) :: input_path,mechanism_hash
  type(nasa7_species), intent(in) :: species(:)
  type(elementary_reaction), intent(in) :: reactions(:)
  type(gas_transport_species), intent(in) :: transport(:)
  type(spray_liquid) :: liquids(2)
  type(spray_options) :: options
  character(len=8) :: breakup_model='none'
  logical :: khrt_children=.false.
  real(dp) :: khrt_reference_diameter_cubed=0
  real(dp) :: surface_tension=0,liquid_viscosity=0,khrt_b0=.61_dp,khrt_b1=7,khrt_c3=1
  type(spray_size_options) :: size_options,port_size(16)
  type(spray_nozzle_model) :: nozzle,port_nozzle(16)
  real(dp) :: nozzle_radius=0,port_nozzle_radius(16)=0
  type(spray_parcel), allocatable :: parcels(:)
  type(spray_wall_film), allocatable :: checkpoint_films(:)
  real(dp), allocatable :: checkpoint_film_mass(:,:)
  character(len=32) :: boundary_faces(6)='periodic'
  character(len=32), allocatable :: active_faces(:)
  real(dp) :: wall_temperature(6)=600,wall_velocity(3,6)=0
  real(dp) :: boundary_total(5+size(species)),boundary_step(5+size(species))
  real(dp) :: inlet_state(5+size(species),6),inlet_temperature(6)
  logical :: physical_boundaries,periodic_axes(3)
  real(dp), allocatable :: state(:,:,:,:),temperature(:,:,:),masses(:,:)
  integer :: cells(3)=[8,4,4],number=4,max_steps=10000,injection_events=1
  real(dp) :: injection_start=0,injection_interval=1e-8_dp
  integer :: injection_profile_points=0
  real(dp) :: injection_profile_times(64)=0,injection_profile_rates(64)=0,profile_mass
  character(len=4096) :: profile_context
  integer :: port_profile_points(16)=0
  real(dp) :: port_profile_times(64,16)=0,port_profile_rates(64,16)=0
  integer :: port_count=0,port_events(16)=0,port_number(16)=4,schedule_count
  real(dp) :: port_start(16)=0,port_interval(16)=1e-8_dp,port_mass(16)=1e-10_dp
  real(dp) :: port_center(3,16),port_axis(3,16),port_fractions(2,16)=.5_dp
  real(dp) :: port_speed(16)=10,port_half_angle(16)=.1_dp,port_diameter(16)=2e-5_dp
  real(dp) :: port_temperature(16)=330
  real(dp), allocatable :: event_times(:),event_ends(:)
  integer, allocatable :: event_ports(:)
  real(dp) :: spacing(3)=[.001_dp,.001_dp,.001_dp],lower(3)=0,gas_density=1,gas_temperature=800
  real(dp) :: gas_velocity(3)=0
  real(dp), allocatable :: gas_fractions(:)
  real(dp) :: total_mass=1e-10_dp,center(3)=[.004_dp,.002_dp,.002_dp],axis(3)=[1.0_dp,0.0_dp,0.0_dp]
  real(dp) :: speed=10,half_angle=.1_dp,diameter=2e-5_dp,drop_temperature=330,fractions(2)=[.5_dp,.5_dp]
  real(dp) :: liquid_density(2)=[1000.0_dp,800.0_dp],liquid_cp(2)=[4184.0_dp,2000.0_dp]
  real(dp) :: boiling_temperature(2)=[373.15_dp,350.0_dp],latent_heat(2)=[2.257e6_dp,1e6_dp]
  real(dp) :: final_time=2e-8_dp,dt_max=1e-8_dp,cfl=.3_dp,rtol=1e-7_dp,atol=1e-12_dp
  logical :: transport_flags(4)=.true.
  real(dp) :: les_cs=0,les_prt=.7_dp,les_sct=.7_dp
  character(len=24) :: liquid_species(2)=[character(len=24)::'H2O','H2O2'],chemistry='implicit'
  character(len=16) :: solver='hllc',limiter='mc'
  character(len=512) :: output_prefix='mpi_spray',restart_path='',checkpoint_path=''
  character(len=spray_checkpoint_context_capacity) :: checkpoint_context
  character(len=131072) :: input_record
  character(len=6) :: rank_text
  integer :: rank,nranks,ierr,unit,status,history_unit,first(3),extent(3),i,j,k,p,n,step,attempt,substeps,event
  integer(int64) :: last_id=0,local_count,global_count
  real(dp) :: primitive(5+size(species)),gas_cell(5+size(species)),sound,initial_t,time,dt,local_dt,transport_dt,max_diffusivity
  real(dp) :: injected(5+size(species)),local_totals(5+size(species)),totals(5+size(species)), &
    initial_totals(5+size(species)),local_mass
  logical :: ok
  namelist /mpi_spray/ cells,spacing,lower,gas_density,gas_temperature,gas_velocity,gas_fractions, &
    number,total_mass,center,axis,speed,half_angle,diameter,drop_temperature,fractions, &
    liquid_density,liquid_cp,boiling_temperature,latent_heat,liquid_species,final_time,dt_max,cfl, &
    rtol,atol,max_steps,transport_flags,chemistry,solver,limiter,output_prefix,restart_path,checkpoint_path, &
    injection_events,injection_start,injection_interval,injection_profile_points, &
    injection_profile_times,injection_profile_rates,port_count,port_events,port_number,port_start,port_interval, &
    port_mass,port_center,port_axis,port_fractions,port_speed,port_half_angle,port_diameter,port_temperature, &
    port_profile_points,port_profile_times,port_profile_rates,nozzle_radius,port_nozzle_radius, &
    size_options,port_size,nozzle,port_nozzle,breakup_model,surface_tension,liquid_viscosity,khrt_b0,khrt_b1,khrt_c3, &
    khrt_children,khrt_reference_diameter_cubed,les_cs,les_prt,les_sct, &
    boundary_faces,wall_temperature,wall_velocity
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  do p=1,16
    port_center(:,p)=center;port_axis(:,p)=axis
  end do
  allocate(gas_fractions(size(species)));gas_fractions=0
  do n=1,size(species)
    select case(trim(species(n)%name))
    case('H2');gas_fractions(n)=.03_dp
    case('O2');gas_fractions(n)=.23_dp
    case('N2');gas_fractions(n)=.74_dp
    end select
  end do
  if (rank==0) then
    open(newunit=unit,file=trim(input_path),status='old',action='read',iostat=status)
    call require(status==0,'Cannot open input')
    read(unit,nml=mpi_spray,iostat=status)
    call require(status==0,'Invalid mpi_spray namelist')
    close(unit)
    write(input_record,nml=mpi_spray,iostat=status)
    call require(status==0,'Input serialization failed')
  end if
  call MPI_Bcast(input_record,len(input_record),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call require(ierr==MPI_SUCCESS,'Input broadcast failed')
  read(input_record,nml=mpi_spray,iostat=status)
  call require(status==0,'Input deserialization failed')
  call require(all(ieee_is_finite([spacing,lower,gas_density,gas_temperature,gas_velocity,gas_fractions, &
    total_mass,center,axis,speed,half_angle,diameter,drop_temperature,fractions,liquid_density,liquid_cp, &
    boiling_temperature,latent_heat,final_time,dt_max,cfl,rtol,atol,injection_start,injection_interval, &
    injection_profile_times,injection_profile_rates])),'Nonfinite input')
  call require(all(spacing>=1e-12_dp).and.all(spacing<=1e9_dp).and.maxval(abs(lower))<=1e9_dp,'Invalid geometry')
  call require(final_time>=0.and.final_time<=1e6_dp.and.dt_max>0.and.dt_max<=1e6_dp,'Invalid time interval')
  call require(cfl>0.and.cfl<=1.and.max_steps>0.and.min(rtol,atol)>0,'Invalid integration controls')
  call require(number>=0.and.number<=125000.and.total_mass>=0,'Invalid injected population')
  call require(injection_events>=0.and.injection_events<=1000000,'Invalid injection event count')
  call require(injection_start>=0.and.injection_start<=1e6_dp.and.injection_interval>0.and. &
    injection_interval<=1e6_dp,'Invalid injection timing')
  call require(injection_start+real(max(0,injection_events-1),dp)*injection_interval<=1e6_dp, &
    'Injection schedule exceeds time range')
  if (injection_events>1) call require(injection_start+injection_interval>injection_start, &
    'Injection interval is below time resolution')
  call require(injection_profile_points==0.or.(injection_profile_points>=2.and.injection_profile_points<=64), &
    'Invalid injection profile point count')
  call require(all(injection_profile_times(injection_profile_points+1:)==0).and. &
    all(injection_profile_rates(injection_profile_points+1:)==0),'Unused injection profile entries')
  call build_schedule()
  call require(all(port_profile_points>=0).and.all(port_profile_points<=64).and. &
    all(port_profile_points/=1),'Invalid port profile point count')
  call require(all(ieee_is_finite(port_profile_times)).and.all(ieee_is_finite(port_profile_rates)), &
    'Nonfinite port profile')
  call require(all(port_profile_points(port_count+1:)==0),'Profile supplied for inactive port')
  do p=1,16
    n=port_profile_points(p)
    call require(all(port_profile_times(n+1:,p)==0).and.all(port_profile_rates(n+1:,p)==0), &
      'Unused port profile entries')
    if (n==0) cycle
    call require(port_events(p)>0,'A port profile requires scheduled events')
    call integrate_spray_injection_profile(port_profile_times(:n,p),port_profile_rates(:n,p),port_start(p), &
      port_start(p)+real(port_events(p),dp)*port_interval(p),profile_mass,ok)
    call require(ok,'Invalid port profile or incomplete schedule coverage')
  end do
  if (injection_profile_points>0) then
    call require(schedule_count>0,'A profile requires scheduled events')
    ! Validate the shared table even when every active port overrides it.
    call integrate_spray_injection_profile(injection_profile_times(:injection_profile_points), &
      injection_profile_rates(:injection_profile_points),injection_profile_times(1), &
      injection_profile_times(injection_profile_points),profile_mass,ok)
    call require(ok,'Invalid shared profile')
    do p=1,schedule_count
      n=event_ports(p)
      if (n>0) then
        if (port_profile_points(n)>0) cycle
      end if
      call integrate_spray_injection_profile(injection_profile_times(:injection_profile_points), &
        injection_profile_rates(:injection_profile_points),event_times(p),event_ends(p),profile_mass,ok)
      call require(ok,'Incomplete shared profile coverage')
    end do
  end if
  call require(gas_density>0.and.gas_temperature>0.and.all(gas_fractions>=0),'Invalid initial gas')
  call require(abs(sum(gas_fractions)-1)<1e-12_dp,'Gas fractions must sum to one')
  call require(len_trim(output_prefix)>0,'Empty output prefix')
  call spray_cartesian_block(cells,[nranks,1,1],rank,first,extent,ok)
  call require(ok.and.extent(1)>=2.and.all(extent(2:)>=3),'Grid requires two x cells per rank and three y/z cells')
  primitive=0;primitive(1)=gas_density;primitive(2:4)=gas_velocity;primitive(6:)=gas_fractions
  primitive(5)=gas_density*gas_temperature*mixture_specific_gas_constant(species,gas_fractions,ok)
  call require(ok,'Invalid initial composition')
  call reactive_primitive_to_conserved(species,primitive,gas_cell,initial_t,sound,ok)
  call require(ok,'Cannot initialize gas')
  allocate(state(5+size(species),extent(1),extent(2),extent(3)),temperature(extent(1),extent(2),extent(3)))
  do k=1,extent(3)
    do j=1,extent(2)
      do i=1,extent(1)
        state(:,i,j,k)=gas_cell
      end do
    end do
  end do
  temperature=initial_t
  do n=1,2
    call initialize_spray_liquid(species,trim(liquid_species(n)),liquid_density(n),liquid_cp(n), &
      boiling_temperature(n),101325.0_dp,latent_heat(n),liquids(n),ok)
    call require(ok,'Invalid liquid properties')
  end do
  ! Bind all physical/evolution inputs, excluding rank count, paths, stop time
  ! and the maximum-step guard, which may change when continuing a run.
  write(checkpoint_context,'(*(g0,1x))',iostat=status) 'mpi-spray-application-v2', &
    mechanism_hash,mixture_evolution_id,cells,lower,spacing, &
    gas_density,gas_temperature,gas_velocity,gas_fractions,number,total_mass,center,axis,speed,half_angle, &
    diameter,drop_temperature,fractions,liquid_species,liquid_density,liquid_cp,boiling_temperature,latent_heat, &
    dt_max,cfl,rtol,atol,chemistry,solver,limiter,transport_flags,injection_events,injection_start,injection_interval
  call require(status==0,'Checkpoint context serialization failed')
  if (injection_profile_points>0) then
    write(profile_context,'(*(g0,1x))',iostat=status) 'linear-injection-profile-v1',injection_profile_points, &
      injection_profile_times(:injection_profile_points),injection_profile_rates(:injection_profile_points)
    call require(status==0,'Profile context serialization failed')
    call require(len_trim(checkpoint_context)+1+len_trim(profile_context)<=len(checkpoint_context), &
      'Checkpoint profile context too long')
    checkpoint_context=trim(checkpoint_context)//' '//trim(profile_context)
  end if
  if (port_count>0) then
    do p=1,port_count
      write(profile_context,'(*(g0,1x))',iostat=status) 'independent-port-v1',p,port_events(p), &
        port_start(p),port_interval(p),port_number(p),port_mass(p),port_center(:,p),port_axis(:,p), &
        port_fractions(:,p),port_speed(p),port_half_angle(p),port_diameter(p),port_temperature(p)
      call require(status==0,'Port context serialization failed')
      call require(len_trim(checkpoint_context)+1+len_trim(profile_context)<=len(checkpoint_context), &
        'Checkpoint port context too long')
      checkpoint_context=trim(checkpoint_context)//' '//trim(profile_context)
    end do
  end if
  do p=1,port_count
    n=port_profile_points(p)
    if (n==0) cycle
    write(profile_context,'(*(g0,1x))',iostat=status) 'port-linear-profile-v1',p,n, &
      port_profile_times(:n,p),port_profile_rates(:n,p)
    call require(status==0,'Port profile context serialization failed')
    call require(len_trim(checkpoint_context)+1+len_trim(profile_context)<=len(checkpoint_context), &
      'Checkpoint port profile context too long')
    checkpoint_context=trim(checkpoint_context)//' '//trim(profile_context)
  end do
  if (port_count==0) then
    call append_injector_context(0,size_options,nozzle,nozzle_radius)
  else
    do p=1,port_count
      call append_injector_context(p,port_size(p),port_nozzle(p),port_nozzle_radius(p))
    end do
  end if
  options%breakup_model=breakup_model;options%surface_tension=surface_tension;options%liquid_viscosity=liquid_viscosity
  options%khrt_b0=khrt_b0;options%khrt_b1=khrt_b1;options%khrt_c3=khrt_c3
  options%khrt_children=khrt_children;options%khrt_reference_diameter_cubed=khrt_reference_diameter_cubed
  call require(valid_spray_options(options),'Invalid breakup settings')
  if (breakup_model/='none'.or.surface_tension/=0.or.liquid_viscosity/=0.or. &
      khrt_b0/=.61_dp.or.khrt_b1/=7.or.khrt_c3/=1) then
    write(profile_context,'(*(g0,1x))',iostat=status) 'representative-breakup-v1',breakup_model, &
      surface_tension,liquid_viscosity,khrt_b0,khrt_b1,khrt_c3
    call require(status==0,'Breakup context serialization failed')
    call require(len_trim(checkpoint_context)+1+len_trim(profile_context)<=len(checkpoint_context), &
      'Checkpoint breakup context too long')
    checkpoint_context=trim(checkpoint_context)//' '//trim(profile_context)
  end if
  if (khrt_children.or.khrt_reference_diameter_cubed/=0) then
    write(profile_context,'(*(g0,1x))',iostat=status) 'mpi-kh-children-v1',khrt_children,khrt_reference_diameter_cubed
    call require(status==0,'KH child context serialization failed')
    call require(len_trim(checkpoint_context)+1+len_trim(profile_context)<=len(checkpoint_context), &
      'Checkpoint KH child context too long')
    checkpoint_context=trim(checkpoint_context)//' '//trim(profile_context)
  end if
  call require(all(ieee_is_finite([les_cs,les_prt,les_sct])),'Nonfinite SGS settings')
  call require(les_cs>=0.and.les_cs<=2.and.min(les_prt,les_sct)>=1e-6_dp.and. &
    max(les_prt,les_sct)<=1e6_dp,'Invalid SGS settings')
  if (les_cs/=0.or.les_prt/=.7_dp.or.les_sct/=.7_dp) then
    write(profile_context,'(*(g0,1x))',iostat=status) 'mpi-smagorinsky-v1',les_cs,les_prt,les_sct
    call require(status==0,'SGS context serialization failed')
    call require(len_trim(checkpoint_context)+1+len_trim(profile_context)<=len(checkpoint_context), &
      'Checkpoint SGS context too long')
    checkpoint_context=trim(checkpoint_context)//' '//trim(profile_context)
  end if
  call require(valid_reactive_boundaries_3d(boundary_faces),'Invalid gas/particle boundary faces')
  call require(valid_wall_velocity_3d(boundary_faces,wall_velocity),'Invalid tangential wall velocity')
  call require(all(ieee_is_finite(wall_temperature)).and.all(wall_temperature>0),'Invalid wall temperature')
  physical_boundaries=any(boundary_faces/='periodic')
  periodic_axes=[boundary_faces(1)=='periodic',boundary_faces(3)=='periodic',boundary_faces(5)=='periodic']
  inlet_state=spread(gas_cell,2,6);inlet_temperature=initial_t;boundary_total=0
  if (physical_boundaries) then
    allocate(active_faces(6));active_faces=boundary_faces
    if (any(boundary_faces=='isothermal_wall')) call require(transport_flags(2),'Isothermal walls require conduction')
    if (any(wall_velocity/=0)) call require(transport_flags(1),'Moving gas walls require viscosity')
    write(profile_context,'(*(g0,1x))',iostat=status) 'mpi-physical-boundaries-v1',boundary_faces, &
      wall_temperature,wall_velocity
    call require(status==0,'Boundary context serialization failed')
    call require(len_trim(checkpoint_context)+1+len_trim(profile_context)<=len(checkpoint_context), &
      'Checkpoint boundary context too long')
    checkpoint_context=trim(checkpoint_context)//' '//trim(profile_context)
  end if
  allocate(checkpoint_films(0),checkpoint_film_mass(2,0))
  allocate(parcels(0),masses(2,0));local_mass=0
  time=0;step=0;substeps=0;event=0;initial_totals=0;injected=0
  if (len_trim(restart_path)>0) then
    call checkpoint(.false.,trim(restart_path))
    call require(ok,'Restart rejected: invalid data or incompatible model/input context')
    call require(time<=final_time.and.event>=0.and.event<=schedule_count, &
      'Restart time/event is incompatible with this run')
    if (event>0) call require(event_times(event)<=time, &
      'Restart contains future injection events')
    if (event<schedule_count) call require(next_injection_time()>time, &
      'Restart is missing a due injection event')
  else
    call budgets();initial_totals=totals
    call inject_due_events()
  end if
  if (rank==0) then
    open(newunit=history_unit,file=trim(output_prefix)//'.history.csv',status='new',action='write',iostat=status)
    call require(status==0,'Cannot create history (existing files are not overwritten)')
    write(history_unit,'(a)') '# mechanism_sha256='//trim(mechanism_hash)
    write(history_unit,'(a,i0)') '# initial_id_high_watermark=',last_id
    write(history_unit,'(a)',advance='no') 'step,time,parcels,substeps,mass,mx,my,mz,energy,'// &
      'relative_mass_error,relative_energy_error,'// &
      'events,last_id,injected_mass,injected_energy'
    if (physical_boundaries) write(history_unit,'(a)',advance='no') ',boundary_mass,boundary_energy'
    write(history_unit,'(a)') ''
  end if
  call budgets();call history()
  do while(time<final_time)
    call require(step<max_steps,'Maximum time steps reached')
    call compute_reactive_cfl_timestep_3d(species,state,temperature,extent(1),extent(2),extent(3), &
      spacing(1),spacing(2),spacing(3),cfl,local_dt,ok)
    call require(ok,'Hydro time-step calculation failed')
    if (any(transport_flags(:3))) then
      call reactive_transport_timestep_3d(species,transport,state,temperature,extent(1),extent(2),extent(3), &
        spacing(1),spacing(2),spacing(3),cfl,transport_flags(1),transport_flags(2),transport_flags(3), &
        transport_dt,max_diffusivity,ok)
      call require(ok,'Transport time-step calculation failed')
      local_dt=min(local_dt,transport_dt)
    end if
    local_dt=min(local_dt,dt_max,final_time-time)
    if (event<schedule_count) local_dt=min(local_dt,next_injection_time()-time)
    call MPI_Allreduce(local_dt,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
    call require(ierr==MPI_SUCCESS,'Time-step reduction failed')
    ok=.false.
    do attempt=1,20
      if (dt<=0.or.time+dt<=time) exit
      call advance_mpi_mixture_reaction(MPI_COMM_WORLD,species,reactions,transport,liquids,options, &
        state,temperature,parcels,masses,lower,spacing,cells,[nranks,1,1],dt,rtol,atol,chemistry,substeps,ok, &
        solver,limiter,transport_flags,last_id,[les_cs,les_prt,les_sct],active_faces,boundary_step, &
        inlet_state,inlet_temperature,wall_temperature,wall_velocity)
      if (ok) exit
      dt=.5_dp*dt
    end do
    call require(ok,'Coupled update rejected after retries')
    boundary_total=boundary_total+boundary_step
    if (event<schedule_count.and.dt==next_injection_time()-time) then
      time=next_injection_time()
    else
      time=time+dt
    end if
    step=step+1
    call inject_due_events()
    call budgets();call history()
  end do
  if (rank==0) close(history_unit)
  if (len_trim(checkpoint_path)>0) then
    call checkpoint(.true.,trim(checkpoint_path))
    call require(ok,'Checkpoint write rejected (existing destinations are not overwritten)')
  end if
  write(rank_text,'(i6.6)') rank
  open(newunit=unit,file=trim(output_prefix)//'.rank'//rank_text//'.gas.csv',status='new',action='write',iostat=status)
  call require(status==0,'Cannot create gas output')
  write(unit,'(a)',advance='no') 'i,j,k,temperature,rho,mx,my,mz,energy'
  do n=1,size(species)
    write(unit,'(a)',advance='no') ',rhoY_'//trim(species(n)%name)
  end do
  write(unit,*)
  do k=1,extent(3)
    do j=1,extent(2)
      do i=1,extent(1)
        write(unit,'(*(g0,:,","))') first(1)+i,j,k,temperature(i,j,k),state(:,i,j,k)
      end do
    end do
  end do
  close(unit)
  open(newunit=unit,file=trim(output_prefix)//'.rank'//rank_text//'.particles.csv',status='new',action='write',iostat=status)
  call require(status==0,'Cannot create particle output')
  write(unit,'(a)') 'id,x,y,z,u,v,w,mass,temperature,multiplicity,component_mass_1,component_mass_2,'// &
    'deformation,deformation_rate,rt_time,parent_id,khrt_initial_count,khrt_shed_volume'
  do p=1,size(parcels)
    write(unit,'(*(g0,:,","))') parcels(p)%id,parcels(p)%position,parcels(p)%velocity,parcels(p)%mass, &
      parcels(p)%temperature,parcels(p)%multiplicity,masses(:,p), &
      parcels(p)%deformation,parcels(p)%deformation_rate,parcels(p)%rt_time, &
      parcels(p)%parent_id,parcels(p)%khrt_initial_count,parcels(p)%khrt_shed_volume
  end do
  close(unit)
  if (rank==0) print *, 'MPI spray completed: steps=',step,' time=',time
  call MPI_Finalize(ierr)
contains
  subroutine append_injector_context(port,sizes,outlet,radius)
    integer, intent(in) :: port
    type(spray_size_options), intent(in) :: sizes
    type(spray_nozzle_model), intent(in) :: outlet
    real(dp), intent(in) :: radius
    call require(all(ieee_is_finite([sizes%shape,sizes%std,sizes%minimum,sizes%maximum, &
      radius,outlet%swirl,outlet%spread])),'Nonfinite injector options')
    ! Default settings preserve the preceding checkpoint identity.
    if (sizes%distribution=='weibull'.and.sizes%weighting=='mass'.and.sizes%std==0.and.sizes%shape==0.and. &
        sizes%minimum==1e-12_dp.and.sizes%maximum==1.and.radius==0.and.outlet%profile=='legacy'.and. &
        outlet%swirl==0.and.outlet%spread==0) return
    write(profile_context,'(*(g0,1x))',iostat=status) 'injector-options-v1',port,sizes%distribution,sizes%weighting, &
      sizes%std,sizes%shape,sizes%minimum,sizes%maximum,radius,outlet%profile,outlet%swirl,outlet%spread
    call require(status==0,'Injector context serialization failed')
    call require(len_trim(checkpoint_context)+1+len_trim(profile_context)<=len(checkpoint_context), &
      'Checkpoint injector context too long')
    checkpoint_context=trim(checkpoint_context)//' '//trim(profile_context)
  end subroutine
  subroutine build_schedule()
    integer :: completed(16),q,e,selected
    real(dp) :: candidate_time,best_time
    call require(port_count>=0.and.port_count<=16,'Invalid port count')
    if (port_count==0) then
      schedule_count=injection_events
      allocate(event_times(schedule_count),event_ends(schedule_count),event_ports(schedule_count))
      event_ports=0
      do e=1,schedule_count
        event_times(e)=injection_start+real(e-1,dp)*injection_interval
        event_ends(e)=injection_start+real(e,dp)*injection_interval
      end do
      return
    end if
    do q=1,port_count
      call require(all(ieee_is_finite([port_start(q),port_interval(q),port_mass(q),port_center(:,q), &
        port_axis(:,q),port_fractions(:,q),port_speed(q),port_half_angle(q),port_diameter(q), &
        port_temperature(q)])),'Nonfinite port input')
      call require(port_events(q)>=0.and.port_events(q)<=1000000,'Invalid port event count')
      call require(port_start(q)>=0.and.port_start(q)<=1e6_dp.and.port_interval(q)>0.and. &
        port_interval(q)<=1e6_dp,'Invalid port timing')
      call require(port_start(q)+real(max(0,port_events(q)-1),dp)*port_interval(q)<=1e6_dp, &
        'Port schedule exceeds time range')
      if (port_events(q)>1) call require(port_start(q)+port_interval(q)>port_start(q), &
        'Port interval is below time resolution')
      call require(port_number(q)>=0.and.port_number(q)<=125000.and.port_mass(q)>=0, &
        'Invalid port population')
      call require(all(port_fractions(:,q)>=0).and.abs(sum(port_fractions(:,q))-1)<1e-12_dp, &
        'Invalid port composition')
      call require(port_diameter(q)>0.and.port_temperature(q)>0.and.port_speed(q)>=0.and. &
        port_half_angle(q)>=0.and.port_half_angle(q)<.5_dp*acos(-1.0_dp).and. &
        maxval(abs(port_axis(:,q)))>0,'Invalid port geometry or droplet state')
    end do
    schedule_count=sum(port_events(:port_count))
    call require(schedule_count<=1000000,'Too many combined port events')
    allocate(event_times(schedule_count),event_ends(schedule_count),event_ports(schedule_count))
    completed=0
    do e=1,schedule_count
      best_time=huge(1.0_dp);selected=0
      do q=1,port_count
        if (completed(q)>=port_events(q)) cycle
        candidate_time=port_start(q)+real(completed(q),dp)*port_interval(q)
        ! Strict comparison preserves port-index ordering at simultaneous events.
        if (candidate_time<best_time) then
          best_time=candidate_time;selected=q
        end if
      end do
      event_times(e)=best_time;event_ports(e)=selected
      completed(selected)=completed(selected)+1
      event_ends(e)=port_start(selected)+real(completed(selected),dp)*port_interval(selected)
    end do
  end subroutine
  real(dp) function next_injection_time() result(next_time)
    next_time=huge(1.0_dp)
    if (event<schedule_count) next_time=event_times(event+1)
  end function
  subroutine inject_due_events()
    type(spray_parcel), allocatable :: candidate(:)
    real(dp), allocatable :: candidate_mass(:,:)
    real(dp) :: delta(5+size(species)),global_delta(5+size(species))
    integer(int64) :: candidate_id
    real(dp) :: amount,event_time,event_center(3),event_axis(3),event_fractions(2)
    real(dp) :: event_speed,event_angle,event_diameter,event_temperature
    integer :: event_number,port,points
    type(spray_size_options) :: event_size
    type(spray_nozzle_model) :: event_nozzle
    real(dp) :: event_radius
    do while(event<schedule_count)
      if (next_injection_time()>time) exit
      candidate=parcels;candidate_mass=masses;candidate_id=last_id;local_mass=0
      amount=total_mass;event_number=number
      event_size=size_options;event_nozzle=nozzle;event_radius=nozzle_radius
      event_center=center;event_axis=axis;event_fractions=fractions
      event_speed=speed;event_angle=half_angle;event_diameter=diameter;event_temperature=drop_temperature
      port=event_ports(event+1)
      if (port>0) then
        amount=port_mass(port);event_number=port_number(port)
        event_size=port_size(port);event_nozzle=port_nozzle(port);event_radius=port_nozzle_radius(port)
        event_center=port_center(:,port);event_axis=port_axis(:,port);event_fractions=port_fractions(:,port)
        event_speed=port_speed(port);event_angle=port_half_angle(port)
        event_diameter=port_diameter(port);event_temperature=port_temperature(port)
      end if
      points=0
      if (port>0) points=port_profile_points(port)
      if (points>0) then
        call integrate_spray_injection_profile(port_profile_times(:points,port),port_profile_rates(:points,port), &
          next_injection_time(),event_ends(event+1),amount,ok)
        call require(ok,'Port profile event integration failed')
        if (amount==0) event_number=0
      else if (injection_profile_points>0) then
        event_time=next_injection_time()
        call integrate_spray_injection_profile(injection_profile_times(:injection_profile_points), &
          injection_profile_rates(:injection_profile_points),event_time, &
          event_ends(event+1),amount,ok)
        call require(ok,'Profile event integration failed')
        if (amount==0) event_number=0
      end if
      if (rank==0) local_mass=amount
      call inject_mpi_cone_parcels(MPI_COMM_WORLD,liquids(1),local_mass,merge(event_number,0,rank==0), &
        event_center,event_axis,event_speed,event_angle,event_diameter,event_temperature, &
        options%maximum_parcels,candidate_id,candidate,delta,ok,liquids,event_fractions,candidate_mass, &
        event_size%shape,event_size%minimum,event_size%maximum,event_size%weighting,event_size%distribution, &
        event_size%std,event_radius,event_nozzle)
      call require(ok,'Scheduled injection rejected')
      call migrate_spray_parcels(MPI_COMM_WORLD,candidate,lower,spacing,cells,[nranks,1,1], &
        periodic_axes,options%maximum_parcels,ok,candidate_mass)
      call require(ok,'Scheduled particle migration rejected')
      call MPI_Allreduce(delta,global_delta,size(delta),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call require(ierr==MPI_SUCCESS,'Injection budget reduction failed')
      call move_alloc(candidate,parcels);call move_alloc(candidate_mass,masses)
      last_id=candidate_id;injected=injected+global_delta;event=event+1
    end do
  end subroutine
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (condition) return
    write(*,'(a,i0,2a)') 'rank ',rank,': ',trim(message)
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  end subroutine
  subroutine budgets()
    call mixture_system_totals(state,spacing,parcels,liquids,masses,local_totals,ok)
    call require(ok,'Inventory calculation failed')
    call MPI_Allreduce(local_totals,totals,size(totals),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(ierr==MPI_SUCCESS,'Inventory reduction failed')
    local_count=int(size(parcels),int64)
    call MPI_Allreduce(local_count,global_count,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(ierr==MPI_SUCCESS,'Particle count reduction failed')
  end subroutine
  subroutine checkpoint(writing,path)
    logical, intent(in) :: writing
    character(len=*), intent(in) :: path
    if (physical_boundaries) then
      call mpi_spray_checkpoint(MPI_COMM_WORLD,writing,path,trim(checkpoint_context),species, &
        state,temperature,parcels,masses,lower,spacing,cells,time,step,event,last_id,initial_totals,injected,ok, &
        breakup_state=breakup_model/='none',khrt_state=breakup_model=='khrt',shedding_state=khrt_children, &
        wall_films=checkpoint_films,film_mass=checkpoint_film_mass,wall_totals=boundary_total)
      if (ok) ok=size(checkpoint_films)==0
    else
      call mpi_spray_checkpoint(MPI_COMM_WORLD,writing,path,trim(checkpoint_context),species, &
        state,temperature,parcels,masses,lower,spacing,cells,time,step,event,last_id,initial_totals,injected,ok, &
        breakup_state=breakup_model/='none',khrt_state=breakup_model=='khrt',shedding_state=khrt_children)
    end if
  end subroutine

  subroutine history()
    if (rank/=0) return
    write(history_unit,'(*(g0,:,","))',advance='no') step,time,global_count,substeps,totals(:5), &
      (totals(1)-initial_totals(1)-injected(1)-boundary_total(1))/max(1e-30_dp,abs(initial_totals(1)+injected(1))), &
      (totals(5)-initial_totals(5)-injected(5)-boundary_total(5))/max(1e-30_dp,abs(initial_totals(5)+injected(5))), &
      event,last_id,injected(1),injected(5)
    if (physical_boundaries) write(history_unit,'(a,g0,a,g0)',advance='no') ',',boundary_total(1),',',boundary_total(5)
    write(history_unit,'(a)') ''
    flush(history_unit)
  end subroutine
  end subroutine
end module
