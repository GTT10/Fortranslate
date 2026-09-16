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
  use spray_parcel_mod
  use spray_partition_mod, only: spray_cartesian_block
  use spray_coupling_3d_mod, only: mixture_system_totals
  use spray_mpi_checkpoint_mod, only: mpi_spray_checkpoint
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
  type(spray_parcel), allocatable :: parcels(:)
  real(dp), allocatable :: state(:,:,:,:),temperature(:,:,:),masses(:,:)
  integer :: cells(3)=[8,4,4],number=4,max_steps=10000,injection_events=1
  real(dp) :: injection_start=0,injection_interval=1e-8_dp
  integer :: injection_profile_points=0
  real(dp) :: injection_profile_times(64)=0,injection_profile_rates(64)=0,profile_mass
  character(len=4096) :: profile_context
  real(dp) :: spacing(3)=[.001_dp,.001_dp,.001_dp],lower(3)=0,gas_density=1,gas_temperature=800
  real(dp) :: gas_velocity(3)=0
  real(dp), allocatable :: gas_fractions(:)
  real(dp) :: total_mass=1e-10_dp,center(3)=[.004_dp,.002_dp,.002_dp],axis(3)=[1.0_dp,0.0_dp,0.0_dp]
  real(dp) :: speed=10,half_angle=.1_dp,diameter=2e-5_dp,drop_temperature=330,fractions(2)=[.5_dp,.5_dp]
  real(dp) :: liquid_density(2)=[1000.0_dp,800.0_dp],liquid_cp(2)=[4184.0_dp,2000.0_dp]
  real(dp) :: boiling_temperature(2)=[373.15_dp,350.0_dp],latent_heat(2)=[2.257e6_dp,1e6_dp]
  real(dp) :: final_time=2e-8_dp,dt_max=1e-8_dp,cfl=.3_dp,rtol=1e-7_dp,atol=1e-12_dp
  logical :: transport_flags(4)=.true.
  character(len=24) :: liquid_species(2)=[character(len=24)::'H2O','H2O2'],chemistry='implicit'
  character(len=16) :: solver='hllc',limiter='mc'
  character(len=512) :: output_prefix='mpi_spray',restart_path='',checkpoint_path=''
  character(len=16384) :: checkpoint_context
  character(len=65536) :: input_record
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
    injection_profile_times,injection_profile_rates
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
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
  if (injection_profile_points>0) then
    call require(injection_events>0,'A profile requires scheduled events')
    call integrate_spray_injection_profile(injection_profile_times(:injection_profile_points), &
      injection_profile_rates(:injection_profile_points),injection_start, &
      injection_start+real(injection_events,dp)*injection_interval,profile_mass,ok)
    call require(ok,'Invalid profile or incomplete schedule coverage')
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
  allocate(parcels(0),masses(2,0));local_mass=0
  time=0;step=0;substeps=0;event=0;initial_totals=0;injected=0
  if (len_trim(restart_path)>0) then
    call mpi_spray_checkpoint(MPI_COMM_WORLD,.false.,trim(restart_path),trim(checkpoint_context),species, &
      state,temperature,parcels,masses,lower,spacing,cells,time,step,event,last_id,initial_totals,injected,ok)
    call require(ok,'Restart rejected: invalid data or incompatible model/input context')
    call require(time<=final_time.and.event>=0.and.event<=injection_events, &
      'Restart time/event is incompatible with this run')
    if (event>0) call require(injection_start+real(event-1,dp)*injection_interval<=time, &
      'Restart contains future injection events')
    if (event<injection_events) call require(next_injection_time()>time, &
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
    write(history_unit,'(a)') 'step,time,parcels,substeps,mass,mx,my,mz,energy,relative_mass_error,relative_energy_error,'// &
      'events,last_id,injected_mass,injected_energy'
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
    if (event<injection_events) local_dt=min(local_dt,next_injection_time()-time)
    call MPI_Allreduce(local_dt,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
    call require(ierr==MPI_SUCCESS,'Time-step reduction failed')
    ok=.false.
    do attempt=1,20
      if (dt<=0.or.time+dt<=time) exit
      call advance_mpi_mixture_reaction(MPI_COMM_WORLD,species,reactions,transport,liquids,options, &
        state,temperature,parcels,masses,lower,spacing,cells,[nranks,1,1],dt,rtol,atol,chemistry,substeps,ok, &
        solver,limiter,transport_flags)
      if (ok) exit
      dt=.5_dp*dt
    end do
    call require(ok,'Coupled update rejected after retries')
    if (event<injection_events.and.dt==next_injection_time()-time) then
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
    call mpi_spray_checkpoint(MPI_COMM_WORLD,.true.,trim(checkpoint_path),trim(checkpoint_context),species, &
      state,temperature,parcels,masses,lower,spacing,cells,time,step,event,last_id,initial_totals,injected,ok)
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
  write(unit,'(a)') 'id,x,y,z,u,v,w,mass,temperature,multiplicity,component_mass_1,component_mass_2'
  do p=1,size(parcels)
    write(unit,'(*(g0,:,","))') parcels(p)%id,parcels(p)%position,parcels(p)%velocity,parcels(p)%mass, &
      parcels(p)%temperature,parcels(p)%multiplicity,masses(:,p)
  end do
  close(unit)
  if (rank==0) print *, 'MPI spray completed: steps=',step,' time=',time
  call MPI_Finalize(ierr)
contains
  real(dp) function next_injection_time() result(next_time)
    next_time=huge(1.0_dp)
    if (event<injection_events) next_time=injection_start+real(event,dp)*injection_interval
  end function
  subroutine inject_due_events()
    type(spray_parcel), allocatable :: candidate(:)
    real(dp), allocatable :: candidate_mass(:,:)
    real(dp) :: delta(5+size(species)),global_delta(5+size(species))
    integer(int64) :: candidate_id
    real(dp) :: amount,event_time
    integer :: event_number
    do while(event<injection_events)
      if (next_injection_time()>time) exit
      candidate=parcels;candidate_mass=masses;candidate_id=last_id;local_mass=0
      amount=total_mass;event_number=number
      if (injection_profile_points>0) then
        event_time=next_injection_time()
        call integrate_spray_injection_profile(injection_profile_times(:injection_profile_points), &
          injection_profile_rates(:injection_profile_points),event_time, &
          injection_start+real(event+1,dp)*injection_interval,amount,ok)
        call require(ok,'Profile event integration failed')
        if (amount==0) event_number=0
      end if
      if (rank==0) local_mass=amount
      call inject_mpi_cone_parcels(MPI_COMM_WORLD,liquids(1),local_mass,merge(event_number,0,rank==0), &
        center,axis,speed,half_angle,diameter,drop_temperature,options%maximum_parcels,candidate_id,candidate,delta,ok, &
        liquids,fractions,candidate_mass)
      call require(ok,'Scheduled injection rejected')
      call migrate_spray_parcels(MPI_COMM_WORLD,candidate,lower,spacing,cells,[nranks,1,1], &
        [.true.,.true.,.true.],options%maximum_parcels,ok,candidate_mass)
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
  subroutine history()
    if (rank/=0) return
    write(history_unit,'(*(g0,:,","))') step,time,global_count,substeps,totals(:5), &
      (totals(1)-initial_totals(1)-injected(1))/max(1e-30_dp,abs(initial_totals(1)+injected(1))), &
      (totals(5)-initial_totals(5)-injected(5))/max(1e-30_dp,abs(initial_totals(5)+injected(5))), &
      event,last_id,injected(1),injected(5)
    flush(history_unit)
  end subroutine
  end subroutine
end module
