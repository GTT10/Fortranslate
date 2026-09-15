program test_spray_integration
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_full_transport
  use mixture_thermo_mod, only: mixture_specific_gas_constant, mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use spray_wall_impact_mod
  use spray_parcel_mod
  use spray_coupling_3d_mod, only: spray_system_totals,advance_mixture_coupling_3d,mixture_system_totals, &
    advance_mixture_film_coupling_3d
  use spray_reactive_3d_mod
  use spray_checkpoint_mod
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(spray_parcel), allocatable :: drops(:),before_drops(:)
  type(spray_liquid) :: liquid,materials(2)
  type(spray_options) :: options
  type(spray_wall_options) :: wall_options
  type(spray_wall_film), allocatable :: films(:),saved_films(:)
  character(len=32) :: faces(6)
  real(dp) :: state(15,5,5,5),temperature(5,5,5),before(15,5,5,5),before_t(5,5,5)
  real(dp) :: x(10),y(10),primitive(15),r,sound,local_t,totals0(15),totals1(15),injected(15),time
  real(dp) :: checkpoint_state(15,5,5,5),checkpoint_temperature(5,5,5)
  real(dp) :: force(15),accel(3),dt,expected(3),saved_force(15)
  real(dp), allocatable :: component_mass(:,:),saved_mass(:,:),film_mass(:,:),saved_film_mass(:,:)
  integer :: i,j,k,iteration,steps,event,record_status
  character(len=16384) :: checkpoint_record
  logical :: ok
  character(len=32) :: which,method
  character(len=256) :: message
  call get_command_argument(1,which)
  call load_h2o2_full_thermo(species,ok); call require(ok,'thermo')
  call load_h2o2_full_mechanism(reactions,ok); call require(ok,'reactions')
  call load_h2o2_full_transport(transport,ok); call require(ok,'transport')
  x=0.0_dp; x([1,4,10])=[2.0_dp,1.0_dp,3.76_dp]; x=x/sum(x)
  call mass_fractions_from_mole_fractions(species,x,y,ok); call require(ok,'composition')
  r=mixture_specific_gas_constant(species,y,ok); call require(ok,'R')
  primitive=0.0_dp; primitive(1)=101325.0_dp/(r*1400.0_dp); primitive(5)=101325.0_dp; primitive(6:)=y
  do k=1,5
    do j=1,5
      do i=1,5
        call reactive_primitive_to_conserved(species,primitive,state(:,i,j,k),local_t,sound,ok)
        call require(ok,'initial state')
      end do
    end do
  end do
  temperature=1400.0_dp; before=state; before_t=temperature
  call initialize_spray_liquid(species,'H2O',1000.0_dp,4184.0_dp,373.15_dp,101325.0_dp,2.257e6_dp,liquid,ok)
  call require(ok,'liquid')
  call inject_cone_parcels(liquid,1.0e-10_dp,2,1_int64,[0.002_dp,0.002_dp,0.002_dp], &
    [0.0_dp,0.0_dp,1.0_dp],10.0_dp,0.1_dp,2.0e-5_dp,330.0_dp,drops,ok)
  call require(ok,'injection'); before_drops=drops
  call spray_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,liquid,totals0)
  select case(trim(which))
  case('native','cvode')
    method='implicit'
    if (trim(which)=='cvode') method='cvode'
    ! 125 cells exceed the backend context-slot table. Three successive calls
    ! expose contexts not finalized after each independently changing cell.
    do iteration=1,3
      call advance_spray_chemistry_3d(species,reactions,state,temperature,5.0e-7_dp, &
        1.0e-7_dp,1.0e-12_dp,method,ok,message)
      call require(ok,trim(message))
    end do
    call require(all(state(1:5,:,:,:)==before(1:5,:,:,:)),'chemistry preserves conserved hydro variables')
    call require(maxval(abs(state(6:,:,:,:)-before(6:,:,:,:)))>1.0e-10_dp,'chemistry actually advances')
    call require(all(state(6:,:,:,:)>=0.0_dp),'nonnegative chemistry')
    call require(maxval(abs(sum(state(6:,:,:,:),dim=1)-state(1,:,:,:)))<1.0e-12_dp,'species closure')
    call require(maxval(abs(state(15,:,:,:)-before(15,:,:,:)))<1.0e-12_dp,'inert nitrogen unchanged')
  case('wall_coupling')
    allocate(films(0))
    faces='periodic';faces(5:6)='no_slip_wall'
    wall_options%enabled(5)=.true.
    options%drag=.false.;options%heat=.false.;options%evaporation=.false.
    options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
    do iteration=1,2
      state=before;temperature=before_t;deallocate(films);allocate(films(0))
      call inject_cone_parcels(liquid,1e-10_dp,1,1_int64,[.002_dp,.002_dp,5e-7_dp], &
        [0.0_dp,0.0_dp,-1.0_dp],merge(10.0_dp,300.0_dp,iteration==1),0.0_dp,2e-5_dp,330.0_dp,drops,ok)
      call require(ok,'wall injection')
      call spray_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,liquid,totals0,wall_films=films)
      call advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,drops, &
        [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],2e-7_dp,.false.,.false.,'implicit', &
        1e-7_dp,1e-12_dp,0.0_dp,.7_dp,.7_dp,ok,message,boundary_faces=faces,boundary_change=force, &
        wall_options=wall_options,wall_films=films)
      call require(ok,trim(message))
      call require(size(films)==1.and.drops(1)%mass==0,'deposition retained in coupled state')
      if (iteration==1) call require(size(drops)==1,'low speed fully deposits')
      if (iteration==2) then
        call require(size(drops)==5.and.all(drops(2:)%mass>0),'high speed creates four children')
        call require(all(drops(2:)%parent_id==1),'coupled child lineage')
      end if
      call spray_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,liquid,totals1,wall_films=films)
      call require(abs(totals1(1)-totals0(1))<1e-12_dp*totals0(1),'gas/drop/film mass conservation')
      call require(abs(totals1(5)-totals0(5)-force(5))<1e-12_dp*abs(totals0(5)),'coupled film energy ledger')
      call require(maxval(abs(totals1(2:4)-totals0(2:4)-force(2:4)))<1e-18_dp,'coupled wall impulse')
    end do
    saved_films=films;before_drops=drops
    time=2e-7_dp;steps=1;event=0;injected=0
    call write_spray_checkpoint('film.chk','film-context',species,state,temperature,drops,time,steps,event, &
      totals0,injected,ok,message,boundary_totals=force,film_state=.true.,wall_films=films)
    call require(ok,trim(message))
    deallocate(films);allocate(films(0));deallocate(drops);allocate(drops(0))
    call read_spray_checkpoint('film.chk','wrong-context',species,state,temperature,drops,time,steps,event, &
      totals1,injected,ok,message,boundary_totals=force,film_state=.true.,wall_films=films)
    call require(.not.ok.and.size(films)==0.and.size(drops)==0,'film checkpoint rejection transaction')
    call read_spray_checkpoint('film.chk','film-context',species,state,temperature,drops,time,steps,event, &
      totals1,injected,ok,message,boundary_totals=force,film_state=.true.,wall_films=films)
    call require(ok.and.size(films)==1.and.size(drops)==5,'film and splash lineage checkpoint')
    call require(all(films%mass==saved_films%mass).and.all(films%height==saved_films%height), &
      'exact film checkpoint geometry')
    call require(all(drops%parent_id==before_drops%parent_id),'wall-only parent IDs restored')
    allocate(component_mass(2,size(drops)),film_mass(2,size(films)))
    do i=1,size(drops)
      component_mass(:,i)=drops(i)%mass*[.6_dp,.4_dp]
    end do
    do i=1,size(films)
      film_mass(:,i)=films(i)%mass*[.2_dp,.8_dp]
    end do
    saved_mass=component_mass;saved_film_mass=film_mass
    call write_spray_checkpoint('film-mixture.chk','mixed-film-context',species,state,temperature,drops,time,steps,event, &
      totals0,injected,ok,message,boundary_totals=force,film_state=.true.,wall_films=films, &
      component_mass=component_mass,film_component_mass=film_mass)
    call require(ok,trim(message))
    film_mass=0
    call read_spray_checkpoint('film-mixture.chk','mixed-film-context',species,state,temperature,drops,time,steps,event, &
      totals1,injected,ok,message,boundary_totals=force,film_state=.true.,wall_films=films,component_mass=component_mass)
    call require(.not.ok.and.all(film_mass==0),'cannot omit mixed film inventory')
    call read_spray_checkpoint('film-mixture.chk','mixed-film-context',species,state,temperature,drops,time,steps,event, &
      totals1,injected,ok,message,boundary_totals=force,film_state=.true.,wall_films=films, &
      component_mass=component_mass,film_component_mass=film_mass)
    call require(ok.and.all(film_mass==saved_film_mass).and.all(component_mass==saved_mass),'exact film composition restart')
    open(newunit=i,file='film-mixture.chk',status='old');j=0
    do
      read(i,'(a)',iostat=record_status) checkpoint_record
      if (record_status/=0) exit
      j=j+1
    end do
    rewind(i);open(newunit=k,file='film-mixture-bad.chk',status='replace')
    do iteration=1,j
      read(i,'(a)') checkpoint_record
      if (iteration==j-1) checkpoint_record='-1 0'
      write(k,'(a)') trim(checkpoint_record)
    end do
    close(i);close(k);checkpoint_state=state;checkpoint_temperature=temperature
    call read_spray_checkpoint('film-mixture-bad.chk','mixed-film-context',species,state,temperature,drops,time,steps,event, &
      totals1,injected,ok,message,boundary_totals=force,film_state=.true.,wall_films=films, &
      component_mass=component_mass,film_component_mass=film_mass)
    call require(.not.ok.and.all(state==checkpoint_state).and.all(temperature==checkpoint_temperature), &
      'bad film composition gas rollback')
    call require(all(film_mass==saved_film_mass).and.all(component_mass==saved_mass),'bad film composition inventory rollback')
    open(newunit=i,file='film-mixture.chk',status='old');close(i,status='delete')
    open(newunit=i,file='film-mixture-bad.chk',status='old');close(i,status='delete')
    open(newunit=i,file='film.chk',status='old');close(i,status='delete')
    state=before;temperature=before_t;deallocate(films);allocate(films(0))
    call inject_cone_parcels(liquid,1e-10_dp,1,1_int64,[.002_dp,.002_dp,5e-7_dp], &
      [0.0_dp,0.0_dp,-1.0_dp],300.0_dp,0.0_dp,2e-5_dp,330.0_dp,drops,ok)
    options%breakup_model='khrt';options%khrt_children=.true.
    options%khrt_reference_diameter_cubed=1e-25_dp
    drops%khrt_initial_count=1;drops%khrt_shed_volume=1e-18_dp
    before_drops=drops;saved_films=films
    call advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,drops, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],2e-7_dp,.true.,.false.,'invalid', &
      1e-7_dp,1e-12_dp,0.0_dp,.7_dp,.7_dp,ok,message,boundary_faces=faces,boundary_change=force, &
      wall_options=wall_options,wall_films=films)
    call require(.not.ok.and.all(state==before).and.all(temperature==before_t),'late chemistry gas rollback')
    call require(size(films)==0.and.size(drops)==1.and.all(force==0),'film and child rollback')
    call require(all(drops%mass==before_drops%mass),'parent mass rollback')
    call require(all(drops%khrt_initial_count==before_drops%khrt_initial_count).and. &
      all(drops%khrt_shed_volume==before_drops%khrt_shed_volume),'wall/KH inventory rollback')
  case('mixture_film_grid')
    materials(1)=liquid
    call initialize_spray_liquid(species,'H2O2',800.0_dp,2000.0_dp,350.0_dp,101325.0_dp,1e6_dp,materials(2),ok)
    call require(ok,'mixed film material')
    deallocate(drops);allocate(drops(0),component_mass(2,0),films(1),film_mass(2,1))
    films(1)=spray_wall_film([.002_dp,.002_dp,0.0_dp],[0.0_dp,0.0_dp,1.0_dp], &
      1e-8_dp,350.0_dp,1e-3_dp,1e-5_dp,1)
    film_mass(:,1)=films(1)%mass*[.6_dp,.4_dp];saved_films=films;saved_film_mass=film_mass
    faces='slip_wall';wall_options%enabled=.true.;wall_options%evolve=.true.
    wall_options%liquid_conductivity=.6_dp;wall_options%temperature=400
    call mixture_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,materials,component_mass,totals0,ok,films,film_mass)
    call require(ok,'mixed film initial totals')
    call advance_mixture_film_coupling_3d(species,transport,materials,options,wall_options,state,temperature,films,film_mass, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],1e-9_dp,faces,steps,force,ok,message)
    call require(ok,trim(message))
    call require(steps>0.and.films(1)%mass<saved_films(1)%mass,'coupled mixed film evolves')
    call mixture_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,materials,component_mass,totals1,ok,films,film_mass)
    call require(ok,'mixed film final totals')
    call require(abs(totals1(1)-totals0(1))<1e-12_dp*totals0(1),'mixed film grid mass')
    call require(maxval(abs(totals1(6:)-totals0(6:)))<1e-18_dp,'mixed film grid species')
    call require(abs(totals1(5)-totals0(5)-force(5))<1e-12_dp*abs(totals0(5)),'mixed film grid energy')
    state=before;temperature=before_t;films=saved_films;film_mass=saved_film_mass
    options%maximum_substeps=1
    call advance_mixture_film_coupling_3d(species,transport,materials,options,wall_options,state,temperature,films,film_mass, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],1e-4_dp,faces,steps,force,ok,message)
    call require(.not.ok.and.steps==0.and.all(force==0),'mixed film substep rejection')
    call require(all(state==before).and.all(temperature==before_t),'mixed film gas rollback')
    call require(all(film_mass==saved_film_mass).and.all(films%mass==saved_films%mass),'mixed film inventory rollback')
    options%maximum_substeps=100000;wall_options%minimum_height=1e-3_dp
    call advance_mixture_film_coupling_3d(species,transport,materials,options,wall_options,state,temperature,films,film_mass, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],1e-9_dp,faces,steps,force,ok,message)
    call require(.not.ok.and.all(state==before).and.all(film_mass==saved_film_mass), &
      'terminal film energy rejection rolls back')
    films(1)%mass=1e-12_dp;film_mass=film_mass*1e-4_dp
    call mixture_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,materials,component_mass,totals0,ok,films,film_mass)
    call require(ok,'admissible terminal initial totals')
    call advance_mixture_film_coupling_3d(species,transport,materials,options,wall_options,state,temperature,films,film_mass, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],1e-9_dp,faces,steps,force,ok,message)
    call require(ok.and.size(films)==0.and.size(film_mass,2)==0,'terminal mixed film removed with matrix column')
    call mixture_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,materials,component_mass,totals1,ok,films,film_mass)
    call require(ok.and.maxval(abs(totals1-totals0-force))<1e-12_dp,'terminal mixed film grid closure')
  case('mixture_reactive')
    materials(1)=liquid
    call initialize_spray_liquid(species,'H2O2',800.0_dp,2000.0_dp,350.0_dp,101325.0_dp,1e6_dp,materials(2),ok)
    call require(ok,'reactive mixture properties')
    allocate(component_mass(2,size(drops)))
    do i=1,size(drops)
      component_mass(:,i)=drops(i)%mass*[.6_dp,.4_dp]
    end do
    saved_mass=component_mass;accel=[0.0_dp,0.0_dp,-9.81_dp]
    call mixture_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,materials,component_mass,totals0,ok)
    call require(ok,'reactive mixture initial totals')
    call advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,drops, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],2e-7_dp,.true.,.true.,'implicit', &
      1e-7_dp,1e-12_dp,.1_dp,.7_dp,.7_dp,ok,message,acceleration=accel,source_change=force, &
      mixture_components=materials,component_mass=component_mass)
    call require(ok,trim(message))
    call require(all(component_mass<saved_mass),'reactive driver advances component evaporation')
    call require(maxval(abs(state(6,:,:,:)-before(6,:,:,:)))>1e-12_dp,'mixture gas chemistry advances')
    call mixture_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,materials,component_mass,totals1,ok)
    call require(ok,'reactive mixture updated totals')
    call require(abs(totals1(1)-totals0(1))<1e-12_dp*totals0(1),'reactive mixture mass conservation')
    call require(maxval(abs(totals1(2:4)-totals0(2:4)-force(2:4)))<1e-18_dp,'reactive mixture force impulse')
    call require(abs(totals1(5)-totals0(5)-force(5))<1e-12_dp*abs(totals0(5)),'reactive mixture force work')
    time=2e-7_dp;steps=1;event=0;injected=0
    call write_spray_checkpoint('mixture-reactive.chk','mixture-reactive-context',species,state,temperature, &
      drops,time,steps,event,totals0,injected,ok,message,source_totals=force,component_mass=component_mass)
    call require(ok,trim(message))
    call advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,drops, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],2e-7_dp,.true.,.true.,'implicit', &
      1e-7_dp,1e-12_dp,.1_dp,.7_dp,.7_dp,ok,message,acceleration=accel,source_change=force, &
      mixture_components=materials,component_mass=component_mass)
    call require(ok,trim(message))
    before=state;before_t=temperature;before_drops=drops;saved_mass=component_mass;saved_force=force
    call read_spray_checkpoint('mixture-reactive.chk','mixture-reactive-context',species,state,temperature, &
      drops,time,steps,event,totals0,injected,ok,message,source_totals=force,component_mass=component_mass)
    call require(ok,trim(message))
    call advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,drops, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],2e-7_dp,.true.,.true.,'implicit', &
      1e-7_dp,1e-12_dp,.1_dp,.7_dp,.7_dp,ok,message,acceleration=accel,source_change=force, &
      mixture_components=materials,component_mass=component_mass)
    call require(ok,trim(message))
    call require(all(state==before).and.all(temperature==before_t),'reactive mixture restart gas equality')
    call require(all(component_mass==saved_mass).and.all(force==saved_force),'reactive mixture restart composition and force')
    do i=1,size(drops)
      call require(drops(i)%mass==before_drops(i)%mass.and.drops(i)%temperature==before_drops(i)%temperature.and. &
        all(drops(i)%position==before_drops(i)%position).and.all(drops(i)%velocity==before_drops(i)%velocity), &
        'reactive mixture restart parcel equality')
    end do
    open(newunit=i,file='mixture-reactive.chk',status='old');close(i,status='delete')
    call advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,drops, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],2e-7_dp,.true.,.true.,'invalid', &
      1e-7_dp,1e-12_dp,.1_dp,.7_dp,.7_dp,ok,message,acceleration=accel,source_change=force, &
      mixture_components=materials,component_mass=component_mass)
    call require(.not.ok.and.all(state==before).and.all(temperature==before_t),'late mixture chemistry rollback')
    call require(all(component_mass==saved_mass).and.all(force==0),'late mixture component and force rollback')
    do i=1,size(drops)
      call require(drops(i)%mass==before_drops(i)%mass.and.drops(i)%temperature==before_drops(i)%temperature.and. &
        all(drops(i)%position==before_drops(i)%position).and.all(drops(i)%velocity==before_drops(i)%velocity), &
        'late mixture parcel rollback')
    end do
    call advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,drops, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],2e-7_dp,.false.,.false.,'implicit', &
      1e-7_dp,1e-12_dp,0.0_dp,.7_dp,.7_dp,ok,message,component_mass=component_mass)
    call require(.not.ok.and.all(component_mass==saved_mass),'missing mixture property table rejected')
  case('mixture_grid')
    materials(1)=liquid
    call initialize_spray_liquid(species,'H2O2',800.0_dp,2000.0_dp,350.0_dp,101325.0_dp,1e6_dp,materials(2),ok)
    call require(ok,'mixture grid components')
    allocate(component_mass(2,size(drops)))
    do i=1,size(drops)
      component_mass(:,i)=drops(i)%mass*[.6_dp,.4_dp]
    end do
    saved_mass=component_mass
    call mixture_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,materials,component_mass,totals0,ok)
    call require(ok,'initial mixture inventory')
    call advance_mixture_coupling_3d(species,transport,materials,options,state,temperature, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],drops,component_mass,2e-7_dp,steps,ok,message)
    call require(ok,trim(message))
    call require(all(component_mass<saved_mass),'both liquid species evaporate into grid')
    call require(sum(state(5+materials(1)%vapor_index,:,:,:))>0.and. &
      sum(state(5+materials(2)%vapor_index,:,:,:))>0,'two vapor sources reach gas cells')
    call mixture_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,materials,component_mass,totals1,ok)
    call require(ok,'updated mixture inventory')
    call require(abs(totals1(1)-totals0(1))<1e-12_dp*totals0(1),'mixture grid mass conservation')
    call require(maxval(abs(totals1(2:4)-totals0(2:4)))<1e-18_dp,'mixture grid momentum conservation')
    call require(abs(totals1(5)-totals0(5))<1e-12_dp*abs(totals0(5)),'mixture grid energy conservation')
    do i=1,2
      j=5+materials(i)%vapor_index
      call require(abs(totals1(j)-totals0(j))<1e-12_dp*totals0(j),'component inventory conservation')
    end do
    time=2e-7_dp;event=0;injected=0
    call write_spray_checkpoint('mixture.chk','mixture-context',species,state,temperature,drops,time,steps,event, &
      totals0,injected,ok,message,component_mass=component_mass)
    call require(ok,trim(message))
    call advance_mixture_coupling_3d(species,transport,materials,options,state,temperature, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],drops,component_mass,2e-7_dp,steps,ok,message)
    call require(ok,trim(message))
    before=state;before_t=temperature;before_drops=drops;saved_mass=component_mass
    call read_spray_checkpoint('mixture.chk','mixture-context',species,state,temperature,drops,time,steps,event, &
      totals0,injected,ok,message)
    call require(.not.ok.and.all(state==before),'mixture checkpoint cannot silently drop composition')
    open(newunit=i,file='mixture.chk',status='old',action='read')
    open(newunit=j,file='mixture-bad.chk',status='new',action='write')
    k=0
    do
      read(i,'(a)',iostat=record_status) checkpoint_record
      if (record_status/=0) exit
      k=k+1
      if (k==size(temperature)+9) checkpoint_record='0 0' ! First parcel's component mass record.
      write(j,'(a)') trim(checkpoint_record)
    end do
    close(i);close(j)
    call read_spray_checkpoint('mixture-bad.chk','mixture-context',species,state,temperature,drops,time,steps,event, &
      totals0,injected,ok,message,component_mass=component_mass)
    call require(.not.ok.and.all(state==before).and.all(component_mass==saved_mass),'malformed component state rollback')
    open(newunit=i,file='mixture-bad.chk',status='old');close(i,status='delete')
    call read_spray_checkpoint('mixture.chk','mixture-context',species,state,temperature,drops,time,steps,event, &
      totals0,injected,ok,message,component_mass=component_mass)
    call require(ok,trim(message))
    call advance_mixture_coupling_3d(species,transport,materials,options,state,temperature, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],drops,component_mass,2e-7_dp,steps,ok,message)
    call require(ok,trim(message))
    call require(all(state==before).and.all(temperature==before_t),'mixture restart gas equality')
    call require(all(component_mass==saved_mass).and.all(drops%mass==before_drops%mass),'mixture restart composition equality')
    do i=1,size(drops)
      call require(all(drops(i)%position==before_drops(i)%position).and. &
        all(drops(i)%velocity==before_drops(i)%velocity).and. &
        drops(i)%temperature==before_drops(i)%temperature,'mixture restart particle equality')
    end do
    open(newunit=i,file='mixture.chk',status='old');close(i,status='delete')
    options%maximum_substeps=1
    call advance_mixture_coupling_3d(species,transport,materials,options,state,temperature, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],drops,component_mass,1e-4_dp,steps,ok,message, &
      boundary_change=force)
    call require(.not.ok.and.all(state==before).and.all(temperature==before_t),'mixture grid late rollback')
    call require(all(component_mass==saved_mass).and.all(force==0),'mixture component and boundary rollback')
    options%maximum_substeps=100000;options%heat=.false.;options%evaporation=.false.;options%drag=.false.
    faces='periodic';faces(1:2)='outflow'
    do i=1,size(drops)
      drops(i)%position(1)=.004999_dp;drops(i)%velocity=[100.0_dp,0.0_dp,0.0_dp]
    end do
    call mixture_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,materials,component_mass,totals0,ok)
    call require(ok,'outflow mixture inventory')
    call advance_mixture_coupling_3d(species,transport,materials,options,state,temperature, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],drops,component_mass,2e-7_dp,steps,ok,message, &
      boundary_faces=faces,boundary_change=force)
    call require(ok.and.all(drops%mass==0).and.all(component_mass==0),'mixed parcel outflow clears components')
    call mixture_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,materials,component_mass,totals1,ok)
    call require(ok.and.maxval(abs(totals1-totals0-force))<1e-12_dp,'mixture outflow ledger')
    drops=before_drops;component_mass=saved_mass;state=before;temperature=before_t
    options%breakup_model='khrt';options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
    options%maximum_multiplicity=10
    do i=1,size(drops)
      drops(i)%velocity=[1000.0_dp,0.0_dp,0.0_dp]
    end do
    call mixture_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,materials,component_mass,totals0,ok)
    call require(ok,'breakup mixture inventory')
    call advance_mixture_coupling_3d(species,transport,materials,options,state,temperature, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],drops,component_mass,2e-7_dp,steps,ok,message)
    call require(ok,trim(message))
    call require(size(drops)>size(before_drops).and.all(drops%multiplicity<=10),'mixture computational splitting')
    call require(size(component_mass,2)==size(drops),'split component columns follow parcel lineage')
    call mixture_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,materials,component_mass,totals1,ok)
    call require(ok.and.maxval(abs(totals1-totals0))<1e-12_dp,'mixture breakup/splitting conservation')
    drops=before_drops;component_mass=saved_mass;state=before;temperature=before_t
    options%khrt_children=.true.;options%khrt_reference_diameter_cubed=1e-25_dp;options%maximum_multiplicity=0
    do i=1,size(drops)
      drops(i)%velocity=[1000.0_dp,0.0_dp,0.0_dp]
    end do
    before_drops=drops
    call mixture_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,materials,component_mass,totals0,ok)
    call require(ok,'mixed KH initial totals')
    call advance_mixture_coupling_3d(species,transport,materials,options,state,temperature, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],drops,component_mass,2e-7_dp,steps,ok,message)
    call require(ok.and.size(drops)>size(before_drops),'mixed KH children created')
    call mixture_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,materials,component_mass,totals1,ok)
    call require(ok.and.maxval(abs(totals1-totals0))<1e-12_dp,'mixed KH gas/child conservation')
    drops=before_drops;component_mass=saved_mass;state=before;temperature=before_t
    options%maximum_parcels=size(drops)
    call advance_mixture_coupling_3d(species,transport,materials,options,state,temperature, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],drops,component_mass,2e-7_dp,steps,ok,message)
    call require(.not.ok.and.all(state==before).and.all(temperature==before_t),'mixed KH capacity gas rollback')
    call require(size(drops)==size(before_drops).and.all(component_mass==saved_mass),'mixed KH capacity composition rollback')
    call require(all(drops%khrt_shed_volume==0).and.all(drops%mass==before_drops%mass),'mixed KH parent rollback')
  case('materials')
    materials(1)=liquid
    ! Illustrative second pure liquid: tests species routing, not peroxide property accuracy.
    call initialize_spray_liquid(species,'H2O2',1400.0_dp,2600.0_dp,420.0_dp,101325.0_dp,1.5e6_dp,materials(2),ok)
    call require(ok,'second liquid')
    drops(2)%material_id=2
    before_drops=drops
    call spray_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,liquid,totals0,materials)
    call advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,drops, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],2e-7_dp,.false.,.true.,'implicit', &
      1e-7_dp,1e-12_dp,0.0_dp,.7_dp,.7_dp,ok,message,materials=materials)
    call require(ok,trim(message))
    call spray_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,liquid,totals1,materials)
    call require(abs(totals1(1)-totals0(1))<1e-12_dp*totals0(1),'mixed mass balance')
    call require(abs(totals1(5)-totals0(5))<1e-12_dp*abs(totals0(5)),'mixed energy balance')
    call require(maxval(abs(totals1(2:4)-totals0(2:4)))<1e-18_dp,'mixed momentum balance')
    do i=1,2
      call require(drops(i)%mass<before_drops(i)%mass,'both liquids evaporate')
      j=5+materials(i)%vapor_index
      call require(sum(state(j,:,:,:))>sum(before(j,:,:,:)),'both gas species receive vapor')
      call require(abs(totals1(j)-totals0(j))<1e-12_dp*abs(totals0(j)),'individual vapor inventory')
    end do
    before=state;before_t=temperature
    drops(2)%material_id=3;before_drops=drops
    call advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,drops, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],2e-7_dp,.false.,.false.,'implicit', &
      1e-7_dp,1e-12_dp,0.0_dp,.7_dp,.7_dp,ok,message,materials=materials)
    call require(.not.ok.and.all(state==before).and.all(temperature==before_t),'invalid material rollback')
    call require(all(drops%mass==before_drops%mass).and.drops(2)%material_id==3,'parcel material rollback')
    drops(2)%material_id=2
    injected=0;time=2e-7_dp;steps=1;event=0;force=0
    call write_spray_checkpoint('materials.chk','materials-context',species,state,temperature,drops,time,steps,event, &
      totals0,injected,ok,message)
    call require(.not.ok,'legacy checkpoint rejects nondefault material')
    call write_spray_checkpoint('materials.chk','materials-context',species,state,temperature,drops,time,steps,event, &
      totals0,injected,ok,message,breakup_state=.true.,split_state=.true.,khrt_state=.true., &
      source_totals=force,material_state=.true.)
    call require(ok,trim(message))
    drops%material_id=1
    call read_spray_checkpoint('materials.chk','materials-context',species,state,temperature,drops,time,steps,event, &
      totals1,injected,ok,message,breakup_state=.true.,split_state=.true.,khrt_state=.true.,source_totals=force)
    call require(.not.ok.and.all(drops%material_id==1),'material schema mismatch rolls back')
    call read_spray_checkpoint('materials.chk','materials-context',species,state,temperature,drops,time,steps,event, &
      totals1,injected,ok,message,breakup_state=.true.,split_state=.true.,khrt_state=.true., &
      source_totals=force,material_state=.true.)
    call require(ok.and.drops(1)%material_id==1.and.drops(2)%material_id==2,'material roundtrip with optional state')
    call require(all(state==before).and.all(temperature==before_t),'mixed exact state roundtrip')
    open(newunit=i,file='materials.chk',status='old');close(i,status='delete')
  case('acceleration')
    options%drag=.false.;options%heat=.false.;options%evaporation=.false.
    accel=[0.0_dp,0.0_dp,-9810.0_dp];dt=2e-7_dp
    call advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,drops, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],dt,.false.,.false.,'implicit', &
      1e-7_dp,1e-12_dp,0.0_dp,.7_dp,.7_dp,ok,message,acceleration=accel,source_change=force)
    call require(ok,trim(message))
    do i=1,size(drops)
      expected=before_drops(i)%velocity+accel*dt
      call require(maxval(abs(drops(i)%velocity-expected))<1e-13_dp,'ballistic velocity')
      expected=before_drops(i)%position+before_drops(i)%velocity*dt+0.5_dp*accel*dt**2
      call require(maxval(abs(drops(i)%position-expected))<1e-16_dp,'ballistic position')
    end do
    call require(maxval(abs(state(4,:,:,:)/state(1,:,:,:)-accel(3)*dt))<1e-13_dp,'gas acceleration')
    call require(maxval(abs(temperature-before_t))<1e-8_dp,'force preserves thermal energy')
    call spray_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,liquid,totals1)
    call require(maxval(abs(totals1(2:4)-totals0(2:4)-force(2:4)))<1e-18_dp,'external impulse ledger')
    call require(abs(totals1(5)-totals0(5)-force(5))<1e-12_dp*abs(totals0(5)),'external work ledger')
    before=state;before_t=temperature;before_drops=drops
    call advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,drops, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],dt,.true.,.false.,'invalid', &
      1e-7_dp,1e-12_dp,0.0_dp,.7_dp,.7_dp,ok,message,acceleration=accel,source_change=force)
    call require(.not.ok.and.all(state==before).and.all(force==0),'force kick rolls back after later failure')
    do i=1,size(drops)
      call require(all(drops(i)%velocity==before_drops(i)%velocity),'force parcel rollback')
    end do
  case('coupled')
    call advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,drops, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],2.0e-7_dp,.true.,.true.,'implicit', &
      1.0e-7_dp,1.0e-12_dp,0.1_dp,0.7_dp,0.7_dp,ok,message)
    call require(ok,trim(message))
    call spray_system_totals(state,[.001_dp,.001_dp,.001_dp],drops,liquid,totals1)
    call require(abs(totals1(1)-totals0(1))<1.0e-12_dp*totals0(1),'coupled mass balance')
    call require(maxval(abs(totals1(2:4)-totals0(2:4)))<1.0e-18_dp,'coupled momentum balance')
    call require(abs(totals1(5)-totals0(5))<1.0e-12_dp*abs(totals0(5)),'coupled energy balance')
    call require(sum(drops%mass)<sum(before_drops%mass),'coupled evaporation active')
    call require(maxval(abs(state(6,:,:,:)-before(6,:,:,:)))>1.0e-12_dp,'coupled chemistry active')
  case('rollback')
    ! The first spray/LES stages run before the invalid backend is encountered.
    call advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,drops, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],2.0e-7_dp,.true.,.true.,'invalid', &
      1.0e-7_dp,1.0e-12_dp,0.1_dp,0.7_dp,0.7_dp,ok,message)
    call require(.not.ok,'late chemistry rejection')
    call require(all(state==before) .and. all(temperature==before_t),'gas whole-step rollback')
    do i=1,size(drops)
      call require(all(parcel_totals(drops(i),liquid)==parcel_totals(before_drops(i),liquid)),'parcel rollback')
      call require(all(drops(i)%position==before_drops(i)%position),'position rollback')
    end do
    options%breakup_model='khrt';options%khrt_children=.true.
    options%khrt_reference_diameter_cubed=1e-25_dp
    options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
    options%drag=.false.;options%heat=.false.;options%evaporation=.false.
    do i=1,size(drops)
      drops(i)%velocity=[0.0_dp,0.0_dp,1000.0_dp]
    end do
    before_drops=drops
    ! Establish that the same first spray stage actually creates children.
    call advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,drops, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],2.0e-7_dp,.false.,.false.,'implicit', &
      1.0e-7_dp,1.0e-12_dp,0.0_dp,0.7_dp,0.7_dp,ok,message)
    call require(ok.and.size(drops)>size(before_drops),'KH children before late failure')
    state=before;temperature=before_t;drops=before_drops
    call advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,drops, &
      [0.0_dp,0.0_dp,0.0_dp],[.001_dp,.001_dp,.001_dp],2.0e-7_dp,.true.,.false.,'invalid', &
      1.0e-7_dp,1.0e-12_dp,0.0_dp,0.7_dp,0.7_dp,ok,message)
    call require(.not.ok.and.all(state==before).and.all(temperature==before_t),'KH gas rollback')
    call require(size(drops)==size(before_drops),'KH children rollback')
    call require(all(drops%khrt_initial_count==0).and.all(drops%khrt_shed_volume==0),'KH shedding rollback')
    do i=1,size(drops)
      call require(all(parcel_totals(drops(i),liquid)==parcel_totals(before_drops(i),liquid)),'KH parent rollback')
      call require(all(drops(i)%position==before_drops(i)%position),'KH parent position rollback')
    end do
  case('checkpoint')
    injected=0.0_dp; time=1.0e-6_dp; steps=3; event=0
    call write_spray_checkpoint('roundtrip.chk','integration-context',species,state,temperature,drops,time,steps,event, &
      totals0,injected,ok,message)
    call require(ok,trim(message))
    call write_spray_checkpoint('roundtrip.chk','integration-context',species,state,temperature,drops,time,steps,event, &
      totals0,injected,ok,message)
    call require(.not.ok,'existing checkpoint not overwritten')
    state=0.0_dp; temperature=0.0_dp; time=0.0_dp; steps=0; event=0; totals1=0.0_dp
    deallocate(drops); allocate(drops(0))
    call read_spray_checkpoint('roundtrip.chk','wrong-context',species,state,temperature,drops,time,steps,event, &
      totals1,injected,ok,message)
    call require(.not.ok .and. all(state==0.0_dp) .and. size(drops)==0,'invalid restart transaction')
    call read_spray_checkpoint('roundtrip.chk','integration-context',species,state,temperature,drops,time,steps,event, &
      totals1,injected,ok,message)
    call require(ok,trim(message))
    call require(all(state==before) .and. all(temperature==before_t),'restart exact state')
    call require(all(totals0==totals1) .and. size(drops)==2 .and. steps==3,'restart metadata')
    open(newunit=i,file='roundtrip.chk',status='old'); close(i,status='delete')
    drops%rt_time=3.125e-6_dp
    call write_spray_checkpoint('khrt-roundtrip.chk','khrt-context',species,state,temperature,drops,time,steps,event, &
      totals0,injected,ok,message,breakup_state=.true.,khrt_state=.true.)
    call require(ok,trim(message))
    drops%rt_time=-1
    call read_spray_checkpoint('khrt-roundtrip.chk','khrt-context',species,state,temperature,drops,time,steps,event, &
      totals1,injected,ok,message,breakup_state=.true.,khrt_state=.true.)
    call require(ok.and.all(drops%rt_time==3.125e-6_dp),'nonzero RT timer roundtrip')
    open(newunit=i,file='khrt-roundtrip.chk',status='old'); close(i,status='delete')
  case default
    error stop 'unknown integration test'
  end select
  print *,trim(which),' PASS'
contains
  subroutine require(condition,text)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: text
    if (.not.condition) then
      print *,trim(text)
      error stop 1
    end if
  end subroutine
end program test_spray_integration
