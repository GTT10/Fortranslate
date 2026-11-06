program test_serial_spray_sdc_3d
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_elementary_transport
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use elementary_kinetics_mod, only: elementary_reaction
  use constant_volume_reactor_mod, only: advance_forced_reactor_interval
  use reactive_1d_mod, only: reactive_primitive_to_conserved,reactive_conserved_to_primitive
  use spray_coupling_3d_mod, only: particle_boundary_stencil
  use mixture_thermo_mod, only: mixture_pressure
  use spray_parcel_mod, only: spray_parcel,spray_liquid,spray_options,spray_film,initialize_spray_liquid,parcel_totals, &
    parcel_diameter,inject_cone_parcels,inject_cone_ports
  use spray_mixture_mod, only: build_mixture_film,advance_mixture_frozen_gas,mixture_transfer_rates
  use spray_breakup_mod, only: advance_spray_breakup,khrt_reference_child_count,khrt_reference_uniforms, &
    construct_reference_etab_children,spray_khrt_shedding_state,spray_khrt_release,advance_khrt_shedding, &
    advance_tab_no_birth_khrt
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options,wall_film_totals
  use spray_sdc_3d_mod, only: spray_sdc_stage_3d,advance_spray_sdc_stage_3d
  use spray_reactive_3d_mod, only: advance_spray_reactive_sdc_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction) :: reactions(0)
  type(elementary_reaction), allocatable :: control_reactions(:)
  type(spray_liquid) :: liquids(1)
  type(spray_options) :: options,properties(1)
  type(spray_parcel), allocatable :: drops(:)
  type(spray_parcel) :: original,reference
  type(spray_film) :: film
  real(dp) :: transfer(1),enthalpy(1),breakup_dt
  type(spray_sdc_stage_3d) :: work
  real(dp), allocatable :: masses(:,:)
  real(dp) :: heat(4,4,4),ir(8,4,4,4)
  real(dp) :: gas(12,4,4,4),temp(4,4,4),initial(12,4,4,4),it(4,4,4),q(12),sound
  real(dp) :: rate(12,4,4,4),old_rate(12,4,4,4),wall(5),outside(12),forced(5),before(5),after(5)
  real(dp) :: boundary(12),source(12),exchange(12),inventory(12),h(3),dt,accel(3),saved_wall(5),saved_forced(5)
  character(len=16) :: faces(6)
  integer :: i,j,k,c,mode,hydro_mode,reaction_calls=0
  integer(int64) :: watermark
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok);call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok);call require(ok,'transport')
  call initialize_spray_liquid(species,'H2O',1000._dp,4184._dp,373.15_dp,101325._dp,2257000._dp,liquids(1),ok)
  call require(ok,'liquid')
  q=0;q(1)=.4_dp;q(5)=101325;q(12)=1
  do k=1,4
  do j=1,4
  do i=1,4
    call reactive_primitive_to_conserved(species,q,initial(:,i,j,k),it(i,j,k),sound,ok)
    call require(ok,'initial cell')
  enddo
  enddo
  enddo
  h=.001_dp;dt=2e-8_dp;accel=[100._dp,-200._dp,300._dp];faces='periodic'
  options%evaporation=.false.
  allocate(drops(1),masses(1,1));original%id=1;original%mass=1e-12_dp;original%multiplicity=10
  original%position=.0015_dp;original%velocity=[1._dp,0._dp,0._dp];original%temperature=330
  drops(1)=original;masses(1,1)=original%mass;gas=initial;temp=it
  before=parcel_totals(original,liquids(1))
  call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,0,2, &
    faces,drops,masses,work,rate,wall,outside,forced,ok,accel)
  call require(ok.and.work%next_stage==1,'old half kick')
  call require(all(drops(1)%position==original%position).and.drops(1)%temperature==original%temperature,'old stage privacy')
  call require(norm2(work%midpoint(1)%position-original%position-dt*work%midpoint(1)%velocity)<1e-18_dp, &
    'old half kick followed by full drift')
  old_rate=rate;gas=0;temp=0
  call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,1,2, &
    faces,drops,masses,work,rate,wall,outside,forced,ok,accel)
  call require(ok.and.all(rate==old_rate).and.work%next_stage==2,'middle iteration reuses source without sampling')
  rate=41;saved_wall=wall;saved_forced=forced
  call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,2,2, &
    faces,drops,masses,work,rate,wall,outside,forced,ok,accel)
  call require(.not.ok.and.all(rate==41).and.all(wall==saved_wall).and.all(forced==saved_forced), &
    'failed final kick preserves outputs')
  call require(work%next_stage==2.and.drops(1)%temperature==original%temperature,'failed final kick preserves work and parcel')
  gas=initial;temp=it
  call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,2,2, &
    faces,drops,masses,work,rate,wall,outside,forced,ok,accel)
  call require(ok.and.work%next_stage==0,'final kick commits')
  call require(all(drops(1)%position==work%midpoint(1)%position),'final kick has no drift')
  after=parcel_totals(drops(1),liquids(1));inventory=0
  do c=1,12
    inventory(c)=.5_dp*dt*sum(old_rate(c,:,:,:)+rate(c,:,:,:))*product(h)
  enddo
  call require(maxval(abs(inventory(:5)+after-before-forced))<1e-18_dp,'particle half kicks and source balance')
  ! Cached stages bind material properties and reject mutation transactionally.
  drops(1)=original;masses=original%mass;gas=initial;temp=it;work=spray_sdc_stage_3d()
  properties(1)=options;rate=41;wall=42;outside=43;forced=44
  call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,0,2, &
    faces,drops,masses,work,rate,wall,outside,forced,ok,breakup_properties=properties)
  call require(ok,'material properties old stage')
  old_rate=rate;saved_wall=wall;saved_forced=forced
  properties(1)%surface_tension=properties(1)%surface_tension+.01_dp
  call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,1,2, &
    faces,drops,masses,work,rate,wall,outside,forced,ok,breakup_properties=properties)
  call require(.not.ok.and.work%next_stage==1.and.all(rate==old_rate).and.all(wall==saved_wall).and. &
    all(forced==saved_forced).and.drops(1)%mass==original%mass,'material mutation rollback')
  properties(1)=options
  call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,1,2, &
    faces,drops,masses,work,rate,wall,outside,forced,ok,breakup_properties=properties)
  call require(ok.and.work%next_stage==2,'unchanged material properties resume')
  do hydro_mode=1,2
    if(allocated(control_reactions)) deallocate(control_reactions)
    if(hydro_mode==2) then
      call load_h2o2_elementary_mechanism(control_reactions,ok);call require(ok,"MOL rollback mechanism")
    else
      control_reactions=reactions
    endif
  drops(1)=original;masses=original%mass;gas=initial;temp=it
  call advance_spray_reactive_sdc_3d(species,control_reactions,transport,liquids,options,gas,temp,drops,masses, &
    [0._dp,0._dp,0._dp],h,dt,0._dp,2,.false.,[.false.,.false.,.false.,.false.], &
    [10._dp,-20._dp,30._dp],0._dp,.7_dp,.7_dp,1e-6_dp,1e-14_dp,faces, &
    boundary,source,exchange,wall,outside,forced,ok,accel,heat_release=heat,mol_mode=hydro_mode==2)
  call require(ok,'coupled serial gas particle SDC')
  call require(all(heat==0),'nonreacting particle and force exchange has zero chemical heat')
  do c=1,12
    inventory(c)=sum(gas(c,:,:,:)-initial(c,:,:,:))*product(h)
  enddo
  after=parcel_totals(drops(1),liquids(1))
  call require(maxval(abs(inventory-boundary-source-exchange))<1e-14_dp,'coupled gas exchange ledger')
  call require(maxval(abs(exchange(:5)+after-before-forced))<1e-18_dp,'coupled particle force exchange ledger')
  call require(all(wall==0).and.all(outside==0),'periodic particle boundaries')
  drops(1)=original;masses=original%mass;gas=initial;temp=it
  boundary=41;source=42;exchange=43;wall=44;outside=45;forced=46;reaction_calls=0
  heat=47;ir=0
  options%maximum_multiplicity=3;watermark=37
  options%breakup_model='etab';options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
  call advance_spray_reactive_sdc_3d(species,control_reactions,transport,liquids,options,gas,temp,drops,masses, &
    [0._dp,0._dp,0._dp],h,dt,0._dp,2,.true.,[.false.,.false.,.false.,.false.], &
    [10._dp,-20._dp,30._dp],0._dp,.7_dp,.7_dp,1e-6_dp,1e-14_dp,faces, &
    boundary,source,exchange,wall,outside,forced,ok,accel,late_reaction_failure,last_id=watermark, &
    reaction_source=ir,heat_release=heat,mol_mode=hydro_mode==2)
  call require(.not.ok.and.reaction_calls==65,'reaction failure after final private particle kick')
  call require(all(gas==initial).and.all(temp==it),'late coupled gas temperature rollback')
  call require(size(drops)==1.and.all(drops(1)%position==original%position).and. &
    all(drops(1)%velocity==original%velocity).and.drops(1)%temperature==original%temperature.and. &
    all(masses==original%mass),'late coupled particle rollback')
  call require(all(boundary==41).and.all(source==42).and.all(exchange==43).and. &
    all(wall==44).and.all(outside==45).and.all(forced==46),'late coupled all ledgers rollback')
  call require(watermark==37,'late coupled split watermark rollback')
  call require(all(heat==47).and.all(ir==0),'late coupled reaction rates and heat rollback')
  options%maximum_multiplicity=0;options%breakup_model='none'
  enddo
  do mode=1,4
    options%mixture_time_scheme='euler';options%fixed_parts=mode==2;faces='periodic'
    if(mode==3) faces(1:2)='outflow'
    if(mode==4) faces(1:2)='slip_wall'
    if(size(drops)/=1) then
      deallocate(drops,masses);allocate(drops(1),masses(1,1))
    endif
    drops(1)=original
    if(mode>=3) then
      drops(1)%position(1)=.0039995_dp;drops(1)%velocity(1)=100
    endif
    before=parcel_totals(drops(1),liquids(1));masses=original%mass
    work=spray_sdc_stage_3d();gas=initial;temp=it
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,0,1, &
      faces,drops,masses,work,rate,wall,outside,forced,ok,accel)
    call require(ok,'Euler old kick')
    old_rate=rate
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,1,1, &
      faces,drops,masses,work,rate,wall,outside,forced,ok,accel)
    call require(ok,'Euler final kick')
    after=0
    do c=1,size(drops)
      after=after+parcel_totals(drops(c),liquids(1))
    enddo
    do c=1,12
      inventory(c)=.5_dp*dt*sum(old_rate(c,:,:,:)+rate(c,:,:,:))*product(h)
    enddo
    call require(maxval(abs(inventory(:5)+after-before-forced-wall-outside(:5)))<1e-18_dp, &
      'Euler particle sources fixed constraint and outflow budget')
    if(mode==2) call require(all(drops(1)%position==original%position).and. &
      all(drops(1)%velocity==0),'fixed particle position and zero velocity')
    if(mode==3) call require(size(drops)==0.and.size(masses,2)==0.and.outside(1)<0,'outflow particle retired')
    if(mode==4) call require(size(drops)==1.and.drops(1)%velocity(1)<0.and.wall(2)<0, &
      'reflection particle wall momentum')
  enddo
  drops(1)=original;masses=original%mass;gas=initial;temp=it;faces='periodic';watermark=37
  options%maximum_multiplicity=3
  options%breakup_model='etab'
  before=parcel_totals(original,liquids(1))
  call advance_spray_reactive_sdc_3d(species,reactions,transport,liquids,options,gas,temp,drops,masses, &
    [0._dp,0._dp,0._dp],h,dt,0._dp,2,.false.,[.false.,.false.,.false.,.false.], &
    [10._dp,-20._dp,30._dp],0._dp,.7_dp,.7_dp,1e-6_dp,1e-14_dp,faces, &
    boundary,source,exchange,wall,outside,forced,ok,accel,last_id=watermark)
  call require(ok.and.size(drops)==4.and.watermark==40,'coupled split preserves watermark gap')
  call require(all(drops%multiplicity<=3).and.all(drops%id==[1_int64,38_int64,39_int64,40_int64]), &
    'bounded multiplicity and unique split IDs')
  after=0
  do c=1,size(drops)
    after=after+parcel_totals(drops(c),liquids(1))
  enddo
  call require(abs(sum(masses(1,:)*drops%multiplicity)-original%mass*original%multiplicity)<1e-25_dp, &
    'split component mass conservation')
  call require(maxval(abs(exchange(:5)+after-before-forced))<1e-18_dp,'split particle exchange conservation')
  options%maximum_multiplicity=0;options%drag=.false.;options%heat=.false.
  do mode=1,2
    options%breakup_model='etab'
    if(mode==2) options%breakup_model='khrt'
    deallocate(drops,masses);allocate(drops(1),masses(1,1))
    drops(1)=original;drops(1)%velocity=[1000._dp,0._dp,0._dp];drops(1)%deformation=1.1_dp
    reference=drops(1);breakup_dt=dt
    if(mode==2) breakup_dt=.5_dp*dt
    call build_mixture_film(species,transport,q,it(1,1,1),reference,liquids,[1._dp],options, &
      film,transfer,enthalpy,ok)
    call require(ok,'reference breakup skin')
    call advance_spray_breakup(reference,liquids(1),options,q(1),q(2:4),film%density,film%viscosity,breakup_dt,ok)
    call require(ok,'reference model-specific breakup clock')
    masses=original%mass;gas=initial;temp=it;work=spray_sdc_stage_3d()
    before=parcel_totals(drops(1),liquids(1));watermark=37
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,0,2, &
      faces,drops,masses,work,rate,wall,outside,forced,ok,last_id=watermark)
    call require(ok.and.work%midpoint(1)%mass<original%mass,'breakup changes private half-kick drop size')
    call require(work%midpoint(1)%mass==reference%mass.and. &
      work%midpoint(1)%deformation==reference%deformation.and. &
      work%midpoint(1)%deformation_rate==reference%deformation_rate,'TAB flow interval and KHRT half interval clocks')
    call require(drops(1)%mass==original%mass.and.watermark==37,'breakup old stage stays private')
    old_rate=rate
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,1,2, &
      faces,drops,masses,work,rate,wall,outside,forced,ok,last_id=watermark)
    call require(ok,'breakup cached middle stage')
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,2,2, &
      faces,drops,masses,work,rate,wall,outside,forced,ok,last_id=watermark)
    call require(ok.and.size(drops)==1.and.watermark==37,'nonspawning breakup final stage')
    call require(abs(masses(1,1)-drops(1)%mass)<1e-25_dp,'breakup per-drop component inventory')
    after=parcel_totals(drops(1),liquids(1))
    call require(maxval(abs(after-before))<1e-18_dp,'breakup preserves physical mass momentum energy')
    call require(all(rate==0).and.all(old_rate==0),'nonspawning breakup has no gas exchange')
  enddo
  call check_etab_replacement()
  call check_khrt_release()
  call check_source_cells()
  call check_primitive_sampling()
  call check_flow_subcycles()
  call check_transfer_reynolds()
  call check_tab_khrt_shared_state()
  call check_serial_films()
  print *, 'PASS: serial particle SDC stage privacy, half kicks, drift, rollback and coupled budgets'
contains
  subroutine check_serial_films()
    type(spray_wall_film), allocatable :: patches(:)
    type(spray_wall_film) :: original_film
    type(spray_wall_options) :: settings
    type(spray_options) :: controls
    type(spray_parcel), allocatable :: particles(:)
    real(dp), allocatable :: inventory_mass(:,:),particle_mass(:,:)
    real(dp) :: film_before(5),film_after(5),gas_delta(12),original_mass
    character(len=16) :: film_faces(6)
    allocate(patches(1),inventory_mass(1,1),particles(0),particle_mass(1,0))
    original_mass=1e-10_dp
    original_film%position=[0._dp,.0015_dp,.0015_dp];original_film%normal=[1._dp,0._dp,0._dp]
    original_film%mass=original_mass;original_film%temperature=330;original_film%height=1e-6_dp
    original_film%diameter=sqrt(4*original_mass/(1000._dp*acos(-1._dp)*original_film%height))
    settings%enabled(1)=.true.;settings%evolve=.true.;settings%liquid_conductivity=.6_dp
    settings%temperature=350;film_faces='periodic';film_faces(:2)='slip_wall'
    film_before=wall_film_totals(original_film,liquids(1))
    patches(1)=original_film;inventory_mass=original_mass;gas=initial;temp=it
    call advance_spray_reactive_sdc_3d(species,reactions,transport,liquids,controls,gas,temp,particles,particle_mass, &
      [0._dp,0._dp,0._dp],h,dt,0._dp,2,.false.,[.false.,.false.,.false.,.false.], &
      [0._dp,0._dp,0._dp],0._dp,.7_dp,.7_dp,1e-6_dp,1e-14_dp,film_faces, &
      boundary,source,exchange,wall,outside,forced,ok, &
      wall_options=settings,wall_films=patches,film_component_mass=inventory_mass)
    call require(ok.and.size(patches)==1,'serial film SDC commits')
    call require(patches(1)%mass<original_mass.and.patches(1)%temperature/=330,'film evaporation and heating')
    call require(inventory_mass(1,1)==patches(1)%mass,'film component inventory')
    film_after=wall_film_totals(patches(1),liquids(1))
    do c=1,12
      gas_delta(c)=sum(gas(c,:,:,:)-initial(c,:,:,:))*product(h)
    enddo
    call require(maxval(abs(gas_delta-boundary-source-exchange))<1e-14_dp,'film gas exchange ledger')
    call require(maxval(abs(exchange(:5)+film_after-film_before-wall))<1e-18_dp,'film and wall heat budget')
    call require(abs(exchange(1)-(original_mass-patches(1)%mass))<1e-25_dp,'film evaporated species mass')
    patches(1)=original_film;inventory_mass=original_mass;gas=initial;temp=it
    boundary=41;source=42;exchange=43;wall=44;outside=45;forced=46;reaction_calls=0
    call advance_spray_reactive_sdc_3d(species,reactions,transport,liquids,controls,gas,temp,particles,particle_mass, &
      [0._dp,0._dp,0._dp],h,dt,0._dp,2,.true.,[.false.,.false.,.false.,.false.], &
      [0._dp,0._dp,0._dp],0._dp,.7_dp,.7_dp,1e-6_dp,1e-14_dp,film_faces, &
      boundary,source,exchange,wall,outside,forced,ok,evaluate_reaction=late_reaction_failure, &
      wall_options=settings,wall_films=patches,film_component_mass=inventory_mass)
    call require(.not.ok.and.reaction_calls==65,'late film coupled failure')
    call require(all(gas==initial).and.all(temp==it).and.all(inventory_mass==original_mass).and. &
      patches(1)%mass==original_mass.and.patches(1)%temperature==330,'film interval state rollback')
    call require(all(boundary==41).and.all(source==42).and.all(exchange==43).and. &
      all(wall==44).and.all(outside==45).and.all(forced==46),'film interval ledger rollback')
    deallocate(particles,particle_mass,patches,inventory_mass)
    allocate(particles(1),particle_mass(1,1),patches(0),inventory_mass(1,0))
    particles(1)=original;particles(1)%position=[5e-7_dp,.0015_dp,.0015_dp]
    particles(1)%velocity=[-300._dp,0._dp,0._dp];particle_mass=original%mass
    controls%surface_tension=.072_dp;controls%liquid_viscosity=.001_dp
    controls%drag=.false.;controls%heat=.false.;controls%evaporation=.false.
    gas=initial;temp=it;watermark=37;reaction_calls=0
    boundary=41;source=42;exchange=43;wall=44;outside=45;forced=46
    call advance_spray_reactive_sdc_3d(species,reactions,transport,liquids,controls,gas,temp,particles,particle_mass, &
      [0._dp,0._dp,0._dp],h,dt,0._dp,2,.true.,[.false.,.false.,.false.,.false.], &
      [0._dp,0._dp,0._dp],0._dp,.7_dp,.7_dp,1e-6_dp,1e-14_dp,film_faces, &
      boundary,source,exchange,wall,outside,forced,ok,evaluate_reaction=late_reaction_failure,last_id=watermark, &
      wall_options=settings,wall_films=patches,film_component_mass=inventory_mass)
    call require(.not.ok.and.reaction_calls==65,'late impact coupled failure')
    call require(size(particles)==1.and.size(patches)==0.and.size(inventory_mass,2)==0.and.watermark==37.and. &
      all(gas==initial).and.all(temp==it).and.all(particle_mass==original%mass).and. &
      all(particles(1)%position==[5e-7_dp,.0015_dp,.0015_dp]),'impact birth and ID rollback')
    call require(all(boundary==41).and.all(source==42).and.all(exchange==43).and. &
      all(wall==44).and.all(outside==45).and.all(forced==46),'impact birth ledger rollback')
  end subroutine check_serial_films
  subroutine check_tab_khrt_shared_state()
    real(dp) :: radius,cd,a,momdiff,lambda,expected_mass,z,we,tau,frequency,new_radius,pi
    type(spray_parcel) :: trial,saved
    type(spray_khrt_shedding_state) :: pending
    type(spray_khrt_release) :: emitted
    integer :: policy
    type(spray_parcel), allocatable :: born(:)
    call inject_cone_parcels(liquids(1),1e-10_dp,3,1_int64,[0._dp,0._dp,0._dp],[0._dp,0._dp,1._dp], &
      10._dp,0._dp,2e-5_dp,330._dp,born,ok,diameter_shape=2._dp,diameter_weighting='number', &
      reference_breakup=.true.)
    call require(ok.and.all(born%khrt_initial_count==born%multiplicity),'number-weighted injection finalizes native N0')
    call inject_cone_ports(liquids(1),1e-10_dp,2,1_int64,reshape([0._dp,0._dp,0._dp],[3,1]), &
      reshape([0._dp,0._dp,1._dp],[3,1]),[1._dp],10._dp,0._dp,2e-5_dp,330._dp,born,ok, &
      reference_breakup=.true.)
    call require(ok.and.all(born%khrt_initial_count==born%multiplicity),'multiport forwards native N0 initialization')
    options=spray_options();options%breakup_model='etab';options%etab_parcel_factor=1
    options%maximum_multiplicity=3;options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
    trial=original;trial%mass=1e-9_dp;trial%multiplicity=.01_dp;trial%velocity=[1000._dp,0._dp,0._dp]
    trial%deformation=.4_dp;trial%deformation_rate=.01_dp;trial%khrt_initial_count=.01_dp
    saved=trial;pi=acos(-1._dp);radius=.5_dp*parcel_diameter(trial,liquids(1))
    cd=24._dp/100*(1+100._dp**(2._dp/3._dp)/6)
    a=.375_dp*cd*.4_dp*1000**2/(liquids(1)%density*radius)
    momdiff=a*(liquids(1)%density-.4_dp);lambda=options%khrt_c3*pi/sqrt(momdiff/(3*options%surface_tension))
    expected_mass=trial%mass*(lambda/radius)**3
    call advance_tab_no_birth_khrt(trial,liquids(1),options,.4_dp,[0._dp,0._dp,0._dp],100._dp,1e-8_dp,ok)
    call require(ok.and.abs(trial%mass-expected_mass)<1e-25_dp.and.trial%deformation_rate==0.and. &
      trial%deformation==saved%deformation,'TAB BM2 serves as RT clock and BM1 survives RT breakup')
    call require(trial%id==saved%id.and.trial%rt_time==saved%rt_time.and. &
      trial%khrt_initial_count==saved%khrt_initial_count,'TAB fallback preserves ID N0 and KHRT-only clock')
    trial=saved;trial%velocity=[200._dp,0._dp,0._dp];trial%deformation=-1e-12_dp;trial%deformation_rate=-7
    saved=trial
    we=.4_dp*200**2*radius/options%surface_tension
    z=options%liquid_viscosity/sqrt(liquids(1)%density*radius*options%surface_tension)
    lambda=radius*9.02_dp*(1+.45_dp*sqrt(z))*(1+.4_dp*(z*sqrt(we))**.7_dp)/ &
      (1+.865_dp*we**1.67_dp)**.6_dp
    frequency=(.34_dp+.385_dp*we*sqrt(we))/((1+z)*(1+1.4_dp*(z*sqrt(we))**.6_dp))* &
      sqrt(options%surface_tension/(liquids(1)%density*radius**3))
    tau=3.726_dp*options%khrt_b1*radius/(lambda*frequency)
    new_radius=(dt/tau*options%khrt_b0*lambda+radius)/(1+dt/tau)
    call advance_tab_no_birth_khrt(trial,liquids(1),options,.4_dp,[0._dp,0._dp,0._dp],100._dp,dt,ok)
    call require(ok.and.abs(trial%mass-saved%mass*(new_radius/radius)**3)<1e-25_dp.and. &
      trial%deformation_rate==0,'negative TAB BM2 skips RT and takes KH branch')
    call require(trial%deformation==saved%deformation+saved%khrt_initial_count*(radius**3-new_radius**3), &
      'TAB fallback retains signed BM1 and adds native N0 volume increment')
    saved=trial
    call advance_tab_no_birth_khrt(trial,liquids(1),options,.4_dp,[0._dp,0._dp,0._dp],-1._dp,dt,ok)
    call require(.not.ok.and.trial%mass==saved%mass.and.trial%deformation==saved%deformation.and. &
      trial%deformation_rate==saved%deformation_rate,'TAB fallback rejects transactionally')
    options%breakup_model='khrt';options%etab_parcel_factor=-1
    options%khrt_children=.true.;options%khrt_parcel_factor=0;options%khrt_reference_diameter_cubed=1e-25_dp
    do policy=1,2
      trial=original;trial%mass=1e-9_dp;trial%multiplicity=.01_dp
      trial%velocity=[200._dp,0._dp,0._dp];saved=trial
      pending=spray_khrt_shedding_state()
      call advance_khrt_shedding(trial,liquids(1),options,.4_dp,[0._dp,0._dp,0._dp],100._dp,dt, &
        options%khrt_reference_diameter_cubed,pending,emitted,ok,make_new_drops=policy==1)
      call require(ok.and.trial%mass<saved%mass.and.trial%mass>0.and.trial%rt_time==0, &
        'reference zero N0 still updates KH radius and enables RT')
      call require(pending%initial_count==0.and.pending%shed_volume==0.and.emitted%mass==0, &
        'reference zero N0 is retained with births allowed or suppressed')
      saved=trial
      call advance_khrt_shedding(trial,liquids(1),options,.4_dp,[0._dp,0._dp,0._dp],100._dp,dt, &
        options%khrt_reference_diameter_cubed,pending,emitted,ok,make_new_drops=policy==1)
      call require(ok.and.pending%initial_count==0.and.pending%shed_volume==0.and.emitted%mass==0, &
        'repeated reference KH updates never infer N0 from multiplicity')
    enddo
  end subroutine check_tab_khrt_shared_state
  subroutine check_transfer_reynolds()
    type(spray_parcel) :: pending,start_drop
    type(spray_parcel), allocatable :: expected_children(:)
    real(dp) :: refq(12),reft,cm(1),saved_cm(1),delta(12),initial_re,last_re,returned_re,post_re
    real(dp) :: invtau,heat,rates(1),hv(1),duration,position(3),tangent,uniforms(1),compensation(4),constraint(5)
    integer :: step
    options=spray_options();options%evaporation=.false.;options%heat=.false.;options%mixture_time_scheme='euler'
    options%breakup_model='etab';options%etab_parcel_factor=1;options%maximum_multiplicity=3
    options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
    start_drop=original;start_drop%mass=1e-12_dp;start_drop%multiplicity=10
    start_drop%velocity=[1000._dp,0._dp,0._dp];start_drop%position=.00175_dp;start_drop%deformation=1.1_dp
    call reference_gas_at(start_drop%position,refq,reft)
    call mixture_transfer_rates(species,transport,refq,reft,start_drop,liquids,[1._dp],options, &
      invtau,heat,rates,hv,ok,initial_re)
    call require(ok.and.invtau>0,'transfer initial Reynolds and timescale')
    duration=1.5_dp/invtau;pending=start_drop;cm=pending%mass
    do step=1,2
      position=pending%position
      call build_mixture_film(species,transport,refq,reft,pending,liquids,[1._dp],options,film,transfer,enthalpy,ok)
      call require(ok,'last internal transfer reference skin')
      last_re=film%density*norm2(refq(2:4)-pending%velocity)* &
        (6*pending%mass/(acos(-1._dp)*liquids(1)%density))**(1._dp/3._dp)/film%viscosity
      call advance_mixture_frozen_gas(species,transport,refq,reft,pending,liquids,cm,options,.5_dp*duration,delta,ok, &
        time_modifier=.5_dp)
      call require(ok,'two prescribed internal transfer steps')
      pending%position=position
    enddo
    reference=start_drop;saved_cm=reference%mass;returned_re=-41
    call advance_mixture_frozen_gas(species,transport,refq,reft,reference,liquids,saved_cm,options,duration,delta,ok, &
      time_modifier=.5_dp,reynolds_number=returned_re)
    call require(ok.and.all(reference%velocity==pending%velocity).and.returned_re==last_re, &
      'frozen transfer returns last internal step pre-update Reynolds')
    call build_mixture_film(species,transport,refq,reft,reference,liquids,[1._dp],options,film,transfer,enthalpy,ok)
    call require(ok,'post-transfer skin')
    post_re=film%density*norm2(refq(2:4)-reference%velocity)* &
      (6*reference%mass/(acos(-1._dp)*liquids(1)%density))**(1._dp/3._dp)/film%viscosity
    call require(returned_re<initial_re.and.returned_re>post_re,'returned Reynolds differs from initial and rebuilt final')
    reference=start_drop;saved_cm=reference%mass;returned_re=-41;options%maximum_substeps=1
    call advance_mixture_frozen_gas(species,transport,refq,reft,reference,liquids,saved_cm,options,duration,delta,ok, &
      time_modifier=.5_dp,reynolds_number=returned_re)
    call require(.not.ok.and.returned_re==-41.and.all(reference%velocity==start_drop%velocity), &
      'failed transfer preserves Reynolds output and particle')
    options%maximum_substeps=100000;options%maximum_multiplicity=0;options%etab_parcel_factor=-1
    options%breakup_model='khrt';start_drop%rt_time=1
    drops(1)=start_drop;masses=start_drop%mass;gas=initial;temp=it;work=spray_sdc_stage_3d()
    reference=start_drop;cm=reference%mass;returned_re=-41
    call advance_mixture_frozen_gas(species,transport,refq,reft,reference,liquids,cm,options,dt,delta,ok, &
      time_modifier=.5_dp,reynolds_number=returned_re)
    call require(ok,'serial KHRT reference transfer')
    reference%position=start_drop%position
    call build_mixture_film(species,transport,refq,reft,reference,liquids,[1._dp],options,film,transfer,enthalpy,ok)
    call require(ok,'serial KHRT post-transfer properties')
    call advance_spray_breakup(reference,liquids(1),options,refq(1),refq(2:4), &
      film%density,film%viscosity,.5_dp*dt,ok,reynolds_number=returned_re)
    call require(ok,'serial KHRT supplied transfer Reynolds')
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,0,1, &
      faces,drops,masses,work,rate,wall,outside,forced,ok)
    call require(ok.and.work%midpoint(1)%mass==reference%mass.and. &
      work%midpoint(1)%rt_time==reference%rt_time,'serial KHRT receives actual transfer Reynolds')
    options%breakup_model='etab';options%maximum_multiplicity=3;options%etab_parcel_factor=1
    drops(1)=start_drop;masses=start_drop%mass;reference=start_drop;cm=reference%mass;watermark=37
    call advance_mixture_frozen_gas(species,transport,refq,reft,reference,liquids,cm,options,dt,delta,ok, &
      time_modifier=.5_dp,reynolds_number=returned_re)
    call require(ok,'serial TAB reference transfer')
    reference%position=start_drop%position
    call build_mixture_film(species,transport,refq,reft,reference,liquids,[1._dp],options,film,transfer,enthalpy,ok)
    call require(ok,'serial TAB post-transfer properties')
    call advance_spray_breakup(reference,liquids(1),options,refq(1),refq(2:4), &
      film%density,film%viscosity,dt,ok,tangent,reynolds_number=returned_re)
    call require(ok.and.tangent>0,'serial TAB uses supplied transfer Reynolds for tangential speed')
    call khrt_reference_uniforms(reference%id,watermark,uniforms)
    call construct_reference_etab_children(reference,reference%velocity,1._dp,tangent,uniforms, &
      watermark+1_int64,expected_children,compensation,ok,drift_dt=dt)
    call require(ok,'serial TAB reference child with supplied Reynolds')
    work=spray_sdc_stage_3d()
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,0,1, &
      faces,drops,masses,work,rate,wall,outside,forced,ok,last_id=watermark)
    call require(ok.and.size(work%midpoint)==1.and.all(work%midpoint(1)%velocity==expected_children(1)%velocity).and. &
      work%midpoint(1)%mass==expected_children(1)%mass,'serial TAB child receives actual transfer Reynolds')
    call require(maxval(abs(.5_dp*dt*product(h)*rate(2:5,2,2,2)-delta(2:5)-compensation))<1e-18_dp, &
      'serial TAB transfer and supplied-Re child feedback combine')
    options%fixed_parts=.true.;options%breakup_model='none';options%etab_parcel_factor=-1
    options%maximum_multiplicity=0;refq(2:4)=[70._dp,-20._dp,10._dp]
    reference=start_drop;pending=start_drop;pending%velocity=0
    call build_mixture_film(species,transport,refq,reft,pending,liquids,[1._dp],options,film,transfer,enthalpy,ok)
    call require(ok,'fixed stationary reference skin')
    last_re=film%density*norm2(refq(2:4))* &
      (6*pending%mass/(acos(-1._dp)*liquids(1)%density))**(1._dp/3._dp)/film%viscosity
    returned_re=-41
    call mixture_transfer_rates(species,transport,refq,reft,reference,liquids,[1._dp],options, &
      invtau,heat,rates,hv,ok,returned_re)
    call require(ok.and.returned_re==last_re,'fixed transfer Reynolds ignores stored moving velocity')
    before=parcel_totals(reference,liquids(1));cm=reference%mass
    call advance_mixture_frozen_gas(species,transport,refq,reft,reference,liquids,cm,options,dt,delta,ok, &
      time_modifier=.5_dp,force_change=constraint,reynolds_number=returned_re)
    call require(ok.and.returned_re==last_re.and.all(reference%velocity==0).and. &
      all(reference%position==start_drop%position),'fixed half kick returns stationary Reynolds and constrains velocity')
    after=parcel_totals(reference,liquids(1))
    call require(maxval(abs(after-before+delta(:5)-constraint))<1e-18_dp,'fixed Reynolds handoff preserves constraint budget')
  end subroutine check_transfer_reynolds
  subroutine check_flow_subcycles()
    type(spray_parcel), allocatable :: expected(:)
    real(dp), allocatable :: uniforms(:)
    real(dp) :: duration,subdt,refq(12),reft,cm(1),delta(12),theta,tangent_total,compensation(4),start(3)
    real(dp) :: expected_exchange(12,4,4,4),saved_temp
    integer :: step,count,p,x,y,z,site(3)
    options=spray_options();options%evaporation=.false.;options%heat=.false.;options%drag=.false.
    options%breakup_model='etab';options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
    options%maximum_multiplicity=.001_dp;options%etab_parcel_factor=1
    drops(1)=original;drops(1)%mass=1e-9_dp;drops(1)%multiplicity=.01_dp
    drops(1)%velocity=[100._dp,0._dp,0._dp];drops(1)%position=.00175_dp;drops(1)%deformation=1.1_dp
    masses=drops(1)%mass;gas=initial;temp=it;faces='periodic';watermark=37;work=spray_sdc_stage_3d()
    duration=1.2e-5_dp;subdt=duration/3;reference=drops(1);cm=reference%mass;tangent_total=0
    before=parcel_totals(reference,liquids(1))
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,duration,0,2, &
      faces,drops,masses,work,rate,wall,outside,forced,ok,last_id=watermark)
    call require(ok.and.work%flow_subcycles==1,'pinned PeleC default CFL disables flow subdivision')
    work=spray_sdc_stage_3d()

    do step=1,3
      start=reference%position
      call reference_gas_at(start,refq,reft)
      call advance_mixture_frozen_gas(species,transport,refq,reft,reference,liquids,cm,options,.5_dp*subdt,delta,ok)
      call require(ok,'subcycle reference thermal half kick')
      reference%position=start
      call build_mixture_film(species,transport,refq,reft,reference,liquids,[1._dp],options, &
        film,transfer,enthalpy,ok)
      call require(ok,'subcycle reference skin')
      call advance_spray_breakup(reference,liquids(1),options,refq(1),refq(2:4), &
        film%density,film%viscosity,subdt,ok,theta)
      call require(ok,'TAB per-flow-subcycle update')
      tangent_total=tangent_total+theta;cm=reference%mass
      if(step<3) reference%position=reference%position+subdt*reference%velocity
    enddo
    call require(tangent_total>0,'subcycle TAB accumulates tangent')
    call khrt_reference_child_count(reference%multiplicity,1._dp,1000000,count,ok)
    call require(ok,'subcycle reference child count')
    allocate(uniforms(count));call khrt_reference_uniforms(reference%id,watermark,uniforms)
    call construct_reference_etab_children(reference,reference%velocity,1._dp,tangent_total,uniforms, &
      watermark+1_int64,expected,compensation,ok,drift_dt=subdt)
    call require(ok,'subcycle reference children move with final subinterval')
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,duration,0,2, &
      faces,drops,masses,work,rate,wall,outside,forced,ok,last_id=watermark,native_flow_cfl=1.2_dp)
    call require(ok.and.work%flow_subcycles==3.and.size(work%midpoint)==count,'old stage selects three flow subcycles')
    after=0
    do p=1,count
      call require(work%midpoint(p)%mass==expected(p)%mass.and. &
        all(work%midpoint(p)%velocity==expected(p)%velocity).and. &
        norm2(work%midpoint(p)%position-expected(p)%position)<1e-18_dp, &
        'TAB repeated updates accumulated tangent and final child drift match')
      after=after+parcel_totals(work%midpoint(p),liquids(1))
    enddo
    do c=1,12
      inventory(c)=.5_dp*duration*sum(rate(c,:,:,:))*product(h)
    enddo
    call require(maxval(abs(inventory(:5)+after-before))<1e-18_dp,'subcycle exchange compensation inventory')
    rate(:,3,2,2)=0
    call require(all(rate==0),'replacement compensation uses last subcycle parent site')
    call require(size(drops)==1.and.watermark==37,'subcycles and births stay private')
    work=spray_sdc_stage_3d();options%maximum_substeps=2;rate=41
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,duration,0,2, &
      faces,drops,masses,work,rate,wall,outside,forced,ok,last_id=watermark,native_flow_cfl=1.2_dp)
    call require(.not.ok.and.work%next_stage==0.and.all(rate==41).and.watermark==37, &
      'flow subcycle limit preserves stage and physical state')
    options=spray_options();options%evaporation=.false.;options%drag=.false.;options%mixture_time_scheme='euler'
    do z=1,4
    do y=1,4
    do x=1,4
      refq=0;refq(1)=.4_dp;refq(12)=1;theta=600._dp+100._dp*x
      refq(5)=mixture_pressure(species,refq(6:),refq(1),theta,ok)
      call require(ok,'moving sample field pressure')
      call reactive_primitive_to_conserved(species,refq,gas(:,x,y,z),temp(x,y,z),sound,ok)
      call require(ok,'moving sample field cell')
    enddo
    enddo
    enddo
    drops(1)=original;drops(1)%mass=1e-12_dp;drops(1)%multiplicity=10
    drops(1)%velocity=[100._dp,0._dp,0._dp];drops(1)%position=.00175_dp;drops(1)%deformation=0
    masses=drops(1)%mass;reference=drops(1);cm=reference%mass;expected_exchange=0
    before=parcel_totals(reference,liquids(1))
    do step=1,3
      start=reference%position;site=floor(start/h)+1
      call reference_gas_at(start,refq,reft,gas,temp)
      call advance_mixture_frozen_gas(species,transport,refq,reft,reference,liquids,cm,options,subdt,delta,ok, &
        time_modifier=.5_dp)
      call require(ok,'moving sample analytic thermal kick')
      reference%position=start+subdt*reference%velocity
      expected_exchange(:,site(1),site(2),site(3))=expected_exchange(:,site(1),site(2),site(3))+delta
    enddo
    work=spray_sdc_stage_3d()
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,duration,0,1, &
      faces,drops,masses,work,rate,wall,outside,forced,ok,native_flow_cfl=1.2_dp)
    call require(ok.and.work%flow_subcycles==3.and.work%midpoint(1)%temperature==reference%temperature.and. &
      all(work%midpoint(1)%position==reference%position),'flow subcycles resample varying temperature along trajectory')
    call require(maxval(abs(.5_dp*duration*product(h)*rate-expected_exchange))<1e-18_dp, &
      'each subcycle deposits exchange into its current containing cell')
    call reference_gas_at(reference%position,refq,reft,gas,temp)
    start=reference%position
    call advance_mixture_frozen_gas(species,transport,refq,reft,reference,liquids,cm,options,duration,delta,ok, &
      time_modifier=.5_dp)
    call require(ok,'final MK analytic full flow interval')
    reference%position=start
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,duration,1,1, &
      faces,drops,masses,work,rate,wall,outside,forced,ok,native_flow_cfl=1.2_dp)
    call require(ok.and.drops(1)%temperature==reference%temperature.and.all(drops(1)%position==start), &
      'final MK is one kick without flow subcycles or drift')
    saved_temp=drops(1)%temperature
    drops(1)%position=.00175_dp;masses=drops(1)%mass;work=spray_sdc_stage_3d();rate=41
    gas(:,4,:,:)=0;temp(4,:,:)=0
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,duration,0,1, &
      faces,drops,masses,work,rate,wall,outside,forced,ok,native_flow_cfl=1.2_dp)
    call require(.not.ok.and.work%next_stage==0.and.all(rate==41).and. &
      drops(1)%temperature==saved_temp.and.all(drops(1)%position==.00175_dp), &
      'late flow-subcycle sample failure preserves entire interval')
    options=spray_options();options%evaporation=.false.;options%heat=.false.;options%drag=.false.
    options%breakup_model='khrt';options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
    gas=initial;temp=it;drops(1)=original;drops(1)%mass=1e-9_dp;drops(1)%multiplicity=.01_dp
    drops(1)%velocity=[100._dp,0._dp,0._dp];drops(1)%position=.00175_dp;drops(1)%rt_time=-1
    masses=drops(1)%mass;reference=drops(1);work=spray_sdc_stage_3d()
    reference%position=reference%position+2*subdt*reference%velocity
    call reference_gas_at(reference%position,refq,reft)
    call build_mixture_film(species,transport,refq,reft,reference,liquids,[1._dp],options,film,transfer,enthalpy,ok)
    call require(ok,'KHRT last-subcycle reference skin')
    call advance_spray_breakup(reference,liquids(1),options,refq(1),refq(2:4), &
      film%density,film%viscosity,.5_dp*duration,ok)
    call require(ok,'KHRT reference uses full half-interval at last subcycle')
    reference%position=reference%position+subdt*reference%velocity
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,duration,0,1, &
      faces,drops,masses,work,rate,wall,outside,forced,ok,native_flow_cfl=1.2_dp)
    call require(ok.and.work%flow_subcycles==3.and.work%midpoint(1)%mass==reference%mass.and. &
      work%midpoint(1)%rt_time==reference%rt_time.and. &
      norm2(work%midpoint(1)%position-reference%position)<1e-18_dp, &
      'KHRT only updates on last subcycle using full half interval')
  end subroutine check_flow_subcycles
  subroutine check_primitive_sampling()
    real(dp) :: cellq(12),expectedq(12),delta(12),cm(1),theta
    integer :: x,y,z
    do z=1,4
    do y=1,4
    do x=1,4
      cellq=0;cellq(1)=.2_dp+.2_dp*x;cellq(2)=10._dp*x;cellq(12)=1
      theta=600._dp+100._dp*y
      cellq(5)=mixture_pressure(species,cellq(6:),cellq(1),theta,ok)
      call require(ok,'nonuniform reference pressure')
      call reactive_primitive_to_conserved(species,cellq,gas(:,x,y,z),temp(x,y,z),sound,ok)
      call require(ok,'nonuniform cell state')
    enddo
    enddo
    enddo
    expectedq=0;expectedq(1)=.65_dp;expectedq(2)=22.5_dp;expectedq(12)=1;theta=825
    expectedq(5)=mixture_pressure(species,expectedq(6:),expectedq(1),theta,ok)
    call require(ok,'analytic interpolated pressure')
    options=spray_options();options%evaporation=.false.;faces='periodic'
    drops(1)=original;drops(1)%position=.00175_dp;drops(1)%velocity=[1._dp,0._dp,0._dp]
    drops(1)%deformation=0;masses=drops(1)%mass;reference=drops(1);cm=reference%mass
    call advance_mixture_frozen_gas(species,transport,expectedq,theta,reference,liquids,cm,options,.5_dp*dt,delta,ok)
    call require(ok,'analytic primitive transfer')
    work=spray_sdc_stage_3d()
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,0,1, &
      faces,drops,masses,work,rate,wall,outside,forced,ok)
    call require(ok.and.norm2(work%midpoint(1)%velocity-reference%velocity)<1e-12_dp.and. &
      abs(work%midpoint(1)%temperature-reference%temperature)<1e-10_dp, &
      'native primitive interpolation agrees with nonuniform analytic state')
    call require(maxval(abs(.5_dp*dt*rate(:,2,2,2)*product(h)-delta))<1e-18_dp, &
      'nonuniform sampled transfer enters containing cell')
    gas=0;temp=0;gas(:,2,2,2)=initial(:,2,2,2);temp(2,2,2)=it(2,2,2)
    drops(1)=original;drops(1)%position=.0015_dp;masses=drops(1)%mass;work=spray_sdc_stage_3d()
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,0,1, &
      faces,drops,masses,work,rate,wall,outside,forced,ok)
    call require(ok,'zero interpolation weight does not read invalid unused gas cells')
  end subroutine check_primitive_sampling
  subroutine check_source_cells()
    real(dp) :: expected_position(3)
    integer :: site(3),next_site(3),variant
    options=spray_options();options%evaporation=.false.;options%drag=.false.
    do variant=1,3
      faces='periodic'
      if(variant==3) faces='slip_wall'
      deallocate(drops,masses);allocate(drops(1),masses(1,1))
      drops(1)=original;drops(1)%mass=1e-12_dp;drops(1)%multiplicity=10
      drops(1)%velocity=[1._dp,0._dp,0._dp];drops(1)%deformation=0
      drops(1)%position=[.0019999999_dp,.00175_dp,.00175_dp];site=[2,2,2]
      if(variant==2) then
        drops(1)%position(1)=.0039999999_dp;site(1)=4
      else if(variant==3) then
        drops(1)%position(1)=.0001_dp;site(1)=1
      endif
      masses=drops(1)%mass;gas=initial;temp=it;work=spray_sdc_stage_3d()
      before=parcel_totals(drops(1),liquids(1))
      call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,0,1, &
        faces,drops,masses,work,rate,wall,outside,forced,ok)
      call require(ok.and.any(rate/=0),'containing-cell old source exists')
      old_rate=rate
      rate(:,site(1),site(2),site(3))=0
      call require(all(rate==0),'old exchange deposits solely into containing cell')
      expected_position=work%midpoint(1)%position
      next_site=min(4,floor(expected_position/h)+1)
      if(variant==1) call require(next_site(1)==3,'old drift crosses interior source cell')
      if(variant==2) call require(next_site(1)==1,'old drift wraps source cell periodically')
      call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,1,1, &
        faces,drops,masses,work,rate,wall,outside,forced,ok)
      call require(ok.and.any(rate/=0),'containing-cell final source exists')
      do c=1,12
        inventory(c)=.5_dp*dt*sum(old_rate(c,:,:,:)+rate(c,:,:,:))*product(h)
      enddo
      rate(:,next_site(1),next_site(2),next_site(3))=0
      call require(all(rate==0),'final exchange deposits solely into post-drift containing cell')
      after=parcel_totals(drops(1),liquids(1))
      call require(maxval(abs(inventory(:5)+after-before-forced-wall-outside(:5)))<1e-18_dp, &
        'containing-cell source preserves coupled inventories')
    enddo
  end subroutine check_source_cells
  subroutine check_khrt_release()
    type(spray_parcel) :: parent,expected
    type(spray_khrt_shedding_state) :: shedding
    type(spray_khrt_release) :: released
    real(dp) :: radius,reynolds,refq(12),reft
    integer :: variant,p,births
    options%breakup_model='khrt';options%etab_parcel_factor=-1;options%maximum_multiplicity=0
    options%khrt_children=.true.;options%khrt_reference_diameter_cubed=1e-25_dp
    ! Keep physical inventory at the existing absolute conservation gate's scale.
    parent=original;parent%mass=1e-9_dp;parent%multiplicity=.01_dp;parent%velocity=[200._dp,0._dp,0._dp]
    parent%position=.00175_dp
    radius=(3*parent%mass/(4*acos(-1._dp)*liquids(1)%density))**(1._dp/3._dp)
    parent%khrt_initial_count=parent%multiplicity;parent%khrt_shed_volume=.95_dp*parent%multiplicity*radius**3
    parent%rt_time=-1
    do variant=1,2
      options%khrt_parcel_factor=-1
      if(variant==2) options%khrt_parcel_factor=1
      deallocate(drops,masses);allocate(drops(1),masses(1,1))
      drops(1)=parent;masses=parent%mass;gas=initial;temp=it;work=spray_sdc_stage_3d();watermark=37
      before=parcel_totals(parent,liquids(1))
      call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,0,2, &
        faces,drops,masses,work,rate,wall,outside,forced,ok,last_id=watermark)
      call require(ok.and.size(work%midpoint)>1.and.work%midpoint(1)%mass==0,'KHRT complete parent release')
      births=size(work%midpoint)-1
      call require((variant==1.and.births==4).or.(variant==2.and.births==1),'KHRT configured child count')
      call require(size(drops)==1.and.drops(1)%mass==parent%mass.and.watermark==37,'KHRT old birth privacy')
      after=0
      do p=1,size(work%midpoint)
        after=after+parcel_totals(work%midpoint(p),liquids(1))
        if(p==1) cycle
        call require(work%midpoint(p)%id==37+p-1.and.work%midpoint(p)%parent_id==parent%id, &
          'KHRT child ID reservation and parent')
        call require(norm2(work%midpoint(p)%position-parent%position-dt*work%midpoint(p)%velocity)<1e-18_dp, &
          'KHRT children drift once after birth')
      enddo
      do c=1,12
        inventory(c)=.5_dp*dt*sum(rate(c,:,:,:))*product(h)
      enddo
      if(maxval(abs(inventory(:5)+after-before))>=1e-18_dp.or.inventory(5)>=0) &
        print *, 'KHRT old balance diagnostic',variant,inventory(:5)+after-before,inventory(5)
      call require(maxval(abs(inventory(:5)+after-before))<1e-18_dp.and.inventory(5)<0, &
        'KHRT old source closes released particle inventory')
      old_rate=rate;gas=0;temp=0
      rate(:,2,2,2)=0
      call require(all(rate==0),'KHRT compensation stays in parent containing cell')
      call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,1,2, &
        faces,drops,masses,work,rate,wall,outside,forced,ok,last_id=watermark)
      call require(ok.and.all(rate==old_rate).and.watermark==37,'KHRT middle source cache')
      gas=initial;temp=it
      call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,2,2, &
        faces,drops,masses,work,rate,wall,outside,forced,ok,last_id=watermark)
      call require(ok.and.size(drops)==births.and.watermark==37+births,'KHRT final retires parent without births')
      after=0
      do p=1,size(drops)
        after=after+parcel_totals(drops(p),liquids(1))
        call require(all(drops(p)%position==work%midpoint(p+1)%position),'KHRT final has no drift')
      enddo
      do c=1,12
        inventory(c)=.5_dp*dt*sum(old_rate(c,:,:,:)+rate(c,:,:,:))*product(h)
      enddo
      call require(maxval(abs(inventory(:5)+after-before))<1e-18_dp,'KHRT complete stage physical conservation')
      call require(abs(sum(masses(1,:)*drops%multiplicity)-parent%mass*parent%multiplicity)<1e-22_dp, &
        'KHRT complete stage component conservation')
      deallocate(drops,masses);allocate(drops(1),masses(1,1))
      drops(1)=parent;masses=parent%mass;work=spray_sdc_stage_3d();watermark=37;options%maximum_parcels=1;rate=41
      call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,0,2, &
        faces,drops,masses,work,rate,wall,outside,forced,ok,last_id=watermark)
      call require(.not.ok.and.work%next_stage==0.and.watermark==37.and.all(rate==41).and. &
        drops(1)%mass==parent%mass,'KHRT birth capacity rollback')
      options%maximum_parcels=1000000
    enddo
    watermark=huge(watermark);rate=41
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,0,2, &
      faces,drops,masses,work,rate,wall,outside,forced,ok,last_id=watermark)
    call require(.not.ok.and.work%next_stage==0.and.watermark==huge(watermark).and.all(rate==41), &
      'KHRT birth ID overflow rollback')
    watermark=37
    call advance_spray_reactive_sdc_3d(species,reactions,transport,liquids,options,gas,temp,drops,masses, &
      [0._dp,0._dp,0._dp],h,dt,0._dp,2,.false.,[.false.,.false.,.false.,.false.], &
      [0._dp,0._dp,0._dp],0._dp,.7_dp,.7_dp,1e-6_dp,1e-14_dp,faces, &
      boundary,source,exchange,wall,outside,forced,ok,last_id=watermark)
    call require(ok.and.size(drops)==1.and.watermark==38,'KHRT accepted gas particle coupling')
    after=parcel_totals(drops(1),liquids(1))
    do c=1,12
      inventory(c)=sum(gas(c,:,:,:)-initial(c,:,:,:))*product(h)
    enddo
    call require(maxval(abs(inventory-boundary-source-exchange))<1e-14_dp.and. &
      maxval(abs(exchange(:5)+after-before))<1e-18_dp,'KHRT coupled gas particle exchange conservation')
    drops(1)=parent;masses=parent%mass;gas=initial;temp=it;watermark=37
    boundary=41;source=42;exchange=43;wall=44;outside=45;forced=46;reaction_calls=0
    call advance_spray_reactive_sdc_3d(species,reactions,transport,liquids,options,gas,temp,drops,masses, &
      [0._dp,0._dp,0._dp],h,dt,0._dp,2,.true.,[.false.,.false.,.false.,.false.], &
      [0._dp,0._dp,0._dp],0._dp,.7_dp,.7_dp,1e-6_dp,1e-14_dp,faces, &
      boundary,source,exchange,wall,outside,forced,ok,evaluate_reaction=late_reaction_failure,last_id=watermark)
    call require(.not.ok.and.reaction_calls==65,'KHRT late reactor failure after private births and retirement')
    call require(size(drops)==1.and.watermark==37.and.drops(1)%mass==parent%mass.and. &
      drops(1)%khrt_shed_volume==parent%khrt_shed_volume.and.all(masses==parent%mass).and. &
      all(gas==initial).and.all(temp==it),'KHRT coupled gas particle and ID rollback')
    call require(all(boundary==41).and.all(source==42).and.all(exchange==43).and. &
      all(wall==44).and.all(outside==45).and.all(forced==46),'KHRT coupled ledger rollback')
    expected=parent;shedding=spray_khrt_shedding_state(parent%khrt_initial_count,parent%khrt_shed_volume)
    call build_mixture_film(species,transport,q,it(1,1,1),expected,liquids,[1._dp],options, &
      film,transfer,enthalpy,ok)
    call require(ok,'KHRT creation-disabled skin')
    reynolds=film%density*norm2(expected%velocity)*2*radius/film%viscosity
    call advance_khrt_shedding(expected,liquids(1),options,q(1),q(2:4),reynolds,.5_dp*dt, &
      options%khrt_reference_diameter_cubed,shedding,released,ok,make_new_drops=.false.)
    call require(ok.and.released%mass==0.and.expected%mass>0.and.shedding%shed_volume>parent%khrt_shed_volume, &
      'KHRT creation-disabled branch retains even above-threshold shedding')
    ! The final MK still accumulates shed volume when creation is disabled.
    drops(1)=parent;drops(1)%khrt_shed_volume=0;masses=parent%mass;work=spray_sdc_stage_3d();watermark=37
    options%khrt_reference_diameter_cubed=1e-3_dp
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,0,1, &
      faces,drops,masses,work,rate,wall,outside,forced,ok,last_id=watermark)
    call require(ok.and.size(work%midpoint)==1.and.work%midpoint(1)%khrt_shed_volume>0, &
      'KHRT old accumulates subthreshold shedding')
    expected=work%midpoint(1)
    call reference_gas_at(expected%position,refq,reft)
    call build_mixture_film(species,transport,refq,reft,expected,liquids,[1._dp],options, &
      film,transfer,enthalpy,ok)
    call require(ok,'KHRT final reference skin')
    radius=(3*expected%mass/(4*acos(-1._dp)*liquids(1)%density))**(1._dp/3._dp)
    reynolds=film%density*norm2(expected%velocity-refq(2:4))*2*radius/film%viscosity
    shedding=spray_khrt_shedding_state(expected%khrt_initial_count,expected%khrt_shed_volume)
    call advance_khrt_shedding(expected,liquids(1),options,refq(1),refq(2:4),reynolds,.5_dp*dt, &
      options%khrt_reference_diameter_cubed,shedding,released,ok,make_new_drops=.false.)
    call require(ok.and.released%mass==0,'KHRT final reference disables creation')
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,1,1, &
      faces,drops,masses,work,rate,wall,outside,forced,ok,last_id=watermark)
    call require(ok.and.size(drops)==1.and.watermark==37.and.drops(1)%mass==expected%mass.and. &
      abs(drops(1)%khrt_shed_volume-shedding%shed_volume)<1e-27_dp.and. &
      drops(1)%khrt_shed_volume>work%midpoint(1)%khrt_shed_volume,'KHRT final updates size and accumulated shedding')
  end subroutine check_khrt_release
  subroutine check_etab_replacement()
    type(spray_parcel), allocatable :: expected(:)
    real(dp), allocatable :: uniforms(:)
    real(dp) :: tangent,compensation(4),refq(12),reft
    integer :: count,p
    options%breakup_model='etab';options%etab_parcel_factor=0;options%maximum_multiplicity=3
    options%mixture_time_scheme='conservative';options%fixed_parts=.false.
    deallocate(drops,masses);allocate(drops(1),masses(1,1))
    original%velocity=[1000._dp,0._dp,0._dp];original%deformation=1.1_dp
    original%position=.00175_dp
    drops(1)=original;masses=original%mass;gas=initial;temp=it;watermark=37
    reference=original
    call reference_gas_at(reference%position,refq,reft)
    call build_mixture_film(species,transport,refq,reft,reference,liquids,[1._dp],options, &
      film,transfer,enthalpy,ok)
    call require(ok,'ETAB replacement reference skin')
    call advance_spray_breakup(reference,liquids(1),options,refq(1),refq(2:4), &
      film%density,film%viscosity,dt,ok,tangent)
    call require(ok.and.tangent>0,'ETAB replacement reference tangent')
    call khrt_reference_child_count(reference%multiplicity,options%etab_parcel_factor, &
      options%maximum_parcels,count,ok)
    call require(ok,'ETAB reference count')
    allocate(uniforms(count));call khrt_reference_uniforms(reference%id,watermark,uniforms)
    call construct_reference_etab_children(reference,reference%velocity,options%etab_parcel_factor, &
      tangent,uniforms,watermark+1_int64,expected,compensation,ok,drift_dt=0._dp)
    call require(ok,'ETAB reference children')
    work=spray_sdc_stage_3d();before=parcel_totals(original,liquids(1))
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,0,2, &
      faces,drops,masses,work,rate,wall,outside,forced,ok,last_id=watermark)
    call require(ok.and.size(work%midpoint)==count,'ETAB replaces parent on old stage')
    call require(size(drops)==1.and.watermark==37.and.work%midpoint_watermark==37+count, &
      'ETAB births and IDs stay private')
    after=0
    do p=1,count
      call require(work%midpoint(p)%id==expected(p)%id.and. &
        all(work%midpoint(p)%velocity==expected(p)%velocity).and. &
        work%midpoint(p)%multiplicity==expected(p)%multiplicity,'ETAB constructor matches reserved children')
      call require(norm2(work%midpoint(p)%position-expected(p)%position-dt*expected(p)%velocity)<1e-18_dp, &
        'ETAB children drift once after birth')
      after=after+parcel_totals(work%midpoint(p),liquids(1))
    enddo
    do c=1,12
      inventory(c)=.5_dp*dt*sum(rate(c,:,:,:))*product(h)
    enddo
    call require(maxval(abs(inventory(:5)+after-before))<1e-18_dp,'ETAB old source compensates child momentum energy')
    call require(maxval(abs(inventory(2:5)-compensation))<1e-18_dp.and.inventory(5)<0, &
      'ETAB nonzero tangential energy feedback')
    old_rate=rate;gas=0;temp=0
    rate(:,2,2,2)=0
    call require(all(rate==0),'ETAB compensation stays in parent containing cell')
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,1,2, &
      faces,drops,masses,work,rate,wall,outside,forced,ok,last_id=watermark)
    call require(ok.and.all(rate==old_rate).and.watermark==37,'ETAB middle caches compensated source')
    gas=initial;temp=it
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,2,2, &
      faces,drops,masses,work,rate,wall,outside,forced,ok,last_id=watermark)
    call require(ok.and.size(drops)==count.and.watermark==37+count,'ETAB final commits reserved IDs without births')
    after=0
    do p=1,count
      call require(all(drops(p)%position==work%midpoint(p)%position),'ETAB final has no drift')
      after=after+parcel_totals(drops(p),liquids(1))
    enddo
    do c=1,12
      inventory(c)=.5_dp*dt*sum(old_rate(c,:,:,:)+rate(c,:,:,:))*product(h)
    enddo
    call require(maxval(abs(inventory(:5)+after-before))<1e-18_dp,'ETAB complete stage conservation')
    call require(abs(sum(masses(1,:)*drops%multiplicity)-original%mass*original%multiplicity)<1e-25_dp, &
      'ETAB component inventory conservation')
    deallocate(drops,masses);allocate(drops(1),masses(1,1))
    drops(1)=original;masses=original%mass;gas=initial;temp=it;work=spray_sdc_stage_3d();watermark=37
    options%maximum_parcels=1;rate=41
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,0,2, &
      faces,drops,masses,work,rate,wall,outside,forced,ok,last_id=watermark)
    call require(.not.ok.and.work%next_stage==0.and.all(rate==41).and.watermark==37.and.size(drops)==1, &
      'ETAB capacity failure rollback')
    options%maximum_parcels=1000000;watermark=huge(watermark)
    call advance_spray_sdc_stage_3d(species,transport,liquids,options,gas,temp,[0._dp,0._dp,0._dp],h,0._dp,dt,0,2, &
      faces,drops,masses,work,rate,wall,outside,forced,ok,last_id=watermark)
    call require(.not.ok.and.work%next_stage==0.and.all(rate==41).and.watermark==huge(watermark), &
      'ETAB ID overflow rollback')
    watermark=37
    call advance_spray_reactive_sdc_3d(species,reactions,transport,liquids,options,gas,temp,drops,masses, &
      [0._dp,0._dp,0._dp],h,dt,0._dp,2,.false.,[.false.,.false.,.false.,.false.], &
      [0._dp,0._dp,0._dp],0._dp,.7_dp,.7_dp,1e-6_dp,1e-14_dp,faces, &
      boundary,source,exchange,wall,outside,forced,ok,last_id=watermark)
    call require(ok.and.size(drops)==count.and.watermark==37+count,'ETAB successful gas particle coupling')
    after=0
    do p=1,size(drops)
      after=after+parcel_totals(drops(p),liquids(1))
    enddo
    do c=1,12
      inventory(c)=sum(gas(c,:,:,:)-initial(c,:,:,:))*product(h)
    enddo
    call require(maxval(abs(inventory-boundary-source-exchange))<1e-14_dp,'ETAB coupled gas exchange ledger')
    call require(maxval(abs(exchange(:5)+after-before))<1e-18_dp,'ETAB coupled total inventory balance')
    deallocate(drops,masses);allocate(drops(1),masses(1,1))
    drops(1)=original;masses=original%mass;gas=initial;temp=it
    watermark=37;boundary=41;source=42;exchange=43;wall=44;outside=45;forced=46;reaction_calls=0
    call advance_spray_reactive_sdc_3d(species,reactions,transport,liquids,options,gas,temp,drops,masses, &
      [0._dp,0._dp,0._dp],h,dt,0._dp,2,.true.,[.false.,.false.,.false.,.false.], &
      [0._dp,0._dp,0._dp],0._dp,.7_dp,.7_dp,1e-6_dp,1e-14_dp,faces, &
      boundary,source,exchange,wall,outside,forced,ok,evaluate_reaction=late_reaction_failure,last_id=watermark)
    call require(.not.ok.and.reaction_calls==65,'ETAB coupled failure after final private replacement')
    call require(size(drops)==1.and.watermark==37.and.all(drops(1)%velocity==original%velocity).and. &
      all(masses==original%mass).and.all(gas==initial).and.all(temp==it),'ETAB coupled state and ID rollback')
    call require(all(boundary==41).and.all(source==42).and.all(exchange==43).and. &
      all(wall==44).and.all(outside==45).and.all(forced==46),'ETAB coupled ledger rollback')
  end subroutine check_etab_replacement
  subroutine reference_gas_at(position,primitive,temperature,field,field_temperature)
    real(dp), intent(in) :: position(3)
    real(dp), intent(out) :: primitive(12),temperature
    real(dp), optional, intent(in) :: field(:,:,:,:),field_temperature(:,:,:)
    real(dp) :: weights(8),cellq(12),cellt,cs
    integer :: indices(3,8),cell
    call particle_boundary_stencil(position,[0._dp,0._dp,0._dp],h,[4,4,4],faces,indices,weights,ok)
    call require(ok,'reference trilinear stencil')
    primitive=0;temperature=0
    do cell=1,8
      if(weights(cell)<=0) cycle
      if(present(field).and.present(field_temperature)) then
        call reactive_conserved_to_primitive(species,field(:,indices(1,cell),indices(2,cell),indices(3,cell)), &
          field_temperature(indices(1,cell),indices(2,cell),indices(3,cell)),cellq,cellt,cs,ok)
      else
        call reactive_conserved_to_primitive(species,initial(:,indices(1,cell),indices(2,cell),indices(3,cell)), &
          it(indices(1,cell),indices(2,cell),indices(3,cell)),cellq,cellt,cs,ok)
      endif
      call require(ok,'reference cell primitive')
      primitive=primitive+weights(cell)*cellq;temperature=temperature+weights(cell)*cellt
    enddo
    primitive(5)=mixture_pressure(species,primitive(6:),primitive(1),temperature,ok)
    call require(ok,'reference sampled primitive')
  end subroutine reference_gas_at
  subroutine late_reaction_failure(records,rxns,energy,forcing,energy_forcing,duration,rtol,atol, &
      densities,temperature,steps,good,maximum_steps)
    type(nasa7_species), intent(in) :: records(:)
    type(elementary_reaction), intent(in) :: rxns(:)
    real(dp), intent(in) :: energy,forcing(:),energy_forcing,duration,rtol,atol
    real(dp), intent(inout) :: densities(:),temperature
    integer, intent(out) :: steps
    logical, intent(out) :: good
    integer, optional, intent(in) :: maximum_steps
    reaction_calls=reaction_calls+1;good=.false.;steps=0
    if(reaction_calls==65) return
    call advance_forced_reactor_interval(records,rxns,energy,forcing,energy_forcing,duration,rtol,atol, &
      densities,temperature,steps,good,maximum_steps=maximum_steps)
  end subroutine
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if(.not.condition) then
      print *, 'FAIL: ',message
      error stop 1
    endif
  end subroutine
end program
