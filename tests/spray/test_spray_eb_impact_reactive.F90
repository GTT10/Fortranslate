program test_spray_eb_impact_reactive
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_full_transport
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d
  use spray_parcel_mod
  use spray_mixture_mod, only: liquid_mixture_properties
  use spray_eb_reactive_mod, only: advance_eb_spray_reactive_allocated
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use spray_eb_impact_mod, only: spray_eb_wall_options
  use spray_wall_impact_mod, only: spray_wall_film
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  type(eb_geometry_3d) :: g
  type(elementary_reaction), allocatable :: reactions(:)
  character(len=16) :: solver="rusanov"
  real(dp), allocatable :: boundary(:),wall(:),chemical(:)
  real(dp) :: theta
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(spray_liquid) :: liquids(2)
  type(spray_options) :: options
  type(spray_parcel), allocatable :: drops(:)
  type(spray_parcel) :: saved_drop
  type(spray_wall_film), allocatable :: films(:)
  type(spray_eb_wall_options) :: walls
  real(dp), allocatable :: masses(:,:),film_mass(:,:)
  integer(int64) :: last_id,split_created
  integer :: mode
  real(dp) :: phi(0:4,0:4,0:4),t,sound,r,temperature(4,4,4),saved_t(4,4,4),saved_mass(2,1)
  real(dp) :: ledger(5),interval=1e-8_dp
  real(dp), allocatable :: primitive(:),cell(:),y(:),state(:,:,:,:),saved(:,:,:,:),before(:),after(:),expected(:)
  integer :: i,j,k,nv,nitrogen,steps
  logical :: ok,rejected
  do i=0,4
    phi(i,:,:)=real(i,dp)*.001_dp-.00125_dp
  end do
  call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[.004_dp,.004_dp,.004_dp],g,ok)
  call require(ok,'geometry')
  call load_h2o2_full_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_full_mechanism(reactions,ok)
  call require(ok,'mechanism')
  call load_h2o2_full_transport(transport,ok)
  call require(ok,'transport')
  call initialize_spray_liquid(species,'H2O',1000._dp,4184._dp,373.15_dp,101325._dp,2.257e6_dp,liquids(1),ok)
  call require(ok,'water')
  call initialize_spray_liquid(species,'H2O2',800._dp,2000._dp,350._dp,101325._dp,1e6_dp,liquids(2),ok)
  call require(ok,'peroxide')
  nv=reactive_nvar(size(species))
  allocate(primitive(nv),cell(nv),y(size(species)),state(nv,4,4,4), &
    before(nv),after(nv),expected(nv),boundary(nv),wall(nv),chemical(nv))
  nitrogen=0
  do i=1,size(species)
    if(species(i)%name=='N2') nitrogen=i
  end do
  call require(nitrogen>0,'nitrogen identity')
  y=0;y(nitrogen)=1;r=mixture_specific_gas_constant(species,y,ok)
  call require(ok,'gas constant')
  primitive=0;primitive(1)=1;primitive(2:4)=[1._dp,2._dp,3._dp];primitive(5)=r*1000;primitive(6:)=y
  call reactive_primitive_to_conserved(species,primitive,cell,t,sound,ok)
  call require(ok,'initial gas')
  options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
  walls%evolve=.true.;walls%maximum_films=8
  do mode=1,2
  allocate(drops(1),masses(2,1),films(0),film_mass(2,0));last_id=1
  do k=1,4
  do j=1,4
  do i=1,4
    state(:,i,j,k)=cell
  end do
  end do
  end do
  state(:,1,:,:)=0;temperature=t;temperature(1,:,:)=42
  drops(1)%id=1;drops(1)%position=[.001250001_dp,.0015_dp,.0015_dp]
  drops(1)%velocity=[-20._dp,10._dp,5._dp];drops(1)%mass=1e-12_dp
  if(mode==2) drops(1)%velocity(1)=-200
  drops(1)%temperature=350;drops(1)%multiplicity=100
  masses(:,1)=drops(1)%mass*[.25_dp,.75_dp]
  saved=state;saved_t=temperature;saved_drop=drops(1);saved_mass=masses
  ledger=42;boundary=42;wall=42;chemical=42;theta=42;solver="invalid"
  call advance()
  call require(.not.ok.and.steps==0,'gas failure after first particle half')
  call require(all(state==saved).and.all(temperature==saved_t).and.all(masses==saved_mass),'gas and composition rollback')
  call require(all(drops(1)%position==saved_drop%position).and.all(boundary==42).and.all(wall==42).and. &
    all(chemical==42).and.theta==42,'particle and wall ledger rollback')
  call require(size(films)==0.and.last_id==1,'film and ID rollback after gas failure')
  solver="rusanov"
  if(mode==1) then
    drops(1)%position(1)=.00125015_dp
    walls%maximum_films=0;options%maximum_substeps=8
    call advance()
    call require(.not.ok.and.steps==0,'second particle half failure')
    call require(all(state==saved).and.all(temperature==saved_t),'gas step rollback after second half')
    call require(drops(1)%position(1)==.00125015_dp.and.all(masses==saved_mass),'second half parcel rollback')
    call require(size(films)==0.and.last_id==1.and.all(wall==42).and.theta==42,'second half diagnostics rollback')
    drops(1)=saved_drop;walls%maximum_films=8;options%maximum_substeps=100000
  end if
  call inventory(before)
  call advance()
  call require(ok.and.steps>1,'coupled substeps')
  call inventory(after)
  expected=boundary+wall+chemical
  call require(all(abs(after-before-expected)<=3e-12_dp*max(abs(before),abs(after),1e-30_dp)), &
    'gas plus particles minus wall ledger conservation')
  if(mode==1) then
    call require(size(drops)==0.and.size(films)==1.and.last_id==1,'reactive deposition')
    call require(films(1)%temperature/=350,'film participates in second source half')
  else
    call require(size(drops)==4.and.last_id==5,'reactive splash IDs')
    call require(all(drops%temperature>350),'children heating')
  end if
  call require(wall(2)>0,'wall impulse')
  call require(all(state(:,1,:,:)==0).and.all(temperature(1,:,:)==42),'covered cells unchanged')
  call require(theta>=0.and.theta<=1,'transport limiter')
  call require(any(chemical(6:)/=0),'active chemistry ledger')
  deallocate(drops,masses,films,film_mass)
  end do
  state=saved;temperature=saved_t
  allocate(drops(1),masses(2,1),films(1),film_mass(2,1));last_id=1
  drops(1)=saved_drop;masses=saved_mass
  films(1)%position=[.00125_dp,.0015_dp,.0015_dp];films(1)%normal=[1._dp,0._dp,0._dp]
  films(1)%mass=1e-10_dp;films(1)%temperature=350;films(1)%height=2e-6_dp;films(1)%diameter=1e-4_dp
  film_mass(:,1)=films(1)%mass*[.25_dp,.75_dp]
  walls%evolve=.false.;wall=42;theta=42
  call advance()
  call require(.not.ok.and.rejected.and.steps==0,'inadmissible wet correlation is terminal')
  call require(all(state==saved).and.all(temperature==saved_t).and.all(masses==saved_mass),'model gas rollback')
  call require(films(1)%mass==1e-10_dp.and.last_id==1.and.all(wall==42).and.theta==42,'model film ID rollback')
  ! Split only after the complete reactive update; retired IDs must never recur.
  deallocate(drops,masses,films,film_mass)
  allocate(drops(1),masses(2,1),films(0),film_mass(2,0))
  state=saved;temperature=saved_t;drops(1)=saved_drop;masses=saved_mass
  drops(1)%position=[.002_dp,.0015_dp,.0015_dp];drops(1)%velocity=0
  saved_drop=drops(1)
  options%breakup_model='etab';options%maximum_multiplicity=20
  options%surface_tension=1e3_dp;options%maximum_parcels=1
  options%drag=.false.;options%heat=.false.;options%evaporation=.false.
  last_id=100;boundary=42;wall=42;chemical=42;theta=42
  call advance()
  call require(.not.ok.and.steps==0.and.last_id==100.and.split_created==0,'late split capacity rollback')
  call require(size(drops)==1.and.all(masses==saved_mass),'split population rollback')
  call require(all(drops(1)%position==saved_drop%position),'split position rollback')
  call require(all(state==saved).and.all(temperature==saved_t),'split gas rollback')
  call require(all(boundary==42).and.all(wall==42).and.all(chemical==42).and.theta==42,'split ledger rollback')
  options%maximum_parcels=8
  do mode=1,2
    state=saved;temperature=saved_t
    deallocate(drops,masses);allocate(drops(1),masses(2,1))
    drops(1)=saved_drop;masses=saved_mass;last_id=100
    call inventory(before)
    if(mode==1) then
      call advance()
    else
      call advance_eb_spray_reactive_allocated(species,reactions,transport,liquids,options,g,phi,state,temperature, &
        drops,masses,interval,solver,'state_redist',.true.,.true.,.true.,.true.,.true.,.true., &
        2e-7_dp,1e-12_dp,steps,theta,boundary,wall,chemical,ok,chemistry_integrator='implicit', &
        last_id=last_id,split_created=split_created)
    end if
    call require(ok.and.size(drops)==5.and.last_id==104.and.split_created==4,'split with retired IDs')
    call require(all(drops%id==[1_int64,101_int64,102_int64,103_int64,104_int64]),'split ID sequence')
    call require(all(drops%multiplicity==20).and.all(drops(2:)%parent_id==1),'split weights and lineage')
    do i=1,size(drops)
      call require(all(masses(:,i)==saved_mass(:,1)),'split component copying')
    end do
    call inventory(after)
    call require(all(abs(after-before-boundary-wall-chemical)<=3e-12_dp*max(abs(before),abs(after),1e-30_dp)), &
      'split composite conservation')
  end do
  state=saved;temperature=saved_t
  deallocate(drops,masses);allocate(drops(1),masses(2,1))
  drops(1)=saved_drop;drops(1)%velocity=[1000._dp,0._dp,0._dp];masses=saved_mass
  options%breakup_model='khrt';options%khrt_children=.true.;options%khrt_reference_diameter_cubed=1e-25_dp
  options%surface_tension=.072_dp;options%maximum_multiplicity=0;options%maximum_parcels=100
  last_id=100;boundary=42;wall=42;chemical=42;theta=42;solver='invalid'
  call kh_advance()
  call require(.not.ok.and.steps==0.and.last_id==100,'KH gas failure ID rollback')
  call require(size(drops)==1.and.all(masses==saved_mass),'KH population rollback')
  call require(all(state==saved).and.all(temperature==saved_t),'KH gas rollback')
  call require(all(boundary==42).and.all(wall==42).and.all(chemical==42).and.theta==42,'KH ledger rollback')
  solver='rusanov'
  call inventory(before)
  call kh_advance()
  call require(ok.and.size(drops)>1.and.last_id>100,'KH reactive children')
  call require(modulo(last_id-100_int64,4_int64)==0,'KH reactive four-child IDs')
  call inventory(after)
  call require(all(abs(after-before-boundary-wall-chemical)<=3e-12_dp*max(abs(before),abs(after),1e-30_dp)), &
    'KH reactive composite conservation')
  print *, 'PASS: reactive EB impacts, films, KH children, splitting and rollback'
contains
  subroutine kh_advance()
    call advance_eb_spray_reactive_allocated(species,reactions,transport,liquids,options,g,phi,state,temperature, &
      drops,masses,interval,solver,'state_redist',.true.,.true.,.true.,.true.,.true.,.true., &
      2e-7_dp,1e-12_dp,steps,theta,boundary,wall,chemical,ok,chemistry_integrator='implicit',last_id=last_id)
  end subroutine
  subroutine advance()
    call advance_eb_spray_reactive_allocated(species,reactions,transport,liquids,options,g,phi,state,temperature, &
      drops,masses,interval,solver,'state_redist',.true.,.true.,.true.,.true.,.true.,.true., &
      2e-7_dp,1e-12_dp,steps,theta,boundary,wall,chemical,ok,chemistry_integrator='implicit', &
      wall_options=walls,films=films,film_mass=film_mass,last_id=last_id,model_rejected=rejected,split_created=split_created)
  end subroutine
  subroutine inventory(total)
    real(dp), intent(out) :: total(:)
    real(dp) :: rho,cp,e,amount
    integer :: ii,jj,kk,n,p
    total=0
    do kk=1,4
    do jj=1,4
    do ii=1,4
      total=total+state(:,ii,jj,kk)*g%volume_fraction(ii,jj,kk)*g%dx*g%dy*g%dz
    end do
    end do
    end do
    do p=1,size(drops)
      call liquid_mixture_properties(liquids,masses(:,p)/drops(p)%mass,drops(p)%temperature,rho,cp,e,ok)
      call require(ok,'liquid inventory')
      amount=drops(p)%mass*drops(p)%multiplicity
      total(1)=total(1)+amount;total(2:4)=total(2:4)+amount*drops(p)%velocity
      total(5)=total(5)+amount*(e+.5_dp*sum(drops(p)%velocity**2))
      do n=1,2
        total(5+liquids(n)%vapor_index)=total(5+liquids(n)%vapor_index)+masses(n,p)*drops(p)%multiplicity
      end do
    end do
    do p=1,size(films)
      call liquid_mixture_properties(liquids,film_mass(:,p)/films(p)%mass,films(p)%temperature,rho,cp,e,ok)
      call require(ok,'film inventory')
      total(1)=total(1)+films(p)%mass;total(5)=total(5)+films(p)%mass*e
      do n=1,2
        total(5+liquids(n)%vapor_index)=total(5+liquids(n)%vapor_index)+film_mass(n,p)
      end do
    end do
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
