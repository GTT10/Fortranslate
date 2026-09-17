program test_spray_eb_impact_coupling
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
  use spray_eb_step_mod, only: advance_eb_mixture_coupling_allocated
  use spray_eb_impact_mod, only: spray_eb_wall_options
  use spray_wall_impact_mod, only: spray_wall_film
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  type(eb_geometry_3d) :: g
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(spray_liquid) :: liquids(2)
  type(spray_options) :: options
  type(spray_parcel), allocatable :: drops(:)
  type(spray_parcel) :: saved_drop
  type(spray_eb_wall_options) :: walls
  type(spray_wall_film), allocatable :: films(:)
  real(dp), allocatable :: masses(:,:),film_mass(:,:)
  integer(int64) :: last_id=1
  real(dp) :: phi(0:4,0:4,0:4),t,sound,r,temperature(4,4,4),saved_t(4,4,4),saved_mass(2,1)
  real(dp) :: ledger(5),interval=2e-5_dp
  real(dp), allocatable :: primitive(:),cell(:),y(:),state(:,:,:,:),saved(:,:,:,:),before(:),after(:),expected(:)
  integer :: i,j,k,nv,nitrogen,steps
  logical :: ok
  do i=0,4
    phi(i,:,:)=real(i,dp)*.001_dp-.00125_dp
  end do
  call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[.004_dp,.004_dp,.004_dp],g,ok)
  call require(ok,'geometry')
  call load_h2o2_full_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_full_transport(transport,ok)
  call require(ok,'transport')
  call initialize_spray_liquid(species,'H2O',1000._dp,4184._dp,373.15_dp,101325._dp,2.257e6_dp,liquids(1),ok)
  call require(ok,'water')
  call initialize_spray_liquid(species,'H2O2',800._dp,2000._dp,350._dp,101325._dp,1e6_dp,liquids(2),ok)
  call require(ok,'peroxide')
  nv=reactive_nvar(size(species))
  allocate(primitive(nv),cell(nv),y(size(species)),state(nv,4,4,4),before(nv),after(nv),expected(nv))
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
  do k=1,4
  do j=1,4
  do i=1,4
    state(:,i,j,k)=cell
  end do
  end do
  end do
  state(:,1,:,:)=0;temperature=t;temperature(1,:,:)=42
  allocate(drops(1),masses(2,1),films(0),film_mass(2,0))
  options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp;walls%maximum_films=8
  drops(1)%id=1;drops(1)%position=[.001251_dp,.0015_dp,.0015_dp]
  drops(1)%velocity=[-20._dp,10._dp,5._dp];drops(1)%mass=1e-12_dp
  drops(1)%temperature=350;drops(1)%multiplicity=100
  masses(:,1)=drops(1)%mass*[.25_dp,.75_dp]
  saved=state;saved_t=temperature;saved_drop=drops(1);saved_mass=masses
  ledger=42;options%maximum_substeps=1
  call advance()
  call require(.not.ok.and.steps==0,'whole-interval failure')
  call require(all(state==saved).and.all(temperature==saved_t).and.all(masses==saved_mass),'gas and composition rollback')
  call require(all(drops(1)%position==saved_drop%position).and.all(ledger==42),'particle and wall ledger rollback')
  call require(size(films)==0.and.last_id==1,'film and watermark rollback')
  options%maximum_substeps=100000
  call inventory(before)
  call advance()
  call require(ok.and.steps>1,'coupled substeps')
  call inventory(after)
  expected=0;expected(1:5)=ledger
  call require(all(abs(after-before-expected)<=3e-12_dp*max(abs(before),abs(after),1e-30_dp)), &
    'gas plus particles minus wall ledger conservation')
  call require(size(drops)==0.and.size(films)==1.and.last_id==1,'deposition and unchanged ID watermark')
  call require(films(1)%temperature>350.and.films(1)%mass<1e-10_dp,'heating and evaporation before impact')
  call require(ledger(2)>0,'deposition wall impulse')
  call require(all(state(:,1,:,:)==0).and.all(temperature(1,:,:)==42),'covered cells unchanged')
  state=saved;temperature=saved_t
  deallocate(drops,masses,films,film_mass)
  allocate(drops(1),masses(2,1),films(0),film_mass(2,0))
  drops(1)=saved_drop;drops(1)%velocity(1)=-200;masses=saved_mass
  last_id=1;interval=2e-6_dp
  call inventory(before)
  call advance()
  call require(ok.and.size(drops)==4.and.last_id==5,'coupled splash population')
  call require(steps>1,'children participate in later source substeps')
  call inventory(after)
  expected=0;expected(1:5)=ledger
  call require(all(abs(after-before-expected)<=3e-12_dp*max(abs(before),abs(after),1e-30_dp)), &
    'splash gas film wall conservation')
  print *, 'PASS: EB mixture coupling with heat, evaporation, drag, deposition and splash'
contains
  subroutine advance()
    call advance_eb_mixture_coupling_allocated(species,transport,liquids,options,g,phi,state,temperature, &
      drops,masses,interval,steps,ledger,ok,walls,films,film_mass,last_id)
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
