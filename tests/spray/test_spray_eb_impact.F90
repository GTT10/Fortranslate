program test_spray_eb_impact
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use nasa7_thermo_mod, only: nasa7_species
  use spray_parcel_mod
  use spray_mixture_mod, only: liquid_mixture_properties
  use spray_wall_impact_mod, only: spray_wall_film,impact_deposit,impact_rebound,impact_dry
  use spray_eb_impact_mod, only: impact_mixture_level_set,advance_mixture_level_set_impacts
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(spray_liquid) :: liquids(2)
  type(spray_options) :: options
  type(spray_parcel) :: parent,original
  type(spray_parcel), allocatable :: children(:)
  type(spray_parcel), allocatable :: population(:)
  type(spray_wall_film), allocatable :: patches(:)
  real(dp), allocatable :: population_mass(:,:),patch_mass(:,:)
  integer(int64) :: watermark
  real(dp) :: exterior(7)
  integer :: hits
  type(spray_wall_film) :: film
  real(dp), allocatable :: child_mass(:,:)
  real(dp) :: phi(0:2,0:2,0:2),masses(2),film_mass(2),wall(5),fraction,normal(3),previous(3)
  real(dp) :: y(2),mw(2),rho,cp,e,amount,before(5),after(5),inventory(2),temperature,speed
  integer :: i,j,k,mode,p,regime
  integer(int64) :: first_id=2_int64
  character(len=16) :: faces(6)
  logical :: hit,ok
  call load_h2o2_full_thermo(species,ok)
  call require(ok,'thermo')
  call initialize_spray_liquid(species,'H2O',1000._dp,4184._dp,373.15_dp,101325._dp,2.257e6_dp,liquids(1),ok)
  call require(ok,'water')
  call initialize_spray_liquid(species,'H2O2',800._dp,2000._dp,350._dp,101325._dp,1e6_dp,liquids(2),ok)
  call require(ok,'peroxide')
  options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
  mw=[(species(liquids(i)%vapor_index)%molecular_weight,i=1,2)]
  do k=0,2
    do j=0,2
      do i=0,2
        phi(i,j,k)=.5_dp*(i+j)-1
      end do
    end do
  end do
  normal=[1._dp,1._dp,0._dp]/sqrt(2._dp);previous=[.6_dp,.6_dp,.5_dp];y=[.6_dp,.4_dp]
  call liquid_mixture_properties(liquids,y,350._dp,rho,cp,e,ok)
  call require(ok,'mixture')
  do mode=1,3
    speed=1;temperature=300
    if(mode==2) temperature=600
    if(mode==3) speed=200
    parent=spray_parcel(1_int64,[.4_dp,.4_dp,.5_dp],-speed*normal,1e-12_dp,350._dp,100._dp)
    original=parent;masses=parent%mass*y;wall=42;fraction=43;regime=44
    call impact()
    call require(ok.and.hit,'oblique contact')
    call require(abs(fraction-.5_dp)<1e-14_dp,'first contact fraction')
    if(mode==1) call require(regime==impact_deposit,'deposit')
    if(mode==2) call require(regime==impact_rebound,'rebound')
    if(mode==3) call require(regime==impact_dry.and.size(children)==4,'dry splash')
    amount=original%mass*original%multiplicity
    before=[amount,amount*original%velocity,amount*(e+.5_dp*sum(original%velocity**2))]
    amount=parent%mass*parent%multiplicity
    after=[amount+film%mass,amount*parent%velocity,amount*(e+.5_dp*sum(parent%velocity**2))+film%mass*e]
    inventory=parent%multiplicity*masses+film_mass
    do p=1,size(children)
      amount=children(p)%mass*children(p)%multiplicity
      after=after+[amount,amount*children(p)%velocity,amount*(e+.5_dp*sum(children(p)%velocity**2))]
      inventory=inventory+children(p)%multiplicity*child_mass(:,p)
      call require(dot_product(children(p)%velocity,normal)>=0,'outgoing splash')
    end do
    call require(all(abs(after+wall-before)<=1e-12_dp*max(abs(before),1e-20_dp)),'mass momentum energy budget')
    call require(all(abs(inventory-original%mass*original%multiplicity*y)<1e-23_dp),'component budget')
    if(film%mass>0) call require(maxval(abs(film%normal-normal))<1e-14_dp,'film normal')
    allocate(population(1),population_mass(2,1),patches(0),patch_mass(2,0))
    population(1)=original;population(1)%position=previous;population_mass(:,1)=original%mass*y
    watermark=1;wall=42;hits=43
    ! Fail after contact construction, then retry from the unchanged population.
    call advance_interval(1,0)
    if(mode/=2) then
      call require(.not.ok.and.watermark==1.and.hits==43.and.all(wall==42),'capacity rollback')
      call require(size(population)==1.and.size(patches)==0,'population rollback')
      call require(all(population(1)%position==previous),'position rollback')
    else
      call require(ok.and.hits==1,'rebound interval')
      call require(sum(population(1)%position(:2))>1,'rebound remaining flight')
      population(1)=original;population(1)%position=previous;watermark=1
    end if
    call advance_interval(5,1)
    call require(ok.and.hits==1,'impact interval')
    if(mode==1) call require(size(population)==0.and.size(patches)==1,'deposited population')
    if(mode==3) then
      call require(size(population)==4.and.watermark==5,'splash IDs')
      call require(norm2(population(1)%position-children(1)%position)>1e-3_dp,'splash remaining flight')
    end if
    after=0;inventory=0
    do p=1,size(population)
      amount=population(p)%mass*population(p)%multiplicity
      after=after+[amount,amount*population(p)%velocity,amount*(e+.5_dp*sum(population(p)%velocity**2))]
      inventory=inventory+population(p)%multiplicity*population_mass(:,p)
    end do
    do p=1,size(patches)
      after(1)=after(1)+patches(p)%mass;after(5)=after(5)+patches(p)%mass*e
      inventory=inventory+patch_mass(:,p)
    end do
    call require(all(abs(after+wall-before)<=1e-12_dp*max(abs(before),1e-20_dp)),'interval budget')
    call require(all(abs(inventory-original%mass*original%multiplicity*y)<1e-23_dp),'interval components')
    deallocate(population,population_mass,patches,patch_mass)
  end do
  parent=original;masses=parent%mass*y;parent%position=[.7_dp,.7_dp,.5_dp]
  call impact()
  call require(ok.and..not.hit.and.fraction==1.and.all(wall==0),'fluid miss')
  parent=original;masses=parent%mass*y;wall=42;fraction=43;regime=44
  temperature=-1
  call impact()
  call require(.not.ok.and..not.hit.and.all(wall==42).and.fraction==43.and.regime==44,'diagnostic rollback')
  call require(parent%mass==original%mass.and.all(masses==original%mass*y),'parent rollback')
  temperature=300;first_id=1
  call impact()
  call require(.not.ok.and..not.hit.and.all(wall==42).and.fraction==43.and.regime==44,'splash ID rollback')
  call require(parent%mass==original%mass.and.all(masses==original%mass*y),'late parent rollback')
  call require(size(children)==4.and.children(1)%id==2_int64,'child buffer rollback')
  do i=0,2
    phi(i,:,:)=.1_dp-abs(.5_dp*i-.5_dp)
  end do
  allocate(population(1),population_mass(2,1),patches(0),patch_mass(2,0))
  population(1)=original;population(1)%position=[.55_dp,.5_dp,.5_dp]
  population(1)%velocity=[1._dp,0._dp,0._dp];population_mass(:,1)=original%mass*y
  watermark=1;wall=42;hits=43
  do mode=1,2
    call advance_mixture_level_set_impacts(phi,[0._dp,0._dp,0._dp],[1._dp,1._dp,1._dp], &
      .4_dp,mode,1,0,liquids,options,mw,101325._dp,600._dp,acos(-1._dp)/3, &
      population,population_mass,patches,patch_mass,watermark,wall,hits,ok)
    if(mode==1) then
      call require(.not.ok.and.hits==43.and.all(wall==42),'second contact rollback')
      call require(population(1)%position(1)==.55_dp,'first contact undone')
    else
      call require(ok.and.hits==2,'two contacts')
      call require(population(1)%position(1)>.4_dp.and.population(1)%position(1)<.6_dp,'two-wall endpoint')
    end if
  end do
  ! Splash children may leave the root domain during their remaining flight.
  do k=0,2
  do j=0,2
  do i=0,2
    phi(i,j,k)=.5_dp*(i+j)-1
  end do
  end do
  end do
  deallocate(population,population_mass,patches,patch_mass)
  allocate(population(1),population_mass(2,1),patches(0),patch_mass(2,0))
  population(1)=original;population(1)%position=previous
  population_mass(:,1)=original%mass*y;watermark=1;wall=42;exterior=43;hits=17
  call escaping_splash(1)
  call require(.not.ok.and.watermark==1.and.hits==17,'splash outflow depth rollback')
  call require(size(population)==1.and.size(patches)==0.and.all(wall==42).and.all(exterior==43), &
    'splash outflow allocation ledger rollback')
  call require(all(population(1)%position==previous).and.all(population_mass(:,1)==original%mass*y), &
    'splash outflow particle rollback')
  amount=original%mass*original%multiplicity
  before=[amount,amount*original%velocity,amount*(e+.5_dp*sum(original%velocity**2))]
  call escaping_splash(4)
  call require(ok.and.size(population)==0.and.size(patches)==1.and.watermark==5,'splash children escape')
  call require(hits==5,'impact and four exterior events counted')
  after=0;inventory=0
  do p=1,size(patches)
    after(1)=after(1)+patches(p)%mass;after(5)=after(5)+patches(p)%mass*e
    inventory=inventory+patch_mass(:,p)
  end do
  call require(all(abs(after+wall-before-exterior(:5))<=1e-12_dp*max(abs(before),1e-20_dp)), &
    'splash exterior mass momentum energy conservation')
  call require(all(abs(inventory-original%mass*original%multiplicity*y-exterior(6:))<1e-23_dp), &
    'splash exterior component conservation')
  ! Outer crossing precedes deposition; preserve the remaining flight time.
  do i=0,2
    phi(i,:,:)=.5_dp*i-.5_dp
  end do
  do mode=1,2
    deallocate(population,population_mass,patches,patch_mass)
    allocate(population(1),population_mass(2,1),patches(0),patch_mass(2,0))
    population(1)=original;population(1)%position=[.6_dp,.5_dp,.99_dp]
    population(1)%velocity=[-1._dp,0._dp,1._dp]
    population_mass(:,1)=original%mass*y;watermark=1;wall=42;exterior=43;hits=17
    faces='outflow';faces(5:6)='periodic'
    if(mode==2) faces(5:6)='no_slip_wall'
    amount=original%mass*original%multiplicity
    before=[amount,amount*population(1)%velocity,amount*(e+.5_dp*sum(population(1)%velocity**2))]
    call crossing_deposit(1)
    call require(.not.ok.and.hits==17.and.all(wall==42).and.all(exterior==43),'crossing deposit rollback')
    call require(size(population)==1.and.size(patches)==0.and.watermark==1,'crossing population rollback')
    call crossing_deposit(4)
    call require(ok.and.hits==2.and.size(population)==0.and.size(patches)==1,'crossing then deposition')
    call require(abs(patches(1)%position(3)-merge(.09_dp,.91_dp,mode==1))<1e-14_dp,'remaining flight after crossing')
    after=0;after(1)=patches(1)%mass;after(5)=patches(1)%mass*e
    call require(all(abs(after+wall-before-exterior(:5))<1e-12_dp*max(abs(before),1e-20_dp)), &
      'crossing deposition budget')
    if(mode==1) call require(all(exterior==0),'periodic crossing has no exterior transfer')
  end do
  print *, 'PASS: EB impact remaining flight and exterior loss'
contains
  subroutine crossing_deposit(limit)
    integer, intent(in) :: limit
    call advance_mixture_level_set_impacts(phi,[0._dp,0._dp,0._dp],[1._dp,1._dp,1._dp], &
      .2_dp,limit,5,1,liquids,options,mw,101325._dp,300._dp,acos(-1._dp)/3, &
      population,population_mass,patches,patch_mass,watermark,wall,hits,ok,outflow_change=exterior,outer_faces=faces)
  end subroutine

  subroutine escaping_splash(limit)
    integer, intent(in) :: limit
    call advance_mixture_level_set_impacts(phi,[0._dp,0._dp,0._dp],[1._dp,1._dp,1._dp], &
      .1_dp,limit,5,1,liquids,options,mw,101325._dp,300._dp,acos(-1._dp)/3, &
      population,population_mass,patches,patch_mass,watermark,wall,hits,ok,outflow_change=exterior)
  end subroutine
  subroutine advance_interval(capacity,film_capacity)
    integer, intent(in) :: capacity,film_capacity
    call advance_mixture_level_set_impacts(phi,[0._dp,0._dp,0._dp],[1._dp,1._dp,1._dp], &
      .2_dp*sqrt(2._dp)/speed,4,capacity,film_capacity,liquids,options,mw,101325._dp,temperature, &
      acos(-1._dp)/3,population,population_mass,patches,patch_mass,watermark,wall,hits,ok)
  end subroutine
  subroutine impact()
    call impact_mixture_level_set(phi,[0._dp,0._dp,0._dp],[1._dp,1._dp,1._dp],previous,parent, &
      liquids,masses,options,mw,101325._dp,temperature,0._dp,acos(-1._dp)/3,first_id, &
      children,child_mass,film,film_mass,wall,fraction,regime,hit,ok)
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
