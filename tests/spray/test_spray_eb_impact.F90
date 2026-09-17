program test_spray_eb_impact
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use nasa7_thermo_mod, only: nasa7_species
  use spray_parcel_mod
  use spray_mixture_mod, only: liquid_mixture_properties
  use spray_wall_impact_mod, only: spray_wall_film,impact_deposit,impact_rebound,impact_dry
  use spray_eb_impact_mod, only: impact_mixture_level_set
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(spray_liquid) :: liquids(2)
  type(spray_options) :: options
  type(spray_parcel) :: parent,original
  type(spray_parcel), allocatable :: children(:)
  type(spray_wall_film) :: film
  real(dp), allocatable :: child_mass(:,:)
  real(dp) :: phi(0:2,0:2,0:2),masses(2),film_mass(2),wall(5),fraction,normal(3),previous(3)
  real(dp) :: y(2),mw(2),rho,cp,e,amount,before(5),after(5),inventory(2),temperature,speed
  integer :: i,j,k,mode,p,regime
  integer(int64) :: first_id=2_int64
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
  print *, 'PASS: oblique EB mixture impact regimes and conservation'
contains
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
