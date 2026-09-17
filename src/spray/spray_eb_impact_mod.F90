module spray_eb_impact_mod
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use spray_parcel_mod
  use spray_mixture_mod, only: valid_component_inventory,liquid_mixture_properties
  use spray_wall_impact_mod, only: spray_wall_film,resolve_mixture_wall_impact,mixture_impact_boiling_temperature
  use spray_eb_intersection_mod, only: intersect_spray_level_set
  implicit none
  private
  public :: impact_mixture_level_set
contains
  ! One earliest stationary EB contact; parent%position is the trial endpoint.
  ! Wall budget is received by the wall, as in resolve_mixture_wall_impact.
  ! Caller supplies local wall temperature/thickness and advances remaining flight.
  subroutine impact_mixture_level_set(phi,lower,upper,previous,parent,components,component_mass,properties, &
      molecular_weights,pressure,wall_temperature,film_thickness,contact_angle,first_id, &
      children,child_mass,deposit,film_mass,wall,fraction,regime,hit,ok)
    real(dp), intent(in) :: phi(0:,0:,0:),lower(3),upper(3),previous(3),molecular_weights(:)
    real(dp), intent(in) :: pressure,wall_temperature,film_thickness,contact_angle
    type(spray_parcel), intent(inout) :: parent
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(inout) :: component_mass(:),film_mass(:),wall(5),fraction
    type(spray_options), intent(in) :: properties
    integer(int64), intent(in) :: first_id
    type(spray_parcel), allocatable, intent(inout) :: children(:)
    real(dp), allocatable, intent(inout) :: child_mass(:,:)
    type(spray_wall_film), intent(inout) :: deposit
    integer, intent(inout) :: regime
    logical, intent(out) :: hit,ok
    type(spray_parcel) :: trial
    type(spray_parcel), allocatable :: spawned(:)
    type(spray_wall_film) :: deposited
    type(spray_liquid) :: effective
    real(dp), allocatable :: spawned_mass(:,:)
    real(dp) :: y(size(components)),masses(size(components)),film(size(components)),budget(5)
    real(dp) :: location(3),normal(3),at,rho,cp,e,boil,d,mu,un,we,re,kv,kcrit
    integer :: cell(3),kind
    logical :: crossed,valid
    ok=.false.;hit=.false.
    if(size(components)<1.or.size(components)>16) return
    if(size(component_mass)/=size(components).or.size(film_mass)/=size(components)) return
    if(.not.valid_spray_parcel(parent).or..not.valid_component_inventory(parent,component_mass)) return
    if(parent%mass<=0) return
    if(.not.all(ieee_is_finite([pressure,wall_temperature,film_thickness,contact_angle]))) return
    if(pressure<1.or.pressure>1e10_dp.or.wall_temperature<1.or.wall_temperature>1e5_dp.or.film_thickness<0) return
    at=0;location=0;normal=0;cell=0;crossed=.false.
    call intersect_spray_level_set(phi,lower,upper,previous,parent%position,crossed,at,location,normal,cell,valid)
    if(.not.valid) return
    if(.not.crossed) then
      wall=0;fraction=1;regime=-1;ok=.true.;return
    end if
    ! Geometry detects travel direction; reject a velocity inconsistent with inward impact.
    un=-dot_product(parent%velocity,normal)
    if(un<=0) return
    y=component_mass/parent%mass
    call liquid_mixture_properties(components,y,parent%temperature,rho,cp,e,valid)
    if(.not.valid) return
    call mixture_impact_boiling_temperature(components,y,molecular_weights,pressure,boil,valid)
    if(.not.valid) return
    effective=components(1);effective%density=rho;effective%density_model='constant'
    effective%density_coefficients=0;effective%density_temperature_range=0
    d=parcel_diameter(parent,effective);mu=liquid_viscosity_at_temperature(properties,parent%temperature)
    if(.not.ieee_is_finite(properties%surface_tension)) return
    if(d<=0.or.mu<=0.or.properties%surface_tension<1e-8_dp.or.properties%surface_tension>1e3_dp) return
    ! Same pinned impact parameter conversion as the regular planar dispatcher.
    we=rho*d*un**2/properties%surface_tension;re=un*d*rho/mu
    kv=sqrt(we*sqrt(re))
    kcrit=(26.788_dp*log(pressure*1e-5_dp)+33.546_dp)*exp(-.006_dp*(wall_temperature-260))+36
    if(kcrit<0) return
    trial=parent;masses=component_mass;film=0
    call resolve_mixture_wall_impact(trial,components,masses,location,normal,kv,wall_temperature/boil,kcrit, &
      film_thickness/d,we,re,contact_angle,first_id,spawned,spawned_mass,deposited,film,budget,kind,valid)
    if(.not.valid) return
    parent=trial;component_mass=masses;film_mass=film;deposit=deposited
    call move_alloc(spawned,children);call move_alloc(spawned_mass,child_mass)
    wall=budget;fraction=at;regime=kind;hit=.true.;ok=.true.
  end subroutine
end module
