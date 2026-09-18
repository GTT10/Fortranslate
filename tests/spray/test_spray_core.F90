program test_spray_core
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_full_transport
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use mixture_transport_mod, only: binary_diffusion_coefficients
  use spray_partition_mod
  use spray_mixture_mod
  use spray_parcel_mod
  use spray_wall_impact_mod
  use spray_breakup_mod
  use spray_coupling_3d_mod
  use les_filter_3d_mod
  use les_dynamic_smagorinsky_mod
  use les_smagorinsky_3d_mod
  implicit none
  integer(int64) :: split_watermark
  real(dp) :: port_centers(3,2),port_axes(3,2),port_fractions(2)
  real(dp) :: value,td,omega,wer,initial_mass,z,tkh,wavelength,tau,acceleration
  real(dp), parameter :: pi=acos(-1.0_dp)
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(spray_liquid) :: liquid,constant_liquid,mixture_components(2)
  real(dp) :: mixture_y(2),mixture_rho,mixture_cp,mixture_e,mixture_t,mixture_e0
  real(dp) :: component_mass(2),saved_mass(2),evaporated_mass(2),vapor_h(2),mw(2),surface_x(2),mixture_rates(2)
  real(dp), allocatable :: child_component_mass(:,:),saved_child_mass(:,:)
  real(dp) :: film_component_mass(2),saved_film_mass(2)
  type(spray_parcel) :: p,q
  type(spray_parcel), allocatable :: drops(:),saved_drops(:)
  type(spray_options) :: options,constant_options
  type(spray_khrt_shedding_state) :: shedding,saved_shedding
  type(spray_khrt_release) :: release
  type(spray_film) :: film
  real(dp) :: state(15,4,4,4),temp(4,4,4),saved(15,4,4,4),saved_t(4,4,4),rhs(15,4,4,4)
  real(dp) :: wall_t(6),wall_v(3,6)
  real(dp) :: totals0(15),totals1(15),spacing(3),lower(3),weights(8)
  real(dp) :: primitive(15),binary(10,10),skin_y(10)
  real(dp) :: inverse_tau,heat,masscoef,transfer(5),expected,dt,d,grad(3,3),stress(3,3),mu,diff,ke0,ke1
  real(dp) :: splash_mass,deposited_mass,splash_vel(5),azimuth_mass(4)
  type(spray_wall_film) :: wall_patch,saved_patch
  type(spray_wall_film) :: empty_films(0)
  logical :: hit
  type(spray_wall_options) :: wall_settings
  real(dp) :: child_totals(5),film_budget(5),film_mdot,film_qgas,film_qwall,film_tdot,solid_energy
  real(dp) :: impact_kv(6),impact_tstar(6),impact_film(6)
  integer :: impact_kinds(6),event_index
  integer, allocatable :: partition_counts(:),partition_offsets(:)
  integer :: les_case,jcell
  character(len=16) :: les_model
  integer :: regime
  integer :: indices(3,8),n,steps
  logical :: ok
  character(len=32) :: faces(6)
  character(len=64) :: test
  character(len=256) :: message
  call get_command_argument(1,test)
  call load_h2o2_full_thermo(species,ok)
  call require(ok,'load thermo')
  call load_h2o2_full_transport(transport,ok)
  call require(ok,'load transport')
  call initialize_spray_liquid(species,'H2O',1000.0_dp,4184.0_dp,373.15_dp,101325.0_dp,2.257e6_dp,liquid,ok)
  call require(ok,'initialize liquid')
  p=spray_parcel(1_int64,[0.0015_dp,0.0015_dp,0.0015_dp],[1.0_dp,0.0_dp,0.0_dp], &
    pi*1000.0_dp*(1.0e-4_dp)**3/6.0_dp,350.0_dp,100.0_dp)
  film=spray_film(1.0_dp,1.8e-5_dp,0.03_dp,2.5e-5_dp,1000.0_dp,500.0_dp, &
    [0.0_dp,0.0_dp,0.0_dp],0.2_dp,0.0_dp,liquid_energy(liquid,350.0_dp)+2.3e6_dp)
  options%abramzon_sirignano=.false.
  lower=0.0_dp; spacing=0.001_dp
  select case(trim(test))
  case('mixture_terminal')
    mixture_components(1)=liquid
    call initialize_spray_liquid(species,'H2O2',800.0_dp,2000.0_dp,350.0_dp,101325.0_dp,1e6_dp, &
      mixture_components(2),ok)
    call require(ok,'terminal mixture component')
    p%temperature=300;q=p;options%drag=.false.;options%heat=.false.
    primitive=0;primitive(1)=1;primitive(5)=101325;primitive(15)=1
    mixture_y=[1e-13_dp,1.0_dp-1e-13_dp];component_mass=p%mass*mixture_y;saved_mass=component_mass
    call liquid_mixture_properties(mixture_components,mixture_y,p%temperature,mixture_rho,mixture_cp,mixture_e0,ok)
    call require(ok,'trace initial energy')
    call advance_mixture_frozen_gas(species,transport,primitive,600.0_dp,p,mixture_components,component_mass, &
      options,1e-9_dp,totals0,ok)
    call require(ok.and.component_mass(1)==0.and.component_mass(2)>0,'trace component fully evaporates')
    call require(totals0(5+liquid%vapor_index)==q%multiplicity*saved_mass(1),'trace mass returned to gas species')
    call liquid_mixture_properties(mixture_components,component_mass/sum(component_mass),p%temperature, &
      mixture_rho,mixture_cp,mixture_e,ok)
    call require(ok,'trace remaining energy')
    expected=q%mass*q%multiplicity*(mixture_e0+.5_dp*sum(q%velocity**2))
    value=p%mass*p%multiplicity*(mixture_e+.5_dp*sum(p%velocity**2))
    call require(abs(value+totals0(5)-expected)<1e-12_dp*abs(expected),'trace cleanup energy conservation')
    p=q;component_mass=saved_mass
    primitive(15)=.1_dp;primitive(5+liquid%vapor_index)=.9_dp
    call advance_mixture_frozen_gas(species,transport,primitive,600.0_dp,p,mixture_components,component_mass, &
      options,1e-9_dp,totals0,ok)
    call require(ok.and.component_mass(1)==saved_mass(1),'suppressed trace remains in liquid')
    p=q;component_mass=p%mass*[.6_dp,.4_dp];saved_mass=component_mass
    primitive=0;primitive(1)=1;primitive(5)=101325;primitive(15)=1
    options%minimum_diameter=2e-4_dp
    call liquid_mixture_properties(mixture_components,[.6_dp,.4_dp],p%temperature,mixture_rho,mixture_cp,mixture_e0,ok)
    call require(ok,'terminal initial energy')
    call advance_mixture_frozen_gas(species,transport,primitive,600.0_dp,p,mixture_components,component_mass, &
      options,1e-9_dp,totals0,ok)
    call require(ok.and.p%mass==0.and.all(component_mass==0),'minimum diameter consumes all components')
    expected=q%mass*q%multiplicity*(mixture_e0+.5_dp*sum(q%velocity**2))
    call require(abs(totals0(5)-expected)<1e-12_dp*abs(expected),'terminal full energy returned')
    call require(abs(totals0(1)-q%mass*q%multiplicity)<1e-22_dp,'terminal full mass returned')
    p=q;component_mass=saved_mass;options%evaporation=.false.
    call advance_mixture_frozen_gas(species,transport,primitive,600.0_dp,p,mixture_components,component_mass, &
      options,1e-9_dp,totals0,ok)
    call require(ok.and.all(component_mass==saved_mass).and.p%mass==q%mass,'evaporation disabled preserves small drops')
  case('mixture_rates')
    mixture_components(1)=liquid
    call initialize_spray_liquid(species,'H2O2',800.0_dp,2000.0_dp,350.0_dp,101325.0_dp,1e6_dp, &
      mixture_components(2),ok)
    call require(ok,'mixture rate component')
    primitive=0;primitive(1)=1;primitive(2)=5;primitive(5)=101325;primitive(15)=1
    mixture_y=[1.0_dp,0.0_dp]
    do n=1,2
      options%abramzon_sirignano=n==2
      call build_spray_film(species,transport,primitive,600.0_dp,p,liquid,options,film,ok)
      call require(ok,'pure comparison film')
      call spray_transfer_rates(p,liquid,film,options,inverse_tau,heat,masscoef,ok)
      call require(ok,'pure comparison rates')
      call mixture_transfer_rates(species,transport,primitive,600.0_dp,p,mixture_components,mixture_y,options, &
        value,td,mixture_rates,vapor_h,ok)
      call require(ok,'pure-limit mixture rates')
      call require(abs(value/inverse_tau-1)<2e-13_dp,'pure-limit drag agreement')
      call require(abs(td/heat-1)<2e-13_dp,'pure-limit heat agreement')
      call require(abs(mixture_rates(1)/(masscoef*parcel_diameter(p,liquid))-1)<2e-13_dp,'pure-limit evaporation')
      call require(mixture_rates(2)==0,'inactive component has no rate')
    end do
    mixture_components(1)%density_model='polynomial'
    mixture_components(1)%density_coefficients=[1035.0_dp,-.1_dp,0.0_dp,0.0_dp]
    mixture_components(1)%density_temperature_range=[300.0_dp,600.0_dp]
    mixture_components(1)%heat_capacity_model='polynomial'
    mixture_components(1)%heat_capacity_coefficients=[3500.0_dp,2.0_dp,0.0_dp,0.0_dp]
    mixture_components(1)%heat_capacity_temperature_range=[300.0_dp,600.0_dp]
    mixture_y=[.8_dp,.2_dp];component_mass=p%mass*mixture_y;saved_mass=component_mass;q=p
    call mixture_transfer_rates(species,transport,primitive,600.0_dp,p,mixture_components,mixture_y,options, &
      inverse_tau,heat,mixture_rates,vapor_h,ok)
    call require(ok.and.all(mixture_rates>0),'two vapor components evaporate')
    call require(abs(mixture_rates(1)/sum(mixture_rates)-mixture_y(1))>.01_dp,'preferential evaporation')
    call liquid_mixture_properties(mixture_components,mixture_y,p%temperature,mixture_rho,mixture_cp,mixture_e0,ok)
    call require(ok,'mixture initial energy')
    call advance_mixture_frozen_gas(species,transport,primitive,600.0_dp,p,mixture_components,component_mass, &
      options,0.0_dp,totals0,ok)
    call require(ok.and.all(p%velocity==q%velocity).and.all(p%position==q%position).and. &
      all(component_mass==saved_mass).and.all(totals0==0),'zero interval preserves mixture exactly')
    call advance_mixture_frozen_gas(species,transport,primitive,600.0_dp,p,mixture_components,component_mass, &
      options,1e-7_dp,totals0,ok)
    call require(ok.and.all(component_mass<saved_mass),'frozen gas updates component masses')
    call require(p%velocity(1)>q%velocity(1).and.p%position(1)>q%position(1),'mixture drag and flight')
    call liquid_mixture_properties(mixture_components,component_mass/sum(component_mass),p%temperature, &
      mixture_rho,mixture_cp,mixture_e,ok)
    call require(ok,'mixture final energy')
    expected=q%mass*q%multiplicity*(mixture_e0+.5_dp*sum(q%velocity**2))
    value=p%mass*p%multiplicity*(mixture_e+.5_dp*sum(p%velocity**2))
    call require(abs(value+totals0(5)-expected)<1e-12_dp*abs(expected),'mixture total energy with drag')
    call require(abs(p%mass*p%multiplicity+totals0(1)-q%mass*q%multiplicity)<1e-20_dp,'mixture transported mass')
    call require(maxval(abs(p%mass*p%multiplicity*p%velocity+totals0(2:4)- &
      q%mass*q%multiplicity*q%velocity))<1e-20_dp,'mixture transported momentum')
    call require(abs(sum(totals0(6:))-totals0(1))<1e-25_dp,'species sum equals total mass transfer')
    q=p;saved_mass=component_mass
    call advance_mixture_frozen_gas(species,transport,primitive,600.0_dp,p,mixture_components,component_mass, &
      options,1e3_dp,totals0,ok)
    call require(.not.ok.and.all(component_mass==saved_mass).and.all(p%position==q%position).and. &
      all(p%velocity==q%velocity).and.p%temperature==q%temperature.and.all(totals0==0),'large-step mixture rollback')
    primitive(15)=.1_dp;primitive(5+mixture_components(2)%vapor_index)=.9_dp
    call mixture_transfer_rates(species,transport,primitive,600.0_dp,p,mixture_components,mixture_y,options, &
      inverse_tau,heat,mixture_rates,vapor_h,ok)
    call require(ok.and.mixture_rates(1)>0.and.mixture_rates(2)==0,'ambient vapor suppresses one component')
  case('mixture_exchange')
    mixture_components(1)=liquid
    call initialize_spray_liquid(species,'H2O2',800.0_dp,2000.0_dp,350.0_dp,101325.0_dp,1e6_dp, &
      mixture_components(2),ok)
    call require(ok,'mixture exchange component')
    mw=[species(liquid%vapor_index)%molecular_weight,species(mixture_components(2)%vapor_index)%molecular_weight]
    mixture_y=[.8_dp,.2_dp];surface_x=-1
    ! Pure component 2 is saturated, but its dilute liquid partial pressure is below ambient.
    call liquid_mixture_surface(mixture_components,mixture_y,mw,350.0_dp,101325.0_dp,surface_x,ok)
    call require(ok.and.all(surface_x>0).and.sum(surface_x)<1,'Raoult surface accepts diluted saturation')
    expected=(mixture_y(2)/mw(2))/sum(mixture_y/mw)
    call require(abs(surface_x(2)-expected)<1e-14_dp,'Raoult uses mole rather than mass fraction')
    call saturation_mole_fraction(liquid,mw(1),350.0_dp,101325.0_dp,value,ok)
    call require(ok.and.abs(surface_x(1)-value*(1-expected))<1e-14_dp,'Raoult partial pressure')
    surface_x=-1
    call liquid_mixture_surface(mixture_components,mixture_y,mw,350.0_dp,1.0_dp,surface_x,ok)
    call require(.not.ok.and.all(surface_x==-1),'saturated mixture output rollback')
    p%mass=1e-9_dp;p%multiplicity=10;p%temperature=350;q=p
    component_mass=p%mass*[.3_dp,.7_dp];saved_mass=component_mass
    evaporated_mass=p%mass*[.05_dp,.02_dp]
    do n=1,2
      vapor_h(n)=liquid_energy(mixture_components(n),350.0_dp)+1e6_dp
    end do
    expected=sum([(component_mass(n)*liquid_energy(mixture_components(n),350.0_dp),n=1,2)])
    heat=p%multiplicity*(sum([((component_mass(n)-evaporated_mass(n))* &
      liquid_energy(mixture_components(n),360.0_dp),n=1,2)])-expected+sum(evaporated_mass*vapor_h))
    call advance_mixture_evaporation(p,mixture_components,component_mass,evaporated_mass,vapor_h,heat,totals0,ok)
    call require(ok.and.abs(p%temperature-360)<1e-10_dp,'species loss and nonlinear energy update')
    call require(all(component_mass==saved_mass-evaporated_mass),'individual component losses')
    call require(abs(p%mass-sum(component_mass))<1e-24_dp,'parcel mass equals component sum')
    do n=1,2
      call require(totals0(5+mixture_components(n)%vapor_index)==q%multiplicity*evaporated_mass(n),'gas species routing')
    end do
    call require(abs(totals0(1)+p%mass*p%multiplicity-q%mass*q%multiplicity)<1e-23_dp,'mixture mass balance')
    expected=q%multiplicity*sum(evaporated_mass*vapor_h)-heat+totals0(1)*.5_dp*sum(q%velocity**2)
    call require(abs(totals0(5)-expected)<1e-12_dp*max(abs(expected),1e-12_dp),'mixture energy exchange')
    call require(all(totals0(2:4)==totals0(1)*q%velocity),'evaporated momentum')
    q=p;saved_mass=component_mass;evaporated_mass=component_mass;evaporated_mass(1)=2*component_mass(1)
    call advance_mixture_evaporation(p,mixture_components,component_mass,evaporated_mass,vapor_h,0.0_dp,totals0,ok)
    call require(.not.ok.and.all(component_mass==saved_mass).and.p%mass==q%mass.and.all(totals0==0),'overdraw rollback')
    evaporated_mass=[component_mass(1),0.0_dp]
    heat=p%multiplicity*component_mass(1)*(vapor_h(1)-liquid_energy(mixture_components(1),p%temperature))
    call advance_mixture_evaporation(p,mixture_components,component_mass,evaporated_mass,vapor_h,heat,totals0,ok)
    call require(ok.and.component_mass(1)==0.and.component_mass(2)>0,'one species completely depleted')
    call require(abs(p%temperature-q%temperature)<1e-10_dp,'depletion respects remaining composition energy')
    q=p;evaporated_mass=component_mass
    expected=q%mass*q%multiplicity*(liquid_energy(mixture_components(2),q%temperature)+.5_dp*sum(q%velocity**2))
    call advance_mixture_evaporation(p,mixture_components,component_mass,evaporated_mass,vapor_h,0.0_dp,totals0,ok)
    call require(ok.and.p%mass==0.and.all(component_mass==0),'complete mixture evaporation')
    call require(abs(totals0(5)-expected)<1e-12_dp*max(abs(expected),1e-12_dp),'terminal full inventory transfer')
  case('mixture')
    mixture_components(1)=liquid
    call initialize_spray_liquid(species,'H2O2',800.0_dp,2000.0_dp,350.0_dp,101325.0_dp,1e6_dp, &
      mixture_components(2),ok)
    call require(ok,'illustrative second mixture component')
    mixture_y=[.25_dp,.75_dp];mixture_t=330
    call liquid_mixture_properties(mixture_components,mixture_y,mixture_t,mixture_rho,mixture_cp,mixture_e,ok)
    call require(ok,'mixture properties')
    call require(abs(mixture_rho/(1/(.25_dp/1000+.75_dp/800))-1)<1e-14_dp,'additive liquid volume')
    call require(abs(mixture_cp-(.25_dp*4184+.75_dp*2000))<1e-12_dp,'mass weighted cp')
    expected=.25_dp*liquid_energy(liquid,330.0_dp)+.75_dp*liquid_energy(mixture_components(2),330.0_dp)
    call require(abs(mixture_e-expected)<1e-9_dp,'common formation energy reference')
    mixture_t=1000
    call liquid_mixture_temperature(mixture_components,mixture_y,mixture_e,mixture_t,ok)
    call require(ok.and.abs(mixture_t-330)<1e-10_dp,'mixture energy inversion')
    mixture_components(2)%heat_capacity_model='polynomial'
    mixture_components(2)%heat_capacity_coefficients=[1000.0_dp,2.0_dp,.001_dp,0.0_dp]
    mixture_components(2)%heat_capacity_temperature_range=[300.0_dp,600.0_dp]
    call require(valid_spray_liquid(mixture_components(2)),'polynomial component')
    p%temperature=330;q=p
    call liquid_mixture_properties(mixture_components,mixture_y,p%temperature,mixture_rho,mixture_cp,mixture_e0,ok)
    call require(ok,'initial polynomial mixture')
    call liquid_mixture_properties(mixture_components,mixture_y,450.0_dp,mixture_rho,mixture_cp,mixture_e,ok)
    call require(ok,'target polynomial mixture')
    heat=p%mass*p%multiplicity*(mixture_e-mixture_e0)
    call advance_mixture_heating(p,mixture_components,mixture_y,heat,solid_energy,ok)
    call require(ok.and.abs(p%temperature-450)<1e-10_dp,'composition-aware nonlinear heating')
    call require(abs(solid_energy+heat)<1e-12_dp*abs(heat),'mixture sensible energy closure')
    call require(p%mass==q%mass.and.p%multiplicity==q%multiplicity,'heating preserves liquid mass')
    q=p
    call advance_mixture_heating(p,mixture_components,mixture_y,1e10_dp,solid_energy,ok)
    call require(.not.ok.and.p%temperature==q%temperature.and.solid_energy==0,'out of range heat rollback')
    mixture_t=330
    call liquid_mixture_temperature(mixture_components,mixture_y,1e15_dp,mixture_t,ok)
    call require(.not.ok.and.mixture_t==330,'unreachable energy preserves temperature')
    mixture_rho=-1;mixture_cp=-2;mixture_e=-3
    call liquid_mixture_properties(mixture_components,[.25_dp,.7_dp],330.0_dp,mixture_rho,mixture_cp,mixture_e,ok)
    call require(.not.ok.and.mixture_rho==-1.and.mixture_cp==-2.and.mixture_e==-3,'composition rejection transaction')
    call liquid_mixture_properties(mixture_components,[1.0_dp,0.0_dp],700.0_dp,mixture_rho,mixture_cp,mixture_e,ok)
    call require(ok.and.mixture_rho==1000.and.mixture_cp==4184,'pure limit ignores inactive fit range')
  case('film_evolution')
    faces='no_slip_wall';wall_settings%enabled(5)=.true.;wall_settings%evolve=.true.
    wall_settings%conductivity_model='polynomial'
    wall_settings%conductivity_coefficients=[.1_dp,.001_dp,1e-6_dp,1e-9_dp]
    wall_settings%conductivity_temperature_range=[300.0_dp,600.0_dp]
    call require(valid_wall_options(wall_settings,faces),'valid cubic conductivity')
    call require(abs(wall_liquid_conductivity_at_temperature(wall_settings,1,400.0_dp)-.724_dp)<1e-14_dp,'cubic lambda')
    wall_settings%port_conductivity_models(2)='constant';wall_settings%port_conductivities(2)=.2_dp
    call require(valid_wall_options(wall_settings,faces),'port constant replaces global fit')
    call require(abs(wall_liquid_conductivity_at_temperature(wall_settings,2,400.0_dp)-.2_dp)<1e-14_dp,'port lambda')
    call require(wall_liquid_conductivity_at_temperature(wall_settings,1,299.0_dp)==0,'lambda range rejection')
    call require(wall_liquid_conductivity_at_temperature(wall_settings,0,400.0_dp)==0,'lambda material rejection')
    wall_settings%conductivity_coefficients=[15.5_dp,-.08_dp,.0001_dp,0.0_dp]
    call require(.not.valid_wall_options(wall_settings,faces),'reject negative interior with positive endpoints')

    call construct_wall_film(p,liquid,p%position,[0.0_dp,0.0_dp,1.0_dp],.6_dp,100.0_dp,1000.0_dp,pi/3, &
      wall_patch,film_budget,ok)
    call require(ok,'initial evolving film')
    call wall_film_transfer_rates(wall_patch,liquid,film,.6_dp,400.0_dp,.true.,.false., &
      film_mdot,film_qgas,film_qwall,film_tdot,ok)
    call require(ok.and.film_mdot==0.and.film_qgas>0.and.film_qwall>0,'film receives gas and wall heat')
    dt=1e-8_dp
    call advance_wall_film_frozen(wall_patch,liquid,film,.6_dp,400.0_dp,.true.,.false.,dt,1e-7_dp, &
      saved_patch,transfer,solid_energy,ok)
    call require(ok.and.saved_patch%mass==wall_patch%mass,'pure film heating')
    call require(abs(saved_patch%temperature-wall_patch%temperature-dt*film_tdot)<1e-12_dp,'explicit thermal rate')
    child_totals=wall_film_totals(saved_patch,liquid)-wall_film_totals(wall_patch,liquid)+transfer
    call require(abs(child_totals(5)+solid_energy)<1e-12_dp,'gas film wall energy closure')
    call require(abs(solid_energy+dt*film_qwall)<1e-15_dp,'wall heat ledger')
    call advance_wall_film_frozen(wall_patch,liquid,film,.6_dp,400.0_dp,.true.,.true.,dt,1e-7_dp, &
      saved_patch,transfer,solid_energy,ok)
    call require(ok.and.saved_patch%mass<wall_patch%mass.and.transfer(1)>0,'film evaporation supplies gas mass')
    call require(saved_patch%diameter==wall_patch%diameter,'film footprint fixed')
    call require(abs(saved_patch%mass+transfer(1)-wall_patch%mass)<1e-20_dp,'film mass closure')
    call advance_wall_film_frozen(wall_patch,liquid,film,-1.0_dp,400.0_dp,.true.,.true.,dt,1e-7_dp, &
      saved_patch,transfer,solid_energy,ok)
    call require(.not.ok.and.saved_patch%mass==wall_patch%mass.and.all(transfer==0),'invalid conduction rollback')
  case('partition')
    call require(spray_cartesian_owner([3.0_dp,0.0_dp,0.0_dp],[0.0_dp,0.0_dp,0.0_dp], &
      [1.0_dp,1.0_dp,1.0_dp],[5,3,2],[2,2,1],[.false.,.false.,.false.])==1,'uneven cell partition boundary')
    call require(spray_cartesian_owner([-1.0_dp,2.0_dp,0.0_dp],[0.0_dp,0.0_dp,0.0_dp], &
      [1.0_dp,1.0_dp,1.0_dp],[5,3,2],[2,2,1],[.true.,.false.,.false.])==3,'periodic Cartesian owner')
    allocate(drops(3),child_component_mass(2,3))
    do n=1,3
      drops(n)=p;drops(n)%id=int(n,int64);drops(n)%position=[4.0_dp,2.0_dp,0.0_dp]
      drops(n)%khrt_initial_count=100;drops(n)%khrt_shed_volume=real(n,dp)*1e-15_dp
      child_component_mass(:,n)=p%mass*[.2_dp*n,1-.2_dp*n]
    end do
    drops(2)%position=[0.0_dp,0.0_dp,0.0_dp];saved_drops=drops;saved_child_mass=child_component_mass
    call partition_spray_parcels(drops,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp], &
      [5,3,2],[2,2,1],[.false.,.false.,.false.],partition_counts,partition_offsets,ok,child_component_mass)
    call require(ok.and.all(partition_counts==[1,0,0,2]).and.all(partition_offsets==[0,1,1,1]),'rank send segments')
    call require(all(drops%id==[2,1,3]).and.all(drops%khrt_shed_volume==saved_drops([2,1,3])%khrt_shed_volume), &
      'stable rank ordering preserves KH state')
    call require(all(child_component_mass==saved_child_mass(:,[2,1,3])),'partition preserves component columns')
    drops=saved_drops;child_component_mass=saved_child_mass;drops(3)%position(1)=5
    call partition_spray_parcels(drops,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp], &
      [5,3,2],[2,2,1],[.false.,.false.,.false.],partition_counts,partition_offsets,ok,child_component_mass)
    call require(.not.ok.and.all(drops%id==[1,2,3]).and.all(child_component_mass==saved_child_mass),'invalid owner rollback')
    call require(all(partition_counts==[1,0,0,2]),'failed partition retains caller counts')
  case('mixture_wall_crossing')
    mixture_components(1)=liquid
    call initialize_spray_liquid(species,'H2O2',800.0_dp,2000.0_dp,350.0_dp,101325.0_dp,1e6_dp, &
      mixture_components(2),ok)
    call require(ok,'crossing mixture component')
    mixture_y=[.6_dp,.4_dp];mw=species([liquid%vapor_index,mixture_components(2)%vapor_index])%molecular_weight
    call mixture_impact_boiling_temperature(mixture_components,mixture_y,mw,101325.0_dp,value,ok)
    call require(ok.and.abs(value-(.6_dp*373.15_dp+.4_dp*350.0_dp))<1e-12_dp,'reference weighted impact boiling')
    options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
    wall_settings%enabled=.true.;wall_settings%temperature=300
    q=p;p%position(3)=-1e-5_dp;p%velocity=[0.0_dp,0.0_dp,-.1_dp]
    component_mass=p%mass*mixture_y;saved_mass=component_mass
    call mixture_impact_at_planar_wall(p,q%position,mixture_components,component_mass,options,mw,101325.0_dp, &
      [0.0_dp,0.0_dp,0.0_dp],[.004_dp,.004_dp,.004_dp],wall_settings,empty_films,2_int64, &
      drops,child_component_mass,wall_patch,film_component_mass,transfer,hit,ok)
    call require(ok.and.hit.and.p%mass==0.and.size(drops)==0,'crossing deposits mixture')
    call require(maxval(abs(film_component_mass-q%multiplicity*saved_mass))<1e-22_dp,'crossing film components')
    call require(all(wall_patch%normal==[0.0_dp,0.0_dp,1.0_dp]).and.wall_patch%position(3)==0,'crossing wall geometry')
    saved_patch=wall_patch;saved_film_mass=film_component_mass
    p=q;component_mass=saved_mass
    call mixture_impact_at_planar_wall(p,q%position,mixture_components,component_mass,options,mw,101325.0_dp, &
      [0.0_dp,0.0_dp,0.0_dp],[.004_dp,.004_dp,.004_dp],wall_settings,empty_films,2_int64, &
      drops,child_component_mass,wall_patch,film_component_mass,transfer,hit,ok)
    call require(ok.and..not.hit.and.all(component_mass==saved_mass),'no crossing retains mixture')
    p%position(3)=-1e-5_dp;p%velocity=[0.0_dp,0.0_dp,-.1_dp]
    call mixture_impact_at_planar_wall(p,q%position,mixture_components,component_mass,options,mw,-1.0_dp, &
      [0.0_dp,0.0_dp,0.0_dp],[.004_dp,.004_dp,.004_dp],wall_settings,empty_films,2_int64, &
      drops,child_component_mass,wall_patch,film_component_mass,transfer,hit,ok)
    call require(.not.ok.and..not.hit.and.all(component_mass==saved_mass).and.p%mass==q%mass,'crossing failure rollback')
    call require(all(film_component_mass==saved_film_mass).and.wall_patch%mass==saved_patch%mass,'crossing film rollback')
  case('mixture_wall_rates')
    mixture_components(1)=liquid
    call initialize_spray_liquid(species,'H2O2',800.0_dp,2000.0_dp,350.0_dp,101325.0_dp,1e6_dp, &
      mixture_components(2),ok)
    call require(ok,'film rate component')
    mixture_y=[.6_dp,.4_dp]
    call liquid_mixture_properties(mixture_components,mixture_y,350.0_dp,mixture_rho,mixture_cp,mixture_e0,ok)
    call require(ok,'film rate properties')
    wall_patch=spray_wall_film([0.0_dp,0.0_dp,0.0_dp],[0.0_dp,0.0_dp,1.0_dp],1e-8_dp,350.0_dp, &
      1e-3_dp,1e-8_dp/(mixture_rho*.25_dp*pi*1e-6_dp),1)
    component_mass=wall_patch%mass*mixture_y;saved_mass=component_mass;saved_patch=wall_patch
    primitive=0;primitive(1)=1;primitive(5)=101325;primitive(15)=1
    call mixture_wall_film_rates(wall_patch,mixture_components,component_mass,species,transport,primitive,600.0_dp, &
      options,.6_dp,400.0_dp,mixture_rates,vapor_h,heat,value,td,ok)
    call require(ok.and.all(mixture_rates>0).and.heat>0.and.value>0,'mixed film heat and evaporation rates')
    expected=(heat+value-sum(mixture_rates*(vapor_h-[liquid_energy(mixture_components(1),350.0_dp), &
      liquid_energy(mixture_components(2),350.0_dp)])))/(wall_patch%mass*mixture_cp)
    call require(abs(td-expected)<1e-12_dp*max(1.0_dp,abs(expected)),'composition-aware film cooling')
    dt=1e-10_dp
    call advance_mixture_wall_film_frozen(wall_patch,mixture_components,component_mass,species,transport,primitive,600.0_dp, &
      options,.6_dp,400.0_dp,dt,1e-15_dp,totals0,heat,ok)
    call require(ok.and.wall_patch%mass<saved_patch%mass,'frozen mixed film update')
    call require(maxval(abs(component_mass-(saved_mass-mixture_rates*dt)))<1e-24_dp,'frozen film component rates')
    call liquid_mixture_properties(mixture_components,component_mass/wall_patch%mass,wall_patch%temperature, &
      mixture_rho,mixture_cp,mixture_e,ok)
    call require(ok,'updated film energy')
    call require(abs(wall_patch%mass*mixture_e+totals0(5)+heat-saved_patch%mass*mixture_e0)<1e-14_dp, &
      'frozen film gas wall closure')
    wall_patch=saved_patch;component_mass=saved_mass
    call advance_mixture_wall_film_frozen(wall_patch,mixture_components,component_mass,species,transport,primitive,600.0_dp, &
      options,.6_dp,400.0_dp,1.0_dp,1e-15_dp,totals0,heat,ok)
    call require(.not.ok.and.all(component_mass==saved_mass).and.wall_patch%mass==saved_patch%mass,'film large step rollback')
    mixture_y=[1.0_dp,0.0_dp];component_mass=wall_patch%mass*mixture_y
    call mixture_wall_film_rates(wall_patch,mixture_components,component_mass,species,transport,primitive,600.0_dp, &
      options,.6_dp,400.0_dp,mixture_rates,vapor_h,heat,value,td,ok)
    call require(ok.and.mixture_rates(2)==0,'pure limit inactive component')
    p%temperature=wall_patch%temperature
    call build_spray_film(species,transport,primitive,600.0_dp,p,liquid,options,film,ok,wall_skin=.true.)
    call require(ok,'pure film skin')
    call wall_film_transfer_rates(wall_patch,liquid,film,.6_dp,400.0_dp,.true.,.true.,masscoef, &
      transfer(1),transfer(2),transfer(3),ok)
    call require(ok.and.abs(mixture_rates(1)/masscoef-1)<1e-12_dp,'pure limit film mass rate')
    call require(abs(heat/transfer(1)-1)<1e-12_dp.and.value==transfer(2),'pure limit heat rates')
    wall_patch=saved_patch;mixture_y=[1e-13_dp,1.0_dp-1e-13_dp]
    component_mass=wall_patch%mass*mixture_y;saved_mass=component_mass
    call advance_mixture_wall_film_frozen(wall_patch,mixture_components,component_mass,species,transport,primitive,600.0_dp, &
      options,.6_dp,400.0_dp,0.0_dp,1e-15_dp,totals0,heat,ok)
    call require(ok.and.all(component_mass==saved_mass).and.all(totals0==0),'zero interval preserves trace film')
    call advance_mixture_wall_film_frozen(wall_patch,mixture_components,component_mass,species,transport,primitive,600.0_dp, &
      options,.6_dp,400.0_dp,1e-10_dp,1e-15_dp,totals0,heat,ok)
    call require(ok.and.component_mass(1)==0.and.component_mass(2)>0,'active film trace fully released')
    call require(totals0(5+liquid%vapor_index)==saved_mass(1),'film trace gas species receives whole residual')
    wall_patch=saved_patch;component_mass=saved_mass
    primitive(5+liquid%vapor_index)=.001_dp;primitive(15)=.999_dp
    call advance_mixture_wall_film_frozen(wall_patch,mixture_components,component_mass,species,transport,primitive,600.0_dp, &
      options,.6_dp,400.0_dp,1e-10_dp,1e-15_dp,totals0,heat,ok)
    call require(ok.and.component_mass(1)==saved_mass(1),'ambient-suppressed film trace retained')
    call require(totals0(5+liquid%vapor_index)==0,'suppressed trace has no gas source')
  case('mixture_wall_exchange')
    mixture_components(1)=liquid
    call initialize_spray_liquid(species,'H2O2',800.0_dp,2000.0_dp,350.0_dp,101325.0_dp,1e6_dp, &
      mixture_components(2),ok)
    call require(ok,'film mixture component')
    mixture_y=[.6_dp,.4_dp]
    call liquid_mixture_properties(mixture_components,mixture_y,350.0_dp,mixture_rho,mixture_cp,mixture_e0,ok)
    call require(ok,'initial mixture film properties')
    wall_patch=spray_wall_film([0.0_dp,0.0_dp,0.0_dp],[0.0_dp,0.0_dp,1.0_dp],1e-8_dp,350.0_dp, &
      1e-3_dp,1e-8_dp/(mixture_rho*.25_dp*pi*1e-6_dp),1)
    saved_patch=wall_patch;component_mass=wall_patch%mass*mixture_y;saved_mass=component_mass
    evaporated_mass=component_mass*[.01_dp,.02_dp];vapor_h=[1e6_dp,2e6_dp]
    call advance_mixture_wall_film_exchange(wall_patch,mixture_components,component_mass,evaporated_mass,vapor_h, &
      1e-6_dp,2e-6_dp,1e-15_dp,totals0,heat,ok)
    call require(ok.and.wall_patch%mass<saved_patch%mass,'partial mixed film evaporation')
    call require(all(component_mass==saved_mass-evaporated_mass),'film component loss')
    call require(maxval(abs(totals0([5+liquid%vapor_index,5+mixture_components(2)%vapor_index])-evaporated_mass)) &
      <1e-24_dp,'film gas species gain')
    call liquid_mixture_properties(mixture_components,component_mass/wall_patch%mass,wall_patch%temperature, &
      mixture_rho,mixture_cp,mixture_e,ok)
    call require(ok,'remaining film properties')
    call require(abs(wall_patch%mass*mixture_e+totals0(5)+heat-saved_patch%mass*mixture_e0)<1e-14_dp, &
      'film gas wall energy closure')
    call require(abs(wall_patch%height*mixture_rho*.25_dp*pi*wall_patch%diameter**2/wall_patch%mass-1)<1e-13_dp, &
      'film evolving density geometry')
    wall_patch=saved_patch;component_mass=saved_mass
    call advance_mixture_wall_film_exchange(wall_patch,mixture_components,component_mass,2*saved_mass,vapor_h, &
      0.0_dp,0.0_dp,1e-15_dp,totals0,heat,ok)
    call require(.not.ok.and.all(component_mass==saved_mass).and.wall_patch%mass==saved_patch%mass,'film overdraw rollback')
    call require(all(totals0==0).and.heat==0,'failed film exchange zero budgets')
    do event_index=1,2
      wall_patch=saved_patch;component_mass=saved_mass
      evaporated_mass=merge(saved_mass,.01_dp*saved_mass,event_index==1)
      call advance_mixture_wall_film_exchange(wall_patch,mixture_components,component_mass,evaporated_mass,vapor_h, &
        0.0_dp,2e-6_dp,1e-3_dp,totals0,heat,ok)
      call require(ok.and.wall_patch%mass==0.and.wall_patch%height==0.and.all(component_mass==0),'film terminal inventory')
      call require(abs(totals0(1)-saved_patch%mass)<1e-24_dp,'terminal film gas mass')
      call require(abs(totals0(5)+heat-saved_patch%mass*mixture_e0)<1e-14_dp,'terminal film energy closure')
    end do
  case('mixture_wall_event')
    mixture_components(1)=liquid
    call initialize_spray_liquid(species,'H2O2',800.0_dp,2000.0_dp,350.0_dp,101325.0_dp,1e6_dp, &
      mixture_components(2),ok)
    call require(ok,'wall mixture component')
    mixture_components(1)%heat_capacity_model='polynomial'
    mixture_components(1)%heat_capacity_coefficients=[4000.0_dp,1.0_dp,0.0_dp,0.0_dp]
    mixture_components(1)%heat_capacity_temperature_range=[300.0_dp,600.0_dp]
    p%velocity=[0.0_dp,0.0_dp,-20.0_dp];q=p;mixture_y=[.6_dp,.4_dp]
    call liquid_mixture_properties(mixture_components,mixture_y,p%temperature,mixture_rho,mixture_cp,mixture_e,ok)
    call require(ok,'wall initial thermodynamics')
    impact_kv=[79.0_dp,79.0_dp,150.0_dp,150.0_dp,150.0_dp,300.0_dp]
    impact_tstar=[1.0_dp,1.1_dp,1.0_dp,1.1_dp,1.0_dp,1.0_dp]
    impact_film=[0.0_dp,0.0_dp,0.0_dp,0.0_dp,.2_dp,.03_dp]
    impact_kinds=[impact_deposit,impact_rebound,impact_dry,impact_thermal,impact_spread,impact_wet]
    do event_index=1,6
      p=q;component_mass=q%mass*mixture_y
      call resolve_mixture_wall_impact(p,mixture_components,component_mass,q%position,[0.0_dp,0.0_dp,1.0_dp], &
        impact_kv(event_index),impact_tstar(event_index),80.0_dp,impact_film(event_index),100.0_dp,1000.0_dp, &
        pi/3,2_int64,drops,child_component_mass,wall_patch,film_component_mass,transfer,regime,ok)
      call require(ok.and.regime==impact_kinds(event_index),'mixture impact regime')
      saved_mass=p%multiplicity*component_mass+film_component_mass
      expected=p%mass*p%multiplicity*(mixture_e+.5_dp*sum(p%velocity**2))+wall_patch%mass*mixture_e
      grad(:,1)=p%mass*p%multiplicity*p%velocity
      do n=1,size(drops)
        saved_mass=saved_mass+drops(n)%multiplicity*child_component_mass(:,n)
        expected=expected+drops(n)%mass*drops(n)%multiplicity*(mixture_e+.5_dp*sum(drops(n)%velocity**2))
        grad(:,1)=grad(:,1)+drops(n)%mass*drops(n)%multiplicity*drops(n)%velocity
        call require(maxval(abs(child_component_mass(:,n)/drops(n)%mass-mixture_y))<1e-14_dp,'splash composition')
      end do
      call require(maxval(abs(saved_mass/(q%mass*q%multiplicity)-mixture_y))<1e-13_dp,'wall component closure')
      value=q%mass*q%multiplicity*(mixture_e+.5_dp*sum(q%velocity**2))
      call require(abs(expected+transfer(5)-value)<1e-12_dp*abs(value),'mixture wall energy')
      call require(maxval(abs(grad(:,1)+transfer(2:4)-q%mass*q%multiplicity*q%velocity))<1e-18_dp, &
        'mixture wall momentum')
    end do
    saved_drops=drops;saved_child_mass=child_component_mass;saved_patch=wall_patch;saved_film_mass=film_component_mass
    p=q;component_mass=q%mass*mixture_y;saved_mass=component_mass
    call resolve_mixture_wall_impact(p,mixture_components,component_mass,q%position,[0.0_dp,0.0_dp,1.0_dp], &
      150.0_dp,1.0_dp,80.0_dp,0.0_dp,100.0_dp,0.0_dp,pi/3,2_int64,drops,child_component_mass, &
      wall_patch,film_component_mass,transfer,regime,ok)
    call require(.not.ok.and.regime==-1.and.all(transfer==0),'mixture late impact failure')
    call require(all(component_mass==saved_mass).and.p%mass==q%mass,'mixture parent rollback')
    call require(all(child_component_mass==saved_child_mass).and.all(drops%mass==saved_drops%mass), &
      'mixture children rollback')
    call require(all(film_component_mass==saved_film_mass).and.wall_patch%mass==saved_patch%mass,'mixture film rollback')
  case('wall_event')
    p%velocity=[0.0_dp,0.0_dp,-20.0_dp];p%material_id=2;q=p
    q%khrt_initial_count=q%multiplicity;q%khrt_shed_volume=1e-12_dp
    impact_kv=[79.0_dp,79.0_dp,150.0_dp,150.0_dp,150.0_dp,300.0_dp]
    impact_tstar=[1.0_dp,1.1_dp,1.0_dp,1.1_dp,1.0_dp,1.0_dp]
    impact_film=[0.0_dp,0.0_dp,0.0_dp,0.0_dp,.2_dp,.03_dp]
    impact_kinds=[impact_deposit,impact_rebound,impact_dry,impact_thermal,impact_spread,impact_wet]
    do event_index=1,6
      p=q
      call resolve_wall_impact(p,liquid,q%position,[0.0_dp,0.0_dp,1.0_dp],impact_kv(event_index), &
        impact_tstar(event_index),80.0_dp,impact_film(event_index),100.0_dp,1000.0_dp,pi/3, &
        2_int64,drops,wall_patch,transfer,regime,ok)
      if (.not.ok) print *, 'rejected impact case',event_index
      call require(ok.and.regime==impact_kinds(event_index),'complete impact regime dispatch')
      call require(all(drops%khrt_initial_count==0).and.all(drops%khrt_shed_volume==0),'fresh splash KH state')
      if (regime==impact_rebound) then
        call require(p%khrt_shed_volume==q%khrt_shed_volume,'rebound preserves unspent KH inventory')
      else
        call require(p%khrt_initial_count==0.and.p%khrt_shed_volume==0,'impact consumed/reset KH inventory')
      end if
      child_totals=parcel_totals(p,liquid)
      do n=1,size(drops)
        child_totals=child_totals+parcel_totals(drops(n),liquid)
      end do
      if (wall_patch%mass>0) child_totals=child_totals+wall_film_totals(wall_patch,liquid)
      call require(abs(child_totals(1)/(q%mass*q%multiplicity)-1)<1e-13_dp,'event total mass')
      call require(maxval(abs(child_totals+transfer-parcel_totals(q,liquid)))<1e-12_dp,'event wall budget')
      if (p%mass>0) call require(p%velocity(3)>0,'retained parent leaves wall')
      if (regime==impact_rebound) then
        expected=.823_dp*exp(-.0019835_dp*79.0_dp**1.6_dp)*20
        call require(abs(p%velocity(3)-expected)<1e-13_dp.and.size(drops)==0,'reference rebound restitution')
      end if
      if (regime==impact_thermal) call require(p%mass>0.and.size(drops)==4.and.wall_patch%mass==0,'thermal split')
      if (regime==impact_dry.or.regime==impact_wet) &
        call require(p%mass==0.and.size(drops)==4.and.wall_patch%mass>0,'splash plus deposition')
    end do
    saved_drops=drops;saved_patch=wall_patch;p=q
    call resolve_wall_impact(p,liquid,q%position,[0.0_dp,0.0_dp,1.0_dp],150.0_dp,1.0_dp,80.0_dp, &
      .03_dp,100.0_dp,1000.0_dp,pi/3,2_int64,drops,wall_patch,transfer,regime,ok)
    call require(.not.ok.and.regime==-1.and.all(transfer==0),'reject inward wet-splash child velocity')
    call require(p%mass==q%mass.and.all(drops%mass==saved_drops%mass),'wet velocity failure rollback')
    call resolve_wall_impact(p,liquid,q%position,[0.0_dp,0.0_dp,1.0_dp],150.0_dp,1.0_dp,80.0_dp, &
      0.0_dp,100.0_dp,0.0_dp,pi/3,2_int64,drops,wall_patch,transfer,regime,ok)
    call require(.not.ok.and.regime==-1.and.all(transfer==0),'late deposition failure')
    call require(p%mass==q%mass.and.all(p%velocity==q%velocity),'parent event rollback')
    call require(all(drops%mass==saved_drops%mass).and.wall_patch%mass==saved_patch%mass,'children and film rollback')
  case('wall_state')
    p%velocity=[0.0_dp,0.0_dp,-20.0_dp];p%material_id=2
    call construct_splash_children(p,liquid,p%position,[0.0_dp,0.0_dp,1.0_dp],150.0_dp,0.0_dp,.5_dp, &
      2_int64,drops,transfer,ok)
    call require(ok.and.size(drops)==4,'construct four splash children')
    child_totals=0
    do n=1,4
      call require(drops(n)%id==n+1.and.drops(n)%parent_id==p%id.and.drops(n)%material_id==2,'child identity')
      call require(drops(n)%multiplicity==p%multiplicity.and.drops(n)%temperature==p%temperature,'child population')
      call require(drops(n)%mass==.125_dp*p%mass,'normal impact equal dry mass partition')
      call require(drops(n)%velocity(3)>0.and.drops(n)%position(3)>p%position(3),'outgoing splash')
      child_totals=child_totals+parcel_totals(drops(n),liquid)
    end do
    call require(maxval(abs(child_totals+transfer-.5_dp*parcel_totals(p,liquid)))<1e-12_dp,'splash wall balance')
    saved_drops=drops
    call construct_splash_children(p,liquid,p%position,[0.0_dp,0.0_dp,1.0_dp],150.0_dp,.2_dp,.5_dp, &
      2_int64,drops,transfer,ok)
    call require(.not.ok.and.all(transfer==0).and.all(drops%mass==saved_drops%mass),'negative child fit rolls back')
    call construct_wall_film(p,liquid,p%position,[0.0_dp,0.0_dp,1.0_dp],.6_dp,100.0_dp,1000.0_dp,pi/3, &
      wall_patch,film_budget,ok)
    call require(ok.and.wall_patch%material_id==2,'construct deposited film')
    call require(abs(wall_patch%diameter/.0014595751300262862_dp-1)<1e-13_dp,'pinned film diameter')
    call require(abs(wall_patch%height/1.8776173162247343e-5_dp-1)<1e-13_dp,'pinned film height')
    call require(abs(.25_dp*pi*wall_patch%diameter**2*wall_patch%height*liquid%density/wall_patch%mass-1)<1e-13_dp, &
      'cylindrical volume conserves deposited mass')
    child_totals=wall_film_totals(wall_patch,liquid)
    call require(maxval(abs(child_totals+film_budget-.6_dp*parcel_totals(p,liquid)))<1e-12_dp,'film wall balance')
    call require(film_budget(5)>0.and.film_budget(3)==0,'stationary film transfers incident kinetic energy')
    saved_patch=wall_patch
    call construct_wall_film(p,liquid,p%position,[0.0_dp,0.0_dp,1.0_dp],.6_dp,100.0_dp,0.0_dp,pi/3, &
      wall_patch,film_budget,ok)
    call require(.not.ok.and.wall_patch%mass==saved_patch%mass.and.all(film_budget==0),'invalid film rolls back')
  case('wall_impact')
    call classify_wall_impact(150.0_dp,1.0_dp,pi/6,80.0_dp,0.0_dp,regime,splash_mass,deposited_mass,ok)
    call require(ok.and.regime==impact_dry,'dry splash regime')
    call require(abs(splash_mass-.14206181552235969_dp)<1e-14_dp,'pinned dry splash fraction')
    call require(splash_mass+deposited_mass==1,'dry splash/deposit closure')
    call splash_velocity_components(-20.0_dp,10.0_dp,150.0_dp,0.0_dp,splash_vel,ok)
    call require(ok.and.maxval(abs(splash_vel-[20.701826903833513_dp,8.6602540378443873_dp, &
      -5.7018269038335125_dp,12.701826903833513_dp,3.3142511671214492_dp]))<1e-13_dp,'pinned dry velocities')
    call splash_mass_fractions(pi/6,.5_dp,0.0_dp,azimuth_mass,ok)
    call require(ok.and.maxval(abs(azimuth_mass-[.125_dp,.21160254037844389_dp, &
      .125_dp,.038397459621556107_dp]))<1e-14_dp,'pinned dry azimuth masses')
    call classify_wall_impact(150.0_dp,1.0_dp,pi/6,80.0_dp,.03_dp,regime,splash_mass,deposited_mass,ok)
    call require(ok.and.regime==impact_wet.and.abs(splash_mass-.26793114657822259_dp)<1e-14_dp,'pinned thin wet splash')
    call splash_velocity_components(-20.0_dp,10.0_dp,150.0_dp,.03_dp,splash_vel,ok)
    call require(ok.and.maxval(abs(splash_vel-[3.2235_dp,1.0_dp,.1065_dp,.8235_dp,.2931_dp]))<1e-14_dp, &
      'pinned thin film velocities')
    call splash_velocity_components(-20.0_dp,10.0_dp,150.0_dp,.07_dp,splash_vel,ok)
    call require(ok.and.maxval(abs(splash_vel-[3.556_dp,1.0_dp,-.226_dp,1.156_dp,1.012_dp]))<1e-14_dp, &
      'pinned intermediate film velocities')
    call splash_mass_fractions(pi/6,.5_dp,.03_dp,azimuth_mass,ok)
    call require(ok.and.maxval(abs(azimuth_mass-[.10595_dp,.2505_dp,.10595_dp,.0376_dp]))<1e-14_dp, &
      'pinned thin film azimuth masses')
    call classify_wall_impact(150.0_dp,1.0_dp,pi/6,80.0_dp,.2_dp,regime,splash_mass,deposited_mass,ok)
    call require(ok.and.regime==impact_spread.and.deposited_mass==1,'wet spreading')
    call splash_mass_fractions(pi/6,.5_dp,.2_dp,azimuth_mass,ok)
    call require(.not.ok.and.all(azimuth_mass==0),'reject negative upstream thick-film azimuth fraction')
    call classify_wall_impact(79.0_dp,1.0_dp,pi/6,80.0_dp,0.0_dp,regime,splash_mass,deposited_mass,ok)
    call require(ok.and.regime==impact_deposit.and.deposited_mass==1,'cold deposition')
    call classify_wall_impact(79.0_dp,1.1_dp,pi/6,80.0_dp,0.0_dp,regime,splash_mass,deposited_mass,ok)
    call require(ok.and.regime==impact_rebound.and.splash_mass==0.and.deposited_mass==0,'hot rebound')
    call classify_wall_impact(150.0_dp,1.1_dp,pi/6,80.0_dp,0.0_dp,regime,splash_mass,deposited_mass,ok)
    call require(ok.and.regime==impact_thermal.and.splash_mass>0.and.deposited_mass==0,'thermal breakup')
    call classify_wall_impact(-1.0_dp,1.0_dp,pi/6,80.0_dp,0.0_dp,regime,splash_mass,deposited_mass,ok)
    call require(.not.ok.and.regime==-1,'invalid impact rejected')
  case('heat_capacity')
    liquid%heat_capacity_model='polynomial'
    liquid%heat_capacity_coefficients=[4000.0_dp,4.0_dp,0.0_dp,0.0_dp]
    liquid%heat_capacity_temperature_range=[300.0_dp,600.0_dp]
    call require(valid_spray_liquid(liquid),'valid cp fit')
    call require(liquid_heat_capacity_at_temperature(liquid,350.0_dp)==5400,'cp evaluation')
    expected=4000*(350.0_dp-373.15_dp)+2*(350.0_dp**2-373.15_dp**2)
    call require(abs(liquid_energy(liquid,350.0_dp)-liquid%reference_energy-expected)<1e-8_dp,'integrated cp energy')
    expected=(liquid_energy(liquid,350.001_dp)-liquid_energy(liquid,349.999_dp))/.002_dp
    call require(abs(expected-5400)<1e-5_dp,'energy derivative equals cp')
    options%drag=.false.;options%evaporation=.false.
    call spray_transfer_rates(p,liquid,film,options,inverse_tau,heat,masscoef,ok)
    call require(ok,'cp heat rate')
    dt=1e-5_dp
    expected=film%gas_temperature+(p%temperature-film%gas_temperature)*exp(-heat*dt/(p%mass*5400))
    call advance_parcel_frozen_film(p,liquid,film,options,dt,q,transfer,ok)
    call require(ok.and.abs(q%temperature-expected)<1e-12_dp,'local cp heating update')
    call require(abs(transfer(5)+p%mass*p%multiplicity* &
      (liquid_energy(liquid,q%temperature)-liquid_energy(liquid,p%temperature)))<1e-12_dp,'integrated energy transfer')
    constant_liquid=liquid
    liquid%heat_capacity_coefficients=[202499.0_dp,-900.0_dp,1.0_dp,0.0_dp]
    call require(.not.valid_spray_liquid(liquid),'reject negative interior cp extremum')
    liquid=constant_liquid
    call inject_cone_parcels(liquid,1e-10_dp,1,1_int64,p%position,[0.0_dp,0.0_dp,1.0_dp], &
      1.0_dp,0.0_dp,1e-5_dp,290.0_dp,drops,ok)
    call require(.not.ok,'reject injection outside cp interval')
  case('port_states')
    port_centers(:,1)=[.001_dp,.001_dp,.001_dp];port_centers(:,2)=[.002_dp,.001_dp,.001_dp]
    port_axes(:,1)=[1.0_dp,0.0_dp,0.0_dp];port_axes(:,2)=[-1.0_dp,0.0_dp,0.0_dp]
    port_fractions=[.25_dp,.75_dp]
    call inject_cone_ports(liquid,1e-9_dp,2,1_int64,port_centers,port_axes,port_fractions, &
      1.0_dp,0.0_dp,1e-4_dp,350.0_dp,drops,ok,port_speeds=[10.0_dp,20.0_dp], &
      port_diameters=[1e-5_dp,2e-5_dp],port_temperatures=[300.0_dp,350.0_dp],port_angles=[0.0_dp,0.0_dp])
    call require(ok.and.size(drops)==4,'distinct port states')
    call require(all(drops(1:2)%temperature==300).and.all(drops(3:4)%temperature==350),'port temperatures')
    call require(drops(1)%velocity(1)==10.and.drops(3)%velocity(1)==-20,'port speed and direction')
    call require(abs(parcel_diameter(drops(1),liquid)/1e-5_dp-1)<1e-13_dp.and. &
      abs(parcel_diameter(drops(3),liquid)/2e-5_dp-1)<1e-13_dp,'port diameters')
    transfer=0
    do n=1,4
      transfer=transfer+parcel_totals(drops(n),liquid)
    end do
    expected=2.5e-10_dp*(liquid_energy(liquid,300.0_dp)+50.0_dp)+ &
      7.5e-10_dp*(liquid_energy(liquid,350.0_dp)+200.0_dp)
    call require(abs(transfer(5)/expected-1)<1e-13_dp,'distinct port thermal and kinetic inventory')
    call inject_cone_ports(liquid,1e-9_dp,2,1_int64,port_centers,port_axes,port_fractions, &
      1.0_dp,0.0_dp,1e-4_dp,350.0_dp,drops,ok,port_speeds=[10.0_dp])
    call require(.not.ok.and..not.allocated(drops),'wrong port override dimensions')
  case('ports')
    port_centers(:,1)=[.001_dp,.001_dp,.001_dp];port_centers(:,2)=[.002_dp,.001_dp,.001_dp]
    port_axes(:,1)=[1.0_dp,0.0_dp,0.0_dp];port_axes(:,2)=[-1.0_dp,0.0_dp,0.0_dp]
    port_fractions=[.25_dp,.75_dp]
    call inject_cone_ports(liquid,1e-9_dp,3,11_int64,port_centers,port_axes,port_fractions, &
      10.0_dp,0.0_dp,1e-4_dp,350.0_dp,drops,ok)
    call require(ok.and.size(drops)==6,'multiport count')
    transfer=0
    do n=1,6
      call require(drops(n)%id==int(n+10,int64),'multiport ordered IDs')
      transfer=transfer+parcel_totals(drops(n),liquid)
    end do
    call require(abs(transfer(1)/1e-9_dp-1)<1e-13_dp,'total injection mass')
    call require(abs(transfer(2)/(-5e-9_dp)-1)<1e-13_dp,'weighted momentum')
    call require(abs(sum(drops(1:3)%mass*drops(1:3)%multiplicity)/2.5e-10_dp-1)<1e-13_dp,'first port share')
    call require(all(drops(1)%position==port_centers(:,1)).and.all(drops(4)%position==port_centers(:,2)), &
      'port centers')
    port_axes(:,2)=0
    call inject_cone_ports(liquid,1e-9_dp,3,11_int64,port_centers,port_axes,port_fractions, &
      10.0_dp,0.0_dp,1e-4_dp,350.0_dp,drops,ok)
    call require(.not.ok.and..not.allocated(drops),'invalid later port leaves no partial injection')
    port_fractions=[.25_dp,.5_dp]
    call inject_cone_ports(liquid,1e-9_dp,3,11_int64,port_centers,port_axes,port_fractions, &
      10.0_dp,0.0_dp,1e-4_dp,350.0_dp,drops,ok)
    call require(.not.ok,'port fractions must close')
  case('surface_cap')
    call saturation_mole_fraction(liquid,18.015_dp,400.0_dp,101325.0_dp,value,ok)
    call require(.not.ok,'default boiling rejection')
    liquid%boiling_model='surface-cap';liquid%critical_temperature=647.0_dp
    constant_liquid=liquid;constant_liquid%boiling_model='reject';constant_liquid%critical_temperature=0
    call saturation_mole_fraction(liquid,18.015_dp,400.0_dp,101325.0_dp,value,ok)
    call require(ok,'boiling surface cap')
    call saturation_mole_fraction(constant_liquid,18.015_dp,.99_dp*373.15_dp,101325.0_dp,expected,ok)
    call require(ok.and.abs(value/expected-1)<1e-13_dp,'normal-pressure surface temperature')
    do n=1,3
      z=101325.0_dp*10.0_dp**(n-2)
      td=1/(1/373.15_dp+log(101325.0_dp/z)*8314.46261815324_dp/(2.257e6_dp*18.015_dp))
      call saturation_mole_fraction(liquid,18.015_dp,600.0_dp,z,value,ok)
      call require(ok,'pressure-dependent surface cap')
      call saturation_mole_fraction(constant_liquid,18.015_dp,.99_dp*min(td,647.0_dp),z,expected,ok)
      call require(ok.and.abs(value/expected-1)<1e-12_dp,'pressure-dependent boiling temperature')
    end do
    liquid%critical_temperature=300
    call require(.not.valid_spray_liquid(liquid),'critical temperature consistency')
  case('viscosity')
    options%breakup_model='etab';options%surface_tension=.072_dp
    options%liquid_viscosity_model='inverse-temperature'
    options%liquid_viscosity_coefficients=[0.0002_dp,.2_dp,10.0_dp,100.0_dp]
    options%liquid_viscosity_temperature_range=[250.0_dp,600.0_dp]
    do n=0,7
      td=250.0_dp+50*n
      expected=.0002_dp+.2_dp/td+10.0_dp/td**2+100.0_dp/td**3
      call require(abs(liquid_viscosity_at_temperature(options,td)/expected-1)<1e-13_dp,'viscosity equation')
    end do
    do n=1,2
      if (n==2) options%breakup_model='khrt'
      constant_options=options;constant_options%liquid_viscosity_model='constant'
      constant_options%liquid_viscosity=liquid_viscosity_at_temperature(options,p%temperature)
      constant_options%liquid_viscosity_coefficients=0;constant_options%liquid_viscosity_temperature_range=0
      q=p
      call advance_spray_breakup(p,liquid,options,20.0_dp,[200.0_dp,0.0_dp,0.0_dp],20.0_dp,1e-5_dp,1e-5_dp,ok)
      call require(ok,'temperature-dependent viscosity breakup')
      call advance_spray_breakup(q,liquid,constant_options,20.0_dp,[200.0_dp,0.0_dp,0.0_dp], &
        20.0_dp,1e-5_dp,1e-5_dp,ok)
      call require(ok.and.p%mass==q%mass.and.p%deformation==q%deformation.and.p%multiplicity==q%multiplicity, &
        'breakup uses current liquid viscosity')
    end do
    options%liquid_viscosity_temperature_range=[250.0_dp,300.0_dp]
    q=p
    call advance_spray_breakup(p,liquid,options,20.0_dp,[200.0_dp,0.0_dp,0.0_dp],20.0_dp,1e-5_dp,1e-5_dp,ok)
    call require(.not.ok.and.p%mass==q%mass.and.p%rt_time==q%rt_time,'viscosity range rollback')
    options%liquid_viscosity_temperature_range=[250.0_dp,600.0_dp]
    options%liquid_viscosity_coefficients=-1
    call require(liquid_viscosity_at_temperature(options,350.0_dp)==0,'negative viscosity rejected')
    options%liquid_viscosity_coefficients=ieee_value(0.0_dp,ieee_quiet_nan)
    call require(.not.valid_spray_options(options),'NaN viscosity rejected')
  case('density')
    liquid%density_model='polynomial'
    liquid%density_coefficients=[1200.0_dp,-.5_dp,0.0_dp,0.0_dp]
    liquid%density_temperature_range=[250.0_dp,600.0_dp]
    call require(liquid_density_at_temperature(liquid,350.0_dp)==1025.0_dp,'density polynomial')
    d=parcel_diameter(p,liquid)
    options%drag=.false.;options%evaporation=.false.
    call advance_parcel_frozen_film(p,liquid,film,options,1e-3_dp,q,transfer,ok)
    call require(ok.and.q%temperature>p%temperature.and.q%mass==p%mass,'thermal expansion preserves mass')
    expected=d*(1025.0_dp/liquid_density_at_temperature(liquid,q%temperature))**(1.0_dp/3.0_dp)
    call require(abs(parcel_diameter(q,liquid)/expected-1)<1e-13_dp,'heated radius uses updated density')
    call require(maxval(abs(transfer-(parcel_totals(p,liquid)-parcel_totals(q,liquid))))<1e-18_dp, &
      'thermal expansion energy exchange')
    call inject_cone_parcels(liquid,1e-9_dp,2,1_int64,[.001_dp,.001_dp,.001_dp], &
      [0.0_dp,0.0_dp,1.0_dp],0.0_dp,0.0_dp,1e-4_dp,350.0_dp,drops,ok)
    call require(ok.and.abs(drops(1)%mass/(pi*1025.0_dp*1e-12_dp/6)-1)<1e-13_dp,'injected volume')
    options%breakup_model='etab';options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
    constant_liquid=liquid;constant_liquid%density_model='constant';constant_liquid%density=1025
    constant_liquid%density_coefficients=0;constant_liquid%density_temperature_range=0
    q=p
    call advance_etab_breakup(p,liquid,options,20.0_dp,[200.0_dp,0.0_dp,0.0_dp],1e-5_dp,ok)
    call require(ok,'polynomial ETAB')
    call advance_etab_breakup(q,constant_liquid,options,20.0_dp,[200.0_dp,0.0_dp,0.0_dp],1e-5_dp,ok)
    call require(ok.and.p%mass==q%mass.and.p%deformation==q%deformation,'ETAB density parity')
    options%breakup_model='khrt'
    call advance_khrt_breakup(p,liquid,options,20.0_dp,[200.0_dp,0.0_dp,0.0_dp],100.0_dp,1e-5_dp,ok)
    call require(ok,'polynomial KHRT')
    call advance_khrt_breakup(q,constant_liquid,options,20.0_dp,[200.0_dp,0.0_dp,0.0_dp],100.0_dp,1e-5_dp,ok)
    call require(ok.and.p%mass==q%mass.and.p%multiplicity==q%multiplicity,'KHRT density parity')
    liquid%density_temperature_range=[250.0_dp,350.0_dp]
    options%breakup_model='none';options%evaporation=.false.
    call advance_parcel_frozen_film(p,liquid,film,options,1e-3_dp,q,transfer,ok)
    call require(.not.ok.and.q%temperature==p%temperature.and.all(transfer==0),'density-range rollback')
    liquid%density_coefficients=[1.0_dp,-.5_dp,0.0_dp,0.0_dp]
    call require(liquid_density_at_temperature(liquid,300.0_dp)==0,'negative density rejected')
  case('antoine')
    liquid%vapor_pressure_model='antoine'
    liquid%antoine_coefficients=[5.40221_dp,1838.675_dp,-31.737_dp,1e5_dp]
    liquid%antoine_temperature_range=[273.0_dp,303.0_dp]
    do n=0,6
      td=273.0_dp+5*n
      call saturation_mole_fraction(liquid,18.015_dp,td,101325.0_dp,value,ok)
      expected=1e5_dp*10.0_dp**(5.40221_dp-1838.675_dp/(td-31.737_dp))/101325.0_dp
      call require(ok.and.abs(value/expected-1)<1e-13_dp,'Antoine pressure units and formula')
    end do
    call saturation_mole_fraction(liquid,18.015_dp,304.0_dp,101325.0_dp,value,ok)
    call require(.not.ok,'Antoine range rejection')
    call saturation_mole_fraction(liquid,18.015_dp,290.0_dp,1.0_dp,value,ok)
    call require(.not.ok,'Antoine boiling rejection')
    liquid%antoine_coefficients(3)=-280.0_dp
    call require(.not.valid_spray_liquid(liquid),'Antoine pole inside range rejected')
    liquid%antoine_coefficients(3)=-31.737_dp
    liquid%antoine_coefficients(4)=0
    call require(.not.valid_spray_liquid(liquid),'Antoine zero pressure scale rejected')
    liquid%antoine_coefficients(4)=ieee_value(0.0_dp,ieee_quiet_nan)
    call require(.not.valid_spray_liquid(liquid),'Antoine NaN rejected')
  case('split')
    call make_drops();saved_drops=drops
    drops%multiplicity=10.0_dp
    drops%khrt_initial_count=7.0_dp;drops%khrt_shed_volume=1e-12_dp
    transfer=0
    do n=1,size(drops)
      transfer=transfer+parcel_totals(drops(n),liquid)
    end do
    call split_spray_parcels(drops,3.0_dp,100,ok)
    call require(ok.and.size(drops)==16,'child allocation')
    call require(all(drops%multiplicity<=3.0_dp),'multiplicity cap')
    call require(all(drops(1:4)%parent_id==0).and.all(drops(5:)%parent_id>0),'parent lineage')
    do n=1,size(drops)
      call require(drops(n)%id==int(n,int64),'stable ordered IDs')
      transfer=transfer-parcel_totals(drops(n),liquid)
    end do
    call require(maxval(abs(transfer))<1e-18_dp,'split conserved inventory')
    saved_drops=drops
    call split_spray_parcels(drops,1.0_dp,17,ok)
    call require(.not.ok.and.size(drops)==size(saved_drops).and.all(drops%id==saved_drops%id), &
      'storage failure preserves IDs and allocation')
    call require(all(drops%multiplicity==saved_drops%multiplicity),'storage failure preserves weights')
    call require(all(drops%khrt_initial_count==saved_drops%khrt_initial_count).and. &
      all(drops%khrt_shed_volume==saved_drops%khrt_shed_volume),'storage failure preserves KH inventory')
    call split_spray_parcels(drops,1.0_dp,100,ok)
    call require(ok.and.size(drops)==48,'repeat split with rounding remainder')
    call require(abs(sum(drops%khrt_initial_count)-28.0_dp)<1e-12_dp,'split conserves KH initial count')
    call require(abs(sum(drops%khrt_shed_volume)/4e-12_dp-1)<1e-14_dp,'split conserves pending shed volume')
    do n=1,size(drops)
      call require(abs(drops(n)%khrt_initial_count/drops(n)%multiplicity-.7_dp)<1e-14_dp,'KH count per drop')
      call require(abs(drops(n)%khrt_shed_volume/drops(n)%multiplicity/1e-13_dp-1)<1e-14_dp,'KH volume per drop')
    end do
    allocate(child_component_mass(2,size(drops)))
    child_component_mass(1,:)=.25_dp*drops%mass
    child_component_mass(2,:)=.75_dp*drops%mass
    saved_child_mass=child_component_mass;saved_drops=drops
    split_watermark=huge(split_watermark)
    call split_spray_parcels(drops,.5_dp,100,ok,split_watermark,child_component_mass)
    call require(.not.ok.and.split_watermark==huge(split_watermark),'split ID overflow rejection')
    call require(all(drops%id==saved_drops%id).and.all(child_component_mass==saved_child_mass),'overflow rollback')
    split_watermark=1000
    child_component_mass(1,1)=-1
    call split_spray_parcels(drops,.5_dp,100,ok,split_watermark,child_component_mass)
    call require(.not.ok.and.split_watermark==1000.and.size(drops)==48,'invalid split components rejected')
    child_component_mass=saved_child_mass
    call split_spray_parcels(drops,.5_dp,100,ok,split_watermark,child_component_mass)
    call require(ok.and.size(drops)==96.and.split_watermark==1048,'mixed split shared watermark')
    do n=49,96
      call require(drops(n)%id==1000+n-48,'mixed split avoids retired IDs')
      call require(all(child_component_mass(:,n)==saved_child_mass(:,int(drops(n)%parent_id))), &
        'mixed split component ancestry')
    end do
  case('khrt')
    options%breakup_model='khrt';options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
    call make_drops();p=drops(1);p%velocity=0;q=p
    dt=1e-7_dp;d=.5_dp*parcel_diameter(p,liquid)
    call advance_khrt_breakup(p,liquid,options,10.0_dp,[0.0_dp,0.0_dp,0.0_dp],0.0_dp,dt,ok)
    call require(ok.and.p%mass==q%mass.and.p%rt_time==-1,'zero slip unchanged')
    value=10*300.0_dp**2*d/options%surface_tension
    z=options%liquid_viscosity/sqrt(liquid%density*d*options%surface_tension);tkh=z*sqrt(value)
    wavelength=d*9.02_dp*(1+.45_dp*sqrt(z))*(1+.4_dp*tkh**.7_dp)/(1+.865_dp*value**1.67_dp)**.6_dp
    omega=(.34_dp+.385_dp*value**1.5_dp)/((1+z)*(1+1.4_dp*tkh**.6_dp))* &
      sqrt(options%surface_tension/(liquid%density*d**3))
    tau=3.726_dp*options%khrt_b1*d/(wavelength*omega)
    expected=((dt/tau*options%khrt_b0*wavelength+d)/(1+dt/tau)/d)**3
    call advance_khrt_breakup(p,liquid,options,10.0_dp,[300.0_dp,0.0_dp,0.0_dp],2000.0_dp,dt,ok)
    call require(ok.and.p%rt_time==0,'KH enables RT timer')
    call near(p%mass/q%mass,expected,1e-13_dp,'KH analytical mass')
    call near(p%mass*p%multiplicity,q%mass*q%multiplicity,1e-14_dp,'KH inventory')
    p=q;p%rt_time=1e-3_dp;options%khrt_c3=.1_dp
    acceleration=.375_dp*.424_dp*10*300.0_dp**2/(liquid%density*d)
    wavelength=options%khrt_c3*pi/sqrt(acceleration*(liquid%density-10)/(3*options%surface_tension))
    call require(wavelength<d,'RT unstable wavelength')
    call advance_khrt_breakup(p,liquid,options,10.0_dp,[300.0_dp,0.0_dp,0.0_dp],2000.0_dp,dt,ok)
    call require(ok.and.p%rt_time==0,'RT event resets timer')
    call near(p%mass/q%mass,(wavelength/d)**3,1e-13_dp,'RT analytical mass')
    call near(p%mass*p%multiplicity,q%mass*q%multiplicity,1e-14_dp,'RT inventory')
    call make_drops();p=drops(1);p%velocity=0;q=p;options%khrt_c3=1
    shedding=spray_khrt_shedding_state()
    call advance_khrt_shedding(p,liquid,options,10.0_dp,[300.0_dp,0.0_dp,0.0_dp],2000.0_dp,dt, &
      1e24_dp,shedding,release,ok)
    call require(ok.and.release%mass==0.and.shedding%shed_volume>0,'KH pending shedding')
    call near(p%mass*p%multiplicity,q%mass*q%multiplicity,1e-13_dp,'pending mass stays in parent')
    p=q;shedding=spray_khrt_shedding_state()
    call advance_khrt_shedding(p,liquid,options,10.0_dp,[300.0_dp,0.0_dp,0.0_dp],2000.0_dp,dt, &
      1e-25_dp,shedding,release,ok)
    call require(ok.and.release%mass>0.and.p%mass>0,'KH partial release')
    call require(shedding%shed_volume==0.and.shedding%initial_count==p%multiplicity,'KH resets shed accumulator')
    call near(p%mass*p%multiplicity+release%mass,q%mass*q%multiplicity,1e-13_dp,'partial shedding mass')
    call construct_khrt_children(p,release,[1.0_dp,0.0_dp,0.0_dp],2_int64,saved_drops,solid_energy,ok)
    call require(ok.and.size(saved_drops)==4,'KH child construction')
    child_totals=parcel_totals(p,liquid)-parcel_totals(q,liquid)
    do n=1,4
      child_totals=child_totals+parcel_totals(saved_drops(n),liquid)
      call require(saved_drops(n)%parent_id==q%id,'KH child lineage')
    end do
    child_totals(5)=child_totals(5)+solid_energy
    call require(maxval(abs(child_totals))<1e-12_dp,'KH gas parent children conservation')
    p=q;shedding=spray_khrt_shedding_state()
    call advance_khrt_shedding(p,liquid,options,10.0_dp,[300.0_dp,0.0_dp,0.0_dp],2000.0_dp,1.0_dp, &
      1e24_dp,shedding,release,ok)
    call require(ok.and.p%mass==0,'KH full parent release')
    call near(release%mass,q%mass*q%multiplicity,1e-13_dp,'full shedding mass')
    saved_shedding=shedding;p=q
    call advance_khrt_shedding(p,liquid,options,10.0_dp,[300.0_dp,0.0_dp,0.0_dp],2000.0_dp,dt, &
      -1.0_dp,shedding,release,ok)
    call require(.not.ok.and.p%mass==q%mass.and.release%mass==0.and. &
      shedding%shed_volume==saved_shedding%shed_volume,'invalid shedding rollback')
    q=p;options%khrt_b1=-1
    call advance_khrt_breakup(p,liquid,options,10.0_dp,[300.0_dp,0.0_dp,0.0_dp],2000.0_dp,dt,ok)
    call require(.not.ok.and.p%mass==q%mass.and.p%rt_time==q%rt_time,'invalid KH coefficient rollback')
  case('etab')
    options%breakup_model='etab';options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
    call make_drops(); p=drops(1); p%velocity=0; p%deformation=.1_dp;p%deformation_rate=100
    q=p; d=.5_dp*parcel_diameter(p,liquid); dt=1e-7_dp
    td=liquid%density*d*d/(5*options%liquid_viscosity)
    omega=sqrt(8*options%surface_tension/(liquid%density*d**3)-1/td**2)
    wer=d*100/options%surface_tension/12
    expected=wer+exp(-dt/td)*((p%deformation-wer)*cos(omega*dt)+ &
      (p%deformation_rate+(p%deformation-wer)/td)/omega*sin(omega*dt))
    call advance_etab_breakup(p,liquid,options,1.0_dp,[10.0_dp,0.0_dp,0.0_dp],dt,ok)
    call require(ok,'ETAB oscillator')
    call near(p%deformation,expected,1e-13_dp,'damped forced oscillator')
    call require(p%mass==q%mass.and.p%multiplicity==q%multiplicity,'subcritical no breakup')
    p=q;p%deformation=1.1_dp;p%deformation_rate=0
    value=10.0_dp*300.0_dp**2*d/options%surface_tension
    expected=exp(-(2.0_dp/9.0_dp)*sqrt(value)*acos(1.0_dp-12.0_dp/value))
    call advance_etab_breakup(p,liquid,options,10.0_dp,[300.0_dp,0.0_dp,0.0_dp],1e-10_dp,ok)
    call require(ok,'instant ETAB event')
    call near(p%mass/q%mass,expected**3,1e-12_dp,'ETAB analytical child mass')
    p=q;p%deformation=0;p%deformation_rate=0
    transfer=parcel_totals(p,liquid);initial_mass=p%mass
    call advance_etab_breakup(p,liquid,options,10.0_dp,[300.0_dp,0.0_dp,0.0_dp],5e-5_dp,ok)
    call require(ok.and.p%mass<initial_mass.and.p%multiplicity>q%multiplicity,'supercritical ETAB breakup')
    call require(maxval(abs(parcel_totals(p,liquid)-transfer))<1e-13_dp*maxval(abs(transfer)), &
      'ETAB inventory conservation')
    q=p;options%maximum_substeps=1
    call advance_etab_breakup(p,liquid,options,10.0_dp,[300.0_dp,0.0_dp,0.0_dp],1.0_dp,ok)
    call require(.not.ok.and.p%mass==q%mass.and.p%deformation==q%deformation,'ETAB substep rollback')
  case('boundaries')
    faces='slip_wall'
    call particle_boundary_stencil([0.0_dp,.002_dp,.002_dp],[0.0_dp,0.0_dp,0.0_dp], &
      [.001_dp,.001_dp,.001_dp],[4,4,4],faces,indices,weights,ok)
    call require(ok.and.abs(sum(weights)-1.0_dp)<1e-14_dp,'physical stencil partition')
    do n=1,8
      if (weights(n)>0) call require(indices(1,n)==1,'no opposite-face interpolation')
    end do
    p%id=1; p%position=[.0041_dp,.002_dp,.002_dp]; p%velocity=[10.0_dp,2.0_dp,3.0_dp]
    p%mass=1e-12_dp; p%temperature=330; p%multiplicity=10
    q=p
    call apply_parcel_boundaries(p,liquid,[0.0_dp,0.0_dp,0.0_dp],[.004_dp,.004_dp,.004_dp],faces,transfer,ok)
    call require(ok.and.abs(p%position(1)-.0039_dp)<1e-15_dp.and.p%velocity(1)==-10,'specular wall')
    call require(abs(transfer(2)+2e-10_dp)<1e-24_dp.and.transfer(1)==0.and.transfer(5)==0,'wall impulse and energy')
    faces='outflow'; p=q
    call apply_parcel_boundaries(p,liquid,[0.0_dp,0.0_dp,0.0_dp],[.004_dp,.004_dp,.004_dp],faces,transfer,ok)
    call require(ok.and.p%mass==0,'parcel outflow removal')
    call require(maxval(abs(transfer+parcel_totals(q,liquid)))<1e-24_dp,'outflow inventory')
  case('drag')
    options%heat=.false.; options%evaporation=.false.
    call spray_transfer_rates(p,liquid,film,options,inverse_tau,heat,masscoef,ok)
    call require(ok,'drag rates')
    d=parcel_diameter(p,liquid)
    expected=18.0_dp*film%viscosity/(liquid%density*d*d) * &
      (1.0_dp+0.15_dp*(film%density*d/film%viscosity)**0.687_dp)
    call near(inverse_tau,expected,1.0e-13_dp,'Schiller-Naumann coefficient')
    dt=0.01_dp
    call advance_parcel_frozen_film(p,liquid,film,options,dt,q,transfer,ok)
    call require(ok,'drag step')
    call near(q%velocity(1),exp(-expected*dt),1.0e-13_dp,'exponential drag')
    call require(q%mass == p%mass .and. q%temperature == p%temperature,'drag leaves mass and T')
    call exchange_check(p,q,transfer)
    options%drag_model='pelephysics'
    do n=1,4
      select case(n)
      case(1)
        value=.5_dp;expected=24
      case(2)
        value=100.0_dp;expected=24*(1+100.0_dp**(2.0_dp/3)/6)
      case(3)
        value=2000.0_dp;expected=.424_dp*2000
      case(4)
        value=0;expected=0
      end select
      film%gas_velocity=p%velocity+[value*film%viscosity/(film%density*d),0.0_dp,0.0_dp]
      expected=3*film%viscosity*expected/(4*liquid%density*d*d)
      call spray_transfer_rates(p,liquid,film,options,inverse_tau,heat,masscoef,ok)
      call require(ok,'PelePhysics drag regime')
      call near(inverse_tau,expected,1e-13_dp,'PelePhysics drag rate')
    end do
    film%gas_velocity=0;options%breakup_model='etab';options%surface_tension=.072_dp
    options%liquid_viscosity=.001_dp;p%deformation=.2_dp
    call spray_transfer_rates(p,liquid,film,options,inverse_tau,heat,masscoef,ok)
    call require(ok,'deformed drag')
    p%deformation=0
    call spray_transfer_rates(p,liquid,film,options,value,heat,masscoef,ok)
    call require(ok,'spherical drag')
    call near(inverse_tau,value*(1+2.6432_dp*.2_dp),1e-13_dp,'ETAB drag enhancement')
    p%deformation=-1
    call advance_parcel_frozen_film(p,liquid,film,options,dt,q,transfer,ok)
    call require(.not.ok.and.all(q%velocity==p%velocity).and.all(transfer==0),'negative drag rollback')
  case('heat')
    options%drag=.false.; options%evaporation=.false.; p%velocity=0.0_dp
    call spray_transfer_rates(p,liquid,film,options,inverse_tau,heat,masscoef,ok)
    call require(ok,'heat rates')
    expected=2.0_dp*pi*film%conductivity*parcel_diameter(p,liquid)
    call near(heat,expected,1.0e-13_dp,'Nu=2 stagnant limit')
    dt=0.01_dp
    expected=film%gas_temperature+(p%temperature-film%gas_temperature)*exp(-heat*dt/(p%mass*liquid%heat_capacity))
    call advance_parcel_frozen_film(p,liquid,film,options,dt,q,transfer,ok)
    call require(ok,'heat step')
    call near(q%temperature,expected,1.0e-13_dp,'analytic thermal relaxation')
    call exchange_check(p,q,transfer)
  case('evaporation')
    p%velocity=0.0_dp; options%drag=.false.; options%heat=.false.
    call spray_transfer_rates(p,liquid,film,options,inverse_tau,heat,masscoef,ok)
    call require(ok,'evaporation rates')
    expected=2.0_dp*pi*film%density*film%diffusivity*log(1.0_dp/0.8_dp)
    call near(masscoef,expected,1.0e-13_dp,'Stefan mass transfer')
    dt=1.0e-5_dp; d=parcel_diameter(p,liquid)
    expected=d*d-4.0_dp*masscoef*dt/(pi*liquid%density)
    call advance_parcel_frozen_film(p,liquid,film,options,dt,q,transfer,ok)
    call require(ok,'evaporation step')
    call near(parcel_diameter(q,liquid)**2,expected,1.0e-13_dp,'frozen-film d-squared law')
    call require(q%mass < p%mass .and. q%temperature < p%temperature,'latent cooling')
    call exchange_check(p,q,transfer)
    film%ambient_fraction=0.3_dp
    call spray_transfer_rates(p,liquid,film,options,inverse_tau,heat,masscoef,ok)
    call require(ok .and. masscoef == 0.0_dp,'supersaturation suppresses evaporation')
  case('terminal')
    options%drag=.false.; options%heat=.false.
    call spray_transfer_rates(p,liquid,film,options,inverse_tau,heat,masscoef,ok)
    dt=2.0_dp*parcel_diameter(p,liquid)**2*pi*liquid%density/(4.0_dp*masscoef)
    call advance_parcel_frozen_film(p,liquid,film,options,dt,q,transfer,ok)
    call require(ok .and. q%mass == 0.0_dp,'complete evaporation')
    call exchange_check(p,q,transfer)
    call near(transfer(1),p%mass*p%multiplicity,1.0e-14_dp,'no discarded residue')
  case('abramzon')
    options%abramzon_sirignano=.true.
    call spray_transfer_rates(p,liquid,film,options,inverse_tau,heat,masscoef,ok)
    call require(ok .and. heat > 0.0_dp .and. masscoef > 0.0_dp,'AS converges')
    film%surface_fraction=film%ambient_fraction; p%velocity=0.0_dp
    call spray_transfer_rates(p,liquid,film,options,inverse_tau,heat,masscoef,ok)
    call require(ok .and. masscoef == 0.0_dp,'AS zero-Spalding limit')
    call near(heat,2.0_dp*pi*film%conductivity*parcel_diameter(p,liquid),1.0e-13_dp,'AS heat limit')
  case('film_binary')
    primitive=0.0_dp; primitive(1)=1.0_dp; primitive(5)=101325.0_dp
    primitive(11)=0.1_dp; primitive(15)=0.9_dp
    call build_spray_film(species,transport,primitive,500.0_dp,p,liquid,options,film,ok)
    call require(ok,'binary film')
    call binary_diffusion_coefficients(species,transport,400.0_dp,101325.0_dp,binary,ok)
    call require(ok,'independent binary coefficient')
    ! In a water/nitrogen binary mixture the Stefan coefficient must be D(H2O,N2),
    ! not a molecule-weight-scaled variant of Cantera/PelePhysics conventions.
    call near(film%diffusivity,binary(6,10),1.0e-13_dp,'binary diffusion limit')
    call build_spray_film(species,transport,primitive,500.0_dp,p,liquid,options,film,ok,wall_skin=.true.)
    call require(ok,'wall skin properties')
    call binary_diffusion_coefficients(species,transport,425.0_dp,101325.0_dp,binary,ok)
    call require(ok,'wall midpoint binary coefficient')
    call near(film%diffusivity,binary(6,10),1.0e-13_dp,'wall half-temperature rule')
    skin_y=0;skin_y(6)=.5_dp*(film%surface_fraction+.1_dp);skin_y(10)=1-skin_y(6)
    expected=101325.0_dp/(mixture_specific_gas_constant(species,skin_y,ok)*425.0_dp)
    call require(ok,'wall mixture gas constant')
    call near(film%density,expected,1e-13_dp,'wall half-composition rule')
    call build_spray_film(species,transport,primitive,300.0_dp,p,liquid,options,film,ok,wall_skin=.true.)
    call require(ok,'cold gas wall skin')
    call binary_diffusion_coefficients(species,transport,350.0_dp,101325.0_dp,binary,ok)
    call require(ok,'cold gas binary coefficient')
    call near(film%diffusivity,binary(6,10),1e-13_dp,'wall skin clamps negative gas temperature difference')
  case('injection')
    call make_drops()
    expected=sum(drops%mass*drops%multiplicity)
    call near(expected,1.0e-10_dp,1.0e-14_dp,'injected mass')
    do n=1,size(drops)
      call require(drops(n)%id == int(n,int64),'unique IDs')
      call near(norm2(drops(n)%velocity),10.0_dp,1.0e-13_dp,'injection speed')
      call require(drops(n)%velocity(3)/10.0_dp >= cos(0.1_dp),'cone angle')
    end do
  case('size_distribution')
    call inject_cone_parcels(liquid,1.0e-8_dp,64,1_int64,[0.0_dp,0.0_dp,0.0_dp], &
      [0.0_dp,0.0_dp,1.0_dp],10.0_dp,0.1_dp,20.0e-6_dp,330.0_dp,drops,ok, &
      diameter_shape=3.0_dp,diameter_min=10.0e-6_dp,diameter_max=40.0e-6_dp)
    call require(ok,'truncated mass-CDF injection')
    value=0.0_dp
    do n=1,size(drops)
      value=value+drops(n)%mass*drops(n)%multiplicity
      d=parcel_diameter(drops(n),liquid)
      call require(d>10.0e-6_dp .and. d<40.0e-6_dp,'bounded size quantiles')
      expected=(exp(-0.5_dp**3)-exp(-(d/20.0e-6_dp)**3))/(exp(-0.5_dp**3)-exp(-2.0_dp**3))
      call near(expected,(real(n,dp)-0.5_dp)/64.0_dp,1.0e-12_dp,'inverse mass CDF')
      call near(drops(n)%mass*drops(n)%multiplicity,1.0e-8_dp/64.0_dp,1.0e-13_dp,'equal parcel mass')
    end do
    call near(value,1.0e-8_dp,1.0e-13_dp,'distributed injected mass')
    deallocate(drops)
    call inject_cone_parcels(liquid,1.0e-8_dp,64,1_int64,[0.0_dp,0.0_dp,0.0_dp], &
      [0.0_dp,0.0_dp,1.0_dp],10.0_dp,0.1_dp,20.0e-6_dp,330.0_dp,drops,ok, &
      diameter_shape=3.0_dp,diameter_min=10.0e-6_dp,diameter_max=40.0e-6_dp,diameter_weighting='number')
    call require(ok,'number-CDF injection')
    call require(all(drops%multiplicity==drops(1)%multiplicity),'equal represented droplet counts')
    call require(drops(64)%mass>drops(1)%mass,'number-CDF parcel masses differ')
    call near(sum(drops%mass*drops%multiplicity),1e-8_dp,1e-13_dp,'number-CDF total mass')
    do n=1,64
      d=parcel_diameter(drops(n),liquid)
      expected=(exp(-.5_dp**3)-exp(-(d/20e-6_dp)**3))/(exp(-.5_dp**3)-exp(-2.0_dp**3))
      call near(expected,(real(n,dp)-.5_dp)/64,1e-12_dp,'inverse number CDF')
    end do
    deallocate(drops)
    call inject_cone_parcels(liquid,1.0e-8_dp,2,1_int64,[0.0_dp,0.0_dp,0.0_dp], &
      [0.0_dp,0.0_dp,1.0_dp],10.0_dp,0.1_dp,20.0e-6_dp,330.0_dp,drops,ok,diameter_shape=-1.0_dp)
    call require(.not. ok .and. .not. allocated(drops),'negative size shape rejected')
  case('stencil')
    call periodic_particle_stencil([0.0015_dp,0.0015_dp,0.0015_dp],lower,spacing,[4,4,4],indices,weights,ok)
    call require(ok,'CIC stencil')
    call near(maxval(weights),1.0_dp,1.0e-14_dp,'cell-center weight')
    call periodic_particle_stencil([-0.00025_dp,0.00425_dp,0.0_dp],lower,spacing,[4,4,4],indices,weights,ok)
    call require(ok .and. all(indices >= 1) .and. all(indices <= 4),'periodic corner')
    call near(sum(weights),1.0_dp,1.0e-14_dp,'partition of unity')
  case('coupling','rollback','substeps','zero')
    call make_drops(); call make_gas(.false.)
    call spray_system_totals(state,spacing,drops,liquid,totals0)
    saved=state; saved_t=temp; saved_drops=drops
    dt=1.0e-6_dp
    if (test == 'rollback') drops(4)%temperature=10000.0_dp
    if (test == 'rollback') saved_drops=drops
    if (test == 'substeps') then
      options%maximum_substeps=1; dt=1.0e-3_dp
    end if
    if (test == 'zero') dt=0.0_dp
    call advance_spray_coupling_3d(species,transport,liquid,options,state,temp,lower,spacing,drops,dt,steps,ok,message)
    if (test == 'rollback' .or. test == 'substeps') then
      call require(.not. ok,'reject invalid batch')
      call require(all(state == saved) .and. all(temp == saved_t),'gas rollback')
      call require(all(drops%mass == saved_drops%mass) .and. &
        all(drops%temperature == saved_drops%temperature),'parcel rollback')
      call require(steps == 0,'no committed substeps on failure')
    else
      if (.not. ok) print *,trim(message)
      call require(ok,'coupled update')
      call spray_system_totals(state,spacing,drops,liquid,totals1)
      call near(totals1(1),totals0(1),1.0e-12_dp,'two-phase mass')
      call near(totals1(5),totals0(5),1.0e-12_dp,'two-phase total energy')
      call near(totals1(11),totals0(11),1.0e-12_dp,'water mass')
      do n=2,4
        call require(abs(totals1(n)-totals0(n)) < 1.0e-21_dp,'two-phase momentum')
      end do
      if (test == 'zero') then
        call require(all(state == saved) .and. all(temp == saved_t),'zero interval identity')
      else
        call require(sum(drops%mass) < sum(saved_drops%mass),'evaporation active in grid')
        call require(sum(state(11,:,:,:)) > sum(saved(11,:,:,:)),'vapor deposited')
      end if
    end if
  case('invalid')
    p%mass=ieee_value(0.0_dp,ieee_quiet_nan)
    call advance_parcel_frozen_film(p,liquid,film,options,1.0_dp,q,transfer,ok)
    call require(.not. ok .and. all(transfer == 0.0_dp),'reject NaN parcel')
    call saturation_mole_fraction(liquid,18.015_dp,500.0_dp,101325.0_dp,expected,ok)
    call require(.not. ok,'boiling explicitly rejected')
  case('les_algebra')
    grad=0.0_dp; grad(1,2)=3.0_dp
    call smagorinsky_stress(grad,2.0_dp,[0.1_dp,0.1_dp,0.1_dp],0.2_dp,mu,stress,ok)
    call require(ok,'Smagorinsky')
    call near(mu,2.0_dp*(0.2_dp*0.1_dp)**2*3.0_dp,1.0e-13_dp,'simple shear eddy viscosity')
    call near(stress(1,2),mu*3.0_dp,1.0e-13_dp,'shear stress')
    grad(2,1)=-3.0_dp
    call smagorinsky_stress(grad,2.0_dp,[0.1_dp,0.1_dp,0.1_dp],0.2_dp,mu,stress,ok)
    call require(ok .and. mu == 0.0_dp .and. all(stress == 0.0_dp),'rigid rotation zero stress')
  case('les_dynamic_flux','les_dynamic_step')
    block
      real(dp) :: other(15,4,4,4),reference(15,4,4,4),sound,power,error_coarse,error_fine
      integer :: i,j,k,c,substep
      call make_gas(.false.)
      call compute_dynamic_les_rhs_3d(species,state,temp,spacing,'box',2,rhs,ok)
      call require(ok.and.all(rhs==0),'uniform dynamic flux zero')
      do k=1,4
        do j=1,4
          do i=1,4
            primitive=0;primitive(1)=.3_dp
            primitive(2)=20*sin(2*pi*(j-1)/4);primitive(3)=15*cos(2*pi*(i-1)/4)
            primitive(5)=101325;primitive(15)=1
            if (test=='les_dynamic_step') then
              primitive(1)=.3_dp*(1+.1_dp*sin(2*pi*(i-1)/4+.2_dp))
              primitive(4)=5*sin(2*pi*(k-1)/4+.3_dp)
              primitive(5)=101325*(1+.05_dp*cos(2*pi*(j-1)/4+.4_dp))
            end if
            call reactive_primitive_to_conserved(species,primitive,state(:,i,j,k),temp(i,j,k),sound,ok)
            call require(ok,'dynamic flux fixture')
          end do
        end do
      end do
      saved=state;saved_t=temp
      if (test=='les_dynamic_step') then
        call compute_dynamic_les_rhs_3d(species,state,temp,spacing,'box',2,rhs,ok,maximum_diffusivity=diff)
        call require(ok.and.diff>0,'dynamic diffusion estimate')
        dt=.05_dp/(diff*sum(1/spacing**2))
        call advance_dynamic_les_ssprk2_3d(species,state,temp,spacing,'box',2,0._dp,ok)
        call require(ok.and.all(state==saved).and.all(temp==saved_t),'dynamic zero step')
        ke0=sum(.5_dp*sum(state(2:4,:,:,:)**2,dim=1)/state(1,:,:,:))
        call advance_dynamic_les_ssprk2_3d(species,state,temp,spacing,'box',2,dt,ok)
        call require(ok,'dynamic SSPRK2 step')
        ke1=sum(.5_dp*sum(state(2:4,:,:,:)**2,dim=1)/state(1,:,:,:))
        call require(ieee_is_finite(ke1).and.ke1/=ke0,'dynamic step changes kinetic energy')
        do c=1,15
          call require(abs(sum(state(c,:,:,:))-sum(saved(c,:,:,:))) &
            <=1e-12_dp*max(1._dp,sum(abs(saved(c,:,:,:)))),'dynamic step conservation')
        end do
        call require(all(state(1,:,:,:)==saved(1,:,:,:)).and.all(state(6:,:,:,:)==saved(6:,:,:,:)), &
          'dynamic step preserves species and density')
        other=state;state=saved;temp=saved_t
        do substep=1,2
          call advance_dynamic_les_ssprk2_3d(species,state,temp,spacing,'box',2,dt/2,ok)
          call require(ok,'dynamic half step')
        end do
        reference=state;state=saved;temp=saved_t
        do substep=1,4
          call advance_dynamic_les_ssprk2_3d(species,state,temp,spacing,'box',2,dt/4,ok)
          call require(ok,'dynamic quarter step')
        end do
        error_coarse=maxval(abs(other-reference));error_fine=maxval(abs(reference-state))
        call require(error_coarse>3*error_fine.and.error_coarse<5*error_fine,'dynamic second order refinement')
        state=saved;temp=saved_t
        call advance_dynamic_les_ssprk2_3d(species,state,temp,spacing,'box',2,10*dt,ok)
        call require(.not.ok.and.all(state==saved).and.all(temp==saved_t),'dynamic oversized step rollback')
        call advance_dynamic_les_ssprk2_3d(species,state,temp,spacing,'box',3,0._dp,ok)
        call require(.not.ok.and.all(state==saved).and.all(temp==saved_t),'dynamic invalid model at zero dt')
      end if
      call compute_dynamic_les_rhs_3d(species,state,temp,spacing,'box',2,rhs,ok)
      call require(ok.and.maxval(abs(rhs(2:5,:,:,:)))>0,'dynamic flux active')
      call require(all(rhs(1,:,:,:)==0).and.all(rhs(6:,:,:,:)==0),'dynamic mass and species unchanged')
      do c=1,15
        call require(abs(sum(rhs(c,:,:,:)))<=1e-12_dp*max(1._dp,sum(abs(rhs(c,:,:,:)))), &
          'dynamic periodic flux conservation')
      end do
      power=sum(state(2:4,:,:,:)*rhs(2:4,:,:,:)/spread(state(1,:,:,:),1,3))
      if (test=='les_dynamic_flux') call require(power<0,'dynamic shear kinetic energy dissipation')
      call require(all(state==saved).and.all(temp==saved_t),'dynamic flux inputs preserved')
      state=cshift(saved,1,dim=2);temp=cshift(saved_t,1,dim=1)
      call compute_dynamic_les_rhs_3d(species,state,temp,spacing,'box',2,other,ok)
      reference=cshift(rhs,1,dim=2)
      call require(ok.and.maxval(abs(other-reference))<=1e-12_dp*maxval(abs(rhs)), &
        'dynamic flux translation')
      state=saved;temp=saved_t;reference=rhs;state(1,1,1,1)=-1
      call compute_dynamic_les_rhs_3d(species,state,temp,spacing,'box',2,rhs,ok)
      call require(.not.ok.and.all(rhs==reference),'dynamic flux rollback')
    end block
  case('les_dynamic_field')
    block
      real(dp) :: coeff(3,4,4,4),other(3,4,4,4),reference(3,4,4,4),sound
      integer :: i,j,k
      call make_gas(.false.)
      call compute_dynamic_les_coefficients_3d(species,state,temp,spacing,'box',2,coeff,ok)
      call require(ok.and.maxval(abs(coeff-1e-8_dp))<1e-20_dp,'uniform dynamic field floor')
      do k=1,4
        do j=1,4
          do i=1,4
            primitive=0;primitive(1)=.3_dp*(1+.1_dp*sin(2*pi*(i-1)/4))
            primitive(2)=20*sin(2*pi*(j-1)/4);primitive(3)=15*cos(2*pi*(i-1)/4)
            primitive(4)=5*sin(2*pi*(k-1)/4);primitive(5)=101325*(1+.05_dp*cos(2*pi*(j-1)/4))
            primitive(15)=1
            call reactive_primitive_to_conserved(species,primitive,state(:,i,j,k),temp(i,j,k),sound,ok)
            call require(ok,'dynamic variable-density fixture')
          end do
        end do
      end do
      saved=state;saved_t=temp
      call compute_dynamic_les_coefficients_3d(species,state,temp,spacing,'box',2,coeff,ok)
      call require(ok.and.all(ieee_is_finite(coeff)),'dynamic nonuniform field')
      call require(minval(coeff)>0.and.maxval(coeff(2,:,:,:))<=1,'dynamic coefficient bounds')
      call require(maxval(coeff)>1e-6_dp,'dynamic nontrivial coefficients')
      call require(all(state==saved).and.all(temp==saved_t),'dynamic field inputs preserved')
      state=cshift(saved,1,dim=2);temp=cshift(saved_t,1,dim=1)
      call compute_dynamic_les_coefficients_3d(species,state,temp,spacing,'box',2,other,ok)
      call require(ok.and.all(other==cshift(coeff,1,dim=2)),'dynamic periodic translation')
      state=saved;temp=saved_t
      call compute_dynamic_les_coefficients_3d(species,state,temp,spacing,'gaussian',2,other,ok)
      call require(ok.and.maxval(abs(other-coeff))>1e-8_dp,'dynamic filter selection changes coefficients')
      reference=coeff;state(1,1,1,1)=-1
      call compute_dynamic_les_coefficients_3d(species,state,temp,spacing,'box',2,coeff,ok)
      call require(.not.ok.and.all(coeff==reference),'dynamic invalid state rollback')
      state=saved
      call compute_dynamic_les_coefficients_3d(species,state,temp,spacing,'box',3,coeff,ok)
      call require(.not.ok.and.all(coeff==reference),'dynamic invalid filter rollback')
    end block
  case('les_filter_approximations')
    block
      character(len=16), parameter :: models(8)=[character(len=16)::'box3','box5','gaussian3','gaussian5', &
        'box3opt','box5opt','gaussian3opt','gaussian5opt']
      real(dp), allocatable :: w(:),fallback(:)
      real(dp) :: moment2,moment4,field(1,5,3,2),filtered(1,5,3,2)
      integer :: m,ratio,i,radius
      do m=1,8
        do ratio=1,10
          call les_filter_weights(models(m),ratio,w,ok)
          call require(ok,'approximation weights')
          call require(all(w==w(size(w):1:-1)),'approximation symmetry')
          call near(sum(w),1._dp,1e-12_dp,'approximation normalization')
        end do
        call les_filter_weights(models(m),2,w,ok)
        field=0;field(1,2,2,1)=1
        call filter_les_field_3d(field,w,filtered,ok)
        call require(ok,'signed approximation application')
        call near(sum(filtered),1._dp,1e-12_dp,'signed periodic filter total')
      end do
      do m=1,4
        call les_filter_weights(models(m),3,w,ok)
        moment2=0;moment4=0;radius=size(w)/2
        do i=-radius,radius
          moment2=moment2+w(i+radius+1)*i**2
          moment4=moment4+w(i+radius+1)*i**4
        end do
        call near(moment2,9._dp/12,1e-13_dp,'polynomial second moment')
        if (m==2) call near(moment4,81._dp/80,1e-13_dp,'box fourth moment')
        if (m==4) call near(moment4,81._dp/48,1e-13_dp,'Gaussian fourth moment')
      end do
      call les_filter_weights('box3opt',4,w,ok)
      call near(w(1)/w(2),-2.375_dp,1e-13_dp,'box optimized table ratio')
      call require(w(2)<0,'negative optimized weight retained')
      call les_filter_weights('gaussian3opt',10,w,ok)
      call near(w(1)/w(2),-.581_dp,1e-13_dp,'Gaussian optimized table ratio')
      call les_filter_weights('box5opt',2,w,ok)
      call near(w(1)/w(3),-.013_dp,1e-13_dp,'box optimized outer ratio')
      call near(w(2)/w(3),.3178_dp,1e-13_dp,'box optimized inner ratio')
      call les_filter_weights('gaussian5opt',1,w,ok)
      call near(w(1)/w(3),-.0175_dp,1e-13_dp,'Gaussian optimized outer ratio')
      call near(w(2)/w(3),.0871_dp,1e-13_dp,'Gaussian optimized inner ratio')
      call les_filter_weights('box5opt',12,w,ok)
      call require(ok,'box optimized fallback')
      call les_filter_weights('box',12,fallback,ok)
      call require(all(w==fallback),'box fallback equals standard kernel')
      call les_filter_weights('gaussian3opt',11,w,ok)
      call require(ok,'Gaussian three point fallback')
      call les_filter_weights('gaussian3',11,fallback,ok)
      call require(all(w==fallback),'Gaussian three point fallback equality')
      call les_filter_weights('gaussian5opt',11,w,ok)
      call require(ok,'Gaussian five point fallback')
      call les_filter_weights('gaussian5',11,fallback,ok)
      call require(all(w==fallback),'Gaussian five point fallback equality')
      call les_filter_weights('box3opt',11,w,ok)
      call require(.not.ok,'odd standard box fallback rejected')
    end block
  case('les_filter')
    block
      real(dp), allocatable :: w(:)
      real(dp) :: field(2,4,3,2),filtered(2,4,3,2),halo(2,10,9,8),cropped(2,4,3,2),before(2,4,3,2)
      integer :: i,j,k
      call les_filter_weights('box',2,w,ok)
      call require(ok.and.size(w)==3,'box width')
      call require(all(w==[.25_dp,.5_dp,.25_dp]),'box endpoint weights')
      field=0;field(1,1,1,1)=1;field(2,:,:,:)=7
      call filter_les_field_3d(field,w,filtered,ok)
      call require(ok,'periodic impulse filter')
      call near(filtered(1,1,1,1),.125_dp,1e-13_dp,'tensor impulse center')
      call near(sum(filtered(1,:,:,:)),1._dp,1e-13_dp,'filter conservation')
      call require(all(filtered(2,:,:,:)==7),'filter constant')
      do i=1,4
        field(1,i,:,:)=cos(2*pi*real(i-1,dp)/4)
      end do
      call filter_les_field_3d(field,w,filtered,ok)
      call require(ok.and.maxval(abs(filtered(1,:,:,:)-.5_dp*field(1,:,:,:)))<1e-14_dp,'Fourier attenuation')
      call les_filter_weights('box',6,w,ok)
      call require(ok.and.size(w)==7,'coefficient smoothing width')
      call near(w(1),1._dp/12,1e-13_dp,'width six endpoint')
      call near(w(4),1._dp/6,1e-13_dp,'width six center')
      do k=1,8
        do j=1,9
          do i=1,10
            halo(:,i,j,k)=field(:,modulo(i-4,4)+1,modulo(j-4,3)+1,modulo(k-4,2)+1)
          end do
        end do
      end do
      call filter_les_field_3d(field,w,filtered,ok)
      call require(ok,'wide periodic filter on small grid')
      call filter_les_field_3d(halo,w,cropped,ok,[4,4,4])
      call require(ok.and.all(cropped==filtered),'halo and periodic filters agree')
      before=cropped
      call filter_les_field_3d(halo,w,cropped,ok,[3,4,4])
      call require(.not.ok.and.all(cropped==before),'insufficient halo rollback')
      halo(1,1,1,1)=ieee_value(1._dp,ieee_quiet_nan)
      call filter_les_field_3d(halo,w,cropped,ok,[4,4,4])
      call require(.not.ok.and.all(cropped==before),'nonfinite halo rollback')
      call les_filter_weights('gaussian',2,w,ok)
      call require(ok,'Gaussian weights')
      call near(w(1)/w(2),exp(-1.5_dp),1e-13_dp,'Gaussian analytic weight ratio')
      call near(sum(w),1._dp,1e-13_dp,'Gaussian normalization')
      call filter_les_field_3d(field,w,filtered,ok)
      call require(ok,'Gaussian periodic filter')
      call near(filtered(1,1,1,1),w(2),1e-13_dp,'Gaussian Fourier response')
      call les_filter_weights('none',2,w,ok)
      call require(ok,'identity weights')
      call filter_les_field_3d(field,w,filtered,ok)
      call require(ok.and.all(filtered==field),'identity filter')
      before=filtered;w=huge(1._dp)
      call filter_les_field_3d(field,w,filtered,ok)
      call require(.not.ok.and.all(filtered==before),'oversized weights rollback')
      call les_filter_weights('box',3,w,ok)
      call require(.not.ok.and..not.allocated(w),'odd box ratio rejected')
      call les_filter_weights('gaussian',1,w,ok)
      call require(.not.ok.and..not.allocated(w),'odd Gaussian ratio rejected')
    end block
  case('les_dynamic_algebra')
    block
      real(dp) :: kij(3,3),rut(3),aij(3,3),a,ft(3),bij(3,3),b,test_ft(3),coeff(3),vel(3),leonard(3,3)
      real(dp), parameter :: eps=1e-8_dp
      integer :: i,j
      grad=0;grad(1,2)=3;vel=[1._dp,2._dp,3._dp]
      call dynamic_smagorinsky_quantities(grad,[2._dp,-1._dp,0._dp],2._dp,vel,300._dp, &
        [.1_dp,.2_dp,.4_dp],2,kij,rut,aij,a,ft,ok)
      call require(ok,'dynamic grid quantities')
      ! Delta=0.4, |S|=3, mu/rho=delta^2 |S|.
      call near(aij(1,2),2*.4_dp**2*3*3,1e-13_dp,'dynamic shear tensor')
      call near(a,2*2*.4_dp**2*9,1e-13_dp,'dynamic isotropic quantity')
      call near(ft(1),2*2*.4_dp**2*3,1e-13_dp,'dynamic thermal quantity')
      call near(kij(2,3),12._dp,1e-13_dp,'dynamic density weighted velocity product')
      call near(rut(3),1800._dp,1e-13_dp,'dynamic density weighted temperature product')
      call require(all(aij==transpose(aij)),'dynamic tensor symmetry')
      ! A manufactured Germano tensor includes shear, trace and signed heat flux.
      aij=0;a=0;bij=0;bij(1,1)=2;bij(2,2)=-1;bij(3,3)=-1;bij(1,2)=1;bij(2,1)=1
      leonard=2*bij
      do i=1,3
        leonard(i,i)=leonard(i,i)+3
      end do
      do j=1,3
        do i=1,3
          kij(i,j)=2*vel(i)*vel(j)+leonard(i,j)
        end do
      end do
      b=6;test_ft=0;ft=[1._dp,-2._dp,2._dp];rut=600*vel+4*ft
      call dynamic_smagorinsky_coefficients(2._dp,vel,300._dp,kij,rut,aij,a,ft,bij,b,test_ft,coeff,ok)
      call require(ok,'dynamic coefficient solve')
      call near(coeff(1),16/(8+eps),1e-13_dp,'dynamic double off-diagonal contraction')
      call near(coeff(2),1._dp,1e-13_dp,'dynamic CI upper bound')
      call near(coeff(3),36/(9+eps),1e-13_dp,'dynamic thermal contraction')
      b=18
      call dynamic_smagorinsky_coefficients(2._dp,vel,300._dp,kij,rut,aij,a,ft,bij,b,test_ft,coeff,ok)
      call require(ok,'dynamic interior CI solve')
      call near(coeff(2),9/(18+eps),1e-13_dp,'dynamic unclipped CI')
      kij=0;rut=0;bij=0;b=0;ft=0;vel=0
      call dynamic_smagorinsky_coefficients(2._dp,vel,300._dp,kij,rut,aij,a,ft,bij,b,test_ft,coeff,ok)
      call require(ok.and.all(coeff==eps),'dynamic zero tensors retain pinned floor')
      kij=-1;rut=-1;bij=1;b=1;ft=1
      call dynamic_smagorinsky_coefficients(2._dp,vel,300._dp,kij,rut,aij,a,ft,bij,b,test_ft,coeff,ok)
      call require(ok.and.all(coeff==eps),'dynamic negative coefficients clipped')
      b=-eps
      call dynamic_smagorinsky_coefficients(2._dp,vel,300._dp,kij,rut,aij,a,ft,bij,b,test_ft,coeff,ok)
      call require(.not.ok.and.all(coeff==0),'dynamic singular CI rejects atomically')
      b=ieee_value(1._dp,ieee_quiet_nan)
      call dynamic_smagorinsky_coefficients(2._dp,vel,300._dp,kij,rut,aij,a,ft,bij,b,test_ft,coeff,ok)
      call require(.not.ok.and.all(coeff==0),'dynamic nonfinite coefficient input')
      call dynamic_smagorinsky_quantities(grad,[2._dp,-1._dp,0._dp],2._dp,vel,300._dp, &
        [.1_dp,.2_dp,.4_dp],0,kij,rut,aij,a,ft,ok)
      call require(.not.ok.and.all(kij==0).and.all(rut==0).and.all(aij==0).and.a==0.and.all(ft==0), &
        'dynamic invalid filter width')
    end block
  case('les_alternative_algebra')
    grad=0
    call wale_stress(grad,2._dp,[.1_dp,.1_dp,.1_dp],.5_dp,mu,stress,ok)
    call require(ok.and.mu==0.and.all(stress==0),'WALE zero gradient')
    call vreman_stress(grad,2._dp,[.1_dp,.1_dp,.1_dp],.2_dp,mu,stress,ok)
    call require(ok.and.mu==0.and.all(stress==0),'Vreman zero gradient')
    grad(1,2)=3
    call wale_stress(grad,2._dp,[.1_dp,.1_dp,.1_dp],.5_dp,mu,stress,ok)
    call require(ok.and.mu==0.and.all(stress==0),'WALE simple shear')
    call vreman_stress(grad,2._dp,[.1_dp,.1_dp,.1_dp],.2_dp,mu,stress,ok)
    call require(ok.and.mu==0.and.all(stress==0),'Vreman simple shear')
    grad=0;grad(1,1)=1;grad(2,2)=2;grad(3,3)=3
    call wale_stress(grad,2._dp,[.1_dp,.2_dp,.4_dp],.5_dp,mu,stress,ok)
    call require(ok,'WALE diagonal strain')
    expected=2*(.5_dp*.2_dp)**2*(98._dp/3)**1.5_dp/(14._dp**2.5_dp+(98._dp/3)**1.25_dp)
    call near(mu,expected,1e-13_dp,'WALE analytic invariant')
    call near(stress(1,1),-2*mu,1e-13_dp,'WALE deviatoric stress')
    call near(stress(1,1)+stress(2,2),-stress(3,3),1e-13_dp,'WALE trace')
    call vreman_stress(grad,2._dp,[.1_dp,.2_dp,.4_dp],.2_dp,mu,stress,ok)
    call require(ok,'Vreman diagonal strain')
    call near(mu,2*(.2_dp*.2_dp)**2*2.5_dp*sqrt(3.5_dp),1e-13_dp,'Vreman analytic invariant')
    call near(stress(3,3),2*mu,1e-13_dp,'Vreman deviatoric stress')
    grad=grad*1e29_dp
    call wale_stress(grad,2._dp,[.1_dp,.2_dp,.4_dp],.5_dp,mu,stress,ok)
    call require(ok.and.all(ieee_is_finite(stress)),'WALE scaled large gradient')
    call near(mu/1e29_dp,expected,1e-13_dp,'WALE scaling')
    call vreman_stress(grad,-1._dp,[.1_dp,.2_dp,.4_dp],.2_dp,mu,stress,ok)
    call require(.not.ok.and.mu==0.and.all(stress==0),'Vreman invalid density')
  case('les_wale','les_vreman')
    les_model=test(5:)
    do les_case=1,2
      call make_gas(.true.)
      do jcell=1,4
        state(2,:,jcell,:)=10*sin(2*pi*(real(jcell,dp)-.5_dp)/4)*state(1,:,jcell,:)
      end do
      state(5,:,:,:)=state(5,:,:,:)+.5_dp*state(2,:,:,:)**2/state(1,:,:,:)
      saved=state;saved_t=temp
      faces='periodic'
      if (les_case==2) faces='outflow'
      call compute_les_rhs_3d(species,state,temp,spacing,.17_dp,.7_dp,.7_dp,rhs,diff,ok,faces,totals0,model=les_model)
      call require(ok.and.diff>0,'alternative LES flux')
      do n=1,15
        call require(abs(sum(rhs(n,:,:,:))*product(spacing)-totals0(n)) &
          <=1e-12_dp*max(1._dp,sum(abs(rhs(n,:,:,:))))*product(spacing),'LES flux balance')
      end do
      ke0=sum(.5_dp*sum(state(2:4,:,:,:)**2,dim=1)/state(1,:,:,:))
      call advance_les_ssprk2_3d(species,state,temp,spacing,.17_dp,.7_dp,.7_dp,1e-7_dp,ok,faces,totals1,model=les_model)
      call require(ok,'alternative LES step')
      do n=1,15
        call require(abs((sum(state(n,:,:,:))-sum(saved(n,:,:,:)))*product(spacing)-totals1(n)) &
          <=1e-12_dp*max(1._dp,sum(abs(saved(n,:,:,:))))*product(spacing), 'LES step boundary budget')
      end do
      if (les_case==1) then
        ke1=sum(.5_dp*sum(state(2:4,:,:,:)**2,dim=1)/state(1,:,:,:))
        call require(ke1<ke0,'alternative LES kinetic dissipation')
      end if
      saved=state;saved_t=temp
      call advance_les_ssprk2_3d(species,state,temp,spacing,.17_dp,.7_dp,.7_dp,1e6_dp,ok,faces,model=les_model)
      call require(.not.ok.and.all(state==saved).and.all(temp==saved_t),'alternative LES timestep rollback')
      call advance_les_ssprk2_3d(species,state,temp,spacing,0._dp,.7_dp,.7_dp,0._dp,ok,model='invalid')
      call require(.not.ok.and.all(state==saved).and.all(temp==saved_t),'invalid LES model rollback')
    end do
  case('les_boundary')
    call make_gas(.false.)
    state(2,:,:,:)=20.0_dp*state(1,:,:,:)
    state(5,:,:,:)=state(5,:,:,:)+200.0_dp*state(1,:,:,:)
    saved=state; saved_t=temp
    faces='periodic'; faces(3:4)='no_slip_wall'
    call compute_les_rhs_3d(species,state,temp,spacing,.17_dp,.7_dp,.7_dp,rhs,diff,ok,faces,totals0)
    call require(ok.and.diff>0,'physical SGS flux')
    mu=state(1,1,1,1)*(.17_dp*spacing(1))**2*(40.0_dp/spacing(2))
    expected=-4*mu*20.0_dp/spacing(2)*(4*spacing(1))*(4*spacing(3))
    call near(totals0(2),expected,1e-12_dp,'analytical SGS wall shear')
    call require(abs(totals0(5))<1e-18_dp,'stationary SGS wall no work')
    do n=1,15
      call require(abs(sum(rhs(n,:,:,:))*product(spacing)-totals0(n))< &
        1e-12_dp*max(1e-12_dp,abs(totals0(n))),'SGS boundary divergence budget')
    end do
    call advance_les_ssprk2_3d(species,state,temp,spacing,.17_dp,.7_dp,.7_dp,1e-7_dp,ok,faces,totals1)
    call require(ok,'physical SGS RK update')
    call near(sum(state(2,:,:,:)-saved(2,:,:,:))*product(spacing),totals1(2),1e-10_dp,'RK wall ledger')
    state=saved;temp=saved_t
    call advance_les_ssprk2_3d(species,state,temp,spacing,.17_dp,.7_dp,.7_dp,1e6_dp,ok,faces,totals1)
    call require(.not.ok.and.all(state==saved).and.all(temp==saved_t).and.all(totals1==0),'physical SGS rollback')
    wall_v=0;wall_v(1,3:4)=10;wall_t=temp(1,1,1)
    faces(3:4)='isothermal_wall'
    call compute_les_rhs_3d(species,state,temp,spacing,.17_dp,.7_dp,.7_dp,rhs,diff,ok,faces,totals0, &
      wall_temperature=wall_t,wall_velocity=wall_v)
    call require(ok,'moving isothermal SGS wall')
    call near(totals0(5),10*totals0(2),1e-12_dp,'SGS wall mechanical work')
    wall_t=wall_t+100
    call compute_les_rhs_3d(species,state,temp,spacing,.17_dp,.7_dp,.7_dp,rhs,diff,ok,faces,totals1, &
      wall_temperature=wall_t,wall_velocity=wall_v)
    call require(ok.and.totals1(5)>totals0(5),'hot wall adds SGS heat')
    call near(sum(rhs(5,:,:,:))*product(spacing),totals1(5),1e-12_dp,'thermal SGS wall ledger')
    wall_v(2,3)=1
    call advance_les_ssprk2_3d(species,state,temp,spacing,.17_dp,.7_dp,.7_dp,1e-7_dp,ok,faces,totals1, &
      wall_temperature=wall_t,wall_velocity=wall_v)
    call require(.not.ok.and.all(state==saved).and.all(temp==saved_t),'normal wall motion rejects transactionally')
  case('les_grid','les_zero','les_reject')
    call make_gas(.true.); saved=state; saved_t=temp
    ke0=sum(0.5_dp*sum(state(2:4,:,:,:)**2,dim=1)/state(1,:,:,:))
    call compute_les_rhs_3d(species,state,temp,spacing,0.17_dp,0.7_dp,0.7_dp,rhs,diff,ok)
    call require(ok .and. diff > 0.0_dp,'nonzero SGS in shear')
    do n=1,15
      call require(abs(sum(rhs(n,:,:,:))) <= 1.0e-12_dp*max(1.0_dp,sum(abs(rhs(n,:,:,:)))), &
        'periodic flux conservation')
    end do
    dt=1.0e-6_dp; mu=0.17_dp
    if (test == 'les_zero') mu=0.0_dp
    if (test == 'les_reject') dt=1.0e6_dp
    call advance_les_ssprk2_3d(species,state,temp,spacing,mu,0.7_dp,0.7_dp,dt,ok)
    if (test == 'les_reject') then
      call require(.not. ok .and. all(state == saved) .and. all(temp == saved_t),'SGS timestep rollback')
    else
      call require(ok,'SGS update')
      call near(sum(state(5,:,:,:)),sum(saved(5,:,:,:)),1.0e-13_dp,'SGS total energy')
      if (test == 'les_zero') then
        call require(all(state == saved) .and. all(temp == saved_t),'Cs zero exact identity')
      else
        ke1=sum(0.5_dp*sum(state(2:4,:,:,:)**2,dim=1)/state(1,:,:,:))
        call require(ke1 < ke0,'SGS dissipates resolved kinetic energy')
      end if
    end if
  case default
    error stop 'Unknown spray test'
  end select
  print '(a)',trim(test)//': PASS'
contains
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(len=*),intent(in)::label
    if (.not. condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
  subroutine near(actual,reference,tolerance,label)
    real(dp),intent(in)::actual,reference,tolerance
    character(len=*),intent(in)::label
    if (abs(actual-reference) > tolerance*max(1.0e-25_dp,abs(reference))) then
      print *,label,actual,reference
      error stop 1
    end if
  end subroutine
  subroutine exchange_check(a,b,exchange)
    type(spray_parcel),intent(in)::a,b
    real(dp),intent(in)::exchange(5)
    real(dp)::ta(5),tb(5)
    ta=parcel_totals(a,liquid); tb=parcel_totals(b,liquid)
    call require(maxval(abs(ta-tb-exchange)) <= 1.0e-14_dp*maxval(abs(ta)),'discrete exchange identity')
  end subroutine
  subroutine make_drops()
    call inject_cone_parcels(liquid,1.0e-10_dp,4,1_int64,[0.0015_dp,0.0015_dp,0.0015_dp], &
      [0.0_dp,0.0_dp,1.0_dp],10.0_dp,0.1_dp,2.0e-5_dp,330.0_dp,drops,ok)
    call require(ok,'cone injection')
  end subroutine
  subroutine make_gas(shear)
    logical,intent(in)::shear
    real(dp)::prim(15),cons(15),y(10),gas_r,t,sound
    integer::i,j,k
    y=0.0_dp; y(10)=1.0_dp
    gas_r=mixture_specific_gas_constant(species,y,ok)
    call require(ok,'mixture R')
    prim=0.0_dp;prim(1)=101325.0_dp/(gas_r*1000.0_dp);prim(5)=101325.0_dp;prim(6:)=y
    do k=1,4
      do j=1,4
        do i=1,4
          if (shear) prim(3)=10.0_dp*sin(2.0_dp*pi*(real(i,dp)-0.5_dp)/4.0_dp)
          call reactive_primitive_to_conserved(species,prim,cons,t,sound,ok)
          call require(ok,'initialize gas')
          state(:,i,j,k)=cons;temp(i,j,k)=t
        end do
      end do
    end do
  end subroutine
end program
