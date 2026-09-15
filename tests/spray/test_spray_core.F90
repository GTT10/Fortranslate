program test_spray_core
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_full_transport
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use mixture_transport_mod, only: binary_diffusion_coefficients
  use spray_parcel_mod
  use spray_breakup_mod
  use spray_coupling_3d_mod
  use les_smagorinsky_3d_mod
  implicit none
  real(dp) :: value,td,omega,wer,initial_mass,z,tkh,wavelength,tau,acceleration
  real(dp), parameter :: pi=acos(-1.0_dp)
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(spray_liquid) :: liquid,constant_liquid
  type(spray_parcel) :: p,q
  type(spray_parcel), allocatable :: drops(:),saved_drops(:)
  type(spray_options) :: options,constant_options
  type(spray_film) :: film
  real(dp) :: state(15,4,4,4),temp(4,4,4),saved(15,4,4,4),saved_t(4,4,4),rhs(15,4,4,4)
  real(dp) :: wall_t(6),wall_v(3,6)
  real(dp) :: totals0(15),totals1(15),spacing(3),lower(3),weights(8)
  real(dp) :: primitive(15),binary(10,10)
  real(dp) :: inverse_tau,heat,masscoef,transfer(5),expected,dt,d,grad(3,3),stress(3,3),mu,diff,ke0,ke1
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
