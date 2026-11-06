program test_spray_mixture_euler
  use precision_mod, only: dp
  use spray_parcel_mod, only: spray_liquid,spray_parcel,spray_options,valid_spray_options,liquid_energy
  use spray_mixture_mod, only: advance_mixture_evaporation,advance_mixture_frozen_gas,mixture_transfer_rates
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_full_transport
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  implicit none
  type(spray_liquid) :: liquid(1)
  type(spray_parcel) :: p,saved
  type(spray_options) :: options
  real(dp) :: mass(1),dm(1),hv(1),exchange(6),before,after,saved_mass(1)
  logical :: ok
  liquid(1)%vapor_name='H2O';liquid(1)%vapor_index=1
  liquid(1)%density=1000;liquid(1)%heat_capacity=2000
  liquid(1)%reference_temperature=300;liquid(1)%reference_pressure=101325
  liquid(1)%latent_heat=400000
  p%id=1;p%mass=1e-9_dp;p%temperature=350;p%multiplicity=1;p%velocity=[1._dp,0._dp,0._dp]
  mass=p%mass;dm=1e-11_dp;hv=500000
  before=p%mass*(liquid_energy(liquid(1),p%temperature)+.5_dp*sum(p%velocity**2))
  call advance_mixture_evaporation(p,liquid,mass,dm,hv,1e-5_dp,exchange,ok,euler_temperature=.true.)
  if(.not.ok) error stop 'Euler mixture update rejected'
  ! Old mass*cp is 2e-6 J/K; net heating is 1e-5 - 1e-11*400000 = 6e-6 J.
  if(abs(p%temperature-353._dp)>1e-12_dp) error stop 'Incorrect Euler temperature'
  if(abs(p%mass-.99e-9_dp)>1e-24_dp) error stop 'Incorrect remaining mass'
  if(abs(exchange(1)-1e-11_dp)>1e-25_dp) error stop 'Incorrect gas mass'
  after=p%mass*(liquid_energy(liquid(1),p%temperature)+.5_dp*sum(p%velocity**2))
  if(abs(before-after-exchange(5))>1e-13_dp*abs(before)) error stop 'Energy not conserved'
  saved=p;saved_mass=mass
  call advance_mixture_evaporation(p,liquid,mass,dm,hv,-1e-3_dp,exchange,ok,euler_temperature=.true.)
  if(ok.or.p%mass/=saved%mass.or.p%temperature/=saved%temperature.or.any(mass/=saved_mass)) &
    error stop 'Failed Euler update changed liquid'
  if(any(exchange/=0)) error stop 'Failed Euler update published gas exchange'
  block
    type(spray_liquid) :: components(3)
    real(dp) :: masses(3),losses(3),enthalpies(3),old_masses(3)
    components=liquid(1);p=saved;p%mass=1e-9_dp;p%temperature=350
    old_masses=p%mass*[.1_dp,.4_dp,.5_dp];masses=old_masses
    losses=p%mass*[.6_dp,0._dp,0._dp]
    enthalpies=liquid_energy(liquid(1),p%temperature)
    before=p%mass*(enthalpies(1)+.5_dp*sum(p%velocity**2))
    call advance_mixture_evaporation(p,components,masses,losses,enthalpies,0._dp,exchange,ok,euler_temperature=.true.)
    if(.not.ok) error stop 'Euler component overdraw rejected'
    if(maxval(abs(masses/1e-9_dp-[0._dp,.2_dp,.2_dp]))>2e-14_dp) &
      error stop 'Incorrect clipped Euler component inventory'
    if(abs(exchange(1)+p%mass-1e-9_dp)>1e-24_dp) error stop 'Clipped Euler total mass not conserved'
    after=p%mass*(liquid_energy(liquid(1),p%temperature)+.5_dp*sum(p%velocity**2))
    if(abs(before-after-exchange(5))>1e-13_dp*abs(before)) error stop 'Clipped Euler energy not conserved'
    p=saved;p%mass=1e-9_dp;p%temperature=350;masses=old_masses
    call advance_mixture_evaporation(p,components,masses,losses,enthalpies,0._dp,exchange,ok)
    if(ok.or.any(masses/=old_masses).or.p%mass/=1e-9_dp) error stop 'Conservative overdraw rollback changed'
    masses=old_masses;losses=2*old_masses
    call advance_mixture_evaporation(p,components,masses,losses,enthalpies,0._dp,exchange,ok,euler_temperature=.true.)
    if(.not.ok.or.p%mass/=0.or.any(masses/=0)) error stop 'Euler total depletion failed'
  end block
  options%mixture_time_scheme='euler'
  if(.not.valid_spray_options(options)) error stop 'Euler option rejected'
  block
    type(nasa7_species), allocatable :: species(:)
    type(gas_transport_species), allocatable :: transport(:)
    type(spray_liquid) :: material(1)
    real(dp) :: primitive(15),feedback(15),force(5),accel(3),dt,invtau,heat,rate(1),enthalpy(1)
    call load_h2o2_full_thermo(species,ok)
    if(.not.ok) error stop 'Euler drift thermo'
    call load_h2o2_full_transport(transport,ok)
    if(.not.ok) error stop 'Euler drift transport'
    material=liquid;material(1)%vapor_index=6
    primitive=0;primitive(1)=1;primitive(2:4)=[3._dp,1._dp,-2._dp];primitive(15)=1
    primitive(5)=800*mixture_specific_gas_constant(species,primitive(6:),ok)
    if(.not.ok) error stop 'Euler drift gas constant'
    p=saved;p%position=[.1_dp,.2_dp,.3_dp];mass=p%mass;before=p%position(1)
    saved=p;options%drag=.false.;options%heat=.false.;options%evaporation=.false.
    dt=1e-3_dp;accel=[4._dp,2._dp,-1._dp]
    call advance_mixture_frozen_gas(species,transport,primitive,800._dp,p,material,mass,options,dt,feedback,ok, &
      time_modifier=.5_dp,acceleration=accel,force_change=force)
    if(.not.ok) error stop 'Euler drift with modified kick'
    if(maxval(abs(p%velocity-saved%velocity-.5_dp*dt*accel))>1e-14_dp) error stop 'Euler modified kick velocity'
    if(maxval(abs(p%position-saved%position-dt*p%velocity))>1e-16_dp) error stop 'Euler outer flow drift'
    p=saved;mass=p%mass;options%drag=.true.
    call mixture_transfer_rates(species,transport,primitive,800._dp,p,material,[1._dp],options, &
      invtau,heat,rate,enthalpy,ok)
    if(.not.ok.or.invtau<=0) error stop 'Euler multistep drag rate'
    dt=1.25_dp/invtau
    call advance_mixture_frozen_gas(species,transport,primitive,800._dp,p,material,mass,options,dt,feedback,ok)
    if(.not.ok) error stop 'Euler internally substepped drift'
    if(maxval(abs(p%position-saved%position-dt*p%velocity))>1e-15_dp) &
      error stop 'Internal drag substeps must not drift the particle'
  end block
  options%mixture_time_scheme='unknown'
  if(valid_spray_options(options)) error stop 'Invalid time scheme accepted'
  print '(a)','spray_mixture_euler: PASS'
end program
