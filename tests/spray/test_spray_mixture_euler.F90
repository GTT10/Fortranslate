program test_spray_mixture_euler
  use precision_mod, only: dp
  use spray_parcel_mod, only: spray_liquid,spray_parcel,spray_options,valid_spray_options,liquid_energy
  use spray_mixture_mod, only: advance_mixture_evaporation
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
  options%mixture_time_scheme='euler'
  if(.not.valid_spray_options(options)) error stop 'Euler option rejected'
  options%mixture_time_scheme='unknown'
  if(valid_spray_options(options)) error stop 'Invalid time scheme accepted'
  print '(a)','spray_mixture_euler: PASS'
end program
