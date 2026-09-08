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
  use spray_coupling_3d_mod
  use les_smagorinsky_3d_mod
  implicit none
  real(dp) :: value
  real(dp), parameter :: pi=acos(-1.0_dp)
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(spray_liquid) :: liquid
  type(spray_parcel) :: p,q
  type(spray_parcel), allocatable :: drops(:),saved_drops(:)
  type(spray_options) :: options
  type(spray_film) :: film
  real(dp) :: state(15,4,4,4),temp(4,4,4),saved(15,4,4,4),saved_t(4,4,4),rhs(15,4,4,4)
  real(dp) :: totals0(15),totals1(15),spacing(3),lower(3),weights(8)
  real(dp) :: primitive(15),binary(10,10)
  real(dp) :: inverse_tau,heat,masscoef,transfer(5),expected,dt,d,grad(3,3),stress(3,3),mu,diff,ke0,ke1
  integer :: indices(3,8),n,steps
  logical :: ok
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
