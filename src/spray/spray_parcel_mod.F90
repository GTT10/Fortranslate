module spray_parcel_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species, nasa7_mass_properties, universal_gas_constant
  implicit none
  private
  real(dp), parameter :: pi = acos(-1.0_dp)

  ! SI throughout. The stored liquid energy uses the gas mechanism's formation
  ! reference, neglecting the incompressible p/rho correction and surface energy.
  type, public :: spray_liquid
    character(len=24) :: vapor_name = ''
    integer :: vapor_index = 0
    real(dp) :: density = 0.0_dp, heat_capacity = 0.0_dp
    real(dp) :: reference_temperature = 0.0_dp, reference_pressure = 0.0_dp
    real(dp) :: latent_heat = 0.0_dp, reference_energy = 0.0_dp
    character(len=16) :: vapor_pressure_model='clausius'
    real(dp) :: antoine_coefficients(4)=0.0_dp ! A, B [K], C [K], D [Pa]
    real(dp) :: antoine_temperature_range(2)=0.0_dp
    character(len=16) :: density_model='constant'
    real(dp) :: density_coefficients(4)=0.0_dp, density_temperature_range(2)=0.0_dp
    character(len=16) :: boiling_model='reject'
    real(dp) :: critical_temperature=0.0_dp
  end type
  type, public :: spray_parcel
    integer(int64) :: id = 0_int64
    real(dp) :: position(3) = 0.0_dp, velocity(3) = 0.0_dp
    real(dp) :: mass = 0.0_dp, temperature = 0.0_dp, multiplicity = 0.0_dp
    real(dp) :: deformation=0.0_dp,deformation_rate=0.0_dp
    integer(int64) :: parent_id=0_int64
    real(dp) :: rt_time=-1.0_dp
  end type
  type, public :: spray_options
    real(dp) :: khrt_b0=.61_dp,khrt_b1=7.0_dp,khrt_c3=1.0_dp
    real(dp) :: maximum_multiplicity=0.0_dp
    integer :: maximum_parcels=1000000
    character(len=8) :: breakup_model='none'
    real(dp) :: surface_tension=0.0_dp,liquid_viscosity=0.0_dp
    logical :: drag = .true., heat = .true., evaporation = .true.
    logical :: abramzon_sirignano = .true.
    real(dp) :: max_fractional_change = 0.1_dp, particle_cfl = 0.25_dp
    real(dp) :: minimum_diameter = 1.0e-9_dp
    integer :: maximum_substeps = 100000
    character(len=24) :: liquid_viscosity_model='constant'
    real(dp) :: liquid_viscosity_coefficients(4)=0.0_dp,liquid_viscosity_temperature_range(2)=0.0_dp
  end type
  type, public :: spray_film
    real(dp) :: density = 0.0_dp, viscosity = 0.0_dp
    real(dp) :: conductivity = 0.0_dp, diffusivity = 0.0_dp, cp = 0.0_dp
    real(dp) :: gas_temperature = 0.0_dp, gas_velocity(3) = 0.0_dp
    real(dp) :: surface_fraction = 0.0_dp, ambient_fraction = 0.0_dp
    real(dp) :: vapor_enthalpy = 0.0_dp
  end type
  public :: initialize_spray_liquid, valid_spray_liquid, valid_spray_parcel
  public :: valid_spray_options, parcel_diameter, liquid_energy, parcel_totals
  public :: spray_transfer_rates, advance_parcel_frozen_film, saturation_mole_fraction
  public :: inject_cone_ports
  public :: inject_cone_parcels, liquid_density_at_temperature, liquid_viscosity_at_temperature
contains
  subroutine initialize_spray_liquid(species, name, density, cp, tref, pref, latent, liquid, ok, &
      vapor_pressure_model, antoine_coefficients, antoine_temperature_range)
    type(nasa7_species), intent(in) :: species(:)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: density, cp, tref, pref, latent
    character(len=*), intent(in), optional :: vapor_pressure_model
    real(dp), intent(in), optional :: antoine_coefficients(4),antoine_temperature_range(2)
    type(spray_liquid), intent(out) :: liquid
    logical, intent(out) :: ok
    real(dp) :: vgcp, cv, h, e, s
    integer :: k, index
    liquid = spray_liquid()
    ok = .false.
    if (.not. all(ieee_is_finite([density, cp, tref, pref, latent]))) return
    if (density < 1.0e-3_dp .or. density > 1.0e5_dp) return
    if (cp < 1.0_dp .or. cp > 1.0e8_dp) return
    if (tref < 1.0_dp .or. tref > 1.0e5_dp) return
    if (pref < 1.0e-6_dp .or. pref > 1.0e12_dp) return
    if (latent <= 0.0_dp .or. latent > 1.0e10_dp) return
    index = 0
    do k = 1, size(species)
      if (trim(species(k)%name) /= trim(name)) cycle
      if (index /= 0) return
      index = k
    end do
    if (index == 0) return
    call nasa7_mass_properties(species(index), tref, vgcp, cv, h, e, s, ok)
    if (.not. ok) return
    liquid = spray_liquid(trim(name), index, density, cp, tref, pref, latent, h-latent)
    if (present(vapor_pressure_model)) liquid%vapor_pressure_model=vapor_pressure_model
    if (present(antoine_coefficients)) liquid%antoine_coefficients=antoine_coefficients
    if (present(antoine_temperature_range)) liquid%antoine_temperature_range=antoine_temperature_range
    ok = valid_spray_liquid(liquid)
  end subroutine

  pure logical function valid_spray_liquid(l) result(ok)
    type(spray_liquid), intent(in) :: l
    ok = .false.
    if (.not. all(ieee_is_finite([l%density, l%heat_capacity, l%reference_temperature, &
        l%reference_pressure, l%latent_heat, l%reference_energy]))) return
    if (.not.all(ieee_is_finite(l%antoine_coefficients))) return
    if (.not.all(ieee_is_finite(l%antoine_temperature_range))) return
    select case(trim(l%vapor_pressure_model))
    case('clausius')
      if (any(l%antoine_coefficients/=0).or.any(l%antoine_temperature_range/=0)) return
    case('antoine')
      if (maxval(abs(l%antoine_coefficients(1:3)))>1e10_dp) return
      if (l%antoine_coefficients(2)<=0) return
      if (l%antoine_coefficients(4)<1e-6_dp.or.l%antoine_coefficients(4)>1e12_dp) return
      if (l%antoine_temperature_range(1)<1.or.l%antoine_temperature_range(2)>1e5_dp) return
      if (l%antoine_temperature_range(2)<=l%antoine_temperature_range(1)) return
      if (l%antoine_temperature_range(1)+l%antoine_coefficients(3)<=0) return
    case default
      return
    end select
    if (.not.ieee_is_finite(l%critical_temperature)) return
    select case(trim(l%boiling_model))
    case('reject')
      if (l%critical_temperature/=0) return
    case('surface-cap')
      if (trim(l%vapor_pressure_model)/='clausius') return
      if (l%critical_temperature<=l%reference_temperature.or.l%critical_temperature>1e5_dp) return
    case default
      return
    end select
    if (.not.all(ieee_is_finite(l%density_coefficients))) return
    if (.not.all(ieee_is_finite(l%density_temperature_range))) return
    select case(trim(l%density_model))
    case('constant')
      if (any(l%density_coefficients/=0).or.any(l%density_temperature_range/=0)) return
    case('polynomial')
      if (maxval(abs(l%density_coefficients))>1e10_dp) return
      if (l%density_temperature_range(1)<1.or.l%density_temperature_range(2)>1e5_dp) return
      if (l%density_temperature_range(2)<=l%density_temperature_range(1)) return
    case default
      return
    end select
    ok = len_trim(l%vapor_name) > 0 .and. l%vapor_index > 0 .and. &
      l%density >= 1.0e-3_dp .and. l%density <= 1.0e5_dp .and. &
      l%heat_capacity >= 1.0_dp .and. l%heat_capacity <= 1.0e8_dp .and. &
      l%reference_temperature >= 1.0_dp .and. l%reference_temperature <= 1.0e5_dp .and. &
      l%reference_pressure >= 1.0e-6_dp .and. l%reference_pressure <= 1.0e12_dp .and. &
      l%latent_heat > 0.0_dp .and. l%latent_heat <= 1.0e10_dp .and. &
      abs(l%reference_energy) <= 1.0e10_dp
  end function

  pure logical function valid_spray_parcel(p) result(ok)
    type(spray_parcel), intent(in) :: p
    ok = .false.
    if (.not. all(ieee_is_finite([p%position, p%velocity, p%mass, p%temperature, &
        p%multiplicity,p%deformation,p%deformation_rate,p%rt_time]))) return
    if (p%rt_time < -1.0_dp.or.p%rt_time>1e6_dp) return
    if (p%parent_id<0_int64.or.p%parent_id>=p%id) return
    if (abs(p%deformation)>1e6_dp.or.abs(p%deformation_rate)>1e20_dp) return
    ok = p%id > 0_int64 .and. p%mass >= 0.0_dp .and. p%mass <= 1.0_dp .and. &
      p%temperature >= 1.0_dp .and. p%temperature <= 1.0e5_dp .and. &
      p%multiplicity > 0.0_dp .and. p%multiplicity <= 1.0e18_dp .and. &
      maxval(abs(p%position)) <= 1.0e9_dp .and. maxval(abs(p%velocity)) <= 1.0e7_dp
  end function

  pure logical function valid_spray_options(o) result(ok)
    type(spray_options), intent(in) :: o
    ok = .false.
    if (.not. all(ieee_is_finite([o%max_fractional_change,o%particle_cfl,o%minimum_diameter]))) return
    if (.not.ieee_is_finite(o%maximum_multiplicity)) return
    if (o%maximum_multiplicity<0.or.o%maximum_multiplicity>1e18_dp) return
    if (o%maximum_multiplicity>0.and.o%maximum_multiplicity<1e-12_dp) return
    if (o%maximum_parcels<1.or.o%maximum_parcels>1000000) return
    if (o%maximum_multiplicity>0.and.trim(o%breakup_model)=='none') return
    select case(trim(o%breakup_model))
    case('none')
    case('etab','khrt')
      if (.not.all(ieee_is_finite([o%surface_tension,o%liquid_viscosity]))) return
      if (o%surface_tension<1e-8_dp.or.o%surface_tension>1e3_dp) return
      if (trim(o%liquid_viscosity_model)=='constant') then
        if (o%liquid_viscosity<1e-10_dp.or.o%liquid_viscosity>1e3_dp) return
      end if
    case default
      return
    end select
    if (.not.all(ieee_is_finite(o%liquid_viscosity_coefficients))) return
    if (.not.all(ieee_is_finite(o%liquid_viscosity_temperature_range))) return
    select case(trim(o%liquid_viscosity_model))
    case('constant')
      if (any(o%liquid_viscosity_coefficients/=0).or.any(o%liquid_viscosity_temperature_range/=0)) return
    case('inverse-temperature')
      if (trim(o%breakup_model)=='none') return
      if (maxval(abs(o%liquid_viscosity_coefficients))>1e20_dp) return
      if (o%liquid_viscosity_temperature_range(1)<1.or.o%liquid_viscosity_temperature_range(2)>1e5_dp) return
      if (o%liquid_viscosity_temperature_range(2)<=o%liquid_viscosity_temperature_range(1)) return
    case default
      return
    end select
    if (trim(o%breakup_model)=='khrt') then
      if (.not.all(ieee_is_finite([o%khrt_b0,o%khrt_b1,o%khrt_c3]))) return
      if (min(o%khrt_b0,o%khrt_b1,o%khrt_c3)<1e-6_dp.or.max(o%khrt_b0,o%khrt_b1,o%khrt_c3)>1e3_dp) return
    end if
    ok = o%max_fractional_change > 0.0_dp .and. o%max_fractional_change <= 0.25_dp .and. &
      o%particle_cfl > 0.0_dp .and. o%particle_cfl <= 0.5_dp .and. &
      o%minimum_diameter >= 1.0e-12_dp .and. o%minimum_diameter <= 1.0_dp .and. o%maximum_substeps > 0
  end function

  pure real(dp) function liquid_viscosity_at_temperature(o,temperature) result(mu)
    type(spray_options), intent(in) :: o
    real(dp), intent(in) :: temperature
    real(dp) :: candidate
    mu=0.0_dp
    if (.not.valid_spray_options(o).or..not.ieee_is_finite(temperature)) return
    if (temperature<1.or.temperature>1e5_dp) return
    candidate=o%liquid_viscosity
    if (trim(o%liquid_viscosity_model)=='inverse-temperature') then
      if (temperature<o%liquid_viscosity_temperature_range(1).or. &
          temperature>o%liquid_viscosity_temperature_range(2)) return
      candidate=o%liquid_viscosity_coefficients(1)+((o%liquid_viscosity_coefficients(4)/temperature+ &
        o%liquid_viscosity_coefficients(3))/temperature+o%liquid_viscosity_coefficients(2))/temperature
    end if
    if (.not.ieee_is_finite(candidate)) return
    if (candidate<1e-10_dp.or.candidate>1e3_dp) return
    mu=candidate
  end function

  pure real(dp) function liquid_density_at_temperature(l,temperature) result(rho)
    type(spray_liquid), intent(in) :: l
    real(dp), intent(in) :: temperature
    real(dp) :: candidate
    rho=0.0_dp
    if (.not.valid_spray_liquid(l).or..not.ieee_is_finite(temperature)) return
    if (temperature<1.or.temperature>1e5_dp) return
    candidate=l%density
    if (trim(l%density_model)=='polynomial') then
      if (temperature<l%density_temperature_range(1).or.temperature>l%density_temperature_range(2)) return
      candidate=l%density_coefficients(1)+temperature*(l%density_coefficients(2)+ &
        temperature*(l%density_coefficients(3)+temperature*l%density_coefficients(4)))
    end if
    if (.not.ieee_is_finite(candidate)) return
    if (candidate<1e-3_dp.or.candidate>1e5_dp) return
    rho=candidate
  end function

  pure real(dp) function parcel_diameter(p, liquid) result(d)
    type(spray_parcel), intent(in) :: p
    type(spray_liquid), intent(in) :: liquid
    real(dp) :: rho
    d = 0.0_dp
    if (.not. valid_spray_liquid(liquid) .or. .not. valid_spray_parcel(p)) return
    rho=liquid_density_at_temperature(liquid,p%temperature)
    if (rho<=0) return
    if (p%mass > 0.0_dp) d = (6.0_dp*p%mass/(pi*rho))**(1.0_dp/3.0_dp)
  end function

  pure real(dp) function liquid_energy(l, temperature) result(e)
    type(spray_liquid), intent(in) :: l
    real(dp), intent(in) :: temperature
    e = l%reference_energy + l%heat_capacity*(temperature-l%reference_temperature)
  end function

  pure function parcel_totals(p, l) result(totals)
    type(spray_parcel), intent(in) :: p
    type(spray_liquid), intent(in) :: l
    real(dp) :: totals(5), m
    m = p%mass*p%multiplicity
    totals(1) = m
    totals(2:4) = m*p%velocity
    totals(5) = m*(liquid_energy(l,p%temperature)+0.5_dp*sum(p%velocity**2))
  end function

  subroutine saturation_mole_fraction(l, molecular_weight, temperature, pressure, fraction, ok)
    type(spray_liquid), intent(in) :: l
    real(dp), intent(in) :: molecular_weight, temperature, pressure
    real(dp), intent(out) :: fraction
    logical, intent(out) :: ok
    real(dp) :: log_fraction,surface_temperature,inverse_boiling_temperature,boiling_temperature
    fraction = 0.0_dp
    ok = .false.
    if (.not. valid_spray_liquid(l)) return
    if (.not. all(ieee_is_finite([molecular_weight, temperature, pressure]))) return
    if (molecular_weight < 1.0e-3_dp .or. molecular_weight > 1.0e5_dp) return
    if (temperature < 1.0_dp .or. temperature > 1.0e5_dp) return
    if (pressure < 1.0e-6_dp .or. pressure > 1.0e12_dp) return
    surface_temperature=temperature
    if (trim(l%boiling_model)=='surface-cap') then
      inverse_boiling_temperature=1/l%reference_temperature+ &
        log(l%reference_pressure/pressure)*universal_gas_constant/(l%latent_heat*molecular_weight)
      if (inverse_boiling_temperature<=0) return
      if (inverse_boiling_temperature<1/l%critical_temperature) then
        boiling_temperature=l%critical_temperature
      else
        boiling_temperature=1/inverse_boiling_temperature
      end if
      surface_temperature=min(temperature,.99_dp*boiling_temperature)
      if (surface_temperature<1) return
    end if
    if (trim(l%vapor_pressure_model)=='antoine') then
      if (temperature<l%antoine_temperature_range(1).or.temperature>l%antoine_temperature_range(2)) return
      ! PelePhysics convention: psat = D*10**(A-B/(T+C)), with T in K.
      log_fraction=log(l%antoine_coefficients(4)/pressure)+log(10.0_dp)* &
        (l%antoine_coefficients(1)-l%antoine_coefficients(2)/(temperature+l%antoine_coefficients(3)))
    else
    log_fraction = log(l%reference_pressure/pressure) + l%latent_heat*molecular_weight/universal_gas_constant * &
      (1.0_dp/l%reference_temperature-1.0_dp/surface_temperature)
    end if
    if (.not.ieee_is_finite(log_fraction)) return
    ! Reject saturation unless the explicitly selected surface cap avoided it.
    if (log_fraction >= 0.0_dp) return
    fraction = exp(max(log(tiny(1.0_dp)),log_fraction))
    ok = fraction < 1.0_dp
  end subroutine

  subroutine spray_transfer_rates(p, l, f, o, inverse_tau, heat_coefficient, mass_coefficient, ok)
    type(spray_parcel), intent(in) :: p
    type(spray_liquid), intent(in) :: l
    type(spray_film), intent(in) :: f
    type(spray_options), intent(in) :: o
    real(dp), intent(out) :: inverse_tau, heat_coefficient, mass_coefficient
    logical, intent(out) :: ok
    real(dp) :: rho
    real(dp) :: d, re, cdre, pr, sc, bm, bt, logbm, sh0, nu0, sh, nu, next_nu, phi, factor, f2
    integer :: iteration
    inverse_tau=0.0_dp; heat_coefficient=0.0_dp; mass_coefficient=0.0_dp
    ok = .false.
    if (.not. valid_spray_liquid(l) .or. .not. valid_spray_parcel(p)) return
    if (.not. valid_spray_options(o)) return
    if (.not. all(ieee_is_finite([f%density,f%viscosity,f%conductivity,f%diffusivity,f%cp, &
        f%gas_temperature,f%gas_velocity,f%surface_fraction,f%ambient_fraction,f%vapor_enthalpy]))) return
    if (f%density <= 0.0_dp .or. f%density > 1.0e8_dp) return
    if (f%viscosity < 1.0e-15_dp .or. f%viscosity > 1.0e5_dp) return
    if (f%conductivity < 1.0e-15_dp .or. f%conductivity > 1.0e8_dp) return
    if (f%diffusivity < 1.0e-20_dp .or. f%diffusivity > 1.0e8_dp) return
    if (f%cp <= 0.0_dp .or. f%cp > 1.0e10_dp) return
    if (f%gas_temperature < 1.0_dp .or. f%gas_temperature > 1.0e5_dp) return
    if (maxval(abs(f%gas_velocity)) > 1.0e7_dp .or. abs(f%vapor_enthalpy) > 1.0e12_dp) return
    if (f%surface_fraction < 0.0_dp .or. f%surface_fraction >= 1.0_dp) return
    if (f%ambient_fraction < 0.0_dp .or. f%ambient_fraction > 1.0_dp) return
    if (p%mass == 0.0_dp) then
      ok=.true.; return
    end if
    rho=liquid_density_at_temperature(l,p%temperature)
    if (rho<=0) return
    d=parcel_diameter(p,l)
    if (d < 1.0e-12_dp) return
    re=f%density*d*norm2(f%gas_velocity-p%velocity)/f%viscosity
    cdre=24.0_dp*(1.0_dp+0.15_dp*re**0.687_dp)
    if (re >= 1000.0_dp) cdre=0.44_dp*re
    pr=f%viscosity*f%cp/f%conductivity
    sc=f%viscosity/(f%density*f%diffusivity)
    sh0=2.0_dp+0.6_dp*sqrt(re)*sc**(1.0_dp/3.0_dp)
    nu0=2.0_dp+0.6_dp*sqrt(re)*pr**(1.0_dp/3.0_dp)
    if (o%abramzon_sirignano) then
      f2=max(1.0_dp,min(400.0_dp,re)**0.077_dp)
      sh0=1.0_dp+(1.0_dp+re*sc)**(1.0_dp/3.0_dp)*f2
      nu0=1.0_dp+(1.0_dp+re*pr)**(1.0_dp/3.0_dp)*f2
    end if
    bm=0.0_dp
    if (o%evaporation) bm=max(0.0_dp,(f%surface_fraction-f%ambient_fraction)/(1.0_dp-f%surface_fraction))
    logbm=log(1.0_dp+bm)
    sh=sh0; nu=nu0; factor=1.0_dp
    if (o%abramzon_sirignano .and. bm > 1.0e-12_dp) then
      sh=2.0_dp+(sh0-2.0_dp)/blowing_factor(bm)
      do iteration=1,100
        phi=f%cp*f%density*f%diffusivity*sh/(f%conductivity*nu)
        if (phi*logbm > 100.0_dp) return
        bt=exp(phi*logbm)-1.0_dp
        next_nu=2.0_dp+(nu0-2.0_dp)/blowing_factor(bt)
        if (abs(next_nu-nu) <= 1.0e-10_dp*max(1.0_dp,nu)) exit
        nu=0.5_dp*(nu+next_nu)
      end do
      if (iteration > 100) return
      nu=next_nu
      if (bt > 1.0e-12_dp) factor=log(1.0_dp+bt)/bt
    end if
    if (o%drag) inverse_tau=3.0_dp*f%viscosity*cdre/(4.0_dp*rho*d*d)
    if (o%heat) heat_coefficient=pi*f%conductivity*d*nu*factor
    if (o%evaporation) mass_coefficient=pi*f%density*f%diffusivity*sh*logbm
    ok=.true.
  end subroutine

  pure real(dp) function blowing_factor(b) result(f)
    real(dp), intent(in) :: b
    f=1.0_dp
    if (b > 1.0e-12_dp) f=(1.0_dp+b)**0.7_dp*log(1.0_dp+b)/b
  end function

  subroutine advance_parcel_frozen_film(p, l, f, o, dt, q, transfer, ok)
    type(spray_parcel), intent(in) :: p
    type(spray_liquid), intent(in) :: l
    type(spray_film), intent(in) :: f
    type(spray_options), intent(in) :: o
    real(dp), intent(in) :: dt
    type(spray_parcel), intent(out) :: q
    real(dp), intent(out) :: transfer(5)
    logical, intent(out) :: ok
    real(dp) :: rho
    real(dp) :: invtau, heat, masscoef, d, d2, dm, mbar, decay, latent, rate, equilibrium
    q=p; transfer=0.0_dp; ok=.false.
    if (.not. ieee_is_finite(dt)) return
    if (dt < 0.0_dp .or. dt > 1.0e6_dp) return
    call spray_transfer_rates(p,l,f,o,invtau,heat,masscoef,ok)
    if (.not. ok) return
    if (dt == 0.0_dp .or. p%mass == 0.0_dp) return
    ok=.false.
    rho=liquid_density_at_temperature(l,p%temperature)
    if (rho<=0) return
    d=parcel_diameter(p,l)
    d2=max(0.0_dp,d*d-4.0_dp*masscoef*dt/(pi*rho))
    q%mass=min(p%mass,(pi*rho/6.0_dp)*d2**1.5_dp)
    ! Avoid drift when no mass transfer is requested.
    if (masscoef == 0.0_dp) q%mass=p%mass
    if (o%evaporation .and. masscoef > 0.0_dp .and. sqrt(d2) < o%minimum_diameter) q%mass=0.0_dp
    q%velocity=f%gas_velocity+(p%velocity-f%gas_velocity)*exp(-min(700.0_dp,invtau*dt))
    if (.not. o%drag) q%velocity=p%velocity
    q%position=p%position+0.5_dp*dt*(p%velocity+q%velocity)
    if (q%mass > 0.0_dp) then
      dm=p%mass-q%mass
      mbar=0.5_dp*(p%mass+q%mass)
      latent=f%vapor_enthalpy-liquid_energy(l,p%temperature)
      if (masscoef > 0.0_dp .and. latent <= 0.0_dp) then
        q=p; return
      end if
      rate=dm/dt
      if (heat > 0.0_dp) then
        equilibrium=f%gas_temperature-rate*latent/heat
        decay=exp(-min(700.0_dp,heat*dt/(mbar*l%heat_capacity)))
        q%temperature=equilibrium+(p%temperature-equilibrium)*decay
      else
        q%temperature=p%temperature-dm*latent/(mbar*l%heat_capacity)
      end if
    end if
    if (.not. valid_spray_parcel(q)) then
      q=p; return
    end if
    if (q%mass>0) then
      if (liquid_density_at_temperature(l,q%temperature)<=0) then
        q=p; return
      end if
    end if
    ! Discrete equal-and-opposite exchange includes evaporated momentum and
    ! kinetic energy. Do not add separate latent heat or drag-work terms again.
    transfer=parcel_totals(p,l)-parcel_totals(q,l)
    ok=all(ieee_is_finite(transfer))
    if (.not. ok) then
      q=p; transfer=0.0_dp
    end if
  end subroutine

  subroutine inject_cone_ports(liquid,total_mass,number,first_id,centers,axes,fractions,speed, &
      half_angle,diameter,temperature,parcels,ok,diameter_shape,diameter_min,diameter_max)
    type(spray_liquid), intent(in) :: liquid
    real(dp), intent(in) :: total_mass,centers(:,:),axes(:,:),fractions(:),speed,half_angle,diameter,temperature
    integer, intent(in) :: number
    integer(int64), intent(in) :: first_id
    type(spray_parcel), allocatable, intent(out) :: parcels(:)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: diameter_shape,diameter_min,diameter_max
    type(spray_parcel), allocatable :: candidate(:),port_parcels(:)
    integer :: ports,port,offset,total
    real(dp) :: fraction_sum
    ok=.false.;ports=size(fractions)
    if (ports<1.or.ports>16.or.number<1) return
    if (number>1000000/ports) return
    if (size(centers,1)/=3.or.size(axes,1)/=3) return
    if (size(centers,2)/=ports.or.size(axes,2)/=ports) return
    if (.not.all(ieee_is_finite(fractions))) return
    if (any(fractions<=0).or.any(fractions>1)) return
    fraction_sum=sum(fractions)
    if (abs(fraction_sum-1.0_dp)>5e-12_dp) return
    total=ports*number
    if (first_id<1_int64.or.first_id>huge(first_id)-int(total,int64)) return
    allocate(candidate(total))
    do port=1,ports
      offset=(port-1)*number
      call inject_cone_parcels(liquid,total_mass*(fractions(port)/fraction_sum),number,first_id+int(offset,int64), &
        centers(:,port),axes(:,port),speed,half_angle,diameter,temperature,port_parcels,ok, &
        diameter_shape,diameter_min,diameter_max)
      if (.not.ok) return
      candidate(offset+1:offset+number)=port_parcels
      deallocate(port_parcels)
    end do
    call move_alloc(candidate,parcels)
  end subroutine inject_cone_ports

  subroutine inject_cone_parcels(liquid, total_mass, number, first_id, center, axis, speed, &
      half_angle, diameter, temperature, parcels, ok, diameter_shape, diameter_min, diameter_max)
    type(spray_liquid), intent(in) :: liquid
    real(dp), intent(in) :: total_mass, center(3), axis(3), speed, half_angle, diameter, temperature
    integer, intent(in) :: number
    integer(int64), intent(in) :: first_id
    type(spray_parcel), allocatable, intent(out) :: parcels(:)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: diameter_shape, diameter_min, diameter_max
    real(dp) :: rho
    real(dp) :: e1(3),e2(3),ez(3), helper(3), cosine, sine, angle, m, weight
    real(dp) :: shape,dlo,dhi,survival_low,survival_high,survival,quantile,drop_size
    integer :: i
    ok=.false.
    if (.not. valid_spray_liquid(liquid)) return
    if (.not. all(ieee_is_finite([total_mass,center,axis,speed,half_angle,diameter,temperature]))) return
    if (number < 1 .or. number > 1000000) return
    if (first_id < 1_int64 .or. first_id > huge(first_id)-int(number,int64)) return
    if (total_mass <= 0.0_dp .or. total_mass > 1.0e6_dp) return
    if (maxval(abs(center)) > 1.0e9_dp .or. maxval(abs(axis)) > 1.0e9_dp) return
    if (norm2(axis) <= 1.0e-12_dp .or. speed < 0.0_dp .or. speed > 1.0e7_dp) return
    if (diameter < 1.0e-12_dp .or. diameter > 1.0_dp) return
    if (half_angle < 0.0_dp .or. half_angle > 0.5_dp*pi) return
    if (temperature < 1.0_dp .or. temperature > 1.0e5_dp) return
    rho=liquid_density_at_temperature(liquid,temperature)
    if (rho<=0) return
    shape=0.0_dp; dlo=1.0e-12_dp; dhi=1.0_dp
    if (present(diameter_shape)) shape=diameter_shape
    if (present(diameter_min)) dlo=diameter_min
    if (present(diameter_max)) dhi=diameter_max
    if (.not. all(ieee_is_finite([shape,dlo,dhi]))) return
    if (shape<0.0_dp .or. shape>20.0_dp) return
    if (dlo<1.0e-12_dp .or. dhi>1.0_dp .or. dlo>=dhi) return
    survival_low=1.0_dp; survival_high=0.0_dp
    if (shape>0.0_dp) then
      if (shape<0.25_dp) return
      survival_low=exp(-min(700.0_dp,(dlo/diameter)**shape))
      survival_high=exp(-min(700.0_dp,(dhi/diameter)**shape))
      if (survival_low-survival_high<=32.0_dp*epsilon(1.0_dp)*survival_low) return
    end if
    ez=axis/norm2(axis); helper=[1.0_dp,0.0_dp,0.0_dp]
    if (abs(ez(1)) > 0.9_dp) helper=[0.0_dp,1.0_dp,0.0_dp]
    e1=helper-dot_product(helper,ez)*ez; e1=e1/norm2(e1)
    e2=[ez(2)*e1(3)-ez(3)*e1(2),ez(3)*e1(1)-ez(1)*e1(3),ez(1)*e1(2)-ez(2)*e1(1)]
    allocate(parcels(number))
    do i=1,number
      drop_size=diameter
      if (shape>0.0_dp) then
        ! The specified truncated Weibull CDF is MASS weighted. Equal parcel
        ! mass and varying multiplicity reproduce it, not a number-weighted CDF.
        quantile=(real(i,dp)-0.5_dp)/real(number,dp)
        survival=survival_low-quantile*(survival_low-survival_high)
        drop_size=diameter*(-log(survival))**(1.0_dp/shape)
        if (drop_size<dlo .or. drop_size>dhi) then
          deallocate(parcels); return
        end if
      end if
      m=pi*rho*drop_size**3/6.0_dp
      weight=total_mass/(real(number,dp)*m)
      if (weight>1.0e18_dp .or. m>1.0_dp) then
        deallocate(parcels); return
      end if
      ! Deterministic equal-solid-angle quadrature, not a stochastic turbulence model.
      cosine=1.0_dp-(1.0_dp-cos(half_angle))*(real(i,dp)-0.5_dp)/real(number,dp)
      sine=sqrt(max(0.0_dp,1.0_dp-cosine*cosine))
      angle=2.0_dp*pi*modulo(real(first_id-1_int64+int(i,int64),dp)*0.6180339887498949_dp,1.0_dp)
      parcels(i)=spray_parcel(first_id+int(i-1,int64),center, &
        speed*(cosine*ez+sine*(cos(angle)*e1+sin(angle)*e2)),m,temperature,weight)
      if (.not. valid_spray_parcel(parcels(i))) then
        deallocate(parcels); return
      end if
    end do
    ok=.true.
  end subroutine
end module spray_parcel_mod
