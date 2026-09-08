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
  end type
  type, public :: spray_parcel
    integer(int64) :: id = 0_int64
    real(dp) :: position(3) = 0.0_dp, velocity(3) = 0.0_dp
    real(dp) :: mass = 0.0_dp, temperature = 0.0_dp, multiplicity = 0.0_dp
  end type
  type, public :: spray_options
    logical :: drag = .true., heat = .true., evaporation = .true.
    logical :: abramzon_sirignano = .true.
    real(dp) :: max_fractional_change = 0.1_dp, particle_cfl = 0.25_dp
    real(dp) :: minimum_diameter = 1.0e-9_dp
    integer :: maximum_substeps = 100000
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
  public :: inject_cone_parcels
contains
  subroutine initialize_spray_liquid(species, name, density, cp, tref, pref, latent, liquid, ok)
    type(nasa7_species), intent(in) :: species(:)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: density, cp, tref, pref, latent
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
    ok = valid_spray_liquid(liquid)
  end subroutine

  pure logical function valid_spray_liquid(l) result(ok)
    type(spray_liquid), intent(in) :: l
    ok = .false.
    if (.not. all(ieee_is_finite([l%density, l%heat_capacity, l%reference_temperature, &
        l%reference_pressure, l%latent_heat, l%reference_energy]))) return
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
    if (.not. all(ieee_is_finite([p%position, p%velocity, p%mass, p%temperature, p%multiplicity]))) return
    ok = p%id > 0_int64 .and. p%mass >= 0.0_dp .and. p%mass <= 1.0_dp .and. &
      p%temperature >= 1.0_dp .and. p%temperature <= 1.0e5_dp .and. &
      p%multiplicity > 0.0_dp .and. p%multiplicity <= 1.0e18_dp .and. &
      maxval(abs(p%position)) <= 1.0e9_dp .and. maxval(abs(p%velocity)) <= 1.0e7_dp
  end function

  pure logical function valid_spray_options(o) result(ok)
    type(spray_options), intent(in) :: o
    ok = .false.
    if (.not. all(ieee_is_finite([o%max_fractional_change,o%particle_cfl,o%minimum_diameter]))) return
    ok = o%max_fractional_change > 0.0_dp .and. o%max_fractional_change <= 0.25_dp .and. &
      o%particle_cfl > 0.0_dp .and. o%particle_cfl <= 0.5_dp .and. &
      o%minimum_diameter >= 1.0e-12_dp .and. o%minimum_diameter <= 1.0_dp .and. o%maximum_substeps > 0
  end function

  pure real(dp) function parcel_diameter(p, liquid) result(d)
    type(spray_parcel), intent(in) :: p
    type(spray_liquid), intent(in) :: liquid
    d = 0.0_dp
    if (.not. valid_spray_liquid(liquid) .or. .not. valid_spray_parcel(p)) return
    if (p%mass > 0.0_dp) d = (6.0_dp*p%mass/(pi*liquid%density))**(1.0_dp/3.0_dp)
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
    real(dp) :: log_fraction
    fraction = 0.0_dp
    ok = .false.
    if (.not. valid_spray_liquid(l)) return
    if (.not. all(ieee_is_finite([molecular_weight, temperature, pressure]))) return
    if (molecular_weight < 1.0e-3_dp .or. molecular_weight > 1.0e5_dp) return
    if (temperature < 1.0_dp .or. temperature > 1.0e5_dp) return
    if (pressure < 1.0e-6_dp .or. pressure > 1.0e12_dp) return
    log_fraction = log(l%reference_pressure/pressure) + l%latent_heat*molecular_weight/universal_gas_constant * &
      (1.0_dp/l%reference_temperature-1.0_dp/temperature)
    ! This sub-boiling model rejects boiling/supercritical states, not clamps them.
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
    if (o%drag) inverse_tau=3.0_dp*f%viscosity*cdre/(4.0_dp*l%density*d*d)
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
    real(dp) :: invtau, heat, masscoef, d, d2, dm, mbar, decay, latent, rate, equilibrium
    q=p; transfer=0.0_dp; ok=.false.
    if (.not. ieee_is_finite(dt)) return
    if (dt < 0.0_dp .or. dt > 1.0e6_dp) return
    call spray_transfer_rates(p,l,f,o,invtau,heat,masscoef,ok)
    if (.not. ok) return
    if (dt == 0.0_dp .or. p%mass == 0.0_dp) return
    ok=.false.
    d=parcel_diameter(p,l)
    d2=max(0.0_dp,d*d-4.0_dp*masscoef*dt/(pi*l%density))
    q%mass=min(p%mass,(pi*l%density/6.0_dp)*d2**1.5_dp)
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
    ! Discrete equal-and-opposite exchange includes evaporated momentum and
    ! kinetic energy. Do not add separate latent heat or drag-work terms again.
    transfer=parcel_totals(p,l)-parcel_totals(q,l)
    ok=all(ieee_is_finite(transfer))
    if (.not. ok) then
      q=p; transfer=0.0_dp
    end if
  end subroutine

  subroutine inject_cone_parcels(liquid, total_mass, number, first_id, center, axis, speed, &
      half_angle, diameter, temperature, parcels, ok, diameter_shape, diameter_min, diameter_max)
    type(spray_liquid), intent(in) :: liquid
    real(dp), intent(in) :: total_mass, center(3), axis(3), speed, half_angle, diameter, temperature
    integer, intent(in) :: number
    integer(int64), intent(in) :: first_id
    type(spray_parcel), allocatable, intent(out) :: parcels(:)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: diameter_shape, diameter_min, diameter_max
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
      m=pi*liquid%density*drop_size**3/6.0_dp
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
