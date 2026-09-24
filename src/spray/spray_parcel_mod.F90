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
    character(len=16) :: heat_capacity_model='constant'
    real(dp) :: heat_capacity_coefficients(4)=0.0_dp,heat_capacity_temperature_range(2)=0.0_dp
    real(dp) :: normal_boiling_temperature=0.0_dp
    real(dp) :: molecular_weight=0.0_dp ! zero uses the mapped gas species
  end type
  type, public :: spray_nozzle_model
    character(len=16) :: profile='legacy'
    real(dp) :: swirl=0,spread=0
  end type
  type, public :: spray_size_options
    character(len=16) :: distribution='weibull',weighting='mass'
    real(dp) :: std=0,shape=0,minimum=1e-12_dp,maximum=1
  end type
  type, public :: spray_parcel
    integer(int64) :: id = 0_int64
    real(dp) :: position(3) = 0.0_dp, velocity(3) = 0.0_dp
    real(dp) :: mass = 0.0_dp, temperature = 0.0_dp, multiplicity = 0.0_dp
    real(dp) :: deformation=0.0_dp,deformation_rate=0.0_dp
    integer(int64) :: parent_id=0_int64
    real(dp) :: rt_time=-1.0_dp
    integer :: material_id=1
    real(dp) :: khrt_initial_count=0,khrt_shed_volume=0
  end type
  type, public :: spray_options
    logical :: khrt_children=.false.
    real(dp) :: khrt_reference_diameter_cubed=0
    real(dp) :: khrt_b0=.61_dp,khrt_b1=7.0_dp,khrt_c3=1.0_dp
    real(dp) :: maximum_multiplicity=0.0_dp
    integer :: maximum_parcels=1000000
    character(len=8) :: breakup_model='none'
    real(dp) :: surface_tension=0.0_dp,liquid_viscosity=0.0_dp
    character(len=24) :: drag_model="schiller-naumann"
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
    real(dp) :: vapor_cp = 0.0_dp ! zero retains the supplied-film cp fallback
  end type
  public :: initialize_spray_liquid, valid_spray_liquid, valid_spray_parcel
  public :: valid_spray_options, parcel_diameter, liquid_energy, parcel_totals
  public :: spray_transfer_rates, advance_parcel_frozen_film, saturation_mole_fraction
  public :: liquid_heat_capacity_at_temperature
  public :: watson_latent_heat, pelemp_boiling_temperature,resolved_liquid_molecular_weight
  public :: select_spray_liquid
  public :: inject_cone_ports
  public :: inject_cone_parcels, liquid_density_at_temperature, liquid_viscosity_at_temperature
contains
  pure real(dp) function resolved_liquid_molecular_weight(liquid,gas_molecular_weight) result(weight)
    type(spray_liquid), intent(in) :: liquid
    real(dp), intent(in) :: gas_molecular_weight
    weight=gas_molecular_weight
    if (liquid%molecular_weight/=0) weight=liquid%molecular_weight
  end function

  pure subroutine watson_latent_heat(latent_reference, reference_temperature, critical_temperature, &
      temperature, latent, ok)
    ! PelePhysics MPLiqProps::latent_correlation_i; SI, subcritical states only.
    ! The reference temperature belongs to the supplied latent heat, not to
    ! the independently specified normal boiling point.
    real(dp), intent(in) :: latent_reference,reference_temperature,critical_temperature,temperature
    real(dp), intent(out) :: latent
    logical, intent(out) :: ok
    real(dp) :: candidate
    latent=0;ok=.false.
    if (.not.all(ieee_is_finite([latent_reference,reference_temperature,critical_temperature,temperature]))) return
    if (latent_reference<=0.or.latent_reference>1e10_dp) return
    if (reference_temperature<1.or.reference_temperature>=critical_temperature) return
    if (critical_temperature>1e5_dp) return
    if (temperature<1.or.temperature>=critical_temperature) return
    candidate=latent_reference*((critical_temperature-reference_temperature)/ &
      (critical_temperature-temperature))**(-.38_dp)
    if (.not.ieee_is_finite(candidate)) return
    if (candidate<=0) return
    latent=candidate;ok=.true.
  end subroutine

  pure subroutine pelemp_boiling_temperature(latent_reference, reference_temperature, normal_boiling_temperature, &
      critical_temperature, molecular_weight, pressure, boiling_temperature, ok)
    ! PelePhysics SprayData::calcBoilT with MPLiqProps. Normal boiling point is
    ! at 101325 Pa; molecular weight is kg/kmol and latent heat is J/kg.
    ! Reject nonphysical inverse temperatures instead of returning a negative
    ! boiling point. This does not change stored liquid energy or vapor pressure.
    real(dp), intent(in) :: latent_reference,reference_temperature,normal_boiling_temperature,critical_temperature
    real(dp), intent(in) :: molecular_weight,pressure
    real(dp), intent(out) :: boiling_temperature
    logical, intent(out) :: ok
    real(dp) :: latent,denominator,numerator,inverse_temperature,candidate
    logical :: valid
    boiling_temperature=0;ok=.false.
    if (.not.all(ieee_is_finite([molecular_weight,pressure]))) return
    if (molecular_weight<1e-3_dp.or.molecular_weight>1e5_dp) return
    if (pressure<1e-6_dp.or.pressure>1e12_dp) return
    call watson_latent_heat(latent_reference,reference_temperature,critical_temperature, &
      normal_boiling_temperature,latent,valid)
    if (.not.valid) return
    denominator=latent*molecular_weight
    numerator=log(101325.0_dp/pressure)*universal_gas_constant
    if (denominator<=0) return
    if (denominator<abs(numerator)/huge(1.0_dp)) return
    inverse_temperature=1/normal_boiling_temperature+numerator/denominator
    if (.not.ieee_is_finite(inverse_temperature)) return
    if (inverse_temperature<=0.or.inverse_temperature>1) return
    if (inverse_temperature<1/critical_temperature) then
      candidate=critical_temperature
    else
      candidate=1/inverse_temperature
    end if
    boiling_temperature=candidate;ok=.true.
  end subroutine

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
    real(dp) :: a,b,c,disc,t,roots(2)
    integer :: i
    ok = .false.
    if (.not. all(ieee_is_finite([l%density, l%heat_capacity, l%reference_temperature, &
        l%reference_pressure, l%latent_heat, l%reference_energy]))) return
    if (.not.ieee_is_finite(l%molecular_weight)) return
    if (l%molecular_weight/=0) then
      if (l%molecular_weight<1e-3_dp.or.l%molecular_weight>1e5_dp) return
    end if
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
    if (.not.all(ieee_is_finite([l%critical_temperature,l%normal_boiling_temperature]))) return
    select case(trim(l%boiling_model))
    case('reject')
      if (l%critical_temperature/=0.or.l%normal_boiling_temperature/=0) return
    case('surface-cap')
      if (l%normal_boiling_temperature/=0) return
      if (trim(l%vapor_pressure_model)/='clausius') return
      if (l%critical_temperature<=l%reference_temperature.or.l%critical_temperature>1e5_dp) return
    case('pelemp-cap')
      if (l%critical_temperature<=l%reference_temperature.or.l%critical_temperature>1e5_dp) return
      if (l%normal_boiling_temperature<1.or.l%normal_boiling_temperature>=l%critical_temperature) return
      if (l%reference_pressure/=101325.0_dp) return
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
    if (.not.all(ieee_is_finite(l%heat_capacity_coefficients)).or. &
        .not.all(ieee_is_finite(l%heat_capacity_temperature_range))) return
    select case(trim(l%heat_capacity_model))
    case('constant')
      if (any(l%heat_capacity_coefficients/=0).or.any(l%heat_capacity_temperature_range/=0)) return
    case('polynomial')
      if (maxval(abs(l%heat_capacity_coefficients))>1e10_dp) return
      if (l%heat_capacity_temperature_range(1)<1.or.l%heat_capacity_temperature_range(2)>1e5_dp) return
      if (l%heat_capacity_temperature_range(2)<=l%heat_capacity_temperature_range(1)) return
      if (l%reference_temperature<l%heat_capacity_temperature_range(1).or. &
          l%reference_temperature>l%heat_capacity_temperature_range(2)) return
      do i=1,2
        if (liquid_heat_capacity_at_temperature(l,l%heat_capacity_temperature_range(i))<=0) return
      end do
      ! Check interior extrema as well as endpoints so the integrated energy is monotone.
      a=3*l%heat_capacity_coefficients(4);b=2*l%heat_capacity_coefficients(3);c=l%heat_capacity_coefficients(2)
      roots=l%heat_capacity_temperature_range(1)
      if (a/=0) then
        disc=b*b-4*a*c
        if (disc>=0) then
          t=-0.5_dp*(b+sign(sqrt(disc),b))
          if (abs(t)<=abs(a)*1e5_dp) roots(1)=t/a
          if (t/=0) then
            if (abs(c)<=abs(t)*1e5_dp) roots(2)=c/t
          end if
        end if
      else if (b/=0) then
        if (abs(c)<=abs(b)*1e5_dp) roots(1)=-c/b
      end if
      do i=1,2
        if (roots(i)<l%heat_capacity_temperature_range(1).or.roots(i)>l%heat_capacity_temperature_range(2)) cycle
        if (liquid_heat_capacity_at_temperature(l,roots(i))<=0) return
      end do
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
    if (.not.all(ieee_is_finite([p%khrt_initial_count,p%khrt_shed_volume]))) return
    if (p%khrt_initial_count<0.or.p%khrt_initial_count>1e18_dp.or.p%khrt_shed_volume<0.or.p%khrt_shed_volume>1e24_dp) return
    if (p%khrt_initial_count==0.and.p%khrt_shed_volume/=0) return
    if (p%material_id<1.or.p%material_id>16) return
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
    if (.not.ieee_is_finite(o%khrt_reference_diameter_cubed)) return
    if (o%khrt_children) then
      if (o%breakup_model/='khrt') return
      if (o%khrt_reference_diameter_cubed<=0.or.o%khrt_reference_diameter_cubed>1e24_dp) return
    else
      if (o%khrt_reference_diameter_cubed/=0) return
    end if
    if (o%drag_model/="schiller-naumann".and.o%drag_model/="pelephysics") return
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

  function select_spray_liquid(parcel,liquid,materials) result(material)
    type(spray_parcel), intent(in) :: parcel
    type(spray_liquid), intent(in) :: liquid
    type(spray_liquid), intent(in), optional :: materials(:)
    type(spray_liquid) :: material
    material=spray_liquid()
    if (present(materials)) then
      if (size(materials)<1.or.size(materials)>16) return
      if (parcel%material_id<1.or.parcel%material_id>size(materials)) return
      material=materials(parcel%material_id)
    else
      if (parcel%material_id/=1) return
      material=liquid
    end if
  end function select_spray_liquid

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

  pure real(dp) function liquid_heat_capacity_at_temperature(l,temperature) result(cp)
    type(spray_liquid), intent(in) :: l
    real(dp), intent(in) :: temperature
    cp=0
    if (.not.ieee_is_finite(temperature)) return
    if (temperature<1.or.temperature>1e5_dp) return
    select case(trim(l%heat_capacity_model))
    case('constant')
      cp=l%heat_capacity
    case('polynomial')
      if (temperature<l%heat_capacity_temperature_range(1).or.temperature>l%heat_capacity_temperature_range(2)) return
      cp=l%heat_capacity_coefficients(1)+temperature*(l%heat_capacity_coefficients(2)+ &
        temperature*(l%heat_capacity_coefficients(3)+temperature*l%heat_capacity_coefficients(4)))
    end select
    if (.not.ieee_is_finite(cp)) then
      cp=0
    else if (cp<1.or.cp>1e8_dp) then
      cp=0
    end if
  end function

  pure real(dp) function liquid_energy(l, temperature) result(e)
    type(spray_liquid), intent(in) :: l
    real(dp), intent(in) :: temperature
    real(dp) :: t,r
    e = l%reference_energy + l%heat_capacity*(temperature-l%reference_temperature)
    if (l%heat_capacity_model=='polynomial') then
      t=temperature;r=l%reference_temperature
      e=l%reference_energy+(t-r)*(l%heat_capacity_coefficients(1)+ &
        l%heat_capacity_coefficients(2)*(t+r)/2+ &
        l%heat_capacity_coefficients(3)*(t*t+t*r+r*r)/3+ &
        l%heat_capacity_coefficients(4)*(t+r)*(t*t+r*r)/4)
    end if
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

  subroutine saturation_mole_fraction(l, molecular_weight, temperature, pressure, fraction, ok, &
      liquid_mole_fraction,vapor_enthalpy)
    type(spray_liquid), intent(in) :: l
    real(dp), intent(in) :: molecular_weight, temperature, pressure
    real(dp), intent(out) :: fraction
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: liquid_mole_fraction,vapor_enthalpy
    real(dp) :: log_fraction,surface_temperature,inverse_boiling_temperature,boiling_temperature,liquid_x
    real(dp) :: latent,boiling_reference,pressure_temperature,liquid_weight
    logical :: valid
    fraction = 0.0_dp
    ok = .false.
    if (.not. valid_spray_liquid(l)) return
    if (.not. all(ieee_is_finite([molecular_weight, temperature, pressure]))) return
    if (molecular_weight < 1.0e-3_dp .or. molecular_weight > 1.0e5_dp) return
    if (temperature < 1.0_dp .or. temperature > 1.0e5_dp) return
    if (pressure < 1.0e-6_dp .or. pressure > 1.0e12_dp) return
    liquid_weight=resolved_liquid_molecular_weight(l,molecular_weight)
    liquid_x=1
    if (present(liquid_mole_fraction)) liquid_x=liquid_mole_fraction
    if (.not.ieee_is_finite(liquid_x)) return
    if (liquid_x<0.or.liquid_x>1) return
    if (liquid_x==0) then
      ok=.true.;return
    end if
    surface_temperature=temperature
    latent=l%latent_heat;boiling_reference=l%reference_temperature
    if (trim(l%boiling_model)=='surface-cap') then
      inverse_boiling_temperature=1/l%reference_temperature+ &
        log(l%reference_pressure/pressure)*universal_gas_constant/(l%latent_heat*liquid_weight)
      if (inverse_boiling_temperature<=0) return
      if (inverse_boiling_temperature<1/l%critical_temperature) then
        boiling_temperature=l%critical_temperature
      else
        boiling_temperature=1/inverse_boiling_temperature
      end if
      surface_temperature=min(temperature,.99_dp*boiling_temperature)
      if (surface_temperature<1) return
    end if
    pressure_temperature=temperature
    if (l%boiling_model=='pelemp-cap') then
      call pelemp_boiling_temperature(l%latent_heat,l%reference_temperature,l%normal_boiling_temperature, &
        l%critical_temperature,liquid_weight,pressure,boiling_temperature,valid)
      if (.not.valid) return
      surface_temperature=min(temperature,.99_dp*boiling_temperature)
      if (surface_temperature<1) return
      pressure_temperature=surface_temperature;boiling_reference=l%normal_boiling_temperature
      if (l%vapor_pressure_model=='clausius') then
        if (.not.present(vapor_enthalpy)) return
        if (.not.ieee_is_finite(vapor_enthalpy)) return
        ! Ordinary EOS reference: h_gas(Tdrop) minus liquid enthalpy at Ts.
        if (liquid_heat_capacity_at_temperature(l,surface_temperature)<=0) return
        latent=vapor_enthalpy-liquid_energy(l,surface_temperature)
        if (.not.ieee_is_finite(latent)) return
        if (latent<=0.or.latent>1e10_dp) return
      end if
    end if
    if (trim(l%vapor_pressure_model)=='antoine') then
      if (pressure_temperature<l%antoine_temperature_range(1).or.pressure_temperature>l%antoine_temperature_range(2)) return
      ! PelePhysics convention: psat = D*10**(A-B/(T+C)), with T in K.
      log_fraction=log(l%antoine_coefficients(4)/pressure)+log(10.0_dp)* &
        (l%antoine_coefficients(1)-l%antoine_coefficients(2)/(pressure_temperature+l%antoine_coefficients(3)))
    else
    log_fraction = log(l%reference_pressure/pressure) + latent*liquid_weight/universal_gas_constant * &
      (1.0_dp/boiling_reference-1.0_dp/surface_temperature)
    end if
    if (liquid_x/=1) log_fraction=log_fraction+log(liquid_x)
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
    real(dp) :: d, re, cdre, pr, sc, bm, bt, old_bt, logbm, sh0, nu0, sh, nu, phi, factor, f2,heat_cp
    integer :: iteration
    inverse_tau=0.0_dp; heat_coefficient=0.0_dp; mass_coefficient=0.0_dp
    ok = .false.
    if (.not. valid_spray_liquid(l) .or. .not. valid_spray_parcel(p)) return
    if (.not. valid_spray_options(o)) return
    if (.not. all(ieee_is_finite([f%density,f%viscosity,f%conductivity,f%diffusivity,f%cp, &
        f%gas_temperature,f%gas_velocity,f%surface_fraction,f%ambient_fraction,f%vapor_enthalpy,f%vapor_cp]))) return
    if (f%density <= 0.0_dp .or. f%density > 1.0e8_dp) return
    if (f%viscosity < 1.0e-15_dp .or. f%viscosity > 1.0e5_dp) return
    if (f%conductivity < 1.0e-15_dp .or. f%conductivity > 1.0e8_dp) return
    if (f%diffusivity < 1.0e-20_dp .or. f%diffusivity > 1.0e8_dp) return
    if (f%cp <= 0.0_dp .or. f%cp > 1.0e10_dp) return
    if (f%vapor_cp<0.or.f%vapor_cp>1e10_dp) return
    if (f%gas_temperature < 1.0_dp .or. f%gas_temperature > 1.0e5_dp) return
    if (maxval(abs(f%gas_velocity)) > 1.0e7_dp .or. abs(f%vapor_enthalpy) > 1.0e12_dp) return
    if (f%surface_fraction < 0.0_dp .or. f%surface_fraction >= 1.0_dp) return
    if (f%ambient_fraction < 0.0_dp .or. f%ambient_fraction > 1.0_dp) return
    if (p%mass == 0.0_dp) then
      ok=.true.; return
    end if
    if (liquid_heat_capacity_at_temperature(l,p%temperature)<=0) return
    rho=liquid_density_at_temperature(l,p%temperature)
    if (rho<=0) return
    d=parcel_diameter(p,l)
    if (d < 1.0e-12_dp) return
    re=f%density*d*norm2(f%gas_velocity-p%velocity)/f%viscosity
    cdre=24.0_dp*(1.0_dp+0.15_dp*re**0.687_dp)
    if (re >= 1000.0_dp) cdre=0.44_dp*re
    if (o%drag_model=='pelephysics') then
      ! Pinned Drag.H: store Cd*Re to avoid division by a vanishing Reynolds number.
      cdre=0
      if (re>1000) then
        cdre=.424_dp*re
      else if (re>1) then
        cdre=24*(1+re**(2.0_dp/3)/6)
      else if (re>0) then
        cdre=24
      end if
      if (o%breakup_model=='etab') then
        if (1+2.6432_dp*p%deformation<0) return
        cdre=cdre*(1+2.6432_dp*p%deformation)
      end if
    end if
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
    if (o%abramzon_sirignano) bm=min(20.0_dp,bm)
    heat_cp=f%cp
    if (f%vapor_cp>0) heat_cp=f%vapor_cp
    logbm=log(1.0_dp+bm)
    sh=sh0; nu=nu0; factor=1.0_dp
    if (o%abramzon_sirignano .and. bm > 1.0e-12_dp) then
      sh=2.0_dp+(sh0-2.0_dp)/blowing_factor(bm)
      phi=heat_cp*f%density*f%diffusivity*sh/(f%conductivity*nu0)
      if (phi*logbm>100.0_dp) return
      bt=(1.0_dp+bm)**phi-1.0_dp
      ! Pinned Drag.H calcHeatCoeff: initial correction then at most 100
      ! direct B_T updates, with absolute B_T convergence tolerance 1e-7.
      do iteration=0,100
        old_bt=bt
        nu=2.0_dp+(nu0-2.0_dp)/blowing_factor(bt)
        phi=heat_cp*f%density*f%diffusivity*sh/(f%conductivity*nu)
        if (phi*logbm > 100.0_dp) return
        bt=(1.0_dp+bm)**phi-1.0_dp
        if (.not.ieee_is_finite(bt)) return
        if (abs(bt-old_bt)<=1.0e-7_dp) exit
      end do
      if (iteration > 100) return
      nu=2.0_dp+(nu0-2.0_dp)/blowing_factor(bt)
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
    real(dp) :: invtau, heat, masscoef, d, d2, dm, mbar, decay, latent, rate, equilibrium,liquid_cp
    q=p; transfer=0.0_dp; ok=.false.
    if (.not. ieee_is_finite(dt)) return
    if (dt < 0.0_dp .or. dt > 1.0e6_dp) return
    call spray_transfer_rates(p,l,f,o,invtau,heat,masscoef,ok)
    if (.not. ok) return
    if (dt == 0.0_dp .or. p%mass == 0.0_dp) return
    ok=.false.
    if (liquid_heat_capacity_at_temperature(l,p%temperature)<=0) return
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
      liquid_cp=liquid_heat_capacity_at_temperature(l,p%temperature)
      rate=dm/dt
      if (heat > 0.0_dp) then
        equilibrium=f%gas_temperature-rate*latent/heat
        decay=exp(-min(700.0_dp,heat*dt/(mbar*liquid_cp)))
        q%temperature=equilibrium+(p%temperature-equilibrium)*decay
      else
        q%temperature=p%temperature-dm*latent/(mbar*liquid_cp)
      end if
    end if
    if (.not. valid_spray_parcel(q)) then
      q=p; return
    end if
    if (q%mass>0) then
      if (liquid_density_at_temperature(l,q%temperature)<=0.or. &
          liquid_heat_capacity_at_temperature(l,q%temperature)<=0) then
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
      half_angle,diameter,temperature,parcels,ok,diameter_shape,diameter_min,diameter_max, &
      port_speeds,port_diameters,port_temperatures,port_angles,materials,diameter_weighting, &
      diameter_distribution,diameter_std,port_size_options,nozzle_radius,port_nozzle_radii, &
      nozzle_model,port_nozzle_models)
    type(spray_liquid), intent(in) :: liquid
    type(spray_nozzle_model), intent(in), optional :: nozzle_model,port_nozzle_models(:)
    type(spray_nozzle_model) :: local_nozzle
    real(dp), intent(in), optional :: nozzle_radius,port_nozzle_radii(:)
    real(dp) :: local_radius
    type(spray_size_options), intent(in), optional :: port_size_options(:)
    type(spray_liquid), intent(in), optional :: materials(:)
    type(spray_liquid) :: material
    real(dp), intent(in) :: total_mass,centers(:,:),axes(:,:),fractions(:),speed,half_angle,diameter,temperature
    integer, intent(in) :: number
    integer(int64), intent(in) :: first_id
    type(spray_parcel), allocatable, intent(out) :: parcels(:)
    logical, intent(out) :: ok
    character(len=*), intent(in), optional :: diameter_weighting,diameter_distribution
    real(dp), intent(in), optional :: diameter_std
    real(dp), intent(in), optional :: diameter_shape,diameter_min,diameter_max
    real(dp), intent(in), optional :: port_speeds(:),port_diameters(:),port_temperatures(:),port_angles(:)
    real(dp) :: local_speed,local_diameter,local_temperature,local_angle
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
    if (present(port_speeds)) then
      if (size(port_speeds)/=ports) return
    end if
    if (present(port_diameters)) then
      if (size(port_diameters)/=ports) return
    end if
    if (present(port_temperatures)) then
      if (size(port_temperatures)/=ports) return
    end if
    if (present(port_angles)) then
      if (size(port_angles)/=ports) return
    end if
    if (present(materials)) then
      if (size(materials)/=ports) return
    end if
    if (present(port_size_options)) then
      if (size(port_size_options)/=ports) return
    end if
    if (present(port_nozzle_radii)) then
      if (size(port_nozzle_radii)/=ports) return
    end if
    if (present(port_nozzle_models)) then
      if (size(port_nozzle_models)/=ports) return
    end if
    total=ports*number
    if (first_id<1_int64.or.first_id>huge(first_id)-int(total,int64)) return
    allocate(candidate(total))
    do port=1,ports
      offset=(port-1)*number
      local_speed=speed;local_diameter=diameter;local_temperature=temperature;local_angle=half_angle
      if (present(port_speeds)) local_speed=port_speeds(port)
      if (present(port_diameters)) local_diameter=port_diameters(port)
      if (present(port_temperatures)) local_temperature=port_temperatures(port)
      if (present(port_angles)) local_angle=port_angles(port)
      local_nozzle=spray_nozzle_model()
      if (present(nozzle_model)) local_nozzle=nozzle_model
      if (present(port_nozzle_models)) local_nozzle=port_nozzle_models(port)
      local_radius=0
      if (present(nozzle_radius)) local_radius=nozzle_radius
      if (present(port_nozzle_radii)) local_radius=port_nozzle_radii(port)
      material=liquid
      if (present(materials)) material=materials(port)
      if (present(port_size_options)) then
        call inject_cone_parcels(material,total_mass*(fractions(port)/fraction_sum),number,first_id+int(offset,int64), &
          centers(:,port),axes(:,port),local_speed,local_angle,local_diameter,local_temperature,port_parcels,ok, &
          port_size_options(port)%shape,port_size_options(port)%minimum,port_size_options(port)%maximum, &
          port_size_options(port)%weighting,port_size_options(port)%distribution,port_size_options(port)%std, &
          local_radius,local_nozzle)
      else
      call inject_cone_parcels(material,total_mass*(fractions(port)/fraction_sum),number,first_id+int(offset,int64), &
        centers(:,port),axes(:,port),local_speed,local_angle,local_diameter,local_temperature,port_parcels,ok, &
        diameter_shape,diameter_min,diameter_max,diameter_weighting,diameter_distribution,diameter_std,local_radius,local_nozzle)
      end if
      if (.not.ok) return
      if (present(materials)) port_parcels%material_id=port
      candidate(offset+1:offset+number)=port_parcels
      deallocate(port_parcels)
    end do
    call move_alloc(candidate,parcels)
  end subroutine inject_cone_ports

  subroutine inject_cone_parcels(liquid, total_mass, number, first_id, center, axis, speed, &
      half_angle, diameter, temperature, parcels, ok, diameter_shape, diameter_min, diameter_max, &
      diameter_weighting,diameter_distribution,diameter_std,nozzle_radius,nozzle_model)
    type(spray_liquid), intent(in) :: liquid
    real(dp), intent(in) :: total_mass, center(3), axis(3), speed, half_angle, diameter, temperature
    integer, intent(in) :: number
    integer(int64), intent(in) :: first_id
    type(spray_parcel), allocatable, intent(out) :: parcels(:)
    type(spray_nozzle_model), intent(in), optional :: nozzle_model
    type(spray_nozzle_model) :: model
    real(dp) :: direction_angle,spread_angle
    real(dp), intent(in), optional :: nozzle_radius
    real(dp) :: radius,radial_position,position(3)
    logical, intent(out) :: ok
    character(len=*), intent(in), optional :: diameter_weighting,diameter_distribution
    real(dp), intent(in), optional :: diameter_std
    character(len=8) :: weighting
    character(len=16) :: distribution
    real(dp) :: std,location,width,zlo,zhi,plo,phi,target,left,right,mid,zmid,prob
    real(dp) :: chi_cdf(100),chi_sizes(100),xi,previous_cdf,current_cdf,probability_sum
    integer :: bin
    logical :: use_survival
    integer :: iteration
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
    if (liquid_heat_capacity_at_temperature(liquid,temperature)<=0) return
    rho=liquid_density_at_temperature(liquid,temperature)
    if (rho<=0) return
    radius=0
    if (present(nozzle_radius)) radius=nozzle_radius
    if (.not.ieee_is_finite(radius)) return
    if (radius<0.or.radius>1e9_dp) return
    model=spray_nozzle_model()
    if (present(nozzle_model)) model=nozzle_model
    if (model%profile/='legacy'.and.model%profile/='solid'.and.model%profile/='hollow') return
    if (.not.all(ieee_is_finite([model%swirl,model%spread]))) return
    if (abs(model%swirl)>pi.or.model%spread<0.or.model%spread>pi) return
    if (model%profile/='hollow'.and.model%spread/=0) return
    if (model%profile/='legacy'.and.radius<=0) return
    weighting='mass'
    if (present(diameter_weighting)) then
      if (diameter_weighting/='mass'.and.diameter_weighting/='number') return
      weighting=diameter_weighting
    end if
    shape=0.0_dp; dlo=1.0e-12_dp; dhi=1.0_dp
    if (present(diameter_shape)) shape=diameter_shape
    if (present(diameter_min)) dlo=diameter_min
    if (present(diameter_max)) dhi=diameter_max
    if (.not. all(ieee_is_finite([shape,dlo,dhi]))) return
    if (shape<0.0_dp .or. shape>20.0_dp) return
    if (dlo<1.0e-12_dp .or. dhi>1.0_dp .or. dlo>=dhi) return
    distribution='weibull';std=0
    if (present(diameter_distribution)) distribution=diameter_distribution
    if (present(diameter_distribution)) then
      if (diameter_distribution/='weibull'.and.diameter_distribution/='normal'.and. &
          diameter_distribution/='lognormal'.and.diameter_distribution/='chi-squared') return
    end if
    if (present(diameter_std)) std=diameter_std
    if (.not.ieee_is_finite(std)) return
    if (std<0.or.std>1) return
    if (std>0.and.std<1e-15_dp) return
    if ((distribution=='weibull'.or.distribution=='chi-squared').and.std/=0) return
    if (distribution/='weibull'.and.shape/=0) return
    if (distribution=='normal'.or.distribution=='lognormal') then
      if (std==0) then
        if (diameter<dlo.or.diameter>dhi) return
      else
        width=std;location=diameter;zlo=dlo;zhi=dhi
        if (distribution=='lognormal') then
          width=(std/diameter)**2
          if (width<1e-8_dp) then
            width=sqrt(width*(1-.5_dp*width))
          else
            width=sqrt(log(1+width))
          end if
          if (width==0) return
          location=log(diameter)-.5_dp*width**2;zlo=log(dlo);zhi=log(dhi)
        end if
        zlo=(zlo-location)/width;zhi=(zhi-location)/width
        use_survival=zlo>=0
        if (use_survival) then
          plo=.5_dp*erfc(zlo/sqrt(2.0_dp));phi=.5_dp*erfc(zhi/sqrt(2.0_dp))
        else
          plo=.5_dp*erfc(-zlo/sqrt(2.0_dp));phi=.5_dp*erfc(-zhi/sqrt(2.0_dp))
        end if
        if (abs(phi-plo)<=max(tiny(1.0_dp),32*epsilon(1.0_dp)*max(plo,phi))) return
      end if
    end if
    if (distribution=='chi-squared') then
      ! Pinned ChiSquared uses 100 left-edge bins over xi=[0,12], d=xi*d32/3.
      ! Condition the discrete probabilities on the configured positive diameter bounds.
      previous_cdf=0;probability_sum=0;xi=0
      do bin=1,100
        chi_sizes(bin)=xi*diameter/3
        xi=.12_dp*real(bin,dp)
        current_cdf=1-exp(-xi)*(1+xi*(1+xi*(.5_dp+xi/6)))
        if (chi_sizes(bin)>=dlo.and.chi_sizes(bin)<=dhi) &
          probability_sum=probability_sum+current_cdf-previous_cdf
        chi_cdf(bin)=probability_sum;previous_cdf=current_cdf
      end do
      if (probability_sum<=32*epsilon(1.0_dp)) return
      chi_cdf=chi_cdf/probability_sum
    end if
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
        ! Midpoint quadrature of the truncated CDF; weights below select mass or number basis.
        quantile=(real(i,dp)-0.5_dp)/real(number,dp)
        survival=survival_low-quantile*(survival_low-survival_high)
        drop_size=diameter*(-log(survival))**(1.0_dp/shape)
        if (drop_size<dlo .or. drop_size>dhi) then
          deallocate(parcels); return
        end if
      end if
      if (distribution/='weibull'.and.std>0) then
        target=plo+((real(i,dp)-.5_dp)/real(number,dp))*(phi-plo)
        left=dlo;right=dhi
        do iteration=1,100
          mid=left+.5_dp*(right-left)
          if (mid==left.or.mid==right) exit
          zmid=mid
          if (distribution=='lognormal') zmid=log(mid)
          zmid=(zmid-location)/width
          if (use_survival) then
            prob=.5_dp*erfc(zmid/sqrt(2.0_dp))
          else
            prob=.5_dp*erfc(-zmid/sqrt(2.0_dp))
          end if
          if ((use_survival.and.prob>target).or.(.not.use_survival.and.prob<target)) then
            left=mid
          else
            right=mid
          end if
        end do
        drop_size=left+.5_dp*(right-left)
      end if
      if (distribution=='chi-squared') then
        target=(real(i,dp)-.5_dp)/real(number,dp)
        bin=1
        do while(bin<100)
          if (target<=chi_cdf(bin)) exit
          bin=bin+1
        end do
        drop_size=chi_sizes(bin)
      end if
      m=pi*rho*drop_size**3/6.0_dp
      weight=1
      if (weighting=='mass') weight=total_mass/(real(number,dp)*m)
      if (weight>1.0e18_dp .or. m>1.0_dp) then
        deallocate(parcels); return
      end if
      ! Deterministic equal-solid-angle quadrature, not a stochastic turbulence model.
      cosine=1.0_dp-(1.0_dp-cos(half_angle))*(real(i,dp)-0.5_dp)/real(number,dp)
      sine=sqrt(max(0.0_dp,1.0_dp-cosine*cosine))
      if (model%profile/='legacy') then
        spread_angle=half_angle*(real(i,dp)-.5_dp)/real(number,dp)
        if (model%profile=='hollow') spread_angle=half_angle+ &
          model%spread*((real(i,dp)-.5_dp)/real(number,dp)-.5_dp)
        cosine=cos(spread_angle);sine=sin(spread_angle)
      end if
      angle=2.0_dp*pi*modulo(real(first_id-1_int64+int(i,int64),dp)*0.6180339887498949_dp,1.0_dp)
      position=center
      if (radius>0) then
        ! SprayInjection.H samples radius linearly, not uniformly in disk area.
        radial_position=radius*(real(i,dp)-.5_dp)/real(number,dp)
        if (model%profile=="hollow") radial_position=radius
        position=center+radial_position*(cos(angle)*e1+sin(angle)*e2)
      end if
      direction_angle=angle+model%swirl
      parcels(i)=spray_parcel(first_id+int(i-1,int64),position, &
        speed*(cosine*ez+sine*(cos(direction_angle)*e1+sin(direction_angle)*e2)),m,temperature,weight)
      if (.not. valid_spray_parcel(parcels(i))) then
        deallocate(parcels); return
      end if
    end do
    if (weighting=='number') then
      ! Equal represented droplet counts; normalize their unequal masses to the requested dose.
      weight=total_mass/sum(parcels%mass)
      if (weight<=0.or.weight>1e18_dp) then
        deallocate(parcels);return
      end if
      parcels%multiplicity=weight
    end if
    ok=.true.
  end subroutine
end module spray_parcel_mod
