! Ideal liquid mixing rules adapted from PelePhysics SprayPropertiesPeleMP.H,
! revision 6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b.
! Copyright/redistribution: references/pelephysics-spray-license.txt.
module spray_mixture_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use gas_transport_mod, only: constant_transport_values
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species,nasa7_mass_properties
  use gas_transport_mod, only: gas_transport_species
  use mixture_thermo_mod, only: mixture_mass_properties
  use mixture_transport_mod, only: mixture_transport_coefficients,binary_diffusion_coefficients
  use spray_parcel_mod
  implicit none
  private
  real(dp), parameter, public :: mixture_trace_fraction=1e-12_dp
  character(len=*), parameter, public :: mixture_evolution_id='mixed-liquid-components-v32-film-vapor-gate'
  public :: liquid_mixture_properties, liquid_mixture_temperature, advance_mixture_heating
  public :: liquid_mixture_conductivity,liquid_mixture_viscosity,apply_mixture_breakup_viscosity
  public :: liquid_mixture_critical_limit
  public :: film_fuel_diffusion_coefficients
  public :: liquid_mixture_surface, advance_mixture_evaporation
  public :: build_mixture_film, mixture_transfer_rates, advance_mixture_frozen_gas, valid_component_inventory
contains
  pure subroutine euler_component_remainder(masses,evaporated,remaining,total)
    real(dp), intent(in) :: masses(:),evaporated(:)
    real(dp), intent(out) :: remaining(:),total
    real(dp) :: survivor_mass
    total=max(0._dp,sum(masses)-sum(evaporated))
    remaining=0
    if(total<=0) return
    remaining=min(total,max(0._dp,masses-evaporated))
    if(size(masses)>1) then
      ! Drag.H clips each fraction to [0,1], then removes Y<=1e-12.
      where(remaining/total<=mixture_trace_fraction) remaining=0
    endif
    survivor_mass=sum(remaining)
    if(survivor_mass>0) then
      remaining=remaining*(total/survivor_mass)
    else
      total=0
    endif
  end subroutine

  pure subroutine liquid_mixture_critical_limit(components,y,temperature,limited,ok)
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(in) :: y(:),temperature
    real(dp), intent(out) :: limited
    logical, intent(out) :: ok
    real(dp) :: critical
    limited=temperature;ok=.false.
    if(.not.valid_composition(components,y)) return
    if(.not.ieee_is_finite(temperature)) return
    if(temperature<1.or.temperature>1e5_dp) return
    ! Legacy liquids without critical-temperature inputs retain their model.
    ! All active PeleMP components provide positive critical temperatures.
    if(all(components%critical_temperature>0.or.y==0)) then
      critical=sum(y*components%critical_temperature)
      limited=min(temperature,.999_dp*critical)
    endif
    ok=limited>=1
  end subroutine

  pure logical function valid_component_inventory(parcel,masses) result(ok)
    type(spray_parcel), intent(in) :: parcel
    real(dp), intent(in) :: masses(:)
    ok=.false.
    if (.not.valid_spray_parcel(parcel).or.size(masses)<1.or.size(masses)>16) return
    if (.not.all(ieee_is_finite(masses))) return
    if (any(masses<0).or.any(masses>1)) return
    ok=abs(sum(masses)-parcel%mass)<=32*epsilon(1.0_dp)*size(masses)*parcel%mass
  end function

  pure logical function valid_composition(components,y) result(ok)
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(in) :: y(:)
    integer :: i,j
    ok=.false.
    if (size(y)<1.or.size(y)>16.or.size(components)/=size(y)) return
    if (.not.all(ieee_is_finite(y))) return
    if (any(y<0).or.any(y>1)) return
    if (abs(sum(y)-1)>32*epsilon(1.0_dp)*size(y)) return
    do i=1,size(y)
      if (.not.valid_spray_liquid(components(i))) return
      do j=1,i-1
        if (components(i)%vapor_index==components(j)%vapor_index.and. &
            components(i)%vapor_name/=components(j)%vapor_name) return
      end do
    end do
    ok=.true.
  end function

  pure subroutine apply_mixture_breakup_viscosity(components,y,temperature,settings,ok,pressure,gas_molecular_weights)
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(in) :: y(:),temperature
    type(spray_options), intent(inout) :: settings
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: pressure,gas_molecular_weights(:)
    type(spray_options) :: candidate
    real(dp) :: value
    ok=.false.
    if(.not.valid_composition(components,y)) return
    if(.not.any(components%viscosity_model/='supplied'.and.y>0)) then
      ok=valid_spray_options(settings);return
    endif
    value=0
    call liquid_mixture_viscosity(components,y,temperature,0._dp,value,ok,pressure,gas_molecular_weights)
    if(.not.ok) return
    candidate=settings;candidate%liquid_viscosity_model='constant';candidate%liquid_viscosity=value
    candidate%liquid_viscosity_coefficients=0;candidate%liquid_viscosity_temperature_range=0
    ok=valid_spray_options(candidate)
    if(ok) settings=candidate
  end subroutine

  pure subroutine liquid_mixture_viscosity(components,y,temperature,supplied,value,ok,pressure,gas_molecular_weights)
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(in) :: y(:),temperature,supplied
    real(dp), intent(inout) :: value
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: pressure,gas_molecular_weights(:)
    real(dp) :: candidate,coefficient,mw
    logical :: configured
    integer :: i
    ok=.false.
    if(.not.valid_composition(components,y)) return
    if(.not.ieee_is_finite(temperature)) return
    if(temperature<1.or.temperature>1e5_dp) return
    if(present(pressure)) then
      if(.not.ieee_is_finite(pressure)) return
      if(pressure<1e-6_dp.or.pressure>1e12_dp) return
    endif
    configured=any(components%viscosity_model/='supplied'.and.y>0)
    candidate=0
    if(configured) then
      if(present(gas_molecular_weights)) then
        if(size(gas_molecular_weights)/=size(y)) return
      endif
      do i=1,size(y)
        if(y(i)==0) cycle
        if(components(i)%viscosity_model=='supplied') return
        mw=components(i)%molecular_weight
        if(present(gas_molecular_weights)) mw=gas_molecular_weights(i)
        coefficient=liquid_component_viscosity_at_temperature(components(i),temperature,pressure,mw)
        if(coefficient<=0) return
        candidate=candidate+y(i)*coefficient
      enddo
    else
      candidate=supplied
    endif
    if(.not.ieee_is_finite(candidate)) return
    if(candidate<1e-10_dp.or.candidate>1e3_dp) return
    value=candidate;ok=.true.
  end subroutine

  pure subroutine liquid_mixture_conductivity(components,y,temperature,supplied,value,ok,pressure,gas_molecular_weights)
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(in) :: y(:),temperature,supplied
    real(dp), intent(inout) :: value
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: pressure,gas_molecular_weights(:)
    real(dp) :: candidate,coefficient,mw
    logical :: configured
    integer :: i
    ok=.false.
    if(.not.valid_composition(components,y)) return
    if(.not.ieee_is_finite(temperature)) return
    if(temperature<1.or.temperature>1e5_dp) return
    if(present(pressure)) then
      if(.not.ieee_is_finite(pressure)) return
      if(pressure<1e-6_dp.or.pressure>1e12_dp) return
    endif
    configured=any(components%conductivity_model/='supplied'.and.y>0)
    candidate=0
    if(configured) then
      if(present(gas_molecular_weights)) then
        if(size(gas_molecular_weights)/=size(y)) return
      endif
      do i=1,size(y)
        if(y(i)==0) cycle
        if(components(i)%conductivity_model=='supplied') return
        mw=components(i)%molecular_weight
        if(present(gas_molecular_weights)) mw=gas_molecular_weights(i)
        coefficient=liquid_conductivity_at_temperature(components(i),temperature,pressure,mw)
        if(coefficient<=0) return
        candidate=candidate+y(i)*coefficient
      enddo
    else
      candidate=supplied
    endif
    if(.not.ieee_is_finite(candidate)) return
    if(candidate<=0.or.candidate>1e4_dp) return
    value=candidate;ok=.true.
  end subroutine

  pure subroutine liquid_mixture_properties(components,y,temperature,density,cp,energy,ok,pressure,gas_molecular_weights)
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(in) :: y(:),temperature
    real(dp), intent(inout) :: density,cp,energy
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: pressure,gas_molecular_weights(:)
    real(dp) :: volume,candidate_cp,candidate_e,rhoi,cpi,ei,mw
    logical :: valid
    integer :: i
    ok=.false.
    if (.not.valid_composition(components,y).or..not.ieee_is_finite(temperature)) return
    if (temperature<1.or.temperature>1e5_dp) return
    if (present(pressure)) then
      if (.not.ieee_is_finite(pressure)) return
      if (pressure<1e-6_dp.or.pressure>1e12_dp) return
    end if
    if (present(gas_molecular_weights)) then
      if (size(gas_molecular_weights)/=size(y)) return
      if (.not.all(ieee_is_finite(gas_molecular_weights))) return
      if (any(gas_molecular_weights<1e-3_dp).or.any(gas_molecular_weights>1e5_dp)) return
    end if
    volume=0;candidate_cp=0;candidate_e=0
    do i=1,size(y)
      if (y(i)==0) cycle
      if (present(pressure)) then
        mw=components(i)%molecular_weight
        if (present(gas_molecular_weights)) mw=gas_molecular_weights(i)
        call liquid_transfer_properties(components(i),temperature,pressure,mw,rhoi,cpi,valid)
        if (.not.valid) return
      else
        rhoi=liquid_density_at_temperature(components(i),temperature)
        cpi=liquid_heat_capacity_at_temperature(components(i),temperature)
      end if
      if (rhoi<=0.or.cpi<=0) return
      ei=liquid_energy(components(i),temperature)
      if (.not.ieee_is_finite(ei)) return
      volume=volume+y(i)/rhoi
      candidate_cp=candidate_cp+y(i)*cpi
      candidate_e=candidate_e+y(i)*ei
    end do
    if (.not.all(ieee_is_finite([volume,candidate_cp,candidate_e])).or.volume<=0) return
    density=1/volume;cp=candidate_cp;energy=candidate_e;ok=.true.
  end subroutine

  pure subroutine liquid_mixture_temperature(components,y,energy,temperature,ok)
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(in) :: y(:),energy
    real(dp), intent(inout) :: temperature
    logical, intent(out) :: ok
    real(dp) :: lower,upper,elo,ehi,t,rho,cp,e,step
    integer :: i,iteration
    logical :: valid
    ok=.false.
    if (.not.valid_composition(components,y)) return
    if (.not.all(ieee_is_finite([temperature,energy]))) return
    lower=1;upper=1e5_dp
    do i=1,size(y)
      if (y(i)==0) cycle
      if (components(i)%density_model=='polynomial') then
        lower=max(lower,components(i)%density_temperature_range(1))
        upper=min(upper,components(i)%density_temperature_range(2))
      end if
      if (components(i)%heat_capacity_model=='polynomial') then
        lower=max(lower,components(i)%heat_capacity_temperature_range(1))
        upper=min(upper,components(i)%heat_capacity_temperature_range(2))
      end if
    end do
    if (lower>upper) return
    call liquid_mixture_properties(components,y,lower,rho,cp,elo,valid)
    if (.not.valid) return
    call liquid_mixture_properties(components,y,upper,rho,cp,ehi,valid)
    if (.not.valid.or.energy<elo.or.energy>ehi) return
    if (energy==elo) then
      temperature=lower;ok=.true.;return
    else if (energy==ehi) then
      temperature=upper;ok=.true.;return
    end if
    t=max(lower,min(upper,temperature))
    do iteration=1,100
      call liquid_mixture_properties(components,y,t,rho,cp,e,valid)
      if (.not.valid) return
      if (abs(e-energy)<=8*epsilon(1.0_dp)*max(1.0_dp,abs(e),abs(energy))) then
        temperature=t;ok=.true.;return
      end if
      if (e<energy) then
        lower=t
      else
        upper=t
      end if
      step=t+(energy-e)/cp
      if (step<=lower.or.step>=upper) step=lower+.5_dp*(upper-lower)
      if (step==t) return
      t=step
    end do
  end subroutine

  subroutine film_fuel_diffusion_coefficients(species,transport,y,temperature,pressure,mixture_diffusion,diffusion,ok)
    ! Simple.H returns rho*D with the species molecular weight in its
    ! numerator. Drag.H converts it by Wskin/Wfuel before vapor weighting.
    ! The Chapman-Enskog path retains its carrier-mole harmonic convention.
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: y(:),temperature,pressure,mixture_diffusion(:)
    real(dp), intent(out) :: diffusion(:)
    logical, intent(out) :: ok
    real(dp) :: binary(size(species),size(species)),candidate(size(species)),mw,carrier,resistance
    real(dp) :: constant_values(4)
    logical :: constant_enabled,valid
    integer :: i,j,n
    diffusion=0;ok=.false.;n=size(species)
    if(n<2.or.size(transport)/=n.or.size(y)/=n.or.size(mixture_diffusion)/=n.or.size(diffusion)/=n) return
    if(.not.all(ieee_is_finite([y,mixture_diffusion,species%molecular_weight,temperature,pressure]))) return
    if(any(y<0).or.abs(sum(y)-1)>5e-12_dp.or.any(mixture_diffusion<0)) return
    if(any(species%molecular_weight<1e-3_dp).or.any(species%molecular_weight>1e5_dp)) return
    if(temperature<1.or.temperature>1e5_dp.or.pressure<1e-6_dp.or.pressure>1e12_dp) return
    mw=1/sum(y/species%molecular_weight)
    call constant_transport_values(transport,constant_values,constant_enabled,valid)
    if(.not.valid) return
    if(constant_enabled.and.constant_values(4)>=0) then
      candidate=mixture_diffusion*mw/species%molecular_weight
      if(.not.all(ieee_is_finite(candidate))) return
      diffusion=candidate;ok=.true.;return
    endif
    if(any(mixture_diffusion<=0)) return
    candidate=mixture_diffusion
    if(allocated(transport(1)%simple_binary)) then
      candidate=candidate*mw/species%molecular_weight
    else
      call binary_diffusion_coefficients(species,transport,temperature,pressure,binary,ok)
      if(.not.ok) return
      ok=.false.
      do i=1,n
        carrier=0;resistance=0
        do j=1,n
          if(i==j) cycle
          carrier=carrier+y(j)/species(j)%molecular_weight
          resistance=resistance+y(j)/(species(j)%molecular_weight*binary(i,j))
        enddo
        if(carrier>0.and.resistance>0) candidate(i)=carrier/resistance
      enddo
    endif
    if(.not.all(ieee_is_finite(candidate)).or.any(candidate<=0)) return
    diffusion=candidate;ok=.true.
  end subroutine

  subroutine build_mixture_film(species,transport,primitive,gas_temperature,parcel,components,y,options, &
      film,mass_weights,vapor_h,ok,wall_skin)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_parcel), intent(in) :: parcel
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    real(dp), intent(in) :: primitive(:),gas_temperature,y(:)
    type(spray_film), intent(out) :: film
    real(dp), intent(out) :: mass_weights(:),vapor_h(:)
    logical, optional, intent(in) :: wall_skin
    logical, intent(out) :: ok
    real(dp) :: xs(size(y)),mw_liq(size(y)),diff_weights(size(y)),enthalpy(size(y))
    real(dp) :: yf(size(species)),yg(size(species)),xg(size(species)),ys(size(species)),mw_g(size(species))
    real(dp) :: diffusion(size(species)),fuel_diffusion(size(species)),tf,r,mw,cp,cv,gamma,h,e,entropy
    real(dp) :: ambient_mw,carrier_y,denom,sum_x,cp_i
    real(dp) :: liquid_x(size(y)),group_x(size(species)),ambient_component
    logical :: active(size(species)),fuel_species(size(species)),valid
    type(spray_film) :: candidate
    integer :: i,j,v
    ok=.false.;film=spray_film();mass_weights=0;vapor_h=0
    if (.not.valid_composition(components,y).or..not.valid_spray_parcel(parcel).or..not.valid_spray_options(options)) return
    if (size(primitive)/=5+size(species).or.size(species)<2) return
    if (size(mass_weights)/=size(y).or.size(vapor_h)/=size(y)) return
    if (.not.all(ieee_is_finite(primitive)).or..not.ieee_is_finite(gas_temperature)) return
    if (gas_temperature<1.or.gas_temperature>1e5_dp) return
    yg=primitive(6:);mw_g=species%molecular_weight
    if (any(yg<0).or.abs(sum(yg)-1)>5e-12_dp) return
    if (.not.all(ieee_is_finite(mw_g))) return
    if (any(mw_g<1e-3_dp).or.any(mw_g>1e5_dp)) return
    xg=(yg/mw_g)/sum(yg/mw_g);fuel_species=.false.
    do i=1,size(y)
      v=components(i)%vapor_index
      if (v>size(species)) return
      if (trim(components(i)%vapor_name)/=trim(species(v)%name)) return
      fuel_species(v)=.true.
      mw_liq(i)=resolved_liquid_molecular_weight(components(i),mw_g(v))
      call nasa7_mass_properties(species(v),parcel%temperature,cp_i,cv,h,e,entropy,valid)
      if (.not.valid) return
      enthalpy(i)=h
    end do
    xs=0;active=.false.
    if (options%evaporation) then
      call liquid_mixture_surface(components,y,mw_liq,parcel%temperature,primitive(5),xs,valid,enthalpy)
      if (.not.valid) return
      liquid_x=(y/mw_liq)/sum(y/mw_liq);group_x=0
      do i=1,size(y)
        v=components(i)%vapor_index
        group_x(v)=group_x(v)+liquid_x(i)
      end do
      do i=1,size(y)
        v=components(i)%vapor_index
        ! Pinned Drag.H: distribute ambient mole fraction to liquid species,
        ! compare pure-component psat to its share, then retain Raoult vapor.
        ambient_component=liquid_x(i)*xg(v)
        if (group_x(v)>liquid_x(i)) ambient_component=ambient_component/group_x(v)
        if (xs(i)<=liquid_x(i)*ambient_component) xs(i)=0
        active(v)=active(v).or.xs(i)>0
      end do
    end if
    sum_x=sum(xs);yf=yg;ys=0
    if (sum_x>0) then
      carrier_y=sum(yg,mask=.not.active);denom=sum(yg/mw_g,mask=.not.active)
      if (carrier_y<=0.or.denom<=0) return
      ambient_mw=1/sum(yg/mw_g)
      denom=sum(xs*mw_liq)+(1-sum_x)*ambient_mw
      do i=1,size(y)
        v=components(i)%vapor_index
        ys(v)=ys(v)+xs(i)*mw_g(v)/denom
      end do
      if (sum(ys)>=1) return
      where (active) yf=(2*ys+yg)/3
      where (.not.active) yf=yg*(1-sum(yf,mask=active))/carrier_y
    end if
    tf=(2*parcel%temperature+gas_temperature)/3
    if(options%abramzon_sirignano.and..not.options%evaporation) tf=gas_temperature
    ! Drag.H uses a nonnegative gas-to-drop temperature difference for the
    ! one-third-rule skin when the reference transfer model is selected.
    if (options%abramzon_sirignano.and.options%evaporation) tf=parcel%temperature+ &
      max(0.0_dp,gas_temperature-parcel%temperature)/3
    if (present(wall_skin)) then
      if (wall_skin) then
        tf=parcel%temperature+.5_dp*max(0.0_dp,gas_temperature-parcel%temperature)
        if (sum_x>0) then
          do i=1,size(y)
            v=components(i)%vapor_index
            if (active(v)) yf(v)=.5_dp*(ys(v)+yg(v))
          end do
          where (.not.active) yf=yg*(1-sum(yf,mask=active))/carrier_y
        end if
      end if
    end if
    call mixture_mass_properties(species,yf,tf,mw,r,cp,cv,gamma,h,e,entropy,valid)
    if (.not.valid) return
    candidate%density=primitive(5)/(r*tf);candidate%cp=cp
    if (sum_x>0) then
      ! Pinned Drag.H normalizes skin Y over all mapped fuel species, including
      ! non-evaporating components. Loop over gas species to avoid alias duplication.
      denom=sum(yf,mask=fuel_species)
      if (denom<=0) return
      do j=1,size(species)
        if (.not.fuel_species(j)) cycle
        call nasa7_mass_properties(species(j),tf,cp_i,cv,h,e,entropy,valid)
        if (.not.valid) return
        candidate%vapor_cp=candidate%vapor_cp+yf(j)*cp_i
      end do
      candidate%vapor_cp=candidate%vapor_cp/denom
    end if
    call mixture_transport_coefficients(species,transport,yf,tf,primitive(5), &
      candidate%viscosity,candidate%conductivity,diffusion,valid)
    if (.not.valid) return
    call film_fuel_diffusion_coefficients(species,transport,yf,tf,primitive(5),diffusion,fuel_diffusion,valid)
    if (.not.valid) return
    diff_weights=0
    do i=1,size(y)
      if (xs(i)==0) cycle
      v=components(i)%vapor_index
      diff_weights(i)=xs(i)*fuel_diffusion(v)
    end do
    candidate%diffusivity=fuel_diffusion(components(1)%vapor_index)
    if (sum_x>0) candidate%diffusivity=sum(diff_weights)/sum_x
    candidate%surface_fraction=sum(ys);candidate%ambient_fraction=sum(yg,mask=active)
    candidate%gas_temperature=gas_temperature;candidate%gas_velocity=primitive(2:4)
    candidate%pressure=primitive(5);candidate%liquid_molecular_weight=mw_liq(1)
    candidate%vapor_enthalpy=sum(y*enthalpy)
    if (sum(diff_weights)>0) mass_weights=diff_weights/sum(diff_weights)
    vapor_h=enthalpy;film=candidate;ok=.true.
  end subroutine

  subroutine mixture_transfer_rates(species,transport,primitive,gas_temperature,parcel,components,y,options, &
      inverse_tau,heat_coefficient,mass_rate,vapor_h,ok,reynolds_number)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_parcel), intent(in) :: parcel
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    real(dp), intent(in) :: primitive(:),gas_temperature,y(:)
    real(dp), intent(out) :: inverse_tau,heat_coefficient,mass_rate(:),vapor_h(:)
    logical, intent(out) :: ok
    real(dp), optional, intent(inout) :: reynolds_number
    type(spray_film) :: film
    type(spray_liquid) :: effective
    real(dp) :: weights(size(y)),enthalpy(size(y)),rho,cp,e,mass_coefficient,rate,heat,d
    logical :: valid
    inverse_tau=0;heat_coefficient=0;mass_rate=0;vapor_h=0;ok=.false.
    if (size(mass_rate)/=size(y).or.size(vapor_h)/=size(y)) return
    call build_mixture_film(species,transport,primitive,gas_temperature,parcel,components,y,options, &
      film,weights,enthalpy,valid)
    if (.not.valid) return
    call liquid_mixture_properties(components,y,parcel%temperature,rho,cp,e,valid,primitive(5), &
      species(components%vapor_index)%molecular_weight)
    if (.not.valid) return
    ! Only density and cp enter the common transfer correlations; energy stays component-resolved.
    effective=components(1);effective%density_model='constant';effective%density=rho
    effective%heat_capacity_model='constant';effective%heat_capacity=cp
    effective%density_coefficients=0;effective%density_temperature_range=0
    effective%heat_capacity_coefficients=0;effective%heat_capacity_temperature_range=0
    effective%boiling_model='reject';effective%critical_temperature=0;effective%normal_boiling_temperature=0
    call spray_transfer_rates(parcel,effective,film,options,rate,heat,mass_coefficient,valid)
    if (.not.valid) return
    d=parcel_diameter(parcel,effective)
    if(present(reynolds_number)) then
      if(options%fixed_parts) then
        reynolds_number=film%density*norm2(primitive(2:4))*d/film%viscosity
      else
        reynolds_number=film%density*norm2(primitive(2:4)-parcel%velocity)*d/film%viscosity
      endif
    endif
    mass_rate=mass_coefficient*d*weights;inverse_tau=rate;heat_coefficient=heat;vapor_h=enthalpy;ok=.true.
  end subroutine

  subroutine advance_mixture_frozen_gas(species,transport,primitive,gas_temperature,parcel,components,component_mass, &
      options,dt,gas_change,ok,time_modifier,acceleration,force_change,reynolds_number)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_parcel), intent(inout) :: parcel
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    real(dp), intent(in) :: primitive(:),gas_temperature,dt
    real(dp), intent(inout) :: component_mass(:)
    real(dp), intent(out) :: gas_change(:)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: time_modifier
    real(dp), optional, intent(in) :: acceleration(3)
    real(dp), optional, intent(out) :: force_change(5)
    real(dp), optional, intent(inout) :: reynolds_number
    type(spray_parcel) :: candidate
    real(dp) :: masses(size(component_mass)),y(size(component_mass)),rates(size(component_mass)),hv(size(component_mass))
    real(dp) :: exchange(size(gas_change)),change(size(gas_change)),invtau,heat,rho,cp,e,scale,delta,modifier
    integer :: nsub,step
    real(dp) :: accel(3),force(5),step_force(5),step_reynolds
    logical :: valid
    ok=.false.;gas_change=0
    if(present(force_change)) force_change=0
    accel=0
    if(present(acceleration)) accel=acceleration
    if(.not.all(ieee_is_finite(accel))) return
    if(any(abs(accel)>1e6_dp)) return
    if(any(accel/=0).and..not.present(force_change)) return
    if(any(accel/=0).and..not.options%fixed_parts.and.options%mixture_time_scheme/='euler') return
    if(options%fixed_parts.and..not.present(force_change)) return
    if(.not.valid_spray_options(options)) return
    modifier=1
    if(present(time_modifier)) modifier=time_modifier
    if(.not.ieee_is_finite(modifier)) return
    if(modifier<=0.or.modifier>1) return
    if(options%mixture_time_scheme/='euler'.or.dt==0) then
      call advance_mixture_frozen_gas_step(species,transport,primitive,gas_temperature,parcel,components,component_mass, &
        options,modifier*dt,gas_change,ok,accel,step_force,step_reynolds)
      if(ok.and.present(force_change)) force_change=step_force
      if(ok.and.present(reynolds_number)) reynolds_number=step_reynolds
      return
    endif
    if(.not.ieee_is_finite(dt)) return
    if(dt<0.or.dt>1e6_dp.or.size(gas_change)/=5+size(species)) return
    if(.not.valid_component_inventory(parcel,component_mass).or.parcel%mass<=tiny(1._dp)) return
    if(size(components)/=size(component_mass)) return
    y=component_mass/parcel%mass
    call mixture_transfer_rates(species,transport,primitive,gas_temperature,parcel,components,y,options, &
      invtau,heat,rates,hv,valid)
    if(.not.valid) return
    call liquid_mixture_properties(components,y,parcel%temperature,rho,cp,e,valid,primitive(5), &
      species(components%vapor_index)%molecular_weight)
    if(.not.valid) return
    ! Pinned Drag.H estimates all substeps on the first state, before dtmod.
    scale=dt*invtau
    delta=max(0._dp,gas_temperature-parcel%temperature)
    if(options%evaporation.and.delta>1e-15_dp) then
      scale=max(scale,dt*sum(rates)/(.2_dp*parcel%mass),dt*heat/(parcel%mass*cp))
    endif
    if(.not.ieee_is_finite(scale)) return
    nsub=1+int(min(99._dp,max(0._dp,scale)))
    if(nsub>options%maximum_substeps) return
    candidate=parcel;masses=component_mass;exchange=0;force=0
    do step=1,nsub
      if(candidate%mass<=0) exit
      call advance_mixture_frozen_gas_step(species,transport,primitive,gas_temperature,candidate,components,masses, &
        options,modifier*dt/real(nsub,dp),change,valid,accel,step_force,step_reynolds)
      if(.not.valid) return
      ! Drag.H substeps update velocity/thermodynamics, never particle position.
      candidate%position=parcel%position
      exchange=exchange+change
      force=force+step_force
      if(.not.all(ieee_is_finite(exchange)).or..not.all(ieee_is_finite(force))) return
    enddo
    ! SprayParticles.cpp drifts once with the final velocity and the flow
    ! interval, independently of the kick's dtmod.
    if(.not.options%fixed_parts.and.candidate%mass>0) candidate%position=parcel%position+dt*candidate%velocity
    if(.not.valid_spray_parcel(candidate)) return
    parcel=candidate;component_mass=masses;gas_change=exchange;ok=.true.
    if(present(force_change)) force_change=force
    if(present(reynolds_number)) reynolds_number=step_reynolds
  end subroutine

  subroutine advance_mixture_frozen_gas_step(species,transport,primitive,gas_temperature,parcel,components,component_mass, &
      options,dt,gas_change,ok,acceleration,force_change,reynolds_number)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_parcel), intent(inout) :: parcel
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    real(dp), intent(in) :: primitive(:),gas_temperature,dt
    real(dp), intent(inout) :: component_mass(:)
    real(dp), intent(out) :: gas_change(:)
    logical, intent(out) :: ok
    real(dp), intent(in) :: acceleration(3)
    real(dp), intent(out) :: force_change(5)
    real(dp), intent(out) :: reynolds_number
    type(spray_parcel) :: candidate
    real(dp) :: remaining(size(component_mass)),y(size(component_mass)),rates(size(component_mass)),hv(size(component_mass))
    real(dp) :: exchange(size(gas_change)),drag(5),invtau,heat,mass
    real(dp) :: evaporated(size(component_mass)),remaining_mass,rho,cp,e
    real(dp) :: dv_force(3),external(5)
    integer :: i
    logical :: valid
    ok=.false.;gas_change=0;force_change=0
    if (.not.ieee_is_finite(dt).or.dt<0.or.dt>1e6_dp) return
    if (size(gas_change)/=5+size(species).or..not.valid_spray_parcel(parcel)) return
    if (.not.all(ieee_is_finite(component_mass)).or.any(component_mass<0).or.any(component_mass>1)) return
    mass=sum(component_mass)
    if (mass<=tiny(1.0_dp)) return
    if (abs(mass-parcel%mass)>32*epsilon(1.0_dp)*size(component_mass)*mass) return
    y=component_mass/mass
    call mixture_transfer_rates(species,transport,primitive,gas_temperature,parcel,components,y,options, &
      invtau,heat,rates,hv,valid,reynolds_number)
    if (.not.valid) return
    if (dt==0) then
      ok=.true.;return
    end if
    evaporated=rates*dt
    if(options%mixture_time_scheme=='euler') then
      call euler_component_remainder(component_mass,evaporated,remaining,remaining_mass)
    else
      do i=1,size(component_mass)
      if (rates(i)<=0) cycle
      if (component_mass(i)<=mixture_trace_fraction*mass) then
        evaporated(i)=component_mass(i)
      else
        if (evaporated(i)>component_mass(i)) return
        if (component_mass(i)-evaporated(i)<=mixture_trace_fraction*mass) evaporated(i)=component_mass(i)
      end if
      end do
      remaining=component_mass-evaporated;remaining_mass=sum(remaining)
    endif
    if (any(rates>0).and.remaining_mass>0) then
      call liquid_mixture_properties(components,remaining/remaining_mass,parcel%temperature,rho,cp,e,valid,primitive(5), &
        species(components%vapor_index)%molecular_weight)
      if (.not.valid) return
      if (remaining_mass<acos(-1.0_dp)*rho*options%minimum_diameter**3/6) evaporated=component_mass
    end if
    candidate=parcel;remaining=component_mass
    if (options%drag.and..not.options%fixed_parts) &
      candidate%velocity=primitive(2:4)+(parcel%velocity-primitive(2:4))*exp(-min(700.0_dp,invtau*dt))
    if(options%drag.and..not.options%fixed_parts.and.options%mixture_time_scheme=='euler') then
      if(invtau*dt>1) return
      candidate%velocity=parcel%velocity+dt*invtau*(primitive(2:4)-parcel%velocity)
    end if
    dv_force=dt*acceleration
    candidate%velocity=candidate%velocity+dv_force
    external=0
    external(2:4)=mass*parcel%multiplicity*dv_force
    external(5)=mass*parcel%multiplicity*dot_product(.5_dp*(parcel%velocity+candidate%velocity),dv_force)
    candidate%position=parcel%position+.5_dp*dt*(parcel%velocity+candidate%velocity)
    drag=0;drag(2:4)=mass*parcel%multiplicity*(parcel%velocity-candidate%velocity)
    drag(5)=.5_dp*mass*parcel%multiplicity*(sum(parcel%velocity**2)-sum(candidate%velocity**2))
    ! Particle kinetic change includes external work. Only the complementary
    ! drag exchange belongs to the gas; keep the force ledger independent.
    drag=drag+external
    ! With drag disabled, body-force rounding must not create gas exchange.
    if(.not.options%drag) drag=0
    if(options%fixed_parts) then
      candidate%velocity=0;candidate%position=parcel%position
      drag=0;drag(2:4)=-mass*parcel%multiplicity*dt*invtau*primitive(2:4)
      external=drag
      external(2:4)=external(2:4)-mass*parcel%multiplicity*parcel%velocity
      external(5)=-.5_dp*mass*parcel%multiplicity*sum(parcel%velocity**2)
    end if
    call advance_mixture_evaporation(candidate,components,remaining,evaporated,hv, &
      heat*merge(max(0.0_dp,gas_temperature-parcel%temperature), &
        gas_temperature-parcel%temperature,options%abramzon_sirignano.and.options%evaporation)*dt* &
        parcel%multiplicity,exchange,valid, &
      euler_temperature=options%mixture_time_scheme=='euler',pressure=primitive(5), &
      gas_molecular_weights=species(components%vapor_index)%molecular_weight)
    if (.not.valid) return
    exchange(1:5)=exchange(1:5)+drag
    if (.not.all(ieee_is_finite(exchange))) return
    if(.not.all(ieee_is_finite(external))) return
    parcel=candidate;component_mass=remaining;gas_change=exchange;force_change=external;ok=.true.
  end subroutine

  subroutine liquid_mixture_surface(components,y,molecular_weights,temperature,pressure,surface_x,ok,vapor_enthalpy)
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(in) :: y(:),molecular_weights(:),temperature,pressure
    real(dp), intent(inout) :: surface_x(:)
    real(dp), intent(in), optional :: vapor_enthalpy(:)
    logical, intent(out) :: ok
    real(dp) :: liquid_x(size(y)),candidate(size(y)),inverse_weight,weights(size(y))
    integer :: i
    logical :: valid
    ok=.false.
    if (.not.valid_composition(components,y)) return
    if (size(molecular_weights)/=size(y).or.size(surface_x)/=size(y)) return
    if (.not.all(ieee_is_finite(molecular_weights))) return
    if (any(molecular_weights<1e-3_dp).or.any(molecular_weights>1e5_dp)) return
    if (present(vapor_enthalpy)) then
      if (size(vapor_enthalpy)/=size(y)) return
      if (.not.all(ieee_is_finite(vapor_enthalpy))) return
    end if
    do i=1,size(y)
      weights(i)=resolved_liquid_molecular_weight(components(i),molecular_weights(i))
    end do
    inverse_weight=sum(y/weights);liquid_x=(y/weights)/inverse_weight
    candidate=0
    do i=1,size(y)
      if (y(i)==0) cycle
      if (present(vapor_enthalpy)) then
        call saturation_mole_fraction(components(i),molecular_weights(i),temperature,pressure,candidate(i),valid, &
          liquid_x(i),vapor_enthalpy(i))
      else
        call saturation_mole_fraction(components(i),molecular_weights(i),temperature,pressure,candidate(i),valid,liquid_x(i))
      end if
      if (.not.valid) return
    end do
    if (sum(candidate)>=1) return
    surface_x=candidate;ok=.true.
  end subroutine

  subroutine advance_mixture_evaporation(parcel,components,component_mass,evaporated_mass,vapor_enthalpy,heat,gas_change,ok, &
      euler_temperature,pressure,gas_molecular_weights)
    type(spray_parcel), intent(inout) :: parcel
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(inout) :: component_mass(:)
    real(dp), intent(in) :: evaporated_mass(:),vapor_enthalpy(:),heat
    real(dp), intent(out) :: gas_change(:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: euler_temperature
    real(dp), optional, intent(in) :: pressure,gas_molecular_weights(:)
    real(dp) :: remaining(size(component_mass)),y(size(component_mass)),exchange(size(gas_change))
    real(dp) :: old_mass,new_mass,rho,cp,old_e,new_e,target,t,kinetic,drop_heat,latent_energy,surface_energy,mw
    integer :: i,n
    logical :: valid,use_euler
    type(spray_parcel) :: candidate
    use_euler=.false.
    if(present(euler_temperature)) use_euler=euler_temperature
    ok=.false.;gas_change=0;n=size(component_mass)
    if (n<1.or.n>16.or.size(components)/=n.or.size(evaporated_mass)/=n.or.size(vapor_enthalpy)/=n) return
    if (size(gas_change)<6.or..not.valid_spray_parcel(parcel)) return
    if (.not.all(ieee_is_finite(component_mass)).or..not.all(ieee_is_finite(evaporated_mass))) return
    if (.not.all(ieee_is_finite(vapor_enthalpy)).or..not.ieee_is_finite(heat)) return
    if (any(component_mass<0).or.any(component_mass>1).or.any(evaporated_mass<0)) return
    if(any(evaporated_mass>1e18_dp).or.any(abs(vapor_enthalpy)>1e15_dp)) return
    if(.not.use_euler.and.any(evaporated_mass>component_mass)) return
    old_mass=sum(component_mass)
    if (old_mass<=tiny(1.0_dp).or.parcel%multiplicity<=tiny(1.0_dp)) return
    if (abs(old_mass-parcel%mass)>32*epsilon(1.0_dp)*n*old_mass) return
    if (abs(heat)>min(1e24_dp,1e15_dp*old_mass*parcel%multiplicity)) return
    y=component_mass/old_mass
    call liquid_mixture_properties(components,y,parcel%temperature,rho,cp,old_e,valid,pressure,gas_molecular_weights)
    if (.not.valid) return
    if(use_euler) then
      call euler_component_remainder(component_mass,evaporated_mass,remaining,new_mass)
    else
      remaining=component_mass-evaporated_mass;new_mass=sum(remaining)
    endif
    candidate=parcel;candidate%mass=new_mass;new_e=0
    if (new_mass>0) then
      y=remaining/new_mass;drop_heat=heat/parcel%multiplicity
      if(use_euler) then
        ! Pinned Drag.H evolves temperature with old mass/cp and latent heat.
        ! Charge the gas for the actual accepted liquid-energy change below.
        latent_energy=0
        do i=1,n
          if(evaporated_mass(i)==0) cycle
          mw=components(i)%molecular_weight
          if(present(gas_molecular_weights)) mw=gas_molecular_weights(i)
          call liquid_evaporation_energy(components(i),parcel%temperature,surface_energy,valid,pressure,mw)
          if(.not.valid) return
          latent_energy=latent_energy+evaporated_mass(i)*(vapor_enthalpy(i)-surface_energy)
        end do
        t=parcel%temperature+(drop_heat-latent_energy)/(old_mass*cp)
      else
        target=old_mass*old_e+drop_heat-sum(evaporated_mass*vapor_enthalpy)
        if (abs(target)>1e15_dp*new_mass) return
        target=target/new_mass;t=parcel%temperature
        call liquid_mixture_temperature(components,y,target,t,valid)
        if (.not.valid) return
      end if
      call liquid_mixture_properties(components,y,t,rho,cp,new_e,valid,pressure,gas_molecular_weights)
      if (.not.valid) return
      candidate%temperature=t
    else
      candidate%khrt_initial_count=0;candidate%khrt_shed_volume=0
    end if
    if (.not.valid_spray_parcel(candidate)) return
    exchange=0
    do i=1,n
      if (components(i)%vapor_index>size(exchange)-5) return
      if(use_euler) then
        exchange(5+components(i)%vapor_index)=exchange(5+components(i)%vapor_index)+ &
          parcel%multiplicity*(component_mass(i)-remaining(i))
      else
        exchange(5+components(i)%vapor_index)=exchange(5+components(i)%vapor_index)+ &
          parcel%multiplicity*evaporated_mass(i)
      endif
    end do
    exchange(1)=sum(exchange(6:));exchange(2:4)=exchange(1)*parcel%velocity
    kinetic=.5_dp*sum(parcel%velocity**2)
    exchange(5)=parcel%multiplicity*(old_mass*(old_e+kinetic)-new_mass*(new_e+kinetic))
    if (.not.all(ieee_is_finite(exchange))) return
    parcel=candidate;component_mass=remaining;gas_change=exchange;ok=.true.
  end subroutine

  subroutine advance_mixture_heating(parcel,components,y,heat,gas_energy,ok)
    type(spray_parcel), intent(inout) :: parcel
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(in) :: y(:),heat
    real(dp), intent(out) :: gas_energy
    logical, intent(out) :: ok
    real(dp) :: mass,rho,cp,old_e,new_e,t,target
    logical :: valid
    ok=.false.;gas_energy=0
    if (.not.valid_spray_parcel(parcel).or..not.ieee_is_finite(heat)) return
    mass=parcel%mass*parcel%multiplicity
    if (mass<=tiny(1.0_dp)) return
    if (abs(heat)>min(1e24_dp,1e15_dp*mass)) return
    call liquid_mixture_properties(components,y,parcel%temperature,rho,cp,old_e,valid)
    if (.not.valid) return
    target=old_e+heat/mass;t=parcel%temperature
    call liquid_mixture_temperature(components,y,target,t,valid)
    if (.not.valid) return
    call liquid_mixture_properties(components,y,t,rho,cp,new_e,valid)
    if (.not.valid) return
    ! Debit the actual liquid-energy increment, including inversion roundoff.
    gas_energy=mass*(old_e-new_e);parcel%temperature=t;ok=.true.
  end subroutine
end module
