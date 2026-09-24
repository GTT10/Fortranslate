! Ideal liquid mixing rules adapted from PelePhysics SprayPropertiesPeleMP.H,
! revision 6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b.
! Copyright/redistribution: references/pelephysics-spray-license.txt.
module spray_mixture_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species,nasa7_mass_properties
  use gas_transport_mod, only: gas_transport_species
  use mixture_thermo_mod, only: mixture_mass_properties
  use mixture_transport_mod, only: mixture_transport_coefficients,binary_diffusion_coefficients
  use spray_parcel_mod
  implicit none
  private
  real(dp), parameter, public :: mixture_trace_fraction=1e-12_dp
  character(len=*), parameter, public :: mixture_evolution_id='mixed-liquid-components-v2-terminal'
  public :: liquid_mixture_properties, liquid_mixture_temperature, advance_mixture_heating
  public :: liquid_mixture_surface, advance_mixture_evaporation
  public :: build_mixture_film, mixture_transfer_rates, advance_mixture_frozen_gas, valid_component_inventory
contains
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
        if (components(i)%vapor_index==components(j)%vapor_index) return
      end do
    end do
    ok=.true.
  end function

  pure subroutine liquid_mixture_properties(components,y,temperature,density,cp,energy,ok)
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(in) :: y(:),temperature
    real(dp), intent(inout) :: density,cp,energy
    logical, intent(out) :: ok
    real(dp) :: volume,candidate_cp,candidate_e,rhoi,cpi,ei
    integer :: i
    ok=.false.
    if (.not.valid_composition(components,y).or..not.ieee_is_finite(temperature)) return
    if (temperature<1.or.temperature>1e5_dp) return
    volume=0;candidate_cp=0;candidate_e=0
    do i=1,size(y)
      if (y(i)==0) cycle
      rhoi=liquid_density_at_temperature(components(i),temperature)
      cpi=liquid_heat_capacity_at_temperature(components(i),temperature)
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
    real(dp) :: diffusion(size(species)),binary(size(species),size(species)),tf,r,mw,cp,cv,gamma,h,e,entropy
    real(dp) :: carrier_mw,carrier_y,denom,sum_x,carrier_sum,resistance,di,cp_i
    logical :: active(size(species)),valid
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
    xg=(yg/mw_g)/sum(yg/mw_g)
    do i=1,size(y)
      v=components(i)%vapor_index
      if (v>size(species)) return
      if (trim(components(i)%vapor_name)/=trim(species(v)%name)) return
      mw_liq(i)=mw_g(v)
      call nasa7_mass_properties(species(v),parcel%temperature,cp_i,cv,h,e,entropy,valid)
      if (.not.valid) return
      enthalpy(i)=h
    end do
    xs=0;active=.false.
    if (options%evaporation) then
      call liquid_mixture_surface(components,y,mw_liq,parcel%temperature,primitive(5),xs,valid,enthalpy)
      if (.not.valid) return
      do i=1,size(y)
        v=components(i)%vapor_index
        ! No condensation: retain only components with outward vapor driving force.
        if (xs(i)<=xg(v)) xs(i)=0
        active(v)=xs(i)>0
      end do
    end if
    sum_x=sum(xs);yf=yg;ys=0
    if (sum_x>0) then
      carrier_y=sum(yg,mask=.not.active);denom=sum(yg/mw_g,mask=.not.active)
      if (carrier_y<=0.or.denom<=0) return
      carrier_mw=carrier_y/denom
      denom=sum(xs*mw_liq)+(1-sum_x)*carrier_mw
      do i=1,size(y)
        v=components(i)%vapor_index
        ys(v)=xs(i)*mw_liq(i)/denom
        if (active(v)) yf(v)=(2*ys(v)+yg(v))/3
      end do
      where (.not.active) yf=yg*(1-sum(yf,mask=active))/carrier_y
    end if
    tf=(2*parcel%temperature+gas_temperature)/3
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
    call mixture_transport_coefficients(species,transport,yf,tf,primitive(5), &
      candidate%viscosity,candidate%conductivity,diffusion,valid)
    if (.not.valid) return
    call binary_diffusion_coefficients(species,transport,tf,primitive(5),binary,valid)
    if (.not.valid) return
    diff_weights=0
    do i=1,size(y)
      if (xs(i)==0) cycle
      v=components(i)%vapor_index;carrier_sum=0;resistance=0
      do j=1,size(species)
        if (j==v) cycle
        carrier_sum=carrier_sum+yf(j)/mw_g(j)
        resistance=resistance+yf(j)/(mw_g(j)*binary(v,j))
      end do
      di=diffusion(v)
      if (carrier_sum>0.and.resistance>0) di=carrier_sum/resistance
      diff_weights(i)=xs(i)*di
    end do
    candidate%diffusivity=diffusion(components(1)%vapor_index)
    if (sum_x>0) candidate%diffusivity=sum(diff_weights)/sum_x
    candidate%surface_fraction=sum(ys);candidate%ambient_fraction=sum(yg,mask=active)
    candidate%gas_temperature=gas_temperature;candidate%gas_velocity=primitive(2:4)
    candidate%vapor_enthalpy=sum(y*enthalpy)
    if (sum(diff_weights)>0) mass_weights=diff_weights/sum(diff_weights)
    vapor_h=enthalpy;film=candidate;ok=.true.
  end subroutine

  subroutine mixture_transfer_rates(species,transport,primitive,gas_temperature,parcel,components,y,options, &
      inverse_tau,heat_coefficient,mass_rate,vapor_h,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_parcel), intent(in) :: parcel
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    real(dp), intent(in) :: primitive(:),gas_temperature,y(:)
    real(dp), intent(out) :: inverse_tau,heat_coefficient,mass_rate(:),vapor_h(:)
    logical, intent(out) :: ok
    type(spray_film) :: film
    type(spray_liquid) :: effective
    real(dp) :: weights(size(y)),enthalpy(size(y)),rho,cp,e,mass_coefficient,rate,heat,d
    logical :: valid
    inverse_tau=0;heat_coefficient=0;mass_rate=0;vapor_h=0;ok=.false.
    if (size(mass_rate)/=size(y).or.size(vapor_h)/=size(y)) return
    call liquid_mixture_properties(components,y,parcel%temperature,rho,cp,e,valid)
    if (.not.valid) return
    call build_mixture_film(species,transport,primitive,gas_temperature,parcel,components,y,options, &
      film,weights,enthalpy,valid)
    if (.not.valid) return
    ! Only density and cp enter the common transfer correlations; energy stays component-resolved.
    effective=components(1);effective%density_model='constant';effective%density=rho
    effective%heat_capacity_model='constant';effective%heat_capacity=cp
    effective%density_coefficients=0;effective%density_temperature_range=0
    effective%heat_capacity_coefficients=0;effective%heat_capacity_temperature_range=0
    call spray_transfer_rates(parcel,effective,film,options,rate,heat,mass_coefficient,valid)
    if (.not.valid) return
    d=parcel_diameter(parcel,effective)
    mass_rate=mass_coefficient*d*weights;inverse_tau=rate;heat_coefficient=heat;vapor_h=enthalpy;ok=.true.
  end subroutine

  subroutine advance_mixture_frozen_gas(species,transport,primitive,gas_temperature,parcel,components,component_mass, &
      options,dt,gas_change,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_parcel), intent(inout) :: parcel
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    real(dp), intent(in) :: primitive(:),gas_temperature,dt
    real(dp), intent(inout) :: component_mass(:)
    real(dp), intent(out) :: gas_change(:)
    logical, intent(out) :: ok
    type(spray_parcel) :: candidate
    real(dp) :: remaining(size(component_mass)),y(size(component_mass)),rates(size(component_mass)),hv(size(component_mass))
    real(dp) :: exchange(size(gas_change)),drag(5),invtau,heat,mass
    real(dp) :: evaporated(size(component_mass)),remaining_mass,rho,cp,e
    integer :: i
    logical :: valid
    ok=.false.;gas_change=0
    if (.not.ieee_is_finite(dt).or.dt<0.or.dt>1e6_dp) return
    if (size(gas_change)/=5+size(species).or..not.valid_spray_parcel(parcel)) return
    if (.not.all(ieee_is_finite(component_mass)).or.any(component_mass<0).or.any(component_mass>1)) return
    mass=sum(component_mass)
    if (mass<=tiny(1.0_dp)) return
    if (abs(mass-parcel%mass)>32*epsilon(1.0_dp)*size(component_mass)*mass) return
    y=component_mass/mass
    call mixture_transfer_rates(species,transport,primitive,gas_temperature,parcel,components,y,options, &
      invtau,heat,rates,hv,valid)
    if (.not.valid) return
    if (dt==0) then
      ok=.true.;return
    end if
    evaporated=rates*dt
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
    if (any(rates>0).and.remaining_mass>0) then
      call liquid_mixture_properties(components,remaining/remaining_mass,parcel%temperature,rho,cp,e,valid)
      if (.not.valid) return
      if (remaining_mass<acos(-1.0_dp)*rho*options%minimum_diameter**3/6) evaporated=component_mass
    end if
    candidate=parcel;remaining=component_mass
    if (options%drag) candidate%velocity=primitive(2:4)+(parcel%velocity-primitive(2:4))*exp(-min(700.0_dp,invtau*dt))
    candidate%position=parcel%position+.5_dp*dt*(parcel%velocity+candidate%velocity)
    drag=0;drag(2:4)=mass*parcel%multiplicity*(parcel%velocity-candidate%velocity)
    drag(5)=.5_dp*mass*parcel%multiplicity*(sum(parcel%velocity**2)-sum(candidate%velocity**2))
    call advance_mixture_evaporation(candidate,components,remaining,evaporated,hv, &
      heat*(gas_temperature-parcel%temperature)*dt*parcel%multiplicity,exchange,valid)
    if (.not.valid) return
    exchange(1:5)=exchange(1:5)+drag
    if (.not.all(ieee_is_finite(exchange))) return
    parcel=candidate;component_mass=remaining;gas_change=exchange;ok=.true.
  end subroutine

  subroutine liquid_mixture_surface(components,y,molecular_weights,temperature,pressure,surface_x,ok,vapor_enthalpy)
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(in) :: y(:),molecular_weights(:),temperature,pressure
    real(dp), intent(inout) :: surface_x(:)
    real(dp), intent(in), optional :: vapor_enthalpy(:)
    logical, intent(out) :: ok
    real(dp) :: liquid_x(size(y)),candidate(size(y)),inverse_weight
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
    inverse_weight=sum(y/molecular_weights);liquid_x=(y/molecular_weights)/inverse_weight
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

  subroutine advance_mixture_evaporation(parcel,components,component_mass,evaporated_mass,vapor_enthalpy,heat,gas_change,ok)
    type(spray_parcel), intent(inout) :: parcel
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(inout) :: component_mass(:)
    real(dp), intent(in) :: evaporated_mass(:),vapor_enthalpy(:),heat
    real(dp), intent(out) :: gas_change(:)
    logical, intent(out) :: ok
    real(dp) :: remaining(size(component_mass)),y(size(component_mass)),exchange(size(gas_change))
    real(dp) :: old_mass,new_mass,rho,cp,old_e,new_e,target,t,kinetic,drop_heat
    integer :: i,n
    logical :: valid
    type(spray_parcel) :: candidate
    ok=.false.;gas_change=0;n=size(component_mass)
    if (n<1.or.n>16.or.size(components)/=n.or.size(evaporated_mass)/=n.or.size(vapor_enthalpy)/=n) return
    if (size(gas_change)<6.or..not.valid_spray_parcel(parcel)) return
    if (.not.all(ieee_is_finite(component_mass)).or..not.all(ieee_is_finite(evaporated_mass))) return
    if (.not.all(ieee_is_finite(vapor_enthalpy)).or..not.ieee_is_finite(heat)) return
    if (any(component_mass<0).or.any(component_mass>1).or.any(evaporated_mass<0)) return
    if (any(evaporated_mass>component_mass).or.any(abs(vapor_enthalpy)>1e15_dp)) return
    old_mass=sum(component_mass)
    if (old_mass<=tiny(1.0_dp).or.parcel%multiplicity<=tiny(1.0_dp)) return
    if (abs(old_mass-parcel%mass)>32*epsilon(1.0_dp)*n*old_mass) return
    if (abs(heat)>min(1e24_dp,1e15_dp*old_mass*parcel%multiplicity)) return
    y=component_mass/old_mass
    call liquid_mixture_properties(components,y,parcel%temperature,rho,cp,old_e,valid)
    if (.not.valid) return
    remaining=component_mass-evaporated_mass;new_mass=sum(remaining)
    candidate=parcel;candidate%mass=new_mass;new_e=0
    if (new_mass>0) then
      y=remaining/new_mass;drop_heat=heat/parcel%multiplicity
      target=old_mass*old_e+drop_heat-sum(evaporated_mass*vapor_enthalpy)
      if (abs(target)>1e15_dp*new_mass) return
      target=target/new_mass;t=parcel%temperature
      call liquid_mixture_temperature(components,y,target,t,valid)
      if (.not.valid) return
      call liquid_mixture_properties(components,y,t,rho,cp,new_e,valid)
      if (.not.valid) return
      candidate%temperature=t
    else
      candidate%khrt_initial_count=0;candidate%khrt_shed_volume=0
    end if
    if (.not.valid_spray_parcel(candidate)) return
    exchange=0
    do i=1,n
      if (components(i)%vapor_index>size(exchange)-5) return
      exchange(5+components(i)%vapor_index)=parcel%multiplicity*evaporated_mass(i)
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
