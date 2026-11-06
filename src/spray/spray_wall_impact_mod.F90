! Ahamed impact correlations from pinned PelePhysics; see references/spray_wall_impact.json.
module spray_wall_impact_mod
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use spray_mixture_mod, only: liquid_mixture_properties,liquid_mixture_temperature,valid_component_inventory, &
    build_mixture_film,mixture_trace_fraction,liquid_mixture_conductivity,liquid_mixture_critical_limit, &
    apply_mixture_breakup_viscosity
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use spray_parcel_mod
  use nasa7_thermo_mod, only: universal_gas_constant
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer, parameter, public :: impact_rebound=0,impact_deposit=1,impact_dry=2
  integer, parameter, public :: impact_thermal=3,impact_spread=4,impact_wet=5
  public :: classify_wall_impact, splash_velocity_components, splash_mass_fractions
  public :: pelec_splash_critical_parameter
  type, public :: spray_wall_film
    real(dp) :: position(3)=0,normal(3)=0,mass=0,temperature=0,diameter=0,height=0
    integer :: material_id=1
  end type
  type, public :: spray_wall_options
    logical :: enabled(6)=.false.
    real(dp) :: temperature(6)=300,contact_angle=1.0471975511965977_dp
    integer :: maximum_films=1000000
    logical :: evolve=.false.
    real(dp) :: liquid_conductivity=0,minimum_height=1e-7_dp
    character(len=16) :: conductivity_model='constant',port_conductivity_models(16)=''
    real(dp) :: conductivity_coefficients(4)=0,conductivity_temperature_range(2)=0
    real(dp) :: port_conductivities(16)=-1,port_conductivity_coefficients(4,16)=0
    real(dp) :: port_conductivity_temperature_ranges(2,16)=0
  end type
  character(len=*), parameter, public :: mixture_film_evolution_id='mixed-film-v7-critical-limit'
  public :: valid_wall_component_inventory
  public :: mixture_impact_boiling_temperature,mixture_impact_at_planar_wall
  public :: impact_at_planar_wall,valid_wall_options,valid_wall_film
  public :: wall_film_transfer_rates,advance_wall_film_frozen,wall_liquid_conductivity_at_temperature
  public :: mixture_wall_film_rates,advance_mixture_wall_film_frozen
  public :: resolve_wall_impact,resolve_mixture_wall_impact,advance_mixture_wall_film_exchange
  public :: construct_splash_children, construct_wall_film, wall_film_totals
  real(dp), parameter :: pi=acos(-1.0_dp)
  public :: count_periodic_film_images

contains
  pure real(dp) function pelec_splash_critical_parameter(pressure,wall_temperature) result(kcrit)
    real(dp), intent(in) :: pressure,wall_temperature
    ! Pinned AhamedSplash.H uses Pa*1e-6 for PeleC. Pa*1e-5 is
    ! the separate PELELM_USE_SPRAY branch, despite its p_bar variable name.
    kcrit=(26.788_dp*log(pressure*1e-6_dp)+33.546_dp)*exp(-.006_dp*(wall_temperature-260))+36
  end function pelec_splash_critical_parameter

  pure subroutine count_periodic_film_images(displacement,normal,radius,tolerance,length,count,ok)
    real(dp), intent(in) :: displacement(3),normal(3),radius,tolerance,length(3)
    integer, intent(inout) :: count
    logical, intent(out) :: ok
    real(dp) :: offset(3),delta(3)
    integer :: first(3),last(3),axis,i,j,k,found
    integer(int64) :: images
    ok=.false.
    if(.not.all(ieee_is_finite([displacement,normal,radius,tolerance,length]))) return
    if(any(length<0).or.any(length>1e9_dp).or.any(length>0.and.length<1e-12_dp)) return
    if(radius<0.or.tolerance<0) return
    offset=displacement;first=0;last=0;images=1;found=0
    do axis=1,3
      if(length(axis)==0) cycle
      offset(axis)=modulo(offset(axis)+.5_dp*length(axis),length(axis))-.5_dp*length(axis)
      if(radius/length(axis)>500000) return
      first(axis)=ceiling((-radius-offset(axis))/length(axis))
      last(axis)=floor((radius-offset(axis))/length(axis))
      images=images*int(max(0,last(axis)-first(axis)+1),int64)
      if(images>1000000_int64) return
    end do
    do k=first(3),last(3)
    do j=first(2),last(2)
    do i=first(1),last(1)
      delta=offset+real([i,j,k],dp)*length
      if(abs(dot_product(delta,normal))>tolerance.or.norm2(delta)>radius) cycle
      found=found+1
    end do
    end do
    end do
    count=found;ok=.true.
  end subroutine count_periodic_film_images

  real(dp) function wall_liquid_conductivity_at_temperature(o,material_id,t) result(value)
    type(spray_wall_options), intent(in) :: o
    integer, intent(in) :: material_id
    real(dp), intent(in) :: t
    character(len=16) :: model
    real(dp) :: coefficients(4),bounds(2),constant
    value=0
    if (material_id<1.or.material_id>16.or..not.ieee_is_finite(t)) return
    if (t<1.or.t>1e5_dp) return
    model=o%conductivity_model;coefficients=o%conductivity_coefficients
    bounds=o%conductivity_temperature_range;constant=o%liquid_conductivity
    if (o%port_conductivities(material_id)/=-1) constant=o%port_conductivities(material_id)
    if (o%port_conductivity_models(material_id)/='') then
      model=o%port_conductivity_models(material_id)
      coefficients=o%port_conductivity_coefficients(:,material_id)
      bounds=o%port_conductivity_temperature_ranges(:,material_id)
    end if
    select case(model)
    case('constant')
      value=constant
    case('polynomial')
      if (t<bounds(1).or.t>bounds(2)) return
      if (.not.all(ieee_is_finite(coefficients))) return
      if (maxval(abs(coefficients))>1e20_dp) return
      value=coefficients(1)+t*(coefficients(2)+t*(coefficients(3)+t*coefficients(4)))
    end select
    if (.not.ieee_is_finite(value)) value=0
    if (value<=0.or.value>1e4_dp) value=0
  end function

  logical function valid_conductivity_model(model,constant,coefficients,bounds,required) result(ok)
    character(len=*), intent(in) :: model
    real(dp), intent(in) :: constant,coefficients(4),bounds(2)
    logical, intent(in) :: required
    type(spray_wall_options) :: probe
    real(dp) :: a,b,c,disc,q,roots(2)
    integer :: i
    ok=.false.
    if (.not.all(ieee_is_finite([constant,coefficients,bounds]))) return
    select case(model)
    case('constant')
      if (any(coefficients/=0).or.any(bounds/=0)) return
      if (required.and.(constant<=0.or.constant>1e4_dp)) return
    case('polynomial')
      if (maxval(abs(coefficients))>1e20_dp) return
      if (bounds(1)<1.or.bounds(2)>1e5_dp.or.bounds(2)<=bounds(1)) return
      probe%conductivity_model=model;probe%conductivity_coefficients=coefficients
      probe%conductivity_temperature_range=bounds
      do i=1,2
        if (wall_liquid_conductivity_at_temperature(probe,1,bounds(i))<=0) return
      end do
      ! A cubic can become negative between positive endpoints: inspect stationary points.
      a=3*coefficients(4);b=2*coefficients(3);c=coefficients(2);roots=bounds(1)
      if (a/=0) then
        disc=b*b-4*a*c
        if (disc>=0) then
          q=-.5_dp*(b+sign(sqrt(disc),b))
          if (abs(q)<=abs(a)*1e5_dp) roots(1)=q/a
          if (q/=0) then
            if (abs(c)<=abs(q)*1e5_dp) roots(2)=c/q
          end if
        end if
      else if (b/=0) then
        if (abs(c)<=abs(b)*1e5_dp) roots(1)=-c/b
      end if
      do i=1,2
        if (roots(i)<bounds(1).or.roots(i)>bounds(2)) cycle
        if (wall_liquid_conductivity_at_temperature(probe,1,roots(i))<=0) return
      end do
    case default
      return
    end select
    ok=.true.
  end function

  subroutine wall_film_transfer_rates(patch,liquid,gas,conductivity,wall_temperature,heat,evaporation, &
      mass_rate,gas_heat,wall_heat,temperature_rate,ok)
    type(spray_wall_film), intent(in) :: patch
    type(spray_liquid), intent(in) :: liquid
    type(spray_film), intent(in) :: gas
    real(dp), intent(in) :: conductivity,wall_temperature
    logical, intent(in) :: heat,evaporation
    real(dp), intent(out) :: mass_rate,gas_heat,wall_heat,temperature_rate
    logical, intent(out) :: ok
    real(dp) :: area,h,bm,cp,latent,rho,liquid_lambda,surface_energy
    mass_rate=0;gas_heat=0;wall_heat=0;temperature_rate=0;ok=.false.
    if (.not.valid_wall_film(patch).or..not.valid_spray_liquid(liquid)) return
    if (.not.all(ieee_is_finite([conductivity,wall_temperature,gas%conductivity,gas%cp,gas%density, &
        gas%diffusivity,gas%gas_temperature,gas%surface_fraction,gas%ambient_fraction,gas%vapor_enthalpy]))) return
    liquid_lambda=conductivity
    if(liquid%conductivity_model/='supplied') liquid_lambda= &
      liquid_conductivity_at_temperature(liquid,patch%temperature,gas%pressure,gas%liquid_molecular_weight)
    if (liquid_lambda<=0.or.liquid_lambda>1e4_dp.or.wall_temperature<1.or.wall_temperature>1e5_dp) return
    if (gas%conductivity<=0.or.gas%conductivity>1e4_dp.or.gas%cp<1.or.gas%cp>1e8_dp) return
    if (gas%density<=0.or.gas%density>1e8_dp.or.gas%diffusivity<=0.or.gas%diffusivity>1e4_dp) return
    if (gas%gas_temperature<1.or.gas%gas_temperature>1e5_dp) return
    if (gas%surface_fraction<0.or.gas%surface_fraction>=1.or.gas%ambient_fraction<0.or.gas%ambient_fraction>1) return
    call liquid_transfer_properties(liquid,patch%temperature,gas%pressure,gas%liquid_molecular_weight,rho,cp,ok)
    if(.not.ok) return
    ok=.false.
    if (cp<=0.or.patch%height<1e-15_dp.or.patch%diameter>1e9_dp) return
    area=.25_dp*pi*patch%diameter**2;h=gas%conductivity/patch%height
    if (heat) then
      gas_heat=area*h*(gas%gas_temperature-patch%temperature)
      wall_heat=2*area*liquid_lambda*(wall_temperature-patch%temperature)/patch%height
    end if
    if (evaporation) then
      bm=min(20.0_dp,max(0.0_dp,(gas%surface_fraction-gas%ambient_fraction)/(1-gas%surface_fraction)))
      mass_rate=area*h*(gas%conductivity/gas%cp)**(1.0_dp/3)*positive_log1p(bm)/gas%conductivity* &
        (gas%density*gas%diffusivity)**(2.0_dp/3)
    end if
    surface_energy=liquid_energy(liquid,patch%temperature)
    if(liquid%boiling_model=='pelemp-cap') then
      call liquid_evaporation_energy(liquid,patch%temperature,surface_energy,ok,gas%pressure,gas%liquid_molecular_weight)
      if(.not.ok) return
      ok=.false.
    endif
    latent=gas%vapor_enthalpy-surface_energy
    if (mass_rate>0.and.latent<=0) return
    temperature_rate=(gas_heat+wall_heat-mass_rate*latent)/(patch%mass*cp)
    ok=all(ieee_is_finite([mass_rate,gas_heat,wall_heat,temperature_rate]))
  end subroutine

  subroutine mixture_wall_film_rates(patch,components,component_mass,species,transport,primitive,gas_temperature, &
      options,conductivity,wall_temperature,mass_rate,vapor_h,gas_heat,wall_heat,temperature_rate,ok)
    type(spray_wall_film), intent(in) :: patch
    type(spray_liquid), intent(in) :: components(:)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_options), intent(in) :: options
    real(dp), intent(in) :: component_mass(:),primitive(:),gas_temperature,conductivity,wall_temperature
    real(dp), intent(out) :: mass_rate(:),vapor_h(:),gas_heat,wall_heat,temperature_rate
    logical, intent(out) :: ok
    type(spray_parcel) :: identity
    type(spray_film) :: skin
    type(spray_liquid) :: effective
    type(spray_wall_film) :: thermal_patch
    real(dp) :: y(size(components)),weights(size(components)),hv(size(components)),rho,cp,e,mdot,qg,qw,tdot,latent,liquid_lambda
    real(dp) :: component_rates(size(components)),surface_energy
    integer :: i,n
    logical :: valid
    mass_rate=0;vapor_h=0;gas_heat=0;wall_heat=0;temperature_rate=0;ok=.false.;n=size(components)
    if (.not.valid_wall_film(patch).or.size(component_mass)/=n.or.size(mass_rate)/=n.or.size(vapor_h)/=n) return
    if (.not.all(ieee_is_finite(component_mass))) return
    if (any(component_mass<0).or.any(component_mass>1e18_dp)) return
    if (abs(sum(component_mass)-patch%mass)>32*epsilon(1.0_dp)*n*patch%mass) return
    y=component_mass/patch%mass
    identity=spray_parcel(1_int64,patch%position,[0.0_dp,0.0_dp,0.0_dp], &
      min(1.0_dp,patch%mass),patch%temperature,max(1.0_dp,patch%mass))
    thermal_patch=patch
    if(options%mixture_time_scheme=='euler') then
      call liquid_mixture_critical_limit(components,y,patch%temperature,thermal_patch%temperature,valid)
      if(.not.valid) return
      identity%temperature=thermal_patch%temperature
    endif
    call build_mixture_film(species,transport,primitive,gas_temperature,identity,components,y,options, &
      skin,weights,hv,valid,wall_skin=.true.)
    if (.not.valid) return
    call liquid_mixture_properties(components,y,patch%temperature,rho,cp,e,valid,primitive(5), &
      species(components%vapor_index)%molecular_weight)
    if (.not.valid) return
    liquid_lambda=0
    call liquid_mixture_conductivity(components,y,patch%temperature,conductivity,liquid_lambda,valid,primitive(5), &
      species(components%vapor_index)%molecular_weight)
    if(.not.valid) return
    effective=components(1);effective%density=rho;effective%heat_capacity=cp
    effective%reference_temperature=patch%temperature;effective%reference_energy=e
    effective%density_model='constant';effective%density_coefficients=0;effective%density_temperature_range=0
    effective%heat_capacity_model='constant';effective%heat_capacity_coefficients=0
    effective%heat_capacity_temperature_range=0;effective%boiling_model='reject';effective%critical_temperature=0
    effective%normal_boiling_temperature=0
    effective%conductivity_model='supplied';effective%conductivity=0
    effective%conductivity_coefficients=0;effective%conductivity_temperature_range=0
    ! Pinned WallFilm.H encloses both heat terms and evaporation in
    ! sumXVap>0. Zero diffusion weights mean no retained surface vapor.
    if(options%mixture_time_scheme=='euler'.and.options%evaporation.and.sum(weights)==0) then
      vapor_h=hv;ok=.true.;return
    endif
    if (sum(weights)>0) skin%vapor_enthalpy=sum(weights*hv)
    call wall_film_transfer_rates(thermal_patch,effective,skin,liquid_lambda,wall_temperature,options%heat,options%evaporation, &
      mdot,qg,qw,tdot,valid)
    if (.not.valid) return
    ! WallFilm.H applies the 2/3 power after each component's vapor share.
    ! The film builder returns weights proportional to X_v,i * rhoD_i;
    ! multiplying the aggregate rate by weights**(2/3) recovers that rule,
    ! including multiple liquid components mapped to one gas species.
    component_rates=mdot*weights**(2._dp/3)
    latent=0
    do i=1,n
      if (weights(i)==0) cycle
      call liquid_evaporation_energy(components(i),thermal_patch%temperature,surface_energy,valid,primitive(5), &
        species(components(i)%vapor_index)%molecular_weight)
      if(.not.valid) return
      latent=latent+component_rates(i)*(hv(i)-surface_energy)
    end do
    tdot=(qg+qw-latent)/(patch%mass*cp)
    if (.not.ieee_is_finite(tdot)) return
    mass_rate=component_rates;vapor_h=hv;gas_heat=qg;wall_heat=qw;temperature_rate=tdot;ok=.true.
  end subroutine

  subroutine advance_mixture_wall_film_frozen(patch,components,component_mass,species,transport,primitive,gas_temperature, &
      options,conductivity,wall_temperature,dt,minimum_height,gas_change,wall_energy,ok)
    type(spray_wall_film), intent(inout) :: patch
    type(spray_liquid), intent(in) :: components(:)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_options), intent(in) :: options
    real(dp), intent(inout) :: component_mass(:)
    real(dp), intent(in) :: primitive(:),gas_temperature,conductivity,wall_temperature,dt,minimum_height
    real(dp), intent(out) :: gas_change(:),wall_energy
    logical, intent(out) :: ok
    real(dp) :: rates(size(components)),hv(size(components)),qgas,qwall,tdot
    logical :: valid
    gas_change=0;wall_energy=0;ok=.false.
    if (.not.ieee_is_finite(dt)) return
    if (dt<0.or.dt>1e6_dp) return
    call mixture_wall_film_rates(patch,components,component_mass,species,transport,primitive,gas_temperature, &
      options,conductivity,wall_temperature,rates,hv,qgas,qwall,tdot,valid)
    if (.not.valid) return
    if(dt==0) then
      ok=.true.;return
    endif
    call advance_mixture_wall_film_exchange(patch,components,component_mass,rates*dt,hv,qgas*dt,qwall*dt, &
      minimum_height,gas_change,wall_energy,ok,deplete_traces=.true.,pressure=primitive(5), &
      gas_molecular_weights=species(components%vapor_index)%molecular_weight, &
      euler_temperature=options%mixture_time_scheme=='euler')
  end subroutine

  ! Supplied component losses are total kg; heat amounts are J into the patch.
  subroutine advance_mixture_wall_film_exchange(patch,components,component_mass,evaporated_mass,vapor_enthalpy, &
      gas_heat,wall_heat,minimum_height,gas_change,wall_energy,ok,deplete_traces,pressure,gas_molecular_weights, &
      euler_temperature)
    type(spray_wall_film), intent(inout) :: patch
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(inout) :: component_mass(:)
    real(dp), intent(in) :: evaporated_mass(:),vapor_enthalpy(:),gas_heat,wall_heat,minimum_height
    real(dp), intent(out) :: gas_change(:),wall_energy
    logical, intent(out) :: ok
    logical, optional, intent(in) :: deplete_traces,euler_temperature
    real(dp), optional, intent(in) :: pressure,gas_molecular_weights(:)
    type(spray_wall_film) :: next
    real(dp) :: remaining(size(component_mass)),loss(size(component_mass)),exchange(size(gas_change))
    real(dp) :: rho,cp,old_e,new_e,target,t,area,total,liquid_lambda,old_rho,latent,surface_energy,mw,base_temperature
    integer :: n,i
    logical :: valid,use_euler
    gas_change=0;wall_energy=0;ok=.false.;n=size(components)
    use_euler=.false.
    if(present(euler_temperature)) use_euler=euler_temperature
    if (n<1.or.n>16.or.size(component_mass)/=n.or.size(evaporated_mass)/=n.or.size(vapor_enthalpy)/=n) return
    if (.not.valid_wall_film(patch).or.size(gas_change)<6) return
    if (.not.all(ieee_is_finite(component_mass)).or..not.all(ieee_is_finite(evaporated_mass))) return
    if (.not.all(ieee_is_finite(vapor_enthalpy)).or. &
        .not.all(ieee_is_finite([gas_heat,wall_heat,minimum_height]))) return
    if (any(component_mass<0).or.any(component_mass>1e18_dp).or.any(evaporated_mass<0)) return
    if (any(evaporated_mass>1e18_dp)) return
    if (any(abs(vapor_enthalpy)>1e15_dp)) return
    if (max(abs(gas_heat),abs(wall_heat))>min(1e24_dp,1e15_dp*patch%mass)) return
    if (minimum_height<1e-15_dp.or.minimum_height>1) return
    total=sum(component_mass)
    if (abs(total-patch%mass)>32*epsilon(1.0_dp)*n*patch%mass) return
    call liquid_mixture_properties(components,component_mass/patch%mass,patch%temperature,rho,cp,old_e,valid, &
      pressure,gas_molecular_weights)
    if (.not.valid) return
    old_rho=rho
    base_temperature=patch%temperature
    if(use_euler) then
      call liquid_mixture_critical_limit(components,component_mass/patch%mass,patch%temperature,base_temperature,valid)
      if(.not.valid) return
    endif
    do i=1,n
      if (components(i)%vapor_index>size(gas_change)-5) return
    end do
    if (all(evaporated_mass==0).and.gas_heat==0.and.wall_heat==0.and.base_temperature==patch%temperature) then
      ok=.true.;return
    end if
    loss=evaporated_mass
    if (present(deplete_traces)) then
      if (deplete_traces.and..not.use_euler) then
        do i=1,n
          if (loss(i)<=0) cycle
          if (component_mass(i)<=mixture_trace_fraction*patch%mass) then
            loss(i)=component_mass(i)
          else if (loss(i)<=component_mass(i)) then
            if (component_mass(i)-loss(i)<=mixture_trace_fraction*patch%mass) loss(i)=component_mass(i)
          end if
        end do
      end if
    end if
    next=patch;new_e=0
    if(use_euler) then
      ! WallFilm.H first advances total mass, then clips component fractions
      ! and normalizes survivors. Preserve that total, including trace loss.
      next%mass=max(0._dp,patch%mass-sum(loss))
      remaining=max(0._dp,component_mass-loss)
      if(next%mass>0) then
        remaining=min(next%mass,remaining)
        where(remaining/next%mass<mixture_trace_fraction) remaining=0
        total=sum(remaining)
        if(total>0) then
          remaining=remaining*(next%mass/total)
        else
          next%mass=0
        endif
      endif
      if(next%mass==0) remaining=0
    else
      if (any(loss>component_mass)) return
      remaining=component_mass-loss;next%mass=sum(remaining)
    endif
    if (next%mass>0) then
      if(use_euler) then
        latent=0
        do i=1,n
          if(loss(i)==0) cycle
          mw=components(i)%molecular_weight
          if(present(gas_molecular_weights)) mw=gas_molecular_weights(i)
          call liquid_evaporation_energy(components(i),base_temperature,surface_energy,valid,pressure,mw)
          if(.not.valid) return
          latent=latent+loss(i)*(vapor_enthalpy(i)-surface_energy)
        enddo
        t=base_temperature+(gas_heat+wall_heat-latent)/(patch%mass*cp)
      else
        target=patch%mass*old_e+gas_heat+wall_heat-sum(loss*vapor_enthalpy)
        if (abs(target)>1e15_dp*next%mass) return
        t=patch%temperature
        call liquid_mixture_temperature(components,remaining/next%mass,target/next%mass,t,valid)
        if (.not.valid) return
      endif
      call liquid_mixture_properties(components,remaining/next%mass,t,rho,cp,new_e,valid,pressure,gas_molecular_weights)
      if (.not.valid) return
      if(any(components%conductivity_model/='supplied'.and.remaining>0)) then
        liquid_lambda=0
        call liquid_mixture_conductivity(components,remaining/next%mass,t,1._dp,liquid_lambda,valid, &
          pressure,gas_molecular_weights)
        if(.not.valid) return
      endif
      area=.25_dp*pi*patch%diameter**2
      if (.not.ieee_is_finite(area).or.area<=0) return
      next%temperature=t;next%height=next%mass/rho/area
      if(use_euler) next%height=next%mass/old_rho/area
      if (any(evaporated_mass>0).and.next%height<=minimum_height) then
        next%mass=0;remaining=0;new_e=0
      end if
    end if
    if (next%mass==0) then
      next%height=0;next%temperature=patch%temperature
    else
      if (.not.valid_wall_film(next)) return
    end if
    loss=component_mass-remaining;exchange=0
    do i=1,n
      exchange(5+components(i)%vapor_index)=exchange(5+components(i)%vapor_index)+loss(i)
    end do
    exchange(1)=sum(loss);exchange(5)=patch%mass*old_e-next%mass*new_e+wall_heat
    if (.not.all(ieee_is_finite(exchange))) return
    patch=next;component_mass=remaining;gas_change=exchange;wall_energy=-wall_heat;ok=.true.
  end subroutine

  subroutine advance_wall_film_frozen(patch,liquid,gas,conductivity,wall_temperature,heat,evaporation,dt, &
      minimum_height,next,gas_change,wall_energy,ok)
    type(spray_wall_film), intent(in) :: patch
    type(spray_liquid), intent(in) :: liquid
    type(spray_film), intent(in) :: gas
    real(dp), intent(in) :: conductivity,wall_temperature,dt,minimum_height
    logical, intent(in) :: heat,evaporation
    type(spray_wall_film), intent(out) :: next
    real(dp), intent(out) :: gas_change(5),wall_energy
    logical, intent(out) :: ok
    type(spray_wall_film) :: candidate
    real(dp) :: mdot,qgas,qwall,tdot,rho,area,cp
    next=patch;gas_change=0;wall_energy=0;ok=.false.
    if (.not.all(ieee_is_finite([dt,minimum_height]))) return
    if (dt<0.or.dt>1e6_dp.or.minimum_height<1e-15_dp.or.minimum_height>1) return
    call wall_film_transfer_rates(patch,liquid,gas,conductivity,wall_temperature,heat,evaporation,mdot,qgas,qwall,tdot,ok)
    if (.not.ok.or.dt==0) return
    ok=.false.;candidate=patch
    candidate%mass=max(0.0_dp,patch%mass-dt*mdot)
    if (candidate%mass>0) then
      candidate%temperature=patch%temperature+dt*tdot
      call liquid_transfer_properties(liquid,candidate%temperature,gas%pressure,gas%liquid_molecular_weight,rho,cp,ok)
      if(.not.ok) return
      ok=.false.
      area=.25_dp*pi*patch%diameter**2
      candidate%height=candidate%mass/(rho*area)
      if (mdot>0.and.candidate%height<=minimum_height) candidate%mass=0
    end if
    if (candidate%mass==0) then
      candidate%height=0;candidate%temperature=patch%temperature
    else
      if (.not.valid_wall_film(candidate)) return
    end if
    gas_change=wall_film_totals(patch,liquid)-wall_film_totals(candidate,liquid)
    gas_change(5)=gas_change(5)+qwall*dt
    wall_energy=-qwall*dt
    if (.not.all(ieee_is_finite([gas_change,wall_energy]))) then
      gas_change=0;wall_energy=0;return
    end if
    next=candidate;ok=.true.
  end subroutine

  logical function valid_wall_options(o,faces) result(ok)
    type(spray_wall_options), intent(in) :: o
    character(len=*), intent(in) :: faces(6)
    integer :: f
    real(dp) :: constant
    ok=.false.
    if (.not.all(ieee_is_finite([o%temperature,o%contact_angle]))) return
    if (any(o%temperature<1).or.any(o%temperature>1e5_dp)) return
    if (o%contact_angle<1e-6_dp.or.o%contact_angle>pi-1e-6_dp) return
    if (o%maximum_films<1.or.o%maximum_films>1000000) return
    if (.not.all(ieee_is_finite([o%liquid_conductivity,o%minimum_height]))) return
    if (o%minimum_height<1e-15_dp.or.o%minimum_height>1) return
    if (o%evolve) then
      if (.not.any(o%enabled)) return
    end if
    if (.not.valid_conductivity_model(o%conductivity_model,o%liquid_conductivity, &
        o%conductivity_coefficients,o%conductivity_temperature_range,o%evolve)) return
    do f=1,16
      constant=o%port_conductivities(f)
      if (.not.ieee_is_finite(constant)) return
      if (constant/=-1.and.(constant<=0.or.constant>1e4_dp)) return
      if (constant==-1) constant=o%liquid_conductivity
      if (o%port_conductivity_models(f)=='') then
        if (any(o%port_conductivity_coefficients(:,f)/=0).or. &
            any(o%port_conductivity_temperature_ranges(:,f)/=0)) return
        if (o%port_conductivities(f)/=-1.and.o%conductivity_model/='constant') return
      else
        if (.not.valid_conductivity_model(o%port_conductivity_models(f),constant, &
            o%port_conductivity_coefficients(:,f),o%port_conductivity_temperature_ranges(:,f),o%evolve)) return
      end if
    end do
    do f=1,6
      if (.not.o%enabled(f)) cycle
      if (faces(f)/='slip_wall'.and.faces(f)/='no_slip_wall'.and.faces(f)/='isothermal_wall') return
    end do
    ok=.true.
  end function

  logical function valid_wall_component_inventory(patch,masses) result(ok)
    type(spray_wall_film), intent(in) :: patch
    real(dp), intent(in) :: masses(:)
    integer :: n
    ok=.false.;n=size(masses)
    if (n<1.or.n>16.or..not.valid_wall_film(patch)) return
    if (.not.all(ieee_is_finite(masses))) return
    if (any(masses<0).or.any(masses>1e18_dp)) return
    ok=abs(sum(masses)-patch%mass)<=32*epsilon(1.0_dp)*n*patch%mass
  end function

  logical function valid_wall_film(f) result(ok)
    type(spray_wall_film), intent(in) :: f
    ok=.false.
    if (.not.all(ieee_is_finite([f%position,f%normal,f%mass,f%temperature,f%diameter,f%height]))) return
    if (maxval(abs(f%position))>1e9_dp.or.maxval(abs(f%normal))>1) return
    ok=abs(sum(f%normal**2)-1)<1e-12_dp.and.f%mass>0.and.f%mass<=1e18_dp.and. &
      f%temperature>=1.and.f%temperature<=1e5_dp.and.f%diameter>0.and.f%height>0.and. &
      f%material_id>=1.and.f%material_id<=16
  end function

  ! AhamedSplash uses the mass-weighted pure-component boiling temperatures.
  subroutine mixture_impact_boiling_temperature(components,y,molecular_weights,pressure,temperature,ok)
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(in) :: y(:),molecular_weights(:),pressure
    real(dp), intent(out) :: temperature
    logical, intent(out) :: ok
    real(dp) :: denom,boil,total
    logical :: valid
    integer :: i
    temperature=0;ok=.false.
    if (size(molecular_weights)/=size(components).or.size(y)/=size(components)) return
    if (.not.all(ieee_is_finite(molecular_weights)).or..not.ieee_is_finite(pressure)) return
    if (any(molecular_weights<1e-3_dp).or.any(molecular_weights>1e5_dp).or.pressure<1.or.pressure>1e10_dp) return
    ! Validate composition/property records without requiring overlapping fit temperatures.
    if (size(y)<1.or.size(y)>16.or..not.all(ieee_is_finite(y))) return
    if (any(y<0).or.abs(sum(y)-1)>1e-12_dp) return
    total=0
    do i=1,size(y)
      if (.not.valid_spray_liquid(components(i))) return
      if (y(i)==0) cycle
      if (components(i)%boiling_model=='pelemp-cap') then
        call pelemp_boiling_temperature(components(i)%latent_heat,components(i)%reference_temperature, &
          components(i)%normal_boiling_temperature,components(i)%critical_temperature, &
          resolved_liquid_molecular_weight(components(i),molecular_weights(i)),pressure,boil,valid)
        if (.not.valid) return
      else if (components(i)%vapor_pressure_model=='antoine') then
        denom=components(i)%antoine_coefficients(1)-log10(pressure/components(i)%antoine_coefficients(4))
        if (denom<=0) return
        boil=components(i)%antoine_coefficients(2)/denom-components(i)%antoine_coefficients(3)
        if (boil<components(i)%antoine_temperature_range(1).or.boil>components(i)%antoine_temperature_range(2)) return
      else
        denom=1/components(i)%reference_temperature-log(pressure/components(i)%reference_pressure)* &
          universal_gas_constant/(components(i)%latent_heat*resolved_liquid_molecular_weight(components(i),molecular_weights(i)))
        if (denom<=0) return
        boil=1/denom
      end if
      if (.not.ieee_is_finite(boil)) return
      if (boil<1.or.boil>1e5_dp) return
      total=total+y(i)*boil
    end do
    temperature=total;ok=.true.
  end subroutine

  subroutine mixture_impact_at_planar_wall(parent,previous,components,component_mass,properties,molecular_weights, &
      pressure,lower,upper,settings,films,first_id,children,child_mass,deposit,film_mass,wall,hit,ok)
    type(spray_parcel), intent(inout) :: parent
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(inout) :: component_mass(:),film_mass(:)
    real(dp), intent(in) :: previous(3),molecular_weights(:),pressure,lower(3),upper(3)
    type(spray_options), intent(in) :: properties
    type(spray_wall_options), intent(in) :: settings
    type(spray_wall_film), intent(in) :: films(:)
    integer(int64), intent(in) :: first_id
    type(spray_parcel), allocatable, intent(inout) :: children(:)
    real(dp), allocatable, intent(inout) :: child_mass(:,:)
    type(spray_wall_film), intent(inout) :: deposit
    real(dp), intent(out) :: wall(5)
    logical, intent(out) :: hit,ok
    type(spray_options) :: local_properties
    type(spray_liquid) :: effective
    type(spray_parcel) :: next
    type(spray_parcel), allocatable :: spawned(:)
    type(spray_wall_film) :: deposited
    real(dp), allocatable :: spawned_mass(:,:)
    real(dp) :: y(size(components)),rho,cp,e,boil,budget(5)
    logical :: valid,crossed
    integer :: n
    wall=0;hit=.false.;ok=.false.
    if (size(component_mass)/=size(components).or.size(film_mass)/=size(components)) return
    if (.not.valid_component_inventory(parent,component_mass).or.parent%mass<=0) return
    y=component_mass/parent%mass
    call liquid_mixture_properties(components,y,parent%temperature,rho,cp,e,valid,pressure,molecular_weights)
    if (.not.valid) return
    call mixture_impact_boiling_temperature(components,y,molecular_weights,pressure,boil,valid)
    if (.not.valid) return
    effective=components(1);effective%density=rho;effective%heat_capacity=cp
    effective%reference_temperature=parent%temperature;effective%reference_energy=e
    effective%density_model='constant';effective%density_coefficients=0;effective%density_temperature_range=0
    effective%heat_capacity_model='constant';effective%heat_capacity_coefficients=0
    effective%heat_capacity_temperature_range=0;effective%boiling_model='reject';effective%critical_temperature=0
    effective%normal_boiling_temperature=0
    local_properties=properties
    call apply_mixture_breakup_viscosity(components,y,parent%temperature,local_properties,valid,pressure,molecular_weights)
    if(.not.valid) return
    next=parent
    call impact_at_planar_wall(next,previous,effective,local_properties,molecular_weights(1),pressure,lower,upper,settings, &
      films,first_id,spawned,deposited,budget,crossed,valid,boiling_temperature=boil)
    if (.not.valid) return
    if (.not.crossed) then
      ok=.true.;return
    end if
    allocate(spawned_mass(size(components),size(spawned)))
    do n=1,size(spawned)
      spawned_mass(:,n)=spawned(n)%mass*y
    end do
    component_mass=next%mass*y;film_mass=deposited%mass*y
    parent=next;deposit=deposited;wall=budget;hit=.true.;ok=.true.
    call move_alloc(spawned,children);call move_alloc(spawned_mass,child_mass)
  end subroutine

  subroutine impact_at_planar_wall(parent,previous,liquid,properties,molecular_weight,pressure,lower,upper, &
      settings,films,first_id,children,deposit,wall,hit,ok,boiling_temperature)
    type(spray_parcel), intent(inout) :: parent
    real(dp), intent(in) :: previous(3),molecular_weight,pressure,lower(3),upper(3)
    type(spray_liquid), intent(in) :: liquid
    type(spray_options), intent(in) :: properties
    type(spray_wall_options), intent(in) :: settings
    type(spray_wall_film), intent(in) :: films(:)
    integer(int64), intent(in) :: first_id
    type(spray_parcel), allocatable, intent(inout) :: children(:)
    type(spray_wall_film), intent(out) :: deposit
    real(dp), intent(out) :: wall(5)
    real(dp), optional, intent(in) :: boiling_temperature
    logical, intent(out) :: hit,ok
    real(dp) :: fraction,earliest,location(3),normal(3),delta(3),boundary,boil,denom,rho,mu,d,un,we,re,kv,kcrit,thickness,cp
    type(spray_liquid) :: effective
    integer :: a,face,f,regime
    logical :: valid
    wall=0;hit=.false.;ok=.false.;deposit=spray_wall_film()
    if (.not.all(ieee_is_finite([previous,lower,upper,pressure,molecular_weight]))) return
    if (any(upper<=lower).or.pressure<1.or.pressure>1e10_dp.or.molecular_weight<=0) return
    if (any(previous<lower).or.any(previous>=upper)) return
    delta=parent%position-previous;earliest=2;face=0
    do a=1,3
      f=0;boundary=0
      if (parent%position(a)<lower(a)) then
        f=2*a-1;boundary=lower(a)
      else if (parent%position(a)>=upper(a)) then
        f=2*a;boundary=upper(a)
      end if
      if (f==0.or.delta(a)==0) cycle
      fraction=(boundary-previous(a))/delta(a)
      if (fraction<earliest) then
        earliest=fraction;face=f
      end if
    end do
    if (face==0) then
      ok=.true.;return
    end if
    if (.not.settings%enabled(face).or.parent%mass<=0) then
      ok=.true.;return
    end if
    if (.not.valid_spray_liquid(liquid)) return
    location=previous+earliest*delta;normal=0;a=(face+1)/2
    normal(a)=merge(1.0_dp,-1.0_dp,modulo(face,2)==1)
    location(a)=merge(lower(a),upper(a),modulo(face,2)==1)
    if (present(boiling_temperature)) then
      if (.not.ieee_is_finite(boiling_temperature)) return
      boil=boiling_temperature
    else if (liquid%boiling_model=='pelemp-cap') then
      call pelemp_boiling_temperature(liquid%latent_heat,liquid%reference_temperature,liquid%normal_boiling_temperature, &
        liquid%critical_temperature,resolved_liquid_molecular_weight(liquid,molecular_weight),pressure,boil,valid)
      if (.not.valid) return
    else if (liquid%vapor_pressure_model=='antoine') then
      denom=liquid%antoine_coefficients(1)-log10(pressure/liquid%antoine_coefficients(4))
      if (denom<=0) return
      boil=liquid%antoine_coefficients(2)/denom-liquid%antoine_coefficients(3)
      if (boil<liquid%antoine_temperature_range(1).or.boil>liquid%antoine_temperature_range(2)) return
    else
      denom=1/liquid%reference_temperature-log(pressure/liquid%reference_pressure)* &
        universal_gas_constant/(liquid%latent_heat*resolved_liquid_molecular_weight(liquid,molecular_weight))
      if (denom<=0) return
      boil=1/denom
    end if
    if (boil<1.or.boil>1e5_dp) return
    if (.not.ieee_is_finite(properties%surface_tension)) return
    call liquid_transfer_properties(liquid,parent%temperature,pressure,molecular_weight,rho,cp,valid)
    if(.not.valid) return
    mu=liquid_viscosity_at_temperature(properties,parent%temperature)
    d=parcel_diameter(parent,liquid,pressure,molecular_weight)
    if (rho<=0.or.mu<=0.or.d<=0.or.properties%surface_tension<1e-8_dp.or.properties%surface_tension>1e3_dp) return
    un=abs(dot_product(parent%velocity,normal))
    we=rho*d*un**2/properties%surface_tension;re=un*d*rho/mu
    kv=sqrt(we*sqrt(re))
    kcrit=pelec_splash_critical_parameter(pressure,settings%temperature(face))
    thickness=0
    do f=1,size(films)
      if (maxval(abs(films(f)%normal-normal))>1e-12_dp) cycle
      delta=location-films(f)%position
      if (abs(dot_product(delta,normal))>1e-10_dp*maxval(upper-lower)) cycle
      if (norm2(delta)<=.5_dp*films(f)%diameter) thickness=thickness+films(f)%height
    end do
    effective=liquid;effective%density=rho;effective%heat_capacity=cp
    effective%reference_energy=liquid_energy(liquid,parent%temperature);effective%reference_temperature=parent%temperature
    effective%density_model='constant';effective%density_coefficients=0;effective%density_temperature_range=0
    effective%heat_capacity_model='constant';effective%heat_capacity_coefficients=0
    effective%heat_capacity_temperature_range=0;effective%boiling_model='reject';effective%critical_temperature=0
    effective%normal_boiling_temperature=0
    call resolve_wall_impact(parent,effective,location,normal,kv,settings%temperature(face)/boil,kcrit,thickness/d, &
      we,re,settings%contact_angle,first_id,children,deposit,wall,regime,ok, &
      reference_parent_position=properties%mixture_time_scheme=='euler')
    hit=ok
  end subroutine

  pure function cross(a,b) result(c)
    real(dp), intent(in) :: a(3),b(3)
    real(dp) :: c(3)
    c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
  end function

  ! Component masses are per physical drop for parcels and total kg for the film.
  ! Dimensionless impact parameters are supplied by the caller, as in the pure kernel.
  subroutine resolve_mixture_wall_impact(parent,components,component_mass,location,normal,kv,tstar,kcrit,film_ratio, &
      weber,reynolds,contact_angle,first_id,children,child_mass,film,film_mass,wall,regime,ok,pressure,gas_molecular_weights, &
      reference_parent_position)
    type(spray_parcel), intent(inout) :: parent
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(inout) :: component_mass(:),film_mass(:)
    real(dp), intent(in) :: location(3),normal(3),kv,tstar,kcrit,film_ratio,weber,reynolds,contact_angle
    integer(int64), intent(in) :: first_id
    type(spray_parcel), allocatable, intent(inout) :: children(:)
    real(dp), allocatable, intent(inout) :: child_mass(:,:)
    type(spray_wall_film), intent(inout) :: film
    real(dp), intent(out) :: wall(5)
    integer, intent(out) :: regime
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: pressure,gas_molecular_weights(:)
    logical, optional, intent(in) :: reference_parent_position
    type(spray_liquid) :: effective
    type(spray_parcel) :: candidate
    type(spray_parcel), allocatable :: spawned(:)
    type(spray_wall_film) :: deposited
    real(dp), allocatable :: spawned_mass(:,:)
    real(dp) :: fractions(size(components)),rho,cp,e,budget(5)
    integer :: n,kind
    logical :: valid
    wall=0;regime=-1;ok=.false.
    if (size(components)/=size(component_mass).or.size(components)/=size(film_mass)) return
    if (.not.valid_component_inventory(parent,component_mass).or.parent%mass<=0) return
    fractions=component_mass/parent%mass
    call liquid_mixture_properties(components,fractions,parent%temperature,rho,cp,e,valid,pressure,gas_molecular_weights)
    if (.not.valid) return
    effective=components(1);effective%density=rho;effective%heat_capacity=cp
    effective%reference_temperature=parent%temperature;effective%reference_energy=e
    effective%density_model='constant';effective%density_coefficients=0;effective%density_temperature_range=0
    effective%heat_capacity_model='constant';effective%heat_capacity_coefficients=0
    effective%heat_capacity_temperature_range=0;effective%boiling_model='reject';effective%critical_temperature=0
    effective%normal_boiling_temperature=0
    candidate=parent
    call resolve_wall_impact(candidate,effective,location,normal,kv,tstar,kcrit,film_ratio,weber,reynolds, &
      contact_angle,first_id,spawned,deposited,budget,kind,valid,reference_parent_position)
    if (.not.valid) return
    allocate(spawned_mass(size(components),size(spawned)))
    do n=1,size(spawned)
      spawned_mass(:,n)=spawned(n)%mass*fractions
    end do
    component_mass=candidate%mass*fractions;film_mass=deposited%mass*fractions
    parent=candidate;film=deposited;wall=budget;regime=kind
    call move_alloc(spawned,children);call move_alloc(spawned_mass,child_mass)
    ok=.true.
  end subroutine

  subroutine resolve_wall_impact(parent,liquid,location,normal,kv,tstar,kcrit,film_ratio,weber,reynolds, &
      contact_angle,first_id,children,film,wall,regime,ok,reference_parent_position)
    type(spray_parcel), intent(inout) :: parent
    type(spray_liquid), intent(in) :: liquid
    real(dp), intent(in) :: location(3),normal(3),kv,tstar,kcrit,film_ratio,weber,reynolds,contact_angle
    integer(int64), intent(in) :: first_id
    type(spray_parcel), allocatable, intent(inout) :: children(:)
    type(spray_wall_film), intent(inout) :: film
    real(dp), intent(out) :: wall(5)
    integer, intent(out) :: regime
    logical, intent(out) :: ok
    logical, optional, intent(in) :: reference_parent_position
    type(spray_parcel) :: next
    type(spray_parcel), allocatable :: spawned(:)
    type(spray_wall_film) :: deposited
    real(dp) :: un,ut,speed,alpha,ms,dm,restitution,remaining,psi(3),beta(3),helper(3),budget(5),total(5)
    integer :: kind,i,axis
    logical :: valid
    wall=0;regime=-1;ok=.false.
    if (.not.valid_spray_parcel(parent).or..not.valid_spray_liquid(liquid)) return
    if (.not.all(ieee_is_finite([location,normal,weber,reynolds,contact_angle]))) return
    if (maxval(abs(location))>1e9_dp.or.maxval(abs(normal))>1) return
    if (abs(sum(normal**2)-1)>1e-12_dp.or.parent%mass<=0) return
    if (liquid_heat_capacity_at_temperature(liquid,parent%temperature)<=0.or. &
        liquid_density_at_temperature(liquid,parent%temperature)<=0) return
    un=dot_product(parent%velocity,normal);speed=norm2(parent%velocity)
    if (un>=0.or.speed<=0) return
    ut=norm2(parent%velocity-un*normal);alpha=max(pi/6,asin(min(1.0_dp,abs(un)/speed)))
    call classify_wall_impact(kv,tstar,alpha,kcrit,film_ratio,kind,ms,dm,valid)
    if (.not.valid) return
    ! The reference wet branch leaves dm unset. Close its incident-mass budget
    ! by depositing the unsplashed fraction; no entrainment from existing film.
    if (kind==impact_wet) dm=1-ms
    next=parent;allocate(spawned(0))
    if (ms>0) then
      call construct_splash_children(parent,liquid,location,normal,kv,film_ratio,ms,first_id,spawned,budget,valid)
      if (.not.valid) return
    end if
    if (dm>0) then
      call construct_wall_film(parent,liquid,location,normal,dm,weber,reynolds,contact_angle,deposited,budget,valid)
      if (.not.valid) return
    end if
    remaining=max(0.0_dp,1-ms-dm);next%mass=parent%mass*remaining
    next%position=location
    if (remaining>0) then
      select case(kind)
      case(impact_rebound)
        restitution=.823_dp*exp(-.0019835_dp*kv**1.6_dp)
        next%velocity=parent%velocity-(1+restitution)*un*normal
      case(impact_thermal)
        psi=cross(-parent%velocity,normal)
        if (norm2(psi)<1e-5_dp) psi=cross(-parent%velocity+1.0_dp,normal)
        if (norm2(psi)<1e-12_dp) then
          axis=minloc(abs(normal),dim=1);helper=0;helper(axis)=1;psi=cross(helper,normal)
        end if
        psi=psi/norm2(psi);beta=cross(psi,normal);beta=beta/norm2(beta)
        next%velocity=(1.5_dp*sqrt(3.0_dp)*ut/pi)*beta+ &
          (.0028_dp*kv*exp(-.0062_dp*kv)*speed)*normal
        next%deformation=0;next%deformation_rate=0;next%rt_time=-1
        next%khrt_initial_count=0;next%khrt_shed_volume=0
      case default
        return
      end select
      next%position=location+parcel_diameter(next,liquid)*normal
      if(present(reference_parent_position)) then
        if(reference_parent_position) then
          if(kind==impact_rebound) then
            next%position=parent%position-(1+restitution)*dot_product(parent%position-location,normal)*normal
          else if(kind==impact_thermal) then
            next%position=location
          endif
        endif
      endif
    end if
    if (next%mass==0) then
      next%khrt_initial_count=0;next%khrt_shed_volume=0
    end if
    if (.not.valid_spray_parcel(next)) return
    total=parcel_totals(next,liquid)
    do i=1,size(spawned)
      total=total+parcel_totals(spawned(i),liquid)
    end do
    if (deposited%mass>0) total=total+wall_film_totals(deposited,liquid)
    budget=parcel_totals(parent,liquid)-total
    if (.not.all(ieee_is_finite(budget))) return
    if (abs(budget(1))>1e-12_dp*parent%mass*parent%multiplicity) return
    budget(1)=0
    parent=next;film=deposited;call move_alloc(spawned,children)
    wall=budget;regime=kind;ok=.true.
  end subroutine

  subroutine construct_splash_children(parent,liquid,location,normal,kv,film_ratio,fraction,first_id,children,wall,ok)
    type(spray_parcel), intent(in) :: parent
    type(spray_liquid), intent(in) :: liquid
    real(dp), intent(in) :: location(3),normal(3),kv,film_ratio,fraction
    integer(int64), intent(in) :: first_id
    type(spray_parcel), allocatable, intent(inout) :: children(:)
    real(dp), intent(out) :: wall(5)
    logical, intent(out) :: ok
    type(spray_parcel), allocatable :: candidate(:)
    real(dp) :: un,ut,speed,alpha,psi(3),beta(3),helper(3),masses(4),v(5),ub,up,angle,total(5)
    integer :: i,axis
    wall=0;ok=.false.
    if (.not.valid_spray_parcel(parent).or..not.valid_spray_liquid(liquid)) return
    if (.not.all(ieee_is_finite([location,normal,fraction]))) return
    if (maxval(abs(location))>1e9_dp.or.maxval(abs(normal))>1.0_dp) return
    if (abs(sum(normal**2)-1)>1e-12_dp.or.fraction<=0.or.fraction>1.or.parent%mass<=0) return
    if (first_id<=parent%id.or.first_id>huge(first_id)-4_int64) return
    if (liquid_heat_capacity_at_temperature(liquid,parent%temperature)<=0) return
    if (liquid_density_at_temperature(liquid,parent%temperature)<=0) return
    un=dot_product(parent%velocity,normal);speed=norm2(parent%velocity)
    if (un>=0.or.speed<=0) return
    ut=norm2(parent%velocity-un*normal)
    alpha=max(pi/6,asin(min(1.0_dp,abs(un)/speed)))
    psi=cross(-parent%velocity,normal)
    if (norm2(psi)<1e-5_dp) then
      ! Reference tangent fallback; use a coordinate axis if its [1,1,1] offset is parallel too.
      psi=cross(-parent%velocity+1.0_dp,normal)
      if (norm2(psi)<1e-12_dp) then
        axis=minloc(abs(normal),dim=1);helper=0;helper(axis)=1
        psi=cross(helper,normal)
      end if
    end if
    psi=psi/norm2(psi);beta=cross(psi,normal);beta=beta/norm2(beta)
    call splash_velocity_components(un,ut,kv,film_ratio,v,ok)
    if (.not.ok) return
    ok=.false.
    if (v(5)<0) return
    call splash_mass_fractions(alpha,fraction,film_ratio,masses,ok)
    if (.not.ok) return
    ok=.false.;allocate(candidate(4));total=0
    do i=1,4
      candidate(i)=parent
      candidate(i)%id=first_id+int(i-1,int64);candidate(i)%parent_id=parent%id
      candidate(i)%mass=parent%mass*masses(i)
      candidate(i)%deformation=0;candidate(i)%deformation_rate=0;candidate(i)%rt_time=-1
      ! Pinned SpraySB.cpp initializes a splash drop's N0 to its own
      ! number density, even when the incident parent used a different N0.
      candidate(i)%khrt_initial_count=candidate(i)%multiplicity;candidate(i)%khrt_shed_volume=0
      ub=v(2)
      if (i==2) ub=v(1)
      if (i==4) ub=v(3)
      angle=.5_dp*pi*real(i-2,dp);up=v(4)*sin(angle)
      candidate(i)%velocity=v(5)*normal+ub*beta+up*psi
      candidate(i)%position=location+parcel_diameter(candidate(i),liquid)*normal
      if (.not.valid_spray_parcel(candidate(i))) return
      total=total+parcel_totals(candidate(i),liquid)
    end do
    wall=fraction*parcel_totals(parent,liquid)-total
    ! Only impulse/work goes to the wall in this operation; mass stays in children.
    if (abs(wall(1))>1e-12_dp*parent%mass*parent%multiplicity) then
      wall=0;return
    end if
    wall(1)=0
    if (.not.all(ieee_is_finite(wall))) then
      wall=0;return
    end if
    call move_alloc(candidate,children);ok=.true.
  end subroutine

  pure function wall_film_totals(film,liquid) result(totals)
    type(spray_wall_film), intent(in) :: film
    type(spray_liquid), intent(in) :: liquid
    real(dp) :: totals(5)
    totals=0;totals(1)=film%mass
    totals(5)=film%mass*liquid_energy(liquid,film%temperature)
  end function

  subroutine construct_wall_film(parent,liquid,location,normal,fraction,weber,reynolds,contact_angle,film,wall,ok)
    type(spray_parcel), intent(in) :: parent
    type(spray_liquid), intent(in) :: liquid
    real(dp), intent(in) :: location(3),normal(3),fraction,weber,reynolds,contact_angle
    type(spray_wall_film), intent(inout) :: film
    real(dp), intent(out) :: wall(5)
    logical, intent(out) :: ok
    type(spray_wall_film) :: candidate
    real(dp) :: rho,volume,combined_diameter,combined_we,combined_re,denominator,drop_diameter
    wall=0;ok=.false.
    if (.not.valid_spray_parcel(parent).or..not.valid_spray_liquid(liquid)) return
    if (.not.all(ieee_is_finite([location,normal,fraction,weber,reynolds,contact_angle]))) return
    if (maxval(abs(location))>1e9_dp.or.maxval(abs(normal))>1) return
    if (abs(sum(normal**2)-1)>1e-12_dp.or.parent%mass<=0.or.fraction<=0.or.fraction>1) return
    if (weber<0.or.weber>1e12_dp.or.reynolds<=0.or.reynolds>1e12_dp) return
    if (contact_angle<1e-6_dp.or.contact_angle>pi-1e-6_dp) return
    rho=liquid_density_at_temperature(liquid,parent%temperature)
    if (rho<=0.or.liquid_heat_capacity_at_temperature(liquid,parent%temperature)<=0) return
    candidate%mass=fraction*parent%mass*parent%multiplicity;volume=candidate%mass/rho
    candidate%position=location;candidate%normal=normal;candidate%temperature=parent%temperature
    candidate%material_id=parent%material_id
    combined_diameter=(6*volume/pi)**(1.0_dp/3)
    drop_diameter=parcel_diameter(parent,liquid)
    combined_we=weber*combined_diameter/drop_diameter;combined_re=reynolds*combined_diameter/drop_diameter
    ! 6(1-cos(theta))/sin(theta)^2 = 6/(1+cos(theta)).
    denominator=6/(1+cos(contact_angle))+4*combined_we/sqrt(combined_re)
    candidate%diameter=combined_diameter*sqrt((12+combined_we)/denominator)
    candidate%height=volume/(.25_dp*pi*candidate%diameter**2)
    if (.not.all(ieee_is_finite([candidate%diameter,candidate%height]))) return
    if (candidate%diameter<=0.or.candidate%height<=0) return
    wall=fraction*parcel_totals(parent,liquid)-wall_film_totals(candidate,liquid);wall(1)=0
    if (.not.all(ieee_is_finite(wall))) then
      wall=0;return
    end if
    film=candidate;ok=.true.
  end subroutine

  subroutine classify_wall_impact(kv,tstar,alpha,kcrit,film_ratio,regime,splashed,deposited,ok)
    real(dp), intent(in) :: kv,tstar,alpha,kcrit,film_ratio
    integer, intent(out) :: regime
    real(dp), intent(out) :: splashed,deposited
    logical, intent(out) :: ok
    real(dp) :: wet_critical,a,b,base,angle
    regime=-1;splashed=0;deposited=0;ok=.false.
    if (.not.all(ieee_is_finite([kv,tstar,alpha,kcrit,film_ratio]))) return
    if (kv<0.or.kv>1e8_dp.or.tstar<=0.or.tstar>1e5_dp.or.abs(kcrit)>1e8_dp) return
    if (film_ratio<0.or.film_ratio>1e5_dp.or.alpha<pi/6.or.alpha>pi/2) return
    if (film_ratio<.01_dp) then
      if (tstar<1.1_dp) then
        if (kv<kcrit) then
          regime=impact_deposit;deposited=1
        else
          regime=impact_dry
          splashed=.75_dp*(1-exp(-.003_dp*(kv-kcrit)));deposited=1-splashed
        end if
      else
        if (kv<kcrit) then
          regime=impact_rebound
        else
          regime=impact_thermal
          splashed=.75_dp*(1-exp(-.0175_dp*(kv-kcrit)))
        end if
      end if
    else
      wet_critical=kcrit+226*film_ratio**.72_dp
      if (kv>wet_critical) then
        regime=impact_wet;a=4.5_dp;b=.2_dp
        if (film_ratio<.5_dp) then
          a=10.42_dp*film_ratio+.32_dp;b=8.75_dp*film_ratio-.107_dp
        end if
        base=a*(1-exp(-.003_dp*(kv-wet_critical)))+b
        angle=alpha*180/pi
        if (film_ratio>.1_dp) then
          if (angle>45) then
            splashed=base*(-.0305_dp*angle+3.7496_dp)
          else
            splashed=base*.0525_dp*angle
          end if
        else
          if (angle>45) then
            splashed=base*(-.0072_dp*angle+1.6389_dp)
          else
            splashed=base*(-.00046_dp*angle**2+.05_dp*angle)
          end if
        end if
        splashed=min(splashed,1.0_dp)
      else
        regime=impact_spread;deposited=1
      end if
    end if
    ! Preserve correlations, but do not pass a negative mass fraction to a caller.
    ok=splashed>=0.and.splashed<=1.and.deposited>=0.and.splashed+deposited<=1
    if (.not.ok) then
      regime=-1;splashed=0;deposited=0
    end if
  end subroutine

  subroutine splash_velocity_components(normal_speed,tangent_speed,kv,film_ratio,velocity,ok)
    real(dp), intent(in) :: normal_speed,tangent_speed,kv,film_ratio
    ! beta(0), beta(pi/2), beta(pi), psi coefficient, outward normal speed.
    real(dp), intent(out) :: velocity(5)
    logical, intent(out) :: ok
    real(dp) :: un,ut
    velocity=0;ok=.false.
    if (.not.all(ieee_is_finite([normal_speed,tangent_speed,kv,film_ratio]))) return
    if (abs(normal_speed)>1e7_dp.or.tangent_speed<0.or.tangent_speed>1e7_dp) return
    if (kv<0.or.kv>1e8_dp.or.film_ratio<0.or.film_ratio>1e5_dp) return
    un=abs(normal_speed)
    if (film_ratio<.01_dp) then
      ut=.0065_dp*kv*exp(-.004_dp*kv)*un
      velocity=[tangent_speed+ut,cos(pi/6)*tangent_speed,cos(pi/3)*tangent_speed-ut, &
        .2_dp*tangent_speed+ut,.0028_dp*kv*exp(-.0062_dp*kv)*un]
    else if (film_ratio<.1_dp) then
      if (film_ratio>.05_dp) then
        ut=(5e-5_dp*kv+.0203_dp)*un
        velocity(5)=.03_dp*tangent_speed+(1.32e-4_dp*kv+.0158_dp)*un
      else
        ut=(3.45e-5_dp*kv+.006_dp)*un
        velocity(5)=.03_dp*tangent_speed+(3.57e-5_dp*kv-.0057_dp)*un
      end if
      velocity(1:4)=[.3_dp*tangent_speed+ut,.1_dp*tangent_speed, &
        .033_dp*tangent_speed-ut,.06_dp*tangent_speed+ut]
    else
      velocity=[.1_dp*tangent_speed,.05_dp*tangent_speed,.01_dp*tangent_speed, &
        .02_dp*tangent_speed,.03_dp*tangent_speed]
    end if
    ok=all(ieee_is_finite(velocity))
  end subroutine

  subroutine splash_mass_fractions(alpha,splashed,film_ratio,fractions,ok)
    real(dp), intent(in) :: alpha,splashed,film_ratio
    ! Azimuths -pi/2, 0, pi/2, pi; fractions of the original incident mass.
    real(dp), intent(out) :: fractions(4)
    logical, intent(out) :: ok
    real(dp) :: angle,forward,side,tangent
    fractions=0;ok=.false.
    if (.not.all(ieee_is_finite([alpha,splashed,film_ratio]))) return
    if (alpha<pi/6.or.alpha>pi/2.or.splashed<0.or.splashed>1.or.film_ratio<0.or.film_ratio>1e5_dp) return
    if (film_ratio>0) then
      angle=alpha*180/pi
      if (film_ratio>.1_dp) then
        forward=(1.9e-4_dp*angle**2+.0202_dp*angle-.0107_dp)*splashed
        side=(6.4e-5_dp*angle**2-.0076_dp*angle+.4176_dp)*splashed
      else
        forward=(-2.3e-4_dp*angle**2+.0236_dp*angle)*splashed
        side=(6e-5_dp*angle**2-.0068_dp*angle+.3619_dp)*splashed
      end if
      fractions=[side,forward,side,splashed-forward-2*side]
    else
      tangent=.4_dp*tan(pi/2-alpha)
      fractions=.25_dp*splashed*[1.0_dp,1+tangent,1.0_dp,1-tangent]
    end if
    ok=all(fractions>=0).and.abs(sum(fractions)-splashed)<=32*epsilon(1.0_dp)
    if (.not.ok) fractions=0
  end subroutine
end module
