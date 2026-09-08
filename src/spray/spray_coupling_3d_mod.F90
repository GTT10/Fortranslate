module spray_coupling_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species, nasa7_mass_properties
  use gas_transport_mod, only: gas_transport_species, compatible_transport_database
  use mixture_thermo_mod, only: mixture_mass_properties
  use mixture_transport_mod, only: mixture_transport_coefficients, binary_diffusion_coefficients
  use reactive_1d_mod, only: reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use spray_parcel_mod
  implicit none
  private
  public :: periodic_particle_stencil, build_spray_film, advance_spray_coupling_3d
  public :: spray_system_totals
contains
  subroutine periodic_particle_stencil(position, lower, spacing, cells, indices, weights, ok)
    real(dp), intent(in) :: position(3), lower(3), spacing(3)
    integer, intent(in) :: cells(3)
    integer, intent(out) :: indices(3,8)
    real(dp), intent(out) :: weights(8)
    logical, intent(out) :: ok
    real(dp) :: coordinate(3), fraction(3), factors(3)
    integer :: base(3), offset(3), a,b,c,n
    indices=0; weights=0.0_dp; ok=.false.
    if (any(cells < 2)) return
    if (.not. all(ieee_is_finite([position,lower,spacing]))) return
    if (any(spacing < 1.0e-12_dp) .or. any(spacing > 1.0e9_dp)) return
    if (maxval(abs(position)) > 1.0e9_dp .or. maxval(abs(lower)) > 1.0e9_dp) return
    coordinate=modulo(position-lower,spacing*real(cells,dp))/spacing-0.5_dp
    base=floor(coordinate)
    fraction=coordinate-real(base,dp)
    n=0
    do c=0,1
      do b=0,1
        do a=0,1
          n=n+1; offset=[a,b,c]
          indices(:,n)=modulo(base+offset,cells)+1
          factors=merge(fraction,1.0_dp-fraction,offset==1)
          weights(n)=product(factors)
        end do
      end do
    end do
    ok=abs(sum(weights)-1.0_dp) < 1.0e-14_dp .and. all(weights >= 0.0_dp)
  end subroutine

  subroutine build_spray_film(species, transport, primitive, gas_temperature, parcel, liquid, options, film, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: primitive(:), gas_temperature
    type(spray_parcel), intent(in) :: parcel
    type(spray_liquid), intent(in) :: liquid
    type(spray_options), intent(in) :: options
    type(spray_film), intent(out) :: film
    logical, intent(out) :: ok
    real(dp) :: y(size(species)), yf(size(species)), diffusion(size(species)), tf
    real(dp) :: binary(size(species),size(species)),carrier_sum,resistance
    real(dp) :: xs,ys,carrier_weight,carrier_inverse,cv,h,e,s,mw,r,cp,gamma,species_cp
    integer :: v,k
    film=spray_film(); ok=.false.
    v=liquid%vapor_index
    if (.not. valid_spray_liquid(liquid) .or. .not. valid_spray_parcel(parcel)) return
    if (v > size(species) .or. v < 1) return
    if (trim(species(v)%name) /= trim(liquid%vapor_name)) return
    if (size(primitive) /= reactive_nprim(size(species))) return
    if (.not. all(ieee_is_finite(primitive))) return
    if (.not. ieee_is_finite(gas_temperature)) return
    if (gas_temperature < 1.0_dp .or. gas_temperature > 1.0e5_dp) return
    if (any(primitive(6:) < 0.0_dp)) return
    y=primitive(6:)
    if (abs(sum(y)-1.0_dp) > 5.0e-12_dp) return
    ys=y(v)
    if (options%evaporation .and. y(v) < 1.0_dp-1.0e-12_dp) then
      call saturation_mole_fraction(liquid,species(v)%molecular_weight,parcel%temperature,primitive(5),xs,ok)
      if (.not. ok) return
      carrier_inverse=0.0_dp
      do k=1,size(species)
        if (k == v) cycle
        carrier_inverse=carrier_inverse+y(k)/species(k)%molecular_weight
      end do
      carrier_weight=(1.0_dp-y(v))/carrier_inverse
      ys=xs*species(v)%molecular_weight/(xs*species(v)%molecular_weight+(1.0_dp-xs)*carrier_weight)
    end if
    yf=y
    if (ys > y(v)) then
      yf=y*(1.0_dp-(2.0_dp*ys+y(v))/3.0_dp)/(1.0_dp-y(v))
      yf(v)=(2.0_dp*ys+y(v))/3.0_dp
    end if
    tf=(2.0_dp*parcel%temperature+gas_temperature)/3.0_dp
    call mixture_mass_properties(species,yf,tf,mw,r,cp,cv,gamma,h,e,s,ok)
    if (.not. ok) return
    film%density=primitive(5)/(r*tf)
    film%cp=cp
    call mixture_transport_coefficients(species,transport,yf,tf,primitive(5), &
      film%viscosity,film%conductivity,diffusion,ok)
    if (.not. ok) return
    ! Direct carrier-mole harmonic average: sum(Xj)/sum(Xj/Dvj), j /= v.
    ! Our mixture D convention differs from PelePhysics' rho-D convention;
    ! multiplying it by Wmix/Wv fails even in a two-species mixture.
    call binary_diffusion_coefficients(species,transport,tf,primitive(5),binary,ok)
    if (.not. ok) return
    carrier_sum=0.0_dp; resistance=0.0_dp
    do k=1,size(species)
      if (k==v) cycle
      carrier_sum=carrier_sum+yf(k)/species(k)%molecular_weight
      resistance=resistance+yf(k)/(species(k)%molecular_weight*binary(v,k))
    end do
    film%diffusivity=diffusion(v)
    if (carrier_sum>0.0_dp .and. resistance>0.0_dp) film%diffusivity=carrier_sum/resistance
    ! At pure vapor the mass flux is zero; diffusivity has no mass-transfer role.
    call nasa7_mass_properties(species(v),parcel%temperature,species_cp,cv,h,e,s,ok)
    if (.not. ok) return
    film%vapor_enthalpy=h
    film%gas_velocity=primitive(2:4)
    film%gas_temperature=gas_temperature
    film%ambient_fraction=y(v)
    film%surface_fraction=ys
    ! In a pure-vapor atmosphere the no-condensation model has zero mass flux.
    if (ys >= 1.0_dp) film%surface_fraction=0.0_dp
    ok=.true.
  end subroutine

  subroutine advance_spray_coupling_3d(species, transport, liquid, options, state, temperature, &
      lower, spacing, parcels, interval, substeps, ok, message)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquid
    type(spray_options), intent(in) :: options
    real(dp), intent(inout) :: state(:,:,:,:), temperature(:,:,:)
    real(dp), intent(in) :: lower(3), spacing(3), interval
    type(spray_parcel), intent(inout) :: parcels(:)
    integer, intent(out) :: substeps
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    real(dp), allocatable :: candidate(:,:,:,:), temp(:,:,:)
    type(spray_parcel), allocatable :: drops(:)
    type(spray_parcel) :: next
    type(spray_film) :: film
    real(dp) :: weights(8), old_cells(size(state,1),8), old_t(8), sampled(size(state,1))
    real(dp) :: primitive(reactive_nprim(size(species))), tg, sound, guess, dt, elapsed, volume
    real(dp) :: invtau,heat,masscoef,d,mdot,tdot,speed,transfer(5),increase(size(state,1)), recovered
    real(dp) :: sampled_pressure,cp,cv,h,e,entropy,xsat
    integer :: shape3(3),indices(3,8),p,n,i,j,k,attempt,total_steps
    logical :: local_ok, accepted
    substeps=0; ok=.false.; message='Invalid spray coupling input'
    if (size(state,1) /= reactive_nvar(size(species)) .or. size(species) < 2) return
    shape3=shape(temperature)
    if (any(shape3 < 2)) return
    if (any(shape(state(1,:,:,:)) /= shape3)) return
    if (.not. valid_spray_liquid(liquid) .or. .not. valid_spray_options(options)) return
    if (liquid%vapor_index > size(species)) return
    if (.not. compatible_transport_database(species,transport)) return
    if (.not. all(ieee_is_finite([lower,spacing,interval]))) return
    if (any(spacing < 1.0e-12_dp) .or. any(spacing > 1.0e9_dp)) return
    if (maxval(abs(lower)) > 1.0e9_dp .or. interval < 0.0_dp .or. interval > 1.0e6_dp) return
    do p=1,size(parcels)
      if (.not. valid_spray_parcel(parcels(p))) return
      if (p > 1) then
        if (parcels(p)%id <= parcels(p-1)%id) then
          message='Parcel IDs must be strictly increasing'; return
        end if
      end if
    end do
    allocate(candidate,source=state); allocate(temp,mold=temperature); allocate(drops,source=parcels)
    call recover_reactive_temperatures_3d(species,state,temperature,shape3(1),shape3(2),shape3(3),temp,local_ok)
    if (.not. local_ok) return
    if (interval == 0.0_dp) then
      ok=.true.; message=''; return
    end if
    volume=product(spacing); total_steps=0
    do p=1,size(drops)
      elapsed=0.0_dp
      do while (elapsed < interval .and. drops(p)%mass > 0.0_dp)
        if (total_steps >= options%maximum_substeps) then
          message='Spray substep limit exceeded; complete update rolled back'; return
        end if
        call periodic_particle_stencil(drops(p)%position,lower,spacing,shape3,indices,weights,local_ok)
        if (.not. local_ok) return
        sampled=0.0_dp; guess=0.0_dp
        do n=1,8
          i=indices(1,n); j=indices(2,n); k=indices(3,n)
          old_cells(:,n)=candidate(:,i,j,k); old_t(n)=temp(i,j,k)
          sampled=sampled+weights(n)*old_cells(:,n)
          guess=guess+weights(n)*old_t(n)
        end do
        call reactive_conserved_to_primitive(species,sampled,guess,primitive,tg,sound,local_ok)
        if (.not. local_ok) return
        call build_spray_film(species,transport,primitive,tg,drops(p),liquid,options,film,local_ok)
        if (.not. local_ok) then
          message='Unsupported film state (including boiling) or invalid liquid/vapor data'; return
        end if
        sampled_pressure=primitive(5)
        call spray_transfer_rates(drops(p),liquid,film,options,invtau,heat,masscoef,local_ok)
        if (.not. local_ok) return
        d=parcel_diameter(drops(p),liquid); mdot=masscoef*d
        tdot=(heat*(tg-drops(p)%temperature)-mdot*(film%vapor_enthalpy-liquid_energy(liquid,drops(p)%temperature))) &
          /(drops(p)%mass*liquid%heat_capacity)
        dt=interval-elapsed
        if (invtau > 0.0_dp) dt=min(dt,options%max_fractional_change/invtau)
        if (mdot > 0.0_dp .and. d > 2.0_dp*options%minimum_diameter) &
          dt=min(dt,options%max_fractional_change*drops(p)%mass/mdot)
        if (abs(tdot) > 0.0_dp) &
          dt=min(dt,min(2.0_dp,options%max_fractional_change*drops(p)%temperature)/abs(tdot))
        speed=max(maxval(abs(drops(p)%velocity)),maxval(abs(film%gas_velocity)))
        if (speed > 0.0_dp) dt=min(dt,options%particle_cfl*minval(spacing)/speed)
        accepted=.false.
        do attempt=1,30
          if (dt <= 0.0_dp .or. elapsed+dt <= elapsed) exit
          call advance_parcel_frozen_film(drops(p),liquid,film,options,dt,next,transfer,local_ok)
          if (local_ok .and. next%mass>0.0_dp) then
            call nasa7_mass_properties(species(liquid%vapor_index),next%temperature,cp,cv,h,e,entropy,local_ok)
            if (local_ok .and. options%evaporation) then
              call saturation_mole_fraction(liquid,species(liquid%vapor_index)%molecular_weight, &
                next%temperature,sampled_pressure,xsat,local_ok)
            end if
          end if
          if (local_ok) then
            increase=0.0_dp; increase(1:5)=transfer; increase(5+liquid%vapor_index)=transfer(1)
            do n=1,8
              i=indices(1,n); j=indices(2,n); k=indices(3,n)
              candidate(:,i,j,k)=old_cells(:,n)+weights(n)*increase/volume
              call reactive_conserved_to_primitive(species,candidate(:,i,j,k),old_t(n),primitive,recovered,sound,local_ok)
              if (.not. local_ok) exit
              temp(i,j,k)=recovered
            end do
          end if
          if (local_ok) then
            accepted=.true.; exit
          end if
          do n=1,8
            i=indices(1,n); j=indices(2,n); k=indices(3,n)
            candidate(:,i,j,k)=old_cells(:,n); temp(i,j,k)=old_t(n)
          end do
          dt=0.5_dp*dt
        end do
        if (.not. accepted) then
          message='Spray feedback cannot preserve gas/liquid admissibility; complete update rolled back'; return
        end if
        next%position=lower+modulo(next%position-lower,spacing*real(shape3,dp))
        drops(p)=next
        elapsed=min(interval,elapsed+dt); total_steps=total_steps+1
      end do
    end do
    state=candidate; temperature=temp; parcels=drops
    substeps=total_steps; ok=.true.; message=''
  end subroutine

  subroutine spray_system_totals(state, spacing, parcels, liquid, totals)
    real(dp), intent(in) :: state(:,:,:,:), spacing(3)
    type(spray_parcel), intent(in) :: parcels(:)
    type(spray_liquid), intent(in) :: liquid
    real(dp), intent(out) :: totals(size(state,1))
    real(dp) :: drop(5)
    integer :: v,p
    do v=1,size(state,1)
      totals(v)=sum(state(v,:,:,:))*product(spacing)
    end do
    do p=1,size(parcels)
      drop=parcel_totals(parcels(p),liquid)
      totals(1:5)=totals(1:5)+drop
      totals(5+liquid%vapor_index)=totals(5+liquid%vapor_index)+drop(1)
    end do
  end subroutine
end module spray_coupling_3d_mod
