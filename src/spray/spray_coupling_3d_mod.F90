module spray_coupling_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use spray_wall_impact_mod
  use nasa7_thermo_mod, only: nasa7_species, nasa7_mass_properties
  use gas_transport_mod, only: gas_transport_species, compatible_transport_database
  use mixture_thermo_mod, only: mixture_mass_properties
  use mixture_transport_mod, only: mixture_transport_coefficients, binary_diffusion_coefficients
  use reactive_1d_mod, only: reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use spray_parcel_mod
  use spray_mixture_mod
  use spray_breakup_mod
  implicit none
  private
  public :: advance_mixture_film_coupling_3d
  public :: periodic_particle_stencil, build_spray_film, advance_spray_coupling_3d
  public :: spray_system_totals, valid_spray_boundary_faces, apply_parcel_boundaries
  public :: particle_boundary_stencil, advance_mixture_coupling_3d, mixture_system_totals
contains
  pure logical function valid_spray_boundary_faces(faces) result(ok)
    character(len=*), intent(in) :: faces(6)
    integer :: a
    ok=.false.
    do a=1,6
      select case(trim(faces(a)))
      case('periodic','outflow','inflow','slip_wall','no_slip_wall','isothermal_wall')
      case default
        return
      end select
    end do
    do a=1,3
      if ((faces(2*a-1)=='periodic').neqv.(faces(2*a)=='periodic')) return
    end do
    ok=.true.
  end function

  subroutine apply_parcel_boundaries(parcel,liquid,lower,upper,faces,change,ok)
    type(spray_parcel), intent(inout) :: parcel
    type(spray_liquid), intent(in) :: liquid
    real(dp), intent(in) :: lower(3),upper(3)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(out) :: change(5)
    logical, intent(out) :: ok
    type(spray_parcel) :: candidate
    real(dp) :: before(5)
    integer :: a,face,attempt
    ok=.false.; change=0
    if (.not.valid_spray_boundary_faces(faces).or..not.valid_spray_parcel(parcel)) return
    if (.not.all(ieee_is_finite([lower,upper])).or.any(upper<=lower)) return
    if (.not.valid_spray_liquid(liquid)) return
    if (maxval(abs([lower,upper]))>1e9_dp) return
    candidate=parcel; before=parcel_totals(parcel,liquid)
    do a=1,3
      if (faces(2*a)=='periodic') then
        candidate%position(a)=lower(a)+modulo(candidate%position(a)-lower(a),upper(a)-lower(a))
        cycle
      end if
      do attempt=1,32
        face=0
        if (candidate%position(a)<lower(a)) face=2*a-1
        if (candidate%position(a)>=upper(a)) face=2*a
        if (face==0) exit
        if (faces(face)=='outflow'.or.faces(face)=='inflow') then
          candidate%mass=0
          candidate%position=max(lower,min(nearest(upper,-1.0_dp),candidate%position))
          exit
        end if
        if (modulo(face,2)==0) then
          candidate%position(a)=2*upper(a)-candidate%position(a)
          candidate%position(a)=min(nearest(upper(a),-1.0_dp),candidate%position(a))
        else
          candidate%position(a)=2*lower(a)-candidate%position(a)
        end if
        candidate%velocity(a)=-candidate%velocity(a)
      end do
      if (attempt>32) return
      if (candidate%mass==0) exit
    end do
    change=parcel_totals(candidate,liquid)-before
    parcel=candidate; ok=.true.
  end subroutine

  subroutine particle_boundary_stencil(position,lower,spacing,cells,faces,indices,weights,ok)
    real(dp), intent(in) :: position(3),lower(3),spacing(3)
    integer, intent(in) :: cells(3)
    character(len=*), intent(in) :: faces(6)
    integer, intent(out) :: indices(3,8)
    real(dp), intent(out) :: weights(8)
    logical, intent(out) :: ok
    real(dp) :: coordinate(3),fraction(3)
    integer :: base(3),offset(3),a,b,c,n,d
    indices=0; weights=0; ok=.false.
    if (.not.valid_spray_boundary_faces(faces)) return
    if (all(faces=='periodic')) then
      call periodic_particle_stencil(position,lower,spacing,cells,indices,weights,ok)
      return
    end if
    if (any(cells<2).or..not.all(ieee_is_finite([position,lower,spacing]))) return
    if (any(spacing<1e-12_dp).or.any(spacing>1e9_dp)) return
    if (maxval(abs(position))>1e9_dp.or.maxval(abs(lower))>1e9_dp) return
    do d=1,3
      if (faces(2*d)=='periodic') then
        coordinate(d)=modulo(position(d)-lower(d),spacing(d)*cells(d))/spacing(d)-0.5_dp
        base(d)=floor(coordinate(d))
      else
        if (position(d)<lower(d).or.position(d)>lower(d)+spacing(d)*cells(d)) return
        coordinate(d)=max(0.0_dp,min(real(cells(d)-1,dp),(position(d)-lower(d))/spacing(d)-0.5_dp))
        base(d)=min(cells(d)-2,floor(coordinate(d)))
      end if
    end do
    fraction=coordinate-real(base,dp); n=0
    do c=0,1
      do b=0,1
        do a=0,1
          n=n+1; offset=[a,b,c]
          indices(:,n)=modulo(base+offset,cells)+1
          weights(n)=product(merge(fraction,1.0_dp-fraction,offset==1))
        end do
      end do
    end do
    ok=abs(sum(weights)-1.0_dp)<1e-14_dp.and.all(weights>=0)
  end subroutine

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

  subroutine build_spray_film(species, transport, primitive, gas_temperature, parcel, liquid, options, film, ok,wall_skin)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: primitive(:), gas_temperature
    type(spray_parcel), intent(in) :: parcel
    type(spray_liquid), intent(in) :: liquid
    type(spray_options), intent(in) :: options
    type(spray_film), intent(out) :: film
    logical, intent(out) :: ok
    logical, intent(in), optional :: wall_skin
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
    if (present(wall_skin)) then
      if (wall_skin) then
        ! WallFilm.H declares rule=1/2 despite its stale one-third-rule comment.
        tf=parcel%temperature+.5_dp*max(0.0_dp,gas_temperature-parcel%temperature)
        if (ys>y(v)) then
          yf(v)=.5_dp*(ys+y(v))
          do k=1,size(species)
            if (k/=v) yf(k)=y(k)*(1-yf(v))/(1-y(v))
          end do
        end if
      end if
    end if
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
      lower, spacing, parcels, interval, substeps, ok, message, boundary_faces,boundary_change,materials,breakup_options, &
      wall_options,wall_films)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquid
    type(spray_liquid), intent(in), optional :: materials(:)
    type(spray_liquid) :: material
    type(spray_options), intent(in) :: options
    type(spray_options), intent(in), optional :: breakup_options(:)
    type(spray_khrt_shedding_state) :: shedding
    type(spray_khrt_release) :: release
    real(dp) :: reynolds,dispersion_energy,release_axis(3)
    type(spray_options) :: local_breakup
    type(spray_wall_options), intent(in), optional :: wall_options
    type(spray_wall_film), allocatable, intent(inout), optional :: wall_films(:)
    type(spray_wall_options) :: impact_settings
    type(spray_wall_film), allocatable :: patches(:)
    type(spray_wall_film) :: deposit,next_patch
    type(spray_parcel), allocatable :: spawned(:)
    type(spray_parcel) :: identity
    real(dp), allocatable :: starts(:),film_starts(:)
    real(dp) :: film_conductivity
    real(dp) :: wall_exchange(5),wall_energy,film_mass_rate,film_gas_heat,film_wall_heat,film_temperature_rate
    integer :: patch_index,wall_face,axis
    logical :: hit
    real(dp), intent(inout) :: state(:,:,:,:), temperature(:,:,:)
    real(dp), intent(in) :: lower(3), spacing(3), interval
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    integer, intent(out) :: substeps
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    character(len=*), intent(in), optional :: boundary_faces(6)
    real(dp), intent(out), optional :: boundary_change(:)
    character(len=32) :: faces(6)
    real(dp) :: ledger(size(state,1)),exchange(5)
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
    impact_settings=spray_wall_options()
    if (present(wall_options)) impact_settings=wall_options
    faces='periodic'
    if (present(boundary_faces)) faces=boundary_faces
    ledger=0
    if (present(boundary_change)) boundary_change=0
    substeps=0; ok=.false.; message='Invalid spray coupling input'
    if (.not.allocated(parcels)) return
    if (.not.valid_spray_boundary_faces(faces)) return
    if (any(faces/='periodic').and..not.present(boundary_change)) return
    if (present(boundary_change)) then
      if (size(boundary_change)/=size(state,1)) return
    end if
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
    if (present(breakup_options)) then
      if (.not.present(materials)) return
      if (size(breakup_options)/=size(materials)) return
      do p=1,size(breakup_options)
        if (.not.valid_spray_options(breakup_options(p))) return
        if (breakup_options(p)%breakup_model/=options%breakup_model) return
      end do
    end if
    if (.not.valid_wall_options(impact_settings,faces)) return
    if (any(impact_settings%enabled).and..not.present(wall_films)) return
    allocate(patches(0))
    if (present(wall_films)) then
      if (.not.allocated(wall_films)) return
      patches=wall_films
    end if
    if (size(patches)>impact_settings%maximum_films) return
    allocate(film_starts(size(patches)));film_starts=0
    do p=1,size(patches)
      if (.not.valid_wall_film(patches(p))) return
      identity=spray_parcel();identity%material_id=patches(p)%material_id
      material=select_spray_liquid(identity,liquid,materials)
      if (.not.valid_spray_liquid(material)) return
      if (material%vapor_index>size(species)) return
      if (trim(species(material%vapor_index)%name)/=trim(material%vapor_name)) return
      if (liquid_heat_capacity_at_temperature(material,patches(p)%temperature)<=0) return
    end do
    do p=1,size(parcels)
      if (.not.options%khrt_children.and.(parcels(p)%khrt_initial_count/=0.or.parcels(p)%khrt_shed_volume/=0)) return
      if (.not. valid_spray_parcel(parcels(p))) return
      material=select_spray_liquid(parcels(p),liquid,materials)
      if (.not.valid_spray_liquid(material)) return
      if (material%vapor_index>size(species)) return
      if (trim(species(material%vapor_index)%name)/=trim(material%vapor_name)) return
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
    allocate(starts(size(drops)));starts=0;p=1
    do while(p<=size(drops))
      material=select_spray_liquid(drops(p),liquid,materials)
      local_breakup=options
      if (present(breakup_options)) then
        local_breakup=breakup_options(drops(p)%material_id)
        local_breakup%minimum_diameter=options%minimum_diameter
        local_breakup%maximum_substeps=options%maximum_substeps
      end if
      elapsed=starts(p)
      do while (elapsed < interval .and. drops(p)%mass > 0.0_dp)
        if (total_steps >= options%maximum_substeps) then
          message='Spray substep limit exceeded; complete update rolled back'; return
        end if
        call particle_boundary_stencil(drops(p)%position,lower,spacing,shape3,faces,indices,weights,local_ok)
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
        call build_spray_film(species,transport,primitive,tg,drops(p),material,options,film,local_ok)
        if (.not. local_ok) then
          message='Unsupported film state (including boiling) or invalid material/vapor data'; return
        end if
        sampled_pressure=primitive(5)
        call spray_transfer_rates(drops(p),material,film,options,invtau,heat,masscoef,local_ok)
        if (.not. local_ok) return
        d=parcel_diameter(drops(p),material); mdot=masscoef*d
        tdot=(heat*(tg-drops(p)%temperature)-mdot*(film%vapor_enthalpy-liquid_energy(material,drops(p)%temperature))) &
          /(drops(p)%mass*liquid_heat_capacity_at_temperature(material,drops(p)%temperature))
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
          call advance_parcel_frozen_film(drops(p),material,film,options,dt,next,transfer,local_ok)
          if (local_ok .and. next%mass>0.0_dp) then
            call nasa7_mass_properties(species(material%vapor_index),next%temperature,cp,cv,h,e,entropy,local_ok)
            if (local_ok .and. options%evaporation) then
              call saturation_mole_fraction(material,species(material%vapor_index)%molecular_weight, &
                next%temperature,sampled_pressure,xsat,local_ok)
            end if
          end if
          if (local_ok) then
            increase=0.0_dp; increase(1:5)=transfer; increase(5+material%vapor_index)=transfer(1)
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
          message='Spray feedback cannot preserve gas/material admissibility; complete update rolled back'; return
        end if
        release=spray_khrt_release();hit=.false.
        if (options%khrt_children) then
          ! A crossed wall owns the incident mass before any new KH release.
          call process_wall_event(accepted)
          if (.not.accepted) return
        end if
        if (options%khrt_children.and..not.hit) then
          shedding=spray_khrt_shedding_state(next%khrt_initial_count,next%khrt_shed_volume)
          reynolds=film%density*norm2(film%gas_velocity-next%velocity)*parcel_diameter(next,material)/film%viscosity
          call advance_khrt_shedding(next,material,local_breakup,sampled(1),film%gas_velocity,reynolds,dt, &
            options%khrt_reference_diameter_cubed,shedding,release,local_ok)
          if (local_ok) then
            next%khrt_initial_count=shedding%initial_count;next%khrt_shed_volume=shedding%shed_volume
          end if
        else if (.not.options%khrt_children) then
          call advance_spray_breakup(next,material,local_breakup,sampled(1),film%gas_velocity, &
            film%density,film%viscosity,dt,local_ok)
        end if
        if (.not.local_ok) then
          message='Spray breakup rejected; complete coupling step rolled back'; return
        end if
        if (release%mass>0) then
          if (size(drops)>options%maximum_parcels-4) return
          if (drops(size(drops))%id>huge(1_int64)-4_int64) return
          release_axis=next%velocity
          if (norm2(release_axis)<1e-12_dp) release_axis=film%gas_velocity-next%velocity
          call construct_khrt_children(next,release,release_axis,drops(size(drops))%id+1_int64, &
            spawned,dispersion_energy,local_ok)
          if (.not.local_ok) return
          do n=1,8
            i=indices(1,n);j=indices(2,n);k=indices(3,n)
            candidate(5,i,j,k)=candidate(5,i,j,k)+weights(n)*dispersion_energy/volume
            call reactive_conserved_to_primitive(species,candidate(:,i,j,k),temp(i,j,k),primitive,recovered,sound,local_ok)
            if (.not.local_ok) then
              message='KH child dispersion energy rejected; complete step rolled back';return
            end if
            temp(i,j,k)=recovered
          end do
          do n=1,size(spawned)
            call apply_parcel_boundaries(spawned(n),material,lower,lower+spacing*real(shape3,dp),faces,exchange,local_ok)
            if (.not.local_ok) return
            ledger(1:5)=ledger(1:5)+exchange
            ledger(5+material%vapor_index)=ledger(5+material%vapor_index)+exchange(1)
          end do
          starts=[starts,spread(min(interval,elapsed+dt),1,size(spawned))]
          drops=[drops,spawned]
        end if
        if (.not.options%khrt_children) then
          call process_wall_event(accepted)
          if (.not.accepted) return
        end if
        call apply_parcel_boundaries(next,material,lower,lower+spacing*real(shape3,dp),faces,exchange,local_ok)
        if (.not.local_ok) return
        ledger(1:5)=ledger(1:5)+exchange
        ledger(5+material%vapor_index)=ledger(5+material%vapor_index)+exchange(1)
        drops(p)=next
        elapsed=min(interval,elapsed+dt); total_steps=total_steps+1
      end do
      p=p+1
    end do
    if (impact_settings%evolve) then
      do patch_index=1,size(patches)
        identity=spray_parcel();identity%id=1;identity%material_id=patches(patch_index)%material_id
        identity%position=patches(patch_index)%position;identity%multiplicity=1
        material=select_spray_liquid(identity,liquid,materials)
        axis=maxloc(abs(patches(patch_index)%normal),dim=1)
        wall_face=2*axis-merge(1,0,patches(patch_index)%normal(axis)>0)
        if (.not.impact_settings%enabled(wall_face)) return
        elapsed=film_starts(patch_index)
        do while(elapsed<interval.and.patches(patch_index)%mass>0)
          if (total_steps>=options%maximum_substeps) return
          ! This proxy supplies only temperature/composition to gas-film property evaluation.
          identity%mass=min(patches(patch_index)%mass,1.0_dp);identity%temperature=patches(patch_index)%temperature
          call particle_boundary_stencil(identity%position,lower,spacing,shape3,faces,indices,weights,local_ok)
          if (.not.local_ok) return
          sampled=0;guess=0
          do n=1,8
            i=indices(1,n);j=indices(2,n);k=indices(3,n)
            old_cells(:,n)=candidate(:,i,j,k);old_t(n)=temp(i,j,k)
            sampled=sampled+weights(n)*old_cells(:,n);guess=guess+weights(n)*old_t(n)
          end do
          call reactive_conserved_to_primitive(species,sampled,guess,primitive,tg,sound,local_ok)
          if (.not.local_ok) return
          call build_spray_film(species,transport,primitive,tg,identity,material,options,film,local_ok,wall_skin=.true.)
          if (.not.local_ok) return
          film_conductivity=wall_liquid_conductivity_at_temperature(impact_settings, &
            patches(patch_index)%material_id,patches(patch_index)%temperature)
          call wall_film_transfer_rates(patches(patch_index),material,film,film_conductivity, &
            impact_settings%temperature(wall_face),options%heat,options%evaporation,film_mass_rate,film_gas_heat, &
            film_wall_heat,film_temperature_rate,local_ok)
          if (.not.local_ok) return
          dt=interval-elapsed
          if (film_mass_rate>0) dt=min(dt,options%max_fractional_change*patches(patch_index)%mass/film_mass_rate)
          if (abs(film_temperature_rate)>0) dt=min(dt, &
            min(2.0_dp,options%max_fractional_change*identity%temperature)/abs(film_temperature_rate))
          accepted=.false.
          do attempt=1,30
            if (dt<=0.or.elapsed+dt<=elapsed) exit
            call advance_wall_film_frozen(patches(patch_index),material,film,film_conductivity, &
              impact_settings%temperature(wall_face),options%heat,options%evaporation,dt,impact_settings%minimum_height, &
              next_patch,transfer,wall_energy,local_ok)
            if (local_ok.and.next_patch%mass>0) local_ok= &
              wall_liquid_conductivity_at_temperature(impact_settings,next_patch%material_id,next_patch%temperature)>0
            if (local_ok) then
              increase=0;increase(1:5)=transfer;increase(5+material%vapor_index)=transfer(1)
              do n=1,8
                i=indices(1,n);j=indices(2,n);k=indices(3,n)
                candidate(:,i,j,k)=old_cells(:,n)+weights(n)*increase/volume
                call reactive_conserved_to_primitive(species,candidate(:,i,j,k),old_t(n),primitive,recovered,sound,local_ok)
                if (.not.local_ok) exit
                temp(i,j,k)=recovered
              end do
            end if
            if (local_ok) then
              accepted=.true.;exit
            end if
            do n=1,8
              i=indices(1,n);j=indices(2,n);k=indices(3,n)
              candidate(:,i,j,k)=old_cells(:,n);temp(i,j,k)=old_t(n)
            end do
            dt=.5_dp*dt
          end do
          if (.not.accepted) then
            message='Film evolution rejected; complete coupling step rolled back';return
          end if
          patches(patch_index)=next_patch;ledger(5)=ledger(5)-wall_energy
          elapsed=min(interval,elapsed+dt);total_steps=total_steps+1
        end do
      end do
      patches=pack(patches,patches%mass>0)
    end if
    call split_spray_parcels(drops,options%maximum_multiplicity,options%maximum_parcels,local_ok)
    if (.not.local_ok) then
      message='Parcel splitting exceeds storage or ID bounds; complete step rolled back';return
    end if
    state=candidate; temperature=temp; parcels=drops
    if (present(wall_films)) call move_alloc(patches,wall_films)
    if (present(boundary_change)) boundary_change=ledger
    substeps=total_steps; ok=.true.; message=''
  contains
    subroutine process_wall_event(event_ok)
      logical, intent(out) :: event_ok
      event_ok=.false.;hit=.false.
      if (any(impact_settings%enabled).and.next%mass>0) then
        if (drops(size(drops))%id>huge(1_int64)-4_int64) return
        call impact_at_planar_wall(next,drops(p)%position,material,local_breakup, &
          species(material%vapor_index)%molecular_weight,sampled_pressure,lower,lower+spacing*real(shape3,dp), &
          impact_settings,patches,drops(size(drops))%id+1_int64,spawned,deposit,wall_exchange,hit,local_ok)
        if (.not.local_ok) then
          message='Wall impact rejected; complete coupling step rolled back';return
        end if
        if (hit) then
          if (size(drops)+size(spawned)>options%maximum_parcels) return
          do n=1,size(spawned)
            if (any(spawned(n)%position<lower).or.any(spawned(n)%position>=lower+spacing*real(shape3,dp))) return
          end do
          if (size(spawned)>0) then
            starts=[starts,spread(min(interval,elapsed+dt),1,size(spawned))]
            drops=[drops,spawned]
          end if
          if (deposit%mass>0) then
            if (size(patches)>=impact_settings%maximum_films) return
            patches=[patches,deposit]
            film_starts=[film_starts,min(interval,elapsed+dt)]
          end if
          ledger(1:5)=ledger(1:5)-wall_exchange
        end if
      end if
      event_ok=.true.
    end subroutine
  end subroutine

  subroutine advance_mixture_coupling_3d(species,transport,components,options,state,temperature,lower,spacing, &
      parcels,component_mass,interval,substeps,ok,message,boundary_faces,boundary_change, &
      wall_options,wall_films,film_component_mass)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(spray_wall_options), optional, intent(in) :: wall_options
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_component_mass(:,:)
    type(spray_wall_options) :: impact_settings
    type(spray_options) :: film_options
    type(spray_wall_film), allocatable :: patches(:)
    type(spray_wall_film) :: deposit
    real(dp), allocatable :: patch_mass(:,:),child_mass(:,:),film_starts(:)
    integer :: film_steps
    logical :: hit
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    real(dp) :: deposited_mass(size(components)),wall_exchange(5),film_ledger(size(state,1)),mw(size(components))
    real(dp), intent(in) :: lower(3),spacing(3),interval
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    integer, intent(out) :: substeps
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(out) :: boundary_change(:)
    real(dp), allocatable :: gas(:,:,:,:),temp(:,:,:),masses(:,:),split_masses(:,:)
    type(spray_parcel), allocatable :: drops(:),children(:)
    type(spray_khrt_shedding_state) :: shedding
    type(spray_khrt_release) :: release
    real(dp), allocatable :: starts(:)
    real(dp) :: reynolds,dispersion_energy,release_axis(3)
    type(spray_parcel) :: next
    type(spray_liquid) :: effective
    type(spray_film) :: mixture_film
    real(dp) :: y(size(components)),rates(size(components)),hv(size(components)),next_mass(size(components))
    real(dp) :: before_boundary(size(components)),rho,cp,e,heat,invtau,tdot,elapsed,dt,speed,old_mass,volume
    real(dp) :: sampled(size(state,1)),primitive(size(state,1)),old_cells(size(state,1),8),old_t(8)
    real(dp) :: change(size(state,1)),ledger(size(state,1)),weights(8),exchange(5),tg,recovered,sound,guess
    integer :: shape3(3),indices(3,8),p,n,i,j,k,attempt,count,old_count,parent,v
    character(len=32) :: faces(6)
    logical :: valid,accepted
    ok=.false.;substeps=0;message='Invalid mixture coupling state'
    if (present(boundary_change)) boundary_change=0
    if (.not.allocated(parcels).or..not.allocated(component_mass)) return
    if (.not.valid_spray_options(options)) return
    if (size(components)<1.or.size(components)>min(16,size(species))) return
    if (any(shape(component_mass)/=[size(components),size(parcels)])) return
    if (size(state,1)/=5+size(species).or.size(parcels)>options%maximum_parcels) return
    shape3=shape(temperature)
    if (any(shape3<2).or.any(shape(state(1,:,:,:))/=shape3)) return
    if (.not.all(ieee_is_finite([lower,spacing,interval]))) return
    if (interval<0.or.interval>1e6_dp.or.any(spacing<1e-12_dp).or.any(spacing>1e9_dp)) return
    faces='periodic'
    if (present(boundary_faces)) faces=boundary_faces
    if (.not.valid_spray_boundary_faces(faces)) return
    impact_settings=spray_wall_options()
    if (present(wall_options)) impact_settings=wall_options
    if (.not.valid_wall_options(impact_settings,faces)) return
    if (present(wall_films).neqv.present(film_component_mass)) return
    if (any(impact_settings%enabled).and..not.present(wall_films)) return
    allocate(patches(0),patch_mass(size(components),0))
    if (present(wall_films)) then
      if (.not.allocated(wall_films).or..not.allocated(film_component_mass)) return
      if (any(shape(film_component_mass)/=[size(components),size(wall_films)])) return
      patches=wall_films;patch_mass=film_component_mass
    end if
    if (size(patches)>impact_settings%maximum_films) return
    do p=1,size(patches)
      if (.not.valid_wall_component_inventory(patches(p),patch_mass(:,p))) return
    end do
    allocate(film_starts(size(patches)));film_starts=0
    do p=1,size(components)
      if (.not.valid_spray_liquid(components(p))) return
      if (components(p)%vapor_index>size(species)) return
      mw(p)=species(components(p)%vapor_index)%molecular_weight
    end do
    if (present(boundary_change)) then
      if (size(boundary_change)/=size(state,1)) return
    end if
    do p=1,size(parcels)
      if (.not.valid_component_inventory(parcels(p),component_mass(:,p))) return
      if (.not.options%khrt_children.and.(parcels(p)%khrt_initial_count/=0.or.parcels(p)%khrt_shed_volume/=0)) return
      if (p>1) then
        if (parcels(p)%id<=parcels(p-1)%id) return
      end if
    end do
    gas=state;temp=temperature
    call recover_reactive_temperatures_3d(species,state,temperature,shape3(1),shape3(2),shape3(3),temp,valid)
    if (.not.valid) return
    if (any(abs(temp-temperature)>1e-9_dp*max(1.0_dp,abs(temperature)))) return
    temp=temperature;drops=parcels;masses=component_mass;ledger=0;count=0;volume=product(spacing)
    allocate(starts(size(drops)));starts=0;p=1
    do while(p<=size(drops))
      elapsed=starts(p)
      do while(elapsed<interval.and.drops(p)%mass>0)
        if (count>=options%maximum_substeps) then
          message='Mixture substep limit; complete update rolled back';return
        end if
        call particle_boundary_stencil(drops(p)%position,lower,spacing,shape3,faces,indices,weights,valid)
        if (.not.valid) return
        sampled=0;guess=0
        do n=1,8
          i=indices(1,n);j=indices(2,n);k=indices(3,n)
          old_cells(:,n)=gas(:,i,j,k);old_t(n)=temp(i,j,k)
          sampled=sampled+weights(n)*old_cells(:,n);guess=guess+weights(n)*old_t(n)
        end do
        call reactive_conserved_to_primitive(species,sampled,guess,primitive,tg,sound,valid)
        if (.not.valid) return
        y=masses(:,p)/sum(masses(:,p))
        call mixture_transfer_rates(species,transport,primitive,tg,drops(p),components,y,options, &
          invtau,heat,rates,hv,valid)
        if (.not.valid) return
        call liquid_mixture_properties(components,y,drops(p)%temperature,rho,cp,e,valid)
        if (.not.valid) return
        dt=interval-elapsed
        if (invtau>0) dt=min(dt,options%max_fractional_change/invtau)
        tdot=heat*(tg-drops(p)%temperature)
        do n=1,size(components)
          if (rates(n)>0.and.y(n)>mixture_trace_fraction) &
            dt=min(dt,options%max_fractional_change*masses(n,p)/rates(n))
          tdot=tdot-rates(n)*(hv(n)-liquid_energy(components(n),drops(p)%temperature))
        end do
        tdot=tdot/(drops(p)%mass*cp)
        if (abs(tdot)>0) dt=min(dt,min(2.0_dp,options%max_fractional_change*drops(p)%temperature)/abs(tdot))
        speed=max(maxval(abs(drops(p)%velocity)),maxval(abs(primitive(2:4))))
        if (speed>0) dt=min(dt,options%particle_cfl*minval(spacing)/speed)
        accepted=.false.
        do attempt=1,30
          if (dt<=0.or.elapsed+dt<=elapsed) exit
          next=drops(p);next_mass=masses(:,p)
          call advance_mixture_frozen_gas(species,transport,primitive,tg,next,components,next_mass,options,dt,change,valid)
          if (valid) then
            do n=1,8
              i=indices(1,n);j=indices(2,n);k=indices(3,n)
              gas(:,i,j,k)=old_cells(:,n)+weights(n)*change/volume
              call reactive_conserved_to_primitive(species,gas(:,i,j,k),old_t(n),primitive,recovered,sound,valid)
              if (.not.valid) exit
              temp(i,j,k)=recovered
            end do
          end if
          if (valid) then
            accepted=.true.;exit
          end if
          do n=1,8
            i=indices(1,n);j=indices(2,n);k=indices(3,n)
            gas(:,i,j,k)=old_cells(:,n);temp(i,j,k)=old_t(n)
          end do
          call reactive_conserved_to_primitive(species,sampled,guess,primitive,tg,sound,valid)
          if (.not.valid) return
          dt=.5_dp*dt
        end do
        if (.not.accepted) then
          message='Mixture feedback rejected; complete update rolled back';return
        end if
        if (next%mass>0) then
          y=next_mass/sum(next_mass)
          call liquid_mixture_properties(components,y,next%temperature,rho,cp,e,valid)
          if (.not.valid) return
        end if
        effective=components(1);effective%density=rho;effective%density_model='constant'
        effective%density_coefficients=0;effective%density_temperature_range=0
        effective%heat_capacity=cp;effective%heat_capacity_model='constant'
        effective%heat_capacity_coefficients=0;effective%heat_capacity_temperature_range=0
        effective%reference_temperature=next%temperature;effective%reference_energy=e
        effective%boiling_model='reject';effective%critical_temperature=0
        old_mass=next%mass
        ! Restore the sampled primitive: the feedback recovery above used it as scratch.
        call reactive_conserved_to_primitive(species,sampled,guess,primitive,tg,sound,valid)
        if (.not.valid) return
        hit=.false.
        if (any(impact_settings%enabled).and.next%mass>0) then
          if (drops(size(drops))%id>huge(1_int64)-4_int64) return
          call mixture_impact_at_planar_wall(next,drops(p)%position,components,next_mass,options,mw,primitive(5), &
            lower,lower+spacing*real(shape3,dp),impact_settings,patches,drops(size(drops))%id+1_int64,children,child_mass, &
            deposit,deposited_mass,wall_exchange,hit,valid)
          if (.not.valid) then
            message='Mixed wall impact rejected; complete update rolled back';return
          end if
          if (hit) then
            if (size(children)>options%maximum_parcels-size(drops)) return
            do n=1,size(children)
              if (any(children(n)%position<lower).or.any(children(n)%position>=lower+spacing*real(shape3,dp))) return
            end do
            if (size(children)>0) then
              old_count=size(drops)
              allocate(split_masses(size(components),old_count+size(children)))
              split_masses(:,:old_count)=masses;split_masses(:,old_count+1:)=child_mass
              call move_alloc(split_masses,masses)
              starts=[starts,spread(min(interval,elapsed+dt),1,size(children))];drops=[drops,children]
            end if
            if (deposit%mass>0) then
              if (size(patches)>=impact_settings%maximum_films) return
              patch_mass=reshape([reshape(patch_mass,[size(patch_mass)]),deposited_mass], &
                [size(components),size(patches)+1])
              patches=[patches,deposit];film_starts=[film_starts,min(interval,elapsed+dt)]
            end if
            ledger(1:5)=ledger(1:5)-wall_exchange
          end if
        end if
        release=spray_khrt_release()
        if (next%mass>0.and.options%breakup_model/='none'.and..not.hit) then
          call build_mixture_film(species,transport,primitive,tg,next,components,y,options,mixture_film,rates,hv,valid)
          if (.not.valid) return
          if (options%khrt_children) then
            shedding=spray_khrt_shedding_state(next%khrt_initial_count,next%khrt_shed_volume)
            reynolds=mixture_film%density*norm2(primitive(2:4)-next%velocity)* &
              parcel_diameter(next,effective)/mixture_film%viscosity
            call advance_khrt_shedding(next,effective,options,sampled(1),primitive(2:4),reynolds,dt, &
              options%khrt_reference_diameter_cubed,shedding,release,valid)
            if (valid) then
              next%khrt_initial_count=shedding%initial_count;next%khrt_shed_volume=shedding%shed_volume
            end if
          else
            call advance_spray_breakup(next,effective,options,sampled(1),primitive(2:4),mixture_film%density, &
              mixture_film%viscosity,dt,valid)
          end if
          if (.not.valid) return
          next_mass=next_mass*(next%mass/old_mass)
        end if
        if (release%mass>0) then
          if (size(drops)>options%maximum_parcels-4) then
            message='Mixture KH children exceed parcel limit; complete update rolled back';return
          end if
          if (drops(size(drops))%id>huge(1_int64)-4_int64) return
          release_axis=next%velocity
          if (norm2(release_axis)<1e-12_dp) release_axis=primitive(2:4)-next%velocity
          call construct_khrt_children(next,release,release_axis,drops(size(drops))%id+1_int64,children,dispersion_energy,valid)
          if (.not.valid) return
          do n=1,8
            i=indices(1,n);j=indices(2,n);k=indices(3,n)
            gas(5,i,j,k)=gas(5,i,j,k)+weights(n)*dispersion_energy/volume
            call reactive_conserved_to_primitive(species,gas(:,i,j,k),temp(i,j,k),primitive,recovered,sound,valid)
            if (.not.valid) then
              message='Mixture KH dispersion energy rejected; complete update rolled back';return
            end if
            temp(i,j,k)=recovered
          end do
          old_count=size(drops)
          allocate(split_masses(size(components),old_count+size(children)));split_masses(:,:old_count)=masses
          do n=1,size(children)
            call apply_parcel_boundaries(children(n),effective,lower,lower+spacing*real(shape3,dp),faces,exchange,valid)
            if (.not.valid) return
            ledger(1:5)=ledger(1:5)+exchange
            do v=1,size(components)
              ledger(5+components(v)%vapor_index)=ledger(5+components(v)%vapor_index)+exchange(1)*y(v)
            end do
            split_masses(:,old_count+n)=children(n)%mass*y
          end do
          call move_alloc(split_masses,masses)
          starts=[starts,spread(min(interval,elapsed+dt),1,size(children))]
          drops=[drops,children]
        end if
        before_boundary=next_mass
        call apply_parcel_boundaries(next,effective,lower,lower+spacing*real(shape3,dp),faces,exchange,valid)
        if (.not.valid) return
        if (next%mass==0) next_mass=0
        ledger(1:5)=ledger(1:5)+exchange
        do n=1,size(components)
          v=components(n)%vapor_index
          ledger(5+v)=ledger(5+v)+next%multiplicity*(next_mass(n)-before_boundary(n))
        end do
        drops(p)=next;masses(:,p)=next_mass;elapsed=min(interval,elapsed+dt);count=count+1
      end do
      p=p+1
    end do
    if (size(patches)>0.and.impact_settings%evolve) then
      film_options=options;film_options%maximum_substeps=max(1,options%maximum_substeps-count)
      call advance_mixture_film_coupling_3d(species,transport,components,film_options,impact_settings,gas,temp, &
        patches,patch_mass,lower,spacing,interval,faces,film_steps,film_ledger,valid,message,film_starts)
      if (.not.valid) return
      if (film_steps>options%maximum_substeps-count) then
        message='Mixed particle/film substep limit; complete update rolled back';return
      end if
      count=count+film_steps;ledger=ledger+film_ledger
    end if
    old_count=size(drops)
    call split_spray_parcels(drops,options%maximum_multiplicity,options%maximum_parcels,valid)
    if (.not.valid) return
    if (size(drops)>old_count) then
      allocate(split_masses(size(components),size(drops)));split_masses(:,:old_count)=masses;parent=1
      do p=old_count+1,size(drops)
        do while(parent<old_count)
          if (drops(parent)%id>=drops(p)%parent_id) exit
          parent=parent+1
        end do
        if (drops(parent)%id/=drops(p)%parent_id) return
        split_masses(:,p)=masses(:,parent)
      end do
      call move_alloc(split_masses,masses)
    end if
    state=gas;temperature=temp;parcels=drops;call move_alloc(masses,component_mass)
    if (present(wall_films)) then
      call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_component_mass)
    end if
    if (present(boundary_change)) boundary_change=ledger
    substeps=count;ok=.true.;message=''
  end subroutine

  subroutine advance_mixture_film_coupling_3d(species,transport,components,options,wall_options,state,temperature, &
      films,film_mass,lower,spacing,interval,faces,substeps,boundary_change,ok,message,start_times)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(spray_wall_options), intent(in) :: wall_options
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    type(spray_wall_film), allocatable, intent(inout) :: films(:)
    real(dp), allocatable, intent(inout) :: film_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3),interval
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in), optional :: start_times(:)
    integer, intent(out) :: substeps
    real(dp), intent(out) :: boundary_change(:)
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    type(spray_wall_film), allocatable :: patches(:)
    type(spray_wall_film) :: next
    real(dp), allocatable :: gas(:,:,:,:),temp(:,:,:),masses(:,:),kept(:,:)
    real(dp) :: sampled(size(state,1)),primitive(size(state,1)),change(size(state,1)),ledger(size(state,1))
    real(dp) :: old_cells(size(state,1),8),old_t(8),weights(8),guess,tg,sound,recovered,volume,elapsed,dt,conductivity
    real(dp) :: rates(size(components)),hv(size(components)),next_mass(size(components)),qg,qw,tdot,wall_energy
    integer :: shape3(3),indices(3,8),p,n,i,j,k,c,axis,face,attempt,total_steps
    logical :: valid,accepted
    ok=.false.;substeps=0;boundary_change=0;message='Invalid mixed film coupling state'
    if (.not.allocated(films).or..not.allocated(film_mass)) return
    if (.not.valid_spray_options(options).or..not.valid_wall_options(wall_options,faces)) return
    if (size(components)<1.or.size(components)>min(16,size(species))) return
    if (any(shape(film_mass)/=[size(components),size(films)]).or.size(films)>wall_options%maximum_films) return
    if (size(state,1)/=5+size(species).or.size(boundary_change)/=size(state,1)) return
    if (.not.all(ieee_is_finite([lower,spacing,interval]))) return
    if (any(spacing<=0).or.any(spacing>1e6_dp).or.interval<0.or.interval>1e6_dp) return
    shape3=shape(temperature)
    if (any(shape3<2).or.any(shape(state(1,:,:,:))/=shape3)) return
    if (present(start_times)) then
      if (size(start_times)/=size(films)) return
      if (.not.all(ieee_is_finite(start_times))) return
      if (any(start_times<0).or.any(start_times>interval)) return
    end if
    do p=1,size(films)
      if (.not.valid_wall_component_inventory(films(p),film_mass(:,p))) return
    end do
    do c=1,size(components)
      if (.not.valid_spray_liquid(components(c))) return
      if (components(c)%vapor_index>size(species)) return
      if (trim(species(components(c)%vapor_index)%name)/=trim(components(c)%vapor_name)) return
    end do
    message='Invalid mixed film gas thermodynamics';temp=temperature
    call recover_reactive_temperatures_3d(species,state,temperature,shape3(1),shape3(2),shape3(3),temp,valid)
    if (.not.valid) return
    if (any(abs(temp-temperature)>1e-9_dp*max(1.0_dp,abs(temperature)))) return
    gas=state;temp=temperature;patches=films;masses=film_mass;ledger=0;total_steps=0;volume=product(spacing)
    if (wall_options%evolve) then
      do p=1,size(patches)
        axis=maxloc(abs(patches(p)%normal),dim=1)
        face=2*axis-merge(1,0,patches(p)%normal(axis)>0)
        if (.not.wall_options%enabled(face)) return
        elapsed=0
        if (present(start_times)) elapsed=start_times(p)
        do while(elapsed<interval.and.patches(p)%mass>0)
          if (total_steps>=options%maximum_substeps) then
            message='Mixed film substep limit; complete update rolled back';return
          end if
          message='Invalid mixed film interpolation stencil'
          call particle_boundary_stencil(patches(p)%position,lower,spacing,shape3,faces,indices,weights,valid)
          if (.not.valid) return
          sampled=0;guess=0
          do n=1,8
            i=indices(1,n);j=indices(2,n);k=indices(3,n)
            old_cells(:,n)=gas(:,i,j,k);old_t(n)=temp(i,j,k)
            sampled=sampled+weights(n)*old_cells(:,n);guess=guess+weights(n)*old_t(n)
          end do
          call reactive_conserved_to_primitive(species,sampled,guess,primitive,tg,sound,valid)
          if (.not.valid) return
          conductivity=wall_liquid_conductivity_at_temperature(wall_options,patches(p)%material_id,patches(p)%temperature)
          message='Invalid mixed film transfer rates'
          call mixture_wall_film_rates(patches(p),components,masses(:,p),species,transport,primitive,tg,options, &
            conductivity,wall_options%temperature(face),rates,hv,qg,qw,tdot,valid)
          if (.not.valid) return
          dt=interval-elapsed
          do c=1,size(components)
            if (rates(c)>0.and.masses(c,p)>mixture_trace_fraction*patches(p)%mass) &
              dt=min(dt,options%max_fractional_change*masses(c,p)/rates(c))
          end do
          if (abs(tdot)>0) dt=min(dt,min(2.0_dp,options%max_fractional_change*patches(p)%temperature)/abs(tdot))
          accepted=.false.
          do attempt=1,30
            if (dt<=0.or.elapsed+dt<=elapsed) exit
            next=patches(p);next_mass=masses(:,p)
            call advance_mixture_wall_film_exchange(next,components,next_mass,rates*dt,hv,qg*dt,qw*dt, &
              wall_options%minimum_height,change,wall_energy,valid,deplete_traces=.true.)
            if (valid.and.next%mass>0) valid= &
              wall_liquid_conductivity_at_temperature(wall_options,next%material_id,next%temperature)>0
            if (valid) then
              do n=1,8
                i=indices(1,n);j=indices(2,n);k=indices(3,n)
                gas(:,i,j,k)=gas(:,i,j,k)+weights(n)*change/volume
              end do
              do n=1,8
                i=indices(1,n);j=indices(2,n);k=indices(3,n)
                call reactive_conserved_to_primitive(species,gas(:,i,j,k),old_t(n),primitive,recovered,sound,valid)
                if (.not.valid) exit
                temp(i,j,k)=recovered
              end do
            end if
            if (valid) then
              accepted=.true.;exit
            end if
            do n=1,8
              i=indices(1,n);j=indices(2,n);k=indices(3,n)
              gas(:,i,j,k)=old_cells(:,n);temp(i,j,k)=old_t(n)
            end do
            dt=.5_dp*dt
          end do
          if (.not.accepted) then
            message='Mixed film feedback rejected; complete update rolled back';return
          end if
          patches(p)=next;masses(:,p)=next_mass;ledger(5)=ledger(5)-wall_energy
          elapsed=min(interval,elapsed+dt);total_steps=total_steps+1
        end do
      end do
    end if
    allocate(kept(size(components),count(patches%mass>0)))
    do c=1,size(components)
      kept(c,:)=pack(masses(c,:),patches%mass>0)
    end do
    patches=pack(patches,patches%mass>0)
    state=gas;temperature=temp;films=patches;call move_alloc(kept,film_mass)
    boundary_change=ledger;substeps=total_steps;ok=.true.;message=''
  end subroutine

  subroutine mixture_system_totals(state,spacing,parcels,components,component_mass,totals,ok,wall_films,film_mass)
    real(dp), intent(in) :: state(:,:,:,:),spacing(3),component_mass(:,:)
    type(spray_parcel), intent(in) :: parcels(:)
    type(spray_wall_film), optional, intent(in) :: wall_films(:)
    real(dp), optional, intent(in) :: film_mass(:,:)
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(inout) :: totals(:)
    logical, intent(out) :: ok
    real(dp) :: candidate(size(state,1)),y(size(components)),rho,cp,e,m
    integer :: p,n,v
    logical :: valid
    ok=.false.
    if (size(totals)/=size(state,1).or.size(components)<1.or.size(components)>16) return
    if (any(shape(component_mass)/=[size(components),size(parcels)])) return
    if (.not.all(ieee_is_finite(spacing)).or.any(spacing<=0)) return
    do n=1,size(components)
      if (.not.valid_spray_liquid(components(n))) return
      if (components(n)%vapor_index>size(state,1)-5) return
    end do
    do n=1,size(candidate)
      candidate(n)=sum(state(n,:,:,:))*product(spacing)
    end do
    do p=1,size(parcels)
      if (.not.valid_component_inventory(parcels(p),component_mass(:,p))) return
      if (parcels(p)%mass==0) cycle
      y=component_mass(:,p)/sum(component_mass(:,p))
      call liquid_mixture_properties(components,y,parcels(p)%temperature,rho,cp,e,valid)
      if (.not.valid) return
      m=parcels(p)%mass*parcels(p)%multiplicity
      candidate(1)=candidate(1)+m;candidate(2:4)=candidate(2:4)+m*parcels(p)%velocity
      candidate(5)=candidate(5)+m*(e+.5_dp*sum(parcels(p)%velocity**2))
      do n=1,size(components)
        v=components(n)%vapor_index
        candidate(5+v)=candidate(5+v)+parcels(p)%multiplicity*component_mass(n,p)
      end do
    end do
    if (present(wall_films).neqv.present(film_mass)) return
    if (present(wall_films)) then
      if (any(shape(film_mass)/=[size(components),size(wall_films)])) return
      do p=1,size(wall_films)
        if (.not.valid_wall_component_inventory(wall_films(p),film_mass(:,p))) return
        call liquid_mixture_properties(components,film_mass(:,p)/wall_films(p)%mass,wall_films(p)%temperature,rho,cp,e,valid)
        if (.not.valid) return
        candidate(1)=candidate(1)+wall_films(p)%mass
        candidate(5)=candidate(5)+wall_films(p)%mass*e
        do n=1,size(components)
          candidate(5+components(n)%vapor_index)=candidate(5+components(n)%vapor_index)+film_mass(n,p)
        end do
      end do
    end if
    if (.not.all(ieee_is_finite(candidate))) return
    totals=candidate;ok=.true.
  end subroutine

  subroutine spray_system_totals(state, spacing, parcels, liquid, totals,materials,wall_films)
    real(dp), intent(in) :: state(:,:,:,:), spacing(3)
    type(spray_parcel), intent(in) :: parcels(:)
    type(spray_liquid), intent(in) :: liquid
    type(spray_liquid), intent(in), optional :: materials(:)
    type(spray_wall_film), intent(in), optional :: wall_films(:)
    type(spray_parcel) :: identity
    type(spray_liquid) :: material
    real(dp), intent(out) :: totals(size(state,1))
    real(dp) :: drop(5)
    integer :: v,p
    do v=1,size(state,1)
      totals(v)=sum(state(v,:,:,:))*product(spacing)
    end do
    do p=1,size(parcels)
      material=select_spray_liquid(parcels(p),liquid,materials)
      drop=parcel_totals(parcels(p),material)
      totals(1:5)=totals(1:5)+drop
      totals(5+material%vapor_index)=totals(5+material%vapor_index)+drop(1)
    end do
    if (present(wall_films)) then
      do p=1,size(wall_films)
        identity=spray_parcel();identity%material_id=wall_films(p)%material_id
        material=select_spray_liquid(identity,liquid,materials)
        drop=wall_film_totals(wall_films(p),material)
        totals(1:5)=totals(1:5)+drop
        totals(5+material%vapor_index)=totals(5+material%vapor_index)+drop(1)
      end do
    end if
  end subroutine
end module spray_coupling_3d_mod
