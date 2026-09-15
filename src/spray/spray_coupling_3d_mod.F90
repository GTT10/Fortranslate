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
  use spray_breakup_mod
  implicit none
  private
  public :: periodic_particle_stencil, build_spray_film, advance_spray_coupling_3d
  public :: spray_system_totals, valid_spray_boundary_faces, apply_parcel_boundaries
  public :: particle_boundary_stencil
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
