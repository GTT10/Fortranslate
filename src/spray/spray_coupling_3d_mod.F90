module spray_coupling_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use spray_wall_impact_mod
  use nasa7_thermo_mod, only: nasa7_species, nasa7_mass_properties
  use gas_transport_mod, only: gas_transport_species, compatible_transport_database
  use mixture_thermo_mod, only: mixture_mass_properties
  use mixture_transport_mod, only: mixture_transport_coefficients
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
    real(dp) :: fuel_diffusion(size(species))
    real(dp) :: xs,ys,ambient_weight,cv,h,e,s,mw,r,cp,gamma,species_cp
    logical :: vapor_active
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
    if (.not.all(ieee_is_finite(species%molecular_weight))) return
    if (any(species%molecular_weight<1e-3_dp).or.any(species%molecular_weight>1e5_dp)) return
    ys=0;vapor_active=.false.
    if (options%evaporation .and. y(v) < 1.0_dp-1.0e-12_dp) then
      h=0
      if (liquid%boiling_model=='pelemp-cap') then
        call nasa7_mass_properties(species(v),parcel%temperature,species_cp,cv,h,e,s,ok)
        if (.not.ok) return
      end if
      call saturation_mole_fraction(liquid,species(v)%molecular_weight,parcel%temperature,primitive(5),xs,ok, &
        vapor_enthalpy=h)
      if (.not. ok) return
      ! Pinned Drag.H uses the complete ambient mixture, including vapor.
      ambient_weight=1/sum(y/species%molecular_weight)
      vapor_active=xs>y(v)*ambient_weight/species(v)%molecular_weight
      if (vapor_active) then
        ys=xs*species(v)%molecular_weight/ &
          (xs*resolved_liquid_molecular_weight(liquid,species(v)%molecular_weight)+(1-xs)*ambient_weight)
        if (ys>=1) then
          ok=.false.;return
        end if
      end if
    end if
    yf=y
    if (vapor_active) then
      yf=y*(1.0_dp-(2.0_dp*ys+y(v))/3.0_dp)/(1.0_dp-y(v))
      yf(v)=(2.0_dp*ys+y(v))/3.0_dp
    end if
    tf=(2.0_dp*parcel%temperature+gas_temperature)/3.0_dp
    if(options%abramzon_sirignano.and..not.options%evaporation) tf=gas_temperature
    if (options%abramzon_sirignano.and.options%evaporation) tf=parcel%temperature+ &
      max(0.0_dp,gas_temperature-parcel%temperature)/3
    if (present(wall_skin)) then
      if (wall_skin) then
        ! WallFilm.H declares rule=1/2 despite its stale one-third-rule comment.
        tf=parcel%temperature+.5_dp*max(0.0_dp,gas_temperature-parcel%temperature)
        if (vapor_active) then
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
    call film_fuel_diffusion_coefficients(species,transport,yf,tf,primitive(5),diffusion,fuel_diffusion,ok)
    if (.not. ok) return
    film%diffusivity=fuel_diffusion(v)
    if (vapor_active) then
      call nasa7_mass_properties(species(v),tf,species_cp,cv,h,e,s,ok)
      if (.not.ok) return
      film%vapor_cp=species_cp
    end if
    ! At pure vapor the mass flux is zero; diffusivity has no mass-transfer role.
    call nasa7_mass_properties(species(v),parcel%temperature,species_cp,cv,h,e,s,ok)
    if (.not. ok) return
    film%vapor_enthalpy=h
    film%gas_velocity=primitive(2:4)
    film%gas_temperature=gas_temperature
    film%pressure=primitive(5)
    film%liquid_molecular_weight=resolved_liquid_molecular_weight(liquid,species(v)%molecular_weight)
    film%ambient_fraction=merge(y(v),0.0_dp,vapor_active)
    film%surface_fraction=ys
    ! Inactive/pure-vapor atmospheres retain ambient transport and zero mass flux.
    ok=.true.
  end subroutine

  subroutine advance_spray_coupling_3d(species, transport, liquid, options, state, temperature, &
      lower, spacing, parcels, interval, substeps, ok, message, boundary_faces,boundary_change,materials,breakup_options, &
      wall_options,wall_films,force_change)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquid
    type(spray_liquid), intent(in), optional :: materials(:)
    type(spray_liquid) :: material,breakup_material
    type(spray_options), intent(in) :: options
    real(dp), optional, intent(inout) :: force_change(5)
    real(dp) :: force_ledger(5),parcel_force(5)
    type(spray_options), intent(in), optional :: breakup_options(:)
    type(spray_khrt_shedding_state) :: shedding
    type(spray_khrt_release) :: release
    real(dp) :: reynolds,child_exchange(4),release_axis(3),etab_delta
    type(spray_options) :: local_breakup
    type(spray_wall_options), intent(in), optional :: wall_options
    type(spray_wall_film), allocatable, intent(inout), optional :: wall_films(:)
    type(spray_wall_options) :: impact_settings
    type(spray_wall_film), allocatable :: patches(:)
    type(spray_wall_film) :: deposit,next_patch
    type(spray_parcel), allocatable :: spawned(:)
    type(spray_parcel) :: identity
    real(dp), allocatable :: starts(:),film_starts(:),etab_tangent(:),etab_last_dt(:),etab_origins(:,:)
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
    type(spray_parcel) :: next,birth_parent
    type(spray_film) :: film
    real(dp) :: weights(8), old_cells(size(state,1),8), old_t(8), sampled(size(state,1))
    real(dp) :: primitive(reactive_nprim(size(species))), tg, sound, guess, dt, elapsed, volume
    real(dp) :: invtau,heat,masscoef,d,mdot,tdot,speed,transfer(5),increase(size(state,1)), recovered
    real(dp) :: sampled_pressure,cp,cv,h,e,entropy,xsat,liquid_rho,liquid_cp
    integer :: shape3(3),indices(3,8),p,n,i,j,k,attempt,total_steps
    logical :: local_ok, accepted
    impact_settings=spray_wall_options()
    if (present(wall_options)) impact_settings=wall_options
    faces='periodic'
    if (present(boundary_faces)) faces=boundary_faces
    ledger=0
    if (present(boundary_change)) boundary_change=0
    substeps=0; ok=.false.; message='Invalid spray coupling input'
    force_ledger=0
    if(options%fixed_parts.and..not.present(force_change)) return
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
      if (.not.options%khrt_children) then
        if(parcels(p)%khrt_shed_volume/=0) return
        ! Initial count also belongs to accepted splash children; only pending
        ! shed volume requires the KH child-history mode.
      endif
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
    if(options%etab_parcel_factor>=0) then
      allocate(etab_tangent(size(drops)),etab_last_dt(size(drops)),etab_origins(3,size(drops)))
      etab_tangent=0;etab_last_dt=0
      do p=1,size(drops)
        etab_origins(:,p)=drops(p)%position
      enddo
    endif
    p=1
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
        call liquid_transfer_properties(material,drops(p)%temperature,film%pressure,film%liquid_molecular_weight, &
          liquid_rho,liquid_cp,local_ok)
        if(.not.local_ok) return
        d=parcel_diameter(drops(p),material,film%pressure,film%liquid_molecular_weight); mdot=masscoef*d
        tdot=heat*(tg-drops(p)%temperature)
        if(options%abramzon_sirignano.and.options%evaporation) tdot=heat*max(0._dp,tg-drops(p)%temperature)
        tdot=(tdot-mdot*(film%vapor_enthalpy-liquid_energy(material,drops(p)%temperature)))/(drops(p)%mass*liquid_cp)
        dt=interval-elapsed
        if (invtau > 0.0_dp) dt=min(dt,options%max_fractional_change/invtau)
        if (mdot > 0.0_dp .and. d > 2.0_dp*options%minimum_diameter) &
          dt=min(dt,options%max_fractional_change*drops(p)%mass/mdot)
        if (abs(tdot) > 0.0_dp) &
          dt=min(dt,min(2.0_dp,options%max_fractional_change*drops(p)%temperature)/abs(tdot))
        speed=max(maxval(abs(drops(p)%velocity)),maxval(abs(film%gas_velocity)))
        if (speed > 0.0_dp.and..not.options%fixed_parts) dt=min(dt,options%particle_cfl*minval(spacing)/speed)
        accepted=.false.
        do attempt=1,30
          if (dt <= 0.0_dp .or. elapsed+dt <= elapsed) exit
          call advance_parcel_frozen_film(drops(p),material,film,options,dt,next,transfer,local_ok,force_change=parcel_force)
          if (local_ok .and. next%mass>0.0_dp) then
            call nasa7_mass_properties(species(material%vapor_index),next%temperature,cp,cv,h,e,entropy,local_ok)
            if (local_ok .and. options%evaporation) then
              call saturation_mole_fraction(material,species(material%vapor_index)%molecular_weight, &
                next%temperature,sampled_pressure,xsat,local_ok,vapor_enthalpy=h)
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
        force_ledger=force_ledger+parcel_force
        if(.not.all(ieee_is_finite(force_ledger))) return
        breakup_material=material
        if(next%mass>0.and.material%boiling_model=='pelemp-cap') then
          call liquid_transfer_properties(material,next%temperature,film%pressure,film%liquid_molecular_weight, &
            liquid_rho,liquid_cp,local_ok)
          if(.not.local_ok) return
          breakup_material%density=liquid_rho;breakup_material%heat_capacity=liquid_cp
          breakup_material%reference_energy=liquid_energy(material,next%temperature)
          breakup_material%reference_temperature=next%temperature
          breakup_material%density_model='constant';breakup_material%density_coefficients=0
          breakup_material%density_temperature_range=0;breakup_material%heat_capacity_model='constant'
          breakup_material%heat_capacity_coefficients=0;breakup_material%heat_capacity_temperature_range=0
          breakup_material%boiling_model='reject';breakup_material%critical_temperature=0
          breakup_material%normal_boiling_temperature=0
        endif
        if(next%mass>0) then
          call apply_mixture_breakup_viscosity([material],[1._dp],next%temperature,local_breakup,local_ok, &
            film%pressure,[film%liquid_molecular_weight])
          if(.not.local_ok) return
        endif
        release=spray_khrt_release();hit=.false.
        if (options%khrt_children) then
          ! A crossed wall owns the incident mass before any new KH release.
          call process_wall_event(accepted)
          if (.not.accepted) return
        end if
        if (options%khrt_children.and..not.hit) then
          shedding=spray_khrt_shedding_state(next%khrt_initial_count,next%khrt_shed_volume)
          reynolds=film%density*norm2(film%gas_velocity-next%velocity)*parcel_diameter(next,breakup_material)/film%viscosity
          call advance_khrt_shedding(next,breakup_material,local_breakup,sampled(1),film%gas_velocity,reynolds,dt, &
            options%khrt_reference_diameter_cubed,shedding,release,local_ok)
          if (local_ok) then
            next%khrt_initial_count=shedding%initial_count;next%khrt_shed_volume=shedding%shed_volume
          end if
        else if (.not.options%khrt_children) then
          if(options%etab_parcel_factor>=0) then
            call advance_spray_breakup(next,breakup_material,local_breakup,sampled(1),film%gas_velocity, &
              film%density,film%viscosity,dt,local_ok,etab_delta)
          else
            call advance_spray_breakup(next,breakup_material,local_breakup,sampled(1),film%gas_velocity, &
              film%density,film%viscosity,dt,local_ok)
          endif
        end if
        if (.not.local_ok) then
          message='Spray breakup rejected; complete coupling step rolled back'; return
        end if
        if (release%mass>0) then
          release_axis=next%velocity
          if (norm2(release_axis)<1e-12_dp.and.local_breakup%khrt_parcel_factor<0) &
            release_axis=film%gas_velocity-next%velocity
          birth_parent=next;birth_parent%position=drops(p)%position
          call construct_configured_khrt_children(local_breakup,birth_parent,release,release_axis,drops(size(drops))%id, &
            options%maximum_parcels-size(drops),spawned,child_exchange,local_ok,drift_dt=dt)
          if (.not.local_ok) return
          do n=1,8
            i=indices(1,n);j=indices(2,n);k=indices(3,n)
            candidate(2:5,i,j,k)=candidate(2:5,i,j,k)+weights(n)*child_exchange/volume
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
        exchange=0;local_ok=.true.
        if(.not.options%fixed_parts) call apply_parcel_boundaries(next,material,lower, &
          lower+spacing*real(shape3,dp),faces,exchange,local_ok)
        if (.not.local_ok) return
        ledger(1:5)=ledger(1:5)+exchange
        ledger(5+material%vapor_index)=ledger(5+material%vapor_index)+exchange(1)
        if(options%etab_parcel_factor>=0) etab_origins(:,p)=drops(p)%position
        drops(p)=next
        if(options%etab_parcel_factor>=0) then
          etab_tangent(p)=etab_tangent(p)+etab_delta
          etab_last_dt(p)=dt
        endif
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
    if(options%etab_parcel_factor>=0) then
      call split_reference_etab_end(local_ok)
    else
      call split_spray_parcels(drops,options%maximum_multiplicity,options%maximum_parcels,local_ok)
    endif
    if (.not.local_ok) then
      message='Parcel splitting exceeds storage or ID bounds; complete step rolled back';return
    end if
    state=candidate; temperature=temp; parcels=drops
    if (present(wall_films)) call move_alloc(patches,wall_films)
    if (present(boundary_change)) boundary_change=ledger
    if(present(force_change)) force_change=force_ledger
    substeps=total_steps; ok=.true.; message=''
  contains
    subroutine split_reference_etab_end(split_ok)
      logical, intent(out) :: split_ok
      type(spray_parcel), allocatable :: replacement(:),children(:)
      type(spray_liquid) :: child_liquid
      integer :: parent,child,number,kept,total,cursor,cell,ii,jj,kk
      integer(int64) :: watermark
      real(dp) :: compensation(4),boundary_exchange(5),cell_weights(8)
      integer :: cell_indices(3,8)
      real(dp) :: cell_primitive(reactive_nprim(size(species))),cell_temperature,cell_sound
      logical :: eligible,valid
      split_ok=.false.;total=0;kept=0
      if(size(etab_tangent)/=size(drops).or.size(etab_last_dt)/=size(drops)) return
      do parent=1,size(drops)
        eligible=drops(parent)%mass>0.and.drops(parent)%multiplicity>options%maximum_multiplicity
        if(.not.eligible) then
          kept=kept+1;total=total+1;cycle
        endif
        number=0
        call khrt_reference_child_count(drops(parent)%multiplicity,options%etab_parcel_factor, &
          options%maximum_parcels-total,number,valid)
        if(.not.valid) then
          message='ETAB child count exceeds parcel capacity; complete step rolled back';return
        endif
        total=total+number
      enddo
      if(total==kept) then
        split_ok=.true.;return
      endif
      if(total>options%maximum_parcels) return
      watermark=drops(size(drops))%id
      if(watermark>huge(watermark)-int(total-kept,int64)) return
      allocate(replacement(total));cursor=0
      do parent=1,size(drops)
        if(drops(parent)%mass>0.and.drops(parent)%multiplicity>options%maximum_multiplicity) cycle
        cursor=cursor+1;replacement(cursor)=drops(parent)
      enddo
      do parent=1,size(drops)
        if(drops(parent)%mass<=0.or.drops(parent)%multiplicity<=options%maximum_multiplicity) cycle
        child_liquid=select_spray_liquid(drops(parent),liquid,materials)
        birth_parent=drops(parent);birth_parent%position=etab_origins(:,parent)
        call construct_configured_etab_children(birth_parent,drops(parent)%velocity, &
          options%etab_parcel_factor,etab_tangent(parent),watermark,total-cursor, &
          children,compensation,valid,etab_last_dt(parent))
        if(.not.valid) then
          message='ETAB child construction rejected; complete step rolled back';return
        endif
        call particle_boundary_stencil(etab_origins(:,parent),lower,spacing,shape3,faces, &
          cell_indices,cell_weights,valid)
        if(.not.valid) return
        do cell=1,8
          ii=cell_indices(1,cell);jj=cell_indices(2,cell);kk=cell_indices(3,cell)
          candidate(2:5,ii,jj,kk)=candidate(2:5,ii,jj,kk)+cell_weights(cell)*compensation/volume
          call reactive_conserved_to_primitive(species,candidate(:,ii,jj,kk),temp(ii,jj,kk), &
            cell_primitive,cell_temperature,cell_sound,valid)
          if(.not.valid) then
            message='ETAB child gas compensation rejected; complete step rolled back';return
          endif
          temp(ii,jj,kk)=cell_temperature
        enddo
        do child=1,size(children)
          call apply_parcel_boundaries(children(child),child_liquid,lower, &
            lower+spacing*real(shape3,dp),faces,boundary_exchange,valid)
          if(.not.valid) return
          ledger(1:5)=ledger(1:5)+boundary_exchange
          ledger(5+child_liquid%vapor_index)=ledger(5+child_liquid%vapor_index)+boundary_exchange(1)
          cursor=cursor+1;replacement(cursor)=children(child)
        enddo
        watermark=watermark+int(size(children),int64)
      enddo
      if(cursor/=total) return
      call move_alloc(replacement,drops)
      split_ok=.true.
    end subroutine

    subroutine process_wall_event(event_ok)
      logical, intent(out) :: event_ok
      event_ok=.false.;hit=.false.
      if (any(impact_settings%enabled).and.next%mass>0.and..not.options%fixed_parts) then
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
            if(options%etab_parcel_factor>=0) then
              etab_tangent=[etab_tangent,spread(0._dp,1,size(spawned))]
              etab_last_dt=[etab_last_dt,spread(0._dp,1,size(spawned))]
              etab_origins=reshape([reshape(etab_origins,[size(etab_origins)]), &
                (spawned(n)%position,n=1,size(spawned))],[3,size(etab_last_dt)])
            endif
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
      wall_options,wall_films,film_component_mass,acceleration,force_change)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    real(dp), optional, intent(in) :: acceleration(3)
    real(dp), optional, intent(inout) :: force_change(5)
    real(dp) :: accel(3),force_ledger(5),parcel_force(5)
    type(spray_wall_options), optional, intent(in) :: wall_options
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_component_mass(:,:)
    type(spray_wall_options) :: impact_settings
    type(spray_options) :: film_options,local_breakup
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
    real(dp), allocatable :: starts(:),etab_tangent(:),etab_last_dt(:),etab_origins(:,:)
    real(dp) :: reynolds,child_exchange(4),release_axis(3),etab_delta
    type(spray_parcel) :: next,birth_parent
    type(spray_liquid) :: effective
    type(spray_film) :: mixture_film
    real(dp) :: y(size(components)),rates(size(components)),hv(size(components)),next_mass(size(components))
    real(dp) :: before_boundary(size(components)),rho,cp,e,heat,invtau,tdot,elapsed,dt,speed,old_mass,volume
    real(dp) :: surface_energy
    real(dp) :: sampled(size(state,1)),primitive(size(state,1)),old_cells(size(state,1),8),old_t(8)
    real(dp) :: change(size(state,1)),ledger(size(state,1)),weights(8),exchange(5),tg,recovered,sound,guess
    integer :: shape3(3),indices(3,8),p,n,i,j,k,attempt,count,old_count,parent,v
    character(len=32) :: faces(6)
    logical :: valid,accepted
    ok=.false.;substeps=0;message='Invalid mixture coupling state'
    accel=0;force_ledger=0
    if(present(acceleration)) accel=acceleration
    if(.not.all(ieee_is_finite(accel)).or.maxval(abs(accel))>1e6_dp) return
    if((options%fixed_parts.or.any(accel/=0)).and..not.present(force_change)) return
    if(any(accel/=0).and.options%mixture_time_scheme/='euler'.and..not.options%fixed_parts) return
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
      if (.not.options%khrt_children) then
        if(parcels(p)%khrt_shed_volume/=0) return
        ! Initial count also belongs to accepted splash children; only pending
        ! shed volume requires the KH child-history mode.
      endif
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
    if(options%etab_parcel_factor>=0) then
      allocate(etab_tangent(size(drops)),etab_last_dt(size(drops)),etab_origins(3,size(drops)))
      etab_tangent=0;etab_last_dt=0
      do p=1,size(drops)
        etab_origins(:,p)=drops(p)%position
      enddo
    endif
    p=1
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
        call liquid_mixture_properties(components,y,drops(p)%temperature,rho,cp,e,valid,primitive(5),mw)
        if (.not.valid) return
        dt=interval-elapsed
        if (invtau>0) dt=min(dt,options%max_fractional_change/invtau)
        tdot=heat*(tg-drops(p)%temperature)
        if(options%abramzon_sirignano.and.options%evaporation) tdot=heat*max(0._dp,tg-drops(p)%temperature)
        do n=1,size(components)
          if (rates(n)>0.and.y(n)>mixture_trace_fraction) &
            dt=min(dt,options%max_fractional_change*masses(n,p)/rates(n))
          if(rates(n)==0) cycle
          surface_energy=liquid_energy(components(n),drops(p)%temperature)
          if(options%mixture_time_scheme=='euler') then
            call liquid_evaporation_energy(components(n),drops(p)%temperature,surface_energy,valid, &
              primitive(5),mw(n))
            if(.not.valid) return
          endif
          tdot=tdot-rates(n)*(hv(n)-surface_energy)
        end do
        tdot=tdot/(drops(p)%mass*cp)
        if (abs(tdot)>0) dt=min(dt,min(2.0_dp,options%max_fractional_change*drops(p)%temperature)/abs(tdot))
        speed=max(maxval(abs(drops(p)%velocity)),maxval(abs(primitive(2:4))))
        if (speed>0.and..not.options%fixed_parts) dt=min(dt,options%particle_cfl*minval(spacing)/speed)
        accepted=.false.
        do attempt=1,30
          if (dt<=0.or.elapsed+dt<=elapsed) exit
          next=drops(p);next_mass=masses(:,p)
          call advance_mixture_frozen_gas(species,transport,primitive,tg,next,components,next_mass,options,dt,change,valid, &
            acceleration=accel,force_change=parcel_force,reynolds_number=reynolds)
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
        force_ledger=force_ledger+parcel_force
        if(.not.all(ieee_is_finite(force_ledger))) return
        ! Feedback recovery used primitive as scratch; recover the sampled gas.
        call reactive_conserved_to_primitive(species,sampled,guess,primitive,tg,sound,valid)
        if (.not.valid) return
        if (next%mass>0) then
          y=next_mass/sum(next_mass)
          call liquid_mixture_properties(components,y,next%temperature,rho,cp,e,valid, &
            primitive(5),mw)
          if (.not.valid) return
        end if
        effective=components(1);effective%density=rho;effective%density_model='constant'
        effective%density_coefficients=0;effective%density_temperature_range=0
        effective%heat_capacity=cp;effective%heat_capacity_model='constant'
        effective%heat_capacity_coefficients=0;effective%heat_capacity_temperature_range=0
        effective%reference_temperature=next%temperature;effective%reference_energy=e
        effective%boiling_model='reject';effective%critical_temperature=0
        effective%normal_boiling_temperature=0
        old_mass=next%mass
        ! Restore the sampled primitive: the feedback recovery above used it as scratch.
        call reactive_conserved_to_primitive(species,sampled,guess,primitive,tg,sound,valid)
        if (.not.valid) return
        hit=.false.
        if (any(impact_settings%enabled).and.next%mass>0.and..not.options%fixed_parts) then
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
              if(options%etab_parcel_factor>=0) then
                etab_tangent=[etab_tangent,spread(0._dp,1,size(children))]
                etab_last_dt=[etab_last_dt,spread(0._dp,1,size(children))]
                etab_origins=reshape([reshape(etab_origins,[size(etab_origins)]), &
                  (children(n)%position,n=1,size(children))],[3,size(etab_last_dt)])
              endif
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
          local_breakup=options
          call apply_mixture_breakup_viscosity(components,y,next%temperature,local_breakup,valid,primitive(5),mw)
          if(.not.valid) return
          call build_mixture_film(species,transport,primitive,tg,next,components,y,options,mixture_film,rates,hv,valid)
          if (.not.valid) return
          if (options%khrt_children) then
            shedding=spray_khrt_shedding_state(next%khrt_initial_count,next%khrt_shed_volume)
            call advance_khrt_shedding(next,effective,local_breakup,sampled(1),primitive(2:4),reynolds,dt, &
              options%khrt_reference_diameter_cubed,shedding,release,valid)
            if (valid) then
              next%khrt_initial_count=shedding%initial_count;next%khrt_shed_volume=shedding%shed_volume
            end if
          else
            if(options%etab_parcel_factor>=0) then
              call advance_spray_breakup(next,effective,local_breakup,sampled(1),primitive(2:4),mixture_film%density, &
                mixture_film%viscosity,dt,valid,etab_delta,reynolds_number=reynolds)
              if(valid) etab_tangent(p)=etab_tangent(p)+etab_delta
            else
              call advance_spray_breakup(next,effective,local_breakup,sampled(1),primitive(2:4),mixture_film%density, &
                mixture_film%viscosity,dt,valid,reynolds_number=reynolds)
            endif
          end if
          if (.not.valid) return
          next_mass=next_mass*(next%mass/old_mass)
        end if
        if (release%mass>0) then
          release_axis=next%velocity
          if (norm2(release_axis)<1e-12_dp.and.options%khrt_parcel_factor<0) &
            release_axis=primitive(2:4)-next%velocity
          birth_parent=next;birth_parent%position=drops(p)%position
          call construct_configured_khrt_children(options,birth_parent,release,release_axis,drops(size(drops))%id, &
            options%maximum_parcels-size(drops),children,child_exchange,valid,drift_dt=dt)
          if (.not.valid) return
          do n=1,8
            i=indices(1,n);j=indices(2,n);k=indices(3,n)
            gas(2:5,i,j,k)=gas(2:5,i,j,k)+weights(n)*child_exchange/volume
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
        exchange=0;valid=.true.
        if(.not.options%fixed_parts) call apply_parcel_boundaries(next,effective,lower, &
          lower+spacing*real(shape3,dp),faces,exchange,valid)
        if (.not.valid) return
        if (next%mass==0) next_mass=0
        ledger(1:5)=ledger(1:5)+exchange
        do n=1,size(components)
          v=components(n)%vapor_index
          ledger(5+v)=ledger(5+v)+next%multiplicity*(next_mass(n)-before_boundary(n))
        end do
        if(options%etab_parcel_factor>=0) etab_origins(:,p)=drops(p)%position
        drops(p)=next;masses(:,p)=next_mass;elapsed=min(interval,elapsed+dt);count=count+1
        if(options%etab_parcel_factor>=0) etab_last_dt(p)=dt
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
    if(options%etab_parcel_factor>=0) then
      call split_reference_mixture_etab(valid)
    else
      call split_spray_parcels(drops,options%maximum_multiplicity,options%maximum_parcels,valid)
    endif
    if (.not.valid) return
    if (options%etab_parcel_factor<0.and.size(drops)>old_count) then
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
    if(present(force_change)) force_change=force_ledger
    substeps=count;ok=.true.;message=''
  contains
    subroutine split_reference_mixture_etab(split_ok)
      logical, intent(out) :: split_ok
      type(spray_parcel), allocatable :: replacement(:),new_children(:)
      real(dp), allocatable :: replacement_mass(:,:)
      type(spray_liquid) :: child_liquid
      integer :: source,child,number,kept,total,cursor,cell,ii,jj,kk,component
      integer(int64) :: watermark
      real(dp) :: compensation(4),boundary_exchange(5),cell_weights(8), &
        cell_primitive(reactive_nprim(size(species))),cell_temperature,cell_sound
      integer :: cell_indices(3,8)
      logical :: eligible,ready
      split_ok=.false.;total=0;kept=0
      if(size(etab_tangent)/=size(drops).or.size(etab_last_dt)/=size(drops)) return
      do source=1,size(drops)
        eligible=drops(source)%mass>0.and.drops(source)%multiplicity>options%maximum_multiplicity
        if(.not.eligible) then
          kept=kept+1;total=total+1;cycle
        endif
        number=0
        call khrt_reference_child_count(drops(source)%multiplicity,options%etab_parcel_factor, &
          options%maximum_parcels-total,number,ready)
        if(.not.ready) then
          message='Mixed ETAB child count exceeds parcel capacity; complete step rolled back';return
        endif
        total=total+number
      enddo
      if(total==kept) then
        split_ok=.true.;return
      endif
      watermark=drops(size(drops))%id
      if(watermark>huge(watermark)-int(total-kept,int64)) return
      allocate(replacement(total),replacement_mass(size(components),total));cursor=0
      do source=1,size(drops)
        if(drops(source)%mass>0.and.drops(source)%multiplicity>options%maximum_multiplicity) cycle
        cursor=cursor+1;replacement(cursor)=drops(source);replacement_mass(:,cursor)=masses(:,source)
      enddo
      do source=1,size(drops)
        if(drops(source)%mass<=0.or.drops(source)%multiplicity<=options%maximum_multiplicity) cycle
        call liquid_mixture_properties(components,masses(:,source)/sum(masses(:,source)), &
          drops(source)%temperature,rho,cp,e,ready)
        if(.not.ready) return
        child_liquid=components(1);child_liquid%density=rho;child_liquid%density_model='constant'
        child_liquid%density_coefficients=0;child_liquid%density_temperature_range=0
        child_liquid%heat_capacity=cp;child_liquid%heat_capacity_model='constant'
        child_liquid%heat_capacity_coefficients=0;child_liquid%heat_capacity_temperature_range=0
        child_liquid%reference_temperature=drops(source)%temperature;child_liquid%reference_energy=e
        birth_parent=drops(source);birth_parent%position=etab_origins(:,source)
        call construct_configured_etab_children(birth_parent,drops(source)%velocity, &
          options%etab_parcel_factor,etab_tangent(source),watermark,total-cursor, &
          new_children,compensation,ready,etab_last_dt(source))
        if(.not.ready) then
          message='Mixed ETAB child construction rejected; complete step rolled back';return
        endif
        call particle_boundary_stencil(etab_origins(:,source),lower,spacing,shape3,faces, &
          cell_indices,cell_weights,ready)
        if(.not.ready) return
        do cell=1,8
          ii=cell_indices(1,cell);jj=cell_indices(2,cell);kk=cell_indices(3,cell)
          gas(2:5,ii,jj,kk)=gas(2:5,ii,jj,kk)+cell_weights(cell)*compensation/volume
          call reactive_conserved_to_primitive(species,gas(:,ii,jj,kk),temp(ii,jj,kk), &
            cell_primitive,cell_temperature,cell_sound,ready)
          if(.not.ready) then
            message='Mixed ETAB gas compensation rejected; complete step rolled back';return
          endif
          temp(ii,jj,kk)=cell_temperature
        enddo
        do child=1,size(new_children)
          call apply_parcel_boundaries(new_children(child),child_liquid,lower, &
            lower+spacing*real(shape3,dp),faces,boundary_exchange,ready)
          if(.not.ready) return
          ledger(1:5)=ledger(1:5)+boundary_exchange
          do component=1,size(components)
            ledger(5+components(component)%vapor_index)= &
              ledger(5+components(component)%vapor_index)+ &
              boundary_exchange(1)*masses(component,source)/drops(source)%mass
          enddo
          cursor=cursor+1;replacement(cursor)=new_children(child)
          replacement_mass(:,cursor)=masses(:,source)
          if(new_children(child)%mass<=0) replacement_mass(:,cursor)=0
        enddo
        watermark=watermark+int(size(new_children),int64)
      enddo
      if(cursor/=total) return
      call move_alloc(replacement,drops);call move_alloc(replacement_mass,masses)
      split_ok=.true.
    end subroutine
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
            call reactive_conserved_to_primitive(species,sampled,guess,primitive,tg,sound,valid)
            if (.not.valid) return
            call advance_mixture_wall_film_exchange(next,components,next_mass,rates*dt,hv,qg*dt,qw*dt, &
              wall_options%minimum_height,change,wall_energy,valid,deplete_traces=.true., &
              pressure=primitive(5),gas_molecular_weights=species(components%vapor_index)%molecular_weight, &
            euler_temperature=options%mixture_time_scheme=='euler')
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
