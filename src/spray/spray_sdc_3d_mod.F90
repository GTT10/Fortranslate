module spray_sdc_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_conserved_to_primitive
  use mixture_thermo_mod, only: mixture_pressure
  use spray_parcel_mod, only: spray_parcel,spray_liquid,spray_options,spray_film,valid_spray_options,parcel_diameter
  use spray_mixture_mod, only: advance_mixture_frozen_gas,valid_component_inventory,liquid_mixture_properties
  use spray_mixture_mod, only: build_mixture_film
  use spray_breakup_mod, only: advance_spray_breakup,split_spray_parcels, &
    khrt_reference_child_count,khrt_reference_uniforms,construct_reference_etab_children, &
    spray_khrt_release,spray_khrt_shedding_state,advance_khrt_shedding, &
    construct_reference_khrt_children,construct_khrt_children,advance_tab_no_birth_khrt
  use spray_coupling_3d_mod, only: particle_boundary_stencil,valid_spray_boundary_faces,apply_parcel_boundaries
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options,valid_wall_options, &
    valid_wall_component_inventory,wall_liquid_conductivity_at_temperature,mixture_wall_film_rates, &
    advance_mixture_wall_film_exchange,mixture_impact_at_planar_wall
  use spray_mixture_mod, only: mixture_trace_fraction,apply_mixture_breakup_viscosity
  implicit none
  private
  ! The outer interval owns this private stage state and discards it on failure.
  type, public :: spray_sdc_stage_3d
    integer :: next_stage=0,iterations=0,flow_subcycles=1
    logical :: has_watermark=.false.
    integer(int64) :: initial_watermark=0,midpoint_watermark=0
    real(dp) :: time=0,interval=0,acceleration(3)=0,native_flow_cfl=-1
    type(spray_options), allocatable :: breakup_properties(:)
    type(spray_parcel), allocatable :: midpoint(:)
    real(dp), allocatable :: midpoint_mass(:,:),old_rate(:,:,:,:)
    real(dp) :: wall_change(5)=0,force_change(5)=0
    real(dp), allocatable :: boundary_change(:)
  end type
  public :: advance_spray_sdc_stage_3d,kick_serial_sdc_films
contains
  ! Every film substep samples the immutable gas stage. Feedback is exported
  ! as a source rate; no intermediate film update changes the sampled gas.
  subroutine kick_serial_sdc_films(species,transport,components,options,settings,state,temperature, &
      lower,spacing,duration,faces,films,film_mass,rate,wall_change,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(spray_wall_options), intent(in) :: settings
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),lower(3),spacing(3),duration
    character(len=*), intent(in) :: faces(6)
    type(spray_wall_film), allocatable, intent(inout) :: films(:)
    real(dp), allocatable, intent(inout) :: film_mass(:,:)
    real(dp), intent(inout) :: rate(:,:,:,:),wall_change(5)
    logical, intent(out) :: ok
    type(spray_wall_film), allocatable :: patches(:)
    type(spray_wall_film) :: trial
    real(dp), allocatable :: masses(:,:),rhs(:,:,:,:),kept(:,:)
    real(dp) :: wall(5),primitive(size(state,1)),cell(size(state,1)),weights(8),tg,tc,sound,volume
    real(dp) :: rates(size(components)),hv(size(components)),trial_mass(size(components)),qg,qw,tdot
    real(dp) :: conductivity,elapsed,dt,change(size(state,1)),wall_energy
    integer :: cells(3),indices(3,8),index(3),p,c,i,j,k,axis,face,total_steps,attempt
    logical :: valid,accepted
    ok=.false.
    if(.not.allocated(films).or..not.allocated(film_mass)) return
    if(.not.valid_spray_options(options).or..not.valid_wall_options(settings,faces)) return
    if(.not.all(ieee_is_finite([lower,spacing,duration]))) return
    if(any(spacing<=0).or.duration<=0.or.duration>1e6_dp) return
    cells=shape(temperature)
    if(any(cells<2).or.any(shape(state)/=[5+size(species),cells]).or.any(shape(rate)/=shape(state))) return
    if(size(components)<1.or.size(components)>16.or.size(films)>settings%maximum_films) return
    if(any(shape(film_mass)/=[size(components),size(films)])) return
    do p=1,size(films)
      if(.not.valid_wall_component_inventory(films(p),film_mass(:,p))) return
    enddo
    volume=product(spacing)
    if(.not.ieee_is_finite(volume).or.volume<=0) return
    patches=films;masses=film_mass;allocate(rhs,mold=rate);rhs=0;wall=0;total_steps=0
    if(settings%evolve) then
      do p=1,size(patches)
        axis=maxloc(abs(patches(p)%normal),dim=1)
        face=2*axis-merge(1,0,patches(p)%normal(axis)>0)
        if(.not.settings%enabled(face)) return
        call particle_boundary_stencil(patches(p)%position,lower,spacing,cells,faces,indices,weights,valid)
        if(.not.valid) return
        primitive=0;tg=0
        do c=1,8
          if(weights(c)<=0) cycle
          i=indices(1,c);j=indices(2,c);k=indices(3,c)
          call reactive_conserved_to_primitive(species,state(:,i,j,k),temperature(i,j,k),cell,tc,sound,valid)
          if(.not.valid) return
          primitive=primitive+weights(c)*cell;tg=tg+weights(c)*tc
        enddo
        primitive(5)=mixture_pressure(species,primitive(6:),primitive(1),tg,valid)
        if(.not.valid) return
        index=max(1,min(cells,floor((patches(p)%position-lower)/spacing)+1))
        elapsed=0
        do while(elapsed<duration.and.patches(p)%mass>0)
          if(total_steps>=options%maximum_substeps) return
          conductivity=wall_liquid_conductivity_at_temperature(settings,patches(p)%material_id,patches(p)%temperature)
          call mixture_wall_film_rates(patches(p),components,masses(:,p),species,transport,primitive,tg,options, &
            conductivity,settings%temperature(face),rates,hv,qg,qw,tdot,valid)
          if(.not.valid) return
          dt=duration-elapsed
          do c=1,size(components)
            if(rates(c)>0.and.masses(c,p)>mixture_trace_fraction*patches(p)%mass) &
              dt=min(dt,options%max_fractional_change*masses(c,p)/rates(c))
          enddo
          if(abs(tdot)>0) dt=min(dt,min(2._dp,options%max_fractional_change*patches(p)%temperature)/abs(tdot))
          accepted=.false.
          do attempt=1,30
            if(dt<=0.or.elapsed+dt<=elapsed) exit
            trial=patches(p);trial_mass=masses(:,p)
            call advance_mixture_wall_film_exchange(trial,components,trial_mass,rates*dt,hv,qg*dt,qw*dt, &
              settings%minimum_height,change,wall_energy,valid,deplete_traces=.true.,pressure=primitive(5), &
              gas_molecular_weights=species(components%vapor_index)%molecular_weight, &
              euler_temperature=options%mixture_time_scheme=='euler')
            if(valid.and.trial%mass>0) valid= &
              wall_liquid_conductivity_at_temperature(settings,trial%material_id,trial%temperature)>0
            if(valid) then
              accepted=.true.;exit
            endif
            dt=.5_dp*dt
          enddo
          if(.not.accepted) return
          patches(p)=trial;masses(:,p)=trial_mass;wall(5)=wall(5)-wall_energy
          rhs(:,index(1),index(2),index(3))=rhs(:,index(1),index(2),index(3))+change/(duration*volume)
          elapsed=min(duration,elapsed+dt);total_steps=total_steps+1
        enddo
      enddo
    endif
    if(.not.all(ieee_is_finite(rhs)).or..not.all(ieee_is_finite(wall))) return
    allocate(kept(size(components),count(patches%mass>0)))
    do c=1,size(components)
      kept(c,:)=pack(masses(c,:),patches%mass>0)
    enddo
    patches=pack(patches,patches%mass>0)
    call move_alloc(patches,films);call move_alloc(kept,film_mass)
    rate=rhs;wall_change=wall;ok=.true.
  end subroutine kick_serial_sdc_films

  subroutine advance_spray_sdc_stage_3d(species,transport,components,options,state,temperature,lower,spacing, &
      time,interval,stage,iterations,faces,parcels,component_mass,work,rate,wall_change,boundary_change,force_change,ok, &
      acceleration,last_id,breakup_properties,wall_options,wall_films,film_component_mass,native_flow_cfl)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    ! Only surface tension and the liquid viscosity model vary by material.
    type(spray_options), optional, intent(in) :: breakup_properties(:)
    type(spray_options) :: parcel_options
    type(spray_wall_options), optional, intent(in) :: wall_options
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_component_mass(:,:)
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: patch_mass(:,:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),lower(3),spacing(3),time,interval
    integer, intent(in) :: stage,iterations
    character(len=*), intent(in) :: faces(6)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    type(spray_sdc_stage_3d), intent(inout) :: work
    real(dp), intent(inout) :: rate(:,:,:,:),wall_change(5),boundary_change(:),force_change(5)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: native_flow_cfl
    real(dp), optional, intent(in) :: acceleration(3)
    integer(int64), optional, intent(inout) :: last_id
    integer(int64) :: watermark
    type(spray_sdc_stage_3d) :: next
    type(spray_parcel), allocatable :: drops(:),kept(:)
    real(dp), allocatable :: masses(:,:),kept_mass(:,:),rhs(:,:,:,:),tangent(:)
    type(spray_khrt_release), allocatable :: releases(:)
    real(dp), allocatable :: release_fractions(:,:),sample_velocity(:,:)
    real(dp) :: accel(3),coupled_accel(3),wall(5),exterior(size(boundary_change)),forced(5),parcel_force(5)
    real(dp) :: cell_primitive(size(state,1)),primitive(size(state,1)),change(size(state,1)),weights(8)
    real(dp) :: before(5),after(5),boundary_delta(5),upper(3),half,volume,cell_temperature,tg,sound
    real(dp) :: flow_cfl,sub_interval,sub_half,kick_position(3),tangent_delta,transfer_reynolds
    integer :: cells(3),indices(3,8),source_index(3),p,c,n,i,j,k,remaining,cycle,cycles
    logical :: valid,hit
    ok=.false.;accel=0;coupled_accel=0;wall=0;exterior=0;forced=0
    flow_cfl=-1
    if(present(native_flow_cfl)) flow_cfl=native_flow_cfl
    if(.not.ieee_is_finite(flow_cfl)) return
    if(flow_cfl< -1) return
    if(present(acceleration)) accel=acceleration
    if(.not.all(ieee_is_finite([time,interval,lower,spacing,accel]))) return
    if(time<0.or.interval<=0.or.interval>1e6_dp.or.any(spacing<=0).or.any(abs(accel)>1e6_dp)) return
    if(iterations<1.or.stage<0.or.stage>iterations.or..not.valid_spray_options(options)) return
    if((options%maximum_multiplicity>0.or.options%khrt_children).and..not.present(last_id)) return
    if(.not.valid_spray_boundary_faces(faces)) return
    cells=shape(temperature);n=5+size(species)
    if(any(cells<2).or.any(shape(state)/=[n,cells]).or.any(shape(rate)/=[n,cells])) return
    if(size(boundary_change)/=n.or.size(components)<1) return
    if(.not.allocated(parcels).or..not.allocated(component_mass)) return
    if(any(shape(component_mass)/=[size(components),size(parcels)])) return
    if(size(parcels)>options%maximum_parcels) return
    watermark=0
    if(present(last_id)) then
      watermark=last_id
      if(watermark<0.or.any(parcels%id>watermark)) return
    endif
    if(present(breakup_properties)) then
      if(size(breakup_properties)<1) return
      do p=1,size(breakup_properties)
        if(.not.valid_spray_options(breakup_properties(p))) return
      enddo
      if(any(parcels%material_id<1).or.any(parcels%material_id>size(breakup_properties))) return
    endif
    if(present(wall_options).neqv.present(wall_films)) return
    if(present(wall_films).neqv.present(film_component_mass)) return
    if(present(wall_films)) then
      if(.not.allocated(wall_films).or..not.allocated(film_component_mass)) return
      if(.not.valid_wall_options(wall_options,faces)) return
      if(any(shape(film_component_mass)/=[size(components),size(wall_films)])) return
      if(size(wall_films)>wall_options%maximum_films) return
      do p=1,size(wall_films)
        if(.not.valid_wall_component_inventory(wall_films(p),film_component_mass(:,p))) return
      enddo
      if(any(wall_options%enabled).and.size(parcels)>0.and..not.present(last_id)) return
      patches=wall_films;patch_mass=film_component_mass
    endif
    if(stage==0) then
      if(work%next_stage/=0) return
    else
      if(work%next_stage/=stage.or.work%iterations/=iterations.or.work%time/=time.or.work%interval/=interval) return
      if(any(work%acceleration/=accel).or.work%native_flow_cfl/=flow_cfl) return
      if(present(breakup_properties).neqv.allocated(work%breakup_properties)) return
      if(present(breakup_properties)) then
        if(size(breakup_properties)/=size(work%breakup_properties)) return
        do p=1,size(breakup_properties)
          if(.not.same_properties(breakup_properties(p),work%breakup_properties(p))) return
        enddo
      endif
      if(present(last_id).neqv.work%has_watermark) return
      if(present(last_id)) then
        if(last_id/=work%initial_watermark) return
      endif
      if(.not.allocated(work%midpoint).or..not.allocated(work%midpoint_mass).or..not.allocated(work%old_rate)) return
    endif
    if(stage>0.and.stage<iterations) then
      rate=work%old_rate;wall_change=work%wall_change;boundary_change=work%boundary_change
      force_change=work%force_change;work%next_stage=stage+1;ok=.true.;return
    endif
    half=.5_dp*interval;volume=product(spacing);upper=lower+spacing*real(cells,dp)
    if(.not.ieee_is_finite(volume).or.volume<=0.or.half<=0) return
    if(stage==0) then
      drops=parcels;masses=component_mass
    else
      drops=work%midpoint;masses=work%midpoint_mass
      watermark=work%midpoint_watermark
    endif
    next=work;allocate(rhs,mold=rate);rhs=0
    allocate(tangent(size(drops)),source=0._dp)
    allocate(releases(size(drops)))
    allocate(release_fractions(size(components),size(drops)),source=0._dp)
    allocate(sample_velocity(3,size(drops)),source=0._dp)
    ! Pinned PeleC uses the default -1 CFL, hence one native flow iteration.
    ! Preserve the explicit PelePhysics subdivision hook for library callers.
    cycles=1
    if(stage==0.and.flow_cfl>.5_dp) then
      if(flow_cfl>.5_dp*real(min(options%maximum_substeps,huge(cycles)-1),dp)) return
      cycles=ceiling(flow_cfl/.5_dp)
    endif
    sub_interval=interval/real(cycles,dp);sub_half=.5_dp*sub_interval
    if(options%mixture_time_scheme=='euler'.or.options%fixed_parts) coupled_accel=accel
    if(stage==0.and.options%mixture_time_scheme/='euler'.and..not.options%fixed_parts) then
      call kick_force();if(.not.valid) return
    endif
    do cycle=1,cycles
    do p=1,size(drops)
      if(.not.valid_component_inventory(drops(p),masses(:,p))) return
      if(drops(p)%mass<=0) cycle
      parcel_options=options
      if(present(breakup_properties)) then
        if(drops(p)%material_id<1.or.drops(p)%material_id>size(breakup_properties)) return
        call select_properties(breakup_properties(drops(p)%material_id),parcel_options)
      endif
      kick_position=drops(p)%position
      call particle_boundary_stencil(drops(p)%position,lower,spacing,cells,faces,indices,weights,valid)
      if(.not.valid) return
      call source_cell(drops(p)%position,source_index,valid);if(.not.valid) return
      primitive=0;tg=0
      do c=1,8
        if(weights(c)<=0) cycle
        i=indices(1,c);j=indices(2,c);k=indices(3,c)
        call reactive_conserved_to_primitive(species,state(:,i,j,k),temperature(i,j,k), &
          cell_primitive,cell_temperature,sound,valid)
        if(.not.valid) return
        ! InterpolateGasPhase converts each cell to rho, velocity, Y and T
        ! before weighting. Recovering EOS from weighted conserved energy
        ! would add unresolved velocity variance to the sampled temperature.
        primitive=primitive+weights(c)*cell_primitive
        tg=tg+weights(c)*cell_temperature
      enddo
      primitive(5)=mixture_pressure(species,primitive(6:),primitive(1),tg,valid)
      if(.not.valid) return
      if(options%mixture_time_scheme=='euler') then
        call advance_mixture_frozen_gas(species,transport,primitive,tg,drops(p),components,masses(:,p), &
          parcel_options,sub_interval, &
          change,valid,time_modifier=.5_dp,acceleration=coupled_accel,force_change=parcel_force, &
          reynolds_number=transfer_reynolds)
      else
        call advance_mixture_frozen_gas(species,transport,primitive,tg,drops(p),components,masses(:,p), &
          parcel_options,sub_half, &
          change,valid,acceleration=coupled_accel,force_change=parcel_force,reynolds_number=transfer_reynolds)
      endif
      if(.not.valid) return
      forced=forced+parcel_force
      if(drops(p)%mass>0.and.options%breakup_model/='none'.and. &
          (options%breakup_model=='etab'.or.cycle==cycles)) then
        block
          type(spray_liquid) :: effective
          type(spray_film) :: film
          type(spray_khrt_shedding_state) :: shedding
          real(dp) :: density,cp,energy,old_mass,fractions(size(components)),breakup_dt
          real(dp) :: transfer(size(components)),enthalpy(size(components))
          old_mass=drops(p)%mass;fractions=masses(:,p)/old_mass
          call liquid_mixture_properties(components,fractions,drops(p)%temperature,density,cp,energy,valid, &
            primitive(5),species(components%vapor_index)%molecular_weight)
          if(.not.valid) return
          effective=components(1);effective%density=density;effective%density_model='constant'
          effective%density_coefficients=0;effective%density_temperature_range=0
          effective%heat_capacity=cp;effective%heat_capacity_model='constant'
          effective%heat_capacity_coefficients=0;effective%heat_capacity_temperature_range=0
          effective%reference_temperature=drops(p)%temperature;effective%reference_energy=energy
          effective%boiling_model='reject';effective%critical_temperature=0;effective%normal_boiling_temperature=0
          call build_mixture_film(species,transport,primitive,tg,drops(p),components,fractions,parcel_options, &
            film,transfer,enthalpy,valid)
          if(.not.valid) return
          call apply_mixture_breakup_viscosity(components,fractions,drops(p)%temperature,parcel_options,valid, &
            primitive(5),species(components%vapor_index)%molecular_weight)
          if(.not.valid) return
          ! Pinned SprayParticles.cpp passes sub_dt to TAB, and explicitly
          ! dtmod*flow_dt to KHRT. A thermal half-kick does not halve TAB's clock.
          breakup_dt=half
          if(options%breakup_model=='etab') breakup_dt=sub_interval
          if(options%khrt_children) then
            shedding=spray_khrt_shedding_state(drops(p)%khrt_initial_count,drops(p)%khrt_shed_volume)
            release_fractions(:,p)=fractions
            sample_velocity(:,p)=primitive(2:4)
            call advance_khrt_shedding(drops(p),effective,parcel_options,primitive(1),primitive(2:4),transfer_reynolds,half, &
              options%khrt_reference_diameter_cubed,shedding,releases(p),valid,make_new_drops=stage==0)
            if(valid) then
              drops(p)%khrt_initial_count=shedding%initial_count;drops(p)%khrt_shed_volume=shedding%shed_volume
            endif
          else if(options%etab_parcel_factor>=0) then
            call advance_spray_breakup(drops(p),effective,parcel_options,primitive(1),primitive(2:4), &
              film%density,film%viscosity,breakup_dt,valid,tangent_delta,reynolds_number=transfer_reynolds)
            if(valid) tangent(p)=tangent(p)+tangent_delta
            if(valid.and.stage==iterations) &
              call advance_tab_no_birth_khrt(drops(p),effective,parcel_options,primitive(1),primitive(2:4), &
                transfer_reynolds,half,valid)
          else
            call advance_spray_breakup(drops(p),effective,parcel_options,primitive(1),primitive(2:4), &
              film%density,film%viscosity,breakup_dt,valid,reynolds_number=transfer_reynolds)
          endif
          if(.not.valid) return
          masses(:,p)=masses(:,p)*(drops(p)%mass/old_mass)
        end block
      endif
      ! Transfer helper motion is replaced by MKD's flow-subcycle drift.
      drops(p)%position=kick_position
      ! SprayParticles.cpp samples with trilinear weights but deposits into
      ! ijkc, the cell containing the particle before this kick's drift.
      i=source_index(1);j=source_index(2);k=source_index(3)
      rhs(:,i,j,k)=rhs(:,i,j,k)+change/(half*volume)
    enddo
    ! Only MKD creates replacement particles. Both retained particles and new
    ! children drift from their source interpolation position afterwards.
    if(stage==0.and.cycle==cycles.and.options%etab_parcel_factor>=0) then
      call replace_etab();if(.not.valid) return
    endif
    if(stage==0.and.cycle==cycles.and.options%khrt_children) then
      call append_khrt();if(.not.valid) return
    endif
    if(stage==0.and.cycle==cycles) then
      ! Breakup replacement can change population size before wall flight.
      ! Its release/tangent payload has already been consumed at this point.
      deallocate(tangent,releases,release_fractions,sample_velocity)
      allocate(tangent(size(drops)),source=0._dp)
      allocate(releases(size(drops)))
      allocate(release_fractions(size(components),size(drops)),source=0._dp)
      allocate(sample_velocity(3,size(drops)),source=0._dp)
    endif
    do p=1,size(drops)
      if(stage/=0.or.options%fixed_parts.or.drops(p)%mass<=0) cycle
      call inventory(drops(p),masses(:,p),before,valid);if(.not.valid) return
      kick_position=drops(p)%position
      drops(p)%position=drops(p)%position+sub_interval*drops(p)%velocity
      hit=.false.
      if(present(wall_options)) then
        if(any(wall_options%enabled)) then
          call planar_impact(p,kick_position,hit)
          if(.not.valid) return
          if(hit) cycle
        endif
      endif
      ! Boundary handling needs mixture energy at the current liquid state.
      block
        type(spray_liquid) :: effective
        real(dp) :: density,cp,energy
        call liquid_mixture_properties(components,masses(:,p)/drops(p)%mass,drops(p)%temperature, &
          density,cp,energy,valid)
        if(.not.valid) return
        effective=components(1);effective%reference_temperature=drops(p)%temperature;effective%reference_energy=energy
        call apply_parcel_boundaries(drops(p),effective,lower,upper,faces,boundary_delta,valid)
      end block
      if(.not.valid) return
      if(drops(p)%mass==0) then
        exterior(:5)=exterior(:5)-before
        do c=1,size(components)
          exterior(5+components(c)%vapor_index)=exterior(5+components(c)%vapor_index)- &
            masses(c,p)*drops(p)%multiplicity
        enddo
        masses(:,p)=0
      else
        call inventory(drops(p),masses(:,p),after,valid);if(.not.valid) return
        wall=wall+after-before
      endif
    enddo
    enddo
    if(stage==iterations.and.options%mixture_time_scheme/='euler'.and..not.options%fixed_parts) then
      call kick_force();if(.not.valid) return
    endif
    if(stage==0) then
      next%midpoint=drops;next%midpoint_mass=masses;next%old_rate=rhs
      next%time=time;next%interval=interval;next%iterations=iterations;next%acceleration=accel
      next%flow_subcycles=cycles;next%native_flow_cfl=flow_cfl
      if(allocated(next%breakup_properties)) deallocate(next%breakup_properties)
      if(present(breakup_properties)) next%breakup_properties=breakup_properties
      next%has_watermark=present(last_id)
      if(present(last_id)) next%initial_watermark=last_id
      next%midpoint_watermark=watermark
    else
      wall=wall+work%wall_change;exterior=exterior+work%boundary_change;forced=forced+work%force_change
    endif
    if(.not.all(ieee_is_finite(rhs)).or..not.all(ieee_is_finite([wall,exterior,forced]))) return
    next%wall_change=wall;next%boundary_change=exterior;next%force_change=forced;next%next_stage=stage+1
    if(stage==iterations) then
      remaining=count(drops%mass>0);allocate(kept(remaining),kept_mass(size(components),remaining));c=0
      do p=1,size(drops)
        if(drops(p)%mass<=0) cycle
        c=c+1;kept(c)=drops(p);kept_mass(:,c)=masses(:,p)
      enddo
      if(options%maximum_multiplicity>0.and.options%etab_parcel_factor<0) then
        call split_spray_parcels(kept,options%maximum_multiplicity,options%maximum_parcels,valid,watermark,kept_mass)
        if(.not.valid) return
      endif
      call move_alloc(kept,parcels);call move_alloc(kept_mass,component_mass);next%next_stage=0
      if(present(last_id)) last_id=watermark
    endif
    if(present(wall_films)) then
      call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_component_mass)
    endif
    rate=rhs;wall_change=wall;boundary_change=exterior;force_change=forced;work=next;ok=.true.
  contains
    subroutine planar_impact(parent,previous,crossed)
      integer, intent(in) :: parent
      real(dp), intent(in) :: previous(3)
      logical, intent(out) :: crossed
      type(spray_parcel), allocatable :: children(:)
      type(spray_wall_film) :: deposit
      type(spray_options) :: properties
      real(dp), allocatable :: child_mass(:,:)
      real(dp) :: deposited_mass(size(components)),budget(5),q(size(state,1)),cell(size(state,1)),tc,tg,cs
      real(dp) :: w(8)
      integer :: ijk(3,8),c,i,j,k
      valid=.false.;crossed=.false.
      if(watermark>huge(watermark)-4_int64) return
      call particle_boundary_stencil(previous,lower,spacing,cells,faces,ijk,w,valid)
      if(.not.valid) return
      q=0;tg=0
      do c=1,8
        if(w(c)<=0) cycle
        i=ijk(1,c);j=ijk(2,c);k=ijk(3,c)
        call reactive_conserved_to_primitive(species,state(:,i,j,k),temperature(i,j,k),cell,tc,cs,valid)
        if(.not.valid) return
        q=q+w(c)*cell;tg=tg+w(c)*tc
      enddo
      q(5)=mixture_pressure(species,q(6:),q(1),tg,valid)
      if(.not.valid) return
      properties=options
      if(present(breakup_properties)) then
        if(drops(parent)%material_id<1.or.drops(parent)%material_id>size(breakup_properties)) then
          valid=.false.;return
        endif
        call select_properties(breakup_properties(drops(parent)%material_id),properties)
      endif
      call mixture_impact_at_planar_wall(drops(parent),previous,components,masses(:,parent),properties, &
        species(components%vapor_index)%molecular_weight,q(5),lower,upper,wall_options,patches,watermark+1_int64, &
        children,child_mass,deposit,deposited_mass,budget,crossed,valid)
      if(.not.valid.or..not.crossed) return
      if(size(children)>options%maximum_parcels-size(drops)) then
        valid=.false.;return
      endif
      if(size(children)>0) then
        masses=reshape([reshape(masses,[size(masses)]),reshape(child_mass,[size(child_mass)])], &
          [size(components),size(drops)+size(children)])
        tangent=[tangent,spread(0._dp,1,size(children))]
        releases=[releases,spread(spray_khrt_release(),1,size(children))]
        release_fractions=reshape([reshape(release_fractions,[size(release_fractions)]), &
          spread(0._dp,1,size(components)*size(children))],[size(components),size(drops)+size(children)])
        sample_velocity=reshape([reshape(sample_velocity,[size(sample_velocity)]), &
          spread(0._dp,1,3*size(children))],[3,size(drops)+size(children)])
        drops=[drops,children];watermark=watermark+int(size(children),int64)
      endif
      if(deposit%mass>0) then
        if(size(patches)>=wall_options%maximum_films) then
          valid=.false.;return
        endif
        patch_mass=reshape([reshape(patch_mass,[size(patch_mass)]),deposited_mass], &
          [size(components),size(patches)+1])
        patches=[patches,deposit]
      endif
      wall=wall-budget
    end subroutine planar_impact
    subroutine select_properties(properties,selected)
      type(spray_options), intent(in) :: properties
      type(spray_options), intent(inout) :: selected
      selected%surface_tension=properties%surface_tension
      selected%liquid_viscosity=properties%liquid_viscosity
      selected%liquid_viscosity_model=properties%liquid_viscosity_model
      selected%liquid_viscosity_coefficients=properties%liquid_viscosity_coefficients
      selected%liquid_viscosity_temperature_range=properties%liquid_viscosity_temperature_range
    end subroutine
    logical function same_properties(a,b)
      type(spray_options), intent(in) :: a,b
      same_properties=a%surface_tension==b%surface_tension.and.a%liquid_viscosity==b%liquid_viscosity.and. &
        a%liquid_viscosity_model==b%liquid_viscosity_model.and. &
        all(a%liquid_viscosity_coefficients==b%liquid_viscosity_coefficients).and. &
        all(a%liquid_viscosity_temperature_range==b%liquid_viscosity_temperature_range)
    end function
    subroutine source_cell(position,index,good)
      real(dp), intent(in) :: position(3)
      integer, intent(out) :: index(3)
      logical, intent(out) :: good
      real(dp) :: coordinate(3)
      integer :: axis
      good=.false.;index=0
      if(.not.all(ieee_is_finite(position))) return
      do axis=1,3
        if(faces(2*axis)=='periodic') then
          coordinate(axis)=modulo(position(axis)-lower(axis),upper(axis)-lower(axis))/spacing(axis)
        else
          if(position(axis)<lower(axis).or.position(axis)>upper(axis)) return
          coordinate(axis)=(position(axis)-lower(axis))/spacing(axis)
        endif
      enddo
      index=max(1,min(cells,floor(coordinate)+1));good=.true.
    end subroutine source_cell
    subroutine append_khrt()
      type(spray_parcel), allocatable :: replacement(:),children(:)
      real(dp), allocatable :: replacement_mass(:,:),uniforms(:)
      integer, allocatable :: numbers(:)
      integer :: parent,child,total,cursor,index(3)
      integer(int64) :: base_id,next_id
      real(dp) :: compensation(4),axis(3)
      valid=.false.;total=size(drops)
      allocate(numbers(size(drops)),source=0)
      do parent=1,size(drops)
        if(releases(parent)%mass<=0) cycle
        numbers(parent)=4
        if(options%khrt_parcel_factor>=0) then
          call khrt_reference_child_count(releases(parent)%multiplicity,options%khrt_parcel_factor, &
            options%maximum_parcels-total,numbers(parent),valid)
          if(.not.valid) return
        endif
        if(numbers(parent)>options%maximum_parcels-total) then
          valid=.false.;return
        endif
        total=total+numbers(parent)
      enddo
      if(watermark>huge(watermark)-int(total-size(drops),int64)) then
        valid=.false.;return
      endif
      allocate(replacement(total),replacement_mass(size(components),total));cursor=size(drops)
      replacement(:cursor)=drops;replacement_mass(:,:cursor)=masses
      base_id=watermark;next_id=watermark
      do parent=1,size(drops)
        if(numbers(parent)==0) cycle
        axis=drops(parent)%velocity;compensation=0
        if(options%khrt_parcel_factor>=0) then
          allocate(uniforms(numbers(parent)))
          call khrt_reference_uniforms(drops(parent)%id,base_id,uniforms)
          call construct_reference_khrt_children(drops(parent),releases(parent),axis,options%khrt_parcel_factor, &
            uniforms,next_id+1_int64,children,compensation,valid,drift_dt=0._dp)
          deallocate(uniforms)
        else
          if(norm2(axis)<1e-12_dp) axis=sample_velocity(:,parent)-axis
          call construct_khrt_children(drops(parent),releases(parent),axis,next_id+1_int64, &
            children,compensation(4),valid,drift_dt=0._dp)
        endif
        if(.not.valid) return
        call source_cell(drops(parent)%position,index,valid)
        if(.not.valid) return
        rhs(2:5,index(1),index(2),index(3))=rhs(2:5,index(1),index(2),index(3))+compensation/(half*volume)
        do child=1,size(children)
          cursor=cursor+1;replacement(cursor)=children(child)
          replacement_mass(:,cursor)=children(child)%mass*release_fractions(:,parent)
        enddo
        next_id=next_id+int(size(children),int64)
      enddo
      valid=cursor==total.and.all(ieee_is_finite(rhs))
      if(.not.valid) return
      call move_alloc(replacement,drops);call move_alloc(replacement_mass,masses);watermark=next_id
    end subroutine append_khrt
    subroutine replace_etab()
      type(spray_parcel), allocatable :: replacement(:),children(:)
      real(dp), allocatable :: replacement_mass(:,:),uniforms(:)
      integer, allocatable :: numbers(:)
      integer :: parent,child,total,births,cursor,index(3)
      integer(int64) :: base_id,next_id
      real(dp) :: compensation(4)
      valid=.false.;total=size(drops);births=0
      allocate(numbers(size(drops)),source=0)
      do parent=1,size(drops)
        if(drops(parent)%mass<=0.or.drops(parent)%multiplicity<=options%maximum_multiplicity) cycle
        call khrt_reference_child_count(drops(parent)%multiplicity,options%etab_parcel_factor, &
          options%maximum_parcels-total+1,numbers(parent),valid)
        if(.not.valid) return
        total=total-1+numbers(parent);births=births+numbers(parent)
      enddo
      if(births==0) then
        valid=.true.;return
      endif
      if(watermark>huge(watermark)-int(births,int64)) then
        valid=.false.;return
      endif
      allocate(replacement(total),replacement_mass(size(components),total));cursor=0
      do parent=1,size(drops)
        if(numbers(parent)>0) cycle
        cursor=cursor+1;replacement(cursor)=drops(parent);replacement_mass(:,cursor)=masses(:,parent)
      enddo
      base_id=watermark;next_id=watermark
      do parent=1,size(drops)
        if(numbers(parent)==0) cycle
        allocate(uniforms(numbers(parent)))
        call khrt_reference_uniforms(drops(parent)%id,base_id,uniforms)
        call construct_reference_etab_children(drops(parent),drops(parent)%velocity,options%etab_parcel_factor, &
          tangent(parent),uniforms,next_id+1_int64,children,compensation,valid,drift_dt=0._dp)
        deallocate(uniforms)
        if(.not.valid) return
        call source_cell(drops(parent)%position,index,valid)
        if(.not.valid) return
        rhs(2:5,index(1),index(2),index(3))=rhs(2:5,index(1),index(2),index(3))+compensation/(half*volume)
        do child=1,size(children)
          cursor=cursor+1;replacement(cursor)=children(child);replacement_mass(:,cursor)=masses(:,parent)
        enddo
        next_id=next_id+int(size(children),int64)
      enddo
      valid=cursor==total.and.all(ieee_is_finite(rhs))
      if(.not.valid) return
      call move_alloc(replacement,drops);call move_alloc(replacement_mass,masses);watermark=next_id
    end subroutine replace_etab
    subroutine inventory(drop,mass,total,good)
      type(spray_parcel), intent(in) :: drop
      real(dp), intent(in) :: mass(:)
      real(dp), intent(out) :: total(5)
      logical, intent(out) :: good
      real(dp) :: density,cp,energy
      call liquid_mixture_properties(components,mass/drop%mass,drop%temperature,density,cp,energy,good)
      if(.not.good) return
      total=drop%mass*drop%multiplicity*[1._dp,drop%velocity,energy+.5_dp*sum(drop%velocity**2)]
      good=all(ieee_is_finite(total))
    end subroutine
    subroutine kick_force()
      real(dp) :: old(5),new(5)
      valid=.true.
      if(all(accel==0)) return
      do p=1,size(drops)
        if(drops(p)%mass<=0) cycle
        call inventory(drops(p),masses(:,p),old,valid);if(.not.valid) return
        drops(p)%velocity=drops(p)%velocity+half*accel
        if(any(abs(drops(p)%velocity)>1e7_dp)) then
          valid=.false.;return
        endif
        call inventory(drops(p),masses(:,p),new,valid);if(.not.valid) return
        forced=forced+new-old
      enddo
    end subroutine
  end subroutine
end module spray_sdc_3d_mod
