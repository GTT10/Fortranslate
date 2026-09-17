module spray_eb_step_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar,reactive_conserved_to_primitive
  use spray_parcel_mod
  use spray_mixture_mod, only: valid_component_inventory,liquid_mixture_properties, &
    advance_mixture_frozen_gas,mixture_transfer_rates,mixture_trace_fraction,build_mixture_film
  use spray_breakup_mod, only: advance_spray_breakup
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use spray_eb_field_mod, only: sample_eb_spray_field,deposit_eb_spray_gas
  use spray_eb_drift_mod, only: drift_spray_level_set
  use, intrinsic :: iso_fortran_env, only: int64
  use spray_eb_impact_mod, only: advance_mixture_level_set_impacts,spray_eb_wall_options
  use spray_wall_impact_mod, only: spray_wall_film,valid_wall_component_inventory
  use spray_eb_film_step_mod, only: advance_eb_film_coupling
  implicit none
  private
  public :: advance_eb_mixture_coupling,advance_eb_mixture_coupling_allocated
contains
  ! Preserve the original fixed-population interface for specular callers.
  subroutine advance_eb_mixture_coupling(species,transport,components,options,g,phi,state,temperature, &
      parcels,component_mass,interval,substeps,boundary_change,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(eb_geometry_3d), intent(in) :: g
    real(dp), intent(in) :: phi(0:,0:,0:),interval
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:),component_mass(:,:),boundary_change(5)
    type(spray_parcel), intent(inout) :: parcels(:)
    integer, intent(out) :: substeps
    logical, intent(out) :: ok
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: masses(:,:)
    drops=parcels;masses=component_mass
    call advance_eb_mixture_coupling_allocated(species,transport,components,options,g,phi,state,temperature, &
      drops,masses,interval,substeps,boundary_change,ok)
    if (.not.ok) return
    parcels=drops;component_mass=masses
  end subroutine
  subroutine advance_eb_mixture_coupling_allocated(species,transport,components,options,g,phi,state,temperature, &
      parcels,component_mass,interval,substeps,boundary_change,ok,wall_options,films,film_mass,last_id,model_rejected)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(eb_geometry_3d), intent(in) :: g
    real(dp), intent(in) :: phi(0:,0:,0:),interval
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:),boundary_change(5)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    type(spray_eb_wall_options), optional, intent(in) :: wall_options
    type(spray_wall_film), allocatable, optional, intent(inout) :: films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    integer(int64), optional, intent(inout) :: last_id
    type(spray_wall_film), allocatable :: patches(:),next_patches(:)
    real(dp), allocatable :: patch_mass(:,:),next_patch_mass(:,:)
    real(dp) :: molecular_weights(size(components))
    integer(int64) :: watermark,next_watermark
    logical :: impacts,rejected
    logical, optional, intent(out) :: model_rejected
    integer, intent(out) :: substeps
    logical, intent(out) :: ok
    type(spray_parcel), allocatable :: drops(:),trial(:)
    type(spray_liquid) :: effective
    type(spray_film) :: film
    real(dp), allocatable :: gas(:,:,:,:),temp(:,:,:),next_gas(:,:,:,:),next_temp(:,:,:),field(:,:,:,:)
    real(dp), allocatable :: masses(:,:),next_mass(:,:),positions(:,:),samples(:,:),primitive(:,:),tg(:),change(:,:)
    real(dp) :: elapsed,dt,sound,invtau,heat,rho,cp,e,tdot,speed,old_mass,h(3),lower(3),upper(3),exchange(5),ledger(5)
    real(dp) :: y(size(components)),rates(size(components)),hv(size(components))
    integer :: np,nc,nvar,p,n,attempt,count_steps,hits,film_steps
    real(dp) :: film_change(5)
    logical :: valid,accepted
    ok=.false.;substeps=0
    if (present(model_rejected)) model_rejected=.false.
    if (.not.g%is_valid().or..not.valid_spray_options(options)) return
    if (.not.ieee_is_finite(interval)) return
    if (interval<0.or.interval>1e6_dp.or.options%khrt_children) return
    if (.not.allocated(parcels).or..not.allocated(component_mass)) return
    impacts=present(wall_options)
    if (present(films).neqv.impacts) return
    if (present(film_mass).neqv.impacts) return
    if (present(last_id).neqv.impacts) return
    if (impacts) then
      if (.not.all(ieee_is_finite([wall_options%temperature,wall_options%contact_angle]))) return
      if (wall_options%temperature<1.or.wall_options%temperature>1e5_dp) return
      if (wall_options%contact_angle<1e-6_dp.or.wall_options%contact_angle>acos(-1._dp)-1e-6_dp) return
      if (wall_options%evolve) then
        if (.not.all(ieee_is_finite([wall_options%conductivity,wall_options%minimum_height]))) return
        if (wall_options%conductivity<=0.or.wall_options%conductivity>1e4_dp) return
        if (wall_options%minimum_height<1e-15_dp.or.wall_options%minimum_height>1) return
      end if
      if (wall_options%maximum_hits<0.or.wall_options%maximum_hits>1024) return
      if (wall_options%maximum_films<0.or.wall_options%maximum_films>1000000) return
      if (last_id<0.or.last_id>huge(last_id)-4_int64) return
      do p=1,size(parcels)
        if (parcels(p)%id>last_id) return
        if (p>1) then
          if (parcels(p)%id<=parcels(p-1)%id) return
        end if
      end do
      if (.not.allocated(films).or..not.allocated(film_mass)) return
      if (size(films)>wall_options%maximum_films) return
      if (any(shape(film_mass)/=[size(components),size(films)])) return
      do p=1,size(films)
        if (.not.valid_wall_component_inventory(films(p),film_mass(:,p))) return
      end do
      patches=films;patch_mass=film_mass;watermark=last_id
    end if
    np=size(parcels);nc=size(components);nvar=reactive_nvar(size(species))
    if (nc<1.or.nc>16.or.np>min(125000,options%maximum_parcels)) return
    if (size(transport)/=size(species).or.any(shape(component_mass)/=[nc,np])) return
    if (any(shape(state)/=[nvar,g%nx,g%ny,g%nz])) return
    if (any(shape(temperature)/=[g%nx,g%ny,g%nz])) return
    do p=1,np
      if (.not.valid_component_inventory(parcels(p),component_mass(:,p))) return
    end do
    do n=1,nc
      if (.not.valid_spray_liquid(components(n))) return
      if (components(n)%vapor_index>size(species)) return
      molecular_weights(n)=species(components(n)%vapor_index)%molecular_weight
    end do
    lower=[g%x_lower,g%y_lower,g%z_lower];upper=[g%x_upper,g%y_upper,g%z_upper];h=[g%dx,g%dy,g%dz]
    gas=state;temp=temperature;drops=parcels;masses=component_mass;ledger=0
    allocate(positions(3,np),change(nvar,np),primitive(nvar,np),tg(np),field(nvar+1,g%nx,g%ny,g%nz))
    do p=1,np
      positions(:,p)=drops(p)%position
    end do
    change=0
    call deposit_eb_spray_gas(species,g,phi,positions,change,gas,temp,valid)
    if (.not.valid) return
    elapsed=0;count_steps=0
    do while(elapsed<interval)
      if (count_steps>=options%maximum_substeps) return
      np=size(drops)
      if (size(positions,2)/=np) then
        deallocate(positions,change,primitive,tg)
        allocate(positions(3,np),change(nvar,np),primitive(nvar,np),tg(np))
      end if
      primitive=0;primitive(5,:)=101325
      field(:nvar,:,:,:)=gas;field(nvar+1,:,:,:)=temp
      do p=1,np
        positions(:,p)=drops(p)%position
      end do
      call sample_eb_spray_field(g,phi,field,positions,samples,valid)
      if (.not.valid) return
      dt=interval-elapsed
      do p=1,np
        if (drops(p)%mass<=0) cycle
        call reactive_conserved_to_primitive(species,samples(:nvar,p),samples(nvar+1,p),primitive(:,p),tg(p),sound,valid)
        if (.not.valid) return
        y=masses(:,p)/drops(p)%mass
        call mixture_transfer_rates(species,transport,primitive(:,p),tg(p),drops(p),components,y,options, &
          invtau,heat,rates,hv,valid)
        if (.not.valid) exit
        call liquid_mixture_properties(components,y,drops(p)%temperature,rho,cp,e,valid)
        if (.not.valid) exit
        if (invtau>0) dt=min(dt,options%max_fractional_change/invtau)
        tdot=heat*(tg(p)-drops(p)%temperature)
        do n=1,nc
          if (rates(n)>0.and.y(n)>mixture_trace_fraction) &
            dt=min(dt,options%max_fractional_change*masses(n,p)/rates(n))
          tdot=tdot-rates(n)*(hv(n)-liquid_energy(components(n),drops(p)%temperature))
        end do
        tdot=tdot/(drops(p)%mass*cp)
        if (abs(tdot)>0) dt=min(dt,min(2.0_dp,options%max_fractional_change*drops(p)%temperature)/abs(tdot))
        speed=max(maxval(abs(drops(p)%velocity)),maxval(abs(primitive(2:4,p))))
        if (speed>0) dt=min(dt,options%particle_cfl*minval(h)/speed)
      end do
      if (.not.valid) return
      accepted=.false.
      do attempt=1,30
        if (dt<=0.or.elapsed+dt<=elapsed) exit
        trial=drops;next_mass=masses;change=0;valid=.true.
        do p=1,np
          if (drops(p)%mass<=0) cycle
          call advance_mixture_frozen_gas(species,transport,primitive(:,p),tg(p),trial(p),components, &
            next_mass(:,p),options,dt,change(:,p),valid)
          if (.not.valid) exit
          if (trial(p)%mass>0.and.options%breakup_model/='none') then
            old_mass=trial(p)%mass;y=next_mass(:,p)/old_mass
            call liquid_mixture_properties(components,y,trial(p)%temperature,rho,cp,e,valid)
            if (.not.valid) exit
            effective=components(1);effective%density=rho;effective%density_model='constant'
            effective%density_coefficients=0;effective%density_temperature_range=0
            effective%heat_capacity=cp;effective%heat_capacity_model='constant'
            effective%heat_capacity_coefficients=0;effective%heat_capacity_temperature_range=0
            effective%reference_temperature=trial(p)%temperature;effective%reference_energy=e
            effective%boiling_model='reject';effective%critical_temperature=0
            call build_mixture_film(species,transport,primitive(:,p),tg(p),trial(p),components,y,options, &
              film,rates,hv,valid)
            if (.not.valid) exit
            call advance_spray_breakup(trial(p),effective,options,primitive(1,p),primitive(2:4,p), &
              film%density,film%viscosity,dt,valid)
            if (.not.valid) exit
            next_mass(:,p)=next_mass(:,p)*(trial(p)%mass/old_mass)
          end if

          ! Source-first split: move at the updated velocity through EB geometry.
          ! Replace the source kernel's free-space trapezoidal position update.
          trial(p)%position=drops(p)%position
        end do
        exchange=0;hits=0;film_change=0
        if (valid) then
          next_gas=gas;next_temp=temp
          call deposit_eb_spray_gas(species,g,phi,positions,change,next_gas,next_temp,valid)
        end if
        if (valid) then
          if (impacts) then
            next_patches=patches;next_patch_mass=patch_mass;next_watermark=watermark
            if (wall_options%evolve) then
              call advance_eb_film_coupling(species,transport,components,options,wall_options,g,phi, &
                next_gas,next_temp,next_patches,next_patch_mass,dt,film_steps,film_change,valid)
            end if
            rejected=.false.
            if (valid) call advance_mixture_level_set_impacts(phi,lower,upper,dt,wall_options%maximum_hits, &
              min(125000,options%maximum_parcels),wall_options%maximum_films,components,options,molecular_weights, &
              101325._dp,wall_options%temperature,wall_options%contact_angle,trial,next_mass, &
              next_patches,next_patch_mass,next_watermark,exchange,hits,valid,primitive(5,:),rejected)
            if (rejected) then
              if (present(model_rejected)) model_rejected=.true.
              return
            end if
            exchange=-exchange
          else
            call drift_spray_level_set(phi,lower,upper,dt,32,trial,exchange,hits,valid)
          end if
        end if
        exchange=exchange+film_change
        if (valid) then
          accepted=.true.;exit
        end if
        dt=.5_dp*dt
      end do
      if (.not.accepted) return
      gas=next_gas;temp=next_temp;drops=trial;masses=next_mass;ledger=ledger+exchange
      if (impacts) then
        patches=next_patches;patch_mass=next_patch_mass;watermark=next_watermark
      end if
      elapsed=elapsed+dt;count_steps=count_steps+1
    end do
    state=gas;temperature=temp;parcels=drops;component_mass=masses
    if (impacts) then
      films=patches;film_mass=patch_mass;last_id=watermark
    end if
    boundary_change=ledger;substeps=count_steps;ok=.true.
  end subroutine
end module
