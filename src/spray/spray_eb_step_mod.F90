module spray_eb_step_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar,reactive_conserved_to_primitive
  use spray_parcel_mod
  use spray_mixture_mod, only: valid_component_inventory,liquid_mixture_properties, &
    advance_mixture_frozen_gas,mixture_transfer_rates,mixture_trace_fraction,build_mixture_film
  use spray_breakup_mod, only: advance_spray_breakup,advance_mixture_khrt_shedding
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use spray_eb_field_mod, only: sample_eb_spray_field,deposit_eb_spray_gas
  use spray_eb_intersection_mod, only: intersect_spray_level_set
  use spray_coupling_3d_mod, only: valid_spray_boundary_faces
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
      parcels,component_mass,interval,substeps,boundary_change,ok,wall_options,films,film_mass,last_id,model_rejected, &
      outer_faces,outer_change)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(eb_geometry_3d), intent(in) :: g
    real(dp), intent(in) :: phi(0:,0:,0:),interval
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:),boundary_change(5)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    character(len=*), optional, intent(in) :: outer_faces(6)
    real(dp), optional, intent(inout) :: outer_change(:)
    real(dp) :: exterior(reactive_nvar(size(species))),trial_exterior(reactive_nvar(size(species)))
    real(dp), allocatable :: liquid_energies(:),impact_exterior(:)
    real(dp) :: clipped_finish(3),preview_change(5),preview_exterior(5)
    type(spray_parcel) :: preview(1)
    logical :: preview_contact(1)
    integer :: preview_hits
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
    real(dp), allocatable :: reynolds(:),child_momentum(:,:),child_energy(:),impact_pressures(:)
    logical, allocatable :: eligible(:)
    real(dp) :: crossing_fraction,crossing_position(3),crossing_normal(3)
    integer :: crossing_cell(3),parent
    logical :: crossing,valid,accepted,periodic(3)
    ok=.false.;substeps=0
    if (present(model_rejected)) model_rejected=.false.
    if (.not.g%is_valid().or..not.valid_spray_options(options)) return
    if (.not.ieee_is_finite(interval)) return
    if (interval<0.or.interval>1e6_dp) return
    if (.not.allocated(parcels).or..not.allocated(component_mass)) return
    impacts=present(wall_options)
    if(present(outer_faces).neqv.present(outer_change)) return
    if(present(outer_faces)) then
      if(.not.valid_spray_boundary_faces(outer_faces)) return
      if(size(outer_change)/=reactive_nvar(size(species))) return
    end if
    periodic=.false.
    if(present(outer_faces)) periodic=outer_faces(1:6:2)=='periodic'
    exterior=0
    if (present(films).neqv.impacts) return
    if (present(film_mass).neqv.impacts) return
    if ((impacts.or.options%khrt_children).and..not.present(last_id)) return
    watermark=0
    if (present(last_id)) then
      watermark=last_id
      if(watermark<0) return
      do p=1,size(parcels)
        if(parcels(p)%id>watermark) return
        if(p>1) then
          if(parcels(p)%id<=parcels(p-1)%id) return
        end if
      end do
    end if
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
    if(impacts.and.present(outer_faces)) allocate(impact_exterior(5+nc))
    lower=[g%x_lower,g%y_lower,g%z_lower];upper=[g%x_upper,g%y_upper,g%z_upper];h=[g%dx,g%dy,g%dz]
    gas=state;temp=temperature;drops=parcels;masses=component_mass;ledger=0
    allocate(positions(3,np),change(nvar,np),primitive(nvar,np),tg(np),field(nvar+1,g%nx,g%ny,g%nz))
    do p=1,np
      positions(:,p)=drops(p)%position
    end do
    change=0
    call deposit_eb_spray_gas(species,g,phi,positions,change,gas,temp,valid,periodic)
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
      call sample_eb_spray_field(g,phi,field,positions,samples,valid,periodic)
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
      if(allocated(reynolds)) deallocate(reynolds,eligible)
      allocate(reynolds(np),eligible(np));reynolds=0;eligible=.true.
      accepted=.false.
      do attempt=1,30
        if (dt<=0.or.elapsed+dt<=elapsed) exit
        trial=drops;next_mass=masses;change=0;valid=.true.;next_watermark=watermark
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
            effective%normal_boiling_temperature=0
            call build_mixture_film(species,transport,primitive(:,p),tg(p),trial(p),components,y,options, &
              film,rates,hv,valid)
            if (.not.valid) exit
            if(options%khrt_children) then
              reynolds(p)=film%density*norm2(primitive(2:4,p)-trial(p)%velocity)* &
                parcel_diameter(trial(p),effective)/film%viscosity
            else
            call advance_spray_breakup(trial(p),effective,options,primitive(1,p),primitive(2:4,p), &
              film%density,film%viscosity,dt,valid)
            if (.not.valid) exit
            next_mass(:,p)=next_mass(:,p)*(trial(p)%mass/old_mass)
            end if
          end if

          ! Source-first split: move at the updated velocity through EB geometry.
          ! Replace the source kernel's free-space trapezoidal position update.
          trial(p)%position=drops(p)%position
        end do
        if(valid.and.options%khrt_children) then
          eligible=.true.
          do p=1,np
            if(trial(p)%mass<=0) cycle
            clipped_finish=trial(p)%position+dt*trial(p)%velocity
            if(present(outer_faces)) then
              ! Follow wrapped/reflecting flight before deciding on KH release.
              ! The preview remains private; actual parents and children drift later.
              preview(1)=trial(p);preview_contact=.false.
              preview_change=0;preview_exterior=0;preview_hits=0
              call drift_spray_level_set(phi,lower,upper,dt,32,preview,preview_change,preview_hits, &
                valid,outer_faces,[0._dp],preview_exterior,preview_contact)
              if(.not.valid) exit
              eligible(p)=.not.preview_contact(1)
              cycle
            end if
            call intersect_spray_level_set(phi,lower,upper,trial(p)%position, &
              clipped_finish,crossing,crossing_fraction, &
              crossing_position,crossing_normal,crossing_cell,valid)
            if(.not.valid) exit
            eligible(p)=.not.crossing
          end do
          if(valid) call advance_mixture_khrt_shedding(components,options,trial,next_mass,primitive(1,:), &
            primitive(2:4,:),reynolds,dt,next_watermark,child_energy,valid,eligible,gas_momentum=child_momentum)
          if(valid) valid=size(trial)<=125000
          if(valid) then
            change(2:4,:)=change(2:4,:)+child_momentum
            change(5,:)=change(5,:)+child_energy
          endif
        end if
        exchange=0;hits=0;film_change=0;trial_exterior=0
        if(allocated(impact_exterior)) impact_exterior=0
        if (valid) then
          next_gas=gas;next_temp=temp
          call deposit_eb_spray_gas(species,g,phi,positions,change,next_gas,next_temp,valid,periodic)
        end if
        if (valid) then
          if (impacts) then
            next_patches=patches;next_patch_mass=patch_mass
            if(allocated(impact_pressures)) deallocate(impact_pressures)
            allocate(impact_pressures(size(trial)));impact_pressures(:np)=primitive(5,:)
            parent=1
            do p=np+1,size(trial)
              do while(parent<np)
                if(trial(parent)%id>=trial(p)%parent_id) exit
                parent=parent+1
              end do
              if(trial(parent)%id/=trial(p)%parent_id) return
              impact_pressures(p)=primitive(5,parent)
            end do
            if (wall_options%evolve) then
              call advance_eb_film_coupling(species,transport,components,options,wall_options,g,phi, &
                next_gas,next_temp,next_patches,next_patch_mass,dt,film_steps,film_change,valid,periodic)
            end if
            rejected=.false.
            if (valid) call advance_mixture_level_set_impacts(phi,lower,upper,dt,wall_options%maximum_hits, &
              min(125000,options%maximum_parcels),wall_options%maximum_films,components,options,molecular_weights, &
              101325._dp,wall_options%temperature,wall_options%contact_angle,trial,next_mass, &
              next_patches,next_patch_mass,next_watermark,exchange,hits,valid,impact_pressures,rejected, &
              impact_exterior,outer_faces)
            if (rejected) then
              if (present(model_rejected)) model_rejected=.true.
              return
            end if
            if(valid.and.allocated(impact_exterior)) then
              trial_exterior(:5)=impact_exterior(:5)
              do n=1,nc
                trial_exterior(5+components(n)%vapor_index)=trial_exterior(5+components(n)%vapor_index)+ &
                  impact_exterior(5+n)
              end do
            end if
            exchange=-exchange
          else if(present(outer_faces)) then
            if(allocated(liquid_energies)) deallocate(liquid_energies)
            allocate(liquid_energies(size(trial)));liquid_energies=0
            do p=1,size(trial)
              if(trial(p)%mass<=0) cycle
              call liquid_mixture_properties(components,next_mass(:,p)/trial(p)%mass,trial(p)%temperature,rho,cp,e,valid)
              if(.not.valid) exit
              liquid_energies(p)=e
            end do
            if(valid) call drift_spray_level_set(phi,lower,upper,dt,32,trial,exchange,hits,valid, &
              outer_faces,liquid_energies,trial_exterior(:5))
            if(valid) then
              do p=1,size(trial)
                if(trial(p)%mass>0) cycle
                do n=1,nc
                  trial_exterior(5+components(n)%vapor_index)=trial_exterior(5+components(n)%vapor_index)- &
                    next_mass(n,p)*trial(p)%multiplicity
                end do
                next_mass(:,p)=0
              end do
            end if
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
      exterior=exterior+trial_exterior
      watermark=next_watermark
      if (impacts) then
        patches=next_patches;patch_mass=next_patch_mass
      end if
      elapsed=elapsed+dt;count_steps=count_steps+1
    end do
    state=gas;temperature=temp;parcels=drops;component_mass=masses
    if (impacts) then
      films=patches;film_mass=patch_mass
    end if
    if(present(last_id)) last_id=watermark
    if(present(outer_change)) outer_change=exterior
    boundary_change=ledger;substeps=count_steps;ok=.true.
  end subroutine
end module
