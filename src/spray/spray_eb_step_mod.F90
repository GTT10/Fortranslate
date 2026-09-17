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
  implicit none
  private
  public :: advance_eb_mixture_coupling
contains
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
    type(spray_parcel), allocatable :: drops(:),trial(:)
    type(spray_liquid) :: effective
    type(spray_film) :: film
    real(dp), allocatable :: gas(:,:,:,:),temp(:,:,:),next_gas(:,:,:,:),next_temp(:,:,:),field(:,:,:,:)
    real(dp), allocatable :: masses(:,:),next_mass(:,:),positions(:,:),samples(:,:),primitive(:,:),tg(:),change(:,:)
    real(dp) :: elapsed,dt,sound,invtau,heat,rho,cp,e,tdot,speed,old_mass,h(3),lower(3),upper(3),exchange(5),ledger(5)
    real(dp) :: y(size(components)),rates(size(components)),hv(size(components))
    integer :: np,nc,nvar,p,n,attempt,count_steps,hits
    logical :: valid,accepted
    ok=.false.;substeps=0
    if (.not.g%is_valid().or..not.valid_spray_options(options)) return
    if (.not.ieee_is_finite(interval)) return
    if (interval<0.or.interval>1e6_dp.or.options%khrt_children) return
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
        exchange=0;hits=0
        if (valid) call drift_spray_level_set(phi,lower,upper,dt,32,trial,exchange,hits,valid)
        if (valid) then
          next_gas=gas;next_temp=temp
          call deposit_eb_spray_gas(species,g,phi,positions,change,next_gas,next_temp,valid)
        end if
        if (valid) then
          accepted=.true.;exit
        end if
        dt=.5_dp*dt
      end do
      if (.not.accepted) return
      gas=next_gas;temp=next_temp;drops=trial;masses=next_mass;ledger=ledger+exchange
      elapsed=elapsed+dt;count_steps=count_steps+1
    end do
    state=gas;temperature=temp;parcels=drops;component_mass=masses
    boundary_change=ledger;substeps=count_steps;ok=.true.
  end subroutine
end module
