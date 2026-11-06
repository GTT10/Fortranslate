module spray_eb_film_step_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar,reactive_conserved_to_primitive
  use spray_parcel_mod, only: spray_liquid,spray_options,valid_spray_options,valid_spray_liquid
  use spray_mixture_mod, only: mixture_trace_fraction
  use spray_wall_impact_mod, only: spray_wall_film,valid_wall_component_inventory, &
    mixture_wall_film_rates,advance_mixture_wall_film_exchange
  use spray_eb_impact_mod, only: spray_eb_wall_options
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use spray_eb_intersection_mod, only: intersect_spray_level_set
  use spray_eb_field_mod, only: sample_eb_spray_field,deposit_eb_spray_gas
  implicit none
  private
  public :: advance_eb_film_coupling
contains
  subroutine advance_eb_film_coupling(species,transport,components,options,walls,g,phi,state,temperature, &
      films,film_mass,interval,substeps,wall_change,ok,periodic)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(spray_eb_wall_options), intent(in) :: walls
    type(eb_geometry_3d), intent(in) :: g
    real(dp), intent(in) :: phi(0:,0:,0:),interval
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:),wall_change(5)
    type(spray_wall_film), allocatable, intent(inout) :: films(:)
    real(dp), allocatable, intent(inout) :: film_mass(:,:)
    integer, intent(out) :: substeps
    logical, intent(out) :: ok
    logical, optional, intent(in) :: periodic(3)
    type(spray_wall_film), allocatable :: patches(:),trial(:),kept(:)
    real(dp), allocatable :: gas(:,:,:,:),temp(:,:,:),next_gas(:,:,:,:),next_temp(:,:,:),field(:,:,:,:)
    real(dp), allocatable :: masses(:,:),next_mass(:,:),kept_mass(:,:),positions(:,:),samples(:,:),change(:,:)
    real(dp), allocatable :: rates(:,:),hv(:,:),qg(:),qw(:)
    real(dp) :: primitive(reactive_nvar(size(species))),tg,sound,tdot,elapsed,dt,wall_energy,ledger,next_ledger
    real(dp) :: lower(3),upper(3),offset,at,point(3),normal(3)
    integer :: np,nc,nvar,p,c,attempt,count_steps,cell(3)
    logical :: valid,hit,accepted
    ok=.false.;substeps=0
    if (.not.allocated(films).or..not.allocated(film_mass)) return
    if (.not.g%is_valid().or..not.valid_spray_options(options)) return
    if (.not.all(ieee_is_finite([interval,walls%temperature,walls%conductivity,walls%minimum_height]))) return
    if (interval<0.or.interval>1e6_dp.or.walls%temperature<1.or.walls%temperature>1e5_dp) return
    if (walls%conductivity<=0.or.walls%conductivity>1e4_dp) return
    if (walls%minimum_height<1e-15_dp.or.walls%minimum_height>1) return
    np=size(films);nc=size(components);nvar=reactive_nvar(size(species))
    if (nc<1.or.nc>16.or.np>min(125000,walls%maximum_films)) return
    if (size(transport)/=size(species).or.any(shape(film_mass)/=[nc,np])) return
    if (any(shape(state)/=[nvar,g%nx,g%ny,g%nz])) return
    if (any(shape(temperature)/=[g%nx,g%ny,g%nz])) return
    do c=1,nc
      if (.not.valid_spray_liquid(components(c))) return
      if (components(c)%vapor_index>size(species)) return
    end do
    lower=[g%x_lower,g%y_lower,g%z_lower];upper=[g%x_upper,g%y_upper,g%z_upper]
    offset=max(1e-8_dp*min(g%dx,g%dy,g%dz),128*epsilon(1._dp)*maxval(abs([lower,upper])))
    allocate(positions(3,np),change(nvar,np),rates(nc,np),hv(nc,np),qg(np),qw(np))
    allocate(field(nvar+1,g%nx,g%ny,g%nz))
    do p=1,np
      if (.not.valid_wall_component_inventory(films(p),film_mass(:,p))) return
      positions(:,p)=films(p)%position+offset*films(p)%normal
      hit=.false.;at=0;point=0;normal=0;cell=0
      call intersect_spray_level_set(phi,lower,upper,positions(:,p), &
        films(p)%position-offset*films(p)%normal,hit,at,point,normal,cell,valid)
      if (.not.valid.or..not.hit) return
      if (norm2(point-films(p)%position)>offset*.01_dp) return
      if (dot_product(normal,films(p)%normal)<1-1e-10_dp) return
    end do
    gas=state;temp=temperature;patches=films;masses=film_mass;change=0
    call deposit_eb_spray_gas(species,g,phi,positions,change,gas,temp,valid,periodic)
    if (.not.valid) return
    elapsed=0;ledger=0;count_steps=0
    do while(walls%evolve.and.elapsed<interval)
      if (count_steps>=options%maximum_substeps) return
      field(:nvar,:,:,:)=gas;field(nvar+1,:,:,:)=temp
      call sample_eb_spray_field(g,phi,field,positions,samples,valid,periodic)
      if (.not.valid) return
      dt=interval-elapsed;rates=0;hv=0;qg=0;qw=0
      do p=1,np
        if (patches(p)%mass<=0) cycle
        call reactive_conserved_to_primitive(species,samples(:nvar,p),samples(nvar+1,p),primitive,tg,sound,valid)
        if (.not.valid) return
        call mixture_wall_film_rates(patches(p),components,masses(:,p),species,transport,primitive,tg,options, &
          walls%conductivity,walls%temperature,rates(:,p),hv(:,p),qg(p),qw(p),tdot,valid)
        if (.not.valid) return
        do c=1,nc
          if (rates(c,p)>0.and.masses(c,p)>mixture_trace_fraction*patches(p)%mass) &
            dt=min(dt,options%max_fractional_change*masses(c,p)/rates(c,p))
        end do
        if (abs(tdot)>0) dt=min(dt,min(2._dp,options%max_fractional_change*patches(p)%temperature)/abs(tdot))
      end do
      accepted=.false.
      do attempt=1,30
        if (dt<=0.or.elapsed+dt<=elapsed) exit
        trial=patches;next_mass=masses;change=0;next_ledger=0;valid=.true.
        do p=1,np
          if (patches(p)%mass<=0) cycle
          call reactive_conserved_to_primitive(species,samples(:nvar,p),samples(nvar+1,p),primitive,tg,sound,valid)
          if (.not.valid) exit
          call advance_mixture_wall_film_exchange(trial(p),components,next_mass(:,p),rates(:,p)*dt,hv(:,p), &
            qg(p)*dt,qw(p)*dt,walls%minimum_height,change(:,p),wall_energy,valid,deplete_traces=.true., &
            pressure=primitive(5),gas_molecular_weights=species(components%vapor_index)%molecular_weight, &
            euler_temperature=options%mixture_time_scheme=='euler')
          if (.not.valid) exit
          next_ledger=next_ledger-wall_energy
        end do
        if (valid) then
          next_gas=gas;next_temp=temp
          call deposit_eb_spray_gas(species,g,phi,positions,change,next_gas,next_temp,valid,periodic)
        end if
        if (valid) then
          accepted=.true.;exit
        end if
        dt=.5_dp*dt
      end do
      if (.not.accepted) return
      gas=next_gas;temp=next_temp;patches=trial;masses=next_mass;ledger=ledger+next_ledger
      elapsed=elapsed+dt;count_steps=count_steps+1
    end do
    kept=pack(patches,patches%mass>0);allocate(kept_mass(nc,size(kept)))
    do c=1,nc
      kept_mass(c,:)=pack(masses(c,:),patches%mass>0)
    end do
    state=gas;temperature=temp
    call move_alloc(kept,films);call move_alloc(kept_mass,film_mass)
    wall_change=0;wall_change(5)=ledger;substeps=count_steps;ok=.true.
  end subroutine
end module
