module spray_mpi_eb_impact_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use spray_parcel_mod
  use spray_mixture_mod, only: valid_component_inventory,liquid_mixture_properties
  use spray_wall_impact_mod, only: spray_wall_film,mixture_impact_boiling_temperature, &
    classify_wall_impact,splash_mass_fractions,splash_velocity_components
  use spray_eb_impact_mod, only: spray_eb_wall_options
  use spray_eb_intersection_mod, only: intersect_spray_level_set
  use spray_mpi_mod, only: sample_mpi_wall_film_thickness,resolve_mpi_wall_impacts
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d
  use spray_mpi_amr_tree_mod, only: migrate_eb_tree_spray_parcels
  use spray_mpi_amr_film_mod, only: migrate_mpi_amr_wall_films
  implicit none
  private
  public :: impact_mpi_level_set,advance_mpi_level_set_impacts,advance_mpi_eb_tree_impacts
contains
  ! One earliest EB contact per original parcel; endpoints are supplied by caller.
  ! Film queries use the pre-event distributed film population. Wall receives budget.
  ! Caller owns remaining flight, compaction and final particle/film migration.
  subroutine impact_mpi_level_set(comm,components,properties,settings,phi,molecular_weights,pressure, &
      previous,lower,upper,parcels,component_mass,films,film_mass,last_id,hit,fractions,wall_change,ok,model_rejected)
    type(MPI_Comm), intent(in) :: comm
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: properties
    type(spray_eb_wall_options), intent(in) :: settings
    real(dp), intent(in) :: phi(0:,0:,0:)
    real(dp), intent(in) :: molecular_weights(:),pressure(:),previous(:,:),lower(3),upper(3)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:),film_mass(:,:)
    type(spray_wall_film), allocatable, intent(inout) :: films(:)
    integer(int64), intent(inout) :: last_id
    logical, allocatable, intent(inout) :: hit(:)
    real(dp), allocatable, intent(inout) :: fractions(:)
    real(dp), intent(inout) :: wall_change(5)
    logical, optional, intent(out) :: model_rejected
    logical, intent(out) :: ok
    integer, allocatable :: indices(:)
    real(dp), allocatable :: location(:,:),normal(:,:),data(:,:),diameters(:),thickness(:)
    logical, allocatable :: crossed(:)
    real(dp), allocatable :: at(:),reference_phi(:,:,:)
    real(dp) :: budget(5),alpha,ms,dm,parts(4),vels(5),ut
    integer :: dims(3),reference_dims(3),cell(3),kind
    logical :: found,rejected,global_rejected
    real(dp) :: bounds(6),low(6),high(6),delta(3),contact_normal(3),fraction,rho,cp,e,boil,mu,un
    real(dp) :: y(size(components)),sigma
    integer :: np,nc,p,ne,n,ierr
    logical :: valid,agreed
    ok=.false.
    if(present(model_rejected)) model_rejected=.false.
    valid=allocated(parcels).and.allocated(component_mass).and.allocated(films).and.allocated(film_mass)
    valid=valid.and.valid_spray_options(properties)
    valid=valid.and.all(ieee_is_finite([lower,upper])).and.all(upper>lower)
    valid=valid.and.maxval(abs([lower,upper]))<=1e9_dp
    valid=valid.and.all(ieee_is_finite([settings%temperature,settings%contact_angle]))
    valid=valid.and.settings%temperature>=1.and.settings%temperature<=1e5_dp
    valid=valid.and.settings%contact_angle>=1e-6_dp.and.settings%contact_angle<=acos(-1._dp)-1e-6_dp
    dims=shape(phi)-1
    valid=valid.and.all(dims>=1).and.all(dims<=1000000)
    if(valid) valid=product(int(dims,int64))<=1000000_int64
    if(valid) valid=all((upper-lower)/real(dims,dp)>=1e-12_dp)
    valid=valid.and.all(ieee_is_finite(phi))
    if(valid) valid=maxval(abs(phi))<=1e100_dp
    call agree()
    if (.not.agreed) return
    reference_dims=dims
    call MPI_Bcast(reference_dims,3,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    valid=all(dims==reference_dims)
    call agree();if(.not.agreed) return
    reference_phi=phi
    call MPI_Bcast(reference_phi,size(reference_phi),MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    valid=all(phi==reference_phi)
    call agree();if(.not.agreed) return
    np=size(parcels);nc=size(components)
    valid=np<=125000.and.nc>=1.and.nc<=16.and.size(molecular_weights)==nc.and.size(pressure)==np
    valid=valid.and.all(shape(previous)==[3,np]).and.all(shape(component_mass)==[nc,np])
    call agree()
    if (.not.agreed) return
    valid=all(ieee_is_finite(previous)).and.all(ieee_is_finite(pressure))
    valid=valid.and.all(pressure>=1).and.all(pressure<=1e10_dp)
    do p=1,np
      valid=valid.and.valid_component_inventory(parcels(p),component_mass(:,p))
      valid=valid.and.all(previous(:,p)>=lower).and.all(previous(:,p)<=upper)
    end do
    call agree()
    if (.not.agreed) return
    bounds=[lower,upper]
    call MPI_Allreduce(bounds,low,6,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(bounds,high,6,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.any(low/=high)) return
    allocate(indices(np),location(3,np),normal(3,np),crossed(np),at(np))
    crossed=.false.;at=1;ne=0
    do p=1,np
      if (parcels(p)%mass<=0) cycle
      found=.false.;fraction=1;delta=0;contact_normal=0;cell=0
      call intersect_spray_level_set(phi,lower,upper,previous(:,p),parcels(p)%position, &
        found,fraction,delta,contact_normal,cell,valid)
      if(.not.valid) exit
      if(.not.found) cycle
      ne=ne+1;indices(ne)=p;crossed(p)=.true.;at(p)=fraction
      location(:,ne)=delta;normal(:,ne)=contact_normal
    end do
    call agree();if(.not.agreed) return
    allocate(data(6,ne),diameters(ne));data=0;sigma=properties%surface_tension
    do n=1,ne
      p=indices(n);y=component_mass(:,p)/parcels(p)%mass
      call liquid_mixture_properties(components,y,parcels(p)%temperature,rho,cp,e,valid)
      if (.not.valid) exit
      call mixture_impact_boiling_temperature(components,y,molecular_weights,pressure(p),boil,valid)
      if (.not.valid) exit
      mu=liquid_viscosity_at_temperature(properties,parcels(p)%temperature)
      valid=mu>0.and.sigma>=1e-8_dp.and.sigma<=1e3_dp
      if (.not.valid) exit
      diameters(n)=(6*parcels(p)%mass/(acos(-1.0_dp)*rho))**(1.0_dp/3.0_dp)
      un=abs(dot_product(parcels(p)%velocity,normal(:,n)))
      data(5,n)=rho*diameters(n)*un**2/sigma;data(6,n)=un*diameters(n)*rho/mu
      data(1,n)=sqrt(data(5,n)*sqrt(data(6,n)));data(2,n)=settings%temperature/boil
      data(3,n)=(26.788_dp*log(pressure(p)*1e-5_dp)+33.546_dp)* &
        exp(-.006_dp*(settings%temperature-260))+36
    end do
    call agree()
    if (.not.agreed) return
    call sample_mpi_wall_film_thickness(comm,films,location(:,:ne),normal(:,:ne), &
      1e-10_dp*maxval(upper-lower),thickness,valid)
    if (.not.valid) return
    data(4,:)=thickness/diameters
    rejected=.false.
    do n=1,ne
      p=indices(n);un=dot_product(parcels(p)%velocity,normal(:,n))
      valid=un<0.and.norm2(parcels(p)%velocity)>0
      if(.not.valid) exit
      alpha=max(acos(-1._dp)/6,asin(min(1._dp,abs(un)/norm2(parcels(p)%velocity))))
      call classify_wall_impact(data(1,n),data(2,n),alpha,data(3,n),data(4,n),kind,ms,dm,valid)
      if(valid.and.ms>0) then
        call splash_mass_fractions(alpha,ms,data(4,n),parts,valid)
        if(valid) then
          ut=norm2(parcels(p)%velocity-un*normal(:,n))
          call splash_velocity_components(un,ut,data(1,n),data(4,n),vels,valid)
          valid=valid.and.vels(5)>=0
        end if
      end if
      if(.not.valid) then
        rejected=.true.;exit
      end if
    end do
    call MPI_Allreduce(rejected,global_rejected,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    if(global_rejected) then
      if(present(model_rejected)) model_rejected=.true.
      return
    end if
    call agree();if(.not.agreed) return
    call resolve_mpi_wall_impacts(comm,components,parcels,component_mass,indices(:ne),location(:,:ne),normal(:,:ne),data, &
      settings%contact_angle,films,film_mass,properties%maximum_parcels,settings%maximum_films,last_id,budget,valid)
    if (.not.valid) return
    call move_alloc(crossed,hit);call move_alloc(at,fractions);wall_change=budget;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine

  ! Collective batches of remaining flight. Film thickness is frozen per batch.
  ! Keep source-rank ownership until the caller performs final migration.
  subroutine advance_mpi_level_set_impacts(comm,components,properties,settings,phi,molecular_weights,pressure, &
      lower,upper,dt,parcels,component_mass,films,film_mass,last_id,wall_change,hits,ok,model_rejected, outflow_change)
    type(MPI_Comm), intent(in) :: comm
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: properties
    type(spray_eb_wall_options), intent(in) :: settings
    real(dp), intent(in) :: phi(0:,0:,0:),molecular_weights(:),pressure(:),lower(3),upper(3),dt
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:),film_mass(:,:)
    type(spray_wall_film), allocatable, intent(inout) :: films(:)
    integer(int64), intent(inout) :: last_id
    real(dp), intent(inout) :: wall_change(5)
    integer, intent(inout) :: hits
    logical, intent(out) :: ok
    real(dp), optional, intent(inout) :: outflow_change(:)
    logical, optional, intent(out) :: model_rejected
    type(spray_parcel), allocatable :: drops(:),kept(:)
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: masses(:,:),patch_mass(:,:),remaining(:),pressures(:),previous(:,:),fractions(:)
    real(dp), allocatable :: next_remaining(:),next_pressure(:),kept_mass(:,:)
    real(dp), allocatable :: exit_fraction(:),exterior(:),global_exterior(:)
    logical, allocatable :: exits(:)
    real(dp) :: delta(3),finish(3),rho,cp,e,amount
    integer :: axis,mode,lowmode,highmode
    integer, allocatable :: depth(:),next_depth(:)
    logical, allocatable :: crossed(:)
    real(dp) :: ledger(5),budget(5),lowdt,highdt
    integer(int64) :: watermark,local_hits,total_hits
    integer :: np,p,q,lo,hi,mid,nc,ierr,limits(2),lowlimits(2),highlimits(2)
    logical :: valid,agreed,rejected,active,global_active
    ok=.false.
    if(present(model_rejected)) model_rejected=.false.
    valid=allocated(parcels).and.allocated(component_mass).and.allocated(films).and.allocated(film_mass)
    valid=valid.and.ieee_is_finite(dt)
    if(valid) valid=dt>=0.and.dt<=1e6_dp
    valid=valid.and.settings%maximum_hits>=0.and.settings%maximum_hits<=1024
    call agree();if(.not.agreed) return
    valid=size(pressure)==size(parcels)
    call agree();if(.not.agreed) return
    mode=merge(1,0,present(outflow_change))
    call MPI_Allreduce(mode,lowmode,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(mode,highmode,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.lowmode/=highmode) return
    if(present(outflow_change)) valid=size(outflow_change)==5+size(components)
    call agree();if(.not.agreed) return
    limits=[settings%maximum_hits,properties%maximum_parcels]
    call MPI_Allreduce(limits,lowlimits,2,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(limits,highlimits,2,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(lowlimits/=highlimits)) return
    call MPI_Allreduce(dt,lowdt,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(dt,highdt,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.lowdt/=highdt) return
    drops=parcels;masses=component_mass;patches=films;patch_mass=film_mass;watermark=last_id
    allocate(remaining(size(drops)),depth(size(drops)));remaining=dt;depth=0;pressures=pressure
    ledger=0;local_hits=0;nc=size(components)
    allocate(exterior(5+nc),global_exterior(5+nc));exterior=0
    do
      np=size(drops)
      if(allocated(previous)) deallocate(previous)
      if(allocated(exit_fraction)) deallocate(exit_fraction,exits)
      allocate(previous(3,np),exit_fraction(np),exits(np));exit_fraction=1;exits=.false.
      do p=1,np
        previous(:,p)=drops(p)%position
        if(drops(p)%mass<=0) cycle
        delta=remaining(p)*drops(p)%velocity;finish=previous(:,p)+delta
        if(present(outflow_change)) then
          do axis=1,3
            if(delta(axis)>0.and.finish(axis)>=upper(axis)) then
              exits(p)=.true.;exit_fraction(p)=min(exit_fraction(p),(upper(axis)-previous(axis,p))/delta(axis))
            else if(delta(axis)<0.and.finish(axis)<=lower(axis)) then
              exits(p)=.true.;exit_fraction(p)=min(exit_fraction(p),(lower(axis)-previous(axis,p))/delta(axis))
            end if
          end do
          if(exits(p)) finish=max(lower,min(upper,previous(:,p)+exit_fraction(p)*delta))
        end if
        drops(p)%position=finish
      end do
      call impact_mpi_level_set(comm,components,properties,settings,phi,molecular_weights,pressures,previous, &
        lower,upper,drops,masses,patches,patch_mass,watermark,crossed,fractions,budget,valid,rejected)
      if(rejected) then
        if(present(model_rejected)) model_rejected=.true.
        return
      end if
      if(.not.valid) return
      allocate(next_remaining(size(drops)),next_pressure(size(drops)),next_depth(size(drops)))
      next_remaining=0;next_pressure(:np)=pressures;next_depth(:np)=depth
      valid=.true.
      do p=1,np
        if(exits(p).and..not.crossed(p).and.drops(p)%mass>0) then
          call liquid_mixture_properties(components,masses(:,p)/drops(p)%mass,drops(p)%temperature,rho,cp,e,valid)
          if(.not.valid) exit
          amount=drops(p)%mass*drops(p)%multiplicity
          exterior(1)=exterior(1)-amount
          exterior(2:4)=exterior(2:4)-amount*drops(p)%velocity
          exterior(5)=exterior(5)-amount*(e+.5_dp*sum(drops(p)%velocity**2))
          exterior(6:)=exterior(6:)-masses(:,p)*drops(p)%multiplicity
          drops(p)%mass=0;masses(:,p)=0
          next_depth(p)=depth(p)+1;local_hits=local_hits+1_int64
        end if
        if(.not.crossed(p)) cycle
        local_hits=local_hits+1_int64
        next_remaining(p)=remaining(p)*(1-exit_fraction(p)*fractions(p));next_depth(p)=depth(p)+1
      end do
      call agree();if(.not.agreed) return
      do p=np+1,size(drops)
        ! Original local IDs are strictly increasing; child ancestry stays local until migration.
        lo=1;hi=np;q=0
        do while(lo<=hi)
          mid=lo+(hi-lo)/2
          if(drops(mid)%id==drops(p)%parent_id) then
            q=mid;exit
          else if(drops(mid)%id<drops(p)%parent_id) then
            lo=mid+1
          else
            hi=mid-1
          end if
        end do
        if(q==0) then
          valid=.false.;exit
        end if
        next_remaining(p)=next_remaining(q);next_pressure(p)=pressures(q);next_depth(p)=next_depth(q)
      end do
      if(valid) valid=all(next_depth<=settings%maximum_hits)
      call agree();if(.not.agreed) return
      call move_alloc(next_remaining,remaining);call move_alloc(next_pressure,pressures);call move_alloc(next_depth,depth)
      ledger=ledger+budget
      active=any(remaining>0.and.drops%mass>0)
      call MPI_Allreduce(active,global_active,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
      if(ierr/=MPI_SUCCESS) return
      if(.not.global_active) exit
    end do
    call MPI_Allreduce(local_hits,total_hits,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.total_hits>int(huge(hits),int64)) return
    call MPI_Allreduce(exterior,global_exterior,5+nc,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    if(present(outflow_change)) outflow_change=global_exterior
    kept=pack(drops,drops%mass>0);allocate(kept_mass(nc,size(kept)))
    do p=1,nc
      kept_mass(p,:)=pack(masses(p,:),drops%mass>0)
    end do
    call move_alloc(kept,parcels);call move_alloc(kept_mass,component_mass)
    call move_alloc(patches,films);call move_alloc(patch_mass,film_mass)
    last_id=watermark;wall_change=ledger;hits=int(total_hits);ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine
  ! Keep pressures attached to their input parcels; migrate only after flight.
  ! Initial and final film migration use private payloads within this transaction.
  subroutine advance_mpi_eb_tree_impacts(comm,tree,owners,components,properties,settings,phi,molecular_weights, &
      pressure,lower,spacing,dt,parcels,component_mass,films,film_mass,last_id,wall_change,hits,ok,model_rejected, &
      outflow_change)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: properties
    type(spray_eb_wall_options), intent(in) :: settings
    real(dp), intent(in) :: phi(0:,0:,0:),molecular_weights(:),pressure(:),lower(3),spacing(3),dt
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:),film_mass(:,:)
    type(spray_wall_film), allocatable, intent(inout) :: films(:)
    integer(int64), intent(inout) :: last_id
    real(dp), intent(inout) :: wall_change(5)
    integer, intent(inout) :: hits
    logical, intent(out) :: ok
    real(dp), optional, intent(inout) :: outflow_change(:)
    logical, optional, intent(out) :: model_rejected
    type(spray_parcel), allocatable :: drops(:)
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: masses(:,:),patch_mass(:,:),exterior(:)
    real(dp) :: upper(3),budget(5)
    integer(int64) :: watermark
    integer :: count,ierr
    logical :: valid,agreed,rejected
    ok=.false.
    if(present(model_rejected)) model_rejected=.false.
    valid=allocated(parcels).and.allocated(component_mass).and.allocated(films).and.allocated(film_mass)
    valid=valid.and.tree%is_valid()
    if(present(outflow_change)) valid=valid.and.size(outflow_change)==5+size(components)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    if(present(outflow_change)) allocate(exterior(size(outflow_change)))
    drops=parcels;masses=component_mass;patches=films;patch_mass=film_mass;watermark=last_id
    call migrate_mpi_amr_wall_films(comm,tree,owners,patches,patch_mass,lower,spacing, &
      settings%maximum_films,valid,phi)
    if(.not.valid) return
    upper=lower+spacing*real(tree%nodes(1)%shape,dp);budget=0;count=0
    call advance_mpi_level_set_impacts(comm,components,properties,settings,phi,molecular_weights,pressure, &
      lower,upper,dt,drops,masses,patches,patch_mass,watermark,budget,count,valid,rejected,exterior)
    if(rejected) then
      if(present(model_rejected)) model_rejected=.true.
      return
    end if
    if(.not.valid) return
    call migrate_eb_tree_spray_parcels(comm,tree,owners,phi,drops,lower,spacing, &
      properties%maximum_parcels,valid,masses)
    if(.not.valid) return
    call migrate_mpi_amr_wall_films(comm,tree,owners,patches,patch_mass,lower,spacing, &
      settings%maximum_films,valid,phi)
    if(.not.valid) return
    call move_alloc(drops,parcels);call move_alloc(masses,component_mass)
    call move_alloc(patches,films);call move_alloc(patch_mass,film_mass)
    if(present(outflow_change)) outflow_change=exterior
    last_id=watermark;wall_change=budget;hits=count;ok=.true.
  end subroutine
end module
