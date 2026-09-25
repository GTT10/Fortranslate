module spray_eb_impact_mod
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use spray_coupling_3d_mod, only: valid_spray_boundary_faces
  use spray_wall_impact_mod, only: count_periodic_film_images
  use spray_parcel_mod
  use spray_mixture_mod, only: valid_component_inventory,liquid_mixture_properties
  use spray_wall_impact_mod, only: spray_wall_film,resolve_mixture_wall_impact,mixture_impact_boiling_temperature
  use spray_wall_impact_mod, only: classify_wall_impact,splash_mass_fractions,splash_velocity_components
  use spray_eb_intersection_mod, only: intersect_spray_level_set
  implicit none
  private
  public :: impact_mixture_level_set,advance_mixture_level_set_impacts
  type, public :: spray_eb_wall_options
    real(dp) :: temperature=300,contact_angle=1.0471975511965977_dp
    integer :: maximum_hits=32,maximum_films=1000000
    logical :: evolve=.false.
    real(dp) :: conductivity=.6_dp,minimum_height=1e-7_dp
  end type
contains
  ! One earliest stationary EB contact; parent%position is the trial endpoint.
  ! Wall budget is received by the wall, as in resolve_mixture_wall_impact.
  ! Caller supplies local wall temperature/thickness and advances remaining flight.
  subroutine impact_mixture_level_set(phi,lower,upper,previous,parent,components,component_mass,properties, &
      molecular_weights,pressure,wall_temperature,film_thickness,contact_angle,first_id, &
      children,child_mass,deposit,film_mass,wall,fraction,regime,hit,ok,model_rejected)
    real(dp), intent(in) :: phi(0:,0:,0:),lower(3),upper(3),previous(3),molecular_weights(:)
    real(dp), intent(in) :: pressure,wall_temperature,film_thickness,contact_angle
    type(spray_parcel), intent(inout) :: parent
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(inout) :: component_mass(:),film_mass(:),wall(5),fraction
    type(spray_options), intent(in) :: properties
    integer(int64), intent(in) :: first_id
    type(spray_parcel), allocatable, intent(inout) :: children(:)
    real(dp), allocatable, intent(inout) :: child_mass(:,:)
    type(spray_wall_film), intent(inout) :: deposit
    integer, intent(inout) :: regime
    logical, intent(out) :: hit,ok
    logical, optional, intent(out) :: model_rejected
    type(spray_parcel) :: trial
    type(spray_parcel), allocatable :: spawned(:)
    type(spray_wall_film) :: deposited
    type(spray_liquid) :: effective
    real(dp), allocatable :: spawned_mass(:,:)
    real(dp) :: y(size(components)),masses(size(components)),film(size(components)),budget(5)
    real(dp) :: location(3),normal(3),at,rho,cp,e,boil,d,mu,un,we,re,kv,kcrit
    real(dp) :: alpha,ms,dm,parts(4),velocities(5),tangent
    integer :: cell(3),kind
    logical :: crossed,valid
    ok=.false.;hit=.false.
    if (present(model_rejected)) model_rejected=.false.
    if(size(components)<1.or.size(components)>16) return
    if(size(component_mass)/=size(components).or.size(film_mass)/=size(components)) return
    if(.not.valid_spray_parcel(parent).or..not.valid_component_inventory(parent,component_mass)) return
    if(parent%mass<=0) return
    if(.not.all(ieee_is_finite([pressure,wall_temperature,film_thickness,contact_angle]))) return
    if(pressure<1.or.pressure>1e10_dp.or.wall_temperature<1.or.wall_temperature>1e5_dp.or.film_thickness<0) return
    at=0;location=0;normal=0;cell=0;crossed=.false.
    call intersect_spray_level_set(phi,lower,upper,previous,parent%position,crossed,at,location,normal,cell,valid)
    if(.not.valid) return
    if(.not.crossed) then
      wall=0;fraction=1;regime=-1;ok=.true.;return
    end if
    ! Geometry detects travel direction; reject a velocity inconsistent with inward impact.
    un=-dot_product(parent%velocity,normal)
    if(un<=0) return
    y=component_mass/parent%mass
    call liquid_mixture_properties(components,y,parent%temperature,rho,cp,e,valid)
    if(.not.valid) return
    call mixture_impact_boiling_temperature(components,y,molecular_weights,pressure,boil,valid)
    if(.not.valid) return
    effective=components(1);effective%density=rho;effective%density_model='constant'
    effective%density_coefficients=0;effective%density_temperature_range=0
    d=parcel_diameter(parent,effective);mu=liquid_viscosity_at_temperature(properties,parent%temperature)
    if(.not.ieee_is_finite(properties%surface_tension)) return
    if(d<=0.or.mu<=0.or.properties%surface_tension<1e-8_dp.or.properties%surface_tension>1e3_dp) return
    ! Same pinned impact parameter conversion as the regular planar dispatcher.
    we=rho*d*un**2/properties%surface_tension;re=un*d*rho/mu
    kv=sqrt(we*sqrt(re))
    kcrit=(26.788_dp*log(pressure*1e-5_dp)+33.546_dp)*exp(-.006_dp*(wall_temperature-260))+36
    if(kcrit<0) return
    ! Identify inadmissible empirical output separately from step-size failures.
    alpha=max(acos(-1._dp)/6,asin(min(1._dp,un/norm2(parent%velocity))))
    call classify_wall_impact(kv,wall_temperature/boil,alpha,kcrit,film_thickness/d,kind,ms,dm,valid)
    if (valid.and.ms>0) then
      call splash_mass_fractions(alpha,ms,film_thickness/d,parts,valid)
      if (valid) then
        tangent=norm2(parent%velocity+un*normal)
        call splash_velocity_components(-un,tangent,kv,film_thickness/d,velocities,valid)
        valid=valid.and.velocities(5)>=0
      end if
    end if
    if (.not.valid) then
      if (present(model_rejected)) model_rejected=.true.
      return
    end if
    trial=parent;masses=component_mass;film=0
    call resolve_mixture_wall_impact(trial,components,masses,location,normal,kv,wall_temperature/boil,kcrit, &
      film_thickness/d,we,re,contact_angle,first_id,spawned,spawned_mass,deposited,film,budget,kind,valid)
    if(.not.valid) return
    parent=trial;component_mass=masses;film_mass=film;deposit=deposited
    call move_alloc(spawned,children);call move_alloc(spawned_mass,child_mass)
    wall=budget;fraction=at;regime=kind;hit=.true.;ok=.true.
  end subroutine
  ! Frozen properties and stationary films; caller owns gas/film exchange.
  ! IDs must be sorted on input. Capacities bound all intermediate births/deposits.
  subroutine advance_mixture_level_set_impacts(phi,lower,upper,dt,maximum_hits,capacity,film_capacity, &
      components,properties,molecular_weights,pressure,wall_temperature,contact_angle, &
      parcels,masses,films,film_masses,last_id,wall,hits,ok,pressures,model_rejected,outflow_change,outer_faces)
    use spray_eb_drift_mod, only: drift_spray_level_set
    use spray_wall_impact_mod, only: valid_wall_component_inventory
    real(dp), intent(in) :: phi(0:,0:,0:),lower(3),upper(3),dt,molecular_weights(:)
    real(dp), intent(in) :: pressure,wall_temperature,contact_angle
    real(dp), optional, intent(in) :: pressures(:)
    integer, intent(in) :: maximum_hits,capacity,film_capacity
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: properties
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: masses(:,:),film_masses(:,:)
    type(spray_wall_film), allocatable, intent(inout) :: films(:)
    integer(int64), intent(inout) :: last_id
    real(dp), optional, intent(inout) :: outflow_change(:)
    real(dp) :: exterior(5+size(components)),exit_fraction,axis_fraction,rho,cp,e,amount
    character(len=*), optional, intent(in) :: outer_faces(6)
    character(len=16) :: faces(6)
    real(dp) :: length(3)
    integer :: axis,exit_face,face,dims(3),image_count
    logical :: exits
    real(dp), intent(inout) :: wall(5)
    integer, intent(inout) :: hits
    logical, intent(out) :: ok
    logical, optional, intent(out) :: model_rejected
    type(spray_parcel), allocatable :: queue(:),children(:),empty(:),result(:)
    type(spray_wall_film), allocatable :: patches(:)
    type(spray_wall_film) :: patch
    real(dp), allocatable :: inventory(:,:),deposits(:,:),child_mass(:,:),remaining(:),local_pressure(:),result_mass(:,:)
    integer, allocatable :: depth(:)
    real(dp) :: previous(3),location(3),normal(3),fraction,thickness,delta(3),budget(5),total(5)
    real(dp) :: film_mass(size(components)),check(5),time_left
    integer :: n,nf,p,f,j,nc,count,cell(3),regime,stat,live
    integer(int64) :: watermark
    logical :: valid,hit,rejected
    ok=.false.
    if (present(model_rejected)) model_rejected=.false.
    if (.not.allocated(parcels).or..not.allocated(masses)) return
    if (.not.allocated(films).or..not.allocated(film_masses)) return
    n=size(parcels);nf=size(films);nc=size(components)
    if(present(outflow_change)) then
      if(size(outflow_change)/=5+nc) return
    end if
    faces='outflow';length=0
    if(present(outer_faces)) then
      if(.not.present(outflow_change)) return
      if(.not.valid_spray_boundary_faces(outer_faces)) return
      faces=outer_faces
      dims=shape(phi)-1
      if(any(dims<1)) return
      if(faces(1)=='periodic') then
        if(any(phi(0,:,:)/=phi(dims(1),:,:))) return
      end if
      if(faces(3)=='periodic') then
        if(any(phi(:,0,:)/=phi(:,dims(2),:))) return
      end if
      if(faces(5)=='periodic') then
        if(any(phi(:,:,0)/=phi(:,:,dims(3)))) return
      end if
      where(faces(1:6:2)=='periodic') length=upper-lower
    end if
    exterior=0
    if (nc<1.or.nc>16.or.capacity<n.or.capacity>1000000.or.capacity<0) return
    if (film_capacity<nf.or.film_capacity>1000000.or.film_capacity<0) return
    if (any(shape(masses)/=[nc,n]).or.any(shape(film_masses)/=[nc,nf])) return
    if (present(pressures)) then
      if (size(pressures)/=n) return
      if (.not.all(ieee_is_finite(pressures))) return
      if (any(pressures<1).or.any(pressures>1e10_dp)) return
    end if
    if (last_id<0.or.last_id>huge(last_id)-4_int64) return
    if (.not.all(ieee_is_finite([pressure,wall_temperature,contact_angle]))) return
    if (pressure<1.or.pressure>1e10_dp.or.wall_temperature<1.or.wall_temperature>1e5_dp) return
    if (contact_angle<1e-6_dp.or.contact_angle>acos(-1._dp)-1e-6_dp) return
    allocate(empty(0));check=0;count=0
    call drift_spray_level_set(phi,lower,upper,dt,maximum_hits,empty,check,count,valid)
    if (.not.valid) return
    do p=1,n
      if (.not.valid_spray_parcel(parcels(p))) return
      if (.not.valid_component_inventory(parcels(p),masses(:,p))) return
      if (parcels(p)%id>last_id) return
      if (p>1) then
        if (parcels(p)%id<=parcels(p-1)%id) return
      end if
    end do
    do f=1,nf
      if (.not.valid_wall_component_inventory(films(f),film_masses(:,f))) return
    end do
    allocate(queue(capacity),inventory(nc,capacity),remaining(capacity),local_pressure(capacity),depth(capacity), &
      patches(film_capacity),deposits(nc,film_capacity),stat=stat)
    if (stat/=0) return
    local_pressure(:n)=pressure
    if (present(pressures)) local_pressure(:n)=pressures
    queue(:n)=parcels;inventory(:,:n)=masses;remaining(:n)=dt;depth(:n)=0
    patches(:nf)=films;deposits(:,:nf)=film_masses;watermark=last_id;total=0;count=0;p=1
    do while(p<=n)
      do while(queue(p)%mass>0)
        previous=queue(p)%position
        queue(p)%position=previous+remaining(p)*queue(p)%velocity
        exits=.false.;exit_fraction=1;exit_face=0
        if(present(outflow_change)) then
          if(any(previous<lower).or.any(previous>upper)) return
          do axis=1,3
            if(remaining(p)*queue(p)%velocity(axis)==0) cycle
            if(queue(p)%velocity(axis)>0.and.queue(p)%position(axis)>=upper(axis)) then
              axis_fraction=(upper(axis)-previous(axis))/(remaining(p)*queue(p)%velocity(axis));face=2*axis
            else if(queue(p)%velocity(axis)<0.and.queue(p)%position(axis)<=lower(axis)) then
              axis_fraction=(lower(axis)-previous(axis))/(remaining(p)*queue(p)%velocity(axis));face=2*axis-1
            else
              cycle
            end if
            if(axis_fraction>exit_fraction) cycle
            if(exits.and.axis_fraction==exit_fraction) cycle
            exits=.true.;exit_fraction=axis_fraction;exit_face=face
          end do
          if(exits) queue(p)%position=max(lower,min(upper,previous+exit_fraction*remaining(p)*queue(p)%velocity))
        end if
        hit=.false.;fraction=0;location=0;normal=0;cell=0
        call intersect_spray_level_set(phi,lower,upper,previous,queue(p)%position, &
          hit,fraction,location,normal,cell,valid)
        if (.not.valid) return
        if (.not.hit) then
          if(exits) then
            if(depth(p)>=maximum_hits) return
            if(faces(exit_face)=='periodic'.or.faces(exit_face)=='slip_wall'.or. &
              faces(exit_face)=='no_slip_wall'.or.faces(exit_face)=='isothermal_wall') then
              axis=(exit_face+1)/2
              if(faces(exit_face)=='periodic') then
                queue(p)%position(axis)=merge(lower(axis),upper(axis),modulo(exit_face,2)==0)
              else
                amount=queue(p)%mass*queue(p)%multiplicity
                exterior(axis+1)=exterior(axis+1)-2*amount*queue(p)%velocity(axis)
                queue(p)%velocity(axis)=-queue(p)%velocity(axis)
              end if
              depth(p)=depth(p)+1;count=count+1;remaining(p)=remaining(p)*(1-exit_fraction)
              if(remaining(p)==0) exit
              cycle
            end if
            call liquid_mixture_properties(components,inventory(:,p)/queue(p)%mass,queue(p)%temperature,rho,cp,e,valid)
            if(.not.valid) return
            amount=queue(p)%mass*queue(p)%multiplicity
            exterior(1)=exterior(1)-amount
            exterior(2:4)=exterior(2:4)-amount*queue(p)%velocity
            exterior(5)=exterior(5)-amount*(e+.5_dp*sum(queue(p)%velocity**2))
            exterior(6:)=exterior(6:)-inventory(:,p)*queue(p)%multiplicity
            queue(p)%mass=0;inventory(:,p)=0;count=count+1
          end if
          exit
        end if
        if (depth(p)>=maximum_hits) return
        thickness=0
        do f=1,nf
          if (maxval(abs(patches(f)%normal-normal))>1e-12_dp) cycle
          delta=location-patches(f)%position
          call count_periodic_film_images(delta,normal,.5_dp*patches(f)%diameter, &
            1e-10_dp*maxval(upper-lower),length,image_count,valid)
          if(.not.valid) return
          if(patches(f)%height*real(image_count,dp)>huge(1._dp)-thickness) return
          thickness=thickness+patches(f)%height*real(image_count,dp)
        end do
        if (watermark>huge(watermark)-4_int64) return
        call impact_mixture_level_set(phi,lower,upper,previous,queue(p),components,inventory(:,p),properties, &
          molecular_weights,local_pressure(p),wall_temperature,thickness,contact_angle,watermark+1_int64, &
          children,child_mass,patch,film_mass,budget,fraction,regime,hit,valid,rejected)
        if (rejected) then
          if (present(model_rejected)) model_rejected=.true.
          return
        end if
        if (.not.valid.or..not.hit) return
        if (size(children)>capacity-n) return
        if (patch%mass>0) then
          if (nf==film_capacity) return
          nf=nf+1;patches(nf)=patch;deposits(:,nf)=film_mass
        end if
        time_left=remaining(p)*(1-exit_fraction*fraction);depth(p)=depth(p)+1
        do j=1,size(children)
          n=n+1;queue(n)=children(j);inventory(:,n)=child_mass(:,j)
          remaining(n)=time_left;depth(n)=depth(p);local_pressure(n)=local_pressure(p)
        end do
        watermark=watermark+int(size(children),int64)
        remaining(p)=time_left;total=total+budget;count=count+1
        if (time_left==0) exit
      end do
      p=p+1
    end do
    ! Remove consumed parents while preserving deterministic increasing IDs.
    live=0
    do p=1,n
      if (queue(p)%mass>0) live=live+1
    end do
    allocate(result(live),result_mass(nc,live),stat=stat)
    if (stat/=0) return
    j=0
    do p=1,n
      if (queue(p)%mass<=0) cycle
      j=j+1;result(j)=queue(p);result_mass(:,j)=inventory(:,p)
    end do
    call move_alloc(result,parcels);call move_alloc(result_mass,masses)
    films=patches(:nf);film_masses=deposits(:,:nf)
    if(present(outflow_change)) outflow_change=exterior
    last_id=watermark;wall=total;hits=count;ok=.true.
  end subroutine
end module
