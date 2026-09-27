! Breakup equations adapted from PelePhysics TABBreakup.H, ReitzKHRT.H and SpraySB.cpp,
! revision 6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b.
! Copyright and redistribution terms: references/pelephysics-spray-license.txt.
module spray_breakup_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use spray_parcel_mod
  use spray_mixture_mod, only: valid_component_inventory,liquid_mixture_properties
  implicit none
  private
  public :: advance_etab_breakup, advance_khrt_breakup, advance_spray_breakup, split_spray_parcels
  type, public :: spray_khrt_shedding_state
    real(dp) :: initial_count=0,shed_volume=0
  end type
  type, public :: spray_khrt_release
    real(dp) :: mass=0,drop_mass=0,multiplicity=0,diameter=0,tangential_speed=0
  end type
  public :: advance_khrt_shedding,construct_khrt_children,advance_mixture_khrt_shedding
  public :: khrt_reference_child_count,construct_reference_khrt_children,khrt_reference_uniforms
contains
  ! Stateless xorshift64 samples keyed by persistent parent ID and pre-release
  ! global ID watermark. No rank-local RNG state is advanced on rejection.
  subroutine khrt_reference_uniforms(parent_id,watermark,uniforms)
    integer(int64), intent(in) :: parent_id,watermark
    real(dp), intent(out) :: uniforms(:)
    integer(int64) :: state
    integer :: n
    state=ieor(ieor(parent_id,ishftc(watermark,29)),88172645463325252_int64)
    if(state==0) state=88172645463325252_int64
    do n=-7,size(uniforms)
      state=ieor(state,shiftl(state,13))
      state=ieor(state,shiftr(state,7))
      state=ieor(state,shiftl(state,17))
      if(n>=1) uniforms(n)=real(iand(state,9007199254740991_int64),dp)/9007199254740992._dp
    enddo
  end subroutine

  subroutine khrt_reference_child_count(multiplicity,exponent,capacity,number,ok)
    real(dp), intent(in) :: multiplicity,exponent
    integer, intent(in) :: capacity
    integer, intent(inout) :: number
    logical, intent(out) :: ok
    real(dp) :: weight,count
    ok=.false.
    if(.not.all(ieee_is_finite([multiplicity,exponent]))) return
    if(multiplicity<=0.or.multiplicity>1e18_dp.or.exponent<0.or.exponent>1.or.capacity<1) return
    ! SpraySB.cpp: N_s = N**factor; N_d = max(1, int(N/N_s)).
    weight=multiplicity**exponent;count=multiplicity/weight
    ! Check before real-to-integer conversion, including truncation at capacity.
    if(count>=real(capacity,dp)+1._dp) return
    number=max(1,int(count));ok=.true.
  end subroutine

  subroutine construct_reference_khrt_children(parent,release,axis,exponent,uniforms,first_id, &
      children,gas_exchange,ok)
    type(spray_parcel), intent(in) :: parent
    type(spray_khrt_release), intent(in) :: release
    real(dp), intent(in) :: axis(3),exponent,uniforms(:)
    integer(int64), intent(in) :: first_id
    type(spray_parcel), allocatable, intent(inout) :: children(:)
    ! Momentum and total-energy compensation; mass/species do not change.
    real(dp), intent(inout) :: gas_exchange(4)
    logical, intent(out) :: ok
    type(spray_parcel), allocatable :: candidate(:)
    real(dp) :: normal(3),testvec(3),psi_axis(3),beta_axis(3),dv(3),angle,exchange(4),weight
    integer :: number,n,stat
    ok=.false.
    if(.not.valid_spray_parcel(parent)) return
    if(.not.all(ieee_is_finite([axis,release%mass,release%drop_mass,release%multiplicity, &
        release%diameter,release%tangential_speed]))) return
    if(maxval(abs(axis))>1e9_dp) return
    if(norm2(axis)<1e-12_dp) return
    if(release%mass<=0.or.release%drop_mass<=0.or.release%drop_mass>1.or.release%diameter<=0) return
    if(release%multiplicity<=0.or.release%multiplicity>1e18_dp) return
    if(release%tangential_speed<0.or.release%tangential_speed>1e7_dp) return
    if(abs(release%mass-release%drop_mass*release%multiplicity)>1e-12_dp*release%mass) return
    number=0
    call khrt_reference_child_count(release%multiplicity,exponent,size(uniforms),number,ok)
    if(.not.ok) return
    ok=.false.
    if(size(uniforms)/=number) return
    if(.not.all(ieee_is_finite(uniforms))) return
    if(any(uniforms<0).or.any(uniforms>=1)) return
    if(first_id<=parent%id.or.first_id>huge(first_id)-int(number-1,int64)) return
    normal=axis/norm2(axis);testvec=[1._dp,0._dp,0._dp]
    psi_axis=cross(testvec,normal)
    if(norm2(psi_axis)<1e-5_dp) then
      testvec=[0._dp,1._dp,0._dp];psi_axis=cross(testvec,normal)
    endif
    psi_axis=psi_axis/norm2(psi_axis)
    beta_axis=cross(psi_axis,normal);beta_axis=beta_axis/norm2(beta_axis)
    allocate(candidate(number),stat=stat)
    if(stat/=0) return
    weight=release%multiplicity/real(number,dp);exchange=0
    do n=1,number
      angle=2*acos(-1._dp)*uniforms(n)
      dv=release%tangential_speed*(sin(angle)*psi_axis+cos(angle)*beta_axis)
      candidate(n)=parent;candidate(n)%id=first_id+int(n-1,int64);candidate(n)%parent_id=parent%id
      candidate(n)%mass=release%drop_mass;candidate(n)%multiplicity=weight
      candidate(n)%velocity=parent%velocity+dv
      candidate(n)%deformation=0;candidate(n)%deformation_rate=0;candidate(n)%rt_time=0
      candidate(n)%khrt_initial_count=weight;candidate(n)%khrt_shed_volume=0
      if(.not.valid_spray_parcel(candidate(n))) return
      exchange(:3)=exchange(:3)-release%drop_mass*weight*dv
      exchange(4)=exchange(4)-release%drop_mass*weight*(dot_product(parent%velocity,dv)+.5_dp*sum(dv**2))
    enddo
    if(.not.all(ieee_is_finite(exchange))) return
    call move_alloc(candidate,children);gas_exchange=exchange;ok=.true.
  contains
    pure function cross(a,b) result(c)
      real(dp), intent(in) :: a(3),b(3)
      real(dp) :: c(3)
      c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
    end function
  end subroutine

  subroutine advance_mixture_khrt_shedding(components,options,parcels,component_mass,gas_density, &
      gas_velocity,reynolds,dt,last_id,gas_energy,ok,eligible,reference_factor,gas_momentum)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:),gas_energy(:)
    real(dp), intent(in) :: gas_density(:),gas_velocity(:,:),reynolds(:),dt
    integer(int64), intent(inout) :: last_id
    logical, intent(out) :: ok
    logical, optional, intent(in) :: eligible(:)
    real(dp), optional, intent(in) :: reference_factor
    real(dp), allocatable, optional, intent(inout) :: gas_momentum(:,:)
    real(dp), allocatable :: momentum(:,:),uniforms(:)
    real(dp) :: exchange(4)
    integer, allocatable :: child_counts(:)
    type(spray_parcel), allocatable :: candidate(:),children(:),result(:)
    type(spray_khrt_release), allocatable :: releases(:)
    type(spray_khrt_shedding_state) :: pending
    type(spray_liquid) :: effective
    real(dp), allocatable :: masses(:,:),fractions(:,:),energy(:),result_mass(:,:)
    real(dp) :: rho,cp,e,axis(3)
    integer :: np,nc,p,n,number,cursor
    integer(int64) :: next_id,first_id
    logical :: valid
    ok=.false.
    valid=allocated(parcels).and.allocated(component_mass).and.valid_spray_options(options)
    valid=valid.and.options%breakup_model=='khrt'.and.options%khrt_children
    valid=valid.and.(present(reference_factor).eqv.present(gas_momentum))
    if(present(reference_factor)) then
      if(ieee_is_finite(reference_factor)) then
        valid=valid.and.reference_factor>=0.and.reference_factor<=1
      else
        valid=.false.
      endif
    endif
    if (.not.valid) return
    np=size(parcels);nc=size(components)
    valid=np<=1000000.and.nc>=1.and.nc<=16.and.all(shape(component_mass)==[nc,np])
    valid=valid.and.size(gas_density)==np.and.size(reynolds)==np.and.all(shape(gas_velocity)==[3,np])
    valid=valid.and.ieee_is_finite(dt).and.dt>=0.and.dt<=1e6_dp
    if(present(eligible)) valid=valid.and.size(eligible)==np
    if (.not.valid) return
    if (last_id<0) return
    candidate=parcels;masses=component_mass
    allocate(releases(np),fractions(nc,np),energy(np),momentum(3,np),child_counts(np))
    energy=0;fractions=0;number=0;momentum=0;child_counts=0
    do p=1,np
      valid=valid_component_inventory(candidate(p),masses(:,p)).and.candidate(p)%id<=last_id
      if (p>1) valid=valid.and.candidate(p)%id>candidate(p-1)%id
      if (.not.valid) exit
      if (candidate(p)%mass<=0) cycle
      if(present(eligible)) then
        if(.not.eligible(p)) cycle
      end if
      fractions(:,p)=masses(:,p)/candidate(p)%mass
      call liquid_mixture_properties(components,fractions(:,p),candidate(p)%temperature,rho,cp,e,valid)
      if (.not.valid) exit
      effective=components(1);effective%density_model='constant';effective%density=rho
      effective%density_coefficients=0;effective%density_temperature_range=0
      effective%heat_capacity_model='constant';effective%heat_capacity=cp
      effective%heat_capacity_coefficients=0;effective%heat_capacity_temperature_range=0
      effective%reference_temperature=candidate(p)%temperature;effective%reference_energy=e
      effective%boiling_model='reject';effective%critical_temperature=0
      effective%normal_boiling_temperature=0
      pending=spray_khrt_shedding_state(candidate(p)%khrt_initial_count,candidate(p)%khrt_shed_volume)
      call advance_khrt_shedding(candidate(p),effective,options,gas_density(p),gas_velocity(:,p),reynolds(p), &
        dt,options%khrt_reference_diameter_cubed,pending,releases(p),valid)
      if (.not.valid) exit
      candidate(p)%khrt_initial_count=pending%initial_count;candidate(p)%khrt_shed_volume=pending%shed_volume
      masses(:,p)=candidate(p)%mass*fractions(:,p)
      if (releases(p)%mass>0) then
        child_counts(p)=4
        if(present(reference_factor)) call khrt_reference_child_count(releases(p)%multiplicity, &
          reference_factor,options%maximum_parcels-np-number,child_counts(p),valid)
        if(.not.valid) exit
        number=number+child_counts(p)
      endif
    end do
    valid=valid.and.number<=options%maximum_parcels-np
    if (.not.valid) return
    if (last_id>huge(last_id)-int(number,int64)) return
    next_id=last_id+int(number,int64)
    first_id=last_id
    if(number>0) first_id=last_id+1_int64
    allocate(result(np+number),result_mass(nc,np+number))
    result(:np)=candidate;result_mass(:,:np)=masses;cursor=np
    do p=1,np
      if (releases(p)%mass<=0) cycle
      axis=candidate(p)%velocity
      if (norm2(axis)<1e-12_dp) axis=gas_velocity(:,p)-candidate(p)%velocity
      if(present(reference_factor)) then
        if(allocated(uniforms)) deallocate(uniforms)
        allocate(uniforms(child_counts(p)))
        call khrt_reference_uniforms(candidate(p)%id,last_id,uniforms)
        call construct_reference_khrt_children(candidate(p),releases(p),axis,reference_factor,uniforms, &
          first_id+int(cursor-np,int64),children,exchange,valid)
        if(valid) then
          momentum(:,p)=exchange(:3);energy(p)=exchange(4)
        endif
      else
        call construct_khrt_children(candidate(p),releases(p),axis, &
          first_id+int(cursor-np,int64),children,energy(p),valid)
      endif
      if (.not.valid) exit
      result(cursor+1:cursor+child_counts(p))=children
      do n=1,child_counts(p)
        result_mass(:,cursor+n)=children(n)%mass*fractions(:,p)
      end do
      cursor=cursor+child_counts(p)
    end do
    if (.not.valid) return
    call move_alloc(result,parcels);call move_alloc(result_mass,component_mass);call move_alloc(energy,gas_energy)
    if(present(gas_momentum)) call move_alloc(momentum,gas_momentum)
    last_id=next_id;ok=.true.
  end subroutine
  subroutine construct_khrt_children(parent,release,axis,first_id,children,gas_energy,ok)
    type(spray_parcel), intent(in) :: parent
    type(spray_khrt_release), intent(in) :: release
    real(dp), intent(in) :: axis(3)
    integer(int64), intent(in) :: first_id
    type(spray_parcel), allocatable, intent(inout) :: children(:)
    real(dp), intent(out) :: gas_energy
    logical, intent(out) :: ok
    type(spray_parcel) :: candidate(4)
    real(dp) :: normal(3),helper(3),tangent(3),other(3),directions(3,4),energy
    integer :: n
    ok=.false.;gas_energy=0
    if (.not.valid_spray_parcel(parent)) return
    if (.not.all(ieee_is_finite([axis,release%mass,release%drop_mass,release%multiplicity, &
        release%diameter,release%tangential_speed]))) return
    if (norm2(axis)<1e-12_dp.or.maxval(abs(axis))>1e9_dp) return
    if (first_id<=parent%id.or.first_id>huge(first_id)-3_int64) return
    if (release%mass<=0.or.release%drop_mass<=0.or.release%drop_mass>1.or.release%multiplicity<=0) return
    if (release%multiplicity>1e18_dp.or.release%tangential_speed<0.or.release%tangential_speed>1e7_dp) return
    if (abs(release%mass-release%drop_mass*release%multiplicity)>1e-12_dp*release%mass) return
    normal=axis/norm2(axis);helper=[1.0_dp,0.0_dp,0.0_dp]
    if (abs(normal(1))>.9_dp) helper=[0.0_dp,1.0_dp,0.0_dp]
    tangent=helper-dot_product(helper,normal)*normal;tangent=tangent/norm2(tangent)
    other=[normal(2)*tangent(3)-normal(3)*tangent(2),normal(3)*tangent(1)-normal(1)*tangent(3), &
      normal(1)*tangent(2)-normal(2)*tangent(1)]
    directions(:,1)=tangent;directions(:,2)=-tangent;directions(:,3)=other;directions(:,4)=-other
    do n=1,4
      candidate(n)=parent;candidate(n)%id=first_id+int(n-1,int64);candidate(n)%parent_id=parent%id
      candidate(n)%mass=release%drop_mass;candidate(n)%multiplicity=release%multiplicity/4
      candidate(n)%velocity=parent%velocity+release%tangential_speed*directions(:,n)
      candidate(n)%deformation=0;candidate(n)%deformation_rate=0;candidate(n)%rt_time=0
      candidate(n)%khrt_initial_count=0;candidate(n)%khrt_shed_volume=0
      if (.not.valid_spray_parcel(candidate(n))) return
    end do
    ! Symmetric azimuth quadrature preserves emitted momentum; gas supplies transverse kinetic energy.
    energy=-.5_dp*release%mass*release%tangential_speed**2
    if (.not.ieee_is_finite(energy)) return
    children=candidate;gas_energy=energy;ok=.true.
  end subroutine

  subroutine advance_khrt_shedding(parcel,liquid,options,gas_density,gas_velocity,reynolds_number,dt, &
      reference_diameter_cubed,state,release,ok)
    type(spray_parcel), intent(inout) :: parcel
    type(spray_liquid), intent(in) :: liquid
    type(spray_options), intent(in) :: options
    real(dp), intent(in) :: gas_density,gas_velocity(3),reynolds_number,dt,reference_diameter_cubed
    type(spray_khrt_shedding_state), intent(inout) :: state
    type(spray_khrt_release), intent(out) :: release
    logical, intent(out) :: ok
    type(spray_parcel) :: candidate
    type(spray_khrt_shedding_state) :: pending
    type(spray_khrt_release) :: emitted
    real(dp) :: old_radius,new_radius,rs,wave,frequency,total,rho,ns,remaining
    release=spray_khrt_release();emitted=release;ok=.false.
    if (.not.all(ieee_is_finite([reference_diameter_cubed,state%initial_count,state%shed_volume]))) return
    if (reference_diameter_cubed<=0.or.reference_diameter_cubed>1e24_dp) return
    if (state%initial_count<0.or.state%initial_count>1e18_dp.or.state%shed_volume<0.or.state%shed_volume>1e24_dp) return
    if (state%initial_count==0.and.state%shed_volume/=0) return
    candidate=parcel;pending=state
    call advance_khrt_breakup(candidate,liquid,options,gas_density,gas_velocity,reynolds_number,dt,ok,wave,frequency)
    if (.not.ok) return
    ok=.false.
    if (wave>0) then
      rho=liquid_density_at_temperature(liquid,parcel%temperature)
      total=parcel%mass*parcel%multiplicity
      old_radius=.5_dp*parcel_diameter(parcel,liquid)
      new_radius=.5_dp*parcel_diameter(candidate,liquid)
      if (pending%initial_count==0) pending%initial_count=parcel%multiplicity
      pending%shed_volume=pending%shed_volume+pending%initial_count*max(0.0_dp,old_radius**3-new_radius**3)
      rs=options%khrt_b0*wave
      ns=pending%shed_volume/rs**3
      if (pending%shed_volume/(parcel%multiplicity*old_radius**3)>.9_dp) then
        rs=max(rs,.5_dp*options%minimum_diameter)
        emitted%mass=total;candidate%mass=0
        pending%shed_volume=0;pending%initial_count=0
      else if (8*pending%shed_volume/max(reference_diameter_cubed,epsilon(1.0_dp)*1e-6_dp)>.03_dp.and. &
          rs>.5_dp*options%minimum_diameter.and.ns>parcel%multiplicity) then
        ! pending volume is N*r^3, not kg or geometric sphere volume.
        emitted%mass=4*acos(-1.0_dp)*rho*pending%shed_volume/3
        remaining=total-emitted%mass
        if (remaining<=0) return
        candidate%multiplicity=remaining/candidate%mass
        pending%initial_count=candidate%multiplicity;pending%shed_volume=0
      end if
      if (emitted%mass>0) then
        emitted%diameter=2*rs;emitted%drop_mass=4*acos(-1.0_dp)*rho*rs**3/3
        emitted%multiplicity=emitted%mass/emitted%drop_mass
        emitted%tangential_speed=.188_dp*frequency*wave
        if (.not.all(ieee_is_finite([emitted%mass,emitted%drop_mass,emitted%multiplicity,emitted%tangential_speed]))) return
        if (emitted%drop_mass<=0.or.emitted%drop_mass>1.or.emitted%multiplicity>1e18_dp) return
        if (emitted%tangential_speed>1e7_dp) return
      end if
    end if
    if (.not.valid_spray_parcel(candidate)) return
    if (.not.all(ieee_is_finite([pending%initial_count,pending%shed_volume]))) return
    if (pending%initial_count>1e18_dp.or.pending%shed_volume>1e24_dp) return
    parcel=candidate;state=pending;release=emitted;ok=.true.
  end subroutine

  subroutine advance_spray_breakup(parcel,liquid,options,gas_density,gas_velocity,film_density,gas_viscosity,dt,ok)
    type(spray_parcel), intent(inout) :: parcel
    type(spray_liquid), intent(in) :: liquid
    type(spray_options), intent(in) :: options
    real(dp), intent(in) :: gas_density,gas_velocity(3),film_density,gas_viscosity,dt
    logical, intent(out) :: ok
    real(dp) :: reynolds
    if (trim(options%breakup_model)=='khrt') then
      ok=.false.
      if (.not.all(ieee_is_finite([film_density,gas_viscosity]))) return
      if (film_density<=0.or.film_density>1e8_dp.or.gas_viscosity<1e-12_dp.or.gas_viscosity>1e3_dp) return
      if (.not.valid_spray_parcel(parcel).or..not.valid_spray_liquid(liquid)) return
      if (.not.all(ieee_is_finite(gas_velocity))) return
      if (maxval(abs(gas_velocity))>1e7_dp) return
      reynolds=film_density*sqrt(sum((gas_velocity-parcel%velocity)**2))*parcel_diameter(parcel,liquid)/gas_viscosity
      call advance_khrt_breakup(parcel,liquid,options,gas_density,gas_velocity,reynolds,dt,ok)
    else
      call advance_etab_breakup(parcel,liquid,options,gas_density,gas_velocity,dt,ok)
    end if
  end subroutine

  ! ReitzKHRT.H representative-drop branch (make_new_drops=false).
  subroutine advance_khrt_breakup(parcel,liquid,options,gas_density,gas_velocity,reynolds_number,dt,ok,kh_wave,kh_frequency)
    type(spray_parcel), intent(inout) :: parcel
    type(spray_liquid), intent(in) :: liquid
    type(spray_options), intent(in) :: options
    real(dp), intent(in) :: gas_density,gas_velocity(3),reynolds_number,dt
    real(dp), intent(out), optional :: kh_wave,kh_frequency
    real(dp) :: release_wave,release_frequency
    logical, intent(out) :: ok
    type(spray_parcel) :: candidate
    real(dp) :: rho,liquid_mu
    real(dp) :: radius,child,speed,re,cd,we,z,t,lambda,omega,tau,acceleration,momdiff,rate,total
    logical :: rt_breakup
    real(dp), parameter :: pi=acos(-1.0_dp)
    if (present(kh_wave)) kh_wave=0
    if (present(kh_frequency)) kh_frequency=0
    release_wave=0;release_frequency=0
    ok=.false.
    if (.not.valid_spray_parcel(parcel).or..not.valid_spray_liquid(liquid).or..not.valid_spray_options(options)) return
    if (trim(options%breakup_model)/='khrt') return
    if (.not.all(ieee_is_finite([gas_density,gas_velocity,reynolds_number,dt]))) return
    if (gas_density<=0.or.gas_density>1e8_dp.or.reynolds_number<0.or.reynolds_number>1e30_dp) return
    if (maxval(abs(gas_velocity))>1e7_dp.or.dt<0.or.dt>1e6_dp) return
    if (parcel%mass<=0.or.dt==0) then
      ok=.true.;return
    end if
    rho=liquid_density_at_temperature(liquid,parcel%temperature)
    if (rho<=0) return
    liquid_mu=liquid_viscosity_at_temperature(options,parcel%temperature)
    if (liquid_mu<=0) return
    candidate=parcel;radius=.5_dp*parcel_diameter(parcel,liquid);child=radius
    if (radius<=.5_dp*options%minimum_diameter) then
      ok=.true.;return
    end if
    speed=sqrt(sum((gas_velocity-parcel%velocity)**2));total=parcel%mass*parcel%multiplicity
    re=reynolds_number
    cd=0
    if (re>1000) then
      cd=.424_dp
    else if (re>1) then
      cd=24/re*(1+re**(2.0_dp/3.0_dp)/6)
    else if (re>0) then
      cd=24/re
    end if
    we=gas_density*speed**2*radius/options%surface_tension
    z=liquid_mu/sqrt(rho*radius*options%surface_tension)
    t=z*sqrt(we);rt_breakup=.false.
    if (candidate%rt_time>=0) then
      acceleration=.375_dp*cd*gas_density*speed**2/(rho*radius)
      momdiff=acceleration*(rho-gas_density)
      if (momdiff>0) then
        omega=sqrt(2/(3*sqrt(3*options%surface_tension))*momdiff*sqrt(momdiff)/(rho+gas_density))
        lambda=options%khrt_c3*pi/sqrt(momdiff/(3*options%surface_tension))
        if (radius<lambda) then
          candidate%rt_time=0
        else
          candidate%rt_time=candidate%rt_time+dt
        end if
        if (candidate%rt_time>1/omega.and.lambda<radius) then
          child=lambda;candidate%rt_time=0;rt_breakup=.true.
        end if
      else
        candidate%rt_time=0
      end if
    end if
    if (.not.rt_breakup) then
      lambda=radius*9.02_dp*(1+.45_dp*sqrt(z))*(1+.4_dp*t**.7_dp)/(1+.865_dp*we**1.67_dp)**.6_dp
      omega=(.34_dp+.385_dp*we*sqrt(we))/((1+z)*(1+1.4_dp*t**.6_dp))* &
        sqrt(options%surface_tension/(rho*radius**3))
      tau=3.726_dp*options%khrt_b1*radius/(lambda*omega)
      if (options%khrt_b0*lambda<radius.and.we>6) then
        release_wave=lambda;release_frequency=omega
        rate=dt/tau
        child=(rate*options%khrt_b0*lambda+radius)/(1+rate)
        candidate%rt_time=0
      end if
    end if
    child=max(.5_dp*options%minimum_diameter,child)
    if (child<radius) then
      candidate%mass=parcel%mass*(child/radius)**3
      candidate%multiplicity=total/candidate%mass
    end if
    if (.not.valid_spray_parcel(candidate)) return
    if (present(kh_wave)) kh_wave=release_wave
    if (present(kh_frequency)) kh_frequency=release_frequency
    parcel=candidate;ok=.true.
  end subroutine

  subroutine split_spray_parcels(parcels,maximum_multiplicity,maximum_parcels,ok,watermark,component_mass)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), intent(in) :: maximum_multiplicity
    integer, intent(in) :: maximum_parcels
    logical, intent(out) :: ok
    integer(int64), optional, intent(inout) :: watermark
    real(dp), allocatable, optional, intent(inout) :: component_mass(:,:)
    real(dp), allocatable :: candidate_mass(:,:)
    type(spray_parcel), allocatable :: candidate(:)
    integer, allocatable :: pieces(:)
    integer :: p,n,total,extra,cursor
    integer(int64) :: last_id
    real(dp) :: weight,count_weight,shed_weight
    ok=.false.
    if (.not.allocated(parcels)) return
    if (.not.ieee_is_finite(maximum_multiplicity).or.maximum_multiplicity<0.or.maximum_multiplicity>1e18_dp) return
    if (maximum_multiplicity>0.and.maximum_multiplicity<1e-12_dp) return
    if (maximum_parcels<1.or.maximum_parcels>1000000.or.size(parcels)>maximum_parcels) return
    if (present(watermark)) then
      if (watermark<0) return
      do p=1,size(parcels)
        if (parcels(p)%id>watermark) return
      end do
    end if
    if (present(component_mass)) then
      if (.not.allocated(component_mass)) return
      if (size(component_mass,2)/=size(parcels)) return
      if (size(component_mass,1)<1.or.size(component_mass,1)>16) return
      do p=1,size(parcels)
        if (.not.valid_component_inventory(parcels(p),component_mass(:,p))) return
      end do
    end if
    if (maximum_multiplicity==0) then
      ok=.true.;return
    end if
    allocate(pieces(size(parcels)));pieces=1;total=size(parcels);last_id=0
    do p=1,size(parcels)
      if (.not.valid_spray_parcel(parcels(p))) return
      if (parcels(p)%id<=last_id) return
      last_id=parcels(p)%id
      if (parcels(p)%mass<=0.or.parcels(p)%multiplicity<=maximum_multiplicity) cycle
      if (parcels(p)%multiplicity/maximum_multiplicity>real(maximum_parcels,dp)) return
      pieces(p)=ceiling(parcels(p)%multiplicity/maximum_multiplicity)
      extra=pieces(p)-1
      if (total>maximum_parcels-extra) return
      total=total+extra
    end do
    if (total==size(parcels)) then
      ok=.true.;return
    end if
    if (present(watermark)) last_id=watermark
    if (last_id>huge(last_id)-int(total-size(parcels),int64)) return
    allocate(candidate(total)); candidate(:size(parcels))=parcels;cursor=size(parcels)
    if (present(component_mass)) then
      allocate(candidate_mass(size(component_mass,1),total))
      candidate_mass(:,:size(parcels))=component_mass
    end if
    do p=1,size(parcels)
      if (pieces(p)==1) cycle
      weight=parcels(p)%multiplicity/real(pieces(p),dp)
      candidate(p)%multiplicity=weight
      ! Both KH accumulators are extensive parcel inventories; duplicating
      ! them would multiply the next release mass by the number of pieces.
      count_weight=parcels(p)%khrt_initial_count/real(pieces(p),dp)
      shed_weight=parcels(p)%khrt_shed_volume/real(pieces(p),dp)
      candidate(p)%khrt_initial_count=count_weight;candidate(p)%khrt_shed_volume=shed_weight
      if (.not.valid_spray_parcel(candidate(p))) return
      do n=2,pieces(p)
        cursor=cursor+1;candidate(cursor)=parcels(p)
        if (present(component_mass)) candidate_mass(:,cursor)=component_mass(:,p)
        last_id=last_id+1;candidate(cursor)%id=last_id
        candidate(cursor)%parent_id=parcels(p)%id
        candidate(cursor)%multiplicity=weight
        candidate(cursor)%khrt_initial_count=count_weight;candidate(cursor)%khrt_shed_volume=shed_weight
        if (n==pieces(p)) then
          candidate(cursor)%multiplicity=parcels(p)%multiplicity-real(n-1,dp)*weight
          candidate(cursor)%khrt_initial_count=parcels(p)%khrt_initial_count-real(n-1,dp)*count_weight
          candidate(cursor)%khrt_shed_volume=parcels(p)%khrt_shed_volume-real(n-1,dp)*shed_weight
        end if
        if (.not.valid_spray_parcel(candidate(cursor))) return
      end do
    end do
    call move_alloc(candidate,parcels)
    if (present(component_mass)) call move_alloc(candidate_mass,component_mass)
    if (present(watermark)) watermark=last_id
    ok=.true.
  end subroutine

  subroutine advance_etab_breakup(parcel,liquid,options,gas_density,gas_velocity,dt,ok)
    type(spray_parcel), intent(inout) :: parcel
    type(spray_liquid), intent(in) :: liquid
    type(spray_options), intent(in) :: options
    real(dp), intent(in) :: gas_density,gas_velocity(3),dt
    logical, intent(out) :: ok
    type(spray_parcel) :: candidate
    real(dp) :: rho,liquid_mu
    real(dp) :: radius,child,td,omega2,omega,we,wer,amplitude,y,ydot,yn,ydn
    real(dp) :: elapsed,step,speed,phi,angle,event_time,kbr,etb,ratio,total_mass
    real(dp), parameter :: pi=acos(-1.0_dp)
    integer :: iteration
    ok=.false.
    if (.not.valid_spray_parcel(parcel).or..not.valid_spray_liquid(liquid).or..not.valid_spray_options(options)) return
    if (.not.all(ieee_is_finite([gas_density,gas_velocity,dt]))) return
    if (gas_density<=0.or.gas_density>1e8_dp.or.maxval(abs(gas_velocity))>1e7_dp.or.dt<0.or.dt>1e6_dp) return
    if (trim(options%breakup_model)=='khrt') return
    if (trim(options%breakup_model)=='none'.or.parcel%mass<=0.or.dt==0) then
      ok=.true.; return
    end if
    rho=liquid_density_at_temperature(liquid,parcel%temperature)
    if (rho<=0) return
    liquid_mu=liquid_viscosity_at_temperature(options,parcel%temperature)
    if (liquid_mu<=0) return
    candidate=parcel; radius=0.5_dp*parcel_diameter(parcel,liquid)
    speed=sqrt(sum((gas_velocity-parcel%velocity)**2))
    total_mass=parcel%mass*parcel%multiplicity
    y=parcel%deformation; ydot=parcel%deformation_rate; elapsed=0
    do iteration=1,options%maximum_substeps
      if (elapsed>=dt.or.radius<=0.5_dp*options%minimum_diameter) exit
      td=rho*radius**2/(5.0_dp*liquid_mu)
      omega2=8.0_dp*options%surface_tension/(rho*radius**3)-1.0_dp/td**2
      if (omega2<=0) then
        y=0; ydot=0; elapsed=dt; exit
      end if
      omega=sqrt(omega2)
      we=gas_density*speed**2*radius/options%surface_tension; wer=we/12.0_dp
      step=dt-elapsed
      if (speed>0) step=min(step,0.1_dp*radius*sqrt(3.0_dp*rho/gas_density)/speed)
      if (step<=0.or.elapsed+step<=elapsed) return
      amplitude=sqrt((y-wer)**2+ydot**2/omega2)
      yn=wer+exp(-step/td)*((y-wer)*cos(omega*step)+(ydot+(y-wer)/td)/omega*sin(omega*step))
      ydn=(wer-yn)/td+exp(-step/td)*((ydot+(y-wer)/td)*cos(omega*step)-omega*(y-wer)*sin(omega*step))
      event_time=huge(1.0_dp)
      if (wer+amplitude>1.0_dp) then
        event_time=0
        if (abs(y)<1.0_dp) then
          phi=atan2(-ydot/omega,y-wer)
          if (phi<0) phi=phi+2*pi
          angle=acos(max(-1.0_dp,min(1.0_dp,(1.0_dp-wer)/amplitude)))
          if (angle<phi) then
            if (2*pi-angle>=phi) angle=-angle
            angle=angle+2*pi
          end if
          event_time=(angle-phi)/omega
        end if
      end if
      if (event_time<step.and.wer>0) then
        kbr=(2.0_dp/9.0_dp)*((sqrt(80.0_dp)-1)*(we/80.0_dp)**4+1)*omega
        if (we>=80) kbr=(2.0_dp/9.0_dp)*omega*sqrt(we)
        etb=acos(max(-1.0_dp,min(1.0_dp,1.0_dp-1.0_dp/wer)))/omega
        child=max(0.5_dp*options%minimum_diameter,radius*exp(-kbr*etb))
        ratio=child/radius
        candidate%mass=candidate%mass*ratio**3
        candidate%multiplicity=total_mass/candidate%mass
        radius=child; y=0; ydot=0
      else
        y=yn; ydot=ydn
      end if
      elapsed=min(dt,elapsed+step)
    end do
    if (iteration>options%maximum_substeps.and.elapsed<dt) return
    candidate%deformation=y; candidate%deformation_rate=ydot
    if (.not.valid_spray_parcel(candidate)) return
    parcel=candidate; ok=.true.
  end subroutine
end module
