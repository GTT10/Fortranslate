! Breakup equations adapted from PelePhysics TABBreakup.H and ReitzKHRT.H,
! revision 6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b.
! Copyright and redistribution terms: references/pelephysics-spray-license.txt.
module spray_breakup_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use spray_parcel_mod
  implicit none
  private
  public :: advance_etab_breakup, advance_khrt_breakup, advance_spray_breakup, split_spray_parcels
contains
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
  subroutine advance_khrt_breakup(parcel,liquid,options,gas_density,gas_velocity,reynolds_number,dt,ok)
    type(spray_parcel), intent(inout) :: parcel
    type(spray_liquid), intent(in) :: liquid
    type(spray_options), intent(in) :: options
    real(dp), intent(in) :: gas_density,gas_velocity(3),reynolds_number,dt
    logical, intent(out) :: ok
    type(spray_parcel) :: candidate
    real(dp) :: radius,child,speed,re,cd,we,z,t,lambda,omega,tau,acceleration,momdiff,rate,total
    logical :: rt_breakup
    real(dp), parameter :: pi=acos(-1.0_dp)
    ok=.false.
    if (.not.valid_spray_parcel(parcel).or..not.valid_spray_liquid(liquid).or..not.valid_spray_options(options)) return
    if (trim(options%breakup_model)/='khrt') return
    if (.not.all(ieee_is_finite([gas_density,gas_velocity,reynolds_number,dt]))) return
    if (gas_density<=0.or.gas_density>1e8_dp.or.reynolds_number<0.or.reynolds_number>1e30_dp) return
    if (maxval(abs(gas_velocity))>1e7_dp.or.dt<0.or.dt>1e6_dp) return
    if (parcel%mass<=0.or.dt==0) then
      ok=.true.;return
    end if
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
    z=options%liquid_viscosity/sqrt(liquid%density*radius*options%surface_tension)
    t=z*sqrt(we);rt_breakup=.false.
    if (candidate%rt_time>=0) then
      acceleration=.375_dp*cd*gas_density*speed**2/(liquid%density*radius)
      momdiff=acceleration*(liquid%density-gas_density)
      if (momdiff>0) then
        omega=sqrt(2/(3*sqrt(3*options%surface_tension))*momdiff*sqrt(momdiff)/(liquid%density+gas_density))
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
        sqrt(options%surface_tension/(liquid%density*radius**3))
      tau=3.726_dp*options%khrt_b1*radius/(lambda*omega)
      if (options%khrt_b0*lambda<radius.and.we>6) then
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
    parcel=candidate;ok=.true.
  end subroutine

  subroutine split_spray_parcels(parcels,maximum_multiplicity,maximum_parcels,ok)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), intent(in) :: maximum_multiplicity
    integer, intent(in) :: maximum_parcels
    logical, intent(out) :: ok
    type(spray_parcel), allocatable :: candidate(:)
    integer, allocatable :: pieces(:)
    integer :: p,n,total,extra,cursor
    integer(int64) :: last_id
    real(dp) :: weight
    ok=.false.
    if (.not.allocated(parcels)) return
    if (.not.ieee_is_finite(maximum_multiplicity).or.maximum_multiplicity<0.or.maximum_multiplicity>1e18_dp) return
    if (maximum_multiplicity>0.and.maximum_multiplicity<1e-12_dp) return
    if (maximum_parcels<1.or.maximum_parcels>1000000.or.size(parcels)>maximum_parcels) return
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
    if (last_id>huge(last_id)-int(total-size(parcels),int64)) return
    allocate(candidate(total)); candidate(:size(parcels))=parcels;cursor=size(parcels)
    do p=1,size(parcels)
      if (pieces(p)==1) cycle
      weight=parcels(p)%multiplicity/real(pieces(p),dp)
      candidate(p)%multiplicity=weight
      do n=2,pieces(p)
        cursor=cursor+1;candidate(cursor)=parcels(p)
        last_id=last_id+1;candidate(cursor)%id=last_id
        candidate(cursor)%parent_id=parcels(p)%id
        candidate(cursor)%multiplicity=weight
        if (n==pieces(p)) candidate(cursor)%multiplicity=parcels(p)%multiplicity-real(n-1,dp)*weight
        if (.not.valid_spray_parcel(candidate(cursor))) return
      end do
    end do
    call move_alloc(candidate,parcels);ok=.true.
  end subroutine

  subroutine advance_etab_breakup(parcel,liquid,options,gas_density,gas_velocity,dt,ok)
    type(spray_parcel), intent(inout) :: parcel
    type(spray_liquid), intent(in) :: liquid
    type(spray_options), intent(in) :: options
    real(dp), intent(in) :: gas_density,gas_velocity(3),dt
    logical, intent(out) :: ok
    type(spray_parcel) :: candidate
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
    candidate=parcel; radius=0.5_dp*parcel_diameter(parcel,liquid)
    speed=sqrt(sum((gas_velocity-parcel%velocity)**2))
    total_mass=parcel%mass*parcel%multiplicity
    y=parcel%deformation; ydot=parcel%deformation_rate; elapsed=0
    do iteration=1,options%maximum_substeps
      if (elapsed>=dt.or.radius<=0.5_dp*options%minimum_diameter) exit
      td=liquid%density*radius**2/(5.0_dp*options%liquid_viscosity)
      omega2=8.0_dp*options%surface_tension/(liquid%density*radius**3)-1.0_dp/td**2
      if (omega2<=0) then
        y=0; ydot=0; elapsed=dt; exit
      end if
      omega=sqrt(omega2)
      we=gas_density*speed**2*radius/options%surface_tension; wer=we/12.0_dp
      step=dt-elapsed
      if (speed>0) step=min(step,0.1_dp*radius*sqrt(3.0_dp*liquid%density/gas_density)/speed)
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
