! Ideal liquid mixing rules adapted from PelePhysics SprayPropertiesPeleMP.H,
! revision 6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b.
! Copyright/redistribution: references/pelephysics-spray-license.txt.
module spray_mixture_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use spray_parcel_mod
  implicit none
  private
  public :: liquid_mixture_properties, liquid_mixture_temperature, advance_mixture_heating
contains
  pure logical function valid_composition(components,y) result(ok)
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(in) :: y(:)
    integer :: i,j
    ok=.false.
    if (size(y)<1.or.size(y)>16.or.size(components)/=size(y)) return
    if (.not.all(ieee_is_finite(y))) return
    if (any(y<0).or.any(y>1)) return
    if (abs(sum(y)-1)>32*epsilon(1.0_dp)*size(y)) return
    do i=1,size(y)
      if (.not.valid_spray_liquid(components(i))) return
      do j=1,i-1
        if (components(i)%vapor_index==components(j)%vapor_index) return
      end do
    end do
    ok=.true.
  end function

  pure subroutine liquid_mixture_properties(components,y,temperature,density,cp,energy,ok)
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(in) :: y(:),temperature
    real(dp), intent(inout) :: density,cp,energy
    logical, intent(out) :: ok
    real(dp) :: volume,candidate_cp,candidate_e,rhoi,cpi,ei
    integer :: i
    ok=.false.
    if (.not.valid_composition(components,y).or..not.ieee_is_finite(temperature)) return
    if (temperature<1.or.temperature>1e5_dp) return
    volume=0;candidate_cp=0;candidate_e=0
    do i=1,size(y)
      if (y(i)==0) cycle
      rhoi=liquid_density_at_temperature(components(i),temperature)
      cpi=liquid_heat_capacity_at_temperature(components(i),temperature)
      if (rhoi<=0.or.cpi<=0) return
      ei=liquid_energy(components(i),temperature)
      if (.not.ieee_is_finite(ei)) return
      volume=volume+y(i)/rhoi
      candidate_cp=candidate_cp+y(i)*cpi
      candidate_e=candidate_e+y(i)*ei
    end do
    if (.not.all(ieee_is_finite([volume,candidate_cp,candidate_e])).or.volume<=0) return
    density=1/volume;cp=candidate_cp;energy=candidate_e;ok=.true.
  end subroutine

  pure subroutine liquid_mixture_temperature(components,y,energy,temperature,ok)
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(in) :: y(:),energy
    real(dp), intent(inout) :: temperature
    logical, intent(out) :: ok
    real(dp) :: lower,upper,elo,ehi,t,rho,cp,e,step
    integer :: i,iteration
    logical :: valid
    ok=.false.
    if (.not.valid_composition(components,y)) return
    if (.not.all(ieee_is_finite([temperature,energy]))) return
    lower=1;upper=1e5_dp
    do i=1,size(y)
      if (y(i)==0) cycle
      if (components(i)%density_model=='polynomial') then
        lower=max(lower,components(i)%density_temperature_range(1))
        upper=min(upper,components(i)%density_temperature_range(2))
      end if
      if (components(i)%heat_capacity_model=='polynomial') then
        lower=max(lower,components(i)%heat_capacity_temperature_range(1))
        upper=min(upper,components(i)%heat_capacity_temperature_range(2))
      end if
    end do
    if (lower>upper) return
    call liquid_mixture_properties(components,y,lower,rho,cp,elo,valid)
    if (.not.valid) return
    call liquid_mixture_properties(components,y,upper,rho,cp,ehi,valid)
    if (.not.valid.or.energy<elo.or.energy>ehi) return
    if (energy==elo) then
      temperature=lower;ok=.true.;return
    else if (energy==ehi) then
      temperature=upper;ok=.true.;return
    end if
    t=max(lower,min(upper,temperature))
    do iteration=1,100
      call liquid_mixture_properties(components,y,t,rho,cp,e,valid)
      if (.not.valid) return
      if (abs(e-energy)<=8*epsilon(1.0_dp)*max(1.0_dp,abs(e),abs(energy))) then
        temperature=t;ok=.true.;return
      end if
      if (e<energy) then
        lower=t
      else
        upper=t
      end if
      step=t+(energy-e)/cp
      if (step<=lower.or.step>=upper) step=lower+.5_dp*(upper-lower)
      if (step==t) return
      t=step
    end do
  end subroutine

  subroutine advance_mixture_heating(parcel,components,y,heat,gas_energy,ok)
    type(spray_parcel), intent(inout) :: parcel
    type(spray_liquid), intent(in) :: components(:)
    real(dp), intent(in) :: y(:),heat
    real(dp), intent(out) :: gas_energy
    logical, intent(out) :: ok
    real(dp) :: mass,rho,cp,old_e,new_e,t,target
    logical :: valid
    ok=.false.;gas_energy=0
    if (.not.valid_spray_parcel(parcel).or..not.ieee_is_finite(heat)) return
    mass=parcel%mass*parcel%multiplicity
    if (mass<=tiny(1.0_dp)) return
    if (abs(heat)>min(1e24_dp,1e15_dp*mass)) return
    call liquid_mixture_properties(components,y,parcel%temperature,rho,cp,old_e,valid)
    if (.not.valid) return
    target=old_e+heat/mass;t=parcel%temperature
    call liquid_mixture_temperature(components,y,target,t,valid)
    if (.not.valid) return
    call liquid_mixture_properties(components,y,t,rho,cp,new_e,valid)
    if (.not.valid) return
    ! Debit the actual liquid-energy increment, including inversion roundoff.
    gas_energy=mass*(old_e-new_e);parcel%temperature=t;ok=.true.
  end subroutine
end module
