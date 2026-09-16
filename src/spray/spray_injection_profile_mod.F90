module spray_injection_profile_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: integrate_spray_injection_profile
contains
  pure subroutine integrate_spray_injection_profile(times,rates,begin_time,end_time,mass,ok)
    real(dp), intent(in) :: times(:),rates(:),begin_time,end_time
    real(dp), intent(out) :: mass
    logical, intent(out) :: ok
    real(dp) :: lo,hi,flo,fhi,width
    integer :: knot
    mass=0;ok=.false.
    if (size(times)<2.or.size(times)>64.or.size(rates)/=size(times)) return
    if (.not.all(ieee_is_finite(times)).or..not.all(ieee_is_finite(rates))) return
    if (.not.all(ieee_is_finite([begin_time,end_time]))) return
    if (any(times<0).or.any(times>1e6_dp).or.any(rates<0).or.any(rates>1)) return
    if (any(times(2:)<=times(:size(times)-1))) return
    if (begin_time<times(1).or.end_time>times(size(times)).or.end_time<begin_time) return
    do knot=1,size(times)-1
      lo=max(begin_time,times(knot));hi=min(end_time,times(knot+1))
      if (hi<=lo) cycle
      width=times(knot+1)-times(knot)
      flo=(lo-times(knot))/width;fhi=(hi-times(knot))/width
      flo=(1-flo)*rates(knot)+flo*rates(knot+1)
      fhi=(1-fhi)*rates(knot)+fhi*rates(knot+1)
      mass=mass+.5_dp*(hi-lo)*(flo+fhi)
    end do
    ok=.true.
  end subroutine
end module
