! Isolated scalar slope experiment; not connected to production EB halos.
! Sekora & Colella, arXiv:0903.4200v2, section 2.1 (C_VL = 1.25).
! No positivity, multidimensional EOS, or fixed-PeleC parity claim.
module eb_extremum_slope_mod
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  public :: extremum_slope,bounded_quartic_slope,continuous_quartic_slope
contains
  pure real(real64) function slope_bound(q,roundoff_safe) result(dl)
    real(real64), intent(in) :: q(-2:2)
    logical, optional, intent(in) :: roundoff_safe
    logical :: extremum
    real(real64) :: dm,dp,dc,dmm,dpp,d2m,d2c,d2p,s2,d2lim,noise
    dm=q(0)-q(-1);dp=q(1)-q(0);dc=0.5_real64*(q(1)-q(-1))
    dmm=q(-1)-q(-2);dpp=q(2)-q(1)
    dl=2*min(abs(dm),abs(dp))
    extremum=min(dm*dp,dmm*dpp)<0
    if(present(roundoff_safe)) then
      if(roundoff_safe) then
        ! Experimental branch guard: signs below the local floating-point
        ! resolution are not evidence of an extremum. Keep the magnitude
        ! bound unchanged; never add a dimensional absolute tolerance.
        noise=8*epsilon(1._real64)*maxval(abs(q))
        extremum=(dm>noise.and.dp< -noise).or.(dm< -noise.and.dp>noise).or. &
          (dmm>noise.and.dpp< -noise).or.(dmm< -noise.and.dpp>noise)
      endif
    endif
    if(extremum) then
      d2m=dm-dmm;d2c=dp-dm;d2p=dpp-dp
      s2=sign(1._real64,d2c)
      d2lim=min(abs(d2c),max(s2*d2m,0._real64),max(s2*d2p,0._real64))
      if(s2*dc<0) then
        dl=min(1.875_real64*d2lim,2*abs(dm))
      else
        dl=min(1.875_real64*d2lim,2*abs(dp))
      endif
    endif
  end function

  pure real(real64) function extremum_slope(q) result(slope)
    real(real64), intent(in) :: q(-2:2)
    real(real64) :: dc
    dc=0.5_real64*(q(1)-q(-1))
    slope=sign(min(abs(dc),slope_bound(q)),dc)
  end function

  ! Experimental adaptation: bound the quartic half-cell difference using
  ! the EP magnitude bound. This is NOT the published MUSCL/PPM algorithm.
  pure real(real64) function bounded_quartic_slope(q) result(slope)
    real(real64), intent(in) :: q(-2:2)
    real(real64) :: dc,high
    dc=0.5_real64*(q(1)-q(-1))
    ! Difference form avoids cancellation of the constant background.
    high=(22*(q(1)-q(-1))-3*(q(2)-q(-2)))/32
    slope=0
    if(high*dc<0) return
    slope=sign(min(abs(high),slope_bound(q,roundoff_safe=.true.)),high)
  end function
  ! Experimental continuous envelope; no production or published-method claim.
  pure real(real64) function continuous_quartic_slope(q) result(slope)
    real(real64), intent(in) :: q(-2:2)
    real(real64) :: dm,dp,dc,dmm,dpp,d2m,d2c,d2p,s2,curvature,monotone,bound,high
    dm=q(0)-q(-1);dp=q(1)-q(0);dc=.5_real64*(q(1)-q(-1))
    dmm=q(-1)-q(-2);dpp=q(2)-q(1)
    monotone=0
    if((dm>0.and.dp>0).or.(dm<0.and.dp<0)) monotone=2*min(abs(dm),abs(dp))
    d2m=dm-dmm;d2c=dp-dm;d2p=dpp-dp;s2=sign(1._real64,d2c)
    curvature=1.875_real64*min(abs(d2c),max(s2*d2m,0._real64),max(s2*d2p,0._real64))
    if(s2*dc<0) then
      curvature=min(curvature,2*abs(dm))
    else
      curvature=min(curvature,2*abs(dp))
    endif
    bound=max(monotone,curvature)
    high=(22*(q(1)-q(-1))-3*(q(2)-q(-2)))/32
    slope=0
    ! A minmod-style cap vanishes continuously when the centered slope does.
    if((high>0.and.dc>0).or.(high<0.and.dc<0)) &
      slope=sign(min(abs(high),2*abs(dc),bound),high)
  end function
end module
