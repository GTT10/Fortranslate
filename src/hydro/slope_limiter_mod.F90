module slope_limiter_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  public :: minmod2
  public :: minmod3
  public :: limited_slope
  public :: amrex_limited_slope
  public :: amrex_extdir_limited_slope

contains

  pure real(dp) function minmod2(first, second) result(value)
    real(dp), intent(in) :: first, second

    if ((first > 0.0_dp .and. second > 0.0_dp) .or. &
        (first < 0.0_dp .and. second < 0.0_dp)) then
      value = sign(min(abs(first), abs(second)), first)
    else
      value = 0.0_dp
    end if
  end function minmod2

  pure real(dp) function minmod3(first, second, third) result(value)
    real(dp), intent(in) :: first, second, third

    value = minmod2(first, minmod2(second, third))
  end function minmod3

  pure subroutine limited_slope(delta_minus, delta_plus, limiter, slope, ok)
    real(dp), intent(in) :: delta_minus, delta_plus
    character(len=*), intent(in) :: limiter
    real(dp), intent(out) :: slope
    logical, intent(out) :: ok

    select case (trim(limiter))
    case ("minmod")
      slope = minmod2(delta_minus, delta_plus)
      ok = .true.
    case ("mc")
      slope = minmod3( &
        0.5_dp * (delta_minus + delta_plus), &
        2.0_dp * delta_minus, &
        2.0_dp * delta_plus)
      ok = .true.
    case default
      slope = 0.0_dp
      ok = .false.
    end select
  end subroutine limited_slope

  pure subroutine amrex_limited_slope( &
      qm2, qm, qc, qp, qp2, flat, order, slope, ok)
    real(dp), intent(in) :: qm2, qm, qc, qp, qp2, flat
    integer, intent(in) :: order
    real(dp), intent(out) :: slope
    logical, intent(out) :: ok

    real(dp) :: dlft, drgt, dcen, dfm, dfp
    real(dp) :: dlim, dsgn, dtemp

    slope = 0.0_dp
    ok = .false.
    if (order /= 2 .and. order /= 4) return
    if (flat < 0.0_dp .or. flat > 1.0_dp) return

    ! AMReX_Slopes_K.H: the order-two branch is centered MC, while
    ! PeleC's characteristic tracing has a different order-two formula.
    if (order == 2) then
      call limited_slope(qc - qm, qp - qc, "mc", slope, ok)
      slope = flat * slope
      return
    end if

    ! Avoid underflow/overflow in sign tests for trace-field gradients.
    dfm = 0.0_dp
    dfp = 0.0_dp
    if (order == 4) then
      dlft = qm - qm2
      drgt = qc - qm
      dcen = 0.5_dp * (dlft + drgt)
      dsgn = sign(1.0_dp, dcen)
      if ((dlft>=0.and.drgt>=0).or.(dlft<=0.and.drgt<=0)) then
        dlim = 2.0_dp * min(abs(dlft), abs(drgt))
      else
        dlim = 0.0_dp
      end if
      dfm = dsgn * min(dlim, abs(dcen))

      dlft = qp - qc
      drgt = qp2 - qp
      dcen = 0.5_dp * (dlft + drgt)
      dsgn = sign(1.0_dp, dcen)
      if ((dlft>=0.and.drgt>=0).or.(dlft<=0.and.drgt<=0)) then
        dlim = 2.0_dp * min(abs(dlft), abs(drgt))
      else
        dlim = 0.0_dp
      end if
      dfp = dsgn * min(dlim, abs(dcen))
    end if

    dlft = qc - qm
    drgt = qp - qc
    dcen = 0.5_dp * (dlft + drgt)
    dsgn = sign(1.0_dp, dcen)
    if ((dlft>=0.and.drgt>=0).or.(dlft<=0.and.drgt<=0)) then
      dlim = 2.0_dp * min(abs(dlft), abs(drgt))
    else
      dlim = 0.0_dp
    end if

    dtemp = (4.0_dp / 3.0_dp) * dcen - (dfp + dfm) / 6.0_dp
    slope = flat * dsgn * min(dlim, abs(dtemp))
    ok = .true.
  end subroutine amrex_limited_slope

  ! Pinned AMReX_Slopes_K.H: ext_dir/hoextrap ghost values represent the
  ! physical boundary face, not a cell center. The order-four correction
  ! also changes the neighboring slope used one cell inside that boundary.
  pure subroutine amrex_extdir_limited_slope(qm2,qm,qc,qp,qp2,order,index, &
      domain_lower,domain_upper,extdir_lower,extdir_upper,slope,ok)
    real(dp), intent(in) :: qm2,qm,qc,qp,qp2
    integer, intent(in) :: order,index,domain_lower,domain_upper
    logical, intent(in) :: extdir_lower,extdir_upper
    real(dp), intent(out) :: slope
    logical, intent(out) :: ok
    real(dp) :: dl,dr,dc,dfm,dfp,dlim,dtemp,dsgn,side_limit

    slope=0;ok=.false.
    if(order/=2.and.order/=4) return
    if(domain_lower>domain_upper.or.index<domain_lower.or.index>domain_upper) return
    if(.not.all(ieee_is_finite([qm2,qm,qc,qp,qp2]))) return
    if(order==2) then
      dl=2*(qc-qm);dr=2*(qp-qc);dc=.5_dp*(qp-qm)
      if(extdir_lower.and.index==domain_lower) then
        dc=(qp+3*qc-4*qm)/3
      else if(extdir_upper.and.index==domain_upper) then
        dc=(4*qp-3*qc-qm)/3
      endif
      slope=min(abs(dl),abs(dc),abs(dr))
      if(.not.((dl>0.and.dr>0).or.(dl<0.and.dr<0))) slope=0
      if(dc<=0) slope=-slope
    else
      dl=qm-qm2;dr=qc-qm;dc=.5_dp*(dl+dr)
      dlim=0
      if((dl>=0.and.dr>=0).or.(dl<=0.and.dr<=0)) dlim=2*min(abs(dl),abs(dr))
      dfm=sign(1._dp,dc)*min(dlim,abs(dc))
      dl=qp-qc;dr=qp2-qp;dc=.5_dp*(dl+dr)
      dlim=0
      if((dl>=0.and.dr>=0).or.(dl<=0.and.dr<=0)) dlim=2*min(abs(dl),abs(dr))
      dfp=sign(1._dp,dc)*min(dlim,abs(dc))
      dl=qc-qm;dr=qp-qc;dc=.5_dp*(dl+dr)
      dsgn=sign(1._dp,dc);dlim=0
      if((dl>=0.and.dr>=0).or.(dl<=0.and.dr<=0)) dlim=2*min(abs(dl),abs(dr))
      dtemp=4._dp/3*dc-(dfp+dfm)/6
      if(extdir_lower.and.index==domain_lower) then
        dtemp=-16._dp/15*qm+.5_dp*qc+2._dp/3*qp-.1_dp*qp2
        dl=2*(qc-qm);dr=2*(qp-qc);dlim=0
        if((dl>=0.and.dr>=0).or.(dl<=0.and.dr<=0)) dlim=min(abs(dl),abs(dr))
        dsgn=sign(1._dp,dtemp)
      else if(extdir_lower.and.index==domain_lower+1) then
        dfm=-16._dp/15*qm2+.5_dp*qm+2._dp/3*qc-.1_dp*qp
        dl=2*(qm-qm2);dr=2*(qc-qm);side_limit=0
        if((dl>=0.and.dr>=0).or.(dl<=0.and.dr<=0)) side_limit=min(abs(dl),abs(dr))
        dfm=sign(1._dp,dfm)*min(side_limit,abs(dfm))
        dtemp=4._dp/3*dc-(dfp+dfm)/6
      endif
      ! The upstream upper-side branch is independent of the lower branch.
      if(extdir_upper.and.index==domain_upper) then
        dtemp=16._dp/15*qp-.5_dp*qc-2._dp/3*qm+.1_dp*qm2
        dl=2*(qc-qm);dr=2*(qp-qc);dlim=0
        if((dl>=0.and.dr>=0).or.(dl<=0.and.dr<=0)) dlim=min(abs(dl),abs(dr))
        dsgn=sign(1._dp,dtemp)
      else if(extdir_upper.and.index==domain_upper-1) then
        dfp=16._dp/15*qp2-.5_dp*qp-2._dp/3*qc+.1_dp*qm
        dl=2*(qp-qc);dr=2*(qp2-qp);side_limit=0
        if((dl>=0.and.dr>=0).or.(dl<=0.and.dr<=0)) side_limit=min(abs(dl),abs(dr))
        dfp=sign(1._dp,dfp)*min(side_limit,abs(dfp))
        dtemp=4._dp/3*dc-(dfp+dfm)/6
      endif
      slope=dsgn*min(dlim,abs(dtemp))
    endif
    if(.not.ieee_is_finite(slope)) then
      slope=0;return
    endif
    ok=.true.
  end subroutine amrex_extdir_limited_slope

end module slope_limiter_mod
