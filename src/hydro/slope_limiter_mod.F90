module slope_limiter_mod
  use precision_mod, only: dp
  implicit none
  private

  public :: minmod2
  public :: minmod3
  public :: limited_slope

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

end module slope_limiter_mod
