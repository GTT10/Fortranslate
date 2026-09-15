! Ahamed impact correlations from pinned PelePhysics; see references/spray_wall_impact.json.
module spray_wall_impact_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer, parameter, public :: impact_rebound=0,impact_deposit=1,impact_dry=2
  integer, parameter, public :: impact_thermal=3,impact_spread=4,impact_wet=5
  public :: classify_wall_impact, splash_velocity_components, splash_mass_fractions
  real(dp), parameter :: pi=acos(-1.0_dp)
contains
  subroutine classify_wall_impact(kv,tstar,alpha,kcrit,film_ratio,regime,splashed,deposited,ok)
    real(dp), intent(in) :: kv,tstar,alpha,kcrit,film_ratio
    integer, intent(out) :: regime
    real(dp), intent(out) :: splashed,deposited
    logical, intent(out) :: ok
    real(dp) :: wet_critical,a,b,base,angle
    regime=-1;splashed=0;deposited=0;ok=.false.
    if (.not.all(ieee_is_finite([kv,tstar,alpha,kcrit,film_ratio]))) return
    if (kv<0.or.kv>1e8_dp.or.tstar<=0.or.tstar>1e5_dp.or.kcrit<0.or.kcrit>1e8_dp) return
    if (film_ratio<0.or.film_ratio>1e5_dp.or.alpha<pi/6.or.alpha>pi/2) return
    if (film_ratio<.01_dp) then
      if (tstar<1.1_dp) then
        if (kv<kcrit) then
          regime=impact_deposit;deposited=1
        else
          regime=impact_dry
          splashed=.75_dp*(1-exp(-.003_dp*(kv-kcrit)));deposited=1-splashed
        end if
      else
        if (kv<kcrit) then
          regime=impact_rebound
        else
          regime=impact_thermal
          splashed=.75_dp*(1-exp(-.0175_dp*(kv-kcrit)))
        end if
      end if
    else
      wet_critical=kcrit+226*film_ratio**.72_dp
      if (kv>wet_critical) then
        regime=impact_wet;a=4.5_dp;b=.2_dp
        if (film_ratio<.5_dp) then
          a=10.42_dp*film_ratio+.32_dp;b=8.75_dp*film_ratio-.107_dp
        end if
        base=a*(1-exp(-.003_dp*(kv-wet_critical)))+b
        angle=alpha*180/pi
        if (film_ratio>.1_dp) then
          if (angle>45) then
            splashed=base*(-.0305_dp*angle+3.7496_dp)
          else
            splashed=base*.0525_dp*angle
          end if
        else
          if (angle>45) then
            splashed=base*(-.0072_dp*angle+1.6389_dp)
          else
            splashed=base*(-.00046_dp*angle**2+.05_dp*angle)
          end if
        end if
        splashed=min(splashed,1.0_dp)
      else
        regime=impact_spread;deposited=1
      end if
    end if
    ! Preserve correlations, but do not pass a negative mass fraction to a caller.
    ok=splashed>=0.and.splashed<=1.and.deposited>=0.and.splashed+deposited<=1
    if (.not.ok) then
      regime=-1;splashed=0;deposited=0
    end if
  end subroutine

  subroutine splash_velocity_components(normal_speed,tangent_speed,kv,film_ratio,velocity,ok)
    real(dp), intent(in) :: normal_speed,tangent_speed,kv,film_ratio
    ! beta(0), beta(pi/2), beta(pi), psi coefficient, outward normal speed.
    real(dp), intent(out) :: velocity(5)
    logical, intent(out) :: ok
    real(dp) :: un,ut
    velocity=0;ok=.false.
    if (.not.all(ieee_is_finite([normal_speed,tangent_speed,kv,film_ratio]))) return
    if (abs(normal_speed)>1e7_dp.or.tangent_speed<0.or.tangent_speed>1e7_dp) return
    if (kv<0.or.kv>1e8_dp.or.film_ratio<0.or.film_ratio>1e5_dp) return
    un=abs(normal_speed)
    if (film_ratio<.01_dp) then
      ut=.0065_dp*kv*exp(-.004_dp*kv)*un
      velocity=[tangent_speed+ut,cos(pi/6)*tangent_speed,cos(pi/3)*tangent_speed-ut, &
        .2_dp*tangent_speed+ut,.0028_dp*kv*exp(-.0062_dp*kv)*un]
    else if (film_ratio<.1_dp) then
      if (film_ratio>.05_dp) then
        ut=(5e-5_dp*kv+.0203_dp)*un
        velocity(5)=.03_dp*tangent_speed+(1.32e-4_dp*kv+.0158_dp)*un
      else
        ut=(3.45e-5_dp*kv+.006_dp)*un
        velocity(5)=.03_dp*tangent_speed+(3.57e-5_dp*kv-.0057_dp)*un
      end if
      velocity(1:4)=[.3_dp*tangent_speed+ut,.1_dp*tangent_speed, &
        .033_dp*tangent_speed-ut,.06_dp*tangent_speed+ut]
    else
      velocity=[.1_dp*tangent_speed,.05_dp*tangent_speed,.01_dp*tangent_speed, &
        .02_dp*tangent_speed,.03_dp*tangent_speed]
    end if
    ok=all(ieee_is_finite(velocity))
  end subroutine

  subroutine splash_mass_fractions(alpha,splashed,film_ratio,fractions,ok)
    real(dp), intent(in) :: alpha,splashed,film_ratio
    ! Azimuths -pi/2, 0, pi/2, pi; fractions of the original incident mass.
    real(dp), intent(out) :: fractions(4)
    logical, intent(out) :: ok
    real(dp) :: angle,forward,side,tangent
    fractions=0;ok=.false.
    if (.not.all(ieee_is_finite([alpha,splashed,film_ratio]))) return
    if (alpha<pi/6.or.alpha>pi/2.or.splashed<0.or.splashed>1.or.film_ratio<0.or.film_ratio>1e5_dp) return
    if (film_ratio>0) then
      angle=alpha*180/pi
      if (film_ratio>.1_dp) then
        forward=(1.9e-4_dp*angle**2+.0202_dp*angle-.0107_dp)*splashed
        side=(6.4e-5_dp*angle**2-.0076_dp*angle+.4176_dp)*splashed
      else
        forward=(-2.3e-4_dp*angle**2+.0236_dp*angle)*splashed
        side=(6e-5_dp*angle**2-.0068_dp*angle+.3619_dp)*splashed
      end if
      fractions=[side,forward,side,splashed-forward-2*side]
    else
      tangent=.4_dp*tan(pi/2-alpha)
      fractions=.25_dp*splashed*[1.0_dp,1+tangent,1.0_dp,1-tangent]
    end if
    ok=all(fractions>=0).and.abs(sum(fractions)-splashed)<=32*epsilon(1.0_dp)
    if (.not.ok) fractions=0
  end subroutine
end module
