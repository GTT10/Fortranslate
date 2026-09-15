! Ahamed impact correlations from pinned PelePhysics; see references/spray_wall_impact.json.
module spray_wall_impact_mod
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use spray_parcel_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer, parameter, public :: impact_rebound=0,impact_deposit=1,impact_dry=2
  integer, parameter, public :: impact_thermal=3,impact_spread=4,impact_wet=5
  public :: classify_wall_impact, splash_velocity_components, splash_mass_fractions
  type, public :: spray_wall_film
    real(dp) :: position(3)=0,normal(3)=0,mass=0,temperature=0,diameter=0,height=0
    integer :: material_id=1
  end type
  public :: resolve_wall_impact
  public :: construct_splash_children, construct_wall_film, wall_film_totals
  real(dp), parameter :: pi=acos(-1.0_dp)
contains
  pure function cross(a,b) result(c)
    real(dp), intent(in) :: a(3),b(3)
    real(dp) :: c(3)
    c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
  end function

  subroutine resolve_wall_impact(parent,liquid,location,normal,kv,tstar,kcrit,film_ratio,weber,reynolds, &
      contact_angle,first_id,children,film,wall,regime,ok)
    type(spray_parcel), intent(inout) :: parent
    type(spray_liquid), intent(in) :: liquid
    real(dp), intent(in) :: location(3),normal(3),kv,tstar,kcrit,film_ratio,weber,reynolds,contact_angle
    integer(int64), intent(in) :: first_id
    type(spray_parcel), allocatable, intent(inout) :: children(:)
    type(spray_wall_film), intent(inout) :: film
    real(dp), intent(out) :: wall(5)
    integer, intent(out) :: regime
    logical, intent(out) :: ok
    type(spray_parcel) :: next
    type(spray_parcel), allocatable :: spawned(:)
    type(spray_wall_film) :: deposited
    real(dp) :: un,ut,speed,alpha,ms,dm,restitution,remaining,psi(3),beta(3),helper(3),budget(5),total(5)
    integer :: kind,i,axis
    logical :: valid
    wall=0;regime=-1;ok=.false.
    if (.not.valid_spray_parcel(parent).or..not.valid_spray_liquid(liquid)) return
    if (.not.all(ieee_is_finite([location,normal,weber,reynolds,contact_angle]))) return
    if (maxval(abs(location))>1e9_dp.or.maxval(abs(normal))>1) return
    if (abs(sum(normal**2)-1)>1e-12_dp.or.parent%mass<=0) return
    if (liquid_heat_capacity_at_temperature(liquid,parent%temperature)<=0.or. &
        liquid_density_at_temperature(liquid,parent%temperature)<=0) return
    un=dot_product(parent%velocity,normal);speed=norm2(parent%velocity)
    if (un>=0.or.speed<=0) return
    ut=norm2(parent%velocity-un*normal);alpha=max(pi/6,asin(min(1.0_dp,abs(un)/speed)))
    call classify_wall_impact(kv,tstar,alpha,kcrit,film_ratio,kind,ms,dm,valid)
    if (.not.valid) return
    ! The reference wet branch leaves dm unset. Close its incident-mass budget
    ! by depositing the unsplashed fraction; no entrainment from existing film.
    if (kind==impact_wet) dm=1-ms
    next=parent;allocate(spawned(0))
    if (ms>0) then
      call construct_splash_children(parent,liquid,location,normal,kv,film_ratio,ms,first_id,spawned,budget,valid)
      if (.not.valid) return
    end if
    if (dm>0) then
      call construct_wall_film(parent,liquid,location,normal,dm,weber,reynolds,contact_angle,deposited,budget,valid)
      if (.not.valid) return
    end if
    remaining=max(0.0_dp,1-ms-dm);next%mass=parent%mass*remaining
    next%position=location
    if (remaining>0) then
      select case(kind)
      case(impact_rebound)
        restitution=.823_dp*exp(-.0019835_dp*kv**1.6_dp)
        next%velocity=parent%velocity-(1+restitution)*un*normal
      case(impact_thermal)
        psi=cross(-parent%velocity,normal)
        if (norm2(psi)<1e-5_dp) psi=cross(-parent%velocity+1.0_dp,normal)
        if (norm2(psi)<1e-12_dp) then
          axis=minloc(abs(normal),dim=1);helper=0;helper(axis)=1;psi=cross(helper,normal)
        end if
        psi=psi/norm2(psi);beta=cross(psi,normal);beta=beta/norm2(beta)
        next%velocity=(1.5_dp*sqrt(3.0_dp)*ut/pi)*beta+ &
          (.0028_dp*kv*exp(-.0062_dp*kv)*speed)*normal
        next%deformation=0;next%deformation_rate=0;next%rt_time=-1
      case default
        return
      end select
      next%position=location+parcel_diameter(next,liquid)*normal
    end if
    if (.not.valid_spray_parcel(next)) return
    total=parcel_totals(next,liquid)
    do i=1,size(spawned)
      total=total+parcel_totals(spawned(i),liquid)
    end do
    if (deposited%mass>0) total=total+wall_film_totals(deposited,liquid)
    budget=parcel_totals(parent,liquid)-total
    if (.not.all(ieee_is_finite(budget))) return
    if (abs(budget(1))>1e-12_dp*parent%mass*parent%multiplicity) return
    budget(1)=0
    parent=next;film=deposited;call move_alloc(spawned,children)
    wall=budget;regime=kind;ok=.true.
  end subroutine

  subroutine construct_splash_children(parent,liquid,location,normal,kv,film_ratio,fraction,first_id,children,wall,ok)
    type(spray_parcel), intent(in) :: parent
    type(spray_liquid), intent(in) :: liquid
    real(dp), intent(in) :: location(3),normal(3),kv,film_ratio,fraction
    integer(int64), intent(in) :: first_id
    type(spray_parcel), allocatable, intent(inout) :: children(:)
    real(dp), intent(out) :: wall(5)
    logical, intent(out) :: ok
    type(spray_parcel), allocatable :: candidate(:)
    real(dp) :: un,ut,speed,alpha,psi(3),beta(3),helper(3),masses(4),v(5),ub,up,angle,total(5)
    integer :: i,axis
    wall=0;ok=.false.
    if (.not.valid_spray_parcel(parent).or..not.valid_spray_liquid(liquid)) return
    if (.not.all(ieee_is_finite([location,normal,fraction]))) return
    if (maxval(abs(location))>1e9_dp.or.maxval(abs(normal))>1.0_dp) return
    if (abs(sum(normal**2)-1)>1e-12_dp.or.fraction<=0.or.fraction>1.or.parent%mass<=0) return
    if (first_id<=parent%id.or.first_id>huge(first_id)-4_int64) return
    if (liquid_heat_capacity_at_temperature(liquid,parent%temperature)<=0) return
    if (liquid_density_at_temperature(liquid,parent%temperature)<=0) return
    un=dot_product(parent%velocity,normal);speed=norm2(parent%velocity)
    if (un>=0.or.speed<=0) return
    ut=norm2(parent%velocity-un*normal)
    alpha=max(pi/6,asin(min(1.0_dp,abs(un)/speed)))
    psi=cross(-parent%velocity,normal)
    if (norm2(psi)<1e-5_dp) then
      ! Reference tangent fallback; use a coordinate axis if its [1,1,1] offset is parallel too.
      psi=cross(-parent%velocity+1.0_dp,normal)
      if (norm2(psi)<1e-12_dp) then
        axis=minloc(abs(normal),dim=1);helper=0;helper(axis)=1
        psi=cross(helper,normal)
      end if
    end if
    psi=psi/norm2(psi);beta=cross(psi,normal);beta=beta/norm2(beta)
    call splash_velocity_components(un,ut,kv,film_ratio,v,ok)
    if (.not.ok) return
    ok=.false.
    if (v(5)<0) return
    call splash_mass_fractions(alpha,fraction,film_ratio,masses,ok)
    if (.not.ok) return
    ok=.false.;allocate(candidate(4));total=0
    do i=1,4
      candidate(i)=parent
      candidate(i)%id=first_id+int(i-1,int64);candidate(i)%parent_id=parent%id
      candidate(i)%mass=parent%mass*masses(i)
      candidate(i)%deformation=0;candidate(i)%deformation_rate=0;candidate(i)%rt_time=-1
      ub=v(2)
      if (i==2) ub=v(1)
      if (i==4) ub=v(3)
      angle=.5_dp*pi*real(i-2,dp);up=v(4)*sin(angle)
      candidate(i)%velocity=v(5)*normal+ub*beta+up*psi
      candidate(i)%position=location+parcel_diameter(candidate(i),liquid)*normal
      if (.not.valid_spray_parcel(candidate(i))) return
      total=total+parcel_totals(candidate(i),liquid)
    end do
    wall=fraction*parcel_totals(parent,liquid)-total
    ! Only impulse/work goes to the wall in this operation; mass stays in children.
    if (abs(wall(1))>1e-12_dp*parent%mass*parent%multiplicity) then
      wall=0;return
    end if
    wall(1)=0
    if (.not.all(ieee_is_finite(wall))) then
      wall=0;return
    end if
    call move_alloc(candidate,children);ok=.true.
  end subroutine

  pure function wall_film_totals(film,liquid) result(totals)
    type(spray_wall_film), intent(in) :: film
    type(spray_liquid), intent(in) :: liquid
    real(dp) :: totals(5)
    totals=0;totals(1)=film%mass
    totals(5)=film%mass*liquid_energy(liquid,film%temperature)
  end function

  subroutine construct_wall_film(parent,liquid,location,normal,fraction,weber,reynolds,contact_angle,film,wall,ok)
    type(spray_parcel), intent(in) :: parent
    type(spray_liquid), intent(in) :: liquid
    real(dp), intent(in) :: location(3),normal(3),fraction,weber,reynolds,contact_angle
    type(spray_wall_film), intent(inout) :: film
    real(dp), intent(out) :: wall(5)
    logical, intent(out) :: ok
    type(spray_wall_film) :: candidate
    real(dp) :: rho,volume,combined_diameter,combined_we,combined_re,denominator,drop_diameter
    wall=0;ok=.false.
    if (.not.valid_spray_parcel(parent).or..not.valid_spray_liquid(liquid)) return
    if (.not.all(ieee_is_finite([location,normal,fraction,weber,reynolds,contact_angle]))) return
    if (maxval(abs(location))>1e9_dp.or.maxval(abs(normal))>1) return
    if (abs(sum(normal**2)-1)>1e-12_dp.or.parent%mass<=0.or.fraction<=0.or.fraction>1) return
    if (weber<0.or.weber>1e12_dp.or.reynolds<=0.or.reynolds>1e12_dp) return
    if (contact_angle<1e-6_dp.or.contact_angle>pi-1e-6_dp) return
    rho=liquid_density_at_temperature(liquid,parent%temperature)
    if (rho<=0.or.liquid_heat_capacity_at_temperature(liquid,parent%temperature)<=0) return
    candidate%mass=fraction*parent%mass*parent%multiplicity;volume=candidate%mass/rho
    candidate%position=location;candidate%normal=normal;candidate%temperature=parent%temperature
    candidate%material_id=parent%material_id
    combined_diameter=(6*volume/pi)**(1.0_dp/3)
    drop_diameter=parcel_diameter(parent,liquid)
    combined_we=weber*combined_diameter/drop_diameter;combined_re=reynolds*combined_diameter/drop_diameter
    ! 6(1-cos(theta))/sin(theta)^2 = 6/(1+cos(theta)).
    denominator=6/(1+cos(contact_angle))+4*combined_we/sqrt(combined_re)
    candidate%diameter=combined_diameter*sqrt((12+combined_we)/denominator)
    candidate%height=volume/(.25_dp*pi*candidate%diameter**2)
    if (.not.all(ieee_is_finite([candidate%diameter,candidate%height]))) return
    if (candidate%diameter<=0.or.candidate%height<=0) return
    wall=fraction*parcel_totals(parent,liquid)-wall_film_totals(candidate,liquid);wall(1)=0
    if (.not.all(ieee_is_finite(wall))) then
      wall=0;return
    end if
    film=candidate;ok=.true.
  end subroutine

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
