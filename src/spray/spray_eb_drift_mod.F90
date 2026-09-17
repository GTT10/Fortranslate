module spray_eb_drift_mod
  use precision_mod, only: dp
  use spray_parcel_mod, only: spray_parcel, valid_spray_parcel
  use spray_eb_intersection_mod, only: intersect_spray_level_set
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: drift_spray_level_set
contains
  ! Force-free drift with stationary, impermeable, specular EB walls.
  ! change is particle inventory after minus before (the wall receives -change).
  ! Publish the entire population and diagnostics only after every path succeeds.
  subroutine drift_spray_level_set(phi,lower,upper,dt,maximum_hits,parcels,change,hits,ok)
    real(dp), intent(in) :: phi(0:,0:,0:),lower(3),upper(3),dt
    integer, intent(in) :: maximum_hits
    type(spray_parcel), intent(inout) :: parcels(:)
    real(dp), intent(inout) :: change(5)
    integer, intent(inout) :: hits
    logical, intent(out) :: ok
    type(spray_parcel), allocatable :: trial(:)
    real(dp) :: remaining,finish(3),fraction,position(3),normal(3),vn,delta(5),weight
    integer :: p,n,cell(3),count,allocation_status
    logical :: hit,valid
    ok=.false.
    if (.not.ieee_is_finite(dt)) return
    if (dt<0.or.dt>1e6_dp.or.maximum_hits<0.or.maximum_hits>1024.or.size(parcels)>1000000) return
    ! Validate geometry independently of fluid location for empty populations too.
    if (any(shape(phi)<2).or.any(shape(phi)>1000001)) return
    if (.not.all(ieee_is_finite(phi)).or..not.all(ieee_is_finite([lower,upper]))) return
    if (product(int(shape(phi)-1,int64))>1000000_int64) return
    if (maxval(abs(phi))>1e100_dp.or.maxval(abs([lower,upper]))>1e9_dp) return
    if (any(upper<=lower)) return
    if (any((upper-lower)/real(shape(phi)-1,dp)<1e-12_dp)) return
    allocate(trial(size(parcels)),stat=allocation_status)
    if (allocation_status/=0) return
    trial=parcels;delta=0;count=0
    hit=.false.;fraction=0;position=0;normal=0;cell=0
    do p=1,size(trial)
      if (.not.valid_spray_parcel(trial(p))) return
      if (trial(p)%mass==0) cycle
      remaining=dt;n=0
      do
        finish=trial(p)%position+remaining*trial(p)%velocity
        call intersect_spray_level_set(phi,lower,upper,trial(p)%position,finish, &
          hit,fraction,position,normal,cell,valid)
        if (.not.valid) return
        trial(p)%position=position
        if (.not.hit) exit
        if (n==maximum_hits) return
        vn=dot_product(trial(p)%velocity,normal)
        if (vn>=0) return
        trial(p)%velocity=trial(p)%velocity-2*vn*normal
        n=n+1;count=count+1
        remaining=remaining*(1-fraction)
        if (remaining==0) exit
      end do
      if (.not.valid_spray_parcel(trial(p))) return
      weight=trial(p)%mass*trial(p)%multiplicity
      delta(2:4)=delta(2:4)+weight*(trial(p)%velocity-parcels(p)%velocity)
      delta(5)=delta(5)+.5_dp*weight*(sum(trial(p)%velocity**2)-sum(parcels(p)%velocity**2))
    end do
    parcels=trial;change=delta;hits=count;ok=.true.
  end subroutine
end module
