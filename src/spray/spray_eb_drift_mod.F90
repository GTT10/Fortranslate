module spray_eb_drift_mod
  use precision_mod, only: dp
  use spray_parcel_mod, only: spray_parcel, valid_spray_parcel
  use spray_coupling_3d_mod, only: valid_spray_boundary_faces
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
  subroutine drift_spray_level_set(phi,lower,upper,dt,maximum_hits,parcels,change,hits,ok, &
      outer_faces,specific_energy,outer_change)
    real(dp), intent(in) :: phi(0:,0:,0:),lower(3),upper(3),dt
    integer, intent(in) :: maximum_hits
    type(spray_parcel), intent(inout) :: parcels(:)
    real(dp), intent(inout) :: change(5)
    integer, intent(inout) :: hits
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: outer_faces(6)
    real(dp), optional, intent(in) :: specific_energy(:)
    real(dp), optional, intent(inout) :: outer_change(5)
    real(dp) :: exterior(5),exit_fraction,test_fraction,old_velocity(3),event_change(5)
    integer :: axis,face,exit_face,dims(3)
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
    if(present(outer_faces).neqv.present(specific_energy)) return
    if(present(outer_faces).neqv.present(outer_change)) return
    if(present(outer_faces)) then
      if(.not.valid_spray_boundary_faces(outer_faces)) return
      if(size(specific_energy)/=size(parcels)) return
      if(.not.all(ieee_is_finite(specific_energy))) return
      if(any(abs(specific_energy)>1e15_dp)) return
      dims=shape(phi)-1
      if(outer_faces(1)=='periodic') then
        if(any(phi(0,:,:)/=phi(dims(1),:,:))) return
      end if
      if(outer_faces(3)=='periodic') then
        if(any(phi(:,0,:)/=phi(:,dims(2),:))) return
      end if
      if(outer_faces(5)=='periodic') then
        if(any(phi(:,:,0)/=phi(:,:,dims(3)))) return
      end if
    end if
    allocate(trial(size(parcels)),stat=allocation_status)
    if (allocation_status/=0) return
    trial=parcels;delta=0;exterior=0;count=0
    hit=.false.;fraction=0;position=0;normal=0;cell=0
    do p=1,size(trial)
      if (.not.valid_spray_parcel(trial(p))) return
      if (trial(p)%mass==0) cycle
      remaining=dt;n=0
      do
        finish=trial(p)%position+remaining*trial(p)%velocity
        exit_face=0;exit_fraction=1
        if(present(outer_faces).and.remaining>0) then
          if(any(trial(p)%position<lower).or.any(trial(p)%position>upper)) return
          do axis=1,3
            if(remaining*trial(p)%velocity(axis)==0) cycle
            if(trial(p)%velocity(axis)>0.and.finish(axis)<upper(axis)) cycle
            if(trial(p)%velocity(axis)<0.and.finish(axis)>lower(axis)) cycle
            face=2*axis
            if(trial(p)%velocity(axis)<0) face=face-1
            if(modulo(face,2)==0) then
              test_fraction=(upper(axis)-trial(p)%position(axis))/(remaining*trial(p)%velocity(axis))
            else
              test_fraction=(lower(axis)-trial(p)%position(axis))/(remaining*trial(p)%velocity(axis))
            end if
            if(test_fraction<0.or.test_fraction>exit_fraction) cycle
            if(exit_face/=0.and.test_fraction==exit_fraction) cycle
            exit_fraction=test_fraction;exit_face=face
          end do
          if(exit_face/=0) finish=max(lower,min(upper,trial(p)%position+ &
            exit_fraction*remaining*trial(p)%velocity))
        end if
        call intersect_spray_level_set(phi,lower,upper,trial(p)%position,finish, &
          hit,fraction,position,normal,cell,valid)
        if (.not.valid) return
        trial(p)%position=position
        if (.not.hit) then
          if(exit_face==0) exit
          if(n==maximum_hits) return
          axis=(exit_face+1)/2
          remaining=remaining*(1-exit_fraction)
          weight=trial(p)%mass*trial(p)%multiplicity
          old_velocity=trial(p)%velocity
          select case(trim(outer_faces(exit_face)))
          case('periodic')
            if(modulo(exit_face,2)==0) then
              trial(p)%position(axis)=lower(axis)
            else
              trial(p)%position(axis)=upper(axis)
            end if
          case('inflow','outflow')
            exterior(1)=exterior(1)-weight
            exterior(2:4)=exterior(2:4)-weight*trial(p)%velocity
            exterior(5)=exterior(5)-weight*(specific_energy(p)+.5_dp*sum(trial(p)%velocity**2))
            trial(p)%mass=0
          case default
            trial(p)%velocity(axis)=-trial(p)%velocity(axis)
            exterior(2:4)=exterior(2:4)+weight*(trial(p)%velocity-old_velocity)
          end select
          n=n+1;count=count+1
          if(remaining==0.or.trial(p)%mass==0) exit
          cycle
        end if
        if (n==maximum_hits) return
        vn=dot_product(trial(p)%velocity,normal)
        if (vn>=0) return
        old_velocity=trial(p)%velocity
        trial(p)%velocity=trial(p)%velocity-2*vn*normal
        if(present(outer_faces)) then
          weight=trial(p)%mass*trial(p)%multiplicity;event_change=0
          event_change(2:4)=weight*(trial(p)%velocity-old_velocity)
          event_change(5)=.5_dp*weight*(sum(trial(p)%velocity**2)-sum(old_velocity**2))
          delta=delta+event_change
        end if
        n=n+1;count=count+1
        remaining=remaining*(1-exit_fraction*fraction)
        if (remaining==0) exit
      end do
      if (.not.valid_spray_parcel(trial(p))) return
      if(.not.present(outer_faces)) then
        weight=trial(p)%mass*trial(p)%multiplicity
        delta(2:4)=delta(2:4)+weight*(trial(p)%velocity-parcels(p)%velocity)
        delta(5)=delta(5)+.5_dp*weight*(sum(trial(p)%velocity**2)-sum(parcels(p)%velocity**2))
      end if
    end do
    if(present(outer_change)) outer_change=exterior
    parcels=trial;change=delta;hits=count;ok=.true.
  end subroutine
end module
