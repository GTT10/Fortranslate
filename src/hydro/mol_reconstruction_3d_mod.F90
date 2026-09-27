module mol_reconstruction_3d_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: mol_characteristic_slope_3d,mol_reconstruct_face_3d,limit_mol_positive_slope_3d
contains
  ! Primitive layout: rho, u, v, w, p, Y(:). Slope layout: two acoustic
  ! amplitudes, two transverse velocities, species contact amplitudes.
  ! Connectivity is supplied by the caller from EB flags; inactive donors
  ! are never read. No characteristic time tracing is applied in MOL.
  subroutine mol_characteristic_slope_3d(left,center,right,sound,axis,connected,slope,ok)
    real(dp), intent(in) :: left(:),center(:),right(:),sound
    integer, intent(in) :: axis
    logical, intent(in) :: connected(2)
    real(dp), intent(inout) :: slope(:)
    logical, intent(out) :: ok
    real(dp) :: dl(size(slope)),dr(size(slope)),candidate(size(slope)),mean,bound
    integer :: velocity(3),n,c
    ok=.false.;n=size(center)-5
    if(n<1.or.axis<1.or.axis>3) return
    if(size(left)/=size(center).or.size(right)/=size(center).or.size(slope)/=n+4) return
    if(.not.all(ieee_is_finite(center)).or..not.ieee_is_finite(sound)) return
    if(center(1)<=0.or.center(5)<=0.or.sound<=0) return
    call velocity_order(axis,velocity)
    dl=0;dr=0
    if(connected(1)) then
      if(.not.all(ieee_is_finite(left))) return
      call differences(left,center,dl)
    end if
    if(connected(2)) then
      if(.not.all(ieee_is_finite(right))) return
      call differences(center,right,dr)
    end if
    if(.not.all(ieee_is_finite(dl)).or..not.all(ieee_is_finite(dr))) return
    do c=1,size(slope)
      mean=.5_dp*(dl(c)+dr(c));bound=0
      ! Sign comparison avoids overflow in the reference dl*dr sign test.
      if((dl(c)>=0.and.dr(c)>=0).or.(dl(c)<=0.and.dr(c)<=0)) bound=2*min(abs(dl(c)),abs(dr(c)))
      candidate(c)=sign(min(bound,abs(mean)),mean)
    end do
    if(.not.all(ieee_is_finite(candidate))) return
    slope=candidate;ok=.true.
  contains
    subroutine differences(a,b,d)
      real(dp), intent(in) :: a(:),b(:)
      real(dp), intent(out) :: d(:)
      real(dp) :: pressure_jump,normal_jump
      pressure_jump=b(5)-a(5);normal_jump=b(velocity(1))-a(velocity(1))
      d(1)=.5_dp*pressure_jump/sound-.5_dp*center(1)*normal_jump
      d(2)=.5_dp*pressure_jump/sound+.5_dp*center(1)*normal_jump
      d(3:4)=b(velocity(2:3))-a(velocity(2:3))
      d(5:)=b(1)*b(6:)-a(1)*a(6:)-center(6:)*pressure_jump/(sound*sound)
    end subroutine
  end subroutine

  ! A single factor bounds both faces of one cell. Independent characteristic
  ! limiting does not ensure nonnegative reconstructed partial densities.
  ! This guard changes slopes only; cell averages and conservative flux assembly
  ! are untouched. It is an additional positivity safeguard beyond pinned MOL.
  subroutine limit_mol_positive_slope_3d(center,sound,axis,slope,ok)
    real(dp), intent(in) :: center(:),sound
    integer, intent(in) :: axis
    real(dp), intent(inout) :: slope(:)
    logical, intent(out) :: ok
    real(dp) :: candidate(size(slope)),partial(size(center)-5),delta(size(center)-5)
    real(dp) :: face(size(center)),theta,pressure_delta,acoustic
    integer :: c,attempt
    logical :: lower_ok,upper_ok
    ok=.false.
    if(size(center)<6.or.size(slope)/=size(center)-1) return
    if(axis<1.or.axis>3) return
    if(.not.all(ieee_is_finite(center)).or..not.all(ieee_is_finite(slope))) return
    if(.not.ieee_is_finite(sound)) return
    if(center(1)<=0.or.center(5)<=0.or.sound<=0.or.any(center(6:)<0)) return
    partial=center(1)*center(6:)
    if(sum(partial)<=0) return
    acoustic=slope(1)+slope(2)
    delta=.5_dp*(slope(5:)+center(6:)*acoustic/sound)
    pressure_delta=.5_dp*acoustic*sound
    if(.not.all(ieee_is_finite(delta)).or..not.ieee_is_finite(pressure_delta)) return
    theta=1
    do c=1,size(partial)
      if(abs(delta(c))>partial(c)) theta=min(theta,partial(c)/abs(delta(c)))
    end do
    if(abs(pressure_delta)>=center(5)) theta=min(theta,center(5)/abs(pressure_delta))
    ! Move inside an active bound, then verify the actual floating-point faces.
    if(theta<1.or.abs(pressure_delta)>=center(5)) theta=theta*(1-16*epsilon(1._dp))
    do attempt=1,55
      candidate=theta*slope
      call mol_reconstruct_face_3d(center,sound,axis,-1,candidate,face,lower_ok)
      call mol_reconstruct_face_3d(center,sound,axis,1,candidate,face,upper_ok)
      if(lower_ok.and.upper_ok) then
        slope=candidate;ok=.true.;return
      end if
      theta=.5_dp*theta
    end do
  end subroutine limit_mol_positive_slope_3d

  ! side = +1 for a cell's upper face, -1 for its lower face.
  ! Reconstruct partial densities, sum them, then normalize mass fractions.
  subroutine mol_reconstruct_face_3d(center,sound,axis,side,slope,face,ok)
    real(dp), intent(in) :: center(:),sound,slope(:)
    integer, intent(in) :: axis,side
    real(dp), intent(inout) :: face(:)
    logical, intent(out) :: ok
    real(dp) :: candidate(size(center)),weight,acoustic
    integer :: velocity(3),n
    ok=.false.;n=size(center)-5
    if(n<1.or.axis<1.or.axis>3.or.abs(side)/=1) return
    if(size(face)/=size(center).or.size(slope)/=n+4) return
    if(.not.all(ieee_is_finite(center)).or..not.all(ieee_is_finite(slope))) return
    if(.not.ieee_is_finite(sound)) return
    if(center(1)<=0.or.center(5)<=0.or.sound<=0) return
    call velocity_order(axis,velocity)
    weight=.5_dp*real(side,dp);acoustic=slope(1)+slope(2)
    candidate=center
    candidate(velocity(1))=center(velocity(1))+weight*(slope(2)-slope(1))/center(1)
    candidate(velocity(2:3))=center(velocity(2:3))+weight*slope(3:4)
    candidate(5)=center(5)+weight*acoustic*sound
    candidate(6:)=center(1)*center(6:)+weight*(slope(5:)+center(6:)*acoustic/sound)
    candidate(1)=sum(candidate(6:))
    if(.not.all(ieee_is_finite(candidate))) return
    if(candidate(1)<=0.or.candidate(5)<=0.or.any(candidate(6:)<0)) return
    candidate(6:)=candidate(6:)/candidate(1)
    face=candidate;ok=.true.
  end subroutine

  pure subroutine velocity_order(axis,indices)
    integer, intent(in) :: axis
    integer, intent(out) :: indices(3)
    select case(axis)
    case(1)
      indices=[2,3,4]
    case(2)
      indices=[3,2,4]
    case(3)
      indices=[4,2,3]
    end select
  end subroutine
end module mol_reconstruction_3d_mod
