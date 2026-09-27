module rotating_frame_3d_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  type, public :: rotating_frame_3d
    ! Fortran axes 1/2/3 correspond to pinned PeleC axes 0/1/2.
    integer :: axis=3
    real(dp) :: angular_speed=0._dp,axis_point(3)=0._dp
  contains
    procedure :: is_valid => valid_frame
  end type
  public :: rotating_potential_3d,rotating_momentum_source_3d,convert_rotating_energy_3d
contains
  pure logical function valid_frame(frame)
    class(rotating_frame_3d), intent(in) :: frame
    valid_frame=frame%axis>=1.and.frame%axis<=3.and. &
      all(ieee_is_finite([frame%angular_speed,frame%axis_point]))
  end function

  ! Positive centrifugal potential at a cell center, in J/kg. Excluding the
  ! axial component directly avoids cancellation in |r|^2-r_axis^2.
  subroutine rotating_potential_3d(frame,position,potential,ok)
    type(rotating_frame_3d), intent(in) :: frame
    real(dp), intent(in) :: position(3)
    real(dp), intent(inout) :: potential
    logical, intent(out) :: ok
    real(dp) :: radius(3),candidate
    ok=.false.
    if(.not.frame%is_valid().or..not.all(ieee_is_finite(position))) return
    radius=position-frame%axis_point;radius(frame%axis)=0
    candidate=.5_dp*frame%angular_speed**2*sum(radius**2)
    if(.not.ieee_is_finite(candidate)) return
    potential=candidate;ok=.true.
  end subroutine

  ! Relative momentum source: -rho [2 Omega x v + Omega x (Omega x r)].
  ! In the rotating total-energy variable below, the reference rotation source
  ! has no mass, species or energy component. Do not add force dot velocity to it.
  subroutine rotating_momentum_source_3d(frame,position,density,momentum,source,ok)
    type(rotating_frame_3d), intent(in) :: frame
    real(dp), intent(in) :: position(3),density,momentum(3)
    real(dp), intent(inout) :: source(3)
    logical, intent(out) :: ok
    real(dp) :: omega(3),radius(3),velocity(3),candidate(3)
    ok=.false.
    if(.not.frame%is_valid().or..not.all(ieee_is_finite([position,density,momentum]))) return
    if(density<=0) return
    omega=0;omega(frame%axis)=frame%angular_speed
    radius=position-frame%axis_point;velocity=momentum/density
    candidate=-density*(2*cross(omega,velocity)+cross(omega,cross(omega,radius)))
    if(.not.all(ieee_is_finite(candidate))) return
    source=candidate;ok=.true.
  end subroutine

  ! Both variables use relative velocity. E_rot = rho(e+|v_rel|^2/2-Phi).
  ! Conversion to rotating storage subtracts rho*Phi; EOS recovery adds it.
  ! This is not a conversion from an inertial-velocity kinetic energy.
  subroutine convert_rotating_energy_3d(frame,position,density,input_energy,to_rotating,output_energy,ok)
    type(rotating_frame_3d), intent(in) :: frame
    real(dp), intent(in) :: position(3),density,input_energy
    logical, intent(in) :: to_rotating
    real(dp), intent(inout) :: output_energy
    logical, intent(out) :: ok
    real(dp) :: potential,candidate
    ok=.false.
    if(.not.all(ieee_is_finite([density,input_energy]))) return
    if(density<=0) return
    call rotating_potential_3d(frame,position,potential,ok)
    if(.not.ok) return
    ok=.false.
    if(to_rotating) then
      candidate=input_energy-density*potential
    else
      candidate=input_energy+density*potential
    end if
    if(.not.ieee_is_finite(candidate)) return
    output_energy=candidate;ok=.true.
  end subroutine

  pure function cross(a,b) result(c)
    real(dp), intent(in) :: a(3),b(3)
    real(dp) :: c(3)
    c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
  end function
end module rotating_frame_3d_mod
