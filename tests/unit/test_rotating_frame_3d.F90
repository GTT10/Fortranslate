program test_rotating_frame_3d
  use precision_mod, only: dp
  use rotating_frame_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  type(rotating_frame_3d) :: frame,shifted
  real(dp) :: position(3),source(3),expected(3),momentum(3),phi,energy,back,other,nan
  integer :: axis,a,b,sign_index,checks
  logical :: ok
  checks=0;nan=ieee_value(0._dp,ieee_quiet_nan)
  do axis=1,3
    a=mod(axis,3)+1;b=mod(axis+1,3)+1
    do sign_index=-1,1,2
      frame%axis=axis;frame%angular_speed=5*real(sign_index,dp);frame%axis_point=[1._dp,-2._dp,3._dp]
      position=frame%axis_point;position(a)=position(a)+3;position(b)=position(b)+4
      position(axis)=position(axis)+1000
      call rotating_potential_3d(frame,position,phi,ok)
      call require(ok.and.phi==312.5_dp,'axis-independent centrifugal potential')
      call rotating_momentum_source_3d(frame,position,2._dp,[0._dp,0._dp,0._dp],source,ok)
      expected=0;expected(a)=150;expected(b)=200
      call require(ok.and.all(source==expected),'outward centrifugal force')
      momentum=0;momentum(a)=6;momentum(b)=8
      call rotating_momentum_source_3d(frame,frame%axis_point,2._dp,momentum,source,ok)
      expected=0;expected(a)=16*frame%angular_speed;expected(b)=-12*frame%angular_speed
      call require(ok.and.all(source==expected),'Coriolis direction and angular-speed sign')
      call require(dot_product(source,momentum)==0,'Coriolis does no work')
      call convert_rotating_energy_3d(frame,position,2._dp,10000._dp,.true.,energy,ok)
      call require(ok.and.energy==9375,'rotating total energy')
      call convert_rotating_energy_3d(frame,position,2._dp,energy,.false.,back,ok)
      call require(ok.and.back==10000,'EOS energy recovery')
      call convert_rotating_energy_3d(frame,position,2._dp,100._dp,.true.,energy,ok)
      call require(ok.and.energy==-525,'negative rotating energy is valid')
      call convert_rotating_energy_3d(frame,position,2._dp,energy,.false.,back,ok)
      call require(ok.and.back==100,'negative rotating energy recovery')
      shifted=frame;shifted%axis_point=frame%axis_point+[2._dp,4._dp,-8._dp]
      call rotating_potential_3d(shifted,position+[2._dp,4._dp,-8._dp],other,ok)
      call require(ok.and.other==phi,'translation invariance')
    end do
  end do
  frame%angular_speed=0
  call rotating_momentum_source_3d(frame,position,2._dp,momentum,source,ok)
  call require(ok.and.all(source==0),'zero rotation source')
  call convert_rotating_energy_3d(frame,position,2._dp,10000._dp,.true.,energy,ok)
  call require(ok.and.energy==10000,'zero rotation energy identity')
  frame%axis=0;source=19;phi=23;energy=29
  call rotating_momentum_source_3d(frame,position,2._dp,momentum,source,ok)
  call require(.not.ok.and.all(source==19),'invalid axis source rollback')
  call rotating_potential_3d(frame,position,phi,ok)
  call require(.not.ok.and.phi==23,'invalid axis potential rollback')
  call convert_rotating_energy_3d(frame,position,2._dp,10000._dp,.true.,energy,ok)
  call require(.not.ok.and.energy==29,'invalid axis energy rollback')
  frame%axis=3
  call rotating_momentum_source_3d(frame,position,-1._dp,momentum,source,ok)
  call require(.not.ok.and.all(source==19),'invalid density rollback')
  frame%angular_speed=nan
  call rotating_potential_3d(frame,position,phi,ok)
  call require(.not.ok.and.phi==23,'nonfinite rotation rollback')
  print *, 'PASS: rotating-frame checks:',checks
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
    checks=checks+1
  end subroutine
end program
