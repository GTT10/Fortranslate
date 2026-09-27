program test_mol_reconstruction_3d
  use precision_mod, only: dp
  use mol_reconstruction_3d_mod, only: mol_characteristic_slope_3d,mol_reconstruct_face_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  real(dp) :: q(7),left(7),right(7),slope(6),expected(6),face(7),target(7),jump(7),nan
  integer :: axis,side,velocity(3),checks
  logical :: ok
  nan=ieee_value(0._dp,ieee_quiet_nan);checks=0
  do axis=1,3
    select case(axis)
    case(1)
      velocity=[2,3,4]
    case(2)
      velocity=[3,2,4]
    case(3)
      velocity=[4,2,3]
    end select
    q=[2._dp,10._dp,20._dp,30._dp,200000._dp,.3_dp,.7_dp]
    jump=0;jump(1)=.1_dp;jump(5)=200;jump(velocity)=[2._dp,3._dp,-4._dp]
    left=q-jump;right=q+jump
    expected=[-1.75_dp,2.25_dp,3._dp,-4._dp,.029625_dp,.069125_dp]
    call mol_characteristic_slope_3d(left,q,right,400._dp,axis,[.true.,.true.],slope,ok)
    call require(ok.and.maxval(abs(slope-expected))<1e-14_dp,'linear characteristic amplitudes')
    do side=-1,1,2
      target=q+.5_dp*side*jump
      call mol_reconstruct_face_3d(q,400._dp,axis,side,slope,face,ok)
      call require(ok.and.maxval(abs(face-target)/max(1._dp,abs(target)))<1e-14_dp,'linear face state')
      call require(abs(sum(face(6:))-1)<1e-15_dp,'normalized face composition')
    end do
    right=nan
    call mol_characteristic_slope_3d(left,q,right,400._dp,axis,[.true.,.false.],slope,ok)
    call require(ok.and.all(slope==0),'closed right donor not accessed')
    left=nan
    call mol_characteristic_slope_3d(left,q,right,400._dp,axis,[.false.,.false.],slope,ok)
    call require(ok.and.all(slope==0),'covered cell zero slopes')
    slope=19
    call mol_characteristic_slope_3d(left,q,right,400._dp,axis,[.true.,.false.],slope,ok)
    call require(.not.ok.and.all(slope==19),'invalid active donor rollback')
    left=q-jump;right=left
    call mol_characteristic_slope_3d(left,q,right,400._dp,axis,[.true.,.true.],slope,ok)
    call require(ok.and.all(slope==0),'extremum limited to zero')
    slope=0;slope(1:2)=-1000000;face=23
    call mol_reconstruct_face_3d(q,400._dp,axis,1,slope,face,ok)
    call require(.not.ok.and.all(face==23),'nonpositive face rollback')
  end do
  print *, 'PASS MOL reconstruction checks:',checks
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL ',label
      error stop 1
    end if
    checks=checks+1
  end subroutine
end program
