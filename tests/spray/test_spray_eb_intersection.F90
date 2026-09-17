program test_spray_eb_intersection
  use precision_mod, only: dp
  use spray_eb_intersection_mod, only: intersect_spray_level_set
  implicit none
  real(dp) :: phi(0:2,0:2,0:2),lower(3)=0,upper(3)=[2._dp,4._dp,6._dp]
  real(dp) :: start(3),finish(3),fraction,position(3),normal(3),expected(3)
  integer :: i,j,k,cell(3)
  logical :: hit,ok
  do k=0,2
  do j=0,2
  do i=0,2
    phi(i,j,k)=real(i,dp)+real(j,dp)+real(k,dp)-2.5_dp
  end do
  end do
  end do
  start=[1.8_dp,3.6_dp,5.4_dp];finish=[.2_dp,.4_dp,.6_dp]
  call run()
  expected=[1._dp,.5_dp,1._dp/3];expected=expected/norm2(expected)
  call require(ok.and.hit,'oblique entry')
  call require(abs(fraction-(5.4_dp-2.5_dp)/4.8_dp)<1e-14_dp,'analytic hit time')
  call require(maxval(abs(normal-expected))<1e-14_dp,'anisotropic solid-to-fluid normal')
  start=position;finish=start+expected*.1_dp
  call run();call require(ok.and..not.hit,'surface outward travel')
  ! Both endpoints in fluid, with a solid slab between them.
  phi(0,:,:)=1;phi(1,:,:)=-1;phi(2,:,:)=1
  start=[.1_dp,1._dp,2._dp];finish=[1.9_dp,1._dp,2._dp]
  call run()
  call require(ok.and.hit.and.abs(position(1)-.5_dp)<1e-14_dp,'first crossing through solid slab')
  call require(all(normal==[-1._dp,0._dp,0._dp]),'slab entry normal')
  start=[.5_dp,1._dp,2._dp];finish=[.7_dp,1._dp,2._dp]
  call run();call require(ok.and.hit.and.fraction==0,'surface inward travel')
  start=[1._dp,1._dp,2._dp];fraction=123;position=456;normal=789;cell=42;hit=.true.
  call run()
  call require(.not.ok.and.fraction==123.and.all(position==456).and.all(normal==789).and.all(cell==42).and.hit, &
    'solid start rejects without publishing')
  phi=1;start=upper;finish=lower
  call run();call require(ok.and..not.hit.and.fraction==1.and.all(position==finish),'closed upper endpoint')
  finish=start
  call run();call require(ok.and..not.hit,'stationary fluid point')
  finish(1)=upper(1)+1
  call run();call require(.not.ok,'outside domain rejects')
  print *, 'PASS: tetrahedral level-set particle segment intersection'
contains
  subroutine run()
    call intersect_spray_level_set(phi,lower,upper,start,finish,hit,fraction,position,normal,cell,ok)
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
