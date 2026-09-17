program test_spray_eb_drift
  use precision_mod, only: dp
  use spray_parcel_mod, only: spray_parcel
  use spray_eb_drift_mod, only: drift_spray_level_set
  implicit none
  type(spray_parcel) :: p(2),saved(2)
  real(dp) :: phi(0:2,0:2,0:2),lo(3)=0,hi(3)=2,change(5),normal(3),v(3),expected(3)
  integer :: i,j,k,hits
  logical :: ok
  ! Oblique plane x+y=1, fluid to the upper right. Tangential z motion persists.
  do k=0,2
  do j=0,2
  do i=0,2
    phi(i,j,k)=real(i+j,dp)-1
  end do
  end do
  end do
  p(1)%id=1;p(1)%mass=.01_dp;p(1)%multiplicity=3;p(1)%temperature=300
  p(1)%position=[.75_dp,.75_dp,.5_dp];p(1)%velocity=[-1._dp,-1._dp,.25_dp]
  p(2)=p(1);p(2)%id=2;p(2)%position=[1.5_dp,1.5_dp,.5_dp];p(2)%velocity=0
  saved=p;normal=[1._dp,1._dp,0._dp]/sqrt(2._dp);v=p(1)%velocity
  expected=v-2*dot_product(v,normal)*normal
  call drift_spray_level_set(phi,lo,hi,.5_dp,8,p,change,hits,ok)
  call require(ok.and.hits==1,'oblique reflection')
  call require(maxval(abs(p(1)%velocity-expected))<1e-14_dp,'normal reversal and tangential preservation')
  call require(maxval(abs(p(1)%position-[.75_dp,.75_dp,.625_dp]))<1e-14_dp,'remaining flight')
  call require(maxval(abs(change(2:4)-.03_dp*(expected-v)))<1e-14_dp,'represented momentum ledger')
  call require(abs(change(5))<1e-14_dp.and.change(1)==0,'elastic energy and mass')
  call require(p(1)%mass==saved(1)%mass.and.p(1)%temperature==saved(1)%temperature,'thermodynamics unchanged')
  ! A later invalid start must roll back an earlier successful reflection.
  p=saved;p(2)%position=.1_dp;saved=p;change=42;hits=17
  call drift_spray_level_set(phi,lo,hi,.5_dp,8,p,change,hits,ok)
  call require(.not.ok.and.all(change==42).and.hits==17,'late failure diagnostics rollback')
  do i=1,2
    call require(all(p(i)%position==saved(i)%position).and.all(p(i)%velocity==saved(i)%velocity),'population rollback')
  end do
  p=saved;p(2)%position=1.5_dp;saved=p
  call drift_spray_level_set(phi,lo,hi,.5_dp,0,p,change,hits,ok)
  call require(.not.ok.and.all(p(1)%position==saved(1)%position).and.all(change==42),'hit budget rollback')
  ! Thin fluid cavity .75<x<1.25: two reflections within a root-contained trial.
  phi(0,:,:)=-3;phi(1,:,:)=1;phi(2,:,:)=-3
  p=saved;p(1)%position=[1._dp,1._dp,1._dp];p(1)%velocity=[1._dp,0._dp,0._dp]
  p(2)%mass=0
  call drift_spray_level_set(phi,lo,hi,.9_dp,8,p,change,hits,ok)
  call require(ok.and.hits==2,'multiple contacts')
  call require(abs(p(1)%position(1)-.9_dp)<1e-14_dp.and.p(1)%velocity(1)==1,'cavity final trajectory')
  call require(maxval(abs(change))<1e-14_dp,'opposite wall impulses cancel')
  call drift_spray_level_set(phi,lo,hi,0._dp,0,p,change,hits,ok)
  call require(ok.and.hits==0.and.all(change==0),'zero interval')
  print *, 'PASS: transactional specular EB parcel drift'
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
