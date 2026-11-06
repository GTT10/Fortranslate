program test_spray_eb_drift
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use spray_parcel_mod, only: spray_parcel
  use spray_eb_drift_mod, only: drift_spray_level_set
  use rotating_frame_3d_mod, only: rotating_frame_3d
  implicit none
  type(spray_parcel) :: p(2),saved(2)
  type(rotating_frame_3d) :: frame
  real(dp) :: phi(0:2,0:2,0:2),lo(3)=0,hi(3)=2,change(5),normal(3),v(3),expected(3)
  real(dp) :: small_phi(0:8,0:8,0:8)
  real(dp) :: exterior(5),energy(2)=10,coordinate_energy
  character(len=16) :: faces(6)
  integer :: i,j,k,hits
  logical :: ok,contact(2)
  ! A dimensional wall contact can round just below zero in the next
  ! tetrahedral interpolation, despite an exactly representable wall position.
  do j=0,8
    small_phi(:,j,:)=real(j,dp)*.0005_dp-.00125_dp
  end do
  p(1)=spray_parcel(1_int64,[.0012_dp,.001250001_dp,.0015_dp], &
    [2.9e-6_dp,-100._dp,8.74e-6_dp],1e-9_dp,350._dp,1._dp)
  p(2)=p(1);p(2)%mass=0
  call drift_spray_level_set(small_phi,[0._dp,0._dp,0._dp],[.004_dp,.004_dp,.004_dp],1e-8_dp,8,p,change,hits,ok)
  call require(ok.and.hits==1,'dimensional wall reflection remains on fluid side')
  call require(p(1)%velocity(2)==100.and.p(1)%position(2)>.00125_dp,'dimensional reflected remainder advances')
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
  ! EB reflection precedes a root-domain outflow in the remaining flight.
  do i=0,2
    phi(i,:,:)=real(i,dp)-.5_dp
  end do
  p=saved;p(1)%position=[.75_dp,1._dp,1._dp];p(1)%velocity=[-1._dp,0._dp,0._dp];p(2)%mass=0
  saved=p;faces='outflow';change=42;exterior=43;hits=17
  call drift_spray_level_set(phi,lo,hi,2._dp,8,p,change,hits,ok,faces,energy,exterior,contact)
  call require(ok.and.hits==2.and.p(1)%mass==0,'EB then exterior escape')
  call require(contact(1).and..not.contact(2),'EB contact mask')
  call require(abs(change(2)-.06_dp)<1e-14_dp.and.abs(exterior(2)+.03_dp)<1e-14_dp,'separate wall and outflow impulse')
  call require(abs(exterior(1)+.03_dp)<1e-14_dp.and.abs(exterior(5)+.315_dp)<1e-14_dp,'outflow mass and liquid energy')
  p=saved;frame%axis=3;frame%angular_speed=2._dp
  call drift_spray_level_set(phi,lo,hi,2._dp,8,p,change,hits,ok,faces,energy,exterior,contact,frame,coordinate_energy)
  call require(ok.and.abs(exterior(5)+.015_dp)<1e-14_dp,'outflow rotating energy at exit point')
  call require(abs(coordinate_energy-.3_dp)<1e-14_dp,'direct exit coordinate conversion')
  ! Same earliest-event ordering with exterior specular reflection.
  p=saved;faces(2)='slip_wall'
  call drift_spray_level_set(phi,lo,hi,2._dp,8,p,change,hits,ok,faces,energy,exterior,contact)
  call require(ok.and.hits==2.and.abs(p(1)%position(1)-1.75_dp)<1e-14_dp,'EB then exterior reflection')
  call require(p(1)%velocity(1)==-1.and.abs(change(2)+exterior(2))<1e-14_dp,'two wall impulses cancel')
  ! Multiple periodic wraps require matching geometry at the seam.
  p=saved;p(1)%position=1;p(1)%velocity=[0._dp,0._dp,5._dp];saved=p
  faces='outflow';faces(5:6)='periodic'
  call drift_spray_level_set(phi,lo,hi,1._dp,8,p,change,hits,ok,faces,energy,exterior,contact)
  call require(ok.and.hits==3.and.abs(p(1)%position(3))<1e-14_dp,'periodic multiple wraps')
  call require(all(change==0).and.all(exterior==0),'periodic inventory unchanged')
  call require(.not.any(contact),'periodic crossings are not wall contacts')
  contact=.true.
  p=saved;change=42;exterior=43;hits=17
  call drift_spray_level_set(phi,lo,hi,1._dp,2,p,change,hits,ok,faces,energy,exterior,contact)
  call require(.not.ok.and.all(p(1)%position==saved(1)%position),'combined event limit rollback')
  call require(all(change==42).and.all(exterior==43).and.hits==17,'combined diagnostics rollback')
  call require(all(contact),'contact mask rollback')
  faces(1:2)='periodic'
  call drift_spray_level_set(phi,lo,hi,1._dp,8,p,change,hits,ok,faces,energy,exterior,contact)
  call require(.not.ok.and.all(exterior==43),'nonperiodic EB seam rejected')
  p=saved;p(1)%position=[.75_dp,1._dp,1.9_dp];p(1)%velocity=[-1._dp,0._dp,2._dp]
  faces='outflow';faces(5:6)='periodic'
  call drift_spray_level_set(phi,lo,hi,.5_dp,8,p,change,hits,ok,faces,energy,exterior,contact)
  call require(ok.and.contact(1).and.p(1)%velocity(1)>0,'EB contact after periodic wrap')
  change=42;exterior=43;hits=17
  ! A later solid start rolls back an earlier outflow.
  faces='outflow';p=saved;p(1)%position(1)=1.9_dp;p(1)%velocity=[1._dp,0._dp,0._dp]
  p(2)=p(1);p(2)%id=2;p(2)%position(1)=.1_dp;saved=p;coordinate_energy=45
  call drift_spray_level_set(phi,lo,hi,1._dp,8,p,change,hits,ok,faces,energy,exterior,contact, &
    frame,coordinate_energy)
  call require(.not.ok.and.coordinate_energy==45.and.p(1)%mass==saved(1)%mass.and.all(exterior==43), &
    'late invalid start preserves escaped parcel')
  print *, 'PASS: transactional EB and exterior parcel drift'
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
