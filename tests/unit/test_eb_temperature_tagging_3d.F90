program test_eb_temperature_tagging_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d, &
    coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_two_level_3d_mod, only: eb_reactive_patch_3d,tag_reactive_eb_temperature_3d
  implicit none
  type(eb_geometry_3d) :: coarse,full
  type(eb_reactive_patch_3d) :: patches(2),empty(0)
  real(dp) :: t(4,4,4),nan
  integer :: lower(3,2),no_lower(3,0),p,ratio
  integer, allocatable :: lo(:,:),hi(:,:),saved_lo(:,:),saved_hi(:,:)
  logical :: ok
  nan=ieee_value(0._dp,ieee_quiet_nan)
  do ratio=2,3
    call build_axis_plane_eb_geometry_3d(4*ratio,4*ratio,4*ratio,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
      'z',.31_dp,full,ok)
    call require(ok,'full geometry')
    call coarsen_eb_geometry_3d(full,ratio,coarse,ok)
    call require(ok,'coarse geometry')
    lower(:,1)=[1,1,2];lower(:,2)=[3,3,3]
    do p=1,2
      call extract_eb_geometry_3d(full,(lower(:,p)-1)*ratio+1,(lower(:,p)+1)*ratio,patches(p)%geometry,ok)
      call require(ok,'patch geometry')
      if (allocated(patches(p)%temperature)) deallocate(patches(p)%temperature)
      allocate(patches(p)%temperature(2*ratio,2*ratio,2*ratio))
      patches(p)%temperature=300._dp
    end do
    t=300._dp
    t(:,:,1)=nan ! Solid cells do not vote or invalidate the fluid state.
    t(1:2,1:2,2:3)=nan ! Covered coarse values may be stale, including invalid.
    t(3:4,3:4,3:4)=2000._dp ! Covered hot parents must not cause refinement.
    patches(1)%temperature(1,1,ratio)=1000._dp ! Inclusive threshold in a fine fluid cell.
    call tag_reactive_eb_temperature_3d(coarse,t,lower,ratio,patches,1000._dp,[2,2,2],0,8,lo,hi,ok)
    call require(ok,'fine hot cell under stale parent')
    call require(size(lo,2)==1,'only fine hot cell votes')
    call require(all(lo(:,1)==[1,1,1]).and.all(hi(:,1)==[2,2,2]),'expected tile')
    saved_lo=lo;saved_hi=hi
    call tag_reactive_eb_temperature_3d(coarse,t,lower(:,2:1:-1),ratio,patches(2:1:-1), &
      1000._dp,[2,2,2],0,8,lo,hi,ok)
    call require(ok.and.all(lo==saved_lo).and.all(hi==saved_hi),'patch-order independence')
    patches(1)%temperature=300._dp
    call tag_reactive_eb_temperature_3d(coarse,t,lower,ratio,patches,1000._dp,[2,2,2],0,8,lo,hi,ok)
    call require(ok.and.size(lo,2)==0.and.size(hi,2)==0,'no hot exposed cells yields no patches')
    t=300._dp;t(:,:,1)=nan;t(4,1,4)=1100._dp
    call tag_reactive_eb_temperature_3d(coarse,t,no_lower,ratio,empty,1000._dp,[1,1,1],1,64,lo,hi,ok)
    call require(ok.and.size(lo,2)==8,'clipped corner buffer')
    call require(all(minval(lo,dim=2)==[3,1,3]).and.all(maxval(hi,dim=2)==[4,2,4]),'buffer bounds')
    saved_lo=lo;saved_hi=hi
    call tag_reactive_eb_temperature_3d(coarse,t,no_lower,ratio,empty,1000._dp,[1,1,1],1,7,lo,hi,ok)
    call require(.not.ok.and.all(lo==saved_lo).and.all(hi==saved_hi),'patch-capacity rollback')
    t(2,2,2)=nan
    call tag_reactive_eb_temperature_3d(coarse,t,no_lower,ratio,empty,1000._dp,[1,1,1],1,64,lo,hi,ok)
    call require(.not.ok.and.all(lo==saved_lo).and.all(hi==saved_hi),'invalid exposed temperature rollback')
    t(2,2,2)=300._dp
    call tag_reactive_eb_temperature_3d(coarse,t,no_lower,ratio,empty,1000._dp,[1,1,1],huge(1),64,lo,hi,ok)
    call require(ok.and.size(lo,2)==48,'large buffer clips safely and excludes all-solid tiles')
    call tag_reactive_eb_temperature_3d(coarse,t,no_lower,ratio,empty,1000._dp,[3,3,3],0,8,lo,hi,ok)
    call require(ok.and.size(lo,2)==1,'partial edge tile')
    call require(all(lo(:,1)==[4,1,4]).and.all(hi(:,1)==[4,3,4]),'partial tile extent')
  end do
  print *, 'test_eb_temperature_tagging_3d: PASS'
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine require
end program test_eb_temperature_tagging_3d
