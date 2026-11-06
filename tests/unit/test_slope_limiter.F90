program test_slope_limiter
  use precision_mod, only: dp
  use slope_limiter_mod, only: minmod2, minmod3, limited_slope, amrex_limited_slope, amrex_extdir_limited_slope
  implicit none

  real(dp), parameter :: tolerance = 1.0e-14_dp
  real(dp) :: slope
  logical :: ok

  call assert_close(minmod2(2.0_dp, 1.0_dp), 1.0_dp, tolerance, "minmod positive")
  call assert_close(minmod2(-2.0_dp, -1.0_dp), -1.0_dp, tolerance, "minmod negative")
  call assert_close(minmod2(-1.0_dp, 1.0_dp), 0.0_dp, tolerance, "minmod extremum")
  call assert_close(minmod3(3.0_dp, 2.0_dp, 1.0_dp), 1.0_dp, tolerance, "minmod3")

  call limited_slope(1.0_dp, 3.0_dp, "minmod", slope, ok)
  call assert_true(ok, "minmod selection")
  call assert_close(slope, 1.0_dp, tolerance, "minmod slope")

  call limited_slope(1.0_dp, 3.0_dp, "mc", slope, ok)
  call assert_true(ok, "MC selection")
  call assert_close(slope, 2.0_dp, tolerance, "MC slope")

  call limited_slope(1.0_dp, -1.0_dp, "mc", slope, ok)
  call assert_true(ok, "MC extremum selection")
  call assert_close(slope, 0.0_dp, tolerance, "MC extremum")

  call limited_slope(1.0_dp, 1.0_dp, "unknown", slope, ok)
  call assert_true(.not. ok, "unknown limiter rejection")

  ! These distinguish AMReX StateRedist from PeleC characteristic slopes.
  call amrex_limited_slope(0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
    1.0_dp, 2, slope, ok)
  call assert_true(ok, "AMReX order-two selection")
  call assert_close(slope, 1.0_dp, tolerance, "AMReX affine order-two slope")
  call amrex_limited_slope(0.0_dp, 1.0_dp, 8.0_dp, 27.0_dp, 64.0_dp, &
    1.0_dp, 2, slope, ok)
  call assert_close(slope, 13.0_dp, tolerance, "AMReX curved order-two slope")
  call amrex_limited_slope(0.0_dp, 1.0_dp, 8.0_dp, 27.0_dp, 64.0_dp, &
    1.0_dp, 4, slope, ok)
  call assert_true(ok, "AMReX order-four selection")
  call assert_close(slope, 37.0_dp / 3.0_dp, tolerance, "AMReX curved order-four slope")
  call amrex_limited_slope(0.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 0.0_dp, &
    1.0_dp, 4, slope, ok)
  call assert_close(slope, 0.0_dp, tolerance, "AMReX order-four extremum")
  call amrex_limited_slope(0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
    1.0_dp, 1, slope, ok)
  call assert_true(.not. ok .and. slope == 0.0_dp, "AMReX invalid order rejection")

  call check_extdir()
  call check_scaled_slopes()

  write(*, '(a)') "test_slope_limiter: PASS"

contains

  subroutine check_scaled_slopes()
    real(dp) :: scale,base,value,q(5),scaled(5),scales(3)
    integer :: order,position,index,which
    logical :: valid
    scales=[1e-200_dp,-1e-200_dp,1e200_dp]
    q=[1._dp,2._dp,4._dp,7._dp,11._dp]
    do order=2,4,2
      do which=1,size(scales)
        scale=scales(which);scaled=scale*q
        call amrex_limited_slope(q(1),q(2),q(3),q(4),q(5),1._dp,order,base,valid)
        call assert_true(valid,'unscaled AMReX slope')
        call amrex_limited_slope(scaled(1),scaled(2),scaled(3),scaled(4),scaled(5), &
          1._dp,order,value,valid)
        call assert_true(valid,'scaled AMReX slope')
        call assert_close(value/scale,base,8*tolerance,'AMReX slope scale invariance')
        do position=1,5
          select case(position)
          case(1);index=0
          case(2);index=1
          case(3);index=5
          case(4);index=9
          case(5);index=10
          end select
          call amrex_extdir_limited_slope(q(1),q(2),q(3),q(4),q(5),order,index,0,10, &
            .true.,.true.,base,valid)
          call assert_true(valid,'unscaled boundary slope')
          call amrex_extdir_limited_slope(scaled(1),scaled(2),scaled(3),scaled(4),scaled(5), &
            order,index,0,10,.true.,.true.,value,valid)
          call assert_true(valid,'scaled boundary slope')
          call assert_close(value/scale,base,8*tolerance,'boundary slope scale invariance')
        enddo
      enddo
    enddo
    q=[1.9669864606043052e-139_dp,2.1381619906164792e-183_dp, &
      2.6355762378475875e-210_dp,4.2460545799637185e-166_dp,2.9072126599529596e-132_dp]
    call amrex_limited_slope(q(1),q(2),q(3),q(4),q(5),1._dp,4,value,valid)
    call assert_true(valid.and.value==0,'AMReX tiny extremum zero slope')
  end subroutine check_scaled_slopes

  subroutine check_extdir()
    real(dp) :: value
    logical :: valid
    integer :: order
    do order=2,4,2
      call amrex_extdir_limited_slope(0._dp,0._dp,.5_dp,1.5_dp,2.5_dp,order,0,0,10, &
        .true.,.false.,value,valid)
      call assert_true(valid,"lower face-sampled affine slope")
      call assert_close(value,1._dp,tolerance,"lower affine slope")
      call amrex_extdir_limited_slope(.5_dp,1.5_dp,2.5_dp,3._dp,3._dp,order,10,0,10, &
        .false.,.true.,value,valid)
      call assert_true(valid,"upper face-sampled affine slope")
      call assert_close(value,1._dp,tolerance,"upper affine slope")
    enddo
    call amrex_extdir_limited_slope(1._dp,1._dp,2._dp,4._dp,7._dp,2,0,0,10, &
      .true.,.false.,value,valid)
    call assert_close(value,2._dp,tolerance,"order-two face correction")
    call amrex_extdir_limited_slope(1._dp,1._dp,2._dp,4._dp,7._dp,4,0,0,10, &
      .true.,.false.,value,valid)
    call assert_close(value,1.9_dp,tolerance,"order-four lower face correction")
    call amrex_extdir_limited_slope(-7._dp,-4._dp,-2._dp,-1._dp,-1._dp,4,10,0,10, &
      .false.,.true.,value,valid)
    call assert_close(value,1.9_dp,tolerance,"order-four upper face correction")
    call amrex_extdir_limited_slope(1._dp,2._dp,4._dp,7._dp,11._dp,4,1,0,10, &
      .true.,.false.,value,valid)
    call assert_close(value,73._dp/30,tolerance,"order-four lower adjacent correction")
    call amrex_extdir_limited_slope(-11._dp,-7._dp,-4._dp,-2._dp,-1._dp,4,9,0,10, &
      .false.,.true.,value,valid)
    call assert_close(value,73._dp/30,tolerance,"order-four upper adjacent correction")
    call amrex_extdir_limited_slope(1._dp,1._dp,2._dp,4._dp,7._dp,4,0,0,10, &
      .false.,.false.,value,valid)
    call assert_close(value,19._dp/12,tolerance,"disabled face correction")
    call amrex_extdir_limited_slope(1._dp,1._dp,2._dp,4._dp,7._dp,4,0,1,10, &
      .true.,.false.,value,valid)
    call assert_true(.not.valid.and.value==0,"invalid domain rejection")
  end subroutine check_extdir

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*, '(a,1x,a)') "FAIL:", trim(label)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected, tol, label)
    real(dp), intent(in) :: actual, expected, tol
    character(len=*), intent(in) :: label
    if (abs(actual - expected) > tol) then
      write(*, '(a,1x,a,2(1x,es24.16))') "FAIL:", trim(label), actual, expected
      error stop 1
    end if
  end subroutine assert_close

end program test_slope_limiter
