program test_amr_hierarchy_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
  use precision_mod, only: dp
  use amr_reactive_3d_mod, only: reflux_physical_coarse_3d,reflux_coarse_3d
  use amr_hierarchy_3d_mod, only: &
    amr_patch_3d, initialize_amr_patch_3d, prolong_pcm_3d, &
    restrict_average_3d, average_down_3d, composite_integrals_amr_3d
  implicit none

  integer, parameter :: nvar = 4, nx = 6, ny = 5, nz = 4
  type(amr_patch_3d) :: patch
  real(dp) :: coarse(nvar, nx, ny, nz), averaged(nvar, nx, ny, nz)
  real(dp), allocatable :: fine(:, :, :, :), restricted(:, :, :, :)
  real(dp) :: composite(nvar), reference(nvar), error
  logical :: ok
  integer :: component, i, j, k

  call initialize_amr_patch_3d( &
    nx, ny, nz, 2, 4, 2, 4, 2, 3, 2, patch, ok)
  call require(ok, "patch initialization")
  call require(patch%is_strictly_interior(), "strictly interior patch")
  call require( &
    patch%fine_nx() == 6 .and. patch%fine_ny() == 6 .and. &
    patch%fine_nz() == 4, "fine extents")

  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        do component = 1, nvar
          coarse(component, i, j, k) = real( &
            1000 * component + 100 * k + 10 * j + i, dp)
        end do
      end do
    end do
  end do
  allocate(fine(nvar, patch%fine_nx(), patch%fine_ny(), patch%fine_nz()))
  allocate(restricted(nvar, 3, 3, 2))
  call prolong_pcm_3d(coarse, patch, fine, ok)
  call require(ok, "PCM prolongation")
  call restrict_average_3d(fine, patch, restricted, ok)
  call require(ok, "restriction")
  error = maxval(abs(restricted - coarse(:, 2:4, 2:4, 2:3)))
  call require(error == 0.0_dp, "PCM restriction identity")

  call composite_integrals_amr_3d( &
    coarse, fine, patch, 0.2_dp, 0.3_dp, 0.4_dp, composite, ok)
  call require(ok, "composite integration")
  do component = 1, nvar
    reference(component) = 0.2_dp * 0.3_dp * 0.4_dp * &
      sum(coarse(component, :, :, :))
  end do
  error = maxval(abs(composite - reference) / max(1.0_dp, abs(reference)))
  call require(error <= 5.0e-15_dp, "composite PCM integral")
  call composite_integrals_amr_3d( &
    coarse, fine, patch, ieee_value(0.0_dp, ieee_quiet_nan), 0.3_dp, &
    0.4_dp, composite, ok)
  call require(.not. ok .and. all(composite == 0.0_dp), &
    "non-finite composite spacing rejection")

  fine = fine + 0.125_dp
  averaged = coarse
  call average_down_3d(averaged, fine, patch, ok)
  call require(ok, "average down")
  error = maxval(abs(averaged(:, 2:4, 2:4, 2:3) - &
    (coarse(:, 2:4, 2:4, 2:3) + 0.125_dp)))
  call require(error <= 5.0e-13_dp, "covered coarse synchronization")
  call require(maxval(abs(averaged(:, 1, :, :) - coarse(:, 1, :, :))) == &
    0.0_dp, "uncovered coarse preservation")

  call initialize_amr_patch_3d( &
    nx, ny, nz, 0, 4, 2, 4, 2, 3, 2, patch, ok)
  call require(.not. ok, "invalid patch rejection")
  call check_physical_reflux()
  write(*, '(a)') "test_amr_hierarchy_3d: PASS"


contains

  subroutine check_physical_reflux()
    integer, parameter :: nc=6,nv=4
    type(amr_patch_3d) :: p
    real(dp) :: c(nv,nc,nc,nc),saved(nv,nc,nc,nc),legacy(nv,nc,nc,nc)
    real(dp) :: fx(nv,0:nc,nc,nc),fy(nv,nc,0:nc,nc),fz(nv,nc,nc,0:nc)
    real(dp), allocatable :: f(:,:,:,:),xl(:,:,:),xh(:,:,:),yl(:,:,:),yh(:,:,:),zl(:,:,:),zh(:,:,:)
    real(dp) :: spacing(3),dt,values(nv),budget(nv),expected(nv),total(nv),rate(nv),largest,old_largest
    character(len=32) :: faces(6)
    integer :: mode,ratio,lo(3),hi(3),width(3),a,v,face
    logical :: valid
    spacing=[.2_dp,.3_dp,.4_dp];dt=.01_dp; values=[1._dp,-2._dp,3._dp,-4._dp]
    do v=1,nv
      fx(v,:,:,:)=.25_dp*values(v);fy(v,:,:,:)=.5_dp*values(v);fz(v,:,:,:)=.75_dp*values(v)
    end do
    do ratio=2,3
      do mode=0,9
        lo=3;hi=4
        if (mode>=1.and.mode<=6) then
          a=(mode+1)/2
          if (modulo(mode,2)==1) lo(a)=1
          if (modulo(mode,2)==0) hi(a)=nc
        else if (mode==7) then
          lo=1;hi=nc
        else if (mode==8) then
          lo=1;hi=3
        else if (mode==9) then
          lo=4;hi=nc
        end if
        width=hi-lo+1;faces='outflow'
        call initialize_amr_patch_3d(nc,nc,nc,lo(1),hi(1),lo(2),hi(2),lo(3),hi(3),ratio,p,valid)
        call require(valid,'physical reflux patch')
        allocate(f(nv,p%fine_nx(),p%fine_ny(),p%fine_nz()))
        allocate(xl(nv,width(2),width(3)),xh(nv,width(2),width(3)))
        allocate(yl(nv,width(1),width(3)),yh(nv,width(1),width(3)))
        allocate(zl(nv,width(1),width(2)),zh(nv,width(1),width(2)))
        do v=1,nv
          xl(v,:,:)=values(v);xh(v,:,:)=2*values(v)
          yl(v,:,:)=3*values(v);yh(v,:,:)=4*values(v)
          zl(v,:,:)=5*values(v);zh(v,:,:)=6*values(v)
        end do
        ! Integrate a linear fine flux field with these six face values.
        rate=-dt*values*sum(1._dp/(real(width,dp)*spacing))
        do v=1,nv
          f(v,:,:,:)=rate(v)
        end do
        c=0;legacy=0;expected=0
        do face=1,6
          a=(face+1)/2
          if (modulo(face,2)==1.and.lo(a)==1) &
            expected=expected+dt*(face-.25_dp*a)*values*product(real(width,dp)*spacing)/(width(a)*spacing(a))
          if (modulo(face,2)==0.and.hi(a)==nc) &
            expected=expected-dt*(face-.25_dp*a)*values*product(real(width,dp)*spacing)/(width(a)*spacing(a))
        end do
        call reflux_physical_coarse_3d(p,c,fx,fy,fz,xl,xh,yl,yh,zl,zh,spacing,dt,faces,budget,largest,valid)
        call require(valid,'physical reflux accepted')
        call require(maxval(abs(budget-expected))<1e-13_dp,'physical surface transfer')
        call composite_integrals_amr_3d(c,f,p,spacing(1),spacing(2),spacing(3),total,valid)
        call require(valid.and.maxval(abs(total-budget))<1e-12_dp,'composite flux conservation')
        call require(all(c(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))==0),'covered cells untouched')
        if (mode==0) then
          call reflux_coarse_3d(p,legacy,fx(:,1:,:,:),fy(:,:,1:,:),fz(:,:,:,1:), &
            xl,xh,yl,yh,zl,zh,spacing(1),spacing(2),spacing(3),dt,old_largest,valid)
          call require(valid.and.all(c==legacy).and.largest==old_largest,'interior legacy equality')
        else if (mode==7) then
          call require(all(c==0).and.largest==0,'whole-domain patch only changes surface ledger')
        end if
        saved=c
        xl(1,1,1)=ieee_value(0._dp,ieee_quiet_nan)
        call reflux_physical_coarse_3d(p,c,fx,fy,fz,xl,xh,yl,yh,zl,zh,spacing,dt,faces,budget,largest,valid)
        call require(.not.valid.and.all(c==saved).and.all(budget==0).and.largest==0,'invalid flux rollback')
        xl(1,1,1)=values(1)
        if (mode/=0) then
          faces='periodic'
          call reflux_physical_coarse_3d(p,c,fx,fy,fz,xl,xh,yl,yh,zl,zh,spacing,dt,faces,budget,largest,valid)
          call require(.not.valid.and.all(c==saved).and.all(budget==0),'periodic seam rejected transactionally')
        end if
        deallocate(f,xl,xh,yl,yh,zl,zh)
      end do
    end do
  end subroutine check_physical_reflux


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
      write(*, '(a)') "FAIL: " // trim(message)
      error stop 1
    end if
  end subroutine require

end program test_amr_hierarchy_3d
