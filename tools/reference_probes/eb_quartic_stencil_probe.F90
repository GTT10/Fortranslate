! Isolated candidate only: ratio-two conservative tensor interpolation from
! pinned AMReX ccquartic_interp. Not used by PeleF production halo filling.
! Source: AMReX 1546c1ccbf41c704a33790f4bf445a2fcee1372e,
! Src/AmrCore/AMReX_Interp_3D_C.H. No positivity/monotonicity claim.
program eb_quartic_stencil_probe
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use eb_quartic_kernel_mod, only: child_average
  implicit none
  real(real64) :: coarse(-2:2,-2:2,-2:2),fine(0:1,0:1,0:1)
  integer :: status,i,j,k
  do
    read(*,*,iostat=status) coarse
    if(status<0) exit
    if(status/=0.or..not.all(ieee_is_finite(coarse))) error stop 'Invalid stencil'
    do k=0,1
    do j=0,1
    do i=0,1
      fine(i,j,k)=child_average(coarse,i,j,k)
    enddo
    enddo
    enddo
    if(.not.all(ieee_is_finite(fine))) error stop 'Nonfinite interpolated child'
    write(*,'(*(es26.17e3,1x))') fine
  enddo
end program eb_quartic_stencil_probe
