program eb_extremum_slope_probe
  use, intrinsic :: iso_fortran_env, only: real64
  use eb_extremum_slope_mod, only: extremum_slope => bounded_quartic_slope
  implicit none
  real(real64) :: q(-2:2),s
  integer :: ios
  do
    read(*,*,iostat=ios) q
    if(ios<0) exit
    if(ios/=0) stop 2
    s=extremum_slope(q)
    write(*,'(3(es25.16e3,1x))') s,q(0)-s/4,q(0)+s/4
  enddo
end program
