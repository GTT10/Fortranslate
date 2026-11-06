program eb_bounded_tensor_probe
  use, intrinsic :: iso_fortran_env, only: real64
  use eb_bounded_tensor_mod, only: bounded_child_average
  implicit none
  real(real64) :: q(-2:2,-2:2,-2:2),out(8)
  integer :: ios,i,j,k
  do
    read(*,*,iostat=ios) q
    if(ios<0) exit
    if(ios/=0) stop 2
    do k=0,1
    do j=0,1
    do i=0,1
      out(1+i+2*j+4*k)=bounded_child_average(q,i,j,k)
    enddo
    enddo
    enddo
    write(*,'(8(es25.16e3,1x))') out
  enddo
end program
