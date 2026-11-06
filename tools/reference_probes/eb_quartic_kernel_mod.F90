! Isolated ratio-two candidate kernel; not a production halo policy.
module eb_quartic_kernel_mod
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  public :: child_average
contains
  pure real(real64) function child_average(q,i,j,k) result(value)
    real(real64), intent(in) :: q(-2:2,-2:2,-2:2)
    integer, intent(in) :: i,j,k
    real(real64), parameter :: c(-2:2)=[-0.01171875_real64,0.0859375_real64, &
      0.5_real64,-0.0859375_real64,0.01171875_real64]
    real(real64) :: plane(-2:2,-2:2),line(-2:2)
    integer :: x,y
    do y=-2,2
    do x=-2,2
      plane(x,y)=2*(c(-2)*q(x,y,-2)+c(-1)*q(x,y,-1)+c(0)*q(x,y,0)+c(1)*q(x,y,1)+c(2)*q(x,y,2))
      if(k==1) plane(x,y)=2*q(x,y,0)-plane(x,y)
    enddo
    enddo
    do x=-2,2
      line(x)=2*(c(-2)*plane(x,-2)+c(-1)*plane(x,-1)+c(0)*plane(x,0)+c(1)*plane(x,1)+c(2)*plane(x,2))
      if(j==1) line(x)=2*plane(x,0)-line(x)
    enddo
    value=2*(c(-2)*line(-2)+c(-1)*line(-1)+c(0)*line(0)+c(1)*line(1)+c(2)*line(2))
    if(i==1) value=2*line(0)-value
  end function child_average
end module eb_quartic_kernel_mod
