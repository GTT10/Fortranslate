! Experimental scalar tensor candidate, not production halo interpolation.
module eb_bounded_tensor_mod
  use, intrinsic :: iso_fortran_env, only: real64
  use eb_extremum_slope_mod, only: bounded_quartic_slope,continuous_quartic_slope
  implicit none
  private
  public :: bounded_child_average
contains
  pure real(real64) function bounded_child_average(q,i,j,k,continuous) result(value)
    real(real64), intent(in) :: q(-2:2,-2:2,-2:2)
    integer, intent(in) :: i,j,k
    logical, optional, intent(in) :: continuous
    integer, parameter :: orders(3,6)=reshape([1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1],[3,6])
    real(real64) :: plane(-2:2,-2:2),line(-2:2),samples(-2:2),offset(3)
    integer :: o,x,y,z,idx(3),a(3)
    offset=.5_real64*real([i,j,k],real64)-.25_real64
    value=0
    ! Nonlinear directional limiting need not commute. Average all six
    ! orders so no coordinate direction receives a privileged ordering.
    do o=1,6
      a=orders(:,o)
      do y=-2,2
      do x=-2,2
        do z=-2,2
          idx(a(1))=x;idx(a(2))=y;idx(a(3))=z
          samples(z)=q(idx(1),idx(2),idx(3))
        enddo
        plane(x,y)=samples(0)+offset(a(3))*selected_slope(samples)
      enddo
      enddo
      do x=-2,2
        line(x)=plane(x,0)+offset(a(2))*selected_slope(plane(x,:))
      enddo
      value=value+(line(0)+offset(a(1))*selected_slope(line))/6
    enddo
  contains
    pure real(real64) function selected_slope(v) result(s)
      real(real64), intent(in) :: v(-2:2)
      if(present(continuous)) then
        if(continuous) then
          s=continuous_quartic_slope(v);return
        endif
      endif
      s=bounded_quartic_slope(v)
    end function
  end function
end module
