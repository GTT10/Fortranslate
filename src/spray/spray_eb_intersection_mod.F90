module spray_eb_intersection_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  private
  public :: intersect_spray_level_set
contains
  ! Match build_level_set_eb_geometry_3d's six tetrahedra per Cartesian cell.
  ! Positive phi is fluid; return the first fluid-to-solid crossing of a segment.
  ! Geometry must be the same nodal field used to construct the gas EB geometry.
  subroutine intersect_spray_level_set(phi,lower,upper,start,finish,hit,fraction,position,normal,cell,ok)
    real(dp), intent(in) :: phi(0:,0:,0:),lower(3),upper(3),start(3),finish(3)
    logical, intent(inout) :: hit
    real(dp), intent(inout) :: fraction,position(3),normal(3)
    integer, intent(inout) :: cell(3)
    logical, intent(out) :: ok
    integer, parameter :: perms(3,6)=reshape([1,2,3,1,3,2,2,1,3,2,3,1,3,1,2,3,2,1],[3,6])
    integer :: dims(3),first(3),last(3),idx(3),bits(3,4),i,j,k,t,a,b,c,q,found_cell(3)
    real(dp) :: h(3),u(3),du(3),w(4),dw(4),v(4),lo,hi,crossing,best,level,slope,gradient(3),found_normal(3)
    logical :: found,inside,checked_start
    ok=.false.;dims=shape(phi)-1
    if (any(dims<1).or.any(dims>1000000)) return
    if (product(int(dims,int64))>1000000_int64) return
    if (.not.all(ieee_is_finite(phi)).or..not.all(ieee_is_finite([lower,upper,start,finish]))) return
    if (maxval(abs(phi))>1e100_dp.or.maxval(abs([lower,upper]))>1e9_dp) return
    if (any(upper<=lower).or.any(start<lower).or.any(start>upper).or.any(finish<lower).or.any(finish>upper)) return
    h=(upper-lower)/real(dims,dp)
    if (any(h<1e-12_dp)) return
    first=max(1,min(dims,floor((min(start,finish)-lower)/h)+1))
    last=max(1,min(dims,floor((max(start,finish)-lower)/h)+1))
    du=(finish-start)/h;found=.false.;checked_start=.false.;best=1;found_normal=0;found_cell=0
    do k=first(3),last(3)
    do j=first(2),last(2)
    do i=first(1),last(1)
      idx=[i,j,k];u=(start-lower)/h-real(idx-1,dp)
      do t=1,6
        a=perms(1,t);b=perms(2,t);c=perms(3,t)
        w=[1-u(a),u(a)-u(b),u(b)-u(c),u(c)]
        dw=[-du(a),du(a)-du(b),du(b)-du(c),du(c)]
        lo=0;hi=1;inside=.true.
        do q=1,4
          if (dw(q)>0) then
            lo=max(lo,-w(q)/dw(q))
          else if (dw(q)<0) then
            hi=min(hi,-w(q)/dw(q))
          else if (w(q)<0) then
            inside=.false.;exit
          end if
        end do
        if (.not.inside.or.lo>hi) cycle
        bits=0;bits(a,2)=1;bits(:,3)=bits(:,2);bits(b,3)=1;bits(:,4)=1
        do q=1,4
          v(q)=phi(i-1+bits(1,q),j-1+bits(2,q),k-1+bits(3,q))
        end do
        level=sum(w*v);slope=sum(dw*v)
        if (lo==0) then
          checked_start=.true.
          if (level<0) return
        end if
        if (slope>=0) cycle
        crossing=-level/slope
        if (crossing<lo.or.crossing>hi) cycle
        if (found.and.crossing>=best) cycle
        gradient=0
        gradient(a)=(v(2)-v(1))/h(a)
        gradient(b)=(v(3)-v(2))/h(b)
        gradient(c)=(v(4)-v(3))/h(c)
        if (norm2(gradient)<=tiny(1._dp)) return
        found=.true.;best=crossing;found_normal=gradient/norm2(gradient);found_cell=idx
      end do
    end do
    end do
    end do
    if (.not.checked_start) return
    hit=found;fraction=best;position=start+best*(finish-start);normal=found_normal;cell=found_cell;ok=.true.
  end subroutine
end module
