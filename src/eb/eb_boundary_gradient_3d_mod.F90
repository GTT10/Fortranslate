module eb_boundary_gradient_3d_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: eb_boundary_gradient_3d,build_eb_quadratic_gradient_3d,apply_eb_boundary_flux_3d
  public :: apply_eb_boundary_viscous_flux_3d
  type :: eb_boundary_gradient_3d
    integer :: offset(3)=0
    real(dp) :: weights(0:2,0:2,0:2)=0,wall_weight=0
  end type
contains
  ! Pinned PeleC type-0 (Johansen-Colella) stencil. Normal has the same
  ! orientation as pc_fill_sv_ebg; centroid is relative to the cell center
  ! in cell widths, and area is normalized by h**2. Weights include area/h.
  subroutine build_eb_quadratic_gradient_3d(normal,centroid,area,spacing,stencil,ok)
    real(dp), intent(in) :: normal(3),centroid(3),area,spacing(3)
    type(eb_boundary_gradient_3d), intent(inout) :: stencil
    logical, intent(out) :: ok
    type(eb_boundary_gradient_3d) :: trial
    real(dp) :: magnitudes(3),b(3),y(2),z(2),d(2),cy(0:2),cz(0:2),w(0:2,0:2,0:2)
    real(dp) :: x(2),factor
    integer :: order(3),sgn(3),base(3),shift(3),index(3),i,j,k,m,axis
    ok=.false.
    if(.not.all(ieee_is_finite([normal,centroid,area,spacing]))) return
    if(any(spacing<=0).or.area<=0.or.any(abs(centroid)>.5_dp)) return
    if(abs(norm2(normal)-1._dp)>128*epsilon(1._dp)) return
    if(maxval(abs(spacing-spacing(1)))>1e-13_dp*spacing(1)) return
    magnitudes=abs(normal)
    do axis=1,3
      order(axis)=maxloc(magnitudes,dim=1)
      magnitudes(order(axis))=-1
    end do
    sgn=int(sign(1._dp,normal(order)))
    b=centroid(order)*sgn
    base=int(sign(1._dp,normal))-1;shift=0;x=[1._dp,2._dp]
    y=b(2)+(x-b(1))*abs(normal(order(2))/normal(order(1)))
    z=b(3)+(x-b(1))*abs(normal(order(3))/normal(order(1)))
    if(any(y<0)) then
      shift(order(2))=-sgn(2);b(2)=b(2)+1
      y=b(2)+(x-b(1))*abs(normal(order(2))/normal(order(1)))
    end if
    if(any(z<0)) then
      shift(order(3))=-sgn(3);b(3)=b(3)+1
      z=b(3)+(x-b(1))*abs(normal(order(3))/normal(order(1)))
    end if
    d=sqrt((x-b(1))**2+(y-b(2))**2+(z-b(3))**2)
    w=0
    do m=1,2
      cy=[.5_dp*(y(m)-1)*(y(m)-2),-y(m)*(y(m)-2),.5_dp*y(m)*(y(m)-1)]
      cz=[.5_dp*(z(m)-1)*(z(m)-2),-z(m)*(z(m)-2),.5_dp*z(m)*(z(m)-1)]
      factor=d(3-m)/(d(m)*(d(3-m)-d(m)))
      do k=0,2
        do j=0,2
          w(m,j,k)=cy(j)*cz(k)*factor
        end do
      end do
    end do
    ! Keep the reference h**2/h ordering for the extensive stencil.
    factor=spacing(1)**2/spacing(1)*area
    do i=0,2
      do j=0,2
        do k=0,2
          index(order)=[i,j,k]*sgn-base(order)
          trial%weights(index(1),index(2),index(3))=factor*w(i,j,k)
        end do
      end do
    end do
    trial%offset=base+shift
    trial%wall_weight=factor*(-(d(1)+d(2))/(d(1)*d(2)))
    if(.not.all(ieee_is_finite(trial%weights)).or..not.ieee_is_finite(trial%wall_weight)) return
    stencil=trial;ok=.true.
  end subroutine

  ! Extensive momentum flux from the pinned normal/tangent construction.
  ! Normal stress uses 4*mu/3+bulk; tangential stress uses mu. This routine
  ! returns traction only; it does not invent an energy flux for a moving wall.
  subroutine apply_eb_boundary_viscous_flux_3d(stencil,velocity,covered,normal,wall_velocity,mu,bulk,flux,ok)
    type(eb_boundary_gradient_3d), intent(in) :: stencil
    real(dp), intent(in) :: velocity(3,0:2,0:2,0:2),normal(3),wall_velocity(3),mu,bulk
    logical, intent(in) :: covered(0:2,0:2,0:2)
    real(dp), intent(inout) :: flux(3)
    logical, intent(out) :: ok
    real(dp) :: n(3),t1(3),t2(3),alpha(3),magnitude,weights(3),basis(3,3),total(3),stress(3),candidate(3)
    integer :: order(3),axis,i,j,k
    ok=.false.
    if(.not.all(ieee_is_finite([normal,wall_velocity,mu,bulk,stencil%wall_weight]))) return
    if(mu<0.or.bulk<0.or..not.all(ieee_is_finite(stencil%weights))) return
    magnitude=norm2(normal)
    if(.not.ieee_is_finite(magnitude).or.magnitude<=tiny(1._dp)) return
    n=normal/magnitude;weights=abs(n)
    do axis=1,3
      order(axis)=maxloc(weights,dim=1);weights(order(axis))=-1
    end do
    alpha=0;alpha(order(3))=1
    t1=alpha-dot_product(n,alpha)*n;t1=t1/norm2(t1)
    t2=[n(2)*t1(3)-n(3)*t1(2),n(3)*t1(1)-n(1)*t1(3),n(1)*t1(2)-n(2)*t1(1)]
    basis(1,:)=n;basis(2,:)=t1;basis(3,:)=t2
    total=0
    do i=0,2
      do j=0,2
        do k=0,2
          if(abs(stencil%weights(i,j,k))<=1e-14_dp) cycle
          if(covered(i,j,k).or..not.all(ieee_is_finite(velocity(:,i,j,k)))) return
          total=total+stencil%weights(i,j,k)*matmul(basis,velocity(:,i,j,k))
        end do
      end do
    end do
    total=total+matmul(basis,wall_velocity)*stencil%wall_weight
    stress=[(4._dp/3._dp*mu+bulk)*total(1),mu*total(2),mu*total(3)]
    candidate=matmul(transpose(basis),stress)
    if(.not.all(ieee_is_finite(candidate))) return
    flux=candidate;ok=.true.
  end subroutine

  ! Values/covered flags are gathered at cell+offset+[0:2,0:2,0:2].
  ! An unused donor may be covered or undefined; an active donor may not.
  ! The result is extensive (W for SI temperature/conductivity inputs).
  subroutine apply_eb_boundary_flux_3d(stencil,values,covered,wall_value,coefficient,flux,ok)
    type(eb_boundary_gradient_3d), intent(in) :: stencil
    real(dp), intent(in) :: values(0:2,0:2,0:2),wall_value,coefficient
    logical, intent(in) :: covered(0:2,0:2,0:2)
    real(dp), intent(inout) :: flux
    logical, intent(out) :: ok
    real(dp) :: total,candidate
    integer :: i,j,k
    ok=.false.
    if(.not.all(ieee_is_finite(stencil%weights))) return
    if(.not.all(ieee_is_finite([stencil%wall_weight,wall_value,coefficient]))) return
    if(coefficient<0) return
    total=0
    do i=0,2
      do j=0,2
        do k=0,2
          if(abs(stencil%weights(i,j,k))<=1e-14_dp) cycle
          if(covered(i,j,k).or..not.ieee_is_finite(values(i,j,k))) return
          total=total+stencil%weights(i,j,k)*values(i,j,k)
        end do
      end do
    end do
    candidate=coefficient*(wall_value*stencil%wall_weight+total)
    if(.not.ieee_is_finite(candidate)) return
    flux=candidate;ok=.true.
  end subroutine
end module
