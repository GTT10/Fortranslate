module eb_polyhedron_clip_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  implicit none
  private
  public :: clip_level_set_triangle,clip_level_set_tetrahedron
contains
  pure function cross(a,b) result(c)
    real(dp), intent(in) :: a(3),b(3)
    real(dp) :: c(3)
    c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
  end function

  pure subroutine clipped_triangle(vertices,phi,polygon,count)
    real(dp), intent(in) :: vertices(3,3),phi(3)
    real(dp), intent(out) :: polygon(3,8)
    integer, intent(out) :: count
    integer :: i,j
    real(dp) :: fraction
    count=0;polygon=0._dp
    if (maxval(phi)<=0._dp) return
    do i=1,3
      j=mod(i,3)+1
      if (phi(i)>=0._dp) then
        count=count+1;polygon(:,count)=vertices(:,i)
      end if
      if ((phi(i)>0._dp.and.phi(j)<0._dp).or.(phi(i)<0._dp.and.phi(j)>0._dp)) then
        fraction=phi(i)/(phi(i)-phi(j))
        count=count+1;polygon(:,count)=vertices(:,i)+fraction*(vertices(:,j)-vertices(:,i))
      end if
    end do
  end subroutine

  pure subroutine polygon_metrics(polygon,count,area,moment)
    real(dp), intent(in) :: polygon(3,8)
    integer, intent(in) :: count
    real(dp), intent(out) :: area,moment(3)
    real(dp) :: weight
    integer :: i
    area=0._dp;moment=0._dp
    do i=2,count-1
      weight=.5_dp*norm2(cross(polygon(:,i)-polygon(:,1),polygon(:,i+1)-polygon(:,1)))
      area=area+weight
      moment=moment+weight*(polygon(:,1)+polygon(:,i)+polygon(:,i+1))/3._dp
    end do
  end subroutine

  pure subroutine clip_level_set_triangle(vertices,values,area,moment,ok)
    real(dp), intent(in) :: vertices(3,3),values(3)
    real(dp), intent(out) :: area,moment(3)
    logical, intent(out) :: ok
    real(dp) :: polygon(3,8),phi(3),scale
    integer :: count
    area=0._dp;moment=0._dp;ok=.false.
    if (.not.all(ieee_is_finite(vertices)).or..not.all(ieee_is_finite(values))) return
    scale=maxval(abs(values))
    if (scale==0._dp) return
    phi=values/scale
    call clipped_triangle(vertices,phi,polygon,count)
    call polygon_metrics(polygon,count,area,moment)
    ok=ieee_is_finite(area).and.all(ieee_is_finite(moment))
  end subroutine

  ! Moments are extensive (volume*centroid and area*centroid). Boundary normals
  ! point toward positive level set, i.e. from solid into fluid.
  pure subroutine clip_level_set_tetrahedron(vertices,values,volume,moment,area,wall_moment,normal_integral,ok)
    real(dp), intent(in) :: vertices(3,4),values(4)
    real(dp), intent(out) :: volume,moment(3),area,wall_moment(3),normal_integral(3)
    logical, intent(out) :: ok
    integer, parameter :: face(3,4)=reshape([1,2,3,1,2,4,1,3,4,2,3,4],[3,4])
    real(dp) :: polygon(3,8,5),phi(4),cap(3,8),angles(8),point(3),center(3),e(3,3)
    real(dp) :: gradient(3),unit(3),u(3),v(3),det,scale,weight,t,tolerance,angle
    integer :: counts(5),i,j,ncap,points
    volume=0._dp;moment=0._dp;area=0._dp;wall_moment=0._dp;normal_integral=0._dp;ok=.false.
    if (.not.all(ieee_is_finite(vertices)).or..not.all(ieee_is_finite(values))) return
    scale=maxval(abs(values))
    if (scale==0._dp) return
    phi=values/scale
    do i=1,3
      e(:,i)=vertices(:,i+1)-vertices(:,1)
    end do
    det=dot_product(e(:,1),cross(e(:,2),e(:,3)))
    if (.not.ieee_is_finite(det).or.abs(det)<=tiny(1._dp)) return
    if (minval(phi)>0._dp) then
      volume=abs(det)/6._dp;moment=volume*sum(vertices,dim=2)/4._dp;ok=.true.;return
    end if
    if (maxval(phi)<=0._dp) then
      ok=.true.;return
    end if
    gradient=((phi(2)-phi(1))*cross(e(:,2),e(:,3))+(phi(3)-phi(1))*cross(e(:,3),e(:,1)) &
      +(phi(4)-phi(1))*cross(e(:,1),e(:,2)))/det
    scale=norm2(gradient)
    if (.not.ieee_is_finite(scale).or.scale<=tiny(1._dp)) return
    unit=gradient/scale
    polygon=0._dp;counts=0
    do i=1,4
      call clipped_triangle(vertices(:,face(:,i)),phi(face(:,i)),polygon(:,:,i),counts(i))
    end do
    ncap=0;cap=0._dp
    tolerance=128._dp*epsilon(1._dp)*maxval(abs(e))
    do i=1,4
      if (phi(i)==0._dp) call add_cap(cap,ncap,vertices(:,i),tolerance)
      do j=i+1,4
        if (.not.((phi(i)>0._dp.and.phi(j)<0._dp).or.(phi(i)<0._dp.and.phi(j)>0._dp))) cycle
        t=phi(i)/(phi(i)-phi(j))
        call add_cap(cap,ncap,vertices(:,i)+t*(vertices(:,j)-vertices(:,i)),tolerance)
      end do
    end do
    if (ncap>=3) then
      center=sum(cap(:,:ncap),dim=2)/real(ncap,dp)
      u=cap(:,1)-center;u=u/norm2(u);v=cross(unit,u)
      do i=1,ncap
        angles(i)=atan2(dot_product(cap(:,i)-center,v),dot_product(cap(:,i)-center,u))
      end do
      do i=2,ncap
        point=cap(:,i);angle=angles(i);j=i-1
        do while(j>=1)
          if (angles(j)<=angle) exit
          cap(:,j+1)=cap(:,j);angles(j+1)=angles(j);j=j-1
        end do
        cap(:,j+1)=point;angles(j+1)=angle
      end do
      polygon(:,:,5)=cap;counts(5)=ncap
      call polygon_metrics(cap,ncap,area,wall_moment)
      normal_integral=area*unit
    end if
    ! Every boundary polygon and its fan triangles form pyramids from an
    ! interior point. Absolute tetra volumes make input orientation immaterial.
    center=0._dp;points=sum(counts)
    if (points==0) return
    do i=1,5
      center=center+sum(polygon(:,:counts(i),i),dim=2)
    end do
    center=center/real(points,dp)
    do i=1,5
      do j=2,counts(i)-1
        weight=abs(dot_product(polygon(:,1,i)-center, &
          cross(polygon(:,j,i)-center,polygon(:,j+1,i)-center)))/6._dp
        volume=volume+weight
        moment=moment+weight*(center+polygon(:,1,i)+polygon(:,j,i)+polygon(:,j+1,i))/4._dp
      end do
    end do
    ok=all(ieee_is_finite([volume,moment,area,wall_moment,normal_integral]))
  contains
    pure subroutine add_cap(cap,ncap,value,tolerance)
      real(dp), intent(inout) :: cap(3,8)
      integer, intent(inout) :: ncap
      real(dp), intent(in) :: value(3),tolerance
      integer :: peer
      do peer=1,ncap
        if (maxval(abs(cap(:,peer)-value))<=tolerance) return
      end do
      ncap=ncap+1;cap(:,ncap)=value
    end subroutine
  end subroutine
end module
