program test_eb_boundary_gradient_3d
  use precision_mod, only: dp
  use eb_boundary_gradient_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  type(eb_boundary_gradient_3d) :: stencil,saved
  real(dp) :: n(3),b(3),h(3),v(0:2,0:2,0:2),r(3),gradient(3),flux,expected
  logical :: covered(0:2,0:2,0:2),ok
  integer :: axis,bits,i,j,k,count
  h=.02_dp;b=[-.31_dp,.12_dp,-.23_dp];gradient=[2._dp,-3._dp,5._dp];count=0
  do axis=1,3
    do bits=0,7
      n=[.2_dp,.3_dp,.4_dp];n(axis)=.8_dp
      do i=1,3
        if(btest(bits,i-1)) n(i)=-n(i)
      end do
      n=n/norm2(n)
      call build_eb_quadratic_gradient_3d(n,b,.7_dp,h,stencil,ok)
      if(.not.ok) error stop 'valid stencil rejected'
      do i=0,2
        do j=0,2
          do k=0,2
            r=h*(real(stencil%offset+[i,j,k],dp)-b)
            v(i,j,k)=400+dot_product(gradient,r)+sum(r*r)+r(1)*r(3)
          end do
        end do
      end do
      covered=.false.;flux=-99
      call apply_eb_boundary_flux_3d(stencil,v,covered,400._dp,3._dp,flux,ok)
      expected=3*.7_dp*h(1)**2*dot_product(gradient,n)
      if(.not.ok.or.abs(flux-expected)>2e-14_dp) error stop 'quadratic normal derivative'
      count=count+1
    end do
  end do
  saved=stencil
  call build_eb_quadratic_gradient_3d(n,b,.7_dp,[.02_dp,.03_dp,.02_dp],stencil,ok)
  if(ok.or.any(stencil%weights/=saved%weights)) error stop 'anisotropic rollback'
  call build_eb_quadratic_gradient_3d([1._dp,0._dp,0._dp],[0._dp,0._dp,0._dp],1._dp,h,stencil,ok)
  if(.not.ok) error stop 'axis aligned construction'
  covered=abs(stencil%weights)<=1e-14_dp;v=400
  where(covered) v=ieee_value(1._dp,ieee_quiet_nan)
  call apply_eb_boundary_flux_3d(stencil,v,covered,400._dp,3._dp,flux,ok)
  if(.not.ok.or.abs(flux)>2e-14_dp) error stop 'inactive covered donors'
  covered=.true.;flux=123
  call apply_eb_boundary_flux_3d(stencil,v,covered,400._dp,3._dp,flux,ok)
  if(ok.or.flux/=123) error stop 'covered donor rollback'
  print *, 'PASS: quadratic fields in',count,'orientations; donor checks and rollback'
end program
