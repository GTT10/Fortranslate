module les_filter_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  implicit none
  private
  public :: les_filter_weights,filter_les_field_3d,les_filter_code
contains
  integer function les_filter_code(model) result(code)
    character(len=*), intent(in) :: model
    code=0
    select case(model)
    case('none');code=1
    case('box');code=2
    case('gaussian');code=3
    case('box3');code=4
    case('box5');code=5
    case('gaussian3');code=6
    case('gaussian5');code=7
    case('box3opt');code=8
    case('box5opt');code=9
    case('gaussian3opt');code=10
    case('gaussian5opt');code=11
    end select
  end function

  ! Standard PelePhysics box and truncated, normalized Gaussian filters.
  recursive subroutine les_filter_weights(model,filter_ratio,weights,ok)
    character(len=*), intent(in) :: model
    integer, intent(in) :: filter_ratio
    real(dp), allocatable, intent(out) :: weights(:)
    logical, intent(out) :: ok
    integer :: radius,i
    real(dp) :: f2,f4,r1,r2
    real(dp), parameter :: box3(10)=[.079_dp,.274_dp,1.377_dp,-2.375_dp,-1._dp, &
      -.779_dp,-.680_dp,-.627_dp,-.596_dp,-.575_dp]
    real(dp), parameter :: gauss3(10)=[.0763_dp,.2527_dp,1.116_dp,-3.144_dp,-1.102_dp, &
      -.809_dp,-.696_dp,-.638_dp,-.604_dp,-.581_dp]
    real(dp), parameter :: box51(10)=[.0886_dp,.3178_dp,1.0237_dp,2.4414_dp,.2949_dp, &
      -.5276_dp,-.6708_dp,-.7003_dp,-.7077_dp,-.6996_dp]
    real(dp), parameter :: box52(10)=[-.0169_dp,-.013_dp,.0368_dp,.5559_dp,.7096_dp, &
      .4437_dp,.3302_dp,.2767_dp,.2532_dp,.2222_dp]
    real(dp), parameter :: gauss51(10)=[.0871_dp,.2596_dp,.474_dp,.1036_dp,-.4252_dp, &
      -.6134_dp,-.6679_dp,-.6836_dp,-.6873_dp,-.687_dp]
    real(dp), parameter :: gauss52(10)=[-.0175_dp,-.0021_dp,.0785_dp,.2611_dp,.3007_dp, &
      .2696_dp,.2419_dp,.2231_dp,.2103_dp,.2014_dp]
    ok=.false.
    if (filter_ratio<1.or.filter_ratio>64) return
    f2=real(filter_ratio,dp)**2;f4=f2*f2
    select case(model)
    case('none')
      allocate(weights(1));weights=1
    case('box')
      if (filter_ratio/=1.and.mod(filter_ratio,2)/=0) return
      radius=filter_ratio/2
      allocate(weights(2*radius+1));weights=1._dp/filter_ratio
      if (radius>0) then
        weights(1)=.5_dp*weights(1);weights(size(weights))=weights(1)
      end if
    case('gaussian')
      if (mod(filter_ratio,2)/=0) return
      radius=filter_ratio/2
      allocate(weights(2*radius+1))
      do i=-radius,radius
        ! The common Gaussian prefactor cancels in normalization.
        weights(i+radius+1)=exp(-6._dp*real(i,dp)**2/real(filter_ratio,dp)**2)
      end do
      weights=weights/sum(weights)
    case('box3','gaussian3')
      allocate(weights(3))
      weights=[f2/24,(12-f2)/12,f2/24]
    case('box5')
      allocate(weights(5))
      weights(1)=(3*f4-20*f2)/5760;weights(2)=(80*f2-3*f4)/1440
      weights(3)=(3*f4-100*f2+960)/960;weights(4)=weights(2);weights(5)=weights(1)
    case('gaussian5')
      allocate(weights(5))
      weights(1)=(f4-4*f2)/1152;weights(2)=(16*f2-f4)/288
      weights(3)=(f4-20*f2+192)/192;weights(4)=weights(2);weights(5)=weights(1)
    case('box3opt','gaussian3opt')
      if (filter_ratio>10) then
        if (model=='box3opt') then
          call les_filter_weights('box',filter_ratio,weights,ok)
        else
          call les_filter_weights('gaussian3',filter_ratio,weights,ok)
        end if
        return
      end if
      r1=box3(filter_ratio)
      if (model=='gaussian3opt') r1=gauss3(filter_ratio)
      allocate(weights(3));weights(1)=r1/(1+2*r1);weights(2)=1-2*weights(1);weights(3)=weights(1)
    case('box5opt','gaussian5opt')
      if (filter_ratio>10) then
        if (model=='box5opt') then
          call les_filter_weights('box',filter_ratio,weights,ok)
        else
          call les_filter_weights('gaussian5',filter_ratio,weights,ok)
        end if
        return
      end if
      r1=box51(filter_ratio);r2=box52(filter_ratio)
      if (model=='gaussian5opt') then
        r1=gauss51(filter_ratio);r2=gauss52(filter_ratio)
      end if
      allocate(weights(5));weights(1)=r2/(1+2*r1+2*r2);weights(2)=r1/r2*weights(1)
      weights(3)=1-2*weights(1)-2*weights(2);weights(4)=weights(2);weights(5)=weights(1)
    case default
      return
    end select
    ok=.true.
  end subroutine

  ! Component-first tensor-product filter. Without first_cell, wrap periodically.
  ! With first_cell, the caller supplies all halos; output shape defines the crop.
  ! Failure preserves filtered. Separable passes avoid cubic stencil work.
  subroutine filter_les_field_3d(field,weights,filtered,ok,first_cell)
    real(dp), intent(in) :: field(:,:,:,:),weights(:)
    real(dp), intent(inout) :: filtered(:,:,:,:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: first_cell(3)
    real(dp), allocatable :: pass_x(:,:,:,:),pass_y(:,:,:,:),candidate(:,:,:,:)
    integer :: dims(3),n(3),first(3),radius,i,j,k,s,index
    logical :: periodic
    ok=.false.
    if (size(field,1)<1.or.size(filtered,1)<1) return
    dims=[size(field,2),size(field,3),size(field,4)]
    n=[size(filtered,2),size(filtered,3),size(filtered,4)]
    if (size(field,1)<1.or.size(filtered,1)/=size(field,1).or.any(dims<1).or.any(n<1)) return
    if (size(weights)<1.or.size(weights)>65.or.mod(size(weights),2)/=1) return
    if (.not.all(ieee_is_finite(weights))) return
    if (maxval(abs(weights))>1e6_dp) return
    if (sum(abs(weights))>1e6_dp.or.abs(sum(weights)-1)>1e-12_dp) return
    if (.not.all(ieee_is_finite(field))) return
    if (maxval(abs(field))>1e100_dp) return
    radius=size(weights)/2;first=1;periodic=.not.present(first_cell)
    if (present(first_cell)) first=first_cell
    if (periodic) then
      if (any(n/=dims)) return
    else
      if (any(first<1+radius).or.any(first>dims-radius)) return
      if (any(n>dims-radius-first+1)) return
    end if
    allocate(pass_x(size(field,1),n(1),dims(2),dims(3)),source=0._dp)
    allocate(pass_y(size(field,1),n(1),n(2),dims(3)),source=0._dp)
    allocate(candidate(size(field,1),n(1),n(2),n(3)),source=0._dp)
    do k=1,dims(3)
      do j=1,dims(2)
        do i=1,n(1)
          do s=-radius,radius
            index=first(1)+i-1+s
            if (periodic) index=modulo(index-1,dims(1))+1
            pass_x(:,i,j,k)=pass_x(:,i,j,k)+weights(s+radius+1)*field(:,index,j,k)
          end do
        end do
      end do
    end do
    do k=1,dims(3)
      do j=1,n(2)
        do i=1,n(1)
          do s=-radius,radius
            index=first(2)+j-1+s
            if (periodic) index=modulo(index-1,dims(2))+1
            pass_y(:,i,j,k)=pass_y(:,i,j,k)+weights(s+radius+1)*pass_x(:,i,index,k)
          end do
        end do
      end do
    end do
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          do s=-radius,radius
            index=first(3)+k-1+s
            if (periodic) index=modulo(index-1,dims(3))+1
            candidate(:,i,j,k)=candidate(:,i,j,k)+weights(s+radius+1)*pass_y(:,i,j,index)
          end do
        end do
      end do
    end do
    filtered=candidate;ok=.true.
  end subroutine
end module
