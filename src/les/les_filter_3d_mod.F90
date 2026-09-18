module les_filter_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  implicit none
  private
  public :: les_filter_weights,filter_les_field_3d
contains
  ! Standard PelePhysics box and truncated, normalized Gaussian filters.
  subroutine les_filter_weights(model,filter_ratio,weights,ok)
    character(len=*), intent(in) :: model
    integer, intent(in) :: filter_ratio
    real(dp), allocatable, intent(out) :: weights(:)
    logical, intent(out) :: ok
    integer :: radius,i
    ok=.false.
    if (filter_ratio<1.or.filter_ratio>64) return
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
