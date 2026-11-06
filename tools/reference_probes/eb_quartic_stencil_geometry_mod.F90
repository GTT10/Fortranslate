! Local stencil preparation only. Global donor/time/periodic mapping is caller owned.
module eb_quartic_stencil_geometry_mod
  use precision_mod, only: dp
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: prepare_connected_quartic_stencil
contains
  subroutine prepare_connected_quartic_stencil(q,geometry,available,prepared,axes,ok)
    real(dp), intent(in) :: q(:,-2:,-2:,-2:)
    type(eb_geometry_3d), intent(in) :: geometry
    logical, intent(in) :: available(-2:2,-2:2,-2:2)
    real(dp), intent(inout) :: prepared(:,-2:,-2:,-2:)
    logical, intent(inout) :: axes(3)
    logical, intent(out) :: ok
    logical :: enabled(3),connected
    real(dp) :: candidate(size(q,1),-2:2,-2:2,-2:2)
    integer :: a,s,i,j,k,idx(3),lo(3),hi(3),mapped(3)
    ok=.false.
    if(size(q,1)<1.or.any(shape(q)/=shape(prepared))) return
    if(any(shape(q(1,:,:,:))/=[5,5,5])) return
    if(any([geometry%nx,geometry%ny,geometry%nz]/=5).or..not.geometry%is_valid()) return
    if(.not.fluid([0,0,0])) return
    if(.not.all(ieee_is_finite(q(:,0,0,0)))) return
    enabled=.true.
    do a=1,3
      do s=-2,2
        idx=0;idx(a)=s
        if(.not.fluid(idx)) enabled(a)=.false.
        if(s<2) then
          if(aperture(a,idx)<=0) enabled(a)=.false.
        endif
      enddo
    enddo
    lo=0;hi=0
    where(enabled) lo=-2
    where(enabled) hi=2
    connected=.true.
    do k=lo(3),hi(3)
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)
      idx=[i,j,k]
      if(.not.fluid(idx)) connected=.false.
      do a=1,3
        if(idx(a)<hi(a)) then
          if(aperture(a,idx)<=0) connected=.false.
        endif
      enddo
    enddo
    enddo
    enddo
    ! Axial visibility cannot bypass a missing retained diagonal or open face.
    if(.not.connected) enabled=.false.
    do k=-2,2
    do j=-2,2
    do i=-2,2
      mapped=0
      where(enabled) mapped=[i,j,k]
      if(.not.all(ieee_is_finite(q(:,mapped(1),mapped(2),mapped(3))))) return
      candidate(:,i,j,k)=q(:,mapped(1),mapped(2),mapped(3))
    enddo
    enddo
    enddo
    prepared=candidate;axes=enabled;ok=.true.
  contains
    logical function fluid(cell)
      integer, intent(in) :: cell(3)
      fluid=available(cell(1),cell(2),cell(3)).and. &
        geometry%volume_fraction(cell(1)+3,cell(2)+3,cell(3)+3)>0
    end function
    real(dp) function aperture(axis,cell) result(value)
      integer, intent(in) :: axis,cell(3)
      select case(axis)
      case(1);value=geometry%x_face_fraction(cell(1)+3,cell(2)+3,cell(3)+3)
      case(2);value=geometry%y_face_fraction(cell(1)+3,cell(2)+3,cell(3)+3)
      case(3);value=geometry%z_face_fraction(cell(1)+3,cell(2)+3,cell(3)+3)
      end select
    end function
  end subroutine
end module eb_quartic_stencil_geometry_mod
