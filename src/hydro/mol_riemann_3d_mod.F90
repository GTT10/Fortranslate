module mol_riemann_3d_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_pelec_primitive_flux_x
  use mol_reconstruction_3d_mod, only: mol_characteristic_slope_3d,mol_reconstruct_face_3d
  implicit none
  private
  public :: mol_riemann_flux_3d
contains
  ! Four cell primitives along axis, with cells 0/1 adjacent to the face.
  ! Sound speeds belong to those cell centers, not the reconstructed states.
  ! Return intensive conservative flux in global momentum order. Caller owns
  ! open-face selection, area multiplication and EB centroid interpolation.
  subroutine mol_riemann_flux_3d(species,q,sound,axis,connected,order,flux,ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: q(:,-1:),sound(2)
    integer, intent(in) :: axis,order
    logical, intent(in) :: connected(2,2)
    real(dp), intent(inout) :: flux(:)
    logical, intent(out) :: ok
    real(dp) :: slopes(size(species)+4,2),left(size(q,1)),right(size(q,1))
    real(dp) :: rotated_left(size(q,1)),rotated_right(size(q,1)),candidate(size(flux)),local(size(flux))
    integer :: velocity(3),cell
    logical :: valid
    ok=.false.
    if(size(q,1)/=reactive_nprim(size(species)).or.size(q,2)/=4) return
    if(size(flux)/=reactive_nvar(size(species)).or.size(species)<1) return
    if(axis<1.or.axis>3.or.(order/=1.and.order/=2)) return
    if(.not.all(ieee_is_finite(sound)).or.any(sound<=0)) return
    slopes=0
    if(order==2) then
      do cell=0,1
        call mol_characteristic_slope_3d(q(:,cell-1),q(:,cell),q(:,cell+1), &
          sound(cell+1),axis,connected(:,cell+1),slopes(:,cell+1),valid)
        if(.not.valid) return
      end do
    end if
    call mol_reconstruct_face_3d(q(:,0),sound(1),axis,1,slopes(:,1),left,valid)
    if(.not.valid) return
    call mol_reconstruct_face_3d(q(:,1),sound(2),axis,-1,slopes(:,2),right,valid)
    if(.not.valid) return
    select case(axis)
    case(1)
      velocity=[2,3,4]
    case(2)
      velocity=[3,2,4]
    case(3)
      velocity=[4,2,3]
    end select
    rotated_left=left;rotated_right=right
    rotated_left(2:4)=left(velocity);rotated_right(2:4)=right(velocity)
    call reactive_pelec_primitive_flux_x(species,rotated_left,rotated_right,.5_dp*(sound(1)+sound(2)), &
      local,valid)
    if(.not.valid.or..not.all(ieee_is_finite(local))) return
    candidate=local;candidate(velocity)=local(2:4)
    flux=candidate;ok=.true.
  end subroutine
end module mol_riemann_3d_mod
