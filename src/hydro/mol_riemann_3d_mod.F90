module mol_riemann_3d_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_pelec_primitive_flux_x, &
    reconstructed_pelec_state
  use mol_reconstruction_3d_mod, only: mol_characteristic_slope_3d,mol_reconstruct_face_3d, &
    limit_mol_positive_slope_3d
  implicit none
  private
  public :: mol_riemann_flux_3d
contains
  ! Four cell primitives along axis, with cells 0/1 adjacent to the face.
  ! Sound speeds belong to those cell centers, not the reconstructed states.
  ! Return intensive conservative flux in global momentum order. Caller owns
  ! open-face selection, area multiplication and EB centroid interpolation.
  subroutine mol_riemann_flux_3d(species,q,sound,axis,connected,order,flux,ok,lax_friedrichs,rotating_potential, &
      internal_energy_flux,reconstructed_pressure,reconstructed_kinematics)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: q(:,-1:),sound(2)
    integer, intent(in) :: axis,order
    logical, intent(in) :: connected(2,2)
    real(dp), intent(inout) :: flux(:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: lax_friedrichs
    ! Cell-center centrifugal potentials for the left/right cells (J/kg).
    real(dp), optional, intent(in) :: rotating_potential(2)
    real(dp), optional, intent(inout) :: internal_energy_flux
    ! Actual left/right pressures passed to the Riemann solver; publish on success only.
    real(dp), optional, intent(inout) :: reconstructed_pressure(2)
    ! Density pair, normal velocity pair, and cell-center sound-speed pair.
    real(dp), optional, intent(inout) :: reconstructed_kinematics(6)
    real(dp) :: auxiliary_flux
    real(dp) :: contact,potential(2),upwind_potential
    real(dp) :: slopes(size(species)+4,2),left(size(q,1)),right(size(q,1))
    real(dp) :: rotated_left(size(q,1)),rotated_right(size(q,1)),candidate(size(flux)),local(size(flux))
    integer :: velocity(3),cell
    logical :: valid,use_lax
    ok=.false.;use_lax=.false.
    if(present(lax_friedrichs)) use_lax=lax_friedrichs
    potential=0
    if(present(rotating_potential)) then
      potential=rotating_potential
      if(.not.all(ieee_is_finite(potential))) return
      if(any(potential<0)) return
    end if
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
        call limit_mol_positive_slope_3d(q(:,cell),sound(cell+1),axis,slopes(:,cell+1),valid)
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
    if(use_lax) then
      call mol_lax_friedrichs_x(species,rotated_left,rotated_right,local,valid,rotating_potential,auxiliary_flux)
    else
      call reactive_pelec_primitive_flux_x(species,rotated_left,rotated_right,.5_dp*(sound(1)+sound(2)), &
        local,valid,contact_velocity=contact,internal_energy_flux=auxiliary_flux)
      if(valid.and.present(rotating_potential)) then
        ! Pinned Riemann.H uses the right potential at zero contact speed.
        upwind_potential=merge(potential(1),potential(2),contact>0)
        local(5)=local(5)-local(1)*upwind_potential
      end if
    end if
    if(.not.valid.or..not.all(ieee_is_finite(local))) return
    candidate=local;candidate(velocity)=local(2:4)
    if(.not.ieee_is_finite(auxiliary_flux)) return
    if(present(internal_energy_flux)) internal_energy_flux=auxiliary_flux
    if(present(reconstructed_pressure)) reconstructed_pressure=[left(5),right(5)]
    if(present(reconstructed_kinematics)) &
      reconstructed_kinematics=[left(1),right(1),rotated_left(2),rotated_right(2),sound]
    flux=candidate;ok=.true.
  end subroutine
  ! Pinned Riemann.H uses acoustic-impedance contact speed plus mean sound
  ! speed here, rather than the maximum of the two one-sided signal speeds.
  subroutine mol_lax_friedrichs_x(species,left,right,flux,ok,rotating_potential,internal_energy_flux)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: left(:),right(:)
    real(dp), intent(out) :: flux(:)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: rotating_potential(2)
    real(dp) :: ul(size(flux)),ur(size(flux)),fl(size(flux)),fr(size(flux))
    real(dp), intent(out) :: internal_energy_flux
    real(dp) :: cl,cr,tl,tr,wl,wr,contact,speed,el,er
    logical :: valid
    ok=.false.
    ! These are reconstructed faces, not accepted cells. Pinned Fuego
    ! evaluates their positive transient temperatures beyond fit endpoints.
    call reconstructed_pelec_state(species,left,ul,tl,cl,valid,internal_energy_density=el)
    if(.not.valid) return
    call reconstructed_pelec_state(species,right,ur,tr,cr,valid,internal_energy_density=er)
    if(.not.valid) return
    wl=max(tiny(1._dp),cl*left(1));wr=max(tiny(1._dp),cr*right(1))
    contact=((wl*left(2)+wr*right(2))+(left(5)-right(5)))/(wl+wr)
    if(abs(contact)<1e-12_dp*.5_dp*(abs(left(2))+abs(right(2)))) contact=0
    speed=abs(contact)+.5_dp*(cl+cr)
    if(.not.ieee_is_finite(speed)) return
    if(present(rotating_potential)) then
      ul(5)=ul(5)-ul(1)*rotating_potential(1)
      ur(5)=ur(5)-ur(1)*rotating_potential(2)
    end if
    fl=ul*left(2);fr=ur*right(2)
    fl(2)=fl(2)+left(5);fr(2)=fr(2)+right(5)
    fl(5)=(ul(5)+left(5))*left(2);fr(5)=(ur(5)+right(5))*right(2)
    internal_energy_flux=.5_dp*(el*left(2)+er*right(2)-speed*(er-el))
    flux=.5_dp*(fl+fr-speed*(ur-ul))
    flux(1)=sum(flux(6:))
    ok=all(ieee_is_finite(flux))
  end subroutine mol_lax_friedrichs_x
end module mol_riemann_3d_mod
