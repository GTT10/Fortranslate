module eb_reactive_transfer_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d
  implicit none
  private
  public :: restrict_reactive_eb_state_3d,prolong_reactive_eb_state_pcm_3d
contains
  ! Full-domain restriction onto the geometry implied by fine and ratio.
  ! Covered destination storage is retained. Publish U and T together only
  ! after all active donor and destination states pass the existing EOS checks.
  subroutine restrict_reactive_eb_state_3d(species,fine,ratio,fine_state,fine_temperature, &
      coarse_state,coarse_temperature,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: fine
    integer, intent(in) :: ratio
    real(dp), intent(in) :: fine_state(:,:,:,:),fine_temperature(:,:,:)
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:)
    logical, intent(out) :: ok
    type(eb_geometry_3d) :: coarse
    real(dp), allocatable :: candidate(:,:,:,:),temps(:,:,:),primitive(:)
    real(dp) :: weight,volume,recovered,sound_speed,guess
    integer :: i,j,k,ii,jj,kk,nvar,allocation_status,lo(3),hi(3)
    logical :: valid
    ok=.false.
    call coarsen_eb_geometry_3d(fine,ratio,coarse,valid)
    if (.not.valid) return
    nvar=reactive_nvar(size(species))
    if (.not.fields_match(fine_state,fine_temperature,fine,nvar)) return
    if (.not.fields_match(coarse_state,coarse_temperature,coarse,nvar)) return
    allocate(candidate(nvar,coarse%nx,coarse%ny,coarse%nz), &
      temps(coarse%nx,coarse%ny,coarse%nz),primitive(reactive_nprim(size(species))),stat=allocation_status)
    if (allocation_status/=0) return
    candidate=coarse_state; temps=coarse_temperature
    do k=1,coarse%nz
      do j=1,coarse%ny
        do i=1,coarse%nx
          if (coarse%volume_fraction(i,j,k)<=0._dp) cycle
          lo=([i,j,k]-1)*ratio+1; hi=[i,j,k]*ratio
          volume=sum(fine%volume_fraction(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
          if (volume<=0._dp) return
          candidate(:,i,j,k)=0._dp; guess=0._dp
          do kk=lo(3),hi(3)
            do jj=lo(2),hi(2)
              do ii=lo(1),hi(1)
                if (fine%volume_fraction(ii,jj,kk)<=0._dp) cycle
                if (.not.ieee_is_finite(fine_temperature(ii,jj,kk))) return
                if (fine_temperature(ii,jj,kk)<=0._dp) return
                call reactive_conserved_to_primitive(species,fine_state(:,ii,jj,kk),fine_temperature(ii,jj,kk), &
                  primitive,recovered,sound_speed,valid)
                if (.not.valid) return
                weight=fine%volume_fraction(ii,jj,kk)/volume
                candidate(:,i,j,k)=candidate(:,i,j,k)+weight*fine_state(:,ii,jj,kk)
                guess=guess+weight*recovered
              end do
            end do
          end do
          call reactive_conserved_to_primitive(species,candidate(:,i,j,k),guess, &
            primitive,temps(i,j,k),sound_speed,valid)
          if (.not.valid) return
        end do
      end do
    end do
    coarse_state=candidate; coarse_temperature=temps
    ok=.true.
  end subroutine restrict_reactive_eb_state_3d

  ! Piecewise constant prolongation preserves the kappa-weighted inventory
  ! because the coarse geometry is restricted from the supplied fine geometry.
  subroutine prolong_reactive_eb_state_pcm_3d(species,fine,ratio,coarse_state,coarse_temperature, &
      fine_state,fine_temperature,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: fine
    integer, intent(in) :: ratio
    real(dp), intent(in) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:)
    real(dp), intent(inout) :: fine_state(:,:,:,:),fine_temperature(:,:,:)
    logical, intent(out) :: ok
    type(eb_geometry_3d) :: coarse
    real(dp), allocatable :: candidate(:,:,:,:),temps(:,:,:),primitive(:)
    real(dp) :: recovered,sound_speed
    integer :: i,j,k,ii,jj,kk,nvar,allocation_status,lo(3),hi(3)
    logical :: valid
    ok=.false.
    call coarsen_eb_geometry_3d(fine,ratio,coarse,valid)
    if (.not.valid) return
    nvar=reactive_nvar(size(species))
    if (.not.fields_match(fine_state,fine_temperature,fine,nvar)) return
    if (.not.fields_match(coarse_state,coarse_temperature,coarse,nvar)) return
    allocate(candidate(nvar,fine%nx,fine%ny,fine%nz),temps(fine%nx,fine%ny,fine%nz), &
      primitive(reactive_nprim(size(species))),stat=allocation_status)
    if (allocation_status/=0) return
    candidate=fine_state; temps=fine_temperature
    do k=1,coarse%nz
      do j=1,coarse%ny
        do i=1,coarse%nx
          if (coarse%volume_fraction(i,j,k)<=0._dp) cycle
          if (.not.ieee_is_finite(coarse_temperature(i,j,k))) return
          if (coarse_temperature(i,j,k)<=0._dp) return
          call reactive_conserved_to_primitive(species,coarse_state(:,i,j,k),coarse_temperature(i,j,k), &
            primitive,recovered,sound_speed,valid)
          if (.not.valid) return
          lo=([i,j,k]-1)*ratio+1; hi=[i,j,k]*ratio
          do kk=lo(3),hi(3)
            do jj=lo(2),hi(2)
              do ii=lo(1),hi(1)
                if (fine%volume_fraction(ii,jj,kk)<=0._dp) cycle
                candidate(:,ii,jj,kk)=coarse_state(:,i,j,k)
                temps(ii,jj,kk)=recovered
              end do
            end do
          end do
        end do
      end do
    end do
    fine_state=candidate; fine_temperature=temps
    ok=.true.
  end subroutine prolong_reactive_eb_state_pcm_3d

  pure logical function fields_match(state,temperature,geometry,nvar) result(valid)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:)
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: nvar
    valid=nvar>0.and.all(shape(state)==[nvar,geometry%nx,geometry%ny,geometry%nz]).and. &
      all(shape(temperature)==[geometry%nx,geometry%ny,geometry%nz])
  end function fields_match
end module eb_reactive_transfer_3d_mod
