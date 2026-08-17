module reconstruction_plm_mod
  use precision_mod, only: dp
  use constants_mod, only: density_floor, pressure_floor
  use state_indices_mod, only: ncons, nprim, qrho, qp
  use state_conversion_mod, only: conserved_to_primitive, primitive_to_conserved
  use slope_limiter_mod, only: limited_slope
  implicit none
  private

  public :: reconstruct_plm_faces

contains

  subroutine reconstruct_plm_faces( &
      conserved, nx, gamma, limiter, boundary_condition, left_faces, right_faces, ok)
    integer, intent(in) :: nx
    real(dp), intent(in) :: conserved(ncons, 0:nx + 1)
    real(dp), intent(in) :: gamma
    character(len=*), intent(in) :: limiter, boundary_condition
    real(dp), intent(out) :: left_faces(ncons, 0:nx)
    real(dp), intent(out) :: right_faces(ncons, 0:nx)
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:, :), slopes(:, :)
    real(dp) :: left_primitive(nprim), right_primitive(nprim)
    real(dp) :: theta
    logical :: cell_ok, limiter_ok, left_ok, right_ok
    integer :: i, component

    allocate(primitive(nprim, 0:nx + 1))
    allocate(slopes(nprim, 0:nx + 1))
    left_faces = 0.0_dp
    right_faces = 0.0_dp
    slopes = 0.0_dp
    ok = .false.

    do i = 0, nx + 1
      call conserved_to_primitive(conserved(:, i), gamma, primitive(:, i), cell_ok)
      if (.not. cell_ok) return
    end do

    do i = 1, nx
      do component = 1, nprim
        call limited_slope( &
          primitive(component, i) - primitive(component, i - 1), &
          primitive(component, i + 1) - primitive(component, i), &
          limiter, slopes(component, i), limiter_ok)
        if (.not. limiter_ok) return
      end do

      theta = min( &
        positivity_scale(primitive(qrho, i), slopes(qrho, i), density_floor), &
        positivity_scale(primitive(qp, i), slopes(qp, i), pressure_floor))
      slopes(:, i) = theta * slopes(:, i)
    end do

    select case (trim(boundary_condition))
    case ("outflow")
      slopes(:, 0) = 0.0_dp
      slopes(:, nx + 1) = 0.0_dp
      slopes(:, 1) = 0.0_dp
      slopes(:, nx) = 0.0_dp
    case ("periodic")
      slopes(:, 0) = slopes(:, nx)
      slopes(:, nx + 1) = slopes(:, 1)
    case default
      return
    end select

    do i = 0, nx
      left_primitive = primitive(:, i) + 0.5_dp * slopes(:, i)
      right_primitive = primitive(:, i + 1) - 0.5_dp * slopes(:, i + 1)

      call primitive_to_conserved(left_primitive, gamma, left_faces(:, i), left_ok)
      if (.not. left_ok) left_faces(:, i) = conserved(:, i)

      call primitive_to_conserved(right_primitive, gamma, right_faces(:, i), right_ok)
      if (.not. right_ok) right_faces(:, i) = conserved(:, i + 1)
    end do

    ok = .true.
  end subroutine reconstruct_plm_faces

  pure real(dp) function positivity_scale(center, slope, lower_bound) result(theta)
    real(dp), intent(in) :: center, slope, lower_bound
    real(dp) :: slope_magnitude

    slope_magnitude = abs(slope)
    if (slope_magnitude <= tiny(1.0_dp)) then
      theta = 1.0_dp
    else if (center - 0.5_dp * slope_magnitude > lower_bound) then
      theta = 1.0_dp
    else
      theta = max(0.0_dp, min(1.0_dp, &
        2.0_dp * (center - lower_bound) / slope_magnitude))
    end if
  end function positivity_scale

end module reconstruction_plm_mod
