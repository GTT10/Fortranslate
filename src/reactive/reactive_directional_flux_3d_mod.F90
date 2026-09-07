module reactive_directional_flux_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use state_indices_mod, only: imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar, reactive_riemann_flux_x
  implicit none
  private

  public :: reactive_riemann_flux_z

contains

  pure subroutine rotate_reactive_conserved_z_to_x(input, output)
    real(dp), intent(in) :: input(:)
    real(dp), intent(out) :: output(:)

    output = input
    if (size(input) /= size(output) .or. size(input) < iet) return
    output(imx) = input(imz)
    output(imy) = input(imy)
    output(imz) = input(imx)
  end subroutine rotate_reactive_conserved_z_to_x

  pure subroutine rotate_reactive_flux_x_to_z(input, output)
    real(dp), intent(in) :: input(:)
    real(dp), intent(out) :: output(:)

    call rotate_reactive_conserved_z_to_x(input, output)
  end subroutine rotate_reactive_flux_x_to_z

  subroutine reactive_riemann_flux_z( &
      species, lower_state, upper_state, lower_temperature_guess, &
      upper_temperature_guess, solver, flux, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: lower_state(:), upper_state(:)
    real(dp), intent(in) :: lower_temperature_guess, upper_temperature_guess
    character(len=*), intent(in) :: solver
    real(dp), intent(out) :: flux(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: lower_rotated(:), upper_rotated(:), flux_rotated(:)
    integer :: nvar

    flux = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (size(lower_state) /= nvar .or. size(upper_state) /= nvar .or. &
        size(flux) /= nvar) return

    allocate(lower_rotated(nvar), upper_rotated(nvar), flux_rotated(nvar))
    call rotate_reactive_conserved_z_to_x(lower_state, lower_rotated)
    call rotate_reactive_conserved_z_to_x(upper_state, upper_rotated)
    call reactive_riemann_flux_x( &
      species, lower_rotated, upper_rotated, lower_temperature_guess, &
      upper_temperature_guess, solver, flux_rotated, ok)
    if (.not. ok) return
    call rotate_reactive_flux_x_to_z(flux_rotated, flux)
    ok = all(ieee_is_finite(flux))
  end subroutine reactive_riemann_flux_z

end module reactive_directional_flux_3d_mod
