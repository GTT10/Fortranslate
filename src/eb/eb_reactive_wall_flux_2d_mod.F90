module eb_reactive_wall_flux_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use state_indices_mod, only: imx, imy
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_cut_cell
  implicit none
  private

  real(dp), parameter :: unit_normal_tolerance = &
    512.0_dp * epsilon(1.0_dp)

  public :: reactive_eb_slip_wall_flux_2d
  public :: reactive_eb_slip_wall_source_2d

contains

  subroutine reactive_eb_slip_wall_flux_2d( &
      species, conserved, temperature_guess, fluid_normal, flux, &
      wall_pressure, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: conserved(:), temperature_guess
    real(dp), intent(in) :: fluid_normal(2)
    real(dp), intent(out) :: flux(:), wall_pressure
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: normal_magnitude, recovered_temperature, sound_speed
    integer :: nspecies

    flux = 0.0_dp
    wall_pressure = 0.0_dp
    ok = .false.
    nspecies = size(species)
    if (size(conserved) /= reactive_nvar(nspecies) .or. &
        size(flux) /= reactive_nvar(nspecies)) return
    if (.not. ieee_is_finite(temperature_guess) .or. &
        temperature_guess <= 0.0_dp .or. &
        any(.not. ieee_is_finite(fluid_normal))) return

    normal_magnitude = sqrt(sum(fluid_normal**2))
    if (.not. ieee_is_finite(normal_magnitude) .or. &
        abs(normal_magnitude - 1.0_dp) > unit_normal_tolerance) return

    allocate(primitive(reactive_nprim(nspecies)))
    call reactive_conserved_to_primitive( &
      species, conserved, temperature_guess, primitive, &
      recovered_temperature, sound_speed, ok)
    if (.not. ok) return

    wall_pressure = primitive(5)
    ! The geometry normal points from solid to fluid.  The outward normal of
    ! the fluid control volume therefore has the opposite sign.
    flux(imx) = -wall_pressure * fluid_normal(1) / normal_magnitude
    flux(imy) = -wall_pressure * fluid_normal(2) / normal_magnitude
    ok = all(ieee_is_finite(flux)) .and. ieee_is_finite(wall_pressure)
  end subroutine reactive_eb_slip_wall_flux_2d

  subroutine reactive_eb_slip_wall_source_2d( &
      species, state, temperature, geometry, source, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), intent(out) :: source(:, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate(:, :, :), wall_flux(:)
    real(dp) :: fluid_volume, wall_pressure
    real(dp) :: fluid_normal(2)
    logical :: local_ok
    integer :: i, j, nvar

    source = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        any(shape(source) /= shape(state)) .or. &
        any(shape(temperature) /= [geometry%nx, geometry%ny])) return

    allocate(candidate(nvar, geometry%nx, geometry%ny), wall_flux(nvar))
    candidate = 0.0_dp
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) /= eb_cut_cell) cycle
        fluid_volume = geometry%volume_fraction(i, j) * &
          geometry%dx * geometry%dy
        if (fluid_volume <= 0.0_dp .or. &
            geometry%boundary_length(i, j) <= 0.0_dp) return
        fluid_normal = [geometry%boundary_normal_x(i, j), &
          geometry%boundary_normal_y(i, j)]
        call reactive_eb_slip_wall_flux_2d( &
          species, state(:, i, j), temperature(i, j), fluid_normal, &
          wall_flux, wall_pressure, local_ok)
        if (.not. local_ok) return
        candidate(:, i, j) = -geometry%boundary_length(i, j) * &
          wall_flux / fluid_volume
      end do
    end do
    if (any(.not. ieee_is_finite(candidate))) return

    source = candidate
    ok = .true.
  end subroutine reactive_eb_slip_wall_source_2d

end module eb_reactive_wall_flux_2d_mod
