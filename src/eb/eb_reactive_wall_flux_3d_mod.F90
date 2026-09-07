module eb_reactive_wall_flux_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use state_indices_mod, only: imx, imy, imz
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d, eb_cut_cell_3d
  implicit none
  private

  real(dp), parameter :: unit_normal_tolerance = &
    512.0_dp * epsilon(1.0_dp)

  public :: reactive_eb_slip_wall_flux_3d
  public :: reactive_eb_slip_wall_source_3d
  public :: reactive_eb_flux_divergence_3d

contains

  subroutine reactive_eb_slip_wall_flux_3d( &
      species, conserved, temperature_guess, fluid_normal, flux, &
      wall_pressure, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: conserved(:), temperature_guess
    real(dp), intent(in) :: fluid_normal(3)
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
    ! The geometry normal points from solid to fluid.  The fluid control
    ! volume has the opposite outward normal at the embedded boundary.
    flux(imx) = -wall_pressure * fluid_normal(1) / normal_magnitude
    flux(imy) = -wall_pressure * fluid_normal(2) / normal_magnitude
    flux(imz) = -wall_pressure * fluid_normal(3) / normal_magnitude
    ok = all(ieee_is_finite(flux)) .and. ieee_is_finite(wall_pressure)
  end subroutine reactive_eb_slip_wall_flux_3d

  subroutine reactive_eb_slip_wall_source_3d( &
      species, state, temperature, geometry, source, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(out) :: source(:, :, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate(:, :, :, :), wall_flux(:)
    real(dp) :: fluid_volume, wall_pressure
    real(dp) :: fluid_normal(3)
    logical :: local_ok
    integer :: i, j, k, nvar

    source = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        size(state, 4) /= geometry%nz .or. &
        any(shape(source) /= shape(state)) .or. &
        any(shape(temperature) /= &
          [geometry%nx, geometry%ny, geometry%nz])) return

    allocate(candidate(nvar, geometry%nx, geometry%ny, geometry%nz))
    allocate(wall_flux(nvar))
    candidate = 0.0_dp
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) /= eb_cut_cell_3d) cycle
          fluid_volume = geometry%volume_fraction(i, j, k) * &
            geometry%dx * geometry%dy * geometry%dz
          if (fluid_volume <= 0.0_dp .or. &
              geometry%boundary_area(i, j, k) <= 0.0_dp) return
          fluid_normal = [ &
            geometry%boundary_normal_x(i, j, k), &
            geometry%boundary_normal_y(i, j, k), &
            geometry%boundary_normal_z(i, j, k)]
          call reactive_eb_slip_wall_flux_3d( &
            species, state(:, i, j, k), temperature(i, j, k), &
            fluid_normal, wall_flux, wall_pressure, local_ok)
          if (.not. local_ok) return
          candidate(:, i, j, k) = &
            -geometry%boundary_area(i, j, k) * wall_flux / fluid_volume
          candidate(imx, i, j, k) = wall_pressure * &
            geometry%boundary_normal_integral_x(i, j, k) / fluid_volume
          candidate(imy, i, j, k) = wall_pressure * &
            geometry%boundary_normal_integral_y(i, j, k) / fluid_volume
          candidate(imz, i, j, k) = wall_pressure * &
            geometry%boundary_normal_integral_z(i, j, k) / fluid_volume
        end do
      end do
    end do
    if (any(.not. ieee_is_finite(candidate))) return

    source = candidate
    ok = .true.
  end subroutine reactive_eb_slip_wall_source_3d

  subroutine reactive_eb_flux_divergence_3d( &
      species, state, temperature, geometry, x_flux, y_flux, z_flux, &
      rhs, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: x_flux(:, 0:, :, :)
    real(dp), intent(in) :: y_flux(:, :, 0:, :)
    real(dp), intent(in) :: z_flux(:, :, :, 0:)
    real(dp), intent(out) :: rhs(:, :, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate(:, :, :, :)
    real(dp), allocatable :: wall_source(:, :, :, :)
    real(dp) :: fluid_volume
    logical :: local_ok
    integer :: i, j, k, nvar

    rhs = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid()) return
    if (size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        size(state, 4) /= geometry%nz .or. &
        any(shape(rhs) /= shape(state)) .or. &
        any(shape(temperature) /= &
          [geometry%nx, geometry%ny, geometry%nz]) .or. &
        size(x_flux, 1) /= nvar .or. &
        size(x_flux, 2) /= geometry%nx + 1 .or. &
        size(x_flux, 3) /= geometry%ny .or. &
        size(x_flux, 4) /= geometry%nz .or. &
        size(y_flux, 1) /= nvar .or. &
        size(y_flux, 2) /= geometry%nx .or. &
        size(y_flux, 3) /= geometry%ny + 1 .or. &
        size(y_flux, 4) /= geometry%nz .or. &
        size(z_flux, 1) /= nvar .or. &
        size(z_flux, 2) /= geometry%nx .or. &
        size(z_flux, 3) /= geometry%ny .or. &
        size(z_flux, 4) /= geometry%nz + 1) return
    if (any(.not. ieee_is_finite(x_flux)) .or. &
        any(.not. ieee_is_finite(y_flux)) .or. &
        any(.not. ieee_is_finite(z_flux))) return

    allocate(candidate(nvar, geometry%nx, geometry%ny, geometry%nz))
    allocate(wall_source(nvar, geometry%nx, geometry%ny, geometry%nz))
    call reactive_eb_slip_wall_source_3d( &
      species, state, temperature, geometry, wall_source, local_ok)
    if (.not. local_ok) return

    candidate = 0.0_dp
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          fluid_volume = geometry%volume_fraction(i, j, k) * &
            geometry%dx * geometry%dy * geometry%dz
          if (fluid_volume <= 0.0_dp) return
          candidate(:, i, j, k) = wall_source(:, i, j, k) - ( &
            geometry%dy * geometry%dz * ( &
              geometry%x_face_fraction(i, j, k) * &
                x_flux(:, i, j, k) - &
              geometry%x_face_fraction(i - 1, j, k) * &
                x_flux(:, i - 1, j, k)) + &
            geometry%dx * geometry%dz * ( &
              geometry%y_face_fraction(i, j, k) * &
                y_flux(:, i, j, k) - &
              geometry%y_face_fraction(i, j - 1, k) * &
                y_flux(:, i, j - 1, k)) + &
            geometry%dx * geometry%dy * ( &
              geometry%z_face_fraction(i, j, k) * &
                z_flux(:, i, j, k) - &
              geometry%z_face_fraction(i, j, k - 1) * &
                z_flux(:, i, j, k - 1))) / fluid_volume
        end do
      end do
    end do
    if (any(.not. ieee_is_finite(candidate))) return

    rhs = candidate
    ok = .true.
  end subroutine reactive_eb_flux_divergence_3d

end module eb_reactive_wall_flux_3d_mod
