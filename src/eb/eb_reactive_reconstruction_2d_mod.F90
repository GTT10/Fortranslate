module eb_reactive_reconstruction_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use constants_mod, only: density_floor, pressure_floor
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved, reactive_conserved_to_primitive, &
    reactive_riemann_flux_x, characteristic_limited_slope, &
    trace_reactive_characteristics
  use reactive_2d_mod, only: reactive_riemann_flux_y
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_covered_cell
  implicit none
  private

  public :: build_reactive_eb_face_center_fluxes_2d
  public :: interpolate_reactive_eb_face_centroid_fluxes_2d

contains

  pure subroutine rotate_primitive_y_to_x(input, output)
    real(dp), intent(in) :: input(:)
    real(dp), intent(out) :: output(:)

    output = input
    if (size(input) /= size(output) .or. size(input) < 5) return
    output(2) = input(3)
    output(3) = input(2)
  end subroutine rotate_primitive_y_to_x

  pure real(dp) function lower_scale(center, slope, lower) result(theta)
    real(dp), intent(in) :: center, slope, lower

    if (abs(slope) <= tiny(1.0_dp) .or. center - abs(slope) > lower) then
      theta = 1.0_dp
    else
      theta = max(0.0_dp, min(1.0_dp, (center - lower) / abs(slope)))
    end if
  end function lower_scale

  pure real(dp) function upper_scale(center, slope, upper) result(theta)
    real(dp), intent(in) :: center, slope, upper

    if (abs(slope) <= tiny(1.0_dp) .or. center + abs(slope) < upper) then
      theta = 1.0_dp
    else
      theta = max(0.0_dp, min(1.0_dp, (upper - center) / abs(slope)))
    end if
  end function upper_scale

  pure real(dp) function primitive_slope_scale(center, slope, nspecies) &
      result(theta)
    real(dp), intent(in) :: center(:), slope(:)
    integer, intent(in) :: nspecies

    integer :: component, k

    theta = 0.0_dp
    if (size(center) /= size(slope)) return
    theta = 1.0_dp
    theta = min(theta, lower_scale(center(1), slope(1), density_floor))
    theta = min(theta, lower_scale(center(5), slope(5), pressure_floor))
    do k = 1, nspecies
      component = reactive_mass_fraction_component(k)
      theta = min(theta, lower_scale(center(component), slope(component), &
        0.0_dp))
      theta = min(theta, upper_scale(center(component), slope(component), &
        1.0_dp))
    end do
  end function primitive_slope_scale

  pure subroutine sanitize_primitive(primitive, fallback, nspecies)
    real(dp), intent(inout) :: primitive(:)
    real(dp), intent(in) :: fallback(:)
    integer, intent(in) :: nspecies

    real(dp) :: total
    integer :: component, k

    if (size(primitive) /= size(fallback) .or. &
        primitive(1) <= density_floor .or. &
        primitive(5) <= pressure_floor .or. &
        any(.not. ieee_is_finite(primitive))) then
      primitive = fallback
      return
    end if
    total = 0.0_dp
    do k = 1, nspecies
      component = reactive_mass_fraction_component(k)
      if (primitive(component) < -1.0e-12_dp) then
        primitive = fallback
        return
      end if
      primitive(component) = max(0.0_dp, primitive(component))
      total = total + primitive(component)
    end do
    if (total <= tiny(1.0_dp)) then
      primitive = fallback
      return
    end if
    do k = 1, nspecies
      component = reactive_mass_fraction_component(k)
      primitive(component) = primitive(component) / total
    end do
  end subroutine sanitize_primitive

  subroutine primitive_face_to_state( &
      species, face_primitive, fallback, state, temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: face_primitive(:), fallback(:)
    real(dp), intent(out) :: state(:), temperature
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate(:)
    real(dp) :: sound_speed

    allocate(candidate(size(face_primitive)))
    candidate = face_primitive
    call sanitize_primitive(candidate, fallback, size(species))
    call reactive_primitive_to_conserved( &
      species, candidate, state, temperature, sound_speed, ok)
    if (ok) return
    call reactive_primitive_to_conserved( &
      species, fallback, state, temperature, sound_speed, ok)
  end subroutine primitive_face_to_state

  subroutine build_reactive_eb_face_center_fluxes_2d( &
      species, state, temperature, geometry, solver, reconstruction, &
      limiter, dt, x_flux, y_flux, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :), temperature(:, :)
    type(eb_geometry_2d), intent(in) :: geometry
    character(len=*), intent(in) :: solver, reconstruction, limiter
    real(dp), intent(in) :: dt
    real(dp), intent(out) :: x_flux(:, 0:, :), y_flux(:, :, 0:)
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate_x(:, :, :), candidate_y(:, :, :)
    real(dp), allocatable :: primitive(:, :, :), sound_speed(:, :)
    real(dp), allocatable :: slope_x(:, :, :), slope_y(:, :, :)
    real(dp), allocatable :: x_minus(:, :, :), x_plus(:, :, :)
    real(dp), allocatable :: y_minus(:, :, :), y_plus(:, :, :)
    real(dp), allocatable :: dl(:), dr(:), rotated_center(:)
    real(dp), allocatable :: rotated_dl(:), rotated_dr(:), rotated_slope(:)
    real(dp), allocatable :: rotated_minus(:), rotated_plus(:)
    real(dp), allocatable :: left_state(:), right_state(:), face_flux(:)
    real(dp) :: local_temperature, left_temperature, right_temperature
    real(dp) :: theta
    logical :: local_ok
    integer :: i, j, left_i, right_i, lower_j, upper_j
    integer :: nvar, nprimitive

    x_flux = 0.0_dp
    y_flux = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    nprimitive = reactive_nprim(size(species))
    if (nvar <= 0 .or. nprimitive <= 0 .or. .not. geometry%is_valid()) return
    if (size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        any(shape(temperature) /= [geometry%nx, geometry%ny]) .or. &
        size(x_flux, 1) /= nvar .or. &
        size(x_flux, 2) /= geometry%nx + 1 .or. &
        size(x_flux, 3) /= geometry%ny .or. &
        size(y_flux, 1) /= nvar .or. &
        size(y_flux, 2) /= geometry%nx .or. &
        size(y_flux, 3) /= geometry%ny + 1 .or. &
        len_trim(solver) == 0 .or. &
        .not. ieee_is_finite(dt) .or. dt < 0.0_dp .or. &
        any(.not. ieee_is_finite(state)) .or. &
        any(.not. ieee_is_finite(temperature))) return
    if (trim(reconstruction) /= "pcm" .and. &
        trim(reconstruction) /= "characteristic_plm") return
    if (trim(reconstruction) == "characteristic_plm" .and. &
        trim(limiter) /= "minmod" .and. trim(limiter) /= "mc") return

    allocate(candidate_x(nvar, 0:geometry%nx, geometry%ny))
    allocate(candidate_y(nvar, geometry%nx, 0:geometry%ny))
    allocate(primitive(nprimitive, geometry%nx, geometry%ny))
    allocate(sound_speed(geometry%nx, geometry%ny))
    allocate(slope_x(nprimitive, geometry%nx, geometry%ny))
    allocate(slope_y(nprimitive, geometry%nx, geometry%ny))
    allocate(x_minus(nprimitive, geometry%nx, geometry%ny))
    allocate(x_plus(nprimitive, geometry%nx, geometry%ny))
    allocate(y_minus(nprimitive, geometry%nx, geometry%ny))
    allocate(y_plus(nprimitive, geometry%nx, geometry%ny))
    allocate(dl(nprimitive), dr(nprimitive), rotated_center(nprimitive))
    allocate(rotated_dl(nprimitive), rotated_dr(nprimitive))
    allocate(rotated_slope(nprimitive), rotated_minus(nprimitive))
    allocate(rotated_plus(nprimitive), left_state(nvar), right_state(nvar))
    allocate(face_flux(nvar))
    candidate_x = 0.0_dp
    candidate_y = 0.0_dp
    primitive = 0.0_dp
    sound_speed = 0.0_dp
    slope_x = 0.0_dp
    slope_y = 0.0_dp
    x_minus = 0.0_dp
    x_plus = 0.0_dp
    y_minus = 0.0_dp
    y_plus = 0.0_dp

    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) == eb_covered_cell) cycle
        call reactive_conserved_to_primitive( &
          species, state(:, i, j), temperature(i, j), primitive(:, i, j), &
          local_temperature, sound_speed(i, j), local_ok)
        if (.not. local_ok) return
      end do
    end do

    if (trim(reconstruction) == "characteristic_plm") then
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j) == eb_covered_cell) cycle
          if (i > 1 .and. i < geometry%nx) then
            if (geometry%cell_type(i - 1, j) /= eb_covered_cell .and. &
                geometry%cell_type(i + 1, j) /= eb_covered_cell) then
              dl = primitive(:, i, j) - primitive(:, i - 1, j)
              dr = primitive(:, i + 1, j) - primitive(:, i, j)
              call characteristic_limited_slope( &
                primitive(:, i, j), dl, dr, sound_speed(i, j), limiter, &
                slope_x(:, i, j), local_ok)
              if (.not. local_ok) return
              theta = primitive_slope_scale( &
                primitive(:, i, j), slope_x(:, i, j), size(species))
              slope_x(:, i, j) = theta * slope_x(:, i, j)
            end if
          end if
          if (j > 1 .and. j < geometry%ny) then
            if (geometry%cell_type(i, j - 1) /= eb_covered_cell .and. &
                geometry%cell_type(i, j + 1) /= eb_covered_cell) then
              call rotate_primitive_y_to_x( &
                primitive(:, i, j), rotated_center)
              call rotate_primitive_y_to_x( &
                primitive(:, i, j) - primitive(:, i, j - 1), rotated_dl)
              call rotate_primitive_y_to_x( &
                primitive(:, i, j + 1) - primitive(:, i, j), rotated_dr)
              call characteristic_limited_slope( &
                rotated_center, rotated_dl, rotated_dr, sound_speed(i, j), &
                limiter, rotated_slope, local_ok)
              if (.not. local_ok) return
              theta = primitive_slope_scale( &
                rotated_center, rotated_slope, size(species))
              rotated_slope = theta * rotated_slope
              call rotate_primitive_y_to_x( &
                rotated_slope, slope_y(:, i, j))
            end if
          end if
        end do
      end do
    end if

    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) == eb_covered_cell) cycle
        if (trim(reconstruction) == "pcm") then
          x_minus(:, i, j) = primitive(:, i, j)
          x_plus(:, i, j) = primitive(:, i, j)
          y_minus(:, i, j) = primitive(:, i, j)
          y_plus(:, i, j) = primitive(:, i, j)
        else
          call trace_reactive_characteristics( &
            primitive(:, i, j), slope_x(:, i, j), sound_speed(i, j), &
            dt / geometry%dx, x_minus(:, i, j), x_plus(:, i, j), local_ok)
          if (.not. local_ok) return
          call sanitize_primitive( &
            x_minus(:, i, j), primitive(:, i, j), size(species))
          call sanitize_primitive( &
            x_plus(:, i, j), primitive(:, i, j), size(species))
          call rotate_primitive_y_to_x( &
            primitive(:, i, j), rotated_center)
          call rotate_primitive_y_to_x(slope_y(:, i, j), rotated_slope)
          call trace_reactive_characteristics( &
            rotated_center, rotated_slope, sound_speed(i, j), &
            dt / geometry%dy, rotated_minus, rotated_plus, local_ok)
          if (.not. local_ok) return
          call rotate_primitive_y_to_x(rotated_minus, y_minus(:, i, j))
          call rotate_primitive_y_to_x(rotated_plus, y_plus(:, i, j))
          call sanitize_primitive( &
            y_minus(:, i, j), primitive(:, i, j), size(species))
          call sanitize_primitive( &
            y_plus(:, i, j), primitive(:, i, j), size(species))
        end if
      end do
    end do

    do j = 1, geometry%ny
      do i = 0, geometry%nx
        if (geometry%x_face_fraction(i, j) <= 0.0_dp) cycle
        left_i = max(1, i)
        right_i = min(geometry%nx, i + 1)
        if (geometry%cell_type(left_i, j) == eb_covered_cell .or. &
            geometry%cell_type(right_i, j) == eb_covered_cell) return
        if (i == 0) then
          call primitive_face_to_state( &
            species, primitive(:, 1, j), primitive(:, 1, j), left_state, &
            left_temperature, local_ok)
          if (.not. local_ok) return
          call primitive_face_to_state( &
            species, x_minus(:, 1, j), primitive(:, 1, j), right_state, &
            right_temperature, local_ok)
        else if (i == geometry%nx) then
          call primitive_face_to_state( &
            species, x_plus(:, geometry%nx, j), &
            primitive(:, geometry%nx, j), left_state, left_temperature, &
            local_ok)
          if (.not. local_ok) return
          call primitive_face_to_state( &
            species, primitive(:, geometry%nx, j), &
            primitive(:, geometry%nx, j), right_state, right_temperature, &
            local_ok)
        else
          call primitive_face_to_state( &
            species, x_plus(:, i, j), primitive(:, i, j), left_state, &
            left_temperature, local_ok)
          if (.not. local_ok) return
          call primitive_face_to_state( &
            species, x_minus(:, i + 1, j), primitive(:, i + 1, j), &
            right_state, right_temperature, local_ok)
        end if
        if (.not. local_ok) return
        call reactive_riemann_flux_x( &
          species, left_state, right_state, left_temperature, &
          right_temperature, solver, face_flux, local_ok)
        if (.not. local_ok) return
        candidate_x(:, i, j) = face_flux
      end do
    end do

    do j = 0, geometry%ny
      do i = 1, geometry%nx
        if (geometry%y_face_fraction(i, j) <= 0.0_dp) cycle
        lower_j = max(1, j)
        upper_j = min(geometry%ny, j + 1)
        if (geometry%cell_type(i, lower_j) == eb_covered_cell .or. &
            geometry%cell_type(i, upper_j) == eb_covered_cell) return
        if (j == 0) then
          call primitive_face_to_state( &
            species, primitive(:, i, 1), primitive(:, i, 1), left_state, &
            left_temperature, local_ok)
          if (.not. local_ok) return
          call primitive_face_to_state( &
            species, y_minus(:, i, 1), primitive(:, i, 1), right_state, &
            right_temperature, local_ok)
        else if (j == geometry%ny) then
          call primitive_face_to_state( &
            species, y_plus(:, i, geometry%ny), &
            primitive(:, i, geometry%ny), left_state, left_temperature, &
            local_ok)
          if (.not. local_ok) return
          call primitive_face_to_state( &
            species, primitive(:, i, geometry%ny), &
            primitive(:, i, geometry%ny), right_state, right_temperature, &
            local_ok)
        else
          call primitive_face_to_state( &
            species, y_plus(:, i, j), primitive(:, i, j), left_state, &
            left_temperature, local_ok)
          if (.not. local_ok) return
          call primitive_face_to_state( &
            species, y_minus(:, i, j + 1), primitive(:, i, j + 1), &
            right_state, right_temperature, local_ok)
        end if
        if (.not. local_ok) return
        call reactive_riemann_flux_y( &
          species, left_state, right_state, left_temperature, &
          right_temperature, solver, face_flux, local_ok)
        if (.not. local_ok) return
        candidate_y(:, i, j) = face_flux
      end do
    end do
    if (any(.not. ieee_is_finite(candidate_x)) .or. &
        any(.not. ieee_is_finite(candidate_y))) return

    x_flux = candidate_x
    y_flux = candidate_y
    ok = .true.
  end subroutine build_reactive_eb_face_center_fluxes_2d

  subroutine interpolate_reactive_eb_face_centroid_fluxes_2d( &
      geometry, center_x_flux, center_y_flux, centroid_x_flux, &
      centroid_y_flux, ok)
    type(eb_geometry_2d), intent(in) :: geometry
    real(dp), intent(in) :: center_x_flux(:, 0:, :), center_y_flux(:, :, 0:)
    real(dp), intent(out) :: centroid_x_flux(:, 0:, :)
    real(dp), intent(out) :: centroid_y_flux(:, :, 0:)
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate_x(:, :, :), candidate_y(:, :, :)
    real(dp) :: offset, weight
    integer :: i, j, neighbor

    centroid_x_flux = 0.0_dp
    centroid_y_flux = 0.0_dp
    ok = .false.
    if (.not. geometry%is_valid()) return
    if (size(center_x_flux, 2) /= geometry%nx + 1 .or. &
        size(center_x_flux, 3) /= geometry%ny .or. &
        size(center_y_flux, 2) /= geometry%nx .or. &
        size(center_y_flux, 3) /= geometry%ny + 1 .or. &
        any(shape(centroid_x_flux) /= shape(center_x_flux)) .or. &
        any(shape(centroid_y_flux) /= shape(center_y_flux)) .or. &
        size(center_x_flux, 1) /= size(center_y_flux, 1) .or. &
        any(.not. ieee_is_finite(center_x_flux)) .or. &
        any(.not. ieee_is_finite(center_y_flux))) return

    allocate(candidate_x( &
      size(center_x_flux, 1), 0:geometry%nx, geometry%ny))
    allocate(candidate_y( &
      size(center_y_flux, 1), geometry%nx, 0:geometry%ny))
    candidate_x = center_x_flux
    candidate_y = center_y_flux
    do j = 1, geometry%ny
      do i = 0, geometry%nx
        if (geometry%x_face_fraction(i, j) <= 0.0_dp) then
          candidate_x(:, i, j) = 0.0_dp
          cycle
        end if
        offset = geometry%x_face_centroid_y(i, j)
        if (offset == 0.0_dp) cycle
        neighbor = j + merge(1, -1, offset > 0.0_dp)
        if (neighbor < 1 .or. neighbor > geometry%ny) cycle
        if (geometry%x_face_fraction(i, neighbor) <= 0.0_dp) cycle
        weight = abs(offset)
        candidate_x(:, i, j) = &
          (1.0_dp - weight) * center_x_flux(:, i, j) + &
          weight * center_x_flux(:, i, neighbor)
      end do
    end do
    do j = 0, geometry%ny
      do i = 1, geometry%nx
        if (geometry%y_face_fraction(i, j) <= 0.0_dp) then
          candidate_y(:, i, j) = 0.0_dp
          cycle
        end if
        offset = geometry%y_face_centroid_x(i, j)
        if (offset == 0.0_dp) cycle
        neighbor = i + merge(1, -1, offset > 0.0_dp)
        if (neighbor < 1 .or. neighbor > geometry%nx) cycle
        if (geometry%y_face_fraction(neighbor, j) <= 0.0_dp) cycle
        weight = abs(offset)
        candidate_y(:, i, j) = &
          (1.0_dp - weight) * center_y_flux(:, i, j) + &
          weight * center_y_flux(:, neighbor, j)
      end do
    end do
    if (any(.not. ieee_is_finite(candidate_x)) .or. &
        any(.not. ieee_is_finite(candidate_y))) return

    centroid_x_flux = candidate_x
    centroid_y_flux = candidate_y
    ok = .true.
  end subroutine interpolate_reactive_eb_face_centroid_fluxes_2d

end module eb_reactive_reconstruction_2d_mod
