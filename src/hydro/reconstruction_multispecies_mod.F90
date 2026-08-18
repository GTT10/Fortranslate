module reconstruction_multispecies_mod
  use precision_mod, only: dp
  use state_indices_mod, only: ncons, nprim, qu
  use state_conversion_mod, only: conserved_to_primitive
  use multispecies_state_mod, only: &
    max_supported_species, multispecies_nvar, mass_fractions_from_state, &
    normalize_mass_fractions, multispecies_state_from_base
  use slope_limiter_mod, only: limited_slope
  use reconstruction_plm_mod, only: reconstruct_plm_faces
  use reconstruction_pelec_plm_mod, only: &
    reconstruct_pelec_plm_faces, pelec_limited_slope, &
    pelec_flattening_coefficient
  implicit none
  private

  public :: reconstruct_multispecies_faces
  public :: trace_passive_mass_fractions

contains

  subroutine reconstruct_multispecies_faces( &
      conserved, nx, nspecies, gamma, reconstruction, limiter, &
      boundary_condition, left_faces, right_faces, ok, dtdx, &
      plm_order, use_flattening)
    integer, intent(in) :: nx, nspecies
    real(dp), intent(in) :: conserved(:, 0:), gamma
    character(len=*), intent(in) :: reconstruction, limiter
    character(len=*), intent(in) :: boundary_condition
    real(dp), intent(out) :: left_faces(:, 0:), right_faces(:, 0:)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: dtdx
    integer, intent(in), optional :: plm_order
    logical, intent(in), optional :: use_flattening

    real(dp), allocatable :: base_state(:, :), base_left(:, :), base_right(:, :)
    real(dp), allocatable :: primitive(:, :), mass_fractions(:, :)
    real(dp), allocatable :: slopes(:, :), cell_left(:, :), cell_right(:, :)
    real(dp) :: center_mass_fractions(max_supported_species)
    real(dp) :: face_mass_fractions(max_supported_species)
    real(dp) :: flat, theta, local_dtdx
    integer :: order, nvar, i, species
    logical :: flattening_enabled, cell_ok, base_ok, slope_ok, face_ok

    ok = .false.
    nvar = multispecies_nvar(nspecies)
    if (nvar == 0 .or. nx < 4) return
    if (size(conserved, 1) /= nvar .or. ubound(conserved, 2) < nx + 1) return
    if (size(left_faces, 1) /= nvar .or. ubound(left_faces, 2) < nx) return
    if (size(right_faces, 1) /= nvar .or. ubound(right_faces, 2) < nx) return

    order = 2
    flattening_enabled = .false.
    local_dtdx = 0.0_dp
    if (present(plm_order)) order = plm_order
    if (present(use_flattening)) flattening_enabled = use_flattening
    if (present(dtdx)) local_dtdx = dtdx
    if (order /= 2 .and. order /= 4) return
    if (local_dtdx < 0.0_dp) return
    if (trim(reconstruction) == "pelec_plm" .and. .not. present(dtdx)) return

    allocate(base_state(ncons, 0:nx + 1))
    allocate(base_left(ncons, 0:nx), base_right(ncons, 0:nx))
    allocate(primitive(nprim, -2:nx + 3))
    allocate(mass_fractions(nspecies, -2:nx + 3))
    allocate(slopes(nspecies, -2:nx + 3))
    allocate(cell_left(nspecies, 0:nx + 1))
    allocate(cell_right(nspecies, 0:nx + 1))

    base_state = conserved(1:ncons, 0:nx + 1)
    primitive = 0.0_dp
    mass_fractions = 0.0_dp
    slopes = 0.0_dp
    left_faces = 0.0_dp
    right_faces = 0.0_dp

    select case (trim(reconstruction))
    case ("pcm")
      do i = 0, nx
        base_left(:, i) = base_state(:, i)
        base_right(:, i) = base_state(:, i + 1)
      end do
    case ("plm")
      call reconstruct_plm_faces( &
        base_state, nx, gamma, limiter, boundary_condition, &
        base_left, base_right, base_ok)
      if (.not. base_ok) return
    case ("pelec_plm")
      call reconstruct_pelec_plm_faces( &
        base_state, nx, gamma, limiter, boundary_condition, local_dtdx, &
        base_left, base_right, base_ok, slope_order=order, &
        use_flattening=flattening_enabled)
      if (.not. base_ok) return
    case default
      return
    end select

    do i = 1, nx
      call conserved_to_primitive(base_state(:, i), gamma, primitive(:, i), cell_ok)
      if (.not. cell_ok) return
      call mass_fractions_from_state( &
        conserved(:, i), nspecies, center_mass_fractions, cell_ok)
      if (.not. cell_ok) return
      mass_fractions(:, i) = center_mass_fractions(1:nspecies)
    end do

    call fill_multispecies_ghosts( &
      primitive, mass_fractions, nx, nspecies, boundary_condition, cell_ok)
    if (.not. cell_ok) return

    if (trim(reconstruction) /= "pcm") then
      do i = 1, nx
        flat = 1.0_dp
        if (flattening_enabled .and. trim(reconstruction) == "pelec_plm") then
          flat = pelec_flattening_coefficient( &
            primitive, nx, i, boundary_condition)
        end if

        do species = 1, nspecies
          if (trim(reconstruction) == "pelec_plm" .and. order == 4) then
            call pelec_limited_slope( &
              mass_fractions(species, i - 2), &
              mass_fractions(species, i - 1), &
              mass_fractions(species, i), &
              mass_fractions(species, i + 1), &
              mass_fractions(species, i + 2), &
              flat, order, slopes(species, i), slope_ok)
          else
            call limited_slope( &
              mass_fractions(species, i) - &
                mass_fractions(species, i - 1), &
              mass_fractions(species, i + 1) - &
                mass_fractions(species, i), &
              limiter, slopes(species, i), slope_ok)
            slopes(species, i) = flat * slopes(species, i)
          end if
          if (.not. slope_ok) return
        end do

        theta = 1.0_dp
        do species = 1, nspecies
          theta = min(theta, bounded_slope_scale( &
            mass_fractions(species, i), slopes(species, i)))
        end do
        slopes(:, i) = theta * slopes(:, i)
      end do
    end if

    call fill_multispecies_slope_ghosts( &
      slopes, nx, nspecies, boundary_condition, cell_ok)
    if (.not. cell_ok) return

    do i = 0, nx + 1
      select case (trim(reconstruction))
      case ("pcm")
        cell_left(:, i) = mass_fractions(:, i)
        cell_right(:, i) = mass_fractions(:, i)
      case ("plm")
        cell_left(:, i) = mass_fractions(:, i) - 0.5_dp * slopes(:, i)
        cell_right(:, i) = mass_fractions(:, i) + 0.5_dp * slopes(:, i)
      case ("pelec_plm")
        call trace_passive_mass_fractions( &
          mass_fractions(:, i), slopes(:, i), primitive(qu, i), &
          local_dtdx, cell_left(:, i), cell_right(:, i), cell_ok)
        if (.not. cell_ok) then
          cell_left(:, i) = mass_fractions(:, i)
          cell_right(:, i) = mass_fractions(:, i)
        end if
      end select

      face_mass_fractions = 0.0_dp
      face_mass_fractions(1:nspecies) = cell_left(:, i)
      call normalize_mass_fractions(face_mass_fractions, nspecies, cell_ok)
      if (.not. cell_ok) then
        face_mass_fractions(1:nspecies) = mass_fractions(:, i)
      end if
      cell_left(:, i) = face_mass_fractions(1:nspecies)

      face_mass_fractions = 0.0_dp
      face_mass_fractions(1:nspecies) = cell_right(:, i)
      call normalize_mass_fractions(face_mass_fractions, nspecies, cell_ok)
      if (.not. cell_ok) then
        face_mass_fractions(1:nspecies) = mass_fractions(:, i)
      end if
      cell_right(:, i) = face_mass_fractions(1:nspecies)
    end do

    do i = 0, nx
      face_mass_fractions = 0.0_dp
      face_mass_fractions(1:nspecies) = cell_right(:, i)
      call multispecies_state_from_base( &
        base_left(:, i), face_mass_fractions, nspecies, gamma, &
        left_faces(:, i), face_ok)
      if (.not. face_ok) left_faces(:, i) = conserved(:, i)

      face_mass_fractions = 0.0_dp
      face_mass_fractions(1:nspecies) = cell_left(:, i + 1)
      call multispecies_state_from_base( &
        base_right(:, i), face_mass_fractions, nspecies, gamma, &
        right_faces(:, i), face_ok)
      if (.not. face_ok) right_faces(:, i) = conserved(:, i + 1)
    end do

    ok = .true.
  end subroutine reconstruct_multispecies_faces

  pure subroutine trace_passive_mass_fractions( &
      center, slope, velocity, dtdx, left_state, right_state, ok)
    real(dp), intent(in) :: center(:), slope(:), velocity, dtdx
    real(dp), intent(out) :: left_state(:), right_state(:)
    logical, intent(out) :: ok

    real(dp) :: speed_minus, speed_plus

    left_state = 0.0_dp
    right_state = 0.0_dp
    ok = .false.
    if (size(center) /= size(slope)) return
    if (size(left_state) /= size(center) .or. &
        size(right_state) /= size(center)) return
    if (dtdx < 0.0_dp) return

    if (velocity > 0.0_dp) then
      speed_minus = -1.0_dp
    else
      speed_minus = velocity * dtdx
    end if
    if (velocity >= 0.0_dp) then
      speed_plus = velocity * dtdx
    else
      speed_plus = 1.0_dp
    end if

    left_state = center + 0.5_dp * (-1.0_dp - speed_minus) * slope
    right_state = center + 0.5_dp * (1.0_dp - speed_plus) * slope
    ok = .true.
  end subroutine trace_passive_mass_fractions

  pure subroutine fill_multispecies_ghosts( &
      primitive, mass_fractions, nx, nspecies, boundary_condition, ok)
    integer, intent(in) :: nx, nspecies
    real(dp), intent(inout) :: primitive(nprim, -2:nx + 3)
    real(dp), intent(inout) :: mass_fractions(nspecies, -2:nx + 3)
    character(len=*), intent(in) :: boundary_condition
    logical, intent(out) :: ok

    integer :: i, wrapped

    select case (trim(boundary_condition))
    case ("outflow")
      do i = -2, 0
        primitive(:, i) = primitive(:, 1)
        mass_fractions(:, i) = mass_fractions(:, 1)
      end do
      do i = nx + 1, nx + 3
        primitive(:, i) = primitive(:, nx)
        mass_fractions(:, i) = mass_fractions(:, nx)
      end do
      ok = .true.
    case ("periodic")
      do i = -2, 0
        wrapped = 1 + modulo(i - 1, nx)
        primitive(:, i) = primitive(:, wrapped)
        mass_fractions(:, i) = mass_fractions(:, wrapped)
      end do
      do i = nx + 1, nx + 3
        wrapped = 1 + modulo(i - 1, nx)
        primitive(:, i) = primitive(:, wrapped)
        mass_fractions(:, i) = mass_fractions(:, wrapped)
      end do
      ok = .true.
    case default
      ok = .false.
    end select
  end subroutine fill_multispecies_ghosts

  pure subroutine fill_multispecies_slope_ghosts( &
      slopes, nx, nspecies, boundary_condition, ok)
    integer, intent(in) :: nx, nspecies
    real(dp), intent(inout) :: slopes(nspecies, -2:nx + 3)
    character(len=*), intent(in) :: boundary_condition
    logical, intent(out) :: ok

    integer :: i, wrapped

    select case (trim(boundary_condition))
    case ("outflow")
      slopes(:, -2:0) = 0.0_dp
      slopes(:, nx + 1:nx + 3) = 0.0_dp
      slopes(:, 1) = 0.0_dp
      slopes(:, nx) = 0.0_dp
      ok = .true.
    case ("periodic")
      do i = -2, 0
        wrapped = 1 + modulo(i - 1, nx)
        slopes(:, i) = slopes(:, wrapped)
      end do
      do i = nx + 1, nx + 3
        wrapped = 1 + modulo(i - 1, nx)
        slopes(:, i) = slopes(:, wrapped)
      end do
      ok = .true.
    case default
      ok = .false.
    end select
  end subroutine fill_multispecies_slope_ghosts

  pure real(dp) function bounded_slope_scale(center, slope) result(theta)
    real(dp), intent(in) :: center, slope

    real(dp) :: maximum_slope

    if (abs(slope) <= tiny(1.0_dp)) then
      theta = 1.0_dp
      return
    end if
    maximum_slope = 2.0_dp * max(0.0_dp, min(center, 1.0_dp - center))
    theta = min(1.0_dp, maximum_slope / abs(slope))
    theta = max(0.0_dp, theta)
  end function bounded_slope_scale

end module reconstruction_multispecies_mod
