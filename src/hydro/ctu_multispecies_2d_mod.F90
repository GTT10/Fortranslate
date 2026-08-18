module ctu_multispecies_2d_mod
  use precision_mod, only: dp
  use constants_mod, only: density_floor
  use state_indices_mod, only: &
    irho, ncons, nbase, nprim, qu, qv
  use state_conversion_mod, only: conserved_to_primitive
  use slope_limiter_mod, only: limited_slope
  use reconstruction_multispecies_mod, only: trace_passive_mass_fractions
  use multispecies_state_mod, only: &
    max_supported_species, species_negative_tolerance, multispecies_nvar, &
    species_component, normalize_mass_fractions, mass_fractions_from_state, &
    synchronize_multispecies_thermodynamics, multispecies_state_is_physical
  use ctu_2d_mod, only: &
    ctu_face_data_2d, compute_cfl_timestep_2d, advance_ctu_2d
  implicit none
  private

  public :: compute_cfl_timestep_multispecies_2d
  public :: advance_ctu_multispecies_2d
  public :: all_cells_multispecies_physical_2d

contains

  pure integer function periodic_index(index, extent) result(wrapped)
    integer, intent(in) :: index, extent

    wrapped = 1 + modulo(index - 1, extent)
  end function periodic_index

  subroutine compute_cfl_timestep_multispecies_2d( &
      conserved, nx, ny, nspecies, dx, dy, gamma, cfl, dt, ok)
    integer, intent(in) :: nx, ny, nspecies
    real(dp), intent(in) :: conserved(:, :, :)
    real(dp), intent(in) :: dx, dy, gamma, cfl
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok

    real(dp), allocatable :: base_state(:, :, :)
    integer :: nvar

    ok = .false.
    dt = 0.0_dp
    nvar = multispecies_nvar(nspecies)
    if (nvar == 0) return
    if (size(conserved, 1) /= nvar .or. &
        size(conserved, 2) /= nx .or. size(conserved, 3) /= ny) return

    allocate(base_state(ncons, nx, ny))
    base_state = conserved(1:ncons, :, :)
    call compute_cfl_timestep_2d( &
      base_state, nx, ny, dx, dy, gamma, cfl, dt, ok)
  end subroutine compute_cfl_timestep_multispecies_2d

  subroutine advance_ctu_multispecies_2d( &
      conserved, nx, ny, nspecies, dx, dy, dt, gamma, limiter, &
      riemann_solver, use_transverse_correction, ok, &
      minimum_transverse_theta, species_face_fallbacks)
    integer, intent(in) :: nx, ny, nspecies
    real(dp), intent(inout) :: conserved(:, :, :)
    real(dp), intent(in) :: dx, dy, dt, gamma
    character(len=*), intent(in) :: limiter, riemann_solver
    logical, intent(in) :: use_transverse_correction
    logical, intent(out) :: ok
    real(dp), intent(out), optional :: minimum_transverse_theta
    integer, intent(out), optional :: species_face_fallbacks

    type(ctu_face_data_2d) :: faces
    real(dp), allocatable :: base_state(:, :, :), new_state(:, :, :)
    real(dp), allocatable :: primitive(:, :, :)
    real(dp), allocatable :: mass_fractions(:, :, :)
    real(dp), allocatable :: slope_x(:, :, :), slope_y(:, :, :)
    real(dp), allocatable :: x_minus(:, :, :), x_plus(:, :, :)
    real(dp), allocatable :: y_minus(:, :, :), y_plus(:, :, :)
    real(dp), allocatable :: x_left_y(:, :, :), x_right_y(:, :, :)
    real(dp), allocatable :: y_lower_y(:, :, :), y_upper_y(:, :, :)
    real(dp), allocatable :: provisional_x_species_flux(:, :, :)
    real(dp), allocatable :: provisional_y_species_flux(:, :, :)
    real(dp), allocatable :: final_x_species_flux(:, :, :)
    real(dp), allocatable :: final_y_species_flux(:, :, :)

    real(dp) :: local_mass_fractions(max_supported_species)
    real(dp) :: corrected(max_supported_species)
    real(dp) :: theta_x, theta_y, minimum_theta
    logical :: cell_ok, limiter_ok, trace_ok, base_ok, correction_ok
    integer :: i, j, species, im, ip, jm, jp, nvar, fallback_count

    ok = .false.
    if (present(minimum_transverse_theta)) minimum_transverse_theta = 0.0_dp
    if (present(species_face_fallbacks)) species_face_fallbacks = 0

    nvar = multispecies_nvar(nspecies)
    if (nvar == 0 .or. nx < 4 .or. ny < 4) return
    if (size(conserved, 1) /= nvar .or. &
        size(conserved, 2) /= nx .or. size(conserved, 3) /= ny) return
    if (dx <= 0.0_dp .or. dy <= 0.0_dp .or. dt <= 0.0_dp) return
    if (.not. all_cells_multispecies_physical_2d( &
        conserved, nx, ny, nspecies, gamma)) return

    allocate(base_state(ncons, nx, ny), new_state(nvar, nx, ny))
    allocate(primitive(nprim, nx, ny))
    allocate(mass_fractions(nspecies, nx, ny))
    allocate(slope_x(nspecies, nx, ny), slope_y(nspecies, nx, ny))
    allocate(x_minus(nspecies, nx, ny), x_plus(nspecies, nx, ny))
    allocate(y_minus(nspecies, nx, ny), y_plus(nspecies, nx, ny))
    allocate(x_left_y(nspecies, nx, ny), x_right_y(nspecies, nx, ny))
    allocate(y_lower_y(nspecies, nx, ny), y_upper_y(nspecies, nx, ny))
    allocate(provisional_x_species_flux(nspecies, nx, ny))
    allocate(provisional_y_species_flux(nspecies, nx, ny))
    allocate(final_x_species_flux(nspecies, nx, ny))
    allocate(final_y_species_flux(nspecies, nx, ny))

    base_state = conserved(1:ncons, :, :)
    call advance_ctu_2d( &
      base_state, nx, ny, dx, dy, dt, gamma, limiter, riemann_solver, &
      use_transverse_correction, base_ok, minimum_theta, face_data=faces)
    if (.not. base_ok) return

    slope_x = 0.0_dp
    slope_y = 0.0_dp
    do j = 1, ny
      do i = 1, nx
        call conserved_to_primitive( &
          conserved(1:ncons, i, j), gamma, primitive(:, i, j), cell_ok)
        if (.not. cell_ok) return
        call mass_fractions_from_state( &
          conserved(:, i, j), nspecies, local_mass_fractions, cell_ok)
        if (.not. cell_ok) return
        mass_fractions(:, i, j) = local_mass_fractions(1:nspecies)
      end do
    end do

    do j = 1, ny
      jm = periodic_index(j - 1, ny)
      jp = periodic_index(j + 1, ny)
      do i = 1, nx
        im = periodic_index(i - 1, nx)
        ip = periodic_index(i + 1, nx)

        do species = 1, nspecies
          call limited_slope( &
            mass_fractions(species, i, j) - &
              mass_fractions(species, im, j), &
            mass_fractions(species, ip, j) - &
              mass_fractions(species, i, j), &
            limiter, slope_x(species, i, j), limiter_ok)
          if (.not. limiter_ok) return

          call limited_slope( &
            mass_fractions(species, i, j) - &
              mass_fractions(species, i, jm), &
            mass_fractions(species, i, jp) - &
              mass_fractions(species, i, j), &
            limiter, slope_y(species, i, j), limiter_ok)
          if (.not. limiter_ok) return
        end do

        theta_x = bounded_vector_slope_scale( &
          mass_fractions(:, i, j), slope_x(:, i, j))
        theta_y = bounded_vector_slope_scale( &
          mass_fractions(:, i, j), slope_y(:, i, j))
        slope_x(:, i, j) = theta_x * slope_x(:, i, j)
        slope_y(:, i, j) = theta_y * slope_y(:, i, j)
      end do
    end do

    do j = 1, ny
      do i = 1, nx
        call trace_passive_mass_fractions( &
          mass_fractions(:, i, j), slope_x(:, i, j), &
          primitive(qu, i, j), dt / dx, x_minus(:, i, j), &
          x_plus(:, i, j), trace_ok)
        if (.not. trace_ok) then
          x_minus(:, i, j) = mass_fractions(:, i, j)
          x_plus(:, i, j) = mass_fractions(:, i, j)
        end if
        call normalize_face_mass_fractions( &
          x_minus(:, i, j), mass_fractions(:, i, j), nspecies)
        call normalize_face_mass_fractions( &
          x_plus(:, i, j), mass_fractions(:, i, j), nspecies)

        call trace_passive_mass_fractions( &
          mass_fractions(:, i, j), slope_y(:, i, j), &
          primitive(qv, i, j), dt / dy, y_minus(:, i, j), &
          y_plus(:, i, j), trace_ok)
        if (.not. trace_ok) then
          y_minus(:, i, j) = mass_fractions(:, i, j)
          y_plus(:, i, j) = mass_fractions(:, i, j)
        end if
        call normalize_face_mass_fractions( &
          y_minus(:, i, j), mass_fractions(:, i, j), nspecies)
        call normalize_face_mass_fractions( &
          y_plus(:, i, j), mass_fractions(:, i, j), nspecies)
      end do
    end do

    do j = 1, ny
      jp = periodic_index(j + 1, ny)
      do i = 1, nx
        ip = periodic_index(i + 1, nx)
        call species_flux_from_mass_flux( &
          faces%provisional_x_flux(irho, i, j), x_plus(:, i, j), &
          x_minus(:, ip, j), provisional_x_species_flux(:, i, j))
        call species_flux_from_mass_flux( &
          faces%provisional_y_flux(irho, i, j), y_plus(:, i, j), &
          y_minus(:, i, jp), provisional_y_species_flux(:, i, j))
      end do
    end do

    fallback_count = 0
    do j = 1, ny
      jm = periodic_index(j - 1, ny)
      jp = periodic_index(j + 1, ny)
      do i = 1, nx
        im = periodic_index(i - 1, nx)
        ip = periodic_index(i + 1, nx)

        if (use_transverse_correction) then
          call correct_species_face( &
            faces%x_left_base(irho, i, j), x_plus(:, i, j), &
            provisional_y_species_flux(:, i, j), &
            provisional_y_species_flux(:, i, jm), 0.5_dp * dt / dy, &
            faces%theta_x_left(i, j), faces%x_left_final(irho, i, j), &
            corrected(1:nspecies), correction_ok)
          if (.not. correction_ok) fallback_count = fallback_count + 1
          if (correction_ok) then
            x_left_y(:, i, j) = corrected(1:nspecies)
          else
            x_left_y(:, i, j) = x_plus(:, i, j)
          end if

          call correct_species_face( &
            faces%x_right_base(irho, i, j), x_minus(:, ip, j), &
            provisional_y_species_flux(:, ip, j), &
            provisional_y_species_flux(:, ip, jm), 0.5_dp * dt / dy, &
            faces%theta_x_right(i, j), faces%x_right_final(irho, i, j), &
            corrected(1:nspecies), correction_ok)
          if (.not. correction_ok) fallback_count = fallback_count + 1
          if (correction_ok) then
            x_right_y(:, i, j) = corrected(1:nspecies)
          else
            x_right_y(:, i, j) = x_minus(:, ip, j)
          end if

          call correct_species_face( &
            faces%y_lower_base(irho, i, j), y_plus(:, i, j), &
            provisional_x_species_flux(:, i, j), &
            provisional_x_species_flux(:, im, j), 0.5_dp * dt / dx, &
            faces%theta_y_lower(i, j), faces%y_lower_final(irho, i, j), &
            corrected(1:nspecies), correction_ok)
          if (.not. correction_ok) fallback_count = fallback_count + 1
          if (correction_ok) then
            y_lower_y(:, i, j) = corrected(1:nspecies)
          else
            y_lower_y(:, i, j) = y_plus(:, i, j)
          end if

          call correct_species_face( &
            faces%y_upper_base(irho, i, j), y_minus(:, i, jp), &
            provisional_x_species_flux(:, i, jp), &
            provisional_x_species_flux(:, im, jp), 0.5_dp * dt / dx, &
            faces%theta_y_upper(i, j), faces%y_upper_final(irho, i, j), &
            corrected(1:nspecies), correction_ok)
          if (.not. correction_ok) fallback_count = fallback_count + 1
          if (correction_ok) then
            y_upper_y(:, i, j) = corrected(1:nspecies)
          else
            y_upper_y(:, i, j) = y_minus(:, i, jp)
          end if
        else
          x_left_y(:, i, j) = x_plus(:, i, j)
          x_right_y(:, i, j) = x_minus(:, ip, j)
          y_lower_y(:, i, j) = y_plus(:, i, j)
          y_upper_y(:, i, j) = y_minus(:, i, jp)
        end if

        call species_flux_from_mass_flux( &
          faces%final_x_flux(irho, i, j), x_left_y(:, i, j), &
          x_right_y(:, i, j), final_x_species_flux(:, i, j))
        call species_flux_from_mass_flux( &
          faces%final_y_flux(irho, i, j), y_lower_y(:, i, j), &
          y_upper_y(:, i, j), final_y_species_flux(:, i, j))
      end do
    end do

    new_state = 0.0_dp
    new_state(1:ncons, :, :) = base_state
    do j = 1, ny
      jm = periodic_index(j - 1, ny)
      do i = 1, nx
        im = periodic_index(i - 1, nx)
        do species = 1, nspecies
          new_state(species_component(species), i, j) = &
            conserved(species_component(species), i, j) - &
            dt / dx * (final_x_species_flux(species, i, j) - &
              final_x_species_flux(species, im, j)) - &
            dt / dy * (final_y_species_flux(species, i, j) - &
              final_y_species_flux(species, i, jm))
        end do
        call synchronize_multispecies_thermodynamics( &
          new_state(:, i, j), gamma, cell_ok)
        if (.not. cell_ok) return
      end do
    end do

    if (.not. all_cells_multispecies_physical_2d( &
        new_state, nx, ny, nspecies, gamma)) return

    conserved = new_state
    if (present(minimum_transverse_theta)) minimum_transverse_theta = minimum_theta
    if (present(species_face_fallbacks)) species_face_fallbacks = fallback_count
    ok = .true.
  end subroutine advance_ctu_multispecies_2d

  pure subroutine species_flux_from_mass_flux( &
      mass_flux, lower_or_left, upper_or_right, species_flux)
    real(dp), intent(in) :: mass_flux
    real(dp), intent(in) :: lower_or_left(:), upper_or_right(:)
    real(dp), intent(out) :: species_flux(:)

    real(dp) :: donor(size(species_flux)), zero_threshold
    integer :: species, nspecies

    nspecies = size(species_flux)
    species_flux = 0.0_dp
    if (size(lower_or_left) /= nspecies .or. &
        size(upper_or_right) /= nspecies) return

    zero_threshold = sqrt(epsilon(1.0_dp)) * max(1.0_dp, abs(mass_flux))
    if (mass_flux > zero_threshold) then
      donor = lower_or_left
    else if (mass_flux < -zero_threshold) then
      donor = upper_or_right
    else
      donor = 0.5_dp * (lower_or_left + upper_or_right)
    end if

    if (nspecies == 1) then
      species_flux(1) = mass_flux
    else
      do species = 1, nspecies - 1
        species_flux(species) = mass_flux * donor(species)
      end do
      species_flux(nspecies) = mass_flux - sum(species_flux(1:nspecies - 1))
    end if
  end subroutine species_flux_from_mass_flux

  pure subroutine correct_species_face( &
      base_density, base_mass_fractions, flux_high, flux_low, scale, theta, &
      corrected_density, corrected_mass_fractions, ok)
    real(dp), intent(in) :: base_density, base_mass_fractions(:)
    real(dp), intent(in) :: flux_high(:), flux_low(:), scale, theta
    real(dp), intent(in) :: corrected_density
    real(dp), intent(out) :: corrected_mass_fractions(:)
    logical, intent(out) :: ok

    real(dp) :: species_density(size(base_mass_fractions))
    real(dp) :: mass_fractions(max_supported_species)
    real(dp) :: closure_scale
    integer :: nspecies

    corrected_mass_fractions = 0.0_dp
    ok = .false.
    nspecies = size(base_mass_fractions)
    if (nspecies < 1 .or. nspecies > max_supported_species) return
    if (size(flux_high) /= nspecies .or. size(flux_low) /= nspecies) return
    if (size(corrected_mass_fractions) /= nspecies) return
    if (base_density <= density_floor .or. corrected_density <= density_floor) return
    if (scale < 0.0_dp .or. theta < 0.0_dp .or. theta > 1.0_dp) return

    species_density = base_density * base_mass_fractions - &
      theta * scale * (flux_high - flux_low)
    closure_scale = max(1.0_dp, corrected_density)
    if (any(species_density < -species_negative_tolerance * closure_scale)) return
    if (abs(sum(species_density) - corrected_density) > &
        1.0e-10_dp * closure_scale) return

    mass_fractions = 0.0_dp
    mass_fractions(1:nspecies) = max(0.0_dp, species_density) / corrected_density
    call normalize_mass_fractions(mass_fractions, nspecies, ok)
    if (ok) corrected_mass_fractions = mass_fractions(1:nspecies)
  end subroutine correct_species_face

  pure subroutine normalize_face_mass_fractions( &
      face, fallback, nspecies)
    integer, intent(in) :: nspecies
    real(dp), intent(inout) :: face(:)
    real(dp), intent(in) :: fallback(:)

    real(dp) :: work(max_supported_species)
    logical :: normalization_ok

    work = 0.0_dp
    work(1:nspecies) = face(1:nspecies)
    call normalize_mass_fractions(work, nspecies, normalization_ok)
    if (normalization_ok) then
      face(1:nspecies) = work(1:nspecies)
    else
      face(1:nspecies) = fallback(1:nspecies)
    end if
  end subroutine normalize_face_mass_fractions

  pure real(dp) function bounded_vector_slope_scale( &
      center, slope) result(theta)
    real(dp), intent(in) :: center(:), slope(:)

    real(dp) :: maximum_slope
    integer :: species

    theta = 1.0_dp
    if (size(center) /= size(slope)) then
      theta = 0.0_dp
      return
    end if

    do species = 1, size(center)
      if (abs(slope(species)) <= tiny(1.0_dp)) cycle
      maximum_slope = 2.0_dp * &
        max(0.0_dp, min(center(species), 1.0_dp - center(species)))
      theta = min(theta, maximum_slope / abs(slope(species)))
    end do
    theta = max(0.0_dp, min(1.0_dp, theta))
  end function bounded_vector_slope_scale

  pure logical function all_cells_multispecies_physical_2d( &
      conserved, nx, ny, nspecies, gamma) result(all_physical)
    integer, intent(in) :: nx, ny, nspecies
    real(dp), intent(in) :: conserved(:, :, :), gamma

    integer :: i, j, nvar

    all_physical = .false.
    nvar = multispecies_nvar(nspecies)
    if (nvar == 0) return
    if (size(conserved, 1) /= nvar .or. &
        size(conserved, 2) /= nx .or. size(conserved, 3) /= ny) return

    do j = 1, ny
      do i = 1, nx
        if (.not. multispecies_state_is_physical( &
            conserved(:, i, j), gamma, nspecies)) return
      end do
    end do
    all_physical = .true.
  end function all_cells_multispecies_physical_2d

end module ctu_multispecies_2d_mod
