program test_multispecies_advection_2d
  use precision_mod, only: dp
  use state_indices_mod, only: ncons, nprim, qrho, qu, qv, qw, qp
  use mesh_2d_mod, only: uniform_cell_centers_2d
  use state_conversion_mod, only: primitive_to_conserved
  use multispecies_state_mod, only: &
    multispecies_nvar, species_component, multispecies_state_from_base, &
    mass_fractions_from_state, species_closure_error
  use ctu_multispecies_2d_mod, only: &
    compute_cfl_timestep_multispecies_2d, advance_ctu_multispecies_2d
  implicit none

  integer, parameter :: nspecies = 2
  integer, parameter :: number_of_resolutions = 3
  integer, parameter :: resolutions(number_of_resolutions) = [20, 40, 80]
  real(dp), parameter :: minimum_order = 1.8_dp
  real(dp), parameter :: conservation_tolerance = 2.0e-11_dp
  real(dp), parameter :: closure_tolerance = 2.0e-11_dp
  real(dp) :: errors(number_of_resolutions), conservation(number_of_resolutions)
  real(dp) :: closure(number_of_resolutions), minimum_theta(number_of_resolutions)
  real(dp) :: order, uncorrected_error, dummy_conservation, dummy_closure
  real(dp) :: dummy_theta
  integer :: case_index, fallback_count, dummy_fallback
  logical :: ok

  do case_index = 1, number_of_resolutions
    call run_case( &
      resolutions(case_index), .true., errors(case_index), &
      conservation(case_index), closure(case_index), minimum_theta(case_index), &
      fallback_count, ok)
    if (.not. ok) error stop "2D multispecies advection run failed"

    write(*, '(a,i0,4(a,es24.16),a,i0)') &
      "nx=", resolutions(case_index), &
      ", Y1 L1=", errors(case_index), &
      ", conservation=", conservation(case_index), &
      ", closure=", closure(case_index), &
      ", min theta=", minimum_theta(case_index), &
      ", fallbacks=", fallback_count

    if (conservation(case_index) > conservation_tolerance) then
      error stop "2D multispecies conservation error exceeds threshold"
    end if
    if (closure(case_index) > closure_tolerance) then
      error stop "2D multispecies closure error exceeds threshold"
    end if
    if (fallback_count /= 0) then
      error stop "smooth 2D multispecies case used a face fallback"
    end if
  end do

  do case_index = 1, number_of_resolutions - 1
    order = log(errors(case_index) / errors(case_index + 1)) / log(2.0_dp)
    write(*, '(a,i0,a,es24.16)') &
      "refinement pair ", case_index, ", order=", order
    if (order < minimum_order) then
      error stop "2D multispecies observed order is below threshold"
    end if
  end do

  call run_case( &
    40, .false., uncorrected_error, dummy_conservation, dummy_closure, &
    dummy_theta, dummy_fallback, ok)
  if (.not. ok) error stop "2D multispecies no-transverse run failed"
  write(*, '(a,es24.16)') &
    "nx=40 no-transverse Y1 L1=", uncorrected_error

  if (errors(2) >= 0.75_dp * uncorrected_error) then
    error stop "species transverse correction did not improve diagonal advection"
  end if

  write(*, '(a)') "test_multispecies_advection_2d: PASS"

contains

  subroutine run_case( &
      nx, use_transverse, mass_fraction_error, conservation_error, &
      maximum_closure_error, minimum_theta, fallback_count, ok)
    integer, intent(in) :: nx
    logical, intent(in) :: use_transverse
    real(dp), intent(out) :: mass_fraction_error, conservation_error
    real(dp), intent(out) :: maximum_closure_error, minimum_theta
    integer, intent(out) :: fallback_count
    logical, intent(out) :: ok

    integer :: ny, nvar, i, j, species, step_fallbacks
    real(dp), parameter :: x_min = 0.0_dp, x_max = 1.0_dp
    real(dp), parameter :: y_min = 0.0_dp, y_max = 1.0_dp
    real(dp), parameter :: final_time = 0.2_dp
    real(dp), parameter :: gamma = 1.4_dp, cfl = 0.4_dp
    real(dp), parameter :: velocity_x = 1.0_dp, velocity_y = 0.5_dp
    real(dp), parameter :: pressure = 1.0_dp
    real(dp), parameter :: amplitude = 0.2_dp
    real(dp), allocatable :: x(:), y(:), conserved(:, :, :)
    real(dp) :: primitive(nprim), base_state(ncons)
    real(dp) :: mass_fractions(nspecies), recovered(nspecies)
    real(dp) :: initial_species_mass(nspecies), final_species_mass(nspecies)
    real(dp) :: dx, dy, time, dt, phase_x, phase_y, exact_y1, step_theta
    logical :: cell_ok, step_ok

    ny = nx
    nvar = multispecies_nvar(nspecies)
    allocate(x(nx), y(ny), conserved(nvar, nx, ny))
    call uniform_cell_centers_2d( &
      nx, ny, x_min, x_max, y_min, y_max, x, y, dx, dy)

    primitive(qrho) = 1.0_dp
    primitive(qu) = velocity_x
    primitive(qv) = velocity_y
    primitive(qw) = 0.0_dp
    primitive(qp) = pressure
    call primitive_to_conserved(primitive, gamma, base_state, cell_ok)
    if (.not. cell_ok) then
      ok = .false.
      return
    end if

    do j = 1, ny
      do i = 1, nx
        mass_fractions(1) = 0.5_dp + amplitude * &
          sin(2.0_dp * acos(-1.0_dp) * x(i)) * &
          sin(2.0_dp * acos(-1.0_dp) * y(j))
        mass_fractions(2) = 1.0_dp - mass_fractions(1)
        call multispecies_state_from_base( &
          base_state, mass_fractions, nspecies, gamma, &
          conserved(:, i, j), cell_ok)
        if (.not. cell_ok) then
          ok = .false.
          return
        end if
      end do
    end do

    initial_species_mass = 0.0_dp
    do species = 1, nspecies
      initial_species_mass(species) = &
        sum(conserved(species_component(species), :, :)) * dx * dy
    end do

    time = 0.0_dp
    minimum_theta = 1.0_dp
    fallback_count = 0
    do while (time < final_time)
      call compute_cfl_timestep_multispecies_2d( &
        conserved, nx, ny, nspecies, dx, dy, gamma, cfl, dt, step_ok)
      if (.not. step_ok) then
        ok = .false.
        return
      end if
      dt = min(dt, final_time - time)
      call advance_ctu_multispecies_2d( &
        conserved, nx, ny, nspecies, dx, dy, dt, gamma, "mc", "pelec", &
        use_transverse, step_ok, step_theta, step_fallbacks)
      if (.not. step_ok) then
        ok = .false.
        return
      end if
      minimum_theta = min(minimum_theta, step_theta)
      fallback_count = fallback_count + step_fallbacks
      time = time + dt
    end do

    final_species_mass = 0.0_dp
    do species = 1, nspecies
      final_species_mass(species) = &
        sum(conserved(species_component(species), :, :)) * dx * dy
    end do
    conservation_error = maxval(abs(final_species_mass - initial_species_mass))

    mass_fraction_error = 0.0_dp
    maximum_closure_error = 0.0_dp
    do j = 1, ny
      do i = 1, nx
        call mass_fractions_from_state( &
          conserved(:, i, j), nspecies, recovered, cell_ok)
        if (.not. cell_ok) then
          ok = .false.
          return
        end if
        phase_x = x_min + modulo( &
          x(i) - x_min - velocity_x * final_time, x_max - x_min)
        phase_y = y_min + modulo( &
          y(j) - y_min - velocity_y * final_time, y_max - y_min)
        exact_y1 = 0.5_dp + amplitude * &
          sin(2.0_dp * acos(-1.0_dp) * phase_x) * &
          sin(2.0_dp * acos(-1.0_dp) * phase_y)
        mass_fraction_error = mass_fraction_error + abs(recovered(1) - exact_y1)
        maximum_closure_error = max(maximum_closure_error, &
          species_closure_error(conserved(:, i, j), nspecies))
      end do
    end do
    mass_fraction_error = mass_fraction_error / real(nx * ny, dp)
    ok = mass_fraction_error > 0.0_dp
  end subroutine run_case

end program test_multispecies_advection_2d
