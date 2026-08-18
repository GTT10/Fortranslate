program test_multispecies_entropy_wave
  use precision_mod, only: dp
  use state_indices_mod, only: ncons, nbase, nprim, qrho, qu, qv, qw, qp
  use mesh_mod, only: uniform_cell_centers
  use state_conversion_mod, only: primitive_to_conserved
  use multispecies_state_mod, only: &
    multispecies_state_from_base, mass_fractions_from_state, species_component
  use boundary_conditions_multispecies_mod, only: &
    apply_multispecies_boundary_conditions
  use time_integrator_multispecies_mod, only: &
    compute_multispecies_cfl_timestep, advance_multispecies_hydro_step, &
    maximum_species_closure_error
  implicit none

  integer, parameter :: nspecies = 2
  integer, parameter :: number_of_resolutions = 3
  integer, parameter :: resolutions(number_of_resolutions) = [40, 80, 160]
  real(dp), parameter :: minimum_order = 1.8_dp
  real(dp) :: errors(number_of_resolutions), order
  logical :: ok
  integer :: case_index

  do case_index = 1, number_of_resolutions
    call run_case(resolutions(case_index), errors(case_index), ok)
    if (.not. ok) error stop "Multispecies entropy-wave run failed"
    write(*, '(a,i0,a,es24.16)') &
      "nx=", resolutions(case_index), ", Y1 L1=", errors(case_index)
  end do

  do case_index = 1, number_of_resolutions - 1
    order = log(errors(case_index) / errors(case_index + 1)) / log(2.0_dp)
    write(*, '(a,i0,a,es24.16)') &
      "refinement pair ", case_index, ", order=", order
    if (order < minimum_order) error stop "Species order below threshold"
  end do

  write(*, '(a)') "test_multispecies_entropy_wave: PASS"

contains

  subroutine run_case(nx, error, ok)
    integer, intent(in) :: nx
    real(dp), intent(out) :: error
    logical, intent(out) :: ok

    real(dp), parameter :: x_min = 0.0_dp, x_max = 1.0_dp
    real(dp), parameter :: final_time = 0.1_dp
    real(dp), parameter :: gamma = 1.4_dp, cfl = 0.40_dp
    real(dp), parameter :: velocity = 1.0_dp
    real(dp), parameter :: amplitude = 0.2_dp
    real(dp), allocatable :: x(:), conserved(:, :)
    real(dp) :: primitive(nprim), base_state(ncons), y(nspecies), recovered(nspecies)
    real(dp) :: dx, time, dt, phase, exact_y1
    real(dp) :: initial_species_mass(nspecies), final_species_mass(nspecies)
    logical :: cell_ok, boundary_ok
    integer :: i, species

    allocate(x(nx))
    allocate(conserved(nbase + nspecies, 0:nx + 1))
    call uniform_cell_centers(nx, x_min, x_max, x, dx)

    primitive(qrho) = 1.0_dp
    primitive(qu) = velocity
    primitive(qv) = 0.0_dp
    primitive(qw) = 0.0_dp
    primitive(qp) = 1.0_dp
    call primitive_to_conserved(primitive, gamma, base_state, cell_ok)
    if (.not. cell_ok) then
      ok = .false.
      error = huge(1.0_dp)
      return
    end if

    do i = 1, nx
      y(1) = 0.5_dp + amplitude * sin(2.0_dp * acos(-1.0_dp) * x(i))
      y(2) = 1.0_dp - y(1)
      call multispecies_state_from_base( &
        base_state, y, nspecies, gamma, conserved(:, i), cell_ok)
      if (.not. cell_ok) then
        ok = .false.
        error = huge(1.0_dp)
        return
      end if
    end do
    call apply_multispecies_boundary_conditions( &
      conserved, nx, "periodic", boundary_ok)
    if (.not. boundary_ok) then
      ok = .false.
      error = huge(1.0_dp)
      return
    end if

    do species = 1, nspecies
      initial_species_mass(species) = &
        sum(conserved(species_component(species), 1:nx)) * dx
    end do

    time = 0.0_dp
    do while (time < final_time)
      call compute_multispecies_cfl_timestep( &
        conserved, nx, nspecies, dx, gamma, cfl, dt, cell_ok)
      if (.not. cell_ok) then
        ok = .false.
        error = huge(1.0_dp)
        return
      end if
      dt = min(dt, final_time - time)
      call advance_multispecies_hydro_step( &
        conserved, nx, nspecies, dx, dt, gamma, cell_ok, &
        reconstruction="pelec_plm", limiter="mc", &
        boundary_condition="periodic", riemann_solver="pelec", &
        plm_order=2, use_flattening=.false.)
      if (.not. cell_ok) then
        ok = .false.
        error = huge(1.0_dp)
        return
      end if
      time = time + dt
    end do

    error = 0.0_dp
    do i = 1, nx
      call mass_fractions_from_state( &
        conserved(:, i), nspecies, recovered, cell_ok)
      if (.not. cell_ok) then
        ok = .false.
        error = huge(1.0_dp)
        return
      end if
      phase = x_min + modulo(x(i) - x_min - velocity * final_time, x_max - x_min)
      exact_y1 = 0.5_dp + amplitude * &
        sin(2.0_dp * acos(-1.0_dp) * phase)
      error = error + abs(recovered(1) - exact_y1)
    end do
    error = error / real(nx, dp)

    do species = 1, nspecies
      final_species_mass(species) = &
        sum(conserved(species_component(species), 1:nx)) * dx
      if (abs(final_species_mass(species) - initial_species_mass(species)) > &
          2.0e-12_dp) then
        ok = .false.
        return
      end if
    end do
    if (maximum_species_closure_error(conserved, nx, nspecies) > 2.0e-12_dp) then
      ok = .false.
      return
    end if
    ok = error > 0.0_dp .and. error < huge(1.0_dp)
  end subroutine run_case

end program test_multispecies_entropy_wave
