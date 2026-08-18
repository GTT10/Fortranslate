program pelef_ms
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
  use state_indices_mod, only: irho, imx, iet, ncons, nprim, qrho, qp
  use mesh_mod, only: uniform_cell_centers
  use multispecies_state_mod, only: &
    max_supported_species, multispecies_nvar, species_component, &
    mass_fractions_from_state
  use state_conversion_mod, only: conserved_to_primitive
  use simulation_config_multispecies_mod, only: &
    multispecies_simulation_config, multispec_sod_config, &
    read_multispecies_configuration
  use multispec_sod_problem_mod, only: initialize_multispec_sod_problem
  use time_integrator_multispecies_mod, only: &
    compute_multispecies_cfl_timestep, advance_multispecies_hydro_step, &
    maximum_species_closure_error
  use csv_io_multispecies_mod, only: write_multispecies_solution_csv
  implicit none

  type(multispecies_simulation_config) :: simulation
  type(multispec_sod_config) :: problem
  character(len=1024) :: input_path, message
  real(dp), allocatable :: x(:), conserved(:, :)
  real(dp), allocatable :: initial_totals(:), final_totals(:)
  real(dp) :: primitive(nprim), mass_fractions(max_supported_species)
  real(dp) :: dx, time, dt, minimum_density, minimum_pressure
  real(dp) :: minimum_mass_fraction, maximum_closure_error
  logical :: ok, cell_ok
  integer :: argument_count, step, i, species, nvar

  argument_count = command_argument_count()
  if (argument_count /= 1) then
    write(*, '(a)') "Usage: pelef_ms <input.nml>"
    error stop 2
  end if

  call get_command_argument(1, input_path)
  call read_multispecies_configuration( &
    trim(input_path), simulation, problem, ok, message)
  if (.not. ok) then
    write(*, '(a)') trim(message)
    error stop 2
  end if

  nvar = multispecies_nvar(simulation%nspecies)
  allocate(x(simulation%nx))
  allocate(conserved(nvar, 0:simulation%nx + 1))
  allocate(initial_totals(nvar), final_totals(nvar))

  call uniform_cell_centers( &
    simulation%nx, simulation%x_min, simulation%x_max, x, dx)
  call initialize_multispec_sod_problem( &
    x, simulation%nx, simulation%nspecies, simulation%gamma, problem, &
    conserved, ok)
  if (.not. ok) error stop "Failed to initialize multispecies Sod problem"

  initial_totals = sum(conserved(:, 1:simulation%nx), dim=2) * dx

  write(*, '(a)') "PeleF " // pelef_version
  write(*, '(a,1x,a)') "Input:", trim(input_path)
  write(*, '(a,1x,a)') "Problem:", trim(simulation%problem)
  write(*, '(a,i0)') "Species: ", simulation%nspecies
  write(*, '(a,i0,a,es12.5)') &
    "Grid: nx=", simulation%nx, ", dx=", dx
  write(*, '(a,1x,a)') "Reconstruction:", trim(simulation%reconstruction)
  write(*, '(a,1x,a)') "Riemann solver:", trim(simulation%riemann_solver)

  time = 0.0_dp
  step = 0
  do while (time < simulation%final_time)
    if (step >= simulation%max_steps) then
      error stop "Maximum step count reached before final_time"
    end if

    call compute_multispecies_cfl_timestep( &
      conserved, simulation%nx, simulation%nspecies, dx, simulation%gamma, &
      simulation%cfl, dt, ok)
    if (.not. ok) error stop "Failed to compute multispecies CFL timestep"

    dt = min(dt, simulation%final_time - time)
    call advance_multispecies_hydro_step( &
      conserved, simulation%nx, simulation%nspecies, dx, dt, &
      simulation%gamma, ok, reconstruction=simulation%reconstruction, &
      limiter=simulation%limiter, &
      boundary_condition=simulation%boundary_condition, &
      riemann_solver=simulation%riemann_solver, &
      plm_order=simulation%plm_order, &
      use_flattening=simulation%use_flattening)
    if (.not. ok) error stop "Non-physical multispecies hydro update"

    time = time + dt
    step = step + 1
  end do

  final_totals = sum(conserved(:, 1:simulation%nx), dim=2) * dx
  minimum_density = huge(1.0_dp)
  minimum_pressure = huge(1.0_dp)
  minimum_mass_fraction = huge(1.0_dp)
  do i = 1, simulation%nx
    call conserved_to_primitive( &
      conserved(1:ncons, i), simulation%gamma, primitive, cell_ok)
    if (.not. cell_ok) error stop "Final base state is not physical"
    call mass_fractions_from_state( &
      conserved(:, i), simulation%nspecies, mass_fractions, cell_ok)
    if (.not. cell_ok) error stop "Final species state is invalid"
    minimum_density = min(minimum_density, primitive(qrho))
    minimum_pressure = min(minimum_pressure, primitive(qp))
    minimum_mass_fraction = min( &
      minimum_mass_fraction, minval(mass_fractions(1:simulation%nspecies)))
  end do
  maximum_closure_error = maximum_species_closure_error( &
    conserved, simulation%nx, simulation%nspecies)

  call write_multispecies_solution_csv( &
    trim(simulation%output_file), x, conserved, simulation%nx, &
    simulation%nspecies, simulation%gamma, ok, message)
  if (.not. ok) then
    write(*, '(a)') trim(message)
    error stop 3
  end if

  write(*, '(a,i0)') "Completed steps: ", step
  write(*, '(a,es24.16)') "Final time: ", time
  write(*, '(a,es24.16)') "Minimum density: ", minimum_density
  write(*, '(a,es24.16)') "Minimum pressure: ", minimum_pressure
  write(*, '(a,es24.16)') "Minimum mass fraction: ", minimum_mass_fraction
  write(*, '(a,es24.16)') &
    "Maximum species closure error: ", maximum_closure_error
  write(*, '(a,es24.16)') &
    "Mass change: ", final_totals(irho) - initial_totals(irho)
  write(*, '(a,es24.16)') &
    "Momentum change: ", final_totals(imx) - initial_totals(imx)
  write(*, '(a,es24.16)') &
    "Energy change: ", final_totals(iet) - initial_totals(iet)
  do species = 1, simulation%nspecies
    write(*, '(a,i0,a,es24.16)') "Species ", species, " mass change: ", &
      final_totals(species_component(species)) - &
      initial_totals(species_component(species))
  end do
  write(*, '(a,1x,a)') "Output:", trim(simulation%output_file)
end program pelef_ms
