program pelef_mpi_reactive_transport_1d
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mpi_f08
  use precision_mod, only: dp
  use state_indices_mod, only: irho
  use mpi_domain_1d_mod, only: &
    mpi_domain_1d, initialize_mpi_domain_1d, global_sum_1d, gather_state_1d
  use mpi_reactive_transport_1d_mod, only: &
    mpi_reactive_transport_timestep, advance_mpi_reactive_transport
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: full_nspecies, load_h2o2_full_thermo
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_full_transport
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_mass_fraction_component, reactive_primitive_to_conserved, &
    reactive_conserved_to_primitive
  implicit none

  integer, parameter :: global_cells = 33
  real(dp), parameter :: domain_length = 1.0e-3_dp
  real(dp), parameter :: final_time = 1.0e-6_dp
  real(dp), parameter :: transport_cfl = 0.05_dp
  real(dp), parameter :: conservation_tolerance = 2.0e-10_dp
  real(dp), parameter :: closure_tolerance = 2.0e-10_dp
  real(dp), parameter :: species_tolerance = 2.0e-12_dp
  integer, parameter :: maximum_steps = 100000
  logical, parameter :: viscosity_enabled = .true.
  logical, parameter :: thermal_conduction_enabled = .true.
  logical, parameter :: species_diffusion_enabled = .true.
  logical, parameter :: barodiffusion_enabled = .true.

  type(mpi_domain_1d) :: domain
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  real(dp), allocatable :: state(:, :), initial_state(:, :)
  real(dp), allocatable :: temperature(:), initial_temperature(:)
  real(dp), allocatable :: local_output(:, :), global_output(:, :)
  real(dp) :: dx, time, dt, maximum_diffusivity, time_tolerance
  real(dp), allocatable :: initial_integrals(:), final_integrals(:)
  real(dp) :: maximum_conservation_error, maximum_closure_error
  real(dp) :: minimum_species, maximum_state_change
  character(len=256) :: output_file
  integer :: ierr, nvar, output_fields, cell, global_cell, steps, output_unit
  logical :: ok

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop 'MPI_Init failed'
  call initialize_mpi_domain_1d( &
    domain, global_cells, MPI_COMM_WORLD, ok)
  if (.not. ok) error stop 'MPI transport decomposition failed'

  if (command_argument_count() >= 1) then
    call get_command_argument(1, output_file)
  else
    output_file = 'mpi_reactive_transport.csv'
  end if

  call load_h2o2_full_thermo(species, ok)
  if (.not. ok) error stop 'Failed to load full H2/O2 thermodynamics'
  call load_h2o2_full_transport(transport, ok)
  if (.not. ok) error stop 'Failed to load full H2/O2 transport data'

  nvar = reactive_nvar(full_nspecies)
  output_fields = nvar + 1
  dx = domain_length / real(global_cells, dp)
  allocate(state(nvar, 0:domain%local_cells + 1))
  allocate(initial_state(nvar, 0:domain%local_cells + 1))
  allocate(temperature(0:domain%local_cells + 1))
  allocate(initial_temperature(0:domain%local_cells + 1))
  allocate(initial_integrals(nvar), final_integrals(nvar))
  state = 0.0_dp
  temperature = 900.0_dp

  do cell = 1, domain%local_cells
    global_cell = domain%global_first + cell - 1
    call initialize_cell( &
      (real(global_cell, dp) - 0.5_dp) * dx, state(:, cell), &
      temperature(cell))
  end do
  initial_state = state
  initial_temperature = temperature
  call validate_local_state( &
    state, temperature, maximum_closure_error, minimum_species)
  call global_integrals(state, initial_integrals)

  time = 0.0_dp
  steps = 0
  time_tolerance = 100.0_dp * epsilon(1.0_dp) * final_time
  do while (time < final_time - time_tolerance)
    if (steps >= maximum_steps) error stop 'MPI transport step limit reached'
    call mpi_reactive_transport_timestep( &
      domain, species, transport, state, temperature, dx, transport_cfl, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, dt, maximum_diffusivity, ok)
    if (.not. ok .or. dt <= 0.0_dp) then
      error stop 'MPI transport timestep failed'
    end if
    dt = min(dt, final_time - time)
    call advance_mpi_reactive_transport( &
      domain, species, transport, state, temperature, dx, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, ok)
    if (.not. ok) error stop 'MPI transport update failed'
    time = time + dt
    steps = steps + 1
  end do

  call validate_local_state( &
    state, temperature, maximum_closure_error, minimum_species)
  call global_integrals(state, final_integrals)
  maximum_conservation_error = maxval(abs(final_integrals - &
    initial_integrals) / max(1.0_dp, abs(initial_integrals)))
  if (maximum_conservation_error > conservation_tolerance) then
    error stop 'MPI transport conservation regression failed'
  end if
  call global_maximum_state_change(maximum_state_change)
  if (maximum_state_change <= 1.0e-10_dp) then
    error stop 'MPI transport response was trivial'
  end if
  if (abs(time - final_time) > 5.0e-15_dp .or. steps < 1 .or. &
      maximum_diffusivity <= 0.0_dp) then
    error stop 'MPI transport scheduler regression failed'
  end if

  allocate(local_output(output_fields, domain%local_cells))
  local_output(1:nvar, :) = state(:, 1:domain%local_cells)
  local_output(output_fields, :) = temperature(1:domain%local_cells)
  call gather_state_1d( &
    domain, local_output, global_output, 0, ok)
  if (.not. ok) error stop 'MPI transport result gather failed'

  if (domain%rank == 0) then
    open(newunit=output_unit, file=trim(output_file), status='replace', &
      action='write')
    write(output_unit, '(a)') &
      'x,rho,rhou,rhov,rhow,rhoE,temperature,' // &
      'rhoY_H2,rhoY_H,rhoY_O,rhoY_O2,rhoY_OH,rhoY_H2O,' // &
      'rhoY_HO2,rhoY_H2O2,rhoY_AR,rhoY_N2'
    do cell = 1, global_cells
      write(output_unit, '(*(es25.16e3,:,","))') &
        (real(cell, dp) - 0.5_dp) * dx, &
        global_output(1:5, cell), global_output(output_fields, cell), &
        global_output(6:nvar, cell)
    end do
    close(output_unit)
    write(*, '(a,i0)') 'MPI ranks: ', domain%nranks
    write(*, '(a,i0)') 'Distributed transport cells: ', global_cells
    write(*, '(a,i0)') 'Accepted SSPRK2 steps: ', steps
    write(*, '(a,es24.16)') 'Final time: ', time
    write(*, '(a,es24.16)') 'Maximum diffusivity: ', maximum_diffusivity
    write(*, '(a,es24.16)') 'Maximum conservation error: ', &
      maximum_conservation_error
    write(*, '(a,es24.16)') 'Maximum state change: ', maximum_state_change
  end if

  call MPI_Finalize(ierr)
  if (ierr /= MPI_SUCCESS) error stop 'MPI_Finalize failed'

contains

  subroutine initialize_cell(position, conserved, cell_temperature)
    real(dp), intent(in) :: position
    real(dp), intent(out) :: conserved(:), cell_temperature

    real(dp) :: mole_fractions(full_nspecies)
    real(dp) :: mass_fractions(full_nspecies)
    real(dp) :: primitive(reactive_nprim(full_nspecies))
    real(dp) :: pressure, density, sound_speed, phase
    integer :: k
    logical :: local_ok

    phase = 2.0_dp * acos(-1.0_dp) * position / domain_length
    mole_fractions = [ &
      2.0_dp + 0.25_dp * sin(phase), &
      0.020_dp + 0.002_dp * cos(phase), &
      0.015_dp + 0.001_dp * sin(2.0_dp * phase), &
      1.0_dp + 0.10_dp * cos(phase), &
      0.020_dp + 0.002_dp * sin(phase), &
      0.20_dp + 0.02_dp * cos(phase), &
      0.010_dp + 0.001_dp * sin(phase), &
      0.010_dp + 0.001_dp * cos(phase), &
      0.20_dp, &
      3.0_dp]
    mole_fractions = mole_fractions / sum(mole_fractions)
    call mass_fractions_from_mole_fractions( &
      species, mole_fractions, mass_fractions, local_ok)
    if (.not. local_ok) error stop 'MPI transport composition failed'

    cell_temperature = 900.0_dp + 80.0_dp * sin(phase) + &
      20.0_dp * cos(2.0_dp * phase)
    pressure = 101325.0_dp * (1.0_dp + 0.01_dp * cos(phase))
    density = mixture_density( &
      species, mass_fractions, pressure, cell_temperature, local_ok)
    if (.not. local_ok) error stop 'MPI transport density failed'

    primitive = 0.0_dp
    primitive(1:5) = [ &
      density, 12.0_dp * sin(phase), 4.0_dp * cos(phase), &
      2.0_dp * sin(2.0_dp * phase), pressure]
    do k = 1, full_nspecies
      primitive(reactive_mass_fraction_component(k)) = mass_fractions(k)
    end do
    call reactive_primitive_to_conserved( &
      species, primitive, conserved, cell_temperature, sound_speed, local_ok)
    if (.not. local_ok) error stop 'MPI transport state initialization failed'
  end subroutine initialize_cell

  subroutine validate_local_state( &
      local_state, local_temperature, closure_error, local_minimum_species)
    real(dp), intent(in) :: local_state(:, 0:), local_temperature(0:)
    real(dp), intent(out) :: closure_error, local_minimum_species

    real(dp) :: primitive(reactive_nprim(full_nspecies))
    real(dp) :: recovered_temperature, sound_speed, rho, local_closure
    real(dp) :: global_values(2), local_values(2)
    logical :: local_ok, reduction_ok
    integer :: i, k

    closure_error = 0.0_dp
    local_minimum_species = huge(1.0_dp)
    local_ok = .true.
    do i = 1, domain%local_cells
      call reactive_conserved_to_primitive( &
        species, local_state(:, i), local_temperature(i), primitive, &
        recovered_temperature, sound_speed, reduction_ok)
      if (.not. reduction_ok .or. .not. ieee_is_finite(recovered_temperature)) then
        local_ok = .false.
        cycle
      end if
      rho = local_state(irho, i)
      local_closure = 0.0_dp
      do k = 1, full_nspecies
        local_closure = local_closure + &
          local_state(reactive_species_component(k), i)
        local_minimum_species = min(local_minimum_species, &
          local_state(reactive_species_component(k), i))
      end do
      closure_error = max(closure_error, abs(local_closure - rho) / rho)
    end do
    local_values = [closure_error, -local_minimum_species]
    call MPI_Allreduce( &
      local_values, global_values, 2, MPI_DOUBLE_PRECISION, MPI_MAX, &
      domain%comm, ierr)
    if (ierr /= MPI_SUCCESS) error stop 'MPI state validation reduction failed'
    closure_error = global_values(1)
    local_minimum_species = -global_values(2)
    if (.not. local_ok .or. closure_error > closure_tolerance .or. &
        local_minimum_species < -species_tolerance) then
      error stop 'MPI transport state validation failed'
    end if
  end subroutine validate_local_state

  subroutine global_integrals(local_state, global_values)
    real(dp), intent(in) :: local_state(:, 0:)
    real(dp), intent(out) :: global_values(:)

    real(dp) :: local_values(size(global_values))
    logical :: reduction_ok

    local_values = sum( &
      local_state(:, 1:domain%local_cells), dim=2) * dx
    call global_sum_1d(domain, local_values, global_values, reduction_ok)
    if (.not. reduction_ok) error stop 'MPI transport integral reduction failed'
  end subroutine global_integrals

  subroutine global_maximum_state_change(global_change)
    real(dp), intent(out) :: global_change

    real(dp) :: local_change
    integer :: reduction_ierr

    local_change = maxval(abs( &
      state(:, 1:domain%local_cells) - &
      initial_state(:, 1:domain%local_cells)) / &
      max(1.0_dp, abs(initial_state(:, 1:domain%local_cells))))
    call MPI_Allreduce( &
      local_change, global_change, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
      domain%comm, reduction_ierr)
    if (reduction_ierr /= MPI_SUCCESS) then
      error stop 'MPI transport change reduction failed'
    end if
  end subroutine global_maximum_state_change

end program pelef_mpi_reactive_transport_1d
