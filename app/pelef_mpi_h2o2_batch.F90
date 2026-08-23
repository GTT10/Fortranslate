program pelef_mpi_h2o2_batch
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use mpi_f08
  use mpi_domain_1d_mod, only: &
    mpi_domain_1d, initialize_mpi_domain_1d, gather_state_1d
  use nasa7_thermo_mod, only: nasa7_species
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density, mixture_pressure, &
    mixture_mass_properties
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_full_thermo_mod, only: full_nspecies, load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use constant_volume_reactor_mod, only: &
    advance_constant_volume_implicit_adaptive
  implicit none

  integer, parameter :: global_reactors = 11
  integer, parameter :: metadata_fields = 14
  integer, parameter :: result_fields = metadata_fields + full_nspecies
  integer, parameter :: i_initial_temperature = 1
  integer, parameter :: i_final_temperature = 2
  integer, parameter :: i_initial_pressure = 3
  integer, parameter :: i_final_pressure = 4
  integer, parameter :: i_density = 5
  integer, parameter :: i_target_energy = 6
  integer, parameter :: i_final_time = 7
  integer, parameter :: i_energy_error = 8
  integer, parameter :: i_closure_error = 9
  integer, parameter :: i_minimum_species = 10
  integer, parameter :: i_species_change = 11
  integer, parameter :: i_steps = 12
  integer, parameter :: i_newton_iterations = 13
  integer, parameter :: i_rejected_steps = 14
  integer, parameter :: first_species = metadata_fields + 1
  real(dp), parameter :: batch_final_time = 5.0e-7_dp
  real(dp), parameter :: initial_time_step = 1.0e-7_dp
  real(dp), parameter :: maximum_time_step = 2.0e-7_dp
  real(dp), parameter :: minimum_time_step = 1.0e-18_dp
  real(dp), parameter :: relative_tolerance = 1.0e-6_dp
  real(dp), parameter :: absolute_tolerance = 1.0e-12_dp
  integer, parameter :: maximum_steps = 10000

  type(mpi_domain_1d) :: domain
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp), allocatable :: local_results(:, :), global_results(:, :)
  character(len=256) :: output_file
  integer :: ierr, cell, global_cell, output_unit
  logical :: ok

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop 'MPI_Init failed'
  call initialize_mpi_domain_1d( &
    domain, global_reactors, MPI_COMM_WORLD, ok)
  if (.not. ok) error stop 'MPI reactor decomposition failed'

  if (command_argument_count() >= 1) then
    call get_command_argument(1, output_file)
  else
    output_file = 'mpi_h2o2_batch.csv'
  end if

  call load_h2o2_full_thermo(species, ok)
  if (.not. ok) error stop 'Failed to load full H2/O2 thermodynamics'
  call load_h2o2_full_mechanism(reactions, ok)
  if (.not. ok) error stop 'Failed to load full H2/O2 mechanism'

  allocate(local_results(result_fields, domain%local_cells))
  do cell = 1, domain%local_cells
    global_cell = domain%global_first + cell - 1
    call integrate_reactor(global_cell, local_results(:, cell))
  end do

  call gather_state_1d( &
    domain, local_results, global_results, 0, ok)
  if (.not. ok) error stop 'MPI chemistry result gather failed'

  if (domain%rank == 0) then
    call validate_batch(global_results)
    open(newunit=output_unit, file=trim(output_file), status='replace', &
      action='write')
    write(output_unit, '(a)') &
      'reactor,initial_temperature,final_temperature,initial_pressure,' // &
      'final_pressure,density,target_internal_energy,final_time,' // &
      'relative_energy_error,closure_error,minimum_mass_fraction,' // &
      'maximum_species_change,accepted_steps,newton_iterations,' // &
      'rejected_trials,Y_H2,Y_H,Y_O,Y_O2,Y_OH,Y_H2O,Y_HO2,Y_H2O2,Y_AR,Y_N2'
    do cell = 1, domain%global_cells
      write(output_unit, '(*(es25.16e3,:,","))') &
        real(cell, dp), global_results(:, cell)
    end do
    close(output_unit)
    write(*, '(a,i0)') 'MPI ranks: ', domain%nranks
    write(*, '(a,i0)') 'Distributed reactors: ', domain%global_cells
    write(*, '(a,i0)') 'Accepted implicit steps: ', &
      nint(sum(global_results(i_steps, :)))
    write(*, '(a,i0)') 'Newton iterations: ', &
      nint(sum(global_results(i_newton_iterations, :)))
    write(*, '(a,es24.16)') 'Maximum species change: ', &
      maxval(global_results(i_species_change, :))
  end if

  call MPI_Finalize(ierr)
  if (ierr /= MPI_SUCCESS) error stop 'MPI_Finalize failed'

contains

  subroutine integrate_reactor(reactor_index, result)
    integer, intent(in) :: reactor_index
    real(dp), intent(out) :: result(result_fields)

    real(dp) :: mole_fractions(full_nspecies)
    real(dp) :: mass_fractions(full_nspecies)
    real(dp) :: initial_mass_fractions(full_nspecies)
    real(dp) :: initial_temperature, temperature
    real(dp) :: initial_pressure, pressure, density
    real(dp) :: target_energy, final_energy
    real(dp) :: molecular_weight, gas_constant, cp, cv, gamma
    real(dp) :: enthalpy, entropy
    real(dp) :: time, requested_step, accepted_step, suggested_step
    real(dp) :: time_tolerance, relative_energy_error, closure_error
    real(dp) :: phase
    integer :: steps, total_newton, total_rejected
    integer :: newton_iterations, rejected_steps
    logical :: reactor_ok

    phase = 2.0_dp * acos(-1.0_dp) * &
      (real(reactor_index, dp) - 0.5_dp) / real(global_reactors, dp)
    mole_fractions = [ &
      2.0_dp + 0.12_dp * sin(phase), &
      1.0e-6_dp * (1.0_dp + 0.20_dp * cos(phase)), &
      1.0e-12_dp, &
      1.0_dp + 0.08_dp * cos(phase), &
      1.0e-10_dp, &
      0.0_dp, &
      0.0_dp, &
      0.0_dp, &
      0.10_dp, &
      3.0_dp]
    mole_fractions = mole_fractions / sum(mole_fractions)
    call mass_fractions_from_mole_fractions( &
      species, mole_fractions, mass_fractions, reactor_ok)
    if (.not. reactor_ok) error stop 'MPI reactor composition failed'

    initial_temperature = 950.0_dp + 300.0_dp * &
      real(reactor_index - 1, dp) / real(global_reactors - 1, dp)
    initial_pressure = 101325.0_dp * &
      (0.80_dp + 0.40_dp * real(reactor_index - 1, dp) / &
      real(global_reactors - 1, dp))
    temperature = initial_temperature
    initial_mass_fractions = mass_fractions
    density = mixture_density( &
      species, mass_fractions, initial_pressure, temperature, reactor_ok)
    if (.not. reactor_ok) error stop 'MPI reactor density failed'
    call mixture_mass_properties( &
      species, mass_fractions, temperature, molecular_weight, gas_constant, &
      cp, cv, gamma, enthalpy, target_energy, entropy, reactor_ok)
    if (.not. reactor_ok) error stop 'MPI reactor initial energy failed'

    time = 0.0_dp
    requested_step = initial_time_step
    steps = 0
    total_newton = 0
    total_rejected = 0
    time_tolerance = 100.0_dp * epsilon(1.0_dp) * batch_final_time

    do while (time < batch_final_time - time_tolerance)
      if (steps >= maximum_steps) error stop 'MPI reactor step limit reached'
      requested_step = min(requested_step, maximum_time_step)
      requested_step = min(requested_step, batch_final_time - time)
      call advance_constant_volume_implicit_adaptive( &
        species, reactions, density, target_energy, requested_step, &
        relative_tolerance, absolute_tolerance, mass_fractions, temperature, &
        accepted_step, suggested_step, newton_iterations, rejected_steps, &
        reactor_ok)
      if (.not. reactor_ok .or. accepted_step < minimum_time_step) then
        error stop 'MPI implicit chemistry step failed'
      end if
      time = time + accepted_step
      requested_step = max(minimum_time_step, &
        min(maximum_time_step, suggested_step))
      steps = steps + 1
      total_newton = total_newton + newton_iterations
      total_rejected = total_rejected + rejected_steps
    end do

    call mixture_mass_properties( &
      species, mass_fractions, temperature, molecular_weight, gas_constant, &
      cp, cv, gamma, enthalpy, final_energy, entropy, reactor_ok)
    if (.not. reactor_ok) error stop 'MPI reactor final energy failed'
    pressure = mixture_pressure( &
      species, mass_fractions, density, temperature, reactor_ok)
    if (.not. reactor_ok) error stop 'MPI reactor final pressure failed'
    relative_energy_error = abs(final_energy - target_energy) / &
      max(1.0_dp, abs(target_energy))
    closure_error = abs(sum(mass_fractions) - 1.0_dp)

    result = 0.0_dp
    result(i_initial_temperature) = initial_temperature
    result(i_final_temperature) = temperature
    result(i_initial_pressure) = initial_pressure
    result(i_final_pressure) = pressure
    result(i_density) = density
    result(i_target_energy) = target_energy
    result(i_final_time) = time
    result(i_energy_error) = relative_energy_error
    result(i_closure_error) = closure_error
    result(i_minimum_species) = minval(mass_fractions)
    result(i_species_change) = &
      maxval(abs(mass_fractions - initial_mass_fractions))
    result(i_steps) = real(steps, dp)
    result(i_newton_iterations) = real(total_newton, dp)
    result(i_rejected_steps) = real(total_rejected, dp)
    result(first_species:result_fields) = mass_fractions
  end subroutine integrate_reactor

  subroutine validate_batch(results)
    real(dp), intent(in) :: results(:, :)

    if (size(results, 1) /= result_fields .or. &
        size(results, 2) /= global_reactors) then
      error stop 'MPI chemistry batch shape regression failed'
    end if
    if (.not. all(ieee_is_finite(results))) then
      error stop 'MPI chemistry batch produced non-finite data'
    end if
    if (maxval(results(i_energy_error, :)) > 5.0e-9_dp) then
      error stop 'MPI chemistry energy regression failed'
    end if
    if (maxval(results(i_closure_error, :)) > 5.0e-12_dp) then
      error stop 'MPI chemistry closure regression failed'
    end if
    if (minval(results(i_minimum_species, :)) < -2.0e-13_dp) then
      error stop 'MPI chemistry positivity regression failed'
    end if
    if (maxval(abs(results(i_final_time, :) - batch_final_time)) > &
        5.0e-15_dp) then
      error stop 'MPI chemistry final-time regression failed'
    end if
    if (minval(results(i_steps, :)) < 1.0_dp .or. &
        maxval(results(i_newton_iterations, :)) < 1.0_dp) then
      error stop 'MPI chemistry scheduler did not execute'
    end if
    if (maxval(results(i_species_change, :)) <= 1.0e-16_dp) then
      error stop 'MPI chemistry response was trivial'
    end if
  end subroutine validate_batch

end program pelef_mpi_h2o2_batch
