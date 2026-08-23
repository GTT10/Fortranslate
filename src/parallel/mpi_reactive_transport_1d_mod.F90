module mpi_reactive_transport_1d_mod
  use mpi_f08
  use precision_mod, only: dp
  use mpi_domain_1d_mod, only: &
    mpi_domain_1d, exchange_periodic_halo_1d, global_minimum_1d
  use nasa7_thermo_mod, only: nasa7_species
  use transport_database_mod, only: gas_transport_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive, &
    reactive_diffusive_flux_x, reactive_transport_timestep
  implicit none
  private

  public :: mpi_reactive_transport_timestep
  public :: advance_mpi_reactive_transport
  public :: advance_mpi_reactive_transport_adaptive

contains

  subroutine mpi_reactive_transport_timestep( &
      domain, species, transport, state, temperature, dx, transport_cfl, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, dt, maximum_diffusivity, ok)
    type(mpi_domain_1d), intent(in) :: domain
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, 0:), temperature(0:)
    real(dp), intent(in) :: dx, transport_cfl
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled
    real(dp), intent(out) :: dt, maximum_diffusivity
    logical, intent(out) :: ok

    real(dp) :: local_dt, local_maximum_diffusivity
    logical :: local_ok, global_ok, reduction_ok
    integer :: ierr

    dt = 0.0_dp
    maximum_diffusivity = 0.0_dp
    ok = valid_local_shape(domain, species, state, temperature)
    call collective_logical_and(domain, ok, global_ok, reduction_ok)
    if (.not. reduction_ok .or. .not. global_ok) then
      ok = .false.
      return
    end if

    call reactive_transport_timestep( &
      species, transport, state, temperature, domain%local_cells, dx, &
      transport_cfl, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, local_dt, local_maximum_diffusivity, &
      local_ok)
    call collective_logical_and(domain, local_ok, global_ok, reduction_ok)
    if (.not. reduction_ok .or. .not. global_ok) then
      ok = .false.
      return
    end if

    call global_minimum_1d(domain, local_dt, dt, reduction_ok)
    if (.not. reduction_ok) then
      ok = .false.
      return
    end if
    call MPI_Allreduce( &
      local_maximum_diffusivity, maximum_diffusivity, 1, &
      MPI_DOUBLE_PRECISION, MPI_MAX, domain%comm, ierr)
    ok = ierr == MPI_SUCCESS .and. dt > 0.0_dp .and. &
      maximum_diffusivity >= 0.0_dp
  end subroutine mpi_reactive_transport_timestep

  subroutine advance_mpi_reactive_transport( &
      domain, species, transport, state, temperature, dx, interval, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, ok)
    type(mpi_domain_1d), intent(in) :: domain
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(inout) :: state(:, 0:), temperature(0:)
    real(dp), intent(in) :: dx, interval
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    logical, intent(out) :: ok

    real(dp), allocatable :: initial_state(:, :), initial_temperature(:)
    real(dp), allocatable :: stage1_state(:, :), stage1_temperature(:)
    real(dp), allocatable :: euler2_state(:, :), euler2_temperature(:)
    real(dp), allocatable :: primitive(:)
    real(dp) :: local_temperature, sound_speed, temperature_guess
    logical :: local_ok, global_ok, reduction_ok
    integer :: i, nvar

    ok = interval >= 0.0_dp .and. dx > 0.0_dp .and. &
      valid_local_shape(domain, species, state, temperature)
    call collective_logical_and(domain, ok, global_ok, reduction_ok)
    if (.not. reduction_ok .or. .not. global_ok) then
      ok = .false.
      return
    end if
    if (interval <= 0.0_dp .or. .not. (viscosity_enabled .or. &
        thermal_conduction_enabled .or. species_diffusion_enabled)) then
      call exchange_transport_halos(domain, state, temperature, ok)
      return
    end if

    nvar = reactive_nvar(size(species))
    allocate(initial_state(nvar, 0:domain%local_cells + 1))
    allocate(initial_temperature(0:domain%local_cells + 1))
    allocate(stage1_state(nvar, 0:domain%local_cells + 1))
    allocate(stage1_temperature(0:domain%local_cells + 1))
    allocate(euler2_state(nvar, 0:domain%local_cells + 1))
    allocate(euler2_temperature(0:domain%local_cells + 1))
    allocate(primitive(reactive_nprim(size(species))))
    initial_state = state
    initial_temperature = temperature

    call distributed_transport_euler_update( &
      domain, species, transport, initial_state, initial_temperature, dx, &
      interval, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, stage1_state, &
      stage1_temperature, ok)
    if (.not. ok) return
    call distributed_transport_euler_update( &
      domain, species, transport, stage1_state, stage1_temperature, dx, &
      interval, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, euler2_state, &
      euler2_temperature, ok)
    if (.not. ok) return

    state = initial_state
    temperature = initial_temperature
    local_ok = .true.
    do i = 1, domain%local_cells
      state(:, i) = 0.5_dp * &
        (initial_state(:, i) + euler2_state(:, i))
      temperature_guess = 0.5_dp * &
        (initial_temperature(i) + euler2_temperature(i))
      call reactive_conserved_to_primitive( &
        species, state(:, i), temperature_guess, primitive, &
        local_temperature, sound_speed, local_ok)
      if (.not. local_ok) exit
      temperature(i) = local_temperature
    end do
    call collective_logical_and(domain, local_ok, global_ok, reduction_ok)
    if (.not. reduction_ok .or. .not. global_ok) then
      ok = .false.
      return
    end if
    call exchange_transport_halos(domain, state, temperature, ok)
  end subroutine advance_mpi_reactive_transport

  subroutine advance_mpi_reactive_transport_adaptive( &
      domain, species, transport, state, temperature, dx, &
      requested_interval, minimum_interval, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, accepted_interval, rejected_trials, ok)
    type(mpi_domain_1d), intent(in) :: domain
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(inout) :: state(:, 0:), temperature(0:)
    real(dp), intent(in) :: dx, requested_interval, minimum_interval
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: accepted_interval
    integer, intent(out) :: rejected_trials
    logical, intent(out) :: ok

    real(dp), allocatable :: trial_state(:, :), trial_temperature(:)
    real(dp) :: trial_interval
    logical :: local_ok, global_ok, reduction_ok
    integer :: nvar

    accepted_interval = 0.0_dp
    rejected_trials = 0
    local_ok = requested_interval >= 0.0_dp .and. &
      minimum_interval > 0.0_dp .and. dx > 0.0_dp .and. &
      valid_local_shape(domain, species, state, temperature)
    call collective_logical_and(domain, local_ok, global_ok, reduction_ok)
    if (.not. reduction_ok .or. .not. global_ok) then
      ok = .false.
      return
    end if
    if (requested_interval <= 0.0_dp) then
      call exchange_transport_halos(domain, state, temperature, ok)
      return
    end if

    nvar = reactive_nvar(size(species))
    allocate(trial_state(nvar, 0:domain%local_cells + 1))
    allocate(trial_temperature(0:domain%local_cells + 1))
    trial_interval = requested_interval
    do
      trial_state = state
      trial_temperature = temperature
      call advance_mpi_reactive_transport( &
        domain, species, transport, trial_state, trial_temperature, dx, &
        trial_interval, viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, local_ok)
      if (local_ok) then
        state = trial_state
        temperature = trial_temperature
        accepted_interval = trial_interval
        ok = .true.
        return
      end if
      rejected_trials = rejected_trials + 1
      if (trial_interval <= minimum_interval) exit
      trial_interval = max(minimum_interval, 0.5_dp * trial_interval)
    end do
    ok = .false.
  end subroutine advance_mpi_reactive_transport_adaptive

  subroutine distributed_transport_euler_update( &
      domain, species, transport, input_state, input_temperature, dx, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, output_state, &
      output_temperature, ok)
    type(mpi_domain_1d), intent(in) :: domain
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: input_state(:, 0:), input_temperature(0:)
    real(dp), intent(in) :: dx, dt
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: output_state(:, 0:), output_temperature(0:)
    logical, intent(out) :: ok

    real(dp), allocatable :: work_state(:, :), work_temperature(:)
    real(dp), allocatable :: flux(:, :), primitive(:)
    real(dp) :: local_temperature, sound_speed
    logical :: local_ok, face_ok, global_ok, reduction_ok
    integer :: i, nvar

    nvar = reactive_nvar(size(species))
    local_ok = dt > 0.0_dp .and. dx > 0.0_dp .and. &
      valid_local_shape(domain, species, input_state, input_temperature) .and. &
      size(output_state, 1) == nvar .and. &
      size(output_state, 2) == domain%local_cells + 2 .and. &
      size(output_temperature) == domain%local_cells + 2
    call collective_logical_and(domain, local_ok, global_ok, reduction_ok)
    if (.not. reduction_ok .or. .not. global_ok) then
      ok = .false.
      return
    end if

    allocate(work_state(nvar, 0:domain%local_cells + 1))
    allocate(work_temperature(0:domain%local_cells + 1))
    allocate(flux(nvar, 0:domain%local_cells))
    allocate(primitive(reactive_nprim(size(species))))
    work_state = input_state
    work_temperature = input_temperature
    call exchange_transport_halos( &
      domain, work_state, work_temperature, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if

    local_ok = .true.
    do i = 0, domain%local_cells
      call reactive_diffusive_flux_x( &
        species, transport, work_state(:, i), work_state(:, i + 1), &
        work_temperature(i), work_temperature(i + 1), dx, &
        viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, flux(:, i), &
        face_ok)
      if (.not. face_ok) then
        local_ok = .false.
        flux(:, i) = 0.0_dp
      end if
    end do
    call collective_logical_and(domain, local_ok, global_ok, reduction_ok)
    if (.not. reduction_ok .or. .not. global_ok) then
      ok = .false.
      return
    end if

    output_state = work_state
    output_temperature = work_temperature
    local_ok = .true.
    do i = 1, domain%local_cells
      output_state(:, i) = work_state(:, i) - dt / dx * &
        (flux(:, i) - flux(:, i - 1))
      call reactive_conserved_to_primitive( &
        species, output_state(:, i), work_temperature(i), primitive, &
        local_temperature, sound_speed, face_ok)
      if (face_ok) then
        output_temperature(i) = local_temperature
      else
        local_ok = .false.
      end if
    end do
    call collective_logical_and(domain, local_ok, global_ok, reduction_ok)
    if (.not. reduction_ok .or. .not. global_ok) then
      ok = .false.
      return
    end if
    call exchange_transport_halos( &
      domain, output_state, output_temperature, ok)
  end subroutine distributed_transport_euler_update

  subroutine exchange_transport_halos(domain, state, temperature, ok)
    type(mpi_domain_1d), intent(in) :: domain
    real(dp), intent(inout) :: state(:, 0:), temperature(0:)
    logical, intent(out) :: ok

    real(dp), allocatable :: temperature_state(:, :)
    logical :: state_ok, temperature_ok

    allocate(temperature_state(1, 0:domain%local_cells + 1))
    call exchange_periodic_halo_1d(domain, state, state_ok)
    temperature_state(1, :) = temperature
    call exchange_periodic_halo_1d( &
      domain, temperature_state, temperature_ok)
    temperature = temperature_state(1, :)
    ok = state_ok .and. temperature_ok
  end subroutine exchange_transport_halos

  logical function valid_local_shape( &
      domain, species, state, temperature) result(valid)
    type(mpi_domain_1d), intent(in) :: domain
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, 0:), temperature(0:)

    valid = domain%local_cells > 0 .and. &
      size(state, 1) == reactive_nvar(size(species)) .and. &
      size(state, 2) == domain%local_cells + 2 .and. &
      size(temperature) == domain%local_cells + 2
  end function valid_local_shape

  subroutine collective_logical_and( &
      domain, local_value, global_value, ok)
    type(mpi_domain_1d), intent(in) :: domain
    logical, intent(in) :: local_value
    logical, intent(out) :: global_value, ok
    integer :: ierr

    call MPI_Allreduce( &
      local_value, global_value, 1, MPI_LOGICAL, MPI_LAND, domain%comm, ierr)
    ok = ierr == MPI_SUCCESS
  end subroutine collective_logical_and

end module mpi_reactive_transport_1d_mod
