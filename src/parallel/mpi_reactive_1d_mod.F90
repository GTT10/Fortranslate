module mpi_reactive_1d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mpi_f08
  use precision_mod, only: dp
  use mpi_domain_1d_mod, only: &
    mpi_domain_1d, exchange_periodic_halo_1d
  use mpi_reactive_transport_1d_mod, only: &
    mpi_reactive_transport_timestep, advance_mpi_reactive_transport
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: gas_transport_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive, &
    reactive_riemann_flux_x, advance_reactive_chemistry
  implicit none
  private

  public :: mpi_reactive_timestep
  public :: advance_mpi_reactive_hydro
  public :: advance_mpi_reactive_chemistry
  public :: advance_mpi_reactive_strang
  public :: advance_mpi_reactive_strang_adaptive

contains

  subroutine mpi_reactive_timestep( &
      domain, species, transport, state, temperature, dx, hydro_cfl, &
      transport_cfl, transport_enabled, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, dt, &
      maximum_wave_speed, maximum_diffusivity, ok)
    type(mpi_domain_1d), intent(in) :: domain
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, 0:), temperature(0:)
    real(dp), intent(in) :: dx, hydro_cfl, transport_cfl
    logical, intent(in) :: transport_enabled, viscosity_enabled
    logical, intent(in) :: thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled
    real(dp), intent(out) :: dt, maximum_wave_speed, maximum_diffusivity
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: local_temperature, sound_speed
    real(dp) :: local_maximum_speed, hydro_dt, transport_dt
    logical :: local_ok, cell_ok, global_ok, reduction_ok
    integer :: i, ierr

    dt = 0.0_dp
    maximum_wave_speed = 0.0_dp
    maximum_diffusivity = 0.0_dp
    local_ok = valid_local_shape(domain, species, state, temperature) .and. &
      dx > 0.0_dp .and. hydro_cfl > 0.0_dp .and. hydro_cfl <= 1.0_dp
    if (transport_enabled) then
      local_ok = local_ok .and. transport_cfl > 0.0_dp .and. &
        transport_cfl <= 0.5_dp
    end if
    call collective_logical_and( &
      domain, local_ok, global_ok, reduction_ok)
    if (.not. reduction_ok .or. .not. global_ok) then
      ok = .false.
      return
    end if

    allocate(primitive(reactive_nprim(size(species))))
    local_maximum_speed = 0.0_dp
    local_ok = .true.
    do i = 1, domain%local_cells
      call reactive_conserved_to_primitive( &
        species, state(:, i), temperature(i), primitive, &
        local_temperature, sound_speed, cell_ok)
      if (cell_ok) then
        local_maximum_speed = max( &
          local_maximum_speed, abs(primitive(2)) + sound_speed)
      else
        local_ok = .false.
      end if
    end do
    call collective_logical_and( &
      domain, local_ok, global_ok, reduction_ok)
    if (.not. reduction_ok .or. .not. global_ok) then
      ok = .false.
      return
    end if

    call MPI_Allreduce( &
      local_maximum_speed, maximum_wave_speed, 1, MPI_DOUBLE_PRECISION, &
      MPI_MAX, domain%comm, ierr)
    if (ierr /= MPI_SUCCESS .or. maximum_wave_speed <= 0.0_dp) then
      ok = .false.
      return
    end if
    hydro_dt = hydro_cfl * dx / maximum_wave_speed

    if (transport_enabled) then
      call mpi_reactive_transport_timestep( &
        domain, species, transport, state, temperature, dx, transport_cfl, &
        viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, transport_dt, maximum_diffusivity, &
        local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
    else
      transport_dt = huge(1.0_dp)
    end if

    dt = min(hydro_dt, transport_dt)
    ok = ieee_is_finite(dt) .and. dt > 0.0_dp .and. &
      ieee_is_finite(maximum_wave_speed) .and. &
      ieee_is_finite(maximum_diffusivity)
  end subroutine mpi_reactive_timestep

  subroutine advance_mpi_reactive_hydro( &
      domain, species, state, temperature, dx, dt, riemann_solver, ok)
    type(mpi_domain_1d), intent(in) :: domain
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:, 0:), temperature(0:)
    real(dp), intent(in) :: dx, dt
    character(len=*), intent(in) :: riemann_solver
    logical, intent(out) :: ok

    real(dp), allocatable :: work_state(:, :), work_temperature(:)
    real(dp), allocatable :: output_state(:, :), output_temperature(:)
    real(dp), allocatable :: flux(:, :), primitive(:)
    real(dp) :: local_temperature, sound_speed
    logical :: local_ok, face_ok, global_ok, reduction_ok
    integer :: i, nvar

    local_ok = valid_local_shape(domain, species, state, temperature) .and. &
      dx > 0.0_dp .and. dt >= 0.0_dp
    call collective_logical_and( &
      domain, local_ok, global_ok, reduction_ok)
    if (.not. reduction_ok .or. .not. global_ok) then
      ok = .false.
      return
    end if
    if (dt <= 0.0_dp) then
      call exchange_state_temperature(domain, state, temperature, ok)
      return
    end if

    nvar = reactive_nvar(size(species))
    allocate(work_state(nvar, 0:domain%local_cells + 1))
    allocate(work_temperature(0:domain%local_cells + 1))
    allocate(output_state(nvar, 0:domain%local_cells + 1))
    allocate(output_temperature(0:domain%local_cells + 1))
    allocate(flux(nvar, 0:domain%local_cells))
    allocate(primitive(reactive_nprim(size(species))))
    work_state = state
    work_temperature = temperature
    call exchange_state_temperature( &
      domain, work_state, work_temperature, local_ok)
    call collective_logical_and( &
      domain, local_ok, global_ok, reduction_ok)
    if (.not. reduction_ok .or. .not. global_ok) then
      ok = .false.
      return
    end if

    local_ok = .true.
    do i = 0, domain%local_cells
      call reactive_riemann_flux_x( &
        species, work_state(:, i), work_state(:, i + 1), &
        work_temperature(i), work_temperature(i + 1), riemann_solver, &
        flux(:, i), face_ok)
      if (.not. face_ok) then
        local_ok = .false.
        flux(:, i) = 0.0_dp
      end if
    end do
    call collective_logical_and( &
      domain, local_ok, global_ok, reduction_ok)
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
    call collective_logical_and( &
      domain, local_ok, global_ok, reduction_ok)
    if (.not. reduction_ok .or. .not. global_ok) then
      ok = .false.
      return
    end if

    call exchange_state_temperature( &
      domain, output_state, output_temperature, local_ok)
    call collective_logical_and( &
      domain, local_ok, global_ok, reduction_ok)
    if (.not. reduction_ok .or. .not. global_ok) then
      ok = .false.
      return
    end if
    state = output_state
    temperature = output_temperature
    ok = .true.
  end subroutine advance_mpi_reactive_hydro

  subroutine advance_mpi_reactive_chemistry( &
      domain, species, reactions, state, temperature, interval, rtol, atol, &
      ok)
    type(mpi_domain_1d), intent(in) :: domain
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(inout) :: state(:, 0:), temperature(0:)
    real(dp), intent(in) :: interval, rtol, atol
    logical, intent(out) :: ok

    real(dp), allocatable :: work_state(:, :), work_temperature(:)
    logical :: local_ok, global_ok, reduction_ok
    integer :: nvar

    local_ok = valid_local_shape(domain, species, state, temperature) .and. &
      interval >= 0.0_dp .and. rtol > 0.0_dp .and. atol > 0.0_dp
    call collective_logical_and( &
      domain, local_ok, global_ok, reduction_ok)
    if (.not. reduction_ok .or. .not. global_ok) then
      ok = .false.
      return
    end if

    nvar = reactive_nvar(size(species))
    allocate(work_state(nvar, 0:domain%local_cells + 1))
    allocate(work_temperature(0:domain%local_cells + 1))
    work_state = state
    work_temperature = temperature
    call advance_reactive_chemistry( &
      species, reactions, work_state, work_temperature, &
      domain%local_cells, interval, rtol, atol, "periodic", local_ok)
    call collective_logical_and( &
      domain, local_ok, global_ok, reduction_ok)
    if (.not. reduction_ok .or. .not. global_ok) then
      ok = .false.
      return
    end if

    call exchange_state_temperature( &
      domain, work_state, work_temperature, local_ok)
    call collective_logical_and( &
      domain, local_ok, global_ok, reduction_ok)
    if (.not. reduction_ok .or. .not. global_ok) then
      ok = .false.
      return
    end if
    state = work_state
    temperature = work_temperature
    ok = .true.
  end subroutine advance_mpi_reactive_chemistry

  subroutine advance_mpi_reactive_strang( &
      domain, species, reactions, transport, state, temperature, dx, dt, &
      chemistry_enabled, rtol, atol, transport_enabled, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, riemann_solver, ok)
    type(mpi_domain_1d), intent(in) :: domain
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(inout) :: state(:, 0:), temperature(0:)
    real(dp), intent(in) :: dx, dt, rtol, atol
    logical, intent(in) :: chemistry_enabled, transport_enabled
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled
    logical, intent(in) :: barodiffusion_enabled
    character(len=*), intent(in) :: riemann_solver
    logical, intent(out) :: ok

    real(dp), allocatable :: work_state(:, :), work_temperature(:)
    logical :: local_ok, global_ok, reduction_ok
    integer :: nvar

    local_ok = valid_local_shape(domain, species, state, temperature) .and. &
      dx > 0.0_dp .and. dt >= 0.0_dp .and. rtol > 0.0_dp .and. &
      atol > 0.0_dp
    call collective_logical_and( &
      domain, local_ok, global_ok, reduction_ok)
    if (.not. reduction_ok .or. .not. global_ok) then
      ok = .false.
      return
    end if

    nvar = reactive_nvar(size(species))
    allocate(work_state(nvar, 0:domain%local_cells + 1))
    allocate(work_temperature(0:domain%local_cells + 1))
    work_state = state
    work_temperature = temperature

    if (chemistry_enabled) then
      call advance_mpi_reactive_chemistry( &
        domain, species, reactions, work_state, work_temperature, &
        0.5_dp * dt, rtol, atol, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
    end if
    if (transport_enabled) then
      call advance_mpi_reactive_transport( &
        domain, species, transport, work_state, work_temperature, dx, &
        0.5_dp * dt, viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
    end if

    call advance_mpi_reactive_hydro( &
      domain, species, work_state, work_temperature, dx, dt, &
      riemann_solver, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if

    if (transport_enabled) then
      call advance_mpi_reactive_transport( &
        domain, species, transport, work_state, work_temperature, dx, &
        0.5_dp * dt, viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
    end if
    if (chemistry_enabled) then
      call advance_mpi_reactive_chemistry( &
        domain, species, reactions, work_state, work_temperature, &
        0.5_dp * dt, rtol, atol, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
    end if

    state = work_state
    temperature = work_temperature
    ok = .true.
  end subroutine advance_mpi_reactive_strang

  subroutine advance_mpi_reactive_strang_adaptive( &
      domain, species, reactions, transport, state, temperature, dx, &
      requested_interval, minimum_interval, chemistry_enabled, rtol, atol, &
      transport_enabled, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, riemann_solver, &
      accepted_interval, rejected_trials, ok)
    type(mpi_domain_1d), intent(in) :: domain
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(inout) :: state(:, 0:), temperature(0:)
    real(dp), intent(in) :: dx, requested_interval, minimum_interval
    real(dp), intent(in) :: rtol, atol
    logical, intent(in) :: chemistry_enabled, transport_enabled
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled
    logical, intent(in) :: barodiffusion_enabled
    character(len=*), intent(in) :: riemann_solver
    real(dp), intent(out) :: accepted_interval
    integer, intent(out) :: rejected_trials
    logical, intent(out) :: ok

    real(dp), allocatable :: trial_state(:, :), trial_temperature(:)
    real(dp) :: trial_interval
    logical :: local_ok, global_ok, reduction_ok
    integer :: nvar

    accepted_interval = 0.0_dp
    rejected_trials = 0
    local_ok = valid_local_shape(domain, species, state, temperature) .and. &
      requested_interval >= 0.0_dp .and. minimum_interval > 0.0_dp .and. &
      dx > 0.0_dp .and. rtol > 0.0_dp .and. atol > 0.0_dp
    if (requested_interval > 0.0_dp) then
      local_ok = local_ok .and. minimum_interval <= requested_interval
    end if
    call collective_logical_and( &
      domain, local_ok, global_ok, reduction_ok)
    if (.not. reduction_ok .or. .not. global_ok) then
      ok = .false.
      return
    end if
    if (requested_interval <= 0.0_dp) then
      call exchange_state_temperature(domain, state, temperature, ok)
      return
    end if

    nvar = reactive_nvar(size(species))
    allocate(trial_state(nvar, 0:domain%local_cells + 1))
    allocate(trial_temperature(0:domain%local_cells + 1))
    trial_interval = requested_interval
    do
      trial_state = state
      trial_temperature = temperature
      call advance_mpi_reactive_strang( &
        domain, species, reactions, transport, trial_state, &
        trial_temperature, dx, trial_interval, chemistry_enabled, rtol, &
        atol, transport_enabled, viscosity_enabled, &
        thermal_conduction_enabled, species_diffusion_enabled, &
        barodiffusion_enabled, riemann_solver, local_ok)
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
  end subroutine advance_mpi_reactive_strang_adaptive

  subroutine exchange_state_temperature(domain, state, temperature, ok)
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
  end subroutine exchange_state_temperature

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

end module mpi_reactive_1d_mod
