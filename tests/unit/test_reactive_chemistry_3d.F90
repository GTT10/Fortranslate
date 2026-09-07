program test_reactive_chemistry_3d
  use, intrinsic :: ieee_arithmetic, only: &
    ieee_get_halting_mode, ieee_invalid, ieee_positive_inf, ieee_quiet_nan, &
    ieee_set_halting_mode, ieee_support_halting, ieee_value
  use precision_mod, only: dp
  use state_indices_mod, only: irho, ncons
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_mass_fraction_component, reactive_primitive_to_conserved, &
    advance_reactive_chemistry, resolve_reactive_chemistry_integrator
  use reactive_3d_mod, only: &
    advance_reactive_chemistry_3d, advance_reactive_euler_ssprk2_3d, &
    advance_reactive_euler_ssprk2_plm_3d, &
    advance_reactive_euler_ssprk2_plm_with_fluxes_3d, &
    advance_reactive_euler_ssprk2_with_fluxes_3d, advance_reactive_full_3d, &
    advance_reactive_strang_3d, compute_reactive_cfl_timestep_3d, &
    compute_reactive_euler_rhs_3d, reactive_integrals_3d, &
    recover_reactive_temperatures_3d
  implicit none

  type(nasa7_species), allocatable :: elementary_species(:), full_species(:)
  type(elementary_reaction), allocatable :: elementary_reactions(:)
  type(elementary_reaction), allocatable :: full_reactions(:)
  type(gas_transport_species), allocatable :: elementary_transport(:)
  logical :: ok

  call load_h2o2_elementary_thermo(elementary_species, ok)
  if (.not. ok) error stop "Failed to load elementary thermodynamics"
  call load_h2o2_elementary_mechanism(elementary_reactions, ok)
  if (.not. ok) error stop "Failed to load elementary mechanism"
  call load_h2o2_elementary_transport(elementary_transport, ok)
  if (.not. ok) error stop "Failed to load elementary transport"
  call load_h2o2_full_thermo(full_species, ok)
  if (.not. ok) error stop "Failed to load full thermodynamics"
  call load_h2o2_full_mechanism(full_reactions, ok)
  if (.not. ok) error stop "Failed to load full mechanism"

  call check_cell_local_reduction( &
    elementary_species, elementary_reactions, 1.0e-6_dp, ok)
  if (.not. ok) error stop "Elementary 3D chemistry reduction failed"
  call check_cell_local_reduction( &
    full_species, full_reactions, 2.0e-7_dp, ok)
  if (.not. ok) error stop "Full 3D chemistry reduction failed"
  call check_chemistry_integrator_metadata( &
    full_species, full_reactions, ok)
  if (.not. ok) error stop "Chemistry integrator metadata contract failed"
  call check_strang_reduction_and_rollback( &
    elementary_species, elementary_reactions, ok)
  if (.not. ok) error stop "Reactive 3D Strang contract failed"
  call check_full_split_contract( &
    elementary_species, elementary_reactions, elementary_transport, ok)
  if (.not. ok) error stop "Reactive 3D full split contract failed"
  call check_plm_full_split_uniform( &
    elementary_species, elementary_reactions, elementary_transport, ok)
  if (.not. ok) error stop "Reactive 3D PLM full split contract failed"
  call check_reactive_3d_finiteness( &
    elementary_species, elementary_reactions, elementary_transport, ok)
  if (.not. ok) error stop "Reactive 3D finiteness contract failed"
  write(*, '(a)') "test_reactive_chemistry_3d: PASS"

contains

  subroutine check_chemistry_integrator_metadata(species, reactions, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: cell_state(:), state(:, :), saved_state(:, :)
    real(dp), allocatable :: temperature(:), saved_temperature(:)
    real(dp) :: cell_temperature
    logical :: local_ok, use_implicit

    ok = .false.
    if (size(species) /= 10) return
    call resolve_reactive_chemistry_integrator( &
      size(species), use_implicit=use_implicit, ok=local_ok)
    if (.not. local_ok .or. .not. use_implicit) return
    call resolve_reactive_chemistry_integrator( &
      size(species), "explicit", use_implicit, local_ok)
    if (.not. local_ok .or. use_implicit) return
    call resolve_reactive_chemistry_integrator( &
      size(species), "implicit", use_implicit, local_ok)
    if (.not. local_ok .or. .not. use_implicit) return
    call resolve_reactive_chemistry_integrator( &
      size(species), "invalid", use_implicit, local_ok)
    if (local_ok) return

    call build_uniform_cell(species, cell_state, cell_temperature, local_ok)
    if (.not. local_ok) return
    allocate(state(size(cell_state), 0:2), saved_state(size(cell_state), 0:2))
    allocate(temperature(0:2), saved_temperature(0:2))
    state = spread(cell_state, 2, 3)
    saved_state = state
    temperature = cell_temperature
    saved_temperature = temperature
    call advance_reactive_chemistry( &
      species, reactions, state, temperature, 1, 2.0e-7_dp, 2.0e-7_dp, &
      1.0e-12_dp, "periodic", local_ok, "invalid")
    ok = .not. local_ok .and. maxval(abs(state - saved_state)) == 0.0_dp .and. &
      maxval(abs(temperature - saved_temperature)) == 0.0_dp
  end subroutine check_chemistry_integrator_metadata

  subroutine check_cell_local_reduction(species, reactions, interval, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: interval
    logical, intent(out) :: ok

    real(dp), allocatable :: cell_state(:), grid_state(:, :, :, :)
    real(dp), allocatable :: reference_state(:, :)
    real(dp), allocatable :: grid_temperature(:, :, :)
    real(dp), allocatable :: reference_temperature(:)
    real(dp) :: cell_temperature, state_error, temperature_error
    real(dp) :: initial_elements(3), final_elements(3), element_error
    real(dp) :: chemistry_change
    logical :: local_ok
    integer :: i, j, k, nvar

    call build_uniform_cell(species, cell_state, cell_temperature, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    nvar = size(cell_state)
    allocate(grid_state(nvar, 2, 2, 2), grid_temperature(2, 2, 2))
    do k = 1, 2
      do j = 1, 2
        do i = 1, 2
          grid_state(:, i, j, k) = cell_state
          grid_temperature(i, j, k) = cell_temperature
        end do
      end do
    end do
    allocate(reference_state(nvar, 0:2), reference_temperature(0:2))
    reference_state(:, 0) = cell_state
    reference_state(:, 1) = cell_state
    reference_state(:, 2) = cell_state
    reference_temperature = cell_temperature
    call elemental_amounts(species, cell_state, initial_elements, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if

    call advance_reactive_chemistry_3d( &
      species, reactions, grid_state, grid_temperature, 2, 2, 2, interval, &
      2.0e-7_dp, 1.0e-12_dp, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call advance_reactive_chemistry( &
      species, reactions, reference_state, reference_temperature, 1, &
      interval, 2.0e-7_dp, 1.0e-12_dp, "periodic", local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if

    state_error = 0.0_dp
    temperature_error = 0.0_dp
    do k = 1, 2
      do j = 1, 2
        do i = 1, 2
          state_error = max(state_error, maxval(abs( &
            grid_state(:, i, j, k) - reference_state(:, 1))))
          temperature_error = max(temperature_error, abs( &
            grid_temperature(i, j, k) - reference_temperature(1)))
        end do
      end do
    end do
    chemistry_change = maxval(abs( &
      grid_state(ncons + 1:nvar, 1, 1, 1) - cell_state(ncons + 1:nvar)))
    call elemental_amounts( &
      species, grid_state(:, 1, 1, 1), final_elements, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    element_error = maxval( &
      abs(final_elements - initial_elements) / max(1.0e-30_dp, &
        abs(initial_elements)))
    write(*, '(a,i0,4(a,es24.16))') &
      "nspecies=", size(species), ", state error=", state_error, &
      ", temperature error=", temperature_error, &
      ", chemistry change=", chemistry_change, &
      ", element error=", element_error
    ok = state_error <= 5.0e-13_dp .and. &
      temperature_error <= 5.0e-9_dp .and. chemistry_change > 1.0e-12_dp .and. &
      element_error <= 2.0e-10_dp .and. maxval(abs( &
        grid_state(1:ncons, 1, 1, 1) - cell_state(1:ncons))) == 0.0_dp
  end subroutine check_cell_local_reduction

  subroutine check_strang_reduction_and_rollback(species, reactions, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    logical, intent(out) :: ok

    real(dp), parameter :: dt = 1.0e-6_dp
    real(dp), allocatable :: cell_state(:), grid_state(:, :, :, :)
    real(dp), allocatable :: saved_state(:, :, :, :), reference_state(:, :)
    real(dp), allocatable :: grid_temperature(:, :, :)
    real(dp), allocatable :: saved_temperature(:, :, :)
    real(dp), allocatable :: reference_temperature(:)
    real(dp) :: cell_temperature, state_error, temperature_error
    logical :: local_ok
    integer :: i, j, k, nvar

    call build_uniform_cell(species, cell_state, cell_temperature, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    nvar = reactive_nvar(size(species))
    allocate(grid_state(nvar, 2, 2, 2), grid_temperature(2, 2, 2))
    do k = 1, 2
      do j = 1, 2
        do i = 1, 2
          grid_state(:, i, j, k) = cell_state
          grid_temperature(i, j, k) = cell_temperature
        end do
      end do
    end do
    allocate(reference_state(nvar, 0:2), reference_temperature(0:2))
    reference_state(:, 0) = cell_state
    reference_state(:, 1) = cell_state
    reference_state(:, 2) = cell_state
    reference_temperature = cell_temperature
    call advance_reactive_chemistry( &
      species, reactions, reference_state, reference_temperature, 1, &
      0.5_dp * dt, 2.0e-7_dp, 1.0e-12_dp, "periodic", local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call advance_reactive_chemistry( &
      species, reactions, reference_state, reference_temperature, 1, &
      0.5_dp * dt, 2.0e-7_dp, 1.0e-12_dp, "periodic", local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call advance_reactive_strang_3d( &
      species, reactions, grid_state, grid_temperature, 2, 2, 2, &
      0.01_dp, 0.01_dp, 0.01_dp, dt, "rusanov", .true., &
      2.0e-7_dp, 1.0e-12_dp, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    state_error = maxval(abs(grid_state(:, 1, 1, 1) - reference_state(:, 1)))
    temperature_error = abs(grid_temperature(1, 1, 1) - &
      reference_temperature(1))
    if (state_error > 5.0e-12_dp .or. temperature_error > 5.0e-8_dp) then
      ok = .false.
      return
    end if

    saved_state = grid_state
    saved_temperature = grid_temperature
    call advance_reactive_strang_3d( &
      species, reactions, grid_state, grid_temperature, 2, 2, 2, &
      0.01_dp, 0.01_dp, 0.01_dp, dt, "invalid", .true., &
      2.0e-7_dp, 1.0e-12_dp, local_ok)
    ok = .not. local_ok .and. &
      maxval(abs(grid_state - saved_state)) == 0.0_dp .and. &
      maxval(abs(grid_temperature - saved_temperature)) == 0.0_dp
  end subroutine check_strang_reduction_and_rollback

  subroutine check_full_split_contract(species, reactions, transport, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    logical, intent(out) :: ok

    real(dp), parameter :: dt = 1.0e-6_dp
    real(dp), allocatable :: cell_state(:)
    real(dp), allocatable :: reference_state(:, :, :, :)
    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: saved_state(:, :, :, :)
    real(dp), allocatable :: reference_temperature(:, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    real(dp), allocatable :: saved_temperature(:, :, :)
    real(dp) :: cell_temperature, state_error, temperature_error, theta
    logical :: local_ok
    integer :: i, j, k, nvar

    ok = .false.
    call build_uniform_cell(species, cell_state, cell_temperature, local_ok)
    if (.not. local_ok) return
    nvar = size(cell_state)
    allocate(reference_state(nvar, 2, 2, 2))
    allocate(candidate_state(nvar, 2, 2, 2))
    allocate(reference_temperature(2, 2, 2))
    allocate(candidate_temperature(2, 2, 2))
    do k = 1, 2
      do j = 1, 2
        do i = 1, 2
          reference_state(:, i, j, k) = cell_state
          reference_temperature(i, j, k) = cell_temperature
        end do
      end do
    end do
    candidate_state = reference_state
    candidate_temperature = reference_temperature

    call advance_reactive_strang_3d( &
      species, reactions, reference_state, reference_temperature, 2, 2, 2, &
      0.01_dp, 0.01_dp, 0.01_dp, dt, "rusanov", .true., &
      2.0e-7_dp, 1.0e-12_dp, local_ok)
    if (.not. local_ok) return
    call advance_reactive_full_3d( &
      species, reactions, transport, candidate_state, candidate_temperature, &
      2, 2, 2, 0.01_dp, 0.01_dp, 0.01_dp, dt, "rusanov", .true., &
      2.0e-7_dp, 1.0e-12_dp, .false., .true., .true., .true., .true., &
      theta, local_ok)
    if (.not. local_ok) return
    state_error = maxval(abs(candidate_state - reference_state))
    temperature_error = maxval(abs( &
      candidate_temperature - reference_temperature))
    if (state_error /= 0.0_dp .or. temperature_error /= 0.0_dp .or. &
        theta /= 1.0_dp) return

    call advance_reactive_strang_3d( &
      species, reactions, reference_state, reference_temperature, 2, 2, 2, &
      0.01_dp, 0.01_dp, 0.01_dp, dt, "rusanov", .false., &
      2.0e-7_dp, 1.0e-12_dp, local_ok)
    if (.not. local_ok) return
    call advance_reactive_full_3d( &
      species, reactions, transport, candidate_state, candidate_temperature, &
      2, 2, 2, 0.01_dp, 0.01_dp, 0.01_dp, dt, "rusanov", .false., &
      2.0e-7_dp, 1.0e-12_dp, .true., .true., .true., .true., .true., &
      theta, local_ok)
    if (.not. local_ok) return
    state_error = maxval(abs(candidate_state - reference_state))
    temperature_error = maxval(abs( &
      candidate_temperature - reference_temperature))
    if (state_error > 5.0e-13_dp .or. temperature_error > 5.0e-9_dp .or. &
        theta /= 1.0_dp) return

    saved_state = candidate_state
    saved_temperature = candidate_temperature
    call advance_reactive_full_3d( &
      species, reactions, transport, candidate_state, candidate_temperature, &
      2, 2, 2, 0.01_dp, 0.01_dp, 0.01_dp, dt, "invalid", .false., &
      2.0e-7_dp, 1.0e-12_dp, .true., .true., .true., .true., .true., &
      theta, local_ok)
    ok = .not. local_ok .and. &
      maxval(abs(candidate_state - saved_state)) == 0.0_dp .and. &
      maxval(abs(candidate_temperature - saved_temperature)) == 0.0_dp
    write(*, '(3(a,es24.16))') &
      "full split state error=", state_error, &
      ", temperature error=", temperature_error, ", theta=", theta
  end subroutine check_full_split_contract

  subroutine check_plm_full_split_uniform(species, reactions, transport, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    logical, intent(out) :: ok

    real(dp), parameter :: dt = 1.0e-6_dp
    real(dp), allocatable :: cell_state(:)
    real(dp), allocatable :: reference_state(:, :, :, :)
    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: saved_state(:, :, :, :)
    real(dp), allocatable :: reference_temperature(:, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    real(dp), allocatable :: saved_temperature(:, :, :)
    real(dp) :: cell_temperature, state_error, temperature_error, theta
    logical :: local_ok
    integer :: i, j, k, nvar

    ok = .false.
    call build_uniform_cell(species, cell_state, cell_temperature, local_ok)
    if (.not. local_ok) return
    nvar = size(cell_state)
    allocate(reference_state(nvar, 3, 3, 3))
    allocate(candidate_state(nvar, 3, 3, 3))
    allocate(reference_temperature(3, 3, 3))
    allocate(candidate_temperature(3, 3, 3))
    do k = 1, 3
      do j = 1, 3
        do i = 1, 3
          reference_state(:, i, j, k) = cell_state
          reference_temperature(i, j, k) = cell_temperature
        end do
      end do
    end do
    candidate_state = reference_state
    candidate_temperature = reference_temperature
    call advance_reactive_strang_3d( &
      species, reactions, reference_state, reference_temperature, 3, 3, 3, &
      0.01_dp, 0.01_dp, 0.01_dp, dt, "rusanov", .true., &
      2.0e-7_dp, 1.0e-12_dp, local_ok)
    if (.not. local_ok) return
    call advance_reactive_full_3d( &
      species, reactions, transport, candidate_state, candidate_temperature, &
      3, 3, 3, 0.01_dp, 0.01_dp, 0.01_dp, dt, "rusanov", .true., &
      2.0e-7_dp, 1.0e-12_dp, .true., .true., .true., .true., .true., &
      theta, local_ok, "characteristic_plm", "minmod")
    if (.not. local_ok) return
    state_error = maxval(abs(candidate_state - reference_state))
    temperature_error = maxval(abs( &
      candidate_temperature - reference_temperature))
    if (state_error > 5.0e-12_dp .or. temperature_error > 5.0e-8_dp .or. &
        theta /= 1.0_dp) return

    saved_state = candidate_state
    saved_temperature = candidate_temperature
    call advance_reactive_full_3d( &
      species, reactions, transport, candidate_state, candidate_temperature, &
      3, 3, 3, 0.01_dp, 0.01_dp, 0.01_dp, dt, "rusanov", .true., &
      2.0e-7_dp, 1.0e-12_dp, .true., .true., .true., .true., .true., &
      theta, local_ok, "characteristic_plm", "invalid")
    ok = .not. local_ok .and. &
      maxval(abs(candidate_state - saved_state)) == 0.0_dp .and. &
      maxval(abs(candidate_temperature - saved_temperature)) == 0.0_dp .and. &
      theta == 1.0_dp
    write(*, '(3(a,es24.16))') &
      "PLM full split state error=", state_error, &
      ", temperature error=", temperature_error, ", theta=", theta
  end subroutine check_plm_full_split_uniform

  subroutine check_reactive_3d_finiteness(species, reactions, transport, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    logical, intent(out) :: ok

    integer, parameter :: nx = 2, ny = 2, nz = 2
    integer, parameter :: plm_nx = 3, plm_ny = 3, plm_nz = 3
    real(dp), parameter :: dx = 0.01_dp, dy = 0.01_dp, dz = 0.01_dp
    real(dp), parameter :: valid_dt = 1.0e-8_dp
    real(dp), parameter :: valid_interval = 1.0e-8_dp
    real(dp), parameter :: valid_rtol = 2.0e-7_dp
    real(dp), parameter :: valid_atol = 1.0e-12_dp
    real(dp), allocatable :: cell_state(:)
    real(dp), allocatable :: state(:, :, :, :), saved_state(:, :, :, :)
    real(dp), allocatable :: plm_state(:, :, :, :), saved_plm_state(:, :, :, :)
    real(dp), allocatable :: temperature(:, :, :), saved_temperature(:, :, :)
    real(dp), allocatable :: plm_temperature(:, :, :), saved_plm_temperature(:, :, :)
    real(dp), allocatable :: rhs(:, :, :, :)
    real(dp), allocatable :: flux_x(:, :, :, :), flux_y(:, :, :, :)
    real(dp), allocatable :: flux_z(:, :, :, :)
    real(dp), allocatable :: plm_flux_x(:, :, :, :), plm_flux_y(:, :, :, :)
    real(dp), allocatable :: plm_flux_z(:, :, :, :)
    real(dp), allocatable :: recovered_temperature(:, :, :), integrals(:)
    real(dp) :: cell_temperature, bad_values(2), bad_value
    real(dp) :: dt, dx_test, dy_test, dz_test, cfl_test
    real(dp) :: interval_test, rtol_test, atol_test, theta
    logical :: local_ok, halting_supported, saved_halting
    integer :: i, j, k, nvar, bad_index, spacing_index, tolerance_index

    ok = .false.
    call build_uniform_cell(species, cell_state, cell_temperature, local_ok)
    if (.not. local_ok) return
    nvar = size(cell_state)
    allocate(state(nvar, nx, ny, nz), saved_state(nvar, nx, ny, nz))
    allocate(temperature(nx, ny, nz), saved_temperature(nx, ny, nz))
    allocate(plm_state(nvar, plm_nx, plm_ny, plm_nz))
    allocate(saved_plm_state(nvar, plm_nx, plm_ny, plm_nz))
    allocate(plm_temperature(plm_nx, plm_ny, plm_nz))
    allocate(saved_plm_temperature(plm_nx, plm_ny, plm_nz))
    allocate(rhs(nvar, nx, ny, nz))
    allocate(flux_x(nvar, nx, ny, nz), flux_y(nvar, nx, ny, nz))
    allocate(flux_z(nvar, nx, ny, nz))
    allocate(plm_flux_x(nvar, plm_nx, plm_ny, plm_nz))
    allocate(plm_flux_y(nvar, plm_nx, plm_ny, plm_nz))
    allocate(plm_flux_z(nvar, plm_nx, plm_ny, plm_nz))
    allocate(recovered_temperature(nx, ny, nz), integrals(nvar))
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          state(:, i, j, k) = cell_state
          temperature(i, j, k) = cell_temperature
        end do
      end do
    end do
    do k = 1, plm_nz
      do j = 1, plm_ny
        do i = 1, plm_nx
          plm_state(:, i, j, k) = cell_state
          plm_temperature(i, j, k) = cell_temperature
        end do
      end do
    end do
    saved_state = state
    saved_temperature = temperature
    saved_plm_state = plm_state
    saved_plm_temperature = plm_temperature
    bad_values(1) = ieee_value(0.0_dp, ieee_quiet_nan)
    bad_values(2) = ieee_value(0.0_dp, ieee_positive_inf)
    halting_supported = ieee_support_halting(ieee_invalid)
    saved_halting = .false.
    if (halting_supported) then
      call ieee_get_halting_mode(ieee_invalid, saved_halting)
    end if

    do bad_index = 1, 2
      cfl_test = bad_values(bad_index)
      dt = -1.0_dp
      if (bad_index == 1) then
        if (halting_supported) call ieee_set_halting_mode(ieee_invalid, .true.)
      end if
      call compute_reactive_cfl_timestep_3d( &
        species, state, temperature, nx, ny, nz, dx, dy, dz, cfl_test, dt, &
        local_ok)
      if (bad_index == 1) then
        if (halting_supported) call ieee_set_halting_mode(ieee_invalid, saved_halting)
      end if
      if (local_ok .or. dt /= 0.0_dp) return
    end do

    do spacing_index = 1, 3
      do bad_index = 1, 2
        dx_test = dx
        dy_test = dy
        dz_test = dz
        bad_value = bad_values(bad_index)
        select case (spacing_index)
        case (1)
          dx_test = bad_value
        case (2)
          dy_test = bad_value
        case (3)
          dz_test = bad_value
        end select
        rhs = 7.0_dp
        call compute_reactive_euler_rhs_3d( &
          species, state, temperature, nx, ny, nz, dx_test, dy_test, dz_test, &
          "rusanov", rhs, local_ok)
        if (local_ok .or. .not. all(rhs == 0.0_dp)) return
      end do
    end do

    do bad_index = 1, 2
      bad_value = bad_values(bad_index)
      state = saved_state
      temperature = saved_temperature
      call advance_reactive_euler_ssprk2_3d( &
        species, state, temperature, nx, ny, nz, dx, dy, dz, bad_value, &
        "rusanov", local_ok)
      if (local_ok .or. .not. all(state == saved_state) .or. &
          .not. all(temperature == saved_temperature)) return

      state = saved_state
      temperature = saved_temperature
      flux_x = 11.0_dp
      flux_y = 13.0_dp
      flux_z = 17.0_dp
      call advance_reactive_euler_ssprk2_with_fluxes_3d( &
        species, state, temperature, nx, ny, nz, dx, dy, dz, bad_value, &
        "rusanov", flux_x, flux_y, flux_z, local_ok)
      if (local_ok .or. .not. all(state == saved_state) .or. &
          .not. all(temperature == saved_temperature) .or. &
          .not. all(flux_x == 0.0_dp) .or. .not. all(flux_y == 0.0_dp) .or. &
          .not. all(flux_z == 0.0_dp)) return

      plm_state = saved_plm_state
      plm_temperature = saved_plm_temperature
      call advance_reactive_euler_ssprk2_plm_3d( &
        species, plm_state, plm_temperature, plm_nx, plm_ny, plm_nz, dx, dy, dz, &
        bad_value, "minmod", "rusanov", local_ok)
      if (local_ok .or. .not. all(plm_state == saved_plm_state) .or. &
          .not. all(plm_temperature == saved_plm_temperature)) return

      plm_state = saved_plm_state
      plm_temperature = saved_plm_temperature
      plm_flux_x = 19.0_dp
      plm_flux_y = 23.0_dp
      plm_flux_z = 29.0_dp
      call advance_reactive_euler_ssprk2_plm_with_fluxes_3d( &
        species, plm_state, plm_temperature, plm_nx, plm_ny, plm_nz, dx, dy, dz, &
        bad_value, "minmod", "rusanov", plm_flux_x, plm_flux_y, plm_flux_z, &
        local_ok)
      if (local_ok .or. .not. all(plm_state == saved_plm_state) .or. &
          .not. all(plm_temperature == saved_plm_temperature) .or. &
          .not. all(plm_flux_x == 0.0_dp) .or. &
          .not. all(plm_flux_y == 0.0_dp) .or. &
          .not. all(plm_flux_z == 0.0_dp)) return
    end do

    do tolerance_index = 1, 3
      do bad_index = 1, 2
        interval_test = valid_interval
        rtol_test = valid_rtol
        atol_test = valid_atol
        bad_value = bad_values(bad_index)
        select case (tolerance_index)
        case (1)
          interval_test = bad_value
        case (2)
          rtol_test = bad_value
        case (3)
          atol_test = bad_value
        end select
        state = saved_state
        temperature = saved_temperature
        call advance_reactive_chemistry_3d( &
          species, reactions, state, temperature, nx, ny, nz, interval_test, &
          rtol_test, atol_test, local_ok)
        if (local_ok .or. .not. all(state == saved_state) .or. &
            .not. all(temperature == saved_temperature)) return
      end do
    end do

    do bad_index = 1, 2
      bad_value = bad_values(bad_index)
      state = saved_state
      temperature = saved_temperature
      call advance_reactive_strang_3d( &
        species, reactions, state, temperature, nx, ny, nz, dx, dy, dz, &
        bad_value, "rusanov", .true., valid_rtol, valid_atol, local_ok)
      if (local_ok .or. .not. all(state == saved_state) .or. &
          .not. all(temperature == saved_temperature)) return

      state = saved_state
      temperature = saved_temperature
      call advance_reactive_full_3d( &
        species, reactions, transport, state, temperature, nx, ny, nz, &
        dx, dy, dz, bad_value, "rusanov", .true., valid_rtol, valid_atol, &
        .false., .false., .false., .false., .false., theta, local_ok)
      if (local_ok .or. .not. all(state == saved_state) .or. &
          .not. all(temperature == saved_temperature) .or. theta /= 1.0_dp) return
    end do

    do tolerance_index = 2, 3
      do bad_index = 1, 2
        rtol_test = valid_rtol
        atol_test = valid_atol
        bad_value = bad_values(bad_index)
        if (tolerance_index == 2) then
          rtol_test = bad_value
        else
          atol_test = bad_value
        end if
        state = saved_state
        temperature = saved_temperature
        call advance_reactive_strang_3d( &
          species, reactions, state, temperature, nx, ny, nz, dx, dy, dz, &
          valid_dt, "rusanov", .true., rtol_test, atol_test, local_ok)
        if (local_ok .or. .not. all(state == saved_state) .or. &
            .not. all(temperature == saved_temperature)) return

        state = saved_state
        temperature = saved_temperature
        call advance_reactive_full_3d( &
          species, reactions, transport, state, temperature, nx, ny, nz, &
          dx, dy, dz, valid_dt, "rusanov", .true., rtol_test, atol_test, &
          .false., .false., .false., .false., .false., theta, local_ok)
        if (local_ok .or. .not. all(state == saved_state) .or. &
            .not. all(temperature == saved_temperature) .or. theta /= 1.0_dp) return
      end do
    end do

    do bad_index = 1, 2
      bad_value = bad_values(bad_index)
      state = saved_state
      state(irho, 1, 1, 1) = bad_value
      recovered_temperature = 42.0_dp
      call recover_reactive_temperatures_3d( &
        species, state, saved_temperature, nx, ny, nz, recovered_temperature, &
        local_ok)
      if (local_ok .or. .not. all(recovered_temperature == 0.0_dp)) return
    end do

    do spacing_index = 1, 3
      do bad_index = 1, 2
        dx_test = dx
        dy_test = dy
        dz_test = dz
        bad_value = bad_values(bad_index)
        select case (spacing_index)
        case (1)
          dx_test = bad_value
        case (2)
          dy_test = bad_value
        case (3)
          dz_test = bad_value
        end select
        integrals = 37.0_dp
        call reactive_integrals_3d( &
          saved_state, nx, ny, nz, dx_test, dy_test, dz_test, integrals, local_ok)
        if (local_ok .or. .not. all(integrals == 0.0_dp)) return
      end do
    end do
    do bad_index = 1, 2
      state = saved_state
      state(irho, 1, 1, 1) = bad_values(bad_index)
      integrals = 41.0_dp
      call reactive_integrals_3d( &
        state, nx, ny, nz, dx, dy, dz, integrals, local_ok)
      if (local_ok .or. .not. all(integrals == 0.0_dp)) return
    end do

    ok = .true.
    write(*, '(a)') "reactive 3D finite-input rollback checks: PASS"
  end subroutine check_reactive_3d_finiteness

  subroutine build_uniform_cell(species, state, temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), allocatable, intent(out) :: state(:)
    real(dp), intent(out) :: temperature
    logical, intent(out) :: ok

    real(dp), allocatable :: mole_fractions(:), mass_fractions(:), primitive(:)
    real(dp) :: density, sound_speed
    logical :: local_ok
    integer :: species_index, nspecies

    nspecies = size(species)
    allocate(mole_fractions(nspecies), mass_fractions(nspecies))
    allocate(primitive(reactive_nprim(nspecies)), state(reactive_nvar(nspecies)))
    if (nspecies == 7) then
      mole_fractions = [ &
        0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
        1.0e-5_dp, 0.0_dp, 0.55643_dp]
    else if (nspecies == 10) then
      mole_fractions = [ &
        0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
        1.0e-5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.55643_dp]
    else
      ok = .false.
      return
    end if
    call mass_fractions_from_mole_fractions( &
      species, mole_fractions, mass_fractions, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    density = mixture_density( &
      species, mass_fractions, 101325.0_dp, 1200.0_dp, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    primitive(1:5) = [density, 0.0_dp, 0.0_dp, 0.0_dp, 101325.0_dp]
    do species_index = 1, nspecies
      primitive(reactive_mass_fraction_component(species_index)) = &
        mass_fractions(species_index)
    end do
    call reactive_primitive_to_conserved( &
      species, primitive, state, temperature, sound_speed, ok)
  end subroutine build_uniform_cell

  subroutine elemental_amounts(species, state, elements, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:)
    real(dp), intent(out) :: elements(3)
    logical, intent(out) :: ok

    real(dp) :: atoms(3), amount
    integer :: species_index

    elements = 0.0_dp
    ok = size(state) == reactive_nvar(size(species))
    if (.not. ok) return
    do species_index = 1, size(species)
      atoms = 0.0_dp
      select case (trim(species(species_index)%name))
      case ("H2")
        atoms = [2.0_dp, 0.0_dp, 0.0_dp]
      case ("H")
        atoms = [1.0_dp, 0.0_dp, 0.0_dp]
      case ("O")
        atoms = [0.0_dp, 1.0_dp, 0.0_dp]
      case ("O2")
        atoms = [0.0_dp, 2.0_dp, 0.0_dp]
      case ("OH")
        atoms = [1.0_dp, 1.0_dp, 0.0_dp]
      case ("H2O")
        atoms = [2.0_dp, 1.0_dp, 0.0_dp]
      case ("HO2")
        atoms = [1.0_dp, 2.0_dp, 0.0_dp]
      case ("H2O2")
        atoms = [2.0_dp, 2.0_dp, 0.0_dp]
      case ("N2")
        atoms = [0.0_dp, 0.0_dp, 2.0_dp]
      case ("AR")
        cycle
      case default
        ok = .false.
        return
      end select
      amount = state(reactive_species_component(species_index)) / &
        species(species_index)%molecular_weight
      elements = elements + atoms * amount
    end do
  end subroutine elemental_amounts

end program test_reactive_chemistry_3d
