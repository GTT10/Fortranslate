program test_mpi_amr_reactive_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, &
    ieee_value
  use, intrinsic :: iso_fortran_env, only: int64
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use gas_transport_mod, only: gas_transport_species
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use transport_database_mod, only: load_h2o2_elementary_transport
  use mesh_3d_mod, only: uniform_cell_centers_3d
  use reactive_1d_mod, only: reactive_nvar, reactive_species_component
  use state_indices_mod, only: irho
  use simulation_config_reactive_3d_mod, only: reactive_3d_config
  use reactive_entropy_wave_3d_problem_mod, only: &
    initialize_reactive_problem_3d
  use reactive_3d_mod, only: &
    recover_reactive_temperatures_3d, &
    advance_reactive_euler_ssprk2_plm_with_fluxes_3d
  use amr_hierarchy_3d_mod, only: &
    amr_patch_3d, initialize_amr_patch_3d, average_down_3d, &
    composite_integrals_amr_3d
  use amr_reactive_3d_mod, only: &
    compute_amr_reactive_cfl_timestep_3d, advance_amr_reactive_hydro_3d
  use amr_reactive_transport_3d_mod, only: &
    compute_amr_reactive_transport_timestep_3d, &
    advance_amr_reactive_transport_3d, advance_amr_reactive_full_3d
  use mpi_amr_reactive_3d_mod, only: &
    mpi_amr_slab_distribution_3d, &
    initialize_mpi_amr_slab_distribution_3d, &
    compute_mpi_amr_reactive_cfl_timestep_3d, &
    advance_mpi_amr_reactive_hydro_3d, &
    broadcast_mpi_amr_reactive_hierarchy_3d
  use mpi_amr_sparse_reactive_3d_mod, only: &
    mpi_amr_sparse_distribution_3d, mpi_amr_sparse_hierarchy_3d, &
    initialize_mpi_amr_sparse_distribution_3d, &
    scatter_mpi_amr_sparse_hierarchy_3d, &
    gather_mpi_amr_sparse_hierarchy_3d, &
    compute_mpi_amr_sparse_cfl_timestep_3d, &
    build_mpi_amr_sparse_periodic_halos_3d, &
    advance_mpi_amr_sparse_periodic_hydro_3d, &
    advance_mpi_amr_sparse_hydro_3d, &
    compute_mpi_amr_sparse_transport_timestep_3d, &
    advance_mpi_amr_sparse_transport_3d, advance_mpi_amr_sparse_full_3d
  implicit none

  integer, parameter :: n = 8, ratio = 2
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(gas_transport_species), allocatable :: mismatched_transport(:)
  type(nasa7_species), allocatable :: mismatched_species(:)
  type(reactive_3d_config) :: config, fine_config
  type(amr_patch_3d) :: patch, mismatched_patch
  type(mpi_amr_slab_distribution_3d) :: distribution
  type(mpi_amr_sparse_distribution_3d) :: sparse_distribution
  type(mpi_amr_sparse_hierarchy_3d) :: sparse_hierarchy
  type(mpi_amr_sparse_hierarchy_3d) :: saved_sparse_hierarchy
  type(mpi_amr_sparse_hierarchy_3d) :: rejected_sparse_hierarchy
  real(dp), allocatable :: coarse_state(:, :, :, :)
  real(dp), allocatable :: fine_state(:, :, :, :)
  real(dp), allocatable :: coarse_temperature(:, :, :)
  real(dp), allocatable :: fine_temperature(:, :, :)
  real(dp), allocatable :: recovered_temperature(:, :, :)
  real(dp), allocatable :: serial_coarse(:, :, :, :)
  real(dp), allocatable :: serial_fine(:, :, :, :)
  real(dp), allocatable :: serial_coarse_temperature(:, :, :)
  real(dp), allocatable :: serial_fine_temperature(:, :, :)
  real(dp), allocatable :: saved_coarse(:, :, :, :)
  real(dp), allocatable :: saved_fine(:, :, :, :)
  real(dp), allocatable :: saved_coarse_temperature(:, :, :)
  real(dp), allocatable :: saved_fine_temperature(:, :, :)
  real(dp), allocatable :: saved_integrals(:)
  real(dp), allocatable :: gathered_coarse(:, :, :, :)
  real(dp), allocatable :: gathered_fine(:, :, :, :)
  real(dp), allocatable :: gathered_coarse_temperature(:, :, :)
  real(dp), allocatable :: gathered_fine_temperature(:, :, :)
  real(dp), allocatable :: sparse_halo_state(:, :, :, :)
  real(dp), allocatable :: sparse_halo_temperature(:, :, :)
  real(dp), allocatable :: sparse_flux_x(:, :, :, :)
  real(dp), allocatable :: sparse_flux_y(:, :, :, :)
  real(dp), allocatable :: sparse_flux_z(:, :, :, :)
  real(dp), allocatable :: serial_flux_x(:, :, :, :)
  real(dp), allocatable :: serial_flux_y(:, :, :, :)
  real(dp), allocatable :: serial_flux_z(:, :, :, :)
  real(dp), allocatable :: mismatched_fine(:, :, :, :)
  real(dp), allocatable :: mismatched_fine_temperature(:, :, :)
  real(dp), allocatable :: saved_mismatched_fine(:, :, :, :)
  real(dp), allocatable :: saved_mismatched_fine_temperature(:, :, :)
  real(dp), allocatable :: initial_integrals(:), mass_fractions(:)
  real(dp), allocatable :: x(:), y(:), z(:), xf(:), yf(:), zf(:)
  real(dp) :: dx, dy, dz, density, fine_density, serial_dt, distributed_dt
  real(dp) :: sparse_dt
  real(dp) :: serial_transport_dt, sparse_transport_dt
  real(dp) :: serial_maximum_diffusivity, sparse_maximum_diffusivity
  real(dp) :: serial_transport_theta, sparse_transport_theta
  real(dp) :: serial_transport_reflux, sparse_transport_reflux
  real(dp) :: serial_full_theta, sparse_full_theta
  real(dp) :: serial_full_reflux, sparse_full_reflux
  real(dp) :: serial_reflux, distributed_reflux, time, maximum_reflux
  real(dp) :: sparse_reflux
  real(dp) :: saved_time, saved_maximum_reflux
  logical :: ok
  integer :: ierr, nvar, steps, saved_steps, inconsistent_root
  integer :: mismatched_fine_nx, global_sparse_value_count
  integer :: candidate_rank, expected_previous_fine_rank
  integer :: expected_next_fine_rank
  character(len=32) :: sparse_reconstruction

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Init failed"
  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "elementary thermodynamics load")
  call load_h2o2_elementary_mechanism(reactions, ok)
  call require(ok, "elementary mechanism load")
  call load_h2o2_elementary_transport(transport, ok)
  call require(ok, "elementary transport load")
  call initialize_amr_patch_3d( &
    n, n, n, 3, 6, 3, 6, 3, 6, ratio, patch, ok)
  call require(ok .and. patch%is_strictly_interior(), "MPI AMR patch")
  call initialize_mpi_amr_slab_distribution_3d( &
    patch, MPI_COMM_WORLD, distribution, ok)
  call require(ok .and. distribution%is_valid(patch), &
    "MPI AMR slab distribution")

  config = reactive_3d_config()
  config%nx = n
  config%ny = n
  config%nz = n
  config%thermo_model = "elementary"
  config%riemann_solver = "pelec"
  config%initial_velocity_x = 240.0_dp
  config%initial_velocity_y = -90.0_dp
  config%initial_velocity_z = 60.0_dp
  config%wave_number_x = 1
  config%wave_number_y = 2
  config%wave_number_z = 1
  nvar = reactive_nvar(size(species))
  allocate(coarse_state(nvar, n, n, n), coarse_temperature(n, n, n))
  allocate(fine_state(nvar, patch%fine_nx(), patch%fine_ny(), &
    patch%fine_nz()))
  allocate(fine_temperature( &
    patch%fine_nx(), patch%fine_ny(), patch%fine_nz()))
  allocate(recovered_temperature(n, n, n))
  allocate(initial_integrals(nvar), mass_fractions(size(species)))
  allocate(x(n), y(n), z(n), xf(patch%fine_nx()), &
    yf(patch%fine_ny()), zf(patch%fine_nz()))
  call uniform_cell_centers_3d( &
    n, n, n, config%x_lower, config%x_upper, &
    config%y_lower, config%y_upper, config%z_lower, config%z_upper, &
    x, y, z, dx, dy, dz)
  call fine_patch_centers(patch, config, dx, dy, dz, xf, yf, zf)
  call initialize_reactive_problem_3d( &
    species, config, x, y, z, coarse_state, coarse_temperature, &
    density, mass_fractions, ok)
  call require(ok .and. density > 0.0_dp, "MPI coarse initialization")
  fine_config = config
  fine_config%nx = patch%fine_nx()
  fine_config%ny = patch%fine_ny()
  fine_config%nz = patch%fine_nz()
  call initialize_reactive_problem_3d( &
    species, fine_config, xf, yf, zf, fine_state, fine_temperature, &
    fine_density, mass_fractions, ok)
  call require(ok .and. fine_density > 0.0_dp, "MPI fine initialization")
  call average_down_3d(coarse_state, fine_state, patch, ok)
  call require(ok, "MPI initial average-down")
  call recover_reactive_temperatures_3d( &
    species, coarse_state, coarse_temperature, n, n, n, &
    recovered_temperature, ok)
  call require(ok, "MPI initial temperature recovery")
  coarse_temperature = recovered_temperature
  call composite_integrals_amr_3d( &
    coarse_state, fine_state, patch, dx, dy, dz, initial_integrals, ok)
  call require(ok, "MPI initial composite integrals")

  time = 0.0_dp
  steps = 0
  maximum_reflux = 0.0_dp
  if (distribution%rank /= 0) then
    coarse_state = -1.0_dp
    coarse_temperature = -1.0_dp
    fine_state = -1.0_dp
    fine_temperature = -1.0_dp
    initial_integrals = -1.0_dp
  end if
  call broadcast_mpi_amr_reactive_hierarchy_3d( &
    distribution, patch, 0, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, time, steps, initial_integrals, &
    maximum_reflux, ok)
  call require(ok, "transactional hierarchy broadcast")

  call initialize_mpi_amr_sparse_distribution_3d( &
    patch, MPI_COMM_WORLD, sparse_distribution, ok)
  call require(ok .and. sparse_distribution%is_valid(patch), &
    "rank-local sparse AMR distribution")
  expected_previous_fine_rank = MPI_PROC_NULL
  expected_next_fine_rank = MPI_PROC_NULL
  if (sparse_distribution%fine_count > 0) then
    do candidate_rank = sparse_distribution%rank - 1, 0, -1
      if (sparse_distribution%fine_counts(candidate_rank + 1) > 0) then
        expected_previous_fine_rank = candidate_rank
        exit
      end if
    end do
    do candidate_rank = sparse_distribution%rank + 1, &
        sparse_distribution%nranks - 1
      if (sparse_distribution%fine_counts(candidate_rank + 1) > 0) then
        expected_next_fine_rank = candidate_rank
        exit
      end if
    end do
  end if
  call require( &
    sparse_distribution%previous_fine_rank == &
      expected_previous_fine_rank .and. &
    sparse_distribution%next_fine_rank == expected_next_fine_rank, &
    "parent-aligned fine neighbor ranks")
  call require(sum(sparse_distribution%fine_counts) == patch%fine_nx(), &
    "parent-aligned fine ownership coverage")
  if (sparse_distribution%nranks > 1) then
    inconsistent_root = merge(0, 1, sparse_distribution%rank == 0)
    call scatter_mpi_amr_sparse_hierarchy_3d( &
      sparse_distribution, species, patch, inconsistent_root, coarse_state, &
      coarse_temperature, fine_state, fine_temperature, &
      rejected_sparse_hierarchy, ok)
    call require(.not. ok, "inconsistent sparse scatter root rejection")
  end if
  call scatter_mpi_amr_sparse_hierarchy_3d( &
    sparse_distribution, species, patch, 0, coarse_state, &
    coarse_temperature, fine_state, fine_temperature, sparse_hierarchy, ok)
  call require(ok .and. sparse_hierarchy%is_valid( &
      sparse_distribution, patch), "rank-local sparse hierarchy scatter")
  call require(size(sparse_hierarchy%coarse_state, 2) == &
      sparse_distribution%coarse_count .and. &
    size(sparse_hierarchy%fine_state, 2) == &
      sparse_distribution%fine_count, &
    "rank-local sparse storage extents")
  global_sparse_value_count = (nvar + 1) * ( &
    n * n * n + patch%fine_nx() * patch%fine_ny() * patch%fine_nz())
  if (sparse_distribution%nranks > 1) then
    call require(sparse_hierarchy%local_value_count() < &
      global_sparse_value_count, "rank-local sparse storage reduction")
  end if
  if (sparse_distribution%nranks == 8) then
    call require(count(sparse_distribution%fine_counts == 0) == 4, &
      "zero-fine-cell ranks remain valid")
  end if

  ! The sparse transport path must use the same elementary transport records,
  ! hierarchy geometry, and two-level state as the serial AMR reference.  Keep
  ! this gate before the existing hydro mutations so both references start from
  ! the identical initialized hierarchy on every rank.
  serial_coarse = coarse_state
  serial_coarse_temperature = coarse_temperature
  serial_fine = fine_state
  serial_fine_temperature = fine_temperature
  call compute_amr_reactive_transport_timestep_3d( &
    species, transport, patch, serial_coarse, serial_coarse_temperature, &
    serial_fine, serial_fine_temperature, dx, dy, dz, 0.35_dp, &
    .true., .true., .true., serial_transport_dt, &
    serial_maximum_diffusivity, ok)
  call require(ok, "serial sparse transport timestep reference")
  call compute_mpi_amr_sparse_transport_timestep_3d( &
    sparse_distribution, species, transport, patch, sparse_hierarchy, &
    dx, dy, dz, 0.35_dp, .true., .true., .true., sparse_transport_dt, &
    sparse_maximum_diffusivity, ok)
  call require(ok .and. same_real_bits( &
      serial_transport_dt, sparse_transport_dt) .and. same_real_bits( &
      serial_maximum_diffusivity, sparse_maximum_diffusivity), &
    "rank-local sparse transport timestep bit parity")

  serial_coarse = coarse_state
  serial_coarse_temperature = coarse_temperature
  serial_fine = fine_state
  serial_fine_temperature = fine_temperature
  call advance_amr_reactive_transport_3d( &
    species, transport, patch, serial_coarse, serial_coarse_temperature, &
    serial_fine, serial_fine_temperature, dx, dy, dz, 1.0e-8_dp, &
    .true., .true., .true., .true., "mc", serial_transport_theta, &
    serial_transport_reflux, ok)
  call require(ok, "serial sparse transport reference step")
  call advance_mpi_amr_sparse_transport_3d( &
    sparse_distribution, species, transport, patch, sparse_hierarchy, &
    dx, dy, dz, 1.0e-8_dp, .true., .true., .true., .true., "mc", &
    sparse_transport_theta, sparse_transport_reflux, ok)
  call require(ok, "rank-local sparse transport step")
  call gather_mpi_amr_sparse_hierarchy_3d( &
    sparse_distribution, patch, 0, sparse_hierarchy, gathered_coarse, &
    gathered_coarse_temperature, gathered_fine, gathered_fine_temperature, &
    ok)
  call require(ok, "rank-local sparse transport gather")
  call require(same_real_bits( &
      serial_transport_theta, sparse_transport_theta) .and. &
    same_real_bits(serial_transport_reflux, sparse_transport_reflux), &
    "rank-local sparse transport metadata bit parity")
  if (sparse_distribution%rank == 0) then
    call require(rank4_bits_match(gathered_coarse, serial_coarse) .and. &
      rank3_bits_match(gathered_coarse_temperature, &
        serial_coarse_temperature) .and. &
      rank4_bits_match(gathered_fine, serial_fine) .and. &
      rank3_bits_match(gathered_fine_temperature, serial_fine_temperature), &
      "rank-local sparse transport state bit parity")
  end if

  ! Invalid barodiffusion is a collective contract failure.  It must leave
  ! both sparse levels and both output metadata values bitwise untouched.
  saved_sparse_hierarchy = sparse_hierarchy
  call advance_mpi_amr_sparse_transport_3d( &
    sparse_distribution, species, transport, patch, sparse_hierarchy, &
    dx, dy, dz, 1.0e-8_dp, .true., .true., .false., .true., "mc", &
    sparse_transport_theta, sparse_transport_reflux, ok)
  call require(.not. ok .and. same_real_bits( &
      sparse_transport_theta, 1.0_dp) .and. same_real_bits( &
      sparse_transport_reflux, 0.0_dp) .and. sparse_hierarchies_bits_match( &
      sparse_hierarchy, saved_sparse_hierarchy), &
    "sparse invalid barodiffusion transaction")

  ! A valid transport record changed on one rank is still a collective
  ! contract violation.  Exercise both new public entry points and verify that
  ! neither output nor hierarchy is published after rejection.
  if (sparse_distribution%nranks > 1) then
    mismatched_transport = transport
    if (sparse_distribution%rank == sparse_distribution%nranks - 1) then
      mismatched_transport(1)%well_depth = nearest( &
        mismatched_transport(1)%well_depth, 1.0_dp)
    end if
    call compute_mpi_amr_sparse_transport_timestep_3d( &
      sparse_distribution, species, mismatched_transport, patch, &
      sparse_hierarchy, dx, dy, dz, 0.35_dp, .true., .true., .true., &
      sparse_transport_dt, sparse_maximum_diffusivity, ok)
    call require(.not. ok .and. same_real_bits( &
        sparse_transport_dt, 0.0_dp) .and. same_real_bits( &
        sparse_maximum_diffusivity, 0.0_dp), &
      "rank-dependent sparse transport metadata rejection")
    saved_sparse_hierarchy = sparse_hierarchy
    call advance_mpi_amr_sparse_transport_3d( &
      sparse_distribution, species, mismatched_transport, patch, &
      sparse_hierarchy, dx, dy, dz, 1.0e-8_dp, .true., .true., .true., &
      .true., "mc", sparse_transport_theta, sparse_transport_reflux, ok)
    call require(.not. ok .and. same_real_bits( &
        sparse_transport_theta, 1.0_dp) .and. same_real_bits( &
        sparse_transport_reflux, 0.0_dp) .and. &
      sparse_hierarchies_bits_match(sparse_hierarchy, saved_sparse_hierarchy), &
      "rank-dependent sparse transport step rollback")
  end if

  call run_adversarial_sparse_transport_tests()
  call run_sparse_full_transport_test()

  call scatter_mpi_amr_sparse_hierarchy_3d( &
    sparse_distribution, species, patch, 0, coarse_state, &
    coarse_temperature, fine_state, fine_temperature, sparse_hierarchy, ok)
  call require(ok, "restore sparse hierarchy after transport parity")
  call run_sparse_nonfinite_guard_tests()

  call compute_amr_reactive_cfl_timestep_3d( &
    species, patch, coarse_state, coarse_temperature, fine_state, &
    fine_temperature, dx, dy, dz, 0.30_dp, serial_dt, ok)
  call require(ok, "sparse CFL serial reference")
  call compute_mpi_amr_sparse_cfl_timestep_3d( &
    sparse_distribution, species, patch, sparse_hierarchy, &
    dx, dy, dz, 0.30_dp, sparse_dt, ok)
  call require(ok .and. same_real_bits(serial_dt, sparse_dt), &
    "rank-local sparse CFL bit parity")
  if (any(sparse_distribution%fine_counts == 0)) then
    rejected_sparse_hierarchy = sparse_hierarchy
    if (sparse_distribution%fine_count == 0) &
      rejected_sparse_hierarchy%coarse_state(1, 1, 1, 1) = -1.0_dp
    call compute_mpi_amr_sparse_cfl_timestep_3d( &
      sparse_distribution, species, patch, rejected_sparse_hierarchy, &
      dx, dy, dz, 0.30_dp, sparse_dt, ok)
    call require(.not. ok .and. same_real_bits(sparse_dt, 0.0_dp), &
      "zero-fine-rank invalid coarse CFL rejection")
  end if
  call build_mpi_amr_sparse_periodic_halos_3d( &
    sparse_distribution, sparse_hierarchy%coarse_state, &
    sparse_hierarchy%coarse_temperature, sparse_halo_state, &
    sparse_halo_temperature, ok)
  call require(ok .and. sparse_periodic_halos_match( &
      sparse_distribution, coarse_state, coarse_temperature, &
      sparse_halo_state, sparse_halo_temperature), &
    "two-plane sparse periodic halo bit parity")
  serial_coarse = coarse_state
  serial_coarse_temperature = coarse_temperature
  allocate(serial_flux_x, mold=serial_coarse)
  allocate(serial_flux_y, mold=serial_coarse)
  allocate(serial_flux_z, mold=serial_coarse)
  call advance_reactive_euler_ssprk2_plm_with_fluxes_3d( &
    species, serial_coarse, serial_coarse_temperature, n, n, n, &
    dx, dy, dz, 1.0e-7_dp, "mc", config%riemann_solver, &
    serial_flux_x, serial_flux_y, serial_flux_z, ok)
  call require(ok, "sparse periodic PLM serial reference")
  call advance_mpi_amr_sparse_periodic_hydro_3d( &
    sparse_distribution, species, patch, sparse_hierarchy%coarse_state, &
    sparse_hierarchy%coarse_temperature, dx, dy, dz, 1.0e-7_dp, &
    config%riemann_solver, "characteristic_plm", "mc", sparse_flux_x, &
    sparse_flux_y, sparse_flux_z, ok)
  call require(ok, "rank-local sparse periodic PLM step")
  call gather_mpi_amr_sparse_hierarchy_3d( &
    sparse_distribution, patch, 0, sparse_hierarchy, gathered_coarse, &
    gathered_coarse_temperature, gathered_fine, &
    gathered_fine_temperature, ok)
  call require(ok, "sparse periodic PLM gather")
  if (sparse_distribution%rank == 0) then
    call require(rank4_bits_match(gathered_coarse, serial_coarse) .and. &
      rank3_bits_match( &
        gathered_coarse_temperature, serial_coarse_temperature), &
      "rank-local sparse periodic PLM bit parity")
  end if
  call scatter_mpi_amr_sparse_hierarchy_3d( &
    sparse_distribution, species, patch, 0, coarse_state, &
    coarse_temperature, fine_state, fine_temperature, sparse_hierarchy, ok)
  call require(ok, "restore sparse hierarchy after periodic PLM test")
  serial_coarse = coarse_state
  serial_coarse_temperature = coarse_temperature
  serial_fine = fine_state
  serial_fine_temperature = fine_temperature
  call advance_amr_reactive_hydro_3d( &
    species, patch, serial_coarse, serial_coarse_temperature, &
    serial_fine, serial_fine_temperature, dx, dy, dz, 1.0e-7_dp, &
    config%riemann_solver, serial_reflux, ok, "characteristic_plm", "mc")
  call require(ok, "rank-local sparse PLM AMR serial reference")
  call advance_mpi_amr_sparse_hydro_3d( &
    sparse_distribution, species, patch, sparse_hierarchy, &
    dx, dy, dz, 1.0e-7_dp, config%riemann_solver, &
    "characteristic_plm", "mc", sparse_reflux, ok)
  call require(ok, "rank-local sparse PLM AMR step")
  call gather_mpi_amr_sparse_hierarchy_3d( &
    sparse_distribution, patch, 0, sparse_hierarchy, gathered_coarse, &
    gathered_coarse_temperature, gathered_fine, &
    gathered_fine_temperature, ok)
  call require(ok, "rank-local sparse PLM AMR gather")
  if (sparse_distribution%rank == 0) then
    call require(rank4_bits_match(gathered_coarse, serial_coarse) .and. &
      rank3_bits_match( &
        gathered_coarse_temperature, serial_coarse_temperature) .and. &
      rank4_bits_match(gathered_fine, serial_fine) .and. &
      rank3_bits_match(gathered_fine_temperature, serial_fine_temperature) &
      .and. same_real_bits(sparse_reflux, serial_reflux), &
      "rank-local sparse PLM AMR bit parity")
  end if
  saved_sparse_hierarchy = sparse_hierarchy
  if (sparse_distribution%nranks > 1) then
    sparse_reconstruction = "characteristic_plm"
    if (sparse_distribution%rank == sparse_distribution%nranks - 1) &
      sparse_reconstruction = "pcm"
    call advance_mpi_amr_sparse_periodic_hydro_3d( &
      sparse_distribution, species, patch, sparse_hierarchy%coarse_state, &
      sparse_hierarchy%coarse_temperature, dx, dy, dz, 1.0e-7_dp, &
      config%riemann_solver, sparse_reconstruction, "mc", sparse_flux_x, &
      sparse_flux_y, sparse_flux_z, ok)
    call require(.not. ok .and. sparse_hierarchies_bits_match( &
        sparse_hierarchy, saved_sparse_hierarchy), &
      "sparse periodic rank-dependent metadata rejection")
  end if
  call advance_mpi_amr_sparse_hydro_3d( &
    sparse_distribution, species, patch, sparse_hierarchy, &
    dx, dy, dz, 1.0e-7_dp, config%riemann_solver, &
    "invalid", "mc", sparse_reflux, ok)
  call require(.not. ok .and. same_real_bits(sparse_reflux, 0.0_dp) .and. &
    sparse_hierarchies_bits_match( &
      sparse_hierarchy, saved_sparse_hierarchy), &
    "sparse invalid-reconstruction transaction")
  if (sparse_distribution%nranks > 1) then
    sparse_reconstruction = "characteristic_plm"
    if (sparse_distribution%rank == sparse_distribution%nranks - 1) &
      sparse_reconstruction = "pcm"
    call advance_mpi_amr_sparse_hydro_3d( &
      sparse_distribution, species, patch, sparse_hierarchy, &
      dx, dy, dz, 1.0e-7_dp, config%riemann_solver, &
      sparse_reconstruction, "mc", sparse_reflux, ok)
    call require(.not. ok .and. same_real_bits(sparse_reflux, 0.0_dp) .and. &
      sparse_hierarchies_bits_match( &
        sparse_hierarchy, saved_sparse_hierarchy), &
      "sparse rank-dependent reconstruction rejection")
    mismatched_species = species
    if (sparse_distribution%rank == sparse_distribution%nranks - 1) then
      mismatched_species(1)%low_coefficients(1) = &
        nearest(mismatched_species(1)%low_coefficients(1), 1.0_dp)
    end if
    call advance_mpi_amr_sparse_hydro_3d( &
      sparse_distribution, mismatched_species, patch, sparse_hierarchy, &
      dx, dy, dz, 1.0e-7_dp, config%riemann_solver, &
      "characteristic_plm", "mc", sparse_reflux, ok)
    call require(.not. ok .and. same_real_bits(sparse_reflux, 0.0_dp) .and. &
      sparse_hierarchies_bits_match( &
        sparse_hierarchy, saved_sparse_hierarchy), &
      "sparse thermodynamics metadata rejection")
  end if
  call scatter_mpi_amr_sparse_hierarchy_3d( &
    sparse_distribution, species, patch, 0, coarse_state, &
    coarse_temperature, fine_state, fine_temperature, sparse_hierarchy, ok)
  call require(ok, "restore sparse hierarchy after AMR PLM test")
  if (sparse_distribution%nranks > 1) then
    mismatched_patch = patch
    if (sparse_distribution%rank == sparse_distribution%nranks - 1) then
      mismatched_patch%coarse_j_lower = patch%coarse_j_lower + 1
      mismatched_patch%coarse_j_upper = patch%coarse_j_upper + 1
    end if
    call gather_mpi_amr_sparse_hierarchy_3d( &
      sparse_distribution, mismatched_patch, 0, sparse_hierarchy, &
      gathered_coarse, gathered_coarse_temperature, gathered_fine, &
      gathered_fine_temperature, ok)
    call require(.not. ok, "rank-dependent sparse gather layout rejection")
  end if
  call gather_mpi_amr_sparse_hierarchy_3d( &
    sparse_distribution, patch, 0, sparse_hierarchy, gathered_coarse, &
    gathered_coarse_temperature, gathered_fine, &
    gathered_fine_temperature, ok)
  call require(ok, "rank-neutral sparse hierarchy gather")
  if (sparse_distribution%rank == 0) then
    call require(rank4_bits_match(gathered_coarse, coarse_state) .and. &
      rank3_bits_match(gathered_coarse_temperature, coarse_temperature) .and. &
      rank4_bits_match(gathered_fine, fine_state) .and. &
      rank3_bits_match(gathered_fine_temperature, fine_temperature), &
      "sparse scatter-gather bit round trip")
  else
    call require(size(gathered_coarse) == 0 .and. &
      size(gathered_coarse_temperature) == 0 .and. &
      size(gathered_fine) == 0 .and. &
      size(gathered_fine_temperature) == 0, &
      "sparse gather remains root-local")
  end if

  saved_coarse = coarse_state
  saved_coarse_temperature = coarse_temperature
  saved_fine = fine_state
  saved_fine_temperature = fine_temperature
  saved_integrals = initial_integrals
  saved_time = time
  saved_steps = steps
  saved_maximum_reflux = maximum_reflux
  if (distribution%nranks > 1) then
    inconsistent_root = merge(0, 1, distribution%rank == 0)
    call broadcast_mpi_amr_reactive_hierarchy_3d( &
      distribution, patch, inconsistent_root, coarse_state, &
      coarse_temperature, fine_state, fine_temperature, time, steps, &
      initial_integrals, maximum_reflux, ok)
    call require(.not. ok .and. &
      rank4_bits_match(coarse_state, saved_coarse) .and. &
      rank3_bits_match(coarse_temperature, saved_coarse_temperature) .and. &
      rank4_bits_match(fine_state, saved_fine) .and. &
      rank3_bits_match(fine_temperature, saved_fine_temperature) .and. &
      rank1_bits_match(initial_integrals, saved_integrals) .and. &
      same_real_bits(time, saved_time) .and. steps == saved_steps .and. &
      same_real_bits(maximum_reflux, saved_maximum_reflux), &
      "inconsistent broadcast root transaction")

    mismatched_fine_nx = patch%fine_nx()
    if (distribution%rank == distribution%nranks - 1) then
      mismatched_fine_nx = mismatched_fine_nx - 1
    end if
    allocate(mismatched_fine( &
      nvar, mismatched_fine_nx, patch%fine_ny(), patch%fine_nz()))
    allocate(mismatched_fine_temperature( &
      mismatched_fine_nx, patch%fine_ny(), patch%fine_nz()))
    mismatched_fine = fine_state(:, :mismatched_fine_nx, :, :)
    mismatched_fine_temperature = &
      fine_temperature(:mismatched_fine_nx, :, :)
    saved_mismatched_fine = mismatched_fine
    saved_mismatched_fine_temperature = mismatched_fine_temperature
    call broadcast_mpi_amr_reactive_hierarchy_3d( &
      distribution, patch, 0, coarse_state, coarse_temperature, &
      mismatched_fine, mismatched_fine_temperature, time, steps, &
      initial_integrals, maximum_reflux, ok)
    call require(.not. ok .and. &
      rank4_bits_match(coarse_state, saved_coarse) .and. &
      rank3_bits_match(coarse_temperature, saved_coarse_temperature) .and. &
      rank4_bits_match(mismatched_fine, saved_mismatched_fine) .and. &
      rank3_bits_match( &
        mismatched_fine_temperature, saved_mismatched_fine_temperature) .and. &
      rank1_bits_match(initial_integrals, saved_integrals) .and. &
      same_real_bits(time, saved_time) .and. steps == saved_steps .and. &
      same_real_bits(maximum_reflux, saved_maximum_reflux), &
      "mismatched broadcast shape transaction")
    deallocate( &
      mismatched_fine, mismatched_fine_temperature, &
      saved_mismatched_fine, saved_mismatched_fine_temperature)

    mismatched_species = species
    if (distribution%rank == distribution%nranks - 1) then
      mismatched_species(1)%low_coefficients(1) = &
        nearest(mismatched_species(1)%low_coefficients(1), 1.0_dp)
    end if
    call compute_mpi_amr_reactive_cfl_timestep_3d( &
      distribution, mismatched_species, patch, coarse_state, &
      coarse_temperature, fine_state, fine_temperature, dx, dy, dz, &
      0.30_dp, distributed_dt, ok)
    call require(.not. ok .and. same_real_bits(distributed_dt, 0.0_dp), &
      "mismatched thermodynamics rejection")
  end if

  serial_coarse = coarse_state
  serial_coarse_temperature = coarse_temperature
  serial_fine = fine_state
  serial_fine_temperature = fine_temperature
  call compute_amr_reactive_cfl_timestep_3d( &
    species, patch, serial_coarse, serial_coarse_temperature, &
    serial_fine, serial_fine_temperature, dx, dy, dz, 0.30_dp, &
    serial_dt, ok)
  call require(ok, "serial AMR timestep")
  call compute_mpi_amr_reactive_cfl_timestep_3d( &
    distribution, species, patch, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, dx, dy, dz, 0.30_dp, &
    distributed_dt, ok)
  call require(ok .and. same_real_bits(serial_dt, distributed_dt), &
    "distributed AMR timestep bit parity")
  distributed_dt = min(distributed_dt, 1.0e-7_dp)
  serial_dt = distributed_dt
  call advance_amr_reactive_hydro_3d( &
    species, patch, serial_coarse, serial_coarse_temperature, &
    serial_fine, serial_fine_temperature, dx, dy, dz, serial_dt, &
    config%riemann_solver, serial_reflux, ok)
  call require(ok, "serial AMR reference step")
  call advance_mpi_amr_reactive_hydro_3d( &
    distribution, species, patch, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, dx, dy, dz, distributed_dt, &
    config%riemann_solver, distributed_reflux, ok)
  call require(ok, "distributed AMR step")
  call require(rank4_bits_match(coarse_state, serial_coarse), &
    "distributed coarse state bit parity")
  call require(rank3_bits_match( &
    coarse_temperature, serial_coarse_temperature), &
    "distributed coarse temperature bit parity")
  call require(rank4_bits_match(fine_state, serial_fine), &
    "distributed fine state bit parity")
  call require(rank3_bits_match(fine_temperature, serial_fine_temperature), &
    "distributed fine temperature bit parity")
  call require(same_real_bits(serial_reflux, distributed_reflux), &
    "distributed reflux bit parity")

  serial_coarse = coarse_state
  serial_coarse_temperature = coarse_temperature
  serial_fine = fine_state
  serial_fine_temperature = fine_temperature
  call advance_amr_reactive_hydro_3d( &
    species, patch, serial_coarse, serial_coarse_temperature, &
    serial_fine, serial_fine_temperature, dx, dy, dz, serial_dt, &
    config%riemann_solver, serial_reflux, ok, "characteristic_plm", "mc")
  call require(ok, "serial characteristic-PLM AMR reference step")
  call advance_mpi_amr_reactive_hydro_3d( &
    distribution, species, patch, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, dx, dy, dz, distributed_dt, &
    config%riemann_solver, distributed_reflux, ok, &
    "characteristic_plm", "mc")
  call require(ok, "distributed characteristic-PLM AMR step")
  call require(rank4_bits_match(coarse_state, serial_coarse), &
    "distributed PLM coarse state bit parity")
  call require(rank3_bits_match( &
    coarse_temperature, serial_coarse_temperature), &
    "distributed PLM coarse temperature bit parity")
  call require(rank4_bits_match(fine_state, serial_fine), &
    "distributed PLM fine state bit parity")
  call require(rank3_bits_match(fine_temperature, serial_fine_temperature), &
    "distributed PLM fine temperature bit parity")
  call require(same_real_bits(serial_reflux, distributed_reflux), &
    "distributed PLM reflux bit parity")

  saved_coarse = coarse_state
  saved_coarse_temperature = coarse_temperature
  saved_fine = fine_state
  saved_fine_temperature = fine_temperature
  call advance_mpi_amr_reactive_hydro_3d( &
    distribution, species, patch, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, dx, dy, dz, distributed_dt, &
    "invalid", distributed_reflux, ok)
  call require(.not. ok .and. &
    rank4_bits_match(coarse_state, saved_coarse) .and. &
    rank3_bits_match(coarse_temperature, saved_coarse_temperature) .and. &
    rank4_bits_match(fine_state, saved_fine) .and. &
    rank3_bits_match(fine_temperature, saved_fine_temperature) .and. &
    same_real_bits(distributed_reflux, 0.0_dp), &
    "distributed invalid-solver transaction")
  call advance_mpi_amr_reactive_hydro_3d( &
    distribution, species, patch, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, dx, dy, dz, distributed_dt, &
    config%riemann_solver, distributed_reflux, ok, "invalid", "mc")
  call require(.not. ok .and. &
    rank4_bits_match(coarse_state, saved_coarse) .and. &
    rank3_bits_match(coarse_temperature, saved_coarse_temperature) .and. &
    rank4_bits_match(fine_state, saved_fine) .and. &
    rank3_bits_match(fine_temperature, saved_fine_temperature) .and. &
    same_real_bits(distributed_reflux, 0.0_dp), &
    "distributed invalid-reconstruction transaction")
  if (distribution%rank == 0) then
    write(*, '(a,i0)') "MPI AMR ranks: ", distribution%nranks
    write(*, '(a,i0,a,i0)') "Coarse slab planes: min=", &
      minval(distribution%coarse_counts), ", max=", &
      maxval(distribution%coarse_counts)
    write(*, '(a,i0,a,i0)') "Fine slab planes: min=", &
      minval(distribution%fine_counts), ", max=", &
      maxval(distribution%fine_counts)
    write(*, '(a)') "test_mpi_amr_reactive_3d: PASS"
  end if
  call MPI_Finalize(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Finalize failed"

contains

  subroutine run_adversarial_sparse_transport_tests()
    character(len=8), parameter :: orientation_name(6) = [ &
      character(len=8) :: "x-lower", "x-upper", "y-lower", "y-upper", &
      "z-lower", "z-upper"]
    real(dp), parameter :: adversarial_interval = 3.75e-5_dp
    real(dp), parameter :: adversarial_low_fraction = 1.0e-12_dp
    real(dp), parameter :: adversarial_high_fraction = 0.70_dp
    real(dp), allocatable :: case_coarse(:, :, :, :)
    real(dp), allocatable :: case_coarse_temperature(:, :, :)
    real(dp), allocatable :: case_fine(:, :, :, :)
    real(dp), allocatable :: case_fine_temperature(:, :, :)
    real(dp), allocatable :: saved_case_coarse(:, :, :, :)
    real(dp), allocatable :: saved_case_coarse_temperature(:, :, :)
    real(dp), allocatable :: saved_case_fine(:, :, :, :)
    real(dp), allocatable :: saved_case_fine_temperature(:, :, :)
    integer :: orientation, first_species, last_species
    logical :: local_ok

    first_species = reactive_species_component(1)
    last_species = reactive_species_component(size(species))
    allocate(case_coarse, mold=coarse_state)
    allocate(case_coarse_temperature, mold=coarse_temperature)
    allocate(case_fine, mold=fine_state)
    allocate(case_fine_temperature, mold=fine_temperature)
    allocate(saved_case_coarse, mold=coarse_state)
    allocate(saved_case_coarse_temperature, mold=coarse_temperature)
    allocate(saved_case_fine, mold=fine_state)
    allocate(saved_case_fine_temperature, mold=fine_temperature)

    do orientation = 1, 6
      call initialize_uniform_adversarial_hierarchy( &
        case_coarse, case_coarse_temperature, case_fine, &
        case_fine_temperature, local_ok)
      call require(local_ok, "adversarial uniform hierarchy " // &
        trim(orientation_name(orientation)))
      call configure_adversarial_interface_case( &
        species, orientation, case_coarse, case_coarse_temperature, &
        case_fine, case_fine_temperature, adversarial_low_fraction, &
        adversarial_high_fraction, local_ok)
      call require(local_ok .and. all(ieee_is_finite(case_coarse)) .and. &
        all(ieee_is_finite(case_fine)) .and. &
        species_nonnegative_4d(case_coarse, first_species, last_species) .and. &
        species_nonnegative_4d(case_fine, first_species, last_species), &
        "adversarial initial state " // trim(orientation_name(orientation)))
      saved_case_coarse = case_coarse
      saved_case_coarse_temperature = case_coarse_temperature
      saved_case_fine = case_fine
      saved_case_fine_temperature = case_fine_temperature

      call advance_amr_reactive_transport_3d( &
        species, transport, patch, case_coarse, case_coarse_temperature, &
        case_fine, case_fine_temperature, dx, dy, dz, &
        adversarial_interval, .false., .false., .true., .false., "mc", &
        serial_transport_theta, serial_transport_reflux, local_ok)
      call require(local_ok .and. serial_transport_theta >= 0.0_dp .and. &
        serial_transport_theta < 1.0_dp .and. &
        species_nonnegative_4d(case_coarse, first_species, last_species) .and. &
        species_nonnegative_4d(case_fine, first_species, last_species) .and. &
        .not. (rank4_bits_match(case_coarse, saved_case_coarse) .and. &
          rank4_bits_match(case_fine, saved_case_fine)), &
        "serial adversarial transport " // trim(orientation_name(orientation)))

      call scatter_mpi_amr_sparse_hierarchy_3d( &
        sparse_distribution, species, patch, 0, saved_case_coarse, &
        saved_case_coarse_temperature, saved_case_fine, &
        saved_case_fine_temperature, sparse_hierarchy, ok)
      call require(ok, "sparse adversarial scatter " // &
        trim(orientation_name(orientation)))
      saved_sparse_hierarchy = sparse_hierarchy
      call advance_mpi_amr_sparse_transport_3d( &
        sparse_distribution, species, transport, patch, sparse_hierarchy, &
        dx, dy, dz, adversarial_interval, .false., .false., .true., &
        .false., "mc", sparse_transport_theta, sparse_transport_reflux, ok)
      call require(ok .and. sparse_transport_theta >= 0.0_dp .and. &
        sparse_transport_theta < 1.0_dp .and. same_real_bits( &
        serial_transport_theta, sparse_transport_theta) .and. &
        same_real_bits(serial_transport_reflux, sparse_transport_reflux) .and. &
        species_nonnegative_4d(sparse_hierarchy%coarse_state, &
          first_species, last_species) .and. &
        species_nonnegative_4d(sparse_hierarchy%fine_state, &
          first_species, last_species) .and. &
        sparse_hierarchy%is_valid(sparse_distribution, patch), &
        "sparse adversarial transport " // trim(orientation_name(orientation)))
      call gather_mpi_amr_sparse_hierarchy_3d( &
        sparse_distribution, patch, 0, sparse_hierarchy, gathered_coarse, &
        gathered_coarse_temperature, gathered_fine, gathered_fine_temperature, &
        ok)
      call require(ok, "sparse adversarial gather " // &
        trim(orientation_name(orientation)))
      if (sparse_distribution%rank == 0) then
        call require(rank4_bits_match(gathered_coarse, case_coarse) .and. &
          rank3_bits_match(gathered_coarse_temperature, &
            case_coarse_temperature) .and. &
          rank4_bits_match(gathered_fine, case_fine) .and. &
          rank3_bits_match(gathered_fine_temperature, case_fine_temperature) &
          .and. species_nonnegative_4d(gathered_coarse, &
            first_species, last_species) .and. &
          species_nonnegative_4d(gathered_fine, first_species, last_species) &
          .and. .not. (rank4_bits_match(gathered_coarse, &
            saved_case_coarse) .and. rank4_bits_match(gathered_fine, &
            saved_case_fine)), &
          "sparse adversarial serial bit parity " // &
          trim(orientation_name(orientation)))
        write(*, '(a,1x,a,1x,a,es24.16,1x,a,es24.16)') &
          "adversarial sparse", trim(orientation_name(orientation)), &
          "theta=", sparse_transport_theta, "reflux=", sparse_transport_reflux
      end if
    end do
  end subroutine run_adversarial_sparse_transport_tests

  subroutine run_sparse_full_transport_test()
    real(dp), parameter :: full_interval = 1.0e-8_dp
    real(dp), parameter :: full_rtol = 2.0e-7_dp
    real(dp), parameter :: full_atol = 1.0e-12_dp
    logical :: local_ok, local_condition

    call scatter_mpi_amr_sparse_hierarchy_3d( &
      sparse_distribution, species, patch, 0, coarse_state, &
      coarse_temperature, fine_state, fine_temperature, sparse_hierarchy, ok)
    call require(ok, "sparse full-operator initial scatter")
    if (sparse_distribution%rank == 0) then
      serial_coarse = coarse_state
      serial_coarse_temperature = coarse_temperature
      serial_fine = fine_state
      serial_fine_temperature = fine_temperature
      call advance_amr_reactive_full_3d( &
        species, reactions, transport, patch, serial_coarse, &
        serial_coarse_temperature, serial_fine, serial_fine_temperature, &
        dx, dy, dz, full_interval, "rusanov", .true., full_rtol, full_atol, &
        .true., .true., .true., .true., .true., serial_full_theta, &
        serial_full_reflux, local_ok, "characteristic_plm", "mc")
    else
      local_ok = .true.
      serial_full_theta = 1.0_dp
      serial_full_reflux = 0.0_dp
    end if
    call require(local_ok, "serial active R-T-H-T-R reference")
    call advance_mpi_amr_sparse_full_3d( &
      sparse_distribution, species, reactions, transport, patch, &
      sparse_hierarchy, dx, dy, dz, full_interval, "rusanov", &
      "characteristic_plm", "mc", .true., full_rtol, full_atol, .true., &
      .true., .true., .true., .true., sparse_full_theta, &
      sparse_full_reflux, ok)
    call require(ok, "sparse active R-T-H-T-R step")
    call gather_mpi_amr_sparse_hierarchy_3d( &
      sparse_distribution, patch, 0, sparse_hierarchy, gathered_coarse, &
      gathered_coarse_temperature, gathered_fine, gathered_fine_temperature, &
      ok)
    call require(ok, "sparse active R-T-H-T-R gather")
    local_condition = .true.
    if (sparse_distribution%rank == 0) then
      local_condition = rank4_bits_match(gathered_coarse, serial_coarse) .and. &
        rank3_bits_match(gathered_coarse_temperature, &
          serial_coarse_temperature) .and. &
        rank4_bits_match(gathered_fine, serial_fine) .and. &
        rank3_bits_match(gathered_fine_temperature, serial_fine_temperature) &
        .and. same_real_bits(sparse_full_theta, serial_full_theta) .and. &
        same_real_bits(sparse_full_reflux, serial_full_reflux)
    end if
    call require(local_condition, "serial/sparse active R-T-H-T-R parity")

    saved_sparse_hierarchy = sparse_hierarchy
    call advance_mpi_amr_sparse_full_3d( &
      sparse_distribution, species, reactions, transport, patch, &
      sparse_hierarchy, dx, dy, dz, full_interval, "invalid", &
      "characteristic_plm", "mc", .true., full_rtol, full_atol, .true., &
      .true., .true., .true., .true., sparse_full_theta, &
      sparse_full_reflux, ok)
    call require(.not. ok .and. same_real_bits(sparse_full_theta, 1.0_dp) &
      .and. same_real_bits(sparse_full_reflux, 0.0_dp) .and. &
      sparse_hierarchies_bits_match(sparse_hierarchy, saved_sparse_hierarchy), &
      "sparse full invalid-middle-hydro rollback")

    if (sparse_distribution%nranks > 1) then
      mismatched_transport = transport
      if (sparse_distribution%rank == sparse_distribution%nranks - 1) then
        mismatched_transport(1)%well_depth = nearest( &
          mismatched_transport(1)%well_depth, 1.0_dp)
      end if
      saved_sparse_hierarchy = sparse_hierarchy
      call advance_mpi_amr_sparse_full_3d( &
        sparse_distribution, species, reactions, mismatched_transport, patch, &
        sparse_hierarchy, dx, dy, dz, full_interval, "rusanov", &
        "characteristic_plm", "mc", .true., full_rtol, full_atol, .true., &
        .true., .true., .true., .true., sparse_full_theta, &
        sparse_full_reflux, ok)
      call require(.not. ok .and. same_real_bits(sparse_full_theta, 1.0_dp) &
        .and. same_real_bits(sparse_full_reflux, 0.0_dp) .and. &
        sparse_hierarchies_bits_match( &
          sparse_hierarchy, saved_sparse_hierarchy), &
        "sparse full transport metadata rollback")
    end if
  end subroutine run_sparse_full_transport_test

  subroutine run_sparse_nonfinite_guard_tests()
    type(mpi_amr_sparse_hierarchy_3d) :: baseline, candidate
    real(dp), allocatable :: rejected_halo_state(:, :, :, :)
    real(dp), allocatable :: rejected_halo_temperature(:, :, :)
    real(dp), allocatable :: rejected_flux_x(:, :, :, :)
    real(dp), allocatable :: rejected_flux_y(:, :, :, :)
    real(dp), allocatable :: rejected_flux_z(:, :, :, :)
    real(dp) :: nan_value

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)

    call compute_mpi_amr_sparse_transport_timestep_3d( &
      sparse_distribution, species, transport, patch, sparse_hierarchy, &
      nan_value, dy, dz, 0.35_dp, .true., .true., .true., sparse_dt, &
      sparse_maximum_diffusivity, ok)
    call require(.not. ok .and. same_real_bits(sparse_dt, 0.0_dp) .and. &
      same_real_bits(sparse_maximum_diffusivity, 0.0_dp), &
      "sparse transport timestep nonfinite-input rejection")

    call compute_mpi_amr_sparse_cfl_timestep_3d( &
      sparse_distribution, species, patch, sparse_hierarchy, &
      dx, dy, nan_value, 0.30_dp, sparse_dt, ok)
    call require(.not. ok .and. same_real_bits(sparse_dt, 0.0_dp), &
      "sparse CFL timestep nonfinite-input rejection")

    candidate = sparse_hierarchy
    candidate%coarse_temperature(1, 1, 1) = nan_value
    call build_mpi_amr_sparse_periodic_halos_3d( &
      sparse_distribution, candidate%coarse_state, &
      candidate%coarse_temperature, rejected_halo_state, &
      rejected_halo_temperature, ok)
    call require(.not. ok, "sparse periodic halo nonfinite rejection")

    baseline = sparse_hierarchy
    candidate = sparse_hierarchy
    call advance_mpi_amr_sparse_periodic_hydro_3d( &
      sparse_distribution, species, patch, candidate%coarse_state, &
      candidate%coarse_temperature, nan_value, dy, dz, 1.0e-7_dp, &
      config%riemann_solver, "characteristic_plm", "mc", rejected_flux_x, &
      rejected_flux_y, rejected_flux_z, ok)
    call require(.not. ok .and. sparse_hierarchies_bits_match( &
      candidate, baseline), "sparse periodic hydro nonfinite rollback")

    candidate = sparse_hierarchy
    baseline = sparse_hierarchy
    call advance_mpi_amr_sparse_transport_3d( &
      sparse_distribution, species, transport, patch, candidate, nan_value, &
      dy, dz, 1.0e-8_dp, .true., .true., .true., .true., "mc", &
      sparse_transport_theta, sparse_transport_reflux, ok)
    call require(.not. ok .and. same_real_bits( &
      sparse_transport_theta, 1.0_dp) .and. same_real_bits( &
      sparse_transport_reflux, 0.0_dp) .and. &
      sparse_hierarchies_bits_match(candidate, baseline), &
      "sparse transport nonfinite rollback")
  end subroutine run_sparse_nonfinite_guard_tests

  subroutine initialize_uniform_adversarial_hierarchy( &
      case_coarse, case_coarse_temperature, case_fine, case_fine_temperature, &
      ok_out)
    real(dp), allocatable, intent(out) :: case_coarse(:, :, :, :)
    real(dp), allocatable, intent(out) :: case_coarse_temperature(:, :, :)
    real(dp), allocatable, intent(out) :: case_fine(:, :, :, :)
    real(dp), allocatable, intent(out) :: case_fine_temperature(:, :, :)
    logical, intent(out) :: ok_out

    type(reactive_3d_config) :: uniform_config, fine_config
    real(dp), allocatable :: mass_fractions(:)
    real(dp), allocatable :: recovered_temperature(:, :, :)
    real(dp) :: density, fine_density
    logical :: local_ok

    ok_out = .false.
    uniform_config = reactive_3d_config()
    uniform_config%nx = n
    uniform_config%ny = n
    uniform_config%nz = n
    uniform_config%x_upper = real(n, dp) * dx
    uniform_config%y_upper = real(n, dp) * dy
    uniform_config%z_upper = real(n, dp) * dz
    uniform_config%problem = "uniform_reactor"
    uniform_config%thermo_model = "elementary"
    uniform_config%initial_temperature = 1200.0_dp
    uniform_config%initial_velocity_x = 15.0_dp
    uniform_config%initial_velocity_y = -7.0_dp
    uniform_config%initial_velocity_z = 3.0_dp

    fine_config = uniform_config
    fine_config%nx = patch%fine_nx()
    fine_config%ny = patch%fine_ny()
    fine_config%nz = patch%fine_nz()
    fine_config%x_lower = real(patch%coarse_i_lower - 1, dp) * dx
    fine_config%x_upper = real(patch%coarse_i_upper, dp) * dx
    fine_config%y_lower = real(patch%coarse_j_lower - 1, dp) * dy
    fine_config%y_upper = real(patch%coarse_j_upper, dp) * dy
    fine_config%z_lower = real(patch%coarse_k_lower - 1, dp) * dz
    fine_config%z_upper = real(patch%coarse_k_upper, dp) * dz

    allocate(case_coarse(reactive_nvar(size(species)), n, n, n))
    allocate(case_coarse_temperature(n, n, n))
    allocate(case_fine(reactive_nvar(size(species)), patch%fine_nx(), &
      patch%fine_ny(), patch%fine_nz()))
    allocate(case_fine_temperature( &
      patch%fine_nx(), patch%fine_ny(), patch%fine_nz()))
    allocate(mass_fractions(size(species)))
    allocate(recovered_temperature(n, n, n))
    call initialize_reactive_problem_3d( &
      species, uniform_config, x, y, z, case_coarse, &
      case_coarse_temperature, density, mass_fractions, local_ok)
    if (.not. local_ok) return
    call initialize_reactive_problem_3d( &
      species, fine_config, xf, yf, zf, case_fine, case_fine_temperature, &
      fine_density, mass_fractions, local_ok)
    if (.not. local_ok) return
    call average_down_3d(case_coarse, case_fine, patch, local_ok)
    if (.not. local_ok) return
    call recover_reactive_temperatures_3d( &
      species, case_coarse, case_coarse_temperature, n, n, n, &
      recovered_temperature, local_ok)
    if (.not. local_ok) return
    case_coarse_temperature = recovered_temperature
    ok_out = density > 0.0_dp .and. fine_density > 0.0_dp
  end subroutine initialize_uniform_adversarial_hierarchy

  subroutine configure_adversarial_interface_case( &
      active_species, orientation, case_coarse, case_coarse_temperature, &
      case_fine, case_fine_temperature, low_fraction, high_fraction, ok_out)
    type(nasa7_species), intent(in) :: active_species(:)
    integer, intent(in) :: orientation
    real(dp), intent(inout) :: case_coarse(:, :, :, :)
    real(dp), intent(inout) :: case_coarse_temperature(:, :, :)
    real(dp), intent(inout) :: case_fine(:, :, :, :)
    real(dp), intent(inout) :: case_fine_temperature(:, :, :)
    real(dp), intent(in) :: low_fraction, high_fraction
    logical, intent(out) :: ok_out

    real(dp), allocatable :: recovered_temperature(:, :, :)
    real(dp) :: fraction
    logical :: local_ok
    integer :: i, j, k, coarse_i, coarse_j, coarse_k, ratio

    ok_out = .false.
    if (orientation < 1 .or. orientation > 6 .or. &
        low_fraction <= 0.0_dp .or. high_fraction <= low_fraction .or. &
        high_fraction >= 1.0_dp) return
    if (any(shape(case_coarse) /= [reactive_nvar(size(active_species)), &
        n, n, n]) .or. any(shape(case_coarse_temperature) /= [n, n, n]) .or. &
        any(shape(case_fine) /= [reactive_nvar(size(active_species)), &
          patch%fine_nx(), patch%fine_ny(), patch%fine_nz()]) .or. &
        any(shape(case_fine_temperature) /= &
          [patch%fine_nx(), patch%fine_ny(), patch%fine_nz()])) return

    do k = 1, n
      do j = 1, n
        do i = 1, n
          fraction = high_fraction
          if (adversarial_low_cell(orientation, i, j, k)) &
            fraction = low_fraction
          call set_adversarial_species_fraction( &
            case_coarse(:, i, j, k), size(active_species), fraction, local_ok)
          if (.not. local_ok) return
        end do
      end do
    end do

    ratio = patch%refinement_ratio
    do k = 1, patch%fine_nz()
      coarse_k = patch%coarse_k_lower + (k - 1) / ratio
      do j = 1, patch%fine_ny()
        coarse_j = patch%coarse_j_lower + (j - 1) / ratio
        do i = 1, patch%fine_nx()
          coarse_i = patch%coarse_i_lower + (i - 1) / ratio
          case_fine(:, i, j, k) = case_coarse(:, coarse_i, coarse_j, coarse_k)
        end do
      end do
    end do

    allocate(recovered_temperature, mold=case_coarse_temperature)
    call recover_reactive_temperatures_3d( &
      active_species, case_coarse, case_coarse_temperature, n, n, n, &
      recovered_temperature, local_ok)
    if (.not. local_ok) return
    case_coarse_temperature = recovered_temperature
    deallocate(recovered_temperature)
    allocate(recovered_temperature, mold=case_fine_temperature)
    call recover_reactive_temperatures_3d( &
      active_species, case_fine, case_fine_temperature, patch%fine_nx(), &
      patch%fine_ny(), patch%fine_nz(), recovered_temperature, local_ok)
    if (.not. local_ok) return
    case_fine_temperature = recovered_temperature
    ok_out = all(ieee_is_finite(case_coarse_temperature))
    if (ok_out) ok_out = all(ieee_is_finite(case_fine_temperature))
    if (ok_out) ok_out = minval(case_coarse_temperature) > 0.0_dp
    if (ok_out) ok_out = minval(case_fine_temperature) > 0.0_dp
  end subroutine configure_adversarial_interface_case

  logical function adversarial_low_cell(orientation, i, j, k) result(low)
    integer, intent(in) :: orientation, i, j, k

    low = .false.
    select case (orientation)
    case (1)
      low = i == patch%coarse_i_lower - 1
    case (2)
      low = i == patch%coarse_i_upper + 1
    case (3)
      low = j == patch%coarse_j_lower - 1
    case (4)
      low = j == patch%coarse_j_upper + 1
    case (5)
      low = k == patch%coarse_k_lower - 1
    case (6)
      low = k == patch%coarse_k_upper + 1
    end select
  end function adversarial_low_cell

  subroutine set_adversarial_species_fraction( &
      state, nspecies, fraction, ok_out)
    real(dp), intent(inout) :: state(:)
    integer, intent(in) :: nspecies
    real(dp), intent(in) :: fraction
    logical, intent(out) :: ok_out

    integer :: target_component, reservoir_component
    real(dp) :: rho, other_species_mass, reservoir_mass

    ok_out = .false.
    target_component = reactive_species_component(1)
    reservoir_component = reactive_species_component(nspecies)
    if (size(state) /= reactive_nvar(nspecies)) return
    if (.not. ieee_is_finite(fraction)) return
    if (fraction <= 0.0_dp .or. fraction >= 1.0_dp) return
    if (.not. all(ieee_is_finite(state))) return
    rho = state(irho)
    other_species_mass = rho - state(target_component) - &
      state(reservoir_component)
    reservoir_mass = rho * (1.0_dp - fraction) - other_species_mass
    if (.not. ieee_is_finite(rho)) return
    if (rho <= 0.0_dp) return
    if (.not. ieee_is_finite(reservoir_mass)) return
    if (reservoir_mass < 0.0_dp) return
    state(target_component) = rho * fraction
    state(reservoir_component) = reservoir_mass
    ok_out = all(ieee_is_finite(state))
    if (ok_out) ok_out = minval( &
      state(target_component:reservoir_component)) >= 0.0_dp
  end subroutine set_adversarial_species_fraction

  logical function species_nonnegative_4d(state, first_component, &
      last_component) result(nonnegative)
    real(dp), intent(in) :: state(:, :, :, :)
    integer, intent(in) :: first_component, last_component

    nonnegative = first_component >= 1 .and. &
      last_component <= size(state, 1) .and. first_component <= last_component
    if (.not. nonnegative) return
    if (size(state, 2) == 0) return
    nonnegative = all(ieee_is_finite( &
      state(first_component:last_component, :, :, :)))
    if (nonnegative) nonnegative = minval( &
      state(first_component:last_component, :, :, :)) >= -1.0e-14_dp
  end function species_nonnegative_4d

  subroutine fine_patch_centers(patch, config, dx, dy, dz, x, y, z)
    type(amr_patch_3d), intent(in) :: patch
    type(reactive_3d_config), intent(in) :: config
    real(dp), intent(in) :: dx, dy, dz
    real(dp), intent(out) :: x(:), y(:), z(:)
    integer :: i

    do i = 1, size(x)
      x(i) = config%x_lower + &
        real(patch%coarse_i_lower - 1, dp) * dx + &
        (real(i, dp) - 0.5_dp) * dx / real(patch%refinement_ratio, dp)
    end do
    do i = 1, size(y)
      y(i) = config%y_lower + &
        real(patch%coarse_j_lower - 1, dp) * dy + &
        (real(i, dp) - 0.5_dp) * dy / real(patch%refinement_ratio, dp)
    end do
    do i = 1, size(z)
      z(i) = config%z_lower + &
        real(patch%coarse_k_lower - 1, dp) * dz + &
        (real(i, dp) - 0.5_dp) * dz / real(patch%refinement_ratio, dp)
    end do
  end subroutine fine_patch_centers

  pure logical function same_real_bits(left, right) result(same)
    real(dp), intent(in) :: left, right

    same = transfer(left, 0_int64) == transfer(right, 0_int64)
  end function same_real_bits

  pure logical function rank4_bits_match(left, right) result(matches)
    real(dp), intent(in) :: left(:, :, :, :), right(:, :, :, :)
    integer :: i, j, k, component

    matches = all(shape(left) == shape(right))
    if (.not. matches) return
    do k = 1, size(left, 4)
      do j = 1, size(left, 3)
        do i = 1, size(left, 2)
          do component = 1, size(left, 1)
            if (.not. same_real_bits( &
                left(component, i, j, k), right(component, i, j, k))) then
              matches = .false.
              return
            end if
          end do
        end do
      end do
    end do
  end function rank4_bits_match

  pure logical function rank3_bits_match(left, right) result(matches)
    real(dp), intent(in) :: left(:, :, :), right(:, :, :)
    integer :: i, j, k

    matches = all(shape(left) == shape(right))
    if (.not. matches) return
    do k = 1, size(left, 3)
      do j = 1, size(left, 2)
        do i = 1, size(left, 1)
          if (.not. same_real_bits(left(i, j, k), right(i, j, k))) then
            matches = .false.
            return
          end if
        end do
      end do
    end do
  end function rank3_bits_match

  pure logical function rank1_bits_match(left, right) result(matches)
    real(dp), intent(in) :: left(:), right(:)
    integer :: i

    matches = size(left) == size(right)
    if (.not. matches) return
    do i = 1, size(left)
      if (.not. same_real_bits(left(i), right(i))) then
        matches = .false.
        return
      end if
    end do
  end function rank1_bits_match

  pure logical function sparse_periodic_halos_match( &
      sparse_distribution, global_state, global_temperature, halo_state, &
      halo_temperature) result(matches)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: sparse_distribution
    real(dp), intent(in) :: global_state(:, :, :, :)
    real(dp), intent(in) :: global_temperature(:, :, :)
    real(dp), intent(in) :: halo_state(:, :, :, :)
    real(dp), intent(in) :: halo_temperature(:, :, :)
    integer :: component, global_i, halo_i, j, k

    matches = size(halo_state, 2) == sparse_distribution%coarse_count + 4
    if (.not. matches) return
    do k = 1, size(global_state, 4)
      do j = 1, size(global_state, 3)
        do halo_i = 1, size(halo_state, 2)
          global_i = 1 + modulo( &
            sparse_distribution%coarse_first + halo_i - 4, &
            size(global_state, 2))
          if (.not. same_real_bits( &
              halo_temperature(halo_i, j, k), &
              global_temperature(global_i, j, k))) then
            matches = .false.
            return
          end if
          do component = 1, size(global_state, 1)
            if (.not. same_real_bits( &
                halo_state(component, halo_i, j, k), &
                global_state(component, global_i, j, k))) then
              matches = .false.
              return
            end if
          end do
        end do
      end do
    end do
  end function sparse_periodic_halos_match

  pure logical function sparse_hierarchies_bits_match( &
      left, right) result(matches)
    type(mpi_amr_sparse_hierarchy_3d), intent(in) :: left, right

    matches = left%nvar == right%nvar .and. &
      rank4_bits_match(left%coarse_state, right%coarse_state) .and. &
      rank3_bits_match( &
        left%coarse_temperature, right%coarse_temperature) .and. &
      rank4_bits_match(left%fine_state, right%fine_state) .and. &
      rank3_bits_match(left%fine_temperature, right%fine_temperature)
  end function sparse_hierarchies_bits_match

  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      write(*, '(a,i0,a,a)') "FAIL rank ", distribution%rank, ": ", &
        trim(label)
      error stop 1
    end if
  end subroutine require

end program test_mpi_amr_reactive_3d
