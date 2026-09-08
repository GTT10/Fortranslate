program test_amr_reactive_3d_checkpoint
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_elementary_transport
  use mesh_3d_mod, only: uniform_cell_centers_3d
  use reactive_1d_mod, only: reactive_nvar
  use simulation_config_reactive_3d_mod, only: reactive_3d_config
  use simulation_config_amr_reactive_3d_mod, only: amr_reactive_3d_config
  use reactive_entropy_wave_3d_problem_mod, only: &
    initialize_reactive_problem_3d
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use amr_hierarchy_3d_mod, only: &
    amr_patch_3d, initialize_amr_patch_3d, average_down_3d, &
    composite_integrals_amr_3d
  use amr_reactive_3d_checkpoint_mod, only: &
    write_amr_reactive_3d_checkpoint, read_amr_reactive_3d_checkpoint
  implicit none

  integer, parameter :: n = 4, ratio = 2
  character(len=*), parameter :: checkpoint_path = &
    "test_amr_reactive_3d_checkpoint.chk"
  character(len=*), parameter :: selected_checkpoint_path = &
    "test_selected_amr_reactive_3d_checkpoint.chk"
  character(len=*), parameter :: transport_checkpoint_path = &
    "test_transport_amr_reactive_3d_checkpoint.chk"
  character(len=*), parameter :: selected_transport_checkpoint_path = &
    "test_selected_transport_amr_reactive_3d_checkpoint.chk"
  character(len=*), parameter :: invalid_checkpoint_path = &
    "test_invalid_selected_amr_reactive_3d_checkpoint.chk"
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(elementary_reaction), allocatable :: incompatible_reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(gas_transport_species), allocatable :: incompatible_transport(:)
  type(reactive_3d_config) :: config, fine_config, incompatible_config
  type(reactive_3d_config) :: selected_config, selected_transport_config
  type(reactive_3d_config) :: transport_config
  type(amr_reactive_3d_config) :: amr_config
  type(amr_patch_3d) :: patch
  real(dp), allocatable :: coarse_state(:, :, :, :)
  real(dp), allocatable :: fine_state(:, :, :, :)
  real(dp), allocatable :: coarse_temperature(:, :, :)
  real(dp), allocatable :: fine_temperature(:, :, :)
  real(dp), allocatable :: recovered_temperature(:, :, :)
  real(dp), allocatable :: saved_coarse_state(:, :, :, :)
  real(dp), allocatable :: saved_fine_state(:, :, :, :)
  real(dp), allocatable :: saved_coarse_temperature(:, :, :)
  real(dp), allocatable :: saved_fine_temperature(:, :, :)
  real(dp), allocatable :: initial_integrals(:), saved_integrals(:)
  real(dp), allocatable :: mass_fractions(:)
  real(dp), allocatable :: selected_composition(:)
  real(dp), allocatable :: incompatible_composition(:)
  real(dp), allocatable :: x(:), y(:), z(:), xf(:), yf(:), zf(:)
  real(dp) :: dx, dy, dz, coarse_density, fine_density
  real(dp) :: time, maximum_reflux, saved_time, saved_maximum_reflux
  real(dp) :: transport_diffusivity, transport_theta
  real(dp) :: saved_transport_diffusivity, saved_transport_theta
  logical :: ok
  integer :: steps, saved_steps, nvar
  character(len=1024) :: message
  character(len=64) :: bundle_sha256, incompatible_bundle_sha256

  call delete_checkpoint(checkpoint_path)
  call delete_checkpoint(selected_checkpoint_path)
  call delete_checkpoint(transport_checkpoint_path)
  call delete_checkpoint(selected_transport_checkpoint_path)
  call delete_checkpoint(invalid_checkpoint_path)
  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "elementary thermodynamics load")
  call load_h2o2_elementary_mechanism(reactions, ok)
  call require(ok, "elementary mechanism load")
  call load_h2o2_elementary_transport(transport, ok)
  call require(ok, "elementary transport load")
  config = reactive_3d_config()
  config%nx = n
  config%ny = n
  config%nz = n
  config%thermo_model = "elementary"
  config%riemann_solver = "pelec"
  config%reconstruction = "characteristic_plm"
  config%limiter = "mc"
  config%maximum_steps = 100
  config%final_time = 1.0e-6_dp
  amr_config = amr_reactive_3d_config()
  amr_config%coarse_i_lower = 2
  amr_config%coarse_i_upper = 3
  amr_config%coarse_j_lower = 2
  amr_config%coarse_j_upper = 3
  amr_config%coarse_k_lower = 2
  amr_config%coarse_k_upper = 3
  call initialize_amr_patch_3d( &
    n, n, n, 2, 3, 2, 3, 2, 3, ratio, patch, ok)
  call require(ok .and. patch%is_strictly_interior(), &
    "checkpoint test patch")

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
    coarse_density, mass_fractions, ok)
  call require(ok .and. coarse_density > 0.0_dp, &
    "coarse checkpoint state initialization")
  fine_config = config
  fine_config%nx = patch%fine_nx()
  fine_config%ny = patch%fine_ny()
  fine_config%nz = patch%fine_nz()
  call initialize_reactive_problem_3d( &
    species, fine_config, xf, yf, zf, fine_state, fine_temperature, &
    fine_density, mass_fractions, ok)
  call require(ok .and. fine_density > 0.0_dp, &
    "fine checkpoint state initialization")
  call average_down_3d(coarse_state, fine_state, patch, ok)
  call require(ok, "checkpoint average-down")
  call recover_reactive_temperatures_3d( &
    species, coarse_state, coarse_temperature, n, n, n, &
    recovered_temperature, ok)
  call require(ok, "checkpoint coarse temperature recovery")
  coarse_temperature = recovered_temperature
  call composite_integrals_amr_3d( &
    coarse_state, fine_state, patch, dx, dy, dz, initial_integrals, ok)
  call require(ok, "checkpoint composite integral")

  time = 2.5e-7_dp
  steps = 2
  maximum_reflux = 3.25_dp
  call write_amr_reactive_3d_checkpoint( &
    checkpoint_path, species, config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message)
  call require(ok, "checkpoint write: " // trim(message))
  saved_coarse_state = coarse_state
  saved_fine_state = fine_state
  saved_coarse_temperature = coarse_temperature
  saved_fine_temperature = fine_temperature
  saved_integrals = initial_integrals
  saved_time = time
  saved_steps = steps
  saved_maximum_reflux = maximum_reflux

  coarse_state = -1.0_dp
  fine_state = -2.0_dp
  coarse_temperature = 300.0_dp
  fine_temperature = 400.0_dp
  initial_integrals = -3.0_dp
  time = 0.0_dp
  steps = 0
  maximum_reflux = 0.0_dp
  call read_amr_reactive_3d_checkpoint( &
    checkpoint_path, species, config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message)
  call require(ok, "checkpoint read: " // trim(message))
  call require(maxval(abs(coarse_state - saved_coarse_state)) == 0.0_dp, &
    "coarse state exact round trip")
  call require(maxval(abs(fine_state - saved_fine_state)) == 0.0_dp, &
    "fine state exact round trip")
  call require(maxval(abs(coarse_temperature - &
    saved_coarse_temperature)) <= 64.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, maxval(abs(saved_coarse_temperature))), &
    "coarse temperature round trip")
  call require(maxval(abs(fine_temperature - &
    saved_fine_temperature)) <= 64.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, maxval(abs(saved_fine_temperature))), &
    "fine temperature round trip")
  call require(all(initial_integrals == saved_integrals) .and. &
    time == saved_time .and. steps == saved_steps .and. &
    maximum_reflux == saved_maximum_reflux, &
    "checkpoint metadata exact round trip")

  transport_config = config
  transport_config%transport_enabled = .true.
  transport_diffusivity = 7.25_dp
  transport_theta = 0.375_dp
  saved_transport_diffusivity = transport_diffusivity
  saved_transport_theta = transport_theta
  call write_amr_reactive_3d_checkpoint( &
    transport_checkpoint_path, species, transport_config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message, &
    transport=transport, &
    maximum_transport_diffusivity=transport_diffusivity, &
    minimum_transport_theta=transport_theta)
  call require(ok, "fixed transport checkpoint write: " // trim(message))
  call require_checkpoint_schema( &
    transport_checkpoint_path, 4, "fixed transport schema four")

  coarse_state = -1.0_dp
  fine_state = -2.0_dp
  coarse_temperature = 300.0_dp
  fine_temperature = 400.0_dp
  initial_integrals = -3.0_dp
  time = 0.0_dp
  steps = 0
  maximum_reflux = 0.0_dp
  transport_diffusivity = 0.0_dp
  transport_theta = 1.0_dp
  call read_amr_reactive_3d_checkpoint( &
    transport_checkpoint_path, species, transport_config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message, &
    transport=transport, &
    maximum_transport_diffusivity=transport_diffusivity, &
    minimum_transport_theta=transport_theta)
  call require(ok, "fixed transport checkpoint read: " // trim(message))
  call require_checkpoint_unchanged( &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, &
    saved_coarse_state, saved_coarse_temperature, &
    saved_fine_state, saved_fine_temperature, &
    saved_time, saved_steps, saved_integrals, saved_maximum_reflux, &
    "fixed transport checkpoint exact round trip", &
    saved_maximum_diffusivity=saved_transport_diffusivity, &
    saved_minimum_theta=saved_transport_theta, &
    maximum_diffusivity=transport_diffusivity, minimum_theta=transport_theta)

  incompatible_transport = transport
  incompatible_transport(1)%name = "O2"
  call read_amr_reactive_3d_checkpoint( &
    transport_checkpoint_path, species, transport_config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message, &
    transport=incompatible_transport, &
    maximum_transport_diffusivity=transport_diffusivity, &
    minimum_transport_theta=transport_theta)
  call require(.not. ok, "transport database/name mismatch rejection")
  call require_checkpoint_unchanged( &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, &
    saved_coarse_state, saved_coarse_temperature, &
    saved_fine_state, saved_fine_temperature, &
    saved_time, saved_steps, saved_integrals, saved_maximum_reflux, &
    "transport name mismatch transaction", &
    saved_maximum_diffusivity=saved_transport_diffusivity, &
    saved_minimum_theta=saved_transport_theta, &
    maximum_diffusivity=transport_diffusivity, minimum_theta=transport_theta)

  incompatible_transport = transport
  incompatible_transport(1)%geometry = 2
  call read_amr_reactive_3d_checkpoint( &
    transport_checkpoint_path, species, transport_config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message, &
    transport=incompatible_transport, &
    maximum_transport_diffusivity=transport_diffusivity, &
    minimum_transport_theta=transport_theta)
  call require(.not. ok, "transport geometry mismatch rejection")
  call require_checkpoint_unchanged( &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, &
    saved_coarse_state, saved_coarse_temperature, &
    saved_fine_state, saved_fine_temperature, &
    saved_time, saved_steps, saved_integrals, saved_maximum_reflux, &
    "transport geometry mismatch transaction", &
    saved_maximum_diffusivity=saved_transport_diffusivity, &
    saved_minimum_theta=saved_transport_theta, &
    maximum_diffusivity=transport_diffusivity, minimum_theta=transport_theta)

  incompatible_transport = transport
  incompatible_transport(1)%well_depth = &
    incompatible_transport(1)%well_depth + 1.0_dp
  call read_amr_reactive_3d_checkpoint( &
    transport_checkpoint_path, species, transport_config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message, &
    transport=incompatible_transport, &
    maximum_transport_diffusivity=transport_diffusivity, &
    minimum_transport_theta=transport_theta)
  call require(.not. ok, "transport value mismatch rejection")
  call require_checkpoint_unchanged( &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, &
    saved_coarse_state, saved_coarse_temperature, &
    saved_fine_state, saved_fine_temperature, &
    saved_time, saved_steps, saved_integrals, saved_maximum_reflux, &
    "transport value mismatch transaction", &
    saved_maximum_diffusivity=saved_transport_diffusivity, &
    saved_minimum_theta=saved_transport_theta, &
    maximum_diffusivity=transport_diffusivity, minimum_theta=transport_theta)

  incompatible_config = transport_config
  incompatible_config%viscosity_enabled = .false.
  call read_amr_reactive_3d_checkpoint( &
    transport_checkpoint_path, species, incompatible_config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message, &
    transport=transport, &
    maximum_transport_diffusivity=transport_diffusivity, &
    minimum_transport_theta=transport_theta)
  call require(.not. ok, "transport policy mismatch rejection")
  call require_checkpoint_unchanged( &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, &
    saved_coarse_state, saved_coarse_temperature, &
    saved_fine_state, saved_fine_temperature, &
    saved_time, saved_steps, saved_integrals, saved_maximum_reflux, &
    "transport policy mismatch transaction", &
    saved_maximum_diffusivity=saved_transport_diffusivity, &
    saved_minimum_theta=saved_transport_theta, &
    maximum_diffusivity=transport_diffusivity, minimum_theta=transport_theta)

  incompatible_config = transport_config
  incompatible_config%transport_cfl = 0.175_dp
  call read_amr_reactive_3d_checkpoint( &
    transport_checkpoint_path, species, incompatible_config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message, &
    transport=transport, &
    maximum_transport_diffusivity=transport_diffusivity, &
    minimum_transport_theta=transport_theta)
  call require(.not. ok, "transport CFL mismatch rejection")
  call require_checkpoint_unchanged( &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, &
    saved_coarse_state, saved_coarse_temperature, &
    saved_fine_state, saved_fine_temperature, &
    saved_time, saved_steps, saved_integrals, saved_maximum_reflux, &
    "transport CFL mismatch transaction", &
    saved_maximum_diffusivity=saved_transport_diffusivity, &
    saved_minimum_theta=saved_transport_theta, &
    maximum_diffusivity=transport_diffusivity, minimum_theta=transport_theta)

  call delete_checkpoint(invalid_checkpoint_path)
  call write_amr_reactive_3d_checkpoint( &
    invalid_checkpoint_path, species, transport_config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message, &
    transport=transport, maximum_transport_diffusivity=transport_diffusivity)
  call require(.not. ok, "incomplete transport context rejection")
  call require_file_absent( &
    invalid_checkpoint_path, "incomplete transport context creates no checkpoint")

  incompatible_config = config
  incompatible_config%riemann_solver = "rusanov"
  call read_amr_reactive_3d_checkpoint( &
    checkpoint_path, species, incompatible_config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message)
  call require(.not. ok, "solver fingerprint mismatch rejection")
  call require_checkpoint_unchanged( &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, &
    saved_coarse_state, saved_coarse_temperature, &
    saved_fine_state, saved_fine_temperature, &
    saved_time, saved_steps, saved_integrals, saved_maximum_reflux, &
    "solver mismatch transaction")

  incompatible_config = config
  incompatible_config%reconstruction = "pcm"
  call read_amr_reactive_3d_checkpoint( &
    checkpoint_path, species, incompatible_config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message)
  call require(.not. ok, "reconstruction fingerprint mismatch rejection")
  call require_checkpoint_unchanged( &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, &
    saved_coarse_state, saved_coarse_temperature, &
    saved_fine_state, saved_fine_temperature, &
    saved_time, saved_steps, saved_integrals, saved_maximum_reflux, &
    "reconstruction mismatch transaction")

  selected_config = config
  selected_config%thermo_model = "selected"
  selected_config%chemistry_enabled = .true.
  bundle_sha256 = repeat("a", len(bundle_sha256))
  allocate(selected_composition(size(species)))
  selected_composition = 1.0_dp / real(size(species), dp)
  call write_amr_reactive_3d_checkpoint( &
    selected_checkpoint_path, species, selected_config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message, &
    bundle_sha256=bundle_sha256, chemistry_integrator="explicit", &
    base_mole_fractions=selected_composition, reactions=reactions)
  call require(ok, "selected checkpoint write: " // trim(message))

  coarse_state = -1.0_dp
  fine_state = -2.0_dp
  coarse_temperature = 300.0_dp
  fine_temperature = 400.0_dp
  initial_integrals = -3.0_dp
  time = 0.0_dp
  steps = 0
  maximum_reflux = 0.0_dp
  call read_amr_reactive_3d_checkpoint( &
    selected_checkpoint_path, species, selected_config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message, &
    bundle_sha256=bundle_sha256, chemistry_integrator="explicit", &
    base_mole_fractions=selected_composition, reactions=reactions)
  call require(ok, "selected checkpoint read: " // trim(message))
  call require_checkpoint_unchanged( &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, &
    saved_coarse_state, saved_coarse_temperature, &
    saved_fine_state, saved_fine_temperature, &
    saved_time, saved_steps, saved_integrals, saved_maximum_reflux, &
    "selected checkpoint exact round trip")

  call read_amr_reactive_3d_checkpoint( &
    selected_checkpoint_path, species, config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message)
  call require(.not. ok, "fixed reader rejects selected schema")
  call require_checkpoint_unchanged( &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, &
    saved_coarse_state, saved_coarse_temperature, &
    saved_fine_state, saved_fine_temperature, &
    saved_time, saved_steps, saved_integrals, saved_maximum_reflux, &
    "fixed reader selected-schema transaction")

  call read_amr_reactive_3d_checkpoint( &
    checkpoint_path, species, selected_config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message, &
    bundle_sha256=bundle_sha256, chemistry_integrator="explicit", &
    base_mole_fractions=selected_composition, reactions=reactions)
  call require(.not. ok, "selected reader rejects fixed schema")
  call require_checkpoint_unchanged( &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, &
    saved_coarse_state, saved_coarse_temperature, &
    saved_fine_state, saved_fine_temperature, &
    saved_time, saved_steps, saved_integrals, saved_maximum_reflux, &
    "selected reader fixed-schema transaction")

  incompatible_bundle_sha256 = bundle_sha256
  incompatible_bundle_sha256(64:64) = "b"
  call read_amr_reactive_3d_checkpoint( &
    selected_checkpoint_path, species, selected_config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message, &
    bundle_sha256=incompatible_bundle_sha256, &
    chemistry_integrator="explicit", &
    base_mole_fractions=selected_composition, reactions=reactions)
  call require(.not. ok .and. index(message, "bundle SHA-256 mismatch") > 0, &
    "selected bundle mismatch rejection")
  call require_checkpoint_unchanged( &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, &
    saved_coarse_state, saved_coarse_temperature, &
    saved_fine_state, saved_fine_temperature, &
    saved_time, saved_steps, saved_integrals, saved_maximum_reflux, &
    "selected bundle mismatch transaction")

  call read_amr_reactive_3d_checkpoint( &
    selected_checkpoint_path, species, selected_config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message, &
    bundle_sha256=bundle_sha256, chemistry_integrator="implicit", &
    base_mole_fractions=selected_composition, reactions=reactions)
  call require(.not. ok .and. &
    index(message, "chemistry integrator mismatch") > 0, &
    "selected integrator mismatch rejection")
  call require_checkpoint_unchanged( &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, &
    saved_coarse_state, saved_coarse_temperature, &
    saved_fine_state, saved_fine_temperature, &
    saved_time, saved_steps, saved_integrals, saved_maximum_reflux, &
    "selected integrator mismatch transaction")

  incompatible_composition = selected_composition
  incompatible_composition(1) = incompatible_composition(1) + 1.0e-6_dp
  incompatible_composition(2) = incompatible_composition(2) - 1.0e-6_dp
  call read_amr_reactive_3d_checkpoint( &
    selected_checkpoint_path, species, selected_config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message, &
    bundle_sha256=bundle_sha256, chemistry_integrator="explicit", &
    base_mole_fractions=incompatible_composition, reactions=reactions)
  call require(.not. ok .and. index(message, "composition mismatch") > 0, &
    "selected composition mismatch rejection")
  call require_checkpoint_unchanged( &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, &
    saved_coarse_state, saved_coarse_temperature, &
    saved_fine_state, saved_fine_temperature, &
    saved_time, saved_steps, saved_integrals, saved_maximum_reflux, &
    "selected composition mismatch transaction")

  incompatible_reactions = reactions
  incompatible_reactions(1)%forward_rate%pre_exponential = &
    1.000001_dp * &
      incompatible_reactions(1)%forward_rate%pre_exponential
  call read_amr_reactive_3d_checkpoint( &
    selected_checkpoint_path, species, selected_config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message, &
    bundle_sha256=bundle_sha256, chemistry_integrator="explicit", &
    base_mole_fractions=selected_composition, &
    reactions=incompatible_reactions)
  call require(.not. ok .and. index(message, "mechanism mismatch") > 0, &
    "selected mechanism mismatch rejection")
  call require_checkpoint_unchanged( &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, &
    saved_coarse_state, saved_coarse_temperature, &
    saved_fine_state, saved_fine_temperature, &
    saved_time, saved_steps, saved_integrals, saved_maximum_reflux, &
    "selected mechanism mismatch transaction")

  incompatible_config = selected_config
  incompatible_config%chemistry_relative_tolerance = &
    2.0_dp * selected_config%chemistry_relative_tolerance
  call read_amr_reactive_3d_checkpoint( &
    selected_checkpoint_path, species, incompatible_config, amr_config, &
    patch, coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message, &
    bundle_sha256=bundle_sha256, chemistry_integrator="explicit", &
    base_mole_fractions=selected_composition, reactions=reactions)
  call require(.not. ok .and. &
    index(message, "chemistry controls mismatch") > 0, &
    "selected chemistry-control mismatch rejection")
  call require_checkpoint_unchanged( &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, &
    saved_coarse_state, saved_coarse_temperature, &
    saved_fine_state, saved_fine_temperature, &
    saved_time, saved_steps, saved_integrals, saved_maximum_reflux, &
    "selected chemistry-control mismatch transaction")

  call write_amr_reactive_3d_checkpoint( &
    invalid_checkpoint_path, species, selected_config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message, &
    bundle_sha256=bundle_sha256)
  call require(.not. ok, "incomplete selected context rejection")
  call require_file_absent( &
    invalid_checkpoint_path, "incomplete context creates no checkpoint")

  incompatible_reactions = reactions
  allocate(incompatible_reactions(1)%third_body_efficiencies(1))
  incompatible_reactions(1)%third_body_efficiencies = 1.0_dp
  call write_amr_reactive_3d_checkpoint( &
    invalid_checkpoint_path, species, selected_config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message, &
    bundle_sha256=bundle_sha256, chemistry_integrator="explicit", &
    base_mole_fractions=selected_composition, &
    reactions=incompatible_reactions)
  call require(.not. ok, "malformed selected mechanism rejection")
  call require_file_absent( &
    invalid_checkpoint_path, "malformed context creates no checkpoint")

  selected_transport_config = selected_config
  selected_transport_config%transport_enabled = .true.
  call write_amr_reactive_3d_checkpoint( &
    selected_transport_checkpoint_path, species, selected_transport_config, &
    amr_config, patch, coarse_state, coarse_temperature, fine_state, &
    fine_temperature, time, steps, initial_integrals, maximum_reflux, ok, &
    message, bundle_sha256=bundle_sha256, chemistry_integrator="explicit", &
    base_mole_fractions=selected_composition, reactions=reactions, &
    transport=transport, maximum_transport_diffusivity=transport_diffusivity, &
    minimum_transport_theta=transport_theta)
  call require(ok, "selected transport checkpoint write: " // trim(message))
  call require_checkpoint_schema( &
    selected_transport_checkpoint_path, 5, "selected transport schema five")

  coarse_state = -1.0_dp
  fine_state = -2.0_dp
  coarse_temperature = 300.0_dp
  fine_temperature = 400.0_dp
  initial_integrals = -3.0_dp
  time = 0.0_dp
  steps = 0
  maximum_reflux = 0.0_dp
  transport_diffusivity = 0.0_dp
  transport_theta = 1.0_dp
  call read_amr_reactive_3d_checkpoint( &
    selected_transport_checkpoint_path, species, selected_transport_config, &
    amr_config, patch, coarse_state, coarse_temperature, fine_state, &
    fine_temperature, time, steps, initial_integrals, maximum_reflux, ok, &
    message, bundle_sha256=bundle_sha256, chemistry_integrator="explicit", &
    base_mole_fractions=selected_composition, reactions=reactions, &
    transport=transport, maximum_transport_diffusivity=transport_diffusivity, &
    minimum_transport_theta=transport_theta)
  call require(ok, "selected transport checkpoint read: " // trim(message))
  call require_checkpoint_unchanged( &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, &
    saved_coarse_state, saved_coarse_temperature, &
    saved_fine_state, saved_fine_temperature, &
    saved_time, saved_steps, saved_integrals, saved_maximum_reflux, &
    "selected transport checkpoint exact round trip", &
    saved_maximum_diffusivity=saved_transport_diffusivity, &
    saved_minimum_theta=saved_transport_theta, &
    maximum_diffusivity=transport_diffusivity, minimum_theta=transport_theta)

  call delete_checkpoint(invalid_checkpoint_path)
  call write_amr_reactive_3d_checkpoint( &
    invalid_checkpoint_path, species, selected_transport_config, amr_config, &
    patch, coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message, &
    bundle_sha256=bundle_sha256, chemistry_integrator="explicit", &
    base_mole_fractions=selected_composition, reactions=reactions, &
    transport=transport, maximum_transport_diffusivity=transport_diffusivity)
  call require(.not. ok, "incomplete selected transport context rejection")
  call require_file_absent( &
    invalid_checkpoint_path, &
    "incomplete selected transport context creates no checkpoint")

  call exercise_transport_corruption( &
    transport_checkpoint_path, "operator", .false.)
  call exercise_transport_corruption( &
    transport_checkpoint_path, "convention", .false.)
  call exercise_transport_corruption( &
    transport_checkpoint_path, "phase", .false.)
  call exercise_transport_corruption( &
    transport_checkpoint_path, "diagnostic_nan", .false.)
  call exercise_transport_corruption( &
    transport_checkpoint_path, "diagnostic_negative", .false.)
  call exercise_transport_corruption( &
    transport_checkpoint_path, "trailing", .false.)
  call exercise_transport_corruption( &
    selected_transport_checkpoint_path, "operator", .true.)
  call exercise_transport_corruption( &
    selected_transport_checkpoint_path, "convention", .true.)
  call exercise_transport_corruption( &
    selected_transport_checkpoint_path, "phase", .true.)
  call exercise_transport_corruption( &
    selected_transport_checkpoint_path, "diagnostic_nan", .true.)
  call exercise_transport_corruption( &
    selected_transport_checkpoint_path, "diagnostic_negative", .true.)
  call exercise_transport_corruption( &
    selected_transport_checkpoint_path, "trailing", .true.)

  call write_truncated_checkpoint(checkpoint_path)
  call read_amr_reactive_3d_checkpoint( &
    checkpoint_path, species, config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message)
  call require(.not. ok, "truncated checkpoint rejection")
  call require_checkpoint_unchanged( &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, &
    saved_coarse_state, saved_coarse_temperature, &
    saved_fine_state, saved_fine_temperature, &
    saved_time, saved_steps, saved_integrals, saved_maximum_reflux, &
    "truncated checkpoint transaction")

  coarse_state(1, 2, 2, 2) = coarse_state(1, 2, 2, 2) + 1.0_dp
  call write_amr_reactive_3d_checkpoint( &
    checkpoint_path, species, config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message)
  call require(.not. ok, "unsynchronized checkpoint write rejection")
  coarse_state = saved_coarse_state
  coarse_temperature(1, 1, 1) = coarse_temperature(1, 1, 1) + 10.0_dp
  call write_amr_reactive_3d_checkpoint( &
    checkpoint_path, species, config, amr_config, patch, &
    coarse_state, coarse_temperature, fine_state, fine_temperature, &
    time, steps, initial_integrals, maximum_reflux, ok, message)
  call require(.not. ok, "inconsistent checkpoint temperature rejection")
  call delete_checkpoint(checkpoint_path)
  call delete_checkpoint(selected_checkpoint_path)
  call delete_checkpoint(transport_checkpoint_path)
  call delete_checkpoint(selected_transport_checkpoint_path)
  call delete_checkpoint(invalid_checkpoint_path)
  write(*, '(a)') "test_amr_reactive_3d_checkpoint: PASS"

contains

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

  subroutine exercise_transport_corruption(path, corruption, selected_context)
    character(len=*), intent(in) :: path, corruption
    logical, intent(in) :: selected_context
    type(reactive_3d_config) :: test_config

    test_config = transport_config
    if (selected_context) test_config = selected_transport_config
    call delete_checkpoint(path)
    if (selected_context) then
      call write_amr_reactive_3d_checkpoint( &
        path, species, test_config, amr_config, patch, coarse_state, &
        coarse_temperature, fine_state, fine_temperature, time, steps, &
        initial_integrals, maximum_reflux, ok, message, &
        bundle_sha256=bundle_sha256, chemistry_integrator="explicit", &
        base_mole_fractions=selected_composition, reactions=reactions, &
        transport=transport, maximum_transport_diffusivity=transport_diffusivity, &
        minimum_transport_theta=transport_theta)
    else
      call write_amr_reactive_3d_checkpoint( &
        path, species, test_config, amr_config, patch, coarse_state, &
        coarse_temperature, fine_state, fine_temperature, time, steps, &
        initial_integrals, maximum_reflux, ok, message, transport=transport, &
        maximum_transport_diffusivity=transport_diffusivity, &
        minimum_transport_theta=transport_theta)
    end if
    call require(ok, "transport corruption fixture write: " // trim(corruption))

    select case (trim(corruption))
    case ("operator")
      call replace_checkpoint_line( &
        path, "STATIC_AMR_3D_R_T_H_T_R_V1", "BROKEN_OPERATOR")
    case ("convention")
      call replace_checkpoint_line( &
        path, "EPSILON_OVER_K_K;SIGMA_ANGSTROM;DIPOLE_DEBYE;" // &
        "POLARIZABILITY_ANGSTROM3;ROT_RELAX_DIMENSIONLESS", &
        "BROKEN_CONVENTION")
    case ("phase")
      call replace_checkpoint_line( &
        path, "POST_ACCEPTED_COARSE_STEP", "BROKEN_PHASE")
    case ("diagnostic_nan")
      call replace_checkpoint_line_after( &
        path, "TRANSPORT_DIAGNOSTICS", "NaN 0.37500000000000000e+00")
    case ("diagnostic_negative")
      call replace_checkpoint_line_after( &
        path, "TRANSPORT_DIAGNOSTICS", "-1.00000000000000000e+00 0.37500000000000000e+00")
    case ("trailing")
      call append_checkpoint_line(path, "TRAILING_DATA")
    case default
      error stop "Unknown transport checkpoint corruption"
    end select

    if (selected_context) then
      call read_amr_reactive_3d_checkpoint( &
        path, species, test_config, amr_config, patch, coarse_state, &
        coarse_temperature, fine_state, fine_temperature, time, steps, &
        initial_integrals, maximum_reflux, ok, message, &
        bundle_sha256=bundle_sha256, chemistry_integrator="explicit", &
        base_mole_fractions=selected_composition, reactions=reactions, &
        transport=transport, maximum_transport_diffusivity=transport_diffusivity, &
        minimum_transport_theta=transport_theta)
    else
      call read_amr_reactive_3d_checkpoint( &
        path, species, test_config, amr_config, patch, coarse_state, &
        coarse_temperature, fine_state, fine_temperature, time, steps, &
        initial_integrals, maximum_reflux, ok, message, transport=transport, &
        maximum_transport_diffusivity=transport_diffusivity, &
        minimum_transport_theta=transport_theta)
    end if
    call require(.not. ok, &
      "transport " // trim(corruption) // " corruption rejection")
    call require_checkpoint_unchanged( &
      coarse_state, coarse_temperature, fine_state, fine_temperature, &
      time, steps, initial_integrals, maximum_reflux, &
      saved_coarse_state, saved_coarse_temperature, &
      saved_fine_state, saved_fine_temperature, &
      saved_time, saved_steps, saved_integrals, saved_maximum_reflux, &
      "transport " // trim(corruption) // " corruption transaction", &
      saved_maximum_diffusivity=saved_transport_diffusivity, &
      saved_minimum_theta=saved_transport_theta, &
      maximum_diffusivity=transport_diffusivity, minimum_theta=transport_theta)
  end subroutine exercise_transport_corruption

  subroutine replace_checkpoint_line(path, target, replacement)
    character(len=*), intent(in) :: path, target, replacement
    character(len=4096), allocatable :: lines(:)
    character(len=4096) :: line
    integer :: unit, status, line_count, index, matches

    allocate(lines(2048))
    line_count = 0
    open(newunit=unit, file=trim(path), status="old", action="read", &
      form="formatted", iostat=status)
    call require(status == 0, "open checkpoint for corruption")
    do
      read(unit, '(a)', iostat=status) line
      if (is_iostat_end(status)) exit
      call require(status == 0, "read checkpoint for corruption")
      line_count = line_count + 1
      call require(line_count <= size(lines), "checkpoint corruption buffer size")
      lines(line_count) = line
    end do
    close(unit, iostat=status)
    call require(status == 0, "close checkpoint after corruption read")
    matches = 0
    do index = 1, line_count
      if (trim(lines(index)) == trim(target)) then
        lines(index) = trim(replacement)
        matches = matches + 1
      end if
    end do
    call require(matches == 1, "checkpoint corruption target exists once")
    call rewrite_checkpoint(path, lines, line_count)
    deallocate(lines)
  end subroutine replace_checkpoint_line

  subroutine replace_checkpoint_line_after(path, marker, replacement)
    character(len=*), intent(in) :: path, marker, replacement
    character(len=4096), allocatable :: lines(:)
    character(len=4096) :: line
    integer :: unit, status, line_count, index, marker_index

    allocate(lines(2048))
    line_count = 0
    open(newunit=unit, file=trim(path), status="old", action="read", &
      form="formatted", iostat=status)
    call require(status == 0, "open checkpoint for diagnostic corruption")
    do
      read(unit, '(a)', iostat=status) line
      if (is_iostat_end(status)) exit
      call require(status == 0, "read checkpoint for diagnostic corruption")
      line_count = line_count + 1
      call require(line_count <= size(lines), "checkpoint diagnostic buffer size")
      lines(line_count) = line
    end do
    close(unit, iostat=status)
    call require(status == 0, "close checkpoint after diagnostic read")
    marker_index = 0
    do index = 1, line_count
      if (trim(lines(index)) == trim(marker)) then
        marker_index = index
        exit
      end if
    end do
    call require(marker_index > 0 .and. marker_index < line_count, &
      "checkpoint diagnostic marker exists")
    lines(marker_index + 1) = trim(replacement)
    call rewrite_checkpoint(path, lines, line_count)
    deallocate(lines)
  end subroutine replace_checkpoint_line_after

  subroutine rewrite_checkpoint(path, lines, line_count)
    character(len=*), intent(in) :: path
    character(len=*), intent(in) :: lines(:)
    integer, intent(in) :: line_count
    integer :: unit, status, index

    open(newunit=unit, file=trim(path), status="replace", action="write", &
      form="formatted", iostat=status)
    call require(status == 0, "open checkpoint for corruption rewrite")
    do index = 1, line_count
      write(unit, '(a)', iostat=status) trim(lines(index))
      call require(status == 0, "write checkpoint corruption rewrite")
    end do
    close(unit, iostat=status)
    call require(status == 0, "close checkpoint corruption rewrite")
  end subroutine rewrite_checkpoint

  subroutine append_checkpoint_line(path, line)
    character(len=*), intent(in) :: path, line
    integer :: unit, status

    open(newunit=unit, file=trim(path), status="old", position="append", &
      action="write", form="formatted", iostat=status)
    call require(status == 0, "open checkpoint for trailing-data corruption")
    write(unit, '(a)', iostat=status) trim(line)
    call require(status == 0, "write checkpoint trailing data")
    close(unit, iostat=status)
    call require(status == 0, "close checkpoint trailing-data corruption")
  end subroutine append_checkpoint_line

  subroutine require_checkpoint_unchanged( &
      coarse_state, coarse_temperature, fine_state, fine_temperature, &
      time, steps, integrals, maximum_reflux, &
      saved_coarse_state, saved_coarse_temperature, &
      saved_fine_state, saved_fine_temperature, &
      saved_time, saved_steps, saved_integrals, saved_maximum_reflux, label, &
      saved_maximum_diffusivity, saved_minimum_theta, maximum_diffusivity, &
      minimum_theta)
    real(dp), intent(in) :: coarse_state(:, :, :, :)
    real(dp), intent(in) :: coarse_temperature(:, :, :)
    real(dp), intent(in) :: fine_state(:, :, :, :)
    real(dp), intent(in) :: fine_temperature(:, :, :)
    real(dp), intent(in) :: time, integrals(:), maximum_reflux
    real(dp), intent(in) :: saved_coarse_state(:, :, :, :)
    real(dp), intent(in) :: saved_coarse_temperature(:, :, :)
    real(dp), intent(in) :: saved_fine_state(:, :, :, :)
    real(dp), intent(in) :: saved_fine_temperature(:, :, :)
    real(dp), intent(in) :: saved_time, saved_integrals(:)
    real(dp), intent(in) :: saved_maximum_reflux
    integer, intent(in) :: steps, saved_steps
    character(len=*), intent(in) :: label
    real(dp), intent(in), optional :: saved_maximum_diffusivity
    real(dp), intent(in), optional :: saved_minimum_theta
    real(dp), intent(in), optional :: maximum_diffusivity
    real(dp), intent(in), optional :: minimum_theta

    call require( &
      maxval(abs(coarse_state - saved_coarse_state)) == 0.0_dp .and. &
      maxval(abs(coarse_temperature - saved_coarse_temperature)) == &
        0.0_dp .and. &
      maxval(abs(fine_state - saved_fine_state)) == 0.0_dp .and. &
      maxval(abs(fine_temperature - saved_fine_temperature)) == &
        0.0_dp .and. &
      time == saved_time .and. steps == saved_steps .and. &
      all(integrals == saved_integrals) .and. &
      maximum_reflux == saved_maximum_reflux, label)
    if (present(saved_maximum_diffusivity) .or. &
        present(saved_minimum_theta) .or. present(maximum_diffusivity) .or. &
        present(minimum_theta)) then
      call require(present(saved_maximum_diffusivity) .and. &
        present(saved_minimum_theta) .and. present(maximum_diffusivity) .and. &
        present(minimum_theta), label // " transport diagnostic arguments")
      call require(maximum_diffusivity == saved_maximum_diffusivity .and. &
        minimum_theta == saved_minimum_theta, label // &
        " transport diagnostics unchanged")
    end if
  end subroutine require_checkpoint_unchanged

  subroutine require_checkpoint_schema(path, expected_schema, label)
    character(len=*), intent(in) :: path, label
    integer, intent(in) :: expected_schema
    character(len=256) :: magic
    integer :: unit, status, schema, species_count, variable_count

    open(newunit=unit, file=trim(path), status="old", action="read", &
      iostat=status)
    call require(status == 0, label // " opens checkpoint")
    read(unit, '(a)', iostat=status) magic
    call require(status == 0 .and. trim(magic) == &
      "PELEF_AMR_REACTIVE_3D_CHECKPOINT", label // " magic")
    read(unit, *, iostat=status) schema, species_count, variable_count
    call require(status == 0 .and. schema == expected_schema .and. &
      species_count > 0 .and. variable_count == species_count + 5, &
      label // " header")
    close(unit, iostat=status)
    call require(status == 0, label // " closes checkpoint")
  end subroutine require_checkpoint_schema

  subroutine write_truncated_checkpoint(path)
    character(len=*), intent(in) :: path
    integer :: unit, status

    open(newunit=unit, file=trim(path), status="replace", action="write", &
      iostat=status)
    if (status /= 0) error stop "Could not create truncated checkpoint"
    write(unit, '(a)', iostat=status) &
      "PELEF_AMR_REACTIVE_3D_CHECKPOINT"
    if (status /= 0) error stop "Could not write truncated checkpoint"
    close(unit, iostat=status)
    if (status /= 0) error stop "Could not close truncated checkpoint"
  end subroutine write_truncated_checkpoint

  subroutine delete_checkpoint(path)
    character(len=*), intent(in) :: path
    logical :: exists
    integer :: unit, status

    inquire(file=trim(path), exist=exists)
    if (.not. exists) return
    open(newunit=unit, file=trim(path), status="old", action="read", &
      iostat=status)
    if (status /= 0) error stop "Could not open checkpoint for deletion"
    close(unit, status="delete", iostat=status)
    if (status /= 0) error stop "Could not delete checkpoint"
  end subroutine delete_checkpoint

  subroutine require_file_absent(path, label)
    character(len=*), intent(in) :: path, label
    logical :: exists

    inquire(file=trim(path), exist=exists)
    call require(.not. exists, label)
  end subroutine require_file_absent

  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      write(*, '(a)') "FAIL: " // trim(label)
      error stop 1
    end if
  end subroutine require

end program test_amr_reactive_3d_checkpoint
