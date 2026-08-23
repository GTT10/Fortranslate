program test_amr_multipatch_reactive_physics_1d
  use precision_mod, only: dp
  use state_indices_mod, only: irho, iet
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_conserved_to_primitive
  use amr_hierarchy_1d_mod, only: restrict_average_1d
  use amr_multipatch_reactive_1d_mod, only: &
    amr_multipatch_reactive_solution_1d, &
    initialize_multipatch_reactive_1d, &
    multipatch_reactive_timestep_1d, advance_multipatch_reactive_1d, &
    multipatch_reactive_integrals_1d
  implicit none

  integer, parameter :: patch_count = 2
  integer, parameter :: patch_lower(patch_count) = [3, 15]
  integer, parameter :: patch_upper(patch_count) = [8, 20]
  real(dp), parameter :: conservation_tolerance = 2.0e-9_dp
  real(dp), parameter :: synchronization_tolerance = 8.0e-13_dp
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(reactive_1d_config) :: config
  type(amr_multipatch_reactive_solution_1d) :: solution
  real(dp), allocatable :: initial_integral(:), final_integral(:)
  real(dp), allocatable :: restricted(:, :), q(:)
  real(dp), allocatable :: initial_patch_state(:, :)
  real(dp) :: dt, conservation_error, synchronization_error
  real(dp) :: local_temperature, sound_speed, closure
  logical :: ok
  integer :: patch, cell, fine_cells, covered_cells, species_index

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "multipatch physics thermodynamics load")
  call load_h2o2_elementary_mechanism(reactions, ok)
  call require(ok, "multipatch physics mechanism load")
  call load_h2o2_elementary_transport(transport, ok)
  call require(ok, "multipatch physics transport load")
  call configure_case(config)
  call initialize_multipatch_reactive_1d( &
    species, config, patch_lower, patch_upper, solution, ok)
  call require(ok .and. solution%is_valid(), &
    "multipatch reacting transport initialization")

  allocate(initial_integral(reactive_nvar(size(species))))
  allocate(final_integral(reactive_nvar(size(species))))
  call multipatch_reactive_integrals_1d(solution, initial_integral, ok)
  call require(ok, "initial multipatch physics integral")
  fine_cells = solution%hierarchy%patches(1)%fine%cell_count()
  initial_patch_state = solution%patches(1)%state(:, 1:fine_cells)

  call multipatch_reactive_timestep_1d( &
    species, config, solution, dt, ok)
  call require(.not. ok, "transport timestep requires its database")
  call multipatch_reactive_timestep_1d( &
    species, config, solution, dt, ok, transport)
  call require(ok .and. dt > 0.0_dp, &
    "multipatch hydro-transport stable timestep")
  dt = min(dt, 1.0e-10_dp)
  call advance_multipatch_reactive_1d( &
    species, reactions, config, dt, solution, ok, transport)
  call require(ok, "multipatch chemistry-transport split step")
  call require(solution%steps == 1 .and. &
    abs(solution%time - dt) <= 32.0_dp * epsilon(1.0_dp) * dt, &
    "multipatch physics time and step accounting")
  call require(maxval(abs(solution%patches(1)%state(:, 1:fine_cells) - &
    initial_patch_state)) > 1.0e-12_dp, &
    "multipatch physics advances the first fine patch")

  call multipatch_reactive_integrals_1d(solution, final_integral, ok)
  call require(ok, "final multipatch physics integral")
  conservation_error = maxval(abs( &
    final_integral(irho:iet) - initial_integral(irho:iet)) / &
    max(1.0_dp, abs(initial_integral(irho:iet))))
  call require(conservation_error < conservation_tolerance, &
    "multipatch reacting transport conserved hydro quantities")

  synchronization_error = 0.0_dp
  do patch = 1, patch_count
    fine_cells = solution%hierarchy%patches(patch)%fine%cell_count()
    covered_cells = solution%hierarchy%patches(patch)%covered_coarse_cells()
    if (allocated(restricted)) deallocate(restricted)
    allocate(restricted(size(solution%coarse, 1), covered_cells))
    call restrict_average_1d( &
      solution%patches(patch)%state(:, 1:fine_cells), &
      solution%hierarchy%patches(patch), restricted, ok)
    call require(ok, "multipatch physics restriction")
    synchronization_error = max(synchronization_error, maxval(abs( &
      restricted - solution%coarse(:, &
        solution%hierarchy%patches(patch)%fine_coarse_lower: &
        solution%hierarchy%patches(patch)%fine_coarse_upper))) / &
      max(1.0_dp, maxval(abs(restricted))))
  end do
  call require(synchronization_error < synchronization_tolerance, &
    "chemistry and transport leave every patch synchronized")

  allocate(q(reactive_nprim(size(species))))
  do cell = 1, config%nx
    call inspect_cell( &
      solution%coarse(:, cell), solution%coarse_temperature(cell))
  end do
  do patch = 1, patch_count
    fine_cells = solution%hierarchy%patches(patch)%fine%cell_count()
    do cell = 1, fine_cells
      call inspect_cell( &
        solution%patches(patch)%state(:, cell), &
        solution%patches(patch)%temperature(cell))
    end do
  end do

  initial_integral = final_integral
  call advance_multipatch_reactive_1d( &
    species, reactions, config, dt, solution, ok)
  call require(.not. ok, "transport advance requires its database")
  call multipatch_reactive_integrals_1d(solution, final_integral, ok)
  call require(ok .and. maxval(abs(final_integral - initial_integral)) == &
    0.0_dp .and. solution%steps == 1, &
    "rejected multipatch physics step is transactional")

  write(*, '(a,1x,es16.8)') &
    "Multipatch reacting transport conservation error:", conservation_error
  write(*, '(a,1x,es16.8)') &
    "Multipatch reacting transport synchronization error:", &
    synchronization_error
  write(*, '(a)') "test_amr_multipatch_reactive_physics_1d: PASS"

contains

  subroutine configure_case(local_config)
    type(reactive_1d_config), intent(out) :: local_config

    local_config = reactive_1d_config()
    local_config%nx = 24
    local_config%x_lower = 0.0_dp
    local_config%x_upper = 0.012_dp
    local_config%cfl = 0.20_dp
    local_config%problem = "reactive_hotspot"
    local_config%riemann_solver = "rusanov"
    local_config%limiter = "mc"
    local_config%boundary_condition = "periodic"
    local_config%chemistry_enabled = .true.
    local_config%transport_enabled = .true.
    local_config%viscosity_enabled = .true.
    local_config%thermal_conduction_enabled = .true.
    local_config%species_diffusion_enabled = .true.
    local_config%barodiffusion_enabled = .true.
    local_config%transport_cfl = 0.30_dp
    local_config%initial_temperature = 1200.0_dp
    local_config%initial_pressure = 101325.0_dp
    local_config%initial_velocity = 0.0_dp
    local_config%hotspot_temperature_rise = 200.0_dp
    local_config%hotspot_center = 0.006_dp
    local_config%hotspot_width = 0.0012_dp
    local_config%amr_enabled = .true.
    local_config%amr_reconstruction = "characteristic_ppm"
    local_config%amr_hybrid_weno = .true.
    local_config%amr_weno_scheme = 2
    local_config%amr_refinement_ratio = 2
  end subroutine configure_case

  subroutine inspect_cell(state, temperature)
    real(dp), intent(in) :: state(:), temperature

    call reactive_conserved_to_primitive( &
      species, state, temperature, q, local_temperature, sound_speed, ok)
    call require(ok .and. local_temperature > 0.0_dp .and. &
      sound_speed > 0.0_dp, "valid multipatch physics cell")
    closure = 0.0_dp
    do species_index = 1, size(species)
      call require(state(reactive_species_component(species_index)) >= &
        -1.0e-13_dp * max(1.0_dp, state(irho)), &
        "multipatch physics nonnegative species")
      closure = closure + &
        state(reactive_species_component(species_index)) / state(irho)
    end do
    call require(abs(closure - 1.0_dp) < 3.0e-10_dp, &
      "multipatch physics species closure")
  end subroutine inspect_cell

  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      write(*, '(a,1x,a)') "FAIL:", trim(label)
      error stop 1
    end if
  end subroutine require

end program test_amr_multipatch_reactive_physics_1d
