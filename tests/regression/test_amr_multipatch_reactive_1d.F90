program test_amr_multipatch_reactive_1d
  use precision_mod, only: dp
  use state_indices_mod, only: irho
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_conserved_to_primitive
  use amr_hierarchy_1d_mod, only: restrict_average_1d
  use amr_multipatch_reactive_1d_mod, only: &
    amr_multipatch_reactive_solution_1d, &
    initialize_multipatch_reactive_1d, &
    multipatch_reactive_timestep_1d, &
    advance_multipatch_reactive_hydro_1d, &
    multipatch_reactive_integrals_1d
  implicit none

  integer, parameter :: patch_count = 2
  integer, parameter :: patch_lower(patch_count) = [8, 40]
  integer, parameter :: patch_upper(patch_count) = [20, 52]
  real(dp), parameter :: conservation_tolerance = 3.0e-10_dp
  real(dp), parameter :: synchronization_tolerance = 5.0e-13_dp
  type(nasa7_species), allocatable :: species(:)
  type(reactive_1d_config) :: config
  type(amr_multipatch_reactive_solution_1d) :: solution
  real(dp), allocatable :: initial_integral(:), final_integral(:)
  real(dp), allocatable :: patch_one_initial(:, :)
  real(dp), allocatable :: patch_two_initial(:, :)
  real(dp), allocatable :: restricted(:, :), q(:)
  real(dp) :: dt, conservation_error, synchronization_error
  real(dp) :: minimum_temperature, maximum_closure_error
  real(dp) :: local_temperature, sound_speed, closure
  logical :: ok
  integer :: patch, cell, fine_cells, covered_cells, species_index, step

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "multipatch thermodynamics load")
  call configure_case(config)
  call initialize_multipatch_reactive_1d( &
    species, config, patch_lower, patch_upper, solution, ok)
  call require(ok .and. solution%is_valid(), &
    "multipatch reactive initialization")
  call require(solution%patch_count() == patch_count, &
    "reactive solution owns both patches")

  allocate(initial_integral(reactive_nvar(size(species))))
  allocate(final_integral(reactive_nvar(size(species))))
  call multipatch_reactive_integrals_1d(solution, initial_integral, ok)
  call require(ok, "initial multipatch reactive integral")
  fine_cells = solution%hierarchy%patches(1)%fine%cell_count()
  patch_one_initial = solution%patches(1)%state(:, 1:fine_cells)
  fine_cells = solution%hierarchy%patches(2)%fine%cell_count()
  patch_two_initial = solution%patches(2)%state(:, 1:fine_cells)
  call multipatch_reactive_timestep_1d( &
    species, config, solution, dt, ok)
  call require(ok .and. dt > 0.0_dp, "multipatch stable time step")
  dt = min(dt, 2.0e-8_dp)

  do step = 1, 2
    call advance_multipatch_reactive_hydro_1d( &
      species, config, dt, solution, ok)
    call require(ok, "multipatch WENO7-Z hydro step")
  end do
  call require(solution%steps == 2 .and. &
    abs(solution%time - 2.0_dp * dt) <= &
      32.0_dp * epsilon(1.0_dp) * dt, &
    "multipatch time and step accounting")
  fine_cells = solution%hierarchy%patches(1)%fine%cell_count()
  call require(maxval(abs(solution%patches(1)%state(:, 1:fine_cells) - &
    patch_one_initial)) > 1.0e-12_dp, &
    "first fine patch advanced")
  fine_cells = solution%hierarchy%patches(2)%fine%cell_count()
  call require(maxval(abs(solution%patches(2)%state(:, 1:fine_cells) - &
    patch_two_initial)) > 1.0e-12_dp, &
    "second fine patch advanced")

  call multipatch_reactive_integrals_1d(solution, final_integral, ok)
  call require(ok, "final multipatch reactive integral")
  conservation_error = maxval(abs(final_integral - initial_integral) / &
    max(1.0_dp, abs(initial_integral)))
  call require(conservation_error < conservation_tolerance, &
    "multipatch WENO7-Z reflux conservation")

  synchronization_error = 0.0_dp
  do patch = 1, patch_count
    fine_cells = solution%hierarchy%patches(patch)%fine%cell_count()
    covered_cells = solution%hierarchy%patches(patch)%covered_coarse_cells()
    if (allocated(restricted)) deallocate(restricted)
    allocate(restricted(size(solution%coarse, 1), covered_cells))
    call restrict_average_1d( &
      solution%patches(patch)%state(:, 1:fine_cells), &
      solution%hierarchy%patches(patch), restricted, ok)
    call require(ok, "multipatch final restriction")
    synchronization_error = max(synchronization_error, maxval(abs( &
      restricted - solution%coarse(:, &
        solution%hierarchy%patches(patch)%fine_coarse_lower: &
        solution%hierarchy%patches(patch)%fine_coarse_upper))) / &
      max(1.0_dp, maxval(abs(restricted))))
  end do
  call require(synchronization_error < synchronization_tolerance, &
    "both fine patches synchronized")

  allocate(q(reactive_nprim(size(species))))
  minimum_temperature = huge(1.0_dp)
  maximum_closure_error = 0.0_dp
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
  call require(minimum_temperature > 0.0_dp, &
    "multipatch thermodynamic positivity")
  call require(maximum_closure_error < 3.0e-10_dp, &
    "multipatch species closure")

  write(*, '(a,1x,es16.8)') &
    "Multipatch AMR conservation error:", conservation_error
  write(*, '(a,1x,es16.8)') &
    "Multipatch AMR synchronization error:", synchronization_error
  write(*, '(a)') "test_amr_multipatch_reactive_1d: PASS"

contains

  subroutine configure_case(local_config)
    type(reactive_1d_config), intent(out) :: local_config

    local_config = reactive_1d_config()
    local_config%nx = 64
    local_config%x_lower = 0.0_dp
    local_config%x_upper = 0.012_dp
    local_config%cfl = 0.20_dp
    local_config%problem = "entropy_wave"
    local_config%riemann_solver = "rusanov"
    local_config%limiter = "mc"
    local_config%boundary_condition = "periodic"
    local_config%chemistry_enabled = .false.
    local_config%transport_enabled = .false.
    local_config%initial_temperature = 1200.0_dp
    local_config%initial_pressure = 101325.0_dp
    local_config%initial_velocity = 25.0_dp
    local_config%density_wave_amplitude = 0.08_dp
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
      sound_speed > 0.0_dp, "valid multipatch reactive cell")
    minimum_temperature = min(minimum_temperature, local_temperature)
    closure = 0.0_dp
    do species_index = 1, size(species)
      closure = closure + &
        state(reactive_species_component(species_index)) / state(irho)
    end do
    maximum_closure_error = max(maximum_closure_error, abs(closure - 1.0_dp))
  end subroutine inspect_cell

  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      write(*, '(a,1x,a)') "FAIL:", trim(label)
      error stop 1
    end if
  end subroutine require

end program test_amr_multipatch_reactive_1d
