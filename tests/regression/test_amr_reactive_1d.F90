program test_amr_reactive_1d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: &
    reactive_nprim, reactive_mass_fraction_component, &
    reactive_conserved_to_primitive
  use amr_hierarchy_1d_mod, only: restrict_average_1d
  use amr_reactive_1d_mod, only: &
    amr_reactive_solution_1d, simulate_amr_reactive_1d
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(reactive_1d_config) :: config
  type(amr_reactive_solution_1d) :: solution
  real(dp), allocatable :: restricted(:, :), q(:)
  real(dp) :: initial_integrals(5), final_integrals(5)
  real(dp) :: conservation_error, synchronization_error
  real(dp) :: minimum_temperature, maximum_temperature
  logical :: ok
  integer :: cell, fine_cells, covered_cells

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "Failed to load AMR thermodynamics"
  call load_h2o2_elementary_mechanism(reactions, ok)
  if (.not. ok) error stop "Failed to load AMR mechanism"

  config = reactive_1d_config()
  config%nx = 24
  config%x_lower = 0.0_dp
  config%x_upper = 0.012_dp
  config%final_time = 5.0e-7_dp
  config%cfl = 0.25_dp
  config%maximum_steps = 1000
  config%problem = "reactive_hotspot"
  config%reconstruction = "pcm"
  config%riemann_solver = "rusanov"
  config%boundary_condition = "periodic"
  config%chemistry_enabled = .true.
  config%chemistry_relative_tolerance = 2.0e-7_dp
  config%chemistry_absolute_tolerance = 1.0e-12_dp
  config%initial_temperature = 1200.0_dp
  config%initial_pressure = 101325.0_dp
  config%initial_velocity = 0.0_dp
  config%hotspot_temperature_rise = 250.0_dp
  config%hotspot_center = 0.006_dp
  config%hotspot_width = 0.0012_dp
  config%amr_enabled = .true.
  config%amr_reconstruction = "plm"
  config%amr_refinement_ratio = 2
  config%amr_regrid_interval = 1
  config%amr_tag_component = 1
  config%amr_buffer_cells = 1
  config%amr_minimum_patch_cells = 6
  config%amr_relative_gradient_threshold = 0.005_dp
  config%amr_absolute_gradient_threshold = 0.0_dp
  config%amr_scale_floor = 1.0e-12_dp

  call simulate_amr_reactive_1d( &
    species, reactions, config, solution, initial_integrals, &
    final_integrals, ok)
  if (.not. ok) error stop "AMR reactive simulation failed"
  call assert_true(solution%steps >= 2, "multiple coarse steps")
  call assert_true(solution%regrid_evaluations >= 2, &
    "dynamic regrid evaluations")
  call assert_true(solution%regrids >= 1, "fine hierarchy created")
  call assert_true(solution%fine_active(), "fine hierarchy remains active")
  call assert_true(solution%hierarchy%fine_coarse_lower > 1 .and. &
    solution%hierarchy%fine_coarse_upper < config%nx, &
    "fine patch remains strictly interior")

  conservation_error = maxval(abs(final_integrals - initial_integrals) / &
    max(1.0_dp, abs(initial_integrals)))
  call assert_true(conservation_error < 2.0e-10_dp, &
    "composite mass momentum energy conservation")

  fine_cells = solution%hierarchy%fine%cell_count()
  covered_cells = solution%hierarchy%covered_coarse_cells()
  allocate(restricted(size(solution%coarse, 1), covered_cells))
  call restrict_average_1d( &
    solution%fine(:, 1:fine_cells), solution%hierarchy, restricted, ok)
  if (.not. ok) error stop "AMR final restriction failed"
  synchronization_error = maxval(abs( &
    restricted - solution%coarse(:, &
      solution%hierarchy%fine_coarse_lower: &
      solution%hierarchy%fine_coarse_upper))) / &
    max(1.0_dp, maxval(abs(restricted)))
  call assert_true(synchronization_error < 5.0e-13_dp, &
    "covered coarse state synchronized")

  allocate(q(reactive_nprim(size(species))))
  minimum_temperature = huge(1.0_dp)
  maximum_temperature = 0.0_dp
  do cell = 1, solution%hierarchy%fine_coarse_lower - 1
    call inspect_cell(solution%coarse(:, cell), &
      solution%coarse_temperature(cell))
  end do
  do cell = 1, fine_cells
    call inspect_cell(solution%fine(:, cell), solution%fine_temperature(cell))
  end do
  do cell = solution%hierarchy%fine_coarse_upper + 1, config%nx
    call inspect_cell(solution%coarse(:, cell), &
      solution%coarse_temperature(cell))
  end do
  call assert_true(minimum_temperature > 0.0_dp, &
    "composite thermodynamic positivity")
  call assert_true(maximum_temperature - minimum_temperature > 20.0_dp, &
    "hotspot structure retained")

  write(*, '(a,1x,es16.8)') "AMR conservation error:", conservation_error
  write(*, '(a,1x,es16.8)') &
    "AMR synchronization error:", synchronization_error
  write(*, '(a)') "test_amr_reactive_1d: PASS"

contains

  subroutine inspect_cell(state, temperature_guess)
    real(dp), intent(in) :: state(:), temperature_guess
    real(dp) :: temperature, sound_speed, closure
    logical :: local_ok
    integer :: k

    call reactive_conserved_to_primitive( &
      species, state, temperature_guess, q, temperature, sound_speed, &
      local_ok)
    if (.not. local_ok) error stop "Invalid AMR composite cell"
    closure = 0.0_dp
    do k = 1, size(species)
      closure = closure + q(reactive_mass_fraction_component(k))
    end do
    if (abs(closure - 1.0_dp) > 2.0e-10_dp) &
      error stop "AMR species closure failure"
    minimum_temperature = min(minimum_temperature, temperature)
    maximum_temperature = max(maximum_temperature, temperature)
  end subroutine inspect_cell

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*, '(a,1x,a)') "FAIL:", trim(label)
      error stop 1
    end if
  end subroutine assert_true

end program test_amr_reactive_1d
