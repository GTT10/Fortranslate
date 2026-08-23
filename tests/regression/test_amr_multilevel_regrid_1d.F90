program test_amr_multilevel_regrid_1d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: reactive_nvar
  use amr_multilevel_reactive_1d_mod, only: &
    amr_multilevel_reactive_solution_1d, &
    initialize_tagged_multilevel_reactive_1d, &
    regrid_multilevel_reactive_1d, simulate_multilevel_reactive_1d, &
    multilevel_reactive_integrals_1d, write_multilevel_reactive_1d_csv
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(reactive_1d_config) :: config
  type(amr_multilevel_reactive_solution_1d) :: solution, old_solution
  type(amr_multilevel_reactive_solution_1d) :: simulated
  real(dp), allocatable :: before(:), after(:)
  real(dp) :: initial_integrals(5), final_integrals(5)
  real(dp) :: regrid_error, simulation_error, overlap_error
  real(dp) :: old_lower, old_upper, new_lower, new_upper, dx
  logical :: ok, changed
  integer :: target_cell, lower, upper, old_first, new_first, overlap_cells

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamics load")
  call load_h2o2_elementary_mechanism(reactions, ok)
  call require(ok, "mechanism load")
  call configure_case(config)
  call initialize_tagged_multilevel_reactive_1d( &
    species, config, solution, ok)
  call require(ok, "tag-driven multilevel initialization")
  call require(solution%level_count() == 3, &
    "tag-driven initialization reaches configured depth")
  call require(solution%regrid_evaluations == 1, &
    "initial hierarchy evaluation counted")

  lower = solution%hierarchy%interfaces(1)%fine_coarse_lower
  upper = solution%hierarchy%interfaces(1)%fine_coarse_upper
  target_cell = max(3, lower - 2)
  if (target_cell >= lower) target_cell = min(config%nx - 2, upper + 2)
  call require(target_cell < lower .or. target_cell > upper, &
    "uncovered root cell available for hierarchy movement")
  solution%levels(1)%state(:, target_cell) = &
    1.20_dp * solution%levels(1)%state(:, target_cell)
  allocate(before(reactive_nvar(size(species))))
  allocate(after(reactive_nvar(size(species))))
  call multilevel_reactive_integrals_1d(solution, before, ok)
  call require(ok, "pre-regrid composite integral")
  old_solution = solution
  call regrid_multilevel_reactive_1d( &
    species, config, solution, changed, ok)
  call require(ok .and. changed, "tag-driven hierarchy changes")
  call require(solution%regrid_evaluations == 2, &
    "dynamic hierarchy evaluation counted")
  call require(solution%regrids >= 2, "dynamic hierarchy change counted")
  call require(solution%level_count() == 3, &
    "changed hierarchy retains overlap test depth")
  call old_solution%hierarchy%level_bounds( &
    2, old_lower, old_upper, ok)
  call require(ok, "old deepest bounds")
  call solution%hierarchy%level_bounds( &
    2, new_lower, new_upper, ok)
  call require(ok, "new deepest bounds")
  dx = solution%hierarchy%level_dx(2)
  old_first = nint((max(old_lower, new_lower) - old_lower) / dx) + 1
  new_first = nint((max(old_lower, new_lower) - new_lower) / dx) + 1
  overlap_cells = nint((min(old_upper, new_upper) - &
    max(old_lower, new_lower)) / dx)
  call require(overlap_cells > 0, "deepest old/new overlap exists")
  overlap_error = maxval(abs( &
    solution%levels(3)%state(:, &
      new_first:new_first + overlap_cells - 1) - &
    old_solution%levels(3)%state(:, &
      old_first:old_first + overlap_cells - 1)))
  call require(overlap_error == 0.0_dp, &
    "deepest overlapping fine state retained exactly")
  call require(solution%overlap_cells_transferred >= overlap_cells, &
    "overlap transfer count includes deepest cells")
  call multilevel_reactive_integrals_1d(solution, after, ok)
  call require(ok, "post-regrid composite integral")
  regrid_error = maxval(abs(after - before) / max(1.0_dp, abs(before)))
  call require(regrid_error < 5.0e-13_dp, &
    "multilevel hierarchy rebuild conserves the composite state")

  call simulate_multilevel_reactive_1d( &
    species, reactions, config, simulated, initial_integrals, &
    final_integrals, ok)
  call require(ok, "tag-driven multilevel simulation")
  call require(simulated%level_count() == 3, &
    "simulation retains three active levels")
  call require(simulated%regrid_evaluations >= 2, &
    "simulation performs dynamic regrid evaluation")
  simulation_error = maxval(abs(final_integrals - initial_integrals) / &
    max(1.0_dp, abs(initial_integrals)))
  call require(simulation_error < 3.0e-10_dp, &
    "dynamic multilevel hydro conserves the composite state")
  call write_multilevel_reactive_1d_csv( &
    "amr_multilevel_regrid.csv", species, simulated, ok)
  call require(ok, "ordered multilevel composite output")

  write(*, '(a,1x,es16.8)') "Regrid conservation:", regrid_error
  write(*, '(a,1x,es16.8)') "Simulation conservation:", simulation_error
  write(*, '(a,1x,es16.8)') "Overlap retention:", overlap_error
  write(*, '(a,i0)') "Regrid evaluations: ", &
    simulated%regrid_evaluations
  write(*, '(a)') "test_amr_multilevel_regrid_1d: PASS"

contains

  subroutine configure_case(local_config)
    type(reactive_1d_config), intent(out) :: local_config

    local_config = reactive_1d_config()
    local_config%nx = 24
    local_config%x_lower = 0.0_dp
    local_config%x_upper = 0.012_dp
    local_config%final_time = 2.0e-7_dp
    local_config%cfl = 0.20_dp
    local_config%maximum_steps = 100
    local_config%problem = "reactive_hotspot"
    local_config%riemann_solver = "rusanov"
    local_config%limiter = "mc"
    local_config%boundary_condition = "periodic"
    local_config%chemistry_enabled = .false.
    local_config%transport_enabled = .false.
    local_config%initial_temperature = 1200.0_dp
    local_config%initial_pressure = 101325.0_dp
    local_config%initial_velocity = 50.0_dp
    local_config%hotspot_temperature_rise = 250.0_dp
    local_config%hotspot_center = 0.006_dp
    local_config%hotspot_width = 0.0012_dp
    local_config%amr_enabled = .true.
    local_config%amr_reconstruction = "plm"
    local_config%amr_refinement_ratio = 2
    local_config%amr_max_levels = 3
    local_config%amr_regrid_interval = 1
    local_config%amr_tag_component = 1
    local_config%amr_buffer_cells = 1
    local_config%amr_minimum_patch_cells = 6
    local_config%amr_relative_gradient_threshold = 0.005_dp
    local_config%amr_absolute_gradient_threshold = 0.0_dp
    local_config%amr_scale_floor = 1.0e-12_dp
  end subroutine configure_case

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
      write(*, '(a)') "FAILED: " // trim(message)
      error stop 1
    end if
  end subroutine require

end program test_amr_multilevel_regrid_1d
