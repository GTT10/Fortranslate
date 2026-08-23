program test_amr_reactive_transport_1d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_species_component, reactive_conserved_to_primitive
  use amr_hierarchy_1d_mod, only: &
    composite_integral_1d, restrict_average_1d
  use amr_reactive_1d_mod, only: &
    amr_reactive_solution_1d, initialize_amr_reactive_1d, &
    advance_amr_reactive_1d, simulate_amr_reactive_1d
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(reactive_1d_config) :: config, reference_config
  type(amr_reactive_solution_1d) :: initial_solution
  type(amr_reactive_solution_1d) :: transport_step, reference_step
  type(amr_reactive_solution_1d) :: solution
  real(dp), allocatable :: initial_all(:), final_all(:)
  real(dp), allocatable :: restricted(:, :), q(:)
  real(dp) :: initial_integrals(5), final_integrals(5)
  real(dp) :: transport_span, reference_span
  real(dp) :: conservation_error, synchronization_error
  real(dp) :: minimum_temperature, species_error, closure
  real(dp), parameter :: comparison_dt = 1.0e-7_dp
  logical :: ok
  integer :: component, cell, fine_cells, covered_cells, k

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamics load")
  call load_h2o2_elementary_mechanism(reactions, ok)
  call require(ok, "mechanism load")
  call load_h2o2_elementary_transport(transport, ok)
  call require(ok, "transport database load")

  call configure_case(config)
  call initialize_amr_reactive_1d(species, config, initial_solution, ok)
  call require(ok, "AMR transport initialization")
  call require(initial_solution%fine_active(), "initial fine level active")
  call composite_all(initial_solution, initial_all, ok)
  call require(ok, "initial composite integral")

  transport_step = initial_solution
  reference_step = initial_solution
  config%viscosity_enabled = .false.
  config%species_diffusion_enabled = .false.
  config%barodiffusion_enabled = .false.
  call advance_amr_reactive_1d( &
    species, reactions, config, comparison_dt, transport_step, ok, transport)
  call require(ok, "single AMR conduction step")
  reference_config = config
  reference_config%transport_enabled = .false.
  call advance_amr_reactive_1d( &
    species, reactions, reference_config, comparison_dt, reference_step, ok)
  call require(ok, "single AMR reference step")
  transport_span = composite_temperature_span(transport_step)
  reference_span = composite_temperature_span(reference_step)
  call require(transport_span < reference_span, &
    "AMR conduction reduces the temperature span relative to inviscid AMR")

  call configure_case(config)
  call simulate_amr_reactive_1d( &
    species, reactions, config, solution, initial_integrals, &
    final_integrals, ok, transport)
  call require(ok, "dynamic AMR transport simulation")
  call require(solution%steps >= 2, "multiple coarse transport steps")
  call require(solution%regrid_evaluations >= 2, &
    "dynamic transport regrid evaluations")
  call require(solution%fine_active(), "final fine level active")

  call composite_all(solution, final_all, ok)
  call require(ok, "final composite integral")
  conservation_error = maxval(abs(final_all - initial_all) / &
    max(1.0_dp, abs(initial_all)))
  call require(conservation_error < 2.0e-10_dp, &
    "all AMR transport conserved components")
  species_error = 0.0_dp
  do k = 1, size(species)
    component = reactive_species_component(k)
    species_error = max(species_error, &
      abs(final_all(component) - initial_all(component)) / &
      max(1.0e-20_dp, abs(initial_all(component))))
  end do
  call require(species_error < 2.0e-10_dp, &
    "each periodic species mass is conserved")

  fine_cells = solution%hierarchy%fine%cell_count()
  covered_cells = solution%hierarchy%covered_coarse_cells()
  allocate(restricted(size(solution%coarse, 1), covered_cells))
  call restrict_average_1d( &
    solution%fine(:, 1:fine_cells), solution%hierarchy, restricted, ok)
  call require(ok, "final fine restriction")
  synchronization_error = maxval(abs( &
    restricted - solution%coarse(:, &
      solution%hierarchy%fine_coarse_lower: &
      solution%hierarchy%fine_coarse_upper))) / &
    max(1.0_dp, maxval(abs(restricted)))
  call require(synchronization_error < 5.0e-13_dp, &
    "covered coarse transport state synchronized")

  allocate(q(reactive_nprim(size(species))))
  minimum_temperature = huge(1.0_dp)
  do cell = 1, solution%hierarchy%fine_coarse_lower - 1
    call inspect_cell( &
      solution%coarse(:, cell), solution%coarse_temperature(cell))
  end do
  do cell = 1, fine_cells
    call inspect_cell(solution%fine(:, cell), solution%fine_temperature(cell))
  end do
  do cell = solution%hierarchy%fine_coarse_upper + 1, config%nx
    call inspect_cell( &
      solution%coarse(:, cell), solution%coarse_temperature(cell))
  end do
  call require(minimum_temperature > 0.0_dp, &
    "AMR transport thermodynamic positivity")

  write(*, '(a,1x,es16.8)') &
    "Transport/reference temperature spans:", transport_span
  write(*, '(a,1x,es16.8)') "Reference span:", reference_span
  write(*, '(a,1x,es16.8)') "All-component conservation:", &
    conservation_error
  write(*, '(a,1x,es16.8)') "Per-species conservation:", species_error
  write(*, '(a,1x,es16.8)') "Coarse/fine synchronization:", &
    synchronization_error
  write(*, '(a)') "test_amr_reactive_transport_1d: PASS"

contains

  subroutine configure_case(local_config)
    type(reactive_1d_config), intent(out) :: local_config

    local_config = reactive_1d_config()
    local_config%nx = 24
    local_config%x_lower = 0.0_dp
    local_config%x_upper = 0.012_dp
    local_config%final_time = 1.0e-6_dp
    local_config%cfl = 0.25_dp
    local_config%maximum_steps = 1000
    local_config%problem = "reactive_hotspot"
    local_config%reconstruction = "pcm"
    local_config%riemann_solver = "rusanov"
    local_config%boundary_condition = "periodic"
    local_config%chemistry_enabled = .false.
    local_config%transport_enabled = .true.
    local_config%viscosity_enabled = .true.
    local_config%thermal_conduction_enabled = .true.
    local_config%species_diffusion_enabled = .true.
    local_config%barodiffusion_enabled = .true.
    local_config%transport_cfl = 0.35_dp
    local_config%initial_temperature = 1200.0_dp
    local_config%initial_pressure = 101325.0_dp
    local_config%initial_velocity = 0.0_dp
    local_config%hotspot_temperature_rise = 250.0_dp
    local_config%hotspot_center = 0.006_dp
    local_config%hotspot_width = 0.0012_dp
    local_config%amr_enabled = .true.
    local_config%amr_reconstruction = "plm"
    local_config%amr_refinement_ratio = 2
    local_config%amr_regrid_interval = 1
    local_config%amr_tag_component = 1
    local_config%amr_buffer_cells = 1
    local_config%amr_minimum_patch_cells = 6
    local_config%amr_relative_gradient_threshold = 0.005_dp
    local_config%amr_absolute_gradient_threshold = 0.0_dp
    local_config%amr_scale_floor = 1.0e-12_dp
  end subroutine configure_case

  subroutine composite_all(local_solution, integral, local_ok)
    type(amr_reactive_solution_1d), intent(in) :: local_solution
    real(dp), allocatable, intent(out) :: integral(:)
    logical, intent(out) :: local_ok
    integer :: nx, nvar

    nx = size(local_solution%coarse, 2) - 2
    nvar = reactive_nvar(size(species))
    allocate(integral(nvar))
    if (local_solution%fine_active()) then
      call composite_integral_1d( &
        local_solution%coarse(:, 1:nx), &
        local_solution%fine(:, &
          1:local_solution%hierarchy%fine%cell_count()), &
        local_solution%hierarchy, integral, local_ok)
    else
      integral = local_solution%coarse_dx * &
        sum(local_solution%coarse(:, 1:nx), dim=2)
      local_ok = .true.
    end if
  end subroutine composite_all

  real(dp) function composite_temperature_span(local_solution) result(span)
    type(amr_reactive_solution_1d), intent(in) :: local_solution
    real(dp) :: minimum, maximum
    integer :: local_cell, local_fine_cells, nx

    nx = size(local_solution%coarse, 2) - 2
    minimum = huge(1.0_dp)
    maximum = -huge(1.0_dp)
    if (.not. local_solution%fine_active()) then
      minimum = minval(local_solution%coarse_temperature(1:nx))
      maximum = maxval(local_solution%coarse_temperature(1:nx))
    else
      do local_cell = 1, local_solution%hierarchy%fine_coarse_lower - 1
        minimum = min(minimum, &
          local_solution%coarse_temperature(local_cell))
        maximum = max(maximum, &
          local_solution%coarse_temperature(local_cell))
      end do
      local_fine_cells = local_solution%hierarchy%fine%cell_count()
      minimum = min(minimum, &
        minval(local_solution%fine_temperature(1:local_fine_cells)))
      maximum = max(maximum, &
        maxval(local_solution%fine_temperature(1:local_fine_cells)))
      do local_cell = local_solution%hierarchy%fine_coarse_upper + 1, nx
        minimum = min(minimum, &
          local_solution%coarse_temperature(local_cell))
        maximum = max(maximum, &
          local_solution%coarse_temperature(local_cell))
      end do
    end if
    span = maximum - minimum
  end function composite_temperature_span

  subroutine inspect_cell(state, temperature_guess)
    real(dp), intent(in) :: state(:), temperature_guess
    real(dp) :: temperature, sound_speed
    logical :: local_ok

    call reactive_conserved_to_primitive( &
      species, state, temperature_guess, q, temperature, sound_speed, local_ok)
    call require(local_ok, "valid final composite cell")
    closure = 0.0_dp
    do k = 1, size(species)
      closure = closure + q(reactive_mass_fraction_component(k))
      call require(state(reactive_species_component(k)) >= -1.0e-14_dp, &
        "AMR transport species positivity")
    end do
    call require(abs(closure - 1.0_dp) < 2.0e-10_dp, &
      "AMR transport species closure")
    minimum_temperature = min(minimum_temperature, temperature)
  end subroutine inspect_cell

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
      write(*, '(a)') "FAILED: " // trim(message)
      error stop 1
    end if
  end subroutine require

end program test_amr_reactive_transport_1d
