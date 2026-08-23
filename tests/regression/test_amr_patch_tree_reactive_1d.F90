program test_amr_patch_tree_reactive_1d
  use precision_mod, only: dp
  use state_indices_mod, only: irho, qp
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
  use amr_hierarchy_1d_mod, only: &
    amr_two_level_hierarchy_1d, restrict_average_1d
  use amr_patch_tree_1d_mod, only: amr_patch_level_plan_1d
  use amr_patch_tree_reactive_1d_mod, only: &
    amr_patch_tree_reactive_solution_1d, &
    initialize_patch_tree_reactive_1d, &
    patch_tree_reactive_timestep_1d, &
    advance_patch_tree_reactive_1d, &
    advance_patch_tree_reactive_hydro_1d, &
    patch_tree_reactive_integrals_1d
  implicit none

  real(dp), parameter :: conservation_tolerance = 3.0e-10_dp
  real(dp), parameter :: transport_conservation_tolerance = 2.0e-9_dp
  real(dp), parameter :: synchronization_tolerance = 5.0e-13_dp
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(reactive_1d_config) :: config, chemistry_config
  type(reactive_1d_config) :: transport_config, split_control_config
  type(amr_patch_level_plan_1d), allocatable :: plans(:)
  type(amr_patch_tree_reactive_solution_1d) :: solution
  type(amr_patch_tree_reactive_solution_1d) :: hydro_control
  type(amr_patch_tree_reactive_solution_1d) :: chemistry_solution
  type(amr_patch_tree_reactive_solution_1d) :: transport_solution
  type(amr_patch_tree_reactive_solution_1d) :: transport_control
  type(amr_patch_tree_reactive_solution_1d) :: rejected_solution
  real(dp), allocatable :: initial_integral(:), final_integral(:)
  real(dp), allocatable :: chemistry_initial_integral(:)
  real(dp), allocatable :: chemistry_final_integral(:)
  real(dp), allocatable :: transport_initial_integral(:)
  real(dp), allocatable :: transport_final_integral(:)
  real(dp), allocatable :: q(:)
  real(dp) :: dt, conservation_error, synchronization_error
  real(dp) :: chemistry_dt, chemistry_difference
  real(dp) :: chemistry_conservation_error
  real(dp) :: transport_dt, transport_difference
  real(dp) :: transport_conservation_error
  real(dp) :: minimum_temperature, minimum_pressure
  real(dp) :: maximum_closure_error
  logical :: ok

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "patch-tree thermodynamics load")
  call load_h2o2_elementary_mechanism(reactions, ok)
  call require(ok, "patch-tree chemistry mechanism load")
  call load_h2o2_elementary_transport(transport, ok)
  call require(ok, "patch-tree transport database load")
  call configure_case(config)
  call configure_plans(plans)
  call initialize_patch_tree_reactive_1d( &
    species, config, plans, solution, ok)
  call require(ok .and. solution%is_valid(), &
    "four-level reactive patch-tree initialization")
  call require(solution%level_count() == 4, &
    "reactive patch-tree level count")
  call require(all([ &
    size(solution%levels(1)%patches), size(solution%levels(2)%patches), &
    size(solution%levels(3)%patches), size(solution%levels(4)%patches)] == &
    [1, 2, 3, 2]), "reactive patch-tree branching")

  allocate(initial_integral(reactive_nvar(size(species))))
  allocate(final_integral(reactive_nvar(size(species))))
  call patch_tree_reactive_integrals_1d(solution, initial_integral, ok)
  call require(ok, "initial patch-tree composite integral")
  call patch_tree_reactive_timestep_1d(species, config, solution, dt, ok)
  call require(ok .and. dt > 0.0_dp, "all-patch stable time step")
  dt = min(0.10_dp * dt, 2.0e-8_dp)

  call advance_patch_tree_reactive_hydro_1d( &
    species, config, dt, solution, ok)
  call require(ok .and. solution%is_valid(), &
    "recursive patch-tree hydro step")
  call require(solution%steps == 1 .and. &
    abs(solution%time - dt) <= 16.0_dp * epsilon(1.0_dp) * dt, &
    "patch-tree time and step accounting")
  call require(all(solution%level_advances == [1, 4, 12, 16]), &
    "recursive subcycling counts")

  call patch_tree_reactive_integrals_1d(solution, final_integral, ok)
  call require(ok, "final patch-tree composite integral")
  conservation_error = maxval(abs(final_integral - initial_integral) / &
    max(1.0_dp, abs(initial_integral)))
  call require(conservation_error < conservation_tolerance, &
    "branched patch-tree reflux conservation")
  call check_synchronization(solution, synchronization_error)
  call require(synchronization_error < synchronization_tolerance, &
    "every patch-tree parent-child relation synchronized")
  call inspect_solution( &
    solution, minimum_temperature, minimum_pressure, maximum_closure_error)
  call require(minimum_temperature > 0.0_dp, &
    "patch-tree thermodynamic temperature positivity")
  call require(minimum_pressure > 0.0_dp, &
    "patch-tree thermodynamic pressure positivity")
  call require(maximum_closure_error < 3.0e-10_dp, &
    "patch-tree species closure")

  allocate(chemistry_initial_integral(reactive_nvar(size(species))))
  allocate(chemistry_final_integral(reactive_nvar(size(species))))
  call patch_tree_reactive_integrals_1d( &
    solution, chemistry_initial_integral, ok)
  call require(ok, "pre-chemistry patch-tree composite integral")
  hydro_control = solution
  chemistry_solution = solution
  chemistry_config = config
  chemistry_config%chemistry_enabled = .true.
  chemistry_config%chemistry_relative_tolerance = 1.0e-8_dp
  chemistry_config%chemistry_absolute_tolerance = 1.0e-14_dp
  chemistry_dt = min(dt, 1.0e-10_dp)
  call advance_patch_tree_reactive_1d( &
    species, reactions, config, chemistry_dt, hydro_control, ok)
  call require(ok, "patch-tree hydro control Strang step")
  call advance_patch_tree_reactive_1d( &
    species, reactions, chemistry_config, chemistry_dt, &
    chemistry_solution, ok)
  call require(ok .and. chemistry_solution%is_valid(), &
    "patch-tree chemistry Strang step")
  chemistry_difference = reactive_species_difference( &
    chemistry_solution, hydro_control)
  call require(chemistry_difference > 100.0_dp * epsilon(1.0_dp), &
    "chemistry changes the recursive patch-tree species state")
  call require(chemistry_solution%steps == 2 .and. &
    all(chemistry_solution%level_advances == [2, 8, 24, 32]) .and. &
    abs(chemistry_solution%time - dt - chemistry_dt) <= &
      32.0_dp * epsilon(1.0_dp) * max(dt, chemistry_dt), &
    "chemistry composition preserves hydro advance accounting")
  call patch_tree_reactive_integrals_1d( &
    chemistry_solution, chemistry_final_integral, ok)
  call require(ok, "post-chemistry patch-tree composite integral")
  chemistry_conservation_error = maxval(abs( &
    chemistry_final_integral(1:5) - chemistry_initial_integral(1:5)) / &
    max(1.0_dp, abs(chemistry_initial_integral(1:5))))
  call require(chemistry_conservation_error < conservation_tolerance, &
    "patch-tree chemistry mass momentum energy conservation")
  call check_synchronization(chemistry_solution, synchronization_error)
  call require(synchronization_error < synchronization_tolerance, &
    "post-chemistry patch-tree synchronization")
  call inspect_solution( &
    chemistry_solution, minimum_temperature, minimum_pressure, &
    maximum_closure_error)
  call require(minimum_temperature > 0.0_dp .and. &
    minimum_pressure > 0.0_dp .and. maximum_closure_error < 3.0e-10_dp, &
    "post-chemistry patch-tree physical state")

  call configure_transport_case(config, transport_config)
  call initialize_patch_tree_reactive_1d( &
    species, transport_config, plans, transport_solution, ok)
  call require(ok .and. transport_solution%is_valid(), &
    "patch-tree reacting transport initialization")
  call patch_tree_reactive_timestep_1d( &
    species, transport_config, transport_solution, transport_dt, ok)
  call require(.not. ok, "patch-tree transport timestep requires database")
  call patch_tree_reactive_timestep_1d( &
    species, transport_config, transport_solution, transport_dt, ok, &
    transport)
  call require(ok .and. transport_dt > 0.0_dp, &
    "all-patch hydro-transport stable timestep")
  transport_dt = min(transport_dt, 1.0e-10_dp)
  allocate(transport_initial_integral(reactive_nvar(size(species))))
  allocate(transport_final_integral(reactive_nvar(size(species))))
  call patch_tree_reactive_integrals_1d( &
    transport_solution, transport_initial_integral, ok)
  call require(ok, "pre-transport patch-tree composite integral")
  transport_control = transport_solution
  split_control_config = transport_config
  split_control_config%transport_enabled = .false.
  call advance_patch_tree_reactive_1d( &
    species, reactions, split_control_config, transport_dt, &
    transport_control, ok)
  call require(ok, "patch-tree chemistry-hydro control step")
  call advance_patch_tree_reactive_1d( &
    species, reactions, transport_config, transport_dt, &
    transport_solution, ok, transport)
  call require(ok .and. transport_solution%is_valid(), &
    "patch-tree reacting transport split step")
  transport_difference = patch_tree_state_difference( &
    transport_solution, transport_control)
  call require(transport_difference > 100.0_dp * epsilon(1.0_dp), &
    "transport changes the recursive patch-tree state")
  call require(all(transport_solution%level_advances == [1, 4, 12, 16]) &
    .and. all(transport_solution%transport_level_advances == &
      [2, 16, 96, 256]), "recursive parabolic subcycling counts")
  call patch_tree_reactive_integrals_1d( &
    transport_solution, transport_final_integral, ok)
  call require(ok, "post-transport patch-tree composite integral")
  transport_conservation_error = maxval(abs( &
    transport_final_integral(1:5) - transport_initial_integral(1:5)) / &
    max(1.0_dp, abs(transport_initial_integral(1:5))))
  call require( &
    transport_conservation_error < transport_conservation_tolerance, &
    "patch-tree reacting transport conservation")
  call check_synchronization(transport_solution, synchronization_error)
  call require(synchronization_error < 8.0e-13_dp, &
    "post-transport patch-tree synchronization")
  call inspect_solution( &
    transport_solution, minimum_temperature, minimum_pressure, &
    maximum_closure_error)
  call require(minimum_temperature > 0.0_dp .and. &
    minimum_pressure > 0.0_dp .and. maximum_closure_error < 3.0e-10_dp, &
    "post-transport patch-tree physical state")

  rejected_solution = transport_solution
  call advance_patch_tree_reactive_1d( &
    species, reactions, transport_config, transport_dt, &
    rejected_solution, ok)
  call require(.not. ok .and. &
    patch_tree_state_difference(rejected_solution, transport_solution) == &
      0.0_dp .and. rejected_solution%steps == transport_solution%steps .and. &
    all(rejected_solution%transport_level_advances == &
      transport_solution%transport_level_advances), &
    "missing transport database leaves patch tree unchanged")

  write(*, '(a,1x,es16.8)') &
    "Patch-tree conservation error:", conservation_error
  write(*, '(a,1x,es16.8)') &
    "Patch-tree synchronization error:", synchronization_error
  write(*, '(a,1x,es16.8)') &
    "Patch-tree chemistry difference:", chemistry_difference
  write(*, '(a,1x,es16.8)') &
    "Patch-tree chemistry conservation:", chemistry_conservation_error
  write(*, '(a,1x,es16.8)') &
    "Patch-tree transport difference:", transport_difference
  write(*, '(a,1x,es16.8)') &
    "Patch-tree transport conservation:", transport_conservation_error
  write(*, '(a)') "test_amr_patch_tree_reactive_1d: PASS"

contains

  subroutine configure_case(local_config)
    type(reactive_1d_config), intent(out) :: local_config

    local_config = reactive_1d_config()
    local_config%nx = 32
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
    local_config%amr_reconstruction = "pcm"
  end subroutine configure_case

  subroutine configure_transport_case(base_config, local_config)
    type(reactive_1d_config), intent(in) :: base_config
    type(reactive_1d_config), intent(out) :: local_config

    local_config = base_config
    local_config%problem = "reactive_hotspot"
    local_config%chemistry_enabled = .true.
    local_config%transport_enabled = .true.
    local_config%viscosity_enabled = .true.
    local_config%thermal_conduction_enabled = .true.
    local_config%species_diffusion_enabled = .true.
    local_config%barodiffusion_enabled = .true.
    local_config%transport_cfl = 0.30_dp
    local_config%initial_velocity = 0.0_dp
    local_config%hotspot_temperature_rise = 200.0_dp
    local_config%hotspot_center = 0.006_dp
    local_config%hotspot_width = 0.0012_dp
  end subroutine configure_transport_case

  subroutine configure_plans(local_plans)
    type(amr_patch_level_plan_1d), allocatable, intent(out) :: local_plans(:)

    allocate(local_plans(3))

    local_plans(1)%refinement_ratio = 2
    allocate(local_plans(1)%patches(2))
    local_plans(1)%patches(1)%parent_patch = 1
    local_plans(1)%patches(1)%lower = 4
    local_plans(1)%patches(1)%upper = 11
    local_plans(1)%patches(2)%parent_patch = 1
    local_plans(1)%patches(2)%lower = 20
    local_plans(1)%patches(2)%upper = 27

    local_plans(2)%refinement_ratio = 2
    allocate(local_plans(2)%patches(3))
    local_plans(2)%patches(1)%parent_patch = 1
    local_plans(2)%patches(1)%lower = 3
    local_plans(2)%patches(1)%upper = 8
    local_plans(2)%patches(2)%parent_patch = 1
    local_plans(2)%patches(2)%lower = 11
    local_plans(2)%patches(2)%upper = 14
    local_plans(2)%patches(3)%parent_patch = 2
    local_plans(2)%patches(3)%lower = 5
    local_plans(2)%patches(3)%upper = 12

    local_plans(3)%refinement_ratio = 2
    allocate(local_plans(3)%patches(2))
    local_plans(3)%patches(1)%parent_patch = 1
    local_plans(3)%patches(1)%lower = 3
    local_plans(3)%patches(1)%upper = 10
    local_plans(3)%patches(2)%parent_patch = 3
    local_plans(3)%patches(2)%lower = 4
    local_plans(3)%patches(2)%upper = 13
  end subroutine configure_plans

  subroutine check_synchronization(local_solution, maximum_error)
    type(amr_patch_tree_reactive_solution_1d), intent(in) :: local_solution
    real(dp), intent(out) :: maximum_error

    type(amr_two_level_hierarchy_1d) :: geometry
    real(dp), allocatable :: restricted(:, :)
    real(dp) :: error
    logical :: local_ok
    integer :: relation, parent, child, global_child
    integer :: fine_cells, covered_cells, nvar

    maximum_error = 0.0_dp
    nvar = size(local_solution%levels(1)%patches(1)%state, 1)
    do relation = 1, size(local_solution%hierarchy%relations)
      do parent = 1, local_solution%hierarchy%relations(relation)% &
          parent_patch_count()
        do child = 1, local_solution%hierarchy%relations(relation)% &
            child_sets(parent)%patch_count()
          geometry = local_solution%hierarchy%relations(relation)% &
            child_sets(parent)%patches(child)
          global_child = local_solution%hierarchy%relations(relation)% &
            child_index(parent, child)
          fine_cells = geometry%fine%cell_count()
          covered_cells = geometry%covered_coarse_cells()
          allocate(restricted(nvar, covered_cells))
          call restrict_average_1d( &
            local_solution%levels(relation + 1)%patches(global_child)% &
              state(:, 1:fine_cells), geometry, restricted, local_ok)
          call require(local_ok, "patch-tree relation restriction")
          error = maxval(abs(restricted - &
            local_solution%levels(relation)%patches(parent)%state(:, &
              geometry%fine_coarse_lower:geometry%fine_coarse_upper))) / &
            max(1.0_dp, maxval(abs(restricted)))
          maximum_error = max(maximum_error, error)
          deallocate(restricted)
        end do
      end do
    end do
  end subroutine check_synchronization

  subroutine inspect_solution( &
      local_solution, minimum_t, minimum_p, maximum_closure)
    type(amr_patch_tree_reactive_solution_1d), intent(in) :: local_solution
    real(dp), intent(out) :: minimum_t, minimum_p, maximum_closure

    real(dp) :: temperature, sound_speed, closure, density
    logical :: local_ok
    integer :: level, patch, cell, nx, species_index

    allocate(q(reactive_nprim(size(species))))
    minimum_t = huge(1.0_dp)
    minimum_p = huge(1.0_dp)
    maximum_closure = 0.0_dp
    do level = 1, local_solution%level_count()
      do patch = 1, size(local_solution%levels(level)%patches)
        nx = size(local_solution%levels(level)%patches(patch)%state, 2) - 2
        do cell = 1, nx
          call reactive_conserved_to_primitive( &
            species, local_solution%levels(level)%patches(patch)% &
              state(:, cell), &
            local_solution%levels(level)%patches(patch)%temperature(cell), &
            q, temperature, sound_speed, local_ok)
          call require(local_ok .and. sound_speed > 0.0_dp, &
            "valid patch-tree reactive cell")
          density = local_solution%levels(level)%patches(patch)% &
            state(irho, cell)
          closure = 0.0_dp
          do species_index = 1, size(species)
            call require(local_solution%levels(level)%patches(patch)% &
              state(reactive_species_component(species_index), cell) >= &
              -1.0e-14_dp * max(1.0_dp, density), &
              "patch-tree species positivity")
            closure = closure + local_solution%levels(level)%patches(patch)% &
              state(reactive_species_component(species_index), cell) / &
              density
          end do
          minimum_t = min(minimum_t, temperature)
          minimum_p = min(minimum_p, q(qp))
          maximum_closure = max(maximum_closure, abs(closure - 1.0_dp))
        end do
      end do
    end do
    deallocate(q)
  end subroutine inspect_solution

  real(dp) function reactive_species_difference( &
      first, second) result(difference)
    type(amr_patch_tree_reactive_solution_1d), intent(in) :: first, second

    integer :: level, patch, cell, nx, species_index, component

    difference = 0.0_dp
    do level = 1, first%level_count()
      do patch = 1, size(first%levels(level)%patches)
        nx = size(first%levels(level)%patches(patch)%state, 2) - 2
        do cell = 1, nx
          do species_index = 1, size(species)
            component = reactive_species_component(species_index)
            difference = max(difference, abs( &
              first%levels(level)%patches(patch)%state(component, cell) - &
              second%levels(level)%patches(patch)%state(component, cell)))
          end do
        end do
      end do
    end do
  end function reactive_species_difference

  real(dp) function patch_tree_state_difference( &
      first, second) result(difference)
    type(amr_patch_tree_reactive_solution_1d), intent(in) :: first, second

    integer :: level, patch, nx

    difference = 0.0_dp
    do level = 1, first%level_count()
      do patch = 1, size(first%levels(level)%patches)
        nx = size(first%levels(level)%patches(patch)%state, 2) - 2
        difference = max(difference, maxval(abs( &
          first%levels(level)%patches(patch)%state(:, 1:nx) - &
          second%levels(level)%patches(patch)%state(:, 1:nx))))
        difference = max(difference, maxval(abs( &
          first%levels(level)%patches(patch)%temperature(1:nx) - &
          second%levels(level)%patches(patch)%temperature(1:nx))))
      end do
    end do
  end function patch_tree_state_difference

  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      write(*, '(a,1x,a)') "FAIL:", trim(label)
      error stop 1
    end if
  end subroutine require

end program test_amr_patch_tree_reactive_1d
