program test_amr_multilevel_reactive_1d
  use precision_mod, only: dp
  use state_indices_mod, only: irho
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
  use amr_hierarchy_1d_mod, only: restrict_average_1d
  use amr_multilevel_reactive_1d_mod, only: &
    amr_multilevel_reactive_solution_1d, &
    initialize_multilevel_reactive_1d, &
    multilevel_reactive_timestep_1d, advance_multilevel_reactive_1d, &
    multilevel_reactive_integrals_1d
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(reactive_1d_config) :: config, chemistry_config, ppm_config
  type(amr_multilevel_reactive_solution_1d) :: initial_solution
  type(amr_multilevel_reactive_solution_1d) :: solution, chemistry_solution
  type(amr_multilevel_reactive_solution_1d) :: ppm_solution
  type(amr_multilevel_reactive_solution_1d) :: weno_js_solution
  type(amr_multilevel_reactive_solution_1d) :: weno_z_solution
  real(dp), allocatable :: initial_integral(:), final_integral(:)
  real(dp), allocatable :: q(:)
  real(dp) :: dt, conservation_error, ppm_conservation_error
  real(dp) :: weno_js_conservation_error, weno_z_conservation_error
  real(dp) :: minimum_temperature
  real(dp) :: maximum_closure_error
  integer, parameter :: patch_lower(2) = [4, 5]
  integer, parameter :: patch_upper(2) = [21, 32]
  integer, parameter :: ratios(2) = [2, 2]
  logical :: ok
  integer :: step

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamics load")
  call load_h2o2_elementary_mechanism(reactions, ok)
  call require(ok, "mechanism load")
  call load_h2o2_elementary_transport(transport, ok)
  call require(ok, "transport database load")
  call configure_case(config)
  call initialize_multilevel_reactive_1d( &
    species, config, patch_lower, patch_upper, ratios, initial_solution, ok)
  call require(ok, "three-level reactive initialization")
  call require(initial_solution%level_count() == 3, &
    "runtime hierarchy owns three reactive levels")
  call require(initial_solution%hierarchy%level_cell_count(0) == 24, &
    "root reactive cell count")
  call require(initial_solution%hierarchy%level_cell_count(1) == 36, &
    "level-one reactive cell count")
  call require(initial_solution%hierarchy%level_cell_count(2) == 56, &
    "level-two reactive cell count")

  allocate(initial_integral(reactive_nvar(size(species))))
  allocate(final_integral(reactive_nvar(size(species))))
  call multilevel_reactive_integrals_1d( &
    initial_solution, initial_integral, ok)
  call require(ok, "initial multilevel composite integral")
  call multilevel_reactive_timestep_1d( &
    species, config, initial_solution, dt, ok, transport)
  call require(ok .and. dt > 0.0_dp, "all-level stable time step")
  dt = min(dt, 1.0e-8_dp)

  solution = initial_solution
  do step = 1, 2
    call advance_multilevel_reactive_1d( &
      species, reactions, config, dt, solution, ok, transport)
    call require(ok, "recursive reactive transport step")
  end do
  call require(solution%steps == 2, "coarse step count")
  call require(abs(solution%time - 2.0_dp * dt) <= 16.0_dp * &
    epsilon(1.0_dp) * dt, "coarse time accumulation")
  call multilevel_reactive_integrals_1d(solution, final_integral, ok)
  call require(ok, "final multilevel composite integral")
  conservation_error = maxval(abs(final_integral - initial_integral) / &
    max(1.0_dp, abs(initial_integral)))
  call require(conservation_error < 3.0e-10_dp, &
    "recursive hydro and transport conservation")
  call check_synchronization(solution)
  call inspect_solution(solution, minimum_temperature, maximum_closure_error)
  call require(minimum_temperature > 0.0_dp, &
    "multilevel thermodynamic positivity")
  call require(maximum_closure_error < 3.0e-10_dp, &
    "multilevel species closure")

  ppm_config = config
  ppm_config%amr_reconstruction = "characteristic_ppm"
  ppm_config%transport_enabled = .false.
  ppm_solution = initial_solution
  do step = 1, 2
    call advance_multilevel_reactive_1d( &
      species, reactions, ppm_config, dt, ppm_solution, ok)
    call require(ok, "recursive characteristic PPM hydro step")
  end do
  call multilevel_reactive_integrals_1d(ppm_solution, final_integral, ok)
  call require(ok, "characteristic PPM composite integral")
  ppm_conservation_error = maxval(abs(final_integral - initial_integral) / &
    max(1.0_dp, abs(initial_integral)))
  call require(ppm_conservation_error < 3.0e-10_dp, &
    "characteristic PPM reflux conservation")
  call check_synchronization(ppm_solution)
  call inspect_solution( &
    ppm_solution, minimum_temperature, maximum_closure_error)
  call require(minimum_temperature > 0.0_dp, &
    "characteristic PPM thermodynamic positivity")
  call require(maximum_closure_error < 3.0e-10_dp, &
    "characteristic PPM species closure")

  call run_weno_case(0, weno_js_solution, weno_js_conservation_error)
  call run_weno_case(1, weno_z_solution, weno_z_conservation_error)

  chemistry_config = config
  chemistry_config%chemistry_enabled = .true.
  chemistry_config%transport_enabled = .false.
  chemistry_solution = initial_solution
  call advance_multilevel_reactive_1d( &
    species, reactions, chemistry_config, 1.0e-10_dp, &
    chemistry_solution, ok)
  call require(ok, "recursive hierarchy chemistry split step")
  call check_synchronization(chemistry_solution)
  call inspect_solution( &
    chemistry_solution, minimum_temperature, maximum_closure_error)
  call require(minimum_temperature > 0.0_dp, &
    "chemistry hierarchy thermodynamic positivity")
  call require(maximum_closure_error < 3.0e-10_dp, &
    "chemistry hierarchy species closure")

  write(*, '(a,1x,es16.8)') "Three-level stable dt:", dt
  write(*, '(a,1x,es16.8)') "Composite conservation:", &
    conservation_error
  write(*, '(a,1x,es16.8)') "PPM composite conservation:", &
    ppm_conservation_error
  write(*, '(a,2(1x,es16.8))') "WENO5 JS/Z composite conservation:", &
    weno_js_conservation_error, weno_z_conservation_error
  write(*, '(a,1x,es16.8)') "Minimum temperature:", &
    minimum_temperature
  write(*, '(a)') "test_amr_multilevel_reactive_1d: PASS"

contains

  subroutine run_weno_case(scheme, local_solution, local_error)
    integer, intent(in) :: scheme
    type(amr_multilevel_reactive_solution_1d), intent(out) :: local_solution
    real(dp), intent(out) :: local_error

    type(reactive_1d_config) :: local_config
    integer :: local_step

    local_config = config
    local_config%amr_reconstruction = "characteristic_ppm"
    local_config%amr_hybrid_weno = .true.
    local_config%amr_weno_scheme = scheme
    local_config%transport_enabled = .false.
    local_solution = initial_solution
    do local_step = 1, 2
      call advance_multilevel_reactive_1d( &
        species, reactions, local_config, dt, local_solution, ok)
      call require(ok, "recursive hybrid WENO hydro step")
    end do
    call multilevel_reactive_integrals_1d( &
      local_solution, final_integral, ok)
    call require(ok, "hybrid WENO composite integral")
    local_error = maxval(abs(final_integral - initial_integral) / &
      max(1.0_dp, abs(initial_integral)))
    call require(local_error < 3.0e-10_dp, &
      "hybrid WENO reflux conservation")
    call check_synchronization(local_solution)
    call inspect_solution( &
      local_solution, minimum_temperature, maximum_closure_error)
    call require(minimum_temperature > 0.0_dp, &
      "hybrid WENO thermodynamic positivity")
    call require(maximum_closure_error < 3.0e-10_dp, &
      "hybrid WENO species closure")
  end subroutine run_weno_case

  subroutine configure_case(local_config)
    type(reactive_1d_config), intent(out) :: local_config

    local_config = reactive_1d_config()
    local_config%nx = 24
    local_config%x_lower = 0.0_dp
    local_config%x_upper = 0.012_dp
    local_config%final_time = 2.0e-8_dp
    local_config%cfl = 0.20_dp
    local_config%maximum_steps = 100
    local_config%problem = "reactive_hotspot"
    local_config%reconstruction = "pcm"
    local_config%riemann_solver = "rusanov"
    local_config%limiter = "mc"
    local_config%boundary_condition = "periodic"
    local_config%chemistry_enabled = .false.
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
    local_config%amr_reconstruction = "plm"
  end subroutine configure_case

  subroutine check_synchronization(local_solution)
    type(amr_multilevel_reactive_solution_1d), intent(in) :: local_solution

    real(dp), allocatable :: restricted(:, :)
    real(dp) :: error
    integer :: relation, child_cells, covered_cells, nvar

    nvar = size(local_solution%levels(1)%state, 1)
    do relation = 1, size(local_solution%hierarchy%interfaces)
      child_cells = local_solution%hierarchy%level_cell_count(relation)
      covered_cells = &
        local_solution%hierarchy%interfaces(relation)%covered_coarse_cells()
      allocate(restricted(nvar, covered_cells))
      call restrict_average_1d( &
        local_solution%levels(relation + 1)%state(:, 1:child_cells), &
        local_solution%hierarchy%interfaces(relation), restricted, ok)
      call require(ok, "relation restriction")
      error = maxval(abs(restricted - &
        local_solution%levels(relation)%state(:, &
          local_solution%hierarchy%interfaces(relation)%fine_coarse_lower: &
          local_solution%hierarchy%interfaces(relation)%fine_coarse_upper))) / &
        max(1.0_dp, maxval(abs(restricted)))
      call require(error < 5.0e-13_dp, &
        "every parent-child relation synchronized")
      deallocate(restricted)
    end do
  end subroutine check_synchronization

  subroutine inspect_solution( &
      local_solution, minimum, maximum_closure)
    type(amr_multilevel_reactive_solution_1d), intent(in) :: local_solution
    real(dp), intent(out) :: minimum, maximum_closure

    real(dp) :: temperature, sound_speed, closure
    logical :: local_ok
    integer :: level, cell, nx, k

    allocate(q(reactive_nprim(size(species))))
    minimum = huge(1.0_dp)
    maximum_closure = 0.0_dp
    do level = 1, local_solution%level_count()
      nx = local_solution%hierarchy%level_cell_count(level - 1)
      do cell = 1, nx
        call reactive_conserved_to_primitive( &
          species, local_solution%levels(level)%state(:, cell), &
          local_solution%levels(level)%temperature(cell), q, temperature, &
          sound_speed, local_ok)
        call require(local_ok, "valid multilevel reactive cell")
        closure = 0.0_dp
        do k = 1, size(species)
          closure = closure + q(reactive_mass_fraction_component(k))
          call require( &
            local_solution%levels(level)%state( &
              reactive_species_component(k), cell) >= -1.0e-14_dp * &
              max(1.0_dp, &
                local_solution%levels(level)%state(irho, cell)), &
            "multilevel species positivity")
        end do
        minimum = min(minimum, temperature)
        maximum_closure = max(maximum_closure, abs(closure - 1.0_dp))
      end do
    end do
    deallocate(q)
  end subroutine inspect_solution

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
      write(*, '(a)') "FAILED: " // trim(message)
      error stop 1
    end if
  end subroutine require

end program test_amr_multilevel_reactive_1d
