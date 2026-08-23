module amr_reactive_1d_mod
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_conserved_to_primitive, reactive_riemann_flux_x, &
    reactive_cfl_timestep, initialize_reactive_1d, &
    advance_reactive_chemistry
  use amr_hierarchy_1d_mod, only: &
    amr_two_level_hierarchy_1d, amr_flux_register_1d, &
    initialize_flux_register_1d, accumulate_coarse_flux_1d, &
    accumulate_fine_flux_1d, reflux_1d, average_down_1d, &
    composite_integral_1d
  use amr_regrid_1d_mod, only: &
    amr_tagging_criteria_1d, amr_regrid_plan_1d, &
    plan_gradient_regrid_1d, regrid_two_level_state_1d
  implicit none
  private

  type, public :: amr_reactive_solution_1d
    type(amr_two_level_hierarchy_1d) :: hierarchy
    real(dp), allocatable :: coarse(:, :)
    real(dp), allocatable :: coarse_temperature(:)
    real(dp), allocatable :: fine(:, :)
    real(dp), allocatable :: fine_temperature(:)
    real(dp) :: x_lower = 0.0_dp
    real(dp) :: x_upper = 0.0_dp
    real(dp) :: coarse_dx = 0.0_dp
    real(dp) :: time = 0.0_dp
    integer :: steps = 0
    integer :: regrid_evaluations = 0
    integer :: regrids = 0
  contains
    procedure :: fine_active => amr_reactive_fine_active
  end type amr_reactive_solution_1d

  public :: initialize_amr_reactive_1d
  public :: amr_reactive_timestep_1d
  public :: advance_amr_reactive_1d
  public :: regrid_amr_reactive_1d
  public :: simulate_amr_reactive_1d
  public :: amr_reactive_integrals_1d
  public :: write_amr_reactive_1d_csv

contains

  pure logical function amr_reactive_fine_active(self) result(active)
    class(amr_reactive_solution_1d), intent(in) :: self

    active = self%hierarchy%is_valid() .and. allocated(self%fine) .and. &
      allocated(self%fine_temperature)
  end function amr_reactive_fine_active

  subroutine initialize_amr_reactive_1d( &
      species, config, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_reactive_solution_1d), intent(out) :: solution
    logical, intent(out) :: ok

    type(amr_tagging_criteria_1d) :: criteria
    type(amr_regrid_plan_1d) :: plan
    logical, allocatable :: tags(:)
    logical :: local_ok, changed

    ok = .false.
    if (.not. valid_amr_configuration(config, reactive_nvar(size(species)))) &
      return
    call initialize_reactive_1d( &
      species, config, solution%coarse, solution%coarse_temperature, &
      solution%coarse_dx, local_ok)
    if (.not. local_ok) return
    solution%x_lower = config%x_lower
    solution%x_upper = config%x_upper
    call criteria_from_config(config, criteria)
    allocate(tags(config%nx))
    call plan_gradient_regrid_1d( &
      solution%coarse(:, 1:config%nx), criteria, tags, plan, local_ok)
    if (.not. local_ok) return
    call apply_regrid_plan_1d( &
      species, config, plan, solution, changed, local_ok)
    if (.not. local_ok) return
    solution%time = 0.0_dp
    solution%steps = 0
    solution%regrid_evaluations = 1
    if (changed) solution%regrids = 1
    ok = .true.
  end subroutine initialize_amr_reactive_1d

  subroutine amr_reactive_timestep_1d( &
      species, config, solution, dt, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_reactive_solution_1d), intent(in) :: solution
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok

    real(dp) :: coarse_dt, fine_dt
    logical :: local_ok

    dt = 0.0_dp
    ok = .false.
    call reactive_cfl_timestep( &
      species, solution%coarse, solution%coarse_temperature, config%nx, &
      solution%coarse_dx, config%cfl, coarse_dt, local_ok)
    if (.not. local_ok) return
    dt = coarse_dt
    if (solution%fine_active()) then
      call reactive_cfl_timestep( &
        species, solution%fine, solution%fine_temperature, &
        solution%hierarchy%fine%cell_count(), solution%hierarchy%fine_dx, &
        config%cfl, fine_dt, local_ok)
      if (.not. local_ok) return
      dt = min(dt, &
        real(solution%hierarchy%refinement_ratio, dp) * fine_dt)
    end if
    ok = dt > 0.0_dp
  end subroutine amr_reactive_timestep_1d

  subroutine advance_amr_reactive_1d( &
      species, reactions, config, dt, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: dt
    type(amr_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok

    type(amr_reactive_solution_1d) :: backup
    type(amr_flux_register_1d) :: flux_register
    real(dp), allocatable :: coarse_start(:, :), coarse_end(:, :)
    real(dp), allocatable :: coarse_flux(:, :), fine_flux(:, :)
    real(dp) :: fine_dt, alpha
    logical :: local_ok
    integer :: substep, ratio, fine_cells

    ok = .false.
    if (dt <= 0.0_dp .or. .not. &
        valid_amr_configuration(config, reactive_nvar(size(species)))) return
    backup = solution

    if (config%chemistry_enabled) then
      call advance_reactive_chemistry( &
        species, reactions, solution%coarse, solution%coarse_temperature, &
        config%nx, 0.5_dp * dt, config%chemistry_relative_tolerance, &
        config%chemistry_absolute_tolerance, config%boundary_condition, &
        local_ok)
      if (.not. local_ok) then
        solution = backup
        return
      end if
      if (solution%fine_active()) then
        fine_cells = solution%hierarchy%fine%cell_count()
        call advance_reactive_chemistry( &
          species, reactions, solution%fine, solution%fine_temperature, &
          fine_cells, 0.5_dp * dt, config%chemistry_relative_tolerance, &
          config%chemistry_absolute_tolerance, "outflow", local_ok)
        if (.not. local_ok) then
          solution = backup
          return
        end if
      end if
    end if

    allocate(coarse_start( &
      size(solution%coarse, 1), 0:config%nx + 1))
    coarse_start = solution%coarse
    allocate(coarse_flux(reactive_nvar(size(species)), 0:config%nx))
    call advance_pcm_level_1d( &
      species, solution%coarse, solution%coarse_temperature, config%nx, &
      solution%coarse_dx, dt, config%riemann_solver, .true., &
      config%boundary_condition, coarse_flux, local_ok)
    if (.not. local_ok) then
      solution = backup
      return
    end if
    allocate(coarse_end( &
      size(solution%coarse, 1), 0:config%nx + 1))
    coarse_end = solution%coarse

    if (solution%fine_active()) then
      ratio = solution%hierarchy%refinement_ratio
      fine_cells = solution%hierarchy%fine%cell_count()
      fine_dt = dt / real(ratio, dp)
      allocate(fine_flux(reactive_nvar(size(species)), 0:fine_cells))
      call initialize_flux_register_1d( &
        flux_register, reactive_nvar(size(species)), local_ok)
      if (.not. local_ok) then
        solution = backup
        return
      end if
      call accumulate_coarse_flux_1d( &
        flux_register, &
        coarse_flux(:, solution%hierarchy%fine_coarse_lower - 1), &
        coarse_flux(:, solution%hierarchy%fine_coarse_upper), dt, local_ok)
      if (.not. local_ok) then
        solution = backup
        return
      end if

      do substep = 1, ratio
        alpha = real(substep - 1, dp) / real(ratio, dp)
        call fill_fine_ghosts_1d( &
          species, solution%hierarchy, coarse_start, coarse_end, alpha, &
          solution%fine, solution%fine_temperature, local_ok)
        if (.not. local_ok) then
          solution = backup
          return
        end if
        call advance_pcm_level_1d( &
          species, solution%fine, solution%fine_temperature, fine_cells, &
          solution%hierarchy%fine_dx, fine_dt, config%riemann_solver, &
          .false., "coarse_fine", fine_flux, local_ok)
        if (.not. local_ok) then
          solution = backup
          return
        end if
        call accumulate_fine_flux_1d( &
          flux_register, fine_flux(:, 0), fine_flux(:, fine_cells), &
          fine_dt, local_ok)
        if (.not. local_ok) then
          solution = backup
          return
        end if
      end do
      call reflux_1d( &
        solution%coarse(:, 1:config%nx), solution%hierarchy, &
        flux_register, local_ok)
      if (.not. local_ok) then
        solution = backup
        return
      end if
      call synchronize_levels_1d(species, config, solution, local_ok)
      if (.not. local_ok) then
        solution = backup
        return
      end if
    end if

    if (config%chemistry_enabled) then
      call advance_reactive_chemistry( &
        species, reactions, solution%coarse, solution%coarse_temperature, &
        config%nx, 0.5_dp * dt, config%chemistry_relative_tolerance, &
        config%chemistry_absolute_tolerance, config%boundary_condition, &
        local_ok)
      if (.not. local_ok) then
        solution = backup
        return
      end if
      if (solution%fine_active()) then
        fine_cells = solution%hierarchy%fine%cell_count()
        call fill_fine_ghosts_1d( &
          species, solution%hierarchy, solution%coarse, solution%coarse, &
          1.0_dp, solution%fine, solution%fine_temperature, local_ok)
        if (.not. local_ok) then
          solution = backup
          return
        end if
        call advance_reactive_chemistry( &
          species, reactions, solution%fine, solution%fine_temperature, &
          fine_cells, 0.5_dp * dt, config%chemistry_relative_tolerance, &
          config%chemistry_absolute_tolerance, "outflow", local_ok)
        if (.not. local_ok) then
          solution = backup
          return
        end if
        call synchronize_levels_1d(species, config, solution, local_ok)
        if (.not. local_ok) then
          solution = backup
          return
        end if
      end if
    end if

    call fill_physical_ghosts_1d( &
      solution%coarse, solution%coarse_temperature, config%nx, &
      config%boundary_condition, local_ok)
    if (.not. local_ok) then
      solution = backup
      return
    end if
    if (solution%fine_active()) then
      call fill_fine_ghosts_1d( &
        species, solution%hierarchy, solution%coarse, solution%coarse, &
        1.0_dp, solution%fine, solution%fine_temperature, local_ok)
      if (.not. local_ok) then
        solution = backup
        return
      end if
    end if
    solution%time = solution%time + dt
    solution%steps = solution%steps + 1
    ok = .true.
  end subroutine advance_amr_reactive_1d

  subroutine regrid_amr_reactive_1d(species, config, solution, changed, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: changed, ok

    type(amr_tagging_criteria_1d) :: criteria
    type(amr_regrid_plan_1d) :: plan
    logical, allocatable :: tags(:)
    logical :: local_ok

    changed = .false.
    ok = .false.
    call criteria_from_config(config, criteria)
    allocate(tags(config%nx))
    call plan_gradient_regrid_1d( &
      solution%coarse(:, 1:config%nx), criteria, tags, plan, local_ok)
    solution%regrid_evaluations = solution%regrid_evaluations + 1
    if (.not. local_ok) return
    call apply_regrid_plan_1d( &
      species, config, plan, solution, changed, local_ok)
    if (.not. local_ok) return
    if (changed) solution%regrids = solution%regrids + 1
    ok = .true.
  end subroutine regrid_amr_reactive_1d

  subroutine simulate_amr_reactive_1d( &
      species, reactions, config, solution, initial_integrals, &
      final_integrals, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_reactive_solution_1d), intent(out) :: solution
    real(dp), intent(out) :: initial_integrals(5), final_integrals(5)
    logical, intent(out) :: ok

    real(dp) :: dt, tolerance
    logical :: local_ok, changed

    initial_integrals = 0.0_dp
    final_integrals = 0.0_dp
    call initialize_amr_reactive_1d(species, config, solution, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call amr_reactive_integrals_1d(solution, initial_integrals, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    tolerance = 50.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, config%final_time)
    do while (solution%time < config%final_time - tolerance)
      if (solution%steps >= config%maximum_steps) then
        ok = .false.
        return
      end if
      call amr_reactive_timestep_1d(species, config, solution, dt, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      dt = min(dt, config%final_time - solution%time)
      call advance_amr_reactive_1d( &
        species, reactions, config, dt, solution, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      if (mod(solution%steps, config%amr_regrid_interval) == 0) then
        call regrid_amr_reactive_1d( &
          species, config, solution, changed, local_ok)
        if (.not. local_ok) then
          ok = .false.
          return
        end if
      end if
    end do
    solution%time = config%final_time
    call amr_reactive_integrals_1d(solution, final_integrals, ok)
  end subroutine simulate_amr_reactive_1d

  subroutine amr_reactive_integrals_1d(solution, integrals, ok)
    type(amr_reactive_solution_1d), intent(in) :: solution
    real(dp), intent(out) :: integrals(5)
    logical, intent(out) :: ok

    real(dp), allocatable :: all_integrals(:)
    integer :: nx, nvar

    integrals = 0.0_dp
    ok = allocated(solution%coarse) .and. solution%coarse_dx > 0.0_dp
    if (.not. ok) return
    nx = size(solution%coarse, 2) - 2
    nvar = size(solution%coarse, 1)
    allocate(all_integrals(nvar))
    if (solution%fine_active()) then
      call composite_integral_1d( &
        solution%coarse(:, 1:nx), &
        solution%fine(:, 1:solution%hierarchy%fine%cell_count()), &
        solution%hierarchy, all_integrals, ok)
      if (.not. ok) return
    else
      all_integrals = solution%coarse_dx * &
        sum(solution%coarse(:, 1:nx), dim=2)
    end if
    integrals = all_integrals([irho, imx, imy, imz, iet])
    ok = .true.
  end subroutine amr_reactive_integrals_1d

  subroutine write_amr_reactive_1d_csv(path, species, solution, ok)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(amr_reactive_solution_1d), intent(in) :: solution
    logical, intent(out) :: ok

    real(dp), allocatable :: q(:)
    real(dp) :: x
    logical :: local_ok
    integer :: unit, status, cell, nx, fine_cells, global_fine, k

    ok = .false.
    nx = size(solution%coarse, 2) - 2
    allocate(q(reactive_nprim(size(species))))
    open(newunit=unit, file=trim(path), status="replace", action="write", &
      iostat=status)
    if (status /= 0) return
    write(unit, '(a)', advance='no') &
      "level,cell_dx,time,x,rho,u,v,w,pressure,temperature,rhoE"
    do k = 1, size(species)
      write(unit, '(a)', advance='no') ",Y_" // trim(species(k)%name)
    end do
    write(unit, '(a)') ""

    if (.not. solution%fine_active()) then
      do cell = 1, nx
        x = solution%x_lower + &
          (real(cell, dp) - 0.5_dp) * solution%coarse_dx
        call write_amr_cell( &
          unit, 0, solution%coarse_dx, solution%time, x, species, &
          solution%coarse(:, cell), solution%coarse_temperature(cell), q, &
          local_ok)
        if (.not. local_ok) then
          close(unit)
          return
        end if
      end do
    else
      do cell = 1, solution%hierarchy%fine_coarse_lower - 1
        x = solution%x_lower + &
          (real(cell, dp) - 0.5_dp) * solution%coarse_dx
        call write_amr_cell( &
          unit, 0, solution%coarse_dx, solution%time, x, species, &
          solution%coarse(:, cell), solution%coarse_temperature(cell), q, &
          local_ok)
        if (.not. local_ok) then
          close(unit)
          return
        end if
      end do
      fine_cells = solution%hierarchy%fine%cell_count()
      do cell = 1, fine_cells
        global_fine = solution%hierarchy%fine%lower + cell - 1
        x = solution%x_lower + &
          (real(global_fine, dp) - 0.5_dp) * solution%hierarchy%fine_dx
        call write_amr_cell( &
          unit, 1, solution%hierarchy%fine_dx, solution%time, x, species, &
          solution%fine(:, cell), solution%fine_temperature(cell), q, &
          local_ok)
        if (.not. local_ok) then
          close(unit)
          return
        end if
      end do
      do cell = solution%hierarchy%fine_coarse_upper + 1, nx
        x = solution%x_lower + &
          (real(cell, dp) - 0.5_dp) * solution%coarse_dx
        call write_amr_cell( &
          unit, 0, solution%coarse_dx, solution%time, x, species, &
          solution%coarse(:, cell), solution%coarse_temperature(cell), q, &
          local_ok)
        if (.not. local_ok) then
          close(unit)
          return
        end if
      end do
    end if
    close(unit)
    ok = .true.
  end subroutine write_amr_reactive_1d_csv

  pure logical function valid_amr_configuration(config, nvar) result(valid)
    type(reactive_1d_config), intent(in) :: config
    integer, intent(in) :: nvar

    valid = config%amr_enabled .and. .not. config%transport_enabled .and. &
      trim(config%reconstruction) == "pcm" .and. config%nx >= 8 .and. &
      config%amr_refinement_ratio >= 2 .and. &
      config%amr_regrid_interval >= 1 .and. &
      config%amr_tag_component >= 1 .and. &
      config%amr_tag_component <= nvar .and. &
      config%amr_buffer_cells >= 0 .and. &
      config%amr_minimum_patch_cells >= 1 .and. &
      config%amr_minimum_patch_cells <= config%nx - 2
  end function valid_amr_configuration

  pure subroutine criteria_from_config(config, criteria)
    type(reactive_1d_config), intent(in) :: config
    type(amr_tagging_criteria_1d), intent(out) :: criteria

    criteria%component = config%amr_tag_component
    criteria%relative_gradient_threshold = &
      config%amr_relative_gradient_threshold
    criteria%absolute_gradient_threshold = &
      config%amr_absolute_gradient_threshold
    criteria%scale_floor = config%amr_scale_floor
    criteria%buffer_cells = config%amr_buffer_cells
    criteria%minimum_patch_cells = config%amr_minimum_patch_cells
  end subroutine criteria_from_config

  subroutine apply_regrid_plan_1d( &
      species, config, plan, solution, changed, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_regrid_plan_1d), intent(in) :: plan
    type(amr_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: changed, ok

    type(amr_two_level_hierarchy_1d) :: new_hierarchy
    real(dp), allocatable :: coarse_interior(:, :), old_fine(:, :)
    real(dp), allocatable :: new_fine(:, :)
    logical :: old_active, local_ok
    integer :: nvar, fine_cells

    changed = .false.
    ok = .false.
    old_active = solution%fine_active()
    if (old_active .eqv. plan%active) then
      if (.not. old_active) then
        ok = .true.
        return
      end if
      if (solution%hierarchy%fine_coarse_lower == plan%patch_lower .and. &
          solution%hierarchy%fine_coarse_upper == plan%patch_upper) then
        ok = .true.
        return
      end if
    end if

    nvar = size(solution%coarse, 1)
    allocate(coarse_interior(nvar, config%nx))
    coarse_interior = solution%coarse(:, 1:config%nx)
    if (old_active) then
      fine_cells = solution%hierarchy%fine%cell_count()
      allocate(old_fine(nvar, fine_cells))
      old_fine = solution%fine(:, 1:fine_cells)
    end if
    call regrid_two_level_state_1d( &
      coarse_interior, solution%hierarchy, old_fine, plan, &
      config%amr_refinement_ratio, config%x_lower, config%x_upper, &
      new_hierarchy, new_fine, local_ok)
    if (.not. local_ok) return
    solution%coarse(:, 1:config%nx) = coarse_interior
    if (allocated(solution%fine)) deallocate(solution%fine)
    if (allocated(solution%fine_temperature)) &
      deallocate(solution%fine_temperature)
    solution%hierarchy = new_hierarchy
    if (plan%active) then
      fine_cells = new_hierarchy%fine%cell_count()
      allocate(solution%fine(nvar, 0:fine_cells + 1))
      allocate(solution%fine_temperature(0:fine_cells + 1))
      solution%fine = 0.0_dp
      solution%fine(:, 1:fine_cells) = new_fine
      solution%fine_temperature = config%initial_temperature
      call recover_level_temperatures_1d( &
        species, solution%fine, solution%fine_temperature, fine_cells, &
        local_ok)
      if (.not. local_ok) return
    end if
    call recover_level_temperatures_1d( &
      species, solution%coarse, solution%coarse_temperature, config%nx, &
      local_ok)
    if (.not. local_ok) return
    call fill_physical_ghosts_1d( &
      solution%coarse, solution%coarse_temperature, config%nx, &
      config%boundary_condition, local_ok)
    if (.not. local_ok) return
    if (solution%fine_active()) then
      call fill_fine_ghosts_1d( &
        species, solution%hierarchy, solution%coarse, solution%coarse, &
        0.0_dp, solution%fine, solution%fine_temperature, local_ok)
      if (.not. local_ok) return
    end if
    changed = .true.
    ok = .true.
  end subroutine apply_regrid_plan_1d

  subroutine advance_pcm_level_1d( &
      species, state, temperature, nx, dx, dt, riemann_solver, &
      physical_boundary, boundary, flux, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, dt
    character(len=*), intent(in) :: riemann_solver
    logical, intent(in) :: physical_boundary
    character(len=*), intent(in) :: boundary
    real(dp), intent(out) :: flux(:, 0:)
    logical, intent(out) :: ok

    real(dp), allocatable :: old_state(:, :), old_temperature(:), q(:)
    real(dp) :: local_temperature, sound_speed
    logical :: local_ok
    integer :: face, cell

    ok = .false.
    if (nx < 1 .or. dx <= 0.0_dp .or. dt <= 0.0_dp .or. &
        size(state, 2) /= nx + 2 .or. size(temperature) /= nx + 2 .or. &
        size(flux, 1) /= size(state, 1) .or. &
        size(flux, 2) /= nx + 1) return
    if (physical_boundary) then
      call fill_physical_ghosts_1d( &
        state, temperature, nx, boundary, local_ok)
      if (.not. local_ok) return
    end if
    allocate(old_state(size(state, 1), 0:nx + 1))
    allocate(old_temperature(0:nx + 1))
    old_state = state
    old_temperature = temperature
    allocate(q(reactive_nprim(size(species))))
    do face = 0, nx
      call reactive_riemann_flux_x( &
        species, old_state(:, face), old_state(:, face + 1), &
        old_temperature(face), old_temperature(face + 1), riemann_solver, &
        flux(:, face), local_ok)
      if (.not. local_ok) return
    end do
    do cell = 1, nx
      state(:, cell) = old_state(:, cell) - dt / dx * &
        (flux(:, cell) - flux(:, cell - 1))
      call reactive_conserved_to_primitive( &
        species, state(:, cell), old_temperature(cell), q, &
        local_temperature, sound_speed, local_ok)
      if (.not. local_ok) return
      temperature(cell) = local_temperature
    end do
    if (physical_boundary) then
      call fill_physical_ghosts_1d( &
        state, temperature, nx, boundary, local_ok)
      if (.not. local_ok) return
    end if
    ok = .true.
  end subroutine advance_pcm_level_1d

  subroutine synchronize_levels_1d(species, config, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok

    logical :: local_ok
    integer :: fine_cells

    ok = .false.
    if (.not. solution%fine_active()) then
      ok = .true.
      return
    end if
    fine_cells = solution%hierarchy%fine%cell_count()
    call average_down_1d( &
      solution%fine(:, 1:fine_cells), solution%hierarchy, &
      solution%coarse(:, 1:config%nx), local_ok)
    if (.not. local_ok) return
    call recover_level_temperatures_1d( &
      species, solution%coarse, solution%coarse_temperature, config%nx, &
      local_ok)
    if (.not. local_ok) return
    call fill_physical_ghosts_1d( &
      solution%coarse, solution%coarse_temperature, config%nx, &
      config%boundary_condition, local_ok)
    if (.not. local_ok) return
    ok = .true.
  end subroutine synchronize_levels_1d

  subroutine recover_level_temperatures_1d( &
      species, state, temperature, nx, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, 0:)
    real(dp), intent(inout) :: temperature(0:)
    integer, intent(in) :: nx
    logical, intent(out) :: ok

    real(dp), allocatable :: q(:)
    real(dp) :: local_temperature, sound_speed, guess
    logical :: local_ok
    integer :: cell

    ok = .false.
    allocate(q(reactive_nprim(size(species))))
    do cell = 1, nx
      guess = temperature(cell)
      if (guess <= 0.0_dp) guess = 1000.0_dp
      call reactive_conserved_to_primitive( &
        species, state(:, cell), guess, q, local_temperature, sound_speed, &
        local_ok)
      if (.not. local_ok) return
      temperature(cell) = local_temperature
    end do
    ok = .true.
  end subroutine recover_level_temperatures_1d

  pure subroutine fill_physical_ghosts_1d( &
      state, temperature, nx, boundary, ok)
    real(dp), intent(inout) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    character(len=*), intent(in) :: boundary
    logical, intent(out) :: ok

    ok = .true.
    select case (trim(boundary))
    case ("periodic")
      state(:, 0) = state(:, nx)
      state(:, nx + 1) = state(:, 1)
      temperature(0) = temperature(nx)
      temperature(nx + 1) = temperature(1)
    case ("outflow")
      state(:, 0) = state(:, 1)
      state(:, nx + 1) = state(:, nx)
      temperature(0) = temperature(1)
      temperature(nx + 1) = temperature(nx)
    case default
      ok = .false.
    end select
  end subroutine fill_physical_ghosts_1d

  subroutine fill_fine_ghosts_1d( &
      species, hierarchy, coarse_start, coarse_end, alpha, fine, &
      fine_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_two_level_hierarchy_1d), intent(in) :: hierarchy
    real(dp), intent(in) :: coarse_start(:, 0:), coarse_end(:, 0:)
    real(dp), intent(in) :: alpha
    real(dp), intent(inout) :: fine(:, 0:), fine_temperature(0:)
    logical, intent(out) :: ok

    real(dp), allocatable :: q(:)
    real(dp) :: guess, local_temperature, sound_speed
    logical :: local_ok
    integer :: fine_cells, coarse_cell, side

    ok = .false.
    if (.not. hierarchy%is_valid() .or. alpha < 0.0_dp .or. &
        alpha > 1.0_dp) return
    fine_cells = hierarchy%fine%cell_count()
    allocate(q(reactive_nprim(size(species))))
    do side = 1, 2
      if (side == 1) then
        coarse_cell = hierarchy%fine_coarse_lower - 1
        fine(:, 0) = (1.0_dp - alpha) * &
          coarse_start(:, coarse_cell) + alpha * coarse_end(:, coarse_cell)
        guess = fine_temperature(1)
        call reactive_conserved_to_primitive( &
          species, fine(:, 0), guess, q, local_temperature, sound_speed, &
          local_ok)
        if (.not. local_ok) return
        fine_temperature(0) = local_temperature
      else
        coarse_cell = hierarchy%fine_coarse_upper + 1
        fine(:, fine_cells + 1) = (1.0_dp - alpha) * &
          coarse_start(:, coarse_cell) + alpha * coarse_end(:, coarse_cell)
        guess = fine_temperature(fine_cells)
        call reactive_conserved_to_primitive( &
          species, fine(:, fine_cells + 1), guess, q, local_temperature, &
          sound_speed, local_ok)
        if (.not. local_ok) return
        fine_temperature(fine_cells + 1) = local_temperature
      end if
    end do
    ok = .true.
  end subroutine fill_fine_ghosts_1d

  subroutine write_amr_cell( &
      unit, level, cell_dx, time, x, species, state, temperature_guess, q, ok)
    integer, intent(in) :: unit, level
    real(dp), intent(in) :: cell_dx, time, x
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:), temperature_guess
    real(dp), intent(out) :: q(:)
    logical, intent(out) :: ok

    real(dp) :: temperature, sound_speed
    integer :: k

    call reactive_conserved_to_primitive( &
      species, state, temperature_guess, q, temperature, sound_speed, ok)
    if (.not. ok) return
    write(unit, '(i0,",",*(es25.16e3,:,","))') &
      level, cell_dx, time, x, state(irho), q(2), q(3), q(4), q(5), &
      temperature, state(iet), &
      (q(reactive_mass_fraction_component(k)), k = 1, size(species))
  end subroutine write_amr_cell

end module amr_reactive_1d_mod
