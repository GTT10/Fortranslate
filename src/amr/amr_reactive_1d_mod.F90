module amr_reactive_1d_mod
  use precision_mod, only: dp
  use constants_mod, only: density_floor, pressure_floor
  use state_indices_mod, only: irho, imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: gas_transport_species
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_conserved_to_primitive, reactive_primitive_to_conserved, &
    reactive_riemann_flux_x, &
    reconstruct_ppm_faces, reconstruct_characteristic_ppm_faces, &
    reactive_cfl_timestep, initialize_reactive_1d, &
    advance_reactive_chemistry, reactive_diffusive_flux_x, &
    reactive_transport_timestep
  use amr_hierarchy_1d_mod, only: &
    amr_two_level_hierarchy_1d, amr_flux_register_1d, &
    initialize_flux_register_1d, accumulate_coarse_flux_1d, &
    accumulate_fine_flux_1d, reflux_1d, average_down_1d, &
    composite_integral_1d
  use amr_regrid_1d_mod, only: &
    amr_tagging_criteria_1d, amr_regrid_plan_1d, &
    plan_gradient_regrid_1d, regrid_two_level_state_1d
  use slope_limiter_mod, only: limited_slope
  implicit none
  private

  ! Reconstructing both cells adjacent to an AMR interface requires the
  ! four exterior cell averages consumed by the widest characteristic stencil.
  integer, parameter, public :: amr_ppm_ghost_width = 4

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
  public :: advance_amr_level_1d
  public :: advance_transport_level_1d
  public :: recover_level_temperatures_1d
  public :: fill_physical_ghosts_1d
  public :: fill_fine_ghosts_1d
  public :: fill_physical_wide_ghosts_1d
  public :: fill_fine_wide_ghosts_1d
  public :: write_amr_cell

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
      species, config, solution, dt, ok, transport)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_reactive_solution_1d), intent(in) :: solution
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok
    type(gas_transport_species), intent(in), optional :: transport(:)

    real(dp) :: coarse_dt, fine_dt, transport_dt, fine_transport_dt
    real(dp) :: maximum_diffusivity
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
    if (config%transport_enabled) then
      if (.not. present(transport)) return
      call reactive_transport_timestep( &
        species, transport, solution%coarse, solution%coarse_temperature, &
        config%nx, solution%coarse_dx, config%transport_cfl, &
        config%viscosity_enabled, config%thermal_conduction_enabled, &
        config%species_diffusion_enabled, transport_dt, &
        maximum_diffusivity, local_ok)
      if (.not. local_ok) return
      dt = min(dt, transport_dt)
      if (solution%fine_active()) then
        call reactive_transport_timestep( &
          species, transport, solution%fine, solution%fine_temperature, &
          solution%hierarchy%fine%cell_count(), &
          solution%hierarchy%fine_dx, config%transport_cfl, &
          config%viscosity_enabled, config%thermal_conduction_enabled, &
          config%species_diffusion_enabled, fine_transport_dt, &
          maximum_diffusivity, local_ok)
        if (.not. local_ok) return
        dt = min(dt, real(solution%hierarchy%refinement_ratio**2, dp) * &
          fine_transport_dt)
      end if
    end if
    ok = dt > 0.0_dp
  end subroutine amr_reactive_timestep_1d

  subroutine advance_amr_reactive_1d( &
      species, reactions, config, dt, solution, ok, transport)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: dt
    type(amr_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok
    type(gas_transport_species), intent(in), optional :: transport(:)

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
    if (config%transport_enabled .and. .not. present(transport)) return
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

    if (config%transport_enabled) then
      call advance_amr_transport_1d( &
        species, transport, config, 0.5_dp * dt, solution, local_ok)
      if (.not. local_ok) then
        solution = backup
        return
      end if
    end if

    allocate(coarse_start( &
      size(solution%coarse, 1), 0:config%nx + 1))
    coarse_start = solution%coarse
    allocate(coarse_flux(reactive_nvar(size(species)), 0:config%nx))
    call advance_amr_level_1d( &
      species, solution%coarse, solution%coarse_temperature, config%nx, &
      solution%coarse_dx, dt, config%amr_reconstruction, config%limiter, &
      config%riemann_solver, .true., config%boundary_condition, &
      coarse_flux, local_ok)
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
        if (trim(config%amr_reconstruction) == "plm") then
          alpha = (real(substep, dp) - 0.5_dp) / real(ratio, dp)
        else
          alpha = real(substep - 1, dp) / real(ratio, dp)
        end if
        call fill_fine_ghosts_1d( &
          species, solution%hierarchy, coarse_start, coarse_end, alpha, &
          solution%fine, solution%fine_temperature, local_ok)
        if (.not. local_ok) then
          solution = backup
          return
        end if
        call advance_amr_level_1d( &
          species, solution%fine, solution%fine_temperature, fine_cells, &
          solution%hierarchy%fine_dx, fine_dt, &
          config%amr_reconstruction, config%limiter, &
          config%riemann_solver, .false., "coarse_fine", fine_flux, &
          local_ok)
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

    if (config%transport_enabled) then
      call advance_amr_transport_1d( &
        species, transport, config, 0.5_dp * dt, solution, local_ok)
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
      final_integrals, ok, transport)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_reactive_solution_1d), intent(out) :: solution
    real(dp), intent(out) :: initial_integrals(5), final_integrals(5)
    logical, intent(out) :: ok
    type(gas_transport_species), intent(in), optional :: transport(:)

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
      if (config%transport_enabled) then
        if (.not. present(transport)) then
          ok = .false.
          return
        end if
        call amr_reactive_timestep_1d( &
          species, config, solution, dt, local_ok, transport)
      else
        call amr_reactive_timestep_1d( &
          species, config, solution, dt, local_ok)
      end if
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      dt = min(dt, config%final_time - solution%time)
      if (config%transport_enabled) then
        call advance_amr_reactive_1d( &
          species, reactions, config, dt, solution, local_ok, transport)
      else
        call advance_amr_reactive_1d( &
          species, reactions, config, dt, solution, local_ok)
      end if
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

    valid = config%amr_enabled .and. config%nx >= 8 .and. &
      (trim(config%amr_reconstruction) == "pcm" .or. &
        trim(config%amr_reconstruction) == "plm") .and. &
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

  subroutine advance_amr_transport_1d( &
      species, transport, config, interval, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: interval
    type(amr_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok

    type(amr_flux_register_1d) :: flux_register
    real(dp), allocatable :: coarse_start(:, :), coarse_end(:, :)
    real(dp), allocatable :: coarse_flux(:, :), fine_flux(:, :)
    real(dp) :: fine_interval, alpha, interface_distance
    logical :: local_ok
    integer :: fine_cells, ratio, subcycles, substep, nvar

    ok = .false.
    if (interval < 0.0_dp .or. size(transport) /= size(species)) return
    if (interval == 0.0_dp .or. .not. &
        (config%viscosity_enabled .or. &
          config%thermal_conduction_enabled .or. &
          config%species_diffusion_enabled)) then
      ok = .true.
      return
    end if

    nvar = reactive_nvar(size(species))
    allocate(coarse_start(nvar, 0:config%nx + 1))
    allocate(coarse_end(nvar, 0:config%nx + 1))
    allocate(coarse_flux(nvar, 0:config%nx))
    coarse_start = solution%coarse
    call advance_transport_level_1d( &
      species, transport, solution%coarse, solution%coarse_temperature, &
      config%nx, solution%coarse_dx, interval, solution%coarse_dx, &
      config, .true., config%boundary_condition, coarse_flux, local_ok)
    if (.not. local_ok) return
    coarse_end = solution%coarse

    if (solution%fine_active()) then
      ratio = solution%hierarchy%refinement_ratio
      subcycles = ratio * ratio
      fine_cells = solution%hierarchy%fine%cell_count()
      fine_interval = interval / real(subcycles, dp)
      interface_distance = 0.5_dp * &
        (solution%coarse_dx + solution%hierarchy%fine_dx)
      allocate(fine_flux(nvar, 0:fine_cells))
      call initialize_flux_register_1d(flux_register, nvar, local_ok)
      if (.not. local_ok) return
      call accumulate_coarse_flux_1d( &
        flux_register, &
        coarse_flux(:, solution%hierarchy%fine_coarse_lower - 1), &
        coarse_flux(:, solution%hierarchy%fine_coarse_upper), interval, &
        local_ok)
      if (.not. local_ok) return

      do substep = 1, subcycles
        alpha = (real(substep, dp) - 0.5_dp) / real(subcycles, dp)
        call fill_fine_ghosts_1d( &
          species, solution%hierarchy, coarse_start, coarse_end, alpha, &
          solution%fine, solution%fine_temperature, local_ok)
        if (.not. local_ok) return
        call advance_transport_level_1d( &
          species, transport, solution%fine, solution%fine_temperature, &
          fine_cells, solution%hierarchy%fine_dx, fine_interval, &
          interface_distance, config, .false., "coarse_fine", fine_flux, &
          local_ok)
        if (.not. local_ok) return
        call accumulate_fine_flux_1d( &
          flux_register, fine_flux(:, 0), fine_flux(:, fine_cells), &
          fine_interval, local_ok)
        if (.not. local_ok) return
      end do
      call reflux_1d( &
        solution%coarse(:, 1:config%nx), solution%hierarchy, &
        flux_register, local_ok)
      if (.not. local_ok) return
      call synchronize_levels_1d(species, config, solution, local_ok)
      if (.not. local_ok) return
      call fill_fine_ghosts_1d( &
        species, solution%hierarchy, solution%coarse, solution%coarse, &
        1.0_dp, solution%fine, solution%fine_temperature, local_ok)
      if (.not. local_ok) return
    end if
    ok = .true.
  end subroutine advance_amr_transport_1d

  subroutine advance_transport_level_1d( &
      species, transport, state, temperature, nx, dx, interval, &
      boundary_distance, config, physical_boundary, boundary, flux, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(inout) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, interval, boundary_distance
    type(reactive_1d_config), intent(in) :: config
    logical, intent(in) :: physical_boundary
    character(len=*), intent(in) :: boundary
    real(dp), intent(out) :: flux(:, 0:)
    logical, intent(out) :: ok

    real(dp), allocatable :: initial_state(:, :), initial_temperature(:)
    real(dp), allocatable :: stage1_state(:, :), stage1_temperature(:)
    real(dp), allocatable :: euler2_state(:, :), euler2_temperature(:)
    real(dp), allocatable :: first_flux(:, :), second_flux(:, :), q(:)
    real(dp) :: guess, local_temperature, sound_speed
    logical :: local_ok
    integer :: nvar, cell

    ok = .false.
    if (interval <= 0.0_dp .or. boundary_distance <= 0.0_dp) return
    nvar = size(state, 1)
    allocate(initial_state(nvar, 0:nx + 1))
    allocate(initial_temperature(0:nx + 1))
    allocate(stage1_state(nvar, 0:nx + 1))
    allocate(stage1_temperature(0:nx + 1))
    allocate(euler2_state(nvar, 0:nx + 1))
    allocate(euler2_temperature(0:nx + 1))
    allocate(first_flux(nvar, 0:nx), second_flux(nvar, 0:nx))
    allocate(q(reactive_nprim(size(species))))
    initial_state = state
    initial_temperature = temperature
    if (physical_boundary) then
      call fill_physical_ghosts_1d( &
        initial_state, initial_temperature, nx, boundary, local_ok)
      if (.not. local_ok) return
    end if

    call transport_euler_update_1d( &
      species, transport, initial_state, initial_temperature, nx, dx, &
      interval, boundary_distance, config, physical_boundary, stage1_state, &
      stage1_temperature, first_flux, local_ok)
    if (.not. local_ok) return
    if (physical_boundary) then
      call fill_physical_ghosts_1d( &
        stage1_state, stage1_temperature, nx, boundary, local_ok)
      if (.not. local_ok) return
    end if
    call transport_euler_update_1d( &
      species, transport, stage1_state, stage1_temperature, nx, dx, &
      interval, boundary_distance, config, physical_boundary, euler2_state, &
      euler2_temperature, second_flux, local_ok)
    if (.not. local_ok) return

    state = initial_state
    temperature = initial_temperature
    do cell = 1, nx
      state(:, cell) = &
        0.5_dp * (initial_state(:, cell) + euler2_state(:, cell))
      guess = 0.5_dp * &
        (initial_temperature(cell) + euler2_temperature(cell))
      call reactive_conserved_to_primitive( &
        species, state(:, cell), guess, q, local_temperature, sound_speed, &
        local_ok)
      if (.not. local_ok) return
      temperature(cell) = local_temperature
    end do
    flux = 0.5_dp * (first_flux + second_flux)
    if (physical_boundary) then
      call fill_physical_ghosts_1d( &
        state, temperature, nx, boundary, local_ok)
      if (.not. local_ok) return
    end if
    ok = .true.
  end subroutine advance_transport_level_1d

  subroutine transport_euler_update_1d( &
      species, transport, state, temperature, nx, dx, interval, &
      boundary_distance, config, physical_boundary, output_state, &
      output_temperature, flux, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, interval, boundary_distance
    type(reactive_1d_config), intent(in) :: config
    logical, intent(in) :: physical_boundary
    real(dp), intent(out) :: output_state(:, 0:), output_temperature(0:)
    real(dp), intent(out) :: flux(:, 0:)
    logical, intent(out) :: ok

    real(dp), allocatable :: q(:)
    real(dp) :: face_distance, local_temperature, sound_speed
    logical :: local_ok
    integer :: face, cell

    ok = .false.
    allocate(q(reactive_nprim(size(species))))
    do face = 0, nx
      face_distance = dx
      if (.not. physical_boundary .and. &
          (face == 0 .or. face == nx)) face_distance = boundary_distance
      call reactive_diffusive_flux_x( &
        species, transport, state(:, face), state(:, face + 1), &
        temperature(face), temperature(face + 1), face_distance, &
        config%viscosity_enabled, config%thermal_conduction_enabled, &
        config%species_diffusion_enabled, config%barodiffusion_enabled, &
        flux(:, face), local_ok)
      if (.not. local_ok) return
    end do
    output_state = state
    output_temperature = temperature
    do cell = 1, nx
      output_state(:, cell) = state(:, cell) - interval / dx * &
        (flux(:, cell) - flux(:, cell - 1))
      call reactive_conserved_to_primitive( &
        species, output_state(:, cell), temperature(cell), q, &
        local_temperature, sound_speed, local_ok)
      if (.not. local_ok) return
      output_temperature(cell) = local_temperature
    end do
    ok = .true.
  end subroutine transport_euler_update_1d

  subroutine advance_amr_level_1d( &
      species, state, temperature, nx, dx, dt, reconstruction, limiter, &
      riemann_solver, physical_boundary, boundary, flux, ok, &
      left_ghost_state, right_ghost_state, left_ghost_temperature, &
      right_ghost_temperature, ppm_contact_steepening, &
      ppm_shock_flattening)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, dt
    character(len=*), intent(in) :: reconstruction, limiter
    character(len=*), intent(in) :: riemann_solver
    logical, intent(in) :: physical_boundary
    character(len=*), intent(in) :: boundary
    real(dp), intent(out) :: flux(:, 0:)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: left_ghost_state(:, :)
    real(dp), intent(in), optional :: right_ghost_state(:, :)
    real(dp), intent(in), optional :: left_ghost_temperature(:)
    real(dp), intent(in), optional :: right_ghost_temperature(:)
    logical, intent(in), optional :: ppm_contact_steepening
    logical, intent(in), optional :: ppm_shock_flattening

    logical :: use_contact_steepening, use_shock_flattening, wide_ghosts

    use_contact_steepening = .false.
    use_shock_flattening = .false.
    if (present(ppm_contact_steepening)) &
      use_contact_steepening = ppm_contact_steepening
    if (present(ppm_shock_flattening)) &
      use_shock_flattening = ppm_shock_flattening
    wide_ghosts = present(left_ghost_state) .and. &
      present(right_ghost_state) .and. &
      present(left_ghost_temperature) .and. &
      present(right_ghost_temperature)
    if (wide_ghosts .neqv. (present(left_ghost_state) .or. &
        present(right_ghost_state) .or. &
        present(left_ghost_temperature) .or. &
        present(right_ghost_temperature))) then
      flux = 0.0_dp
      ok = .false.
      return
    end if

    select case (trim(reconstruction))
    case ("pcm")
      call advance_pcm_level_1d( &
        species, state, temperature, nx, dx, dt, riemann_solver, &
        physical_boundary, boundary, flux, ok)
    case ("plm")
      call advance_plm_level_1d( &
        species, state, temperature, nx, dx, dt, limiter, riemann_solver, &
        physical_boundary, boundary, flux, ok)
    case ("ppm", "characteristic_ppm")
      if (physical_boundary) then
        call advance_ppm_level_1d( &
          species, state, temperature, nx, dx, dt, reconstruction, &
          riemann_solver, boundary, use_contact_steepening, &
          use_shock_flattening, flux, ok)
      else if (wide_ghosts) then
        call advance_ppm_level_1d( &
          species, state, temperature, nx, dx, dt, reconstruction, &
          riemann_solver, boundary, use_contact_steepening, &
          use_shock_flattening, flux, ok, left_ghost_state, &
          right_ghost_state, left_ghost_temperature, &
          right_ghost_temperature)
      else
        flux = 0.0_dp
        ok = .false.
      end if
    case default
      flux = 0.0_dp
      ok = .false.
    end select
  end subroutine advance_amr_level_1d

  subroutine advance_ppm_level_1d( &
      species, state, temperature, nx, dx, dt, reconstruction, &
      riemann_solver, boundary, use_contact_steepening, &
      use_shock_flattening, flux, ok, left_ghost_state, &
      right_ghost_state, left_ghost_temperature, right_ghost_temperature)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, dt
    character(len=*), intent(in) :: reconstruction, riemann_solver, boundary
    logical, intent(in) :: use_contact_steepening, use_shock_flattening
    real(dp), intent(out) :: flux(:, 0:)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: left_ghost_state(:, :)
    real(dp), intent(in), optional :: right_ghost_state(:, :)
    real(dp), intent(in), optional :: left_ghost_temperature(:)
    real(dp), intent(in), optional :: right_ghost_temperature(:)

    real(dp), allocatable :: initial_state(:, :), initial_temperature(:)
    real(dp), allocatable :: stage1_state(:, :), stage1_temperature(:)
    real(dp), allocatable :: euler2_state(:, :), euler2_temperature(:)
    real(dp), allocatable :: stage2_state(:, :), stage2_temperature(:)
    real(dp), allocatable :: euler3_state(:, :), euler3_temperature(:)
    real(dp), allocatable :: first_flux(:, :), second_flux(:, :)
    real(dp), allocatable :: third_flux(:, :), q(:)
    real(dp), allocatable :: left_wide(:, :), right_wide(:, :)
    real(dp), allocatable :: left_wide_temperature(:)
    real(dp), allocatable :: right_wide_temperature(:)
    real(dp) :: guess, local_temperature, sound_speed
    logical :: local_ok, external_ghosts
    integer :: nvar, cell

    ok = .false.
    external_ghosts = present(left_ghost_state) .and. &
      present(right_ghost_state) .and. &
      present(left_ghost_temperature) .and. &
      present(right_ghost_temperature)
    if (external_ghosts .neqv. (present(left_ghost_state) .or. &
        present(right_ghost_state) .or. &
        present(left_ghost_temperature) .or. &
        present(right_ghost_temperature))) return
    nvar = size(state, 1)
    allocate(initial_state(nvar, 0:nx + 1))
    allocate(initial_temperature(0:nx + 1))
    allocate(stage1_state(nvar, 0:nx + 1))
    allocate(stage1_temperature(0:nx + 1))
    allocate(euler2_state(nvar, 0:nx + 1))
    allocate(euler2_temperature(0:nx + 1))
    allocate(stage2_state(nvar, 0:nx + 1))
    allocate(stage2_temperature(0:nx + 1))
    allocate(euler3_state(nvar, 0:nx + 1))
    allocate(euler3_temperature(0:nx + 1))
    allocate(first_flux(nvar, 0:nx), second_flux(nvar, 0:nx))
    allocate(third_flux(nvar, 0:nx), q(reactive_nprim(size(species))))
    allocate(left_wide(nvar, amr_ppm_ghost_width))
    allocate(right_wide(nvar, amr_ppm_ghost_width))
    allocate(left_wide_temperature(amr_ppm_ghost_width))
    allocate(right_wide_temperature(amr_ppm_ghost_width))
    initial_state = state
    initial_temperature = temperature

    call prepare_ppm_ghosts( &
      initial_state, initial_temperature, nx, boundary, external_ghosts, &
      left_wide, right_wide, left_wide_temperature, &
      right_wide_temperature, local_ok, left_ghost_state, &
      right_ghost_state, left_ghost_temperature, right_ghost_temperature)
    if (.not. local_ok) return
    call ppm_euler_update_1d( &
      species, initial_state, initial_temperature, nx, dx, dt, &
      reconstruction, riemann_solver, boundary, use_contact_steepening, &
      use_shock_flattening, left_wide, right_wide, &
      left_wide_temperature, right_wide_temperature, stage1_state, &
      stage1_temperature, first_flux, local_ok)
    if (.not. local_ok) return

    call prepare_ppm_ghosts( &
      stage1_state, stage1_temperature, nx, boundary, external_ghosts, &
      left_wide, right_wide, left_wide_temperature, &
      right_wide_temperature, local_ok, left_ghost_state, &
      right_ghost_state, left_ghost_temperature, right_ghost_temperature)
    if (.not. local_ok) return
    call ppm_euler_update_1d( &
      species, stage1_state, stage1_temperature, nx, dx, dt, &
      reconstruction, riemann_solver, boundary, use_contact_steepening, &
      use_shock_flattening, left_wide, right_wide, &
      left_wide_temperature, right_wide_temperature, euler2_state, &
      euler2_temperature, second_flux, local_ok)
    if (.not. local_ok) return
    stage2_state = initial_state
    stage2_temperature = initial_temperature
    do cell = 1, nx
      stage2_state(:, cell) = &
        0.75_dp * initial_state(:, cell) + &
        0.25_dp * euler2_state(:, cell)
      guess = 0.75_dp * initial_temperature(cell) + &
        0.25_dp * euler2_temperature(cell)
      call reactive_conserved_to_primitive( &
        species, stage2_state(:, cell), guess, q, local_temperature, &
        sound_speed, local_ok)
      if (.not. local_ok) return
      stage2_temperature(cell) = local_temperature
    end do

    call prepare_ppm_ghosts( &
      stage2_state, stage2_temperature, nx, boundary, external_ghosts, &
      left_wide, right_wide, left_wide_temperature, &
      right_wide_temperature, local_ok, left_ghost_state, &
      right_ghost_state, left_ghost_temperature, right_ghost_temperature)
    if (.not. local_ok) return
    call ppm_euler_update_1d( &
      species, stage2_state, stage2_temperature, nx, dx, dt, &
      reconstruction, riemann_solver, boundary, use_contact_steepening, &
      use_shock_flattening, left_wide, right_wide, &
      left_wide_temperature, right_wide_temperature, euler3_state, &
      euler3_temperature, third_flux, local_ok)
    if (.not. local_ok) return

    state = initial_state
    temperature = initial_temperature
    do cell = 1, nx
      state(:, cell) = initial_state(:, cell) / 3.0_dp + &
        2.0_dp * euler3_state(:, cell) / 3.0_dp
      guess = initial_temperature(cell) / 3.0_dp + &
        2.0_dp * euler3_temperature(cell) / 3.0_dp
      call reactive_conserved_to_primitive( &
        species, state(:, cell), guess, q, local_temperature, sound_speed, &
        local_ok)
      if (.not. local_ok) return
      temperature(cell) = local_temperature
    end do
    ! This is the conservative flux combination represented by SSPRK3:
    ! U^(n+1) = U^n - dt*div((F0 + F1)/6 + 2*F2/3).
    flux = (first_flux + second_flux) / 6.0_dp + &
      2.0_dp * third_flux / 3.0_dp
    if (.not. external_ghosts) then
      call fill_physical_ghosts_1d( &
        state, temperature, nx, boundary, local_ok)
      if (.not. local_ok) return
    end if
    ok = .true.
  end subroutine advance_ppm_level_1d

  subroutine prepare_ppm_ghosts( &
      state, temperature, nx, boundary, external_ghosts, left_wide, &
      right_wide, left_wide_temperature, right_wide_temperature, ok, &
      external_left, external_right, external_left_temperature, &
      external_right_temperature)
    real(dp), intent(in) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    character(len=*), intent(in) :: boundary
    logical, intent(in) :: external_ghosts
    real(dp), intent(out) :: left_wide(:, :), right_wide(:, :)
    real(dp), intent(out) :: left_wide_temperature(:)
    real(dp), intent(out) :: right_wide_temperature(:)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: external_left(:, :)
    real(dp), intent(in), optional :: external_right(:, :)
    real(dp), intent(in), optional :: external_left_temperature(:)
    real(dp), intent(in), optional :: external_right_temperature(:)

    ok = .false.
    if (external_ghosts) then
      if (.not. present(external_left) .or. &
          .not. present(external_right) .or. &
          .not. present(external_left_temperature) .or. &
          .not. present(external_right_temperature)) return
      if (size(external_left, 1) /= size(left_wide, 1) .or. &
          size(external_right, 1) /= size(right_wide, 1) .or. &
          size(external_left, 2) /= amr_ppm_ghost_width .or. &
          size(external_right, 2) /= amr_ppm_ghost_width .or. &
          size(external_left_temperature) /= amr_ppm_ghost_width .or. &
          size(external_right_temperature) /= amr_ppm_ghost_width) return
      left_wide = external_left
      right_wide = external_right
      left_wide_temperature = external_left_temperature
      right_wide_temperature = external_right_temperature
      ok = .true.
      return
    end if
    call fill_physical_wide_ghosts_1d( &
      state, temperature, nx, boundary, left_wide, right_wide, &
      left_wide_temperature, right_wide_temperature, ok)
  end subroutine prepare_ppm_ghosts

  subroutine ppm_euler_update_1d( &
      species, state, temperature, nx, dx, dt, reconstruction, &
      riemann_solver, boundary, use_contact_steepening, &
      use_shock_flattening, left_wide, right_wide, &
      left_wide_temperature, right_wide_temperature, output_state, &
      output_temperature, flux, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, dt
    character(len=*), intent(in) :: reconstruction, riemann_solver, boundary
    logical, intent(in) :: use_contact_steepening, use_shock_flattening
    real(dp), intent(in) :: left_wide(:, :), right_wide(:, :)
    real(dp), intent(in) :: left_wide_temperature(:)
    real(dp), intent(in) :: right_wide_temperature(:)
    real(dp), intent(out) :: output_state(:, 0:), output_temperature(0:)
    real(dp), intent(out) :: flux(:, 0:)
    logical, intent(out) :: ok

    real(dp), allocatable :: left_face(:, :), right_face(:, :)
    real(dp), allocatable :: left_t(:), right_t(:), q(:)
    real(dp) :: local_temperature, sound_speed
    logical :: local_ok
    integer :: nvar, face, cell

    ok = .false.
    nvar = size(state, 1)
    allocate(left_face(nvar, 0:nx), right_face(nvar, 0:nx))
    allocate(left_t(0:nx), right_t(0:nx))
    allocate(q(reactive_nprim(size(species))))
    select case (trim(reconstruction))
    case ("ppm")
      call reconstruct_ppm_faces( &
        species, state, temperature, nx, boundary, left_face, right_face, &
        left_t, right_t, local_ok, left_wide, right_wide, &
        left_wide_temperature, right_wide_temperature)
    case ("characteristic_ppm")
      call reconstruct_characteristic_ppm_faces( &
        species, state, temperature, nx, boundary, dt / dx, &
        use_contact_steepening, use_shock_flattening, left_face, &
        right_face, left_t, right_t, local_ok, left_wide, right_wide, &
        left_wide_temperature, right_wide_temperature)
    case default
      return
    end select
    if (.not. local_ok) return
    do face = 0, nx
      call reactive_riemann_flux_x( &
        species, left_face(:, face), right_face(:, face), &
        left_t(face), right_t(face), riemann_solver, flux(:, face), local_ok)
      if (.not. local_ok) return
    end do
    output_state = state
    output_temperature = temperature
    do cell = 1, nx
      output_state(:, cell) = state(:, cell) - dt / dx * &
        (flux(:, cell) - flux(:, cell - 1))
      call reactive_conserved_to_primitive( &
        species, output_state(:, cell), temperature(cell), q, &
        local_temperature, sound_speed, local_ok)
      if (.not. local_ok) return
      output_temperature(cell) = local_temperature
    end do
    ok = .true.
  end subroutine ppm_euler_update_1d

  subroutine advance_plm_level_1d( &
      species, state, temperature, nx, dx, dt, limiter, riemann_solver, &
      physical_boundary, boundary, flux, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, dt
    character(len=*), intent(in) :: limiter, riemann_solver
    logical, intent(in) :: physical_boundary
    character(len=*), intent(in) :: boundary
    real(dp), intent(out) :: flux(:, 0:)
    logical, intent(out) :: ok

    real(dp), allocatable :: initial_state(:, :), initial_temperature(:)
    real(dp), allocatable :: stage1_state(:, :), stage1_temperature(:)
    real(dp), allocatable :: euler2_state(:, :), euler2_temperature(:)
    real(dp), allocatable :: first_flux(:, :), second_flux(:, :), q(:)
    real(dp) :: temperature_guess, local_temperature, sound_speed
    logical :: local_ok
    integer :: nvar, cell

    ok = .false.
    nvar = size(state, 1)
    allocate(initial_state(nvar, 0:nx + 1))
    allocate(initial_temperature(0:nx + 1))
    allocate(stage1_state(nvar, 0:nx + 1))
    allocate(stage1_temperature(0:nx + 1))
    allocate(euler2_state(nvar, 0:nx + 1))
    allocate(euler2_temperature(0:nx + 1))
    allocate(first_flux(nvar, 0:nx), second_flux(nvar, 0:nx))
    allocate(q(reactive_nprim(size(species))))
    initial_state = state
    initial_temperature = temperature
    if (physical_boundary) then
      call fill_physical_ghosts_1d( &
        initial_state, initial_temperature, nx, boundary, local_ok)
      if (.not. local_ok) return
    end if

    call plm_euler_update_1d( &
      species, initial_state, initial_temperature, nx, dx, dt, limiter, &
      riemann_solver, physical_boundary, boundary, stage1_state, &
      stage1_temperature, first_flux, local_ok)
    if (.not. local_ok) return
    if (physical_boundary) then
      call fill_physical_ghosts_1d( &
        stage1_state, stage1_temperature, nx, boundary, local_ok)
      if (.not. local_ok) return
    end if
    call plm_euler_update_1d( &
      species, stage1_state, stage1_temperature, nx, dx, dt, limiter, &
      riemann_solver, physical_boundary, boundary, euler2_state, &
      euler2_temperature, second_flux, local_ok)
    if (.not. local_ok) return

    state = initial_state
    temperature = initial_temperature
    do cell = 1, nx
      state(:, cell) = &
        0.5_dp * (initial_state(:, cell) + euler2_state(:, cell))
      temperature_guess = 0.5_dp * &
        (initial_temperature(cell) + euler2_temperature(cell))
      call reactive_conserved_to_primitive( &
        species, state(:, cell), temperature_guess, q, local_temperature, &
        sound_speed, local_ok)
      if (.not. local_ok) return
      temperature(cell) = local_temperature
    end do
    flux = 0.5_dp * (first_flux + second_flux)
    if (physical_boundary) then
      call fill_physical_ghosts_1d( &
        state, temperature, nx, boundary, local_ok)
      if (.not. local_ok) return
    end if
    ok = .true.
  end subroutine advance_plm_level_1d

  subroutine plm_euler_update_1d( &
      species, state, temperature, nx, dx, dt, limiter, riemann_solver, &
      physical_boundary, boundary, output_state, output_temperature, &
      flux, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, dt
    character(len=*), intent(in) :: limiter, riemann_solver
    logical, intent(in) :: physical_boundary
    character(len=*), intent(in) :: boundary
    real(dp), intent(out) :: output_state(:, 0:), output_temperature(0:)
    real(dp), intent(out) :: flux(:, 0:)
    logical, intent(out) :: ok

    real(dp), allocatable :: left_state(:, :), right_state(:, :)
    real(dp), allocatable :: left_temperature(:), right_temperature(:)
    real(dp), allocatable :: q(:)
    real(dp) :: local_temperature, sound_speed
    logical :: local_ok
    integer :: nvar, face, cell

    ok = .false.
    nvar = size(state, 1)
    allocate(left_state(nvar, nx), right_state(nvar, nx))
    allocate(left_temperature(nx), right_temperature(nx))
    allocate(q(reactive_nprim(size(species))))
    call build_plm_cell_edges_1d( &
      species, state, temperature, nx, limiter, left_state, right_state, &
      left_temperature, right_temperature, local_ok)
    if (.not. local_ok) return

    do face = 1, nx - 1
      call reactive_riemann_flux_x( &
        species, right_state(:, face), left_state(:, face + 1), &
        right_temperature(face), left_temperature(face + 1), &
        riemann_solver, flux(:, face), local_ok)
      if (.not. local_ok) return
    end do
    if (physical_boundary .and. trim(boundary) == "periodic") then
      call reactive_riemann_flux_x( &
        species, right_state(:, nx), left_state(:, 1), &
        right_temperature(nx), left_temperature(1), riemann_solver, &
        flux(:, 0), local_ok)
      if (.not. local_ok) return
      flux(:, nx) = flux(:, 0)
    else
      call reactive_riemann_flux_x( &
        species, state(:, 0), left_state(:, 1), temperature(0), &
        left_temperature(1), riemann_solver, flux(:, 0), local_ok)
      if (.not. local_ok) return
      call reactive_riemann_flux_x( &
        species, right_state(:, nx), state(:, nx + 1), &
        right_temperature(nx), temperature(nx + 1), riemann_solver, &
        flux(:, nx), local_ok)
      if (.not. local_ok) return
    end if

    output_state = state
    output_temperature = temperature
    do cell = 1, nx
      output_state(:, cell) = state(:, cell) - dt / dx * &
        (flux(:, cell) - flux(:, cell - 1))
      call reactive_conserved_to_primitive( &
        species, output_state(:, cell), temperature(cell), q, &
        local_temperature, sound_speed, local_ok)
      if (.not. local_ok) return
      output_temperature(cell) = local_temperature
    end do
    ok = .true.
  end subroutine plm_euler_update_1d

  subroutine build_plm_cell_edges_1d( &
      species, state, temperature, nx, limiter, left_state, right_state, &
      left_temperature, right_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    character(len=*), intent(in) :: limiter
    real(dp), intent(out) :: left_state(:, :), right_state(:, :)
    real(dp), intent(out) :: left_temperature(:), right_temperature(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:, :), slope(:)
    real(dp), allocatable :: left_primitive(:), right_primitive(:)
    real(dp) :: local_temperature, sound_speed
    logical :: local_ok, slope_ok
    integer :: nprimitive, cell, component

    ok = .false.
    nprimitive = reactive_nprim(size(species))
    allocate(primitive(nprimitive, 0:nx + 1), slope(nprimitive))
    allocate(left_primitive(nprimitive), right_primitive(nprimitive))
    do cell = 0, nx + 1
      call reactive_conserved_to_primitive( &
        species, state(:, cell), temperature(cell), primitive(:, cell), &
        local_temperature, sound_speed, local_ok)
      if (.not. local_ok) return
    end do
    do cell = 1, nx
      do component = 1, nprimitive
        call limited_slope( &
          primitive(component, cell) - primitive(component, cell - 1), &
          primitive(component, cell + 1) - primitive(component, cell), &
          limiter, slope(component), slope_ok)
        if (.not. slope_ok) return
      end do
      left_primitive = primitive(:, cell) - 0.5_dp * slope
      right_primitive = primitive(:, cell) + 0.5_dp * slope
      call sanitize_amr_primitive_1d( &
        left_primitive, primitive(:, cell), size(species))
      call sanitize_amr_primitive_1d( &
        right_primitive, primitive(:, cell), size(species))
      call reactive_primitive_to_conserved( &
        species, left_primitive, left_state(:, cell), &
        left_temperature(cell), sound_speed, local_ok)
      if (.not. local_ok) then
        call reactive_primitive_to_conserved( &
          species, primitive(:, cell), left_state(:, cell), &
          left_temperature(cell), sound_speed, local_ok)
        if (.not. local_ok) return
      end if
      call reactive_primitive_to_conserved( &
        species, right_primitive, right_state(:, cell), &
        right_temperature(cell), sound_speed, local_ok)
      if (.not. local_ok) then
        call reactive_primitive_to_conserved( &
          species, primitive(:, cell), right_state(:, cell), &
          right_temperature(cell), sound_speed, local_ok)
        if (.not. local_ok) return
      end if
    end do
    ok = .true.
  end subroutine build_plm_cell_edges_1d

  pure subroutine sanitize_amr_primitive_1d(primitive, fallback, nspecies)
    real(dp), intent(inout) :: primitive(:)
    real(dp), intent(in) :: fallback(:)
    integer, intent(in) :: nspecies

    real(dp) :: total
    integer :: component, k

    if (primitive(1) <= density_floor .or. &
        primitive(5) <= pressure_floor) then
      primitive = fallback
      return
    end if
    total = 0.0_dp
    do k = 1, nspecies
      component = reactive_mass_fraction_component(k)
      primitive(component) = max(0.0_dp, primitive(component))
      total = total + primitive(component)
    end do
    if (total <= tiny(1.0_dp)) then
      primitive = fallback
      return
    end if
    do k = 1, nspecies
      component = reactive_mass_fraction_component(k)
      primitive(component) = primitive(component) / total
    end do
  end subroutine sanitize_amr_primitive_1d

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

  pure subroutine fill_physical_wide_ghosts_1d( &
      state, temperature, nx, boundary, left_ghost, right_ghost, &
      left_temperature, right_temperature, ok)
    real(dp), intent(in) :: state(:, 0:), temperature(0:)
    integer, intent(in) :: nx
    character(len=*), intent(in) :: boundary
    real(dp), intent(out) :: left_ghost(:, :), right_ghost(:, :)
    real(dp), intent(out) :: left_temperature(:), right_temperature(:)
    logical, intent(out) :: ok

    integer :: layer, left_source, right_source

    left_ghost = 0.0_dp
    right_ghost = 0.0_dp
    left_temperature = 0.0_dp
    right_temperature = 0.0_dp
    ok = nx >= amr_ppm_ghost_width .and. &
      size(left_ghost, 1) == size(state, 1) .and. &
      size(right_ghost, 1) == size(state, 1) .and. &
      size(left_ghost, 2) == amr_ppm_ghost_width .and. &
      size(right_ghost, 2) == amr_ppm_ghost_width .and. &
      size(left_temperature) == amr_ppm_ghost_width .and. &
      size(right_temperature) == amr_ppm_ghost_width
    if (.not. ok) return
    do layer = 1, amr_ppm_ghost_width
      select case (trim(boundary))
      case ("periodic")
        left_source = modulo(nx - layer, nx) + 1
        right_source = modulo(layer - 1, nx) + 1
      case ("outflow")
        left_source = 1
        right_source = nx
      case default
        ok = .false.
        return
      end select
      left_ghost(:, layer) = state(:, left_source)
      right_ghost(:, layer) = state(:, right_source)
      left_temperature(layer) = temperature(left_source)
      right_temperature(layer) = temperature(right_source)
    end do
    ok = .true.
  end subroutine fill_physical_wide_ghosts_1d

  subroutine fill_fine_wide_ghosts_1d( &
      species, hierarchy, coarse_start, coarse_end, alpha, left_ghost, &
      right_ghost, left_temperature, right_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_two_level_hierarchy_1d), intent(in) :: hierarchy
    real(dp), intent(in) :: coarse_start(:, 0:), coarse_end(:, 0:)
    real(dp), intent(in) :: alpha
    real(dp), intent(out) :: left_ghost(:, :), right_ghost(:, :)
    real(dp), intent(out) :: left_temperature(:), right_temperature(:)
    logical, intent(out) :: ok

    logical :: local_ok
    integer :: layer, global_fine

    left_ghost = 0.0_dp
    right_ghost = 0.0_dp
    left_temperature = 0.0_dp
    right_temperature = 0.0_dp
    ok = hierarchy%is_valid() .and. alpha >= 0.0_dp .and. alpha <= 1.0_dp .and. &
      size(coarse_start, 1) == size(coarse_end, 1) .and. &
      size(left_ghost, 1) == size(coarse_start, 1) .and. &
      size(right_ghost, 1) == size(coarse_start, 1) .and. &
      size(left_ghost, 2) == amr_ppm_ghost_width .and. &
      size(right_ghost, 2) == amr_ppm_ghost_width .and. &
      size(left_temperature) == amr_ppm_ghost_width .and. &
      size(right_temperature) == amr_ppm_ghost_width
    if (.not. ok) return
    do layer = 1, amr_ppm_ghost_width
      ! Layer one is adjacent to the fine patch; subsequent layers move away.
      global_fine = hierarchy%fine%lower - layer
      call interpolate_parent_fine_cell( &
        species, hierarchy, coarse_start, coarse_end, alpha, global_fine, &
        left_ghost(:, layer), left_temperature(layer), local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      global_fine = hierarchy%fine%upper + layer
      call interpolate_parent_fine_cell( &
        species, hierarchy, coarse_start, coarse_end, alpha, global_fine, &
        right_ghost(:, layer), right_temperature(layer), local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
    end do
    ok = .true.
  end subroutine fill_fine_wide_ghosts_1d

  subroutine interpolate_parent_fine_cell( &
      species, hierarchy, coarse_start, coarse_end, alpha, global_fine, &
      sampled_state, sampled_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_two_level_hierarchy_1d), intent(in) :: hierarchy
    real(dp), intent(in) :: coarse_start(:, 0:), coarse_end(:, 0:)
    real(dp), intent(in) :: alpha
    integer, intent(in) :: global_fine
    real(dp), intent(out) :: sampled_state(:), sampled_temperature
    logical, intent(out) :: ok

    real(dp), allocatable :: center(:), left(:), right(:), slope(:), q(:)
    real(dp) :: offset, local_temperature, sound_speed
    logical :: local_ok, slope_ok
    integer :: parent_cell, child, component, ratio, parent_cells

    sampled_state = 0.0_dp
    sampled_temperature = 0.0_dp
    ok = .false.
    ratio = hierarchy%refinement_ratio
    parent_cells = hierarchy%coarse%cell_count()
    if (global_fine < 1 .or. global_fine > parent_cells * ratio .or. &
        size(sampled_state) /= size(coarse_start, 1)) return
    parent_cell = (global_fine - 1) / ratio + 1
    child = modulo(global_fine - 1, ratio) + 1
    allocate(center(size(sampled_state)), left(size(sampled_state)))
    allocate(right(size(sampled_state)), slope(size(sampled_state)))
    allocate(q(reactive_nprim(size(species))))
    ! Interpolate the parent conserved averages in time first, then form a
    ! conservative MC-limited subcell value at the requested refined center.
    center = (1.0_dp - alpha) * coarse_start(:, parent_cell) + &
      alpha * coarse_end(:, parent_cell)
    left = (1.0_dp - alpha) * coarse_start(:, parent_cell - 1) + &
      alpha * coarse_end(:, parent_cell - 1)
    right = (1.0_dp - alpha) * coarse_start(:, parent_cell + 1) + &
      alpha * coarse_end(:, parent_cell + 1)
    do component = 1, size(sampled_state)
      call limited_slope( &
        center(component) - left(component), &
        right(component) - center(component), "mc", slope(component), &
        slope_ok)
      if (.not. slope_ok) return
    end do
    offset = (real(child, dp) - 0.5_dp) / real(ratio, dp) - 0.5_dp
    sampled_state = center + offset * slope
    call reactive_conserved_to_primitive( &
      species, sampled_state, 1000.0_dp, q, local_temperature, sound_speed, &
      local_ok)
    if (.not. local_ok) then
      sampled_state = center
      call reactive_conserved_to_primitive( &
        species, sampled_state, 1000.0_dp, q, local_temperature, &
        sound_speed, local_ok)
      if (.not. local_ok) return
    end if
    sampled_temperature = local_temperature
    ok = .true.
  end subroutine interpolate_parent_fine_cell

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
