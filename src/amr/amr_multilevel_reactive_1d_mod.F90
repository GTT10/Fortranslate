module amr_multilevel_reactive_1d_mod
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: gas_transport_species
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_cfl_timestep, &
    reactive_transport_timestep, &
    initialize_reactive_1d, advance_reactive_chemistry
  use amr_hierarchy_1d_mod, only: &
    amr_two_level_hierarchy_1d, amr_multilevel_hierarchy_1d, &
    amr_level_field_1d, &
    amr_flux_register_1d, initialize_multilevel_hierarchy_1d, &
    initialize_two_level_hierarchy_1d, prolong_conservative_1d, &
    prolong_multilevel_1d, initialize_flux_register_1d, &
    accumulate_coarse_flux_1d, accumulate_fine_flux_1d, reflux_1d, &
    average_down_1d, composite_integral_multilevel_1d
  use amr_reactive_1d_mod, only: &
    amr_ppm_ghost_width, &
    advance_amr_level_1d, advance_transport_level_1d, &
    recover_level_temperatures_1d, fill_physical_ghosts_1d, &
    fill_fine_ghosts_1d, fill_physical_wide_ghosts_1d, &
    fill_fine_wide_ghosts_1d, write_amr_cell
  use amr_regrid_1d_mod, only: &
    amr_tagging_criteria_1d, amr_regrid_plan_1d, &
    tag_gradient_1d, build_regrid_plan_1d
  implicit none
  private

  type, public :: amr_reactive_level_1d
    real(dp), allocatable :: state(:, :)
    real(dp), allocatable :: temperature(:)
    real(dp), allocatable :: left_ghost_state(:, :)
    real(dp), allocatable :: right_ghost_state(:, :)
    real(dp), allocatable :: left_ghost_temperature(:)
    real(dp), allocatable :: right_ghost_temperature(:)
  end type amr_reactive_level_1d

  type, public :: amr_multilevel_reactive_solution_1d
    type(amr_multilevel_hierarchy_1d) :: hierarchy
    type(amr_reactive_level_1d), allocatable :: levels(:)
    real(dp) :: time = 0.0_dp
    integer :: steps = 0
    integer :: regrid_evaluations = 0
    integer :: regrids = 0
    integer :: overlap_cells_transferred = 0
  contains
    procedure :: level_count => multilevel_reactive_level_count
    procedure :: is_valid => multilevel_reactive_is_valid
  end type amr_multilevel_reactive_solution_1d

  public :: initialize_multilevel_reactive_1d
  public :: initialize_tagged_multilevel_reactive_1d
  public :: multilevel_reactive_timestep_1d
  public :: advance_multilevel_reactive_1d
  public :: regrid_multilevel_reactive_1d
  public :: simulate_multilevel_reactive_1d
  public :: multilevel_reactive_integrals_1d
  public :: write_multilevel_reactive_1d_csv

contains

  pure integer function multilevel_reactive_level_count(self) result(count)
    class(amr_multilevel_reactive_solution_1d), intent(in) :: self

    count = 0
    if (allocated(self%levels)) count = size(self%levels)
  end function multilevel_reactive_level_count

  pure logical function multilevel_reactive_is_valid(self) result(valid)
    class(amr_multilevel_reactive_solution_1d), intent(in) :: self

    integer :: level, nx, nvar

    valid = self%hierarchy%is_valid() .and. allocated(self%levels)
    if (.not. valid) return
    valid = size(self%levels) == self%hierarchy%level_count()
    if (.not. valid .or. size(self%levels) < 1) return
    valid = allocated(self%levels(1)%state) .and. &
      allocated(self%levels(1)%temperature)
    if (.not. valid) return
    nvar = size(self%levels(1)%state, 1)
    valid = nvar >= 1
    if (.not. valid) return
    do level = 1, size(self%levels)
      valid = allocated(self%levels(level)%state) .and. &
        allocated(self%levels(level)%temperature) .and. &
        allocated(self%levels(level)%left_ghost_state) .and. &
        allocated(self%levels(level)%right_ghost_state) .and. &
        allocated(self%levels(level)%left_ghost_temperature) .and. &
        allocated(self%levels(level)%right_ghost_temperature)
      if (.not. valid) return
      nx = self%hierarchy%level_cell_count(level - 1)
      valid = lbound(self%levels(level)%state, 2) == 0 .and. &
        ubound(self%levels(level)%state, 2) == nx + 1 .and. &
        size(self%levels(level)%state, 1) == nvar .and. &
        lbound(self%levels(level)%temperature, 1) == 0 .and. &
        ubound(self%levels(level)%temperature, 1) == nx + 1 .and. &
        size(self%levels(level)%left_ghost_state, 1) == nvar .and. &
        size(self%levels(level)%right_ghost_state, 1) == nvar .and. &
        size(self%levels(level)%left_ghost_state, 2) == &
          amr_ppm_ghost_width .and. &
        size(self%levels(level)%right_ghost_state, 2) == &
          amr_ppm_ghost_width .and. &
        size(self%levels(level)%left_ghost_temperature) == &
          amr_ppm_ghost_width .and. &
        size(self%levels(level)%right_ghost_temperature) == &
          amr_ppm_ghost_width
      if (.not. valid) return
    end do
  end function multilevel_reactive_is_valid

  subroutine initialize_multilevel_reactive_1d( &
      species, config, patch_parent_lower, patch_parent_upper, &
      refinement_ratios, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    integer, intent(in) :: patch_parent_lower(:), patch_parent_upper(:)
    integer, intent(in) :: refinement_ratios(:)
    type(amr_multilevel_reactive_solution_1d), intent(out) :: solution
    logical, intent(out) :: ok

    type(amr_level_field_1d), allocatable :: fields(:)
    real(dp), allocatable :: root_state(:, :), root_temperature(:)
    real(dp) :: root_dx
    logical :: local_ok
    integer :: level, nx, nvar

    ok = .false.
    if (size(species) < 1) return
    call initialize_reactive_1d( &
      species, config, root_state, root_temperature, root_dx, local_ok)
    if (.not. local_ok) return
    call initialize_multilevel_hierarchy_1d( &
      config%nx, patch_parent_lower, patch_parent_upper, refinement_ratios, &
      config%x_lower, config%x_upper, solution%hierarchy, local_ok)
    if (.not. local_ok) return
    call prolong_multilevel_1d( &
      root_state(:, 1:config%nx), solution%hierarchy, fields, local_ok)
    if (.not. local_ok) return

    nvar = reactive_nvar(size(species))
    allocate(solution%levels(solution%hierarchy%level_count()))
    do level = 1, size(solution%levels)
      nx = solution%hierarchy%level_cell_count(level - 1)
      allocate(solution%levels(level)%state(nvar, 0:nx + 1))
      allocate(solution%levels(level)%temperature(0:nx + 1))
      allocate(solution%levels(level)%left_ghost_state( &
        nvar, amr_ppm_ghost_width))
      allocate(solution%levels(level)%right_ghost_state( &
        nvar, amr_ppm_ghost_width))
      allocate(solution%levels(level)%left_ghost_temperature( &
        amr_ppm_ghost_width))
      allocate(solution%levels(level)%right_ghost_temperature( &
        amr_ppm_ghost_width))
      solution%levels(level)%state = 0.0_dp
      solution%levels(level)%temperature = 0.0_dp
      solution%levels(level)%left_ghost_state = 0.0_dp
      solution%levels(level)%right_ghost_state = 0.0_dp
      solution%levels(level)%left_ghost_temperature = 0.0_dp
      solution%levels(level)%right_ghost_temperature = 0.0_dp
      solution%levels(level)%state(:, 1:nx) = fields(level)%values
      if (level == 1) then
        solution%levels(level)%temperature = root_temperature
      else
        call recover_level_temperatures_1d( &
          species, solution%levels(level)%state, &
          solution%levels(level)%temperature, nx, local_ok)
        if (.not. local_ok) return
      end if
    end do
    call refresh_multilevel_ghosts( &
      species, config, solution, local_ok)
    if (.not. local_ok) return
    solution%time = 0.0_dp
    solution%steps = 0
    solution%regrid_evaluations = 0
    solution%regrids = 0
    solution%overlap_cells_transferred = 0
    ok = solution%is_valid() .and. root_dx > 0.0_dp
  end subroutine initialize_multilevel_reactive_1d

  subroutine initialize_tagged_multilevel_reactive_1d( &
      species, config, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_multilevel_reactive_solution_1d), intent(out) :: solution
    logical, intent(out) :: ok

    real(dp), allocatable :: root_state(:, :), root_temperature(:)
    real(dp) :: root_dx
    logical :: local_ok

    ok = .false.
    if (.not. config%amr_enabled .or. config%amr_max_levels < 2) return
    call initialize_reactive_1d( &
      species, config, root_state, root_temperature, root_dx, local_ok)
    if (.not. local_ok .or. root_dx <= 0.0_dp) return
    call build_tagged_solution_from_root( &
      species, config, root_state, root_temperature, solution, local_ok)
    if (.not. local_ok) return
    solution%time = 0.0_dp
    solution%steps = 0
    solution%regrid_evaluations = 1
    solution%regrids = merge(1, 0, solution%level_count() > 1)
    solution%overlap_cells_transferred = 0
    ok = .true.
  end subroutine initialize_tagged_multilevel_reactive_1d

  subroutine build_tagged_solution_from_root( &
      species, config, root_state, root_temperature, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: root_state(:, 0:), root_temperature(0:)
    type(amr_multilevel_reactive_solution_1d), intent(out) :: solution
    logical, intent(out) :: ok

    type(amr_level_field_1d), allocatable :: candidate(:)
    type(amr_two_level_hierarchy_1d) :: relation_geometry
    type(amr_tagging_criteria_1d) :: criteria
    type(amr_regrid_plan_1d) :: plan
    integer, allocatable :: patch_lower(:), patch_upper(:), ratios(:)
    logical, allocatable :: tags(:)
    real(dp) :: parent_lower, parent_upper, child_lower, child_upper
    logical :: local_ok
    integer :: level, relation_count, nx, nvar, maximum_relations
    integer :: parent_buffer, allowed_lower, allowed_upper

    ok = .false.
    nvar = reactive_nvar(size(species))
    if (.not. config%amr_enabled .or. config%amr_max_levels < 2 .or. &
        size(root_state, 1) /= nvar .or. &
        ubound(root_state, 2) /= config%nx + 1 .or. &
        ubound(root_temperature, 1) /= config%nx + 1) return
    maximum_relations = config%amr_max_levels - 1
    allocate(candidate(config%amr_max_levels))
    allocate(patch_lower(maximum_relations), patch_upper(maximum_relations))
    allocate(ratios(maximum_relations))
    allocate(candidate(1)%values(nvar, config%nx))
    candidate(1)%values = root_state(:, 1:config%nx)
    parent_lower = config%x_lower
    parent_upper = config%x_upper
    relation_count = 0
    call criteria_from_config(config, criteria)

    do level = 1, maximum_relations
      nx = size(candidate(level)%values, 2)
      criteria%minimum_patch_cells = &
        min(config%amr_minimum_patch_cells, nx - 2)
      if (criteria%minimum_patch_cells < 1) exit
      allocate(tags(nx))
      call tag_gradient_1d( &
        candidate(level)%values, criteria, tags, local_ok)
      if (.not. local_ok) return
      if (uses_ppm_reconstruction(config)) then
        parent_buffer = (amr_ppm_ghost_width + &
          config%amr_refinement_ratio - 1) / &
          config%amr_refinement_ratio + 1
        if (2 * parent_buffer >= nx) then
          tags = .false.
        else
          tags(1:parent_buffer) = .false.
          tags(nx - parent_buffer + 1:nx) = .false.
        end if
      else if (level > 1) then
        tags(1) = .false.
        tags(nx) = .false.
      end if
      call build_regrid_plan_1d( &
        tags, criteria%buffer_cells, criteria%minimum_patch_cells, &
        plan, local_ok)
      deallocate(tags)
      if (.not. local_ok) return
      if (.not. plan%active) exit
      if (uses_ppm_reconstruction(config)) then
        allowed_lower = parent_buffer + 1
        allowed_upper = nx - parent_buffer
        plan%patch_lower = max(plan%patch_lower, allowed_lower)
        plan%patch_upper = min(plan%patch_upper, allowed_upper)
        do while (plan%patch_upper - plan%patch_lower + 1 < &
            criteria%minimum_patch_cells)
          if (plan%patch_lower > allowed_lower) &
            plan%patch_lower = plan%patch_lower - 1
          if (plan%patch_upper - plan%patch_lower + 1 >= &
              criteria%minimum_patch_cells) exit
          if (plan%patch_upper < allowed_upper) &
            plan%patch_upper = plan%patch_upper + 1
          if (plan%patch_lower == allowed_lower .and. &
              plan%patch_upper == allowed_upper) exit
        end do
        if (plan%patch_upper - plan%patch_lower + 1 < &
            criteria%minimum_patch_cells) exit
      end if
      relation_count = relation_count + 1
      patch_lower(relation_count) = plan%patch_lower
      patch_upper(relation_count) = plan%patch_upper
      ratios(relation_count) = config%amr_refinement_ratio
      call initialize_two_level_hierarchy_1d( &
        nx, plan%patch_lower, plan%patch_upper, &
        config%amr_refinement_ratio, parent_lower, parent_upper, &
        relation_geometry, local_ok, level - 1)
      if (.not. local_ok) return
      allocate(candidate(level + 1)%values( &
        nvar, relation_geometry%fine%cell_count()))
      call prolong_conservative_1d( &
        candidate(level)%values, relation_geometry, &
        candidate(level + 1)%values, local_ok)
      if (.not. local_ok) return
      child_lower = parent_lower + &
        real(plan%patch_lower - 1, dp) * relation_geometry%coarse_dx
      child_upper = parent_lower + &
        real(plan%patch_upper, dp) * relation_geometry%coarse_dx
      parent_lower = child_lower
      parent_upper = child_upper
    end do

    call initialize_multilevel_hierarchy_1d( &
      config%nx, patch_lower(1:relation_count), &
      patch_upper(1:relation_count), ratios(1:relation_count), &
      config%x_lower, config%x_upper, solution%hierarchy, local_ok)
    if (.not. local_ok) return
    allocate(solution%levels(relation_count + 1))
    do level = 1, relation_count + 1
      nx = solution%hierarchy%level_cell_count(level - 1)
      allocate(solution%levels(level)%state(nvar, 0:nx + 1))
      allocate(solution%levels(level)%temperature(0:nx + 1))
      allocate(solution%levels(level)%left_ghost_state( &
        nvar, amr_ppm_ghost_width))
      allocate(solution%levels(level)%right_ghost_state( &
        nvar, amr_ppm_ghost_width))
      allocate(solution%levels(level)%left_ghost_temperature( &
        amr_ppm_ghost_width))
      allocate(solution%levels(level)%right_ghost_temperature( &
        amr_ppm_ghost_width))
      solution%levels(level)%state = 0.0_dp
      solution%levels(level)%temperature = 0.0_dp
      solution%levels(level)%left_ghost_state = 0.0_dp
      solution%levels(level)%right_ghost_state = 0.0_dp
      solution%levels(level)%left_ghost_temperature = 0.0_dp
      solution%levels(level)%right_ghost_temperature = 0.0_dp
      solution%levels(level)%state(:, 1:nx) = candidate(level)%values
      if (level == 1) then
        solution%levels(level)%state(:, 0) = root_state(:, 0)
        solution%levels(level)%state(:, nx + 1) = root_state(:, nx + 1)
        solution%levels(level)%temperature = root_temperature
      else
        call recover_level_temperatures_1d( &
          species, solution%levels(level)%state, &
          solution%levels(level)%temperature, nx, local_ok)
        if (.not. local_ok) return
      end if
    end do
    call refresh_multilevel_ghosts(species, config, solution, local_ok)
    if (.not. local_ok) return
    ok = solution%is_valid()
  end subroutine build_tagged_solution_from_root

  subroutine multilevel_reactive_timestep_1d( &
      species, config, solution, dt, ok, transport)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_multilevel_reactive_solution_1d), intent(in) :: solution
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok
    type(gas_transport_species), intent(in), optional :: transport(:)

    real(dp) :: local_dt, transport_dt, maximum_diffusivity
    real(dp) :: hydro_scale, transport_scale, dx
    logical :: local_ok
    integer :: level, nx

    dt = huge(1.0_dp)
    ok = solution%is_valid() .and. size(species) >= 1
    if (.not. ok) return
    if (config%transport_enabled .and. .not. present(transport)) then
      ok = .false.
      return
    end if
    hydro_scale = 1.0_dp
    do level = 1, solution%level_count()
      if (level > 1) hydro_scale = hydro_scale * real( &
        solution%hierarchy%interfaces(level - 1)%refinement_ratio, dp)
      nx = solution%hierarchy%level_cell_count(level - 1)
      dx = solution%hierarchy%level_dx(level - 1)
      call reactive_cfl_timestep( &
        species, solution%levels(level)%state, &
        solution%levels(level)%temperature, nx, dx, config%cfl, &
        local_dt, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      dt = min(dt, hydro_scale * local_dt)
      if (config%transport_enabled) then
        call reactive_transport_timestep( &
          species, transport, solution%levels(level)%state, &
          solution%levels(level)%temperature, nx, dx, &
          config%transport_cfl, config%viscosity_enabled, &
          config%thermal_conduction_enabled, &
          config%species_diffusion_enabled, transport_dt, &
          maximum_diffusivity, local_ok)
        if (.not. local_ok) then
          ok = .false.
          return
        end if
        transport_scale = hydro_scale * hydro_scale
        dt = min(dt, transport_scale * transport_dt)
      end if
    end do
    ok = dt > 0.0_dp .and. dt < huge(1.0_dp)
  end subroutine multilevel_reactive_timestep_1d

  subroutine advance_multilevel_reactive_1d( &
      species, reactions, config, dt, solution, ok, transport)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: dt
    type(amr_multilevel_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok
    type(gas_transport_species), intent(in), optional :: transport(:)

    type(amr_multilevel_reactive_solution_1d) :: backup
    real(dp), allocatable :: left_integral(:), right_integral(:)
    logical :: local_ok
    integer :: nvar

    ok = .false.
    if (dt <= 0.0_dp .or. .not. solution%is_valid()) return
    if (config%transport_enabled .and. .not. present(transport)) return
    backup = solution
    nvar = size(solution%levels(1)%state, 1)
    allocate(left_integral(nvar), right_integral(nvar))

    if (config%chemistry_enabled) then
      call advance_chemistry_all_levels( &
        species, reactions, config, 0.5_dp * dt, solution, local_ok)
      if (.not. local_ok) then
        solution = backup
        return
      end if
    end if
    if (config%transport_enabled) then
      call advance_transport_recursive( &
        species, transport, config, solution, 1, 0.5_dp * dt, &
        left_integral, right_integral, local_ok)
      if (.not. local_ok) then
        solution = backup
        return
      end if
      call refresh_multilevel_ghosts(species, config, solution, local_ok)
      if (.not. local_ok) then
        solution = backup
        return
      end if
    end if
    call advance_hydro_recursive( &
      species, config, solution, 1, dt, left_integral, right_integral, &
      local_ok)
    if (.not. local_ok) then
      solution = backup
      return
    end if
    call refresh_multilevel_ghosts(species, config, solution, local_ok)
    if (.not. local_ok) then
      solution = backup
      return
    end if
    if (config%transport_enabled) then
      call advance_transport_recursive( &
        species, transport, config, solution, 1, 0.5_dp * dt, &
        left_integral, right_integral, local_ok)
      if (.not. local_ok) then
        solution = backup
        return
      end if
      call refresh_multilevel_ghosts(species, config, solution, local_ok)
      if (.not. local_ok) then
        solution = backup
        return
      end if
    end if
    if (config%chemistry_enabled) then
      call advance_chemistry_all_levels( &
        species, reactions, config, 0.5_dp * dt, solution, local_ok)
      if (.not. local_ok) then
        solution = backup
        return
      end if
    end if
    call refresh_multilevel_ghosts(species, config, solution, local_ok)
    if (.not. local_ok) then
      solution = backup
      return
    end if
    solution%time = solution%time + dt
    solution%steps = solution%steps + 1
    ok = .true.
  end subroutine advance_multilevel_reactive_1d

  subroutine regrid_multilevel_reactive_1d( &
      species, config, solution, changed, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_multilevel_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: changed, ok

    type(amr_multilevel_reactive_solution_1d) :: backup, rebuilt
    real(dp), allocatable :: root_state(:, :), root_temperature(:)
    logical :: local_ok
    integer :: transferred_cells

    changed = .false.
    ok = .false.
    if (.not. solution%is_valid() .or. config%amr_max_levels < 2) return
    backup = solution
    call average_down_all_levels(species, config, solution, local_ok)
    if (.not. local_ok) then
      solution = backup
      return
    end if
    root_state = solution%levels(1)%state
    root_temperature = solution%levels(1)%temperature
    call build_tagged_solution_from_root( &
      species, config, root_state, root_temperature, rebuilt, local_ok)
    if (.not. local_ok) then
      solution = backup
      return
    end if
    changed = .not. same_hierarchy(backup%hierarchy, rebuilt%hierarchy)
    if (.not. changed) then
      solution = backup
      solution%regrid_evaluations = backup%regrid_evaluations + 1
      ok = .true.
      return
    end if
    rebuilt%time = backup%time
    rebuilt%steps = backup%steps
    rebuilt%regrid_evaluations = backup%regrid_evaluations + 1
    rebuilt%regrids = backup%regrids + 1
    call transfer_multilevel_overlap( &
      backup, rebuilt, transferred_cells, local_ok)
    if (.not. local_ok) then
      solution = backup
      return
    end if
    call average_down_all_levels(species, config, rebuilt, local_ok)
    if (.not. local_ok) then
      solution = backup
      return
    end if
    rebuilt%overlap_cells_transferred = &
      backup%overlap_cells_transferred + transferred_cells
    solution = rebuilt
    ok = .true.
  end subroutine regrid_multilevel_reactive_1d

  subroutine transfer_multilevel_overlap( &
      old_solution, new_solution, transferred_cells, ok)
    type(amr_multilevel_reactive_solution_1d), intent(in) :: old_solution
    type(amr_multilevel_reactive_solution_1d), intent(inout) :: new_solution
    integer, intent(out) :: transferred_cells
    logical, intent(out) :: ok

    real(dp) :: old_lower, old_upper, new_lower, new_upper
    real(dp) :: overlap_lower, overlap_upper, old_dx, new_dx, tolerance
    real(dp) :: old_offset, new_offset
    logical :: old_bounds_ok, new_bounds_ok
    integer :: level, common_levels, old_first, new_first, cell_count

    transferred_cells = 0
    ok = old_solution%is_valid() .and. new_solution%is_valid()
    if (.not. ok) return
    common_levels = min( &
      old_solution%level_count(), new_solution%level_count())
    do level = 2, common_levels
      old_dx = old_solution%hierarchy%level_dx(level - 1)
      new_dx = new_solution%hierarchy%level_dx(level - 1)
      tolerance = 128.0_dp * epsilon(1.0_dp) * &
        max(1.0_dp, abs(old_dx), abs(new_dx))
      if (abs(old_dx - new_dx) > tolerance) cycle
      call old_solution%hierarchy%level_bounds( &
        level - 1, old_lower, old_upper, old_bounds_ok)
      call new_solution%hierarchy%level_bounds( &
        level - 1, new_lower, new_upper, new_bounds_ok)
      if (.not. old_bounds_ok .or. .not. new_bounds_ok) then
        ok = .false.
        return
      end if
      overlap_lower = max(old_lower, new_lower)
      overlap_upper = min(old_upper, new_upper)
      if (overlap_upper <= overlap_lower + tolerance) cycle
      old_offset = (overlap_lower - old_lower) / old_dx
      new_offset = (overlap_lower - new_lower) / new_dx
      if (abs(old_offset - real(nint(old_offset), dp)) > tolerance / old_dx .or. &
          abs(new_offset - real(nint(new_offset), dp)) > &
            tolerance / new_dx) then
        ok = .false.
        return
      end if
      old_first = nint(old_offset) + 1
      new_first = nint(new_offset) + 1
      cell_count = nint((overlap_upper - overlap_lower) / old_dx)
      if (cell_count < 1) cycle
      new_solution%levels(level)%state(:, &
        new_first:new_first + cell_count - 1) = &
        old_solution%levels(level)%state(:, &
          old_first:old_first + cell_count - 1)
      new_solution%levels(level)%temperature( &
        new_first:new_first + cell_count - 1) = &
        old_solution%levels(level)%temperature( &
          old_first:old_first + cell_count - 1)
      transferred_cells = transferred_cells + cell_count
    end do
    ok = .true.
  end subroutine transfer_multilevel_overlap

  subroutine simulate_multilevel_reactive_1d( &
      species, reactions, config, solution, initial_integrals, &
      final_integrals, ok, transport)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_multilevel_reactive_solution_1d), intent(out) :: solution
    real(dp), intent(out) :: initial_integrals(5), final_integrals(5)
    logical, intent(out) :: ok
    type(gas_transport_species), intent(in), optional :: transport(:)

    real(dp), allocatable :: all_integrals(:)
    real(dp) :: dt, tolerance
    logical :: local_ok, changed
    integer :: nvar

    initial_integrals = 0.0_dp
    final_integrals = 0.0_dp
    ok = .false.
    if (config%transport_enabled .and. .not. present(transport)) return
    call initialize_tagged_multilevel_reactive_1d( &
      species, config, solution, local_ok)
    if (.not. local_ok) return
    nvar = reactive_nvar(size(species))
    allocate(all_integrals(nvar))
    call multilevel_reactive_integrals_1d( &
      solution, all_integrals, local_ok)
    if (.not. local_ok) return
    initial_integrals = all_integrals([irho, imx, imy, imz, iet])
    tolerance = 50.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, config%final_time)
    do while (solution%time < config%final_time - tolerance)
      if (solution%steps >= config%maximum_steps) return
      if (config%transport_enabled) then
        call multilevel_reactive_timestep_1d( &
          species, config, solution, dt, local_ok, transport)
      else
        call multilevel_reactive_timestep_1d( &
          species, config, solution, dt, local_ok)
      end if
      if (.not. local_ok) return
      dt = min(dt, config%final_time - solution%time)
      if (config%transport_enabled) then
        call advance_multilevel_reactive_1d( &
          species, reactions, config, dt, solution, local_ok, transport)
      else
        call advance_multilevel_reactive_1d( &
          species, reactions, config, dt, solution, local_ok)
      end if
      if (.not. local_ok) return
      if (mod(solution%steps, config%amr_regrid_interval) == 0) then
        call regrid_multilevel_reactive_1d( &
          species, config, solution, changed, local_ok)
        if (.not. local_ok) return
      end if
    end do
    solution%time = config%final_time
    call multilevel_reactive_integrals_1d(solution, all_integrals, local_ok)
    if (.not. local_ok) return
    final_integrals = all_integrals([irho, imx, imy, imz, iet])
    ok = .true.
  end subroutine simulate_multilevel_reactive_1d

  recursive subroutine advance_hydro_recursive( &
      species, config, solution, level, interval, left_integral, &
      right_integral, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_multilevel_reactive_solution_1d), intent(inout) :: solution
    integer, intent(in) :: level
    real(dp), intent(in) :: interval
    real(dp), intent(out) :: left_integral(:), right_integral(:)
    logical, intent(out) :: ok

    type(amr_flux_register_1d) :: flux_register
    real(dp), allocatable :: state_start(:, :), state_end(:, :), flux(:, :)
    real(dp), allocatable :: child_left(:), child_right(:)
    real(dp) :: child_interval, alpha, dx
    logical :: local_ok, physical_boundary
    integer :: nx, nvar, ratio, substep

    ok = .false.
    left_integral = 0.0_dp
    right_integral = 0.0_dp
    if (interval <= 0.0_dp) return
    nx = solution%hierarchy%level_cell_count(level - 1)
    dx = solution%hierarchy%level_dx(level - 1)
    nvar = size(solution%levels(level)%state, 1)
    allocate(state_start(nvar, 0:nx + 1), state_end(nvar, 0:nx + 1))
    allocate(flux(nvar, 0:nx))
    state_start = solution%levels(level)%state
    physical_boundary = level == 1
    if (physical_boundary .or. .not. uses_ppm_reconstruction(config)) then
      call advance_amr_level_1d( &
        species, solution%levels(level)%state, &
        solution%levels(level)%temperature, nx, dx, interval, &
        config%amr_reconstruction, config%limiter, config%riemann_solver, &
        physical_boundary, level_boundary(config, level), flux, local_ok, &
        ppm_contact_steepening=config%ppm_contact_steepening, &
        ppm_shock_flattening=config%ppm_shock_flattening)
    else
      call advance_amr_level_1d( &
        species, solution%levels(level)%state, &
        solution%levels(level)%temperature, nx, dx, interval, &
        config%amr_reconstruction, config%limiter, config%riemann_solver, &
        physical_boundary, level_boundary(config, level), flux, local_ok, &
        solution%levels(level)%left_ghost_state, &
        solution%levels(level)%right_ghost_state, &
        solution%levels(level)%left_ghost_temperature, &
        solution%levels(level)%right_ghost_temperature, &
        config%ppm_contact_steepening, config%ppm_shock_flattening)
    end if
    if (.not. local_ok) return
    state_end = solution%levels(level)%state
    left_integral = interval * flux(:, 0)
    right_integral = interval * flux(:, nx)
    if (level == solution%level_count()) then
      ok = .true.
      return
    end if

    ratio = solution%hierarchy%interfaces(level)%refinement_ratio
    child_interval = interval / real(ratio, dp)
    allocate(child_left(nvar), child_right(nvar))
    call initialize_flux_register_1d(flux_register, nvar, local_ok)
    if (.not. local_ok) return
    call accumulate_coarse_flux_1d( &
      flux_register, &
      flux(:, solution%hierarchy%interfaces(level)%fine_coarse_lower - 1), &
      flux(:, solution%hierarchy%interfaces(level)%fine_coarse_upper), &
      interval, local_ok)
    if (.not. local_ok) return
    do substep = 1, ratio
      if (trim(config%amr_reconstruction) /= "pcm") then
        alpha = (real(substep, dp) - 0.5_dp) / real(ratio, dp)
      else
        alpha = real(substep - 1, dp) / real(ratio, dp)
      end if
      call fill_fine_ghosts_1d( &
        species, solution%hierarchy%interfaces(level), state_start, &
        state_end, alpha, solution%levels(level + 1)%state, &
        solution%levels(level + 1)%temperature, local_ok)
      if (.not. local_ok) return
      if (uses_ppm_reconstruction(config)) then
        call fill_fine_wide_ghosts_1d( &
          species, solution%hierarchy%interfaces(level), state_start, &
          state_end, alpha, &
          solution%levels(level + 1)%left_ghost_state, &
          solution%levels(level + 1)%right_ghost_state, &
          solution%levels(level + 1)%left_ghost_temperature, &
          solution%levels(level + 1)%right_ghost_temperature, local_ok)
        if (.not. local_ok) return
      end if
      call advance_hydro_recursive( &
        species, config, solution, level + 1, child_interval, &
        child_left, child_right, local_ok)
      if (.not. local_ok) return
      call accumulate_fine_flux_1d( &
        flux_register, child_left / child_interval, &
        child_right / child_interval, child_interval, local_ok)
      if (.not. local_ok) return
    end do
    call synchronize_relation( &
      species, solution, level, flux_register, local_ok)
    if (.not. local_ok) return
    ok = .true.
  end subroutine advance_hydro_recursive

  recursive subroutine advance_transport_recursive( &
      species, transport, config, solution, level, interval, &
      left_integral, right_integral, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_multilevel_reactive_solution_1d), intent(inout) :: solution
    integer, intent(in) :: level
    real(dp), intent(in) :: interval
    real(dp), intent(out) :: left_integral(:), right_integral(:)
    logical, intent(out) :: ok

    type(amr_flux_register_1d) :: flux_register
    real(dp), allocatable :: state_start(:, :), state_end(:, :), flux(:, :)
    real(dp), allocatable :: child_left(:), child_right(:)
    real(dp) :: child_interval, alpha, dx, boundary_distance
    logical :: local_ok, physical_boundary
    integer :: nx, nvar, ratio, subcycles, substep

    ok = .false.
    left_integral = 0.0_dp
    right_integral = 0.0_dp
    if (interval <= 0.0_dp) return
    nx = solution%hierarchy%level_cell_count(level - 1)
    dx = solution%hierarchy%level_dx(level - 1)
    nvar = size(solution%levels(level)%state, 1)
    allocate(state_start(nvar, 0:nx + 1), state_end(nvar, 0:nx + 1))
    allocate(flux(nvar, 0:nx))
    state_start = solution%levels(level)%state
    physical_boundary = level == 1
    boundary_distance = dx
    if (.not. physical_boundary) boundary_distance = 0.5_dp * ( &
      solution%hierarchy%interfaces(level - 1)%coarse_dx + dx)
    call advance_transport_level_1d( &
      species, transport, solution%levels(level)%state, &
      solution%levels(level)%temperature, nx, dx, interval, &
      boundary_distance, config, physical_boundary, &
      level_boundary(config, level), flux, local_ok)
    if (.not. local_ok) return
    state_end = solution%levels(level)%state
    left_integral = interval * flux(:, 0)
    right_integral = interval * flux(:, nx)
    if (level == solution%level_count()) then
      ok = .true.
      return
    end if

    ratio = solution%hierarchy%interfaces(level)%refinement_ratio
    subcycles = ratio * ratio
    child_interval = interval / real(subcycles, dp)
    allocate(child_left(nvar), child_right(nvar))
    call initialize_flux_register_1d(flux_register, nvar, local_ok)
    if (.not. local_ok) return
    call accumulate_coarse_flux_1d( &
      flux_register, &
      flux(:, solution%hierarchy%interfaces(level)%fine_coarse_lower - 1), &
      flux(:, solution%hierarchy%interfaces(level)%fine_coarse_upper), &
      interval, local_ok)
    if (.not. local_ok) return
    do substep = 1, subcycles
      alpha = (real(substep, dp) - 0.5_dp) / real(subcycles, dp)
      call fill_fine_ghosts_1d( &
        species, solution%hierarchy%interfaces(level), state_start, &
        state_end, alpha, solution%levels(level + 1)%state, &
        solution%levels(level + 1)%temperature, local_ok)
      if (.not. local_ok) return
      call advance_transport_recursive( &
        species, transport, config, solution, level + 1, child_interval, &
        child_left, child_right, local_ok)
      if (.not. local_ok) return
      call accumulate_fine_flux_1d( &
        flux_register, child_left / child_interval, &
        child_right / child_interval, child_interval, local_ok)
      if (.not. local_ok) return
    end do
    call synchronize_relation( &
      species, solution, level, flux_register, local_ok)
    if (.not. local_ok) return
    ok = .true.
  end subroutine advance_transport_recursive

  subroutine advance_chemistry_all_levels( &
      species, reactions, config, interval, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: interval
    type(amr_multilevel_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok

    logical :: local_ok
    integer :: level, nx

    ok = .false.
    do level = 1, solution%level_count()
      nx = solution%hierarchy%level_cell_count(level - 1)
      call advance_reactive_chemistry( &
        species, reactions, solution%levels(level)%state, &
        solution%levels(level)%temperature, nx, interval, &
        config%chemistry_relative_tolerance, &
        config%chemistry_absolute_tolerance, &
        chemistry_boundary(config, level), local_ok)
      if (.not. local_ok) return
    end do
    call average_down_all_levels(species, config, solution, local_ok)
    if (.not. local_ok) return
    ok = .true.
  end subroutine advance_chemistry_all_levels

  subroutine synchronize_relation(species, solution, relation, register, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_multilevel_reactive_solution_1d), intent(inout) :: solution
    integer, intent(in) :: relation
    type(amr_flux_register_1d), intent(inout) :: register
    logical, intent(out) :: ok

    logical :: local_ok
    integer :: parent_cells, child_cells

    ok = .false.
    parent_cells = solution%hierarchy%level_cell_count(relation - 1)
    child_cells = solution%hierarchy%level_cell_count(relation)
    call reflux_1d( &
      solution%levels(relation)%state(:, 1:parent_cells), &
      solution%hierarchy%interfaces(relation), register, local_ok)
    if (.not. local_ok) return
    call average_down_1d( &
      solution%levels(relation + 1)%state(:, 1:child_cells), &
      solution%hierarchy%interfaces(relation), &
      solution%levels(relation)%state(:, 1:parent_cells), local_ok)
    if (.not. local_ok) return
    call recover_level_temperatures_1d( &
      species, solution%levels(relation)%state, &
      solution%levels(relation)%temperature, parent_cells, local_ok)
    if (.not. local_ok) return
    ok = .true.
  end subroutine synchronize_relation

  subroutine average_down_all_levels(species, config, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_multilevel_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok

    logical :: local_ok
    integer :: relation, parent_cells, child_cells

    ok = .false.
    do relation = size(solution%hierarchy%interfaces), 1, -1
      parent_cells = solution%hierarchy%level_cell_count(relation - 1)
      child_cells = solution%hierarchy%level_cell_count(relation)
      call average_down_1d( &
        solution%levels(relation + 1)%state(:, 1:child_cells), &
        solution%hierarchy%interfaces(relation), &
        solution%levels(relation)%state(:, 1:parent_cells), local_ok)
      if (.not. local_ok) return
      call recover_level_temperatures_1d( &
        species, solution%levels(relation)%state, &
        solution%levels(relation)%temperature, parent_cells, local_ok)
      if (.not. local_ok) return
    end do
    call refresh_multilevel_ghosts(species, config, solution, local_ok)
    if (.not. local_ok) return
    ok = .true.
  end subroutine average_down_all_levels

  subroutine refresh_multilevel_ghosts(species, config, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_multilevel_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok

    logical :: local_ok
    integer :: level, nx

    ok = .false.
    nx = solution%hierarchy%base_cells
    call fill_physical_ghosts_1d( &
      solution%levels(1)%state, solution%levels(1)%temperature, nx, &
      config%boundary_condition, local_ok)
    if (.not. local_ok) return
    if (uses_ppm_reconstruction(config)) then
      call fill_physical_wide_ghosts_1d( &
        solution%levels(1)%state, solution%levels(1)%temperature, nx, &
        config%boundary_condition, &
        solution%levels(1)%left_ghost_state, &
        solution%levels(1)%right_ghost_state, &
        solution%levels(1)%left_ghost_temperature, &
        solution%levels(1)%right_ghost_temperature, local_ok)
      if (.not. local_ok) return
    end if
    do level = 2, solution%level_count()
      call fill_fine_ghosts_1d( &
        species, solution%hierarchy%interfaces(level - 1), &
        solution%levels(level - 1)%state, &
        solution%levels(level - 1)%state, 1.0_dp, &
        solution%levels(level)%state, &
        solution%levels(level)%temperature, local_ok)
      if (.not. local_ok) return
      if (uses_ppm_reconstruction(config)) then
        call fill_fine_wide_ghosts_1d( &
          species, solution%hierarchy%interfaces(level - 1), &
          solution%levels(level - 1)%state, &
          solution%levels(level - 1)%state, 1.0_dp, &
          solution%levels(level)%left_ghost_state, &
          solution%levels(level)%right_ghost_state, &
          solution%levels(level)%left_ghost_temperature, &
          solution%levels(level)%right_ghost_temperature, local_ok)
        if (.not. local_ok) return
      end if
    end do
    ok = .true.
  end subroutine refresh_multilevel_ghosts

  pure logical function uses_ppm_reconstruction(config) result(enabled)
    type(reactive_1d_config), intent(in) :: config

    enabled = trim(config%amr_reconstruction) == "ppm" .or. &
      trim(config%amr_reconstruction) == "characteristic_ppm"
  end function uses_ppm_reconstruction

  subroutine multilevel_reactive_integrals_1d(solution, integral, ok)
    type(amr_multilevel_reactive_solution_1d), intent(in) :: solution
    real(dp), intent(out) :: integral(:)
    logical, intent(out) :: ok

    type(amr_level_field_1d), allocatable :: fields(:)
    integer :: level, nx, nvar

    integral = 0.0_dp
    ok = solution%is_valid()
    if (.not. ok) return
    nvar = size(solution%levels(1)%state, 1)
    if (size(integral) /= nvar) then
      ok = .false.
      return
    end if
    allocate(fields(solution%level_count()))
    do level = 1, solution%level_count()
      nx = solution%hierarchy%level_cell_count(level - 1)
      allocate(fields(level)%values(nvar, nx))
      fields(level)%values = solution%levels(level)%state(:, 1:nx)
    end do
    call composite_integral_multilevel_1d( &
      fields, solution%hierarchy, integral, ok)
  end subroutine multilevel_reactive_integrals_1d

  subroutine write_multilevel_reactive_1d_csv(path, species, solution, ok)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(amr_multilevel_reactive_solution_1d), intent(in) :: solution
    logical, intent(out) :: ok

    real(dp), allocatable :: q(:)
    logical :: local_ok
    integer :: unit, status, k

    ok = .false.
    if (.not. solution%is_valid()) return
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
    call write_composite_level( &
      unit, 1, species, solution, q, local_ok)
    close(unit)
    if (.not. local_ok) return
    ok = .true.
  end subroutine write_multilevel_reactive_1d_csv

  recursive subroutine write_composite_level( &
      unit, level, species, solution, q, ok)
    integer, intent(in) :: unit, level
    type(nasa7_species), intent(in) :: species(:)
    type(amr_multilevel_reactive_solution_1d), intent(in) :: solution
    real(dp), intent(out) :: q(:)
    logical, intent(out) :: ok

    logical :: local_ok
    integer :: cell, first_after_child, last_before_child, nx

    ok = .false.
    nx = solution%hierarchy%level_cell_count(level - 1)
    if (level == solution%level_count()) then
      do cell = 1, nx
        call write_multilevel_cell( &
          unit, level, cell, species, solution, q, local_ok)
        if (.not. local_ok) return
      end do
      ok = .true.
      return
    end if
    last_before_child = &
      solution%hierarchy%interfaces(level)%fine_coarse_lower - 1
    first_after_child = &
      solution%hierarchy%interfaces(level)%fine_coarse_upper + 1
    do cell = 1, last_before_child
      call write_multilevel_cell( &
        unit, level, cell, species, solution, q, local_ok)
      if (.not. local_ok) return
    end do
    call write_composite_level( &
      unit, level + 1, species, solution, q, local_ok)
    if (.not. local_ok) return
    do cell = first_after_child, nx
      call write_multilevel_cell( &
        unit, level, cell, species, solution, q, local_ok)
      if (.not. local_ok) return
    end do
    ok = .true.
  end subroutine write_composite_level

  subroutine write_multilevel_cell( &
      unit, level, cell, species, solution, q, ok)
    integer, intent(in) :: unit, level, cell
    type(nasa7_species), intent(in) :: species(:)
    type(amr_multilevel_reactive_solution_1d), intent(in) :: solution
    real(dp), intent(out) :: q(:)
    logical, intent(out) :: ok

    real(dp) :: x_lower, x_upper, dx, x
    logical :: bounds_ok

    call solution%hierarchy%level_bounds( &
      level - 1, x_lower, x_upper, bounds_ok)
    if (.not. bounds_ok) then
      ok = .false.
      return
    end if
    dx = solution%hierarchy%level_dx(level - 1)
    x = x_lower + (real(cell, dp) - 0.5_dp) * dx
    if (x < x_lower .or. x > x_upper) then
      ok = .false.
      return
    end if
    call write_amr_cell( &
      unit, level - 1, dx, solution%time, x, species, &
      solution%levels(level)%state(:, cell), &
      solution%levels(level)%temperature(cell), q, ok)
  end subroutine write_multilevel_cell

  pure logical function same_hierarchy(left, right) result(same)
    type(amr_multilevel_hierarchy_1d), intent(in) :: left, right

    integer :: relation

    same = left%is_valid() .and. right%is_valid()
    if (.not. same) return
    same = left%base_cells == right%base_cells .and. &
      size(left%interfaces) == size(right%interfaces)
    if (.not. same) return
    do relation = 1, size(left%interfaces)
      same = &
        left%interfaces(relation)%fine_coarse_lower == &
          right%interfaces(relation)%fine_coarse_lower .and. &
        left%interfaces(relation)%fine_coarse_upper == &
          right%interfaces(relation)%fine_coarse_upper .and. &
        left%interfaces(relation)%refinement_ratio == &
          right%interfaces(relation)%refinement_ratio
      if (.not. same) return
    end do
  end function same_hierarchy

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

  pure function level_boundary(config, level) result(boundary)
    type(reactive_1d_config), intent(in) :: config
    integer, intent(in) :: level
    character(len=32) :: boundary

    if (level == 1) then
      boundary = config%boundary_condition
    else
      boundary = "coarse_fine"
    end if
  end function level_boundary

  pure function chemistry_boundary(config, level) result(boundary)
    type(reactive_1d_config), intent(in) :: config
    integer, intent(in) :: level
    character(len=32) :: boundary

    if (level == 1) then
      boundary = config%boundary_condition
    else
      boundary = "outflow"
    end if
  end function chemistry_boundary

end module amr_multilevel_reactive_1d_mod
