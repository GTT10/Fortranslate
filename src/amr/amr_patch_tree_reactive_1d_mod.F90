module amr_patch_tree_reactive_1d_mod
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: gas_transport_species
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_cfl_timestep, reactive_transport_timestep, &
    initialize_reactive_1d, advance_reactive_chemistry
  use amr_hierarchy_1d_mod, only: &
    amr_two_level_hierarchy_1d, amr_level_field_1d, &
    accumulate_coarse_flux_1d, accumulate_fine_flux_1d
  use amr_multipatch_1d_mod, only: synchronize_patch_set_1d
  use amr_patch_tree_1d_mod, only: &
    amr_patch_level_plan_1d, amr_patch_tree_hierarchy_1d, &
    amr_patch_tree_level_fields_1d, &
    amr_patch_tree_relation_flux_registers_1d, &
    initialize_patch_tree_1d, prolong_patch_tree_1d, &
    average_down_patch_tree_1d, &
    initialize_patch_tree_flux_registers_1d, &
    composite_integral_patch_tree_1d, patch_tree_child_geometry_1d
  use amr_reactive_1d_mod, only: &
    amr_ppm_ghost_width, advance_amr_level_1d, advance_transport_level_1d, &
    recover_level_temperatures_1d, fill_physical_ghosts_1d, &
    fill_fine_ghosts_1d, fill_fine_wide_ghosts_1d
  implicit none
  private

  type, public :: amr_patch_tree_reactive_patch_1d
    real(dp), allocatable :: state(:, :)
    real(dp), allocatable :: temperature(:)
    real(dp), allocatable :: left_ghost_state(:, :)
    real(dp), allocatable :: right_ghost_state(:, :)
    real(dp), allocatable :: left_ghost_temperature(:)
    real(dp), allocatable :: right_ghost_temperature(:)
  end type amr_patch_tree_reactive_patch_1d

  type, public :: amr_patch_tree_reactive_level_1d
    type(amr_patch_tree_reactive_patch_1d), allocatable :: patches(:)
  end type amr_patch_tree_reactive_level_1d

  type, public :: amr_patch_tree_reactive_solution_1d
    type(amr_patch_tree_hierarchy_1d) :: hierarchy
    type(amr_patch_tree_reactive_level_1d), allocatable :: levels(:)
    integer, allocatable :: level_advances(:)
    integer, allocatable :: transport_level_advances(:)
    real(dp) :: time = 0.0_dp
    integer :: steps = 0
    integer :: regrid_evaluations = 0
    integer :: regrids = 0
    integer :: overlap_cells_transferred = 0
  contains
    procedure :: level_count => patch_tree_reactive_level_count
    procedure :: is_valid => patch_tree_reactive_is_valid
  end type amr_patch_tree_reactive_solution_1d

  public :: initialize_patch_tree_reactive_1d
  public :: patch_tree_reactive_timestep_1d
  public :: advance_patch_tree_reactive_1d
  public :: advance_patch_tree_reactive_hydro_1d
  public :: regrid_patch_tree_reactive_1d
  public :: patch_tree_reactive_integrals_1d

contains

  pure integer function patch_tree_reactive_level_count(self) result(count)
    class(amr_patch_tree_reactive_solution_1d), intent(in) :: self

    count = 0
    if (allocated(self%levels)) count = size(self%levels)
  end function patch_tree_reactive_level_count

  pure logical function patch_tree_reactive_is_valid(self) result(valid)
    class(amr_patch_tree_reactive_solution_1d), intent(in) :: self

    type(amr_two_level_hierarchy_1d) :: geometry
    logical :: local_ok
    integer :: level, patch, nx, nvar

    valid = self%hierarchy%is_valid() .and. allocated(self%levels) .and. &
      allocated(self%level_advances) .and. &
      allocated(self%transport_level_advances)
    if (.not. valid) return
    valid = size(self%levels) == self%hierarchy%level_count() .and. &
      size(self%level_advances) == size(self%levels) .and. &
      size(self%transport_level_advances) == size(self%levels) .and. &
      all(self%level_advances >= 0) .and. self%time >= 0.0_dp .and. &
      all(self%transport_level_advances >= 0) .and. self%steps >= 0 .and. &
      self%regrid_evaluations >= 0 .and. self%regrids >= 0 .and. &
      self%overlap_cells_transferred >= 0
    if (.not. valid .or. size(self%levels) < 1) return
    valid = allocated(self%levels(1)%patches) .and. &
      size(self%levels(1)%patches) == 1 .and. &
      allocated(self%levels(1)%patches(1)%state)
    if (.not. valid) return
    nvar = size(self%levels(1)%patches(1)%state, 1)
    valid = nvar >= 1
    if (.not. valid) return

    do level = 1, size(self%levels)
      valid = allocated(self%levels(level)%patches) .and. &
        size(self%levels(level)%patches) == &
          self%hierarchy%level_patch_count(level - 1)
      if (.not. valid) return
      do patch = 1, size(self%levels(level)%patches)
        if (level == 1) then
          nx = self%hierarchy%base_cells
        else
          call patch_tree_child_geometry_1d( &
            self%hierarchy%relations(level - 1), patch, geometry, local_ok)
          if (.not. local_ok) then
            valid = .false.
            return
          end if
          nx = geometry%fine%cell_count()
        end if
        valid = reactive_patch_is_valid( &
          self%levels(level)%patches(patch), nvar, nx)
        if (.not. valid) return
      end do
    end do
  end function patch_tree_reactive_is_valid

  pure logical function reactive_patch_is_valid( &
      patch, nvar, nx) result(valid)
    type(amr_patch_tree_reactive_patch_1d), intent(in) :: patch
    integer, intent(in) :: nvar, nx

    valid = nvar >= 1 .and. nx >= 1 .and. allocated(patch%state) .and. &
      allocated(patch%temperature) .and. &
      allocated(patch%left_ghost_state) .and. &
      allocated(patch%right_ghost_state) .and. &
      allocated(patch%left_ghost_temperature) .and. &
      allocated(patch%right_ghost_temperature)
    if (.not. valid) return
    valid = size(patch%state, 1) == nvar .and. &
      lbound(patch%state, 2) == 0 .and. ubound(patch%state, 2) == nx + 1 .and. &
      lbound(patch%temperature, 1) == 0 .and. &
      ubound(patch%temperature, 1) == nx + 1 .and. &
      size(patch%left_ghost_state, 1) == nvar .and. &
      size(patch%right_ghost_state, 1) == nvar .and. &
      size(patch%left_ghost_state, 2) == amr_ppm_ghost_width .and. &
      size(patch%right_ghost_state, 2) == amr_ppm_ghost_width .and. &
      size(patch%left_ghost_temperature) == amr_ppm_ghost_width .and. &
      size(patch%right_ghost_temperature) == amr_ppm_ghost_width
  end function reactive_patch_is_valid

  subroutine initialize_patch_tree_reactive_1d( &
      species, config, plans, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_patch_level_plan_1d), intent(in) :: plans(:)
    type(amr_patch_tree_reactive_solution_1d), intent(out) :: solution
    logical, intent(out) :: ok

    type(amr_patch_tree_level_fields_1d), allocatable :: fields(:)
    type(amr_two_level_hierarchy_1d) :: geometry
    real(dp), allocatable :: root_state(:, :), root_temperature(:)
    real(dp) :: root_dx
    logical :: local_ok
    integer :: level, patch, nx, nvar

    ok = .false.
    if (size(species) < 1 .or. .not. config%amr_enabled) return
    call initialize_reactive_1d( &
      species, config, root_state, root_temperature, root_dx, local_ok)
    if (.not. local_ok .or. root_dx <= 0.0_dp) return
    call initialize_patch_tree_1d( &
      config%nx, config%x_lower, config%x_upper, plans, &
      solution%hierarchy, local_ok)
    if (.not. local_ok) return
    if (.not. patch_tree_children_are_interior(solution%hierarchy)) return
    call prolong_patch_tree_1d( &
      root_state(:, 1:config%nx), solution%hierarchy, fields, local_ok)
    if (.not. local_ok) return

    nvar = reactive_nvar(size(species))
    allocate(solution%levels(solution%hierarchy%level_count()))
    allocate(solution%level_advances(solution%hierarchy%level_count()))
    allocate(solution%transport_level_advances( &
      solution%hierarchy%level_count()))
    solution%level_advances = 0
    solution%transport_level_advances = 0
    do level = 1, solution%level_count()
      allocate(solution%levels(level)%patches( &
        solution%hierarchy%level_patch_count(level - 1)))
      do patch = 1, size(solution%levels(level)%patches)
        if (level == 1) then
          nx = config%nx
        else
          call patch_tree_child_geometry_1d( &
            solution%hierarchy%relations(level - 1), patch, &
            geometry, local_ok)
          if (.not. local_ok) return
          nx = geometry%fine%cell_count()
        end if
        call allocate_reactive_patch( &
          solution%levels(level)%patches(patch), nvar, nx)
        solution%levels(level)%patches(patch)%state(:, 1:nx) = &
          fields(level)%patches(patch)%values
        if (level == 1) then
          solution%levels(level)%patches(patch)%state = root_state
          solution%levels(level)%patches(patch)%temperature = root_temperature
        else
          call recover_level_temperatures_1d( &
            species, solution%levels(level)%patches(patch)%state, &
            solution%levels(level)%patches(patch)%temperature, nx, local_ok)
          if (.not. local_ok) return
        end if
      end do
    end do
    call refresh_patch_tree_ghosts(species, config, solution, local_ok)
    if (.not. local_ok) return
    solution%time = 0.0_dp
    solution%steps = 0
    solution%regrid_evaluations = 0
    solution%regrids = 0
    solution%overlap_cells_transferred = 0
    ok = solution%is_valid()
  end subroutine initialize_patch_tree_reactive_1d

  subroutine allocate_reactive_patch(patch, nvar, nx)
    type(amr_patch_tree_reactive_patch_1d), intent(out) :: patch
    integer, intent(in) :: nvar, nx

    allocate(patch%state(nvar, 0:nx + 1))
    allocate(patch%temperature(0:nx + 1))
    allocate(patch%left_ghost_state(nvar, amr_ppm_ghost_width))
    allocate(patch%right_ghost_state(nvar, amr_ppm_ghost_width))
    allocate(patch%left_ghost_temperature(amr_ppm_ghost_width))
    allocate(patch%right_ghost_temperature(amr_ppm_ghost_width))
    patch%state = 0.0_dp
    patch%temperature = 0.0_dp
    patch%left_ghost_state = 0.0_dp
    patch%right_ghost_state = 0.0_dp
    patch%left_ghost_temperature = 0.0_dp
    patch%right_ghost_temperature = 0.0_dp
  end subroutine allocate_reactive_patch

  subroutine patch_tree_reactive_timestep_1d( &
      species, config, solution, dt, ok, transport)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_patch_tree_reactive_solution_1d), intent(in) :: solution
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok
    type(gas_transport_species), intent(in), optional :: transport(:)

    real(dp) :: local_dt, transport_dt, maximum_diffusivity
    real(dp) :: scale, transport_scale, dx
    logical :: local_ok
    integer :: level, patch, nx

    dt = huge(1.0_dp)
    ok = solution%is_valid() .and. size(species) >= 1
    if (.not. ok) return
    if (config%transport_enabled .and. .not. present(transport)) then
      ok = .false.
      return
    end if
    scale = 1.0_dp
    do level = 1, solution%level_count()
      if (level > 1) scale = scale * real( &
        solution%hierarchy%relations(level - 1)%refinement_ratio, dp)
      dx = solution%hierarchy%level_dx(level - 1)
      do patch = 1, size(solution%levels(level)%patches)
        nx = size(solution%levels(level)%patches(patch)%state, 2) - 2
        call reactive_cfl_timestep( &
          species, solution%levels(level)%patches(patch)%state, &
          solution%levels(level)%patches(patch)%temperature, nx, dx, &
          config%cfl, local_dt, local_ok)
        if (.not. local_ok) then
          ok = .false.
          return
        end if
        dt = min(dt, scale * local_dt)
        if (config%transport_enabled) then
          call reactive_transport_timestep( &
            species, transport, &
            solution%levels(level)%patches(patch)%state, &
            solution%levels(level)%patches(patch)%temperature, nx, dx, &
            config%transport_cfl, config%viscosity_enabled, &
            config%thermal_conduction_enabled, &
            config%species_diffusion_enabled, transport_dt, &
            maximum_diffusivity, local_ok)
          if (.not. local_ok) then
            ok = .false.
            return
          end if
          transport_scale = scale * scale
          dt = min(dt, transport_scale * transport_dt)
        end if
      end do
    end do
    ok = dt > 0.0_dp .and. dt < huge(1.0_dp)
  end subroutine patch_tree_reactive_timestep_1d

  subroutine advance_patch_tree_reactive_1d( &
      species, reactions, config, dt, solution, ok, transport)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: dt
    type(amr_patch_tree_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok
    type(gas_transport_species), intent(in), optional :: transport(:)

    type(amr_patch_tree_reactive_solution_1d) :: backup
    logical :: local_ok

    ok = .false.
    if (dt <= 0.0_dp .or. .not. solution%is_valid()) return
    if (config%transport_enabled .and. .not. present(transport)) return
    backup = solution
    if (config%chemistry_enabled) then
      call advance_patch_tree_chemistry( &
        species, reactions, config, 0.5_dp * dt, solution, local_ok)
      if (.not. local_ok) then
        solution = backup
        return
      end if
    end if
    if (config%transport_enabled) then
      call advance_patch_tree_transport( &
        species, transport, config, 0.5_dp * dt, solution, local_ok)
      if (.not. local_ok) then
        solution = backup
        return
      end if
    end if
    call advance_patch_tree_reactive_hydro_1d( &
      species, config, dt, solution, local_ok)
    if (.not. local_ok) then
      solution = backup
      return
    end if
    if (config%transport_enabled) then
      call advance_patch_tree_transport( &
        species, transport, config, 0.5_dp * dt, solution, local_ok)
      if (.not. local_ok) then
        solution = backup
        return
      end if
    end if
    if (config%chemistry_enabled) then
      call advance_patch_tree_chemistry( &
        species, reactions, config, 0.5_dp * dt, solution, local_ok)
      if (.not. local_ok) then
        solution = backup
        return
      end if
    end if
    call refresh_patch_tree_ghosts(species, config, solution, local_ok)
    if (.not. local_ok) then
      solution = backup
      return
    end if
    ok = solution%is_valid()
    if (.not. ok) solution = backup
  end subroutine advance_patch_tree_reactive_1d

  subroutine regrid_patch_tree_reactive_1d( &
      species, config, plans, solution, changed, transferred_cells, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_patch_level_plan_1d), intent(in) :: plans(:)
    type(amr_patch_tree_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: changed
    integer, intent(out) :: transferred_cells
    logical, intent(out) :: ok

    type(amr_patch_tree_reactive_solution_1d) :: backup, rebuilt
    type(amr_patch_tree_level_fields_1d), allocatable :: fields(:)
    real(dp) :: tolerance
    logical :: local_ok
    integer :: common_levels

    changed = .false.
    transferred_cells = 0
    ok = solution%is_valid() .and. config%amr_enabled .and. &
      config%nx == solution%hierarchy%base_cells
    if (.not. ok) return
    tolerance = 128.0_dp * epsilon(1.0_dp) * max( &
      1.0_dp, abs(config%x_lower), abs(config%x_upper), &
      abs(solution%hierarchy%x_lower), abs(solution%hierarchy%x_upper))
    ok = abs(config%x_lower - solution%hierarchy%x_lower) <= tolerance .and. &
      abs(config%x_upper - solution%hierarchy%x_upper) <= tolerance
    if (.not. ok) return
    backup = solution
    call initialize_patch_tree_reactive_1d( &
      species, config, plans, rebuilt, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    changed = .not. same_patch_tree_hierarchy( &
      backup%hierarchy, rebuilt%hierarchy)
    if (.not. changed) then
      solution%regrid_evaluations = backup%regrid_evaluations + 1
      ok = .true.
      return
    end if

    call average_down_patch_tree_solution(species, config, solution, local_ok)
    if (.not. local_ok) then
      solution = backup
      ok = .false.
      return
    end if
    call prolong_patch_tree_1d( &
      solution%levels(1)%patches(1)%state(:, 1:config%nx), &
      rebuilt%hierarchy, fields, local_ok)
    if (.not. local_ok) then
      solution = backup
      ok = .false.
      return
    end if
    call install_patch_tree_fields(species, fields, rebuilt, local_ok)
    if (.not. local_ok) then
      solution = backup
      ok = .false.
      return
    end if
    call transfer_patch_tree_overlap( &
      solution, rebuilt, transferred_cells, local_ok)
    if (.not. local_ok) then
      solution = backup
      ok = .false.
      return
    end if
    call average_down_patch_tree_solution( &
      species, config, rebuilt, local_ok)
    if (.not. local_ok) then
      solution = backup
      ok = .false.
      return
    end if

    rebuilt%time = backup%time
    rebuilt%steps = backup%steps
    common_levels = min( &
      backup%level_count(), rebuilt%level_count())
    rebuilt%level_advances(1:common_levels) = &
      backup%level_advances(1:common_levels)
    rebuilt%transport_level_advances(1:common_levels) = &
      backup%transport_level_advances(1:common_levels)
    rebuilt%regrid_evaluations = backup%regrid_evaluations + 1
    rebuilt%regrids = backup%regrids + 1
    rebuilt%overlap_cells_transferred = &
      backup%overlap_cells_transferred + transferred_cells
    solution = rebuilt
    ok = solution%is_valid()
    if (.not. ok) solution = backup
  end subroutine regrid_patch_tree_reactive_1d

  subroutine average_down_patch_tree_solution(species, config, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_patch_tree_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok

    type(amr_patch_tree_level_fields_1d), allocatable :: fields(:)
    logical :: local_ok
    integer :: level, patch, nx

    call extract_patch_tree_fields(solution, fields, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call average_down_patch_tree_1d(fields, solution%hierarchy, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    do level = 1, solution%level_count()
      do patch = 1, size(solution%levels(level)%patches)
        nx = size(solution%levels(level)%patches(patch)%state, 2) - 2
        solution%levels(level)%patches(patch)%state(:, 1:nx) = &
          fields(level)%patches(patch)%values
        if (level > size(solution%hierarchy%relations)) cycle
        if (solution%hierarchy%relations(level)% &
            child_sets(patch)%patch_count() == 0) cycle
        call recover_level_temperatures_1d( &
          species, solution%levels(level)%patches(patch)%state, &
          solution%levels(level)%patches(patch)%temperature, nx, local_ok)
        if (.not. local_ok) then
          ok = .false.
          return
        end if
      end do
    end do
    call refresh_patch_tree_ghosts(species, config, solution, local_ok)
    ok = local_ok
  end subroutine average_down_patch_tree_solution

  subroutine install_patch_tree_fields(species, fields, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_level_fields_1d), intent(in) :: fields(:)
    type(amr_patch_tree_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok

    logical :: local_ok
    integer :: level, patch, nx

    ok = size(fields) == solution%level_count()
    if (.not. ok) return
    do level = 1, solution%level_count()
      if (.not. allocated(fields(level)%patches)) then
        ok = .false.
        return
      end if
      if (size(fields(level)%patches) /= &
          size(solution%levels(level)%patches)) then
        ok = .false.
        return
      end if
      do patch = 1, size(solution%levels(level)%patches)
        nx = size(solution%levels(level)%patches(patch)%state, 2) - 2
        if (.not. allocated(fields(level)%patches(patch)%values)) then
          ok = .false.
          return
        end if
        if (size(fields(level)%patches(patch)%values, 1) /= &
              size(solution%levels(level)%patches(patch)%state, 1) .or. &
            size(fields(level)%patches(patch)%values, 2) /= nx) then
          ok = .false.
          return
        end if
        solution%levels(level)%patches(patch)%state(:, 1:nx) = &
          fields(level)%patches(patch)%values
        call recover_level_temperatures_1d( &
          species, solution%levels(level)%patches(patch)%state, &
          solution%levels(level)%patches(patch)%temperature, nx, local_ok)
        if (.not. local_ok) then
          ok = .false.
          return
        end if
      end do
    end do
    ok = .true.
  end subroutine install_patch_tree_fields

  subroutine transfer_patch_tree_overlap( &
      old_solution, new_solution, transferred_cells, ok)
    type(amr_patch_tree_reactive_solution_1d), intent(in) :: old_solution
    type(amr_patch_tree_reactive_solution_1d), intent(inout) :: new_solution
    integer, intent(out) :: transferred_cells
    logical, intent(out) :: ok

    type(amr_two_level_hierarchy_1d) :: old_geometry, new_geometry
    real(dp) :: old_lower, old_upper, new_lower, new_upper
    real(dp) :: overlap_lower, overlap_upper, old_dx, new_dx, tolerance
    real(dp) :: old_offset, new_offset
    logical :: local_ok
    integer :: level, old_patch, new_patch, common_levels
    integer :: old_first, new_first, cell_count

    transferred_cells = 0
    ok = old_solution%is_valid() .and. new_solution%is_valid()
    if (.not. ok) return
    common_levels = min(old_solution%level_count(), new_solution%level_count())
    do level = 2, common_levels
      old_dx = old_solution%hierarchy%level_dx(level - 1)
      new_dx = new_solution%hierarchy%level_dx(level - 1)
      tolerance = 128.0_dp * epsilon(1.0_dp) * &
        max(1.0_dp, abs(old_dx), abs(new_dx))
      if (abs(old_dx - new_dx) > tolerance) cycle
      do old_patch = 1, size(old_solution%levels(level)%patches)
        call patch_tree_child_geometry_1d( &
          old_solution%hierarchy%relations(level - 1), old_patch, &
          old_geometry, local_ok)
        if (.not. local_ok) then
          ok = .false.
          return
        end if
        call patch_physical_bounds( &
          old_geometry, old_lower, old_upper)
        do new_patch = 1, size(new_solution%levels(level)%patches)
          call patch_tree_child_geometry_1d( &
            new_solution%hierarchy%relations(level - 1), new_patch, &
            new_geometry, local_ok)
          if (.not. local_ok) then
            ok = .false.
            return
          end if
          call patch_physical_bounds( &
            new_geometry, new_lower, new_upper)
          overlap_lower = max(old_lower, new_lower)
          overlap_upper = min(old_upper, new_upper)
          if (overlap_upper <= overlap_lower + tolerance) cycle
          old_offset = (overlap_lower - old_lower) / old_dx
          new_offset = (overlap_lower - new_lower) / new_dx
          if (abs(old_offset - real(nint(old_offset), dp)) > &
                tolerance / old_dx .or. &
              abs(new_offset - real(nint(new_offset), dp)) > &
                tolerance / new_dx) then
            ok = .false.
            return
          end if
          old_first = nint(old_offset) + 1
          new_first = nint(new_offset) + 1
          cell_count = nint((overlap_upper - overlap_lower) / old_dx)
          if (cell_count < 1) cycle
          new_solution%levels(level)%patches(new_patch)%state(:, &
            new_first:new_first + cell_count - 1) = &
            old_solution%levels(level)%patches(old_patch)%state(:, &
              old_first:old_first + cell_count - 1)
          new_solution%levels(level)%patches(new_patch)%temperature( &
            new_first:new_first + cell_count - 1) = &
            old_solution%levels(level)%patches(old_patch)%temperature( &
              old_first:old_first + cell_count - 1)
          transferred_cells = transferred_cells + cell_count
        end do
      end do
    end do
    ok = .true.
  end subroutine transfer_patch_tree_overlap

  pure logical function same_patch_tree_hierarchy(first, second) result(same)
    type(amr_patch_tree_hierarchy_1d), intent(in) :: first, second

    real(dp) :: tolerance
    integer :: relation, parent, child

    same = first%is_valid() .and. second%is_valid()
    if (.not. same) return
    tolerance = 128.0_dp * epsilon(1.0_dp) * max( &
      1.0_dp, abs(first%x_lower), abs(first%x_upper), &
      abs(second%x_lower), abs(second%x_upper))
    same = first%base_cells == second%base_cells .and. &
      size(first%relations) == size(second%relations) .and. &
      abs(first%x_lower - second%x_lower) <= tolerance .and. &
      abs(first%x_upper - second%x_upper) <= tolerance
    if (.not. same) return
    do relation = 1, size(first%relations)
      same = first%relations(relation)%refinement_ratio == &
          second%relations(relation)%refinement_ratio .and. &
        first%relations(relation)%parent_patch_count() == &
          second%relations(relation)%parent_patch_count()
      if (.not. same) return
      do parent = 1, first%relations(relation)%parent_patch_count()
        same = first%relations(relation)%child_sets(parent)%patch_count() == &
          second%relations(relation)%child_sets(parent)%patch_count()
        if (.not. same) return
        do child = 1, first%relations(relation)% &
            child_sets(parent)%patch_count()
          same = first%relations(relation)%child_sets(parent)% &
              patches(child)%fine_coarse_lower == &
                second%relations(relation)%child_sets(parent)% &
                  patches(child)%fine_coarse_lower .and. &
            first%relations(relation)%child_sets(parent)% &
              patches(child)%fine_coarse_upper == &
                second%relations(relation)%child_sets(parent)% &
                  patches(child)%fine_coarse_upper
          if (.not. same) return
        end do
      end do
    end do
  end function same_patch_tree_hierarchy

  pure subroutine patch_physical_bounds(geometry, lower, upper)
    type(amr_two_level_hierarchy_1d), intent(in) :: geometry
    real(dp), intent(out) :: lower, upper

    lower = geometry%x_lower + &
      real(geometry%fine_coarse_lower - 1, dp) * geometry%coarse_dx
    upper = geometry%x_lower + &
      real(geometry%fine_coarse_upper, dp) * geometry%coarse_dx
  end subroutine patch_physical_bounds

  subroutine advance_patch_tree_reactive_hydro_1d( &
      species, config, dt, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: dt
    type(amr_patch_tree_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok

    type(amr_patch_tree_reactive_solution_1d) :: backup
    type(amr_patch_tree_relation_flux_registers_1d), allocatable :: registers(:)
    real(dp), allocatable :: left_integral(:), right_integral(:)
    logical :: local_ok
    integer :: nvar

    ok = dt > 0.0_dp .and. solution%is_valid()
    if (.not. ok) return
    backup = solution
    nvar = size(solution%levels(1)%patches(1)%state, 1)
    allocate(left_integral(nvar), right_integral(nvar))
    call initialize_patch_tree_flux_registers_1d( &
      solution%hierarchy, nvar, registers, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call advance_patch_recursive( &
      species, config, solution, registers, 1, 1, dt, &
      left_integral, right_integral, local_ok)
    if (.not. local_ok) then
      solution = backup
      ok = .false.
      return
    end if
    call refresh_patch_tree_ghosts(species, config, solution, local_ok)
    if (.not. local_ok) then
      solution = backup
      ok = .false.
      return
    end if
    solution%time = solution%time + dt
    solution%steps = solution%steps + 1
    ok = solution%is_valid()
    if (.not. ok) solution = backup
  end subroutine advance_patch_tree_reactive_hydro_1d

  subroutine advance_patch_tree_chemistry( &
      species, reactions, config, interval, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: interval
    type(amr_patch_tree_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok

    type(amr_patch_tree_level_fields_1d), allocatable :: fields(:)
    logical :: local_ok
    integer :: level, patch, nx

    ok = solution%is_valid() .and. interval >= 0.0_dp
    if (.not. ok) return
    do level = 1, solution%level_count()
      do patch = 1, size(solution%levels(level)%patches)
        nx = size(solution%levels(level)%patches(patch)%state, 2) - 2
        call advance_reactive_chemistry( &
          species, reactions, &
          solution%levels(level)%patches(patch)%state, &
          solution%levels(level)%patches(patch)%temperature, nx, interval, &
          config%chemistry_relative_tolerance, &
          config%chemistry_absolute_tolerance, &
          chemistry_boundary(config, level), local_ok)
        if (.not. local_ok) then
          ok = .false.
          return
        end if
      end do
    end do
    call extract_patch_tree_fields(solution, fields, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call average_down_patch_tree_1d( &
      fields, solution%hierarchy, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    do level = 1, solution%level_count()
      do patch = 1, size(solution%levels(level)%patches)
        nx = size(solution%levels(level)%patches(patch)%state, 2) - 2
        solution%levels(level)%patches(patch)%state(:, 1:nx) = &
          fields(level)%patches(patch)%values
        call recover_level_temperatures_1d( &
          species, solution%levels(level)%patches(patch)%state, &
          solution%levels(level)%patches(patch)%temperature, nx, local_ok)
        if (.not. local_ok) then
          ok = .false.
          return
        end if
      end do
    end do
    call refresh_patch_tree_ghosts(species, config, solution, local_ok)
    ok = local_ok
  end subroutine advance_patch_tree_chemistry

  subroutine advance_patch_tree_transport( &
      species, transport, config, interval, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: interval
    type(amr_patch_tree_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok

    type(amr_patch_tree_relation_flux_registers_1d), allocatable :: registers(:)
    real(dp), allocatable :: left_integral(:), right_integral(:)
    logical :: local_ok
    integer :: nvar

    ok = interval > 0.0_dp .and. solution%is_valid()
    if (.not. ok) return
    nvar = size(solution%levels(1)%patches(1)%state, 1)
    allocate(left_integral(nvar), right_integral(nvar))
    call initialize_patch_tree_flux_registers_1d( &
      solution%hierarchy, nvar, registers, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call advance_patch_transport_recursive( &
      species, transport, config, solution, registers, 1, 1, interval, &
      left_integral, right_integral, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call refresh_patch_tree_ghosts(species, config, solution, local_ok)
    ok = local_ok
  end subroutine advance_patch_tree_transport

  recursive subroutine advance_patch_transport_recursive( &
      species, transport, config, solution, registers, level, &
      parent_patch, interval, left_integral, right_integral, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_patch_tree_reactive_solution_1d), intent(inout) :: solution
    type(amr_patch_tree_relation_flux_registers_1d), &
      intent(inout) :: registers(:)
    integer, intent(in) :: level, parent_patch
    real(dp), intent(in) :: interval
    real(dp), intent(out) :: left_integral(:), right_integral(:)
    logical, intent(out) :: ok

    type(amr_level_field_1d), allocatable :: child_fields(:)
    real(dp), allocatable :: state_start(:, :), state_end(:, :), flux(:, :)
    real(dp), allocatable :: child_left(:), child_right(:)
    real(dp) :: child_interval, alpha, dx, boundary_distance
    logical :: local_ok, physical_boundary
    integer :: nx, nvar, ratio, subcycles, substep
    integer :: child, global_child, child_count

    ok = .false.
    left_integral = 0.0_dp
    right_integral = 0.0_dp
    if (interval <= 0.0_dp .or. level < 1 .or. &
        level > solution%level_count()) return
    if (parent_patch < 1 .or. &
        parent_patch > size(solution%levels(level)%patches)) return
    nx = size(solution%levels(level)%patches(parent_patch)%state, 2) - 2
    nvar = size(solution%levels(level)%patches(parent_patch)%state, 1)
    dx = solution%hierarchy%level_dx(level - 1)
    allocate(state_start(nvar, 0:nx + 1), state_end(nvar, 0:nx + 1))
    allocate(flux(nvar, 0:nx))
    state_start = solution%levels(level)%patches(parent_patch)%state
    physical_boundary = level == 1
    boundary_distance = dx
    if (.not. physical_boundary) boundary_distance = 0.5_dp * ( &
      solution%hierarchy%level_dx(level - 2) + dx)
    call advance_transport_level_1d( &
      species, transport, &
      solution%levels(level)%patches(parent_patch)%state, &
      solution%levels(level)%patches(parent_patch)%temperature, nx, dx, &
      interval, boundary_distance, config, physical_boundary, &
      patch_boundary(config, level), flux, local_ok)
    if (.not. local_ok) return
    solution%transport_level_advances(level) = &
      solution%transport_level_advances(level) + 1
    state_end = solution%levels(level)%patches(parent_patch)%state
    left_integral = interval * flux(:, 0)
    right_integral = interval * flux(:, nx)
    if (level > size(solution%hierarchy%relations)) then
      ok = .true.
      return
    end if

    child_count = solution%hierarchy%relations(level)% &
      child_sets(parent_patch)%patch_count()
    if (child_count == 0) then
      ok = .true.
      return
    end if
    ratio = solution%hierarchy%relations(level)%refinement_ratio
    subcycles = ratio * ratio
    child_interval = interval / real(subcycles, dp)
    allocate(child_left(nvar), child_right(nvar))
    do child = 1, child_count
      call accumulate_coarse_flux_1d( &
        registers(level)%parents(parent_patch)%children(child), &
        flux(:, solution%hierarchy%relations(level)% &
          child_sets(parent_patch)%patches(child)%fine_coarse_lower - 1), &
        flux(:, solution%hierarchy%relations(level)% &
          child_sets(parent_patch)%patches(child)%fine_coarse_upper), &
        interval, local_ok)
      if (.not. local_ok) return
    end do

    do substep = 1, subcycles
      alpha = (real(substep, dp) - 0.5_dp) / real(subcycles, dp)
      do child = 1, child_count
        global_child = solution%hierarchy%relations(level)% &
          child_index(parent_patch, child)
        call fill_one_child_ghosts( &
          species, config, solution%hierarchy%relations(level)% &
            child_sets(parent_patch)%patches(child), state_start, state_end, &
          alpha, solution%levels(level + 1)%patches(global_child), local_ok)
        if (.not. local_ok) return
      end do
      do child = 1, child_count
        global_child = solution%hierarchy%relations(level)% &
          child_index(parent_patch, child)
        call advance_patch_transport_recursive( &
          species, transport, config, solution, registers, level + 1, &
          global_child, child_interval, child_left, child_right, local_ok)
        if (.not. local_ok) return
        call accumulate_fine_flux_1d( &
          registers(level)%parents(parent_patch)%children(child), &
          child_left / child_interval, child_right / child_interval, &
          child_interval, local_ok)
        if (.not. local_ok) return
      end do
    end do

    allocate(child_fields(child_count))
    do child = 1, child_count
      global_child = solution%hierarchy%relations(level)% &
        child_index(parent_patch, child)
      nx = size(solution%levels(level + 1)% &
        patches(global_child)%state, 2) - 2
      child_fields(child)%values = solution%levels(level + 1)% &
        patches(global_child)%state(:, 1:nx)
    end do
    nx = size(solution%levels(level)%patches(parent_patch)%state, 2) - 2
    call synchronize_patch_set_1d( &
      solution%levels(level)%patches(parent_patch)%state(:, 1:nx), &
      child_fields, &
      solution%hierarchy%relations(level)%child_sets(parent_patch), &
      registers(level)%parents(parent_patch)%children, local_ok)
    if (.not. local_ok) return
    call recover_level_temperatures_1d( &
      species, solution%levels(level)%patches(parent_patch)%state, &
      solution%levels(level)%patches(parent_patch)%temperature, nx, local_ok)
    if (.not. local_ok) return
    ok = .true.
  end subroutine advance_patch_transport_recursive

  recursive subroutine advance_patch_recursive( &
      species, config, solution, registers, level, parent_patch, interval, &
      left_integral, right_integral, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_patch_tree_reactive_solution_1d), intent(inout) :: solution
    type(amr_patch_tree_relation_flux_registers_1d), &
      intent(inout) :: registers(:)
    integer, intent(in) :: level, parent_patch
    real(dp), intent(in) :: interval
    real(dp), intent(out) :: left_integral(:), right_integral(:)
    logical, intent(out) :: ok

    type(amr_level_field_1d), allocatable :: child_fields(:)
    real(dp), allocatable :: state_start(:, :), state_end(:, :), flux(:, :)
    real(dp), allocatable :: child_left(:), child_right(:)
    real(dp) :: child_interval, alpha, dx
    logical :: local_ok, physical_boundary
    integer :: nx, nvar, ratio, substep, child, global_child, child_count

    ok = .false.
    left_integral = 0.0_dp
    right_integral = 0.0_dp
    if (interval <= 0.0_dp .or. level < 1 .or. &
        level > solution%level_count()) return
    if (parent_patch < 1 .or. &
        parent_patch > size(solution%levels(level)%patches)) return
    nx = size(solution%levels(level)%patches(parent_patch)%state, 2) - 2
    nvar = size(solution%levels(level)%patches(parent_patch)%state, 1)
    dx = solution%hierarchy%level_dx(level - 1)
    allocate(state_start(nvar, 0:nx + 1), state_end(nvar, 0:nx + 1))
    allocate(flux(nvar, 0:nx))
    state_start = solution%levels(level)%patches(parent_patch)%state
    physical_boundary = level == 1
    if (physical_boundary .or. .not. uses_ppm_reconstruction(config)) then
      call advance_amr_level_1d( &
        species, solution%levels(level)%patches(parent_patch)%state, &
        solution%levels(level)%patches(parent_patch)%temperature, nx, dx, &
        interval, config%amr_reconstruction, config%limiter, &
        config%riemann_solver, physical_boundary, &
        patch_boundary(config, level), flux, local_ok, &
        ppm_contact_steepening=config%ppm_contact_steepening, &
        ppm_shock_flattening=config%ppm_shock_flattening, &
        amr_hybrid_weno=config%amr_hybrid_weno, &
        amr_weno_scheme=config%amr_weno_scheme)
    else
      call advance_amr_level_1d( &
        species, solution%levels(level)%patches(parent_patch)%state, &
        solution%levels(level)%patches(parent_patch)%temperature, nx, dx, &
        interval, config%amr_reconstruction, config%limiter, &
        config%riemann_solver, physical_boundary, &
        patch_boundary(config, level), flux, local_ok, &
        solution%levels(level)%patches(parent_patch)%left_ghost_state, &
        solution%levels(level)%patches(parent_patch)%right_ghost_state, &
        solution%levels(level)%patches(parent_patch)% &
          left_ghost_temperature, &
        solution%levels(level)%patches(parent_patch)% &
          right_ghost_temperature, config%ppm_contact_steepening, &
        config%ppm_shock_flattening, config%amr_hybrid_weno, &
        config%amr_weno_scheme)
    end if
    if (.not. local_ok) return
    solution%level_advances(level) = solution%level_advances(level) + 1
    state_end = solution%levels(level)%patches(parent_patch)%state
    left_integral = interval * flux(:, 0)
    right_integral = interval * flux(:, nx)
    if (level > size(solution%hierarchy%relations)) then
      ok = .true.
      return
    end if

    child_count = solution%hierarchy%relations(level)% &
      child_sets(parent_patch)%patch_count()
    if (child_count == 0) then
      ok = .true.
      return
    end if
    ratio = solution%hierarchy%relations(level)%refinement_ratio
    child_interval = interval / real(ratio, dp)
    allocate(child_left(nvar), child_right(nvar))
    do child = 1, child_count
      call accumulate_coarse_flux_1d( &
        registers(level)%parents(parent_patch)%children(child), &
        flux(:, solution%hierarchy%relations(level)% &
          child_sets(parent_patch)%patches(child)%fine_coarse_lower - 1), &
        flux(:, solution%hierarchy%relations(level)% &
          child_sets(parent_patch)%patches(child)%fine_coarse_upper), &
        interval, local_ok)
      if (.not. local_ok) return
    end do

    do substep = 1, ratio
      if (trim(config%amr_reconstruction) == "pcm") then
        alpha = real(substep - 1, dp) / real(ratio, dp)
      else
        alpha = (real(substep, dp) - 0.5_dp) / real(ratio, dp)
      end if
      do child = 1, child_count
        global_child = solution%hierarchy%relations(level)% &
          child_index(parent_patch, child)
        call fill_one_child_ghosts( &
          species, config, solution%hierarchy%relations(level)% &
            child_sets(parent_patch)%patches(child), state_start, state_end, &
          alpha, solution%levels(level + 1)%patches(global_child), local_ok)
        if (.not. local_ok) return
      end do
      do child = 1, child_count
        global_child = solution%hierarchy%relations(level)% &
          child_index(parent_patch, child)
        call advance_patch_recursive( &
          species, config, solution, registers, level + 1, global_child, &
          child_interval, child_left, child_right, local_ok)
        if (.not. local_ok) return
        call accumulate_fine_flux_1d( &
          registers(level)%parents(parent_patch)%children(child), &
          child_left / child_interval, child_right / child_interval, &
          child_interval, local_ok)
        if (.not. local_ok) return
      end do
    end do

    allocate(child_fields(child_count))
    do child = 1, child_count
      global_child = solution%hierarchy%relations(level)% &
        child_index(parent_patch, child)
      nx = size(solution%levels(level + 1)% &
        patches(global_child)%state, 2) - 2
      child_fields(child)%values = solution%levels(level + 1)% &
        patches(global_child)%state(:, 1:nx)
    end do
    nx = size(solution%levels(level)%patches(parent_patch)%state, 2) - 2
    call synchronize_patch_set_1d( &
      solution%levels(level)%patches(parent_patch)%state(:, 1:nx), &
      child_fields, &
      solution%hierarchy%relations(level)%child_sets(parent_patch), &
      registers(level)%parents(parent_patch)%children, local_ok)
    if (.not. local_ok) return
    call recover_level_temperatures_1d( &
      species, solution%levels(level)%patches(parent_patch)%state, &
      solution%levels(level)%patches(parent_patch)%temperature, nx, local_ok)
    if (.not. local_ok) return
    ok = .true.
  end subroutine advance_patch_recursive

  subroutine patch_tree_reactive_integrals_1d(solution, integral, ok)
    type(amr_patch_tree_reactive_solution_1d), intent(in) :: solution
    real(dp), intent(out) :: integral(:)
    logical, intent(out) :: ok

    type(amr_patch_tree_level_fields_1d), allocatable :: fields(:)

    call extract_patch_tree_fields(solution, fields, ok)
    if (.not. ok) return
    call composite_integral_patch_tree_1d( &
      fields, solution%hierarchy, integral, ok)
  end subroutine patch_tree_reactive_integrals_1d

  subroutine extract_patch_tree_fields(solution, fields, ok)
    type(amr_patch_tree_reactive_solution_1d), intent(in) :: solution
    type(amr_patch_tree_level_fields_1d), allocatable, intent(out) :: fields(:)
    logical, intent(out) :: ok

    integer :: level, patch, nx

    ok = solution%is_valid()
    if (.not. ok) return
    allocate(fields(solution%level_count()))
    do level = 1, solution%level_count()
      allocate(fields(level)%patches(size(solution%levels(level)%patches)))
      do patch = 1, size(solution%levels(level)%patches)
        nx = size(solution%levels(level)%patches(patch)%state, 2) - 2
        fields(level)%patches(patch)%values = &
          solution%levels(level)%patches(patch)%state(:, 1:nx)
      end do
    end do
    ok = .true.
  end subroutine extract_patch_tree_fields

  subroutine refresh_patch_tree_ghosts(species, config, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_patch_tree_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok

    logical :: local_ok
    integer :: relation, parent, child, child_index, nx

    nx = solution%hierarchy%base_cells
    call fill_physical_ghosts_1d( &
      solution%levels(1)%patches(1)%state, &
      solution%levels(1)%patches(1)%temperature, nx, &
      config%boundary_condition, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    do relation = 1, size(solution%hierarchy%relations)
      do parent = 1, &
          solution%hierarchy%relations(relation)%parent_patch_count()
        do child = 1, solution%hierarchy%relations(relation)% &
            child_sets(parent)%patch_count()
          child_index = solution%hierarchy%relations(relation)% &
            child_index(parent, child)
          call fill_one_child_ghosts( &
            species, config, solution%hierarchy%relations(relation)% &
              child_sets(parent)%patches(child), &
            solution%levels(relation)%patches(parent)%state, &
            solution%levels(relation)%patches(parent)%state, 1.0_dp, &
            solution%levels(relation + 1)%patches(child_index), local_ok)
          if (.not. local_ok) then
            ok = .false.
            return
          end if
        end do
      end do
    end do
    ok = .true.
  end subroutine refresh_patch_tree_ghosts

  subroutine fill_one_child_ghosts( &
      species, config, geometry, parent_start, parent_end, alpha, child, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_two_level_hierarchy_1d), intent(in) :: geometry
    real(dp), intent(in) :: parent_start(:, 0:), parent_end(:, 0:)
    real(dp), intent(in) :: alpha
    type(amr_patch_tree_reactive_patch_1d), intent(inout) :: child
    logical, intent(out) :: ok

    logical :: local_ok

    call fill_fine_ghosts_1d( &
      species, geometry, parent_start, parent_end, alpha, child%state, &
      child%temperature, local_ok, config%boundary_condition)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    if (uses_ppm_reconstruction(config)) then
      call fill_fine_wide_ghosts_1d( &
        species, geometry, parent_start, parent_end, alpha, &
        child%left_ghost_state, child%right_ghost_state, &
        child%left_ghost_temperature, child%right_ghost_temperature, &
        local_ok, child%state, child%temperature, config%boundary_condition)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
    end if
    ok = .true.
  end subroutine fill_one_child_ghosts

  pure logical function patch_tree_children_are_interior(hierarchy) &
      result(interior)
    type(amr_patch_tree_hierarchy_1d), intent(in) :: hierarchy

    integer :: relation, parent, child

    interior = hierarchy%is_valid()
    if (.not. interior) return
    do relation = 1, size(hierarchy%relations)
      do parent = 1, hierarchy%relations(relation)%parent_patch_count()
        do child = 1, hierarchy%relations(relation)% &
            child_sets(parent)%patch_count()
          interior = .not. hierarchy%relations(relation)% &
            child_sets(parent)%patches(child)%touches_left_boundary() .and. &
            .not. hierarchy%relations(relation)% &
            child_sets(parent)%patches(child)%touches_right_boundary()
          if (.not. interior) return
        end do
      end do
    end do
  end function patch_tree_children_are_interior

  pure logical function uses_ppm_reconstruction(config) result(enabled)
    type(reactive_1d_config), intent(in) :: config

    enabled = trim(config%amr_reconstruction) == "ppm" .or. &
      trim(config%amr_reconstruction) == "characteristic_ppm"
  end function uses_ppm_reconstruction

  pure function patch_boundary(config, level) result(boundary)
    type(reactive_1d_config), intent(in) :: config
    integer, intent(in) :: level
    character(len=32) :: boundary

    boundary = "coarse_fine"
    if (level == 1) boundary = config%boundary_condition
  end function patch_boundary

  pure function chemistry_boundary(config, level) result(boundary)
    type(reactive_1d_config), intent(in) :: config
    integer, intent(in) :: level
    character(len=32) :: boundary

    boundary = "outflow"
    if (level == 1) boundary = config%boundary_condition
  end function chemistry_boundary

end module amr_patch_tree_reactive_1d_mod
