module amr_multilevel_reactive_1d_mod
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: gas_transport_species
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_cfl_timestep, reactive_transport_timestep, &
    initialize_reactive_1d, advance_reactive_chemistry
  use amr_hierarchy_1d_mod, only: &
    amr_multilevel_hierarchy_1d, amr_level_field_1d, &
    amr_flux_register_1d, initialize_multilevel_hierarchy_1d, &
    prolong_multilevel_1d, initialize_flux_register_1d, &
    accumulate_coarse_flux_1d, accumulate_fine_flux_1d, reflux_1d, &
    average_down_1d, composite_integral_multilevel_1d
  use amr_reactive_1d_mod, only: &
    advance_amr_level_1d, advance_transport_level_1d, &
    recover_level_temperatures_1d, fill_physical_ghosts_1d, &
    fill_fine_ghosts_1d
  implicit none
  private

  type, public :: amr_reactive_level_1d
    real(dp), allocatable :: state(:, :)
    real(dp), allocatable :: temperature(:)
  end type amr_reactive_level_1d

  type, public :: amr_multilevel_reactive_solution_1d
    type(amr_multilevel_hierarchy_1d) :: hierarchy
    type(amr_reactive_level_1d), allocatable :: levels(:)
    real(dp) :: time = 0.0_dp
    integer :: steps = 0
  contains
    procedure :: level_count => multilevel_reactive_level_count
    procedure :: is_valid => multilevel_reactive_is_valid
  end type amr_multilevel_reactive_solution_1d

  public :: initialize_multilevel_reactive_1d
  public :: multilevel_reactive_timestep_1d
  public :: advance_multilevel_reactive_1d
  public :: multilevel_reactive_integrals_1d

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
        allocated(self%levels(level)%temperature)
      if (.not. valid) return
      nx = self%hierarchy%level_cell_count(level - 1)
      valid = lbound(self%levels(level)%state, 2) == 0 .and. &
        ubound(self%levels(level)%state, 2) == nx + 1 .and. &
        size(self%levels(level)%state, 1) == nvar .and. &
        lbound(self%levels(level)%temperature, 1) == 0 .and. &
        ubound(self%levels(level)%temperature, 1) == nx + 1
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
      solution%levels(level)%state = 0.0_dp
      solution%levels(level)%temperature = 0.0_dp
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
    ok = solution%is_valid() .and. root_dx > 0.0_dp
  end subroutine initialize_multilevel_reactive_1d

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
    call advance_amr_level_1d( &
      species, solution%levels(level)%state, &
      solution%levels(level)%temperature, nx, dx, interval, &
      config%amr_reconstruction, config%limiter, config%riemann_solver, &
      physical_boundary, level_boundary(config, level), flux, local_ok)
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
      if (trim(config%amr_reconstruction) == "plm") then
        alpha = (real(substep, dp) - 0.5_dp) / real(ratio, dp)
      else
        alpha = real(substep - 1, dp) / real(ratio, dp)
      end if
      call fill_fine_ghosts_1d( &
        species, solution%hierarchy%interfaces(level), state_start, &
        state_end, alpha, solution%levels(level + 1)%state, &
        solution%levels(level + 1)%temperature, local_ok)
      if (.not. local_ok) return
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
    do level = 2, solution%level_count()
      call fill_fine_ghosts_1d( &
        species, solution%hierarchy%interfaces(level - 1), &
        solution%levels(level - 1)%state, &
        solution%levels(level - 1)%state, 1.0_dp, &
        solution%levels(level)%state, &
        solution%levels(level)%temperature, local_ok)
      if (.not. local_ok) return
    end do
    ok = .true.
  end subroutine refresh_multilevel_ghosts

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
