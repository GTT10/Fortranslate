module mpi_amr_sparse_patch_1d_mod
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: gas_transport_species
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: advance_reactive_chemistry
  use amr_hierarchy_1d_mod, only: &
    amr_level_field_1d, accumulate_coarse_flux_1d, &
    accumulate_fine_flux_1d
  use amr_multipatch_1d_mod, only: &
    average_down_patch_set_1d, synchronize_patch_set_1d
  use amr_reactive_1d_mod, only: &
    recover_level_temperatures_1d, fill_physical_ghosts_1d, &
    advance_amr_level_1d, advance_transport_level_1d
  use amr_patch_tree_1d_mod, only: &
    amr_patch_tree_hierarchy_1d, &
    amr_patch_tree_relation_flux_registers_1d, &
    initialize_patch_tree_flux_registers_1d
  use amr_patch_tree_reactive_1d_mod, only: &
    amr_patch_tree_reactive_patch_1d, &
    amr_patch_tree_reactive_solution_1d, fill_one_child_ghosts
  use mpi_amr_patch_1d_mod, only: &
    mpi_amr_patch_distribution_1d, &
    mpi_amr_distribution_matches_hierarchy_1d, &
    synchronize_owned_patch_tree_reactive_1d
  implicit none
  private

  type, public :: mpi_amr_sparse_reactive_level_1d
    type(amr_patch_tree_reactive_patch_1d), allocatable :: patches(:)
    logical, allocatable :: is_local(:)
  end type mpi_amr_sparse_reactive_level_1d

  type, public :: mpi_amr_sparse_reactive_solution_1d
    type(amr_patch_tree_hierarchy_1d) :: hierarchy
    type(mpi_amr_sparse_reactive_level_1d), allocatable :: levels(:)
    integer, allocatable :: level_advances(:)
    integer, allocatable :: transport_level_advances(:)
    integer :: rank = -1
    integer :: nranks = 0
    integer :: nvar = 0
    integer :: ghost_width = 0
    real(dp) :: time = 0.0_dp
    integer :: steps = 0
    integer :: regrid_evaluations = 0
    integer :: regrids = 0
    integer :: overlap_cells_transferred = 0
  contains
    procedure :: is_valid => mpi_amr_sparse_reactive_is_valid
    procedure :: local_patch_count => mpi_amr_sparse_local_patch_count
    procedure :: local_cell_count => mpi_amr_sparse_local_cell_count
    procedure :: local_value_count => mpi_amr_sparse_local_value_count
  end type mpi_amr_sparse_reactive_solution_1d

  public :: scatter_owned_patch_tree_reactive_1d
  public :: gather_owned_patch_tree_reactive_1d
  public :: migrate_owned_patch_tree_reactive_1d
  public :: advance_sparse_patch_tree_chemistry_1d
  public :: advance_sparse_patch_tree_hydro_1d
  public :: advance_sparse_patch_tree_transport_1d

contains

  subroutine advance_sparse_patch_tree_chemistry_1d( &
      species, reactions, config, interval, distribution, solution, ok, &
      local_patch_advances)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: interval
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_patch_advances

    type(mpi_amr_sparse_reactive_solution_1d) :: backup
    character(len=32) :: boundary
    logical :: local_ok, accepted, mpi_ok
    integer :: level, patch, nx, advances

    ok = .false.
    advances = 0
    if (present(local_patch_advances)) local_patch_advances = 0
    local_ok = interval >= 0.0_dp .and. size(species) >= 1 .and. &
      solution%is_valid(distribution)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
    backup = solution

    do level = 1, size(solution%levels)
      boundary = "outflow"
      if (level == 1) boundary = config%boundary_condition
      do patch = 1, size(solution%levels(level)%patches)
        local_ok = .true.
        if (solution%levels(level)%is_local(patch)) then
          nx = size(solution%levels(level)%patches(patch)%state, 2) - 2
          call advance_reactive_chemistry( &
            species, reactions, &
            solution%levels(level)%patches(patch)%state, &
            solution%levels(level)%patches(patch)%temperature, nx, interval, &
            config%chemistry_relative_tolerance, &
            config%chemistry_absolute_tolerance, boundary, local_ok)
          if (local_ok) advances = advances + 1
        end if
        call all_ranks_accept_sparse_1d( &
          distribution, local_ok, accepted, mpi_ok)
        if (.not. mpi_ok .or. .not. accepted) then
          solution = backup
          advances = 0
          return
        end if
      end do
    end do

    call average_down_sparse_reactive_solution_1d( &
      species, distribution, solution, local_ok)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) then
      solution = backup
      advances = 0
      return
    end if
    call refresh_sparse_reactive_ghosts_1d( &
      species, config, distribution, solution, local_ok)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) then
      solution = backup
      advances = 0
      return
    end if
    local_ok = solution%is_valid(distribution)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    ok = mpi_ok .and. accepted
    if (.not. ok) then
      solution = backup
      advances = 0
      return
    end if
    if (present(local_patch_advances)) local_patch_advances = advances
  end subroutine advance_sparse_patch_tree_chemistry_1d

  subroutine advance_sparse_patch_tree_hydro_1d( &
      species, config, dt, distribution, solution, ok, &
      local_patch_advances)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: dt
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_patch_advances

    type(mpi_amr_sparse_reactive_solution_1d) :: backup
    type(amr_patch_tree_relation_flux_registers_1d), allocatable :: registers(:)
    real(dp), allocatable :: left_integral(:), right_integral(:)
    logical :: local_ok, accepted, mpi_ok
    integer :: advances

    ok = .false.
    advances = 0
    if (present(local_patch_advances)) local_patch_advances = 0
    local_ok = dt > 0.0_dp .and. size(species) >= 1 .and. &
      solution%is_valid(distribution)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
    backup = solution

    allocate(left_integral(solution%nvar), right_integral(solution%nvar))
    call initialize_patch_tree_flux_registers_1d( &
      solution%hierarchy, solution%nvar, registers, local_ok)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) then
      solution = backup
      return
    end if
    call advance_sparse_patch_hydro_recursive_1d( &
      species, config, distribution, solution, registers, 1, 1, dt, &
      left_integral, right_integral, advances, local_ok)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) then
      solution = backup
      advances = 0
      return
    end if
    call refresh_sparse_reactive_ghosts_1d( &
      species, config, distribution, solution, local_ok)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) then
      solution = backup
      advances = 0
      return
    end if
    solution%time = solution%time + dt
    solution%steps = solution%steps + 1
    local_ok = solution%is_valid(distribution)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    ok = mpi_ok .and. accepted
    if (.not. ok) then
      solution = backup
      advances = 0
      return
    end if
    if (present(local_patch_advances)) local_patch_advances = advances
  end subroutine advance_sparse_patch_tree_hydro_1d

  recursive subroutine advance_sparse_patch_hydro_recursive_1d( &
      species, config, distribution, solution, registers, level, &
      parent_patch, interval, left_integral, right_integral, &
      local_advances, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    type(amr_patch_tree_relation_flux_registers_1d), &
      intent(inout) :: registers(:)
    integer, intent(in) :: level, parent_patch
    real(dp), intent(in) :: interval
    real(dp), intent(out) :: left_integral(:), right_integral(:)
    integer, intent(inout) :: local_advances
    logical, intent(out) :: ok

    real(dp), allocatable :: state_start(:, :), state_end(:, :), flux(:, :)
    real(dp), allocatable :: child_left(:, :), child_right(:, :)
    real(dp) :: child_interval, alpha, dx
    logical :: local_ok, patch_ok, accepted, mpi_ok, physical_boundary
    integer :: nx, owner, ratio, substep, child, child_count, global_child
    integer :: ierr
    character(len=32) :: boundary

    ok = .false.
    left_integral = 0.0_dp
    right_integral = 0.0_dp
    if (interval <= 0.0_dp .or. level < 1 .or. &
        level > size(solution%levels)) return
    if (parent_patch < 1 .or. parent_patch > &
        size(solution%levels(level)%patches)) return
    if (size(left_integral) /= solution%nvar .or. &
        size(right_integral) /= solution%nvar) return
    nx = distribution%levels(level)%cell_counts(parent_patch)
    owner = distribution%owner_of(level - 1, parent_patch)
    allocate(state_start(solution%nvar, 0:nx + 1))
    allocate(state_end(solution%nvar, 0:nx + 1))
    allocate(flux(solution%nvar, 0:nx))
    patch_ok = .true.
    if (distribution%rank == owner) then
      state_start = solution%levels(level)%patches(parent_patch)%state
      dx = solution%hierarchy%level_dx(level - 1)
      physical_boundary = level == 1
      boundary = "outflow"
      if (physical_boundary) boundary = config%boundary_condition
      if (physical_boundary .or. .not. sparse_uses_wide_ghosts(config)) then
        call advance_amr_level_1d( &
          species, solution%levels(level)%patches(parent_patch)%state, &
          solution%levels(level)%patches(parent_patch)%temperature, nx, dx, &
          interval, config%amr_reconstruction, config%limiter, &
          config%riemann_solver, physical_boundary, boundary, flux, &
          patch_ok, ppm_contact_steepening=config%ppm_contact_steepening, &
          ppm_shock_flattening=config%ppm_shock_flattening, &
          amr_hybrid_weno=config%amr_hybrid_weno, &
          amr_weno_scheme=config%amr_weno_scheme)
      else
        call advance_amr_level_1d( &
          species, solution%levels(level)%patches(parent_patch)%state, &
          solution%levels(level)%patches(parent_patch)%temperature, nx, dx, &
          interval, config%amr_reconstruction, config%limiter, &
          config%riemann_solver, physical_boundary, boundary, flux, &
          patch_ok, &
          solution%levels(level)%patches(parent_patch)%left_ghost_state, &
          solution%levels(level)%patches(parent_patch)%right_ghost_state, &
          solution%levels(level)%patches(parent_patch)% &
            left_ghost_temperature, &
          solution%levels(level)%patches(parent_patch)% &
            right_ghost_temperature, &
          config%ppm_contact_steepening, config%ppm_shock_flattening, &
          config%amr_hybrid_weno, config%amr_weno_scheme)
      end if
      if (patch_ok) then
        state_end = solution%levels(level)%patches(parent_patch)%state
        solution%level_advances(level) = &
          solution%level_advances(level) + 1
        local_advances = local_advances + 1
      end if
    end if
    call all_ranks_accept_sparse_1d( &
      distribution, patch_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
    call MPI_Bcast(state_start, size(state_start), MPI_DOUBLE_PRECISION, &
      owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Bcast(state_end, size(state_end), MPI_DOUBLE_PRECISION, &
      owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Bcast(flux, size(flux), MPI_DOUBLE_PRECISION, owner, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Bcast(solution%level_advances(level), 1, MPI_INTEGER, owner, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
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
    allocate(child_left(solution%nvar, child_count))
    allocate(child_right(solution%nvar, child_count))
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
      local_ok = .true.
      do child = 1, child_count
        global_child = solution%hierarchy%relations(level)% &
          child_index(parent_patch, child)
        if (.not. solution%levels(level + 1)%is_local(global_child)) cycle
        call fill_one_child_ghosts( &
          species, config, solution%hierarchy%relations(level)% &
            child_sets(parent_patch)%patches(child), state_start, state_end, &
          alpha, solution%levels(level + 1)%patches(global_child), local_ok)
        if (.not. local_ok) exit
      end do
      call all_ranks_accept_sparse_1d( &
        distribution, local_ok, accepted, mpi_ok)
      if (.not. mpi_ok .or. .not. accepted) return
      call exchange_sparse_adjacent_child_ghosts_1d( &
        config, distribution, solution, level, parent_patch, local_ok)
      call all_ranks_accept_sparse_1d( &
        distribution, local_ok, accepted, mpi_ok)
      if (.not. mpi_ok .or. .not. accepted) return
      do child = 1, child_count
        global_child = solution%hierarchy%relations(level)% &
          child_index(parent_patch, child)
        call advance_sparse_patch_hydro_recursive_1d( &
          species, config, distribution, solution, registers, level + 1, &
          global_child, child_interval, child_left(:, child), &
          child_right(:, child), local_advances, local_ok)
        if (.not. local_ok) return
      end do
      call reconcile_sparse_adjacent_child_fluxes_1d( &
        species, distribution, solution, level, parent_patch, child_left, &
        child_right, local_ok)
      call all_ranks_accept_sparse_1d( &
        distribution, local_ok, accepted, mpi_ok)
      if (.not. mpi_ok .or. .not. accepted) return
      do child = 1, child_count
        call accumulate_fine_flux_1d( &
          registers(level)%parents(parent_patch)%children(child), &
          child_left(:, child) / child_interval, &
          child_right(:, child) / child_interval, child_interval, local_ok)
        if (.not. local_ok) return
      end do
    end do

    call synchronize_sparse_parent_1d( &
      species, distribution, solution, registers, level, parent_patch, &
      local_ok)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    ok = mpi_ok .and. accepted
  end subroutine advance_sparse_patch_hydro_recursive_1d

  subroutine advance_sparse_patch_tree_transport_1d( &
      species, transport, config, interval, distribution, solution, ok, &
      local_patch_advances)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: interval
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_patch_advances

    type(mpi_amr_sparse_reactive_solution_1d) :: backup
    type(amr_patch_tree_relation_flux_registers_1d), allocatable :: registers(:)
    real(dp), allocatable :: left_integral(:), right_integral(:)
    logical :: local_ok, accepted, mpi_ok
    integer :: advances

    ok = .false.
    advances = 0
    if (present(local_patch_advances)) local_patch_advances = 0
    local_ok = interval > 0.0_dp .and. size(species) >= 1 .and. &
      size(transport) == size(species) .and. &
      solution%is_valid(distribution)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
    backup = solution

    allocate(left_integral(solution%nvar), right_integral(solution%nvar))
    call initialize_patch_tree_flux_registers_1d( &
      solution%hierarchy, solution%nvar, registers, local_ok)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) then
      solution = backup
      return
    end if
    call advance_sparse_patch_transport_recursive_1d( &
      species, transport, config, distribution, solution, registers, 1, 1, &
      interval, left_integral, right_integral, advances, local_ok)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) then
      solution = backup
      advances = 0
      return
    end if
    call refresh_sparse_reactive_ghosts_1d( &
      species, config, distribution, solution, local_ok)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) then
      solution = backup
      advances = 0
      return
    end if
    local_ok = solution%is_valid(distribution)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    ok = mpi_ok .and. accepted
    if (.not. ok) then
      solution = backup
      advances = 0
      return
    end if
    if (present(local_patch_advances)) local_patch_advances = advances
  end subroutine advance_sparse_patch_tree_transport_1d

  recursive subroutine advance_sparse_patch_transport_recursive_1d( &
      species, transport, config, distribution, solution, registers, level, &
      parent_patch, interval, left_integral, right_integral, &
      local_advances, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_1d_config), intent(in) :: config
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    type(amr_patch_tree_relation_flux_registers_1d), &
      intent(inout) :: registers(:)
    integer, intent(in) :: level, parent_patch
    real(dp), intent(in) :: interval
    real(dp), intent(out) :: left_integral(:), right_integral(:)
    integer, intent(inout) :: local_advances
    logical, intent(out) :: ok

    real(dp), allocatable :: state_start(:, :), state_end(:, :), flux(:, :)
    real(dp), allocatable :: child_left(:, :), child_right(:, :)
    real(dp) :: child_interval, alpha, dx, boundary_distance
    logical :: local_ok, patch_ok, accepted, mpi_ok, physical_boundary
    integer :: nx, owner, ratio, subcycles, substep, child, child_count
    integer :: global_child, ierr
    character(len=32) :: boundary

    ok = .false.
    left_integral = 0.0_dp
    right_integral = 0.0_dp
    if (interval <= 0.0_dp .or. level < 1 .or. &
        level > size(solution%levels)) return
    if (parent_patch < 1 .or. parent_patch > &
        size(solution%levels(level)%patches)) return
    if (size(left_integral) /= solution%nvar .or. &
        size(right_integral) /= solution%nvar) return
    nx = distribution%levels(level)%cell_counts(parent_patch)
    owner = distribution%owner_of(level - 1, parent_patch)
    allocate(state_start(solution%nvar, 0:nx + 1))
    allocate(state_end(solution%nvar, 0:nx + 1))
    allocate(flux(solution%nvar, 0:nx))
    patch_ok = .true.
    if (distribution%rank == owner) then
      state_start = solution%levels(level)%patches(parent_patch)%state
      dx = solution%hierarchy%level_dx(level - 1)
      physical_boundary = level == 1
      boundary_distance = dx
      if (.not. physical_boundary) boundary_distance = 0.5_dp * ( &
        solution%hierarchy%level_dx(level - 2) + dx)
      boundary = "outflow"
      if (physical_boundary) boundary = config%boundary_condition
      call advance_transport_level_1d( &
        species, transport, &
        solution%levels(level)%patches(parent_patch)%state, &
        solution%levels(level)%patches(parent_patch)%temperature, nx, dx, &
        interval, boundary_distance, config, physical_boundary, boundary, &
        flux, patch_ok)
      if (patch_ok) then
        state_end = solution%levels(level)%patches(parent_patch)%state
        solution%transport_level_advances(level) = &
          solution%transport_level_advances(level) + 1
        local_advances = local_advances + 1
      end if
    end if
    call all_ranks_accept_sparse_1d( &
      distribution, patch_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
    call MPI_Bcast(state_start, size(state_start), MPI_DOUBLE_PRECISION, &
      owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Bcast(state_end, size(state_end), MPI_DOUBLE_PRECISION, &
      owner, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Bcast(flux, size(flux), MPI_DOUBLE_PRECISION, owner, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Bcast( &
      solution%transport_level_advances(level), 1, MPI_INTEGER, owner, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
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
    allocate(child_left(solution%nvar, child_count))
    allocate(child_right(solution%nvar, child_count))
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
      local_ok = .true.
      do child = 1, child_count
        global_child = solution%hierarchy%relations(level)% &
          child_index(parent_patch, child)
        if (.not. solution%levels(level + 1)%is_local(global_child)) cycle
        call fill_one_child_ghosts( &
          species, config, solution%hierarchy%relations(level)% &
            child_sets(parent_patch)%patches(child), state_start, state_end, &
          alpha, solution%levels(level + 1)%patches(global_child), local_ok)
        if (.not. local_ok) exit
      end do
      call all_ranks_accept_sparse_1d( &
        distribution, local_ok, accepted, mpi_ok)
      if (.not. mpi_ok .or. .not. accepted) return
      call exchange_sparse_adjacent_child_ghosts_1d( &
        config, distribution, solution, level, parent_patch, local_ok)
      call all_ranks_accept_sparse_1d( &
        distribution, local_ok, accepted, mpi_ok)
      if (.not. mpi_ok .or. .not. accepted) return
      do child = 1, child_count
        global_child = solution%hierarchy%relations(level)% &
          child_index(parent_patch, child)
        call advance_sparse_patch_transport_recursive_1d( &
          species, transport, config, distribution, solution, registers, &
          level + 1, global_child, child_interval, child_left(:, child), &
          child_right(:, child), local_advances, local_ok)
        if (.not. local_ok) return
      end do
      call reconcile_sparse_adjacent_child_fluxes_1d( &
        species, distribution, solution, level, parent_patch, child_left, &
        child_right, local_ok)
      call all_ranks_accept_sparse_1d( &
        distribution, local_ok, accepted, mpi_ok)
      if (.not. mpi_ok .or. .not. accepted) return
      do child = 1, child_count
        call accumulate_fine_flux_1d( &
          registers(level)%parents(parent_patch)%children(child), &
          child_left(:, child) / child_interval, &
          child_right(:, child) / child_interval, child_interval, local_ok)
        if (.not. local_ok) return
      end do
    end do

    call synchronize_sparse_parent_1d( &
      species, distribution, solution, registers, level, parent_patch, &
      local_ok)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    ok = mpi_ok .and. accepted
  end subroutine advance_sparse_patch_transport_recursive_1d

  subroutine synchronize_sparse_parent_1d( &
      species, distribution, solution, registers, level, parent_patch, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    type(amr_patch_tree_relation_flux_registers_1d), &
      intent(inout) :: registers(:)
    integer, intent(in) :: level, parent_patch
    logical, intent(out) :: ok

    type(amr_level_field_1d), allocatable :: child_fields(:)
    real(dp), allocatable :: child_state(:, :)
    logical :: local_ok
    integer :: child, child_count, child_index, child_owner, child_nx
    integer :: parent_owner, parent_nx, ierr

    ok = .false.
    child_count = solution%hierarchy%relations(level)% &
      child_sets(parent_patch)%patch_count()
    parent_owner = distribution%owner_of(level - 1, parent_patch)
    allocate(child_fields(child_count))
    do child = 1, child_count
      child_index = solution%hierarchy%relations(level)% &
        child_index(parent_patch, child)
      child_owner = distribution%owner_of(level, child_index)
      child_nx = distribution%levels(level + 1)%cell_counts(child_index)
      allocate(child_state(solution%nvar, child_nx))
      if (distribution%rank == child_owner) &
        child_state = solution%levels(level + 1)% &
          patches(child_index)%state(:, 1:child_nx)
      call MPI_Bcast(child_state, size(child_state), MPI_DOUBLE_PRECISION, &
        child_owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
      if (distribution%rank == parent_owner) &
        child_fields(child)%values = child_state
      deallocate(child_state)
    end do
    local_ok = .true.
    if (distribution%rank == parent_owner) then
      parent_nx = distribution%levels(level)%cell_counts(parent_patch)
      call synchronize_patch_set_1d( &
        solution%levels(level)%patches(parent_patch)%state(:, 1:parent_nx), &
        child_fields, solution%hierarchy%relations(level)% &
          child_sets(parent_patch), &
        registers(level)%parents(parent_patch)%children, local_ok)
      if (local_ok) call recover_level_temperatures_1d( &
        species, solution%levels(level)%patches(parent_patch)%state, &
        solution%levels(level)%patches(parent_patch)%temperature, &
        parent_nx, local_ok)
    end if
    ok = local_ok
  end subroutine synchronize_sparse_parent_1d

  subroutine reconcile_sparse_adjacent_child_fluxes_1d( &
      species, distribution, solution, relation, parent, left_integrals, &
      right_integrals, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    integer, intent(in) :: relation, parent
    real(dp), intent(inout) :: left_integrals(:, :), right_integrals(:, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: shared_integral(:)
    logical, allocatable :: touched(:)
    logical :: local_ok
    real(dp) :: dx
    integer :: child, child_count, left_index, right_index, left_nx

    ok = .false.
    child_count = solution%hierarchy%relations(relation)% &
      child_sets(parent)%patch_count()
    if (size(left_integrals, 1) /= solution%nvar .or. &
        size(right_integrals, 1) /= solution%nvar .or. &
        size(left_integrals, 2) /= child_count .or. &
        size(right_integrals, 2) /= child_count) return
    allocate(shared_integral(solution%nvar), touched(child_count))
    touched = .false.
    dx = solution%hierarchy%level_dx(relation)
    if (dx <= 0.0_dp) return
    do child = 1, child_count - 1
      if (solution%hierarchy%relations(relation)%child_sets(parent)% &
            patches(child)%fine_coarse_upper + 1 /= &
          solution%hierarchy%relations(relation)%child_sets(parent)% &
            patches(child + 1)%fine_coarse_lower) cycle
      left_index = solution%hierarchy%relations(relation)% &
        child_index(parent, child)
      right_index = solution%hierarchy%relations(relation)% &
        child_index(parent, child + 1)
      shared_integral = 0.5_dp * ( &
        right_integrals(:, child) + left_integrals(:, child + 1))
      if (solution%levels(relation + 1)%is_local(left_index)) then
        left_nx = distribution%levels(relation + 1)%cell_counts(left_index)
        solution%levels(relation + 1)%patches(left_index)% &
          state(:, left_nx) = &
          solution%levels(relation + 1)%patches(left_index)% &
            state(:, left_nx) - &
          (shared_integral - right_integrals(:, child)) / dx
      end if
      if (solution%levels(relation + 1)%is_local(right_index)) &
        solution%levels(relation + 1)%patches(right_index)%state(:, 1) = &
          solution%levels(relation + 1)%patches(right_index)%state(:, 1) + &
          (shared_integral - left_integrals(:, child + 1)) / dx
      right_integrals(:, child) = shared_integral
      left_integrals(:, child + 1) = shared_integral
      touched(child) = .true.
      touched(child + 1) = .true.
    end do
    local_ok = .true.
    do child = 1, child_count
      if (.not. touched(child)) cycle
      left_index = solution%hierarchy%relations(relation)% &
        child_index(parent, child)
      if (.not. solution%levels(relation + 1)%is_local(left_index)) cycle
      left_nx = distribution%levels(relation + 1)%cell_counts(left_index)
      call recover_level_temperatures_1d( &
        species, solution%levels(relation + 1)%patches(left_index)%state, &
        solution%levels(relation + 1)%patches(left_index)%temperature, &
        left_nx, local_ok)
      if (.not. local_ok) return
    end do
    ok = .true.
  end subroutine reconcile_sparse_adjacent_child_fluxes_1d

  subroutine average_down_sparse_reactive_solution_1d( &
      species, distribution, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok

    type(amr_level_field_1d), allocatable :: children(:)
    real(dp), allocatable :: child_state(:, :)
    logical :: local_ok, accepted, mpi_ok
    integer :: relation, parent, child, child_index, child_count
    integer :: parent_owner, child_owner, parent_nx, child_nx, ierr

    ok = .false.
    do relation = size(solution%hierarchy%relations), 1, -1
      do parent = 1, solution%hierarchy%relations(relation)% &
          parent_patch_count()
        child_count = solution%hierarchy%relations(relation)% &
          child_sets(parent)%patch_count()
        if (child_count == 0) cycle
        parent_owner = distribution%owner_of(relation - 1, parent)
        allocate(children(child_count))
        do child = 1, child_count
          child_index = solution%hierarchy%relations(relation)% &
            child_index(parent, child)
          child_nx = distribution%levels(relation + 1)% &
            cell_counts(child_index)
          allocate(child_state(solution%nvar, child_nx))
          child_owner = distribution%owner_of(relation, child_index)
          if (distribution%rank == child_owner) &
            child_state = solution%levels(relation + 1)% &
              patches(child_index)%state(:, 1:child_nx)
          call MPI_Bcast(child_state, size(child_state), &
            MPI_DOUBLE_PRECISION, child_owner, distribution%comm, ierr)
          if (ierr /= MPI_SUCCESS) return
          if (distribution%rank == parent_owner) &
            children(child)%values = child_state
          deallocate(child_state)
        end do
        local_ok = .true.
        if (distribution%rank == parent_owner) then
          parent_nx = distribution%levels(relation)%cell_counts(parent)
          call average_down_patch_set_1d( &
            children, solution%hierarchy%relations(relation)% &
              child_sets(parent), &
            solution%levels(relation)%patches(parent)% &
              state(:, 1:parent_nx), local_ok)
          if (local_ok) call recover_level_temperatures_1d( &
            species, solution%levels(relation)%patches(parent)%state, &
            solution%levels(relation)%patches(parent)%temperature, &
            parent_nx, local_ok)
        end if
        deallocate(children)
        call all_ranks_accept_sparse_1d( &
          distribution, local_ok, accepted, mpi_ok)
        if (.not. mpi_ok .or. .not. accepted) return
      end do
    end do
    ok = .true.
  end subroutine average_down_sparse_reactive_solution_1d

  subroutine refresh_sparse_reactive_ghosts_1d( &
      species, config, distribution, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok

    real(dp), allocatable :: parent_state(:, :)
    logical :: local_ok, accepted, mpi_ok
    integer :: relation, parent, child, child_index, child_count
    integer :: parent_owner, parent_nx, ierr

    ok = .false.
    local_ok = .true.
    parent_owner = distribution%owner_of(0, 1)
    if (distribution%rank == parent_owner) then
      parent_nx = distribution%levels(1)%cell_counts(1)
      call fill_physical_ghosts_1d( &
        solution%levels(1)%patches(1)%state, &
        solution%levels(1)%patches(1)%temperature, parent_nx, &
        config%boundary_condition, local_ok)
    end if
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return

    do relation = 1, size(solution%hierarchy%relations)
      do parent = 1, solution%hierarchy%relations(relation)% &
          parent_patch_count()
        child_count = solution%hierarchy%relations(relation)% &
          child_sets(parent)%patch_count()
        if (child_count == 0) cycle
        parent_owner = distribution%owner_of(relation - 1, parent)
        parent_nx = distribution%levels(relation)%cell_counts(parent)
        allocate(parent_state(solution%nvar, 0:parent_nx + 1))
        if (distribution%rank == parent_owner) &
          parent_state = solution%levels(relation)%patches(parent)%state
        call MPI_Bcast(parent_state, size(parent_state), &
          MPI_DOUBLE_PRECISION, parent_owner, distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS) return
        local_ok = .true.
        do child = 1, child_count
          child_index = solution%hierarchy%relations(relation)% &
            child_index(parent, child)
          if (.not. solution%levels(relation + 1)% &
              is_local(child_index)) cycle
          call fill_one_child_ghosts( &
            species, config, solution%hierarchy%relations(relation)% &
              child_sets(parent)%patches(child), &
            parent_state, parent_state, 1.0_dp, &
            solution%levels(relation + 1)%patches(child_index), local_ok)
          if (.not. local_ok) exit
        end do
        deallocate(parent_state)
        call all_ranks_accept_sparse_1d( &
          distribution, local_ok, accepted, mpi_ok)
        if (.not. mpi_ok .or. .not. accepted) return
        call exchange_sparse_adjacent_child_ghosts_1d( &
          config, distribution, solution, relation, parent, local_ok)
        call all_ranks_accept_sparse_1d( &
          distribution, local_ok, accepted, mpi_ok)
        if (.not. mpi_ok .or. .not. accepted) return
      end do
    end do
    ok = .true.
  end subroutine refresh_sparse_reactive_ghosts_1d

  subroutine exchange_sparse_adjacent_child_ghosts_1d( &
      config, distribution, solution, relation, parent, ok)
    type(reactive_1d_config), intent(in) :: config
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    integer, intent(in) :: relation, parent
    logical, intent(out) :: ok

    real(dp), allocatable :: source_state(:, :), source_temperature(:)
    integer :: source, source_index, source_owner, source_nx, child_count
    integer :: ierr

    ok = .false.
    child_count = solution%hierarchy%relations(relation)% &
      child_sets(parent)%patch_count()
    do source = 1, child_count
      source_index = solution%hierarchy%relations(relation)% &
        child_index(parent, source)
      source_owner = distribution%owner_of(relation, source_index)
      source_nx = distribution%levels(relation + 1)% &
        cell_counts(source_index)
      allocate(source_state(solution%nvar, 0:source_nx + 1))
      allocate(source_temperature(0:source_nx + 1))
      if (distribution%rank == source_owner) then
        source_state = solution%levels(relation + 1)% &
          patches(source_index)%state
        source_temperature = solution%levels(relation + 1)% &
          patches(source_index)%temperature
      end if
      call MPI_Bcast(source_state, size(source_state), MPI_DOUBLE_PRECISION, &
        source_owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
      call MPI_Bcast(source_temperature, size(source_temperature), &
        MPI_DOUBLE_PRECISION, source_owner, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
      call apply_sparse_adjacent_source_1d( &
        config, solution, relation, parent, source, source_state, &
        source_temperature)
      deallocate(source_state, source_temperature)
    end do
    ok = .true.
  end subroutine exchange_sparse_adjacent_child_ghosts_1d

  subroutine apply_sparse_adjacent_source_1d( &
      config, solution, relation, parent, source, source_state, &
      source_temperature)
    type(reactive_1d_config), intent(in) :: config
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    integer, intent(in) :: relation, parent, source
    real(dp), intent(in) :: source_state(:, 0:), source_temperature(0:)

    integer :: target, target_index, target_nx, source_cell, global_fine
    integer :: source_lower, source_upper, layer, child_count

    child_count = solution%hierarchy%relations(relation)% &
      child_sets(parent)%patch_count()
    source_lower = solution%hierarchy%relations(relation)% &
      child_sets(parent)%patches(source)%fine%lower
    source_upper = solution%hierarchy%relations(relation)% &
      child_sets(parent)%patches(source)%fine%upper
    do target = 1, child_count
      if (target == source) cycle
      target_index = solution%hierarchy%relations(relation)% &
        child_index(parent, target)
      if (.not. solution%levels(relation + 1)%is_local(target_index)) cycle
      target_nx = size(solution%levels(relation + 1)% &
        patches(target_index)%state, 2) - 2
      global_fine = solution%hierarchy%relations(relation)% &
        child_sets(parent)%patches(target)%fine%lower - 1
      if (global_fine >= source_lower .and. global_fine <= source_upper) then
        source_cell = global_fine - source_lower + 1
        solution%levels(relation + 1)%patches(target_index)%state(:, 0) = &
          source_state(:, source_cell)
        solution%levels(relation + 1)% &
          patches(target_index)%temperature(0) = &
            source_temperature(source_cell)
      end if
      global_fine = solution%hierarchy%relations(relation)% &
        child_sets(parent)%patches(target)%fine%upper + 1
      if (global_fine >= source_lower .and. global_fine <= source_upper) then
        source_cell = global_fine - source_lower + 1
        solution%levels(relation + 1)% &
          patches(target_index)%state(:, target_nx + 1) = &
            source_state(:, source_cell)
        solution%levels(relation + 1)% &
          patches(target_index)%temperature(target_nx + 1) = &
            source_temperature(source_cell)
      end if
      if (.not. sparse_uses_wide_ghosts(config)) cycle
      do layer = 1, solution%ghost_width
        global_fine = solution%hierarchy%relations(relation)% &
          child_sets(parent)%patches(target)%fine%lower - layer
        if (global_fine >= source_lower .and. &
            global_fine <= source_upper) then
          source_cell = global_fine - source_lower + 1
          solution%levels(relation + 1)%patches(target_index)% &
            left_ghost_state(:, layer) = source_state(:, source_cell)
          solution%levels(relation + 1)%patches(target_index)% &
            left_ghost_temperature(layer) = &
              source_temperature(source_cell)
        end if
        global_fine = solution%hierarchy%relations(relation)% &
          child_sets(parent)%patches(target)%fine%upper + layer
        if (global_fine >= source_lower .and. &
            global_fine <= source_upper) then
          source_cell = global_fine - source_lower + 1
          solution%levels(relation + 1)%patches(target_index)% &
            right_ghost_state(:, layer) = source_state(:, source_cell)
          solution%levels(relation + 1)%patches(target_index)% &
            right_ghost_temperature(layer) = &
              source_temperature(source_cell)
        end if
      end do
    end do
  end subroutine apply_sparse_adjacent_source_1d

  pure logical function sparse_uses_wide_ghosts(config) result(enabled)
    type(reactive_1d_config), intent(in) :: config

    enabled = trim(config%amr_reconstruction) == "ppm" .or. &
      trim(config%amr_reconstruction) == "characteristic_ppm"
  end function sparse_uses_wide_ghosts

  logical function mpi_amr_sparse_reactive_is_valid( &
      self, distribution) result(valid)
    class(mpi_amr_sparse_reactive_solution_1d), intent(in) :: self
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution

    logical :: local
    integer :: level, patch, nx

    valid = self%rank == distribution%rank .and. &
      self%nranks == distribution%nranks .and. self%nvar >= 1 .and. &
      self%ghost_width >= 1 .and. allocated(self%levels) .and. &
      allocated(self%level_advances) .and. &
      allocated(self%transport_level_advances) .and. &
      mpi_amr_distribution_matches_hierarchy_1d( &
        distribution, self%hierarchy)
    if (.not. valid) return
    valid = size(self%levels) == self%hierarchy%level_count() .and. &
      size(self%level_advances) == size(self%levels) .and. &
      size(self%transport_level_advances) == size(self%levels) .and. &
      all(self%level_advances >= 0) .and. &
      all(self%transport_level_advances >= 0) .and. self%time >= 0.0_dp .and. &
      self%steps >= 0 .and. self%regrid_evaluations >= 0 .and. &
      self%regrids >= 0 .and. self%overlap_cells_transferred >= 0
    if (.not. valid) return
    do level = 1, size(self%levels)
      valid = allocated(self%levels(level)%patches) .and. &
        allocated(self%levels(level)%is_local)
      if (.not. valid) return
      valid = size(self%levels(level)%patches) == &
        self%hierarchy%level_patch_count(level - 1) .and. &
        size(self%levels(level)%is_local) == &
          self%hierarchy%level_patch_count(level - 1)
      if (.not. valid) return
      do patch = 1, size(self%levels(level)%patches)
        local = distribution%is_local(level - 1, patch)
        valid = self%levels(level)%is_local(patch) .eqv. local
        if (.not. valid) return
        nx = distribution%levels(level)%cell_counts(patch)
        if (local) then
          valid = sparse_patch_has_shape( &
            self%levels(level)%patches(patch), self%nvar, nx, &
            self%ghost_width)
        else
          valid = sparse_patch_is_empty(self%levels(level)%patches(patch))
        end if
        if (.not. valid) return
      end do
    end do
    valid = self%local_patch_count() == &
      distribution%rank_patch_counts(self%rank + 1) .and. &
      self%local_cell_count() == &
        distribution%rank_cell_counts(self%rank + 1)
  end function mpi_amr_sparse_reactive_is_valid

  pure integer function mpi_amr_sparse_local_patch_count(self) result(patch_count)
    class(mpi_amr_sparse_reactive_solution_1d), intent(in) :: self

    integer :: level

    patch_count = 0
    if (.not. allocated(self%levels)) return
    do level = 1, size(self%levels)
      if (allocated(self%levels(level)%is_local)) &
        patch_count = patch_count + count(self%levels(level)%is_local)
    end do
  end function mpi_amr_sparse_local_patch_count

  pure integer function mpi_amr_sparse_local_cell_count(self) result(count)
    class(mpi_amr_sparse_reactive_solution_1d), intent(in) :: self

    integer :: level, patch

    count = 0
    if (.not. allocated(self%levels)) return
    do level = 1, size(self%levels)
      if (.not. allocated(self%levels(level)%patches)) cycle
      do patch = 1, size(self%levels(level)%patches)
        if (.not. allocated(self%levels(level)%patches(patch)%state)) cycle
        count = count + size(self%levels(level)%patches(patch)%state, 2) - 2
      end do
    end do
  end function mpi_amr_sparse_local_cell_count

  pure integer function mpi_amr_sparse_local_value_count(self) result(count)
    class(mpi_amr_sparse_reactive_solution_1d), intent(in) :: self

    integer :: level, patch

    count = 0
    if (.not. allocated(self%levels)) return
    do level = 1, size(self%levels)
      if (.not. allocated(self%levels(level)%patches)) cycle
      do patch = 1, size(self%levels(level)%patches)
        if (.not. allocated(self%levels(level)%patches(patch)%state)) cycle
        count = count + sparse_patch_value_count( &
          self%levels(level)%patches(patch))
      end do
    end do
  end function mpi_amr_sparse_local_value_count

  subroutine scatter_owned_patch_tree_reactive_1d( &
      distribution, replicated, sparse, ok)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(amr_patch_tree_reactive_solution_1d), intent(inout) :: replicated
    type(mpi_amr_sparse_reactive_solution_1d), intent(out) :: sparse
    logical, intent(out) :: ok

    logical :: local_ok, accepted, mpi_ok
    integer :: level, patch

    ok = .false.
    local_ok = replicated%is_valid() .and. &
      mpi_amr_distribution_matches_hierarchy_1d( &
        distribution, replicated%hierarchy)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
    call synchronize_owned_patch_tree_reactive_1d( &
      distribution, replicated, local_ok)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return

    sparse%hierarchy = replicated%hierarchy
    sparse%rank = distribution%rank
    sparse%nranks = distribution%nranks
    sparse%nvar = size(replicated%levels(1)%patches(1)%state, 1)
    sparse%ghost_width = &
      size(replicated%levels(1)%patches(1)%left_ghost_state, 2)
    sparse%level_advances = replicated%level_advances
    sparse%transport_level_advances = replicated%transport_level_advances
    sparse%time = replicated%time
    sparse%steps = replicated%steps
    sparse%regrid_evaluations = replicated%regrid_evaluations
    sparse%regrids = replicated%regrids
    sparse%overlap_cells_transferred = replicated%overlap_cells_transferred
    allocate(sparse%levels(replicated%level_count()))
    do level = 1, size(sparse%levels)
      allocate(sparse%levels(level)%patches( &
        replicated%hierarchy%level_patch_count(level - 1)))
      allocate(sparse%levels(level)%is_local( &
        replicated%hierarchy%level_patch_count(level - 1)))
      do patch = 1, size(sparse%levels(level)%patches)
        sparse%levels(level)%is_local(patch) = &
          distribution%is_local(level - 1, patch)
        if (sparse%levels(level)%is_local(patch)) &
          sparse%levels(level)%patches(patch) = &
            replicated%levels(level)%patches(patch)
      end do
    end do
    local_ok = sparse%is_valid(distribution)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    ok = mpi_ok .and. accepted
    if (.not. ok) sparse = mpi_amr_sparse_reactive_solution_1d()
  end subroutine scatter_owned_patch_tree_reactive_1d

  subroutine gather_owned_patch_tree_reactive_1d( &
      distribution, sparse, replicated, ok)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(in) :: sparse
    type(amr_patch_tree_reactive_solution_1d), intent(inout) :: replicated
    logical, intent(out) :: ok

    logical :: local_ok, accepted, mpi_ok
    integer :: level, patch, root_owner

    ok = .false.
    local_ok = sparse%is_valid(distribution) .and. replicated%is_valid() .and. &
      mpi_amr_distribution_matches_hierarchy_1d( &
        distribution, replicated%hierarchy)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
    do level = 1, sparse%hierarchy%level_count()
      do patch = 1, sparse%hierarchy%level_patch_count(level - 1)
        if (distribution%is_local(level - 1, patch)) &
          replicated%levels(level)%patches(patch) = &
            sparse%levels(level)%patches(patch)
      end do
    end do
    root_owner = distribution%owner_of(0, 1)
    if (distribution%rank == root_owner) then
      replicated%level_advances = sparse%level_advances
      replicated%transport_level_advances = sparse%transport_level_advances
      replicated%time = sparse%time
      replicated%steps = sparse%steps
      replicated%regrid_evaluations = sparse%regrid_evaluations
      replicated%regrids = sparse%regrids
      replicated%overlap_cells_transferred = sparse%overlap_cells_transferred
    end if
    call synchronize_owned_patch_tree_reactive_1d( &
      distribution, replicated, local_ok)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    ok = mpi_ok .and. accepted
  end subroutine gather_owned_patch_tree_reactive_1d

  subroutine migrate_owned_patch_tree_reactive_1d( &
      old_distribution, new_distribution, old_sparse, new_sparse, ok)
    type(mpi_amr_patch_distribution_1d), intent(in) :: old_distribution
    type(mpi_amr_patch_distribution_1d), intent(in) :: new_distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(in) :: old_sparse
    type(mpi_amr_sparse_reactive_solution_1d), intent(out) :: new_sparse
    logical, intent(out) :: ok

    logical :: local_ok, accepted, mpi_ok
    integer :: level, patch, nx, old_owner, new_owner

    ok = .false.
    local_ok = old_sparse%is_valid(old_distribution) .and. &
      old_distribution%rank == new_distribution%rank .and. &
      old_distribution%nranks == new_distribution%nranks .and. &
      mpi_amr_distribution_matches_hierarchy_1d( &
        new_distribution, old_sparse%hierarchy)
    call all_ranks_accept_sparse_1d( &
      old_distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
    call copy_sparse_metadata(old_sparse, new_distribution, new_sparse)
    do level = 1, old_sparse%hierarchy%level_count()
      do patch = 1, old_sparse%hierarchy%level_patch_count(level - 1)
        nx = old_distribution%levels(level)%cell_counts(patch)
        old_owner = old_distribution%owner_of(level - 1, patch)
        new_owner = new_distribution%owner_of(level - 1, patch)
        call migrate_one_patch_1d( &
          old_distribution%comm, old_distribution%rank, old_owner, new_owner, &
          old_sparse%nvar, nx, old_sparse%ghost_width, &
          old_sparse%levels(level)%patches(patch), &
          new_sparse%levels(level)%patches(patch), local_ok)
        if (.not. local_ok) return
      end do
    end do
    local_ok = new_sparse%is_valid(new_distribution)
    call all_ranks_accept_sparse_1d( &
      old_distribution, local_ok, accepted, mpi_ok)
    ok = mpi_ok .and. accepted
    if (.not. ok) new_sparse = mpi_amr_sparse_reactive_solution_1d()
  end subroutine migrate_owned_patch_tree_reactive_1d

  subroutine copy_sparse_metadata(source, distribution, target)
    type(mpi_amr_sparse_reactive_solution_1d), intent(in) :: source
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(out) :: target

    integer :: level, patch_count

    target%hierarchy = source%hierarchy
    target%rank = distribution%rank
    target%nranks = distribution%nranks
    target%nvar = source%nvar
    target%ghost_width = source%ghost_width
    target%level_advances = source%level_advances
    target%transport_level_advances = source%transport_level_advances
    target%time = source%time
    target%steps = source%steps
    target%regrid_evaluations = source%regrid_evaluations
    target%regrids = source%regrids
    target%overlap_cells_transferred = source%overlap_cells_transferred
    allocate(target%levels(source%hierarchy%level_count()))
    do level = 1, size(target%levels)
      patch_count = source%hierarchy%level_patch_count(level - 1)
      allocate(target%levels(level)%patches(patch_count))
      allocate(target%levels(level)%is_local(patch_count))
      target%levels(level)%is_local = distribution%levels(level)%owners == &
        distribution%rank
    end do
  end subroutine copy_sparse_metadata

  subroutine migrate_one_patch_1d( &
      comm, rank, old_owner, new_owner, nvar, nx, ghost_width, source, &
      destination, ok)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: rank, old_owner, new_owner, nvar, nx, ghost_width
    type(amr_patch_tree_reactive_patch_1d), intent(in) :: source
    type(amr_patch_tree_reactive_patch_1d), intent(inout) :: destination
    logical, intent(out) :: ok

    type(amr_patch_tree_reactive_patch_1d) :: work
    integer :: ierr

    call allocate_sparse_patch(work, nvar, nx, ghost_width)
    if (rank == old_owner) call copy_sparse_patch_values(source, work)
    call MPI_Bcast(work%state, size(work%state), MPI_DOUBLE_PRECISION, &
      old_owner, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast(work%temperature, size(work%temperature), &
      MPI_DOUBLE_PRECISION, old_owner, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast(work%left_ghost_state, size(work%left_ghost_state), &
      MPI_DOUBLE_PRECISION, old_owner, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast(work%right_ghost_state, size(work%right_ghost_state), &
      MPI_DOUBLE_PRECISION, old_owner, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast( &
      work%left_ghost_temperature, size(work%left_ghost_temperature), &
      MPI_DOUBLE_PRECISION, old_owner, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast( &
      work%right_ghost_temperature, size(work%right_ghost_temperature), &
      MPI_DOUBLE_PRECISION, old_owner, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    if (rank == new_owner) destination = work
    ok = .true.
  end subroutine migrate_one_patch_1d

  subroutine allocate_sparse_patch(patch, nvar, nx, ghost_width)
    type(amr_patch_tree_reactive_patch_1d), intent(out) :: patch
    integer, intent(in) :: nvar, nx, ghost_width

    allocate(patch%state(nvar, 0:nx + 1))
    allocate(patch%temperature(0:nx + 1))
    allocate(patch%left_ghost_state(nvar, ghost_width))
    allocate(patch%right_ghost_state(nvar, ghost_width))
    allocate(patch%left_ghost_temperature(ghost_width))
    allocate(patch%right_ghost_temperature(ghost_width))
    patch%state = 0.0_dp
    patch%temperature = 0.0_dp
    patch%left_ghost_state = 0.0_dp
    patch%right_ghost_state = 0.0_dp
    patch%left_ghost_temperature = 0.0_dp
    patch%right_ghost_temperature = 0.0_dp
  end subroutine allocate_sparse_patch

  subroutine copy_sparse_patch_values(source, destination)
    type(amr_patch_tree_reactive_patch_1d), intent(in) :: source
    type(amr_patch_tree_reactive_patch_1d), intent(inout) :: destination

    destination%state = source%state
    destination%temperature = source%temperature
    destination%left_ghost_state = source%left_ghost_state
    destination%right_ghost_state = source%right_ghost_state
    destination%left_ghost_temperature = source%left_ghost_temperature
    destination%right_ghost_temperature = source%right_ghost_temperature
  end subroutine copy_sparse_patch_values

  pure logical function sparse_patch_has_shape( &
      patch, nvar, nx, ghost_width) result(valid)
    type(amr_patch_tree_reactive_patch_1d), intent(in) :: patch
    integer, intent(in) :: nvar, nx, ghost_width

    valid = allocated(patch%state) .and. allocated(patch%temperature) .and. &
      allocated(patch%left_ghost_state) .and. &
      allocated(patch%right_ghost_state) .and. &
      allocated(patch%left_ghost_temperature) .and. &
      allocated(patch%right_ghost_temperature)
    if (.not. valid) return
    valid = size(patch%state, 1) == nvar .and. &
      lbound(patch%state, 2) == 0 .and. &
      ubound(patch%state, 2) == nx + 1 .and. &
      lbound(patch%temperature, 1) == 0 .and. &
      ubound(patch%temperature, 1) == nx + 1 .and. &
      size(patch%left_ghost_state, 1) == nvar .and. &
      size(patch%right_ghost_state, 1) == nvar .and. &
      size(patch%left_ghost_state, 2) == ghost_width .and. &
      size(patch%right_ghost_state, 2) == ghost_width .and. &
      size(patch%left_ghost_temperature) == ghost_width .and. &
      size(patch%right_ghost_temperature) == ghost_width
  end function sparse_patch_has_shape

  pure logical function sparse_patch_is_empty(patch) result(empty)
    type(amr_patch_tree_reactive_patch_1d), intent(in) :: patch

    empty = .not. allocated(patch%state) .and. &
      .not. allocated(patch%temperature) .and. &
      .not. allocated(patch%left_ghost_state) .and. &
      .not. allocated(patch%right_ghost_state) .and. &
      .not. allocated(patch%left_ghost_temperature) .and. &
      .not. allocated(patch%right_ghost_temperature)
  end function sparse_patch_is_empty

  pure integer function sparse_patch_value_count(patch) result(count)
    type(amr_patch_tree_reactive_patch_1d), intent(in) :: patch

    count = size(patch%state) + size(patch%temperature) + &
      size(patch%left_ghost_state) + size(patch%right_ghost_state) + &
      size(patch%left_ghost_temperature) + &
      size(patch%right_ghost_temperature)
  end function sparse_patch_value_count

  subroutine all_ranks_accept_sparse_1d( &
      distribution, local_value, global_value, ok)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    logical, intent(in) :: local_value
    logical, intent(out) :: global_value, ok

    integer :: ierr

    call MPI_Allreduce( &
      local_value, global_value, 1, MPI_LOGICAL, MPI_LAND, &
      distribution%comm, ierr)
    ok = ierr == MPI_SUCCESS
    if (.not. ok) global_value = .false.
  end subroutine all_ranks_accept_sparse_1d

end module mpi_amr_sparse_patch_1d_mod
